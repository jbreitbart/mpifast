/*************************************************************************
** Copyright 2009 by Virginia Polytechnic Institute and State
** University. All rights reserved. Virginia Polytechnic Institute and
** State University (Virginia Tech) owns the mpiBLAST software and its
** associated documentation ("Software"). You should carefully read the
** following terms and conditions before using this software. Your use
** of this Software indicates your acceptance of this license agreement
** and all terms and conditions.
** 
** The parallel input/output (PIO) portion of this work was based in
** part on published research that was supported by the North Carolina
** State University and Oak Ridge National Laboratory. The actual
** implementation was completed at Virginia Tech in the summer of 2007.
** 
** Additionally, portions of this work are derivative works of the NCBI
** C Toolkit. Although, the NCBI C Toolkit is released freely without
** restriction under the terms of their license, the following files
** listed below, have been modified by Virginia Tech, and as such are
** redistributed under the terms of this license.
** 
** ncbi/api/txalign.c
** ncbi/corelib/ncbifile.c
** ncbi/object/objalign.c
** ncbi/tools/blast.c
** ncbi/tools/blastdef.h
** ncbi/tools/blastool.c
** ncbi/tools/blastutl.c
** ncbi/tools/blfmtutl.c
** ncbi/tools/ncbisam.c
** ncbi/tools/readdb.c
** ncbi/tools/readdb.h
** 
** License:
** 
** This file is part of mpiBLAST.
** 
** mpiBLAST is free software: you can redistribute it and/or modify it 
** under the terms of the GNU General Public License version 2 as published 
** by the Free Software Foundation. 
** 
** Accordingly, mpiBLAST is distributed in the hope that it will be
** useful, but WITHOUT ANY WARRANTY; without even the implied warranty
** of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
** General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with mpiBLAST. If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/
#include "blastjob.hpp"
#include "mpiblast.hpp"
#include <algorithm>
using namespace std;

BlastJob::BlastJob( int frag_count, FragmentDistribution& frag_dist ) :
fragment_count( frag_count ),
frag_distribution( frag_dist )
{	
	// init copy_fragments and unassigned_fragments
	uint fragmentI;
	for( fragmentI = 0; fragmentI < (uint)fragment_count; fragmentI++ ){
		unassigned_fragments.insert( fragmentI );
	}
}

BlastJob::BlastJob( const BlastJob& bj ) :
fragment_count( bj.fragment_count ),
frag_distribution( bj.frag_distribution ),
unassigned_fragments( bj.unassigned_fragments )
{
}

BlastJob& BlastJob::operator=( const BlastJob& bj ){
	BlastJob( bj.fragment_count, bj.frag_distribution );
	unassigned_fragments = bj.unassigned_fragments;
	return *this;
}

void BlastJob::getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses, bool& should_copy){
	if( debug_msg )
		LOG_MSG << "Computing next assignment for node " << nodeI << endl;

	// first try assigning the least distributed unassigned fragment that
	// already exists on nodeI
	// intersect unassigned with the set of frags on nodeI
	set< int > nodeI_unassigned;
	insert_iterator< set< int > > node_ins(nodeI_unassigned, nodeI_unassigned.begin());
	set_intersection( frag_distribution.getNodeFrags( nodeI ).begin(), frag_distribution.getNodeFrags( nodeI ).end(),
			unassigned_fragments.begin(), unassigned_fragments.end(), node_ins );
	
	// print some debugging info
	if( debug_msg ){
		LOG_MSG << "Node " << nodeI << " fragments: ";
		frag_distribution.printNodeFrags( nodeI );
		CONT_LOG_MSG << endl;
		
		LOG_MSG << "Unassigned fragments: ";
		set< int >::const_iterator iter = unassigned_fragments.begin();
		for( ; iter != unassigned_fragments.end(); iter++ ){
			CONT_LOG_MSG << " " << *iter;
		}
		CONT_LOG_MSG << endl;

		LOG_MSG << "Intersection fragments: ";
		iter = nodeI_unassigned.begin();
		for( ; iter != nodeI_unassigned.end(); iter++ ){
			CONT_LOG_MSG << " " << *iter;
		}
		CONT_LOG_MSG << endl;
	}


	// get the least distributed
	int search_fragI = frag_distribution.getLeastDistributed( nodeI_unassigned );	
	if( debug_msg )
		LOG_MSG << "search_fragI is " << search_fragI << endl;

	if( search_fragI != -1 ){
		unassigned_fragments.erase( search_fragI );
		operation = SEARCH_FRAGMENT;
		fragment_id = search_fragI;
		return;
	}
	
	// if -1 was returned (the intersected set was empty) then we need to copy a frag.
	// find the least distributed unassigned fragment to copy
	// take the set difference of unassigned and copy_fragments
	
	// WARNING!  This code only permits one node to copy a particular DB fragment at a time
	// This behavior is fine so long as query segmentation is not implemented
	set< int > uncopy;
	insert_iterator< set< int > > unass_ins(uncopy, uncopy.begin());
	set_difference( unassigned_fragments.begin(), unassigned_fragments.end(),
			frag_distribution.getCopyFrags().begin(), frag_distribution.getCopyFrags().end(),
			unass_ins );

	int copy_fragI = frag_distribution.getLeastDistributed( uncopy );
	
	if( copy_fragI != -1 ){
		operation = COPY_FRAGMENT;
		fragment_id = copy_fragI;
		// JA -- Once a node copies a fragment, it immediately starts searching it.
		//       This bypasses the Worker asking for an assignment, which we know will be SEARCH_FRAGMENT
//		unassigned_fragments.erase( fragment_id );
		frag_distribution.copyFragment( copy_fragI );
		return;
	}
	
	// if all of the fragments waiting to be assigned are already being copied
	// then this worker can't do any work on this query and should move on to the next
	// query.  return a state of DO_NOTHING
	fragment_id = 0;
	if( frag_distribution.getCopyFrags().size() > 0 )
		operation = DO_NOTHING;
	else
		// if -1 was returned (the intersected set was empty) then the search must be complete
		operation = SEARCH_COMPLETE;

}

int PrecopySchedulerPolicy::db_copy_number = 1;
int PrecopySchedulerPolicy::node_copy_number = 0;

PrecopySchedulerPolicy::PrecopySchedulerPolicy( int fragment_count, FragmentDistribution& frag_distribution )
: BlastJob( fragment_count, frag_distribution )
{
}

PrecopySchedulerPolicy::PrecopySchedulerPolicy( const PrecopySchedulerPolicy& bj )
: BlastJob( bj ) {}

PrecopySchedulerPolicy& PrecopySchedulerPolicy::operator=( const PrecopySchedulerPolicy& bj )
{
	BlastJob::operator =( bj );
	return *this;
}

void PrecopySchedulerPolicy::getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses, bool& should_copy )
{
	static set< int > all_frags;
	if(all_frags.empty()) {
		for( uint fragI = 0; fragI < fragment_count; fragI++ )
			all_frags.insert( fragI );
	}

	set< int > cur_frags = all_frags;
	const set< int >& cf = frag_distribution.getCopyFrags();
	// if this worker has no fragment, tell it to copy the least distributed
	// fragment
	if( frag_distribution.getNodeFrags( nodeI ).size() == 0 ){
		int copy_fragI = -1;
		int least_count = -1;
		while( copy_fragI == -1 ){
			if(cur_frags.empty()) 
				break;

			int least_frag = frag_distribution.getLeastDistributed( cur_frags );
			if( least_count == -1 )
				least_count = frag_distribution.getFragNodes(least_frag).size();
			if( least_count < frag_distribution.getFragNodes(least_frag).size() )
				break;	// failed to make an assignment
			if( cf.find( least_frag ) != cf.end() ){
					cur_frags.erase( least_frag );
					continue;	// this one is being copied already.  keep looking.
			}
			copy_fragI = least_frag;
		}
		if( copy_fragI != -1 ){
			operation = COPY_FRAGMENT;
			fragment_id = copy_fragI;
			frag_distribution.copyFragment( copy_fragI );
			num_copying++;
			return;
		}
		// if the minimally distributed fragments are already being copied
		// tell this worker to do nothing for the time being
		fragment_id = 0;
		if( frag_distribution.getCopyFrags().size() > 0 )
			operation = DO_NOTHING;
		return;
	}

	// check if there's any fragment to search first
	// from here on out, assign the least distributed fragment to get searched

	// intersect unassigned with the set of frags on nodeI
	set< int > nodeI_unassigned;
	insert_iterator< set< int > > node_ins(nodeI_unassigned, nodeI_unassigned.begin());
	set_intersection( frag_distribution.getNodeFrags( nodeI ).begin(), frag_distribution.getNodeFrags( nodeI ).end(),
			unassigned_fragments.begin(), unassigned_fragments.end(), node_ins );

	// get the least distributed
	int search_fragI = frag_distribution.getLeastDistributed( nodeI_unassigned );	

	if( search_fragI != -1 ){
		unassigned_fragments.erase( search_fragI );
		operation = SEARCH_FRAGMENT;
		fragment_id = search_fragI;
		return;
	}
	
	// ensure the fragment distribution is complete
	// if an undistributed fragment is not already being copied
	// then tell this worker to copy it!

	cur_frags = all_frags;
	int copy_fragI = frag_distribution.getLeastDistributed( cur_frags );
	while( 
		frag_distribution.getFragNodes( copy_fragI ).size() < db_copy_number &&
		(cf.find( copy_fragI ) != cf.end() || 
		frag_distribution.getFragNodes( copy_fragI ).find( nodeI ) !=
		frag_distribution.getFragNodes( copy_fragI ).end() )
		)
	{
			cur_frags.erase( copy_fragI );
			if( cur_frags.size() == 0 ){
				copy_fragI = -1;
				break;
			}
			copy_fragI = frag_distribution.getLeastDistributed( cur_frags );
			continue;	// this one is being copied already.  keep looking.
	}
	if( copy_fragI != -1 && frag_distribution.getFragNodes( copy_fragI ).size() < db_copy_number ){
		bool do_copy = false;
		// node_copy_number == 0 means unset, do the copy
		if(node_copy_number == 0 || frag_distribution.getNodeFrags( nodeI ).size() < node_copy_number) { 
			do_copy = true;
		}

		if(do_copy) {
			if(num_copying < concurrent_accesses) {
				num_copying++;
				// this guy needs to get copied
				operation = COPY_FRAGMENT;
				fragment_id = copy_fragI;
				frag_distribution.copyFragment( copy_fragI );
				return;
			} else {
				operation = DO_NOTHING; // avoid have workers' searching query backward
				should_copy = true;
				return;
			}
		}
	}

	fragment_id = 0;
	if( unassigned_fragments.size() == 0 )
		operation = SEARCH_COMPLETE;
	else
		operation = DO_NOTHING;
	return;
}

FragmentDistribution::FragmentDistribution( uint fragment_count, uint node_count ) {

	set< int > frag_set;
	fragment_nodes = vector< set< int > >( fragment_count );
	node_fragments = vector< set< int > >( node_count );
	
}

FragmentDistribution::FragmentDistribution( const FragmentDistribution& fd ) {
	*this = fd;
}

FragmentDistribution& FragmentDistribution::operator=( const FragmentDistribution& fd ) {
	node_fragments = fd.node_fragments;
	fragment_nodes = fd.fragment_nodes;
	copy_fragments = fd.copy_fragments;
	return *this;
}

void FragmentDistribution::update( int nodeI, const std::set< int >& fragment_list ){
	node_fragments[ nodeI ] = fragment_list;
	
	// remove this node from all of the fragment lists
	for( int fragI = 0; fragI < fragment_nodes.size(); fragI++ ){
		fragment_nodes[ fragI ].erase( nodeI );
	}
	
	// add this node to all of the correct fragment node lists
	set< int >::const_iterator frag_iter = fragment_list.begin();
	for( ; frag_iter != fragment_list.end(); frag_iter++ ){
		fragment_nodes[ *frag_iter ].insert( nodeI );
	}
}

void FragmentDistribution::printNodeFrags( int nodeI ) const {
	set< int >::const_iterator iter = node_fragments[ nodeI ].begin();
	for( ; iter != node_fragments[ nodeI ].end(); iter++ ){
		CONT_LOG_MSG << " " << *iter;
	}
}

const set<int>& FragmentDistribution::getNodeFrags( int nodeI ) const {
	return node_fragments[ nodeI ];
}

const set<int>& FragmentDistribution::getFragNodes( int fragmentI ) const {
	return fragment_nodes[ fragmentI ];
}

const set<int>& FragmentDistribution::getCopyFrags() const {
	return copy_fragments;
}

void FragmentDistribution::copyFragment( int fragmentI ) {
	copy_fragments.insert( fragmentI );
}

void FragmentDistribution::copyCompleted( int nodeI, int fragment_id ) {
	fragment_nodes[ fragment_id ].insert( nodeI );
	copy_fragments.erase( fragment_id );
	node_fragments[ nodeI ].insert( fragment_id );
	
	if(debug_msg){
		LOG_MSG << "node " << nodeI << " copied fragment " << fragment_id << endl;
	}
}


int FragmentDistribution::getLeastDistributed( const set< int >& possible_fragments ) const {
	// find the lowest replication fragment in possible_fragments

	int best_fragment = -1;
	int lowest_replication = -1;
	set< int >::const_iterator frag_iter = possible_fragments.begin();
	
	for(; frag_iter != possible_fragments.end(); frag_iter++ ){
		int fragmentI = *frag_iter;
		if( (fragment_nodes[ fragmentI ].size() < lowest_replication) || 
		   (lowest_replication == -1) ){
			best_fragment = fragmentI;
			lowest_replication = fragment_nodes[ fragmentI ].size();
		}
	}
	return best_fragment;
}


