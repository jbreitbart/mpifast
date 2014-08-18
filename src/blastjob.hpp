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
#ifndef __BLASTJOB_HPP__
#define __BLASTJOB_HPP__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mpi.h"

#include <string>
#include <set>
#include <vector>

#include "mpiblast_types.h"

/**
 * Tracks which nodes have which fragments.  A single instance of this
 * object should be created and used to track the DB distribution 
 * among all participating workers
 */
class FragmentDistribution {
public:
	FragmentDistribution(){};
	/**
	 * Construct an empty fragment distribution with a given number of nodes and fragments
	 */
	FragmentDistribution( uint fragment_count, uint node_count );
	FragmentDistribution( const FragmentDistribution& fd );
	FragmentDistribution& operator=( const FragmentDistribution& fd );
	
	/** Sets the fragments on a node's local storage */
	void update( int nodeI, const std::set<int>& fragment_list );
	
	/** Prints the node_fragments to cout */
	void printNodeFrags( int nodeI ) const;
	
	/** Returns the set of fragments on a given node */
	const std::set<int>& getNodeFrags( int nodeI ) const;

	/** Returns the set of nodes that have a given fragment */
	const std::set<int>& getFragNodes( int fragmentI ) const;
	
	/** Return the set of fragments currently being copied */
	const std::set<int>& getCopyFrags() const;

	/** Register that a fragment is being copied */
	void copyFragment( int fragmentI );

	/** Unregister a fragment copy that has completed, add it to the node's list */
	void copyCompleted( int nodeI, int fragment_id );

	/**
	 * Given a set of fragments, find the least distributed fragment
	 */
	int getLeastDistributed( const std::set< int >& possible_fragments ) const;
protected:
	std::vector< std::set< int > > node_fragments;	/**< The set of fragments residing on each node -- index[node][fragment] */
	std::vector< std::set< int > > fragment_nodes;	/**< The set of nodes each fragment resides on -- index[fragment][node] */
	std::set< int > copy_fragments; /**< Fragments that are currently being copied somewhere */
};


/**
 * Manages each Blast job by allocating work to worker nodes according to a policy
 * A "blast job" can be a single query or a set of queries.
 * Each BlastJob contains a reference to a single instance of the FragmentDistribution
 * object for all participating workers.  The BlastJob class can then make fragment
 * search assignments to workers based on knowledge of the database fragment 
 * distribution.
 * 
 */
class BlastJob{
	public:
		BlastJob( int fragment_count, FragmentDistribution& frag_distribution );
		BlastJob( const BlastJob& bj );
		BlastJob& operator=( const BlastJob& bj );
		
		/**
		 * Calculates the best assignment for a node based on the current
		 * distribution of database fragments
		 */
		// virtual void getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses );
		virtual void getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses, bool& should_copy);
		int GetNumUnassigned() { return unassigned_fragments.size(); }

	protected:
		int fragment_count;	/**< The total number of fragments for this database */
		std::set< int > unassigned_fragments;	/**< Fragments that haven't been assigned yet */
		FragmentDistribution& frag_distribution;	/**< Tracks which nodes have which DB fragments */
};

/** 
 * A scheduler policy that precopies a minimal database 
 * distribution prior to making search assignments
 */
class PrecopySchedulerPolicy : public BlastJob
{
public:
	PrecopySchedulerPolicy( int fragment_count, FragmentDistribution& frag_distribution );
	PrecopySchedulerPolicy( const PrecopySchedulerPolicy& bj );
	PrecopySchedulerPolicy& operator=( const PrecopySchedulerPolicy& bj );
	//virtual void getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses );
	virtual void getAssignment( int nodeI, int& operation, int& fragment_id, int& num_copying, int concurrent_accesses, bool& should_copy);
	static int db_copy_number;	/**< The number of complete copies of the database that will be created across cluster nodes */
	static int node_copy_number;	/**< The number of database fragments that will cached by a node */
};

#endif	//__BlastJob_h__
