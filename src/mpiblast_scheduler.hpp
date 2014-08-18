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
#ifndef __MPIBLAST_SCHEDULER_HPP__
#define __MPIBLAST_SCHEDULER_HPP__

#include "mpiblast_querydata.hpp"
#include "blastjob.hpp"
#include "mpiblast_writer.hpp"
#include <vector>
#include <queue>
using namespace std;

class Scheduler {
public:
	Scheduler(unsigned int atotal_frags, int anode_count, int aquery_count, int aconcurrent_accesses, int acopy_via);
	void RecvFragmentLists( void );
	int HandleMessages(MPI_Status &event_status, int* event_array);
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status);
	void AddPendingSend(CommSendStruct* cssp) {
		_send_ops_list.push_back(cssp);
	}
	void CheckPendingSends();
	void WaitPendingSends();
	void SetWriterMaster(WriterMaster* writer) { writer_m = writer; }
	int GetNumWorkingQueries() { return working_queries.size(); }
	inline int IssuePrefetchRequest();
	inline int CheckRecvQuerySegment(bool wait);
	int AddAssignment(int segment_id, int start_query, int end_query);
	void Finalize(void);

private:
	int finished_nodes;
	// int current_queryI;
	int current_working_queryI;
	int query_count;
	int my_node_count;
	int num_copying;
	int writer_index;

	int num_copying_db; // JA -- number of workes currently copying DB fragments, value must be less then concurrent_accesses
	int output_queryI;	/**< Stores the index of the next query to output */
	unsigned int total_frags;
	int concurrent_accesses;	/**< number of concurrent accesses to shared storage */
	int copy_via;			/**< how to copy the fragments */
	bool no_more_queries;
	bool prefetching;
	CommRecvStruct* crsp_prefetch;
	WriterMaster* writer_m;

	// vector< QueryData > query_data;
	map <int, QueryData*> query_data; // key: query_id, data: query_data
	queue< pair<int, int> > db_access_queue; // JA DB access queue of workers waiting to access the database, pair values are <worker rank, fragment id>
	FragmentDistribution frag_dist;

	list <CommSendStruct*> _send_ops_list; // send operations list 
	vector <int> working_queries; // queries that have been assigned to this group
	map <int, int> query_map; // key: query_id, data: segment_id
	vector <int> worker_last_query;
};

#endif
