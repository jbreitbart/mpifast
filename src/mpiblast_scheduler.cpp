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
#include "mpiblast_scheduler.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_tags.hpp"
#include "file_util.hpp"
#include "query_manager.hpp"

static WriterMaster* writer_m = NULL;

Scheduler::Scheduler(unsigned int atotal_frags, int anode_count, int aquery_count, int aconcurrent_accesses, int acopy_via) 
	: frag_dist(atotal_frags, anode_count) {
	total_frags = atotal_frags;
	my_node_count = anode_count;
	query_count = aquery_count;
	finished_nodes = 0;
	current_working_queryI = 0;
	output_queryI = 0;
	num_copying_db = 0;
	no_more_queries = false;
	concurrent_accesses = aconcurrent_accesses;
	copy_via = acopy_via;
	prefetching = false;
	writer_m = NULL;
	num_copying = 0;

	worker_last_query = vector<int>(my_node_count, -1);

// initialize a FragmentDistribution which tracks the distribution of
// fragments on worker nodes
//	FragmentDistribution frag_dist( total_frags, my_node_count );

	crsp_prefetch = new CommRecvStruct(5 * int_size);
}

void Scheduler::RecvFragmentLists( void ) {
	//
	// refresh each worker's local storage fragment list
	//
	for( uint workerI = 0; workerI < my_node_count; workerI++ ){
		if(!IsWorker(workerI)) {
			continue;
		}

		vector< int > fvec;
		RecvIntVec(group_comm, fvec, workerI, DB_FRAGMENTS_TAG);
		set< int > fset = set< int >( fvec.begin(), fvec.end() );
		if( debug_msg ){
			LOG_MSG << "Node " << workerI;
			if( fvec.size() > 0 )
				CONT_LOG_MSG << " has fragments ";
			else
				CONT_LOG_MSG << " has no fragments";
			for( int fI = 0; fI < (int)fvec.size(); fI++ ){
				CONT_LOG_MSG << " " << fvec[ fI ];
			}
			CONT_LOG_MSG << endl;
		}
		// Update the FragmentDistribution
		frag_dist.update( workerI, fset );

		if( debug_msg ){
			LOG_MSG << "Received fragment list update from node " << workerI << endl;
		}
	}
}

void Scheduler::CheckPendingSends() {
	list <CommSendStruct*>::iterator it;
	list <CommSendStruct*>::iterator tmp_it;
	
	it=_send_ops_list.begin();
	while(it!=_send_ops_list.end()) {
		if((*it)->TestIsend()) {
			delete(*it);
			tmp_it = it;
			it++;
			_send_ops_list.erase(tmp_it);
		} else {
			it++;
		}
	}
}

// wait for non-blocking sends to be finished
void Scheduler::WaitPendingSends() {
	list <CommSendStruct*>::iterator it;
	list <CommSendStruct*>::iterator tmp_it;
	
	it=_send_ops_list.begin();
	while(it!=_send_ops_list.end()) {
		(*it)->WaitIsend();
		delete(*it);
		tmp_it = it;
		it++;
		_send_ops_list.erase(tmp_it);
	}
}

void Scheduler::Finalize() {
	WaitPendingSends();
	delete crsp_prefetch;

	int assign_array[5];
	assign_array[0] = WORKER_QUIT;
	
	for(int i=0; i<my_node_count; i++) {
		if(IsWorker(i)) {
			MPI_Send( &assign_array, 5, MPI_INT, i, ASSIGNMENT_TYPE, group_comm);
		}
	}
}

int Scheduler::HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) {
	scheduler_called++;
	
	CheckPendingSends();
	static int complete_nodes = 0;

	int worker_id = recv_status.MPI_SOURCE;

	IssuePrefetchRequest();
	
	if(recv_status.MPI_TAG == FRAGMENT_COPY_COMPLETE) {
		int new_fragment = -1;
		crsp->ExtractData(&new_fragment, int_size);
		frag_dist.copyCompleted( worker_id, new_fragment );
		num_copying--;
	} else if (recv_status.MPI_TAG == WORKER_IDLE) {
		int assign_array[5];
		assign_array[0] = DO_NOTHING;

		// loop until we get an assignment or until we 
		// exhaust the list of unsearched queries
		int working_queryI = current_working_queryI;
		while( true ) {
			
			if(!no_more_queries) {
				if(working_queryI >= working_queries.size()) {
					CheckRecvQuerySegment(true);
				} else {
					CheckRecvQuerySegment(false);
				}
			}

			if(working_queryI >= working_queries.size()) {
				break;
			}

			bool should_copy = false;
			if(working_queryI < working_queries.size()) {
//				query_data[ working_queries[working_queryI] ].blast_job->getAssignment( worker_id, assign_array[0], assign_array[2], num_copying, concurrent_accesses );
				QueryData* qdp = query_data[ working_queries[working_queryI] ];
				qdp->blast_job->getAssignment( worker_id, assign_array[0], assign_array[2], num_copying, concurrent_accesses, should_copy);
			}

			// if current_queryI has been completely assigned, advance it to the next query
			if( working_queryI == current_working_queryI && assign_array[0] == SEARCH_COMPLETE ) {
				delete query_data[current_working_queryI];
				query_data.erase(current_working_queryI);
				query_map.erase(current_working_queryI);
				current_working_queryI++;
			}

			if(assign_array[0] != DO_NOTHING && assign_array[0] != SEARCH_COMPLETE){
				// get an real assignment
				break;
			}

			if(assign_array[0] == DO_NOTHING && should_copy) {
				// this worker should copy but exceeding the concurrent copy limit, break to avoid backward assignment
				break;
			}
			
			working_queryI++;
		}
		
		if( working_queryI < working_queries.size()) {
			assign_array[3] = working_queries[working_queryI]; // query_id
			assign_array[1] = query_map[assign_array[3]]; // segment_id
		} else {
			assign_array[3] = -1;
			assign_array[1] = -1;
		}

		if( no_more_queries && current_working_queryI >= working_queries.size() ) {
			assign_array[0] = SEARCH_COMPLETE;
			complete_nodes++;
		}
			
		if( debug_msg ) {
			LOG_MSG << "Giving " << worker_id << " assignment " << assign_array[0] << " on frag " << assign_array[2] << " query_id " << assign_array[3] << " segment_id " << assign_array[1] << endl;
		}
		
		if (assign_array[0] == COPY_FRAGMENT){
			if (debug_msg)
				LOG_MSG << "Telling node " << worker_id << " to copy fragment " << assign_array[2] << endl;
		}

		assign_array[4] = -1; // the length of next message of query data

		int query_len = 0;
		char* query_buf = NULL;
		if(!query_in_mem && assign_array[0] == SEARCH_FRAGMENT) {
			if(assign_array[3] > worker_last_query[worker_id]) { // worker does not have this query
				query_len = QueryM::Instance()->GetQueryLen(assign_array[3]);
				query_buf = new char[query_len];
				QueryM::Instance()->CloneQueryData(assign_array[3], query_buf);
				
				// need to send the query data
				assign_array[4] = query_len;
				worker_last_query[worker_id] = assign_array[3]; 
			}
		} 

		// use isend to avoid blocking scheduler
		CommSendStruct* cssp_assign = new CommSendStruct(5 * int_size);
		cssp_assign->AddData(&assign_array, 5 * int_size);
		cssp_assign->IsendData(group_comm, worker_id, ASSIGNMENT_TYPE);
		AddPendingSend(cssp_assign);

		if(query_len > 0) {
			// send message to worker
			CommSendStruct* cssp_query = new CommSendStruct(query_len);
			cssp_query->AddData(query_buf, query_len);
			cssp_query->IsendData(group_comm, worker_id, QUERY_DATA);
			AddPendingSend(cssp_query);
		}

		if(query_buf != NULL) delete query_buf;

	} else {
		throw __FILE__ "Scheduler::HandleMessages(): Illegal tag";
	} 

	delete crsp;

	IssuePrefetchRequest();

	if(complete_nodes == GetNumWorkers()) {
		if(debug_msg) {
			LOG_MSG << "Scheduler end" << endl;
		}
		return 1;
	}

	return 0;
}

int Scheduler::IssuePrefetchRequest() {
	if(prefetching || no_more_queries) {
		return 0;
	}
	
	int total_unassigned = 0;
	for(int i=current_working_queryI; i<working_queries.size(); i++) {
		total_unassigned += query_data[ working_queries[i] ]->blast_job->GetNumUnassigned();
	}

	// working: add prefetch
	if((total_unassigned < 2 * GetNumWorkers()) && !no_more_queries) {
		if(debug_msg) {
			LOG_MSG << "Issue prefetch of a query segment" << endl;
		}
		
		int event_array[2];
		event_array[0] = QUERY_SEGMENT_REQUEST;
		event_array[1] = -1;
		CommSendStruct* cssp_prefetch = new CommSendStruct(2 * int_size);
		cssp_prefetch->AddData(&event_array, 2 * int_size);
		cssp_prefetch->IsendData(MPI_COMM_WORLD, super_master_process, EVENT_TYPE);
		AddPendingSend(cssp_prefetch);
		
		prefetching = true;

		crsp_prefetch->IrecvData(MPI_COMM_WORLD, super_master_process, ASSIGNMENT_TYPE);
	}

	return 0;
}

int Scheduler::CheckRecvQuerySegment(bool wait) {
	if(!prefetching) return 0;

	bool received = false;

	if(wait) {
		double start_t = MPI_Wtime();
		crsp_prefetch->WaitIrecv();
		received = true;
		acquire_segement_time += MPI_Wtime() - start_t;
	} else {
		if(crsp_prefetch->TestIrecv()) {
			received = true;
		}
	}

	if(received) {
		int command = -1;
		int segment_id = -1;
		int start_query = -1;
		int end_query = -1;

		crsp_prefetch->ExtractData(&command, int_size);
		crsp_prefetch->ExtractData(&segment_id, int_size);
		crsp_prefetch->ExtractData(&start_query, int_size);
		crsp_prefetch->ExtractData(&end_query, int_size);
		crsp_prefetch->ResetDataBuf();  
		
		if(command == SEARCH_SEGMENT) {
			AddAssignment(segment_id, start_query, end_query);
			
			if(debug_msg) {
				LOG_MSG << "Got segment " << segment_id << " from supermaster" << endl;
			}
		} else if (command == NO_MORE_SEGMENT) {
			no_more_queries = true;
			writer_m->SetNoMoreQueries();
		} else {
			throw __FILE__ "Scheduler::CheckRecvQuerySegment - received invalid message from supermaster";
		}

		prefetching = false;
	}

	return 0;
}

int Scheduler::AddAssignment(int segment_id, int start_query, int end_query) {
	if(debug_msg) {
		LOG_MSG << "Adding assignments: start_query " << start_query << " end_query " << end_query << endl;
	}
	for(int query_id = start_query; query_id < end_query; query_id++) {
		working_queries.push_back(query_id);
		query_map[query_id] = segment_id;

		query_data[query_id] = new QueryData(NULL, NULL);
		query_data[query_id]->blast_job = new PrecopySchedulerPolicy( total_frags, frag_dist ); 
		QueryM::Instance()->GetQueryEntry(query_id);
	}

	// add working queries for writer
	writer_m->AddWorkingQueries(segment_id, start_query, end_query);

	// load queries
	if(!query_in_mem) {
		QueryM::Instance()->LoadQueries(start_query, end_query);
	}

	return 0;
}
