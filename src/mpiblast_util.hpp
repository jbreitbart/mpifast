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
#ifndef __MPI_BLAST_UTIL
#define __MPI_BLAST_UTIL

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <iterator>
#include <fcntl.h>
#include <assert.h>
#include "mpi.h"
#include "mpiblast_types.h"

// Logging for debugging and profiling: All variables are declared in mpiblast_util.cpp
extern bool debug_msg;		/* Default, set to true with --debug */
extern std::ostream* log_stream; 		/**< Output stream for logging information */
extern double realstarttime ;		/* Given value immediately on first line of main, before MPI_Init */
extern double prog_start;		/* Given value after MPI_Init and MPE defs */
extern double prog_end ;		/* Given value right before MPI_Finalize */
extern int rank;			/**< Rank of the current MPI process */
extern int node_count;		/**< Number of MPI processes allocated for the job */
extern int group_rank;  // rank within a group
extern int group_node_count; // how manay procs in my group
extern int writer_process;
extern int scheduler_process;
extern int super_master_process;
extern bool use_master_write;
extern int max_write_buffer;
extern bool fast_evalue_approximation;
extern int parallel_evalue_adjust;
extern int output_strategy;
extern int io_function;
extern int ncpus;
extern int concurrent_write;
extern int query_segment_size;
extern int max_segment_size;
extern double init_assign_percent;
extern bool dump_raw_output;
extern bool disable_posix_lock;
extern bool sync_comm;
extern bool use_real_dblen;
extern int use_query_map;
extern int query_in_mem;
extern int preload_query;
extern int num_pending_writes;
extern int num_query_segments;
extern int use_brief_report;
// extern int use_virtual_frags; /* 0-false; 1-true */
extern MPI_Comm group_write_comm;
extern MPI_Comm group_comm;

//extern MPI_File file_handle;
extern uint64 max_data_to_send;	/**< max data to send when using MPI to copy fragments */

// debug
extern FILE* dbgfp;
extern double search_time;
extern double copy_time;
extern double result_time;
extern double run_blast_time;
extern double write_time;
extern double output_info_time;
extern double query_broadcast_time;
extern double idle_time;
extern double master_idle_time;
extern double scheduler_idle_time;
extern double writer_master_idle_time;
extern double scheduler_handle_time;
extern double writer_master_handle_time;
extern double writer_worker_handle_time;
extern double results_gather_time;
extern double wait_result_time;
extern double wait_write_time;
extern double receive_msg_time;
extern double receive_msg_size;
extern double db_distribute_time;
extern double acquire_segement_time;
extern double scheduler_called;
extern double load_queries_time;
extern double default_output_time;
extern double process_output_time;
extern double curr_process_output_time;

extern int max_num_write;
extern int peak_pending_offsets;  

extern int max_pending_offsets; 

extern int num_threads;
extern int num_reserved_events;
extern int max_query_load;

extern const int int_size;
extern const int float_size;
extern const int offset_size;

#ifdef USING_MPI
#include <iomanip>
#define LOG_MSG (*log_stream) << "[" << rank << "]\t" << setiosflags(ios::fixed) << setw(20) << setprecision(6) << MPI_Wtime() - prog_start << '\t'
#else
#define LOG_MSG (*log_stream)
#endif // USING_MPI

#define CONT_LOG_MSG (*log_stream) 

/**
 * Initializes the NCBI library with a particular vector of options
 */
void initNCBI( std::vector< std::string >& ncbi_opts );
void cleanupNCBI();

/**
 * Add a command line option to an option vector
 */
void addOpt( std::vector< std::string >& opt_vector, int opt_character, const char*
	opt_argument );

void SendIntVec(MPI_Comm comm, std::vector<int> &vec, int dest, int tag);

void RecvIntVec(MPI_Comm comm, std::vector<int> &vec, int src, int tag);

void printProgress (uint64 progress_sofar, float64 total);

int getlock (int fd, int cmd, int type);


void broadcastFile( std::string& b_filename );

void recvBroadcastFile( std::string& b_filename );


void sendFile( const std::string& filename, int dest, int tag );

void recvFile( const std::string& filename, int src, int tag );

int IsWorker(int src);

int GetNumWorkers();

class GroupManager {
	public:
		GroupManager(int group_size, int num_frags);

		static GroupManager* Instance(int group_size, int num_frags) {
			if(_instance == NULL) {
				_instance = new GroupManager(group_size, num_frags);
			}

			return _instance;
		}
		
		static GroupManager* Instance() {
			return _instance;
		}

		int GetGroupId(int world_rank) { 
			return _proc_group[world_rank]; 
		}
		int Rank2GroupRank(int world_rank); // need to add one map
		int GroupRank2Rank(int group_rank, int group_id); 
		int GetNumGroups() { return _num_groups; }
		int GetGroupSize(int group_id) { return _actual_group_size[group_id]; }
		void PrintGroupMap(std::ostream& os);

	private:
		int _group_size;
		int _num_groups;
		std::vector <int> _actual_group_size;
		std::map <int, int> _proc_group; // map proc to group_id
		static GroupManager* _instance;
};


#endif // __MPI_BLAST_UTIL
