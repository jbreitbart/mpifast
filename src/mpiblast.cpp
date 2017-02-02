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
#include "mpiblast.hpp"

extern "C" {
#include "getopt.h"
#include "mpi_util.h"
#include "blast_hooks.h"
#include "readdb.h"
#include "distributed_bioseq.h"
}
#include "virtual_files.hpp"
#include "query_manager.hpp"
#include <iterator>

#ifdef WIN32
	const char* bitbucket = "NUL";
#else
	const char* bitbucket = "/dev/null";
#endif

using namespace std;

#ifdef MPE
#include <mpe.h>
#endif
int cpstart, cpend, blaststart, blastend, mergestart, mergeend; //MPE ids
int outputstart, outputend, recvbioseqstart, recvbioseqend; //more MPE ids

// These MPE ids were added to replace the profile_msg functionality
int init_blast_start, init_blast_end, load_queries_start, load_queries_end;
int cleanup_blast_start, cleanup_blast_end, send_results_start, send_results_end;
int send_frag_list_start, send_frag_list_end, file_cleanup_start, file_cleanup_end;
int bcast_filter_start, bcast_filter_end;
int bcast_query_start, bcast_query_end, bcast_qadj_start, bcast_qadj_end;
int receive_results_start, receive_results_end, worker_file_setup_start, worker_file_setup_end;
int init_ncbi_start, init_ncbi_end;

#ifndef WIN32
#include <signal.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

// File local variables functions
string localPath;	/**< Used by cleanupDB() to give the local storage path */
static MpiBlast* mpiblast_ptr = NULL;
static bool running_blast = false;	/**< set to true when executing code in the NCBI library */
static bool sameSharedLocal = false;	/**< set to true when shared and local storage are identical */
static bool replica_count_set = false;
static string global_output_file;
static string profile_filename;
static int partition_size;
static int replica_group_size = 0;
static int frags_per_worker = 0;
static int query_file_len = 0;
static bool predistribute_db = false;
static bool query_segment_size_set = false;
static bool increment_db_copy = true;
static double system_init_time = 0;
static string hpm_filename;
bool collect_profile = false;
void exitCleanly();
void terminateProgram( int sig );
static void CollectSearchProfile(const string& file_name, double total_time);

// DEBUG --  Electric fence global variables
//extern int EF_ALLOW_MALLOC_0;
//extern int EF_PROTECT_BELOW;

int main( int argc, char* argv[] )
{
	if(init_intercept_lib())  { // has to be called before any fwrite
		cerr<<"Error calling init_intercept_lib()";
		exit(1);
	}

	// DEBUG --  Electric fence global variables
	//EF_ALLOW_MALLOC_0 = 1;
	//EF_PROTECT_BELOW = 1;

	// Perform a quick version check (and exit the program if
	// the version flag is found). If we don't check now, there's no way for
	// the user to query the version with out having mpi up and running.
	for(int i = 1;i < argc;i++){
		if(strcmp(argv[i], "--version") == 0){
			cerr << argv[0] << " version "
			     << VERSION << endl;
			exit(1);
		}
	}
	if (argc < 6) {
		cerr << "mpiBLAST requires the following options: -d [database] -i [query file] -p [blast program name]\n";
		exit(1);
	}

#ifdef MPE
	char mpelogfmt[20] = "MPE_LOG_FORMAT=SLOG";
	putenv(mpelogfmt);
#endif
	double realstarttime = clock();

	/* Start up MPI */

	// Original version (no MPI thread support)
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &node_count);

	super_master_process = node_count - 1;

	double profile_time = MPI_Wtime();

	//Check for -p, -i, -d -- the required args
	int ArgsOK = 0;
	for(int i = 1;i < argc;i++){
		if(strcmp(argv[i], "-d") == 0)
			ArgsOK++;
		if(strcmp(argv[i], "-i") == 0)
			ArgsOK++;
		if(strcmp(argv[i], "-p") == 0)
			ArgsOK++;
		if ( (ArgsOK != 3) && (i == argc-1)){
			cerr << "mpiBLAST requires the following options: -d [database] -i [query file] -p [blast program name]\n";
			terminateProgram(-1);
		}
	}
	// set an error handler that won't abort
	// MPI_Errhandler_set( MPI_COMM_WORLD, MPI_ERRORS_RETURN );

	prog_start = MPI_Wtime();
	double progstarttime = clock();

	log_stream = &cerr;
	mpiblast_ptr = new MpiBlast();

	if(node_count < 2){
		//sorry, we need at least 2 processes.
		terminateProgram(0);
	}
#ifndef WIN32
    signal( SIGINT, terminateProgram );
    signal( SIGTERM, terminateProgram );
    signal( SIGSEGV, terminateProgram );
#endif

	// register an exit function to abort MPI if blast exits
	// without an error code
	atexit( &exitCleanly );

	// mpiblast_ptr = &mpb;

	// int retcode = mpb.main( argc, argv );
	int retcode = mpiblast_ptr->main( argc, argv );

	if( debug_msg ){
		LOG_MSG << "rank " << my_rank << " exiting successfully\n";
		LOG_MSG << "MPI startup  time is " << (double)((progstarttime - realstarttime)/CLOCKS_PER_SEC) << endl;
		LOG_MSG << "Copy time: " << copy_time << endl;
		LOG_MSG << "Search time: " << search_time << endl;
		LOG_MSG << "Results processing time: " << result_time << endl;
		LOG_MSG << "Write time: " << write_time << endl;
		LOG_MSG << "Run blast time: " << run_blast_time << endl;
		LOG_MSG << "Results gathering time: " << results_gather_time << endl;
	}

	// mpiblast_ptr = NULL;

	MPI_Barrier(MPI_COMM_WORLD);
	profile_time = MPI_Wtime() - profile_time;

	if(my_rank == super_master_process) {
		cerr<<"Total Execution Time: "<<profile_time<<endl;
	}

	try{
		if(collect_profile) {
			CollectSearchProfile(profile_filename, profile_time);
		}
		delete mpiblast_ptr;
	} catch(const char *error) {
		cerr << "Error collecting profiling info: " << error << endl;
	} catch(exception& e) {
		cerr << "Fatal exeption!" << endl;
	}

	MPI_Finalize();

	fin_intercept_lib();

	return retcode;
}


void deleteRegisteredFiles();
void deleteRegisteredFiles() {
	// don't be a slob, clean up after yourself:
	// delete any files that are laying around
	vector< string >& del_files = registerFileToDelete();
	for( int fileI = 0; fileI < del_files.size(); fileI++ )
		removeFile( del_files[ fileI ] );
	del_files.clear();	// clear the deleted files from the list
}

/**
 * This function is registered with atexit() to run before mpiBLAST terminates.
 * It deletes files that need to be deleted and calls terminateProgram() if
 * the NCBI library tried to abort execution
 */
void exitCleanly(){
	deleteRegisteredFiles();
	// abort MPI if blast crashed
	if( running_blast )
		terminateProgram( -1 );
}

/**
 * Aborts the mpiBLAST program using MPI_Abort()  On unix this function gets
 * registered to handle interrupt, terminate, and segfault signals.
 * @param sig	A unix signal to die with.  Special case is 0, indicating that mpiblast is being run with too few processes.
 */
void terminateProgram( int sig ){
	if(sig==0)
		cerr << "Sorry, mpiBLAST must be run on 3 or more nodes\n";
	else {
		cerr << my_rank << "\t" << (MPI_Wtime()-prog_start)
		     << "\tBailing out with signal "<< sig << endl;
		if( mpiblast_ptr != NULL )
			if( mpiblast_ptr->remove_db ){
				deleteMpiBlastFiles( localPath );
				removeDirectory( localPath );
			}
	}
	MPI_Abort( MPI_COMM_WORLD, 0 );
}

void CollectSearchProfile(const string& file_name, double total_time) {
	double max_search_time = 0;
	double min_search_time = FLT_MAX;
	double avg_search_time = 0;

	double max_result_time = 0;
	double min_result_time = FLT_MAX;
	double avg_result_time = 0;

	double max_process_output_time = 0;
	double min_process_output_time = FLT_MAX;
	double avg_process_output_time = 0;

	double max_run_blast_time = 0;
	double min_run_blast_time = FLT_MAX;
	double avg_run_blast_time = 0;

	double max_write_time = 0;
	double min_write_time = FLT_MAX;
	double avg_write_time = 0;

	double max_wait_write_time = 0;
	double min_wait_write_time = FLT_MAX;
	double avg_wait_write_time = 0;

	double max_wait_result_time = 0;
	double min_wait_result_time = FLT_MAX;
	double avg_wait_result_time = 0;

	double max_receive_message_time = 0;
	double min_receive_message_time = FLT_MAX;
	double avg_receive_message_time = 0;

	int max_peak_pending_offsets = 0;
	int min_peak_pending_offsets = INT_MAX;
	int avg_peak_pending_offsets = 0;

	if(group_rank == scheduler_process && my_rank != super_master_process) {
		// receive profiling data from writer master
		MPI_Status probe_status, status;

		if(scheduler_process != writer_process) {
			MPI_Probe(writer_process, PROFILE_TAG, group_comm, &probe_status);
			int msg_size;
			MPI_Get_count(&probe_status, MPI_BYTE, &msg_size);
			char* recv_buf = (char*)malloc(msg_size+1);
			CHECK_NULPTR(recv_buf);
			MPI_Recv(recv_buf, msg_size, MPI_BYTE, writer_process, PROFILE_TAG, group_comm, &status);

			istringstream is_recv;
			is_recv.str(recv_buf);
			free(recv_buf);

			is_recv >> output_info_time;
			is_recv >> writer_master_idle_time;
			is_recv >> writer_master_handle_time;
			is_recv >> write_time;
			is_recv >> receive_msg_time;
			is_recv >> receive_msg_size;
			is_recv >> max_num_write;
		}

		// receive profiling data from workers
		vector <double> copy_time_vec;
		vector <double> search_time_vec;
		vector <double> result_time_vec;
		vector <double> process_output_time_vec;
		vector <double> write_time_vec;
		vector <double> run_blast_time_vec;
		vector <double> idle_time_vec;
		vector <double> handle_time_vec;
		vector <double> wait_result_time_vec;
		vector <double> wait_write_time_vec;
		vector <double> receive_message_time_vec;
		vector <int> peak_pending_offsets_vec;

		for(int workerI=0; workerI<group_node_count; workerI++) {
			if(!IsWorker(workerI)) {
				continue;
			}

			MPI_Status probe_status, status;
			MPI_Probe(workerI, PROFILE_TAG, group_comm, &probe_status);
			int msg_size;
			MPI_Get_count(&probe_status, MPI_BYTE, &msg_size);
			char* recv_buf = (char*)malloc(msg_size+1);
			CHECK_NULPTR(recv_buf);
			MPI_Recv(recv_buf, msg_size, MPI_BYTE, workerI, PROFILE_TAG, group_comm, &status);
			recv_buf[msg_size] = 0;
			double wrk_copy_time;
			double wrk_search_time;
			double wrk_result_time;
			double wrk_process_output_time;
			double wrk_write_time;
			double wrk_run_blast_time;
			double wrk_idle_time;
			double wrk_writer_handle_time;
			double wrk_wait_result_time;
			double wrk_wait_write_time;
			double wrk_receive_message_time;
			int wrk_peak_pending_offsets;

			istringstream is_recv;
			is_recv.str(recv_buf);
			free(recv_buf);

			is_recv >> wrk_copy_time;
			is_recv >> wrk_search_time;
			is_recv >> wrk_result_time;
			is_recv >> wrk_process_output_time;
			is_recv >> wrk_write_time;
			is_recv >> wrk_run_blast_time;
			is_recv >> wrk_idle_time;
			is_recv >> wrk_writer_handle_time;
			is_recv >> wrk_wait_result_time;
			is_recv >> wrk_wait_write_time;
			is_recv >> wrk_receive_message_time;
			is_recv >> wrk_peak_pending_offsets;

			if(wrk_search_time > max_search_time) {
				max_search_time = wrk_search_time;
			}
			if(wrk_search_time < min_search_time) {
				min_search_time = wrk_search_time;
			}
			avg_search_time += wrk_search_time;

			if(wrk_result_time > max_result_time) {
				max_result_time = wrk_result_time;
			}
			if(wrk_result_time < min_result_time) {
				min_result_time = wrk_result_time;
			}
			avg_result_time += wrk_result_time;
			if(wrk_process_output_time > max_process_output_time) {
				max_process_output_time = wrk_process_output_time;
			}
			if(wrk_process_output_time < min_process_output_time) {
				min_process_output_time = wrk_process_output_time;
			}
			avg_process_output_time += wrk_process_output_time;


			if(wrk_run_blast_time > max_run_blast_time) {
				max_run_blast_time = wrk_run_blast_time;
			}
			if(wrk_run_blast_time < min_run_blast_time) {
				min_run_blast_time = wrk_run_blast_time;
			}
			avg_run_blast_time += wrk_run_blast_time;

			if(wrk_write_time > max_write_time) {
				max_write_time = wrk_write_time;
			}
			if(wrk_write_time < min_write_time) {
				min_write_time = wrk_write_time;
			}
			avg_write_time += wrk_write_time;

			if(wrk_wait_write_time > max_wait_write_time) {
				max_wait_write_time = wrk_wait_write_time;
			}
			if(wrk_wait_write_time < min_wait_write_time) {
				min_wait_write_time = wrk_wait_write_time;
			}
			avg_wait_write_time += wrk_wait_write_time;

			if(wrk_wait_result_time > max_wait_result_time) {
				max_wait_result_time = wrk_wait_result_time;
			}
			if(wrk_wait_result_time < min_wait_result_time) {
				min_wait_result_time = wrk_wait_result_time;
			}
			avg_wait_result_time += wrk_wait_result_time;

			if(wrk_receive_message_time > max_receive_message_time) {
				max_receive_message_time = wrk_receive_message_time;
			}
			if(wrk_receive_message_time < min_receive_message_time) {
				min_receive_message_time = wrk_receive_message_time;
			}
			avg_receive_message_time += wrk_receive_message_time;

			if(wrk_peak_pending_offsets > max_peak_pending_offsets) {
				max_peak_pending_offsets = wrk_peak_pending_offsets;
			}
			if(wrk_peak_pending_offsets < min_peak_pending_offsets) {
				min_peak_pending_offsets = wrk_peak_pending_offsets;
			}
			avg_peak_pending_offsets += wrk_peak_pending_offsets;

			copy_time_vec.push_back(wrk_copy_time);
			search_time_vec.push_back(wrk_search_time);
			result_time_vec.push_back(wrk_result_time);
			process_output_time_vec.push_back(wrk_process_output_time);
			write_time_vec.push_back(wrk_write_time);
			run_blast_time_vec.push_back(wrk_run_blast_time);
			idle_time_vec.push_back(wrk_idle_time);
			handle_time_vec.push_back(wrk_writer_handle_time);
			wait_result_time_vec.push_back(wrk_wait_result_time);
			wait_write_time_vec.push_back(wrk_wait_write_time);
			receive_message_time_vec.push_back(wrk_receive_message_time);
			peak_pending_offsets_vec.push_back(wrk_peak_pending_offsets);
		}

		avg_search_time = avg_search_time / search_time_vec.size();
		avg_result_time = avg_result_time / result_time_vec.size();
		avg_process_output_time = avg_process_output_time / process_output_time_vec.size();
		avg_run_blast_time = avg_run_blast_time / run_blast_time_vec.size();
		avg_write_time = avg_write_time / write_time_vec.size();
		avg_wait_write_time = avg_wait_write_time / wait_write_time_vec.size();
		avg_wait_result_time = avg_wait_result_time / wait_result_time_vec.size();
		avg_receive_message_time = avg_receive_message_time / receive_message_time_vec.size();
		avg_peak_pending_offsets = avg_peak_pending_offsets / peak_pending_offsets_vec.size();

		int group_id = GroupManager::Instance()->GetGroupId(my_rank);
		ostringstream os;
		os << file_name << "." << group_id;
		ofstream ofs(os.str().c_str());
		if(!ofs.is_open()) {
			throw __FILE__ "cannot open profile output file";
		}

		ofs << "worker profile: copy_time, search_time, result_time, process_output_time, write_time, run_blast_time, idle_time, writer_worker_handle_time, wait_result_time, wait_write_time, receive_message_time, peak_pending_offsets" << endl;
		for(int idx=0; idx<copy_time_vec.size(); idx++) {
			ofs << copy_time_vec[idx] << ", " << search_time_vec[idx] << ", " << result_time_vec[idx] << ", " << process_output_time_vec[idx] << ", " << write_time_vec[idx] << ", " << run_blast_time_vec[idx] << ", " << idle_time_vec[idx] << ", " << handle_time_vec[idx] << ", " << wait_result_time_vec[idx] << ", " << wait_write_time_vec[idx] << ", " <<  receive_message_time_vec[idx] << peak_pending_offsets_vec[idx] << endl;
		}

		ofs << endl;

		// statistics of scheduler
		ofs << "master broadcast time: " << query_broadcast_time << endl;
		ofs << "scheduler idle time: " << scheduler_idle_time << endl;
		ofs << "scheduler handle time: " << scheduler_handle_time << endl;
		ofs << "aquire query segments time: " << acquire_segement_time << endl;
		ofs << "load queries time: " << load_queries_time << endl;
		ofs << "number of scheduler execution: " << scheduler_called << endl;
		ofs << endl;

		// statistics of writer master
		ofs << "get output info time: " << output_info_time << endl;
		ofs << "writer master idle time: " << writer_master_idle_time << endl;
		ofs << "writer master handle time: " << writer_master_handle_time << endl;
		ofs << "master write time: " << write_time << endl;
		ofs << "master receive message time: " << receive_msg_time << endl;
		ofs << "master receive message size (bytes): " << receive_msg_size << endl;
		ofs << "max concurrent writing: " << max_num_write << endl;
		ofs << "default output time: " << default_output_time << endl;
		ofs << endl;

		// statistics for workers
		ofs << "max worker search time: " << max_search_time << endl;
		ofs << "min worker search time: " << min_search_time << endl;
		ofs << "avg worker search time: " << avg_search_time << endl;
		ofs << "max worker result time: " << max_result_time << endl;
		ofs << "min worker result time: " << min_result_time << endl;
		ofs << "avg worker result time: " << avg_result_time << endl;
		ofs << "max worker process_output time: " << max_process_output_time << endl;
		ofs << "min worker process_output time: " << min_process_output_time << endl;
		ofs << "avg worker process_output time: " << avg_process_output_time << endl;
		ofs << "max worker run_blast time: " << max_run_blast_time << endl;
		ofs << "min worker run_blast time: " << min_run_blast_time << endl;
		ofs << "avg worker run_blast time: " << avg_run_blast_time << endl;
		ofs << "max worker write time: " << max_write_time << endl;
		ofs << "min worker write time: " << min_write_time << endl;
		ofs << "avg worker write time: " << avg_write_time << endl;
		ofs << "max worker wait result time: " << max_wait_result_time << endl;
		ofs << "min worker wait result time: " << min_wait_result_time << endl;
		ofs << "avg worker wait result time: " << avg_wait_result_time << endl;
		ofs << "max worker wait write time: " << max_wait_write_time << endl;
		ofs << "min worker wait write time: " << min_wait_write_time << endl;
		ofs << "avg worker wait write time: " << avg_wait_write_time << endl;
		ofs << "max worker receive message time: " << max_receive_message_time << endl;
		ofs << "min worker receive message time: " << min_receive_message_time << endl;
		ofs << "avg worker receive message time: " << avg_receive_message_time << endl;
		ofs << "max peak pending offsets: " << max_peak_pending_offsets << endl;
		ofs << "min peak pending offsets: " << min_peak_pending_offsets << endl;
		ofs << "avg peak pending offsets: " << avg_peak_pending_offsets << endl;
		ofs.close();

	} else if (group_rank == writer_process && my_rank != super_master_process) {
		ostringstream os_send;
		os_send << output_info_time << " ";
		os_send << writer_master_idle_time << " ";
		os_send << writer_master_handle_time << " ";
		os_send << write_time << " ";
		os_send << receive_msg_time << " ";
		os_send << receive_msg_size << " ";
		os_send << max_num_write << " ";

		string send_buf = os_send.str();
		MPI_Send((void*)(send_buf.c_str()), send_buf.length(), MPI_BYTE, scheduler_process, PROFILE_TAG, group_comm);
	} else if (my_rank != super_master_process) {
		ostringstream os_send;

		os_send << copy_time << " ";
		os_send << search_time << " ";
		os_send << result_time << " ";
		os_send << process_output_time << " ";
		os_send << write_time << " ";
		os_send << run_blast_time << " ";
		os_send << idle_time << " ";
		os_send << writer_worker_handle_time << " ";
		os_send << wait_result_time << " ";
		os_send << wait_write_time << " ";
		os_send << receive_msg_time << " ";
		os_send << peak_pending_offsets << " ";

		string send_buf = os_send.str();
		MPI_Send((void*)(send_buf.c_str()), send_buf.length(), MPI_BYTE, scheduler_process, PROFILE_TAG, group_comm);
	}

	// collect statistics from groups
	if(my_rank == super_master_process) {

		int num_groups = GroupManager::Instance()->GetNumGroups();

		vector <double> query_broadcast_time_vec(num_groups, 0);
		vector <double> scheduler_handle_time_vec(num_groups, 0);
		vector <double> acquire_segment_time_vec(num_groups, 0);
		vector <double> load_queries_time_vec(num_groups, 0);
		vector <double> writer_master_handle_time(num_groups, 0);
		vector <double> master_write_time_vec(num_groups, 0);

		vector <double> max_search_time_vec(num_groups, 0);
		vector <double> min_search_time_vec(num_groups, 0);
		vector <double> avg_search_time_vec(num_groups, 0);
		vector <double> max_result_time_vec(num_groups, 0);
		vector <double> min_result_time_vec(num_groups, 0);
		vector <double> avg_result_time_vec(num_groups, 0);
		vector <double> max_run_blast_time_vec(num_groups, 0);
		vector <double> min_run_blast_time_vec(num_groups, 0);
		vector <double> avg_run_blast_time_vec(num_groups, 0);
		vector <double> max_write_time_vec(num_groups, 0);
		vector <double> min_write_time_vec(num_groups, 0);
		vector <double> avg_write_time_vec(num_groups, 0);

		MPI_Status probe_status, status;
		for(int i=0; i<num_groups; i++) {
			MPI_Probe(MPI_ANY_SOURCE, PROFILE_TAG, MPI_COMM_WORLD, &probe_status);
			int msg_size;
			MPI_Get_count(&probe_status, MPI_BYTE, &msg_size);
			char* recv_buf = (char*)malloc(msg_size+1);
			CHECK_NULPTR(recv_buf);
			MPI_Recv(recv_buf, msg_size, MPI_BYTE, probe_status.MPI_SOURCE, PROFILE_TAG, MPI_COMM_WORLD, &status);
			recv_buf[msg_size] = 0;

			istringstream is_recv;
			is_recv.str(recv_buf);
			free(recv_buf);

			int group_id = -1;
			is_recv >> group_id;

			double grp_query_broadcast_time = 0;
			double grp_scheduler_handle_time = 0;
			double grp_acquire_segment_time = 0;
			double grp_load_queries_time = 0;
			double grp_writer_master_handle_time = 0;
			double grp_master_write_time = 0;

			double grp_max_search_time = 0;
			double grp_min_search_time = 0;
			double grp_avg_search_time = 0;
			double grp_max_result_time = 0;
			double grp_min_result_time = 0;
			double grp_avg_result_time = 0;
			double grp_max_run_blast_time = 0;
			double grp_min_run_blast_time = 0;
			double grp_avg_run_blast_time = 0;
			double grp_max_write_time = 0;
			double grp_min_write_time = 0;
			double grp_avg_write_time = 0;

			is_recv >> grp_query_broadcast_time;
			is_recv >> grp_scheduler_handle_time;
			is_recv >> grp_acquire_segment_time;
			is_recv >> grp_load_queries_time;
			is_recv >> grp_writer_master_handle_time;
			is_recv >> grp_master_write_time;

			is_recv >> grp_max_search_time;
			is_recv >> grp_min_search_time;
			is_recv >> grp_avg_search_time;
			is_recv >> grp_max_result_time;
			is_recv >> grp_min_result_time;
			is_recv >> grp_avg_result_time;
			is_recv >> grp_max_run_blast_time;
			is_recv >> grp_min_run_blast_time;
			is_recv >> grp_avg_run_blast_time;
			is_recv >> grp_max_write_time;
			is_recv >> grp_min_write_time;
			is_recv >> grp_avg_write_time;

			query_broadcast_time_vec[group_id] = grp_query_broadcast_time;
			scheduler_handle_time_vec[group_id] = grp_scheduler_handle_time;
			acquire_segment_time_vec[group_id] = grp_acquire_segment_time;
			load_queries_time_vec[group_id] = load_queries_time;
			writer_master_handle_time[group_id] = grp_writer_master_handle_time;
			master_write_time_vec[group_id] = grp_master_write_time;

			max_search_time_vec[group_id] = grp_max_search_time;
			min_search_time_vec[group_id] = grp_min_search_time;
			avg_search_time_vec[group_id] = grp_avg_search_time;
			max_result_time_vec[group_id] = grp_max_result_time;
			min_result_time_vec[group_id] = grp_min_result_time;
			avg_result_time_vec[group_id] = grp_avg_result_time;
			max_run_blast_time_vec[group_id] = grp_max_run_blast_time;
			min_run_blast_time_vec[group_id] = grp_min_run_blast_time;
			avg_run_blast_time_vec[group_id] = grp_avg_run_blast_time;
			max_write_time_vec[group_id] = grp_max_write_time;
			min_write_time_vec[group_id] = grp_min_write_time;
			avg_write_time_vec[group_id] = grp_avg_write_time;
		}

		ofstream ofs(file_name.c_str());
		if(!ofs.is_open()) {
			throw __FILE__ "cannot open profile output file";
		}
		ofs << "group profile: query_broadcast, scheduler_handle, acquire_segment_time, load_queries_time, writer_master_handle, master_write, max_search, min_search, avg_search, max_result, min_result, avg_result, max_run_blast, min_run_blast, avg_run_blast, max_write, min_write, avg_write" << endl;
		for(int group_id=0; group_id<num_groups; group_id++) {
			ofs << query_broadcast_time_vec[group_id] << ", "
				<< scheduler_handle_time_vec[group_id] << ", "
				<< acquire_segment_time_vec[group_id] << ", "
				<< load_queries_time_vec[group_id] << ", "
				<< writer_master_handle_time[group_id] << ", "
				<< master_write_time_vec[group_id] << ", "
				<< max_search_time_vec[group_id] << ", "
				<< min_search_time_vec[group_id] << ", "
				<< avg_search_time_vec[group_id] << ", "
				<< max_result_time_vec[group_id] << ", "
				<< min_result_time_vec[group_id] << ", "
				<< avg_result_time_vec[group_id] << ", "
				<< max_run_blast_time_vec[group_id] << ", "
				<< min_run_blast_time_vec[group_id] << ", "
				<< avg_run_blast_time_vec[group_id] << ", "
				<< max_write_time_vec[group_id] << ", "
				<< min_write_time_vec[group_id] << ", "
				<< avg_write_time_vec[group_id] << endl;
		}

		ofs << "system initialization time: " << system_init_time << endl;
		ofs << "global query broadcast time: " << query_broadcast_time << endl;
		ofs << "db pre-distribution time: " << db_distribute_time << endl;
		ofs << "overall execution time: " << total_time << endl;
		ofs << endl;
		ofs.close();

	} else if (group_rank == 0) {

		ostringstream os_send;
		int group_id = GroupManager::Instance()->GetGroupId(my_rank);

		os_send << group_id << " ";

		os_send << query_broadcast_time << " ";
		os_send << scheduler_handle_time << " ";
		os_send << acquire_segement_time << " ";
		os_send << load_queries_time << " ";
		os_send << writer_master_handle_time << " ";
		os_send << write_time << " "; // master write time

		os_send << max_search_time << " ";
		os_send << min_search_time << " ";
		os_send << avg_search_time << " ";
		os_send << max_result_time << " ";
		os_send << min_result_time << " ";
		os_send << avg_result_time << " ";
		os_send << max_run_blast_time << " ";
		os_send << min_run_blast_time << " ";
		os_send << avg_run_blast_time << " ";
		os_send << max_write_time << " ";
		os_send << min_write_time << " ";
		os_send << avg_write_time << " ";

		string send_buf = os_send.str();
		MPI_Send((void*)(send_buf.c_str()), send_buf.length(), MPI_BYTE, super_master_process, PROFILE_TAG, MPI_COMM_WORLD);
	}
}

MpiBlast::~MpiBlast(){
	log_stream = &cerr;
	if(myscheduler!=NULL) delete myscheduler;
	if(writer_m!=NULL) delete writer_m;
	if(writer_w!=NULL) delete writer_w;
	if(ncbiblast!=NULL) delete ncbiblast;
}


/**
 *  Scans through the list of SeqAlignPtrs separating it into several lists,
 *  each of which corresponds to a single BLAST database hit.
 */
void BreakUpResults( vector< pair< int, SeqAlignPtr > >& results_vec, SeqAlignPtr sap, int searched_frag_id ){
	SeqAlignPtr cur_a = sap;
	results_vec.clear();
	results_vec.push_back( pair< int, SeqAlignPtr >( searched_frag_id, cur_a ) );
	int cur_count = 1;
	//	cerr << cur_count << "\t" << "cur_a: " << cur_a << endl;
	while( cur_a->next != NULL ){
		//		cerr << cur_count << "\t" << "cur_a: " << cur_a->next << endl;
		SeqIdPtr cur_id = TxGetSubjectIdFromSeqAlign( cur_a );
		SeqIdPtr next_id = TxGetSubjectIdFromSeqAlign( cur_a->next );
		if( SeqIdComp( cur_id, next_id ) != SIC_YES ){
			// break the chain at cur_a, add cur_a->next
			SeqAlignPtr tmp_next = cur_a->next;
			cur_a->next = NULL;
			//			cerr << "tmp_next\t" << tmp_next << endl;
			results_vec.push_back( pair< int, SeqAlignPtr >( searched_frag_id, tmp_next ) );
			cur_a = tmp_next;
			//			cerr << "Segment count is: " << cur_count << endl;
			cur_count = 1;
		}else{
			cur_a = cur_a->next;
			cur_count++;
		}
	}
}


/**
 * Mergesorts new results with existing results for a single query
 */
void mergeSeqAlignResults( const vector< pair< int, SeqAlignPtr > >& results, const vector< pair< int, SeqAlignPtr > >& new_results,
			  vector< pair< int, SeqAlignPtr > >& merged_results ){
	uint resultI = 0;
	uint resultJ = 0;
	while( resultI != results.size() && resultJ != new_results.size() ){

		int a_score = 0, b_score = 0, a_number, b_number;
		Nlm_FloatHi a_bit_score = 0, b_bit_score = 0, a_evalue, b_evalue;

		SeqAlignPtr sapI = results[ resultI ].second;
		SeqAlignPtr sapJ = new_results[ resultJ ].second;

		GetScoreAndEvalue( sapI, &a_score, &a_bit_score, &a_evalue, &a_number );
		GetScoreAndEvalue( sapJ, &b_score, &b_bit_score, &b_evalue, &b_number );
//		cerr << "a_score: " << a_bit_score << "\t b_score: " << b_bit_score << endl;
		if( a_evalue < b_evalue ||
			a_evalue == b_evalue && a_bit_score > b_bit_score ||
			a_evalue == b_evalue && a_bit_score == b_bit_score && results[ resultI ].first > new_results[ resultJ ].first)
		{
			merged_results.push_back( results[ resultI ] );
			resultI++;
		}else{
			merged_results.push_back( new_results[ resultJ ] );
			resultJ++;
		}
	}
	for( ; resultI != results.size(); resultI++ ){
		merged_results.push_back( results[ resultI ] );
	}
	for( ; resultJ != new_results.size(); resultJ++ ){
		merged_results.push_back( new_results[ resultJ ] );
	}
}


void mergeResults( map< int, vector< pair< int, SeqAlignPtr > > >& query_map, vector< SeqAlignPtr >& new_results, int searched_frag_id, int* result_list )
{
	for( uint resultI = 0; resultI < new_results.size(); resultI++ ){

		int query_id = result_list[ resultI ];
		vector< pair< int, SeqAlignPtr > > cur_results;
		BreakUpResults( cur_results, new_results[ resultI ], searched_frag_id );

		map< int, vector< pair< int, SeqAlignPtr > > >::iterator query_iter = query_map.find( query_id );

		if( query_iter == query_map.end() ){
			// no results for this query yet.  add it.
			query_map.insert( map< int, vector< pair< int, SeqAlignPtr > > >::value_type( query_id, cur_results ) );
		}else{
			// merge these results with the existing query results.
			vector< pair< int, SeqAlignPtr > > merged_results;
			mergeSeqAlignResults( query_iter->second, cur_results, merged_results );
			query_map.erase( query_iter );
			query_map.insert( map< int, vector< pair< int, SeqAlignPtr > > >::value_type( query_id, merged_results ) );
		}
	}
}

// JA -- send fragments over MPI
bool MpiBlast::sendFragsViaMPI(int fragment_id, int worker_id){

	const vector<string>& extensions = FragmentExtensions(db_type);

	char fragid[8];
	memset( fragid, 0, 8 );
	sprintf( fragid, "%03d", fragment_id );

	string copy_fragment_name;
	string copy_src;

	vector< string > copy_files;
	vector< string > ext;

	// determine which frags to copy and how many
	for( uint extI = 0; extI < extensions.size(); extI++ ){
		copy_fragment_name = database_name + "." + fragid + extensions[ extI ];
		copy_src = config.sharedPath() + copy_fragment_name;
		if( doesFileExist( copy_src ) ){
			ext.push_back( extensions[ extI ] );
			copy_files.push_back( copy_src );
		}
	}

	// tell the worker how many files it's getting
	int file_count = copy_files.size();
	MPI_Send( &file_count, 1, MPI_INT, worker_id, DB_DATATOSEND_TAG, MPI_COMM_WORLD );

	// send each filename extension and file
	for( uint fileI = 0; fileI < copy_files.size(); fileI++ ){
		char ext_name[5];
		strncpy( ext_name, ext[ fileI ].c_str(), 5 );
		MPI_Send( ext_name, ext[ fileI ].length() + 1, MPI_CHAR, worker_id, DB_EXTENSION_DATA_TAG, MPI_COMM_WORLD );
		sendFile( copy_files[ fileI ], worker_id, DB_DATATOSEND_TAG );
	}

	return true;
}


// JA -- receive fragments over MPI
bool MpiBlast::recvFragsViaMPI( int fragment_id ){

	char fragid[8];
	memset( fragid, 0, 8 );
	sprintf( fragid, "%03d", fragment_id );

	string copy_fragment_name;
	string copy_dest;

	if( debug_msg )	LOG_MSG << "Requesting fragment " << fragment_id << " from rank 0\n";
	// tell rank 0 which fragment we want to copy
	int f_tmp = fragment_id;
	int event_array[2];
	event_array[0] = DB_DATATOSEND_TAG;
	event_array[1] = int_size;
	MPI_Send( &event_array, 2, MPI_INT, super_master_process, EVENT_TYPE, MPI_COMM_WORLD );
	MPI_Send( &f_tmp, 1, MPI_INT, super_master_process, DB_DATATOSEND_TAG, MPI_COMM_WORLD );

	if( debug_msg )	LOG_MSG << "Getting frag file count for " << fragment_id << " from rank 0\n";
	// get file count from rank 0
	int file_count;
	MPI_Recv( &file_count, 1, MPI_INT, super_master_process, DB_DATATOSEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

	if( debug_msg )	LOG_MSG << "Got frag file count " << file_count << " from rank 0\n";

	for( uint fileI = 0; fileI < file_count; fileI++ ){
		char ext[5];
		MPI_Recv( ext, 5, MPI_CHAR, super_master_process, DB_EXTENSION_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		if( debug_msg )	LOG_MSG << "Got frag extension name " << ext << " from rank 0\n";
		copy_fragment_name = database_name + "." + fragid;
		copy_fragment_name += ext;
		copy_dest = config.localPath() + copy_fragment_name;
		recvFile( copy_dest, super_master_process, DB_DATATOSEND_TAG );
	}

	return true;
}

void MpiBlast::super_master() {
	// scheduling query segments to different groups
	// limit the concurrent copy activity

	// while loop waiting for query segment request
	// distribute query batch in form of <start_query, end_query>, where queries in range [start_query, end_query-1] will be searched
	int curr_query_id = 0;
	int curr_segment_id = 0;
	int complete_groups = 0;
	int event_array[2];

	queue< int > db_access_queue; // JA DB access queue of workers waiting to access the database

	int num_copying_db = 0;
	int num_groups = GroupManager::Instance()->GetNumGroups();
	int num_first_assign = CalculateFirstAssignment();

	curr_segment_id = num_groups;
	curr_query_id = num_first_assign * num_groups;
	if(curr_query_id > query_count) {
		curr_query_id = query_count;
	}

	IsendMgr ism;

	while(true) {

		ism.CheckPendingSends();

		MPI_Status probe_status, status;

		MPI_Probe( MPI_ANY_SOURCE, EVENT_TYPE, MPI_COMM_WORLD, &probe_status );
		MPI_Recv(event_array, 2, MPI_INT, probe_status.MPI_SOURCE, EVENT_TYPE, MPI_COMM_WORLD, &status);

		int master_rank = probe_status.MPI_SOURCE;
		int group_id = GroupManager::Instance()->GetGroupId(master_rank);

		if(debug_msg) {
			LOG_MSG << "Supermaster got tag " << event_array[0] << " from " << probe_status.MPI_SOURCE << endl;
		}

		if( event_array[0] == QUERY_SEGMENT_REQUEST ) {
			int assign_array[5];
			memset(assign_array, 0, 5);
			if(curr_query_id < query_count) {
				assign_array[0] = SEARCH_SEGMENT;
				assign_array[1] = curr_segment_id;
				assign_array[2] = curr_query_id;

				int end_query = 0;
				int curr_assign = query_segment_size;

				end_query = curr_query_id + curr_assign;
				if(end_query > query_count) {
					end_query = query_count;
				}
				assign_array[3] = end_query;

				if(assign_array[3] - assign_array[2] <= 0) {
					throw __FILE__ " Empty assignment!";
				}

				curr_query_id = end_query;
				curr_segment_id++;
			} else {
				// tell group all queries have finished
				assign_array[0] = NO_MORE_SEGMENT;
			}

			CommSendStruct* cssp = new CommSendStruct(5 * int_size);
			cssp->AddData(&assign_array[0], 5 * int_size);
			cssp->IsendData(MPI_COMM_WORLD, probe_status.MPI_SOURCE, ASSIGNMENT_TYPE);
			ism.AddPendingSend(cssp);

			if(debug_msg) {
				if(assign_array[0] != NO_MORE_SEGMENT) {
					LOG_MSG << "Give " << probe_status.MPI_SOURCE << " segment " << assign_array[1] << ", start_query " << assign_array[2] << ", end_query " << assign_array[3] << endl;
				} else {
					LOG_MSG << "No more query segments to assign" << endl;
				}
			}

		} else if (event_array[0] == GROUP_FINISHED) {
			if(++complete_groups == num_groups) {
				break;
			}
		} else if (event_array[0] == DB_DATATOSEND_TAG) {
			if( debug_msg ) LOG_MSG << "Worker " << probe_status.MPI_SOURCE << " requested fragment ";
			int fragment_id;
			MPI_Recv( &fragment_id, 1, MPI_INT, probe_status.MPI_SOURCE, DB_DATATOSEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			if( debug_msg ) CONT_LOG_MSG << fragment_id << " for MPI-copy" << endl;
			if (!sendFragsViaMPI( fragment_id, probe_status.MPI_SOURCE ))
				LOG_MSG << "***ERROR_2*** Could not send fragment " << fragment_id << " to Worker " << probe_status.MPI_SOURCE << endl;
		} else if (event_array[0] == REQUEST_FRAGMENT_COPY) {
			if (num_copying_db >= concurrent_accesses && copy_via != COPY_VIA_NONE){
				// JA -- put on Queue and wait turn for DB access
				db_access_queue.push( probe_status.MPI_SOURCE);
				if (debug_msg)
					LOG_MSG << "Node " << probe_status.MPI_SOURCE << " Waiting turn for DB access" << endl;
			} else {
				num_copying_db++;
				// JA -- send message so as to start copying
				if (debug_msg)
					LOG_MSG << "Node " << probe_status.MPI_SOURCE << " turn for DB access is now" << endl;

				int copy_command = CONFIRM_COPY_FRAGMENT;
				MPI_Send(&copy_command, 1, MPI_INT, probe_status.MPI_SOURCE, ASSIGNMENT_TYPE, MPI_COMM_WORLD);
			}
		} else if (event_array[0] == FRAGMENT_COPY_COMPLETE) {
			if(!db_access_queue.empty()) {
				if (debug_msg){
					LOG_MSG << "Node " << db_access_queue.front() << " has turn for DB access now" << endl;
				}

				int worker_to_copy = db_access_queue.front();

				int copy_command = CONFIRM_COPY_FRAGMENT;
				MPI_Send(&copy_command, 1, MPI_INT, worker_to_copy, ASSIGNMENT_TYPE, MPI_COMM_WORLD);

				db_access_queue.pop();
			} else {
				num_copying_db--;
			}
		}
	}

	ism.WaitPendingSends();
}

void MpiBlast::master() {

	myscheduler = new Scheduler(total_frags, group_node_count, query_count, concurrent_accesses, copy_via);

	// set up the first assignment
	int num_first_assign = CalculateFirstAssignment();
	int group_id = GroupManager::Instance()->GetGroupId(my_rank);
	int start_query = num_first_assign * group_id;
	int end_query = num_first_assign * (group_id + 1);
	if(end_query > query_count) {
		end_query = query_count;
	}

	addOpt( writer_opts, 'o', bitbucket );

	running_blast = true;
	ncbiblast->InitNCBI( writer_opts );
	ncbiblast->Init(MASTER_MODE);
	running_blast = false;

	double load_time = MPI_Wtime();

	if(query_in_mem) {
		QueryM::Instance()->ParseQuerydata();
		VFM::Instance()->Erase(query_filename);
	} else {
		QueryM::Instance()->EnableAutoLoad();
	}

	if(debug_msg) {
		LOG_MSG << "Loading queries took " << MPI_Wtime() - load_time << " seconds." << endl;
	}

	if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_SPLIT) {
		writer_m = WriterMasterIndividual::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	} else if (output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		writer_m = WriterMasterStreamline::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	} else if (output_strategy == WORKER_COLLECTIVE) {
		writer_m = WriterMasterCollective::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	}

	myscheduler->SetWriterMaster(writer_m);
	myscheduler->AddAssignment(group_id, start_query, end_query);

	if(use_brief_report) {
		if(debug_msg) {
			LOG_MSG << "Start preparing default output report" << endl;
		}
		double track_time = MPI_Wtime();
		writer_m->PrepareDefaultReport();
		default_output_time = MPI_Wtime() - track_time;
		if(debug_msg) {
			LOG_MSG << "End preparing default output report" << endl;
		}
	}

	myscheduler->RecvFragmentLists();

	// add irecv for event
	for(int i=0; i<num_reserved_events; i++) {
		CommRecvStruct* crsp = new CommRecvStruct(2 * int_size);
		crsp->IrecvData(group_comm, MPI_ANY_SOURCE, EVENT_TYPE);
		_irm->AddPendingRecv(crsp);
	}

	int scheduler_tags = SCHEDULER_TAGS;
	int writer_tags = WRITER_TAGS;

	bool scheduler_end = false;
	bool writer_end = false;
	static double handle_track_time = 0;

    int pool_count = 0;
	while( true ){
		CommRecvStruct* crsp_finish;
		MPI_Status recv_status;

		writer_m->ProcessAsyncOutput();
        if(QueryM::Instance()->GetNumWorkingQueries() == 0) {
            scheduler_end = true;
        }

		if(scheduler_end && writer_end) {
			break;
		}

		if(!_irm->TestAnyRecv(crsp_finish, recv_status)) {
/*            pool_count++;
            if(pool_count > 100) {
                usleep(1000);
                pool_count = 0;
            }
*/
			continue;
		}

		// _irm->WaitAnyRecv(crsp_finish, recv_status);

		if(recv_status.MPI_TAG == EVENT_TYPE) {
			int event = -1;
			int next_msg_size = -1;

			crsp_finish->ExtractData(&event, int_size);
			crsp_finish->ExtractData(&next_msg_size, int_size);

			if(next_msg_size > 0) {
				CommRecvStruct* crsp_next = new CommRecvStruct(next_msg_size);
				crsp_next->IrecvData(group_comm, recv_status.MPI_SOURCE, event);
				_irm->AddPendingRecv(crsp_next);
			}

			// start receiving a new event
			crsp_finish->ResetDataBuf();
			crsp_finish->IrecvData(group_comm, MPI_ANY_SOURCE, EVENT_TYPE);
			_irm->AddPendingRecv(crsp_finish);
		} else if (recv_status.MPI_TAG & scheduler_tags) {
			handle_track_time = MPI_Wtime();
			myscheduler->HandleMessages(crsp_finish, recv_status);

			scheduler_handle_time += MPI_Wtime() - handle_track_time;
		} else if ( recv_status.MPI_TAG & writer_tags ) {
			handle_track_time = MPI_Wtime();
			if(writer_m->HandleMessages(crsp_finish, recv_status)) {
				writer_end = true;
			}

			writer_master_handle_time += MPI_Wtime() - handle_track_time;
		} else {
			// Catch any out of place messages
			LOG_MSG << "Error: Received msg " << recv_status.MPI_TAG << " from " << recv_status.MPI_SOURCE << endl;
			throw __FILE__ "master(): Illegal message in main loop!";
		}

	} // end while (true)

	writer_m->Finalize();
	myscheduler->Finalize();

	// tell super master we finished
	int event_array[2];
	memset(event_array, 0 , 2);
	event_array[0] = GROUP_FINISHED;
	event_array[1] = -1;
	MPI_Send(event_array, 2, MPI_INT, super_master_process, EVENT_TYPE, MPI_COMM_WORLD);

	if(_irm->GetNumPendingRecvs() != num_reserved_events) {
		throw __FILE__ " MpiBlast::master - there are unprocessed pending messages. ";
	}

	delete _irm;
}

void MpiBlast::worker(){

	double track_time;

	meta_MPE_Log_event(worker_file_setup_start, 0,"setting up worker files");
	//
	// create blast alias file
	//
	string alias_basename = config.localPath() + database_name + "XXXXXX";
	string alias_filename;

	if(!use_virtual_frags) {
		getTempFileName( alias_basename );
		alias_filename = alias_basename + "." + db_type + "al";

		// do we really need this?
		moveFile( alias_basename, alias_filename );
		// be sure the alias file gets deleted before exiting
		registerFileToDelete( alias_filename );
	}

	stringstream rank_str;
	rank_str << my_rank;

	meta_MPE_Log_event(worker_file_setup_end, 0,"done setting up worker files");

	//
	// initialize the NCBI Toolbox, create blast command line
	//
	meta_MPE_Log_event(init_ncbi_start, 0,"start initNCBI()");

	addOpt( worker_opts, 'i', query_filename.c_str() );

	/** are we using a gi filter? */
	if ( filter_filename != "" )
	        addOpt( worker_opts, 'l', filter_filename.c_str() );

//	if( use_binary_asn )
//		addOpt( worker_opts, 'm', "11" );	/** ask for binary SeqAlign output */
//	else
//		addOpt( worker_opts, 'm', "10" );	/** ask for text SeqAlign output */
	addOpt( worker_opts, 'o', bitbucket );	/** since results get written to memory we can use /dev/null for any other output */
	addOpt( worker_opts, 'd', alias_basename.c_str() );

	if(num_threads > 0) {
		ostringstream os;
		os << num_threads;
		addOpt( worker_opts, 'a', os.str().c_str() );
	}

	// Init the NCBI library with the contrived blast command line
	running_blast = true;
	ncbiblast->InitNCBI( worker_opts );
	running_blast = false;

	meta_MPE_Log_event(init_ncbi_end, 0,"end initNCBI()");

	// construct writer
	// init blast to get the max show size
	ncbiblast->Init(WORKER_MODE);

	if(query_in_mem) {
		QueryM::Instance()->ParseQuerydata();
		VFM::Instance()->Erase(query_filename);
	}

	if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_SPLIT) {
		writer_w = WriterWorkerIndividual::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	} else if (output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		writer_w = WriterWorkerStreamline::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	} else if (output_strategy == WORKER_COLLECTIVE) {
		writer_w = WriterWorkerCollective::Instance(ncbiblast, query_count, total_frags, output_filename, 1);
	}

	ncbiblast->Cleanup();

	MPI_Status status;

	// These variables manage the database fragment list (for sending, receiving and writing
	// to disk).
	vector<int> fragment_list;
	int frag_id;	/**< The fragment currently being searched */

	// construct the fragment list file name
	string fragment_filename = config.localPath() + database_name + FRAG_LIST_EXTENSION;

	int query_id;	/**< The query that is currently being worked on */
	int segment_id;

	//
	// tell the scheduler node what fragments this node has
	//
	meta_MPE_Log_event(send_frag_list_start, 0,"begin send frag list to scheduler");

	if(use_virtual_frags) {
		frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type, true);
	} else {
		frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type );
	}
	//frag_list_file.setTimestamps( all_dates, total_datelen );
	frag_list_file.setTimestamps( ncbiblast->GetAllDates(), ncbiblast->GetTotalDateLen() );
	frag_list_file.checkTimestamps();
	frag_list_file.SendList( scheduler_process );
	if( debug_msg ){
		LOG_MSG << "Fragment list sent." << endl;
	}
	meta_MPE_Log_event(send_frag_list_end, 0,"end send frag list to scheduler");

	// loop until the assignment is SEARCH_COMPLETE or WORKER_QUIT
	int* assign_array = NULL;	/**< will store the assignment */
	int event_array[2];

	bool send_idle_message = true;
	bool pending_idle_messag = false;

	while( true ){

		double idle_track_time = MPI_Wtime();

		if(send_idle_message) {
			//
			// tell scheduler we are idle
			//
			event_array[0] = WORKER_IDLE;
			event_array[1] = int_size;

			CommSendStruct* cssp_event = new CommSendStruct(2 * int_size);
			cssp_event->AddData(&event_array[0], 2 * int_size);
			cssp_event->IsendData(group_comm, scheduler_process, EVENT_TYPE);
			writer_w->AddCommSend(cssp_event);

			int idle_msg = WORKER_IDLE;

			CommSendStruct* cssp_data = new CommSendStruct(int_size);
			cssp_data->AddData(&idle_msg, int_size);
			cssp_data->IsendData(group_comm, scheduler_process, WORKER_IDLE);
			writer_w->AddCommSend(cssp_data);
			pending_idle_messag = true;

			if( debug_msg ){
					LOG_MSG << "Idle message sent" << endl;
			}
		}

		if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			writer_w->ProcessAsyncOutput();
		}
		writer_w->CheckCommSends();
		writer_w->CheckResultSends();

		// check arriving messages
		int array_size;
		memset( &status, 0, sizeof(MPI_Status) );
		// MPI_Probe( MPI_ANY_SOURCE, ASSIGNMENT_TYPE, group_comm, &status );

		int flag = 0;
		int pool_count = 0;
		while (!flag) {
			MPI_Iprobe( MPI_ANY_SOURCE, ASSIGNMENT_TYPE, group_comm, &flag, &status );
			//pool_count++;

			writer_w->CheckCommSends();
			writer_w->CheckResultSends();
			if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
				writer_w->ProcessAsyncOutput();
			}

/*			if(pool_count > 50) {
				usleep(1000);
				pool_count = 0;
			}*/
		}

		idle_time += MPI_Wtime() - idle_track_time;
		if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			writer_w->ProcessAsyncOutput();
		}
		writer_w->CheckCommSends();
		writer_w->CheckResultSends();

		MPI_Get_count( &status, MPI_INT, &array_size );
		assign_array = (int*)malloc(int_size*array_size);
		CHECK_NULPTR(assign_array);

		MPI_Status recv_status;
		MPI_Recv( assign_array, array_size, MPI_INT, status.MPI_SOURCE, ASSIGNMENT_TYPE, group_comm, &recv_status );

		if( debug_msg ){
    			LOG_MSG << "Assignment " << assign_array[0] << " is " << array_size << " ints long from rank " << status.MPI_SOURCE << endl;
   		}

		//
		// we have an assignment, now act on it!
		//
		// should not receive WORKER_QUIT here
		if( assign_array[0] == WORKER_QUIT ) {
			throw __FILE__ "MpiBlast::worker -- unexpected WORKER_QUIT tag";
		}

		// exit the worker loop if there's no more work to do
		if( assign_array[0] == SEARCH_COMPLETE ) {
			free(assign_array);
			assign_array = NULL;
			break;
		}

		send_idle_message = true;

		if( assign_array[0] == WRITER_EVENT) {
			double handle_time = MPI_Wtime();
			writer_w->HandleMessages(status, &(assign_array[1]));
			writer_worker_handle_time += MPI_Wtime() - handle_time;

			if(!pending_idle_messag && (writer_w->GetNumPendingOffsets() < max_pending_offsets)) {
				send_idle_message = true;
				if(debug_msg) {
					LOG_MSG << "The number of pending offsets dropped, start fetching new assignment" << endl;
				}
			} else {
				send_idle_message = false;
			}
		} else {
			pending_idle_messag = false;

			if (assign_array[0] != DO_NOTHING ) {
				query_id = assign_array[3];
				segment_id = assign_array[1];
				//end_query_id = assign_array[4];
				frag_id = assign_array[2];

				// JA -- If assignment is SEARCH_FRAGMENT, worker will immediately receive these messages
				//       If assignment is COPY_FRAGMENT, worker is put in queue and will receive messages when at front of queue

				if( assign_array[0] == COPY_FRAGMENT ) {

					track_time = MPI_Wtime();

					// first check whether another worker has already copied this
					// fragment to the local storage directory
					// note that the scheduler mustn't allow two workers to copy
					// the same fragment simultaneously otherwise they could overwrite
					// each other's files when multiple workers share local storage
					if(use_virtual_frags) {
						frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type, true);
					} else {
						frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type );
					}
					//frag_list_file.setTimestamps( all_dates, total_datelen );
					frag_list_file.setTimestamps( ncbiblast->GetAllDates(), ncbiblast->GetTotalDateLen() );
					frag_list_file.checkTimestamps();
					if( !frag_list_file.contains( frag_id ) ) {
						// copy the database fragments
						// .nhr, .nin, .nsq, .nnd, .nni, .nsd, .nsi for nucleotide
						// .phr, .pin, .psq, .pnd, .pni, .psd, .psi for amino acids
						bool copy_success = true;
						const vector<string>& extensions = FragmentExtensions(db_type);

						char fragid[8];
						memset( fragid, 0, 8 );
						// always support 3 digit fragment identifiers
						sprintf( fragid, "%03d", frag_id );
						if( debug_msg ){
							LOG_MSG << "copying fragment " << fragid << endl;
						}

						meta_MPE_Log_event(cpstart,0,"start cMasteropy");
						if (copy_via != COPY_VIA_MPI) { // copy via cp/rcp/scp

	/*
							// request copy fragment
							event_array[0] = REQUEST_FRAGMENT_COPY;
							event_array[1] = -1;
							MPI_Send(&event_array, 2, MPI_INT, super_master_process, EVENT_TYPE, MPI_COMM_WORLD);

							int copy_command = -1;
							MPI_Status recv_status;
							MPI_Recv( &copy_command, 1, MPI_INT, super_master_process, ASSIGNMENT_TYPE, MPI_COMM_WORLD, &recv_status);
							if(copy_command != CONFIRM_COPY_FRAGMENT) {
								throw __FILE__ "message out of order while waiting for copy command";
							}
	*/

							for( uint extI = 0; extI < extensions.size(); extI++ ){
								string copy_fragment_name = database_name + ".";
								copy_fragment_name += fragid + extensions[ extI ];
								string copy_src = config.sharedPath() + copy_fragment_name;
								string copy_dest = config.localPath() + copy_fragment_name;

								if(use_virtual_frags) {
									int ret = VFM::Instance()->LoadVirtualFile(copy_src, copy_dest);
									if(ret != 0) {
										copy_success = false;
									}
								} else {
									if( copyFile( copy_src, copy_dest, copy_via ) != 0 )
										copy_success = false;
								}
							}
						} else { // copy via mpi
							copy_success = recvFragsViaMPI( frag_id );
						}
						meta_MPE_Log_event(cpend,0,"end copy");

						if( copy_success ){
							frag_list_file.addFragment( frag_id );
						}
						else{
							LOG_MSG << "(" << my_rank << ") unable to copy fragment!" << endl;
							MPI_Abort( MPI_COMM_WORLD, -1 );
						}
					}	// end "if( !fragment_list_file.contains( frag_id ) )"

					copy_time += MPI_Wtime() - track_time;

					//
					// Tell the scheduler that the fragment copy completed successfully
					//
					event_array[0] = FRAGMENT_COPY_COMPLETE;
					event_array[1] = int_size;

					CommSendStruct* cssp_event = new CommSendStruct(2 * int_size);
					cssp_event->AddData(&event_array[0], 2 * int_size);
					cssp_event->IsendData(group_comm, scheduler_process, EVENT_TYPE);
					writer_w->AddCommSend(cssp_event);

					CommSendStruct* cssp_data = new CommSendStruct(int_size);
					cssp_data->AddData(&frag_id, int_size);
					cssp_data->IsendData(group_comm, scheduler_process, FRAGMENT_COPY_COMPLETE);
					writer_w->AddCommSend(cssp_data);

					// also inform supermaster
					CommSendStruct* cssp_event1 = new CommSendStruct(2 * int_size);
					cssp_event1->AddData(&event_array[0], 2 * int_size);
					cssp_event1->IsendData(MPI_COMM_WORLD, super_master_process, EVENT_TYPE);
					writer_w->AddCommSend(cssp_event1);

					free(assign_array);
					continue;
				} // end "if (assignment == COPY_FRAGMENT)"

				// Assignment is SEARCH_FRAGMENT
				// fill alias file with the names of fragments to search.
				//
				if( debug_msg )
					LOG_MSG << "Searching query " << assign_array[3] << " vs frag " << assign_array[2] << endl;

				// working receive query data from the scheduler and put it into query map
				int query_len = assign_array[4];
				if(query_len > 0 && !query_in_mem) {
					MPI_Status data_status;
					char* query_buf = new char[query_len + 1]; // will be freed in QueryM
					MPI_Recv(query_buf, query_len, MPI_BYTE, status.MPI_SOURCE, QUERY_DATA, group_comm, &data_status);
					query_buf[query_len] = NULLB;
					QueryM::Instance()->AddQueryData(query_id, query_buf);
				}

				// writer_w->AddWorkingQueries(segment_id, query_id, query_id + 1);
				writer_w->InitQueryOutput(segment_id, query_id);

				if(!preload_query) {
					QueryM::Instance()->UnloadOldQueries(query_id);
				}

				string frag_to_search;
				if(use_virtual_frags) {
					ostringstream tmp_os;
					char sufix[8];
					sprintf(sufix, "%03d", frag_id);
					tmp_os << config.localPath() << database_name << "." << sufix;
					frag_to_search = tmp_os.str();
				} else {
					vector<int> single_frag_vector(1, frag_id);
					WriteAliasFile(alias_filename, single_frag_vector);
				}

				//
				// initialize BLAST with the current data set
				//
				if(debug_msg) {
					LOG_MSG << "initBLAST() start" << endl;
				}
				running_blast = true;
				ncbiblast->Init(WORKER_MODE);
				if(debug_msg) {
					LOG_MSG << "initBLAST() end" << endl;
				}
				running_blast = false;

				if(use_virtual_frags) {
					updateBlastDB(frag_to_search.c_str(), frag_to_search.size());
				}

				//
				// execute blast
				//
				if( debug_msg ){
					LOG_MSG << "node " << my_rank << " Executing BLAST search\n";
				}
				meta_MPE_Log_event(blaststart,0,"start blastall");
				running_blast = true;
				//int retcode = ncbiblast->Run( SEARCH_MODE, query_id, query_id );
				writer_w->SetStartQuery(query_id);

				track_time = MPI_Wtime();

				int retcode = ncbiblast->RunPIO( SEARCH_MODE, query_id, query_id, frag_id );

				run_blast_time += MPI_Wtime() - track_time;

				running_blast = false;
				meta_MPE_Log_event(blastend,0,"end blastall");
				if( retcode != 0 ){
					LOG_MSG << "blastall exited with code: " << retcode << endl;
					MPI_Abort( MPI_COMM_WORLD, -1 );
				}

				meta_MPE_Log_event(cleanup_blast_start,0,"cleanupBLAST() start");
				running_blast = true;
				ncbiblast->Cleanup();
				running_blast = false;
				meta_MPE_Log_event(cleanup_blast_end,0,"cleanupBLAST() end");

				if( debug_msg ){
					LOG_MSG << "node " << my_rank << " blast retcode= " << retcode << endl;
				}

				writer_w->ProcessResults(true, frag_id);

				if(writer_w->GetNumPendingOffsets() >= max_pending_offsets) {
					send_idle_message = false;
					if(debug_msg) {
						LOG_MSG << "Reached the max number of pending offsets " << max_pending_offsets << ", stop fetching new assignments." << endl;
					}
				}

			}else{
				// assignment was DO_NOTHING...
				// in the future the worker shouldn't ever see this message.
				// for now, just have it sleep briefly to avoid swamping the
				// scheduler with requests
#ifndef WIN32
				// sleep(1);
				usleep(200000);
#else
				Sleep(1000);
#endif
			}
		}
		free(assign_array);

	} // end while (true) worker loop

	// search is now complete
	// finish any pending send operations
	// finish output all queries
	writer_w->Finalize();

//	while(true){
//		int array_size = 0;
//		MPI_Probe( scheduler_process, ASSIGNMENT_TYPE, group_comm, &status );
//		MPI_Get_count( &status, MPI_INT, &array_size );
//		assign_array = (int*)malloc(int_size * array_size);
//		CHECK_NULPTR(assign_array);
//		MPI_Recv( assign_array, array_size, MPI_INT, scheduler_process, ASSIGNMENT_TYPE, group_comm, &status );
//
//		if(assign_array[0] == WORKER_QUIT) {
//			break;
//		}
//
//		free(assign_array);
//	}

    int flag = 0;
    int pool_count = 0;
	int array_size = 0;
    while (true) {
        MPI_Iprobe( MPI_ANY_SOURCE, ASSIGNMENT_TYPE, group_comm, &flag, &status );

        if(flag) {
            MPI_Get_count( &status, MPI_INT, &array_size );
            assign_array = (int*)malloc(int_size * array_size);
            CHECK_NULPTR(assign_array);
            MPI_Recv( assign_array, array_size, MPI_INT, scheduler_process, ASSIGNMENT_TYPE, group_comm, &status );

            if(assign_array[0] == WORKER_QUIT) {
                break;
            }

            if(debug_msg) {
                LOG_MSG << "Warning: received assignment " << assign_array[0] << endl;
            }
            free(assign_array);
        }

//        if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
//            writer_w->ProcessAsyncOutput();
//        }

/*        pool_count++;
        if(pool_count > 50) {
            usleep(50000);
            pool_count = 0;
		}*/
    }

    free(assign_array);

    // free query sequence memory
	if(use_query_map) {
		QueryM::Instance()->Destroy();
	} else {
		cleanupQueries();
	}
}

void  MpiBlast::WriteAliasFile( const string& alias_filename, const vector< int >& fragment_list)
{
	if(fragment_list.size() == 0){
		return;	// don't write a file
		// this isn't an error condition when # workers > # db frags
//		throw __FILE__ "(WriteAliasFile) Empty fragment_list";
	}

	ofstream alias_file( alias_filename.c_str() );

	if( !alias_file.is_open() ){
		LOG_MSG << "Error opening " << alias_filename << endl;
		throw __FILE__ "(MpiBlast::worker): Unable to open alias file";
	}

	alias_file << "TITLE " << config.localPath() << database_name << endl;
	alias_file << "DBLIST";

	for( uint iter = 0; iter != fragment_list.size(); iter++ ){
		alias_file << " " << database_name << ".";
		char fragid[8];

		memset(fragid, 0, 8);

		// always use 3 digit fragment identifiers
		sprintf( fragid, "%03d", fragment_list[ iter ] );

		alias_file << fragid;
		if( debug_msg ){
			LOG_MSG << "Will search "<<database_name << "." << fragid << endl;
		}
	}
	alias_file << endl;
	alias_file.close();

}

void MpiBlast::DistributeDB() {

	// form distribute communication
	int color = -1;
	if(group_rank == 0) {
		color = node_count; // make sure it has different color than others
	} else {
		color = (group_rank - 1) % replica_group_size;
	}
	MPI_Comm distribute_comm;
	int distribute_rank;

    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &distribute_comm);
    MPI_Comm_rank(distribute_comm, &distribute_rank);

	if(group_rank == 0) {   // skip the first rank (supermaster and master)
		return;
	}

	if(debug_msg) {
		LOG_MSG << "rank=" << my_rank << ", group_rank=" << group_rank << ", distribute_color=" << color << ", distribute_rank=" << distribute_rank << endl;
	}

	vector <string> local_frags;

	// assign frags among replica group
	for(int i=0; i<total_frags; i++) {
		char fragid[8];
		memset( fragid, 0, 8);
		sprintf( fragid, "%03d", i);
		if(i % replica_group_size == color) {
			local_frags.push_back(fragid);
		}
	}

	if(debug_msg) {
		ostringstream tmp_os;
		copy(local_frags.begin(), local_frags.end(), ostream_iterator<string>(tmp_os, " "));
		LOG_MSG << "Assigned local frags: " << tmp_os.str() << endl;
	}

	const vector<string>& extensions = FragmentExtensions(db_type);

	for(int i=0; i<local_frags.size(); i++) {
		for(int j=0; j<extensions.size(); j++) {
			// note file_len is integer here, cannot support fragment larger than 2GB for now
			int file_len = 0;
			string file_name = database_name + "." + local_frags[i] + extensions[j];
			string shared_file_path = config.sharedPath() + file_name;
			ifstream ifs;

			if(distribute_rank == 0) {
				file_len = statFileSize(shared_file_path.c_str());

				if(debug_msg) {
					LOG_MSG << "Distributing " << file_name << endl;
				}

				ifs.open(shared_file_path.c_str());
				if(!ifs.is_open()) {
					throw __FILE__ "MpiBlast::DistributeDB() - cannot open input file";
				}
			}
			MPI_Bcast(&file_len, 1, MPI_INT, 0, distribute_comm);

			if(file_len == 0) {
				continue;
			}

			// store file in virtual file
			// create virtual file
			string local_file_path = config.localPath() + file_name;
			NlmMFILEPtr mfp = VFM::Instance()->InsertMapFile(local_file_path, MAPVIRTUAL, file_len);

			char* ptr = (char*)(mfp->mmp_begin);
			int total_read = 0; // not supporting file > 2GB

			while(total_read < file_len) {
				int curr_read = 0;

				if(distribute_rank == 0) {
					if(total_read + max_data_to_send <= file_len) {
						curr_read = max_data_to_send;
					} else {
						curr_read = file_len - total_read;
					}
					ifs.read(ptr, curr_read);
				}

				MPI_Bcast(&curr_read, 1, MPI_INT, 0, distribute_comm);
				MPI_Bcast(ptr, curr_read, MPI_BYTE, 0, distribute_comm);

				ptr += curr_read;
				total_read += curr_read;
			}

			// if use local storage, dump vitual file into local and delete it
			if(!use_virtual_frags) {
				ofstream ofs(local_file_path.c_str());
				if(!ofs.is_open()) {
					throw __FILE__ "MpiBlast::DistributeDB() - cannot open output file";
				}

				ofs.write((char*)(mfp->mmp_begin), file_len);
			}
		}
	}

	// creat mbf file
	for(int i=0; i<local_frags.size(); i++) {
		int frag_id = atoi(local_frags[i].c_str());
		frag_list_file.addFragment(frag_id);
	}
}

void MpiBlast::ReadDBSpecs(const string& spec_file) {
	ifstream ifs(spec_file.c_str());

	if(!ifs.is_open()) {
		throw __FILE__ " MpiBlast::ReadDBSpecs -- error opening of database specification file. Please reformat your database with the updated version of mpiformatdb.";
	}

	ifs >> global_db_len;
	ifs >> global_dbseq_num;
	ifs.close();

	use_real_dblen = true;

	return;
}

void MpiBlast::CreateGroupComm() {
	// create group communicators
	int color = -1;
	color = ( my_rank == super_master_process? node_count : GroupManager::Instance()->GetGroupId(my_rank)); // make sure supermaster has different color
	MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &group_comm);
	MPI_Comm_rank(group_comm, &group_rank);
	MPI_Comm_size(group_comm, &group_node_count);

	// initialize write group

	if(scheduler_process == writer_process) {
		group_write_comm = group_comm;
	}

	if(replica_group_size == 0) {
		replica_group_size = total_frags;
	}

	if(replica_group_size > partition_size - 1) {
		replica_group_size = partition_size - 1;
	}

	if(replica_group_size > total_frags) {
		replica_group_size = total_frags;
	}

	// set default replica count
	int num_copies = GetNumWorkers() / replica_group_size;
	int remain = GetNumWorkers() % replica_group_size;

	if(remain > 0) {
		num_copies += 1;
	}

	if(num_copies == 0) {
		num_copies = 1;
	}

	if(my_rank != super_master_process) {
		PrecopySchedulerPolicy::db_copy_number = num_copies;
		if(frags_per_worker > 0) {
			PrecopySchedulerPolicy::node_copy_number = frags_per_worker;
		} else {
			if(!predistribute_db || increment_db_copy) {
				int max_frags = (int)ceil((double)num_copies * (double)total_frags / (double)(group_node_count - 1));
				if(use_virtual_frags) {
					PrecopySchedulerPolicy::node_copy_number = max_frags > 0? max_frags : 0;
				}
			} else {
				PrecopySchedulerPolicy::node_copy_number = total_frags / replica_group_size;
			}
		}

		if(debug_msg) {
			LOG_MSG << "Number of cached fragments per node: " << PrecopySchedulerPolicy::node_copy_number << endl;
		}
	}

	if(debug_msg) {
		LOG_MSG << "Group Info: rank=" << my_rank << ", group_id=" << color << ", group_rank=" << group_rank << ", group_node_count=" << group_node_count << endl;
	}

	return;
}

void MpiBlast::ParseArguments (int argc, char* argv[]) {
	//
	// must read -d, -i, and -o options
	//
	const char* short_opts = "-d:i:o:p:m:O:l:";
	int ac = argc;
	char** av = argv;
	int opt_code, opt_type = 0;
	opterr = 0;
//	config_f_name = MpiBlastConfig::defaultConfigFileName();
	config_f_name = ".ncbirc";
	int config_opt = 0;
	int long_index = 0;
	int niceness = 0;
	int db_replica_count = 0;
	bool lock_set = false;
	bool io_set = false;
	bool query_in_file_set = false;
	bool replica_size_set = false;
	partition_size = node_count - 1; // default one group
	hpm_filename = "mpiblast";

	writer_opts.push_back( "blastall" );
	scheduler_opts.push_back( "blastall" );
	worker_opts.push_back( "blastall" );
	struct option long_opts[] = {
		{"config-file", required_argument, &config_opt, 1 },
		{"debug", optional_argument, &config_opt, 2 },
		{"pro-phile", required_argument, &config_opt, 3 },
		{"removedb", no_argument, &config_opt, 4 },
		{"version", no_argument, &config_opt, 5 },
		{"disable-mpi-db", no_argument, &config_opt, 6 },
		{"nice", required_argument, &config_opt, 7 },
		{"concurrent", required_argument, &config_opt, 8},
		{"mpi-size", required_argument, &config_opt, 9},
		{"copy-via", required_argument, &config_opt, 10},
		{"lock", optional_argument, &config_opt, 11},
		{"resume-run", no_argument, &config_opt, 12},
		{"obselete-arg1", required_argument, &config_opt, 13},
		{"altschul-reference", no_argument, &config_opt, 14},
		{"db-replicate-count", required_argument, &config_opt, 15},
		{"disable-posix-lock", optional_argument, &config_opt, 16 },
		{"dump-raw-output", no_argument, &config_opt, 17 },
		{"max-write-buffer", required_argument, &config_opt, 18 },
		{"fast-evalue-approximation", no_argument, &config_opt, 19},
		{"time-profile", required_argument, &config_opt, 20 },
		{"output-search-stats", no_argument, &config_opt, 21},
		{"output-strategy", required_argument, &config_opt, 22},
		{"io", required_argument, &config_opt, 23},
		{"concurrent-write", required_argument, &config_opt, 24},
		{"partition-size", required_argument, &config_opt, 25},
		{"sync-comm", no_argument, &config_opt, 26},
		{"use-virtual-frags", no_argument, &config_opt, 27},
		{"use-evalue-dryrun", no_argument, &config_opt, 28},
		{"query-segment-size", required_argument, &config_opt, 29},
		{"replica-group-size", required_argument, &config_opt, 30},
		{"predistribute-db", no_argument, &config_opt, 31},
		{"query-in-file", no_argument, &config_opt, 32},
		{"num-pending-writes", required_argument, &config_opt, 33},
		{"use-parallel-write", no_argument, &config_opt, 34},
		{"real-db-len", required_argument, &config_opt, 35},
		{"real-dbseq-num", required_argument, &config_opt, 36},
		{"init-assign-percent", required_argument, &config_opt, 37},
		{"frags-per-worker", required_argument, &config_opt, 38},
		{"use-query-map", no_argument, &config_opt, 39},
		{"lazy-load-query", no_argument, &config_opt, 40},
		{"obselete-arg2", required_argument, &config_opt, 41},
		{"increment-db-copy", no_argument, &config_opt, 42},
		{"hpm-filename", required_argument, &config_opt, 43},
		{"query-in-mem", no_argument, &config_opt, 44},
		{"query-load-size", required_argument, &config_opt, 45},
		{0,0,0,0}	// signifies termination of option list
	};

	while((opt_code = getopt_long( ac, av, short_opts, long_opts, &long_index ) ) != EOF ){
		switch( opt_code ){
			case 0:
				switch( config_opt ){
					case 1:	// --config-file
						config_f_name = optarg;
						break;
					case 2:	// --debug
						if( optarg != NULL )
							log_file_name = optarg;
						debug = true;
						debug_msg = true;
						debug_bsfetch = true;
						debug_bslookup = true;
						break;
					case 3:	// was pro-phile -- now deprecated
						cerr << "WARNING: --pro-phile is no longer supported\n";
						break;
					case 4:	// --removedb
						remove_db = true;
						break;
					case 5:	// --version
						cout << PACKAGE << " version " << VERSION << endl;
						break;
					case 6:	// --disable-mpi-db
						// This option is not used in mpiblast/PIO.
						//disable_mpi_db = true;
						break;
					case 7:	// --nice
						niceness = atoi(optarg);
						break;
					case 8:	// --concurrent
						concurrent_accesses = atoi(optarg);
						if( concurrent_accesses > 1 && !lock_set )
							lock_mbf = true;
						break;
					case 9:	// --mpi-size
						max_data_to_send = (uint64) atoi(optarg);
						break;
					case 10:	// --copy-via
						//options = cp, rcp, scp, none/pvfs/gpfs, mpi
						copy_via_set = true;
						if (!strcmp(optarg,"cp")){
							copy_via = COPY_VIA_CP;
						} else if (!strcmp(optarg,"rcp")){
							copy_via = COPY_VIA_RCP;
						} else if (!strcmp(optarg, "scp")){
							copy_via = COPY_VIA_SCP;
						} else if (!strcmp(optarg, "none")){
							copy_via = COPY_VIA_NONE;
						} else {
							copy_via = COPY_VIA_MPI;
						}
						break;
					case 11:	// --lock
						lock_mbf = true;
						if( optarg != NULL && strcmp( optarg, "off" ) )
							lock_mbf = false;
						lock_set = true;
						break;
					case 12:	// --resume-run
						cerr << "\"--resume-run\" currently not supported!" <<endl;
						exit(1);
						//resume_run = true;
						break;
					case 13:	// --scheduler-rank
						// cerr << "\"--scheduler-rank\" currently not supported!" <<endl;
						// exit(1);
						//scheduler_process=atoi(optarg);
						//writer_process=scheduler_process;
						break;
					case 14:	// --altschul-reference
						strict_output_conformance = TRUE;
						break;
					case 15:	// --db-replicate-count
						// PrecopySchedulerPolicy::db_copy_number = atoi(optarg);
						db_replica_count = atoi(optarg);
						replica_count_set = true;
						break;
					case 16:
						disable_posix_lock = true;
						break;
					case 17:
						dump_raw_output = true;
						break;
					case 18:
						max_write_buffer = atoi(optarg);
						break;
					case 19:
						fast_evalue_approximation = true;
						break;
					case 20:
						collect_profile = true;
						profile_filename = optarg;
						break;
					case 21:
						use_brief_report = 0;
						break;
					case 22:
						if(optarg == string("master_streamline")) {
							output_strategy = MASTER_STREAMLINE;
						} else if (optarg == string("worker_split")) {
							output_strategy = WORKER_SPLIT;
						} else if (optarg == string("worker_individual")) {
							output_strategy = WORKER_INDIVIDUAL;
						} else if (optarg == string("worker_collective")) {
							output_strategy = WORKER_COLLECTIVE;
						} else if (optarg == string("worker_streamline")) {
							output_strategy = WORKER_STREAMLINE;
						} else {
							cerr << "Invalid output strategy!!!";
							exit(-1);
						}
						break;
					case 23:
						if(optarg == string("mpi")) {
							io_function = MPI_IO_FUNC;
						} else if (optarg == string("posix")) {
							io_function = POSIX_IO_FUNC;
						}

						io_set = true;
						break;
					case 24:
						concurrent_write = atoi(optarg);
						break;
					case 25:
						partition_size = atoi(optarg);
						if(partition_size < 1) {
							cerr << "Invalid processor group size" << endl;
							exit(-1);
						}
						if(partition_size > node_count - 1) {
							partition_size = node_count - 1;
						}
						break;
					case 26:
						sync_comm = true;
						break;
					case 27:
						use_virtual_frags = 1;
						break;
					case 28:
						parallel_evalue_adjust = 0;
						break;
					case 29:
						query_segment_size = atoi(optarg);
						query_segment_size_set = true;
						break;
					case 30:
						replica_group_size = atoi(optarg);
						replica_size_set = true;
						break;
					case 31:
						predistribute_db = true;
						break;
					case 32:
						query_in_file_set = true;
						break;
					case 33:
						num_pending_writes = atoi(optarg);
						break;
					case 34:
						output_strategy = WORKER_STREAMLINE;
						break;
					case 35:
						global_db_len = atol(optarg);
						use_real_dblen = true;
						break;
					case 36:
						global_dbseq_num = atoi(optarg);
						break;
					case 37:
						init_assign_percent = atof(optarg);
						if(init_assign_percent > 1.0) {
							init_assign_percent = 1.0;
						}
						break;
					case 38:
						frags_per_worker = atoi(optarg);
						break;
					case 39:
						use_query_map = 1;
						break;
					case 40:
						preload_query = 0;
						break;
					case 41:
						break;
					case 42:
						increment_db_copy = true;
						break;
					case 43:
						hpm_filename = optarg;
						break;
					case 44:
						query_in_mem = 1;
						break;
					case 45:
						max_query_load = atoi(optarg);
						break;
				}
				break;

			case 'i':
				query_filename = optarg;
				addOpt( writer_opts, opt_code, optarg );
				break;
			case 'o':
				output_filename = optarg;
				global_output_file = optarg;
				//addOpt( writer_opts, opt_code, optarg );
				break;
			case 'd':
				database_name = optarg;
				break;
			case 'p':
				addOpt( writer_opts, opt_code, optarg );
				addOpt( worker_opts, opt_code, optarg );
				blast_type = optarg;
				break;
			case 'm':
				// worker always uses -m 11
				output_format = atoi( optarg );
				addOpt( writer_opts, opt_code, optarg );
				addOpt( worker_opts, opt_code, optarg );
				break;
			case 'O':
				// Only the master is allowed to see this
				// flag since workers must use -m 11. Don't let the
				// -O force ASCII ASN.1 output on the workers!
				addOpt( writer_opts, opt_code, optarg );
				break;
			case 'l':
				// this specifies a gi list to use when searching
				// for now just assume it's located on shared storage
				filter_filename = optarg;
				addOpt( writer_opts, opt_code, optarg );
				break;
			case 1:
				worker_opts.push_back( optarg );
				writer_opts.push_back( optarg );
				break;
			case '?':
				addOpt( worker_opts, optopt, optarg );
				addOpt( writer_opts, optopt, optarg );
				opt_type = optopt;
				break;
			default:
				addOpt( worker_opts, opt_type, optarg );
				addOpt( writer_opts, opt_type, optarg );
				opt_type = optopt;
				break;
		}
	}

//	if(use_virtual_frags && !query_in_file_set) {
//		use_query_map = true;
//		preload_query = false;
//	}
//
	if(use_virtual_frags) {
		predistribute_db = true;
	}

	if(copy_via == COPY_VIA_MPI && use_virtual_frags) {
		cerr << "copy-via cannot be mpi if use-virtual-frags is specified" << endl;
		exit(-1);
	}

	if(!io_set) {
		if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_COLLECTIVE) {
			io_function = MPI_IO_FUNC;
		} else {
			io_function = POSIX_IO_FUNC;
		}
	}

	if(replica_count_set && replica_size_set) {
		cerr << __FILE__ " db-replicate-count and replica-group-size cannot be both specified." << endl;
		exit(-1);
	}

	if(replica_count_set) { // convert to db_replica_count to equivalent replica_group_size
		replica_group_size = (partition_size - 1) / db_replica_count;
	}

	if(output_strategy == MASTER_STREAMLINE) {
		disable_posix_lock = true;
	}

#ifndef WIN32
	//Set the niceness if given as cmd line arg & we're in unix
	if (niceness > 0)
		setpriority(PRIO_PROCESS, 0, niceness);
#endif

}

int MpiBlast::CalculateFirstAssignment() {
/*
	int num_first_assign = 0;
	int num_groups = GroupManager::Instance()->GetNumGroups();

	if(num_groups == 1) {
		num_first_assign = query_count;
	} else {
		if(init_assign_percent > 0.0) {
			num_first_assign = (int)ceil((double)query_count * (double)init_assign_percent / (double)num_groups);
		} else {
			num_first_assign = (int)ceil((double)(partition_size - 1) / (double)total_frags);
		}
	}

	if(num_first_assign < query_segment_size) {
		num_first_assign = query_segment_size;
	}

	// to do: adjust assignment of last group
	return num_first_assign;
*/
	return query_segment_size;
}

// need getpid for debugging under win32
//#include <process.h>

int MpiBlast::main( int argc, char* argv[] )
{
#if defined(WIN32) && defined(_DEBUG)
	// for debugging under win32, wait so the process can be attached.
//	LOG_MSG << "Process id: " << _getpid() << endl;
	cerr << "Gentlemen, start your debuggers!\n";
	Sleep( 7000 );
#endif

	if( argc == 0 ){
		cerr << "Incorrect usage\n";
		MPI_Abort( MPI_COMM_WORLD, -1 );
		throw __FILE__ "(MpiBlast::main): Incorrect usage";
	}
	exec_path = argv[0];
	exec_path = getPath( exec_path );

	ParseArguments(argc, argv);

	try{
		if(my_rank == super_master_process) {
			int buf_len = 0;
			config = MpiBlastConfig( config_f_name );
			string local = config.localPath();

			buf_len = local.size();
			MPI_Bcast( &buf_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);
			MPI_Bcast( (char*)(local.c_str()), buf_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);

			string shared = config.sharedPath();
			buf_len = shared.size();
			MPI_Bcast( &buf_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);
			MPI_Bcast( (char*)(shared.c_str()), buf_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);
		} else {
			int buf_len = 0;
			char* buf;

			MPI_Bcast( &buf_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);
			buf = new char[buf_len+1];
			MPI_Bcast( buf, buf_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);
			buf[buf_len] = 0;
			string local = buf;
			delete buf;

			MPI_Bcast( &buf_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);
			buf = new char[buf_len+1];
			MPI_Bcast( buf, buf_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);
			buf[buf_len] = 0;
			string shared = buf;
			delete buf;

			config = MpiBlastConfig( local, shared );
		}

		// if the user wants a "clean" run, create a temp
		// directory on the local storage device for the job
		if( remove_db )
			config.createLocalTempDir();

		localPath = config.localPath();
		if (config.localPath() == config.sharedPath() ) {
			// if the user didn't specify copy_via and the config
			// file has identical shared and local then default
			// to --copy_via=none for ease of use
			if( !copy_via_set )
				copy_via = COPY_VIA_NONE;
			if( copy_via != COPY_VIA_NONE ){
				// can't have same shared/local and not have copy_via=none
				cerr << "Error: You must set --copy_via=none when Shared and Local storage are identical\n";
				terminateProgram(1);
			}
			sameSharedLocal = 1;
		} else {
			if ( copy_via == COPY_VIA_NONE){
				// can't have different shared/local when copying via none
				cerr << "Error: Shared and Local storage must be identical when --copy_via=none\n";
				terminateProgram(1);
			}
			sameSharedLocal = 0;
		}

		db_type = "n"; //Default to nucleotide
		if( blast_type == "blastp" || blast_type == "blastx" ){
			db_type = "p";
		}

		ncbiblast->SetDBType(db_type);

		// set location of substitution matrices if not set by the user
		if( getenv( "BLASTMAT" ) == NULL ){
			string env_str = "BLASTMAT=" + config.sharedPath();
			setEnvironmentVariable( env_str.c_str() );
		}

		//
		// open the log file
		//
		if( log_file_name != "" ){
			ostringstream log_file_str;
			log_file_str << log_file_name << "." << my_rank;
			log_file.open( log_file_str.str().c_str() );

			if( !log_file.is_open() ){
				cerr << "Error opening log file \"" << log_file_name << "\"\n";
				throw __FILE__ "(main): Unable to open log file";
			}

			if( debug_msg ){
				cout << "logging to " << log_file_str.str() << endl;
			}

			log_stream = &log_file;
		}

		if(debug_msg) {
			LOG_MSG << "log file created" << endl;
			char proc_name[MPI_MAX_PROCESSOR_NAME];
			int len;
			MPI_Get_processor_name(proc_name, &len);
			LOG_MSG << "Running on: " << proc_name << endl;
		}

		//
		// add database name for master
		//
		// The old version
		// string shared_db_name = config.sharedPath() + database_name;
		// Assume that the database resides on either the master, or a valid nfs mount
		string shared_db_name;
		if( !isRemotePath( config.sharedPath() ) ){
			// NFS mount
			shared_db_name = config.sharedPath() + database_name;
		}
		else{
			// On master
			string::size_type pos = config.sharedPath().find(":", 0);
			pos ++;
			shared_db_name = config.sharedPath().substr(pos, config.sharedPath().length() - pos)
				+ database_name;
		}

        if(debug_msg) {
            LOG_MSG << "Shared = " << config.sharedPath() << endl;
            LOG_MSG << "Local = " << config.localPath() << endl;
            LOG_MSG << "shared_db_name = " << shared_db_name << endl;
        }

		ncbiblast->SetSharedDBName(shared_db_name);

		// have rank 0 count the number of fragments
		// and send the number to workers
		if( my_rank == super_master_process ){
			// MPI_File_delete((char*)output_filename.c_str(), MPI_INFO_NULL);

			string mbf_filename = shared_db_name + FRAG_LIST_EXTENSION;
			FragmentListFile tmp_frag_list_file( mbf_filename, config, database_name, db_type );
			total_frags = tmp_frag_list_file.fragmentCount();

			if(total_frags == 0) {
				throw __FILE__ "Error reading database mbf file";
			}

			MPI_Bcast( &total_frags, 1, MPI_UNSIGNED, super_master_process, MPI_COMM_WORLD );
		}else{
			MPI_Bcast( &total_frags, 1, MPI_UNSIGNED, super_master_process, MPI_COMM_WORLD );
		}

		addOpt( writer_opts, 'd', shared_db_name.c_str() );

		query_file = getFilename( query_filename );
		local_query_filename = config.localPath() + query_file;

		// for gi list filters:
		filter_file = getFilename( filter_filename );
		local_filter_filename = config.localPath() + filter_file;

/*
		//read the output file, and return the number of the last query w/output -- X
		//read off the first X queries, and place the remaining ones into a new query file
		if(resume_run){
			if(!doesFileExist(output_filename) || (output_format == 7) ||(output_format == 8) || (output_format == 10) || (output_format == 11) ){
				LOG_MSG << "To resume a run the outputfile must already exist\n";
				LOG_MSG << "Only -m 1-6 & 9 are supported for --resume\n";
				terminateProgram(-1);
			}
			unsigned int lastquerynum = determineLastQueryNum(output_filename);
			string new_qf = config.localPath() + query_file + "XXXXXX";
			makeNewQuery(new_qf,query_filename, lastquerynum);
			// Overwrite the query value in the writer's argument list
			addOpt( writer_opts, 'i', new_qf.c_str() );
			registerFileToDelete( new_qf );
			query_filename = new_qf ;
		}
*/

		GroupManager::Instance(partition_size, total_frags);
		if(debug_msg) {
			if(my_rank == super_master_process) {
				LOG_MSG << "Printing group map: " << endl;
				GroupManager::Instance()->PrintGroupMap(*log_stream);
				LOG_MSG << endl;
			}
		}

		CreateGroupComm();

		if(group_rank != 0 && my_rank != super_master_process) {
			string fragment_filename = config.localPath() + database_name + FRAG_LIST_EXTENSION;
			if(use_virtual_frags) {
				frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type, true );
			} else {
				frag_list_file = FragmentListFile( fragment_filename, config, database_name, db_type );
			}
		} else {
			use_virtual_frags = 0; // disable virtual frags otherwise
		}

		if(predistribute_db) {
			double distribute_start_time = MPI_Wtime();

			if(debug_msg) LOG_MSG << "Start distributing db" << endl;

			DistributeDB();

			if(debug_msg) LOG_MSG << "End distributing db" << endl;
			MPI_Barrier(MPI_COMM_WORLD);

			db_distribute_time = MPI_Wtime() - distribute_start_time;
		}

		// Have the supermaster process send the query and gi filter file to all the other processes
		if( my_rank == super_master_process ){

			if( !isRemotePath( query_filename ) || copy_via == COPY_VIA_MPI){
				// The query is either (a) on the masters local disk or (b) visible to the master
				// via an NFS mount. No path name adjustment is required.
				// JA -- or (c) when using MPI to copy, rank 0 must have direct access to query
			}
			else{
				// The query resides on another host (not local to the master). Use rcp to copy the query to the
				// masters local disk.
				if( copyFile( query_filename, local_query_filename, copy_via ) != 0 ){
					LOG_MSG << "Unable to copy " << query_filename << " to "
					     << local_query_filename << endl;

					throw __FILE__ "(main): Unable to copy query to master";
				}

				// Make sure the master now reads from the local copy
				query_filename = local_query_filename;

				// Overwrite the query value in the writer's argument list
				addOpt( writer_opts, 'i', local_query_filename.c_str() );
				registerFileToDelete( local_query_filename );
			}

			//
			// broadcast the query to every node
			//
			double broadcast_track_time = MPI_Wtime();

			if(debug_msg) {
				LOG_MSG << "Start broadcasting query file" << endl;
			}

			ncbiblast->SetQueryFile(query_filename);
			QueryM::Instance(query_filename, query_in_mem);

			if(query_in_mem) {
				query_file_len = statFileSize(query_filename.c_str()) + 1;

				MPI_Bcast(&query_file_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);

				NlmMFILEPtr mfp = VFM::Instance()->InsertMapFile(query_filename, MAPVIRTUAL, query_file_len);
				char* ptr = (char*)(mfp->mmp_begin);
				ifstream ifs(query_filename.c_str());
				if(!ifs.is_open()) {
					throw __FILE__ "cannot open query file";
				}
				ifs.read(ptr, query_file_len - 1);

				ptr[query_file_len - 1] = NULLB;

				MPI_Bcast(ptr, query_file_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);

			}

			if(debug_msg) {
				LOG_MSG << "End broadcasting query file" << endl;
			}

			query_broadcast_time = MPI_Wtime() - broadcast_track_time;

			// filter file handling
			if ( filter_filename != "" ){
				// MC - This is crudely cut and pasted from the query file broadcast.
				// It ought to be unified and refactored
				// Have the rank 0 process send the filter file to all the other processes
				if( !isRemotePath( filter_filename ) || copy_via == COPY_VIA_MPI ){
					// The filter is either (a) on the masters local disk or (b) visible to the master
					// via an NFS mount. No path name adjustment is required.
					// JA -- or (c) when using MPI to copy, rank 0 must have direct access to filter
				}
				else{
					// The filter resides on another host (not local to the master).
					// Use rcp to copy the filter to the
					// masters local disk.
					if( copyFile( filter_filename, local_filter_filename, copy_via ) != 0 ){
						LOG_MSG << "Unable to copy " << filter_filename << " to "
							<< local_filter_filename << endl;

						throw __FILE__ "(main): Unable to copy filter to master";
					}

					// Make sure the master now reads from the local copy
					filter_filename = local_filter_filename;

					// Overwrite the query value in the writer's argument list
					addOpt( writer_opts, 'l', local_filter_filename.c_str() );
					registerFileToDelete( local_filter_filename );
				}

				//
				// broadcast the filter to every node
				//
				meta_MPE_Log_event(bcast_filter_start,0,"broadcast filter start");
				broadcastFile( filter_filename );
				meta_MPE_Log_event(bcast_query_end,0,"broadcast filter end");
			}

			if(debug_msg) {
				LOG_MSG << "Start initNCBI()" << endl;
			}

			running_blast = true;
			addOpt( writer_opts, 'o', bitbucket );
			ncbiblast->InitNCBI( writer_opts );
			ncbiblast->Init(MASTER_MODE);
			running_blast = false;

			if(debug_msg) {
				LOG_MSG << "End initNCBI()" << endl;
			}

			QueryM::Instance(query_filename, query_in_mem);
			query_count = QueryM::Instance()->IndexQueries();

			QueryM::Instance()->BcastIndexes(super_master_process);

			ncbiblast->SetQueryCount(query_count);

		// broadcast the array size
			MPI_Bcast( &query_count, 1, MPI_INT, super_master_process, MPI_COMM_WORLD );

		// collect and distribute evalue adjustments
			// read global DB info from the DBS file
			ReadDBSpecs(shared_db_name + ".dbs");

			if(debug_msg) {
				LOG_MSG << "bcast query adjustments start" << endl;
			}
			ncbiblast->MasterDistributeEvalAdjusts();
			if(debug_msg) {
				LOG_MSG << "bcast query adjustments end" << endl;
			}
		}else{
			//
			// receive the query file broadcast
			//
			if(debug_msg) {
				LOG_MSG << "Start broadcasting query file" << endl;
			}

			if(query_in_mem) {
				query_filename = local_query_filename + "XXXXXX";
			}

//			if(!query_in_mem) {
//				// create local query tempfile
//				getTempFileName( query_filename );
//			}

			ncbiblast->SetQueryFile(query_filename);
			QueryM::Instance(query_filename, query_in_mem);

			if(query_in_mem) {
				query_file_len = 0;
				MPI_Bcast(&query_file_len, 1, MPI_INT, super_master_process, MPI_COMM_WORLD);

				NlmMFILEPtr mfp = VFM::Instance()->InsertMapFile(query_filename, MAPVIRTUAL, query_file_len);
				char* ptr = (char*)(mfp->mmp_begin);
				MPI_Bcast(ptr, query_file_len, MPI_BYTE, super_master_process, MPI_COMM_WORLD);
			}

			if( debug_msg ){
				LOG_MSG << "End broadcasting query file" << endl;
			}

			if ( filter_filename != "" ){
				//
				// receive the gi filter file broadcast
				//
				meta_MPE_Log_event(bcast_filter_start,0,"broadcast filter start");
				// create local filter tempfile
				filter_filename = local_filter_filename + "XXXXXX";
				getTempFileName( filter_filename );
				// receive the filter file broadcast
				recvBroadcastFile( filter_filename );
				// make sure the filter file gets deleted before exiting
				registerFileToDelete( filter_filename );
				if( debug_msg ){
						LOG_MSG << "Filter file received as " << local_filter_filename << endl;
				}
				meta_MPE_Log_event(bcast_filter_end,0,"broadcast filter end");
			}

			QueryM::Instance()->BcastIndexes(super_master_process);

			// receive the array size
			MPI_Bcast( &query_count, 1, MPI_INT, super_master_process, MPI_COMM_WORLD );

			ncbiblast->SetQueryCount(query_count);

			if(debug_msg) {
				LOG_MSG << "bcast query adjustments start" << endl;
			}
			ncbiblast->WorkerRecvEvalAdjusts();
			if(debug_msg) {
				LOG_MSG << "bcast query adjustments end" << endl;
			}

			if( debug_msg ){
				LOG_MSG << query_count << " adjustments and effective db sizes received successfully\n";
			}
		}

		if(!query_segment_size_set && partition_size == node_count - 1) {
			// if only one group, use one query segment
			query_segment_size = max_segment_size > query_count? query_count : max_segment_size;
		}

		//
		// have rank super_master_process send database fragment timestamps to all nodes
		//
		if(debug_msg) {
			LOG_MSG << "bcast fragment times start" << endl;
		}
		if( my_rank == super_master_process ){
			ncbiblast->MasterDistributeFragTimes();
		}else{
			ncbiblast->WorkerRecvFragTimes();
		}
		if(debug_msg) {
			LOG_MSG << "bcast fragment times end" << endl;
		}

		// open a per process output file to dump debug info
		if(dump_raw_output) {
			char tmpfile[100];
			sprintf(tmpfile, "%d.data", my_rank);
			dbgfp=fopen(tmpfile, "w");
		}

		// write one file for a group for now
		ostringstream tmp_os;
		int partition_num = GroupManager::Instance()->GetGroupId(my_rank);
		if(partition_size != node_count - 1) {
			tmp_os << output_filename << ".part" << GroupManager::Instance()->GetGroupId(my_rank);
			output_filename = tmp_os.str();
		}

		// if(rank != super_master_process && group_rank == scheduler_process) {
			// MPI_File_delete((char*)output_filename.c_str(), MPI_INFO_NULL);
		// }

		// prevent deadloop when some group has nothing to do
		num_query_segments = query_count / query_segment_size;

		int num_groups = GroupManager::Instance()->GetNumGroups();

		if(num_query_segments < num_groups) {
			query_segment_size = query_count / num_groups;
			num_query_segments = query_count / query_segment_size;
			if(query_count % query_segment_size > 0) {
				num_query_segments += 1;
			}
		}

		if(num_query_segments == 0) {
			num_query_segments = 1;
		}

		if(num_query_segments < num_groups) {
			throw __FILE__ " The number of query segments are less than the number of partitions. Please try to specify larger partition-size or smaller query-segment-size";
		}

		system_init_time = MPI_Wtime() - prog_start;

		if(my_rank == super_master_process) {
			super_master();
		} else {
			if( group_rank == scheduler_process ) {
				master();
			} else {
				worker();
			}
		}

		if(dump_raw_output) {
			fclose(dbgfp);
		}

	}
	catch(const char *error){
		cerr << "Fatal Error:" << endl;
		cerr << error << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );
		return -1;
	}catch(exception& e ){
		cerr << "Fatal Exception!" << endl;
		cerr << e.what() << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );
		return -1;
	}
/*	catch(...){
		cerr << "Fatal Unhandled Exception!" << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );

		return -1;
	}
*/
	// clean up
	cleanupNCBI();

	VFM::Instance()->Destroy();

	deleteRegisteredFiles();
	if (remove_db) {
		MPI_Barrier(MPI_COMM_WORLD);
		string local_no_path_sep = config.localPath().substr(0, config.localPath().size() - 1 );
		deleteMpiBlastFiles( config.localPath() );
		if( !removeDirectory( config.localPath() ) )
			perror( "rmdir" );
	}

	return 0;
}

unsigned int determineLastQueryNum(const std::string outputfile){
	string line;
	int lastquerynum = 0;
	streampos lastQpos;
	ifstream outputfile_in(outputfile.c_str());
	while( getline( outputfile_in, line)){
		if ((line.substr(0,9) == "# Query: ") ||(line.substr(0,7) == "Query= ") || (line.substr(0,9) == "<?xml ver") ) {
			lastquerynum++;
			lastQpos = outputfile_in.tellg();
			lastQpos -= outputfile_in.gcount();
		}
	}
	outputfile_in.close();
	ofstream outputfile_out(outputfile.c_str(),ios::app);
	outputfile_out.seekp(lastQpos);
	outputfile_out << "";
	outputfile_out.close();

	return(lastquerynum);
}

void makeNewQuery(std::string& new_query, std::string& old_query, unsigned int lastquerynum){
	string line;
	int success = 0;
	int queryN = 0;
	getTempFileName(new_query);
	ifstream qfile_in(old_query.c_str());
	ofstream qfile_out(new_query.c_str());
	while (getline(qfile_in, line)){
		if (line.substr(0,1) == ">")
			queryN++;
		if (queryN > lastquerynum){
			qfile_out << line << endl;
			success++;
		}
	}
	if (success == 0) {
		cerr << "Looks like we don't want to be resuming the run. I think we completed\n";
		terminateProgram(-1);
	}
}

