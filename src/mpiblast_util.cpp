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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"
#include "mpi_util.h"

#include "mpiblast_util.hpp"
#include "fragment_list.hpp"
#include "mpiblast_types.h"

extern "C"{
#include "ncbi.h"
#include "objseq.h"
#include "blastdef.h"
#include "blast_hooks.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"
}

#include <fstream>
#include <iomanip>

using namespace std;

/* Logging for debugging and profiling */
bool debug_msg = false;		/* Default, set to true with --debug */
std::ostream* log_stream; 		/**< Output stream for logging information */
int rank;			/**< Rank of the current MPI process */
int node_count;		/**< Number of MPI processes allocated for the job */
int group_rank;  // rank within a group
int group_node_count; // how manay procs in my group
double realstarttime ;		/* Given value immediately on first line of main, before MPI_Init */
double prog_start;		/* Given value after MPI_Init and MPE defs */
double prog_end ;		/* Given value right before MPI_Finalize */

// process roles:
int scheduler_process = 0;
int writer_process = 0;
int super_master_process = 0;

int max_write_buffer = 16 * 1024 * 1024; // maximum data volume allowed for a streamlined write operation
bool fast_evalue_approximation = false;
int parallel_evalue_adjust = 1;  // 0-false, 1-true
int output_strategy = 0; // MASTER_STREAMLINE
int io_function = 0;
int ncpus = 1; // specify how many processors per compute node
int concurrent_write = 0;
int query_segment_size = 5;
int max_segment_size = 1000;
double init_assign_percent = -1.0;
bool dump_raw_output = false;
int use_brief_report = 1;
bool disable_posix_lock = false;
bool sync_comm = false;
bool use_real_dblen = false;
int use_query_map = 1;
int query_in_mem = 0;
int preload_query = 0;
int num_query_segments = 0;
int num_pending_writes = 25;
int num_threads = 0;
int num_reserved_events = 10;
int max_query_load = 0;
MPI_Comm group_write_comm = MPI_COMM_WORLD;
MPI_Comm group_comm;

const int int_size = sizeof(int);
const int float_size = sizeof(Nlm_FloatHi);
const int offset_size = sizeof(MPI_Offset);

// debug
FILE* dbgfp;
double search_time = 0;
double copy_time = 0;
double result_time = 0;
double run_blast_time = 0;
double write_time = 0;
double output_info_time = 0;
double query_broadcast_time = 0;
double idle_time = 0;
double master_idle_time = 0;
double scheduler_idle_time = 0;
double writer_master_idle_time = 0;
double scheduler_handle_time = 0;
double writer_master_handle_time = 0;
double writer_worker_handle_time = 0;
double results_gather_time = 0;
double wait_result_time = 0;
double wait_write_time = 0;
double receive_msg_time = 0;
double receive_msg_size = 0;
double db_distribute_time = 0;
double acquire_segement_time = 0;
double load_queries_time = 0;
double default_output_time = 0;
double process_output_time = 0;
double curr_process_output_time = 0;

double scheduler_called = 0;
int max_num_write = 0;
int peak_pending_offsets = 0; // the peak number of pending offsets, only used by workers

int max_pending_offsets = 50; // the maximum number of pending offsets, only used by workers

GroupManager* GroupManager::_instance = NULL;

int IsWorker(int src) {
	if(src == scheduler_process || src == writer_process) {
		return 0;
	} 
	
	return 1;
}

int GetNumWorkers() {
	if(scheduler_process == writer_process) {
		return group_node_count - 1;
	} else {
		return group_node_count - 2;
	}
}

void initNCBI( vector< string >& ncbi_opts ){
	//
	// initialize ncbi library
	//
	size_t ncbi_argc = ncbi_opts.size();
	//char** ncbi_argv = new char*[ ncbi_opts.size() ];
	char** ncbi_argv = NULL;
	
	if( debug_msg ) {
		LOG_MSG << "initializing ncbi ...";
	}
	
	if(ncbi_argc > 0){
		ncbi_argv = (char**)calloc(ncbi_argc, sizeof(char*));
	}
	
	for(uint optI = 0; optI < (uint)ncbi_argc; optI++ ){
		if( debug_msg )
			CONT_LOG_MSG << ncbi_opts[ optI ] << " ";
		//ncbi_argv[ optI ] = new char[ ncbi_opts[ optI ].size() + 1 ];
		ncbi_argv[ optI ] =(char*)calloc((ncbi_opts[ optI ].size() + 1), sizeof(char));
		strcpy( ncbi_argv[ optI ], ncbi_opts[ optI ].c_str() );
	}

	if( debug_msg )
		CONT_LOG_MSG << endl;
		
	Nlm_SetupArguments( ncbi_argc, ncbi_argv );
	
	// Note! Do not clean up the memory allocated in this function!
	// (i.e. ncbi_argv and ncbi_argv[i]). It is now owned by the
	// ncbi library. It is not clear if it needs to be deleted with
	// a call to Nlm_FreeCmdLineArguments(char** argv) in cleanupNCBI().
	// Since this would require saving the pointer, ncbi_argv, for
	// future reference, just forget about it (not a lot of memory
	// is lost anyway ...).
	
#ifdef MSC_VIRT
	if ( !_vheapinit(0, 1250, _VM_ALLSWAP) )
	{
		ErrPost( CTX_NCBIOBJ, 1, "Can't open virtual memory" );
		return 1;
	}
#endif

	/* Initialize connection library's logger, registry and lock */
	CONNECT_Init( 0 );
	
	if( debug_msg ) {
		LOG_MSG << "\n(" << rank << ") done initializing ncbi." << endl;
	}
}


void cleanupNCBI(){
	//
	// cleanup ncbi library
	//
	NlmThreadJoinAll();

	Nlm_FreeConfigStruct();
	ErrSetLogfile( NULL, 0 );
	Nlm_ReleaseAppContext();

#ifdef MSC_VIRT
	_vheapterm();
#endif
	NlmThreadDestroyAll();
}


void addOpt(vector< string >& opt_vector, int opt_character, const char* opt_argument)
{
	char opt_char = opt_character;
	string opt_string = "-";
	opt_string += opt_char;
	opt_vector.push_back( opt_string );
	if( opt_argument != NULL ){
		opt_vector.push_back( opt_argument );
	}
}

void SendIntVec(MPI_Comm comm, vector<int> &vec, int dest, int tag)
{
	size_t count = vec.size();
	
	MPI_Send(&count, 1, MPI_INT, dest, tag, comm);

	if(count > 0){
		int *tmp = (int*)calloc(count, sizeof(int));
		
		if(tmp == NULL){
			throw __FILE__ "SendIntVec: Unable to allocate memory for buffer";
		}
		
		for(size_t i = 0;i < count;i++){
			tmp[i] = vec[i];
		}

		MPI_Send(tmp, count, MPI_INT, dest, tag, comm);
			
		free(tmp);
	}

}

void RecvIntVec(MPI_Comm comm, vector<int> &vec, int src, int tag)
{
	int count;
	
	MPI_Recv(&count, 1, MPI_INT, src, tag, comm, MPI_STATUS_IGNORE);
	
	vec = vector<int>(count);
	
	if(count > 0){
		int *tmp = (int*)calloc(count, sizeof(int));
		
		if(tmp == NULL){
			throw __FILE__ "RecvIntVec: Unable to allocate memory for buffer";
		}
		
		MPI_Recv(tmp, count, MPI_INT, src, tag, comm, MPI_STATUS_IGNORE);
		
		for(int i = 0;i < count;i++){
			vec[i] = tmp[i];
		}
			
		free(tmp);
	}
}

void broadcastFile( string& b_filename ){
	ifstream b_file;
	char* b_buf;
	
	b_file.open( b_filename.c_str(), ios::binary );
	
	if( !b_file.is_open() ){
		cerr << "Error opening \"" << b_filename << "\"\n";
		throw __FILE__ "(MpiBlast::broadcastFile): Unable to open file";
	}
	
	b_file.seekg( 0, ios::end );
	
	uint b_filesize = b_file.tellg();
	
	if( debug_msg ){
		LOG_MSG << "broadcasting file size of " << b_filesize << endl;
	}
	
	MPI_Bcast( &b_filesize, 1, MPI_UNSIGNED, super_master_process, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "file size broadcasted\n";
	}
	
	//b_buf = new char[ b_filesize ];
	b_buf = (char*)calloc(b_filesize, sizeof(unsigned char));
	
	if(b_buf == NULL){
		throw __FILE__ "(MpiBlast::broadcastFile): Unable to allocate memory";
	}
	
	b_file.seekg( 0, ios::beg );
	
	if( b_file.read( b_buf, b_filesize ) ){
		uint b_read = b_file.gcount();
		if( b_read != b_filesize ){
		
			free(b_buf);
			
			cerr << "Only read " << b_read << " of " << b_filesize << " bytes\n";
			cerr << "Error completely reading \"" << b_filename << "\"\n";
			
			throw __FILE__ "(MpiBlast::broadcastFile): Error completely reading file";
		}
		
		if( debug_msg ){
			LOG_MSG << "broadcasting file\n";
		}
		
		MPI_Bcast( b_buf, b_read, MPI_BYTE, super_master_process, MPI_COMM_WORLD );
		
		if( debug_msg ){
			LOG_MSG << "file broadcasted\n";
		}
	}else{
		free(b_buf);
		
		cerr << "Error reading\"" << b_filename << "\"\n";
		
		throw __FILE__ "(MpiBlast::broadcastFile):Error reading file";
	}
	
	b_file.close();
	
	free(b_buf);
}


void recvBroadcastFile( string& b_filename ){
	ofstream b_file;
	char* b_buf;

	if( debug_msg ){
		LOG_MSG << "waiting for file size broadcast\n";
	}
	
	uint b_filesize;
	
	MPI_Bcast( &b_filesize, 1, MPI_UNSIGNED, super_master_process, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "received file size broadcast of " << b_filesize << endl;
	}
	
	//b_buf = new char[ b_filesize ];
	b_buf = (char*)calloc(b_filesize, sizeof(unsigned char));
	
	if(b_buf == NULL){
		throw __FILE__ "(MpiBlast::recvBroadcastFile): Unable to allocate memory";
	}
	
	if( debug_msg ){
		LOG_MSG << "opening receive file " << b_filename << endl;
	}
	
	b_file.open( b_filename.c_str() );

	if(!b_file.is_open()){
		free(b_buf);
		cerr << "Error opening file: " << b_filename.c_str() << endl;
		throw __FILE__ "(MpiBlast::recvBroadcastFile): Unable to open file";
	}
	
	if( debug_msg ){
		LOG_MSG << "receiving file to " << b_filename << endl;
	}
	
	MPI_Bcast( b_buf, b_filesize, MPI_BYTE, super_master_process, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "received file broadcast\n";
	}
	
	b_file.write( b_buf, b_filesize );
	b_file.close();
	
	free(b_buf);
}


uint64 	max_data_to_send = 1024*1024;	/**< max data to send when using MPI to copy fragments */

void sendFile( const string& filename, int dest, int tag ){
	ifstream file_input( filename.c_str() );
	if( !file_input.is_open() ){
		LOG_MSG << "Unable to open \"" << filename << "\"\n";
		throw "Unable to open file for send.";
	}
	
	// send the file size
	file_input.seekg( 0, ios::end );
	uint filesize = file_input.tellg();	
	MPI_Send( &filesize, 1, MPI_UNSIGNED, dest, tag, MPI_COMM_WORLD );
	file_input.seekg( 0, ios::beg );

	// read and send the file
	char* buffer;
	int bufsize = max_data_to_send;
	meta_MPI_Alloc_mem( bufsize, MPI_INFO_NULL, &buffer );

	while( true ){
		file_input.read( buffer, bufsize );
		int count = file_input.gcount();
		if( count <= 0 )
			break;
		MPI_Send( buffer, count, MPI_BYTE, dest, tag, MPI_COMM_WORLD );
	}
	meta_MPI_Free_mem( buffer );
}

void recvFile( const string& filename, int src, int tag ){
	ofstream output_file;
	char* buffer;
	int bufsize = max_data_to_send;
	
	output_file.open( filename.c_str() );
	if( !output_file.is_open() ){
		throw "Unable to open file";
	}

	uint b_filesize;
	
	MPI_Recv( &b_filesize, 1, MPI_UNSIGNED, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	
	meta_MPI_Alloc_mem( bufsize, MPI_INFO_NULL, &buffer );
	
	while( output_file.tellp() < b_filesize ){
		int remaining = b_filesize - output_file.tellp() < bufsize ? b_filesize - output_file.tellp() : bufsize;
		MPI_Recv( buffer, remaining, MPI_BYTE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		output_file.write( buffer, remaining );
	}
	output_file.close();
	meta_MPI_Free_mem( buffer );
}

void printProgress (uint64 progress_sofar, float64 total) {
	int temp_ratio = (int)((progress_sofar/total)*10000);
	float ratio = (float)temp_ratio / 100.0 ;
	cout << setprecision(2);
	cout << "\b\b\b\b\b" << ratio << "%" ;
	if (progress_sofar == total) 
		cout << "\b\b\b\b\b\b" ;
	cout.flush();
}

//Lock or unlock an entire file
#ifndef WIN32
int getlock (int fd, int cmd, int type) {
	struct flock lock;
	lock.l_type = type;  //F_RDLCK , F_WRLCK, F_UNLCK
	lock.l_start = 0;
	lock.l_whence = SEEK_SET;
	lock.l_len = 0;
	// only lock the file if locking was enabled--
	// some NFS systems don't behave well with locking
	// enabled so we only want to use it when absolutely
	// necessary
	if( lock_mbf )
		return( fcntl(fd,cmd,&lock));
	else
		return 0;
}
#endif

GroupManager::GroupManager(int group_size, int num_frags) {
	// initialize groups
	_group_size = group_size;

	std::vector <int> proc_pool;

	for(int proc_id=0; proc_id<node_count; proc_id++) {
		if(proc_id == super_master_process) {
			continue;
		}
		proc_pool.push_back(proc_id);
	}
	
	_num_groups = (node_count - 1) / _group_size;
	int remain = (node_count - 1) % _group_size;

	int i = 0;
	for(int group_id=0; group_id<_num_groups; group_id++) {
		_actual_group_size.push_back(0);
		for(int j=0; j<group_size; j++) {
			_proc_group[proc_pool[i++]] = group_id;
			_actual_group_size[group_id]++;
		}
	}
	
	if( remain > 0 ) {
		int last_group = -1; 

		if(remain >= num_frags) {
			_num_groups++;
			_actual_group_size.push_back(0);
		}

		last_group = _num_groups - 1;
		
		for(; i<proc_pool.size(); i++) {
			_proc_group[proc_pool[i]] = last_group;
			_actual_group_size[last_group]++;
		}
	}
}

void GroupManager::PrintGroupMap(std::ostream& os) {
	std::vector < std::vector <int> > group_procs(_num_groups, std::vector<int>());

	std::map <int, int>::iterator mit;
	for(mit=_proc_group.begin(); mit!=_proc_group.end(); mit++) {
		group_procs[mit->second].push_back(mit->first);
	}

	for(int i=0; i<_num_groups; i++) {
		os << "Partition " << i << " has " <<  _actual_group_size[i] << " processes: ";
		copy(group_procs[i].begin(), group_procs[i].end(), std::ostream_iterator<int>(os, " "));
		os << std::endl;
	}

	return;
}

