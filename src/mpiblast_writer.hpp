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
#ifndef __MPIBLAST_WRITER_HPP__
#define __MPIBLAST_WRITER_HPP__

#ifdef __cplusplus
#include <blast_wrapper.hpp>
#include <pio_intercept.h>
#include <mpiblast_util.hpp>
#include <mpiblast_querydata.hpp>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <ncbi.h>
#include <mpi.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define SHARED_FS	"NFS"
#endif

#define MAX_SINGLE_WRITE		2*1000*1024*1024

// output strategies
#define MASTER_STREAMLINE				0
#define WORKER_SPLIT					1
#define WORKER_INDIVIDUAL				2
#define WORKER_COLLECTIVE				3
#define WORKER_STREAMLINE				4

// io functions
#define MPI_IO_FUNC						0
#define POSIX_IO_FUNC					1

using namespace std;

// Manage a data buffer
class DataStruct {
public:
	DataStruct() : _size(0), _buffer(NULL) { _ptr = _buffer; }
	DataStruct(MPI_Offset buffer_size) : _size(buffer_size) {
		_buffer = (char*)malloc(_size);
		CHECK_NULPTR(_buffer);
		_ptr = _buffer;
	}
	~DataStruct() { 
		if(_buffer!=NULL) free(_buffer); 
	}
	
	// Append data to current location
	void AddData(const void* inptr, int data_size) {
#ifdef DEBUG
		if(_ptr - _buffer + data_size > _size) {
			throw __FILE__ "Buffer over flow";
		}
#endif
		memcpy(_ptr, inptr, data_size);
		_ptr += data_size;
	}

	// Extract data from current location
	void ExtractData(void* inptr, int data_size) {
#ifdef DEBUG
		if(_ptr - _buffer + data_size > _size) {
			throw __FILE__ "Buffer over flow";
		}
#endif
		memcpy(inptr, _ptr, data_size);
		_ptr += data_size;
	}
	
	void ResetDataBuf() {
		memset(_buffer, _size, 0);
		_ptr = _buffer;
	}

	// notice caller must take care of deallocation of the outptr
	void ReapData(char** outptr) {
		*outptr = _buffer;
		_buffer = NULL;
	}

	int GetSize() { return _size; }

	MPI_Request* GetReqAddress() { return &_req; }
	MPI_Request GetReq() { return _req; }

protected:
	MPI_Request _req; // request object, used by non-blocking communication or I/O
	MPI_Offset _size;		 // size of the buffer
	char* _ptr;  	// data processing pointer 
	char* _buffer;	// buffer
};

// send communication data structure
class CommSendStruct : public DataStruct {
public:
	CommSendStruct() : DataStruct() {}
	CommSendStruct(MPI_Offset buffer_size) : DataStruct(buffer_size) {}
	
	void SendData(MPI_Comm comm, int dst, int tag) {
		_dst = dst;
		MPI_Send(_buffer, _size, MPI_BYTE, dst, tag, comm);
	}

	void IsendData(MPI_Comm comm, int dst, int tag) {
		_dst = dst;
		MPI_Isend(_buffer, _size, MPI_BYTE, dst, tag, comm, &_req);
	}
	
	// test isend
	int TestIsend() {
		MPI_Status status;
		int flag;
		
		MPI_Test(&_req, &flag, &status);
		return flag;
	}

	void WaitIsend() {
		MPI_Status status;
		
		MPI_Wait(&_req, &status);
	}

	int GetDest() { return _dst; }
protected:
	int _dst;
};

// recv communication data structure
class CommRecvStruct : public DataStruct {
public:
	CommRecvStruct() : DataStruct() {}
	CommRecvStruct(MPI_Offset buffer_size) : DataStruct(buffer_size) {}
	CommRecvStruct(MPI_Offset buffer_size, int src, int tag, MPI_Comm comm) : DataStruct(buffer_size) {
		_comm = comm;
		_src = src;
		_tag = tag;
	}
	
	void RecvData(MPI_Comm comm, int src, int tag) {
		_src = src;

		double profile_time = MPI_Wtime();

		MPI_Status status;
		MPI_Recv(_buffer, _size, MPI_BYTE, src, tag, comm, &status);
	
		receive_msg_time += MPI_Wtime() - profile_time;
		receive_msg_size += _size;
	}

	void IrecvData(MPI_Comm comm, int src, int tag) {
		_src = src;
		_comm = comm;
		_tag = tag;
		MPI_Irecv(_buffer, _size, MPI_BYTE, src, tag, comm, &_req);
	}

	void IrecvData() {
		MPI_Irecv(_buffer, _size, MPI_BYTE, _src, _tag, _comm, &_req);
	}

	// test irecv
	int TestIrecv() {
		MPI_Status status;
		int flag;
		
		MPI_Test(&_req, &flag, &status);
		return flag;
	}

	void WaitIrecv() {
		MPI_Status status;
		
		MPI_Wait(&_req, &status);
	}

	void Cancel() { MPI_Cancel(&_req); }

	int GetSource() { return _src; }
	int GetTag() { return _tag; }
protected:
	int _src; // source of the receive call
	int _tag; // tag of the receive call
	MPI_Comm _comm;
};

// class to streamline output
class Streamliner;

// Wrapper class to manage OutputRecord
class COutputRecord {
public: 
	// If auto_clean is set to true, the internal data structure will be freed 
	// in the desctructor. The destructor won't do cleanup if auto_clean is set
	// to false. The purpose of this is to save memory copy operations. 
	// Different COutputRecord object could point to the same OutputRecordPtr, 
	// but only one one will free up the actual OutputRecord. 
	COutputRecord(OutputRecordPtr orp, bool auto_clean) {
		assert(orp);
		_orp = orp;
		_auto_clean = auto_clean;
	}
	
	// the copy always set auto_clean to true 
	COutputRecord(const COutputRecord& cor) {
		_orp = cor._orp;
		_auto_clean = true;
	}
	
	// used in calculating global output offsets 
	void UpdateDescriptionOffset(MPI_Offset& out_off, bool update_off=false) {
		if(_orp->des_size==0) {
			return;
		}
		_orp->des_offset = out_off;
		
		if(update_off) {
			// update output offset
			out_off += _orp->des_size;
		}
	}
	
	// used in calculating global output offsets 
	void UpdateAlignOffset(MPI_Offset& out_off, bool update_off=false) {
		if(_orp->aln_size==0) {
			return;
		}
		_orp->aln_offset = out_off;
		
		if(update_off) {
			// update output offset
			out_off += _orp->aln_size;
		}
	}
	
	// compare two records
	bool operator >= (const COutputRecord& cor) {
		return (OutputRecordCmp(_orp, cor._orp)>=0);
	}
	
	void DebugPrint(FILE *fp) {
		print_output_record(fp, _orp);
	}
	
	~COutputRecord() {
		// free up the left data
		if(_auto_clean && _orp != NULL) {
			if(_orp->des_data!=NULL) {
				free(_orp->des_data);
				_orp->des_data = NULL;
			}
			if(_orp->aln_data!=NULL) {
				free(_orp->aln_data);
				_orp->aln_data = NULL;
			}
			free(_orp);
			_orp = NULL;
		}
	}
	
	OutputRecordPtr _orp;
	bool _auto_clean;
};

// Class responsible for write buffered output to the file system
class FileWriter {
public:
	FileWriter(const string& output_file, bool split_write=false) {
		_file_name = output_file;
		_initialized = true;
		_split_write = split_write;
	}
	virtual ~FileWriter() {};
	
	void AddWriteEntry(MPI_Offset& out_off, string str_out_data, bool update_off=false) {
		int data_size = str_out_data.length();
		char* data_buf = (char*)malloc(data_size); // will be free after write
		CHECK_NULPTR(data_buf);
		memcpy(data_buf, str_out_data.c_str(), data_size);
		
		AddWriteEntry(out_off, data_size, data_buf, update_off);
	}
	
	void AddWriteEntry(MPI_Offset& out_off, int out_size, char* out_data, bool update_off=false) {
		_out_off_vec.push_back(out_off);
		_out_size_vec.push_back(out_size);
		_out_stuff_vec.push_back(out_data);

		if(update_off) {
			//update current output position
			out_off += out_size;
		}
	}
	
	virtual void Write();	
	virtual void CollWrite() {}	
	virtual void SingleWrite(int start_idx, int end_idx) = 0;
	virtual void SplitWrite(int start_idx, int end_idx) = 0;
	
	void DumpWriteRequests(FILE* fp);

	void DebugPrint(FILE* fp) {
		fprintf(fp, "\n\ndebug print\n");
		fprintf(fp, "_out_off_vec.size()=%d\n", _out_off_vec.size());
		fprintf(fp, "_out_size_vec.size()=%d\n", _out_size_vec.size());
		fprintf(fp, "_out_stuff_vec.size()=%d\n", _out_stuff_vec.size());
		fflush(fp);
		int num_elements = _out_size_vec.size();
		for(int i=0; i<num_elements; i++) {
			fprintf(fp, "off=%d\n", _out_off_vec[i]);
			fprintf(fp, "size=%d\n", _out_size_vec[i]);
			fflush(fp);
			fwrite(_out_stuff_vec[i], _out_size_vec[i], 1, fp);
			fflush(fp);
		}
	}

protected:
	string _file_name;
	vector <MPI_Offset> _out_off_vec; // output offsets
	vector <int> _out_size_vec; // output sizes
	vector <char*> _out_stuff_vec; // output stuff
	bool _initialized;
	bool _split_write;
};

class MPIFileWriter : public FileWriter {
public:
	MPIFileWriter(const string& output_file, bool shared=false, bool split_write=false, bool use_file_view=true) : FileWriter(output_file) {
		MPI_Comm comm;
		if(shared) {
			comm = group_comm;
		} else {
			comm = MPI_COMM_SELF;
		}

		_split_write = split_write;
		_use_file_view = use_file_view;

		if(strcasecmp(SHARED_FS, "PVFS") == 0) {
			MPI_Info info;
			
			MPI_Info_create(&info);
			MPI_Info_set(info, "romio_pvfs2_posix_write", "disable");
			MPI_Info_set(info, "romio_pvfs2_dtype_write", "enable");
			MPI_Info_set(info, "romio_cb_write", "disable");
			
			MPI_File_open(comm, (char*)_file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &_fh);
			MPI_Info_free(&info);
		} else {
			MPI_File_open(comm, (char*)_file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &_fh);
		}
#ifdef MPI_ATOMICITY
		MPI_File_set_atomicity(_fh, TRUE);
#endif
	}

	~MPIFileWriter() {
		if(_out_stuff_vec.size() > 0) {
			throw __FILE__ "There is unwritten data";
		}

		MPI_File_sync(_fh);
		MPI_File_close(&_fh);
	}
	
	void SingleWrite(int start_idx, int end_idx);
	void SplitWrite(int start_idx, int end_idx);
	
private:
	MPI_File _fh;
	bool _use_file_view;
};

class CollFileWriter : public FileWriter {
public:
	CollFileWriter(const string& output_file) : FileWriter(output_file) {
		MPI_File_open(group_write_comm, (char*)_file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &_fh);
	}

	~CollFileWriter() {
		if(_out_stuff_vec.size() > 0) {
			throw __FILE__ "There is unwritten data";
		}

		MPI_File_sync(_fh);
		MPI_File_close(&_fh);
	}
	
	void Write() {}
	void CollWrite();

	void SingleWrite(int start_idx, int end_idx) {}
	void SplitWrite(int start_idx, int end_idx) {}
	
private:
	MPI_File _fh;
	bool _use_file_view;
};

class PosixFileWriter : public FileWriter {
public:
	PosixFileWriter(const string& output_file, bool init, bool split_write=false, bool use_lock=true) : FileWriter(output_file) {
		if(init) {
#ifdef O_LARGEFILE
			_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT|O_LARGEFILE, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#else
			_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#endif
			if( _fd < 0 ) {
				perror("file open error:");
				throw __FILE__ "PosixFileWriter(): cannot open the output file";
			}	
		} else {
			_initialized = false;
		}
		
		_split_write = split_write;
		_use_lock = use_lock;
	}	

	~PosixFileWriter() {
		if(_out_stuff_vec.size() > 0) {
			throw __FILE__ "There is unwritten data";
		}

		if(_initialized) {
			close(_fd);
		}
	}
	
	void SingleWrite(int start_idx, int end_idx);
	void SplitWrite(int start_idx, int end_idx); 
	
private:
	int _fd;
	bool _use_lock;
};

// Output data of a query
class QueryOutput {
public:
	QueryOutput(int max_size);
	//QueryOutput(ValNodePtr orplist, QueryOutputInfoPtr qoip);
	~QueryOutput();
	
	void AddOutputRecordsToList(ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr, bool add_orp=true);
	void PruneOutputList();
	void CalculateOutputOffsets(BLAST* blastp, FileWriter* fwp, MPI_Offset& curr_output_pos);
	void AdjustOutputOffsets(MPI_Offset& curr_output_pos);
	void PrepareOutputData(BLAST* blastp, FileWriter* fwp);
	void UpdateInfo(QueryOutputInfoPtr qoip);
    void CleanupOutputs();
	
	// for debug
	void DebugPrint(FILE* fp);
	void DebugWriteAsn(FILE* fp);
	
	QueryOutputInfoPtr _qoip;
	unsigned char* _stats_buf;
	int _stats_buf_size;
	bool _is_ready;
    bool _is_offsets_sent;
	bool _outputs_freed;

	list < COutputRecord > _sorted_output_list;
	int _max_size;
	Streamliner* _stliner;

	int _ready_workers;
};

class IsendMgr {
	public:
		void AddPendingSend(CommSendStruct* cssp) { _send_ops_list.push_back(cssp); }
		void CheckPendingSends();
		void WaitPendingSends();
	protected:
		list <CommSendStruct*> _send_ops_list; // send operations list 
};

class IrecvMgr {
    public:
        ~IrecvMgr();
        void AddPendingRecv(CommRecvStruct* crsp) { _recv_ops_list.push_back(crsp); }
        int WaitAnyRecv(CommRecvStruct* &crsp, MPI_Status& status);
        int TestAnyRecv(CommRecvStruct* &crsp, MPI_Status& status);
		int GetNumPendingRecvs() { return _recv_ops_list.size(); }
    protected:
        list <CommRecvStruct*> _recv_ops_list;
};

// Base class for output management
class Writer {
public:
	Writer(BLAST* blast, int query_count, int total_frags) {
		_query_count = query_count;
		_total_frags = total_frags;
		_blast = blast;
		_offsets_to_send = 0;
		_last_output_query = 0;
		_last_written_query = -1;
//		_query_to_write = -1;
		_is_finished = false;
		_output_drp = NULL;
		_no_more_queries = false;
		_assigned_queries = 0;
		
	/*
		for(int i=0; i<query_count; i++) {
			// _output_vec.push_back(new QueryOutput(_blast->GetMaxShow()));
			_output_vec.push_back(NULL);
	*/		
	}
	
	void PrintOutputsList(FILE* fp);
	void AddOutputRecords(int query_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr, bool add_orp=true);
	void UpdateQueryOutputInfo(int query_id, QueryOutputInfoPtr qoip);
	void UpdateDefaultReportInfo(QueryOutputInfoPtr qoip);
	void PrepareDefaultReport() { _blast->GetDefaultReportInfo(); }
	void AddCommSend(CommSendStruct* cssp) { _comm_send.AddPendingSend(cssp); }
	void CheckCommSends() { _comm_send.CheckPendingSends(); }
	void WaitCommSends() { _comm_send.WaitPendingSends(); }
	void ProcessAsyncOutput();
	void AddWriteOps(int query_id, int num_write_segments) {
		for(int i=0; i<num_write_segments; i++) {
			_writing_ops.push_back(pair<int,int>(query_id, i));
		}
	}
	void WaitAllWriteOps();
	void SendWriteAck(int query_id, int dst);
	void AddWorkingQueries(int segment_id, int start_query, int end_query);
	void AddWorkingQueries(int segment_id, int start_query, int end_query, int writer_id);
	void AddWorkingQueries(int query_id);
	void InitQueryOutput(int segment_id, int query_id);
	void SetNoMoreQueries() { _no_more_queries = true; }
	void UnloadQuery(int query_id);
	int GetNumWorkingQueries() { return _working_queries.size(); }
	virtual int HandleMessages(MPI_Status &probe_status, int* event_array) = 0;
	virtual int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) = 0;
	virtual void Finalize() = 0;
	virtual ~Writer() {
		if(_output_drp != NULL) {
			free_default_reportp(_output_drp);
		}

		if(!_output_vec.empty()) {
			// if(debug_msg) {
			//	LOG_MSG << "WARNING: Some queries have not been unloaded, possibly they don't have output." << endl;
			// }
			// throw __FILE__ " Invalid writer output structure state";
		}
		if(!_pending_offsets.empty()) {
			throw __FILE__ " There is unprocessed results";
		}
		if(!_working_queries.empty()) {
			throw __FILE__ " There is still working queries left";

		}
	};

protected:
	// vector < QueryOutput* > _output_vec; // sorted output records list
	map < int, QueryOutput* > _output_vec; // sorted output records list
	
	queue <int> _lead_queries; // queries that are ready for output, currently queries are output in order of query_id
	list < pair<int, int> > _writing_ops;  // pending writing ops: query_id, segment_id
	
	BLAST* _blast;
	FileWriter* _fw;
	IsendMgr _comm_send; //common send ops
	DefaultReportPtr _output_drp;
	int _query_count;
	int _total_frags;
	int _offsets_to_send;
	int _last_output_query;
	int _last_written_query;
	int _assigned_queries;
	bool _is_finished;

	bool _no_more_queries;

	// vector <int> _working_queries; // queries that have been assigned to this group
	list <int> _working_queries; // queries that have been assigned to this group
	map <int, int> _query_map; // key: query_id, data: segment_id
	set <int> _pending_offsets;
};

// Master side output management class
class WriterMaster : public Writer {
public:
	WriterMaster(BLAST* blast, int query_count, int total_frags, int output_batch) 
		: Writer(blast, query_count, total_frags), _output_batch(output_batch),
		_curr_output_pos(0), _complete_worker(0), _num_writing(0), _last_sent_offset(0),
	    _num_written_queries(0)	{

			if(output_strategy ==  WORKER_COLLECTIVE) {
				_curr_output_pos += group_node_count;
			}
	
#ifdef COLL_SEARCH_PROF
			for(int i=0; i<query_count; i++) {
				_search_end_t.push_back(vector <double>(0, total_frags));
				_merge_t.push_back(vector <double>(0, total_frags));
			}
#endif

		}
	void RecvPartialOutputs(CommRecvStruct* crsp, int src);
	void CheckUnsentOffsets();
	void SendOutputOffsets(int start_query, int end_query);
	bool WriterWorkerComplete() { return _complete_worker == GetNumWorkers(); }
	void ProcessWriteRequest(int src);
	void ProcessWriteComplete(int src);
	void TellWriterWorkerQuit();
	virtual ~WriterMaster() {}
    virtual int CheckEmptyOutput() {
		throw __FILE__ "a prototype of virtual function, should not be call at this level";
	};

#ifdef COLL_SEARCH_PROF
	void DumpSearchProfile(string& file_name);
#endif

protected:
	map < int, QueryData* > query_data;
	queue <int> _write_request_queue;
	
	int _output_batch;
	int _complete_worker;
	int _num_writing;
	int _last_sent_offset;
	int _num_written_queries;
	
	MPI_Offset _curr_output_pos;  // current global output pos

#ifdef COLL_SEARCH_PROF
	vector < vector <double> > _search_end_t;
	vector < vector <double> > _merge_t;
#endif
};

/* WriterMaster with individual output strategy */
/* make it a singleton class */
class WriterMasterIndividual : public WriterMaster {
public:
	static WriterMasterIndividual* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int output_batch) { return NULL; }
	static WriterMasterIndividual* Instance(void) { return NULL; }

	int HandleMessages(MPI_Status &probe_status, int* event_array) { return 0; }
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) { return 0; }
	void Finalize(void) { return; }
	
protected:
	WriterMasterIndividual(BLAST *blast, int query_count, int total_frags, const string& output_file, int output_batch=1) 
		: WriterMaster(blast, query_count, total_frags, output_batch) {}
	
	WriterMasterIndividual &operator=(WriterMasterIndividual&);

private:
	static WriterMasterIndividual* _instance;
};

class WriterMasterCollective : public WriterMaster {
public:
	static WriterMasterCollective* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int output_batch) { return NULL; }
	static WriterMasterCollective* Instance(void) { return NULL; }

	int HandleMessages(MPI_Status &probe_status, int* event_array) { return 0; }
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) { return 0; }
	void ProcessCollRequest(CommRecvStruct* crsp) {}
	void Finalize(void) {} 
	
protected:
	WriterMasterCollective(BLAST *blast, int query_count, int total_frags, const string& output_file, int output_batch=1) 
		: WriterMaster(blast, query_count, total_frags, output_batch) {}
	
	WriterMasterCollective &operator=(WriterMasterCollective&) {}

private:
	static WriterMasterCollective* _instance;
};

/* WriterMaster with streamline output strategies */
/* make it a singleton class */
class WriterMasterStreamline : public WriterMaster {
public:
	static WriterMasterStreamline* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int output_batch);
	static WriterMasterStreamline* Instance(void);

	int HandleMessages(MPI_Status &probe_status, int* event_array);
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status);
	void RecvMetaDataToWrite(CommRecvStruct* crsp, int src, FileWriter* fwp);
	void FinishQuery(CommRecvStruct* crsp);
	void ProcessWaitWriteNotification(CommRecvStruct* crsp);
	virtual int CheckEmptyOutput();
	void Finalize(void);
	
protected:
	WriterMasterStreamline(BLAST *blast, int query_count, int total_frags, const string& output_file, int output_batch=1) 
		: WriterMaster(blast, query_count, total_frags, output_batch) {
		// install FileWriter in this level so master and worker can have different write strategy
		if( io_function == POSIX_IO_FUNC) {
			if(output_strategy == MASTER_STREAMLINE) {
				_fw = new PosixFileWriter(output_file, true, false, !disable_posix_lock);
			} else if (output_strategy == WORKER_STREAMLINE) {
				_fw = new PosixFileWriter(output_file, false, false, !disable_posix_lock);
			}
		} else if (io_function == MPI_IO_FUNC) {
			_fw = new MPIFileWriter(output_file, true, false, false);
		}

		// output html header
		char* buffer = (char*) malloc (PRINT_BUFFER_LENGTH);
		CHECK_NULPTR(buffer);
		memset(buffer, 0, PRINT_BUFFER_LENGTH);
		_blast->GetHtmlHeader(buffer);
		int print_len = strlen(buffer);
		if(print_len>0) {
			_fw->AddWriteEntry(_curr_output_pos, print_len, buffer, true); // buffer will be freed by FileWriter
			_fw->Write();
		} else {
			free(buffer);
		}
	}
	
	WriterMasterStreamline &operator=(WriterMasterStreamline&);

private:
	static WriterMasterStreamline* _instance;
};

// Worker side output management
class WriterWorker : public Writer {
public:
	WriterWorker(BLAST* blast, int query_count, int total_frags, int send_batch) 
		: Writer(blast, query_count, total_frags), _tmp_output_vec(query_count), 
		_queries_to_processed(0), _last_processed_query(0), _send_batch(send_batch) {}
	void AddLocalOutputs(int query_id, int frag_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr);
	void SetStartQuery(int start_query); // set the start query for this run
	void SendPartialOutputs(int start_query, int end_query, int pack_type, int frag_id);
	void ProcessResults(bool is_search_end, int frag_id);
	void ReceiveOutputOffsets(CommRecvStruct* crsp);
	void CheckResultSends() { _results_send.CheckPendingSends(); }
//	virtual void ProcessResults(bool is_search_end, int frag_id) {
//		throw __FILE__ "a prototype of virtual function, should not be call at this level";
//	};
	virtual ~WriterWorker() {}; // singleton won't work without definition of destructor
	int GetNumPendingOffsets() { return _pending_offsets.size(); }
	
protected:
	vector <vector <OutputRecordPtr> > _tmp_output_vec; // vector to store tempory results

	IsendMgr _results_send;
	int _queries_to_processed; // number of queries having results to be processed 
	int _last_processed_query; // the last query whose results have been processed
	int _send_batch;  // number of query related data that will be sent in group
};

/* WriterWorker with individual output strategy */
/* make it a singleton class */
class WriterWorkerIndividual : public WriterWorker {
public:
	static WriterWorkerIndividual* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int send_batch) {} 
	static WriterWorkerIndividual* Instance(void) {}
	void Finalize(void) {}
	
	int HandleMessages(MPI_Status &event_status, int* event_array);
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) {}
	
protected:
	WriterWorkerIndividual(BLAST *blast, int query_count, int total_frags, const string& output_file, int send_batch=1) 
		: WriterWorker(blast, query_count, total_frags, send_batch) {}
	
	WriterWorkerIndividual (WriterWorkerIndividual&);
	WriterWorkerIndividual &operator=(WriterWorkerIndividual&);
private:
	static WriterWorkerIndividual* _instance;
};

/* WriterWorker with collective output strategy */
/* make it a singleton class */
class WriterWorkerCollective : public WriterWorker {
public:
	static WriterWorkerCollective* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int send_batch) { return NULL; }
	static WriterWorkerCollective* Instance(void) { return NULL; }
	void Finalize(void);
	
	int HandleMessages(MPI_Status &event_status, int* event_array);
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) {}
	
protected:
	WriterWorkerCollective(BLAST *blast, int query_count, int total_frags, const string& output_file, int send_batch=1) 
		: WriterWorker(blast, query_count, total_frags, send_batch) {}
	
	WriterWorkerCollective (WriterWorkerCollective&);
	WriterWorkerCollective &operator=(WriterWorkerCollective&) {}
private:
	static WriterWorkerCollective* _instance;
};

/* WriterWorker with streamline output strategy */
/* make it a singleton class */
class WriterWorkerStreamline : public WriterWorker {
public:
	static WriterWorkerStreamline* Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int send_batch);
	static WriterWorkerStreamline* Instance(void);
	void Finalize(void);
	
	int HandleMessages(MPI_Status &event_status, int* event_array);
	int HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) {}
	void ReceiveLeaderData(CommRecvStruct* crsp);
	void RecvPeerMetaData(CommRecvStruct* crsp, int src);
	void OutputReadyQueries();
	//void ProcessAsyncOutput();
	
protected:
	WriterWorkerStreamline(BLAST *blast, int query_count, int total_frags, const string& output_file, int send_batch=1) 
		: WriterWorker(blast, query_count, total_frags, send_batch){
		// install FileWriter in this level so master and worker can have different write strategy
		if( io_function == POSIX_IO_FUNC ) {
			_fw = new PosixFileWriter(output_file, false, false, !disable_posix_lock);
		} else if ( io_function == MPI_IO_FUNC ) {
			_fw = new MPIFileWriter(output_file, true, false, false);
		}
	}
	
	WriterWorkerStreamline (WriterWorkerStreamline&);
	WriterWorkerStreamline &operator=(WriterWorkerStreamline&);
private:
	static WriterWorkerStreamline* _instance;
};

inline void Writer::UpdateQueryOutputInfo(int query_id, QueryOutputInfoPtr qoip) {
	_output_vec[query_id]->UpdateInfo(qoip);
}

inline void Writer::AddOutputRecords(int query_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr, bool add_orp) {
	_output_vec[query_id]->AddOutputRecordsToList(orplist, qoip, bs_ptr, add_orp);

	// HL-debug-open:
	// if(dump_raw_output && group_rank == 0) {
	//	_output_vec[query_id]->DebugPrint(dbgfp);
	//}
}

class OutputEntry {
	public:
		OutputEntry(MPI_Offset out_off, int out_size, char* out_data) {
			_out_off = out_off;
			_out_size = out_size;
			_out_data = out_data;
		}
		MPI_Offset _out_off;
		int _out_size; // output sizes
		char* _out_data; // output stuff
};

// Streamline the output data to a single processor for serialized writing
class Streamliner {
public:
	Streamliner(int amax_write) {
		_max_write = amax_write;
		_is_written = false;
		_is_processed = false;
		_is_ready = false;
		_is_empty = false;
		_assigned_leader = false;
		_write_leader = -1;
		_output_len = 0;
		_begin_pos = 0;
		_curr_write_segment = -1;
		_num_write_segments = -1;
		_complete_irecv = -1;
		_query_id = -1;
	}
	
	~Streamliner() {
		map <MPI_Offset, OutputEntry*>::iterator mit;
		for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
			delete mit->second;
		}
		_output_entries.clear();
		_submit_workers.clear();
	}


	MPI_Offset AddOutputEntry(MPI_Offset& out_off, int out_size, char* out_data, bool update_off=false, bool mark_write=false);
	void SetLastWritePos (MPI_Offset pos);
	MPI_Offset GetOutputLength() { return _output_len; }
	MPI_Offset GetBeginPos() { return _begin_pos; }
	void AppendWriteLeader(vector <MPI_Offset> &offset_buf_vec);
	void SendOutputData(int dst, int query_id, Writer* wtrp);
	void ExtractOutputEntries(CommRecvStruct* crsp);
	bool IsWritten();
	bool IsReady();
	bool IsProcessed();
	void AddWorkerWOutputs(int worker) { _workers_with_output.insert(worker); }
	void AddWorkerWOutputs(const set<int>& workers);
	void ExtractWritePeers(CommRecvStruct* crsp);
	void ExtractSubmitters(CommRecvStruct* crsp);
	void GatherOutputs(FileWriter* fwp);
	void GatherOutputsForAWrite(FileWriter* fwp);
	void WriteCurrentSegment(FileWriter *fwp);
	void WriteASegment(int segment_id, FileWriter *fwp);
	void WriteOutputEntries(FileWriter *fwp, MPI_Offset write_end);
	void SendDataForAWrite(int dst, int query_id, int segment_id);
	void SendDataForAWriteEx(int dst, int query_id, Writer* wtrp, MPI_Offset write_end, vector<CommSendStruct*>& msg_send_list);
	int AddPendingRecvs(int src);
	void CheckEmpty();
	bool IsEmpty() { return _is_empty; }
	void SetOutputLength(MPI_Offset output_len) { 
		_output_len = output_len; 
		MPI_Offset num_write_segments = _output_len / _max_write;
		MPI_Offset remain = _output_len % _max_write;
		if(remain > 0) {
			num_write_segments += 1;
		}

		_num_write_segments = num_write_segments;
		_num_ready_workers = vector <int> (_num_write_segments, 0);
	} 
	void SetBeginPos(MPI_Offset begin_pos) { _begin_pos = begin_pos; } ; 
	void ShiftLeader(int query_id, Writer* wtrp);
	void SetWriterLeader(int leader) { _write_leader = leader; }
	int GetWriterLeader() { return _write_leader; }
	void AddSearchedWorker(int worker) { _submit_workers.insert(worker); }
	bool HasWorkerResult(int worker);
	void SetAssignLeader() { _assigned_leader = true; }
	bool IsLeader() { return _assigned_leader; }
	bool LeaderCheckReady() {
		if(_assigned_leader && 
			_pend_recv_nodes.size() == _workers_with_output.size()) {
			_is_ready = true;
		}

		return _is_ready;
	}
	void PrepareCommWriteData();
	void AddMessageSizes(int src, const vector<int>& msg_sizes) {
		_data_msg_sizes[src] =  msg_sizes;
/*
		// debug !!!! just compare sizes
		if(msg_sizes.size() != _data_msg_sizes[src].size()) {
			throw __FILE__ "Message size mismatched";
		}

		for(int i=0; i<msg_sizes.size(); i++) {
			if(msg_sizes[i] != _data_msg_sizes[src][i]) {
				throw __FILE__ "Message size mismatched";
			}
		}
*/
	}

	void InitAWrite(int query_id, int segment_id);
	int CheckAWrite(int segment_id, FileWriter* fwp);
	int WaitAWrite(int segment_id, FileWriter* fwp);

	void AlignOutputEntries();
	void InitAccessLocation(int worker) {
		_access_locations[worker] = list< pair<MPI_Offset, int> >();
	}
	void AddAccessLocation(int worker, MPI_Offset offset, int size) {
		_access_locations[worker].push_back(pair<MPI_Offset, int>(offset, size));
	}
	void CalculateMessageSizes();
	int GetNumWriteSegments() { return _num_write_segments; }
	int HandleWaitWriteNotify(int write_segment_id); 
	bool PeersReady( int write_segment_id ) { 
		return _num_ready_workers[write_segment_id] == _workers_with_output.size(); 
	}

	// debug
	void DumpOutputEntries(FILE* fp);
	void DumpMessageSizes(FILE* fp);

private:
	
	map <MPI_Offset, OutputEntry*> _output_entries;
	map <int, list< pair <MPI_Offset, int> > > _access_locations;
	map <int, vector <int> > _data_msg_sizes; // key: rank, data: data message sizes

	set <int> _submit_workers; // workers that have submit results for this query
	set <int> _workers_with_output; // workers that have output data for this query
	
	vector < list <CommRecvStruct*> > _irecv_list; // recv operations list 
	vector < list <CommSendStruct*> > _isend_list;
	vector < int > _write_status; // indicuate write status of each write op: -1 - uninitialized; 0 - issued isend/irecv; 1 - complete isend/irecv

	vector <int> _pend_recv_nodes;
	vector <int> _num_ready_workers; // number of ready workers for a write segment
	//MPI_Offset _last_write_pos;
	MPI_Offset _output_len; // the output length for this query
	MPI_Offset _begin_pos; // the begin output position for this query
	int _curr_write_segment; // current write segment for which gathering results
	int _num_write_segments; // number of write segments
	int _complete_irecv; // num of completed irecv for current write segment
	MPI_Offset _max_write; // max length of each write
	bool _is_written;
	bool _is_processed;
	bool _is_ready;
	bool _is_empty;
	bool _assigned_leader;
	int _write_leader;
	int _query_id;
};

extern "C" {
#endif // __cplusplus

// Wrapper functions for c caller
void cProcessLocalOutputs(int query_id, int frag_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr);
void cUpdateQueryOutputInfo(int query_id, QueryOutputInfoPtr qoip);
void cUpdateDefaultReportInfo(QueryOutputInfoPtr qoip);

#ifdef __cplusplus
} // extern "C"
#endif	// __cplusplus

#endif	// __MPIBLAST_WRITER_HPP__
