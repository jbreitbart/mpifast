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
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <mpiblast_tags.hpp>
#include <mpiblast_writer.hpp>
#include <mpiblast.hpp>
#include <query_manager.hpp>
#include <iostream>
#include <iterator>
using namespace std;

// get send size for a output record according to the send type 
static int get_send_record_size(int send_type, OutputRecordPtr orp) {
	int i;
	int bit;
	int record_size=0;
	
	bit = RANK;
	if(bit & send_type) { /* send rank */
		record_size += int_size;
	}
	
	bit = FRAG_ID;
	if(bit & send_type) { /* send frag_id */
		record_size += int_size;
	}
	
	bit = ID;
	if(bit & send_type) { /* send seq id */
		/* not used currently */
	}
	
	bit = BIT_SCORE;
	if(bit & send_type) { /* send bit_score */
		record_size += float_size;
	}
	
	bit = EVALUE;
	if(bit & send_type) { /* send e_value */
		record_size += float_size;
	}
	
	bit = DES_SIZE;
	if(bit & send_type) { /* send des_size */
		record_size += int_size;
	}
	
	bit = DES_DATA;
	if(bit & send_type) { /* send des_data */
		record_size += orp->des_size;
	}
	
	bit = ALN_SIZE;
	if(bit & send_type) { /* send aln_size */
		record_size += int_size;
	}
	
	bit = ALN_DATA;
	if(bit & send_type) { /* send aln_data */
		record_size += orp->aln_size;
	}
	
	bit = DES_OFFSET;
	if(bit & send_type) { /* send des_offset */
		record_size += offset_size;
	}

	bit = ALN_OFFSET;
	if(bit & send_type) { /* send algin_offset */
		record_size += offset_size;
	}

	return record_size;	
}

// pack output record into the send buffer 
static int pack_record(int send_type, OutputRecordPtr orp, CommSendStruct* cssp) {
	int i;
	int bit;
	int record_size=0;
	
	bit = RANK;
	if(bit & send_type) { /* send rank */
		cssp->AddData(&(orp->rank), int_size);
	}
	
	bit = FRAG_ID;
	if(bit & send_type) { /* send frag_id */
		cssp->AddData(&(orp->frag_id), int_size);
	}
	
	bit = ID;
	if(bit & send_type) { /* send seq id */
		/* not used currently */
	}
	
	bit = BIT_SCORE;
	if(bit & send_type) { /* send bit_score */
		cssp->AddData(&(orp->bit_score), float_size);
	}
	
	bit = EVALUE;
	if(bit & send_type) { /* send e_value */
		cssp->AddData(&(orp->evalue), float_size);
	}
	
	bit = DES_SIZE;
	if(bit & send_type) { /* send des_size */
		cssp->AddData(&(orp->des_size), int_size);
	}
	
	bit = DES_DATA;
	if(orp->des_size > 0) {
		if(bit & send_type) { /* send des_data */
			cssp->AddData(orp->des_data, orp->des_size);
		}
	}
	
	bit = ALN_SIZE;
	if(bit & send_type) { /* send aln_size */
		cssp->AddData(&(orp->aln_size), int_size);
	}
	
	bit = ALN_DATA;
	if(orp->aln_size > 0) {
		if(bit & send_type) { /* send aln_data */
			cssp->AddData(orp->aln_data, orp->aln_size);
		}
	}
	
	bit = DES_OFFSET;
	if(bit & send_type) { /* send des_offset */
		cssp->AddData(&(orp->des_offset), offset_size);
	}

	bit = ALN_OFFSET;
	if(bit & send_type) { /* send algin_offset */
		cssp->AddData(&(orp->aln_offset), offset_size);
	}
	
	return 0;
}

// unpack output record from the receiving buffer 
static int unpack_record(int send_type, OutputRecordPtr orp, CommRecvStruct* crsp) {
	int i;
	int bit;
	int record_size=0;
	
	bit = RANK;
	if(bit & send_type) { /* send rank */
		crsp->ExtractData(&(orp->rank), int_size);
	}
	
	bit = FRAG_ID;
	if(bit & send_type) { /* send frag_id */
		crsp->ExtractData(&(orp->frag_id), int_size);
	}
	
	bit = ID;
	if(bit & send_type) { /* send seq id */
		/* not used currently */
	}
	
	bit = BIT_SCORE;
	if(bit & send_type) { /* send bit_score */
		crsp->ExtractData(&(orp->bit_score), float_size);
	}
	
	bit = EVALUE;
	if(bit & send_type) { /* send e_value */
		crsp->ExtractData(&(orp->evalue), float_size);
	}
	
	bit = DES_SIZE;
	if(bit & send_type) { /* send des_size */
		crsp->ExtractData(&(orp->des_size), int_size);
	}
	
	bit = DES_DATA;
	if(bit & send_type) { /* send des_data */
		if(orp->des_size > 0) {
			orp->des_data = (char *)malloc(orp->des_size);
			CHECK_NULPTR(orp->des_data);
			crsp->ExtractData(orp->des_data, orp->des_size);
		}
	}
	
	bit = ALN_SIZE;
	if(bit & send_type) { /* send aln_size */
		crsp->ExtractData(&(orp->aln_size), int_size);
	}
	
	bit = ALN_DATA;
	if(bit & send_type) { /* send aln_data */
		if(orp->aln_size > 0) {
			orp->aln_data = (char *)malloc(orp->aln_size);
			CHECK_NULPTR(orp->aln_data);
			crsp->ExtractData(orp->aln_data, orp->aln_size);
		}
	}
	bit = DES_OFFSET;
	if(bit & send_type) { /* send des_offset */
		crsp->ExtractData(&(orp->des_offset), offset_size);
	}

	bit = ALN_OFFSET;
	if(bit & send_type) { /* send algin_offset */
		crsp->ExtractData(&(orp->aln_offset), offset_size);
	}

	// HL-debug-open:
	if(debug_msg) {
		if(dump_raw_output) {
			print_output_record(dbgfp, orp);
		}
		// if(orp->rank < 0 || orp->rank > node_count) {
		if(orp->rank != crsp->GetSource()) {
			LOG_MSG << "Unpacking records: invalid worker rank" << orp->rank << endl;	
			throw __FILE__ "Error unpacking records";
		}
	}
	
	return 0;
}
	
// write all buffered data
void FileWriter::Write() {

	double write_track_time = 0;

	write_track_time = MPI_Wtime();
	
	if(_out_off_vec.size() == 0) {
		return;
	}

	// add writer synchronization here, notice it only limit the number of concurrent worker write, master can write at any time
	if(concurrent_write > 0 && group_rank != writer_process) {

		if(debug_msg) {
			LOG_MSG << "Send write request to writer master" << endl;
		}

		// ask master for write permission
		int event_array[2];
		event_array[0] = WRITE_REQUEST;
		event_array[1] = int_size;
		CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
		cssp_event->AddData(&event_array[0], 2 * int_size);
		cssp_event->SendData(group_comm, writer_process, EVENT_TYPE);
		delete cssp_event;

		// only need a tag here, the content is not important
		MPI_Send(event_array, 1, MPI_INT, writer_process, WRITE_REQUEST, group_comm);

		// wait for write permission
		int command;
		MPI_Status status;
		MPI_Recv(&command, 1, MPI_INT, writer_process, WORKER_WRITE, group_comm, &status);

		if(debug_msg) {
			LOG_MSG << "Receive write permission from writer master" << endl;
		}
	}

	if(_split_write) {
		SplitWrite(0, _out_off_vec.size());
	} else {
		SingleWrite(0, _out_off_vec.size());
	}

	if(concurrent_write > 0 && group_rank != writer_process) {
		// ask master for write permission
		int event_array[2];
		event_array[0] = WRITE_COMPLETE;
		event_array[1] = int_size;
		CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
		cssp_event->AddData(&event_array[0], 2*int_size);
		cssp_event->SendData(group_comm, writer_process, EVENT_TYPE);
		delete cssp_event;

		if(debug_msg) {
			LOG_MSG << "Send write completion to writer master" << endl;
		}

		// only need a tag here, the content is not important
		MPI_Send(event_array, 1, MPI_INT, writer_process, WRITE_COMPLETE, group_comm);
	}
	
	_out_off_vec.clear();
	_out_size_vec.clear();
	_out_stuff_vec.clear();

	write_time += MPI_Wtime() - write_track_time;
}

void FileWriter::DumpWriteRequests(FILE* fp) {
	// number of elements
	int num_elements = _out_off_vec.size();

	fwrite(&num_elements, sizeof(int), 1, fp);
	
	if(num_elements == 0) {
		return;
	}

	MPI_Offset buf_size = 0;
	// dumping sizes
	for(int i=0; i<num_elements; i++) {
		fwrite(&(_out_size_vec[i]), sizeof(int), 1, fp);
		buf_size += _out_size_vec[i];
	}

	// dumping offsets
	for(int i=0; i<num_elements; i++) {
		fwrite(&(_out_off_vec[i]), sizeof(MPI_Offset), 1, fp);
	}

	// buffer size
	fwrite(&buf_size, sizeof(MPI_Offset), 1, fp);

	// dumping writing stuff
	for(int i=0; i<num_elements; i++) {
		fwrite(_out_stuff_vec[i], sizeof(char), _out_size_vec[i], fp);
	}

	fflush(fp);
}

// perform a single mpi write operation
void MPIFileWriter::SingleWrite(int start_idx, int end_idx) {
	// construct file view
	int *elmoffset = 0;
	int *elmsize = 0;
	int num_elements = 0;
	static int int_size = sizeof(int);

	if(debug_msg) {
		LOG_MSG << "Start SingleWrite: from index " << start_idx <<" to "<<end_idx<<endl;
	}

	num_elements = end_idx - start_idx;
	if(num_elements == 0) {
		return;
	}
	
	elmoffset=(int *)calloc(num_elements, int_size);
	CHECK_NULPTR(elmoffset);
	elmsize=(int *)calloc(num_elements, int_size);
	CHECK_NULPTR(elmsize);
	
	// calculate buffer size;
	MPI_Offset buf_size = 0;
	for(int i=start_idx; i<end_idx; i++) {
		buf_size += _out_size_vec[i];
	}

	DataStruct* dsp = new DataStruct(buf_size);
	
	// adjust offsets based on the displacement
	MPI_Offset disp = _out_off_vec[start_idx];
		
	for(int i=0; i<num_elements; i++) {
		elmoffset[i] = _out_off_vec[start_idx + i] - disp;
		elmsize[i] = _out_size_vec[start_idx + i];
		dsp->AddData(_out_stuff_vec[start_idx + i], _out_size_vec[start_idx + i]);
		free(_out_stuff_vec[start_idx + i]);
	}
	
	char* data_2write;

	// remember to free data_2write
	dsp->ReapData(&data_2write);

	if(_use_file_view) {
		MPI_Info info;
		MPI_Datatype file_type;
		MPI_Type_indexed(num_elements, elmsize, elmoffset, MPI_BYTE, &file_type);
		if(strcasecmp(SHARED_FS, "PVFS") == 0) {
			MPI_Info_create(&info);
			MPI_Info_set(info, "romio_pvfs2_posix_write", "disable");
			MPI_Info_set(info, "romio_pvfs2_dtype_write", "enable");
			MPI_Info_set(info, "romio_cb_write", "disable");
			MPI_File_set_view(_fh, disp, MPI_BYTE, file_type, "native", info);
			MPI_Info_free(&info);
		} else {
			MPI_File_set_view(_fh, disp, MPI_BYTE, file_type, "native", MPI_INFO_NULL);
		}

		MPI_File_write(_fh, data_2write, buf_size, MPI_BYTE, MPI_STATUS_IGNORE);
		MPI_Type_free(&file_type);
	} else {
		MPI_File_write_at(_fh, disp, data_2write, buf_size, MPI_BYTE, MPI_STATUS_IGNORE);
	}

	free(data_2write);
	delete dsp;

	// free related data structure
	free(elmoffset);
	free(elmsize);

	if(debug_msg) {
		LOG_MSG << "End SingleWrite" <<endl;
	}
}

void MPIFileWriter::SplitWrite(int start_idx, int end_idx) {
	int num_elements = 0;
	static int int_size = sizeof(int);

	if(debug_msg) {
		LOG_MSG << "Start SplitWrite: from index " << start_idx <<" to "<<end_idx<<endl;
	}

	num_elements = end_idx - start_idx;
	if(num_elements == 0) {
		return;
	}

	for(int i=start_idx; i<end_idx; i++) {
		MPI_File_write_at(_fh, _out_off_vec[i], _out_stuff_vec[i], _out_size_vec[i], MPI_BYTE, MPI_STATUS_IGNORE);
		free(_out_stuff_vec[i]);
	}

	if(debug_msg) {
		LOG_MSG << "End SplitWrite" <<endl;
	}
}

// perform a single posix write operation
void PosixFileWriter::SingleWrite(int start_idx, int end_idx) {
	int num_elements = 0;
	static int int_size = sizeof(int);

	if(debug_msg) {
		LOG_MSG << "Start SingleWrite: from index " << start_idx <<" to "<<end_idx<<endl;
	}

	num_elements = end_idx - start_idx;
	if(num_elements == 0) {
		return;
	}

	if(!_initialized) {
#ifdef O_LARGEFILE
		_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT|O_LARGEFILE, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#else
		_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#endif
		if( _fd < 0 ) {
			perror("file open error:");
			throw __FILE__ "PosixFileWriter(): cannot open the output file";
		}	
	}

	// calculate buffer size;
	MPI_Offset buf_size = 0;
	for(int i=start_idx; i<end_idx; i++) {
		buf_size += _out_size_vec[i];
	}
	off_t disp = _out_off_vec[start_idx];
	
	if(_use_lock) {
		Posix_Set_lock(_fd, F_SETLKW, F_WRLCK, disp, SEEK_SET, buf_size);
	}

	lseek(_fd, disp, SEEK_SET);
	
/* !!! seems write_v is causing some problems on certain platforms, use write as more portable solution for now
	struct iovec iov[num_elements];
		
	for(int i=0; i<num_elements; i++) {
		iov[i].iov_base = _out_stuff_vec[start_idx + i];
		iov[i].iov_len = _out_size_vec[start_idx + i];
	}

	int write_len = writev(_fd, iov, num_elements);
	if (write_len != buf_size) {
		perror("writev error:");
		throw __FILE__ "write incomplet";
	}
*/

	DataStruct* dsp = new DataStruct(buf_size);
	
	for(int i=0; i<num_elements; i++) {
		dsp->AddData(_out_stuff_vec[start_idx + i], _out_size_vec[start_idx + i]);
		free(_out_stuff_vec[start_idx + i]);
	}
	
	char* data_2write;

	// remember to free data_2write
	dsp->ReapData(&data_2write);

	int write_len = write(_fd, data_2write, buf_size);
	if (write_len != buf_size) {
		perror("write error:");
		throw __FILE__ "write incomplet";
	}

	free(data_2write);
	delete dsp;

	if(_use_lock) {
		Posix_Set_lock(_fd, F_SETLK, F_UNLCK, disp, SEEK_SET, buf_size);
	}

/*
	// clean up buffered output
	for(int i=0; i<num_elements; i++) {
		free(_out_stuff_vec[start_idx + i]);
	}
*/

	if(!_initialized) {
		close(_fd);
	}

	if(debug_msg) {
		LOG_MSG << "End SingleWrite" <<endl;
	}
}

void PosixFileWriter::SplitWrite(int start_idx, int end_idx) {
	int num_elements = 0;
	static int int_size = sizeof(int);

	if(debug_msg) {
		LOG_MSG << "Start SplitWrite: from index " << start_idx <<" to "<<end_idx<<endl;
	}

	num_elements = end_idx - start_idx;
	if(num_elements == 0) {
		return;
	}

	// open file if the file is not open
	if(!_initialized) {
#ifdef O_LARGEFILE
		_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT|O_LARGEFILE, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#else
		_fd = open(_file_name.c_str(), O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#endif

		if( _fd < 0 ) {
			perror("file open error:");
			throw __FILE__ "PosixFileWriter(): cannot open the output file";
		}	
	}
	
	for(int i=start_idx; i<end_idx; i++) {
		off_t disp = _out_off_vec[i];

		if(_use_lock) {
			Posix_Set_lock(_fd, F_SETLKW, F_WRLCK, disp, SEEK_SET, _out_off_vec[i]);
		}
		lseek(_fd, disp, SEEK_SET);

		int write_len = write(_fd, _out_stuff_vec[i], _out_size_vec[i]);
		if (write_len != _out_size_vec[i]) {
			perror("write error:");
			throw __FILE__ "write incomplet";
		}

		if(_use_lock) {
			Posix_Set_lock(_fd, F_SETLK, F_UNLCK, disp, SEEK_SET, _out_off_vec[i]);
		}
		free(_out_stuff_vec[i]);
	}

	if(!_initialized) {
		close(_fd);
	}
	
	if(debug_msg) {
		LOG_MSG << "End SplitWrite" <<endl;
	}
}

void CollFileWriter::CollWrite() {

	double write_track_time = MPI_Wtime();

	int start_idx = 0;
	int end_idx = _out_off_vec.size();

	int *elmoffset = 0;
	int *elmsize = 0;

	int num_elements = end_idx - start_idx;

	if(num_elements == 0) {
		char buf[1];
		buf[0] = ' ';
		MPI_Offset disp = group_rank;
		MPI_File_set_view(_fh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		MPI_File_write_all(_fh, buf, 1, MPI_BYTE, MPI_STATUS_IGNORE);
	} else {
		elmoffset=(int *)calloc(num_elements, int_size);
		CHECK_NULPTR(elmoffset);
		elmsize=(int *)calloc(num_elements, int_size);
		CHECK_NULPTR(elmsize);
		
		// calculate buffer size;
		MPI_Offset buf_size = 0;
		for(int i=start_idx; i<end_idx; i++) {
			buf_size += _out_size_vec[i];
		}

		DataStruct* dsp = new DataStruct(buf_size);
		
		// adjust offsets based on the displacement
		MPI_Offset disp = _out_off_vec[0];
			
		for(int i=0; i<num_elements; i++) {
			elmoffset[i] = _out_off_vec[start_idx + i] - disp;
			elmsize[i] = _out_size_vec[start_idx + i];
			dsp->AddData(_out_stuff_vec[start_idx + i], _out_size_vec[start_idx + i]);
			free(_out_stuff_vec[start_idx + i]);
		}
		
		char* data_2write;

		// remember to free data_2write
		dsp->ReapData(&data_2write);

		MPI_Datatype file_type;
		MPI_Type_indexed(num_elements, elmsize, elmoffset, MPI_BYTE, &file_type);
		MPI_File_set_view(_fh, disp, MPI_BYTE, file_type, "native", MPI_INFO_NULL);
		MPI_File_write_all(_fh, data_2write, buf_size, MPI_BYTE, MPI_STATUS_IGNORE);
		MPI_Type_free(&file_type);
		
		free(data_2write);
		delete dsp;
		// free related data structure
		free(elmoffset);
		free(elmsize);
	}

	_out_off_vec.clear();
	_out_size_vec.clear();
	_out_stuff_vec.clear();

	write_time += MPI_Wtime() - write_track_time;
}

inline QueryOutput::QueryOutput(int max_size) :
    _max_size(max_size) {
    _stliner = new Streamliner(max_write_buffer);
    _qoip = new_query_output_infop();
    _is_ready = false;
    _is_offsets_sent = false;
    _outputs_freed = false;
	_stats_buf_size = 0;
	_stats_buf = NULL;

	_ready_workers = 0;
}

inline QueryOutput::~QueryOutput() {
    CleanupOutputs();

	if(_stats_buf != NULL) {
		Nlm_MemFree(_stats_buf);
	}

    delete _stliner;
}

// update QueryOutput::_qoip structure
void QueryOutput::UpdateInfo(QueryOutputInfoPtr qoip) {
	if(qoip == NULL) {
		return;
	}
	
	if(_qoip->header_size==0 && qoip->header_size!=0) {
		_qoip->header_size = qoip->header_size;
		_qoip->header = (char*)malloc(qoip->header_size);
		CHECK_NULPTR(_qoip->header);
		memcpy(_qoip->header, qoip->header, qoip->header_size);
	}
	
	if(_qoip->des_hdr_size==0 && qoip->des_hdr_size!=0) {
		_qoip->des_hdr_size = qoip->des_hdr_size;
		_qoip->des_hdr = (char*)malloc(qoip->des_hdr_size);
		CHECK_NULPTR(_qoip->des_hdr);
		memcpy(_qoip->des_hdr, qoip->des_hdr, qoip->des_hdr_size);
	}
	
	if(_qoip->des_ftr_size==0 && qoip->des_ftr_size!=0) {
		_qoip->des_ftr_size = qoip->des_ftr_size;
		_qoip->des_ftr = (char*)malloc(qoip->des_ftr_size);
		CHECK_NULPTR(_qoip->des_ftr);
		memcpy(_qoip->des_ftr, qoip->des_ftr, qoip->des_ftr_size);
	}
	
	if(_qoip->no_hits_size==0 && qoip->no_hits_size!=0) {
		_qoip->no_hits_size = qoip->no_hits_size;
		_qoip->no_hits = (char*)malloc(qoip->no_hits_size);
		CHECK_NULPTR(_qoip->no_hits);
		memcpy(_qoip->no_hits, qoip->no_hits, qoip->no_hits_size);
	}
	
	if(_qoip->stat_size==0 && qoip->stat_size!=0) {
		_qoip->stat_size = qoip->stat_size;
		_qoip->stat = (char*)malloc(qoip->stat_size);
		CHECK_NULPTR(_qoip->stat);
		memcpy(_qoip->stat, qoip->stat, qoip->stat_size);
	}
	
	if(_qoip->footer_size==0 && qoip->footer_size!=0) {
		_qoip->footer_size = qoip->footer_size;
		_qoip->footer = (char*)malloc(qoip->footer_size);
		CHECK_NULPTR(_qoip->footer);
		memcpy(_qoip->footer, qoip->footer, qoip->footer_size);
	}
	
}

// Add a list of output records to the buffered output list
void QueryOutput::AddOutputRecordsToList(ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr, bool add_orp) {
	ValNodePtr curr_node;
	list <COutputRecord>::iterator it;

	UpdateInfo(qoip);
	// add stats data here
	if(bs_ptr != NULL) {
		_stats_buf_size = BSLen(bs_ptr);
		_stats_buf = (unsigned char*)malloc(_stats_buf_size);
		CHECK_NULPTR(_stats_buf);
		BSMerge(bs_ptr, _stats_buf);
		BSFree(bs_ptr);
	}

	if(orplist == NULL || !add_orp) {
		return;
	}

	if(debug_msg) {
		LOG_MSG << "Start adding output records" <<endl;
	}
	
	// computation of result processing here can be improved as following
	// 1. check if the list if full
	// 2. if the list if full, mark the score of last record as threshold
	// 3. only add data above threshold
	// 4. this improvement could be aplied during receiving and unpacking the output records
	
	it = _sorted_output_list.begin();
	bool is_end = false;
	for(curr_node=orplist; curr_node!=NULL; curr_node=curr_node->next) {
		// !!! curr_cor serves as a temporary holder of _orp here 
		// set auto_clean to false to prevent _orp being cleaned up when curr_cor is distroyed 
		// (see COutputRecord declaration)
		
		COutputRecord curr_cor((OutputRecordPtr)(curr_node->data.ptrvalue), false); 
		if(!is_end) {
			while((it != _sorted_output_list.end()) && (curr_cor >= (*it))) {
				it++;
			}
			
			if(it == _sorted_output_list.end()) {
				is_end = true;
			}
		}
		
		// the COutputRecord added to the list will automatically clean up _orp
		if(is_end) {
			_sorted_output_list.push_back(curr_cor);  
		} else {
			_sorted_output_list.insert(it, curr_cor);
		}
	}
	
	if(debug_msg) {
		LOG_MSG << "End adding output records" <<endl;
	}
}

void QueryOutput::PruneOutputList() {
	list <COutputRecord>::iterator it;

	// simple result pruning
	if(_sorted_output_list.size() > _max_size) {
		int count = 0;
		it = _sorted_output_list.begin();
		while(++count <= _max_size) {
			it++;
		}
		_sorted_output_list.erase(it, _sorted_output_list.end());
	}
}

// calculate output offsets and add master writing stuff to FileWriter
// only called in master side
void QueryOutput::CalculateOutputOffsets(BLAST* blastp, FileWriter* fwp, MPI_Offset& curr_output_pos) {
	static int blast_align_view = blastp->GetAlignView();
	static bool blast_html = blastp->GetHtml();
	static int blast_num_desc = blastp->GetNumDes(); // number of descriptions to be showed
	static int blast_num_align = blastp->GetNumAlign(); // number of alignments to be showed
	static int max_show = blastp->GetMaxShow();
	int max_des = 0;
	int max_align = 0;
	list <COutputRecord>::iterator it;
	int count;

	MPI_Offset begin_pos = curr_output_pos;
	
	switch (blast_align_view) {
	case 0: // pairwise outupt
		max_des = blast_num_desc;
		max_align = blast_num_align;
		break;
	case 7:
		max_des = 0;
		max_align = max_show;
		break;
	case 8:
	case 9:
		max_des = 0;
		max_align = blast_num_align;
		break;
	case 10:
	case 11:
		max_des = 0;
		max_align = max_show;
		break;
	}

	set <int> workers_with_output;
	
	if(blast_align_view != 8) {
		// output header
		if(_qoip->header!=NULL) {
			_stliner->AddOutputEntry(curr_output_pos, _qoip->header_size, _qoip->header, true, true);

			_qoip->header = NULL; _qoip->header_size=0;
		}
		// output description header
		if(_qoip->des_hdr!=NULL) {
			// calculate offsets for one line descriptions
			_stliner->AddOutputEntry(curr_output_pos, _qoip->des_hdr_size, _qoip->des_hdr, true, true);

			_qoip->des_hdr = NULL; _qoip->des_hdr_size = 0;
		}
	}

	if(debug_msg) {
		LOG_MSG << "There are " << _sorted_output_list.size() << " elements in the output list" << endl;

		for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
			if((*it)._orp->rank < 0 || (*it)._orp->rank > node_count) {
				LOG_MSG << "Error: invalid worker rank " << (*it)._orp->rank << endl;
				throw __FILE__ "Ivalid worker rank!";
			}
		}
	}
	
	if(_sorted_output_list.size() > 0) {

		if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
				workers_with_output.insert(((*it)._orp)->rank);
			}

			set <int>::iterator sit;
			for(sit=workers_with_output.begin(); sit!=workers_with_output.end(); sit++) {
				_stliner->InitAccessLocation(*sit);
			}
		}
			
		if(blast_align_view < 8) {
			count = 0;
			for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
				if(++count > max_des) break;
				
				(*it).UpdateDescriptionOffset(curr_output_pos, true);

				if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
					_stliner->AddAccessLocation(((*it)._orp)->rank, ((*it)._orp)->des_offset, ((*it)._orp)->des_size);
				}
				
			}
		}

		if(blast_align_view != 7) {
			if(_qoip->des_ftr!=NULL) {
				// output description footer
				_stliner->AddOutputEntry(curr_output_pos, _qoip->des_ftr_size, _qoip->des_ftr, true, true);

				_qoip->des_ftr = NULL; _qoip->des_ftr_size = 0;
			}
		}

		count = 0;		
		// calculate offsets for alignments
		for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {

			if(++count > max_align) break;

			if(blast_align_view == 7) {
				// print hit number
				char* tmp_buf = (char*)malloc(100); // this buffer will be freed in FileWriter
				CHECK_NULPTR(tmp_buf);
				char* ptr = tmp_buf;
				sprintf(ptr, "        <Hit>\n");
				ptr += strlen(tmp_buf);
				sprintf(ptr, "          <Hit_num>%d</Hit_num>\n", count);
				_stliner->AddOutputEntry(curr_output_pos, strlen(tmp_buf), tmp_buf, true, true);      
			}
			
			if(blast_align_view == 10) {
				char* tmp_buf = (char*)malloc(10); // this buffer will be freed in FileWriter
				CHECK_NULPTR(tmp_buf);
				if(it == _sorted_output_list.begin()) {
					sprintf(tmp_buf, "{"); // open struct for the first hit
				} else {
					sprintf(tmp_buf, " ,\n      {"); // open struct for the rest of hits
				}
				_stliner->AddOutputEntry(curr_output_pos, strlen(tmp_buf), tmp_buf, true, true);      
			}
			
			(*it).UpdateAlignOffset(curr_output_pos, true);
			if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
				_stliner->AddAccessLocation(((*it)._orp)->rank, ((*it)._orp)->aln_offset, ((*it)._orp)->aln_size);
			}
		}
		
		if(blast_align_view == 7) {
			if(_qoip->des_ftr!=NULL) {
				// output description footer
				_stliner->AddOutputEntry(curr_output_pos, _qoip->des_ftr_size, _qoip->des_ftr, true, true);

				_qoip->des_ftr = NULL; _qoip->des_ftr_size = 0;
			}
		}
		
	} else {
		if(blast_align_view != 7) {
			if(_qoip->no_hits!=NULL) {
				_stliner->AddOutputEntry(curr_output_pos, _qoip->no_hits_size, _qoip->no_hits, true, true);

				_qoip->no_hits = NULL; _qoip->no_hits_size = 0;
			}
		}
	}

	// output statistics
	if(_qoip->stat!=NULL) {
		_stliner->AddOutputEntry(curr_output_pos, _qoip->stat_size, _qoip->stat, true, true);

		_qoip->stat = NULL; _qoip->stat_size = 0;
	}

	if(_sorted_output_list.size() == 0 && blast_align_view == 7) {
		if(_qoip->no_hits!=NULL) {
			_stliner->AddOutputEntry(curr_output_pos, _qoip->no_hits_size, _qoip->no_hits, true, true);
			
			_qoip->no_hits = NULL; _qoip->no_hits_size = 0;
		}
	}
	
	// output footer
	if(_qoip->footer!=NULL) {
		_stliner->AddOutputEntry(curr_output_pos, _qoip->footer_size, _qoip->footer, true, true);

		_qoip->footer = NULL;
	}

	if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		
		if(debug_msg) {
			LOG_MSG << "There are " << workers_with_output.size() << " workers with output" << endl;
		}

		_stliner->AddWorkerWOutputs(workers_with_output);
	}

	_stliner->SetOutputLength(curr_output_pos - begin_pos);
}

// adding data to FileWriter
void QueryOutput::PrepareOutputData(BLAST *blastp, FileWriter* fwp) {
	static int blast_align_view = blastp->GetAlignView();
	static bool blast_html = blastp->GetHtml();
	static int blast_num_desc = blastp->GetNumDes(); // number of descriptions to be showed
	static int blast_num_align = blastp->GetNumAlign(); // number of alignments to be showed
	list <COutputRecord>::iterator it;

	if(_sorted_output_list.size() == 0) { // no results for this query
		return;
	}

	MPI_Offset query_begin = _stliner->GetBeginPos();
		
	switch (blast_align_view) {
	case 0: // pairwis output
		// _sorted_output_list should have been pruned
		
		// output descriptions

		for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
			OutputRecordPtr curr_orp = (*it)._orp;
			if(curr_orp->des_offset != -1) {
                if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
                    _stliner->AddOutputEntry(curr_orp->des_offset, curr_orp->des_size, curr_orp->des_data, false);
                } else {
					MPI_Offset real_off = query_begin + curr_orp->des_offset;
                    fwp->AddWriteEntry(real_off, curr_orp->des_size, curr_orp->des_data, false);
                }

				curr_orp->des_data = NULL; // set to NULL so it won't be cleanup by ~COutputRecord()
			} else {
				break;
			}
		}
		
	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
		// output alignments
		for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
			OutputRecordPtr curr_orp = (*it)._orp;
			if(curr_orp->aln_offset != -1) {
                if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
                    _stliner->AddOutputEntry(curr_orp->aln_offset, curr_orp->aln_size, curr_orp->aln_data, false);
                } else {
					MPI_Offset real_off = query_begin + curr_orp->aln_offset;
                    fwp->AddWriteEntry(real_off, curr_orp->aln_size, curr_orp->aln_data, false);
                }

				curr_orp->aln_data = NULL; // set to NULL so it won't be cleanup by ~COutputRecord()
			} else {
				break;
			}
		}
	}
}

void QueryOutput::AdjustOutputOffsets(MPI_Offset& curr_output_pos) {
	
	_stliner->SetBeginPos(curr_output_pos);	
	curr_output_pos += _stliner->GetOutputLength();
}

// Debug function, print out buffered output records
void QueryOutput::DebugPrint(FILE* fp) {
	fprintf(fp, "Query Output Info\n");
	fflush(fp);
	fprintf(fp, "_qoip->header=");
	fwrite(_qoip->header, _qoip->header_size, 1, fp);
	fprintf(fp, "_qoip->des_hdr=");
	fwrite(_qoip->des_hdr, _qoip->des_hdr_size, 1, fp);
	fprintf(fp, "_qoip->des_ftr=");
	fwrite(_qoip->des_ftr, _qoip->des_ftr_size, 1, fp);
	fprintf(fp, "_qoip->no_hits=");
	fwrite(_qoip->no_hits, _qoip->no_hits_size, 1, fp);
	fprintf(fp, "_qoip->stat=");
	fwrite(_qoip->stat, _qoip->stat_size, 1, fp);
	fprintf(fp, "_qoip->footer=");
	fwrite(_qoip->footer, _qoip->footer_size, 1, fp);
	fprintf(fp, "\nNumber of Output Records: %d\n", _sorted_output_list.size());
	fflush(fp);

	list < COutputRecord >::iterator it;
	for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
		(*it).DebugPrint(fp);
	}

	fprintf(fp, "\n");
	fflush(fp);
}

void QueryOutput::DebugWriteAsn(FILE* fp) {
	// header
	fwrite(_qoip->header, _qoip->header_size, 1, fp);
	
	// alignments
	list < COutputRecord >::iterator it;
	for(it=_sorted_output_list.begin(); it!=_sorted_output_list.end(); it++) {
		fwrite((*it)._orp->aln_data, (*it)._orp->aln_size, 1, fp);
	}
	
	// footer
	fwrite(_qoip->footer, _qoip->footer_size, 1, fp);
}

void QueryOutput::CleanupOutputs() {
    if(!_outputs_freed) {
        _outputs_freed = true;
    } else {
        return;
    }

    _sorted_output_list.clear();
    free_query_output_infop(_qoip);
}

// check whether non-blocking sends have been finished
void IsendMgr::CheckPendingSends() {
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
void IsendMgr::WaitPendingSends() {
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

int IrecvMgr::WaitAnyRecv(CommRecvStruct* &crsp, MPI_Status& status) {
	// wait a finished recv, remove it from the list
	int num_reqs = _recv_ops_list.size();

	list <CommRecvStruct*>::iterator lit;
	
	MPI_Request* req_array = new MPI_Request[num_reqs];
	int i=0;

	for(lit=_recv_ops_list.begin(); lit!=_recv_ops_list.end(); lit++) {
		req_array[i++] = (*lit)->GetReq();
	}

	int index = -1;
	MPI_Waitany(num_reqs, req_array, &index, &status);
	
	assert(index >= 0);

	lit=_recv_ops_list.begin(); 
	for(i=0; i<index; i++) {
		lit++;
	}

	crsp = *lit;
	_recv_ops_list.erase(lit);
	
	delete req_array;	

	return 0;
}

int IrecvMgr::TestAnyRecv(CommRecvStruct* &crsp, MPI_Status& status) {
	// wait a finished recv, remove it from the list
	int num_reqs = _recv_ops_list.size();

	list <CommRecvStruct*>::iterator lit;
	
	MPI_Request* req_array = new MPI_Request[num_reqs];
	int i=0;

	for(lit=_recv_ops_list.begin(); lit!=_recv_ops_list.end(); lit++) {
		req_array[i++] = (*lit)->GetReq();
	}

	int index = -1;
	int flag = 0;
	MPI_Testany(num_reqs, req_array, &index, &flag, &status);
	
	delete req_array;	

	if(flag) {
		lit=_recv_ops_list.begin(); 
		for(i=0; i<index; i++) {
			lit++;
		}

		crsp = *lit;
		_recv_ops_list.erase(lit);

		return 1;
	}

	return 0;
}

IrecvMgr::~IrecvMgr() {
	list <CommRecvStruct*>::iterator lit;
	for(lit=_recv_ops_list.begin(); lit!=_recv_ops_list.end(); lit++) {
		(*lit)->Cancel();
		delete(*lit);
	}
	_recv_ops_list.clear();
}

// debug function, print out a list of output records
void Writer::PrintOutputsList(FILE *fp) {
	list <OutputRecordPtr>::iterator it;
	
	for(int i=0; i<_query_count; i++) {
		fflush(fp);
		fprintf(fp, "\nQUERY = %d\n", i);
		fflush(fp);
		_output_vec[i]->DebugPrint(fp);
	}
}

void Writer::WaitAllWriteOps() {
	list < pair<int, int> >::iterator lit;

	while(!_writing_ops.empty()) {
		lit = _writing_ops.begin();
		int query_id = (*lit).first;
		int write_segment_id = (*lit).second;

		_output_vec[query_id]->_stliner->InitAWrite(query_id, write_segment_id);

		if(debug_msg) {
			LOG_MSG << "Before wait for writing query " << query_id << " segment " << write_segment_id << endl;
		}

		if(output_strategy == MASTER_STREAMLINE && group_rank != writer_process) {
			if(!_output_vec[query_id]->_stliner->IsEmpty()) {
				// notify WriterMaster
				int event_array[2];
				event_array[0] = WAIT_WRITE_COMPLETE;
				event_array[1] = 2*int_size;
				CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
				cssp_event->AddData(&event_array[0], 2*int_size);
				cssp_event->IsendData(group_comm, writer_process, EVENT_TYPE);
				_comm_send.AddPendingSend(cssp_event);

				CommSendStruct* cssp_data = new CommSendStruct(2*int_size);
				cssp_data->AddData(&query_id, int_size);
				cssp_data->AddData(&write_segment_id, int_size);
				cssp_data->IsendData(group_comm, writer_process, WAIT_WRITE_COMPLETE);
				_comm_send.AddPendingSend(cssp_data);
			}
		}
	
		_output_vec[query_id]->_stliner->WaitAWrite(write_segment_id, _fw);

		if(debug_msg) {
			LOG_MSG << "After wait for writing query " << query_id << " segment " << write_segment_id << endl;
		}

		// check if the query has been finished
		if(_output_vec[query_id]->_stliner->IsWritten()) {

			if(output_strategy == WORKER_STREAMLINE && group_rank == _output_vec[query_id]->_stliner->GetWriterLeader()) {
				SendWriteAck(query_id, writer_process);
			}

			if( output_strategy == WORKER_STREAMLINE ||
					(output_strategy == MASTER_STREAMLINE && write_segment_id == _output_vec[query_id]->_stliner->GetNumWriteSegments() - 1)) {
				// delete _output_vec[query_id];
				// _output_vec[query_id] = NULL;
				UnloadQuery(query_id);
				if(_last_written_query < query_id + 1) {
					_last_written_query = query_id + 1;
				}
				if(_last_output_query < query_id + 1) {
					_last_output_query = query_id + 1;
				}
			}
		}
	
		_writing_ops.erase(lit);
	}
}

void Writer::ProcessAsyncOutput() {
	
	if(_writing_ops.empty()) {
		return;
	}

	list < pair<int, int> >::iterator lit, tmp_lit;

	lit = _writing_ops.begin();
	while(lit != _writing_ops.end()) {
		int query_id = (*lit).first;
		int write_segment_id = (*lit).second;

		bool is_finished = false;
		
		_output_vec[query_id]->_stliner->InitAWrite(query_id, write_segment_id);
		int ret = _output_vec[query_id]->_stliner->CheckAWrite(write_segment_id, _fw);

		tmp_lit = lit;
		lit++;

		if(ret == 1) {
			is_finished = true;

			if(output_strategy == WORKER_STREAMLINE) {
				_writing_ops.erase(tmp_lit);

				if(_output_vec[query_id]->_stliner->IsWritten()) {

					if(group_rank == _output_vec[query_id]->_stliner->GetWriterLeader()) {
						SendWriteAck(query_id, writer_process);
					}

					// delete _output_vec[query_id];
					// _output_vec[query_id] = NULL;
					UnloadQuery(query_id);
					
					_last_written_query = query_id + 1;
					_last_output_query = query_id + 1;
				}
			} 
		}

		if( output_strategy == MASTER_STREAMLINE ) {
			if (group_rank == writer_process && 
				_output_vec[query_id]->_stliner->PeersReady(write_segment_id)) {
				if(debug_msg) {
					LOG_MSG << "Before waiting write for query " << query_id << " segment " << write_segment_id << endl;
				}
				
				_output_vec[query_id]->_stliner->WaitAWrite(write_segment_id, _fw);

				if(debug_msg) {
					LOG_MSG << "After waiting write for query " << query_id << " segment " << write_segment_id << endl;
				}
				
				_writing_ops.erase(tmp_lit);

				if(_output_vec[query_id]->_stliner->IsWritten() && write_segment_id == _output_vec[query_id]->_stliner->GetNumWriteSegments() - 1) {

					// delete _output_vec[query_id];
					// _output_vec[query_id] = NULL;
					UnloadQuery(query_id);

					if(_last_written_query < query_id + 1) {
						_last_written_query = query_id + 1;
					}
					if(_last_output_query < query_id + 1) {
						_last_output_query = query_id + 1;
					}
				}

				is_finished = true;
			}
		}
		
		if(!is_finished) {	
			// avoid dead loop
			break;
		}
	}
}

void Writer::SendWriteAck(int query_id, int dst) {

	if(debug_msg) {
		LOG_MSG << "Finish writing of query " << query_id << endl;
	}

	int event_array[2];
	event_array[0] = QUERY_WRITE_ACK;
	event_array[1] = int_size;
	CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
	cssp_event->AddData(&event_array[0], 2 * int_size);
	cssp_event->IsendData(group_comm, writer_process, EVENT_TYPE);
	_comm_send.AddPendingSend(cssp_event);

	CommSendStruct* cssp_data = new CommSendStruct(int_size);
	cssp_data->AddData(&query_id, int_size);
	cssp_data->IsendData(group_comm, writer_process, QUERY_WRITE_ACK);
	_comm_send.AddPendingSend(cssp_data);
}

void Writer::AddWorkingQueries(int segment_id, int start_query, int end_query) {
	for(int i=start_query; i<end_query; i++) {
		_assigned_queries++;
		_working_queries.push_back(i);
		InitQueryOutput(segment_id, i);
	}
}

void Writer::InitQueryOutput(int segment_id, int query_id) {
	if(_output_vec.find(query_id) == _output_vec.end()) {
		if(debug_msg) {
			LOG_MSG << "Init output structure for query " << query_id << endl;
		}
		_query_map[query_id] = segment_id;
		_output_vec[query_id] = new QueryOutput(_blast->GetMaxShow());
	}
}

void Writer::UnloadQuery(int query_id) {
	if(_output_vec.find(query_id) == _output_vec.end()) {
		throw __FILE__ " Writer::UnloadQuery -- the query has not been initialized";
	}

	if(debug_msg) {
		LOG_MSG << "Unload query " << query_id << endl;
	}
	
	_query_map.erase(query_id);
	delete _output_vec[query_id];
	_output_vec.erase(query_id);
}

void Writer::UpdateDefaultReportInfo(QueryOutputInfoPtr qoip) {
	assert(qoip != NULL);
	assert(_output_drp == NULL); // this function can only be called once

	_output_drp = new_default_reportp();

	char* find_ptr;
	int align_view = _blast->GetAlignView();

	// !!!! check what will happen for all output format
	
	int copy_len = 0;
	if(align_view == 0) { 
		if(qoip->header_size > 0) {
			// reference
			find_ptr = strstr(qoip->header, "Query=");
			copy_len = find_ptr - qoip->header;
			assert(copy_len > 0);
			_output_drp->ref_size = copy_len;
			_output_drp->ref = (char*)malloc(copy_len);
			CHECK_NULPTR(_output_drp->ref);
			memcpy(_output_drp->ref, qoip->header, copy_len);
		
			// db info
			find_ptr = strstr(qoip->header, "Database:");
			int start_pos = find_ptr - qoip->header;
			copy_len = qoip->header_size - start_pos;
			_output_drp->db_info_size = copy_len;
			_output_drp->db_info = (char*)malloc(copy_len);
			CHECK_NULPTR(_output_drp->db_info);
			memcpy(_output_drp->db_info, find_ptr, copy_len);
		}

		// no hits
		if(qoip->no_hits_size > 0) {
			_output_drp->no_hits_size = qoip->no_hits_size;
			_output_drp->no_hits = (char*)malloc(qoip->no_hits_size);
			CHECK_NULPTR(_output_drp->no_hits);
			memcpy(_output_drp->no_hits, qoip->no_hits, qoip->no_hits_size);
		}

		// db report
		if(qoip->stat_size > 0) {
			assert(qoip->stat != NULL);
			find_ptr = strstr(qoip->stat, "Lambda");
			copy_len = find_ptr - qoip->stat;
			assert(copy_len > 0);
			_output_drp->db_report_size = copy_len;
			_output_drp->db_report = (char*)malloc(copy_len);
			CHECK_NULPTR(_output_drp->db_report);
			memcpy(_output_drp->db_report, qoip->stat, copy_len);
		}

		// footer
		if(qoip->footer_size > 0) {
			_output_drp->footer_size = qoip->footer_size;
			_output_drp->footer = (char*)malloc(qoip->footer_size);
			CHECK_NULPTR(_output_drp->footer);
			memcpy(_output_drp->footer, qoip->footer, qoip->footer_size);
		}
	}
}

// sigleton pattern, initialize _instance to NULL
WriterMasterIndividual* WriterMasterIndividual::_instance = NULL;

WriterMasterStreamline* WriterMasterStreamline::_instance = NULL;

// singleton pattern, get instance for this class
WriterMasterStreamline* WriterMasterStreamline::Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int send_batch) {
	if(_instance == NULL) {
		_instance = new WriterMasterStreamline(blast, query_count, total_frags, output_file, send_batch);
	} 
	
	return _instance;
}

WriterMasterStreamline* WriterMasterStreamline::Instance() {
	return _instance;
}

WriterMasterCollective* WriterMasterCollective::_instance = NULL;

// receive partial output records from workers
/* The format of packed message:
	- packing type
	- fragment id
	- number of queries results
	- query_id, number of records
	- description header length
	- description header data
	- description footer length
	- description footer data
	- actual records 
	- query_id, number of records
	- actual records 
	- ......
 */ 
void WriterMaster::RecvPartialOutputs(CommRecvStruct* crsp, int src) {
	ValNodePtr orplist, curr_node, last_node;
	int query_id;
	bool last_query_ready = false;
	static int blast_align_view = _blast->GetAlignView();
	static bool blast_html = _blast->GetHtml();
	
	// unpack the partial outputs to a ValNode list
	int pack_type;
	crsp->ExtractData(&pack_type, int_size);
	int frag_id;
	crsp->ExtractData(&frag_id, int_size);
	int num_queries;
	crsp->ExtractData(&num_queries, int_size);
	
	for(int i=0; i<num_queries; i++) {
		crsp->ExtractData(&query_id, int_size);
		if(debug_msg) {
			LOG_MSG << "Start receiving partial outputs: query_id " << query_id << " frag_id " << frag_id << " from " << " src " << src << endl;
		}
		int num_orps;
		crsp->ExtractData(&num_orps, int_size);
		
		QueryOutputInfoPtr qoip = new_query_output_infop();
		
		int tmp_size;
		if(blast_align_view < 10) { // not asn
			// receive description header, description footer
			crsp->ExtractData(&tmp_size, int_size);
			if(tmp_size > 0) {
				qoip->des_hdr_size = tmp_size;
				qoip->des_hdr = (char*)malloc(tmp_size);
				CHECK_NULPTR(qoip->des_hdr);
				crsp->ExtractData(qoip->des_hdr, tmp_size);
			}
			crsp->ExtractData(&tmp_size, int_size);
			if(tmp_size > 0) {
				qoip->des_ftr_size = tmp_size;
				qoip->des_ftr = (char*)malloc(tmp_size);
				CHECK_NULPTR(qoip->des_ftr);
				crsp->ExtractData(qoip->des_ftr, tmp_size);
			}
		} else { // asn
			// receive header and footer
			crsp->ExtractData(&tmp_size, int_size);
			if(tmp_size > 0) {
				qoip->header_size = tmp_size;
				qoip->header = (char*)malloc(tmp_size);
				CHECK_NULPTR(qoip->header);
				crsp->ExtractData(qoip->header, tmp_size);
			}
			crsp->ExtractData(&tmp_size, int_size);
			if(tmp_size > 0) {
				qoip->footer_size = tmp_size;
				qoip->footer = (char*)malloc(tmp_size);
				CHECK_NULPTR(qoip->footer);
				crsp->ExtractData(qoip->footer, tmp_size);
			}
		}

		// extract stats buffer
		crsp->ExtractData(&tmp_size, int_size);
		if(tmp_size > 0) {
			unsigned char* stats_buf = (unsigned char*)malloc(tmp_size);
			CHECK_NULPTR(stats_buf);
			crsp->ExtractData(stats_buf, tmp_size);

			_blast->ExtractStats(query_id, stats_buf, tmp_size);	

			free(stats_buf);
		}	
		
		orplist = NULL;
		last_node = NULL;
		curr_node = NULL;
		if(num_orps > 0) {
			for(int j=0; j<num_orps; j++) {
				OutputRecordPtr curr_orp;
				curr_orp = (OutputRecordPtr)malloc(sizeof(OutputRecord));
				CHECK_NULPTR(curr_orp);
	
				curr_orp->rank = -1;
				curr_orp->frag_id = -1;
				curr_orp->id = NULL;
				curr_orp->des_size = 0;
				curr_orp->des_data = NULL;
				curr_orp->aln_size = 0;
				curr_orp->aln_data = NULL;
				curr_orp->des_offset = -1;
				curr_orp->aln_offset = -1;
				
				unpack_record(pack_type, curr_orp, crsp);
				curr_node=ValNodeAddPointer(&last_node, 0, curr_orp);
				if(j==0) { // first record of the query
					orplist = last_node;
				}
				last_node = curr_node;
			}
		}
		
		// insert to the list
		AddOutputRecords(query_id, orplist, qoip, NULL);

		_output_vec[query_id]->PruneOutputList();
		_output_vec[query_id]->_stliner->AddSearchedWorker(src);
		
		free_query_output_infop(qoip);
		// free up ValNode list
		orplist = ValNodeFree(orplist);
		
		if(query_data.find(query_id) == query_data.end()) {
			query_data[query_id] = new QueryData(NULL, NULL);
		}
		query_data[query_id]->searched_frags++;
		if(query_data[query_id]->searched_frags == _total_frags) { //received all outputs for this query
			delete query_data[query_id];
			query_data.erase(query_id);
			
			_output_vec[query_id]->_is_ready = true; // mark this query ready for writing

			double track_time = MPI_Wtime();
			if(use_brief_report && _blast->GetAlignView() != 7) {
				_blast->SynthesizeReport(query_id, _output_drp, _output_vec[query_id]->_qoip);
			} else {
				_blast->GetQueryOutputInfo(query_id);
			}
			output_info_time += MPI_Wtime() - track_time;

			// unload query data
			QueryM::Instance()->RemoveQuery(query_id);

			MPI_Offset zero_pos = 0;
			
			if(debug_msg) {
				LOG_MSG << "Start calculating output offsets for query "<< query_id <<endl;
			}
			
			_output_vec[query_id]->CalculateOutputOffsets(_blast, _fw, zero_pos);
			
			if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
				_output_vec[query_id]->_stliner->CalculateMessageSizes();
				//_output_vec[query_id]->_stliner->DumpMessageSizes(stderr);
			}

			if(debug_msg) {
				LOG_MSG << "End calculating output offsets for query "<< query_id <<endl;
			}
			
			_offsets_to_send++;
		}
	}
			
	if(debug_msg) {
		LOG_MSG << "End receiving partial outputs" << endl;
	}

	if(_offsets_to_send > 0) {
		CheckUnsentOffsets();
	}
}

void WriterMaster::CheckUnsentOffsets() {

	if(_offsets_to_send > 0) {
		if(debug_msg) {
			LOG_MSG << "Start checking unsent offsets" <<endl;
		}
	
		list <int>::iterator lit;
		// for(int wrk_idx=_last_sent_offset; wrk_idx<_working_queries.size(); wrk_idx++) {
		for(lit=_working_queries.begin(); lit!=_working_queries.end(); lit++) {

			// int curr_query_id = _working_queries[wrk_idx];
			int curr_query_id = *lit;
			
			if(_output_vec.find(curr_query_id) == _output_vec.end())
				break;
	
			if(!_output_vec[curr_query_id]->_is_ready) 
				break;

			if(_output_vec[curr_query_id]->_is_offsets_sent) {
				throw __FILE__ "Offsets have been sent for this query";
			}
		
			// currently only process one query at a time		
            int start_query = curr_query_id;
            int end_query = curr_query_id + 1;

			_output_vec[start_query]->AdjustOutputOffsets(_curr_output_pos);
			if(output_strategy == MASTER_STREAMLINE) {

				//_output_vec[start_query]->_stliner->SetWriterLeader(writer_process);
				_output_vec[start_query]->_stliner->SetWriterLeader(group_rank);
				_output_vec[start_query]->_stliner->SetAssignLeader();
				_output_vec[start_query]->_stliner->PrepareCommWriteData();

				int num_write_segments = _output_vec[start_query]->_stliner->GetNumWriteSegments();
				AddWriteOps(start_query, num_write_segments);
				_lead_queries.push(start_query);
				_output_vec[start_query]->_stliner->AlignOutputEntries();
				
				// Write will be initialized in ProcessAsyncOutput
				if(_writing_ops.empty() && num_write_segments > 0 ) {
					_output_vec[start_query]->_stliner->InitAWrite(start_query, 0);
				}
			}

			SendOutputOffsets(start_query, end_query);
			_output_vec[start_query]->_is_offsets_sent = true;

			QueryM::Instance()->UnloadOldQueries(start_query);

            if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_COLLECTIVE) {
                _output_vec[curr_query_id]->CleanupOutputs();
//                _query_to_write = start_query;
            } else {
				if(output_strategy == WORKER_SPLIT || output_strategy == WORKER_INDIVIDUAL) {
					// master print
					_output_vec[curr_query_id]->_stliner->WriteOutputEntries(_fw, _curr_output_pos);

					if(dump_raw_output) {
						_fw->DumpWriteRequests(dbgfp);
					}
	
					_fw->Write();
				} 

                // clean up query output data
                // delete _output_vec[curr_query_id];
				// _output_vec[curr_query_id] = NULL;
				UnloadQuery(curr_query_id);
            }

			_last_sent_offset++;
			_offsets_to_send--;
		}

		_working_queries.erase(_working_queries.begin(), lit);

		if(debug_msg) {
			LOG_MSG << "After checking unsent offsets" <<endl;
		}
	}
}

void WriterMaster::TellWriterWorkerQuit() {
	int command = WRITER_QUIT;

	for(int i=0; i<group_node_count; i++) {
		if(!IsWorker(i)) {
			continue;
		}

		// send event
		int assign_array[3];
		assign_array[0] = WRITER_EVENT;
		assign_array[1] = WRITER_QUIT;
		assign_array[2] = -1;
		CommSendStruct* cssp_event = new CommSendStruct(3*int_size);
		cssp_event->AddData(&assign_array[0], 3 * int_size);
		cssp_event->SendData(group_comm, i, ASSIGNMENT_TYPE);
		delete cssp_event;

		MPI_Send(&command, 1, MPI_INT, i, WRITER_QUIT, group_comm);
	}
}

// handle writer related messages that received from worker side
int WriterMasterStreamline::HandleMessages(MPI_Status &event_status, int* event_array) {
    _comm_send.CheckPendingSends();

	MPI_Status probe_status;
	MPI_Probe( event_status.MPI_SOURCE, event_array[0], group_comm, &probe_status );

	int msg_size;
	MPI_Get_count(&probe_status, MPI_BYTE, &msg_size);
	CommRecvStruct* crsp = new CommRecvStruct(msg_size);
	crsp->RecvData(group_comm, probe_status.MPI_SOURCE, probe_status.MPI_TAG);

	switch (probe_status.MPI_TAG) {
		case PARTIAL_OUTPUTS_DATA:
			//receive partial outputs
			RecvPartialOutputs(crsp, event_status.MPI_SOURCE);
			break;
		case WRITE_REQUEST:
			ProcessWriteRequest(probe_status.MPI_SOURCE);
			break;
		case WRITE_COMPLETE:
			ProcessWriteComplete(probe_status.MPI_SOURCE);
			break;
		case WRITER_FINALIZED:
			_complete_worker++;
			break;
		case QUERY_WRITE_ACK:
			FinishQuery(crsp);
			break;
		case WAIT_WRITE_COMPLETE:
			ProcessWaitWriteNotification(crsp);
			break;
		default:
			throw __FILE__ "Don't know how to process this tag";
	}

	if(crsp!=NULL) { // allow crsp to be freed earlier elsewhere to save memory
		delete crsp;
	}

	if(!_writing_ops.empty()) {
		ProcessAsyncOutput();
	}
	
	if(WriterWorkerComplete()) {
		return 1;
	}

    return 0;
}

int WriterMasterStreamline::HandleMessages(CommRecvStruct* crsp, MPI_Status &recv_status) {
    _comm_send.CheckPendingSends();

	switch (recv_status.MPI_TAG) {
		case PARTIAL_OUTPUTS_DATA:
			//receive partial outputs
			RecvPartialOutputs(crsp, recv_status.MPI_SOURCE);
			break;
		case WRITE_REQUEST:
			ProcessWriteRequest(recv_status.MPI_SOURCE);
			break;
		case WRITE_COMPLETE:
			ProcessWriteComplete(recv_status.MPI_SOURCE);
			break;
		case WRITER_FINALIZED:
			_complete_worker++;
			break;
		case QUERY_WRITE_ACK:
			FinishQuery(crsp);
			break;
		case WAIT_WRITE_COMPLETE:
			ProcessWaitWriteNotification(crsp);
			break;
		default:
			throw __FILE__ "Don't know how to process this tag";
	}

	if(crsp!=NULL) { // allow crsp to be freed earlier elsewhere to save memory
		delete crsp;
	}

	if(!_writing_ops.empty()) {
		ProcessAsyncOutput();
	}
	
	if(WriterWorkerComplete()) {
		return 1;
	}

    return 0;
}

void WriterMasterStreamline::RecvMetaDataToWrite(CommRecvStruct* crsp, int src, FileWriter* fwp) {
    int query_id;

    crsp->ExtractData(&query_id, int_size);

    int ret_code = _output_vec[query_id]->_stliner->AddPendingRecvs(src);

	int num_data_msgs;
    crsp->ExtractData(&num_data_msgs, int_size);

	// receive msg sizes
	vector <int> msg_sizes;

	for(int i=0; i<num_data_msgs; i++) {
		int msg_size; 
		crsp->ExtractData(&msg_size, int_size);

		msg_sizes.push_back(msg_size);
	}
	_output_vec[query_id]->_stliner->AddMessageSizes(src, msg_sizes);
	
	if(debug_msg) {
		LOG_MSG << "Received MetaData of query " << query_id << " from " << src << endl;
	}
}

void WriterMasterStreamline::FinishQuery(CommRecvStruct* crsp) {
	int query_id;
	crsp->ExtractData(&query_id, int_size);

	_num_written_queries++;

	// if( _no_more_queries && _num_written_queries == _working_queries.size() ) {
	if( _no_more_queries && _num_written_queries == _assigned_queries ) {
		TellWriterWorkerQuit();
	}
}

void WriterMasterStreamline::ProcessWaitWriteNotification(CommRecvStruct* crsp) {
	int wait_query_id = -1;
	int wait_write_segment_id = -1;

	crsp->ExtractData(&wait_query_id, int_size);
	crsp->ExtractData(&wait_write_segment_id, int_size);

	if(debug_msg) {
		LOG_MSG << "Received waiting for write notification of query " << wait_query_id << " segment " << wait_write_segment_id << endl;
	}
	
	list < pair<int, int> >::iterator lit, tmp_lit;

	lit = _writing_ops.begin();
	while(lit != _writing_ops.end()) {
		int query_id = lit->first;
		int write_segment_id = lit->second;

		if(query_id == wait_query_id && write_segment_id == wait_write_segment_id) {
			break;
		}
		lit++;
	}

	if(lit == _writing_ops.end()) {
		throw __FILE__ " WriterMasterStreamline::ProcessWaitWriteNotification() - cannot find the write Op";
	}

	_output_vec[wait_query_id]->_stliner->HandleWaitWriteNotify(wait_write_segment_id);
}

int WriterWorkerStreamline::HandleMessages(MPI_Status &event_status, int* event_array) {
    
	_comm_send.CheckPendingSends();
	_results_send.CheckPendingSends();

	MPI_Status probe_status;
	MPI_Probe( event_status.MPI_SOURCE, event_array[0], group_comm, &probe_status );

	int msg_size;
	MPI_Get_count(&probe_status, MPI_BYTE, &msg_size);
	CommRecvStruct* crsp = new CommRecvStruct(msg_size);
	crsp->RecvData(group_comm, probe_status.MPI_SOURCE, probe_status.MPI_TAG);
	
	switch (probe_status.MPI_TAG) {
		case OUTPUT_OFFSETS:
			// receive output offsets for local buffered records
			ReceiveOutputOffsets(crsp);
			break;
		case SHIFT_WRITE_LEADER: // only for WORKER_STREAMLINE
			ReceiveLeaderData(crsp);
			break;
		case META_DATA_TO_WRITE: // only for WORKER_STREAMLINE 
			RecvPeerMetaData(crsp, probe_status.MPI_SOURCE);
			break;
		case WRITER_QUIT: // only for WORKER_STREAMLINE
			_is_finished = true;
			break;
		default:
			throw __FILE__ " WriterWorkerStreamline::HandleMessages -- don't know how to process this tag";
	}

	delete crsp;

	if(_writing_ops.size() > 0) {
		ProcessAsyncOutput();
	}
	
	return 0;
}

void WriterWorkerStreamline::ReceiveLeaderData(CommRecvStruct* crsp) {
	if(debug_msg) {
		LOG_MSG << "Start receiving lead data" <<endl;
	}

	int query_id = -1;
	crsp->ExtractData(&query_id, int_size);
	
	MPI_Offset begin_pos = -1;
	crsp->ExtractData(&begin_pos, offset_size);

	_output_vec[query_id]->_stliner->SetWriterLeader(group_rank);
	_output_vec[query_id]->_stliner->SetBeginPos(begin_pos);

	MPI_Offset output_len = -1;
	crsp->ExtractData(&output_len, offset_size);
	_output_vec[query_id]->_stliner->SetOutputLength(output_len);

	_output_vec[query_id]->_stliner->ExtractOutputEntries(crsp);
	_output_vec[query_id]->_stliner->ExtractSubmitters(crsp);
	_output_vec[query_id]->_stliner->ExtractWritePeers(crsp);

	_lead_queries.push(query_id);

	_output_vec[query_id]->_stliner->SetAssignLeader();
	// _output_vec[query_id]->_stliner->LeaderCheckReady();
	
	// _output_vec[query_id]->_stliner->PrepareCommWriteData();

	if(debug_msg) {
		LOG_MSG << "End receiving lead data of " << query_id <<endl;
	}
}

void WriterWorkerStreamline::RecvPeerMetaData(CommRecvStruct* crsp, int src) {
    int query_id;

    crsp->ExtractData(&query_id, int_size);

    int ret_code = _output_vec[query_id]->_stliner->AddPendingRecvs(src);

	int num_data_msgs;
    crsp->ExtractData(&num_data_msgs, int_size);

	// receive msg sizes
	vector <int> msg_sizes;

	for(int i=0; i<num_data_msgs; i++) {
		int msg_size; 
		crsp->ExtractData(&msg_size, int_size);

		msg_sizes.push_back(msg_size);
	}
	_output_vec[query_id]->_stliner->AddMessageSizes(src, msg_sizes);
	
	if(debug_msg) {
		LOG_MSG << "Received MetaData of query " << query_id << " from " << src << endl;
	}
}

void WriterWorkerStreamline::OutputReadyQueries() {
	while(!_lead_queries.empty()) {
		int query_id = _lead_queries.front();

		if(_output_vec[query_id]->_stliner->IsReady()) {
			if(debug_msg) {
				LOG_MSG << "Start gathering outputs for " << query_id <<endl;
			}

			_output_vec[query_id]->_stliner->GatherOutputs(_fw);
			// send write acknowledgement to master
			SendWriteAck(query_id, writer_process);

			if(debug_msg) {
				LOG_MSG << "End gathering outputs for " << query_id <<endl;
			}

			// delete _output_vec[query_id];
			// _output_vec[query_id] = NULL;
			UnloadQuery(query_id);

			_lead_queries.pop();
		} else {
			break;
		}
	}
}

void WriterMasterStreamline::Finalize(void) {
	_comm_send.WaitPendingSends();
	
	WaitAllWriteOps();

	// output html footer
	char* buffer = (char*) malloc (PRINT_BUFFER_LENGTH);
	CHECK_NULPTR(buffer);
	memset(buffer, 0, PRINT_BUFFER_LENGTH);
	_blast->GetHtmlFooter(buffer);
	int print_len = strlen(buffer);
	if(print_len>0) {
		_fw->AddWriteEntry(_curr_output_pos, print_len, buffer, true); // buffer will be freed by FileWriter
	} else {
		free(buffer);
	}
	_fw->Write();
	
	delete _fw;
}

/* Send offsets to all the workers
 * Packing format
   - number of queries
   - query_id
   - number of records
   - description offset, alignment offset
   - (opt) end output position
   - ...
*/
void WriterMaster::SendOutputOffsets(int start_query, int end_query) {
	// group offsets according to the node rank
	vector < vector <MPI_Offset> > offsets_buf(group_node_count);

#ifdef DEBUG
	if(end_query - start_query <=0) {
		throw __FILE__ "try to send empty set of output offsets";
	}
#endif

	if(debug_msg) {
		LOG_MSG << "Start sending output offsets of query " << start_query << endl;
	}
	
	// !!!! only support sending results for one query for now
	
	// number of queries
	for(int k=0; k<group_node_count; k++) {
		offsets_buf[k].push_back(end_query-start_query);
	}
	
	// group offsets according to the node rank
	for(int i=start_query; i<end_query; i++) {
		// index of the position storing number of records
		int index[group_node_count]; 
		// begin pos of the query
		MPI_Offset query_begin = _output_vec[i]->_stliner->GetBeginPos();

		for(int j=0; j<group_node_count; j++) {
			if(!IsWorker(j)) {
				continue;
			}

			offsets_buf[j].push_back(i); // add query_id
			index[j] = offsets_buf[j].size(); // record the position of records count
			offsets_buf[j].push_back(0); // records count
			offsets_buf[j].push_back(query_begin); 
		}
		
		list <COutputRecord>::iterator it;
		// adding records
		for(it=(_output_vec[i]->_sorted_output_list).begin(); it!=(_output_vec[i]->_sorted_output_list).end(); it++) {
			OutputRecordPtr curr_orp = (*it)._orp;
			offsets_buf[curr_orp->rank].push_back(curr_orp->des_offset);
			offsets_buf[curr_orp->rank].push_back(curr_orp->aln_offset);
			offsets_buf[curr_orp->rank][index[curr_orp->rank]]++;
		}

        if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
            for(int j=0; j<group_node_count; j++) {
				if(!IsWorker(j)) {
					continue;
				}

				MPI_Offset output_len = _output_vec[i]->_stliner->GetOutputLength();
				offsets_buf[j].push_back(output_len);
            }
			_output_vec[i]->_stliner->CheckEmpty();
        }
	}
	
	if(output_strategy == WORKER_STREAMLINE) {
		int leader = -1;
		int max_results = -1;
		// pick the worker with most results as write leader
		for(int k=0; k<group_node_count; k++) {
			if(!IsWorker(k)) {
				continue;
			}

			if(_output_vec[start_query]->_stliner->HasWorkerResult(k)) {
				int results_len = offsets_buf[k].size();
				if( results_len > max_results) {
					max_results = offsets_buf[k].size();
					leader = k;
				}
			}
		}

		if(debug_msg) {
			LOG_MSG << "Leader for query " << start_query << " is " << leader << endl;
		}

		_output_vec[start_query]->_stliner->SetWriterLeader(leader);
		_output_vec[start_query]->_stliner->ShiftLeader(start_query, this);
	}
	
	for(int k=0; k<group_node_count; k++) {
		if(!IsWorker(k)) {
			continue;
		}

		if(_output_vec[start_query]->_stliner->HasWorkerResult(k) 
				|| output_strategy == WORKER_COLLECTIVE) {
			
			if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
				_output_vec[start_query]->_stliner->AppendWriteLeader(offsets_buf[k]);
			}

			MPI_Offset msg_size = offsets_buf[k].size()*offset_size;	

			// send event message first
			int assign_array[3];
			assign_array[0] = WRITER_EVENT;
			assign_array[1] = OUTPUT_OFFSETS;
			assign_array[2] = (int)msg_size;
			CommSendStruct* cssp_event = new CommSendStruct(3*int_size);
			cssp_event->AddData(&assign_array[0], 3 * int_size);
			cssp_event->IsendData(group_comm, k, ASSIGNMENT_TYPE);
			_comm_send.AddPendingSend(cssp_event);
		
			CommSendStruct* cssp = new CommSendStruct(msg_size);
			for(int l=0; l<offsets_buf[k].size(); l++) {
				cssp->AddData(&(offsets_buf[k][l]), offset_size);
			}

			cssp->IsendData(group_comm, k, OUTPUT_OFFSETS);
			_comm_send.AddPendingSend(cssp);
		} 
	}

	if(debug_msg) {
		LOG_MSG << "End sending output offsets of query " << start_query << endl;
	}
}

void WriterMaster::ProcessWriteRequest(int src) {
	if(_num_writing >= concurrent_write) {
		// enqueue request
		_write_request_queue.push(src);
	} else {
		_num_writing++;
		
		if(_num_writing > max_num_write) {
			max_num_write = _num_writing;
		}
		
		int command = WORKER_WRITE;
		
		CommSendStruct* cssp_assign = new CommSendStruct(int_size);
		cssp_assign->AddData(&command, int_size);
		cssp_assign->IsendData(group_comm, src, WORKER_WRITE);
		_comm_send.AddPendingSend(cssp_assign);
	}
}

void WriterMaster::ProcessWriteComplete(int src) {
	if(!_write_request_queue.empty()) {
		int dst = _write_request_queue.front();
		
		int command = WORKER_WRITE;
		CommSendStruct* cssp_assign = new CommSendStruct(int_size);
		cssp_assign->AddData(&command, int_size);
		cssp_assign->IsendData(group_comm, dst, WORKER_WRITE);
		_comm_send.AddPendingSend(cssp_assign);

		_write_request_queue.pop();
	} else {
		_num_writing--;
	}
}

int WriterMasterStreamline::CheckEmptyOutput() {

	throw __FILE__ "WriterMasterStreamline::CheckEmptyOutput -- should never be called now";
	
/*
    if(_last_written_query == -1) {
        _last_written_query = 0;
    }

    for(int i=_last_written_query; i<=_query_to_write; i++) {
        if(!_output_vec[i]->_stliner->IsEmpty()) {
            break;
        }

        if(debug_msg) {
            LOG_MSG << "Start gathering outputs for query " << i <<endl;
        }

        _output_vec[i]->_stliner->GatherOutputs(_fw);
        delete _output_vec[i];
        _output_vec[i] = NULL;

        if(debug_msg) {
            LOG_MSG << "End gathering outputs for query " << i <<endl;
        }
        _last_written_query++;
    }

    if(_last_written_query == _query_count) {
        return 1;
    } else {
        return 0;
    }
*/
}

void WriterWorker::AddLocalOutputs(int query_id, int frag_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr) {
	ValNodePtr curr_node = orplist;
	
	// HL-debug-open:
	if(dump_raw_output) {
		fprintf(dbgfp, "Adding output for query %d, fragment %d, num_records=%d\n", query_id, frag_id, _tmp_output_vec[query_id].size());
	}
	
	// add results to tmporary vec, which buffers the output records to be sent to master writer
	while(curr_node) {
		((OutputRecordPtr)curr_node->data.ptrvalue)->rank = group_rank;
		((OutputRecordPtr)curr_node->data.ptrvalue)->frag_id = frag_id;
		_tmp_output_vec[query_id].push_back((OutputRecordPtr)curr_node->data.ptrvalue);
		
		// HL-debug-open:
		if(dump_raw_output) {
			print_output_record(dbgfp, (OutputRecordPtr)curr_node->data.ptrvalue);
		}
		
		curr_node = curr_node->next;
	}
	
	AddOutputRecords(query_id, orplist, qoip, bs_ptr, true);

	// HL-debug-open:
	if(dump_raw_output) {
	    fprintf(dbgfp, "Finish adding output number of records = %d\n", _tmp_output_vec[query_id].size());
    }
	
	_queries_to_processed++;
}

// set the start query for this run
void WriterWorker::SetStartQuery(int start_query) {
	_last_processed_query = start_query;
	_queries_to_processed = 0;
}

// sigleton pattern, initialize _instance to NULL
WriterWorkerIndividual* WriterWorkerIndividual::_instance = NULL;

WriterWorkerStreamline* WriterWorkerStreamline::_instance = NULL;

// singleton pattern, get instance for this class
WriterWorkerStreamline* WriterWorkerStreamline::Instance(BLAST* blast, int query_count, int total_frags, const string& output_file, int send_batch) {
	if(_instance == NULL) {
		_instance = new WriterWorkerStreamline(blast, query_count, total_frags, output_file, send_batch);
	} 
	
	return _instance;
}

WriterWorkerStreamline* WriterWorkerStreamline::Instance() {
	return _instance;
}

void WriterWorkerStreamline::Finalize() {
	if(debug_msg) {
		LOG_MSG << "WriterWorker start finalizaing" <<endl;
	}	
	
	if(_writing_ops.empty()) {
		 _results_send.WaitPendingSends();
	}

	MPI_Status probe_status, recv_status;
	
	int assign_array[3];

	if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		while(!_is_finished) {
			_comm_send.CheckPendingSends();

            int flag = 0;
			int pool_count = 0;

            // use Iproble here to prevent isend from freezing
            MPI_Iprobe(MPI_ANY_SOURCE, ASSIGNMENT_TYPE, group_comm, &flag, &probe_status);
            if(flag) {
				MPI_Recv(assign_array, 3, MPI_INT, probe_status.MPI_SOURCE, ASSIGNMENT_TYPE, group_comm, &recv_status);

                if(assign_array[0] != WRITER_EVENT) {
                    throw __FILE__ "WriterWorkerStreamline::Finalize -- illegal assignment";
                }

                HandleMessages(probe_status, &(assign_array[1]));
            } else {
				pool_count++;

				if(pool_count > 100) {
					sleep(1);
					pool_count = 0;
				}
            }

            if(!_writing_ops.empty()) {
                ProcessAsyncOutput();
				WaitAllWriteOps();
            }

			if(output_strategy == MASTER_STREAMLINE) {
				if( _last_processed_query == 0 || _last_output_query == _last_processed_query) {
					_is_finished = true;
				}
			}
		}
	}
	
	_comm_send.WaitPendingSends();

	if(!_writing_ops.empty()) {
		throw __FILE__ "WriterWorkerStreamline::Finalize -- unwritten queries left";
	}

	int event_array[2];
	event_array[0] = WRITER_FINALIZED;
	event_array[1] = int_size;
	CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
	cssp_event->AddData(&event_array[0], 2*int_size);
	cssp_event->SendData(group_comm, writer_process, EVENT_TYPE);
	delete cssp_event;

	// only need a tag here, the content is not important
	MPI_Send(event_array, 1, MPI_INT, writer_process, WRITER_FINALIZED, group_comm);

	delete _fw;

	if(debug_msg) {
		LOG_MSG << "WriterWorker end finalizaing" <<endl;
	}	
}

/* send scores and sizes to the master 
 * The format of packed message:
	- packing type
	- fragment id
	- number of queries results
	- query_id, number of records
	- description header length
	- description header data
	- description footer length
	- description footer data
	- stats data size
	- stats data 
	- actual records 
	- query_id, number of records
	- actual records 
	- ......
 */ 
void WriterWorker::SendPartialOutputs(int start_query, int end_query, int pack_type, int frag_id) {
    if(debug_msg) {
        LOG_MSG << "Start sending partial outputs frag " << frag_id << " from query "<<start_query<<" to query "<<end_query<<endl;
    }

	static int blast_align_view = _blast->GetAlignView();
	static bool blast_html = _blast->GetHtml();
	
	MPI_Offset msg_size=0;
	// calculate the size of message
		// header size
	msg_size += int_size  // pack type
			  + int_size  // fragment id
			  + int_size; // number of queries results

	// HL-debug-open:
	if(dump_raw_output) {
		fprintf(dbgfp, "Sending output for query %d, fragment %d, num_records=%d\n", start_query, frag_id, _tmp_output_vec[start_query].size());
	}

	for(int i=start_query; i<end_query; i++) {

		// if(dump_raw_output) { // print accumulated results for this query
		//	_output_vec[i]->DebugPrint(dbgfp);
		//}

		if(blast_align_view < 10) { // not asn
			msg_size += int_size  // query_id
					+ int_size  // number of records
					+ int_size  // description header size
					+ _output_vec[i]->_qoip->des_hdr_size
					+ int_size  // description footer size
					+ _output_vec[i]->_qoip->des_ftr_size;
		} else { // asn
			msg_size += int_size  // query_id
					+ int_size  // number of records
					+ int_size  // header size
					+ _output_vec[i]->_qoip->header_size
					+ int_size  // footer size
					+ _output_vec[i]->_qoip->footer_size;
		}
		
		msg_size += int_size 	// _stats_buf_size
				+ _output_vec[i]->_stats_buf_size; 
		
		for(int j=0; j<_tmp_output_vec[i].size(); j++) {
			msg_size += get_send_record_size(pack_type, _tmp_output_vec[i][j]);
		}
	}

	// package message
	CommSendStruct* cssp=new CommSendStruct(msg_size);
	cssp->AddData(&pack_type, int_size);
	cssp->AddData(&frag_id, int_size);
	int num_queries = end_query - start_query;
	cssp->AddData(&num_queries, int_size);

	for(int i=start_query; i<end_query; i++) {
		cssp->AddData(&i, int_size); // query_id
		int num_records = _tmp_output_vec[i].size();
		cssp->AddData(&num_records, int_size); // number of records
		
		_pending_offsets.insert(start_query); // duplication will be handled by set insert

		if(_pending_offsets.size() > peak_pending_offsets) {
			peak_pending_offsets = _pending_offsets.size();
		}

		int tmp_size=0;
		if(blast_align_view < 10) {// not asn
			//description header and footer are stored in _output_vec
			tmp_size = _output_vec[i]->_qoip->des_hdr_size;
			cssp->AddData(&tmp_size, int_size);
			if(tmp_size > 0) {
				cssp->AddData(_output_vec[i]->_qoip->des_hdr, tmp_size);
			}
			tmp_size = _output_vec[i]->_qoip->des_ftr_size;
			cssp->AddData(&tmp_size, int_size);
			if(tmp_size > 0) {
				cssp->AddData(_output_vec[i]->_qoip->des_ftr, tmp_size);
			}
		} else {
			//header and footer 
			tmp_size = _output_vec[i]->_qoip->header_size;
			cssp->AddData(&tmp_size, int_size);
			if(tmp_size > 0) {
				cssp->AddData(_output_vec[i]->_qoip->header, tmp_size);
			}
			tmp_size = _output_vec[i]->_qoip->footer_size;
			cssp->AddData(&tmp_size, int_size);
			if(tmp_size > 0) {
				cssp->AddData(_output_vec[i]->_qoip->footer, tmp_size);
			}
		}

		// stats buf
		tmp_size = _output_vec[i]->_stats_buf_size;
		cssp->AddData(&tmp_size, int_size);
		if(tmp_size > 0) {
			cssp->AddData(_output_vec[i]->_stats_buf, tmp_size);
			Nlm_MemFree(_output_vec[i]->_stats_buf);
			_output_vec[i]->_stats_buf = NULL;
		}
		
		for(int j=0; j<_tmp_output_vec[i].size(); j++) {
			pack_record(pack_type, _tmp_output_vec[i][j], cssp);

			// HL-debug-open:
			if(dump_raw_output) {
				print_output_record(dbgfp, _tmp_output_vec[i][j]);
			}
		}
	}

	if(output_strategy == MASTER_STREAMLINE) {
		WaitAllWriteOps();
	}
	
	int tag = PARTIAL_OUTPUTS_DATA;
	// send event first
	int event_array[2];
	event_array[0] = tag;
	event_array[1] = (int)(cssp->GetSize());
	CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
	cssp_event->AddData(&event_array[0], 2*int_size);

	if(sync_comm) {
		double track_time = MPI_Wtime();
		cssp_event->SendData(group_comm, writer_process, EVENT_TYPE);
		cssp->SendData(group_comm, writer_process, tag);
		delete cssp_event;
		delete cssp;
		wait_result_time += MPI_Wtime() - track_time;
	} else {
		// wait previous sending finished
		if(debug_msg) {
			LOG_MSG << "Before wait previous sent results" << endl;
		}
		
		double track_time = MPI_Wtime();
		_results_send.WaitPendingSends();
		wait_result_time += MPI_Wtime() - track_time;

		if(debug_msg) {
			LOG_MSG << "After wait previous sent results" << endl;
		}

		cssp_event->IsendData(group_comm, writer_process, EVENT_TYPE);
		_results_send.AddPendingSend(cssp_event);
	
		cssp->IsendData(group_comm, writer_process, tag);
		_results_send.AddPendingSend(cssp);
	}

	// prune unqualified output for this query
    for(int i=start_query; i<end_query; i++) {
		_output_vec[i]->PruneOutputList();
	}

    // clean up tmp send vector
    for(int i=start_query; i<end_query; i++) {
        _tmp_output_vec[i].clear();
    }

    if(debug_msg) {
        LOG_MSG << "End sending partial outputs from query "<<start_query<<" to query "<<end_query<<endl;
    }
}

WriterWorkerCollective* WriterWorkerCollective::_instance = NULL;

// process query results, send size and scores to master writer
void WriterWorker::ProcessResults(bool is_search_end, int frag_id) {
	int start_query, end_query;
	int pack_type;
	
	if(_queries_to_processed > 0 && ((_queries_to_processed >= _send_batch) || is_search_end)) {
	
		start_query = _last_processed_query;
		end_query = _last_processed_query + _queries_to_processed;
		
		static int blast_align_view = _blast->GetAlignView();
		static bool blast_html = _blast->GetHtml();
		
		switch(blast_align_view) {
		case 0:
			pack_type = SEND_TYPE1;
			break;
		case 7:
		case 8:
		case 9:
		case 10:
		case 11:
			pack_type = SEND_TYPE2;
			break;
		}
		
		SendPartialOutputs(start_query, end_query, pack_type, frag_id);
		_last_processed_query = end_query;
		_queries_to_processed = 0;
	}
}

// receive offsets from the master
// after offsets are received, the data to be outputted are passed to the FileWriter for writing
/* Packing format
   - number of queries
   - query_id
   - number of records
   - description offset, alignment offset
   - ...
*/
void WriterWorker::ReceiveOutputOffsets(CommRecvStruct* crsp) {
	if(debug_msg) {
		LOG_MSG << "Start receiving output offsets" <<endl;
	}
	
	// get number of queries
	MPI_Offset num_queries;
	crsp->ExtractData(&num_queries, offset_size);
	
	for(int i=0; i<num_queries; i++) {
		MPI_Offset query_id;
		crsp->ExtractData(&query_id, offset_size);

		if(debug_msg) {
			LOG_MSG << "Output offsets of query " << query_id <<endl;
		}

		InitQueryOutput(0, query_id);

		_pending_offsets.erase(query_id);

		MPI_Offset num_records;
		crsp->ExtractData(&num_records, offset_size);
		MPI_Offset query_begin = -1;
		crsp->ExtractData(&query_begin, offset_size);

		_output_vec[query_id]->_stliner->SetBeginPos(query_begin);
		
		list <COutputRecord>::iterator it;
		it=_output_vec[query_id]->_sorted_output_list.begin();
		for(int j=0; j<num_records; j++, it++) {
			OutputRecordPtr curr_orp = (*it)._orp;
			crsp->ExtractData(&(curr_orp->des_offset), offset_size);
			crsp->ExtractData(&(curr_orp->aln_offset), offset_size);

			if(curr_orp->des_offset == -1 && curr_orp->aln_offset == -1) {
				throw __FILE__ "WriterWorker::ReceiveOutputOffsets: received illegal offset values";
			}
		}

        if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			MPI_Offset output_len = -1;
			crsp->ExtractData(&output_len, offset_size);
			_output_vec[query_id]->_stliner->SetOutputLength(output_len);

			if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
				MPI_Offset leader; // notice it is offset type
				crsp->ExtractData(&leader, offset_size);
				int int_leader = leader; // should use type cast?
				_output_vec[query_id]->_stliner->SetWriterLeader(int_leader);
			}
        }

		// prune the unqualified local outputs
		_output_vec[query_id]->_sorted_output_list.erase(it, _output_vec[query_id]->_sorted_output_list.end());

		if(debug_msg) {
			LOG_MSG << "Start preparing output data for query " << query_id <<endl;
		}
		
		_output_vec[query_id]->PrepareOutputData(_blast, _fw);

		if(debug_msg) {
			LOG_MSG << "End preparing output data" <<endl;
		}
		
        if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			// wait for previous write to finish
			if( output_strategy == MASTER_STREAMLINE || 
					(output_strategy == WORKER_STREAMLINE && _writing_ops.size() > num_pending_writes)) {
				WaitAllWriteOps();
			}

			_output_vec[query_id]->_stliner->AlignOutputEntries();

			int num_write_segments = _output_vec[query_id]->_stliner->GetNumWriteSegments();
			
			AddWriteOps(query_id, num_write_segments);
			_output_vec[query_id]->_stliner->PrepareCommWriteData();

			// send local output to the write leader
			_output_vec[query_id]->_stliner->SendOutputData(-1, query_id, this);

			if(num_write_segments == 0) {
				if(output_strategy == WORKER_STREAMLINE) {
					if(_output_vec[query_id]->_stliner->IsLeader()) {
						SendWriteAck(query_id, writer_process);
					}			
				}
				
				_last_output_query = query_id + 1;
				_last_written_query = query_id + 1;
			}

		} else if (output_strategy == WORKER_COLLECTIVE) {
			// tell writer master I am ready to write
			int event_array[2];
			event_array[0] = COLLECTIVE_READY;
			event_array[1] = int_size;
			CommSendStruct* cssp_event = new CommSendStruct(2*int_size);
			cssp_event->AddData(&event_array[0], 2*int_size);
			cssp_event->SendData(group_comm, writer_process, EVENT_TYPE);
			delete cssp_event;

			// don't send query_id, it has MPI_Offset type
			event_array[0] = query_id;
			MPI_Send(event_array, 1, MPI_INT, writer_process, COLLECTIVE_READY, group_comm);

			if(dump_raw_output) {
				_fw->DumpWriteRequests(dbgfp);
			}
			
			// collectively write data
			_fw->CollWrite();

			_last_written_query = query_id + 1;
		}

		if(output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
			_output_vec[query_id]->CleanupOutputs();

			// !!!! this if block seems never to be true 
			if(_output_vec[query_id]->_stliner->IsProcessed()) {
				// clean up the buffered output for this query
				// delete _output_vec[query_id];
				// _output_vec[query_id] = NULL;
				UnloadQuery(query_id);
				
				_last_output_query = query_id + 1;
			}
		} else {
			// clean up the buffered output for this query
			//  delete _output_vec[query_id];
			// _output_vec[query_id] = NULL;
			UnloadQuery(query_id);

			_last_output_query = query_id + 1;
		}
	}
	
	if(debug_msg) {
		LOG_MSG << "End receiving output offsets" <<endl;
	}
}

void Streamliner::SendOutputData(int dst, int query_id, Writer* wtrp) {

	if(dst == -1) {
		dst = _write_leader;
	}

	if(dst == group_rank) { // write leader is myself
		if(_workers_with_output.find(group_rank) != _workers_with_output.end()) {
			_workers_with_output.erase(group_rank);
		}

		return;
	}

	if(debug_msg) {
		LOG_MSG << "Start sending output data of query "<<query_id << " to " << dst <<endl;
	}

	for(int i=0; i<_num_write_segments; i++) {
		if(!_is_empty) {
			SendDataForAWrite(dst, query_id, i);
		}
		_write_status[i] = 0;
	}

	if(debug_msg) {
		LOG_MSG << "End sending output data of query "<<query_id<<endl;
	}
}


void Streamliner::SendDataForAWrite(int dst, int query_id, int write_segment_id) {
	MPI_Offset buf_size = 0;
	int num_records = 0;

	MPI_Offset write_end = (write_segment_id + 1) * _max_write;

	// calculate buffer size
	buf_size += int_size; // query_id
	buf_size += int_size; // num_records

	map <MPI_Offset, OutputEntry*>::iterator mit;

	if(_output_entries.size() > 0) {
		for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
			MPI_Offset begin_pos = mit->first;
			OutputEntry* oep = mit->second;

			if(begin_pos + oep->_out_size > write_end) {
				break;
			}

			buf_size += offset_size;
			buf_size += int_size;
			buf_size += oep->_out_size;	
			num_records++;
		}	
	}

	CommSendStruct* cssp;	
	cssp = new CommSendStruct(buf_size);
	cssp->AddData(&query_id, int_size);
	cssp->AddData(&num_records, int_size);
	
	if(num_records > 0) {
		for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
			MPI_Offset begin_pos = mit->first;
			OutputEntry* oep = mit->second;

			if(begin_pos + oep->_out_size > write_end) {
				break;
			}

			cssp->AddData(&(oep->_out_off), offset_size);
			cssp->AddData(&(oep->_out_size), int_size);
			cssp->AddData(oep->_out_data, oep->_out_size);
			free(oep->_out_data);
			delete oep;
		}	
		_output_entries.erase(_output_entries.begin(), mit);
	}

	cssp->IsendData(group_comm, dst, DATA_TO_WRITE);
	_isend_list[write_segment_id].push_back(cssp);
}

void Streamliner::SendDataForAWriteEx(int dst, int query_id, Writer* wtrp, MPI_Offset write_end, vector<CommSendStruct*>& msg_send_list) {
	MPI_Offset buf_size = 0;
	int num_records = 0;

	// calculate buffer size
	buf_size += int_size; // query_id
	buf_size += int_size; // num_records

	map <MPI_Offset, OutputEntry*>::iterator mit;

	if(_output_entries.size() > 0) {
		for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
			MPI_Offset begin_pos = mit->first;
			OutputEntry* oep = mit->second;

			if(begin_pos + oep->_out_size > write_end) {
				break;
			}

			buf_size += offset_size;
			buf_size += int_size;
			buf_size += oep->_out_size;	
			num_records++;
		}	
	}

	CommSendStruct* cssp2;	
	cssp2 = new CommSendStruct(buf_size);
	cssp2->AddData(&query_id, int_size);
	cssp2->AddData(&num_records, int_size);
	
	if(num_records > 0) {
		for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
			MPI_Offset begin_pos = mit->first;
			OutputEntry* oep = mit->second;

			if(begin_pos + oep->_out_size > write_end) {
				break;
			}

			cssp2->AddData(&(oep->_out_off), offset_size);
			cssp2->AddData(&(oep->_out_size), int_size);
			cssp2->AddData(oep->_out_data, oep->_out_size);
			free(oep->_out_data);
			delete oep;
		}	
		_output_entries.erase(_output_entries.begin(), mit);
	}

	msg_send_list.push_back(cssp2);
}

void Streamliner::ExtractOutputEntries(CommRecvStruct* crsp) {
	MPI_Offset output_key;
	int num_records;
	int record_size;

	crsp->ExtractData(&num_records, int_size);

	if(num_records == 0) {
		return;
	}

	//debug
	//fprintf(dbgfp, "num_records = %ld\n", num_records);
	//fflush(dbgfp);
	
	for(int i=0; i<num_records; i++) {
		MPI_Offset out_off = -1;
		crsp->ExtractData(&out_off, offset_size);
		int out_size = -1;
		crsp->ExtractData(&out_size, int_size);
		char* out_data = (char*)malloc(out_size); // this memory will be freed in FileWriter
		CHECK_NULPTR(out_data);
		crsp->ExtractData(out_data, out_size);

		AddOutputEntry(out_off, out_size, out_data, false, false);
	}		
}

MPI_Offset Streamliner::AddOutputEntry(MPI_Offset& out_off, int out_size, char* out_data, bool update_off, bool mark_write) {

	_output_entries[out_off] = new OutputEntry(out_off, out_size, out_data);
	if(update_off) {
		//update current output position
		out_off += out_size;
	}

	return 0;
}

// append write leader to the offsets buffer to be sent to workers
void Streamliner::AppendWriteLeader(vector <MPI_Offset> &offset_buf_vec) {
	offset_buf_vec.push_back(_write_leader);
}

inline bool Streamliner::IsWritten() {
	return _is_written;
}

inline bool Streamliner::IsReady() {
	return _is_ready;
}

inline bool Streamliner::IsProcessed() {
	return _is_processed;
}

inline void Streamliner::AddWorkerWOutputs(const set<int>& workers) {
	if(debug_msg) {
		ostringstream os;
		copy(workers.begin(), workers.end(), ostream_iterator<int>(os<<"Workers with output: ", " "));
		LOG_MSG << os.str() << endl;
	}

	_workers_with_output.insert(workers.begin(), workers.end());
}

void Streamliner::ExtractWritePeers(CommRecvStruct* crsp) {
	int num_write_peers;
	crsp->ExtractData(&num_write_peers, int_size);

	vector <int> peers;
	for(int i=0; i<num_write_peers; i++) {
		int worker_id;
		crsp->ExtractData(&worker_id, int_size);
		_workers_with_output.insert(worker_id);
		peers.push_back(worker_id);
	}


	// data messages
	crsp->ExtractData(&_num_write_segments, int_size);

	if(num_write_peers > 0) {
		for(int i=0; i<peers.size(); i++) {
			for(int j=0; j<_num_write_segments; j++) {
				int msg_size;
				crsp->ExtractData(&msg_size, int_size);
				_data_msg_sizes[peers[i]].push_back(msg_size);
			}
		}
	}
}

void Streamliner::ExtractSubmitters(CommRecvStruct* crsp) {
	int num_submitters;
	crsp->ExtractData(&num_submitters, int_size);

	for(int i=0; i<num_submitters; i++) {
		int worker_id;
		crsp->ExtractData(&worker_id, int_size);
		_submit_workers.insert(worker_id);
	}
}

void Streamliner::WriteCurrentSegment(FileWriter *fwp) {
	MPI_Offset curr_write_end = (_curr_write_segment + 1) * _max_write;
	WriteOutputEntries(fwp, curr_write_end);
	fwp->Write();

	_curr_write_segment++;
	_is_ready = false;
	_complete_irecv = 0;

	if(_curr_write_segment == _num_write_segments) {
		// finish output for this query
		_is_written = true;
	}
}
	
void Streamliner::WriteASegment(int write_segment_id, FileWriter *fwp) {
	MPI_Offset curr_write_end = (write_segment_id + 1) * _max_write;
	WriteOutputEntries(fwp, curr_write_end);
	fwp->Write();

	_curr_write_segment = write_segment_id + 1;
	_is_ready = false;
	_complete_irecv = 0;
}

void Streamliner::WriteOutputEntries(FileWriter *fwp, MPI_Offset write_end) {
	map <MPI_Offset, OutputEntry*>::iterator mit;
	for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
		OutputEntry* oep = mit->second;
		MPI_Offset entry_begin = mit->first;
		MPI_Offset end_pos = entry_begin + oep->_out_size;
		
		if(end_pos <= write_end) {
			MPI_Offset real_begin = entry_begin + _begin_pos;
			fwp->AddWriteEntry(real_begin, oep->_out_size, oep->_out_data, false);
			delete oep;
		} else {
			break;
		}
	}

	_output_entries.erase(_output_entries.begin(), mit);
}

void Streamliner::GatherOutputs(FileWriter* fwp) {

	double track_time = MPI_Wtime();
	
	if(_is_written) {
		throw __FILE__ "Streamliner::GatherOutputs -- query has been written";
	}

	MPI_Offset curr_write_end = 0;
	while(curr_write_end < _output_len) {
		curr_write_end += _max_write;
		GatherOutputsForAWrite(fwp);
		WriteOutputEntries(fwp, curr_write_end);	
		fwp->Write();
	}

	_is_written = true;
	
	results_gather_time += MPI_Wtime() - track_time;
}

void Streamliner::GatherOutputsForAWrite(FileWriter* fwp) {
    int num_nodes;
    num_nodes = _pend_recv_nodes.size();

    int buf_size;
    for(int i=0; i<num_nodes; i++) {
        // receive actual data to write
        MPI_Status probe_status;
        MPI_Probe( _pend_recv_nodes[i], DATA_TO_WRITE, group_comm, &probe_status );
        MPI_Get_count(&probe_status, MPI_BYTE, &buf_size);

        CommRecvStruct* crsp = new CommRecvStruct(buf_size);
        crsp->RecvData(group_comm, _pend_recv_nodes[i], DATA_TO_WRITE);
        int curr_query;
        crsp->ExtractData(&curr_query, int_size);

        // debug
        //fprintf(dbgfp, "received from nodes %d\n", _pend_recv_nodes[i]);
        //fprintf(dbgfp, "current_query = %d\n", curr_query);
        //fprintf(dbgfp, "buf_size = %d\n", buf_size);

        ExtractOutputEntries(crsp);
        delete crsp;
    }
}

int Streamliner::AddPendingRecvs(int src) {
	_pend_recv_nodes.push_back(src);

	return 0;
}

void Streamliner::CheckEmpty() {
	if(_workers_with_output.size() == 0) {
		_is_ready = true;
		_is_empty = true;
	}
}

void Streamliner::ShiftLeader(int query_id, Writer* wtrp) {
	
	int dst = _write_leader;

	// send all output entries, workers that has results
	// should implemented as a method of streamliner
	
	// format
	//   query_id
	//   write end posistion
	//   number of output entries
	//     offset, size, buffer
	//     offset, size, buffer
	//     ....
	//   number of workers with outputs
	//     rank1
	//     rank2
	//     ....
	
	// send event message first
	// calculate buffer size
	MPI_Offset buf_size = 0;
	int num_entries = 0;
	
	buf_size += int_size; // query_id
	buf_size += offset_size; // write begin position
	buf_size += offset_size; // output length of this query
	buf_size += int_size; // num entries

	// expected message sizes

	map <MPI_Offset, OutputEntry*>::iterator mit;
	for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
		MPI_Offset begin_pos = mit->first;
		OutputEntry* oep = mit->second;

		buf_size += offset_size;
		buf_size += int_size;
		buf_size += oep->_out_size;	
		num_entries++;
	}	

	buf_size += int_size; // num workers submitting results
	buf_size += int_size * _submit_workers.size();
	
	buf_size += int_size; // num workers with data
	buf_size += int_size * _workers_with_output.size();

	// data msg sizes
	int num_write_segments = _num_write_segments;
	
	buf_size += int_size; // num write segments

	int num_workers = _workers_with_output.size();
	buf_size += num_workers * num_write_segments * int_size;

	int assign_array[3];
	assign_array[0] = WRITER_EVENT;
	assign_array[1] = SHIFT_WRITE_LEADER;
	assign_array[2] = (int)buf_size;
	CommSendStruct* cssp_event = new CommSendStruct(3*int_size);
	cssp_event->AddData(&assign_array[0], 3 * int_size);
	cssp_event->IsendData(group_comm, dst, ASSIGNMENT_TYPE);
	wtrp->AddCommSend(cssp_event);
	
	CommSendStruct* cssp;	
	cssp = new CommSendStruct(buf_size);
	cssp->AddData(&query_id, int_size);
	cssp->AddData(&_begin_pos, offset_size);
	cssp->AddData(&_output_len, offset_size);
	cssp->AddData(&num_entries, int_size);
	
	for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
		MPI_Offset begin_pos = mit->first;
		OutputEntry* oep = mit->second;

		cssp->AddData(&(oep->_out_off), offset_size);
		cssp->AddData(&(oep->_out_size), int_size);
		cssp->AddData(oep->_out_data, oep->_out_size);
		free(oep->_out_data);
		delete oep;
	}	

	set <int>::iterator sit;

	int num_submit = _submit_workers.size();
	cssp->AddData(&num_submit, int_size);
	for(sit=_submit_workers.begin(); sit!=_submit_workers.end(); sit++) {
		cssp->AddData(&(*sit), int_size);
	}

	cssp->AddData(&num_workers, int_size);
	for(sit=_workers_with_output.begin(); sit!=_workers_with_output.end(); sit++) {
		cssp->AddData(&(*sit), int_size);
	}

	cssp->AddData(&num_write_segments, int_size);

	// data message sizes
	if(num_workers > 0) {
		for(sit=_workers_with_output.begin(); sit!=_workers_with_output.end(); sit++) {
			int wrk_id = *sit;

			for(int i=0; i<num_write_segments; i++) {
				cssp->AddData(&(_data_msg_sizes[wrk_id][i]), int_size);
			}
		}
	}

	cssp->IsendData(group_comm, dst, SHIFT_WRITE_LEADER);
	wtrp->AddCommSend(cssp);

	_output_entries.clear();
}

bool Streamliner::HasWorkerResult(int src) {
	set <int>::iterator sit;
	sit = _submit_workers.find(src);
	if (sit != _submit_workers.end()) {
		return true;
	}
	return false;
}

/*
void Streamliner::UpdateRealOffsets(MPI_Offset& curr_output_pos) {
	
	map <MPI_Offset, OutputEntry*>::iterator mit;
	map <MPI_Offset, OutputEntry*> tmp_entries;
	
	for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
		MPI_Offset curr_off = mit->first;
		OutputEntry* oep = mit->second;
		curr_off += curr_output_pos;
		oep->_out_off += curr_output_pos;
		tmp_entries[curr_off] = oep;
	}

	_output_entries.clear();
	_output_entries = tmp_entries;
	
	curr_output_pos += _output_len;
}
*/

void Streamliner::PrepareCommWriteData() {
	// we are actually compute ceiling here

	if(_output_entries.empty()) {
		_is_empty = true;
	}

	for(int i=0; i<_num_write_segments; i++) {
		_irecv_list.push_back(list <CommRecvStruct*>());
		_write_status.push_back(-1);
		_isend_list.push_back(list <CommSendStruct*>());
	}

	_curr_write_segment = 0;
	_complete_irecv = 0;
}

void Streamliner::InitAWrite(int query_id, int write_segment_id) {
	_query_id = query_id;
	
	if(_write_status[write_segment_id] > -1) { // already initialized
		return;
	}

	if(_write_leader == group_rank) { // myself is leader, issue irecvs
		set <int>::iterator sit;

		if(debug_msg) {
			LOG_MSG << "Initialize write for query " << query_id << " segment " << write_segment_id << " with " << _workers_with_output.size() << " irecvs" << endl;
		}

		for(sit=_workers_with_output.begin(); sit!=_workers_with_output.end(); sit++) {
			int wrk_id = *sit;
			
			if(wrk_id == group_rank) {
				continue;
			}

			if(_data_msg_sizes.find(wrk_id) == _data_msg_sizes.end()) {
				throw __FILE__ "Streamliner::IssueDataIrecv -- invalid message size";
			}
			if(_data_msg_sizes[wrk_id][_curr_write_segment] > -1) {
				// issue irecv
				CommRecvStruct* crsp = new CommRecvStruct(_data_msg_sizes[wrk_id][_curr_write_segment]);
				crsp->IrecvData(group_comm, wrk_id, DATA_TO_WRITE);
				_irecv_list[write_segment_id].push_back(crsp);

				// mark irecv has been issued
				_data_msg_sizes[wrk_id][_curr_write_segment] = -1;
			}
		}
	}

	if(group_rank != writer_process) {
		if(_write_leader == group_rank) { // write leader is myself
			if(_workers_with_output.find(group_rank) != _workers_with_output.end()) {
				_workers_with_output.erase(group_rank);
			}
		} else {
			if(!_is_empty) {
				SendDataForAWrite(_write_leader, query_id, write_segment_id);
			}
		}
	}

	_write_status[write_segment_id] = 0;
}

// return 1 when current write is finished
// return 0 otherwise
int Streamliner::CheckAWrite(int write_segment_id, FileWriter* fwp) {

	// check status
	if(_write_status[write_segment_id] == -1) {
		throw __FILE__ "Streamliner::CheckAWrite - segment has not been initialized";
	} 

	if(_write_status[write_segment_id] == 1) {
		return 1;
	}

	if(!_isend_list[write_segment_id].empty()) {
		// check isend
		list <CommSendStruct*>::iterator it;
		list <CommSendStruct*>::iterator tmp_it;
		
		it=_isend_list[write_segment_id].begin();
		while(it!=_isend_list[write_segment_id].end()) {
			if((*it)->TestIsend()) {
				delete(*it);
				tmp_it = it;
				it++;
				_isend_list[write_segment_id].erase(tmp_it);
			} else {
				it++;
			}
		}
	}

	if(!_irecv_list[write_segment_id].empty()) {
		// check irecv
		list <CommRecvStruct*>::iterator it;
		list <CommRecvStruct*>::iterator tmp_it;
		
		it=_irecv_list[write_segment_id].begin();
		while(it!=_irecv_list[write_segment_id].end()) {
			if((*it)->TestIrecv()) {
				int query_id;
				(*it)->ExtractData(&query_id, int_size);
				ExtractOutputEntries(*it);

				delete(*it);
				tmp_it = it;
				it++;
				_irecv_list[write_segment_id].erase(tmp_it);
			} else {
				it++;
			}
		}
	} 

	if( _irecv_list[write_segment_id].empty() && _isend_list[write_segment_id].empty()) {
		if(_write_leader == group_rank) {
			// received all output data
			// ready for a write
			_is_ready = true;

			if(debug_msg) {
				LOG_MSG << "Writing a segment" << endl;
			}

			WriteASegment(write_segment_id, fwp);
		} 

		_write_status[write_segment_id] = 1; // finished processing

		if(write_segment_id == _num_write_segments - 1) {
			_is_written = true;
		}

		return 1;
	}

	return 0;
}

int Streamliner::WaitAWrite(int write_segment_id, FileWriter* fwp) {
	
	double track_time = MPI_Wtime();

	if(_write_status[write_segment_id] == 1) { // finished
		return 0;
	}

	vector <MPI_Request> requests;

	if(!_isend_list[write_segment_id].empty()) {

		// debug !!!!
		//if(_isend_list[write_segment_id].size() > 1) {
		//	throw __FILE__ "Streamliner::WaitAWrite - invalid send list";
		//}

		// check isend
		list <CommSendStruct*>::iterator lit;
		for(lit=_isend_list[write_segment_id].begin(); lit!=_isend_list[write_segment_id].end(); lit++) {
			requests.push_back((*lit)->GetReq());
		}
	}

	if(!_irecv_list[write_segment_id].empty()) {
		// check irecv
		list <CommRecvStruct*>::iterator lit;
		for(lit=_irecv_list[write_segment_id].begin(); lit!=_irecv_list[write_segment_id].end(); lit++) {
			requests.push_back((*lit)->GetReq());
		}
	}

	int count = requests.size();

	// count can be zero for queries with empty output
	if(count > 0) {
		MPI_Request *req_array;
		req_array=(MPI_Request *)malloc(count * sizeof(MPI_Request));
		CHECK_NULPTR(req_array);

		copy(requests.begin(), requests.end(), req_array);

		if(debug_msg) {
			LOG_MSG << "Wait for " << count << " pending non-blocking MPI communication ops" << endl;
		}

		MPI_Waitall(count, req_array, MPI_STATUSES_IGNORE);

		if(!_isend_list[write_segment_id].empty()) {
			// check isend
			list <CommSendStruct*>::iterator it;
			
			it=_isend_list[write_segment_id].begin();
			while(it!=_isend_list[write_segment_id].end()) {
				delete(*it);
				it++;
			}
			_isend_list[write_segment_id].clear();
		}

		if(!_irecv_list[write_segment_id].empty()) {
			// check irecv
			list <CommRecvStruct*>::iterator it;
			
			it=_irecv_list[write_segment_id].begin();
			while(it!=_irecv_list[write_segment_id].end()) {
				int query_id;
				(*it)->ExtractData(&query_id, int_size);
				ExtractOutputEntries(*it);

				delete(*it);
				it++;
			}
			_irecv_list[write_segment_id].clear();
		} 
	}
		
	if(_write_leader == group_rank) {
		// received all output data
		// ready for a write
		_is_ready = true;

		if(debug_msg) {
			LOG_MSG << "Writing a segment" << endl;
		}

		WriteASegment(write_segment_id, fwp);
	}

	_write_status[write_segment_id] = 1; // finished processing

	if(write_segment_id == _num_write_segments - 1) {
		_is_written = true;
	}

	wait_write_time += MPI_Wtime() - track_time;

	return 1;
}

// split output entries if accross writer buffer boundary
void Streamliner::AlignOutputEntries() {

	if(_output_entries.empty()) {
		return;
	}

	MPI_Offset write_end = 0;

	MPI_Offset last_pos = 0;
	while(write_end <= _output_len) {

		write_end += _max_write;

		map <MPI_Offset, OutputEntry*>::iterator mit;
		
		if(last_pos > 0) {
			mit = _output_entries.find(last_pos);
			if(mit == _output_entries.end()) {
				throw __FILE__ "Streamliner::AlignOutputEntries - invalid last pos";
			}
		} else {
			mit = _output_entries.begin();
		}

		for(; mit!=_output_entries.end(); mit++) {
			MPI_Offset begin_pos = mit->first;
			OutputEntry* oep = mit->second;
			last_pos = begin_pos;
			
			if(begin_pos + oep->_out_size > write_end) {
				if(begin_pos < write_end) { // accross boundary

					// split the entry
					MPI_Offset offset1, offset2;
					int size1, size2;
					char* out_data1;
				    char* out_data2;

					offset1 = begin_pos;
					offset2 = write_end;

					size1 = write_end - offset1;
					size2 = oep->_out_size - size1;

					out_data1 = (char*)malloc(size1);
					CHECK_NULPTR(out_data1);
					out_data2 = (char*)malloc(size2);
					CHECK_NULPTR(out_data2);
					
					memcpy(out_data1, oep->_out_data, size1);
					memcpy(out_data2, oep->_out_data + size1, size2);

					free(oep->_out_data);
					delete oep;
					_output_entries.erase(mit);

					_output_entries[offset1] = new OutputEntry(offset1, size1, out_data1);
					_output_entries[offset2] = new OutputEntry(offset2, size2, out_data2);
				}
				break;
			}
		}
	}
}

void Streamliner::CalculateMessageSizes() {

	if(_workers_with_output.empty()) {
		return;
	}

	int num_msgs = _num_write_segments;
	assert(num_msgs > 0);

	map <int, list< pair <MPI_Offset, int> > >::iterator mit;

	for(mit=_access_locations.begin(); mit!=_access_locations.end(); mit++) {
		int worker = mit->first;

		MPI_Offset write_end = 0;
		while(write_end <= _output_len) {
			
			write_end += _max_write;

			int msg_size = 0;
			msg_size += int_size; // query_id
			msg_size += int_size; // num_records
			
			list< pair<MPI_Offset, int> >::iterator lit, tmp_lit;
			
			lit = _access_locations[worker].begin();
			while(lit != _access_locations[worker].end()) {
				MPI_Offset offset = lit->first;
				int size = lit->second;

				if(offset + size <= write_end) {
					// count in msg_size
					msg_size += offset_size;
					msg_size += int_size;
					msg_size += size;

					tmp_lit = lit;
					lit++;

					_access_locations[worker].erase(tmp_lit);
				} else {
					// check if it needs to be splitted
					if(offset < write_end) {
						MPI_Offset offset1, offset2;
						int size1, size2;

						offset1 = offset;
						offset2 = write_end;

						size1 = write_end - offset1;
						size2 = size - size1;

						msg_size += offset_size;
						msg_size += int_size;
						msg_size += size1;

						tmp_lit = lit;
						lit++;
						_access_locations[worker].erase(tmp_lit);

						_access_locations[worker].push_front(pair<MPI_Offset, int>(offset2, size2));
					}
					break;
				}
			}

			_data_msg_sizes[worker].push_back(msg_size);
		}

		if(_data_msg_sizes[worker].size() != num_msgs) {
			throw __FILE__ "Streamliner::CalculateMessageSizes - miscalculate message sizes";
		}
	}
}

int Streamliner::HandleWaitWriteNotify(int write_segment_id) {
	_num_ready_workers[write_segment_id]++;
	
	if(_num_ready_workers[write_segment_id] == _workers_with_output.size()) {
		return 1;
	}

	return 0;
}

void Streamliner::DumpOutputEntries(FILE* fp) {
	map <MPI_Offset, OutputEntry*>::iterator mit;

	fprintf(fp, "num_entries=%ld\n", _output_entries.size()); 
	for(mit=_output_entries.begin(); mit!=_output_entries.end(); mit++) {
		OutputEntry* oep = mit->second;

		fprintf(fp, "offset=%ld\n", oep->_out_off);
		fprintf(fp, "size=%ld\n", oep->_out_size);
		fwrite(oep->_out_data, 1, oep->_out_size, fp);
	}
	fflush(fp);
}

void Streamliner::DumpMessageSizes(FILE* fp) {
	map<int, vector<int> >::iterator mit;

	fprintf(fp, "Dumping data message sizes\n");

	for(mit=_data_msg_sizes.begin(); mit!=_data_msg_sizes.end(); mit++) {
		int worker = mit->first;
		vector<int>& msg_sizes = mit->second;

		fprintf(fp, "Worker %d:", worker);
		for(int i=0; i<msg_sizes.size(); i++) {
			fprintf(fp, " %d", msg_sizes[i]);
		}
		fprintf(fp, "\n");
	}
}

// provide a C wrapper function for result processing so it can be called from blast_hooks.c
void cProcessLocalOutputs(int query_id, int frag_id, ValNodePtr orplist, QueryOutputInfoPtr qoip, ByteStorePtr bs_ptr) {
	double track_time = MPI_Wtime();
	WriterWorker* writer_w;
	if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_SPLIT) {
		writer_w = WriterWorkerIndividual::Instance();
	} else if (output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		writer_w = WriterWorkerStreamline::Instance();
	} else if (output_strategy == WORKER_COLLECTIVE) {
		writer_w = WriterWorkerCollective::Instance();
	}

	writer_w->AddLocalOutputs(query_id, frag_id, orplist, qoip, bs_ptr);
	
	writer_w->ProcessResults(false, frag_id);
	process_output_time += MPI_Wtime() - track_time;
}

// C wrapper function for WriterMaster::UpdateQueryOutputInfo
void cUpdateQueryOutputInfo(int query_id, QueryOutputInfoPtr qoip) {
	WriterMaster* writer_m;
	if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_SPLIT) {
		writer_m = WriterMasterIndividual::Instance();
	} else if (output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		writer_m = WriterMasterStreamline::Instance();
	} else if (output_strategy == WORKER_COLLECTIVE) {
		writer_m = WriterMasterCollective::Instance();
	}

	writer_m->UpdateQueryOutputInfo(query_id, qoip);
}

// C wrapper function for Writer::UpdateDefaultReportInfo
void cUpdateDefaultReportInfo(QueryOutputInfoPtr qoip) {
	WriterMaster* writer_m;
	if(output_strategy == WORKER_INDIVIDUAL || output_strategy == WORKER_SPLIT) {
		writer_m = WriterMasterIndividual::Instance();
	} else if (output_strategy == MASTER_STREAMLINE || output_strategy == WORKER_STREAMLINE) {
		writer_m = WriterMasterStreamline::Instance();
	} else if (output_strategy == WORKER_COLLECTIVE) {
		writer_m = WriterMasterCollective::Instance();
	}

	writer_m->UpdateDefaultReportInfo(qoip);
}
