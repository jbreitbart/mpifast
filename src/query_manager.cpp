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

#include <query_manager.hpp>
#include <mpiblast_util.hpp>
#include <virtual_files.hpp>
#include <file_util.hpp>
using namespace std;

extern Int8** stats_array;

extern "C" {
extern SeqEntryPtr bufferToQueryEntry( char* buffer );
}

QueryM* QueryM::_instance=NULL;

QueryM* QueryM::Instance(const std::string& filename, bool use_vfile) {
	if(_instance == NULL) {
		_instance = new QueryM(filename, use_vfile);
	}

	return _instance;
}

QueryM* QueryM::Instance() {
	return _instance;
}

QueryM::QueryM(const std::string& filename, bool use_vfile) {
	_filename = filename;
	_use_vfile = use_vfile;
	_auto_load = false;
	_query_count = 0;

	if(max_query_load > 0) {
		_max_buf_size = max_query_load * 1024; 
	} else {
		_max_buf_size = 4 * 1024 * 1024; 
	}
}

void QueryM::AddQueryEntry(int query_id, SeqEntryPtr ptr) {
	if(_query_entryM.find(query_id) != _query_entryM.end()) {
		throw __FILE__ " QueryM::AddQuery -- duplicated query";
	}

	_query_entryM[query_id] = ptr;
	return;
}

void QueryM::AddQueryData(int query_id, char* buffer) {
	if(_query_dataM.find(query_id) != _query_dataM.end()) {
		throw __FILE__ " QueryM::AddQuery -- duplicated query";
	}

	if(stats_array[query_id] == NULL) {
		stats_array[query_id] = (Int8*)MemNew( sizeof(Int8)*12 );
	}

	_query_dataM[query_id] = buffer;
	return;
}

void QueryM::RemoveQuery(int query_id) {
	if(_query_entryM.find(query_id) != _query_entryM.end()) {
		SeqEntryFree(_query_entryM[query_id]);
		_query_entryM.erase(query_id);
	}

	if(_query_dataM.find(query_id) != _query_dataM.end()) {
		delete _query_dataM[query_id];
		_query_dataM.erase(query_id);
	}

	if(stats_array[query_id] != NULL) {
		MemFree(stats_array[query_id]);
		stats_array[query_id] = NULL;
	}

	return;
}

int QueryM::IndexQueries() {
	if(_use_vfile) {
		NlmMFILEPtr mfp=(NlmMFILEPtr)(VFM::Instance()->GetFilePtr(_filename));
		char* buffer = (char*)(mfp->mmp_begin);
		int len = VFM::Instance()->GetFileLength(_filename);

		for(int i=0; i<len; i++) {
			if(buffer[i] == '>') {
				_start_posV.push_back(i);
				int query_id = _query_count;
				_query_count++;

				if(query_id > 0) {
					_query_lenV.push_back(i - _start_posV[query_id - 1]);
				}
			}
		}
		_query_lenV.push_back(len - 1 - _start_posV[_query_count - 1]); // last position is NULLB

		return _query_count;

	} else {
		uint64 file_len = statFileSize(_filename.c_str());

		// read in all files for now, but really need to change to incremental reads
		int len = file_len;
		char* buffer = new char[len];

		ifstream ifs;
		ifs.open(_filename.c_str());
		if(!ifs.is_open()) {
			throw __FILE__ "QueryM::IndexQueries -- Cannot open query for read";
		}
		ifs.read(buffer, len);
		
		ifs.close();

		for(int i=0; i<len; i++) {
			if(buffer[i] == '>') {
				_start_posV.push_back(i);
				int query_id = _query_count;
				_query_count++;

				if(query_id > 0) {
					_query_lenV.push_back(i - _start_posV[query_id - 1]);
				}
			}
		}
		_query_lenV.push_back(len - _start_posV[_query_count - 1]); // last position is NULLB

		delete buffer;

		return _query_count;
	}

	return 0;
}

void QueryM::ParseQuerydata() {
	if((_query_count != _start_posV.size()) || (_query_count != _query_lenV.size())) {
		throw __FILE__ "BLAST::ParseQuerydata -- invalid query indexes.";
	}

	if(!query_in_mem) {
		throw __FILE__ "BLAST::ParseQuerydata -- query parsing has not been implemented for file I/O.";
	}

	NlmMFILEPtr mfp=(NlmMFILEPtr)(VFM::Instance()->GetFilePtr(_filename));
	char* buffer = (char*)(mfp->mmp_begin);

	for(int i=0; i<_query_count; i++) {
		char* data = new char[_query_lenV[i] + 1];
		memcpy(data, buffer + _start_posV[i], _query_lenV[i]);
		data[_query_lenV[i]] = NULLB;
		AddQueryData(i, data);
	}

	return;
}

void QueryM::UnloadOldQueries(int query_id) {
	if(query_id == 0) {
		return;
	}

	if(_query_dataM.empty() && _query_entryM.empty()) {
		return;
	}

	int min_query_id1 = _query_dataM.begin()->first;
	int min_query_id2 = _query_entryM.begin()->first;

	int min_query_id = min(min_query_id1, min_query_id2);

	for(int i=min_query_id; i<query_id; i++) {
		RemoveQuery(i);
	}

	return;
}

SeqEntryPtr QueryM::GetQueryEntry(int query_id) {
	bool found = false;

	if(_query_entryM.find(query_id) != _query_entryM.end()) {
		found = true;
	} else {
		bool loaded = false;
		if(_query_dataM.find(query_id) == _query_dataM.end()) {
			if(_auto_load && !query_in_mem) {
				LoadQueries(query_id, query_id + 1);
				loaded = true;
			}
		} else {
			loaded = true;
		}
		
		if(loaded) {
			SeqEntryPtr cur_sep = bufferToQueryEntry(_query_dataM[query_id]);
			if(cur_sep != NULL) {
				_query_entryM[query_id] = cur_sep;
				found = true;
			}
		}
	}

	if(found) {
		return _query_entryM[query_id];
	} else {
		throw __FILE__ "QueryM::GetQueryEntry -- cannot find query entry.";
		return NULL;
	}
}

void cQueryMAddQueryEntry(int query_id, SeqEntryPtr ptr) {
	QueryM::Instance()->AddQueryEntry(query_id, ptr);
}

int QueryM::LoadQueries(int start_query, int end_query) {
	// decide which query has not been loaded
	int load_start = -1;
	for(int i=start_query; i<end_query; i++) {
		if(_query_dataM.find(i) == _query_dataM.end()) {
			load_start = i;
			break;
		}
	}

	if(load_start == -1) {
		return 0;
	}
	
	double read_start = MPI_Wtime();
	// fetch a chunk of data from the query file, parse and load it
	int read_len = 0;
	int load_stop = -1;
	for(int i=load_start; i<_query_count; i++) {
		if((read_len >= _max_buf_size) || (_query_dataM.find(i) != _query_dataM.end())) {
			break;
		}
		read_len += _query_lenV[i];		
		load_stop = i + 1;
	}

	if(debug_msg) {
		LOG_MSG << "Loading queries from " << load_start << " to " << load_stop << " read_len " << read_len << endl;
	}
	
	if(read_len > 0) {
		char* buffer = new char[read_len];
		ifstream ifs(_filename.c_str());
		ifs.seekg(_start_posV[load_start], ios::beg);
		ifs.read(buffer, read_len);
		ifs.close();

		char* ptr = buffer;
		// parsing queries
		for(int i=load_start; i<load_stop; i++) {
			char * data = new char[_query_lenV[i] + 1];
			memcpy(data, ptr, _query_lenV[i]);
			data[_query_lenV[i]] = NULLB;
			AddQueryData(i, data);

			ptr += _query_lenV[i];
		}

		delete buffer;
	}
	load_queries_time += MPI_Wtime() - read_start;

	return 0;
}

void QueryM::Destroy() {
	map <int, SeqEntryPtr>::iterator mit;

	for(mit=_query_entryM.begin(); mit!=_query_entryM.end(); mit++) {
		SeqEntryFree(mit->second);
	}
	_query_entryM.clear();

	map <int, char*>::iterator mit_data;
	for(mit_data=_query_dataM.begin(); mit_data!=_query_dataM.end(); mit_data++) {
		delete mit_data->second;
	}
	_query_dataM.clear();
}

void QueryM::BcastIndexes(int brank) {
	if(debug_msg) {
		LOG_MSG << "Begin broadcasting query indexes" << endl;
	}
	
	MPI_Bcast(&_query_count, 1, MPI_INT, brank, MPI_COMM_WORLD);

	int* buf = new int[_query_count];
	if(rank == brank) {
		for(int i=0; i<_query_count; i++) {
			buf[i]=_start_posV[i];
		}
	}
	
	MPI_Bcast(buf, _query_count, MPI_INT, brank, MPI_COMM_WORLD);

	if(rank != brank) {
		for(int i=0; i<_query_count; i++) {
			_start_posV.push_back(buf[i]);
		}
	}

	if(rank == brank) {
		for(int i=0; i<_query_count; i++) {
			buf[i]=_query_lenV[i];
		}
	}
	
	MPI_Bcast(buf, _query_count, MPI_INT, brank, MPI_COMM_WORLD);

	if(rank != brank) {
		for(int i=0; i<_query_count; i++) {
			_query_lenV.push_back(buf[i]);
		}
	}
	
	delete buf;
	if(debug_msg) {
		LOG_MSG << "After broadcasting query indexes" << endl;
	}
}

void QueryM::CloneQueryData(int query_id, char* buffer) {
	if(_query_dataM.find(query_id) == _query_dataM.end()) {
		throw __FILE__ " QueryM::CloneQueryData -- query data has not been loaded";	
	}
	memcpy(buffer, _query_dataM[query_id], _query_lenV[query_id]);
}

int QueryM::GetFirstQuery() {
	if(_query_dataM.empty()) {
		GetQueryEntry(0);
	}
	int query_id = _query_dataM.begin()->first;
	return query_id;
}

SeqEntryPtr cQueryMGetQueryEntry(int query_id) {
	return QueryM::Instance()->GetQueryEntry(query_id);
}
