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

#ifndef __QUERY_MANAGER_HPP__
#define __QUERY_MANAGER_HPP__

#include <tofasta.h>

#ifdef __cplusplus
#include <string>
#include <map>
#include <vector>

class QueryM {
	public:
		static QueryM* Instance(const std::string& filename, bool use_vfile);
		static QueryM* Instance();
		SeqEntryPtr GetQueryEntry(int query_id);
		void AddQueryEntry(int query_id, SeqEntryPtr ptr);
		void AddQueryData(int query_id, char* buffer);
		void RemoveQuery(int query_id);
		int IndexQueries();
		void ParseQuerydata();
		void UnloadOldQueries(int query_id);
		void BcastIndexes(int brank);
		int LoadQueries(int start_query, int end_query);
		int GetQueryLen(int query_id) { return _query_lenV[query_id]; }
		int GetNumWorkingQueries() { return _query_entryM.size(); }
		void CloneQueryData(int query_id, char* buffer);
		int GetFirstQuery();
		void Destroy();

		void EnableAutoLoad() { _auto_load = true; }

	private:	
		static QueryM* _instance;
		std::string _filename;
		std::map <int, SeqEntryPtr> _query_entryM;
		bool _use_vfile; // use virtual files?
		int _query_count;
		std::vector <int> _start_posV; /** start positions of each query sequence */
		std::vector <int> _query_lenV; /** lengths of query sequences */
		std::map <int, char*> _query_dataM; /** query data */
		bool _auto_load;
		int _max_buf_size; // 
		
		QueryM(const std::string& filename, bool use_vfile);
		QueryM(const QueryM&);
		QueryM& operator=(QueryM &);
};

#endif

#ifdef __cplusplus
extern "C" {
#endif
void cQueryMAddQueryEntry(int query_id, SeqEntryPtr);
SeqEntryPtr cQueryMGetQueryEntry(int query_id);
#ifdef __cplusplus
}
#endif

#endif
