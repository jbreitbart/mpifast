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
#ifndef __BLAST_WRAPPER_HPP__
#define __BLAST_WRAPPER_HPP__
#include <vector>
#include <string>
#include <stdlib.h>
#include <ncbi.h>
#include <objseq.h>
#include "pio_intercept.h"

using namespace std;

class BLAST {
public:
	BLAST() {
		_query_count=-1;
		_all_dates = NULL;
		_total_datelen = 0;
		_max_show = 500;
	}

	~BLAST() {
		if(_all_dates!=NULL) {
			free(_all_dates);
		}
	}

	char* GetAllDates(void);
	int GetTotalDateLen();
	int GetNumDes();
	int GetNumAlign();
	int GetMaxShow();
	int GetAlignView();
	bool GetHtml();
	void GetHtmlHeader(char* buffer);
	void GetHtmlFooter(char* buffer);
	void SetSharedDBName(string& adb_name);
	void SetDBType(string& adb_type);
	void SetQueryCount(int aquery_count);
	void SetQueryFile(const string& filename);
	
	void InitNCBI(vector< string > &opts);
	void Init(int mpi_mode);
	short Run(int mode, int first_query, int last_query);
	short RunPIO(int mode, int first_query, int last_query, int frag_id);
	void GetQueryOutputInfo(int query_id);
	void GetDefaultReportInfo();
	void SynthesizeReport(int query_id, DefaultReportPtr drp, QueryOutputInfoPtr qoip);
	void Cleanup(void);
	int LoadQueries(void);
	int LoadQueriesFromBuffer( char* buffer );
	
	void AllocSearchStats(void);
	void MasterDistributeEvalAdjusts();
	void WorkerRecvEvalAdjusts(void);
	void MasterDistributeFragTimes();
	void WorkerRecvFragTimes();
	void ExtractStats(int queryI, unsigned char* stats_buf, int stats_buf_size);
	
	// for debug
	void PrintSearchStatics(void);
private:
	string _shared_db_name;
	string _db_type;

	int _total_datelen;	/**< size of storage for fragment datestamps */
	char* _all_dates;	/**< storage for fragment datestamps */
	int _query_count;	/**< The number of queries in the query file */
	int _num_desc;		/** -v option, number of descriptions to be showed*/
	int _num_align;		/** -b option, number of alignments to be showed*/
	int _max_show;		/** MAX(_num_desc, _num_align) */
	int _align_view;	/** output format */
	bool _html;			/** using html */
	string _query_filename; /** name of the query file */
};

inline char* BLAST::GetAllDates(void) {
	return _all_dates;
}

inline int BLAST::GetTotalDateLen() {
	return _total_datelen;
}

inline void BLAST::SetSharedDBName(string& adb_name) {
	_shared_db_name = adb_name;
}

inline void BLAST::SetDBType(string& adb_type) {
	_db_type = adb_type;
}

inline void BLAST::SetQueryCount(int aquery_count) {
	_query_count = aquery_count;
}

inline void BLAST::SetQueryFile(const string& filename) {
	_query_filename = filename;
}

inline int BLAST::GetNumDes() {
	return _num_desc;
}

inline int BLAST::GetNumAlign() {
	return _num_align;
}

inline int BLAST::GetMaxShow() {
	return _max_show;
}

inline int BLAST::GetAlignView() {
	return _align_view;
}

inline bool BLAST::GetHtml() {
	return _html;
}

#endif
