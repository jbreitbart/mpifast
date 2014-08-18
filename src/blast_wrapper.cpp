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
//#include <ncbi.h>
//#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <xmlblast.h>
#include <mblast.h>
#include <blfmtutl.h>
#include <blast_hooks.h>
#include <mpiblast_util.hpp>
#include <blast_wrapper.hpp>
#include <pio_intercept.h>
#include <virtual_files.hpp>
#include <mpi.h>
#include <query_manager.hpp>

void BLAST::InitNCBI(vector< string > &opts) {
	initNCBI(opts);
}

void BLAST::Init(int mpi_mode) {
	int rval = initBLAST(mpi_mode);

	if( rval != 0 ){
		// bail out on error for now.
		LOG_MSG << "Error: got initBLAST exit code " << rval << endl;
		throw __FILE__ "(BLAST::Init): Received non-zero initBLAST exit code";
	}

	_num_desc=GetNumOfShowDes();
	_num_align=GetNumOfShowAlign();

	_align_view=GetFormatAlignView();

	if(_align_view == 8 || _align_view == 9) {
		_max_show=_num_align;
	} else {
		_max_show=MAX(_num_desc, _num_align);
	}

	if(GetFormatHtml()==TRUE) {
		_html = true;
	} else {
		_html = false;
	}
}

short BLAST::Run( int mode, int first_query, int last_query ) {
	return runBLAST(mode, first_query, last_query);
}

short BLAST::RunPIO( int mode, int first_query, int last_query, int frag_id ) {
	return runBLASTPIO(mode, first_query, last_query, frag_id);
}

// do a dry run of outputResultsPIO, get query output info and update related QueryOutput
void BLAST::GetQueryOutputInfo(int query_id) {
	outputResultsPIO(COLLECT_INFO_MODE, NULL, query_id, -1, NULL); 
}

void BLAST::SynthesizeReport(int query_id, DefaultReportPtr drp, QueryOutputInfoPtr qoip) {
	int query_ack_size;
	char* query_ack = NULL;
	
	if(_align_view < 7) {
		query_ack_size = GetQueryAcknowledge(query_id, &query_ack);
		if(query_ack_size == 0) {
			throw __FILE__ "BLAST::SynthesizeReport - fail to get query acknowledge";
		}

		assert(qoip->header_size == 0);
		qoip->header_size = drp->ref_size + query_ack_size + drp->db_info_size; 
		qoip->header = (char*)malloc(qoip->header_size);
		CHECK_NULPTR(qoip->header);
		char* ptr = qoip->header;
		memcpy(ptr, drp->ref, drp->ref_size);
		ptr += drp->ref_size;
		memcpy(ptr, query_ack, query_ack_size);
		ptr += query_ack_size;
		memcpy(ptr, drp->db_info, drp->db_info_size);

		if(drp->no_hits_size > 0) {
			qoip->no_hits_size = drp->no_hits_size;
			qoip->no_hits = (char*)malloc(drp->no_hits_size);
			CHECK_NULPTR(qoip->no_hits);
			memcpy(qoip->no_hits, drp->no_hits, drp->no_hits_size);
		}

		if(drp->db_report_size > 0) {
			qoip->stat_size = drp->db_report_size;
			qoip->stat = (char*)malloc(drp->db_report_size);
			CHECK_NULPTR(qoip->stat);
			memcpy(qoip->stat, drp->db_report, drp->db_report_size);
		}

		if(query_ack != NULL) {
			free(query_ack);
		}
	}
}

void BLAST::GetDefaultReportInfo() {
	if(_align_view > 0) {
		return;
	}

	hijack_default = 1;
	int query_id = QueryM::Instance()->GetFirstQuery();
	if(debug_msg) {
		LOG_MSG << "First query is " << query_id << endl;
	}
	outputResultsPIO(COLLECT_INFO_MODE, NULL, query_id, -1, NULL); 
	hijack_default = 0;
}

int BLAST::LoadQueries( void ) {
	if(query_in_mem) {
		NlmMFILEPtr mfp=(NlmMFILEPtr)(VFM::Instance()->GetFilePtr(_query_filename));
		char* buffer = (char*)(mfp->mmp_begin);
		
		if(use_query_map) {
			return loadQueriesFromBufferEx(buffer);
		} else {
			return loadQueriesFromBuffer(buffer);
		}
	} else {
		if(use_query_map) {
			return loadQueriesEx();
		} else {
			return loadQueries();
		}
	}
}

int BLAST::LoadQueriesFromBuffer( char* buffer ) {
	if(use_query_map) {
		return loadQueriesFromBufferEx(buffer);
	} else {
		return loadQueriesFromBuffer(buffer);
	}
}

void BLAST::GetHtmlHeader(char* buffer) {
	outputHtmlHeaderPIO(buffer);
}

void BLAST::GetHtmlFooter(char* buffer) {
	outputHtmlFooterPIO(buffer);
}

void BLAST::AllocSearchStats( void ) {
	// allocate space for search statistics
	stats_array = (Int8**)MemNew( sizeof(Int8*)*_query_count );
//	stats_array[0] = (Int8*)MemNew( sizeof(Int8)*12*_query_count );
//	memset( stats_array[0], 0, sizeof(Int8)*12*_query_count );
//	for( int qI = 0; qI < _query_count; qI++ )
//		stats_array[qI] = &stats_array[0][qI*12];

//	for( int qI = 0; qI < _query_count; qI++ )
//		stats_array[qI] = NULL;
}

void BLAST::MasterDistributeEvalAdjusts() {
	AllocSearchStats();

	if(parallel_evalue_adjust) {
		if(!use_real_dblen) {
			Run( COLLECT_STATS_MODE, 0, 0 );
		}

		MPI_Bcast( &global_db_len, 8, MPI_BYTE, super_master_process, MPI_COMM_WORLD); 
		MPI_Bcast( &global_dbseq_num, 4, MPI_BYTE, super_master_process, MPI_COMM_WORLD); 

		if(debug_msg) {
			LOG_MSG << "global_db_len=" << global_db_len << ", global_dbseq_num=" << global_dbseq_num << endl;
		}
		
	} else {
		// allocate the eff. query and db length arrays
		options->query_adjustments = (Int8*)MemNew( sizeof( Int8 ) * _query_count );
		if(options->query_adjustments == NULL) {
			throw __FILE__ "BLAST::AllocMasterEvalAdjusts: cannot allocate enough memory";
		}
		options->effective_db_lengths = (Int8*)MemNew( sizeof( Int8 ) * _query_count );
		if(options->effective_db_lengths == NULL) {
			throw __FILE__ "BLAST::AllocMasterEvalAdjusts: cannot allocate enough memory";
		}
		
		if(fast_evalue_approximation) {
			Run( COLLECT_STATS_MODE, 0, 0 );
			for( int qI = 1; qI < _query_count; qI++ )
			{
				options->query_adjustments[qI] = options->query_adjustments[0];
				options->effective_db_lengths[qI] = options->effective_db_lengths[0];
			}
		} else {
			Run( COLLECT_STATS_MODE, 0, _query_count - 1 );
		}
		MPI_Bcast( options->query_adjustments, _query_count, MPI_LONG_LONG_INT, super_master_process, MPI_COMM_WORLD );
		MPI_Bcast( options->effective_db_lengths, _query_count, MPI_LONG_LONG_INT, super_master_process, MPI_COMM_WORLD );
	}
}

void BLAST::WorkerRecvEvalAdjusts(void) {
	AllocSearchStats();

	if(parallel_evalue_adjust) {
		MPI_Bcast( &global_db_len, 8, MPI_BYTE, super_master_process, MPI_COMM_WORLD); 
		MPI_Bcast( &global_dbseq_num, 4, MPI_BYTE, super_master_process, MPI_COMM_WORLD); 
	} else {
		// allocate space and receive the array
		query_adj_array = (Int8*)MemNew( sizeof( Int8 ) * _query_count );
		if(query_adj_array == NULL) {
			throw __FILE__ "BLAST::AllocWorkerEvalAdjusts: cannot allocate enough memory";
		}
		db_len_array = (Int8*)MemNew( sizeof( Int8 ) * _query_count );
		if(db_len_array == NULL) {
			throw __FILE__ "BLAST::AllocWorkerEvalAdjusts: cannot allocate enough memory";
		}

		MPI_Bcast( query_adj_array, _query_count, MPI_LONG_LONG_INT, super_master_process, MPI_COMM_WORLD );
		MPI_Bcast( db_len_array, _query_count, MPI_LONG_LONG_INT, super_master_process, MPI_COMM_WORLD );
	}	
}

void BLAST::MasterDistributeFragTimes( void ) {
	// open the database info files
	char* m_dbname = (char*)malloc( _shared_db_name.size() + 1 );
	CHECK_NULPTR(m_dbname);
	strcpy( m_dbname, _shared_db_name.c_str() );
	ReadDBFILEPtr rdfp = readdb_new_ex2( m_dbname, _db_type == "p", READDB_NEW_DO_REPORT, NULL, NULL );

	free( m_dbname );
	ReadDBFILEPtr cur_rdfp = rdfp;
	_total_datelen = 0;
	while( cur_rdfp != NULL ){
		char* date_str = readdb_get_date( cur_rdfp );
		_total_datelen += strlen( date_str ) + 1;
		cur_rdfp = cur_rdfp->next;
	}
	_all_dates = (char*)malloc( _total_datelen );
	CHECK_NULPTR(_all_dates);
	char* cur_date = _all_dates;
	cur_rdfp = rdfp;
	while( cur_rdfp != NULL ){
		char* date_str = readdb_get_date( cur_rdfp );
		strcpy( cur_date, date_str );
		cur_date += strlen( date_str ) + 1;
		cur_rdfp = cur_rdfp->next;
	}
	MPI_Bcast( &_total_datelen, 1, MPI_INT, super_master_process, MPI_COMM_WORLD );
	MPI_Bcast( _all_dates, _total_datelen, MPI_CHAR, super_master_process, MPI_COMM_WORLD );
	rdfp = readdb_destruct( rdfp );
	if( debug_msg )	LOG_MSG << "First date: " << _all_dates << endl;
}

void BLAST::WorkerRecvFragTimes( void ) {
	MPI_Bcast( &_total_datelen, 1, MPI_INT, super_master_process, MPI_COMM_WORLD );
	_all_dates = (char*)malloc( _total_datelen );
	CHECK_NULPTR(_all_dates);
	MPI_Bcast( _all_dates, _total_datelen, MPI_CHAR, super_master_process, MPI_COMM_WORLD );
	if( debug_msg )	LOG_MSG << "First date: " << _all_dates << endl;
}

void BLAST::Cleanup(void) {
	cleanupBLAST();
}

void BLAST::PrintSearchStatics(void) {
	for(int i=0; i<_query_count; i++) {
		cerr<<"query "<<i<<" query_adjustment= "<<options->query_adjustments[i]<<" , "
			<<" query_effective_db_lengths= "<<options->effective_db_lengths[i]<<endl;
	}
}

void BLAST::ExtractStats(int queryI, unsigned char* stats_buf, int stats_buf_size) {
	int ret = extractStats(queryI, stats_buf, stats_buf_size);
	if(ret != 0) {
		throw __FILE__ "BLAST::ExtractStats() errors in extracing statistics";
	}
}
