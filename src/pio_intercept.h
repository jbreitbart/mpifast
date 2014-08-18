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
#ifndef __pio_intercept_h__
#define __pio_intercept_h__

#include <mpiblast_types.h>

#ifdef __cplusplus
#include <mpi.h>
extern "C"{
#endif

//#include <stdio.h>
#include <txalign.h>

/* marcros are duplicated in ncbi files to avoid header file conflict */
#define CHECK_NULPTR(x) if(x==NULL) fprintf(stderr, "%s, line%d: cannot allocate enough memory\n", __FILE__, __LINE__);

#ifndef NO_MPI
typedef struct _output_record {
	int rank;     /* 1) rank of the worker who submit it */
	int frag_id;  /* 2) fragment id */
	SeqIdPtr id;  /* 3) not used now */
	Nlm_FloatHi bit_score; /* 4) */ 
	Nlm_FloatHi evalue;    /* 5) */
	int des_size;   /* 6) size of the description output data */
	char* des_data; /* 7) description output data */
	int aln_size;   /* 8) size of the alignment output data */
	char* aln_data; /* 9) segment output data */
	MPI_Offset des_offset; /* 10) offset of the description data */
	MPI_Offset aln_offset; /* 11) offset of the segment data */
} OutputRecord, *OutputRecordPtr;
#endif

// Macros for send type, use bits to indicate what field in the OutputRecord should be sent
#define RANK 		1<<0
#define FRAG_ID 	1<<1
#define ID			1<<2
#define BIT_SCORE	1<<3
#define EVALUE		1<<4
#define DES_SIZE	1<<5
#define DES_DATA	1<<6
#define ALN_SIZE	1<<7
#define ALN_DATA	1<<8
#define DES_OFFSET	1<<9
#define ALN_OFFSET	1<<10
#define TXT_ID		1<<11

#define SEND_TYPE1	0|RANK|FRAG_ID|BIT_SCORE|EVALUE|DES_SIZE|ALN_SIZE
#define SEND_TYPE2	0|RANK|FRAG_ID|BIT_SCORE|EVALUE|ALN_SIZE
#define SEND_TYPE3	0|RANK|FRAG_ID|BIT_SCORE|EVALUE|DES_SIZE|DES_DATA|ALN_SIZE|ALN_DATA
#define SEND_TYPE4	0|RANK|FRAG_ID|BIT_SCORE|EVALUE|ALN_SIZE|ALN_DATA

typedef struct _query_output_info {
	char* header;	// header info, collected in master
	int header_size;
	char* des_hdr;	// description header, collected in worker
	int des_hdr_size;
	char* des_ftr;	// description footer, collected in worker
	int des_ftr_size;
	char* no_hits;	// no hits prompt, collected in master
	int no_hits_size;
	char* stat;		// statistics output, collected in master
	int stat_size;
	char* footer;	// footer, collected in master
	int footer_size;
} QueryOutputInfo, *QueryOutputInfoPtr;

typedef struct _default_report {
	char* ref; // reference
	int ref_size;
	char* db_info; // database info
	int db_info_size;
	char* no_hits;	// no hits prompt
	int no_hits_size;
	char* db_report; // database report
	int db_report_size;
	char* footer; 
	int footer_size;
} DefaultReport, *DefaultReportPtr;

#define INT_FPRINTF		1<<0	// intercept fprintf call
#define INT_FWRITE		1<<1	// intercept fwrite call

#define PRINT_BUFFER_LENGTH  2048

#define INT_INFO				1<<0	// only intercept output info
#define INT_RECORDS				1<<1	// intercept outut records
#define INT_INFO_RECORDS		1<<2	// intercept outut info and records

extern int hijack_default;

int intercept_output_asn(int tag);
int intercept_output_pairwise(int tag);
int intercept_output_tabular(int tag);

//FILE *fopen(const char *path, const char *mode);
QueryOutputInfoPtr new_query_output_infop(void);
int free_query_output_infop(QueryOutputInfoPtr qoip);

DefaultReportPtr new_default_reportp(void);
int free_default_reportp(DefaultReportPtr drp);

int init_intercept_lib();
int fin_intercept_lib();

int create_output_records(SeqAlignPtr seqalign, ValNodePtr PNTR orplist, int maxnum);
//int start_intercept(FILE *stream, int format, Boolean HTML, int mode, int function); 
// Assume NCBI BLAST either output with fprintf or fwrite, if both, add function parameter
int start_intercept(FILE *stream, int format, Boolean HTML, int mode);  
int set_results_info(ValNodePtr orplist, int effective_len);
int end_intercept(int query_id, int frag_id, ByteStorePtr bs_ptr);

int start_comm_intercept(FILE* stream);
int end_comm_intercept(char** buffer);

int parse_output_xml();

#ifndef NO_MPI
void print_output_records(FILE* fp, ValNodePtr orplist);
void print_output_record(FILE* fp, OutputRecordPtr orp);
int OutputRecordCmp(OutputRecordPtr lorp, OutputRecordPtr rorp);
#endif

#ifdef __cplusplus
}
#endif
#endif /* __pio_intercept_h__ */
