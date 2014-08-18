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
#ifndef __mpiblastall_h__
#define __mpiblastall_h__

#ifdef __cplusplus
extern "C"{
#endif

//#include <stdlib.h>

#ifndef Int2
#define Int2 short
#endif

#ifndef Int4
#define Int4 int
#endif

#ifndef Nlm_FloatHi
#define	Nlm_FloatHi double
#endif

#ifndef Boolean
#define	Boolean  unsigned char
#endif

#ifndef Uint1
#define	Uint1  unsigned char
#endif

#ifdef USE_NCBI_ASCII
const static int use_binary_asn = FALSE;
#else
const static int use_binary_asn = TRUE;
#endif

#define MASTER_MODE		0
#define WORKER_MODE		1

	/** Reads command line and initializes BLAST data structures */
Int2 initBLAST (int mpi_mode);
/** Outputs an HTML header IF the HTML command line option was set */
void outputHtmlHeader();
void outputHtmlHeaderPIO(char * buffer);
/** Outputs the results in <code>sappy</code> */
Int2 outputResults( SeqAlignPtr sappy, int queryI );
Int2 outputResultsPIO( int mode, SeqAlignPtr sappy, int queryI, int frag_id, ByteStorePtr bs_ptr );
/** Outputs an HTML footer IF the HTML command line option was set */
void outputHTMLfooter();
void outputHtmlFooterPIO(char *buffer);
/** Cleans up BLAST data structures */
void cleanupBLAST();


/** Load query sequences from the multi-FastA file */
int loadQueries( void );
int loadQueriesEx( void );
int loadQueriesFromBuffer( char* buffer );
int loadQueriesFromBufferEx( char* buffer );
int countQueriesInBuffer( char* buffer, int buf_size );
int updateBlastDB( const char* db_name, int len );
/** Free memory used to store the query sequences */
void cleanupQueries( void );

void MpiBlastEnableBioseqFetch();
void MpiBlastDisableBioseqFetch();

SeqEntryPtr bufferToQueryEntry(char* buffer);

/** Extrace stats data sent from worker */
int extractStats(int queryI, unsigned char* stats_buf, int stats_buf_size);

/** 
 * runBLAST executes BLAST in one of two modes:
 * SEARCH_MODE Executes BLAST and stores successful result ids in result_id_list 
 * COLLECT_STATS_MODE Prepares all data structures and performs calculations related
 *                    to BLAST without actually performing the search
 * When first_query == last_query only a single query is processed
 * to process all queries, set first_query = 0 and last_query = query_count - 1
 * @param mode	The mode to run BLAST in, as described above
 * @param first_query	The first query to process
 * @param last_query	The last query to process
 */
Int2 runBLAST( int mode, int first_query, int last_query );
Int2 runBLASTPIO( int mode, int first_query, int last_query, int frag_id );
extern BLAST_OptionsBlkPtr options;
extern ValNodePtr result_id_list;
extern ValNodePtr result_bsp_list;	/**< List of ByteStorePtrs, each of which contains all results (SeqAligns+Bioseqs) for a query */
extern int debug;			/**< print debugging messages when set to true */
extern int resume_run;  /**<choose append to outputfile when true */ 
extern Int8* query_adj_array;  /**< Array of query length adjustments (one per query) */
extern Int8* db_len_array;     /**< Array of effective database lengths (one per query) */
extern Int8** stats_array;		/**< Array of per-query search statistics */
extern Boolean strict_output_conformance;	/**< Set to true if the Altschul paper reference must be given instead of the mpiBLAST paper reference */
const static int SEARCH_MODE = 0;
const static int COLLECT_STATS_MODE = 1;

const static int OUTPUT_MODE = 0;
const static int COLLECT_INFO_MODE = 1;

extern Int4 global_dbseq_num;
extern Int8 global_db_len;

/** stub for Bioseq* functions.  Can't directly include edutil.h
  * for some reason.  */
NLM_EXTERN Boolean LIBCALL BioseqDelete( SeqIdPtr target, Int4 from, Int4 to, 
					   Boolean do_feat, Boolean do_split );

NLM_EXTERN Boolean LIBCALL BioseqInsert( SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, 
					   Boolean from_feat, Boolean to_feat, Boolean do_split );

NLM_EXTERN BioseqPtr LIBCALL BioseqCopy (SeqIdPtr newid, SeqIdPtr sourceid, Int4 from, Int4 to,
                               Uint1 strand, Boolean do_feat);

int GetNumOfShowDes();
int GetNumOfShowAlign();
Uint1 GetFormatAlignView();
Boolean GetFormatHtml();
int GetQueryAcknowledge(int query_id, char** query_ack);

static AsnType biginttype = {BIGINT_TYPE, "BIGINT" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL};
static AsnType set_type = {313, "SET" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL}; 
static AsnType len_type = {0, "len" ,128,6,0,0,0,0,0,0,NULL,&biginttype,NULL,0,NULL};
static AsnType lenset_type = {402, "Subject-len" ,1,0,0,0,0,1,0,0,NULL,&set_type,&len_type,0,NULL};
static AsnType statdata_type = {0, "stat-data" ,128,6,0,0,0,0,0,0,NULL,&biginttype,NULL,0,NULL};
static AsnType statset_type = {402, "Stat-data" ,1,0,0,0,0,1,0,0,NULL,&set_type,&statdata_type,0,NULL};

#ifdef __cplusplus
}
#endif

#endif /* __mpiblastall_h__ */
