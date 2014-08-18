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
/**************************************************************************
*																		  *
*							  COPYRIGHT	NOTICE							  *
*																		  *
* This software/database is	categorized	as "United States Government	  *
* Work"	under the terms	of the United States Copyright Act.	 It	was		  *
* produced as part of the author's official	duties as a	Government		  *
* employee and thus	can	not	be copyrighted.	 This software/database	is	  *
* freely available to the public for use without a copyright notice.	  *
* Restrictions can not be placed on	its	present	or future use.			  *
*																		  *
* Although all reasonable efforts have been	taken to ensure	the	accuracy  *
* and reliability of the software and data,	the	National Library of		  *
* Medicine (NLM) and the U.S. Government do	not	and	can	not	warrant	the	  *
* performance or results that may be obtained by using this	software,	  *
* data,	or derivative works	thereof.  The NLM and the U.S. Government	  *
* disclaim any and all warranties, expressed or	implied, as	to the		  *
* performance, merchantability or fitness for any particular purpose or	  *
* use.																	  *
*																		  *
* In any work or product derived from this material, proper	attribution	  *
* of the author(s) as the source of	the	software or	data would be		  *
* appreciated.															  *
*																		  *
************************************************************************** 
*/
// #include <mpi.h>

#include <ncbi.h>
#include <objseq.h>
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
#include <algo/blast/composition_adjustment/composition_constants.h>
#ifdef BLAST_CS_API
#include <objblst3.h>
#include <netblap3.h>
#endif

#include "blast_hooks.h"
#include "mpi_util.h"
#include "intercept_def.h"
#include "query_manager.hpp"

#define NO_MPI		1
#include "pio_intercept.h"

/* Amount to relax the evalue threshold for preliminary alignments
 *  * when compositionally adjusted score matrices are used. */
#define EVALUE_EXPAND 1

extern int pio_hijacking;
extern int use_brief_report;
extern double result_time, search_time, process_output_time, curr_process_time;
extern int query_in_mem;
extern int use_query_map;
extern FILE* dbgfp;

/* Local function definition */
void writePartialSubjectBioseq( SeqAlignPtr sap_p, AsnIoBSPtr aibp );

SeqEntryPtr	getQuery( int queryI );
SeqEntryPtr	getQueryEx( int queryI );
extern SeqEntryPtr cQueryMGetQueryEntry(int query_id);

/* Used	by the callback	function. */
FILE *global_fp=NULL;
/*
Callback to print out ticks, in UNIX only due to file systems
portability issues.
*/

tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{
#ifdef OS_UNIX
    /* #ifndef BLAST_CS_API */
    fprintf(global_fp, "%s", ".");
    fflush(global_fp);
    /* #endif */
#endif
    return 0;
}

/* FIXME: better name, move to API directory?? */
static Boolean
sGetLoc(char* location, Int4* start, Int4* end)
{
        CharPtr delimiters = " ,;";
       
        if (start == NULL || end == NULL)
           return FALSE;
      
        *start = 0;
        *end = 0;
      
        if (location == NULL)
           return TRUE;
                       
        *start =  atoi(StringTokMT(location, delimiters, &location));
        *end = atoi(location);
                   
        return TRUE;  
}                  

static Int2
BlastGetMaskingLoc(FILE *infp, FILE *outfp, CharPtr instructions)
{
	BioseqPtr bsp;
	Char buffer[50];
	SeqEntryPtr sep;
	SeqLocPtr slp, slp_start, tmp_slp;

	if (infp == NULL || outfp == NULL || instructions == NULL)
		return 1;

	while ((sep=FastaToSeqEntryEx(infp, TRUE, NULL, TRUE)) != NULL) 
	{
		bsp = NULL;
		SeqEntryExplore(sep, &bsp, FindNuc);

		if (bsp == NULL)
		{
			ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
			return 2;
		}
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, 50);
		fprintf(outfp, ">%s\n", buffer);
		slp_start = slp = BlastBioseqFilter(bsp, instructions);
		while (slp)
		{
			tmp_slp=NULL;
			while((tmp_slp = SeqLocFindNext(slp, tmp_slp))!=NULL)
			{
				fprintf(outfp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp)), (long) (1+SeqLocStop(tmp_slp)));
			}
			slp = slp->next;
		}

		slp_start = SeqLocSetFree(slp_start);
		sep = SeqEntryFree(sep);
	}

	return 0;
}

typedef enum {
ARG_PROGRAM = 0,
ARG_DB,
ARG_QUERY,
ARG_EVALUE,
ARG_FORMAT,
ARG_OUT,
ARG_FILTER,
ARG_GAPOPEN,
ARG_GAPEXT,
ARG_XDROP,
ARG_SHOWGIS,
ARG_MISMATCH,
ARG_MATCH,
ARG_DESCRIPTIONS,
ARG_ALIGNMENTS,
ARG_THRESHOLD,
ARG_GAPPED,
ARG_QGENETIC_CODE,
ARG_DBGENCODE,
ARG_THREADS, 
ARG_ASNOUT,
ARG_BELIEVEQUERY,
ARG_MATRIX,
ARG_WORDSIZE,
ARG_DBSIZE,
ARG_BESTHITS,
ARG_MULTIPLEHITS,
ARG_SEARCHSP,
ARG_STRAND,
ARG_HTML,
#ifdef BLAST_CS_API
ARG_ENTREZQ,
#else
ARG_GILIST,
#endif
ARG_LCASE,
ARG_XDROP_UNGAPPED,
ARG_XDROP_FINAL,
#ifdef BLAST_CS_API
ARG_RPSBLAST,
#else
ARG_PSITCHKPNT,
#endif
ARG_USEMEGABLAST,
ARG_QUERYLOC,
ARG_WINDOW,
ARG_FRAMESHIFT,
ARG_INTRON,
#ifndef BLAST_CS_API
ARG_NUMQUERIES,
#ifndef BLASTALL_TOOLS_ONLY
ARG_FORCE_OLD,
#endif
#endif
ARG_COMP_BASED_STATS,
ARG_SMITH_WATERMAN,
#ifdef ALLOW_FULL_SMITH_WATERMAN
ARG_SMITH_WATERMAN_ALL
#endif
} BlastArguments;

#define	NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs[] = {
    { "Program Name",           
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},    /* ARG_PROGRAM */
    { "Database",               
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},    /* ARG_DB */
    { "Query File",            
      "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL}, /* ARG_QUERY */
    { "Expectation value (E)",  
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},    /* ARG_EVALUE */
    { "alignment view options:\n0 = pairwise,\n1 = query-anchored showing identities,\n2 = query-anchored no identities,\n3 = flat query-anchored, show identities,\n4 = flat query-anchored, no identities,\n5 = query-anchored no identities and blunt ends,\n6 = flat query-anchored, no identities and blunt ends,\n7 = XML Blast output,\n8 = tabular, \n9 tabular with comment lines\n10 ASN, text\n11 ASN, binary", /* 4 */
      "0", "0", "11", FALSE, 'm', ARG_INT, 0.0, 0, NULL},         /* ARG_FORMAT */
    { "BLAST report Output File", 
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL}, /* ARG_OUT */
    { "Filter query sequence (DUST with blastn, SEG with others)", 
      "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},       /* ARG_FILTER */
    { "Cost to open a gap (-1 invokes default behavior)", 
      "-1", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},          /* ARG_GAPOPEN */
    { "Cost to extend a gap (-1 invokes default behavior)", 
      "-1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},          /* ARG_GAPEXT */
    { "X dropoff value for gapped alignment (in bits) (zero invokes default "
      "behavior)\n      blastn 30, megablast 20, tblastx 0, all others 15", 
      "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},          /* ARG_XDROP */
    { "Show GI's in deflines",  /* 10 */
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},      /* ARG_SHOWGIS */
    { "Penalty for a nucleotide mismatch (blastn only)", 
      "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},         /* ARG_MISMATCH */
    { "Reward for a nucleotide match (blastn only)", 
      "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},          /* ARG_MATCH */
    { "Number of database sequences to show one-line descriptions for (V)", 
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},         /*  ARG_DESCRIPTIONS */
    { "Number of database sequence to show alignments for (B)", 
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},        /* ARG_ALIGNMENTS */
    { "Threshold for extending hits, default if zero\n" 
      "      blastp 11, blastn 0, blastx 12, tblastn 13\n"
      "      tblastx 13, megablast 0",
      "0", NULL, NULL, FALSE, 'f', ARG_FLOAT, 0.0, 0, NULL},           /* ARG_THRESHOLD */
    { "Perform gapped alignment (not available with tblastx)", 
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},     /* ARG_GAPPED */
    { "Query Genetic code to use", /* 17 */
      "1", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},           /* ARG_QGENETIC_CODE */
    { "DB Genetic code (for tblast[nx] only)", /* 18 */
      "1", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},           /* ARG_DBGENCODE */
    { "Number of processors to use", /* 19 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},           /* ARG_THREADS */
    { "SeqAlign file",          /* 20 */
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},      /* ARG_ASNOUT */
    { "Believe the query defline", /* 21 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},        /* ARG_BELIEVEQUERY */
    { "Matrix",                 /* 22 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},  /* ARG_MATRIX */
    { "Word size, default if zero (blastn 11, megablast 28, "
        "all others 3)", /* 23 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},            /* ARG_WORDSIZE */
    { "Effective length of the database (use zero for the real size)", 
      "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},          /* ARG_DBSIZE */
    { "Number of best hits from a region to keep. Off by default.\nIf used a value of 100 is recommended.  Very high values of -v or -b is also suggested", 
      "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},            /* ARG_BESTHITS */
    { "0 for multiple hit, 1 for single hit (does not apply to blastn)",
       "0",  NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},           /* ARG_MULTIPLEHITS */
    { "Effective length of the search space (use zero for the real size)", 
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},           /* ARG_SEARCHSP */
    { "Query strands to search against database (for blast[nx], and tblastx)\n"
      "       3 is both, 1 is top, 2 is bottom", 
      "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},             /* ARG_STRAND */
    { "Produce HTML output",    /* 29 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},         /* ARG_HTML */
#ifdef BLAST_CS_API
    { "Restrict search of database to results of Entrez2 lookup", 
      NULL, NULL, NULL, TRUE, 'u', ARG_STRING, 0.0, 0, NULL},          /* ARG_ENTREZQ */
#else
    { "Restrict search of database to list of GI's",             
      NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},          /* ARG_GILIST */
#endif
    {"Use lower case filtering of FASTA sequence", 
     NULL, NULL, NULL, TRUE, 'U', ARG_BOOLEAN, 0.0, 0, NULL},          /* ARG_LCASE */
    { "X dropoff value for ungapped extensions in bits (0.0 invokes default "
      "behavior)\n      blastn 20, megablast 10, all others 7", 
      "0.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},         /* ARG_XDROP_UNGAPPED */       
    { "X dropoff value for final gapped alignment in bits " 
      "(0.0 invokes default behavior)\n"
      "      blastn/megablast 100, tblastx 0, all others 25",
      "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},             /* ARG_XDROP_FINAL */
#ifdef BLAST_CS_API
    { "RPS Blast search",            /* 34 */
      "F", NULL, NULL, FALSE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},          /* ARG_RPSBLAST */
#else
    { "PSI-TBLASTN checkpoint file", /* 34 */
      NULL, NULL, NULL, TRUE, 'R', ARG_FILE_IN, 0.0, 0, NULL},         /* ARG_PSITCHKPNT */
#endif
    { "MegaBlast search",       /* 35 */
      "F", NULL, NULL, FALSE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},         /* ARG_USEMEGABLAST */
    { "Location on query sequence",/* 36 */
      NULL, NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},          /* ARG_QUERYLOC */
    { "Multiple Hits window size, default if zero (blastn/megablast 0, "
        "all others 40", /* 37 */
      "0", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},             /* ARG_WINDOW */
    { "Frame shift penalty (OOF algorithm for blastx)", 
      "0", NULL, NULL, FALSE, 'w', ARG_INT, 0.0, 0, NULL},             /* ARG_FRAMESHIFT */
    { "Length of the largest intron allowed in a translated nucleotide "
      "sequence when "
      "linking multiple distinct alignments. (0 invokes default behavior; a "
      "negative value disables linking.)", 
      "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL},             /* ARG_INTRON */
/*--KM
   seems ok to add another param b/c NUMARG is defined based on 
    sizeof(myargs) itself
   made optional=TRUE but this may change?
*/
#ifndef BLAST_CS_API
    { "Number of concatenated queries, for blastn and tblastn", 
      "0", NULL, NULL, TRUE, 'B', ARG_INT, 0.0, 0, NULL},               /* ARG_NUMQUERIES */
#ifndef BLASTALL_TOOLS_ONLY
    { "Force use of the legacy BLAST engine", 
      "F", NULL, NULL, TRUE, 'V', ARG_BOOLEAN, 0.0, 0, NULL},              /* ARG_FORCE_OLD */
#endif  /* BLASTALL_TOOLS_ONLY */
#endif
    { "Use composition-based score adjustments for blastp or tblastn:\n"                /* ARG_COMP_BASED_STATS */
      "      As first character:\n"
      "      D or d: default (equivalent to T)\n"
      "      0 or F or f: no composition-based statistics\n"
      "      2 or T or t: Composition-based score adjustments as in "
      "Bioinformatics 21:902-911,\n"
      "      1: Composition-based statistics as in "
      "NAR 29:2994-3005, 2001\n"
      "          2005, conditioned on sequence properties\n"
      "      3: Composition-based score adjustment as in "
      "Bioinformatics 21:902-911,\n"
      "          2005, unconditionally\n"
      "      For programs other than tblastn, must either be absent "
      "or be D, F or 0.\n     "
      "      As second character, if first character is "
      "equivalent to 1, 2, or 3:\n"
      "      U or u: unified p-value combining alignment p-value "
      "and compositional p-value in round 1 only\n",
      "D", NULL, NULL, FALSE, 'C', ARG_STRING, 0.0, 0, NULL},
    { "Compute locally optimal Smith-Waterman alignments "
        "(This option is only\n"
      "      available for gapped tblastn.)",                          /* ARG_SMITH_WATERMAN */
      "F", NULL, NULL, FALSE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
#ifdef ALLOW_FULL_SMITH_WATERMAN
    { "Compute only Smith-Waterman alignments (new engine only)",
      "F", NULL, NULL, FALSE, 'h', ARG_BOOLEAN, 0.0, 0, NULL},         /* ARG_SMITH_WATERMAN_ALL */
#endif
};

/* Needed for Mega BLAST only */
#define	MAX_NUM_QUERIES	16383 /* ==	1/2	INT2_MAX */

/*
* these variables were local to blastall's main function,
* they were made global so the main function can be split up.
*/
AsnIoPtr aip, xml_aip;
BioseqPtr fake_bsp = NULL, query_bsp, bsp;
BioSourcePtr source;
BLAST_MatrixPtr	matrix;
Int4Ptr	PNTR txmatrix;
BLAST_OptionsBlkPtr	options;
BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
BlastPruneSapStructPtr prune;
Boolean	db_is_na, query_is_na, show_gi,	believe_query=FALSE;
Boolean	html = FALSE;
CharPtr	params_buffer=NULL;
Int4 number_of_descriptions, number_of_alignments;
SeqAlignPtr	 seqalign;
SeqAnnotPtr	seqannot = NULL;
SeqEntryPtr	sep;
TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
Uint1 align_type, align_view, err_ticket;
Uint4 align_options, print_options;
ValNodePtr mask_loc, mask_loc_start	= NULL,	vnp, next_mask_loc = NULL;
ValNodePtr other_returns, error_returns;
CharPtr	blast_program, blast_database, blast_inputfile,	blast_outputfile;
FILE *infp=NULL, *outfp=NULL;
Char buf[256] =	{ '\0' } ;
/* Mega	BLAST related variables	*/
SeqAlignPtr	sap, next_seqalign,	PNTR seqalignp;
Int4 num_bsps, m_index;
SeqLocPtr last_mask, mask_slp, slp = NULL, tmp_slp;
Int2 ctr = 1;
Char prefix[2];
Boolean	done = TRUE;
int	(LIBCALLBACK *handle_results)(VoidPtr srch);	   
Int4 from =	0, to =	-1;
Uint4 num_queries;		   /*--KM for concatenated queries in blastn, tblastn */
Uint4 num_iters;
Uint4 sap_iter;
SeqAlignPtr	curr_seqalign;
SeqAlignPtrArray sap_array;				   /*--KM for separating seqaligns to test concat printing,	temporary?*/
SeqAnnotPtr	curr_seqannot;
SeqAnnotPtrArray seq_annot_arr;
Uint4 bsp_iter;
BspArray fake_bsp_arr;	   /*--KM the array	of fake_bsps for indiv.	queries	*/
SeqLocPtr PNTR lcase_mask_arr =	NULL;	  /* AM: information about lower case masked parts of queries */
Boolean	concat_done, nuc_concat;
QueriesPtr mult_queries	= NULL;	   /*--KM, AM: stores information related to
								   query multipolexing, to put in search */
BioseqPtr curr_bsp;

/* AM: Support for query multiplexing. */
Uint4 num_spacers;
ValNodePtr orig_mask_loc = NULL;

ValNodePtr sep_list;	/**< a list	of all sep data	structures */
ValNodePtr cur_sep = NULL;		/**< tracks	the	current	position in	sep_list */
ValNodePtr result_id_list;		/**< a list	of query indices that had successful results */
ValNodePtr result_bsp_list;	/**< List of ByteStorePtrs, each of which contains all results (SeqAligns+Bioseqs) for a query */
int	debug;			/**< print debugging messages when set to true */
int	resume_run = 0;	  /**<choose append	to outputfile when true	*/ 
Int8* query_adj_array;	/**< Array of query	length adjustments (one	per	query) */
Int8* db_len_array;	/**< Array of effective	database lengths (one per query) */
Int8** stats_array;		/**< Array of per-query search statistics */
Boolean strict_output_conformance = TRUE;	/**< Set to true if the Altschul paper reference must be given instead of the mpiBLAST paper reference */

Int4 global_dbseq_num = 0;
Int8 global_db_len = 0;

int fakebioseq_start = -1;
int fakebioseq_end = -1;
int blastengine_start = -1;
int blastengine_end = -1;
int addaligninfo_start = -1;
int addaligninfo_end = -1;
int asnwrite_start = -1;
int asnwrite_end = -1;

/* role=0 - master; role=1 - worker */
Int2 initBLAST (int mpi_mode) 

{
	StringCpy(buf, "blastall ");
	StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf));
	if (! GetArgs (buf, NUMARG, myargs)) {
		return (1);
	}

	UseLocalAsnloadDataAndErrMsg ();

	if (! SeqEntryLoad())
		return 1;

	ErrSetMessageLevel(SEV_WARNING);

	blast_program = myargs[ARG_PROGRAM].strvalue;

#ifdef BLAST_CS_API
	/* For RPS Blast - anything not "blastp" - is "tblastn" */    
	if(myargs[ARG_RPSBLAST].intvalue) {
		if(StringICmp(blast_program, "blastp")) {
			StringCpy(blast_program, "blastx");
		}
	}
#endif

	blast_database = myargs[ARG_DB].strvalue;
	blast_inputfile = myargs[ARG_QUERY].strvalue;
	blast_outputfile = myargs[ARG_OUT].strvalue;

	if (myargs[ARG_HTML].intvalue)
		html = TRUE;

	if(!query_in_mem) {
		if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
			ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", blast_inputfile);
			return (1);
		}
	}

	align_view = (Int1) myargs[ARG_FORMAT].intvalue;
	outfp = NULL;
	if (align_view != 7 && align_view != 10 && align_view != 11 && blast_outputfile != NULL) {
		const char* mode = (resume_run == 0) ? "w" : "a";
		if ((outfp = FileOpen(blast_outputfile, (char*)mode) ) == NULL) {
			ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", blast_outputfile);
			return (1);
		}
	}

	if (StringCmp("filter", blast_program) == 0) {
		BlastGetMaskingLoc(infp, outfp, myargs[ARG_FILTER].strvalue);
		FileClose(outfp);
		FileClose(infp);	
		return 0;
	}

	align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);

	if(align_view < 7) {
		if (StringICmp("blastx", blast_program) == 0) {
			if (align_view != 0) {
				ErrPostEx(SEV_FATAL, 1, 0, "This option is not available with blastx");
				return 1;
			}
		} else if (StringICmp("tblastx", blast_program) == 0) {
			if (align_view != 0) {
				ErrPostEx(SEV_FATAL, 1, 0, "This option is not available with tblastx");
				return 1;
			}
		}
	}

	believe_query = FALSE;
	if (myargs[ARG_BELIEVEQUERY].intvalue != 0)
		believe_query = TRUE;

//	if (believe_query == FALSE && (myargs[ARG_ASNOUT].strvalue || align_view == 10 || align_view ==11)) {
//		ErrPostEx(SEV_FATAL, 1, 0, "-J option must be TRUE to produce a SeqAlign file");
//	}

    options = BLASTOptionNewEx(blast_program, (Boolean) myargs[ARG_GAPPED].intvalue, (Boolean) myargs[ARG_USEMEGABLAST].intvalue);
	if (options == NULL)
		return 3;

#ifdef BLAST_CS_API
	if(myargs[ARG_RPSBLAST].intvalue) 
		options->is_rps_blast = TRUE;
#endif

//	if (align_view == 8 && options->is_megablast_search) {
//		options->output = (VoidPtr) outfp;
//		handle_results = MegaBlastPrintAlignInfo;
//	} else 
//		handle_results = NULL;

    handle_results = NULL;

    BLASTOptionSetGapParams(options, myargs[ARG_MATRIX].strvalue, 0, 0); 
    options->kappa_expect_value =
        options->expect_value  = (Nlm_FloatHi) myargs[ARG_EVALUE].floatvalue;
    number_of_descriptions = myargs[ARG_DESCRIPTIONS].intvalue;	
    number_of_alignments = myargs[ARG_ALIGNMENTS].intvalue;	
    options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);

    if (StringICmp("blastn", blast_program) == 0) {
        options->penalty = myargs[ARG_MISMATCH].intvalue;
        options->reward = myargs[ARG_MATCH].intvalue;
        if (options->reward > 1) {
           /* Scale the default values for gap costs; will be overridden
              later, if command line values are non-zero */
           options->gap_open *= options->reward;
           options->gap_extend *= options->reward;
        }
    } else {
        if ((Int4)myargs[ARG_THRESHOLD].floatvalue != 0) {
            options->threshold_second = (Int4)myargs[ARG_THRESHOLD].floatvalue;
        }
    }
    
    if (myargs[ARG_GAPOPEN].intvalue >= 0)
        options->gap_open = myargs[ARG_GAPOPEN].intvalue;
    if (myargs[ARG_GAPEXT].intvalue >= 0)
        options->gap_extend = myargs[ARG_GAPEXT].intvalue;
    if (myargs[ARG_XDROP].intvalue != 0)
        options->gap_x_dropoff = myargs[ARG_XDROP].intvalue;

    /* use one-hit if specified or it's a blastn search */
    if ( (myargs[ARG_MULTIPLEHITS].intvalue == 1) || (StringICmp("blastn", blast_program) == 0 ) )
      {
        options->two_pass_method  = FALSE;
        options->multiple_hits_only  = FALSE;
      }
    /* otherwise, use two-hit */
    else
      { 
        /* all other inputs, including the default 0 use 2-hit method */
        options->two_pass_method  = FALSE;
        options->multiple_hits_only  = TRUE;
      }
    
    if(myargs[ARG_XDROP_FINAL].intvalue != 0) 
        options->gap_x_dropoff_final = myargs[ARG_XDROP_FINAL].intvalue;

    if (StringICmp(myargs[ARG_FILTER].strvalue, "T") == 0) {
        if (StringICmp("blastn", blast_program) == 0)
            options->filter_string = StringSave("D");
        else
            options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[ARG_FILTER].strvalue);
    }
    
    show_gi = (Boolean) myargs[ARG_SHOWGIS].intvalue;

    options->genetic_code = myargs[ARG_QGENETIC_CODE].intvalue;
    options->db_genetic_code = myargs[ARG_DBGENCODE].intvalue;
    options->number_of_cpus = myargs[ARG_THREADS].intvalue;
    if (myargs[ARG_WORDSIZE].intvalue != 0) {
        options->wordsize = myargs[ARG_WORDSIZE].intvalue;
    }
    
    if (options->is_megablast_search) {
       options->cutoff_s2 = options->wordsize*options->reward;
    }

    options->db_length = (Int8) myargs[ARG_DBSIZE].floatvalue;
    
    options->hsp_range_max  = myargs[ARG_BESTHITS].intvalue;
    if (options->hsp_range_max != 0)
        options->perform_culling = TRUE;
    if (myargs[ARG_SEARCHSP].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[ARG_SEARCHSP].floatvalue;

    if ((0 != StringICmp("tblastn", blast_program) &&
         0 != StringICmp("blastp", blast_program)) ||
        !options->gapped_calculation) {
        /* Set some gapped tblastn-specific options to the correct
         * defaults for non-tblastn or non-gapped modes of operation.
         */
        options->tweak_parameters = eNoCompositionBasedStats;
        options->smith_waterman = 0;
        options->unified_p = 0;
        
        switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
        case '0': case 'D': case 'd': case 'F': case 'f':
            options->tweak_parameters = eNoCompositionBasedStats;
            break;
        default:
            ErrPostEx(SEV_FATAL, 1, 0,
               "Invalid option -C: only gapped blastp or gapped tblastn "
               "may use composition based statistics.");
            break;
        }
        if(myargs[ARG_SMITH_WATERMAN].intvalue) {
            ErrPostEx(SEV_FATAL, 1, 0,
               "Invalid option -s: Smith-Waterman alignments are only "
               "available for gapped blastp or gapped tblastn.");
        }
    } else {
        /* Set options specific to gapped tblastn and blastp */
        switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
        case '0': case 'F': case 'f':
            options->tweak_parameters = eNoCompositionBasedStats;
            break;
        case 'D': case 'd':
        case '1': case 'T': case 't':
            options->tweak_parameters = eCompositionBasedStats;
            break;
        case '2':
            ErrPostEx(SEV_WARNING, 1, 0, "the -C 2 argument "
                      "is currently experimental\n");
            options->tweak_parameters = eCompositionMatrixAdjust;
            break;
        case '3':
            ErrPostEx(SEV_WARNING, 1, 0, "the -C 3 argument "
                      "is currently experimental\n");
            options->tweak_parameters = eCompoForceFullMatrixAdjust;
        break;
        default:
            ErrPostEx(SEV_FATAL, 1, 0, "invalid argument for composition-"
                      "based statistics; see -C options\n");
            break;
        }
	if (options->tweak_parameters > 0) {
            switch (myargs[ARG_COMP_BASED_STATS].strvalue[1]) {
            case 'U':
            case 'u': 
                if (0 == StringICmp("blastp", blast_program)) {
                    options->unified_p = 1;
                    ErrPostEx(SEV_WARNING, 1, 0, "unified p-values "
                              "are currently experimental\n");
                } else {
                    ErrPostEx(SEV_FATAL, 1, 0, "unified p-values "
                              "are currently only available for blastp\n");
                }
                break;
	  case '\0':
            break;
	  default:
            ErrPostEx(SEV_WARNING, 1, 0, "unrecognized second character"
                      "in value of -t, ignoring it\n");
            break;
	  }
	}
        options->smith_waterman =
            (Boolean) myargs[ARG_SMITH_WATERMAN].intvalue;
    }
    if (options->tweak_parameters > 1) {
        /* Compositionally adjusted score matrices are being used, and
         * these can improve evalue, so relax the evalue cutoff for
         * the preliminary alignments.  (Note that traditional
         * composition based statistics can only make evalues larger.)
         */
        options->expect_value *= EVALUE_EXPAND;
    }

    options->strand_option = myargs[ARG_STRAND].intvalue;

    if(myargs[ARG_XDROP_UNGAPPED].floatvalue != 0.0) {
        options->dropoff_2nd_pass  = myargs[ARG_XDROP_UNGAPPED].floatvalue;
        if(options->dropoff_1st_pass > options->dropoff_2nd_pass)
            options->dropoff_1st_pass = options->dropoff_2nd_pass;
    }

    if (myargs[ARG_WINDOW].intvalue != 0)
        options->window_size = myargs[ARG_WINDOW].intvalue;

    print_options = 0;
    align_options = 0;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    if (StringICmp("blastx", blast_program) == 0) {
        align_options += TXALIGN_BLASTX_SPECIAL;
    }
    if (show_gi) {
        align_options += TXALIGN_SHOW_GI;
        print_options += TXALIGN_SHOW_GI;
    }
    if (myargs[ARG_GAPPED].intvalue == 0 || StringICmp("tblastx", blast_program) == 0)
        print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    if (align_view) {
        align_options += TXALIGN_MASTER;
        if (align_view == 1 || align_view == 3)
            align_options += TXALIGN_MISMATCH;
        if (align_view == 3 || align_view == 4 || align_view == 6)
            align_options += TXALIGN_FLAT_INS;
        if (align_view == 5 || align_view == 6)
            align_options += TXALIGN_BLUNT_END;
    } else {
        align_options += TXALIGN_MATRIX_VAL;
        align_options += TXALIGN_SHOW_QS;
    }
    
    if (html) {
        align_options += TXALIGN_HTML;
        print_options += TXALIGN_HTML;
    }

#ifdef BLAST_CS_API
    if(myargs[ARG_ENTREZQ].strvalue)
        options->entrez_query = StringSave(myargs[ARG_ENTREZQ].strvalue);
#else    
    if (myargs[ARG_GILIST].strvalue) {
        options->gifile = StringSave(myargs[ARG_GILIST].strvalue);
    }
#endif
    
    /* 
       Out-of-frame option is valid only for blastx, tblastn and 
       psitblastnsearches
    */

    if(myargs[ARG_FRAMESHIFT].intvalue > 0) {
        if (!StringICmp("blastx", blast_program) || 
            !StringICmp("tblastn", blast_program)||
	    !StringICmp("psitblastn", blast_program)) {
           if (!StringICmp("blastx", blast_program)) {
              options->is_ooframe = TRUE;
              options->shift_pen = myargs[ARG_FRAMESHIFT].intvalue;
           }
        }
    }
        
    /* Input longest intron length is in nucleotide scale; in the lower level
       code it will be used in protein scale */
    options->longest_intron = myargs[ARG_INTRON].intvalue;

    aip = NULL;
    if (myargs[ARG_ASNOUT].strvalue != NULL) {
        if ((aip = AsnIoOpen (myargs[ARG_ASNOUT].strvalue,"w")) == NULL) {
                ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[ARG_ASNOUT].strvalue);
                return 1;
        }
    }
    else if (align_view == 10 || align_view == 11) 
    {
        const char* mode = (align_view == 10) ? "w" : "wb";
        if ((aip = AsnIoOpen (blast_outputfile, (char*) mode)) == NULL) {
                ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[ARG_ASNOUT].strvalue);
                return 1;
        }
    }

    if(align_view < 7) {
       if (html) {
          fprintf(outfp, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
          fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                  "VLINK=\"#660099\" ALINK=\"#660099\">\n");
          fprintf(outfp, "<PRE>\n");
       }
    } else if (align_view == 7 ) {
        xml_aip = AsnIoOpen(blast_outputfile, "wx");
    }

#ifndef BLAST_CS_API
    if(align_view >= 7 && myargs[ARG_NUMQUERIES].intvalue > 1)
    {
      ErrPostEx(SEV_FATAL, 1, 0, 
                 "blast: Query concatenation is currently not supported with -m > 7");
      return 1;
    }
#endif


                  /* Futamura: Setting up the psitblastn options */
#ifndef BLAST_CS_API
    if (NULL != myargs[ARG_PSITCHKPNT].strvalue) {
          options->recoverCheckpoint = TRUE;
          options->freqCheckpoint = TRUE;
    }
    options->CheckpointFileName=myargs[ARG_PSITCHKPNT].strvalue;
#endif

#ifdef BLAST_CS_API
    if (align_view < 7)
    	bl3hp = BNETInitializeBlast(blast_database, blast_program, outfp, 
                                db_is_na, options->is_rps_blast, html, TRUE);
    else
    	bl3hp = BNETInitializeBlast(blast_database, blast_program, outfp, 
                                db_is_na, options->is_rps_blast, html, FALSE);
#endif

    /*--KM get number of queries for concatenated blastn/tblastn queries */

#ifndef BLAST_CS_API
    options->NumQueries=myargs[ARG_NUMQUERIES].intvalue;  
#endif

    num_queries = options->NumQueries;
    if (num_queries>0 && 
	!( (StringICmp("blastn",  blast_program) == 0) || 
	   (StringICmp("tblastn", blast_program) == 0)   ) ) {

	ErrPostEx(SEV_FATAL, 1, 0, "blast: Can't concat with program %s\n", myargs[ARG_PROGRAM].strvalue);
       return 1;
    }
    
    /* AM: Query concatenation is not consistent with ungapped search */
    if( num_queries > 0 && !myargs[ARG_GAPPED].intvalue )
    {
      ErrPostEx(SEV_FATAL, 1, 0, 
                 "blast: Query concatenation is inconsistent with ungapped search\n" );
      return 1;
    }
    if( !myargs[ARG_GAPPED].intvalue &&
        0 == StringCmp("psitblastn", blast_program ) ) {
      ErrPostEx(SEV_FATAL, 1, 0,"blast: Ungapped alignment is not appropriate "
                "for PSI-tBLASTn.\n" );
    }

    /* --KM set bool value if DNA and concat needed, need for Fasta->seq functions */
    if (num_queries>0 && query_is_na == TRUE) {
        nuc_concat = TRUE;
    } else {
        nuc_concat = FALSE;
    }
 
    /* --- Main loop over all FASTA entries in the input file ---- */

    concat_done = FALSE;	/*--KM */

	/* mpiBLAST */
	/* set the effective query and db lengths */
	options->query_adjustments = query_adj_array;
	options->effective_db_lengths = db_len_array;

	return 0;
}


/**
* function to parse the FastA format query file, load the queries, and count the number loaded
*/
int	loadQueries( void )	
{
	int query_count = 0;
	while (TRUE) {
		if (!options->is_megablast_search) {
			if((Boolean)myargs[ARG_LCASE].intvalue) {
				sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, 
					NULL, NULL, &options->query_lcase_mask);
			} else {
				sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query );
			}

			if(sep == NULL){
				break; /* no more queries, can go to finish with next break */
			}

			ValNodeAddPointer( &sep_list, 0, sep );
			query_count++;
		}
	}
	return query_count;
}

int	loadQueriesEx( void )	
{
	int query_count = 0;
	while (TRUE) {
		if (!options->is_megablast_search) {
			if((Boolean)myargs[ARG_LCASE].intvalue) {
				sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, 
					NULL, NULL, &options->query_lcase_mask);
			} else {
				sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query );
			}

			if(sep == NULL){
				break; /* no more queries, can go to finish with next break */
			}

			cQueryMAddQueryEntry(query_count, sep);
			query_count++;
		}
	}
	return query_count;
}

int	loadQueriesFromBuffer( char* buffer )	
{
	CharPtr last_char;
	int query_count = 0;

	last_char = buffer;
	while (TRUE) {
		if (!options->is_megablast_search) {
			if((Boolean)myargs[ARG_LCASE].intvalue) {
				sep = FastaToSeqBuffForDb (last_char, &last_char, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
			} else {
				sep = FastaToSeqBuffEx(last_char, &last_char, query_is_na, NULL, believe_query );
			}

			if(sep == NULL){
				break; /* no more queries, can go to finish with next break */
			}

			ValNodeAddPointer( &sep_list, 0, sep );
			query_count++;
		}
	}
	return query_count;
}

int	loadQueriesFromBufferEx( char* buffer )	
{
	CharPtr last_char;
	int query_count = 0;

	last_char = buffer;
	while (TRUE) {
		if (!options->is_megablast_search) {
			if((Boolean)myargs[ARG_LCASE].intvalue) {
				sep = FastaToSeqBuffForDb (last_char, &last_char, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
			} else {
				sep = FastaToSeqBuffEx(last_char, &last_char, query_is_na, NULL, believe_query );
			}

			if(sep == NULL){
				break; /* no more queries, can go to finish with next break */
			}

			cQueryMAddQueryEntry(query_count, sep);
			query_count++;
		}
	}
	return query_count;
}

SeqEntryPtr	bufferToQueryEntry( char* buffer )	
{
	CharPtr last_char;

	last_char = buffer;
	sep = NULL;
	if (!options->is_megablast_search) {
		if((Boolean)myargs[ARG_LCASE].intvalue) {
			sep = FastaToSeqBuffForDb (last_char, &last_char, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
		} else {
			sep = FastaToSeqBuffEx(last_char, &last_char, query_is_na, NULL, believe_query );
		}
	}

	return sep;
}

int countQueriesInBuffer( char *buffer, int buf_size ) {
	int i=0;
	int count = 0;
	for(i=0; i<buf_size; i++) {
		if(buffer[i] == '>') {
			count++;
		}
	}

	return count;
}

int updateBlastDB(const char* db_name, int len) {
	if(blast_database != NULL) {
		MemFree(blast_database);
	}

	blast_database = (char*)MemNew(len + 1);
	memcpy(blast_database, db_name, len);
	blast_database[len] = 0;
	
	return 0;
}

void cleanupQueries( void )	{
	ValNodePtr cur_sep = sep_list;
	for( ; cur_sep != NULL; cur_sep = cur_sep->next )
		cur_sep->data.ptrvalue = SeqEntryFree( (SeqEntryPtr)cur_sep->data.ptrvalue );
	sep_list = ValNodeFree( sep_list );
}

/**
* function to select a certain query sequence given the order of the query which were 
* loaded by loadQueries function above.
*/
SeqEntryPtr	getQuery( int queryI )
{
	int qI = 0;

	if(use_query_map) {
		return getQueryEx(queryI);
	}
	
	if(sep_list != NULL){
		ValNodePtr curr = sep_list;
		for( ; qI < queryI; qI++ ){
			curr = curr->next;
		}
		/* set the query adjustment to the current query */
		options->current_queryI = queryI;
		return curr->data.ptrvalue;
	}
	return NULL;
}

SeqEntryPtr	getQueryEx( int queryI ) {
	options->current_queryI = queryI;
	return cQueryMGetQueryEntry(queryI);
}

/**
* function to load values into fake_bsp for the current bioseq
* queryI is only valid if preloaded is true.
* preloading only works with standard blast searches, not megablast
* returns 0 on success
* returns > 0 on error condition
* returns -1 if there are no more bioseqs to load
*/
int	getFakeBioseq( Boolean preloaded, int queryI ){
	if (options->is_megablast_search) {
		StrCpy(prefix, "");
		slp = NULL;
		num_bsps = 0;
		done = TRUE;
		SeqMgrHoldIndexing(TRUE);
		mask_slp = last_mask = NULL;
		while ((sep=FastaToSeqEntryForDb(infp, query_is_na, NULL,
			believe_query, prefix, &ctr, 
			&mask_slp)) != NULL) {
				if ((Boolean)myargs[ARG_LCASE].intvalue) {
					if (mask_slp) {
						if (!last_mask)
							options->query_lcase_mask = last_mask = mask_slp;
						else {
							last_mask->next = mask_slp;
							last_mask = last_mask->next;
						}
						mask_slp = NULL;
					}
				} else {
					mask_slp = SeqLocSetFree(mask_slp);
				}
				query_bsp = NULL;
				if (query_is_na) 
					SeqEntryExplore(sep, &query_bsp, FindNuc);
				else
					SeqEntryExplore(sep, &query_bsp, FindProt);

				if (query_bsp == NULL) {
					ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
					return 2;
				}

				/* Only for the first query */
				if (num_bsps == 0) {
					to = MIN(to, query_bsp->length - 1);

					/* -1 means end of sequence */
					if (to < 0)
						to = query_bsp->length - 1;
					if (from >= query_bsp->length || to < 0) {
						ErrPostEx(SEV_FATAL, 1, 0, 
							"Location outside of the query sequence range\n");
						return 3;
					}
					slp = SeqLocIntNew(from, to, options->strand_option, 
						SeqIdFindBest(query_bsp->id, SEQID_GI));
				} else 
					ValNodeAddPointer(&slp, SEQLOC_WHOLE,
					SeqIdDup(SeqIdFindBest(query_bsp->id,
					SEQID_GI)));
				num_bsps++;
				if (num_bsps >= MAX_NUM_QUERIES) {
					done = FALSE;
					break;
				}
				/*sep = MemFree(sep);*/ /* Do not free the underlying Bioseq */
		}
		SeqMgrHoldIndexing(FALSE);
		if (num_bsps == 0) 
			return -1;
	} else {

		/* not megablast */

		/*--KM make array of fake_bsp's if concat. query */
		if (concat_done)
			return -1;
		if (num_queries > 0)  {
			fake_bsp_arr = (BspArray) MemNew(sizeof(BioseqPtr)*num_queries);

			if( myargs[ARG_LCASE].intvalue )
				lcase_mask_arr = (SeqLocPtr PNTR)MemNew( sizeof( SeqLocPtr )*num_queries );
		}
		num_iters = (num_queries>0) ? num_queries : 1;
		for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {

			if( preloaded ){
				/*					if(cur_sep == NULL){
				sep = NULL;
				}else{
				sep = (SeqEntryPtr)cur_sep->data.ptrvalue;
				cur_sep = cur_sep->next;
				}
				*/
				sep = getQuery( queryI++ );	/** increment for query concatenation */
			}else{
				if(myargs[ARG_LCASE].intvalue) {	  
					/* AM: query multiplexing */
					if( !num_queries )
						sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
					else
						sep = FastaToSeqEntryInternalEx( infp, FASTA_FILE_IO, NULL, query_is_na, NULL, believe_query,
						NULL, NULL, NULL, lcase_mask_arr + bsp_iter );

				} else {	  
					sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query);
				}
			}

			/* if concat and num_queries has not been reached and sep is NULL, crap out */
			if (sep == NULL && bsp_iter < num_queries) {   /* implies num_queries>0 */
				ErrPostEx(SEV_FATAL, 1, 0, "blast: Only %d queries found!\n", bsp_iter); 
				return (1);
			}
			if(sep == NULL)
				return -1;	/* go to finish, all bioseqs have been parsed */

			query_bsp = NULL;
			if (query_is_na) {
				SeqEntryExplore(sep, &query_bsp, FindNuc);
			} else {
				SeqEntryExplore(sep, &query_bsp, FindProt);
			}

			if (query_bsp == NULL) {
				ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
				return 2;
			}

			if (num_queries>0) {
				*(fake_bsp_arr + bsp_iter) = query_bsp;
			}
		}
		if ( (sep == NULL && num_queries ==0) || (num_queries>0 && concat_done) )
			return -1;  /* go to finish */

		/* --KM */
		if (num_queries>0) {
			concat_done = TRUE;   /* --KM to prevent futher looping */

			/* AM: Determine the number of query separators. */
			num_spacers = GetNumSpacers( options, believe_query, fake_bsp_arr );

			if( num_spacers%2 ) ++num_spacers;

			/* --KM make the concatenated fake_bsp */
			/* AM: Added num_spacers. */
			if( query_is_na )
				fake_bsp = (BioseqPtr)
				BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers);
			else
				fake_bsp = (BioseqPtr)
				BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers);

			/* construct the MultQueries struct here*/
			mult_queries = (QueriesPtr) BlastMakeMultQueries(fake_bsp_arr, num_queries, query_is_na, num_spacers, lcase_mask_arr);
		} else {

			if(believe_query){
				fake_bsp = query_bsp;
			}
			else{ 
				fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);
			}
		}
		err_ticket = BlastSetUserErrorString(NULL, query_bsp->id, believe_query);

		/* If fake_bsp created mask should be updated to use it's id */
		/* AM: query multiplexing */
		if( !mult_queries )
			BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, fake_bsp->id);
		else for( bsp_iter = 0; bsp_iter < num_iters; ++bsp_iter )
			if( mult_queries->LCaseMasks )
				BLASTUpdateSeqIdInSeqInt( mult_queries->LCaseMasks[bsp_iter],
				mult_queries->FakeBsps[bsp_iter]->id );

		source = BioSourceNew();
		source->org = OrgRefNew();
		source->org->orgname = OrgNameNew();
		source->org->orgname->gcode = options->genetic_code;
		ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
	}
	return 0;
}

int	performBLAST(){

	global_fp = outfp;

	other_returns = NULL;
	error_returns = NULL;

	if (options->is_megablast_search) {
		seqalignp = BioseqMegaBlastEngineByLoc(slp, blast_program,
			blast_database, options, &other_returns, 
			&error_returns, 
			align_view < 7 ? tick_callback : NULL,
			NULL, NULL, 0, handle_results);
		seqalign = NULL;
		for (m_index=0; m_index<num_bsps; m_index++) { 
			if (seqalignp && seqalignp[m_index]) {
				if (seqalign == NULL) 
					sap = seqalign = seqalignp[m_index];
				else
					sap->next = seqalignp[m_index];
				while (sap->next != NULL)
					sap = sap->next;
			}
		}
	   seqalignp = MemFree(seqalignp);

	} else if (!myargs[ARG_QUERYLOC].strvalue) {       
//		if(!hijack_default) {
			/* KM added mult_queries param */
			seqalign = BioseqBlastEngineWithCallbackMult(fake_bsp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, handle_results, mult_queries);
//		}
	} else { /* Location on query provided */
		to = MIN(to, fake_bsp->length - 1);

		/* -1 means end of sequence */
		if (to < 0)
			to = fake_bsp->length - 1;
		if (from >= fake_bsp->length || to < 0) {
			ErrPostEx(SEV_FATAL, 1, 0, 
				"Location outside of the query sequence range\n");
			return 3;
		}
		slp = SeqLocIntNew(from, to, options->strand_option, 
			fake_bsp->id);
		seqalign = BioseqBlastEngineByLocWithCallbackMult(slp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, NULL, NULL, 0, handle_results, mult_queries);

	}

	return 0;
}


void writePartialSubjectBioseq( SeqAlignPtr sap_p, AsnIoBSPtr aibp );
void writeSearchStatistics( AsnIoBSPtr aibp );


/**
* This function takes results from a blast search and writes
* them into a memory buffer (a ByteStore) as binary ASN.1 results
* Specifically, results are written as:
* 1.  All the SeqAlign structs associated with the results
* 2.  All the partial Bioseq structs used in alignments described
*     in the SeqAligns
* 3.  The array of search statistics -- 12 long long int values
*/
Int2 createBinASNResults( SeqAlignPtr sappy, ByteStorePtr bs_ptr, int queryI	)
{	
	seqalign = sappy;
	global_fp = outfp;

	ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

	tmp_slp = slp;
	if (slp)
		query_bsp = NULL;

	if (getenv("POST_BLAST_CLUSTER_HITS") != NULL)
		BlastClusterHitsFromSeqAlign(seqalign, blast_program, blast_database, 
		options, 0.9, 1.6, 0.5, TRUE);

	if (mask_loc) {
		mask_loc_start = mask_loc;
	}	
	else
	{	/* Could have become non-NUll for last query. */
		mask_loc_start = NULL;
	}

	if (seqalign) {
		if (num_queries > 0) { /* AM: Support for query multiplexing. */
			sap_array = mult_queries->sap_array_data->sap_array;
		}

		while (seqalign) {
			/*if( debug )fprintf( stderr, "Processing next seqalign\n" );*/

			if (!options->is_megablast_search){
				next_seqalign = NULL;
			} else {
				SeqIdPtr sip, next_sip = NULL;

				sap = seqalign;
				sip = TxGetQueryIdFromSeqAlign(seqalign);

				while (sap != NULL) { 
					if (sap->next != NULL) {
						next_sip = TxGetQueryIdFromSeqAlign(sap->next);

						if (SeqIdComp(sip, next_sip) != SIC_YES) {
							next_seqalign = sap->next;
							sap->next = NULL;
						}
					} else{
						next_seqalign = NULL;
					}
					sap = sap->next;
				}

				while (tmp_slp && SeqIdComp(sip, SeqLocId(tmp_slp)) != SIC_YES)
					tmp_slp = tmp_slp->next;
				if (tmp_slp == NULL) /* Should never happen */
					break;
				/* Separate the mask locations list for this query */
				if (!mask_loc && next_mask_loc) {
					mask_loc = next_mask_loc;
					next_mask_loc = NULL;
				}
				if (mask_loc) {
					if (next_mask_loc) {
						mask_loc->next = next_mask_loc;
						mask_loc = next_mask_loc;
					}
					mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
					next_mask_loc = mask_loc;
					while (SeqIdComp(SeqLocId(mask_slp), sip) != SIC_YES) {
						mask_loc = mask_loc->next;
						if (!mask_loc)
							break;
						mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
					}
					if (mask_loc) {
						next_mask_loc = mask_loc->next;
						mask_loc->next = NULL;
					}
				}

			}

			/* create the array of SeqAnnotPtrs, if necessary */
			/*if( debug )fprintf( stderr, "Creating array of SeqAnnotPtrs\n" );*/

			num_iters = (num_queries > 0) ? num_queries : 1;
			for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
				curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
				if ( (num_queries > 0) && (sap_iter == 0) ) {
					seq_annot_arr = (SeqAnnotPtrArray) MemNew(sizeof(SeqAnnotPtr)*num_queries);
				}
				seqannot = SeqAnnotNew();
				seqannot->type = 2;
				/*if( debug )fprintf( stderr, "AddAlignInfoToSeqAnnot\n" );*/
				AddAlignInfoToSeqAnnot(seqannot, align_type);
				seqannot->data = curr_seqalign;

				if (num_queries > 0) {
					*(seq_annot_arr + sap_iter) = seqannot;
				}
			} /* make seqannots over the sap_iters from concat, or the single seqalign */

			{
				/* Write the SeqAlign to a temporary byte store */
				AsnIoBSPtr aibp;
				if( use_binary_asn )
					aibp = AsnIoBSOpen("wb", bs_ptr);
				else
					aibp = AsnIoBSOpen("w", bs_ptr);

				/* Write results to binary ASN.1 objects in memory and add 
				* the pointer to the asn_results_list.  This functionality
				* allows the worker process to directly send the results
				* to the master without the intermediate step of writing them
				* to a file.
				*/
				for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
					/*if( debug )fprintf( stderr, "ASN result output sap_iter: %d\n", sap_iter );*/
					curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;

					if( !SeqAnnotAsnWrite(curr_seqannot, aibp->aip, NULL) ){
						fprintf(stderr, "packSeqAnnot: Error writing SeqAnnot to ByteStore\n");
						return -1;
					}
				}

				/* now write all the bioseqs to the ByteStore */
				for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
					curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;
					writePartialSubjectBioseq( curr_seqannot->data, aibp );
				}
				writeSearchStatistics( aibp );
				AsnIoBSClose( aibp );
			}
			/*if( debug )fprintf( stderr, "deallocate\n" );*/
			for (sap_iter=0; sap_iter < num_queries; sap_iter++) {
				/* upper bound is num_queries, take care not to do this unless concat */
				*(seq_annot_arr + sap_iter) = SeqAnnotFree(*(seq_annot_arr + sap_iter));
			}
			/*--KM free seqalign array and all seqaligns?? */

			if (options->is_megablast_search)
				tmp_slp = tmp_slp->next;

			/* --KM watch for memory leaks */
			if (seqannot &&	num_queries	== 0)
				seqannot = SeqAnnotFree(seqannot);

			if( aip == NULL ){
				seqalign = SeqAlignSetFree(seqalign);
			}

			seqalign = next_seqalign;
		} /* End of loop on all seqaligns */
		if (mask_loc && next_mask_loc)
			mask_loc->next = next_mask_loc;

	} else {         
		if (error_returns != NULL) {
			for (vnp = error_returns; vnp; vnp = vnp->next) {
				BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
			}
			ValNodeFree(error_returns);
		}
	}

	return 0;
}


/**
 * Given a SeqAlign, extract the aligned portion of the subject 
 * Bioseq and write it to a ByteStore.  Also write the original
 * length of the bioseq
 * @param sap_p		The SeqAlignPtr where blast results are stored
 * @param aibp		The ASN I/O structure for the ByteStore
 */
void writePartialSubjectBioseq( SeqAlignPtr sap_p, AsnIoBSPtr aibp )
{
	DataVal dv;
	SeqIdPtr tmp_id;
	while(sap_p !=NULL) {
		Int4 start = -1;
		Int4 stop = -1;
		Uint1 strand = 0;
		BioseqPtr subject_bsp;
		BioseqPtr partial_seq;
		SeqIdPtr bioseq_id = TxGetSubjectIdFromSeqAlign( sap_p );
		get_align_ends( sap_p, bioseq_id, &start, &stop, &strand );
		subject_bsp = BioseqLockById( bioseq_id );

		// copy the used portion of the bioseq
		if( strand == 2 )
			strand = 1;	// always copy and send the forward strand
		partial_seq = BioseqCopy(NULL, bioseq_id, start, stop, strand, FALSE);

		// temporarily set the description and id then
		// write the bioseq to the ByteStore
		partial_seq->descr = subject_bsp->descr;
		tmp_id = partial_seq->id;
		partial_seq->id = subject_bsp->id;
		if( !BioseqAsnWrite(partial_seq, aibp->aip, NULL) ){
			fprintf(stderr, "writePartialSubjectBioseq: Error writing Bioseq to ByteStore\n");
			return;
		}
		partial_seq->id = tmp_id;
		partial_seq->descr = NULL;
		
		// now write the original bioseq length
		dv.bigintvalue = subject_bsp->length;
		AsnOpenStruct( aibp->aip, &lenset_type, &dv );
		AsnWrite( aibp->aip, &len_type, &dv );
		AsnCloseStruct( aibp->aip, &lenset_type, &dv );

		// free the partial bioseq
		BioseqFree( partial_seq );
		
		// unlock and free the subject bioseq
		BioseqUnlock( subject_bsp );
		BioseqFree( subject_bsp );
		
		// move to the next alignment
		sap_p = sap_p->next;
	}
}

/**
 * Write the array of search statistics
 * @param aibp		The ASN I/O structure for the ByteStore
 */
void writeSearchStatistics( AsnIoBSPtr aibp )
{
	AsnTypePtr atp;
	DataVal dv;
	SeqIdPtr anp = NULL;
	int statI;
	for( statI = 0; statI < 12; statI++ ){
		AsnOpenStruct( aibp->aip, &statset_type, &dv );
		dv.bigintvalue = options->stats[statI];
		AsnWrite( aibp->aip, &statdata_type, &dv );
		AsnCloseStruct( aibp->aip, &statset_type, &dv );
	}
}

/**	
* process the data returned by the BLAST search 
* by copying it to local variables 
*/
void processReturnData(){
    BlastErrorPrint(error_returns);

	dbinfo = NULL;
	ka_params = NULL;
	ka_params_gap = NULL;
	params_buffer = NULL;
	mask_loc = NULL;
	matrix = NULL;
	txmatrix = NULL;
	for (vnp=other_returns; vnp; vnp = vnp->next) {
		switch (vnp->choice) {
		case TXDBINFO:
			dbinfo = vnp->data.ptrvalue;
			break;
		case TXKABLK_NOGAP:
			ka_params = vnp->data.ptrvalue;
			break;
		case TXKABLK_GAP:
			ka_params_gap = vnp->data.ptrvalue;
			break;
		case TXPARAMETERS:
			params_buffer = vnp->data.ptrvalue;
			break;
		case TXMATRIX:
			matrix = vnp->data.ptrvalue;
			if (matrix)
			   txmatrix = BlastMatrixToTxMatrix(matrix);
			break;
		case SEQLOC_MASKING_NOTSET:
		case SEQLOC_MASKING_PLUS1:
		case SEQLOC_MASKING_PLUS2:
		case SEQLOC_MASKING_PLUS3:
		case SEQLOC_MASKING_MINUS1:
		case SEQLOC_MASKING_MINUS2:
		case SEQLOC_MASKING_MINUS3:
			ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
			break;
		default:
			break;
		}
	}	
}

void cleanupTheMess(){
	/**	start 'cleanup the mess' code */

	slp = SeqLocSetFree(slp);
	matrix = BLAST_MatrixDestruct(matrix);
	if (txmatrix)
		txmatrix = TxMatrixDestruct(txmatrix);

	init_buff_ex(85);
	dbinfo_head = dbinfo;

	dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);

	if (ka_params) {
		MemFree(ka_params);
	}

	if (ka_params_gap) {
		MemFree(ka_params_gap);
	}

	MemFree(params_buffer);
	free_buff();
	mask_loc = mask_loc_start;
	while (mask_loc) {
		SeqLocSetFree((ValNode*)(mask_loc->data.ptrvalue));
		mask_loc = mask_loc->next;
	}
	ValNodeFree(mask_loc_start);

	if(!believe_query)
		fake_bsp = BlastDeleteFakeBioseq(fake_bsp);

	other_returns = ValNodeFree(other_returns);

#ifndef	BLAST_CS_API
	/* This is freed earlier in client-server case */
	options->query_lcase_mask = SeqLocSetFree(options->query_lcase_mask);
	ReadDBBioseqFetchDisable();
#endif

	if (!options->is_megablast_search) 
		BlastDeleteUserErrorString(err_ticket);

	ObjMgrFreeCache(0);

	/**	end	cleanup	code */
}


/**
* run the blast algorithm in one of several modes:
*  0	Standard BLAST search as performed on a worker node
*  1	Collection of search space and effective HSP length
*/
Int2 runBLAST( int mode, int first_query, int last_query ) {
	/* 
	* Use a 500 kb buffer for the initial byte store. Byte stores
	* grown dynamically and there may be some opportunity for optimizing this
	* initial size.
	*/
	const unsigned init_bsp_size = 512000;
	ByteStorePtr bs_ptr;

	int queryI = first_query;
	int query_increment = num_queries <= 0 ? 1 : num_queries;
	cur_sep = sep_list;

	if( mode == COLLECT_STATS_MODE ){
		options->calculate_statistics_and_exit = TRUE;
	}else
		options->calculate_statistics_and_exit = FALSE;

	/* --- Main loop over all FASTA entries in the input file ---- */

    sGetLoc(myargs[ARG_QUERYLOC].strvalue, &from, &to);

//	if (myargs[ARG_QUERYLOC].strvalue) {       
//		CharPtr delimiters = " ,;";
//		CharPtr location;
//		location = myargs[ARG_QUERYLOC].strvalue;
//		from = atoi(StringTokMT(location, delimiters, &location)) - 1;
//		to = atoi(location) - 1;
//		from = MAX(from, 0);
//	}
	/* if results will be generated a ByteStore needs to
	 * be allocated to store them */
	if( mode == SEARCH_MODE ){
		bs_ptr = BSNew(init_bsp_size);
		/* Put the result ByteStore where the c++ code can find it */
		ValNodeAddPointer( &result_bsp_list, 0, (Nlm_VoidPtr)bs_ptr );
	}

	/* loop over all requested queries, incrementing by the number 
	* of queries concatenated in an individual BLAST search */
	for( ; queryI <= last_query; queryI += query_increment ){
		int fake_bsp_rval = 0;
		meta_MPE_Log_event(fakebioseq_start,0,"get fakebioseq start");
		/*if( debug )fprintf( stderr, "Getting fake bioseq\n" );*/
		fake_bsp_rval = getFakeBioseq( TRUE, queryI );
		options->current_queryI = queryI;
		meta_MPE_Log_event(fakebioseq_end,0,"get fakebioseq end");
		if( fake_bsp_rval > 0 )
			return fake_bsp_rval;
		else if( fake_bsp_rval == -1 )
			break;

		/*if( debug )fprintf( stderr, "Performing BLAST\n" );*/
		fake_bsp_rval = performBLAST();
		if( fake_bsp_rval != 0 )
			return fake_bsp_rval;
		/*if( debug )fprintf( stderr, "Processing return data\n" );*/
		processReturnData();
		if( mode == SEARCH_MODE && seqalign != NULL ){
			/* this query had results, add its identifier to the list... */
			ValNodeAddInt( &result_id_list, 0, queryI );
			/* createBinASNResults will add all ASN.1 SeqAlign data to the asn_results_list
			* and standard SeqAlign data to the seqalign_results_list */
			/*if( debug )fprintf( stderr, "Creating ASN.1 results\n" );*/
			MpiBlastEnableBioseqFetch();
			createBinASNResults( seqalign, bs_ptr, queryI );
			MpiBlastDisableBioseqFetch();
		}else if( mode == SEARCH_MODE ){
			/* need to write search statistics even if no hits were found */
			AsnIoBSPtr aibp;
			if( use_binary_asn )
				aibp = AsnIoBSOpen("wb", bs_ptr);
			else
				aibp = AsnIoBSOpen("w", bs_ptr);
			writeSearchStatistics( aibp );
			AsnIoBSClose( aibp );
		}
		/*if( debug )fprintf( stderr, "Cleaning up the mess\n" );*/
		cleanupTheMess();

	} /* while(TRUE)  - main loop of the program over all FASTA entries */
	return 0;

}

/**
* run the blast algorithm in one of several modes:
*  0	Standard BLAST search as performed on a worker node
*  1	Collection of search space and effective HSP length
*/
Int2 runBLASTPIO( int mode, int first_query, int last_query, int frag_id ) {
	/* 
	* Use a 500 kb buffer for the initial byte store. Byte stores
	* grown dynamically and there may be some opportunity for optimizing this
	* initial size.
	*/
	const unsigned init_bsp_size = 512000;
	Uint1 tmp_align_view;
 	ByteStorePtr bs_ptr;

	int queryI = first_query;
	int query_increment = num_queries <= 0 ? 1 : num_queries;
	cur_sep = sep_list;

	/* used for profiling */
	double search_track_time = 0;
	double curr_search_time = 0;
	double result_track_time = 0;
	double curr_result_time = 0;

    SeqAlignPtr sap_filter, sap_cur;
    SeqIdPtr subject_sip;
    BioseqPtr tmpbsp;

	search_track_time = meta_MPI_Wtime();

	if( mode == COLLECT_STATS_MODE ){
		options->calculate_statistics_and_exit = TRUE;
	}else
		options->calculate_statistics_and_exit = FALSE;

	/* --- Main loop over all FASTA entries in the input file ---- */

	if (myargs[ARG_QUERYLOC].strvalue) {       
		CharPtr delimiters = " ,;";
		CharPtr location;
		location = myargs[ARG_QUERYLOC].strvalue;
		from = atoi(StringTokMT(location, delimiters, &location)) - 1;
		to = atoi(location) - 1;
		from = MAX(from, 0);
	}
	/* if results will be generated a ByteStore needs to
	 * be allocated to store them */
 	if( mode == SEARCH_MODE){
 		bs_ptr = BSNew(init_bsp_size);
		/* bs_ptr will be freed at QueryOutput::AddOutputRecordsToList()" */
 	}

	/* loop over all requested queries, incrementing by the number 
	* of queries concatenated in an individual BLAST search */
//	for( ; queryI <= last_query; queryI += query_increment ){
	for( ; queryI <= last_query; queryI += 1 ){
		int fake_bsp_rval = 0;
		meta_MPE_Log_event(fakebioseq_start,0,"get fakebioseq start");
		/*if( debug )fprintf( stderr, "Getting fake bioseq\n" );*/
		fake_bsp_rval = getFakeBioseq( TRUE, queryI );
		options->current_queryI = queryI;
		meta_MPE_Log_event(fakebioseq_end,0,"get fakebioseq end");
		if( fake_bsp_rval > 0 )
			return fake_bsp_rval;
		else if( fake_bsp_rval == -1 )
			break;

		if( mode == SEARCH_MODE){
			tmp_align_view = align_view;
			align_view = 11;	/* speed up search ?? */
		}
		
		fake_bsp_rval = performBLAST();
		
		if( mode == SEARCH_MODE){
			align_view = tmp_align_view;
		}
		
		if( fake_bsp_rval != 0 )
			return fake_bsp_rval;
		/*if( debug )fprintf( stderr, "Processing return data\n" );*/
		processReturnData();

		if( mode == SEARCH_MODE){
			/* write stats data here */
			AsnIoBSPtr aibp;
			if( use_binary_asn )
				aibp = AsnIoBSOpen("wb", bs_ptr);
			else
				aibp = AsnIoBSOpen("w", bs_ptr);
			writeSearchStatistics( aibp );
			AsnIoBSClose( aibp );
		}

		result_track_time = meta_MPI_Wtime();	

//        if(align_view == 7 || align_view == 8 || align_view == 9) { /* filter unknown subject sequences, may have little memory leak in the rare case*/
//            sap_filter = NULL;
//
//            while (seqalign) {
//                subject_sip = TxGetSubjectIdFromSeqAlign(seqalign);
//                tmpbsp = BioseqLockById(subject_sip);
//
//                if (tmpbsp) { 
//                    if(sap_filter == NULL) {
//                        sap_filter = seqalign;
//                        sap_cur = sap_filter;
//                    } else {
//                        sap_cur->next = seqalign;
//                        sap_cur = sap_cur->next;
//                    }
//                }
//                seqalign = seqalign->next;
//            }
//
//            seqalign = sap_filter;
//        }

		if( mode == SEARCH_MODE && seqalign != NULL ){
			/* this query had results, add its identifier to the list... */
			// ValNodeAddInt( &result_id_list, 0, queryI );
			
			MpiBlastEnableBioseqFetch();
			/* print and cache local results */

			outputResultsPIO(OUTPUT_MODE, seqalign, queryI, frag_id, bs_ptr);
		
			//createBinASNResults( seqalign, bs_ptr, queryI );
			MpiBlastDisableBioseqFetch();
			
		}else if( mode == SEARCH_MODE ){
			/* need to write search statistics even if no hits were found */
/*			AsnIoBSPtr aibp;
			if( use_binary_asn )
				aibp = AsnIoBSOpen("wb", bs_ptr);
			else
				aibp = AsnIoBSOpen("w", bs_ptr);
			writeSearchStatistics( aibp );
			AsnIoBSClose( aibp );*/
			
			if (mask_loc) {
				mask_loc_start = mask_loc;
			}	
			else
			{	/* Could have become non-NUll for last query. */
				mask_loc_start = NULL;
			}
			cProcessLocalOutputs(queryI, frag_id, NULL, NULL, bs_ptr);
		} else if( mode == COLLECT_STATS_MODE ){
			outputResultsPIO(COLLECT_INFO_MODE, NULL, queryI, -1, NULL); 
		}

		curr_result_time = meta_MPI_Wtime() - result_track_time;
		
		/*if( debug )fprintf( stderr, "Cleaning up the mess\n" );*/
		cleanupTheMess();

	} /* while(TRUE)  - main loop of the program over all FASTA entries */

	curr_search_time = meta_MPI_Wtime() - search_track_time;

	search_time += curr_search_time - curr_result_time;

	return 0;

}

/*
* Printing out html header here 
*/
void outputHtmlHeader(){
	if(align_view < 7) {
		if (html) {
			fprintf(outfp, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
			fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
				"VLINK=\"#660099\" ALINK=\"#660099\">\n");
			fprintf(outfp, "<PRE>\n");
		}
	} else if (align_view == 7) { 
		xml_aip = AsnIoOpen(blast_outputfile, "wx");
	}
#ifndef	BLAST_CS_API
	if(align_view >= 7 && myargs[ARG_NUMQUERIES].intvalue > 1)
	{
		ErrPostEx( SEV_FATAL, 1, 0,
			"blast: Query concatenation is currently not supported with -m > 7");
	}
#endif
}

void outputHtmlHeaderPIO(char* buffer){
	char *ptr = buffer;
	int print_len;
	
	if(align_view < 7) {
		if (html) {
			sprintf(ptr, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
			print_len = strlen(ptr);
			ptr += print_len;
			
			sprintf(ptr, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
				"VLINK=\"#660099\" ALINK=\"#660099\">\n");
			print_len = strlen(ptr);
			ptr += print_len;
			
			sprintf(ptr, "<PRE>\n");
			print_len = strlen(ptr);
			ptr += print_len;
		}
	} 
	
	return; 
}

void MpiBlastPrintReference( int html, Int4	line_length, FILE* outfp ) {
	if (outfp == NULL)
		return;

	if (html) {
		fprintf( outfp, "<b><a href=\"http://mpiblast.lanl.gov\">Reference</a>:</b><br>" );
		fprintf( outfp, "Aaron E. Darling, Lucas Carey, and Wu-chun Feng,<br> " );
		fprintf( outfp, "\"The design, implementation, and evaluation of mpiBLAST\".<br>" );
		fprintf( outfp, "In Proceedings of <a href=\"http://www.linuxclustersinstitute.org\">4th International Conference on Linux Clusters: The HPC Revolution 2003</a>, <br>June 24-26 2003, San Jose, CA<br>" );
		fprintf( outfp, "<br>Heshan Lin, Xiaosong Ma, Praveen Chandramohan, Al Geist, and Nagiza Samatova,<br> " );
		fprintf( outfp, "\"Efficient Data Access for Parallel BLAST\".<br>" );
		fprintf( outfp, "In Proceedings of <a href=\"http://www.ipdps.org/\">19th International Parallel & Distributed Processing Symposium 2005</a>, <br>April 3-8 2005, Denver, CO<br>" );
	} else {
		fprintf( outfp, "\nReference: \nAaron E. Darling, Lucas Carey, and Wu-chun Feng,\n" );
		fprintf( outfp, "\"The design, implementation, and evaluation of mpiBLAST\".\n" );
		fprintf( outfp, "In Proceedings of 4th International Conference on Linux Clusters: The HPC Revolution 2003, \nJune 24-26 2003, San Jose, CA\n\n" );
		
		fprintf( outfp, "Heshan Lin, Xiaosong Ma, Praveen Chandramohan, Al Geist, and Nagiza Samatova,\n" );
		fprintf( outfp, "\"Efficient Data Access for Parallel BLAST\".\n" );
		fprintf( outfp, "In Proceedings of 19th International Parallel & Distributed Processing Symposium 2005, \nApril 3-8 2005, Denver, CO\n\n" );
	}
}

/** a list of bioseqs that were added to the SeqMgr index */
ValNodePtr indexed_seqs = NULL;
/** 
 * Adds all subject bioseqs in a list of SeqAligns to the SeqMgr's index 
 *  needed because some NCBI output functions don't use BioseqLockById() but
 *  instead call BioseqFindCore()
 */
void indexSubjectBioseqs( SeqAlignPtr sap_p )
{
	SeqAlignPtr cur_sap = sap_p;
	while( cur_sap ){
		SeqIdPtr bioseq_id = TxGetSubjectIdFromSeqAlign( cur_sap );
		SeqMgrPtr smp = SeqMgrWriteLock();
		BioseqPtr subject_bsp = smp->bsfetch( bioseq_id, 0 );
		SeqMgrAddToBioseqIndex( subject_bsp );
		SeqMgrUnlock();
		indexed_seqs = ValNodeAddPointer( &indexed_seqs, 0, subject_bsp );
		cur_sap = cur_sap->next;
	}
}
/** remove previously indexed bioseqs from the index */
void unindexSubjectBioseqs()
{
	ValNodePtr cur_seq = indexed_seqs;
	while( cur_seq ){
		SeqMgrWriteLock();
		SeqMgrDeleteFromBioseqIndex( (BioseqPtr)(cur_seq->data.ptrvalue) );
		SeqMgrUnlock();
		cur_seq = cur_seq->next;
	}
	indexed_seqs = ValNodeFree( indexed_seqs );
}

/**
* output is done here
*/
Int2 outputResults(	SeqAlignPtr	sappy, int queryI )
{
	int fake_bsp_rval = 0;

	/* NOTE: if this function is called more times than the number of queries
	then it will re-cycle through the query list!
	*/
	if( cur_sep == NULL )
		cur_sep = sep_list;

	/* bioseqs are now pre-loaded here */
	fake_bsp_rval = getFakeBioseq( TRUE, queryI );

	if( fake_bsp_rval > 0 )
		return fake_bsp_rval;
	else if( fake_bsp_rval == -1 )
		return 0;

	/* do a 'dry run' to get the BLAST statistic data structures and db info */
	options->calculate_statistics_and_exit = TRUE;
	memcpy( options->stats, stats_array[queryI], sizeof( options->stats ) );
	fake_bsp_rval = performBLAST();
	if( fake_bsp_rval != 0 )
		return fake_bsp_rval;
	processReturnData();

	seqalign = sappy;
	global_fp = outfp;

	if(align_view < 7) {
		init_buff_ex(90);
		BlastPrintVersionInfo(blast_program, html, outfp);
		fprintf(outfp, "\n");
		if(strict_output_conformance)
			BlastPrintReference(html, 90, outfp);
		else
			MpiBlastPrintReference(html, 90, outfp);
		fprintf(outfp, "\n");
		if (!options->is_megablast_search) {
			/* KM added loop here for concat case */
			num_iters = (num_queries>0) ? num_queries : 1;
			for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {
				curr_bsp = (num_queries>0) ? *(fake_bsp_arr + bsp_iter) : query_bsp;
				AcknowledgeBlastQuery(curr_bsp, 70, outfp, believe_query, html);
			}
		}

		/* Here we first check, that database do no exists */

		if(!PrintDbInformation(blast_database, !db_is_na, 70, outfp, html))
			return 1;
		free_buff();
		if (options->is_ooframe)
			ErrPostEx(SEV_WARNING, 1, 0, "Out-of-frame option selected, Expect values are only approximate and calculated not assuming out-of-frame alignments");
	}
#ifdef OS_UNIX
	/*		  if(align_view	< 7) {
	fprintf(global_fp, "%s", "Searching");
	}
	*/
#endif		

#ifndef	BLAST_CS_API
	/*
	The function call to ReadDBBioseqFetchEnable has been moved into 
	mpiblast.cpp (in the outer loop).
	ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
	*/
#endif

	ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

	tmp_slp = slp;
	if (slp)
		query_bsp = NULL;

	if (getenv("POST_BLAST_CLUSTER_HITS") != NULL)
		BlastClusterHitsFromSeqAlign(seqalign, blast_program, blast_database, 
		options, 0.9, 1.6, 0.5, TRUE);

	if (mask_loc) {
		mask_loc_start = mask_loc;
	}	
	else
	{	/* Could have become non-NUll for last query. */
		mask_loc_start = NULL;
	}

	/* Print header in any case */
	if (align_view == 9) {
		PrintTabularOutputHeader(blast_database, query_bsp, slp, 
			blast_program, 0, believe_query,
			global_fp);
	}

	if (seqalign) {
		if (num_queries > 0) { /* AM: Support for query multiplexing. */
			sap_array = mult_queries->sap_array_data->sap_array;
		}

		if (align_view == 8 || align_view == 9) {
			/* --KM	need to	put	a loop around this.	seqaligns already broken up
			note the method for looping if num_aligns > 0 - reuse this method everywhere */
			num_iters = (num_queries>0) ? num_queries : 1;
			for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
				curr_seqalign = (num_queries>0) ? *(sap_array + sap_iter) : seqalign;
				BlastPrintTabularResults(curr_seqalign, query_bsp, slp, 
					number_of_alignments, blast_program, 
					!options->gapped_calculation, options->is_ooframe,
					believe_query, 0, 0, global_fp, NULL, (align_view == 9));

				/*				seqalign = NULL;
				*/				SeqAlignSetFree(curr_seqalign);
			}
		} else {
			while (seqalign) {

				if (!options->is_megablast_search){
					next_seqalign = NULL;
				} else {
					SeqIdPtr sip, next_sip = NULL;

					sap = seqalign;
					sip = TxGetQueryIdFromSeqAlign(seqalign);

					while (sap != NULL) { 
						if (sap->next != NULL) {
							next_sip = TxGetQueryIdFromSeqAlign(sap->next);

							if (SeqIdComp(sip, next_sip) != SIC_YES) {
								next_seqalign = sap->next;
								sap->next = NULL;
							}
						} else{
							next_seqalign = NULL;
						}
						sap = sap->next;
				 }

					while (tmp_slp && SeqIdComp(sip, SeqLocId(tmp_slp)) != SIC_YES)
						tmp_slp = tmp_slp->next;
					if (tmp_slp == NULL) /* Should never happen */
						break;
					/* Separate the mask locations list for this query */
					if (!mask_loc && next_mask_loc) {
						mask_loc = next_mask_loc;
						next_mask_loc = NULL;
				 }
					if (mask_loc) {
						if (next_mask_loc) {
							mask_loc->next = next_mask_loc;
							mask_loc = next_mask_loc;
						}
						mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
						next_mask_loc = mask_loc;
						while (SeqIdComp(SeqLocId(mask_slp), sip) != SIC_YES) {
							mask_loc = mask_loc->next;
							if (!mask_loc)
								break;
							mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
						}
						if (mask_loc) {
							next_mask_loc = mask_loc->next;
							mask_loc->next = NULL;
						}
				 }

					if (align_view < 7) {
						bsp = BioseqLockById(SeqLocId(tmp_slp));
						init_buff_ex(85);
						fprintf(outfp, "\n");
						AcknowledgeBlastQuery(bsp, 70, outfp, believe_query, html);
						free_buff();
						BioseqUnlock(bsp);
				 }
				}

				if(align_view == 7 && !options->is_ooframe) {
					if (options->is_megablast_search) {
						bsp = BioseqLockById(SeqLocId(tmp_slp));
						BXMLPrintOutput(xml_aip, seqalign, 
							options, blast_program, blast_database, 
							bsp, other_returns, 0, NULL, mask_loc);
						BioseqUnlock(bsp);
						AsnIoReset(xml_aip);
						seqalign = SeqAlignSetFree(seqalign);

					} else {
						num_iters = (num_queries>0) ? num_queries : 1;
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
							BXMLPrintOutput(xml_aip, curr_seqalign,
								options, blast_program, blast_database, 
								fake_bsp, other_returns, 0, NULL, mask_loc);

							AsnIoReset(xml_aip);
							curr_seqalign = SeqAlignSetFree(curr_seqalign);
							if( num_queries == 0 )
								seqalign = curr_seqalign;

						}  /* for loop over sap-array (concat) */
				 }
				} else {

					/* create the array of SeqAnnotPtrs, if necessary */

					num_iters = (num_queries > 0) ? num_queries : 1;
					for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
						curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
						if ( (num_queries > 0) && (sap_iter == 0) ) {
							seq_annot_arr = (SeqAnnotPtrArray) MemNew(sizeof(SeqAnnotPtr)*num_queries);
						}
						seqannot = SeqAnnotNew();
						seqannot->type = 2;
						AddAlignInfoToSeqAnnot(seqannot, align_type);
						seqannot->data = curr_seqalign;

						if (num_queries > 0) {
							*(seq_annot_arr + sap_iter) = seqannot;
						}
					} /* make seqannots over the sap_iters from concat, or the single seqalign */


					if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
						ObjMgrSetHold();
						/* print deflines */
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							intercept_output_pairwise(INT_DES_BEGIN);
							
							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;

							init_buff_ex(85);

							PrintDefLinesFromSeqAlignEx2(curr_seqalign, 80, outfp,
								print_options, FIRST_PASS, NULL,
								number_of_descriptions, NULL, NULL);

							free_buff();
							intercept_output_pairwise(INT_DES_END);
						} /* print deflines, looped if concat */
						
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							/* AM: Query concatenation. */
							if( mult_queries && mask_loc )
							{
								orig_mask_loc = mask_loc;

								if( !mask_loc->data.ptrvalue ) mask_loc = NULL;
							}

							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
							curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;
							/* mpiBLAST: some NCBI functions do not use the 
							 * SeqMgr's bioseq fetch function (e.g. MuskSeqIdWrite)
							 * In order to allow such functions to look up our bioseqs
							 * we need to preload them into the SeqMgr
							 * This should probably get reported as a bug to NCBI...
							 */
							prune = BlastPruneHitsFromSeqAlign(curr_seqalign,
								number_of_alignments, NULL);

							curr_seqannot->data = prune->sap;

							if( align_view >= 1 && align_view <= 6 )
								indexSubjectBioseqs(curr_seqannot->data);
							
							intercept_output_pairwise(INT_ALN_BEGIN);
							
							if(options->is_ooframe) {
								OOFShowBlastAlignment(curr_seqalign, /*mask*/ NULL,
									outfp, align_options, txmatrix);
							} else {
								if (align_view != 0)
									ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
									align_options, txmatrix, mask_loc, NULL);
								else
									ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
									align_options, txmatrix, mask_loc,
									FormatScoreFunc);
							}
							
							intercept_output_pairwise(INT_ALN_END);

							if( align_view >= 1 && align_view <= 6 )
								unindexSubjectBioseqs(curr_seqannot->data);
							if(aip == NULL){
								curr_seqannot->data = NULL;  
							}
							else{
								/* 
								* Restore the SeqAlign data that was temporarily replaced
								* by the prune->sap data. The SeqAlign data will be needed
								* since the user is writing a ASN.1 of SeqAnnot data.
								*/
								curr_seqannot->data = curr_seqalign;  
							}

							prune = BlastPruneSapStructDestruct(prune);

							/* AM: Query concatenation. */
							if( mult_queries && orig_mask_loc ) 
							{
								mask_loc = orig_mask_loc;
								mask_loc = mask_loc->next;
							}
						} /* show text align, loop over seqalign/seqannots for concat */

						ObjMgrClearHold();

					} /* if outfp */


					if (aip) {
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							curr_seqannot = 
								(num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;

							SeqAnnotAsnWrite((SeqAnnotPtr) curr_seqannot, aip, NULL);
							AsnIoReset(aip);
						}
					} /* if aip */

					for (sap_iter=0; sap_iter < num_queries; sap_iter++) {
						/* upper bound is num_queries, take care not to do this unless concat */
						*(seq_annot_arr + sap_iter) = SeqAnnotFree(*(seq_annot_arr + sap_iter));
					}
					/*--KM free seqalign array and all seqaligns?? */

				} /* end of else (not XML Printing) */

				if (options->is_megablast_search)
					tmp_slp = tmp_slp->next;

				/* --KM watch for memory leaks */
				if (seqannot && num_queries == 0)
					seqannot = SeqAnnotFree(seqannot);

				if( aip == NULL ){
					seqalign = SeqAlignSetFree(seqalign);
				}

				seqalign = next_seqalign;
			} /* End of loop on all seqaligns */
			if (mask_loc && next_mask_loc)
				mask_loc->next = next_mask_loc;

		} /* end of align_view not tabular case */
	} else {         /* seqalign is NULL */
		if(align_view == 7 && !options->is_ooframe) {
			BlastErrorMsgPtr error_msg;
			CharPtr message;

			if (error_returns == NULL) {
				message = "No hits found";
			} else {
				error_msg = error_returns->data.ptrvalue;
				message = error_msg->msg;
			}
			if (options->is_megablast_search) {
				bsp = BioseqLockById(SeqLocId(tmp_slp));
				BXMLPrintOutput(xml_aip, seqalign, 
					options, blast_program, blast_database, 
					bsp, other_returns, 0, NULL, mask_loc);
				BioseqUnlock(bsp);
			} else {
				BXMLPrintOutput(xml_aip, NULL, 
					options, blast_program, blast_database, 
					fake_bsp, other_returns, 0, message, mask_loc);
			}
			AsnIoReset(xml_aip);
		} else if (align_view < 8) {
			fprintf(outfp, "\n\n ***** No hits found ******\n\n");
		}
		if (error_returns != NULL) {
			for (vnp = error_returns; vnp; vnp = vnp->next) {
				BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
			}
			ValNodeFree(error_returns);
		}
	}
	
	if(html) {
		fprintf(outfp, "<PRE>\n");
	}

	init_buff_ex(85);
	dbinfo_head = dbinfo;

	if(align_view < 7 && done) {
		while (dbinfo) {
			PrintDbReport(dbinfo, 70, outfp);
			dbinfo = dbinfo->next;
		}
	}

	if (ka_params) {
		if(align_view < 7 && done) {
			PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
		}
	}

	if (ka_params_gap) {
		if(align_view < 7 && done) {
			PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
		}
	}

	if(align_view < 7 && done) {
		PrintTildeSepLines(params_buffer, 70, outfp);
	}

	free_buff();
	if (html)
		fprintf(outfp, "</PRE>\n<P><HR><BR>\n<PRE>");
	/*
	if (done) 
	sep = SeqEntryFree(sep);
	*/
	/* mpiBLAST: make sure the entire query results get written */
	if( outfp ){
		fflush( outfp );
	}
	if( aip ){
		AsnIoFlush( aip );
	}
	cleanupTheMess();

	return 0;
}

/**
* output is done here
*/
Int2 outputResultsPIO(int mode, SeqAlignPtr sappy, int queryI, int frag_id, ByteStorePtr bs_ptr )
{
	int fake_bsp_rval = 0;
	ValNodePtr orplist;
	int intercept_mode;
	int tag = -1;

	double result_track_time = 0;
	double curr_result_time = 0;
	double process_track_time = 0;
	double curr_process_time = 0;

	result_track_time = meta_MPI_Wtime();

	if(mode == COLLECT_INFO_MODE) {
		/* !!! in OUTPUT_MODE fake_bsp has been fetched in runBLASTPIO, no need to fetch again */
		/* NOTE: if this function is called more times than the number of queries
		then it will re-cycle through the query list!
		*/
		if( cur_sep == NULL )
			cur_sep = sep_list;
	
		/* bioseqs are now pre-loaded here */
		fake_bsp_rval = getFakeBioseq( TRUE, queryI );
	
		if( fake_bsp_rval > 0 )
			return fake_bsp_rval;
		else if( fake_bsp_rval == -1 )
			return 0;
	
		if(!(align_view == 7 && use_brief_report)) {
			/* do a 'dry run' to get the BLAST statistic data structures and db info */
			options->calculate_statistics_and_exit = TRUE;
			memcpy( options->stats, stats_array[queryI], sizeof( options->stats ) );
			fake_bsp_rval = performBLAST();
			if( fake_bsp_rval != 0 )
				return fake_bsp_rval;
			processReturnData();
		}
	}
	seqalign = sappy;
	global_fp = outfp;

	if(mode == OUTPUT_MODE) {
		intercept_mode =  INT_INFO_RECORDS;
	} else if (mode == COLLECT_INFO_MODE) {
		intercept_mode = INT_INFO;
	} else {
		return -1;
	}
	
	if(outfp) {
		start_intercept(outfp, align_view, html, intercept_mode);
	} else if (aip) {
		start_intercept(aip->fp, align_view, html, intercept_mode);
	} else if (xml_aip) {
		start_intercept(xml_aip->fp, align_view, html, intercept_mode);
	} else {
		fprintf(stderr, "outputResultsPIO(): output stream has not been opened\n");
		return 0;
	}

	if(align_view==10 || align_view==11) {
		intercept_output_asn(INT_HDR_BEGIN);
	}
	
	if(mode == OUTPUT_MODE) {
		create_output_records(seqalign, &orplist, options->hitlist_size);
		set_results_info(orplist, 0);
	}

	if(outfp) {
		switch (align_view) {
		case 0:
			intercept_output_pairwise(INT_HDR_BEGIN);
			break;
		case 8:
		case 9:
			intercept_output_tabular(INT_HDR_BEGIN);
			break;
		}
	} 
	
	if(align_view < 7) {
		init_buff_ex(90);
		BlastPrintVersionInfo(blast_program, html, outfp);
		fprintf(outfp, "\n");
		if(strict_output_conformance)
			BlastPrintReference(html, 90, outfp);
		else
			MpiBlastPrintReference(html, 90, outfp);
		fprintf(outfp, "\n");
		if (!options->is_megablast_search) {
			/* KM added loop here for concat case */
			num_iters = (num_queries>0) ? num_queries : 1;
			for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {
				curr_bsp = (num_queries>0) ? *(fake_bsp_arr + bsp_iter) : query_bsp;
				AcknowledgeBlastQuery(curr_bsp, 70, outfp, believe_query, html);
			}
		}

		/* Here we first check, that database do no exists */

		if(!PrintDbInformation(blast_database, !db_is_na, 70, outfp, html))
			return 1;
		free_buff();
		if (options->is_ooframe)
			ErrPostEx(SEV_WARNING, 1, 0, "Out-of-frame option selected, Expect values are only approximate and calculated not assuming out-of-frame alignments");
	}
#ifdef OS_UNIX
	/*		  if(align_view	< 7) {
	fprintf(global_fp, "%s", "Searching");
	}
	*/
#endif		

#ifndef	BLAST_CS_API
	/*
	The function call to ReadDBBioseqFetchEnable has been moved into 
	mpiblast.cpp (in the outer loop).
	ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
	*/
#endif

	ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

	tmp_slp = slp;
	if (slp)
		query_bsp = NULL;

	if (getenv("POST_BLAST_CLUSTER_HITS") != NULL)
		BlastClusterHitsFromSeqAlign(seqalign, blast_program, blast_database, 
		options, 0.9, 1.6, 0.5, TRUE);

	if (mask_loc) {
		mask_loc_start = mask_loc;
	}	
	else
	{	/* Could have become non-NUll for last query. */
		mask_loc_start = NULL;
	}

	/* Print header in any case */
	if (align_view == 9) {
		PrintTabularOutputHeader(blast_database, query_bsp, slp, 
			blast_program, 0, believe_query,
			global_fp);
	}

	if(outfp) {
		switch (align_view) {
		case 0:
			intercept_output_pairwise(INT_HDR_END);
			break;
		case 8:
		case 9:
			intercept_output_tabular(INT_HDR_END);
			break;
		}
	} 
	
	if (seqalign) {
		if (num_queries > 0) { /* AM: Support for query multiplexing. */
			sap_array = mult_queries->sap_array_data->sap_array;
		}

		if (align_view == 8 || align_view == 9) {
			
			/* --KM	need to	put	a loop around this.	seqaligns already broken up
			note the method for looping if num_aligns > 0 - reuse this method everywhere */
			num_iters = (num_queries>0) ? num_queries : 1;
			for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
				intercept_output_tabular(INT_ALN_BEGIN);

				curr_seqalign = (num_queries>0) ? *(sap_array + sap_iter) : seqalign;
				BlastPrintTabularResults(curr_seqalign, query_bsp, slp, 
					number_of_alignments, blast_program, 
					!options->gapped_calculation, options->is_ooframe,
					believe_query, 0, 0, global_fp, NULL, (align_view == 9));

				/*				seqalign = NULL;
				*/				SeqAlignSetFree(curr_seqalign);
			
				intercept_output_tabular(INT_ALN_END);
			}
			
		} else {
			while (seqalign) {

				if (!options->is_megablast_search){
					next_seqalign = NULL;
				} else {
					SeqIdPtr sip, next_sip = NULL;

					sap = seqalign;
					sip = TxGetQueryIdFromSeqAlign(seqalign);

					while (sap != NULL) { 
						if (sap->next != NULL) {
							next_sip = TxGetQueryIdFromSeqAlign(sap->next);

							if (SeqIdComp(sip, next_sip) != SIC_YES) {
								next_seqalign = sap->next;
								sap->next = NULL;
							}
						} else{
							next_seqalign = NULL;
						}
						sap = sap->next;
				 }

					while (tmp_slp && SeqIdComp(sip, SeqLocId(tmp_slp)) != SIC_YES)
						tmp_slp = tmp_slp->next;
					if (tmp_slp == NULL) /* Should never happen */
						break;
					/* Separate the mask locations list for this query */
					if (!mask_loc && next_mask_loc) {
						mask_loc = next_mask_loc;
						next_mask_loc = NULL;
				 }
					if (mask_loc) {
						if (next_mask_loc) {
							mask_loc->next = next_mask_loc;
							mask_loc = next_mask_loc;
						}
						mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
						next_mask_loc = mask_loc;
						while (SeqIdComp(SeqLocId(mask_slp), sip) != SIC_YES) {
							mask_loc = mask_loc->next;
							if (!mask_loc)
								break;
							mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
						}
						if (mask_loc) {
							next_mask_loc = mask_loc->next;
							mask_loc->next = NULL;
						}
				 }

					if (align_view < 7) {
						bsp = BioseqLockById(SeqLocId(tmp_slp));
						init_buff_ex(85);
						fprintf(outfp, "\n");
						AcknowledgeBlastQuery(bsp, 70, outfp, believe_query, html);
						free_buff();
						BioseqUnlock(bsp);
				 }
				}

				if(align_view == 7 && !options->is_ooframe) {
					if (options->is_megablast_search) {
						bsp = BioseqLockById(SeqLocId(tmp_slp));
						BXMLPrintOutput(xml_aip, seqalign, 
							options, blast_program, blast_database, 
							bsp, other_returns, 0, NULL, mask_loc);
						BioseqUnlock(bsp);
						AsnIoReset(xml_aip);
						seqalign = SeqAlignSetFree(seqalign);

					} else {
						num_iters = (num_queries>0) ? num_queries : 1;
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
							BXMLPrintOutput(xml_aip, curr_seqalign,
								options, blast_program, blast_database, 
								fake_bsp, other_returns, 0, NULL, mask_loc);

							AsnIoReset(xml_aip);
							curr_seqalign = SeqAlignSetFree(curr_seqalign);
							if( num_queries == 0 )
								seqalign = curr_seqalign;

						}  /* for loop over sap-array (concat) */
				 }
				} else {

					/* create the array of SeqAnnotPtrs, if necessary */

					num_iters = (num_queries > 0) ? num_queries : 1;
					for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
						curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
						if ( (num_queries > 0) && (sap_iter == 0) ) {
							seq_annot_arr = (SeqAnnotPtrArray) MemNew(sizeof(SeqAnnotPtr)*num_queries);
						}
						seqannot = SeqAnnotNew();
						seqannot->type = 2;
						AddAlignInfoToSeqAnnot(seqannot, align_type);
						seqannot->data = curr_seqalign;

						if (num_queries > 0) {
							*(seq_annot_arr + sap_iter) = seqannot;
						}
					} /* make seqannots over the sap_iters from concat, or the single seqalign */


					if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
						ObjMgrSetHold();
						
						/* print deflines */
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
 							intercept_output_pairwise(INT_DES_BEGIN);

							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;

							init_buff_ex(85);

							PrintDefLinesFromSeqAlignEx2(curr_seqalign, 80, outfp,
								print_options, FIRST_PASS, NULL,
								number_of_descriptions, NULL, NULL);

							free_buff();

 							intercept_output_pairwise(INT_DES_END);
						} /* print deflines, looped if concat */
						
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							/* AM: Query concatenation. */
							if( mult_queries && mask_loc )
							{
								orig_mask_loc = mask_loc;

								if( !mask_loc->data.ptrvalue ) mask_loc = NULL;
							}

							curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
							curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;
							/* mpiBLAST: some NCBI functions do not use the 
							 * SeqMgr's bioseq fetch function (e.g. MuskSeqIdWrite)
							 * In order to allow such functions to look up our bioseqs
							 * we need to preload them into the SeqMgr
							 * This should probably get reported as a bug to NCBI...
							 */
							prune = BlastPruneHitsFromSeqAlign(curr_seqalign,
								number_of_alignments, NULL);

							curr_seqannot->data = prune->sap;

							//indexSubjectBioseqs(curr_seqannot->data); not needed with pio
							
 							intercept_output_pairwise(INT_ALN_BEGIN);
							if(options->is_ooframe) {
								OOFShowBlastAlignment(curr_seqalign, /*mask*/ NULL,
									outfp, align_options, txmatrix);
							} else {
								if (align_view != 0)
									ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
									align_options, txmatrix, mask_loc, NULL);
								else
									ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
									align_options, txmatrix, mask_loc,
									FormatScoreFunc);
							}
 							intercept_output_pairwise(INT_ALN_END);

							//unindexSubjectBioseqs(curr_seqannot->data); not needed with pio
							if(aip == NULL){
								curr_seqannot->data = NULL;  
							}
							else{
								/* 
								* Restore the SeqAlign data that was temporarily replaced
								* by the prune->sap data. The SeqAlign data will be needed
								* since the user is writing a ASN.1 of SeqAnnot data.
								*/
								curr_seqannot->data = curr_seqalign;  
							}

							prune = BlastPruneSapStructDestruct(prune);

							/* AM: Query concatenation. */
							if( mult_queries && orig_mask_loc ) 
							{
								mask_loc = orig_mask_loc;
								mask_loc = mask_loc->next;
							}
						} /* show text align, loop over seqalign/seqannots for concat */

						ObjMgrClearHold();

					} /* if outfp */


					if (aip) {
						for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
							curr_seqannot = 
								(num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;

							SeqAnnotAsnWrite((SeqAnnotPtr) curr_seqannot, aip, NULL);
							AsnIoReset(aip);
						}
					} /* if aip */

					for (sap_iter=0; sap_iter < num_queries; sap_iter++) {
						/* upper bound is num_queries, take care not to do this unless concat */
						*(seq_annot_arr + sap_iter) = SeqAnnotFree(*(seq_annot_arr + sap_iter));
					}
					/*--KM free seqalign array and all seqaligns?? */

				} /* end of else (not XML Printing) */

				if (options->is_megablast_search)
					tmp_slp = tmp_slp->next;

				/* --KM watch for memory leaks */
				if (seqannot && num_queries == 0)
					seqannot = SeqAnnotFree(seqannot);

				if( aip == NULL ){
					seqalign = SeqAlignSetFree(seqalign);
				}

				seqalign = next_seqalign;
			} /* End of loop on all seqaligns */
			if (mask_loc && next_mask_loc)
				mask_loc->next = next_mask_loc;

		} /* end of align_view not tabular case */
	} else {         /* seqalign is NULL */
		if(align_view == 7 && !options->is_ooframe) {
			BlastErrorMsgPtr error_msg;
			CharPtr message;

			if (error_returns == NULL) {
				message = "No hits found";
			} else {
				error_msg = error_returns->data.ptrvalue;
				message = error_msg->msg;
			}
			if (options->is_megablast_search) {
				bsp = BioseqLockById(SeqLocId(tmp_slp));
				BXMLPrintOutput(xml_aip, seqalign, 
					options, blast_program, blast_database, 
					bsp, other_returns, 0, NULL, mask_loc);
				BioseqUnlock(bsp);
			} else {
				BXMLPrintOutput(xml_aip, NULL, 
					options, blast_program, blast_database, 
					fake_bsp, other_returns, 0, message, mask_loc);
			}
			AsnIoReset(xml_aip);
		} else if (align_view < 8) {
			fprintf(outfp, "\n\n ***** No hits found ******\n\n");
		}
		if (error_returns != NULL) {
			for (vnp = error_returns; vnp; vnp = vnp->next) {
				BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
			}
			ValNodeFree(error_returns);
		}
	}
	
	if(outfp) {
		if(align_view == 0) {
			intercept_output_pairwise(INT_STA_BEGIN);
		}
	}

	if(html) {
		fprintf(outfp, "<PRE>\n");
	}

	init_buff_ex(85);
	dbinfo_head = dbinfo;

	if(align_view < 7 && done) {
		while (dbinfo) {
			PrintDbReport(dbinfo, 70, outfp);
			dbinfo = dbinfo->next;
		}
	}

	if (ka_params) {
		if(align_view < 7 && done) {
			PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
		}
	}

	if (ka_params_gap) {
		if(align_view < 7 && done) {
			PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
		}
	}
	
	if(align_view < 7 && done) {
		PrintTildeSepLines(params_buffer, 70, outfp);
	}

	if(outfp) {
		if(align_view == 0) {
			intercept_output_pairwise(INT_STA_END);
		}
	}

	free_buff();
	
	if(align_view==10 || align_view==11) {
		AsnIoFlush(aip);
		intercept_output_asn(INT_FTR_END);
	}
	
	if(outfp) {
		switch (align_view) {
		case 0:
			intercept_output_pairwise(INT_FTR_BEGIN);
			break;
		case 8:
		case 9:
			intercept_output_tabular(INT_FTR_BEGIN);
			break;
		}
	}
	
	if (html)
		fprintf(outfp, "</PRE>\n<P><HR><BR>\n<PRE>");
	
	if(outfp) {
		switch (align_view) {
		case 0:
			intercept_output_pairwise(INT_FTR_END);
			break;
		case 8:
		case 9:
			intercept_output_tabular(INT_FTR_END);
			break;
		}
	}
	/*
	if (done) 
	sep = SeqEntryFree(sep);
	*/
	/* mpiBLAST: make sure the entire query results get written */
	if( outfp ){
		fflush( outfp );
	}
	if( aip ){
		AsnIoFlush( aip );
	}
	
	if(mode == COLLECT_INFO_MODE) {
		cleanupTheMess();
	}
	
	process_track_time = meta_MPI_Wtime();
	//print_output_records(stdout, orplist);
	end_intercept(queryI, frag_id, bs_ptr);
	curr_process_time = meta_MPI_Wtime() - process_track_time;
	
	curr_result_time = meta_MPI_Wtime() - result_track_time;

	result_time += curr_result_time - curr_process_time;

	return 0;
}

void outputHTMLfooter(){
	if(align_view < 7) {
		if (html) {
			fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
		}
	}
}

void outputHtmlFooterPIO(char *buffer){
	if(align_view < 7) {
		if (html) {
			sprintf(buffer, "</PRE>\n</BODY>\n</HTML>\n");
		}
	}
}

void cleanupBLAST(){
	aip = AsnIoClose(aip);

	if (align_view == 7)
		xml_aip = AsnIoClose(xml_aip);

	options = BLASTOptionDelete(options);
	FileClose(infp);
	
	if(outfp != NULL) {
		FileClose(outfp);
	}
}

void MpiBlastEnableBioseqFetch()
{
#ifndef BLAST_CS_API
	ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
#endif
}

void MpiBlastDisableBioseqFetch()
{
#ifndef BLAST_CS_API
	ReadDBBioseqFetchDisable();
#endif
}

int GetNumOfShowDes() {
	return myargs[ARG_DESCRIPTIONS].intvalue;	
}

int GetNumOfShowAlign() {
	return myargs[ARG_ALIGNMENTS].intvalue;	
}

Uint1 GetFormatAlignView() {
	return align_view;
}

Boolean GetFormatHtml() {
	return html;
}

int extractStats(int queryI, unsigned char* stats_buf, int stats_buf_size) {
	AsnIoMemPtr aimp;
	AsnTypePtr atp;
	DataVal dv;
	int statI;
	int num_des;
	int num_align;
	int max_show;

	if( use_binary_asn )
		aimp = AsnIoMemOpen( "rb", stats_buf, stats_buf_size);
	else
		aimp = AsnIoMemOpen( "r", stats_buf, stats_buf_size);

	if(aimp == NULL) {
		return -1;
	}
	
	for( statI = 0; statI < 12; statI++ ){
		// read "Stat-data" struct opening
		atp = AsnReadId( aimp->aip, NULL, &statset_type );
		AsnReadVal( aimp->aip, atp, &dv );
		// read each stat-data value
//			at.type = &at;
		atp = AsnReadId( aimp->aip, NULL, &statdata_type );
		AsnReadVal( aimp->aip, atp, &dv );
		stats_array[queryI][statI] += dv.bigintvalue;
		// read "Stat-data" struct close
		atp = AsnReadId( aimp->aip, NULL, &statset_type );
		AsnReadVal( aimp->aip, atp, &dv );
	}

	AsnIoMemClose(aimp);
	return 0;
}

int GetQueryAcknowledge(int query_id, char** query_ack) {
	int fake_bsp_rval = 0;
	int buf_size = 0;

	start_comm_intercept(outfp);

	fake_bsp_rval = getFakeBioseq( TRUE, query_id);

	// get query info
	if (!options->is_megablast_search) {
		num_iters = (num_queries>0) ? num_queries : 1;
		for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {
			curr_bsp = (num_queries>0) ? *(fake_bsp_arr + bsp_iter) : query_bsp;
			AcknowledgeBlastQuery(curr_bsp, 70, outfp, believe_query, html);
		}
	}

	buf_size = end_comm_intercept(query_ack);

	return buf_size;
}
