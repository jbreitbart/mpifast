#ifndef API_ACEREAD__H
#define API_ACEREAD__H

/*
 * $Id: aceread.h,v 1.12 2008/12/02 18:58:24 bollin Exp $
 *
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Colleen Bollin
 *
 */

#include <util/creaders/creaders_export.h>

#ifdef __cplusplus
extern "C" {
#endif

/* defines from ncbistd.h */
#ifndef FAR
#define FAR
#endif
#ifndef PASCAL
#define PASCAL
#endif
#ifndef EXPORT
#define EXPORT
#endif

#ifndef PASCAL
#define PASCAL
#endif
#ifndef EXPORT
#define EXPORT
#endif

#if defined (WIN32)
#    define ASSEMBLY_CALLBACK __stdcall
#else
#    define ASSEMBLY_CALLBACK
#endif

typedef struct gapinfo {
    int num_gaps;
    int *gap_offsets;
} SGapInfo, * TGapInfoPtr;

extern TGapInfoPtr GapInfoNew (void);
extern void GapInfoFree (TGapInfoPtr g);
extern TGapInfoPtr GapInfoFromSequenceString (char *seq_str, char *gap_chars);
extern void RemoveGapCharsFromSequenceString (char *seq_str, char *gap_chars);
extern int SeqPosFromTilingPos (int tiling_pos, TGapInfoPtr gap_info);
extern int TilingPosFromSeqPos (int seq_pos, TGapInfoPtr gap_info);

typedef struct SContigRead {
    char * read_id;
    int    ti;
    char * srr;
    char * read_seq;
    int    read_len;
    char   is_complement;
    int    cons_start;
    int    cons_stop;
    int    read_start;
    int    read_stop;
    int    read_assem_start;
    int    read_assem_stop;
    int    tiling_start;
    int    tiling_stop;
    TGapInfoPtr gaps;
    int    valid;
    int    local;
    char * tag; /* notes, comments, annotation for the read */
    /* quality scores - these are optional, used when recalculating consensus sequence */
    int  * qual_scores;
    int    num_qual_scores;
} SContigRead, * TContigReadPtr;

extern TContigReadPtr ContigReadNew (void);
extern void ContigReadFree (TContigReadPtr r);

typedef struct SConsensusReadAln {
    int numseg;
    int *cons_starts;
    int *read_starts;
    int *lens;
    char is_complement;
} SConsensusReadAln, * TConsensusReadAlnPtr;

extern TConsensusReadAlnPtr ConsensusReadAlnNew (int numseg);
extern TConsensusReadAlnPtr ConsensusReadAlnFree (TConsensusReadAlnPtr a);
extern TConsensusReadAlnPtr GetConsensusReadAln (char *consensus_seq, TContigReadPtr read);


typedef struct SBaseSeg {
    char * read_id;
    int    cons_start;
    int    cons_stop;
} SBaseSeg, * TBaseSegPtr;

extern TBaseSegPtr BaseSegNew (void);
extern void BaseSegFree (TBaseSegPtr b);

typedef struct SContig {
    char  * consensus_id;
    char  * consensus_seq;
    int     consensus_assem_len;
    int     consensus_seq_len;
    char    is_complement;
    int     num_qual_scores;
    int   * qual_scores;
    TGapInfoPtr gaps;
    int     num_reads;
    TContigReadPtr * reads;
    int     num_base_segs;
    TBaseSegPtr *base_segs;
    char  * tag; /* notes, comments, annotation for the contig */
} SContig, * TContigPtr;

extern TContigPtr ContigNew (void);
extern void ContigFree (TContigPtr c);
   
typedef struct SACEFile {
    int        num_contigs;
    TContigPtr * contigs;
} SACEFile, * TACEFilePtr;

extern NCBI_CREADERS_EXPORT TACEFilePtr ACEFileNew (void);
extern NCBI_CREADERS_EXPORT void ACEFileFree (TACEFilePtr afp);

extern NCBI_CREADERS_EXPORT TACEFilePtr ReadACEFile (
  FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
  void *               fileuserdata,  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
  char                 make_qual_scores, /* false if ignoring 
                                          * known-bad qual scores
                                          */
  char *               has_errors        /* starts false if errors have already been reported
                                          * set to true if errors are encountered
                                          */
);


extern NCBI_CREADERS_EXPORT TACEFilePtr ReadMAQFile (
 FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
);


extern void WriteACEFile (FILE *fp, TACEFilePtr afp);

extern TAlignmentFilePtr AlignmentFileFromContig (TContigPtr contig);

extern char * TraceArchiveGapStringFromACESequence (char *seq_str);

extern TContigReadPtr 
ReadContigFromString 
(char  *str,
 char **consensus_id,
 int    id_col,
 int    seq_col, 
 int    contig_id_col,
 int    strand_col,
 int    start_col,
 int    interpret_n_col);
extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromMAQString (char *str, char **consensus_id);
extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandMostCompressed (char *str, char **consensus_id);
extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandSanger (char *str, char **consensus_id);
extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandStandalone (char *str, char **consensus_id);

typedef TContigReadPtr (ASSEMBLY_CALLBACK *FReadFromStringFunction) (char *str, char **consensus_id);
extern TACEFilePtr ReadAssemblyFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata,  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
 FReadFromStringFunction makeread_func); /* function to transform a string into a read */

extern TACEFilePtr ReadMAQFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata);  /* data to be passed back each time
                                       * readfunc is invoked
                                       */

extern TACEFilePtr ReadElandStandaloneFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata);  /* data to be passed back each time
                                       * readfunc is invoked
                                       */


extern void 
WriteTraceAssemblyFromAceFile 
(TACEFilePtr afp,
 char      * subref,
 char      * center_name, 
 int         taxid,
 char      * description,
 FILE      * fp);

extern void
WriteTraceAssemblyHeader
(char * assembly_type,
 char * subref,
 char * center_name,
 int    taxid,
 char * description,
 char * assembly,
 int    num_contigs,
 unsigned int    num_conbases,
 int    num_reads,
 unsigned int    num_readbases,
 FILE * fp);

extern void WriteTraceAssemblyTrailer (FILE *fp);


extern void WriteTraceAssemblyFromContig (TContigPtr contig, FILE *fp);

extern void WriteTraceArchiveRead (FILE *fp, TContigReadPtr read);

extern void
WriteFASTAFromAceFile
(TACEFilePtr afp,
 FILE        *fp);

extern void PrintACEFormatErrorXMLStart (char *id, char *has_errors);
extern void PrintACEFormatErrorXMLEnd (void);
extern void PrintACEFormatErrorXML (char *msg, char *id, char *has_errors);

extern int AddReadQualScores (TACEFilePtr afp, FReadLineFunction readfunc, void *userdata, FReadLineFunction fasta_readfunc, void *fasta_userdata);

extern int ReplaceConsensusSequenceFromTraces (TContigPtr contig, char only_ns);
extern void RecalculateConsensusSequences (TACEFilePtr ace_file, char only_ns);

extern void WriteFASTAFromContig (TContigPtr contig, FILE *fp);
extern void WriteContigQualScores (TContigPtr contig, FILE *out);

typedef char (*ProcessContigFunc) (TContigPtr, void *);

extern char
ProcessLargeACEFileForContigFastaAndQualScores
(FReadLineFunction    readfunc,
 void *               userdata,
 char                 make_qual_scores,
 char *               has_errors,
 ProcessContigFunc    process_func,
 void *               process_data);


#ifdef __cplusplus
}
#endif

/*
 * ==========================================================================
 *
 * $Log: aceread.h,v $
 * Revision 1.12  2008/12/02 18:58:24  bollin
 * Added argument to WriteTraceAssemblyHeader for assembly type.
 *
 * Revision 1.11  2008/12/02 18:41:39  bollin
 * Checking in unfinished work on creating pairwise denseg alignment for consensus-read comparison.  Unfinished.
 *
 * Revision 1.10  2008/11/26 18:30:02  bollin
 * Changes to make aceread_tst more efficient when handling large ACE files,
 * added TSA field tags for assembly and taxid.
 *
 * Revision 1.9  2008/11/19 15:21:48  bollin
 * Changes for handling large files.
 *
 * Revision 1.8  2008/11/14 20:16:12  bollin
 * Allow correction of just Ns in consensus sequences.
 *
 * Revision 1.7  2008/11/07 18:28:00  bollin
 * Added functions for reading read FASTA files and quality scores for the read
 * sequences.
 * Also added functions for recalculating the consensus sequence and consensus
 * sequence quality scores based on the reads and read quality scores.
 *
 * Revision 1.6  2008/08/13 15:35:30  bollin
 * Added wrapping header for XML errors during ACE read.
 *
 * Revision 1.5  2008/08/13 14:37:23  bollin
 * Changed error messages to use XML format, removed some unused functions.
 *
 * Revision 1.4  2008/08/13 12:30:01  bollin
 * Changes to allow use of srr numbers in XML and suppress lookups.  Also fixes segfault.
 *
 * Revision 1.3  2008/07/22 19:40:25  kans
 * use brackets instead of quotes on includes, put void in parentheses for no argument prototypes
 *
 * Revision 1.2  2008/07/22 18:45:09  bollin
 * Added function declarations
 *
 * Revision 1.1  2008/07/22 18:10:33  bollin
 * New files for parsing ACE format files.
 *
 *
 * ==========================================================================
 */

#endif /* UTIL_CREADERS___ACEREAD__H */
