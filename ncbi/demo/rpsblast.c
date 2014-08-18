static char const rcsid[] = "$Id: rpsblast.c,v 6.93 2008/07/23 14:06:57 madden Exp $";

/* $Id: rpsblast.c,v 6.93 2008/07/23 14:06:57 madden Exp $
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
* File Name:  $RCSfile: rpsblast.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 12/14/1999
*
* $Revision: 6.93 $
*
* File Description:
*         Main file for RPS BLAST program
*
* $Log: rpsblast.c,v $
* Revision 6.93  2008/07/23 14:06:57  madden
* Fix ASN.1 output (JIRA SB-89)
*
* Revision 6.92  2007/08/21 20:07:01  kans
* include gencode_singleton.h, cast first argument to BlastFormattingInfoNew to fix CodeWarrior complaint
*
* Revision 6.91  2007/03/20 14:56:58  camacho
* Call GeneticCodeSingletonInit/GeneticCodeSingletonFini
*
* Revision 6.90  2007/03/15 14:24:35  coulouri
* added call to FreeSeqLocSetComponents to free query sequences referenced by query_slp list
*
* Revision 6.89  2007/03/12 16:14:23  madden
* Pass a NULL Blast_PsiCheckpointLoc* to Blast_DatabaseSearch [from Mike Gertz].
*
* Revision 6.88  2007/03/05 14:54:39  camacho
* - Call Blast_FindRepeatFilterSeqLoc with a NULL pointer for a PSI-BLAST
*   checkpoint file.
*
* Revision 6.87  2007/02/08 17:07:22  papadopo
* change signature of FillInitialWordOptions; ungapped extensions are always turned on now by default
*
* Revision 6.86  2007/02/01 16:42:44  papadopo
* remove legacy blast engine code
*
* Revision 6.85  2006/08/29 18:12:09  papadopo
* reduce input sequence batch size
*
* Revision 6.84  2006/07/28 21:09:18  papadopo
* allow database length to override the real value when using the rewritten blast engine
*
* Revision 6.83  2006/05/18 16:29:13  papadopo
* do not set search space field directly
*
* Revision 6.82  2006/04/26 12:47:48  madden
* Use SBlastMessage in place of Blast_Message
*
* Revision 6.81  2006/04/21 14:34:50  madden
* BLAST_GetQuerySeqLoc prototype change
*
* Revision 6.80  2006/04/20 15:32:37  papadopo
* if query IDs are actually used, verify that there are no duplicate IDs
*
* Revision 6.79  2006/01/23 16:44:06  papadopo
* change signature of FillHitSavingOptions
*
* Revision 6.78  2006/01/10 20:44:10  madden
* Use SBlastSeqalignArray
*
* Revision 6.77  2006/01/06 16:36:14  papadopo
* 1. By default, log messages go to stderr and not a logfile
* 2. Do not continue execution after a fatal error is encountered
*
* Revision 6.76  2005/12/22 14:22:19  papadopo
* change signature of BLAST_FillLookupTableOptions
*
* Revision 6.75  2005/08/29 14:45:34  camacho
* From Ilya Dondoshansky:
* Retrieve mask_at_hash option from the SBlastOptions structure instead of
* passing as argument in search API calls
*
* Revision 6.74  2005/08/09 21:36:54  dondosha
* Added extra argument in calls to BLAST_GetQuerySeqLoc
*
* Revision 6.73  2005/06/09 17:32:29  dondosha
* Added options validation at the end of s_FillOptions
*
* Revision 6.72  2005/06/08 20:24:48  dondosha
* Previous change wrong - rolled back
*
* Revision 6.70  2005/06/02 20:43:38  dondosha
* 1. Added loop over concatenated query sets;
* 2. Removed unused variables and other clean-up;
* 3. Use BlastFormattingInfo structure for formatting.
*
* Revision 6.69  2005/05/20 18:57:51  camacho
* Update to use new signature to BLAST_FillLookupTableOptions
*
* Revision 6.68  2005/05/02 17:00:28  coulouri
* change default to new engine
*
* Revision 6.67  2005/04/27 14:55:09  papadopo
* change signature of BlastFillHitSavingOptions
*
* Revision 6.66  2005/02/03 18:02:07  dondosha
* Pass summary returns to BLAST_FormatResults, needed for XML output
*
* Revision 6.65  2005/02/02 19:01:36  dondosha
* Use new high level API for performing the search
*
* Revision 6.64  2005/01/10 18:36:47  dondosha
* BlastMaskLocDNAToProtein and BlastMaskLocProteinToDNA moved to core with changed signatures
*
* Revision 6.63  2005/01/10 13:48:20  madden
* Change to BLAST_FillInitialWordOptions prototype
*
* Revision 6.62  2004/12/22 15:38:34  kans
* removed include of blast_tback.h
*
* Revision 6.61  2004/12/21 18:03:04  dondosha
* BLAST_SearchEngine renamed to Blast_RunFullSearch
*
* Revision 6.60  2004/12/07 16:52:57  papadopo
* fix check for range of nucleotide query
*
* Revision 6.59  2004/12/02 20:46:38  dondosha
* Corrected off by 1 error in handling of query location argument
*
* Revision 6.58  2004/12/01 21:15:41  papadopo
* update handling of query range restriction when using the new engine
*
* Revision 6.57  2004/11/23 20:43:06  papadopo
* 1. Add support for the new engine (controlled by -V flag)
* 2. Use enum values instead of hardcoded offsets into myargs[]
* 3. Remove the option for ungapped RPS blast
*
* Revision 6.56  2004/06/30 12:33:30  madden
* Add include for blfmtutl.h
*
* Revision 6.55  2004/04/30 15:33:05  dondosha
* Added argument in call to BXMLPrintOutput
*
* Revision 6.54  2004/02/09 23:14:44  camacho
* Fix to prevent printing multiple Seq-annots in ASN.1 output
*
* Revision 6.53  2004/02/06 05:32:36  camacho
* Properly close ASN.1 stream
*
* Revision 6.52  2003/05/30 17:31:10  coulouri
* add rcsid
*
* Revision 6.51  2003/05/13 16:02:42  coulouri
* make ErrPostEx(SEV_FATAL, ...) exit with nonzero status
*
* Revision 6.50  2003/03/20 14:47:16  madden
* StringSave on asn1_mode
*
* Revision 6.49  2003/03/20 13:44:24  madden
* Fix -m 10/11 output to make them SeqAnnots
*
* Revision 6.48  2002/12/31 22:47:16  boemker
* Added support for printing output as ASN (text, with -m 10, or binary, with
* -m 11).
*
* Revision 6.47  2002/11/29 20:13:13  camacho
* Fix incorrect parameter to FastaToSeqEntryForDb
*
* Revision 6.46  2002/10/17 20:36:00  camacho
* Disallow -L option for tblastn
*
* Revision 6.45  2002/10/13 22:43:51  camacho
* Minor correction
*
* Revision 6.44  2002/10/10 14:49:43  camacho
*
* 1. Removed irrelevant options: -E, -G, -S, -H
* 2. Added -L option to provide range restriction on query sequence
*
* Revision 6.43  2002/09/18 20:34:30  camacho
* Restored -P option
*
* Revision 6.42  2002/08/20 15:17:42  camacho
* Fixed small memory leak
*
* Revision 6.41  2002/08/09 19:41:25  camacho
* 1) Added blast version number to command-line options
* 2) Added explanations for some default parameters
*
* Revision 6.40  2002/06/19 22:50:17  dondosha
* Added all queries information for tabular output with multiple queries
*
* Revision 6.39  2002/04/29 19:55:25  madden
* Use ARG_FLOAT for db length
*
* Revision 6.38  2001/08/28 17:45:01  madden
* Add -m 9 as tabular output with comments
*
* Revision 6.36  2001/06/27 16:20:00  dondosha
* Enabled tabular output for RPS Blast
*
* Revision 6.35  2001/04/13 14:19:08  madden
* Do not print verison banner for XML
*
* Revision 6.34  2001/03/26 14:28:53  madden
* Fix XML problems
*
* Revision 6.33  2000/11/07 18:32:38  shavirin
* Added program version to the output.
*
* Revision 6.32  2000/11/01 20:03:15  shavirin
* Removed not-used option -f for threshold.
*
* Revision 6.31  2000/10/27 19:14:41  madden
* Change description of -b option
*
* Revision 6.30  2000/10/23 19:58:21  dondosha
* Open and close AsnIo outside of call(s) to BXMLPrintOutput
*
* Revision 6.29  2000/10/18 20:55:16  shavirin
* Removed unused command-line parameters.
*
* Revision 6.28  2000/10/02 16:40:53  shavirin
* Fixed setting of TXALIGN_BLASTX_SPECIAL; option.
*
* Revision 6.27  2000/09/29 19:04:57  shavirin
* Fixed warnings and minor errors detected on Windows NT.
*
* Revision 6.26  2000/09/28 18:51:16  shavirin
* Adopted to new parameter BioseqPtr in print results callback.
*
* Revision 6.25  2000/09/27 19:09:04  shavirin
* Significantly redesigned external interface to RPS Blast.
*
* Revision 6.24  2000/08/29 16:54:55  shavirin
* Added option (m = 7) to print XML output.
*
* Revision 6.23  2000/08/21 21:19:16  shavirin
* Removed absolete variable MaxThreadsSem and related code.
*
* Revision 6.22  2000/08/21 19:24:28  madden
* Fix problem writing ASN.1 when multi-threading
*
* Revision 6.21  2000/08/18 19:38:31  madden
* Set believe_query and html flags in AcknowledgeBlastQuery
*
* Revision 6.20  2000/08/04 16:36:05  madden
* Concatenate rather than overwrite SeqAnnot
*
* Revision 6.19  2000/06/28 17:12:29  shavirin
* Fixed problem with -U T option: NULL-ed slp between different sequences.
*
* Revision 6.18  2000/06/27 15:25:19  madden
* Changed master-slave to query-anchored
*
* Revision 6.17  2000/06/20 15:49:35  shavirin
* Added BLAST database title to SeqAnnot output.
*
* Revision 6.16  2000/05/02 18:01:32  shavirin
* Adjusted to changes in RPSInti() function.
*
* Revision 6.15  2000/04/13 18:50:55  shavirin
* Fixed serious memory leaks.
*
* Revision 6.14  2000/04/12 14:15:41  shavirin
* Added back ObjMgrFreeCache together wirg seqmgr changes, those removed
* deadlock.
*
* Revision 6.13  2000/03/28 20:33:44  shavirin
* Changed logic of processing MT - multiple FASTA files.
*
* Revision 6.12  2000/03/10 20:00:15  shavirin
* Added multi-thread support for multi-FASTA files.X
*
* Revision 6.11  2000/03/02 21:06:09  shavirin
* Added -U option, that allows to consider low characters in FASTA files
* as filtered regions (for blastn, blastp and tblastn).
*
* Revision 6.10  2000/02/23 21:03:36  shavirin
* Fixed -z and -Y options in rpsblast.
*
* Revision 6.9  2000/02/17 21:28:55  shavirin
* Added option is_rps_blast = TRUE.
*
* Revision 6.8  2000/02/15 16:14:04  shavirin
* Minor changes.
*
* Revision 6.7  2000/02/11 22:05:01  shavirin
* Oprion do_sum_stats set to FALSE.
*
* Revision 6.6  2000/02/11 20:51:01  shavirin
* Added possibility to search PSSM database against DNA sequences.
*
* Revision 6.5  2000/02/08 17:39:08  shavirin
* Empty log message.
*
* Revision 6.4  2000/02/01 17:22:48  shavirin
* Updated function RPSViewSeqAlign().
*
* Revision 6.3  2000/01/07 22:34:05  shavirin
* Added printing of SeqAlignment if necessary.
*
* Revision 6.2  1999/12/30 18:37:53  shavirin
* Added NCBI header and Log information.
*
*
* ==========================================================================
*/

#include <readdb.h>
#include <algo/blast/api/blast_input.h>
#include <algo/blast/api/blast_format.h>
#include <algo/blast/api/blast_api.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/blast_message_api.h>
#include <algo/blast/core/gencode_singleton.h>

#ifdef PURIFY
#include "/am/purew/solaris2/new/../purify/purify-4.5-solaris2/purify.h"
#endif

enum {
    OPT_QUERY_FILE = 0,
    OPT_DB,
    OPT_PROT_QUERY,
    OPT_EVALUE,
    OPT_FORMAT,
    OPT_OUTPUT_FILE,
    OPT_XDROP_UNGAPPED,
    OPT_UNGAPPED_HITS,
    OPT_FILTER,
    OPT_XDROP_GAPPED,
    OPT_GAP_TRIGGER,
    OPT_NUMTHREADS,
    OPT_SHOW_GI,
    OPT_BELIEVE_QUERY,
    OPT_XDROP_GAPPED_FINAL,
    OPT_ASNOUT,
    OPT_NUM_DESC,
    OPT_NUM_RESULTS,
    OPT_DBLENGTH,
    OPT_SEARCHSP,
    OPT_HTML,
    OPT_LOGFILE,
    OPT_LCASE,
    OPT_RANGE,
    NUM_ARGS    /* must come last */
};
    
static Args myargs[] = {
    /* OPT_QUERY_FILE */
    {"Input query sequence (this parameter must be set)",
     "stdin", NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    /* OPT_DB */
    {"RPS BLAST Database",
     NULL, NULL,NULL,FALSE,'d',ARG_FILE_IN, 0.0,0,NULL},
    /* OPT_PROT_QUERY */
    {"Query sequence is protein ",
     "T", NULL,NULL,TRUE, 'p', ARG_BOOLEAN, 0.0,0,NULL},
    /* OPT_EVALUE */
    { "Expectation value (E)",
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
    /* OPT_FORMAT*/
    { "alignment view options:\n"
"0 = pairwise,\n"
"1 = query-anchored showing identities,\n"
"2 = query-anchored no identities,\n"
"3 = flat query-anchored, show identities,\n"
"4 = flat query-anchored, no identities,\n"
"5 = query-anchored no identities and blunt ends,\n"
"6 = flat query-anchored, no identities and blunt ends,\n"
"7 = XML Blast output,\n8 = tabular output, \n"
"9 = tabular output with comments",
      "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
    /* OPT_OUTPUT_FILE */
    { "Output File for Alignment",
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
    /* OPT_XDROP_UNGAPPED */
    { "Dropoff (X) for blast extensions in bits (default if zero)",
      "7.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},
    /* OPT_UNGAPPED_HITS */
    { "0 for multiple hit, 1 for single hit",
       "0",  NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},
    /* OPT_FILTER */
    { "Filter query sequence with SEG",
      "F", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
    /* OPT_XDROP_GAPPED */
    { "X dropoff value for gapped alignment (in bits)",
      "15", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
    /* OPT_GAP_TRIGGER */
    { "Number of bits to trigger gapping",
      "22.0", NULL, NULL, FALSE, 'N', ARG_FLOAT, 0.0, 0, NULL},
    /* OPT_NUMTHREADS */
    { "Number of processors to use",
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
    /* OPT_SHOW_GI */
    { "Show GI's in deflines",
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
    /* OPT_BELIEVE_QUERY */
    { "Believe the query defline",
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
    /* OPT_XDROP_GAPPED_FINAL */
    { "X dropoff value for final gapped alignment (in bits)",
      "25", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
    /* OPT_ASNOUT */
    { "SeqAlign file ('Believe the query defline' must be TRUE)",
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
    /* OPT_NUM_DESC */
    { "Number of database sequences to show one-line descriptions for (V)",
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
    /* OPT_NUM_RESULTS */
    { "Number of database sequence to show alignments for (B)",
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
    /* OPT_DBLENGTH */
    { "Effective length of the database (use zero for the real size)",
      "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},
    /* OPT_SEARCHSP */
    { "Effective length of the search space (use zero for the real size)",
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
    /* OPT_HTML */
    { "Produce HTML output",
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
    /* OPT_LOGIFLE */
    {"Logfile name ",
     "stderr", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    /* OPT_LCASE */
    {"Use lower case filtering of FASTA sequence",
     "F", NULL,NULL,TRUE,'U',ARG_BOOLEAN, 0.0,0,NULL},
    /* OPT_RANGE */
    { "Range restriction on query sequence (Format: start,stop) blastp only\n"
      "      0 in 'start' refers to the beginning of the sequence\n"
      "      0 in 'stop' refers to the end of the sequence",
      "0,0", NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},
};


/* This is very similar to the routine of the same name
   in blast_driver.c; unneeded setup is removed */

static Int2 
s_FillOptions(SBlastOptions* options, Blast_SummaryReturn* sum_returns)
{
   Int2 status;
   const EBlastProgramType kProgram = options->program;
   LookupTableOptions* lookup_options = options->lookup_options;
   QuerySetUpOptions* query_setup_options = options->query_options; 
   BlastInitialWordOptions* word_options = options->word_options;
   BlastExtensionOptions* ext_options = options->ext_options ;
   BlastHitSavingOptions* hit_options = options->hit_options;
   BlastScoringOptions* score_options = options->score_options;
   BlastEffectiveLengthsOptions* eff_len_options = options->eff_len_options;
   BlastDatabaseOptions* db_options = options->db_options;
   Blast_Message *core_msg = NULL;

   BLAST_FillLookupTableOptions(lookup_options, kProgram, 
                                FALSE, /* megablast */
                                0,     /* default threshold */
                                0);    /* default wordsize */

   BLAST_FillQuerySetUpOptions(query_setup_options, kProgram, 
                              myargs[OPT_FILTER].strvalue, 0);

   BLAST_FillInitialWordOptions(word_options, kProgram,
                                0, myargs[OPT_XDROP_UNGAPPED].intvalue);

   /* set one-hit extensions if desired */
   if (myargs[OPT_UNGAPPED_HITS].intvalue == 1)
       word_options->window_size = 0;

   BLAST_FillExtensionOptions(ext_options, kProgram, FALSE, 
                              myargs[OPT_XDROP_GAPPED].intvalue, 
                              myargs[OPT_XDROP_GAPPED_FINAL].intvalue);

   BLAST_FillScoringOptions(score_options, kProgram, FALSE,
                0, 0,  /* nucleotide match/mismatch scores */
                NULL, 0, 0); /* Matrix, gap opening and gap extension penalties
                                will be set from the RPS database, so there is
                                no need to set them here. */

   score_options->gapped_calculation = TRUE;

   BLAST_FillHitSavingOptions(hit_options, 
                      myargs[OPT_EVALUE].floatvalue, 
                      MAX(myargs[OPT_NUM_DESC].intvalue, 
                          myargs[OPT_NUM_RESULTS].intvalue),
                      TRUE,
                      0,        /* turn off culling */
                      0);       /* min diag separation */

   if (myargs[OPT_SEARCHSP].floatvalue != 0 ||
       myargs[OPT_DBLENGTH].floatvalue != 0) {
      Int8 searchsp = (Int8) myargs[OPT_SEARCHSP].floatvalue; 
      Int8 dblength = (Int8) myargs[OPT_DBLENGTH].floatvalue; 
      BLAST_FillEffectiveLengthsOptions(eff_len_options, 0, 
                                        dblength, &searchsp, 1);
   }

   if (db_options && kProgram == eBlastTypeRpsTblastn) {
       Uint1* gc = NULL;
       BLAST_GeneticCodeFind(db_options->genetic_code, &gc);
       GenCodeSingletonAdd(db_options->genetic_code, gc);
       free(gc);
   }

   /* Validate the options. */
   status = BLAST_ValidateOptions(kProgram, ext_options, score_options, 
                                  lookup_options, word_options, hit_options, 
                                  &core_msg);
   sum_returns->error = Blast_MessageToSBlastMessage(core_msg, NULL, NULL, FALSE);

   return status;
}

static Int2
s_ParseIntervalLocationArgument(char* arg, Int4* from_ptr, Int4* to_ptr)
{
   const char* delimiters = " ,;";
   Int4 from = 0, to = 0;

   if (!arg)
      return 0;
   
   from = atoi(StringTokMT(arg, delimiters, &arg));
   to = atoi(arg);
      
   from = MAX(from, 0);
   to = MAX(to, 0);

   *from_ptr = from;
   *to_ptr = to;

   if (from > to)
      return -1;
   else
      return 0;
}

/* Similar to the Main() of blast_driver.c; major differences are
        - only one thread supported
        - no two-sequences capability
        - on-the-fly tabular output not supported
*/
Int2 Main(void)
{
   Int2 status;
   Boolean query_is_na;
   char* blast_program = NULL;
   EBlastProgramType program_number;
   char* dbname = NULL;
   FILE *infp;
   Int4 ctr = 1;
   Int4 query_from = 0, query_to = 0;
   SBlastOptions* options = NULL;
   BlastFormattingInfo* format_info = NULL;
   Blast_SummaryReturn* sum_returns=Blast_SummaryReturnNew();
   Int4 num_queries = 0;
   SeqLoc* lcase_mask = NULL;
   const int kMaxConcatLength = 10000;
   Blast_SummaryReturn* full_sum_returns = NULL;
   Boolean believe_query = (Boolean) myargs[OPT_BELIEVE_QUERY].intvalue;
   Char buf[256] = { '\0' };
   BlastFormattingInfo* asn_format_info = NULL;
   GeneticCodeSingletonInit();

   StringCpy(buf, "rpsblast ");
   StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf)-1);
   if (!GetArgs(buf, NUM_ARGS, myargs))
       return 1;
   
   if (strcmp(myargs[OPT_LOGFILE].strvalue, "stderr") != 0) {
      if ( !ErrSetLog (myargs[OPT_LOGFILE].strvalue) ) { /* Logfile */
         ErrShow();
         return 1;
      }
   }

   /* select protein or nucleotide query, and choose
      the appropriate program type */

   if (myargs[OPT_PROT_QUERY].intvalue == FALSE) {
       query_is_na = TRUE;
       blast_program = "rpstblastn";
   } else {
       query_is_na = FALSE;
       blast_program = "rpsblast";
   }

   status = SBlastOptionsNew(blast_program, &options, sum_returns);

   if (!status)
       s_FillOptions(options, sum_returns);

   if (status) {
       options = SBlastOptionsFree(options);
       if (sum_returns->error) {
           SBlastMessageErrPost(sum_returns->error);
           sum_returns = Blast_SummaryReturnFree(sum_returns);
       }
       return -1;
   }

   SBlastOptionsSetBelieveQuery(options, believe_query);

   /* initialize the database and then the RPS_specific data files */

   dbname = myargs[OPT_DB].strvalue;

   program_number = options->program;

   /* Get the query */

   if ((infp = FileOpen(myargs[OPT_QUERY_FILE].strvalue, "r")) == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "Unable to open input file '%s'\n", 
                myargs[OPT_QUERY_FILE].strvalue);
      return (1);
   }

   if (s_ParseIntervalLocationArgument(myargs[OPT_RANGE].strvalue,
                                        &query_from, &query_to)) {
      ErrPostEx(SEV_FATAL, 1, 0, "Invalid subject sequence location\n");
      return -1;                   
   }

   if ((query_from != 0 || query_to != 0) && 
       (Boolean)myargs[OPT_PROT_QUERY].intvalue == FALSE) {
      ErrPostEx(SEV_FATAL, 1, 0, "No query range allowed for nucleotide query");
      return -1;                   
   }

   BlastFormattingInfoNew((EAlignView) myargs[OPT_FORMAT].intvalue, options, blast_program,
                          dbname, myargs[OPT_OUTPUT_FILE].strvalue,
                          &format_info);

   /* This is not megablast, so pass FALSE for the respective arguments. */
   BlastFormattingInfoSetUpOptions(format_info, myargs[OPT_NUM_DESC].intvalue, 
                                   myargs[OPT_NUM_RESULTS].intvalue,
                                   (Boolean) myargs[OPT_HTML].intvalue,
                                   FALSE, (Boolean) myargs[OPT_SHOW_GI].intvalue,
                                   believe_query);

   BLAST_PrintOutputHeader(format_info);
   if (myargs[OPT_ASNOUT].strvalue) {
               /* This just prints out the ASN.1 to a secondary file. */
               BlastFormattingInfoNew(eAlignViewAsnText, options,
                            blast_program, dbname,
                            myargs[OPT_ASNOUT].strvalue, &asn_format_info);

               BlastFormattingInfoSetUpOptions(asn_format_info,
                            myargs[OPT_NUM_DESC].intvalue,
                            myargs[OPT_NUM_DESC].intvalue,
                            FALSE,
                            FALSE,
                            FALSE,
                            TRUE);
   }

   /* Loop over sets of queries. */
   while (1) {
       SBlastSeqalignArray* seqalign_arr=NULL;
       SeqLoc* query_slp = NULL;
       SeqLoc* filter_loc = NULL;
       Int4 letters_read;

       if ((Boolean)myargs[OPT_LCASE].intvalue) {
           letters_read = 
               BLAST_GetQuerySeqLoc(infp, query_is_na, 0, kMaxConcatLength, 
                                    query_from, query_to, &lcase_mask,
                                    &query_slp, &ctr, &num_queries, 
                                    believe_query, 0);
       } else {
           letters_read = 
               BLAST_GetQuerySeqLoc(infp, query_is_na, 0, kMaxConcatLength,
                                    query_from, query_to, NULL, 
                                    &query_slp, &ctr, &num_queries, 
                                    believe_query, 0);
       }
       
       if (letters_read <= 0)
           break;
       
       if (believe_query && BlastSeqlocsHaveDuplicateIDs(query_slp)) {
          ErrPostEx(SEV_FATAL, 1, 0, 
                  "Duplicate IDs detected; please ensure that "
                  "all query sequence identifiers are unique");
       }

       /* Call database search function. Pass NULL for tabular formatting 
        * structure pointer, because on-the-fly tabular formatting is not
        * allowed for RPS BLAST.  Pass NULL for a PSIBLAST checkpoint file
        * as well, since this parameter is not relevant to RPS BLAST.
        */
       status =
           Blast_DatabaseSearch(query_slp, (Blast_PsiCheckpointLoc *) NULL,
                                dbname, lcase_mask, options,
                                (BlastTabularFormatData*) NULL,
                                &seqalign_arr, &filter_loc,
                                sum_returns);

       /* Free the lower case mask in SeqLoc form. */
       lcase_mask = Blast_ValNodeMaskListFree(lcase_mask);
       
       /* If masking was done for lookup table only, free the masking locations,
          because they will not be used for formatting. */
       if (SBlastOptionsGetMaskAtHash(options))
           filter_loc = Blast_ValNodeMaskListFree(filter_loc);
       
       /* Post warning or error messages, no matter what the search status 
          was. */
       SBlastMessageErrPost(sum_returns->error);

       if (status != 0) {
           ErrPostEx(SEV_FATAL, 1, 0, "BLAST search failed");
           return status;
       }
       
       /* do initial cleanup */
       
       /* format results */
       
       if (myargs[OPT_ASNOUT].strvalue) {
                   status =
                       BLAST_FormatResults(seqalign_arr, num_queries, query_slp,
                                   NULL, asn_format_info, sum_returns);
       }
       
       status = 
           BLAST_FormatResults(seqalign_arr, num_queries, query_slp, filter_loc,
                               format_info, sum_returns);
       
       /* finish cleanup */
       filter_loc = Blast_ValNodeMaskListFree(filter_loc);
       FreeSeqLocSetComponents(query_slp);
       query_slp = SeqLocSetFree(query_slp);
       seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
       
       /* Update the cumulative summary returns structure and clean the returns
          substructures for the current search iteration. */
       Blast_SummaryReturnUpdate(sum_returns, &full_sum_returns);
       Blast_SummaryReturnClean(sum_returns);
   }

   if (infp)
      FileClose(infp);
    if (asn_format_info)
      asn_format_info = BlastFormattingInfoFree(asn_format_info);
   
   /* Print the footer with summary information. */
   Blast_PrintOutputFooter(format_info, full_sum_returns);

   sum_returns = Blast_SummaryReturnFree(sum_returns);
   full_sum_returns = Blast_SummaryReturnFree(full_sum_returns);

   format_info = BlastFormattingInfoFree(format_info);
   options = SBlastOptionsFree(options);
   GeneticCodeSingletonFini();
   
   return status;
}
