#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = "$Id: blast_input.c,v 1.32 2008/06/09 17:28:41 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */
/* ===========================================================================
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
* Author: Ilya Dondoshansky
*
*/

/** @file blast_input.c
 * Reading FASTA sequences for BLAST
 */

#include <objloc.h>
#include <tofasta.h>
#include <algo/blast/api/blast_input.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_posit.h>
#include <objscoremat.h>
#include <asn.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** Maximal number of queries allowed in one SeqLoc chain (== 1/2 INT2_MAX) */
#define MAX_NUM_QUERIES 16383 
/** Maximal total length of queries in a single SeqLoc chain. */
#define MAX_TOTAL_LENGTH 2000000

Int4
BLAST_GetQuerySeqLoc(FILE *infp, Boolean query_is_na, Uint1 strand, 
                     Int4 max_total_length, Int4 start, Int4 end, 
                     SeqLoc** lcase_mask, SeqLocPtr* query_slp, 
                     Int4Ptr ctr, Int4* num_queries, Boolean believe_query,
                     Int4 genetic_code)
{
   Int4 total_length=0; /* total number of letters read this call, also final 
                           return value. */
   SeqEntryPtr sep;
   SeqLocPtr mask_slp, last_slp;
   char prefix[4];     /* for FastaToSeqEntryForDb */
   Int4 query_index = 0;  /* number of sequences read. */
   ValNodePtr vnp=NULL; /* used to keep lower-case masking SeqLoc's */
   Int2 query_count;    /* Query count on this call for FastaToSeqEntryForDb */

   if (!query_slp)
   {
      ErrPostEx(SEV_FATAL, 0, 0, "NULL query_slp obtained in BLAST_GetQuerySeqLoc");
      return -1;
   }

   if (!ctr)  /* Not providing this can cause problems if multiple calls to this
                 function are made. */
   {
      ErrPostEx(SEV_FATAL, 0, 0, "ctr must be non-NULL in BLAST_GetQuerySeqLoc");
      return -1;
   }

   if (max_total_length <= 0)
     max_total_length = MAX_TOTAL_LENGTH;

   *query_slp = NULL;
   last_slp = NULL;

   if (query_is_na && strand == Seq_strand_unknown)
      strand = Seq_strand_both;

   SeqMgrHoldIndexing(TRUE);
   mask_slp = NULL;
   if (lcase_mask) /* Make sure we don't get old (possibly freed) locations. */
     *lcase_mask = NULL;
   
   /* This is a workaround on the Int2 ctr input for FastaToSeqEntryForDb */
   sprintf(prefix, "%ld", (long) (*ctr/MAX_NUM_QUERIES));
   query_count = *ctr%MAX_NUM_QUERIES;
   
   while ((sep=FastaToSeqEntryForDb(infp, query_is_na, NULL, believe_query, prefix, 
                                    &query_count, (lcase_mask ? &mask_slp : NULL))) != NULL)
   {
      BioseqPtr query_bsp;
      Int4 from, to;

      if (lcase_mask)  /* Only keep if lcase masking is being read in. */
         ValNodeAddPointer(&vnp, 0, mask_slp);
      ++query_index;
      mask_slp = NULL;
      
      query_bsp = NULL;
      if (query_is_na) {
         SeqEntryExplore(sep, &query_bsp, FindNuc);
      } else {
         SeqEntryExplore(sep, &query_bsp, FindProt);
      }

      if (query_bsp == NULL) {
         ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq for sequence number %ld", (long) query_count-1);
         *ctr += query_index;
         return -1;
      }

      if (query_bsp->length <= 0) {
         ErrPostEx(SEV_WARNING, 0, 0, "Sequence number %ld had length %ld", (long) query_count-1, (long) query_bsp->length);
         *ctr += query_index;
      }


      
      /* Original from and to are 1-offsets, except when they are 0's,
         in which case they are start and end of sequence respectively */
      from = ((start > 0) ? start - 1 : 0);
      to = ((end > 0) ? end - 1 : query_bsp->length - 1);

      to = MIN(to, query_bsp->length - 1);

      /* If location starting offset is larger than sequence length, skip this
         sequence. */
      if (from > to) 
         continue;

      /* Fill the query genetic code option. */
      if (query_is_na && genetic_code > 0) {
          BioSourcePtr source;
          source = BioSourceNew();
          source->org = OrgRefNew();
          source->org->orgname = OrgNameNew();
          source->org->orgname->gcode = genetic_code;
          ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
      }

      if ((strand == Seq_strand_plus) || (strand == Seq_strand_minus) ||
          (from > 0) || (to < query_bsp->length - 1))
      {
         SeqLocPtr new_slp = SeqLocIntNew(from, to, strand, 
                                SeqIdFindBest(query_bsp->id, SEQID_GI)); 
         if (last_slp) {
            last_slp->next = new_slp;
            last_slp = last_slp->next;
         } else {
            *query_slp = last_slp = new_slp;
         }
      } else {
         last_slp = ValNodeAddPointer(&last_slp, SEQLOC_WHOLE, 
                       SeqIdDup(SeqIdFindBest(query_bsp->id, SEQID_GI)));
         if (*query_slp == NULL)
            *query_slp = last_slp;
      }

      total_length += query_bsp->length;
      if (total_length > max_total_length || query_index >= MAX_NUM_QUERIES) {
         break;  /* Read maximum allowed amount of data. */
      }
   }

   if (lcase_mask)
       *lcase_mask = vnp;

   SeqMgrHoldIndexing(FALSE);

   if (num_queries)
      *num_queries = query_index;

   *ctr += query_index;

   return total_length;
}


/**
 * Read the query data stored in a PSI-BLAST checkpoint file.
 *
 * @param *pold_query_length   length of the old query, as stored in the
 *                             checkpoint file
 * @param old_query            a buffer of length at least query_length
 *                             to hold the query data from the checkpoint
 *                             file
 * @param query_length         the expected query length, and the length
 *                             of the old_query buffer
 * @param file                 the checkpoint file
 * @return                     the number of query characters read,
 *                             <= the value of the query_length parameter.
 */
static int
s_GetOldQueryFromCheckpoint(Int4 * pold_query_length, Uint1 * old_query,
                            int query_length, FILE * file)
{
    int i;
    int num_read;    /* number of query characters read from disk */

    ASSERT(query_length > 0);
    if (1 != fread(pold_query_length, sizeof(Int4), 1, file)) {
        *pold_query_length = 0;
        return 0;
    }
    /* Read at most query_length bytes, no matter what the value of
     * *pold_query_length is */
    num_read = (int)fread(old_query, sizeof(Uint1), query_length, file);
    for (i = 0;  i < num_read;  i++) {
        old_query[i] = AMINOACID_TO_NCBISTDAA[old_query[i]];
    }
    return num_read;
}


/**
 * Compare the query to another sequence, typically a query stored in
 * a PSI-BLAST checkpoint file.  Warn if they don't match.
 * @param query          one query to be compared
 * @param query_length   the length of query and old_query
 * @param old_query      the other query
 * @return 0 if the sequences match, -1 otherwise */
static Int4
s_ValidateOldQuery(const Uint1 query[], int query_length,
                   const Uint1 old_query[], int old_query_length)
{
    int query_index;
    /* Value for the X ambiguity character */
    enum { eXchar = 21 };

    if (query_length != old_query_length) {
        ErrPostEx(SEV_WARNING, 0, 0, "Invalid usage of checkpoint recovery; "
                  "old query has length %ld, new query has length %ld",
                  (long) old_query_length,  (long) query_length);
        return -1;
    }
    for (query_index = 0;  query_index < query_length; query_index++) {
        if (old_query[query_index] != query[query_index]) {
            char old_char = NCBISTDAA_TO_AMINOACID[old_query[query_index]];
            char new_char = NCBISTDAA_TO_AMINOACID[query[query_index]];

            if (old_query[query_index] != eXchar) {
                if (query[query_index] == eXchar) {
                    ErrPostEx(SEV_WARNING, 0, 0,
                              "\nStored query has a %c at position "
                              "%d, while new query has a %c there.\n%c "
                              "appears in query sequence: The query "
                              "could be filtered. Run with \"-F F\" "
                              "option to turn the filter off.",
                              old_char, query_index, new_char, new_char);
                } else {
                    ErrPostEx(SEV_WARNING, 0, 0,
                              "Stored query has a %c at position %d, "
                              "while new query has a %c there.",
                              old_char, query_index, new_char);
                }
                return -1;
            } else { /* old_query[c] == eXchar */
                ErrPostEx(SEV_WARNING, 0, 0,
                          "Stored query has a %c at position %d, "
                          "while new query has a %c there\n%c appears "
                          "in the stored query: The stored query may be "
                          "filtered.  Run blastpgp with \"-F F\" option "
                          "to turn the filter off.",
                          old_char, query_index, new_char, old_char);
                /* This is only a warning */
            }
        }
    }
    return 0;
}


/**
 * Read frequency ratios from a PSI-BLAST checkpoint file; read the
 * query data for the file before calling this routine.
 *
 * @param freq_ratios     the frequency ratios [out]
 * @param query_length    the length of the query
 * @param file            file to read
 *
 * @return 0 on sucess, nonzero on error */
static int
s_PosReadFreqRatios(double ** freq_ratios,
                    int qlength,
                    FILE * file)
{
    int query_index;  /* Query position (column) in the PSSM */
    int stdaa_index;  /* Index in the ncbi_stdaa alphabet */
    int trueaa_index; /* Index in the ARND... alphabet of true
                         amino acids */
    /* conversion from 28 letter NCBIstdaa alphabet to 20 letter order
     * for true amino acids: ARNDCQEGHILKMFPSTWYV. */
    static int alphaConvert[BLASTAA_SIZE] =
        {(-1), 0, (-1),  4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
         16, 19, 17, (-1), 18, (-1), (-1), (-1), (-1), (-1)};
    /* Buffer to hold frequency data in ARND... order */
    double trueaa_buffer[PRO_TRUE_ALPHABET_SIZE];
    int num_read;   /* Number of items retrieved by a call to read */

    for (query_index = 0;  query_index < qlength;  query_index++) {
        num_read = (int)fread(trueaa_buffer, sizeof(double),
                              PRO_TRUE_ALPHABET_SIZE, file);
        if (num_read < PRO_TRUE_ALPHABET_SIZE)
            return -1;
        for (stdaa_index = 0;  stdaa_index < BLASTAA_SIZE;  stdaa_index++) {
            trueaa_index = alphaConvert[stdaa_index];
            if (trueaa_index < 0) {
                freq_ratios[query_index][stdaa_index] = 0.0;
            } else {
                freq_ratios[query_index][stdaa_index] =
                    trueaa_buffer[trueaa_index];
            }
        }
    }
    return 0;
}


/**
 * Read frequency ratios from a standard format PSI-BLAST checkpoint file.
 *
 * @param freq_ratios   the frequency ratios
 * @param query_length  the length of the query, and second dimension of
 *                      freq_ratios
 * @param query         query sequence data
 * @param file          an open file to be read
 * @param blast_msg     a pointer to hold BLAST warnings.
 *
 * @return 0 on success, nonzero otherwise
 */
static int
s_PosReadStdCheckpointFile(double ** freq_ratios,
                           int query_length,
                           const Uint1 * query,
                           FILE * file,
                           Blast_Message* *blast_msg)
{
    int chkpt_query_length = 0;  /* Length of the query saved in the
                                    checkpoint file */
    int num_read = 0;  /* number if query characters actually read from
                          the checkpoint file */
    int status = 0;    /* error status */
    /* Buffer to hold the query from the checkpoint file */
    Uint1 * chkpt_query = NULL;

    chkpt_query = calloc(query_length, sizeof(Uint1));
    if (NULL != chkpt_query) {
        num_read = s_GetOldQueryFromCheckpoint(&chkpt_query_length,
                                               chkpt_query, query_length,
                                               file);
    }
    if (NULL == chkpt_query || num_read < chkpt_query_length) {
        Blast_MessageWrite(blast_msg, eBlastSevFatal,
                           kBlastMessageNoContext,
                           "Blast_PosReadCheckpoint: "
                           "Failed to reconstruct previous query\n");
        goto error_return;
    }
    status = s_ValidateOldQuery(query, query_length,
                                chkpt_query, chkpt_query_length);
    free(chkpt_query);

    if (0 != status)
        goto error_return;

    status = s_PosReadFreqRatios(freq_ratios, query_length, file);
    if (0 != status)
        goto error_return;

    return 0;
error_return:
    Blast_MessageWrite(blast_msg, eBlastSevFatal,
                       kBlastMessageNoContext,
                       "Blast_PosReadCheckpoint: "
                       "Failed to recover data\n");
    return -1;
}


/**
 * Read frequency ratios from a standard format PSI-BLAST checkpoint file
 * of the given name.
 *
 * @param freq_ratios   the frequency ratios
 * @param query_length  the length of the query, and second dimension of
 *                      freq_ratios
 * @param query         query sequence data
 * @param file          the name of the file to be read
 * @param blast_msg     a pointer to hold BLAST warnings.
 *
 * @return 0 on success, nonzero otherwise
 */
static int
s_PosReadStdCheckpoint(Nlm_FloatHi ** freq_ratios,
                       int qlength,
                       const Uint1 * query,
                       const char * filename,
                       Blast_Message* *blast_msg)
{
    FILE * file = fopen(filename, "rb");  
    int status, close_status;

    if (!file) {
     Blast_MessageWrite(blast_msg, eBlastSevFatal,
                        kBlastMessageNoContext,
                        "Blast_PosReadCheckpointFile: "
                        "Could not open checkpoint file\n");
        return -1;
    }
    status = s_PosReadStdCheckpointFile(freq_ratios, qlength, query,
                                        file, blast_msg);
    close_status = fclose(file);

    return (0 == status) ? close_status : status;
}


static int
s_PosReadAsnCheckpoint(double ** freq_ratios,
                       int query_length,
                       const Uint1 query[],
                       char fileName[],
                       int is_ascii_scoremat)
{
    AsnIoPtr infile = NULL;
    PssmWithParametersPtr scoremat = NULL;
    PssmPtr pssm = NULL;
    PssmIntermediateDataPtr freqs = NULL;
    ValNodePtr freq_list;
    Bioseq *bsp;
    int i, j, c;
    enum { eXchar = 21 };

    if (is_ascii_scoremat) {
        infile = AsnIoOpen(fileName, "r");
    } else {
        infile = AsnIoOpen(fileName, "rb");
    }
    if (infile == NULL) {
        ErrPostEx(SEV_WARNING, 0, 0,"Could not open scoremat file\n");
        return -1;
    }

    scoremat = PssmWithParametersAsnRead(infile, NULL);
    AsnIoClose(infile);
    if (scoremat == NULL) {
        ErrPostEx(SEV_WARNING, 0, 0,
                  "Could not read scoremat from input file\n");
        return 1;
    }
    pssm = scoremat->pssm;
    if (pssm == NULL) {
        ErrPostEx(SEV_WARNING, 0, 0,"Scoremat is empty\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    freqs = pssm->intermediateData;
    if (freqs == NULL) {
        ErrPostEx(SEV_WARNING, 0, 0,
                  "Scoremat doesn't contain intermediate data\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (freqs->freqRatios == NULL) {
        ErrPostEx(SEV_WARNING, 0, 0,
                  "Scoremat does not contain frequency ratios\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (pssm->numRows != BLASTAA_SIZE) {
        ErrPostEx(SEV_WARNING, 0, 0, "Wrong alphabet size of %d in "
                  "input scoremat\n", pssm->numRows);
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (!pssm->query || !pssm->query->data.ptrvalue) {
        ErrPostEx(SEV_WARNING, 0, 0, 
                  "Missing sequence data in input scoremat\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    bsp = (Bioseq *)(pssm->query->data.ptrvalue);
    if (pssm->numColumns != bsp->length) {
        ErrPostEx(SEV_WARNING, 0, 0, "Different sequence lengths "
                  "(%d and %d) in input scoremat\n", pssm->numColumns,
                  bsp->length);
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (pssm->numColumns != query_length) {
        ErrPostEx(SEV_WARNING, 0, 0, "Scoremat sequence length "
                  "(%d) does not match query length (%d)\n",
                  pssm->numColumns, query_length);
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (!bsp->seq_data || !ISA_aa(bsp->mol)) {
        ErrPostEx(SEV_WARNING, 0, 0,
                  "Sequence within checkpoint file has no data or is "
                  "not protein\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }

    if (bsp->seq_data_type == Seq_code_gap) {
        ErrPostEx(SEV_WARNING, 0, 0,"Seq_code_gap passed to s_PosReadAsnCheckpoint\n");
        return -1;
    }

    BSSeek((ByteStorePtr) bsp->seq_data, 0, SEEK_SET);

    /* Convert sequence data into Seq_code_ncbistdaa */
    if (bsp->seq_data_type != Seq_code_ncbistdaa) {

        ByteStore* new_byte_store = BSConvertSeq((ByteStorePtr) bsp->seq_data,
                                                 Seq_code_ncbistdaa,
                                                 bsp->seq_data_type,
                                                 bsp->length);
        if ( !new_byte_store ) {
            ErrPostEx(SEV_FATAL, 1, 0, "Failed to convert Bioseq in ASN.1"
                      " PSSM to Seq_code_ncbistdaa");
        }
        bsp->seq_data = (SeqDataPtr) new_byte_store;
        bsp->seq_data_type = Seq_code_ncbistdaa;
        BSSeek((ByteStorePtr) bsp->seq_data, 0, SEEK_SET);
    }

    /* verify the input query is the same as the sequence
       within the checkpoint file */

    for (i = 0; i < query_length; i++) {
        c = BSGetByte((ByteStorePtr) bsp->seq_data);
        if (c == EOF) {
            ErrPostEx(SEV_WARNING, 0, 0, "Premature end of sequence data\n");
            PssmWithParametersFree(scoremat);
            return -1;
        }
        if (c != query[i]) {
           char old_char = NCBISTDAA_TO_AMINOACID[query[i]];
           char new_char = NCBISTDAA_TO_AMINOACID[c];
           if (query[i] == eXchar) {
                ErrPostEx(SEV_WARNING, 0, 0,
                          "Query sequence contains '%c' at position %d; "
                          "if filtering was used, rerun the search with "
                          "filtering turned off ('-F F')\n", old_char, i);
            }
            else {
                ErrPostEx(SEV_WARNING, 0, 0,
                          "Query sequence contains '%c' at position %d, "
                          "while sequence withing checkpoint file contains "
                          "'%c' at this position\n",
                          old_char, i, new_char);
            }
            PssmWithParametersFree(scoremat);
            return -1;
        }
    }

    /* Read in the frequency ratios, verify they fall
       in the correct range, and verify that the linked list
       of residue frequencies is exactly as long as it should be */

    freq_list = freqs->freqRatios;
    /* This is bad */
    if (pssm->byRow == FALSE) {
        j = 0;
        for (i = 0; i < pssm->numColumns; i++) {
            for (j = 0; j < pssm->numRows; j++) {
                if (freq_list == NULL)
                    break;
                freq_ratios[i][j] = freq_list->data.realvalue;

                if (freq_ratios[i][j] < 0.0) {
                    ErrPostEx(SEV_WARNING, 0, 0, "position frequency (%d,%d) "
                              "out of bounds\n", i, j);
                    PssmWithParametersFree(scoremat);
                    return -1;
                }
                freq_list = freq_list->next;
            }
            if (j < pssm->numRows)
                break;
        }
    } else {
        i = 0;
        for (j = 0; j < pssm->numRows; j++) {
            for (i = 0; i < pssm->numColumns; i++) {
                if (freq_list == NULL)
                    break;
                freq_ratios[i][j] = freq_list->data.realvalue;

                if (freq_ratios[i][j] < 0.0) {
                    ErrPostEx(SEV_WARNING, 0, 0, "position frequency (%d,%d) "
                              "out of bounds\n", i, j);
                    PssmWithParametersFree(scoremat);
                    return -1;
                }
                freq_list = freq_list->next;
            }
            if (i < pssm->numColumns)
                break;
        }
    }
    if (i < pssm->numColumns || j < pssm->numRows) {
        ErrPostEx(SEV_WARNING, 0, 0, "Not enough frequency "
                  "ratios in input scoremat\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    if (freq_list != NULL) {
        ErrPostEx(SEV_WARNING, 0, 0, "Too many frequency "
                  "ratios in input scoremat\n");
        PssmWithParametersFree(scoremat);
        return -1;
    }
    PssmWithParametersFree(scoremat);
    return 0;
}


Blast_PsiCheckpointLoc *
Blast_PsiCheckpointLocNew(EPsiCheckpointType checkpoint_type,
                          char * filename)
{
    Blast_PsiCheckpointLoc * psi_checkpoint =
        malloc(sizeof(Blast_PsiCheckpointLoc));
    if (psi_checkpoint) {
        size_t length_filename = strlen(filename);
        psi_checkpoint->filename = calloc(length_filename + 1, sizeof(char));
        if (!psi_checkpoint->filename) {
            free(psi_checkpoint);
            psi_checkpoint = NULL;
        }
        memcpy(psi_checkpoint->filename, filename, length_filename + 1);
        psi_checkpoint->checkpoint_type = checkpoint_type;
    }
    return psi_checkpoint;
}


void
Blast_PsiCheckpointLocFree(Blast_PsiCheckpointLoc ** ppsi_checkpoint)
{
    Blast_PsiCheckpointLoc * psi_checkpoint = *ppsi_checkpoint;
    if (psi_checkpoint) {
        if (psi_checkpoint->filename) {
            free(psi_checkpoint->filename);
        }
        free(psi_checkpoint);
    }
    *ppsi_checkpoint = NULL;
}


int
Blast_PosReadCheckpoint(double ** freq_ratios,
                        int query_length,
                        const Uint1 * query,
                        Blast_PsiCheckpointLoc * psi_checkpoint,
                        Blast_Message* *blast_msg)
{
    int status;

    switch(psi_checkpoint->checkpoint_type) {
        case eStandardCheckpoint:
            status = s_PosReadStdCheckpoint(freq_ratios, query_length,
                                            query,
                                            psi_checkpoint->filename,
                                            blast_msg);
            break;
        case eAsnTextCheckpoint:
        case eAsnBinaryCheckpoint:
            {{
                int is_ascii_checkpoint =
                    psi_checkpoint->checkpoint_type == eAsnTextCheckpoint;
                status = s_PosReadAsnCheckpoint(freq_ratios,
                                                query_length, query,
                                                psi_checkpoint->filename,
                                                is_ascii_checkpoint);
            }}
            break;
        default:
            ASSERT(0 && "Impossible type of checkpoint file");
            status = -1;
            break;
    }
    return status;
}
/* @} */

