/*
 * $Id: acerdapi.c,v 1.15 2008/12/02 17:13:14 bollin Exp $
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


#include <stdlib.h>
#include <ncbi.h>
#include <ncbistr.h>
#include <seqport.h>
#include <sqnutils.h>
#include <gather.h>
#include <pmfapi.h>
#include <alignmgr2.h>
#include <explore.h>
#include <aceread.h>
#include <acerdapi.h>


/* This constructs an ASN.1 SeqGraph that contains the quality scores from the consensus sequence */
static SeqGraphPtr SeqGraphFromContig (TContigPtr contig, BioseqPtr bsp)
{
  SeqGraphPtr       sgp;
  ByteStorePtr      bs;
  Uint1             bytes[128]; 
  Int2              max = INT2_MIN;
  Int2              min = INT2_MAX;
  Int4              q_pos, b_pos;
  SeqIntPtr         sintp;

  if (contig == NULL || contig->num_qual_scores == 0 || contig->qual_scores == NULL
      || bsp == NULL) {
    return NULL;
  }

  sgp = SeqGraphNew ();
  bs = BSNew (1000);
  q_pos = 0;
  while (q_pos < contig->num_qual_scores) {
    b_pos = 0;
    while (b_pos < sizeof (bytes) && q_pos < contig->num_qual_scores) {
      max = MAX (max, (Int2) contig->qual_scores[q_pos]);
      min = MIN (min, (Int2) contig->qual_scores[q_pos]);
      bytes[b_pos++] = (Uint1) contig->qual_scores[q_pos++];
    }
    BSWrite (bs, (Pointer) bytes, (Int4) b_pos);
  }
  sgp->numval = BSLen (bs);
  BSPutByte (bs, EOF);
  sgp->title = StringSave ("Phrap Quality");
  sgp->flags [0] = 0;
  sgp->compr = 1;
  sgp->flags [1] = 0;
  sgp->flags [2] = 3;
  sgp->axis.intvalue = 0;
  sgp->min.intvalue = min;
  sgp->max.intvalue = max;
  sgp->a = 1.0;
  sgp->b = 0;
  sgp->values = (Pointer) bs;

  sintp = SeqIntNew ();
  sintp->from = 0;
  sintp->to = bsp->length - 1;
  sintp->id = SeqIdDup (bsp->id);
  ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);

  return sgp;
}


NLM_EXTERN SeqEntryPtr MakeSeqEntryFromRead (TContigReadPtr read)
{
  CharPtr seq_data;
  SeqIdPtr sip;
  SeqEntryPtr sep = NULL;
  BioseqPtr   bsp;
  SeqDescrPtr sdp;
  MolInfoPtr  mip;

  if (read == NULL) {
    return NULL;
  }

  seq_data = AlignmentStringToSequenceString (read->read_seq, Seq_mol_na);
  sip = MakeSeqID (read->read_id);
  sep = SequenceStringToSeqEntry (seq_data, sip, Seq_mol_na);
  if (sep != NULL && IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    bsp->mol = Seq_mol_rna;
    if (read->is_complement) {
      BioseqRevComp (bsp);
    }
    /* add molinfo */
    sdp = bsp->descr;
    while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
      sdp = sdp->next;
    }
    if (sdp == NULL) {
      sdp = SeqDescrNew (bsp->descr);
      if (bsp->descr == NULL) {
        bsp->descr = sdp;
      }
      sdp->choice = Seq_descr_molinfo;
      mip = MolInfoNew ();
      mip->biomol = MOLECULE_TYPE_MRNA;
      sdp->data.ptrvalue = mip;
    } else {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
    }
    mip->tech = MI_TECH_tsa;
  }
  return sep;
}


NLM_EXTERN SeqEntryPtr MakeSeqEntryFromContig (TContigPtr contig)
{
  CharPtr seq_data;
  SeqIdPtr sip;
  SeqEntryPtr sep = NULL;
  BioseqPtr   bsp;
  SeqGraphPtr sgp;
  SeqAnnotPtr sap;
  SeqDescrPtr sdp;
  MolInfoPtr  mip;

  if (contig == NULL) {
    return NULL;
  }

  seq_data = AlignmentStringToSequenceString (contig->consensus_seq, Seq_mol_na);
  sip = MakeSeqID (contig->consensus_id);
  sep = SequenceStringToSeqEntry (seq_data, sip, Seq_mol_na);
  if (sep != NULL && IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    bsp->mol = Seq_mol_rna;
    /* add molinfo */
    sdp = bsp->descr;
    while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
      sdp = sdp->next;
    }
    if (sdp == NULL) {
      sdp = SeqDescrNew (bsp->descr);
      if (bsp->descr == NULL) {
        bsp->descr = sdp;
      }
      sdp->choice = Seq_descr_molinfo;
      mip = MolInfoNew ();
      mip->biomol = MOLECULE_TYPE_MRNA;
      sdp->data.ptrvalue = mip;
    } else {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
    }
    mip->tech = MI_TECH_tsa;

    sgp = SeqGraphFromContig (contig, bsp);
    if (sgp != NULL) {
      sap = SeqAnnotNew ();
      sap->type = 3;
      sap->data = sgp;
      sap->next = bsp->annot;
      bsp->annot = sap;
    }
  }
  return sep;
}


/* This function compares a string of nucleotide characters to an existing Bioseq */
static Boolean DoesSeqStringMatchBsp (CharPtr seq_str, BioseqPtr bsp, Uint1 strand)
{
  Char buf[51];
  CharPtr cp_s, cp_b;
  Int4 ctr, pos = 0, i, len, seq_len;

  if (seq_str == NULL || bsp == NULL) return FALSE;
  cp_s = seq_str;
  len = sizeof (buf) - 1;
  seq_len = StringLen (seq_str);
  
  while (pos < bsp->length) {
    if (strand == Seq_strand_minus) {
      ctr = SeqPortStreamInt (bsp, MAX (0,  bsp->length - pos - len), bsp->length - pos - 1, Seq_strand_minus,
                            STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                            (Pointer) buf, NULL);
    } else {
      ctr = SeqPortStreamInt (bsp, pos, MIN(pos + len - 1, bsp->length - 1), Seq_strand_plus,
                            STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                            (Pointer) buf, NULL);
    }

    for (i = 0, cp_b = buf; i < ctr && *cp_s != 0; i++, cp_b++) { 
      while (*cp_s == '*') cp_s++;
      if (*cp_s != *cp_b) return FALSE;
      cp_s++;
    }
    if (ctr < len) {
      return TRUE;
    } else {
      pos = pos + len;
    }
  }
  if (*cp_s != 0) return FALSE;
  return TRUE;                          
}


static Int4 GetTraceID (SeqIdPtr sip)
{
  DbtagPtr dbtag;

  if (sip == NULL || sip->choice != SEQID_GENERAL) return 0;
  dbtag = (DbtagPtr) sip->data.ptrvalue;
  if (dbtag == NULL || StringCmp (dbtag->db, "ti") != 0 || dbtag->tag == NULL) {
    return 0;
  } 
  return dbtag->tag->id;
}


static Int4 GetTraceIDFromIdList (SeqIdPtr sip)
{
  Int4 ti = 0;

  while (sip != NULL && ti == 0) {
    ti = GetTraceID (sip);
    sip = sip->next;
  }
  return ti;
}

 

/* This function retrieves a sequence.  It would be better to use BioseqLockById. */
static SeqEntryPtr FetchRead (SeqIdPtr sip)
{
  Uint4       tid = 0;
  Int4        uid = 0;
  SeqEntryPtr sep = NULL;

  if (sip == NULL) return NULL;

  tid = GetTraceID (sip);
  if (tid > 0) {
    sep = PubSeqSynchronousQueryTI (tid, 0, -1);
  } else {
    uid = GetGIForSeqId (sip);
    if (uid > 0) {
      sep = PubSeqSynchronousQuery (uid, 0, -1);
    }
  }
  
  return sep;
}



static SeqIdPairPtr SeqIdPairNew ()
{
  SeqIdPairPtr pair;

  pair = (SeqIdPairPtr) MemNew (sizeof (SeqIdPairData));
  pair->sip_find = NULL;
  pair->sip_replace = NULL;
  return pair;
}


static SeqIdPairPtr SeqIdPairFree (SeqIdPairPtr pair)
{
  if (pair != NULL) {
    pair->sip_find = SeqIdFree (pair->sip_find);
    pair->sip_replace = SeqIdFree (pair->sip_replace);
    pair = MemFree (pair);
  }
  return pair;
}


static int SeqIdPairCompare (SeqIdPairPtr sp1, SeqIdPairPtr sp2)
{
  if (sp1 == NULL || sp2 == NULL) {
    return 0;
  }
  return StringICmp (sp1->buf_find, sp2->buf_find);
}


static int LIBCALLBACK SortSeqIdPairList (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      rval = SeqIdPairCompare (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return rval;
}


static ValNodePtr SeqIdPairListFree (ValNodePtr pair_list)
{
  ValNodePtr vnp_next;

  while (pair_list != NULL) {
    vnp_next = pair_list->next;
    pair_list->data.ptrvalue = SeqIdPairFree (pair_list->data.ptrvalue);
    pair_list->next = NULL;
    pair_list = ValNodeFree (pair_list);
    pair_list = vnp_next;
  }
  return pair_list;
}


static SeqIdReplaceListPtr SeqIdReplaceListNew (ValNodePtr id_list)
{
  SeqIdReplaceListPtr replace_list;
  SeqIdPairPtr        pair;
  Int4                i;

  replace_list = (SeqIdReplaceListPtr) MemNew (sizeof (SeqIdReplaceListData));
  replace_list->num_ids = ValNodeLen (id_list);
  replace_list->list = (SeqIdPairPtr) MemNew (sizeof (SeqIdPairData) * replace_list->num_ids);
  for (i = 0; id_list != NULL; id_list = id_list->next, i++) {
    pair = (SeqIdPairPtr) id_list->data.ptrvalue;
    replace_list->list[i].sip_find = SeqIdDup (pair->sip_find);
    StringCpy (replace_list->list[i].buf_find, pair->buf_find);
    replace_list->list[i].sip_replace = SeqIdDup (pair->sip_replace);
    replace_list->list[i].is_complement = pair->is_complement;
    replace_list->list[i].trim5 = pair->trim5;
    replace_list->list[i].trim3 = pair->trim3;
    replace_list->list[i].is_consensus = pair->is_consensus;
    replace_list->list[i].ti = pair->ti;
  }
  return replace_list;
}


NLM_EXTERN SeqIdReplaceListPtr SeqIdReplaceListFree (SeqIdReplaceListPtr replace_list)
{
  Int4 i;
  if (replace_list != NULL) {
    for (i = 0; i < replace_list->num_ids; i++) {
      replace_list->list[i].sip_find = SeqIdFree (replace_list->list[i].sip_find);
      replace_list->list[i].sip_replace = SeqIdFree (replace_list->list[i].sip_replace);
    }
    replace_list->list = MemFree (replace_list->list);
    replace_list = MemFree (replace_list);
  }
  return replace_list;
}


NLM_EXTERN SeqIdReplaceListPtr ReadSeqIdPairListFromFile (FILE *fp)
{
  ReadBufferData rbd;
  CharPtr        linestring, cp, id2, buf = NULL;
  Int4           len, buf_len = 0;
  SeqIdPairPtr   pair;
  ValNodePtr     pair_list = NULL, last = NULL, vnp;
  SeqIdReplaceListPtr replace_list = NULL;
  
  if (fp == NULL) return NULL;

  rbd.fp = fp;
  rbd.current_data = NULL;

  linestring = AbstractReadFunction (&rbd);
  while (linestring != NULL && linestring[0] != EOF) {
    cp = linestring + StringSpn (linestring, " \t");
    if (*cp != 0) {
      len = StringCSpn (cp, " \t");
      id2 = cp + len + StringSpn (cp + len, " \t");
      if (*id2 != 0) {
        if (len + 1 > buf_len) {
          buf = MemFree (buf);
          buf_len = len + 1;
          buf = (CharPtr) MemNew (sizeof (Char) * buf_len);
        }
        StringNCpy (buf, cp, len);
        buf[len] = 0;
        pair = SeqIdPairNew ();
        pair->sip_find = MakeSeqID (buf);
        SeqIdWrite (pair->sip_find, pair->buf_find, PRINTID_REPORT, sizeof (pair->buf_find) - 1);
        pair->sip_replace = MakeSeqID (id2);
        vnp = ValNodeNew (NULL);
        vnp->data.ptrvalue = pair;
        if (last == NULL) {
          pair_list = vnp;
        } else {
          last->next = vnp;
        }
        last = vnp;
      }
    }
    free (linestring);
    linestring = AbstractReadFunction (&rbd);     
  }
  pair_list = ValNodeSort (pair_list, SortSeqIdPairList);

  replace_list = SeqIdReplaceListNew (pair_list);
  pair_list = SeqIdPairListFree (pair_list);

  return replace_list;
}


static SeqIdPairPtr FindReplacementInSeqIdReplaceList (SeqIdPtr sip, SeqIdReplaceListPtr pair_list)
{
  Int4         l, r, m;
  Char         buf_find[100];
  int          cmp;

  if (sip == NULL || pair_list == NULL) return NULL;

  SeqIdWrite (sip, buf_find, PRINTID_REPORT, sizeof (buf_find) - 1);
  l = 0;
  r = pair_list->num_ids - 1;
  m = (r + l) / 2;

  while ((cmp = StringICmp (buf_find, pair_list->list[m].buf_find)) != 0 && l <= r) {
    if (cmp < 0) {
      r = m - 1;
    } else {
      l = m + 1;
    }
    m = (r + l) / 2;
  }
  if (cmp == 0) {
    return pair_list->list + m;
  } else {
    return NULL;
  }
}



static void ReportInvalidReplacement (SeqIdPtr sip, CharPtr reason, char *has_errors)
{
  Char         buf[128];

  SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
  PrintACEFormatErrorXMLStart (buf, has_errors);
  printf ("%s", reason);
  PrintACEFormatErrorXMLEnd ();
}


static Boolean OkToReplaceId (SeqIdPairPtr pair, CharPtr seq_str, char *has_errors)
{
  Boolean rval = FALSE;
  SeqEntryPtr fetched_sep, old_scope;
  BioseqPtr   bsp_replace;

  if (StringHasNoText (seq_str)) {
    rval = FALSE;
  }

  if (pair == NULL || pair->sip_replace == NULL) {
    rval = FALSE;
  } else if ((fetched_sep = FetchRead (pair->sip_replace)) == NULL) {
    rval = FALSE;
    ReportInvalidReplacement (pair->sip_replace, "Unable to fetch far sequence", has_errors);
  } else {
    old_scope = SeqEntrySetScope (fetched_sep);
    bsp_replace = BioseqFind (pair->sip_replace);
    SeqEntrySetScope (old_scope);
    if (bsp_replace == NULL) {
      rval = FALSE;
      ReportInvalidReplacement (pair->sip_replace, "Unable to locate far sequence after fetch", has_errors);
    } else if (DoesSeqStringMatchBsp (seq_str, bsp_replace, Seq_strand_plus)) {
      /* matches */
      rval = TRUE;
      pair->ti = GetTraceIDFromIdList (bsp_replace->id);
    } else if (DoesSeqStringMatchBsp (seq_str, bsp_replace, Seq_strand_minus)) {
      /* matches on complement */
      pair->is_complement = TRUE;
      rval = TRUE;
      pair->ti = GetTraceIDFromIdList (bsp_replace->id);
    } else {
      /* later, are we going to try to find trim lengths? */
      rval = FALSE;
      ReportInvalidReplacement (pair->sip_replace, "Replacement does not match local", has_errors);
    }
    SeqEntryFree (fetched_sep);
  }
  return rval;
}


static Boolean UpdateContigReadId (TContigReadPtr read, SeqIdReplaceListPtr pair_list, Boolean no_lookup, Boolean is_srr, char *has_errors)
{
  SeqIdPairPtr pair;
  SeqIdPtr     sip_find;
  Char         id_buf[255];
  Boolean      rval = TRUE;

  if (read == NULL || StringHasNoText (read->read_id)) {
    rval = FALSE;
  } else {
    sip_find = MakeSeqID (read->read_id);
    pair = FindReplacementInSeqIdReplaceList (sip_find, pair_list);
    if (pair != NULL && (no_lookup || OkToReplaceId (pair, read->read_seq, has_errors))) {
      if (pair->is_complement) {
        if (read->is_complement) {
          read->is_complement = FALSE;
        } else {
          read->is_complement = TRUE;
        }
      }
      if (pair->ti > 0) {
        read->ti = pair->ti;
      } else {
        if (pair->sip_replace->choice == SEQID_LOCAL) {
          SeqIdWrite (pair->sip_replace, id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
        } else {
          SeqIdWrite (pair->sip_replace, id_buf, PRINTID_FASTA_LONG, sizeof (id_buf) - 1);
        }
        if (is_srr) {
          if (read->srr != NULL) {
            free (read->srr);
          }
          read->srr = malloc (sizeof (Char) * (StringLen (id_buf) + 1));
          sprintf (read->srr, "%s", id_buf);
          free (read->read_id);
          read->read_id = NULL;
        } else {
          free (read->read_id);
          read->read_id = malloc (sizeof (Char) * (StringLen (id_buf) + 1));
          sprintf (read->read_id, "%s", id_buf);
        }
      }
      read->local = FALSE;
    }
    sip_find = SeqIdFree (sip_find);
  }
  return rval;
}


NLM_EXTERN Boolean UpdateContigIds (TContigPtr contig, SeqIdReplaceListPtr pair_list, Boolean no_lookup, Boolean is_srr, char *has_errors)
{
  Int4 i;
  SeqIdPairPtr pair;
  SeqIdPtr     sip_find;
  Char         id_buf[255];
  Boolean      rval = TRUE;

  if (contig == NULL) return FALSE;
  if (pair_list == NULL) return TRUE;

  if (contig->consensus_id != NULL) {
    sip_find = MakeSeqID (contig->consensus_id);
    pair = FindReplacementInSeqIdReplaceList (sip_find, pair_list);
    if (pair != NULL && (no_lookup || OkToReplaceId (pair, contig->consensus_seq, has_errors))) {
      if (pair->is_complement) {
        if (contig->is_complement) {
          contig->is_complement = FALSE;
        } else {
          contig->is_complement = TRUE;
        }
      }
      SeqIdWrite (pair->sip_replace, id_buf, PRINTID_FASTA_LONG, sizeof (id_buf) - 1);
      free (contig->consensus_id);
      contig->consensus_id = malloc (sizeof (Char) * (StringLen (id_buf) + 1));
      sprintf (contig->consensus_id, "%s", id_buf);
    } else {
      rval = FALSE;
    }
    sip_find = SeqIdFree (sip_find);
  }
  for (i = 0; i < contig->num_reads; i++) {
    rval &= UpdateContigReadId (contig->reads[i], pair_list, no_lookup, is_srr, has_errors);
  }
  return rval;
}


NLM_EXTERN Boolean UpdateAceFileIds (TACEFilePtr afp, FILE *id_file, Boolean no_lookup, Boolean is_srr, char *has_errors)
{
  Boolean    rval = TRUE;
  SeqIdReplaceListPtr pair_list;
  SeqEntryPtr old_scope;
  Int4        i;

  if (afp == NULL || id_file == NULL) return FALSE;
  old_scope = SeqEntrySetScope (NULL);
  pair_list = ReadSeqIdPairListFromFile (id_file);
  for (i = 0; i < afp->num_contigs; i++) {
    rval &= UpdateContigIds (afp->contigs[i], pair_list, no_lookup, is_srr, has_errors);
  }  
  
  pair_list = SeqIdReplaceListFree (pair_list);
  SeqEntrySetScope (old_scope);
  return rval; 
}


static Boolean ValidateContigReadId (TContigReadPtr read, char *has_errors)
{
  SeqIdPairData pair;
  Char          id_buf[255];
  Boolean       rval = TRUE;

  if (read == NULL || StringHasNoText (read->read_id)) {
    rval = FALSE;
  } else if (!read->local) {
    rval = TRUE;
  } else {
    pair.sip_find = NULL;
    pair.is_complement = FALSE;
    pair.is_consensus = FALSE;
    pair.trim3 = 0;
    pair.trim5 = 0;
    pair.sip_replace = MakeSeqID (read->read_id);
    pair.ti = 0;
    if (OkToReplaceId (&pair, read->read_seq, has_errors)) {
      if (pair.is_complement && !read->is_complement) {
        read->is_complement = TRUE;
      } else if (!pair.is_complement && read->is_complement) {
        read->is_complement = FALSE;
      }
      if (pair.ti > 0) {
        read->ti = pair.ti;
      } else {
        SeqIdWrite (pair.sip_replace, id_buf, PRINTID_FASTA_LONG, sizeof (id_buf) - 1);
        free (read->read_id);
        read->read_id = malloc (sizeof (Char) * (StringLen (id_buf) + 1));
        sprintf (read->read_id, "%s", id_buf);
      }
      read->local = FALSE;
    }
    pair.sip_replace = SeqIdFree (pair.sip_replace);
  }
  return rval;
}


static Boolean ValidateContigIds (TContigPtr contig, char *has_errors)
{
  Int4 i;
  Boolean      rval = TRUE;

  if (contig == NULL) return FALSE;

  if (contig->consensus_id != NULL) {
    /* check consensus later... */
  }
  for (i = 0; i < contig->num_reads; i++) {
    rval &= ValidateContigReadId (contig->reads[i], has_errors);
  }
  return rval;
}


NLM_EXTERN Boolean ValidateAceFileIds (TACEFilePtr afp, char *has_errors)
{
  Boolean    rval = TRUE;
  SeqEntryPtr old_scope;
  Int4        i;

  if (afp == NULL) return FALSE;
  old_scope = SeqEntrySetScope (NULL);
  for (i = 0; i < afp->num_contigs; i++) {
    rval &= ValidateContigIds (afp->contigs[i], has_errors);
  }  
  
  SeqEntrySetScope (old_scope);
  return rval; 
}


NLM_EXTERN ValNodePtr GetTransitionsFromGapInfo (TGapInfoPtr gaps, Int4 offset, Int4 seq_offset, Int4 seq_len)
{
  ValNodePtr list = NULL;
  Int4 i = 0, tiling_pos = offset, seq_pos = 0, diff = 0;
  Boolean added_gap = FALSE;

  /* add a transition to the list for where a sequence "begins" in the alignment, if not at 0 */
  if (seq_offset == 0) {
    if (tiling_pos > 0) {
      ValNodeAddInt (&list, 0, tiling_pos);
    }
  } else {
    /* if seq_offset causes sequence to "start" in the middle of a between-gap interval, add a transition for it */
    if (gaps == NULL || gaps->num_gaps == 0) {
      ValNodeAddInt (&list, 0, tiling_pos + seq_offset);
    } else {
      while (seq_pos < seq_offset && i < gaps->num_gaps && !added_gap) {
        if (seq_pos + gaps->gap_offsets[i] <= seq_offset) {
          tiling_pos += gaps->gap_offsets[i] + 1;
          seq_pos += gaps->gap_offsets[i];
          diff += gaps->gap_offsets[i];
          i++;
        } else {
          ValNodeAddInt (&list, 0, tiling_pos + seq_offset);
          added_gap = TRUE;
        }
      }
    }
  }

  if (gaps != NULL) {
    while (i < gaps->num_gaps) {
      seq_pos += gaps->gap_offsets[i];
      if (gaps->gap_offsets[i] > 0) {
        tiling_pos += gaps->gap_offsets[i];
        ValNodeAddInt (&list, 0, tiling_pos);
      }
      tiling_pos++;
      if (gaps->num_gaps == i + 1
          || gaps->gap_offsets[i + 1] > 0) {
        ValNodeAddInt (&list, 0, tiling_pos);
      }
      i++;
    }
  }
  if (seq_pos < seq_len) {
    ValNodeAddInt (&list, 0, tiling_pos + seq_len - seq_pos);
  }
  return list;
}


static Boolean ValidateContigAgainstSeqEntry (TContigPtr contig, SeqEntryPtr sep, char *has_errors)
{
  CharPtr      seq_data = NULL;
  SeqIdPtr     sip;
  BioseqPtr    bsp;
  Boolean      rval = FALSE;

  if (contig == NULL || sep == NULL) {
    return FALSE;
  }

  seq_data = AlignmentStringToSequenceString (contig->consensus_seq, Seq_mol_na);
  sip = MakeSeqID (contig->consensus_id);
  
  bsp = BioseqFind (sip);
  if (bsp == NULL) {
    PrintACEFormatErrorXML ("not found in supplied SeqEntry", contig->consensus_id, has_errors);
  } else if (!DoesSeqStringMatchBsp (seq_data, bsp, Seq_strand_plus)) {
    PrintACEFormatErrorXML ("does not match sequence in supplied SeqEntry", contig->consensus_id, has_errors);
  } else {
    rval = TRUE;
  }
  seq_data = MemFree (seq_data);
  return rval;
}


NLM_EXTERN Boolean ValidateACEFileAgainstSeqEntry (TACEFilePtr ace_file, SeqEntryPtr sep, char *has_errors)
{
  Boolean rval = TRUE;
  Int4    i;
  SeqEntryPtr oldscope;

  if (ace_file == NULL || sep == NULL) {
    return FALSE;
  }

  oldscope = SeqEntrySetScope (sep);

  for (i = 0; i < ace_file->num_contigs; i++) {
    rval |= ValidateContigAgainstSeqEntry (ace_file->contigs[i], sep, has_errors);
  }
  SeqEntrySetScope (oldscope);
  return rval;
}

