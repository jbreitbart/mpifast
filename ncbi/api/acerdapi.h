#ifndef API_ACERDAPI__H
#define API_ACERDAPI__H

/*
 * $Id: acerdapi.h,v 1.8 2008/12/02 16:16:01 bollin Exp $
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


#ifdef __cplusplus
extern "C" {
#endif

/* for doing ID replacements */
typedef struct seqidpair {
  SeqIdPtr sip_find;
  Char     buf_find[100];
  SeqIdPtr sip_replace;
  Boolean  is_complement;
  Int4     trim5;
  Int4     trim3;
  Boolean  is_consensus;
  Int4     ti;
} SeqIdPairData, PNTR SeqIdPairPtr;

typedef struct seqidreplacelist {
  Int4 num_ids;
  SeqIdPairPtr list;
} SeqIdReplaceListData, PNTR SeqIdReplaceListPtr;

NLM_EXTERN SeqIdReplaceListPtr SeqIdReplaceListFree (SeqIdReplaceListPtr replace_list);
NLM_EXTERN SeqIdReplaceListPtr ReadSeqIdPairListFromFile (FILE *fp);


NLM_EXTERN Boolean UpdateContigIds (TContigPtr contig, SeqIdReplaceListPtr pair_list, Boolean no_lookup, Boolean is_srr, char *has_errors);
NLM_EXTERN Boolean UpdateAceFileIds (TACEFilePtr afp, FILE *id_file, Boolean no_lookup, Boolean is_srr, char *has_errors);
NLM_EXTERN Boolean ValidateAceFileIds (TACEFilePtr afp, char *has_errors);

NLM_EXTERN ValNodePtr GetTransitionsFromGapInfo (TGapInfoPtr gaps, Int4 offset, Int4 seq_offset, Int4 tiling_len);

NLM_EXTERN SeqEntryPtr MakeSeqEntryFromRead (TContigReadPtr read);
NLM_EXTERN SeqEntryPtr MakeSeqEntryFromContig (TContigPtr contig);

NLM_EXTERN Boolean ValidateACEFileAgainstSeqEntry (TACEFilePtr ace_file, SeqEntryPtr sep, char *has_errors);

#ifdef __cplusplus
}
#endif

/*
 * ==========================================================================
 *
 * $Log: acerdapi.h,v $
 * Revision 1.8  2008/12/02 16:16:01  bollin
 * Added function for making a SeqEntry from a read, fixed bug in seqid lookup.
 *
 * Revision 1.7  2008/11/26 18:30:02  bollin
 * Changes to make aceread_tst more efficient when handling large ACE files,
 * added TSA field tags for assembly and taxid.
 *
 * Revision 1.6  2008/08/13 18:18:28  bollin
 * Improved aceread_tst - adds molinfo and SeqGraph with quality scores to generated ASN.1
 *
 * Revision 1.5  2008/08/13 15:35:30  bollin
 * Added wrapping header for XML errors during ACE read.
 *
 * Revision 1.4  2008/08/13 12:30:01  bollin
 * Changes to allow use of srr numbers in XML and suppress lookups.  Also fixes segfault.
 *
 * Revision 1.3  2008/07/22 20:25:49  bollin
 * Added functiosn for validating consensus sequences against a supplied SeqEntry
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

#endif

