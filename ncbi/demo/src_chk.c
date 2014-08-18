/*   src_chk.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  src_chk.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   4/12/07
*
* $Revision: 1.10 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <pmfapi.h>
#define NLM_GENERATED_CODE_PROTO
#include <asnmacro.h>
#include <objmacro.h>
#include <macroapi.h>

#define SRC_CHK_APP_VER "1.0"

CharPtr SRC_CHK_APPLICATION = SRC_CHK_APP_VER;


static ValNodePtr CollectFieldList(BioseqPtr bsp)
{
  BioSourcePtr biop;
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  ValNodePtr list = NULL, vnp;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    vnp = GetSourceQualFieldListFromBioSource (biop);
    ValNodeLink (&list, vnp);
  }
  return list;
}


static void PrintHeader (FILE *fp, ValNodePtr field_list)
{
  CharPtr txt;

  if (fp == NULL || field_list == NULL) {
    return;
  }
  /* first field accession, second field GI, third field tax ID */
  fprintf (fp, "\t\tTaxID");
  while (field_list != NULL) {
    txt = SummarizeFieldType (field_list);
    fprintf (fp, "\t%s", txt);
    txt = MemFree (txt);
    field_list = field_list->next;
  }
  fprintf (fp, "\n");
}


static Int4 GetTaxIdFromOrgRef (OrgRefPtr orp)
{
  Int4       tax_id = -1;
  ValNodePtr vnp;
  DbtagPtr   d;

  if (orp != NULL)
  {
    for (vnp = orp->db; vnp != NULL; vnp = vnp->next) 
    {
      d = (DbtagPtr) vnp->data.ptrvalue;
      if (StringCmp(d->db, "taxon") == 0) 
      {
        tax_id = d->tag->id;
        break;
      }
    }
  }
  return tax_id;
}


static void PrintBioSourceLine (FILE *fp, BioSourcePtr biop, ValNodePtr field_list)
{
  CharPtr txt;

  if (fp == NULL || biop == NULL || field_list == NULL) {
    return;
  }

  fprintf (fp, "\t%d", GetTaxIdFromOrgRef(biop->org));
 
  while (field_list != NULL) {
    txt = GetSourceQualFromBioSource (biop, field_list->data.ptrvalue, NULL);
    fprintf (fp, "\t%s", txt == NULL ? "" : txt);
    txt = MemFree (txt);
    field_list = field_list->next;
  }
}


static void PrintBioseqLines (FILE *fp, BioseqPtr bsp, ValNodePtr field_list)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Char              id_txt[255], id_txt2[255];
  SeqIdPtr          sip, sip_gi = NULL, sip_gb = NULL;

  if (fp == NULL || bsp == NULL || field_list == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK
        || (sip->choice == SEQID_EMBL && sip_gb == NULL)
        || (sip->choice == SEQID_SWISSPROT && sip_gb == NULL)
        || (sip->choice == SEQID_DDBJ && sip_gb == NULL)
        || (sip->choice == SEQID_PIR && sip_gb == NULL)) {
      sip_gb = sip;
    } else if (sip->choice == SEQID_GI) {
      sip_gi = sip;
    }
  }

  if (sip_gb == NULL && sip_gi == NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    id_txt2[0] = 0;
  } else {
    if (sip_gb == NULL) {
      id_txt[0] = 0;
    } else {
      SeqIdWrite (sip_gb, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    }
    if (sip_gi == NULL) {
      id_txt2[0] = 0;
    } else {
      SeqIdWrite (sip_gi, id_txt2, PRINTID_REPORT, sizeof (id_txt2) - 1);
    }
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
    fprintf (fp, "%s\t%s", id_txt, id_txt2);
    PrintBioSourceLine (fp, sdp->data.ptrvalue, field_list);
    fprintf (fp, "\n");
  }
}


static void PrintBioseqErrorLine (FILE *fp, SeqIdPtr sip)
{
  Char              id_txt[255];

  if (fp == NULL || sip == NULL) {
    return;
  }

  SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);

  if (sip->choice == SEQID_GI) {
    fprintf (fp, "\t%s\n", id_txt);
  } else {
    fprintf (fp, "%s\t\n", id_txt);
  }
}


static Boolean IsAllDigits (CharPtr str)
{
  CharPtr cp;

  if (StringHasNoText (str)) return FALSE;

  cp = str;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (*cp == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static SeqIdPtr SmartGuessMakeId (CharPtr str)
{
  CharPtr id_txt;
  SeqIdPtr sip = NULL;

  if (StringHasNoText (str)) {
    return NULL;
  } else if (StringChr (str, '|') != NULL) {
    sip = MakeSeqID (str);
  } else if (IsAllDigits (str)) {
    id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 4));
    sprintf (id_txt, "gi|%s", str);
    sip = MakeSeqID (id_txt);
    id_txt = MemFree (id_txt);
  } else {
    id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 4));
    sprintf (id_txt, "gb|%s", str);
    sip = MakeSeqID (id_txt);
    id_txt = MemFree (id_txt);
  }
  return sip;
}


/* Args structure contains command-line arguments */

#define i_argInputFile         0
#define o_argOutputFile        1

Args myargs [] = {
  {"Input File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL}
};


static void SortFieldListForSrcChk (ValNodePtr PNTR field_list)
{
  ValNodePtr vnp, vnp_s, vnp_prev = NULL;

  if (field_list == NULL || *field_list == NULL) return;

  SortUniqueFieldTypeList (field_list);

  /* move taxname to front of list */
  for (vnp = *field_list; vnp != NULL; vnp_prev = vnp, vnp = vnp->next) {
    if (vnp->choice == FieldType_source_qual) {
      vnp_s = vnp->data.ptrvalue;
      if (vnp_s != NULL
          && vnp_s->choice == SourceQualChoice_textqual
          && vnp_s->data.intvalue == Source_qual_taxname) {
        /* only need to move if not already at front of list */
        if (vnp_prev != NULL) {
          vnp_prev->next = vnp->next;
          vnp->next = *field_list;
          *field_list = vnp;
        }
        break;
      }
    }
  }       


}


Int2 Main(void)
{
  Char             app [64];
  Int4             rval = 0;
  CharPtr          id_file, line;
  ReadBufferData   rbd;
  ValNodePtr       field_list = NULL;
  SeqIdPtr         sip;
  ValNodePtr       bsp_list = NULL, vnp;
  BioseqPtr        bsp;
  FILE *fp;


  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  PubSeqFetchEnable ();

  /* process command line arguments */

  sprintf (app, "src_chk %s", SRC_CHK_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  id_file = (CharPtr) myargs [i_argInputFile].strvalue;

  rbd.fp = FileOpen (id_file, "r");
  if (rbd.fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", (CharPtr) myargs [i_argInputFile].strvalue);
    return 1;
  }
  rbd.current_data = NULL;
  line = AbstractReadFunction (&rbd);  
  while (line != NULL && line[0] != EOF) {
    if (!StringHasNoText (line)) {

      sip = SmartGuessMakeId (line);
      bsp = BioseqLockById (sip);
      if (bsp == NULL) {
        printf ("Unable to download Bioseq for %s\n", line);
      } else {
        ValNodeLink (&field_list, CollectFieldList (bsp));
        BioseqUnlock (bsp);
      }
      ValNodeAddPointer (&bsp_list, 0, sip);
    }
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }

  FileClose (rbd.fp);

  SortFieldListForSrcChk (&field_list);

  fp = FileOpen ((CharPtr) myargs [o_argOutputFile].strvalue, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", (CharPtr) myargs [o_argOutputFile].strvalue);
    rval = 1;
  } else {
    PrintHeader (fp, field_list);
    for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
      bsp = BioseqLockById (vnp->data.ptrvalue);
      if (bsp == NULL) {
        PrintBioseqErrorLine (fp, vnp->data.ptrvalue);
      } else {
        PrintBioseqLines (fp, bsp, field_list);
      }
      BioseqUnlock (bsp);
      vnp->data.ptrvalue = SeqIdFree (vnp->data.ptrvalue);
    }
  }
  FileClose (fp);
  bsp_list = ValNodeFree (bsp_list);
  field_list = FieldTypeListFree (field_list);
  return rval;
}
