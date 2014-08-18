/*   sequin3.c
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
* File Name:  sequin3.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.960 $
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

#include "sequin.h"
#include <document.h>
#include <sequtil.h>
#include <biosrc.h>
#include <import.h>
#include <gather.h>
#include <asn2gnbk.h>
#include <asn2gnbp.h>
#include <edutil.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <subutil.h>    /* TOPOLOGY_xxx definitions */
#include <salutil.h>
#include <valid.h>
#include <vsm.h>
#include <bspview.h>
#include <toasn3.h>
#include <salstruc.h>
#include <explore.h>
#include <utilpub.h>
#include <tofasta.h>
#include <rpsutil.h>
#include <findrepl.h>
#include <explore.h>
#include <seqpanel.h>
#include <tax3api.h>
#include <saledit.h>
#include <alignmgr2.h>
#include <suggslp.h>
#include <asn2gnbi.h>
#include <alignval.h>
#include <cdrgn.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macrodlg.h>
#include <macroapi.h>

/* For converting primer names to primer seqs and vice versa */
static void PrimerSeqToPrimerNameCallback (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr ssp;

  if (biop == NULL) return;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      ssp->subtype = SUBSRC_fwd_primer_name;
    } else if (ssp->subtype == SUBSRC_rev_primer_seq) {
      ssp->subtype = SUBSRC_rev_primer_name;
    }
  }
}


static void PrimerNameToPrimerSeqCallback (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr ssp;

  if (biop == NULL) return;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_fwd_primer_name) {
      ssp->subtype = SUBSRC_fwd_primer_seq;
    } else if (ssp->subtype == SUBSRC_rev_primer_name) {
      ssp->subtype = SUBSRC_rev_primer_seq;
    }
  }
}


static void SwapPrimerNamesAndSeq (IteM i, Boolean seq_to_name)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (seq_to_name) {
    VisitBioSourcesInSep (sep, NULL, PrimerSeqToPrimerNameCallback);
  } else {
    VisitBioSourcesInSep (sep, NULL, PrimerNameToPrimerSeqCallback);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void ConvertPrimerSeqToName (IteM i)
{
  SwapPrimerNamesAndSeq (i, TRUE);
}


static void ConvertPrimerNameToSeq (IteM i)
{
  SwapPrimerNamesAndSeq (i, FALSE);
}


static void TSATableCallback (Pointer data)
{
  ProcessExternalEvent ();
}

static ValNodePtr ApplyTranscriptomeIdListWithProgress (ValNodePtr ids_list, LocalAlignFunc aln_func)
{
  WindoW     w;
  GrouP      h;
  PrompT     ppt;
  CharPtr    prompt_fmt = "%d out of %d";
  Char       prompt_str[40];
  Int4       max, num = 1;
  ValNodePtr err_list = NULL, vnp;

  if (ids_list == NULL || aln_func == NULL) return NULL;

  WatchCursor();
  Update();
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  max = ValNodeLen (ids_list);
  sprintf (prompt_str, prompt_fmt, max, max);
  ppt = StaticPrompt (h, prompt_str, 0, dialogTextHeight, programFont, 'c');
  sprintf (prompt_str, prompt_fmt, 1, max);
  SetTitle (ppt, prompt_str);  
  
  Show(w); 
  Select (w);
  for (num = 1, vnp = ids_list; vnp != NULL; vnp = vnp->next, num++) {
    sprintf (prompt_str, prompt_fmt, num, max);
    SetTitle (ppt, prompt_str);  
    ValNodeLink (&err_list, MakeTranscriptomeAssemblySeqHist (vnp->data.ptrvalue, aln_func, TSATableCallback, NULL));
  }

  Remove (w);
  ArrowCursor();
  Update();
  return err_list;
}


static void AddTSATableToBioseq (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ValNodePtr   err_list = NULL, ids_list, vnp;
  LogInfoPtr   lip;
  Char         path [PATH_MAX];
  FILE        *fp;
  MsgAnswer    ans;
  ValNodePtr   coverage_report;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }

  ids_list = GetTranscriptomeIdsList (fp, sep, &err_list);
  FileClose (fp);

  if (err_list != NULL) {
    lip = OpenLog ("Problems Reading TSA Table File");
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", vnp->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
    err_list = ValNodeFreeData (err_list);
    if (ANS_NO == Message (MSG_YN, "Problems reading input file - continue creating TSA tables?")) {
      ids_list = TranscriptomeIdsListFree (ids_list);
      return;
    }
  }

  if (HasExistingSeqHistAssembly(ids_list)) {
    ans = Message (MSG_YNC, "Some of the bioseqs listed in the file already have Seq-Hist assemblies.  Should these be replaced?");
    if (ans == ANS_YES) {
      DeleteSeqHistAssembliesForList(ids_list);
    } else if (ans == ANS_CANCEL) {
      ids_list = TranscriptomeIdsListFree (ids_list);
      return;
    }
  }

  
  ValNodeLink (&err_list, ApplyTranscriptomeIdListWithProgress (ids_list, GetSeqAlignTSA));

  coverage_report = ReportCoverageForTranscriptomeIdsListSeqHist (ids_list);

  ids_list = TranscriptomeIdsListFree (ids_list);

  ValNodeLink (&coverage_report, err_list);
  err_list = coverage_report;

  if (err_list != NULL) {
    lip = OpenLog ("TSA Table Problems");
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", vnp->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
    err_list = ValNodeFreeData (err_list);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void RefreshTSATables (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ValNodePtr   err_list = NULL, ids_list, vnp;
  LogInfoPtr   lip;
  ValNodePtr   coverage_report;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;


  ids_list = GetExistingTSATableIds (sep);

  DeleteSeqHistAssembliesForList(ids_list);
  
  ValNodeLink (&err_list, ApplyTranscriptomeIdListWithProgress (ids_list, GetSeqAlignTSA));

  coverage_report = ReportCoverageForTranscriptomeIdsListSeqHist (ids_list);

  ids_list = TranscriptomeIdsListFree (ids_list);

  ValNodeLink (&coverage_report, err_list);
  err_list = coverage_report;

  if (err_list != NULL) {
    lip = OpenLog ("TSA Table Problems");
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", vnp->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
    err_list = ValNodeFreeData (err_list);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void CreateTSAIDsFromLocalIDs (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ConvertLocalIdsToTSAIds (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void AddmRNASequenceCallback (SeqDescrPtr sdp, Pointer data)
{
  Int4 len, len_add;
  CharPtr add = ", mRNA sequence.", tmp, title;

  if (sdp != NULL && sdp->choice == Seq_descr_title) {
    title = (CharPtr) sdp->data.ptrvalue;
    len_add = StringLen (add);
    len = StringLen (title);
    if (len < len_add || StringCmp (title + (len - len_add), add) != 0) {
      /* remove trailing period if present */
      if (title[len - 1] == '.') {
        title[len - 1] = 0;
      }
      tmp = (CharPtr) MemNew (sizeof (Char) * (len + len_add + 1));
      sprintf (tmp, "%s%s", title, add);
      sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
      sdp->data.ptrvalue = tmp;
    }
  }
}


static void AddmRNASequenceToDeflines (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  AutoDefEntityIDNoOptions (bfp->input_entityID, TRUE);
  VisitDescriptorsInSep (sep, NULL, AddmRNASequenceCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void SplitPubsByTable (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  FILE *fp;
  ValNodePtr   tab_table, split_list;
  ValNodePtr   err_list = NULL;
  Char         path [PATH_MAX];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;
  tab_table = ReadTabTableFromFile (fp);
  FileClose (fp);
  split_list = MakeSplitPubListFromTabList (&tab_table, sep, &err_list);
  tab_table = FreeTabTable (tab_table);
  SplitPubsByList (split_list);
  SplitPubListFree (split_list);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update();  
}


#define DDBUFFERMAX 32768

typedef struct fndeldata {
  CharPtr     buffer;
  CharPtr     curr;
  Int4        seqlen;
  Int4        gap_count;
  Int4        bsp_length;
  Boolean     in_gap;
  Boolean     unk_gap;
  ValNodePtr  head;
  ValNodePtr  last;
} DelData, PNTR DelDataPtr;

static void AddNextBlock (
  DelDataPtr ddp
)

{
  ByteStorePtr  bs;
  IntFuzzPtr    ifp;
  SeqLitPtr     slitp;
  CharPtr       str;
  ValNodePtr    vnp;

  if (ddp == NULL) return;
  if (ddp->in_gap) {
    slitp = (SeqLitPtr) MemNew (sizeof (SeqLit));
    if (slitp != NULL) {
      slitp->length = ddp->gap_count;
      ddp->bsp_length += slitp->length;
      if (ddp->unk_gap) {
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          slitp->fuzz = ifp;
        }
      }
      vnp = ValNodeAddPointer (&(ddp->last), 2, (Pointer) slitp);
      if (ddp->head == NULL) {
        ddp->head = vnp;
      }
      ddp->last = vnp;
    }
  } else if (ddp->seqlen > 0) {
    str = ddp->curr;
    if (str != NULL) {
      *str = '\0';
    }
    bs = BSNew (ddp->seqlen);
    if (bs != NULL) {
      BSWrite (bs, (VoidPtr) ddp->buffer, ddp->seqlen);
      slitp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slitp != NULL) {
        slitp->length = ddp->seqlen;
        slitp->seq_data_type = Seq_code_iupacna;
        slitp->seq_data = (SeqDataPtr) bs;
        ddp->bsp_length += slitp->length;
        vnp = ValNodeAddPointer (&(ddp->last), 2, (Pointer) slitp);
        if (ddp->head == NULL) {
          ddp->head = vnp;
        }
        ddp->last = vnp;
      }
    }
  }
}

static void LIBCALLBACK FarToNearCallback (
  CharPtr sequence,
  Pointer userdata
)

{
  Char        ch;
  DelDataPtr  ddp;
  CharPtr     ptr;
  CharPtr     str;

  if (sequence == NULL || userdata == NULL) return;
  ddp = (DelDataPtr) userdata;

  ptr = sequence;
  ch = *ptr;
  str = ddp->curr;

  while (ch != '\0') {
    if (ch != '-' && ch != '+') {
      if (ddp->in_gap) {
        AddNextBlock (ddp);
        ddp->in_gap = FALSE;
        ddp->unk_gap = FALSE;
        ddp->gap_count = 0;
        ddp->curr = ddp->buffer;
        ddp->seqlen = 0;
        str = ddp->curr;
      }
      *str = ch;
      str++;
      (ddp->seqlen)++;
      if (ddp->seqlen >= DDBUFFERMAX) {
        *str = '\0';
        ddp->curr = str;
        AddNextBlock (ddp);
        ddp->curr = ddp->buffer;
        ddp->seqlen = 0;
        str = ddp->curr;
      }
    } else {
      if (! ddp->in_gap) {
        AddNextBlock (ddp);
        ddp->in_gap = TRUE;
        ddp->unk_gap = (Boolean) (ch == '-');
        ddp->gap_count = 0;
        ddp->curr = ddp->buffer;
        ddp->seqlen = 0;
        str = ddp->curr;
      }
      (ddp->gap_count)++;
    }

    ptr++;
    ch = *ptr;
  }
  ddp->curr = str;
}

static void FarToNearProc (BioseqPtr bsp, Pointer userdata)

{
  DelData  dd;
  Pointer  olddelta;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;
  MemSet ((Pointer) &dd, 0, sizeof (DelData));
  dd.buffer = MemNew (sizeof (Char) * (DDBUFFERMAX + 2));
  if (dd.buffer == NULL) return;
  dd.curr = dd.buffer;
  dd.seqlen = 0;
  SeqPortStream (bsp, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, (Pointer) &dd, FarToNearCallback);
  AddNextBlock (&dd);
  MemFree (dd.buffer);
  olddelta = bsp->seq_ext;
  bsp->seq_ext = (Pointer) dd.head;
  if (dd.bsp_length != bsp->length) {
    Message (MSG_OK, "Old length %ld differs from new length %ld",
             (long) bsp->length, (long) dd.bsp_length);
  }
}

static void FarDeltaToNear (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, FarToNearProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Message (MSG_OK, "Some manual desktop manipulations remain");
}

static CharPtr HIVDBStart = "##HIVDataBaseData-START##";
static CharPtr HIVDBEnd = "##HIVDataBaseData-END##";

static CharPtr HIVStart = "##HIVData-START##";
static CharPtr HIVEnd = "##HIVData-END##";

static CharPtr FluStart = "##FluData-START##";
static CharPtr FluEnd = "##FluData-END##";


static void RemoveStructuredCommentFromString (CharPtr str, CharPtr prefix, CharPtr suffix)
{
  CharPtr pPrefix, pSuffix, pAfter;

  if (StringHasNoText (str) || StringHasNoText (prefix) || StringHasNoText (suffix)) return;

  pPrefix = StringSearch (str, prefix);
  pSuffix = StringSearch (str, suffix);
  if (pPrefix == NULL || pSuffix == NULL || pPrefix >= pSuffix) return;

  pAfter = pSuffix + StringLen (suffix);
  StringCpy (pPrefix, pAfter);
}


static void FindStructComment (SeqDescrPtr sdp, Pointer userdata)

{
  CharPtr        str;
  UserObjectPtr  uop = NULL;
  UserObjectPtr  PNTR uopp;
  CharPtr        prefix = "", suffix = "";

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;
  uopp = (UserObjectPtr PNTR) userdata;
  if (uopp == NULL) return;

  if (StringStr (str, HIVDBStart) != NULL) {
    prefix = HIVDBStart;
    suffix = HIVDBEnd;
  } else if (StringStr (str, HIVStart) != NULL) {
    prefix = HIVStart;
    suffix = HIVEnd;
  } else if (StringStr (str, FluStart) != NULL) {
    prefix = FluStart;
    prefix = FluEnd;
  }

  uop = ParseStringIntoStructuredComment (NULL, str, prefix, suffix);

  if (uop != NULL) {
    *uopp = uop;
    RemoveStructuredCommentFromString (str, prefix, suffix);
    if (StringHasNoText (str) && sdp->extended != 0) {
      ((ObjValNodePtr)sdp)->idx.deleteme = TRUE;
    }
  }
}

static void DoStructCommentBsp (BioseqPtr bsp, Pointer userdata)

{
  SeqDescrPtr    sdp;
  UserObjectPtr  uop = NULL;

  VisitDescriptorsOnBsp (bsp, (Pointer) &uop, FindStructComment);
  if (uop == NULL) return;
  sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = uop;
}

static void ParseCommentIntoStructuredObject (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, DoStructCommentBsp);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void DoHivOrFluTag (SeqDescrPtr sdp, Pointer userdata, CharPtr prefix, CharPtr suffix)

{
  UserFieldPtr   curr;
  Boolean        hasPrefix = FALSE;
  Boolean        hasSuffix = FALSE;
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "StructuredComment") != 0) return;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "StructuredCommentPrefix") == 0) {
      hasPrefix = TRUE;
      if (curr->choice == 1) {
        MemFree (curr->data.ptrvalue);
        curr->data.ptrvalue = (Pointer) StringSave (prefix);
      }
    }
  }
  if (! hasPrefix) {
    AddItemStructuredCommentUserObject (uop, "StructuredCommentPrefix", prefix);
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "StructuredCommentSuffix") == 0) {
      hasSuffix = TRUE;
      if (curr->choice == 1) {
        MemFree (curr->data.ptrvalue);
        curr->data.ptrvalue = (Pointer) StringSave (suffix);
      }
    }
  }
  if (! hasSuffix) {
    AddItemStructuredCommentUserObject (uop, "StructuredCommentSuffix", suffix);
  }
}

static void DoHivTag (SeqDescrPtr sdp, Pointer userdata)

{
  DoHivOrFluTag (sdp, userdata, "##HIVData-START##", "##HIVData-END##");
}

static void DoFluTag (SeqDescrPtr sdp, Pointer userdata)

{
  DoHivOrFluTag (sdp, userdata, "##FluData-START##", "##FluData-END##");
}

static void AddStructuredCommentHivTag (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitDescriptorsInSep (sep, NULL, DoHivTag);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void AddStructuredCommentFluTag (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitDescriptorsInSep (sep, NULL, DoFluTag);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void ConvertBadInfProc (SeqFeatPtr sfp, Pointer userdata)

{
  GBQualPtr  gbq;
  size_t     len;
  CharPtr    ptr, str;

  if (sfp == NULL) return;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    if (ValidateInferenceQualifier (gbq->val, FALSE) != VALID_INFERENCE) {
      if (StringNICmp (gbq->val, "similar to ", 11) == 0) {
        ptr = StringChr (gbq->val, ':');
        if (ptr != NULL) {
          ptr++;
          if (StringDoesHaveText (ptr)) {
            len = StringLen ("similar to ") + StringLen (ptr);
            str = MemNew (len + 5);
            if (str != NULL) {
              StringCpy (str, "similar to ");
              StringCat (str, ptr);
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave (str);
              MemFree (str);
            }
          }
        }
      }
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("note");
    }
  }
}

static void ConvertBadInferenceToNote (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  VisitFeaturesInSep (sep, NULL, ConvertBadInfProc);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update (); 
}

typedef struct {
  Int2    level;
  Boolean prevWasModified;
} VecToNData, PNTR VecToNDataPtr;

typedef struct {
  Boolean dirtyFlag;
  Int2    number;
} CdsDelStruct, PNTR CdsDelStructPtr;

typedef struct {
  FORM_MESSAGE_BLOCK
  Int2  number;
} CdsDelForm, PNTR CdsDelFormPtr;

static void AddOrgNameToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, 0, 0, NULL);
}

static void AddStrainToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, ORGMOD_strain, 0, NULL);
}

static void AddCloneToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, 0, SUBSRC_clone, NULL);
}

static void AddIsolateToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, ORGMOD_isolate, 0, NULL);
}

static void AddHaplotypeToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, 0, SUBSRC_haplotype, NULL);
}

static void AddCultivarToDefLines (IteM i)

{
  CommonAddOrgOrModsToDefLines (i, ORGMOD_cultivar, 0, NULL);
}

static Boolean FindBestCitSubCallback (GatherContextPtr gcp)

{
  CitSubPtr   best;
  CitSubPtr   PNTR bestp;
  CitSubPtr   csp;
  PubdescPtr  pdp;
  ValNodePtr  sdp;
  ValNodePtr  vnp;

  if (gcp == NULL) return TRUE;
  bestp = (CitSubPtr PNTR) gcp->userdata;
  if (bestp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQDESC) return TRUE;
  sdp = (ValNodePtr) gcp->thisitem;
  if (sdp == NULL || sdp->choice != Seq_descr_pub) return TRUE;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return TRUE;
  vnp = pdp->pub;
  if (vnp == NULL || vnp->choice != PUB_Sub) return TRUE;
  csp = (CitSubPtr) vnp->data.ptrvalue;
  if (csp == NULL) return TRUE;
  if (*bestp == NULL) {
    *bestp = csp;
    return TRUE;
  }
  best = *bestp;
  if (DateMatch (best->date, csp->date, FALSE) == -1) {
    *bestp = csp;
    return TRUE;
  }
  return TRUE;
}

CharPtr kSubmitterUpdateText = "Sequence update by submitter";
CharPtr kIndexerUpdateVecScreenText = "Sequence update by database staff to remove vector contamination";

extern Boolean CreateUpdateCitSubFromBestTemplate (
  SeqEntryPtr top_sep,
  SeqEntryPtr upd_sep,
  CharPtr     update_txt
)
{
  CitSubPtr    best;
  CitSubPtr    csp;
  DatePtr      dp;
  GatherScope  gs;
  PubdescPtr   pdp;
  ValNodePtr   sdp;
  ValNodePtr   vnp;

  best = NULL;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.scope = top_sep;
  GatherSeqEntry (top_sep, (Pointer) (&best), FindBestCitSubCallback, &gs);
  if (best == NULL) {
    Message (MSG_OK, "There is no earlier cit-sub template");
    return FALSE;
  }
  dp = DateCurr ();
  if (dp != NULL) {
    if (DateMatch (best->date, dp, FALSE) == 0) {
      DateFree (dp);
	    dp = NULL;
      if (upd_sep == top_sep
          && StringICmp (best->descr, update_txt) == 0)
      {
        Message (MSG_OK, "There already exists an update on today's date");
        Update ();
        return FALSE;
      } else if (best->descr == NULL) {
        best->descr = StringSave (update_txt);
        Message (MSG_OK, "Adding update indication to existing cit-sub");
        return TRUE;
      }
    }
    DateFree (dp);
  }
  sdp = CreateNewDescriptor (upd_sep, Seq_descr_pub);
  if (sdp == NULL) return FALSE;
  pdp = PubdescNew ();
  if (pdp == NULL) return FALSE;
  sdp->data.ptrvalue = (Pointer) pdp;
  vnp = ValNodeNew (NULL);
  csp = AsnIoMemCopy ((Pointer) best,
                      (AsnReadFunc) CitSubAsnRead,
                      (AsnWriteFunc) CitSubAsnWrite);
  pdp->pub = vnp;
  vnp->choice = PUB_Sub;
  vnp->data.ptrvalue = csp;
  csp->date = DateFree (csp->date);
  csp->date = DateCurr ();
  if (!StringHasNoText (update_txt)) {
    csp->descr = MemFree (csp->descr);
    csp->descr = StringSave (update_txt);
  }

  if (top_sep == upd_sep)
  {
    Message (MSG_OK, "The update Cit-sub has been placed on the top Seq-entry");
  }
  return TRUE;
}                               

static void AddCitSubForUpdateProc (BaseFormPtr bfp)

{
  SeqEntryPtr  sep;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if ( CreateUpdateCitSubFromBestTemplate (sep, sep, kSubmitterUpdateText))
  {
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    Update ();
  }
  return;
}

static void AddCitSubForUpdate (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  AddCitSubForUpdateProc (bfp);
}

/* FixLocus is copied from taxutil, then modified to remove embl_code use */
static void FixLocus (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr bsp;
	ValNodePtr vnp;
    TextSeqIdPtr tsip;
    Char tbuf[40];
	CharPtr ptr, embl_code;

	if (! IS_Bioseq(sep))
		return;

	bsp = (BioseqPtr)(sep->data.ptrvalue);
    embl_code = (CharPtr) data;
	/* if (! embl_code) return; */

    for (vnp = bsp->id; vnp != NULL; vnp = vnp->next) {
		if ((vnp->choice == SEQID_GENBANK || vnp->choice == SEQID_TPG) && vnp->data.ptrvalue != NULL) {
			tsip = (TextSeqIdPtr) (vnp->data.ptrvalue);
			if (tsip->accession != NULL)
			{
				tsip->name = MemFree (tsip->name);
				if (bsp->repr != Seq_repr_seg) {
					/* ptr = StringMove (tbuf, embl_code); */
					ptr = tbuf;
					StringMove (ptr, tsip->accession);
					tsip->name = StringSave (tbuf);
					SeqMgrReplaceInBioseqIndex (bsp);
				}
				return;
			}
		}
	}
	
	return;
}

extern void DoFixupLocus (SeqEntryPtr sep);
extern void DoFixupLocus (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  CharPtr       embl_code;
  OrgRefPtr     orp;
  ValNodePtr    sdp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        DoFixupLocus (sep);
      }
      return;
    }
  }
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  /* embl_code = get_embl_code (orp); */
  embl_code = NULL;
  SeqEntryExplore (sep, (Pointer) embl_code, FixLocus);
  MemFree (embl_code);
}

static TextSeqIdPtr GetTSIPforBSP (BioseqPtr bsp)

{
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip = NULL;

  if (bsp == NULL) return NULL;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_TPG) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    }
  }

  return tsip;
}

static void DoFixPartLocus (SeqEntryPtr sep, CharPtr prefix, Int2 segment, Int2 digits)

{
  BioseqPtr     bsp;
  Char          ch;
  Char          digs [16];
  CharPtr       ptr;
  Char          str [32];
  TextSeqIdPtr  tsip;

  if (sep == NULL || prefix == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_const) return;
  tsip = GetTSIPforBSP (bsp);
  if (tsip == NULL) return;
  tsip->name = MemFree (tsip->name);
  sprintf (digs, "%*d", (int) digits, (int) segment);
  ptr = digs;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == ' ') {
      *ptr = '0';
    }
    ptr++;
    ch = *ptr;
  }
  sprintf (str, "%sS%s", prefix, digs);
  tsip->name = StringSave (str);
}

extern void DoFixupSegSet (SeqEntryPtr sep);
extern void DoFixupSegSet (SeqEntryPtr sep)

{
  CharPtr       accn;
  Uint1         chs;
  Int2          digits;
  BioseqPtr     nbsp;
  BioseqSetPtr  bssp;
  size_t        len;
  SeqEntryPtr   nsep;
  Int2          numsegs;
  BioseqSetPtr  pbssp;
  CharPtr       prefix;
  SeqEntryPtr   psep;
  BioseqPtr     sbsp;
  Int2          segment;
  SeqIdPtr      sip;
  Char          str [32];
  Char          tmp [128];
  TextSeqIdPtr  tsip;

  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      DoFixupSegSet (sep);
    }
    return;
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;
  if (! IS_Bioseq (nsep)) return;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  if (nbsp == NULL) return;
  if (nbsp->repr != Seq_repr_seg) return;
  sbsp = NULL;
  tsip = GetTSIPforBSP (nbsp);
  if (tsip == NULL || StringHasNoText (tsip->accession)) {
    sbsp = nbsp;
    nsep = FindNthSeqEntry (sep, 4);
    if (nsep == NULL) return;
    if (! IS_Bioseq (nsep)) return;
    nbsp = (BioseqPtr) nsep->data.ptrvalue;
    if (nbsp == NULL) return;
    if (nbsp->repr != Seq_repr_raw && nbsp->repr != Seq_repr_const) return;
    tsip = GetTSIPforBSP (nbsp);
  }
  if (tsip == NULL) return;
  chs = SEQID_GENBANK;
  for (sip = nbsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_TPG) {
      chs = sip->choice;
    }
  }
  accn = tsip->accession;
  if (StringHasNoText (accn)) return;
  psep = FindBioseqSetByClass (sep, BioseqseqSet_class_parts);
  if (psep == NULL) return;
  pbssp = (BioseqSetPtr) psep->data.ptrvalue;
  if (pbssp == NULL) return;
  for (sep = pbssp->seq_set, numsegs = 0; sep != NULL; sep = sep->next, numsegs++) continue;
  StringNCpy_0 (tmp, accn, sizeof (tmp));
  prefix = tmp;
  len = StringLen (prefix);
  if (numsegs > 999) {
    digits = 4;
  } else if (numsegs > 99) {
    digits = 3;
  } else if (numsegs > 9) {
    digits = 2;
  } else {
    digits = 1;
  }
  if (digits + 1 + len > 16) {
    prefix += digits + 1 + len - 16;
  }
  
  for (sep = pbssp->seq_set, segment = 1; sep != NULL; sep = sep->next, segment++) {
    DoFixPartLocus (sep, prefix, segment, digits);
  }
  
  if (sbsp == NULL) return;
  tsip = NULL;
  for (sip = sbsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_TPG) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    }
  }
  if (tsip == NULL) {
    tsip = TextSeqIdNew ();
    if (tsip == NULL) return;
    ValNodeAddPointer (&(sbsp->id), chs, (Pointer) tsip);
  }
  tsip->name = MemFree (tsip->name);
  sprintf (str, "SEG_%sS", prefix);
  tsip->name = StringSave (str);
  SeqMgrReplaceInBioseqIndex (sbsp);
}

static void ForceLocusFixup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ErrSev       sev;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sev = ErrSetMessageLevel (SEV_FATAL);
  WatchCursor ();
  Update ();
  DoFixupLocus (sep);
  DoFixupSegSet (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
}

extern void GetRidOfRedundantSourceNotes (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);

static void CleanupGenbankBlockCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GBBlockPtr     gbp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    empty = FALSE;
    if (sdp->choice == Seq_descr_genbank && sdp->data.ptrvalue == NULL) {
      empty = TRUE;
    } else if (sdp->choice == Seq_descr_genbank && sdp->data.ptrvalue != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      gbp->source = MemFree (gbp->source); /* remove within Sequin */
      gbp->origin = MemFree (gbp->origin);
      gbp->taxonomy = MemFree (gbp->taxonomy);
      if (gbp->extra_accessions == NULL && gbp->source == NULL &&
          gbp->keywords == NULL && gbp->origin == NULL &&
          gbp->date == NULL && gbp->entry_date == NULL &&
          gbp->div == NULL && gbp->taxonomy == NULL) {
        empty = TRUE;
      }
    }
    if (empty) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static void CheckBioSourceForEnvironmentalSample (BioSourcePtr biop, Pointer userdata)
{
  BoolPtr bp;
  
  if (biop == NULL || userdata == NULL) return;
  
  bp = (BoolPtr) userdata;
  if (! *bp) return;
  
  if (biop->org == NULL || biop->org->orgname == NULL)
  {
  	*bp = FALSE;
  	return;
  }
  if (StringSearch (biop->org->orgname->lineage, "environmental sample") == NULL)
  {
  	*bp = FALSE;
  }
}

static Boolean CheckForEnvironmentalSample (BioseqSetPtr bssp)
{
  Boolean all_enviro = TRUE;
  
  if (bssp == NULL) return FALSE;
  
  VisitBioSourcesInSet (bssp, (Pointer) &all_enviro, CheckBioSourceForEnvironmentalSample);
   
  return all_enviro; 
}


static void ConvertToEcoSets (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  sep_entry;
  
  if (sep == NULL) return;
  
  if (!IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  
  if (bssp->_class == BioseqseqSet_class_mut_set
      || bssp->_class == BioseqseqSet_class_pop_set 
      || bssp->_class == BioseqseqSet_class_phy_set)
  {
  	if (CheckForEnvironmentalSample (bssp))
  	{
  	  bssp->_class = BioseqseqSet_class_eco_set;
  	}
  }
  for (sep_entry = bssp->seq_set; sep_entry != NULL; sep_entry = sep_entry->next)
  {
  	if (IS_Bioseq_set (sep_entry))
  	{
  	  ConvertToEcoSets (sep_entry);
  	}
  }
}

static void AddQualsFromLineageAndDiv (BioSourcePtr biop, Pointer userdata)
{
  ValNode      vn;
  
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return;

  if (StringISearch (biop->org->orgname->lineage, "environmental sample") != NULL
      || StringCmp (biop->org->orgname->div, "ENV") == 0)
  {
    vn.choice = SourceQualChoice_textqual;
    vn.data.intvalue = Source_qual_environmental_sample;
    vn.next = NULL;
    SetSourceQualInBioSource (biop, &vn, NULL, "", ExistingTextOption_replace_old);
  }
  if (StringISearch (biop->org->orgname->lineage, "metagenomes") != NULL) {
    vn.choice = SourceQualChoice_textqual;
    vn.data.intvalue = Source_qual_metagenomic;
    vn.next = NULL;
    SetSourceQualInBioSource (biop, &vn, NULL, "", ExistingTextOption_replace_old);
  }
}

extern void ForceCleanupEntityID (Uint2 entityID)

{
  SeqEntryPtr   sep;
  ModTextFixPtr tfp;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  SeqMgrClearFeatureIndexes (entityID, NULL);
  SeqEntryExplore (sep, NULL, CleanupGenbankBlockCallback);
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  MySeqEntryToAsn3 (sep, TRUE, FALSE, TRUE);
  CorrectGenCodes (sep, entityID);
  CleanUpPseudoProducts (entityID, sep);
  SeqEntryExplore (sep, NULL, GetRidOfRedundantSourceNotes);
  RenormalizeNucProtSets (sep, TRUE);
  CdCheck (sep, NULL);
  ConvertToEcoSets (sep);
  VisitBioSourcesInSep (sep, NULL, AddQualsFromLineageAndDiv);
  tfp = ModTextFixNew ();
  VisitBioSourcesInSep (sep, tfp, RemoveTextFromTextFreeSubSourceModifiers);
  MemFree (tfp);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
}

static void ForceCleanupBtn (IteM i, ButtoN b, Boolean validate)

{
  BaseFormPtr   bfp;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  ForceCleanupEntityID (bfp->input_entityID);
  if (validate) {
    ValSeqEntryForm (bfp->form);
  }
}


static void ForceCleanup (IteM i)

{
  ForceCleanupBtn (i, NULL, TRUE);
}

static void MapWholeToInt (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr  bsp;
  SeqLocPtr  location;
  SeqIntPtr  sint;
  SeqIdPtr   sip;

  if (sfp == NULL) return;
  location = sfp->location;
  if (location == NULL) return;
  if (location->choice != SEQLOC_WHOLE) return;
  sip = SeqLocId (location);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sint = SeqIntNew ();
  if (sint == NULL) return;
  sint->from = 0;
  sint->to = bsp->length - 1;
  sint->id = SeqIdDup (sip);
  location = ValNodeAddPointer (NULL, SEQLOC_INT, (Pointer) sint);
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = location;
}

static void ConvertWholeToInterval (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, MapWholeToInt);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
}

static void CorrectBothStrands (SeqFeatPtr sfp, Pointer userdata)

{
  SeqLocPtr  location, slp = NULL;
  SeqIntPtr  sint;
  SeqPntPtr  spp;

  if (sfp == NULL) return;
  location = sfp->location;
  if (location == NULL) return;

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    switch (slp->choice) {
      case SEQLOC_INT :
        sint = (SeqIntPtr) slp->data.ptrvalue;
        if (sint != NULL) {
          if (sint->strand == Seq_strand_both) {
            sint->strand = Seq_strand_plus;
          } else if (sint->strand == Seq_strand_both_rev) {
            sint->strand = Seq_strand_minus;
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          if (spp->strand == Seq_strand_both) {
            spp->strand = Seq_strand_plus;
          } else if (spp->strand == Seq_strand_both_rev) {
            spp->strand = Seq_strand_minus;
          }
        }
        break;
      default :
        break;
    }
  }
}

static void ConvertBothStrands (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, CorrectBothStrands);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
}

static void AssignFeatIDs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  AssignFeatureIDs (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void ReassignFeatIDs (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ans = Message (MSG_OKC, "Are you sure you want to reassign feature identifiers?");
  if (ans == ANS_CANCEL) return;

  ReassignFeatureIDs (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void UniqifyFeatIDs (IteM i)

{
  MsgAnswer     ans;
  BaseFormPtr   bfp;
  BioseqSetPtr  bssp;
  Int4          count;
  Int4          offset = 0;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_genbank) return;

  ans = Message (MSG_OKC, "Are you sure you want to make merged feature identifiers unique?");
  if (ans == ANS_CANCEL) return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    count = FindHighestFeatureID (sep);
    OffsetFeatureIDs (sep, offset);
    OffsetFeatureIDXrefs (sep, offset);
    offset += count;
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void ClearFeatIDsAndLinks (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ans = Message (MSG_YN, "Are you sure you want to remove feature IDs and links?");
  if (ans == ANS_NO) return;

  ClearFeatureIDs (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void LinkByOverlap (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  LinkCDSmRNAbyOverlap (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void LinkByProduct (IteM i)

{
  BaseFormPtr  bfp;
  ValNodePtr   bsplist;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  bsplist = LockFarComponentsEx (sep, FALSE, TRUE, TRUE, NULL);
  LinkCDSmRNAbyProduct (sep);
  bsplist = UnlockFarComponents (bsplist);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void LinkByLabel (IteM i)

{
  BaseFormPtr  bfp;
  ValNodePtr   bsplist;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  bsplist = LockFarComponentsEx (sep, FALSE, TRUE, TRUE, NULL);
  LinkCDSmRNAbyLabel (sep);
  bsplist = UnlockFarComponents (bsplist);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}


extern void RemoveFeatureLink (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  SeqFeatXrefPtr  xref, next, PNTR prevlink;
  ObjectIdPtr     oip;
  SeqFeatPtr      link_sfp;
  Char            buf [32];
  CharPtr         str = NULL;

  if (sfp1 == NULL) return;

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp1->xref);
  xref = sfp1->xref;
  while (xref != NULL) {
    next = xref->next;
    link_sfp = NULL;

    if (xref->id.choice == 3) {
      oip = (ObjectIdPtr) xref->id.value.ptrvalue;
      if (oip != NULL) {
        if (StringDoesHaveText (oip->str)) {
          str = oip->str;
        } else {
          sprintf (buf, "%ld", (long) oip->id);
          str = buf;
        }
        link_sfp = SeqMgrGetFeatureByFeatID (sfp1->idx.entityID, NULL, str, NULL, NULL);
      }
    }
    if (link_sfp == sfp2) {
      *prevlink = xref->next;
      xref->next = NULL;
      MemFree (xref);
    } else {
      prevlink = (SeqFeatXrefPtr PNTR) &(xref->next);
    }

    xref = next;
  }
}


extern void LinkTwoFeatures (SeqFeatPtr dst, SeqFeatPtr sfp)

{
  ChoicePtr       cp;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref, prev_xref, next_xref;
  SeqFeatPtr      old_match;

  if (dst == NULL || sfp == NULL) return;

  cp = &(dst->id);
  if (cp == NULL) return;
  if (cp->choice == 3) {
    /* don't create a duplicate xref, remove links to other features */
    xref = sfp->xref;
    prev_xref = NULL;
    while (xref != NULL) {
      next_xref = xref->next;
      if (xref->id.choice == 3 && xref->id.value.ptrvalue != NULL) {
        if (ObjectIdMatch (cp->value.ptrvalue, xref->id.value.ptrvalue)) {
          /* already have this xref */
          return;
        } else {
          old_match = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
          RemoveFeatureLink (sfp, old_match);
          RemoveFeatureLink (old_match, sfp);
        }
      } else {
        prev_xref = xref;
      }
      xref = next_xref;
    }

    oip = (ObjectIdPtr) cp->value.ptrvalue;
    if (oip != NULL) {
      oip = AsnIoMemCopy (oip, (AsnReadFunc) ObjectIdAsnRead,
                          (AsnWriteFunc) ObjectIdAsnWrite);
      if (oip != NULL) {
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->id.choice = 3;
          xref->id.value.ptrvalue = (Pointer) oip;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
    }
  }
}

static void LinkSelected (IteM i)

{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp1, sfp2;
  SelStructPtr  ssp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  AssignFeatureIDs (sep);

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
    Message (MSG_OK, "No features selected");
    return;
  }
  if (ssp->itemtype != OBJ_SEQFEAT) {
    Message (MSG_OK, "Something other than a feature is selected");
    return;
  }
  sfp1 = SeqMgrGetDesiredFeature (ssp->entityID, NULL, ssp->itemID, 0, NULL, NULL);
  if (sfp1 == NULL) {
    Message (MSG_OK, "Unable to find selected feature");
    return;
  }
  ssp = ssp->next;
  if (ssp == NULL) {
    Message (MSG_OK, "Only one feature was selected, two must be selected");
    return;
  }
  if (ssp->itemtype != OBJ_SEQFEAT) {
    Message (MSG_OK, "Something other than a feature is selected");
    return;
  }
  sfp2 = SeqMgrGetDesiredFeature (ssp->entityID, NULL, ssp->itemID, 0, NULL, NULL);
  if (ssp->next != NULL) {
    Message (MSG_OK, "Only two features must be selected");
    return;
  }
  if (sfp2 == NULL) {
    Message (MSG_OK, "Unable to find selected feature");
    return;
  }

  LinkTwoFeatures (sfp1, sfp2);
  LinkTwoFeatures (sfp2, sfp1);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void SelCDSmRNALink (IteM i)

{
  Char            buf [32];
  ObjectIdPtr     oip;
  SeqFeatPtr      sfp;
  SelStructPtr    ssp;
  CharPtr         str = NULL;
  SeqFeatXrefPtr  xref;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
    Message (MSG_OK, "No feature selected");
    return;
  }
  if (ssp->next != NULL) {
    Message (MSG_OK, "Only one feature must be selected");
    return;
  }
  if (ssp->itemtype != OBJ_SEQFEAT) {
    Message (MSG_OK, "Something other than a feature is selected");
    return;
  }

  sfp = SeqMgrGetDesiredFeature (ssp->entityID, NULL, ssp->itemID, 0, NULL, NULL);
  if (sfp == NULL) {
    Message (MSG_OK, "Unable to find selected feature");
    return;
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        str = oip->str;
      } else {
        sprintf (buf, "%ld", (long) oip->id);
        str = buf;
      }
    }
  }

  if (str == NULL) {
    Message (MSG_OK, "Unable to extract feature ID xref");
    return;
  }

  sfp = SeqMgrGetFeatureByFeatID (ssp->entityID, NULL, str, NULL, NULL);
  if (sfp == NULL) {
    Message (MSG_OK, "Unable to find referenced feature");
    return;
  }

  ObjMgrAlsoSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
  ObjMgrSendMsg (OM_MSG_UPDATE, sfp->idx.entityID, 0, 0);
  Update ();
}

extern void ForceTaxonFixupBtn (IteM i, ButtoN b)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  /* SeqEntryExplore (sep, (Pointer) bfp, CleanupGenbankBlockCallback); */
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  MySeqEntryToAsn3 (sep, TRUE, FALSE, TRUE);
  CorrectGenCodes (sep, bfp->input_entityID);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ForceTaxonFixup (IteM i)

{
  ForceTaxonFixupBtn (i, NULL);
}

static void JustTaxonFixup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct hasdupfeatsdata {
  Uint1 featdef;
  Int4  num_feat;
  Int4  num_cds_with_product;
  ValNodePtr delete_list;
  Boolean case_sensitive;
  Boolean ignore_partials;
} HasDupFeatsData, PNTR HasDupFeatsPtr;

typedef struct dupfeatsdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG feat_type;
  ButtoN ignore_partials;
  Uint2  entityID;
  WindoW window;
} DupFeatsDlgData, PNTR DupFeatsDlgPtr;

static Pointer DupFeatsDlgToData (DialoG d)
{
  DupFeatsDlgPtr dlg;
  HasDupFeatsPtr data;
  ValNodePtr     vnp;

  dlg = (DupFeatsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  data = (HasDupFeatsPtr) MemNew (sizeof (HasDupFeatsData));
  data->num_feat = 0;
  data->num_cds_with_product = 0;
  data->delete_list = NULL;
  data->case_sensitive = TRUE;

  vnp = DialogToPointer (dlg->feat_type);
  if (vnp == NULL) {
    data->featdef = 0;
  } else {
    if (vnp->choice == Feature_type_any) {
      data->featdef = 0;
    } else {
      data->featdef = GetFeatdefFromFeatureType (vnp->choice);
    }
    vnp = ValNodeFreeData (vnp);
  }

  data->ignore_partials = GetStatus (dlg->ignore_partials);

  return data;
}


static DialoG DupFeatsDialog (GrouP g, Uint2 entityID, WindoW w)
{
  DupFeatsDlgPtr dlg;
  GrouP          p;
  ValNode        vn;
  ValNodePtr feature_list = NULL;

  dlg = (DupFeatsDlgPtr) MemNew (sizeof (DupFeatsDlgData));

  p = HiddenGroup (g, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = DupFeatsDlgToData;

  dlg->entityID = entityID;
  dlg->window = w;

  feature_list = ValNodeNew (NULL);
  feature_list->choice = Feature_type_any;
  feature_list->data.ptrvalue = StringSave ("All");
  AddAllFeaturesToChoiceList (&feature_list);
  dlg->feat_type = ValNodeSelectionDialog (p, feature_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                NULL, NULL, FALSE);

  vn.choice = Feature_type_any;
  vn.data.ptrvalue = NULL;
  vn.next = NULL;
  PointerToDialog (dlg->feat_type, &vn);
  dlg->ignore_partials = CheckBox (p, "Ignore partials", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feat_type, (HANDLE) dlg->ignore_partials, NULL);

  return (DialoG) p;
}


static void HasDupFeatsBioseqCallback (BioseqPtr bsp, Pointer data)
{
  HasDupFeatsPtr h;
  SeqFeatPtr     sfp1, sfp2;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || data == NULL) {
    return;
  }

  h = (HasDupFeatsPtr) data;

  sfp1 = SeqMgrGetNextFeature (bsp, NULL, 0,h->featdef, &fcontext);
  while (sfp1 != NULL) {
    sfp2 = SeqMgrGetNextFeature (bsp, sfp1, 0, h->featdef, &fcontext);
    if (DoFeaturesMatch (sfp1, sfp2, FALSE, h->case_sensitive, h->ignore_partials)) {
      h->num_feat++;
      if (sfp2->data.choice == SEQFEAT_CDREGION && sfp2->product != NULL) {
        h->num_cds_with_product++;
      }
      ValNodeAddPointer (&(h->delete_list), OBJ_SEQFEAT, sfp2);
    }
    sfp1 = sfp2;
  }
}


static Boolean RemoveDuplicateFeatsAction (Uint2 entityID, HasDupFeatsPtr data)
{
  SeqEntryPtr sep;
  Boolean     rval = FALSE;
  MsgAnswer       ans;
  Boolean         remove_proteins = FALSE;
  ValNodePtr      vnp;
  SeqFeatPtr      sfp;
  BioseqPtr       protbsp;

  sep = GetTopSeqEntryForEntityID (entityID);

  VisitBioseqsInSep (sep, data, HasDupFeatsBioseqCallback);
  if (data->num_feat == 0) {
    Message (MSG_ERROR, "No duplicates found");
    return rval;
  }

  if (data->num_cds_with_product > 0) {
    ans = Message (MSG_YNC, "%d duplicate features found.  %d coding regions have protein products.  Remove protein products?",
                   data->num_feat, data->num_cds_with_product);
    if (ans == ANS_CANCEL) {
      data->delete_list = ValNodeFree (data->delete_list);
      return rval;
    } else if (ans == ANS_YES) {
      remove_proteins = TRUE;
    }
  } 

  for (vnp = data->delete_list; vnp != NULL; vnp = vnp->next) {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL) {
      if (remove_proteins && sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) {
        protbsp = BioseqFindFromSeqLoc (sfp->product);
        if (protbsp != NULL) {
          protbsp->idx.deleteme = TRUE;
        }
      }
      sfp->idx.deleteme = TRUE;
    }
  }

  DeleteMarkedObjects (entityID, 0, NULL);
  if (remove_proteins) {
    RenormalizeNucProtSets (sep, TRUE);   	
  }

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  rval = TRUE;
  return rval;
}


static void DoRemoveDuplicateFeats (ButtoN b)
{
  DialoG         d;
  DupFeatsDlgPtr dlg;
  HasDupFeatsPtr data;

  d = (DialoG) GetObjectExtra (b);
  dlg = (DupFeatsDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) {
    return;
  }

  data = (HasDupFeatsPtr) DialogToPointer (dlg->dialog);

  if (RemoveDuplicateFeatsAction (dlg->entityID, data)) {
    Remove (dlg->window);
  }
  data = MemFree (data);
}


static void RemoveDuplicateFeatsWithOptions (IteM i)
{
  BaseFormPtr  bfp;
  WindoW       w;
  GrouP        h, c;
  ButtoN       b;
  DialoG       d;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Remove Duplicate Features", StdCloseWindowProc);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  d = DupFeatsDialog (h, bfp->input_entityID, w);
  c = HiddenGroup (h, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveDuplicateFeats);
  SetObjectExtra (b, d, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) d, (HANDLE) c, NULL);

  Show (w);
}


static void NewRemoveDuplicateFeats (IteM i, Uint2 featdef)
{
  BaseFormPtr  bfp;
  HasDupFeatsData h;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  h.featdef = (Uint1) featdef;
  h.num_cds_with_product = 0;
  h.num_feat = 0;
  h.delete_list = NULL;
  h.case_sensitive = TRUE;
  h.ignore_partials = FALSE;

  RemoveDuplicateFeatsAction (bfp->input_entityID, &h);
}


static void NewRemoveDuplicateFeatsAll (IteM i)
{
  NewRemoveDuplicateFeats (i, 0);
}


static void NewRemoveDuplicateFeatsGenes (IteM i)
{
  NewRemoveDuplicateFeats (i, FEATDEF_GENE);
}



/*
static Boolean LIBCALLBACK RemoveGeneXrefsOnBioseqs (BioseqPtr bsp, SeqMgrBioseqContextPtr bcontext)

{
  SeqFeatXrefPtr     curr;
  SeqMgrFeatContext  fcontext;
  GeneRefPtr         grp;
  GeneRefPtr         grpx;
  SeqFeatXrefPtr     PNTR last;
  SeqFeatXrefPtr     next;
  Boolean            redundantgenexref;
  SeqFeatPtr         sfp;
  SeqFeatPtr         sfpx;
  CharPtr            syn1;
  CharPtr            syn2;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice != SEQFEAT_GENE) {
      grp = SeqMgrGetGeneXref (sfp);
      if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
        sfpx = SeqMgrGetOverlappingGene (sfp->location, NULL);
        if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
          grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          if (grpx != NULL) {
            redundantgenexref = FALSE;
            if ((! StringHasNoText (grp->locus)) && (! StringHasNoText (grpx->locus))) {
              if ((StringICmp (grp->locus, grpx->locus) == 0)) {
                redundantgenexref = TRUE;
              }
            } else if ((! StringHasNoText (grp->locus_tag)) && (! StringHasNoText (grpx->locus_tag))) {
              if ((StringICmp (grp->locus_tag, grpx->locus_tag) == 0)) {
                redundantgenexref = TRUE;
              }
            } else if (grp->syn != NULL && grpx->syn != NULL) {
              syn1 = (CharPtr) grp->syn->data.ptrvalue;
              syn2 = (CharPtr) grpx->syn->data.ptrvalue;
              if ((! StringHasNoText (syn1)) && (! StringHasNoText (syn2))) {
                if ((StringICmp (syn1, syn2) == 0)) {
                  redundantgenexref = TRUE;
                }
              }
            }
            if (redundantgenexref) {
              last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
              curr = sfp->xref;
              while (curr != NULL) {
                next = curr->next;
                if (curr->data.choice == SEQFEAT_GENE) {
                  *last = next;
                  curr->next = NULL;
                  SeqFeatXrefFree (curr);
                } else {
                  last = &(curr->next);
                }
                curr = next;
              }
            }
          }
        }
      } 
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  return TRUE;
}
*/

typedef struct dummysmfedata {
  Int4  max;
  Int4  num_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK SQNDummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  Int4          len;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
  }

  return TRUE;
}

static Boolean LIBCALLBACK RemoveGeneXrefsOnBioseqs (BioseqPtr bsp, SeqMgrBioseqContextPtr bcontext)

{
  Int2               count;
  SeqFeatXrefPtr     curr;
  DummySmfeData      dsd;
  SeqMgrFeatContext  fcontext;
  GeneRefPtr         grp;
  GeneRefPtr         grpx;
  SeqFeatXrefPtr     PNTR last;
  SeqFeatXrefPtr     next;
  Boolean            redundantgenexref;
  SeqFeatPtr         sfp;
  SeqFeatPtr         sfpx;
  CharPtr            syn1, syn2;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice != SEQFEAT_GENE) {
      grp = SeqMgrGetGeneXref (sfp);
      if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
        sfpx = SeqMgrGetOverlappingGene (sfp->location, NULL);
        if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
          grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          if (grpx != NULL) {

            redundantgenexref = TRUE;
            if ((!StringHasNoText (grp->locus)) && (!StringHasNoText (grpx->locus))) {
              if ((StringICmp (grp->locus, grpx->locus) != 0)) {
                redundantgenexref = FALSE;
              }
            } else if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grp->locus_tag)) {
              if ((StringICmp (grp->locus_tag, grpx->locus_tag) != 0)) {
                redundantgenexref = FALSE;
              }
            } else if (grp->syn != NULL && grpx->syn != NULL) {
              syn1 = (CharPtr) grp->syn->data.ptrvalue;
              syn2 = (CharPtr) grpx->syn->data.ptrvalue;
              if ((!StringHasNoText (syn1)) && (!StringHasNoText (syn2))) {
                if ((StringICmp (syn1, syn2) != 0)) {
                  redundantgenexref = FALSE;
                }
              }
            }

            if (redundantgenexref) {
              MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
              dsd.max = INT4_MAX;
              dsd.num_at_max = 0;
              count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                                       LOCATION_SUBSET, (Pointer) &dsd, SQNDummySMFEProc);

              if (dsd.num_at_max < 2) {
                last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
                curr = sfp->xref;
                while (curr != NULL) {
                  next = curr->next;
                  if (curr->data.choice == SEQFEAT_GENE) {
                    *last = next;
                    curr->next = NULL;
                    SeqFeatXrefFree (curr);
                  } else {
                    last = &(curr->next);
                  }
                  curr = next;
                }
              }
            }
          }
        }
      } 
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  return TRUE;
}

static void RemoveGeneXrefs (IteM i)

{
  BaseFormPtr  bfp;
  Uint2        entityID;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  entityID = SeqMgrIndexFeatures (bfp->input_entityID, NULL);
  if (entityID < 1) return;
  SeqMgrExploreBioseqs (entityID, 0, NULL, RemoveGeneXrefsOnBioseqs, TRUE, TRUE, TRUE);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveGeneXNoLC (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr       bsp;
  GeneRefPtr      grp;
  SeqFeatXrefPtr  curr;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  if (StringHasNoText (grp->locus_tag)) return;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  if (SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, NULL) != NULL) return;
  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  curr = sfp->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      *last = next;
      curr->next = NULL;
      SeqFeatXrefFree (curr);
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
}

static void RemoveGeneXrefsNoGene (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, RemoveGeneXNoLC);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveGeneXNoLT (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr       bsp;
  GeneRefPtr      grp;
  SeqFeatXrefPtr  curr;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  if (StringHasNoText (grp->locus_tag)) return;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  if (SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, NULL) != NULL) return;
  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  curr = sfp->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      *last = next;
      curr->next = NULL;
      SeqFeatXrefFree (curr);
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
}

static void RemoveGeneXrefsNoLocTag (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, RemoveGeneXNoLT);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoRemoveNSGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatXrefPtr  curr;
  GeneRefPtr      grp;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;

  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  curr = sfp->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) curr->data.value.ptrvalue;
      if (SeqMgrGeneIsSuppressed (grp)) {
        last = &(curr->next);
      } else {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      }
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
}

static void RemoveNSGeneXrefs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, DoRemoveNSGeneXrefs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoRemoveGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatXrefPtr  curr;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;

  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  curr = sfp->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      *last = next;
      curr->next = NULL;
      SeqFeatXrefFree (curr);
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
}

static void RemoveAllGeneXrefs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, DoRemoveGeneXrefs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct removegenexref
{
  FORM_MESSAGE_BLOCK
  DialoG  feature_select;
  DialoG  constraints;
  DialoG  accept_cancel;
  
  ButtoN  suppressing;
  ButtoN  nonsuppressing;
  ButtoN  unnecessary;
    
  Boolean do_suppressing;
  Boolean do_nonsuppressing;
  Boolean do_unnecessary;
} RemoveGeneXrefData, PNTR RemoveGeneXrefPtr;

static void RemoveGeneXrefsChangeNotify (Pointer userdata)
{
  RemoveGeneXrefPtr dlg;
  ValNodePtr        vnp;
  Boolean           do_enable = FALSE;
  
  dlg = (RemoveGeneXrefPtr) userdata;
  if (dlg == NULL) 
  {
    return;
  }
  
  vnp = (ValNodePtr) DialogToPointer (dlg->feature_select);
  if (vnp != NULL)
  {
    ValNodeFree (vnp);
    if (GetStatus (dlg->suppressing)
        || GetStatus (dlg->nonsuppressing)
        || GetStatus (dlg->unnecessary))
    {
      do_enable = TRUE;
    }
  }
  if (do_enable)
  {
    EnableAcceptCancelDialogAccept (dlg->accept_cancel);
  }
  else
  {
    DisableAcceptCancelDialogAccept (dlg->accept_cancel);
  }
} 

static void RemoveGeneXrefsChangeBtn (ButtoN b)
{
  RemoveGeneXrefPtr dlg;

  dlg = (RemoveGeneXrefPtr) GetObjectExtra (b);
  RemoveGeneXrefsChangeNotify (dlg);
}

static void RemoveGeneXrefsClear (Pointer data)
{
  RemoveGeneXrefPtr dlg;

  dlg = (RemoveGeneXrefPtr) data;
  if (dlg == NULL) return;
 
  PointerToDialog (dlg->feature_select, NULL);
  PointerToDialog (dlg->constraints, NULL);
}

static void RemoveGeneXrefsClearText (Pointer data)
{
  RemoveGeneXrefPtr    dlg;
  FilterSetPtr          fsp;

  dlg = (RemoveGeneXrefPtr) data;
  if (dlg == NULL) return;
 
  fsp = DialogToPointer (dlg->constraints);
  FilterSetClearText (fsp);
  PointerToDialog (dlg->constraints, fsp);
  FilterSetFree (fsp);
}

static Boolean IsUnnecessaryGeneXref (SeqFeatPtr sfp, GeneRefPtr grp)
{
  SeqFeatPtr sfpx;
  GeneRefPtr grpx;
  Boolean    redundantgenexref = FALSE;
  CharPtr    syn1;
  CharPtr    syn2;
  
  if (sfp == NULL || grp == NULL || SeqMgrGeneIsSuppressed (grp))
  {
    return FALSE;
  }
  
  sfpx = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE)
  {
    return FALSE;
  }

  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL)
  {
    return FALSE;
  }
  
  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
      redundantgenexref = TRUE;
    }
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if (StringICmp (grp->locus, grpx->locus) == 0) {
      redundantgenexref = TRUE;
    }
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) == 0) {
        redundantgenexref = TRUE;
      }
    }
  }
  return redundantgenexref;
}

static void RemoveGeneXrefsCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RemoveGeneXrefPtr    dlg;
  SeqFeatXrefPtr       curr;
  GeneRefPtr           grp;
  SeqFeatXrefPtr       PNTR last;
  SeqFeatXrefPtr       next;
  Boolean              is_suppressed;

  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  dlg = (RemoveGeneXrefPtr) userdata;

  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  curr = sfp->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) curr->data.value.ptrvalue;
      is_suppressed = SeqMgrGeneIsSuppressed (grp);
      
      if ((dlg->do_suppressing && is_suppressed)
          || (dlg->do_nonsuppressing && !is_suppressed)
          || (dlg->do_unnecessary && IsUnnecessaryGeneXref (sfp, grp)))
      {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      }
      else 
      {
        last = &(curr->next);
      }
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
 
}

static Boolean RemoveGeneXrefsAction (Pointer userdata)
{
  RemoveGeneXrefPtr    dlg;
  FilterSetPtr         fsp;
  SeqEntryPtr          sep;
  ValNodePtr           feature_type_list, vnp;
  Uint1                 feat_def_choice;
  
  if (userdata == NULL) return FALSE;
  
  dlg = (RemoveGeneXrefPtr) userdata;
  
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  if (sep == NULL) return FALSE;
  
  feature_type_list = (ValNodePtr) DialogToPointer (dlg->feature_select);
  
  if (feature_type_list == NULL)
  {
    return FALSE;
  }
  
  fsp = (FilterSetPtr) DialogToPointer (dlg->constraints);
  
  dlg->do_suppressing = GetStatus (dlg->suppressing);
  dlg->do_nonsuppressing = GetStatus (dlg->nonsuppressing);
  dlg->do_unnecessary = GetStatus (dlg->unnecessary);
    
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    feat_def_choice = vnp->choice;
    if (feat_def_choice == 255)
    {
      feat_def_choice = 0;
    }
    OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                         RemoveGeneXrefsCallback,
                                         NULL, 0, feat_def_choice, 0, dlg);
  }
  
  ValNodeFree (feature_type_list);
  FilterSetFree (fsp);
  
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);  
  Update ();
  return TRUE;
}

static void NewRemoveGeneXrefs (IteM i)

{
  BaseFormPtr       bfp;
  RemoveGeneXrefPtr dlg;
  WindoW            w;
  GrouP             h, k;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (RemoveGeneXrefPtr) MemNew (sizeof (RemoveGeneXrefData));
  if (dlg == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove Gene Xrefs", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupExtraProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  dlg->feature_select =  FeatureSelectionDialogEx (h, TRUE, sep,
                                                 RemoveGeneXrefsChangeNotify, 
                                                 dlg);
  
  k = NormalGroup (h, 0, 3, "Remove Gene Xrefs that are", programFont, NULL);
  dlg->suppressing = CheckBox (k, "Suppressing", RemoveGeneXrefsChangeBtn);
  SetObjectExtra (dlg->suppressing, dlg, NULL);
  dlg->nonsuppressing = CheckBox (k, "Non-Suppressing", RemoveGeneXrefsChangeBtn);
  SetObjectExtra (dlg->nonsuppressing, dlg, NULL);
  dlg->unnecessary = CheckBox (k, "Unnecessary", RemoveGeneXrefsChangeBtn);
  SetObjectExtra (dlg->unnecessary, dlg, NULL);
  
  dlg->constraints = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text");
  dlg->accept_cancel = AcceptCancelDialog (h, RemoveGeneXrefsAction, NULL, RemoveGeneXrefsClear, RemoveGeneXrefsClearText, (Pointer)dlg, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_select,
                              (HANDLE) k,
                              (HANDLE) dlg->constraints,
                              (HANDLE) dlg->accept_cancel, NULL);
                                
  Show (w);
}

static void DoRefreshGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatXrefPtr    curr;
  GeneRefPtr        grp, grpfeat;
  SeqFeatPtr        gene;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp;

  if (sfp == NULL) return;

  for (curr = sfp->xref; curr != NULL; curr = curr->next)
  {
    if (curr->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) curr->data.value.ptrvalue;
      if (grp != NULL)
      {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        gene = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
        if (gene != NULL && gene->data.choice == SEQFEAT_GENE) {
          grpfeat = (GeneRefPtr) gene->data.value.ptrvalue;
          if (grpfeat != NULL) {
            GeneRefFree (grp);
            grp = GeneRefDup (grpfeat);
            curr->data.value.ptrvalue = grp;
          }
        }
      }
    }
  }
}

static void RefreshGeneXRefs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, DoRefreshGeneXrefs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static ValNodePtr RemoveDbxrefList (ValNodePtr vnp)

{
  ValNodePtr  next;

  while (vnp != NULL) {
    next = vnp->next;
    DbtagFree ((DbtagPtr) vnp->data.ptrvalue);
    MemFree (vnp);
    vnp = next;
  }
  return NULL;
}

static void DoRemoveDbxrefs (SeqFeatPtr sfp, Pointer userdata)

{
  GeneRefPtr  grp;
  Uint1       feattype;

  if (sfp == NULL) return;

  if (userdata != NULL && (feattype = *(Uint1 PNTR)userdata) != sfp->data.choice) return;
  sfp->dbxref = RemoveDbxrefList (sfp->dbxref);
  if (sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
  grp->db = RemoveDbxrefList (grp->db);
}

static void DoRemoveBioSourceDbxrefs (BioSourcePtr biop, Pointer userdata)
{
  if (biop == NULL || biop->org == NULL) return;
  
  biop->org->db = RemoveDbxrefList (biop->org->db);
}

static void RemoveDbxrefsCommon (IteM i, Uint1 PNTR subtype, Boolean remove_source, Boolean remove_features)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (remove_features) {
    VisitFeaturesInSep (sep, subtype, DoRemoveDbxrefs);
  }
  if (subtype == NULL && remove_source)
  {
    VisitBioSourcesInSep (sep, NULL, DoRemoveBioSourceDbxrefs);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveGeneDbxrefs (IteM i)
{
  Uint1        subtype = SEQFEAT_GENE;
  RemoveDbxrefsCommon (i, &subtype, FALSE, TRUE);
}

static void RemoveRNADbxrefs (IteM i)

{
  Uint1        subtype = SEQFEAT_RNA;

  RemoveDbxrefsCommon (i, &subtype, FALSE, TRUE);
}

static void RemoveCDSDbxrefs (IteM i)

{
  Uint1        subtype = FEATDEF_CDS;
  RemoveDbxrefsCommon (i, &subtype, FALSE, TRUE);
}

static void RemoveAllFeatureDbxrefs (IteM i)
{
  RemoveDbxrefsCommon (i, NULL, FALSE, TRUE);
}

static void RemoveAllDbxrefs (IteM i)

{
  RemoveDbxrefsCommon (i, NULL, TRUE, TRUE);
}


static void RemoveAllBioSourceDbxrefs (IteM i)

{
  RemoveDbxrefsCommon (i, NULL, TRUE, FALSE);
}


static void ClearPubFig (PubdescPtr pdp, Pointer userdata)

{
  if (pdp == NULL) return;
  pdp->fig = MemFree (pdp->fig);
  pdp->num = NumberingFree (pdp->num);
}

static void RemovePubFig (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitPubdescsInSep (sep, NULL, ClearPubFig);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveSelFeats (IteM i)

{
  MsgAnswer          ans;
  BaseFormPtr        bfp;
  BioseqPtr          cdna;
  SeqMgrFeatContext  fcontext;
  BioseqPtr          prot;
  Boolean            remove_asked;
  Boolean            remove_mrnas = FALSE;
  Boolean            remove_prots = FALSE;
  SelStructPtr       sel;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  SelStructPtr       ssp;
  Boolean            unremoved_feats = FALSE;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) return;

  remove_asked = FALSE;
  for (sel = ssp; sel != NULL; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) 
    {
      unremoved_feats = TRUE;
      continue;      
    }
    if (sel->itemtype == OBJ_SEQFEAT) {
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        sfp->idx.deleteme = TRUE;
        if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) {
          prot = BioseqFindFromSeqLoc (sfp->product);
          if (prot != NULL) {
            if (! remove_asked) {
              ans = Message (MSG_YN, "Remove protein products?");
              if (ans == ANS_YES) {
                remove_prots = TRUE;
              }
              remove_asked = TRUE;
            }
            if (remove_prots) {
              prot->idx.deleteme = TRUE;
            }
          }
        }
      }
    }
  }

  remove_asked = FALSE;
  for (sel = ssp; sel != NULL; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) 
    {
      unremoved_feats = TRUE;
      continue;      
    }
    if (sel->itemtype == OBJ_SEQFEAT) {
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        sfp->idx.deleteme = TRUE;
        if (sfp->data.choice == SEQFEAT_RNA && sfp->product != NULL) {
          cdna = BioseqFindFromSeqLoc (sfp->product);
          if (cdna != NULL) {
            if (! remove_asked) {
              ans = Message (MSG_YN, "Remove mRNA products?");
              if (ans == ANS_YES) {
                remove_mrnas = TRUE;
              }
              remove_asked = TRUE;
            }
            if (remove_mrnas) {
              cdna->idx.deleteme = TRUE;
            }
          }
        }
      }
    }
  }

  ObjMgrSelect (0, 0, 0, 0, NULL);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  if (unremoved_feats)
  {
    Message (MSG_ERROR, "Warning!  Features mapped to far sequences cannot be deleted!");
  }
}

static void RemoveUnselFeats (IteM i)

{
  MsgAnswer          ans;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  BioseqPtr          cdna;
  SeqMgrFeatContext  fcontext;
  BioseqPtr          prot;
  Boolean            remove_asked;
  Boolean            remove_mrnas = FALSE;
  Boolean            remove_prots = FALSE;
  SelStructPtr       sel;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  SelStructPtr       ssp;
  SeqEntryPtr        scope;
  Boolean            sel_on_local = FALSE;
  Boolean            sel_on_far = FALSE;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) return;

  scope = SeqEntrySetScope (NULL);
  bsp = NULL;
  for (sel = ssp; sel != NULL && bsp == NULL; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) 
    {
      if (sel->itemtype == OBJ_SEQFEAT)
      {
        sel_on_far = TRUE;
      }
      continue;
    }    
    if (sel->itemtype == OBJ_SEQFEAT) {
      sel_on_local = TRUE;
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
      }
    }
  }
  SeqEntrySetScope (scope);
  if (bsp == NULL) 
  {
    if (sel_on_far && ! sel_on_local)
    {
      Message (MSG_ERROR, "Warning!  Features mapped to far sequences cannot be deleted!");
    }
    return;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    sfp->idx.deleteme = TRUE;
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  for (sel = ssp; sel != NULL; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) continue;
    if (sel->itemtype == OBJ_SEQFEAT) {
      sfp = SeqMgrGetDesiredFeature (0, bsp, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        sfp->idx.deleteme = FALSE;
      }
    }
  }

  remove_asked = FALSE;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.deleteme) {
      if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) {
        prot = BioseqFindFromSeqLoc (sfp->product);
        if (prot != NULL) {
          if (! remove_asked) {
            ans = Message (MSG_YN, "Remove protein products?");
            if (ans == ANS_YES) {
              remove_prots = TRUE;
            }
            remove_asked = TRUE;
          }
          if (remove_prots) {
            prot->idx.deleteme = TRUE;
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  remove_asked = FALSE;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.deleteme) {
      if (sfp->data.choice == SEQFEAT_RNA && sfp->product != NULL) {
        cdna = BioseqFindFromSeqLoc (sfp->product);
        if (cdna != NULL) {
          if (! remove_asked) {
            ans = Message (MSG_YN, "Remove mRNA products?");
            if (ans == ANS_YES) {
              remove_mrnas = TRUE;
            }
            remove_asked = TRUE;
          }
          if (remove_mrnas) {
            cdna->idx.deleteme = TRUE;
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  ObjMgrSelect (0, 0, 0, 0, NULL);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveFeatsOutOfSelRange (IteM i)

{
  MsgAnswer          ans;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  BioseqPtr          cdna;
  SeqMgrFeatContext  fcontext;
  BioseqPtr          prot;
  Boolean            remove_asked;
  Boolean            remove_mrnas = FALSE;
  Boolean            remove_prots = FALSE;
  SelStructPtr       sel;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  SeqLocPtr          slp = NULL;
  SelStructPtr       ssp;
  SeqEntryPtr        scope;
  Boolean            sel_on_local = FALSE;
  Boolean            sel_on_far = FALSE;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
    Message (MSG_ERROR, "Warning!  No features are selected - please select one");
    return;
  }
  if (ssp->next != NULL) {
    Message (MSG_ERROR, "Warning!  Multiple features are selected - please select only one");
    return;
  }

  scope = SeqEntrySetScope (NULL);
  bsp = NULL;
  for (sel = ssp; sel != NULL && bsp == NULL; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) 
    {
      if (sel->itemtype == OBJ_SEQFEAT)
      {
        sel_on_far = TRUE;
      }
      continue;
    }    
    if (sel->itemtype == OBJ_SEQFEAT) {
      sel_on_local = TRUE;
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        slp = sfp->location;
        bsp = BioseqFindFromSeqLoc (slp);
      }
    }
  }
  SeqEntrySetScope (scope);
  if (bsp == NULL) 
  {
    if (sel_on_far && ! sel_on_local)
    {
      Message (MSG_ERROR, "Warning!  Please select only a feature on a near sequence");
    }
    return;
  }
  if (slp == NULL) {
    Message (MSG_ERROR, "Warning!  Failure to get location from selected feature");
    return;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (SeqLocAinB (sfp->location, slp) >= 0) {
      sfp->idx.deleteme = FALSE;
    } else {
      sfp->idx.deleteme = TRUE;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  remove_asked = FALSE;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.deleteme) {
      if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) {
        prot = BioseqFindFromSeqLoc (sfp->product);
        if (prot != NULL) {
          if (! remove_asked) {
            ans = Message (MSG_YN, "Remove protein products?");
            if (ans == ANS_YES) {
              remove_prots = TRUE;
            }
            remove_asked = TRUE;
          }
          if (remove_prots) {
            prot->idx.deleteme = TRUE;
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  remove_asked = FALSE;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.deleteme) {
      if (sfp->data.choice == SEQFEAT_RNA && sfp->product != NULL) {
        cdna = BioseqFindFromSeqLoc (sfp->product);
        if (cdna != NULL) {
          if (! remove_asked) {
            ans = Message (MSG_YN, "Remove mRNA products?");
            if (ans == ANS_YES) {
              remove_mrnas = TRUE;
            }
            remove_asked = TRUE;
          }
          if (remove_mrnas) {
            cdna->idx.deleteme = TRUE;
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  ObjMgrSelect (0, 0, 0, 0, NULL);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveCDDRegions (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  FreeCDDRegions (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveCDDAligns (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  FreeCDDAligns (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveCDDDups (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  RemoveDuplicateCDDs (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RetainBestCDD (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  LeaveBestCDD (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void LastDitchLookup (BioSourcePtr biop)

{
  DbtagPtr      dbt;
  Char          orgname [256];
  OrgRefPtr     orp;
  CharPtr       ptr;
  CharPtr       binomial_end = NULL, trinomial_end = NULL;
  SubSourcePtr  ssp;
  Int4          taxID;
  ValNodePtr    vnp;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  /* if this organism has already been assigned a taxon ID, skip it */
  if (orp->db != NULL) {
    for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
      dbt = (DbtagPtr) vnp->data.ptrvalue;
      if (dbt != NULL) {
        if (StringICmp (dbt->db, "taxon") == 0) return;
      }
    }
  }
  /* we are looking for a name that might become a recognizable name if we
   * truncated it after the first four words (if the third word is subsp.),
   * after the first three words (if the third word is not subsp.),
   * or after the first two words if the above two attempts failed.
   */
 
  StringNCpy_0 (orgname, orp->taxname, sizeof (orgname));
  if (StringHasNoText (orgname)) return;
  ptr = StringChr (orgname, ' ');
  if (ptr == NULL) return;
  /* skip over the first word and the spaces after it. */
  while (*ptr == ' ') 
  {
    ptr++;
  }
  ptr = StringChr (ptr, ' ');
  /* if there are only two words, give up. */
  if (ptr == NULL) 
  {
    return;
  }
  binomial_end = ptr;
  while (*ptr == ' ')
  {
    ptr++;
  }
  if (StringNCmp (ptr, "subsp.", 6) == 0)
  {
    ptr += 6;
    while (*ptr == ' ' )
    {
      ptr++;
    }
  }
  
  ptr = StringChr (ptr, ' ');
  if (ptr != NULL)
  {
    trinomial_end = ptr;
  }
  
  /* see if trinomial produces a tax server hit */
  taxID = 0;
  if (trinomial_end != NULL)
  {
    *trinomial_end = '\0';
    taxID = Taxon3GetTaxIdByName (orgname);
  }
  if (taxID < 1)
  {
    /* see if binomial produces a tax server hit */
    *binomial_end = '\0';
    taxID = Taxon3GetTaxIdByName (orgname);
  }
  
  if (taxID < 1) return;
  
  /* create a note that contains the original name of the organism, and truncate
   * the actual organism after the binomial or trinomial.
   */
  ssp = SubSourceNew ();
  if (ssp != NULL) {
    ssp->subtype = 255;
    ssp->name = orp->taxname;
	orp->taxname = NULL;
    ssp->next = biop->subtype;
    biop->subtype = ssp;
    SetTaxNameAndRemoveTaxRef (orp, StringSave (orgname));
  }
}

static void LastDitchTaxonFixup (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)

{
	BioseqPtr bsp = NULL;
	BioseqSetPtr bssp = NULL;
	ValNodePtr descr, vnp;
	BioSourcePtr biosp = NULL;
	SeqAnnotPtr annot, sap;
	SeqFeatPtr sfp;
	
	if (IS_Bioseq(sep)) {
		bsp = sep->data.ptrvalue;
		descr = bsp->descr;
		annot = bsp->annot;
	} else {
		bssp = sep->data.ptrvalue;
		descr = bssp->descr;
		annot = bssp->annot;
	}
	for (vnp = descr; vnp; vnp = vnp->next) {
		if (vnp->choice == Seq_descr_source) {
			biosp = vnp->data.ptrvalue;
			LastDitchLookup (biosp); 
		}
	}
	for (sap = annot; sap != NULL; sap = sap->next) {
		if (sap->type == 1) {
			for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
				if (sfp->data.choice == SEQFEAT_BIOSRC) {
					biosp = sfp->data.value.ptrvalue;
					LastDitchLookup (biosp);
				}
			}
		}
	}
}

static void GenSpecTaxonFixup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, LastDitchTaxonFixup);
  ForceTaxonFixupBtn (i, NULL);
}


static void MergeToPartsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  BoolPtr ordered;
  
  if (sfp == NULL) return;
  ordered = (BoolPtr) userdata;
  
  MergeFeatureIntervalsToParts (sfp, *ordered);
}

static void MergeToParts (IteM i, Boolean ordered)
{
  BaseFormPtr       bfp;
  SelStructPtr      ssp;
  Boolean           some_sel_not_feat = FALSE;
  Boolean           some_sel_feat = FALSE;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  SeqEntryPtr       sep;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
    if (bfp == NULL) {
      Message (MSG_ERROR, "MergeToParts error");
      return;
    }
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep == NULL) return;
    VisitFeaturesInSep (sep, (Pointer) &ordered, MergeToPartsCallback);
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  }
  else
  {
    while (ssp != NULL)
    {
      if (ssp->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (ssp->entityID, NULL, ssp->itemID, 0, NULL, &fcontext);
        if (sfp != NULL)
        {
          MergeFeatureIntervalsToParts (sfp, ordered);
          some_sel_feat = TRUE;
          ObjMgrSetDirtyFlag (ssp->entityID, TRUE);
          ObjMgrSendMsg (OM_MSG_UPDATE, ssp->entityID, 0, 0);
        }
      }
      else
      {
        some_sel_not_feat = TRUE;        
      }
      ssp = ssp->next;
    }
    if (!some_sel_feat)
    {
      Message (MSG_ERROR, "No features selected!");
    }
    else if (some_sel_not_feat)
    {
      Message (MSG_ERROR, "Some selected items were not features!");
    }
  }
  Update ();
}

extern void MergeToPartsJoin (IteM i)
{
  MergeToParts (i, FALSE);
}

extern void MergeToPartsOrdered (IteM i)
{
  MergeToParts (i, TRUE);
}


static Boolean MergeSegSeqCallback (GatherContextPtr gcp)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Boolean       noLeft;
  Boolean       noRight;
  RnaRefPtr     rrp;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
  SeqLocPtr     slp;
  tRNAPtr       trna;

  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->location == NULL) return TRUE;
  bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  if (bsp == NULL) return TRUE;
  if (ISA_aa (bsp->mol)) return TRUE;
  if (bsp->repr != Seq_repr_seg) {
    sep = GetBestTopParentForData (gcp->entityID, bsp);
    if (sep == NULL) return TRUE;
    sep = FindNucSeqEntry (sep);
    if (sep == NULL || sep->choice != 1) return TRUE;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return TRUE;
  }
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
  if (slp == NULL) return TRUE;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  FreeAllFuzz (sfp->location);
  SetSeqLocPartial (sfp->location, noLeft, noRight);
  sfp->partial = (sfp->partial || noLeft || noRight);
  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          slp = SeqLocMerge (bsp, cbp->loc, NULL, FALSE, TRUE, FALSE);
          if (slp != NULL) {
            cbp->loc = SeqLocFree (cbp->loc);
            cbp->loc = slp;
            FreeAllFuzz (cbp->loc);
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) {
        trna = rrp->ext.value.ptrvalue;
        if (trna != NULL && trna->anticodon != NULL) {
          slp = SeqLocMerge (bsp, trna->anticodon, NULL, FALSE, TRUE, FALSE);
          if (slp != NULL) {
            trna->anticodon = SeqLocFree (trna->anticodon);
            trna->anticodon = slp;
            FreeAllFuzz (trna->anticodon);
          }
        }
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static void MergeToSegSeq (IteM i)

{
  BaseFormPtr   bfp;
  GatherScope   gs;
  SelStructPtr  ssp;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
    if (bfp == NULL) {
      Message (MSG_ERROR, "MergeToSegSeq error");
      return;
    }
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    GatherEntity (bfp->input_entityID, NULL, MergeSegSeqCallback, &gs);
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    /*
    Message (MSG_ERROR, "No item selected");
    */
    return;
  }
  if (ssp->itemtype != OBJ_SEQFEAT) {
    Message (MSG_ERROR, "Feature must be selected");
    return;
  }
  GatherItem (ssp->entityID, ssp->itemID, ssp->itemtype, NULL, MergeSegSeqCallback);
  ObjMgrSetDirtyFlag (ssp->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ssp->entityID, 0, 0);
}

static void SegToRawCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || (bsp->repr != Seq_repr_seg && bsp->repr != Seq_repr_delta)) return;
  SegOrDeltaBioseqToRaw (bsp);
}

static void SegSeqToRawSeq (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, (Pointer) bfp, SegToRawCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Message (MSG_OK, "Some manual desktop manipulations remain");
}

static void CollectBioseqsForConversion (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR list;
  
  if (bsp == NULL || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)) return;
  if (userdata == NULL)
  {
    return;
  }
  list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (list, 0, bsp);
}

static ValNodePtr GetAnnotListForBspList (ValNodePtr bsp_list)
{
  BioseqPtr   bsp;
  SeqAnnotPtr sanp;
  ValNodePtr  complete_annot_list = NULL, align_annot_list;
  ValNodePtr  annot_vnp, vnp, bsp_vnp;
  Boolean     found;

  for (bsp_vnp = bsp_list; bsp_vnp != NULL; bsp_vnp = bsp_vnp->next)
  {
    bsp = (BioseqPtr)(bsp_vnp->data.ptrvalue);
    align_annot_list = FindAlignSeqAnnotsForBioseq (bsp);
    for (annot_vnp = align_annot_list; annot_vnp != NULL; annot_vnp = annot_vnp->next)
    {
      sanp = (SeqAnnotPtr) annot_vnp->data.ptrvalue;
      if (sanp->type == 2 && sanp->data != NULL)
      {
        for (vnp = complete_annot_list, found = FALSE;
             vnp != NULL && !found;
             vnp = vnp->next)
        {
          if (vnp->data.ptrvalue == sanp)
          {
            found = TRUE;
          }
        }
        if (!found)
        {
          ValNodeAddPointer (&complete_annot_list, 0, sanp);
        }
      }
    }
    align_annot_list = ValNodeFree (align_annot_list);
  }
  return complete_annot_list;
}


static void 
CleanupAlignmentsAfterConversion 
(ValNodePtr adjusted_bsp_list,
 Uint2 entityID)
{
  ValNodePtr  align_annot_list, annot_vnp;
  SeqAnnotPtr sanp;
  
  align_annot_list = GetAnnotListForBspList (adjusted_bsp_list);
  for (annot_vnp = align_annot_list;
       annot_vnp != NULL; 
       annot_vnp = annot_vnp->next)
  {
    sanp = annot_vnp->data.ptrvalue;
    if (sanp != NULL && sanp->type == 2)
    {
      ConsolidateSegmentsOverKnownLengthGaps (sanp->data);
      FixOneAlignmentOverGaps (sanp->data, entityID);
    }
  }
  align_annot_list = ValNodeFree (align_annot_list);
  DeleteMarkedObjects (entityID, 0, NULL);
  
}

static void CorrectAlignmentsAfterNToGapConversion (SeqEntryPtr sep, Uint2 entityID)
{
  ValNodePtr adjusted_bsp_list = NULL;
  
  VisitBioseqsInSep (sep, &adjusted_bsp_list, CollectBioseqsForConversion);
  CleanupAlignmentsAfterConversion (adjusted_bsp_list, entityID);
  adjusted_bsp_list = ValNodeFree (adjusted_bsp_list);
}

static void ConvertNsToGapsWithSizeList (Uint2 entityID, Int4Ptr gap_sizes, Boolean adjust_cds)
{
  ValNodePtr  adjusted_bsp_list = NULL;
  SeqEntryPtr sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, &adjusted_bsp_list, CollectBioseqsForConversion);

  if (adjust_cds)
  {
    /* remove the gap locations from the coding regions first */
    VisitBioseqsInSep (sep, gap_sizes, PrepareCodingRegionLocationsForDeltaConversionCallback);
  }
  
  VisitBioseqsInSep (sep, gap_sizes, ConvertNsToGaps);
  CleanupAlignmentsAfterConversion (adjusted_bsp_list, entityID);
  adjusted_bsp_list = ValNodeFree (adjusted_bsp_list);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}


static void RawSeqToDeltaSeq (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  ConvertNsToGapsWithSizeList (bfp->input_entityID, NULL, FALSE);
}

static void RawSeqToDeltaSeqUnknownLengthGaps (IteM i)

{
  BaseFormPtr  bfp;
  Int4         gap_sizes[2];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  gap_sizes [0] = 100;
  gap_sizes [1] = -1;
  ConvertNsToGapsWithSizeList (bfp->input_entityID, gap_sizes, FALSE);
}

static void RawSeqToDeltaSeqUnknown100LengthGaps (IteM i)

{
  BaseFormPtr  bfp;
  Int4         gap_sizes[2];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  gap_sizes [0] = -1;
  gap_sizes [1] = 0;
  ConvertNsToGapsWithSizeList (bfp->input_entityID, gap_sizes, FALSE);
}

#if 0
static void ConvertSegSetToDeltaItem (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  Int4         gap_sizes[2];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ConvertSegSetsToDeltaSequences (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}
#endif


typedef struct gapconversiondata 
{
  FEATURE_FORM_BLOCK

  GrouP  unknown_op;
  TexT   unknown_val_txt;
  GrouP  known_op;
  TexT   known_val_txt;
  ButtoN acceptBtn;
  DoC    explanation;
  ButtoN adjust_CDS_locations;
  
} GapConversionData, PNTR GapConversionPtr;

static ParData faParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData faColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static void SetConvertGapsAcceptAndText (GapConversionPtr gcp)
{
  Char             str[15];
  Int4 unknown_val, known_val;
  Int4 unknown_op, known_op;
  Char explanation[300];
  	
  if (gcp == NULL || gcp->explanation == NULL) return;
  GetTitle (gcp->unknown_val_txt, str, sizeof (str));
  unknown_val = atoi (str);
  unknown_op = GetValue (gcp->unknown_op);
  GetTitle (gcp->known_val_txt, str, sizeof (str));
  known_val = atoi (str);
  known_op = GetValue (gcp->known_op);

  Reset (gcp->explanation);
  explanation [0] = 0;
  
  if (unknown_val < 0 || known_val < 0)
  {
    sprintf (explanation, "%s", "Negative values are not permitted.");
  	Disable (gcp->acceptBtn);
  }
  else if (known_val == unknown_val)
  {
    sprintf (explanation, 
  	         "You must specify different values for known and unknown gap conversion sizes.");
  	Disable (gcp->acceptBtn);
  }  
  else if (unknown_val == 0 && unknown_op == 1 && known_val == 0 && known_op == 1)
  {
    sprintf (explanation, "%s", "This combination will have no effect.");
    Disable (gcp->acceptBtn);
  }
  else if (unknown_val == 0 && unknown_op == 1 && (known_val == 0 || known_val == 1) && known_op == 2)
  {
    sprintf (explanation, "%s", "All sequences of Ns will be converted to gaps of known length.");
    Enable (gcp->acceptBtn);
  }
  else if ((unknown_val == 0 || unknown_val == 1) && unknown_op == 2 && known_val == 0 && known_op == 1)
  {
    sprintf (explanation, "%s", "All sequences of Ns will be converted to gaps of unknown length (size 100).");
    Enable (gcp->acceptBtn);
  }
  else if (unknown_val == 0 && unknown_op == 1)
  {  	
  	if (known_op == 1)
  	{
  	  sprintf (explanation, 
  	           "All sequences of exactly %d Ns will be converted to gaps of known length.\n"
  	           "All other sequences of Ns will remain as Ns.", known_val);
   	  Enable (gcp->acceptBtn);
 	}
 	else
 	{
	  sprintf (explanation, 
  	             "All sequences of %d or more Ns will be converted to gaps of known length.\n"
  	             "All sequences of less than %d Ns will remain as Ns.", known_val, known_val);
  	  Enable (gcp->acceptBtn);
 	}
  }
  else if (known_val == 0 && known_op == 1)
  {
  	if (unknown_op == 1)
  	{
  	  sprintf (explanation, 
  	           "All sequences of exactly %d Ns will be converted to gaps of unknown length (size 100).\n"
  	           "All other sequences of Ns will remain as Ns.", unknown_val);
  	  Enable (gcp->acceptBtn);
  	}
  	else if (unknown_op == 2)
  	{
  	  sprintf (explanation, 
  	           "All sequences of %d or more Ns will be converted to gaps of unknown length.\n"
  	           "All sequences of less than %d Ns will remain as Ns.", unknown_val, unknown_val);
  	  Enable (gcp->acceptBtn);
  	}
  }
  else if (known_op == 1 && unknown_op == 1)
  {
  	sprintf (explanation, 
  	         "All sequences of exactly %d Ns will be converted to gaps of known length.\n"
  	         "All sequences of exactly %d Ns will be converted to gaps of unknown length (size 100).\n"
  	         "All sequences of Ns with lengths other than %d or %d will remain as Ns.", 
  	         known_val, unknown_val, known_val, unknown_val);
  	Enable (gcp->acceptBtn);
  }
  else if (known_op == 2 && unknown_op == 2)
  {
    /* handle cases where both operators are greater than/equal to */
  	if (known_val > unknown_val)
  	{
  	  if (unknown_val == 0 || unknown_val == 1)
  	  {
  	    sprintf (explanation, 
  	             "All sequences of Ns with lengths >= 1 and <= %d will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths >= %d will be converted to gaps of known length.\n",
  	             known_val -1, known_val);         
  	  }
  	  else
  	  {
  	  	sprintf (explanation,
  	             "All sequences of Ns with lengths >= %d and <= %d will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths >= %d will be converted to gaps of known length.\n"
  	             "All sequences of Ns with lengths <= %d will remain as Ns.",
  	             unknown_val, known_val -1, known_val, unknown_val - 1);         
  	  }
  	}
  	else
  	{
  	  if (known_val == 0 || known_val == 1)
  	  {
  	    sprintf (explanation, 
  	             "All sequences of Ns with lengths >= 1 and <= %d will be converted to gaps of known length.\n"
  	             "All sequences of Ns with lengths >= %d will be converted to gaps of unknown length (size 100).\n",
  	             unknown_val -1, unknown_val);         
  	  }
  	  else
  	  {
  	  	sprintf (explanation,
  	             "All sequences of Ns with lengths >= %d and <= %d will be converted to gaps of known length.\n"
  	             "All sequences of Ns with lengths >= %d will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths <= %d will remain as Ns.",
  	             known_val, unknown_val -1, unknown_val, known_val - 1);         
  	  }
  	}
  }
  else if (known_val > unknown_val)
  {
  	if (known_op == 1)
  	{
  	  /* we know that unknown_op is 2 */
      if (unknown_val == 0 || unknown_val == 1)
  	  {
        sprintf (explanation, 
  	           "All sequences of exactly %d Ns will be converted to gaps of known length.\n"
  	           "All sequences of Ns with lengths >= 1 and <= %d or lengths >= %d will be converted to gaps of unknown length (size 100).\n",
  	           known_val, known_val - 1, known_val + 1);
  	  }
  	  else
  	  {
        sprintf (explanation, 
  	            "All sequences of exactly %d Ns will be converted to gaps of known length.\n"
  	             "All sequences of Ns with lengths >= %d and <= %d or lengths >= %d will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths <= %d will remain as Ns.", 
  	             known_val, unknown_val, known_val - 1, known_val + 1, unknown_val - 1);
  	  }
  	  Enable (gcp->acceptBtn);
  	}
  	else if (known_op == 2)
  	{
  	  /* we know that unknown_op is 1 */
  	  if (unknown_val == 1)
  	  {
  	  	sprintf (explanation,
  	           "All sequences of exactly 1 Ns will be converted to gaps of unknown length (size 100).\n"
  	           "All sequences with lengths >= %d will be converted to gaps of known length.\n"
  	           "All sequences of Ns with lengths >= 2 and <= %d will remain as Ns",
  	  	        known_val, known_val - 1);
  	  }
  	  else
  	  {
        sprintf (explanation, 
  	           "All sequences of exactly %d Ns will be converted to gaps of unknown length (size 100).\n"
  	           "All sequences of Ns with lengths >= %d will be converted to gaps of known length.\n"
  	           "All sequences of Ns with lengths >= %d or lengths >= %d and <= %d will remain as Ns.  ",
  	           unknown_val, known_val, unknown_val - 1, unknown_val + 1, known_val - 1);
  	  }
  	  Enable (gcp->acceptBtn);
  	}
  }
  else
  {
  	/* unknown_val > known_val */
    if (known_op == 1)
    {
      /* we know that unknown_op is 2 */
  	  sprintf (explanation, 
  	           "All sequences of exactly %d Ns will be converted to gaps of known length.\n"
  	           "All sequences of Ns with lengths >= %d will be converted to gaps of unknown length (size 100).\n"
  	           "All sequences of Ns with lengths <= %d or lengths >= %d and <= %d will remain as Ns.\n",
  	            known_val, unknown_val, known_val - 1, known_val + 1, unknown_val - 1);
    }
    else
    {
      /* we know that unknown_op is 1, known_op is 2 */
      if (known_val == 0 || known_val == 1)
      {
      	sprintf (explanation,
  	             "All sequences of exactly %d Ns will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths >= 1 and <= %d or lengths >= %d will be converted to gaps of known length.\n",
      	         unknown_val, unknown_val - 1, unknown_val + 1);
      }
      else 
      {
      	sprintf (explanation,
  	             "All sequences of exactly %d Ns will be converted to gaps of unknown length (size 100).\n"
  	             "All sequences of Ns with lengths >= %d and <= %d or lengths >= %d will be converted to gaps of known length.\n"
  	             "All sequences of Ns with lengths <= %d will remain as Ns.",
      	         unknown_val, known_val, unknown_val - 1, unknown_val + 1, known_val);
      	
      }
    }
    Enable (gcp->acceptBtn);
  }  
  AppendText (gcp->explanation, explanation, &faParFmt, &faColFmt, programFont);
  UpdateDocument (gcp->explanation, 0, 0);
  Update ();	
}

static void SetGapsConvertAcceptButtonAndTextGroup (GrouP g)
{
  GapConversionPtr gcp;
  
  gcp = (GapConversionPtr) GetObjectExtra (g);
  SetConvertGapsAcceptAndText (gcp);
}

static void SetGapsConvertAcceptButtonAndTextText (TexT nText)
{
  GapConversionPtr gcp;
  
  gcp = (GapConversionPtr) GetObjectExtra (nText);
  SetConvertGapsAcceptAndText (gcp);
}

static void ConvertGaps (ButtoN b)
{
  GapConversionPtr gcp;
  SeqEntryPtr      sep;
  Int4             gap_sizes[2];
  Char             str[15];
  Int4             unknown_val, known_val;
  Int4             unknown_op, known_op;
  
  if (b == NULL) return;
  gcp = (GapConversionPtr) GetObjectExtra (b);
  if (gcp == NULL) return;
  sep = GetTopSeqEntryForEntityID (gcp->input_entityID);
  if (sep == NULL) return;
  
  GetTitle (gcp->unknown_val_txt, str, sizeof (str));
  unknown_val = atoi (str);
  unknown_op = GetValue (gcp->unknown_op);
  if (unknown_op == 1)
  {
  	gap_sizes [0] = unknown_val;
  }
  else
  {
    if (unknown_val == 0)
      unknown_val = 1;
  	gap_sizes[0] = 0 - unknown_val;
  }
  GetTitle (gcp->known_val_txt, str, sizeof (str));
  known_val = atoi (str);
  known_op = GetValue (gcp->known_op);
  if (known_op == 1)
  {
  	gap_sizes [1] = known_val;
  }
  else
  {
    if (known_val == 0)
      known_val = 1;
  	gap_sizes[1] = 0 - known_val;
  }
  
  
  Hide (gcp->form);
  ConvertNsToGapsWithSizeList (gcp->input_entityID, gap_sizes, GetStatus (gcp->adjust_CDS_locations));

  Remove (gcp->form);  
}

static void ListRawSequencesInAlignments (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR raw_list;
  Char            id_str [128];
  
  if (bsp == NULL || bsp->repr != Seq_repr_raw || userdata == NULL)
  {
  	return;
  }
  
  if (!IsBioseqInAnyAlignment (bsp, bsp->idx.entityID))
  {
    return;
  }
  
  raw_list = (ValNodePtr PNTR) userdata;
  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
  ValNodeAddPointer (raw_list, 0, StringSave (id_str));
}

static Boolean ContinueWithSequencesInAlignments (SeqEntryPtr sep)
{
  ValNodePtr raw_in_aln = NULL;
  CharPtr    msg;
  MsgAnswer  ans;
  
  VisitBioseqsInSep (sep, &raw_in_aln, ListRawSequencesInAlignments);
  if (raw_in_aln != NULL)
  {
  	msg = CreateListMessage ("Sequence", 
  	                         raw_in_aln->next == NULL 
  	                              ? " is in an alignment.  The alignment will be automatically adjusted after conversion.  Do you want to continue?" 
  	                              : " are in alignments.  The alignment will be automatically adjusted after conversion.  Do you want to continue?",
  	                              raw_in_aln);
  	raw_in_aln = ValNodeFreeData (raw_in_aln);
    ans = Message (MSG_YN, msg);
    msg = MemFree (msg);
    if (ans == ANS_NO)
    {
      return FALSE;
    }
  }
  return TRUE;
}

static void ListDeltaSequences (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR delta_list;
  Char            id_str [128];
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta || userdata == NULL)
  {
  	return;
  }
  
  delta_list = (ValNodePtr PNTR) userdata;
  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
  ValNodeAddPointer (delta_list, 0, StringSave (id_str));
}

static Boolean ContinueWithDeltaSequences (SeqEntryPtr sep)
{
  ValNodePtr delta_list = NULL;
  CharPtr    msg;
  MsgAnswer  ans;
  
  VisitBioseqsInSep (sep, &delta_list, ListDeltaSequences);
  if (delta_list != NULL)
  {
  	msg = CreateListMessage ("Sequence", 
  	                         delta_list->next == NULL 
  	                              ? " is already a delta sequence and gaps will not be added.  Do you want to continue?" 
  	                              : " are already delta sequences and gaps will not be added.  Do you want to continue?",
  	                              delta_list);
  	delta_list = ValNodeFreeData (delta_list);
    ans = Message (MSG_YN, msg);
    msg = MemFree (msg);
    if (ans == ANS_NO)
    {
      return FALSE;
    }
  }  
  return TRUE;
}


static void RawSeqToDeltaSeqUnknownWithUnknownLengthGaps (IteM i)
{
  BaseFormPtr      bfp;
  GapConversionPtr gcp;
  WindoW           w;
  GrouP            h, l, g, c;
  RecT             r;
  PrompT           p1, p2;
  SeqEntryPtr      sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (!ContinueWithDeltaSequences (sep))
  {
    return;
  }
  if (!ContinueWithSequencesInAlignments (sep))
  {
    return;
  }

  gcp = (GapConversionPtr) MemNew (sizeof (GapConversionData));
  if (gcp == NULL) return;
  
  gcp->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Convert Ns to Gaps of Unknown Length", StdCloseWindowProc);
  SetObjectExtra (w, gcp, StdCleanupFormProc);
  gcp->form = (ForM) w;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  p1 = StaticPrompt (h, "Choose size of length of Ns to convert to gaps of unknown length", 0, dialogTextHeight, programFont, 'c');
  g = HiddenGroup (h, 2, 0, NULL);
  gcp->unknown_op = HiddenGroup (g, 0, 2, SetGapsConvertAcceptButtonAndTextGroup);
  SetObjectExtra (gcp->unknown_op, gcp, NULL);
  RadioButton (gcp->unknown_op, "=");
  RadioButton (gcp->unknown_op, ">=");
  SetValue (gcp->unknown_op, 1);
  gcp->unknown_val_txt = DialogText (g, "100", 14, SetGapsConvertAcceptButtonAndTextText);
  SetObjectExtra (gcp->unknown_val_txt, gcp, NULL);
  p2 = StaticPrompt (h, "Choose size of length of Ns to convert to gaps of known length", 0, dialogTextHeight, programFont, 'c');
  l = HiddenGroup (h, 2, 0, NULL);
  gcp->known_op = HiddenGroup (l, 0, 2, SetGapsConvertAcceptButtonAndTextGroup);
  SetObjectExtra (gcp->known_op, gcp, NULL);
  RadioButton (gcp->known_op, "=");
  RadioButton (gcp->known_op, ">=");
  SetValue (gcp->known_op, 1);
  gcp->known_val_txt = DialogText (l, "0", 14, SetGapsConvertAcceptButtonAndTextText);
  SetObjectExtra (gcp->known_val_txt, gcp, NULL);
  
  /* status text */
  gcp->explanation = DocumentPanel (h, stdCharWidth * 27, stdLineHeight * 8);
  ObjectRect (gcp->explanation, &r);
  InsetRect (&r, 4, 4);
  faColFmt.pixWidth = r.right - r.left;
  
  gcp->adjust_CDS_locations = CheckBox (h, "Adjust CDS locations for gaps", NULL);
  SetStatus (gcp->adjust_CDS_locations, TRUE);

  c = HiddenGroup (h, 4, 0, NULL);
  gcp->acceptBtn = PushButton (c, "Accept", ConvertGaps);
  SetObjectExtra (gcp->acceptBtn, gcp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p1,
                              (HANDLE) g, 
                              (HANDLE) p2, 
                              (HANDLE) l,
                              (HANDLE) gcp->explanation,
                              (HANDLE) gcp->adjust_CDS_locations, 
                              (HANDLE) c, 
                              NULL);
  
  SetConvertGapsAcceptAndText (gcp);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static Uint2 deltaconversionedittypes [] = 
{
  TAGLIST_TEXT, TAGLIST_POPUP, TAGLIST_TEXT, TAGLIST_POPUP
};

ENUM_ALIST(deltaconversionedit_alist)
  {"Unknown length",  0},
  {"Known Length",    1},
END_ENUM_ALIST

ENUM_ALIST(deltaconversionedit_type_alist)
  {"Insert",  0},
  {"Replace",    1},
END_ENUM_ALIST

EnumFieldAssocPtr deltaconversionedit_alists [] = {
  NULL, deltaconversionedit_alist, NULL, deltaconversionedit_type_alist
};

static Uint2 deltaconversionedit_widths [] = {
  5, 5, 5, 5
};

typedef struct deltaconversion 
{
  FEATURE_FORM_BLOCK

  DialoG gap_locations;
  GrouP  coord_grp;
  ButtoN adjust_CDS_locations;
  
} DeltaConversionData, PNTR DeltaConversionPtr;


static SeqAlignPtr 
AdjustOneAlignmentForOneBioseqBasedOnGapLocations 
(SeqAlignPtr salp_to_adjust,
 BioseqPtr   bsp,
 ValNodePtr  gaps,
 SeqAlignPtr gap_salp)
{
  Int4          aln_index;
  SeqIdPtr      sip;
  ValNodePtr    gap_vnp;
  GapLocInfoPtr glip;
  SeqLocPtr     slp;
  Int4          gap_start;
  Int4          cumulative_offset = 0;
  DenseSegPtr   dsp;
  
  if (salp_to_adjust == NULL || salp_to_adjust->segtype != SAS_DENSEG
      || bsp == NULL || gaps == NULL)
  {
    return salp_to_adjust;
  }

  dsp = (DenseSegPtr) salp_to_adjust->segs;
  if (dsp == NULL) {
     return salp_to_adjust;
  }
  sip = bsp->id;
  aln_index = 0;
  while (sip != NULL && aln_index == 0)
  {
    aln_index = SeqIdOrderInBioseqIdList (sip, dsp->ids);
    if (aln_index == 0)
    {
      sip = sip->next;
    }
  }
  if (aln_index == 0) {
     /* bioseq not in alignment */
     return salp_to_adjust;
  }
  
  for (gap_vnp = gaps; gap_vnp != NULL; gap_vnp = gap_vnp->next)
  {
    glip = (GapLocInfoPtr) gap_vnp->data.ptrvalue;
    if (glip == NULL || glip->replace)
    {
      continue;
    }
    gap_start = glip->start_pos;
    if (gap_salp != NULL)
    {
      gap_start = AlnMgr2MapSeqAlignToBioseq(gap_salp, gap_start, aln_index);
    }
    if (gap_start < 1)
    {
      continue;
    }
    gap_start += cumulative_offset;
    if (gap_start > bsp->length)
    {
      continue;
    }
    slp = SeqLocIntNew (gap_start, gap_start + glip->length - 1, Seq_strand_plus, sip);
    salp_to_adjust = SeqAlignInsertByLoc (slp, salp_to_adjust);
    cumulative_offset += glip->length;
  }
  
  return salp_to_adjust;
}


static void ShiftAlignmentForGapLocations (SeqAnnotPtr sanp, ValNodePtr gaps)
{
  SeqAlignPtr   salp;
  DenseSegPtr   dsp;
  Int4          alignment_position = 0, k;
  Int4          cumulative_offset = 0;
  BoolPtr       shift_this_row;
  Int4          seg_num = 0, shift_seg;
  ValNodePtr    gap_vnp;
  GapLocInfoPtr glip;
  
  if (sanp == NULL || sanp->data == NULL || sanp->type != 2 || gaps == NULL)
  {
    return;
  }
  
  salp = (SeqAlignPtr) sanp->data;
  if (salp->segtype != SAS_DENSEG || salp->segs == NULL)
  {
    return;
  }
  
  dsp = (DenseSegPtr) salp->segs;
  
  shift_this_row = (BoolPtr) MemNew (sizeof (Boolean) * dsp->dim);
  
  for (gap_vnp = gaps; gap_vnp != NULL; gap_vnp = gap_vnp->next)
  {
    glip = (GapLocInfoPtr) gap_vnp->data.ptrvalue;
    if (glip == NULL || glip->replace)
    {
      continue;
    }
    while (seg_num < dsp->numseg 
           && alignment_position + dsp->lens [seg_num] < glip->start_pos + cumulative_offset)
           
    {
      alignment_position += dsp->lens [seg_num];
      seg_num++;
    }

    if (seg_num < dsp->numseg)
    {
      dsp->lens [seg_num] += glip->length;
      for (k = 0; k < dsp->dim; k++)
      {
        if (dsp->starts [seg_num * dsp->dim + k] == -1)
        {
          shift_this_row[k] = FALSE;
        }
        else
        {
          shift_this_row[k] = TRUE;
        }
      }
      
      /* shift the rows for which a gap was inserted */
      for (k = 0; k < dsp->dim; k++)
      {
        if (!shift_this_row[k])
        {
          continue;
        }
        /* shift start for minus strand rows for seg_num only */
        if (dsp->strands != NULL 
            && dsp->strands [seg_num * dsp->dim + k] == Seq_strand_minus
            && dsp->starts [seg_num * dsp->dim + k] != -1)
        {
          dsp->starts [seg_num * dsp->dim + k] -= glip->length;
        }
        /* shift starts for rows after seg_num */
        for (shift_seg = seg_num + 1;
             shift_seg < dsp->numseg;
             shift_seg++)
        {
          if (dsp->starts [shift_seg * dsp->dim + k] != -1)
          {
            if (dsp->strands == NULL
                || dsp->strands [shift_seg * dsp->dim + k] != Seq_strand_minus)
            {
              dsp->starts [shift_seg * dsp->dim + k] += glip->length;
            }
            else
            {
              dsp->starts [shift_seg * dsp->dim + k] -= glip->length;
            }
          }
        }
      }
    }
    cumulative_offset += glip->length;
  }  
  shift_this_row = MemFree (shift_this_row);
}

static void 
AdjustAlignmentsBasedOnGapLocations 
(ValNodePtr  adjusted_bsp_list,
 ValNodePtr  gaps,
 SeqAlignPtr salp,
 Uint2       entityID)
{
  ValNodePtr  align_annot_list, annot_vnp, bsp_vnp;
  SeqAnnotPtr sanp, adjustment_sanp = NULL;
  BioseqPtr   bsp;
  
  for (bsp_vnp = adjusted_bsp_list; bsp_vnp != NULL; bsp_vnp = bsp_vnp->next)
  {
    bsp = (BioseqPtr)(bsp_vnp->data.ptrvalue);
    align_annot_list = FindAlignSeqAnnotsForBioseq (bsp);
    for (annot_vnp = align_annot_list;
         annot_vnp != NULL; 
         annot_vnp = annot_vnp->next)
    {
      sanp = annot_vnp->data.ptrvalue;
      if (sanp != NULL && sanp->type == 2)
      {
        if (sanp->data == salp)
        {
          adjustment_sanp = sanp;
        }
        else
        {
          sanp->data = AdjustOneAlignmentForOneBioseqBasedOnGapLocations 
                              (sanp->data, bsp, gaps, salp);
        }
      }
    }
    align_annot_list = ValNodeFree (align_annot_list);
  }

  if (adjustment_sanp != NULL)
  {
    ShiftAlignmentForGapLocations (adjustment_sanp, gaps);
  }
    
  CleanupAlignmentsAfterConversion (adjusted_bsp_list, entityID);
}


static void 
ConvertRawBioseqToDelta 
(BioseqPtr bsp,
 ValNodePtr gaps,
 SeqAlignPtr salp,
 Int4        aln_row) 

{
  CharPtr       bases;
  Int4          len;
  ValNodePtr    seq_ext;
  SeqLitPtr     slp;
  IntFuzzPtr    ifp;
  ValNodePtr    gap_vnp;
  GapLocInfoPtr glip;
  Int4          orig_seq_offset;
  Char          tmp_ch;
  Int4          gap_start;
  SeqEntryPtr   sep;

  if (bsp == NULL || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)
      || gaps == NULL) 
  {
    return;
  }
  if (salp != NULL && aln_row < 1)
  {
    return;
  }

  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) return;

  seq_ext = NULL;
  len = 0;
  orig_seq_offset = 0;

  for (gap_vnp = gaps; gap_vnp != NULL; gap_vnp = gap_vnp->next)
  {
    glip = (GapLocInfoPtr) gap_vnp->data.ptrvalue;
    if (glip == NULL)
    {
      continue;
    }
    gap_start = glip->start_pos;
    /* remap for alignment coordinates if desired */
    if (salp != NULL)
    {
      gap_start = AlnMgr2MapSeqAlignToBioseq(salp, gap_start, aln_row);
    }
    if (gap_start < 1 || gap_start > bsp->length)
    {
      continue;
    }
    
    /* add data since last gap */
    if (gap_start - orig_seq_offset > 0) {
      slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slp != NULL) 
      {
        slp->length = gap_start - orig_seq_offset;
        ValNodeAddPointer (&(seq_ext), (Int2) 2, (Pointer) slp);
        slp->seq_data = (SeqDataPtr) BSNew (slp->length);
        slp->seq_data_type = Seq_code_iupacna;
        tmp_ch = bases [gap_start];
        bases [gap_start] = 0;
        AddBasesToByteStore ((ByteStorePtr) slp->seq_data, bases + orig_seq_offset);
        bases [gap_start] = tmp_ch;
        len += slp->length;
        orig_seq_offset += slp->length;
      }
    }
    
    /* add gap */
    slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
    if (slp != NULL) 
    {
      ValNodeAddPointer ((ValNodePtr PNTR) &(seq_ext), (Int2) 2, (Pointer) slp);
      if (glip->is_known)
      {
        slp->length = glip->length;
      }
      else
      {
        ifp = IntFuzzNew ();
        ifp->choice = 4;
        slp->fuzz = ifp;
        slp->length = 100;
      }
      len += slp->length;
    }
    if (glip->replace)
    {
      orig_seq_offset += glip->length;
    }
  }
  
  /* add remaining data after last gap to end */
  if (bsp->length - orig_seq_offset > 0) {
    slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
    if (slp != NULL) 
    {
      slp->length = bsp->length - orig_seq_offset;
      ValNodeAddPointer (&(seq_ext), (Int2) 2, (Pointer) slp);
      slp->seq_data = (SeqDataPtr) BSNew (slp->length);
      slp->seq_data_type = Seq_code_iupacna;
      AddBasesToByteStore ((ByteStorePtr) slp->seq_data, bases + orig_seq_offset);
      len += slp->length;
    }
  }
  
  MemFree (bases);

  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data_type = 0;
  bsp->repr = Seq_repr_delta;
  bsp->seq_ext_type = 4;
  bsp->seq_ext = seq_ext;
  bsp->length = len;

  BioseqPack (bsp);
  
  /* now adjust features for insertion */
  orig_seq_offset = 0;
  for (gap_vnp = gaps; gap_vnp != NULL; gap_vnp = gap_vnp->next)
  {
    glip = (GapLocInfoPtr) gap_vnp->data.ptrvalue;
    if (glip == NULL)
    {
      continue;
    }
    gap_start = glip->start_pos;
    /* remap for alignment coordinates if desired */
    if (salp != NULL)
    {
      gap_start = AlnMgr2MapSeqAlignToBioseq(salp, gap_start, aln_row);
    }
    if (gap_start < 1 || gap_start > bsp->length)
    {
      continue;
    }
    if (!glip->replace) {
      AdjustFeaturesForInsertion (bsp, bsp->id, 
                                gap_start + orig_seq_offset,
                                glip->length,
                                FALSE);
    } else if (glip->length != 100 && !glip->is_known) {
      AdjustFeaturesForInsertion (bsp, bsp->id, gap_start + orig_seq_offset,
                                  100 - glip->length,
                                  FALSE);
    }
    orig_seq_offset += glip->length;
  }

  sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);

  VisitFeaturesInSep (sep, (Pointer) Sequin_GlobalAlign2Seq, AdjustCDSLocationsForUnknownGapsCallback);
  
}

typedef struct converttodelta 
{
  ValNodePtr  location_list;
  ValNodePtr  affected_bioseq_list;
  SeqAlignPtr salp;
} ConvertToDeltaData, PNTR ConvertToDeltaPtr;

static void ConvertBioseqToDeltaWithSequenceGapList (BioseqPtr bsp, Pointer userdata)
{
  ConvertToDeltaPtr ctdp;
  
  if (bsp == NULL || bsp->id == NULL 
      || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)
      || userdata == NULL)
  {
    return;
  }
  ctdp = (ConvertToDeltaPtr) userdata;
  if (ctdp->location_list == NULL)
  {
    return;
  }
  ValNodeAddPointer (&(ctdp->affected_bioseq_list), 0, bsp);
  ConvertRawBioseqToDelta (bsp, ctdp->location_list, NULL, -1);
}

static void ConvertBioseqToDeltaWithAlignmentGapList (BioseqPtr bsp, Pointer userdata)
{
  Int4              aln_row = -1;
  SeqIdPtr          sip, sip_next;
  SeqAlignPtr       salp;
  ConvertToDeltaPtr ctdp;
  ValNodePtr        annot_list;
  SeqAnnotPtr       sanp;
  
  if (bsp == NULL || bsp->id == NULL 
      || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)
      || userdata == NULL)
  {
    return;
  }
  ctdp = (ConvertToDeltaPtr) userdata;
  if (ctdp->location_list == NULL)
  {
    return;
  }
  
  annot_list = FindAlignSeqAnnotsForBioseq (bsp);
  if (annot_list == NULL || annot_list->data.ptrvalue == NULL)
  {
    return;
  }
  sanp = (SeqAnnotPtr) annot_list->data.ptrvalue;
  if (sanp->type != 2 || sanp->data == NULL)
  {
    return;
  }

  salp = sanp->data;
  ctdp->salp = salp;
  
  /* Make sure alignment is indexed.
   */
  AlnMgr2IndexSeqAlignEx(salp, FALSE);
  
  sip = bsp->id;
  
  while (sip != NULL && aln_row == -1)
  {
    sip_next = sip->next;
    sip->next = NULL;
    aln_row = AlnMgr2GetFirstNForSip (salp, sip);
    sip->next = sip_next;
    sip = sip_next;
  }

  ValNodeAddPointer (&(ctdp->affected_bioseq_list), 0, bsp);
  ConvertRawBioseqToDelta (bsp, ctdp->location_list, salp, aln_row);
}

static int LIBCALLBACK SortGapLocations (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  GapLocInfoPtr glip1, glip2;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      glip1 = (GapLocInfoPtr) vnp1->data.ptrvalue;
      glip2 = (GapLocInfoPtr) vnp2->data.ptrvalue;
      if (glip1 != NULL && glip2 != NULL) {
        if (glip1->start_pos < glip2->start_pos)
        {
          return -1;
        }
        else if (glip1->start_pos > glip2->start_pos)
        {
          return 1;
        }
        else if (glip1->length < glip2->length)
        {
          return -1;
        }
        else if (glip1->length > glip2->length)
        {
          return 1;
        }
        else
        {
          return 0;
        }
      }
    }
  }
  return 0;
}

static Pointer DeltaLocToData (DialoG d)
{
  TagListPtr    tlp;
  ValNodePtr    result_list = NULL, vnp;
  CharPtr       str;
  Int4          start_pos, len;
  Boolean       is_known, replace;
  GapLocInfoPtr glip;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL) return NULL;
  
  for (vnp = tlp->vnp;
       vnp != NULL;
       vnp = vnp->next)
  {
    str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    if (StringHasNoText (str))
    {
      str = MemFree (str);
    }
    else
    {
      start_pos = atoi (str) - 1;
      str = MemFree (str);
      if (start_pos > 0)
      {
        str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 3);
        if (StringICmp (str, "1") == 0)
        {
          replace = TRUE;
        }
        else
        {
          replace = FALSE;
        }
        str = MemFree (str);

        str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
        if (StringICmp (str, "1") == 0)
        {
          str = MemFree (str);
          is_known = TRUE;
          len = -1;
          str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 2);
          if (!StringHasNoText (str))
          {
            len = atoi (str);
          }
          str = MemFree (str);
          if (len < 1)
          {
            Message (MSG_ERROR, "Must supply a length greater than zero for gaps of known length!");
            result_list = ValNodeFreeData (result_list);
            return NULL;
          }
        }
        else
        {
          str = MemFree (str);
          is_known = FALSE;
          if (replace)
          {
            str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 2);
            if (StringHasNoText (str))
            {
              len = 100;
            }
            else
            {
              len = atoi (str);
              if (len < 1)
              {
                len = 100;
              }
            }
            str = MemFree (str);
          }
          else
          {
            len = 100;
          }
        }
        
        glip = (GapLocInfoPtr) MemNew (sizeof (GapLocInfoData));
        if (glip != NULL)
        {
          glip->start_pos = start_pos;
          glip->is_known = is_known;
          glip->length = len;
          glip->replace = replace;
          ValNodeAddPointer (&result_list, 0, glip);
        }
      }
    }
  }
  
  result_list = ValNodeSort (result_list, SortGapLocations);
  return result_list;
}

static void DoConvertRawToDeltaWithGapLocations (ButtoN b)
{
  DeltaConversionPtr dcp;
  SeqEntryPtr        sep;
  ConvertToDeltaData ctdd;

  dcp = (DeltaConversionPtr) GetObjectExtra (b);
  if (dcp == NULL)
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (dcp->input_entityID);
  if (sep == NULL)
  {
    return;
  }

  ctdd.location_list = DialogToPointer (dcp->gap_locations);
  if (ctdd.location_list == NULL)
  {
    Message (MSG_ERROR, "Must supply valid gap locations!");
    return;
  }
  ctdd.affected_bioseq_list = NULL;
  ctdd.salp = NULL;
  
  Hide (dcp->form);
  WatchCursor ();
  Update (); 
  if (dcp->coord_grp == NULL || GetValue (dcp->coord_grp) == 1)
  {
    VisitBioseqsInSep (sep, &ctdd, ConvertBioseqToDeltaWithSequenceGapList);
  }
  else
  {
    VisitBioseqsInSep (sep, &ctdd, ConvertBioseqToDeltaWithAlignmentGapList);
  }
  
  AdjustAlignmentsBasedOnGapLocations (ctdd.affected_bioseq_list, 
                                       ctdd.location_list, ctdd.salp, 
                                       dcp->input_entityID);                                      

  ctdd.affected_bioseq_list = ValNodeFree (ctdd.affected_bioseq_list);
  ctdd.location_list = ValNodeFreeData (ctdd.location_list);
  
  if (GetStatus (dcp->adjust_CDS_locations))
  {
    VisitFeaturesInSep (sep, (Pointer)Sequin_GlobalAlign2Seq, AdjustCDSLocationsForUnknownGapsCallback);
    DeleteMarkedObjects (dcp->input_entityID, 0, NULL);
  }  
  
  ObjMgrSetDirtyFlag (dcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dcp->input_entityID, 0, 0);
  Remove (dcp->form);    
  ArrowCursor ();
  Update (); 
}

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

static void ConvertRawToDeltaWithGapLocations (IteM i)
{
  BaseFormPtr        bfp;
  DeltaConversionPtr dcp;
  WindoW             w;
  GrouP              h, p, c;
  ButtoN             b;
  PrompT             p1, p2, p3;
  SeqEntryPtr        sep;
  RecT               r1, r2;
  TagListPtr         tlp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dcp = (DeltaConversionPtr) MemNew (sizeof (DeltaConversionData));
  if (dcp == NULL) return;
  
  dcp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID (dcp->input_entityID);
  if (sep == NULL)
  {
    return;
  }
  
  if (!ContinueWithDeltaSequences (sep))
  {
    return;
  }
  if (!ContinueWithSequencesInAlignments (sep))
  {
    return;
  }

  w = FixedWindow (-50, -33, -10, -10, "Convert Raw Sequence to Delta Sequence", StdCloseWindowProc);
  SetObjectExtra (w, dcp, StdCleanupFormProc);
  dcp->form = (ForM) w;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
 
  p = HiddenGroup (h, 3, 0, NULL);
  p1 = StaticPrompt (p, "Start", 0, dialogTextHeight, programFont, 'c');
  p2 = StaticPrompt (p, "Type", 0, dialogTextHeight, programFont, 'c');
  p3 = StaticPrompt (p, "Length", 0, dialogTextHeight, programFont, 'c');
  dcp->gap_locations = CreateTagListDialogEx (h, 4, 4, 2,
                                              deltaconversionedittypes, 
                                              deltaconversionedit_widths,
                                              deltaconversionedit_alists,
                                              TRUE, FALSE, NULL, DeltaLocToData);
  dcp->coord_grp = NormalGroup (h, 2, 0, "Coordinates", programFont, NULL);
  RadioButton (dcp->coord_grp, "Sequence");
  RadioButton (dcp->coord_grp, "Alignment");
  SetValue (dcp->coord_grp, 1);
  if (!SeqEntryHasAligns (dcp->input_entityID, sep))
  {
    Disable (dcp->coord_grp);
  }
  
  dcp->adjust_CDS_locations = CheckBox (h, "Adjust CDS locations for gaps", NULL);
  SetStatus (dcp->adjust_CDS_locations, TRUE);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoConvertRawToDeltaWithGapLocations);
  SetObjectExtra (b, dcp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dcp->gap_locations,
                              (HANDLE) dcp->coord_grp,
                              (HANDLE) dcp->adjust_CDS_locations,
                              (HANDLE) c,
                              NULL);

  tlp = GetObjectExtra (dcp->gap_locations);
  if (tlp != NULL)
  {
    ObjectRect (tlp->control [0], &r1);
    ObjectRect (p1, &r2);
    r2.left = r1.left;
    r2.right = r1.right;
    SetPosition (p1, &r2);
    
    ObjectRect (tlp->control [1], &r1);
    ObjectRect (p2, &r2);
    r2.left = r1.left;
    r2.right = r1.right;
    SetPosition (p2, &r2);
    
    ObjectRect (tlp->control [2], &r1);
    ObjectRect (p3, &r2);
    r2.left = r1.left;
    r2.right = r1.right;
    SetPosition (p3, &r2);
  }
  RealizeWindow (w);
  
  Show (w);
  Update ();  
}

typedef struct xrefgenedata {
  GeneRefPtr  grp;
  SeqLocPtr   slp;
  Boolean     pseudo;
} XrefGeneData, PNTR XrefGenePtr;

static Boolean XrefToGeneCallback (GatherContextPtr gcp)

{
  GeneRefPtr           grp;
  SeqFeatXrefPtr PNTR  last;
  SeqFeatXrefPtr       next;
  SeqFeatPtr           sfp;
  ValNodePtr PNTR      vnpp;
  XrefGenePtr          xgp;
  SeqFeatXrefPtr       xref;

  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->location == NULL) return TRUE;
  if (sfp->xref == NULL) return TRUE;
  grp = NULL;
  last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;
    if (xref->data.choice == SEQFEAT_GENE) {
      *last = next;
      xref->next = NULL;
      grp = (GeneRefPtr) xref->data.value.ptrvalue;
      xref->data.value.ptrvalue = NULL;
      SeqFeatXrefFree (xref);
    } else {
      last = &(xref->next);
    }
    xref = next;
  }
  if (grp == NULL) return TRUE;
  xgp = MemNew (sizeof (XrefGeneData));
  if (xgp == NULL) return TRUE;
  xgp->grp = grp;
  xgp->slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead,
                           (AsnWriteFunc) SeqLocAsnWrite);
  /* if feature is pseudo, gene created from xref on feature should also be pseudo */
  xgp->pseudo = sfp->pseudo;

  ValNodeAddPointer (vnpp, 0, xgp);
  return TRUE;
}

static void XrefToGene (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  Uint2         entityID = 0;
  GatherScope   gs;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  SeqFeatPtr    sfp;
  SelStructPtr  ssp;
  ValNodePtr    vnp;
  XrefGenePtr   xgp;

  head = NULL;
  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
    if (bfp == NULL) {
      Message (MSG_ERROR, "XrefToGene error");
      return;
    }
    entityID = bfp->input_entityID;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    GatherEntity (bfp->input_entityID, (Pointer) &head, XrefToGeneCallback, &gs);
  } else {
    entityID = ssp->entityID;
    if (ssp->itemtype != OBJ_SEQFEAT) {
      Message (MSG_ERROR, "Feature must be selected");
      return;
    }
    GatherItem (ssp->entityID, ssp->itemID, ssp->itemtype, (Pointer) &head, XrefToGeneCallback);
  }
  if (head == NULL || entityID == 0) return;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    xgp = (XrefGenePtr) vnp->data.ptrvalue;
    if (xgp != NULL) {
      bsp = GetBioseqGivenSeqLoc (xgp->slp, entityID);
      if (bsp != NULL) {
        nsep = SeqMgrGetSeqEntryForData (bsp);
        if (! ExtendGene (xgp->grp, nsep, xgp->slp)) {
          sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
          if (sfp != NULL) {
            sfp->data.value.ptrvalue = (Pointer) xgp->grp;
            xgp->grp = NULL;
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = xgp->slp;
            xgp->slp = NULL;
            sfp->pseudo = xgp->pseudo;
          }
        }
      }
      GeneRefFree (xgp->grp);
      SeqLocFree (xgp->slp);
    }
  }
  ValNodeFreeData (head);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void GeneToXrefCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BaseFormPtr     bfp;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  SeqFeatPtr      sfp;
  SeqAnnotPtr     sap;
  SeqFeatXrefPtr  xref;

  if (sep == NULL) return;
  bfp = (BaseFormPtr) mydata;
  if (bfp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice != SEQFEAT_GENE) {
          FindGeneAndProtForCDS (bfp->input_entityID, sfp, &gene, NULL);
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (grp != NULL) {
              xref = SeqFeatXrefNew ();
              if (xref != NULL) {
                xref->data.choice = SEQFEAT_GENE;
                xref->data.value.ptrvalue = AsnIoMemCopy ((Pointer) grp,
                                                          (AsnReadFunc) GeneRefAsnRead,
                                                          (AsnWriteFunc) GeneRefAsnWrite);
                xref->next = sfp->xref;
                sfp->xref = xref;
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void GeneToXref (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    Message (MSG_ERROR, "GeneToXref error");
    return;
  }
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, (Pointer) bfp, GeneToXrefCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct featurefieldforgenechoice
{
  Int4    first_choice;
  Int4    second_choice;
  CharPtr field_txt;  
} FeatureFieldForGeneChoiceData, PNTR FeatureFieldForGeneChoicePtr;

static FeatureFieldForGeneChoicePtr FeatureFieldForGeneChoiceFree (FeatureFieldForGeneChoicePtr fcp)
{
  if (fcp != NULL)
  {
    fcp->field_txt = MemFree (fcp->field_txt);
    fcp = MemFree (fcp);
  }
  return fcp;
}

typedef struct featuretogene
{
  FEATURE_FORM_BLOCK

  DialoG *gene_src_dlg_list;
  Int4   num_choices;
  ValNodePtr feature_choices;
  
  TexT   label_txt;
  PopuP  genechoice;
  GrouP  qual_caps_grp;
  DialoG feature_dlg;
  DialoG filter_grp;
  DialoG accept_cancel;
  
  ButtoN single_interval_btn;
  ButtoN selected_features_only_btn;
  
  Boolean single_interval;
  FeatureFieldForGeneChoicePtr fcp;
  
} FeatureToGeneData, PNTR FeatureToGenePtr;

static Int4 GetGeneSrcDlgIndex (FeatureToGenePtr fgp)
{
  ValNodePtr vnp, check_vnp;
  Int4       rval = -1, i;
  
  if (fgp == NULL) return -1;

  vnp = DialogToPointer (fgp->feature_dlg);
  if (vnp == NULL)
  {
    return -1;
  }

  check_vnp = fgp->feature_choices;
  i = 0;
  while (check_vnp != NULL && check_vnp->choice != vnp->choice)
  {
    check_vnp = check_vnp->next;
    i++;    
  }
  if (check_vnp != NULL)
  {
    rval = i;
  }
  ValNodeFreeData (vnp);
  return rval;
}

static void EnableFeatureToGeneControls (Pointer userdata)
{
  Int4             i;
  FeatureToGenePtr fgp;
  FeatureFieldForGeneChoicePtr fcp;
  
  fgp = (FeatureToGenePtr) userdata;
  if (fgp == NULL) return;
  Disable (fgp->qual_caps_grp);
  Disable (fgp->genechoice);
  
  i = GetGeneSrcDlgIndex (fgp);
  if (i >= 0)
  {
    fcp = (FeatureFieldForGeneChoicePtr) DialogToPointer (fgp->gene_src_dlg_list [i]);
    if (fcp != NULL)
    {
      if (fcp->first_choice > 1 || !StringHasNoText (fcp->field_txt))
      {
        Enable (fgp->qual_caps_grp);
        Enable (fgp->genechoice);
      }     
      fcp = FeatureFieldForGeneChoiceFree (fcp);
    }
  }
}

static void EnableFeatureToGeneControlsPopup (PopuP p)
{
  FeatureToGenePtr fgp;
  
  fgp = (FeatureToGenePtr) GetObjectExtra (p);
  EnableFeatureToGeneControls (fgp);
}

static void EnableFeatureToGeneControlsText (TexT t)
{
  FeatureToGenePtr fgp;
  
  fgp = (FeatureToGenePtr) GetObjectExtra (t);
  EnableFeatureToGeneControls (fgp);
  
}

typedef struct featurefieldforgenedlg 
{
  DIALOG_MESSAGE_BLOCK
  PopuP first_choice_popup;
  PopuP second_choice_popup;
  TexT  label_txt;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} FeatureFieldForGeneDlgData, PNTR FeatureFieldForGeneDlgPtr;

static void ResetFeatureFieldForGeneDlg (FeatureFieldForGeneDlgPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  
  SetValue (dlg->first_choice_popup, 1);
  SafeSetValue (dlg->second_choice_popup, 1);
  SetTitle (dlg->label_txt, "");
}

static void FeatureFieldForGeneToDialog (DialoG d, Pointer userdata)
{
  FeatureFieldForGeneDlgPtr    dlg;
  FeatureFieldForGeneChoicePtr data;

  dlg = (FeatureFieldForGeneDlgPtr) GetObjectExtra (d);
  data = (FeatureFieldForGeneChoicePtr) userdata;
  if (dlg == NULL)
  {
    return;
  }
  ResetFeatureFieldForGeneDlg (dlg);
  if (data == NULL)
  {
    return;
  }
  if (data->first_choice > 0)
  {
    SetValue (dlg->first_choice_popup, data->first_choice);
  }
  if (data->second_choice > 0)
  {
    SafeSetValue (dlg->second_choice_popup, data->second_choice);
  }
  if (!StringHasNoText (data->field_txt))
  {
    SetTitle (dlg->label_txt, data->field_txt);
  }
}

static Pointer FeatureFieldForGeneFromDialog (DialoG d)
{
  FeatureFieldForGeneDlgPtr    dlg;
  FeatureFieldForGeneChoicePtr data;

  dlg = (FeatureFieldForGeneDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  data = (FeatureFieldForGeneChoicePtr) MemNew (sizeof (FeatureFieldForGeneChoiceData));
  if (data != NULL)
  {
    data->first_choice = GetValue (dlg->first_choice_popup);
    if (dlg->second_choice_popup != NULL)
    {
      data->second_choice = GetValue (dlg->second_choice_popup);
    }
    else
    {
      data->second_choice = 0;
    }
    if (TextHasNoText (dlg->label_txt))
    {
      data->field_txt = NULL;
    }
    else
    {
      data->field_txt = SaveStringFromText (dlg->label_txt);
    }
  }
  return data;
}

static void FeatureFieldForGeneChange (FeatureFieldForGeneDlgPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  if (GetValue (dlg->first_choice_popup) > 1)
  {
    SafeEnable (dlg->second_choice_popup);
  }
  else
  {
    SafeDisable (dlg->second_choice_popup);
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void FeatureFieldForGeneChangePopup (PopuP p)
{
  FeatureFieldForGeneDlgPtr  dlg;

  dlg = (FeatureFieldForGeneDlgPtr) GetObjectExtra (p);
  FeatureFieldForGeneChange (dlg);  
}

static void FeatureFieldForGeneChangeText (TexT t)
{
  FeatureFieldForGeneDlgPtr  dlg;

  dlg = (FeatureFieldForGeneDlgPtr) GetObjectExtra (t);
  FeatureFieldForGeneChange (dlg);  
}

static PopuP MakeFieldChoicePopup (GrouP g, Int4 featdef_choice, Pointer extradata)
{
  PopuP p;
  
  p = PopupList (g, TRUE, FeatureFieldForGeneChangePopup);
  SetObjectExtra (p, extradata, NULL);
  
  if (featdef_choice == FEATDEF_CDS
      || featdef_choice == FEATDEF_tRNA
      || featdef_choice == FEATDEF_rRNA
      || featdef_choice == FEATDEF_misc_RNA
      || featdef_choice == FEATDEF_mRNA)
  {  
    PopupItem (p, "None");
    PopupItem (p, "Comment");
    PopupItem (p, "Product");
  }
  else
  {
    PopupItem (p, "None");
    PopupItem (p, "Comment");
  }
  SetValue (p, 1);
  return p;	
}

static DialoG 
FeatureFieldForGeneDialog 
(GrouP parent, 
 Int4  featdef_choice,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  FeatureFieldForGeneDlgPtr  dlg;
  GrouP                     p;
  
  dlg = (FeatureFieldForGeneDlgPtr) MemNew (sizeof (FeatureFieldForGeneDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = NormalGroup (parent, 2, 0, "Select qualifier to use in gene", NULL, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatureFieldForGeneToDialog;
  dlg->fromdialog = FeatureFieldForGeneFromDialog;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  StaticPrompt (p, "1st Choice", 0, 0, programFont, 'c');
  dlg->first_choice_popup = MakeFieldChoicePopup (p, featdef_choice, dlg);
  
  if (featdef_choice == FEATDEF_CDS
      || featdef_choice == FEATDEF_tRNA
      || featdef_choice == FEATDEF_rRNA
      || featdef_choice == FEATDEF_misc_RNA
      || featdef_choice == FEATDEF_mRNA)
  {
    StaticPrompt (p, "2nd Choice", 0, 0, programFont, 'c');
    dlg->second_choice_popup = MakeFieldChoicePopup (p, featdef_choice, dlg);
    Disable (dlg->second_choice_popup);
  }
  else
  {
    dlg->second_choice_popup = NULL;
  }

  StaticPrompt (p, "Use this string:", 0, 0, programFont, 'c');
  dlg->label_txt = DialogText (p, "", 10, FeatureFieldForGeneChangeText);
  SetObjectExtra (dlg->label_txt, dlg, NULL);
    
  return (DialoG) p;  
}

static CharPtr GetGeneSrcChoice (SeqFeatPtr sfp, Int4 choice)
{
  CharPtr           src_txt = NULL;
  RnaRefPtr         rrp;
  SeqMgrFeatContext fcontext;
  
  if (sfp == NULL || choice < 2)
  {
    return NULL;
  }
  
  if (choice == 2 && !StringHasNoText (sfp->comment))
  {
    src_txt = StringSave (sfp->comment);
  }
  else if (choice == 3)
  {
    if (sfp->idx.subtype == FEATDEF_tRNA)
    {
      sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL,
                                     0, 0, sfp, &fcontext);
      if (!StringHasNoText (fcontext.label))
      {
        src_txt = (CharPtr) MemNew (StringLen (fcontext.label) + 6);
        if (src_txt != NULL)
        {
          sprintf (src_txt, "tRNA-%s", fcontext.label);
        }
      }
    }
    else if (sfp->idx.subtype == FEATDEF_CDS
             || sfp->idx.subtype == FEATDEF_mRNA)
    {
      sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL,
                                     0, 0, sfp, &fcontext);
      if (!StringHasNoText (fcontext.label))
      {
        src_txt = StringSave (fcontext.label);
      }
    }
    else if (sfp->data.choice == SEQFEAT_RNA)
    {
      rrp = (RnaRefPtr) (sfp->data.value.ptrvalue);
      if (rrp != NULL && rrp->ext.choice == 1
          && !StringHasNoText (rrp->ext.value.ptrvalue))
      {
        src_txt = StringSave (rrp->ext.value.ptrvalue);
      }
    }
  }
  return src_txt;  
}

static CharPtr GetGeneSrc (SeqFeatPtr sfp, FeatureFieldForGeneChoicePtr fcp)
{
  CharPtr src_txt = NULL;
  
  if (sfp == NULL || fcp == NULL)
  {
    return NULL;
  }
  
  src_txt = GetGeneSrcChoice (sfp, fcp->first_choice);
  if (src_txt == NULL)
  {
    src_txt = GetGeneSrcChoice (sfp, fcp->second_choice);
  }
  if (src_txt == NULL && !StringHasNoText (fcp->field_txt))
  {
    src_txt = StringSave (fcp->field_txt);
  }
  return src_txt;  
}

static Boolean IsBioseqmRNA (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mip;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL || sdp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol == 3)
  {
    return TRUE;
  }
  return FALSE;
}

static void FeatureToGeneCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  SeqFeatPtr         gene, overlap_gene;
  GeneRefPtr         grp;
  Boolean            partial5, partial3;
  FeatureToGenePtr   fgp;
  CharPtr            gene_val = NULL;
  Int4               i;
  BioseqPtr          bsp;
  CharPtr            cp;
  SeqLocPtr          gene_location;

  if (sfp == NULL || userdata == NULL) return;
  fgp = (FeatureToGenePtr) userdata;
  
  if (SeqMgrGetGeneXref (sfp) != NULL) return;
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  
  gene_location = SeqLocMerge (bsp, sfp->location, NULL, fgp->single_interval, FALSE, FALSE);
  if (gene_location == NULL) return;
  
  overlap_gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (overlap_gene != NULL && SeqLocCompare (gene_location, overlap_gene->location) == SLC_A_EQ_B)
  {
    gene_location = SeqLocFree (gene_location);
    return;
  }
  
  if (IsBioseqmRNA (bsp))
  {
    gene_location = SeqLocFree (gene_location);
    gene_location = SeqLocIntNew (0, bsp->length - 1, 
                                  SeqLocStrand (sfp->location),
                                  SeqIdFindBest (bsp->id, 0));
  }
  
  gene_val = GetGeneSrc (sfp, fgp->fcp);

  if (gene_val != NULL)
  {
  	/* apply capitalization */
	  switch (GetValue (fgp->qual_caps_grp))
  	{
  	  case 1:
  	    /* do nothing, leave capitalization as is */
  	    break;
  	  case 2:
  	    /* capitalize first letter */
  	    gene_val [0] = toupper (gene_val[0]);
  	    break;
  	  case 3:
  	  	/* capitalize all letters */
  	  	for (cp = gene_val; *cp != 0; cp++)
  	  	{
  	  	  *cp = toupper (*cp);
  	  	}
  	  	break;
  	}
  }

  gene = SeqFeatNew ();
  if (gene == NULL) return;
    
  grp = GeneRefNew ();
  if (grp == NULL) return;
  i = GetValue (fgp->genechoice);
  switch (i)
  {
  	case 1:
  	  grp->locus = gene_val;
  	  break;
  	case 2:
  	  grp->locus_tag = gene_val;
  	  break;
  	case 3:
  	  grp->desc = gene_val;
  	  break;
  	case 4:
  	  grp->allele = gene_val;
  	  break;
  	default:
  	  gene->comment = gene_val;
  	  break;
  }
    
  gene->data.choice = SEQFEAT_GENE;
  gene->data.value.ptrvalue = (Pointer) grp;
  gene->location = gene_location;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  SetSeqLocPartial (gene->location, partial5, partial3);
  gene->partial = sfp->partial;
  gene->next = sfp->next;
  sfp->next = gene;
}

static Boolean FeatureToGeneAccept (Pointer userdata)
{
  FeatureToGenePtr  fgp;
  SeqEntryPtr       sep;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  FilterSetPtr      fsp;
  Int4              i, featdef_choice;
  ValNodePtr        vnp;
  
  fgp = (FeatureToGenePtr) userdata;
  if (fgp == NULL) return FALSE;
  
  vnp = DialogToPointer (fgp->feature_dlg);
  if (vnp == NULL)
  {
    return FALSE;
  }
  featdef_choice = vnp->choice;
  vnp = ValNodeFreeData (vnp);
  
  i = GetGeneSrcDlgIndex (fgp);
  if (i < 0)
  {
    return FALSE;
  }

  fgp->fcp = DialogToPointer (fgp->gene_src_dlg_list [i]);
  
  fsp = DialogToPointer (fgp->filter_grp);
  
  sep = GetTopSeqEntryForEntityID (fgp->input_entityID);
  if (sep == NULL) return FALSE;
  fgp->single_interval = GetStatus (fgp->single_interval_btn);
  if (GetStatus (fgp->selected_features_only_btn))
  {
    sel = ObjMgrGetSelected ();
    while (sel != NULL)
    {
      if (sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if (sfp != NULL)
        {
          FeatureToGeneCallback (sfp, fgp, NULL);
        }
      }
      sel = sel->next;
    }
  }
  else
  {
    OperateOnSeqEntryConstrainedObjects (sep, fsp, FeatureToGeneCallback,
                                         NULL, 0, featdef_choice, 0, fgp);
  }
  
  FilterSetFree (fsp);
  fgp->fcp = FeatureFieldForGeneChoiceFree (fgp->fcp);
  
  ObjMgrSetDirtyFlag (fgp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, fgp->input_entityID, 0, 0);
  return TRUE;
}

static void ChangeGeneFeatureSelection (Pointer userdata)
{
  FeatureToGenePtr  fgp;
  Int4              i;
  
  fgp = (FeatureToGenePtr) userdata;
  if (fgp == NULL)
  {
    return;
  }
  
  for (i = 0; i < fgp->num_choices; i++)
  {
    Hide (fgp->gene_src_dlg_list [i]);
  }

  i = GetGeneSrcDlgIndex (fgp);
  
  if (i < 0)
  {
    DisableAcceptCancelDialogAccept (fgp->accept_cancel);
  }
  else
  {
    EnableAcceptCancelDialogAccept (fgp->accept_cancel);
    Show (fgp->gene_src_dlg_list [i]);
  }
  EnableFeatureToGeneControls (fgp);
}

static void CleanupFeatureToGeneForm (GraphiC g, VoidPtr data)
{
  FeatureToGenePtr fgp;

  fgp = (FeatureToGenePtr) data;
  if (fgp != NULL) 
  {
    fgp->feature_choices = ValNodeFreeData (fgp->feature_choices);
    MemFree (fgp);
  }
}

static void FeatureToGeneClearText (Pointer userdata)
{
  FeatureToGenePtr             fgp;
  Int4                         j;
  FilterSetPtr                 fsp;
  FeatureFieldForGeneChoicePtr fcp;

  fgp = (FeatureToGenePtr) userdata;
  if (fgp == NULL)
  {
    return;
  }
  
  for (j = 0; j < fgp->num_choices; j++)
  {
    fcp = (FeatureFieldForGeneChoicePtr) DialogToPointer (fgp->gene_src_dlg_list [j]);
    if (fcp != NULL)
    {
      fcp->field_txt = MemFree (fcp->field_txt);
    }
    PointerToDialog (fgp->gene_src_dlg_list [j], fcp);
    FeatureFieldForGeneChoiceFree (fcp);
  }
  
  fsp = (FilterSetPtr) DialogToPointer (fgp->filter_grp);
  FilterSetClearText (fsp);
  PointerToDialog (fgp->filter_grp, fsp);
  FilterSetFree (fsp);
  EnableFeatureToGeneControls (fgp);
}

static void FeatureToGeneClear (Pointer userdata)
{
  FeatureToGenePtr  fgp;
  Int4              j;

  fgp = (FeatureToGenePtr) userdata;
  if (fgp == NULL)
  {
    return;
  }
  
  PointerToDialog (fgp->feature_dlg, NULL);
  for (j = 0; j < fgp->num_choices; j++)
  {
    PointerToDialog (fgp->gene_src_dlg_list [j], NULL);
  }
  PointerToDialog (fgp->filter_grp, NULL);
  
  SetValue (fgp->qual_caps_grp, 1);
  SetValue (fgp->genechoice, 1);
  SetStatus (fgp->single_interval_btn, TRUE);
  SetStatus (fgp->selected_features_only_btn, FALSE);

  ChangeGeneFeatureSelection (fgp);
  EnableFeatureToGeneControls (fgp);
}

static void FeatureToGene (IteM i)
{
  BaseFormPtr       bfp;
  FeatureToGenePtr  fgp;
  WindoW            w;
  Int4              j;
  GrouP             h, k, m, n;
  ValNodePtr        vnp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    Message (MSG_ERROR, "Feature to Gene error");
    return;
  }
  
  fgp = (FeatureToGenePtr) MemNew (sizeof (FeatureToGeneData));
  if (fgp == NULL) return;
  fgp->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Feature to Gene", StdCloseWindowProc);
  SetObjectExtra (w, fgp, CleanupFeatureToGeneForm);
  fgp->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  fgp->feature_dlg = FeatureSelectionDialogEx (h, FALSE, sep, ChangeGeneFeatureSelection, fgp);
  k = HiddenGroup (h, 0, 0, NULL);

  fgp->feature_choices = BuildFeatureDialogList (TRUE, sep);
  fgp->num_choices = ValNodeLen (fgp->feature_choices);
  fgp->gene_src_dlg_list = (DialoG *) MemNew (fgp->num_choices * sizeof (DialoG));
  
  for (j = 0, vnp = fgp->feature_choices;
       j < fgp->num_choices && vnp != NULL;
       j++, vnp = vnp->next)
  {
    fgp->gene_src_dlg_list [j] = FeatureFieldForGeneDialog (k, vnp->choice,
                                                            EnableFeatureToGeneControls,
                                                            fgp);
  }
  
  fgp->qual_caps_grp = NormalGroup (h, 3, 0,
                     "Capitalization for gene", NULL, NULL);
  RadioButton (fgp->qual_caps_grp, "As is");
  RadioButton (fgp->qual_caps_grp, "Capitalize first initial");
  RadioButton (fgp->qual_caps_grp, "Capitalize all");
  SetValue (fgp->qual_caps_grp, 1);
  Disable (fgp->qual_caps_grp);
  m = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (m, "Select gene qualifier to populate", 0,0, programFont, 'c');
  fgp->genechoice = PopupList (m, TRUE, NULL);
  PopupItem (fgp->genechoice, "locus");
  PopupItem (fgp->genechoice, "locus_tag");
  PopupItem (fgp->genechoice, "gene description");
  PopupItem (fgp->genechoice, "allele");
  PopupItem (fgp->genechoice, "gene comment");
  SetValue (fgp->genechoice, 1);
  Disable (fgp->genechoice);
  
  n = HiddenGroup (h, 2, 0, NULL);
  fgp->single_interval_btn = CheckBox (n, "Create gene with single interval location", NULL);
  SetStatus (fgp->single_interval_btn, TRUE);
  fgp->selected_features_only_btn = CheckBox (n, "Only create genes for selected features", NULL);
  SetStatus (fgp->selected_features_only_btn, FALSE);
  
  fgp->filter_grp = FilterGroup (h, TRUE, FALSE, FALSE, FALSE, FALSE, "Where feature field contains");

  fgp->accept_cancel = AcceptCancelDialog (h, FeatureToGeneAccept, NULL, FeatureToGeneClear, FeatureToGeneClearText, (Pointer)fgp, w);

  AlignObjects (ALIGN_CENTER, (HANDLE) fgp->feature_dlg,
                              (HANDLE) k, 
                              (HANDLE) fgp->qual_caps_grp, 
                              (HANDLE) m,
                              (HANDLE) n,
                              (HANDLE) fgp->filter_grp,
                              (HANDLE) fgp->accept_cancel, 
                              NULL);
                
  RealizeWindow (w);
  Show (w);
  Update ();
  ChangeGeneFeatureSelection (fgp);
}


typedef struct featuretocdsform
{
  FEATURE_FORM_BLOCK
  PopuP   feature_choice;
  
  DialoG  gene_field;
  DialoG  exon_field;
  DialoG  mrna_field;
  GrouP   qual_caps_grp;
  TexT    append_text;
  ButtoN  accept_btn;
  ButtoN  fuse_multiple_btn;
  
  Int4       featdef_choice;
  ValNodePtr gene_field_list;
  ValNodePtr mrna_field_list;
  ValNodePtr exon_field_list;
  CharPtr    append_string;
  Int4       caps_choice;
  Boolean    fuse_multiple;
  Boolean    makemRNA;
  SeqEntryPtr sep;
  LogInfoData lid;
  
} FeatureToCDSFormData, PNTR FeatureToCDSFormPtr;

#define FEATURE_TO_CDS_GENE 1
#define FEATURE_TO_CDS_MRNA 2
#define FEATURE_TO_CDS_EXON 3

static void EnableFeatureToCDSControls (Pointer data)
{
  FeatureToCDSFormPtr    fp;
  Int4                   feature_choice;
  Boolean                need_name_string = FALSE;
  ValNodePtr             missing = NULL;
  CharPtr                cp;

  fp = (FeatureToCDSFormPtr) data;
  if (fp == NULL)
  {
    return;
  }
  Enable (fp->accept_btn);
  feature_choice = GetValue (fp->feature_choice);
  if (feature_choice == FEATURE_TO_CDS_GENE)
  {
    missing = TestDialog (fp->gene_field);
    if (missing != NULL)
    {
      need_name_string = TRUE;
      ValNodeFree (missing);
    }
  }
  else if (feature_choice == FEATURE_TO_CDS_EXON)
  {
    missing = DialogToPointer (fp->exon_field);
    if (missing == NULL)
    {
      need_name_string = TRUE;
    }
    ValNodeFree (missing);
  }
  else if (feature_choice == FEATURE_TO_CDS_MRNA)
  {
    missing = DialogToPointer (fp->mrna_field);
    if (missing == NULL)
    {
      need_name_string = TRUE;
    }
    ValNodeFree (missing);
  }
  
  if (need_name_string)
  {
    Disable (fp->qual_caps_grp);
    cp = SaveStringFromText (fp->append_text);
    if (StringHasNoText (cp))
    {
      Disable (fp->accept_btn);
    }
    MemFree (cp);
  }
  else
  {
    Enable (fp->qual_caps_grp);
  }
}

static void EnableFeatureToCDSControlsText (TexT t)
{
  FeatureToCDSFormPtr    fp;

  fp = (FeatureToCDSFormPtr) GetObjectExtra (t);
  EnableFeatureToCDSControls (fp);
}

static void ChangeFeatureToCDS (PopuP p)
{
  FeatureToCDSFormPtr    fp;
  Int4                   feature_choice;

  fp = (FeatureToCDSFormPtr) GetObjectExtra (p);
  if (fp != NULL)
  {
    feature_choice = GetValue (fp->feature_choice);
    if (feature_choice == FEATURE_TO_CDS_GENE)
    {
      Show (fp->gene_field);
      Hide (fp->mrna_field);
      Hide (fp->exon_field);
    }
    else if (feature_choice == FEATURE_TO_CDS_MRNA)
    {
      Show (fp->mrna_field);
      Hide (fp->gene_field);
      Hide (fp->exon_field);
    }
    else
    {
      Show (fp->exon_field);
      Hide (fp->mrna_field);
      Hide (fp->gene_field);
    }
  }
  EnableFeatureToCDSControls (fp);
}

static void CleanupFeatureToCDSPage (GraphiC g, VoidPtr data)
{
  FeatureToCDSFormPtr fp;

  fp = (FeatureToCDSFormPtr) data;
  if (fp != NULL) 
  {
    fp->append_string = MemFree (fp->append_string);
    fp->gene_field_list = ValNodeFree (fp->gene_field_list);
    fp->mrna_field_list = ValNodeFree (fp->mrna_field_list);
    fp->exon_field_list = ValNodeFree (fp->exon_field_list);
    MemFree (fp);
  }
}

static void LogFeatureToCDSNoProtein (FeatureToCDSFormPtr fp, CharPtr feature_label, SeqIdPtr sip)
{
  Char       id [41];

  if (fp == NULL || fp->lid.fp == NULL || sip == NULL) return;

  id [0] = '\0';
  SeqIdWrite (SeqIdFindBest (sip, 0), id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  
  if (StringHasNoText (feature_label))
  {
    fprintf (fp->lid.fp, "No protein added for unlabeled feature on sequence %s\n", id);
  }
  else
  {
    fprintf (fp->lid.fp, "No protein added for %s on sequence %s\n", feature_label, id);
  }
  fp->lid.data_in_log = TRUE;
}

static CharPtr GetNewProteinName (SeqFeatPtr sfp, FeatureToCDSFormPtr fp)
{
  SeqMgrFeatContext context;
  CharPtr           prot_name_src = NULL;
  CharPtr           protName = NULL;
  Int4              prot_len;
  CharPtr           cp;
  
  if (sfp == NULL || fp == NULL || sfp->idx.subtype != fp->featdef_choice)
  {
    return NULL;
  }

  if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL,
                               0, 0, sfp, &context) == NULL) return NULL;
                               
  /* determine the new protein name */
  if (fp->featdef_choice == FEATDEF_GENE)
  {
    prot_name_src = GetGeneFieldString (sfp, fp->gene_field_list, NULL);
  }
  else if (fp->featdef_choice == FEATDEF_mRNA)
  {
    prot_name_src = GetmRNAFieldString (sfp, fp->mrna_field_list, NULL);
  }
  else if (fp->featdef_choice == FEATDEF_exon)
  {
    prot_name_src = GetExonFieldString (sfp, fp->exon_field_list);
  }
  else if (fp->featdef_choice == FEATDEF_CDS)
  {
    prot_name_src = GetCDSFieldString(sfp, fp->mrna_field_list, NULL);
  }
  
  if (StringHasNoText (prot_name_src) && StringHasNoText (fp->append_string))
  {
  	return NULL;
  }

  prot_len = StringLen (prot_name_src) + StringLen (fp->append_string) + 2;
  
  /* add one to length for terminating null */
  protName = (CharPtr) MemNew (prot_len * sizeof (Char));
  if (protName == NULL) return NULL;
  
  if (!StringHasNoText (prot_name_src))
  {
    StringCpy (protName, prot_name_src);
  	/* apply capitalization */
    switch (fp->caps_choice)
  	{
  	  case 1:
  	    /* do nothing, leave capitalization as is */
  	    break;
  	  case 2:
  	    /* capitalize first letter */
  	    protName [0] = toupper (protName[0]);
  	    break;
  	  case 3:
  	  	/* capitalize all letters */
  	  	for (cp = protName; *cp != 0; cp++)
  	  	{
  	  	  *cp = toupper (*cp);
  	  	}
  	  	break;
  	}
  	if (!StringHasNoText (fp->append_string))
  	{
  	  StringCat (protName, " ");
  	}
  }
  prot_name_src = MemFree (prot_name_src);
  if (!StringHasNoText (fp->append_string))
  {
    StringCat (protName, fp->append_string);
  }
  return protName;  
}

static SeqFeatPtr MakeCDS(SeqLocPtr location, CharPtr name_str, Uint2 entityID, SeqEntryPtr sep, LogInfoPtr lip)
{
  SeqFeatPtr        cds, sfp;
  Int2              genCode;
  CdRegionPtr       crp;
  ByteStorePtr      bs;
  BioseqPtr         bsp;
  SeqEntryPtr       psep, nsep, old, topsep;
  MolInfoPtr        mip;
  ProtRefPtr        prp;
  Char              ch;
  ValNodePtr        descr;
  Int2              i;
  Boolean           partial5, partial3;
  CharPtr           prot;
  CharPtr           ptr;
  ValNodePtr        vnp;
  
  cds = SeqFeatNew ();
  
  CheckSeqLocForPartial (location, &partial5, &partial3);

  if (cds == NULL) return NULL;
  cds->data.choice = SEQFEAT_CDREGION;

  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (crp == NULL) {
    return NULL;
  }
 
  cds->data.value.ptrvalue = (Pointer) crp;
  cds->location = SeqLocFree (cds->location);
  cds->location = AsnIoMemCopy ((Pointer) location,
                                (AsnReadFunc) SeqLocAsnRead,
                                (AsnWriteFunc) SeqLocAsnWrite);
  cds->partial = partial5 | partial3;                                
  
  if (!SetBestFrameByLocation (cds)) {
    LogCDSAmbiguousFrame (lip, cds);
  }
 
  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
  if (bs != NULL) {
    prot = BSMerge (bs, NULL);
    bs = BSFree (bs);
    if (prot != NULL) {
      ptr = prot;
      ch = *ptr;
      while (ch != '\0') {
        *ptr = TO_UPPER (ch);
        ptr++;
        ch = *ptr;
      }
      i = (Int2) StringLen (prot);
      if (i > 0 && prot [i - 1] == '*') {
        prot [i - 1] = '\0';
      }
      bs = BSNew (1000);
	  if (bs != NULL) {
        ptr = prot;
        BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
      }
    }
  }
  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (NULL);
  bsp->id = MakeNewProteinSeqId (cds->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);
  psep = SeqEntryNew ();
  if (psep != NULL) {
    psep->choice = 1;
    psep->data.ptrvalue = (Pointer) bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 8;
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    descr = ExtractBioSourceAndPubs (sep);
    if (IS_Bioseq_set (sep))
    {
      topsep = sep;
    }
    else
    {
      topsep = GetBestTopParentForData (entityID, sep->data.ptrvalue);
    }
    AddSeqEntryToSeqEntry (topsep, psep, TRUE);
    nsep = FindNucSeqEntry (sep);
    ReplaceBioSourceAndPubs (sep, descr);
    SetSeqFeatProduct (cds, bsp);
    prp = CreateNewProtRef (name_str, NULL, NULL, NULL);
    sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (sfp != NULL) {
       sfp->data.value.ptrvalue = (Pointer) prp;
    }
  }
  return cds;

}


static SeqFeatPtr MakemRNA (SeqLocPtr location, CharPtr protName, SeqEntryPtr sep)
{
  SeqFeatPtr sfp;
  RnaRefPtr  rrp;
  Boolean    partial5, partial3;
  
  sfp = SeqFeatNew();
  rrp = RnaRefNew();
  rrp->type = RNA_TYPE_mRNA;
  rrp->ext.choice = 1;
  rrp->ext.value.ptrvalue = protName;
  sfp->data.value.ptrvalue = rrp;
  sfp->data.choice = SEQFEAT_RNA;
  
  CheckSeqLocForPartial (location, &partial5, &partial3);
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = AsnIoMemCopy ((Pointer) location,
                                (AsnReadFunc) SeqLocAsnRead,
                                (AsnWriteFunc) SeqLocAsnWrite);
  sfp->partial = partial5 | partial3;
  
  return sfp;
}


static SeqFeatPtr FeatureToCDS (SeqFeatPtr sfp, FeatureToCDSFormPtr fp)

{
  SeqFeatPtr        cds;
  SeqMgrFeatContext context;
  CharPtr           protName = NULL;
  RnaRefPtr         rrp;

  
  if (sfp == NULL || fp == NULL || sfp->idx.subtype != fp->featdef_choice)
  {
    return NULL;
  }

  if (fp->makemRNA) {
    cds = SeqMgrGetOverlappingmRNA (sfp->location, &context);
  } else  {
    cds = SeqMgrGetOverlappingCDS (sfp->location, &context);
  }
  if (cds != NULL) 
  {
    return NULL;	
  }

  if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL,
                               0, 0, sfp, &context) == NULL) 
  {
    return NULL;
  }
                               
  protName = GetNewProteinName (sfp, fp);
  if (protName == NULL)
  {
    LogFeatureToCDSNoProtein (fp, context.label, SeqLocId (sfp->location));
    return NULL;
  }

  /* make sure mRNA product name matches CDS protein name */
  if (fp->featdef_choice == FEATDEF_mRNA && !fp->makemRNA)
  {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL)
    {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave (protName);
    }    
  }
  
  /* need to remove string */
  if (fp->featdef_choice == FEATDEF_GENE)
  {
    RemoveGeneFieldString (sfp, fp->gene_field_list);
  }
  else if (fp->featdef_choice == FEATDEF_mRNA)
  {
    RemovemRNAFieldString (sfp, fp->mrna_field_list);
  }
  else if (fp->featdef_choice == FEATDEF_exon)
  {
    RemoveExonFieldString (sfp, fp->exon_field_list);
  }
  else if (fp->featdef_choice == FEATDEF_CDS)
  {
    RemoveCDSFieldString (sfp, fp->mrna_field_list);
  }
  
  if (fp->makemRNA) {
    cds = MakemRNA (sfp->location, protName, fp->sep);
  } else {
    cds = MakeCDS(sfp->location, protName, fp->input_entityID, fp->sep, &(fp->lid));
  }
  
  cds->partial = sfp->partial;
  cds->next = sfp->next;
  sfp->next = cds;

  return cds;
}

static void 
LogFeatureToCDSMismatchProtein 
(FeatureToCDSFormPtr fp, 
 SeqIdPtr            sip)
{
  Char       id [41];

  if (fp == NULL || fp->lid.fp == NULL || sip == NULL) return;

  id [0] = '\0';
  SeqIdWrite (SeqIdFindBest (sip, 0), id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  
  fprintf (fp->lid.fp, "The protein names generated from features on %s do not match\n", id);
  
  fp->lid.data_in_log = TRUE;
}

extern void LogCDSAmbiguousFrame (LogInfoPtr lip, SeqFeatPtr sfp)
{
  CharPtr loc_str;
  
  if (lip == NULL || lip->fp == NULL || sfp == NULL) {
    return;
  }
  
  loc_str = SeqLocPrintUseBestID (sfp->location);
  fprintf (lip->fp, "Ambiguous frames detected for coding region at %s\n", loc_str);
  loc_str = MemFree (loc_str);
  lip->data_in_log = TRUE;
}

static void FeatureToCDSBioseqCheckCallback (BioseqPtr bsp, Pointer userdata)

{
  FeatureToCDSFormPtr fp;
  CharPtr             last_prot_name = NULL, new_prot_name;
  SeqMgrFeatContext   feature_context, cds_context;
  SeqFeatPtr          sfp, cds;

  if (bsp == NULL || userdata == NULL) return;
  fp = (FeatureToCDSFormPtr) userdata;
  fp->sep = SeqMgrGetSeqEntryForData (bsp);
  
  /* if we are fusing features, find all the features on this Bioseq for
   * our type.
   */
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &feature_context);
  while (sfp != NULL)
  {
    if (sfp->idx.subtype == fp->featdef_choice)
    {
      if (fp->makemRNA) {
        cds = SeqMgrGetOverlappingmRNA (sfp->location, &cds_context);
      } else {
        cds = SeqMgrGetOverlappingCDS (sfp->location, &cds_context);
      }
      if (cds == NULL) 
      {
        new_prot_name = GetNewProteinName (sfp, fp);
        if (!fp->fuse_multiple)
        {
          if (new_prot_name == NULL)
          {
            LogFeatureToCDSNoProtein (fp, feature_context.label, SeqLocId (sfp->location));
          }
        }
        else if (last_prot_name == NULL && new_prot_name != NULL)
        {
          if (new_prot_name == NULL)
          {
            LogFeatureToCDSNoProtein (fp, feature_context.label, SeqLocId (sfp->location));
          }
          else
          {
            last_prot_name = new_prot_name;
          }
        }
        else if (new_prot_name == NULL)
        {
          /* do nothing - will use protein name from prior CDS */ 
        }
        else if (StringCmp (last_prot_name, new_prot_name) == 0)
        {
          /* no conflict */
          MemFree (new_prot_name);
        }
        else
        {
          LogFeatureToCDSMismatchProtein (fp, SeqLocId (sfp->location));
          MemFree (last_prot_name);
          last_prot_name = new_prot_name;
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, fp->featdef_choice, &feature_context);
  }
}

static void FeatureToCDSBioseqActCallback (BioseqPtr bsp, Pointer userdata)

{
  FeatureToCDSFormPtr fp;
  SeqFeatPtr          first = NULL;
  CharPtr             new_prot_name;
  SeqMgrFeatContext   feature_context, cds_context;
  SeqFeatPtr          sfp, cds;
  SeqLocPtr           new_cds_location = NULL, slp, tmp_loc;
  Boolean             first_partial5, first_partial3;
  Int4                first_start_pos, first_stop_pos;
  Boolean             sfp_partial5, sfp_partial3;
  Int4                sfp_start_pos, sfp_stop_pos;

  if (bsp == NULL || userdata == NULL) return;
  fp = (FeatureToCDSFormPtr) userdata;
  fp->sep = SeqMgrGetSeqEntryForData (bsp);
  
  /* if we are fusing features, find all the features on this Bioseq for
   * our type.
   */
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &feature_context);
  while (sfp != NULL)
  {
    if (sfp->idx.subtype == fp->featdef_choice)
    {
      if (fp->makemRNA) {
        cds = SeqMgrGetOverlappingmRNA (sfp->location, &cds_context);
      } else {
        cds = SeqMgrGetOverlappingCDS (sfp->location, &cds_context);
      }
      if (cds == NULL) 
      {
        if (fp->fuse_multiple)
        {
          if (first == NULL)
          {
            new_prot_name = GetNewProteinName (sfp, fp);
            if (new_prot_name != NULL)
            {
              first = sfp;
              new_cds_location = (SeqLocPtr) AsnIoMemCopy ((Pointer) sfp->location,
                                                     (AsnReadFunc) SeqLocAsnRead,
                                                     (AsnWriteFunc) SeqLocAsnWrite);
              MemFree (new_prot_name);
            }
          }
          else
          {
            /* preserve partialness of ends */
            CheckSeqLocForPartial (new_cds_location, &first_partial5, &first_partial3);
            first_start_pos = SeqLocStart (new_cds_location);
            first_stop_pos = SeqLocStop (new_cds_location);
            CheckSeqLocForPartial (sfp->location, &sfp_partial5, &sfp_partial3);
            sfp_start_pos = SeqLocStart (sfp->location);
            sfp_stop_pos = SeqLocStop (sfp->location);
            if (first_start_pos > sfp_start_pos)
            {
              first_partial5 = sfp_partial5;
            }
            if (first_stop_pos < sfp_stop_pos)
            {
              first_partial3 = sfp_partial3;
            }
            
            slp = SeqLocMerge (bsp, sfp->location, new_cds_location, FALSE, TRUE, FALSE);
            new_cds_location = SeqLocFree (new_cds_location);
            new_cds_location = slp;
            SetSeqLocPartial (new_cds_location, first_partial5, first_partial3);            
            first_start_pos = SeqLocStart (new_cds_location);
            first_stop_pos = SeqLocStop (new_cds_location);
          }
        }
        else
        {
          FeatureToCDS (sfp, fp);
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, fp->featdef_choice, &feature_context);
  }
  if (fp->fuse_multiple && first != NULL)
  {
    tmp_loc = first->location;
    first->location = new_cds_location;
    cds = FeatureToCDS (first, fp);
    if (cds != NULL)
    {
      cds->partial = CheckSeqLocForPartial (cds->location, &first_partial5, &first_partial3);
    }
    first->location = tmp_loc;
    new_cds_location = SeqLocFree (new_cds_location);
  }
}

static void FeatureToCDSAccept (ButtoN b)
{
  FeatureToCDSFormPtr    fp;
  Int4                   feature_choice;
  MsgAnswer              ans = ANS_YES;

  fp = (FeatureToCDSFormPtr) GetObjectExtra (b);
  if (fp == NULL)
  {
    return;
  }
  
  feature_choice = GetValue (fp->feature_choice);
  if (feature_choice == FEATURE_TO_CDS_GENE)
  {
    fp->featdef_choice = FEATDEF_GENE;
    fp->gene_field_list = DialogToPointer (fp->gene_field);
  }
  else if (feature_choice == FEATURE_TO_CDS_MRNA)
  {
    if (fp->makemRNA) {
      fp->featdef_choice = FEATDEF_CDS;
    } else {
      fp->featdef_choice = FEATDEF_mRNA;
    }
    fp->mrna_field_list = DialogToPointer (fp->mrna_field);
  }
  else if (feature_choice == FEATURE_TO_CDS_EXON)
  {
    fp->featdef_choice = FEATDEF_exon;
    fp->exon_field_list = DialogToPointer (fp->exon_field);
  }
  else
  {
    return;
  }

  TmpNam (fp->lid.path);
  fp->lid.fp = FileOpen (fp->lid.path, "wb");
  if (fp->lid.fp == NULL) return;
  fp->lid.display_title = "Name Problems";
  
  Hide (fp->form);
  WatchCursor ();
  Update ();

  fp->sep = GetTopSeqEntryForEntityID (fp->input_entityID);
  if (fp->sep == NULL) return;
  
  fp->caps_choice = GetValue (fp->qual_caps_grp);
  
  fp->append_string = JustSaveStringFromText (fp->append_text);
  
  fp->fuse_multiple = GetStatus (fp->fuse_multiple_btn);

  VisitBioseqsInSep (fp->sep, (Pointer) fp, FeatureToCDSBioseqCheckCallback);

  if (ans == ANS_YES)
  {
    fp->sep = GetTopSeqEntryForEntityID (fp->input_entityID);
    VisitBioseqsInSep (fp->sep, (Pointer) fp, FeatureToCDSBioseqActCallback);
  }
  
  CloseLog (&(fp->lid));
    
  ArrowCursor ();
  Update ();

  ObjMgrSetDirtyFlag (fp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, fp->input_entityID, 0, 0);  
  Remove (fp->form); 
}

static void FeatureToCDSormRNA (IteM i, Boolean makemRNA)

{
  BaseFormPtr            bfp;
  FeatureToCDSFormPtr    fp;
  WindoW                 w;
  GrouP                  h, g, m, c;
  CharPtr                title, grp_title;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    Message (MSG_ERROR, "FeatureToCDS error");
    return;
  }
  
  fp = (FeatureToCDSFormPtr) MemNew (sizeof (FeatureToCDSFormData));
  if (fp == NULL) return;
  fp->input_entityID = bfp->input_entityID;
  fp->input_itemID = bfp->input_itemID;
  fp->input_itemtype = bfp->input_itemtype;
  
  fp->makemRNA = makemRNA;
  
  if (makemRNA) {
    title = "Feature to mRNA";
  } else {
    title = "Feature to CDS";
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, fp, CleanupFeatureToCDSPage);
  fp->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  fp->feature_choice = PopupList (h, TRUE, ChangeFeatureToCDS);
  PopupItem (fp->feature_choice, "Gene");
  if (makemRNA) {
    PopupItem (fp->feature_choice, "CDS");
  } else {
    PopupItem (fp->feature_choice, "mRNA");
  }
  PopupItem (fp->feature_choice, "exon");
  SetValue (fp->feature_choice, FEATURE_TO_CDS_GENE);
  SetObjectExtra (fp->feature_choice, fp, NULL);

  if (makemRNA) {
    grp_title = "Source for new mRNA name";
  } else {
    grp_title = "Source for new CDS protein name";
  }
  g = NormalGroup (h, 0, 0, grp_title, NULL, NULL);
  fp->gene_field = FeatureFieldChoiceDialog (g, GeneFieldSelectionDialog, TRUE, EnableFeatureToCDSControls, fp);
  fp->mrna_field = FeatureFieldChoiceDialog (g, MRNAFieldSelectionDialog, TRUE, EnableFeatureToCDSControls, fp);
  fp->exon_field = FeatureFieldChoiceDialog (g, ExonFieldSelectionDialog, TRUE, EnableFeatureToCDSControls, fp);
  Hide (fp->mrna_field);
  Hide (fp->exon_field);
  AlignObjects (ALIGN_CENTER, (HANDLE) fp->gene_field,
                              (HANDLE) fp->mrna_field, 
                              (HANDLE) fp->exon_field, 
                              NULL);
  
  fp->qual_caps_grp = NormalGroup (h, 3, 0,
                     "Capitalization for field in name", NULL, NULL);
  RadioButton (fp->qual_caps_grp, "As is");
  RadioButton (fp->qual_caps_grp, "Capitalize first initial");
  RadioButton (fp->qual_caps_grp, "Capitalize all");
  SetValue (fp->qual_caps_grp, 1);
  
  m = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (m, "Append text to name", 0, 0, programFont, 'c');
  fp->append_text = DialogText (m, "", 14, EnableFeatureToCDSControlsText); 
  SetObjectExtra (fp->append_text, fp, NULL);
  
  if(makemRNA) {
    fp->fuse_multiple_btn = CheckBox (h, "Fuse multiple intervals for new mRNA", NULL);
  } else {
    fp->fuse_multiple_btn = CheckBox (h, "Fuse multiple intervals for new CDS", NULL);
  }
   
  c = HiddenGroup (h, 4, 0, NULL); 
  fp->accept_btn = PushButton (c, "Accept", FeatureToCDSAccept);
  SetObjectExtra (fp->accept_btn, fp, NULL);
  EnableFeatureToCDSControls (fp);
  
  PushButton (c, "Cancel", StdCancelButtonProc);


  AlignObjects (ALIGN_CENTER, (HANDLE) fp->feature_choice,
                              (HANDLE) g,
                              (HANDLE) fp->qual_caps_grp,
                              (HANDLE) m, 
                              (HANDLE) fp->fuse_multiple_btn, 
                              (HANDLE) c, NULL);
                
  RealizeWindow (w);
  Show (w);
  Update ();
  
}

static void FeatureToCds (IteM i)
{
  FeatureToCDSormRNA (i, FALSE);
}


static void FeatureTomRNA (IteM i)
{
  FeatureToCDSormRNA (i, TRUE);
}


static void AddFeatureToBioseq (SeqFeatPtr sfp, BioseqPtr bsp)

{
  SeqFeatPtr   prev;
  SeqAnnotPtr  sap;

  if (sfp == NULL || bsp == NULL) return;
  sap = bsp->annot;
  while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
    sap = sap->next;
  }
  if (sap == NULL) {
    sap = SeqAnnotNew ();
    if (sap != NULL) {
      sap->type = 1;
      sap->next = bsp->annot;
      bsp->annot = sap;
    }
  }
  sap = bsp->annot;
  if (sap != NULL) {
    if (sap->data != NULL) {
      prev = sap->data;
      while (prev->next != NULL) {
        prev = prev->next;
      }
      prev->next = sfp;
    } else {
      sap->data = (Pointer) sfp;
    }
  }
}

static void PackageFeatureCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BaseFormPtr    bfp;
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp = NULL;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  BioseqPtr      target;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  bfp = (BaseFormPtr) mydata;
  if (bfp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        target = GetBioseqGivenSeqLoc (sfp->location, bfp->input_entityID);
        if (target != bsp) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          AddFeatureToBioseq (sfp, target);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void PackageOnParts (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  SeqEntryExplore (sep, bfp, PackageFeatureCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ForcePackageFeatureCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BaseFormPtr    bfp;
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp = NULL;
  SeqAnnotPtr    firstsap;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  BioseqPtr      target;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  bfp = (BaseFormPtr) mydata;
  if (bfp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  firstsap = sap;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        target = GetBioseqGivenSeqLoc (sfp->location, bfp->input_entityID);
        if (target != bsp || sap != firstsap) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          AddFeatureToBioseq (sfp, target);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void ForceRepackage (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  SeqEntryExplore (sep, bfp, ForcePackageFeatureCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static BioseqSetPtr GetParentNPS (BioseqPtr bsp)
{
  BioseqSetPtr  bssp;

  if (bsp == NULL) return NULL;
  if (bsp->idx.parenttype != OBJ_BIOSEQSET) return NULL;
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot && bssp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bssp->idx.parentptr;
  }
  if (bssp->_class == BioseqseqSet_class_nuc_prot) return bssp;
  return NULL;
}

static Boolean NucAndProtNotInSameNPS (BioseqPtr nuc, BioseqPtr prt)
{
  BioseqSetPtr    bssp;

  if (nuc == NULL || prt == NULL) return FALSE;
  bssp = GetParentNPS (nuc);
  if (bssp == NULL) return TRUE;
  if (GetParentNPS (prt) != bssp) return TRUE;
  return FALSE;
}

typedef struct nppack {
  BioseqPtr  nuc;
  BioseqPtr  prt;
} NpPack, PNTR NpPackPtr;

static void PackageProteinsCallback (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqSetPtr       bssp;
  SeqMgrDescContext  context;
  ValNodePtr PNTR    headp;
  MolInfoPtr         mip;
  NpPackPtr          npp;
  BioseqPtr          nuc, prt;
  SeqDescrPtr        sdp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->product == NULL) return;
  nuc = BioseqFindFromSeqLoc (sfp->location);
  prt = BioseqFindFromSeqLoc (sfp->product);
  if (nuc == NULL || prt == NULL) return;

  headp = (ValNodePtr PNTR) userdata;

  /* if CDS location is on genomic nucleotide in genomic product set, bail */

  if (nuc->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) nuc->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      sdp = SeqMgrGetNextDescriptor (nuc, NULL, Seq_descr_molinfo, &context);
      if (sdp != NULL) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip != NULL) {
          if (mip->biomol != MOLECULE_TYPE_MRNA) return;
        }
      }
    }
  }

  if (NucAndProtNotInSameNPS (nuc, prt)) {
    npp = (NpPackPtr) MemNew (sizeof (NpPack));
    if (npp == NULL) return;
    npp->nuc = nuc;
    npp->prt = (BioseqPtr) AsnIoMemCopy ((Pointer) prt,
                                         (AsnReadFunc) BioseqAsnRead,
                                         (AsnWriteFunc) BioseqAsnWrite);
    ValNodeAddPointer (headp, 0, (Pointer) npp);
    prt->idx.deleteme = TRUE;
  }
}

static void PackageOnNucs (IteM i)

{
  BaseFormPtr    bfp;
  ValNodePtr     descr;
  ValNodePtr     head = NULL, vnp;
  NpPackPtr      npp;
  BioseqSetPtr   nps;
  SeqEntryPtr    nsep, psep;
  BioseqPtr      nuc, prt;
  SeqEntryPtr    sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();

  VisitFeaturesInSep (sep, &head, PackageProteinsCallback);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    npp = (NpPackPtr) vnp->data.ptrvalue;
    if (npp == NULL) continue;
    nuc = npp->nuc;
    prt = npp->prt;
    if (nuc == NULL || prt == NULL) continue;
    nsep = nuc->seqentry;
    nps = GetParentNPS (nuc);
    if (nps != NULL) {
      nsep = nps->seqentry;
    }
    psep = SeqEntryNew ();
    if (nsep == NULL || psep == NULL) continue;
    psep->choice = 1;
    psep->data.ptrvalue = (Pointer) prt;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) prt, psep);
    descr = ExtractBioSourceAndPubs (nsep);
    AddSeqEntryToSeqEntry (nsep, psep, TRUE);
    ReplaceBioSourceAndPubs (nsep, descr);
    AssignIDsInEntity (bfp->input_entityID, 0, NULL);
  }
  ValNodeFreeData (head);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern EnumFieldAssoc  enum_bond_alist [];
extern EnumFieldAssoc  enum_site_alist [];

static Boolean            alistBoxUp;
static Boolean            alistBoxAccept;
static PopuP              alistBoxPopup;

static void AcceptAlistMessage (ButtoN b)

{
  alistBoxAccept = TRUE;
  alistBoxUp = FALSE;
}

static void CancelAlistMessage (ButtoN b)

{
  alistBoxAccept = FALSE;
  alistBoxUp = FALSE;
}

extern Boolean AlistMessage (EnumFieldAssocPtr al, UIEnumPtr val, UIEnum dflt, CharPtr mssg)

{
  GrouP   c;
  GrouP   g;
  PrompT  ppt;
  WindoW  w;

  if (al == NULL || val == NULL) return FALSE;
  alistBoxUp = TRUE;
  alistBoxAccept = FALSE;
  w = ModalWindow (-50, -20, -20, -20, NULL);
  if (w != NULL) {
    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 10, 10);
    ppt = StaticPrompt (g, mssg, 0, 0, programFont, 'c');
    alistBoxPopup = PopupList (g, TRUE, NULL);
    InitEnumPopup (alistBoxPopup, al, NULL);
    SetEnumPopup (alistBoxPopup, al, dflt);
    c = HiddenGroup (g, 2, 0, NULL);
    DefaultButton (c, "Accept", AcceptAlistMessage);
    PushButton (c, "Cancel", CancelAlistMessage);
    AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) alistBoxPopup, (HANDLE) c, NULL);
    RealizeWindow (w);
    Show (w);
    Select (w);
    while (alistBoxUp) {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (! GetEnumPopup (alistBoxPopup, al, val)) {
      alistBoxAccept = FALSE;
    }
    Remove (w);
  }
  return alistBoxAccept;
} 

static void DefaultMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

typedef struct multtransformdata {
  FEATURE_FORM_BLOCK

  Uint1          type;
  LisT           fromfeatlist;
  PopuP          siteorbondlist;
  Uint1          feattype;
  Int4           siteorbondtype;
} MultTransFormData, PNTR MultTransFormPtr;

static void ProcessMultiTrans (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  SeqFeatPtr        cds;
  MultTransFormPtr  mfp;
  SeqFeatPtr        newsfp;
  BioseqPtr         pbsp;
  SeqEntryPtr       psep;
  SeqAnnotPtr       sap;
  SeqBondPtr        sbp;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;
  SeqLocPtr         slp;
  SeqPntPtr         spp;
  Uint1             subtype;

  mfp = (MultTransFormPtr) mydata;
  if (sep == NULL || mfp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        subtype = FindFeatDefType (sfp);
        if (subtype == mfp->feattype) {
          cds = SeqMgrGetOverlappingCDS (sfp->location, NULL);
          if (cds != NULL) {
            slp = dnaLoc_to_aaLoc (cds, sfp->location, TRUE, NULL, FALSE);
            if (slp != NULL) {
              pbsp = BioseqFindFromSeqLoc (slp);
              if (pbsp != NULL) {
                psep = SeqMgrGetSeqEntryForData (pbsp);
                if (psep != NULL) {
                  newsfp = SeqFeatNew ();
                  if (newsfp != NULL) {
                    newsfp->data.choice = mfp->type;
                    if (mfp->type == SEQFEAT_BOND || mfp->type == SEQFEAT_SITE) {
                      newsfp->data.value.intvalue = (Int4) mfp->siteorbondtype;
                    } else if (mfp->type == SEQFEAT_REGION) {
                      newsfp->data.value.ptrvalue = sfp->comment;
                      sfp->comment = NULL;
                    }
                    CreateNewFeature (psep, NULL, mfp->type, newsfp);
                    newsfp->location = slp;
                    newsfp->comment = sfp->comment;
                    sfp->comment = NULL;
                    SeqFeatDataFree (&sfp->data);            /* free duplicate feature data */
                    sfp->data.choice = SEQFEAT_GENE;         /* fake as empty gene, to be removed */
                    sfp->data.value.ptrvalue = GeneRefNew ();
                    if (mfp->type == SEQFEAT_BOND) {
                      if (slp->choice != SEQLOC_BOND) {
                        sip = SeqLocId (slp);
                        if (sip != NULL) {
                          sbp = SeqBondNew ();
                          if (sbp != NULL) {
                            slp = ValNodeNew (NULL);
                            if (slp != NULL) {
                              slp->choice = SEQLOC_BOND;
                              slp->data.ptrvalue = (Pointer) sbp;
                              spp = SeqPntNew ();
                              if (spp != NULL) {
                                spp->strand = SeqLocStrand (newsfp->location);
                                spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
                                spp->point = SeqLocStart (newsfp->location);
                                sbp->a = spp;
                              }
                              spp = SeqPntNew ();
                              if (spp != NULL) {
                                spp->strand = SeqLocStrand (newsfp->location);
                                spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
                                spp->point = SeqLocStop (newsfp->location);
                                sbp->b = spp;
                              }
                              newsfp->location = SeqLocFree (newsfp->location);
                              newsfp->location = slp;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void DoMultiTransform (ButtoN b)

{
  UIEnum            enumval;
  Int2              feattype;
  MultTransFormPtr  mfp;
  SeqEntryPtr       sep;
  Int2              siteorbondtype = 0;
  UIEnumPtr         uptr;

  mfp = (MultTransFormPtr) GetObjectExtra (b);
  if (mfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (mfp->input_entityID);
  if (sep == NULL) return;
  Hide (mfp->form);
  WatchCursor ();
  Update ();
  uptr = (UIEnumPtr) DialogToPointer ((DialoG) mfp->fromfeatlist);
  if (uptr != NULL) {
    enumval = *uptr;
    feattype = (Int2) enumval;
    uptr = (UIEnumPtr) DialogToPointer ((DialoG) mfp->siteorbondlist);
    if (uptr != NULL) {
      enumval = *uptr;
      siteorbondtype = (Int2) enumval;
    }
    if (feattype > 0) {
      mfp->feattype = (Uint1) feattype;
      mfp->siteorbondtype = (Int4) siteorbondtype;
      if (siteorbondtype > 0 || mfp->type == SEQFEAT_REGION) {
        SeqEntryExplore (sep, (Pointer) mfp, ProcessMultiTrans);
      }
    }
  }
  ArrowCursor ();
  Update ();
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  ObjMgrSetDirtyFlag (mfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mfp->input_entityID, 0, 0);
  Remove (mfp->form);
}

static void MultiTransformProc (IteM i, Uint1 targettype)

{
  EnumFieldAssocPtr  ap;
  ButtoN             b;
  BaseFormPtr        bfp;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  GrouP              k;
  Int2               listHeight;
  MultTransFormPtr   mfp;
  EnumFieldAssocPtr  realalist;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title = "?";
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  mfp = (MultTransFormPtr) MemNew (sizeof (MultTransFormData));
  if (mfp == NULL) return;

  if (targettype == SEQFEAT_BOND) {
    title = "Transform to Bond";
  } else if (targettype == SEQFEAT_SITE) {
    title = "Transform to Site";
  } else if (targettype == SEQFEAT_REGION) {
    title = "Transform to Region";
  }

  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, mfp, StdCleanupFormProc);
  mfp->form = (ForM) w;
  mfp->formmessage = DefaultMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    mfp->appmessage = sepp->handleMessages;
  }

  mfp->type = targettype;

  mfp->input_entityID = bfp->input_entityID;
  mfp->input_itemID = bfp->input_itemID;
  mfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);

  ap = import_featdef_alist (TRUE, FALSE, TRUE);
  realalist = ap;
  ap++;
  StaticPrompt (g, "From Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  mfp->fromfeatlist = CreateEnumListDialog (g, 10, listHeight, NULL, ap, 0, NULL);
  FreeEnumFieldAlist (realalist);

  k = NULL;
  if (targettype != SEQFEAT_REGION) {
    k = HiddenGroup (h, 0, 2, NULL);

    if (targettype == SEQFEAT_BOND) {
      StaticPrompt (k, "Select type of bond", 0, 0, programFont, 'c');
      mfp->siteorbondlist = CreateEnumPopupDialog (k, TRUE, NULL, enum_bond_alist, 0, NULL);
    } else if (targettype == SEQFEAT_SITE) {
      StaticPrompt (k, "Select type of site", 0, 0, programFont, 'c');
      mfp->siteorbondlist = CreateEnumPopupDialog (k, TRUE, NULL, enum_site_alist, 0, NULL);
    }
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoMultiTransform);
  SetObjectExtra (b, mfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) k, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static Boolean UnlinkPubDesc (GatherContextPtr gcp)

{
  ValNodePtr       sdp;
  ValNodePtr PNTR  vnpp;

  if (gcp->thistype != OBJ_SEQDESC) return TRUE;
  sdp = (ValNodePtr) gcp->thisitem;
  if (sdp == NULL || sdp->choice != Seq_descr_pub) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;
  *(gcp->prevlink) = sdp->next;
  sdp->next = NULL;
  *vnpp = sdp;
  return TRUE;
}

static void AddConvertedDescFeat (SeqEntryPtr sep, ValNodePtr vnp)
{
  BioseqSetPtr  bssp;
  BioseqPtr     bsp;
  SeqFeatPtr    sfp;

  if (sep == NULL || vnp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
	for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      AddConvertedDescFeat (sep, vnp);
	}
  } else {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_aa (bsp->mol)) return;
    sfp = CreateNewFeature (sep, NULL, SEQFEAT_PUB, NULL);
    if (sfp != NULL) {
      sfp->data.value.ptrvalue = AsnIoMemCopy ((Pointer) vnp->data.ptrvalue,
                                               (AsnReadFunc) PubdescAsnRead,
                                               (AsnWriteFunc) PubdescAsnWrite);
    }
  }
}

static void CommonConvertDescToFeat (BaseFormPtr bfp, Boolean pub, Boolean src, Boolean com, CharPtr findstring)

{
  SeqEntryPtr   sep;
  SelStructPtr  sel;
  ValNodePtr    vnp;
  Boolean       asked_about_prop = FALSE;
  Boolean       propagate_descriptors = FALSE;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sel = ObjMgrGetSelected ();
  if (pub && sel != NULL) {
    if (sel->entityID == bfp->input_entityID &&
        sel->next == NULL && sel->itemtype == OBJ_SEQDESC) {
      vnp = NULL;
      sep = GetBestTopParentForItemID (sel->entityID, sel->itemID, sel->itemtype);
      if (sep != NULL) {
        /* unlink changes itemID, so sep must be determined first */
        GatherItem (sel->entityID, sel->itemID, sel->itemtype, (Pointer) &vnp, UnlinkPubDesc);
        if (vnp != NULL) {
          AddConvertedDescFeat (sep, vnp);
          ValNodeFree (vnp);
        }
      }
    }
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    return;
  }
  if (ConvertPubSrcComDescsToFeats (sep, pub, src, com, FALSE, &asked_about_prop, &propagate_descriptors, findstring)) {
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  }
}


static void CommonConvertDescToFeatMenuItem (IteM i, Boolean pub, Boolean src, Boolean com)
{
  BaseFormPtr   bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  CommonConvertDescToFeat (bfp, pub, src, com, NULL);
}

static void ConvertPubDescToFeat (IteM i)

{
  CommonConvertDescToFeatMenuItem (i, TRUE, FALSE, FALSE);
}

static void ConvertSrcDescToFeat (IteM i)

{
  CommonConvertDescToFeatMenuItem (i, FALSE, TRUE, FALSE);
}

static void ConvertComDescToFeat (IteM i)

{
  CommonConvertDescToFeatMenuItem (i, FALSE, FALSE, TRUE);
}

typedef struct convertpubdescdata {
  DESCRIPTOR_FORM_BLOCK
  BaseFormPtr bfp;
  TexT    findthis;
  Char    findString [255];
} ConvertPubDescData, PNTR ConvertPubDescPtr;


static void DoConvertPubDescStringConstraint (ButtoN b)
{
  ConvertPubDescPtr cpdp;

  cpdp = (ConvertPubDescPtr) GetObjectExtra (b);
  if (cpdp == NULL) return;
  Hide (cpdp->form);

  GetTitle (cpdp->findthis, cpdp->findString, sizeof (cpdp->findString) - 1);
  if (StringHasNoText (cpdp->findString)) {
    CommonConvertDescToFeat (cpdp->bfp, TRUE, FALSE, FALSE, NULL);
  } else {
	CommonConvertDescToFeat (cpdp->bfp, TRUE, FALSE, FALSE, cpdp->findString);
  }

  ObjMgrSetDirtyFlag (cpdp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cpdp->input_entityID, 0, 0);
  Remove (cpdp->form);
}

static void CreateConvertPubDescStringConstraintDialogX (IteM i)
{
  BaseFormPtr       bfp;
  ConvertPubDescPtr cpdp;
  WindoW            w;
  GrouP             h, g, c;
  ButtoN            b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cpdp = (ConvertPubDescPtr) MemNew (sizeof (ConvertPubDescData));
  if (cpdp == NULL) return;
  cpdp->bfp = bfp;

  w = FixedWindow (-50, -33, -10, -10, "Convert Pub Descriptors", StdCloseWindowProc);
  SetObjectExtra (w, cpdp, StdCleanupFormProc);
  cpdp->form = (ForM) w;
  cpdp->formmessage = NULL;

  cpdp->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, 0, 2, NULL);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
  cpdp->findthis = DialogText (g, "", 14, NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoConvertPubDescStringConstraint);
  SetObjectExtra (b, cpdp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct feattodesc {
  Boolean pub;
  Boolean src;
  Boolean com;
} FeatToDescData, PNTR FeatToDescPtr;

static void ConvertOneFeatToDesc (SeqFeatPtr sfp, Pointer userdata)
{
  SeqDescPtr sdp;
  SeqEntryPtr sep;
  SeqIdPtr    sip;
  BioseqPtr   bsp;
  ValNode     vn;
  FeatToDescPtr fdp;
  BioSourcePtr  biop;
  SubSourcePtr  ssp, last_ssp = NULL;
  CharPtr       new_note;
  
  if (sfp == NULL || (fdp = (FeatToDescPtr) userdata) == NULL
      || (sfp->data.choice == SEQFEAT_PUB && !fdp->pub)
      || (sfp->data.choice == SEQFEAT_BIOSRC && !fdp->src)
      || (sfp->data.choice == SEQFEAT_COMMENT && !fdp->com)
      || (sfp->data.choice != SEQFEAT_PUB
          && sfp->data.choice != SEQFEAT_BIOSRC
          && sfp->data.choice != SEQFEAT_COMMENT)) {
    return;
  }

  /* perform only if feature is full-length */  
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return;
  
  sep = GetBestTopParentForItemID (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT);
  if (sep == NULL) return;
  
  if (sfp->data.choice == SEQFEAT_PUB) {
    sdp = CreateNewDescriptor (sep, Seq_descr_pub);
    sdp->data.ptrvalue = AsnIoMemCopy ((Pointer) sfp->data.value.ptrvalue,
                                       (AsnReadFunc) PubdescAsnRead,
                                       (AsnWriteFunc) PubdescAsnWrite);
  } else if (sfp->data.choice == SEQFEAT_BIOSRC) {
    sdp = CreateNewDescriptor (sep, Seq_descr_source);
    sdp->data.ptrvalue = AsnIoMemCopy ((Pointer) sfp->data.value.ptrvalue,
                                       (AsnReadFunc) BioSourceAsnRead,
                                       (AsnWriteFunc) BioSourceAsnWrite);
    if (!StringHasNoText (sfp->comment)) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != 255) {
        last_ssp = ssp;
        ssp = ssp->next;
      }
      if (ssp == NULL) {
        ssp = SubSourceNew ();
        ssp->subtype = 255;
        ssp->name = StringSave (sfp->comment);
        if (last_ssp == NULL) {
          biop->subtype = ssp;
        } else {
          last_ssp->next = ssp;
        }
      } else if (StringHasNoText (ssp->name)) {
        ssp->name = MemFree (ssp->name);
        ssp->name = StringSave (sfp->comment);
      } else {
        new_note = (CharPtr) MemNew (sizeof(Char) * (StringLen (ssp->name) + StringLen (sfp->comment) + 2));
        sprintf (new_note, "%s;%s", ssp->name, sfp->comment);
        ssp->name = MemFree (ssp->name);
        ssp->name = new_note;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_COMMENT) {
    sdp = CreateNewDescriptor (sep, Seq_descr_comment);
    sdp->data.ptrvalue = StringSave (sfp->comment);
  }
  
  sfp->idx.deleteme = TRUE;  
}


static void CommonConvertFeatToDesc (IteM i, Boolean pub, Boolean src, Boolean com)

{
  SeqEntryPtr   sep;
  SelStructPtr  sel;
  SeqFeatPtr    sfp;
  FeatToDescData fdd;
  BaseFormPtr   bfp;
  SeqMgrFeatContext fcontext;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  fdd.pub = pub;
  fdd.src = src;
  fdd.com = com;
  
  sel = ObjMgrGetSelected ();  
  if (sel != NULL) {
    while (sel != NULL) {    
      if (sel->entityID == bfp->input_entityID && sel->itemtype == OBJ_SEQFEAT) {
        sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        ConvertOneFeatToDesc (sfp, &fdd);
      }
      sel = sel->next;
    }
  } else {
    VisitFeaturesInSep (sep, &fdd, ConvertOneFeatToDesc);
  }
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void ConvertPubFeatToDesc (IteM i)
{
  CommonConvertFeatToDesc (i, TRUE, FALSE, FALSE);
}


static void ConvertSrcFeatToDesc (IteM i)
{
  CommonConvertFeatToDesc (i, FALSE, TRUE, FALSE);
}


static void ConvertCommentFeatToDesc (IteM i)
{
  CommonConvertFeatToDesc (i, FALSE, FALSE, TRUE);
}



static void MoveSecondProtName (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr  nxt;
  ProtRefPtr  prp;
  ValNodePtr  vnp;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL || prp->desc != NULL) return;
  vnp = prp->name;
  if (vnp == NULL) return;
  nxt = vnp->next;
  if (nxt == NULL || nxt->next != NULL) return;
  vnp->next = NULL;
  prp->desc = nxt->data.ptrvalue;
  nxt->data.ptrvalue = NULL;
  ValNodeFree (nxt);
}

static void SecondProtNameToDesc (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, MoveSecondProtName);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void MoveProtDesc (SeqFeatPtr sfp, Pointer userdata)

{
  ProtRefPtr  prp;
  CharPtr     str;
  ValNodePtr  vnp;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL || prp->desc == NULL) return;
  vnp = prp->name;
  if (vnp == NULL) return;
  str = prp->desc;
  prp->desc = NULL;
  ValNodeAddStr (&(prp->name), 0, str);
}

static void ProtDescToSecondName (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, MoveProtDesc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct ecrepdata {
  CharPtr  before;
  CharPtr  after;
} EcRepData, PNTR EcRepPtr;

static ValNodePtr     ec_rep_list = NULL;
static EcRepPtr PNTR  ec_rep_data = NULL;
static Int4           ec_rep_len = 0;

static int LIBCALLBACK SortVnpByEcBefore (VoidPtr ptr1, VoidPtr ptr2)

{
  EcRepPtr    erp1, erp2;
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  erp1 = (EcRepPtr) vnp1->data.ptrvalue;
  erp2 = (EcRepPtr) vnp2->data.ptrvalue;
  if (erp1 == NULL || erp2 == NULL) return 0;
  str1 = erp1->before;
  str2 = erp2->before;
  if (str1 == NULL || str2 == NULL) return 0;
  return StringCmp (str1, str2);
}

static void SetupECReplacementTable (CharPtr file)

{
  EcRepPtr    erp;
  FileCache   fc;
  FILE        *fp = NULL;
  Int4        i;
  ValNodePtr  last = NULL;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  ValNodePtr  vnp;

  if (ec_rep_data != NULL) return;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            erp = (EcRepPtr) MemNew (sizeof (EcRepData));
            if (erp != NULL) {
              erp->before = StringSave (str);
              erp->after = StringSave (ptr);
              vnp = ValNodeAddPointer (&last, 0, (Pointer) erp);
              if (ec_rep_list == NULL) {
                ec_rep_list = vnp;
              }
              last = vnp;
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
      ec_rep_len = ValNodeLen (ec_rep_list);
      if (ec_rep_len > 0) {
        ec_rep_list = ValNodeSort (ec_rep_list, SortVnpByEcBefore);
        ec_rep_data = (EcRepPtr PNTR) MemNew (sizeof (EcRepPtr) * (ec_rep_len + 1));
        if (ec_rep_data != NULL) {
          for (vnp = ec_rep_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
            erp = (EcRepPtr) vnp->data.ptrvalue;
            ec_rep_data [i] = erp;
          }
        }
      }
    }
  }
}

static void UpdateProtEC (SeqFeatPtr sfp, Pointer userdata)

{
  EcRepPtr    erp;
  Int4        L, R, mid;
  ProtRefPtr  prp;
  CharPtr     str;
  ValNodePtr  vnp;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL || prp->ec == NULL) return;
  for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    L = 0;
    R = ec_rep_len - 1;
    while (L < R) {
      mid = (L + R) / 2;
      erp = ec_rep_data [(int) mid];
      if (erp != NULL && StringCmp (erp->before, str) < 0) {
        L = mid + 1;
      } else {
        R = mid;
      }
    }
    erp = ec_rep_data [(int) R];
    if (erp != NULL && StringCmp (erp->before, str) == 0) {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = StringSave (erp->after);
    }
  }
}

static void UpdateECnumbers (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SetupECReplacementTable ("ecnum_replaced.txt");
  if (ec_rep_data != NULL && ec_rep_len > 0) {
    VisitFeaturesInSep (sep, NULL, UpdateProtEC);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct icrepdata {
  CharPtr  code;
  CharPtr  name;
} IcRepData, PNTR IcRepPtr;

static ValNodePtr     ic_rep_list = NULL;
static IcRepPtr PNTR  ic_rep_data = NULL;
static Int4           ic_rep_len = 0;

static int LIBCALLBACK SortVnpByInstituteName (VoidPtr ptr1, VoidPtr ptr2)

{
  int         compare;
  IcRepPtr    irp1, irp2;
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  irp1 = (IcRepPtr) vnp1->data.ptrvalue;
  irp2 = (IcRepPtr) vnp2->data.ptrvalue;
  if (irp1 == NULL || irp2 == NULL) return 0;
  str1 = irp1->name;
  str2 = irp2->name;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringCmp (str1, str2);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }
  str1 = irp1->code;
  str2 = irp2->code;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringCmp (str1, str2);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }
  return 0;
}

static void SetupInstituteCodeNameTable (void)

{
  FileCache   fc;
  CharPtr     file = "institution_codes.txt";
  FILE        *fp = NULL;
  Int4        i;
  IcRepPtr    irp;
  ValNodePtr  last = NULL;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  ValNodePtr  vnp;

  if (ic_rep_data != NULL) return;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            ptr = StringChr (ptr, '\t');
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
              irp = (IcRepPtr) MemNew (sizeof (IcRepData));
              if (irp != NULL) {
                TrimSpacesAroundString (str);
                TrimSpacesAroundString (ptr);
                irp->code = StringSave (str);
                irp->name = StringSave (ptr);
                vnp = ValNodeAddPointer (&last, 0, (Pointer) irp);
                if (ic_rep_list == NULL) {
                  ic_rep_list = vnp;
                }
                last = vnp;
              }
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
      ic_rep_len = ValNodeLen (ic_rep_list);
      if (ic_rep_len > 0) {
        ic_rep_list = ValNodeSort (ic_rep_list, SortVnpByInstituteName);
        ic_rep_data = (IcRepPtr PNTR) MemNew (sizeof (IcRepPtr) * (ic_rep_len + 1));
        if (ic_rep_data != NULL) {
          for (vnp = ic_rep_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
            irp = (IcRepPtr) vnp->data.ptrvalue;
            ic_rep_data [i] = irp;
          }
        }
      }
    }
  }
}

static CharPtr CheckForInstitutionFullName (CharPtr name, BoolPtr ambigP)

{
  CharPtr   code = NULL;
  IcRepPtr  irp;
  Int4      L, R, mid;

  if (ambigP != NULL) {
    *ambigP = FALSE;
  }
  if (StringHasNoText (name)) return NULL;
  L = 0;
  R = ic_rep_len - 1;
  while (L < R) {
    mid = (L + R) / 2;
    irp = ic_rep_data [(int) mid];
    if (irp != NULL && StringCmp (irp->name, name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  irp = ic_rep_data [(int) R];
  if (irp != NULL && StringCmp (irp->name, name) == 0) {
    code = irp->code;
    /* need to check for uniqueness */
    if (R < ic_rep_len - 1) {
      irp = ic_rep_data [(int) R + 1];
      if (irp != NULL && StringCmp (irp->name, name) == 0) {
        if (ambigP != NULL) {
          *ambigP = TRUE;
        }
        return NULL;
      }
    }
  }
  return code;
}

static void ConvertInstitutionName (OrgModPtr mod)

{
  Boolean  ambig = FALSE;
  Char     buf [512];
  CharPtr  code = NULL, inst = NULL, id = NULL, coll = NULL;

  if (mod == NULL) return;

  StringNCpy_0 (buf, mod->subname, sizeof (buf));
  if (StringChr (buf, ':') == NULL) return;
  if (! ParseStructuredVoucher (buf, &inst, &id)) return;
  if (inst == NULL) return;

  if (VoucherInstitutionIsValid (inst)) return;

  /* ignore personal collections */
  if (StringNICmp (inst, "personal", 8) == 0) return;

  coll = StringChr (inst, ':');
  if (coll != NULL) {
    *coll = '\0';
    if (VoucherInstitutionIsValid (inst)) {
      /* DNA is a valid collection for any institution (using bio_material) */
      if (StringCmp (coll, "DNA") != 0) return;
    }
    *coll = ':';
  }

  code = CheckForInstitutionFullName (inst, &ambig);
  if (StringHasNoText (code)) return;
  if (ambig) {
    Message (MSG_POST, "Institution '%s' maps to multiple codes", inst);
    return;
  }

  mod->subname = MemFree (mod->subname);
  mod->subname = StringSave (code);
}

static void InstitutionNameToCode (BioSourcePtr biop, Pointer userdata)

{
  OrgNamePtr  onp;
  OrgModPtr   omp;
  OrgRefPtr   orp;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_specimen_voucher ||
        omp->subtype == ORGMOD_culture_collection ||
        omp->subtype == ORGMOD_bio_material) {
      ConvertInstitutionName (omp);
    }
  }
}

static void UpdateVoucherName (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SetupInstituteCodeNameTable ();
  if (ic_rep_data != NULL && ic_rep_len > 0) {
    VisitBioSourcesInSep (sep, NULL, InstitutionNameToCode);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void ParseGoTermsFromUserFieldsItem (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ParseGoTermsFromFields (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* HasInternalStops () -- Checks to see if a given protein has         */
/*                        any internal stop codons.                    */
/*                                                                     */
/*   Returns : TRUE -- If internal stop codons are found.              */
/*             FALSE -- If internal stop codons are NOT found.         */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean HasInternalStops (SeqFeatPtr sfp)
{
  ByteStorePtr  bs;
  CharPtr       protStr;
  CharPtr       stopStr;
  Boolean       partial5, partial3;
  Boolean       rval = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }
  /* Get the protein sequence in ASCII chars */

  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (NULL == bs)
    return TRUE;
  protStr = BSMerge (bs, NULL);
  bs = BSFree (bs);

  /* If the is a "*" in the str in any position  */
  /* except at the end, then we have an internal */
  /* stop codon.                                 */
  
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  stopStr = StringStr (protStr, "*");
  if (stopStr != NULL && (partial3 || StringLen (stopStr) > 1)) {
    rval = TRUE;
  }
  MemFree (protStr);
  return rval;
}


static Boolean HasBadStartCodon (SeqFeatPtr sfp)
{
  Boolean partial5, partial3;
  ByteStorePtr  bs;
  CharPtr       protStr;
  Boolean       rval = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5) {
    return FALSE;
  }
 
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs == NULL) {
    rval = FALSE;
  } else {
    protStr = BSMerge (bs, NULL);
    bs = BSFree (bs);
    if (StringNCmp (protStr, "M", 1) != 0) {
      rval = TRUE;
    }
    protStr = MemFree (protStr);
  }
  return rval;
}


static Boolean HasBadStopCodon (SeqFeatPtr sfp)
{
  Boolean partial5, partial3;
  ByteStorePtr  bs;
  CharPtr       protStr;
  Boolean       rval = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial3) {
    return FALSE;
  }
 
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs == NULL) {
    rval = FALSE;
  } else {
    protStr = BSMerge (bs, NULL);
    bs = BSFree (bs);
    if (protStr[StringLen (protStr) - 1] != '*') {
      rval = TRUE;
    }
    protStr = MemFree (protStr);
  }
  return rval;
}


extern Boolean 
ConvertOneCDSToMiscFeat 
(SeqFeatPtr sfp,
 Boolean viral,
 Boolean must_have_stops,
 CDSConversionOptsPtr opts)
{
  CdRegionPtr          cdrp;
  ImpFeatPtr           ifp;
  CharPtr              noteStr;
  BioseqPtr            protBsp;
  SeqMgrFeatContext    protContext, mRNAcontext;
  CharPtr              protName = NULL;
  SeqFeatPtr           protSfp, mRNA, gene;
  ProtRefPtr           prp;
  ValNodePtr           vnp;
  Int4                 note_len = 0;
  CharPtr              viral_fmt = "nonfunctional %s due to mutation";
  CharPtr              similar_fmt = "similar to %s";

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->product == NULL
      || must_have_stops && !HasInternalStops (sfp))
  {
    return FALSE;
  }
      
  /* Get the CD region part of the feature, and */
  /* the associated protein bioseq.             */

  cdrp = (CdRegionPtr) sfp->data.value.ptrvalue;
  protBsp = BioseqFindFromSeqLoc (sfp->product);

  if (protBsp == NULL) return FALSE;

  /* remove the overlapping mRNA feature if requested */
  if (opts != NULL && opts->remove_mrna) {
    mRNA = SeqMgrGetOverlappingmRNA (sfp->location, &mRNAcontext);
    if (mRNA != NULL) {
      mRNA->idx.deleteme = TRUE;
    }
  }

  /* remove the overlapping gene feature if requested */
  if (opts != NULL && opts->remove_gene) {
    gene = SeqMgrGetOverlappingGene (sfp->location, &mRNAcontext);
    if (gene != NULL) {
      gene->idx.deleteme = TRUE;
    }
  }


  /* Convert the CDS feature to a misc_feat */

  CdRegionFree (cdrp);
  sfp->data.value.ptrvalue = NULL;

  ifp = ImpFeatNew ();
  if (NULL == ifp) return FALSE;
  ifp->key = StringSave ("misc_feature");

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  if (! IsFeatInGPS (sfp) || (opts != NULL && opts->remove_transcript_id))
  {
    sfp->product = SeqLocFree (sfp->product);
  }

  /* Add a name key to the misc_feature */

  protSfp = SeqMgrGetBestProteinFeature (protBsp, &protContext);
  if (protSfp != NULL)
  {
    prp = (ProtRefPtr) protSfp->data.value.ptrvalue;

    if (prp != NULL) 
    {
      note_len = StringLen (sfp->comment) + StringLen (prp->desc) + 5;
          
      vnp = prp->name;
      if (NULL != vnp)
      {
        protName = (CharPtr) vnp->data.ptrvalue;
        if (NULL != protName) 
        {
          if (viral) {
            note_len += StringLen (viral_fmt) + StringLen (protName);
          } else {   
            note_len += StringLen (similar_fmt) + StringLen (protName);
          }
        }
      }  
      noteStr = (CharPtr) MemNew (sizeof (Char) * note_len);
        
      if (NULL != protName) {
        if (viral) {
          sprintf (noteStr, viral_fmt, protName);
        } else {
          sprintf (noteStr, similar_fmt, protName);
        }
      }
      if (!StringHasNoText (prp->desc)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, prp->desc);
      }
      if (!StringHasNoText (sfp->comment)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, sfp->comment);
      }
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = noteStr;
    }
  }

  /* Delete the protein Bioseq that */
  /* the CDS points to.             */

  protBsp->idx.deleteme = TRUE;

  /* set the subtype to zero so that it will be reindexed */
  sfp->idx.subtype = 0;
  return TRUE;
}


extern void ConvertCDSToMiscFeat (SeqFeatPtr sfp, Pointer userdata)
{
  CDStoMiscFeatPtr     cmfp;
  Boolean              viral = FALSE;
  Boolean              must_have_stops = TRUE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->product == NULL)
  {
    return;
  }
  
  if (userdata != NULL) {
    cmfp = (CDStoMiscFeatPtr) userdata;
    viral = cmfp->viral;
    must_have_stops = cmfp->must_have_stops;
  }

  ConvertOneCDSToMiscFeat (sfp, viral, must_have_stops, cmfp->opts);
}

static void ChangeOnlyPseudoOption (ButtoN b)
{
  GrouP nonpseudo_opts;

  nonpseudo_opts = (GrouP) GetObjectExtra (b);
  if (nonpseudo_opts != NULL) {
    if (GetStatus (b)) {
      Disable (nonpseudo_opts);
    } else {
      Enable (nonpseudo_opts);
    }
  }
}


extern CDSConversionOptsPtr 
GetCDSConversionOptions (Boolean all_are_pseudo, Boolean any_pseudo, Boolean any_gps, BoolPtr cancel)
{
  ModalAcceptCancelData acd;
  WindoW                w;
  GrouP                 h, g, c, nonpseudo_opts;
  ButtoN                b, remove_mrna_btn, remove_gene_btn, remove_transcript_btn = NULL, only_pseudo_btn = NULL;
  CDSConversionOptsPtr opts;

  if (cancel == NULL) return NULL;

  *cancel = FALSE;
  opts = (CDSConversionOptsPtr) MemNew (sizeof(CDSConversionOptsData));

  opts->all_are_pseudo = all_are_pseudo;  

  if (opts->all_are_pseudo)
  {
    opts->remove_mrna = FALSE;
    opts->remove_gene = FALSE;
    opts->only_pseudo = FALSE;
    opts->remove_transcript_id = FALSE;
    return opts;
  }

  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = NormalGroup (h, 0, 10, "CDS Conversion Options", programFont, NULL);
  SetGroupSpacing (g, 10, 10);
  
  if (any_pseudo) 
  {
    only_pseudo_btn = CheckBox (g, "Only Convert Pseudo CDS", ChangeOnlyPseudoOption);
    SetStatus (only_pseudo_btn, FALSE);
  }

  nonpseudo_opts = HiddenGroup (g, 0, 3, NULL);
  remove_mrna_btn = CheckBox (nonpseudo_opts, "Remove Overlapping mRNAs", NULL);
  SetStatus (remove_mrna_btn, FALSE);
  remove_gene_btn = CheckBox (nonpseudo_opts, "Remove Overlapping Genes", NULL);
  SetStatus (remove_gene_btn, FALSE);
  if (any_gps)
  {
    remove_transcript_btn = CheckBox (nonpseudo_opts, "Remove Transcript ID", NULL);
    SetStatus (remove_transcript_btn, FALSE);
  }

  if (only_pseudo_btn != NULL) 
  {
    SetObjectExtra (only_pseudo_btn, nonpseudo_opts, NULL);
  }
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = DefaultButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) nonpseudo_opts, (HANDLE) c, (HANDLE) only_pseudo_btn, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.accepted)
  {
    if (only_pseudo_btn == NULL || !GetStatus (only_pseudo_btn)) 
    {
      opts->only_pseudo = FALSE;
      opts->remove_mrna = GetStatus (remove_mrna_btn);
      opts->remove_gene = GetStatus (remove_gene_btn);
      opts->remove_transcript_id = GetStatus (remove_transcript_btn);
    }
    else
    {
      opts->only_pseudo = TRUE;
      opts->remove_mrna = FALSE;
      opts->remove_gene = FALSE;
      opts->remove_transcript_id = FALSE;
    }    
  }
  else
  {
    *cancel = TRUE;
    opts = MemFree (opts);
  }

  Remove (w);
  return opts;
}


typedef struct convertcdsoptions {
  Uint2 featdef_to;
  Boolean keep_cds;
  Boolean require_internal_stop;
  Boolean require_bad_start;
  Boolean require_bad_stop;
  Boolean require_pseudo;
  Boolean require_or;
  CharPtr note_fmt;
  Boolean remove_mrna;
  Boolean remove_gene;
  Boolean remove_transcript_id;
  StringConstraintPtr string_constraint;
  StringConstraintPtr id_constraint;
} ConvertCDSOptionsData, PNTR ConvertCDSOptionsPtr;

static const CharPtr kViralNoteFmt = "nonfunctional %s due to mutation";
static const CharPtr kSimilarNoteFmt = "similar to %s";

static ConvertCDSOptionsPtr ConvertCDSOptionsNew ()
{
  ConvertCDSOptionsPtr c;

  c = (ConvertCDSOptionsPtr) MemNew (sizeof (ConvertCDSOptionsData));
  c->featdef_to = FEATDEF_misc_feature;
  c->keep_cds = FALSE;
  c->require_internal_stop = FALSE;
  c->require_bad_start = FALSE;
  c->require_bad_stop = FALSE;
  c->require_pseudo = FALSE;
  c->require_or = FALSE; 
  c->note_fmt = StringSave (kSimilarNoteFmt);
  c->remove_mrna = FALSE;
  c->remove_gene = FALSE;
  c->remove_transcript_id = FALSE;
  c->string_constraint = NULL;
  c->id_constraint = NULL;
  return c;
}


static ConvertCDSOptionsPtr ConvertCDSOptionsFree (ConvertCDSOptionsPtr c)
{
  if (c != NULL) {
    c->note_fmt = MemFree (c->note_fmt);
    c->string_constraint = StringConstraintFree (c->string_constraint);
    c->id_constraint = StringConstraintFree (c->id_constraint);
    c = MemFree (c);
  }
  return c;
}


static ValNodePtr ConstraintChoiceSetFromStringConstraintAndIDStringConstraint (StringConstraintPtr scp, StringConstraintPtr id)
{
  ValNodePtr constraint = NULL;
  SequenceConstraintPtr seq;

  if (!IsStringConstraintEmpty (scp)) {
    ValNodeAddPointer (&constraint, ConstraintChoice_string, AsnIoMemCopy (scp, (AsnReadFunc) StringConstraintAsnRead, (AsnWriteFunc) StringConstraintAsnWrite));
  }
  if (!IsStringConstraintEmpty (id)) {
    seq = SequenceConstraintNew ();
    seq->id = AsnIoMemCopy (id, (AsnReadFunc) StringConstraintAsnRead, (AsnWriteFunc) StringConstraintAsnWrite);
    ValNodeAddPointer (&constraint, ConstraintChoice_sequence, seq);
  }
  return constraint;
}


static void ConvertOneCDS (SeqFeatPtr sfp, Pointer data)
{
  ConvertCDSOptionsPtr c;
  SeqFeatPtr           gene = NULL;
  SeqMgrFeatContext    context;
  ValNodePtr           constraint;
  SeqFeatPtr           new_sfp = NULL, mrna, protSfp;
  GeneRefPtr           grp;
  ProtRefPtr           prp = NULL;
  ImpFeatPtr           ifp;
  CharPtr              prot_name = NULL, prot_desc = NULL, new_comment = NULL, prot_list = NULL;
  Int4                 comment_len;
  BioseqPtr            protBsp = NULL;
  ValNodePtr           vnp;
  Boolean              any_false = FALSE, any_true = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || data == NULL) {
    return;
  }

  c = (ConvertCDSOptionsPtr) data;

  /* determine if original feature meets requirements */
  
  if (!c->require_pseudo && !c->require_internal_stop && !c->require_bad_start && !c->require_bad_stop) {
    any_true = TRUE;
  } else {
    if (c->require_pseudo) {
      if (!sfp->pseudo 
        && ((gene = SeqMgrGetOverlappingGene (sfp->location, &context)) == NULL || !gene->pseudo)) {
        any_false = TRUE;
      } else {
        any_true = TRUE;
      }
    }
    if (c->require_internal_stop) {
      if (HasInternalStops(sfp)) {
        any_true = TRUE;
      } else {
        any_false = TRUE;
      }
    }
    if (c->require_bad_start) {
      if (HasBadStartCodon(sfp)) {
        any_true = TRUE;
      } else {
        any_false = TRUE;
      }
    }
    if (c->require_bad_stop) {
      if (HasBadStopCodon (sfp)) {
        any_true = TRUE;
      } else {
        any_false = TRUE;
      }
    }
  }
  if (!any_true) {
    return;
  } else if (!c->require_or && any_false) {
    return;
  }

  constraint = ConstraintChoiceSetFromStringConstraintAndIDStringConstraint (c->string_constraint, c->id_constraint);
  if (constraint != NULL && !DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, constraint)) {
    constraint = ConstraintChoiceSetFree (constraint);
    return;
  }
  constraint = ConstraintChoiceSetFree (constraint);

  /* remove mrna if requested */
  if (c->remove_mrna) {
    mrna = SeqMgrGetOverlappingmRNA (sfp->location, &context);
    if (mrna != NULL) {
      mrna->idx.deleteme = TRUE;
    }
  }
  
  /* remove gene if requested */
  if (c->remove_gene) {
    if (gene == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &context);
    }
    if (gene != NULL) {
      gene->idx.deleteme = TRUE;
    }
  }

  /* get Prot-ref before removing transcript ID/marking protein Bioseq for deletion */
  protBsp = BioseqFindFromSeqLoc (sfp->product);
  protSfp = SeqMgrGetBestProteinFeature (protBsp, &context);
  if (protSfp != NULL) {
    prp = protSfp->data.value.ptrvalue;
  }

  if (c->featdef_to == FEATDEF_misc_feature) {
    if (prp != NULL) {
      if (prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue)) {
        prot_name = prp->name->data.ptrvalue;
      }
      if (!StringHasNoText (prp->desc)) {
        prot_desc = prp->desc;
      }
    }

    if (StringHasNoText (sfp->comment)) {
      comment_len = 1;
    } else {
      comment_len = StringLen (sfp->comment) + 3;
    }
    
    if (prot_name != NULL) {
      if (comment_len > 1) {
        comment_len += 2;
      }
      comment_len += StringLen (c->note_fmt) + StringLen (prot_name);
    }
    if (prot_desc != NULL) {
      if (comment_len > 1) {
        comment_len += 2;
      }
      comment_len += StringLen (prot_desc);
    }

    if (comment_len > 1) {
      new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
      new_comment[0] = 0;

      if (prot_name != NULL) {
        sprintf (new_comment, c->note_fmt, prot_name);
      }
      if (!StringHasNoText (sfp->comment)) {
        if (new_comment[0] != 0) {
          StringCat (new_comment, "; ");
        }
        StringCat (new_comment, sfp->comment);
      }
      if (prot_desc != NULL) {
        if (new_comment[0] != 0) {
          StringCat (new_comment, "; ");
        }
        StringCat (new_comment, prot_desc);
      }
    }

    if (c->keep_cds) {
      new_sfp = CreateNewFeatureOnBioseq (BioseqFindFromSeqLoc (sfp->location), SEQFEAT_IMP, sfp->location);
      if (c->remove_transcript_id) {
        if (protBsp != NULL) {
          protBsp->idx.deleteme = TRUE;
        }
        sfp->product = SeqLocFree (sfp->product);
      }
    } else {
      new_sfp = sfp;
      new_sfp->data.value.ptrvalue = CdRegionFree (new_sfp->data.value.ptrvalue);
      new_sfp->data.choice = SEQFEAT_IMP;
      /* if we don't keep the CDS, we need to remove the product */
      if (protBsp != NULL) {
        protBsp->idx.deleteme = TRUE;
      }
      if (! IsFeatInGPS (sfp) || c->remove_transcript_id) {
        sfp->product = SeqLocFree (sfp->product);
      }     
    }
    ifp = ImpFeatNew ();
    new_sfp->data.value.ptrvalue = ifp;
    ifp->key = StringSave ("misc_feature");

    new_sfp->comment = MemFree (new_sfp->comment);
    new_sfp->comment = new_comment;

    /* set the subtype to zero so that it will be reindexed */
    sfp->idx.subtype = 0;
  } else {
    /* conversion to pseudogene */

    /* find the gene associated with the coding region */
    if (gene == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &context);
    }

    if (c->keep_cds) {
      if (protSfp != NULL) {
        /* append protein comment to CDS comment */
        if (!StringHasNoText (protSfp->comment)) {
          SetStringValue (&(sfp->comment), protSfp->comment, ExistingTextOption_append_semi);
        }

        if (prp != NULL) {
          /* prepend protein description (if any) to CDS comment */
          if (!StringHasNoText (prp->desc)) {
            SetStringValue (&(sfp->comment), prp->desc, ExistingTextOption_prefix_semi);
          }

          /* Get the protein names and prepend them to the coding region comment or use them to create a gene if one didn't exist */
          if (gene == NULL) {
            grp = GeneRefNew ();
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              SetStringValue (&(grp->desc), vnp->data.ptrvalue, ExistingTextOption_append_semi);
            }
            gene = CreateNewFeatureOnBioseq (BioseqFindFromSeqLoc (sfp->location), SEQFEAT_GENE, sfp->location);
            gene->data.value.ptrvalue = (Pointer) grp;
          } else {
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              SetStringValue (&(sfp->comment), vnp->data.ptrvalue, ExistingTextOption_prefix_semi);
            }
          }
        }
      }
      if (protBsp != NULL) {
        /* Delete the protein Bioseq that */
        /* the CDS points to.             */
        protBsp->idx.deleteme = TRUE;
      }
      /* Mark the CDS as Pseudo */
      sfp->pseudo = TRUE;
    
      if (! IsFeatInGPS (sfp) || c->remove_transcript_id) {
        sfp->product = SeqLocFree (sfp->product);
      }     

      /* Set the subtype to zero so that it will be reindexed */
      sfp->idx.subtype = 0;

      /* Mark the overlapping gene as pseudo */
      if (gene != NULL) {
        gene->pseudo = TRUE;
      }
    } else {
      /* remove CDS */
      if (gene == NULL) {
        /* need to create gene */
        gene = CreateNewFeatureOnBioseq (BioseqFindFromSeqLoc (sfp->location), SEQFEAT_GENE, sfp->location);
      }
      /* if gene doesn't have generef, create one */
      grp = gene->data.value.ptrvalue;
      if (grp == NULL) {
        grp = GeneRefNew ();
        gene->data.value.ptrvalue = grp;
      }

      /* mark the gene as pseudo */
      gene->pseudo = TRUE;

      /* prepend CDS comment to gene comment */
      if (!StringHasNoText (sfp->comment)) {
        SetStringValue (&(gene->comment), sfp->comment, ExistingTextOption_prefix_semi);
      }

      if (protSfp != NULL) {
        /* prepend protein comment to gene comment */
        if (!StringHasNoText (protSfp->comment)) {
          SetStringValue (&(gene->comment), protSfp->comment, ExistingTextOption_prefix_semi);
        }
      }

      if (prp != NULL) {
        /* prepend protein description to gene comment */
        if (!StringHasNoText (prp->desc)) {
          SetStringValue (&(gene->comment), prp->desc, ExistingTextOption_prefix_semi);
        }
        if (StringHasNoText (grp->locus)) {
          /* use protein product names for gene locus if gene doesn't already have one */
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            SetStringValue (&(grp->locus), vnp->data.ptrvalue, ExistingTextOption_append_semi);
          }
        } else if (StringHasNoText (grp->desc)) {
          /* use protein product names for gene description if gene doesn't already have one */
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            SetStringValue (&(grp->desc), vnp->data.ptrvalue, ExistingTextOption_append_semi);
          }
        } else {
          /* prepend list of protein product names to gene comment */
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            SetStringValue (&prot_list, vnp->data.ptrvalue, ExistingTextOption_append_semi);
          }
          SetStringValue (&(gene->comment), prot_list, ExistingTextOption_prefix_semi);
          prot_list = MemFree (prot_list);
        }
      }      

      /* Delete the protein Bioseq and */
      /* the CDS that points to it.    */
      sfp->idx.deleteme = TRUE;
      if (protBsp != NULL) {
        protBsp->idx.deleteme = TRUE;
      }
    }
  }
}  


typedef struct convertcdsoptionsdlg {
  DIALOG_MESSAGE_BLOCK
  
  GrouP dest_feat_type;
  ButtoN keep_cds;
  ButtoN remove_mrna;
  ButtoN remove_gene;
  ButtoN remove_transcript_id;
  GrouP  all_or_conditional;
  GrouP  conditions;
  ButtoN require_internal_stop;
  ButtoN require_bad_start;
  ButtoN require_bad_stop;
  ButtoN require_pseudo;
  GrouP  require_or;
  GrouP note_fmt;
  DialoG string_constraint;
  DialoG id_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;  
} ConvertCDSOptionsDlgData, PNTR ConvertCDSOptionsDlgPtr;


static void DataToConvertCDSOptionsDlg (DialoG d, Pointer data)
{
  ConvertCDSOptionsDlgPtr dlg;
  ConvertCDSOptionsPtr    c;

  dlg = (ConvertCDSOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  c = (ConvertCDSOptionsPtr) data;

  if (c == NULL) {
    c = ConvertCDSOptionsNew ();
    DataToConvertCDSOptionsDlg (d, c);
    c = ConvertCDSOptionsFree (c);
  } else {
    if (c->featdef_to == FEATDEF_misc_feature) {      
      SetValue (dlg->dest_feat_type, 1);
    } else {
      SetValue (dlg->dest_feat_type, 2);
    }
    SetStatus (dlg->keep_cds, c->keep_cds);
    SetStatus (dlg->remove_mrna, c->remove_mrna);
    SetStatus (dlg->remove_gene, c->remove_gene);
    SetStatus (dlg->remove_transcript_id, c->remove_transcript_id);

    if (!c->require_internal_stop && !c->require_bad_start && !c->require_bad_stop && !c->require_pseudo) {
      SetValue (dlg->all_or_conditional, 1);
      Disable (dlg->conditions);
    } else {
      SetValue (dlg->all_or_conditional, 2);
      Enable (dlg->conditions);
    }
    if (c->require_or) {
      SetValue (dlg->require_or, 1);
    } else {
      SetValue (dlg->require_or, 2);
    }

    SetStatus (dlg->require_internal_stop, c->require_internal_stop);
    SetStatus (dlg->require_bad_start, c->require_bad_start);
    SetStatus (dlg->require_bad_stop, c->require_bad_stop);
    SetStatus (dlg->require_pseudo, c->require_pseudo);

    if (StringCmp (c->note_fmt, kSimilarNoteFmt) == 0) {
      SetValue (dlg->note_fmt, 1);
    } else {
      SetValue (dlg->note_fmt, 2);
    }
    PointerToDialog (dlg->string_constraint, c->string_constraint);
    PointerToDialog (dlg->id_constraint, c->id_constraint);
  }
  if (dlg->change_notify) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer DataFromConvertCDSOptionsDlg (DialoG d)
{
  ConvertCDSOptionsDlgPtr dlg;
  ConvertCDSOptionsPtr    c;

  dlg = (ConvertCDSOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  c = ConvertCDSOptionsNew();

  if (GetValue (dlg->dest_feat_type) == 1) {
    c->featdef_to = FEATDEF_misc_feature;
  } else {
    c->featdef_to = FEATDEF_GENE;
  }
  c->keep_cds = GetStatus (dlg->keep_cds);
  c->remove_mrna = GetStatus (dlg->remove_mrna);
  c->remove_gene = GetStatus (dlg->remove_gene);
  c->remove_transcript_id = GetStatus (dlg->remove_transcript_id);

  if (GetValue (dlg->all_or_conditional) == 1) {
    c->require_internal_stop = FALSE;
    c->require_bad_start = FALSE;
    c->require_bad_stop = FALSE;
    c->require_pseudo = FALSE;
    c->require_or = TRUE;
  } else {
    c->require_internal_stop = GetStatus (dlg->require_internal_stop);
    c->require_bad_start = GetStatus (dlg->require_bad_start);
    c->require_bad_stop = GetStatus (dlg->require_bad_stop);
    c->require_pseudo = GetStatus (dlg->require_pseudo);
    if (GetValue (dlg->require_or) == 1) {
      c->require_or = TRUE;
    } else {
      c->require_or = FALSE;
    }
  }

  if (GetValue (dlg->note_fmt) == 1) {
    c->note_fmt = MemFree (c->note_fmt);
    c->note_fmt = StringSave (kSimilarNoteFmt);
  } else {
    c->note_fmt = MemFree (c->note_fmt);
    c->note_fmt = StringSave (kViralNoteFmt);
  }

  c->string_constraint = DialogToPointer (dlg->string_constraint);
  c->id_constraint = DialogToPointer (dlg->id_constraint);
  return c;
}


static ValNodePtr TestConvertCDSOptionsDlg (DialoG d)
{
  ConvertCDSOptionsDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (ConvertCDSOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  if (GetValue (dlg->all_or_conditional) == 2) {
    if (!GetStatus (dlg->require_internal_stop)
        && !GetStatus (dlg->require_bad_start)
        && !GetStatus (dlg->require_bad_stop)
        && !GetStatus (dlg->require_pseudo)) {
      ValNodeAddPointer (&err_list, 0, "no condition");
    }
  }
  return err_list;
}


static void ChangeAllOrConditional (GrouP g)
{
  ConvertCDSOptionsDlgPtr dlg;

  dlg = (ConvertCDSOptionsDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) {
    return;
  }

  if (GetValue (dlg->all_or_conditional) == 1) {
    Disable (dlg->conditions);
  } else {
    Enable (dlg->conditions);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeCDSOptionsButton (ButtoN b)
{
  ConvertCDSOptionsDlgPtr dlg;

  dlg = (ConvertCDSOptionsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG ConvertCDSOptionsDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ConvertCDSOptionsDlgPtr dlg;
  GrouP                   p, g, g2, opts;
  
  dlg = (ConvertCDSOptionsDlgPtr) MemNew (sizeof (ConvertCDSOptionsDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToConvertCDSOptionsDlg;
  dlg->fromdialog = DataFromConvertCDSOptionsDlg;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestConvertCDSOptionsDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->dest_feat_type = NormalGroup (p, 2, 0, "Convert to", NULL, NULL);
  RadioButton (dlg->dest_feat_type, "misc_feat");
  RadioButton (dlg->dest_feat_type, "pseudogene");
  SetValue (dlg->dest_feat_type, 1);

  opts = HiddenGroup (p, 0, 4, NULL);
  dlg->keep_cds = CheckBox (opts, "Keep original CDS", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->keep_cds, dlg, NULL);
  dlg->remove_mrna = CheckBox (opts, "Remove overlapping mRNA", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->remove_mrna, dlg, NULL);
  dlg->remove_gene = CheckBox (opts, "Remove overlapping gene", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->remove_gene, dlg, NULL);
  dlg->remove_transcript_id = CheckBox (opts, "Remove transcript ID", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->remove_transcript_id, dlg, NULL);

  g = NormalGroup (p, 0, 2, "Conditions", NULL, NULL);
  dlg->all_or_conditional = HiddenGroup (g, 0, 2, ChangeAllOrConditional);
  RadioButton (dlg->all_or_conditional, "All");
  RadioButton (dlg->all_or_conditional, "Where");
  SetObjectExtra (dlg->all_or_conditional, dlg, NULL);
  SetValue (dlg->all_or_conditional, 1);
  g2 = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (g2, "", 10, dialogTextHeight, systemFont, 'l');
  dlg->conditions = HiddenGroup (g2, 0, 5, NULL);
  dlg->require_internal_stop = CheckBox (dlg->conditions, "CDS has internal stop codon", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->require_internal_stop, dlg, NULL);
  dlg->require_bad_start = CheckBox (dlg->conditions, "CDS has bad start codon", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->require_bad_start, dlg, NULL);
  dlg->require_bad_stop = CheckBox (dlg->conditions, "CDS has bad stop codon", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->require_bad_stop, dlg, NULL);
  dlg->require_pseudo = CheckBox (dlg->conditions, "CDS is pseudo", ChangeCDSOptionsButton);
  SetObjectExtra (dlg->require_pseudo, dlg, NULL);
  dlg->require_or = HiddenGroup (dlg->conditions, 2, 0, NULL);
  RadioButton (dlg->require_or, "OR");
  RadioButton (dlg->require_or, "AND");
  SetValue (dlg->require_or, 1);
  Disable (dlg->conditions);

  dlg->note_fmt = NormalGroup (p, 0, 2, "Comment Format", NULL, NULL);
  RadioButton (dlg->note_fmt, "similar to X");
  RadioButton (dlg->note_fmt, "nonfunctional X due to mutation");
  SetValue (dlg->note_fmt, 1);

  dlg->string_constraint = StringConstraintDialog (p, "Where CDS feature", FALSE, change_notify, change_userdata);
  dlg->id_constraint = StringConstraintDialog (p, "Where sequence ID", FALSE, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->dest_feat_type,
                              (HANDLE) opts, 
                              (HANDLE) g,
                              (HANDLE) dlg->note_fmt,
                              (HANDLE) dlg->string_constraint,
                              (HANDLE) dlg->id_constraint,
                              NULL);

  return (DialoG) p;
}


typedef struct convertcdsfrm {
  FORM_MESSAGE_BLOCK
  DialoG options;
  ButtoN accept_btn;
  ButtoN leave_dlg_up;
} ConvertCDSFrmData, PNTR ConvertCDSFrmPtr;


static void EnableCDSConversionAccept (Pointer data)
{
  ConvertCDSFrmPtr frm;
  ValNodePtr       err_list = NULL;

  frm = (ConvertCDSFrmPtr) data;
  if (frm == NULL) {
    return;
  }
  err_list = TestDialog (frm->options);
  if (err_list == NULL) {
    Enable (frm->accept_btn);
  } else {
    Disable (frm->accept_btn);
  }
  err_list = ValNodeFree (err_list);
}


static void DoCDSConversion (ButtoN b)
{
  ConvertCDSFrmPtr frm;
  SeqEntryPtr      sep;
  ConvertCDSOptionsPtr c;

  frm = (ConvertCDSFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  c = DialogToPointer (frm->options);

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  VisitFeaturesInSep (sep, c, ConvertOneCDS);

  c = ConvertCDSOptionsFree (c);

  DeleteMarkedObjects (frm->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);

  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  if (!GetStatus (frm->leave_dlg_up)) {
    Remove (frm->form);
  }
  Update ();
}


static void ConvertCDSMenuItem (IteM i)
{
  BaseFormPtr  bfp;
  ConvertCDSFrmPtr frm;
  WindoW           w;
  GrouP            h, c;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
    
  frm = (ConvertCDSFrmPtr) MemNew (sizeof (ConvertCDSFrmData));
  if (frm == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Convert CDS", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  frm->options = ConvertCDSOptionsDialog (h, EnableCDSConversionAccept, frm);

  c = HiddenGroup (h, 3, 0, NULL);
  frm->accept_btn = DefaultButton (c, "Accept", DoCDSConversion);
  SetObjectExtra (frm->accept_btn, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  frm->leave_dlg_up = CheckBox (c, "Leave dialog up", NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->options,
                              (HANDLE) c,
                              NULL);
  EnableCDSConversionAccept (frm);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


static void NewConvertCDSWithInternalStopsToMiscFeat (IteM i)
{
  BaseFormPtr  bfp;
  ConvertCDSOptionsPtr o;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  o = ConvertCDSOptionsNew ();
  o->require_internal_stop = TRUE;
  o->note_fmt = MemFree (o->note_fmt);
  o->note_fmt = StringSave (kSimilarNoteFmt);
  VisitFeaturesInSep (sep, o, ConvertOneCDS);
  o = ConvertCDSOptionsFree (o);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void NewConvertCDSWithInternalStopsToMiscFeatViral (IteM i)
{
  BaseFormPtr  bfp;
  ConvertCDSOptionsPtr o;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  o = ConvertCDSOptionsNew ();
  o->require_internal_stop = TRUE;
  o->note_fmt = MemFree (o->note_fmt);
  o->note_fmt = StringSave (kViralNoteFmt);
  VisitFeaturesInSep (sep, o, ConvertOneCDS);
  o = ConvertCDSOptionsFree (o);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}


static void RemoveTaxonXrefsFromList (ValNodePtr PNTR list)
{
  ValNodePtr vnp, next;
  DbtagPtr   dbt;
  
  if (list == NULL) return;
  vnp = *list;
  while (vnp != NULL)
  {
  	next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringICmp ((CharPtr) dbt->db, "taxon") == 0) {
      *list = vnp->next;
      vnp->next = NULL;
      DbtagFree (dbt);
      ValNodeFree (vnp);
    } else {
      list = (ValNodePtr PNTR) &(vnp->next);
    }    
    vnp = next;
  }
}

static Boolean RemoveTaxonProc (GatherObjectPtr gop)

{
  BioSourcePtr     biop = NULL;
  OrgRefPtr        orp;
  ObjValNodePtr    ovn;
  SeqDescPtr       sdp;
  SeqFeatPtr       sfp;
  ValNodePtr       vnp;
  BoolPtr          remove_from_biosources;

  switch (gop->itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gop->dataptr;
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_BIOSRC) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      }
      RemoveTaxonXrefsFromList (&(sfp->dbxref));
      break;
    case OBJ_SEQDESC :
      sdp = (SeqDescPtr) gop->dataptr;
      ovn = (ObjValNodePtr) sdp;
      if (sdp != NULL && sdp->extended != 0 && ovn->idx.subtype == Seq_descr_source ) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
      break;
    default :
      return TRUE;
  }
  remove_from_biosources = (BoolPtr) (gop->userdata);
  if (!*remove_from_biosources) return TRUE;
  if (biop == NULL) return TRUE;
  orp = biop->org;
  if (orp == NULL) return TRUE;
  vnp = orp->db;
  if (vnp == NULL) return TRUE;
  RemoveTaxonXrefsFromList (&(orp->db));

  return TRUE;
}

static void RemoveTaxonXrefs (IteM i, Boolean remove_from_biosources)

{
  BaseFormPtr  bfp;
  Boolean      objMgrFilter [OBJ_MAX];
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
  objMgrFilter [OBJ_SEQFEAT] = TRUE;
  objMgrFilter [OBJ_SEQDESC] = TRUE;
  GatherObjectsInEntity (bfp->input_entityID, 0, NULL,
                         RemoveTaxonProc, (Pointer) &remove_from_biosources, objMgrFilter);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void RemoveTaxonXrefsFromFeatures (IteM i)
{
  RemoveTaxonXrefs (i, FALSE);
}

static void RemoveTaxonXrefsFromFeaturesAndBioSources (IteM i)
{
  RemoveTaxonXrefs (i, TRUE);
}

/*=========================================================================*/
/*                                                                         */
/* SuppressError_Callback () --                                            */
/*                                                                         */
/*=========================================================================*/

static void SuppressError_Callback (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr        grp = NULL;
  SeqFeatXrefPtr    xref;
  SeqFeatPtr        gene;
  SeqMgrFeatContext geneContext;

  /* If we already have a gene xref, then we're OK */

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
    return;

  /* If we have a completely overlapping gene, we're OK */

  gene = SeqMgrGetOverlappingGene (sfp->location, &geneContext);
  if (gene != NULL)
    return;

  /* If we have a partially overlapping gene */
  /* then it is a problem.                   */

  gene = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL,
				      0, NULL, SIMPLE_OVERLAP, &geneContext);
  if (gene == NULL)
    return;

  /* If we got to here, then we have no gene xref and */
  /* and a partially overlapping gene, so add on a    */
  /* blank gene xref to suppress the error message.   */

  grp = GeneRefNew ();
  if (grp != NULL) {
    xref = SeqFeatXrefNew ();
    xref->data.choice = SEQFEAT_GENE;
    xref->data.value.ptrvalue = grp;
    xref->next = sfp->xref;
    sfp->xref = xref;
  }

  /* Return successfully */

  return;
}

/*=========================================================================*/
/*                                                                         */
/* SuppressCDSGeneRangeError () --                                         */
/*                                                                         */
/*=========================================================================*/

static void SuppressCDSGeneRangeError (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  Uint2        entityID;
  SeqEntryPtr  sep;

  /* Get the Bioseq */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  bsp =  GetBioseqGivenIDs (bfp->input_entityID, 1, OBJ_BIOSEQ);
  sep = SeqMgrGetSeqEntryForData (bsp);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  sep = GetTopSeqEntryForEntityID (entityID);

  /* Explore the Bioseq's CDS features, searching for */
  /* ones that would generate a CDSGeneRange error.   */

  VisitFeaturesInSep (sep, NULL, SuppressError_Callback);

  /* Force an update and redraw */

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  ArrowCursor ();
  Update ();

  /* Return successfully */

  return;
}

/*=========================================================================*/
/*                                                                         */
/* RestoreError_Callback () --                                             */
/*                                                                         */
/*=========================================================================*/

static void RestoreError_Callback (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr        grp = NULL;
  SeqFeatPtr        gene;
  SeqMgrFeatContext geneContext;

  /* Check for a suppressed gene xref */

  grp = SeqMgrGetGeneXref (sfp);
  if (NULL == grp)
    return;

  if (FALSE == SeqMgrGeneIsSuppressed (grp))
    return;

  /* Check to make sure we don't have a fully overlapping gene */

  gene = SeqMgrGetOverlappingGene (sfp->location, &geneContext);
  if (NULL != gene)
    return;

  /* If we have a partially overlapping gene */
  /* then we need to remove the suppression. */

  gene = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL,
				      0, NULL, SIMPLE_OVERLAP, &geneContext);
  if (gene == NULL)
    return;

  /* If we got to here, then we have a suppressed xref */
  /* and a partially overlapping gene, so remove the   */
  /* suppression.                                      */

  GeneRefFree (grp);
  sfp->xref = SeqFeatXrefFree (sfp->xref);

  /* Return successfully */

  return;
}

/*=========================================================================*/
/*                                                                         */
/* RestoreCDSGeneRangeError () --                                          */
/*                                                                         */
/*=========================================================================*/

static void RestoreCDSGeneRangeError (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  Uint2        entityID;
  SeqEntryPtr  sep;

  /* Get the Bioseq */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  bsp =  GetBioseqGivenIDs (bfp->input_entityID, 1, OBJ_BIOSEQ);
  sep = SeqMgrGetSeqEntryForData (bsp);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  sep = GetTopSeqEntryForEntityID (entityID);

  /* Explore the Bioseq's CDS features, searching for */
  /* ones that would generate a CDSGeneRange error.   */

  VisitFeaturesInSep (sep, NULL, RestoreError_Callback);

  /* Force an update and redraw */

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  ArrowCursor ();
  Update ();

  /* Return successfully */

  return;
}

/*-------------------------------------------------------------------------*/
/*                                                                         */
/* SuppressFeatureGeneXref () -- Suppresses any gene xref on the feature   */
/*                               that is passed to it.                     */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static void SuppressFeatureGeneXref (SeqFeatPtr sfp)
{
  GeneRefPtr      grp = NULL;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL)
  {
    return;
  }
  /* If there is a gene xref, then change it */
  /* to a suppression gene xref.             */

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
  {
    if (SeqMgrGeneIsSuppressed (grp) == FALSE)
	  {
	    if (NULL != grp->locus)
	    {
	      MemFree(grp->locus);
        grp->locus = NULL;
	    }
	    if (NULL != grp->allele)
	    {
	      MemFree(grp->allele);
	      grp->allele = NULL;
	    }
	    if (NULL != grp->desc)
	    {
	      MemFree(grp->desc);
	      grp->desc = NULL;
	    }
	    if (NULL != grp->maploc)
      {
	      MemFree (grp->maploc);
	      grp->maploc = NULL;
      }
      if (NULL != grp->locus_tag)
      {
        MemFree (grp->locus_tag);
        grp->locus_tag = NULL;
      }
      grp->db  = ValNodeFreeData (grp->db);
      grp->syn = ValNodeFreeData (grp->syn);
	  }
	}    
  /* Otherwise, if there is an overlapping gene, add */
  /* a suppression xref for it.                      */
  else if (SeqMgrGetOverlappingGene (sfp->location, NULL) != NULL)
  {
    grp = GeneRefNew ();
    if (grp != NULL) 
    {
	    xref = SeqFeatXrefNew ();
	    xref->data.choice = SEQFEAT_GENE;
	    xref->data.value.ptrvalue = grp;
	    xref->next = sfp->xref;
	    sfp->xref = xref;
    }
  }
}

/*=============================================================================*/
/*                                                                             */
/* SuppressGenesOnFeatures () -- Suppress gene xrefs on selected feature types */
/*                                                                             */
/*=============================================================================*/

typedef struct suppressgenes
{
  FORM_MESSAGE_BLOCK
  DialoG feature_choice_dlg;
  DialoG constraint_dlg;
  ButtoN accept_btn;
  
} SuppressGenesData, PNTR SuppressGenesPtr;

static void SetSuppressGenesAccept (Pointer userdata)
{
  SuppressGenesPtr dlg;
  ValNodePtr       feature_type_list;

  dlg = (SuppressGenesPtr) userdata;
  if (dlg == NULL)
  {
    return;
  }
  feature_type_list = DialogToPointer (dlg->feature_choice_dlg);
  if (feature_type_list == NULL)
  {
    Disable (dlg->accept_btn);
  }
  else
  {
    Enable (dlg->accept_btn);
  }
  feature_type_list = ValNodeFree (feature_type_list);
}

static void SuppressOneGeneOnFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  SuppressFeatureGeneXref (sfp);  
}

static void DoSuppressGenesOnFeatures (ButtoN b)
{
  SuppressGenesPtr dlg;
  ValNodePtr       feature_type_list, vnp;
  FilterSetPtr     fsp;
  Uint1            feat_def_choice;
  SeqEntryPtr      sep;

  dlg = (SuppressGenesPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  if (sep == NULL)
  {
    return;
  }
  
  feature_type_list = DialogToPointer (dlg->feature_choice_dlg);
  if (feature_type_list == NULL)
  {
    return;
  }
  
  fsp = (FilterSetPtr) DialogToPointer (dlg->constraint_dlg);
  
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    feat_def_choice = vnp->choice;
    if (feat_def_choice == 255)
    {
      feat_def_choice = 0;
    }
    OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                         SuppressOneGeneOnFeature,
                                         NULL, 0, feat_def_choice, 0, dlg);
  }
  
  ValNodeFree (feature_type_list);
  FilterSetFree (fsp);
  
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  Remove (dlg->form);
  Update ();
}

static void SuppressGenesOnFeatures (IteM i)
{
  BaseFormPtr      bfp;
  SuppressGenesPtr dlg;
  WindoW           w;
  GrouP            h, c;
  ButtoN           b;
  SeqEntryPtr      sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (SuppressGenesPtr) MemNew (sizeof (SuppressGenesData));
  if (dlg == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Suppress Genes on Features", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupExtraProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  dlg->feature_choice_dlg = FeatureSelectionDialogEx (h, TRUE, sep, SetSuppressGenesAccept, dlg);
  dlg->constraint_dlg = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text");
  
  c = HiddenGroup (h, 2, 0, NULL);
  dlg->accept_btn = PushButton (c, "Accept", DoSuppressGenesOnFeatures);
  SetObjectExtra (dlg->accept_btn, dlg, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_choice_dlg,
                              (HANDLE) dlg->constraint_dlg,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


/*=======================================================================*/
/*                                                                       */
/*  CopyGeneRef () -                                                     */
/*                                                                       */
/*=======================================================================*/

static GeneRefPtr CopyGeneRef (GeneRefPtr srcGrp)
{
  GeneRefPtr   destGrp;
  DbtagPtr     destDbt;
  DbtagPtr     srcDbt;
  ValNodePtr   vnp;

  destGrp = GeneRefNew ();

  /* Copy the string fields from source to destination */
  
  if (srcGrp->locus != NULL)
    destGrp->locus = StringSave (srcGrp->locus);

  if (srcGrp->allele != NULL)
    destGrp->allele = StringSave (srcGrp->allele);

  if (srcGrp->desc != NULL)
    destGrp->desc = StringSave (srcGrp->desc);

  if (srcGrp->maploc != NULL)
    destGrp->maploc = StringSave (srcGrp->maploc);

  if (srcGrp->locus_tag != NULL)
    destGrp->locus_tag = StringSave (srcGrp->locus_tag);

  /* Copy the DB references */

  destGrp->db = NULL;

  if (srcGrp->db != NULL)
    {
      vnp = srcGrp->db;
      while (vnp != NULL)
	{
	  ValNodeAdd (&(destGrp->db));
	  destDbt = DbtagNew ();
	  srcDbt  = (DbtagPtr) vnp->data.ptrvalue;
	  if (srcDbt != NULL)
	    {
	      if (srcDbt->db != NULL)
		{
		  destDbt->db = (CharPtr) MemNew (sizeof(srcDbt->db));
		  StringCpy (destDbt->db, srcDbt->db);
		}
	      if (srcDbt->tag != NULL)
		{
		  destDbt->tag = ObjectIdNew ();
		  destDbt->tag->id = srcDbt->tag->id;
		  destDbt->tag->str = (CharPtr)MemNew(sizeof(srcDbt->tag->str));
		  StringCpy (destDbt->tag->str, srcDbt->tag->str);
		}
	    }
	  vnp = vnp->next;
	}
    }

  /* Copy the synonyms */

  if (srcGrp->syn != NULL)
    {
      vnp = srcGrp->syn;
      while (vnp != NULL)
	{
	  ValNodeCopyStr (&(destGrp->syn), vnp->choice,
			  (CharPtr) vnp->data.ptrvalue);
	  vnp = vnp->next;
	}
    }

  /* Return the new gene reference */

  return destGrp;

}

static void ZapRBSGene (SeqFeatPtr sfp, Pointer userdata)

{
  Int4               diff;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;

  if (sfp->idx.subtype != FEATDEF_RBS) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) return;
  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;
  diff = SeqLocAinB (sfp->location, gene->location);
  if (diff == 0) {
    gene->idx.deleteme = TRUE;
  } else if (diff > 0) {
    SeqLocSubtract (gene->location, sfp->location);
  }
}

static void RemoveRBSGenes (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  WatchCursor ();
  Update ();

  VisitFeaturesInSep (sep, NULL, ZapRBSGene);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

/*=========================================================================*/
/*                                                                         */
/* MapRBSsToDownstreamGene () -- Creates a xref to the correct gene for    */
/*                               all RBS features.                         */
/*                                                                         */
/*=========================================================================*/

static void MapRBSsToDownstreamGene (IteM i)
{
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  Uint2              entityID;
  Boolean            found;
  SeqMgrFeatContext  cdsContext;
  Int4               cdsLeft;
  Int4               cdsRight;
  SeqFeatPtr         cdsSfp = NULL;
  SeqLocPtr          cdsSlp;
  GeneRefPtr         geneGrp;
  SeqFeatPtr         geneSfp = NULL;
  Boolean            overshot = FALSE;
  SeqFeatPtr         prevCdsSfp = NULL;
  SeqFeatXrefPtr     prevXref = NULL;
  SeqMgrFeatContext  rbsContext, geneContext;
  GeneRefPtr         rbsGrp;
  Int4               rbsLeft;
  Int4               rbsRight;
  SeqFeatPtr         rbsSfp = NULL;
  SeqLocPtr          rbsSlp;
  Uint1              strand;
  SeqFeatXrefPtr     xref = NULL;

  /* Get the Bioseq */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  WatchCursor ();
  Update ();

  bsp =  GetBioseqGivenIDs (bfp->input_entityID,
			   bfp->input_itemID,
			   bfp->input_itemtype);

  rbsSfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_RBS, &rbsContext);

  /* Loop through all RBS features */

  while (NULL != rbsSfp)
  {
    /* If the feature has an existing gene xref */
    /* then we need to delete it.               */
  
    xref = rbsSfp->xref;
    while (xref != NULL && xref->data.choice != SEQFEAT_GENE)
	  {
	    prevXref = xref;
	    xref = xref->next;
	  }
      
    if (xref != NULL)
	  {
	    if (prevXref != NULL)
	      prevXref->next = xref->next;
	    else
	      rbsSfp->xref = xref->next;
	    xref->next = NULL;
	    SeqFeatXrefFree (xref);
	    xref = NULL;
	  }

    /* Search for a 'downstream' gene xref */
  
    if (overshot)
	    overshot = FALSE;
    else
    {
	    prevCdsSfp = cdsSfp;
	    cdsSfp = SeqMgrGetNextFeature (bsp, cdsSfp, 0, FEATDEF_CDS,
					                           &cdsContext);
    }
      
    rbsSlp  = rbsSfp->location;
    strand  = SeqLocStrand (rbsSlp);

    if (strand == Seq_strand_plus)
    {
	    /* If the gene starts at or 'after' the  */
	    /* RBS, then the RBS probably belongs to */
	    /* it, so make that gene a xref for the  */
	    /* RBS feature.                          */
	  
	    found = FALSE;
	    while ((cdsSfp != NULL) && (!found))
	    {
	      cdsSlp  = cdsSfp->location;
	      cdsLeft = GetOffsetInBioseq (cdsSlp, bsp,
					   SEQLOC_LEFT_END);
	      rbsLeft = GetOffsetInBioseq (rbsSlp , bsp,
					   SEQLOC_LEFT_END);
	      if (rbsLeft <= cdsLeft)
        {
          found = TRUE;
          geneSfp = SeqMgrGetOverlappingGene (cdsSlp, &geneContext);
          if (geneSfp != NULL) {
            geneGrp = (GeneRefPtr) geneSfp->data.value.ptrvalue;
            rbsGrp = CopyGeneRef (geneGrp);
            if (rbsGrp != NULL)
            {
		          xref = SeqFeatXrefNew ();
		          xref->data.choice = SEQFEAT_GENE;
		          xref->data.value.ptrvalue = rbsGrp;
		          xref->next = rbsSfp->xref;
		          rbsSfp->xref = xref;
		        }
		        
		        /* enlarge gene to cover RBS */
		        if (geneContext.left > rbsContext.left)
		        {
              ExtendSeqLocToPosition (geneSfp->location, TRUE, rbsContext.left);
		        }
		        
          }
        }
	      else
        {
          if (overshot)
            overshot = FALSE;
          else
          {
            prevCdsSfp = cdsSfp;
            cdsSfp = SeqMgrGetNextFeature (bsp, cdsSfp, 0,
                                           FEATDEF_CDS,
                                           &cdsContext);
          }
        }
	    }
    }
    else if (strand == Seq_strand_minus)
    {
	    found = FALSE;
	    while ((cdsSfp != NULL) && (!found))
      {
        cdsSlp   = cdsSfp->location;
        cdsRight = GetOffsetInBioseq (cdsSlp, bsp,
                                      SEQLOC_RIGHT_END);
        rbsRight = GetOffsetInBioseq (rbsSlp , bsp,
                                      SEQLOC_RIGHT_END);
        if (cdsRight > rbsRight)
        {
          found = TRUE;
          if (prevCdsSfp != NULL)
          {
            overshot = TRUE;
            geneSfp = SeqMgrGetOverlappingGene (prevCdsSfp->location,
                                                NULL);
            if (geneSfp != NULL) {
              geneGrp = (GeneRefPtr) geneSfp->data.value.ptrvalue;
              rbsGrp = CopyGeneRef (geneGrp);
              if (rbsGrp != NULL)
              {
                xref = SeqFeatXrefNew ();
                xref->data.choice = SEQFEAT_GENE;
                xref->data.value.ptrvalue = rbsGrp;
                xref->next = rbsSfp->xref;
                rbsSfp->xref = xref;
              }
              
		          /* enlarge gene to cover RBS */
		          if (geneContext.right < rbsContext.right)
		          {
                ExtendSeqLocToPosition (geneSfp->location, TRUE, rbsContext.right);
		          }
              
            }
			
          }
        }
        else
        {
          if (overshot)
            overshot = FALSE;
          else
          {
            prevCdsSfp = cdsSfp;
            cdsSfp = SeqMgrGetNextFeature (bsp, cdsSfp, 0,
                                           FEATDEF_CDS,
                                           &cdsContext);
          }
        }
	      
      }
    }
      
    rbsSfp  = SeqMgrGetNextFeature (bsp, rbsSfp, 0, FEATDEF_RBS, &rbsContext);
  }

  /* Force an update and redraw */

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  ArrowCursor ();
  Update ();

  /* Return successfully */

  return;

}

static Boolean MarkBankitComments (GatherObjectPtr gop)

{
  ObjectIdPtr    oip;
  ObjValNodePtr  ovn;
  SeqDescrPtr    sdp;
  UserObjectPtr  uop;

  if (gop == NULL || gop->itemtype != OBJ_SEQDESC) return TRUE;
  sdp = (SeqDescrPtr) gop->dataptr;
  if (sdp == NULL || sdp->choice != Seq_descr_user || sdp->extended == 0) return TRUE;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return TRUE;
  oip = uop->type;
  if (oip == NULL || StringCmp (oip->str, "Submission") != 0) return TRUE;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
  return TRUE;
}

static void RemoveBankitComments (IteM i)

{
  BaseFormPtr  bfp;
  Boolean      objMgrFilter [OBJ_MAX];
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
  objMgrFilter [OBJ_SEQDESC] = TRUE;
  GatherObjectsInEntity (bfp->input_entityID, 0, NULL,
                         MarkBankitComments, NULL, objMgrFilter);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoCopyBankItComments (SeqDescrPtr sdp, Pointer userdata)

{
  ObjectIdPtr    oip;
  CharPtr        str;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;
  SeqDescrPtr    vnp;

  if (sdp == NULL || sdp->choice != Seq_descr_user || sdp->extended == 0) return;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringCmp (oip->str, "Submission") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
      str = ufp->data.ptrvalue;
      if (! StringHasNoText (str)) {
        vnp = SeqDescrNew (NULL);
        if (vnp != NULL) {
          vnp->choice = Seq_descr_comment;
          vnp->data.ptrvalue = StringSave (str);
          vnp->next = sdp->next;
          sdp->next = vnp;
        }
      }
    }
  }

}

static void CopyBankitComments (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitDescriptorsInSep (sep, NULL, DoCopyBankItComments);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static ENUM_ALIST(molinfo_biomol_alist)
  {"Any",                  254},
  {" ",                      0},
  {"Genomic DNA or RNA",     1},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Peptide",                8},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"cRNA",                  11},
  {"Transcribed RNA",       13},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_tech_alist)
  {"Any",               254},
  {" ",                   0},
  {"Standard",            1},
  {"EST",                 2},
  {"STS",                 3},
  {"Survey",              4},
  {"Genetic Map",         5},
  {"Physical Map",        6},
  {"Derived",             7},
  {"Concept-Trans",       8},
  {"Seq-Pept",            9},
  {"Both",               10},
  {"Seq-Pept-Overlap",   11},
  {"Seq-Pept-Homol",     12},
  {"Concept-Trans-A",    13},
  {"HTGS 0",             18},
  {"HTGS 1",             14},
  {"HTGS 2",             15},
  {"HTGS 3",             16},
  {"FLI_cDNA",           17},
  {"HTC",                19},
  {"WGS",                20},
  {"Barcode",            21},
  {"Composite-WGS-HTGS", 22},
  {"TSA",                MI_TECH_tsa},
  {"Other:",            255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_complete_alist)
  {"Any",     254},
  {" ",         0},
  {"Complete",  1},
  {"Partial",   2},
  {"No Left",   3},
  {"No Right",  4},
  {"No Ends",   5},
  {"Other",   255},
END_ENUM_ALIST

static ENUM_ALIST(mol_alist)
{"Any",             Seq_mol_other - 1},
{" ",               0},
{"DNA",             Seq_mol_dna},
{"RNA",             Seq_mol_rna},
{"Protein",         Seq_mol_aa},
{"Nucleotide",      Seq_mol_na},
{"Other",           Seq_mol_other},
END_ENUM_ALIST

static ENUM_ALIST(topology_alist)
{"Any",             254},
{" ",               0},
{"Linear",          TOPOLOGY_LINEAR},
{"Circular",        TOPOLOGY_CIRCULAR},
{"Tandem",          TOPOLOGY_TANDEM},
{"Other",           255},
END_ENUM_ALIST

static ENUM_ALIST(strand_alist)
{"Any",             254},
{" ",               Seq_strand_unknown},
{"Single",          Seq_strand_plus},
{"Double",          Seq_strand_minus},
{"Mixed",           Seq_strand_both},
{"Mixed Rev",       Seq_strand_both_rev},
{"Other",           Seq_strand_other},
END_ENUM_ALIST

static ValNodePtr MakeValNodeListFromEnum ( EnumFieldAssocPtr al)
{
  EnumFieldAssocPtr efap;
  ValNodePtr        list;

  efap = al;
  list = NULL;
  while (efap->name != NULL)
  {
    ValNodeAddStr (&list, efap->value, StringSave (efap->name));
    efap ++;
  }
  return list;
}

static void ReplaceValNodeString (ValNodePtr list, CharPtr find, CharPtr repl)
{
  ValNodePtr vnp;

  if (list == NULL || find == NULL || repl == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (StringCmp (vnp->data.ptrvalue, find) == 0)
    {
      MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = StringSave (repl);
    }
  }
}

extern void InitValNodePopup (ValNodePtr list, PopuP p)
{
  ValNodePtr vnp;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      PopupItem (p, vnp->data.ptrvalue);
    }
  }
}

extern Int2 GetValNodePopup (PopuP p, ValNodePtr list)
{
  ValNodePtr vnp;
  Int2       popupval;

  popupval = GetValue (p);
  if (popupval == 0)
  {
    return -1;
  }
  for (vnp = list; vnp != NULL && popupval > 1; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      popupval --;
    }
  }
  if (popupval > 1 || vnp == NULL) return -1;
  return vnp->choice;
}

extern void SetValNodePopupValue (ValNodePtr list, PopuP p, CharPtr val)
{
  ValNodePtr vnp;
  Int2       popupval;

  popupval = 1;
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      if (StringCmp (vnp->data.ptrvalue, val) == 0)
      {
        SetValue (p, popupval);
      }
      popupval ++;
    }
  }
}

typedef struct molinfoblockx {
  PopuP           moltype;
  ValNodePtr      moltype_values;
  PopuP           technique;
  ValNodePtr      technique_values;
  PopuP           complete;
  ValNodePtr      complete_values;
  PopuP           molPopup;
  ValNodePtr      mol_values;
  PopuP           topologyPopup;
  ValNodePtr      topology_values;
  PopuP           strandPopup;
  ValNodePtr      strand_values;
} MolInfoBlock, PNTR MolInfoBlockPtr;

static void FreeMolInfoBlockData (MolInfoBlockPtr mibp)
{
  if (mibp == NULL) return;
  ValNodeFreeData (mibp->moltype_values);
  mibp->moltype_values = NULL;
  ValNodeFreeData (mibp->technique_values);
  mibp->technique_values = NULL;
  ValNodeFreeData (mibp->complete_values);
  mibp->complete_values = NULL;
  ValNodeFreeData (mibp->mol_values);
  mibp->mol_values = NULL;
  ValNodeFreeData (mibp->topology_values);
  mibp->topology_values = NULL;
  ValNodeFreeData (mibp->strand_values);
  mibp->strand_values = NULL;
}


extern SequenceConstraintXPtr SequenceConstraintXFree (SequenceConstraintXPtr scp)
{
  if (scp != NULL)
  {
    scp->string_constraint = StringConstraintXFree (scp->string_constraint);
    scp->source_constraint = ChoiceConstraintFree (scp->source_constraint);
    scp->feature_list = ValNodeFree (scp->feature_list);
    scp = MemFree (scp);
  }
  return scp;
}


typedef struct sequenceconstraintdlg 
{
  DIALOG_MESSAGE_BLOCK
  GrouP  seq_type_constraint;
  GrouP  other_constraint_type;
  DialoG source_constraint;
  DialoG feature_list;
  DialoG feature_string_constraint;
  DialoG id_string_constraint;
  
} SequenceConstraintXDlgData, PNTR SequenceConstraintXDlgPtr;


enum sequenceconstraintothertype 
{
  SEQ_CONSTRAINT_ANY = 0,
  SEQ_CONSTRAINT_SOURCE,
  SEQ_CONSTRAINT_FEATURE_TEXT,
  SEQ_CONSTRAINT_ID
};


static void ChangeSequenceConstraintXType (GrouP g)
{
  SequenceConstraintXDlgPtr scdp;
  Int4                     other_constraint_type;
  
  scdp = (SequenceConstraintXDlgPtr) GetObjectExtra (g);
  if (scdp == NULL)
  {
    return;
  }
  other_constraint_type = (GetValue (scdp->other_constraint_type) - 1) / 2;
  Disable (scdp->source_constraint);
  Disable (scdp->feature_string_constraint);
  Disable (scdp->id_string_constraint);
  if (other_constraint_type == SEQ_CONSTRAINT_SOURCE)
  {
    Enable (scdp->source_constraint);
  }
  else if (other_constraint_type == SEQ_CONSTRAINT_FEATURE_TEXT)
  {
    Enable (scdp->feature_string_constraint);
  }
  else if (other_constraint_type == SEQ_CONSTRAINT_ID)
  {
    Enable (scdp->id_string_constraint);
  }
}

static void ResetSequenceConstraintXDialog (DialoG d)
{
  SequenceConstraintXDlgPtr scdp;
  
  scdp = (SequenceConstraintXDlgPtr) GetObjectExtra (d);
  if (scdp == NULL)
  {
    return;
  }
  
  SetValue (scdp->seq_type_constraint, 1);
  SetValue (scdp->other_constraint_type, SEQ_CONSTRAINT_ANY + 1);
  ChangeSequenceConstraintXType (scdp->other_constraint_type);
}

static void SequenceConstraintXToDialog (DialoG d, Pointer data)
{
  SequenceConstraintXPtr    scp;
  SequenceConstraintXDlgPtr scdp;
  
  scdp = (SequenceConstraintXDlgPtr) GetObjectExtra (d);
  if (scdp == NULL)
  {
    return;
  }
  
  scp = (SequenceConstraintXPtr) data;
  
  if (scp == NULL)
  {
    ResetSequenceConstraintXDialog (d);
  }
  else
  {
    if (scp->nucs_ok && scp->prots_ok)
    {
      SetValue (scdp->seq_type_constraint, 1);
    }
    else if (scp->nucs_ok)
    {
      SetValue (scdp->seq_type_constraint, 2);
    }
    else
    {
      SetValue (scdp->seq_type_constraint, 3);
    }
    SetValue (scdp->other_constraint_type, scp->other_constraint_type + 1);
    switch (scp->other_constraint_type)
    {
      case SEQ_CONSTRAINT_SOURCE:
          PointerToDialog (scdp->source_constraint, scp->source_constraint);
          break;
      case SEQ_CONSTRAINT_FEATURE_TEXT:
          PointerToDialog (scdp->feature_string_constraint, scp->string_constraint);
          PointerToDialog (scdp->feature_list, scp->feature_list);
          break;
      case SEQ_CONSTRAINT_ID:
          PointerToDialog (scdp->id_string_constraint, scp->string_constraint);
          break;
    }
    ChangeSequenceConstraintXType (scdp->other_constraint_type);
  }
}


static Pointer DialogToSequenceConstraintX (DialoG d)
{
  SequenceConstraintXPtr    scp;
  SequenceConstraintXDlgPtr scdp;
  
  scdp = (SequenceConstraintXDlgPtr) GetObjectExtra (d);
  if (scdp == NULL)
  {
    return NULL;
  }
  
  scp = (SequenceConstraintXPtr) MemNew (sizeof (SequenceConstraintXData));
  if (scp != NULL)
  {
    switch (GetValue (scdp->seq_type_constraint))
    {
      case 1:
        scp->nucs_ok = TRUE;
        scp->prots_ok = TRUE;
        break;
      case 2:
        scp->nucs_ok = TRUE;
        break;
      case 3:
        scp->prots_ok = TRUE;
        break;
    }
    
    scp->other_constraint_type = (GetValue (scdp->other_constraint_type) - 1)/2;
    switch (scp->other_constraint_type)
    {
      case SEQ_CONSTRAINT_SOURCE:
        scp->source_constraint = DialogToPointer (scdp->source_constraint);
        break;
      case SEQ_CONSTRAINT_FEATURE_TEXT:
        scp->string_constraint = DialogToPointer (scdp->feature_string_constraint);
        scp->feature_list = DialogToPointer (scdp->feature_list);
        break;
      case SEQ_CONSTRAINT_ID:
        scp->string_constraint = DialogToPointer (scdp->id_string_constraint);
        break;
    }
  }
  return scp;
}


static void SequenceConstraintXMessage (DialoG d, Int2 mssg)

{
  SequenceConstraintXDlgPtr scdp;

  scdp = (SequenceConstraintXDlgPtr) GetObjectExtra (d);
  if (scdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetSequenceConstraintXDialog (d);        
        break;
      case VIB_MSG_ENTER :
        Select (scdp->other_constraint_type);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSequenceConstraintXDialog (DialoG d)

{
  return NULL;
}


extern DialoG SequenceConstraintXDialog (GrouP g)
{
  SequenceConstraintXDlgPtr dlg;
  GrouP p, k;
  
  dlg = (SequenceConstraintXDlgPtr) MemNew (sizeof (SequenceConstraintXDlgData));
  if (dlg == NULL) return NULL;

  p = HiddenGroup (g, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceConstraintXToDialog;
  dlg->fromdialog = DialogToSequenceConstraintX;
  dlg->dialogmessage = SequenceConstraintXMessage;
  dlg->testdialog = TestSequenceConstraintXDialog;

  dlg->seq_type_constraint = HiddenGroup (p, 4, 0, NULL);
  RadioButton (dlg->seq_type_constraint, "All Sequences");
  RadioButton (dlg->seq_type_constraint, "Nucleotides");
  RadioButton (dlg->seq_type_constraint, "Proteins");
  SetValue (dlg->seq_type_constraint, 2);

  dlg->other_constraint_type = HiddenGroup (p, 2, 0, ChangeSequenceConstraintXType);
  RadioButton (dlg->other_constraint_type, "Any sequence");
  StaticPrompt (dlg->other_constraint_type, "", 0, dialogTextHeight, systemFont, 'l');  
  
  RadioButton (dlg->other_constraint_type, "Where sequence source");
  dlg->source_constraint = SourceConstraintDialogX (dlg->other_constraint_type, FALSE);
    
  RadioButton (dlg->other_constraint_type, "Where feature text");
  k = HiddenGroup (dlg->other_constraint_type, 2, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  dlg->feature_list = FeatureSelectionDialog (k, TRUE, NULL, NULL);
  dlg->feature_string_constraint = StringConstraintDialogX (k, NULL, FALSE);
    
  RadioButton (dlg->other_constraint_type, "Where sequence ID");
  dlg->id_string_constraint = StringConstraintDialogX (dlg->other_constraint_type, NULL, FALSE);
  
  SetValue (dlg->other_constraint_type, CHOICE_CONSTRAINT_ANY);
  SetObjectExtra (dlg->other_constraint_type, dlg, NULL);

  Disable (dlg->source_constraint);
  Disable (dlg->feature_string_constraint);
  Disable (dlg->id_string_constraint);
    
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->seq_type_constraint, 
                              (HANDLE) dlg->other_constraint_type, NULL);
   
  return (DialoG) p;
    
}


extern Boolean DoesIDListMeetStringConstraint (SeqIdPtr sip, StringConstraintXPtr string_constraint)
{
  Char       id [41];
  CharPtr    cp, cp_dst;
  SeqIdPtr   tmp;
  Boolean    match, changed;

  if (sip == NULL) 
  {
    return FALSE;
  }
  if (string_constraint == NULL)
  {
    return TRUE;
  }
  
  while (sip != NULL)
  {
    /* temporary disconnect ID from list */
    tmp = sip->next;
    sip->next = NULL;
    id [0] = '\0';
    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id) - 1);
    match = DoesStringMatchConstraintX (id, string_constraint);
    if (!match) 
    {
      changed = FALSE;
      /* remove terminating pipe character */
      if (id[StringLen(id) - 1] == '|') 
      {
        id[StringLen(id) - 1] = 0;
        changed = TRUE;
      }
      /* remove leading pipe identifier */
      cp = StringChr (id, '|');
      if (cp != NULL)
      {
        changed = TRUE;
        cp++;
        cp_dst = id;
        while (*cp != 0) 
        {
          *cp_dst = *cp;
          cp_dst++;
          cp++;
        }
        *cp_dst = 0;
      }  
      if (changed) 
      {
        match = DoesStringMatchConstraintX (id, string_constraint);
      }

      /* if search text doesn't have ., try ID without version */
      if (!match && StringChr (string_constraint->match_text, '.') == NULL) 
      {
        cp = StringChr (id, '.');
        if (cp != NULL) 
        {
          *cp = 0;
          match = DoesStringMatchConstraintX (id, string_constraint);
        }
      }       
    }
    sip->next = tmp;

    if (match)
    {
      if (string_constraint->not_present)
      {
        return FALSE;
      }
      else
      {
        return TRUE;
      }
    }
    sip = sip->next;
  }
  if (string_constraint->not_present)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


extern Boolean DoesSequenceMatchSequenceConstraintX (BioseqPtr bsp, SequenceConstraintXPtr scp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  
  if (bsp == NULL) return FALSE;
  if (scp == NULL) return TRUE;
  if (!scp->nucs_ok && ISA_na (bsp->mol)) return FALSE;
  if (!scp->prots_ok && ISA_aa (bsp->mol)) return FALSE;
  switch (scp->other_constraint_type)
  {
    case SEQ_CONSTRAINT_ANY:
      return TRUE;
      break;
    case SEQ_CONSTRAINT_SOURCE:
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
      if (sdp != NULL)
      {
        return DoesOneSourceMatchConstraint(sdp->data.ptrvalue, scp->source_constraint);
      }
      break;
    case SEQ_CONSTRAINT_FEATURE_TEXT:
      return DoBioseqFeaturesMatchSequenceConstraintX (bsp, scp->feature_list, scp->string_constraint);
      break;
    case SEQ_CONSTRAINT_ID:
      return DoesIDListMeetStringConstraint(bsp->id, scp->string_constraint);
      break;
  }
  return FALSE; 
}

typedef struct molinfoedit {
  DESCRIPTOR_FORM_BLOCK
  SeqEntryPtr     sep;
  MolInfoBlock    from;
  MolInfoBlock    to;
  DialoG          sequence_constraint;
  SequenceConstraintXPtr scp;
  Boolean               allow_existing;
} MolInfoEdit, PNTR MolInfoEditPtr;

static CharPtr  labels [] = {
  "Molecule", "Technique", "Completedness",
  "Class", "Topology", "Strand", NULL
};

static void SetupMolInfoBlockPopup (GrouP g, PopuP PNTR p, ValNodePtr PNTR v,
                                    EnumFieldAssocPtr al, Boolean IsTo)
{
  *p = PopupList (g, TRUE, NULL);
  *v = MakeValNodeListFromEnum ( al );
  if (IsTo)
  {
    ReplaceValNodeString (*v, "Any", "No change");
  }
  InitValNodePopup (*v, *p);
  SetValue (*p, 1);
}

static void CreateMolInfoBlock (MolInfoBlockPtr mibp, GrouP h, Boolean IsTo)

{
  GrouP  g;
  Int2   wid;

  if (mibp == NULL || h == NULL) return;
  SelectFont (programFont);
  wid = MaxStringWidths (labels);
  SelectFont (systemFont);

  g = HiddenGroup (h, -2, 0, NULL);

  StaticPrompt (g, "Molecule", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->moltype), &(mibp->moltype_values),
                          molinfo_biomol_alist, IsTo);

  StaticPrompt (g, "Technique", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->technique), &(mibp->technique_values),
                          molinfo_tech_alist, IsTo);

  StaticPrompt (g, "Completedness", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->complete), &(mibp->complete_values),
                          molinfo_complete_alist, IsTo);

  StaticPrompt (g, "Class", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->molPopup), &(mibp->mol_values),
                          mol_alist, IsTo);

  StaticPrompt (g, "Topology", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->topologyPopup), &(mibp->topology_values),
                          topology_alist, IsTo);

  StaticPrompt (g, "Strand", wid, popupMenuHeight, programFont, 'l');
  SetupMolInfoBlockPopup (g, &(mibp->strandPopup), &(mibp->strand_values),
                          strand_alist, IsTo);
}

static void EditMolInfoCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  MolInfoEditPtr  miep;
  MolInfoPtr      mip;
  ValNodePtr      sdp;
  Uint1           new_mol;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  miep = (MolInfoEditPtr) mydata;
  if (miep == NULL) return;
  sdp = NULL;
  bsp = NULL;
  bssp = NULL;

  if (IS_Bioseq (sep)) {

    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (DoesSequenceMatchSequenceConstraintX (bsp, miep->scp)) {
      if ((Uint1)GetValNodePopup (miep->from.molPopup,
                                  miep->from.mol_values) == bsp->mol
        || GetValue (miep->from.molPopup) == 1)
      {
        if (GetValue (miep->to.molPopup) != 1)
        {
          new_mol = (Uint1) GetValNodePopup (miep->to.molPopup,
                                             miep->to.mol_values);
          if (ISA_na (bsp->mol) && ISA_aa (new_mol))
          {
            BioseqConvert (bsp, Seq_code_ncbieaa);
          }
          else if (ISA_aa (bsp->mol) && ISA_na (new_mol))
          {
            BioseqConvert (bsp, Seq_code_ncbi2na);
          }
          bsp->mol = new_mol;                                              
        }
      }
      if ((Uint1)GetValNodePopup (miep->from.strandPopup,
                                  miep->to.strand_values) == bsp->strand
        || GetValue (miep->from.strandPopup) == 1)
      {
        if (GetValue (miep->to.strandPopup) != 1)
        {
          bsp->strand = (Uint1) GetValNodePopup (miep->to.strandPopup,
                                                 miep->to.strand_values);
        }
      }
      if ((Uint1) GetValNodePopup (miep->from.topologyPopup,
                                   miep->from.topology_values) == bsp->topology
        || GetValue (miep->from.topologyPopup) == 1)
      {
        if (GetValue (miep->to.topologyPopup) != 1)
        {
          bsp->topology = (Uint1) GetValNodePopup (miep->to.topologyPopup,
                                                   miep->to.topology_values);
        }
      }
      sdp = bsp->descr;
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == 2 || bssp->_class == 4) {
      if (miep->scp != NULL && ! miep->scp->nucs_ok) return;
    }
    sdp = bssp->descr;
  } else {
    return;
  }

  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_molinfo && sdp->data.ptrvalue != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if ((Uint2) GetValNodePopup (miep->from.moltype,
                                   miep->from.moltype_values) == mip->biomol
        || GetValue (miep->from.moltype) == 1)
      {
        if ( GetValue (miep->to.moltype) != 1)
        {
          mip->biomol = (Uint1) GetValNodePopup (miep->to.moltype,
                                                 miep->to.moltype_values);
        }
      }
      if ((Uint1) GetValNodePopup (miep->from.technique,
                                   miep->from.technique_values) == mip->tech
        || GetValue (miep->from.technique) == 1)
      {
        if ( GetValue (miep->to.technique) != 1)
        {
          mip->tech = (Uint1) GetValNodePopup (miep->to.technique,
                                     miep->to.technique_values);
        }
      }
      if ((Uint1) GetValNodePopup (miep->from.complete,
                              miep->from.complete_values) == mip->completeness
        || GetValue (miep->from.complete) == 1)
      {
        if ( GetValue (miep->to.complete) != 1)
        {
          mip->completeness = (Uint1) GetValNodePopup (miep->to.complete,
                                                   miep->to.complete_values);
        }
      }
    }
    sdp = sdp->next;
  }
}

static void DoProcessEditMolInfo (ButtoN b)

{
  MolInfoEditPtr  miep;

  miep = (MolInfoEditPtr) GetObjectExtra (b);
  if (miep == NULL) return;
  Hide (miep->form);
  WatchCursor ();
  Update ();
  
  miep->scp = DialogToPointer (miep->sequence_constraint);
  
  SeqEntryExplore (miep->sep, (Pointer) miep, EditMolInfoCallback);
  FreeMolInfoBlockData (&(miep->from));
  FreeMolInfoBlockData (&(miep->to));
  miep->scp = SequenceConstraintXFree (miep->scp);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (miep->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, miep->input_entityID, 0, 0);
  Remove (miep->form);
}

static void EditMolInfoFields (IteM i)

{
  ButtoN          b;
  BaseFormPtr     bfp;
  GrouP           c;
  GrouP           g;
  GrouP           h;
  MolInfoEditPtr  miep;
  SeqEntryPtr     sep;
  WindoW          w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  miep = (MolInfoEditPtr) MemNew (sizeof (MolInfoEdit));
  if (miep == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Molecule Information Editor", StdCloseWindowProc);
  SetObjectExtra (w, miep, StdCleanupFormProc);
  miep->form = (ForM) w;
  miep->formmessage = NULL;

  miep->sep = sep;
  miep->input_entityID = bfp->input_entityID;

  g = HiddenGroup (w, 2, 0, NULL);
  h = HiddenGroup (g, 0, 2, NULL);
  StaticPrompt (h, "From", 0, 0, programFont, 'c');
  CreateMolInfoBlock (&miep->from, h, FALSE);
  h = HiddenGroup (g, 0, 2, NULL);
  StaticPrompt (h, "To", 0, 0, programFont, 'c');
  CreateMolInfoBlock (&miep->to, h, TRUE);

  miep->sequence_constraint = SequenceConstraintXDialog (w);

  c = HiddenGroup (w, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoProcessEditMolInfo);
  SetObjectExtra (b, miep, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) miep->sequence_constraint, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct molinfofrm {
  FORM_MESSAGE_BLOCK
  DialoG d;
  ButtoN leave_dlg_up;
} MolInfoFrmData, PNTR MolInfoFrmPtr;

static void DoMolInfoEdit (ButtoN b)
{
  MolInfoFrmPtr   frm;
  SeqEntryPtr sep;
  MolinfoBlockPtr mib;

  frm = (MolInfoFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }


  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  mib = DialogToPointer (frm->d);
  ApplyMolinfoBlockToSeqEntry (sep, mib); 
  mib = MolinfoBlockFree (mib);

  Update ();
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  if (!GetStatus (frm->leave_dlg_up)) {
    Remove (frm->form);
  }
}

  
static void NewMolInfoFields (IteM i, Boolean edit)
{
  ButtoN          b;
  BaseFormPtr     bfp;
  GrouP           c;
  GrouP           h;
  MolInfoFrmPtr   frm;
  WindoW          w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  frm = (MolInfoFrmPtr) MemNew (sizeof (MolInfoFrmData));
  if (frm == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Molecule Information Editor", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->form = (ForM) w;
  frm->formmessage = NULL;

  frm->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, -1, 0, NULL);

  frm->d = MolInfoBlockDialog (h, edit, NULL, NULL);
  PointerToDialog (frm->d, NULL);

  c = HiddenGroup (w, 3, 0, NULL);
  b = DefaultButton (c, "Accept", DoMolInfoEdit);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  frm->leave_dlg_up = CheckBox (c, "Leave dialog up", NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->d, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


static void NewMolInfoFieldsApply (IteM i)
{
  NewMolInfoFields (i, FALSE);
}


static void NewMolInfoFieldsEdit (IteM i)
{
  NewMolInfoFields (i, TRUE);
}



static void RibosomalRNAToGenomicDNAVisitFunc (
  BioseqPtr bsp,
  Pointer userdata
)
{
  SeqDescrPtr sdp;
  MolInfoPtr mip;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL 
    || sdp->choice != Seq_descr_molinfo
    || sdp->data.ptrvalue == NULL)
  {
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol == MOLECULE_TYPE_RRNA)
  {
    mip->biomol = MOLECULE_TYPE_GENOMIC;
    bsp->mol = MOLECULE_CLASS_DNA;
    bsp->strand = 0;
    bsp->topology = 0;
  }
}

static void RibosomalRNAToGenomicDNA (BaseFormPtr bfp)
{
  SeqEntryPtr     sep;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, RibosomalRNAToGenomicDNAVisitFunc);

  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RibosomalRNAToGenomicDNAMenuItem (IteM i)
{
  BaseFormPtr     bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  RibosomalRNAToGenomicDNA (bfp);
}

static void RibosomalRNAToGenomicDNAToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  
  RibosomalRNAToGenomicDNA (bfp);
}

static Boolean MolInfoOkForSeqEntry (SeqEntryPtr sep, SequenceConstraintXPtr scp)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  seq_set;

  if (IS_Bioseq (sep)) 
  {
    return DoesSequenceMatchSequenceConstraintX(sep->data.ptrvalue, scp);
  } 
  else if (IS_Bioseq_set (sep)) 
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL 
        || bssp->seq_set == NULL)
    {
      return FALSE;
    }
    else if (bssp->_class == BioseqseqSet_class_segset)
    {
      for (seq_set = bssp->seq_set; seq_set != NULL; seq_set = seq_set->next)
      {
        if (!MolInfoOkForSeqEntry (seq_set, scp))
        {
          return FALSE;
        }
      }
      return TRUE;
    }
    else if (bssp->_class == BioseqseqSet_class_parts)
    {
      for (seq_set = bssp->seq_set; seq_set != NULL; seq_set = seq_set->next)
      {
        if (!MolInfoOkForSeqEntry (seq_set, scp))
        {
          return FALSE;
        }
      }
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  else
  {
    return FALSE;
  }
}

static void DoApplyMolInfoToBioseq (BioseqPtr bsp, MolInfoEditPtr miep)
{
  if (bsp != NULL && miep != NULL) 
  {
    if ( GetValue (miep->to.molPopup) != 1)
    {
      bsp->mol = (Uint1) GetValNodePopup (miep->to.molPopup,
                                          miep->to.mol_values);
    }
    if ( GetValue (miep->to.strandPopup) != 1)
    {
      bsp->strand = (Uint1) GetValNodePopup (miep->to.strandPopup,
                                             miep->to.strand_values);
    }
    if (GetValue (miep->to.topologyPopup) != 1)
    {
      bsp->topology = (Uint1) GetValNodePopup (miep->to.topologyPopup,
                                               miep->to.topology_values);
    }
  }
}

static void DoApplyMolInfoToBioseqSet (BioseqSetPtr bssp, MolInfoEditPtr miep)
{
  SeqEntryPtr sep;
  if (bssp != NULL && miep != NULL) 
  {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
      if (IS_Bioseq (sep)) 
      {
        DoApplyMolInfoToBioseq (sep->data.ptrvalue, miep);
      }
      else if (IS_Bioseq_set (sep)) 
      {
        DoApplyMolInfoToBioseqSet (sep->data.ptrvalue, miep);
      }
    }
  }
}


static void DoApplyMolInfo (SeqEntryPtr sep, MolInfoEditPtr miep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  MolInfoPtr    mip = NULL;
  SeqDescrPtr   sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL || miep == NULL) return;
  if (IS_Bioseq_set (sep)) 
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
/* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        DoApplyMolInfo (sep, miep);
      }
      return;
    }
  }
  else if (!MolInfoOkForSeqEntry(sep, miep->scp))
  {
    return;
  }

  bsp = NULL;
  bssp = NULL;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  if (!miep->allow_existing && sdp != NULL) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (DoesSequenceMatchSequenceConstraintX(bsp, miep->scp)) {
      DoApplyMolInfoToBioseq (bsp, miep);
    } else return;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    DoApplyMolInfoToBioseqSet (bssp, miep);
  } else return;
  if (sdp != NULL) 
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }
  if (mip == NULL) 
  {
    mip = MolInfoNew ();
  }
  if (mip == NULL) return;
  if (sdp == NULL) 
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
  }
  if (sdp == NULL) return;
  sdp->data.ptrvalue = mip;

  if ( GetValue (miep->to.moltype) != 1)
  {
    mip->biomol = (Uint1) GetValNodePopup (miep->to.moltype,
                                           miep->to.moltype_values);
  }
  if ( GetValue (miep->to.technique) != 1)
  {
    mip->tech = (Uint1) GetValNodePopup (miep->to.technique,
                                         miep->to.technique_values);
  }
  if ( GetValue (miep->to.complete) != 1)
  {
    mip->completeness = (Uint1) GetValNodePopup (miep->to.complete,
                                                 miep->to.complete_values);
  }
}

static Boolean WillEditMolInfo (MolInfoEditPtr miep)
{
  if (miep == NULL) 
  {
    return FALSE;
  }
  else if (GetValue (miep->to.moltype) != 1
           || GetValue (miep->to.technique) != 1
           || GetValue (miep->to.complete) != 1) 
  {
    return TRUE;
  } 
  else
  {
    return FALSE;
  }
}

typedef struct molinfocheck {
  SequenceConstraintXPtr scp;
  Int4                  num_with;
  Int4                  num_without;
} MolInfoCheckData, PNTR MolInfoCheckPtr;

static void CheckForExistingMolInfo (SeqEntryPtr sep, MolInfoCheckPtr micp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqDescrPtr   sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL || micp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
/* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CheckForExistingMolInfo (sep, micp);
      }
      return;
    }
  }
  bsp = NULL;
  bssp = NULL;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  if (MolInfoOkForSeqEntry (sep, micp->scp)) 
  {
    if (sdp == NULL) 
    {
      micp->num_without++;
    }
    else
    {
      micp->num_with++;
    }
  }
}


static Int4 
ThreeOptionsDlg
(CharPtr title_txt,
 CharPtr explain_txt,
 CharPtr opt1,
 CharPtr opt2,
 CharPtr opt3)
{
  WindoW w;
  GrouP  h, c;
  ButtoN b;
  GrouP  p;
  ModalAcceptCancelData acd;
  
  
  w = MovableModalWindow (-20, -13, -10, -10, title_txt, NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  p = MultiLinePrompt (h, explain_txt, 300, programFont);
    
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, opt1, ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, opt2, ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  if (!StringHasNoText (opt3)) {
    b = PushButton (c, opt3, ModalThirdOptionButton);
    SetObjectExtra (b, &acd, NULL);
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p,
                              (HANDLE) c, 
                              NULL);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  acd.third_option = FALSE;
  acd.accepted = FALSE;
  while (!acd.accepted && ! acd.cancelled && !acd.third_option)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  
  if (acd.accepted)
  {
    return 1;
  }
  else if (acd.cancelled)
  {
    return 2;
  }
  else
  {
    return 3;
  }
}

static void DoProcessApplyMolInfo (ButtoN b)

{
  MolInfoEditPtr   miep;
  MolInfoCheckData micd;
  Int4             opt;
  Char             err_msg[1000];
  CharPtr          err_fmt = "%d sequences that match your constraint already have MolInfo descriptors.";

  miep = (MolInfoEditPtr) GetObjectExtra (b);
  if (miep == NULL) return;
  Hide (miep->form);
  WatchCursor ();
  Update ();
  miep->scp = DialogToPointer (miep->sequence_constraint);
  micd.scp = miep->scp;
  micd.num_with = 0;
  micd.num_without = 0;
  
  CheckForExistingMolInfo (miep->sep, &micd);
  if (micd.num_with == 0 && micd.num_without == 0) 
  {
    Message (MSG_ERROR, "No sequences meet your constraint!");
    return;
  }

  if (WillEditMolInfo (miep))
  {
    miep->allow_existing = FALSE;
    if (micd.num_without == 0) 
    {
      opt = ThreeOptionsDlg ("MolInfo Descriptors Already Exist",
                             "All sequences that match your constraint already have MolInfo descriptors.",
                             "Edit existing MolInfo descriptors",
                             "Cancel",
                             NULL);
      if (opt == 2) 
      {
        return;
      }
      else
      {
        miep->allow_existing = TRUE;
      }
    } 
    else
    {                      
      sprintf (err_msg, err_fmt, micd.num_with);
      opt = ThreeOptionsDlg ("MolInfo Descriptors Already Exist", err_msg, 
                             "Continue and edit existing MolInfo descriptors",
                             "Continue, but skip existing MolInfo descriptors",
                             "Cancel");
      if (opt == 3) 
      {
        return;
      } 
      else if (opt == 1)
      {
        miep->allow_existing = TRUE;
      }
    }
  }
  else
  {
    miep->allow_existing = TRUE;
  }

  DoApplyMolInfo (miep->sep, miep);
  FreeMolInfoBlockData (&miep->to);
  miep->scp = SequenceConstraintXFree (miep->scp);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (miep->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, miep->input_entityID, 0, 0);
  Remove (miep->form);
}

static void ApplyMolInfo (IteM i)

{
  ButtoN          b;
  BaseFormPtr     bfp;
  GrouP           c;
  GrouP           g;
  GrouP           h;
  MolInfoEditPtr  miep;
  SeqEntryPtr     sep;
  WindoW          w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  miep = (MolInfoEditPtr) MemNew (sizeof (MolInfoEdit));
  if (miep == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Apply Molecule Information", StdCloseWindowProc);
  SetObjectExtra (w, miep, StdCleanupFormProc);
  miep->form = (ForM) w;
  miep->formmessage = NULL;

  miep->sep = sep;
  miep->input_entityID = bfp->input_entityID;

  g = HiddenGroup (w, 2, 0, NULL);
  h = HiddenGroup (g, 0, 2, NULL);
  CreateMolInfoBlock (&miep->to, h, TRUE);

  miep->sequence_constraint = SequenceConstraintXDialog (w);

  c = HiddenGroup (w, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoProcessApplyMolInfo);
  SetObjectExtra (b, miep, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) miep->sequence_constraint, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

typedef struct objstringdata 
{
  CharPtr match;
  Boolean found;
  Boolean insensitive;
  Boolean whole_word;
} ObjStringData, PNTR ObjStringPtr;


static Boolean StringMatchesConstraint (CharPtr pchSource, ObjStringPtr osp)
{
  CharPtr pchFind;
  CharPtr pFound;
  
  if (pchSource == NULL || osp == NULL) return FALSE;
  
	pchFind = osp->match;
	  
	if (osp->insensitive)
	{
	  pFound = StringISearch (pchSource, pchFind);
	}
	else
	{
	  pFound = StringSearch (pchSource, pchFind);
	}
	  
	if (pFound != NULL)
	{
	  if (!osp->whole_word)
	  {
	    return TRUE;
	  }
	  else if ((pFound == pchSource || (! isalpha ((Int4)(*(pFound - 1))) && ! isdigit ((Int4)(*(pFound - 1)))))
	           && ! isalpha ((Int4)(*(pFound + StringLen (pchFind))))
	           && ! isdigit ((Int4)(*(pFound + StringLen (pchFind)))))
	  {
	    return TRUE;
	  }
	}
	return FALSE;
}

static void LIBCALLBACK AsnWriteSearchCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr        pchSource;
  ObjStringPtr   osp;

  osp = (ObjStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) 
  {
	  pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	  osp->found |= StringMatchesConstraint (pchSource, osp);
  }
}


static Boolean ObjectHasSubstring (ObjMgrTypePtr omtp, AsnIoPtr aip, Pointer ptr, ObjStringPtr osp)

{
  osp->found = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  return osp->found;
}

extern void OperateOnBioseqFeaturesWithText 
(BioseqPtr         bsp,
 Pointer           userdata)
{
  AsnExpOptPtr            aeop;
  AsnIoPtr                aip;
  ObjStringData           osd;
  SeqFeatPtr              sfp;
  ObjMgrPtr               omp;
  ObjMgrTypePtr           omtp;
  SeqMgrFeatContext       fcontext;
  FeaturesWithTextPtr     fdp;
  
  fdp = (FeaturesWithTextPtr) userdata;
  if (bsp == NULL || fdp == NULL || fdp->callback == NULL || fdp->search_text == NULL) return; 
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp == NULL) return;

  aip = AsnIoNullOpen ();
  osd.insensitive = fdp->case_insensitive;
  osd.whole_word = fdp->whole_word;
  
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteSearchCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &osd;
  }
  osd.match = fdp->search_text;

  sfp = NULL;
  while ((sfp = SeqMgrGetNextFeature (bsp, sfp, fdp->seqFeatChoice, fdp->featDefChoice, &fcontext)) != NULL)
  {
    if (fdp->no_text || ObjectHasSubstring (omtp, aip, (Pointer) sfp, &osd)) 
    {
      fdp->callback (sfp, fdp->userdata);  
    }
    else if (StringMatchesConstraint(fcontext.label, &osd))
    {
      fdp->callback (sfp, fdp->userdata);  
    }
  }
  AsnIoClose (aip);
}

typedef struct seqentryitemswithtext 
{
  FeaturesWithTextPtr    fdp;
  DescriptorsWithTextPtr ddp;
  ObjMgrTypePtr          omtp;
  AsnIoPtr               aip;
  Uint2                  entityID;
  ObjStringData          osd;
} SeqEntryItemsWithTextData, PNTR SeqEntryItemsWithTextPtr;

static void SeqEntryFeaturesWithTextCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqEntryItemsWithTextPtr sp;
  
  if (sfp == NULL || userdata == NULL) return;
  sp = (SeqEntryItemsWithTextPtr) userdata;
  sfp = SeqMgrGetDesiredFeature (sp->entityID, NULL, 0, 0, sfp, &fcontext);
  if (sfp == NULL) return;
  
  if (sp->fdp->seqFeatChoice != 0 && sp->fdp->seqFeatChoice != sfp->data.choice) return;
  if (sp->fdp->featDefChoice != 0 && sp->fdp->featDefChoice != sfp->idx.subtype) return;
  
  if (sp->fdp->no_text
      || (! sp->fdp->act_when_string_not_present &&
         (ObjectHasSubstring (sp->omtp, sp->aip, (Pointer) sfp, &(sp->osd))
          || StringMatchesConstraint (fcontext.label, &(sp->osd))))
      || (sp->fdp->act_when_string_not_present
          && ! ObjectHasSubstring (sp->omtp, sp->aip, (Pointer) sfp, &(sp->osd))
          && ! StringMatchesConstraint (fcontext.label, &(sp->osd))))
  {
    sp->fdp->callback (sfp, sp->fdp->userdata);  
  }
}

extern void OperateOnSeqEntryFeaturesWithText (SeqEntryPtr sep, FeaturesWithTextPtr fdp)
{
  SeqEntryItemsWithTextData sd;
  ObjMgrPtr                 omp;
  AsnExpOptPtr              aeop;
  
  if (sep == NULL || fdp == NULL) return;
  
  sd.fdp = fdp;
  sd.entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  sd.omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
  if (sd.omtp == NULL) return;

  sd.aip = AsnIoNullOpen ();
  sd.osd.insensitive = fdp->case_insensitive;
  sd.osd.whole_word = fdp->whole_word;
  aeop = AsnExpOptNew (sd.aip, NULL, NULL, AsnWriteSearchCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &(sd.osd);
  }
  sd.osd.match = fdp->search_text;
  
  VisitFeaturesInSep (sep, &sd, SeqEntryFeaturesWithTextCallback);
  AsnIoClose (sd.aip);
}


static void SeqEntryDescriptorsWithTextCallback (SeqDescrPtr sdp, Pointer userdata)
{
  SeqEntryItemsWithTextPtr sp;
  
  if (sdp == NULL || userdata == NULL) return;
  sp = (SeqEntryItemsWithTextPtr) userdata;
  if (sp->ddp->no_text 
      || (! sp->ddp->act_when_string_not_present && ObjectHasSubstring (sp->omtp, sp->aip, (Pointer) sdp, &(sp->osd)))
      || (sp->ddp->act_when_string_not_present && ! ObjectHasSubstring (sp->omtp, sp->aip, (Pointer) sdp, &(sp->osd))))
  {
    sp->ddp->callback (sdp, sp->ddp->userdata);  
  }
}


extern void OperateOnSeqEntryDescriptorsWithText (SeqEntryPtr sep, DescriptorsWithTextPtr ddp)
{
  SeqEntryItemsWithTextData sd;
  ObjMgrPtr                 omp;
  AsnExpOptPtr              aeop;
  
  if (sep == NULL || ddp == NULL) return;
  
  sd.ddp = ddp;
  sd.entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  sd.omtp = ObjMgrTypeFind (omp, OBJ_SEQDESC, NULL, NULL);
  if (sd.omtp == NULL) return;

  sd.aip = AsnIoNullOpen ();
  sd.osd.insensitive = ddp->case_insensitive;
  sd.osd.whole_word = ddp->whole_word;
  aeop = AsnExpOptNew (sd.aip, NULL, NULL, AsnWriteSearchCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &(sd.osd);
  }
  sd.osd.match = ddp->search_text;
  
  VisitDescriptorsInSep (sep, &sd, SeqEntryDescriptorsWithTextCallback);
  AsnIoClose (sd.aip);
}


static void GenePseudoOn (SeqFeatPtr sfp, Pointer userdata)

{
  GeneRefPtr  grp;

  if (sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
  grp->pseudo = TRUE;
}

static void GenePseudoOff (SeqFeatPtr sfp, Pointer userdata)

{
  GeneRefPtr  grp;

  if (sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
  grp->pseudo = FALSE;
}

typedef struct genepseudodata 
{
  FORM_MESSAGE_BLOCK
  Boolean set_on;
  TexT    search_text;  
  ButtoN  case_insensitive;
  ButtoN  when_string_not_present;
} GenePseudoData, PNTR GenePseudoPtr;

static void DoGenePseudoConstraint (ButtoN b)
{
  GenePseudoPtr        gsp;
  Char                 search_text[255];
  FeaturesWithTextData fd;
  SeqEntryPtr          sep;
  
  gsp = (GenePseudoPtr) GetObjectExtra (b);
  if (gsp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (gsp->input_entityID);
  if (sep == NULL) return;
  
  GetTitle (gsp->search_text, search_text, sizeof (search_text) - 1);
  
  /* set up text operation */
  fd.seqFeatChoice = SEQFEAT_GENE;
  fd.featDefChoice = 0;
  fd.search_text = search_text;
  fd.no_text = StringHasNoText (search_text);
  fd.case_insensitive = GetStatus (gsp->case_insensitive);
  fd.whole_word = FALSE;
  fd.act_when_string_not_present = GetStatus (gsp->when_string_not_present);
  fd.userdata = NULL;
  if (gsp->set_on)
  {
    fd.callback = GenePseudoOn;
  }
  else 
  {
  	fd.callback = GenePseudoOff;
  }
  
  OperateOnSeqEntryFeaturesWithText (sep, &fd);
  ObjMgrSetDirtyFlag (gsp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, gsp->input_entityID, 0, 0);  
  Update ();
  Remove (gsp->form);
}

static void GenePseudoConstraint (IteM i, Boolean set_on)
{
  BaseFormPtr        bfp;
  WindoW             w;
  GenePseudoPtr      gsp;
  GrouP              h, g, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  gsp = (GenePseudoPtr) MemNew (sizeof (GenePseudoData));
  if (gsp == NULL) return;

  gsp->input_entityID = bfp->input_entityID;
  gsp->set_on = set_on;
  if (gsp->set_on)
  {
    w = FixedWindow (-50, -33, -10, -10, "Set Gene Pseudo", StdCloseWindowProc);
  }
  else
  {
    w = FixedWindow (-50, -33, -10, -10, "Clear Gene Pseudo", StdCloseWindowProc);
  }
  SetObjectExtra (w, gsp, StdCleanupFormProc);
  gsp->form = (ForM) w;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 4, NULL);
  StaticPrompt (g, "String Constraint", 0, 0, programFont, 'c');
  gsp->search_text = DialogText (g, "", 20, NULL);
  gsp->case_insensitive = CheckBox (g, "Case Insensitive", NULL);
  gsp->when_string_not_present = CheckBox (g, "When String Not Present", NULL);
  
  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoGenePseudoConstraint);
  SetObjectExtra (b, gsp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  RealizeWindow (w);
  Show (w);
  Select (gsp->search_text);
  Update ();
}

static void ClearGenePseudoConstraint (IteM i)
{
  GenePseudoConstraint (i, FALSE);	
}


static void AddGenePseudo (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, GenePseudoOn);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ClearGenePseudo (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, GenePseudoOff);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean ClearOrfCallback (GatherContextPtr gcp)

{
  CdRegionPtr  crp;
  SeqFeatPtr   sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;
  crp->orf = FALSE;
  return TRUE;
}

static void ClearOrfFlagInCDS (IteM i)

{
  BaseFormPtr  bfp;
  GatherScope  gs;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, NULL, ClearOrfCallback, &gs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean ClearTranslExceptCallback (GatherContextPtr gcp)

{
  CdRegionPtr  crp;
  SeqFeatPtr   sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;
  crp->code_break = CodeBreakFree (crp->code_break);
  return TRUE;
}

static void ClearTranslExceptInCDS (IteM i)

{
  BaseFormPtr  bfp;
  GatherScope  gs;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, NULL, ClearTranslExceptCallback, &gs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean ConvertTranslExceptCallback (GatherContextPtr gcp)

{
  CdRegionPtr  crp;
  SeqFeatPtr   sfp;
  CharPtr      text = "RNA editing";

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL || crp->code_break == NULL) return TRUE;
  crp->code_break = CodeBreakFree (crp->code_break);
  sfp->excpt = TRUE;
  if (sfp->except_text != NULL)
  {
  	MemFree (sfp->except_text);
  }
  sfp->except_text = (CharPtr) MemNew (sizeof (Char) * (StringLen (text) + 1));
  if (sfp->except_text != NULL)
  {
  	StringCpy (sfp->except_text, text);
  }
  return TRUE;
}

static void ConvertTranslExceptToRNAEditingException (IteM i)

{
  BaseFormPtr  bfp;
  GatherScope  gs;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, NULL, ConvertTranslExceptCallback, &gs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);	
}

typedef struct bioseqsetform {
  FEATURE_FORM_BLOCK
  BioseqSetPtr       bssp;
  Uint2              entityID;
  PopuP              class_control;
} BioseqSetForm, PNTR BioseqSetFormPtr;

static ENUM_ALIST(bioseqset_class_alist)
  {"  ",                     0},
  {"Nuc-prot",               1},
  {"Segset",                 2},
  {"Conset",                 3},
  {"Parts",                  4},
  {"Gibb",                   5},
  {"GI",                     6},
  {"Genbank",                7},
  {"PIR",                    8},
  {"Pubset",                 9},
  {"Equiv",                 10},
  {"Swissprot",             11},
  {"PDB-entry",             12},
  {"Mut-set",               13},
  {"Pop-set",               14},
  {"Phy-set",               15},
  {"Eco-set",               16},
  {"Gen-prod-set",          17},
  {"WGS-set",               18},
  {"Other",                255},
END_ENUM_ALIST

static void AcceptBioseqSetEditProc (ButtoN b)

{
  BioseqSetFormPtr  bsfp;
  BioseqSetPtr      bssp;
  UIEnum            val;

  bsfp = (BioseqSetFormPtr) GetObjectExtra (b);
  if (bsfp != NULL) {
    Hide (bsfp->form);
    bssp = bsfp->bssp;
    if (bssp == NULL && bsfp->entityID == 0) {
      bssp = BioseqSetNew ();
    }
    if (bssp != NULL) {
      GetEnumPopup (bsfp->class_control, bioseqset_class_alist, &val);
      bssp->_class = (Uint1) val;
      if (bsfp->entityID == 0) {
        if (! ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp)) {
          Message (MSG_ERROR, "ObjMgrRegister failed");
        }
      }
    }
    Remove (bsfp->form);
    ObjMgrSetDirtyFlag (bsfp->entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bsfp->entityID, 0, 0);
  }
}

static void CreateBioseqSetMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

static ForM CreateBioseqSetEditForm (BioseqSetPtr bssp, Uint2 entityID)

{
  ButtoN            b;
  BioseqSetFormPtr  bsfp;
  GrouP             c;
  GrouP             h;
  WindoW            w;

  bsfp = (BioseqSetFormPtr) MemNew (sizeof (BioseqSetForm));
  if (bsfp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Bioseq Set", NULL);
    SetObjectExtra (w, bsfp, StdCleanupFormProc);
    bsfp->form = (ForM) w;
    bsfp->formmessage = CreateBioseqSetMessageProc;

    bsfp->bssp = bssp;
    bsfp->entityID = entityID;

    h = HiddenGroup (w, -2, 0, NULL);
    StaticPrompt (h, "Class", 0, popupMenuHeight, programFont, 'l');
    bsfp->class_control = PopupList (h, TRUE, NULL);
    SetObjectExtra (bsfp->class_control, bsfp, NULL);
    InitEnumPopup (bsfp->class_control, bioseqset_class_alist, NULL);
    if (bssp != NULL) {
      SetEnumPopup (bsfp->class_control, bioseqset_class_alist, (UIEnum) bssp->_class);
    }

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", AcceptBioseqSetEditProc);
    SetObjectExtra (b, bsfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK BioseqSetEditFunc (Pointer data)

{
  BioseqSetPtr      bssp;
  OMProcControlPtr  ompcp;
  ForM              w;

  ompcp = (OMProcControlPtr) data;
  bssp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQSET :
      bssp = (BioseqSetPtr) ompcp->input_data;
      break;
   case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  /* if (bssp == NULL) return OM_MSG_RET_ERROR; */

  w = CreateBioseqSetEditForm (bssp, ompcp->input_entityID);
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}

/*#ifdef EXTRA_SERVICES*/

static int LIBCALLBACK SortVnpByClickableItemDesc (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr cip1, cip2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      cip1 = (ClickableItemPtr) vnp1->data.ptrvalue;
      cip2 = (ClickableItemPtr) vnp2->data.ptrvalue;
      if (cip1 != NULL && cip2 != NULL) {
        str1 = cip1->description;
        str2 = cip2->description;
        if (str1 != NULL && str2 != NULL) {
          return StringCmp (str1, str2);
        }
      }
    }
  }
  return 0;
}

static int LIBCALLBACK SortVnpByChoiceAndClickableItemDesc (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr cip1, cip2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      cip1 = (ClickableItemPtr) vnp1->data.ptrvalue;
      cip2 = (ClickableItemPtr) vnp2->data.ptrvalue;
      if (cip1 != NULL && cip2 != NULL) {
        str1 = cip1->description;
        str2 = cip2->description;
        if (str1 != NULL && str2 != NULL) {
          if (vnp1->choice > vnp2->choice) {
            return 1;
          } else if (vnp1->choice < vnp2->choice) {
            return -1;
          }
          return StringCmp (str1, str2);
        }
      }
    }
  }
  return 0;
}

static ValNodePtr UniqueAndCountValNodeCS (ValNodePtr list)

{
  Int4          count;
  ValNodePtr    curr;
  size_t        len;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       tmp;
  ValNodePtr    vnp;
  ClickableItemPtr cip, cip_last = NULL;

  if (list == NULL) return NULL;
  cip_last = (ClickableItemPtr) list->data.ptrvalue;
  vnp = list->next;
  if (vnp == NULL) return list;
  prev = (Pointer PNTR) &(list->next);
  count = 1;
  curr = list;
  while (vnp != NULL) {
    next = vnp->next;
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (StringCmp (cip_last->description, cip->description) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeLink (&(cip_last->item_list), cip->item_list);
      cip->item_list = NULL;
      ValNodeLink (&(cip_last->subcategories), cip->subcategories);
      cip->subcategories = NULL;
      FreeClickableList (vnp); 
      count++;
    } else {
      len = StringLen (cip_last->description) + 20;
      tmp = (CharPtr) MemNew (len);
      if (tmp != NULL) {
        sprintf (tmp, "%6ld   %s", (long) count, cip_last->description);
        cip_last->description = MemFree (cip_last->description);
        cip_last->description = tmp;
        cip_last->subcategories = ValNodeSort (cip_last->subcategories, SortVnpByClickableItemDesc);
        cip_last->subcategories = UniqueAndCountValNodeCS (cip_last->subcategories);
      }
      cip_last = cip;
      prev = (Pointer PNTR) &(vnp->next);
      count = 1;
      curr = vnp;
    }
    vnp = next;
  }
  len = StringLen (cip_last->description) + 20;
  tmp = (CharPtr) MemNew (len);
  if (tmp != NULL) {
    sprintf (tmp, "%6ld   %s", (long) count, cip_last->description);
    cip_last->description = MemFree (cip_last->description);
    cip_last->description = tmp;
    cip_last->subcategories = ValNodeSort (cip_last->subcategories, SortVnpByClickableItemDesc);
    cip_last->subcategories = UniqueAndCountValNodeCS (cip_last->subcategories);
  }

  return list;
}

static ValNodePtr SortFlatFile (ValNodePtr head, Boolean reverse, Boolean byblock)

{
  ValNodePtr  next;
  ValNodePtr  tail = NULL;
  ValNodePtr  vnp;

  if (head == NULL) return NULL;

  if (byblock) {
    head = ValNodeSort (head, SortVnpByChoiceAndClickableItemDesc);
  } else {
    head = ValNodeSort (head, SortVnpByClickableItemDesc);
  }
  if (reverse) {
    for (vnp = head; vnp != NULL; vnp = next) {
      next = vnp->next;
      vnp->next = tail;
      tail = vnp;
    }
    head = tail;
  } else {
    head = UniqueAndCountValNodeCS (head);
  }

  return head;
}

static CharPtr StartsWithQualOrFeat (CharPtr str)
{
  Int4 space_before_qual, qual_len, space_after_feat;
  CharPtr qual_name = NULL;

  if (StringHasNoText (str)) return NULL;

  space_before_qual = StringSpn (str, " \t");
  /* qual is name between slash and equals sign */
  if (str[space_before_qual] == '/') {
    qual_len = StringCSpn (str + space_before_qual, "=");
    if (qual_len != 0 && qual_len != StringLen (str + space_before_qual)) {
      qual_name = (CharPtr) MemNew ((qual_len + 1) * sizeof(Char));
      StringNCpy (qual_name, str + space_before_qual, qual_len);
      qual_name[qual_len] = 0;
    }
  } else {
    qual_len = StringCSpn (str + space_before_qual, " \t");
    space_after_feat = StringSpn (str + space_before_qual + qual_len, " \t");
    /* look for location  after feature name */
    if (space_before_qual == 5
        && qual_len + space_after_feat == 16
        && (isdigit(str[space_before_qual + qual_len + space_after_feat])
        || StringNCmp (str + space_before_qual + qual_len + space_after_feat, "complement", 10) == 0)) {
      qual_name = (CharPtr) MemNew ((qual_len + 1) * sizeof(Char));
      StringNCpy (qual_name, str + space_before_qual, qual_len);
      qual_name[qual_len] = 0;
    }
  }

  return qual_name;
}


static void CaptureFFLineEx (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  Boolean include_sequence,
  Boolean byqual
)

{
  Char             ch;
  CharPtr          copy;
  ValNodePtr PNTR  head;
  CharPtr          ptr;
  CharPtr          tmp;
  ValNodePtr       vnp;
  ClickableItemPtr cip, subcip;
  ValNodePtr       item_list = NULL;
  BioseqPtr        bsp;
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  CharPtr           qual_name;

  if (!include_sequence && blocktype == SEQUENCE_BLOCK) return;

  head = (ValNodePtr PNTR) userdata;
  copy = StringSaveNoNull (str);
  if (copy == NULL) return;

  ptr = copy;
  tmp = StringChr (ptr, '\n');
  while (tmp != NULL) {
    ch = *tmp;
    *tmp = '\0';
    if (!StringHasNoText (ptr)) {
      item_list = NULL;
      if (itemtype == OBJ_BIOSEQ) {
        bsp =  GetBioseqGivenIDs (entityID, itemID, itemtype);
        if (bsp != NULL) {
          item_list = ValNodeNew (NULL);
          item_list->choice = OBJ_BIOSEQ;
          item_list->data.ptrvalue = bsp;
        }
      } else if (itemtype == OBJ_SEQFEAT) {
        sfp = SeqMgrGetDesiredFeature (entityID, NULL, itemID, 0, NULL, &fcontext);
        if (sfp != NULL) {
          item_list = ValNodeNew (NULL);
          item_list->choice = OBJ_SEQFEAT;
          item_list->data.ptrvalue = sfp;
        }
      } else if (itemtype == OBJ_SEQDESC) {
        sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, itemID, 0, NULL, &dcontext);
        item_list = ValNodeNew (NULL);
        item_list->choice = OBJ_SEQDESC;
        item_list->data.ptrvalue = sdp;
      }

      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->clickable_item_type = blocktype;

      if (byqual && (qual_name = StartsWithQualOrFeat (ptr)) != NULL) {
        cip->description = qual_name;
        subcip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (subcip, 0, sizeof (ClickableItemData));
        subcip->clickable_item_type = blocktype;
        subcip->description = StringSave (ptr);
        if (item_list != NULL) {
          ValNodeAddPointer (&(subcip->item_list), item_list->choice, item_list->data.ptrvalue);
        }
        ValNodeAddPointer (&(cip->subcategories), 0, subcip);
        /* iterate to add the rest of the lines of the qual */
        while (ch != 0 && tmp != NULL && (qual_name = StartsWithQualOrFeat (tmp + 1)) == NULL) {
          *tmp = ch;
          tmp++;
          ptr = tmp;
          tmp = StringChr (ptr, '\n');
          if (tmp != NULL) {
            ch = *tmp;
            *tmp = '\0';
            subcip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
            MemSet (subcip, 0, sizeof (ClickableItemData));
            subcip->clickable_item_type = blocktype;
            subcip->description = StringSave (ptr);
            if (item_list != NULL) {
              ValNodeAddPointer (&(subcip->item_list), item_list->choice, item_list->data.ptrvalue);
            }
            ValNodeAddPointer (&(cip->subcategories), 0, subcip);
          }
        }
        qual_name = MemFree (qual_name);
      } else {
        cip->description = StringSave (ptr);
      }
      cip->item_list = item_list;     
      vnp = ValNodeNew(NULL);
      vnp->choice = blocktype;
      vnp->data.ptrvalue = cip;
      if (*head == NULL) {
        *head = vnp;
      } else {
        vnp->next = *head;
        *head = vnp;
      }
    }
    /* tmp may have become NULL while processing quals */
    if (tmp != NULL) {
      *tmp = ch;
      tmp++;
      ptr = tmp;
      tmp = StringChr (ptr, '\n');
    }
  }

  MemFree (copy);
}

static void CaptureFFLine (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID
)

{
  CaptureFFLineEx (str, userdata, blocktype, entityID, itemtype, itemID, TRUE, FALSE);
}

static void CaptureFFLineNoSequence (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID
)

{
  CaptureFFLineEx (str, userdata, blocktype, entityID, itemtype, itemID, FALSE, FALSE);
}

static void CaptureFFLineNoSequenceByQual (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID
)

{
  CaptureFFLineEx (str, userdata, blocktype, entityID, itemtype, itemID, FALSE, TRUE);
}

static void CaptureFFLineByQual (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID
)

{
  CaptureFFLineEx (str, userdata, blocktype, entityID, itemtype, itemID, TRUE, TRUE);
}

extern ValNodePtr GetSUCCommonList (SeqEntryPtr sep, Boolean reverse, Boolean byblock, Boolean showsequence, Boolean byqual)
{
  XtraBlock       xtra;
  ValNodePtr      head = NULL;
  ErrSev          level;
  Boolean         okay;
  SeqEntryPtr     oldscope;
  Uint2           entityID;

  if (sep == NULL) return NULL;
  
  MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
  if (showsequence)
  {
    if (byqual) 
    {
      xtra.ffwrite = CaptureFFLineByQual;
    } 
    else
    {
      xtra.ffwrite = CaptureFFLine;  	
    }
  }
  else
  {
    if (byqual)
    {
      xtra.ffwrite = CaptureFFLineNoSequenceByQual;
    }
    else 
    {
  	  xtra.ffwrite = CaptureFFLineNoSequence;
    }
  }
  xtra.userdata = (Pointer) &head;
  level = ErrSetMessageLevel (SEV_MAX);
  oldscope = SeqEntrySetScope (sep);
  okay = SeqEntryToGnbk (sep, NULL, GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE,
                         SHOW_CONTIG_FEATURES, 0, 0, &xtra, NULL);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  SeqEntrySetScope (oldscope);
  ErrSetMessageLevel (level);
  if (okay) {
    head = SortFlatFile (head, FALSE, byblock);
    if (reverse) {
      head = SortFlatFile (head, TRUE, FALSE);
    }
  }
  return head;
}

static CharPtr block_names[] = {
  "Unknown Block Type",
  "HEAD",
  "LOCUS",
  "DEFLINE",
  "ACCESSION",
  "VERSION",
  "PROJECT",
  "PID",
  "DBSOURCE",
  "DATE",
  "KEYWORDS",
  "SEGMENT",
  "SOURCE",
  "ORGANISM",
  "REF_STATS",
  "REFERENCE",
  "PRIMARY",
  "COMMENT",
  "FEAT_STATS",
  "FEATHEADER",
  "SOURCEFEAT",
  "FEATURE",
  "BASECOUNT",
  "ORIGIN",
  "SEQUENCE",
  "CONTIG",
  "WGS",
  "GENOME",
  "SLASH",
  "TAIL"
};

extern ValNodePtr CategorizeSUCBlocks (ValNodePtr head)
{
  ValNodePtr vnp, next_vnp, category_list = NULL, new_list = NULL;
  ClickableItemPtr cip, first_cip;
  
  vnp = head;
  category_list = head;
  while (vnp != NULL && vnp->next != NULL) {
    next_vnp = vnp->next;
    if (vnp->choice != next_vnp->choice) {
      if (vnp->choice == FEATHEADER_BLOCK
          || vnp->choice == SLASH_BLOCK) {
        /* remove these blocktypes */
        vnp->next = NULL;
        category_list = FreeClickableList (category_list);
      } else {
        /* create a category to hold this blocktypes */
        first_cip = category_list->data.ptrvalue;
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (cip, 0, sizeof (ClickableItemData));
        cip->clickable_item_type = first_cip->clickable_item_type;
        cip->description = StringSave (block_names[first_cip->clickable_item_type]);
        vnp->next = NULL;
        cip->subcategories = category_list;
        ValNodeAddPointer (&new_list, vnp->choice, cip);
      }
      category_list = next_vnp;
    }
    vnp = next_vnp;
  }      
  if (category_list != NULL) {
    if (category_list->choice == FEATHEADER_BLOCK
        || category_list->choice == SLASH_BLOCK) {
      category_list = FreeClickableList (category_list);
    } else {
      first_cip = category_list->data.ptrvalue;
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->clickable_item_type = first_cip->clickable_item_type;
      cip->description = StringSave (block_names[cip->clickable_item_type]);
      vnp->next = NULL;
      cip->subcategories = category_list;
      category_list = NULL;
      ValNodeAddPointer (&new_list, vnp->choice, cip);
    }
  }
  return new_list;
}


typedef struct resucoptions
{
  Boolean     reverse;
  Boolean     byblock;
  Boolean     showsequence;
  Uint2       entityID;
} ReSUCOptionsData, PNTR ReSUCOptionsPtr;

static Pointer ReSUCFree (Pointer userdata)
{
  return MemFree (userdata);  
}

static CharPtr ReSUCCallback (Pointer userdata)
{
  ReSUCOptionsPtr rp;
  SeqEntryPtr     sep;
  FILE            *fp;
  ValNodePtr      head = NULL;
  Char            path [PATH_MAX];
  CharPtr         str;
  ValNodePtr      vnp;
  ClickableItemPtr cip;
  
  if (userdata == NULL)
  {
    return NULL;
  }
  
  rp = (ReSUCOptionsPtr) userdata;
  
  sep = GetTopSeqEntryForEntityID (rp->entityID);
  if (sep == NULL) return NULL;
  
  head = GetSUCCommonList (sep, rp->reverse, rp->byblock, rp->showsequence, FALSE);

  if (head != NULL) {  
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        cip = (ClickableItemPtr) vnp->data.ptrvalue;        
        str = cip->description;
        fprintf (fp, "%s\n", str);
      }
      FileClose (fp);
      return StringSave (path);
    }
  }
  return NULL;
}

static void SUCCommonProcEx (IteM i, Boolean reverse, Boolean byblock, 
                             Boolean showsequence, ButtoN b, Boolean use_new_suc)

{
  BaseFormPtr  bfp;
  FILE         *fp;
  ValNodePtr   head = NULL;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  CharPtr      str;
  ValNodePtr   vnp;
  ReSUCOptionsPtr rp = NULL;
  ClickableItemPtr cip;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  head = GetSUCCommonList (sep, reverse, byblock, showsequence, use_new_suc);
  
  if (head != NULL) {  
    if (use_new_suc) {    
      NewSUC (head, bfp->input_entityID, reverse, byblock, showsequence);
    } else {
      TmpNam (path);
      fp = FileOpen (path, "w");
      if (fp != NULL) {
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          cip = (ClickableItemPtr) vnp->data.ptrvalue;        
          str = cip->description;
          fprintf (fp, "%s\n", str);
        }
        FileClose (fp);
      }

      rp = (ReSUCOptionsPtr) MemNew (sizeof (ReSUCOptionsData));
      if (rp != NULL)
      {
        rp->reverse = reverse;
        rp->byblock = byblock;
        rp->showsequence = showsequence;
        rp->entityID = bfp->input_entityID;
        LaunchGeneralTextViewerWithRepopulate (path, "Sort Unique Count", 
                                               ReSUCCallback, rp, ReSUCFree);
      }
      else
      {
        LaunchGeneralTextViewer (path, "Sort Unique Count");
      }
    
      FileRemove (path);
      ValNodeFreeData (head);
    }    
  }
  ArrowCursor ();
  Update ();
}

static void SUCCommonProc (IteM i, Boolean reverse, Boolean byblock, 
                           Boolean showsequence, ButtoN b)
{
  SUCCommonProcEx (i, reverse, byblock, showsequence, b, FALSE);
}

static void SUCProc (IteM i)

{
  SUCCommonProc (i, FALSE, FALSE, TRUE, NULL);
}

static void SUCRProc (IteM i)

{
  SUCCommonProc (i, TRUE, FALSE, TRUE, NULL);
}

static void SUCBProc (IteM i)

{
  SUCCommonProc (i, FALSE, TRUE, TRUE, NULL);
}

static void SUCBNoSequenceProc (IteM i)

{
  SUCCommonProc (i, FALSE, TRUE, FALSE, NULL);
}

extern void SUCSubmitterProc (IteM i)

{
  SUCCommonProc (i, FALSE, TRUE, TRUE, NULL);
}


static void SUCClickableProc (IteM i)

{
  SUCCommonProcEx (i, FALSE, TRUE, TRUE, NULL, TRUE);
}

static void SUCClickableNoSequenceProc (IteM i)

{
  SUCCommonProcEx (i, FALSE, TRUE, FALSE, NULL, TRUE);
}

/*#endif*/

#if 0
#ifdef OS_UNIX
extern void SUCCommonProc (IteM i, Boolean reverse, ButtoN b);
extern void SUCCommonProc (IteM i, Boolean reverse, ButtoN b)

{
  BaseFormPtr  bfp;
  Char         cmmd [256];
  FILE         *fp;
  ErrSev       level;
  Boolean      okay;
  SeqEntryPtr  oldscope;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return;
  WatchCursor ();
  Update ();
  level = ErrSetMessageLevel (SEV_MAX);
  oldscope = SeqEntrySetScope (sep);
  okay = SeqEntryToGnbk (sep, NULL, GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  SeqEntrySetScope (oldscope);
  ErrSetMessageLevel (level);
  FileClose (fp);
  if (okay) {
    if (reverse) {
      sprintf (cmmd, "sort %s | uniq -c | sort -nr > %s.suc; rm %s", path, path, path);
    } else {
      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path, path, path);
    }
    system (cmmd); /* removes original flat file */
    StringCat (path, ".suc");
    LaunchGeneralTextViewer (path, "Sort Unique Count");
    FileRemove (path); /* removes sorted/uniqued/counted file */
  } else {
    FileRemove (path);
  }
  ArrowCursor ();
  Update ();
}

static void SUCProc (IteM i)

{
  SUCCommonProc (i, FALSE, NULL);
}

static void SUCRProc (IteM i)

{
  SUCCommonProc (i, TRUE, NULL);
}
#endif
#endif

/*#ifdef INTERNAL_NCBI_SEQUIN*/
/*#ifdef NEW_TAXON_SERVICE*/
/*
*  Get top 900 organisms, remove all ending with 'sp.',
*  all starting with unknown, unidentified, uncultured,
*  unclassified, cloning vector, synthetic construct,
*  mixed EST library, also remove Rattus rattus,
*  trim to remaining 800 organisms, remove any counts
*
*  Process output in UNIX with
*    'sort -f +0 sorttax | uniq > taxlist.txt'
*  and
*    'sort -n +0 sortlin | uniq > lineages.txt'
*
*  Then add version number to top of each file
*/

static Boolean WriteTaxNode (CharPtr sci, CharPtr com, CharPtr syn,
                             Int4 gc, Int4 mgc, CharPtr div,
                             Int4 taxID, FILE *fout)

{
  Char  str [256];

  if (fout == NULL) return FALSE;
  if (sci != NULL && sci [0] != '\0') {
    if (div == NULL || div [0] == '\0') div = "???";
    if (com != NULL && com [0] != '\0') {
      sprintf (str, "%s\t%s\t%ld\t%ld\t%s\t%ld\n", sci, com,
               (long) gc, (long) mgc, div, (long) taxID);
      fprintf (fout, "%s", str);
    } else {
      sprintf (str, "%s\t\t%ld\t%ld\t%s\t%ld\n", sci,
               (long) gc, (long) mgc, div, (long) taxID);
      fprintf (fout, "%s", str);
    }
    return TRUE;
  }
  return FALSE;
}

static void WriteLineage (Int4 taxID, CharPtr lineage, FILE *fout)

{
  if (fout == NULL || lineage == NULL) return;
  fprintf (fout, "%ld\t%s\n", (long) taxID, lineage);
}

static void ProcessTaxNode (OrgRefPtr orp, FILE *fout1, FILE *fout2)

{
  Int4           gc;
  Int4           mgc;
  OrgModPtr      omp;
  OrgNamePtr     onp;
  Int4           taxID = -1;
  DbtagPtr       d;
  ValNodePtr     vnp;
  CharPtr        div = NULL;

  if (orp == NULL || fout1 == NULL || fout2 == NULL) return;
  
  for (vnp = orp->db; vnp != NULL; vnp = vnp->next) 
  {
    d = (DbtagPtr) vnp->data.ptrvalue;
    if (StringCmp(d->db, "taxon") == 0) 
    {
      taxID = d->tag->id;
      break;
    }
  }
  if (taxID > -1)
  {
    gc = 0;
    mgc = 0;
    onp = orp->orgname;
    if (onp != NULL) 
    {
      gc = onp->gcode;
      mgc = onp->mgcode;
      div = onp->div;
    }
    if (WriteTaxNode (orp->taxname, orp->common, NULL, gc, mgc, div, taxID, fout1))
    {
      if (onp != NULL) 
      {
        WriteLineage (taxID, onp->lineage, fout2);
        for (omp = onp->mod; omp != NULL; omp = omp->next) 
        {
          if (omp->subtype == ORGMOD_gb_anamorph ||
            omp->subtype == ORGMOD_gb_synonym) 
          {
            WriteTaxNode (omp->subname, orp->common, NULL, gc, mgc, div, taxID, fout1);
          }
        }
      }
    }
  }
}

static void PrepareTaxListProc (IteM i)

{
  Char     ch;
  FILE     *fin;
  FILE     *fout1;
  FILE     *fout2;
  Boolean  goOn;
  Char     orgname [256];
  CharPtr  ptr;
  ValNodePtr org_list = NULL;
  ValNodePtr response_list = NULL;
  ValNodePtr vnp;
  OrgRefPtr  orp;

  orgname [0] = '\0';
  goOn = TRUE;
  Message (MSG_POST, "Finding File");
  fin = FileOpen ("orglist", "r");
  if (fin != NULL) {
    fout1 = FileOpen ("sorttax", "w");
    if (fout1 != NULL) {
      fout2 = FileOpen ("sortlin", "w");
      if (fout2 != NULL) {
        Message (MSG_POST, "Processing");
        while (goOn) {
          goOn = (FileGets (orgname, sizeof (orgname), fin) != NULL);
          if (goOn) {
            ptr = orgname;
            ch = *ptr;
            while (ch != '\n' && ch != '\r' && ch != '\0') {
              ptr++;
              ch = *ptr;
            }
            *ptr = '\0';
            Message (MSG_POST, "Organism '%s'", orgname);
            ValNodeCopyStr (&org_list, 2, orgname);
          }
        }
        Message (MSG_POST, "Retrieving from server");
        response_list =  Taxon3GetOrgRefList (org_list);
        for (vnp = response_list; vnp != NULL; vnp = vnp->next)
        {
          orp = (OrgRefPtr) vnp->data.ptrvalue;
          ProcessTaxNode (orp, fout1, fout2);
          vnp->data.ptrvalue = OrgRefFree (orp);
        }
        response_list = ValNodeFree (response_list);
        org_list = ValNodeFreeData (org_list);

        FileClose (fout2);
        
      }
      FileClose (fout1);
    }
    FileClose (fin);
  } else {
    Message (MSG_OK, "Could not find orglist file");
  }
  
  Message (MSG_POST, "Finished");
}
/*#endif*/
/*#endif*/


#define REMOVE_PUB   1

static ENUM_ALIST(pub_field_alist)
  {" ",                    0},
  {"Remark",               1},
END_ENUM_ALIST

typedef struct pubformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  PopuP          fromfield;
  Int2           fromval;
  ButtoN         removeIncompletePubs;
} PubFormData, PNTR PubFormPtr;

static Boolean ProcessEachPubFunc (GatherContextPtr gcp)

{
  PubdescPtr  pdp;
  PubFormPtr  pfp;
  ValNodePtr  sdp;
  SeqFeatPtr  sfp;

  if (gcp == NULL) return TRUE;
  pfp = (PubFormPtr) gcp->userdata;
  if (pfp == NULL) return TRUE;
  pdp = NULL;
  if (gcp->thistype == OBJ_SEQDESC) {
    sdp = (ValNodePtr) gcp->thisitem;
    if (sdp != NULL && sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
    }
  } else if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    }
  }
  if (pdp == NULL) return TRUE;
  switch (pfp->fromval) {
    case 1 :
      pdp->comment = MemFree (pdp->comment);
      break;
    default :
      break;
  }
  return TRUE;
}

/* PubdescIsIncomplete modified from ValidatePubdesc in valid.c */
static Boolean HasNoText (CharPtr str)

{
  Char  ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static Boolean HasNoName (ValNodePtr name)

{
	AuthorPtr  ap;
	NameStdPtr  nsp;
	PersonIdPtr  pid;

	if (name != NULL) {
		ap = name->data.ptrvalue;
		if (ap != NULL) {
			pid = ap->name;
			if (pid != NULL) {
				if (pid->choice == 2) {
					nsp = pid->data;
					if (nsp != NULL) {
						if (! HasNoText (nsp->names [0])) {
							return FALSE;
						}
					}
				}
			}
		}
	}
	return TRUE;
}

static Boolean PubdescIsIncomplete (PubdescPtr pdp)

{
	AuthListPtr  alp;
	CitArtPtr  cap;
	Boolean  hasName, hasTitle;
	ValNodePtr  name;
	ValNodePtr  title;
	ValNodePtr  vnp;

	if (pdp == NULL) return TRUE;
	for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
		switch (vnp->choice) {
			case PUB_Article :
				cap = (CitArtPtr) vnp->data.ptrvalue;
				hasName = FALSE;
				hasTitle = FALSE;
				if (cap != NULL) {
					for (title = cap->title; title != NULL; title = title->next) {
						if (! HasNoText ((CharPtr) title->data.ptrvalue)) {
							hasTitle = TRUE;
						}
					}
					if (! hasTitle) {
						return TRUE;
					}
					alp = cap->authors;
					if (alp != NULL) {
						if (alp->choice == 1) {
							for (name = alp->names; name != NULL; name = name->next) {
								if (! HasNoName (name)) {
									hasName = TRUE;
								}
							}
						} else if (alp->choice == 2 || alp->choice == 3) {
							for (name = alp->names; name != NULL; name = name->next) {
								if (! HasNoText ((CharPtr) name->data.ptrvalue)) {
									hasName = TRUE;
								}
							}
						}
					}
					if (! hasName) {
						return TRUE;
					}
				}
				break;
			default :
				break;
		}
	}
	return FALSE;
}

static void RemoveIncompletePubs (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Boolean       empty;
  SeqAnnotPtr   nextsap;
  ValNodePtr    nextsdp;
  SeqFeatPtr    nextsfp;
  PubdescPtr    pdp;
  Pointer PNTR  prevsap;
  Pointer PNTR  prevsdp;
  Pointer PNTR  prevsfp;
  SeqAnnotPtr   sap;
  ValNodePtr    sdp;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        empty = FALSE;
        if (sfp->data.choice == SEQFEAT_PUB && sfp->data.value.ptrvalue != NULL) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          empty = PubdescIsIncomplete (pdp);
        }
        if (empty) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          SeqFeatFree (sfp);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
  while (sdp != NULL) {
    nextsdp = sdp->next;
    empty = FALSE;
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      empty = PubdescIsIncomplete (pdp);
    }
    if (empty) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static void DoProcessPub (ButtoN b)

{
  PubFormPtr   pfp;
  GatherScope  gs;
  SeqEntryPtr  sep;
  UIEnum       val;

  pfp = (PubFormPtr) GetObjectExtra (b);
  if (pfp == NULL || pfp->input_entityID == 0) return;
  Hide (pfp->form);
  WatchCursor ();
  Update ();
  pfp->fromval = 0;
  if (GetEnumPopup (pfp->fromfield, pub_field_alist, &val)) {
    pfp->fromval = (Int2) val;
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  GatherEntity (pfp->input_entityID, (Pointer) pfp, ProcessEachPubFunc, &gs);
  if (GetStatus (pfp->removeIncompletePubs)) {
    sep = GetTopSeqEntryForEntityID (pfp->input_entityID);
    if (sep != NULL) {
      SeqEntryExplore (sep, NULL, RemoveIncompletePubs);
    }
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (pfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, pfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (pfp->form);
}

static void PubMessageProc (ForM f, Int2 mssg)

{
  PubFormPtr  pfp;

  pfp = (PubFormPtr) GetObjectExtra (f);
  if (pfp != NULL) {
    if (pfp->appmessage != NULL) {
      pfp->appmessage (f, mssg);
    }
  }
}

static void ProcessPub (IteM i, Int2 type)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  PubFormPtr         pfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  pfp = (PubFormPtr) MemNew (sizeof (PubFormData));
  if (pfp == NULL) return;
  pfp->type = type;
  switch (type) {
    case REMOVE_PUB :
      title = "Pub Removal";
      break;
    default :
      title = "?";
      break;
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, pfp, StdCleanupFormProc);
  pfp->form = (ForM) w;
  pfp->formmessage = PubMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    pfp->appmessage = sepp->handleMessages;
  }

  pfp->input_entityID = bfp->input_entityID;
  pfp->input_itemID = bfp->input_itemID;
  pfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 4, 0, NULL);
  switch (type) {
    case REMOVE_PUB :
      StaticPrompt (g, "Remove", 0, popupMenuHeight, programFont, 'l');
      pfp->fromfield = PopupList (g, TRUE, NULL);
      SetObjectExtra (pfp->fromfield, pfp, NULL);
      InitEnumPopup (pfp->fromfield, pub_field_alist, NULL);
      SetEnumPopup (pfp->fromfield, pub_field_alist, 0);
      break;
    default :
      break;
  }

  pfp->removeIncompletePubs = NULL;
  if (type == REMOVE_PUB) {
    pfp->removeIncompletePubs = CheckBox (h, "Remove Incomplete Publications", NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoProcessPub);
  SetObjectExtra (b, pfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) pfp->removeIncompletePubs, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static void RemovePub (IteM i)

{
  ProcessPub (i, REMOVE_PUB);
}

static void SequinSetPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  AuthListPtr    alp;
  Char           buf [2048];
  CitGenPtr      cgp;
  CharPtr        consortium, str, tmp;
  ValNodePtr     vnp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    alp = GetAuthListPtr (pdp, NULL);
    if (alp != NULL) {
      consortium = NULL;
      str = GetAuthorsString (GENBANK_FMT, alp, &consortium, NULL, NULL);
      if (StringDoesHaveText (buf) || StringDoesHaveText (str) || StringDoesHaveText (consortium)) {
        tmp = MemNew (StringLen (buf) + StringLen (str) + StringLen (consortium) + 10);
        if (tmp != NULL) {
          StringCpy (tmp, buf);
          if (StringDoesHaveText (str)) {
            StringCat (tmp, "; ");
            StringCat (tmp, str);
          }
          if (StringDoesHaveText (consortium)) {
            StringCat (tmp, "; ");
            StringCat (tmp, consortium);
          }
          ovp->idx.scratch = tmp;
        }
      }
      MemFree (str);
      MemFree (consortium);
    }
  }
}

static void SequinClearPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  ovp->idx.scratch = MemFree (ovp->idx.scratch);
}

static void MarkCollidingPubs (BioseqPtr bsp, Pointer userdata)

{
  SeqMgrDescContext  dcontext;
  ValNodePtr         head = NULL, last = NULL, vnp, nxt;
  ObjValNodePtr      ovp, ovp1, ovp2;
  PubdescPtr         pdp;
  SeqDescrPtr        sdp;
  CharPtr            str1, str2;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  if (SeqMgrGetParentOfPart (bsp, NULL) != NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.scratch != NULL) {
        vnp = ValNodeAddPointer (&last, 0, (Pointer) ovp);
        if (head == NULL) {
          head = vnp;
        }
        last = vnp;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    ovp1 = (ObjValNodePtr) vnp->data.ptrvalue;
    if (ovp1 == NULL) continue;
    str1 = (CharPtr) ovp1->idx.scratch;
    if (StringHasNoText (str1)) continue;
    for (nxt = vnp->next; nxt != NULL; nxt = nxt->next) {
      ovp2 = (ObjValNodePtr) nxt->data.ptrvalue;
      if (ovp2 == NULL) continue;
      str2 = (CharPtr) ovp2->idx.scratch;
      if (StringHasNoText (str2)) continue;
      if (StringICmp (str1, str2) == 0) {
        ovp1->idx.deleteme = TRUE;
      }
    }
  }

  ValNodeFree (head);
}

static void RemoveCollidingPublications (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitDescriptorsInSep (sep, NULL, SequinSetPubScratchData);
  VisitBioseqsInSep (sep, NULL, MarkCollidingPubs);
  VisitDescriptorsInSep (sep, NULL, SequinClearPubScratchData);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void StripSerials (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  EntryStripSerialNumber (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveLocalIDsProc (SeqEntryPtr sep, Boolean do_nuc, Boolean do_prt)

{
  BioseqPtr     bsp;
  SeqIdPtr      nextsip;
  Pointer PNTR  prevsip;
  Boolean       replaced;
  SeqIdPtr      sip;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (ISA_na (bsp->mol) && (! do_nuc)) return;
  if (ISA_aa (bsp->mol) && (! do_prt)) return;
  sip = bsp->id;
  prevsip = (Pointer PNTR) &(bsp->id);
  replaced = FALSE;
  while (sip != NULL) {
    nextsip = sip->next;
    if (sip->choice == SEQID_LOCAL) {
      (*prevsip) = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
      replaced = TRUE;
    } else {
      prevsip = (Pointer PNTR) &(sip->next);
    }
    sip = nextsip;
  }
  if (replaced) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }
}

static void RemoveAllLocalIDsProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  RemoveLocalIDsProc (sep, TRUE, TRUE);
}

static void RemoveNucLocalIDsProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  RemoveLocalIDsProc (sep, TRUE, FALSE);
}

static void RemovePrtLocalIDsProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  RemoveLocalIDsProc (sep, FALSE, TRUE);
}

static void RemoveLocalIDs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveAllLocalIDsProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveNucLocalIDs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveNucLocalIDsProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemovePrtLocalIDs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemovePrtLocalIDsProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveLocusProc (BioseqPtr bsp, Pointer userdata)

{
  Boolean       reindex = FALSE;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_OTHER :
      case SEQID_DDBJ :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL) {
          tsip->name = MemFree (tsip->name);
          reindex = TRUE;
        }
        break;
      default :
        break;
    }
  }
  if (reindex) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }
}

static void RemoveLocusFromParts (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (Message (MSG_OKC, "Are you SURE you want to do this?") == ANS_CANCEL) return;
  VisitSequencesInSep (sep, NULL, VISIT_PARTS, RemoveLocusProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct gosstruc {
  CharPtr       term;
  Char          goid [32];
  CharPtr       evidence;
  Int4          pmid;
  UserFieldPtr  ufp;
} GosStruc, PNTR GosStrucPtr;

static int LIBCALLBACK SortVnpByGssp (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GosStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GosStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GosStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}

static UserFieldPtr SortGoTerms (
  UserFieldPtr entryhead
)

{
  UserFieldPtr  entry, topufp, ufp, lastufp;
  CharPtr       evidence, goid, textstr;
  Char          gid [32];
  GosStrucPtr   gsp, lastgsp;
  ValNodePtr    head = NULL, vnp;
  Int2          j;
  ObjectIdPtr   oip;
  Int4          pmid;

  if (entryhead == NULL) return entryhead;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    textstr = NULL;
    evidence = NULL;
    goid = NULL;
    pmid = 0;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; goFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, goFieldType [j]) == 0) break;
      }
      if (goFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }

    if (StringDoesHaveText (textstr)) {
      gsp = (GosStrucPtr) MemNew (sizeof (GosStruc));
      if (gsp != NULL) {
        gsp->term = textstr;
        StringNCpy_0 (gsp->goid, goid, sizeof (gsp->goid));
        gsp->evidence = evidence;
        gsp->pmid = pmid;
        gsp->ufp = entry;
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL) return entryhead;
  head = ValNodeSort (head, SortVnpByGssp);

  entryhead = NULL;
  lastgsp = NULL;
  lastufp = NULL;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GosStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL || gsp->ufp == NULL) continue;
    if (lastgsp != NULL &&
        (StringICmp (gsp->term, lastgsp->term) == 0 || StringICmp (gsp->goid, lastgsp->goid) == 0) &&
        (gsp->pmid == lastgsp->pmid && StringICmp (gsp->evidence, lastgsp->evidence) == 0)) {
      gsp->ufp->next = NULL;
      UserFieldFree (gsp->ufp);
    } else {
      if (lastufp != NULL) {
        lastufp->next = gsp->ufp;
      } else {
        entryhead = gsp->ufp;
      }
      lastufp = gsp->ufp;
      lastufp->next = NULL;
    }
    lastgsp = gsp;
  }

  ValNodeFreeData (head);

  return entryhead;
}

static void SortGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
 
  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; goQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, goQualType [i]) == 0) break;
  }
  if (goQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  ufp->data.ptrvalue = SortGoTerms (entry);
}

static void SortGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, SortGoTermsUfp);
  }
}

static void LookForGoTerms (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp != NULL && sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, NULL, SortGoTermsSfp);
  }
}

static void RemoveDuplicateGeneOntology (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, LookForGoTerms);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void LookForDupPCR (BioSourcePtr biop, Pointer userdata)

{
  SubSourcePtr       newssp = NULL;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  SubSourcePtr       ssp;
  ValNodePtr         pset;

  if (biop == NULL) return;

  pset = ParsePCRSet (biop);
  if (pset != NULL) {
    pset = ValNodeSort (pset, SortVnpByPCRSetSeq);
    pset = UniqueVnpByPCRSetSeq (pset);
    pset = ValNodeSort (pset, SortVnpByPCRSetOrder);
    newssp = WritePCRSet (pset);
    FreePCRSet (pset);
  }

  prev = &(biop->subtype);
  ssp = biop->subtype;
  while (ssp != NULL) {
    next = ssp->next;
    if (ssp->subtype == SUBSRC_fwd_primer_seq ||
        ssp->subtype == SUBSRC_rev_primer_seq ||
        ssp->subtype == SUBSRC_fwd_primer_name ||
        ssp->subtype == SUBSRC_rev_primer_name) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }

  *prev = newssp;
}

static void RemoveDuplicatePCRPrimers (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioSourcesInSep (sep, NULL, LookForDupPCR);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveAssemblyProc (BioseqPtr bsp, Pointer userdata)

{
  SeqHistPtr  hist;

  hist = bsp->hist;
  if (hist == NULL) return;
  hist->assembly = SeqAlignFree (hist->assembly);
  if (hist->assembly != NULL || hist->replace_date != NULL ||
      hist->replace_ids != NULL || hist->replaced_by_date != NULL ||
      hist->replaced_by_ids != NULL || hist->deleted_date != NULL) return;
  bsp->hist = SeqHistFree (bsp->hist);
}

static void RemoveSeqHistAssembly (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, RemoveAssemblyProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void LeaveFirstAssemblyProc (BioseqPtr bsp, Pointer userdata)

{
  SeqHistPtr   hist;
  SeqAlignPtr  next;
  SeqAlignPtr  sap;

  hist = bsp->hist;
  if (hist == NULL) return;
  sap = hist->assembly;
  if (sap == NULL) return;
  next = sap->next;
  if (next == NULL) return;
  sap->next = NULL;
  SeqAlignFree (next);
}

static void LeaveFirstSeqHistAssembly (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, LeaveFirstAssemblyProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void FindTPAAccns (SeqDescrPtr sdp, Pointer userdata)

{
  UserFieldPtr     curr;
  ValNodePtr PNTR  head;
  ObjectIdPtr      oip;
  CharPtr          str;
  UserFieldPtr     ufp;
  UserObjectPtr    uop;
  ValNodePtr       vnp;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "TpaAssembly") != 0) return;
  head = (ValNodePtr PNTR) userdata;
  if (head == NULL) return;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      vnp = ValNodeCopyStr (NULL, 0, str);
      vnp->next = *head;
      *head = vnp;
    }
  }
}

typedef struct streamfsa {
  FILE  *fp;
  Char  buf [128];
  Int2  idx;
  Int2  maxlen;
} StreamFsa, PNTR StreamFsaPtr;

static void LIBCALLBACK FsaStreamProc (
  CharPtr sequence,
  Pointer userdata
)

{
  Char          ch;
  StreamFsaPtr  sfp;

  if (StringHasNoText (sequence) || userdata == NULL) return;
  sfp = (StreamFsaPtr) userdata;

  ch = *sequence;
  while (ch != '\0') {
    sfp->buf [sfp->idx] = ch;
    (sfp->idx)++;

    if (sfp->idx >= sfp->maxlen) {
      sfp->buf [sfp->idx] = '\0';
      fprintf (sfp->fp, "%s\n", sfp->buf);
      sfp->idx = 0;
    }

    sequence++;
    ch = *sequence;
  }
}

static void StreamFastaForBsp (BioseqPtr bsp, FILE *fp)

{
  Char       id [41];
  StreamFsa  sf;

  if (bsp == NULL || fp == NULL) return;

  id [0] = '\0';
  SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  fprintf (fp, ">%s\n", id);

  MemSet ((Pointer) &sf, 0, sizeof (StreamFsa));
  sf.fp = fp;
  sf.idx = 0;
  sf.maxlen = 60;
  SeqPortStream (bsp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) &sf, FsaStreamProc);

  if (sf.idx > 0) {
    sf.buf [sf.idx] = '\0';
    fprintf (fp, "%s\n", sf.buf);
  }
}

static void CreateTPABlastDB (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  FILE         *fp;
  ValNodePtr   head = NULL;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  CharPtr      str;
  ValNodePtr   vnp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitDescriptorsInSep (sep, (Pointer) &head, FindTPAAccns);
  if (head == NULL) return;
  head = ValNodeSort (head, SortVnpByString);
  head = UniqueValNode (head);
  if (! GetOutputFileName (path, sizeof (path), NULL)) return;
  fp = FileOpen (path, "w");
  if (fp == NULL) return;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    Message (MSG_POST, "Making DB for %s", str);
    sip = SeqIdFromAccessionDotVersion (str);
    bsp = BioseqLockById (sip);
    StreamFastaForBsp (bsp, fp);
    BioseqUnlock (bsp);
  }
  FileClose (fp);
  ValNodeFreeData (head);
}

static void FindBioseqForHist (SeqIdPtr sip, Pointer userdata)

{
  BioseqPtr       bsp;
  BioseqPtr PNTR  bspp;

  bsp = BioseqFindCore (sip);
  if (bsp == NULL) return;
  bspp = (BioseqPtr PNTR) userdata;
  *bspp = bsp;
}

static void AttachSapToHistory (SeqAlignPtr sap)

{
  BioseqPtr   bsp = NULL;
  SeqHistPtr  shp;

  VisitSeqIdsInSeqAlign (sap, (Pointer) &bsp, FindBioseqForHist);
  if (bsp == NULL) return;
  shp = SeqHistNew ();
  shp->assembly = sap;
  bsp->hist = shp;
}

static void ReadSeqHistAssembly (IteM i)

{
  BaseFormPtr  bfp;
  Pointer      dataptr;
  Uint2        datatype;
  FILE         *fp;
  Char         path [PATH_MAX];
  SeqAlignPtr  salp;
  SeqAnnotPtr  sap;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
    if (datatype == OBJ_SEQANNOT) {
      sap = (SeqAnnotPtr) dataptr;
      if (sap != NULL && sap->type == 2) {
        salp = (SeqAlignPtr) sap->data;
        sap->data = NULL;
        AttachSapToHistory (salp);
        SeqAnnotFree (sap);
      }
    }
  }

  FileClose (fp);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct recondata {
  BioseqPtr   prod;
  SeqFeatPtr  cds;
  SeqFeatPtr  prt;
  Boolean     notunique;
} ReconData, PNTR ReconDataPtr;

static Boolean GetReconFunc (GatherContextPtr gcp)

{
  BioseqPtr     bsp;
  ReconDataPtr  rdp;
  SeqFeatPtr    sfp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  rdp = (ReconDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        if (rdp->cds != NULL) {
          rdp->notunique = TRUE;
        } else {
          rdp->cds = sfp;
        }
      } else if (sfp->data.choice == SEQFEAT_PROT) {
        if (rdp->prt == NULL) {
          rdp->prt = sfp; /* gets first protein, not largest, since location wrong */
        }
      }
      break;
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) gcp->thisitem;
      if (ISA_aa (bsp->mol)) {
        if (rdp->prod != NULL) {
          rdp->notunique = TRUE;
        } else {
          rdp->prod = bsp;
        }
      }
      break;
    default :
      break;
  }
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  return TRUE;
}

static void ReconnectCDSProc (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  GatherScope   gs;
  MolInfoPtr    mip;
  Boolean       partial5;
  Boolean       partial3;
  ReconData     rd;
  SeqLocPtr     slp;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ReconnectCDSProc (entityID, sep);
      }
      return;
    }
  }

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  rd.prod = NULL;
  rd.cds = NULL;
  rd.prt = NULL;
  rd.notunique = FALSE;
  partial5 = FALSE;
  partial3 = FALSE;
  GatherEntity (entityID, (Pointer) (&rd), GetReconFunc, &gs);
  if (rd.notunique) return;
  if (rd.prod != NULL && rd.cds != NULL) {
    slp = SeqLocFindNext (rd.cds->location, NULL);
    if (slp != NULL) {
      CheckSeqLocForPartial (slp, &partial5, &partial3);
    }
    sep = SeqMgrGetSeqEntryForData (rd.prod);
    if (sep != NULL) {
      SetSeqFeatProduct (rd.cds, rd.prod);
      if (rd.prt != NULL) {
        rd.prt->location = SeqLocFree (rd.prt->location);
        rd.prt->location = CreateWholeInterval (sep);
        SetSeqLocPartial (rd.prt->location, partial5, partial3);
        rd.prt->partial = (partial5 || partial3);
        vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
        if (vnp == NULL) {
          vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
          if (vnp != NULL) {
            mip = MolInfoNew ();
            vnp->data.ptrvalue = (Pointer) mip;
            if (mip != NULL) {
              mip->biomol = 8;
              mip->tech = 13;
            }
          }
        }
        if (vnp != NULL) {
          mip = (MolInfoPtr) vnp->data.ptrvalue;
          if (mip != NULL) {
            if (partial5 && partial3) {
              mip->completeness = 5;
            } else if (partial5) {
              mip->completeness = 3;
            } else if (partial3) {
              mip->completeness = 4;
            /*
            } else if (partial) {
              mip->completeness = 2;
            */
            } else {
              mip->completeness = 0;
            }
          }
        }
      }
    }
  }
}

static void ReconnectCDSProduct (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ReconnectCDSProc (bfp->input_entityID, sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoAddDeltaGap (BioseqPtr bsp, Pointer userdata)

{
  DeltaSeqPtr        del, dsp, next;
  SeqMgrDescContext  dcontext;
  Int4               len = 0;
  Int4Ptr            lenp;
  SeqLitPtr          litp;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;
  lenp = (Int4Ptr) userdata;
  if (lenp != NULL) {
    len = *lenp;
  }

  /* if has molinfo, but tech is not htgs1 or htgs2, bail */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->tech != MI_TECH_htgs_1 && mip->tech != MI_TECH_htgs_2) return;
    }
  }

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL && dsp->next != NULL; dsp = next) {
    next = dsp->next;
    if (dsp->choice != 2) continue;
    litp = (SeqLitPtr) dsp->data.ptrvalue;
    if (litp == NULL || IsDeltaSeqGap (dsp) || litp->length == 0) continue;

    if (next->choice == 2) {
      litp = (SeqLitPtr) next->data.ptrvalue;
      if (litp != NULL && IsDeltaSeqGap (next) && litp->length > 0) continue;
    }

    del = ValNodeNew (NULL);
    if (del == NULL) continue;
    litp = SeqLitNew ();
    if (litp == NULL) continue;
    litp->length = len;
    del->choice = 2;
    del->data.ptrvalue = litp;

    del->next = next;
    dsp->next = del;

    bsp->length += len;
  }
}

static void AddDelta100Gaps (IteM i)

{
  BaseFormPtr  bfp;
  Int4         len = 100;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, (Pointer) &len, DoAddDeltaGap);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void AddDelta0Gaps (IteM i)

{
  BaseFormPtr  bfp;
  Int4         len = 0;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, (Pointer) &len, DoAddDeltaGap);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoDeleteDeltaEndGap (BioseqPtr bsp, Pointer userdata)

{
  DeltaSeqPtr  dsp, last, prev;
  SeqLitPtr    litp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;

  dsp = (DeltaSeqPtr) bsp->seq_ext;
  if (dsp == NULL) return;

  /* if first DeltaSeq is gap, remove it */
  if (IsDeltaSeqGap (dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp)) {
    litp = (SeqLitPtr) dsp->data.ptrvalue;
    if (litp != NULL) {
      bsp->length -= litp->length;
      SeqLitFree (litp);
    }
    bsp->seq_ext = dsp->next;
    MemFree (dsp);
  }

  prev = (DeltaSeqPtr) bsp->seq_ext;
  if (prev == NULL) return;
  last = prev->next;
  if (last == NULL) return;
  while (last->next != NULL) {
    prev = last;
    last = last->next;
  }

  /* if last DeltaSeq is gap, remove it */
  if (IsDeltaSeqGap (last) && !DoesDeltaSeqHaveGapTypeOrLinkage(last)) {
    litp = (SeqLitPtr) last->data.ptrvalue;
    if (litp != NULL) {
      bsp->length -= litp->length;
      SeqLitFree (litp);
    }
    prev->next = NULL;
    MemFree (last);
  }

  /* combine adjacent gaps */
  CombineAdjacentGapsOnBioseq (bsp, NULL);
}

static void DeleteDeltaEndGaps (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, DoDeleteDeltaEndGap);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void AtccStrainToCultureCollection (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ValNodePtr   object_list, vnp;
  AECRSamplePtr       sample;
  SourceQualPairPtr   pair;
  ValNodePtr          field_from, field_to;
  AECRParseActionPtr  parse;
  AECRActionPtr       action;
  CharPtr             str1, str2, cp, new_str;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  parse = AECRParseActionNew ();

  parse->fields = ValNodeNew (NULL);
  parse->fields->choice = FieldPairType_source_qual;
  pair = SourceQualPairNew ();
  pair->field_from = Source_qual_strain;
  pair->field_to = Source_qual_culture_collection;
  parse->fields->data.ptrvalue = pair;

  parse->portion = TextPortionNew ();
  parse->portion->left_text = StringSave ("ATCC ");
  parse->portion->include_left = FALSE;
  parse->portion->right_text = StringSave ("");
  parse->portion->include_right = FALSE;
  parse->portion->inside = TRUE;
  parse->portion->case_sensitive = FALSE;
  parse->portion->whole_word = FALSE;

  parse->remove_from_parsed = FALSE;
  parse->remove_left = FALSE;
  parse->remove_right = FALSE;
  parse->existing_text = ExistingTextOption_add_qual;

  action = AECRActionNew ();
  action->action = ValNodeNew (NULL);
  action->action->choice = ActionChoice_parse;
  action->action->data.ptrvalue = parse;

  object_list = GetObjectListForAECRAction (sep, action);

  field_to = GetToFieldFromFieldPair (parse->fields);

  sample = GetAECRSampleFromObjectList (object_list, field_to);
  if (sample != NULL && sample->num_found > 0) {
    parse->existing_text = TwoStepExistingText (sample->num_found, FALSE, TRUE);
  }
  sample = AECRSampleFree (sample);
  if (parse->existing_text != 0) {
    WatchCursor ();
    Update ();
    field_from = GetFromFieldFromFieldPair (parse->fields);
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
      str1 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, NULL);
      str2 = GetTextPortionFromString (str1, parse->portion);
      if (str2 != NULL) {
        cp = StringChr (str2, ';');
        if (cp != NULL) {
          *cp = 0;
        }
        new_str = (CharPtr) MemNew (sizeof (Char) * (5 + StringLen (str2) + 1));
        sprintf (new_str, "ATCC:%s", str2);
        SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, new_str, parse->existing_text);
        new_str = MemFree (new_str);
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    }
    field_from = FieldTypeFree (field_from);
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    ArrowCursor ();
    Update ();   
  }
  field_to = FieldTypeFree (field_to);
  object_list = FreeObjectList (object_list);
  action = AECRActionFree (action);
}


typedef struct extendgeneformglobalitemdata {
  Int4		iFeatDef;
  Boolean	item_value;
  Boolean	upstream;
  CharPtr	label;
} ExtendGeneFormGlobalItemData, PNTR ExtendGeneFormGlobalItemPtr;

typedef struct extendgeneformlocalitemdata {
  ButtoN	item_button;
  Boolean	item_value;
  SeqFeatPtr	plus_sfp;
  SeqFeatPtr	minus_sfp;
} ExtendGeneFormLocalItemData, PNTR ExtendGeneFormLocalItemPtr;

static ExtendGeneFormGlobalItemData list_of_regulatory_items[] = {
	{ FEATDEF_attenuator, FALSE, FALSE, "attenuator" },
	{ FEATDEF_attenuator, FALSE, TRUE, "attenuator" },
	{ FEATDEF_CAAT_signal, FALSE, TRUE, "CAAT signal" },
	{ FEATDEF_enhancer, FALSE, TRUE, "enhancer" },
	{ FEATDEF_GC_signal, FALSE, TRUE, "GC signal" },
	{ FEATDEF_misc_binding, FALSE, TRUE, "misc binding" },
	{ FEATDEF_misc_feature, FALSE, TRUE, "misc feature" },
	{ FEATDEF_misc_signal, FALSE, TRUE, "misc signal" },
	{ FEATDEF_misc_structure, FALSE, TRUE, "misc structure" },
	{ FEATDEF_promoter, FALSE, TRUE, "promoter" },
	{ FEATDEF_protein_bind, FALSE, TRUE, "protein bind" },
	{ FEATDEF_RBS, TRUE, TRUE, "RBS" },
	{ FEATDEF_stem_loop, FALSE, TRUE, "stem loop" },
	{ FEATDEF_TATA_signal, FALSE, TRUE, "TATA signal" },
	{ FEATDEF_terminator, FALSE, FALSE, "terminator" },
	{ FEATDEF_5UTR, TRUE, TRUE, "5'UTR" },
	{ FEATDEF_3UTR, TRUE, FALSE, "3'UTR" },
	{ FEATDEF_10_signal, FALSE, TRUE, "-10 signal" },
	{ FEATDEF_35_signal, FALSE, TRUE, "-35 signal" },
};

const Int4 NUM_REGULATORY_ITEMS =
	sizeof(list_of_regulatory_items) / sizeof(ExtendGeneFormGlobalItemData); 

typedef struct extendgeneformdata {
  FEATURE_FORM_BLOCK

  ButtoN	PlusStrand;
  ButtoN	MinusStrand;
  Boolean	doPlusStrand;
  Boolean	doMinusStrand;
  ButtoN	ResetGenes;
  ButtoN	StealFeatures;
  Boolean	doResetGenes;
  Boolean	doStealFeatures;
  ButtoN	LogEvents;
  ButtoN	LogUnextendedEvents;
  Boolean	doLogEvents;
  Boolean	doLogUnextendedEvents;
  ExtendGeneFormLocalItemPtr	item_list;
  FILE		*fp;
} ExtendGeneFormData, PNTR ExtendGeneFormPtr;

static CharPtr GetGeneName(GeneRefPtr grp)
{
  if (!HasNoText (grp->locus)) return grp->locus;
  if (!HasNoText (grp->locus_tag)) return grp->locus_tag;
  if (!HasNoText (grp->desc)) return grp->desc;
  return NULL;
}

static void LogLocation (SeqMgrFeatContextPtr context, FILE *fp)
{
  if (context->strand == Seq_strand_minus)
  {
    fprintf (fp, "at complement(%d..%d)", context->left, context->right);
  }
  else
  {
    fprintf (fp, "at %d..%d", context->left, context->right);
  }
}

typedef struct correctgeneformdata {
  FEATURE_FORM_BLOCK

  ButtoN      setIntervalBtn;
  ButtoN      setStrandBtn;
  ButtoN      logChangesBtn;
  ButtoN      correct_selected_pair;
  GrouP       feature_choice_grp;
  Uint1       featdeftype;
  Boolean     setInterval;
  Boolean     setStrand;
  SeqEntryPtr sep;
  FILE *      fp;
} CorrectGeneFormData, PNTR CorrectGeneFormPtr;

typedef struct genesforcdsdata {
  SeqLocPtr      test_loc;
  ValNodePtr     gene_list;
  Boolean        same_strand_only;
} GenesForCDSData, PNTR GenesForCDSPtr;

static void FindGenesForCDS (SeqFeatPtr sfp, Pointer userdata)
{
  GenesForCDSPtr gfc;
  ValNodePtr     vnp;
  Uint1          strand_cds, strand_gene;
  
  gfc = (GenesForCDSPtr) userdata;

  if (gfc == NULL || sfp == NULL || sfp->idx.subtype != FEATDEF_GENE) 
  {
    return;
  }

  /* unless we plan to correct to same strand, only look at locations on same strand */
  if (gfc->same_strand_only)
  {
    strand_cds = SeqLocStrand (gfc->test_loc);
    strand_gene = SeqLocStrand (sfp->location);
    if (strand_cds == Seq_strand_minus && strand_gene != Seq_strand_minus)
    {
      return;
    }
    else if (strand_cds != Seq_strand_minus && strand_gene == Seq_strand_minus)
    {
      return;
    }
  }

  if (SeqLocAinB (gfc->test_loc, sfp->location) > -1
    || SeqLocAinB (sfp->location, gfc->test_loc) > -1)
  {
    vnp = ValNodeNew (gfc->gene_list);
    if (vnp == NULL) return;
    vnp->data.ptrvalue = sfp;
    if (gfc->gene_list == NULL)
    {
      gfc->gene_list = vnp;
    }
  }
}

extern void SetSeqLocStrand (SeqLocPtr location, Uint1 strand)
{
  SeqLocPtr slp;
  SeqIntPtr sip;
  SeqPntPtr spp;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL)
  {
    if (slp == NULL || slp->data.ptrvalue == NULL) 
    {
      slp = SeqLocFindNext (location, slp);
      continue;
    }
    if (slp->choice == SEQLOC_INT)
    {
      sip = (SeqIntPtr)slp->data.ptrvalue;
      sip->strand = strand;
    }
    else if (slp->choice == SEQLOC_PNT)
    {
      spp = (SeqPntPtr)slp->data.ptrvalue;
      spp->strand = strand;
    }
    slp = SeqLocFindNext (location, slp);
  }
}

static void 
CorrectOneCDSGenePair 
(SeqFeatPtr cds, SeqFeatPtr gene, BioseqPtr bsp, CorrectGeneFormPtr cgp)
{
  Boolean   need_change;
  SeqLocPtr slp;
  Boolean   oldpartial5, oldpartial3, partial5, partial3;
  CharPtr   cds_label;
  CharPtr            genename = NULL;
  SeqMgrFeatContext  fcontext;
  Uint1              strand_cds, strand_gene;
  Char               id_txt [128];
  
  if (cds == NULL || gene == NULL || cgp == NULL)
  {
    return;
  }
  
  if (bsp == NULL)
  {
    bsp = BioseqFindFromSeqLoc (cds->location);
    if (bsp == NULL)
    {
      return;
    }
  }
  
  if (cgp->setInterval)
  {
    need_change = FALSE;
    slp = SeqLocMerge (bsp, cds->location, NULL, TRUE, TRUE, FALSE);
    if (SeqLocCompare (slp, gene->location) != SLC_A_EQ_B)
    {
      need_change = TRUE;
    }
    CheckSeqLocForPartial (gene->location, &oldpartial5, &oldpartial3);
    CheckSeqLocForPartial (cds->location, &partial5, &partial3);
    if (oldpartial5 != partial5 || oldpartial3 != partial3)
    {
      need_change = TRUE;
    }
    if (need_change)
    {
      SetSeqLocPartial (slp, partial5, partial3);
      SeqLocFree (gene->location);
      gene->location = slp;
      gene->partial = partial5 || partial3;
      if (cgp->fp != NULL)
      {
        SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
        if (!StringHasNoText (id_txt))
        {
          fprintf (cgp->fp, "Sequence %s: ", id_txt);
        }
        genename = GetGeneName ((GeneRefPtr) gene->data.value.ptrvalue);
        cds = SeqMgrGetDesiredFeature (cds->idx.entityID, NULL, 0, 0, cds, &fcontext);
        if (genename == NULL)
        {
          fprintf (cgp->fp, "Unnamed gene ");
          LogLocation (&fcontext, cgp->fp);
        }
        else
        {
          fprintf (cgp->fp, "Gene %s", genename);
        }
        cds_label = StringSave (fcontext.label);
        if (StringLen (cds_label) > 20)
        {
          cds_label[20] = 0;
        }
        fprintf (cgp->fp, " reset to %s interval (%s)\n", cgp->featdeftype == FEATDEF_CDS ? "CDS" : "mRNA", cds_label);
        cds_label = MemFree (cds_label);
      }
    }
    else
    {
      SeqLocFree (slp);      
    }              
  }
  if (cgp->setStrand)
  {
    need_change = FALSE;
    strand_cds = SeqLocStrand (cds->location);
    strand_gene = SeqLocStrand (gene->location);
    if (strand_cds == Seq_strand_minus && strand_gene != Seq_strand_minus)
    {
      need_change = TRUE;
    }
    else if (strand_cds != Seq_strand_minus && strand_gene == Seq_strand_minus)
    {
      need_change = TRUE;
    }
    if (need_change)
    {
      SetSeqLocStrand (gene->location, strand_cds);
      if (cgp->fp != NULL)
      {
        cds = SeqMgrGetDesiredFeature (cds->idx.entityID, NULL, 0, 0, cds, &fcontext);
        genename = GetGeneName ((GeneRefPtr) gene->data.value.ptrvalue);
        if (genename == NULL)
        {
          fprintf (cgp->fp, "Unnamed gene ");
          LogLocation (&fcontext, cgp->fp);
        }
        else
        {
          fprintf (cgp->fp, "Gene %s", genename);
        }
        cds_label = StringSave (fcontext.label);
        if (StringLen (cds_label) > 20)
        {
          cds_label[20] = 0;
        }
        fprintf (cgp->fp, " strand set to match CDS (%s)\n", cds_label);
        cds_label = MemFree (cds_label);
      }
    }
  }
}
  
static void CorrectOneCDSGene (SeqFeatPtr sfp, Pointer userdata)
{
  CorrectGeneFormPtr cgp;
  GenesForCDSData    gfc;
  SeqFeatPtr         gene;
  SeqFeatXrefPtr     gene_xref;
  GeneRefPtr         grp;
  SeqMgrFeatContext  fcontext;
  CharPtr            cds_label;
  BioseqPtr          bsp;

  if (sfp == NULL 
      || (cgp = (CorrectGeneFormPtr)userdata) == NULL
      || sfp->idx.subtype != cgp->featdeftype)
  {
    return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  gfc.test_loc = SeqLocMerge (bsp, sfp->location, NULL, TRUE, TRUE, FALSE);

  gfc.gene_list = NULL;
  gfc.same_strand_only = ! cgp->setStrand;

  gene_xref = sfp->xref;
  while (gene_xref != NULL)
  {
    if (gene_xref->data.choice == SEQFEAT_GENE && gene_xref->data.value.ptrvalue != NULL)
    {
      grp = (GeneRefPtr) gene_xref->data.value.ptrvalue;
      if (!StringHasNoText (grp->locus))
      {
        gene = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
        if (gene != NULL)
        {
          ValNodeAddPointer (&(gfc.gene_list), 0, gene);
        }
      }
    }
    gene_xref = gene_xref->next;
  }

  if (gfc.gene_list == NULL)
  {
    VisitFeaturesInSep (cgp->sep, &gfc, FindGenesForCDS);
  }

  if (gfc.gene_list == NULL || gfc.gene_list->data.ptrvalue == NULL
      || gfc.gene_list->next != NULL)
  {
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &fcontext);
    cds_label = StringSave (fcontext.label);
    if (StringLen (cds_label) > 20)
    {
      cds_label[20] = 0;
    }
    if (cgp->fp != NULL)
    {
      if (gfc.gene_list == NULL || gfc.gene_list->data.ptrvalue == NULL)
      {
        fprintf (cgp->fp, "No gene found for %s (%s)\n", cgp->featdeftype == FEATDEF_CDS ? "CDS" : "mRNA", cds_label);
      }
      else
      {
        fprintf (cgp->fp, "Found more than one gene for %s (%s)\n", cgp->featdeftype == FEATDEF_CDS ? "CDS" : "mRNA", cds_label);
      }
    }
    cds_label = MemFree (cds_label);
  }
  else
  {
    gene = gfc.gene_list->data.ptrvalue;
    CorrectOneCDSGenePair (sfp, gene, bsp, cgp);
  }

  gfc.gene_list = ValNodeFree (gfc.gene_list);
  gfc.test_loc = SeqLocFree (gfc.test_loc);
}

static void CorrectGenesForCDSs (ButtoN b)
{
  CorrectGeneFormPtr cgp;
  Char			         path [PATH_MAX];
  SelStructPtr       ssp;
  SeqFeatPtr         gene = NULL, cds = NULL, sfp;
  SeqMgrFeatContext  fcontext;
  Boolean            cancel = FALSE;
  
  cgp = GetObjectExtra (b);
  if (cgp == NULL) return;

  Hide (cgp->form);
  WatchCursor ();
  Update ();

  cgp->setStrand = GetStatus (cgp->setStrandBtn);
  cgp->setInterval = GetStatus (cgp->setIntervalBtn);
  cgp->sep = GetTopSeqEntryForEntityID (cgp->input_entityID);

  if (GetValue (cgp->feature_choice_grp) == 1) {
    cgp->featdeftype = FEATDEF_CDS;
  } else {
    cgp->featdeftype = FEATDEF_mRNA;
  }

  /* set up tmp file if a log has been requested */  
  if (GetStatus (cgp->logChangesBtn))
  {
    TmpNam (path);
    cgp->fp = FileOpen (path, "wb");
  }
  else
  {
    cgp->fp = NULL;
  }
  
  if (GetStatus (cgp->correct_selected_pair))
  {
    ssp = ObjMgrGetSelected ();
    while (ssp != NULL && ! cancel)
    {
      if(ssp->itemtype == OBJ_SEQFEAT) 
      {
        sfp = SeqMgrGetDesiredFeature (cgp->input_entityID, NULL, ssp->itemID, 0, NULL, &fcontext);
        if(sfp != NULL && sfp->idx.subtype == cgp->featdeftype)
        {
          if (cds == NULL)
          {
            cds = sfp;
          }
          else
          {
            if (cgp->featdeftype == FEATDEF_CDS) {
              Message (MSG_ERROR, "You have selected more than one CDS!");
            } else {
              Message (MSG_ERROR, "You have selected more than one mRNA!");
            }
            cancel = TRUE;
          }
        }
        else if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE)
        {
          if (gene == NULL)
          {
            gene = sfp;
          }
          else
          {
            Message (MSG_ERROR, "You have selected more than one gene!");
            cancel = TRUE;
          }
        }
      }
      ssp = ssp->next;
    }
    
    if (!cancel && cds != NULL && gene != NULL)
    {
      CorrectOneCDSGenePair (cds, gene, NULL, cgp);
    }
    else if (!cancel)
    {
      if (cgp->featdeftype == FEATDEF_CDS) {
        Message (MSG_ERROR, "You need to select one CDS and one gene");
      } else {
        Message (MSG_ERROR, "You need to select one mRNA and one gene");
      }
    }
  }
  else
  {
    VisitFeaturesInSep (cgp->sep, cgp, CorrectOneCDSGene);
  }
  
  /* close, display, and remove log file */
  if (cgp->fp != NULL)
  {
    FileClose (cgp->fp);
    cgp->fp = NULL;
    LaunchGeneralTextViewer (path, "Gene Change Log");
    FileRemove (path);
  }
  
  ObjMgrSetDirtyFlag (cgp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cgp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  
  if (GetStatus (cgp->leave_dlg_up))
  {
    Show (cgp->form);
  }
  else
  {
    Remove (cgp->form);
  }

}

static void CorrectGenes (IteM i)
{
  BaseFormPtr        bfp;
  CorrectGeneFormPtr cgp;
  WindoW             w;
  GrouP              g, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cgp = MemNew (sizeof (CorrectGeneFormData));
  if (cgp == NULL) return;
  
  w = FixedWindow (-50, -33, -20, -10, "Correct Genes",
                   StdCloseWindowProc);
  SetObjectExtra (w, cgp, StdCleanupFormProc);
  cgp->form = (ForM) w;

  cgp->input_entityID = bfp->input_entityID;
  cgp->input_itemID = bfp->input_itemID;
  cgp->input_itemtype = bfp->input_itemtype;

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  cgp->setIntervalBtn = CheckBox (g, "Set gene interval to match", NULL);
  SetStatus (cgp->setIntervalBtn, TRUE);
  cgp->setStrandBtn = CheckBox (g, "Set gene strand to match", NULL);
  SetStatus (cgp->setStrandBtn, TRUE);
  cgp->logChangesBtn = CheckBox (g, "Log gene changes", NULL);
  SetStatus (cgp->logChangesBtn, FALSE);
  
  cgp->correct_selected_pair = CheckBox (g, "Correct only selected feature-gene pair", NULL);

  cgp->feature_choice_grp = HiddenGroup (g, 2, 0, NULL);
  RadioButton (cgp->feature_choice_grp, "CDS");
  RadioButton (cgp->feature_choice_grp, "mRNA");
  SetValue (cgp->feature_choice_grp, 1);

  c = HiddenGroup (g, 3, 0, NULL);
  b = DefaultButton(c, "Accept", CorrectGenesForCDSs);
  SetObjectExtra(b, cgp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cgp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER,
                (HANDLE) cgp->setIntervalBtn,
                (HANDLE) cgp->setStrandBtn,
                (HANDLE) cgp->correct_selected_pair,
                (HANDLE) cgp->feature_choice_grp,
	        (HANDLE) c, NULL);
  RealizeWindow(w);
  Show(w);
  Update();
}

static void ResetAllGenes (SeqEntryPtr sep,
			Boolean doMinus, Boolean doPlus, FILE *fp)
{
  SeqEntryPtr		sep_index;
  SeqFeatPtr		gene;
  SeqFeatPtr		cds;
  SeqMgrFeatContext  fcontext;
  SeqMgrFeatContext  cdscontext;
  CharPtr		genename;
  BioseqSetPtr		bssp;
  BioseqPtr		bsp;
  SeqLocPtr		slp;
  SeqLocPtr		tmpslp;
  Boolean               cds_partial5, cds_partial3;
  Boolean               partial5, partial3, partial_flag;

  if (sep == NULL) return;
  if ( IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep_index = bssp->seq_set;
	sep_index != NULL;
	sep_index = sep_index->next)
    {
      ResetAllGenes ( sep_index, doMinus, doPlus, fp);
    }
    return;
  }
  else if ( IS_Bioseq(sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
    while (gene != NULL)
    {
      if (fcontext.strand == Seq_strand_minus && !doMinus)
        continue;
      if (fcontext.strand != Seq_strand_minus && !doPlus)
        continue;
      slp = NULL;
      partial3 = FALSE;
      partial5 = FALSE;
      partial_flag = FALSE;
      /* look for next feature, starting right after the gene */
      /* an overlapping CDS won't occur before the gene in the indexing */
      cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &cdscontext);
      while (cds != NULL)
      {
        if (((fcontext.strand == Seq_strand_minus
             && cdscontext.strand == Seq_strand_minus)
            || (fcontext.strand != Seq_strand_minus
             && cdscontext.strand != Seq_strand_minus))
            && SeqLocAinB (cds->location, gene->location) > -1)
        {
          CheckSeqLocForPartial (cds->location, &cds_partial5, &cds_partial3);
          partial5 |= cds_partial5;
          partial3 |= cds_partial3;
          partial_flag |= cds->partial;
          tmpslp = SeqLocMerge (bsp, cds->location, slp, TRUE, TRUE, FALSE);
          if (tmpslp == NULL) return;
          SeqLocFree (slp);
          slp = tmpslp;
        }
        cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0,
					&cdscontext);
      }

      if (slp != NULL)
      {
        SetSeqLocPartial (slp, partial5, partial3);
        if (fp != NULL
          && (SeqLocStart (gene->location) != SeqLocStart (slp)
            || SeqLocStop (gene->location) != SeqLocStop (slp)))
        {
          genename = GetGeneName ((GeneRefPtr) gene->data.value.ptrvalue);
          if (genename == NULL)
          {
            fprintf (fp, "Unnamed gene ");
            LogLocation (&fcontext, fp);
          }
          else
          {
            fprintf (fp, "Gene %s", genename);
          }
          fprintf (fp, " reset to CDS interval\n");
        }
        SeqLocFree (gene->location);
        gene->location = slp;
        gene->partial = partial_flag;
      }
      gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &fcontext);
    }
  }
}

static void ExtendGeneMessageProc (ForM f, Int2 mssg)
{
  ExtendGeneFormPtr	egfp;

  egfp = (ExtendGeneFormPtr) GetObjectExtra(f);
  if (egfp != NULL) {
    if (egfp->appmessage != NULL) {
      egfp->appmessage (f, mssg);
    }
  }
}

static void CleanupExtendGenePage (GraphiC g, VoidPtr data)

{
  ExtendGeneFormPtr  egfp;

  egfp = (ExtendGeneFormPtr) data;
  if (egfp != NULL) {
    if (egfp->item_list != NULL)
      MemFree (egfp->item_list);
  }
  StdCleanupFormProc (g, data);
}

static void makeGeneXref (SeqFeatPtr sfpr, SeqFeatPtr sfpg)
{
  SeqFeatXrefPtr  curr;
  GeneRefPtr      grp;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;
  SeqFeatXrefPtr  xref;

  if (sfpr == NULL || sfpg == NULL || sfpg->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfpg->data.value.ptrvalue;
  if (grp == NULL) return;
  last = (SeqFeatXrefPtr PNTR) &(sfpr->xref);
  curr = sfpr->xref;
  while (curr != NULL) {
    next = curr->next;
    if (curr->data.choice == SEQFEAT_GENE) {
      *last = next;
      curr->next = NULL;
      SeqFeatXrefFree (curr);
    } else {
      last = &(curr->next);
    }
    curr = next;
  }
  grp = GeneRefDup (grp);
  if (grp == NULL) return;
  xref = SeqFeatXrefNew ();
  sfpr->xref = xref;
  if (xref != NULL) {
    xref->data.choice = SEQFEAT_GENE;
    xref->data.value.ptrvalue = (Pointer) grp;
  }
}

static void extendGeneInStreamOnStrand(ExtendGeneFormPtr egfp,
					BioseqPtr bsp,
					SeqFeatPtr GeneToExtend,
					SeqFeatPtr FarthestRegulator,
					Boolean stream,
					Boolean isPlusStrand)
{
  SeqLocPtr		slp;
  Int2			item_index;

  /* extend previous minus gene to last minus_upstream_regulator */
  slp = SeqLocMerge (bsp,
  	GeneToExtend->location, FarthestRegulator->location,
	TRUE, TRUE, FALSE);
  if (slp == NULL) return;
  GeneToExtend->location = SeqLocFree (GeneToExtend->location);
  GeneToExtend->location = slp;
  for (item_index = 0 ; item_index < NUM_REGULATORY_ITEMS; item_index++)
  {
    if (list_of_regulatory_items[item_index].upstream == stream)
    {
      if(isPlusStrand)
      {
        if (egfp->item_list[item_index].plus_sfp != NULL)
        {
          /* Make XRefs for all stream plus regulators */
          makeGeneXref (egfp->item_list[item_index].plus_sfp, GeneToExtend);

          /* Then reset the stream plus pointers */
          egfp->item_list[item_index].plus_sfp = NULL;
        }
      }
      else
      {
        if (egfp->item_list[item_index].minus_sfp != NULL)
        {
          /* Make XRefs for all stream minus regulators */
          makeGeneXref (egfp->item_list[item_index].minus_sfp, GeneToExtend);

          /* Then reset the stream minus pointers */
          egfp->item_list[item_index].minus_sfp = NULL;
        }
      }
    }
  }
}

typedef struct featureloginfo {
  SeqMgrFeatContext	context;
  CharPtr		featurename;
} FeatureLogInfo, FAR *FeatureLogInfoPtr;
 
typedef struct strandloginfo {
  FeatureLogInfo	upstream;
  FeatureLogInfo	downstream;
} StrandLogInfo, FAR *StrandLogInfoPtr;

static void LogExtend (SeqFeatPtr gene, SeqMgrFeatContextPtr genecontext,
		StrandLogInfoPtr log_info, Boolean doLogUnextended, FILE *fp)
{
  GeneRefPtr grp;
  CharPtr genename;

  if (fp == NULL) return;

  if (gene == NULL) return;

  grp = (GeneRefPtr) gene->data.value.ptrvalue;
  if (grp == NULL) return;

  if (!doLogUnextended && log_info->upstream.featurename == NULL && log_info->downstream.featurename == NULL)
    return;

  fprintf (fp, "Gene ");
  genename = GetGeneName (grp);
  if (genename != NULL)
  {
    fprintf (fp, "%s", genename);
  }
  else
  {
    LogLocation (genecontext, fp);
  }

  if (log_info->upstream.featurename == NULL && log_info->downstream.featurename == NULL)
  {
    fprintf (fp, " not extended\n");
  }
  else if (log_info->upstream.featurename != NULL)
  {
    fprintf (fp, " extended to %s ", log_info->upstream.featurename);
    LogLocation (&(log_info->upstream.context), fp);
    if (log_info->downstream.featurename != NULL)
    {
      fprintf (fp, ", %s ", log_info->downstream.featurename);
      LogLocation (&(log_info->downstream.context), fp);
    }
    fprintf (fp, "\n");
  }
  else
  {
    fprintf (fp, " extended to %s ", log_info->downstream.featurename);
    LogLocation (&(log_info->downstream.context), fp);
    fprintf (fp, "\n");
  }  
}

static void ExtendGenesOnSequence (SeqEntryPtr sep, ExtendGeneFormPtr egfp)
{
  SeqEntryPtr		sep_index;
  BioseqSetPtr		bssp;
  BioseqPtr		bsp;
  SeqFeatPtr		sfp;
  SeqMgrFeatContext	fcontext;
  SeqFeatPtr		plus_upstream_regulator;
  SeqFeatPtr		plus_downstream_regulator;
  SeqFeatPtr		plus_gene;
  SeqFeatPtr		minus_upstream_regulator;
  SeqFeatPtr		minus_downstream_regulator;
  SeqFeatPtr		minus_gene;
  Int2			item_index;
  GeneRefPtr		grp;
  SeqFeatPtr		overlapping_gene;
  StrandLogInfo		previous_plus_log;
  StrandLogInfo		next_plus_log;
  StrandLogInfo		previous_minus_log;
  StrandLogInfo		next_minus_log;

  if (sep == NULL) return;
  if ( IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep_index = bssp->seq_set;
	sep_index != NULL;
	sep_index = sep_index->next)
    {
      ExtendGenesOnSequence ( sep_index, egfp);
    }
    return;
  }
  else if ( ! IS_Bioseq(sep))
  {
    return;
  }
 
  bsp = (BioseqPtr) sep->data.ptrvalue;

  /* Walk through gene sequence, extending genes to specified regulatory */
  /* elements, respecting strand direction */
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);

  /* Loop through all features */

  plus_upstream_regulator = NULL;
  plus_downstream_regulator = NULL;
  plus_gene = NULL;
  minus_upstream_regulator = NULL;
  minus_downstream_regulator = NULL;
  minus_gene = NULL;
  previous_plus_log.upstream.featurename = NULL;
  previous_plus_log.downstream.featurename = NULL;
  previous_minus_log.upstream.featurename = NULL;
  previous_minus_log.downstream.featurename = NULL;
  next_plus_log.upstream.featurename = NULL;
  next_plus_log.downstream.featurename = NULL;
  next_minus_log.upstream.featurename = NULL;
  next_minus_log.downstream.featurename = NULL;

  while (NULL != sfp)
  {
    /* determine whether item has gene reference */
    overlapping_gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    grp = SeqMgrGetGeneXref (sfp);

    /* process strand items separately */
    if (fcontext.strand == Seq_strand_minus)
    {
      if (egfp->doMinusStrand)
      {
        if (fcontext.featdeftype == FEATDEF_GENE)
        {
          if (minus_gene != NULL && minus_upstream_regulator != NULL)
          {
            /* extend previous minus gene to last minus_upstream_regulator */
            extendGeneInStreamOnStrand(egfp, bsp, minus_gene,
					minus_upstream_regulator,
					TRUE, FALSE);
          }
          LogExtend (minus_gene, &fcontext, &previous_minus_log,
			egfp->doLogUnextendedEvents, egfp->fp);
          previous_minus_log = next_minus_log;
          next_minus_log.upstream.featurename = NULL;
          next_minus_log.downstream.featurename = NULL;
          minus_gene = sfp;
          if (minus_gene != NULL && minus_downstream_regulator != NULL)
          {
            /* extend this minus gene to minus_downstream_regulator */
            extendGeneInStreamOnStrand(egfp, bsp, minus_gene,
					minus_downstream_regulator,
					FALSE, FALSE);
          }
          minus_upstream_regulator = NULL;
          minus_downstream_regulator = NULL;
        }
        else
        {
          for (item_index = 0;
  		       item_index < NUM_REGULATORY_ITEMS;
  		       item_index++) 
  		  {
  		    if (list_of_regulatory_items[item_index].iFeatDef == fcontext.featdeftype)
            {
              if (egfp->item_list[item_index].item_value &&
		          ((grp == NULL && overlapping_gene == NULL) || egfp->doStealFeatures))
              {
                if (list_of_regulatory_items[item_index].upstream)
                {
                  minus_upstream_regulator = sfp;
                  previous_minus_log.upstream.featurename = list_of_regulatory_items[item_index].label;
                  previous_minus_log.upstream.context = fcontext;
                }
                else
                {
                  /* found downstream regulator */
                  if (minus_downstream_regulator == NULL)
                  {
                    minus_downstream_regulator = sfp;
                    next_minus_log.downstream.featurename = list_of_regulatory_items[item_index].label;
                    next_minus_log.downstream.context = fcontext;
                  }
                }
                egfp->item_list[item_index].minus_sfp = sfp;
              }
            }
          }
        }
      }
    }
    else /* treat all other conditions as plus strand */
    {
      if (egfp->doPlusStrand)
      {
        if (fcontext.featdeftype == FEATDEF_GENE)
        {
          if (plus_gene != NULL && plus_downstream_regulator != NULL)
          {
            /* extend previous plus gene to last plus_downstream_regulator */
            extendGeneInStreamOnStrand(egfp, bsp, plus_gene,
					plus_downstream_regulator,
					FALSE, TRUE);
          }
          LogExtend (plus_gene, &fcontext, &previous_plus_log, 
			egfp->doLogUnextendedEvents,egfp->fp);
          plus_gene = sfp;
          if (plus_gene != NULL && plus_upstream_regulator != NULL)
          {
            /* extend previous plus gene to first plus_upstream_regulator */
            extendGeneInStreamOnStrand(egfp, bsp, plus_gene,
					plus_upstream_regulator,
					TRUE, TRUE);
          }
          plus_upstream_regulator = NULL;
          plus_downstream_regulator = NULL;
          previous_plus_log = next_plus_log;
          next_plus_log.upstream.featurename = NULL;
          next_plus_log.downstream.featurename = NULL;
        }
        else
        {
          for (item_index = 0;
		       item_index < NUM_REGULATORY_ITEMS;
		       item_index++ )
		  {
		    if (list_of_regulatory_items[item_index].iFeatDef == fcontext.featdeftype)
            {
              if (egfp->item_list[item_index].item_value &&
                  ((grp == NULL && overlapping_gene == NULL) || egfp->doStealFeatures))
              {
                if (list_of_regulatory_items[item_index].upstream)
                {
                  if (plus_upstream_regulator == NULL)
                  {
                    plus_upstream_regulator = sfp;
                    next_plus_log.upstream.featurename = list_of_regulatory_items[item_index].label;
                    next_plus_log.upstream.context = fcontext;
                  }
                }
                else
                {
                  /* found downstream regulator */
                  plus_downstream_regulator = sfp;
                  previous_plus_log.downstream.featurename = list_of_regulatory_items[item_index].label;
                  previous_plus_log.downstream.context = fcontext;
                }
                egfp->item_list[item_index].plus_sfp = sfp;
              }
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  /* extend final minus strand gene */
  if (minus_gene != NULL && minus_upstream_regulator != NULL && egfp->doMinusStrand)
  {
    /* extend this minus gene to minus_upstream_regulator */
    extendGeneInStreamOnStrand(egfp, bsp, minus_gene,
					minus_upstream_regulator,
					TRUE, FALSE);
  }
  LogExtend (minus_gene, &fcontext, &previous_minus_log, 
			egfp->doLogUnextendedEvents,egfp->fp);

  /* extend final plus strand gene */
  if (plus_gene != NULL && plus_downstream_regulator != NULL && egfp->doPlusStrand)
  {
    /* extend this plus gene to plus_downstream_regulator */
    extendGeneInStreamOnStrand(egfp, bsp, plus_gene,
					plus_downstream_regulator,
					FALSE, TRUE);
  }
  LogExtend (plus_gene, &fcontext, &previous_plus_log, 
			egfp->doLogUnextendedEvents,egfp->fp);

}
  
static void DoExtendGene (ButtoN b)
{
  SeqEntryPtr		sep;
  ExtendGeneFormPtr	egfp;
  Char			path [PATH_MAX];
  Int2			item_index;

  egfp = GetObjectExtra (b);
  if (egfp == NULL) return;
  Hide (egfp->form);
  WatchCursor ();
  Update ();

  egfp->doLogEvents = GetStatus (egfp->LogEvents);
  egfp->doLogUnextendedEvents = GetStatus (egfp->LogUnextendedEvents);
  if (egfp->doLogEvents)
  {
    TmpNam (path);
    egfp->fp = FileOpen (path, "wb");
  }
  else
  {
    egfp->fp = NULL;
  }
  
  /* get button statuses */
  for (item_index = 0; item_index < NUM_REGULATORY_ITEMS; item_index ++ )
  {
    egfp->item_list[item_index].item_value = GetStatus (egfp->item_list[item_index].item_button);
  }

  egfp->doPlusStrand = GetStatus (egfp->PlusStrand);
  egfp->doMinusStrand = GetStatus (egfp->MinusStrand);
  egfp->doResetGenes = GetStatus (egfp->ResetGenes);
  egfp->doStealFeatures = GetStatus (egfp->StealFeatures);

  sep = GetTopSeqEntryForEntityID (egfp->input_entityID);

  if (egfp->doResetGenes)
  {
    ResetAllGenes (sep, egfp->doMinusStrand, egfp->doPlusStrand, egfp->fp);
    /* reindex features */
    SeqMgrIndexFeatures (egfp->input_entityID, NULL);
  }

  ExtendGenesOnSequence (sep, egfp);

  /* reindex features */
  SeqMgrIndexFeatures (egfp->input_entityID, NULL);

  /* remove redundant x-refs */
  SeqMgrExploreBioseqs (egfp->input_entityID, 0, NULL, RemoveGeneXrefsOnBioseqs, TRUE, TRUE, TRUE);
  
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (egfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, egfp->input_entityID, 0, 0);
  Remove (egfp->form);
  if (egfp->doLogEvents)
  {
    FileClose (egfp->fp);
    LaunchGeneralTextViewer (path, "Extend Gene Log");
    FileRemove (path);
  }
}

static void DoGeneReset (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr		sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  WatchCursor ();
  Update ();
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  ResetAllGenes (sep, TRUE, TRUE, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

static void DisableLogOption (ButtoN b)
{
  ExtendGeneFormPtr	egfp;
  Boolean	DoLogEvents;

  egfp = (ExtendGeneFormPtr) GetObjectExtra (b);
  DoLogEvents = GetStatus (egfp->LogEvents);
  if (DoLogEvents)
  {
    Enable (egfp->LogUnextendedEvents);
  }
  else
  {
    SetStatus (egfp->LogUnextendedEvents, FALSE);
    Disable (egfp->LogUnextendedEvents);
  }
}

static void ExtendGeneReg (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ExtendGeneFormPtr	egfp;
  WindoW	w;
  GrouP		g;
  GrouP		h;
  GrouP		j;
  GrouP		k;
  GrouP		UpstreamGroup;
  GrouP		DownstreamGroup;
  GrouP		m;
  GrouP		c;
  ButtoN	b;
  Int2		item_index;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  egfp = (ExtendGeneFormPtr) MemNew (sizeof (ExtendGeneFormData));
  if (egfp == NULL) return;
  w = FixedWindow (-50, -33, -20, -10, "Extend Gene", StdCloseWindowProc);
  SetObjectExtra (w, egfp, CleanupExtendGenePage);
  egfp->form = (ForM) w;
  egfp->formmessage = ExtendGeneMessageProc;

  egfp->input_entityID = bfp->input_entityID;
  egfp->input_itemID = bfp->input_itemID;
  egfp->input_itemtype = bfp->input_itemtype;

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  h = HiddenGroup (g, 0, 1, NULL);
  egfp->PlusStrand = CheckBox(h, "Extend on Plus Strand", NULL);
  SetStatus(egfp->PlusStrand, TRUE);
  egfp->MinusStrand = CheckBox(h, "Extend on Minus Strand", NULL);
  SetStatus(egfp->MinusStrand, TRUE);
  k = HiddenGroup (g, 0, 2, NULL);
  egfp->ResetGenes = CheckBox (k, "Reset Genes to CDS Locations before Extending", NULL);
  SetStatus (egfp->ResetGenes, FALSE);
  egfp->StealFeatures = CheckBox (k, "Reassign Features from Other Genes", NULL);
  SetStatus (egfp->StealFeatures, FALSE);
 
  egfp->item_list = MemNew (NUM_REGULATORY_ITEMS * sizeof (ExtendGeneFormLocalItemData));
  if (egfp->item_list == NULL) return;
  for (item_index = 0; item_index < NUM_REGULATORY_ITEMS; item_index ++)
  {
    egfp->item_list[item_index].item_button = NULL;
    egfp->item_list[item_index].item_value =
		list_of_regulatory_items[item_index].item_value;
    egfp->item_list[item_index].plus_sfp = NULL;
    egfp->item_list[item_index].minus_sfp = NULL;
  }   

  j = HiddenGroup (g, 1, 0, NULL);
  SetGroupSpacing (j, 10, 10);
  UpstreamGroup = NormalGroup (j, 0, 4, "Upstream", NULL, NULL);
  for(item_index = 0; item_index < NUM_REGULATORY_ITEMS; item_index ++)
  {
    if(list_of_regulatory_items[item_index].upstream)
    {
      egfp->item_list[item_index].item_button = CheckBox(UpstreamGroup, list_of_regulatory_items[item_index].label, NULL);
      SetStatus(egfp->item_list[item_index].item_button,
              list_of_regulatory_items[item_index].item_value);
    }
  }
  DownstreamGroup = NormalGroup (j, 4, 0, "Downstream", NULL, NULL);
  for(item_index = 0; item_index < NUM_REGULATORY_ITEMS; item_index ++)
  {
    if(! list_of_regulatory_items[item_index].upstream)
    {
      egfp->item_list[item_index].item_button = CheckBox(DownstreamGroup, list_of_regulatory_items[item_index].label, NULL);
      SetStatus(egfp->item_list[item_index].item_button,
              list_of_regulatory_items[item_index].item_value);
    }
  }
  m = HiddenGroup (g, 2, 0, NULL);
  egfp->LogEvents = CheckBox (m, "Log Changes", DisableLogOption);
  SetObjectExtra(egfp->LogEvents, egfp, NULL);
  SetStatus (egfp->LogEvents, FALSE);
  egfp->LogUnextendedEvents = CheckBox (m, "Log Unextended Genes", NULL);
  SetStatus (egfp->LogUnextendedEvents, FALSE);
  Disable (egfp->LogUnextendedEvents);

  c = HiddenGroup (g, 2, 0, NULL);
  b = DefaultButton(c, "Accept", DoExtendGene);
  SetObjectExtra(b, egfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) k,
	(HANDLE) j, (HANDLE) m, (HANDLE) c, NULL);
  RealizeWindow(w);
  Show(w);
  Update();
}

static void ResynchronizeCDSPartials (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  ResynchCodingRegionPartials (sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ResynchronizeMRNAPartials (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  ResynchMessengerRNAPartials (sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ResynchronizePeptidePartials (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  ResynchProteinPartials (sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void TrimProtsCallback (BioseqPtr bsp, Pointer userdata)

{
  SeqFeatPtr    sfp;
  SeqIntPtr     sintp;
  SeqLocPtr     slp;
  SeqMgrFeatContext context;

  if (bsp == NULL || !ISA_aa (bsp->mol)) return;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PROT, 0, &context)) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL && slp->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL) {
        if (sintp->from == 0 && sintp->to > bsp->length - 1) {
          sintp->to = bsp->length - 1;
        }
      }
    }
  }
}

static void TrimProtFeatLengthsEx (Uint2 entityID)

{
  SeqEntryPtr  sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, TrimProtsCallback);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void TrimProtFeatLengths (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  TrimProtFeatLengthsEx (bfp->input_entityID);
}

static void UnsetPartialForTruncatedProteins (SeqFeatPtr sfp, Pointer userdata)
{
  Boolean partial5, partial3;
  ProtRefPtr prp;
  CharPtr    new_name;
  BoolPtr    change_name;

  if (sfp == NULL || userdata == NULL) return;
  change_name = (BoolPtr) userdata;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  SetSeqLocPartial (sfp->location, partial5, FALSE);
  sfp->partial = partial5;
  if (*change_name) {
    prp = sfp->data.value.ptrvalue;
    if (prp != NULL && prp->name != NULL 
        && prp->name->data.ptrvalue != NULL
        && StringNCmp (prp->name->data.ptrvalue, "truncated ", 10) != 0) {
      new_name = (CharPtr) MemNew (StringLen (prp->name->data.ptrvalue)
                                   + 11);
      if (new_name != NULL) {
        sprintf (new_name, "truncated %s", (CharPtr) prp->name->data.ptrvalue);
        MemFree (prp->name->data.ptrvalue);
        prp->name->data.ptrvalue = new_name;
      }
    }
  }
}

static void UnsetCDSPartialForTruncatedProteins (SeqFeatPtr sfp, Pointer userdata) 
{
  Boolean     partial5, partial3;
  BioseqPtr   product_bsp, prot_bsp;
  CdRegionPtr crp;
  Int4        cds_len;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->location == NULL) return;
  if ((prot_bsp = (BioseqPtr)userdata) == NULL) return;

  product_bsp = BioseqFindFromSeqLoc (sfp->product);
  if (product_bsp == prot_bsp) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    SetSeqLocPartial (sfp->location, partial5, FALSE);
    sfp->partial = partial5;

    cds_len = (prot_bsp->length + 1) * 3 ;
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL)
    {
      if (crp->frame == 2)
      {
        cds_len += 1;
      }
      else if (crp->frame == 3)
      {
        cds_len += 2;
      }
    }
    sfp->location = TruncateLocation (sfp->location, cds_len);
  }
}

typedef struct truncprotsdata {
  Uint2       entityID;
  WindoW      w;
  ButtoN      reset_genes_btn;
  ButtoN      trim_prot_feats_btn;
  ButtoN      change_name_btn;
  ButtoN      retranslate_cds_btn;
  ButtoN      truncate_mrna_btn;
  
  SeqEntryPtr sep;
  Boolean     change_name;
  Boolean     trim_prot_feats;
  Boolean     retranslate_cds;
  Boolean     reset_genes;
  Boolean     truncate_mrna;
} TruncProtsData, PNTR TruncProtsPtr;

static void TruncProtsCallback (BioseqPtr bsp, Pointer userdata)

{
  ByteStorePtr  bs;
  BioseqPtr     cds_bsp = NULL;
  Int4          i;
  Int2          residue;
  SeqPortPtr    spp;
  Boolean       found_stop;
  TruncProtsPtr tpp;
  SeqFeatPtr    cds;
  SeqFeatPtr    gene = NULL, mRNA = NULL;
  GeneRefPtr    xref;
  SeqMgrFeatContext context;
  Boolean           partial5 = FALSE, partial3 = FALSE;
  Boolean           oldpartial5 = FALSE, oldpartial3 = FALSE;
  Boolean           need_change = FALSE;
  SeqLocPtr         slp;
  RecompData        rd;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_const) return;

  if ((tpp = (TruncProtsPtr)userdata) == NULL) return;

  /* get coding region for this protein */
  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
 
  if (cds != NULL)
  {
  
    cds_bsp = BioseqFindFromSeqLoc (cds->location);
  
    CheckSeqLocForPartial (cds->location, &partial5, &partial3);
    /* get gene for coding region */
    xref = SeqMgrGetGeneXref (cds);
    if (xref == NULL)
    {
      gene = SeqMgrGetOverlappingGene (cds->location, &context);
    }
    
    /* get mRNA for coding region */
    mRNA = SeqMgrGetOverlappingmRNA (cds->location, &context);
  }

  bs = BSNew (1000);
  if (bs == NULL) return;
  spp = SeqPortNew (bsp, 0, -1, 0, Seq_code_ncbieaa);
  if (spp == NULL) return;

  i = 0;
  residue = SeqPortGetResidue (spp);
  while (i < bsp->length && residue != '*') {
    BSPutByte (bs, residue);
    i++;
    residue = SeqPortGetResidue (spp);
  }
  if (residue == '*') {
    found_stop = TRUE;
  } else {
    found_stop = FALSE;
  }
  
  if (i >= bsp->length || !found_stop)
  {
    bs = BSFree (bs);
    SeqPortFree (spp);
    return;
  }

  SeqPortFree (spp);
  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bsp->seq_data_type = Seq_code_ncbieaa;
  
  VisitFeaturesOnBsp (bsp, (Pointer) &(tpp->change_name), 
                      UnsetPartialForTruncatedProteins);
  VisitFeaturesInSep (tpp->sep, bsp,  UnsetCDSPartialForTruncatedProteins);
  /* trim protein features */
  TrimProtsCallback (bsp, NULL);
  
  if (cds != NULL)
  {
    rd.no_stop_at_end_of_complete_cds = FALSE;
    rd.include_stop = FALSE;
    /* retranslate coding region associated with this protein, if requested */
    if (tpp->retranslate_cds) 
    {
      RetranslateOneCDS (cds, tpp->entityID, rd.include_stop, rd.no_stop_at_end_of_complete_cds);
    }
  }
  
  /* reset gene location if requested */
  if (tpp->reset_genes && gene != NULL && (mRNA == NULL || tpp->truncate_mrna))
  {
    need_change = FALSE;
    slp = SeqLocMerge (cds_bsp, cds->location, NULL, TRUE, TRUE, FALSE);
    if (SeqLocCompare (slp, gene->location) != SLC_A_EQ_B)
    {
      need_change = TRUE;
    }
    CheckSeqLocForPartial (gene->location, &oldpartial5, &oldpartial3);
    CheckSeqLocForPartial (cds->location, &partial5, &partial3);
    if (oldpartial5 != partial5 || oldpartial3 != partial3)
    {
      need_change = TRUE;
    }
    if (need_change)
    {
      SetSeqLocPartial (slp, partial5, partial3);
      SeqLocFree (gene->location);
      gene->location = slp;
      gene->partial = partial5 || partial3;
    }
    else
    {
      slp = SeqLocFree (slp);
    }
  }


  if (tpp->truncate_mrna && mRNA != NULL)
  {
    mRNA->location = TruncateLocation (mRNA->location, SeqLocLen (cds->location));    
  }
}

static void 
TruncProtsAtStopsExEx 
(Uint2 entityID,
 Boolean change_name,
 Boolean trim_prot_feats,
 Boolean retranslate_cds)

{
  SeqEntryPtr    sep;
  TruncProtsData tpd;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  tpd.sep = sep;
  tpd.change_name = change_name;
  tpd.trim_prot_feats = trim_prot_feats;
  tpd.retranslate_cds = retranslate_cds;
  VisitBioseqsInSep (sep, &tpd, TruncProtsCallback);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void DoTruncateProteins (ButtoN b)
{
  TruncProtsPtr     tpp;

  tpp = (TruncProtsPtr) GetObjectExtra (b);
  if (tpp == NULL)
  {
    return;
  }
  
  tpp->sep = GetTopSeqEntryForEntityID (tpp->entityID);
  if (tpp->sep == NULL) return;
  
  
  tpp->change_name = GetStatus (tpp->change_name_btn);
  tpp->trim_prot_feats = GetStatus (tpp->trim_prot_feats_btn);
  tpp->retranslate_cds = GetStatus (tpp->retranslate_cds_btn);
  tpp->reset_genes = GetStatus (tpp->reset_genes_btn);
  tpp->truncate_mrna = GetStatus (tpp->truncate_mrna_btn);
    
  VisitBioseqsInSep (tpp->sep, tpp, TruncProtsCallback);

  ObjMgrSetDirtyFlag (tpp->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tpp->entityID, 0, 0);
  Remove (tpp->w);
}

static void TruncateProteins (IteM i)
{
  BaseFormPtr       bfp;
  TruncProtsPtr     tpp;
  GrouP             g, k, c;
  ButtoN            b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  tpp = (TruncProtsPtr) MemNew (sizeof (TruncProtsData));
  if (tpp == NULL) return;
  
  tpp->w = FixedWindow (-50, -33, -20, -10, "Truncate Proteins at Stop Codons",
                   StdCloseWindowProc);
  SetObjectExtra (tpp->w, tpp, StdCleanupExtraProc);
  tpp->entityID = bfp->input_entityID;

  g = HiddenGroup (tpp->w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  k = HiddenGroup (g, 0, 6, NULL);
  SetGroupSpacing (g, 10, 10);
  tpp->trim_prot_feats_btn = CheckBox (k, "Trim Protein Features", NULL);
  SetStatus (tpp->trim_prot_feats_btn, TRUE);
  tpp->retranslate_cds_btn = CheckBox (k, "Retranslate Coding Regions for Truncated Proteins", NULL);
  SetStatus (tpp->retranslate_cds_btn, TRUE);
  tpp->truncate_mrna_btn = CheckBox (k, "Truncate Associated mRNA Features", NULL);
  SetStatus (tpp->truncate_mrna_btn, FALSE);
  tpp->reset_genes_btn = CheckBox (k, "Reset Genes", NULL);
  SetStatus (tpp->reset_genes_btn, TRUE);
  tpp->change_name_btn = CheckBox (k, "Prepend 'Truncated' to Protein Name", NULL);
  SetStatus (tpp->change_name_btn, TRUE);

  c = HiddenGroup (g, 2, 0, NULL);
  b = DefaultButton(c, "Accept", DoTruncateProteins);
  SetObjectExtra(b, tpp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow(tpp->w);
  Show(tpp->w);
  Update();
}

static void ExtendProtsForBioseq (BioseqPtr bsp, Int4 orig_len)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  SeqLocPtr         slp;
  SeqIntPtr         sintp;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PROT, 0, &context))
  {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL && slp->choice == SEQLOC_INT) 
    {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL) 
      {
        if (sintp->to == orig_len) 
        {
          sintp->to = bsp->length - 1;
        }
      }
    }
  }
}

static void ExtendProtsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr         bsp, cds_bsp = NULL;
  TruncProtsPtr     tpp;
  SeqFeatPtr        cds;
  SeqFeatPtr        gene = NULL, mRNA = NULL;
  GeneRefPtr        xref;
  SeqMgrFeatContext context;
  Boolean           partial5 = FALSE, partial3 = FALSE;
  Boolean           oldpartial5 = FALSE, oldpartial3 = FALSE;
  Boolean           need_change = FALSE;
  SeqLocPtr         slp;
  CharPtr           new_prot = NULL;
  Int4              orig_len, orig_cds_len;
  Int4              start, stop, len_extended;
  Uint1             strand;

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_const) return;

  if ((tpp = (TruncProtsPtr)mydata) == NULL) return;

  orig_len = bsp->length;
  /* get coding region for this protein */
  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds == NULL)
  {
    return;
  }
  
  orig_cds_len = SeqLocLen (cds->location);
  
  cds_bsp = BioseqFindFromSeqLoc (cds->location);
  
  /* get gene for coding region */
  xref = SeqMgrGetGeneXref (cds);
  if (xref == NULL)
  {
    gene = SeqMgrGetOverlappingGene (cds->location, &context);
  }
    
  /* get mRNA for coding region */
  mRNA = SeqMgrGetOverlappingmRNA (cds->location, &context);

  new_prot = ExtendProtein3 (cds, cds->idx.entityID, FALSE);
  len_extended = SeqLocLen (cds->location) - orig_cds_len;
  if (new_prot == NULL 
      || (StringLen (new_prot) == bsp->length && len_extended == 0))
  {
    new_prot = MemFree (new_prot);
    return;
  }
 
  CheckSeqLocForPartial (cds->location, &partial5, &partial3);
  if (!partial3) {
    VisitFeaturesOnBsp (bsp, (Pointer) &(tpp->change_name), 
                        UnsetPartialForTruncatedProteins);
  }

  /* retranslate coding region associated with this protein, if requested */
  if (tpp->retranslate_cds)
  {
    RetranslateOneCDS (cds, tpp->entityID, FALSE, FALSE); 
  
    /* extend protein features if requested */
    if (tpp->trim_prot_feats)
    {
      ExtendProtsForBioseq (bsp, orig_len);
    }
  }

  /* reset gene location if requested */
  if (tpp->reset_genes && gene != NULL && (mRNA == NULL || tpp->truncate_mrna))
  {
    need_change = FALSE;
    slp = SeqLocMerge (cds_bsp, cds->location, NULL, TRUE, TRUE, FALSE);
    if (SeqLocCompare (slp, gene->location) != SLC_A_EQ_B)
    {
      need_change = TRUE;
    }
    CheckSeqLocForPartial (gene->location, &oldpartial5, &oldpartial3);
    CheckSeqLocForPartial (cds->location, &partial5, &partial3);
    if (oldpartial5 != partial5 || oldpartial3 != partial3)
    {
      need_change = TRUE;
    }
    if (need_change)
    {
      SetSeqLocPartial (slp, partial5, partial3);
      SeqLocFree (gene->location);
      gene->location = slp;
    }
    else
    {
      slp = SeqLocFree (slp);
    }
  }

  if (tpp->truncate_mrna && mRNA != NULL)
  {
    start = GetOffsetInBioseq (mRNA->location, cds_bsp, SEQLOC_START);
    stop = GetOffsetInBioseq (mRNA->location, cds_bsp, SEQLOC_STOP);
    strand = SeqLocStrand (mRNA->location);
  
    mRNA->location = ExpandSeqLoc (start, stop + len_extended,
                                   strand, cds_bsp, mRNA->location);   
 
  }
}

static void DoExtendProteins (ButtoN b)
{
  TruncProtsPtr     tpp;

  tpp = (TruncProtsPtr) GetObjectExtra (b);
  if (tpp == NULL)
  {
    return;
  }
  
  tpp->sep = GetTopSeqEntryForEntityID (tpp->entityID);
  if (tpp->sep == NULL) return;
  
  
  tpp->change_name = FALSE;
  tpp->trim_prot_feats = GetStatus (tpp->trim_prot_feats_btn);
  tpp->retranslate_cds = GetStatus (tpp->retranslate_cds_btn);
  tpp->reset_genes = GetStatus (tpp->reset_genes_btn);
  tpp->truncate_mrna = GetStatus (tpp->truncate_mrna_btn);
    
  SeqEntryExplore (tpp->sep, tpp, ExtendProtsCallback);

  ObjMgrSetDirtyFlag (tpp->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tpp->entityID, 0, 0);
  Remove (tpp->w);
}

static void ExtendProteins (IteM i)
{
  BaseFormPtr       bfp;
  TruncProtsPtr     tpp;
  GrouP             g, k, c;
  ButtoN            b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  tpp = (TruncProtsPtr) MemNew (sizeof (TruncProtsData));
  if (tpp == NULL) return;
  
  tpp->w = FixedWindow (-50, -33, -20, -10, "Extend Proteins to Stop Codons",
                   StdCloseWindowProc);
  SetObjectExtra (tpp->w, tpp, StdCleanupExtraProc);
  tpp->entityID = bfp->input_entityID;

  g = HiddenGroup (tpp->w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  k = HiddenGroup (g, 0, 6, NULL);
  SetGroupSpacing (g, 10, 10);
  tpp->trim_prot_feats_btn = CheckBox (k, "Extend Protein Features", NULL);
  SetStatus (tpp->trim_prot_feats_btn, TRUE);
  tpp->retranslate_cds_btn = CheckBox (k, "Retranslate Coding Regions for Extended Proteins", NULL);
  SetStatus (tpp->retranslate_cds_btn, TRUE);
  tpp->truncate_mrna_btn = CheckBox (k, "Extend Associated mRNA Features", NULL);
  SetStatus (tpp->truncate_mrna_btn, FALSE);
  tpp->reset_genes_btn = CheckBox (k, "Reset Genes", NULL);
  SetStatus (tpp->reset_genes_btn, TRUE);
  tpp->change_name_btn = NULL;

  c = HiddenGroup (g, 2, 0, NULL);
  b = DefaultButton(c, "Accept", DoExtendProteins);
  SetObjectExtra(b, tpp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow(tpp->w);
  Show(tpp->w);
  Update();
  
}

typedef struct trimnandlog
{
  SeqEntryPtr top_sep;
  LogInfoPtr  lip;
} TrimNAndLogData, PNTR TrimNAndLogPtr;

extern void LogTrimmedLocation (LogInfoPtr lip, SeqLocPtr slp)
{
  CharPtr loc_str;
  
  if (lip != NULL && lip->fp != NULL && slp != NULL)
  {
    loc_str = SeqLocPrintUseBestID(slp);
    fprintf (lip->fp, "%s\n", loc_str);
		MemFree(loc_str);
		lip->data_in_log = TRUE;
  }
    
}


static void BioseqTrimNAndLog (BioseqPtr bsp, Pointer userdata)
{
  TrimNAndLogPtr tnalp;
  SeqIdPtr       sip;
  SeqLocPtr      slp1 = NULL,
                 slp2 = NULL;
  CharPtr       str;
  Int4          j, lens;
    
  tnalp = (TrimNAndLogPtr) userdata;
  if (bsp == NULL || ! ISA_na (bsp->mol) || tnalp == NULL)
  {
    return;
  }
  
  sip = SeqIdFindBest (bsp->id, 0);
  str = load_seq_data (sip, -1, -1, FALSE, &lens);
  if (str != NULL) 
  {
     j = lens-1;
     while (j>0) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j--;
     }
     if (j<lens-1) 
     {
        slp1 = SeqLocIntNew (j+1, lens-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp1, TRUE, FALSE);
        TrimQualityScores (bsp, lens - 1 - j, FALSE);        
     }
     j=0;
     while (j<lens) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j++;
     }
     if (j>0) {
        slp2 = SeqLocIntNew (0, j-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp2, TRUE, FALSE);
        TrimQualityScores (bsp, j, TRUE);        
     }
     if (slp1!=NULL) {
        LogTrimmedLocation (tnalp->lip, slp1);
        if (tnalp->top_sep!=NULL)
           SeqEntryExplore (tnalp->top_sep, (Pointer)slp1, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp1);
     }
     if (slp2!=NULL) {
        LogTrimmedLocation (tnalp->lip, slp2);
        if (tnalp->top_sep!=NULL)
           SeqEntryExplore (tnalp->top_sep, (Pointer)slp2, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp2);
     }
  }
}


static void TrimNsFromNucsCommon (Uint2 entityID)

{
  TrimNAndLogData tnald;

  tnald.top_sep = GetTopSeqEntryForEntityID (entityID);
  if (tnald.top_sep == NULL) return;
  
  tnald.lip = OpenLog ("Trimmed Locations");
  VisitBioseqsInSep (tnald.top_sep, &tnald, BioseqTrimNAndLog);
  
  CloseLog (tnald.lip);
  tnald.lip = FreeLog (tnald.lip);
  
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void TrimNsFromNucs (IteM i)
{
  BaseFormPtr     bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  TrimNsFromNucsCommon (bfp->input_entityID);
  
}

static void TrimNsFromNucsToolBtn (ButtoN b)

{
  BaseFormPtr  bfp;

  if (b == NULL) return;
  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  TrimNsFromNucsCommon (bfp->input_entityID);
}

static Boolean s_NextQualifies (SeqFeatPtr sfp,
				BioseqPtr  bsp,
				Int2       level,
				SeqMgrFeatContextPtr context)
{
  SeqFeatPtr        nextSfp;
  SeqMgrFeatContext fContext;
  GBQualPtr         gbq;
  Boolean           isvector;

  /* Get the next feature */

  fContext = *context;
  nextSfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature,
				  &fContext);
  if (nextSfp == NULL)
    return FALSE;

  /* Make sure that it's a Vector Contamination misc feat */

  for (isvector = FALSE, gbq = nextSfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "standard_name") == 0) {
      if (StringCmp (gbq->val, "Vector Contamination") == 0) {
        isvector = TRUE;
      }
    }
  }
  if (! isvector)
    return FALSE;

  /* Check the match level against the level that we're looking for */

  for (isvector = FALSE, gbq = nextSfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "phenotype") == 0) {
      if (StringCmp (gbq->val, "Strong match") == 0) {
        isvector = TRUE;
      } else if (StringCmp (gbq->val, "Moderate match") == 0 && level > 1) {
        isvector = TRUE;
      } else if (StringCmp (gbq->val, "Weak match") == 0 && level > 2) {
        isvector = TRUE;
      }
    }
  }
  if (! isvector)
    return FALSE;

  /* If we've made it this far, then the next feature */
  /* qualifies to be overwritten with N's.            */

  return TRUE;
}

static Boolean LIBCALLBACK VecToN (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)

{
  BioseqPtr      bsp;
  GBQualPtr      gbq;
  Int2           level;
  Boolean        isvector;
  SeqLocPtr      slp, subslp;
  Int4           start, stop, swap;
  VecToNDataPtr  vDataPtr;

  if (sfp == NULL || context == NULL) return TRUE;
  if (context->featdeftype != FEATDEF_misc_feature) return TRUE;
  vDataPtr = (VecToNDataPtr) context->userdata;
  if (vDataPtr == NULL) return TRUE;
  level = vDataPtr->level;

  /* Make sure that it's a Vector Contamination misc feat */

  for (isvector = FALSE, gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "standard_name") == 0) {
      if (StringCmp (gbq->val, "Vector Contamination") == 0) {
        isvector = TRUE;
      }
    }
  }
  if (! isvector) {
    vDataPtr->prevWasModified = FALSE;
    return TRUE;
  }

  /* Check the match level against the level that we're looking for */

  for (isvector = FALSE, gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "phenotype") == 0) {
      if (StringCmp (gbq->val, "Strong match") == 0) {
        isvector = TRUE;
      } else if (StringCmp (gbq->val, "Moderate match") == 0 && level > 1) {
        isvector = TRUE;
      } else if (StringCmp (gbq->val, "Weak match") == 0 && level > 2) {
        isvector = TRUE;
      } else if (StringCmp (gbq->val, "Suspect origin") == 0) {
	if (vDataPtr->prevWasModified ||
	    s_NextQualifies (sfp, context->bsp, level,
			     context))
	  isvector = TRUE;
      }
    }
  }
  if (! isvector) {
    vDataPtr->prevWasModified = FALSE;
    return TRUE;
  }

  /* based on BioseqOverwrite */

  bsp = context->bsp;
  if (bsp == NULL || (! ISA_na (bsp->mol))) return TRUE;
  if (bsp->repr == Seq_repr_seg) {
    slp = SegLocToPartsEx (bsp, sfp->location, FALSE);
    bsp = BioseqFind(SeqLocId(slp));
    if (bsp == NULL) return TRUE;
  } else {
    slp = sfp->location;
  }
  if (bsp->repr != Seq_repr_raw) {
    if (slp != sfp->location) {
      slp = SeqLocFree (slp);
    }
    return TRUE;
  }
  if (bsp->seq_data_type == Seq_code_gap) return TRUE;

  if (bsp->seq_data_type != Seq_code_iupacna) {
    BioseqRawConvert (bsp, Seq_code_iupacna);
  }
  
  subslp = SeqLocFindNext (slp, NULL);
  while (subslp != NULL) {
    start = SeqLocStart (subslp);
    stop = SeqLocStop (subslp);
    if (start > stop) {
      swap = start;
      start = stop;
      stop = swap;
    }
    if (start <= stop) {
      BSSeek ((ByteStorePtr) bsp->seq_data, start, SEEK_SET);
      while (start <= stop) {
        BSPutByte ((ByteStorePtr) bsp->seq_data, (Int2) 'N');
        start++;
      }
    }
    subslp = SeqLocFindNext (slp, subslp);
  }
  vDataPtr->prevWasModified = TRUE;
  if (slp != sfp->location) {
    slp = SeqLocFree (slp);
  }
  return TRUE;
}

static void VecScreenBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  VecToNDataPtr     vDataPtr;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || userdata == NULL) {
    return;
  }

  vDataPtr = (VecToNDataPtr) userdata;
    
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext)) {
       fcontext.userdata = vDataPtr;
       VecToN (sfp, &fcontext);
  }

}

static void VecScreenToNs (IteM i, Int2 level)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  VecToNData   vData;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  WatchCursor ();
  Update ();
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  vData.level = level;
  vData.prevWasModified = FALSE;

  bsp = GetBioseqGivenIDs (bfp->input_entityID,
			   bfp->input_itemID,
			   bfp->input_itemtype);
  if (NULL == bsp)
    VisitBioseqsInSep (sep, &vData, VecScreenBioseqCallback);
  else
  {
    VecScreenBioseqCallback (bsp, &vData);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();
}

static void VecScreenToNsStrong (IteM i)

{
  VecScreenToNs (i, 1);
}

static void VecScreenToNsModerate (IteM i)

{
  VecScreenToNs (i, 2);
}

static void VecScreenToNsWeak (IteM i)

{
  VecScreenToNs (i, 3);
}

static void AddGenomesUserObject (IteM i)

{
  BaseFormPtr    bfp;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjectIdPtr    oip;
  SeqEntryPtr    sep;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;
  ValNodePtr     vnp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL || sep->data.ptrvalue == NULL) return;

  oip = ObjectIdNew ();
  oip->str = StringSave ("info");
  uop = UserObjectNew ();
  uop->type = oip;
  uop->_class = StringSave ("Genomes");

  oip = ObjectIdNew ();
  oip->id = 0;
  ufp = UserFieldNew ();
  ufp->choice = 2;
  ufp->data.intvalue = 0;
  ufp->label = oip;

  uop->data = ufp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    vnp = SeqDescrNew (NULL);
    vnp->choice = Seq_descr_user;
    vnp->data.ptrvalue = (Pointer) uop;
    vnp->next = bsp->descr;
    bsp->descr = vnp;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    vnp = SeqDescrNew (NULL);
    vnp->choice = Seq_descr_user;
    vnp->data.ptrvalue = (Pointer) uop;
    vnp->next = bssp->descr;
    bssp->descr = vnp;
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean AddBlastTag (GatherContextPtr gcp)

{
  Uint1        align_type;
  SeqAnnotPtr  sap;

  align_type = (Uint1) (Int4) gcp->userdata;
  if (gcp->thistype != OBJ_SEQANNOT) return TRUE;
  sap = (SeqAnnotPtr) gcp->thisitem;
  if (sap == NULL || sap->type != 2) return TRUE;
  AddAlignInfoToSeqAnnot (sap, align_type);
  return TRUE;
}

static void CommonAddBlastToSeqAnnot (IteM i, Int2 align_type)

{
  BaseFormPtr   bfp;
  SelStructPtr  sel;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sel = ObjMgrGetSelected ();
  if (sel != NULL && sel->next == NULL &&
      sel->entityID == bfp->input_entityID &&
      sel->itemtype == OBJ_SEQANNOT) {
    GatherItem (sel->entityID, sel->itemID, sel->itemtype, (Pointer) (Int4) align_type, AddBlastTag);
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void AddBlastNToSeqAnnot (IteM i)

{
  CommonAddBlastToSeqAnnot (i, ALIGN_BLASTN);
}

static void AddBlastPToSeqAnnot (IteM i)

{
  CommonAddBlastToSeqAnnot (i, ALIGN_BLASTP);
}

static void AddBlastXToSeqAnnot (IteM i)

{
  CommonAddBlastToSeqAnnot (i, ALIGN_BLASTX);
}

static void AddTBlastNToSeqAnnot (IteM i)

{
  CommonAddBlastToSeqAnnot (i, ALIGN_TBLASTN);
}

static void ToolBtn1 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  SaveSeqSubmitProc (bfp, FALSE);
}

static void ToolBtn2 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  AutoDefBaseFormCommon (bfp, FALSE, TRUE);
}

/*
static void ToolBtn3 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  ForceTaxonFixupBtn (NULL, b);
}
*/

static void ToolBtn3 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  ForceCleanupBtn (NULL, b, FALSE);
}

static void ToolBtn4 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  CommonAddOrgOrModsToDefLines (NULL, 0, 0, b);
}

static void ToolBtn5 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  WatchCursor ();
  Update ();
  GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, bfp->input_itemID,
                    bfp->input_itemtype, OBJ_SEQFEAT, FEATDEF_CDS, 0, 0);
  ArrowCursor ();
  Update ();
}

static void ToolBtn6 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;


  WatchCursor ();
  Update ();
  GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, bfp->input_itemID,
                    bfp->input_itemtype, OBJ_SEQFEAT, FEATDEF_rRNA, 0, 0);
  ArrowCursor ();
  Update ();
}


static void ToolBtnGene (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;


  WatchCursor ();
  Update ();
  GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, bfp->input_itemID,
                    bfp->input_itemtype, OBJ_SEQFEAT, FEATDEF_GENE, 0, 0);
  ArrowCursor ();
  Update ();
}


static void ToolBtn7 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  MRnaFromCdsProc (bfp->input_entityID);
}

static void ToolBtn8 (ButtoN b)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  LaunchOrfViewer (bsp, bfp->input_entityID, bfp->input_itemID, FALSE);
}

static void ToolBtn9 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;


  WatchCursor ();
  Update ();
  GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, bfp->input_itemID,
                    bfp->input_itemtype, OBJ_SEQFEAT, FEATDEF_misc_feature, 0, 0);
  ArrowCursor ();
  Update ();
}

static void ToolBtn10 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  SUCCommonProc (NULL, FALSE, TRUE, TRUE, b);
}

static void ToolBtn11 (ButtoN b)

{
  BaseFormPtr  bfp;
  Int2         mssgsub;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  mssgsub = RegisterFormMenuItemName ("SequinEditSubmitterItem");
  SendMessageToForm (bfp->form, mssgsub);
}

static Boolean MakePubOnTopSep (GatherContextPtr gcp)

{
  Pointer  ptrvalue;

  if (gcp == NULL) return TRUE;
  ptrvalue = (Pointer) gcp->userdata;
  if (ptrvalue == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ || gcp->thistype == OBJ_BIOSEQSET) {
    if (ptrvalue == gcp->thisitem) {

      WatchCursor ();
      Update ();
      GatherProcLaunch (OMPROC_EDIT, FALSE, gcp->entityID, gcp->itemID,
                        gcp->thistype, OBJ_SEQDESC, Seq_descr_pub, 0, 0);
      ArrowCursor ();
      Update ();
      return FALSE;
    }
  }
  return TRUE;
}

static void ToolBtn12 (ButtoN b)

{
  BaseFormPtr  bfp;
  GatherScope  gs;
  SeqEntryPtr  sep;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL || sep->data.ptrvalue == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQSET] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  GatherEntity (bfp->input_entityID, (Pointer) sep->data.ptrvalue, MakePubOnTopSep, &gs);
}

static void ToolBtn13 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  EditPubsEx (bfp);
}

static void ToolBtn14 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  ValSeqEntryForm (bfp->form);
}

static void ToolBtn15 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  VSeqMgrShow ();
}

static void ToolBtn16 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  AddCitSubForUpdateProc (bfp);
}


typedef struct addsetform {
  FORM_MESSAGE_BLOCK
  PopuP class_control;
} AddSetFormData, PNTR AddSetFormPtr;

static void InsertSetBtn (ButtoN b) {
  AddSetFormPtr frm;
  SeqEntryPtr   sep;
  UIEnum        val;
  BioseqSetPtr  bssp, newbssp;
  ObjMgrDataPtr omdptop;
  ObjMgrData    omdata;
  Uint2         parenttype;
  Pointer       parentptr;


  frm = (AddSetFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }
  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  if (sep == NULL) {
    return;
  }

  WatchCursor ();
  Update ();

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  RemoveDupGenBankSets (sep);
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      GetEnumPopup (frm->class_control, bioseqset_class_alist, &val);
      newbssp = InsetNewSet (bssp, (Uint1) val);
      if (newbssp != NULL) {
        newbssp->descr = bssp->descr;
        bssp->descr = NULL;
        newbssp->annot = bssp->annot;
        bssp->annot = NULL;
      }
    }
  }
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);

  ArrowCursor ();
  Remove (frm->form);
  Update ();
}


static void AddSetToolBtn (ButtoN b)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  WindoW       w;
  GrouP        h, c;
  AddSetFormPtr frm;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) {
    return;
  }

  /* get class for new set */
  frm = (AddSetFormPtr) MemNew (sizeof (AddSetFormData));
  frm->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Bioseq Set", NULL);
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -2, 0, NULL);
  StaticPrompt (h, "Class", 0, popupMenuHeight, programFont, 'l');
  frm->class_control = PopupList (h, TRUE, NULL);
  InitEnumPopup (frm->class_control, bioseqset_class_alist, NULL);
  SetEnumPopup (frm->class_control, bioseqset_class_alist, (UIEnum) BioseqseqSet_class_phy_set);

  c = HiddenGroup (w, 2, 0, NULL);
  b = DefaultButton (c, "Accept", InsertSetBtn);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
}


static void MarkHypotheticalProteinTitles (SeqDescrPtr sdp, Pointer userdata)
{
  CharPtr pDefLine;
  CharPtr pCh;
  CharPtr pMatch = "hypothetical protein";
  size_t match_len = StringLen(pMatch);
  ObjValNodePtr ovn;

  if (sdp == NULL || sdp->extended == 0) return;
  if (sdp->choice != Seq_descr_title) return;
  pDefLine = sdp->data.ptrvalue;
  if (pDefLine == NULL) return;
  if(StringNCmp(pDefLine, pMatch, match_len) != 0)
  {
    return;
  }
  pCh = pDefLine + match_len;
  if(*pCh == 0
	|| (*pCh == '.' && *(pCh + 1) == 0)
	|| (*pCh == ' ' && *(pCh + 1) == '['))
  {
    ovn = (ObjValNodePtr) sdp;
    ovn->idx.deleteme = TRUE;
  }
}

static void RemoveBioseqHypotheticalProteinTitles (BioseqPtr bsp, Pointer userdata)
{
  if (bsp == NULL) return;
  if (! ISA_aa(bsp->mol)) return;
  VisitDescriptorsOnBsp(bsp, userdata, MarkHypotheticalProteinTitles);
}

static void RemoveHypotheticalProteinTitles (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, RemoveBioseqHypotheticalProteinTitles);
  DeleteMarkedObjects(bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void RemoveProtTitlesProc (BaseFormPtr bfp)

{
  BioseqPtr    bsp;
  SeqEntryPtr  fsep;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;

  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && IS_Bioseq (fsep)) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_OTHER) {
          Message (MSG_OK, "Use Remove ALL Protein Titles for RefSeq records");
          break;
        }
      }
    }
  }

  ClearProteinTitlesInNucProts (bfp->input_entityID, NULL);

  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoRemoveProtTitles (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  RemoveProtTitlesProc (bfp);
}

static void ToolBtn17 (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  RemoveProtTitlesProc (bfp);
}

static void RemoveBioseqInconsistentProteinTitles (BioseqPtr bsp, Pointer userdata)
{
  BioseqSetPtr       bssp;
  Char               buf [1001];
  ObjValNodePtr      ovn;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) {
    vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
    if (vnp != NULL) {
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bsp->idx.parentptr;
        while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot) {
          if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bssp->idx.parentptr;
          } else {
            bssp = NULL;
          }
        }
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
          title = (CharPtr) vnp->data.ptrvalue;
          if (NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), TRUE, FALSE)) {
            if (StringICmp (buf, title) != 0) {
              ovn = (ObjValNodePtr) vnp;
              ovn->idx.deleteme = TRUE;
            }
          }
        }
      }
    }
  }
}

static void RemoveInconsistentProteinTitles (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, RemoveBioseqInconsistentProteinTitles);
  DeleteMarkedObjects(bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static ButtoN SqnPushButton (GrouP prnt, CharPtr title, BtnActnProc actn, BaseFormPtr bfp)

{
  ButtoN  b;

  b = PushButton (prnt, title, actn);
  SetObjectExtra (b, (Pointer) bfp, NULL);
  return b;
}

static void FixCapitalizationToolBtn (ButtoN b);

extern void BioseqViewFormToolBar (GrouP h)

{
  BaseFormPtr  bfp;
  GrouP        g, g2;
  ButtoN       b1, b2, b3, b4, b5, b_top;

  bfp = (BaseFormPtr) GetObjectExtra (h);
  if (bfp == NULL) return;
  g = HiddenGroup (h, 1, 0, NULL);
  if (indexerVersion)
  {
    SqnPushButton (g, "Auto_Def", AutoDefToolBtn, bfp);
    b_top = SqnPushButton (g, "Auto_Def Options", AutoDefOptionsToolBtn, bfp);
    g2 = HiddenGroup (g, 2, 0, NULL);
    SqnPushButton (g2, "  Ad_ID  ", AutoDefStrainToolBtn, bfp);
    b1 = SqnPushButton (g2, "Ad_misc", AutoDefMiscFeatToolBtn, bfp);

    if (useTaxon) {
      SqnPushButton (g, "Tax_Fix/Clean_Up", ToolBtn3, bfp);
    }
    g2 = HiddenGroup (g, 2, 0, NULL);
    SqnPushButton (g2, "   CDS   ", ToolBtn5, bfp);
    b2 = SqnPushButton (g2, "rRNA", ToolBtn6, bfp);
    SqnPushButton (g2, "misc_feat", ToolBtn9, bfp);
    b3 = SqnPushButton (g2, "gene", ToolBtnGene, bfp);
    SqnPushButton (g, "mRna_CDS", ToolBtn7, bfp);
    SqnPushButton (g, "ORF_Find", ToolBtn8, bfp);
    SqnPushButton (g, "SUC", ToolBtn10, bfp);
    SqnPushButton (g, "rRNA->DNA", RibosomalRNAToGenomicDNAToolBtn, bfp);
    SqnPushButton (g, "Add Set", AddSetToolBtn, bfp);
    SqnPushButton (g, "Remove DefLines", RemoveDefLinesToolBtn, bfp);
    SqnPushButton (g, "Trim Ns", TrimNsFromNucsToolBtn, bfp);
    g2 = HiddenGroup (g, 2, 0, NULL);
    SqnPushButton (g2, "Find ASN.1", FindStringProcToolBtn, bfp);
    b4 = SqnPushButton (g2, "Find FF", FindFlatfileProcToolBtn, bfp);
    SqnPushButton (g, "Local ID->Src", ParseLocalIDToSourceQual, bfp);
    SqnPushButton (g, "Fix Local IDs", ResolveExistingLocalIDsToolBtn, bfp);
    SqnPushButton (g, "Group Explode", GroupExplodeToolBtn, bfp);
    SqnPushButton (g, "cit-sub-upd", ToolBtn16, bfp);
    g2 = HiddenGroup (g, 2, 0, NULL);
    SqnPushButton (g2, "Fix Caps", FixCapitalizationToolBtn, bfp);
    b5 = SqnPushButton (g2, "Edit Pubs", ToolBtn13, bfp);
    SqnPushButton (g, "rem_prot_titles", ToolBtn17, bfp);
    SqnPushButton (g, "Validate", ToolBtn14, bfp);
    SqnPushButton (g, "Desktop", ToolBtn15, bfp);
    AlignObjects (ALIGN_RIGHT, (HANDLE) b1, (HANDLE) b2, (HANDLE) b3, (HANDLE) b4, (HANDLE) b5, (HANDLE) b_top, NULL);

  }
  else
  {
    SqnPushButton (g, "Save", ToolBtn1, bfp);
    SqnPushButton (g, "Auto_Def", ToolBtn2, bfp);
    if (useTaxon) {
      /* SqnPushButton (g, "Tax_Fix", ToolBtn3, bfp); */
      SqnPushButton (g, "Tax_Fix/Clean_Up", ToolBtn3, bfp);
    }
    SqnPushButton (g, "Def_Org", ToolBtn4, bfp);
    SqnPushButton (g, "CDS", ToolBtn5, bfp);
    SqnPushButton (g, "rRNA", ToolBtn6, bfp);
    SqnPushButton (g, "mRna_CDS", ToolBtn7, bfp);
    SqnPushButton (g, "ORF_Find", ToolBtn8, bfp);
    SqnPushButton (g, "misc_feat", ToolBtn9, bfp);
    SqnPushButton (g, "SUC", ToolBtn10, bfp);
    SqnPushButton (g, "sub_affil", ToolBtn11, bfp);
    SqnPushButton (g, "sub_add", ToolBtn12, bfp);
    SqnPushButton (g, "cit-sub-upd", ToolBtn16, bfp);
    SqnPushButton (g, "Edit Pubs", ToolBtn13, bfp);
    SqnPushButton (g, "rem_prot_titles", ToolBtn17, bfp);
    SqnPushButton (g, "Validate", ToolBtn14, bfp);
    SqnPushButton (g, "Desktop", ToolBtn15, bfp);
  }
}

static void MakeToolBarWindow (IteM i)

{
  BaseFormPtr  bfp;
  ForM         f;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  f = MakeToolFormForBioseqView (bfp, BioseqViewFormToolBar);
  Show (f);
  Select (f);
}

static void DoNCCleanup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  NC_Cleanup (bfp->input_entityID, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void DoSSECleanup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  SeriousSeqEntryCleanup (sep, NULL, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void DoNormalizeDescs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  NormalizeDescriptorOrder (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void MarkProtTitlesProc (BioseqPtr bsp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      if (sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }
}

static void DoRemoveAllProtTitles (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, MarkProtTitlesProc);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void MarkNucTitlesProc (BioseqPtr bsp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      if (sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }
}

static void DoRemoveAllNucTitles (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, MarkNucTitlesProc);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void MarkMrnaTitlesProc (BioseqPtr bsp, Pointer userdata)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  ObjValNodePtr      ovp;
  SeqDescrPtr        sdp;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  if (mip->biomol != MOLECULE_TYPE_MRNA) return;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      if (sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }
}

static void DoRemoveMrnaTitles (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, MarkMrnaTitlesProc);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void MarkNPSTitlesProc (BioseqSetPtr bssp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;

  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      if (sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }
}

static void DoRemoveNPSTitles (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitSetsInSep (sep, NULL, MarkNPSTitlesProc);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void DoPTCleanup (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  InstantiateProteinTitles (bfp->input_entityID, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void TakeProteinsFromGPS (BioseqPtr bsp, Pointer userdata)

{
  SeqEntryPtr PNTR  lastp;
  SeqEntryPtr       sep;

  if (bsp == NULL || (! ISA_aa (bsp->mol))) return;
  lastp = (SeqEntryPtr PNTR) userdata;
  if (lastp == NULL) return;

  /* unlink from existing chain */

  sep = bsp->seqentry;
  if (sep != NULL) {
    sep->data.ptrvalue = NULL;
  }

  /* link after genomic sequence */

  sep = ValNodeAddPointer (lastp, 1, (Pointer) bsp);
  *lastp = sep;
}

static void StripMrnaProducts (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->product == NULL) return;

  sfp->product = SeqLocFree (sfp->product);
}

static void GPStoNPS (IteM i)

{
  BaseFormPtr   bfp;
  BioseqSetPtr  bssp;
  BioseqSet     dum;
  SeqEntryPtr   last, sep, top;
  Uint2         parenttype;
  Pointer       parentptr;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  top = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (top == NULL || top->choice != 2) return;
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
    sep = bssp->seq_set;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
    }
  }
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_gen_prod_set) return;

  GetSeqEntryParent (top, &parentptr, &parenttype);

  sep = bssp->seq_set;
  if (sep == NULL || sep->choice != 1) return;

  /* unlink nuc-prot sets, etc., from genomic Bioseq */

  MemSet ((Pointer) &dum, 0, sizeof (BioseqSet));
  dum._class = 1;
  dum.seq_set = sep->next;
  sep->next = NULL;

  last = sep;
  VisitBioseqsInSet (&dum, (Pointer) &last, TakeProteinsFromGPS);

  /* leave dum.seq_set dangling for now */

  bssp->_class = BioseqseqSet_class_nuc_prot;

  SeqMgrLinkSeqEntry (top, parenttype, parentptr);

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  VisitFeaturesInSet (bssp, NULL, StripMrnaProducts);

  move_cds (top);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void CDStoNPS (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  move_cds (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void AddNewMrnaTitle (BioseqPtr bsp, Pointer userdata)

{
  BioSourcePtr       biop;
  SeqMgrFeatContext  ccontext;
  CharPtr            cdslabel = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  gcontext;
  CharPtr            genelabel = NULL;
  size_t             len;
  MolInfoPtr         mip;
  CharPtr            organism;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            str;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  if (mip->biomol != MOLECULE_TYPE_MRNA) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || sdp->choice != Seq_descr_source) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  organism = orp->taxname;

  if (BioseqGetTitle (bsp) != NULL) return;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (sfp != NULL) {
    genelabel = gcontext.label;
    if (StringHasNoText (genelabel)) {
      genelabel = NULL;
    }
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  if (sfp != NULL) {
    cdslabel = ccontext.label;
    if (StringHasNoText (cdslabel)) {
      cdslabel = NULL;
    }
  }
  len = StringLen (organism) + StringLen (genelabel) + StringLen (cdslabel) +
        StringLen (" mRNA, complete cds.") + 10;
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str == NULL) return;
  str [0] = '\0';

  if (! StringHasNoText (organism)) {
    StringCat (str, organism);
  }
  if (cdslabel != NULL) {
    StringCat (str, " ");
    StringCat (str, cdslabel);
  }
  if (genelabel != NULL) {
      StringCat (str, " (");
      StringCat (str, genelabel);
      StringCat (str, ")");
  }
  if (cdslabel != NULL && genelabel != NULL) {
    if (ccontext.partialL || ccontext.partialR) {
      StringCat (str, " mRNA, partial cds.");
    } else {
      StringCat (str, " mRNA, complete cds.");
    }
  } else if (genelabel != NULL) {
    StringCat (str, " mRNA.");
  }
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
}

static void MrnaTitlesToGPS (IteM i)

{
  BaseFormPtr   bfp;
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep, top;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  top = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (top == NULL || top->choice != 2) return;
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
    sep = bssp->seq_set;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
    }
  }
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_gen_prod_set) return;

  VisitBioseqsInSep (top, NULL, AddNewMrnaTitle);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

typedef struct npsseqs {
  BioseqPtr  nuc;
  BioseqPtr  prot;
} NpsSeqs, PNTR NpsSeqsPtr;

static void FindNucProtSeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  NpsSeqsPtr  nsp;

  if (bsp == NULL) return;
  nsp = (NpsSeqsPtr) userdata;
  if (nsp == NULL) return;

  if (ISA_na (bsp->mol)) {
    nsp->nuc = bsp;
  } else if (ISA_aa (bsp->mol)) {
    nsp->prot = bsp;
  }
}

static void SqnMakeNucProtCDS (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CodeBreakPtr  cbp;
  SeqFeatPtr    cds;
  CdRegionPtr   crp;
  SeqFeatPtr    mrna;
  NpsSeqs       ns;
  Boolean       partial5, partial3;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Int4          start, stop;
  SeqFeatPtr    temp;

  ns.nuc = NULL;
  ns.prot = NULL;
  if (VisitBioseqsInSet (bssp, (Pointer) &ns, FindNucProtSeqs) != 2) return;
  if (ns.nuc == NULL || ns.prot == NULL) return;

  cds = SeqMgrGetCDSgivenProduct (ns.prot, NULL);
  mrna = SeqMgrGetRNAgivenProduct (ns.nuc, NULL);
  if (cds == NULL || mrna == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  start = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_START);
  stop = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_STOP);

  if (start < 0 || start >= ns.nuc->length ||
      stop < 0 || stop >= ns.nuc->length) return;

  sip = SeqIdFindBest (ns.nuc->id, 0);
  if (sip == NULL) return;

  /* copy cds feature fields to paste into new cds feature */
  temp = AsnIoMemCopy (cds,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* clear feature ID xrefs on cds */
  ClearFeatIDXrefs (temp);

  sfp = CreateNewFeatureOnBioseq (ns.nuc, SEQFEAT_CDREGION, NULL);
  if (sfp == NULL) return;

  sfp->location = SeqLocFree (sfp->location);
  if (StringISearch (cds->except_text, "ribosomal slippage") == NULL &&
      StringISearch (cds->except_text, "ribosome slippage") == NULL &&
      StringISearch (cds->except_text, "trans splicing") == NULL &&
      StringISearch (cds->except_text, "trans-splicing") == NULL &&
      StringISearch (cds->except_text, "artificial frameshift") == NULL) {
    sfp->location = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);
  } else {
    slp = SeqLocFindNext (cds->location, NULL);
    while (slp != NULL) {
      start = GetOffsetInLoc (slp, mrna->location, SEQLOC_START);
      stop = GetOffsetInLoc (slp, mrna->location, SEQLOC_STOP);
      sfp->location = AddIntervalToLocation (sfp->location, sip, start, stop, partial5, partial3);
      slp = SeqLocFindNext (cds->location, slp);
    }
    sfp->location = SeqLocMergeEx (ns.nuc, sfp->location, NULL, FALSE, TRUE, FALSE, FALSE);
  }
  SetSeqFeatProduct (sfp, ns.prot);

  /* paste fields from temp copy of original cds */
  crp = (CdRegionPtr) temp->data.value.ptrvalue;
  sfp->data.value.ptrvalue = (Pointer) crp;

  sfp->partial = temp->partial;
  sfp->excpt = temp->excpt;
  sfp->comment = temp->comment;
  sfp->qual = temp->qual;
  sfp->title = temp->title;
  sfp->ext = temp->ext;
  sfp->cit = temp->cit;
  sfp->exp_ev = temp->exp_ev;
  sfp->xref = temp->xref;
  sfp->dbxref = temp->dbxref;
  sfp->pseudo = temp->pseudo;
  sfp->except_text = temp->except_text;

  MemFree (temp); /* do not SeqFeatFree */

  /* update code break locations */
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
    CheckSeqLocForPartial (cbp->loc, &partial5, &partial3);
    start = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_START);
    stop = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_STOP);
    if (start < 0 || start >= ns.nuc->length ||
        stop < 0 || stop >= ns.nuc->length) continue;
	cbp->loc = SeqLocFree (cbp->loc);
	cbp->loc = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);;
  }
}

/* copy gene from contig onto nuc-prot, single interval on cdna bioseq */

static void SqnCopyGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene, copy, temp;
  GeneRefPtr         grp, xref;
  Boolean            partial5, partial3;

  /* input mrna features are multi-interval on contig */

  if (sfp->data.choice != SEQFEAT_RNA) return;

  /* find cdna product of mrna */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  /* check for gene xref */

  xref = SeqMgrGetGeneXref (sfp);
  if (xref != NULL) {
    if (SeqMgrGeneIsSuppressed (xref)) return;

    /* copy gene xref for new gene feature */

    grp = AsnIoMemCopy (xref,
                        (AsnReadFunc) GeneRefAsnRead,
                        (AsnWriteFunc) GeneRefAsnWrite);
    if (grp == NULL) return;

    /* make new gene feature on full-length of cdna */

    copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
    if (copy == NULL) return;

    copy->data.value.ptrvalue = grp;
    return;
  }

  /* overlapping gene should be single interval on contig */

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;

  CheckSeqLocForPartial (gene->location, &partial5, &partial3);

  /* copy gene feature fields to paste into new gene feature */

  temp = AsnIoMemCopy (gene,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* make new gene feature on full-length of cdna */

  copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
  if (copy == NULL) {
    SeqFeatFree (temp);
    return;
  }

  /* paste fields from temp copy of original gene */

  copy->data.value.ptrvalue = temp->data.value.ptrvalue;
  copy->partial = temp->partial;
  copy->excpt = temp->excpt;
  copy->comment = temp->comment;
  copy->qual = temp->qual;
  copy->title = temp->title;
  copy->ext = temp->ext;
  copy->cit = temp->cit;
  copy->exp_ev = temp->exp_ev;
  copy->xref = temp->xref;
  copy->dbxref = temp->dbxref;
  copy->pseudo = temp->pseudo;
  copy->except_text = temp->except_text;

  SetSeqLocPartial (copy->location, partial5, partial3);

  SeqLocFree (temp->location);
  MemFree (temp); /* do not SeqFeatFree */
}

static void SqnAddMrnaTitles (
  SeqLocPtr slp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  CharPtr            cdslabel = NULL;
  SeqMgrFeatContext  gcontext;
  CharPtr            genelabel = NULL;
  size_t             len;
  CharPtr            organism;
  SeqFeatPtr         sfp;
  CharPtr            str;

  if (slp == NULL) return;
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  organism = (CharPtr) userdata;
  if (BioseqGetTitle (bsp) != NULL) return;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (sfp != NULL) {
    genelabel = gcontext.label;
    if (StringHasNoText (genelabel)) {
      genelabel = NULL;
    }
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  if (sfp != NULL) {
    cdslabel = ccontext.label;
    if (StringHasNoText (cdslabel)) {
      cdslabel = NULL;
    }
  }
  len = StringLen (organism) + StringLen (genelabel) + StringLen (cdslabel) +
        StringLen (" mRNA, complete cds.") + 10;
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str == NULL) return;
  str [0] = '\0';

  if (! StringHasNoText (organism)) {
    StringCat (str, organism);
  }
  if (cdslabel != NULL) {
    StringCat (str, " ");
    StringCat (str, cdslabel);
  }
  if (genelabel != NULL) {
      StringCat (str, " (");
      StringCat (str, genelabel);
      StringCat (str, ")");
  }
  if (cdslabel != NULL && genelabel != NULL) {
    if (ccontext.partialL || ccontext.partialR) {
      StringCat (str, " mRNA, partial cds.");
    } else {
      StringCat (str, " mRNA, complete cds.");
    }
  } else if (genelabel != NULL) {
    StringCat (str, " mRNA.");
  }
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
}

static Boolean LIBCALLBACK DummyACMProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  return TRUE;
}

static Boolean AmbiguousCdsMrnaByOverlap (BioseqPtr bsp)

{
  Int2               count;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;

  if (bsp == NULL) return TRUE;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                             CHECK_INTERVALS, NULL, DummyACMProc);
    if (count > 1) return TRUE;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  return FALSE;
}

static SeqFeatPtr GetmRNAByFeatureID (SeqFeatPtr cds)
{
  SeqFeatPtr         mRNA = NULL;
  SeqFeatXrefPtr     xref;
  ObjectIdPtr        oip;
  Char               buf [32];

  if (cds == NULL) return NULL;
  
  xref = cds->xref;
  while (xref != NULL && mRNA == NULL) {
    if (xref->id.choice == 3) {
      oip = (ObjectIdPtr) xref->id.value.ptrvalue;
      if (oip != NULL) {
        if (StringDoesHaveText (oip->str)) {
          mRNA = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, oip->str, NULL, NULL);
        } else {
          sprintf (buf, "%ld", (long) oip->id);
          mRNA = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, buf, NULL, NULL);
        }
      }
    }
    xref = xref->next;
  }
  return mRNA;
}

static Boolean AmbiguousCdsMrnaByFeatureID (BioseqPtr bsp)

{
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;
  
  if (bsp == NULL) return TRUE;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    if (GetmRNAByFeatureID(sfp) == NULL) return TRUE;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  return FALSE;
}


static SeqIdPtr MakeIdFromLocusTag (SeqFeatPtr mrna)

{
  Char        buf [64];
  SeqFeatPtr  gene;
  GeneRefPtr  grp;

  if (mrna == NULL) return NULL;
  grp = SeqMgrGetGeneXref (mrna);
  if (grp != NULL) {
    if (SeqMgrGeneIsSuppressed (grp)) return NULL;
  }
  if (grp == NULL) {
    gene = SeqMgrGetOverlappingGene (mrna->location, NULL);
    if (gene != NULL) {
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
    }
  }
  if (grp != NULL) {
    if (StringDoesHaveText (grp->locus_tag)) {
      /* StringCpy (buf, "lcl|"); */
      StringCpy (buf, "gnl|MTRACK|");
      StringCat (buf, grp->locus_tag);
      return MakeSeqID (buf);
    }
  }
  return NULL;
}

static void InstantiateMrnaIntoProt (SeqFeatPtr cds, SeqFeatPtr mrna, Int2Ptr ctrp, Boolean useLocusTag)

{
  ByteStorePtr  bs;
  Int4          len;
  BioseqPtr     mbsp, pbsp;
  MolInfoPtr    mip;
  SeqEntryPtr   msep, psep;
  Boolean       partial5, partial3;
  CharPtr       rnaseq;
  ValNodePtr    vnp;

  if (cds == NULL || mrna == NULL) return;
  if (cds->product == NULL || mrna->product != NULL) return;

  pbsp = BioseqFindFromSeqLoc (cds->product);
  if (pbsp  == NULL) return;
  psep = SeqMgrGetSeqEntryForData (pbsp);
  if (psep == NULL) return;

  rnaseq = GetSequenceByFeature (mrna);
  if (rnaseq == NULL) return;
  len = (Int4) StringLen (rnaseq);

  bs = BSNew (len + 2);
  if (bs == NULL) return;
  BSWrite (bs, (VoidPtr) rnaseq, len);
  MemFree (rnaseq);

  mbsp = BioseqNew ();
  if (mbsp == NULL) return;
  mbsp->repr = Seq_repr_raw;
  mbsp->mol = Seq_mol_rna;
  mbsp->seq_data_type = Seq_code_iupacna;
  mbsp->seq_data = (SeqDataPtr) bs;
  mbsp->length = BSLen (bs);
  BioseqPack (mbsp);

  if (useLocusTag) {
    mbsp->id = MakeIdFromLocusTag (mrna);
    if (mbsp->id == NULL) {
      mbsp->id = MakeNewProteinSeqIdEx (mrna->location, NULL, NULL, ctrp);
    }
  } else {
    mbsp->id = MakeNewProteinSeqIdEx (mrna->location, NULL, NULL, ctrp);
  }
  CheckSeqLocForPartial (mrna->location, &partial5, &partial3);
  SeqMgrAddToBioseqIndex (mbsp);

  msep = SeqEntryNew ();
  if (msep == NULL) return;
  msep->choice = 1;
  msep->data.ptrvalue = (Pointer) mbsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mbsp, msep);

  mip = MolInfoNew ();
  if (mip == NULL) return;
  mip->biomol = MOLECULE_TYPE_MRNA;
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (msep, Seq_descr_molinfo);
  if (vnp != NULL) {
    vnp->data.ptrvalue = (Pointer) mip;
  }

  SetSeqFeatProduct (mrna, mbsp);

  AddSeqEntryToSeqEntry (psep, msep, FALSE);
}

static void NPStoGPS (IteM i, Boolean useLocusTag, Boolean useFeatureIDs)

{
  BaseFormPtr        bfp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Int2               ctr = 1;
  SeqMgrFeatContext  fcontext, mcontext;
  SeqFeatPtr         mrna, sfp, lastsfp;
  SeqEntryPtr        old, sep, top;
  CharPtr            organism;
  OrgRefPtr          orp;
  Uint2              parenttype;
  Pointer            parentptr;
  SeqAnnotPtr        sap, annot, lastsap;
  SeqDescrPtr        sdp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  top = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (top == NULL || top->choice != 2) return;
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
    sep = bssp->seq_set;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
    }
  }
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_nuc_prot) return;

  bsp = FindNucBioseq (top);
  if (bsp == NULL) {
    Message (MSG_OK, "Unable to find nucleotide Bioseq");
    return;
  }

  if (useFeatureIDs) {
    if (AmbiguousCdsMrnaByFeatureID(bsp)) {
      Message (MSG_OK, "Unable to proceed because of unresolved CDS-mRNA ambiguity");
      return;
    }
  } else if (AmbiguousCdsMrnaByOverlap (bsp)) {
    Message (MSG_OK, "Unable to proceed because of unresolved CDS-mRNA ambiguity");
    return;
  }

  GetSeqEntryParent (top, &parentptr, &parenttype);

  bssp->_class = BioseqseqSet_class_gen_prod_set;

  sap = bssp->annot;
  bssp->annot = NULL;
  annot = bsp->annot;
  if (annot == NULL) {
    bsp->annot = sap;
  } else if (sap->next == NULL && annot->next == NULL &&
             sap->type == 1 && annot->type == 1 &&
             sap->data != NULL && annot->data != NULL) {
    lastsfp = (SeqFeatPtr) annot->data;
    while (lastsfp->next != NULL) {
      lastsfp = lastsfp->next;
    }
    sfp = (SeqFeatPtr) sap->data;
    lastsfp->next = sfp;
    sap->data = NULL;
    SeqAnnotFree (sap);
  } else {
    lastsap = bsp->annot;
    while (lastsap->next != NULL) {
      lastsap = lastsap->next;
    }
    lastsap->next = sap;
  }

  old = SeqEntrySetScope (top);

  ctr = (Int2) VisitBioseqsInSep (top, NULL, NULL) + 1;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    if (useFeatureIDs) {
      mrna = GetmRNAByFeatureID (sfp);
    } else {
      mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0,
                                          NULL, CHECK_INTERVALS, &mcontext);
    }
    if (mrna != NULL) {
      InstantiateMrnaIntoProt (sfp, mrna, &ctr, useLocusTag);
    }    
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  SeqMgrLinkSeqEntry (top, parenttype, parentptr);

  SeqEntrySetScope (old);

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  /* need to reindex to get mRNA and CDS features from cDNA and protein */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);
  VisitSetsInSet (bssp, NULL, SqnMakeNucProtCDS);

  /* need to reindex before copying genes, instantiating protein titles */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);

  VisitFeaturesInSep (top, NULL, SqnCopyGene);

  /* need to reindex before instantiating mRNA titles */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);

  organism = NULL;
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_source) continue;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) continue;
    orp = biop->org;
    if (orp == NULL) continue;
    if (! StringHasNoText (orp->taxname)) {
      organism = orp->taxname;
    }
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &mcontext);
  while (sfp != NULL) {
    SqnAddMrnaTitles (sfp->product, organism);
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &mcontext);
  }

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  move_cds (top);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

static void NPStoGPSLocusTag (IteM i)

{
  NPStoGPS (i, TRUE, FALSE);
}

static void NPStoGPSArbitrary (IteM i)

{
  NPStoGPS (i, FALSE, FALSE);
}

static void NPStoGPSLocusTagFeatID (IteM i)

{
  NPStoGPS (i, TRUE, TRUE);
}

static void NPStoGPSArbitraryFeatID (IteM i)

{
  NPStoGPS (i, FALSE, TRUE);
}

static void MakePhrap (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  SeqAnnotPtr  sap;
  SeqEntryPtr  sep, nsep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;
  bsp = (BioseqPtr) nsep->data.ptrvalue;
  if (bsp == NULL) return;
  sap = PhrapGraphForContig (bsp);
  if (sap == NULL) return;
  sap->next = bsp->annot;
  bsp->annot = sap;
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct marksegfeatsdata {
  SeqLocPtr seg_loc;
  CdsDelStructPtr   cdsDelInfoPtr;
} MarkSegFeatData, PNTR MarkSegFeatPtr;
  
static Boolean LIBCALLBACK MarkSegmentFeatures (SeqFeatPtr sfp,
					       SeqMgrFeatContextPtr scontext)
{
  MarkSegFeatPtr msfp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_CDS) return TRUE;
  msfp = scontext->userdata;
  if (msfp == NULL || msfp->seg_loc == NULL
    || msfp->cdsDelInfoPtr == NULL
    || sfp->location == NULL)
  {
    return TRUE;
  }
  if (SeqLocAinB (sfp->location, msfp->seg_loc) >=0)
  {
    sfp->idx.deleteme = TRUE;
    msfp->cdsDelInfoPtr->dirtyFlag = TRUE;
  }
  return TRUE;
} 

/*=====================================================================*/
/*                                                                     */
/* RemoveCDS_Callback () - Called for each segment explored in         */
/*                         RemoveNthCDSfromSegSets_Callback(), this    */
/*                         removes the first CDS of a segment.         */
/*                                                                     */
/*=====================================================================*/

static Boolean LIBCALLBACK RemoveCDS_Callback (SeqLocPtr slp,
					       SeqMgrSegmentContextPtr scontext)
{
  BioseqPtr         bsp;
  MarkSegFeatData   msf;

  /* If not the Nth segment, skip */

  msf.cdsDelInfoPtr = (CdsDelStructPtr) scontext->userdata;
  if (scontext->index != msf.cdsDelInfoPtr->number)
    return TRUE;

  /* Else, delete its CDS */

  else {

    /* Mark the coding regions for this segment */
    bsp = BioseqFindFromSeqLoc (slp);
    msf.seg_loc = slp;
    SeqMgrExploreFeatures (bsp, (Pointer) &msf, MarkSegmentFeatures, NULL, NULL, NULL);
    
    /* We only want to remove a CDS from the Nth   */
    /* segment so stop the exploring at this point */
    
    return FALSE;
  }
}

/*=====================================================================*/
/*                                                                     */
/* RemoveNthCDS() - Remove the CDS from the first segment of a given   */
/*                  segmented sequence.                                */
/*                                                                     */
/*=====================================================================*/

static void RemoveNthCDS (SeqEntryPtr      sep,
			  CdsDelStructPtr  cdsDelInfoPtr)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          segCount;

  /* If the SeqEntry is a Bioseq, then remove */
  /* the Nth CDS from each segment.           */

  if (IS_Bioseq(sep))
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp->repr == Seq_repr_seg) 
	segCount = SeqMgrExploreSegments (bsp, cdsDelInfoPtr,
					  RemoveCDS_Callback);
    }

  /* If we have a Bioseq set, then recursively call */
  /* this function for each member of the set.      */

  else
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
    	RemoveNthCDS (sep, cdsDelInfoPtr);
    }

}

/*=====================================================================*/
/*                                                                     */
/* RemoveNthCDS_AcceptCallback () --                                   */
/*                                                                     */
/*=====================================================================*/

static void RemoveNthCDS_AcceptCallback (ButtoN b)
{
  CdsDelFormPtr cdfp;
  SeqEntryPtr   sep;
  Uint2         entityID;
  CdsDelStruct  cdsDelInfo;

  /* Make sure that conditions are right */

  cdfp = GetObjectExtra (b);
  if (NULL == cdfp) {
    Remove (ParentWindow (b));
    Update ();
    return;
  }

  if (cdfp->number == 0) {
    Remove (ParentWindow (b));
    MemFree (cdfp);
    Update ();
    return;
  }

  /* Get the top level SeqEntry */

  sep = GetTopSeqEntryForEntityID (cdfp->input_entityID);
  if (sep == NULL) {
    Remove (ParentWindow (b));
    MemFree (cdfp);
    Update ();
    return;
  }

  /* Remove CDS from segments */

  WatchCursor ();
  Update ();

  cdsDelInfo.dirtyFlag = FALSE;
  cdsDelInfo.number    = cdfp->number;

  RemoveNthCDS (sep, &cdsDelInfo);

  /* If any CDS were marked for deletion */
  /* then delete them and update.        */

  if (TRUE == cdsDelInfo.dirtyFlag)
    {
      entityID = ObjMgrGetEntityIDForChoice (sep);
      DeleteMarkedObjects (entityID, 0, NULL);
      Update ();
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }

  /* Update the display */

  Remove (ParentWindow (b));
  MemFree (cdfp);
  ArrowCursor ();
  Update ();
}

/*=====================================================================*/
/*                                                                     */
/* GetN_Callback () --                                                 */
/*                                                                     */
/*=====================================================================*/

static void GetN_Callback (TexT nText)

{
  CdsDelFormPtr cdfp;
  CharPtr       nStr;

  cdfp = (CdsDelFormPtr) GetObjectExtra (nText);
  if (cdfp == NULL)
    return;

  nStr = SaveStringFromText (nText);
  cdfp->number = atoi (nStr);

}

/*=====================================================================*/
/*                                                                     */
/* CreateGetNWindow () --                                              */
/*                                                                     */
/*=====================================================================*/

static ForM CreateGetNWindow (CdsDelFormPtr cdfp)
{
  /*
  FindFormPtr        ffp;
  GrouP              q = NULL;
  */

  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              j;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  TexT               nthText;

  w = ModalWindow (-50, -30, -10, -10, StdCloseWindowProc);
  cdfp->form = (ForM) w;

#ifndef WIN_MAC
  CreateStdEditorFormMenus (w);
#endif

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
  }
  
  j = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (j, 10, 10);
  
  g = HiddenGroup (j, 2, 0, NULL);
  StaticPrompt (g, "CDS to Delete", 0, dialogTextHeight, programFont, 'l');
  nthText = DialogText (g, "", 6, (TxtActnProc) GetN_Callback);
  SetObjectExtra (nthText, cdfp, NULL);
  
  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  b = DefaultButton (c, "Accept", RemoveNthCDS_AcceptCallback);
  SetObjectExtra (b, cdfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) j, (HANDLE) c, NULL);
  
  RealizeWindow (w);
  
  return (ForM) w;
}

/*=====================================================================*/
/*                                                                     */
/* RemoveNthCDSFromSegSets_Callback() - Remove the CDS from the first  */
/*                                      segment of all Bioseqs.        */
/*                                                                     */
/*=====================================================================*/

static void RemoveNthCDSFromSegSets_Callback (IteM i)
{
  ForM          w;
  CdsDelFormPtr cdfp;
  BaseFormPtr   bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  cdfp = (CdsDelFormPtr) MemNew (sizeof (CdsDelForm));

  /* Find out which segments to delete the CDS from */

  cdfp->input_entityID = bfp->input_entityID;
  w = CreateGetNWindow (cdfp);
  Show (w);
  Select (w);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* RemoveSeqIdName_FeatureCallback () --                               */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK RemoveSeqIdName_FeatureCallback
                               (SeqFeatPtr sfp,
				SeqMgrFeatContextPtr fcontext)
{
  SeqIdPtr     sip;
  TextSeqIdPtr tsip;

  /* If we don't have a product then */
  /* go on to the next seq feat.     */

  if (sfp->product == NULL)
    return TRUE;
  
  /* Otherwise check the id name field */

  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
      switch (sip->choice) {
        case SEQID_LOCAL :
	  break;
        case SEQID_GENBANK :
        case SEQID_EMBL :
        case SEQID_DDBJ :
        case SEQID_OTHER :
        case SEQID_TPG :
        case SEQID_TPE :
        case SEQID_TPD :
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if ((tsip != NULL) &&
	      (tsip->accession == NULL) &&
	      (! StringHasNoText (tsip->name)) &&
	      (IsNuclAcc (tsip->name)))
	    tsip->name[0] = '\0';
          break;
        default :
          break;
      }
    }

  /* Return TRUE to continue on to the next feature */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* RemoveSeqIdName_BioseqCallback () --                                */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK RemoveSeqIdName_BioseqCallback (BioseqPtr bsp,
					 SeqMgrBioseqContextPtr bcontext)
{
  TextSeqIdPtr tsip;

  /* If it is a protein Bioseq then */
  /* check the ID name field.       */

  if (ISA_aa (bsp->mol))
    {
      if (bsp->id != NULL) {
	switch (bsp->id->choice)
	  {
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_OTHER :
          case SEQID_TPG :
	  case SEQID_TPE :
	  case SEQID_TPD :
	    tsip = (TextSeqIdPtr) bsp->id->data.ptrvalue;
	    if ((tsip != NULL) &&
		(tsip->accession == NULL) &&
		(! StringHasNoText (tsip->name)) &&
		(IsNuclAcc (tsip->name)))
	      {
		tsip->name[0] = '\0';
	      }
	    break;
	  default :
	    break;
	  }
      }
    }

  /* Else for nucl bioseqs check the name */
  /* field of the protein sfps.           */

  else
    SeqMgrExploreFeatures (bsp, NULL, RemoveSeqIdName_FeatureCallback,
			   NULL, NULL, NULL);

  /* Return TRUE to continue on to the next Bioseq */

  return TRUE;
}

/*=====================================================================*/
/*                                                                     */
/* RemoveSeqIdName_Callback() - Remove the contents of the seq-id name */
/*                              field in the cases where there is no   */
/*                              corresponding GI number.               */
/*                                                                     */
/*=====================================================================*/

static void RemoveSeqIdName_Callback (IteM i)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;

  SeqMgrExploreBioseqs (bfp->input_entityID, NULL, bfp,
			RemoveSeqIdName_BioseqCallback,
			TRUE, TRUE, TRUE);

  /* Force an update and redraw */

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();

  /* Return successfully */

  return;

}

static void ReverseNameOrderInAuthor (AuthorPtr pAuthor)
{
  NameStdPtr pNameStandard;
  NameStdPtr pNewName;
  CharPtr	first_name;
  CharPtr	middle_initial;
  CharPtr	last_name;
  CharPtr	suffix;
  Char		str[128];
  CharPtr	pTmp;

  if (pAuthor == NULL)
    return;
  else if(pAuthor->name->choice != 2)
    return;
  pNameStandard = pAuthor->name->data;
  if (pNameStandard != NULL)
  {
    pTmp = NameStdPtrToAuthorSpreadsheetString(pNameStandard);
    first_name = ExtractTagListColumn(pTmp, 2);
    middle_initial = ExtractTagListColumn(pTmp, 1);
    last_name = ExtractTagListColumn(pTmp, 0);
    suffix = ExtractTagListColumn(pTmp, 3);
    sprintf(str, "%s\t%s\t%s\t%s\n",
	first_name, middle_initial, last_name, suffix);
    MemFree(first_name);
    MemFree(middle_initial);
    MemFree(last_name);
    MemFree(suffix);
    pNewName = AuthorSpreadsheetStringToNameStdPtr(str);
    NameStdFree(pNameStandard);
    pAuthor->name->data = pNewName;
  }

}

extern void 
FixCapitalizationInTitle 
(CharPtr PNTR pTitle,
 Boolean      first_is_upper,
 ValNodePtr   org_names)
{
  if (pTitle == NULL) return;
  ResetCapitalization (first_is_upper, *pTitle);
  FixAbbreviationsInElement (pTitle);
  FixOrgNamesInString (*pTitle, org_names);
}


static void FixCapitalizationInAuthor (AuthorPtr pAuthor)
{
  NameStdPtr pNameStandard;
  CharPtr    cp;
  
  if (pAuthor == NULL)
    return;
  else if(pAuthor->name->choice != 2)
    return;
  pNameStandard = pAuthor->name->data;
  if (pNameStandard != NULL)
  {
    FixCapitalizationInElement (&(pNameStandard->names[0]), FALSE, FALSE, TRUE);
    FixCapitalizationInElement (&(pNameStandard->names[1]), FALSE, FALSE, FALSE);
    /* Set initials to all caps */
    for (cp = pNameStandard->names[4]; cp != NULL && *cp != 0; cp++)
    {
      *cp = toupper (*cp);
    }
  }
}

static void StripSuffixFromAuthor (AuthorPtr pAuthor)
{
  NameStdPtr pNameStandard;
  
  if (pAuthor == NULL)
    return;
  else if(pAuthor->name->choice != 2)
    return;
  pNameStandard = pAuthor->name->data;
  if (pNameStandard != NULL && pNameStandard->names[5] != NULL)
  {
    pNameStandard->names[5][0] = 0;
  }
}

static void TruncateAuthorMiddleInitials (AuthorPtr pAuthor)
{
  NameStdPtr pNameStandard;
  CharPtr cp;
  
  if (pAuthor == NULL)
    return;
  else if(pAuthor->name->choice != 2)
    return;
  pNameStandard = pAuthor->name->data;
  if (pNameStandard != NULL)
  {
    cp = StringChr (pNameStandard->names[4], '.');
    if (cp == NULL || StringChr (cp + 1, '.') == NULL) {
      if (StringLen (pNameStandard->names[4]) > 3) 
      {
        pNameStandard->names[4][3] = 0;
        pNameStandard->names[4][2] = '.';
      }
    } else if (StringLen (pNameStandard->names[4]) > 4) {
        pNameStandard->names[4][4] = 0;
        pNameStandard->names[4][3] = '.';
    }
  }
}

static void ChangeAuthorLastNameToConsortium (AuthorPtr pAuthor, StringConstraintXPtr author_scp)
{
  NameStdPtr pNameStandard;
  CharPtr    consortium;
  Int4       len;
  
  if (pAuthor == NULL || pAuthor->name == NULL || pAuthor->name->choice == 5)
  {
  	return;
  }
  if (pAuthor->name->choice == 2)
  {
    pNameStandard = pAuthor->name->data;
    if (pNameStandard != NULL)
    {
      len = StringLen (pNameStandard->names[0]);
      if (len > 0)
      {
        if (DoesStringMatchConstraintX (pNameStandard->names[0], author_scp))
        {
    	    consortium = (CharPtr) MemNew ((len + 1) * sizeof (Char));
    	    if (consortium != NULL)
  	      {
  	        StringCpy (consortium, pNameStandard->names[0]);
  	        NameStdFree (pNameStandard);
  	        pAuthor->name->choice = 5;
  	        pAuthor->name->data = consortium;
  	      }
  	    }
      }
    }  
  }
}


#define FIX_AUTHOR_NAME_ORDER	          1
#define FIX_PUB_AUTHOR_CAPITALIZATION	  2
#define FIX_PUB_TITLE_CAPITALIZATION	  4
#define FIX_PUB_AFFIL_CAPITALIZATION	  8
#define FIX_SELECTED			              16
#define FIX_ALL                         32
#define STRIP_AUTHOR_SUFFIX             64
#define FIX_PUB_SWAP_NAME_CONSORTIUM    128
#define TRUNCATE_AUTHOR_MIDDLE_INITIALS 256

typedef struct fixpubdesc 
{
  Int4                iType;
  StringConstraintXPtr author_scp;
  ValNodePtr          org_names;
} FixPubdescData, PNTR FixPubdescPtr;

static void FixPubdesc (PubdescPtr pdp, Pointer userdata)
{
  ValNodePtr	  vnp;
  AuthListPtr 	alp;
  CitArtPtr     cap;
  CitBookPtr    cbp;
  CitGenPtr     cgp;
  CitPatPtr     cpp;
  CitSubPtr     csp;
  ValNodePtr	  names;
  AuthorPtr	    ap;
  CharPtr PNTR	title;
  Int4	        iType;
  AffilPtr	    affil;
  FixPubdescPtr fpp;
  
  fpp = (FixPubdescPtr) userdata;
  if (fpp == NULL)
  {
    return;
  }
  
  iType = fpp->iType;

  if (pdp == NULL) return;
  /* search for PUB_PMid or PUB_Muid - if found, do not edit this pub */
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  }
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) continue;
    if (vnp->data.ptrvalue == NULL) continue;
    alp = NULL;
    title = NULL;
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;   
        if( !(cgp->cit != NULL && StringCmp(cgp->cit, "unpublished") != 0)
  	      && !((iType & FIX_SELECTED) || (iType & FIX_ALL)) ) continue;
        alp = cgp->authors;
        title = &(cgp->title);
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        alp = csp->authors;
        title = NULL;
        break;
      case PUB_Article :
        if ( !((iType & FIX_SELECTED) || (iType & FIX_ALL)) ) continue;
        cap = (CitArtPtr) vnp->data.ptrvalue;
        alp = cap->authors;
        if(cap->title != NULL) {
          title = (CharPtr PNTR) &(cap->title->data.ptrvalue);
        }
        break;
      case PUB_Book :
      case PUB_Man :
        if ( !((iType & FIX_SELECTED) || (iType & FIX_ALL)) ) continue;
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        alp = cbp->authors;
        if(cbp->title != NULL) {
          title = (CharPtr PNTR) &(cbp->title->data.ptrvalue);
        }
        break;
      case PUB_Patent :
        if ( !((iType & FIX_SELECTED) || (iType & FIX_ALL)) ) continue;
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        alp = cpp->authors;
        title = &(cpp->title);
        break;
      default :
        break;
    }

    if(iType & FIX_PUB_TITLE_CAPITALIZATION ) {
      FixCapitalizationInTitle (title, TRUE, fpp->org_names);
    }

    if(alp == NULL) continue;

    /* Loop through author list */
    for (names = alp->names; names != NULL; names = names->next) { 
      ap = names->data.ptrvalue;
      if(iType & FIX_PUB_AUTHOR_CAPITALIZATION ) {
        FixCapitalizationInAuthor (ap);
      }
      if(iType & FIX_AUTHOR_NAME_ORDER ) {
        ReverseNameOrderInAuthor (ap);
      }
      if (iType & TRUNCATE_AUTHOR_MIDDLE_INITIALS)
      {
        TruncateAuthorMiddleInitials (ap);
      }
      if (iType & STRIP_AUTHOR_SUFFIX) {
        StripSuffixFromAuthor (ap);
      }
      if (iType & FIX_PUB_SWAP_NAME_CONSORTIUM)
      {
      	ChangeAuthorLastNameToConsortium (ap, fpp->author_scp);
      }
    }

    if (iType & FIX_PUB_AFFIL_CAPITALIZATION ) {
      affil = alp->affil;
      if (affil == NULL) continue;
      FixCapitalizationInElement (&(affil->affil), TRUE, TRUE, FALSE);
      FixCapitalizationInElement (&(affil->div), TRUE, TRUE, FALSE);
      FixCapitalizationInElement (&(affil->city), FALSE, TRUE, FALSE);

      /* special handling for states */
      if (affil->sub != NULL && StringLen (affil->sub) == 2
	      && isalpha((Int4)(affil->sub[0]))	&& isalpha((Int4)(affil->sub[1])))
      {
        affil->sub[0] = toupper(affil->sub[0]);
        affil->sub[1] = toupper(affil->sub[1]);
      } else {
        FixCapitalizationInElement (&(affil->sub), FALSE, TRUE, FALSE);
      }
      FixCapitalizationInElement (&(affil->country), TRUE, TRUE, FALSE);
      FixCapitalizationInElement (&(affil->street), FALSE, TRUE, FALSE);
    }
  }
}


extern void GetOrgNamesInRecordCallback (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr PNTR org_names;
  
  if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)
      || userdata == NULL)
  {
    return;
  }
  
  org_names = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (org_names, 0, biop->org->taxname);
}


/* This function is used to apply fixes to citations */
static void FixPubs (Uint2 entityID, Int4 iType, StringConstraintXPtr author_scp)
{
  SeqEntryPtr  sep;
  SelStructPtr	sel;
  SeqFeatPtr	sfp;
  SeqMgrFeatContext fcontext;
  SeqDescPtr	sdp;
  SeqMgrDescContext dcontext;
  PubdescPtr	pdp;
  FixPubdescData fpd;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (iType & FIX_ALL)
  {
    sel = NULL;
  }
  else
  {
    sel = ObjMgrGetSelected ();
  }
  
  fpd.iType = iType;
  fpd.author_scp = author_scp;
  
  fpd.org_names = NULL;
  if (iType & FIX_PUB_TITLE_CAPITALIZATION)
  {
    VisitBioSourcesInSep (sep, &(fpd.org_names), GetOrgNamesInRecordCallback);
  }
  
  if(sel == NULL)
  {
    VisitPubdescsInSep (sep, (Pointer)&fpd, FixPubdesc);
  }
  else
  {
    fpd.iType |= FIX_SELECTED;
    while( sel != NULL )
    {
      pdp = NULL;
      if(sel->entityID == entityID)
      {
        if(sel->itemtype == OBJ_SEQFEAT) 
        {
          sfp = SeqMgrGetDesiredFeature (entityID, NULL, sel->itemID, 0, NULL, &fcontext);
          if(sfp != NULL && sfp->data.choice == SEQFEAT_PUB)
          {
            pdp = sfp->data.value.ptrvalue;
          }
        }
        else if(sel->itemtype == OBJ_SEQDESC)
        {
          sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, sel->itemID, 0, NULL, &dcontext);
          if(sdp != NULL && sdp->choice == Seq_descr_pub)
          {
            pdp = sdp->data.ptrvalue;
          }
        }
      } 
      if (pdp != NULL) FixPubdesc (pdp, (Pointer) &fpd);
      sel = sel->next;      
    }
  }
  
  fpd.org_names = ValNodeFree (fpd.org_names);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  Update ();
}

static void FixPubsMenuItem (IteM i, Int4 iType)
{
  BaseFormPtr bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  FixPubs (bfp->input_entityID, iType, NULL);
}

static void FixCapitalizationToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  FixPubs (bfp->input_entityID, FIX_PUB_AUTHOR_CAPITALIZATION
              | FIX_PUB_TITLE_CAPITALIZATION
              | FIX_PUB_AFFIL_CAPITALIZATION
              | FIX_ALL,
              NULL);
}

static void FixNameOrder (IteM i)
{
  FixPubsMenuItem (i, FIX_AUTHOR_NAME_ORDER);
}

static void FixAllPubCapitalization (IteM i)
{
  FixPubsMenuItem (i, FIX_PUB_AUTHOR_CAPITALIZATION | FIX_PUB_TITLE_CAPITALIZATION | FIX_PUB_AFFIL_CAPITALIZATION);
}

static void FixAuthorCapitalization (IteM i)
{
  FixPubsMenuItem (i, FIX_PUB_AUTHOR_CAPITALIZATION);
}

static void TruncateAuthorMiddle (IteM i)
{
  FixPubsMenuItem (i, TRUNCATE_AUTHOR_MIDDLE_INITIALS);
}

static void FixTitleCapitalization (IteM i)
{
  FixPubsMenuItem (i, FIX_PUB_TITLE_CAPITALIZATION);
}

static void FixAffiliationCapitalization (IteM i)
{
  FixPubsMenuItem (i, FIX_PUB_AFFIL_CAPITALIZATION);
}

static void StripAuthorSuffix (IteM i)
{
  FixPubsMenuItem (i, STRIP_AUTHOR_SUFFIX);
}

static void ChangeAllAuthorNameToConsortium (IteM i)
{
  FixPubsMenuItem (i, FIX_PUB_SWAP_NAME_CONSORTIUM);
}

static void ChangeAuthorNameWithConsortiumToConsortium (IteM i)
{
  BaseFormPtr bfp;
  StringConstraintXPtr scp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  
  scp = (StringConstraintXPtr) MemNew (sizeof (StringConstraintData));
  if (scp == NULL)
  {
    return;
  }
  
  scp->match_text = StringSave ("consortium");
  scp->match_location = eStringConstraintContains;
  scp->insensitive = TRUE;
  scp->whole_word = TRUE;
  scp->not_present = FALSE;   
  
  FixPubs (bfp->input_entityID, FIX_PUB_SWAP_NAME_CONSORTIUM, scp);
  scp = StringConstraintXFree (scp);
}

typedef struct authornamestringconstraintform 
{
  FORM_MESSAGE_BLOCK
  DialoG author_name_constraint_dlg;
} AuthorNameStringConstraintFormData, PNTR AuthorNameStringConstraintFormPtr;

static void DoChangeAuthorNameToConsortiumWithConstraint (ButtoN b)
{
  AuthorNameStringConstraintFormPtr frm;
  StringConstraintXPtr               scp;
  
  frm = (AuthorNameStringConstraintFormPtr) GetObjectExtra (b);
  if (frm == NULL)
  {
    return;
  }
  
  scp = (StringConstraintXPtr) DialogToPointer (frm->author_name_constraint_dlg);
  
  Hide (frm->form);
  WatchCursor ();
  Update ();
  
  FixPubs (frm->input_entityID, FIX_PUB_SWAP_NAME_CONSORTIUM, scp);
  scp = StringConstraintXFree (scp);
  Remove ((WindoW)frm->form);  
  ArrowCursor ();
  Update ();
}

static void ChangeAuthorNameToConsortiumWithConstraint (IteM i)
{
  AuthorNameStringConstraintFormPtr frm;
  BaseFormPtr bfp;
  WindoW      w;
  GrouP       h, c;
  ButtoN      b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  
  frm = (AuthorNameStringConstraintFormPtr) MemNew (sizeof (AuthorNameStringConstraintFormData));
  if (frm == NULL)
  {
    return;
  }
  
  w = FixedWindow (-50, -33, -10, -10, "Change Author Last Name to Consortium", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  frm->author_name_constraint_dlg =  StringConstraintDialogX (h, "Where author last name", TRUE);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  b = PushButton (c, "Accept", DoChangeAuthorNameToConsortiumWithConstraint);
  SetObjectExtra (b, frm, NULL);
  
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->author_name_constraint_dlg,
                              (HANDLE) c, 
                              NULL);
                                
  Show (w);  
}

static void FixProductCapitalizationCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ProtRefPtr prp;
  RnaRefPtr  rrp;
  ValNodePtr vnp;
  CharPtr    product_name;
  
  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_PROT && sfp->data.choice != SEQFEAT_RNA) return;
  
  if (sfp->data.choice == SEQFEAT_PROT)
  {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp != NULL)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        product_name = (CharPtr)(vnp->data.ptrvalue);
        FixCapitalizationInTitle (&product_name, FALSE, NULL);
        vnp->data.ptrvalue = product_name;
      }
    }
  }
  else if (sfp->data.choice == SEQFEAT_RNA)
  {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 1)
    {
      product_name = (CharPtr)(rrp->ext.value.ptrvalue);
      FixCapitalizationInTitle (&product_name, FALSE, NULL);
      rrp->ext.value.ptrvalue = product_name;
    }
  }
}

static void FixProductCapitalization (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  SelStructPtr      sel;
  Boolean           any_found = FALSE;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  RnaRefPtr         rrp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  sel = ObjMgrGetSelected ();
  if(sel == NULL)
  {
    VisitFeaturesInSep (sep, NULL, FixProductCapitalizationCallback);
  }
  else
  {
    while( sel != NULL )
    {
      if(sel->entityID == bfp->input_entityID && sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if(sfp != NULL)
        {
          if (sfp->data.choice == SEQFEAT_PROT)
          {
          	FixProductCapitalizationCallback (sfp, NULL);
          	any_found = TRUE;
          }
          else if (sfp->data.choice == SEQFEAT_RNA)
          {
            rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          	if (rrp != NULL && rrp->ext.choice == 1)
          	{
          	  FixProductCapitalizationCallback (sfp, NULL);
          	  any_found = TRUE;
          	}
          }
        }
      } 
      sel = sel->next;      
    }
    if (!any_found)
    {
      VisitFeaturesInSep (sep, NULL, FixProductCapitalizationCallback);      
    }
  }
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();	
}

static void FixDeltaSeqDataLenCallback (BioseqPtr bsp, Pointer userdata)
{
  DeltaSeqPtr dsp;
  SeqLocPtr   slp;
  SeqLitPtr   slip;
  Int4        bsp_len = 0;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4)
  {
    return;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    if (dsp->choice == 2)
    {
      slip = (SeqLitPtr) dsp->data.ptrvalue;
      bsp_len += slip->length;
    }
    else if (dsp->choice == 1)
    {
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      bsp_len += SeqLocLen (slp);
    }
    dsp = dsp->next;
  }
  bsp->length = bsp_len;
}

static void FixDeltaSeqDataLen (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, FixDeltaSeqDataLenCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();	  
}

static void InstantiateGapFeatCallback (BioseqPtr bsp, Pointer userdata)

{
  Char        buf [32];
  Int4        currpos = 0;
  IntFuzzPtr  fuzz;
  ImpFeatPtr  ifp;
  SeqLitPtr   litp;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  SeqLocPtr   slp;
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  /* suppress on far delta contigs for now */
  if (! DeltaLitOnly (bsp)) return;

  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      if (slp == NULL) continue;
      currpos += SeqLocLen (slp);
    }
    if (vnp->choice == 2) {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL) continue;
      if (IsDeltaSeqGap (vnp) && litp->length > 0) {
        ifp = ImpFeatNew ();
        if (ifp == NULL) continue;
        ifp->key = StringSave ("gap");
        fuzz = litp->fuzz;
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, NULL);
        if (sfp == NULL) continue;
        sfp->data.choice = SEQFEAT_IMP;
        sfp->data.value.ptrvalue = (Pointer) ifp;
        if (IsDeltaSeqUnknownGap (vnp)) {
          AddQualifierToFeature (sfp, "estimated_length", "unknown");
          sfp->location = SeqLocFree (sfp->location);
          sfp->location = AddIntervalToLocation (NULL, sip, currpos, currpos + litp->length - 1, FALSE, FALSE);
        } else {
          sprintf (buf, "%ld", (long) litp->length);
          AddQualifierToFeature (sfp, "estimated_length", buf);
          sfp->location = SeqLocFree (sfp->location);
          sfp->location = AddIntervalToLocation (NULL, sip, currpos, currpos + litp->length - 1, FALSE, FALSE);
        }
      }
      currpos += litp->length;
    }
  }
}

static void InstantiateGapFeatures (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, InstantiateGapFeatCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();	  
}

static void RemoveGapFeatCallback (SeqFeatPtr sfp, Pointer userdata)

{
  GBQualPtr  gbq;

  if (sfp == NULL) return;
  if (sfp->idx.subtype != FEATDEF_gap) return;
  if (StringDoesHaveText (sfp->comment)) return;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "estimated_length") != 0) return;
  }
  sfp->idx.deleteme = TRUE;
}

static void RemoveUnnecessaryGapFeatures (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitFeaturesInSep (sep, NULL, RemoveGapFeatCallback);

  ObjMgrSelect (0, 0, 0, 0, NULL);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();	  
}


static void RenormalizeNucProtSetsMenuItem (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  /* Get the current data */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor();
  Update();
  RenormalizeNucProtSets (sep, TRUE);   	
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();

  /* Return successfully */

  return;
}


static void AddAllDbxrefsBioseqCallback (BioseqPtr bsp, Pointer data)
{
  AddAllDbxrefsToBioseq (bsp);
}


static void AddAllDbxrefsMenuItem (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
 
  /* Get the current data */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitBioseqsInSep (sep, NULL, AddAllDbxrefsBioseqCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern void LookupAllPubs (IteM i);
extern void ResolveExistingLocalIDs (IteM i);
extern void PromoteAlignsToBestIDProc (IteM i);
extern void SetSourceFocus (IteM i);
extern void ClearSourceFocus (IteM i);
extern void ExtraAccToHistByPos (IteM i);
extern void ClearCdsProducts (IteM i);
extern void ClearMrnaProducts (IteM i);
extern void MakeRedundantGPSwithProp (IteM i);
extern void MakeRedundantGPSnoProp (IteM i);
extern void MakeRedundantGPSjustXref (IteM i);
extern void FuseSlpJoins (IteM i);
extern void CreateSeqHistTPA (IteM i);
extern void CreateSeqHistTPADetailed (IteM i);
extern void CreateSeqHistDelta (IteM i);
extern void AddGlobalCodeBreak (IteM i);
extern void FixCdsAfterPropagate (IteM i);
extern void ParseTrinomial (IteM i);
extern void ReprocessPeptideProducts (IteM i);
extern void ReprocessmRNAProducts (IteM i);
extern void EditSeqEndsProc (IteM i);

static void NewAligmentViewer (IteM i)
{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  ValNodePtr   seq_annot_list;
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp = NULL;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  
  bsp = GetBioseqGivenIDs (bfp->input_entityID, 
                           bfp->input_itemID,
                           bfp->input_itemtype);
  if (bsp == NULL) return;
  seq_annot_list = FindAlignSeqAnnotsForBioseq (bsp);
  if (seq_annot_list != NULL && seq_annot_list->data.ptrvalue != NULL)
  {
    sap = (SeqAnnotPtr) seq_annot_list->data.ptrvalue;
    if (sap->type == 2)
    {
      salp = (SeqAlignPtr) sap->data;
    }
  }
  if (salp == NULL)
  {
    Message (MSG_ERROR, "No alignments found for this Bioseq!");
    return;
  }
  
  OpenNewAlignmentEditor (salp, bfp->input_entityID);

}

static void ShowDiscrepancyReport (IteM i)
{
  BaseFormPtr  bfp;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  
  CreateReportWindow (eReportTypeDiscrepancy);
}

static void ShowOnCallerTool (IteM i)
{
  BaseFormPtr  bfp;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;
  
  CreateReportWindow (eReportTypeOnCaller);
}

static void TestTax3 (IteM i, Int4 taxid, CharPtr name)

{
  BaseFormPtr       bfp;
  OrgNamePtr        onp;
  OrgRefPtr         orp;
  T3DataPtr         tdp;
  T3ErrorPtr        tep;
  Taxon3RequestPtr  t3rq;
  Taxon3ReplyPtr    t3ry;
  T3ReplyPtr        trp;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  t3rq = CreateTaxon3Request (taxid, name, NULL);
  if (t3rq == NULL) return;
  t3ry = Tax3SynchronousQuery (t3rq);
  Taxon3RequestFree (t3rq);
  if (t3ry != NULL) {
    for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
      switch (trp->choice) {
        case T3Reply_error :
          tep = (T3ErrorPtr) trp->data.ptrvalue;
          if (tep != NULL) {
            Message (MSG_POST, "T3 error %s", tep->message);
          }
          break;
        case T3Reply_data :
          tdp = (T3DataPtr) trp->data.ptrvalue;
          if (tdp != NULL) {
            orp = (OrgRefPtr) tdp->org;
            if (orp != NULL) {
              onp = orp->orgname;
              if (onp != NULL) {
                Message (MSG_POST, "T3 data %s code %d", orp->taxname, (int) onp->gcode);
              }
            }
          }
          break;
        default :
          break;
      }
    }
    Taxon3ReplyFree (t3ry);
  }
}

static void TestTax3A (IteM i)

{
  TestTax3 (i, 9606, NULL);
}

static void TestTax3B (IteM i)

{
  TestTax3 (i, 0, "Escherichia coli");
}

static void MakeFarProtRefXref (SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  Uint2              entityID;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prot;
  ProtRefPtr         prp;
  SeqIdPtr           sip;
  SeqFeatXrefPtr     xref;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice == SEQFEAT_PROT) return;
  }

  sip = SeqLocIdForProduct (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqLockById (sip);
  if (bsp == NULL) return;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  if (entityID != 0 && SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  prot = SeqMgrGetBestProteinFeature (bsp, &pcontext);
  if (prot != NULL && prot->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) prot->data.value.ptrvalue;
    if (prp != NULL) {
      xref = SeqFeatXrefNew ();
      if (xref != NULL) {
        xref->data.choice = SEQFEAT_PROT;
        xref->data.value.ptrvalue = AsnIoMemCopy ((Pointer) prp,
                                                  (AsnReadFunc) ProtRefAsnRead,
                                                  (AsnWriteFunc) ProtRefAsnWrite);
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
    }
  }

  BioseqUnlock (bsp);
}

static void MakeRedundantGenomicProtXref (
  IteM i
)

{
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  /* find genomic sequence */
  bsp = FindNucBioseq (sep);
  if (bsp != NULL) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      MakeFarProtRefXref (sfp);
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void MarkOrphanedProteins (BioseqPtr bsp, Pointer userdata)

{
  SeqFeatPtr         cds;
  SeqMgrFeatContext  context;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  cds = SeqMgrGetCDSgivenProduct (bsp, &context);
  if (cds == NULL) {
    bsp->idx.deleteme = TRUE;
  }
}

static void RemoveOrphanedProteins (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, MarkOrphanedProteins);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ViewSortedProteinMenuItem(IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ViewSortedProteins(sep);
}


static void WeightedAlignmentPercentIdentityCallback (SeqAlignPtr salp, Pointer userdata)
{
  ValNodePtr PNTR pct_list;
  Uint2           wp;
  
  if (salp == NULL || userdata == NULL) return;
  pct_list = (ValNodePtr PNTR) userdata;

  wp = WeightedAlignmentPercentIdentity (salp, FALSE);
  ValNodeAddInt (pct_list, 0, wp);
}

static void WeightedAlignmentPercentIdentityItem (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ValNodePtr   pct_list = NULL, vnp;
  FILE         *fp;
  Char         path [PATH_MAX];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitAlignmentsInSep (sep, &pct_list, WeightedAlignmentPercentIdentityCallback);

  if (pct_list == NULL) {
    Message (MSG_OK, "No alignments!");
  } else if (pct_list->next == NULL) {
    Message (MSG_OK, "Weighted Alignment Percent Identity: %d", pct_list->data.intvalue);
  } else {
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      fprintf (fp, "This record contains %d alignments, with the following weighted alignment percent identity scores:\n", ValNodeLen (pct_list));
      for (vnp = pct_list; vnp != NULL; vnp = vnp->next) {
        fprintf (fp, "%d\n", vnp->data.intvalue);
      }
      FileClose (fp);
      LaunchGeneralTextViewer (path, "Weighted Alignment Percent Identity");
      FileRemove (path);
    }    
  }
  pct_list = ValNodeFree (pct_list);
}


static void ApplyBarcodeDbxrefs (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, ApplyBarcodeDbxrefsToBioseq);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update();
}


static void InstantiateMatPeptideProductsMenuItem (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  InstantiateMatPeptideProducts (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update();
}

static void ModPCRCallback (BioSourcePtr biop, Pointer userdata)

{
  ModernizePCRPrimers (biop);
}

static void ModernizePCRs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioSourcesInSep (sep, NULL, ModPCRCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ModRNACallback (SeqFeatPtr sfp, Pointer userdata)

{
  ModernizeRNAFields (sfp);
}

static void ModernizeRNAs (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ModRNACallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ModGeneCallback (SeqFeatPtr sfp, Pointer userdata)

{
  ModernizeGeneFields (sfp);
}

static void ModernizeGenes (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ModGeneCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void GI2dbSTSCallback (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr    bsp;
  Char         buf [64];
  DbtagPtr     dbt;
  ImpFeatPtr   ifp;
  MolInfoPtr   mip;
  ObjectIdPtr  oip;
  SeqDescrPtr  sdp;
  SeqId        si;
  SeqIdPtr     sip;
  Int4         uid;
  ValNodePtr   vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP) return;
  ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
  if (ifp == NULL || StringICmp (ifp->key, "STS") != 0) return;
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringICmp (dbt->db, "GI") != 0) continue;
    oip = dbt->tag;
    if (oip == NULL) continue;
    if (oip->str != NULL) continue;
    uid = oip->id;
    if (uid < 1) continue;
    MemSet ((Pointer) &si, 0, sizeof (SeqId));
    si.choice = SEQID_GI;
    si.data.intvalue = uid;
    bsp = BioseqLockById (&si);
    if (bsp == NULL) continue;
    sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
    if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (mip->tech == MI_TECH_sts) {
          sip = GetSeqIdForGI (uid);
          if (sip != NULL) {
            SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_ONLY, sizeof (buf) - 1);
            SeqIdFree (sip);
            if (StringDoesHaveText (buf)) {
              dbt->db = MemFree (dbt->db);
              dbt->db = StringSave ("dbSTS");
              oip->id = 0;
              oip->str = StringSave (buf);
            }
          }
        } else {
          Message (MSG_POST, "Tech for gi %ld is %d", (long) uid, (int) mip->tech);
        }
      }
    }
    BioseqUnlock (bsp);
  }
}

static void DoGI2dbSTSConversion (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, GI2dbSTSCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void SplitSemicolonedCallback (BioSourcePtr biop, Pointer userdata)

{
  CharPtr     ptr, str;
  OrgModPtr   nwomp, omp;
  OrgNamePtr  onp;
  OrgRefPtr   orp;

  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_specimen_voucher ||
        omp->subtype == ORGMOD_culture_collection ||
        omp->subtype == ORGMOD_bio_material) {
      str = omp->subname;
      ptr = StringChr (str, ';');
      if (ptr == NULL) continue;
      *ptr = '\0';
      ptr++;
      nwomp = OrgModNew ();
      nwomp->subtype = omp->subtype;
      nwomp->subname = StringSave (ptr);
      nwomp->next = omp->next;
      omp->next = nwomp;
    }
  }
}

static void SplitSemicolonedVouchers (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioSourcesInSep (sep, NULL, SplitSemicolonedCallback);

  BasicSeqEntryCleanup (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void PrintDiscrepancyTestListMenuItem (IteM i)

{
  Char         path [PATH_MAX];
  FILE *fp;

  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
  PrintDiscrepancyTestList (fp);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Discrepancy Tests List");
  FileRemove (path);
}


static void TestNewParse (IteM i)
{
  BaseFormPtr  bfp;
  ForM         f;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  f = SingleParseAction (bfp->input_entityID);
  if (f != NULL) {
    Show (f);
  }
}


static void MakeSpecialOrganismMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;
  MenU  x;
  
  s = SubMenu (m, "Organism/ G");
  i = CommandItem (s, "Parse Text After Bionomial", ParseTaxnameToSourceQual);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Parse ATCC Strain to Culture Collection", AtccStrainToCultureCollection);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Trim Organism Name", TrimOrganismName);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Genus-Species Fixup", GenSpecTaxonFixup);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Prefix Authority with Organism", PrefixAuthorityWithOrganism);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Append Modifier to Organism", AddModToOrg);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Parse Trinomial into Subspecies", ParseTrinomial);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Country Fixup");
  i = CommandItem (x, "Do Not Fix Capitalization After Colon", CountryLookupWithoutCapFix);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Fix Capitalization After Colon", CountryLookupWithCapFix);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Lat-lon Tool", LatLonTool);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Lat-Lon Country Conflict Tool", LatLonCountryTool);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Fix Specific-Host Values", FixSpecificHostValues);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Make Bad Specific-Host Table", MakeBadSpecificHostValueTable);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "List Failed Taxonomy Lookups", ListFailedTaxonomyLookups);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "TaxFix Tool", TaxFixTool);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Abbreviate Voucher Institution Codes", UpdateVoucherName);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  x = SubMenu (s, "Parse CollectionDate formats");
  i = CommandItem (x, "Month First", ParseCollectionDateMonthFirst);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Day First", ParseCollectionDateDayFirst);
  SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
  x = SubMenu (s, "Source Focus");
  i = CommandItem (x, "Set", SetSourceFocus);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Clear", ClearSourceFocus);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Set on Source Desc when Source Feat Present");
  i = CommandItem (x, "Transgenic", SetTransgenicOnSourceDescWhenSourceFeatPresent);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Focus", SetFocusOnSourceDescWhenSourceFeatPresent);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  i = CommandItem (s, "Consolidate Organism Notes", ConsolidateOrganismNotes);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Consolidate Like Modifiers");
  i = CommandItem (x, "With semicolons", ConsolidateLikeModifiersWithSemicolons);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Without semicolons", ConsolidateLikeModifiersWithoutSemicolons);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu(s, "Parse Organism Modifiers");
  i = CommandItem (x, "Load From File", LoadOrganismModifierTable);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Load Tax Consult File", LoadTaxConsult);
  SetObjectExtra (i, bfp, NULL);
/*  i = CommandItem (x, "New Source Table Reader", NewLoadSourceQualifierTable);
  SetObjectExtra (i, bfp, NULL); */
  i = CommandItem (x, "Export To File", ExportOrganismTable);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Export Last Lineage Table", ExportLastLineage);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Load Tax Table", LoadTaxTableReader);
  SetObjectExtra (i, bfp, NULL);  
  i = CommandItem (s, "Parse File to Source", ParseFileToSource);
  SetObjectExtra (i, bfp, NULL);  
  SeparatorItem (s);
  x = SubMenu (s, "Add Type Strain Comments");
  i = CommandItem (x, "To All", AddTypeStrainCommentsToAll);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Constraints", AddTypeStrainCommentsWithConstraint);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Influenza Virus Names");
  i = CommandItem (x, "Parse Strain,Serotype from Names",
                   ParseInfluenzaAVirusNames);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Add Strain,Serotype to Names",
                   AddStrainAndSerotypeToInfluenzaAVirusNames);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Fixup Organism Names", FixupInfluenzaAVirusNames);
  SetObjectExtra (i, bfp, NULL);

  x = SubMenu (s, "Convert Primer Quals");
  i = CommandItem (x, "Seq to Name", ConvertPrimerSeqToName);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Name to Seq", ConvertPrimerNameToSeq);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Trim Junk in Primer Seqs", TrimPrimerSeqJunk);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  i = CommandItem (s, "Merge Multiple BioSources On Sequence", MergeBiosources);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  i = CommandItem (s, "Split Structured Collections at Semicolon", SplitSemicolonedVouchers);
  SetObjectExtra (i, bfp, NULL);
}

static void MakeSpecialProjectsMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;
  MenU  x;

  s = SubMenu (m, "Projects/ J");
  if (indexerVersion) {
    i = CommandItem (s, "Create Seq-Hist for TPA", CreateSeqHistTPA);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Create Seq-Hist for TPA - Detailed", 
                     CreateSeqHistTPADetailed);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Load TPA Accession Numbers From File",
                     LoadTPAAccessionNumbersFromFile);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Import Single Interval Alignment for TPA Seq-hist", ImportAlignmentForSeqHistInterval);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Advanced Assembly Alignment Editor", AdvancedAssemblyAlignmentEditor);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Assembly Alignment Interval Resolution", AssemblyAlignmentIntervalResolution);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);    
    i = CommandItem (s, "Create TPA BLAST Database", CreateTPABlastDB);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Leave First TPA Assembly", LeaveFirstSeqHistAssembly);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Read Seq-Hist Assembly BLASTs", ReadSeqHistAssembly);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Remove Seq-Hist Assembly", RemoveSeqHistAssembly);
    SetObjectExtra (i, bfp, NULL);    
    i = CommandItem (s, "Infer Seq-Hist from Delta", CreateSeqHistDelta);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Add Genomes User Object", AddGenomesUserObject);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Add Genomic Phrap Graphs", MakePhrap);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    x = SubMenu (s, "Add Gaps to HTGS Deltas");
    i = CommandItem (x, "Add Gaps of 100", AddDelta100Gaps);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "Add Gaps of 0", AddDelta0Gaps);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Remove Gaps at ends of HTGS Deltas", DeleteDeltaEndGaps);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Parse Structured Comment into User Object", ParseCommentIntoStructuredObject);
    SetObjectExtra (i, bfp, NULL);
    x = SubMenu (s, "Add Structured Comment Prefix-Suffix");
    i = CommandItem (x, "HIVData", AddStructuredCommentHivTag);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "FluData", AddStructuredCommentFluTag);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "BARCODE Discrepancy Tool", BarcodeTestTool);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Apply BARCODE Dbxrefs", ApplyBarcodeDbxrefs);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Add Note to Coding Regions in Source Features", ApplyTagToCodingRegionsInSourceFeatures);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Add TSA Table", AddTSATableToBioseq);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Refresh TSA Tables", RefreshTSATables);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Create TSA IDs", CreateTSAIDsFromLocalIDs); 
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Edit TSA Assembly", EditTSAAssembly);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Add mRNA sequence to deflines", AddmRNASequenceToDeflines);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Create RefSeq Protein IDs", CreateRefSeqProteinIDs);
    SetObjectExtra (i, bfp, NULL); 
    SeparatorItem (s);
    i = CommandItem (s, "Add Flu Comments", AddFluComments);
    SetObjectExtra (i, bfp, NULL); 
    SeparatorItem (s);
    i = CommandItem (s, "Create Structured Comments from Table", CreateStructuredCommentsItem);
    SetObjectExtra (i, bfp, NULL); 

    /* remove after retro */
    SeparatorItem (s);
    i = CommandItem (s, "GI to dbSTS Conversion", DoGI2dbSTSConversion);
    SetObjectExtra (i, bfp, NULL);
  }
}

static void MakeSpecialLocusTagMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;

  s = SubMenu (m, "Locus_Tags");
  i = CommandItem (s, "Insert Gene Locus Tag Prefix", InsertGeneLocusTagPrefix);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Collect Locus for Locus Tag Gene Xrefs", FixLocusTagGeneXrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Convert Gene Locus Tag to Old Locus Tag", ConvertLocusTagToOldLocusTag);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Copy Gene Locus to Empty Gene Locus Tag", CopyLocusToLocusTag);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Replace Repeat Region Locus Tags with FLYBASE Dbxref",
                   ReplaceRepeatRegionLocusTagWithDbxref);
  SetObjectExtra (i, bfp, NULL);  
  
}


static void MakeSpecialApplyMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s, x;

  s = SubMenu (m, "Apply/ A");
  i = CommandItem (s, "Add CDS", ApplyCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add RNA", ApplyRRNA);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Add Named rRNA");
  i = CommandItem (x, "12S", Apply12SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "16S", Apply16SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "23S", Apply23SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "18S", Apply18SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "28S", Apply28SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "26S", Apply26SRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Small", ApplySmallRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Large", ApplyLargeRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Apply rRNA_ITS", ApplyRNA_ITS);
  SetObjectExtra (i, bfp, NULL); 

  i = CommandItem (s, "Add other Feature", ApplyImpFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add Qualifier", MacroApplyGBQual);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Add Source Qual", MacroApplySourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add CDS-Gene-Prot-mRNA Qual", MacroApplyCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add RNA Qual", MacroApplyRNAQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Apply Structured Comment Field", MacroApplyStructuredComment);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Add Molecule Information", ApplyMolInfo);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Add Gene Pseudo");
  i = CommandItem (x, "To All", AddGenePseudo);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Where Feature is Pseudo", MarkGenesWithPseudoFeaturesPseudo);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add Global Code Break", AddGlobalCodeBreak);
  SetObjectExtra (i, bfp, NULL);  
  if (indexerVersion)
  {
    x = SubMenu (s, "Add Keyword");
    i = CommandItem (x, "GDS", ApplyGDSKeyword);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "TPA:inferential", ApplyTPAInferentialKeyword);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "TPA:experimental", ApplyTPAExperimentalKeyword);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "TPA:reassembly", ApplyTPAReassemblyKeyword);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "With Constraint", ApplyKeywordWithStringConstraint);
    SetObjectExtra (i, bfp, NULL);
  }
  else
  {
    i = CommandItem (s, "Add GDS Keyword", ApplyGDSKeyword);
    SetObjectExtra (i, bfp, NULL);
  }
  SeparatorItem(s);  
  if (indexerVersion) {
    SetupEditSecondary (s, bfp);
  }
  x = SubMenu (s, "Load Secondary Accessions from File");
  i = CommandItem (x, "Extra Accessions",
                   LoadSecondaryAccessionNumbersFromFile);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "History Takeover",
                   LoadHistoryAccessionNumbersFromFile);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Load Genome Project IDs from File", LoadGenomeProjectIDsFromFile);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "Load Feature Qualifiers from File",
                   LoadFeatureQualifierTable);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Load Feature Fields from File", LoadFeatureFieldTable);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Table Reader", NewLoadFeatureQualifierTable);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Export Table", ExportQualifiers);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Export Bankit Comments", ExportBankitComments);
  SetObjectExtra (i, bfp, NULL);
}


static void MakeSpecialRemoveMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s, x, y;

  s = SubMenu (m, "Remove/ R");
  i = CommandItem (s, "Remove Descriptors", RemoveDescriptor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Features", FeatureRemove);
  SetObjectExtra (i, bfp, NULL);

  x = SubMenu (s, "Remove Features by Selection");
  i = CommandItem (x, "Selected", RemoveSelFeats);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Unselected", RemoveUnselFeats);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "All Features Outside Range of Selected Feature", RemoveFeatsOutOfSelRange);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Remove Unindexed Features", RemoveUnindexedFeatures);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Remove Genes by Features", RemoveGeneByUnderlyingFeatureType);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Remove RBSs Genes", RemoveRBSGenes);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "Remove Qualifiers", MacroRemoveGBQual);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Remove Source Qualifier", MacroRemoveSourceQual);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Remove CDS-Gene-Prot-mRNA Qual", MacroRemoveCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Remove RNA Qual", MacroRemoveRNAQual);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Remove Structured Comment Field", MacroRemoveStructuredComment);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);

  x = SubMenu (s, "Remove Proteins");
  i = CommandItem (x, "Just Remove Proteins", RemoveProteins);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "And Renormalize Nuc-Prot Sets", RemoveProteinsAndRenormalize);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Orphaned Proteins", RemoveOrphanedProteins);
  SetObjectExtra (i, bfp, NULL);

  x = SubMenu (s, "Remove CDD");
  i = CommandItem (x, "Regions", RemoveCDDRegions);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Alignments", RemoveCDDAligns);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Duplicates", RemoveCDDDups);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Overlapping (Retain Best)", RetainBestCDD);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Nth Seg's CDS for Each SegSet",
		   RemoveNthCDSFromSegSets_Callback);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);

  i = CommandItem (s, "Remove BankIt Comments", RemoveBankitComments);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);

  x = SubMenu (s, "Remove Text");
  i = CommandItem (x, "Inside String", RemoveTextInsideString);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Outside String", RemoveTextOutsideString);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  
  x = SubMenu (s, "Remove Duplicate");
  i = CommandItem (x, "All Feat Types", NewRemoveDuplicateFeatsAll);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Feats With Options", RemoveDuplicateFeatsWithOptions);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Gene Ontology Terms", RemoveDuplicateGeneOntology);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "PCR Primers", RemoveDuplicatePCRPrimers);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Remove Redundant Proprotein Misc Feats",
                   RemoveRedundantProproteinMiscFeats);
  SetObjectExtra (i, bfp, NULL);
  
  SeparatorItem (s);

  x = SubMenu (s, "Remove Sequences");
  i = CommandItem (x, "From Alignments", RemoveSequencesFromAlignment);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "From Record", RemoveSequencesFromRecord);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Remove Alignments", RemoveAlignment);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Graphs", RemoveGraph);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  
  x = SubMenu (s, "Remove Xrefs");
  y = SubMenu (x, "Gene Xrefs");
  i = CommandItem (y, "Choose Types and Constraint", NewRemoveGeneXrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Unnecessary Gene Xrefs", RemoveGeneXrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Non-Suppressing Gene Xrefs", RemoveNSGeneXrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Gene Xrefs with Orphan Gene Locus", RemoveGeneXrefsNoGene);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Gene Xrefs with Orphan Locus_tag", RemoveGeneXrefsNoLocTag);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "ALL Gene Xrefs", RemoveAllGeneXrefs);
  SetObjectExtra (i, bfp, NULL);
  
  y = SubMenu (x, "Taxon Xrefs");
  i = CommandItem (y, "From Features", RemoveTaxonXrefsFromFeatures);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From Features and BioSources", 
                   RemoveTaxonXrefsFromFeaturesAndBioSources);
  SetObjectExtra (i, bfp, NULL);
  
  y = SubMenu (x, "Db_xrefs");
  i = CommandItem (y, "From Genes", RemoveGeneDbxrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From RNAs", RemoveRNADbxrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From CDSs", RemoveCDSDbxrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From All Features", RemoveAllFeatureDbxrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From All BioSources", RemoveAllBioSourceDbxrefs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "From All Features and BioSources", RemoveAllDbxrefs);
  SetObjectExtra (i, bfp, NULL);
  
  SeparatorItem (s);
  y = SubMenu (s, "Remove Titles");
  x = SubMenu (y, "Protein");
  i = CommandItem (x, "All but RefSeq", DoRemoveProtTitles);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "All", DoRemoveAllProtTitles);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Hypothetical", RemoveHypotheticalProteinTitles);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Inconsistent", RemoveInconsistentProteinTitles);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (y, "Nucleotide");
  i = CommandItem (x, "All", DoRemoveAllNucTitles);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "mRNA", DoRemoveMrnaTitles);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "NPS", DoRemoveNPSTitles);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  if (indexerVersion) {
    y = SubMenu (s, "Remove IDs, SeqIDs, GI_IDs, LOCUS");
    x = SubMenu (y, "Local");
    i = CommandItem (x, "From Nucleotides", RemoveNucLocalIDs);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "From Proteins", RemovePrtLocalIDs);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "From All Bioseqs", RemoveLocalIDs);
    SetObjectExtra (i, bfp, NULL);

    x = SubMenu (y, "GenBank");    
    i = CommandItem (x, "From Proteins", RemoveGBIDsFromProteins);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "From All Bioseqs", RemoveGBIDsFromBioseqs);
    SetObjectExtra (i, bfp, NULL);
    
    i = CommandItem (y, "Seq-ID Name From Prot Feats", RemoveSeqIdName_Callback);
    SetObjectExtra (i, bfp, NULL);
    
    i = CommandItem (y, "GI_IDs from Bioseqs", RemoveGIsFromBioseqs);
    SetObjectExtra (i, bfp, NULL);

    i = CommandItem (y, "LOCUS from Parts", RemoveLocusFromParts);
    SetObjectExtra (i, bfp, NULL);    
    SeparatorItem (s);
  }
  
  i = CommandItem (s, "Remove Empty Genome Project ID Descriptors", RemoveEmptyGenomeProjectIDs);
  SetObjectExtra (i, bfp, NULL);
}

static void MakeSpecialConvertMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s, x, y, z;

  s = SubMenu (m, "Convert/ C");
  i = CommandItem (s, "Convert Features", ConvertFeatures);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Convert Qualifiers", MacroConvertGBQual);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  i = CommandItem (s, "Convert Source Qualifier", MacroConvertSourceQual);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Convert CDS-Gene-Prot-mRNA Qual", MacroConvertCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Convert RNA Qual", MacroConvertRNAQual);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  i = CommandItem (s, "Convert CDS to Misc_feat or Pseudogene", ConvertCDSMenuItem);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Convert CDS with Internal Stop Codon");
  i = CommandItem (x, "To Misc Feat", NewConvertCDSWithInternalStopsToMiscFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "To Misc Feat (Viral)", NewConvertCDSWithInternalStopsToMiscFeatViral);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Convert Bad Inference Qualifiers to Notes", ConvertBadInferenceToNote);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Convert transl-except in CDSs to RNA Editing Exception", ConvertTranslExceptToRNAEditingException);
  SetObjectExtra (i, bfp, NULL);  
  
  SeparatorItem (s);

  i = CommandItem (s, "Swap Source Qualifiers", MacroSwapSourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Swap CDS-Gene-Prot-mRNA Qualifiers", MacroSwapCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Swap RNA Qualifiers", MacroSwapRNAQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Swap Qualifiers", MacroSwapGBQual);
  SetObjectExtra (i, bfp, NULL);



  SeparatorItem (s);

  i = CommandItem (s, "Parse Text", ParseFlatfileToSourceQual);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  x = SubMenu (s, "Parse from tRNA Comment");
  i = CommandItem (x, "Codons", ParseCodonsFromtRNAComment);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Anticodons", ParseAntiCodonsFromtRNAComment);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Parse Codon Quals to Code Breaks", ParseCodonQualToCodeBreak);
  SetObjectExtra (i, bfp, NULL);
  
  SeparatorItem (s);

  i = CommandItem (s, "Extra Accession to History by Position", ExtraAccToHistByPos);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  x = SubMenu (s, "Convert Descriptors and Features");

  y = SubMenu (x, "Descriptor to Feature");
  i = CommandItem (y, "Source", ConvertSrcDescToFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Comment", ConvertComDescToFeat);
  SetObjectExtra (i, bfp, NULL);
  z = SubMenu (y, "Publication");
  i = CommandItem (z, "All", ConvertPubDescToFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (z, "By Constraint", CreateConvertPubDescStringConstraintDialogX);
  SetObjectExtra (i, bfp, NULL);

  y = SubMenu (x, "Convert Full Length Feature to Descriptor");
  i = CommandItem (y, "Source", ConvertSrcFeatToDesc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Comment", ConvertCommentFeatToDesc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (y, "Publication", ConvertPubFeatToDesc);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (x, "BankIt Comment to Comment Descriptor", CopyBankitComments);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  x = SubMenu (s, "Convert Second Protein Name and Description");
  i = CommandItem (x, "Second Protein Name to Description", SecondProtNameToDesc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Protein Description to Second Name", ProtDescToSecondName);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Whole Location to Interval", ConvertWholeToInterval);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Both Strands to Proper Strand", ConvertBothStrands);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Refresh Gene Xrefs", RefreshGeneXRefs);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Convert Biosource Dbxrefs to Feature Dbxrefs", ConvertBioSourceDbxrefToFeatureDbxref);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
#if 0
  i = CommandItem (s, "Convert SegSets to Delta Sequences", ConvertSegSetToDeltaItem);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
#endif
  if (indexerVersion) {
    x = SubMenu (s, "Convert IDs, Accessions, or GenBank.names");
    i = CommandItem (x, "LocalIDs to GeneralIDs", ConvertLocalToGeneral);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "Convert GeneralIDs to LocalIDs", ConvertGeneralIdToLocalID);
    SetObjectExtra (i, bfp, NULL);
    y = SubMenu (x, "Accession to LocalIDs");
    i = CommandItem (y, "For All Sequences", ConvertToLocalProcAll);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (y, "For Nucleotide Sequences", ConvertToLocalProcOnlyNucs);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (y, "For Protein Sequences", ConvertToLocalProcOnlyProts);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "GenBank.name to LocalID", ChangeGenBankNameToLocal);
    SetObjectExtra (i, bfp, NULL);
  } 
}


static void MakeSpecialEditMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s, x;

  s = SubMenu (m, "Edit/ E");
  i = CommandItem (s, "Edit Qualifiers", MacroEditGBQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit Source Qual", MacroEditSourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit CDS-Gene-Prot-mRNA Qual", MacroEditCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit RNA Qual", MacroEditRNAQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit Structured Comment Field", MacroEditStructuredComment);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit Coding Region Frame", ApplyCDSFrame);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Edit Molecule Information", EditMolInfoFields);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "New Edit Molecule Information", NewMolInfoFieldsEdit);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Edit Feature");
  i = CommandItem (x, "Evidence", FeatureEvidenceEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Exception", FeatureExceptionEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Location", FeaturePartialEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Strand", FeatureStrandEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Citation", FeatureCitationEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Experiment", FeatureExperimentEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Inference", FeatureInferenceEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Pseudo", FeaturePseudoEditor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Reverse Location Interval Order", ReverseFeatureIntervals);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Resolve Intersecting Feature Locations", ResolveFeatureOverlaps);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Update Features", UpdateFeatures);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  x = SubMenu (s, "Bulk Edit");
  i = CommandItem (x, "CDS", BulkEditCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Gene", BulkEditGene);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "RNA", BulkEditRNA);
  SetObjectExtra (i, bfp, NULL);


  SeparatorItem (s);
  i = CommandItem (s, "Lowercase Qualifiers", FixCaseByField);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  x = SubMenu (s, "Extend Partial Features");
  i = CommandItem (x, "All", ExtendPartialFeatures);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Constraint", ExtendPartialFeaturesWithConstraint);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Extend Gene", ExtendGeneReg);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Correct Genes for CDSs or mRNAs", CorrectGenes);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Correct Strand for primer_bind Pairs", SetPrimerBindPairStrands);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion) {
    SeparatorItem (s);
    i = CommandItem (s, "Edit Locus", EditLocusProc);
    SetObjectExtra (i, bfp, NULL);
    
    SeparatorItem (s);
    i = CommandItem (s, "Edit Sequence Ends", EditSeqEndsProc);
    SetObjectExtra (i, bfp, NULL);
  }

  i = CommandItem (s, "Trim Ns from Bioseqs", TrimNsFromNucs);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion) {
    x = SubMenu (s, "Put Ns under VecScreen Hits");
    i = CommandItem (x, "Strong", VecScreenToNsStrong);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "Strong + Moderate", VecScreenToNsModerate);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "All", VecScreenToNsWeak);
    SetObjectExtra (i, bfp, NULL);
  }
}

static void MakeSpecialClearMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s, x;

  s = SubMenu (m, "Clear/ L");
  i = CommandItem (s, "Clear CDS Products", ClearCdsProducts);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Clear mRNA Products", ClearMrnaProducts);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);

  i = CommandItem (s, "Clear ORF Flag in CDSs", ClearOrfFlagInCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Clear transl-except in CDSs", ClearTranslExceptInCDS);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Clear Gene Pseudo");
  i = CommandItem (x, "All", ClearGenePseudo);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Constraint", ClearGenePseudoConstraint);
  SetObjectExtra (i, bfp, NULL);

  SeparatorItem (s);
  
  i = CommandItem (s, "Clear GenBank Components", (ItmActnProc) EditGenbankElements);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion)
  {
    i = CommandItem (s, "Clear Keywords", RemoveKeywordWithStringConstraint);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Clear Nomenclature", RemoveNomenclature);
    SetObjectExtra (i, bfp, NULL);    
  }

  SeparatorItem (s);

  i = CommandItem (s, "Clear Serial Numbers", StripSerials);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Clear Figure from Pubs", RemovePubFig);
  SetObjectExtra (i, bfp, NULL);
}


static void MakeSpecialGapSequencesMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;
  MenU  x;

  /* Gap Sequences Menu */
  s = SubMenu (m, "Gap Sequences/ G");
  x = SubMenu (s, "Raw Sequence to Delta");
  i = CommandItem (x, "By Ns", 
                   RawSeqToDeltaSeqUnknownWithUnknownLengthGaps);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "By Location",
                   ConvertRawToDeltaWithGapLocations);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Delta Sequence to Raw", SegSeqToRawSeq);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Far Delta to Near Delta", FarDeltaToNear);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Convert Selected Gaps");
  i = CommandItem (x, "To Known Length", ConvertSelectedGapFeaturesToKnown);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "To Unknown Length", ConvertSelectedGapFeaturesToUnknown);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Convert Known to Unknown (by Size)", ConvertGapFeaturesToUnknown);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit Length of Selected Known Gaps", ChangeKnownGapLength);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Expand Known Gaps to Include Flanking Ns", AddFlankingNsToKnownLengthGaps);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Combine Adjacent Gaps", CombineAdjacentGaps);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Adjust Features for Gaps", AdjustFeaturesForGaps);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Trim CDS-mRNA for Known Gaps", AdjustCodingRegionsEndingInGap);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Convert CDS with Internal Gaps to Misc_feat", ConvertCodingRegionsWithInternalKnownGapToMiscFeat);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Fix Delta SeqDataLen", FixDeltaSeqDataLen);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "List Delta Sequences", DeltaReport);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Instantiate Gap Features from Delta Instructions", InstantiateGapFeatures);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Unnecessary Gap Features", RemoveUnnecessaryGapFeatures);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Flip Intervals Between Gaps", FlipSequenceIntervals);
  SetObjectExtra (i, bfp, NULL);
}

static void MakeSpecialSelectMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;

  s = SubMenu (m, "Select/ S");
  i = CommandItem (s, "Select Descriptors", SelectDescriptor);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Select Features", SelectFeatures);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Select Sequences", SelectBioseq);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Select Publications", SelectPubs);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion) {
    i = CommandItem (s, "Find Contig for Features", FindContig);
    SetObjectExtra (i, bfp, NULL);
  }  
}


static void MakeSpecialLinkMenu (MenU m, BaseFormPtr bfp)
{
  IteM  i;
  MenU  s;
  MenU  x;

  s = SubMenu (m, "Link/ K");
  x = SubMenu (s, "Feature IDs");
  i = CommandItem (x, "Assign", AssignFeatIDs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Reassign", ReassignFeatIDs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Unique in Merged Set", UniqifyFeatIDs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Clear IDs and Links", ClearFeatIDsAndLinks);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Create CDS-mRNA Links");
  i = CommandItem (x, "By Overlap", LinkByOverlap);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "By Product", LinkByProduct);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "By Label", LinkByLabel);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Tool", CDSmRNALinkTool);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Link Selected CDS and mRNA Pair", LinkSelected);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Show Linked CDS or mRNA Feature", SelCDSmRNALink);
  SetObjectExtra (i, bfp, NULL);
}


static void MakeSpecialDesktopMenu (MenU m, BaseFormPtr bfp)
{
  MenU s;
  IteM i;

  s = SubMenu (m, "Desktop");
  i = CommandItem (s, "Reorder by ID", ReorderSetByAccessionMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Repackage Parts", RepackagePartsMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Descriptor Propagate", DescriptorPropagateMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Segregate Sets", NewSegregateBioseqSetMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Sequester Sets", SequesterSequencesMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Get Rid of Seg Gaps", GetRidOfSegGapMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Create SeqAlign", GenerateSeqAlignMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Update SeqAlign", UpdateSeqAlignMenuItem);
  SetObjectExtra (i, bfp, NULL);
#if defined(OS_UNIX) || defined(OS_MSWIN) 
  i = CommandItem (s, "Correct RNA strandedness (Use SMART)", CorrectRNAStrandednessMenuItem);
  SetObjectExtra (i, bfp, NULL);
#endif 
  i = CommandItem (s, "Bioseq RevComp By ID", BioseqRevCompByIDMenuItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Sets in Set", RemoveSetsInSetMenuItem);
  SetObjectExtra (i, bfp, NULL);


}


extern void SetupSpecialMenu (MenU m, BaseFormPtr bfp)

{
  IteM  i;
  MenU  s;
  MenU  x, y;

  s = SubMenu (m, "Def Line/ D");
  x = SubMenu (s, "Automatic Def Line");
  i = CommandItem (x, "Default Options", AutoDef);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "No Modifiers", AutoDefWithoutModifiers);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Select Options...", AutoDefWithOptions);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Add Definition Line", ApplyTitle);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Prefix Def Line With");
  i = CommandItem (x, "Organism", AddOrgNameToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Strain", AddStrainToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Clone", AddCloneToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Isolate", AddIsolateToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Haplotype", AddHaplotypeToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Cultivar", AddCultivarToDefLines);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "List", PrefixDefLines);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Parse Text", ParseDefLineToSourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Parse Modifiers from Defline", ParseModifiersFromDefline);
  SetObjectExtra (i, bfp, NULL);

  MakeSpecialOrganismMenu (m, bfp);

  MakeSpecialApplyMenu (m, bfp);

  MakeSpecialRemoveMenu (m, bfp);

  MakeSpecialConvertMenu (m, bfp);

  MakeSpecialEditMenu (m, bfp);

  MakeSpecialClearMenu (m, bfp);

  s = SubMenu (m, "CDS, Protein, and RNA");
  i = CommandItem (s, "Fix Product Name Capitalization", FixProductCapitalization);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Retranslate Coding Regions");
  i = CommandItem (x, "Obey Stop Codon", RetranslateCdRegionsNoStop);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Ignore Stop Codon", RetranslateCdRegionsDoStop);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Ignore Stop Codon Except at End of Complete CDS", RetranslateCdRegionsNoStopExceptEndCompleteCDS);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Recompute Suggested Intervals");
  i = CommandItem (x, "And Update Genes", RecomputeSuggestFixGenes);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Do Not Update Genes", RecomputeSuggest);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Retranscribe mRNA Products", ReprocessmRNAProducts);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Regenerate Peptide Products", ReprocessPeptideProducts);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Resynchronize Partials");
  i = CommandItem (x, "CDS", ResynchronizeCDSPartials);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "mRNA", ResynchronizeMRNAPartials);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Peptide", ResynchronizePeptidePartials);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Cleanup CDS partials after propagation", FixCdsAfterPropagate);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Truncate Proteins and CDSs at Stops", TruncateProteins);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Extend Proteins to Stops", ExtendProteins);
  SetObjectExtra (i, bfp, NULL);
  
  i = CommandItem (s, "Trim Protein Feature Lengths", TrimProtFeatLengths);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Update Proteins From CDS", UpdateProteinsFromCDS);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Adjust Features for Introns");
  i = CommandItem (x, "CDS", RemoveIntronLocationsFromCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "rRNA", RemoveIntronLocationsFromrRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "tRNA", RemoveIntronLocationsFromtRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "mRNA", RemoveIntronLocationsFrommRNA);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Adjust Last Exons for Stop Codons");
  i = CommandItem (x, "Do Not Change Last Exon Partial", FixLastExonLocNoPartial);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Make Last Exon 3' Partial", FixLastExonLocMakePartial);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Convert CDS to mat_peptide", CombineMultipleCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Extract Protein Features From CDS comments", ExtractProteinFeaturesFromNote);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Add Transl Excepts");
  i = CommandItem (x, "Without Comment", GlobalAddTranslExcept);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Comment", AddTranslExceptWithComment);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Instantiate Mat-peptide Products", InstantiateMatPeptideProductsMenuItem);
  SetObjectExtra (i, bfp, NULL);
    
  /* Authors and Publications */
  s = SubMenu (m, "Authors and Publications");
  i = CommandItem (s, "Edit Publications", EditPubs);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Fix Capitalization");
  i = CommandItem (x, "All", FixAllPubCapitalization);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Authors", FixAuthorCapitalization);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Titles", FixTitleCapitalization);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Affiliations", FixAffiliationCapitalization);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  x = SubMenu (s, "Remove Publications");
  i = CommandItem (x, "Unpublished", RemoveUnpublishedPublications);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Published", RemovePublishedPublications);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Colliding", RemoveCollidingPublications);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "All", RemovePub);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Remove Author Consortiums", RemovePubConsortiums);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Reverse Author Order", FixNameOrder);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Strip Author Suffixes", StripAuthorSuffix);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Truncate Author Middle Initials", TruncateAuthorMiddle);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Convert Author to Consortium");
  i = CommandItem (x, "All", ChangeAllAuthorNameToConsortium);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With String Constraint", ChangeAuthorNameToConsortiumWithConstraint);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Where Last Name Contains 'Consortium'", ChangeAuthorNameWithConsortiumToConsortium);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Reload Publications", LookupAllPubs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add Cit-sub for update", AddCitSubForUpdate);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Fix USA and States", AbbreviateCitSubAffilStates);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Load Individual Author Affils from File", SplitPubsByTable);
  SetObjectExtra (i, bfp, NULL); 
  
  MakeSpecialGapSequencesMenu (m, bfp);

  s = SubMenu (m, "Make Features and Fuse");
  i = CommandItem (s, "Gene Features from Xrefs", XrefToGene);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Gene Xrefs from Features", GeneToXref);
  SetObjectExtra (i, bfp, NULL);
  i  = CommandItem (s, "Gene Features From Other Features", FeatureToGene);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "CDS Features from Gene, mRNA, or exon", FeatureToCds);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Make AltSplice CDSs for mRNAs", MakeCDSmRNAPairs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "mRNA Features from Gene, CDS, or exon", FeatureTomRNA);
  SetObjectExtra (i, bfp, NULL);  
  x = SubMenu (s, "Exon features from");
  i = CommandItem (x, "CDS Intervals", MakeExonsFromCDSIntervals);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "mRNA Intervals", MakeExonsFromMRNAIntervals);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  i = CommandItem (s, "Fuse Features", FuseFeature);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion) {
    i = CommandItem (s, "Fuse Joins in Locations", FuseSlpJoins);
    SetObjectExtra (i, bfp, NULL);
  }
  SeparatorItem (s);
  i = CommandItem (s, "Combine Selected Genes into Pseudogene", CombineToCreatePseudoGene);
  SetObjectExtra (i, bfp, NULL);
  
  s = SubMenu (m, "Normalize, Map, Package");
  i = CommandItem (s, "Renormalize Nuc-Prot Sets", RenormalizeNucProtSetsMenuItem);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "Package Features on Parts", PackageOnParts);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Package Proteins on Nucleotides", PackageOnNucs);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  x = SubMenu (s, "Map Feature Intervals to Parts");
  i = CommandItem (x, "Joined Location", MergeToPartsJoin);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Ordered Location", MergeToPartsOrdered);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Map Intervals to Segmented Sequence", MergeToSegSeq);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Map Features to Protein Sequences", MapFeaturesToProteinSequence);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  if (indexerVersion) {
    i = CommandItem (s, "Connect CDS to Closest Protein", ReconnectCDSProduct);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
  }
  
  i = CommandItem (s, "Remap Feature Locations Based on Alignment Coordinates", FixFeatureIntervals);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);

  i = CommandItem (s, "Map RBSs to Downstream Gene", MapRBSsToDownstreamGene);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "Suppress Genes On Features", SuppressGenesOnFeatures);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "CDSGeneRange Error");
  i = CommandItem (x, "Suppress", SuppressCDSGeneRangeError);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Restore", RestoreCDSGeneRangeError);
  SetObjectExtra (i, bfp, NULL);

  i = CommandItem (s, "Suppress genes inside mobile elements", SuppressGenesOnFeaturesInsideMobileElements);
  SetObjectExtra (i, bfp, NULL);


  MakeSpecialProjectsMenu (m, bfp);

  s = SubMenu (m, "GenProdSets and SeqAnnots");
  if (indexerVersion) {
    x = SubMenu (s, "Promote NucProtSet to GenProdSet");
    y = SubMenu (x, "By Overlap");
    i = CommandItem (y, "Locus Tag ID", NPStoGPSLocusTag);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (y, "Arbitrary ID", NPStoGPSArbitrary);
    SetObjectExtra (i, bfp, NULL);
    y = SubMenu (x, "By Feature ID");
    i = CommandItem (y, "Locus Tag ID", NPStoGPSLocusTagFeatID);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (y, "Arbitrary ID", NPStoGPSArbitraryFeatID);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Demote GenProdSet To NucProtSet", GPStoNPS);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Move CDS To NucProtSet", CDStoNPS);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Add mRNA Titles to GenProdSet", MrnaTitlesToGPS);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
  }
  
  x = SubMenu (s, "Generate GenProdSet Redundancy");
  i = CommandItem (x, "With Descriptor Propagation", MakeRedundantGPSwithProp);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Without Descriptor Propagation", MakeRedundantGPSnoProp);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "Just Copy CDS Prot Xref", MakeRedundantGPSjustXref);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Copy Far Prot Xref to CDS", MakeRedundantGenomicProtXref);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  
  i = CommandItem (s, "Force Features into one SeqAnnot", ForceRepackage);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
 
  if (indexerVersion) {
    x = SubMenu (s, "Add BLAST tag to Seq-annot");
    i = CommandItem (x, "BLASTN", AddBlastNToSeqAnnot);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "BLASTP", AddBlastPToSeqAnnot);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "BLASTX", AddBlastXToSeqAnnot);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "TBLASTN", AddTBlastNToSeqAnnot);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);

    x = SubMenu (s, "Remove Seq-annot");
    i = CommandItem (x, "ids", RemoveSeqAnnotIDs);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "locs", RemoveSeqAnnotLOCs);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);

    x = SubMenu (s, "Promote Features to ID by Quality");
    i = CommandItem (x, "Best ID", PromoteToBestIDProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (x, "Worst ID", PromoteToWorstIDProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Promote Alignments to Best ID", PromoteAlignsToBestIDProc);
    SetObjectExtra (i, bfp, NULL);
  }
  
  s = SubMenu (m, "Cleanup and Fixup");
  i = CommandItem (s, "Clean Up Record", ForceCleanup);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Force Taxonomy Fixup", ForceTaxonFixup);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Force Locus Fixup", ForceLocusFixup);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  /*
  i = CommandItem (s, "Just Taxonomy Fixup", JustTaxonFixup);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  */
  i = CommandItem (s, "Correct CDS Genetic Codes", CorrectCDSGenCodes);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Update Replaced EC_numbers", UpdateECnumbers);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Parse GO Terms from User Fields", ParseGoTermsFromUserFieldsItem);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "NC_Cleanup", DoNCCleanup);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "PT_Cleanup", DoPTCleanup);
  SetObjectExtra (i, bfp, NULL);
  if (indexerVersion) {
    SeparatorItem (s);
    i = CommandItem (s, "SeriousSeqEntryCleanup", DoSSECleanup);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (s);
    i = CommandItem (s, "Normalize Descriptor Order", DoNormalizeDescs);
    SetObjectExtra (i, bfp, NULL);
  }  

  MakeSpecialLocusTagMenu (m, bfp);

  s = SubMenu (m, "Process/ P");
  i = CommandItem (s, "Process FASTA Proteins", ParseInMoreProteins);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Process Additional FASTA cDNAs", ParseInMoreMRNAs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Process FASTA Nucleotide Updates", ParseInNucUpdates);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Process Oligonucleotide PCR Primers", ParseInOligoPrimers);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Process Pseudo Misc_feats", ProcessPseudoMiscFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Process Scanned tRNA Anticodons", tRNAScanUpdate);
  SetObjectExtra (i, bfp, NULL);
    
  s = SubMenu (m, "Display");
  i = CommandItem (s, "View Sorted Protein FASTA", ViewSortedProteinMenuItem);
  SetObjectExtra (i, bfp, NULL);
  if (extraServices) {
    i = CommandItem (s, "Make ToolBar Window/ T", MakeToolBarWindow);
    SetObjectExtra (i, bfp, NULL);
  }
  i = CommandItem (s, "Alignment Summary", ViewAlignmentSummary);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Weighted Alignment Percent Identity", WeightedAlignmentPercentIdentityItem);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "VecScreenTool", VecScreenTool);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "On-Caller Tool", ShowOnCallerTool);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Print Test List", PrintDiscrepancyTestListMenuItem);
  SetObjectExtra (i, bfp, NULL);

  MakeSpecialLinkMenu (m, bfp);

  MakeSpecialSelectMenu (m, bfp);

  MakeSpecialDesktopMenu (m, bfp);

  SeparatorItem (m);
  s = SubMenu (m, "Sort Unique Count/ U");
  x = SubMenu (s, "By Group");
  i = CommandItem (x, "Without Sequence", SUCBNoSequenceProc);
  SetObjectExtra (i, bfp, NULL);  
  i = CommandItem (x, "With Sequence", SUCBProc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Original (by position)", SUCProc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Resorted (by frequency)/ N", SUCRProc);
  SetObjectExtra (i, bfp, NULL);
  x = SubMenu (s, "Clickable");
  i = CommandItem (x, "Without Sequence", SUCClickableNoSequenceProc);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (x, "With Sequence", SUCClickableProc);
  SetObjectExtra (i, bfp, NULL);

  if (indexerVersion) {
    SeparatorItem (m);
    s = SubMenu (m, "Export Alignment");
    i = CommandItem (s, "Interleave", ExportAlignmentInterleave);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (s, "Contiguous", ExportAlignmentContiguous);
    SetObjectExtra (i, bfp, NULL);

    SeparatorItem (m);
    i = CommandItem (m, "Test New Parse", TestNewParse);
    SetObjectExtra (i, bfp, NULL);

#if 0        
    SeparatorItem (m);
    s = SubMenu (m, "Test Taxon3 Server");
    i = CommandItem (s, "9606", TestTax3A);
    i = CommandItem (s, "E. coli", TestTax3B);
    SetObjectExtra (i, bfp, NULL);
#endif
    
#if 0
    SeparatorItem (m);
    i = CommandItem (m, "Add All Dbxrefs", AddAllDbxrefsMenuItem);
    SetObjectExtra (i, bfp, NULL);
#endif
    
  }
  SeparatorItem (m);
  i = CommandItem (m, "Discrepancy Report", ShowDiscrepancyReport);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (m);
  i = CommandItem (m, "Macro Editor", LaunchMacroEditor);
  SetObjectExtra (i, bfp, NULL);

  /*
  SeparatorItem (m);
  i = CommandItem (m, "Modernize PCR Primers", ModernizePCRs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (m, "Modernize Misc RNAs", ModernizeRNAs);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (m, "Modernize Gene Nomenclature", ModernizeGenes);
  SetObjectExtra (i, bfp, NULL);
  */


/*#ifdef INTERNAL_NCBI_SEQUIN*/
/*#ifdef NEW_TAXON_SERVICE*/
/*
  SeparatorItem (m);
  i = CommandItem (m, "Prepare TaxList", PrepareTaxListProc);
*/
/*#endif*/
/*#endif*/
}

extern void SetupBatchApplyMenu (MenU s, BaseFormPtr bfp)

{
  IteM  i;

  i = CommandItem (s, "Add CDS", ApplyCDS);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add RNA", ApplyRRNA);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add other Feature", ApplyImpFeat);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add Qualifier", PublicApplyGBQual);
  SetObjectExtra (i, bfp, NULL);
  SeparatorItem (s);
  i = CommandItem (s, "Add Source Qual", PublicApplySourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add CDS-Gene-Prot-mRNA Qual", PublicApplyCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Add RNA Qual", PublicApplyRNAQual);
  SetObjectExtra (i, bfp, NULL);
}

extern void SetupBatchEditMenu (MenU s, BaseFormPtr bfp)

{
  IteM  i;

  i = CommandItem (s, "Edit Qualifiers", PublicEditGBQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit Source Qual", PublicEditSourceQual);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit CDS-Gene-Prot Qual", PublicEditCDSGeneProt);
  SetObjectExtra (i, bfp, NULL);
  i = CommandItem (s, "Edit RNA Qual", PublicEditRNAQual);
  SetObjectExtra (i, bfp, NULL);
}

