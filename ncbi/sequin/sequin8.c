/*   sequin8.c
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
* File Name:  sequin8.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   2/3/98
*
* $Revision: 6.524 $
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
#include <objsub.h>
#include <valid.h>
#include <cdrgn.h>
#include <suggslp.h>
#include <toasn3.h>
#include <subutil.h>
#include <explore.h>
#include <medarch.h>
#include <medutil.h>
#include <tofasta.h>
#include <asn2gnbk.h>
#include <alignmgr2.h>
#include <spidey.h>
#include <blast.h>
#include <salpanel.h>
#include <seqpanel.h>
#include <edutil.h>
#include <tax3api.h>
#include <asn2gnbp.h> /* included for discrepancy report */
#include <algo/blast/api/twoseq_api.h>
#include <algo/blast/api/blast_seqalign.h>

#include <macrodlg.h>
#include <macroapi.h>
#include <alignval.h>

#define DEFLINE_MAX_LEN          380
#define TEXT_MAX_LEN             64
#define DEFLINE_MAX_GENENAME_LEN 64
#define ALL_FEATURES             255

typedef struct evidenceformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  PopuP          evdence;
  Uint2          exp_ev;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  ButtoN         case_insensitive;
  ButtoN         when_string_not_present;
} EvidenceFormData, PNTR 

EvidenceFormPtr;

typedef struct codebreakformdata {
  FEATURE_FORM_BLOCK
  PopuP  aminoAcidPopup;
  Char   currentCodonStr [4];
  TexT   codonText;
  ButtoN acceptButton;
} CodeBreakFormData, PNTR CodeBreakFormPtr;

static Boolean IsRealImpFeat (Uint2 subtype)

{
  if (subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref) return TRUE;
  if (subtype == FEATDEF_oriT) return TRUE;
  return FALSE;
}



static void BreakIntoAGroup (BioseqSetPtr parent, Uint1 _class, SeqEntryPtr list)

{
  BioseqSetPtr  bssp;
  Int2          count;
  SeqEntryPtr   sep;
  SeqEntryPtr   tmp;

  while (list != NULL) {
    bssp = BioseqSetNew ();
    if (bssp == NULL) return;
    bssp->_class = _class;
    sep = SeqEntryNew ();
    if (sep == NULL) return;
    sep->choice = 2;
    sep->data.ptrvalue = (Pointer) bssp;
    if (parent->seq_set != NULL) {
      tmp = parent->seq_set;
      while (tmp->next != NULL) {
        tmp = tmp->next;
      }
      tmp->next = sep;
    } else {
      parent->seq_set = sep;
    }
    bssp->seq_set = list;
    for (tmp = list, count = 0; tmp != NULL && count < 99; tmp = tmp->next, count++) continue;
    if (tmp != NULL) {
      list = tmp->next;
      tmp->next = NULL;
    } else {
      list = NULL;
    }
  }
}


extern Int2 LIBCALLBACK MakeGroupsOf200 (Pointer data)

{
  BioseqSetPtr      bssp;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  OMProcControlPtr  ompcp;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       sep;
  AsnIoPtr       aip;
  Uint1          _class;
  Int2           count;
  Char           file [FILENAME_MAX];
  SeqEntryPtr    list;
  SeqEntryPtr    next;
  Char           output [PATH_MAX];
  Char           path [PATH_MAX];
  CharPtr        ptr;
  SeqSubmitPtr   ssp;
  Char           str [FILENAME_MAX];
  SeqEntryPtr    tmp;
#ifdef WIN_MAC
  FILE           *f;
#endif

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  if (ompcp->input_itemtype != OBJ_BIOSEQSET) {
    Message (MSG_ERROR, "Must select Bioseq-set!");
    return OM_MSG_RET_ERROR;
  }
  bssp = (BioseqSetPtr) ompcp->input_data; 
  if (bssp == NULL) return OM_MSG_RET_ERROR;
  sep = SeqMgrGetSeqEntryForData (bssp);
  
  _class = bssp->_class;
  if (_class != 7 && _class != 13 && _class != 14 &&
      _class != 15 && _class != 16 && _class != 18) { 
    Message (MSG_ERROR, "Can only use this for GenBank, Mut, Pop, Phy, Eco, and WGS sets!");
    return OM_MSG_RET_ERROR;
  }

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  list = bssp->seq_set;
  bssp->seq_set = NULL;
  bssp->_class = 7;
  BreakIntoAGroup (bssp, _class, list);

  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  PropagateFromGenBankBioseqSet (sep, TRUE);

  if (parenttype == OBJ_SEQSUB) {
    if (GetOutputFileName (path, sizeof (path), "")) {
      ssp = (SeqSubmitPtr) parentptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sep = (SeqEntryPtr) ssp->data;
        ptr = StringRChr (path, DIRDELIMCHR);
        if (ptr != NULL) {
          ptr++;
          StringNCpy_0 (file, ptr, sizeof (file));
          *ptr = '\0';
          tmp = bssp->seq_set;
          count = 0;
          while (tmp != NULL) {
            next = tmp->next;
            tmp->next = NULL;
            ssp->data = (Pointer) tmp;
            StringCpy (output, path);
            count++;
            if (count < 10) {
              sprintf (str, "%s0%1d", file, (int) count);
            } else {
              sprintf (str, "%s%2d", file, (int) count);
            }
            FileBuildPath (output, NULL, str);
#ifdef WIN_MAC
            f = FileOpen (output, "r");
            if (f != NULL) {
              FileClose (f);
            } else {
              FileCreate (output, "TEXT", "ttxt");
            }
#endif
            aip = AsnIoOpen (output, "w");
            if (aip != NULL) {
              SeqSubmitAsnWrite (ssp, aip, NULL);
            }
            AsnIoClose (aip);
            tmp->next = next;
            tmp = next;
          }
          ssp->data = (Pointer) sep;
        }
      }
    }
  }

  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  Update ();  
  return OM_MSG_RET_DONE;
}

extern void ParseInNucUpdates (IteM i)

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
  Message (MSG_OK, "Not yet implemented");
}

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

static Int2 GeneticCodeFromCrp (CdRegionPtr crp)

{
  Int2            code;
  GeneticCodePtr  gcp;
  Char            name [256];
  ValNodePtr      tmp;

  code = 0;
  name [0] = '\0';
  gcp = crp->genetic_code;
  if (gcp != NULL) {
    tmp = (ValNodePtr) gcp->data.ptrvalue;
    for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      switch (tmp->choice) {
        case 1 :
          if (name [0] == '\0') {
            StringNCpy_0 (name, (CharPtr) tmp->data.ptrvalue, sizeof (name));
          }
          break;
        case 2 :
          code = tmp->data.intvalue;
          break;
        default :
          break;
      }
    }
    if (code == 0) {
      gcp = GeneticCodeFind (code, name);
      if (gcp != NULL) {
        for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
          switch (tmp->choice) {
            case 2 :
              code = tmp->data.intvalue;
              break;
            default :
              break;
          }
        }
      }
    }
  }
  return code;
}

extern void ExtendSeqLocToPosition (SeqLocPtr slp, Boolean end5, Int4 pos)
{
  Uint1          strand;
  SeqLocPtr      slp_to_change, slp_index;
  Int4           extent_to_change;
  Int4           start, stop;
  SeqIdPtr       sip;
  BioseqPtr      bsp;
  
  if (slp == NULL || pos < 0) return;
  
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) return;

  slp_to_change = NULL;
  strand = SeqLocStrand (slp);
  switch (slp->choice)
  {
    case SEQLOC_INT:
      slp_to_change = slp;
      break;
    case SEQLOC_MIX:
  	case SEQLOC_PACKED_INT:
      sip = SeqLocId (slp);
      if (sip == NULL) return; /* can only process if all on one bioseq */
      slp_to_change = NULL;
      if ((strand == Seq_strand_minus && end5)
        || (strand != Seq_strand_minus && !end5))
      {
        extent_to_change = 0;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          stop = GetOffsetInBioseq (slp_index, bsp, SEQLOC_STOP);
          if (stop > extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = stop;
          }
        }
      }
      else
      {
        extent_to_change = bsp->length;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          start = GetOffsetInBioseq (slp_index, bsp, SEQLOC_START);
          if (start < extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = start;
          }
        }
      }
      break;
  }

  if (slp_to_change != NULL)
  {
    if ((strand == Seq_strand_minus && end5)
      || (strand != Seq_strand_minus && !end5))
    {
      start = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_START);
      stop = pos;
    }
    else
    {
      start = pos;
      stop = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_STOP);
    }
    if (start < 0 
        || stop > bsp->length - 1
        || start > stop)
    {
      return;
    }
    expand_seq_loc (start, stop, strand, slp_to_change);
  }
}

static void ExtendOnePartialFeatureEx (SeqFeatPtr sfp, Boolean extend5, Boolean extend3)
{
  BioseqPtr   bsp;
  Boolean     partial3, partial5;
  Int4        start_diff;
  CdRegionPtr crp;

  if (sfp == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5 && extend5)
  {
    start_diff = ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
    if (start_diff > 0 && sfp->data.choice == SEQFEAT_CDREGION) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
          if (crp->frame == 0) {
              crp->frame = 1;
          }
          crp->frame = (crp->frame + start_diff - 1) % 3 + 1;
      }
    }
  }
  if (partial3 && extend3)
  {
    ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
  }
}

static void ExtendOnePartialFeature (SeqFeatPtr sfp, Pointer userdata)
{
  ExtendOnePartialFeatureEx (sfp, TRUE, TRUE);
}

extern void ExtendPartialFeatures (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep, old_scope;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sel = ObjMgrGetSelected ();
  WatchCursor ();
  Update ();
  old_scope = SeqEntrySetScope(sep);
  if (sel == NULL)
  {
    VisitFeaturesInSep (sep, NULL, ExtendOnePartialFeature);
  }
  else
  {
    while (sel != NULL)
    {
      if (sel->entityID == bfp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if (sfp != NULL)
        {
          ExtendOnePartialFeature (sfp, NULL);
        }
      }
      sel = sel->next;
    }
  }
  SeqEntrySetScope(old_scope);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct extendpartialfeaturesform {
  FORM_MESSAGE_BLOCK
  DialoG feature_type;
  ButtoN extend5;
  ButtoN extend3;
  DialoG string_constraint;
  ButtoN leave_dlg_up;
} ExtendPartialFeaturesFormData, PNTR ExtendPartialFeaturesFormPtr;

static void DoExtendPartialFeatures (ButtoN b)
{
  ExtendPartialFeaturesFormPtr f;
  SeqEntryPtr sep;
  AECRActionPtr action;
  ApplyActionPtr fake_apply;
  FeatureFieldPtr feature_field;
  ValNodePtr vnp;
  StringConstraintPtr scp;
  ValNodePtr object_list;
  Boolean    extend5, extend3;

  f = (ExtendPartialFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  sep = GetTopSeqEntryForEntityID (f->input_entityID);
  if (sep == NULL) return;

  /* note - this is a small amount of hackery, designed to put off adding "extend partial features" as a macro language action.
   * if this is ever added, this code should be moved to macroapi.c.
   */
  feature_field = FeatureFieldNew ();
  vnp = DialogToPointer (f->feature_type);
  if (vnp == NULL) {
    feature_field->type = Feature_type_any;
  } else {
    feature_field->type = vnp->choice;
    vnp = ValNodeFree (vnp);
  }
  fake_apply = ApplyActionNew ();
  ValNodeAddPointer (&(fake_apply->field), FieldType_feature_field, feature_field);
  action = AECRActionNew ();
  ValNodeAddPointer (&(action->action), ActionChoice_apply, fake_apply);
  
  scp = DialogToPointer (f->string_constraint);
  if (scp != NULL) {
    ValNodeAddPointer (&(action->constraint), ConstraintChoice_string, scp);
  }
  object_list = GetObjectListForAECRAction (sep, action);
  action = AECRActionFree (action);

  extend5 = GetStatus (f->extend5);
  extend3 = GetStatus (f->extend3);
  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      ExtendOnePartialFeatureEx (vnp->data.ptrvalue, extend5, extend3);
    }
  }
  object_list = ValNodeFree (object_list);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (f->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, f->input_entityID, 0, 0);
  if (!GetStatus (f->leave_dlg_up)) {
    Remove (f->form);
  }
}


extern void ExtendPartialFeaturesWithConstraint (IteM i)
{
  BaseFormPtr       bfp;
  ExtendPartialFeaturesFormPtr f;
  ValNodePtr          feature_list = NULL;
  ValNode             vn;
  WindoW              w;
  GrouP               h, g, c;
  PrompT              p1, p2;
  ButtoN              b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  f = (ExtendPartialFeaturesFormPtr) MemNew (sizeof (ExtendPartialFeaturesFormData));
  if (f == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Extend Partial Features", StdCloseWindowProc);
  SetObjectExtra (w, f, StdCleanupExtraProc);
  f->form = (ForM) w;
  f->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  p1 = StaticPrompt (h, "Feature Type to Extend", 0, dialogTextHeight, programFont, 'c');
  ValNodeAddPointer (&feature_list, Feature_type_any, StringSave ("Any"));
  AddAllFeaturesToChoiceList (&feature_list);
  
  f->feature_type = ValNodeSelectionDialog (h, feature_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                NULL, NULL, FALSE);
  vn.choice = Feature_type_any;
  vn.data.ptrvalue = NULL;
  vn.next = NULL;
  PointerToDialog (f->feature_type, &vn);

  g = HiddenGroup (h, 2, 0, NULL);
  f->extend5 = CheckBox (g, "Extend partial 5'", NULL);
  SetStatus (f->extend5, TRUE);
  f->extend3 = CheckBox (g, "Extend partial 3'", NULL);
  SetStatus (f->extend3, TRUE);
  
  p2 = StaticPrompt (h, "Optional Constraint", 0, dialogTextHeight, programFont, 'c');
  f->string_constraint = StringConstraintDialog (h, "Where feature text", FALSE, NULL, NULL);
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", DoExtendPartialFeatures);
  SetObjectExtra (b, f, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  f->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) p1,
                              (HANDLE) f->feature_type,
                              (HANDLE) g,
                              (HANDLE) p2,
                              (HANDLE) f->string_constraint,
                              (HANDLE) c,
                              NULL);
  Show (w);
}


static Boolean HasValidStartCodon (SeqFeatPtr cds)
{
  ByteStorePtr  bs;
  CharPtr       prot;

  if (cds == NULL) return FALSE;

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);

  if (bs == NULL) return FALSE;
  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot == NULL) return FALSE;
  if (prot [0] != 'M') return FALSE;
  return TRUE;
}

static void FixReadingFrame (
  CdRegionPtr crp,
  SeqLocPtr slp,
  BioseqPtr bsp,
  Int4 start
)
{
  Int4 offset;

  if (crp == NULL || slp == NULL) return;

  if (SeqLocStrand (slp) == Seq_strand_minus)
  {
    offset = bsp->length - start;
  }
  else
  {
    offset = start;
  }
  start = offset % 3;
  start = start + 1;
  crp->frame = start;
}

extern void RecomputeSuggestedIntervalsForCDS 
(Uint2          entityID,
 BioseqPtr PNTR batchbsp,
 Int4Ptr        count,
 MonitorPtr     mon,
 SeqFeatPtr     sfp)
{
  Int2           code;
  CdRegionPtr    crp;
  BioseqPtr      nucbsp;
  BioseqPtr      protbsp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  Char           str [256];
  Char           tmp [256];
  Boolean        partial3, partial5;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return;

  code = GeneticCodeFromCrp (crp);

  nucbsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
  if (nucbsp != NULL && batchbsp != NULL && *batchbsp != NULL 
      && nucbsp != *batchbsp) {
    ClearBatchSuggestNucleotide ();
    *batchbsp = nucbsp;
    SetBatchSuggestNucleotide (*batchbsp, code);
/*    Message (MSG_POSTERR, "Recompute Suggest is reverting to slower processing"); */
  }
  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    protbsp = BioseqFind (sip);
    if (nucbsp != NULL && protbsp != NULL &&
        ISA_na (nucbsp->mol) && ISA_aa (protbsp->mol) &&
        nucbsp->length > 0 && protbsp->length > 0) {
      str [0] = '\0';
      tmp [0] = '\0';
      sip = SeqIdFindWorst (protbsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      if (count != NULL)
      {
        (*count) ++;
        if (mon != NULL)
        {
          sprintf (str, "Processing sequence %d [%s]", *count, tmp);
          MonitorStrValue (mon, str);
          Update ();
        }
      }
      slp = PredictCodingRegion (nucbsp, protbsp, code);
      if (slp == NULL) return;

      /* correct for partial conditions */
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp;

      /* if no valid start codon, cds is 5' partial */
      if (! HasValidStartCodon (sfp))
      {
        partial5 = TRUE;
      }

      if (partial5 || partial3)
      {
#if 0 
        /* removed, per DeAnne's request */     
        if (partial5)
        {
          start = GetOffsetInBioseq (sfp->location, nucbsp, SEQLOC_START);
          FixReadingFrame (crp, sfp->location, nucbsp, start);
          ExtendSeqLocToEnd (sfp->location, nucbsp, TRUE);
        }
#endif        
        if (partial3)
        {
          ExtendSeqLocToEnd (sfp->location, nucbsp, FALSE);
        }
        SetSeqLocPartial (sfp->location, partial5, partial3);
      }
      sfp->partial = LocationHasNullsBetween (sfp->location);
      sfp->partial |= partial5 || partial3;
    }
  }
}

extern void RecomputeIntervalsForOneCDS (SeqFeatPtr sfp, RecompDataPtr rdp)
{
  SeqFeatPtr gene_to_update = NULL;
  SeqLocPtr      orig_loc = NULL;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || rdp == NULL)
  {
    return;
  }

  if (rdp->fix_genes)
  {
    gene_to_update = SeqMgrGetOverlappingGene (sfp->location, NULL);
    orig_loc = (SeqLocPtr) AsnIoMemCopy (sfp->location, 
                                         (AsnReadFunc) SeqLocAsnRead,
                                         (AsnWriteFunc) SeqLocAsnWrite);
  }
  RecomputeSuggestedIntervalsForCDS (rdp->entityID, &(rdp->batchbsp),
                                     &(rdp->count), rdp->mon, sfp);

  if (gene_to_update != NULL)
  {
    UpdateGeneLocation (gene_to_update, orig_loc,
                        sfp->location, rdp->entityID);
  }
  orig_loc = SeqLocFree (orig_loc);
  
}

static Boolean RecomputeSuggCallback (GatherContextPtr gcp)
{
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;

  RecomputeIntervalsForOneCDS (sfp, rdp);
  return TRUE;
}

extern void RecomputeSuggestEx (Uint2 entityID, Boolean fix_genes, Boolean recompute_all)

{
  Int2              code;
  GatherScope       gs;
  SeqEntryPtr       nucsep;
  RecompData        rd;
  SeqEntryPtr       sep;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  sel = ObjMgrGetSelected ();
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  rd.batchbsp = NULL;
  rd.no_stop_at_end_of_complete_cds = FALSE;
  rd.fix_genes = fix_genes;
  rd.entityID = entityID;
  nucsep = FindNucSeqEntry (sep);
  if (nucsep != NULL && IS_Bioseq (nucsep)) {
    rd.batchbsp = (BioseqPtr) nucsep->data.ptrvalue;
  }
  if (rd.batchbsp != NULL) {
    code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
    SetBatchSuggestNucleotide (rd.batchbsp, code);
  }
  if (sel == NULL || recompute_all)
  {
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    GatherEntity (entityID, (Pointer) (&rd), RecomputeSuggCallback, &gs);
  }
  else
  {
    while (sel != NULL)
    {
      if (sel->entityID == entityID
        && sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_CDS)
        {
          RecomputeIntervalsForOneCDS (sfp, &rd);
        }
      }
      sel = sel->next;
    }
  }
  MonitorFree (rd.mon);
  if (rd.batchbsp != NULL) {
    ClearBatchSuggestNucleotide ();
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

extern void RecomputeSuggest (IteM i)
{
  BaseFormPtr       bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
  {
    return;
  }
  RecomputeSuggestEx (bfp->input_entityID, FALSE, FALSE);
}

extern void RecomputeSuggestFixGenes (IteM i)
{
  BaseFormPtr       bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
  {
    return;
  }
  RecomputeSuggestEx (bfp->input_entityID, TRUE, FALSE);
}


static Boolean RetranslateCDSCallback (GatherContextPtr gcp)

{
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  return RetranslateOneCDS (sfp, gcp->entityID, rdp->include_stop,
                            rdp->no_stop_at_end_of_complete_cds);
}

extern void RetranslateCdRegionsEx (
  Uint2   entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds )

{
  GatherScope  gs;
  RecompData   rd;
  SeqEntryPtr  sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  rd.include_stop = include_stop;
  rd.no_stop_at_end_of_complete_cds = no_stop_at_end_of_complete_cds;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (entityID, (Pointer) (&rd), RetranslateCDSCallback, &gs);
  MonitorFree (rd.mon);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void RetranslateCdRegions (
  IteM i,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds )

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  RetranslateCdRegionsEx (bfp->input_entityID, include_stop, no_stop_at_end_of_complete_cds);
}

extern void RetranslateCdRegionsNoStop (IteM i)

{
  RetranslateCdRegions (i, FALSE, FALSE);
}

extern void RetranslateCdRegionsDoStop (IteM i)

{
  RetranslateCdRegions (i, TRUE, FALSE);
}

extern void RetranslateCdRegionsNoStopExceptEndCompleteCDS (IteM i)
{
  RetranslateCdRegions (i, TRUE, TRUE);
}


static void DoReprocessPeptides (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  MolInfoPtr    mip;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  ProtRefPtr    prp;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;

  if (sfp->data.choice != SEQFEAT_PROT) return;
  if (sfp->product == NULL) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;
  if (prp->processed < 1 || prp->processed > 4) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
    bestprot = FindBestProtein (sfp->idx.entityID, sfp->product);
    prot = GetSequenceByFeature (sfp);
    if (prot == NULL) return;
    ptr = prot;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
    bs = BSNew (1000);
    if (bs != NULL) {
      ptr = prot;
      BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
    }
    bsp->repr = Seq_repr_raw;
    bsp->mol = Seq_mol_aa;
    bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
    bsp->seq_data = (SeqDataPtr) bs;
    bsp->seq_data_type = Seq_code_ncbieaa;
    bsp->length = BSLen (bs);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    if (bestprot != NULL) {
      bestprot->location = SeqLocFree (bestprot->location);
      bestprot->location = CreateWholeInterval (sep);
      SetSeqLocPartial (bestprot->location, partial5, partial3);
      bestprot->partial = (partial5 || partial3);
    }
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

extern void ReprocessPeptideProducts (IteM i);
extern void ReprocessPeptideProducts (IteM i)

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
  VisitFeaturesInSep (sep, NULL, DoReprocessPeptides);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoReprocessMrnas (SeqFeatPtr sfp, Pointer userdata)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  MolInfoPtr    mip;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->product == NULL) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;
  if (rrp->type != 2) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_na (bsp->mol) && bsp->repr == Seq_repr_raw) {
    prot = GetSequenceByFeature (sfp);
    if (prot == NULL) return;
    ptr = prot;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
    bs = BSNew (1000);
    if (bs != NULL) {
      ptr = prot;
      BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
    }
    bsp->repr = Seq_repr_raw;
    bsp->mol = Seq_mol_na;
    bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
    bsp->seq_data = (SeqDataPtr) bs;
    bsp->seq_data_type = Seq_code_iupacna;
    bsp->length = BSLen (bs);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
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

extern void ReprocessmRNAProducts (IteM i);
extern void ReprocessmRNAProducts (IteM i)

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
  VisitFeaturesInSep (sep, NULL, DoReprocessMrnas);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

/*
 * ApplyCodeBreakToCDS
 */
 
static Boolean ApplyCodeBreakToCDS (SeqFeatPtr sfp, CharPtr codonStr, Uint1 aaNum)
{
  Uint1            aaChar;
  SeqLocPtr        aaSlp;
  Int4             aaPosition;
  SeqPntPtr        aaSpp;
  CharPtr          basePtr;
  CharPtr          bases;
  CodeBreakPtr     cbp;
  CdRegionPtr      crp;
  Int4             dnaLen;
  SeqLocPtr        dnaSlp;
  CodeBreakPtr     lastCbp;
  SeqCodeTablePtr  sctp;
  Int4             total;
  Boolean          added_code_breaks = FALSE;
  
  if (sfp == NULL || codonStr == NULL)
  {
    return FALSE;
  }

  /* Get the nucleotide sequence */

  dnaLen = SeqLocLen (sfp->location);
  if (dnaLen < 1)
    return FALSE;

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;

  bases = ReadCodingRegionBases (sfp->location, dnaLen, crp->frame, &total);

  /* Search for the selected codon in the */
  /* nucleotide sequence.  If found, add  */
  /* it as a codebreak.                   */

  basePtr = bases;
  aaPosition = 0;
  while (basePtr[0] != '\0') {
    if (StringNCmp (basePtr, codonStr, 3) == 0) {

      /* Create a new seq point object with the aa location */

      aaSpp = SeqPntNew ();
      aaSpp->point  = aaPosition;
      aaSpp->strand = Seq_strand_plus;
      aaSpp->id = SeqLocId (sfp->product);
      aaSpp->fuzz   = NULL;

      /* Make a SeqLoc using the seq point */

      aaSlp = (SeqLocPtr) ValNodeNew (NULL);
      aaSlp->choice = SEQLOC_PNT;
      aaSlp->data.ptrvalue = (Pointer) aaSpp;

      /* Convert the seqloc to a DNA location */

      dnaSlp = aaLoc_to_dnaLoc (sfp, aaSlp);

      /* Create the code break using the DNA location */

      cbp = CodeBreakNew ();
      cbp->loc = dnaSlp;
      sctp = SeqCodeTableFind (Seq_code_ncbieaa);
      aaChar = (Uint1) GetSymbolForResidue (sctp, aaNum);
      cbp->aa.value.intvalue = aaChar;
      cbp->aa.choice = 1; /* ncbieaa */

      /* Insert the code break into the CDS's */
      /* existing list of code breaks.        */

      lastCbp = crp->code_break;
      if (lastCbp == NULL)
      {
	      crp->code_break = cbp;        
      }
      else 
      {
        while (lastCbp->next != NULL)
        {
          lastCbp = lastCbp->next;          
        }
	      lastCbp->next = cbp;
	      cbp->next = NULL;
      }

      added_code_breaks = TRUE;
    }
    basePtr += 3;
    aaPosition++;
  }

  return added_code_breaks;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyCodeBreak_FeatureCallback () -- Called for each CDS feature in */
/*                                      a Bioseq.  Checks for any      */
/*                                      nucleotide triplets that match */
/*                                      the one in the given code      */
/*                                      break and sets a code break    */
/*                                      for each one that is found.    */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK ApplyCodeBreak_FeatureCallback (SeqFeatPtr sfp,
					SeqMgrFeatContextPtr fcontext)
{
  Uint1            aaNum;
  CodeBreakFormPtr cbfp;
  Char             codonStr [4];
  Int2             i;

  cbfp = (CodeBreakFormPtr) fcontext->userdata;


  /* Get the selected Amino Acid and codon triplet */

  GetTitle (cbfp->codonText, codonStr, sizeof (codonStr));
  for (i = 0; i < 3; i++)
    codonStr [i] = TO_UPPER (codonStr [i]);

  aaNum = (Uint1) GetValue (cbfp->aminoAcidPopup);
  aaNum += 63;

  /*
  if (aaNum >= 74)
  {
  	aaNum++;
  }
  if (aaNum >= 79)
  {
  	aaNum++;
  }
  */

  if (ApplyCodeBreakToCDS (sfp, codonStr, aaNum))
  {
    /* Retranslate the CDS */

    RetranslateOneCDS (sfp, fcontext->entityID, TRUE, FALSE);
    
  }
  
  /* Return TRUE to continue on to the next CDS feature */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyCodeBreak_BioseqCallback () -- Called by SeqMgrExploreBioseqs  */
/*                                     for each Bioseq.  Searches the  */
/*                                     Bioseq for CDS features and adds*/
/*                                     the given code break to any     */
/*                                     found.                          */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK ApplyCodeBreak_BioseqCallback (BioseqPtr bsp,
					 SeqMgrBioseqContextPtr bcontext)
{
  Boolean featureFilterArray [SEQFEAT_MAX];

  /* Set up to explore only CDS features */

  MemSet ((Pointer) (featureFilterArray),
	  (int) FALSE,
	  SEQFEAT_MAX);

  featureFilterArray[SEQFEAT_CDREGION] = TRUE;

  /* Explore the Bioseq's CDS features, marking the */
  /* ones with internal stop codons as pseudo.      */

  SeqMgrExploreFeatures (bsp, bcontext->userdata,
			 ApplyCodeBreak_FeatureCallback, NULL,
			 featureFilterArray, NULL);

  /* Return TRUE to continue on to the next Bioseq */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* DoAddCodeBreak_Callback () -- Called when the 'Apply' button is     */
/*                               pressed in the "Add Code Break"       */
/*                               window.  Adds the entered code break  */
/*                               to all CDS features.                  */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void DoAddCodeBreak_Callback (ButtoN b)
{
  CodeBreakFormPtr  cbfp;

  cbfp = (CodeBreakFormPtr) GetObjectExtra (b);

  /* Change to the "working" cursor */

  Hide (cbfp->form);
  WatchCursor ();
  Update ();

  /* Visit all the Bioseqs, where we will */
  /* then explore their CDS features.     */

  SeqMgrExploreBioseqs (cbfp->input_entityID, NULL, (Pointer) cbfp,
			ApplyCodeBreak_BioseqCallback, TRUE, FALSE, TRUE);

  /* Restore the cursor and force an update */

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (cbfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cbfp->input_entityID, 0, 0);
  Remove (cbfp->form);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* PopulateAAPopup () -- Creates a popup list of amino acids.          */
/*                                                                     */
/*     NOTE : This function is identical to (and identically named as) */
/*            a function in cdrgn.c                                    */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void PopulateAAPopup (PopuP AAitem)

{
  Char             ch;
  Uint1            first;
  Uint1            i;
  Char             item [77];
  Uint1            last;
  SeqCodeTablePtr  sctp;
  CharPtr          str;

  sctp = SeqCodeTableFind (Seq_code_ncbieaa);
  first = FirstResidueInCode (sctp);
  last = LastResidueInCode (sctp);
  PopupItem (AAitem, " ");
  for (i = 65; i <= last; i++) {
    /*
    if (i == 74 || i == 79) {
      continue;
    }
    */
    ch = GetSymbolForResidue (sctp, i);
    str = (CharPtr) GetNameForResidue (sctp, i);
    sprintf (item, "%c    %s", ch, str);
    PopupItem (AAitem, item);
  } 
  SetValue (AAitem, 1); 
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* IsLegalCodon () - Determines if a three base string is a legal      */
/*                   codon.                                            */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean IsLegalCodon (CharPtr codonStr)
{
  Int2 i;
  Char baseChar;

  /* Only allow three characters */

  if (StringLen (codonStr) > 3)
    return FALSE;

  /* Allow only the character A,C,G,T,U and */
  /* convert the U to a T.                  */

  i = 0;
  while (i < 3) {

    if (codonStr [i] == '\0')
      break;

    baseChar = codonStr [i];
    
    if (StringChr ("acgtuACGTU", baseChar) == NULL)
      return FALSE;
    if ('U' == baseChar)
      codonStr [i] = 'T';
    else if ('u' == baseChar)
      codonStr [i] = 't';
    
    i++;
  }

  /* If we made it this far, it's a valid codon */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* CodonText_Callback () -- Called whenever a keystoke is entered in   */
/*                          Codon text field.  Validates to see if the */
/*                          keystroke should be allowed.               */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void CodonText_Callback (TexT codonText)

{
  CodeBreakFormPtr  cbfp;
  Int2              aaNum;
  Char              newCodonStr [5];

  /* Get the currect code break data */

  cbfp = (CodeBreakFormPtr) GetObjectExtra (codonText);
  if (cbfp == NULL)
    return;

  /* If the new codon string is not legal */
  /* then reset to the previous text.     */

  GetTitle (codonText, newCodonStr, sizeof (newCodonStr));

  if (!IsLegalCodon (newCodonStr))
    StringCpy (newCodonStr, cbfp->currentCodonStr);
  else
    StringCpy (cbfp->currentCodonStr, newCodonStr);

  SafeSetTitle (cbfp->codonText, newCodonStr);

  /* Only enable the accept button if */
  /* we have a full codon.            */

  if (StringLen (newCodonStr) != 3) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* See if an amino acid has been selected yet */

  aaNum = GetValue (cbfp->aminoAcidPopup);
  if (aaNum <= 1) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* If we made it this far then we have both a codon and */
  /* an amino acid, so enable the accept button.          */
 
  SafeEnable (cbfp->acceptButton);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* SelectAminoAcid_Callback () -- Called whenever a new amino acid is  */
/*                                selected in the Amino Acid Popup.    */
/*                                Toggles 'Accept' button base on      */
/*                                current state.                       */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void SelectAminoAcid_Callback (PopuP p)
{
  CodeBreakFormPtr  cbfp;
  Char              codonStr [4];
  Int2              aaNum;

  /* Get the currect code break data */

  cbfp = (CodeBreakFormPtr) GetObjectExtra (p);
  if (cbfp == NULL)
    return;

  /* Only enable the accept button if */
  /* we have a full codon.            */

  GetTitle (cbfp->codonText, codonStr, sizeof (codonStr));
  if (StringLen (codonStr) != 3) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* Get the newly selected amino acid */

  aaNum = GetValue (cbfp->aminoAcidPopup);

  /* If an amino acid is selected then */
  /* enable the accept button.         */

  if (aaNum > 1)
    SafeEnable (cbfp->acceptButton);
  else
    SafeDisable (cbfp->acceptButton);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* AddGlobalCodeBreak () -- Gets a nucleotide triplet and an amino     */
/*                          acid from the user and adds them as        */
/*                          codebreaks for all CDS features.           */
/*                                                                     */
/*---------------------------------------------------------------------*/

extern void AddGlobalCodeBreak (IteM i);
extern void AddGlobalCodeBreak (IteM i)

{
  BaseFormPtr      bfp;
  WindoW           breakWin;
  GrouP            mainGroup;
  GrouP            buttGroup;
  CodeBreakFormPtr cbfp;

  /* Get the current state of things */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  cbfp = (CodeBreakFormPtr) MemNew (sizeof (EvidenceFormData));

  /* Create a window to get the codon and */
  /* the Amino acid from the user.        */

  breakWin = FixedWindow (-50, -33, -10, -10,
			  "Add Code Break", StdCloseWindowProc);
  SetObjectExtra (breakWin, cbfp, StdCleanupFormProc);
  cbfp->form = (ForM) breakWin;
  cbfp->formmessage = NULL;
  cbfp->input_entityID = bfp->input_entityID;

  mainGroup = HiddenGroup (breakWin, -2, 0, NULL);

  /* Create a text entry box for the nucl. codon */

  StaticPrompt (mainGroup, "Triplet Codon", 0, popupMenuHeight,
		programFont, 'l');
  cbfp->codonText = DialogText (mainGroup, "", 3, CodonText_Callback);
  SetObjectExtra (cbfp->codonText, cbfp, NULL);
  cbfp->currentCodonStr [0] = '\0';

  /* Add a Popup list of Amino Acids */

  StaticPrompt (mainGroup, "Amino Acid", 0, popupMenuHeight,
		programFont, 'l');
  cbfp->aminoAcidPopup = PopupList (mainGroup, TRUE,
				    SelectAminoAcid_Callback);
  PopulateAAPopup (cbfp->aminoAcidPopup);
  SetObjectExtra (cbfp->aminoAcidPopup, cbfp, NULL);

  /* Add Accept and Cancel buttons */

  buttGroup = HiddenGroup (breakWin, 2, 0, NULL);
  cbfp->acceptButton = DefaultButton (buttGroup, "Accept",
				   DoAddCodeBreak_Callback);
  SetObjectExtra (cbfp->acceptButton, cbfp, NULL);
  SafeDisable (cbfp->acceptButton);
  PushButton (buttGroup, "Cancel", StdCancelButtonProc);

  /* Line things up and display the window */

  AlignObjects (ALIGN_CENTER, (HANDLE) mainGroup, (HANDLE) buttGroup, NULL);
  RealizeWindow (breakWin);
  Show (breakWin);
  Update ();

}

static void ParseCodonQualToCodeBreakCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqCodeTablePtr sctp;
  GBQualPtr       gqual, prev_qual = NULL, next_qual;
  CharPtr         cp;
  Char            codon_text[4];
  Char            symbol_text [4];
  Uint1           aaNum;
  Int4            i;
  Boolean         converted_qual;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL)
  {
    return;
  }
  
  sctp = (SeqCodeTablePtr) userdata;
  
  for (gqual = sfp->qual; gqual != NULL; gqual = next_qual)
  {
    next_qual = gqual->next;
    converted_qual = FALSE;
    if (StringCmp (gqual->qual, "codon") == 0)
    {
      cp = StringSearch (gqual->val, "seq:\"");
      if (cp != NULL)
      {
        cp += 5;
        StringNCpy (codon_text, cp, 3);
        codon_text [3] = 0;
        for (i = 0; i < 3; i++)
          codon_text [i] = TO_UPPER (codon_text [i]);

        cp = StrChr (cp, ':');
        if (cp != NULL)
        {
          cp++;
          StringNCpy (symbol_text, cp, 3);
          symbol_text [3] = 0;
          aaNum = FindResidueByName (symbol_text, sctp);         
          if (ApplyCodeBreakToCDS (sfp, codon_text, aaNum))
          {
            /* Retranslate the CDS */

            RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, FALSE);
            
            /* remove the codon qual */
            if (prev_qual == NULL)
            {
              sfp->qual = gqual->next;
            }
            else
            {
              prev_qual->next = gqual->next;
            }
            gqual->next = NULL;
            GBQualFree (gqual);
            converted_qual = TRUE;
          }
        }
      }
    }
    if (!converted_qual)
    {
      prev_qual = gqual;
    }
  }
}

extern void ParseCodonQualToCodeBreak (IteM i)
{
  BaseFormPtr      bfp;
  SeqEntryPtr      sep;
  SeqCodeTablePtr  sctp;
  
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
  
  sctp = SeqCodeTableFind (Seq_code_ncbieaa);
  
  VisitFeaturesInSep (sep, sctp, ParseCodonQualToCodeBreakCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

static void CorrectGenCodeIndexedCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  Int2Ptr         pGenCode;
  ValNodePtr      vnp;
  Boolean         need_replacement = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION 
      || sfp->data.value.ptrvalue == NULL
      || userdata == NULL) return;
 
  pGenCode = (Int2Ptr) userdata;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp->genetic_code != NULL
      && crp->genetic_code->choice == 254) {
    if (crp->genetic_code->data.ptrvalue == NULL) {
      vnp = ValNodeNew (NULL);
      vnp->choice = 2;
      vnp->data.intvalue = (Int4) *pGenCode;
    } else {
      vnp = crp->genetic_code->data.ptrvalue;
      if (vnp->next == NULL && vnp->choice == 2) {
        vnp->data.intvalue = (Int4) *pGenCode;
      } else {
        need_replacement = TRUE;
      }
    }
  } else {
    need_replacement = TRUE;
  }
  if (need_replacement) {
    gc = GeneticCodeNew ();
    if (gc == NULL) return;
    crp->genetic_code = GeneticCodeFree (crp->genetic_code);
    vnp = ValNodeNew (NULL);
    gc->data.ptrvalue = vnp;
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = (Int4) *pGenCode;
    }
    crp->genetic_code = gc;
  }
}

static void CorrectGenCodesBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;

  if (bsp == NULL || userdata == NULL) return;
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext)) {
    CorrectGenCodeIndexedCallback (sfp, userdata);
  }

}

typedef struct gencodescan {
  Boolean mito;
  Boolean plastid;
  Int2    nuclCode;
  Int2    mitoCode;
  Boolean already_found;
} GenCodeScanData, PNTR GenCodeScanPtr;

static void JustGetGenCodeFromOrgRef (OrgRefPtr orp, GenCodeScanPtr gp)
{
  if (orp == NULL || orp->orgname == NULL || gp == NULL || gp->already_found) return;

  gp->nuclCode = orp->orgname->gcode;
  gp->mitoCode = orp->orgname->mgcode;
}

static void JustGetGenCodeFromBiop (BioSourcePtr biop, GenCodeScanPtr gp)
{
  if (biop == NULL || gp == NULL) return;
  if (gp->already_found && !biop->is_focus) return;

  gp->mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                        biop->genome == GENOME_mitochondrion ||
                        biop->genome == GENOME_hydrogenosome);

  gp->plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                                biop->genome == GENOME_chromoplast ||
                                biop->genome == GENOME_plastid ||
                                biop->genome == GENOME_cyanelle ||
                                biop->genome == GENOME_apicoplast ||
                                biop->genome == GENOME_leucoplast ||
                                biop->genome == GENOME_proplastid);

  JustGetGenCodeFromOrgRef (biop->org, gp);
  gp->already_found = TRUE;
}


static void JustGetGenCodeFromFeat (SeqFeatPtr sfp, Pointer userdata) 
{
  GenCodeScanPtr gp;

  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;

  gp = (GenCodeScanPtr) userdata;

  JustGetGenCodeFromBiop (sfp->data.value.ptrvalue, gp);
}

static void JustGetGenCodeFromDesc (SeqDescrPtr sdp, Pointer userdata)
{
  GenCodeScanPtr gp;

  if (sdp == NULL || userdata == NULL || sdp->choice != Seq_descr_source) return;

  gp = (GenCodeScanPtr) userdata;

  JustGetGenCodeFromBiop (sdp->data.ptrvalue, gp);
}

static Int2 JustGetGenCodeForSeqEntry (SeqEntryPtr sep) 
{
  GenCodeScanData gd;

  gd.already_found = FALSE;
  gd.mito = FALSE;
  gd.mitoCode = 0;
  gd.nuclCode = 0;
  gd.plastid = FALSE;

  VisitDescriptorsInSep (sep, &gd, JustGetGenCodeFromDesc);
  VisitFeaturesInSep (sep, &gd, JustGetGenCodeFromFeat);

  if (gd.plastid) {
    return 11;
  } else if (gd.mito) {
    return gd.mitoCode;
  } else {
    return gd.nuclCode;
  }
}


extern void CorrectGenCodes (SeqEntryPtr sep, Uint2 entityID)

{
  BioseqSetPtr  bssp;
  Int2          genCode;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CorrectGenCodes (sep, entityID);
      }
      return;
    }
  }

  genCode = JustGetGenCodeForSeqEntry(sep);
  VisitFeaturesInSep (sep, &genCode, CorrectGenCodeIndexedCallback);
  VisitBioseqsInSep (sep, &genCode, CorrectGenCodesBioseqCallback);
}

extern void CorrectCDSGenCodes (IteM i)

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
  CorrectGenCodes (sep, bfp->input_entityID);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static BioseqPtr SqnGetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    tmp;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else {
    tmp = SeqLocFindNext (slp, NULL);
    if (tmp != NULL) {
      sip = SeqLocId (tmp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = SeqMgrGetSeqEntryForData (bsp);
          entityID = ObjMgrGetEntityIDForChoice (sep);
          bsp = GetBioseqGivenSeqLoc (slp, entityID);
        }
      }
    }
  }
  return bsp;
}

static BioseqPtr GetBioseqReferencedByAnnot (SeqAnnotPtr sap, Uint2 entityID)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqFeatPtr    feat;
  SeqGraphPtr   graph;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (sap == NULL) return NULL;
  switch (sap->type) {
    case 1 :
      feat = (SeqFeatPtr) sap->data;
      while (feat != NULL) {
        slp = feat->location;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        feat = feat->next;
      }
      break;
    case 2 :
      align = (SeqAlignPtr) sap->data;
      while (align != NULL) {
        if (align->segtype == 1) {
          ddp = (DenseDiagPtr) align->segs;
          if (ddp != NULL) {
            for (sip = ddp->id; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 2) {
          dsp = (DenseSegPtr) align->segs;
          if (dsp != NULL) {
            for (sip = dsp->ids; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 3) {
          ssp = (StdSegPtr) align->segs;
          if (ssp != NULL && ssp->loc != NULL) {
            for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
              bsp = BioseqFind (SeqLocId (tloc));
              if (bsp != NULL) return bsp;
            }
          }
        }
        align = align->next;
      }
      break;
    case 3 :
      graph = (SeqGraphPtr) sap->data;
      while (graph != NULL) {
        slp = graph->loc;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        graph = graph->next;
      }
      break;
    default :
      break;
  }
  return NULL;
}

static Int4 GetScore (ScorePtr score)

{
  ObjectIdPtr  id;

  while (score != NULL) {
    id = score->id;
    if (id != NULL) {
      if (StringICmp (id->str, "score") == 0) {
        if (score->choice == 1) {
          return (score->value.intvalue);
        }
      }
    }
    score = score->next;
  }
  return 0;
}

static Int4 FindScore (SeqAlignPtr align)

{
  if (align == NULL) return 0;
  if (align->score != NULL) {
    return GetScore (align->score);
  }
  return 0;
}

static int LIBCALLBACK SortByScoreCallback (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqAlignPtr   sap1;
  SeqAlignPtr   sap2;
  Int4          score1;
  Int4          score2;

  if (ptr1 != NULL && ptr2 != NULL) {
    sap1 = *((SeqAlignPtr PNTR) ptr1);
    sap2 = *((SeqAlignPtr PNTR) ptr2);
    if (sap1 != NULL && sap2 != NULL) {
      score1 = FindScore (sap1);
      score2 = FindScore (sap2);
      if (score1 < score2) {
        return 1;
      } else if (score1 > score2) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static SeqAlignPtr SortBySeqAlignScore (SeqAlignPtr list)

{
  SeqAlignPtr  align;
  Int4         count, i;
  SeqAlignPtr  PNTR head;

  if (list == NULL) return 0;
  count = 0;
  for (align = list; align != NULL; align = align->next) {
    count++;
  }
  head = MemNew (sizeof (SeqAlignPtr) * (size_t) (count + 1));
  if (head == NULL) return 0;
  for (align = list, i = 0; align != NULL && i < count; i++) {
    head [i] = align;
    align = align->next;
  }
  HeapSort (head, (size_t) count, sizeof (SeqAlignPtr), SortByScoreCallback);
  for (i = 0; i < count; i++) {
    align = head [i];
    align->next = head [i + 1];
  }
  list = head [0];
  MemFree (head);
  return list;
}

static void TakeTop10Alignments (SeqAnnotPtr sap)

{
  SeqAlignPtr  align;
  MsgAnswer    ans;
  Int2         count;
  SeqAlignPtr  next;

  if (sap == NULL || sap->type != 2 || sap->data == NULL) return;
  count = 0;
  for (align = (SeqAlignPtr) sap->data; align != NULL; align = align->next) {
    count++;
  }
  if (count <= 10) return;
  ans = Message (MSG_YN, "Do you want to take only the top 10 (out of %d) alignments?", (int) count);
  if (ans == ANS_NO) return;
  sap->data = SortBySeqAlignScore ((SeqAlignPtr) sap->data);
  for (align = (SeqAlignPtr) sap->data, count = 0; align != NULL && count < 10; align = align->next) {
    count++;
  }
  next = align->next;
  align->next = NULL;
  align = next;
  while (align != NULL) {
    next = align->next;
    align->next = NULL;
    SeqAlignFree (align);
    align = next;
  }
}

static void DoOnePub (PubdescPtr pdp)

{
  ValNodePtr    citartptr = NULL;
  Int4          muid = 0;
  Int4          pmid = 0;
  ValNodePtr    tmp = NULL;
  ValNodePtr    vnp;

  if (pdp != NULL) {
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Muid) {
        muid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_PMid) {
        pmid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_Article) {
        citartptr = vnp;
      }
    }
    if (pmid != 0) {
      tmp = MedArchGetPubPmId (pmid);
      muid = MedArchPm2Mu (pmid);
    } else if (muid != 0) {
      tmp = MedArchGetPub (muid);
      pmid = MedArchMu2Pm (muid);
    } else if (citartptr != NULL) {
      muid = MedArchCitMatch (citartptr);
      if (muid != 0) {
        tmp = MedArchGetPub (muid);
        pmid = MedArchMu2Pm (muid);
      }
    }
    if (tmp != NULL) {
      MedlineToISO (tmp);
      if (pmid != 0) {
        ValNodeAddInt (&tmp, PUB_PMid, pmid);
      }
      if (muid != 0) {
        ValNodeAddInt (&tmp, PUB_Muid, muid);
      }
      pdp->pub = PubEquivFree (pdp->pub);
      pdp->pub = tmp;
    }
  }
}

static void DoLookupPub (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  PubdescPtr    pdp;
  SeqAnnotPtr   sap;
  ValNodePtr    sdp;
  SeqFeatPtr    sfp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          DoOnePub (pdp);
        }
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      DoOnePub (pdp);
    }
    sdp = sdp->next;
  }
}

extern void LookupAllPubs (IteM i);
extern void LookupAllPubs (IteM i)

{
  BaseFormPtr  bfp;
  MonitorPtr   mon = NULL;
  SeqEntryPtr  sep;
  ErrSev       sev;


  if (! useMedarch) return;
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
  mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
  MonitorStrValue (mon, "Connecting to MedArch");
  Update ();
  if (! MedArchInit ()) {
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
    Message (MSG_POST, "Unable to connect to MedArch");
    return;
  }
  SeqEntryExplore (sep, NULL, DoLookupPub);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  MonitorStrValue (mon, "Closing MedArch");
  Update ();
  MedArchFini ();
  MonitorFree (mon);
  ArrowCursor ();
  Update ();
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
}

static void LookupPublications (SeqAnnotPtr sap)

{
  MonitorPtr  mon = NULL;
  PubdescPtr  pdp;
  SeqFeatPtr  sfp;
  ValNodePtr  tmp;
  Int4        uid;
  Boolean     usingMedarch = FALSE;
  ValNodePtr  vnp;

  if (! useMedarch) return;
  if (sap == NULL || sap->type != 1) return;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      if (pdp != NULL) {
        vnp = pdp->pub;
        if (vnp != NULL && vnp->next == NULL) {
          if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
            if (! usingMedarch) {
              WatchCursor ();
              mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
              MonitorStrValue (mon, "Connecting to MedArch");
              Update ();
              if (MedArchInit ()) {
                usingMedarch = TRUE;
              } else {
                MonitorFree (mon);
                ArrowCursor ();
                Update ();
                Message (MSG_POST, "Unable to connect to MedArch");
                return;
              }
            }
          }
          tmp = NULL;
          if (vnp->choice == PUB_Muid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPub (uid);
          } else if (vnp->choice == PUB_PMid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPubPmId (uid);
          }
          if (tmp != NULL) {
            MedlineToISO (tmp);
            tmp->next = vnp;
            pdp->pub = tmp;
          }
        }
      }
    }
  }
  if (usingMedarch) {
    MonitorStrValue (mon, "Closing MedArch");
    Update ();
    MedArchFini ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
  }
}

static void PromotePubs (SeqFeatPtr first, BioseqPtr bsp, Uint2 entityID)

{
  MsgAnswer    ans;
  Boolean      asked = FALSE;
  PubdescPtr   pdp;
  SeqDescrPtr  sdp;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;
  ValNode      vn;

  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
  vn.next = NULL;

  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        if (! asked) {
          ans = Message (MSG_YN, "Do you wish to convert full-length publication features to descriptors?");
          if (ans == ANS_NO) return;
          asked = TRUE;
        }
      }
    }
  }

  sep = GetBestTopParentForData (entityID, bsp);
  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        sfp->idx.deleteme = TRUE;
        sfp->data.choice = SEQFEAT_COMMENT;
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        sfp->data.value.ptrvalue = NULL;
        sdp = CreateNewDescriptor (sep, Seq_descr_pub);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) pdp;
        }
      }
    }
  }

  DeleteMarkedObjects (entityID, 0, NULL);
}

extern Uint2 SmartAttachSeqAnnotToSeqEntry (Uint2 entityID, SeqAnnotPtr sap, ValNodePtr PNTR err_list)

{
  BioseqPtr      bsp;
  Int2           genCode;
  SeqEntryPtr    oldscope;
  OMProcControl  ompc;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp = NULL;

  if (sap == NULL) return entityID;
  bsp = GetBioseqReferencedByAnnot (sap, entityID);
  if (bsp == NULL) {
    oldscope = SeqEntrySetScope (NULL);
    if (oldscope != NULL) {
      bsp = GetBioseqReferencedByAnnot (sap, entityID);
      SeqEntrySetScope (oldscope);
    }
  }
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      sep = GetBestTopParentForData (entityID, bsp);
      genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
      SetEmptyGeneticCodes (sap, genCode);
      LookupPublications (sap);
    } else if (sap->type == 2) {
      TakeTop10Alignments (sap);
    }
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = entityID;
    ompc.input_itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    ompc.input_itemtype = OBJ_BIOSEQ;
    ompc.output_itemtype = OBJ_SEQANNOT;
    ompc.output_data = (Pointer) sap;
    if (! AttachDataForProc (&ompc, FALSE)) {
      if (err_list == NULL) {
        Message (MSG_ERROR, "SmartAttachSeqAnnotToSeqEntry failed");
      } else {
        ValNodeAddPointer  (err_list, 0, StringSave ("SmartAttachSeqAnnotToSeqEntry failed"));
      }
    } else if (sfp != NULL) {
      PromoteXrefs (sfp, bsp, entityID);
      PromotePubs (sfp, bsp, entityID);
    }
  } else {
    if (err_list == NULL) {
      Message (MSG_ERROR, "Feature table identifiers do not match record");
    } else {
      ValNodeAddPointer (err_list, 0, StringSave ("Feature table identifiers do not match record"));
    }
  }
  return entityID;
}

typedef struct removeformdata {
  FEATURE_FORM_BLOCK

  Boolean        is_feature;
  LisT           objlist;
  TexT           findthis;
  TexT           fromTxt;
  TexT           toTxt;
  Uint2          itemtype;
  Uint2          subtype;
  CharPtr        extra_string;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  Boolean        take_action_when_string_present;
  GrouP          string_constraint_type;
  ButtoN         case_insensitive;
  Int4           from;
  Int4           to;
  ValNodePtr     bsplist;
  ValNodePtr     bssplist;
} RemoveFormData, PNTR RemoveFormPtr;

static Boolean ObjectInRange (SeqFeatPtr sfp, Int4 from, Int4 to)

{
  SeqMgrFeatContext  context;

  if (sfp == NULL || from < 0 || to < 0) return TRUE;
  if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) == sfp) {
    if (context.left > to) return FALSE;
    if (context.right < from) return FALSE;
  }
  return TRUE;
}


static void RemoveFeatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  RemoveFormPtr rfp;
  SeqIdPtr      sip;
  BioseqPtr     productbsp, productcdna;
  BioseqSetPtr  productnps;
  
  if (sfp == NULL || userdata == NULL) return;

  rfp = (RemoveFormPtr) userdata;
  if (rfp == NULL) return;
  if (sfp->idx.subtype == rfp->subtype ||
      (rfp->subtype == FEATDEF_IMP && IsRealImpFeat (sfp->idx.subtype)) ||
      rfp->subtype == ALL_FEATURES) 
  {
    if ((rfp->from == -1 && rfp->to == -1) || ObjectInRange (sfp, rfp->from, rfp->to)) 
    {
      if (sfp->data.choice == SEQFEAT_CDREGION) 
      {
        if (sfp->product != NULL) 
        {
          sip = SeqLocId (sfp->product);
          if (sip != NULL) 
          {
            productbsp = BioseqFind (sip);
            if (productbsp != NULL) 
            {
              ValNodeAddPointer (&(rfp->bsplist), 0, (Pointer) productbsp);
            }
          }
        }
      } 
      else if (sfp->data.choice == SEQFEAT_RNA) 
      {
        if (sfp->product != NULL) 
        {
          sip = SeqLocId (sfp->product);
          if (sip != NULL) 
          {
            productcdna = BioseqFind (sip);
            if (productcdna != NULL && productcdna->idx.parenttype == OBJ_BIOSEQSET) 
            {
              productnps = (BioseqSetPtr) productcdna->idx.parentptr;
              if (productnps != NULL && productnps->_class == BioseqseqSet_class_nuc_prot) 
              {
                ValNodeAddPointer (&(rfp->bssplist), 0, (Pointer) productnps);
              }
            }
          }
        }
      }
      sfp->idx.deleteme = TRUE;  
    }
  }
}

static void RemoveFeatures (SeqEntryPtr sep, RemoveFormPtr rfp)
{
  FeaturesWithTextData fd;
  Char           str [32];
  Int4           swap;
  long int       val;
  
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  fd.search_text = rfp->findStr;
  fd.no_text = StringHasNoText (rfp->findStr);
  fd.seqFeatChoice = 0;
  fd.featDefChoice = 0;
  fd.case_insensitive = GetStatus (rfp->case_insensitive);
  fd.whole_word = FALSE;
  fd.act_when_string_not_present = ! rfp->take_action_when_string_present;
  fd.userdata = rfp;
  fd.callback = RemoveFeatureCallback;
  GetTitle (rfp->fromTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->from = (Int4) val;
  } else {
    rfp->from = -1;
  }
  GetTitle (rfp->toTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->to = (Int4) val;
  } else {
    rfp->to = -1;
  }
  if (rfp->from > rfp->to) {
    swap = rfp->from;
    rfp->from = rfp->to;
    rfp->to = swap;
  }
  OperateOnSeqEntryFeaturesWithText (sep, &fd);
  DeleteMarkedObjects (rfp->input_entityID, OBJ_SEQENTRY, (Pointer) sep);
}


static Boolean IsUserObjectType (SeqDescrPtr sdp, CharPtr string)
{
  UserObjectPtr uop;
  ObjectIdPtr   oip;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return FALSE;
  if (StringHasNoText (string)) return TRUE;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return FALSE;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, string) != 0) return FALSE;
  return TRUE;
}


static void RemoveDescriptorCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovp;
  RemoveFormPtr rfp;
  
  if (sdp == NULL || userdata == NULL || sdp->extended == 0) return;
  rfp = (RemoveFormPtr) userdata;
  if (rfp == NULL) return;
  
  ovp = (ObjValNodePtr) sdp;
    
  if (sdp->choice == rfp->subtype) 
  {
    if (sdp->choice != Seq_descr_user || IsUserObjectType (sdp, rfp->extra_string)) {
      ovp->idx.deleteme = TRUE;	
    }
  }
}

static void RemoveDescriptors (SeqEntryPtr sep, RemoveFormPtr rfp)
{
  DescriptorsWithTextData dd;
  Char                    str [32];
  Int4                    swap;
  long int                val;
  
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  dd.search_text = rfp->findStr;
  dd.no_text = StringHasNoText (rfp->findStr);
  dd.case_insensitive = GetStatus (rfp->case_insensitive);
  dd.whole_word = FALSE;
  dd.act_when_string_not_present = ! rfp->take_action_when_string_present;
  dd.userdata = rfp;
  dd.callback = RemoveDescriptorCallback;
  GetTitle (rfp->fromTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->from = (Int4) val;
  } else {
    rfp->from = -1;
  }
  GetTitle (rfp->toTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->to = (Int4) val;
  } else {
    rfp->to = -1;
  }
  if (rfp->from > rfp->to) {
    swap = rfp->from;
    rfp->from = rfp->to;
    rfp->to = swap;
  }
  OperateOnSeqEntryDescriptorsWithText (sep, &dd);
  DeleteMarkedObjects (rfp->input_entityID, OBJ_SEQENTRY, (Pointer) sep);
}

static void DoRemoveAsnObject (ButtoN b)

{
  MsgAnswer      ans;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Uint4          itemID;
  OMProcControl  ompc;
  RemoveFormPtr  rfp;
  SeqEntryPtr    sep;
  ValNodePtr     tmp;
  Int2           val;
  ValNodePtr     vnp;
  Boolean        removed_some_features;

  rfp = GetObjectExtra (b);
  if (rfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (rfp->input_entityID);
  if (sep == NULL) return;
  Hide (rfp->form);
  WatchCursor ();
  Update ();
  if (rfp->is_feature) {
    rfp->itemtype = OBJ_SEQFEAT;
  } else {
    rfp->itemtype = OBJ_SEQDESC;
  }

  if (rfp->itemtype == 0) return;

  removed_some_features = FALSE;
  
  if (GetValue (rfp->string_constraint_type) == 1)
  {
  	rfp->take_action_when_string_present = TRUE;
  }
  else
  {
    rfp->take_action_when_string_present = FALSE;
  }

  val = 1;
  for (vnp = rfp->head; vnp != NULL; vnp = vnp->next)
  {
    if (GetItemStatus (rfp->objlist, val))
    {
      rfp->subtype = vnp->choice;
      rfp->extra_string = NULL;
      if (rfp->subtype != 0) {
        if (rfp->is_feature) {
          RemoveFeatures (sep, rfp);
          removed_some_features = TRUE;
        } else {
          if (rfp->subtype == Seq_descr_user && StringCmp (vnp->data.ptrvalue, "User") != 0) {
            rfp->extra_string = vnp->data.ptrvalue;
          }
          RemoveDescriptors (sep, rfp);
        }
      }
    }
    val ++;
  }

  if (removed_some_features) {
    if (rfp->bsplist != NULL) {
      ans = Message (MSG_YN, "Remove protein products?");
      if (ans == ANS_YES) {
        for (tmp = rfp->bsplist; tmp != NULL; tmp = tmp->next) {
          bsp = (BioseqPtr) tmp->data.ptrvalue;
          itemID = GetItemIDGivenPointer (rfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
          if (itemID > 0) {
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = rfp->input_entityID;
            ompc.input_itemID = itemID;
            ompc.input_itemtype = OBJ_BIOSEQ;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_POSTERR, "DetachDataForProc failed");
            }
            SeqMgrDeleteFromBioseqIndex (bsp);
          }
        }
        ans = Message (MSG_YN, "Renormalize Nuc-Prot sets?");
        if (ans == ANS_YES)
        {
          RemoveOrphanProteins (rfp->input_entityID, sep);
          RenormalizeNucProtSets (sep, TRUE);   	
        }
      }
    }
    if (rfp->bssplist != NULL) {
      ans = Message (MSG_YN, "Remove cDNA nuc-prot products?");
      if (ans == ANS_YES) {
        for (tmp = rfp->bssplist; tmp != NULL; tmp = tmp->next) {
          bssp = (BioseqSetPtr) tmp->data.ptrvalue;
          itemID = GetItemIDGivenPointer (rfp->input_entityID, OBJ_BIOSEQSET, (Pointer) bssp);
          if (itemID > 0) {
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = rfp->input_entityID;
            ompc.input_itemID = itemID;
            ompc.input_itemtype = OBJ_BIOSEQSET;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_POSTERR, "DetachDataForProc failed");
            }
          }
        }
      }
    }
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (rfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (rfp->form);
}

static void RemoveDefLinesCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
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
    if (sdp->choice == Seq_descr_title)
    {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

extern void RemoveDefLinesToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  WatchCursor ();
  Update ();

  SeqEntryExplore (sep, NULL, RemoveDefLinesCallback);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  CommonApplyToAllProc (bfp, ADD_TITLE);
}

int LIBCALLBACK SortByVnpChoice (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice > vnp2->choice) {
        return 1;
      } else if (vnp1->choice < vnp2->choice) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static void RemoveMessageProc (ForM f, Int2 mssg)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    if (rfp->appmessage != NULL) {
      rfp->appmessage (f, mssg);
    }
  }
}

static void CleanupRemovePage (GraphiC g, VoidPtr data)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) data;
  if (rfp != NULL) {
    ValNodeFreeData (rfp->head);
    ValNodeFree (rfp->bsplist);
    ValNodeFree (rfp->bssplist);
  }
  StdCleanupFormProc (g, data);
}

static CharPtr descNames [] = {
  " ", " ", " ", " ", "Name",
  "Title", " ", "Comment", "Numbering",
  "MapLoc", "PIR", "GenBank", "Publication",
  "Region", "User", "SWISS-PROT", "dbXREF",
  "EMBL", "Create Date", "Update Date", "PRF",
  "PDB", "Heterogen", "BioSource", "MolInfo", NULL
};

/*
#ifdef INTERNAL_NCBI_SEQUIN
#define LISTHEIGHT 16
#else
#define LISTHEIGHT 8
#endif
*/

CharPtr MostUsedDescriptorList[] = { "Title" };

static Boolean isMostUsedDescriptor (CharPtr descname)
{
  Int2 i;

  if (descname == NULL) return FALSE;

  for (i=0; i < sizeof (MostUsedDescriptorList) / sizeof (CharPtr); i++)
  {
    if (StringCmp (descname, MostUsedDescriptorList[i]) == 0)
      return TRUE;
  }
  return FALSE;
}

static int LIBCALLBACK SortMostUsedDescriptorsFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;
  CharPtr      str1;
  CharPtr      str2;
  Boolean      str1_is_most_used;
  Boolean      str2_is_most_used;

  /* Check parameters */

  if ((NULL == ptr1) || (NULL == ptr2))
    return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if ((NULL == vnp1) || (NULL == vnp2))
    return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if ((NULL == str1) || (NULL == str2))
    return 0;

  str1_is_most_used = isMostUsedDescriptor (str1);
  str2_is_most_used = isMostUsedDescriptor (str2);

  if ((str1_is_most_used && str2_is_most_used)
    || (!str1_is_most_used && !str2_is_most_used))
  {
    return SortVnpByString (ptr1, ptr2);
  }
  else if (str1_is_most_used)
  {
    return -1;
  }
  else
  {
    return 1;
  }
}

extern ValNodePtr BuildDescriptorValNodeList (void)
{
  Int4 j;
  ValNodePtr vnp;
  ValNodePtr head = NULL;
  
  for (j = 1; descNames [j] != NULL; j++) {
    if (StringHasNoText (descNames [j])) continue;
    vnp = ValNodeNew (head);
    if (head == NULL) {
      head = vnp;
    }
    if (vnp != NULL) {
      vnp->choice = j;
      vnp->data.ptrvalue = StringSave (descNames [j]);
    }
  }
  head = SortValNode (head, SortMostUsedDescriptorsFirst);
  return head;
}


static void RemoveAsnObject (IteM i, Boolean feature)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  GrouP              k;
  Int2               listHeight;
  GrouP              m;
  RemoveFormPtr      rfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  rfp = (RemoveFormPtr) MemNew (sizeof (RemoveFormData));
  if (rfp == NULL) return;
  if (feature) {
    title = "Feature Removal";
  } else {
    title = "Descriptor Removal";
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, rfp, CleanupRemovePage);
  rfp->form = (ForM) w;
  rfp->formmessage = RemoveMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    rfp->appmessage = sepp->handleMessages;
  }

  rfp->input_entityID = bfp->input_entityID;
  rfp->input_itemID = bfp->input_itemID;
  rfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  rfp->is_feature = feature;
  if (feature) {
    StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  } else {
    StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  rfp->objlist = MultiList (g, 16, listHeight, NULL);
  head = NULL;
  if (feature) {
    head = BuildFeatureValNodeList (TRUE, "All", ALL_FEATURES, TRUE, FALSE);
  } else {
    head = BuildDescriptorValNodeList();
    vnp = ValNodeNew (NULL);
    vnp->choice = Seq_descr_user;
    vnp->data.ptrvalue = StringSave ("StructuredComment");
    ValNodeInsert (&(head->next), vnp, SortVnpByString);
  }
  if (head != NULL) {

    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (rfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  rfp->head = head;
  rfp->bsplist = NULL;
  rfp->bssplist = NULL;

  k = NormalGroup (h, 0, 3, "Optional string constraint", NULL, NULL);
  rfp->string_constraint_type = HiddenGroup (k, 0, 2, NULL);
  RadioButton (rfp->string_constraint_type, "Remove when text is present");
  RadioButton (rfp->string_constraint_type, "Remove when text is not present");
  SetValue (rfp->string_constraint_type, 1);
  rfp->findthis = DialogText (k, "", 14, NULL);
  rfp->case_insensitive = CheckBox (k, "Case Insensitive", NULL);

  m = NULL;
  if (feature) {
    m = HiddenGroup (h, 4, 0, NULL);
    StaticPrompt (m, "From", 0, dialogTextHeight, programFont, 'l');
    rfp->fromTxt = DialogText (m, "", 6, NULL);
    StaticPrompt (m, "To", 0, dialogTextHeight, programFont, 'l');
    rfp->toTxt = DialogText (m, "", 6, NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveAsnObject);
  SetObjectExtra (b, rfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) c, (HANDLE) m, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void RemoveDescriptor (IteM i)

{
  RemoveAsnObject (i, FALSE);
}

#define SLCT_FEAT    1
#define SLCT_DESC    2
#define SLCT_BIOSEQ  3
#define SLCT_PUB     4

typedef struct selectformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  ButtoN         when_string_not_present;
  ButtoN         case_insensitive;
} SelectFormData, PNTR SelectFormPtr;

static void FeatureSelectCallback (SeqFeatPtr sfp, Pointer userdata)
{
  Uint1Ptr subtype;
  if (sfp == NULL) return;
  
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != sfp->idx.subtype) return;
  }
  ObjMgrAlsoSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
}

static void DescriptorSelectCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovp;
  Uint1Ptr subtype;
  
  if (sdp == NULL || sdp->extended == 0) return;
  
  ovp = (ObjValNodePtr) sdp;
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != ovp->idx.subtype) return;
  }
    
  ObjMgrAlsoSelect (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, 0, NULL);
}

static void BioseqSelectCallback (BioseqPtr bsp, Pointer userdata)
{
  Uint1Ptr subtype;
  
  if (bsp == NULL) return;
  
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != bsp->idx.subtype) return;
  }
    
  ObjMgrAlsoSelect (bsp->idx.entityID, bsp->idx.itemID, bsp->idx.itemtype, 0, NULL);
}


static void DoSelectAsnObject (ButtoN b)

{
  SelectFormPtr           selfp;
  SeqEntryPtr             sep;
  Int2                    val;
  ValNodePtr              vnp;
  FeaturesWithTextData    fd;
  DescriptorsWithTextData dd;
  Uint1                   bioseq_choice = Seq_repr_raw;
  Uint1                   pub_choice = FEATDEF_PUB;

  selfp = GetObjectExtra (b);
  if (selfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (selfp->input_entityID);
  if (sep == NULL) return;
  Hide (selfp->form);
  
  vnp = NULL;
  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC)
  {
    val = GetValue (selfp->objlist);
    if (val > 0) {
      vnp = selfp->head;
      while (vnp != NULL && val > 1) {
        val--;
        vnp = vnp->next;
      }
    }
  }

  switch (selfp->type) {
    case SLCT_FEAT :
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      fd.search_text = selfp->findStr;
      fd.no_text = StringHasNoText (selfp->findStr);
      fd.seqFeatChoice = 0;
      fd.featDefChoice = 0;
      fd.case_insensitive = GetStatus (selfp->case_insensitive);
      fd.whole_word = FALSE;
      fd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      fd.callback = FeatureSelectCallback;
      if (vnp == NULL)
      {
      	fd.userdata = NULL;
      }
      else
      {
        fd.userdata = (Pointer) &(vnp->choice);
      }
      OperateOnSeqEntryFeaturesWithText (sep, &fd); 	
      break;
    case SLCT_DESC :
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      dd.search_text = selfp->findStr;
      dd.no_text = StringHasNoText (selfp->findStr);
      dd.case_insensitive = GetStatus (selfp->case_insensitive);
      dd.whole_word = FALSE;
      dd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      dd.callback = DescriptorSelectCallback;
      if (vnp == NULL)
      {
      	dd.userdata = NULL;
      }
      else
      {
        dd.userdata = (Pointer) &(vnp->choice);
      }
      OperateOnSeqEntryDescriptorsWithText (sep, &dd); 	
      break;
    case SLCT_BIOSEQ :
 	  VisitBioseqsInSep (sep, (Pointer) &bioseq_choice, BioseqSelectCallback);	
      break;
  	case SLCT_PUB:
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      fd.search_text = selfp->findStr;
      fd.no_text = StringHasNoText (selfp->findStr);
      fd.seqFeatChoice = 0;
      fd.featDefChoice = 0;
      fd.case_insensitive = GetStatus (selfp->case_insensitive);
      fd.whole_word = FALSE;
      fd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      fd.callback = FeatureSelectCallback;
      fd.userdata = (Pointer) &pub_choice;
      OperateOnSeqEntryFeaturesWithText (sep, &fd); 
      dd.search_text = fd.search_text;
      dd.no_text = fd.no_text;
      dd.case_insensitive = fd.case_insensitive;
      dd.whole_word = fd.whole_word;
      dd.act_when_string_not_present = fd.act_when_string_not_present;
      dd.callback = DescriptorSelectCallback;
      dd.userdata = fd.userdata;	
      OperateOnSeqEntryDescriptorsWithText (sep, &dd); 	
      break;  	
    default :
      Remove (selfp->form);
      Update ();
      return;
  }
  WatchCursor ();
  Update ();

  ArrowCursor ();
  Update ();
  /* ObjMgrSendMsg (OM_MSG_UPDATE, selfp->input_entityID, 0, 0); */
  Remove (selfp->form);
}

static void SelectMessageProc (ForM f, Int2 mssg)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) GetObjectExtra (f);
  if (selfp != NULL) {
    if (selfp->appmessage != NULL) {
      selfp->appmessage (f, mssg);
    }
  }
}

static void CleanupSelectPage (GraphiC g, VoidPtr data)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) data;
  if (selfp != NULL) {
    ValNodeFreeData (selfp->head);
  }
  StdCleanupFormProc (g, data);
}

static void SelectAsnObject (IteM i, Int2 type)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  GrouP              k, m;
  ValNodePtr         head;
  Uint1              j;
  Int2               listHeight;
  SelectFormPtr      selfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  selfp = (SelectFormPtr) MemNew (sizeof (SelectFormData));
  if (selfp == NULL) return;
  switch (type) {
    case SLCT_FEAT :
      title = "Feature Selection";
      break;
    case SLCT_DESC :
      title = "Descriptor Selection";
      break;
    case SLCT_BIOSEQ :
      title = "Sequence Selection";
      break;
  	case SLCT_PUB:
  	  title = "Publication Selection";
  	  break;
    default :
      title = "? Selection";
      break;
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, selfp, CleanupSelectPage);
  selfp->form = (ForM) w;
  selfp->formmessage = SelectMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    selfp->appmessage = sepp->handleMessages;
  }

  selfp->input_entityID = bfp->input_entityID;
  selfp->input_itemID = bfp->input_itemID;
  selfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  selfp->type = type;
  switch (type) {
    case SLCT_FEAT :
      StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
      break;
    case SLCT_DESC :
      StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
      break;
    case SLCT_BIOSEQ :
      StaticPrompt (g, "Sequence", 0, 0, programFont, 'c');
      break;
  	case SLCT_PUB:
  	  StaticPrompt (g, "Publication", 0, 0, programFont, 'c');
  	  break;
    default :
      break;
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  
  if (type != SLCT_PUB)
  {
    selfp->objlist = SingleList (g, 16, listHeight, NULL);
  }
  head = NULL;
  if (type == SLCT_FEAT) {
    head = BuildFeatureValNodeList (TRUE, NULL, ALL_FEATURES, TRUE, FALSE);
  } else if (type == SLCT_DESC) {
    for (j = 1; descNames [j] != NULL; j++) {
      if (StringHasNoText (descNames [j])) continue;
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = j;
        vnp->data.ptrvalue = StringSave (descNames [j]);
      }
    }
  }
  if (head != NULL) {
    if (type != SLCT_FEAT) {
      head = SortValNode (head, SortByVnpChoice);
    }
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (selfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  selfp->head = head;

  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC || selfp->type == SLCT_PUB)
  {
    k = HiddenGroup (h, 0, 3, NULL);
    StaticPrompt (k, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
    selfp->findthis = DialogText (k, "", 14, NULL);
    m = HiddenGroup (k, 2, 0, NULL);
    selfp->case_insensitive = CheckBox (m, "Case Insensitive", NULL);
    selfp->when_string_not_present = CheckBox (m, "When String Not Present", NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoSelectAsnObject);
  SetObjectExtra (b, selfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) c, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  }
  RealizeWindow (w);
  if (type == SLCT_BIOSEQ) {
    DoSelectAsnObject (b);
    Update ();
    return;
  }
  Show (w);
  Update ();
}

extern void SelectDescriptor (IteM i)

{
  SelectAsnObject (i, SLCT_DESC);
}

extern void SelectBioseq (IteM i)

{
  SelectAsnObject (i, SLCT_BIOSEQ);
}

extern void SelectPubs (IteM i)

{
  SelectAsnObject (i, SLCT_PUB);
}

typedef struct fuseformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  Uint2          subtype;
  ValNodePtr     head;
} FuseFormData, PNTR FuseFormPtr;

static SeqLocPtr FuseTwoLocations (Uint2 entityID, SeqLocPtr slp1, SeqLocPtr slp2)
{
  Boolean           slp1_partial5, slp1_partial3;
  Boolean           slp2_partial5, slp2_partial3;
  Int4              slp1_start_pos, slp1_stop_pos;
  Int4              slp2_start_pos, slp2_stop_pos;
  SeqLocPtr         slp_result = NULL, slp2_copy, slp1_copy, slp_list;
  BioseqPtr         bsp1, bsp2;

  if (slp1 == NULL || slp2 == NULL)
  {
    return NULL;
  }
  
  bsp1 = GetBioseqGivenSeqLoc (slp1, entityID);
  bsp2 = GetBioseqGivenSeqLoc (slp2, entityID);
  
  /* preserve partialness of ends */
  CheckSeqLocForPartial (slp1, &slp1_partial5, &slp1_partial3);
  slp1_start_pos = SeqLocStart (slp1);
  slp1_stop_pos = SeqLocStop (slp1);
  CheckSeqLocForPartial (slp2, &slp2_partial5, &slp2_partial3);
  slp2_start_pos = SeqLocStart (slp2);
  slp2_stop_pos = SeqLocStop (slp2);
  if (slp1_start_pos > slp2_start_pos)
  {
    slp1_partial5 = slp2_partial5;
  }
  if (slp1_stop_pos < slp2_stop_pos)
  {
    slp1_partial3 = slp2_partial3;
  }

  if (bsp1 == bsp2)
  {
    slp_result = SeqLocMerge (bsp1, slp2, slp1, FALSE, TRUE, FALSE);
    SetSeqLocPartial (slp_result, slp1_partial5, slp1_partial3);
  }
  else 
  {
    /* we are dealing with a segmented set, and both locations are segments */
    slp1_copy = (SeqLocPtr) AsnIoMemCopy (slp1, (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
    slp2_copy = (SeqLocPtr) AsnIoMemCopy (slp2, (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
    if (slp1_copy != NULL && slp2_copy != NULL)
    {
      /* if the second location is already a mix, don't create a nested
       * mixed location.
       */
      if (slp2_copy->choice == SEQLOC_MIX)
      {
        slp_list = (SeqLocPtr) slp2_copy->data.ptrvalue;
        slp2_copy->data.ptrvalue = NULL;
        slp2_copy = SeqLocFree (slp2_copy);
        slp2_copy = slp_list;
      }
    
      if (slp1_copy->choice == SEQLOC_MIX)
      {
        slp_list = slp1_copy->data.ptrvalue;
        while (slp_list != NULL && slp_list->next != NULL)
        {
          slp_list = slp_list->next;
        }
        
        
        if (slp_list == NULL)
        {
          slp1_copy->data.ptrvalue = slp2_copy;
        }
        else
        {
          slp_list->next = slp2_copy;
        }
        slp_result = slp1_copy;
      }
      else
      {
        slp_result = ValNodeNew (NULL);
        if (slp_result != NULL)
        {
          slp_result->choice = SEQLOC_MIX;
          slp_result->data.ptrvalue = slp1_copy;
          slp1_copy->next = slp2_copy;
        }
      }
    }
  }
  
  return slp_result;
}

static void CombineProductFeatures (BioseqPtr pbsp1, BioseqPtr pbsp2, Int4 old_len)
{
  SeqAnnotPtr sap1, sap2, last_sap1 = NULL;
  SeqFeatPtr  sfp2, sfp_new, last_sfp1 = NULL, main_prot = NULL;
  Boolean     partial5_orig = TRUE, partial3_orig = TRUE;
  Boolean     partial5_new, partial3_new;
  SeqEntryPtr psep;
  ProtRefPtr  prp;
  
  if (pbsp1 == NULL || pbsp2 == NULL)
  {
    return;
  }
  
  sap2 = pbsp2->annot;
  while (sap2 != NULL && sap2->type != 1)
  {
    sap2 = sap2->next;
  }
  if (sap2 == NULL || sap2->data == NULL)
  {
    /* second sequence has no features */
    return;
  }
  
  sap1 = pbsp1->annot;
  while (sap1 != NULL && sap1->type != 1)
  {
    last_sap1 = sap1;
    sap1 = sap1->next;
  }
  if (sap1 == NULL)
  {
    sap1 = SeqAnnotNew();
    if (sap1 == NULL)
    {
      return;
    }
    sap1->type = 1;
    sap1->data = NULL;
    if (last_sap1 == NULL)
    {
      pbsp1->annot = sap1;
    }
    else
    {
      last_sap1->next = sap1;
    }
  }
  
  last_sfp1 = sap1->data;
  if (last_sfp1 != NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
  {
    main_prot = last_sfp1;
    CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
  }
  while (last_sfp1 != NULL && last_sfp1->next != NULL)
  {
    if (main_prot == NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
    {
      main_prot = last_sfp1;
      CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
    }
    last_sfp1 = last_sfp1->next;
  }
  if (last_sfp1 != NULL && main_prot == NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
  {
    main_prot = last_sfp1;
    CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
  }
  
  partial3_new = partial3_orig;
  
  while (sap2 != NULL)
  {
    if (sap2->type == 1)
    {
      for (sfp2 = sap2->data; sfp2 != NULL; sfp2 = sfp2->next)
      {  
        if (sfp2->idx.subtype == FEATDEF_PROT)
        {
          CheckSeqLocForPartial (sfp2->location, &partial5_new, &partial3_new);
          /* do not create additional full-length protein features */
          continue;
        }
        sfp_new = (SeqFeatPtr) AsnIoMemCopy (sfp2, (AsnReadFunc) SeqFeatAsnRead,
                                               (AsnWriteFunc) SeqFeatAsnWrite);
        if (sfp_new != NULL)
        {
          OffsetLocation (sfp_new->location, old_len, pbsp1->id);

          if (last_sfp1 == NULL)
          {
            sap1->data = sfp_new;
          }
          else
          {
            last_sfp1->next = sfp_new;
          }
          last_sfp1 = sfp_new;
        }
      }
    }
    sap2 = sap2->next;
  }

  /* make sure there is one full-length protein feature */
  if (main_prot == NULL)
  {
    psep = SeqMgrGetSeqEntryForData (pbsp1);
    main_prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (main_prot != NULL) 
    {
      prp = ProtRefNew ();
      main_prot->data.value.ptrvalue = (Pointer) prp;
    }
  }
  if (main_prot != NULL)
  {
    if (main_prot->location == NULL || main_prot->location->choice != SEQLOC_INT)
    {
      main_prot->location = SeqLocFree (main_prot->location);
      main_prot->location = SeqLocIntNew (0, pbsp1->length - 1, 
                                          Seq_strand_plus,
                                          SeqIdDup (pbsp1->id));
    }
    SetSeqLocPartial (main_prot->location, partial5_orig, partial3_new);
    main_prot->partial = (Boolean) (partial5_orig || partial3_new);
  }
}

static void FuseTwoProducts (SeqFeatPtr sfp1, SeqFeatPtr sfp2, Uint2 entityID)
{
  BioseqPtr    pbsp1, pbsp2;
  CharPtr      pstr1, pstr2;
  ByteStorePtr byte_store, bs2 = NULL;
  Int4         old_length, added_length;
  SeqIdPtr     sip1;
  
  if (sfp1 == NULL || sfp2 == NULL 
      || sfp1->idx.subtype != FEATDEF_CDS
      || sfp2->idx.subtype != FEATDEF_CDS)
  {
    return;
  }
  pbsp1 = BioseqFindFromSeqLoc (sfp1->product);
  pbsp2 = BioseqFindFromSeqLoc (sfp2->product);
  
  if (pbsp1 == NULL)
  {
    sip1 = SeqLocId (sfp1->product);
    pbsp1 = BioseqFind (sip1);
    if (pbsp1 == NULL)
    {
      RetranslateOneCDS (sfp1, entityID, TRUE, FALSE);
      pbsp1 = BioseqFindFromSeqLoc (sfp1->product);
      if (pbsp1 == NULL)
      {
        sip1 = SeqLocId (sfp1->product);
        pbsp1 = BioseqFind (sip1);
      }
    }
  }
  if (pbsp1 == NULL)
  {
    return;
  }
  pstr1 = BSMerge ((ByteStorePtr)(pbsp1->seq_data), NULL);
  old_length = pbsp1->length;
  
  if (pbsp2 == NULL)
  {
    bs2 = ProteinFromCdRegionEx (sfp2, TRUE, FALSE);
    pstr2 = BSMerge (bs2, NULL);
    added_length = BSLen (bs2);
  }
  else
  {
    pstr2 = BSMerge ((ByteStorePtr)(pbsp2->seq_data), NULL);
    added_length = pbsp2->length;
  }
  
  byte_store = BSNew (old_length + added_length);
  if (byte_store != NULL)
  {
    BSWrite (byte_store, pstr1, StringLen (pstr1));
    BSWrite (byte_store, pstr2, StringLen (pstr2));
    pbsp1->seq_data = SeqDataFree (pbsp1->seq_data, pbsp1->seq_data_type);
    pbsp1->seq_data = (SeqDataPtr) byte_store;
    pbsp1->length += added_length;
      
    /* now copy features from the second protein to the first */
    CombineProductFeatures (pbsp1, pbsp2, old_length);     
      
    /* remove unused protein */
    if (pbsp2 != NULL)
    {
      pbsp2->idx.deleteme = TRUE;
    }
  }
  bs2 = BSFree (bs2);
}

static void FuseFeatureCallback (BioseqPtr bsp, Pointer userdata)
{
  FuseFormPtr       ffp;
  SeqFeatPtr        first = NULL, sfp = NULL;
  SeqMgrFeatContext context;
  SeqLocPtr         slp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  ffp = (FuseFormPtr) userdata;
  
  sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  while (sfp != NULL)
  {
    if (sfp->idx.subtype == ffp->subtype ||
           (ffp->subtype == FEATDEF_IMP &&
            IsRealImpFeat (sfp->idx.subtype))) 
    {
      if (first == NULL)
      {
        first = sfp;
      }
      else
      {
        slp = FuseTwoLocations (ffp->input_entityID, first->location, sfp->location);
        first->location = SeqLocFree (first->location);
        first->location = slp;
        first->partial = CheckSeqLocForPartial (slp, NULL, NULL);
        sfp->idx.deleteme = TRUE;
        if (sfp->idx.subtype == FEATDEF_CDS)
        {
          FuseTwoProducts (first, sfp, ffp->input_entityID);
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }  
}


static void DoFuseFeature (ButtoN b)

{
  FuseFormPtr  ffp;
  SeqEntryPtr  sep;
  Int2         val;
  ValNodePtr   vnp;

  ffp = (FuseFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
  if (sep == NULL) return;
  Hide (ffp->form);
  WatchCursor ();
  Update ();

  vnp = NULL;
  val = GetValue (ffp->objlist);
  if (val > 0) {
    vnp = ffp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    ffp->subtype = vnp->choice;
    VisitBioseqsInSep (sep, ffp, FuseFeatureCallback);
    DeleteMarkedObjects (ffp->input_entityID, 0, NULL);
  }

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (ffp->form);
}

static void FuseMessageProc (ForM f, Int2 mssg)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    if (ffp->appmessage != NULL) {
      ffp->appmessage (f, mssg);
    }
  }
}

static void CleanupFusePage (GraphiC g, VoidPtr data)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) data;
  if (ffp != NULL) {
    ValNodeFreeData (ffp->head);
  }
  StdCleanupFormProc (g, data);
}

extern void FuseFeature (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  FuseFormPtr        ffp;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint1              subtype;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ffp = (FuseFormPtr) MemNew (sizeof (FuseFormData));
  if (ffp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Fuse Feature", StdCloseWindowProc);
  SetObjectExtra (w, ffp, CleanupFusePage);
  ffp->form = (ForM) w;
  ffp->formmessage = FuseMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ffp->appmessage = sepp->handleMessages;
  }

  ffp->input_entityID = bfp->input_entityID;
  ffp->input_itemID = bfp->input_itemID;
  ffp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  ffp->objlist = SingleList (g, 16, listHeight, NULL);
  head = NULL;
  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD) {
      subtype = curr->featdef_key;
      if (subtype != FEATDEF_misc_RNA &&
          subtype != FEATDEF_precursor_RNA &&
          subtype != FEATDEF_mat_peptide &&
          subtype != FEATDEF_sig_peptide &&
          subtype != FEATDEF_transit_peptide &&
          subtype != FEATDEF_Imp_CDS &&
          !IsUnwantedFeatureType(subtype)) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (curr->typelabel);
        }
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  if (head != NULL) {
    head = SortValNode (head, SortByVnpChoice);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (ffp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  ffp->head = head;

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoFuseFeature);
  SetObjectExtra (b, ffp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}





typedef struct genomeprojiduserdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        ids;
} GenomeprojidUserDialog, PNTR GenomeprojidUserDialogPtr;

typedef struct genomeprojiduserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} GenomeprojidUserForm, PNTR GenomeprojidUserFormPtr;

static void UserObjectPtrToGenomeprojidDialog (DialoG d, Pointer data)

{
  Char                       buf [64];
  UserFieldPtr               curr;
  ValNodePtr                 head = NULL;
  ObjectIdPtr                oip;
  Int4                       parentID;
  Int4                       projectID;
  GenomeprojidUserDialogPtr  rdp;
  UserObjectPtr              uop;
  Int4                       val;

  rdp = (GenomeprojidUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return;
  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "GenomeProjectsDB") != 0) {
    PointerToDialog (rdp->ids, NULL);
    return;
  }
  projectID = 0;
  parentID = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip == NULL) continue;
    if (StringICmp (oip->str, "ProjectID") == 0) {
      if (curr->choice == 2) {
        val = (Int4) curr->data.intvalue;
        if (projectID > 0) {
          sprintf (buf, "%ld\t%ld", (long) projectID, (long) parentID);
          ValNodeCopyStr (&head, 0, buf);
          parentID = 0;
        }
        projectID = val;
      }
    } else if (StringICmp (oip->str, "ParentID") == 0) {
      if (curr->choice == 2) {
        val = (Int4) curr->data.intvalue;
        parentID = val;
      }
    }
  }
  if (projectID > 0) {
    sprintf (buf, "%ld\t%ld", (long) projectID, (long) parentID);
    ValNodeCopyStr (&head, 0, buf);
  }

  PointerToDialog (rdp->ids, (Pointer) head);
  ValNodeFreeData (head);
}

static Pointer GenomeprojidDialogToUserObjectPtr (DialoG d)

{
  Char                       buf [64];
  ValNodePtr                 head;
  Int4                       parentID;
  Int4                       projectID;
  CharPtr                    ptr1;
  CharPtr                    ptr2;
  GenomeprojidUserDialogPtr  rdp;
  CharPtr                    str;
  UserObjectPtr              uop;
  long int                   val;
  ValNodePtr                 vnp;

  rdp = (GenomeprojidUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return NULL;

  uop = CreateGenomeProjectsDBUserObject ();
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (rdp->ids);
  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    projectID = 0;
    parentID = 0;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    StringNCpy_0 (buf, str, sizeof (buf));
    ptr1 = StringChr (buf, '\t');
    if (ptr1 != NULL) {
      *ptr1 = '\0';
      ptr1++;
      ptr2 = StringChr (ptr1, '\n');
      if (ptr2 == NULL) {
        ptr2 = StringChr (ptr1, '\t');
      }
      if (ptr2 != NULL) {
        *ptr2 = '\0';
      }
      if (sscanf (buf, "%ld", &val) == 1 && val > 0) {
        projectID = (Int4) val;
        if (sscanf (ptr1, "%ld", &val) == 1 && val > 0) {
          parentID = (Int4) val;
        }
        AddIDsToGenomeProjectsDBUserObject (uop, projectID, parentID);
      }
    }
  }

  ValNodeFreeData (head);

  return uop;
}

static void ValNodePtrToGenomeprojidDialog (DialoG d, Pointer data)

{
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   list;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (ValNodePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        str = MemNew (StringLen ((CharPtr) list->data.ptrvalue) + 3);
        if (str != NULL) {
          StringCpy (str, (CharPtr) list->data.ptrvalue);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer GenomeprojidDialogToValNodePtr (DialoG d)

{
  Char         ch;
  ValNodePtr   head;
  Int2         j;
  Int2         len;
  ValNodePtr   list;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    list = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        list = ValNodeNew (list);
        if (head == NULL) {
          head = list;
        }
        if (list != NULL) {
          list->choice = 0;
          list->data.ptrvalue = StringSave (str);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 genproj_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 genproj_widths [] = {
  10, 10, 0
};

static DialoG CreateGenomeProjectsDBDialog (GrouP g)

{
  GrouP                      p;
  GenomeprojidUserDialogPtr  rdp;
  GrouP                      x;
  GrouP                      y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rdp = (GenomeprojidUserDialogPtr) MemNew (sizeof (GenomeprojidUserDialog));
  if (rdp == NULL) return NULL;

  SetObjectExtra (p, rdp, NULL);
  rdp->dialog = (DialoG) p;
  rdp->todialog = UserObjectPtrToGenomeprojidDialog;
  rdp->fromdialog = GenomeprojidDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  StaticPrompt (y, "Project ID", 10 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Parent ID", 10 * stdCharWidth, 0, programFont, 'c');

  rdp->ids = CreateTagListDialog (x, 3, 2, -1,
                                  genproj_types, genproj_widths, NULL,
                                  ValNodePtrToGenomeprojidDialog,
                                  GenomeprojidDialogToValNodePtr);

  return (DialoG) p;
}

static void GenomeProjectsDBUserFormMessage (ForM f, Int2 mssg)

{
  GenomeprojidUserFormPtr  rfp;

  rfp = (GenomeprojidUserFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateGenomeProjectsDBDescForm (Int2 left, Int2 top, Int2 width,
                                           Int2 height, CharPtr title, ValNodePtr sdp,
                                           SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN                   b;
  GrouP                    c;
  GrouP                    g;
  GenomeprojidUserFormPtr  rfp;
  StdEditorProcsPtr        sepp;
  WindoW                   w;

  w = NULL;
  rfp = (GenomeprojidUserFormPtr) MemNew (sizeof (GenomeprojidUserForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdDescFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->formmessage = GenomeProjectsDBUserFormMessage;

    rfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      rfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    rfp->data = CreateGenomeProjectsDBDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK GenomeProjectsDBUserGenFunc (Pointer data);
extern Int2 LIBCALLBACK GenomeProjectsDBUserGenFunc (Pointer data)

{
  ObjectIdPtr              oip;
  OMProcControlPtr         ompcp;
  OMUserDataPtr            omudp;
  ObjMgrProcPtr            proc;
  GenomeprojidUserFormPtr  rfp;
  ValNodePtr               sdp;
  SeqEntryPtr              sep;
  UserObjectPtr            uop;
  WindoW                   w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit GenomeProjectsDB User Desc") == 0) {
      rfp = (GenomeprojidUserFormPtr) omudp->userdata.ptrvalue;
      if (rfp != NULL) {
        Select (rfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "GenomeProjectsDB") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateGenomeProjectsDBDescForm (-50, -33, -10, -10,
                                               "Genome Projects DB", sdp, sep,
                                               StdDescFormActnProc);
  rfp = (GenomeprojidUserFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQDESC;
    rfp->this_subtype = Seq_descr_user;
    rfp->procid = ompcp->proc->procid;
    rfp->proctype = ompcp->proc->proctype;
    rfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, rfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) rfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (rfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct dblinkdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG  traceassm;
  DialoG  biosample;
} DblinkDialog, PNTR DblinkDialogPtr;

typedef struct dblinkform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} DblinkForm, PNTR DblinkFormPtr;

static void UserObjectPtrToDblinkDialog (
  DialoG d,
  Pointer data
)

{
  Char             buf [32];
  CharPtr PNTR     cpp;
  UserFieldPtr     curr;
  DblinkDialogPtr  ddp;
  ValNodePtr       head;
  Int4             i, val;
  Int4Ptr          ip;
  Int4             num;
  ObjectIdPtr      oip;
  CharPtr          str;
  UserObjectPtr    uop;
 
  ddp = (DblinkDialogPtr) GetObjectExtra (d);
  if (ddp == NULL) return;

  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "DBLink") != 0) {
    PointerToDialog (ddp->traceassm, NULL);
    PointerToDialog (ddp->biosample, NULL);
    return;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip == NULL) continue;
    if (StringICmp (oip->str, "Trace Assembly Archive") == 0) {
      if (curr->choice == 8) {
        num = curr->num;
        ip = (Int4Ptr) curr->data.ptrvalue;
        if (num > 0 && ip != NULL) {
          head = NULL;
          for (i = 0; i < num; i++) {
            val = ip [i];
            if (val > 0) {
              sprintf (buf, "%ld", (long) val);
              ValNodeCopyStr (&head, 0, buf);
            }
          }
          if (head != NULL) {
            PointerToDialog (ddp->traceassm, (Pointer) head);
          }
          ValNodeFreeData (head);
        }
      }
    } else if (StringICmp (oip->str, "Bio Sample") == 0) {
      if (curr->choice == 7) {
        num = curr->num;
        cpp = (CharPtr PNTR) curr->data.ptrvalue;
        if (num > 0 && cpp != NULL) {
          head = NULL;
          for (i = 0; i < num; i++) {
            str = cpp [i];
            if (StringDoesHaveText (str)) {
              ValNodeCopyStr (&head, 0, str);
            }
          }
          if (head != NULL) {
            PointerToDialog (ddp->biosample, (Pointer) head);
          }
          ValNodeFreeData (head);
        }
      }
    }
  }
}

static Pointer DblinkDialogToUserObjectPtr (
  DialoG d
)

{
  CharPtr PNTR     cpp;
  DblinkDialogPtr  ddp;
  ValNodePtr       head, vnp;
  Int4             i, num;
  Int4Ptr          ipp;
  Boolean          okay = FALSE;
  CharPtr          str;
  UserObjectPtr    uop;
  long int         val;

  ddp = (DblinkDialogPtr) GetObjectExtra (d);
  if (ddp == NULL) return NULL;

  uop = CreateDBLinkUserObject ();
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (ddp->traceassm);
  if (head != NULL) {
    num = 0;
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      num++;
    }
    if (num > 0) {
      ipp = (Int4Ptr) MemNew (sizeof (Int4) * num);
      if (ipp != NULL) {
        i = 0;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (sscanf (str, "%ld", &val) == 1) {
            ipp [i] = (Int4) val;
            i++;
          }
        }
        if (i > 0) {
          AddTraceAssemblyIDsToDBLinkUserObject (uop, i, ipp);
          okay = TRUE;
        }
      }
    }
  }
  ValNodeFreeData (head);

  head = (ValNodePtr) DialogToPointer (ddp->biosample);
  if (head != NULL) {
    num = 0;
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      num++;
    }
    if (num > 0) {
      cpp = (CharPtr PNTR) MemNew (sizeof (CharPtr) * num);
      if (cpp != NULL) {
        i = 0;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          cpp [i] = str;
          i++;
        }
        if (i > 0) {
          AddBioSampleIDsToDBLinkUserObject (uop, i, cpp);
          okay = TRUE;
        }
      }
    }
  }
  ValNodeFreeData (head);

  if (! okay) {
    uop = UserObjectFree (uop);
  }

  return uop;
}

static DialoG CreateDblinkDialog (
  GrouP g
)

{
  DblinkDialogPtr  ddp;
  GrouP            p, x;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ddp = (DblinkDialogPtr) MemNew (sizeof (DblinkDialog));
  if (ddp == NULL) return NULL;

  SetObjectExtra (p, ddp, NULL);
  ddp->dialog = (DialoG) p;
  ddp->todialog = UserObjectPtrToDblinkDialog;
  ddp->fromdialog = DblinkDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 8, NULL);

  StaticPrompt (x, "Trace Assembly", 10 * stdCharWidth, 0, programFont, 'c');
  ddp->traceassm = CreateVisibleStringDialog (x, 3, -1, 15);

  StaticPrompt (x, "Bio Sample", 10 * stdCharWidth, 0, programFont, 'c');
  ddp->biosample = CreateVisibleStringDialog (x, 3, -1, 15);

  return (DialoG) p;
}

static void DblinkFormMessage (
  ForM f,
  Int2 mssg
)

{
  DblinkFormPtr  dfp;

  dfp = (DblinkFormPtr) GetObjectExtra (f);
  if (dfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (dfp->appmessage != NULL) {
          dfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateDblinkDescForm (
  Int2 left,
  Int2 top,
  Int2 width,
  Int2 height,
  CharPtr title,
  ValNodePtr sdp,
  SeqEntryPtr sep,
  FormActnFunc actproc
)

{
  ButtoN             b;
  DblinkFormPtr      dfp;
  GrouP              c, g;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  dfp = (DblinkFormPtr) MemNew (sizeof (DblinkForm));

  if (dfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, dfp, StdDescFormCleanupProc);
    dfp->form = (ForM) w;
    dfp->actproc = actproc;
    dfp->formmessage = DblinkFormMessage;

    dfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      dfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    dfp->data = CreateDblinkDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, dfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }

  return (ForM) w;
}

extern Int2 LIBCALLBACK DBlinkUserGenFunc (Pointer data);
extern Int2 LIBCALLBACK DBlinkUserGenFunc (Pointer data)

{
  ObjectIdPtr              oip;
  OMProcControlPtr         ompcp;
  OMUserDataPtr            omudp;
  ObjMgrProcPtr            proc;
  DblinkFormPtr            dfp;
  ValNodePtr               sdp;
  SeqEntryPtr              sep;
  UserObjectPtr            uop;
  WindoW                   w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit DBLink User Desc") == 0) {
      dfp = (DblinkFormPtr) omudp->userdata.ptrvalue;
      if (dfp != NULL) {
        Select (dfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "DBLink") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateDblinkDescForm (-50, -33, -10, -10,
                                     "DBLink", sdp, sep,
                                     StdDescFormActnProc);
  dfp = (DblinkFormPtr) GetObjectExtra (w);
  if (dfp != NULL) {
    dfp->input_entityID = ompcp->input_entityID;
    dfp->input_itemID = ompcp->input_itemID;
    dfp->input_itemtype = ompcp->input_itemtype;
    dfp->this_itemtype = OBJ_SEQDESC;
    dfp->this_subtype = Seq_descr_user;
    dfp->procid = ompcp->proc->procid;
    dfp->proctype = ompcp->proc->proctype;
    dfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, dfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) dfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (dfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (dfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) dfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data);

#define REFGENE_ASSEMBLY   1
#define REFGENE_RELATED    2
#define REFGENE_SPLICEVAR  3
#define REFGENE_RELDREK    4
#define REFGENE_REJECT     5
#define REFGENE_UNKNOWN    6

typedef struct refgeneuserdialog {
  DIALOG_MESSAGE_BLOCK
  GrouP         status;
  ButtoN        generated;
  TexT          curator;
  TexT          url;
  TexT          source;
  Int2          indexer;
  DialoG        fields;
  ButtoN        pipebtn;
} RefgeneUserDialog, PNTR RefgeneUserDialogPtr;

typedef struct refgeneuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} RefgeneUserForm, PNTR RefgeneUserFormPtr;

static ENUM_ALIST(refgene_alist)
  {" ",          0},
  {"Assembly",    REFGENE_ASSEMBLY},
  {"Related",     REFGENE_RELATED},
  {"SpliceVar",   REFGENE_SPLICEVAR},
  {"RelatedDrek", REFGENE_RELDREK},
  {"Reject",      REFGENE_REJECT},
  {"Unknown",     REFGENE_UNKNOWN},
END_ENUM_ALIST

static Uint2 refgene_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP
};

static Uint2 refgene_widths [] = {
  9, 7, 7, 7, 8, 0
};

static EnumFieldAssocPtr refgene_popups [] = {
  NULL, NULL, NULL, NULL, NULL, refgene_alist
};

static CharPtr refgene_labels [] = {
  "", "Assembly", "Related", "SpliceVar", "RelatedDrek", "Reject", "Unknown", NULL
};

static CharPtr refgene_fields [] = {
  "Accession", "GI", "Start", "Stop", "Comment", "Type", NULL
};

static void AccessionUserFieldPtrToVisStringDialog (DialoG d, Pointer data)

{
  CharPtr       accession;
  CharPtr       comment;
  UserFieldPtr  curr;
  UserFieldPtr  entry;
  Int2          field;
  Char          fm [16];
  Int4          from;
  Int4          gi;
  ValNodePtr    head;
  Int2          i;
  Int2          j;
  CharPtr       name;
  ObjectIdPtr   oip;
  CharPtr       str;
  TagListPtr    tlp;
  Int4          to;
  Char          tu [16];
  UserFieldPtr  ufp;
  ValNodePtr    vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;
  str = MemNew (sizeof (Char) * 1024);
  head = NULL;
  curr = (UserFieldPtr) data;
  while (curr != NULL) {
    oip = curr->label;
    if (oip != NULL) {
      field = 0;
      for (i = REFGENE_ASSEMBLY; i <= REFGENE_UNKNOWN; i++) {
        if (StringICmp (oip->str, refgene_labels [i]) == 0 && curr->choice == 11) {
          field = i;
        }
      }
      if (field > 0) {
        entry = (UserFieldPtr) curr->data.ptrvalue;
        while (entry != NULL && entry->choice == 11) {
          accession = NULL;
          comment = NULL;
          name = NULL;
          gi = 0;
          from = 0;
          to = 0;
          ufp = (UserFieldPtr) entry->data.ptrvalue;
          while (ufp != NULL) {
            oip = ufp->label;
            if (oip != NULL && oip->str != NULL) {
              if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                accession = (CharPtr) ufp->data.ptrvalue;
              } else if (StringICmp (oip->str, "gi") == 0 && ufp->choice == 2) {
                gi = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
                from = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
                to = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "comment") == 0 && ufp->choice == 1) {
                comment = (CharPtr) ufp->data.ptrvalue;
              } else if (StringICmp (oip->str, "name") == 0 && ufp->choice == 1) {
                name = (CharPtr) ufp->data.ptrvalue;
              }
            }
            ufp = ufp->next;
          }
          if (accession != NULL) {
            if (comment == NULL) {
              comment = "";
            }
            fm [0] = '\0';
            tu [0] = '\0';
            if (from > 0 && to > 0) {
              sprintf (fm, "%ld", (long) from);
              sprintf (tu, "%ld", (long) to);
            }
            sprintf (str, "%s\t%ld\t%s\t%s\t%s\t%d\n", accession,
                     (long) gi, fm, tu, comment, (int) field);
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = StringSave (str);
            }
          } else if (name != NULL) {
            sprintf (str, "\t\t\t\t%s\t%d\n", name, (int) field);
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = StringSave (str);
            }
          }
          entry = entry->next;
        }
      }
    }
    curr = curr->next;
  }
  MemFree (str);
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

static Pointer VisStringDialogToUserFieldPtr (DialoG d)

{
  return NULL;
}

static void UserObjectPtrToRefGeneDialog (DialoG d, Pointer data)

{
  UserFieldPtr          curr;
  Boolean               gen;
  ObjectIdPtr           oip;
  RefgeneUserDialogPtr  rdp;
  Int2                  status = 0;
  CharPtr               str;
  UserObjectPtr         uop;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return;
  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "RefGeneTracking") != 0) {
    SetValue (rdp->status, 0);
    PointerToDialog (rdp->fields, NULL);
    return;
  }
  PointerToDialog (rdp->fields, uop->data);
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Status") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    if (StringICmp (str, "Inferred") == 0) {
      status = 1;
    } else if (StringICmp (str, "Predicted") == 0) {
      status = 2;
    } else if (StringICmp (str, "Provisional") == 0) {
      status = 3;
    } else if (StringICmp (str, "Validated") == 0) {
      status = 4;
    } else if (StringICmp (str, "Reviewed") == 0) {
      status = 5;
    } else if (StringICmp (str, "Model") == 0) {
      status = 6;
    } else if (StringICmp (str, "WGS") == 0) {
      status = 7;
    } else if (StringICmp (str, "Pipeline") == 0) {
      status = 8;
      SafeEnable (rdp->pipebtn);
    }
  }
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Generated") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 4) {
    gen = curr->data.boolvalue;
    SetStatus (rdp->generated, gen);
  }
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Collaborator") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    SetTitle (rdp->curator, str);
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "CollaboratorURL") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    SetTitle (rdp->url, str);
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "GenomicSource") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    SetTitle (rdp->source, str);
  }
  SetValue (rdp->status, status);
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Indexer") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 2) {
    rdp->indexer = (Int2) curr->data.intvalue;
  }
}

static void AddIndexerToRefGeneTrackUserObject (UserObjectPtr uop, Int2 indexer)

{
  UserFieldPtr  curr;
  ObjectIdPtr   oip;

  if (uop == NULL || indexer < 1) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "RefGeneTracking") != 0) return;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Indexer") == 0) {
      break;
    }
  }

  if (curr == NULL) {
    curr = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave ("Indexer");
    curr->label = oip;
    curr->choice = 2; /* integer */

    /* link indexer at beginning of list */

    curr->next = uop->data;
    uop->data = curr;
  }

  if (curr == NULL || curr->choice != 2) return;

  /* replace any existing indexer indication */

  curr->data.intvalue = (Int4) indexer;
}

static Pointer RefGeneDialogToUserObjectPtr (DialoG d)

{
  Char                  ch;
  Char                  curator [256];
  Int2                  i;
  Uint2                 j;
  size_t                len;
  Int4                  num [6];
  Boolean               okay;
  RefgeneUserDialogPtr  rdp;
  Char                  source [64];
  Int2                  status;
  CharPtr               str;
  TagListPtr            tlp;
  CharPtr               txt [6];
  UserObjectPtr         uop;
  Char                  url [512];
  long int              val;
  ValNodePtr            vnp;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return NULL;

  uop = CreateRefGeneTrackUserObject ();
  if (uop == NULL) return NULL;

  status = GetValue (rdp->status);
  if (status == 1) {
    AddStatusToRefGeneTrackUserObject (uop, "Inferred");
  } else if (status == 2) {
    AddStatusToRefGeneTrackUserObject (uop, "Predicted");
  } else if (status == 3) {
    AddStatusToRefGeneTrackUserObject (uop, "Provisional");
  } else if (status == 4) {
    AddStatusToRefGeneTrackUserObject (uop, "Validated");
  } else if (status == 5) {
    AddStatusToRefGeneTrackUserObject (uop, "Reviewed");
  } else if (status == 6) {
    AddStatusToRefGeneTrackUserObject (uop, "Model");
  } else if (status == 7) {
    AddStatusToRefGeneTrackUserObject (uop, "WGS");
  } else if (status == 8) {
    AddStatusToRefGeneTrackUserObject (uop, "Pipeline");
  }

  GetTitle (rdp->source, source, sizeof (source));
  if (! StringHasNoText (source)) {
    AddSourceToRefGeneTrackUserObject (uop, source);
  }

  if (GetStatus (rdp->generated)) {
    AddGeneratedToRefGeneTrackUserObject (uop, TRUE);
  }

  GetTitle (rdp->curator, curator, sizeof (curator));
  if (! StringHasNoText (curator)) {
    AddCuratorToRefGeneTrackUserObject (uop, curator);
  }

  GetTitle (rdp->url, url, sizeof (url));
  if (! StringHasNoText (url)) {
    AddCuratorURLToRefGeneTrackUserObject (uop, url);
  }

  if (rdp->indexer > 0) {
    AddIndexerToRefGeneTrackUserObject (uop, rdp->indexer);
  }

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL && tlp->vnp != NULL) {
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        for (j = 0; j < 6; j++) {
          txt [j] = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, j);
          num [j] = 0;
        }
        for (j = 1; j < 4; j++) {
          num [j] = 0;
          if (txt [j] != NULL && sscanf (txt [j], "%ld", &val) == 1) {
            num [j] = val;
          }
        }
        if (txt [5] != NULL && sscanf (txt [5], "%ld", &val) == 1) {
          num [5] = val;
        }
        i = num [5];
        if (i >= REFGENE_ASSEMBLY && i <= REFGENE_UNKNOWN) {
          if (! StringHasNoText (txt [0])) {
            AddAccessionToRefGeneTrackUserObject (uop, refgene_labels [i],
                                                  txt [0], num [1], num [2],
                                                  num [3], txt [4]);
          } else if (! StringHasNoText (txt [4])) {
            /* comment by itself goes into name */
            AddAccessionToRefGeneTrackUserObject (uop, refgene_labels [i],
                                                  NULL, num [1], num [2],
                                                  num [3], txt [4]);
          }
        }
        for (j = 0; j < 6; j++) {
          txt [j] = MemFree (txt [j]);
        }
      }
    }
  }

  return uop;
}

static DialoG CreateRefGeneDialog (GrouP g)

{
  Int2                  i;
  PrompT                lastppt;
  GrouP                 p;
  PrompT                ppt;
  GrouP                 q;
  RefgeneUserDialogPtr  rdp;
  TagListPtr            tlp;
  GrouP                 x;
  GrouP                 y;
  GrouP                 z;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rdp = (RefgeneUserDialogPtr) MemNew (sizeof (RefgeneUserDialog));
  if (rdp == NULL) return NULL;

  SetObjectExtra (p, rdp, NULL);
  rdp->dialog = (DialoG) p;
  rdp->todialog = UserObjectPtrToRefGeneDialog;
  rdp->fromdialog = RefGeneDialogToUserObjectPtr;

  x = HiddenGroup (p, 4, 0, NULL);
  /* StaticPrompt (x, "Status", 0, stdLineHeight, programFont, 'l'); */
  rdp->status = HiddenGroup (x, 8, 0, NULL);
  SetObjectExtra (rdp->status, rdp, NULL);
  RadioButton (rdp->status, "Inferred");
  RadioButton (rdp->status, "Predicted");
  RadioButton (rdp->status, "Provisional");
  RadioButton (rdp->status, "Validated");
  RadioButton (rdp->status, "Reviewed");
  RadioButton (rdp->status, "Model");
  RadioButton (rdp->status, "WGS");
  rdp->pipebtn = RadioButton (rdp->status, "Pipeline");
  Disable (rdp->pipebtn);

  y = HiddenGroup (p, 6, 0, NULL);
  rdp->generated = CheckBox (y, "Generated", NULL);
  z = HiddenGroup (y, 2, 0, NULL);
  StaticPrompt (z, "Curator", 0, dialogTextHeight, programFont, 'l');
  rdp->curator = DialogText (z, "", 14, NULL);
  StaticPrompt (z, "URL", 0, dialogTextHeight, programFont, 'r');
  rdp->url = DialogText (z, "", 14, NULL);
  StaticPrompt (y, "Genomic Source", 0, dialogTextHeight, programFont, 'l');
  rdp->source = DialogText (y, "", 7, NULL);

  rdp->indexer = 0;

  q = HiddenGroup (p, -7, 0, NULL);
  lastppt = NULL;
  ppt = NULL;
  for (i = 0; i < 6; i++) {
    lastppt = ppt;
    ppt = StaticPrompt (q, refgene_fields [i], refgene_widths [i] * stdCharWidth, 0, systemFont, 'c');
  }
  rdp->fields = CreateTagListDialog (p, 6, 6, STD_TAG_SPACING,
                                     refgene_types, refgene_widths, refgene_popups,
                                     AccessionUserFieldPtrToVisStringDialog,
                                     VisStringDialogToUserFieldPtr);

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL) {
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [4], (HANDLE) lastppt, NULL);
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [5], (HANDLE) ppt, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) q, (HANDLE) rdp->fields, NULL);
  return (DialoG) p;
}

static void RefgeneUserFormMessage (ForM f, Int2 mssg)

{
  RefgeneUserFormPtr  rfp;

  rfp = (RefgeneUserFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateRefGeneDescForm (Int2 left, Int2 top, Int2 width,
                                   Int2 height, CharPtr title, ValNodePtr sdp,
                                   SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN              b;
  GrouP               c;
  GrouP               g;
  RefgeneUserFormPtr  rfp;
  StdEditorProcsPtr   sepp;
  WindoW              w;

  w = NULL;
  rfp = (RefgeneUserFormPtr) MemNew (sizeof (RefgeneUserForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdDescFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->formmessage = RefgeneUserFormMessage;

    rfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      rfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    rfp->data = CreateRefGeneDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data)

{
  ObjectIdPtr         oip;
  OMProcControlPtr    ompcp;
  OMUserDataPtr       omudp;
  ObjMgrProcPtr       proc;
  RefgeneUserFormPtr  rfp;
  ValNodePtr          sdp;
  SeqEntryPtr         sep;
  UserObjectPtr       uop;
  WindoW              w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit RefGene UserTrack Desc") == 0) {
      rfp = (RefgeneUserFormPtr) omudp->userdata.ptrvalue;
      if (rfp != NULL) {
        Select (rfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "RefGeneTracking") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateRefGeneDescForm (-50, -33, -10, -10,
                                      "Reference Gene Tracking", sdp, sep,
                                      StdDescFormActnProc);
  rfp = (RefgeneUserFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQDESC;
    rfp->this_subtype = Seq_descr_user;
    rfp->procid = ompcp->proc->procid;
    rfp->proctype = ompcp->proc->proctype;
    rfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, rfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) rfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (rfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}

extern Int2 LIBCALLBACK StruCommUserGenFunc (Pointer data);

typedef struct strucommuserdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        fields;
} StruCommUserDialog, PNTR StruCommUserDialogPtr;

typedef struct strucommuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} StruCommUserForm, PNTR StruCommUserFormPtr;

static void UserObjectPtrToStruCommDialog (DialoG d, Pointer data)

{
  UserFieldPtr           curr;
  CharPtr                field;
  ValNodePtr             head = NULL;
  ObjectIdPtr            oip;
  StruCommUserDialogPtr  sdp;
  CharPtr                str;
  CharPtr                tmp;
  UserObjectPtr          uop;

  sdp = (StruCommUserDialogPtr) GetObjectExtra (d);
  if (sdp == NULL) return;

  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    PointerToDialog (sdp->fields, NULL);
    return;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
   if (curr->choice != 1) continue;
    oip = curr->label;
    if (oip == NULL) continue;
    field = oip->str;
    if (StringHasNoText (field)) continue;
    str = (CharPtr) curr->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    tmp = MemNew (StringLen (field) + StringLen (str) + 5);
    if (tmp == NULL) continue;
    sprintf (tmp, "%s\t%s", field, str);
    ValNodeAddStr (&head, 0, (Pointer) tmp);
  }

  PointerToDialog (sdp->fields, (Pointer) head);
  ValNodeFreeData (head);
}

static void FixSpecialCharactersInStructuredCommentUserObject (UserObjectPtr uop, BoolPtr changed)
{
  UserFieldPtr           curr;
  CharPtr                field;
  ValNodePtr             find_list = NULL;
  ObjectIdPtr            oip;
  CharPtr                str;
 

  if (changed != NULL) {
    *changed = FALSE;
  }
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
   if (curr->choice != 1) continue;
    oip = curr->label;
    if (oip == NULL) continue;
    field = oip->str;
    if (StringHasNoText (field)) continue;
    str = (CharPtr) curr->data.ptrvalue;
    if (StringHasNoText (str)) continue;

    SpecialCharFindWithContext ((CharPtr PNTR) (&(oip->str)), &find_list, NULL, NULL);
    SpecialCharFindWithContext ((CharPtr PNTR) (&(curr->data.ptrvalue)), &find_list, NULL, NULL);
  }
  FixSpecialCharactersForStringsInList (find_list, "Special characters are not permitted.", TRUE);  
  if (find_list != NULL)
  {
    if (changed != NULL) {
      *changed = TRUE;
    }
    find_list = FreeContextList (find_list);
  }
}


static Pointer StruCommDialogToUserObjectPtr (DialoG d)

{
  CharPtr                field;
  ValNodePtr             head;
  CharPtr                item;
  UserObjectPtr          uop;
  StruCommUserDialogPtr  sdp;
  CharPtr                str;
  ValNodePtr             vnp;
  Boolean                fixed_special = FALSE;

  sdp = (StruCommUserDialogPtr) GetObjectExtra (d);
  if (sdp == NULL) return NULL;

  uop = CreateStructuredCommentUserObject (NULL, NULL);
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (sdp->fields);
  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    field = ExtractTagListColumn (str, 0);
    item = ExtractTagListColumn (str, 1);
    if (StringDoesHaveText (field) && StringDoesHaveText (item)) {
      AddItemStructuredCommentUserObject (uop, field, item);
    }
    MemFree (field);
    MemFree (item);
  }

  ValNodeFreeData (head);

  FixSpecialCharactersInStructuredCommentUserObject (uop, &fixed_special);
  if (fixed_special) {
    PointerToDialog (d, uop);
  }

  return uop;
}

static void ValNodePtrToStruCommDialog (DialoG d, Pointer data)

{
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   list;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (ValNodePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        str = MemNew (StringLen ((CharPtr) list->data.ptrvalue) + 3);
        if (str != NULL) {
          StringCpy (str, (CharPtr) list->data.ptrvalue);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer StruCommDialogToValNodePtr (DialoG d)

{
  Char         ch;
  ValNodePtr   head;
  Int2         j;
  Int2         len;
  ValNodePtr   list;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    list = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        list = ValNodeNew (list);
        if (head == NULL) {
          head = list;
        }
        if (list != NULL) {
          list->choice = 0;
          list->data.ptrvalue = StringSave ((CharPtr) vnp->data.ptrvalue);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 strccmm_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 strccmm_widths [] = {
  16, 16, 0
};

static DialoG CreateStruCommDialog (GrouP g)

{
  StruCommUserDialogPtr  sdp;
  GrouP                  p;
  GrouP                  x;
  GrouP                  y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  sdp = (StruCommUserDialogPtr) MemNew (sizeof (StruCommUserDialog));
  if (sdp == NULL) return NULL;

  SetObjectExtra (p, sdp, NULL);
  sdp->dialog = (DialoG) p;
  sdp->todialog = UserObjectPtrToStruCommDialog;
  sdp->fromdialog = StruCommDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  StaticPrompt (y, "Field", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Text", 16 * stdCharWidth, 0, programFont, 'c');
  sdp->fields = CreateTagListDialog (x, 8, 2, -1,
                                     strccmm_types, strccmm_widths, NULL,
                                     ValNodePtrToStruCommDialog,
                                     StruCommDialogToValNodePtr);

  return (DialoG) p;
}

static void StruCommUserFormMessage (ForM f, Int2 mssg)

{
  StruCommUserFormPtr  sfp;

  sfp = (StruCommUserFormPtr) GetObjectExtra (f);
  if (sfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (sfp->appmessage != NULL) {
          sfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateStruCommDescForm (Int2 left, Int2 top, Int2 width,
                                     Int2 height, CharPtr title, ValNodePtr sdp,
                                     SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN               b;
  GrouP                c;
  GrouP                g;
  StdEditorProcsPtr    sepp;
  StruCommUserFormPtr  sfp;
  WindoW               w;

  w = NULL;
  sfp = (StruCommUserFormPtr) MemNew (sizeof (StruCommUserForm));
  if (sfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, sfp, StdDescFormCleanupProc);
    sfp->form = (ForM) w;
    sfp->actproc = actproc;
    sfp->formmessage = StruCommUserFormMessage;

    sfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      sfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    sfp->data = CreateStruCommDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, sfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK StruCommUserGenFunc (Pointer data)

{
  ObjectIdPtr          oip;
  OMProcControlPtr     ompcp;
  OMUserDataPtr        omudp;
  ObjMgrProcPtr        proc;
  ValNodePtr           sdp;
  SeqEntryPtr          sep;
  StruCommUserFormPtr  sfp;
  UserObjectPtr        uop;
  WindoW               w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit StructuredComment User Desc") == 0) {
      sfp = (StruCommUserFormPtr) omudp->userdata.ptrvalue;
      if (sfp != NULL) {
        Select (sfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "StructuredComment") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateStruCommDescForm (-50, -33, -10, -10,
                                       "Structured Comment", sdp, sep,
                                       StdDescFormActnProc);
  sfp = (StruCommUserFormPtr) GetObjectExtra (w);
  if (sfp != NULL) {
    sfp->input_entityID = ompcp->input_entityID;
    sfp->input_itemID = ompcp->input_itemID;
    sfp->input_itemtype = ompcp->input_itemtype;
    sfp->this_itemtype = OBJ_SEQDESC;
    sfp->this_subtype = Seq_descr_user;
    sfp->procid = ompcp->proc->procid;
    sfp->proctype = ompcp->proc->proctype;
    sfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, sfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) sfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (sfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (sfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) sfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


/*
static void TestGeneRefStuff (void)

{
  UserObjectPtr uop;
  ValNodePtr    sdp;

  uop = CreateRefGeneTrackUserObject ();
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "U12345", 57, 29, 1995);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "L97531", 142, 66, 963);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "M66778", 823, 7677, 343);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "P34345", 445, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Reject", "S19635", 1765, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "Q14884", 664, 35, 97);
  sdp = ValNodeNew (NULL);
  sdp->choice = Seq_descr_user;
  sdp->data.ptrvalue = (Pointer) uop;
  if (! ObjMgrRegister (OBJ_SEQDESC, (Pointer) sdp)) {
     ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
  }
}
*/

#define CKA_GAPLEN  50 /* max allowed unaligned gap size */

typedef struct cka_acc {
   CharPtr      accession;
   SeqIdPtr     sip_whole;
   SeqAlignPtr  sap;
   Int4         start_acc;
   Int4         stop_acc;
   Int4         start_seq;
   Int4         stop_seq;
   Uint1        strand;
   Int4         num;
   struct cka_acc PNTR next;
} CKA_Acc, PNTR CKA_AccPtr;

static Int4     CKA_blast_wordsize;
static FloatHi  CKA_blast_expect_value;
static Boolean  CKA_blast_allow_repeats;
static Int4     CKA_blast_detailed_wordsize;
static FloatHi  CKA_blast_detailed_expect_value;
static Boolean  CKA_blast_detailed_allow_repeats;

static SeqAlignPtr CKA_MakeAlign(BioseqPtr bsp, CKA_AccPtr acc_head, LogInfoPtr lip);

static Boolean SPI_GetAccessionFromSeqId(SeqIdPtr sip, Int4Ptr gi, CharPtr PNTR id)
{
   Boolean numeric_id_type = FALSE;
   Int2 id_len;
   GiimPtr gip;
   ObjectIdPtr oip;
   TextSeqIdPtr textsip;
   DbtagPtr dbtag;
   PatentSeqIdPtr psip;
   PDBSeqIdPtr pdbsip;

   *id = NULL;
   *gi = 0;

   switch (sip->choice) {
   case SEQID_GI: case SEQID_GIBBSQ: case SEQID_GIBBMT:
      *gi = sip->data.intvalue;
      numeric_id_type = TRUE;
      break;
   case SEQID_GIIM:
      gip = (GiimPtr) sip->data.ptrvalue;
      *gi = gip->id;
      numeric_id_type = TRUE;
      break;
   case SEQID_LOCAL:
      oip = (ObjectIdPtr) sip->data.ptrvalue;

      if (oip->str) {
         id_len = StringLen(oip->str);
         *id = (CharPtr) MemNew(id_len+1);
         sprintf(*id, "%s", oip->str);
      } else {
         *id = (CharPtr) MemNew(6);
         sprintf(*id, "%d", oip->id);
      }
      break;
   case SEQID_GENBANK: case SEQID_EMBL: case SEQID_PIR: case SEQID_TPG: case SEQID_TPE: case SEQID_TPD:
   case SEQID_SWISSPROT: case SEQID_DDBJ: case SEQID_PRF:
   case SEQID_OTHER: case SEQID_GPIPE:
      textsip = (TextSeqIdPtr)sip->data.ptrvalue;
      id_len = StringLen(textsip->accession);
      *id = (CharPtr) MemNew(id_len+1);
      if (textsip->version > 0)
         sprintf(*id, "%s.%d", textsip->accession, textsip->version);
      else
         sprintf(*id, "%s", textsip->accession);
      break;
   case SEQID_GENERAL:
      dbtag = (DbtagPtr) sip->data.ptrvalue;
      if (dbtag->tag->str == NULL) {
         numeric_id_type = TRUE;
         *gi = dbtag->tag->id;
      } else {
         id_len = StringLen(dbtag->tag->str);
         *id = (CharPtr) MemNew(id_len+1);
         sprintf(*id, "%s", dbtag->tag->str);
      }
      break;
   case SEQID_PATENT:
      psip = (PatentSeqIdPtr) sip->data.ptrvalue;
      *gi = (Int4) psip->seqid;
      numeric_id_type = TRUE;
      break;
   case SEQID_PDB:
      pdbsip = (PDBSeqIdPtr) sip->data.ptrvalue;
      id_len = StringLen(pdbsip->mol);
      *id = (CharPtr) MemNew(id_len+4);
      sprintf(*id, "%s%d", pdbsip->mol, pdbsip->chain);
      break;
   default: break;
   }

   return numeric_id_type;
}

static void CKA_FindAllTpaDescr(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
   CKA_AccPtr         acc;
   CKA_AccPtr         PNTR acc_head;
   CKA_AccPtr         acc_prev;
   BioseqPtr          bsp;
   SeqMgrDescContext  context;
   UserFieldPtr       curr;
   ObjectIdPtr        oip;
   SeqDescrPtr        sdp;
   UserFieldPtr       ufp;
   UserObjectPtr      uop;

   acc_head = (CKA_AccPtr PNTR)data;
   acc_prev = *acc_head;
   while (acc_prev != NULL && acc_prev->next != NULL)
   {
      acc_prev = acc_prev->next;
   }
   sdp = NULL;
   if (IS_Bioseq(sep))
   {
      bsp = (BioseqPtr)sep->data.ptrvalue;
      if (ISA_na(bsp->mol))
      {
         while ((sdp = SeqMgrGetNextDescriptor(bsp, sdp, Seq_descr_user, &context)) != NULL)
         {
            uop = (UserObjectPtr)sdp->data.ptrvalue;
            if (!StringICmp(uop->type->str, "TpaAssembly"))
            {
               for (curr = uop->data; curr != NULL; curr = curr->next)
               {
                  if (curr->choice != 11) continue;
                  
                  acc = (CKA_AccPtr)MemNew(sizeof(CKA_Acc));
                  acc->sip_whole = SeqIdSetDup(bsp->id);
                  /* will use these to mark the span for blast2seq */
                  acc->start_acc = acc->stop_acc = -1;
                  if (acc_prev == NULL)
                    *acc_head = acc_prev = acc;
                  else {
                    acc_prev->next = acc;
                    acc_prev = acc;
                  }
                  
                  for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
                    oip = ufp->label;
                    if (oip == NULL) continue;
                    if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                      acc->accession = StringSave((CharPtr)ufp->data.ptrvalue);
                    } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
                      acc->start_acc = (Int4) ufp->data.intvalue;
                    } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
                      acc->stop_acc = (Int4) ufp->data.intvalue;
                    }
                  }
               }
            }
         }
      }
   }
}

static int LIBCALLBACK CKA_SortAccs(VoidPtr ptr1, VoidPtr ptr2)
{
   CKA_AccPtr  acc1;
   CKA_AccPtr  acc2;

   acc1 = *((CKA_AccPtr PNTR)ptr1);
   acc2 = *((CKA_AccPtr PNTR)ptr2);
   if (acc1->start_seq < acc2->start_seq)
      return -1;
   else if (acc1->start_seq > acc2->start_seq)
      return 1;
   else if (acc1->stop_seq < acc2->stop_seq)
      return -1;
   else if (acc1->stop_seq > acc2->stop_seq)
      return 1;
   else
      return 0; /* no alignment */
}

static SeqIdPtr SqnSeqIdFindBestAccession (SeqIdPtr sip)
{
	Uint1 order[NUM_SEQID];

	if (sip == NULL)
		return NULL;
	SeqIdBestRank(order, NUM_SEQID);
        order[SEQID_GI]=order[SEQID_LOCAL]+2;
        order[SEQID_PATENT]=order[SEQID_LOCAL]+1;
	return SeqIdSelect (sip, order, NUM_SEQID);
}


static Boolean ValidateTPAHistAlign (BioseqPtr bsp, ValNodePtr PNTR errors)
{
  ValNodePtr new_errors;
  Boolean    retval = TRUE;
  SeqAlignPtr salp;

  if (bsp == NULL || bsp->hist == NULL || bsp->hist->assembly == NULL) {
    return FALSE;
  }

  for (salp = bsp->hist->assembly; salp != NULL; salp = salp->next) {
    AlnMgr2IndexSingleChildSeqAlign(salp);
  }

  new_errors = ReportCoverageForBioseqSeqHist (bsp);
  if (new_errors != NULL) {
    ValNodeLink (errors, new_errors);
    retval = FALSE;
  }
  return retval;
}

 
static Boolean CKA_ValidateSeqAlign(SeqAlignPtr sap, CKA_AccPtr acc_head, Int4 bioseqlen, ValNodePtr PNTR errors)
{
   CKA_AccPtr        acc;
   CKA_AccPtr        PNTR accarray;
   AMAlignIndex2Ptr  amaip;
   Int4              first, first_align;
   Boolean           found;
   Int4              gi;
   Int4              i;
   Int4              j;
   Int4              k;
   Int4              last;
   Int4              longest;
   Int4              max;
   Int4              n;
   Int4              prev;
   Boolean           retval = TRUE;
   CharPtr           textid;
   Char              textid2[42];
   CharPtr           err_msg;
   CharPtr           no_cover_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s aligns to %d-%d but next aln is %s to %d-%d\n"; 
   CharPtr           no_cover_ok_gap_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s aligns to %d-%d but the next aln is %s to %d-%d;\n the gap is less than %d and is acceptable.\n";
   CharPtr           bad_start_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s (the first aln) starts at position %d\n";
   CharPtr           bad_start_ok_gap_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s (the first alignment) starts at position %d, but the gap is less than %d and is acceptable.\n";
   CharPtr           bad_end_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s (the last aln) goes to %d, bioseq length is %d\n";
   CharPtr           bad_end_ok_gap_fmt = "Primary accessions do not completely cover the bioseq %s:\n %s (the last alignment) goes to %d, bioseq length is %d, but the gap is less than %d and is acceptable.\n";

   if (sap == NULL || sap->saip == NULL || sap->saip->indextype != INDEX_PARENT || errors == NULL)
      return FALSE;

   amaip = (AMAlignIndex2Ptr)(sap->saip);
   for (i=0; i<amaip->numsaps; i++)
   {
      acc = acc_head;
      found = FALSE;
      while (acc != NULL && !found)
      {
         if (amaip->saps[i] == acc->sap)
            found = TRUE;
         if (!found)
            acc = acc->next;
      }
      if (!found) /* big error */
         return FALSE;
      acc->num = i+1;
      AlnMgr2GetNthSeqRangeInSA(amaip->saps[i], 1, &acc->start_seq, &acc->stop_seq);
      AlnMgr2GetNthSeqRangeInSA(amaip->saps[i], 2, &acc->start_acc, &acc->stop_acc);
      acc->strand = AlnMgr2GetNthStrand(amaip->saps[i], 2);
      acc->start_seq++;
      acc->stop_seq++;
      acc->start_acc++;
      acc->stop_acc++;
   }
   acc = acc_head;
   i = 0;
   while (acc != NULL)
   {
      if (acc->start_seq == 0 && acc->stop_seq == 0)
      {
         AlnMgr2GetNthSeqRangeInSA(acc->sap, 1, &acc->start_seq, &acc->stop_seq);
         AlnMgr2GetNthSeqRangeInSA(acc->sap, 2, &acc->start_acc, &acc->stop_acc);
         acc->strand = AlnMgr2GetNthStrand(acc->sap, 2);
         acc->start_seq++;
         acc->stop_seq++;
         acc->start_acc++;
         acc->stop_acc++;
      }
      if (acc->num == 0)
         acc->num = amaip->numsaps; /* sort these guys all to the end */
      i++;
      acc = acc->next;
   }
   accarray = (CKA_AccPtr PNTR)MemNew(i*sizeof(CKA_AccPtr));
   i = 0;
   acc = acc_head;
   while (acc != NULL)
   {
      accarray[i] = acc;
      i++;
      acc = acc->next;
   }
   HeapSort(accarray, i, sizeof(CKA_AccPtr), CKA_SortAccs);
   n=0;
   while (accarray[n]->sap == NULL && n < i)
   {
      n++;
   }
   SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(accarray[0]->sip_whole), &gi, &textid);
   if (textid == NULL)
   {
      sprintf(textid2, "%d", gi);
      textid = textid2;
   }
   first = last = -1;
   prev = -1;
   retval = TRUE;
   for (j=0; j<i /*&& first <=0*/ ; j++)
   {
      acc = accarray[j];
      if (acc->sap != NULL)
      {
        if (first == -1) {
          first = acc->start_seq;
          first_align = j;
        }
        last = MAX(last, acc->stop_seq);
      } else {
        continue;
      }
      if (prev != -1)
      {
         if (acc->start_seq > prev + CKA_GAPLEN)
         {
            err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_cover_fmt) + StringLen (textid) + StringLen (accarray[j-1]->accession)
                                                         + StringLen (acc->accession) + 60));
            sprintf (err_msg, no_cover_fmt, 
                     textid, accarray[j-1]->accession, accarray[j-1]->start_seq, accarray[j-1]->stop_seq, acc->accession, acc->start_seq, acc->stop_seq);
            ValNodeAddPointer (errors, 0, err_msg);
            retval = FALSE;
         } 
         else if (acc->start_seq > prev)
         {
            err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_cover_ok_gap_fmt) + StringLen (textid)
                                                         + StringLen (accarray[j-1]->accession)
                                                         + StringLen (acc->accession)
                                                         + 75));
            sprintf (err_msg, no_cover_ok_gap_fmt,
                     textid, accarray[j-1]->accession, accarray[j-1]->start_seq, accarray[j-1]->stop_seq, acc->accession, acc->start_seq, acc->stop_seq, CKA_GAPLEN);
            ValNodeAddPointer (errors, 0, err_msg);
            retval = FALSE;
         }
      }
      prev = acc->stop_seq+1;
   }
   if (first != 1 || last != bioseqlen)
   {
      if (first > CKA_GAPLEN)
      {
         err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_start_fmt) + StringLen (textid)
                                                      + StringLen (accarray[first_align]->accession)
                                                     + 15));
         sprintf (err_msg, bad_start_fmt,
                          textid, accarray[first_align]->accession, accarray[first_align]->start_seq);
         ValNodeAddPointer (errors, 0, err_msg);
         retval = FALSE;
      } 
      else if (first != 1)
      {
         err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_start_ok_gap_fmt) 
                                                      + StringLen (textid)
                                                      + StringLen (accarray[first_align]->accession)
                                                      + 30));
         sprintf (err_msg, bad_start_ok_gap_fmt, 
                  textid, accarray[first_align]->accession, accarray[first_align]->start_seq, CKA_GAPLEN);
         ValNodeAddPointer (errors, 0, err_msg);
      }
      max = 0;
      for (k=0; k<i; k++)
      {
         if (accarray[k]->stop_seq > max)
         {
            max = accarray[k]->stop_seq;
            longest = k;
         }
      }
      if (accarray[longest]->stop_seq < bioseqlen-CKA_GAPLEN)
      {
         err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_end_fmt) 
                                                      + StringLen (textid)
                                                      + StringLen (accarray[longest]->accession)
                                                      + 30));
         sprintf (err_msg, bad_end_fmt,
                          textid, accarray[longest]->accession, accarray[longest]->stop_seq, bioseqlen);
         ValNodeAddPointer (errors, 0, err_msg);
         retval = FALSE;
      } 
      else if (accarray[longest]->stop_seq < bioseqlen)
      {
         err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_end_ok_gap_fmt)
                                                      + StringLen (textid)
                                                      + StringLen (accarray[longest]->accession)
                                                      + 45));
         sprintf (err_msg, bad_end_ok_gap_fmt,
                           textid, accarray[longest]->accession, accarray[longest]->stop_seq, bioseqlen, CKA_GAPLEN);
         ValNodeAddPointer (errors, 0, err_msg);
      }
      MemFree(accarray);
      accarray = NULL;
   }
   MemFree(accarray);
   return retval;
}

static void FrameVwr (
  VieweR vwr,
  SegmenT pict
)

{
  RecT  r;

  ResetClip ();
  ObjectRect (vwr, &r);
  FrameRect (&r);
}

static void CKA_ShowAln(SeqAlignPtr sap, CKA_AccPtr acc_head)
{
   CKA_AccPtr   acc;
   BioseqPtr    bsp;
   DenseSegPtr  dsp;
   Boolean      found;
   GrouP        g;
   Int4         gi;
   Int4         i;
   Int4         len;
   Int4         numsaps;
   SegmenT      picture;
   SeqAlignPtr  salp;
   SeqIdPtr     sip;
   Int4         start;
   Int4         start_r;
   Int4         stop;
   Int4         stop_r;
   Char         tmp[42];
   CharPtr      textid;
   Char         textid2[42];
   VieweR       v;
   WindoW       w;

   w = FixedWindow(-1, -1, -1, -1, "TPA display", StdCloseWindowProc);
   g = HiddenGroup(w, 1, 0, NULL);
   v = CreateViewer(g, 750, 300, TRUE, TRUE);
   picture = CreatePicture();
   salp = (SeqAlignPtr)(sap->segs);
   numsaps = 0;
   while (salp != NULL)
   {
      numsaps++;
      salp = salp->next;
   }
   salp = (SeqAlignPtr)(sap->segs);
   numsaps++;
   dsp = (DenseSegPtr)(salp->segs);
   sip = dsp->ids;
   SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(sip), &gi, &textid);
   if (textid == NULL)
   {
      sprintf(textid2, "%d", gi);
      textid = textid2;
   }
   bsp = BioseqLockById(sip);
   len = bsp->length;
   AddRectangle(picture, 0, numsaps*10, (bsp->length*680)/len, numsaps*10-7, 0, TRUE, 0);
   sprintf(tmp, "1");
   AddLabel(picture, 0-10, numsaps*10-3, tmp, 0, 0, MIDDLE_LEFT, 0);
   sprintf(tmp, "%d  %s", bsp->length, textid);
   AddLabel(picture, ((bsp->length+10)*680)/len, numsaps*10-3, tmp, 0, 0, MIDDLE_RIGHT, 0);
   BioseqUnlock(bsp);
   i = numsaps-1;
   while (salp != NULL)
   {
      acc = acc_head;
      found = FALSE;
      while (acc != NULL && !found)
      {
         if (acc->sap == salp)
            found = TRUE;
         else
            acc = acc->next;
      }
      AlnMgr2GetNthSeqRangeInSA(salp, 1, &start, &stop);
      start_r = (start*680)/len;
      stop_r = (stop*680)/len;
      AddRectangle(picture, start_r, i*10, stop_r, i*10-7, 0, TRUE, 0);
      dsp = (DenseSegPtr)(salp->segs);
      sprintf(tmp, "%d", start+1);
      AddLabel(picture, start_r-10, i*10-3, tmp, 0, 0, MIDDLE_LEFT, 0);
      sprintf(tmp, "%d  %s", stop+1, acc->accession);
      AddLabel(picture, stop_r+10, i*10-3, tmp, 0, 0, MIDDLE_RIGHT, 0);
      salp = salp->next;
      i--;
   }
   AttachPicture(v, picture, 0, 0, UPPER_LEFT, 1, 1, FrameVwr);
   Show(w);
}


static void PrintTPAHistErrors (LogInfoPtr lip, ValNodePtr errors)
{
  ValNodePtr vnp;

  if (lip == NULL || lip->fp == NULL || errors == NULL) return;

  for (vnp = errors; vnp != NULL; vnp = vnp->next)
  {
    fprintf (lip->fp, "%s\n", vnp->data.ptrvalue);
    lip->data_in_log = TRUE;
  }
  fprintf (lip->fp, "\n\n");
}


static void CKA_RunChecker(SeqEntryPtr sep)
{
   CKA_AccPtr   acc;
   CKA_AccPtr   acc_head;
   CKA_AccPtr   acc_head_next;
   CKA_AccPtr   acc_head_prev;
   CKA_AccPtr   acc_head_real;
   CKA_AccPtr   acc_head_tmp;
   BioseqPtr    bsp;
   Boolean      found;
   Int4         gi;
   SeqIdPtr     lastid;
   SeqAlignPtr  sap;
   SeqHistPtr   shp;
   CharPtr      textid;
   Char         textid2[42];
   LogInfoPtr   lip;
   ValNodePtr   err_list;
   CharPtr      err_msg;
   CharPtr      no_align_fmt = "Accession %s does not align to the bioseq %s.\n";

   if (sep == NULL)
   {
      Message(MSG_ERROR, "Null SeqEntry passed to CKA_RunChecker");
      return;
   }
   acc_head = NULL;
   SeqEntryExplore(sep, &acc_head, CKA_FindAllTpaDescr);
   lastid = NULL;
   if (acc_head == NULL)
   {
      Message(MSG_ERROR, "No Tpa features found in SeqEntry.");
      return;
   }

   lip = OpenLog ("TPA Alignment Assembly Problems");

   acc_head_real = acc_head;
   while (acc_head != NULL)
   {
      lastid = acc_head->sip_whole;
      acc_head_prev = acc_head;
      acc_head_tmp = acc_head->next;
      found = FALSE;
      while (!found && acc_head_tmp != NULL)
      {
         if (SeqIdComp(lastid, acc_head_tmp->sip_whole) != SIC_YES)
            found = TRUE;
         else
         {
            acc_head_prev = acc_head_tmp;
            acc_head_tmp = acc_head_tmp->next;
         }
      }
      acc_head_next = acc_head_prev->next;
      acc_head_prev->next = NULL;
      bsp = BioseqLockById(acc_head->sip_whole);
      if (ISA_na(bsp->mol))
      {
         sap = CKA_MakeAlign(bsp, acc_head, lip);
         acc = acc_head;
         while (acc != NULL && acc->sap == NULL)
         {
            acc = acc->next;
         }
         SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(acc_head->sip_whole), &gi, &textid);
         if (textid == NULL)
         {
            sprintf(textid2, "%d", gi);
            textid = textid2;
         }

         err_list = NULL;

         /* report each accession that does not align to the bioseq */
         acc = acc_head;
         while (acc != NULL) 
         {
           if (acc->sap == NULL) 
           {
             err_msg = (CharPtr) MemNew (sizeof(Char) * (StringLen (no_align_fmt) + StringLen (acc->accession) + StringLen (textid)));
             sprintf (err_msg, no_align_fmt, acc->accession, textid);
             ValNodeAddPointer (&err_list, 0, err_msg);
           }
           acc = acc->next;
         }

         if (sap != NULL) {
            AlnMgr2IndexLite(sap);
            AlnMgr2SortAlnSetByNthRowPos(sap, 1);
	          /* make seq-hist and add it to record */
	          if (bsp->hist != NULL)
	          {
	              shp = bsp->hist;
	              if (shp->assembly != NULL)
	                SeqAlignSetFree(shp->assembly);
	              shp->assembly = (SeqAlignPtr)(sap->segs);
	          }
            else
	          {
	              shp = SeqHistNew();
	              shp->assembly = (SeqAlignPtr)(sap->segs);
	              bsp->hist = shp;
	          }
         }
      
         if (ValidateTPAHistAlign(bsp, &err_list)) 
         {
            fprintf (lip->fp, "Alignments were successfully created and are being added to %s.\n", textid);
            lip->data_in_log = TRUE;
            PrintTPAHistErrors (lip, err_list);
            err_list = ValNodeFreeData (err_list);
         }
	       else if (sap != NULL)
         {
            fprintf (lip->fp, "Alignments were created but are not valid. They are being added to %s for review.\n", textid);
            lip->data_in_log = TRUE;
            PrintTPAHistErrors (lip, err_list);
            err_list = ValNodeFreeData (err_list);
         }
         else if (sap == NULL) 
         {
            fprintf (lip->fp, "No alignments could be created for %s.\n", textid);
            lip->data_in_log = TRUE;
            PrintTPAHistErrors (lip, err_list);
            err_list = ValNodeFreeData (err_list);
            acc_head = acc_head_next;
            BioseqUnlock (bsp);
            continue;
         } else {
            PrintTPAHistErrors (lip, err_list);
            err_list = ValNodeFreeData (err_list);
         }

         if (sap != NULL) {
            sap->segs = NULL;
            SeqAlignFree(sap);
         }
	     } 
     else
     {
        fprintf (lip->fp, "%s is annotated on a non-nucleotide bioseq.\n", acc_head->accession);
        lip->data_in_log = TRUE;
     }
     BioseqUnlock (bsp);
     acc_head = acc_head_next;
  }
  /*CKA_ShowAln(sap, acc_head_real);*/
  while (acc_head_real != NULL)
  {
     acc_head_tmp = acc_head_real->next;
     MemFree(acc_head_real->accession);
     SeqIdFree(acc_head_real->sip_whole);
     MemFree(acc_head_real);
     acc_head_real = acc_head_tmp;
  }
  CloseLog (lip);
  lip = FreeLog (lip);
}


typedef struct seqalignrow {
  Int4 start;
  Int4 stop;
  Uint1 strand;
} SeqAlignRowData, PNTR SeqAlignRowPtr;

static SeqAlignRowPtr SeqAlignRowNew (Int4 start, Int4 stop, Uint1 strand)
{
  SeqAlignRowPtr r;
  Int4 tmp;

  r = (SeqAlignRowPtr) MemNew (sizeof (SeqAlignRowData));
  r->start = start;
  r->stop = stop;

  if (r->start > r->stop) {
    tmp = r->start;
    r->start = r->stop;
    r->stop = tmp;
  }
  r->strand = strand;
  return r;
}


static SeqAlignRowPtr SeqAlignRowCopy (SeqAlignRowPtr orig)
{
  SeqAlignRowPtr r = NULL;

  if (orig != NULL) {
    r = SeqAlignRowNew (orig->start, orig->stop, orig->strand);
  }
  return r;
}


static SeqAlignRowPtr SeqAlignRowFree (SeqAlignRowPtr r)
{
  r = MemFree (r);
  return r;
}


static Int4 RowDiff (SeqAlignRowPtr r1, SeqAlignRowPtr r2)
{
  Int4 diff = 0;

  if (r1 == NULL || r2 == NULL) {
    return -1;
  }

  diff = ABS(r1->start - r2->start) + ABS (r1->stop - r2->stop);
  return diff;
}


static Int4 SeqAlignRowLen (SeqAlignRowPtr r) 
{
  Int4 len = 0;

  if (r != NULL) {
    len = r->stop - r->start + 1;
  }
  return len;
}


typedef struct seqalignsort {
  SeqAlignRowPtr row1;
  SeqAlignRowPtr row2;
  SeqAlignPtr salp;
} SeqAlignSortData, PNTR SeqAlignSortPtr;


static SeqAlignSortPtr SeqAlignSortNew (SeqAlignPtr salp)
{
  SeqAlignSortPtr s;

  if (salp == NULL) {
    return NULL;
  }

  s = (SeqAlignSortPtr) MemNew (sizeof (SeqAlignSortData));
  s->salp = salp;

  AlnMgr2IndexSingleChildSeqAlign(salp);

  s->row1 = SeqAlignRowNew (SeqAlignStart (salp, 0), SeqAlignStop (salp, 0), SeqAlignStrand (salp, 0));
  s->row2 = SeqAlignRowNew (SeqAlignStart (salp, 1), SeqAlignStop (salp, 1), SeqAlignStrand (salp, 1));

  return s;
}


static SeqAlignSortPtr SeqAlignSortFree (SeqAlignSortPtr s)
{
  if (s != NULL) {
    s->row1 = SeqAlignRowFree (s->row1);
    s->row2 = SeqAlignRowFree (s->row2);
    s = MemFree (s);
  }
  return s;
}


static ValNodePtr SeqAlignSortListNew (SeqAlignPtr salp)
{
  ValNodePtr list = NULL;
  SeqAlignPtr salp_next;

  while (salp != NULL) {
    salp_next = salp->next;
    salp->next = NULL;
    ValNodeAddPointer (&list, 0, SeqAlignSortNew (salp));
    salp = salp_next;
  }
  return list;
}


static ValNodePtr SeqAlignSortListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = SeqAlignSortFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static SeqAlignRowPtr SeqAlignRowFromSeqAlignSort (SeqAlignSortPtr s, Int4 row)
{
  if (s == NULL) {
    return NULL;
  } else if (row == 1) {
    return s->row1;
  } else {
    return s->row2;
  }
}


static Uint1 SeqAlignSortRowStrand (SeqAlignSortPtr s, Int4 row)
{
  Uint1 strand = Seq_strand_plus;
  SeqAlignRowPtr r;

  r = SeqAlignRowFromSeqAlignSort (s, row);
  if (r != NULL) {
    strand = r->strand;
  }
  return strand;
}


static Uint1 SeqAlignSortListFindBestStrand (ValNodePtr vnp, Int4 row)
{
  Int4 num_plus = 0, num_minus = 0;
  SeqAlignSortPtr s;

  while (vnp != NULL) {
    s = (SeqAlignSortPtr) vnp->data.ptrvalue;
    if (s != NULL) {
      if (SeqAlignSortRowStrand(s, row) == Seq_strand_minus) {
        num_minus++;
      } else {
        num_plus++;
      }
    }
    vnp = vnp->next;
  }

  if (num_minus > num_plus) {
    return Seq_strand_minus;
  } else {
    return Seq_strand_plus;
  }
}


static void SeqAlignSortListMarkStrand (ValNodePtr vnp, Int4 row, Uint1 strand)
{
  while (vnp != NULL) {
    if (SeqAlignSortRowStrand (vnp->data.ptrvalue, row) == strand) {
      vnp->choice = 1;
    }
    vnp = vnp->next;
  }
}


static ValNodePtr SeqAlignSortListRemoveAll (ValNodePtr list)
{
  ValNodePtr vnp;
  SeqAlignSortPtr s;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    s = vnp->data.ptrvalue;
    if (s != NULL && s->salp != NULL) {
      s->salp->next = NULL;
      s->salp = SeqAlignFree (s->salp);
    }
  }

  list = SeqAlignSortListFree (list);
  return list;
}


static void SeqAlignSortListRemoveMarked (ValNodePtr PNTR list)
{
  ValNodePtr remove_list;

  if (list == NULL) {
    return;
  }

  remove_list = ValNodeExtractList (list, 1);
  remove_list = SeqAlignSortListRemoveAll (remove_list);
}


static Uint1 SeqAlignSortListRemoveConflictingStrands (ValNodePtr PNTR list, Int4 row)
{
  Uint1 strand;

  if (list == NULL) {
    return Seq_strand_plus;
  }

  strand = SeqAlignSortListFindBestStrand (*list, row);
  if (strand == Seq_strand_plus) {
    SeqAlignSortListMarkStrand (*list, row, Seq_strand_minus);
  } else {
    SeqAlignSortListMarkStrand (*list, row, Seq_strand_plus);
  }

  SeqAlignSortListRemoveMarked (list);
  return strand;
}


static SeqAlignPtr SeqAlignFromSeqAlignSortList (ValNodePtr vnp)
{
  SeqAlignPtr salp_list = NULL, salp_prev = NULL;
  SeqAlignSortPtr s;

  while (vnp != NULL) {
    s = (SeqAlignSortPtr) vnp->data.ptrvalue;
    if (s != NULL && s->salp != NULL) {
      s->salp->next = NULL;
      if (salp_prev == NULL) {
        salp_list = s->salp;
      } else {
        salp_prev->next = s->salp;
      }
      salp_prev = s->salp;
      s->salp->next = NULL;
    }
    vnp = vnp->next;
  }
  return salp_list;
}


static int CompareSeqAlignRow (SeqAlignRowPtr r1, SeqAlignRowPtr r2)
{
  int rval = 0;

  if (r1 == NULL && r2 == NULL) {
    rval = 0;
  } else if (r1 == NULL) {
    rval = -1;
  } else if (r2 == NULL) {
    rval = 1;
  } else if (r1->start < r2->start) {
    rval = -1;
  } else if (r1->start > r2->start) {
    rval = 1;
  } else if (r1->stop < r2->stop) {
    rval = -1;
  } else if (r1->stop > r2->stop) {
    rval = 1;
  } else if (r1->strand < r2->strand) {
    rval = -1;
  } else if (r1->strand > r2->strand) {
    rval = 1;
  }
  return rval;
}


static int CompareSeqAlignSortPreferRow1 (SeqAlignSortPtr s1, SeqAlignSortPtr s2)
{
  int rval = 0;
  if (s1 == NULL && s2 == NULL) {
    rval = 0;
  } else if (s1 == NULL) {
    rval = -1;
  } else if (s2 == NULL) {
    rval = 1;
  } else if ((rval = CompareSeqAlignRow (s1->row1,s2->row1)) == 0) {
    rval = CompareSeqAlignRow (s1->row2, s2->row2);
  }
  return rval;
}


static int CompareSeqAlignSortPreferRow2 (SeqAlignSortPtr s1, SeqAlignSortPtr s2)
{
  int rval = 0;
  if (s1 == NULL && s2 == NULL) {
    rval = 0;
  } else if (s1 == NULL) {
    rval = -1;
  } else if (s2 == NULL) {
    rval = 1;
  } else if ((rval = CompareSeqAlignRow (s1->row2,s2->row2)) == 0) {
    rval = CompareSeqAlignRow (s1->row1, s2->row1);
  }
  return rval;
}



static int LIBCALLBACK SortVnpBySeqAlignSortRow1 (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return CompareSeqAlignSortPreferRow1 (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static int LIBCALLBACK SortVnpBySeqAlignSortRow2 (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return CompareSeqAlignSortPreferRow2 (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static ValNodePtr SeqAlignSortListExtractRepeats (ValNodePtr PNTR list, Int4 row, Int4 fuzz)
{
  ValNodePtr repeat_start, repeat_prev = NULL, vnp_prev = NULL, vnp;
  ValNodePtr repeat_list = NULL;
  SeqAlignSortPtr s1, s2;
  SeqAlignRowPtr  interval, r2;
  Int4            diff;
  Boolean         is_repeat;

  if (list == NULL || *list == NULL || (*list)->next == NULL) {
    return NULL;
  }

  if (row == 1) {
    *list = ValNodeSort (*list, SortVnpBySeqAlignSortRow1);
  } else {
    *list = ValNodeSort (*list, SortVnpBySeqAlignSortRow2);
  }

  repeat_start = *list;
  s1 = repeat_start->data.ptrvalue;
  interval = SeqAlignRowCopy (SeqAlignRowFromSeqAlignSort(s1, row));
  vnp_prev = *list;
  for (vnp = (*list)->next; vnp != NULL; vnp = vnp->next) {
    s2 = vnp->data.ptrvalue;
    is_repeat = FALSE;
    r2 = SeqAlignRowFromSeqAlignSort (s2, row);

    if (interval->start <= r2->start && interval->stop >= r2->stop) {
      /* contained */
      is_repeat = TRUE;
    } else if (r2->start <= interval->start && r2->stop >= interval->stop) {
      /* contained */
      is_repeat = TRUE;
    } else if ((diff = RowDiff (interval, r2)) > -1 && diff < fuzz) {
      is_repeat = TRUE;
    }
    if (is_repeat) {
      vnp_prev = vnp;
      if (interval->start > r2->start) {
        interval->start = r2->start;
      }
      if (interval->stop < r2->stop) {
        interval->stop = r2->stop;
      }
    } else {
      if (repeat_start->next == vnp) {
        repeat_prev = vnp_prev;
        vnp_prev = vnp;
      } else {
        if (repeat_prev == NULL) {
          *list = vnp;
        } else {
          repeat_prev->next = vnp;
        }
        if (vnp_prev != NULL) {
          vnp_prev->next = NULL;
        }
        ValNodeAddPointer (&repeat_list, 0, repeat_start);
        vnp_prev = vnp;
      }
      repeat_start = vnp;
      s1 = vnp->data.ptrvalue;
      interval = SeqAlignRowFree (interval);
      interval = SeqAlignRowCopy (SeqAlignRowFromSeqAlignSort(s1, row));
    }
  }

  if (repeat_start->next != NULL) {
    if (repeat_prev == NULL) {
      *list = NULL;
    } else {
      repeat_prev->next = NULL;
    }
    ValNodeAddPointer (&repeat_list, 0, repeat_start);
  }

  interval = SeqAlignRowFree (interval);
  return repeat_list;
}


static int SeqAlignRowFuzzyCompare (SeqAlignRowPtr r1, SeqAlignRowPtr r2, Int4 fuzz)
{
  if (r1 == NULL && r2 == NULL) {
    return 0;
  } else if (r1 == NULL) {
    return -1;
  } else if (r2 == NULL) {
    return 1;
  }

  if (r1->stop < r2->start || r1->stop - r2->start < fuzz) {
    return -1;
  } else if (r2->stop < r1->start || r2->stop - r1->start < fuzz) {
    return 1;
  } else {
    return 0;
  }
}


static int SeqAlignSortFuzzyCompare (SeqAlignSortPtr s1, SeqAlignSortPtr s2, Int4 row, Int4 fuzz)
{
  if (s1 == NULL && s2 == NULL) {
    return 0;
  } else if (s1 == NULL) {
    return -1;
  } else if (s2 == NULL) {
    return 1;
  } else if (row == 1) {
    return SeqAlignRowFuzzyCompare (s1->row1, s2->row1, fuzz);
  } else {
    return SeqAlignRowFuzzyCompare (s1->row2, s2->row2, fuzz);
  }
}


static void SeqAlignSortListRemoveIntervalsOutOfOrder (ValNodePtr PNTR list, Int4 row, Int4 fuzz)
{
  ValNodePtr vnp, vnp_prev = NULL;

  if (list == NULL) {
    return;
  }

  for (vnp = *list; vnp != NULL; vnp = vnp->next) {
    if (vnp_prev != NULL && SeqAlignSortFuzzyCompare (vnp_prev->data.ptrvalue, vnp->data.ptrvalue, row, fuzz) != -1) {
      vnp->choice = 1;
    } else if (vnp->next != NULL && SeqAlignSortFuzzyCompare (vnp->data.ptrvalue, vnp->next->data.ptrvalue, row, fuzz) != -1) {
      if (vnp->next->next != NULL 
          && SeqAlignSortFuzzyCompare (vnp->data.ptrvalue, vnp->next->next->data.ptrvalue, row, fuzz) == -1
          && SeqAlignSortFuzzyCompare (vnp->next->data.ptrvalue, vnp->next->next->data.ptrvalue, row, fuzz) != -1) {
        /* ok to keep this one, we'll toss the next one */
      } else {
        vnp->choice = 1;
      }
    }
    if (vnp->choice == 0) {
      vnp_prev = vnp;
    }
  }

  SeqAlignSortListRemoveMarked (list);
}


static Int4 GetRepeatIntervalFuzz (SeqAlignSortPtr s_repeat, SeqAlignSortPtr s_before, SeqAlignSortPtr s_after, Int4 row)
{
  Int4 start_fuzz = 0, end_fuzz = 0;

  if (s_repeat == NULL) {
    return -1;
  }

  if (s_before != NULL) {
    if (row == 1) {
      start_fuzz = ABS (s_repeat->row1->start - s_before->row1->stop);
    } else {
      start_fuzz = ABS (s_repeat->row2->start - s_before->row2->stop);
    }
  }
  if (s_after != NULL) {
    if (row == 1) {
      end_fuzz = ABS (s_after->row1->start - s_repeat->row1->stop);
    } else {
      end_fuzz = ABS (s_after->row2->start - s_repeat->row2->stop);
    }
  }

  return start_fuzz + end_fuzz;
}


static int StrandedSeqAlignSortRowCompare (SeqAlignSortPtr s1, SeqAlignSortPtr s2, Int4 row, Int4 fuzz)
{
  SeqAlignRowPtr r1 = NULL, r2 = NULL;
  int rval = 0;
  Uint1 strand = Seq_strand_plus;

  if (s1 == NULL || s2 == NULL) {
    return 0;
  } 

  r1 = SeqAlignRowFromSeqAlignSort (s1, row);
  r2 = SeqAlignRowFromSeqAlignSort (s2, row);
  strand = r1->strand;

  if (strand == Seq_strand_minus) {
    if (r1->start < r2->start - fuzz) {
      rval = 1;
    } else if (r1->start >r2->start + fuzz) {
      rval = -1;
    }
  } else {
    if (r1->start > r2->start + fuzz) {
      rval = 1;
    } else if (r1->stop < r2->stop - fuzz) {
      rval = -1;
    } 
  }

  return rval;
}


static Boolean FindSeqAlignSortWithPoint (ValNodePtr list, Int4 point, Int4 row)
{
  ValNodePtr vnp;
  SeqAlignRowPtr r;
  SeqAlignSortPtr s;
  Boolean found = FALSE;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    s = vnp->data.ptrvalue;
    r = SeqAlignRowFromSeqAlignSort (s, row);
    if (point >= r->start && point <= r->stop) {
      found = TRUE;
    }
  }
  return found;
}


static ValNodePtr FindBestRepeat (ValNodePtr PNTR repeat_list, SeqAlignSortPtr s_before, SeqAlignSortPtr s_after, Int4 row, Int4 fuzz)
{
  Int4       best_diff = -1, diff;
  ValNodePtr vnp, vnp_best = NULL, vnp_best_prev = NULL, vnp_prev = NULL;
  SeqAlignSortPtr s_this;

  if (repeat_list == NULL || *repeat_list == NULL) {
    return NULL;
  }

  for (vnp = *repeat_list; vnp != NULL; vnp = vnp->next) {
    s_this = vnp->data.ptrvalue;

    if (StrandedSeqAlignSortRowCompare (s_this, s_before, row, fuzz) < 0
      || StrandedSeqAlignSortRowCompare (s_this, s_after, row, fuzz) > 0) {
      /* skip - already out of order */
    } else {
      diff = GetRepeatIntervalFuzz (vnp->data.ptrvalue, s_before, s_after, row);
      if (diff > -1 && (best_diff < 0 || best_diff > diff)) {
        vnp_best = vnp;
        vnp_best_prev = vnp_prev;
        best_diff = diff;
      }
    }
    vnp_prev = vnp;
  }

  if (vnp_best != NULL) {
    if (vnp_best_prev == NULL) {
      *repeat_list = vnp_best->next;
    } else {
      vnp_best_prev->next = vnp_best->next;
    }
    vnp_best->next = NULL;
  }

  return vnp_best;
}


static void RemoveRepeatsCoincidingWithBest (ValNodePtr PNTR repeat_list, ValNodePtr best, Int4 row, Int4 fuzz)
{
  ValNodePtr vnp;
  SeqAlignSortPtr s_best, s;
  SeqAlignRowPtr r_best, r2;
  Boolean        is_repeat;
  Int4           diff;

  if (repeat_list == NULL || *repeat_list == NULL || best == NULL) {
    return;
  }

  s_best = best->data.ptrvalue;
  r_best = SeqAlignRowFromSeqAlignSort (s_best, row);

  for (vnp = *repeat_list; vnp != NULL; vnp = vnp->next) {
    s = vnp->data.ptrvalue;
    r2 = SeqAlignRowFromSeqAlignSort (s, row);

    is_repeat = FALSE;
    if (r_best->start <= r2->start && r_best->stop >= r2->stop) {
      /* contained */
      is_repeat = TRUE;
    } else if (r2->start <= r_best->start && r2->stop >= r_best->stop) {
      /* contained */
      is_repeat = TRUE;
    } else if (r2->stop < r_best->stop || r2->stop - r_best->stop < fuzz) {
      is_repeat = TRUE;
    } else if ((diff = RowDiff (r_best, r2)) > -1 && diff < fuzz) {
      is_repeat = TRUE;
    }

    if (is_repeat) {
      vnp->choice = 1;
    }
  }
  SeqAlignSortListRemoveMarked (repeat_list);
}


static ValNodePtr ExtractLongestSeqAlignRow (ValNodePtr PNTR list, Int4 row)
{
  ValNodePtr vnp, vnp_prev = NULL, rval = NULL;
  Int4       longest = 0, len;
  Boolean    found = FALSE;

  if (list == NULL || *list == NULL) {
    return NULL;
  }

  for (vnp = *list; vnp != NULL; vnp = vnp->next) {
    len = SeqAlignRowLen (SeqAlignRowFromSeqAlignSort (vnp->data.ptrvalue, row));
    if (longest < len) {
      longest = len;
    }
  }

  for (vnp = *list; vnp != NULL && !found; vnp = vnp->next) {
    len = SeqAlignRowLen (SeqAlignRowFromSeqAlignSort (vnp->data.ptrvalue, row));
    if (len == longest) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      rval = vnp;
      found = TRUE;
    }
    vnp_prev = vnp;
  }
  return rval;
}


static void InsertBestRepeat (ValNodePtr repeat_list, ValNodePtr PNTR sorted_list, Int4 row, Int4 fuzz)
{
  ValNodePtr vnp, vnp_prev = NULL, vnp_new;
  SeqAlignSortPtr s_repeat, s, s_before = NULL;
  Boolean         found = TRUE;
  SeqAlignRowPtr  r1, r2;
  Int4 other_row;

  if (repeat_list == NULL || sorted_list == NULL) {
    return;
  }

  if (row == 1) {
    other_row = 2;
  } else {
    other_row = 1;
  }
  if (sorted_list == NULL || *sorted_list == NULL) {
    /* keep longest, mark others for removal */
    vnp = ExtractLongestSeqAlignRow (&repeat_list, row);
    ValNodeLink (sorted_list,vnp);
    repeat_list = SeqAlignSortListRemoveAll (repeat_list);
  } else {
    s_repeat = repeat_list->data.ptrvalue;
    found = FALSE;
    vnp = *sorted_list;

    /* find first entry that is after this repeat, and insert before that */
    while (vnp != NULL && !found) {
      s = vnp->data.ptrvalue;
      r1 = SeqAlignRowFromSeqAlignSort (s, row);
      r2 = SeqAlignRowFromSeqAlignSort (s_repeat, row);

      if (r1->start > r2->start || r2->start - r1->start < fuzz) {
        while (repeat_list != NULL) {
          /* extract best repeat */
          vnp_new = FindBestRepeat (&repeat_list, s_before, vnp->data.ptrvalue, other_row, fuzz);
          if (vnp_new == NULL) {
            repeat_list = SeqAlignSortListRemoveAll (repeat_list);
          } else {
            RemoveRepeatsCoincidingWithBest (&repeat_list, vnp_new, row, fuzz);
            vnp_new->next = vnp;
            if (vnp_prev == NULL) {
              *sorted_list = vnp_new;
            } else {
              vnp_prev->next = vnp_new;
            }
            vnp_prev = vnp_new;
            s_before = vnp_new->data.ptrvalue;
          }
        }
        found = TRUE;
      }
      if (!found) {
        s_before = vnp->data.ptrvalue;
        vnp_prev = vnp;
        vnp = vnp->next;
      }
    }
    if (!found) {
      while (repeat_list != NULL) {
        /* extract best repeat */
        vnp_new = FindBestRepeat (&repeat_list, s_before, NULL, other_row, fuzz);
        if (vnp_new == NULL) {
          repeat_list = SeqAlignSortListRemoveAll (repeat_list);
        } else {
          RemoveRepeatsCoincidingWithBest (&repeat_list, vnp_new, row, fuzz);
          vnp_new->next = NULL;
          if (vnp_prev == NULL) {
            *sorted_list = vnp_new;
          } else {
            vnp_prev->next = vnp_new;
          }
          vnp_prev = vnp_new;
          s_before = vnp_new->data.ptrvalue;
        }
      }
    }
  }

  SeqAlignSortListRemoveMarked (&repeat_list);
  repeat_list = ValNodeFree (repeat_list);
}


static void SelectBestRepeatsFromList (SeqAlignPtr PNTR salp)
{
  ValNodePtr list, vnp;
  ValNodePtr row1_repeats, row2_repeats;
  Uint1 strand1, strand2;
  SeqAlignPtr    tmp_salp;
  Int4           fuzz = 15;

  Int4           missing = 600;

  if (salp == NULL || *salp == NULL || (*salp)->next == NULL) {
    return;
  }

  list = SeqAlignSortListNew (*salp);

  FindSeqAlignSortWithPoint (list, missing, 1);

  /* remove conflicting strands for row 1 */
  strand1 = SeqAlignSortListRemoveConflictingStrands (&list, 1);

  /* remove conflicting strands for row 1 */
  strand2 = SeqAlignSortListRemoveConflictingStrands (&list, 2);

  FindSeqAlignSortWithPoint (list, missing, 1);

  if (list != NULL && list->next != NULL) {
    row1_repeats = SeqAlignSortListExtractRepeats (&list, 1, fuzz);
    row2_repeats = SeqAlignSortListExtractRepeats (&list, 2, fuzz);

    FindSeqAlignSortWithPoint (list, missing, 1);

    /* remove scaffold intervals that are out of order */
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow1);
    SeqAlignSortListRemoveIntervalsOutOfOrder (&list, 1, fuzz);
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow2);
    SeqAlignSortListRemoveIntervalsOutOfOrder (&list, 2, fuzz);

    FindSeqAlignSortWithPoint (list, missing, 1);

    /* Remove overlaps.*/
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow1);
    tmp_salp = SeqAlignFromSeqAlignSortList (list);
    list = SeqAlignSortListFree (list);
    ACT_RemoveInconsistentAlnsFromSet (tmp_salp, 1, 1);
    list = SeqAlignSortListNew (tmp_salp);

    FindSeqAlignSortWithPoint (list, missing, 1);

    /* for each repeat on row 1, we want to pick the most consistent interval for row 2 */
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow1);
    for (vnp = row1_repeats; vnp != NULL; vnp = vnp->next) {   
      InsertBestRepeat (vnp->data.ptrvalue, &list, 1, fuzz);
    }
    row1_repeats = ValNodeFree (row1_repeats);

    /* for each repeat on row 2, we want to pick the most consistent interval for row 1 */
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow2);
    for (vnp = row2_repeats; vnp != NULL; vnp = vnp->next) {   
      InsertBestRepeat (vnp->data.ptrvalue, &list, 2, fuzz);
    }
    row2_repeats = ValNodeFree (row2_repeats);
  }

  list = ValNodeSort (list, SortVnpBySeqAlignSortRow1);
  *salp = SeqAlignFromSeqAlignSortList (list);

  list = SeqAlignSortListFree (list);
}


static void amconssetfree(AMConsSetPtr acp)
{
   AMConsSetPtr  acp_next;

   while (acp != NULL)
   {
      acp_next = acp->next;
      MemFree(acp->starts);
      MemFree(acp->stops);
      MemFree(acp->strands);
      MemFree(acp);
      acp = acp_next;
   }
}

static int LIBCALLBACK CKA_SortForConsistent(VoidPtr ptr1, VoidPtr ptr2)
{
   AMConsSetPtr  acp1;
   AMConsSetPtr  acp2;
   FloatHi       bitscore;
   FloatHi       evalue;
   Int4          number;
   SAIndex2Ptr   saip1;
   SAIndex2Ptr   saip2;

   acp1 = *((AMConsSetPtr PNTR)ptr1);
   acp2 = *((AMConsSetPtr PNTR)ptr2);
   saip1 = (SAIndex2Ptr)(acp1->sap->saip);
   saip2 = (SAIndex2Ptr)(acp2->sap->saip);
   if (saip1->score == 0)
      GetScoreAndEvalue(acp1->sap, &saip1->score, &bitscore, &evalue, &number);
   if (saip2->score == 0)
      GetScoreAndEvalue(acp2->sap, &saip2->score, &bitscore, &evalue, &number);
   if (saip1->score > saip2->score)
      return -1;
   else if (saip1->score < saip2->score)
      return 1;
   else
      return 0;
}


static void CKA_RemoveInconsistentAlnsFromSet(SeqAlignPtr sap_head, Int4 fuzz)
{
   AMConsSetPtr  acp;
   AMConsSetPtr  acp_head;
   AMConsSetPtr  acp_prev;
   AMConsSetPtr  PNTR acparray;
   DenseSegPtr   dsp;
   Int4          i;
   Int4          j;
   Int4          k;
   Int4          lfuzz;
   SeqAlignPtr   newsap;
   Int4          numrows;
   Int4          numsaps;
   Int4          orientation;
   Int4          row;
   SAIndex2Ptr   saip;
   SeqAlignPtr   salp_head;
   SeqAlignPtr   salp_prev;
   SeqAlignPtr   sap;
   SeqAlignPtr   sapnext;
   Int4          score;
   SeqIdPtr      sip;
   SeqIdPtr      sip_head;
   Uint1         strand;

   lfuzz = fuzz;
   if (fuzz < 0)
      fuzz = 1;
   sap = (SeqAlignPtr)(sap_head->segs);
   if (sap->next == NULL)
      return;
   dsp = (DenseSegPtr)(sap->segs);
   sip_head = dsp->ids;
   numrows = AlnMgr2GetNumRows(sap);
   acp_head = NULL;
   strand = AlnMgr2GetNthStrand(sap, 1);
   numsaps = 0;
   while (sap != NULL)
   {
      if (AlnMgr2GetNumRows(sap) != numrows)
      {
         amconssetfree(acp_head);
         return;
      }
      numsaps++;
      acp = (AMConsSetPtr)MemNew(sizeof(AMConsSet));
      acp->starts = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->stops = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->strands = (Uint1Ptr)MemNew(numrows*sizeof(Uint1));
      acp->which = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->sap = sap;
      if (acp_head != NULL)
      {
         acp_prev->next = acp;
         acp_prev = acp;
      } else
         acp_head = acp_prev = acp;
      sip = sip_head;
      row = AlnMgr2GetFirstNForSip(sap, sip);
      if (row <= 0)
      {
         amconssetfree(acp_head);
         return;
      }
      if (acp->strands[row] != strand)
      {
         sapnext = acp->sap->next;
         acp->sap->next = NULL;
         score = ((SAIndex2Ptr)(acp->sap->saip))->score;
         SeqAlignListReverseStrand(acp->sap);
         AMAlignIndexFreeEitherIndex(acp->sap);
         AlnMgr2IndexSingleChildSeqAlign(acp->sap);
         saip = (SAIndex2Ptr)(acp->sap->saip);
         saip->score = score;
         acp->strands[row] = strand;
         acp->sap->next = sapnext;
      }
      for (i=0; i<numrows; i++)
      {
         acp->which[i] = row;
         AlnMgr2GetNthSeqRangeInSA(sap, i+1, &acp->starts[i], &acp->stops[i]);
         acp->strands[i] = AlnMgr2GetNthStrand(sap, i+1);
      }
      sap = sap->next;
   }
   acparray = (AMConsSetPtr PNTR)MemNew(numsaps*sizeof(AMConsSetPtr));
   acp = acp_head;
   i = 0;
   while (acp != NULL)
   {
      acparray[i] = acp;
      acp = acp->next;
      i++;
   }
   HeapSort(acparray, numsaps, sizeof(AMConsSetPtr), CKA_SortForConsistent);
   /* orientation -1 means that ith is before jth in ALL rows, 1 means ith is after jth in ALL rows */
   for (i=0; i<numsaps; i++)
   {
      if (acparray[i]->used != -1)
      {
         for (j=i+1; j<numsaps; j++)
         {
            orientation = 0;
            for (k=0; acparray[j]->used != -1 && k<numrows; k++)
            {
               if (acparray[i]->starts[k] - fuzz < acparray[j]->starts[k])
               {
                  if (acparray[i]->stops[k] - fuzz < acparray[j]->starts[k])
                  {
                     if (orientation == 0)
                     {
                        if (acparray[i]->strands[k] == Seq_strand_minus)
                           orientation = 1;
                        else
                           orientation = -1;
                     }
                  } else
                  {
                     if (lfuzz >= 0) /* just mark it for deletion */
                        acparray[j]->used = -1;
                     else /* truncate it */
                     {
                        if (acparray[j]->stops[k] >
                            acparray[i]->stops[k] + CKA_blast_wordsize)
                        {
                           newsap = AlnMgr2GetSubAlign(acparray[j]->sap, acparray[i]->stops[k]+1,
 acparray[j]->stops[k], k+1, TRUE);
                           AlnMgr2IndexSingleChildSeqAlign(newsap);
                           SeqAlignFree(acparray[j]->sap);
                           acparray[j]->sap = newsap;
                           acparray[j]->starts[k] = acparray[i]->stops[k]+1;
                        } else
                           acparray[j]->used = -1;
                     }
                  }
               } else if (acparray[i]->starts[k] - fuzz > acparray[j]->starts[k])
               {
                 if (acparray[i]->starts[k] + fuzz > acparray[j]->stops[k])
                  {
                     if (orientation == 0)
                     {
                        if (acparray[i]->strands[k] == Seq_strand_minus)
                           orientation = -1;
                        else
                           orientation = 1;
                     }
                  } else
                  {
                     if (lfuzz >= 0) /* mark for deletion */
                        acparray[j]->used = -1;
                     else /* truncate */
                     {
                        if (acparray[j]->starts[k] <
                            acparray[i]->starts[k] - CKA_blast_wordsize)
                        {
                           newsap = AlnMgr2GetSubAlign(acparray[j]->sap, acparray[j]->starts[k], acparray[i]->starts[k]-1, k+1, TRUE);
                           AlnMgr2IndexSingleChildSeqAlign(newsap);
                           SeqAlignFree(acparray[j]->sap);
                           acparray[j]->sap = newsap;
                           acparray[j]->stops[k] = acparray[i]->starts[k]-1;
                        } else
                           acparray[j]->used = -1;
                     }
                  }
               } else
                  acparray[j]->used = -1;
            }
         }
      }
   }
   /* now free all the unused ones, stick the rest back together, reindex, and return */
   salp_head = salp_prev = NULL;
   for (i=0; i<numsaps; i++)
   {
      if (acparray[i]->used == -1)
      {
         SeqAlignFree(acparray[i]->sap);
         acparray[i]->sap = NULL;
      } else
      {
         if (salp_head != NULL)
         {
            salp_prev->next = acparray[i]->sap;
            salp_prev = acparray[i]->sap;
            salp_prev->next = NULL;
         } else
         {
            salp_head = salp_prev = acparray[i]->sap;
            salp_prev->next = NULL;
         }
      }
   }
   amconssetfree(acp_head);
   MemFree(acparray);
   sap_head->segs = (Pointer)(salp_head);
   AMAlignIndex2Free2(sap_head->saip);
   AlnMgr2IndexLite(sap_head);
}


static BioseqPtr ReadFromTraceDb (CharPtr number)

{
  BioseqPtr    bsp = NULL;
  CONN         conn;
  time_t       currtime, starttime;
  FILE         *fp;
  time_t       max = 0;
  size_t       n_written;
  Char         path [PATH_MAX];
  Char         query [64];
  SeqEntryPtr  sep = NULL;
  EIO_Status   status;
  STimeout     timeout;
  long int     val;

  if (StringHasNoText (number)) return NULL;
  if (sscanf (number, "%ld", &val) != 1) return NULL;
  sprintf (query, "cmd=raw&query=retrieve+fasta+%ld", (long) val);
  conn = QUERY_OpenUrlQuery ("www.ncbi.nlm.nih.gov", 80, "/Traces/trace.cgi",
                             query, "Sequin", 30, eMIME_T_NcbiData,
                             eMIME_Fasta, eENCOD_None, 0);
  if (conn == NULL) return NULL;
  status = CONN_Write (conn, (const void *) query, StringLen (query),
                       &n_written, eIO_WritePersist);
  if (status != eIO_Success) return NULL;
  QUERY_SendQuery (conn);

#ifdef OS_MAC 
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 100;
  timeout.usec = 0;
#endif

  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) != eIO_Success && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
  }

  if (status == eIO_Success) {
    TmpNam (path);
    fp = FileOpen (path, "w");
    QUERY_CopyResultsToFile (conn, fp);
    FileClose (fp);
    /*
    LaunchGeneralTextViewer (path, "QueueFastaQueryToURL results");
    */
    fp = FileOpen (path, "r");
    sep = FastaToSeqEntry (fp, TRUE);
    FileClose (fp);
    FileRemove (path);
    if (sep != NULL) {
      bsp = FindNucBioseq (sep);
    }
  }
  CONN_Close (conn);

  return bsp;
}

static SeqAlignPtr GetNewBlastTPAHistAlignPiece (BioseqPtr bsp1, BioseqPtr bsp2)
{
   BLAST_SummaryOptions *options = NULL;
   SeqAlignPtr          salp = NULL;

   BLAST_SummaryOptionsInit(&options);
   options->program = eBlastn;
   options->use_megablast = TRUE;
   options->word_size = CKA_blast_wordsize;
   options->cutoff_evalue = CKA_blast_expect_value;
   options->hint = eNone;
   options->gap_x_dropoff = 30;
   options->gap_open = -1;
   options->gap_extend = -1;
   options->filter_string = StringSave ("F");

   BLAST_TwoSequencesSearch(options, bsp1, bsp2, &salp);
   BLAST_SummaryOptionsFree(options);

   return salp;
}


static SeqAlignPtr GetOldBlastTPAHistAlignPiece (BioseqPtr bsp1, BioseqPtr bsp2)
{
  BLAST_OptionsBlkPtr  options;
  SeqAlignPtr          salp;

  options = BLASTOptionNew("blastn", TRUE);
  options->is_megablast_search = TRUE;
  options->gap_open = options->gap_extend = 0;
  options->wordsize = CKA_blast_wordsize;
  options->expect_value = CKA_blast_expect_value;

  salp = BlastTwoSequences(bsp1, bsp2, "blastn", options);
  BLASTOptionDelete(options);

  return salp;
}


static Boolean IsHUPIDAccession (BioseqPtr bsp)
{
  Int4 j, num;
  ObjMgrDataPtr omdp;
  OMUserDataPtr omudp;
  ObjMgrProcPtr ompp = NULL;
  ObjMgrPtr omp;
  ObjMgrDataPtr PNTR omdpp;
  Boolean           rval = FALSE;

  ompp = NULL;
  omp = ObjMgrReadLock();

  omdpp = omp->datalist;
  if (omdpp != NULL) {
    num = omp->currobj;
    for (j = 0; j < num && ompp == NULL; j++) {
      if (omdpp[j] != NULL && omdpp[j]->datatype == OBJ_BIOSEQ 
          && omdpp[j]->dataptr == bsp) {

        omdp = ObjMgrFindTop (omp, omdpp[j]);
        if (omdp != NULL) {
          for (omudp = omdp->userdata; omudp != NULL && ompp == NULL; omudp = omudp->next)
          {
            if (omudp->proctype == OMPROC_FETCH)  /* caching function */
            {
              ompp = ObjMgrProcFind(omp, omudp->procid, NULL, 0);
            }
          }
        }
      }
    }
  }

  if (ompp != NULL && StringCmp (ompp->procname, "HUPBioseqFetch") == 0) 
  {
    rval = TRUE;
  }

  ObjMgrUnlock();
  return rval;
}


static SeqAlignPtr CKA_MakeAlign(BioseqPtr bsp, CKA_AccPtr acc_head, LogInfoPtr lip)
{
   CKA_AccPtr           acc;
   CKA_AccPtr           acc_new;
   CKA_AccPtr           acc_new_head_head;
   CKA_AccPtr           acc_new_head;
   CKA_AccPtr           acc_new_prev;
   SeqAlignPtr          allsap;
   SeqAlignPtr          allsap_prev;
   AMAlignIndex2Ptr     amaip;
   BioseqPtr            bsp_tmp = NULL;
   Int4                 i;
   BLAST_SummaryOptions *options = NULL;
   SBlastSeqalignArray * seqalign_arr=NULL;
   SeqAlignPtr          sap_new;
   SeqAlignPtr          sap_tmp;
   SeqAlignPtr          sap_tmp_next;
   SeqIdPtr             sip;
   Uint1                strand;
   Boolean              need_to_unlock = FALSE;

   if (bsp == NULL || acc_head == NULL)
      return NULL;
   acc = acc_head;
   allsap = NULL;
   allsap_prev = NULL;
   acc_new_head_head = acc_new_prev = NULL;
   BLAST_SummaryOptionsInit(&options);
   options->program = eBlastn;
   options->use_megablast = TRUE;
   options->word_size = CKA_blast_wordsize;
   options->cutoff_evalue = CKA_blast_expect_value;
   options->hint = eNone;
   options->gap_open = -1;
   options->gap_extend = -1;
   options->filter_string = StringSave ("F");
   while (acc != NULL)
   {

      if (need_to_unlock)
      {
        BioseqUnlock (bsp_tmp);
        need_to_unlock = FALSE;
      }

      bsp_tmp = NULL;
      if (StringNICmp (acc->accession, "ti", 2) == 0) {
        bsp_tmp = ReadFromTraceDb (acc->accession + 2);
      } else {
        sip = SeqIdFromAccessionDotVersion(acc->accession);
        bsp_tmp = BioseqLockById(sip);
        if (bsp_tmp != NULL) {
          need_to_unlock = TRUE;
          if (lip != NULL && lip->fp != NULL && IsHUPIDAccession (bsp_tmp)) {
            fprintf (lip->fp, "%s is unreleased accession in TPA Seq-Hist.\n", acc->accession);
            lip->data_in_log = TRUE;
          }
        }
      }
      if (bsp_tmp == NULL) {
        fprintf (lip->fp, "Unable to load %s", acc->accession);
        lip->data_in_log = TRUE;
        break;
      }
      if (bsp_tmp->id->next) {
        /* find the best accession */
        SeqIdPtr sip = SeqIdDup(SqnSeqIdFindBestAccession(bsp_tmp->id));
        bsp_tmp->id = SeqIdSetFree(bsp_tmp->id);
        bsp_tmp->id = sip;
      }
      if (!ISA_na(bsp_tmp->mol))
      {
         if (need_to_unlock)
         {
           BioseqUnlock(bsp_tmp);
           need_to_unlock = FALSE;
         }
         Message(MSG_ERROR, "%s is not a nucleotide bioseq.", acc->accession);
         break;
      }
      WatchCursor();
      if (acc->start_acc >=0 && acc->stop_acc >=0 &&
          acc->start_acc < bsp_tmp->length &&
          acc->start_acc < bsp_tmp->length) {
        SeqLocPtr slp1, slp2;
        if (acc->start_acc <= acc->stop_acc) {
          slp1 = SeqLocIntNew
            (acc->start_acc, acc->stop_acc, Seq_strand_plus, bsp_tmp->id);
        } else {
          slp1 = SeqLocIntNew
            (acc->stop_acc, acc->start_acc, Seq_strand_minus, bsp_tmp->id);
        }
        slp2 = SeqLocIntNew(0, bsp->length-1, Seq_strand_plus, bsp->id);
        acc->sap = NULL;
        seqalign_arr = NULL;
        BLAST_TwoSeqLocSets (options, slp1, slp2, NULL, &seqalign_arr, NULL, NULL, NULL);
        if (seqalign_arr != NULL)
        {
          acc->sap = seqalign_arr->array[0];
          seqalign_arr->array[0] = NULL;
          seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
        }
        SeqLocFree(slp1);
        SeqLocFree(slp2);
      } else {
        acc->sap = NULL;
        BLAST_TwoSequencesSearch(options, bsp_tmp, bsp, &(acc->sap));
      }
      ArrowCursor();
      acc->start_acc = acc->stop_acc = 0; /* reset, for later usage */
      if (acc->sap != NULL)
         SPI_flip_sa_list(acc->sap);
      acc_new_head = NULL;
      if (acc->sap != NULL && acc->sap->next != NULL)
      {
         if (!CKA_blast_allow_repeats) {
           SelectBestRepeatsFromList (&(acc->sap));
         }
         AlnMgr2IndexLite(acc->sap);
         if (!CKA_blast_allow_repeats) {
           CKA_RemoveInconsistentAlnsFromSet(acc->sap, -1);
         }
         sap_tmp = acc->sap;
         acc->sap = (SeqAlignPtr)(acc->sap->segs);
         sap_tmp->segs = NULL;
         SeqAlignFree(sap_tmp);
         sap_tmp = acc->sap->next;
         acc->sap->next = NULL;
         while (sap_tmp != NULL)
         {
            AlnMgr2IndexSingleChildSeqAlign(sap_tmp);
            sap_tmp_next = sap_tmp->next;
            sap_tmp->next = NULL;
            acc_new = (CKA_AccPtr)MemNew(sizeof(CKA_Acc));
            acc_new->accession = StringSave(acc->accession);
            acc_new->sip_whole = SeqIdDup(acc->sip_whole);
            acc_new->sap = sap_tmp;
            sap_tmp = sap_tmp_next;
            if (!acc_new_head) {
              acc_new_head = acc_new;
            }
            if (acc_new_prev) {
              acc_new_prev->next = acc_new;
            }
            acc_new_prev = acc_new;
         }
      } else if (acc->sap != NULL)
         AlnMgr2IndexSingleChildSeqAlign(acc->sap);
      if (acc->sap != NULL)
      {
         strand = AlnMgr2GetNthStrand(acc->sap, 1);
         if (strand == Seq_strand_minus)
         {
            SeqAlignListReverseStrand(acc->sap);
            SAIndex2Free2(acc->sap->saip);
            acc->sap->saip = NULL;
         }
      }
      if (allsap != NULL && acc->sap != NULL)
      {
         allsap_prev->next = acc->sap;
         allsap_prev = allsap_prev->next;;
      } else if (acc->sap != NULL)
         allsap_prev = allsap = (acc->sap);
      acc_new = acc_new_head;
      while (acc_new != NULL)
      {
         strand = AlnMgr2GetNthStrand(acc_new->sap, 1);
         if (strand == Seq_strand_minus)
         {
            SeqAlignListReverseStrand(acc_new->sap);
            SAIndex2Free2(acc_new->sap->saip);
            acc_new->sap->saip = NULL;
         }
         if (allsap != NULL)
         {
            allsap_prev->next = acc_new->sap;
            allsap_prev = allsap_prev->next;;
         } else
            allsap_prev = allsap = acc_new->sap;
         acc_new = acc_new->next;
      }
      if (allsap_prev != NULL)
      {
         while (allsap_prev->next != NULL)
         {
            allsap_prev = allsap_prev->next;
         }
      }  
      if (need_to_unlock) 
      {
         BioseqUnlock(bsp_tmp);
         need_to_unlock = FALSE;
      }
      acc = acc->next;
      if (!acc_new_head_head) {
        acc_new_head_head = acc_new_head;
      }
   }

   BLAST_SummaryOptionsFree(options);
   if (need_to_unlock)
   {
      BioseqUnlock (bsp_tmp);
      need_to_unlock = FALSE;
   }

   acc = acc_head;
   while (acc->next != NULL)
   {
      acc = acc->next;
   }
   acc->next = acc_new_head_head;
   if (allsap == NULL)
      return NULL;
   sap_new = SeqAlignNew();
   sap_new->segtype = SAS_DISC;
   sap_new->segs = (Pointer)(allsap);
   allsap = sap_new;
   AlnMgr2IndexLite(allsap);
   AlnMgr2SortAlnSetByNthRowPos(allsap, 1);
   amaip = (AMAlignIndex2Ptr)(allsap->saip);
   for (i=0; i<amaip->numsaps-1; i++)
   {
      amaip->saps[i]->next = amaip->saps[i+1];
   }
   amaip->saps[amaip->numsaps-1]->next = NULL;
   allsap->segs = (Pointer)(amaip->saps[0]);
   return allsap;
}

static void DoCreateSeqHistTPA (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  CKA_RunChecker(sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern void CreateSeqHistTPA (IteM i);
extern void CreateSeqHistTPA (IteM i)
{
  CKA_blast_wordsize = 28;
  CKA_blast_expect_value = 0.000001;
  CKA_blast_allow_repeats = FALSE;
  DoCreateSeqHistTPA(i);
}

static SeqAlignPtr DeltaSeq2SeqAlign(BioseqPtr bsp)
{
  DeltaSeqPtr  deltasp;
  SeqLocPtr    slp;
  SeqAlignPtr  sap, sap_head = NULL;
  DenseSegPtr  dsp;
  SeqIntPtr    intp;
  SeqLitPtr    litp;
  int          curr_start = 0;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) {
    return NULL;
  }
  if (!(deltasp = (DeltaSeqPtr) bsp->seq_ext)) {
    return NULL;
  }

  while (deltasp) {
    if (deltasp->choice == 1) { 
      slp = (SeqLocPtr) deltasp->data.ptrvalue;
      if (sap_head) {
        sap = sap->next = SeqAlignNew();
      } else {
        sap_head = sap = SeqAlignNew();
      }
      dsp = DenseSegNew();
      
      sap->type = SAT_PARTIAL;
      sap->segtype = SAS_DENSEG;
      sap->dim = 2;
      sap->segs = (Pointer)(dsp);
      dsp->dim = 2;
      dsp->numseg = 1;
      dsp->lens = (Int4Ptr)MemNew((dsp->numseg)*sizeof(Int4));
      dsp->starts = (Int4Ptr)MemNew((dsp->numseg)*(dsp->dim)*sizeof(Int4));
      dsp->strands = (Uint1Ptr)MemNew((dsp->numseg)*(dsp->dim)*sizeof(Int4));
      
      dsp->ids = SeqIdDup(bsp->id);
      if (dsp->ids->next) {
        /* Dense-seg ids do not support lists, only 1 id per sequence */
        SeqIdFree(dsp->ids->next);
        dsp->ids->next = NULL;
      }
      switch (slp->choice) {
      case SEQLOC_INT:
        intp = (SeqIntPtr) slp->data.ptrvalue;
        dsp->starts[0] = curr_start;
        dsp->starts[1] = intp->from;
        curr_start += dsp->lens[0] = intp->to - intp->from + 1;
        dsp->strands[0] = Seq_strand_plus;
        dsp->strands[1] = intp->strand;
        dsp->ids->next = SeqIdDup(intp->id);
        break;
      default:
        /* exception */
        break;
      }
    } else if (deltasp->choice == 2) { 
      litp = (SeqLitPtr) deltasp->data.ptrvalue;
      if (litp != NULL) {
        curr_start += litp->length;
      }
    }
    deltasp = deltasp->next;
  }
  return sap_head;
}

static void DoDeltaHist (BioseqPtr bsp, Pointer userdata)

{
  SeqAlignPtr  salp;
  SeqHistPtr   shp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;
  shp = bsp->hist;
  if (shp != NULL && shp->assembly != NULL) return;
  salp = DeltaSeq2SeqAlign (bsp);
  if (salp == NULL) return;
  if (shp == NULL) {
    shp = SeqHistNew ();
    bsp->hist = shp;
  }
  if (shp == NULL) return;
  shp->assembly = salp;
}

extern void CreateSeqHistDelta (IteM i);
extern void CreateSeqHistDelta (IteM i)
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

  VisitBioseqsInSep (sep, NULL, DoDeltaHist);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static TexT blast_wordsize_text = NULL;
static TexT blast_expect_value_text = NULL;
static ButtoN blast_allow_repeats_button = NULL;
static IteM blast_i;

static void DoAcceptBlastOptions (ButtoN b)

{
   Char    buf [64];
   long    val1;
   FloatHi val2;

   GetTitle (blast_wordsize_text, buf, sizeof (buf));
   if (sscanf (buf, "%ld", &val1) == 1) {
     CKA_blast_wordsize = CKA_blast_detailed_wordsize = (Int4) val1;
   }
   GetTitle (blast_expect_value_text, buf, sizeof (buf));
   if (sscanf (buf, "%lf", &val2) == 1) {
     CKA_blast_expect_value = CKA_blast_detailed_expect_value = (FloatHi) val2;
   }
   CKA_blast_allow_repeats = CKA_blast_detailed_allow_repeats =
     (Boolean) GetStatus (blast_allow_repeats_button);
   Remove (ParentWindow (b));
   DoCreateSeqHistTPA(blast_i);
}

extern void CreateSeqHistTPADetailed (IteM i);
extern void CreateSeqHistTPADetailed (IteM i)
{
   GrouP   c;
   GrouP   g;
   GrouP   h;
   WindoW  w;
   Char    buf[64];

   blast_i = i;

   w = FixedWindow (-50, -33, -10, -10, "Blast Options", NULL);
   h = HiddenGroup (w, -1, 0, NULL);
   SetGroupSpacing (h, 3, 2);
   g = HiddenGroup (h, 2, 0, NULL);

   StaticPrompt (g, "Word Size", 0, dialogTextHeight, programFont, 'l');
   if (CKA_blast_detailed_wordsize <= 0) {
     CKA_blast_detailed_wordsize = 14;
   }
   CKA_blast_wordsize = CKA_blast_detailed_wordsize;
   sprintf(buf, "%d", CKA_blast_wordsize);
   blast_wordsize_text = DialogText (g, buf, 10, NULL);

   StaticPrompt (g, "Expect Value", 0, dialogTextHeight, programFont, 'l');
   if (CKA_blast_detailed_expect_value <= 0.0) {
     CKA_blast_detailed_expect_value = 0.001;
   }
   CKA_blast_expect_value = CKA_blast_detailed_expect_value;
   sprintf(buf, "%f", CKA_blast_expect_value);
   blast_expect_value_text = DialogText (g, buf, 10, NULL);

   blast_allow_repeats_button = CheckBox (g, "Allow Repeats", NULL);
   if (CKA_blast_detailed_allow_repeats) {
     SetStatus(blast_allow_repeats_button, TRUE);
   }

   c = HiddenGroup (w, 2, 0, NULL);
   SetGroupSpacing (c, 5, 5);
   DefaultButton (c, "Accept", DoAcceptBlastOptions);
   PushButton (c, "Cancel", StdCancelButtonProc);
   AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);
   RealizeWindow (w);
   Show (w);
   Select (w);
   Select (blast_wordsize_text);
}

extern Int2 LIBCALLBACK AssemblyUserGenFunc (Pointer data);

typedef struct assemblyuserdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        accns;
} AssemblyUserDialog, PNTR AssemblyUserDialogPtr;

typedef struct assemblyuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  SeqDescrPtr   orig_sdp;
} AssemblyUserForm, PNTR AssemblyUserFormPtr;

static void UserObjectPtrToAssemblyDialog (DialoG d, Pointer data)

{
  AssemblyUserDialogPtr  adp;
  Char                   buf [16];
  UserFieldPtr           curr;
  Int4                   from;
  ValNodePtr             head = NULL;
  ObjectIdPtr            oip;
  CharPtr                str;
  CharPtr                tmp;
  Int4                   to;
  UserFieldPtr           ufp;
  UserObjectPtr          uop;

  adp = (AssemblyUserDialogPtr) GetObjectExtra (d);
  if (adp == NULL) return;

  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "TpaAssembly") != 0) {
    PointerToDialog (adp->accns, NULL);
    return;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    str = NULL;
    from = 0;
    to = 0;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
      } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
        from = (Int4) ufp->data.intvalue;
      } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
        to = (Int4) ufp->data.intvalue;
      }
    }
    if (StringHasNoText (str)) continue;
    tmp = MemNew (StringLen (str) + 32);
    StringCpy (tmp, str);
    StringCat (tmp, "\t");
    if (from > 0 || to > 0) {
      sprintf (buf, "%ld", (long) (from + 1));
      StringCat (tmp, buf);
      StringCat (tmp, "\t");
      sprintf (buf, "%ld", (long) (to + 1));
      StringCat (tmp, buf);
      StringCat (tmp, "\t");
    } else {
      StringCat (tmp, "\t\t");
    }
    ValNodeAddStr (&head, 0, (Pointer) tmp);
  }

  PointerToDialog (adp->accns, (Pointer) head);
  ValNodeFreeData (head);
}

static void DoAddAccessionToTpa (UserObjectPtr uop, CharPtr last)

{
  Int4      from = 0;
  CharPtr   ptr1, ptr2, ptr3;
  Int4      to = 0;
  long int  val;

  ptr1 = StringChr (last, '\t');
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = StringChr (ptr1, '\t');
    if (ptr2 != NULL) {
      *ptr2 = '\0';
      ptr2++;
      ptr3 = StringChr (ptr2, '\t');
      if (ptr3 != NULL) {
        *ptr3 = '\0';
      }
    }
    if (sscanf (ptr1, "%ld", &val) == 1 && val > 0) {
      from = val - 1;
      if (sscanf (ptr2, "%ld", &val) == 1 && val > 0) {
        to = val - 1;
      } else {
        from = 0;
        to = 0;
      }
    }
  }
  AddAccessionToTpaAssemblyUserObject (uop, last, from, to);
}

static Pointer AssemblyDialogToUserObjectPtr (DialoG d)

{
  AssemblyUserDialogPtr  adp;
  Char                   ch;
  ValNodePtr             head;
  CharPtr                last;
  UserObjectPtr          uop;
  CharPtr                ptr;
  CharPtr                str;
  CharPtr                tmp;
  ValNodePtr             vnp;

  adp = (AssemblyUserDialogPtr) GetObjectExtra (d);
  if (adp == NULL) return NULL;

  uop = CreateTpaAssemblyUserObject ();
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (adp->accns);
  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',' || ch == ' ') {
        *ptr = '\0';
        TrimSpacesAroundString (last);
        DoAddAccessionToTpa (uop, last);
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      DoAddAccessionToTpa (uop, last);
    }
    MemFree (tmp);
  }

  ValNodeFreeData (head);

  return uop;
}

static void ValNodePtrToAssemblyDialog (DialoG d, Pointer data)

{
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   list;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (ValNodePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        str = MemNew (StringLen ((CharPtr) list->data.ptrvalue) + 3);
        if (str != NULL) {
          StringCpy (str, (CharPtr) list->data.ptrvalue);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer AssemblyDialogToValNodePtr (DialoG d)

{
  Char         ch;
  ValNodePtr   head;
  Int2         j;
  Int2         len;
  ValNodePtr   list;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    list = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        list = ValNodeNew (list);
        if (head == NULL) {
          head = list;
        }
        if (list != NULL) {
          list->choice = 0;
          list->data.ptrvalue = StringSave ((CharPtr) vnp->data.ptrvalue);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 assmbly_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 assmbly_widths [] = {
  16, 8, 8, 0
};

static DialoG CreateAssemblyDialog (GrouP g)

{
  AssemblyUserDialogPtr  adp;
  GrouP                  p;
  GrouP                  x;
  GrouP                  y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  adp = (AssemblyUserDialogPtr) MemNew (sizeof (AssemblyUserDialog));
  if (adp == NULL) return NULL;

  SetObjectExtra (p, adp, NULL);
  adp->dialog = (DialoG) p;
  adp->todialog = UserObjectPtrToAssemblyDialog;
  adp->fromdialog = AssemblyDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  StaticPrompt (y, "Accessions", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "From", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "To", 8 * stdCharWidth, 0, programFont, 'c');
  adp->accns = CreateTagListDialog (x, 3, 3, -1,
                                    assmbly_types, assmbly_widths, NULL,
                                    ValNodePtrToAssemblyDialog,
                                    AssemblyDialogToValNodePtr);

  return (DialoG) p;
}

static void AssemblyUserFormMessage (ForM f, Int2 mssg)

{
  AssemblyUserFormPtr  afp;

  afp = (AssemblyUserFormPtr) GetObjectExtra (f);
  if (afp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (afp->appmessage != NULL) {
          afp->appmessage (f, mssg);
        }
        break;
    }
  }
}


static Int4 CountNumberSequence (CharPtr str)
{
  if (StringHasNoText (str)) return 0;
  else return StringSpn (str, "0123456789");
}

static Int4 CountCapitalLetterSequence (CharPtr str)
{
  if (StringHasNoText (str)) return 0;
  else return StringSpn (str, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
}

static Boolean LooksLikeTPAAccession (CharPtr str)
{
  Int4 len, num_len, letter_len;

  if (StringHasNoText (str)) return TRUE;

  len = StringLen (str);

  if (StringNICmp (str, "ti", 2) == 0 && CountNumberSequence (str + 2) == len - 2)
  {
    return TRUE;
  }
  else if (StringChr (str, '_') != NULL)
  {
    return FALSE;
  }
  letter_len = CountCapitalLetterSequence (str);
  num_len = CountNumberSequence(str + letter_len);

  if ((letter_len == 1 && num_len == 5)
      || (letter_len == 2 && num_len == 6)
      || (letter_len == 4 && num_len == 8))
  {
    if (str[letter_len + num_len] == 0
        || (str[letter_len + num_len] == '.' && letter_len + num_len + 1 + CountNumberSequence (str + letter_len + num_len + 1) == len))
    {
      return TRUE;
    }
  }
  return FALSE;
}


static void TPAAssemblyFormAccept (ButtoN b)

{
  AssemblyUserFormPtr afp;
  UserObjectPtr     uop;
  UserFieldPtr      curr, ufp;
  LogInfoPtr        lip;
  CharPtr           str;
  Boolean           some_bad = FALSE;
  ObjectIdPtr       oip;

  afp = (AssemblyUserFormPtr) GetObjectExtra (b);

  uop = (UserObjectPtr) DialogToPointer (afp->data);
  if (uop == NULL) 
  {
    StdAcceptFormButtonProc (b);
    return;
  }

  lip = OpenLog ("Possible Problem Accessions");

  /* put explanatory text here, but do not set data_in_log flag until
   * suspect accession is actually found
   */
  fprintf (lip->fp, "The information you have input does not appear to be a GenBank accession number.\n");
  fprintf (lip->fp, "Please confirm that you have identified the primary sequences used to assemble \n"
                    "or derive your TPA sequences with GenBank accession numbers in the correct format.\n");

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    str = NULL;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (!LooksLikeTPAAccession (str))
        {
          fprintf (lip->fp, "%s\n", str);
          lip->data_in_log = TRUE;
          some_bad = TRUE;
        }
      }
    }
  }

  CloseLog (lip);
  lip = FreeLog (lip);
  uop = UserObjectFree (uop);

  if (some_bad && ANS_NO == Message (MSG_YN, "Some of your accession numbers may not be in the correct format - continue anyway?"))
  {
    Select (afp->form);    
    return;
  }

  StdAcceptFormButtonProc (b);
}


static void PopulateAssemblyIntervals (ButtoN b)
{
  AssemblyUserFormPtr  afp;
  BioseqPtr            bsp;
  UserObjectPtr        uop;
  SeqAlignPtr          salp;
  Char                 id_txt[355];
  SeqIdPtr             sip;
  Int4                 primary_start, primary_stop;

  afp = (AssemblyUserFormPtr) GetObjectExtra (b);
  if (afp == NULL) {
    return;
  }

  bsp = GetSequenceForObject (OBJ_SEQDESC, afp->orig_sdp);
  if (bsp != NULL && bsp->hist != NULL && bsp->hist->assembly != NULL) {
    uop = CreateTpaAssemblyUserObject  ();
    /* populate user object with intervals */
    for (salp = bsp->hist->assembly; salp != NULL; salp = salp->next) {
      AlnMgr2IndexSingleChildSeqAlign (salp); 
      sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
      SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
      AlnMgr2GetNthSeqRangeInSA (salp, 2, &primary_start, &primary_stop);
      AddAccessionToTpaAssemblyUserObject (uop, id_txt, primary_start, primary_stop);
      sip = SeqIdFree (sip);
    }

    PointerToDialog (afp->data, uop);
    uop = UserObjectFree (uop);
  }
}


static ForM CreateAssemblyDescForm (Int2 left, Int2 top, Int2 width,
                                   Int2 height, CharPtr title, ValNodePtr sdp,
                                   SeqEntryPtr sep, FormActnFunc actproc)

{
  AssemblyUserFormPtr  afp;
  ButtoN               b, pop_btn = NULL;
  GrouP                c;
  GrouP                g;
  StdEditorProcsPtr    sepp;
  WindoW               w;
  BioseqPtr            bsp;

  w = NULL;
  afp = (AssemblyUserFormPtr) MemNew (sizeof (AssemblyUserForm));
  if (afp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, afp, StdDescFormCleanupProc);
    afp->form = (ForM) w;
    afp->actproc = actproc;
    afp->formmessage = AssemblyUserFormMessage;

    afp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      afp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    afp->data = CreateAssemblyDialog (g);

    if (sdp != NULL) {
      bsp = GetSequenceForObject (OBJ_SEQDESC, sdp);
      if (bsp != NULL && bsp->hist != NULL && bsp->hist->assembly != NULL) {
        pop_btn = PushButton (g, "Populate Intervals from Assembly Alignment", PopulateAssemblyIntervals);
        SetObjectExtra (pop_btn, afp, NULL);
      }
    }

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", TPAAssemblyFormAccept);
    SetObjectExtra (b, afp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) pop_btn, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}


extern Int2 LIBCALLBACK AssemblyUserGenFunc (Pointer data)

{
  AssemblyUserFormPtr  afp;
  ObjectIdPtr          oip;
  OMProcControlPtr     ompcp;
  OMUserDataPtr        omudp;
  ObjMgrProcPtr        proc;
  ValNodePtr           sdp;
  SeqEntryPtr          sep;
  UserObjectPtr        uop;
  WindoW               w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit Assembly User Desc") == 0) {
      afp = (AssemblyUserFormPtr) omudp->userdata.ptrvalue;
      if (afp != NULL) {
        Select (afp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "TpaAssembly") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateAssemblyDescForm (-50, -33, -10, -10,
                                       "Assembly Tracking", sdp, sep,
                                       StdDescFormActnProc);
  afp = (AssemblyUserFormPtr) GetObjectExtra (w);
  if (afp != NULL) {
    afp->input_entityID = ompcp->input_entityID;
    afp->input_itemID = ompcp->input_itemID;
    afp->input_itemtype = ompcp->input_itemtype;
    afp->this_itemtype = OBJ_SEQDESC;
    afp->this_subtype = Seq_descr_user;
    afp->procid = ompcp->proc->procid;
    afp->proctype = ompcp->proc->proctype;
    afp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, afp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) afp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (afp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (afp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) afp);
      afp->orig_sdp = sdp;
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


/* advanced editor for Seq-hist assembly alignment */
typedef struct assemblyalignmentdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG intervals_dialog;

  
} AssemblyAlignmentDlgData, PNTR AssemblyAlignmentDlgPtr;

Uint2 assmbly_aln_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_PROMPT, TAGLIST_POPUP
};

Uint2 assmbly_aln_widths [] = {
  16, 8, 8, 8, 8, 8, 0
};

ENUM_ALIST(assmbly_aln_strand_alist)
  {"Plus",  0},
  {"Minus",    1},
END_ENUM_ALIST



static EnumFieldAssocPtr assmbly_aln_alists[] = {
  NULL, NULL, NULL, NULL, NULL, assmbly_aln_strand_alist
};


typedef struct assemblyalignmentinterval {
  CharPtr prim_accession;
  Int4    tpa_from;
  Int4    tpa_to;
  Int4    prim_from;
  Uint1   prim_strand;
  Uint2   percent_identity;
} AssemblyAlignmentIntervalData, PNTR AssemblyAlignmentIntervalPtr;


static AssemblyAlignmentIntervalPtr AssemblyAlignmentIntervalFree (AssemblyAlignmentIntervalPtr interval)
{
  if (interval != NULL) {
    interval->prim_accession = MemFree (interval->prim_accession);
    interval = MemFree (interval);
  }
  return interval;
}


static AssemblyAlignmentIntervalPtr AssemblyAlignmentIntervalFromTagListString (CharPtr str)
{
  AssemblyAlignmentIntervalPtr interval;
  CharPtr cp;
  Int4    len, val;

  if (StringHasNoText (str)) {
    return NULL;
  }

  interval = (AssemblyAlignmentIntervalPtr) MemNew (sizeof (AssemblyAlignmentIntervalData));
  MemSet (interval, 0, sizeof (AssemblyAlignmentIntervalData));

  cp = StringChr (str, '\t');
  if (cp == NULL) {
    interval->prim_accession = StringSave (str);
  } else {
    len = cp - str + 1;
    interval->prim_accession = (CharPtr) MemNew (sizeof (Char) * len);
    StringNCpy (interval->prim_accession, str, len - 1);
    interval->prim_accession[len - 1] = 0;
    str = cp + 1;
    cp = StringChr (str, '\t');
    interval->tpa_from = atoi (str);
    if (cp != NULL) {
      str = cp + 1;
      cp = StringChr (str, '\t');
      interval->tpa_to = atoi (str);
      if (cp != NULL) {
        str = cp + 1;
        cp = StringChr (str, '\t');
        interval->prim_from = atoi (str);
        if (cp != NULL) {
          cp = StringChr (cp + 1, '\t');
          if (cp != NULL) {
            val = atoi (cp + 1);
            if (val == 1) {
              interval->prim_strand = Seq_strand_minus;
            } else {
              interval->prim_strand = Seq_strand_plus;
            }
          }
        }
      }
    }
  }
  return interval;
}


static CharPtr TagListStringFromAssemblyAlignmentInterval (AssemblyAlignmentIntervalPtr interval)
{
  CharPtr str, str_fmt = "%s\t%d\t%d\t%d\t%d (%d)\t%d\n";

  if (interval == NULL) {
    return NULL;
  }

  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (str_fmt) + StringLen (interval->prim_accession) + 61));
  sprintf (str, str_fmt, interval->prim_accession == NULL ? "" : interval->prim_accession,
                         interval->tpa_from,
                         interval->tpa_to,
                         interval->prim_from,
                         interval->prim_from + interval->tpa_to - interval->tpa_from,
                         interval->percent_identity,
                         interval->prim_strand == Seq_strand_minus ? 1 : 0);
  return str;
}


static AssemblyAlignmentIntervalPtr AssemblyAlignmentIntervalFromSeqAlign (SeqAlignPtr salp)
{
  AssemblyAlignmentIntervalPtr interval;
  DenseSegPtr dsp;
  Char     id_txt[200];

  if (salp == NULL || salp->dim != 2 || salp->segtype != SAS_DENSEG) {
    return NULL;
  }

  dsp = (DenseSegPtr) salp->segs;
  if (dsp == NULL || dsp->numseg != 1) {
    return NULL;
  }

  interval = (AssemblyAlignmentIntervalPtr) MemNew (sizeof (AssemblyAlignmentIntervalData));
  MemSet (interval, 0, sizeof (AssemblyAlignmentIntervalData));

  /* first row is TPA, second row is primary */  
  SeqIdWrite (dsp->ids->next, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  interval->prim_accession = StringSave (id_txt);

  interval->tpa_from = dsp->starts[0] + 1;
  interval->tpa_to = dsp->starts[0] + dsp->lens[0];
  interval->prim_from = dsp->starts[1] + 1;
  if (dsp->strands == NULL) {
    interval->prim_strand = Seq_strand_plus;
  } else {
    interval->prim_strand = dsp->strands[1];
  }

  interval->percent_identity = AlignmentPercentIdentity (salp, FALSE);

  return interval;
}


static SeqAlignPtr SeqAlignFromAssemblyAlignmentInterval (AssemblyAlignmentIntervalPtr interval)
{
  SeqAlignPtr salp;
  DenseSegPtr dsp;

  if (interval == NULL) {
    return NULL;
  }
  dsp = DenseSegNew ();
  dsp->dim = 2;
  dsp->numseg = 1;
  dsp->starts = (Int4Ptr) MemNew (sizeof (Int4) * dsp->dim * dsp->numseg);
  dsp->lens = (Int4Ptr) MemNew (sizeof (Int4) * dsp->numseg);
  dsp->strands = (Uint1Ptr) MemNew (sizeof (Uint1) * dsp->dim * dsp->numseg);

  dsp->ids = ValNodeNew (NULL);
  dsp->ids->next = SeqIdFromAccessionDotVersion(interval->prim_accession);

  dsp->starts[0] = interval->tpa_from - 1;
  dsp->starts[1] = interval->prim_from - 1;
  dsp->lens[0] = interval->tpa_to - interval->tpa_from + 1;
  dsp->strands[0] = Seq_strand_plus;
  dsp->strands[1] = interval->prim_strand;

  salp = SeqAlignNew ();
  salp->dim = 2;
  salp->segtype = SAS_DENSEG;
  salp->segs = dsp;
  salp->type = SAT_PARTIAL;
  return salp;
}


static void SeqAlignToAssemblyAlignmentDialog (DialoG d, Pointer data) 
{
  AssemblyAlignmentDlgPtr dlg;

  dlg = (AssemblyAlignmentDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  PointerToDialog (dlg->intervals_dialog, data);
}


static Pointer AssemblyAlignmentDialogToSeqAlign (DialoG d)
{
  AssemblyAlignmentDlgPtr dlg;
  SeqAlignPtr salp = NULL;

  dlg = (AssemblyAlignmentDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  salp = DialogToPointer (dlg->intervals_dialog);
  
  return salp;
}


static void SeqAlignToTagDlg (DialoG d, Pointer data)
{
  TagListPtr tlp;
  SeqAlignPtr salp;
  AssemblyAlignmentIntervalPtr interval;

  tlp =(TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) {
    return;
  }

  tlp->vnp = ValNodeFreeData (tlp->vnp);
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  salp = (SeqAlignPtr) data;

  while (salp != NULL) {
    interval = AssemblyAlignmentIntervalFromSeqAlign (salp);
    if (interval != NULL) {
      ValNodeAddPointer (&(tlp->vnp), 0, TagListStringFromAssemblyAlignmentInterval (interval));
      interval = AssemblyAlignmentIntervalFree (interval);
    }
    salp = salp->next;
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
}


static Pointer SeqAlignFromTagDlg (DialoG d)
{
  TagListPtr tlp;
  SeqAlignPtr salp = NULL, salp_last = NULL, salp_tmp;
  ValNodePtr vnp;
  AssemblyAlignmentIntervalPtr interval;

  tlp =(TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) {
    return NULL;
  }

  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    interval = AssemblyAlignmentIntervalFromTagListString (vnp->data.ptrvalue);
    if (interval != NULL) {
      salp_tmp = SeqAlignFromAssemblyAlignmentInterval (interval);
      if (salp_tmp != NULL) {
        if (salp_last == NULL) {
          salp = salp_tmp;
        } else {
          salp_last->next = salp_tmp;
        }
        salp_last = salp_tmp;
      }
      interval = AssemblyAlignmentIntervalFree (interval);
    }
  }
  return salp;
}


static DialoG CreateAssemblyAlignmentDialog (GrouP g)

{
  AssemblyAlignmentDlgPtr dlg;
  GrouP                   p;
  GrouP                   x;
  GrouP                   y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (AssemblyAlignmentDlgPtr) MemNew (sizeof (AssemblyAlignmentDlgData));
  if (dlg == NULL) return NULL;

  SetObjectExtra (p, dlg, NULL);
  dlg->dialog = (DialoG) p;
  dlg->todialog = SeqAlignToAssemblyAlignmentDialog;
  dlg->fromdialog = AssemblyAlignmentDialogToSeqAlign;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 6, 0, NULL);
  StaticPrompt (y, "", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Primary To", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Accessions", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "TPA From", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "TPA To", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Primary From", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "(Percent Identity)", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Primary Strand", 8 * stdCharWidth, 0, programFont, 'c');
  dlg->intervals_dialog = CreateTagListDialogExEx (x, 6, 6, -1, assmbly_aln_types, assmbly_aln_widths, assmbly_aln_alists,
                                                   TRUE, FALSE, SeqAlignToTagDlg, SeqAlignFromTagDlg,
                                                   NULL, NULL, FALSE);


  return (DialoG) p;
}


static void AddTPAIdToAlignment (SeqAlignPtr salp, BioseqPtr bsp)
{
  DenseSegPtr dsp;
  SeqIdPtr    sip;

  if (salp == NULL || salp->dim != 2 || salp->segtype != SAS_DENSEG || bsp == NULL) {
    return;
  }

  dsp = salp->segs;
  if (dsp == NULL) {
    return;
  }

  sip = SeqIdFindWorst (bsp->id);
  sip = SeqIdDup (sip);

  sip->next = dsp->ids->next;
  dsp->ids->next = NULL;
  dsp->ids = SeqIdFree (dsp->ids);
  dsp->ids = sip;
}


static void AddTPAIdToAlignmentList (SeqAlignPtr salp, BioseqPtr bsp)
{
  while (salp != NULL) {
    AddTPAIdToAlignment (salp, bsp);
    salp = salp->next;
  }
}

typedef struct assemblyalignmentform {
  FORM_MESSAGE_BLOCK
  DialoG dlg;
  BioseqPtr bsp;
} AssemblyAlignmentFormData, PNTR AssemblyAlignmentFormPtr;


static void CheckCoverageWithAddedIntervals (LogInfoPtr lip, BioseqPtr bsp, SeqAlignPtr salp)
{
  ValNodePtr err_list = NULL;
  SeqAlignPtr salp_orig, salp_prev = NULL;

  if (bsp == NULL) {
    return;
  }

  if (bsp->hist == NULL) {
    bsp->hist = SeqHistNew ();
  }

  salp_orig = bsp->hist->assembly;
  while (salp_orig != NULL) {
    salp_prev = salp_orig;
    salp_orig = salp_orig->next;
  }

  if (salp_prev == NULL) {
    bsp->hist->assembly = salp;
  } else {
    salp_prev->next = salp;
  }

  ValidateTPAHistAlign (bsp, &err_list);
  if (err_list != NULL) {
    fprintf (lip->fp, "Projected Coverage Problems\n");
    PrintTPAHistErrors (lip, err_list);
    lip->data_in_log = TRUE;
    err_list = ValNodeFreeData (err_list);
  }

  if (salp_prev == NULL) {
    bsp->hist->assembly = NULL;
  } else {
    salp_prev->next = NULL;
  }
}


static Boolean ReportAssemblyIntervalProblems (BioseqPtr bsp, SeqAlignPtr salp)
{
  LogInfoPtr lip;
  AssemblyAlignmentIntervalPtr interval;
  Boolean    has_errors = FALSE;
  SeqAlignPtr salp_tmp;

  lip = OpenLog ("Assembly Alignment Interval Problems");
  fprintf (lip->fp, "Primary Accession\tTPA From\tTPA To\tPrimary From\tPrimary To\tStrand\tPercent Identity\n");
  for (salp_tmp = salp; salp_tmp != NULL; salp_tmp = salp_tmp->next) {
    interval = AssemblyAlignmentIntervalFromSeqAlign (salp_tmp);
    fprintf (lip->fp, "%s\t%d\t%d\t%d\t%d\t%s\t%d%s\n",
             interval->prim_accession,
             interval->tpa_from,
             interval->tpa_to,
             interval->prim_from,
             interval->prim_from + interval->tpa_to - interval->tpa_from,
             interval->prim_strand == Seq_strand_minus ? "c" : "",
             interval->percent_identity,
             interval->percent_identity < 75 ? "(Suspiciously low percent identity!)" : "");
    if (interval->percent_identity < 75) {
      lip->data_in_log = TRUE;
    }
    interval = AssemblyAlignmentIntervalFree (interval);
  }
  fprintf (lip->fp, "\n");
  CheckCoverageWithAddedIntervals (lip, bsp, salp);
  CloseLog (lip);
  has_errors = lip->data_in_log;
  return has_errors;
}


static void AcceptAssemblyAlignment (ButtoN b)
{
  AssemblyAlignmentFormPtr frm;
  SeqAlignPtr              salp, salp_next;
  ValNodePtr               list;
  MsgAnswer                ans = ANS_OK;

  frm = (AssemblyAlignmentFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }
  salp = DialogToPointer (frm->dlg);
  if (salp == NULL) {
    Message (MSG_ERROR, "No intervals specified");
    return;
  }
  AddTPAIdToAlignmentList (salp, frm->bsp);
  if (ReportAssemblyIntervalProblems (frm->bsp, salp)) {
    ans = Message (MSG_OKC, "Continue with errors?");
  }

  if (ans == ANS_OK) {
    if (frm->bsp->hist == NULL) {
      frm->bsp->hist = SeqHistNew ();
    }

    /* something is wrong here, alignment is not being sorted */
    salp_next = salp;
    while (salp_next->next != NULL) {
      salp_next = salp_next->next;
    }
    salp_next->next = frm->bsp->hist->assembly;

    list = SeqAlignSortListNew (salp);
    list = ValNodeSort (list, SortVnpBySeqAlignSortRow1);
    salp = SeqAlignFromSeqAlignSortList (list);
    list = SeqAlignSortListFree (list);
    frm->bsp->hist->assembly = salp;
    ObjMgrSetDirtyFlag (frm->bsp->idx.entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, frm->bsp->idx.entityID, 0, 0);
    Remove (frm->form);
  } else {
    while (salp != NULL) {
      salp_next = salp->next;
      salp->next = NULL;
      salp = SeqAlignFree (salp);
      salp = salp_next;
    }
  }
}


static void CheckAssemblyAlignment (ButtoN b)
{
  AssemblyAlignmentFormPtr frm;
  SeqAlignPtr              salp, salp_next;

  frm = (AssemblyAlignmentFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  salp = DialogToPointer (frm->dlg);
  if (salp == NULL) {
    Message (MSG_ERROR, "No intervals specified");
    return;
  }
  AddTPAIdToAlignmentList (salp, frm->bsp);
  PointerToDialog (frm->dlg, salp);
  ReportAssemblyIntervalProblems (frm->bsp, salp);
  
  while (salp != NULL) {
    salp_next = salp->next;
    salp->next = NULL;
    salp = SeqAlignFree (salp);
    salp = salp_next;
  }
}


extern void AdvancedAssemblyAlignmentEditor (IteM i)
{
  BaseFormPtr        bfp;
  BioseqPtr   bsp;
  WindoW      w;
  GrouP       h, c;
  AssemblyAlignmentFormPtr frm;
  ButtoN                   b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) {
    Message (MSG_ERROR, "Must select single Bioseq!");
    return;
  }

  frm = (AssemblyAlignmentFormPtr) MemNew (sizeof (AssemblyAlignmentFormData));
  frm->bsp = bsp;

  w = FixedWindow (-50, -33, -10, -10, "Add Intervals to Assembly Alignment", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  frm->dlg = CreateAssemblyAlignmentDialog(h);

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", AcceptAssemblyAlignment);
  SetObjectExtra (b, frm, NULL);
  b = PushButton (c, "Check", CheckAssemblyAlignment);
  SetObjectExtra (b, frm, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->dlg, (HANDLE) c, NULL);

  Show (w);
  Update ();
}

typedef enum {
  eAssemblyIntervalInfo_NoAction = 0,
  eAssemblyIntervalInfo_Remove,
  eAssemblyIntervalInfo_Truncate_Left,
  eAssemblyIntervalInfo_Truncate_Right,
  eAssemblyIntervalInfo_Truncate_Both
} EAssemblyIntervalInfoAction;


typedef struct assemblyintervalinfo {
  SeqAlignPtr salp;
  Int4        prim_left;
  Int4        prim_right;
  Uint1       prim_strand;
  Int4        tpa_left;
  Int4        tpa_right;
  CharPtr     prim_id;
  ValNodePtr  conflict_list;
  EAssemblyIntervalInfoAction action;
} AssemblyIntervalInfoData, PNTR AssemblyIntervalInfoPtr;

typedef enum {
  eIntervalConflict_none,
  eIntervalConflict_a_contains_b,
  eIntervalConflict_a_contained_in_b,
  eIntervalConflict_a_overlaps_b_on_5,
  eIntervalConflict_a_overlaps_b_on_3
} EIntervalConflict;


static Int4 FindIntervalConflict (Int4 left1, Int4 right1, Int4 left2, Int4 right2, Int4 overlap)
{
  if (right1 < left2 + overlap || left1 > right2 - overlap) {
    return eIntervalConflict_none;
  } else if (left1 <= left2 && right1 < right2 && right1 - left2  + 1> overlap) {
    return eIntervalConflict_a_overlaps_b_on_5;
  } else if (left1 > left2 && right1 >= right2 && right2 - left1  + 1> overlap) {
    return eIntervalConflict_a_overlaps_b_on_3;
  } else if (left1 <= left2 && right1 >= right2) {
    return eIntervalConflict_a_contains_b;
  } else if (left2 <= left1 && right2 >= right2) {
    return eIntervalConflict_a_contained_in_b;
  } else {
    Message (MSG_ERROR, "Conflict calculation failed");
    return eIntervalConflict_none;
  }
}



typedef struct intervalconflictinfo {
  AssemblyIntervalInfoPtr conflict_interval;
  EIntervalConflict prim_conflict;
  EIntervalConflict tpa_conflict;
} IntervalConflictInfoData, PNTR IntervalConflictInfoPtr;


static IntervalConflictInfoPtr IntervalConflictInfoNew (AssemblyIntervalInfoPtr conflict_interval, EIntervalConflict prim_conflict, EIntervalConflict tpa_conflict)
{
  IntervalConflictInfoPtr ip;

  ip = (IntervalConflictInfoPtr) MemNew (sizeof (IntervalConflictInfoData));
  ip->conflict_interval = conflict_interval;
  ip->prim_conflict = prim_conflict;
  ip->tpa_conflict = tpa_conflict;
  return ip;
}


static IntervalConflictInfoPtr IntervalConflictInfoFree (IntervalConflictInfoPtr ip)
{
  if (ip != NULL) {
    ip = MemFree (ip);
  }
  return ip;
}


static ValNodePtr IntervalConflictInfoListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->data.ptrvalue = IntervalConflictInfoFree (vnp->data.ptrvalue);
    vnp->next = NULL;
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static AssemblyIntervalInfoPtr AssemblyIntervalInfoFree (AssemblyIntervalInfoPtr ip)
{
  if (ip != NULL) {
    ip->prim_id = MemFree (ip->prim_id);
    ip->conflict_list = IntervalConflictInfoListFree (ip->conflict_list);
    ip = MemFree (ip);
  }
  return ip;
}


static ValNodePtr AssemblyIntervalInfoListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->data.ptrvalue = AssemblyIntervalInfoFree (vnp->data.ptrvalue);
    vnp->next = NULL;
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void TruncateForConflictsOnLeft (AssemblyIntervalInfoPtr ai)
{
  IntervalConflictInfoPtr ip;
  ValNodePtr vnp;
  Int4       prim_change, tpa_change;

  if (ai == NULL || ai->conflict_list == NULL) {
    return;
  }

  for (vnp = ai->conflict_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      ip = (IntervalConflictInfoPtr) vnp->data.ptrvalue;
      prim_change = 0;
      tpa_change = 0;
      if (ip->tpa_conflict == eIntervalConflict_a_overlaps_b_on_5) {
        tpa_change = ip->conflict_interval->tpa_right - ai->tpa_left;
        if (tpa_change > 0) {
          ai->tpa_left += tpa_change;
          if (ai->prim_strand == Seq_strand_minus) {
            ai->prim_right -= tpa_change;
          } else {
            ai->prim_left += tpa_change;
          }
        }
      }
      if (ip->prim_conflict == eIntervalConflict_a_overlaps_b_on_5) {
        prim_change = ip->conflict_interval->prim_right - ai->prim_left; 
        if (prim_change > 0) {
          ai->prim_left += prim_change;
          if (ai->prim_strand == Seq_strand_minus) {
            ai->tpa_right -= prim_change;
          } else {
            ai->tpa_left += prim_change;
          }
        }
      }
      vnp->choice = 0;
    }
  }
}


static void TruncateForConflictsOnRight (AssemblyIntervalInfoPtr ai)
{
  IntervalConflictInfoPtr ip;
  ValNodePtr vnp;
  Int4       prim_change, tpa_change;

  if (ai == NULL || ai->conflict_list == NULL) {
    return;
  }

  for (vnp = ai->conflict_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      ip = (IntervalConflictInfoPtr) vnp->data.ptrvalue;
      prim_change = 0;
      tpa_change = 0;
      if (ip->tpa_conflict == eIntervalConflict_a_overlaps_b_on_3) {
        tpa_change = ai->tpa_right - ip->conflict_interval->tpa_left;
        if (tpa_change > 0) {
          ai->tpa_right -= tpa_change;
          if (ai->prim_strand == Seq_strand_minus) {
            ai->prim_left += tpa_change;
          } else {
            ai->prim_right -= tpa_change;
          }
        }
      }
      if (ip->prim_conflict == eIntervalConflict_a_overlaps_b_on_3) {
        prim_change = ai->prim_right - ip->conflict_interval->prim_left; 
        if (prim_change > 0) {
          ai->prim_right -= prim_change;
          if (ai->prim_strand == Seq_strand_minus) {
            ai->tpa_left += prim_change;
          } else {
            ai->tpa_right -= prim_change;
          }
        }
      }
      vnp->choice = 0;
    }
  }
}


static void ReevaluateConflicts (AssemblyIntervalInfoPtr ip, Int4 overlap)
{
  ValNodePtr vnp;
  IntervalConflictInfoPtr cp;

  if (ip == NULL) {
    return;
  }

  for (vnp = ip->conflict_list; vnp != NULL; vnp = vnp->next) {
    cp = vnp->data.ptrvalue;
    if (cp->conflict_interval->action == eAssemblyIntervalInfo_Remove) {
      vnp->choice = 0;
    } else if (FindIntervalConflict (ip->tpa_left, ip->tpa_right, 
                                     cp->conflict_interval->tpa_left, cp->conflict_interval->tpa_right, overlap) == eIntervalConflict_none
               && FindIntervalConflict (ip->prim_left, ip->prim_right,
                                        cp->conflict_interval->tpa_left, cp->conflict_interval->tpa_right, overlap) == eIntervalConflict_none) {
      vnp->choice = 0;
    } else {
      vnp->choice = 1;
    }
  }
}


static void RecalculateAssemblyIntervalInfoEndpoints (AssemblyIntervalInfoPtr ip, Int4 overlap)
{
  IntervalConflictInfoPtr cp;
  ValNodePtr vnp;

  if (ip == NULL) {
    return;
  }

  /* calculate endpoints */
  AlnMgr2IndexSingleChildSeqAlign (ip->salp);
  ip->prim_strand = SeqAlignStrand (ip->salp, 2);
  AlnMgr2GetNthSeqRangeInSA (ip->salp, 1, &(ip->tpa_left), &(ip->tpa_right));
  AlnMgr2GetNthSeqRangeInSA (ip->salp, 2, &(ip->prim_left), &(ip->prim_right));
  
  /* apply changes for conflicts and actions */
  if (ip->action == eAssemblyIntervalInfo_Truncate_Left) {
    TruncateForConflictsOnLeft (ip);
  } else if (ip->action == eAssemblyIntervalInfo_Truncate_Right) {
    TruncateForConflictsOnRight (ip);
  } else if (ip->action == eAssemblyIntervalInfo_Truncate_Both) {
    TruncateForConflictsOnLeft (ip);
    TruncateForConflictsOnRight (ip);
  } 

  ReevaluateConflicts (ip, overlap);
  for (vnp = ip->conflict_list; vnp != NULL; vnp = vnp->next) {
    cp = vnp->data.ptrvalue;
    ReevaluateConflicts (cp->conflict_interval, overlap);
  }
}


static AssemblyIntervalInfoPtr AssemblyIntervalInfoNew (SeqAlignPtr salp, Int4 overlap)
{
  AssemblyIntervalInfoPtr ip;
  Char        id_txt[255];
  SeqIdPtr    sip;
  BioseqPtr   bsp;

  if (salp == NULL) {
    return NULL;
  }

  ip = (AssemblyIntervalInfoPtr) MemNew (sizeof (AssemblyIntervalInfoData));
  ip->salp = salp;
  AlnMgr2IndexSingleChildSeqAlign (ip->salp);

  sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sip = SeqIdFree (sip);
    sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
    BioseqUnlock (bsp);
  }

  SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  sip = SeqIdFree (sip);

  ip->prim_id = StringSave (id_txt);
  RecalculateAssemblyIntervalInfoEndpoints (ip, overlap);
  return ip;
}


static CharPtr SummarizeAssemblyInterval (AssemblyIntervalInfoPtr ip)
{
  CharPtr summary = NULL;
  CharPtr fmt = "%d-%d %s %d-%d%s";
  Int4    len;

  if (ip == NULL) {
    return NULL;
  }

  len = StringLen (fmt) + StringLen (ip->prim_id) + 60;
  if (ip->prim_strand == Seq_strand_minus) {
    len += 3;
  }
  summary = (CharPtr) MemNew (sizeof (Char) + len);
  sprintf (summary, fmt, ip->tpa_left + 1, ip->tpa_right + 1,
                         ip->prim_id, 
                         ip->prim_left + 1, ip->prim_right + 1,
                         ip->prim_strand == Seq_strand_minus ? "(c)" : "");
  return summary;
}


static CharPtr SummarizeConflictList (ValNodePtr list)
{
  CharPtr summary = NULL, str, conflict = "Conflict with ";
  ValNodePtr vnp, strings = NULL;
  IntervalConflictInfoPtr ip;
  Int4 len = 0;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1 && vnp->data.ptrvalue != NULL) {
      ip = (IntervalConflictInfoPtr) vnp->data.ptrvalue;
      str = SummarizeAssemblyInterval (ip->conflict_interval);
      ValNodeAddPointer (&strings, 0, str);
      len += StringLen (str) + 3;
    }
  }
  if (strings == NULL) {
    summary = StringSave ("No conflicts");
  } else {
    summary = (CharPtr) MemNew (sizeof (Char) * (StringLen (conflict) + len));
    StringCpy (summary, conflict);
    for (vnp = strings; vnp != NULL; vnp = vnp->next) {
      StringCat (summary, vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        StringCat (summary, ", ");
      }
    }
    strings = ValNodeFreeData (strings);
  }
  return summary;
}


static void RemoveConflictLists (ValNodePtr list)
{
  AssemblyIntervalInfoPtr ip;

  while (list != NULL) {
    ip = (AssemblyIntervalInfoPtr) list->data.ptrvalue;
    ip->conflict_list = IntervalConflictInfoListFree (ip->conflict_list);
    list = list->next;
  }
}


static void BuildConflictLists (ValNodePtr list, Int4 overlap)
{
  ValNodePtr vnp1, vnp2;
  AssemblyIntervalInfoPtr ip1, ip2;
  EIntervalConflict prim_conflict, tpa_conflict;
  IntervalConflictInfoPtr conflict;

  if (list == NULL || list->next == NULL) {
    return;
  }

  RemoveConflictLists (list);

  for (vnp1 = list; vnp1->next != NULL; vnp1 = vnp1->next) {
    ip1 = (AssemblyIntervalInfoPtr) vnp1->data.ptrvalue;
    for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next) {
      ip2 = (AssemblyIntervalInfoPtr) vnp2->data.ptrvalue;
      if (StringCmp (ip1->prim_id, ip2->prim_id) == 0) {
        tpa_conflict = FindIntervalConflict (ip1->tpa_left, ip1->tpa_right, ip2->tpa_left, ip2->tpa_right, overlap);
        prim_conflict = FindIntervalConflict (ip1->prim_left, ip1->prim_right, ip2->prim_left, ip2->prim_right, overlap);
        if (tpa_conflict != eIntervalConflict_none || prim_conflict != prim_conflict) {
          conflict = IntervalConflictInfoNew (ip2, prim_conflict, tpa_conflict);
          ValNodeAddPointer (&(ip1->conflict_list), 1, conflict);
        }
        tpa_conflict = FindIntervalConflict (ip2->tpa_left, ip2->tpa_right, ip1->tpa_left, ip1->tpa_right, overlap);
        prim_conflict = FindIntervalConflict (ip2->prim_left, ip2->prim_right, ip1->prim_left, ip1->prim_right, overlap);
        if (tpa_conflict != eIntervalConflict_none || prim_conflict != prim_conflict) {
          conflict = IntervalConflictInfoNew (ip1, prim_conflict, tpa_conflict);
          ValNodeAddPointer (&(ip2->conflict_list), 1, conflict);
        }
      }
    }
  }
}


static ValNodePtr AssemblyIntervalInfoListFromSeqAlign (SeqAlignPtr salp, Int4 overlap)
{
  ValNodePtr list = NULL;

  while (salp != NULL) {
    ValNodeAddPointer (&list, 0, AssemblyIntervalInfoNew (salp, overlap));
    salp = salp->next;
  }

  BuildConflictLists (list, overlap);
  return list;
}


typedef struct assemblyalignmentintervalresolutiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG intervals_dialog;
  TexT   overlap;
  ValNodePtr list;
} AssemblyAlignmentIntervalResolutionDlgData, PNTR AssemblyAlignmentIntervalResolutionDlgPtr;

CharPtr assmbly_aln_int_res_labels [] = {
  "Action",
  "Alignment",
  "Conflicts With" };

Uint2 assmbly_aln_int_res_types [] = {
  TAGLIST_POPUP, TAGLIST_PROMPT, TAGLIST_PROMPT
};

Uint2 assmbly_aln_int_res_widths [] = {
  10, 15, 30, 0
};

ENUM_ALIST(assmbly_aln_int_res_action_alist)
  {"No change",  eAssemblyIntervalInfo_NoAction},
  {"Remove",     eAssemblyIntervalInfo_Remove},
  {"Truncate Left", eAssemblyIntervalInfo_Truncate_Left},
  {"Truncate Right", eAssemblyIntervalInfo_Truncate_Right},
  {"Truncate Both", eAssemblyIntervalInfo_Truncate_Both},
END_ENUM_ALIST


static EnumFieldAssocPtr assmbly_aln_int_res_alists[] = {
  assmbly_aln_int_res_action_alist, NULL, NULL
};


static CharPtr SummarizeOneIntervalResolutionRow (AssemblyIntervalInfoPtr ip) 
{
  CharPtr str, int_str, conf_str, fmt = "%d\t%s\t%s\n";

  if (ip == NULL) {
    return NULL;
  }

  int_str = SummarizeAssemblyInterval (ip);
  conf_str = SummarizeConflictList (ip->conflict_list);
  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 15 + StringLen (int_str) + StringLen (conf_str)));
  sprintf (str, fmt, ip->action, int_str, conf_str);
  int_str = MemFree (int_str);
  conf_str = MemFree (conf_str);

  return str;
}

static void UpdateConflicts (Pointer userdata)
{
  AssemblyAlignmentIntervalResolutionDlgPtr dlg;
  AssemblyIntervalInfoPtr ip;
  TagListPtr tlp;
  ValNodePtr vnp, vnp_t;
  Int4 action;
  Boolean any_change = FALSE;
  CharPtr str;
  Int4    overlap = 0;
  
  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) userdata;
  
  if (dlg == NULL) {
    return;
  }

  tlp = (TagListPtr) GetObjectExtra (dlg->intervals_dialog);
  if (tlp == NULL) {
    return;
  }

  if (!TextHasNoText (dlg->overlap)) {
    str = SaveStringFromText (dlg->overlap);
    overlap = atoi (str);
    str = MemFree (str);
    if (overlap < 0) {
      overlap = 0;
    }
  }

  /* now update conflict lists */
  for (vnp = dlg->list, vnp_t = tlp->vnp; vnp != NULL && tlp->vnp != NULL; vnp = vnp->next, vnp_t = vnp_t->next) {
    ip = (AssemblyIntervalInfoPtr) vnp->data.ptrvalue;
    str = ExtractTagListColumn ((CharPtr) vnp_t->data.ptrvalue, 0);
    action = atoi (str);
    str = MemFree (str);
    if (action != ip->action) {
      ip->action = action;
      RecalculateAssemblyIntervalInfoEndpoints (ip, overlap);
      any_change = TRUE;
    }
  }

  if (any_change) {
    for (vnp = dlg->list, vnp_t = tlp->vnp; vnp != NULL && tlp->vnp != NULL; vnp = vnp->next, vnp_t = vnp_t->next) {
      ip = (AssemblyIntervalInfoPtr) vnp->data.ptrvalue;
      vnp_t->data.ptrvalue = MemFree (vnp_t->data.ptrvalue);
      vnp_t->data.ptrvalue = SummarizeOneIntervalResolutionRow (ip);
    }
    
    /* update dialog */
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);
  }
}  

static TaglistCallback assmbly_int_res_callback_list[3] = 
 { UpdateConflicts, UpdateConflicts, UpdateConflicts };

static void SeqAlignToAssemblyAlignmentIntervalResolutionDialog (DialoG d, Pointer data)
{
  AssemblyAlignmentIntervalResolutionDlgPtr dlg;
  SeqAlignPtr salp;
  ValNodePtr  vnp;
  CharPtr     str;
  Int4        overlap = 0;
  TagListPtr  tlp;

  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) GetObjectExtra (d);
  salp = (SeqAlignPtr) data;
  if (dlg == NULL) {
    return;
  }
  tlp = (TagListPtr) GetObjectExtra (dlg->intervals_dialog);
  if (tlp == NULL) {
    return;
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = ValNodeFreeData (tlp->vnp);

  if (!TextHasNoText (dlg->overlap)) {
    str = SaveStringFromText (dlg->overlap);
    overlap = atoi (str);
    str = MemFree (str);
    if (overlap < 0) {
      overlap = 0;
    }
  }

  dlg->list = AssemblyIntervalInfoListFree (dlg->list);
  dlg->list = AssemblyIntervalInfoListFromSeqAlign (salp, overlap);

  for (vnp = dlg->list; vnp != NULL; vnp = vnp->next) {
    str = SummarizeOneIntervalResolutionRow (vnp->data.ptrvalue);
    ValNodeAddPointer (&(tlp->vnp), 0, str);
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (ValNodeLen (tlp->vnp) - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
  } else {
    SafeHide (tlp->bar);
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);

}


static Pointer AssemblyAlignmentIntervalResolutionDialogToSeqAlign (DialoG d)
{
  AssemblyAlignmentIntervalResolutionDlgPtr dlg;
  AssemblyIntervalInfoPtr ai;
  SeqAlignPtr salp_list = NULL, salp_prev = NULL, salp;
  ValNodePtr vnp;

  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  for (vnp = dlg->list; vnp != NULL; vnp = vnp->next) {
    ai = vnp->data.ptrvalue;
    if (ai->action != eAssemblyIntervalInfo_Remove) {
      salp = (SeqAlignPtr) AsnIoMemCopy (ai->salp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);
      /* TODO: truncate alignments */

      /* add to list */
      if (salp_prev == NULL) {
        salp_list = salp;
      } else {
        salp_prev->next = salp;
      }
      salp_prev = salp;
    }
  }

  return (Pointer) salp_list;
}


static void ChangeAssemblyAlignmentIntervalResolutionOverlap (TexT t)
{
  AssemblyAlignmentIntervalResolutionDlgPtr  dlg;
  AssemblyIntervalInfoPtr ip;
  TagListPtr tlp;
  ValNodePtr vnp, vnp_t;
  CharPtr    str;
  Int4       overlap = 0;

  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  tlp = (TagListPtr) GetObjectExtra (dlg->intervals_dialog);
  if (tlp == NULL) {
    return;
  }

  if (!TextHasNoText (dlg->overlap)) {
    str = SaveStringFromText (dlg->overlap);
    overlap = atoi (str);
    str = MemFree (str);
    if (overlap < 0) {
      overlap = 0;
    }
  }

  BuildConflictLists (dlg->list, overlap);

  /* change text */
  for (vnp = dlg->list, vnp_t = tlp->vnp; vnp != NULL && tlp->vnp != NULL; vnp = vnp->next, vnp_t = vnp_t->next) {
    ip = (AssemblyIntervalInfoPtr) vnp->data.ptrvalue;
    RecalculateAssemblyIntervalInfoEndpoints (ip, overlap);
    vnp_t->data.ptrvalue = MemFree (vnp_t->data.ptrvalue);
    vnp_t->data.ptrvalue = SummarizeOneIntervalResolutionRow (ip);
  }
  
  /* update dialog */
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);

}


static void CleanupAssemblyAlignmentIntervalResolutionDialog (GraphiC g, VoidPtr data)

{
  AssemblyAlignmentIntervalResolutionDlgPtr  dlg;

  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) data;
  if (dlg != NULL) {
    dlg->list =  AssemblyIntervalInfoListFree (dlg->list);
    dlg = MemFree (dlg);
  }
}


static DialoG CreateAssemblyAlignmentIntervalResolutionDialog (GrouP g)

{
  AssemblyAlignmentIntervalResolutionDlgPtr dlg;
  GrouP                   p;
  GrouP                   x, y, z;
  Int4                    i;
  Int4                    num_columns = sizeof (assmbly_aln_int_res_labels) / sizeof (CharPtr);

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (AssemblyAlignmentIntervalResolutionDlgPtr) MemNew (sizeof (AssemblyAlignmentIntervalResolutionDlgData));
  if (dlg == NULL) return NULL;

  SetObjectExtra (p, dlg, CleanupAssemblyAlignmentIntervalResolutionDialog);
  dlg->dialog = (DialoG) p;
  
  dlg->todialog = SeqAlignToAssemblyAlignmentIntervalResolutionDialog;
  dlg->fromdialog = AssemblyAlignmentIntervalResolutionDialogToSeqAlign;

  z = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (z, "Allowable Overlap", 0, 0, programFont, 'r');
  dlg->overlap = DialogText (z, "5", 5, ChangeAssemblyAlignmentIntervalResolutionOverlap);
  SetObjectExtra (dlg->overlap, dlg, NULL);

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, num_columns, 0, NULL);
  for (i = 0; i < num_columns; i++) {
      StaticPrompt (y, assmbly_aln_int_res_labels[i], assmbly_aln_int_res_widths[i] * stdCharWidth, 0, programFont, 'l');
  }
  dlg->intervals_dialog = CreateTagListDialogExEx (x, 6, num_columns, -1, assmbly_aln_int_res_types, 
                                                   assmbly_aln_int_res_widths, assmbly_aln_int_res_alists,
                                                   TRUE, TRUE, NULL, NULL,
                                                   assmbly_int_res_callback_list, dlg, FALSE);

  AlignObjects (ALIGN_CENTER, (HANDLE) z, (HANDLE) x, NULL);
  return (DialoG) p;
}


typedef struct assemblyalignmentintervalresolutionfrm {
  FORM_MESSAGE_BLOCK
  DialoG dlg;

  BioseqPtr bsp;
} AssemblyAlignmentIntervalResolutionFrmData, PNTR AssemblyAlignmentIntervalResolutionFrmPtr;


static void AcceptAssemblyAlignmentIntervalResolution (ButtoN b)
{
  AssemblyAlignmentIntervalResolutionFrmPtr frm;
  SeqAlignPtr new_assem, salp, salp_next;

  frm = (AssemblyAlignmentIntervalResolutionFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  /* TODO: check for problems before accepting */

  new_assem = DialogToPointer (frm->dlg);
  if (new_assem == NULL) {
    Message (MSG_ERROR, "No intervals left in assembly!");
    return;
  }

  /* remove old assembly */
  salp = frm->bsp->hist->assembly;
  frm->bsp->hist->assembly = NULL;
  while (salp != NULL) {
    salp_next = salp->next;
    salp->next = NULL;
    salp = SeqAlignFree (salp);
    salp = salp_next;
  }

  /* assign new assembly */
  frm->bsp->hist->assembly = new_assem;

  ObjMgrSetDirtyFlag (frm->bsp->idx.entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->bsp->idx.entityID, 0, 0);
  Remove (frm->form);
}


extern void AssemblyAlignmentIntervalResolution (IteM i)
{
  BaseFormPtr        bfp;
  BioseqPtr   bsp;
  WindoW      w;
  GrouP       h, c;
  AssemblyAlignmentIntervalResolutionFrmPtr frm;
  ButtoN                   b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) {
    Message (MSG_ERROR, "Must select single Bioseq!");
    return;
  }
  if (bsp->hist == NULL || bsp->hist->assembly == NULL) {
    Message (MSG_ERROR, "No assembly alignment!");
    return;
  }
  frm = (AssemblyAlignmentIntervalResolutionFrmPtr) MemNew (sizeof (AssemblyAlignmentIntervalResolutionFrmData));
  frm->bsp = bsp;

  w = FixedWindow (-50, -33, -10, -10, "Resolve Assembly Alignment Intervals", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  frm->dlg = CreateAssemblyAlignmentIntervalResolutionDialog(h);
  PointerToDialog (frm->dlg, bsp->hist->assembly);

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", AcceptAssemblyAlignmentIntervalResolution);
  SetObjectExtra (b, frm, NULL);
/*  b = PushButton (c, "Check", CheckAssemblyAlignmentIntervalResolution);
  SetObjectExtra (b, frm, NULL); */
  b = PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->dlg, (HANDLE) c, NULL);
  Show (w);
  Update ();
}



typedef struct historyformdata {
  FEATURE_FORM_BLOCK

  BioseqPtr      bsp;
  DialoG         replace_date;
  DialoG         replace_ids;
  ButtoN         secondary_on_part;
  DialoG         replaced_by_date;
  DialoG         replaced_by_ids;
  ButtoN         deleted;
  DialoG         deleted_date;
} HistoryFormData, PNTR HistoryFormPtr;

static SeqIdPtr VisStrDialogToSeqIdSet (DialoG d)

{
  long        gi;
  SeqIdPtr    head = NULL;
  ValNodePtr  list;
  SeqIdPtr    sip;
  CharPtr     str;
  ValNodePtr  vnp;

  if (d == NULL) return NULL;
  list = DialogToPointer (d);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      if (sscanf (str, "%ld", &gi)) {
        /*
        ValNodeAddInt (&head, SEQID_GI, (Int4) gi);
        */
      } else {
        sip = SeqIdFromAccessionDotVersion (str);
        if (sip != NULL) {
          ValNodeLink (&head, sip);
        }
      }
    }
  }
  ValNodeFreeData (list);
  return head;
}

static int LIBCALLBACK SortByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static ValNodePtr GetStringsForSeqIDs (SeqIdPtr sip)

{
  Char          buf [40];
  ValNodePtr    head = NULL;
  TextSeqIdPtr  tsip;

  if (sip == NULL) return NULL;
  while (sip != NULL) {
    buf [0] = '\0';
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && (! StringHasNoText (tsip->accession))) {
          StringNCpy_0 (buf, tsip->accession, sizeof (buf));
        }
        break;
      case SEQID_GI :
        /*
        gi = sip->data.intvalue;
        if (gi > 0) {
          sprintf (buf, "%ld", (long) gi);
        }
        */
        break;
      default :
        break;
    }
    if (! StringHasNoText (buf)) {
      ValNodeCopyStr (&head, 0, buf);
    }
    sip = sip->next;
  }
  return head;
}

static void AddGenBankBlockToBioseq (BioseqPtr bsp, ValNodePtr head1, ValNodePtr head2)

{
  GBBlockPtr       gbp = NULL;
  CharPtr          last = NULL;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       sdp;
  SeqEntryPtr      sep;
  CharPtr          str1, str2;
  ValNodePtr       vnp, vnp1, vnp2;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
        str1 = (CharPtr) vnp1->data.ptrvalue;
        if (str1 != NULL) {
          for (vnp2 = gbp->extra_accessions; vnp2 != NULL; vnp2 = vnp2->next) {
            str2 = (CharPtr) vnp2->data.ptrvalue;
            if (StringICmp (str1, str2) == 0) {
              vnp2->data.ptrvalue = MemFree (vnp2->data.ptrvalue);
            }
          }
        }
      }
    }
  }
  if (sdp == NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    if (sdp != NULL) {
      sdp->data.ptrvalue = GBBlockNew ();
    }
  }
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      while (head2 != NULL) {
        ValNodeCopyStr (&(gbp->extra_accessions), 0, (CharPtr) head2->data.ptrvalue);
        head2 = head2->next;
      }
      /*
      ValNodeLink (&(gbp->extra_accessions), head2);
      head2 = NULL;
      */
      gbp->extra_accessions = SortValNode (gbp->extra_accessions, SortByName);
      prev = &(gbp->extra_accessions);
      vnp = gbp->extra_accessions;
      last = NULL;
      while (vnp != NULL) {
        next = vnp->next;
        str2 = (CharPtr) vnp->data.ptrvalue;
        if (str2 == NULL || StringHasNoText (str2) || StringICmp (last, str2) == 0) {
          *prev = next;
          vnp->next = NULL;
          MemFree (vnp);
          vnp = next;
        } else {
          last = str2;
          prev = &(vnp->next);
          vnp = next;
        }
      }
    }
  }

}

static void DoChangeHistory (ButtoN b)

{
  MsgAnswer       ans;
  BioseqPtr       bsp;
  ValNodePtr      head1 = NULL, head2 = NULL;
  HistoryFormPtr  hfp;
  SeqHistPtr      hist;
  BioseqPtr       pbsp;
  SeqEntryPtr     sep;
  SeqLocPtr       slp;
  CharPtr         str1, str2;
  ValNodePtr      vnp1, vnp2;

  hfp = (HistoryFormPtr) GetObjectExtra (b);
  if (hfp == NULL) return;
  ans = Message (MSG_OKC, "Are you sure you want to edit the history?");
  if (ans == ANS_CANCEL) {
    return;
  }
  Hide (hfp->form);
  Update ();
  bsp = hfp->bsp;
  hist = bsp->hist;
  if (hist == NULL) {
    hist = SeqHistNew ();
    bsp->hist = hist;
  }
  if (hist != NULL) {

    hist->replace_date = DateFree (hist->replace_date);
    hist->replace_date = DialogToPointer (hfp->replace_date);
    head1 = GetStringsForSeqIDs (hist->replace_ids);
    hist->replace_ids = SeqIdSetFree (hist->replace_ids);
    hist->replace_ids = VisStrDialogToSeqIdSet (hfp->replace_ids);
    head2 = GetStringsForSeqIDs (hist->replace_ids);

    hist->replaced_by_date = DateFree (hist->replaced_by_date);
    hist->replaced_by_date = DialogToPointer (hfp->replaced_by_date);
    hist->replaced_by_ids = SeqIdSetFree (hist->replaced_by_ids);
    hist->replaced_by_ids = VisStrDialogToSeqIdSet (hfp->replaced_by_ids);

    hist->deleted = GetStatus (hfp->deleted);
    hist->deleted_date = DateFree (hist->deleted_date);
    hist->deleted_date = DialogToPointer (hfp->deleted_date);
  }

  if (hist->assembly == NULL &&
      hist->replace_date == NULL && hist->replace_ids == NULL &&
      hist->replaced_by_date == NULL && hist->replaced_by_ids == NULL &&
      (! hist->deleted) && hist->deleted_date == NULL) {
    bsp->hist = SeqHistFree (bsp->hist);
  }

  head1 = SortValNode (head1, SortByName);
  head2 = SortValNode (head2, SortByName);
  for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    for (vnp2 = head2; vnp2 != NULL; vnp2 = vnp2->next) {
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (StringICmp (str1, str2) == 0) {
        vnp1->data.ptrvalue = MemFree (vnp1->data.ptrvalue);
      }
    }
  }

  AddGenBankBlockToBioseq (bsp, head1, head2);

  if (GetStatus (hfp->secondary_on_part)) {
    if (bsp->repr == Seq_repr_seg) {
      for (slp = (SeqLocPtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
        pbsp = BioseqFind (SeqLocId (slp));
        if (pbsp != NULL) {
          AddGenBankBlockToBioseq (pbsp, head1, head2);
        }
      }
    }
  }

  ValNodeFreeData (head1);
  ValNodeFreeData (head2);

  sep = GetTopSeqEntryForEntityID (hfp->input_entityID);
  EntryCheckGBBlock (sep);

  Update ();
  ObjMgrSetDirtyFlag (hfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, hfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (hfp->form);
}

static void SeqIdSetToVisStrDialog (DialoG d, SeqIdPtr sip)

{
  ValNodePtr  head = NULL;

  if (d == NULL || sip == NULL) return;
  head = GetStringsForSeqIDs (sip);
  PointerToDialog (d, head);
  ValNodeFreeData (head);
}

static void HistoryMessageProc (ForM f, Int2 mssg)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) GetObjectExtra (f);
  if (hfp != NULL) {
    if (hfp->appmessage != NULL) {
      hfp->appmessage (f, mssg);
    }
  }
}

static void CleanupHistoryPage (GraphiC g, VoidPtr data)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) data;
  if (hfp != NULL) {
  }
  StdCleanupFormProc (g, data);
}

extern void EditSequenceHistory (IteM i)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  HistoryFormPtr     hfp;
  SeqHistPtr         hist;
  GrouP              j;
  GrouP              k;
  PrompT             ppt1, ppt2, ppt3;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  hist = bsp->hist;

  hfp = (HistoryFormPtr) MemNew (sizeof (HistoryFormData));
  if (hfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Sequence History", StdCloseWindowProc);
  SetObjectExtra (w, hfp, CleanupHistoryPage);
  hfp->form = (ForM) w;
  hfp->formmessage = HistoryMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    hfp->appmessage = sepp->handleMessages;
  }

  hfp->input_entityID = bfp->input_entityID;
  hfp->input_itemID = bfp->input_itemID;
  hfp->input_itemtype = bfp->input_itemtype;

  hfp->bsp = bsp;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt1 = StaticPrompt (h, "Replaces", 0, 0, systemFont, 'c');

  g = HiddenGroup (h, -1, 0, NULL);
  hfp->replace_ids = CreateVisibleStringDialog (g, 4, -1, 10);
  hfp->secondary_on_part = CheckBox (g, "Secondary on Parts", NULL);
  hfp->replace_date = CreateDateDialog (g, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replace_ids, (HANDLE) hfp->secondary_on_part, (HANDLE) hfp->replace_date, NULL);

  ppt2 = StaticPrompt (h, "Replaced By", 0, 0, systemFont, 'c');

  j = HiddenGroup (h, -1, 0, NULL);
  hfp->replaced_by_ids = CreateVisibleStringDialog (j, 4, -1, 10);
  hfp->replaced_by_date = CreateDateDialog (j, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replaced_by_ids, (HANDLE) hfp->replaced_by_date, NULL);

  ppt3 = StaticPrompt (h, "Status", 0, 0, systemFont, 'c');

  k = HiddenGroup (h, -1, 0, NULL);
  hfp->deleted = CheckBox (k, "Deleted", NULL);
  hfp->deleted_date = CreateDateDialog (k, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->deleted, (HANDLE) hfp->deleted_date, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoChangeHistory);
  SetObjectExtra (b, hfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) j, (HANDLE) k, (HANDLE) c,
                (HANDLE) ppt1, (HANDLE) ppt2, (HANDLE) ppt3, NULL);
  RealizeWindow (w);

  if (bsp->repr != Seq_repr_seg) {
    Disable (hfp->secondary_on_part);
  } else {
    SetStatus (hfp->secondary_on_part, TRUE);
  }

  if (hist != NULL) {
    PointerToDialog (hfp->replace_date, hist->replace_date);
    SeqIdSetToVisStrDialog (hfp->replace_ids, hist->replace_ids);
    PointerToDialog (hfp->replaced_by_date, hist->replaced_by_date);
    SeqIdSetToVisStrDialog (hfp->replaced_by_ids, hist->replaced_by_ids);
    PointerToDialog (hfp->deleted_date, hist->deleted_date);
    SetStatus (hfp->deleted, hist->deleted);
  }

  Show (w);
  Update ();
}


typedef struct nw {
  Int4 score;
  Int4 traceback_pos;
} NWData, PNTR NWPtr;


static Int4 SegmentsNeededForAlignment (CharPtr buf1, CharPtr buf2)
{
  CharPtr cp1, cp2;
  Boolean gap1, gap2, change = FALSE;
  Int4    num_segs = 1;

  cp1 = buf1;
  cp2 = buf2;
  if (*cp1 == '-') {
    gap1 = TRUE;
  } else {
    gap1 = FALSE;
  }
  cp1++;
  if (*cp2 == '-') {
    gap2 = TRUE;
  } else {
    gap2 = FALSE;
  }
  cp2++;
  while (*cp1 != 0 && *cp2 != 0) {
    change = FALSE;
    if (*cp1 == '-') {
      if (!gap1) {
        gap1 = TRUE;
        change = TRUE;
      }
    } else if (gap1) {
      gap1 = FALSE;
      change = TRUE;
    }
    if (*cp2 == '-') {
      if (!gap2) {
        gap2 = TRUE;
        change = TRUE;
      }
    } else if (gap2) {
      gap2 = FALSE;
      change = TRUE;
    }
    if (change) {
      num_segs ++;
    }
    cp1++;
    cp2++;
  }
  return num_segs;
}

static void 
FillInStartsAndLensForAlignment 
(DenseSegPtr dsp,
 CharPtr     buf1,
 CharPtr     buf2)
{
  CharPtr cp1, cp2;
  Boolean gap1, gap2, change = FALSE;
  Int4    num_segs = 0;
  Int4    pos1 = 0, pos2 = 0;

  cp1 = buf1;
  cp2 = buf2;
  if (*cp1 == '-') {
    dsp->starts[dsp->dim * num_segs] = -1;
    gap1 = TRUE;
  } else {
    dsp->starts[dsp->dim * num_segs] = pos1;
    gap1 = FALSE;
    pos1++;
  }
  cp1++;
  if (*cp2 == '-') {
    gap2 = TRUE;
    dsp->starts[dsp->dim * num_segs + 1] = -1;
  } else {
    dsp->starts[dsp->dim * num_segs + 1] = pos2;
    gap2 = FALSE;
    pos2++;
  }
  cp2++;
  dsp->lens[num_segs] = 1;

  while (*cp1 != 0 && *cp2 != 0) {
    change = FALSE;
    if (*cp1 == '-') {
      if (!gap1) {
        gap1 = TRUE;
        change = TRUE;
      }
    } else {
      if (gap1) {
        gap1 = FALSE;
        change = TRUE;
      }
    }
    if (*cp2 == '-') {
      if (!gap2) {
        gap2 = TRUE;
        change = TRUE;
      }
    } else {
      if (gap2) {
        gap2 = FALSE;
        change = TRUE;
      }
    }
    if (change) {
      num_segs ++;
      if (gap1) {
        dsp->starts[dsp->dim * num_segs] = -1;
      } else {
        dsp->starts[dsp->dim * num_segs] = pos1;
      }
      if (gap2) {
        dsp->starts[dsp->dim * num_segs + 1] = -1;
      } else {
        dsp->starts[dsp->dim * num_segs + 1] = pos2;
      }
      dsp->lens[num_segs] = 1;
    } else {
      dsp->lens[num_segs]++;
    }
    cp1++;
    cp2++;
    if (!gap1) {
      pos1++;
    }
    if (!gap2) {
      pos2++;
    }
  }
}


static void AdjustStringAlingmentForOffsetAndStrand (DenseSegPtr dsp, Int4 start1, Int4 start2, Int4 stop2, Uint1 strand2)
{
  Int4 i;

  if (dsp == NULL) {
    return;
  }

  for (i = 0; i < dsp->numseg; i++) {
    if (dsp->starts[2 * i] != -1) {
      dsp->starts[2 * i] += start1;
    }
    if (dsp->starts[2 * i + 1] != -1) {
      if (strand2 == Seq_strand_plus) {
        dsp->starts[2 * i + 1] += start2;
      } else {
        dsp->starts[2 * i + 1] = stop2 - dsp->starts[2 * i + 1] - dsp->lens[i];
      }
    }
  }
}


/* assumption - first interval always on plus strand */
static SeqAlignPtr NWAlignmentForInterval (SeqIdPtr sip1, SeqIdPtr sip2, Int4 start1, Int4 stop1, Int4 start2, Int4 stop2)
{
  BioseqPtr bsp1, bsp2;
  Int4      len1, len2, tmp, i, row, col, back_row, back_col;
  Uint1     strand2 = Seq_strand_plus;
  CharPtr   buf1, buf2;
  CharPtr   alnbuf1, alnbuf2, cp1, cp2;
  NWPtr     matrix;
  Int4      gap_penalty = -1;
  Int4      mismatch_penalty = -1;
  Int4      match_score = 1;
  Int4      left, up, diag;
  Int4      num_segs;
  DenseSegPtr dsp;
  SeqAlignPtr salp = NULL;

  if (sip1 == NULL || sip2 == NULL) {
    return NULL;
  }

  bsp1 = BioseqLockById (sip1);
  bsp2 = BioseqLockById (sip2);
  
  if (bsp1 != NULL && bsp2 != NULL) {
    if (stop2 < start2) {
      strand2 = Seq_strand_minus;
      tmp = start2;
      start2 = stop2;
      stop2 = tmp;
    }
    len1 = stop1 - start1 + 1;
    len2 = stop2 - start2 + 1;
    buf1 = MemNew (sizeof (Char) * (len1 + 1));
    buf2 = MemNew (sizeof (Char) * (len2 + 1));
    SeqPortStreamInt (bsp1, start1, stop1, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf1, NULL);
    SeqPortStreamInt (bsp2, start2, stop2, strand2, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);

    matrix = (NWPtr) MemNew (sizeof (NWData) * (len1 + 1) * (len2 + 1));
    /* initalize matrix */
    MemSet (matrix, 0, sizeof (NWData) * (len1 + 1) * (len2 + 1));
    matrix[0].score = 0;
    matrix[0].traceback_pos = 0;
    row = 0;
    for (col = 1; col <= len1; col++) {
      matrix[(row * (len1  + 1)) + col].score = matrix[(row * (len1  + 1)) + col - 1].score + gap_penalty;
      matrix[(row * (len1  + 1)) + col].traceback_pos = (row * (len1  + 1)) + col - 1;
    }
    col = 0;
    for (row = 1; row <= len2; row++) {
      matrix[(row * (len1  + 1)) + col].score = matrix[((row - 1) * (len1  + 1)) + col].score + gap_penalty;
      matrix[(row * (len1  + 1)) + col].traceback_pos = ((row - 1) * (len1  + 1)) + col;
    }

    /* fill in scores */
    for (row = 1; row <= len2; row++) {
      for (col = 1; col <= len1; col++) {
        /* diagonal */
        diag = matrix[((row - 1) * (len1  + 1)) + col - 1].score;
        if (buf1[col - 1] == buf2[row - 1]) {
          diag += match_score;
        } else {
          diag += mismatch_penalty;
        }
        left = matrix[((row) * (len1  + 1)) + col - 1].score + gap_penalty;
        up = matrix[((row - 1) * (len1  + 1)) + col].score + gap_penalty;

        /* choose best */
        if (left > diag && left > up) {
          matrix[((row) * (len1 + 1)) + col].score = left + gap_penalty;
          matrix[((row) * (len1 + 1)) + col].traceback_pos = ((row) * (len1  + 1)) + col - 1;
        } else if (up > diag && up > left) {
          matrix[((row) * (len1 + 1)) + col].score = up + gap_penalty;
          matrix[((row) * (len1 + 1)) + col].traceback_pos = ((row - 1) * (len1  + 1)) + col;
        } else {
          matrix[((row) * (len1 + 1)) + col].score = diag;
          matrix[((row) * (len1 + 1)) + col].traceback_pos = ((row - 1) * (len1  + 1)) + col - 1;
        }
        
      }
    }

    /* trace back, create alignment strings */
    alnbuf1 = (CharPtr) MemNew (sizeof (Char) * (len1 + len2 + 1));
    alnbuf2 = (CharPtr) MemNew (sizeof (Char) * (len1 + len2 + 1));
    cp1 = alnbuf1 + len1 + len2;
    cp2 = alnbuf2 + len1 + len2;
    *cp1 = 0;
    *cp2 = 0;
    cp1--;
    cp2--;
    row = len2;
    col = len1;
    while (row > 0 || col > 0) {
      back_row = matrix[(row * (len1 + 1)) + col].traceback_pos / (len1 + 1);
      back_col = matrix[(row * (len1 + 1)) + col].traceback_pos % (len1 + 1);
      if (row == back_row) {
        *cp1 = buf1[col - 1];
        *cp2 = '-';
      } else if (col == back_col) {
        *cp1 = '-';
        *cp2 = buf2[row - 1];
      } else {
        *cp1 = buf1[col - 1];
        *cp2 = buf2[row - 1];
      }
      cp1--;
      cp2--;
      row = back_row;
      col = back_col;
    }

    /* no longer need matrix or original sequence buffers */
    matrix = MemFree (matrix);
    buf1 = MemFree (buf1);
    buf2 = MemFree (buf2);

    /* count number of segments needed */
    num_segs = SegmentsNeededForAlignment (cp1 + 1, cp2 + 1);
        
    /* create DenseSeg */
    dsp = DenseSegNew ();
    dsp->dim = 2;
    dsp->ids = SeqIdDup (sip1);
    dsp->ids->next = SeqIdDup (sip2);
    dsp->numseg = num_segs;
    dsp->lens = (Int4Ptr)MemNew (sizeof (Int4) * dsp->numseg);
    dsp->starts = (Int4Ptr)MemNew (sizeof (Int4) * dsp->dim * dsp->numseg);
    dsp->strands = (Uint1Ptr)MemNew (sizeof (Uint1) * dsp->dim * dsp->numseg);

    /* fill in strands */
    for (i = 0; i < dsp->numseg; i++) {
      dsp->strands[2 * i] = Seq_strand_plus;
      dsp->strands[2 * i + 1] = strand2;
    }
    /* fill in starts and lens */
    FillInStartsAndLensForAlignment (dsp, cp1 + 1, cp2 + 1);

    /* no longer need FASTA+GAP alignment strings */
    alnbuf1 = MemFree (alnbuf1);
    alnbuf2 = MemFree (alnbuf2);

    /* adjust for real sequence position and strand */
    AdjustStringAlingmentForOffsetAndStrand (dsp, start1, start2, stop2, strand2);

    salp = SeqAlignNew ();
    salp->segs = dsp;
    salp->segtype = SAS_DENSEG;
    salp->dim = 2;
  }
  BioseqUnlock (bsp1);
  BioseqUnlock (bsp2);
  AlnMgr2IndexSingleChildSeqAlign (salp);  
  return salp;
}


/* Assumptions:
 * first interval always on plus strand
 */
static SeqAlignPtr AlignmentForInterval (SeqIdPtr sip1, SeqIdPtr sip2, Int4 start1, Int4 stop1, Int4 start2, Int4 stop2)
{
  DenseSegPtr dsp;
  SeqAlignPtr salp = NULL;
  Int4        len1, len2, i;
  Uint1       strand = Seq_strand_plus;

  if (sip1 == NULL || sip2 == NULL) {
    return NULL;
  }

  salp = NWAlignmentForInterval(sip1, sip2, start1, stop1, start2, stop2);
  if (salp != NULL) {
    return salp;
  }

  dsp = DenseSegNew ();
  dsp->dim = 2;
  dsp->ids = SeqIdDup (sip1);
  dsp->ids->next = SeqIdDup (sip2);
  len1 = stop1 - start1 + 1;

  if (stop2 > start2) {
    len2 = stop2 - start2 + 1;
  } else {
    len2 = start2 - stop2 + 1;
    strand = Seq_strand_minus;
  }

  if (len1 == len2) {
    dsp->numseg = 1;
  } else {
    dsp->numseg = 2;
  }

  dsp->starts = (Int4Ptr) MemNew (sizeof (Int4) * dsp->dim * dsp->numseg);
  dsp->strands = (Uint1Ptr) MemNew (sizeof (Uint1) * dsp->dim * dsp->numseg);
  dsp->lens = (Int4Ptr) MemNew (sizeof (Int4) * dsp->dim);
  
  if (len1 == len2) {
    dsp->lens[0] = len1;
  } else if (len1 > len2) {
    dsp->lens[0] = len2;
    dsp->lens[1] = len1 - len2;
  } else {
    dsp->lens[0] = len1;
    dsp->lens[1] = len2 - len1;
  }

  dsp->starts[0] = start1;
  if (strand == Seq_strand_minus) {
    dsp->starts[1] = stop2;
  } else {
    dsp->starts[1] = start2;
  }  

  for (i = 0; i < dsp->numseg; i++) {
    dsp->strands[2 * i] = Seq_strand_plus;
    dsp->strands[2 * i + 1] = strand;
  }

  if (dsp->numseg > 1) {
    if (len1 > len2) {
      dsp->starts[2] = dsp->starts[0] + dsp->lens[0];
      dsp->starts[3] = -1;
    } else {
      dsp->starts[2] = -1;
      if (strand == Seq_strand_minus) {
        dsp->starts[3] = dsp->starts[1] + dsp->lens[0] + dsp->lens[1] - 1;
      } else {
        dsp->starts[3] = dsp->starts[1] + dsp->lens[0];
      }
    }
  }

  salp = SeqAlignNew ();
  salp->segtype = SAS_DENSEG;
  salp->segs = dsp;
  
  AlnMgr2IndexSingleChildSeqAlign (salp);  
  return salp;
}


static void ReportCreatedAlignment (LogInfoPtr lip, SeqAlignPtr salp)
{
  Int4        from_1, to_1, from_2, to_2;
  Uint1       strand;
  Char        id1[255], id2[255];
  SeqIdPtr    sip;
  BioseqPtr   bsp;

  if (lip == NULL || lip->fp == NULL || salp == NULL) {
    return;
  }

  sip = AlnMgr2GetNthSeqIdPtr (salp, 1);
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sip = SeqIdFree (sip);
    sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
    BioseqUnlock (bsp);
  }

  SeqIdWrite (sip, id1, PRINTID_REPORT, sizeof (id1) - 1);
  sip = SeqIdFree (sip);
  sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sip = SeqIdFree (sip);
    sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
    BioseqUnlock (bsp);
  }
  SeqIdWrite (sip, id2, PRINTID_REPORT, sizeof (id2) - 1);
  sip = SeqIdFree (sip);
  
  strand = SeqAlignStrand (salp, 2);
  AlnMgr2GetNthSeqRangeInSA (salp, 1, &from_1, &to_1);
  AlnMgr2GetNthSeqRangeInSA (salp, 2, &from_2, &to_2);
  fprintf (lip->fp, "Created alignment to cover space between local alignments: %s:%d-%d, %s:%d-%d%s\n",
                     id1, from_1, to_1, id2, from_2, to_2,
                     strand == Seq_strand_minus ? "(c)" : "");

  WriteAlignmentInterleaveToFileEx (salp, lip->fp, 40, FALSE, TRUE);
}


static void FillInAlignmentHoles (SeqAlignPtr salp_list, LogInfoPtr lip)
{
  SeqAlignPtr salp, salp_new;
  Int4        start1_this, start1_next;
  Int4        stop1_this, stop1_next;
  Int4        start2_this, start2_next;
  Int4        stop2_this, stop2_next;
  Uint1       strand;
  SeqIdPtr    sip1, sip2;

  if (salp_list == NULL || salp_list->next == NULL) {
    return;
  }

  sip1 = AlnMgr2GetNthSeqIdPtr (salp_list, 1);
  sip2 = AlnMgr2GetNthSeqIdPtr (salp_list, 2);

  /* note - unlike the other functions, SeqAlignStrand uses 0-based index */
  strand = SeqAlignStrand (salp_list, 1);

  salp = salp_list;
  AlnMgr2GetNthSeqRangeInSA (salp, 1, &start1_this, &stop1_this);
  AlnMgr2GetNthSeqRangeInSA (salp, 2, &start2_this, &stop2_this);
  
  while (salp->next != NULL) {
    AlnMgr2GetNthSeqRangeInSA (salp->next, 1, &start1_next, &stop1_next);
    AlnMgr2GetNthSeqRangeInSA (salp->next, 2, &start2_next, &stop2_next);
    if (start1_next > stop1_this + 1
        || (strand == Seq_strand_minus && start2_next < stop2_this - 1)
        || (strand != Seq_strand_minus && start2_next > stop2_this + 1)) {
      if (strand == Seq_strand_minus) {
        salp_new = AlignmentForInterval (sip1, sip2, stop1_this + 1, start1_next - 1, stop2_this - 1, start2_next + 1);
      } else {
        salp_new = AlignmentForInterval (sip1, sip2, stop1_this + 1, start1_next - 1, stop2_this + 1, start2_next - 1);
      }
      ReportCreatedAlignment (lip, salp_new);
      salp_new->next = salp->next;
      salp->next = salp_new;
      salp = salp_new;
    }
    start1_this = start1_next;
    stop1_this = stop1_next;
    start2_this = start2_next;
    stop2_this = stop2_next;
    salp = salp->next;
  }

}


static SeqAlignPtr MergeAlignments (SeqAlignPtr salp_list)
{
  SeqAlignPtr salp_new = NULL, salp, salp_next;
  DenseSegPtr dsp, dsp_new;
  Int4        seg_num, k;

  if (salp_list == NULL || salp_list->next == NULL) {
    return salp_list;
  }

  dsp_new = DenseSegNew ();
  dsp_new->dim = 2;
  dsp_new->ids = AlnMgr2GetNthSeqIdPtr (salp_list, 1);
  dsp_new->ids->next = AlnMgr2GetNthSeqIdPtr (salp_list, 2);

  /* get total number of segments */
  for (salp = salp_list; salp != NULL; salp = salp->next) {
    dsp = (DenseSegPtr) salp->segs;
    dsp_new->numseg += dsp->numseg;
  }

  dsp_new->starts = (Int4Ptr) MemNew (sizeof (Int4) * dsp_new->dim * dsp_new->numseg);
  dsp_new->strands = (Uint1Ptr) MemNew (sizeof (Uint1) * dsp_new->dim * dsp_new->numseg);
  dsp_new->lens = (Int4Ptr) MemNew (sizeof (Int4) * dsp_new->numseg);
  
  seg_num = 0;
  for (salp = salp_list; salp != NULL; salp = salp_next) {
    salp_next = salp->next;
    dsp = (DenseSegPtr) salp->segs;
    for (k = 0; k < dsp->numseg; k++) {
      dsp_new->lens[seg_num] = dsp->lens[k];
      dsp_new->starts[2 * seg_num] = dsp->starts[2 * k];
      dsp_new->starts[2 * seg_num + 1] = dsp->starts[2 * k + 1];
      dsp_new->strands[2 * seg_num] = dsp->strands[2 * k];
      dsp_new->strands[2 * seg_num + 1] = dsp->strands[2 * k + 1];
      seg_num++;
    }
    salp->next = NULL;
    salp = SeqAlignFree (salp);
  }

  salp_new = SeqAlignNew ();
  salp_new->segtype = SAS_DENSEG;
  salp_new->segs = dsp_new;
  salp_new->dim = 2;

  return salp_new;
}


static Boolean MatchWithAmbiguity (Char ch1, Char ch2)
{
  Boolean rval = FALSE;

  ch1 = toupper (ch1);
  ch2 = toupper (ch2);

  if (ch1 == ch2) {
    return TRUE;
  }
  if (ch1 == 'X' || ch2 == 'X') {
    return TRUE;
  }
  switch (ch1) {
    case 'A':
      if (ch2 == 'M' || ch2 == 'R' || ch2 == 'W' || ch2 == 'V' || ch2 == 'H' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'T':
      if (ch2 == 'K' || ch2 == 'Y' || ch2 == 'W' || ch2 == 'B' || ch2 == 'H' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'G':
      if (ch2 == 'K' || ch2 == 'R' || ch2 == 'S' || ch2 == 'B' || ch2 == 'V' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'C':
      if (ch2 == 'M' || ch2 == 'Y' || ch2 == 'S' || ch2 == 'B' || ch2 == 'V' || ch2 == 'H') {
        rval = TRUE;
      }
      break;
    case 'K':
      if (ch2 == 'G' || ch2 == 'T' || ch2 == 'B' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'M':
      if (ch2 == 'A' || ch2 == 'C' || ch2 == 'V' || ch2 == 'H') {
        rval = TRUE;
      }
      break;
    case 'R':
      if (ch2 == 'A' || ch2 == 'G' || ch2 == 'V' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'Y':
      if (ch2 == 'C' || ch2 == 'T' || ch2 == 'B' || ch2 == 'H') {
        rval = TRUE;
      }
      break;
    case 'S':
      if (ch2 == 'C' || ch2 == 'G' || ch2 == 'B' || ch2 == 'V') {
        rval = TRUE;
      }
      break;
    case 'W':
      if (ch2 == 'A' || ch2 == 'T' || ch2 == 'H' || ch2 == 'D') {
        rval = TRUE;
      }
      break;
    case 'B':
      if (ch2 == 'C' || ch2 == 'G' || ch2 == 'T' || ch2 == 'K' || ch2 == 'Y' || ch2 == 'S') {
        rval = TRUE;
      }
      break;
    case 'V':
      if (ch2 == 'A' || ch2 == 'C' || ch2 == 'G' || ch2 == 'M' || ch2 == 'R' || ch2 == 'S') {
        rval = TRUE;
      }
      break;
    case 'H':
      if (ch2 == 'A' || ch2 == 'C' || ch2 == 'T' || ch2 == 'M' || ch2 == 'Y' || ch2 == 'W') {
        rval = TRUE;
      }
      break;
    case 'D':
      if (ch2 == 'A' || ch2 == 'G' || ch2 == 'T' || ch2 == 'K' || ch2 == 'R' || ch2 == 'W') {
        rval = TRUE;
      }
      break;
  }
  return rval;
}  
  

/* expand for ambiguity characters and poly-A tail */
static SeqAlignPtr ExtendAlignmentList (SeqAlignPtr salp_list)
{
  Int4        from_1, to_1, from_2, to_2, len1, len2, len_check, len_extend;
  SeqAlignPtr salp_tmp;
  BioseqPtr   bsp1 = NULL;
  BioseqPtr   bsp2 = NULL;
  Uint1       strand;
  CharPtr     buf1, buf2;
  Int2        ctr;
  SeqIdPtr    sip1, sip2;
  DenseSegPtr dsp;

  if (salp_list == NULL) {
    return salp_list;
  }

  strand = SeqAlignStrand (salp_list, 1);
  AlnMgr2IndexSingleChildSeqAlign (salp_list);

  sip1 = AlnMgr2GetNthSeqIdPtr (salp_list, 1);
  bsp1 = BioseqLockById (sip1);
  sip2 = AlnMgr2GetNthSeqIdPtr (salp_list, 2);
  bsp2 = BioseqLockById (sip2);

  AlnMgr2GetNthSeqRangeInSA (salp_list, 1, &from_1, &to_1);  
  AlnMgr2GetNthSeqRangeInSA (salp_list, 2, &from_2, &to_2);  

  if (from_1 > 0 
      && ((strand == Seq_strand_plus && from_2 > 0)
          || (strand == Seq_strand_minus && to_2 < bsp2->length - 1))) {
    len1 = from_1;
    if (strand == Seq_strand_plus) {
      len2 = from_2;
    } else {
      len2 = bsp2->length - to_2 - 1;
    }
    if (len1 > len2) {
      len_check = len2;
    } else {
      len_check = len1;
    }
    buf1 = (CharPtr) MemNew (sizeof (Char) * (len_check  + 1));
    buf2 = (CharPtr) MemNew (sizeof (Char) * (len_check  + 1));
    ctr = SeqPortStreamInt (bsp1, from_1 - len_check, from_1 - 1, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf1, NULL);
    buf1[ctr] = 0;
    if (strand == Seq_strand_plus) {
      ctr = SeqPortStreamInt (bsp2, from_2 - len_check, from_2 - 1, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);
    } else {
      ctr = SeqPortStreamInt (bsp2, to_2 + 1, to_2 + len_check, Seq_strand_minus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);
    }
    buf2[ctr] = 0;

    len_extend = 0;
    while (len_extend < len_check 
           && MatchWithAmbiguity (buf1[len_check - len_extend - 1], buf2[len_check - len_extend - 1])) {
      len_extend++;
    }
    buf1 = MemFree (buf1);
    buf2 = MemFree (buf2);
    if (len_extend > 0) {
      dsp = (DenseSegPtr) salp_list->segs;
      dsp->lens[0] += len_extend;
      dsp->starts[0] -= len_extend;
      if (strand == Seq_strand_plus) {
        dsp->starts[1] -= len_extend;
      }
      SeqAlignIndexFree(salp_list->saip);
      salp_list->saip = NULL;
      AlnMgr2IndexSingleChildSeqAlign (salp_list);
    }
  }

  /* extend at other end */
  salp_tmp = salp_list;
  while (salp_tmp->next != NULL) {
    salp_tmp = salp_tmp->next;
  }

  AlnMgr2IndexSingleChildSeqAlign (salp_tmp);

  AlnMgr2GetNthSeqRangeInSA (salp_tmp, 1, &from_1, &to_1);  
  AlnMgr2GetNthSeqRangeInSA (salp_tmp, 2, &from_2, &to_2);
  if (to_1 < bsp1->length - 1
      && ((strand == Seq_strand_plus && to_2 < bsp2->length - 1)
          || (strand == Seq_strand_minus && from_2 > 0))) {
    len1 = bsp1->length - to_1 - 1;
    if (strand == Seq_strand_plus) {
      len2 = bsp2->length - to_2 - 1;
    } else {
      len2 = from_2;
    }
    if (len1 > len2) {
      len_check = len2;
    } else {
      len_check = len1;
    }
    if (len_check > 0) {
      buf1 = (CharPtr) MemNew (sizeof (Char) * (len_check  + 1));
      buf2 = (CharPtr) MemNew (sizeof (Char) * (len_check  + 1));
      ctr = SeqPortStreamInt (bsp1, to_1 + 1, to_1 + len_check, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf1, NULL);
      buf1[ctr] = 0;
      if (strand == Seq_strand_plus) {
        ctr = SeqPortStreamInt (bsp2, to_2 + 1, to_2 + len_check, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);
      } else {
        ctr = SeqPortStreamInt (bsp2, from_2 - len_check, from_2 - 1, Seq_strand_minus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);
      }
      buf2[ctr] = 0;

      len_extend = 0;
      while (len_extend < len_check 
            && MatchWithAmbiguity (buf1[len_extend], buf2[len_extend])) {
        len_extend++;
      }
      buf1 = MemFree (buf1);
      buf2 = MemFree (buf2);
      if (len_extend > 0) {
        dsp = (DenseSegPtr) salp_tmp->segs;
        dsp->lens[dsp->numseg - 1] += len_extend;
        if (strand == Seq_strand_minus) {
          dsp->starts[(dsp->numseg - 1) * dsp->dim + 1] -= len_extend;
        }
        SeqAlignIndexFree(salp_tmp->saip);
        salp_tmp->saip = NULL;
        AlnMgr2IndexSingleChildSeqAlign (salp_tmp);
      }
    }
  }   

  BioseqUnlock (bsp1);
  BioseqUnlock (bsp2);
  sip1 = SeqIdFree (sip1);  
  sip2 = SeqIdFree (sip2);  

  return salp_list;
}


static void ReportInitialBlastResults (LogInfoPtr lip, SeqAlignPtr salp_list)
{
  Int4        from_1, to_1, from_2, to_2;
  Uint1       strand;
  Char        id1[255], id2[255];
  SeqIdPtr    sip;
  BioseqPtr   bsp;

  if (lip == NULL || lip->fp == NULL || salp_list == NULL) {
    return;
  }

  AlnMgr2IndexSingleChildSeqAlign (salp_list);
  sip = AlnMgr2GetNthSeqIdPtr (salp_list, 1);
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sip = SeqIdFree (sip);
    sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
    BioseqUnlock (bsp);
  }
  SeqIdWrite (sip, id1, PRINTID_REPORT, sizeof (id1) - 1);
  sip = SeqIdFree (sip);
  sip = AlnMgr2GetNthSeqIdPtr (salp_list, 2);
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sip = SeqIdFree (sip);
    sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
    BioseqUnlock (bsp);
  }
  SeqIdWrite (sip, id2, PRINTID_REPORT, sizeof (id2) - 1);
  sip = SeqIdFree (sip);
  fprintf (lip->fp, "Initial BLAST results\n");
  while (salp_list != NULL) {
    AlnMgr2IndexSingleChildSeqAlign (salp_list);
    strand = SeqAlignStrand (salp_list, 2);
    AlnMgr2GetNthSeqRangeInSA (salp_list, 1, &from_1, &to_1);
    AlnMgr2GetNthSeqRangeInSA (salp_list, 2, &from_2, &to_2);
    fprintf (lip->fp, "%s:%d-%d, %s:%d-%d%s\n", id1, from_1, to_1, id2, from_2, to_2,
             strand == Seq_strand_minus ? "(c)" : "");
    salp_list = salp_list->next;
  }
  lip->data_in_log = TRUE;
}


static void ReportForRemoval (LogInfoPtr lip, SeqAlignPtr salp, CharPtr reason)
{
  Int4        from_1, to_1, from_2, to_2;
  Uint1       strand;
  Char        id1[255], id2[255];
  SeqIdPtr    sip;

  if (lip == NULL || lip->fp == NULL || salp == NULL) {
    return;
  }

  sip = AlnMgr2GetNthSeqIdPtr (salp, 1);
  SeqIdWrite (sip, id1, PRINTID_REPORT, sizeof (id1) - 1);
  sip = SeqIdFree (sip);
  sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
  SeqIdWrite (sip, id2, PRINTID_REPORT, sizeof (id2) - 1);
  sip = SeqIdFree (sip);
  AlnMgr2IndexSingleChildSeqAlign (salp);
  strand = SeqAlignStrand (salp, 2);
  AlnMgr2GetNthSeqRangeInSA (salp, 1, &from_1, &to_1);
  AlnMgr2GetNthSeqRangeInSA (salp, 2, &from_2, &to_2);
  fprintf (lip->fp, "Removed alignment %s:%d-%d, %s:%d-%d%s %s\n",
                     id1, from_1, to_1, id2, from_2, to_2,
                     strand == Seq_strand_minus ? "(c)" : "", reason);
}


/* Assume earlier alignments are better.  
 * Remove all alignments that are not the same strand as the first.
 * Remove all alignments that overlap on the first sequence.
 * Remove all alignments that overlap on the second sequence.
 * use Needleman-Wunsch to fill in holes between alignments.
 */
static SeqAlignPtr CombineTSAAlignments (SeqAlignPtr salp)
{
  SeqAlignPtr salp_tmp, salp_match, salp_prev, salp_next;
  Uint1       strand;
  Int4        from_1, to_1, from_2, to_2;
  LogInfoPtr  lip = NULL;

  if (salp == NULL) {
    return salp;
  }

  for (salp_tmp = salp; salp_tmp != NULL; salp_tmp = salp_tmp->next) {
    /* if alignment is to minus strand for first sequence, flip alignment */
    if (SeqAlignStrand (salp_tmp, 0) != Seq_strand_plus) {
      FlipAlignment (salp_tmp);
    }
  }

  if (salp->next != NULL) {
    lip = OpenLog ("TSA Alignment Adjustments");
  }

  ReportInitialBlastResults (lip, salp);

  if (salp->next != NULL) {
    /* remove alignments that are not the same strand as the initial alignment */
    strand = SeqAlignStrand (salp, 1);
    salp_prev = salp;
    AlnMgr2IndexSingleChildSeqAlign (salp);
    salp_tmp = salp->next;
    while (salp_tmp != NULL) {
      AlnMgr2IndexSingleChildSeqAlign (salp_tmp);
      salp_next = salp_tmp->next;
      if (SeqAlignStrand (salp_tmp, 1) != strand) {
        ReportForRemoval (lip, salp_tmp, "because strands do not match.");
        salp_prev->next = salp_tmp->next;
        salp_tmp->next = NULL;
        salp_tmp = SeqAlignFree (salp_tmp);
      }
      salp_tmp = salp_next;
    }
  }
  
  /* remove alignments that overlap on the first sequence */
  salp_match = salp;
  while (salp_match != NULL) {
    salp_prev = salp_match;
    salp_tmp = salp_match->next;
    AlnMgr2GetNthSeqRangeInSA (salp_match, 1, &from_1, &to_1);
    while (salp_tmp != NULL) {
      salp_next = salp_tmp->next;
      AlnMgr2GetNthSeqRangeInSA (salp_tmp, 1, &from_2, &to_2);
      if ((from_2 >= from_1 && from_2 <= to_1) || (to_2 >= from_1 && to_2 <= to_1)) { 
        ReportForRemoval (lip, salp_tmp, "because alignments overlap for first sequence.");
        salp_prev->next = salp_tmp->next;
        salp_tmp->next = NULL;
        salp_tmp = SeqAlignFree (salp_tmp);
      }
      salp_tmp = salp_next;
    }
    salp_match = salp_match->next;
  }

  if (salp->next != NULL) {
    /* remove alignments that overlap on the second sequence */
    salp_match = salp;
    while (salp_match != NULL) {
      salp_prev = salp_match;
      salp_tmp = salp_match->next;
      AlnMgr2GetNthSeqRangeInSA (salp_match, 2, &from_1, &to_1);
      while (salp_tmp != NULL) {
        salp_next = salp_tmp->next;
        AlnMgr2GetNthSeqRangeInSA (salp_tmp, 2, &from_2, &to_2);
        if ((from_2 >= from_1 && from_2 <= to_1) || (to_2 >= from_1 && to_2 <= to_1)) { 
          ReportForRemoval (lip, salp_tmp, "because alignments overlap for second sequence.");
          salp_prev->next = salp_tmp->next;
          salp_tmp->next = NULL;
          salp_tmp = SeqAlignFree (salp_tmp);
        }
        salp_tmp = salp_next;
      }
      salp_match = salp_match->next;
    }
  }

  /* sort remaining alignments by start position on the first sequence */
  salp = SortPairwiseAlignmentsByFirstSeqRange (salp);

  /* temporary hack.  only interested in "unusual" errors. */
  if (lip != NULL) {
    lip->data_in_log = FALSE;
  }

  if (salp->next != NULL) {
    /* remove alignments that are out of order on the second sequence */
    salp_prev = salp;
    AlnMgr2GetNthSeqRangeInSA (salp_prev, 2, &from_1, &to_1);
    salp_tmp = salp->next;
    while (salp_tmp != NULL) {
      AlnMgr2GetNthSeqRangeInSA (salp_tmp, 2, &from_2, &to_2);
      if (from_2 < from_1) {
        ReportForRemoval (lip, salp_tmp, "because alignments are out of order for second sequence.");
        salp_prev->next = salp_tmp->next;
        salp_tmp->next = NULL;
        salp_tmp = SeqAlignFree (salp_tmp);
      } else {
        salp_prev = salp_tmp;
        from_1 = from_2;
        to_1 = to_2;
      }
      salp_tmp = salp_prev->next;
    }
  }

  /* fill in holes */
  FillInAlignmentHoles (salp, lip);

  /* extend for good matches */
  salp = ExtendAlignmentList (salp);

  /* make new alignment by stringing together local alignments */
  salp = MergeAlignments (salp);
  CloseLog (lip);
  lip = FreeLog (lip);
  return salp;
}


NLM_EXTERN SeqAlignPtr GetSeqAlignTSA (BioseqPtr bsp1, BioseqPtr bsp2)
{
   BLAST_SummaryOptions *options = NULL;
   SeqAlignPtr           salp = NULL;
   
   if (bsp1 == NULL || bsp2 == NULL) return NULL;

   BLAST_SummaryOptionsInit(&options);
   options->filter_string = StringSave ("m L");
   options->word_size = 20;
   options->cutoff_evalue = act_get_eval (1);
   options->hint = eNone;
   if (ISA_na (bsp1->mol))
   {
     options->program = eBlastn;
   }
   else
   {
     options->program = eBlastp;
   }

   BLAST_TwoSequencesSearch(options, bsp1, bsp2, &salp);

   /* if there were no alignments from the first search, try again
    * with low-complexity masking turned off.
    */
   if (salp == NULL) {
     options->filter_string = MemFree (options->filter_string);
     options->filter_string = StringSave ("m F");
     BLAST_TwoSequencesSearch(options, bsp1, bsp2, &salp);
   }

   BLAST_SummaryOptionsFree(options);

   if (salp != NULL) {
     salp = CombineTSAAlignments (salp);
     AlnMgr2IndexSeqAlign(salp);
   }
   return salp;
}


/* code for editing TSA assembly */
typedef struct tsaassemblydialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        intervals;
  BioseqPtr     consensus_bsp;
} TSAAssemblyDialog, PNTR TSAAssemblyDialogPtr;


Uint2 tsa_assembly_types [] = {
  TAGLIST_TEXT, TAGLIST_PROMPT, TAGLIST_PROMPT
};


Uint2 tsa_assembly_widths [] = {
  16, 16, 16, 0
};


static void SetTSAAssemblyDialogConsensusBioseq (DialoG d, BioseqPtr bsp)
{
  TSAAssemblyDialogPtr dlg;

  dlg = (TSAAssemblyDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    dlg->consensus_bsp = bsp;
  }
}


static void SeqAlignsToTSAAssemblyDialog (DialoG d, Pointer data)
{
  TSAAssemblyDialogPtr  dlg;
  TagListPtr            tlp;
  SeqAlignPtr           salp_list, salp, salp_tmp;
  ValNodePtr            data_list = NULL;
  Int4                  tsa_from, tsa_to, primary_from, primary_to;
  SeqIdPtr              sip;
  Char                  id_buf[255];
  CharPtr               str;
  CharPtr               fmt = "%s\t%d-%d%s\t%d-%d\n";
  Uint1                 strand;
  BioseqPtr             bsp;

  dlg = (TSAAssemblyDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  tlp = (TagListPtr) GetObjectExtra (dlg->intervals);
  if (tlp == NULL) {
    return;
  }

  salp_list = (SeqAlignPtr) data;

  salp = salp_list;
  while (salp != NULL) {
    salp_tmp = salp->next;
    salp->next = NULL;
    AlnMgr2IndexSeqAlign (salp);
    AlnMgr2GetNthSeqRangeInSA(salp, 1, &tsa_from, &tsa_to);
    AlnMgr2GetNthSeqRangeInSA(salp, 2, &primary_from, &primary_to);
    sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
    bsp = BioseqLockById (sip);
    if (bsp != NULL) {
      sip = SeqIdFree (sip);
      sip = SeqIdDup (SeqIdFindBest (bsp->id, SEQID_GENBANK));
      BioseqUnlock (bsp);
    }
    SeqIdWrite (sip, id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
    sip = SeqIdFree (sip);
    strand = AlnMgr2GetNthStrand (salp, 2);
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (id_buf) + 70));
    sprintf (str, fmt, id_buf, primary_from + 1, primary_to + 1, strand == Seq_strand_minus ? "(c)" : "", tsa_from + 1, tsa_to + 1);
    ValNodeAddPointer (&data_list, 0, str);
    salp->next = salp_tmp;
    salp = salp->next;
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = data_list;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (ValNodeLen (data_list) - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
  } else {
    SafeHide (tlp->bar);
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);
}


static void TSATableCallback (Pointer data)
{
  ProcessExternalEvent ();
}


static Pointer TSAAssemblyDialogToTranscriptomeIdsList (DialoG d)
{
  TSAAssemblyDialogPtr  dlg;
  TagListPtr            tlp;
  ValNodePtr            vnp;
  CharPtr               txt;
  TranscriptomeIdsPtr   t;
  ValNodePtr            token_list = NULL;

  dlg = (TSAAssemblyDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  tlp = (TagListPtr) GetObjectExtra (dlg->intervals);
  if (tlp == NULL) {
    return NULL;
  }

  if (dlg->consensus_bsp == NULL) {
    return NULL;
  }

  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    txt = ExtractTagListColumn (vnp->data.ptrvalue, 0);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    } else {
      ValNodeAddPointer (&token_list, 0, txt);
    }
  }
  t = TranscriptomeIdsNew (dlg->consensus_bsp, token_list);
  vnp = ValNodeNew (NULL);
  vnp->choice = 0;
  vnp->data.ptrvalue = t;
  return vnp;
}


static DialoG CreateTSAAssemblyDialog (GrouP g)

{
  TSAAssemblyDialogPtr  dlg;
  GrouP                  p;
  GrouP                  x;
  GrouP                  y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (TSAAssemblyDialogPtr) MemNew (sizeof (TSAAssemblyDialog));
  if (dlg == NULL) return NULL;

  SetObjectExtra (p, dlg, NULL);
  dlg->dialog = (DialoG) p;
  dlg->todialog = SeqAlignsToTSAAssemblyDialog;
  dlg->fromdialog = TSAAssemblyDialogToTranscriptomeIdsList;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  /* first line */
  StaticPrompt (y, "PrimaryID", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Primary Interval", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "TSA Interval", 16 * stdCharWidth, 0, programFont, 'c');
  dlg->intervals = CreateTagListDialogEx3 (x, 10, 3, 1, tsa_assembly_types, tsa_assembly_widths, NULL, TRUE, FALSE, NULL, NULL, 
                                           NULL, NULL, FALSE, TRUE);

  return (DialoG) p;
}


static void ShowAssemblyAlignment (SeqAlignPtr salp)
{
  LogInfoPtr lip;

  if (salp == NULL) {
    return;
  }

  lip = OpenLog ("Assembly Alignments");
  while (salp != NULL) {
    WriteAlignmentInterleaveToFileEx (salp, lip->fp, 40, FALSE, TRUE);
    lip->data_in_log = TRUE;
    salp = salp->next;
  }
  CloseLog (lip);
  lip = FreeLog (lip);
}


typedef struct tsaassemblyform {
  FORM_MESSAGE_BLOCK
  DialoG intervals;

} TSAAssemblyFormData, PNTR TSAAssemblyFormPtr;


static void AcceptTSAAssembly (ButtoN b)
{
  TSAAssemblyFormPtr frm;
  BioseqPtr          bsp;
  ValNodePtr         err_list, coverage_report, vnp, ids_list, match_errs;
  SeqAlignPtr        salp, salp_next;
  LogInfoPtr         lip;

  frm = (TSAAssemblyFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }
  bsp = GetBioseqGivenIDs (frm->input_entityID, frm->input_itemID, frm->input_itemtype);
  if (bsp == NULL) {
    return;
  }

  WatchCursor();
  Update();
  SetTSAAssemblyDialogConsensusBioseq (frm->intervals, bsp);
  
  ids_list = DialogToPointer (frm->intervals);

  /* remove existing assembly */
  if (bsp->hist != NULL) {
    salp = bsp->hist->assembly;
    bsp->hist->assembly = NULL;
    while (salp != NULL) {
      salp_next = salp->next;
      salp->next = NULL;
      salp = SeqAlignFree (salp);
      salp = salp_next;
    }
  }

  if (ids_list == NULL) {
    Message (MSG_ERROR, "TSA assembly removed");
  } else {
    err_list = MakeTranscriptomeAssemblySeqHist (ids_list->data.ptrvalue, GetSeqAlignTSA, TSATableCallback, NULL);
    coverage_report = ReportCoverageForTranscriptomeIdsListSeqHist (ids_list);
    match_errs = ReportConsensusMatchForBioseqSeqHist (bsp);
    ValNodeLink (&coverage_report, match_errs);
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
  }

  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  Remove (frm->form);
  ArrowCursor();
  Update();
}


NLM_EXTERN void EditTSAAssembly (IteM i)
{
  BaseFormPtr        bfp;
  WindoW             w;
  TSAAssemblyFormPtr frm;
  BioseqPtr          bsp;
  GrouP              h, c;
  ButtoN             b;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) {
    Message (MSG_ERROR, "Must select sequence for editing TSA Assembly");
    return;
  }
  frm = (TSAAssemblyFormPtr) MemNew (sizeof (TSAAssemblyFormData));
  if (frm == NULL) return;
  frm->input_entityID = bfp->input_entityID;
  frm->input_itemID = bfp->input_itemID;
  frm->input_itemtype = bfp->input_itemtype;

  w = FixedWindow (-50, -33, -10, -10, "TSA Assembly", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->form = (ForM) w;
  h = HiddenGroup (w, -1, 0, NULL);
  frm->intervals = CreateTSAAssemblyDialog (h);

  if (bsp->hist == NULL) {
    PointerToDialog (frm->intervals, NULL);
  } else {
    PointerToDialog (frm->intervals, bsp->hist->assembly);
  }


  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", AcceptTSAAssembly);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->intervals, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();    

}


/* automatic defline generator */

typedef struct deffeats {
  SeqFeatPtr  sfp;
  SeqFeatPtr  gene;
  SeqFeatPtr  prot;
  CharPtr     genename;
  CharPtr     allelename;
  CharPtr     protname;
  Boolean     alreadyTrimmed;
  Uint2       entityID;
  Uint4       itemID;
  Uint2       subtype;
  Boolean     isDNA;
  Boolean     isAlleleGroup;
  Boolean     lastInString;
  Boolean     lastInGroup;
  Boolean     lastInType;
  Boolean     lastInPenultimate;
  Boolean     pseudo;
  Boolean     ignore;
  Boolean     suppressprefix;
  Int2        altSplices;
  Int2        numUnknown;
} DefFeatsData, PNTR DefFeatsPtr;

static Boolean GetMolBioFeatsGatherFunc (GatherContextPtr gcp, Boolean getGene, Boolean getSnoRNA)

{
  DefFeatsPtr  dfp;
  RnaRefPtr    rrp;
  SeqFeatPtr   sfp;
  CharPtr      str;
  Uint1        type;
  ValNodePtr   PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL)
    return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      if (getGene) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = FEATDEF_GENE;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    case SEQFEAT_CDREGION :
      dfp = MemNew (sizeof (DefFeatsData));
      if (dfp == NULL) return TRUE;
      dfp->entityID = gcp->entityID;
      dfp->itemID = gcp->itemID;
      dfp->sfp = sfp;
      dfp->subtype = FEATDEF_CDS;
      dfp->altSplices = 1;
      ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return TRUE;
      switch (rrp->type) {
        case 3 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_tRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 4 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_rRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
         case 5 :
          if (getSnoRNA) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_snRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
       case 7 :
          if (getSnoRNA) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_snoRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
        case 255 :
          if (rrp->ext.choice == 1) {
            str = (CharPtr) rrp->ext.value.ptrvalue;
            if (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                StringICmp (str, "internal transcribed spacer 2") == 0 ||
                StringICmp (str, "internal transcribed spacer 3") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                StringICmp (str, "ITS1") == 0 ||
                StringICmp (str, "ITS2") == 0 ||
                StringICmp (str, "ITS3") == 0) {
              dfp = MemNew (sizeof (DefFeatsData));
              if (dfp == NULL) return TRUE;
              dfp->entityID = gcp->entityID;
              dfp->itemID = gcp->itemID;
              dfp->sfp = sfp;
              dfp->subtype = FEATDEF_otherRNA;
              ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
            }
          }
          break;
        default :
          break;
      }
      break;
    case SEQFEAT_IMP :
      type = FindFeatDefType (sfp);
      if (type == FEATDEF_LTR || type == FEATDEF_exon) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = type;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetCDStRNArRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, FALSE, FALSE);
}


extern void StringToLower (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    *str = TO_LOWER (ch);
    str++;
    ch = *str;
  }
}

static CharPtr molinfo_tech_list [] = {
  "?", "standard", "EST", "STS", "survey", "genetic map", "physical map",
  "derived", "concept-trans", "seq-pept", "both", "seq-pept-overlap",
  "seq-pept-homol", "concept-trans-a", "htgs 1", "htgs 2", "htgs 3",
  "fli cDNA", "htgs 0", "htc", "wgs", "barcode", "composite-wgs-htgs", NULL
};


static void AddSubSourceValuesToNucTitle (BioSourcePtr biop, CharPtr str)
{
  CharPtr       ssp_name;
  SubSourcePtr  ssp;
  Char          text [256];

  if (biop == NULL || str == NULL) return;
  ssp = biop->subtype;
  while (ssp != NULL) {
    StringCpy (text, "[");
    ssp_name = GetSubsourceQualName (ssp->subtype);
    if (StringHasNoText (ssp_name)) {
      StringCat (text, "subsource");
    } else {
      StringCat (text, ssp_name);
    }
    StringToLower (text);
    StringCat (text, "=");
    StringCat (text, ssp->name);
    StringCat (text, "] ");
    StringCat (str, text);
    ssp = ssp->next;
  }
}

static void AddOrgModValuesToNucTitle (BioSourcePtr biop, CharPtr str)
{
  CharPtr    mod_name;
  OrgModPtr  mod;
  Char       text [256];

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || str == NULL) return;
  mod = biop->org->orgname->mod;
  while (mod != NULL) {
    StringCpy (text, "[");
    mod_name = GetOrgModQualName (mod->subtype);
    StringCat (text, mod_name);
    StringToLower (text);
    StringCat (text, "=");
    StringCat (text, mod->subname);
    StringCat (text, "] ");
    StringCat (str, text);
    mod = mod->next;
  }
}

static void MakeNucleotideTitlesInSequinStyle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqContextPtr   bcp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  MolInfoPtr         mip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            str;
  Uint1              tech = 0;
  Char               text [256];
  ValNodePtr         ttl;
  ValNodePtr         vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
  BioseqContextFree (bcp);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  if (bsp->descr != NULL) {
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
    vnp = ValNodeFreeData (vnp);
  }
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->tech) {
        case MI_TECH_est :
        case MI_TECH_sts :
        case MI_TECH_survey :
        case MI_TECH_htgs_1 :
        case MI_TECH_htgs_2 :
        case MI_TECH_htgs_3 :
        case MI_TECH_fli_cdna :
        case MI_TECH_htgs_0 :
        case MI_TECH_htc :
        case MI_TECH_wgs :
          tech = mip->tech;
          break;
        default :
          break;
      }
    }
  }
  BioseqContextFree (bcp);
  str = MemNew (2000);

  orp = biop->org;
  if (orp != NULL) {
    StringCpy (text, "[organism=");
    StringCat (text, orp->taxname);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  AddSubSourceValuesToNucTitle (biop, str);

  AddOrgModValuesToNucTitle (biop, str);

  if (tech > 0) {
    StringCpy (text, "[tech=");
    StringCat (text, molinfo_tech_list [tech]);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  TrimSpacesAroundString (str);
  if (! StringHasNoText (str)) {
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
    if (ttl != NULL) {
      ttl->data.ptrvalue = StringSave (str);
    }
  }
  MemFree (str);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data);

extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SeqEntryExplore (sep, NULL, MakeNucleotideTitlesInSequinStyle);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static Boolean StringHasEqualSignOrBrackets (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;
  ch = *str;
  while (ch != '\0') {
    if (ch == '=' || ch == '[' || ch == ']') return TRUE;
    str++;
    ch = *str;
  }
  return FALSE;
}

static void MakeProteinTitlesInSequinStyle (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  Char          quot [4];
  GeneRefPtr    grp;
  GatherScope   gs;
  Boolean       has_equal_or_brackets;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  Char          str [256];
  Char          text [256];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 ||
        (IsPopPhyEtcSet (bssp->_class))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeProteinTitlesInSequinStyle (entityID, sep);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetCDStRNArRNAGatherFunc, &gs);
  /* head = SortValNode (head, SortCDStRNArRNAByLocation); */
  if (head == NULL) return;

  quot [0] = '"';
  quot [1] = '\0';

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL && dfp->subtype == FEATDEF_CDS) {
      FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
      bsp = GetBioseqGivenSeqLoc (dfp->sfp->product, entityID);
      if (bsp != NULL) {
        str [0] = '\0';
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (text, (CharPtr) grp->locus, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[gene=");
              StringCat (str, text);
              StringCat (str, "]");
            }
            if (grp->syn != NULL) {
              StringNCpy_0 (text, (CharPtr) grp->syn->data.ptrvalue, sizeof (text));
              if (! StringHasNoText (text)) {
                if (str [0] != '\0') {
                  StringCat (str, " ");
                }
                StringCat (str, "[gene_syn=");
                StringCat (str, text);
                StringCat (str, "]");
              }
            }
          }
        }
        if (dfp->prot != NULL) {
          prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
          if (prp != NULL) {
            if (prp->name != NULL) {
              StringNCpy_0 (text, (CharPtr) prp->name->data.ptrvalue, sizeof (text));
              if (! StringHasNoText (text)) {
                if (str [0] != '\0') {
                  StringCat (str, " ");
                }
                StringCat (str, "[protein=");
                has_equal_or_brackets = StringHasEqualSignOrBrackets (text);
                if (has_equal_or_brackets) {
                  StringCat (str, quot);
                }
                StringCat (str, text);
                if (has_equal_or_brackets) {
                  StringCat (str, quot);
                }
                StringCat (str, "]");
              }
            }
            StringNCpy_0 (text, (CharPtr) prp->desc, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[prot_desc=");
              has_equal_or_brackets = StringHasEqualSignOrBrackets (text);
              if (has_equal_or_brackets) {
                StringCat (str, quot);
              }
              StringCat (str, text);
              if (has_equal_or_brackets) {
                StringCat (str, quot);
              }
              StringCat (str, "]");
            }
          }
        }
        if (! StringHasNoText (str)) {
          psep = SeqMgrGetSeqEntryForData (bsp);
          if (psep != NULL) {
            ttl = CreateNewDescriptor (psep, Seq_descr_title);
            if (ttl != NULL) {
              ttl->data.ptrvalue = StringSave (str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  MakeProteinTitlesInSequinStyle (ompcp->input_entityID, sep);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static Boolean LIBCALLBACK SequinFTableBioseq (BioseqPtr bsp, SeqMgrBioseqContextPtr context)

{
  FILE      *fp;

  if (bsp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  BioseqToGnbk (bsp, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  return TRUE;
}

extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data)

{
  FILE              *fp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
  SeqMgrExploreBioseqs (ompcp->input_entityID, NULL, (Pointer) fp, SequinFTableBioseq, TRUE, FALSE, FALSE);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Gene - CDS Feature Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

static void InsertGeneLocusTagPrefixCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CharPtr    prefix;
  GeneRefPtr grp;
  CharPtr    new_locus_tag = NULL;
  Int4       new_locus_tag_len = 0;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_GENE)
  {
    return;
  }
  
  prefix = (CharPtr) userdata;
  
  if (StringHasNoText (prefix))
  {
    return;
  }
  
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL)
  {
    grp = GeneRefNew();
    sfp->data.value.ptrvalue = grp;
  }
  if (grp == NULL)
  {
    return;
  }
  
  if (StringHasNoText (grp->locus_tag))
  {
    grp->locus_tag = MemFree (grp->locus_tag);
    grp->locus_tag = StringSave (prefix);
  }
  else
  {
    new_locus_tag_len = StringLen (prefix) + StringLen (grp->locus_tag) + 1;
    new_locus_tag = (CharPtr) MemNew (sizeof (Char) * new_locus_tag_len);
    if (new_locus_tag != NULL)
    {
      StringCpy (new_locus_tag, prefix);
      StringCat (new_locus_tag, grp->locus_tag);
      grp->locus_tag = MemFree (grp->locus_tag);
      grp->locus_tag = new_locus_tag;
    }
  }
} 

typedef struct locustagprefix
{
  FEATURE_FORM_BLOCK

  TexT prefix_txt;  
} LocusTagPrefixData, PNTR LocusTagPrefixPtr;

static void InsertGeneLocusTagPrefixButton (ButtoN b)
{
  LocusTagPrefixPtr ltpp;
  CharPtr           prefix;
  SeqEntryPtr       sep;
  
  ltpp = (LocusTagPrefixPtr) GetObjectExtra (b);
  if (ltpp == NULL)
  {
    return;
  }
  
  prefix = SaveStringFromText (ltpp->prefix_txt);
  if (!StringHasNoText (prefix))
  {
    sep = GetTopSeqEntryForEntityID (ltpp->input_entityID);
    if (sep != NULL)
    {
      VisitFeaturesInSep (sep, prefix, InsertGeneLocusTagPrefixCallback); 
    }
    
  }
  prefix = MemFree (prefix);
  ObjMgrSetDirtyFlag (ltpp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ltpp->input_entityID, 0, 0);
  Remove (ltpp->form);
  Update ();
}

extern void InsertGeneLocusTagPrefix (IteM i)
{
  BaseFormPtr        bfp;
  WindoW             w;
  LocusTagPrefixPtr  ltpp;
  GrouP              h, g, c;
  ButtoN             b;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  ltpp = (LocusTagPrefixPtr) MemNew (sizeof (LocusTagPrefixData));
  if (ltpp == NULL) return;
  ltpp->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Feature Evidence", StdCloseWindowProc);
  SetObjectExtra (w, ltpp, StdCleanupFormProc);
  ltpp->form = (ForM) w;
  h = HiddenGroup (w, -1, 0, NULL);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Prefix for Gene Locus Tag", 0, 0, programFont, 'c');
  ltpp->prefix_txt = DialogText (g, "", 14, NULL);
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", InsertGeneLocusTagPrefixButton);
  SetObjectExtra (b, ltpp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();    
}


static CharPtr MakeFlybaseTagString (CharPtr locus_tag)
{
  CharPtr new_str;
  CharPtr cp;
  
  if (locus_tag == NULL) return NULL;
  cp = locus_tag;
  while (*cp != 0 && !isdigit(*cp))
  {
    cp++;
  }
  if (StringLen (cp) > 7)
  {
    cp += StringLen (cp) - 7;
  }
  
  new_str = (CharPtr) MemNew (12 * sizeof (Char));
  if (new_str != NULL)
  {
    sprintf (new_str, "FBti0000000");
    StringCpy (new_str + 4 + 7 - StringLen (cp), cp);
  }
  return new_str;
}

static void 
ReplaceRepeatRegionLocusTagWithDbxrefCallback
(SeqFeatPtr sfp,
 Pointer userdata)
{
  GBQualPtr       gbqual, prev_qual = NULL, next_qual;
  ValNodePtr      vnp;
  Boolean         has_dbxref = FALSE;
  DbtagPtr        tag;
  CharPtr         new_string = NULL;
  CharPtr         new_comment;
	SeqFeatXrefPtr 	xrp, prev_xrp = NULL, next_xrp;
	GeneRefPtr      grp;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region)
  {
    return;
  }
  
	for (xrp = sfp->xref; xrp; xrp = next_xrp) 
	{
	  next_xrp = xrp->next;
		if (xrp->data.choice == SEQFEAT_GENE) 
		{
			grp = (GeneRefPtr) xrp->data.value.ptrvalue;
			if (grp != NULL && !StringHasNoText (grp->locus_tag))
			{
			  new_string = StringSave (grp->locus_tag);
			  if (prev_xrp == NULL)
			  {
			    sfp->xref = xrp->next;
			  }
			  else
			  {
			    prev_xrp->next = xrp->next;
			  }
			  xrp->next = NULL;
			  SeqFeatXrefFree (xrp);
			}
			else
			{
			  prev_xrp = xrp;
			}
		}
		else
		{
		  prev_xrp = xrp;
		}
	}
  
  for (gbqual = sfp->qual; 
       gbqual != NULL && new_string == NULL;
       gbqual = next_qual)
  {
    next_qual = gbqual->next;
    if (StringCmp (gbqual->qual, "locus_tag") == 0)
    {
      new_string = StringSave (gbqual->val);
      if (prev_qual == NULL)
      {
        sfp->qual = gbqual->next;
      }
      else
      {
        prev_qual->next = gbqual->next;
      }
      gbqual->next = NULL;
      GBQualFree (gbqual);
    }
    else
    {
      prev_qual = gbqual;
    }
  }
  if (new_string == NULL)
  {
    return;
  }
  
  for (vnp = sfp->dbxref; vnp != NULL && !has_dbxref; vnp = vnp->next)
  {
    tag = (DbtagPtr) vnp->data.ptrvalue;
    if (tag != NULL && StringCmp (tag->db, "FLYBASE") == 0)
    {
      has_dbxref = TRUE;
    }
  }
  if (!has_dbxref)
  {
    tag = DbtagNew ();
    if (tag != NULL)
    {
      tag->db = StringSave ("FLYBASE");
      tag->tag = ObjectIdNew ();
      tag->tag->str = MakeFlybaseTagString (new_string);
      ValNodeAddPointer (&(sfp->dbxref), 0, tag);
    }
  }
  
  if (StringHasNoText (sfp->comment))
  {
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = StringSave (new_string);
  }
  else
  {
    new_comment = (CharPtr) MemNew (sizeof (Char) * (StringLen (sfp->comment) + StringLen (new_string) + 2));
    if (new_comment != NULL)
    {
      StringCpy (new_comment, new_string);
      StringCat (new_comment, ";");
      StringCat (new_comment, sfp->comment);
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = new_comment;
    }
  }
  
  MemFree (new_string);
}

extern void ReplaceRepeatRegionLocusTagWithDbxref (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ReplaceRepeatRegionLocusTagWithDbxrefCallback);
    
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update (); 
}

static void SetPrimerBindPairStrandsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr        primer_1 = NULL, primer_2 = NULL, sfp;
  SeqMgrFeatContext context;
  Uint1             first_strand, second_strand;
  
  if (bsp == NULL)
  {
    return;
  }
  
  /* must have exactly two primer_bind features */
  primer_1 = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_primer_bind, &context);
  if (primer_1 == NULL)
  {
    return;
  }
  primer_2 = SeqMgrGetNextFeature (bsp, primer_1, 0, FEATDEF_primer_bind, &context);
  if (primer_2 == NULL)
  {
    return;
  }
  
  /* if there are three, must abandon */
  sfp = SeqMgrGetNextFeature (bsp, primer_1, 0, FEATDEF_primer_bind, &context);
  if (sfp != NULL)
  {
    return;
  }
  
  first_strand = SeqLocStrand (primer_1->location);
  second_strand = SeqLocStrand (primer_2->location);
  
  if (first_strand == Seq_strand_minus)
  {
    if (second_strand == Seq_strand_minus)
    {
      SetSeqLocStrand (primer_2->location, Seq_strand_plus);
    }
  }
  else
  {
    if (second_strand != Seq_strand_minus)
    {
      SetSeqLocStrand (primer_2->location, Seq_strand_minus);
    }
  }
}

extern void SetPrimerBindPairStrands (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, SetPrimerBindPairStrandsCallback);
    
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update (); 
}

/* This function returns the length change, so that the start
 * position for subsequent features can be adjusted.
 */
static Int4 
FindAndConvertGapFeat 
(BioseqPtr  bsp,
 Int4       start,
 Boolean    make_known)
{
  ValNodePtr vnp;
  SeqLocPtr  slp;
  SeqLitPtr  litp;
  Int4       currpos = 0, len_diff;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta
      || start < 0)
  {
    return 0;
  }
  for (vnp = (ValNodePtr)(bsp->seq_ext); 
       vnp != NULL && currpos < start; 
       vnp = vnp->next) 
  {
    if (vnp->choice == 1) 
    {
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      if (slp == NULL) continue;
      currpos += SeqLocLen (slp);
    }
    else if (vnp->choice == 2) 
    {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL) continue;
      currpos += litp->length;
    }
  }
  
  if (currpos < start || vnp == NULL || vnp->choice != 2)
  {
    return 0;
  }
  
  litp = (SeqLitPtr) vnp->data.ptrvalue;
  if (litp == NULL)
  {
    return 0;
  }
  
  if (make_known)
  {
    litp->fuzz = IntFuzzFree (litp->fuzz);
    len_diff = 0;
  }
  else if (litp->fuzz == NULL)
  {
    litp->fuzz = IntFuzzNew ();
    litp->fuzz->choice = 4;
    len_diff = litp->length - 100;
    litp->length = 100;
    bsp->length -= len_diff;
  }
  return len_diff;
}

typedef struct gapadjust 
{
  Int4 adjust_start;
  Int4 adjust_len;
} GapAdjustData, PNTR GapAdjustPtr;

static Int4 GetAdjustedGapStart (Int4 feat_left, ValNodePtr prev_adjust)
{
  GapAdjustPtr p;
  Int4         new_feat_left;
  
  new_feat_left = feat_left;
  while (prev_adjust != NULL)
  {
    p = (GapAdjustPtr) prev_adjust->data.ptrvalue;
    if (p != NULL && feat_left > p->adjust_start)
    {
      new_feat_left -= p->adjust_len;
    }
    prev_adjust = prev_adjust->next;
  }
  return new_feat_left;
}

static void ConvertSelectedGapFeatures (IteM i, Boolean to_known)
{
  BaseFormPtr       bfp;
  SelStructPtr      sel;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp;
  SeqFeatPtr        sfp;
  Int4              adjusted_start, len_diff;
  ValNodePtr        adjustment_list = NULL;
  GapAdjustPtr      p;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sel = ObjMgrGetSelected ();
  if (sel == NULL)
  {
    Message (MSG_ERROR, "Must select gap features to convert!");
    return;
  }
  WatchCursor ();
  Update ();
  while (sel != NULL)
  {
    if (sel->entityID == bfp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
    {
      sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_gap)
      {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL && bsp->repr == Seq_repr_delta)
        {
          adjusted_start = GetAdjustedGapStart (fcontext.left, adjustment_list);
          len_diff = FindAndConvertGapFeat (bsp, adjusted_start, to_known);
          if (len_diff != 0)
          {
            p = (GapAdjustPtr) MemNew (sizeof (GapAdjustData));
            if (p != NULL)
            {
              p->adjust_start = fcontext.left;
              p->adjust_len = len_diff;
              ValNodeAddPointer (&adjustment_list, 0, p);
            }
          }
        }
      }
    }
    sel = sel->next;
  }
  
  adjustment_list = ValNodeFreeData (adjustment_list);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


typedef struct convertgaptounknown
{
  FEATURE_FORM_BLOCK
  
  TexT start_unknown_txt;
  Int4 start_unknown;
} ConvertGapToUnknownData, PNTR ConvertGapToUnknownPtr;


static void FixDeltaFeatures (BioseqPtr bsp, Int4 offset, Int4 len_diff)
{
  SeqFeatPtr             sfp;
  SeqMgrFeatContext      fcontext;
  
  if (len_diff == 0)
  {
    return;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext))
  {
    AdjustFeatureForGapChange (sfp, bsp, offset, len_diff);
  }  
}


static void ConvertGapFeaturesToUnknownCallback (BioseqPtr bsp, Pointer userdata)
{
  ConvertGapToUnknownPtr cgtup;
  DeltaSeqPtr            dsp;
  SeqLitPtr              slip;
  Int4                   len_diff;
  Int4                   offset = 0, add_len;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL
      || userdata == NULL)
  {
    return;
  }
  
  cgtup = (ConvertGapToUnknownPtr) userdata;
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    add_len = GetDeltaSeqLen (dsp);
    if (IsDeltaSeqGap (dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp) && add_len >= cgtup->start_unknown)
    {
      slip = (SeqLitPtr) (dsp->data.ptrvalue);
      len_diff = slip->length - 100;
      slip->length = 100;
      if (slip->fuzz != NULL)
      {
        slip->fuzz = IntFuzzFree (slip->fuzz);
      }
      slip->fuzz = IntFuzzNew();
      slip->fuzz->choice = 4;
        
      if (len_diff > 0) {
        FixDeltaFeatures (bsp, offset, len_diff);
      }
      add_len -= len_diff;
    }
    offset += add_len;
    dsp = dsp->next;
  }
  bsp->length = offset;
}

static void ConvertGapFeaturesToUnknownButton (ButtoN b)
{
  ConvertGapToUnknownPtr cgtup;
  SeqEntryPtr            sep;
  CharPtr                str;
  
  cgtup = (ConvertGapToUnknownPtr) GetObjectExtra (b);
  if (cgtup == NULL)
  {
    return;
  }
  
  str = SaveStringFromText (cgtup->start_unknown_txt);
  if (StringHasNoText (str))
  {
    str = MemFree (str);
    return;
  }
  
  cgtup->start_unknown = atoi (str);
  str = MemFree (str);
  if (cgtup->start_unknown <= 0)
  {
    return;
  }
  
  WatchCursor();
  Update();
  
  sep = GetTopSeqEntryForEntityID (cgtup->input_entityID);
  
  VisitBioseqsInSep (sep, cgtup, ConvertGapFeaturesToUnknownCallback);
  
  ObjMgrSetDirtyFlag (cgtup->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cgtup->input_entityID, 0, 0);
  Remove (cgtup->form);
  ArrowCursor ();
  Update (); 
}

extern void ConvertGapFeaturesToUnknown (IteM i)
{
  BaseFormPtr            bfp;
  ConvertGapToUnknownPtr cgtup;
  WindoW                 w;
  GrouP                  h, g, c;
  ButtoN                 b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cgtup = (ConvertGapToUnknownPtr) MemNew (sizeof (ConvertGapToUnknownData));
  if (cgtup == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Convert Known Length Gaps to Unknown", StdCloseWindowProc);
  SetObjectExtra (w, cgtup, StdCleanupFormProc);
  cgtup->form = (ForM) w;
  cgtup->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Convert gaps longer or equal to", 0, popupMenuHeight, programFont, 'r');
  cgtup->start_unknown_txt = DialogText (g, "100", 10, NULL);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ConvertGapFeaturesToUnknownButton);
  SetObjectExtra (b, cgtup, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  
  Show (w);
}


typedef struct changegaplen
{
  FEATURE_FORM_BLOCK
  
  TexT   length_txt;
  Int4   length;
} ChangeGapLenData, PNTR ChangeGapLenPtr;


static void ChangeOneGapLength (SeqFeatPtr sfp, Int4 new_length)
{
  BioseqPtr   bsp;
  DeltaSeqPtr dsp;
  SeqLitPtr              slip;
  SeqLocPtr              slp;
  Int4                   len_diff = 0;
  Int4                   offset = 0, gap_start;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_gap || new_length < 0
      || new_length == SeqLocLen (sfp->location))
  {
    return;
  }
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return;
  }
  
  gap_start = SeqLocStart (sfp->location);
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && offset < gap_start)
  {
    if (dsp->choice == 1 && dsp->data.ptrvalue != NULL)
    {
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      offset += SeqLocLen (slp);
    }
    else if (dsp->choice == 2 && dsp->data.ptrvalue != NULL)
    {
      slip = (SeqLitPtr) (dsp->data.ptrvalue);
      offset += slip->length;
    }
    dsp = dsp->next;
  }
  
  if (offset == gap_start && dsp != NULL && dsp->choice == 2 && dsp->data.ptrvalue != NULL)
  {
    slip = (SeqLitPtr) (dsp->data.ptrvalue);
    if (IsDeltaSeqKnownGap (dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp))
    {
      len_diff = slip->length - new_length;
      slip->length = new_length;
      FixDeltaFeatures (bsp, offset, len_diff);
    }
  }
  
  
  bsp->length -= len_diff;
}


static void ChangeGapLength (ButtoN b)
{
  ChangeGapLenPtr   cglp;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  CharPtr           str;

  cglp = (ChangeGapLenPtr) GetObjectExtra (b);
  if (cglp == NULL)
  {
    return;
  }
  
  cglp->length = 0;
  str = SaveStringFromText (cglp->length_txt);
  if (!StringHasNoText (str))
  {
    cglp->length = atoi (str);
  }
  str = MemFree (str);
  if (cglp->length < 1)
  {
    Message (MSG_ERROR, "Must select a gap size greater than zero!");
    return;
  }
  
  sel = ObjMgrGetSelected ();
  if (sel == NULL)
  {
    Message (MSG_ERROR, "No gaps selected!");
    return;
  }
  WatchCursor ();
  Update ();
  while (sel != NULL)
  {
    if (sel->entityID == cglp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
    {
      sfp = SeqMgrGetDesiredFeature (cglp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_gap)
      {
        ChangeOneGapLength (sfp, cglp->length);
      }
    }
    sel = sel->next;
  }
  
  ObjMgrSetDirtyFlag (cglp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cglp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
  
  if (!GetStatus (cglp->leave_dlg_up))
  {
    Remove (cglp->form);
  }
}

extern void ChangeKnownGapLength (IteM i)
{
  BaseFormPtr     bfp;
  ChangeGapLenPtr cglp;
  WindoW          w;
  GrouP           h, g, c;
  ButtoN          b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cglp = (ChangeGapLenPtr) MemNew (sizeof (ChangeGapLenData));
  if (cglp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Change Length of Selected Known Length Gaps", StdCloseWindowProc);
  SetObjectExtra (w, cglp, StdCleanupFormProc);
  cglp->form = (ForM) w;
  cglp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Change length to", 0, popupMenuHeight, programFont, 'r');
  cglp->length_txt = DialogText (g, "100", 10, NULL);
  
  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ChangeGapLength);
  SetObjectExtra (b, cglp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cglp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  
  Show (w);
}

extern void ConvertSelectedGapFeaturesToKnown (IteM i)
{
  ConvertSelectedGapFeatures (i, TRUE);  
}

extern void ConvertSelectedGapFeaturesToUnknown (IteM i)
{
  ConvertSelectedGapFeatures (i, FALSE);
}


static Int4 CountNsAtEndOfSeqLit (SeqLitPtr slip, Uint2 which_end)
{
  Int2         residue;
  Int4         change_len = 0;
  
  if (slip == NULL || slip->seq_data == NULL || slip->seq_data_type == Seq_code_gap)
  {
    return 0;
  }

  if (slip->seq_data_type != Seq_code_iupacna)
  {
    slip->seq_data = (SeqDataPtr) BSConvertSeq((ByteStorePtr) slip->seq_data, Seq_code_iupacna, 
                          slip->seq_data_type, 
                          slip->length);
    slip->seq_data_type = Seq_code_iupacna;
  }
  
  if (which_end == SEQLOC_LEFT_END)
  {
    BSSeek((ByteStorePtr) slip->seq_data, 0, SEEK_SET);
  }
  else
  {
    BSSeek((ByteStorePtr) slip->seq_data, slip->length - 1, SEEK_SET);    
  }
  
  residue = BSGetByte((ByteStorePtr) slip->seq_data);
  while (change_len < slip->length && residue == 'N')
  {
    change_len++;
    if (which_end != SEQLOC_LEFT_END)
    {
      BSSeek ((ByteStorePtr) slip->seq_data, slip->length - change_len - 1, SEEK_SET);
    }
    residue = BSGetByte ((ByteStorePtr) slip->seq_data);
  }

  return change_len;
}

static void RemoveSeqLitEnd (SeqLitPtr slip, Int4 change_len, Uint2 which_end)
{
  ByteStorePtr bs_new;

  if (slip == NULL || change_len < 1)
  {
    return;
  }

  if (slip->seq_data_type == Seq_code_gap) return;

  if (slip->seq_data_type != Seq_code_iupacna)
  {
    slip->seq_data = (SeqDataPtr) BSConvertSeq((ByteStorePtr) slip->seq_data, Seq_code_iupacna, 
                                  slip->seq_data_type, 
                                  slip->length);
    slip->seq_data_type = Seq_code_iupacna;
  }
  bs_new = BSNew (slip->length - change_len);
  if (which_end == SEQLOC_LEFT_END)
  {
    BSSeek ((ByteStorePtr) slip->seq_data, change_len, SEEK_SET);
  }
  else
  {
    BSSeek((ByteStorePtr) slip->seq_data, 0, SEEK_SET);
  }
  BSInsertFromBS (bs_new, (ByteStorePtr) slip->seq_data, slip->length - change_len);
  slip->seq_data = SeqDataFree (slip->seq_data, slip->seq_data_type);
  slip->seq_data = (SeqDataPtr) bs_new;
  slip->length -= change_len;
}

static void ExpandGapsToIncludeFlankingNs (BioseqPtr bsp, Pointer userdata)
{
  DeltaSeqPtr  dsp, prev_dsp = NULL, next_dsp = NULL, prev_prev_dsp = NULL;
  Int4         change_len;
  SeqLitPtr    slip, prev_slip = NULL, next_slip;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return;
  }
  
  dsp = bsp->seq_ext;
  while (dsp != NULL) 
  {
    next_dsp = dsp->next;
    /* look for gap of known length */
    if (IsDeltaSeqKnownGap(dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp)) 
    {
      /* check for Ns before gap of known length */
      if (prev_dsp != NULL && prev_dsp->choice == 2 
          && prev_dsp->data.ptrvalue != NULL
          && !IsDeltaSeqGap (prev_dsp))
      {
        slip = (SeqLitPtr) dsp->data.ptrvalue;
        prev_slip = (SeqLitPtr) prev_dsp->data.ptrvalue;
            
        change_len = CountNsAtEndOfSeqLit (prev_slip, SEQLOC_RIGHT_END);
        if (change_len > 0)
        {
          RemoveSeqLitEnd (prev_slip, change_len, SEQLOC_RIGHT_END);
          slip->length += change_len;
          if (prev_slip->length == 0) {
            if (prev_prev_dsp == NULL) {
              bsp->seq_ext = dsp;
            } else {
              prev_prev_dsp->next = dsp;
            }
            prev_dsp->next = NULL;
            prev_dsp = DeltaSeqFree (prev_dsp);
            prev_dsp = dsp;
          }
        }
      } else {
        prev_prev_dsp = prev_dsp;
        prev_dsp = dsp;
      }
      /* check for Ns after gap of known length */
      if (dsp->next != NULL && dsp->next->choice == 2
          && dsp->next->data.ptrvalue != NULL
          && !IsDeltaSeqGap (dsp->next))
      {
        slip = (SeqLitPtr) dsp->data.ptrvalue;
        next_slip = (SeqLitPtr) dsp->next->data.ptrvalue;
        change_len = CountNsAtEndOfSeqLit (next_slip, SEQLOC_LEFT_END);
            
        if (change_len < next_slip->length)
        {
          RemoveSeqLitEnd (next_slip, change_len, SEQLOC_LEFT_END);
          slip->length += change_len;
        }
        else 
        {
          dsp->next = next_dsp->next;
          next_dsp->next = NULL;
          next_dsp = DeltaSeqFree (next_dsp);
          next_dsp = dsp->next;
        }
      }
    }
    dsp = next_dsp;
  }    
  BioseqPack (bsp);
}

extern void AddFlankingNsToKnownLengthGaps (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitBioseqsInSep (sep, NULL, ExpandGapsToIncludeFlankingNs);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static Boolean CanCombineDeltaSeq (DeltaSeqPtr dsp1, DeltaSeqPtr dsp2)
{

  if (dsp1 == NULL || dsp2 == NULL || dsp1->choice != 2 || dsp2->choice != 2
      || dsp1->data.ptrvalue == NULL || dsp2->data.ptrvalue == NULL
      || DoesDeltaSeqHaveGapTypeOrLinkage(dsp1) || DoesDeltaSeqHaveGapTypeOrLinkage(dsp2)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


extern void CombineAdjacentGapsOnBioseq (BioseqPtr bsp, Pointer userdata)
{
  SeqLitPtr  litp, pitp;
  Boolean    l_unknown, p_unknown;
  Int4       offset = 0;
  Int4       len_diff;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  DeltaSeqPtr prev, last;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta)
  {
    return;
  }

  /* combine adjacent gaps */
  prev = (DeltaSeqPtr) bsp->seq_ext;
  if (prev == NULL) return;
  last = prev->next;
  while (last != NULL) {
    if (IsDeltaSeqGap (prev) && IsDeltaSeqGap (last) && CanCombineDeltaSeq (prev, last)) {
      pitp = (SeqLitPtr) prev->data.ptrvalue;
      litp = (SeqLitPtr) last->data.ptrvalue;
      p_unknown = IsDeltaSeqUnknownGap (prev);
      l_unknown = IsDeltaSeqUnknownGap (last);
      len_diff = 0;
      if (p_unknown && l_unknown) {
        /* combine two unknown length gaps to make one unknown length gap */
        len_diff = litp->length;
      } else if (!p_unknown && !l_unknown) {
        /* combine two known length gaps into one known length gap */
        pitp->length += litp->length;
      } else {
        if (l_unknown) {
          /* swap gaps to put unknown length gap into prev */
          prev->data.ptrvalue = litp;
          last->data.ptrvalue = pitp;
          pitp = (SeqLitPtr) prev->data.ptrvalue;
          litp = (SeqLitPtr) last->data.ptrvalue;
        }
        /* remove length of known gap */
        len_diff = litp->length;
      } 
      prev->next = last->next;
      last->next = NULL;
      SeqLitFree (litp);
      MemFree (last);
      last = prev;
      
      if (len_diff > 0) {
        FixDeltaFeatures (bsp, offset, len_diff);
        bsp->length -= len_diff;
      }
    } else {
      offset += GetDeltaSeqLen (prev);
    }
    prev = last;
    last = last->next;
  }

  /* adjust coding region locations for unknown gaps */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext))
  {
    AdjustCDSLocationsForUnknownGapsCallback (sfp, NULL);
  }
}

extern void CombineAdjacentGaps (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitBioseqsInSep (sep, NULL, CombineAdjacentGapsOnBioseq);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}

static void MarkPseudoGenesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqFeatPtr gene;
  
  if (sfp == NULL || ! sfp->pseudo 
      || sfp->data.choice == SEQFEAT_GENE
      || SeqMgrGetGeneXref (sfp) != NULL)
  {
    return;
  }
  
  gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene != NULL) 
  {
    gene->pseudo = TRUE;
  }
}

extern void MarkGenesWithPseudoFeaturesPseudo (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitFeaturesInSep (sep, NULL, MarkPseudoGenesCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
  
}

static void RemoveOneNomenclature (UserObjectPtr PNTR puop)
{
  UserObjectPtr  uop, obj, prev_obj, next_obj;
  UserFieldPtr   prev_ufp, next_ufp, ufp;
  ObjectIdPtr    oip;

  if (puop == NULL || *puop == NULL) return;
  uop = *puop;
  
  for (ufp = uop->data, prev_ufp = NULL; 
       ufp != NULL; 
       ufp = next_ufp) {
    next_ufp = ufp->next;
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      RemoveOneNomenclature (&obj);
      ufp->data.ptrvalue = obj;
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue, prev_obj = NULL;
           obj != NULL;
           obj = next_obj) {
        next_obj = obj->next;
        RemoveOneNomenclature (&obj);
        if (obj == NULL)
        {
          if (prev_obj == NULL)
          {
            ufp->data.ptrvalue = next_obj;
          }
          else
          {
            prev_obj->next = next_obj;
          }
          obj = UserObjectFree (obj);
        }
        else
        {
          prev_obj = obj;
        }
      }
    }
    if ((ufp->choice == 6 || ufp->choice == 12) && ufp->data.ptrvalue == NULL)
    {
      if (prev_ufp == NULL)
      {
        uop->data = ufp->next;
      }
      else
      {
        prev_ufp->next = ufp->next;
      }
      ufp = UserFieldFree (ufp);
    }
    else
    {
      prev_ufp = ufp;
    }
  }
  
  oip = uop->type;
  if (oip != NULL && StringCmp (oip->str, "OfficialNomenclature") == 0)
  {
    uop = UserObjectFree (uop);
    *puop = uop;
  }
}

static void RemoveNomenclatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  UserObjectPtr uop;
  
  if (sfp != NULL)
  {
    uop = sfp->ext;
    RemoveOneNomenclature (&uop);
    sfp->ext = uop;
  }
}

extern void RemoveNomenclature (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  VisitFeaturesInSep (sep, NULL, RemoveNomenclatureCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static void RemoveUnindexedFeaturesInSeqEntry (SeqEntryPtr sep, Uint2 entityID)
{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  SeqAnnotPtr       sap = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;

  if (sep == NULL || sep->data.ptrvalue == NULL)
  {
	return;
  }

  if (IS_Bioseq(sep))
  {
	bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;

	sep = bssp->seq_set;
	while (sep != NULL)
	{
      RemoveUnindexedFeaturesInSeqEntry (sep, entityID);
	  sep = sep->next;
	}
    sap = bssp->annot;
  }

  while (sap != NULL)
  {
    if (sap->type == 1)
	{
      sfp = (SeqFeatPtr) sap->data;
	  while (sfp != NULL)
	  {
	    if (SeqMgrGetDesiredFeature (entityID, NULL, 0, 0, sfp, &context) == NULL)
		{
		  sfp->idx.deleteme = TRUE;
		}
		sfp = sfp->next;
	  }
	}
	sap = sap->next;
  }
}


extern void RemoveUnindexedFeatures (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  RemoveUnindexedFeaturesInSeqEntry (sep, bfp->input_entityID);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static void CopyLocusToLocusTagCallback (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr grp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (!StringHasNoText (grp->locus) && StringHasNoText (grp->locus_tag))
  {
    grp->locus_tag = MemFree (grp->locus_tag);
    grp->locus_tag = StringSave (grp->locus);
  }
}

extern void CopyLocusToLocusTag (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  VisitFeaturesInSep (sep, NULL, CopyLocusToLocusTagCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


/* data structure and functions for a generic form displaying a clickable list */

#define CLICKABLE_LIST_FORM_BLOCK   \
  FORM_MESSAGE_BLOCK                \
  ValNodePtr      clickable_list_data; \
  DialoG          clickable_list_dlg;  \
  ButtoN          recheck_btn;         \
  
typedef struct clickablelistform {
  CLICKABLE_LIST_FORM_BLOCK
} ClickableListFormData, PNTR ClickableListFormPtr;

static void CleanupClickableListForm (GraphiC g, VoidPtr data)

{
  ClickableListFormPtr clfp;

  clfp = (ClickableListFormPtr) data;
  if (clfp != NULL) {
    clfp->clickable_list_data = FreeClickableList (clfp->clickable_list_data);
    ObjMgrFreeUserData (clfp->input_entityID, clfp->procid, clfp->proctype, clfp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static Int2 LIBCALLBACK ClickableListFormMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  ClickableListFormPtr     clfp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  clfp = (ClickableListFormPtr) omudp->userdata.ptrvalue;
  if (clfp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (clfp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          clfp->clickable_list_data = FreeClickableList (clfp->clickable_list_data);
          PointerToDialog (clfp->clickable_list_dlg, NULL);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          clfp->clickable_list_data = FreeClickableList (clfp->clickable_list_data);
          PointerToDialog (clfp->clickable_list_dlg, NULL);
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}


static void ClickableListFormMessage (ForM f, Int2 mssg)

{
  ClickableListFormPtr drfp;

  drfp = (ClickableListFormPtr) GetObjectExtra (f);
  if (drfp != NULL) {
    switch (mssg) {
      case VIB_MSG_EXPORT :
        if (drfp->exportform != NULL) {
          (drfp->exportform) (f, NULL);
        }
        break;
      case VIB_MSG_PRINT :
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
      case VIB_MSG_COPY :
        SendMessageToDialog (drfp->clickable_list_dlg, VIB_MSG_COPY);
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        drfp->clickable_list_data = FreeClickableList (drfp->clickable_list_data);
        PointerToDialog (drfp->clickable_list_dlg, NULL);
        break;
      default :
        if (drfp->appmessage != NULL) {
          drfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

/* Discrepancy Report */

/* There will only be one Discrepancy Report window at a time */
static WindoW discrepancyReportWindow = NULL;
static WindoW oncallerReportWindow = NULL;

static WindoW GetWindowForReportType (EDiscrepancyReportType report_type)
{
  WindoW w = NULL;

  switch (report_type) {
    case eReportTypeDiscrepancy:
      w = discrepancyReportWindow;
      break;
    case eReportTypeOnCaller:
      w = oncallerReportWindow;
      break;
  }
  return w;
}


static void ClearWindowForReportType (WindoW w)
{
  if (discrepancyReportWindow == w) {
    discrepancyReportWindow = NULL;
  }
  if (oncallerReportWindow == w) {
    oncallerReportWindow = NULL;
  }
}


static void SetWindowForReportType (WindoW w, EDiscrepancyReportType report_type)
{
  switch (report_type) {
    case eReportTypeDiscrepancy:
      discrepancyReportWindow = w;
      break;
    case eReportTypeOnCaller:
      oncallerReportWindow = w;
      break;
  }
}


typedef void (*DiscrepancyCallback) (ValNodePtr item_list, Pointer userdata);
typedef void (*DiscrepancyCallbackDataFree) (Pointer userdata);




typedef struct discrepancyreportform 
{
  CLICKABLE_LIST_FORM_BLOCK

  DiscrepancyConfigPtr dcp;
} DiscrepancyReportFormData, PNTR DiscrepancyReportFormPtr;

static void CleanupDiscrepancyReportForm (GraphiC g, VoidPtr data)

{
  DiscrepancyReportFormPtr drfp;
  ValNodePtr               vnp;
  ClickableItemPtr         cip;

  drfp = (DiscrepancyReportFormPtr) data;
  if (drfp != NULL) {
    /* find whether source qual report is open or closed */
    for (vnp = drfp->clickable_list_data; vnp != NULL; vnp = vnp->next) {
      cip = vnp->data.ptrvalue;
      if (cip != NULL && cip->clickable_item_type == DISC_SRC_QUAL_PROBLEM) {
        if (cip->expanded) {
          SetAppParam ("SEQUINCUSTOM", "ONCALLERTOOL", "EXPAND_SRCQUAL_REPORT", "TRUE");
        } else {
          SetAppParam ("SEQUINCUSTOM", "ONCALLERTOOL", "EXPAND_SRCQUAL_REPORT", "FALSE");
        }
        break;
      }
    }
    drfp->clickable_list_data = FreeClickableList (drfp->clickable_list_data);
    drfp->dcp = DiscrepancyConfigFree (drfp->dcp);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
    ClearWindowForReportType ((WindoW) drfp->form);
  }
  StdCleanupFormProc (g, data);
}


static void UnselectDiscrepancyList(ButtoN b)
{
  ButtoN *test_options;
  Int4    i;
  
  test_options = (ButtoN *) GetObjectExtra (b);
  if (test_options != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      SetStatus (test_options[i], FALSE);
    }
  }
}

static void SelectDiscrepancyList(ButtoN b)
{
  ButtoN *test_options;
  Int4    i;
  
  test_options = (ButtoN *) GetObjectExtra (b);
  if (test_options != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      SetStatus (test_options[i], TRUE);
    }
  }
}


/* This function returns TRUE if there was a change to the discrepancy config,
 * FALSE otherwise.
 */
static Boolean EditDiscrepancyConfig (DiscrepancyConfigPtr dcp, EDiscrepancyReportType report_type)
{
  WindoW                w;
  GrouP                 h, g, k, c;
  ButtoN                b, use_feature_table_format_btn;
  ModalAcceptCancelData acd;
  Int4                  i;
  ButtoN                test_options[MAX_DISC_TYPE];
  Boolean               rval = FALSE;
  
  if (dcp == NULL)
  {
    return rval;
  }
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = NormalGroup (h, 0, 14, "Discrepancy Tests to Run", programFont, NULL);
  SetGroupSpacing (g, 10, 10);
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if (IsTestTypeAppropriateForReportType (i, report_type)) {
      test_options[i] = CheckBox (g, GetDiscrepancyTestConfName ((DiscrepancyType) i), NULL);
      SetStatus (test_options[i], dcp->conf_list[i]);
    } else {
      test_options[i] = NULL;
    }
  }
  
  use_feature_table_format_btn = CheckBox (h, "Use feature table format for features in report", NULL);
  SetStatus (use_feature_table_format_btn, dcp->use_feature_table_format);
  
  k = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (k, "Select All", SelectDiscrepancyList);
  SetObjectExtra (b, test_options, NULL);
  b = PushButton (k, "Unselect All", UnselectDiscrepancyList);
  SetObjectExtra (b, test_options, NULL);  
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) use_feature_table_format_btn, (HANDLE) k, (HANDLE) c, NULL);
  
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
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      dcp->conf_list [i] = GetStatus (test_options[i]);
    }
    dcp->use_feature_table_format = GetStatus (use_feature_table_format_btn);
    rval = TRUE;
    SaveDiscrepancyConfig (dcp);
  }

  Remove (w);
  return rval;
}

static void VisitFeaturesInViewedBioseqs (Pointer userdata, VisitFeaturesFunc callback)
{
  ValNodePtr  base_form_list, vnp;
  BaseFormPtr bfp;
  SeqEntryPtr sep;
  SeqEntryPtr orig_scope;
  
  orig_scope = SeqEntryGetScope ();
  base_form_list = GetBaseFormList();
  for (vnp = base_form_list; vnp != NULL; vnp = vnp->next) {
    bfp = (BaseFormPtr) vnp->data.ptrvalue;
    if (bfp != NULL && bfp->input_entityID != 0) {
      sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
      SeqEntrySetScope (sep);
      VisitFeaturesInSep (sep, userdata, callback);
    }
  }
  base_form_list = ValNodeFree (base_form_list);
  SeqEntrySetScope (orig_scope);
}


extern Int4 CountChosenDiscrepancies (ValNodePtr discrepancy_list, Boolean count_all)
{
  Int4               num_chosen = 0;
  ClickableItemPtr dip;
  
  while (discrepancy_list != NULL)
  {
    dip = (ClickableItemPtr) discrepancy_list->data.ptrvalue;
    if (dip != NULL)
    {
      if (dip->chosen || count_all)
      {
        if (dip->expanded && dip->subcategories != NULL)
        {
          num_chosen += CountChosenDiscrepancies (dip->subcategories, TRUE);
        }
        else
        {
          num_chosen ++;
        }
      }
      else if (dip->expanded)
      {
        num_chosen += CountChosenDiscrepancies (dip->subcategories, FALSE);
      }
    }
    discrepancy_list = discrepancy_list->next;
  }
  return num_chosen;
}

extern void EditDiscrepancyItem (ValNodePtr vnp, Pointer userdata)
{
  SeqFeatPtr sfp, cds;
  BioseqPtr  bsp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  
  if (vnp == NULL)
  {
    return;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (sfp->idx.subtype == FEATDEF_PROT) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            sfp = cds;
          }
        }
      }
      GatherProcLaunch (OMPROC_EDIT, FALSE, sfp->idx.entityID, sfp->idx.itemID,
                        OBJ_SEQFEAT, 0, 0, OBJ_SEQFEAT, 0);
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      GatherProcLaunch (OMPROC_EDIT, FALSE, bsp->idx.entityID, bsp->idx.itemID,
                         OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
    }
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) (vnp->data.ptrvalue);
    if (sdp != NULL && sdp->extended != 0)
    {
      ovp = (ObjValNodePtr) sdp;
      GatherProcLaunch (OMPROC_EDIT, FALSE, ovp->idx.entityID, ovp->idx.itemID,
                         OBJ_SEQDESC, 0, 0, OBJ_SEQDESC, 0);
    }
  }

}


extern void BulkEditDiscrepancy (ValNodePtr vnp, Pointer userdata)
{
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  ValNodePtr       feat_list, desc_list;
  ValNodePtr       feat_list_list = NULL, desc_list_list = NULL, list_vnp;
  Uint2            entityID = 0;
  
  if (vnp == NULL)
  {
    return;
  }
  while (vnp != NULL) 
  {
    if (vnp->choice == OBJ_SEQFEAT) 
    {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp != NULL)
      {
        for (list_vnp = feat_list_list;
             list_vnp != NULL && list_vnp->choice != sfp->idx.entityID;
             list_vnp = list_vnp->next) 
        {
        }
        
        if (list_vnp == NULL) {
          list_vnp = ValNodeAddPointer (&feat_list_list, sfp->idx.entityID, NULL);
        }
        feat_list = (ValNodePtr) list_vnp->data.ptrvalue;
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
        list_vnp->data.ptrvalue = feat_list;
        entityID = sfp->idx.entityID;
      }
    }
    else if (vnp->choice == OBJ_SEQDESC) 
    {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp != NULL) 
      {
        for (list_vnp = desc_list_list;
             list_vnp != NULL && list_vnp->choice != sdp->choice;
             list_vnp = list_vnp->next) 
        {
        }
        
        if (list_vnp == NULL) {
          list_vnp = ValNodeAddPointer (&desc_list_list, sdp->choice, NULL);
        }
        desc_list = (ValNodePtr) list_vnp->data.ptrvalue;
        ValNodeAddPointer (&desc_list, OBJ_SEQDESC, sdp);
        list_vnp->data.ptrvalue = desc_list;
        if (sdp->extended > 0) {
          entityID = ((ObjValNodePtr)sdp)->idx.entityID;
        }
      }
    }
    vnp = vnp->next;
  }

  for (list_vnp = feat_list_list; list_vnp != NULL; list_vnp = list_vnp->next) {
    BulkEditorFeatList (entityID, list_vnp->data.ptrvalue);
  }
  feat_list_list = ValNodeFree (feat_list_list);

  for (list_vnp = desc_list_list; list_vnp != NULL; list_vnp = list_vnp->next) {
    BulkEditorDescrList (entityID, list_vnp->data.ptrvalue);
  }
  desc_list_list = ValNodeFree (desc_list_list);

}


static void EditCDStRNAOverlap (ValNodePtr item_list);

static void EditCDStRNAOverlapCallback (ValNodePtr vnp, Pointer userdata)
{
  EditCDStRNAOverlap (vnp); 
}


static void ExtendPartialsToEndOrGapCallback (ValNodePtr list, Pointer userdata)
{
  ValNodePtr vnp, entityID_list = NULL;
  SeqFeatPtr sfp;
  LogInfoPtr lip;
  CharPtr    orig_location, new_location;

  if (Message (MSG_OKC, "Extend partial ends of features to gaps/end of sequence if within two nucleotides?") == ANS_CANCEL) {
    return;
  }
  WatchCursor();
  Update();
  lip = OpenLog ("Extended Features");
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && (sfp = vnp->data.ptrvalue) != NULL) {
      orig_location = SeqLocPrintUseBestID (sfp->location);
      if (ExtendPartialsToEndOrGap (sfp)) {
        ValNodeAddInt (&entityID_list, 0, sfp->idx.entityID);
        new_location = SeqLocPrintUseBestID (sfp->location);
        fprintf (lip->fp, "Extended %s to %s\n", orig_location, new_location);
        new_location = MemFree (new_location);
        lip->data_in_log = TRUE;
      }
      orig_location = MemFree (orig_location);
    }
  }
  entityID_list = ValNodeSort (entityID_list, SortByIntvalue);
  ValNodeUnique (&entityID_list, SortByIntvalue, ValNodeFree);
  for (vnp = entityID_list; vnp != NULL; vnp = vnp->next) {
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);
  }
  ArrowCursor();
  Update();
  CloseLog (lip);
  lip = FreeLog (lip);
}


extern void SetBioseqViewTargetByBioseq (BaseFormPtr bfp, BioseqPtr bsp)
{
  Char       id_text [41];
  
  if (bsp != NULL && bfp != NULL)
  {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_text, PRINTID_REPORT, sizeof (id_text));
    SetBioseqViewTarget (bfp, id_text);
  }
}


static BioseqPtr GetFirstBioseqInSeqEntry (SeqEntryPtr sep)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else if (IS_Bioseq (sep))
  {
    bsp = sep->data.ptrvalue;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL && bsp == NULL; sep = sep->next)
    {
      bsp = GetFirstBioseqInSeqEntry (sep);
    }
  }
  return bsp;
}


static BioseqPtr GetBioseqForDescriptor (ObjValNodePtr ovp)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;
  
  if (ovp == NULL || ovp->idx.parentptr == NULL)
  {
    return NULL;
  }
  else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
    bsp = (BioseqPtr) ovp->idx.parentptr;
  } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) ovp->idx.parentptr;
    for (sep = bssp->seq_set; sep != NULL && bsp == NULL; sep = sep->next)
    {
      bsp = GetFirstBioseqInSeqEntry (sep);
    }  
  }
  return bsp;
}


extern BaseFormPtr GetBaseFormForEntityID (Uint2 entityID)
{
  ValNodePtr base_form_list, vnp;
  BaseFormPtr bfp = NULL;
  
  base_form_list = GetBaseFormList();
  vnp = base_form_list;
  while (vnp != NULL && bfp == NULL) {
    bfp = (BaseFormPtr) vnp->data.ptrvalue;
    if (bfp != NULL && bfp->input_entityID != entityID) {
      bfp = NULL;
    }
    vnp = vnp->next;
  }
  base_form_list = ValNodeFree (base_form_list);
  return bfp;
}

extern void ScrollToDiscrepancyItem (ValNodePtr vnp, Pointer userdata)
{
  SeqFeatPtr    sfp, cds;
  BioseqPtr     bsp;
  SeqDescrPtr   sdp;
  ObjValNodePtr ovp;
  BaseFormPtr   bfp;
    
  if (vnp == NULL)
  {
    return;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (sfp->idx.subtype == FEATDEF_PROT) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            sfp = cds;
          }
        }
      }
    
      /* need to scroll to item */
      bsp = BioseqFindFromSeqLoc (sfp->location);      
      bfp = GetBaseFormForEntityID (sfp->idx.entityID);
      if (bfp != NULL) {
        Select (bfp->form);
        SetBioseqViewTargetByBioseq (bfp, bsp);
        ObjMgrSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
      }
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    bfp = GetBaseFormForEntityID (bsp->idx.entityID);
    if (bfp != NULL) {
      Select (bfp->form);
      SetBioseqViewTargetByBioseq (bfp, bsp);
    }
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) (vnp->data.ptrvalue);
    if (sdp != NULL && sdp->extended != 0)
    {
      ovp = (ObjValNodePtr) sdp;      
      bsp = GetBioseqForDescriptor (ovp);
      bfp = GetBaseFormForEntityID (bsp->idx.entityID);
      if (bfp != NULL) {
        Select (bfp->form);
        SetBioseqViewTargetByBioseq (bfp, bsp);
        ObjMgrSelect (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, 0, NULL);
      }
    }
  }
}

static Uint2 GetEntityIDFromItem (ValNodePtr vnp);
static void ApplyTagToCodingRegionsInEntityID (Uint2 entityID);

static void ApplyTagToCodingRegionsCallback (ValNodePtr item_list, Pointer userdata)
{
  ValNodePtr vnp;
  ValNodePtr entityIDList = NULL;
  Uint2       entityID;
  SeqEntryPtr sep;
  ValNodePtr  cds_list;
  Boolean     found_any = FALSE;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next)
  {
    entityID = GetEntityIDFromItem(vnp);
    if (!EntityIDAlreadyInList(entityID, entityIDList))
    {
      ValNodeAddInt (&entityIDList, 0, entityID);
    }
  }
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next)
  {
    sep = GetTopSeqEntryForEntityID (vnp->data.intvalue);
    cds_list = ListCodingRegionsContainedInSourceFeatures (sep);
    if (cds_list != NULL)
    {
      found_any = TRUE;
      cds_list = ValNodeFree (cds_list);
      ApplyTagToCodingRegionsInEntityID (vnp->data.intvalue);
    }
  }
  if (!found_any)
  {
    Message (MSG_ERROR, "No coding regions found in source features!  Try editing the definition lines.");
  }  
  entityIDList = ValNodeFree (entityIDList);  
}

static void AddBulkEditing (ValNodePtr clickable_list)
{
  ClickableItemPtr cip;
  Uint1            subtype;

  while (clickable_list != NULL) {
    cip = (ClickableItemPtr) clickable_list->data.ptrvalue;
    if (cip->callback_func == NULL && cip->item_list != NULL) {
      if (cip->clickable_item_type == DISC_CDS_OVERLAP_TRNA) {
        cip->callback_func = EditCDStRNAOverlapCallback;    
      } else if (cip->clickable_item_type == DISC_INCONSISTENT_BIOSRC_DEFLINE) {
        cip->callback_func = ApplyTagToCodingRegionsCallback; 
      } else if (cip->clickable_item_type == DISC_BACTERIAL_PARTIAL_PROBLEMS) {
        cip->callback_func = ExtendPartialsToEndOrGapCallback;
      } else {
        subtype = GetSubtypeForBulkEdit (cip->item_list);
        /* Note - using FEATDEF_rRNA to represent all editable RNA features */
        if (subtype == FEATDEF_CDS || subtype == FEATDEF_GENE || subtype == FEATDEF_rRNA) {
          cip->callback_func = BulkEditDiscrepancy;
        }
      }
    }
    AddBulkEditing (cip->subcategories);
    clickable_list = clickable_list->next;
  }
}

static void RecheckDiscrepancyProc (ButtoN b)
{
  DiscrepancyReportFormPtr drfp;
  ValNodePtr               sep_list;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (drfp != NULL)
  {
    WatchCursor();
    Update();
    drfp->clickable_list_data = FreeClickableList (drfp->clickable_list_data);
    
    sep_list = GetViewedSeqEntryList ();
    drfp->clickable_list_data = CollectDiscrepancies (drfp->dcp, sep_list, CheckTaxNamesAgainstTaxDatabase);

    /* add bulk editing where appropriate */
    AddBulkEditing (drfp->clickable_list_data);
    
    PointerToDialog (drfp->clickable_list_dlg, drfp->clickable_list_data);
    ArrowCursor();
    Update();
  }
}


extern void 
WriteClickableListReport 
(FILE       *fp,
 ValNodePtr discrepancy_list, 
 Boolean    show_all,
 Boolean    use_feature_table_fmt)
{
  ClickableItemPtr       dip;
  ValNodePtr               vnp;
  Int4                     num_chosen;

  if (fp == NULL || discrepancy_list == NULL)
  {
    return;
  }
  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next)
  {
    dip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (dip != NULL)
    {
      if (dip->expanded)
      {
        num_chosen = CountChosenDiscrepancies (dip->subcategories, show_all | dip->chosen);
        if (num_chosen > 0)
        {
          if (dip->chosen || show_all)
          {
            fprintf (fp, "%s\n", dip->description);
          }
          WriteClickableListReport (fp, dip->subcategories, show_all | dip->chosen, use_feature_table_fmt);          
        }
      }
      else if (dip->chosen || show_all)
      {
        WriteDiscrepancy (fp, dip, use_feature_table_fmt);
      }
    }
  }
}

static Boolean DiscrepancyReportExportProc (ForM f, CharPtr filename)

{
  FILE           *fp;
  Char           path [PATH_MAX];
  DiscrepancyReportFormPtr drfp;
  Int4                     num_disc = 0;
  Boolean                  show_all = FALSE;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (f);
  if (drfp == NULL) 
  {
    return FALSE;
  }
  
  num_disc = CountChosenDiscrepancies (drfp->clickable_list_data, FALSE);

  if (num_disc == 0) 
  {
    if (ANS_CANCEL == Message (MSG_OKC, "No discrepancies selected!  Export all?"))
    {
      return FALSE;
    }
    else
    {
      show_all = TRUE;
    }
  }
  
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      FileClose (fp);
    } else {
      FileCreate (path, "TEXT", "ttxt");
    }
#endif
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      WriteClickableListReport (fp, drfp->clickable_list_data, show_all, 
                              (Boolean)(drfp->dcp != NULL && drfp->dcp->use_feature_table_format));
      FileClose (fp);
      return TRUE;
    }
  }
  return FALSE;
}


static void GenerateDiscrepancyReport (ButtoN b)
{
  DiscrepancyReportFormPtr drfp;
  Char                     path [PATH_MAX];

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (drfp == NULL)
  {
    return;
  }

  TmpNam (path);  
  if (DiscrepancyReportExportProc (drfp->form, path))
  {
    LaunchGeneralTextViewer (path, "Discrepancy Report");
  }
  FileRemove (path);  
}


static void ReactivateDiscrepancyReport (EDiscrepancyReportType report_type)
{
  DiscrepancyReportFormPtr drfp;
  WindoW                   w;

  w = GetWindowForReportType (report_type);
  if (w == NULL) 
  {
    CreateReportWindow (report_type);
  }
  
  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (w);
  if (drfp == NULL)
  {
    Remove (w);
    ClearWindowForReportType (w);
    CreateReportWindow (report_type);
    w = GetWindowForReportType (report_type);
  }
    
  /* populate discrepancy lists */
  RecheckDiscrepancyProc (drfp->recheck_btn);
  Show (w);  
  Select (w);
}



static void EditReportConfigBtn (EDiscrepancyReportType report_type)
{
  DiscrepancyReportFormPtr drfp;
  WindoW                   w;

  w = GetWindowForReportType (report_type);
  
  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (w);
  if (drfp == NULL) return;
  
  if (EditDiscrepancyConfig (drfp->dcp, report_type))
  {
    RecheckDiscrepancyProc (drfp->recheck_btn);
  }
}


static void EditDiscrepancyConfigBtn (ButtoN b)
{
  EditReportConfigBtn (eReportTypeDiscrepancy);
}


static void EditOnCallerConfigBtn (ButtoN b)
{
  EditReportConfigBtn (eReportTypeOnCaller);
}


static Nlm_BtnActnProc GetReportEditButtonProc (EDiscrepancyReportType report_type)
{
  Nlm_BtnActnProc proc = NULL;

  switch (report_type) {
    case eReportTypeDiscrepancy:
      proc = EditDiscrepancyConfigBtn;
      break;
    case eReportTypeOnCaller:
      proc = EditOnCallerConfigBtn;
      break;
  }
  return proc;
}


#ifndef WIN_MAC
extern void CreateStdValidatorFormMenus (WindoW w);
#endif


static CharPtr GetReportName (EDiscrepancyReportType report_type)
{
  CharPtr report_name = "";

  switch (report_type) {
    case eReportTypeDiscrepancy:
      report_name = "Discrepancy Report";
      break;
    case eReportTypeOnCaller:
      report_name = "On Caller Tool";
      break;
  }
  return report_name;
}


static CharPtr GetReportConfigName (EDiscrepancyReportType report_type)
{
  CharPtr report_name = "";

  switch (report_type) {
    case eReportTypeDiscrepancy:
      report_name = "DISCREPANCY_REPORT";
      break;
    case eReportTypeOnCaller:
      report_name = "ON_CALLER_TOOL";
      break;
  }
  return report_name;
}


static void AdjustConfigForReportType (EDiscrepancyReportType report_type, DiscrepancyConfigPtr dcp)
{
  Int4 i;

  if (dcp == NULL) {
    return;
  }

  for (i = 0; i < MAX_DISC_TYPE; i++) {
    if (!IsTestTypeAppropriateForReportType (i, report_type)) {
      dcp->conf_list[i] = FALSE;
    }
  }

}


static void ExpandAllDiscReportItems (ButtoN b)
{
  DiscrepancyReportFormPtr d;

  d = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (d == NULL) {
    return;
  }

  ExpandClickableItemList (d->clickable_list_data);

  PointerToDialog (d->clickable_list_dlg, d->clickable_list_data);
}


static void ContractAllDiscReportItems (ButtoN b)
{
  DiscrepancyReportFormPtr d;

  d = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (d == NULL) {
    return;
  }
  ContractClickableItemList (d->clickable_list_data);

  PointerToDialog (d->clickable_list_dlg, d->clickable_list_data);
}


extern void CreateReportWindow (EDiscrepancyReportType report_type)
{
  DiscrepancyReportFormPtr drfp;
  GrouP                    h;
  ButtoN                   b;
  GrouP                    c;
  GrouP                    c1;
  WindoW                   w;
  OMUserDataPtr            omudp;

  if (GetWindowForReportType(report_type) != NULL)
  {
    ReactivateDiscrepancyReport (report_type);
    return; 
  }
  
  drfp = (DiscrepancyReportFormPtr) MemNew (sizeof (DiscrepancyReportFormData));
  if (drfp == NULL)
  {
    return;
  }
  
  w = FixedWindow (-50, -33, -10, -10, GetReportName (report_type), StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupDiscrepancyReportForm);
  drfp->form = (ForM) w;
  drfp->formmessage = ClickableListFormMessage;
  drfp->exportform = DiscrepancyReportExportProc;
  
  /* read in config file */
  drfp->dcp = ReadDiscrepancyConfigEx(GetReportConfigName (report_type));

  /* adjust for report type */
  AdjustConfigForReportType (report_type, drfp->dcp);
  
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (0, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = ClickableListFormMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list_dlg = CreateClickableListDialog (h, "Discrepancies", "Affected Items",
                                                    ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                    GetDiscrepancyItemText);
                              
  c1 = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c1, 10, 10);
  b = PushButton (c1, "Expand All", ExpandAllDiscReportItems);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c1, "Contract All", ContractAllDiscReportItems);
  SetObjectExtra (b, drfp, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Generate Report", GenerateDiscrepancyReport);
  SetObjectExtra (b, drfp, NULL);
  drfp->recheck_btn = PushButton (c, "Recheck", RecheckDiscrepancyProc);
  SetObjectExtra (drfp->recheck_btn, drfp, NULL);
  
  b = PushButton (c, "Configure", GetReportEditButtonProc (report_type));
  SetObjectExtra (b, drfp, NULL);
  
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list_dlg, (HANDLE) c1, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  /* populate discrepancy lists */
  RecheckDiscrepancyProc (drfp->recheck_btn);
  Show (w);
  SetWindowForReportType (w, report_type);
}


typedef struct sucform 
{
  CLICKABLE_LIST_FORM_BLOCK

  Boolean reverse;
  Boolean byblock;
  Boolean showsequence;
} SUCFormData, PNTR SUCFormPtr;

static void ReSUCProc (ButtoN b)
{
  SUCFormPtr sfp;
  SeqEntryPtr sep;
  sfp = (SUCFormPtr) GetObjectExtra (b);
  if (sfp != NULL)
  {
    sep = GetTopSeqEntryForEntityID (sfp->input_entityID);
    WatchCursor();
    Update();    
    sfp->clickable_list_data = FreeClickableList (sfp->clickable_list_data);
    sfp->clickable_list_data = GetSUCCommonList (sep, sfp->reverse, sfp->byblock, sfp->showsequence, TRUE);    
    if (sfp->byblock) {
      sfp->clickable_list_data = CategorizeSUCBlocks (sfp->clickable_list_data);
    }
    AddBulkEditing (sfp->clickable_list_data);  
    PointerToDialog (sfp->clickable_list_dlg, sfp->clickable_list_data);
    ArrowCursor();
    Update();
  }
}

static void SetClickableItemExpanded (ValNodePtr vnp, Boolean expand)
{
  ClickableItemPtr cip;

  while (vnp != NULL) {
    cip = vnp->data.ptrvalue;
    if (cip != NULL) {
      cip->expanded = expand;
      SetClickableItemExpanded (cip->subcategories, expand);
    }
    vnp = vnp->next;
  }
}

static void SUCExpand (ButtoN b, Boolean expand)
{
  SUCFormPtr       sfp;

  sfp = (SUCFormPtr) GetObjectExtra (b);
  if (sfp != NULL)
  {
    SetClickableItemExpanded (sfp->clickable_list_data, expand);
    PointerToDialog (sfp->clickable_list_dlg, sfp->clickable_list_data);
  }
}

static void SUCExpandAll (ButtoN b)
{
  SUCExpand (b, TRUE);
}

static void SUCCollapseAll (ButtoN b)
{
  SUCExpand (b, FALSE);
}


extern void NewSUC (ValNodePtr suc_list, Uint2 entityID, Boolean reverse, Boolean byblock, Boolean showsequence)
{
  SUCFormPtr    sfp;
  GrouP         h;
  GrouP         c;
  WindoW        w;
  ButtoN        b;
  OMUserDataPtr omudp;

  sfp = (SUCFormPtr) MemNew (sizeof (SUCFormData));
  if (sfp == NULL)
  {
    return;
  }
  
  w = FixedWindow (-50, -33, -10, -10, "SUC", StdCloseWindowProc);
  SetObjectExtra (w, sfp, CleanupClickableListForm);
  sfp->form = (ForM) w;
  sfp->formmessage = ClickableListFormMessage;
  sfp->exportform = DiscrepancyReportExportProc;
    
  /* register to receive update messages */
  sfp->input_entityID = entityID;
  sfp->userkey = OMGetNextUserKey ();
  sfp->procid = 0;
  sfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (sfp->input_entityID, sfp->procid, sfp->proctype, sfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) sfp;
    omudp->messagefunc = ClickableListFormMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif

  /* set up parameters for re-SUC */
  sfp->reverse = reverse;
  sfp->byblock = byblock; 
  sfp->showsequence = showsequence;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  sfp->clickable_list_dlg = CreateClickableListDialogEx (h, "Text", "Affected Items", NULL, NULL,
                                                         ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                         GetDiscrepancyItemText,
                                                         stdCharWidth * 55,
                                                         stdCharWidth * 55,
                                                         FALSE, TRUE);
  if (byblock) {
    suc_list = CategorizeSUCBlocks (suc_list);
    SetClickableItemExpanded (suc_list, FALSE);
  }
  AddBulkEditing (suc_list);  

  sfp->clickable_list_data = suc_list;                                                  
  PointerToDialog (sfp->clickable_list_dlg, sfp->clickable_list_data);
                                                    
                                                    
  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  sfp->recheck_btn = PushButton (c, "Recheck", ReSUCProc);
  SetObjectExtra (sfp->recheck_btn, sfp, NULL);

  b = PushButton (c, "Expand All", SUCExpandAll);
  SetObjectExtra (b, sfp, NULL);
  b = PushButton (c, "Collapse All", SUCCollapseAll);
  SetObjectExtra (b, sfp, NULL);
    
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) sfp->clickable_list_dlg, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);
}


static void ConvertLocalToGeneralCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqIdPtr    sip, sip_new;
  ObjectIdPtr oip;
  CharPtr     db_end;
  DbtagPtr    dbtag;
  Boolean     found = TRUE;
  
  if (bsp == NULL) {
    return;
  }
  
  while (found) {
    found = FALSE;
    for (sip = bsp->id; sip != NULL && !found; sip = sip->next) {
      if (sip->choice == SEQID_LOCAL) {
        oip = (ObjectIdPtr) sip->data.ptrvalue;
        if (oip->str) {
          db_end = StringChr (oip->str, ':');
            if (db_end != NULL && !StringHasNoText (db_end + 1)) {
            dbtag = DbtagNew();
            dbtag->tag = ObjectIdNew();
            dbtag->tag->str = StringSave (db_end + 1);
            *db_end = 0;
            dbtag->db = StringSave (oip->str);
            sip_new = ValNodeNew(NULL);
            sip_new->choice = SEQID_GENERAL;
            sip_new->data.ptrvalue = dbtag;
            BioseqReplaceID (bsp, sip_new);
            sip_new = SeqIdFree (sip_new);
            found = TRUE;
          }
        }
      }
    }
  }
}

extern void ConvertLocalToGeneral(IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

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
  VisitBioseqsInSep (sep, NULL, ConvertLocalToGeneralCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void SetFeatureID (SeqFeatPtr sfp, Int4 id)
{
  SeqFeatXrefPtr    xref;
  ObjectIdPtr       oip;

  if (sfp == NULL) return;

  for (xref = sfp->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
  if (xref != NULL) {
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str != NULL) {
        oip->str = MemFree (oip->str);
      }
      oip->id = id;
    }
  } else {
    xref = SeqFeatXrefNew ();
    if (xref != NULL) {
      oip = ObjectIdNew ();
      if (oip != NULL) {
        oip->id = id;
        xref->id.choice = 3;
        xref->id.value.ptrvalue = (Pointer) oip;
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
    }
  }
}

static void SetCDSProductName (SeqFeatPtr cds, CharPtr product_name)
{
  BioseqPtr         prot_bsp;
  SeqFeatPtr        prot_feat;
  SeqMgrFeatContext context;
  Boolean           partial3, partial5;
  SeqEntryPtr       prot_sep;
  ProtRefPtr        prp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION
      || cds->product == NULL) return;

  prot_bsp = BioseqFindFromSeqLoc (cds->product);
  if (prot_bsp == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &context);
  if (prot_feat == NULL) {
    prot_sep = SeqMgrGetSeqEntryForData (prot_bsp);
    prot_feat = CreateNewFeature (prot_sep, NULL, SEQFEAT_PROT, NULL);
    prp = ProtRefNew ();
    ValNodeAddPointer (&(prp->name), 0, StringSave (product_name));
    prot_feat->data.value.ptrvalue = prp;
    prot_feat->location = SeqLocFree (prot_feat->location);
    prot_feat->location = CreateWholeInterval (prot_sep);
    SetSeqLocPartial (prot_feat->location, partial5, partial3);
    prot_feat->partial = (partial5 || partial3);
  } else {
    prp = prot_feat->data.value.ptrvalue;
    prp->name = ValNodeFreeData (prp->name);
    ValNodeAddPointer (&(prp->name), 0, StringSave (product_name));
  }
}


static CharPtr GetmRNAProductName (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return NULL;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp != NULL && rrp->type == RNA_TYPE_mRNA && rrp->ext.choice == 1) {
    return rrp->ext.value.ptrvalue;
  } else {
    return NULL;
  }
}


static void MakeCDSmRNAPairsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext context1, context2;
  SeqFeatPtr        cds, mRNA, first_mRNA, new_cds;
  Int2              loc_compare;
  SeqEntryPtr       sep, top_sep;
  ObjectIdPtr       oip_cds, oip_mrna;
  ValNodePtr        pair_list = NULL, vnp;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;

  for (cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context1);
       cds != NULL;
       cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, FEATDEF_CDS, &context1)) {
    first_mRNA = NULL;

    for (mRNA = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_mRNA, &context2);
         mRNA != NULL;
         mRNA = SeqMgrGetNextFeature (bsp, mRNA, SEQFEAT_RNA, FEATDEF_mRNA, &context2)) {
      if (context2.left > context1.right) {
        /* out of range */
        break;
      }
      if (StringNCmp (context2.label, context1.label, StringLen (context1.label)) != 0) {
        /* need root product name to match */
        continue;
      }
      loc_compare = SeqLocCompare (cds->location, mRNA->location);
      if (loc_compare != SLC_A_IN_B && loc_compare != SLC_A_EQ_B) {
        /* CDS must be equal to or inside mRNA */
        continue;
      }

      if (first_mRNA == NULL) {
        first_mRNA = mRNA;
        new_cds = cds;
      } else {
        new_cds = (SeqFeatPtr) AsnIoMemCopy (cds, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
        new_cds->product = SeqLocFree (new_cds->product);
        new_cds->xref = SeqFeatXrefFree (new_cds->xref);
        SeqFeatIdFree(&new_cds->id);
        CreateNewFeature (sep, NULL, SEQFEAT_CDREGION, new_cds);
        RetranslateOneCDS (new_cds, cds->idx.entityID, FALSE, FALSE);
      }
      ValNodeAddPointer (&pair_list, SEQFEAT_RNA, mRNA);
      ValNodeAddPointer (&pair_list, SEQFEAT_CDREGION, new_cds);
    }
  }
  top_sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  SeqMgrIndexFeatures (bsp->idx.entityID, NULL);
  AssignFeatureIDs (top_sep);
  vnp = pair_list;
  while (vnp != NULL) {
    mRNA = vnp->data.ptrvalue;
    vnp = vnp->next;
    if (vnp != NULL) {
      new_cds = vnp->data.ptrvalue;
      SetCDSProductName (new_cds, GetmRNAProductName(mRNA));
      oip_mrna = (ObjectIdPtr) mRNA->id.value.ptrvalue;
      oip_cds = (ObjectIdPtr) new_cds->id.value.ptrvalue;
      if (oip_mrna != NULL && oip_cds != NULL
          && oip_mrna->id > 0 && oip_cds->id > 0) {
        /* add xref to cds to point to mRNA */
        SetFeatureID (new_cds, oip_mrna->id);
        /* add xref to mrna to point to cds */
        SetFeatureID (mRNA, oip_cds->id);
      }
      vnp = vnp->next;
    }
  }
  pair_list = ValNodeFree (pair_list);
}


extern void MakeCDSmRNAPairs (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

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
  VisitBioseqsInSep (sep, NULL, MakeCDSmRNAPairsCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void ApplyNameRNACallback (BioseqPtr bsp, Pointer userdata)
{
  CharPtr    rnaName = (CharPtr) userdata;
  RnaRefPtr  rrp;
  SeqFeatPtr sfp;
  SeqEntryPtr sep;

  if (bsp == NULL || ISA_aa (bsp->mol) || StringHasNoText (rnaName)) {
    return;
  }
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;

  rrp = RnaRefNew ();
  rrp->type = RNA_TYPE_rRNA;
  rrp->ext.choice = 1;
  rrp->ext.value.ptrvalue = StringSave (rnaName);
  sfp = SeqFeatNew ();
  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = rrp;
  sfp->location = CreateWholeInterval (sep);
  SetSeqLocPartial (sfp->location, TRUE, TRUE);
  sfp->partial = TRUE;
  CreateNewFeature (sep, NULL, SEQFEAT_RNA, sfp);
}


static void ApplyNamedRNA (IteM i, CharPtr rnaName)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

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
  VisitBioseqsInSep (sep, rnaName, ApplyNameRNACallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


extern void Apply16SRNA (IteM i)
{
  ApplyNamedRNA (i, "16S ribosomal RNA");
}


extern void Apply23SRNA (IteM i)
{
  ApplyNamedRNA (i, "23S ribosomal RNA");
}


extern void Apply18SRNA (IteM i)
{
  ApplyNamedRNA (i, "18S ribosomal RNA");
}


extern void Apply28SRNA (IteM i)
{
  ApplyNamedRNA (i, "28S ribosomal RNA");
}


extern void Apply26SRNA (IteM i)
{
  ApplyNamedRNA (i, "26S ribosomal RNA");
}


extern void Apply12SRNA (IteM i)
{
  ApplyNamedRNA (i, "12S ribosomal RNA");
}


extern void ApplySmallRNA (IteM i)
{
  ApplyNamedRNA (i, "small subunit ribosomal RNA");
}


extern void ApplyLargeRNA (IteM i)
{
  ApplyNamedRNA (i, "large subunit ribosomal RNA");
}


typedef struct rnaits {
  FORM_MESSAGE_BLOCK
  ButtoN       partial_5;
  ButtoN       partial_3;
  DialoG       rna_list;
  ButtoN       leave_dlg_up;
  ButtoN       accept;
  /* for alignment coordinates */
  ButtoN       use_aln;
  PopuP        this_aln;

  BaseFormPtr       bfp;
  TaglistCallback   callbacks[3];
  Uint2             interval_types[3]; 
  Uint2             interval_widths[3];
  EnumFieldAssocPtr alists[3];

  ValNodePtr        prev_list;
  /* for alignment coordinates */
  SeqAlignPtr       salp_list;
  SeqAlignPtr       chosen_salp;

  /* for propagation */
  ValNodePtr        feature_list;
} RNAITSData, PNTR RNAITSPtr;

typedef enum {
  eRNA_ITS_INVALID = 0,
  eRNA_ITS_18S,
  eRNA_ITS_18S_SMALL,
  eRNA_ITS_SMALL,
  eRNA_ITS_ITS1,
  eRNA_ITS_58S,
  eRNA_ITS_ITS2,
  eRNA_ITS_28S,
  eRNA_ITS_28S_LARGE,
  eRNA_ITS_LARGE,
  eRNA_ITS_26S,
  eRNA_ITS_26S_LARGE,
  eRNA_ITS_25S,
  eRNA_ITS_25S_LARGE
} ERNA_ITS_type;

static ENUM_ALIST(rna_its_alist)
{" ",                                 eRNA_ITS_INVALID},
{"18S ribosomal RNA",                 eRNA_ITS_18S},
{"18S small subunit ribosomal RNA",   eRNA_ITS_18S_SMALL},
{"small subunit ribosomal RNA",       eRNA_ITS_SMALL},

{"internal transcribed spacer 1",     eRNA_ITS_ITS1},

{"5.8S ribosomal RNA",                eRNA_ITS_58S},

{"internal transcribed spacer 2",     eRNA_ITS_ITS2},

{"28S ribosomal RNA",                 eRNA_ITS_28S},
{"28S large subunit ribosomal RNA",   eRNA_ITS_28S_LARGE},
{"large subunit ribosomal RNA",       eRNA_ITS_LARGE},
{"26S ribosomal RNA",                 eRNA_ITS_26S},
{"26S large subunit ribosomal RNA",   eRNA_ITS_26S_LARGE},
{"25S ribosomal RNA",                 eRNA_ITS_25S},
{"25S large subunit ribosomal RNA",   eRNA_ITS_25S_LARGE},
END_ENUM_ALIST

static ENUM_ALIST(rna_its_popup_alist)
{" ",                                 eRNA_ITS_INVALID},
{"18S rRNA",                          eRNA_ITS_18S},
{"small subunit",                     eRNA_ITS_SMALL},
{"ITS1",                              eRNA_ITS_ITS1},
{"5.8S rRNA",                         eRNA_ITS_58S},
{"ITS2",                              eRNA_ITS_ITS2},
{"28S rRNA",                          eRNA_ITS_28S},
{"26S rRNA",                          eRNA_ITS_26S},
{"25S rRNA",                          eRNA_ITS_25S},
{"large subunit rRNA",                eRNA_ITS_LARGE},
{"18S small subunit rRNA",            eRNA_ITS_18S_SMALL},
{"28S large subunit rRNA",            eRNA_ITS_28S_LARGE},
{"26S large subunit rRNA",            eRNA_ITS_26S_LARGE},
{"25S large subunit rRNA",            eRNA_ITS_25S_LARGE},
END_ENUM_ALIST

static CharPtr RNANameFromRNAType (ERNA_ITS_type rna_type)
{
  EnumFieldAssocPtr eap = rna_its_alist;
  while (eap != NULL && eap->name != NULL && eap->value != rna_type) 
  {
    eap++;
  }
  return eap->name;
}

typedef struct rna_its_item {
  ERNA_ITS_type rna_type;
  Int4          left;
  Int4          right;
} RNA_ITS_ItemData, PNTR RNA_ITS_ItemPtr;

static void MarkRNAITSDiffs (ValNodePtr vnp_new, ValNodePtr vnp_old)
{
  RNA_ITS_ItemPtr rip_new, rip_old;
  while (vnp_new != NULL) {
    if (vnp_old == NULL) {
      vnp_new->choice = 1;
    } else {
      rip_new = vnp_new->data.ptrvalue;
      rip_old = vnp_old->data.ptrvalue;
      if (rip_old == NULL && rip_new == NULL) {
        /* no change */
      } else if (rip_old == NULL || rip_new == NULL) {
        vnp_new->choice = 1;
      } else if (rip_new->rna_type != rip_old->rna_type
                 || rip_new->left != rip_old->left
                 || rip_new->right != rip_old->right) {
        vnp_new->choice = 1;
      }
      vnp_old = vnp_old->next;
    }
    vnp_new = vnp_new->next;
  }
}

static RNA_ITS_ItemPtr TagListStringToRNAITSItem (CharPtr line)
{
  RNA_ITS_ItemPtr rip;
  CharPtr         txt;
  Int4            tmp;

  rip = (RNA_ITS_ItemPtr) MemNew (sizeof (RNA_ITS_ItemData));

  /* get interval values */
  txt = ExtractTagListColumn (line, 1);
  if (StringHasNoText (txt)) {
    rip->left = -1;
  } else {
    StrToLong (txt, &rip->left);
  }
  txt = MemFree (txt);
  txt = ExtractTagListColumn (line, 2);
  if (StringHasNoText (txt)) {
    rip->right = -1;
  } else {
    StrToLong (txt, &rip->right);
  }
  txt = MemFree (txt);

  txt = ExtractTagListColumn (line, 0);
  if (StringHasNoText (txt)) {
    rip->rna_type = eRNA_ITS_INVALID;
  } else {
    StrToLong (txt, &tmp);
    rip->rna_type = (ERNA_ITS_type) tmp;
  }
  txt = MemFree (txt);

  return rip;
}

static CharPtr RNAITSItemToTagListString (RNA_ITS_ItemPtr rip)
{
  Char new_line[100];
  Char left_str[15];
  Char right_str[15];

  if (rip == NULL) return NULL;
  if (rip->left > -1) {
    sprintf (left_str, "%d", rip->left);
  } else {
    sprintf (left_str, "");
  }
  if (rip->right > -1) {
    sprintf (right_str, "%d", rip->right);
  } else {
    sprintf (right_str, "");
  }
  sprintf (new_line, "%d\t%s\t%s\n", rip->rna_type, left_str, right_str);
  return StringSave (new_line);
}

static ERNA_ITS_type NextRNAITSChoice (ERNA_ITS_type rna_type)
{
  ERNA_ITS_type next_type = eRNA_ITS_INVALID;

  switch (rna_type) {
    case eRNA_ITS_INVALID:
      next_type = eRNA_ITS_INVALID;
      break;
    case eRNA_ITS_18S:
    case eRNA_ITS_18S_SMALL:
    case eRNA_ITS_SMALL:
      next_type = eRNA_ITS_ITS1;
      break;
    case eRNA_ITS_ITS1:
      next_type = eRNA_ITS_58S;
      break;
    case eRNA_ITS_58S:
      next_type = eRNA_ITS_ITS2;
      break;
    case eRNA_ITS_ITS2:
      next_type = eRNA_ITS_28S;
      break;
    case eRNA_ITS_28S:
    case eRNA_ITS_28S_LARGE:
    case eRNA_ITS_LARGE:
    case eRNA_ITS_26S:
    case eRNA_ITS_26S_LARGE:
    case eRNA_ITS_25S:
    case eRNA_ITS_25S_LARGE:
      next_type = eRNA_ITS_INVALID;
      break;
  }
  return next_type;
}

static Boolean OkNextRNAITSChoice (ERNA_ITS_type rna_type, ERNA_ITS_type next_rna_type)
{
  Boolean is_ok = FALSE;

  switch (rna_type) {
    case eRNA_ITS_INVALID:
      is_ok = (next_rna_type == eRNA_ITS_INVALID);
      break;
    case eRNA_ITS_18S:
    case eRNA_ITS_18S_SMALL:
    case eRNA_ITS_SMALL:
      is_ok = (next_rna_type == eRNA_ITS_ITS1);
      break;
    case eRNA_ITS_ITS1:
      is_ok = (next_rna_type == eRNA_ITS_58S);
      break;
    case eRNA_ITS_58S:
      is_ok = (next_rna_type == eRNA_ITS_ITS2);
      break;
    case eRNA_ITS_ITS2:
      if (next_rna_type == eRNA_ITS_28S
          || next_rna_type == eRNA_ITS_28S_LARGE
          || next_rna_type == eRNA_ITS_LARGE
          || next_rna_type == eRNA_ITS_26S
          || next_rna_type == eRNA_ITS_26S_LARGE
          || next_rna_type == eRNA_ITS_25S
          || next_rna_type == eRNA_ITS_25S_LARGE) {
        is_ok = TRUE;
      } else {
        is_ok = FALSE;
      }
      break;
    case eRNA_ITS_28S:
    case eRNA_ITS_28S_LARGE:
    case eRNA_ITS_LARGE:
    case eRNA_ITS_26S:
    case eRNA_ITS_26S_LARGE:
    case eRNA_ITS_25S:
    case eRNA_ITS_25S_LARGE:
      is_ok = (next_rna_type == eRNA_ITS_INVALID);
      break;
  }
  return is_ok;
}

static Pointer TagListToRNAITSList (DialoG d)
{
  TagListPtr  tlp;
  ValNodePtr  rip_list = NULL, vnp;
  RNA_ITS_ItemPtr rip;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return NULL;
  vnp = tlp->vnp;
  while (vnp != NULL) {
    rip = TagListStringToRNAITSItem (vnp->data.ptrvalue);
    if (rip != NULL) {
      ValNodeAddPointer (&rip_list, 0, rip);
    }
    vnp = vnp->next;
  }
  return (Pointer) rip_list;
}

static void RNAITSListToTagList (DialoG d, Pointer data)
{
  TagListPtr  tlp;
  ValNodePtr  rip_list, vnp;
  RNA_ITS_ItemPtr rip;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;
  rip_list = (ValNodePtr) data;
  tlp->vnp = ValNodeFreeData (tlp->vnp);
  vnp = rip_list;
  while (vnp != NULL) {
    rip = (vnp->data.ptrvalue);
    if (rip != NULL) {
      ValNodeAddPointer (&tlp->vnp, 0, RNAITSItemToTagListString (rip));
    }
    vnp = vnp->next;
  }
  SendMessageToDialog (d, VIB_MSG_REDRAW);
}

static void RNA_ITS_ChoiceCallback (Pointer userdata)
{
  RNAITSPtr   rp;
  TagListPtr  tlp;
  ValNodePtr  vnp, prev_vnp;
  RNA_ITS_ItemPtr rip, rip_prev = NULL;
  Boolean         any_changes = FALSE, any_changes_this;

  rp = (RNAITSPtr) userdata;
  if (rp == NULL) return;

  tlp = (TagListPtr) GetObjectExtra (rp->rna_list);
  if (tlp == NULL) {
    return;
  }

  for (vnp = tlp->vnp, prev_vnp = NULL;
       vnp != NULL; 
       prev_vnp = vnp, vnp = vnp->next) {
    rip = TagListStringToRNAITSItem (vnp->data.ptrvalue);
    if (rip_prev == NULL) {
      if (rip->left == -1) {
        any_changes = TRUE;
        rip->left = 1;
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        vnp->data.ptrvalue = RNAITSItemToTagListString (rip);
      }
    } else {
      any_changes_this = FALSE;
      if (!OkNextRNAITSChoice(rip_prev->rna_type, rip->rna_type)) {
        /* fix type for this line */
        any_changes_this = TRUE;
        rip->rna_type = NextRNAITSChoice (rip_prev->rna_type);
      }
      /* fix interval for this line */
      if (rip->left == -1 && rip_prev->right != -1) {
        any_changes_this = TRUE;
        rip->left = rip_prev->right + 1;
      }
      if (any_changes_this) {
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        vnp->data.ptrvalue = RNAITSItemToTagListString (rip);
        any_changes = TRUE;
      }
      /* fix interval for previous line */
      if (rip_prev->right == -1 && rip->left > 0) {
        any_changes = TRUE;
        rip_prev->right = rip->left - 1;
        prev_vnp->data.ptrvalue = MemFree (prev_vnp->data.ptrvalue);
        prev_vnp->data.ptrvalue = RNAITSItemToTagListString (rip_prev);
      }
    }
    rip_prev = MemFree (rip_prev);
    rip_prev = rip;
  }
  /* fill out the chain */
  while (rip_prev != NULL && NextRNAITSChoice(rip_prev->rna_type) != eRNA_ITS_INVALID) {
    rip = (RNA_ITS_ItemPtr) MemNew (sizeof (RNA_ITS_ItemData));
    rip->rna_type = NextRNAITSChoice(rip_prev->rna_type);
    if (rip_prev->right > -1) {
      rip->left = rip_prev->right + 1;
    } else {
      rip->left = -1;
    }
    rip->right = -1;
    ValNodeAddPointer (&tlp->vnp, 0, RNAITSItemToTagListString (rip));
    any_changes = TRUE;
    rip_prev = MemFree (rip_prev);
    rip_prev = rip;
  }
  rip_prev = MemFree (rip_prev);
  if (any_changes) {
    SendMessageToDialog (rp->rna_list, VIB_MSG_REDRAW);
  }
  rp->prev_list = ValNodeFreeData (rp->prev_list);
  rp->prev_list = TagListToRNAITSList (rp->rna_list);
}

extern void RNA_ITS_IntervalCallback (Pointer data)
{
  RNAITSPtr   rp;
  TagListPtr  tlp;
  ValNodePtr  vnp, prev_vnp, new_list;
  RNA_ITS_ItemPtr rip, rip_prev = NULL;
  Boolean         any_changes = FALSE;

  rp = (RNAITSPtr) data;
  if (rp == NULL) return;

  tlp = (TagListPtr) GetObjectExtra (rp->rna_list);
  if (tlp == NULL) {
    return;
  }
  new_list = TagListToRNAITSList (rp->rna_list);
  MarkRNAITSDiffs (new_list, rp->prev_list);
  for (vnp = new_list, prev_vnp = NULL;
       vnp != NULL;
       prev_vnp = vnp, vnp = vnp->next) {
    rip = vnp->data.ptrvalue;
    if (rip_prev == NULL) {
      if (rip->left == -1) {
        rip->left = 0;
        any_changes = TRUE;
      }
    } else {
      /* fix interval for this line */
      if (rip->left == -1 && rip_prev->right != -1) {
        any_changes = TRUE;
        rip->left = rip_prev->right + 1;
      }
      /* fix interval for previous line */
      if (rip_prev->right == -1 && rip->left > 0) {
        any_changes = TRUE;
        rip_prev->right = rip->left - 1;
      }
      if (vnp->choice == 1) {
        if (rip->left != -1 && rip_prev->right != rip->left - 1) {
          any_changes = TRUE;
          rip_prev->right = rip->left - 1;
        }
      } else if (prev_vnp->choice == 1) {
        if (rip_prev->right != -1 && rip->left != rip_prev->right + 1) {
          any_changes = TRUE;
          rip->left = rip_prev->right + 1;
        }
      }
    }
    rip_prev = rip;
  }
  if (any_changes) {
    RNAITSListToTagList (rp->rna_list, new_list);
    SendMessageToDialog (rp->rna_list, VIB_MSG_REDRAW);
  }
  new_list = ValNodeFreeData (new_list);
  rp->prev_list = ValNodeFreeData (rp->prev_list);
  rp->prev_list = TagListToRNAITSList (rp->rna_list);
}

static Int4 GetAlnRowForBsp (BioseqPtr bsp, SeqAlignPtr salp)
{
  Int4     aln_row;
  SeqIdPtr sip;
  
  if (bsp == NULL || salp == NULL)
  {
    return 0;
  }
  
  for (aln_row = 1; aln_row <= salp->dim; aln_row++) 
  {
    sip = AlnMgr2GetNthSeqIdPtr(salp, aln_row);
    if (SeqIdIn (sip, bsp->id))
    {
      return aln_row;
    }
  }
  return 0;
}

static Int4 ChooseEndpoint (Int4 bsp_len, Int4 pos, Boolean is_left, SeqAlignPtr salp, Int4 aln_row)
{
  Int4 tmp;
  Int4 aln_len;

  aln_len = SeqAlignLength (salp);
  if (pos < 0) {
    if (is_left) {
      pos = 0;
    } else {
      pos = bsp_len - 1;
    }
  } else if (salp == NULL) { 
    if (pos > bsp_len - 1) {
      if (is_left) {
        pos = -1;
      } else {
        pos = bsp_len - 1;
      }
    }
  } else if (salp != NULL) {      
    tmp = AlnMgr2MapSeqAlignToBioseq (salp, pos, aln_row);
    if (is_left) {
      while (tmp < 0 && pos < aln_len) {
        pos++;
        tmp = AlnMgr2MapSeqAlignToBioseq (salp, pos, aln_row);
      }
      if (tmp < 0) {
        pos = -1;
      } else {
        pos = tmp;
      }
    } else {
      while (tmp < 0 && pos > -1) {
        pos--;
        tmp = AlnMgr2MapSeqAlignToBioseq (salp, pos, aln_row);
      }
      if (tmp < 0) {
        pos = -1;
      } else {
        pos = tmp;
      }
    }
  }
  return pos;
}


static void MakeRNA_ITSChainForBioseq (BioseqPtr bsp, Pointer userdata)
{
  RNAITSPtr          rp;
  ValNodePtr         vnp;
  RNA_ITS_ItemPtr    rip;
  Int4               last_untranslated_right = 0, this_left, this_right;
  RnaRefPtr          rrp;
  SeqFeatPtr         sfp;
  SeqEntryPtr        sep;
  Boolean            partial5, partial3 = FALSE;
  Int4               aln_row = -1;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) return;
  rp = (RNAITSPtr) userdata;

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;

  if (rp->chosen_salp != NULL) {
    aln_row = GetAlnRowForBsp (bsp, rp->chosen_salp);
    if (aln_row < 1) {
      return;
    }
  }

  partial5 = GetStatus (rp->partial_5);
  last_untranslated_right = 0;
  for (vnp = rp->prev_list; vnp != NULL; vnp = vnp->next) {
    rip = vnp->data.ptrvalue;

    /* stop if we are at the end of the chain */
    if (rip->rna_type == eRNA_ITS_INVALID) {
      break;
    }

    if (rip->left < 0) {      
      this_left = last_untranslated_right + 1;
    } else {
      this_left = rip->left - 1;
    } 

    this_left = ChooseEndpoint (bsp->length, this_left, TRUE, rp->chosen_salp, aln_row);

    if (rip->right < 0) {
      if (rp->chosen_salp == NULL) {
        this_right = bsp->length - 1;
      } else {
        this_right = AlnMgr2GetAlnLength (rp->chosen_salp, TRUE) - 1;
      }    
    } else {
      this_right = rip->right - 1;
    }
    last_untranslated_right = this_right;
    this_right = ChooseEndpoint (bsp->length, this_right, FALSE, rp->chosen_salp, aln_row);

    /* if using alignment coordinates and left is not yet available, continue */
    if (this_right < 0 && this_left > -1 && rp->chosen_salp != NULL) {
      continue;
    }

    /* stop if the chain extends past the length of this Bioseq */
    if (this_left > bsp->length - 1 || this_left < 0 || this_right < 0 || this_right <= this_left) {
      break;
    }

    rrp = RnaRefNew ();
    if (rip->rna_type == eRNA_ITS_ITS1 || rip->rna_type == eRNA_ITS_ITS2) {
      rrp->type = RNA_TYPE_other;
    } else {
      rrp->type = RNA_TYPE_rRNA;
    }
    rrp->ext.choice = 1;

    rrp->ext.value.ptrvalue = StringSave (RNANameFromRNAType(rip->rna_type));
    sfp = SeqFeatNew ();
    sfp->data.choice = SEQFEAT_RNA;
    sfp->data.value.ptrvalue = rrp;
    /* if last endpoint not set, use end of sequence.
     * if last endpoint past end of sequence, use end of sequence.
     */
    sfp->location = SeqLocIntNew (this_left, this_right, Seq_strand_plus,
                                  SeqIdDup (SeqIdFindBest (bsp->id, 0)));
    if ((vnp->next == NULL || this_right == bsp->length - 1) && GetStatus (rp->partial_3)) {
      partial3 = TRUE;
    }
    SetSeqLocPartial (sfp->location, partial5, partial3);
    sfp->partial = partial5 || partial3;
    CreateNewFeature (sep, NULL, SEQFEAT_RNA, sfp);
    ValNodeAddPointer (&(rp->feature_list), OBJ_SEQFEAT, sfp);
    partial5 = FALSE;
  }
}

/* This function will propagate the RNA features pointed to by feature_list
 * from aln_row in salp to the other rows in salp. 
*/
static void PropagateRNAList (ValNodePtr feature_list, Int4 aln_row, SeqAlignPtr salp)
{
  ValNodePtr seq_for_prop = NULL, vnp;
  Int4       i;
  SeqIdPtr   sip;
  Boolean    warned_about_master = FALSE;

  if (feature_list == NULL || aln_row < 1 || salp == NULL) return;

  /* get list of sequences to propagate to */
  for (i = 1; i <= salp->dim; i++) {
    if (i != aln_row) {
      ValNodeAddPointer (&seq_for_prop, 0, AlnMgr2GetNthSeqIdPtr(salp, i));
    }
  }

  /* propagate each feature */
  for (vnp = feature_list; vnp != NULL; vnp = vnp->next) {
    PropagateOneFeat (vnp->data.ptrvalue, TRUE, FALSE, FALSE, FALSE, FALSE, seq_for_prop, &warned_about_master);
  }

  /* free sequence IDs in seq_for_prop */
  for (vnp = seq_for_prop; vnp != NULL; vnp = vnp->next) {
    sip = (SeqIdPtr) vnp->data.ptrvalue;
    sip = SeqIdFree (sip);
  }
  seq_for_prop = ValNodeFree (seq_for_prop);

}


static void MakeRNA_ITSChain (ButtoN b)
{
  RNAITSPtr          rp;
  SeqEntryPtr        sep;
  BioseqPtr          bsp;
  Boolean            apply_to_all = FALSE, propagate = FALSE;
  SeqAlignPtr        next_salp, propagate_salp = NULL;
  Int4               salp_index, aln_row = -1;

  rp = (RNAITSPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  sep = GetTopSeqEntryForEntityID (rp->input_entityID);
  rp->prev_list = ValNodeFreeData (rp->prev_list);
  rp->prev_list = DialogToPointer (rp->rna_list);

  bsp = GetBioseqViewTarget (rp->bfp);
  if (GetStatus (rp->use_aln)) {
    apply_to_all = TRUE;    
    if (bsp != NULL && Message (MSG_YN, "Apply to all Bioseqs?") == ANS_NO) {
      apply_to_all = FALSE;
    }
  } else if (bsp == NULL) {
    if (Message (MSG_OKC, "You are not viewing a single Bioseq.  Apply to all Bioseqs?") == ANS_CANCEL) {
      return;
    } else {
      apply_to_all = TRUE;
    }
  } else {
    propagate_salp = rp->salp_list;
    while (aln_row < 0 && propagate_salp != NULL) {
      aln_row = GetAlnRowForBsp (bsp, propagate_salp);
      if (aln_row < 0) {
        propagate_salp = propagate_salp->next;
      }
    }
    if (aln_row > 0  && propagate_salp != NULL && Message (MSG_YN, "Propagate to other sequences?") == ANS_YES) {
      propagate = TRUE;
    }
  }

  rp->chosen_salp = NULL;
  next_salp = NULL;
  if (GetStatus (rp->use_aln)) {
    rp->chosen_salp = rp->salp_list;
    if (rp->this_aln != NULL) {
      salp_index = GetValue (rp->this_aln);
      while (salp_index > 1 && rp->chosen_salp != NULL) {
        rp->chosen_salp = rp->chosen_salp->next;
        salp_index--;
      }
      if (salp_index > 1) {
        rp->chosen_salp = NULL;
      }
    }

    /* if choosing only one sequence, make sure it's in alignment */
    if (bsp != NULL && !apply_to_all && GetAlnRowForBsp (bsp, rp->chosen_salp) < 1) {
      Message (MSG_ERROR, "This sequence is not in the alignment!");
      return;
    }

    /* temporarily break alignment chain */
    if (rp->chosen_salp != NULL) {
      next_salp = rp->chosen_salp->next;
      rp->chosen_salp->next = NULL;
    }
  }


  WatchCursor ();
  Update ();

  /* rp->feature list will contain a list of all the features created */
  rp->feature_list = NULL;
  if (apply_to_all) {
    VisitBioseqsInSep (sep, rp, MakeRNA_ITSChainForBioseq);
  } else {
    MakeRNA_ITSChainForBioseq (bsp, rp);
    if (propagate) {
      PropagateRNAList (rp->feature_list, aln_row, propagate_salp);
    }
  }
  rp->feature_list = ValNodeFree (rp->feature_list);


  /* restore alignment chain if it was broken */
  if (rp->chosen_salp != NULL) {
    rp->chosen_salp->next = next_salp;
  }
  
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (rp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rp->input_entityID, 0, 0);
  if (!GetStatus (rp->leave_dlg_up)) {
    Remove (rp->form);
  }
}

static WindoW RNAITSWindow = NULL;

static void CleanupRNA_ITS (GraphiC g, VoidPtr data)

{
  RNAITSPtr          rp;
  SeqAlignPtr        salp;

  rp = (RNAITSPtr) data;
  if (rp != NULL) {
    ValNodeFree (rp->prev_list);
    /* free the alignment list we created for alignment coordinates */
    while (rp->salp_list != NULL) {
      salp = rp->salp_list->next;
      rp->salp_list->next = NULL;
      rp->salp_list = SeqAlignFree (rp->salp_list);
      rp->salp_list = salp;
    }
    ObjMgrFreeUserData (rp->input_entityID, rp->procid, rp->proctype, rp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static Int2 LIBCALLBACK RNAITSFormMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW         currentport,
                 temport;
  OMUserDataPtr  omudp;
  RNAITSPtr      rp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  rp = (RNAITSPtr) omudp->userdata.ptrvalue;
  if (rp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (rp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (rp->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (rp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (rp->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}

static void UseAlnCallback (ButtoN b)
{
  RNAITSPtr      rp = NULL;

  rp = (RNAITSPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  if (rp->use_aln != NULL && rp->this_aln != NULL) {
    if (GetStatus (rp->use_aln)) {
      Enable (rp->this_aln);
    } else {
      Disable (rp->this_aln);
    }
  }
}

extern void ApplyRNA_ITS (IteM i)
{
  GrouP              c;
  RNAITSPtr          rp;
  GrouP              h, g, aln_group = NULL;
  WindoW             w;
  BaseFormPtr        bfp;
  OMUserDataPtr      omudp;
  SeqEntryPtr        sep;
  SeqAlignPtr        salp;
  Char               str[31];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  rp = (RNAITSPtr) MemNew (sizeof (RNAITSData));
  if (rp == NULL)
    return;

  w = FixedWindow (-50, -33, -10, -10, "Apply RNA-ITS",
		   StdCloseWindowProc);
  SetObjectExtra (w, rp, CleanupRNA_ITS);
  rp->form = (ForM) w;
  rp->input_entityID = bfp->input_entityID;
  rp->bfp = bfp;

  /* register to receive update messages */
  rp->userkey = OMGetNextUserKey ();
  rp->procid = 0;
  rp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (rp->input_entityID, rp->procid, rp->proctype, rp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) rp;
    omudp->messagefunc = RNAITSFormMsgFunc;
  }


  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  rp->partial_5 = CheckBox (g, "5' Partial", NULL);
  rp->partial_3 = CheckBox (g, "3' Partial", NULL);
  SetStatus (rp->partial_5, TRUE);
  SetStatus (rp->partial_3, TRUE);

  /* if there are alignments, allow the user to specify alignment coordinates */
  sep = GetTopSeqEntryForEntityID (rp->input_entityID);
  VisitAnnotsInSep (sep, &(rp->salp_list), GetAlignmentsInSeqEntryCallback);
  if (rp->salp_list == NULL) {
    rp->use_aln = NULL;
    rp->this_aln = NULL;
  } else {
    aln_group = HiddenGroup (h, 2, 0, NULL);
    rp->use_aln = CheckBox (aln_group, "Use alignment coordinates", UseAlnCallback);
    SetObjectExtra (rp->use_aln, rp, NULL);
    if (rp->salp_list->next != NULL) {
      rp->this_aln = PopupList (aln_group, TRUE, NULL);
      salp = rp->salp_list;
      while (salp != NULL)
      {
        CreateSeqAlignLabel (salp, str, sizeof (str) - 1);
        PopupItem (rp->this_aln, str);
        salp = salp->next;
      }
      SetValue (rp->this_aln, 1); 
      Disable (rp->this_aln);
    } else {
      rp->this_aln = NULL;
    }
  }


  /* set up callbacks */
  rp->callbacks[0] = RNA_ITS_ChoiceCallback;
  rp->callbacks[1] = RNA_ITS_IntervalCallback;
  rp->callbacks[2] = RNA_ITS_IntervalCallback;

  /* set up inteval types */
  rp->interval_types[0] = TAGLIST_POPUP;
  rp->interval_types[1] = TAGLIST_TEXT;
  rp->interval_types[2] = TAGLIST_TEXT;

  /* set up interval widths */
  rp->interval_widths[0] = 0;
  rp->interval_widths[1] = 5;
  rp->interval_widths[2] = 5;

  /* set up alist for RNA type */
  rp->alists[0] = rna_its_popup_alist;
  rp->alists[1] = NULL;
  rp->alists[2] = NULL;

  rp->rna_list = CreateTagListDialogExEx (h, 4, 3, 2,
                                          rp->interval_types, rp->interval_widths, rp->alists,
                                          TRUE, FALSE, RNAITSListToTagList, TagListToRNAITSList,
                                          rp->callbacks, rp, FALSE);

  rp->prev_list = NULL;
  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 3, 0, NULL);
  rp->accept = DefaultButton (c, "Accept", MakeRNA_ITSChain);
  SetObjectExtra (rp->accept, rp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  rp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) g,
                              (HANDLE) rp->rna_list,
                              (HANDLE) c, NULL);


  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (rp->accept);
  Update ();
}


/* CDS-mRNA Link Tool */
typedef struct cdsmrnalinktool {
  FORM_MESSAGE_BLOCK
  GrouP edit_grp;
  GrouP report_grp;
  DoC left_doc;
  DoC right_doc;
  DoC report_doc;
  DoC report_title_doc;

  BaseFormPtr       bfp;
  SeqFeatPtr        selected_cds;
  ValNodePtr        cds_list;
  ValNodePtr        mrna_list;
  Int2       lineheight;  

  Nlm_ParData cdsParFmt;
  Nlm_ParData mrnaParFmt;
  Nlm_ParData reportParFmt;
  Nlm_ColData cdsColFmt [2];
  Nlm_ColData mrnaColFmt [3];
  Nlm_ColData reportColFmt [4];
} CDSmRNALinkToolData, PNTR CDSmRNALinkToolPtr;

static void CleanupCDSmRNALinkTool (GraphiC g, VoidPtr data)

{
  CDSmRNALinkToolPtr          tp;

  tp = (CDSmRNALinkToolPtr) data;
  if (tp != NULL) {
    ObjMgrFreeUserData (tp->input_entityID, tp->procid, tp->proctype, tp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static Int2 LIBCALLBACK CDSmRNALinkToolFormMsgFunc (OMMsgStructPtr ommsp)
{
  OMUserDataPtr  omudp;
  CDSmRNALinkToolPtr      tp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  tp = (CDSmRNALinkToolPtr) omudp->userdata.ptrvalue;
  if (tp == NULL) return OM_MSG_RET_ERROR;

  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (tp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (tp->form);	
          break;
      default:
          break;
  }
  return OM_MSG_RET_OK;
}

static ValNodePtr GetLinkedFeatureList (SeqFeatPtr sfp)
{
  Char            buf [32];
  ObjectIdPtr     oip;
  SeqFeatPtr      link_sfp;
  CharPtr         str = NULL;
  SeqFeatXrefPtr  xref;
  ValNodePtr      link_list = NULL;

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
      link_sfp = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, str, NULL, NULL);
      if (link_sfp != NULL) {
        ValNodeAddPointer (&link_list, OBJ_SEQFEAT, link_sfp);
      }
    }
  }

  return link_list;
}

static int LIBCALLBACK SortVnpByLinkStatus (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  SeqFeatPtr  cds1, cds2;
  ValNodePtr  link1, link2;
  Int4        rval = 0;
  SeqMgrFeatContext fcontext1, fcontext2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;
  
  cds1 = vnp1->data.ptrvalue;
  cds2 = vnp2->data.ptrvalue;

  link1 = GetLinkedFeatureList (cds1);
  link2 = GetLinkedFeatureList (cds2);

  if (link1 != NULL && link2 == NULL) {
    rval = 1;
  } else if (link1 == NULL && link2 != NULL) {
    rval = -1;
  } else {
    SeqMgrGetDesiredFeature (cds1->idx.entityID, NULL, cds1->idx.itemID, 0, cds1, &fcontext1);
    SeqMgrGetDesiredFeature (cds2->idx.entityID, NULL, cds2->idx.itemID, 0, cds2, &fcontext2);
    if (fcontext1.left > fcontext2.left) {
      rval = 1;
    } else if (fcontext1.left < fcontext2.left) {
      rval = -1;
    } else if (fcontext1.itemID > fcontext2.itemID) {
      rval = 1;
    } else if (fcontext2.itemID < fcontext2.itemID) {
      rval = -1;
    } else {
      rval = 0;
    }
  }
  link1 = ValNodeFree (link1);
  link2 = ValNodeFree (link2);

  return rval;
}

static int LIBCALLBACK SortVnpBymRNAPos (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  SeqFeatPtr  mrna1, mrna2;
  Int4        rval = 0;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;
  
  mrna1 = vnp1->data.ptrvalue;
  mrna2 = vnp2->data.ptrvalue;

  if (vnp1->choice > 0 && vnp2->choice == 0) {
    rval = -1;
  } else if (vnp1->choice == 0 && vnp2->choice > 0) {
    rval = 1;
  } else {
    rval = SortVnpByLinkStatus (ptr1, ptr2);
  }

  return rval;
}

static CharPtr GetLinkLabel (SeqFeatPtr sfp)
{
  CharPtr location, label;
  SeqMgrFeatContext context;

  if (sfp == NULL) return NULL;

  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &context);
  location = SeqLocPrintUseBestID (sfp->location);

  label = (CharPtr) MemNew (sizeof (Char) * (StringLen (context.label) + StringLen (location) + 3));
  sprintf (label, "%s\t%s\n", context.label, location);
  return label;
}

static void PopulateLinkReportDoc (CDSmRNALinkToolPtr tp)
{
  CharPtr cds_desc, mrna_desc, whole_line, tmp_desc;
  ValNodePtr vnp, vnp_m, mrnas;
  CharPtr    cp;
  CharPtr    no_mrna = "\t\n";
  Int4       col_width;
  RecT       r;

  if (tp == NULL) return;

  Reset (tp->report_doc);
  Reset (tp->report_title_doc);
  ObjectRect (tp->report_doc, &r);
  InsetRect (&r, 4, 4);

  col_width = (r.right - r.left) / 4;
  tp->reportColFmt[0].pixWidth = col_width;
  tp->reportColFmt[1].pixWidth = col_width;
  tp->reportColFmt[2].pixWidth = col_width;
  tp->reportColFmt[3].pixWidth = col_width;

  AppendText (tp->report_title_doc, "CDS product\tCDS location\tLinked mRNA product\tLinked mRNA Location\n",
              &(tp->reportParFmt), tp->reportColFmt, programFont);

  tp->cds_list = ValNodeSort (tp->cds_list, SortVnpByLinkStatus);

  for (vnp = tp->cds_list; vnp != NULL; vnp = vnp->next) {
    cds_desc = GetLinkLabel(vnp->data.ptrvalue);
    /* replace carriage return with tab */
    cp = StringChr (cds_desc, '\n');
    if (cp != NULL) {
      *cp = '\t';
    }
    mrnas = GetLinkedFeatureList (vnp->data.ptrvalue);
    if (mrnas == NULL) {
      whole_line = (CharPtr) MemNew (sizeof (Char) * (StringLen (cds_desc) + StringLen (no_mrna) + 1));
      sprintf (whole_line, "%s%s", cds_desc, no_mrna);
      AppendText (tp->report_doc, whole_line, &(tp->reportParFmt), tp->reportColFmt, programFont);
      whole_line = MemFree (whole_line);
    } else {
      tmp_desc = cds_desc;
      for (vnp_m = mrnas; vnp_m != NULL; vnp_m = vnp_m->next) {
        mrna_desc = GetLinkLabel(vnp_m->data.ptrvalue);
        whole_line = (CharPtr) MemNew (sizeof (Char) * (StringLen (cds_desc) + StringLen (mrna_desc) + 1));
        sprintf (whole_line, "%s%s", tmp_desc, mrna_desc);
        AppendText (tp->report_doc, whole_line, &(tp->reportParFmt), tp->reportColFmt, programFont);
        whole_line = MemFree (whole_line);
        mrna_desc = MemFree (mrna_desc);
        tmp_desc = "\t\t";
      }
    }
    cds_desc = MemFree (cds_desc);
  }
  InvalRect (&r);
  ObjectRect (tp->report_title_doc, &r);
  InvalRect (&r);
}

static void PopulateLinkLeftDoc (CDSmRNALinkToolPtr tp)
{
  RecT       r;
  ValNodePtr vnp;
  CharPtr    desc;
  Int4       col_width;
  Int2       i = 1, item_num;
  Int4       startsAt;

  if (tp == NULL) return;

  /* first, reorder the coding regions so that coding regions with links
   * appear last.
   */
  tp->cds_list = ValNodeSort (tp->cds_list, SortVnpByLinkStatus);

  ObjectRect (tp->left_doc, &r);
  InsetRect (&r, 4, 4);
  
  col_width = (r.right - r.left) / 2;
  tp->cdsColFmt[0].pixWidth = col_width;
  tp->cdsColFmt[1].pixWidth = col_width;

  Reset (tp->left_doc);

  for (vnp = tp->cds_list, i = 1; vnp != NULL; vnp = vnp->next, i++) {

    desc = GetLinkLabel(vnp->data.ptrvalue);
    AppendText (tp->left_doc, desc, &(tp->cdsParFmt), tp->cdsColFmt, programFont);
    desc = MemFree (desc);
    if (vnp->data.ptrvalue == tp->selected_cds) {
      item_num = i;
    }
  }

  GetItemParams4 (tp->left_doc, item_num, &startsAt, NULL, NULL, NULL, NULL);
  SetScrlParams4 (tp->left_doc, startsAt);

  InvalRect (&r);  
}

static void PopulateLinkRightDoc (CDSmRNALinkToolPtr tp)
{
  RecT       r;
  ValNodePtr vnp;
  CharPtr    desc, tmp_desc;
  Uint1      old_choice;
  Int4       col_width;

  if (tp == NULL) return;

  Reset (tp->right_doc);
  ObjectRect (tp->right_doc, &r);
  InsetRect (&r, 4, 4);
  col_width = (r.right - r.left - tp->lineheight) / 2;
  tp->mrnaColFmt[0].pixWidth = tp->lineheight;
  tp->mrnaColFmt[1].pixWidth = col_width;
  tp->mrnaColFmt[2].pixWidth = col_width;

  for (vnp = tp->mrna_list; vnp != NULL; vnp = vnp->next) {
    old_choice = vnp->choice;
    vnp->choice = OBJ_SEQFEAT;
    desc = GetLinkLabel(vnp->data.ptrvalue);
    tmp_desc = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc) + 2));
    sprintf (tmp_desc, "\t%s", desc);
    AppendText (tp->right_doc, tmp_desc, &(tp->mrnaParFmt), tp->mrnaColFmt, programFont);
    desc = MemFree (desc);
    tmp_desc = MemFree (tmp_desc);
    vnp->choice = old_choice;
  }
  InvalRect (&r);  
}


static void PopulateLinkTool (CDSmRNALinkToolPtr tp)
{
  if (tp == NULL) return;

  /* first, reorder the coding regions so that coding regions with links
   * appear last.
   */
  tp->cds_list = ValNodeSort (tp->cds_list, SortVnpByLinkStatus);

  ResetClip();
  if (tp->selected_cds == NULL) {
    PopulateLinkReportDoc (tp);
    Show (tp->report_grp);
    Hide (tp->edit_grp);
  } else {
    PopulateLinkLeftDoc (tp);
    PopulateLinkRightDoc (tp);
    Show (tp->edit_grp);
    Hide (tp->report_grp);
  }
  Update ();
}

static void RefreshLinkBtn (ButtoN b)
{
  PopulateLinkTool ((CDSmRNALinkToolPtr)GetObjectExtra (b));
}

static Boolean HighlightSelectedCDS (DoC doc, Int2 item, Int2 row, Int2 col)
{
  CDSmRNALinkToolPtr tp;
  ValNodePtr         vnp;
  Int4               cds_num;
  
  tp = (CDSmRNALinkToolPtr) GetObjectExtra (doc);
  if (tp == NULL) return FALSE;
  
  for (vnp = tp->cds_list, cds_num = 1; vnp != NULL && cds_num != item; vnp = vnp->next, cds_num++)
  {
  }
  if (vnp != NULL && cds_num == item && vnp->data.ptrvalue == tp->selected_cds) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}

static SeqFeatPtr CDSFromLinkDocItem (CDSmRNALinkToolPtr tp, Int2 item)
{                           
  ValNodePtr vnp, link_list;
  Int4       cds_num = 1;
  SeqFeatPtr cds = NULL;

  if (item == 0 || tp == NULL) return NULL;

  vnp = tp->cds_list;
  while (vnp != NULL && cds_num < item) {
    if (tp->selected_cds == NULL) {
      /* we're showing the report view */
      link_list = GetLinkedFeatureList (vnp->data.ptrvalue);
      if (link_list != NULL && link_list->next != NULL) {
        cds_num += ValNodeLen (link_list->next);
      }
      link_list = ValNodeFree (link_list);
    }
    if (cds_num < item) {
      vnp = vnp->next;
      cds_num++;
    }
  }
  if (vnp != NULL) {
    cds = vnp->data.ptrvalue;
  }
  return cds;
}


static Boolean GrayLinkedCDS (DoC doc, Int2 item, Int2 row, Int2 col)
{
  CDSmRNALinkToolPtr tp;
  SeqFeatPtr         cds;
  ValNodePtr         link_list;
  Boolean            rval = FALSE;
  
  tp = (CDSmRNALinkToolPtr) GetObjectExtra (doc);
  if (tp == NULL) return FALSE;
  
  cds = CDSFromLinkDocItem (tp, item);
  if (cds != NULL && cds != tp->selected_cds) 
  {
    link_list = GetLinkedFeatureList (cds);
    if (link_list != NULL) {
      rval = TRUE;
    }
    link_list = ValNodeFree (link_list);
  } 

  return rval;
}

static Boolean GrayLinkedmRNA (DoC doc, Int2 item, Int2 row, Int2 col)
{
  CDSmRNALinkToolPtr tp;
  ValNodePtr         vnp, link_list;
  Int4               num;
  Boolean            rval = FALSE, is_checked = FALSE;
  
  tp = (CDSmRNALinkToolPtr) GetObjectExtra (doc);
  if (tp == NULL) return FALSE;
  
  for (vnp = tp->mrna_list, num = 1; vnp != NULL && num != item; vnp = vnp->next, num++)
  {
  }
  if (vnp != NULL && num == item) 
  {
    if (vnp->choice > 0) {
      is_checked = TRUE;
    }
    link_list = GetLinkedFeatureList (vnp->data.ptrvalue);
    if (link_list != NULL) {
      if (is_checked) {
        if (tp->selected_cds != NULL) {
          vnp = link_list;
          while (vnp != NULL && !rval) {
            if (vnp->data.ptrvalue != tp->selected_cds) {
              rval = TRUE;
            }
            vnp = vnp->next;
          }
        }
      } else {
        rval = TRUE;
      }
    }
    link_list = ValNodeFree (link_list);
  } 

  return rval;
}

static void DrawmRNA (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  CDSmRNALinkToolPtr tp;
  RecT         rct;
  RecT         doc_rect;
  Int4         mrna_num;
  ValNodePtr   vnp;

  tp = (CDSmRNALinkToolPtr) GetObjectExtra (d);
  
  if (tp == NULL || tp->selected_cds == NULL || r == NULL 
      || item < 1 
      || firstLine != 0)
  {
    return;
  }

  for (vnp = tp->mrna_list, mrna_num = 1; vnp != NULL && mrna_num != item; vnp = vnp->next, mrna_num++) {};
  if (vnp == NULL || mrna_num != item)
  {
    /* don't draw box for empty mRNA field */
    return;
  }
  rct = *r;
  rct.left ++;
  rct.right = rct.left + tp->lineheight;
  rct.bottom = rct.top + (rct.right - rct.left);
  
  /* make sure we don't draw a box where we aren't drawing text */
  ObjectRect (tp->right_doc, &doc_rect);
  InsetRect (&doc_rect, 4, 4);
  if (rct.bottom > doc_rect.bottom)
  {
    return;
  }
  
  FrameRect (&rct);
  
  if (vnp->choice > 0) {
    MoveTo (rct.left, rct.top);
    LineTo (rct.right - 1, rct.bottom - 1);
    MoveTo (rct.left, rct.bottom - 1);
    LineTo (rct.right - 1, rct.top);
  }
}

static void DrawCDSDivide (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  CDSmRNALinkToolPtr tp;
  SeqFeatPtr   this_cds, prev_cds;
  ValNodePtr   this_link, prev_link;

  tp = (CDSmRNALinkToolPtr) GetObjectExtra (d);
  
  if (tp == NULL || r == NULL 
      || item < 2 
      || firstLine != 0)
  {
    return;
  }

  this_cds = CDSFromLinkDocItem (tp, item);
  prev_cds = CDSFromLinkDocItem (tp, item - 1);

  this_link = GetLinkedFeatureList (this_cds);
  if (this_link == NULL) return;
  prev_link = GetLinkedFeatureList (prev_cds);
  if (prev_link == NULL) {
    /* draw line */
    MoveTo (r->left, r->top);
    LineTo (r->right, r->top);
  }
  this_link = ValNodeFree (this_link);
  prev_link = ValNodeFree (prev_link);

}

static void ReleasemRNA (DoC d, PoinT pt)

{
  Int2            col;
  CDSmRNALinkToolPtr tp;
  Int2            item;
  RecT            rct;
  Int2            row;
  Int4            mrna_num;
  ValNodePtr      vnp, vnp_l;
  SeqEntryPtr     sep;
  SeqFeatPtr      mrna;
  ValNodePtr      link_list;
  SeqFeatPtr      linked_to_other = NULL, cds_link_to_other = NULL;
  Boolean         linked_to_selected = FALSE;
  SeqMgrFeatContext fcontext;

  tp = (CDSmRNALinkToolPtr) GetObjectExtra (d);
  if (tp == NULL || tp->selected_cds == NULL) return;

  MapDocPoint (d, pt, &item, &row, &col, &rct);
  rct.left += 1;
  rct.right = rct.left + tp->lineheight;
  rct.bottom = rct.top + (rct.right - rct.left);
  if (row != 1 || col != 1 || !PtInRect (pt, &rct)) {
    /* didn't click on a box */
    return;
  }

  for (vnp = tp->mrna_list, mrna_num = 1;
       vnp != NULL && mrna_num != item;
       vnp = vnp->next, mrna_num++)
  {
  }
  if (vnp == NULL || mrna_num != item) {
    /* clicked beyond end of list */
    return;
  }

  /* check to see if this mRNA is already linked to this CDS - 
   * always want to let the user uncheck and break the link
   */
  mrna = vnp->data.ptrvalue;

  link_list = GetLinkedFeatureList (tp->selected_cds);
  for (vnp_l = link_list; 
       vnp_l != NULL && !linked_to_selected; 
       vnp_l = vnp_l->next) {
    if (vnp_l->data.ptrvalue == mrna) {
      linked_to_selected = TRUE;
    } else {
      cds_link_to_other = vnp_l->data.ptrvalue;
    }
  }
  link_list = ValNodeFree (link_list);

  if (!linked_to_selected) {
    if (cds_link_to_other != NULL) {
      Message (MSG_ERROR, "This CDS is already linked to an mRNA - you must unlink the other mRNA before you can link to this mRNA.");
      return;
    } else {
      link_list = GetLinkedFeatureList (mrna);
      for (vnp_l = link_list;
           vnp_l != NULL && linked_to_other == NULL;
           vnp_l = vnp_l->next) {
        if (vnp_l->data.ptrvalue != tp->selected_cds) {
          linked_to_other = vnp_l->data.ptrvalue;
        }
      }
      link_list = ValNodeFree (link_list);
      if (linked_to_other != NULL) {
        linked_to_other = SeqMgrGetDesiredFeature (linked_to_other->idx.entityID, 
                                                   NULL, 
                                                   linked_to_other->idx.itemID, 
                                                   0,
                                                   linked_to_other, 
                                                   &fcontext);
        Message (MSG_ERROR, "This mRNA is linked to another CDS (%s).  You must unlink it before you can link it to this CDS.", fcontext.label);
        return;
      }
    }
  }

  sep = GetTopSeqEntryForEntityID (tp->selected_cds->idx.entityID);
  AssignFeatureIDs (sep);
  if (vnp->choice > 0) {
    vnp->choice = 0;
    /* remove this link */
    RemoveFeatureLink (tp->selected_cds, mrna);
    RemoveFeatureLink (mrna, tp->selected_cds);
  } else {
    vnp->choice = OBJ_SEQFEAT;
    /* add this link */
    LinkTwoFeatures (tp->selected_cds, mrna);
    LinkTwoFeatures (mrna, tp->selected_cds);
  }
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  InsetRect (&rct, -1, -1);
  InvalRect (&rct);
  ObjMgrSetDirtyFlag (tp->selected_cds->idx.entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tp->selected_cds->idx.entityID, 0, 0);
  Update ();
}


static void ReleaseCDS (DoC d, PoinT pt)

{
  Int2            col;
  CDSmRNALinkToolPtr tp;
  Int2            item;
  RecT            rct;
  Int2            row;
  ValNodePtr      vnp, vnp_m, linked_mrnas;
  Boolean         found;
  Boolean         was_selected = FALSE;
  SeqFeatPtr      cds;

  tp = (CDSmRNALinkToolPtr) GetObjectExtra (d);
  if (tp == NULL) return;

  if (tp->selected_cds != NULL) {
    was_selected = TRUE;
  }

  MapDocPoint (d, pt, &item, &row, &col, &rct);
  cds = CDSFromLinkDocItem (tp, item);

  if (cds == NULL || cds == tp->selected_cds) {
    tp->selected_cds = NULL;
  } else {
    tp->selected_cds = cds;
  }

  if (tp->selected_cds != NULL) {
    /* need to label mRNAs - those already linked to CDS should be 
     * listed at the end, those not linked should be unchecked.
     */
    linked_mrnas = GetLinkedFeatureList (tp->selected_cds);

    for (vnp = tp->mrna_list; vnp != NULL; vnp = vnp->next) {
      for (vnp_m = linked_mrnas, found = FALSE; vnp_m != NULL && !found; vnp_m = vnp_m->next) {
        if (vnp_m->data.ptrvalue == vnp->data.ptrvalue) {
          found = TRUE;
        }
      }
      if (found) {
        vnp->choice = OBJ_SEQFEAT;
      } else {
        vnp->choice = 0;
      }
    }
    linked_mrnas = ValNodeFree (linked_mrnas);

    tp->mrna_list = ValNodeSort (tp->mrna_list, SortVnpBymRNAPos);
      
  }

  if (tp->selected_cds == NULL || !was_selected) {
    PopulateLinkTool (tp);
  } else {
    ResetClip();
    ObjectRect (tp->left_doc, &rct);
    InvalRect (&rct);

    PopulateLinkRightDoc (tp);
    Show (tp->edit_grp);
    Hide (tp->report_grp);
    Update ();
  }
}

static void CollectCDSCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || data == NULL) return;
  ValNodeAddPointer (data, OBJ_SEQFEAT, sfp);
}

static void CollectmRNACallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA || data == NULL) return;
  ValNodeAddPointer (data, OBJ_SEQFEAT, sfp);
}

static void FormatLinkColumn (Nlm_ColPtr col)
{
  if (col == NULL) return;

  col->pixWidth    = 0;
  col->pixInset    = 0;
  col->charWidth   = 10;
  col->charInset   = 0;
  col->font        = NULL;
  col->just        = 'l';
  col->wrap        = TRUE;
  col->bar         = FALSE;
  col->underline   = FALSE;
  col->left        = FALSE;
  col->last        = FALSE;
}

extern void CDSmRNALinkTool (IteM i)
{
  CDSmRNALinkToolPtr tp;
  GrouP              h, k, button_grp;
  ButtoN             b;
  WindoW             w;
  BaseFormPtr        bfp;
  OMUserDataPtr      omudp;
  SeqEntryPtr        sep;
  Int4               col_num;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  tp = (CDSmRNALinkToolPtr) MemNew (sizeof (CDSmRNALinkToolData));
  if (tp == NULL)
    return;

  w = FixedWindow (-50, -33, -10, -10, "CDS mRNA Link Tool",
		   StdCloseWindowProc);
  SetObjectExtra (w, tp, CleanupCDSmRNALinkTool);
  tp->form = (ForM) w;
  tp->input_entityID = bfp->input_entityID;
  tp->bfp = bfp;

  /* register to receive update messages */
  tp->userkey = OMGetNextUserKey ();
  tp->procid = 0;
  tp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (tp->input_entityID, tp->procid, tp->proctype, tp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) tp;
    omudp->messagefunc = CDSmRNALinkToolFormMsgFunc;
  }

  tp->cdsParFmt.openSpace    = FALSE;
  tp->cdsParFmt.keepWithNext = FALSE;
  tp->cdsParFmt.keepTogether = FALSE;
  tp->cdsParFmt.newPage      = FALSE;
  tp->cdsParFmt.tabStops     = FALSE;
  tp->cdsParFmt.minLines     = 0;
  tp->cdsParFmt.minHeight    = 0;

  tp->mrnaParFmt.openSpace    = FALSE;
  tp->mrnaParFmt.keepWithNext = FALSE;
  tp->mrnaParFmt.keepTogether = FALSE;
  tp->mrnaParFmt.newPage      = FALSE;
  tp->mrnaParFmt.tabStops     = FALSE;
  tp->mrnaParFmt.minLines     = 0;
  tp->mrnaParFmt.minHeight    = 0;

  tp->reportParFmt.openSpace    = FALSE;
  tp->reportParFmt.keepWithNext = FALSE;
  tp->reportParFmt.keepTogether = FALSE;
  tp->reportParFmt.newPage      = FALSE;
  tp->reportParFmt.tabStops     = FALSE;
  tp->reportParFmt.minLines     = 0;
  tp->reportParFmt.minHeight    = 0;

  for (col_num = 0; col_num < 2; col_num++) {
    FormatLinkColumn (tp->cdsColFmt + col_num);
    FormatLinkColumn (tp->mrnaColFmt + col_num);
    FormatLinkColumn (tp->reportColFmt + col_num);
    FormatLinkColumn (tp->reportColFmt + col_num + 2);
  }
  FormatLinkColumn (tp->mrnaColFmt + 2);
  tp->cdsColFmt[1].last       = TRUE;
  tp->mrnaColFmt[2].last      = TRUE;
  tp->reportColFmt[3].last    = TRUE;

  SelectFont (programFont);
  tp->lineheight = LineHeight ();

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  k = HiddenGroup (w, 0, 0, NULL);
  tp->edit_grp = HiddenGroup (k, 2, 0, NULL);
  StaticPrompt (tp->edit_grp, "CDS", 0, dialogTextHeight, programFont, 'c');
  StaticPrompt (tp->edit_grp, "mRNA", 0, dialogTextHeight, programFont, 'c');
  tp->left_doc = DocumentPanel (tp->edit_grp, stdCharWidth * 30 + 5, tp->lineheight * 40);
  SetObjectExtra (tp->left_doc, tp, NULL);
  SetDocProcs (tp->left_doc, NULL, NULL, ReleaseCDS, NULL); 
  SetDocShade (tp->left_doc, NULL, GrayLinkedCDS, HighlightSelectedCDS, NULL);

  tp->right_doc = DocumentPanel (tp->edit_grp, stdCharWidth * 30 + 5 + tp->lineheight, tp->lineheight * 40);
  SetObjectExtra (tp->right_doc, tp, NULL);
  SetDocProcs (tp->right_doc, NULL, NULL, ReleasemRNA, NULL); 
  SetDocShade (tp->right_doc, DrawmRNA, GrayLinkedmRNA, NULL, NULL);

  tp->report_grp = HiddenGroup (k, 0, 2, NULL);
  tp->report_title_doc = DocumentPanel (tp->report_grp, stdCharWidth * 60 + 10, tp->lineheight);
  tp->report_doc = DocumentPanel (tp->report_grp, stdCharWidth * 60 + 10, tp->lineheight * 40);
  SetObjectExtra (tp->report_doc, tp, NULL);
  SetDocProcs (tp->report_doc, NULL, NULL, ReleaseCDS, NULL); 
  SetDocShade (tp->report_doc, DrawCDSDivide, GrayLinkedCDS, NULL, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) tp->edit_grp, (HANDLE) tp->report_grp, NULL);

  button_grp = HiddenGroup (w, 2, 0, NULL);
  b = PushButton (button_grp, "Refresh", RefreshLinkBtn);
  SetObjectExtra (b, tp, NULL);
  PushButton (button_grp, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) button_grp, NULL);

  sep = GetTopSeqEntryForEntityID (tp->input_entityID);
  /* collect CDSs */
  VisitFeaturesInSep (sep, &(tp->cds_list), CollectCDSCallback);
  /* collect mRNAs */
  VisitFeaturesInSep (sep, &(tp->mrna_list), CollectmRNACallback);

  PopulateLinkTool (tp);
  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


typedef struct setqual {
  SourceQualDescPtr qual;
  CharPtr           val;
} SetQualData, PNTR SetQualPtr;

static void SetQualOnSourceDescWhenSourceFeatPresentCallback (BioseqPtr bsp, Pointer userdata)
{
  SetQualPtr        sqp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  OrgModPtr         mod;
  SubSourcePtr      ssp;

  if (bsp == NULL || userdata == NULL) return;
  sqp = (SetQualPtr) userdata;
  if (sqp->qual == NULL) return;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  if (sfp == NULL) return;

  if (sqp->qual->isOrgMod)
  {
    if (biop->org == NULL) 
    {
      biop->org = OrgRefNew();
    }
    if (biop->org->orgname == NULL) 
    {
      biop->org->orgname = OrgNameNew();
    }
    mod = biop->org->orgname->mod;
    while (mod != NULL && mod->subtype != sqp->qual->subtype)
    {
      mod = mod->next;
    }
    if (mod == NULL)
    {
      mod = OrgModNew ();
      mod->subtype = sqp->qual->subtype;
      mod->next = biop->org->orgname->mod;
      biop->org->orgname->mod = mod;
    }
    if (IsNonTextModifier (sqp->qual->name))
    {
      if (mod->subname == NULL)
      {
        mod->subname = StringSave ("");
      }
    }
    else
    {
      mod->subname = MemFree (mod->subname);
      mod->subname = StringSave (sqp->val);
    }
  }
  else
  {
    ssp = biop->subtype;
    while (ssp != NULL && ssp->subtype != sqp->qual->subtype)
    {
      ssp = ssp->next;
    }
    if (ssp == NULL)
    {
      ssp = SubSourceNew();
      ssp->subtype = sqp->qual->subtype;
      ssp->next = biop->subtype; 
      biop->subtype = ssp;
    }
    if (IsNonTextModifier (sqp->qual->name))
    {
      if (ssp->name == NULL)
      {
        ssp->name = StringSave ("");
      }
    }
    else
    {
      ssp->name = MemFree (ssp->name);
      ssp->name = StringSave (sqp->val);
    }
  } 

}

extern void SetTransgenicOnSourceDescWhenSourceFeatPresent (IteM i)
{
  BaseFormPtr         bfp;
  SeqEntryPtr         sep;
  SetQualData         sqd;
  SourceQualDescData  sqdd;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sqdd.isOrgMod = FALSE;
  sqdd.name = "transgenic";
  sqdd.subtype = SUBSRC_transgenic;

  sqd.qual = &sqdd;
  sqd.val = NULL;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitBioseqsInSep (sep, &sqd, SetQualOnSourceDescWhenSourceFeatPresentCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();  
}


static void SetFocusOnSourceDescWhenSourceFeatPresentCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;

  if (bsp == NULL) return;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  if (sfp == NULL) return;

  biop->is_focus = TRUE;
}

extern void SetFocusOnSourceDescWhenSourceFeatPresent (IteM i)
{
  BaseFormPtr         bfp;
  SeqEntryPtr         sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitBioseqsInSep (sep, NULL, SetFocusOnSourceDescWhenSourceFeatPresentCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();  
}


typedef struct cdstrnaoverlap {
  SeqFeatPtr cds;
  ValNodePtr trna_list;
  Int4       overlap_len;
  CharPtr    end_str;
} CDStRNAOverlapData, PNTR CDStRNAOverlapPtr;

static void TrimCDSFortRNAOverlap (ValNodePtr item_list)
{
   
}

static Int4 GetOverlapLen (SeqLocPtr slp1, SeqLocPtr slp2)
{
  Int4  start1, stop1, start2, stop2, tmp, a, b, len = 0;

  if (slp1 == NULL || slp2 == NULL) return 0;

  start1 = SeqLocStart (slp1);
  stop1 = SeqLocStop (slp1);
  if (start1 > stop1) {
    tmp = start1;
    start1 = stop1;
    stop1 = tmp;
  }
  start2 = SeqLocStart (slp2);
  stop2 = SeqLocStop (slp2);
  if (start2 > stop2) {
    tmp = start2;
    start2 = stop2;
    stop2 = tmp;
  }
  
  a = MAX(start1, start2);
  b = MIN (stop1, stop2);
  len = b - a + 1;
  return len; 
}


static CDStRNAOverlapPtr CDStRNAOverlapFree (CDStRNAOverlapPtr p)
{
  if (p != NULL) {
    p->trna_list = ValNodeFree (p->trna_list);
    p->end_str = MemFree (p->end_str);
    p = MemFree (p);
  }
  return p;
}


static ValNodePtr FreeCDStRNAOverlapList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = CDStRNAOverlapFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}

static Boolean IsBadEndStr (CharPtr end_str)
{
  if (end_str == NULL) return TRUE;
  if (end_str[0] != 'T') return TRUE;
  if (end_str[1] == 0) return FALSE;
  if (end_str[1] == 'A' && end_str[2] == 0) return FALSE;
  return TRUE;
}


static ClickableItemPtr ClickableItemFromCDStRNAOverlapList (CDStRNAOverlapPtr p)
{
  ClickableItemPtr cip;
  CharPtr          fmt1 = "CDS is overlapped by more than one tRNA, longest overlap is %d";
  CharPtr          fmt2 = "Overlap is too long: %d";

  if (p == NULL || p->cds == NULL || p->trna_list == NULL) return NULL;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, p->cds);

  /* potential problems - too many tRNAs, too much overlap, wrong string */
  /* only report wrong string if only one tRNA and overlap correct length */

  if (p->overlap_len > 2) {
    if (p->trna_list->next != NULL) {
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt1) + 15));
      sprintf (cip->description, fmt1, p->overlap_len);
    } else {
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt2) + 15));
      sprintf (cip->description, fmt2, p->overlap_len);
    }
  } else if (p->trna_list->next != NULL) {
    cip->description = StringSave ("CDS is overlapped by more than one tRNA!");
  } else if (IsBadEndStr(p->end_str)) {
    cip->description = StringSave ("Base pairs in partial codon before trim are neither 'T' nor 'TA'");
  } else {
    cip->description = StringSave ("Expected overlap");
  }

  return cip;
}


static CharPtr GetEndStrForOverlap (SeqLocPtr slp, Int4 overlap_len)
{
  Uint1 strand;
  Int4  end, start, stop;
  SeqPortPtr spp;
  CharPtr    end_str;

  if (slp == NULL || overlap_len < 0 || overlap_len > 2) {
    return NULL;
  }
  strand = SeqLocStrand (slp);
  if (strand == Seq_strand_minus) {
    end = SeqLocStart (slp);
    start = end + overlap_len;
    stop = end + 2;
  } else {
    end = SeqLocStop (slp);
    start = end - 2;
    stop = end - overlap_len + 1;
  }
  spp = SeqPortNew (BioseqFindFromSeqLoc (slp), start, stop, strand, Seq_code_iupacna);
  end_str = (CharPtr) MemNew (sizeof (Char) + (3 - overlap_len + 1));
  SeqPortRead (spp, (Uint1Ptr)end_str, 3 - overlap_len);
  spp = SeqPortFree (spp);
  return end_str;
}

static ValNodePtr CDStRNAOverlapListFromItemList (ValNodePtr item_list)
{
  ValNodePtr new_item_list = NULL, vnp, vnp_trna;
  ValNodePtr        good_list = NULL, bad_list = NULL;
  SeqFeatPtr        sfp, trna;
  CDStRNAOverlapPtr p = NULL;
  Int4              len;

  if (item_list == NULL) return NULL;

  /* collect overlap groups */
  vnp = item_list;
  while (vnp != NULL) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp->idx.subtype == FEATDEF_CDS) {
        p = (CDStRNAOverlapPtr) MemNew (sizeof (CDStRNAOverlapData));
        p->cds = sfp;
        ValNodeAddPointer (&new_item_list, 0, p);
      } else if (sfp->idx.subtype == FEATDEF_tRNA && p != NULL) {
        ValNodeAddPointer (&(p->trna_list), OBJ_SEQFEAT, sfp);
      }
    }
    vnp = vnp->next;
  }

  for (vnp = new_item_list; vnp != NULL; vnp = vnp->next) {
    p = (CDStRNAOverlapPtr) vnp->data.ptrvalue;
    if (p->cds == NULL || p->trna_list == NULL) continue;
    for (vnp_trna = p->trna_list; vnp_trna != NULL; vnp_trna = vnp_trna->next) {
      if (vnp_trna->choice != OBJ_SEQFEAT || vnp_trna->data.ptrvalue == NULL) continue;
      trna = vnp_trna->data.ptrvalue;
      len = GetOverlapLen (p->cds->location, trna->location);
      if (p->overlap_len < len) {
        p->overlap_len = len;
      }
    }
    if (p->overlap_len > 0) {
      p->end_str = GetEndStrForOverlap (p->cds->location, p->overlap_len);
    }

    if (p->overlap_len == 0 || p->overlap_len > 2)
    {
      ValNodeAddPointer (&bad_list, 0, p);
    }
    else if (IsBadEndStr(p->end_str))
    {
      ValNodeAddPointer (&bad_list, 0, p);
    }
    else 
    {
      ValNodeAddPointer (&good_list, 0, p);
    }
    vnp->data.ptrvalue = NULL;
  }
  new_item_list = ValNodeFree (new_item_list);
  new_item_list = bad_list;
  ValNodeLink (&new_item_list, good_list);
  return new_item_list;
}


static ValNodePtr ClickableListFromCDStRNAOverlapList (ValNodePtr overlap_list)
{
  ClickableItemPtr cip;
  ValNodePtr       vnp, clickable_list = NULL;

  /* create clickable list of coding regions to trim */
  for (vnp = overlap_list; vnp != NULL; vnp = vnp->next)
  {
    cip = ClickableItemFromCDStRNAOverlapList (vnp->data.ptrvalue);
    if (cip != NULL) {
      ValNodeAddPointer (&clickable_list, 0, cip);
    }
  }
  return clickable_list;
}


static void FixCDStRNAOverlaps (ValNodePtr overlap_list)
{
  CDStRNAOverlapPtr p;
  ValNodePtr        vnp;
  SeqFeatPtr        gene;

  for (vnp = overlap_list; vnp != NULL; vnp = vnp->next)
  {
    p = (CDStRNAOverlapPtr) vnp->data.ptrvalue;
    if (p == NULL || p->cds == NULL || p->overlap_len < 1)
    {
      continue;
    }
    p->cds->location = TruncateLocation (p->cds->location, SeqLocLen (p->cds->location) - p->overlap_len);
    gene = GetGeneForFeature (p->cds);
    if (gene != NULL) 
    {
      gene->location = TruncateLocation (gene->location, SeqLocLen (gene->location) - p->overlap_len);
    }

    AddTranslExcept (p->cds, "TAA stop codon is completed by the addition of 3' A residues to the mRNA", FALSE, FALSE, FALSE);
  }
}  


static Uint2 GetEntityIDFromItem (ValNodePtr vnp)
{
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr   bsp;
  BioseqSetPtr bssp;  
  Uint2        entityID = 0;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return 0;
  switch (vnp->choice) {
    case OBJ_SEQFEAT:
      sfp = vnp->data.ptrvalue;
      entityID = sfp->idx.entityID;
      break;
    case OBJ_SEQDESC:
      sdp = vnp->data.ptrvalue;
      entityID = ((ObjValNodePtr) sdp)->idx.entityID;
      break;
    case OBJ_BIOSEQ:
      bsp = vnp->data.ptrvalue;
      entityID = bsp->idx.entityID;
      break;
    case OBJ_BIOSEQSET:
      bssp = vnp->data.ptrvalue;
      entityID = bssp->idx.entityID;
      break;
  }
  return entityID;
}


static Uint2 GetEntityIDFromClickableItem (ClickableItemPtr cip)
{
  ValNodePtr vnp;
  Uint2      entityID = 0;

  if (cip == NULL || cip->item_list == NULL) return 0;

  for (vnp = cip->item_list; vnp != NULL && entityID == 0; vnp = vnp->next)
  {
    entityID = GetEntityIDFromItem (vnp);
  }
  return entityID;
}
  

static void SendUpdatesForClickableList (ValNodePtr clickable_list)
{
  Uint2Ptr   entityIDtable;
  Uint2      this_entityID;
  ValNodePtr vnp;
  Int4       num, i;

  num = ValNodeLen (clickable_list);
  if (num == 0) return;

  entityIDtable = (Uint2Ptr) MemNew (sizeof (Uint2) * num);
  MemSet (entityIDtable, 0, sizeof (Uint2) * num);
  for (vnp = clickable_list; vnp != NULL; vnp = vnp->next) {
    this_entityID = GetEntityIDFromClickableItem (vnp->data.ptrvalue);
    for (i = 0; i < num && entityIDtable[i] != 0 && entityIDtable[i] != this_entityID; i++) {}
    if (i < num && entityIDtable[i] == 0) {
      entityIDtable[i] = this_entityID;
    }
  }
  for (i = 0; i < num && entityIDtable[i] != 0; i++) {
    ObjMgrSetDirtyFlag (entityIDtable[i], TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, entityIDtable[i], 0, 0);
  }
  Update();
  entityIDtable = MemFree (entityIDtable);
}


static ValNodePtr DivideClickableItemListByEntityID (ValNodePtr clickable_list)
{
  ValNodePtr vnp, vnp2, list_list = NULL, this_list;
  Uint2            entityID;

  for (vnp = clickable_list; vnp != NULL; vnp = vnp->next)
  {
    entityID = GetEntityIDFromClickableItem (vnp->data.ptrvalue);
    for (vnp2 = list_list; vnp2 != NULL && vnp2->choice != entityID; vnp2 = vnp2->next)
    {};
    if (vnp2 == NULL) 
    {    
      vnp2 = ValNodeAddPointer (&list_list, (Uint1) entityID, NULL);
    }
    this_list = vnp2->data.ptrvalue;
    ValNodeAddPointer (&this_list, vnp->choice, vnp->data.ptrvalue);
    vnp->data.ptrvalue = NULL;
    vnp2->data.ptrvalue = this_list;
  }
  return list_list;
}


typedef struct cdstrnatool {
  FORM_MESSAGE_BLOCK
  DialoG clickable_list_dlg;
  ValNodePtr clickable_list;
  ValNodePtr overlap_list;
} CDStRNAToolData, PNTR CDStRNAToolPtr;

static void CleanupCDStRNAToolForm (GraphiC g, VoidPtr data)

{
  CDStRNAToolPtr dlg;

  dlg = (CDStRNAToolPtr) data;
  if (dlg != NULL) {
    dlg->clickable_list = FreeClickableList (dlg->clickable_list);
    dlg->overlap_list = FreeCDStRNAOverlapList (dlg->overlap_list);
    ObjMgrFreeUserData (dlg->input_entityID, dlg->procid, dlg->proctype, dlg->userkey);
  }
  StdCleanupFormProc (g, data);
}

static ValNodePtr GetSelectedClickableItems (ValNodePtr list)
{
  ValNodePtr       vnp, selected = NULL;
  ClickableItemPtr cip;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip == NULL) continue;
    if (cip->chosen) {
      ValNodeAddPointer (&selected, 0, cip);
    } else {
      ValNodeLink (&selected, GetSelectedClickableItems (cip->subcategories));
    }
  }
  return selected;
}


static void TrimSelectedCDS (ButtoN b)
{
  CDStRNAToolPtr    dlg;
  ValNodePtr        tmp_list = NULL, overlap_list = NULL, vnp_o, vnp_c;
  ValNodePtr        prev_o = NULL, next_o;
  ClickableItemPtr  cip;
  CDStRNAOverlapPtr p;
  SeqFeatPtr        cds;

  dlg = (CDStRNAToolPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  tmp_list = GetSelectedClickableItems (dlg->clickable_list);

  if (tmp_list == NULL) {
    if (ANS_CANCEL == Message (MSG_OKC, "No features selected!  Trim all?")) {
      return;
    } else {
      tmp_list = dlg->clickable_list;
    }
  }

  vnp_c = tmp_list;
  vnp_o = dlg->overlap_list;
  while (vnp_c != NULL && vnp_o != NULL) {
    cip = (ClickableItemPtr) vnp_c->data.ptrvalue;
    if (cip == NULL || cip->item_list == NULL || cip->item_list->choice != OBJ_SEQFEAT) {
      vnp_c = vnp_c->next;
    } else {
      cds = cip->item_list->data.ptrvalue;
      p = (CDStRNAOverlapPtr) vnp_o->data.ptrvalue;
      if (p == NULL || p->cds != cds) {
        prev_o = vnp_o;
        vnp_o = vnp_o->next;
      } else {
        ValNodeAddPointer (&overlap_list, 0, p);
        /* better not be repeats in overlap list */
        next_o = vnp_o->next;
        if (prev_o == NULL) {
          dlg->overlap_list = next_o;
        } else {
          prev_o->next = next_o;
        }
        vnp_o->next = NULL;
        vnp_o = ValNodeFree (vnp_o);
        vnp_o = next_o;
        /* move on to next item in clickable list */
        vnp_c = vnp_c->next;
      }
    }
  }
  
  FixCDStRNAOverlaps (overlap_list);
  overlap_list = FreeCDStRNAOverlapList (overlap_list);
  SendUpdatesForClickableList (tmp_list);

  if (tmp_list != dlg->clickable_list) {
    tmp_list = ValNodeFree (tmp_list);
  }

  /* close the window if we're done */
  if (dlg->overlap_list == NULL) {
    Remove (dlg->form);
  } else {
    /* update list of overlaps (remove the trimmed ones) */
    PointerToDialog (dlg->clickable_list_dlg, NULL);
    dlg->clickable_list = FreeClickableList (dlg->clickable_list);
    dlg->clickable_list = ClickableListFromCDStRNAOverlapList (dlg->overlap_list);
    PointerToDialog (dlg->clickable_list_dlg, dlg->clickable_list);
  }
}


static void SelectAllCDS (ButtoN b)
{
  CDStRNAToolPtr    dlg;

  dlg = (CDStRNAToolPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  ChooseCategories (dlg->clickable_list, TRUE);
  PointerToDialog (dlg->clickable_list_dlg, dlg->clickable_list);

}


static void EditCDStRNAOverlap (ValNodePtr item_list)
{
  ValNodePtr overlap_list;
  CDStRNAToolPtr dlg;
  WindoW      w;
  GrouP       h, c;
  ButtoN      b;
  OMUserDataPtr            omudp;

  overlap_list = CDStRNAOverlapListFromItemList (item_list);
  if (overlap_list == NULL) {
    Message (MSG_ERROR, "No overlaps to correct!");
    return;
  }

  /* create window to display and correct overlaps */
  dlg = (CDStRNAToolPtr) MemNew (sizeof (CDStRNAToolData));
  dlg->clickable_list = ClickableListFromCDStRNAOverlapList (overlap_list);
  dlg->overlap_list = overlap_list;

  w = FixedWindow (-50, -33, -10, -10, "CDS tRNA Overlaps", StdCloseWindowProc);
  SetObjectExtra (w, dlg, CleanupCDStRNAToolForm);
  dlg->form = (ForM) w;
    
  /* register to receive update messages */
  dlg->userkey = OMGetNextUserKey ();
  dlg->procid = 0;
  dlg->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (0, dlg->procid, dlg->proctype, dlg->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) dlg;
    omudp->messagefunc = ClickableListFormMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  dlg->clickable_list_dlg = CreateClickableListDialog (h, "CDS-tRNA Overlaps", "Features",
                                                        ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                        GetDiscrepancyItemText);
  PointerToDialog (dlg->clickable_list_dlg, dlg->clickable_list);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Trim Selected", TrimSelectedCDS);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (c, "Select All", SelectAllCDS);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->clickable_list_dlg, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


extern GetSamplePtr 
GetSampleForItemList 
(ValNodePtr               item_list,
 ValNodePtr               requested_field,
 GetFeatureFieldString    fieldstring_func,
 GetDescriptorFieldString descrstring_func)
{
  ValNodePtr   vnp;
  GetSamplePtr gsp;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  CharPtr      str = NULL;

  if (item_list == NULL || requested_field == NULL) return NULL;

  gsp = GetSampleNew();

  for (vnp = item_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT && fieldstring_func != NULL)
    {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      str = fieldstring_func (sfp, requested_field, NULL);
    }
    else if (vnp->choice == OBJ_SEQDESC && descrstring_func != NULL)
    {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      str = descrstring_func (sdp, requested_field, NULL);
    }
    if (!StringHasNoText (str))
    {
      gsp->num_found ++;
      if (gsp->sample_text == NULL)
      {
        gsp->sample_text = str;
      }
      else
      {
        if (StringCmp (str, gsp->sample_text) != 0)
        {
          gsp->all_same = FALSE;
        }
      }
    }
    str = MemFree (str);
  }

  if (gsp->num_found == 0)
  {
    gsp = GetSampleFree (gsp);
  }
  return gsp;
}


typedef struct applytagtocdsinsrcfeat {
  FORM_MESSAGE_BLOCK
  TexT       note_text;
  ValNodePtr cds_list;
} ApplyTagToCDSInSrcFeatData, PNTR ApplyTagToCDSInSrcFeatPtr;


static void CleanupApplyTagToCDSInSrcFeatForm (GraphiC g, VoidPtr data)

{
  ApplyTagToCDSInSrcFeatPtr dlg;

  dlg = (ApplyTagToCDSInSrcFeatPtr) data;
  if (dlg != NULL) {
    dlg->cds_list = ValNodeFree (dlg->cds_list);
  }
  StdCleanupFormProc (g, data);
}


static void DoApplyTagToCodingRegionsInSourceFeatures (ButtoN b)
{
  ApplyTagToCDSInSrcFeatPtr dlg;
  ValNodePtr                vnp;
  GetSamplePtr              gsp;
  ValNode                   vn;
  SeqFeatPtr                sfp;
  ApplyValueData            avd;

  dlg = (ApplyTagToCDSInSrcFeatPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (TextHasNoText (dlg->note_text)) {
    Message (MSG_ERROR, "No text supplied for note!");
    return;
  }
  vn.choice = 0;
  vn.next = NULL;
  vn.data.intvalue = 2;
  gsp = GetSampleForItemList (dlg->cds_list, &vn, GetCDSFieldString, NULL);
  avd.etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
  gsp = GetSampleFree (gsp);  
  if (avd.etp != NULL && avd.etp->existing_text_choice == eExistingTextChoiceCancel)
  {
    avd.etp = MemFree (avd.etp);
    return;
  }
  avd.new_text = SaveStringFromText (dlg->note_text);
  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;

  WatchCursor();
  Update();
  for (vnp = dlg->cds_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      sfp->comment = HandleApplyValue (sfp->comment, &avd);
    }
  }
  avd.etp = MemFree (avd.etp);
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  ArrowCursor();
  Remove (dlg->form);
  Update ();
}


static void ApplyTagToCodingRegionsInEntityID (Uint2 entityID)
{
  ApplyTagToCDSInSrcFeatPtr dlg;
  ValNodePtr   cds_list;
  WindoW       w;
  GrouP        h, g, c;
  ButtoN       b;
  SeqEntryPtr  sep;

  sep = GetTopSeqEntryForEntityID (entityID);

  cds_list = ListCodingRegionsContainedInSourceFeatures (sep);

  if (cds_list == NULL) {
    Message (MSG_ERROR, "No coding regions found in source features!");
    return;
  }

  /* create window to display and correct overlaps */
  dlg = (ApplyTagToCDSInSrcFeatPtr) MemNew (sizeof (ApplyTagToCDSInSrcFeatData));
  dlg->input_entityID = entityID;
  dlg->cds_list = cds_list;

  w = FixedWindow (-50, -33, -10, -10, "Apply Note to Coding Regions in Source Features", StdCloseWindowProc);
  SetObjectExtra (w, dlg, CleanupApplyTagToCDSInSrcFeatForm);
  dlg->form = (ForM) w;
    
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Note Text:", 0, dialogTextHeight, programFont, 'r');
  dlg->note_text = DialogText (g, "prophage-encoded protein", 0, NULL);
  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Accept", DoApplyTagToCodingRegionsInSourceFeatures);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


extern void ApplyTagToCodingRegionsInSourceFeatures (IteM i)
{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ApplyTagToCodingRegionsInEntityID (bfp->input_entityID);
}

