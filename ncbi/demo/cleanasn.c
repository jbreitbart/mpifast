/*   cleanasn.c
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
* File Name:  cleanasn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   10/19/99
*
* $Revision: 6.29 $
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
#include <objfdef.h>
#include <objsub.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <tofasta.h>
#include <toasn3.h>
#include <subutil.h>
#include <asn2gnbk.h>
#include <pmfapi.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_CLEANASN
#include <accpubseq.h>
#endif

#define CLEANASN_APP_VER "2.2"

CharPtr CLEANASN_APPLICATION = CLEANASN_APP_VER;

typedef struct cleanflags {
  Char          buf [64];
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Int2          type;
  CharPtr       results;
  CharPtr       outfile;
  CharPtr       report;
  CharPtr       ffdiff;
  ModType       ffmode;
  CharPtr       clean;
  CharPtr       modernize;
  CharPtr       link;
  CharPtr       feat;
  CharPtr       desc;
  CharPtr       mods;
  Boolean       taxon;
  Boolean       pub;
  Int4          okay;
  Int4          bsec;
  Int4          ssec;
  Int4          norm;
  Int4          cumokay;
  Int4          cumbsec;
  Int4          cumssec;
  Int4          cumnorm;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  BioseqSet     bss;
  FILE          *logfp;
} CleanFlagData, PNTR CleanFlagPtr;

static void RemoveFeatUser (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  if (sfp->ext != NULL) {
    sfp->ext = UserObjectFree (sfp->ext);
  }
}

static void RemoveFeatDbxref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  DbtagPtr    dbt;
  ValNodePtr  next, vnp;

  if (sfp == NULL) return;
  for (vnp = sfp->dbxref; vnp != NULL; vnp = next) {
    next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    DbtagFree (dbt);
    MemFree (vnp);
  }
  sfp->dbxref = NULL;
}

typedef struct dummysmfedata {
  Int4  max;
  Int4  num_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK CADummySMFEProc (
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

static void RemoveUnnecGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  DummySmfeData        dsd;
  SeqMgrFeatContext    fcontext;
  SeqFeatXrefPtr PNTR  last;
  GeneRefPtr           grp, grpx;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL) return;

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) != 0) return;
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if (StringICmp (grp->locus, grpx->locus) != 0) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) != 0) return;
    }
  }

  MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
  dsd.max = INT4_MAX;
  dsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE,
                                           NULL, 0, LOCATION_SUBSET,
                                           (Pointer) &dsd, CADummySMFEProc);

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

static void MarkTitles (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

static void AddSpTaxnameToList (
  SeqDescrPtr sdp,
  Pointer userdata
)
{
  BioSourcePtr biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL || biop->org == NULL || !IsSpName (biop->org->taxname)) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, biop->org->taxname);
}


static Boolean ShouldExcludeSp (
  SeqEntryPtr sep
)

{
  ValNodePtr name_list = NULL, vnp1, vnp2;
  Boolean    all_diff = TRUE;

  if (sep == NULL) return TRUE;
  VisitDescriptorsInSep (sep, &name_list, AddSpTaxnameToList);

  name_list = ValNodeSort (name_list, SortVnpByString);

  if (name_list != NULL && name_list->next != NULL) {
    for (vnp1 = name_list; vnp1 != NULL && vnp1->next != NULL && all_diff; vnp1 = vnp1->next) {
      for (vnp2 = vnp1->next; vnp2 != NULL && all_diff; vnp2 = vnp2->next) {
        if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == 0) {
          all_diff = FALSE;
        }
      }
    }
  }

  name_list = ValNodeFree (name_list);

  return all_diff;
}

static void DoAutoDef (
  SeqEntryPtr sep,
  Uint2 entityID
)

{
  ValNodePtr                    defline_clauses = NULL;
  DeflineFeatureRequestList     feature_requests;
  Int4                          index;
  ValNodePtr                    modifier_indices = NULL;
  ModifierItemLocalPtr          modList;
  OrganismDescriptionModifiers  odmp;
  SeqEntryPtr                   oldscope;

  if (sep == NULL) return;
  if (entityID < 1) return;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;

  InitFeatureRequests (&feature_requests);

  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  odmp.use_labels = TRUE;
  odmp.max_mods = -99;
  odmp.keep_paren = TRUE;
  odmp.exclude_sp = ShouldExcludeSp (sep);
  odmp.exclude_cf = FALSE;
  odmp.exclude_aff = FALSE;
  odmp.exclude_nr = FALSE;
  odmp.include_country_extra = FALSE;
  odmp.clone_isolate_HIV_rule_num = clone_isolate_HIV_rule_want_both;
  odmp.use_modifiers = FALSE;
  odmp.allow_semicolon_in_modifier = FALSE;


  RemoveNucProtSetTitles (sep);  
  oldscope = SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if (AreFeatureClausesUnique (defline_clauses)) {
    modifier_indices = GetModifierIndicesFromModList (modList);
  } else {
    modifier_indices = FindBestModifiers (sep, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);
  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL) {
    for (index = 0; index < NumDefLineModifiers (); index++) {
      ValNodeFree (modList [index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);

  SeqEntrySetScope (oldscope);
}

static void LookupPubdesc (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid = 0;
  ValNodePtr       vnp;

  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        /* ignore obsolete muids */
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        /* return on real pub */
        return;
        break;
    }
  }

  if (pmid == 0) return;

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
  }

  PubmedEntryFree (pep);
}

static void ModGenes (SeqFeatPtr sfp, Pointer userdata)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (SeqFeatPtr sfp, Pointer userdata)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (BioSourcePtr biop, Pointer userdata)

{
  ModernizePCRPrimers (biop);
}

static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}

typedef struct chgdata {
  Boolean       rubisco;
  Boolean       rbc;
  Boolean       its;
  Boolean       rnaother;
  Boolean       trnanote;
  Boolean       oldbiomol;
  Int4          protdesc;
  Int4          sfpnote;
  Int4          gbsource;
  Int4          cdsconf;
} ChangeData, PNTR ChangeDataPtr;

static Boolean IsRubisco (
  CharPtr name
)

{
  return (StringICmp (name, "rubisco large subunit") == 0 ||
          StringICmp (name, "rubisco small subunit") == 0);
}

static Boolean IsRbc (
  CharPtr name
)

{
  return (StringICmp (name, "RbcL") == 0 ||
          StringICmp (name, "RbcS") == 0);
}

static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  CharPtr        name;
  ProtRefPtr     prp;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (IsRubisco (str)) {
          cdp->rubisco = TRUE;
        }
        if (IsRbc (str)) {
          cdp->rbc = TRUE;
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "ncRNA") == 0 || StringCmp (name, "tmRNA") == 0) {
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
        }
      }
      break;
    default:
      break;
  }
}

static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
}

static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
}

static void DoASNReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  Boolean     bsec = FALSE, ssec = FALSE, norm = FALSE;
  ChangeData  cdbefore, cdafter;
  CharPtr     str1, str2, str3, str4;

  if (sep == NULL || cfp == NULL) return;

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  CheckForChanges (sep, &cdbefore);

  str1 = Se2Str (sep);
  NormalizeDescriptorOrder (sep);
  str2 = Se2Str (sep);
  if (StringCmp (str1, str2) != 0) {
    norm = TRUE;
  }
  BasicSeqEntryCleanup (sep);
  str3 = Se2Str (sep);
  if (StringCmp (str2, str3) != 0) {
    bsec = TRUE;
  }
  SeriousSeqEntryCleanup (sep, NULL, NULL);
  NormalizeDescriptorOrder (sep);
  str4 = Se2Str (sep);
  if (StringCmp (str3, str4) != 0) {
    ssec = TRUE;
  }

  CheckForChanges (sep, &cdafter);

  if (ssec) {
    (cfp->ssec)++;
    (cfp->cumssec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (bsec) {
    (cfp->bsec)++;
    (cfp->cumbsec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (norm) {
    (cfp->norm)++;
    (cfp->cumnorm)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NORM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else {
    (cfp->okay)++;
    (cfp->cumokay)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "OKAY %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.rubisco) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RUB %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rbc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RBC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.its) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ITS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rnaother) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RNA %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.trnanote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TRN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldbiomol) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MOL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PRT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "COM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SRC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CNF %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  MemFree (str1);
  MemFree (str2);
  MemFree (str3);
  MemFree (str4);
}

static void DoGBFFReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
#ifdef OS_UNIX
  BioseqPtr    bsp;
  Char         cmmd [256];
  FILE         *fp;
  SeqEntryPtr  fsep;
  Char         path1 [PATH_MAX];
  Char         path2 [PATH_MAX];
  CharPtr      rep = "reports";
  SeqIdPtr     sip;

  if (sep == NULL || cfp == NULL) return;

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
            rep = "gbreports";
            break;
          case SEQID_EMBL :
            rep = "ebreports";
            break;
          case SEQID_DDBJ :
            rep = "djreports";
            break;
          case SEQID_OTHER :
            rep = "rfreports";
            break;
          default :
            break;
        }
      }
    }
  }

  TmpNam (path1);
  TmpNam (path2);

  fp = FileOpen (path1, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  }
  FileClose (fp);
  SeriousSeqEntryCleanupBulk (sep);
  fp = FileOpen (path2, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  }
  FileClose (fp);

  sprintf (cmmd, "%s -o %s -n %s -d %s", cfp->ffdiff, path1, path2, rep);
  system (cmmd);

  sprintf (cmmd, "rm %s; rm %s", path1, path2);
  system (cmmd);
#endif
}

static void DoModernizeReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  CharPtr  str1, str2, str3, str4;

  str1 = Se2Str (sep);
  VisitFeaturesInSep (sep, NULL, ModGenes);
  str2 = Se2Str (sep);
  if (StringCmp (str1, str2) != 0) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GEN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  VisitFeaturesInSep (sep, NULL, ModRNAs);
  str3 = Se2Str (sep);
  if (StringCmp (str2, str3) != 0) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  VisitBioSourcesInSep (sep, NULL, ModPCRs);
  str4 = Se2Str (sep);
  if (StringCmp (str3, str4) != 0) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  MemFree (str1);
  MemFree (str2);
  MemFree (str3);
  MemFree (str4);
}

static void DoCleanup (
  SeqEntryPtr sep,
  Uint2 entityID,
  CleanFlagPtr cfp
)

{
  BioseqPtr    bsp;
  SeqEntryPtr  fsep;
  SeqIdPtr     sip, siphead;

  if (sep == NULL || cfp == NULL) return;

  StringCpy (cfp->buf, "");
  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      siphead = SeqIdSetDup (bsp->id);
      for (sip = siphead; sip != NULL; sip = sip->next) {
        SeqIdStripLocus (sip);
      }
      SeqIdWrite (siphead, cfp->buf, PRINTID_FASTA_LONG, sizeof (cfp->buf));
      SeqIdSetFree (siphead);
    }
  }

  if (StringChr (cfp->report, 'r') != NULL) {
    DoASNReport (sep, cfp);
    return;
  }
  if (StringChr (cfp->report, 'g') != NULL) {
    DoGBFFReport (sep, cfp);
    return;
  }
  if (StringChr (cfp->report, 'm') != NULL) {
    DoModernizeReport (sep, cfp);
    return;
  }

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  if (StringChr (cfp->clean, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->clean, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
  }
  if (StringChr (cfp->clean, 'n') != NULL) {
    NormalizeDescriptorOrder (sep);
  }

  if (StringChr (cfp->modernize, 'g') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModGenes);
  }
  if (StringChr (cfp->modernize, 'r') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModRNAs);
  }
  if (StringChr (cfp->modernize, 'p') != NULL) {
    VisitBioSourcesInSep (sep, NULL, ModPCRs);
  }

  if (cfp->taxon) {
    Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  }

  if (cfp->pub) {
    VisitPubdescsInSep (sep, NULL, LookupPubdesc);
  }

  if (StringChr (cfp->link, 'o') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyOverlap (sep);
  }
  if (StringChr (cfp->link, 'p') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyProduct (sep);
  }
  if (StringChr (cfp->link, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    ReassignFeatureIDs (sep);
  }
  if (StringChr (cfp->link, 'c') != NULL) {
    ClearFeatureIDs (sep);
  }

  if (StringChr (cfp->feat, 'u') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatUser);
  }
  if (StringChr (cfp->feat, 'd') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatDbxref);
  }
  if (StringChr (cfp->feat, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
  }

  if (StringChr (cfp->desc, 't') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
  }

  if (StringChr (cfp->mods, 'd') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID);
  }
}

static void CleanupSingleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FILE          *fp;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  SeqEntryPtr   sep;

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  if (cfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", filename);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (cfp->type >= 2 && cfp->type <= 5) {
    aip = AsnIoOpen (filename, cfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    switch (cfp->type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) cfp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {

      path [0] = '\0';
      if (StringDoesHaveText (cfp->outfile)) {

        StringNCpy_0 (path, cfp->outfile, sizeof (path));
      
      } else if (StringDoesHaveText (cfp->results)) {

        ptr = StringRChr (filename, DIRDELIMCHR);
        if (ptr != NULL) {
          StringNCpy_0 (path, cfp->results, sizeof (path));
          ptr++;
          FileBuildPath (path, NULL, ptr);
        }
      }

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL && StringDoesHaveText (path)) {

        DoCleanup (sep, entityID, cfp);

        aop = AsnIoOpen (path, "w");
        if (aop != NULL) {
          if (datatype == OBJ_SEQSUB) {
            SeqSubmitAsnWrite ((SeqSubmitPtr) dataptr, aop, NULL);
          } else {
            SeqEntryAsnWrite (sep, aop, NULL);
          }
          AsnIoFlush (aop);
          AsnIoClose (aop);
        }
      }

      ObjMgrFreeByEntityID (entityID);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

static void CleanupMultipleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr     aip, aop;
  AsnTypePtr   atp;
  DataVal      av;
  Uint2        entityID;
  FILE         *fp;
  size_t       len;
  Char         longest [64];
  Int4         numrecords;
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  time_t       starttime, stoptime, worsttime;
#ifdef OS_UNIX
  Char         cmmd [256];
  CharPtr      gzcatprog;
  int          ret;
  Boolean      usedPopen = FALSE;
#endif

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  path [0] = '\0';
  if (StringDoesHaveText (cfp->outfile)) {

    StringNCpy_0 (path, cfp->outfile, sizeof (path));
      
  } else if (StringDoesHaveText (cfp->results)) {

    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      StringNCpy_0 (path, cfp->results, sizeof (path));
      ptr++;
      if (cfp->compressed) {
        len = StringLen (ptr);
        if (len > 4 && StringCmp (ptr + len - 3, ".gz") == 0) {
          ptr [len - 3] = '\0';
        }
      }
      FileBuildPath (path, NULL, ptr);
    }
  }
  if (StringHasNoText (path)) return;

#ifndef OS_UNIX
  if (cfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

#ifdef OS_UNIX
  if (cfp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* cfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, cfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, cfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (cfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n\n", filename);
    fflush (cfp->logfp);
  }

  longest [0] = '\0';
  worsttime = 0;
  numrecords = 0;

  aop = AsnIoOpen (path, cfp->binary? "wb" : "w");
  if (aop != NULL) {

    AsnOpenStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
    av.intvalue = 7;
    AsnWrite (aop, cfp->atp_bsc, &av);
    AsnOpenStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));

    atp = cfp->atp_bss;

    while ((atp = AsnReadId (aip, cfp->amp, atp)) != NULL) {
      if (atp == cfp->atp_se) {

        sep = SeqEntryAsnRead (aip, atp);
        if (sep != NULL) {

          entityID = ObjMgrGetEntityIDForChoice (sep);

          starttime = GetSecs ();
          DoCleanup (sep, entityID, cfp);
          stoptime = GetSecs ();

          if (stoptime - starttime > worsttime) {
            worsttime = stoptime - starttime;
            StringCpy (longest, cfp->buf);
          }
          numrecords++;

          SeqEntryAsnWrite (sep, aop, cfp->atp_se);

          ObjMgrFreeByEntityID (entityID);
        }

      } else {

        AsnReadVal (aip, atp, NULL);
      }
    }

    AsnCloseStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
    AsnCloseStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
  }

  AsnIoClose (aop);
  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "Total number of records %ld\n", (long) numrecords);
    if (StringDoesHaveText (longest)) {
      fprintf (cfp->logfp, "Longest processing time %ld seconds on %s\n",
               (long) worsttime, longest);
    }
    if (cfp->okay > 0 || cfp->norm > 0 || cfp->bsec > 0 || cfp->ssec > 0) {
      fprintf (cfp->logfp, "%ld OKAY, %ld NORM, %ld BSEC, %ld SSEC\n",
               (long) cfp->okay, (long) cfp->norm, (long) cfp->bsec, (long) cfp->ssec);
    }
    fflush (cfp->logfp);
  }
}

static void CleanupOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;

  if (StringHasNoText (filename)) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  cfp->okay = 0;
  cfp->bsec = 0;
  cfp->ssec = 0;
  cfp->norm = 0;

  if (cfp->batch) {
    CleanupMultipleRecord (filename, cfp);
  } else {
    CleanupSingleRecord (filename, cfp);
  }
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define r_argOutputPath    1
#define i_argInputFile     2
#define o_argOutputFile    3
#define f_argFilter        4
#define x_argSuffix        5
#define a_argType          6
#define b_argBinary        7
#define c_argCompressed    8
#define L_argLogFile       9
#define R_argRemote       10
#define Q_argReport       11
#define q_argFfDiff       12
#define m_argFfMode       13
#define K_argClean        14
#define U_argModernize    15
#define N_argLink         16
#define F_argFeat         17
#define D_argDesc         18
#define M_argMods         19
#define T_argTaxonLookup  20
#define P_argPubLookup    21

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Processing", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Report\n"
   "      r ASN.1 BSEC/SSEC Report\n"
   "      g GenBank SSEC Diff\n"
   "      m Modernize Gene/RNA/PCR", NULL, NULL, NULL,
    TRUE, 'Q', ARG_STRING, 0.0, 0, NULL},
  {"Ffdiff Executable", "/netopt/genbank/subtool/bin/ffdiff", NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"Flatfile Mode\n"
   "      r Release\n"
   "      e Entrez\n"
   "      s Sequin\n"
   "      d Dump\n", NULL, NULL, NULL,
    TRUE, 'm', ARG_STRING, 0.0, 0, NULL},
  {"Cleanup\n"
   "      b BasicSeqEntryCleanup\n"
   "      s SeriousSeqEntryCleanup\n"
   "      n Normalize Descriptor Order", NULL, NULL, NULL,
    TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
  {"Modernize\n"
   "      g Gene\n"
   "      r RNA\n"
   "      p PCR Primers", NULL, NULL, NULL,
    TRUE, 'U', ARG_STRING, 0.0, 0, NULL},
  {"Link\n"
   "      o LinkCDSmRNAbyOverlap\n"
   "      p LinkCDSmRNAbyProduct\n"
   "      r ReassignFeatureIDs\n"
   "      c ClearFeatureIDs", NULL, NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Feature\n"
   "      u Remove User Object\n"
   "      d Remove db_xref\n"
   "      r Remove Redundant Gene xref", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Descriptor\n"
   "      t Remove Title", NULL, NULL, NULL,
    TRUE, 'D', ARG_STRING, 0.0, 0, NULL},
  {"Miscellaneous\n"
   "      d Automatic Definition Line", NULL, NULL, NULL,
    TRUE, 'M', ARG_STRING, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Publication Lookup", "F", NULL, NULL,
    TRUE, 'P', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char           app [64], mode, type;
  CleanFlagData  cfd;
  CharPtr        directory, filter, infile, logfile, outfile, results, str, suffix;
  Boolean        remote;
  time_t         runtime, starttime, stoptime;

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

  /* process command line arguments */

  sprintf (app, "cleanasn %s", CLEANASN_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &cfd, 0, sizeof (CleanFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  cfd.batch = FALSE;
  cfd.binary = (Boolean) myargs [b_argBinary].intvalue;
  cfd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  cfd.type = 1;

  str = myargs [a_argType].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    type = str [0];
  } else {
    type = 'a';
  }

  type = TO_LOWER (type);
  switch (type) {
    case 'a' :
      cfd.type = 1;
      break;
    case 'e' :
      cfd.type = 2;
      break;
    case 'b' :
      cfd.type = 3;
      break;
    case 's' :
      cfd.type = 4;
      break;
    case 'm' :
      cfd.type = 5;
      break;
    case 't' :
      cfd.type = 1;
      cfd.batch = TRUE;
      break;
    default :
      cfd.type = 1;
      break;
  }

  remote = (Boolean) myargs [R_argRemote].intvalue;

  cfd.report = myargs [Q_argReport].strvalue;
  cfd.ffdiff = myargs [q_argFfDiff].strvalue;

  str = myargs [m_argFfMode].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    mode = str [0];
  } else {
    mode = 'e';
  }

  mode = TO_LOWER (mode);
  switch (mode) {
    case 'r' :
      cfd.ffmode = RELEASE_MODE;
      break;
    case 'e' :
      cfd.ffmode = ENTREZ_MODE;
      break;
    case 's' :
      cfd.ffmode = SEQUIN_MODE;
      break;
    case 'd' :
      cfd.ffmode = DUMP_MODE;
      break;
    default :
      cfd.ffmode = ENTREZ_MODE;
      break;
  }

  cfd.clean = myargs [K_argClean].strvalue;
  cfd.modernize = myargs [U_argModernize].strvalue;
  cfd.link = myargs [N_argLink].strvalue;
  cfd.feat = myargs [F_argFeat].strvalue;
  cfd.desc = myargs [D_argDesc].strvalue;
  cfd.mods = myargs [M_argMods].strvalue;
  cfd.taxon = (Boolean) myargs [T_argTaxonLookup].intvalue;
  cfd.pub = (Boolean) myargs [P_argPubLookup].intvalue;

  cfd.amp = AsnAllModPtr ();
  cfd.atp_bss = AsnFind ("Bioseq-set");
  cfd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  cfd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  cfd.atp_bsc = AsnFind ("Bioseq-set.class");
  cfd.bssp_atp = AsnLinkType (NULL, cfd.atp_bss);

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;
  if (StringDoesHaveText (logfile)) {
    cfd.logfp = FileOpen (logfile, "w");
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    if (! PUBSEQBioseqFetchEnable ("cleanasn", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
#else
    PubSeqFetchEnable ();
#endif
  }

  if (remote || cfd.pub) {
    PubMedFetchEnable ();
  }

  starttime = GetSecs ();

  if (StringDoesHaveText (directory)) {
    if (StringCmp (directory, results) == 0) {
      Message (MSG_POSTERR, "-r results path must be different than -p data path");
      if (cfd.logfp != NULL) {
        fprintf (cfd.logfp, "-r results path must be different than -p data path\n");
      }
    } else {

      cfd.results = results;

      DirExplore (directory, NULL, suffix, FALSE, CleanupOneRecord, (Pointer) &cfd);
    }

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    cfd.outfile = outfile;

    CleanupOneRecord (infile, (Pointer) &cfd);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;
  if (cfd.logfp != NULL) {
    fprintf (cfd.logfp, "Finished in %ld seconds\n", (long) runtime);
    if (cfd.cumokay > 0 || cfd.cumnorm > 0 || cfd.cumbsec > 0 || cfd.cumssec > 0) {
      fprintf (cfd.logfp, "Cumulative counts - %ld OKAY, %ld NORM, %ld BSEC, %ld SSEC\n",
               (long) cfd.cumokay, (long) cfd.cumnorm, (long) cfd.cumbsec, (long) cfd.cumssec);
    }
    FileClose (cfd.logfp);
  }

  if (remote || cfd.pub) {
    PubMedFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return 0;
}

