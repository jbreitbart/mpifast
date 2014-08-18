/*   asndisc.c
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
* File Name:  asndisc.c
*
* Author:  Jonathan Kans, adapted from asnval.c by Colleen Bollin
*
* Version Creation Date:   1/23/07
*
* $Revision: 1.25 $
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
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <gather.h>
#include <explore.h>
#include <lsqfetch.h>
#include <valid.h>
#include <pmfapi.h>
#ifdef INTERNAL_NCBI_ASNDISC
#include <accpubseq.h>
#include <tax3api.h>
#endif

#define ASNDISC_APP_VER "1.1"

CharPtr ASNDISC_APPLICATION = ASNDISC_APP_VER;

typedef struct drflags {
  Boolean  farFetchCDSproducts;
  Boolean  batch;
  Boolean  binary;
  Boolean  compressed;
  Boolean  lock;
  Boolean  useThreads;
  Boolean  usePUBSEQ;
  Int2     type;
  Int4     maxcount;
  CharPtr  outpath;
  CharPtr  output_suffix;
  CharPtr  output_dir;
  FILE     *outfp;
  Int4     numrecords;
  ValNodePtr            sep_list;
  ValNodePtr            bsplist;

  GlobalDiscrepReportPtr global_report;
} DRFlagData, PNTR DRFlagPtr;

#ifdef INTERNAL_NCBI_ASNDISC
const PerformDiscrepancyTest taxlookup = CheckTaxNamesAgainstTaxDatabase;
#else
const PerformDiscrepancyTest taxlookup = NULL;
#endif

#ifdef INTERNAL_NCBI_ASNDISC
static CharPtr dirsubfetchproc = "DirSubBioseqFetch";

static CharPtr dirsubfetchcmd = NULL;

extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK DirSubBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean DirSubFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, dirsubfetchproc, dirsubfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  DirSubBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr smartfetchproc = "SmartBioseqFetch";

static CharPtr smartfetchcmd = NULL;

extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK SmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr tpasmartfetchproc = "TPASmartBioseqFetch";

static CharPtr tpasmartfetchcmd = NULL;

extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK TPASmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_TPG) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean TPASmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, tpasmartfetchproc, tpasmartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  TPASmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}
#endif

static ValNodePtr DoLockFarComponents (
  SeqEntryPtr sep,
  DRFlagPtr drfp
)

{
  ValNodePtr  rsult;

#ifdef INTERNAL_NCBI_ASNDISC
  if (drfp->useThreads) {
    Message (MSG_POST, "Threads will not be used in this executable");
    drfp->useThreads = FALSE;;
  }
#endif

  if (NlmThreadsAvailable () && drfp->useThreads) {
    rsult = AdvcLockFarComponents (sep, TRUE, drfp->farFetchCDSproducts, drfp->farFetchCDSproducts, NULL, TRUE);
  } else if (drfp->useThreads) {
    Message (MSG_POST, "Threads not available in this executable");
    rsult = AdvcLockFarComponents (sep, TRUE, drfp->farFetchCDSproducts, drfp->farFetchCDSproducts, NULL, FALSE);
  } else {
    rsult = AdvcLockFarComponents (sep, TRUE, drfp->farFetchCDSproducts, drfp->farFetchCDSproducts, NULL, FALSE);
  }

  return rsult;
}


static void ReleaseDiscrepancyReportSeqEntries (DRFlagPtr drfp)
{
  ValNodePtr vnp;
  SeqEntryPtr sep;
  ObjMgrPtr   omp;

  if (drfp == NULL) {
    return;
  }

  for (vnp = drfp->sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    SeqEntryFree (sep);
    omp = ObjMgrGet ();
    ObjMgrReapOne (omp);
  }
  SeqMgrClearBioseqIndex ();
  ObjMgrFreeCache (0);
  FreeSeqIdGiCache ();
  SeqEntrySetScope (NULL);
  drfp->sep_list = ValNodeFree (drfp->sep_list);
  
  drfp->bsplist = UnlockFarComponents (drfp->bsplist);
}


static void ProcessSeqEntryList (DRFlagPtr drfp, CharPtr filename)
{
  ValNodePtr  discrepancy_list;
  FILE        *ofp = NULL;
  Char        path [PATH_MAX];
  CharPtr     ptr;

  if (drfp == NULL || drfp->sep_list == NULL) return;

  if (StringDoesHaveText (drfp->output_dir)) {
    if (StringLen (drfp->output_dir) > PATH_MAX) {
      Message (MSG_ERROR, "Unable to generate output file - path name is too long");
      return;
    }
    StringCpy (path, drfp->output_dir);
#ifdef OS_WINNT
    ptr = StringRChr (filename, '\\');
    if (path[StringLen(path) - 1] != '\\') {
      StringCat (path, "\\");
    }
#else
    ptr = StringRChr (filename, '/');
    if (path[StringLen(path) - 1] != '/') {
      StringCat (path, "/");
    }
#endif;
    if (ptr == NULL) {
      StringNCat (path, filename, PATH_MAX - StringLen(path) - 1);
    } else {
      StringNCat (path, ptr + 1, PATH_MAX - StringLen(path) - 1);
    }
  } else {
    StringNCpy_0 (path, filename, sizeof (path));
  }
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  if (StringDoesHaveText (drfp->output_suffix)) {
    StringNCat (path, drfp->output_suffix, PATH_MAX - StringLen(path) - 1);
    path[PATH_MAX - 1] = 0;
  } else {
    StringCat (path, ".dr");
  }
  ofp = FileOpen (path, "w");

  discrepancy_list = CollectDiscrepancies (drfp->global_report->test_config, drfp->sep_list, taxlookup);
  WriteAsnDiscReport (discrepancy_list, ofp, drfp->global_report->output_config, TRUE);
  discrepancy_list = FreeClickableList (discrepancy_list);

  FileClose (ofp);
}


static void ProcessSingleRecord (
  CharPtr filename,
  DRFlagPtr drfp
)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  ValNodePtr     bsplist_next = NULL;
  BioseqSetPtr   bssp;
  Char           path [PATH_MAX];
  Pointer        dataptr = NULL;
  Uint2          datatype, entityID = 0;
  FILE           *fp;
  SeqEntryPtr    sep;

  if (StringHasNoText (filename)) return;
  if (drfp == NULL) return;

  if (drfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", path);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (drfp->type >= 2 && drfp->type <= 5) {
    aip = AsnIoOpen (filename, drfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    switch (drfp->type) {
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
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) drfp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (SeqMgrFeaturesAreIndexed(entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
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
      ValNodeAddPointer (&(drfp->sep_list), 0, sep);

      if (drfp->lock) {
        bsplist_next = DoLockFarComponents (sep, drfp);
        ValNodeLink (&(drfp->bsplist), bsplist_next);
      }
    }
  } else {
    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  SeqEntrySetScope (NULL);
}

static void ProcessMultipleRecord (
  CharPtr filename,
  DRFlagPtr drfp
)

{
  AsnIoPtr        aip;
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_bss, atp_desc, atp_sbp, atp_se, atp_ssp;
  ValNodePtr      bsplist_next;
  Int2            maxcount = 0;
  CitSubPtr       csp = NULL;
  FILE            *fp, *ofp = NULL;
  Int4            numrecords = 0;
  SeqEntryPtr     sep;
  ObjValNode      ovn;
  Pubdesc         pd;
  SubmitBlockPtr  sbp = NULL;
  SeqDescrPtr     subcit = NULL;
  ValNode         vn;
#ifdef OS_UNIX
  Char            cmmd [256];
  Boolean         detailed_report = FALSE;
  CharPtr         gzcatprog;
  Boolean         memory_usage = FALSE;
  int             ret;
  Boolean         usedPopen = FALSE;
#endif

  if (StringHasNoText (filename)) return;
  if (drfp == NULL) return;

#ifndef OS_UNIX
  if (drfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_POSTERR, "Unable to load AsnAllModPtr");
    return;
  }

  atp_ssp = AsnFind ("Seq-submit");
  if (atp_ssp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return;
  }

  atp_sbp = AsnFind ("Seq-submit.sub");
  if (atp_sbp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.sub");
    return;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set");
    return;
  }

  atp_desc = AsnFind ("Bioseq-set.descr");
  if (atp_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.descr");
    return;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return;
  }

#ifdef OS_UNIX
  if (getenv ("ASNVAL_LOG_OBJMGR_REPORT") != NULL) {
    detailed_report = TRUE;
  }
  if (getenv ("ASNVAL_LOG_MEMORY_REPORT") != NULL) {
    memory_usage = TRUE;
  }

  if (drfp->compressed) {
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
    fp = popen (cmmd, /* drfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, drfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, drfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (drfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (drfp->type == 4) {
    atp = atp_bss;
  } else if (drfp->type == 5) {
    atp = atp_ssp;
  } else {
    Message (MSG_ERROR, "Batch processing type not set properly");
    return;
  }

  while ((atp = AsnReadId (aip, amp, atp)) != NULL && maxcount < drfp->maxcount) {
    if (atp == atp_se) {
      sep = SeqEntryAsnRead (aip, atp);
      ValNodeAddPointer (&(drfp->sep_list), 0, sep);

      if (drfp->lock) {
        bsplist_next = DoLockFarComponents (sep, drfp);
        ValNodeLink (&(drfp->bsplist), bsplist_next);
      }

      numrecords++;
      maxcount++;
    } else if (atp == atp_sbp) {
      sbp = SubmitBlockAsnRead (aip, atp);
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          MemSet ((Pointer) &ovn, 0, sizeof (ObjValNode));
          MemSet ((Pointer) &pd, 0, sizeof (Pubdesc));
          MemSet ((Pointer) &vn, 0, sizeof (ValNode));
          vn.choice = PUB_Sub;
          vn.data.ptrvalue = (Pointer) csp;
          vn.next = NULL;
          pd.pub = &vn;
          ovn.vn.choice = Seq_descr_pub;
          ovn.vn.data.ptrvalue = (Pointer) &pd;
          ovn.vn.next = NULL;
          ovn.vn.extended = 1;
          subcit = (SeqDescrPtr) &ovn;
        }
      }
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }



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

}


static void ProcessSeqEntryListWithCollation (GlobalDiscrepReportPtr g, ValNodePtr sep_list, CharPtr filename)
{
  ValNodePtr  vnp;
  SeqEntryPtr sep;

  if (g == NULL || sep_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    AddSeqEntryToGlobalDiscrepReport (sep, g, filename);
  }

}


static void ProcessOneRecord (CharPtr filename, Pointer userdata)
{
  DRFlagPtr  drfp;

  drfp = (DRFlagPtr) userdata;
  if (drfp == NULL) return;

  if (drfp->batch) {
    ProcessMultipleRecord (filename, drfp);
  } else {
    ProcessSingleRecord (filename, drfp);
  }

  if (drfp->outfp == NULL) {
    ProcessSeqEntryList (drfp, filename);
  } else {
    ProcessSeqEntryListWithCollation (drfp->global_report, drfp->sep_list, filename);
  }
  ReleaseDiscrepancyReportSeqEntries (drfp);
}


/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  i_argInputFile,
  o_argOutputFile,
  x_argSuffix,
  u_argRecurse,
  f_argUseFT,
  e_argEnableTests,
  d_argDisableTests,
  s_argOutputSuffix,
  r_argOutputDir,
  Z_argRemoteCDS,
  a_argType,
  b_argBinary,
  c_argCompressed,
  R_argRemote,
  k_argLocalFetch,
  I_argAsnIdx,
  l_argLockFar,
  T_argThreads,
  X_argExpandCategories,
  S_argSummaryReport,
  B_argBigSequenceReport,
  C_argMaxCount
} DRFlagNum;

Args myargs [] = {
  {"Path to ASN.1 Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"File Selection Substring", ".sqn", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Feature Table Output Format", "F", NULL, NULL,
    FALSE, 'f', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Enable Tests (comma-delimited list of test names)\n\tMISSING_GENES\n\tEXTRA_GENES\n\tMISSING_LOCUS_TAGS\n\tDUPLICATE_LOCUS_TAGS\n\tBAD_LOCUS_TAG_FORMAT\n"
   "\tINCONSISTENT_LOCUS_TAG_PREFIX\n\tNON_GENE_LOCUS_TAG\n\tMISSING_PROTEIN_ID\n\tINCONSISTENT_PROTEIN_ID\n"
   "\tFEATURE_LOCATION_CONFLICT\n\tGENE_PRODUCT_CONFLICT\n\tDUPLICATE_GENE_LOCUS\n\tEC_NUMBER_NOTE\n\tPSEUDO_MISMATCH\n"
   "\tJOINED_FEATURES\n\tOVERLAPPING_GENES\n\tOVERLAPPING_CDS\n\tSHORT_CONTIG\n\tINCONSISTENT_BIOSOURCE\n\tSUSPECT_PRODUCT_NAMES\n"
   "\tINCONSISTENT_SOURCE_DEFLINE\n\tPARTIAL_CDS_COMPLETE_SEQUENCE\n\tEC_NUMBER_ON_UNKNOWN_PROTEIN\n\tTAX_LOOKUP_MISSING\n"
   "\tTAX_LOOKUP_MISMATCH\n\tSHORT_SEQUENCES\n\tSUSPECT_PHRASES\n", "", NULL, NULL,
    TRUE, 'e', ARG_STRING, 0.0, 0, NULL},
  {"Disable Tests (comma-delimited list of test names)\n\tMISSING_GENES\n\tEXTRA_GENES\n\tMISSING_LOCUS_TAGS\n\tDUPLICATE_LOCUS_TAGS\n\tBAD_LOCUS_TAG_FORMAT\n"
   "\tINCONSISTENT_LOCUS_TAG_PREFIX\n\tNON_GENE_LOCUS_TAG\n\tMISSING_PROTEIN_ID\n\tINCONSISTENT_PROTEIN_ID\n"
   "\tFEATURE_LOCATION_CONFLICT\n\tGENE_PRODUCT_CONFLICT\n\tDUPLICATE_GENE_LOCUS\n\tEC_NUMBER_NOTE\n\tPSEUDO_MISMATCH\n"
   "\tJOINED_FEATURES\n\tOVERLAPPING_GENES\n\tOVERLAPPING_CDS\n\tSHORT_CONTIG\n\tINCONSISTENT_BIOSOURCE\n\tSUSPECT_PRODUCT_NAMES\n"
   "\tINCONSISTENT_SOURCE_DEFLINE\n\tPARTIAL_CDS_COMPLETE_SEQUENCE\n\tEC_NUMBER_ON_UNKNOWN_PROTEIN\n\tTAX_LOOKUP_MISSING\n"
   "\tTAX_LOOKUP_MISMATCH\n\tSHORT_SEQUENCES\n\tSUSPECT_PHRASES\n", "", NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Output File Suffix", ".dr", NULL, NULL,
    TRUE, 's', ARG_STRING, 0.0, 0, NULL},
  {"Output Directory", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Remote CDS Product Fetch", "F", NULL, NULL,
    TRUE, 'Z', ARG_BOOLEAN, 0.0, 0, NULL},
  {"ASN.1 Type (a Any, e Seq-entry, b Bioseq, s Bioseq-set, m Seq-submit, t Batch Bioseq-set, u Batch Seq-submit)", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Batch File is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch File is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Local Fetching", "F", NULL, NULL,
    TRUE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Path to Indexed Binary ASN.1 Data", NULL, NULL, NULL,
    TRUE, 'I', ARG_STRING, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Expand Report Categories (comma-delimited list of test names or ALL)\n\tALL\n\tMISSING_GENES\n\tEXTRA_GENES\n\tMISSING_LOCUS_TAGS\n\tDUPLICATE_LOCUS_TAGS\n\tBAD_LOCUS_TAG_FORMAT\n"
   "\tINCONSISTENT_LOCUS_TAG_PREFIX\n\tNON_GENE_LOCUS_TAG\n\tMISSING_PROTEIN_ID\n\tINCONSISTENT_PROTEIN_ID\n"
   "\tFEATURE_LOCATION_CONFLICT\n\tGENE_PRODUCT_CONFLICT\n\tDUPLICATE_GENE_LOCUS\n\tEC_NUMBER_NOTE\n\tPSEUDO_MISMATCH\n"
   "\tJOINED_FEATURES\n\tOVERLAPPING_GENES\n\tOVERLAPPING_CDS\n\tSHORT_CONTIG\n\tINCONSISTENT_BIOSOURCE\n\tSUSPECT_PRODUCT_NAMES\n"
   "\tINCONSISTENT_SOURCE_DEFLINE\n\tPARTIAL_CDS_COMPLETE_SEQUENCE\n\tEC_NUMBER_ON_UNKNOWN_PROTEIN\n\tTAX_LOOKUP_MISSING\n"
   "\tTAX_LOOKUP_MISMATCH\n\tSHORT_SEQUENCES\n\tSUSPECT_PHRASES\n", "", NULL, NULL,
    TRUE, 'X', ARG_STRING, 0.0, 0, NULL},
  {"Summary Report", "F", NULL, NULL,
    TRUE, 'S', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Big Sequence Report", "F", NULL, NULL,
  TRUE, 'B', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Max Count", "0", NULL, NULL,
    TRUE, 'C', ARG_INT, 0.0, 0, NULL},
};


static CharPtr GetTestNameList (CharPtr intro)
{
  Int4 i, len;
  CharPtr text;

  len = StringLen (intro) + 1;

  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    len += StringLen (GetDiscrepancyTestSettingName (i)) + 2;
  }

  text = (CharPtr) MemNew (sizeof (Char) * len);
  StringCat (text, intro);
  for (i = 0; i < MAX_DISC_TYPE; i++) {
    StringCat (text, "\t");
    StringCat (text, GetDiscrepancyTestSettingName (i));
    StringCat (text, "\n");
  }
  return text;
}


Int2 Main (void)

{
  Char         app [64];
  CharPtr      asnidx, directory, infile, outfile, str, suffix, output_dir;
  CharPtr      enabled_list, disabled_list, err_msg;
  Boolean      batch, binary, compressed, dorecurse,
               indexed, local, lock, remote, usethreads;
  Int2         type = 0;
  DRFlagData   dfd;
  Boolean      big_sequence_report;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  ErrSetOpts (ERR_IGNORE, ERR_LOG_ON);

  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

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

  /* set up help descriptions for enable and disable */
  myargs[e_argEnableTests].prompt = GetTestNameList("Enable Tests (comma-delimited list of test names)\n");
  myargs[d_argDisableTests].prompt = GetTestNameList("Disable Tests (comma-delimited list of test names)\n");
  myargs[X_argExpandCategories].prompt = GetTestNameList("Expand Report Categories (comma-delimited list of test names or ALL)\n");
  /* process command line arguments */

  sprintf (app, "asndisc %s", ASNDISC_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* additional setup modifications */
  MemSet (&dfd, 0, sizeof (DRFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  dfd.output_suffix = (CharPtr) myargs [s_argOutputSuffix].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  output_dir = (CharPtr) myargs [r_argOutputDir].strvalue;
  if (StringDoesHaveText (outfile) && StringDoesHaveText (output_dir)) {
    Message (MSG_FATAL, "-o and -q are incompatible: specify the output file name with the full path.");
    return 1;
  }
  if (StringDoesHaveText (output_dir)) {
    dfd.output_dir = output_dir;
    if (! CreateDir (output_dir)) {
      Message (MSG_FATAL, "Unable to create output directory %s", output_dir);
    }
  }

  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  remote = (Boolean ) myargs [R_argRemote].intvalue;
  local = (Boolean) myargs [k_argLocalFetch].intvalue;

  asnidx = (CharPtr) myargs [I_argAsnIdx].strvalue;
  indexed = (Boolean) StringDoesHaveText (asnidx);
  lock = (Boolean) myargs [l_argLockFar].intvalue;
  usethreads = (Boolean) myargs [T_argThreads].intvalue;
  dfd.farFetchCDSproducts = (Boolean) myargs [Z_argRemoteCDS].intvalue;

  /* set up Discrepancy Report Configuration */
  dfd.global_report = GlobalDiscrepReportNew ();
  dfd.global_report->test_config = DiscrepancyConfigNew();
  DisableTRNATests (dfd.global_report->test_config);

  ExpandDiscrepancyReportTestsFromString ((CharPtr) myargs [X_argExpandCategories].strvalue, TRUE, dfd.global_report->output_config);
  dfd.global_report->output_config->summary_report = (Boolean) myargs [S_argSummaryReport].intvalue;

  big_sequence_report = (Boolean) myargs [B_argBigSequenceReport].intvalue;

  enabled_list = (CharPtr) myargs [e_argEnableTests].strvalue;
  disabled_list = (CharPtr) myargs [d_argDisableTests].strvalue;


#ifdef INTERNAL_NCBI_ASNDISC
  dfd.global_report->taxlookup = CheckTaxNamesAgainstTaxDatabase;
#endif
  
  err_msg = NULL;
  if (StringDoesHaveText (enabled_list) && StringDoesHaveText (disabled_list)) {
    err_msg = StringSave ("Cannot specify both -e and -d.  Choose -e to enable only a few tests and disable the rest, choose -d to disable only a few tests and enable the rest.");
  } else if (StringDoesHaveText (disabled_list)) {
    if (big_sequence_report) {
      ConfigureForBigSequence (dfd.global_report->test_config);
    } else {
      ConfigureForGenomes (dfd.global_report->test_config);
    }

    /* now disable tests from string */
    err_msg = SetDiscrepancyReportTestsFromString (disabled_list, FALSE, dfd.global_report->test_config);
  } else if (StringDoesHaveText (enabled_list)) {
    if (big_sequence_report) {
      ConfigureForBigSequence (dfd.global_report->test_config);
    } else {
      ConfigureForGenomes (dfd.global_report->test_config);
    }
    /* now enable tests from string */
    err_msg = SetDiscrepancyReportTestsFromString (enabled_list, TRUE, dfd.global_report->test_config);
  } else {
    if (big_sequence_report) {
      ConfigureForBigSequence (dfd.global_report->test_config);
    } else {
      ConfigureForGenomes (dfd.global_report->test_config);
    }
  }
  if (err_msg != NULL) {
    Message (MSG_FATAL, err_msg);
    err_msg = MemFree (err_msg);
    return 1;
  }

  if ((Boolean) myargs[f_argUseFT].intvalue) {
    dfd.global_report->test_config->use_feature_table_format = TRUE;
    dfd.global_report->output_config->use_feature_table_format = TRUE;
  }

  dfd.maxcount = (Int4) myargs [C_argMaxCount].intvalue;
  if (dfd.maxcount < 1) {
    dfd.maxcount = INT4_MAX;
  }

  batch = FALSE;
  binary = (Boolean) myargs [b_argBinary].intvalue;
  compressed = (Boolean) myargs [c_argCompressed].intvalue;

  str = myargs [a_argType].strvalue;
  if (StringICmp (str, "a") == 0) {
    type = 1;
  } else if (StringICmp (str, "e") == 0) {
    type = 2;
  } else if (StringICmp (str, "b") == 0) {
    type = 3;
  } else if (StringICmp (str, "s") == 0) {
    type = 4;
  } else if (StringICmp (str, "m") == 0) {
    type = 5;
  } else if (StringICmp (str, "t") == 0) {
    type = 4;
    batch = TRUE;
  } else if (StringICmp (str, "u") == 0) {
    type = 5;
    batch = TRUE;
  } else {
    type = 1;
  }

  if ((binary || compressed) && (! batch)) {
    if (type == 1) {
      Message (MSG_FATAL, "-b or -c cannot be used without -t or -a");
      return 1;
    }
  }

  if (StringHasNoText (directory) && StringHasNoText (infile)) {
    Message (MSG_FATAL, "Input path or input file must be specified");
    return 1;
  }

  /* populate parameter structure */

  dfd.batch = batch;
  dfd.binary = binary;
  dfd.compressed = compressed;
  dfd.lock = lock;
  dfd.useThreads = usethreads;
  dfd.type = type;
  dfd.numrecords = 0;

  if (! StringHasNoText (outfile)) {
    dfd.outpath = outfile;
    dfd.outfp = FileOpen (outfile, "w");
    if (dfd.outfp == NULL) {
      Message (MSG_FATAL, "Unable to open single output file");
      return 1;
    }
  }

  /* register fetch functions */

  if (remote) {
#ifdef INTERNAL_NCBI_ASNDISC

    if (! PUBSEQBioseqFetchEnable ("asnval", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    dfd.usePUBSEQ = TRUE;
    dfd.useThreads = FALSE;
#else
    PubSeqFetchEnable ();
#endif
  }

  if (local) {
    LocalSeqFetchInit (FALSE);
  }

  if (indexed) {
    AsnIndexedLibFetchEnable (asnidx, TRUE);
  }

  if (StringDoesHaveText (directory)) {
    DirExplore (directory, NULL, suffix, dorecurse, ProcessOneRecord, (Pointer) &dfd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, (Pointer) &dfd);
  }
  if (dfd.outfp != NULL) {
    WriteGlobalDiscrepancyReport (dfd.global_report, dfd.outfp);
    FileClose (dfd.outfp);
    dfd.outfp = NULL;
  }

  dfd.global_report = GlobalDiscrepReportFree (dfd.global_report);

  /* close fetch functions */

  if (indexed) {
    AsnIndexedLibFetchDisable ();
  }

  if (local) {
    LocalSeqFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_ASNDISC
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
    SeqMgrSetPreCache (NULL);
    SeqMgrSetSeqIdSetFunc (NULL);
  }

  TransTableFreeAll ();

  ECNumberFSAFreeAll ();

  return 0;
}

