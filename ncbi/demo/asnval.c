/*   asnval.c
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
* File Name:  asnval.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/3/04
*
* $Revision: 1.80 $
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
#ifdef INTERNAL_NCBI_ASN2VAL
#include <accpubseq.h>
#endif

#define ASNVAL_APP_VER "6.2"

CharPtr ASNVAL_APPLICATION = ASNVAL_APP_VER;

typedef struct valflags {
  Int2     severity;
  Int2     lowCutoff;
  Int2     highCutoff;
  CharPtr  errcode;
  Boolean  validateAlignments;
  Boolean  alignFindRemoteBsp;
  Boolean  doSeqHistAssembly;
  Boolean  farIDsInAlignments;
  Boolean  alwaysRequireIsoJTA;
  Boolean  farFetchCDSproducts;
  Boolean  farFetchMRNAproducts;
  Boolean  locusTagGeneralMatch;
  Boolean  validateIDSet;
  Boolean  seqSubmitParent;
  Boolean  ignoreExceptions;
  Boolean  validateExons;
  Boolean  inferenceAccnCheck;
  Boolean  testLatLonSubregion;
  Boolean  strictLatLonCountry;
  Boolean  indexerVersion;
  Boolean  batch;
  Boolean  binary;
  Boolean  compressed;
  Boolean  lock;
  Boolean  useThreads;
  Boolean  usePUBSEQ;
  Boolean  validateBarcode;
  Int2     verbosity;
  Int2     type;
  Int4     skipcount;
  Int4     maxcount;
  CharPtr  outpath;
  FILE     *outfp;
  FILE     *logfp;
  Int4     num_errors;
  Int4     fatal_errors;
  Boolean  has_errors;
  Boolean  io_failure;
  Char     longest [64];
  time_t   worsttime;
  Int4     numrecords;
} ValFlagData, PNTR ValFlagPtr;

#ifdef INTERNAL_NCBI_ASN2VAL
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
  ValFlagPtr vfp
)

{
  Boolean     farFetch;
  ValNodePtr  rsult;
  time_t      start_time, stop_time;

  start_time = GetSecs ();

#ifdef INTERNAL_NCBI_ASN2VAL
  if (vfp->useThreads) {
    Message (MSG_POST, "Threads will not be used in this executable");
    vfp->useThreads = FALSE;;
  }
#endif

  farFetch = (Boolean) (vfp->farFetchCDSproducts);

  if (NlmThreadsAvailable () && vfp->useThreads) {
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, TRUE);
  } else if (vfp->useThreads) {
    Message (MSG_POST, "Threads not available in this executable");
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, FALSE);
  } else {
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, FALSE);
  }

  stop_time = GetSecs ();

  return rsult;
}

static CharPtr severityLabel [] = {
  "NONE", "INFO", "WARNING", "ERROR", "REJECT", "FATAL", "MAX", NULL
};

static CharPtr compatSeverityLabel [] = {
  "NONE", "NOTE: valid", "WARNING: valid", "ERROR: valid", "REJECT: valid", "FATAL: valid", "MAX", NULL
};

typedef struct vcdaa {
  FILE        *ofp;
  Int2        verbosity;
  Int2        lowCutoff;
  Int2        highCutoff;
  CharPtr     errcode;
  ValFlagPtr  vfp;
} VCData, PNTR VCPtr;

static void XmlEncode (CharPtr dst, CharPtr src)

{
  Char  ch;

  if (dst == NULL || src == NULL) return;

  ch = *src;
  while (ch != '\0') {
    if (ch == '<') {
      *dst = '&';
      dst++;
      *dst = 'l';
      dst++;
      *dst = 't';
      dst++;
      *dst = ';';
      dst++;
    } else if (ch == '>') {
      *dst = '&';
      dst++;
      *dst = 'g';
      dst++;
      *dst = 't';
      dst++;
      *dst = ';';
      dst++;
    } else {
      *dst = ch;
      dst++;
    }
    src++;
    ch = *src;
  }
  *dst = '\0';
}


static CharPtr GetXmlHeaderText (ErrSev cutoff)
{
  CharPtr         xml_header = NULL;
  CharPtr         xml_4_fmt = "asnval version=\"%s\" severity_cutoff=\"%s\"";

  xml_header = (CharPtr) MemNew (sizeof (Char) * (10 + StringLen (xml_4_fmt) + 
                StringLen (ASNVAL_APPLICATION) + StringLen (severityLabel[cutoff])));
  sprintf (xml_header, xml_4_fmt, ASNVAL_APPLICATION, severityLabel[cutoff]);
  return xml_header;
}


static void LIBCALLBACK ValidCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  Char        buf [256];
  CharPtr     catname, errname, urlmssg = NULL;
  ErrSev      cutoff;
  FILE        *fp;
  size_t      len;
  VCPtr       vcp;
  ValFlagPtr  vfp;
  CharPtr     xml_header;

  vcp = (VCPtr) userdata;
  if (vcp == NULL) return;
  fp = vcp->ofp;
  if (fp == NULL) return;
  vfp = vcp->vfp;
  if (vfp == NULL) return;

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  if (severity < vcp->lowCutoff || severity > vcp->highCutoff) return;

  catname = GetValidCategoryName (errcode);
  errname = GetValidErrorName (errcode, subcode);

  if (catname == NULL) {
    catname = "?";
  }
  if (errname == NULL) {
    errname = "?";
  }

  if (StringDoesHaveText (vcp->errcode)) {
    if (StringICmp (vcp->errcode, errname) != 0) return;
  }

  if (accession == NULL) {
    accession = "";
  }
  if (message == NULL) {
    message = "";
  }
  if (objtype == NULL) {
    objtype = "";
  }
  if (label == NULL) {
    label = "";
  }

  if (vcp->verbosity == 1) {

    fprintf (fp, "%s [%s.%s] %s %s: %s",
             compatSeverityLabel [severity],
             catname, errname, message, objtype, label);
    if (location != NULL) {
      fprintf (fp, " %s", location);
    }
    if (context != NULL) {
      fprintf (fp, " %s", context);
    }
    if (product != NULL) {
      fprintf (fp, " -> %s", product);
    }
    fprintf (fp, "\n");

  } else if (vcp->verbosity == 2) {

    StringCpy (buf, accession);
    StringCat (buf, "                    ");
    buf [15] = '\0';

    StringCat (buf, severityLabel [severity]);
    StringCat (buf, "                    ");
    buf [30] = '\0';

    StringCat (buf, catname);
    StringCat (buf, "_");
    StringCat (buf, errname);

    fprintf (fp, "%s\n", buf);

  } else if (vcp->verbosity == 3) {

    fprintf (fp, "%s\t%s\t%s_%s\n",
             accession, severityLabel [severity],
             catname, errname);

  } else if (vcp->verbosity == 4) {

    if (! vfp->has_errors) {
      cutoff = (ErrSev) vcp->lowCutoff;
      if (cutoff < SEV_NONE || cutoff > SEV_MAX) {
        cutoff = SEV_MAX;
      }

      xml_header = GetXmlHeaderText (cutoff);
      fprintf (fp, "<%s>\n", xml_header);
      xml_header = MemFree (xml_header);
    }

    len = StringLen (message);
    if (len > 0) {
      urlmssg = MemNew (len * 3 + 2);
      if (urlmssg != NULL) {
        XmlEncode (urlmssg, message);
        fprintf (fp, "  <message severity=\"%s\" seq-id=\"%s\" code=\"%s_%s\">%s</message>\n",
                 severityLabel [severity], accession, catname, errname, urlmssg);
        MemFree (urlmssg);
      }
    }
  }

  vfp->has_errors = TRUE;
}

static void DoValidation (
  SeqEntryPtr sep,
  ValFlagPtr vfp,
  FILE *ofp
)

{
  Int2            i;
  VCData          vcd;
  ValidStructPtr  vsp;
  ErrSev          cutoff;
  CharPtr         xml_header = NULL;

  if (vfp == NULL) return;

  vsp = ValidStructNew ();
  if (vsp == NULL) return;

  MemSet ((Pointer) &vcd, 0, sizeof (VCData));

  vsp->useSeqMgrIndexes = TRUE;

  vsp->cutoff = vfp->lowCutoff;
  vsp->validateAlignments = vfp->validateAlignments;
  vsp->alignFindRemoteBsp = vfp->alignFindRemoteBsp;
  vsp->doSeqHistAssembly = vfp->doSeqHistAssembly;
  vsp->farIDsInAlignments = vfp->farIDsInAlignments;
  vsp->alwaysRequireIsoJTA = vfp->alwaysRequireIsoJTA;
  vsp->farFetchCDSproducts = vfp->farFetchCDSproducts;
  vsp->farFetchMRNAproducts = vfp->farFetchMRNAproducts;
  vsp->locusTagGeneralMatch = vfp->locusTagGeneralMatch;
  vsp->validateIDSet = vfp->validateIDSet;
  vsp->seqSubmitParent = vfp->seqSubmitParent;
  vsp->ignoreExceptions = vfp->ignoreExceptions;
  vsp->validateExons = vfp->validateExons;
  vsp->inferenceAccnCheck = vfp->inferenceAccnCheck;
  vsp->testLatLonSubregion = vfp->testLatLonSubregion;
  vsp->strictLatLonCountry = vfp->strictLatLonCountry;
  vsp->indexerVersion = vfp->indexerVersion;

  if (ofp == NULL && vfp->outfp != NULL) {
    ofp = vfp->outfp;
  }
  if (ofp != NULL) {
    vcd.ofp = ofp;
    vcd.verbosity = vfp->verbosity;
    vcd.lowCutoff = vfp->lowCutoff;
    vcd.highCutoff = vfp->highCutoff;
    vcd.errcode = vfp->errcode;
    vcd.vfp = vfp;
    vsp->errfunc = ValidCallback;
    vsp->userdata = (Pointer) &vcd;
    vsp->convertGiToAccn = FALSE;
  }

  ValidateSeqEntry (sep, vsp);

  for (i = 0; i <= 4; i++) {
    vfp->num_errors += vsp->errors [i];
    if (i >= vfp->severity) {
      vfp->fatal_errors += vsp->errors [i];
    }
  }

  ValidStructFree (vsp);
  if (vfp->validateBarcode) {
    if (vfp->verbosity == 4 && !vfp->has_errors) {
      cutoff = (ErrSev) vfp->lowCutoff;
      if (cutoff < SEV_NONE || cutoff > SEV_MAX) {
        cutoff = SEV_MAX;
      }
      xml_header = GetXmlHeaderText(cutoff);
    }
    if (!BarcodeValidateOneSeqEntry (ofp, sep, FALSE,
                                     vfp->verbosity == 4,
                                     !vfp->has_errors,
                                     xml_header)) {
      vfp->has_errors = TRUE;
    }
    xml_header = MemFree (xml_header);
  }
}

static void ProcessSingleRecord (
  CharPtr filename,
  ValFlagPtr vfp
)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  ValNodePtr     bsplist;
  BioseqSetPtr   bssp;
  Char           buf [64], path [PATH_MAX];
  Pointer        dataptr = NULL;
  Uint2          datatype, entityID = 0;
  FILE           *fp, *ofp = NULL;
  SeqEntryPtr    fsep, sep;
  ObjMgrPtr      omp;
  CharPtr        ptr;
  time_t         starttime, stoptime;

  if (StringHasNoText (filename)) return;
  if (vfp == NULL) return;

  if (vfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", path);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (vfp->type >= 2 && vfp->type <= 5) {
    aip = AsnIoOpen (filename, vfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    switch (vfp->type) {
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
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) vfp->type);
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

      starttime = GetSecs ();
      buf [0] = '\0';

      if (vfp->logfp != NULL) {
        fsep = FindNthBioseq (sep, 1);
        if (fsep != NULL && fsep->choice == 1) {
          bsp = (BioseqPtr) fsep->data.ptrvalue;
          if (bsp != NULL) {
            SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
            fprintf (vfp->logfp, "%s\n", buf);
            fflush (vfp->logfp);
          }
        }
      }

      StringNCpy_0 (path, filename, sizeof (path));
      ptr = StringRChr (path, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      StringCat (path, ".val");

      if (vfp->outpath != NULL) {
        ErrSetLogfile (vfp->outpath, ELOG_APPEND);
      } else if (vfp->verbosity == 0) {
        ErrSetLogfile (path, ELOG_APPEND);
      } else if (vfp->outfp == NULL) {
        ofp = FileOpen (path, "w");
      }

      bsplist = NULL;
    
      if (vfp->inferenceAccnCheck) {
        LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
      }
      if (vfp->lock) {
        bsplist = DoLockFarComponents (sep, vfp);
      }

      DoValidation (sep, vfp, ofp);

      bsplist = UnlockFarComponents (bsplist);

      if (ofp != NULL) {
        if (vfp->has_errors) {
          if (vfp->verbosity == 4) {
            fprintf (ofp, "</asnval>\n");
          }
          vfp->has_errors = FALSE;
        }
        FileClose (ofp);
      }

      stoptime = GetSecs ();
      if (stoptime - starttime > vfp->worsttime && StringDoesHaveText (buf)) {
        vfp->worsttime = stoptime - starttime;
        StringCpy (vfp->longest, buf);
      }
      (vfp->numrecords)++;
    }
  } else {
    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  ObjMgrFree (datatype, dataptr);

  omp = ObjMgrGet ();
  ObjMgrReapOne (omp);
  SeqMgrClearBioseqIndex ();
  ObjMgrFreeCache (0);
  FreeSeqIdGiCache ();

  SeqEntrySetScope (NULL);
}

static void ProcessMultipleRecord (
  CharPtr filename,
  ValFlagPtr vfp
)

{
  AsnIoPtr        aip;
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_bss, atp_desc, atp_sbp, atp_se, atp_ssp;
  BioseqPtr       bsp;
  ValNodePtr      bsplist;
  BioseqSetPtr    bssp;
  Char            buf [64], path [PATH_MAX], longest [64];
  Int2            skipcount = 0, maxcount = 0;
  CitSubPtr       csp = NULL;
  FILE            *fp, *ofp = NULL;
  Int4            numrecords = 0;
  SeqEntryPtr     fsep, sep;
  ObjMgrPtr       omp;
  ObjValNode      ovn;
  Pubdesc         pd;
  CharPtr         ptr;
  SubmitBlockPtr  sbp = NULL;
  time_t          starttime, stoptime, worsttime;
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
  if (vfp == NULL) return;

#ifndef OS_UNIX
  if (vfp->compressed) {
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

  if (vfp->compressed) {
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
    fp = popen (cmmd, /* vfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, vfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, vfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (vfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (vfp->type == 4) {
    atp = atp_bss;
  } else if (vfp->type == 5) {
    atp = atp_ssp;
  } else {
    Message (MSG_POSTERR, "Batch processing type not set properly");
    return;
  }

  longest [0] = '\0';
  worsttime = 0;

  StringNCpy_0 (path, filename, sizeof (path));
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (path, ".val");

  if (vfp->outpath != NULL) {
    ErrSetLogfile (vfp->outpath, ELOG_APPEND);
  } else if (vfp->verbosity == 0) {
    ErrSetLogfile (path, ELOG_APPEND);
  } else if (vfp->outfp == NULL) {
    ofp = FileOpen (path, "w");
  }

  while ((! vfp->io_failure) && maxcount < vfp->maxcount &&
         (atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (aip->io_failure) {
      vfp->io_failure = TRUE;
      aip->io_failure = FALSE;
    }
    if (atp == atp_se) {
      sep = SeqEntryAsnRead (aip, atp);
    
      /* propagate submission citation as descriptor onto each Seq-entry */

      if (subcit != NULL && sep != NULL && sep->data.ptrvalue != NULL) {
        if (sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        } else if (sep->choice == 2) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
      }

      if (sep != NULL) {
        if (skipcount < vfp->skipcount) {
          skipcount++;
        } else {

          starttime = GetSecs ();
          buf [0] = '\0';

          if (vfp->logfp != NULL) {
            fsep = FindNthBioseq (sep, 1);
            if (fsep != NULL && fsep->choice == 1) {
              bsp = (BioseqPtr) fsep->data.ptrvalue;
              if (bsp != NULL) {
                SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
                fprintf (vfp->logfp, "%s\n", buf);
                fflush (vfp->logfp);
              }
            }
          }

          bsplist = NULL;
    
          if (vfp->inferenceAccnCheck) {
            LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
          }
          if (vfp->lock) {
            bsplist = DoLockFarComponents (sep, vfp);
          }

          DoValidation (sep, vfp, ofp);

          bsplist = UnlockFarComponents (bsplist);

          stoptime = GetSecs ();
          if (stoptime - starttime > worsttime && StringDoesHaveText (buf)) {
            worsttime = stoptime - starttime;
            StringCpy (longest, buf);
          }
          numrecords++;
          maxcount++;
        }
      }

      SeqEntryFree (sep);
      omp = ObjMgrGet ();
      ObjMgrReapOne (omp);
      SeqMgrClearBioseqIndex ();
      ObjMgrFreeCache (0);
      FreeSeqIdGiCache ();

      SeqEntrySetScope (NULL);

#ifdef OS_UNIX
      if (detailed_report && vfp->logfp != NULL) {
        ObjMgrReportProc (vfp->logfp);
      }

      if (memory_usage && vfp->logfp != NULL) {
        Char     mbuf [512];
        FILE     *mufp;
        Char     ch;
        Int4     len;
        CharPtr  ptr1, ptr2;
        Int2     spaces;
        uid_t    uid;
        uid = getpid ();
        sprintf (cmmd, "cat /proc/%d/stat", (int) uid);
        mufp = popen (cmmd, "r");
        if (mufp != NULL) {
          len = FileRead ((Pointer) mbuf, sizeof (Char), sizeof (mbuf), mufp);
          if (len > 0) {
            mbuf [(int) len] = '\0';
            ptr1 = mbuf;
            ch = *ptr1;
            spaces = 0;
            while (ch != '\0' && spaces < 22) {
              if (ch == ' ') {
                spaces++;
              }
              ptr1++;
              ch = *ptr1;
            }
            if (ch != '\0') {
              ptr2 = StringChr (ptr1, ' ');
              if (ptr2 != NULL) {
                *ptr2 = '\0';
                fprintf (vfp->logfp, "Memory usage %s\n", ptr1);
              }
            }
          }
          pclose (mufp);
        }
      }
#endif

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

    if (aip->io_failure) {
      vfp->io_failure = TRUE;
      aip->io_failure = FALSE;
    }
  }

  if (aip->io_failure) {
    vfp->io_failure = TRUE;
  }

  if (vfp->io_failure) {
    Message (MSG_POSTERR, "Asn io_failure for input file '%s'", filename);
  }

  if (ofp != NULL) {
    if (vfp->has_errors) {
      if (vfp->verbosity == 4) {
        fprintf (ofp, "</asnval>\n");
      }
      vfp->has_errors = FALSE;
    }
    FileClose (ofp);
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

  if (vfp->logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (vfp->logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (vfp->logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (vfp->logfp);
  }
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  ValFlagPtr  vfp;

  vfp = (ValFlagPtr) userdata;
  if (vfp == NULL) return;

  if (vfp->logfp != NULL) {
    fprintf (vfp->logfp, "%s\n", filename);
    fflush (vfp->logfp);
  }

  if (vfp->batch) {
    ProcessMultipleRecord (filename, vfp);
  } else {
    ProcessSingleRecord (filename, vfp);
  }
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define i_argInputFile     1
#define o_argOutputFile    2
#define x_argSuffix        3
#define u_argRecurse       4
#define R_argSeverity      5
#define Q_argLowCutoff     6
#define P_argHighCutoff    7
#define E_argOnlyThisErr   8
#define A_argAlignments    9
#define J_argIsoJta       10
#define Z_argRemoteCDS    11
#define X_argExonSplice   12
#define G_argInfAccns     13
#define N_argLatLonStrict 14
#define M_argMatchTag     15
#define Y_argCheckOld     16
#define e_argIgnoreExcept 17
#define v_argVerbosity    18
#define a_argType         19
#define b_argBinary       20
#define c_argCompressed   21
#define r_argRemote       22
#define k_argLocalFetch   23
#define d_argAsnIdx       24
#define l_argLockFar      25
#define T_argThreads      26
#define L_argLogFile      27
#define S_argSkipCount    28
#define B_argBarcodeVal   29
#define C_argMaxCount     30
#ifdef INTERNAL_NCBI_ASN2VAL
#define w_argSeqSubParent 31
#define H_argAccessHUP    32
#define y_argAIndexer     33
#endif

#define LAT_LON_STATE    1
#define LAT_LON_STRICT   2

Args myargs [] = {
  {"Path to ASN.1 Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"File Selection Substring", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Severity for Error in Return Code", "4", "0", "4",
    FALSE, 'R', ARG_INT, 0.0, 0, NULL},
  {"Lowest Severity for Error to Show", "3", "0", "4",
    FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
  {"Highest Severity for Error to Show", "4", "0", "4",
    FALSE, 'P', ARG_INT, 0.0, 0, NULL},
  {"Only Error Code to Show", NULL, NULL, NULL,
    TRUE, 'E', ARG_STRING, 0.0, 0, NULL},
  {"Validate Alignments", "F", NULL, NULL,
    TRUE, 'A', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Require ISO-JTA?", "F", NULL, NULL,
    TRUE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote CDS Product Fetch", "F", NULL, NULL,
    TRUE, 'Z', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Exon Splice Check", "F", NULL, NULL,
    TRUE, 'X', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verify Inference Accessions", "F", NULL, NULL,
    TRUE, 'G', ARG_BOOLEAN, 0.0, 0, NULL},
  {"LatLon/Country Flags (1 Test State/Province, 2 Ignore Water Exception)", "0", "0", "3",
    TRUE, 'N', ARG_INT, 0.0, 0, NULL},
  {"Match locus_tag against General ID", "F", NULL, NULL,
    TRUE, 'M', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Check Against Old IDs", "F", NULL, NULL,
    TRUE, 'Y', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Ignore Transcription/Translation Exceptions", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbosity", "0", "0", "4",
    FALSE, 'v', ARG_INT, 0.0, 0, NULL},
  {"ASN.1 Type (a Any, e Seq-entry, b Bioseq, s Bioseq-set, m Seq-submit, t Batch Bioseq-set, u Batch Seq-submit)", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Batch File is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch File is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Local Fetching", "F", NULL, NULL,
    TRUE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Path to Indexed Binary ASN.1 Data", NULL, NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Skip Count", "0", NULL, NULL,
    TRUE, 'S', ARG_INT, 0.0, 0, NULL},
  {"Barcode Validate", "F", NULL, NULL,
    TRUE, 'B', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Max Count", "0", NULL, NULL,
    TRUE, 'C', ARG_INT, 0.0, 0, NULL},
#ifdef INTERNAL_NCBI_ASN2VAL
  {"SeqSubmitParent Flag", "F", NULL, NULL,
    TRUE, 'w', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Internal Access to HUP", "F", NULL, NULL,
    TRUE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Special Indexer Tests", "F", NULL, NULL,
    TRUE, 'y', ARG_BOOLEAN, 0.0, 0, NULL},
#endif
};

Int2 Main (void)

{
  Char         app [64];
  CharPtr      asnidx, directory, infile, logfile, outfile, str, suffix;
  Boolean      batch, binary, compressed, dorecurse,
               indexed, local, lock, remote, usethreads;
#ifdef INTERNAL_NCBI_ASN2VAL
  Boolean      hup = FALSE;
#endif
  time_t       run_time, start_time, stop_time;
  Int2         type = 0, val;
  ValFlagData  vfd;

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

  /* process command line arguments */

  sprintf (app, "asnval %s", ASNVAL_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* additional setup modifications */

  MemSet ((Pointer) &vfd, 0, sizeof (ValFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  remote = (Boolean ) myargs [r_argRemote].intvalue;
  local = (Boolean) myargs [k_argLocalFetch].intvalue;
#ifdef INTERNAL_NCBI_ASN2VAL
  hup = (Boolean) myargs [H_argAccessHUP].intvalue;
#endif
  asnidx = (CharPtr) myargs [d_argAsnIdx].strvalue;
  indexed = (Boolean) StringDoesHaveText (asnidx);
  lock = (Boolean) myargs [l_argLockFar].intvalue;
  usethreads = (Boolean) myargs [T_argThreads].intvalue;

  vfd.severity = (Int2) myargs [R_argSeverity].intvalue;
  vfd.lowCutoff = (Int2) myargs [Q_argLowCutoff].intvalue;
  vfd.highCutoff = (Int2) myargs [P_argHighCutoff].intvalue;
  vfd.errcode = (CharPtr) myargs [E_argOnlyThisErr].strvalue;
  vfd.validateAlignments = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.alignFindRemoteBsp = (Boolean) (vfd.validateAlignments && remote);
  vfd.doSeqHistAssembly = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.farIDsInAlignments = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.alwaysRequireIsoJTA = (Boolean) myargs [J_argIsoJta].intvalue;
  vfd.farFetchCDSproducts = (Boolean) myargs [Z_argRemoteCDS].intvalue;
  vfd.farFetchMRNAproducts = (Boolean) myargs [Z_argRemoteCDS].intvalue;
  vfd.locusTagGeneralMatch = (Boolean) myargs [M_argMatchTag].intvalue;
  vfd.validateIDSet = (Boolean) myargs [Y_argCheckOld].intvalue;
  vfd.ignoreExceptions = (Boolean) myargs [e_argIgnoreExcept].intvalue;
  vfd.validateExons = (Boolean) myargs [X_argExonSplice].intvalue;
  vfd.inferenceAccnCheck = (Boolean) myargs [G_argInfAccns].intvalue;
  vfd.validateBarcode = (Boolean) myargs[B_argBarcodeVal].intvalue;


  val = (Int2) myargs [N_argLatLonStrict].intvalue;
  vfd.testLatLonSubregion = (Boolean) ((val & LAT_LON_STATE) != 0);
  vfd.strictLatLonCountry = (Boolean) ((val & LAT_LON_STRICT) != 0);

  vfd.verbosity = (Int2) myargs [v_argVerbosity].intvalue;

  vfd.skipcount = (Int4) myargs [S_argSkipCount].intvalue;
  vfd.maxcount = (Int4) myargs [C_argMaxCount].intvalue;
  if (vfd.maxcount < 1) {
    vfd.maxcount = INT4_MAX;
  }

#ifdef INTERNAL_NCBI_ASN2VAL
  vfd.seqSubmitParent = (Boolean) myargs [w_argSeqSubParent].intvalue;
  vfd.indexerVersion = (Boolean) myargs [y_argAIndexer].intvalue;
#endif

#ifdef INTERNAL_NCBI_ASN2VAL
  SetAppProperty ("InternalNcbiSequin", (void *) 1024);
#endif

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

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;

  start_time = GetSecs ();

  /* populate parameter structure */

  vfd.batch = batch;
  vfd.binary = binary;
  vfd.compressed = compressed;
  vfd.lock = lock;
  vfd.useThreads = usethreads;
  vfd.type = type;
  vfd.logfp = NULL;
  vfd.num_errors = 0;
  vfd.fatal_errors = 0;
  vfd.has_errors = FALSE;
  vfd.io_failure = FALSE;
  vfd.longest [0] = '\0';
  vfd.worsttime = 0;
  vfd.numrecords = 0;

  if (! StringHasNoText (outfile)) {
    if (vfd.verbosity == 0) {
      vfd.outpath = outfile;
    } else {
      vfd.outfp = FileOpen (outfile, "w");
      if (vfd.outfp == NULL) {
        Message (MSG_FATAL, "Unable to open single output file");
        return 1;
      }
    }
  }

  if (! StringHasNoText (logfile)) {
    vfd.logfp = FileOpen (logfile, "w");
    if (vfd.logfp == NULL) {
      Message (MSG_FATAL, "Unable to open log file");
      return 1;
    }
  }

  /* register fetch functions */

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2VAL
    if (hup) {
      DirSubFetchEnable ();
      SmartFetchEnable ();
      TPASmartFetchEnable ();
    }

    if (! PUBSEQBioseqFetchEnable ("asnval", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    vfd.usePUBSEQ = TRUE;
    vfd.useThreads = FALSE;
#else
    PubSeqFetchEnable ();
#endif
    if (vfd.inferenceAccnCheck) {
      SeqMgrSetPreCache (GiRevHistLookupFarSeqIDs);
    }
    if (vfd.validateIDSet) {
      SeqMgrSetSeqIdSetFunc (GiRevHistLookupSeqIdSet);
    }
  }

  if (local) {
    LocalSeqFetchInit (FALSE);
  }

  if (indexed) {
    AsnIndexedLibFetchEnable (asnidx, TRUE);
  }

  /* recurse through all files within source directory or subdirectories */

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, dorecurse, ProcessOneRecord, (Pointer) &vfd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, (Pointer) &vfd);
  }

  stop_time = GetSecs ();
  run_time = stop_time - start_time;

  if (vfd.outfp != NULL) {
    if (vfd.has_errors) {
      if (vfd.verbosity == 4) {
        fprintf (vfd.outfp, "</asnval>\n");
      }
      vfd.has_errors = FALSE;
    }
    FileClose (vfd.outfp);
  }

  if (vfd.logfp != NULL) {
    fprintf (vfd.logfp, "Finished in %ld seconds\n", (long) run_time);
    if (StringDoesHaveText (vfd.longest)) {
      fprintf (vfd.logfp, "Longest processing time %ld seconds on %s\n",
               (long) vfd.worsttime, vfd.longest);
      fprintf (vfd.logfp, "Total number of records %ld\n", (long) vfd.numrecords);
    }
    FileClose (vfd.logfp);
  }

  /* close fetch functions */

  if (indexed) {
    AsnIndexedLibFetchDisable ();
  }

  if (local) {
    LocalSeqFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2VAL
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
    SeqMgrSetPreCache (NULL);
    SeqMgrSetSeqIdSetFunc (NULL);
  }

  TransTableFreeAll ();

  ECNumberFSAFreeAll ();

  if (vfd.fatal_errors > 0) return 1;
  if (vfd.io_failure) return 1;

  return 0;
}

