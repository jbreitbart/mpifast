/*   raw2delt.c
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
* File Name:  raw2delt.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   4/12/07
*
* $Revision: 1.2 $
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
#include <seqport.h>

#define RAW2DELT_APP_VER "1.0"

CharPtr RAW2DELT_APPLICATION = RAW2DELT_APP_VER;

typedef struct outputstream {
  CharPtr  results_dir;
  CharPtr  base;
  CharPtr  suffix;
  CharPtr  outfile;
  CharPtr  outsuffix;
  AsnIoPtr aip;
  Boolean  is_binary;
} OutputStreamData, PNTR OutputStreamPtr;

typedef struct inputstream {
  CharPtr directory;
  CharPtr base;
  CharPtr suffix;
  Boolean is_binary;
  Boolean is_seqentry;
} InputStreamData, PNTR InputStreamPtr;

typedef struct asnstream {
  AsnModulePtr amp;
  AsnTypePtr   atp_se;
  AsnTypePtr   atp_bss;
  AsnTypePtr   atp_bss_se;
} AsnStreamData, PNTR AsnStreamPtr;

static FILE* OpenOneFile (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix
)

{
  Char  file [FILENAME_MAX], path [PATH_MAX];

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  return FileOpen (path, "r");
}

static AsnIoPtr AsnIoFromInputStream (
  InputStreamPtr isp
)

{
  AsnIoPtr aip;
  Char     file [FILENAME_MAX], path [PATH_MAX];
  CharPtr  read_flag;

  if (isp == NULL) return NULL;

  if (isp->is_binary) {
    read_flag = "rb";
  } else {
    read_flag = "r";
  }

  if (isp->base == NULL) {
    aip = AsnIoOpen ("stdin", read_flag);
  } else {
    StringNCpy_0 (path, isp->directory, sizeof (path));
    sprintf (file, "%s%s", isp->base, isp->suffix);
    FileBuildPath (path, NULL, file);
    aip = AsnIoOpen (path, read_flag);
  }
  return aip;
}


static AsnIoPtr AsnIoFromOutputStream (OutputStreamPtr osp)
{
  AsnIoPtr   aip;
  Char       file [FILENAME_MAX], path [PATH_MAX];
  CharPtr    write_flag;

  if (osp == NULL) return NULL;
  if (osp->aip == NULL) {
    write_flag = osp->is_binary ? "wb" : "w";
    if (StringDoesHaveText (osp->outfile)) {
      StringNCpy_0 (path, osp->outfile, sizeof (path));
    } else {
      if (osp->base == NULL) {
        aip = AsnIoOpen ("stdout", write_flag);
      } else {
        if (osp->outsuffix == NULL) {
          osp->outsuffix = "";
        }
        StringNCpy_0 (path, osp->results_dir, sizeof (path));
        sprintf (file, "%s%s%s", osp->base, osp->suffix, osp->outsuffix);
        FileBuildPath (path, NULL, file);
        aip = AsnIoOpen (path, write_flag);
        if (aip == NULL) {
          Message (MSG_POSTERR, "Unable to write to %s.", path);
        }
      }
    }
  } else {
    aip = osp->aip;
  }
  return aip;
}

static void WriteOneFile (
  OutputStreamPtr osp,
  SeqEntryPtr sep
)

{
  AsnIoPtr   aip;

  aip = AsnIoFromOutputStream (osp);
  if (aip != NULL) {
    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoFlush (aip);
  }
  if (aip != osp->aip) {
    AsnIoClose (aip);
  }
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


static void ProcessSeqEntry (SeqEntryPtr sep, Int4Ptr gap_sizes)
{
  ValNodePtr     bsp_list = NULL, vnp;
  BioseqPtr      bsp;

  if (sep == NULL || gap_sizes == NULL) return;

  VisitBioseqsInSep (sep, &bsp_list, CollectBioseqsForConversion);

  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    ConvertNsToGaps (bsp, gap_sizes);
  }
  bsp_list = ValNodeFree (bsp_list);
}

static Uint2 ProcessOneAsn (
  FILE* fp,
  CharPtr path,
  Int4Ptr gap_sizes
)

{
  Pointer        dataptr;
  Uint2          datatype, entityID = 0;
  SeqEntryPtr    sep;

  if (fp == NULL || gap_sizes == NULL) return 0;

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, TRUE, FALSE, TRUE, FALSE);
  if (dataptr == NULL) {
    Message (MSG_POSTERR, "Unable to read data from %s.", path);
    return 0;
  }

  sep = GetTopSeqEntryForEntityID (entityID);
  ProcessSeqEntry (sep, gap_sizes);

  return entityID;
}

/* return -1 if failure, 0 if success */
static Int4 ProcessOneRecord (
  CharPtr directory,
  OutputStreamPtr osp,
  Int4Ptr   gap_sizes
)

{
  Uint2              entityID;
  FILE               *fp;
  SeqEntryPtr        sep;

  if (osp == NULL || gap_sizes == NULL) return -1;
  fp = OpenOneFile (directory, osp->base, osp->suffix);
  if (fp == NULL) return -1;

  entityID = ProcessOneAsn (fp, osp->base == NULL ? "input stream" : osp->base, gap_sizes);
  
  FileClose (fp);

  if (entityID == 0) return -1;

  /* finish processing */

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL) {
    WriteOneFile (osp, sep);
  }

  ObjMgrFreeByEntityID (entityID);
  return 0;
}

static Int4 ProcessStream (InputStreamPtr isp, OutputStreamPtr osp, AsnStreamPtr asp, Int4Ptr gap_sizes)
{
  AsnTypePtr   atp, atp_srch;
  AsnIoPtr asn_in, asn_out;
  Int4         rval = 0;
  SeqEntryPtr  sep;
  Uint2        entityID;
  DataVal av;
  
  if (isp == NULL || osp == NULL || asp == NULL || gap_sizes == NULL) return 1;

  asn_in = AsnIoFromInputStream (isp);
  asn_out = AsnIoFromOutputStream (osp);

  if (isp->is_seqentry) {
    atp = asp->atp_se;
    atp_srch = asp->atp_se;
  }
  else {
    atp = asp->atp_bss;
    atp_srch = asp->atp_bss_se;
  }

  while ((atp = AsnReadId(asn_in, asp->amp, atp)) != NULL && rval == 0) {
    if (atp != atp_srch) {
      AsnReadVal(asn_in, atp, &av);
      AsnWrite(asn_out, atp, &av);
      AsnKillValue(atp, &av);
      continue;
    }
    if ((sep = SeqEntryAsnRead(asn_in, atp)) == NULL) {
      Message (MSG_POSTERR, "SeqEntryAsnRead failure");
      rval = 1;
    } else  {
      entityID = ObjMgrRegister (OBJ_SEQENTRY, sep);
      ProcessSeqEntry (sep, gap_sizes);

      if (! SeqEntryAsnWrite(sep, asn_out, atp)) {
       Message (MSG_POSTERR, "SeqEntryAsnWrite failure");
       rval = 1;
      }
      AsnIoFlush(asn_out);
      ObjMgrFreeByEntityID (entityID);
    }
  }                             /* Endwhile, AsnReadId */

  AsnIoClose(asn_in);
  if (asn_out != osp->aip) {
    AsnIoClose(asn_out);
  }
  
  return rval;
}

/* return -1 on failure, 0 on success */
static Int4 FileRecurse (
  CharPtr         directory,
  InputStreamPtr  isp,
  OutputStreamPtr osp,
  AsnStreamPtr    asp,
  Int4Ptr         gap_sizes
)

{
  Char        path [PATH_MAX];
  CharPtr     ptr;
  CharPtr     str;
  ValNodePtr  head, vnp;
  CharPtr     orig_dir, orig_base;
  Int4        rval = 0;

  /* get list of all files in source directory */

  head = DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringDoesHaveText (str)) {

        /* does filename have desired substring? */

        ptr = StringStr (str, osp->suffix);

        if (ptr != NULL) {

          /* make sure detected suffix is really at end of filename */

          if (StringCmp (ptr, osp->suffix) == 0) {
            *ptr = '\0';

            /* process file that has desired suffix (usually .fsa) */
            osp->base = str;
            orig_dir = isp->directory;
            isp->directory = directory;
            orig_base = isp->base;
            isp->base = str;
            if (isp->is_binary) {
              rval |= ProcessStream (isp, osp, asp, gap_sizes);
            } else {
              rval |= ProcessOneRecord (directory, osp, gap_sizes);
            }
            isp->directory = orig_dir;
            isp->base = orig_base;
            osp->base = NULL;
          }
        }
      }
    } else if (vnp->choice == 1) {

      /* recurse into subdirectory */

      StringNCpy_0 (path, directory, sizeof (path));
      str = (CharPtr) vnp->data.ptrvalue;
      FileBuildPath (path, str, NULL);
      rval |= FileRecurse (path, isp, osp, asp, gap_sizes);
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);
  return rval;
}

static Boolean SetUpAsnStreamData (AsnStreamPtr asp)

{
  if (asp == NULL) return FALSE;

  if (! SeqSetAsnLoad()) {
    Message (MSG_POSTERR, "Unable to load SeqSet parse tree");
    return FALSE;
  }
  asp->amp = AsnAllModPtr();
  if (asp->amp == NULL) {
    Message (MSG_POSTERR, "Unable to obtain ASN.1 module pointer");
    return FALSE;
  }

  /* Get pointers to ASN.1 types that must be dealt with in asn_in */

  if ( (asp->atp_bss = AsnFind("Bioseq-set")) == NULL) {
    Message (MSG_POSTERR, "could not find type Bioseq-set");
    return FALSE;
  }
  if ( (asp->atp_bss_se = AsnFind("Bioseq-set.seq-set.E")) == NULL) {
    Message (MSG_POSTERR, "AsnFind failure: Bioseq-set.seq-set.E");
    return FALSE;
  }
  if ( (asp->atp_se = AsnFind("Seq-entry")) == NULL) {
    Message (MSG_POSTERR, "AsnFind failure: Seq-entry");
    return FALSE;
  }
  return TRUE;
}

/* Args structure contains command-line arguments */

#define p_argInputPath         0
#define r_argOutputPath        1
#define i_argInputFile         2
#define o_argOutputFile        3
#define x_argSuffix            4
#define s_argOutSuffix         5
#define b_argInputBinary       6
#define e_argInputSeqEntry     7
#define d_argOutputBinary      8
#define u_argEqUnknownGap      9
#define U_argGTEqUnknownGap   10
#define k_argEqUnknownGap     11
#define K_argGtEqUnknownGap   12


Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Suffix", ".sqn", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Suffix for converted files", "", NULL, NULL,
    TRUE, 's', ARG_STRING, 0.0, 0, NULL},
  {"Input is binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Input is Seq-entry", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Output is binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Exact size to convert to unknown length gap", "-1", NULL, NULL,
    TRUE, 'u', ARG_INT, 0.0, 0, NULL},
  {"Size greater than/equal to to convert to unknown length gap", "-1", NULL, NULL,
    TRUE, 'U', ARG_INT, 0.0, 0, NULL},
  {"Exact size to convert to known length gap", "-1", NULL, NULL,
    TRUE, 'k', ARG_INT, 0.0, 0, NULL},
  {"Size greater than/equal to to convert to known length gap", "-1", NULL, NULL,
    TRUE, 'K', ARG_INT, 0.0, 0, NULL},
};

Int2 Main(void)
{
  Char             app [64];
  CharPtr          directory;
  CharPtr          ptr;
  Char             sfx [32];
  OutputStreamData osd;
  InputStreamData  isd;
  AsnStreamData    asd;
  Int4             gap_sizes[2];
  Int4             rval = 0;
  Int4             u_eq = 0, u_gteq = -1, k_eq = 0, k_gteq = -1;

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

  SetUpAsnStreamData (&asd);

  /* initialize OuputStreamData */
  MemSet (&osd, 0, sizeof (osd));

  /* initialize InputStreamData */
  MemSet (&isd, 0, sizeof (isd));

  /* initialize gap_sizes */
  gap_sizes[0] = 0;
  gap_sizes[1] = 0;

  /* process command line arguments */

  sprintf (app, "raw2delt %s", RAW2DELT_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  u_eq = (Int4) myargs [u_argEqUnknownGap].intvalue;
  u_gteq = (Int4) myargs [U_argGTEqUnknownGap].intvalue;
  k_eq = (Int4) myargs [k_argEqUnknownGap].intvalue;
  k_gteq = (Int4) myargs [K_argGtEqUnknownGap].intvalue;

  if (u_eq < 1 && u_gteq < 1 && k_eq < 1 && k_gteq < 1) {
    Message (MSG_FATAL, "Must specify values for at least one of -u, -U, -k, -K");
    return 1;
  } else if (u_eq > -1 && u_gteq > -1) {
    Message (MSG_FATAL, "May only specify value for -u or -U, not both");
    return 1;
  } else if (k_eq > -1 && k_gteq > -1) {
    Message (MSG_FATAL, "May only specify value for -k or -K, not both");
    return 1;
  }
    
  if (u_eq > 0) {
    gap_sizes[0] = u_eq;
  } else if (u_gteq > 0) {
    gap_sizes[0] = 0 - u_gteq;
  }

  if (k_eq > 0) {
    gap_sizes[1] = k_eq;
  } else if (k_gteq > 0) {
    gap_sizes[1] = 0 - k_gteq;
  }

  if (gap_sizes[0] == gap_sizes[1]) {
    Message (MSG_FATAL, "Cannot specify the same size for known and unknown length gaps");
    return 1;
  }

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  osd.results_dir = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (osd.results_dir)) {
    osd.results_dir = NULL;
  }
  osd.suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  osd.outsuffix = (CharPtr) myargs [s_argOutSuffix].strvalue;
  osd.base = (CharPtr) myargs [i_argInputFile].strvalue;
  osd.outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  if (StringHasNoText (osd.outfile)) {
    osd.outfile = NULL;
  }
  osd.is_binary = (Boolean) myargs [d_argOutputBinary].intvalue;

  if (osd.base == "stdin") {
    osd.base = NULL;
  }

  /* if we don't have an output directory or an output file, and the user hasn't provided an
   * output suffix, add a default.
   */
  if (osd.results_dir == NULL && osd.outfile == NULL && StringHasNoText (osd.outsuffix)) {
    osd.outsuffix = ".delta";
  } 

  isd.is_binary = (Boolean) myargs [b_argInputBinary].intvalue;
  isd.is_seqentry = (Boolean) myargs [e_argInputSeqEntry].intvalue;
  isd.directory = directory;
  isd.base = osd.base;
  isd.suffix = osd.suffix;

  if (StringDoesHaveText (osd.outfile)) {
    osd.aip = AsnIoOpen (osd.outfile, "w");
    if (osd.aip == NULL) {
      Message (MSG_FATAL, "Unable to open output file");
      return 1;
    }
  } else {
    if (StringHasNoText (osd.results_dir)) {
      osd.results_dir = directory;
    }
    /* if we're putting the results in a separate directory, strip the directory name from the output base */
    if (!StringHasNoText (osd.results_dir) && !StringHasNoText (osd.base)) {
#ifdef OS_MSWIN
      ptr = StringRChr (osd.base, '\\');
#else
      ptr = StringRChr (osd.base, '/');
#endif;
      if (ptr != NULL) {
        osd.base = ptr + 1;
      }
    }
  }


  if (StringHasNoText(directory) && StringHasNoText(osd.base)) {
    rval = ProcessStream (&isd, &osd, &asd, gap_sizes);
  } else if (StringDoesHaveText (osd.base)) {
    ptr = StringRChr (osd.base, '.');
    sfx[0] = '\0';
    if (ptr != NULL) {
      StringNCpy_0 (sfx, ptr, sizeof (sfx));
      *ptr = '\0';
    }
    osd.suffix = sfx;
    isd.suffix = sfx;
    if (isd.is_binary) {
      rval = ProcessStream (&isd, &osd, &asd, gap_sizes);
    } else {
      rval = ProcessOneRecord (directory, &osd, gap_sizes);
    }
  } else {

    rval = FileRecurse (directory, &isd, &osd, &asd, gap_sizes);
  }

  if (osd.aip != NULL) {
    AsnIoFlush (osd.aip);
    AsnIoClose (osd.aip);
  }
  return rval;
}
