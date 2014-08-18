/*   sugint.c
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
* File Name:  sugint.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   10/31/08
*
* $Revision: 1.1 $
*
* File Description:
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
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
#include <tofasta.h>
#include <gather.h>
#include <explore.h>
#include <suggslp.h>

static SeqEntryPtr ReadSep (
  FILE *fp,
  Boolean forceNuc,
  Boolean forcePrt
)

{
  Pointer  dataptr;
  Uint2    datatype, entityID;

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, forceNuc, forcePrt, TRUE, FALSE);
  if (dataptr == NULL) return NULL;
  entityID = ObjMgrRegister (datatype, dataptr);
  return GetTopSeqEntryForEntityID (entityID);
}

static void ProcessSuggest (
  FILE *nfp,
  FILE *pfp,
  AsnIoPtr ofp,
  Int2 gencode
)

{
  BioseqPtr    nbsp = NULL, pbsp = NULL;
  SeqEntryPtr  nsep, psep, sep;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;
  SeqLocPtr    slp;

  nsep = ReadSep (nfp, TRUE, FALSE);
  psep = ReadSep (pfp, FALSE, TRUE);

  if (nsep != NULL && psep != NULL) {
    sep = FindNthBioseq (nsep, 1);
    if (sep != NULL && IS_Bioseq (sep)) {
      nbsp = (BioseqPtr) sep->data.ptrvalue;
    }
    sep = FindNthBioseq (psep, 1);
    if (sep != NULL && IS_Bioseq (sep)) {
      pbsp = (BioseqPtr) sep->data.ptrvalue;
    }
    if (nbsp != NULL && pbsp != NULL) {
      if (ISA_na (nbsp->mol) && ISA_aa (pbsp->mol)) {
        sap = SuggestCodingRegion (nbsp, pbsp, gencode);

        if (sap != NULL && sap->type == 1) {
          sfp = (SeqFeatPtr) sap->data;
          if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
            slp = sfp->location;
            if (slp != NULL) {
              SeqLocAsnWrite (slp,  ofp, NULL);
            }
          }
        }

        SeqAnnotFree (sap);
      }
    }
  }

  SeqEntryFree (nsep);
  SeqEntryFree (psep);
}

#define n_argNucInputFile  0
#define p_argPrtInputFile  1
#define o_argOutputFile    2
#define g_argGeneticCode   3

Args myargs [] = {
  {"Nucleotide Input File", NULL, NULL, NULL,
    FALSE, 'n', ARG_FILE_IN, 0.0, 0, NULL},
  {"Protein Input File", NULL, NULL, NULL,
    FALSE, 'p', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", NULL, NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Genetic Code", "1", "0", "20",
    TRUE, 'g', ARG_INT, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Int2      gencode;
  FILE      *nfp, *pfp;
  AsnIoPtr  ofp;
  CharPtr   nucfile, prtfile, outfile;

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

  if (! GetArgs ("sugint", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  nucfile = (CharPtr) myargs [n_argNucInputFile].strvalue;
  prtfile = (CharPtr) myargs [p_argPrtInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  gencode = (Int2) myargs [g_argGeneticCode].intvalue;

  nfp = FileOpen (nucfile, "r");
  if (nfp == NULL) {
    Message (MSG_FATAL, "Unable to open nucleotide input file");
    return 1;
  }

  pfp = FileOpen (prtfile, "r");
  if (pfp == NULL) {
    Message (MSG_FATAL, "Unable to open protein input file");
    return 1;
  }

  ofp = AsnIoOpen (outfile, "w");
  if (ofp == NULL) {
    Message (MSG_FATAL, "Unable to open output file");
    return 1;
  }

  ProcessSuggest (nfp, pfp, ofp, gencode);

  AsnIoClose (ofp);
  FileClose (pfp);
  FileClose (nfp);

  return 0;
}

