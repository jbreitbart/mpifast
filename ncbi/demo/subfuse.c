/*   subfuse.c
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
* File Name:  subfuse.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/30/01
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

static SeqSubmitPtr ReadOneSubmission (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix
)

{
  AsnIoPtr      aip;
  Char          file [FILENAME_MAX], path [PATH_MAX];
  SeqSubmitPtr  ssp;

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }
  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  aip = AsnIoOpen (path, "r");
  if (aip == NULL) return NULL;
  ssp = SeqSubmitAsnRead (aip, NULL);
  AsnIoClose (aip);

  return ssp;
}

static void WriteOneSubmission (
  CharPtr path,
  SeqSubmitPtr ssp
)

{
  AsnIoPtr  aip;

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  SeqSubmitAsnWrite (ssp, aip, NULL);

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static void ProcessOneRecord (
  SeqSubmitPtr master,
  BioseqSetPtr bssp,
  CharPtr directory,
  CharPtr base,
  CharPtr suffix
)

{
  SeqEntryPtr   sep;
  SeqSubmitPtr  ssp;

  ssp = ReadOneSubmission (directory, base, suffix);
  if (ssp == NULL || ssp->datatype != 1) return;

  if (master->sub == NULL) {
    master->sub = ssp->sub;
    ssp->sub = NULL;
  }

  sep = (SeqEntryPtr) ssp->data;
  ssp->data = NULL;

  ValNodeLink (&(bssp->seq_set), sep);
}

/* Args structure contains command-line arguments */

#define p_argInputPath  0
#define o_argOutputFile 1
#define x_argSuffix     2

Args myargs [] = {
  {"Path to files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Output file", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Suffix", ".sqn", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  CharPtr       base, directory, outfile, suffix, ptr;
  BioseqSetPtr  bssp;
  ValNodePtr    head, vnp;
  SeqEntryPtr   sep;
  SeqSubmitPtr  ssp;

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

  if (! GetArgs ("subfuse", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  bssp = BioseqSetNew ();
  if (bssp == NULL) return 0;
  bssp->_class = BioseqseqSet_class_genbank;

  sep = SeqEntryNew ();
  if (sep == NULL) return 0;
  sep->choice = 2;
  sep->data.ptrvalue = (Pointer) bssp;

  ssp = SeqSubmitNew ();
  if (ssp == NULL) return 0;
  ssp->datatype = 1;
  ssp->data = (Pointer) sep;

  /* get list of all files in source directory */

  head = DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      base = (CharPtr) vnp->data.ptrvalue;
      if (! StringHasNoText (base)) {
        ptr = StringStr (base, suffix);
        if (ptr != NULL) {
          *ptr = '\0';
          Message (MSG_POST, "Processing %s\n", base);
          ProcessOneRecord (ssp, bssp, directory, base, suffix);
        }
      }
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);

  /* write output file */

  WriteOneSubmission (outfile, ssp);

  return 0;
}

