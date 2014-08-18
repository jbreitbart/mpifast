/*   asn2idx.c
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
* File Name:  asn2idx.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   8/2/04
*
* $Revision: 1.5 $
*
* File Description:
*
*   Generates accession/file offset indices for Bioseq-set release files
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
#include <sqnutils.h>
#include <lsqfetch.h>

static Boolean NotInFilter (
  CharPtr filename,
  ValNodePtr exclude
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  for (vnp = exclude; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringStr (filename, str) != NULL) return FALSE;
  }
  return TRUE;
}

static void FileRecurse (
  CharPtr directory,
  CharPtr results,
  CharPtr subdir,
  CharPtr subfile,
  ValNodePtr exclude,
  Boolean binary,
  Boolean dorecurse
)

{
  Char        path [PATH_MAX];
  CharPtr     str;
  ValNodePtr  head, vnp;

  /* get list of all files in source directory */

  head = DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      if (StringHasNoText (subdir) || StringStr (directory, subdir) != NULL) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (! StringHasNoText (str)) {

          /* does filename have desired substring? */

          if (StringHasNoText (subfile) || StringStr (str, subfile) != NULL) {

            if (NotInFilter (str, exclude)) {

              /* process file that has desired suffix (usually .aso) */

              StringNCpy_0 (path, directory, sizeof (path));
              FileBuildPath (path, NULL, str);

              CreateAsnIndex (path, results, binary);
            }
          }
        }
      }
    } else if (vnp->choice == 1 && dorecurse) {
      StringNCpy_0 (path, directory, sizeof (path));
      str = (CharPtr) vnp->data.ptrvalue;
      FileBuildPath (path, str, NULL);
      FileRecurse (path, results, subdir, subfile, exclude, binary, dorecurse);
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);
}

static ValNodePtr ParseStringByCommas (
  CharPtr str
)

{
  Char        ch;
  ValNodePtr  head = NULL;
  CharPtr     last, ptr, tmp;

  if (StringHasNoText (str)) return NULL;

  tmp = StringSave (str);
  last = tmp;
  ptr = last;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == ',') {
      *ptr = '\0';
      if (! StringHasNoText (last)) {
        TrimSpacesAroundString (last);
        ValNodeCopyStr (&head, 0, last);
      }
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
    ValNodeCopyStr (&head, 0, last);
  }
  MemFree (tmp);

  return head;
}

/* Args structure contains optional command-line arguments */

#define p_argInputPath    0
#define r_argOutputPath   1
#define d_argSubDirName   2
#define x_argFileSelect   3
#define f_argFilter       4
#define b_argBinary       5
#define t_argRecurse      6

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Required Subdirectory", NULL, NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Substring", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Filter", "gbcon,gbest,gbgss,gbhtg,gbsts", NULL, NULL,
    FALSE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-sets are Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 't', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Boolean     binary, dorecurse;
  CharPtr     directory, filter, results, subdir, subfile;
  ValNodePtr  exclude = NULL;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
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

  if (! GetArgs ("asn2idx", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* additional setup modifications */

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = NULL;
  }
  subdir = (CharPtr) myargs [d_argSubDirName].strvalue;
  subfile = (CharPtr) myargs [x_argFileSelect].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  binary = (Boolean) myargs [b_argBinary].intvalue;
  dorecurse = (Boolean) myargs [t_argRecurse].intvalue;

  if (StringHasNoText (directory)) {
    Message (MSG_FATAL, "Required path to files not specified");
    return 1;
  }

  exclude = ParseStringByCommas (filter);

  /* index recursively within specified directory */

  FileRecurse (directory, results, subdir, subfile, exclude, binary, dorecurse);

  /* now create master index combining individual indices */

  CreateMasterAsnIndex (directory, results);

  ValNodeFreeData (exclude);

  return 0;
}

