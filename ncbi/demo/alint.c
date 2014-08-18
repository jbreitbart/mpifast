/*   alint.c
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
* File Name:  alint.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/10/08
*
* $Revision: 1.1 $
*
* File Description:
*
*  Lint for Alignments in FASTA format - upper cases points of exact match
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#include <ncbi.h>
#include <sqnutils.h>

static CharPtr GetSequence (
  CharPtr str,
  Boolean skiptoken
)

{
  Char  ch;

  if (str == NULL) return NULL;

  if (! skiptoken) return str;

  ch = *str;
  while (ch != '\0' && ch != ' ') {
    str++;
    ch = *str;
  }
  if (ch == ' ') {
    str++;
  }

  return str;
}

static void ProcessAlignedFASTA (
  FILE *ifp,
  FILE *ofp,
  Boolean skiptoken
)

{
  CharPtr PNTR  array;
  Char          ch, ch0;
  FileCache     fc;
  ValNodePtr    head = NULL, last = NULL, vnp;
  Int2          i, j, num = 0, len, minlen = INT2_MAX, matches = 0, mismatches = 0;
  Char          line [4096];
  Boolean       match;
  CharPtr       ptr, str;

  FileCacheSetup (&fc, ifp);

  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  if (str == NULL) return;

  while (str != NULL) {
    TrimSpacesAroundString (str);
    if (StringDoesHaveText (str)) {
      vnp = ValNodeCopyStr (&last, 0, str);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
      num++;
      str = GetSequence (str, skiptoken);
      len = (Int2) StringLen (str);
      if (minlen > len) {
        minlen = len;
      }
    }
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  if (num < 1 || minlen < 1) return;

  array = (CharPtr PNTR) MemNew (sizeof (CharPtr) * (num + 1));
  if (array == NULL) return;

  for (vnp = head, i = 0; vnp != NULL; vnp = vnp->next, i++) {
    str = (CharPtr) vnp->data.ptrvalue;
    array [i] = str;
  }

  for (j = 0; j < minlen; j++) {
    ptr = GetSequence (array [0], skiptoken);
    ch0 = ptr [j];
    match = TRUE;

    for (i = 1; i < num; i++) {
      ptr = GetSequence (array [i], skiptoken);
      ch = ptr [j];
      if (ch != ch0) {
        match = FALSE;
      }
    }

    if (match) {
      matches++;
    } else {
      mismatches++;
    }

    for (i = 0; i < num; i++) {
      ptr = GetSequence (array [i], skiptoken);
      ch = ptr [j];
      if (match) {
        ptr [j] = TO_UPPER (ch);
      } else {
        ptr [j] = TO_LOWER (ch);
      }
    }
  }

  for (vnp = head, i = 0; vnp != NULL; vnp = vnp->next, i++) {
    str = (CharPtr) vnp->data.ptrvalue;
    fprintf (ofp, "%s\n", str);
  }

  fprintf (ofp, "\n%d matches, %d mismatches, length %d, %d percent matching\n",
           (int) matches, (int) mismatches, (int) minlen,
           (int) (matches * 100 / minlen));

  MemFree (array);
  ValNodeFreeData (head);
}

#define i_argInputFile    0
#define o_argOutputFile   1
#define s_argSkipToken    2

Args myargs [] = {
  {"Input File", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
 {"Skip First Token", "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  FILE     *ifp, *ofp;
  CharPtr  infile, outfile;
  Boolean  skiptoken;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrPathReset ();

  if (! GetArgs ("alint", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  skiptoken = (Boolean) myargs [s_argSkipToken].intvalue;

  ifp = FileOpen (infile, "r");
  if (ifp == NULL) {
    Message (MSG_FATAL, "Unable to open input file");
    return 1;
  }

  ofp = FileOpen (outfile, "w");
  if (ofp == NULL) {
    Message (MSG_FATAL, "Unable to open output file");
    return 1;
  }

  ProcessAlignedFASTA (ifp, ofp, skiptoken);

  FileClose (ofp);
  FileClose (ifp);

  return 0;
}

