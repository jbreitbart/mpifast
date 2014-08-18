/*   vibmain.h
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
* File Name:  vibmain.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   12/13/06
*
* $Revision: 1.1 $
*
* File Description:
*       Vibrant main functions
*
* ==========================================================================
*/

#ifndef _VIBMAIN_
#define _VIBMAIN_

#ifndef _NCBI_
#include <ncbi.h>
#endif


#ifdef WIN_MAC
#ifdef OS_UNIX_DARWIN
extern int Nlm_VibMainPrelude (int argc, char *argv[]);
extern int Nlm_VibMainFinale (int argc, char *argv[]);
#else
extern void Nlm_VibMainPrelude ();
extern void Nlm_VibMainFinale ();
#endif
#endif


#ifdef WIN_MSWIN
extern int CALLBACK Nlm_VibMainPrelude (
  HINSTANCE hInstance,
  HINSTANCE hPrevInstance,
  LPSTR lpszCmdLine,
  int nCmdShow
);
extern int CALLBACK Nlm_VibMainFinale (
  HINSTANCE hInstance,
  HINSTANCE hPrevInstance,
  LPSTR lpszCmdLine,
  int nCmdShow
);
#endif


#ifdef WIN_MOTIF
#if defined(OS_UNIX) && defined(COMP_SUNPRO)
extern void Nlm_VibMainPrelude (int argc, char *argv[], char *envp[]);
extern void Nlm_VibMainFinale (int argc, char *argv[], char *envp[]);
#else
extern void Nlm_VibMainPrelude (int argc, char *argv[]);
extern void Nlm_VibMainFinale (int argc, char *argv[]);
#endif
#endif

#endif /* _VIBMAIN_ */

