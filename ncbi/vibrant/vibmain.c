/*   vibmain.c
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
* File Name:  vibmain.c
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

#include <vibtypes.h>
#include <vibprocs.h>
#include <vibincld.h>
#include <vibmain.h>



#ifdef WIN_MAC
#ifdef OS_UNIX_DARWIN
int main (int argc, char *argv[]) 
{
  Nlm_VibMainPrelude (argc, argv);

  Nlm_Main ();

  Nlm_VibMainFinale (argc, argv);

  ExitToShell ();
}
# else /* ! OS_UNIX_DARWIN */
void main ()
{
  Nlm_VibMainPrelude ();

  Nlm_Main ();

  Nlm_VibMainFinale ();

  ExitToShell ();
}
#endif
#endif


#ifdef WIN_MSWIN
int CALLBACK WinMain (
  HINSTANCE hInstance,
  HINSTANCE hPrevInstance,
  LPSTR lpszCmdLine,
  int nCmdShow
)
{
  if (! Nlm_VibMainPrelude (hInstance, hPrevInstance, lpszCmdLine, nCmdShow)) return FALSE;

  Nlm_Main ();

  return Nlm_VibMainFinale (hInstance, hPrevInstance, lpszCmdLine, nCmdShow);
}
#endif


#ifdef WIN_MOTIF
#if defined(OS_UNIX) && defined(COMP_SUNPRO)
main (int argc, char *argv[], char *envp[])
{
  Nlm_Int2  retval;

  Nlm_VibMainPrelude (argc, argv, envp);

  retval = Nlm_Main ();

  Nlm_VibMainFinale (argc, argv, envp);

  exit (retval);
}
#else
main (int argc, char *argv[])
{
  Nlm_Int2  retval;

  Nlm_VibMainPrelude (argc, argv);

  retval = Nlm_Main ();

  Nlm_VibMainFinale (argc, argv);

  exit (retval);
}
#endif
#endif


