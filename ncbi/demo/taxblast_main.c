static char const rcsid[] = "$Id: taxblast_main.c,v";

/* $Id: taxblast_main.c,v
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* File Name:  $RCSfile: taxblast_main.c,v $
*
* Authors:  Tom Madden
*
* ==========================================================================
*/

#include <ncbi.h>
#include <sequtil.h>
#include <treemgr.h>
#include <taxext.h>
#include <txclient.h>
#include <objseq.h>
#include <objgen.h>
#include <taxblast.h>


#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    { "Input ASN.1 File (SeqAnnot)",             /* 0 */
      NULL, NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL },
    { "Sequence is DNA",                         /* 1 */
      "F", NULL, NULL, TRUE, 'p', ARG_BOOLEAN, 0.0, 0, NULL },
    { "Database used to get SeqAnnot ASN.1",     /* 2 */
      "nr", NULL, NULL, TRUE, 'd', ARG_STRING, 0.0, 0, NULL },
    { "Output file name",                        /* 3 */
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL }
};

Int2 Main (void)
{
    AsnIoPtr aip;
    SeqAnnotPtr sap;
    Boolean is_na = FALSE;
    FILE *outfile;
    Char ofile[128];
    
    if (!GetArgs("txblast", NUMARG, myargs)) {
        return 1;
    }
    
    if (myargs[1].intvalue) 
        is_na = TRUE;
    
    if((aip = AsnIoOpen(myargs[0].strvalue, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 1,0, "AsnIoOpen failure\n");
        return 1;
    }
    
    if((sap = SeqAnnotAsnRead (aip, NULL)) == NULL) {
        ErrPostEx(SEV_FATAL, 1,0,"SeqAlignAsnRead failure\n");
        return 1;
    }

    if(StringCmp(myargs[3].strvalue, "stdout")) {
        sprintf (ofile, "%s.html", myargs[0].strvalue);
        outfile = FileOpen(ofile, "w");
    } else {
        outfile = FileOpen(myargs[3].strvalue, "w");
    }
    
    TXBHtmlReport((SeqAlignPtr)sap->data, outfile, is_na, is_na, 
                  myargs[2].strvalue, NULL, NULL, FALSE);
    
    FileClose(outfile);
    
    AsnIoClose(aip);
    SeqAnnotFree(sap);
    
    return (0);
}
