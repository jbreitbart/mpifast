/*   cdrgn.h
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
* File Name:  cdrgn.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.13 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: cdrgn.h,v $
* Revision 6.13  2008/01/17 21:18:09  bollin
* Fixes for RNA dialog for ncRNA product/class handling
*
* Revision 6.12  2007/12/21 17:30:27  bollin
* Changes to dialog for selecting ncRNA class and RNA type - allow ncRNA class
* Any if selecting RNA type in a constraint context, test dialog to generate
* error if "other" is the class but no text is specified.
*
* Revision 6.11  2007/11/26 21:16:29  bollin
* Moved CreatencRNAClassDialog proto into cdrgn.h
*
* Revision 6.10  2007/09/21 18:05:23  bollin
* code in place for conversion of old-style misc_RNA, snRNA, scRNA, and snoRNA
* features to new ncRNA features, commented out until changeover.
*
* Revision 6.9  2007/09/10 20:08:51  bollin
* Correction to MatchesRnaType
*
* Revision 6.8  2007/09/10 18:47:11  kans
* prototype for SetRnaSpecificQuals, cast for codon argument to ParseTRnaString
*
* Revision 6.7  2007/09/10 18:33:15  bollin
* Changes for new ncRNA and tmRNA class editor.
*
* Revision 6.6  2007/07/25 13:53:12  bollin
* Removed local copy of TruncateLocation from sequin3.c, made TruncateLocation
* function in desktop/cdrgn.c extern and added prototype to cdrgn.h
*
* Revision 6.5  2003/10/23 16:43:56  kans
* changed operon to import feature
*
* Revision 6.4  2003/10/07 13:52:01  kans
* added gap, operon, oriT features and ecotype, estimated_length and operon qualifiers
*
* Revision 6.3  2000/07/08 20:44:00  vakatov
* Get all "#include" out of the 'extern "C" { }' scope;  other cleanup...
*
* ==========================================================================
*/

#ifndef _CDRGN_
#define _CDRGN_

#include <dlogutil.h>

#ifdef __cplusplus
extern "C" {
#endif

#define REGISTER_CDRGN_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit CdRgn","CDS",OBJ_SEQFEAT,FEATDEF_CDS,OBJ_SEQFEAT,FEATDEF_CDS,NULL,CdRgnGenFunc,PROC_PRIORITY_DEFAULT)

extern ForM CreateCdRgnForm (Int2 left, Int2 top, CharPtr title,
                             SeqFeatPtr sfp, SeqEntryPtr sep,
                             FormActnFunc actproc);
extern Int2 LIBCALLBACK CdRgnGenFunc (Pointer data);

extern void CdRgnFeatFormActnProc (ForM f);
extern void CdRgnTranslateWithFrame (ForM f, Uint1 frame);

extern SeqLocPtr PredictCodingRegion (BioseqPtr nuc, BioseqPtr prot, Int2 genCode);

#define REGISTER_GENE_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit Gene","Gene",OBJ_SEQFEAT,FEATDEF_GENE,OBJ_SEQFEAT,FEATDEF_GENE,NULL,GeneGenFunc,PROC_PRIORITY_DEFAULT)

extern ForM CreateGeneForm (Int2 left, Int2 top, CharPtr title,
                            SeqFeatPtr sfp, SeqEntryPtr sep,
                            FormActnFunc actproc);
extern Int2 LIBCALLBACK GeneGenFunc (Pointer data);

#define REGISTER_PROT_EDIT(PROCNAME,PROCLABEL,SUBTYPE) ObjMgrProcLoad(OMPROC_EDIT,PROCNAME,PROCLABEL,OBJ_SEQFEAT,SUBTYPE,OBJ_SEQFEAT,SUBTYPE,NULL,ProtGenFunc,PROC_PRIORITY_DEFAULT)

extern ForM CreateProtForm (Int2 left, Int2 top, CharPtr title,
                            SeqFeatPtr sfp, SeqEntryPtr sep,
                            FormActnFunc actproc);
extern Int2 LIBCALLBACK ProtGenFunc (Pointer data);

#define REGISTER_RNA_EDIT(PROCNAME,PROCLABEL,SUBTYPE) ObjMgrProcLoad(OMPROC_EDIT,PROCNAME,PROCLABEL,OBJ_SEQFEAT,SUBTYPE,OBJ_SEQFEAT,SUBTYPE,NULL,RnaGenFunc,PROC_PRIORITY_DEFAULT)

extern ForM CreateRnaForm (Int2 left, Int2 top, CharPtr title,
                           SeqFeatPtr sfp, SeqEntryPtr sep,
                           Uint2 subtype, FormActnFunc actproc);
extern Int2 LIBCALLBACK RnaGenFunc (Pointer data);

extern void AddRnaSpecificQuals (SeqFeatPtr sfp, DialoG d);
extern void ConvertProductQualToRnaRefName (SeqFeatPtr sfp);
extern void SetRnaSpecificQuals (SeqFeatPtr sfp, DialoG d);
extern void ConvertToOldRNAFormat (SeqFeatPtr sfp);

extern SeqLocPtr TruncateLocation (SeqLocPtr head, Int4 len);

/* for searching for RNA values of a certain type */
typedef struct rnatype {
  Int4    rna_featdef; /* use FEATDEF_ANY for match any RNA */
  CharPtr ncrna_class; /* value to look for in ncrna_class qual */
} RnaTypeData, PNTR RnaTypePtr;

extern RnaTypePtr RnaTypeFree (RnaTypePtr rtp);
extern Boolean MatchesRnaType (SeqFeatPtr sfp, RnaTypePtr rtp);
extern void ApplyRnaTypeToSeqFeat (SeqFeatPtr sfp, RnaTypePtr rtp);
extern void ApplyProductToRNA (SeqFeatPtr sfp, CharPtr product);
extern void AddToComment (SeqFeatPtr sfp, CharPtr comment);
extern DialoG RnaTypeDialog (GrouP h, Boolean is_constraint, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
extern DialoG CreatencRNAClassDialog (GrouP h, Boolean is_constraint, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

#ifdef __cplusplus
}
#endif

#endif /* ndef _CDRGN_ */

