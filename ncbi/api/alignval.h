/*  alignval.hA
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
* File Name:  alignval.h
*
* Author:  Jian Ye, Colombe Chappey
*
* Version Creation Date:   6/3/99
*
* $Revision: 6.18 $
*
* File Description:
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: alignval.h,v $
* Revision 6.18  2007/07/05 17:50:49  bollin
* Added validation INFO when alignment stops before ends of sequences.
* Added red vertical bar to columns in alignment assistant when weighted
* percent identity for column is less than 50%.
*
* Revision 6.17  2007/02/28 21:07:13  bollin
* Added function for weighted percent identity for alignments
*
* Revision 6.16  2006/10/23 15:11:42  bollin
* Moved AlignmentPercentIdentity and supporting functions here from
* tools/salptool.c to avoid library dependency problems.
*
* Revision 6.15  2006/10/20 13:30:22  bollin
* Added Validator error for alignment percent identity.
*
* Revision 6.14  2003/11/14 18:06:42  kans
* added do_hist_assembly parameter
*
* Revision 6.13  1999/11/23 21:47:31  vakatov
* Fixed for C++ and/or DLL compilation
*
* ==========================================================================
*/

#ifndef ALIGNVAL_H
#define ALIGNVAL_H

#include <ncbi.h>
#include <objall.h>
#include <objseq.h>
#include <objmgr.h>
#include <objfdef.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif


/*call back function for REGISTER_ALIGNVALIDATION defined in sequin4.c.  
Starting point for seqalignment validation if user clicked on SeqalignValidation 
under menu Filer/Alignment.  Either individual alignment or alignment block 
should be highlighted for this validation to work*/

NLM_EXTERN Int2 LIBCALLBACK ValidateSeqAlignFromData (Pointer data);

/*validate each alignment sequentially.  This function will subject the seqalign to all validation functions*/ 
NLM_EXTERN Boolean ValidateSeqAlign (SeqAlignPtr salp, Uint2 entityID, Boolean message,
                         Boolean msg_success, Boolean find_remote_bsp,
                         Boolean delete_bsp, Boolean delete_salp, BoolPtr dirty);

NLM_EXTERN Boolean ValidateSeqAlignInSeqEntry (SeqEntryPtr sep, Boolean message, 
                                 Boolean msg_success, Boolean find_remote_bsp, 
                                 Boolean delete_bsp, Boolean delete_salp,
                                 Boolean do_hist_assembly);

extern Uint2 AlignmentPercentIdentity (SeqAlignPtr salp, Boolean internal_gaps);
extern Uint2 WeightedAlignmentPercentIdentity (SeqAlignPtr salp, Boolean internal_gaps);

extern double *
GetAlignmentColumnPercentIdentities 
(SeqAlignPtr salp,
 Int4    start,
 Int4    stop,
 Boolean internal_gaps,
 Boolean internal_validation);


#define Err_SeqId 1
#define Err_Strand_Rev 2
#define Err_Denseg_Len_Start 3
#define Err_Start_Less_Than_Zero 4
#define Err_Start_More_Than_Biolen 5
#define Err_End_Less_Than_Zero 6
#define Err_End_More_Than_Biolen 7
#define Err_Len_Less_Than_Zero 8
#define Err_Len_More_Than_Biolen 9
#define Err_Sum_Len_Start 10
#define Err_SeqAlign_DimSeqId_Not_Match 11
#define Err_Segs_DimSeqId_Not_Match 12
#define Err_Fastalike 13
#define Err_Null_Segs 14
#define Err_Segment_Gap 15 
#define Err_Segs_Dim_One 16
#define Err_SeqAlign_Dim_One 17
#define Err_Segtype 18
#define Err_Pcnt_ID 20
#define Err_Short_Aln 21

#ifdef __cplusplus
}
#endif


#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
 

