/*  objalign.h
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
* File Name:  objalign.h
*
* Author:  James Ostell
*   
* Version Creation Date: 4/1/91
*
* $Revision: 6.11 $
*
* File Description:  Object manager interface for module NCBI-Seqalign
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: objalign.h,v $
* Revision 6.11  2008/10/09 18:43:56  kans
* splice_5_prime and splice_3_prime changed to acceptor_before_exon and donor_after_exon
*
* Revision 6.10  2007/05/12 21:03:32  kans
* Spliced-seg code recompiled with -Z so product_length 0 is not always written
*
* Revision 6.9  2007/04/09 20:25:47  kans
* added support for sparse alignment type
*
* Revision 6.8  2006/07/28 16:07:48  kans
* added new fields due to spec change
*
* Revision 6.7  2002/01/10 14:35:18  dondosha
* Added GenericSeqAlignSetAsnWrite to allow writing seqalign set without a parent asn type
*
* Revision 6.6  1999/09/28 12:09:56  kans
* added alignID field
*
* Revision 6.5  1999/09/27 17:48:37  kans
* using GatherIndex structure
*
* Revision 6.4  1999/09/24 23:09:23  kans
* adds EXTRA_OBJMGR_FIELDS to several objects
*
* Revision 6.3  1999/09/07 17:00:26  kans
* added entityID, itemID, itemtype fields for new Alignment Indexing functions
*
* Revision 6.2  1999/07/29 15:49:58  ostell
* added pointer and free for a SeqAlignIndex
*
* Revision 6.1  1999/07/26 20:41:01  ostell
* added SAT_ and SAS_ defines, added master to SeqAlign
*
* Revision 6.0  1997/08/25 18:49:14  madden
* Revision changed to 6.0
*
* Revision 4.1  1997/06/19 18:40:41  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 4.0  1995/07/26 13:48:06  ostell
* force revision to 4.0
*
 * Revision 3.5  1995/07/22  21:59:13  ostell
 * added support for ASN.1 spec 4.0
 *
 * Revision 3.4  1995/06/21  17:17:50  epstein
 * *** empty log message ***
 *
 * Revision 3.3  95/05/15  21:22:00  ostell
 * added Log line
 * 
*
*
*
* ==========================================================================
*/

#ifndef _NCBI_Seqalign_
#define _NCBI_Seqalign_

#ifndef _ASNTOOL_
#include <asn.h>
#endif
#ifndef _NCBI_General_
#include <objgen.h>
#endif
#ifndef _NCBI_Seqloc_
#include <objloc.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
*
*   loader
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqAlignAsnLoad PROTO((void));

/*****************************************************************************
*
*   internal structures for NCBI-Seqalign objects
*
*****************************************************************************/

/*****************************************************************************
*
*   Score
*     NOTE: read, write, and free always process GROUPS of scores
*
*****************************************************************************/
typedef struct score {
    ObjectIdPtr id;
    Uint1 choice;          /* 0=not set, 1=int, 2=real */
    DataVal value;
    struct score PNTR next;    /* for sets of scores */
} Score, PNTR ScorePtr;

NLM_EXTERN ScorePtr LIBCALL ScoreNew PROTO((void));
NLM_EXTERN Boolean  LIBCALL ScoreSetAsnWrite PROTO((ScorePtr sp, AsnIoPtr aip, AsnTypePtr settype));
NLM_EXTERN ScorePtr LIBCALL ScoreSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr settype));
NLM_EXTERN ScorePtr LIBCALL ScoreSetFree PROTO((ScorePtr anp));

/****************************************************************************
*
*  SeqAlignIndex
*
*    This structure is the public face of a data structure which is
*    extended in AlignMgr for alignment indexing, alignment features
*    and other utilities.
*
*    It is left as a limited structure "stub" here
*
****************************************************************************/
/** the VoidPtr below should really be a SeqAlignIndexPtr **/
typedef Boolean (LIBCALLBACK * SeqAlignIndexFreeFunc)(VoidPtr);

typedef struct seqalignindex {
	Uint1 indextype;
	SeqAlignIndexFreeFunc freefunc;
} SeqAlignIndex, PNTR SeqAlignIndexPtr;

NLM_EXTERN SeqAlignIndexPtr LIBCALL SeqAlignIndexFree (SeqAlignIndexPtr saip);
    
/*****************************************************************************
*
*   SeqAlign
*   type =  type of alignment
        not-set (0) ,
        global (1) ,
        diags (2) ,
        partial (3) ,           -- mapping pieces together
		disc (4) ,
        other (255) } ,
    segtype = type of segs structure
        not-set 0
        dendiag 1
        denseq 2
        std 3
		packed 4
		disc 5      SeqAlignSet is used
		spliced 6
*   
*
*****************************************************************************/
/** SeqAlign.type values ***/
#define SAT_GLOBAL 1    /* ordered segments, over full length of seqs */
#define SAT_DIAGS 2     /* unordered, possibly overlapping segments */
#define SAT_PARTIAL 3   /* ordered segments, over part of sequence */
#define SAT_MASTERSLAVE 4   /* set of SeqAligns, all of which have one common */
                            /* sequence. Not in ASN.1 yet */

/** SeqAlign.segtype values ***/
#define SAS_DENDIAG 1
#define SAS_DENSEG 2
#define SAS_STD 3
#define SAS_PACKED 4
#define SAS_DISC 5
#define SAS_SPLICED 6
#define SAS_SPARSE 7


typedef struct seqalign {
    Uint1 type,
        segtype;
    Int2 dim;
    ScorePtr score;
    Pointer segs;
    struct seqalign PNTR next;
	SeqLocPtr bounds;      /* sequence of SeqLocPtr */
	ValNodePtr   id;
	struct struct_User_object PNTR   ext;
    SeqIdPtr master;   /* for SAT_MASTERSLAVE */
    SeqAlignIndexPtr saip;  /* for added Alignment Indexing structures */
	GatherIndex idx;      /* internal gather/objmgr tracking fields */
	Uint2 alignID;        /* unique number assigned to alignment */
} SeqAlign, PNTR SeqAlignPtr;

NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignNew PROTO((void));
NLM_EXTERN Boolean     LIBCALL SeqAlignAsnWrite PROTO((SeqAlignPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignFree PROTO((SeqAlignPtr anp));
NLM_EXTERN Int2 LIBCALL SeqAlignLabel PROTO((SeqAlignPtr sap, CharPtr buffer, Int2 buflen, Uint1 content));

/*****************************************************************************
*
*   SeqAlignSet
*
*****************************************************************************/
NLM_EXTERN Boolean     LIBCALL SeqAlignSetAsnWrite PROTO((SeqAlignPtr anp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));
NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));
NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignSetFree PROTO((SeqAlignPtr sap));
NLM_EXTERN SeqAlignPtr LIBCALL SeqAlignSetNew PROTO((void));
NLM_EXTERN SeqAlignPtr LIBCALL SpecialSeqAlignSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr set));
NLM_EXTERN Boolean LIBCALL SpecialSeqAlignSetAsnWrite PROTO((SeqAlignPtr sap, AsnIoPtr aip, AsnTypePtr set));
NLM_EXTERN Boolean LIBCALL GenericSeqAlignSetAsnWrite PROTO((SeqAlignPtr sap, AsnIoPtr aip));


/*****************************************************************************
*
*   DenseDiag
*   
*
*****************************************************************************/
typedef struct dendiag {
    Int2 dim;                   /* this is a convenience, not in asn1 */
    SeqIdPtr id;
    Int4Ptr starts;
    Int4 len;
    Uint1Ptr strands;
    ScorePtr scores;
    struct dendiag PNTR next;
} DenseDiag, PNTR DenseDiagPtr;

NLM_EXTERN DenseDiagPtr LIBCALL DenseDiagNew PROTO((void));
NLM_EXTERN Boolean      LIBCALL DenseDiagAsnWrite PROTO((DenseDiagPtr ddp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DenseDiagPtr LIBCALL DenseDiagAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DenseDiagPtr LIBCALL DenseDiagFree PROTO((DenseDiagPtr ddp));

/*****************************************************************************
*
*   DenseSeg
*   
*
*****************************************************************************/
typedef struct denseg {
    Int2 dim,
        numseg;
    SeqIdPtr ids;           /* dimension is dim */
    Int4Ptr starts;			/* dimension is dim * numseg */
    Int4Ptr lens;			/* dimension is numseg */
    Uint1Ptr strands;		/* dimension is dim * numseg */
    ScorePtr scores;		/* dimension is numseg */
} DenseSeg, PNTR DenseSegPtr;

NLM_EXTERN DenseSegPtr LIBCALL DenseSegNew PROTO((void));
NLM_EXTERN Boolean     LIBCALL DenseSegAsnWrite PROTO((DenseSegPtr dsp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DenseSegPtr LIBCALL DenseSegAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DenseSegPtr LIBCALL DenseSegFree PROTO((DenseSegPtr dsp));

/*****************************************************************************
*
*   PackSeg
*   
*
*****************************************************************************/
typedef struct packseg {
    Int2 dim,
        numseg;
    SeqIdPtr ids;			/* dimension is dim */
    Int4Ptr starts;			/* dimension is dim */
	ByteStorePtr present;	/* dimension is dim * numseg booleans */
    Int4Ptr lens;			/* dimension is numseg */
    Uint1Ptr strands;		/* dimension is dim */
    ScorePtr scores;		/* dimension is numseg */
} PackSeg, PNTR PackSegPtr;

NLM_EXTERN PackSegPtr LIBCALL PackSegNew PROTO((void));
NLM_EXTERN Boolean     LIBCALL PackSegAsnWrite PROTO((PackSegPtr psp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN PackSegPtr LIBCALL PackSegAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN PackSegPtr LIBCALL PackSegFree PROTO((PackSegPtr psp));

/*****************************************************************************
*
*   StdSeg
*   
*
*****************************************************************************/
typedef struct stdseg {
    Int2 dim;
    SeqIdPtr ids;    /* SeqId s */
    SeqLocPtr loc;    /* SeqLoc s */
    ScorePtr scores;
    struct stdseg PNTR next;
} StdSeg, PNTR StdSegPtr;

NLM_EXTERN StdSegPtr LIBCALL StdSegNew PROTO((void));
NLM_EXTERN Boolean   LIBCALL StdSegAsnWrite PROTO((StdSegPtr ssp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN StdSegPtr LIBCALL StdSegAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN StdSegPtr LIBCALL StdSegFree PROTO((StdSegPtr ssp));

/**************************************************
*
*    SplicedSeg
*
**************************************************/
typedef struct struct_Spliced_seg {
   Uint4 OBbits__;
   ValNodePtr   product_id;
   ValNodePtr   genomic_id;
#define OB__Spliced_seg_product_strand 0

   Uint2   product_strand;
#define OB__Spliced_seg_genomic_strand 1

   Uint2   genomic_strand;
   Uint2   product_type;
   /* following #defines are for enumerated type, not used by object loaders */
#define Spliced_seg_product_type_transcript 0
#define Spliced_seg_product_type_protein 1

   struct struct_Spliced_exon PNTR   exons;
#define OB__Spliced_seg_poly_a 2

   Int4   poly_a;
#define OB__Spliced_seg_product_length 3

   Int4   product_length;
   ValNodePtr   modifiers;
} SplicedSeg, PNTR SplicedSegPtr;


NLM_EXTERN SplicedSegPtr LIBCALL SplicedSegFree PROTO ((SplicedSegPtr ));
NLM_EXTERN SplicedSegPtr LIBCALL SplicedSegNew PROTO (( void ));
NLM_EXTERN SplicedSegPtr LIBCALL SplicedSegAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SplicedSegAsnWrite PROTO (( SplicedSegPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SparseSeg
*
**************************************************/
typedef struct struct_Sparse_seg {
   Uint4 OBbits__;
   ValNodePtr   master_id;
   struct struct_Sparse_align PNTR   rows;
   ScorePtr    row_scores;
   struct struct_Sparse_seg_ext PNTR   ext;
} SparseSeg, PNTR SparseSegPtr;


NLM_EXTERN SparseSegPtr LIBCALL SparseSegFree PROTO ((SparseSegPtr ));
NLM_EXTERN SparseSegPtr LIBCALL SparseSegNew PROTO (( void ));
NLM_EXTERN SparseSegPtr LIBCALL SparseSegAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SparseSegAsnWrite PROTO (( SparseSegPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SplicedExon
*
**************************************************/
typedef struct struct_Spliced_exon {
   struct struct_Spliced_exon PNTR next;
   ValNodePtr   product_start;
   ValNodePtr   product_end;
   Int4   genomic_start;
   Int4   genomic_end;
   ValNodePtr   product_id;
   ValNodePtr   genomic_id;
   Uint2   product_strand;
   Uint2   genomic_strand;
   ValNodePtr   parts;
   ScorePtr   scores;
   struct struct_Splice_site PNTR   acceptor_before_exon;
   struct struct_Splice_site PNTR   donor_after_exon;
   Uint1   partial;
   struct struct_User_object PNTR   ext;
} SplicedExon, PNTR SplicedExonPtr;


NLM_EXTERN SplicedExonPtr LIBCALL SplicedExonFree PROTO ((SplicedExonPtr ));
NLM_EXTERN SplicedExonPtr LIBCALL SplicedExonNew PROTO (( void ));
NLM_EXTERN SplicedExonPtr LIBCALL SplicedExonAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SplicedExonAsnWrite PROTO (( SplicedExonPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr SplicedSegModifierPtr;
typedef ValNode SplicedSegModifier;
#define SplicedSegModifier_start_codon_found 1
#define SplicedSegModifier_stop_codon_found 2


NLM_EXTERN SplicedSegModifierPtr LIBCALL SplicedSegModifierFree PROTO ((SplicedSegModifierPtr ));
NLM_EXTERN SplicedSegModifierPtr LIBCALL SplicedSegModifierAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SplicedSegModifierAsnWrite PROTO (( SplicedSegModifierPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ProductPosPtr;
typedef ValNode ProductPos;
#define ProductPos_nucpos 1
#define ProductPos_protpos 2


NLM_EXTERN ProductPosPtr LIBCALL ProductPosFree PROTO ((ProductPosPtr ));
NLM_EXTERN ProductPosPtr LIBCALL ProductPosAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ProductPosAsnWrite PROTO (( ProductPosPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr SplicedExonChunkPtr;
typedef ValNode SplicedExonChunk;
#define SplicedExonChunk_match 1
#define SplicedExonChunk_mismatch 2
#define SplicedExonChunk_diag 3
#define SplicedExonChunk_product_ins 4
#define SplicedExonChunk_genomic_ins 5


NLM_EXTERN SplicedExonChunkPtr LIBCALL SplicedExonChunkFree PROTO ((SplicedExonChunkPtr ));
NLM_EXTERN SplicedExonChunkPtr LIBCALL SplicedExonChunkAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SplicedExonChunkAsnWrite PROTO (( SplicedExonChunkPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SpliceSite
*
**************************************************/
typedef struct struct_Splice_site {
   CharPtr   bases;
} SpliceSite, PNTR SpliceSitePtr;


NLM_EXTERN SpliceSitePtr LIBCALL SpliceSiteFree PROTO ((SpliceSitePtr ));
NLM_EXTERN SpliceSitePtr LIBCALL SpliceSiteNew PROTO (( void ));
NLM_EXTERN SpliceSitePtr LIBCALL SpliceSiteAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SpliceSiteAsnWrite PROTO (( SpliceSitePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ProtPos
*
**************************************************/
typedef struct struct_Prot_pos {
   Int4   amin;
   Int4   frame;
} ProtPos, PNTR ProtPosPtr;


NLM_EXTERN ProtPosPtr LIBCALL ProtPosFree PROTO ((ProtPosPtr ));
NLM_EXTERN ProtPosPtr LIBCALL ProtPosNew PROTO (( void ));
NLM_EXTERN ProtPosPtr LIBCALL ProtPosAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ProtPosAsnWrite PROTO (( ProtPosPtr , AsnIoPtr, AsnTypePtr));


/**************************************************
*
*    SparseAlign
*
**************************************************/
typedef struct struct_Sparse_align {
   struct struct_Sparse_align PNTR next;
   Uint4 OBbits__;
   ValNodePtr   first_id;
   ValNodePtr   second_id;
   Int4   numseg;
   ValNodePtr   first_starts;
   ValNodePtr   second_starts;
   ValNodePtr   lens;
   ValNodePtr   second_strands;
   ScorePtr   seg_scores;
} SparseAlign, PNTR SparseAlignPtr;


NLM_EXTERN SparseAlignPtr LIBCALL SparseAlignFree PROTO ((SparseAlignPtr ));
NLM_EXTERN SparseAlignPtr LIBCALL SparseAlignNew PROTO (( void ));
NLM_EXTERN SparseAlignPtr LIBCALL SparseAlignAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SparseAlignAsnWrite PROTO (( SparseAlignPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SparseSegExt
*
**************************************************/
typedef struct struct_Sparse_seg_ext {
   struct struct_Sparse_seg_ext PNTR next;
   Uint4 OBbits__;
   Int4   index;
} SparseSegExt, PNTR SparseSegExtPtr;


NLM_EXTERN SparseSegExtPtr LIBCALL SparseSegExtFree PROTO ((SparseSegExtPtr ));
NLM_EXTERN SparseSegExtPtr LIBCALL SparseSegExtNew PROTO (( void ));
NLM_EXTERN SparseSegExtPtr LIBCALL SparseSegExtAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SparseSegExtAsnWrite PROTO (( SparseSegExtPtr , AsnIoPtr, AsnTypePtr));

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
