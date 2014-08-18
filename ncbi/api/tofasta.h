/*  tofasta.h
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
* File Name:  tofasta.h
*
* Author:  James Ostell
*   
* Version Creation Date: 7/12/91
*
* $Revision: 6.34 $
*
* File Description:  various sequence objects to fasta output
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_Tofasta_
#define _NCBI_Tofasta_

#include <seqport.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

	                            /* keys returned by FastaWriteFunc */

#define FASTA_ID ((Int2)1)       /* the SeqId string */
#define FASTA_DEFLINE ((Int2)2)	 /* the definition line */
#define FASTA_SEQLINE ((Int2)3)	 /* a line of sequence */
#define FASTA_EOS ((Int2)4)		 /* all the sequence has been printed */
#define FASTA_FORMATDB_AMB ((Int2)5) /* make conversion to standard form for nucleotides */

typedef Boolean (* FastaWriteFunc) PROTO ((BioseqPtr bsp, Int2 key,
                                           CharPtr buf, Uint4 buflen, Pointer data));

typedef struct iteminfo {   /* struct used for defline */
	Uint2	entityID,	   
			itemtype;
    Uint4   itemID;
} ItemInfo, PNTR ItemInfoPtr;

typedef struct descrinfo {   /* struct used for defline */
	BioseqPtr bsp;
	ValNodePtr vnp;
	ItemInfoPtr iip;
	Uint1 choice;
} DescrInfo, PNTR DescrInfoPtr;

typedef struct myfsa {   /* struct used for fasta searches */
    CharPtr buf;         /* buffer for strings (suggest 255 minimum) */
    Uint4 buflen,         /* length of buf */
          seqlen;          /* length of sequence lines (must be buflen-1 or less) */
    Pointer mydata;      /* pointer to your own data */
    FastaWriteFunc myfunc;  /* callback to process parts of fasta file */
    BioseqPtr bsp;       /* Bioseq being processed */
    Boolean bad_asn1; /* set if error in input object like mol not set */
    CharPtr accession;   /* used internally for GenPept def lines */
    CharPtr organism;    /* used internally for GenPept/PRF def lines */
    Uint1 order;          /* used to order def lines for BLAST */
    Boolean do_virtual;   /* if TRUE, instantiate virtual sequences */
    Uint1 tech;           /* for MolInfo.tech */
    Boolean no_sequence;  /* used to disable sequence printing */
    Uint1	code;	/* coding of sequence */
    Boolean	formatdb; /* TRUE, if is used in formatdb */
    Boolean printid_general; /* show gi and gnl for SeqId */
    SeqLocPtr seqloc;  /* sub-sequence of interest */
} MyFsa, PNTR MyFsaPtr;

typedef struct tofasta {
    Boolean is_na;
    Boolean got_one;
    MyFsaPtr mfp;
	Uint1 group_segs;
	Int2 last_indent,
		parts, seg;
} FastaDat, PNTR FastaPtr;
    

NLM_EXTERN Boolean BioseqRawToFastaExtra PROTO((BioseqPtr bsp, FILE *fp, Int2 line_length));
NLM_EXTERN Boolean BioseqRawToFastaExtraEx PROTO((BioseqPtr bsp, FILE *fp, Int2 line_length, SeqLocPtr slp));
NLM_EXTERN Boolean BioseqRawToFasta PROTO((BioseqPtr bsp, FILE * fp, Boolean is_na));
NLM_EXTERN Boolean SeqEntryToFasta PROTO((SeqEntryPtr sep, FILE * fp, Boolean is_na));
NLM_EXTERN Boolean SeqEntryToFastaEx PROTO((SeqEntryPtr sep, FILE * fp, Boolean is_na, Boolean printid_general));
NLM_EXTERN Boolean BioseqToFasta PROTO((BioseqPtr bsp, FILE *fp, Boolean is_na));
NLM_EXTERN Boolean BioseqToFastaDump PROTO((BioseqPtr bsp, FILE *fp, Boolean is_na));



NLM_EXTERN void	SeqEntryFasta PROTO ((SeqEntryPtr sep, Pointer data,
                                     Int4 index, Int2 indent));
   

/*****************************************************************************
*
*   SeqEntrysToFasta(sep, fp, is_na, group_segs)
*
*   	group_segs = 0 ... take only raw Bioseqs
*       group_segs = 1 ... group segmented seqs into single entry.. no parts
*       group_segs = 2 ... show only parts of segmented seqs
*       group_segs = 3 ... like 1, but also instantiate virtual sequences
*   
*****************************************************************************/
NLM_EXTERN Boolean SeqEntrysToFasta PROTO((SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs));

NLM_EXTERN Boolean SeqEntrysToFastaX PROTO((SeqEntryPtr sep, MyFsaPtr mfp, Boolean is_na, Uint1 group_segs));
NLM_EXTERN Boolean SeqEntrysToDefline PROTO((SeqEntryPtr sep, 
                                  FILE *fp, Boolean is_na, Uint1 group_segs));
NLM_EXTERN Boolean BioseqRawToFastaX PROTO((BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na));
NLM_EXTERN Boolean BioseqToFastaX PROTO((BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na));

/*****************************************************************************
*
*   BioseqFastaStream (bsp, fp, flags, linelen, blocklen, grouplen, do_defline)
*   BioseqFastaMemStream (bsp, bs, flags, linelen, blocklen, grouplen, do_defline)
*   SeqLocFastaStream (slp, fp, flags, linelen, blocklen, grouplen)
*   CdRegionFastaStream (sfp, fp, flags, linelen, blocklen, grouplen)
*   TranslationFastaStream (sfp, fp, flags, linelen, blocklen, grouplen)
*   SeqEntryFastaStream (sep, fp, flags, linelen, blocklen, grouplen,
*                        do_na, do_aa, master_style)
*
*   	Rapid FASTA generators using SeqPortStream
*
*****************************************************************************/
NLM_EXTERN Int4 BioseqFastaStream (
  BioseqPtr bsp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
);

NLM_EXTERN Int4 BioseqFastaStreamEx (
  BioseqPtr bsp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  Boolean substitute_ids,
  Boolean sorted_protein
);

NLM_EXTERN Int4 BioseqFastaMemStream (
  BioseqPtr bsp,
  ByteStorePtr bs,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
);

NLM_EXTERN Int4 SeqLocFastaStream (
  SeqLocPtr slp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen
);

NLM_EXTERN Int4 CdRegionFastaStream (
  SeqFeatPtr sfp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
);

NLM_EXTERN Int4 TranslationFastaStream (
  SeqFeatPtr sfp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
);

NLM_EXTERN Int4 SeqEntryFastaStream (
  SeqEntryPtr sep,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_na,
  Boolean do_aa,
  Boolean master_style
);

NLM_EXTERN Int4 SeqEntryFastaStreamEx (
  SeqEntryPtr sep,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_na,
  Boolean do_aa,
  Boolean master_style,
  Boolean substitute_ids,
  Boolean sorted_prot
);

/*****************************************************************************
*
*   FastaFileFunc(key, buf, data)
*   	standard "write to file" callback
*
*****************************************************************************/
NLM_EXTERN Boolean FastaFileFunc PROTO((BioseqPtr bsp, Int2 key,
                                        CharPtr buf, Uint4 buflen, Pointer data));


/*****************************************************************************
*
*   FastaDumpFileFunc(key, buf, data)
*       standard "write to file" callback
*
*       Used for BLAST (FASTA) databases.  If the defline is
*       longer than buflen, then check that an ID is not
*       truncated in the middle.
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDumpFileFunc PROTO((BioseqPtr bsp, Int2 key, CharPtr buf, Uint4 buflen, Pointer data));

/*****************************************************************************
*
*   Reads a Fasta File into a SeqEntry structure
*   Conventions:
*   >name def
*   agaggagagagagag
*   agagagagagagagag
*
*   "name" = string is considered SEQID_LOCAL until first white space
*   "def"  = after first white space, and before first newline will be "title"
*   "agaga.." = sequence follows until EOF. can be upper or lower case IUPAC
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntry PROTO((FILE *fp, Boolean is_na));

/*****************************************************************************
*
*   Reads a Fasta Buffer into a SeqEntry structure
*   Conventions:
*   >name def
*   agaggagagagagag
*   agagagagagagagag
*
*   "name" = string is considered SEQID_LOCAL until first white space
*   "def"  = after first white space, and before first newline will be "title"
*   "agaga.." = sequence follows until EOF. can be upper or lower case IUPAC
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuff PROTO((CharPtr buffer, 
                                  CharPtr PNTR last_char,Boolean is_na));

NLM_EXTERN SeqEntryPtr FastaToSeqBuffEx(CharPtr buffer, CharPtr PNTR last_char, 
                             Boolean is_na, CharPtr PNTR errormsg,
                             Boolean parseSeqId);
NLM_EXTERN SeqEntryPtr FastaToSeqBuffForDb(CharPtr buffer, CharPtr PNTR last_char, 
                             Boolean is_na, CharPtr PNTR errormsg,
                             Boolean parseSeqId, CharPtr prefix, Int2Ptr ctrptr,
			     SeqLocPtr PNTR mask_ptr);
NLM_EXTERN SeqEntryPtr FastaToSeqEntryEx (FILE *fp, Boolean is_na, 
                               CharPtr PNTR errormsg,
                               Boolean parseSeqId); 
NLM_EXTERN SeqEntryPtr FastaToSeqEntryForDb ( FILE *fp, Boolean is_na,
                               CharPtr PNTR errormsg,
                               Boolean parseSeqId,
                               CharPtr prefix, Int2Ptr ctrptr, 
                               SeqLocPtr PNTR mask_ptr);

/********* DEFINES for input type *********/

#define FASTA_MEM_IO  1   /* type of reading from buffer in memory */
#define FASTA_FILE_IO 2   /* type of reading from file */

NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternalEx ( VoidPtr input,
                               Int4 type, CharPtr PNTR next_char,
                               Boolean is_na, CharPtr PNTR errormsg,
                               Boolean parseSeqId, CharPtr special_symbol,
                               CharPtr prefix, Int2Ptr ctrptr,
                               SeqLocPtr PNTR mask_ptr);

/*****************************************************************************
*
*   Boolean FastaReadSequence() - read sequence from file
*
*****************************************************************************/

Boolean FastaReadSequence
(
 FILE *fd,                 /* input file pointer) */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
 );

/*****************************************************************************
*
*   Boolean FastaReadSequenceMem() - read sequence from buffer
*
*****************************************************************************/

Boolean FastaReadSequenceMem
(
 CharPtr buffer,           /* input buffer with sequence */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
);

/*****************************************************************************
*
*   This routines get parts needed to make FASTA format from ASN.1
*
*****************************************************************************/
/*****************************************************************************
*
*   FastaId(bsp, buf, buflen)
*      Makes the string for the id part of fasta format.
*      buf should be at least 40 bytes
*
*****************************************************************************/
NLM_EXTERN Boolean FastaId PROTO((BioseqPtr bsp, CharPtr buf, Uint4 buflen));

/*****************************************************************************
*
*   FastaDefLine(bsp, buf, buflen, accession, organism)
*   	Finds or makes a FASTA format defline (just locates the string)
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*       if (accession != NULL) prefixes defline with (accession)
*          used for translated GenBank records
*       if (organism != NULL) adds [organism] to end
*       if (tech == MI_TECH_phase1 or phase2, adds order comment to defline)
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDefLine PROTO((BioseqPtr bsp, CharPtr buf, Uint4 buflen, CharPtr accession, CharPtr organism, Uint1 tech));

NLM_EXTERN Boolean CreateDefLine PROTO((ItemInfoPtr dip, BioseqPtr bsp, CharPtr buf, Uint4 buflen,
                                        Uint1 tech, CharPtr accession, CharPtr organism));
NLM_EXTERN Boolean CreateDefLineEx (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Uint4 buflen, Uint1 tech,
                                    CharPtr accession, CharPtr organism, Boolean ignoreTitle);
NLM_EXTERN Boolean CreateDefLineExEx (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Uint4 buflen, Uint1 tech,
                                      CharPtr accession, CharPtr organism, Boolean ignoreTitle, Boolean extProtTitle);

/*****************************************************************************
*
*   NewCreateDefLine and NewCreateDefLineBuf replace the CreateDefLine family
*
*****************************************************************************/

NLM_EXTERN CharPtr NewCreateDefLine (
  ItemInfoPtr iip,
  BioseqPtr bsp,
  Boolean ignoreTitle,
  Boolean extProtTitle
);

NLM_EXTERN Boolean NewCreateDefLineBuf (
  ItemInfoPtr iip,
  BioseqPtr bsp,
  CharPtr buf,
  Uint4 buflen,
  Boolean ignoreTitle,
  Boolean extProtTitle
);

/*****************************************************************************
*
*   FastaSeqPort(bsp, is_na, do_virtual)
*   	opens a SeqPort for a fasta output of bsp
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr FastaSeqPort PROTO((BioseqPtr bsp, Boolean is_na,
                                          Boolean do_virtual, Uint1 code));
NLM_EXTERN SeqPortPtr FastaSeqPortEx PROTO((BioseqPtr bsp, Boolean is_na, 
        Boolean do_virtual, Uint1 code, SeqLocPtr slp));

/*****************************************************************************
*
*   FastaSeqLine(spp, buf, linelen)
*     an open seqport is passed in.
*     fills buf with linelen bases
*     assumes buf[linelen] = '\0'
*     returns FALSE when no more residues to print
*     here for backward compatibility. Is the same as calling
*      FastaSeqLineEx() with do_virtual = FALSE;
*
*****************************************************************************/
NLM_EXTERN Boolean FastaSeqLine PROTO((SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na));

/****************************************************************************
*
*  FastaSeqLineEx(spp, buf, linelen, is_na, do_virtual)
*
*    adds a parameter to indicate if virtual sequences should be shown
*     this is necessary in order to avoid showing SEQPORT_VIRT return
*     (gap of 0 length) as '-', when do_virtual is TRUE. 
*
*****************************************************************************/
NLM_EXTERN Boolean FastaSeqLineEx(SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na, Boolean
do_virtual);

Int4 GetOrderBySeqId(Int4 choice, Boolean is_prot);

/*****************************************************************************
*
*   NC_Cleanup (entityID, ptr) and ClearGenBankKeywords (entityID, ptr)
*     internal functions for genome RefSeq processing
*
*****************************************************************************/
NLM_EXTERN void NC_Cleanup (Uint2 entityID, Pointer ptr);
NLM_EXTERN void ClearGenBankKeywords (Uint2 entityID, Pointer ptr);

/*****************************************************************************
*
*   InstantiateProteinTitles (entityID, ptr) and ClearProteinTitles (entityID, ptr)
*     allows proteins titles to be kept as Seq_descr_title rather than always
*     being generated on the fly
*
*****************************************************************************/
NLM_EXTERN void InstantiateProteinTitles (Uint2 entityID, Pointer ptr);
NLM_EXTERN void ClearProteinTitlesInNucProts (Uint2 entityID, Pointer ptr);

NLM_EXTERN CharPtr MakeCompleteChromTitle (BioseqPtr bsp, Uint1 biomol, Uint1 completeness);


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
