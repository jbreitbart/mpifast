/*  valid.h
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
* File Name:  valid.h
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/94
*
* $Revision: 6.50 $
*
* File Description:  Sequence editing utilities
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_Valid_
#define _NCBI_Valid_

/*****************************************************************************
*
*   valid.h
*       values for cutoff are
*       0 INFO
*       1 WARN
*       2 ERROR
*       3 FATAL
*
*****************************************************************************/
#ifndef _NCBI_Seqport_
#include <seqport.h>
#endif

#ifndef _GATHER_
#include <gather.h>
#endif

#ifndef _SQNUTILS_
#include <sqnutils.h>
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

typedef void (*SpellCallBackFunc) (char * str);
typedef int (* SpellCheckFunc) (char *String, SpellCallBackFunc);

/* callback type for finer error reporting */

typedef void (LIBCALLBACK *ValidErrorFunc) (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
);

#define SET_DEPTH 20

#define VALIDATE_ALL 0
#define VALIDATE_INST 1
#define VALIDATE_HIST 2
#define VALIDATE_CONTEXT 3
#define VALIDATE_GRAPH 4
#define VALIDATE_SET 5
#define VALIDATE_FEAT 6
#define VALIDATE_DESC 7

typedef struct validstruct {
    Int2 cutoff;                   /* lowest errmsg to show 0=default */
    Int2 errors[6];
    SeqEntryPtr sep;               /* top level SeqEntryPtr */
    BioseqSetPtr bssp;               /* current bioseqset */
    BioseqPtr bsp;                 /* current bioseq */
    SeqFeatPtr sfp;                /* current feature */
    ValNodePtr descr;              /* current descriptor */
    Uint4 descrs [SET_DEPTH];      /* bit flags set by descriptor type */
    Int2 protcnt, nuccnt, segcnt;
    CharPtr errbuf;
    Boolean patch_seq;             /* repair invalid sequence residues? */
    Boolean non_ascii_chars;       /* non ascii chars found in read? */
    Boolean suppress_no_pubs;      /* internal use for no pub anywhere message */
    Boolean suppress_no_biosrc;    /* internal use for no biosource anywhere message */
    SpellCheckFunc spellfunc;
    SpellCallBackFunc spellcallback;
    GatherContextPtr gcp;          /* used for reporting the errors */
                                   /* this section used for checking Bioseqs */
    Uint2 bsp_partial_val;         /* return from SeqLocPartial on segmented SeqLoc */
    Boolean onlyspell;             /* only do spell check */
    Boolean justwarnonspell;       /* severity WARNING instead of ERROR on spell */
    Boolean useSeqMgrIndexes;      /* new style indexing to speed up validation */
    Boolean suppressContext;       /* suppress context part of message */
    Boolean validateAlignments;    /* call alignval test suite */
    Boolean farIDsInAlignments;    /* fetch to get far IDs in alignments */
    Boolean alignFindRemoteBsp;    /* do remote fetching in alignment validation */
    Boolean doSeqHistAssembly;     /* do alignment validation in Seq-hist.assembly */
    Boolean alwaysRequireIsoJTA;   /* force check for iso_jta */
    Boolean farFetchCDSproducts;   /* lock CDS->products for CdTransCheck, if necessary */
    Boolean farFetchMRNAproducts;  /* lock MRNA->products for MrnaTransCheck, if necessary */
    Boolean locusTagGeneralMatch;  /* expect locus_tag to match Seq-id.general of CDS and mRNA product */
    Boolean validateIDSet;         /* look for gain or loss of general IDs on sequence update */
    Boolean seqSubmitParent;       /* flag from tbl2asn to suppress no pub message */
    Boolean justShowAccession;     /* extremely terse output with accession and error type */
    Boolean ignoreExceptions;      /* report translation and transcription problems even if exception set */
    Boolean validateExons;         /* check splice sites on exon features except at ends of overlapping CDS */
    Boolean inferenceAccnCheck;    /* lookup inference qualifier accession.version reference */
    Boolean testLatLonSubregion;   /* validate coordinates of states and provinces within a country */
    Boolean strictLatLonCountry;   /* bodies of water do not relax country vs. lat_lon mismatch */
    Boolean indexerVersion;        /* special tests for GenBank indexers */
    Int2 validationLimit;          /* limit validation to major classes in Valid1GatherProc */
                                   /* this section used for finer error reporting callback */
    ValidErrorFunc errfunc;
    Pointer userdata;
    Boolean convertGiToAccn;
                                   /* this section used for internal flags */
    TextFsaPtr sourceQualTags;     /* for detecting structured qual tags in notes */
    TextFsaPtr modifiedBases;      /* permitted modified bases in PCR_primer qualifier */
    Boolean is_htg_in_sep;         /* record has technique of htgs 0 through htgs 3 */
    Boolean is_barcode_sep;        /* record has technique barcode */
    Boolean is_refseq_in_sep;      /* record has seqid of type other (refseq) */
    Boolean is_gps_in_sep;         /* record has genomic product set */
    Boolean is_embl_ddbj_in_sep;   /* record has embl or ddbj seqid */
    Boolean is_insd_in_sep;        /* record has genbank/embl/ddbj or tpg/tpe/tpd seqid */
    Boolean only_lcl_gnl_in_sep;   /* record has seqid of only local or general */
    Boolean has_gnl_prot_sep;      /* protein Bioseq has general seqid */
    Boolean is_smupd_in_sep;       /* record in INSD internal processing */
    Boolean feat_loc_has_gi;       /* at least one feature has a gi location reference */
    Boolean feat_prod_has_gi;      /* at least one feature has a gi product reference */
    Boolean far_fetch_failure;     /* a far location or bioseq with no fetch function */
    VoidPtr rrna_array;            /* sorted feature index array of rRNA features */
    VoidPtr trna_array;            /* sorted feature index array of tRNA features */
    Int4 numrrna;                  /* number of rRNA features */
    Int4 numtrna;                  /* number of tRNA features */
} ValidStruct, PNTR ValidStructPtr;

NLM_EXTERN Boolean ValidateSeqEntry PROTO((SeqEntryPtr sep, ValidStructPtr vsp));
NLM_EXTERN void ValidStructClear (ValidStructPtr vsp);  /* 0 out a ValidStruct */
NLM_EXTERN ValidStructPtr ValidStructNew (void);
NLM_EXTERN ValidStructPtr ValidStructFree (ValidStructPtr vsp);
NLM_EXTERN void SpellCallBack (char * str);
NLM_EXTERN Boolean IsNuclAcc (CharPtr name);

NLM_EXTERN CharPtr GetValidCategoryName (int errcode);
NLM_EXTERN CharPtr GetValidErrorName (int errcode, int subcode);
NLM_EXTERN CharPtr GetValidExplanation (int errcode, int subcode);

NLM_EXTERN CharPtr PNTR GetValidCountryList (void);
NLM_EXTERN Boolean CountryIsValid (CharPtr name, BoolPtr old_countryP, BoolPtr bad_capP);
NLM_EXTERN CharPtr GetCorrectedCountryCapitalization (CharPtr name);
NLM_EXTERN Boolean LookForECnumberPattern (CharPtr str);

NLM_EXTERN Boolean IsCountryInLatLonList (CharPtr country);
NLM_EXTERN Boolean TestLatLonForCountry (CharPtr country, FloatHi lat, FloatHi lon);
NLM_EXTERN CharPtr GuessCountryForLatLon (FloatHi lat, FloatHi lon);
NLM_EXTERN Boolean StringContainsBodyOfWater (CharPtr str);

NLM_EXTERN Boolean ParseStructuredVoucher (CharPtr subname, CharPtr PNTR inst, CharPtr PNTR id);
NLM_EXTERN Boolean VoucherInstitutionIsValid (CharPtr inst);

/* EC_number finite state machine persists to avoid expensive reload, should free on program exit */
NLM_EXTERN void ECNumberFSAFreeAll (void);

NLM_EXTERN Boolean HasTpaUserObject (BioseqPtr bsp);

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
