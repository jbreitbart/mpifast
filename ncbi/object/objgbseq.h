#ifndef _objgbseq_ 
#define _objgbseq_ 

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module NCBI-GBSeq
*    Generated using ASNCODE Revision: 6.16 at Jan 15, 2009  2:16 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objgbseqAsnLoad PROTO((void));


/**************************************************
*
*    GBSet
*
**************************************************/
typedef struct struct_GBSeq GBSet;
typedef struct struct_GBSeq PNTR GBSetPtr;
#define GBSetNew() GBSeqNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN GBSetPtr LIBCALL GBSetFree PROTO ((GBSetPtr ));
NLM_EXTERN GBSetPtr LIBCALL GBSetNew PROTO (( void ));
NLM_EXTERN GBSetPtr LIBCALL GBSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBSetAsnWrite PROTO (( GBSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    GBSeq
*
**************************************************/
typedef struct struct_GBSeq {
   struct struct_GBSeq PNTR next;
   Uint4 OBbits__;
   CharPtr   locus;
   Int4   length;
   CharPtr   strandedness;
   CharPtr   moltype;
   CharPtr   topology;
   CharPtr   division;
   CharPtr   update_date;
   CharPtr   create_date;
   CharPtr   update_release;
   CharPtr   create_release;
   CharPtr   definition;
   CharPtr   primary_accession;
   CharPtr   entry_version;
   CharPtr   accession_version;
   ValNodePtr   other_seqids;
   ValNodePtr   secondary_accessions;
   CharPtr   project;
   ValNodePtr   keywords;
   CharPtr   segment;
   CharPtr   source;
   CharPtr   organism;
   CharPtr   taxonomy;
   struct struct_GBReference PNTR   references;
   CharPtr   comment;
   struct struct_GBTagset PNTR   tagset;
   CharPtr   primary;
   CharPtr   source_db;
   CharPtr   database_reference;
   struct struct_GBFeature PNTR   feature_table;
   CharPtr   sequence;
   CharPtr   contig;
} GBSeq, PNTR GBSeqPtr;


NLM_EXTERN GBSeqPtr LIBCALL GBSeqFree PROTO ((GBSeqPtr ));
NLM_EXTERN GBSeqPtr LIBCALL GBSeqNew PROTO (( void ));
NLM_EXTERN GBSeqPtr LIBCALL GBSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBSeqAsnWrite PROTO (( GBSeqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBReference
*
**************************************************/
typedef struct struct_GBReference {
   struct struct_GBReference PNTR next;
   Uint4 OBbits__;
   CharPtr   reference;
   CharPtr   position;
   ValNodePtr   authors;
   CharPtr   consortium;
   CharPtr   title;
   CharPtr   journal;
   struct struct_GBXref PNTR   xref;
#define OB__GBReference_pubmed 0

   Int4   pubmed;
   CharPtr   remark;
} GBReference, PNTR GBReferencePtr;


NLM_EXTERN GBReferencePtr LIBCALL GBReferenceFree PROTO ((GBReferencePtr ));
NLM_EXTERN GBReferencePtr LIBCALL GBReferenceNew PROTO (( void ));
NLM_EXTERN GBReferencePtr LIBCALL GBReferenceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBReferenceAsnWrite PROTO (( GBReferencePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBTagset
*
**************************************************/
typedef struct struct_GBTagset {
   Uint4 OBbits__;
   CharPtr   authority;
   CharPtr   version;
   CharPtr   url;
   struct struct_GBTag PNTR   tags;
} GBTagset, PNTR GBTagsetPtr;


NLM_EXTERN GBTagsetPtr LIBCALL GBTagsetFree PROTO ((GBTagsetPtr ));
NLM_EXTERN GBTagsetPtr LIBCALL GBTagsetNew PROTO (( void ));
NLM_EXTERN GBTagsetPtr LIBCALL GBTagsetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagsetAsnWrite PROTO (( GBTagsetPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBFeature
*
**************************************************/
typedef struct struct_GBFeature {
   struct struct_GBFeature PNTR next;
   Uint4 OBbits__;
   CharPtr   key;
   CharPtr   location;
   struct struct_GBInterval PNTR   intervals;
   CharPtr   operator__;
#define OB__GBFeature_partial5 0

   Uint1   partial5;
#define OB__GBFeature_partial3 1

   Uint1   partial3;
   struct struct_GBQualifier PNTR   quals;
} GBFeature, PNTR GBFeaturePtr;


NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureFree PROTO ((GBFeaturePtr ));
NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureNew PROTO (( void ));
NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBFeatureAsnWrite PROTO (( GBFeaturePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBXref
*
**************************************************/
typedef struct struct_GBXref {
   struct struct_GBXref PNTR next;
   Uint4 OBbits__;
   CharPtr   dbname;
   CharPtr   id;
} GBXref, PNTR GBXrefPtr;


NLM_EXTERN GBXrefPtr LIBCALL GBXrefFree PROTO ((GBXrefPtr ));
NLM_EXTERN GBXrefPtr LIBCALL GBXrefNew PROTO (( void ));
NLM_EXTERN GBXrefPtr LIBCALL GBXrefAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBXrefAsnWrite PROTO (( GBXrefPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBTags
*
**************************************************/
typedef struct struct_GBTag GBTags;
typedef struct struct_GBTag PNTR GBTagsPtr;
#define GBTagsNew() GBTagNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN GBTagsPtr LIBCALL GBTagsFree PROTO ((GBTagsPtr ));
NLM_EXTERN GBTagsPtr LIBCALL GBTagsNew PROTO (( void ));
NLM_EXTERN GBTagsPtr LIBCALL GBTagsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagsAsnWrite PROTO (( GBTagsPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    GBTag
*
**************************************************/
typedef struct struct_GBTag {
   struct struct_GBTag PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
   CharPtr   unit;
} GBTag, PNTR GBTagPtr;


NLM_EXTERN GBTagPtr LIBCALL GBTagFree PROTO ((GBTagPtr ));
NLM_EXTERN GBTagPtr LIBCALL GBTagNew PROTO (( void ));
NLM_EXTERN GBTagPtr LIBCALL GBTagAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagAsnWrite PROTO (( GBTagPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBInterval
*
**************************************************/
typedef struct struct_GBInterval {
   struct struct_GBInterval PNTR next;
   Uint4 OBbits__;
#define OB__GBInterval_from 0

   Int4   from;
#define OB__GBInterval_to 1

   Int4   to;
#define OB__GBInterval_point 2

   Int4   point;
#define OB__GBInterval_iscomp 3

   Uint1   iscomp;
#define OB__GBInterval_interbp 4

   Uint1   interbp;
   CharPtr   accession;
} GBInterval, PNTR GBIntervalPtr;


NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalFree PROTO ((GBIntervalPtr ));
NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalNew PROTO (( void ));
NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBIntervalAsnWrite PROTO (( GBIntervalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBQualifier
*
**************************************************/
typedef struct struct_GBQualifier {
   struct struct_GBQualifier PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
} GBQualifier, PNTR GBQualifierPtr;


NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierFree PROTO ((GBQualifierPtr ));
NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierNew PROTO (( void ));
NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBQualifierAsnWrite PROTO (( GBQualifierPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBTagsetRules
*
**************************************************/
typedef struct struct_GBTagsetRules {
   struct struct_GBTagsetRules PNTR next;
   Uint4 OBbits__;
   CharPtr   authority;
   CharPtr   version;
   ValNodePtr   mandatorytags;
   ValNodePtr   optionaltags;
   ValNodePtr   uniquetags;
#define OB__GBTagsetRules_extensible 0

   Uint1   extensible;
} GBTagsetRules, PNTR GBTagsetRulesPtr;


NLM_EXTERN GBTagsetRulesPtr LIBCALL GBTagsetRulesFree PROTO ((GBTagsetRulesPtr ));
NLM_EXTERN GBTagsetRulesPtr LIBCALL GBTagsetRulesNew PROTO (( void ));
NLM_EXTERN GBTagsetRulesPtr LIBCALL GBTagsetRulesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagsetRulesAsnWrite PROTO (( GBTagsetRulesPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBTagNames
*
**************************************************/
typedef ValNode GBTagNames;
typedef ValNodePtr GBTagNamesPtr;
#define GBTagNamesNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN GBTagNamesPtr LIBCALL GBTagNamesFree PROTO ((GBTagNamesPtr ));
NLM_EXTERN GBTagNamesPtr LIBCALL GBTagNamesNew PROTO (( void ));
NLM_EXTERN GBTagNamesPtr LIBCALL GBTagNamesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagNamesAsnWrite PROTO (( GBTagNamesPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    GBTagsetRuleSet
*
**************************************************/
typedef struct struct_GBTagsetRules GBTagsetRuleSet;
typedef struct struct_GBTagsetRules PNTR GBTagsetRuleSetPtr;
#define GBTagsetRuleSetNew() GBTagsetRulesNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN GBTagsetRuleSetPtr LIBCALL GBTagsetRuleSetFree PROTO ((GBTagsetRuleSetPtr ));
NLM_EXTERN GBTagsetRuleSetPtr LIBCALL GBTagsetRuleSetNew PROTO (( void ));
NLM_EXTERN GBTagsetRuleSetPtr LIBCALL GBTagsetRuleSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBTagsetRuleSetAsnWrite PROTO (( GBTagsetRuleSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objgbseq_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

