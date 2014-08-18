#ifndef _objinsdseq_ 
#define _objinsdseq_ 

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
*    Generated objects for Module INSD-INSDSeq
*    Generated using ASNCODE Revision: 6.16 at Jan 15, 2009  2:16 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objinsdseqAsnLoad PROTO((void));


/**************************************************
*
*    INSDSet
*
**************************************************/
typedef struct struct_INSDSeq INSDSet;
typedef struct struct_INSDSeq PNTR INSDSetPtr;
#define INSDSetNew() INSDSeqNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN INSDSetPtr LIBCALL INSDSetFree PROTO ((INSDSetPtr ));
NLM_EXTERN INSDSetPtr LIBCALL INSDSetNew PROTO (( void ));
NLM_EXTERN INSDSetPtr LIBCALL INSDSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDSetAsnWrite PROTO (( INSDSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    INSDSeq
*
**************************************************/
typedef struct struct_INSDSeq {
   struct struct_INSDSeq PNTR next;
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
   struct struct_INSDReference PNTR   references;
   CharPtr   comment;
   struct struct_INSDTagset PNTR   tagset;
   CharPtr   primary;
   CharPtr   source_db;
   CharPtr   database_reference;
   struct struct_INSDFeature PNTR   feature_table;
   CharPtr   sequence;
   CharPtr   contig;
} INSDSeq, PNTR INSDSeqPtr;


NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqFree PROTO ((INSDSeqPtr ));
NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqNew PROTO (( void ));
NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDSeqAsnWrite PROTO (( INSDSeqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDReference
*
**************************************************/
typedef struct struct_INSDReference {
   struct struct_INSDReference PNTR next;
   Uint4 OBbits__;
   CharPtr   reference;
   CharPtr   position;
   ValNodePtr   authors;
   CharPtr   consortium;
   CharPtr   title;
   CharPtr   journal;
   struct struct_INSDXref PNTR   xref;
#define OB__INSDReference_pubmed 0

   Int4   pubmed;
   CharPtr   remark;
} INSDReference, PNTR INSDReferencePtr;


NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceFree PROTO ((INSDReferencePtr ));
NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceNew PROTO (( void ));
NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDReferenceAsnWrite PROTO (( INSDReferencePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDTagset
*
**************************************************/
typedef struct struct_INSDTagset {
   Uint4 OBbits__;
   CharPtr   authority;
   CharPtr   version;
   CharPtr   url;
   struct struct_INSDTag PNTR   tags;
} INSDTagset, PNTR INSDTagsetPtr;


NLM_EXTERN INSDTagsetPtr LIBCALL INSDTagsetFree PROTO ((INSDTagsetPtr ));
NLM_EXTERN INSDTagsetPtr LIBCALL INSDTagsetNew PROTO (( void ));
NLM_EXTERN INSDTagsetPtr LIBCALL INSDTagsetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagsetAsnWrite PROTO (( INSDTagsetPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDFeature
*
**************************************************/
typedef struct struct_INSDFeature {
   struct struct_INSDFeature PNTR next;
   Uint4 OBbits__;
   CharPtr   key;
   CharPtr   location;
   struct struct_INSDInterval PNTR   intervals;
   CharPtr   operator__;
#define OB__INSDFeature_partial5 0

   Uint1   partial5;
#define OB__INSDFeature_partial3 1

   Uint1   partial3;
   struct struct_INSDQualifier PNTR   quals;
} INSDFeature, PNTR INSDFeaturePtr;


NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureFree PROTO ((INSDFeaturePtr ));
NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureNew PROTO (( void ));
NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDFeatureAsnWrite PROTO (( INSDFeaturePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDXref
*
**************************************************/
typedef struct struct_INSDXref {
   struct struct_INSDXref PNTR next;
   Uint4 OBbits__;
   CharPtr   dbname;
   CharPtr   id;
} INSDXref, PNTR INSDXrefPtr;


NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefFree PROTO ((INSDXrefPtr ));
NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefNew PROTO (( void ));
NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDXrefAsnWrite PROTO (( INSDXrefPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDTags
*
**************************************************/
typedef struct struct_INSDTag INSDTags;
typedef struct struct_INSDTag PNTR INSDTagsPtr;
#define INSDTagsNew() INSDTagNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN INSDTagsPtr LIBCALL INSDTagsFree PROTO ((INSDTagsPtr ));
NLM_EXTERN INSDTagsPtr LIBCALL INSDTagsNew PROTO (( void ));
NLM_EXTERN INSDTagsPtr LIBCALL INSDTagsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagsAsnWrite PROTO (( INSDTagsPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    INSDTag
*
**************************************************/
typedef struct struct_INSDTag {
   struct struct_INSDTag PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
   CharPtr   unit;
} INSDTag, PNTR INSDTagPtr;


NLM_EXTERN INSDTagPtr LIBCALL INSDTagFree PROTO ((INSDTagPtr ));
NLM_EXTERN INSDTagPtr LIBCALL INSDTagNew PROTO (( void ));
NLM_EXTERN INSDTagPtr LIBCALL INSDTagAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagAsnWrite PROTO (( INSDTagPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDInterval
*
**************************************************/
typedef struct struct_INSDInterval {
   struct struct_INSDInterval PNTR next;
   Uint4 OBbits__;
#define OB__INSDInterval_from 0

   Int4   from;
#define OB__INSDInterval_to 1

   Int4   to;
#define OB__INSDInterval_point 2

   Int4   point;
#define OB__INSDInterval_iscomp 3

   Uint1   iscomp;
#define OB__INSDInterval_interbp 4

   Uint1   interbp;
   CharPtr   accession;
} INSDInterval, PNTR INSDIntervalPtr;


NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalFree PROTO ((INSDIntervalPtr ));
NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalNew PROTO (( void ));
NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDIntervalAsnWrite PROTO (( INSDIntervalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDQualifier
*
**************************************************/
typedef struct struct_INSDQualifier {
   struct struct_INSDQualifier PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
} INSDQualifier, PNTR INSDQualifierPtr;


NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierFree PROTO ((INSDQualifierPtr ));
NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierNew PROTO (( void ));
NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDQualifierAsnWrite PROTO (( INSDQualifierPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDTagsetRules
*
**************************************************/
typedef struct struct_INSDTagsetRules {
   struct struct_INSDTagsetRules PNTR next;
   Uint4 OBbits__;
   CharPtr   authority;
   CharPtr   version;
   ValNodePtr   mandatorytags;
   ValNodePtr   optionaltags;
   ValNodePtr   uniquetags;
#define OB__INSDTagsetRules_extensible 0

   Uint1   extensible;
} INSDTagsetRules, PNTR INSDTagsetRulesPtr;


NLM_EXTERN INSDTagsetRulesPtr LIBCALL INSDTagsetRulesFree PROTO ((INSDTagsetRulesPtr ));
NLM_EXTERN INSDTagsetRulesPtr LIBCALL INSDTagsetRulesNew PROTO (( void ));
NLM_EXTERN INSDTagsetRulesPtr LIBCALL INSDTagsetRulesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagsetRulesAsnWrite PROTO (( INSDTagsetRulesPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDTagNames
*
**************************************************/
typedef ValNode INSDTagNames;
typedef ValNodePtr INSDTagNamesPtr;
#define INSDTagNamesNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN INSDTagNamesPtr LIBCALL INSDTagNamesFree PROTO ((INSDTagNamesPtr ));
NLM_EXTERN INSDTagNamesPtr LIBCALL INSDTagNamesNew PROTO (( void ));
NLM_EXTERN INSDTagNamesPtr LIBCALL INSDTagNamesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagNamesAsnWrite PROTO (( INSDTagNamesPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    INSDTagsetRuleSet
*
**************************************************/
typedef struct struct_INSDTagsetRules INSDTagsetRuleSet;
typedef struct struct_INSDTagsetRules PNTR INSDTagsetRuleSetPtr;
#define INSDTagsetRuleSetNew() INSDTagsetRulesNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN INSDTagsetRuleSetPtr LIBCALL INSDTagsetRuleSetFree PROTO ((INSDTagsetRuleSetPtr ));
NLM_EXTERN INSDTagsetRuleSetPtr LIBCALL INSDTagsetRuleSetNew PROTO (( void ));
NLM_EXTERN INSDTagsetRuleSetPtr LIBCALL INSDTagsetRuleSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDTagsetRuleSetAsnWrite PROTO (( INSDTagsetRuleSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objinsdseq_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

