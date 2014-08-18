#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objgbseq.h>

static Boolean loaded = FALSE;

#include <asngbseq.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objgbseqAsnLoad(void)
{

   if ( ! loaded) {
      NLM_EXTERN_LOADS

      if ( ! AsnLoad ())
      return FALSE;
      loaded = TRUE;
   }

   return TRUE;
}



/**************************************************
*    Generated object loaders for Module NCBI-GBSeq
*    Generated using ASNCODE Revision: 6.16 at Jan 15, 2009  2:16 PM
*
**************************************************/


/**************************************************
*
*    GBSetFree()
*
**************************************************/
NLM_EXTERN 
GBSetPtr LIBCALL
GBSetFree(GBSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) GBSeqFree);
   return NULL;
}


/**************************************************
*
*    GBSetAsnRead()
*
**************************************************/
NLM_EXTERN 
GBSetPtr LIBCALL
GBSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBSetPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBSet ::= (self contained) */
      atp = AsnReadId(aip, amp, GBSET);
   } else {
      atp = AsnLinkType(orig, GBSET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBSeqAsnRead, (AsnOptFreeFunc) GBSeqFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBSetAsnWrite(GBSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBSET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) GBSeqAsnWrite, aip, atp, GBSET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBSeqNew()
*
**************************************************/
NLM_EXTERN 
GBSeqPtr LIBCALL
GBSeqNew(void)
{
   GBSeqPtr ptr = MemNew((size_t) sizeof(GBSeq));

   return ptr;

}


/**************************************************
*
*    GBSeqFree()
*
**************************************************/
NLM_EXTERN 
GBSeqPtr LIBCALL
GBSeqFree(GBSeqPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> locus);
   MemFree(ptr -> strandedness);
   MemFree(ptr -> moltype);
   MemFree(ptr -> topology);
   MemFree(ptr -> division);
   MemFree(ptr -> update_date);
   MemFree(ptr -> create_date);
   MemFree(ptr -> update_release);
   MemFree(ptr -> create_release);
   MemFree(ptr -> definition);
   MemFree(ptr -> primary_accession);
   MemFree(ptr -> entry_version);
   MemFree(ptr -> accession_version);
   AsnGenericBaseSeqOfFree(ptr -> other_seqids ,ASNCODE_PTRVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> secondary_accessions ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> project);
   AsnGenericBaseSeqOfFree(ptr -> keywords ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> segment);
   MemFree(ptr -> source);
   MemFree(ptr -> organism);
   MemFree(ptr -> taxonomy);
   AsnGenericUserSeqOfFree(ptr -> references, (AsnOptFreeFunc) GBReferenceFree);
   MemFree(ptr -> comment);
   GBTagsetFree(ptr -> tagset);
   MemFree(ptr -> primary);
   MemFree(ptr -> source_db);
   MemFree(ptr -> database_reference);
   AsnGenericUserSeqOfFree(ptr -> feature_table, (AsnOptFreeFunc) GBFeatureFree);
   MemFree(ptr -> sequence);
   MemFree(ptr -> contig);
   return MemFree(ptr);
}


/**************************************************
*
*    GBSeqAsnRead()
*
**************************************************/
NLM_EXTERN 
GBSeqPtr LIBCALL
GBSeqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBSeqPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBSeq ::= (self contained) */
      atp = AsnReadId(aip, amp, GBSEQ);
   } else {
      atp = AsnLinkType(orig, GBSEQ);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBSeqNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBSEQ_locus) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> locus = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> length = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_strandedness) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strandedness = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_moltype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> moltype = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_topology) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> topology = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_division) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> division = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_update_date) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> update_date = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_create_date) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> create_date = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_update_release) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> update_release = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_create_release) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> create_release = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_definition) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> definition = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_primary_accession) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> primary_accession = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_entry_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> entry_version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_accession_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> accession_version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_other_seqids) {
      ptr -> other_seqids = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> other_seqids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_secondary_accessions) {
      ptr -> secondary_accessions = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> secondary_accessions == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_project) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> project = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_keywords) {
      ptr -> keywords = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> keywords == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_segment) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> segment = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_organism) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> organism = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_taxonomy) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> taxonomy = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_references) {
      ptr -> references = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBReferenceAsnRead, (AsnOptFreeFunc) GBReferenceFree);
      if (isError && ptr -> references == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_comment) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> comment = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_tagset) {
      ptr -> tagset = GBTagsetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_primary) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> primary = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_source_db) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source_db = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_database_reference) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> database_reference = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_feature_table) {
      ptr -> feature_table = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBFeatureAsnRead, (AsnOptFreeFunc) GBFeatureFree);
      if (isError && ptr -> feature_table == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_sequence) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> sequence = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBSEQ_contig) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> contig = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBSeqFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBSeqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBSeqAsnWrite(GBSeqPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBSEQ);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> locus != NULL) {
      av.ptrvalue = ptr -> locus;
      retval = AsnWrite(aip, GBSEQ_locus,  &av);
   }
   av.intvalue = ptr -> length;
   retval = AsnWrite(aip, GBSEQ_length,  &av);
   if (ptr -> strandedness != NULL) {
      av.ptrvalue = ptr -> strandedness;
      retval = AsnWrite(aip, GBSEQ_strandedness,  &av);
   }
   if (ptr -> moltype != NULL) {
      av.ptrvalue = ptr -> moltype;
      retval = AsnWrite(aip, GBSEQ_moltype,  &av);
   }
   if (ptr -> topology != NULL) {
      av.ptrvalue = ptr -> topology;
      retval = AsnWrite(aip, GBSEQ_topology,  &av);
   }
   if (ptr -> division != NULL) {
      av.ptrvalue = ptr -> division;
      retval = AsnWrite(aip, GBSEQ_division,  &av);
   }
   if (ptr -> update_date != NULL) {
      av.ptrvalue = ptr -> update_date;
      retval = AsnWrite(aip, GBSEQ_update_date,  &av);
   }
   if (ptr -> create_date != NULL) {
      av.ptrvalue = ptr -> create_date;
      retval = AsnWrite(aip, GBSEQ_create_date,  &av);
   }
   if (ptr -> update_release != NULL) {
      av.ptrvalue = ptr -> update_release;
      retval = AsnWrite(aip, GBSEQ_update_release,  &av);
   }
   if (ptr -> create_release != NULL) {
      av.ptrvalue = ptr -> create_release;
      retval = AsnWrite(aip, GBSEQ_create_release,  &av);
   }
   if (ptr -> definition != NULL) {
      av.ptrvalue = ptr -> definition;
      retval = AsnWrite(aip, GBSEQ_definition,  &av);
   }
   if (ptr -> primary_accession != NULL) {
      av.ptrvalue = ptr -> primary_accession;
      retval = AsnWrite(aip, GBSEQ_primary_accession,  &av);
   }
   if (ptr -> entry_version != NULL) {
      av.ptrvalue = ptr -> entry_version;
      retval = AsnWrite(aip, GBSEQ_entry_version,  &av);
   }
   if (ptr -> accession_version != NULL) {
      av.ptrvalue = ptr -> accession_version;
      retval = AsnWrite(aip, GBSEQ_accession_version,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> other_seqids ,ASNCODE_PTRVAL_SLOT, aip, GBSEQ_other_seqids, GBSEQ_other_seqids_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> secondary_accessions ,ASNCODE_PTRVAL_SLOT, aip, GBSEQ_secondary_accessions, GBSEQ_secondary_accessions_E);
   if (ptr -> project != NULL) {
      av.ptrvalue = ptr -> project;
      retval = AsnWrite(aip, GBSEQ_project,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> keywords ,ASNCODE_PTRVAL_SLOT, aip, GBSEQ_keywords, GBSEQ_keywords_E);
   if (ptr -> segment != NULL) {
      av.ptrvalue = ptr -> segment;
      retval = AsnWrite(aip, GBSEQ_segment,  &av);
   }
   if (ptr -> source != NULL) {
      av.ptrvalue = ptr -> source;
      retval = AsnWrite(aip, GBSEQ_source,  &av);
   }
   if (ptr -> organism != NULL) {
      av.ptrvalue = ptr -> organism;
      retval = AsnWrite(aip, GBSEQ_organism,  &av);
   }
   if (ptr -> taxonomy != NULL) {
      av.ptrvalue = ptr -> taxonomy;
      retval = AsnWrite(aip, GBSEQ_taxonomy,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> references, (AsnWriteFunc) GBReferenceAsnWrite, aip, GBSEQ_references, GBSEQ_references_E);
   if (ptr -> comment != NULL) {
      av.ptrvalue = ptr -> comment;
      retval = AsnWrite(aip, GBSEQ_comment,  &av);
   }
   if (ptr -> tagset != NULL) {
      if ( ! GBTagsetAsnWrite(ptr -> tagset, aip, GBSEQ_tagset)) {
         goto erret;
      }
   }
   if (ptr -> primary != NULL) {
      av.ptrvalue = ptr -> primary;
      retval = AsnWrite(aip, GBSEQ_primary,  &av);
   }
   if (ptr -> source_db != NULL) {
      av.ptrvalue = ptr -> source_db;
      retval = AsnWrite(aip, GBSEQ_source_db,  &av);
   }
   if (ptr -> database_reference != NULL) {
      av.ptrvalue = ptr -> database_reference;
      retval = AsnWrite(aip, GBSEQ_database_reference,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> feature_table, (AsnWriteFunc) GBFeatureAsnWrite, aip, GBSEQ_feature_table, GBSEQ_feature_table_E);
   if (ptr -> sequence != NULL) {
      av.ptrvalue = ptr -> sequence;
      retval = AsnWrite(aip, GBSEQ_sequence,  &av);
   }
   if (ptr -> contig != NULL) {
      av.ptrvalue = ptr -> contig;
      retval = AsnWrite(aip, GBSEQ_contig,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBReferenceNew()
*
**************************************************/
NLM_EXTERN 
GBReferencePtr LIBCALL
GBReferenceNew(void)
{
   GBReferencePtr ptr = MemNew((size_t) sizeof(GBReference));

   return ptr;

}


/**************************************************
*
*    GBReferenceFree()
*
**************************************************/
NLM_EXTERN 
GBReferencePtr LIBCALL
GBReferenceFree(GBReferencePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> reference);
   MemFree(ptr -> position);
   AsnGenericBaseSeqOfFree(ptr -> authors ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> consortium);
   MemFree(ptr -> title);
   MemFree(ptr -> journal);
   AsnGenericUserSeqOfFree(ptr -> xref, (AsnOptFreeFunc) GBXrefFree);
   MemFree(ptr -> remark);
   return MemFree(ptr);
}


/**************************************************
*
*    GBReferenceAsnRead()
*
**************************************************/
NLM_EXTERN 
GBReferencePtr LIBCALL
GBReferenceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBReferencePtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBReference ::= (self contained) */
      atp = AsnReadId(aip, amp, GBREFERENCE);
   } else {
      atp = AsnLinkType(orig, GBREFERENCE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBReferenceNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBREFERENCE_reference) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> reference = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_position) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> position = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_authors) {
      ptr -> authors = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> authors == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_consortium) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> consortium = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_title) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> title = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_journal) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> journal = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_xref) {
      ptr -> xref = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBXrefAsnRead, (AsnOptFreeFunc) GBXrefFree);
      if (isError && ptr -> xref == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_pubmed) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> pubmed = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBREFERENCE_remark) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remark = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBReferenceFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBReferenceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBReferenceAsnWrite(GBReferencePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBREFERENCE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> reference != NULL) {
      av.ptrvalue = ptr -> reference;
      retval = AsnWrite(aip, GBREFERENCE_reference,  &av);
   }
   if (ptr -> position != NULL) {
      av.ptrvalue = ptr -> position;
      retval = AsnWrite(aip, GBREFERENCE_position,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> authors ,ASNCODE_PTRVAL_SLOT, aip, GBREFERENCE_authors, GBREFERENCE_authors_E);
   if (ptr -> consortium != NULL) {
      av.ptrvalue = ptr -> consortium;
      retval = AsnWrite(aip, GBREFERENCE_consortium,  &av);
   }
   if (ptr -> title != NULL) {
      av.ptrvalue = ptr -> title;
      retval = AsnWrite(aip, GBREFERENCE_title,  &av);
   }
   if (ptr -> journal != NULL) {
      av.ptrvalue = ptr -> journal;
      retval = AsnWrite(aip, GBREFERENCE_journal,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> xref, (AsnWriteFunc) GBXrefAsnWrite, aip, GBREFERENCE_xref, GBREFERENCE_xref_E);
   if (ptr -> pubmed || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> pubmed;
      retval = AsnWrite(aip, GBREFERENCE_pubmed,  &av);
   }
   if (ptr -> remark != NULL) {
      av.ptrvalue = ptr -> remark;
      retval = AsnWrite(aip, GBREFERENCE_remark,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagsetNew()
*
**************************************************/
NLM_EXTERN 
GBTagsetPtr LIBCALL
GBTagsetNew(void)
{
   GBTagsetPtr ptr = MemNew((size_t) sizeof(GBTagset));

   return ptr;

}


/**************************************************
*
*    GBTagsetFree()
*
**************************************************/
NLM_EXTERN 
GBTagsetPtr LIBCALL
GBTagsetFree(GBTagsetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> authority);
   MemFree(ptr -> version);
   MemFree(ptr -> url);
   GBTagsFree(ptr -> tags);
   return MemFree(ptr);
}


/**************************************************
*
*    GBTagsetAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagsetPtr LIBCALL
GBTagsetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagsetPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTagset ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAGSET);
   } else {
      atp = AsnLinkType(orig, GBTAGSET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBTagsetNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBTAGSET_authority) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> authority = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSET_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSET_url) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> url = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSET_tags) {
      ptr -> tags = GBTagsAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagsetFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagsetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagsetAsnWrite(GBTagsetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAGSET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> authority != NULL) {
      av.ptrvalue = ptr -> authority;
      retval = AsnWrite(aip, GBTAGSET_authority,  &av);
   }
   if (ptr -> version != NULL) {
      av.ptrvalue = ptr -> version;
      retval = AsnWrite(aip, GBTAGSET_version,  &av);
   }
   if (ptr -> url != NULL) {
      av.ptrvalue = ptr -> url;
      retval = AsnWrite(aip, GBTAGSET_url,  &av);
   }
   if (ptr -> tags != NULL) {
      if ( ! GBTagsAsnWrite(ptr -> tags, aip, GBTAGSET_tags)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBFeatureNew()
*
**************************************************/
NLM_EXTERN 
GBFeaturePtr LIBCALL
GBFeatureNew(void)
{
   GBFeaturePtr ptr = MemNew((size_t) sizeof(GBFeature));

   return ptr;

}


/**************************************************
*
*    GBFeatureFree()
*
**************************************************/
NLM_EXTERN 
GBFeaturePtr LIBCALL
GBFeatureFree(GBFeaturePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> key);
   MemFree(ptr -> location);
   AsnGenericUserSeqOfFree(ptr -> intervals, (AsnOptFreeFunc) GBIntervalFree);
   MemFree(ptr -> operator__);
   AsnGenericUserSeqOfFree(ptr -> quals, (AsnOptFreeFunc) GBQualifierFree);
   return MemFree(ptr);
}


/**************************************************
*
*    GBFeatureAsnRead()
*
**************************************************/
NLM_EXTERN 
GBFeaturePtr LIBCALL
GBFeatureAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBFeaturePtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBFeature ::= (self contained) */
      atp = AsnReadId(aip, amp, GBFEATURE);
   } else {
      atp = AsnLinkType(orig, GBFEATURE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBFeatureNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBFEATURE_key) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> key = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> location = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_intervals) {
      ptr -> intervals = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBIntervalAsnRead, (AsnOptFreeFunc) GBIntervalFree);
      if (isError && ptr -> intervals == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_operator) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> operator__ = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_partial5) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial5 = av.boolvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_partial3) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial3 = av.boolvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBFEATURE_quals) {
      ptr -> quals = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBQualifierAsnRead, (AsnOptFreeFunc) GBQualifierFree);
      if (isError && ptr -> quals == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBFeatureFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBFeatureAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBFeatureAsnWrite(GBFeaturePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBFEATURE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> key != NULL) {
      av.ptrvalue = ptr -> key;
      retval = AsnWrite(aip, GBFEATURE_key,  &av);
   }
   if (ptr -> location != NULL) {
      av.ptrvalue = ptr -> location;
      retval = AsnWrite(aip, GBFEATURE_location,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> intervals, (AsnWriteFunc) GBIntervalAsnWrite, aip, GBFEATURE_intervals, GBFEATURE_intervals_E);
   if (ptr -> operator__ != NULL) {
      av.ptrvalue = ptr -> operator__;
      retval = AsnWrite(aip, GBFEATURE_operator,  &av);
   }
   if (ptr -> partial5 || (ptr -> OBbits__ & (1<<0) )){   av.boolvalue = ptr -> partial5;
      retval = AsnWrite(aip, GBFEATURE_partial5,  &av);
   }
   if (ptr -> partial3 || (ptr -> OBbits__ & (1<<1) )){   av.boolvalue = ptr -> partial3;
      retval = AsnWrite(aip, GBFEATURE_partial3,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> quals, (AsnWriteFunc) GBQualifierAsnWrite, aip, GBFEATURE_quals, GBFEATURE_quals_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBXrefNew()
*
**************************************************/
NLM_EXTERN 
GBXrefPtr LIBCALL
GBXrefNew(void)
{
   GBXrefPtr ptr = MemNew((size_t) sizeof(GBXref));

   return ptr;

}


/**************************************************
*
*    GBXrefFree()
*
**************************************************/
NLM_EXTERN 
GBXrefPtr LIBCALL
GBXrefFree(GBXrefPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> dbname);
   MemFree(ptr -> id);
   return MemFree(ptr);
}


/**************************************************
*
*    GBXrefAsnRead()
*
**************************************************/
NLM_EXTERN 
GBXrefPtr LIBCALL
GBXrefAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBXrefPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBXref ::= (self contained) */
      atp = AsnReadId(aip, amp, GBXREF);
   } else {
      atp = AsnLinkType(orig, GBXREF);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBXrefNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBXREF_dbname) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> dbname = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBXREF_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> id = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBXrefFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBXrefAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBXrefAsnWrite(GBXrefPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBXREF);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> dbname != NULL) {
      av.ptrvalue = ptr -> dbname;
      retval = AsnWrite(aip, GBXREF_dbname,  &av);
   }
   if (ptr -> id != NULL) {
      av.ptrvalue = ptr -> id;
      retval = AsnWrite(aip, GBXREF_id,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagsFree()
*
**************************************************/
NLM_EXTERN 
GBTagsPtr LIBCALL
GBTagsFree(GBTagsPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) GBTagFree);
   return NULL;
}


/**************************************************
*
*    GBTagsAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagsPtr LIBCALL
GBTagsAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagsPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTags ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAGS);
   } else {
      atp = AsnLinkType(orig, GBTAGS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBTagAsnRead, (AsnOptFreeFunc) GBTagFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagsFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagsAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagsAsnWrite(GBTagsPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAGS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) GBTagAsnWrite, aip, atp, GBTAGS_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagNew()
*
**************************************************/
NLM_EXTERN 
GBTagPtr LIBCALL
GBTagNew(void)
{
   GBTagPtr ptr = MemNew((size_t) sizeof(GBTag));

   return ptr;

}


/**************************************************
*
*    GBTagFree()
*
**************************************************/
NLM_EXTERN 
GBTagPtr LIBCALL
GBTagFree(GBTagPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   MemFree(ptr -> value);
   MemFree(ptr -> unit);
   return MemFree(ptr);
}


/**************************************************
*
*    GBTagAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagPtr LIBCALL
GBTagAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTag ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAG);
   } else {
      atp = AsnLinkType(orig, GBTAG);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBTagNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBTAG_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAG_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAG_unit) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> unit = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagAsnWrite(GBTagPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAG);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, GBTAG_name,  &av);
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, GBTAG_value,  &av);
   }
   if (ptr -> unit != NULL) {
      av.ptrvalue = ptr -> unit;
      retval = AsnWrite(aip, GBTAG_unit,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBIntervalNew()
*
**************************************************/
NLM_EXTERN 
GBIntervalPtr LIBCALL
GBIntervalNew(void)
{
   GBIntervalPtr ptr = MemNew((size_t) sizeof(GBInterval));

   return ptr;

}


/**************************************************
*
*    GBIntervalFree()
*
**************************************************/
NLM_EXTERN 
GBIntervalPtr LIBCALL
GBIntervalFree(GBIntervalPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> accession);
   return MemFree(ptr);
}


/**************************************************
*
*    GBIntervalAsnRead()
*
**************************************************/
NLM_EXTERN 
GBIntervalPtr LIBCALL
GBIntervalAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBIntervalPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBInterval ::= (self contained) */
      atp = AsnReadId(aip, amp, GBINTERVAL);
   } else {
      atp = AsnLinkType(orig, GBINTERVAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBIntervalNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBINTERVAL_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBINTERVAL_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBINTERVAL_point) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> point = av.intvalue;
      ptr -> OBbits__ |= 1<<2;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBINTERVAL_iscomp) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> iscomp = av.boolvalue;
      ptr -> OBbits__ |= 1<<3;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBINTERVAL_interbp) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> interbp = av.boolvalue;
      ptr -> OBbits__ |= 1<<4;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBINTERVAL_accession) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> accession = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBIntervalFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBIntervalAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBIntervalAsnWrite(GBIntervalPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBINTERVAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> from || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> from;
      retval = AsnWrite(aip, GBINTERVAL_from,  &av);
   }
   if (ptr -> to || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> to;
      retval = AsnWrite(aip, GBINTERVAL_to,  &av);
   }
   if (ptr -> point || (ptr -> OBbits__ & (1<<2) )){   av.intvalue = ptr -> point;
      retval = AsnWrite(aip, GBINTERVAL_point,  &av);
   }
   if (ptr -> iscomp || (ptr -> OBbits__ & (1<<3) )){   av.boolvalue = ptr -> iscomp;
      retval = AsnWrite(aip, GBINTERVAL_iscomp,  &av);
   }
   if (ptr -> interbp || (ptr -> OBbits__ & (1<<4) )){   av.boolvalue = ptr -> interbp;
      retval = AsnWrite(aip, GBINTERVAL_interbp,  &av);
   }
   if (ptr -> accession != NULL) {
      av.ptrvalue = ptr -> accession;
      retval = AsnWrite(aip, GBINTERVAL_accession,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBQualifierNew()
*
**************************************************/
NLM_EXTERN 
GBQualifierPtr LIBCALL
GBQualifierNew(void)
{
   GBQualifierPtr ptr = MemNew((size_t) sizeof(GBQualifier));

   return ptr;

}


/**************************************************
*
*    GBQualifierFree()
*
**************************************************/
NLM_EXTERN 
GBQualifierPtr LIBCALL
GBQualifierFree(GBQualifierPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   MemFree(ptr -> value);
   return MemFree(ptr);
}


/**************************************************
*
*    GBQualifierAsnRead()
*
**************************************************/
NLM_EXTERN 
GBQualifierPtr LIBCALL
GBQualifierAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBQualifierPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBQualifier ::= (self contained) */
      atp = AsnReadId(aip, amp, GBQUALIFIER);
   } else {
      atp = AsnLinkType(orig, GBQUALIFIER);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBQualifierNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBQUALIFIER_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBQUALIFIER_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBQualifierFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBQualifierAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBQualifierAsnWrite(GBQualifierPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBQUALIFIER);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, GBQUALIFIER_name,  &av);
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, GBQUALIFIER_value,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagsetRulesNew()
*
**************************************************/
NLM_EXTERN 
GBTagsetRulesPtr LIBCALL
GBTagsetRulesNew(void)
{
   GBTagsetRulesPtr ptr = MemNew((size_t) sizeof(GBTagsetRules));

   return ptr;

}


/**************************************************
*
*    GBTagsetRulesFree()
*
**************************************************/
NLM_EXTERN 
GBTagsetRulesPtr LIBCALL
GBTagsetRulesFree(GBTagsetRulesPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> authority);
   MemFree(ptr -> version);
   GBTagNamesFree(ptr -> mandatorytags);
   GBTagNamesFree(ptr -> optionaltags);
   GBTagNamesFree(ptr -> uniquetags);
   return MemFree(ptr);
}


/**************************************************
*
*    GBTagsetRulesAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagsetRulesPtr LIBCALL
GBTagsetRulesAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagsetRulesPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTagsetRules ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAGSETRULES);
   } else {
      atp = AsnLinkType(orig, GBTAGSETRULES);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = GBTagsetRulesNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GBTAGSETRULES_authority) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> authority = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSETRULES_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSETRULES_mandatorytags) {
      ptr -> mandatorytags = GBTagNamesAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSETRULES_optionaltags) {
      ptr -> optionaltags = GBTagNamesAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSETRULES_uniquetags) {
      ptr -> uniquetags = GBTagNamesAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GBTAGSETRULES_extensible) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> extensible = av.boolvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagsetRulesFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagsetRulesAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagsetRulesAsnWrite(GBTagsetRulesPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAGSETRULES);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> authority != NULL) {
      av.ptrvalue = ptr -> authority;
      retval = AsnWrite(aip, GBTAGSETRULES_authority,  &av);
   }
   if (ptr -> version != NULL) {
      av.ptrvalue = ptr -> version;
      retval = AsnWrite(aip, GBTAGSETRULES_version,  &av);
   }
   if (ptr -> mandatorytags != NULL) {
      if ( ! GBTagNamesAsnWrite(ptr -> mandatorytags, aip, GBTAGSETRULES_mandatorytags)) {
         goto erret;
      }
   }
   if (ptr -> optionaltags != NULL) {
      if ( ! GBTagNamesAsnWrite(ptr -> optionaltags, aip, GBTAGSETRULES_optionaltags)) {
         goto erret;
      }
   }
   if (ptr -> uniquetags != NULL) {
      if ( ! GBTagNamesAsnWrite(ptr -> uniquetags, aip, GBTAGSETRULES_uniquetags)) {
         goto erret;
      }
   }
   if (ptr -> extensible || (ptr -> OBbits__ & (1<<0) )){   av.boolvalue = ptr -> extensible;
      retval = AsnWrite(aip, GBTAGSETRULES_extensible,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagNamesFree()
*
**************************************************/
NLM_EXTERN 
GBTagNamesPtr LIBCALL
GBTagNamesFree(GBTagNamesPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr,ASNCODE_PTRVAL_SLOT);
   return NULL;
}


/**************************************************
*
*    GBTagNamesAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagNamesPtr LIBCALL
GBTagNamesAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagNamesPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTagNames ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAGNAMES);
   } else {
      atp = AsnLinkType(orig, GBTAGNAMES);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagNamesFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagNamesAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagNamesAsnWrite(GBTagNamesPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAGNAMES);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr, ASNCODE_PTRVAL_SLOT, aip, atp, GBTAGNAMES_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    GBTagsetRuleSetFree()
*
**************************************************/
NLM_EXTERN 
GBTagsetRuleSetPtr LIBCALL
GBTagsetRuleSetFree(GBTagsetRuleSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) GBTagsetRulesFree);
   return NULL;
}


/**************************************************
*
*    GBTagsetRuleSetAsnRead()
*
**************************************************/
NLM_EXTERN 
GBTagsetRuleSetPtr LIBCALL
GBTagsetRuleSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   GBTagsetRuleSetPtr ptr;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GBTagsetRuleSet ::= (self contained) */
      atp = AsnReadId(aip, amp, GBTAGSETRULESET);
   } else {
      atp = AsnLinkType(orig, GBTAGSETRULESET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) GBTagsetRulesAsnRead, (AsnOptFreeFunc) GBTagsetRulesFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = GBTagsetRuleSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    GBTagsetRuleSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GBTagsetRuleSetAsnWrite(GBTagsetRuleSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objgbseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GBTAGSETRULESET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) GBTagsetRulesAsnWrite, aip, atp, GBTAGSETRULESET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

