#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objtable.h>

static Boolean loaded = FALSE;

#include <asntable.h>
#include <objloc.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objtableAsnLoad(void)
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
*    Generated object loaders for Module NCBI-SeqTable
*    Generated using ASNCODE Revision: 6.16 at Dec 12, 2007  3:30 PM
*
**************************************************/


/**************************************************
*
*    SeqTableColumnInfoNew()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnInfoPtr LIBCALL
SeqTableColumnInfoNew(void)
{
   SeqTableColumnInfoPtr ptr = MemNew((size_t) sizeof(SeqTableColumnInfo));

   return ptr;

}


/**************************************************
*
*    SeqTableColumnInfoFree()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnInfoPtr LIBCALL
SeqTableColumnInfoFree(SeqTableColumnInfoPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> title);
   MemFree(ptr -> field_name);
   return MemFree(ptr);
}


/**************************************************
*
*    SeqTableColumnInfoAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnInfoPtr LIBCALL
SeqTableColumnInfoAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SeqTableColumnInfoPtr ptr;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTableColumnInfo ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQTABLE_COLUMN_INFO);
   } else {
      atp = AsnLinkType(orig, SEQTABLE_COLUMN_INFO);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SeqTableColumnInfoNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQTABLE_COLUMN_INFO_title) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> title = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_INFO_field_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_id = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_INFO_field_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> field_name = av.ptrvalue;
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
   ptr = SeqTableColumnInfoFree(ptr);
   goto ret;
}



/**************************************************
*
*    SeqTableColumnInfoAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableColumnInfoAsnWrite(SeqTableColumnInfoPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQTABLE_COLUMN_INFO);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> title != NULL) {
      av.ptrvalue = ptr -> title;
      retval = AsnWrite(aip, SEQTABLE_COLUMN_INFO_title,  &av);
   }
   if (ptr -> field_id || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> field_id;
      retval = AsnWrite(aip, SEQTABLE_COLUMN_INFO_field_id,  &av);
   }
   if (ptr -> field_name != NULL) {
      av.ptrvalue = ptr -> field_name;
      retval = AsnWrite(aip, SEQTABLE_COLUMN_INFO_field_name,  &av);
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
*    SeqTableColumnNew()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnPtr LIBCALL
SeqTableColumnNew(void)
{
   SeqTableColumnPtr ptr = MemNew((size_t) sizeof(SeqTableColumn));

   return ptr;

}


/**************************************************
*
*    SeqTableColumnFree()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnPtr LIBCALL
SeqTableColumnFree(SeqTableColumnPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SeqTableColumnInfoFree(ptr -> header);
   SeqTableMultiDataFree(ptr -> data);
   SeqTableSparseIndexFree(ptr -> sparse);
   SeqTableSingleDataFree(ptr -> default__);
   SeqTableSingleDataFree(ptr -> sparse_other);
   return MemFree(ptr);
}


/**************************************************
*
*    SeqTableColumnAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTableColumnPtr LIBCALL
SeqTableColumnAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SeqTableColumnPtr ptr;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTableColumn ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQTABLE_COLUMN);
   } else {
      atp = AsnLinkType(orig, SEQTABLE_COLUMN);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SeqTableColumnNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQTABLE_COLUMN_header) {
      ptr -> header = SeqTableColumnInfoAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_data) {
      ptr -> data = SeqTableMultiDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_sparse) {
      ptr -> sparse = SeqTableSparseIndexAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_default) {
      ptr -> default__ = SeqTableSingleDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQTABLE_COLUMN_sparse_other) {
      ptr -> sparse_other = SeqTableSingleDataAsnRead(aip, atp);
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
   ptr = SeqTableColumnFree(ptr);
   goto ret;
}



/**************************************************
*
*    SeqTableColumnAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableColumnAsnWrite(SeqTableColumnPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQTABLE_COLUMN);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> header != NULL) {
      if ( ! SeqTableColumnInfoAsnWrite(ptr -> header, aip, SEQTABLE_COLUMN_header)) {
         goto erret;
      }
   }
   if (ptr -> data != NULL) {
      if ( ! SeqTableMultiDataAsnWrite(ptr -> data, aip, SEQTABLE_COLUMN_data)) {
         goto erret;
      }
   }
   if (ptr -> sparse != NULL) {
      if ( ! SeqTableSparseIndexAsnWrite(ptr -> sparse, aip, SEQTABLE_COLUMN_sparse)) {
         goto erret;
      }
   }
   if (ptr -> default__ != NULL) {
      if ( ! SeqTableSingleDataAsnWrite(ptr -> default__, aip, SEQTABLE_COLUMN_default)) {
         goto erret;
      }
   }
   if (ptr -> sparse_other != NULL) {
      if ( ! SeqTableSingleDataAsnWrite(ptr -> sparse_other, aip, SEQTABLE_COLUMN_sparse_other)) {
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
*    SeqTableNew()
*
**************************************************/
NLM_EXTERN 
SeqTablePtr LIBCALL
SeqTableNew(void)
{
   SeqTablePtr ptr = MemNew((size_t) sizeof(SeqTable));

   return ptr;

}


/**************************************************
*
*    SeqTableFree()
*
**************************************************/
NLM_EXTERN 
SeqTablePtr LIBCALL
SeqTableFree(SeqTablePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> columns, (AsnOptFreeFunc) SeqTableColumnFree);
   return MemFree(ptr);
}


/**************************************************
*
*    SeqTableAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTablePtr LIBCALL
SeqTableAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SeqTablePtr ptr;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTable ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQ_TABLE);
   } else {
      atp = AsnLinkType(orig, SEQ_TABLE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SeqTableNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQ_TABLE_feat_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feat_type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQ_TABLE_feat_subtype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> feat_subtype = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQ_TABLE_num_rows) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> num_rows = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQ_TABLE_columns) {
      ptr -> columns = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqTableColumnAsnRead, (AsnOptFreeFunc) SeqTableColumnFree);
      if (isError && ptr -> columns == NULL) {
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
   ptr = SeqTableFree(ptr);
   goto ret;
}



/**************************************************
*
*    SeqTableAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableAsnWrite(SeqTablePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQ_TABLE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> feat_type;
   retval = AsnWrite(aip, SEQ_TABLE_feat_type,  &av);
   if (ptr -> feat_subtype || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> feat_subtype;
      retval = AsnWrite(aip, SEQ_TABLE_feat_subtype,  &av);
   }
   av.intvalue = ptr -> num_rows;
   retval = AsnWrite(aip, SEQ_TABLE_num_rows,  &av);
   AsnGenericUserSeqOfAsnWrite(ptr -> columns, (AsnWriteFunc) SeqTableColumnAsnWrite, aip, SEQ_TABLE_columns, SEQ_TABLE_columns_E);
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
*    CommonStringTableNew()
*
**************************************************/
NLM_EXTERN 
CommonStringTablePtr LIBCALL
CommonStringTableNew(void)
{
   CommonStringTablePtr ptr = MemNew((size_t) sizeof(CommonStringTable));

   return ptr;

}


/**************************************************
*
*    CommonStringTableFree()
*
**************************************************/
NLM_EXTERN 
CommonStringTablePtr LIBCALL
CommonStringTableFree(CommonStringTablePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> strings ,ASNCODE_PTRVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> indexes ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    CommonStringTableAsnRead()
*
**************************************************/
NLM_EXTERN 
CommonStringTablePtr LIBCALL
CommonStringTableAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CommonStringTablePtr ptr;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CommonStringTable ::= (self contained) */
      atp = AsnReadId(aip, amp, COMMONSTRING_TABLE);
   } else {
      atp = AsnLinkType(orig, COMMONSTRING_TABLE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CommonStringTableNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == COMMONSTRING_TABLE_strings) {
      ptr -> strings = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> strings == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMONSTRING_TABLE_indexes) {
      ptr -> indexes = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> indexes == NULL) {
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
   ptr = CommonStringTableFree(ptr);
   goto ret;
}



/**************************************************
*
*    CommonStringTableAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CommonStringTableAsnWrite(CommonStringTablePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, COMMONSTRING_TABLE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> strings ,ASNCODE_PTRVAL_SLOT, aip, COMMONSTRING_TABLE_strings, COMMONSTRING_TABLE_strings_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> indexes ,ASNCODE_INTVAL_SLOT, aip, COMMONSTRING_TABLE_indexes, COMMONSTRING_TABLE_indexes_E);
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
*    CommonBytesTableNew()
*
**************************************************/
NLM_EXTERN 
CommonBytesTablePtr LIBCALL
CommonBytesTableNew(void)
{
   CommonBytesTablePtr ptr = MemNew((size_t) sizeof(CommonBytesTable));

   return ptr;

}


/**************************************************
*
*    CommonBytesTableFree()
*
**************************************************/
NLM_EXTERN 
CommonBytesTablePtr LIBCALL
CommonBytesTableFree(CommonBytesTablePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> bytes ,ASNCODE_BYTEVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> indexes ,ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    CommonBytesTableAsnRead()
*
**************************************************/
NLM_EXTERN 
CommonBytesTablePtr LIBCALL
CommonBytesTableAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CommonBytesTablePtr ptr;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CommonBytesTable ::= (self contained) */
      atp = AsnReadId(aip, amp, COMMONBYTES_TABLE);
   } else {
      atp = AsnLinkType(orig, COMMONBYTES_TABLE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CommonBytesTableNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == COMMONBYTES_TABLE_bytes) {
      ptr -> bytes = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_BYTEVAL_SLOT, &isError);
      if (isError && ptr -> bytes == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == COMMONBYTES_TABLE_indexes) {
      ptr -> indexes = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> indexes == NULL) {
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
   ptr = CommonBytesTableFree(ptr);
   goto ret;
}



/**************************************************
*
*    CommonBytesTableAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CommonBytesTableAsnWrite(CommonBytesTablePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, COMMONBYTES_TABLE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> bytes ,ASNCODE_BYTEVAL_SLOT, aip, COMMONBYTES_TABLE_bytes, COMMONBYTES_TABLE_bytes_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> indexes ,ASNCODE_INTVAL_SLOT, aip, COMMONBYTES_TABLE_indexes, COMMONBYTES_TABLE_indexes_E);
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
*    SeqTableMultiDataFree()
*
**************************************************/
NLM_EXTERN 
SeqTableMultiDataPtr LIBCALL
SeqTableMultiDataFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SeqTableMultiData_int__:
      AsnGenericBaseSeqOfFree((ValNodePtr) pnt,ASNCODE_INTVAL_SLOT);
      break;
   case SeqTableMultiData_real:
      AsnGenericBaseSeqOfFree((ValNodePtr) pnt,ASNCODE_REALVAL_SLOT);
      break;
   case SeqTableMultiData_string:
      AsnGenericBaseSeqOfFree((ValNodePtr) pnt,ASNCODE_PTRVAL_SLOT);
      break;
   case SeqTableMultiData_bytes:
      AsnGenericBaseSeqOfFree((ValNodePtr) pnt,ASNCODE_BYTEVAL_SLOT);
      break;
   case SeqTableMultiData_common_string:
      CommonStringTableFree(anp -> data.ptrvalue);
      break;
   case SeqTableMultiData_common_bytes:
      CommonBytesTableFree(anp -> data.ptrvalue);
      break;
   case SeqTableMultiData_bit:
      BSFree(anp -> data.ptrvalue);
      break;
   case SeqTableMultiData_loc:
      AsnGenericChoiceSeqOfFree((Pointer) pnt, (AsnOptFreeFunc) SeqLocFree);
      break;
   case SeqTableMultiData_id:
      AsnGenericChoiceSeqOfFree((Pointer) pnt, (AsnOptFreeFunc) SeqIdFree);
      break;
   case SeqTableMultiData_interval:
      AsnGenericUserSeqOfFree((Pointer) pnt, (AsnOptFreeFunc) SeqIntFree);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SeqTableMultiDataAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTableMultiDataPtr LIBCALL
SeqTableMultiDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTableMultiData ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQTABLE_MULTI_DATA);
   } else {
      atp = AsnLinkType(orig, SEQTABLE_MULTI_DATA);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEQTABLE_MULTI_DATA_int) {
      choice = SeqTableMultiData_int__;
      anp -> data.ptrvalue = 
      AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_MULTI_DATA_real) {
      choice = SeqTableMultiData_real;
      anp -> data.ptrvalue = 
      AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_REALVAL_SLOT, &isError);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_MULTI_DATA_string) {
      choice = SeqTableMultiData_string;
      anp -> data.ptrvalue = 
      AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_MULTI_DATA_bytes) {
      choice = SeqTableMultiData_bytes;
      anp -> data.ptrvalue = 
      AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_BYTEVAL_SLOT, &isError);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == MULTI_DATA_common_string) {
      choice = SeqTableMultiData_common_string;
      func = (AsnReadFunc) CommonStringTableAsnRead;
   }
   else if (atp == MULTI_DATA_common_bytes) {
      choice = SeqTableMultiData_common_bytes;
      func = (AsnReadFunc) CommonBytesTableAsnRead;
   }
   else if (atp == SEQTABLE_MULTI_DATA_bit) {
      choice = SeqTableMultiData_bit;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == SEQTABLE_MULTI_DATA_loc) {
      choice = SeqTableMultiData_loc;
      anp -> data.ptrvalue =
      AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqLocAsnRead,             (AsnOptFreeFunc) SeqLocFree);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_MULTI_DATA_id) {
      choice = SeqTableMultiData_id;
      anp -> data.ptrvalue =
      AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqIdAsnRead,             (AsnOptFreeFunc) SeqIdFree);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_MULTI_DATA_interval) {
      choice = SeqTableMultiData_interval;
      anp -> data.ptrvalue =
      AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqIntAsnRead,             (AsnOptFreeFunc) SeqIntFree);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SeqTableMultiDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableMultiDataAsnWrite(SeqTableMultiDataPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEQTABLE_MULTI_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SeqTableMultiData_int__:
      retval = AsnGenericBaseSeqOfAsnWrite((Pointer) pnt,ASNCODE_INTVAL_SLOT, aip, SEQTABLE_MULTI_DATA_int, SEQTABLE_MULTI_DATA_int_E);            break;
   case SeqTableMultiData_real:
      retval = AsnGenericBaseSeqOfAsnWrite((Pointer) pnt,ASNCODE_REALVAL_SLOT, aip, SEQTABLE_MULTI_DATA_real, SEQTABLE_MULTI_DATA_real_E);            break;
   case SeqTableMultiData_string:
      retval = AsnGenericBaseSeqOfAsnWrite((Pointer) pnt,ASNCODE_PTRVAL_SLOT, aip, SEQTABLE_MULTI_DATA_string, SEQTABLE_MULTI_DATA_string_E);            break;
   case SeqTableMultiData_bytes:
      retval = AsnGenericBaseSeqOfAsnWrite((Pointer) pnt,ASNCODE_BYTEVAL_SLOT, aip, SEQTABLE_MULTI_DATA_bytes, SEQTABLE_MULTI_DATA_bytes_E);            break;
   case SeqTableMultiData_common_string:
      writetype = MULTI_DATA_common_string;
      func = (AsnWriteFunc) CommonStringTableAsnWrite;
      break;
   case SeqTableMultiData_common_bytes:
      writetype = MULTI_DATA_common_bytes;
      func = (AsnWriteFunc) CommonBytesTableAsnWrite;
      break;
   case SeqTableMultiData_bit:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEQTABLE_MULTI_DATA_bit, &av);
      break;
   case SeqTableMultiData_loc:
      retval = AsnGenericChoiceSeqOfAsnWrite((Pointer) pnt, (AsnWriteFunc) SeqLocAsnWrite, aip, SEQTABLE_MULTI_DATA_loc, SEQTABLE_MULTI_DATA_loc_E);
      break;
   case SeqTableMultiData_id:
      retval = AsnGenericChoiceSeqOfAsnWrite((Pointer) pnt, (AsnWriteFunc) SeqIdAsnWrite, aip, SEQTABLE_MULTI_DATA_id, SEQTABLE_MULTI_DATA_id_E);
      break;
   case SeqTableMultiData_interval:
      retval = AsnGenericUserSeqOfAsnWrite((Pointer) pnt, (AsnWriteFunc) SeqIntAsnWrite, aip, SEQTABLE_MULTI_DATA_interval, SEQTABLE_MULTI_DATA_interval_E);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SeqTableSingleDataFree()
*
**************************************************/
NLM_EXTERN 
SeqTableSingleDataPtr LIBCALL
SeqTableSingleDataFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SeqTableSingleData_string:
      MemFree(anp -> data.ptrvalue);
      break;
   case SeqTableSingleData_bytes:
      BSFree(anp -> data.ptrvalue);
      break;
   case SeqTableSingleData_loc:
      SeqLocFree(anp -> data.ptrvalue);
      break;
   case SeqTableSingleData_id:
      SeqIdFree(anp -> data.ptrvalue);
      break;
   case SeqTableSingleData_interval:
      SeqIntFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SeqTableSingleDataAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTableSingleDataPtr LIBCALL
SeqTableSingleDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTableSingleData ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQTABLE_SINGLE_DATA);
   } else {
      atp = AsnLinkType(orig, SEQTABLE_SINGLE_DATA);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEQTABLE_SINGLE_DATA_int) {
      choice = SeqTableSingleData_int__;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.intvalue = av.intvalue;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_real) {
      choice = SeqTableSingleData_real;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.realvalue = av.realvalue;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_string) {
      choice = SeqTableSingleData_string;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_bytes) {
      choice = SeqTableSingleData_bytes;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_bit) {
      choice = SeqTableSingleData_bit;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_loc) {
      choice = SeqTableSingleData_loc;
      func = (AsnReadFunc) SeqLocAsnRead;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_id) {
      choice = SeqTableSingleData_id;
      func = (AsnReadFunc) SeqIdAsnRead;
   }
   else if (atp == SEQTABLE_SINGLE_DATA_interval) {
      choice = SeqTableSingleData_interval;
      func = (AsnReadFunc) SeqIntAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SeqTableSingleDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableSingleDataAsnWrite(SeqTableSingleDataPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEQTABLE_SINGLE_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SeqTableSingleData_int__:
      av.intvalue = anp->data.intvalue;
      retval = AsnWrite(aip, SEQTABLE_SINGLE_DATA_int, &av);
      break;
   case SeqTableSingleData_real:
      av.realvalue = anp->data.realvalue;
      retval = AsnWrite(aip, SEQTABLE_SINGLE_DATA_real, &av);
      break;
   case SeqTableSingleData_string:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEQTABLE_SINGLE_DATA_string, &av);
      break;
   case SeqTableSingleData_bytes:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEQTABLE_SINGLE_DATA_bytes, &av);
      break;
   case SeqTableSingleData_bit:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, SEQTABLE_SINGLE_DATA_bit, &av);
      break;
   case SeqTableSingleData_loc:
      writetype = SEQTABLE_SINGLE_DATA_loc;
      func = (AsnWriteFunc) SeqLocAsnWrite;
      break;
   case SeqTableSingleData_id:
      writetype = SEQTABLE_SINGLE_DATA_id;
      func = (AsnWriteFunc) SeqIdAsnWrite;
      break;
   case SeqTableSingleData_interval:
      writetype = SEQTABLE_SINGLE_DATA_interval;
      func = (AsnWriteFunc) SeqIntAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SeqTableSparseIndexFree()
*
**************************************************/
NLM_EXTERN 
SeqTableSparseIndexPtr LIBCALL
SeqTableSparseIndexFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case SeqTableSparseIndex_indexes:
      AsnGenericBaseSeqOfFree((ValNodePtr) pnt,ASNCODE_INTVAL_SLOT);
      break;
   case SeqTableSparseIndex_bit_set:
      BSFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    SeqTableSparseIndexAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqTableSparseIndexPtr LIBCALL
SeqTableSparseIndexAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean isError = FALSE;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! objtableAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqTableSparseIndex ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQTABLE_SPARSE_INDEX);
   } else {
      atp = AsnLinkType(orig, SEQTABLE_SPARSE_INDEX);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == SEQTABLE_SPARSE_INDEX_indexes) {
      choice = SeqTableSparseIndex_indexes;
      anp -> data.ptrvalue = 
      AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && anp -> data.ptrvalue == NULL) {
         goto erret;
      }
   }
   else if (atp == SEQTABLE_SPARSE_INDEX_bit_set) {
      choice = SeqTableSparseIndex_bit_set;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SeqTableSparseIndexAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqTableSparseIndexAsnWrite(SeqTableSparseIndexPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objtableAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, SEQTABLE_SPARSE_INDEX);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case SeqTableSparseIndex_indexes:
      retval = AsnGenericBaseSeqOfAsnWrite((Pointer) pnt,ASNCODE_INTVAL_SLOT, aip, SEQTABLE_SPARSE_INDEX_indexes, SEQTABLE_SPARSE_INDEX_indexes_E);            break;
   case SeqTableSparseIndex_bit_set:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, SEQTABLE_SPARSE_INDEX_bit_set, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}
