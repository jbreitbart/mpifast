#ifndef _objtable_ 
#define _objtable_ 

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
*    Generated objects for Module NCBI-SeqTable
*    Generated using ASNCODE Revision: 6.16 at Dec 12, 2007  3:30 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objtableAsnLoad PROTO((void));


/**************************************************
*
*    SeqTableColumnInfo
*
**************************************************/
typedef struct struct_SeqTable_column_info {
   Uint4 OBbits__;
   CharPtr   title;
#define OB__SeqTable_column_info_field_id 0

   Int4   field_id;
   CharPtr   field_name;
} SeqTableColumnInfo, PNTR SeqTableColumnInfoPtr;


NLM_EXTERN SeqTableColumnInfoPtr LIBCALL SeqTableColumnInfoFree PROTO ((SeqTableColumnInfoPtr ));
NLM_EXTERN SeqTableColumnInfoPtr LIBCALL SeqTableColumnInfoNew PROTO (( void ));
NLM_EXTERN SeqTableColumnInfoPtr LIBCALL SeqTableColumnInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableColumnInfoAsnWrite PROTO (( SeqTableColumnInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SeqTableColumn
*
**************************************************/
typedef struct struct_SeqTable_column {
   struct struct_SeqTable_column PNTR next;
   Uint4 OBbits__;
   struct struct_SeqTable_column_info PNTR   header;
   ValNodePtr   data;
   ValNodePtr   sparse;
   ValNodePtr   default__;
   ValNodePtr   sparse_other;
} SeqTableColumn, PNTR SeqTableColumnPtr;


NLM_EXTERN SeqTableColumnPtr LIBCALL SeqTableColumnFree PROTO ((SeqTableColumnPtr ));
NLM_EXTERN SeqTableColumnPtr LIBCALL SeqTableColumnNew PROTO (( void ));
NLM_EXTERN SeqTableColumnPtr LIBCALL SeqTableColumnAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableColumnAsnWrite PROTO (( SeqTableColumnPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SeqTable
*
**************************************************/
typedef struct struct_Seq_table {
   Uint4 OBbits__;
   Int4   feat_type;
#define OB__Seq_table_feat_subtype 0

   Int4   feat_subtype;
   Int4   num_rows;
   struct struct_SeqTable_column PNTR   columns;
} SeqTable, PNTR SeqTablePtr;


NLM_EXTERN SeqTablePtr LIBCALL SeqTableFree PROTO ((SeqTablePtr ));
NLM_EXTERN SeqTablePtr LIBCALL SeqTableNew PROTO (( void ));
NLM_EXTERN SeqTablePtr LIBCALL SeqTableAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableAsnWrite PROTO (( SeqTablePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CommonStringTable
*
**************************************************/
typedef struct struct_CommonString_table {
   Uint4 OBbits__;
   ValNodePtr   strings;
   ValNodePtr   indexes;
} CommonStringTable, PNTR CommonStringTablePtr;


NLM_EXTERN CommonStringTablePtr LIBCALL CommonStringTableFree PROTO ((CommonStringTablePtr ));
NLM_EXTERN CommonStringTablePtr LIBCALL CommonStringTableNew PROTO (( void ));
NLM_EXTERN CommonStringTablePtr LIBCALL CommonStringTableAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommonStringTableAsnWrite PROTO (( CommonStringTablePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CommonBytesTable
*
**************************************************/
typedef struct struct_CommonBytes_table {
   Uint4 OBbits__;
   ValNodePtr   bytes;
   ValNodePtr   indexes;
} CommonBytesTable, PNTR CommonBytesTablePtr;


NLM_EXTERN CommonBytesTablePtr LIBCALL CommonBytesTableFree PROTO ((CommonBytesTablePtr ));
NLM_EXTERN CommonBytesTablePtr LIBCALL CommonBytesTableNew PROTO (( void ));
NLM_EXTERN CommonBytesTablePtr LIBCALL CommonBytesTableAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommonBytesTableAsnWrite PROTO (( CommonBytesTablePtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr SeqTableMultiDataPtr;
typedef ValNode SeqTableMultiData;
#define SeqTableMultiData_int__ 1
#define SeqTableMultiData_real 2
#define SeqTableMultiData_string 3
#define SeqTableMultiData_bytes 4
#define SeqTableMultiData_common_string 5
#define SeqTableMultiData_common_bytes 6
#define SeqTableMultiData_bit 7
#define SeqTableMultiData_loc 8
#define SeqTableMultiData_id 9
#define SeqTableMultiData_interval 10


NLM_EXTERN SeqTableMultiDataPtr LIBCALL SeqTableMultiDataFree PROTO ((SeqTableMultiDataPtr ));
NLM_EXTERN SeqTableMultiDataPtr LIBCALL SeqTableMultiDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableMultiDataAsnWrite PROTO (( SeqTableMultiDataPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr SeqTableSingleDataPtr;
typedef ValNode SeqTableSingleData;
#define SeqTableSingleData_int__ 1
#define SeqTableSingleData_real 2
#define SeqTableSingleData_string 3
#define SeqTableSingleData_bytes 4
#define SeqTableSingleData_bit 5
#define SeqTableSingleData_loc 6
#define SeqTableSingleData_id 7
#define SeqTableSingleData_interval 8


NLM_EXTERN SeqTableSingleDataPtr LIBCALL SeqTableSingleDataFree PROTO ((SeqTableSingleDataPtr ));
NLM_EXTERN SeqTableSingleDataPtr LIBCALL SeqTableSingleDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableSingleDataAsnWrite PROTO (( SeqTableSingleDataPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr SeqTableSparseIndexPtr;
typedef ValNode SeqTableSparseIndex;
#define SeqTableSparseIndex_indexes 1
#define SeqTableSparseIndex_bit_set 2


NLM_EXTERN SeqTableSparseIndexPtr LIBCALL SeqTableSparseIndexFree PROTO ((SeqTableSparseIndexPtr ));
NLM_EXTERN SeqTableSparseIndexPtr LIBCALL SeqTableSparseIndexAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqTableSparseIndexAsnWrite PROTO (( SeqTableSparseIndexPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objtable_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

