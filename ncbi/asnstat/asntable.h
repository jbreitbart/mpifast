/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asntable.h15";
static AsnValxNode avnx[32] = {
    {20,"location" ,0,0.0,&avnx[1] } ,
    {20,"location-id" ,1,0.0,&avnx[2] } ,
    {20,"location-gi" ,2,0.0,&avnx[3] } ,
    {20,"location-from" ,3,0.0,&avnx[4] } ,
    {20,"location-to" ,4,0.0,&avnx[5] } ,
    {20,"location-strand" ,5,0.0,&avnx[6] } ,
    {20,"location-fuzz-from-lim" ,6,0.0,&avnx[7] } ,
    {20,"location-fuzz-to-lim" ,7,0.0,&avnx[8] } ,
    {20,"product" ,10,0.0,&avnx[9] } ,
    {20,"product-id" ,11,0.0,&avnx[10] } ,
    {20,"product-gi" ,12,0.0,&avnx[11] } ,
    {20,"product-from" ,13,0.0,&avnx[12] } ,
    {20,"product-to" ,14,0.0,&avnx[13] } ,
    {20,"product-strand" ,15,0.0,&avnx[14] } ,
    {20,"product-fuzz-from-lim" ,16,0.0,&avnx[15] } ,
    {20,"product-fuzz-to-lim" ,17,0.0,&avnx[16] } ,
    {20,"id-local" ,20,0.0,&avnx[17] } ,
    {20,"xref-id-local" ,21,0.0,&avnx[18] } ,
    {20,"partial" ,22,0.0,&avnx[19] } ,
    {20,"comment" ,23,0.0,&avnx[20] } ,
    {20,"title" ,24,0.0,&avnx[21] } ,
    {20,"ext" ,25,0.0,&avnx[22] } ,
    {20,"qual" ,26,0.0,&avnx[23] } ,
    {20,"dbxref" ,27,0.0,&avnx[24] } ,
    {20,"data-imp-key" ,30,0.0,&avnx[25] } ,
    {20,"data-region" ,31,0.0,&avnx[26] } ,
    {20,"data-cdregion-frame" ,32,0.0,&avnx[27] } ,
    {20,"ext-type" ,40,0.0,&avnx[28] } ,
    {20,"qual-qual" ,41,0.0,&avnx[29] } ,
    {20,"qual-val" ,42,0.0,&avnx[30] } ,
    {20,"dbxref-db" ,43,0.0,&avnx[31] } ,
    {20,"dbxref-tag" ,44,0.0,NULL } };

static AsnType atx[68] = {
  {401, "SeqTable-column-info" ,1,0,0,0,0,1,0,0,NULL,&atx[6],&atx[1],0,&atx[7]} ,
  {0, "title" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "field-id" ,128,1,0,1,0,0,0,0,NULL,&atx[4],&avnx[0],0,&atx[5]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "field-name" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "SeqTable-column" ,1,0,0,0,0,1,0,0,NULL,&atx[6],&atx[8],0,&atx[62]} ,
  {0, "header" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[9]} ,
  {0, "data" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[45]} ,
  {409, "SeqTable-multi-data" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[11],0,&atx[51]} ,
  {0, "int" ,128,0,0,0,0,0,0,0,NULL,&atx[13],&atx[12],0,&atx[14]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "real" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[15],0,&atx[17]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "string" ,128,2,0,0,0,0,0,0,NULL,&atx[13],&atx[18],0,&atx[19]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "bytes" ,128,3,0,0,0,0,0,0,NULL,&atx[13],&atx[20],0,&atx[22]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "common-string" ,128,4,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[28]} ,
  {407, "CommonString-table" ,1,0,0,0,0,0,0,0,NULL,&atx[6],&atx[24],0,&atx[29]} ,
  {0, "strings" ,128,0,0,0,0,0,0,0,NULL,&atx[13],&atx[25],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "indexes" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[27],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "common-bytes" ,128,5,0,0,0,0,0,0,NULL,&atx[29],NULL,0,&atx[34]} ,
  {408, "CommonBytes-table" ,1,0,0,0,0,0,0,0,NULL,&atx[6],&atx[30],0,&atx[10]} ,
  {0, "bytes" ,128,0,0,0,0,0,0,0,NULL,&atx[13],&atx[31],0,&atx[32]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {0, "indexes" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[33],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "bit" ,128,6,0,0,0,0,0,0,NULL,&atx[21],NULL,0,&atx[35]} ,
  {0, "loc" ,128,7,0,0,0,0,0,0,NULL,&atx[13],&atx[36],0,&atx[38]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[37],NULL,0,NULL} ,
  {405, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[43]} ,
  {0, "id" ,128,8,0,0,0,0,0,0,NULL,&atx[13],&atx[39],0,&atx[41]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[40],NULL,0,NULL} ,
  {404, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[37]} ,
  {0, "interval" ,128,9,0,0,0,0,0,0,NULL,&atx[13],&atx[42],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[43],NULL,0,NULL} ,
  {406, "Seq-interval" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[23]} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "sparse" ,128,2,0,1,0,0,0,0,NULL,&atx[46],NULL,0,&atx[50]} ,
  {411, "SeqTable-sparse-index" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[47],0,NULL} ,
  {0, "indexes" ,128,0,0,0,0,0,0,0,NULL,&atx[13],&atx[48],0,&atx[49]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "bit-set" ,128,1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {0, "default" ,128,3,0,1,0,0,0,0,NULL,&atx[51],NULL,0,&atx[61]} ,
  {410, "SeqTable-single-data" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[52],0,&atx[46]} ,
  {0, "int" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[53]} ,
  {0, "real" ,128,1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[54]} ,
  {0, "string" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[55]} ,
  {0, "bytes" ,128,3,0,0,0,0,0,0,NULL,&atx[21],NULL,0,&atx[56]} ,
  {0, "bit" ,128,4,0,0,0,0,0,0,NULL,&atx[57],NULL,0,&atx[58]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "loc" ,128,5,0,0,0,0,0,0,NULL,&atx[37],NULL,0,&atx[59]} ,
  {0, "id" ,128,6,0,0,0,0,0,0,NULL,&atx[40],NULL,0,&atx[60]} ,
  {0, "interval" ,128,7,0,0,0,0,0,0,NULL,&atx[43],NULL,0,NULL} ,
  {0, "sparse-other" ,128,4,0,1,0,0,0,0,NULL,&atx[51],NULL,0,NULL} ,
  {403, "Seq-table" ,1,0,0,0,0,1,0,0,NULL,&atx[6],&atx[63],0,&atx[40]} ,
  {0, "feat-type" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[64]} ,
  {0, "feat-subtype" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[65]} ,
  {0, "num-rows" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[66]} ,
  {0, "columns" ,128,3,0,0,0,0,0,0,NULL,&atx[13],&atx[67],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-SeqTable" , "asntable.h15",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-SeqTable
*
**************************************************/

#define SEQTABLE_COLUMN_INFO &at[0]
#define SEQTABLE_COLUMN_INFO_title &at[1]
#define SEQTABLE_COLUMN_INFO_field_id &at[3]
#define SEQTABLE_COLUMN_INFO_field_name &at[5]

#define SEQTABLE_COLUMN &at[7]
#define SEQTABLE_COLUMN_header &at[8]
#define SEQTABLE_COLUMN_data &at[9]
#define SEQTABLE_COLUMN_sparse &at[45]
#define SEQTABLE_COLUMN_default &at[50]
#define SEQTABLE_COLUMN_sparse_other &at[61]

#define SEQ_TABLE &at[62]
#define SEQ_TABLE_feat_type &at[63]
#define SEQ_TABLE_feat_subtype &at[64]
#define SEQ_TABLE_num_rows &at[65]
#define SEQ_TABLE_columns &at[66]
#define SEQ_TABLE_columns_E &at[67]

#define COMMONSTRING_TABLE &at[23]
#define COMMONSTRING_TABLE_strings &at[24]
#define COMMONSTRING_TABLE_strings_E &at[25]
#define COMMONSTRING_TABLE_indexes &at[26]
#define COMMONSTRING_TABLE_indexes_E &at[27]

#define COMMONBYTES_TABLE &at[29]
#define COMMONBYTES_TABLE_bytes &at[30]
#define COMMONBYTES_TABLE_bytes_E &at[31]
#define COMMONBYTES_TABLE_indexes &at[32]
#define COMMONBYTES_TABLE_indexes_E &at[33]

#define SEQTABLE_MULTI_DATA &at[10]
#define SEQTABLE_MULTI_DATA_int &at[11]
#define SEQTABLE_MULTI_DATA_int_E &at[12]
#define SEQTABLE_MULTI_DATA_real &at[14]
#define SEQTABLE_MULTI_DATA_real_E &at[15]
#define SEQTABLE_MULTI_DATA_string &at[17]
#define SEQTABLE_MULTI_DATA_string_E &at[18]
#define SEQTABLE_MULTI_DATA_bytes &at[19]
#define SEQTABLE_MULTI_DATA_bytes_E &at[20]
#define MULTI_DATA_common_string &at[22]
#define MULTI_DATA_common_bytes &at[28]
#define SEQTABLE_MULTI_DATA_bit &at[34]
#define SEQTABLE_MULTI_DATA_loc &at[35]
#define SEQTABLE_MULTI_DATA_loc_E &at[36]
#define SEQTABLE_MULTI_DATA_id &at[38]
#define SEQTABLE_MULTI_DATA_id_E &at[39]
#define SEQTABLE_MULTI_DATA_interval &at[41]
#define SEQTABLE_MULTI_DATA_interval_E &at[42]

#define SEQTABLE_SINGLE_DATA &at[51]
#define SEQTABLE_SINGLE_DATA_int &at[52]
#define SEQTABLE_SINGLE_DATA_real &at[53]
#define SEQTABLE_SINGLE_DATA_string &at[54]
#define SEQTABLE_SINGLE_DATA_bytes &at[55]
#define SEQTABLE_SINGLE_DATA_bit &at[56]
#define SEQTABLE_SINGLE_DATA_loc &at[58]
#define SEQTABLE_SINGLE_DATA_id &at[59]
#define SEQTABLE_SINGLE_DATA_interval &at[60]

#define SEQTABLE_SPARSE_INDEX &at[46]
#define SEQTABLE_SPARSE_INDEX_indexes &at[47]
#define SEQTABLE_SPARSE_INDEX_indexes_E &at[48]
#define SEQTABLE_SPARSE_INDEX_bit_set &at[49]
