/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "seqsplit.h19";
static AsnValxNode avnx[4] = {
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } };

static AsnType atx[139] = {
  {401, "ID2S-Chunk-Id" ,1,0,0,0,0,1,0,0,NULL,&atx[1],NULL,0,&atx[2]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "ID2S-Seq-annot-Info" ,1,0,0,0,0,1,0,0,NULL,&atx[15],&atx[3],0,&atx[20]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "align" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "graph" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {0, "feat" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[9],0,&atx[16]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {425, "ID2S-Feat-type-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[11],0,&atx[17]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[12]} ,
  {0, "subtypes" ,128,1,0,1,0,0,0,0,NULL,&atx[14],&atx[13],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seq-loc" ,128,4,0,1,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {426, "ID2S-Seq-loc" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[18],0,&atx[115]} ,
  {0, "whole-gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[19]} ,
  {0, "whole-seq-id" ,128,1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[21]} ,
  {403, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[51]} ,
  {0, "whole-gi-range" ,128,2,0,0,0,0,0,0,NULL,&atx[22],NULL,0,&atx[25]} ,
  {430, "ID2S-Gi-Range" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[23],0,&atx[26]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[24]} ,
  {0, "count" ,128,1,0,0,1,0,0,0,&avnx[0],&atx[1],NULL,0,NULL} ,
  {0, "gi-interval" ,128,3,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[30]} ,
  {431, "ID2S-Gi-Interval" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[27],0,&atx[31]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[28]} ,
  {0, "start" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[29]} ,
  {0, "length" ,128,2,0,0,1,0,0,0,&avnx[1],&atx[1],NULL,0,NULL} ,
  {0, "seq-id-interval" ,128,4,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[35]} ,
  {432, "ID2S-Seq-id-Interval" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[32],0,&atx[36]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[33]} ,
  {0, "start" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[34]} ,
  {0, "length" ,128,2,0,0,1,0,0,0,&avnx[2],&atx[1],NULL,0,NULL} ,
  {0, "gi-ints" ,128,5,0,0,0,0,0,0,NULL,&atx[36],NULL,0,&atx[43]} ,
  {433, "ID2S-Gi-Ints" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[37],0,&atx[44]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[38]} ,
  {0, "ints" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[39],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[40],NULL,0,NULL} ,
  {435, "ID2S-Interval" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[41],0,NULL} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[42]} ,
  {0, "length" ,128,1,0,0,1,0,0,0,&avnx[3],&atx[1],NULL,0,NULL} ,
  {0, "seq-id-ints" ,128,6,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[48]} ,
  {434, "ID2S-Seq-id-Ints" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[45],0,&atx[40]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[46]} ,
  {0, "ints" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[47],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[40],NULL,0,NULL} ,
  {0, "loc-set" ,128,7,0,0,0,0,0,0,NULL,&atx[14],&atx[49],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {404, "Seq-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[52]} ,
  {405, "Bioseq" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[53]} ,
  {406, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[54]} ,
  {407, "Seq-descr" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[55]} ,
  {408, "Seq-literal" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[56]} ,
  {409, "Seq-align" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[57]} ,
  {410, "Feat-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[58]} ,
  {411, "ID2S-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[59],0,&atx[61]} ,
  {0, "bioseqs-info" ,128,0,0,1,0,0,0,0,NULL,&atx[14],&atx[60],0,&atx[73]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {412, "ID2S-Bioseqs-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[62],0,&atx[75]} ,
  {0, "info" ,128,0,0,0,0,0,0,0,NULL,&atx[63],NULL,0,&atx[67]} ,
  {414, "ID2S-Bioseq-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[64],0,&atx[68]} ,
  {0, "gap-count" ,128,0,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[65]} ,
  {0, "seq-map-has-ref" ,128,1,0,1,0,0,0,0,NULL,&atx[66],NULL,0,NULL} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bioseqs" ,128,1,0,0,0,0,0,0,NULL,&atx[68],NULL,0,NULL} ,
  {415, "ID2S-Bioseq-Ids" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[69],0,&atx[79]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],&atx[70],0,NULL} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[71]} ,
  {0, "seq-id" ,128,1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[72]} ,
  {0, "gi-range" ,128,2,0,0,0,0,0,0,NULL,&atx[22],NULL,0,NULL} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[74],0,&atx[114]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[75],NULL,0,NULL} ,
  {413, "ID2S-Chunk-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[76],0,&atx[63]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[77]} ,
  {0, "content" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[78],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[79],NULL,0,NULL} ,
  {416, "ID2S-Chunk-Content" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[80],0,&atx[81]} ,
  {0, "seq-descr" ,128,0,0,0,0,0,0,0,NULL,&atx[81],NULL,0,&atx[87]} ,
  {417, "ID2S-Seq-descr-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[82],0,&atx[89]} ,
  {0, "type-mask" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[83]} ,
  {0, "bioseqs" ,128,1,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[84]} ,
  {0, "bioseq-sets" ,128,2,0,1,0,0,0,0,NULL,&atx[85],NULL,0,NULL} ,
  {424, "ID2S-Bioseq-set-Ids" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[86],0,&atx[10]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "seq-annot" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[88]} ,
  {0, "seq-assembly" ,128,2,0,0,0,0,0,0,NULL,&atx[89],NULL,0,&atx[91]} ,
  {418, "ID2S-Seq-assembly-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[90],0,&atx[92]} ,
  {0, "bioseqs" ,128,0,0,0,0,0,0,0,NULL,&atx[68],NULL,0,NULL} ,
  {0, "seq-map" ,128,3,0,0,0,0,0,0,NULL,&atx[92],NULL,0,&atx[93]} ,
  {419, "ID2S-Seq-map-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[94]} ,
  {0, "seq-data" ,128,4,0,0,0,0,0,0,NULL,&atx[94],NULL,0,&atx[95]} ,
  {420, "ID2S-Seq-data-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[96]} ,
  {0, "seq-annot-place" ,128,5,0,0,0,0,0,0,NULL,&atx[96],NULL,0,&atx[100]} ,
  {421, "ID2S-Seq-annot-place-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[97],0,&atx[102]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[98]} ,
  {0, "bioseqs" ,128,1,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[99]} ,
  {0, "bioseq-sets" ,128,2,0,1,0,0,0,0,NULL,&atx[85],NULL,0,NULL} ,
  {0, "bioseq-place" ,128,6,0,0,0,0,0,0,NULL,&atx[14],&atx[101],0,&atx[105]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[102],NULL,0,NULL} ,
  {422, "ID2S-Bioseq-place-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[103],0,&atx[107]} ,
  {0, "bioseq-set" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[104]} ,
  {0, "seq-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[68],NULL,0,NULL} ,
  {0, "feat-ids" ,128,7,0,0,0,0,0,0,NULL,&atx[14],&atx[106],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[107],NULL,0,NULL} ,
  {423, "ID2S-Seq-feat-Ids-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[108],0,&atx[85]} ,
  {0, "feat-types" ,128,0,0,1,0,0,0,0,NULL,&atx[14],&atx[109],0,&atx[110]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "xref-types" ,128,1,0,1,0,0,0,0,NULL,&atx[14],&atx[111],0,&atx[112]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "local-ids" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[113],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "skeleton" ,128,2,0,1,0,0,0,0,NULL,&atx[51],NULL,0,NULL} ,
  {427, "ID2S-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[116],0,&atx[118]} ,
  {0, "data" ,128,0,0,0,0,0,0,0,NULL,&atx[14],&atx[117],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[118],NULL,0,NULL} ,
  {428, "ID2S-Chunk-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[119],0,&atx[130]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[50],&atx[120],0,&atx[123]} ,
  {0, "bioseq-set" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[121]} ,
  {0, "gi" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[122]} ,
  {0, "seq-id" ,128,2,0,0,0,0,0,0,NULL,&atx[20],NULL,0,NULL} ,
  {0, "descr" ,128,1,0,1,0,0,0,0,NULL,&atx[54],NULL,0,&atx[124]} ,
  {0, "annots" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[125],0,&atx[126]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[53],NULL,0,NULL} ,
  {0, "assembly" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[127],0,&atx[128]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {0, "seq-map" ,128,4,0,1,0,0,0,0,NULL,&atx[134],&atx[129],0,&atx[135]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[130],NULL,0,NULL} ,
  {429, "ID2S-Sequence-Piece" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[131],0,&atx[22]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[132]} ,
  {0, "data" ,128,1,0,0,0,0,0,0,NULL,&atx[134],&atx[133],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[55],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seq-data" ,128,5,0,1,0,0,0,0,NULL,&atx[134],&atx[136],0,&atx[137]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[130],NULL,0,NULL} ,
  {0, "bioseqs" ,128,6,0,1,0,0,0,0,NULL,&atx[14],&atx[138],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[52],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Seq-split" , "seqsplit.h19",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Seq-split
*
**************************************************/

#define ID2S_CHUNK_ID &at[0]

#define ID2S_SEQ_ANNOT_INFO &at[2]
#define ID2S_SEQ_ANNOT_INFO_name &at[3]
#define ID2S_SEQ_ANNOT_INFO_align &at[5]
#define ID2S_SEQ_ANNOT_INFO_graph &at[7]
#define ID2S_SEQ_ANNOT_INFO_feat &at[8]
#define ID2S_SEQ_ANNOT_INFO_feat_E &at[9]
#define ID2S_SEQ_ANNOT_INFO_seq_loc &at[16]

#define ID2S_SPLIT_INFO &at[58]
#define ID2S_SPLIT_INFO_bioseqs_info &at[59]
#define ID2S_SPLIT_INFO_bioseqs_info_E &at[60]
#define ID2S_SPLIT_INFO_chunks &at[73]
#define ID2S_SPLIT_INFO_chunks_E &at[74]
#define ID2S_SPLIT_INFO_skeleton &at[114]

#define ID2S_BIOSEQS_INFO &at[61]
#define ID2S_BIOSEQS_INFO_info &at[62]
#define ID2S_BIOSEQS_INFO_bioseqs &at[67]

#define ID2S_CHUNK_INFO &at[75]
#define ID2S_CHUNK_INFO_id &at[76]
#define ID2S_CHUNK_INFO_content &at[77]
#define ID2S_CHUNK_INFO_content_E &at[78]

#define ID2S_BIOSEQ_INFO &at[63]
#define ID2S_BIOSEQ_INFO_gap_count &at[64]
#define BIOSEQ_INFO_seq_map_has_ref &at[65]

#define ID2S_BIOSEQ_IDS &at[68]
#define ID2S_BIOSEQ_IDS_E &at[69]
#define ID2S_BIOSEQ_IDS_E_gi &at[70]
#define ID2S_BIOSEQ_IDS_E_seq_id &at[71]
#define ID2S_BIOSEQ_IDS_E_gi_range &at[72]

#define ID2S_CHUNK_CONTENT &at[79]
#define ID2S_CHUNK_CONTENT_seq_descr &at[80]
#define ID2S_CHUNK_CONTENT_seq_annot &at[87]
#define ID2S_CHUNK_CONTENT_seq_assembly &at[88]
#define ID2S_CHUNK_CONTENT_seq_map &at[91]
#define ID2S_CHUNK_CONTENT_seq_data &at[93]
#define CHUNK_CONTENT_seq_annot_place &at[95]
#define ID2S_CHUNK_CONTENT_bioseq_place &at[100]
#define CHUNK_CONTENT_bioseq_place_E &at[101]
#define ID2S_CHUNK_CONTENT_feat_ids &at[105]
#define ID2S_CHUNK_CONTENT_feat_ids_E &at[106]

#define ID2S_SEQ_DESCR_INFO &at[81]
#define ID2S_SEQ_DESCR_INFO_type_mask &at[82]
#define ID2S_SEQ_DESCR_INFO_bioseqs &at[83]
#define ID2S_SEQ_DESCR_INFO_bioseq_sets &at[84]

#define ID2S_SEQ_ASSEMBLY_INFO &at[89]
#define ID2S_SEQ_ASSEMBLY_INFO_bioseqs &at[90]

#define ID2S_SEQ_MAP_INFO &at[92]

#define ID2S_SEQ_DATA_INFO &at[94]

#define ID2S_SEQ_ANNOT_PLACE_INFO &at[96]
#define ID2S_SEQ_ANNOT_PLACE_INFO_name &at[97]
#define SEQ_ANNOT_PLACE_INFO_bioseqs &at[98]
#define ANNOT_PLACE_INFO_bioseq_sets &at[99]

#define ID2S_BIOSEQ_PLACE_INFO &at[102]
#define BIOSEQ_PLACE_INFO_bioseq_set &at[103]
#define ID2S_BIOSEQ_PLACE_INFO_seq_ids &at[104]

#define ID2S_SEQ_FEAT_IDS_INFO &at[107]
#define SEQ_FEAT_IDS_INFO_feat_types &at[108]
#define SEQ_FEAT_IDS_INFO_feat_types_E &at[109]
#define SEQ_FEAT_IDS_INFO_xref_types &at[110]
#define SEQ_FEAT_IDS_INFO_xref_types_E &at[111]
#define SEQ_FEAT_IDS_INFO_local_ids &at[112]
#define SEQ_FEAT_IDS_INFO_local_ids_E &at[113]

#define ID2S_BIOSEQ_SET_IDS &at[85]
#define ID2S_BIOSEQ_SET_IDS_E &at[86]

#define ID2S_FEAT_TYPE_INFO &at[10]
#define ID2S_FEAT_TYPE_INFO_type &at[11]
#define ID2S_FEAT_TYPE_INFO_subtypes &at[12]
#define ID2S_FEAT_TYPE_INFO_subtypes_E &at[13]

#define ID2S_SEQ_LOC &at[17]
#define ID2S_SEQ_LOC_whole_gi &at[18]
#define ID2S_SEQ_LOC_whole_seq_id &at[19]
#define ID2S_SEQ_LOC_whole_gi_range &at[21]
#define ID2S_SEQ_LOC_gi_interval &at[25]
#define ID2S_SEQ_LOC_seq_id_interval &at[30]
#define ID2S_SEQ_LOC_gi_ints &at[35]
#define ID2S_SEQ_LOC_seq_id_ints &at[43]
#define ID2S_SEQ_LOC_loc_set &at[48]
#define ID2S_SEQ_LOC_loc_set_E &at[49]

#define ID2S_CHUNK &at[115]
#define ID2S_CHUNK_data &at[116]
#define ID2S_CHUNK_data_E &at[117]

#define ID2S_CHUNK_DATA &at[118]
#define ID2S_CHUNK_DATA_id &at[119]
#define ID2S_CHUNK_DATA_id_bioseq_set &at[120]
#define ID2S_CHUNK_DATA_id_gi &at[121]
#define ID2S_CHUNK_DATA_id_seq_id &at[122]
#define ID2S_CHUNK_DATA_descr &at[123]
#define ID2S_CHUNK_DATA_annots &at[124]
#define ID2S_CHUNK_DATA_annots_E &at[125]
#define ID2S_CHUNK_DATA_assembly &at[126]
#define ID2S_CHUNK_DATA_assembly_E &at[127]
#define ID2S_CHUNK_DATA_seq_map &at[128]
#define ID2S_CHUNK_DATA_seq_map_E &at[129]
#define ID2S_CHUNK_DATA_seq_data &at[135]
#define ID2S_CHUNK_DATA_seq_data_E &at[136]
#define ID2S_CHUNK_DATA_bioseqs &at[137]
#define ID2S_CHUNK_DATA_bioseqs_E &at[138]

#define ID2S_SEQUENCE_PIECE &at[130]
#define ID2S_SEQUENCE_PIECE_start &at[131]
#define ID2S_SEQUENCE_PIECE_data &at[132]
#define ID2S_SEQUENCE_PIECE_data_E &at[133]

#define ID2S_GI_RANGE &at[22]
#define ID2S_GI_RANGE_start &at[23]
#define ID2S_GI_RANGE_count &at[24]

#define ID2S_GI_INTERVAL &at[26]
#define ID2S_GI_INTERVAL_gi &at[27]
#define ID2S_GI_INTERVAL_start &at[28]
#define ID2S_GI_INTERVAL_length &at[29]

#define ID2S_SEQ_ID_INTERVAL &at[31]
#define ID2S_SEQ_ID_INTERVAL_seq_id &at[32]
#define ID2S_SEQ_ID_INTERVAL_start &at[33]
#define ID2S_SEQ_ID_INTERVAL_length &at[34]

#define ID2S_GI_INTS &at[36]
#define ID2S_GI_INTS_gi &at[37]
#define ID2S_GI_INTS_ints &at[38]
#define ID2S_GI_INTS_ints_E &at[39]

#define ID2S_SEQ_ID_INTS &at[44]
#define ID2S_SEQ_ID_INTS_seq_id &at[45]
#define ID2S_SEQ_ID_INTS_ints &at[46]
#define ID2S_SEQ_ID_INTS_ints_E &at[47]

#define ID2S_INTERVAL &at[40]
#define ID2S_INTERVAL_start &at[41]
#define ID2S_INTERVAL_length &at[42]
