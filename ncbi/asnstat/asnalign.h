/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnalign.h64";
static AsnValxNode avnx[13] = {
    {20,"not-set" ,0,0.0,&avnx[1] } ,
    {20,"global" ,1,0.0,&avnx[2] } ,
    {20,"diags" ,2,0.0,&avnx[3] } ,
    {20,"partial" ,3,0.0,&avnx[4] } ,
    {20,"disc" ,4,0.0,&avnx[5] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {20,"transcript" ,0,0.0,&avnx[11] } ,
    {20,"protein" ,1,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } };

static AsnType atx[161] = {
  {401, "Seq-align" ,1,0,0,0,0,1,0,0,NULL,&atx[15],&atx[1],0,&atx[7]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[2],&avnx[0],0,&atx[3]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "dim" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "score" ,128,2,0,1,0,0,0,0,NULL,&atx[16],&atx[6],0,&atx[17]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {402, "Score" ,1,0,0,0,0,1,0,0,NULL,&atx[15],&atx[8],0,&atx[111]} ,
  {0, "id" ,128,0,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[10]} ,
  {409, "Object-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[20]} ,
  {0, "value" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[11],0,NULL} ,
  {0, "real" ,128,0,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[13]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "int" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "segs" ,128,3,0,0,0,0,0,0,NULL,&atx[14],&atx[18],0,&atx[155]} ,
  {0, "dendiag" ,128,0,0,0,0,0,0,0,NULL,&atx[25],&atx[19],0,&atx[34]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,NULL} ,
  {410, "Dense-diag" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[21],0,&atx[35]} ,
  {0, "dim" ,128,0,0,0,1,0,0,0,&avnx[6],&atx[4],NULL,0,&atx[22]} ,
  {0, "ids" ,128,1,0,0,0,0,0,0,NULL,&atx[25],&atx[23],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {405, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[56]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "starts" ,128,2,0,0,0,0,0,0,NULL,&atx[25],&atx[27],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "len" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[29]} ,
  {0, "strands" ,128,4,0,1,0,0,0,0,NULL,&atx[25],&atx[30],0,&atx[32]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,NULL} ,
  {407, "Na-strand" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[122]} ,
  {0, "scores" ,128,5,0,1,0,0,0,0,NULL,&atx[16],&atx[33],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "denseg" ,128,1,0,0,0,0,0,0,NULL,&atx[35],NULL,0,&atx[48]} ,
  {411, "Dense-seg" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[36],0,&atx[50]} ,
  {0, "dim" ,128,0,0,0,1,0,0,0,&avnx[7],&atx[4],NULL,0,&atx[37]} ,
  {0, "numseg" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[38]} ,
  {0, "ids" ,128,2,0,0,0,0,0,0,NULL,&atx[25],&atx[39],0,&atx[40]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "starts" ,128,3,0,0,0,0,0,0,NULL,&atx[25],&atx[41],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "lens" ,128,4,0,0,0,0,0,0,NULL,&atx[25],&atx[43],0,&atx[44]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "strands" ,128,5,0,1,0,0,0,0,NULL,&atx[25],&atx[45],0,&atx[46]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,NULL} ,
  {0, "scores" ,128,6,0,1,0,0,0,0,NULL,&atx[25],&atx[47],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "std" ,128,2,0,0,0,0,0,0,NULL,&atx[25],&atx[49],0,&atx[59]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {412, "Std-seg" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[51],0,&atx[60]} ,
  {0, "dim" ,128,0,0,0,1,0,0,0,&avnx[8],&atx[4],NULL,0,&atx[52]} ,
  {0, "ids" ,128,1,0,1,0,0,0,0,NULL,&atx[25],&atx[53],0,&atx[54]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "loc" ,128,2,0,0,0,0,0,0,NULL,&atx[25],&atx[55],0,&atx[57]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {406, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[31]} ,
  {0, "scores" ,128,3,0,1,0,0,0,0,NULL,&atx[16],&atx[58],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "packed" ,128,3,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[75]} ,
  {413, "Packed-seg" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[61],0,&atx[79]} ,
  {0, "dim" ,128,0,0,0,1,0,0,0,&avnx[9],&atx[4],NULL,0,&atx[62]} ,
  {0, "numseg" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[63]} ,
  {0, "ids" ,128,2,0,0,0,0,0,0,NULL,&atx[25],&atx[64],0,&atx[65]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "starts" ,128,3,0,0,0,0,0,0,NULL,&atx[25],&atx[66],0,&atx[67]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "present" ,128,4,0,0,0,0,0,0,NULL,&atx[68],NULL,0,&atx[69]} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "lens" ,128,5,0,0,0,0,0,0,NULL,&atx[25],&atx[70],0,&atx[71]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "strands" ,128,6,0,1,0,0,0,0,NULL,&atx[25],&atx[72],0,&atx[73]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,NULL} ,
  {0, "scores" ,128,7,0,1,0,0,0,0,NULL,&atx[25],&atx[74],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "disc" ,128,4,0,0,0,0,0,0,NULL,&atx[76],NULL,0,&atx[78]} ,
  {404, "Seq-align-set" ,1,0,0,0,0,1,0,0,NULL,&atx[16],&atx[77],0,&atx[24]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {0, "spliced" ,128,5,0,0,0,0,0,0,NULL,&atx[79],NULL,0,&atx[130]} ,
  {414, "Spliced-seg" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[80],0,&atx[131]} ,
  {0, "product-id" ,128,0,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[81]} ,
  {0, "genomic-id" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[82]} ,
  {0, "product-strand" ,128,2,0,1,0,0,0,0,NULL,&atx[31],NULL,0,&atx[83]} ,
  {0, "genomic-strand" ,128,3,0,1,0,0,0,0,NULL,&atx[31],NULL,0,&atx[84]} ,
  {0, "product-type" ,128,4,0,0,0,0,0,0,NULL,&atx[2],&avnx[10],0,&atx[85]} ,
  {0, "exons" ,128,5,0,0,0,0,0,0,NULL,&atx[25],&atx[86],0,&atx[123]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[87],NULL,0,NULL} ,
  {416, "Spliced-exon" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[88],0,&atx[127]} ,
  {0, "product-start" ,128,0,0,0,0,0,0,0,NULL,&atx[89],NULL,0,&atx[95]} ,
  {418, "Product-pos" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[90],0,&atx[104]} ,
  {0, "nucpos" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[91]} ,
  {0, "protpos" ,128,1,0,0,0,0,0,0,NULL,&atx[92],NULL,0,NULL} ,
  {421, "Prot-pos" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[93],0,&atx[135]} ,
  {0, "amin" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[94]} ,
  {0, "frame" ,128,1,0,0,1,0,0,0,&avnx[12],&atx[4],NULL,0,NULL} ,
  {0, "product-end" ,128,1,0,0,0,0,0,0,NULL,&atx[89],NULL,0,&atx[96]} ,
  {0, "genomic-start" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[97]} ,
  {0, "genomic-end" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[98]} ,
  {0, "product-id" ,128,4,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[99]} ,
  {0, "genomic-id" ,128,5,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[100]} ,
  {0, "product-strand" ,128,6,0,1,0,0,0,0,NULL,&atx[31],NULL,0,&atx[101]} ,
  {0, "genomic-strand" ,128,7,0,1,0,0,0,0,NULL,&atx[31],NULL,0,&atx[102]} ,
  {0, "parts" ,128,8,0,1,0,0,0,0,NULL,&atx[25],&atx[103],0,&atx[110]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[104],NULL,0,NULL} ,
  {419, "Spliced-exon-chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[105],0,&atx[114]} ,
  {0, "match" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[106]} ,
  {0, "mismatch" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[107]} ,
  {0, "diag" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[108]} ,
  {0, "product-ins" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[109]} ,
  {0, "genomic-ins" ,128,4,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "scores" ,128,9,0,1,0,0,0,0,NULL,&atx[111],NULL,0,&atx[113]} ,
  {403, "Score-set" ,1,0,0,0,0,1,0,0,NULL,&atx[16],&atx[112],0,&atx[76]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "acceptor-before-exon" ,128,10,0,1,0,0,0,0,NULL,&atx[114],NULL,0,&atx[117]} ,
  {420, "Splice-site" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[115],0,&atx[92]} ,
  {0, "bases" ,128,0,0,0,0,0,0,0,NULL,&atx[116],NULL,0,NULL} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "donor-after-exon" ,128,11,0,1,0,0,0,0,NULL,&atx[114],NULL,0,&atx[118]} ,
  {0, "partial" ,128,12,0,1,0,0,0,0,NULL,&atx[119],NULL,0,&atx[120]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ext" ,128,13,0,1,0,0,0,0,NULL,&atx[25],&atx[121],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[122],NULL,0,NULL} ,
  {408, "User-object" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[9]} ,
  {0, "poly-a" ,128,6,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[124]} ,
  {0, "product-length" ,128,7,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[125]} ,
  {0, "modifiers" ,128,8,0,1,0,0,0,0,NULL,&atx[16],&atx[126],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[127],NULL,0,NULL} ,
  {417, "Spliced-seg-modifier" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[128],0,&atx[89]} ,
  {0, "start-codon-found" ,128,0,0,0,0,0,0,0,NULL,&atx[119],NULL,0,&atx[129]} ,
  {0, "stop-codon-found" ,128,1,0,0,0,0,0,0,NULL,&atx[119],NULL,0,NULL} ,
  {0, "sparse" ,128,6,0,0,0,0,0,0,NULL,&atx[131],NULL,0,NULL} ,
  {415, "Sparse-seg" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[132],0,&atx[87]} ,
  {0, "master-id" ,128,0,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[133]} ,
  {0, "rows" ,128,1,0,0,0,0,0,0,NULL,&atx[16],&atx[134],0,&atx[149]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[135],NULL,0,NULL} ,
  {422, "Sparse-align" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[136],0,&atx[153]} ,
  {0, "first-id" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[137]} ,
  {0, "second-id" ,128,1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[138]} ,
  {0, "numseg" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[139]} ,
  {0, "first-starts" ,128,3,0,0,0,0,0,0,NULL,&atx[25],&atx[140],0,&atx[141]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "second-starts" ,128,4,0,0,0,0,0,0,NULL,&atx[25],&atx[142],0,&atx[143]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "lens" ,128,5,0,0,0,0,0,0,NULL,&atx[25],&atx[144],0,&atx[145]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "second-strands" ,128,6,0,1,0,0,0,0,NULL,&atx[25],&atx[146],0,&atx[147]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,NULL} ,
  {0, "seg-scores" ,128,7,0,1,0,0,0,0,NULL,&atx[16],&atx[148],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "row-scores" ,128,2,0,1,0,0,0,0,NULL,&atx[16],&atx[150],0,&atx[151]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {0, "ext" ,128,3,0,1,0,0,0,0,NULL,&atx[16],&atx[152],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[153],NULL,0,NULL} ,
  {423, "Sparse-seg-ext" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[154],0,NULL} ,
  {0, "index" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "bounds" ,128,4,0,1,0,0,0,0,NULL,&atx[16],&atx[156],0,&atx[157]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {0, "id" ,128,5,0,1,0,0,0,0,NULL,&atx[25],&atx[158],0,&atx[159]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "ext" ,128,6,0,1,0,0,0,0,NULL,&atx[25],&atx[160],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[122],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Seqalign" , "asnalign.h64",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Seqalign
*
**************************************************/

#define SEQ_ALIGN &at[0]
#define SEQ_ALIGN_type &at[1]
#define SEQ_ALIGN_dim &at[3]
#define SEQ_ALIGN_score &at[5]
#define SEQ_ALIGN_score_E &at[6]
#define SEQ_ALIGN_segs &at[17]
#define SEQ_ALIGN_segs_dendiag &at[18]
#define SEQ_ALIGN_segs_dendiag_E &at[19]
#define SEQ_ALIGN_segs_denseg &at[34]
#define SEQ_ALIGN_segs_std &at[48]
#define SEQ_ALIGN_segs_std_E &at[49]
#define SEQ_ALIGN_segs_packed &at[59]
#define SEQ_ALIGN_segs_disc &at[75]
#define SEQ_ALIGN_segs_spliced &at[78]
#define SEQ_ALIGN_segs_sparse &at[130]
#define SEQ_ALIGN_bounds &at[155]
#define SEQ_ALIGN_bounds_E &at[156]
#define SEQ_ALIGN_id &at[157]
#define SEQ_ALIGN_id_E &at[158]
#define SEQ_ALIGN_ext &at[159]
#define SEQ_ALIGN_ext_E &at[160]

#define SCORE &at[7]
#define SCORE_id &at[8]
#define SCORE_value &at[10]
#define SCORE_value_real &at[11]
#define SCORE_value_int &at[13]

#define SCORE_SET &at[111]
#define SCORE_SET_E &at[112]

#define SEQ_ALIGN_SET &at[76]
#define SEQ_ALIGN_SET_E &at[77]

#define DENSE_DIAG &at[20]
#define DENSE_DIAG_dim &at[21]
#define DENSE_DIAG_ids &at[22]
#define DENSE_DIAG_ids_E &at[23]
#define DENSE_DIAG_starts &at[26]
#define DENSE_DIAG_starts_E &at[27]
#define DENSE_DIAG_len &at[28]
#define DENSE_DIAG_strands &at[29]
#define DENSE_DIAG_strands_E &at[30]
#define DENSE_DIAG_scores &at[32]
#define DENSE_DIAG_scores_E &at[33]

#define DENSE_SEG &at[35]
#define DENSE_SEG_dim &at[36]
#define DENSE_SEG_numseg &at[37]
#define DENSE_SEG_ids &at[38]
#define DENSE_SEG_ids_E &at[39]
#define DENSE_SEG_starts &at[40]
#define DENSE_SEG_starts_E &at[41]
#define DENSE_SEG_lens &at[42]
#define DENSE_SEG_lens_E &at[43]
#define DENSE_SEG_strands &at[44]
#define DENSE_SEG_strands_E &at[45]
#define DENSE_SEG_scores &at[46]
#define DENSE_SEG_scores_E &at[47]

#define STD_SEG &at[50]
#define STD_SEG_dim &at[51]
#define STD_SEG_ids &at[52]
#define STD_SEG_ids_E &at[53]
#define STD_SEG_loc &at[54]
#define STD_SEG_loc_E &at[55]
#define STD_SEG_scores &at[57]
#define STD_SEG_scores_E &at[58]

#define PACKED_SEG &at[60]
#define PACKED_SEG_dim &at[61]
#define PACKED_SEG_numseg &at[62]
#define PACKED_SEG_ids &at[63]
#define PACKED_SEG_ids_E &at[64]
#define PACKED_SEG_starts &at[65]
#define PACKED_SEG_starts_E &at[66]
#define PACKED_SEG_present &at[67]
#define PACKED_SEG_lens &at[69]
#define PACKED_SEG_lens_E &at[70]
#define PACKED_SEG_strands &at[71]
#define PACKED_SEG_strands_E &at[72]
#define PACKED_SEG_scores &at[73]
#define PACKED_SEG_scores_E &at[74]

#define SPLICED_SEG &at[79]
#define SPLICED_SEG_product_id &at[80]
#define SPLICED_SEG_genomic_id &at[81]
#define SPLICED_SEG_product_strand &at[82]
#define SPLICED_SEG_genomic_strand &at[83]
#define SPLICED_SEG_product_type &at[84]
#define SPLICED_SEG_exons &at[85]
#define SPLICED_SEG_exons_E &at[86]
#define SPLICED_SEG_poly_a &at[123]
#define SPLICED_SEG_product_length &at[124]
#define SPLICED_SEG_modifiers &at[125]
#define SPLICED_SEG_modifiers_E &at[126]

#define SPARSE_SEG &at[131]
#define SPARSE_SEG_master_id &at[132]
#define SPARSE_SEG_rows &at[133]
#define SPARSE_SEG_rows_E &at[134]
#define SPARSE_SEG_row_scores &at[149]
#define SPARSE_SEG_row_scores_E &at[150]
#define SPARSE_SEG_ext &at[151]
#define SPARSE_SEG_ext_E &at[152]

#define SPLICED_EXON &at[87]
#define SPLICED_EXON_product_start &at[88]
#define SPLICED_EXON_product_end &at[95]
#define SPLICED_EXON_genomic_start &at[96]
#define SPLICED_EXON_genomic_end &at[97]
#define SPLICED_EXON_product_id &at[98]
#define SPLICED_EXON_genomic_id &at[99]
#define SPLICED_EXON_product_strand &at[100]
#define SPLICED_EXON_genomic_strand &at[101]
#define SPLICED_EXON_parts &at[102]
#define SPLICED_EXON_parts_E &at[103]
#define SPLICED_EXON_scores &at[110]
#define EXON_acceptor_before_exon &at[113]
#define SPLICED_EXON_donor_after_exon &at[117]
#define SPLICED_EXON_partial &at[118]
#define SPLICED_EXON_ext &at[120]
#define SPLICED_EXON_ext_E &at[121]

#define SPLICED_SEG_MODIFIER &at[127]
#define SEG_MODIFIER_start_codon_found &at[128]
#define SEG_MODIFIER_stop_codon_found &at[129]

#define PRODUCT_POS &at[89]
#define PRODUCT_POS_nucpos &at[90]
#define PRODUCT_POS_protpos &at[91]

#define SPLICED_EXON_CHUNK &at[104]
#define SPLICED_EXON_CHUNK_match &at[105]
#define SPLICED_EXON_CHUNK_mismatch &at[106]
#define SPLICED_EXON_CHUNK_diag &at[107]
#define SPLICED_EXON_CHUNK_product_ins &at[108]
#define SPLICED_EXON_CHUNK_genomic_ins &at[109]

#define SPLICE_SITE &at[114]
#define SPLICE_SITE_bases &at[115]

#define PROT_POS &at[92]
#define PROT_POS_amin &at[93]
#define PROT_POS_frame &at[94]

#define SPARSE_ALIGN &at[135]
#define SPARSE_ALIGN_first_id &at[136]
#define SPARSE_ALIGN_second_id &at[137]
#define SPARSE_ALIGN_numseg &at[138]
#define SPARSE_ALIGN_first_starts &at[139]
#define SPARSE_ALIGN_first_starts_E &at[140]
#define SPARSE_ALIGN_second_starts &at[141]
#define SPARSE_ALIGN_second_starts_E &at[142]
#define SPARSE_ALIGN_lens &at[143]
#define SPARSE_ALIGN_lens_E &at[144]
#define SPARSE_ALIGN_second_strands &at[145]
#define SPARSE_ALIGN_second_strands_E &at[146]
#define SPARSE_ALIGN_seg_scores &at[147]
#define SPARSE_ALIGN_seg_scores_E &at[148]

#define SPARSE_SEG_EXT &at[153]
#define SPARSE_SEG_EXT_index &at[154]
