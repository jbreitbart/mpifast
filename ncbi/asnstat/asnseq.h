/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnseq.h20";
static AsnValxNode avnx[164] = {
    {3,NULL,1,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[4] } ,
    {20,"sources" ,1,0.0,&avnx[5] } ,
    {20,"aligns" ,2,0.0,NULL } ,
    {20,"seq" ,0,0.0,&avnx[7] } ,
    {20,"sites" ,1,0.0,&avnx[8] } ,
    {20,"feats" ,2,0.0,&avnx[9] } ,
    {20,"no-target" ,3,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"ref" ,1,0.0,&avnx[12] } ,
    {20,"alt" ,2,0.0,&avnx[13] } ,
    {20,"blocks" ,3,0.0,&avnx[14] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[16] } ,
    {20,"genomic" ,1,0.0,&avnx[17] } ,
    {20,"pre-mRNA" ,2,0.0,&avnx[18] } ,
    {20,"mRNA" ,3,0.0,&avnx[19] } ,
    {20,"rRNA" ,4,0.0,&avnx[20] } ,
    {20,"tRNA" ,5,0.0,&avnx[21] } ,
    {20,"snRNA" ,6,0.0,&avnx[22] } ,
    {20,"scRNA" ,7,0.0,&avnx[23] } ,
    {20,"peptide" ,8,0.0,&avnx[24] } ,
    {20,"other-genetic" ,9,0.0,&avnx[25] } ,
    {20,"genomic-mRNA" ,10,0.0,&avnx[26] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"dna" ,0,0.0,&avnx[28] } ,
    {20,"rna" ,1,0.0,&avnx[29] } ,
    {20,"extrachrom" ,2,0.0,&avnx[30] } ,
    {20,"plasmid" ,3,0.0,&avnx[31] } ,
    {20,"mitochondrial" ,4,0.0,&avnx[32] } ,
    {20,"chloroplast" ,5,0.0,&avnx[33] } ,
    {20,"kinetoplast" ,6,0.0,&avnx[34] } ,
    {20,"cyanelle" ,7,0.0,&avnx[35] } ,
    {20,"synthetic" ,8,0.0,&avnx[36] } ,
    {20,"recombinant" ,9,0.0,&avnx[37] } ,
    {20,"partial" ,10,0.0,&avnx[38] } ,
    {20,"complete" ,11,0.0,&avnx[39] } ,
    {20,"mutagen" ,12,0.0,&avnx[40] } ,
    {20,"natmut" ,13,0.0,&avnx[41] } ,
    {20,"transposon" ,14,0.0,&avnx[42] } ,
    {20,"insertion-seq" ,15,0.0,&avnx[43] } ,
    {20,"no-left" ,16,0.0,&avnx[44] } ,
    {20,"no-right" ,17,0.0,&avnx[45] } ,
    {20,"macronuclear" ,18,0.0,&avnx[46] } ,
    {20,"proviral" ,19,0.0,&avnx[47] } ,
    {20,"est" ,20,0.0,&avnx[48] } ,
    {20,"sts" ,21,0.0,&avnx[49] } ,
    {20,"survey" ,22,0.0,&avnx[50] } ,
    {20,"chromoplast" ,23,0.0,&avnx[51] } ,
    {20,"genemap" ,24,0.0,&avnx[52] } ,
    {20,"restmap" ,25,0.0,&avnx[53] } ,
    {20,"physmap" ,26,0.0,&avnx[54] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"concept-trans" ,1,0.0,&avnx[56] } ,
    {20,"seq-pept" ,2,0.0,&avnx[57] } ,
    {20,"both" ,3,0.0,&avnx[58] } ,
    {20,"seq-pept-overlap" ,4,0.0,&avnx[59] } ,
    {20,"seq-pept-homol" ,5,0.0,&avnx[60] } ,
    {20,"concept-trans-a" ,6,0.0,&avnx[61] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[63] } ,
    {20,"genomic" ,1,0.0,&avnx[64] } ,
    {20,"pre-RNA" ,2,0.0,&avnx[65] } ,
    {20,"mRNA" ,3,0.0,&avnx[66] } ,
    {20,"rRNA" ,4,0.0,&avnx[67] } ,
    {20,"tRNA" ,5,0.0,&avnx[68] } ,
    {20,"snRNA" ,6,0.0,&avnx[69] } ,
    {20,"scRNA" ,7,0.0,&avnx[70] } ,
    {20,"peptide" ,8,0.0,&avnx[71] } ,
    {20,"other-genetic" ,9,0.0,&avnx[72] } ,
    {20,"genomic-mRNA" ,10,0.0,&avnx[73] } ,
    {20,"cRNA" ,11,0.0,&avnx[74] } ,
    {20,"snoRNA" ,12,0.0,&avnx[75] } ,
    {20,"transcribed-RNA" ,13,0.0,&avnx[76] } ,
    {20,"ncRNA" ,14,0.0,&avnx[77] } ,
    {20,"tmRNA" ,15,0.0,&avnx[78] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[81] } ,
    {20,"standard" ,1,0.0,&avnx[82] } ,
    {20,"est" ,2,0.0,&avnx[83] } ,
    {20,"sts" ,3,0.0,&avnx[84] } ,
    {20,"survey" ,4,0.0,&avnx[85] } ,
    {20,"genemap" ,5,0.0,&avnx[86] } ,
    {20,"physmap" ,6,0.0,&avnx[87] } ,
    {20,"derived" ,7,0.0,&avnx[88] } ,
    {20,"concept-trans" ,8,0.0,&avnx[89] } ,
    {20,"seq-pept" ,9,0.0,&avnx[90] } ,
    {20,"both" ,10,0.0,&avnx[91] } ,
    {20,"seq-pept-overlap" ,11,0.0,&avnx[92] } ,
    {20,"seq-pept-homol" ,12,0.0,&avnx[93] } ,
    {20,"concept-trans-a" ,13,0.0,&avnx[94] } ,
    {20,"htgs-1" ,14,0.0,&avnx[95] } ,
    {20,"htgs-2" ,15,0.0,&avnx[96] } ,
    {20,"htgs-3" ,16,0.0,&avnx[97] } ,
    {20,"fli-cdna" ,17,0.0,&avnx[98] } ,
    {20,"htgs-0" ,18,0.0,&avnx[99] } ,
    {20,"htc" ,19,0.0,&avnx[100] } ,
    {20,"wgs" ,20,0.0,&avnx[101] } ,
    {20,"barcode" ,21,0.0,&avnx[102] } ,
    {20,"composite-wgs-htgs" ,22,0.0,&avnx[103] } ,
    {20,"tsa" ,23,0.0,&avnx[104] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[107] } ,
    {20,"complete" ,1,0.0,&avnx[108] } ,
    {20,"partial" ,2,0.0,&avnx[109] } ,
    {20,"no-left" ,3,0.0,&avnx[110] } ,
    {20,"no-right" ,4,0.0,&avnx[111] } ,
    {20,"no-ends" ,5,0.0,&avnx[112] } ,
    {20,"has-left" ,6,0.0,&avnx[113] } ,
    {20,"has-right" ,7,0.0,&avnx[114] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[117] } ,
    {20,"virtual" ,1,0.0,&avnx[118] } ,
    {20,"raw" ,2,0.0,&avnx[119] } ,
    {20,"seg" ,3,0.0,&avnx[120] } ,
    {20,"const" ,4,0.0,&avnx[121] } ,
    {20,"ref" ,5,0.0,&avnx[122] } ,
    {20,"consen" ,6,0.0,&avnx[123] } ,
    {20,"map" ,7,0.0,&avnx[124] } ,
    {20,"delta" ,8,0.0,&avnx[125] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[127] } ,
    {20,"dna" ,1,0.0,&avnx[128] } ,
    {20,"rna" ,2,0.0,&avnx[129] } ,
    {20,"aa" ,3,0.0,&avnx[130] } ,
    {20,"na" ,4,0.0,&avnx[131] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[133] } ,
    {20,"linear" ,1,0.0,&avnx[134] } ,
    {20,"circular" ,2,0.0,&avnx[135] } ,
    {20,"tandem" ,3,0.0,&avnx[136] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[139] } ,
    {20,"ss" ,1,0.0,&avnx[140] } ,
    {20,"ds" ,2,0.0,&avnx[141] } ,
    {20,"mixed" ,3,0.0,&avnx[142] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[144] } ,
    {20,"fragment" ,1,0.0,&avnx[145] } ,
    {20,"clone" ,2,0.0,&avnx[146] } ,
    {20,"short-arm" ,3,0.0,&avnx[147] } ,
    {20,"heterochromatin" ,4,0.0,&avnx[148] } ,
    {20,"centromere" ,5,0.0,&avnx[149] } ,
    {20,"telomere" ,6,0.0,&avnx[150] } ,
    {20,"repeat" ,7,0.0,&avnx[151] } ,
    {20,"contig" ,8,0.0,&avnx[152] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"unlinked" ,0,0.0,&avnx[154] } ,
    {20,"linked" ,1,0.0,&avnx[155] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"genbank" ,1,0.0,&avnx[157] } ,
    {20,"embl" ,2,0.0,&avnx[158] } ,
    {20,"ddbj" ,3,0.0,&avnx[159] } ,
    {20,"pir" ,4,0.0,&avnx[160] } ,
    {20,"sp" ,5,0.0,&avnx[161] } ,
    {20,"bbone" ,6,0.0,&avnx[162] } ,
    {20,"pdb" ,7,0.0,&avnx[163] } ,
    {20,"other" ,255,0.0,NULL } };

static AsnType atx[219] = {
  {401, "Annotdesc" ,1,0,0,0,0,1,0,0,NULL,&atx[39],&atx[1],0,&atx[63]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "title" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[4]} ,
  {0, "comment" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[5]} ,
  {0, "pub" ,128,3,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[48]} ,
  {408, "Pubdesc" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[7],0,&atx[188]} ,
  {0, "pub" ,128,0,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[9]} ,
  {426, "Pub-equiv" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[82]} ,
  {0, "name" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[10]} ,
  {0, "fig" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[11]} ,
  {0, "num" ,128,3,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[40]} ,
  {407, "Numbering" ,1,0,0,0,0,1,0,0,NULL,&atx[39],&atx[13],0,&atx[6]} ,
  {0, "cont" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[21]} ,
  {440, "Num-cont" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[15],0,&atx[22]} ,
  {0, "refnum" ,128,0,0,0,1,0,0,0,&avnx[0],&atx[16],NULL,0,&atx[17]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "has-zero" ,128,1,0,0,1,0,0,0,&avnx[1],&atx[18],NULL,0,&atx[19]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ascending" ,128,2,0,0,1,0,0,0,&avnx[2],&atx[18],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "enum" ,128,1,0,0,0,0,0,0,NULL,&atx[22],NULL,0,&atx[27]} ,
  {441, "Num-enum" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[23],0,&atx[28]} ,
  {0, "num" ,128,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[24]} ,
  {0, "names" ,128,1,0,0,0,0,0,0,NULL,&atx[26],&atx[25],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ref" ,128,2,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[33]} ,
  {442, "Num-ref" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[29],0,&atx[34]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[30],&avnx[3],0,&atx[31]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "aligns" ,128,1,0,1,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {423, "Seq-align" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[162]} ,
  {0, "real" ,128,3,0,0,0,0,0,0,NULL,&atx[34],NULL,0,NULL} ,
  {443, "Num-real" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[35],0,&atx[155]} ,
  {0, "a" ,128,0,0,0,0,0,0,0,NULL,&atx[36],NULL,0,&atx[37]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "b" ,128,1,0,0,0,0,0,0,NULL,&atx[36],NULL,0,&atx[38]} ,
  {0, "units" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "numexc" ,128,4,0,1,0,0,0,0,NULL,&atx[18],NULL,0,&atx[41]} ,
  {0, "poly-a" ,128,5,0,1,0,0,0,0,NULL,&atx[18],NULL,0,&atx[42]} ,
  {0, "maploc" ,128,6,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[43]} ,
  {0, "seq-raw" ,128,7,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[45]} ,
  {351, "StringStore" ,64,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "align-group" ,128,8,0,1,0,0,0,0,NULL,&atx[16],NULL,0,&atx[46]} ,
  {0, "comment" ,128,9,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[47]} ,
  {0, "reftype" ,128,10,0,0,1,0,0,0,&avnx[10],&atx[16],&avnx[6],0,NULL} ,
  {0, "user" ,128,4,0,0,0,0,0,0,NULL,&atx[49],NULL,0,&atx[50]} ,
  {422, "User-object" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[32]} ,
  {0, "create-date" ,128,5,0,0,0,0,0,0,NULL,&atx[51],NULL,0,&atx[52]} ,
  {418, "Date" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[122]} ,
  {0, "update-date" ,128,6,0,0,0,0,0,0,NULL,&atx[51],NULL,0,&atx[53]} ,
  {0, "src" ,128,7,0,0,0,0,0,0,NULL,&atx[54],NULL,0,&atx[55]} ,
  {429, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[62]} ,
  {0, "align" ,128,8,0,0,0,0,0,0,NULL,&atx[56],NULL,0,&atx[61]} ,
  {462, "Align-def" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[57],0,NULL} ,
  {0, "align-type" ,128,0,0,0,0,0,0,0,NULL,&atx[16],&avnx[11],0,&atx[58]} ,
  {0, "ids" ,128,1,0,1,0,0,0,0,NULL,&atx[60],&atx[59],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "region" ,128,9,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} ,
  {430, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[90]} ,
  {402, "Annot-descr" ,1,0,0,0,0,1,0,0,NULL,&atx[60],&atx[64],0,&atx[65]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {403, "Bioseq" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[66],0,&atx[73]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[60],&atx[67],0,&atx[68]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {0, "descr" ,128,1,0,1,0,0,0,0,NULL,&atx[69],NULL,0,&atx[116]} ,
  {412, "Seq-descr" ,1,0,0,0,0,1,0,0,NULL,&atx[60],&atx[70],0,&atx[153]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[71],NULL,0,NULL} ,
  {411, "Seqdesc" ,1,0,0,0,0,1,0,0,NULL,&atx[39],&atx[72],0,&atx[69]} ,
  {0, "mol-type" ,128,0,0,0,0,0,0,0,NULL,&atx[73],NULL,0,&atx[74]} ,
  {404, "GIBB-mol" ,1,0,0,0,0,1,0,0,NULL,&atx[30],&avnx[15],0,&atx[106]} ,
  {0, "modif" ,128,1,0,0,0,0,0,0,NULL,&atx[60],&atx[75],0,&atx[77]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[76],NULL,0,NULL} ,
  {438, "GIBB-mod" ,1,0,0,0,0,0,0,0,NULL,&atx[30],&avnx[27],0,&atx[78]} ,
  {0, "method" ,128,2,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[79]} ,
  {439, "GIBB-method" ,1,0,0,0,0,0,0,0,NULL,&atx[30],&avnx[55],0,&atx[14]} ,
  {0, "name" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[80]} ,
  {0, "title" ,128,4,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[81]} ,
  {0, "org" ,128,5,0,0,0,0,0,0,NULL,&atx[82],NULL,0,&atx[83]} ,
  {427, "Org-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[108]} ,
  {0, "comment" ,128,6,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[84]} ,
  {0, "num" ,128,7,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[85]} ,
  {0, "maploc" ,128,8,0,0,0,0,0,0,NULL,&atx[86],NULL,0,&atx[87]} ,
  {420, "Dbtag" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[193]} ,
  {0, "pir" ,128,9,0,0,0,0,0,0,NULL,&atx[88],NULL,0,&atx[89]} ,
  {432, "PIR-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[98]} ,
  {0, "genbank" ,128,10,0,0,0,0,0,0,NULL,&atx[90],NULL,0,&atx[91]} ,
  {431, "GB-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[88]} ,
  {0, "pub" ,128,11,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[92]} ,
  {0, "region" ,128,12,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[93]} ,
  {0, "user" ,128,13,0,0,0,0,0,0,NULL,&atx[49],NULL,0,&atx[94]} ,
  {0, "sp" ,128,14,0,0,0,0,0,0,NULL,&atx[95],NULL,0,&atx[96]} ,
  {434, "SP-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[102]} ,
  {0, "dbxref" ,128,15,0,0,0,0,0,0,NULL,&atx[86],NULL,0,&atx[97]} ,
  {0, "embl" ,128,16,0,0,0,0,0,0,NULL,&atx[98],NULL,0,&atx[99]} ,
  {433, "EMBL-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[95]} ,
  {0, "create-date" ,128,17,0,0,0,0,0,0,NULL,&atx[51],NULL,0,&atx[100]} ,
  {0, "update-date" ,128,18,0,0,0,0,0,0,NULL,&atx[51],NULL,0,&atx[101]} ,
  {0, "prf" ,128,19,0,0,0,0,0,0,NULL,&atx[102],NULL,0,&atx[103]} ,
  {435, "PRF-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[104]} ,
  {0, "pdb" ,128,20,0,0,0,0,0,0,NULL,&atx[104],NULL,0,&atx[105]} ,
  {436, "PDB-block" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[218]} ,
  {0, "het" ,128,21,0,0,0,0,0,0,NULL,&atx[106],NULL,0,&atx[107]} ,
  {405, "Heterogen" ,1,0,0,0,0,1,0,0,NULL,&atx[2],NULL,0,&atx[110]} ,
  {0, "source" ,128,22,0,0,0,0,0,0,NULL,&atx[108],NULL,0,&atx[109]} ,
  {428, "BioSource" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[54]} ,
  {0, "molinfo" ,128,23,0,0,0,0,0,0,NULL,&atx[110],NULL,0,NULL} ,
  {406, "MolInfo" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[111],0,&atx[12]} ,
  {0, "biomol" ,128,0,0,0,1,0,0,0,&avnx[79],&atx[16],&avnx[62],0,&atx[112]} ,
  {0, "tech" ,128,1,0,0,1,0,0,0,&avnx[105],&atx[16],&avnx[80],0,&atx[113]} ,
  {0, "techexp" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[114]} ,
  {0, "completeness" ,128,3,0,0,1,0,0,0,&avnx[115],&atx[16],&avnx[106],0,&atx[115]} ,
  {0, "gbmoltype" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "inst" ,128,2,0,0,0,0,0,0,NULL,&atx[117],NULL,0,&atx[186]} ,
  {415, "Seq-inst" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[118],0,&atx[169]} ,
  {0, "repr" ,128,0,0,0,0,0,0,0,NULL,&atx[30],&avnx[116],0,&atx[119]} ,
  {0, "mol" ,128,1,0,0,0,0,0,0,NULL,&atx[30],&avnx[126],0,&atx[120]} ,
  {0, "length" ,128,2,0,1,0,0,0,0,NULL,&atx[16],NULL,0,&atx[121]} ,
  {0, "fuzz" ,128,3,0,1,0,0,0,0,NULL,&atx[122],NULL,0,&atx[123]} ,
  {419, "Int-fuzz" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[86]} ,
  {0, "topology" ,128,4,0,0,1,0,0,0,&avnx[137],&atx[30],&avnx[132],0,&atx[124]} ,
  {0, "strand" ,128,5,0,1,0,0,0,0,NULL,&atx[30],&avnx[138],0,&atx[125]} ,
  {0, "seq-data" ,128,6,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[152]} ,
  {410, "Seq-data" ,1,0,0,0,0,1,0,0,NULL,&atx[39],&atx[127],0,&atx[71]} ,
  {0, "iupacna" ,128,0,0,0,0,0,0,0,NULL,&atx[128],NULL,0,&atx[129]} ,
  {449, "IUPACna" ,1,0,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[130]} ,
  {0, "iupacaa" ,128,1,0,0,0,0,0,0,NULL,&atx[130],NULL,0,&atx[131]} ,
  {450, "IUPACaa" ,1,0,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[132]} ,
  {0, "ncbi2na" ,128,2,0,0,0,0,0,0,NULL,&atx[132],NULL,0,&atx[134]} ,
  {451, "NCBI2na" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[135]} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ncbi4na" ,128,3,0,0,0,0,0,0,NULL,&atx[135],NULL,0,&atx[136]} ,
  {452, "NCBI4na" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[137]} ,
  {0, "ncbi8na" ,128,4,0,0,0,0,0,0,NULL,&atx[137],NULL,0,&atx[138]} ,
  {453, "NCBI8na" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[139]} ,
  {0, "ncbipna" ,128,5,0,0,0,0,0,0,NULL,&atx[139],NULL,0,&atx[140]} ,
  {454, "NCBIpna" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[141]} ,
  {0, "ncbi8aa" ,128,6,0,0,0,0,0,0,NULL,&atx[141],NULL,0,&atx[142]} ,
  {455, "NCBI8aa" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[143]} ,
  {0, "ncbieaa" ,128,7,0,0,0,0,0,0,NULL,&atx[143],NULL,0,&atx[144]} ,
  {456, "NCBIeaa" ,1,0,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[145]} ,
  {0, "ncbipaa" ,128,8,0,0,0,0,0,0,NULL,&atx[145],NULL,0,&atx[146]} ,
  {457, "NCBIpaa" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[147]} ,
  {0, "ncbistdaa" ,128,9,0,0,0,0,0,0,NULL,&atx[147],NULL,0,&atx[148]} ,
  {458, "NCBIstdaa" ,1,0,0,0,0,0,0,0,NULL,&atx[133],NULL,0,&atx[149]} ,
  {0, "gap" ,128,10,0,0,0,0,0,0,NULL,&atx[149],NULL,0,NULL} ,
  {459, "Seq-gap" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[150],0,&atx[197]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[16],&avnx[143],0,&atx[151]} ,
  {0, "linkage" ,128,1,0,1,0,0,0,0,NULL,&atx[16],&avnx[153],0,NULL} ,
  {0, "ext" ,128,7,0,1,0,0,0,0,NULL,&atx[153],NULL,0,&atx[173]} ,
  {413, "Seq-ext" ,1,0,0,0,0,1,0,0,NULL,&atx[39],&atx[154],0,&atx[174]} ,
  {0, "seg" ,128,0,0,0,0,0,0,0,NULL,&atx[155],NULL,0,&atx[157]} ,
  {444, "Seg-ext" ,1,0,0,0,0,0,0,0,NULL,&atx[26],&atx[156],0,&atx[158]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} ,
  {0, "ref" ,128,1,0,0,0,0,0,0,NULL,&atx[158],NULL,0,&atx[159]} ,
  {445, "Ref-ext" ,1,0,0,0,0,0,0,0,NULL,&atx[62],NULL,0,&atx[160]} ,
  {0, "map" ,128,2,0,0,0,0,0,0,NULL,&atx[160],NULL,0,&atx[163]} ,
  {446, "Map-ext" ,1,0,0,0,0,0,0,0,NULL,&atx[26],&atx[161],0,&atx[166]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[162],NULL,0,NULL} ,
  {424, "Seq-feat" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[212]} ,
  {0, "delta" ,128,3,0,0,0,0,0,0,NULL,&atx[164],NULL,0,NULL} ,
  {417, "Delta-ext" ,1,0,0,0,0,1,0,0,NULL,&atx[26],&atx[165],0,&atx[51]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[166],NULL,0,NULL} ,
  {447, "Delta-seq" ,1,0,0,0,0,0,0,0,NULL,&atx[39],&atx[167],0,&atx[178]} ,
  {0, "loc" ,128,0,0,0,0,0,0,0,NULL,&atx[62],NULL,0,&atx[168]} ,
  {0, "literal" ,128,1,0,0,0,0,0,0,NULL,&atx[169],NULL,0,NULL} ,
  {416, "Seq-literal" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[170],0,&atx[164]} ,
  {0, "length" ,128,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[171]} ,
  {0, "fuzz" ,128,1,0,1,0,0,0,0,NULL,&atx[122],NULL,0,&atx[172]} ,
  {0, "seq-data" ,128,2,0,1,0,0,0,0,NULL,&atx[126],NULL,0,NULL} ,
  {0, "hist" ,128,8,0,1,0,0,0,0,NULL,&atx[174],NULL,0,NULL} ,
  {414, "Seq-hist" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[175],0,&atx[117]} ,
  {0, "assembly" ,128,0,0,1,0,0,0,0,NULL,&atx[60],&atx[176],0,&atx[177]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {0, "replaces" ,128,1,0,1,0,0,0,0,NULL,&atx[178],NULL,0,&atx[182]} ,
  {448, "Seq-hist-rec" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[179],0,&atx[128]} ,
  {0, "date" ,128,0,0,1,0,0,0,0,NULL,&atx[51],NULL,0,&atx[180]} ,
  {0, "ids" ,128,1,0,0,0,0,0,0,NULL,&atx[60],&atx[181],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {0, "replaced-by" ,128,2,0,1,0,0,0,0,NULL,&atx[178],NULL,0,&atx[183]} ,
  {0, "deleted" ,128,3,0,1,0,0,0,0,NULL,&atx[39],&atx[184],0,NULL} ,
  {0, "bool" ,128,0,0,0,0,0,0,0,NULL,&atx[18],NULL,0,&atx[185]} ,
  {0, "date" ,128,1,0,0,0,0,0,0,NULL,&atx[51],NULL,0,NULL} ,
  {0, "annot" ,128,3,0,1,0,0,0,0,NULL,&atx[60],&atx[187],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[188],NULL,0,NULL} ,
  {409, "Seq-annot" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[189],0,&atx[126]} ,
  {0, "id" ,128,0,0,1,0,0,0,0,NULL,&atx[60],&atx[190],0,&atx[202]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[191],NULL,0,NULL} ,
  {461, "Annot-id" ,1,0,0,0,0,0,0,0,NULL,&atx[39],&atx[192],0,&atx[56]} ,
  {0, "local" ,128,0,0,0,0,0,0,0,NULL,&atx[193],NULL,0,&atx[194]} ,
  {421, "Object-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[49]} ,
  {0, "ncbi" ,128,1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[195]} ,
  {0, "general" ,128,2,0,0,0,0,0,0,NULL,&atx[86],NULL,0,&atx[196]} ,
  {0, "other" ,128,3,0,0,0,0,0,0,NULL,&atx[197],NULL,0,NULL} ,
  {460, "Textannot-id" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[198],0,&atx[191]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[199]} ,
  {0, "accession" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[200]} ,
  {0, "release" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[201]} ,
  {0, "version" ,128,3,0,1,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {0, "db" ,128,1,0,1,0,0,0,0,NULL,&atx[16],&avnx[156],0,&atx[203]} ,
  {0, "name" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[204]} ,
  {0, "desc" ,128,3,0,1,0,0,0,0,NULL,&atx[63],NULL,0,&atx[205]} ,
  {0, "data" ,128,4,0,0,0,0,0,0,NULL,&atx[39],&atx[206],0,NULL} ,
  {0, "ftable" ,128,0,0,0,0,0,0,0,NULL,&atx[60],&atx[207],0,&atx[208]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[162],NULL,0,NULL} ,
  {0, "align" ,128,1,0,0,0,0,0,0,NULL,&atx[60],&atx[209],0,&atx[210]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {0, "graph" ,128,2,0,0,0,0,0,0,NULL,&atx[60],&atx[211],0,&atx[213]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[212],NULL,0,NULL} ,
  {425, "Seq-graph" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[8]} ,
  {0, "ids" ,128,3,0,0,0,0,0,0,NULL,&atx[60],&atx[214],0,&atx[215]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {0, "locs" ,128,4,0,0,0,0,0,0,NULL,&atx[60],&atx[216],0,&atx[217]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} ,
  {0, "seq-table" ,128,5,0,0,0,0,0,0,NULL,&atx[218],NULL,0,NULL} ,
  {437, "Seq-table" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[76]} };

static AsnModule ampx[1] = {
  { "NCBI-Sequence" , "asnseq.h20",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Sequence
*
**************************************************/

#define ANNOTDESC &at[0]
#define ANNOTDESC_name &at[1]
#define ANNOTDESC_title &at[3]
#define ANNOTDESC_comment &at[4]
#define ANNOTDESC_pub &at[5]
#define ANNOTDESC_user &at[48]
#define ANNOTDESC_create_date &at[50]
#define ANNOTDESC_update_date &at[52]
#define ANNOTDESC_src &at[53]
#define ANNOTDESC_align &at[55]
#define ANNOTDESC_region &at[61]

#define ANNOT_DESCR &at[63]
#define ANNOT_DESCR_E &at[64]

#define BIOSEQ &at[65]
#define BIOSEQ_id &at[66]
#define BIOSEQ_id_E &at[67]
#define BIOSEQ_descr &at[68]
#define BIOSEQ_inst &at[116]
#define BIOSEQ_annot &at[186]
#define BIOSEQ_annot_E &at[187]

#define GIBB_MOL &at[73]

#define HETEROGEN &at[106]

#define MOLINFO &at[110]
#define MOLINFO_biomol &at[111]
#define MOLINFO_tech &at[112]
#define MOLINFO_techexp &at[113]
#define MOLINFO_completeness &at[114]
#define MOLINFO_gbmoltype &at[115]

#define NUMBERING &at[12]
#define NUMBERING_cont &at[13]
#define NUMBERING_enum &at[21]
#define NUMBERING_ref &at[27]
#define NUMBERING_real &at[33]

#define PUBDESC &at[6]
#define PUBDESC_pub &at[7]
#define PUBDESC_name &at[9]
#define PUBDESC_fig &at[10]
#define PUBDESC_num &at[11]
#define PUBDESC_numexc &at[40]
#define PUBDESC_poly_a &at[41]
#define PUBDESC_maploc &at[42]
#define PUBDESC_seq_raw &at[43]
#define PUBDESC_align_group &at[45]
#define PUBDESC_comment &at[46]
#define PUBDESC_reftype &at[47]

#define SEQ_ANNOT &at[188]
#define SEQ_ANNOT_id &at[189]
#define SEQ_ANNOT_id_E &at[190]
#define SEQ_ANNOT_db &at[202]
#define SEQ_ANNOT_name &at[203]
#define SEQ_ANNOT_desc &at[204]
#define SEQ_ANNOT_data &at[205]
#define SEQ_ANNOT_data_ftable &at[206]
#define SEQ_ANNOT_data_ftable_E &at[207]
#define SEQ_ANNOT_data_align &at[208]
#define SEQ_ANNOT_data_align_E &at[209]
#define SEQ_ANNOT_data_graph &at[210]
#define SEQ_ANNOT_data_graph_E &at[211]
#define SEQ_ANNOT_data_ids &at[213]
#define SEQ_ANNOT_data_ids_E &at[214]
#define SEQ_ANNOT_data_locs &at[215]
#define SEQ_ANNOT_data_locs_E &at[216]
#define SEQ_ANNOT_data_seq_table &at[217]

#define SEQ_DATA &at[126]
#define SEQ_DATA_iupacna &at[127]
#define SEQ_DATA_iupacaa &at[129]
#define SEQ_DATA_ncbi2na &at[131]
#define SEQ_DATA_ncbi4na &at[134]
#define SEQ_DATA_ncbi8na &at[136]
#define SEQ_DATA_ncbipna &at[138]
#define SEQ_DATA_ncbi8aa &at[140]
#define SEQ_DATA_ncbieaa &at[142]
#define SEQ_DATA_ncbipaa &at[144]
#define SEQ_DATA_ncbistdaa &at[146]
#define SEQ_DATA_gap &at[148]

#define SEQDESC &at[71]
#define SEQDESC_mol_type &at[72]
#define SEQDESC_modif &at[74]
#define SEQDESC_modif_E &at[75]
#define SEQDESC_method &at[77]
#define SEQDESC_name &at[79]
#define SEQDESC_title &at[80]
#define SEQDESC_org &at[81]
#define SEQDESC_comment &at[83]
#define SEQDESC_num &at[84]
#define SEQDESC_maploc &at[85]
#define SEQDESC_pir &at[87]
#define SEQDESC_genbank &at[89]
#define SEQDESC_pub &at[91]
#define SEQDESC_region &at[92]
#define SEQDESC_user &at[93]
#define SEQDESC_sp &at[94]
#define SEQDESC_dbxref &at[96]
#define SEQDESC_embl &at[97]
#define SEQDESC_create_date &at[99]
#define SEQDESC_update_date &at[100]
#define SEQDESC_prf &at[101]
#define SEQDESC_pdb &at[103]
#define SEQDESC_het &at[105]
#define SEQDESC_source &at[107]
#define SEQDESC_molinfo &at[109]

#define SEQ_DESCR &at[69]
#define SEQ_DESCR_E &at[70]

#define SEQ_EXT &at[153]
#define SEQ_EXT_seg &at[154]
#define SEQ_EXT_ref &at[157]
#define SEQ_EXT_map &at[159]
#define SEQ_EXT_delta &at[163]

#define SEQ_HIST &at[174]
#define SEQ_HIST_assembly &at[175]
#define SEQ_HIST_assembly_E &at[176]
#define SEQ_HIST_replaces &at[177]
#define SEQ_HIST_replaced_by &at[182]
#define SEQ_HIST_deleted &at[183]
#define SEQ_HIST_deleted_bool &at[184]
#define SEQ_HIST_deleted_date &at[185]

#define SEQ_INST &at[117]
#define SEQ_INST_repr &at[118]
#define SEQ_INST_mol &at[119]
#define SEQ_INST_length &at[120]
#define SEQ_INST_fuzz &at[121]
#define SEQ_INST_topology &at[123]
#define SEQ_INST_strand &at[124]
#define SEQ_INST_seq_data &at[125]
#define SEQ_INST_ext &at[152]
#define SEQ_INST_hist &at[173]

#define SEQ_LITERAL &at[169]
#define SEQ_LITERAL_length &at[170]
#define SEQ_LITERAL_fuzz &at[171]
#define SEQ_LITERAL_seq_data &at[172]

#define DELTA_EXT &at[164]
#define DELTA_EXT_E &at[165]

#define GIBB_MOD &at[76]

#define GIBB_METHOD &at[78]

#define NUM_CONT &at[14]
#define NUM_CONT_refnum &at[15]
#define NUM_CONT_has_zero &at[17]
#define NUM_CONT_ascending &at[19]

#define NUM_ENUM &at[22]
#define NUM_ENUM_num &at[23]
#define NUM_ENUM_names &at[24]
#define NUM_ENUM_names_E &at[25]

#define NUM_REF &at[28]
#define NUM_REF_type &at[29]
#define NUM_REF_aligns &at[31]

#define NUM_REAL &at[34]
#define NUM_REAL_a &at[35]
#define NUM_REAL_b &at[37]
#define NUM_REAL_units &at[38]

#define SEG_EXT &at[155]
#define SEG_EXT_E &at[156]

#define REF_EXT &at[158]

#define MAP_EXT &at[160]
#define MAP_EXT_E &at[161]

#define DELTA_SEQ &at[166]
#define DELTA_SEQ_loc &at[167]
#define DELTA_SEQ_literal &at[168]

#define SEQ_HIST_REC &at[178]
#define SEQ_HIST_REC_date &at[179]
#define SEQ_HIST_REC_ids &at[180]
#define SEQ_HIST_REC_ids_E &at[181]

#define IUPACNA &at[128]

#define IUPACAA &at[130]

#define NCBI2NA &at[132]

#define NCBI4NA &at[135]

#define NCBI8NA &at[137]

#define NCBIPNA &at[139]

#define NCBI8AA &at[141]

#define NCBIEAA &at[143]

#define NCBIPAA &at[145]

#define NCBISTDAA &at[147]

#define SEQ_GAP &at[149]
#define SEQ_GAP_type &at[150]
#define SEQ_GAP_linkage &at[151]

#define TEXTANNOT_ID &at[197]
#define TEXTANNOT_ID_name &at[198]
#define TEXTANNOT_ID_accession &at[199]
#define TEXTANNOT_ID_release &at[200]
#define TEXTANNOT_ID_version &at[201]

#define ANNOT_ID &at[191]
#define ANNOT_ID_local &at[192]
#define ANNOT_ID_ncbi &at[194]
#define ANNOT_ID_general &at[195]
#define ANNOT_ID_other &at[196]

#define ALIGN_DEF &at[56]
#define ALIGN_DEF_align_type &at[57]
#define ALIGN_DEF_ids &at[58]
#define ALIGN_DEF_ids_E &at[59]
