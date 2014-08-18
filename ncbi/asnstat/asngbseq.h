/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asngbseq.h67";
static AsnType atx[106] = {
  {401, "GBSet" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[1],0,&atx[2]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {402, "GBSeq" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[3],0,&atx[21]} ,
  {0, "locus" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "length" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "strandedness" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[8]} ,
  {0, "moltype" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[9]} ,
  {0, "topology" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[10]} ,
  {0, "division" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[11]} ,
  {0, "update-date" ,128,6,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[12]} ,
  {0, "create-date" ,128,7,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[13]} ,
  {0, "update-release" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[14]} ,
  {0, "create-release" ,128,9,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[15]} ,
  {0, "definition" ,128,10,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[16]} ,
  {0, "primary-accession" ,128,11,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[17]} ,
  {0, "entry-version" ,128,12,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[18]} ,
  {0, "accession-version" ,128,13,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[19]} ,
  {0, "other-seqids" ,128,14,0,1,0,0,0,0,NULL,&atx[22],&atx[20],0,&atx[23]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {403, "GBSeqid" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[25]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "secondary-accessions" ,128,15,0,1,0,0,0,0,NULL,&atx[22],&atx[24],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {404, "GBSecondary-accn" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[29]} ,
  {0, "project" ,128,16,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[27]} ,
  {0, "keywords" ,128,17,0,1,0,0,0,0,NULL,&atx[22],&atx[28],0,&atx[30]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[29],NULL,0,NULL} ,
  {405, "GBKeyword" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[36]} ,
  {0, "segment" ,128,18,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[31]} ,
  {0, "source" ,128,19,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[32]} ,
  {0, "organism" ,128,20,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[33]} ,
  {0, "taxonomy" ,128,21,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[34]} ,
  {0, "references" ,128,22,0,1,0,0,0,0,NULL,&atx[22],&atx[35],0,&atx[54]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[36],NULL,0,NULL} ,
  {406, "GBReference" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[37],0,&atx[56]} ,
  {0, "reference" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[38]} ,
  {0, "position" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[39]} ,
  {0, "authors" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[40],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[41],NULL,0,NULL} ,
  {409, "GBAuthor" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[47]} ,
  {0, "consortium" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[43]} ,
  {0, "title" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[44]} ,
  {0, "journal" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[45]} ,
  {0, "xref" ,128,6,0,1,0,0,0,0,NULL,&atx[51],&atx[46],0,&atx[52]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {410, "GBXref" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[48],0,&atx[61]} ,
  {0, "dbname" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[49]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "pubmed" ,128,7,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[53]} ,
  {0, "remark" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "comment" ,128,23,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[55]} ,
  {0, "tagset" ,128,24,0,1,0,0,0,0,NULL,&atx[56],NULL,0,&atx[67]} ,
  {407, "GBTagset" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[57],0,&atx[72]} ,
  {0, "authority" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[58]} ,
  {0, "version" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[59]} ,
  {0, "url" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[60]} ,
  {0, "tags" ,128,3,0,1,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {411, "GBTags" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[62],0,&atx[63]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[63],NULL,0,NULL} ,
  {412, "GBTag" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[64],0,&atx[77]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[65]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[66]} ,
  {0, "unit" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "primary" ,128,25,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[68]} ,
  {0, "source-db" ,128,26,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[69]} ,
  {0, "database-reference" ,128,27,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[70]} ,
  {0, "feature-table" ,128,28,0,1,0,0,0,0,NULL,&atx[22],&atx[71],0,&atx[93]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {408, "GBFeature" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[73],0,&atx[41]} ,
  {0, "key" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[74]} ,
  {0, "location" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[75]} ,
  {0, "intervals" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[76],0,&atx[85]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[77],NULL,0,NULL} ,
  {413, "GBInterval" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[78],0,&atx[90]} ,
  {0, "from" ,128,0,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[79]} ,
  {0, "to" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[80]} ,
  {0, "point" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[81]} ,
  {0, "iscomp" ,128,3,0,1,0,0,0,0,NULL,&atx[82],NULL,0,&atx[83]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "interbp" ,128,4,0,1,0,0,0,0,NULL,&atx[82],NULL,0,&atx[84]} ,
  {0, "accession" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "operator" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[86]} ,
  {0, "partial5" ,128,4,0,1,0,0,0,0,NULL,&atx[82],NULL,0,&atx[87]} ,
  {0, "partial3" ,128,5,0,1,0,0,0,0,NULL,&atx[82],NULL,0,&atx[88]} ,
  {0, "quals" ,128,6,0,1,0,0,0,0,NULL,&atx[22],&atx[89],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[90],NULL,0,NULL} ,
  {414, "GBQualifier" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[91],0,&atx[95]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[92]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "sequence" ,128,29,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[94]} ,
  {0, "contig" ,128,30,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {415, "GBTagsetRules" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[96],0,&atx[99]} ,
  {0, "authority" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[97]} ,
  {0, "version" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[98]} ,
  {0, "mandatorytags" ,128,2,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[101]} ,
  {416, "GBTagNames" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[100],0,&atx[104]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "optionaltags" ,128,3,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[102]} ,
  {0, "uniquetags" ,128,4,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[103]} ,
  {0, "extensible" ,128,5,0,1,0,0,0,0,NULL,&atx[82],NULL,0,NULL} ,
  {417, "GBTagsetRuleSet" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[105],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[95],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-GBSeq" , "asngbseq.h67",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-GBSeq
*
**************************************************/

#define GBSET &at[0]
#define GBSET_E &at[1]

#define GBSEQ &at[2]
#define GBSEQ_locus &at[3]
#define GBSEQ_length &at[5]
#define GBSEQ_strandedness &at[7]
#define GBSEQ_moltype &at[8]
#define GBSEQ_topology &at[9]
#define GBSEQ_division &at[10]
#define GBSEQ_update_date &at[11]
#define GBSEQ_create_date &at[12]
#define GBSEQ_update_release &at[13]
#define GBSEQ_create_release &at[14]
#define GBSEQ_definition &at[15]
#define GBSEQ_primary_accession &at[16]
#define GBSEQ_entry_version &at[17]
#define GBSEQ_accession_version &at[18]
#define GBSEQ_other_seqids &at[19]
#define GBSEQ_other_seqids_E &at[20]
#define GBSEQ_secondary_accessions &at[23]
#define GBSEQ_secondary_accessions_E &at[24]
#define GBSEQ_project &at[26]
#define GBSEQ_keywords &at[27]
#define GBSEQ_keywords_E &at[28]
#define GBSEQ_segment &at[30]
#define GBSEQ_source &at[31]
#define GBSEQ_organism &at[32]
#define GBSEQ_taxonomy &at[33]
#define GBSEQ_references &at[34]
#define GBSEQ_references_E &at[35]
#define GBSEQ_comment &at[54]
#define GBSEQ_tagset &at[55]
#define GBSEQ_primary &at[67]
#define GBSEQ_source_db &at[68]
#define GBSEQ_database_reference &at[69]
#define GBSEQ_feature_table &at[70]
#define GBSEQ_feature_table_E &at[71]
#define GBSEQ_sequence &at[93]
#define GBSEQ_contig &at[94]

#define GBSEQID &at[21]

#define GBSECONDARY_ACCN &at[25]

#define GBKEYWORD &at[29]

#define GBREFERENCE &at[36]
#define GBREFERENCE_reference &at[37]
#define GBREFERENCE_position &at[38]
#define GBREFERENCE_authors &at[39]
#define GBREFERENCE_authors_E &at[40]
#define GBREFERENCE_consortium &at[42]
#define GBREFERENCE_title &at[43]
#define GBREFERENCE_journal &at[44]
#define GBREFERENCE_xref &at[45]
#define GBREFERENCE_xref_E &at[46]
#define GBREFERENCE_pubmed &at[52]
#define GBREFERENCE_remark &at[53]

#define GBTAGSET &at[56]
#define GBTAGSET_authority &at[57]
#define GBTAGSET_version &at[58]
#define GBTAGSET_url &at[59]
#define GBTAGSET_tags &at[60]

#define GBFEATURE &at[72]
#define GBFEATURE_key &at[73]
#define GBFEATURE_location &at[74]
#define GBFEATURE_intervals &at[75]
#define GBFEATURE_intervals_E &at[76]
#define GBFEATURE_operator &at[85]
#define GBFEATURE_partial5 &at[86]
#define GBFEATURE_partial3 &at[87]
#define GBFEATURE_quals &at[88]
#define GBFEATURE_quals_E &at[89]

#define GBAUTHOR &at[41]

#define GBXREF &at[47]
#define GBXREF_dbname &at[48]
#define GBXREF_id &at[49]

#define GBTAGS &at[61]
#define GBTAGS_E &at[62]

#define GBTAG &at[63]
#define GBTAG_name &at[64]
#define GBTAG_value &at[65]
#define GBTAG_unit &at[66]

#define GBINTERVAL &at[77]
#define GBINTERVAL_from &at[78]
#define GBINTERVAL_to &at[79]
#define GBINTERVAL_point &at[80]
#define GBINTERVAL_iscomp &at[81]
#define GBINTERVAL_interbp &at[83]
#define GBINTERVAL_accession &at[84]

#define GBQUALIFIER &at[90]
#define GBQUALIFIER_name &at[91]
#define GBQUALIFIER_value &at[92]

#define GBTAGSETRULES &at[95]
#define GBTAGSETRULES_authority &at[96]
#define GBTAGSETRULES_version &at[97]
#define GBTAGSETRULES_mandatorytags &at[98]
#define GBTAGSETRULES_optionaltags &at[101]
#define GBTAGSETRULES_uniquetags &at[102]
#define GBTAGSETRULES_extensible &at[103]

#define GBTAGNAMES &at[99]
#define GBTAGNAMES_E &at[100]

#define GBTAGSETRULESET &at[104]
#define GBTAGSETRULESET_E &at[105]
