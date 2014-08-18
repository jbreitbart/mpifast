/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asninsdseq.h17";
static AsnType atx[106] = {
  {401, "INSDSet" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[1],0,&atx[2]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {402, "INSDSeq" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[3],0,&atx[21]} ,
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
  {403, "INSDSeqid" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[25]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "secondary-accessions" ,128,15,0,1,0,0,0,0,NULL,&atx[22],&atx[24],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {404, "INSDSecondary-accn" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[29]} ,
  {0, "project" ,128,16,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[27]} ,
  {0, "keywords" ,128,17,0,1,0,0,0,0,NULL,&atx[22],&atx[28],0,&atx[30]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[29],NULL,0,NULL} ,
  {405, "INSDKeyword" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[36]} ,
  {0, "segment" ,128,18,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[31]} ,
  {0, "source" ,128,19,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[32]} ,
  {0, "organism" ,128,20,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[33]} ,
  {0, "taxonomy" ,128,21,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[34]} ,
  {0, "references" ,128,22,0,1,0,0,0,0,NULL,&atx[22],&atx[35],0,&atx[54]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[36],NULL,0,NULL} ,
  {406, "INSDReference" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[37],0,&atx[56]} ,
  {0, "reference" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[38]} ,
  {0, "position" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[39]} ,
  {0, "authors" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[40],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[41],NULL,0,NULL} ,
  {409, "INSDAuthor" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[47]} ,
  {0, "consortium" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[43]} ,
  {0, "title" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[44]} ,
  {0, "journal" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[45]} ,
  {0, "xref" ,128,6,0,1,0,0,0,0,NULL,&atx[51],&atx[46],0,&atx[52]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {410, "INSDXref" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[48],0,&atx[61]} ,
  {0, "dbname" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[49]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "pubmed" ,128,7,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[53]} ,
  {0, "remark" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "comment" ,128,23,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[55]} ,
  {0, "tagset" ,128,24,0,1,0,0,0,0,NULL,&atx[56],NULL,0,&atx[67]} ,
  {407, "INSDTagset" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[57],0,&atx[72]} ,
  {0, "authority" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[58]} ,
  {0, "version" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[59]} ,
  {0, "url" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[60]} ,
  {0, "tags" ,128,3,0,1,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {411, "INSDTags" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[62],0,&atx[63]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[63],NULL,0,NULL} ,
  {412, "INSDTag" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[64],0,&atx[77]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[65]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[66]} ,
  {0, "unit" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "primary" ,128,25,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[68]} ,
  {0, "source-db" ,128,26,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[69]} ,
  {0, "database-reference" ,128,27,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[70]} ,
  {0, "feature-table" ,128,28,0,1,0,0,0,0,NULL,&atx[22],&atx[71],0,&atx[93]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {408, "INSDFeature" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[73],0,&atx[41]} ,
  {0, "key" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[74]} ,
  {0, "location" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[75]} ,
  {0, "intervals" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[76],0,&atx[85]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[77],NULL,0,NULL} ,
  {413, "INSDInterval" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[78],0,&atx[90]} ,
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
  {414, "INSDQualifier" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[91],0,&atx[95]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[92]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "sequence" ,128,29,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[94]} ,
  {0, "contig" ,128,30,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {415, "INSDTagsetRules" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[96],0,&atx[99]} ,
  {0, "authority" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[97]} ,
  {0, "version" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[98]} ,
  {0, "mandatorytags" ,128,2,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[101]} ,
  {416, "INSDTagNames" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[100],0,&atx[104]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "optionaltags" ,128,3,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[102]} ,
  {0, "uniquetags" ,128,4,0,1,0,0,0,0,NULL,&atx[99],NULL,0,&atx[103]} ,
  {0, "extensible" ,128,5,0,1,0,0,0,0,NULL,&atx[82],NULL,0,NULL} ,
  {417, "INSDTagsetRuleSet" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[105],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[95],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "INSD-INSDSeq" , "asninsdseq.h17",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module INSD-INSDSeq
*
**************************************************/

#define INSDSET &at[0]
#define INSDSET_E &at[1]

#define INSDSEQ &at[2]
#define INSDSEQ_locus &at[3]
#define INSDSEQ_length &at[5]
#define INSDSEQ_strandedness &at[7]
#define INSDSEQ_moltype &at[8]
#define INSDSEQ_topology &at[9]
#define INSDSEQ_division &at[10]
#define INSDSEQ_update_date &at[11]
#define INSDSEQ_create_date &at[12]
#define INSDSEQ_update_release &at[13]
#define INSDSEQ_create_release &at[14]
#define INSDSEQ_definition &at[15]
#define INSDSEQ_primary_accession &at[16]
#define INSDSEQ_entry_version &at[17]
#define INSDSEQ_accession_version &at[18]
#define INSDSEQ_other_seqids &at[19]
#define INSDSEQ_other_seqids_E &at[20]
#define INSDSEQ_secondary_accessions &at[23]
#define INSDSEQ_secondary_accessions_E &at[24]
#define INSDSEQ_project &at[26]
#define INSDSEQ_keywords &at[27]
#define INSDSEQ_keywords_E &at[28]
#define INSDSEQ_segment &at[30]
#define INSDSEQ_source &at[31]
#define INSDSEQ_organism &at[32]
#define INSDSEQ_taxonomy &at[33]
#define INSDSEQ_references &at[34]
#define INSDSEQ_references_E &at[35]
#define INSDSEQ_comment &at[54]
#define INSDSEQ_tagset &at[55]
#define INSDSEQ_primary &at[67]
#define INSDSEQ_source_db &at[68]
#define INSDSEQ_database_reference &at[69]
#define INSDSEQ_feature_table &at[70]
#define INSDSEQ_feature_table_E &at[71]
#define INSDSEQ_sequence &at[93]
#define INSDSEQ_contig &at[94]

#define INSDSEQID &at[21]

#define INSDSECONDARY_ACCN &at[25]

#define INSDKEYWORD &at[29]

#define INSDREFERENCE &at[36]
#define INSDREFERENCE_reference &at[37]
#define INSDREFERENCE_position &at[38]
#define INSDREFERENCE_authors &at[39]
#define INSDREFERENCE_authors_E &at[40]
#define INSDREFERENCE_consortium &at[42]
#define INSDREFERENCE_title &at[43]
#define INSDREFERENCE_journal &at[44]
#define INSDREFERENCE_xref &at[45]
#define INSDREFERENCE_xref_E &at[46]
#define INSDREFERENCE_pubmed &at[52]
#define INSDREFERENCE_remark &at[53]

#define INSDTAGSET &at[56]
#define INSDTAGSET_authority &at[57]
#define INSDTAGSET_version &at[58]
#define INSDTAGSET_url &at[59]
#define INSDTAGSET_tags &at[60]

#define INSDFEATURE &at[72]
#define INSDFEATURE_key &at[73]
#define INSDFEATURE_location &at[74]
#define INSDFEATURE_intervals &at[75]
#define INSDFEATURE_intervals_E &at[76]
#define INSDFEATURE_operator &at[85]
#define INSDFEATURE_partial5 &at[86]
#define INSDFEATURE_partial3 &at[87]
#define INSDFEATURE_quals &at[88]
#define INSDFEATURE_quals_E &at[89]

#define INSDAUTHOR &at[41]

#define INSDXREF &at[47]
#define INSDXREF_dbname &at[48]
#define INSDXREF_id &at[49]

#define INSDTAGS &at[61]
#define INSDTAGS_E &at[62]

#define INSDTAG &at[63]
#define INSDTAG_name &at[64]
#define INSDTAG_value &at[65]
#define INSDTAG_unit &at[66]

#define INSDINTERVAL &at[77]
#define INSDINTERVAL_from &at[78]
#define INSDINTERVAL_to &at[79]
#define INSDINTERVAL_point &at[80]
#define INSDINTERVAL_iscomp &at[81]
#define INSDINTERVAL_interbp &at[83]
#define INSDINTERVAL_accession &at[84]

#define INSDQUALIFIER &at[90]
#define INSDQUALIFIER_name &at[91]
#define INSDQUALIFIER_value &at[92]

#define INSDTAGSETRULES &at[95]
#define INSDTAGSETRULES_authority &at[96]
#define INSDTAGSETRULES_version &at[97]
#define INSDTAGSETRULES_mandatorytags &at[98]
#define INSDTAGSETRULES_optionaltags &at[101]
#define INSDTAGSETRULES_uniquetags &at[102]
#define INSDTAGSETRULES_extensible &at[103]

#define INSDTAGNAMES &at[99]
#define INSDTAGNAMES_E &at[100]

#define INSDTAGSETRULESET &at[104]
#define INSDTAGSETRULESET_E &at[105]
