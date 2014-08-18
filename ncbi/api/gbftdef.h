/*****************************************************************************
*   gbfeatdfn.h:
*   -- GenBank Feature table define file
*
******************************************************************************/
#ifndef _GBFEATDFN_
#define _GBFEATDFN_


#define GBQUAL_allele             0
#define GBQUAL_anticodon          1
#define GBQUAL_bound_moiety       2
#define GBQUAL_cell_line          3
#define GBQUAL_cell_type          4
#define GBQUAL_chromosome         5
#define GBQUAL_chloroplast        6
#define GBQUAL_chromoplast        7
#define GBQUAL_citation           8
#define GBQUAL_clone              9
#define GBQUAL_clone_lib         10
#define GBQUAL_codon             11
#define GBQUAL_codon_start       12
#define GBQUAL_cons_splice       13
#define GBQUAL_cultivar          14
#define GBQUAL_cyanelle          15
#define GBQUAL_db_xref           16
#define GBQUAL_dev_stage         17
#define GBQUAL_direction         18
#define GBQUAL_EC_number         19
#define GBQUAL_evidence          20
#define GBQUAL_exception         21
#define GBQUAL_frequency         22
#define GBQUAL_function          23
#define GBQUAL_gene              24
#define GBQUAL_gdb_xref          25
#define GBQUAL_germline          26
#define GBQUAL_haplotype         27
#define GBQUAL_insertion_seq     28
#define GBQUAL_isolate           29
#define GBQUAL_kinetoplast       30
#define GBQUAL_label             31
#define GBQUAL_lab_host          32
#define GBQUAL_map               33
#define GBQUAL_macronuclear      34
#define GBQUAL_mitochondrion     35
#define GBQUAL_mod_base          36
#define GBQUAL_note              37
#define GBQUAL_number            38
#define GBQUAL_organism          39
#define GBQUAL_partial           40
#define GBQUAL_PCR_conditions    41
#define GBQUAL_pop_variant       42
#define GBQUAL_phenotype         43
#define GBQUAL_plasmid           44
#define GBQUAL_product           45
#define GBQUAL_proviral          46
#define GBQUAL_pseudo            47          
#define GBQUAL_rearranged        48
#define GBQUAL_replace           49
#define GBQUAL_rpt_family        50
#define GBQUAL_rpt_type          51
#define GBQUAL_rpt_unit          52
#define GBQUAL_sex               53
#define GBQUAL_sequenced_mol     54
#define GBQUAL_serotype		     55
#define GBQUAL_specific_host     56
#define GBQUAL_standard_name     57
#define GBQUAL_strain            58
#define GBQUAL_sub_clone         59
#define GBQUAL_sub_species       60
#define GBQUAL_sub_strain        61
#define GBQUAL_tissue_lib        62
#define GBQUAL_tissue_type       63
#define GBQUAL_translation       64
#define GBQUAL_transl_except     65
#define GBQUAL_transl_table      66
#define GBQUAL_transposon        67
#define GBQUAL_usedin            68
#define GBQUAL_variety           69
#define GBQUAL_virion            70
#define GBQUAL_focus             71
#define GBQUAL_specimen_voucher  72
#define GBQUAL_protein_id        73
#define GBQUAL_country           74
#define GBQUAL_organelle         75
#define GBQUAL_transcript_id     76
#define GBQUAL_transgenic        77
#define GBQUAL_environmental_sample 78
#define GBQUAL_isolation_source  79
#define GBQUAL_serovar           80
#define GBQUAL_locus_tag         81
#define GBQUAL_mol_type          82
#define GBQUAL_segment           83
#define GBQUAL_ecotype           84
#define GBQUAL_estimated_length  85
#define GBQUAL_operon            86
#define GBQUAL_old_locus_tag     87
#define GBQUAL_compare           88
#define GBQUAL_experiment        89
#define GBQUAL_inference         90
#define GBQUAL_rpt_unit_seq      91
#define GBQUAL_rpt_unit_range    92
#define GBQUAL_ribosomal_slippage 93
#define GBQUAL_trans_splicing    94
#define GBQUAL_collected_by      95
#define GBQUAL_collection_date   96
#define GBQUAL_identified_by     97
#define GBQUAL_lat_lon           98
#define GBQUAL_PCR_primers       99
#define GBQUAL_mobile_element   100
#define GBQUAL_metagenomic      101
#define GBQUAL_culture_collection 102
#define GBQUAL_bio_material     103
#define GBQUAL_ncRNA_class      104
#define GBQUAL_tag_peptide      105
#define GBQUAL_mating_type      106
#define GBQUAL_satellite        107
#define GBQUAL_gene_synonym     108

#define ParFlat_TOTAL_GBQUAL    109
#define ParFlat_TOTAL_IntOr       3
#define ParFlat_TOTAL_LRB         3
#define ParFlat_TOTAL_Exp         2
#define ParFlat_TOTAL_Rpt         7
#define ParFlat_TOTAL_GBFEAT     69

#define  Class_pos_aa             1
#define  Class_text               2
#define  Class_bracket_int        3
#define  Class_seq_aa             4
#define  Class_int_or             5
#define  Class_site               6
#define  Class_L_R_B              7
#define  Class_ecnum              8
#define  Class_exper              9
#define  Class_none              10
#define  Class_token             11
#define  Class_int               12
#define  Class_rpt               13
#define  Class_flabel_base       14
#define  Class_flabel_dbname     15
#define  Class_note              16
#define  Class_number            17


#define  ParFlat_Stoken_type            1
#define  ParFlat_BracketInt_type        2
#define  ParFlat_Integer_type           3      
#define  ParFlat_Number_type            4      

/*********************************************************************/

typedef struct sematic_gbfeature {
    CharPtr  key;
    Int2     mand_num;
    Int2     mand_qual[5];
    Int2     opt_num;
    Int2     opt_qual[65];
} SematicFeat, PNTR SematicFeatPtr;

typedef struct gbfeat_name {
    CharPtr  name;
    Uint1     gbclass;
} GbFeatName, PNTR GbFeatNamePtr;

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

NLM_EXTERN GbFeatNamePtr x_ParFlat_GBQual_names PROTO((void));
#define ParFlat_GBQual_names x_ParFlat_GBQual_names()

NLM_EXTERN SematicFeatPtr x_ParFlat_GBFeat PROTO((void));
#define ParFlat_GBFeat x_ParFlat_GBFeat()

NLM_EXTERN ValNodePtr Validate_ParFlat_GBFeat (void);

extern CharPtr     ParFlat_IntOrString [ParFlat_TOTAL_IntOr];
extern CharPtr     ParFlat_LRBString   [ParFlat_TOTAL_LRB];
extern CharPtr     ParFlat_ExpString   [ParFlat_TOTAL_Exp];
extern CharPtr     ParFlat_RptString   [ParFlat_TOTAL_Rpt];

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
