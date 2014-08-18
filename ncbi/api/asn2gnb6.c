/*   asn2gnb6.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  asn2gnb6.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* Version Creation Date:   10/21/98
*
* $Revision: 1.187 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objpubme.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbi.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

static CharPtr link_tax = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?";

static CharPtr link_featn = "http://www.ncbi.nlm.nih.gov/nuccore/";
static CharPtr link_featp = "http://www.ncbi.nlm.nih.gov/protein/";

static CharPtr link_seqn = "http://www.ncbi.nlm.nih.gov/nuccore/";
static CharPtr link_seqp = "http://www.ncbi.nlm.nih.gov/protein/";

static CharPtr link_lat_lon = "http://www.ncbi.nlm.nih.gov/projects/Sequin/latlonview.html?";




/* ordering arrays for qualifiers and note components */

static SourceType source_qual_order [] = {
  SCQUAL_organism,

  SCQUAL_organelle,

  SCQUAL_mol_type,

  SCQUAL_strain,
  SCQUAL_sub_strain,
  SCQUAL_variety,
  SCQUAL_serotype,
  SCQUAL_serovar,
  SCQUAL_cultivar,
  SCQUAL_isolate,
  SCQUAL_isolation_source,
  SCQUAL_spec_or_nat_host,
  SCQUAL_sub_species,

  SCQUAL_specimen_voucher,
  SCQUAL_culture_collection,
  SCQUAL_bio_material,

  SCQUAL_db_xref,
  SCQUAL_org_xref,

  SCQUAL_chromosome,

  SCQUAL_segment,

  SCQUAL_map,
  SCQUAL_clone,
  SCQUAL_sub_clone,
  SCQUAL_haplotype,
  SCQUAL_sex,
  SCQUAL_mating_type,
  SCQUAL_cell_line,
  SCQUAL_cell_type,
  SCQUAL_tissue_type,
  SCQUAL_clone_lib,
  SCQUAL_dev_stage,
  SCQUAL_ecotype,
  SCQUAL_frequency,

  SCQUAL_germline,
  SCQUAL_rearranged,
  SCQUAL_transgenic,
  SCQUAL_environmental_sample,

  SCQUAL_lab_host,
  SCQUAL_pop_variant,
  SCQUAL_tissue_lib,

  SCQUAL_plasmid_name,
  SCQUAL_transposon_name,
  SCQUAL_ins_seq_name,

  SCQUAL_country,

  SCQUAL_focus,

  SCQUAL_lat_lon,
  SCQUAL_collection_date,
  SCQUAL_collected_by,
  SCQUAL_identified_by,
  /*
  SCQUAL_fwd_primer_seq,
  SCQUAL_rev_primer_seq,
  SCQUAL_fwd_primer_name,
  SCQUAL_rev_primer_name,
  */
  SCQUAL_PCR_primers,
  SCQUAL_PCR_reaction,

  SCQUAL_note,

  SCQUAL_sequenced_mol,
  SCQUAL_label,
  SCQUAL_usedin,
  SCQUAL_citation,
  (SourceType) 0
};

static SourceType source_desc_note_order [] = {
  SCQUAL_seqfeat_note,
  SCQUAL_orgmod_note,
  SCQUAL_subsource_note,

  SCQUAL_metagenomic,

  SCQUAL_linkage_group,

  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_serogroup,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,

  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,

  SCQUAL_haplogroup,

  SCQUAL_metagenome_source,
  SCQUAL_metagenome_note,

  SCQUAL_genotype,
  SCQUAL_plastid_name,

  SCQUAL_endogenous_virus_name,

  SCQUAL_common_name,

  SCQUAL_PCR_primer_note,
  SCQUAL_PCR_reaction,

  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_zero_subsrc,

  /* SCQUAL_old_lineage, */

  /* SCQUAL_old_name, */
  (SourceType) 0
};

static SourceType source_feat_note_order [] = {
  SCQUAL_unstructured,

  SCQUAL_metagenomic,

  SCQUAL_linkage_group,
  SCQUAL_mating_type,

  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_serogroup,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,

  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,
  
  SCQUAL_haplogroup,
  
  SCQUAL_metagenome_source,
  SCQUAL_metagenome_note,

  SCQUAL_genotype,
  SCQUAL_plastid_name,

  SCQUAL_endogenous_virus_name,

  SCQUAL_seqfeat_note,
  SCQUAL_orgmod_note,
  SCQUAL_subsource_note,

  SCQUAL_common_name,

  SCQUAL_PCR_primer_note,
  SCQUAL_PCR_reaction,

  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_zero_subsrc,

  /* SCQUAL_old_lineage, */

  /* SCQUAL_old_name, */
  (SourceType) 0
};

NLM_EXTERN SourceQual asn2gnbk_source_quals [ASN2GNBK_TOTAL_SOURCE] = {
  { "",                         Qual_class_ignore     },
  { "acronym",                  Qual_class_orgmod     },
  { "anamorph",                 Qual_class_orgmod     },
  { "authority",                Qual_class_orgmod     },
  { "biotype",                  Qual_class_orgmod     },
  { "biovar",                   Qual_class_orgmod     },
  { "bio_material",             Qual_class_voucher    },
  { "breed",                    Qual_class_orgmod     },
  { "cell_line",                Qual_class_subsource  },
  { "cell_type",                Qual_class_subsource  },
  { "chemovar",                 Qual_class_orgmod     },
  { "chromosome",               Qual_class_subsource  },
  { "citation",                 Qual_class_pubset     },
  { "clone",                    Qual_class_subsource  },
  { "clone_lib",                Qual_class_subsource  },
  { "collected_by",             Qual_class_subsource  },
  { "collection_date",          Qual_class_subsource  },
  { "common",                   Qual_class_orgmod     },
  { "common",                   Qual_class_string     },
  { "country",                  Qual_class_subsource  },
  { "cultivar",                 Qual_class_orgmod     },
  { "culture_collection",       Qual_class_voucher    },
  { "db_xref",                  Qual_class_db_xref    },
  { "db_xref",                  Qual_class_db_xref    },
  { "dev_stage",                Qual_class_subsource  },
  { "dosage",                   Qual_class_orgmod     },
  { "ecotype",                  Qual_class_orgmod     },
  { "endogenous_virus",         Qual_class_subsource  },
  { "environmental_sample",     Qual_class_subsource  },
  { "extrachromosomal",         Qual_class_boolean    },
  { "focus",                    Qual_class_boolean    },
  { "forma",                    Qual_class_orgmod     },
  { "forma_specialis",          Qual_class_orgmod     },
  { "frequency",                Qual_class_subsource  },
  { "fwd_primer_name",          Qual_class_subsource  },
  { "fwd_primer_seq",           Qual_class_subsource  },
  { "gb_acronym",               Qual_class_orgmod     },
  { "gb_anamorph",              Qual_class_orgmod     },
  { "gb_synonym",               Qual_class_orgmod     },
  { "genotype",                 Qual_class_subsource  },
  { "germline",                 Qual_class_subsource  },
  { "group",                    Qual_class_orgmod     },
  { "haplogroup",               Qual_class_subsource  },
  { "haplotype",                Qual_class_subsource  },
  { "identified_by",            Qual_class_subsource  },
  { "insertion_seq",            Qual_class_subsource  },
  { "isolate",                  Qual_class_orgmod     },
  { "isolation_source",         Qual_class_subsource  },
  { "lab_host",                 Qual_class_subsource  },
  { "label",                    Qual_class_label      },
  { "lat_lon",                  Qual_class_lat_lon    },
  { "linkage_group",            Qual_class_subsource  },
  { "macronuclear",             Qual_class_boolean    },
  { "map",                      Qual_class_subsource  },
  { "mating_type",              Qual_class_subsource  },
  { "derived from metagenome",  Qual_class_orgmod     },
  { "metagenome_source",        Qual_class_orgmod     },
  { "metagenomic",              Qual_class_subsource  },
  { "mol_type",                 Qual_class_string     },
  { "note",                     Qual_class_note       },
  { "old_lineage",              Qual_class_orgmod     },
  { "old_name",                 Qual_class_orgmod     },
  { "organism",                 Qual_class_string     },
  { "organelle",                Qual_class_organelle  },
  { "orgmod_note",              Qual_class_orgmod     },
  { "pathovar",                 Qual_class_orgmod     },
  { "PCR_primers",              Qual_class_pcr        },
  { "PCR_primers",              Qual_class_pcr        },
  { "PCR_primers",              Qual_class_pcr_react  },
  { "plasmid",                  Qual_class_subsource  },
  { "plastid",                  Qual_class_subsource  },
  { "pop_variant",              Qual_class_subsource  },
  { "rearranged",               Qual_class_subsource  },
  { "rev_primer_name",          Qual_class_subsource  },
  { "rev_primer_seq",           Qual_class_subsource  },
  { "segment",                  Qual_class_subsource  },
  { "seqfeat_note",             Qual_class_string     },
  { "sequenced_mol",            Qual_class_quote      },
  { "serogroup",                Qual_class_orgmod     },
  { "serotype",                 Qual_class_orgmod     },
  { "serovar",                  Qual_class_orgmod     },
  { "sex",                      Qual_class_subsource  },
  { "host",                     Qual_class_orgmod     },
  { "specimen_voucher",         Qual_class_voucher    },
  { "strain",                   Qual_class_orgmod     },
  { "sub_clone",                Qual_class_subsource  },
  { "subgroup",                 Qual_class_orgmod     },
  { "sub_species",              Qual_class_orgmod     },
  { "sub_strain",               Qual_class_orgmod     },
  { "subtype",                  Qual_class_orgmod     },
  { "subsource_note",           Qual_class_subsource  },
  { "synonym",                  Qual_class_orgmod     },
  { "teleomorph",               Qual_class_orgmod     },
  { "tissue_lib",               Qual_class_subsource  },
  { "tissue_type",              Qual_class_subsource  },
  { "transgenic",               Qual_class_subsource  },
  { "transposon",               Qual_class_subsource  },
  { "type",                     Qual_class_orgmod     },
  { "unstructured",             Qual_class_valnode    },
  { "usedin",                   Qual_class_quote      },
  { "variety",                  Qual_class_orgmod     },
  { "?",                        Qual_class_orgmod     },
  { "?",                        Qual_class_orgmod     },
  { "?",                        Qual_class_subsource  }
};

NLM_EXTERN SourceType subSourceToSourceIdx [42] = {
  SCQUAL_zero_subsrc,
  SCQUAL_chromosome,
  SCQUAL_map,
  SCQUAL_clone,
  SCQUAL_sub_clone,
  SCQUAL_haplotype,
  SCQUAL_genotype,
  SCQUAL_sex,
  SCQUAL_cell_line,
  SCQUAL_cell_type,
  SCQUAL_tissue_type,
  SCQUAL_clone_lib,
  SCQUAL_dev_stage,
  SCQUAL_frequency,
  SCQUAL_germline,
  SCQUAL_rearranged,
  SCQUAL_lab_host,
  SCQUAL_pop_variant,
  SCQUAL_tissue_lib,
  SCQUAL_plasmid_name,
  SCQUAL_transposon_name,
  SCQUAL_ins_seq_name,
  SCQUAL_plastid_name,
  SCQUAL_country,
  SCQUAL_segment,
  SCQUAL_endogenous_virus_name,
  SCQUAL_transgenic,
  SCQUAL_environmental_sample,
  SCQUAL_isolation_source,
  SCQUAL_lat_lon,
  SCQUAL_collection_date,
  SCQUAL_collected_by,
  SCQUAL_identified_by,
  SCQUAL_fwd_primer_seq,
  SCQUAL_rev_primer_seq,
  SCQUAL_fwd_primer_name,
  SCQUAL_rev_primer_name,
  SCQUAL_metagenomic,
  SCQUAL_mating_type,
  SCQUAL_linkage_group,
  SCQUAL_haplogroup,
  SCQUAL_subsource_note
};

/* ********************************************************************** */

/* ********************************************************************** */

/* format functions allocate printable string for given paragraph */

/* superset of http://www.ncbi.nlm.nih.gov/collab/db_xref.html and RefSeq db_xrefs */

NLM_EXTERN CharPtr legalDbXrefs [] = {
  "AceView/WormGenes",
  "AFTOL",
  "ApiDB",
  "ApiDB_CryptoDB",
  "ApiDB_PlasmoDB",
  "ApiDB_ToxoDB",
  "ASAP",
  "ATCC",
  "ATCC(in host)",
  "ATCC(dna)",
  "Axeldb",
  "BDGP_EST",
  "BDGP_INS",
  "BEETLEBASE",
  "BioHealthBase",
  "BOLD",
  "CDD",
  "CK",
  "COG",
  "dbClone",
  "dbCloneLib",
  "dbEST",
  "dbProbe",
  "dbSNP",
  "dbSTS",
  "dictyBase",
  "EcoGene",
  "ENSEMBL",
  "ERIC",
  "ESTLIB",
  "FANTOM_DB",
  "FLYBASE",
  "GABI",
  "GDB",
  "GeneDB",
  "GeneID",
  "GO",
  "GOA",
  "GRIN",
  "H-InvDB",
  "HGNC",
  "HSSP",
  "IMGT/GENE-DB",
  "IMGT/HLA",
  "IMGT/LIGM",
  "InterimID",
  "InterPro",
  "ISD",
  "ISFinder",
  "JCM",
  "JGIDB",
  "LocusID",
  "MaizeGDB",
  "MGI",
  "MIM",
  "NBRC",
  "NextDB",
  "niaEST",
  "NMPDR",
  "NRESTdb",
  "Osa1",
  "Pathema",
  "PBmice",
  "PDB",
  "PFAM",
  "PGN",
  "PIR",
  "PSEUDO",
  "PseudoCap",
  "RAP-DB",
  "RATMAP",
  "RFAM",
  "RGD",
  "RiceGenes",
  "RZPD",
  "SEED",
  "SGD",
  "SoyBase",
  "SubtiList",
  "taxon",
  "TIGRFAM",
  "UniGene",
  "UNILIB",
  "UniProtKB/Swiss-Prot",
  "UniProtKB/TrEMBL",
  "UniSTS",
  "UNITE",
  "VBASE2",
  "VectorBase",
  "WorfDB",
  "WormBase",
  "Xenbase",
  "ZFIN",
  NULL
};

NLM_EXTERN CharPtr legalSrcDbXrefs [] = {
  "AFTOL",
  "ATCC",
  "ATCC(dna)",
  "ATCC(in host)",
  "BOLD",
  "FANTOM_DB",
  "FLYBASE",
  "GRIN",
  "IMGT/HLA",
  "IMGT/LIGM",
  "JCM",
  "MGI",
  "NBRC",
  "RZPD",
  "taxon",
  "UNILIB",
  "UNITE",
  NULL
};

NLM_EXTERN CharPtr legalRefSeqDbXrefs [] = {
  "CCDS",
  "CGNC",
  "CloneID",
  "ECOCYC",
  "HPRD",
  "miRBase",
  "PBR",
  "REBASE",
  "TAIR",
  "VBRC",
  NULL
};

static Boolean IsDbxrefInList (
  CharPtr name,
  CharPtr PNTR list,
  size_t num,
  BoolPtr badcapP,
  CharPtr PNTR goodcapP
)

{
  Int2  L, R, mid;

  L = 0;
  R = num;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (list [mid], name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (list [R], name) == 0) {
    if (StringCmp (list [R], name) != 0) {
      if (badcapP != NULL) {
        *badcapP = TRUE;
      }
      if (goodcapP != NULL) {
        *goodcapP = list [R];
      }
    }
    return TRUE;
  }

  return FALSE;
}

NLM_EXTERN Boolean DbxrefIsValid (
  CharPtr name,
  BoolPtr is_refseq_P,
  BoolPtr is_source_P,
  BoolPtr is_badcap_P,
  CharPtr PNTR goodcapP
)

{
  if (is_refseq_P != NULL) {
    *is_refseq_P = FALSE;
  }
  if (is_source_P != NULL) {
    *is_source_P = FALSE;
  }
  if (is_badcap_P != NULL) {
    *is_badcap_P = FALSE;
  }
  if (goodcapP != NULL) {
    *goodcapP = NULL;
  }

  if (StringHasNoText (name)) return FALSE;

  if (IsDbxrefInList (name, legalRefSeqDbXrefs,
                      sizeof (legalRefSeqDbXrefs) / sizeof (legalRefSeqDbXrefs [0]) - 1,
                      is_badcap_P, goodcapP)) {
    if (is_refseq_P != NULL) {
      *is_refseq_P = TRUE;
    }
    return TRUE;
  }

  if (IsDbxrefInList (name, legalSrcDbXrefs,
                      sizeof (legalSrcDbXrefs) / sizeof (legalSrcDbXrefs [0]) - 1,
                      is_badcap_P, goodcapP)) {
    if (is_source_P != NULL) {
      *is_source_P = TRUE;
    }
    return TRUE;
  }

  if (IsDbxrefInList (name, legalDbXrefs,
                      sizeof (legalDbXrefs) / sizeof (legalDbXrefs [0]) - 1,
                      is_badcap_P, goodcapP)) {
    return TRUE;
  }

  return FALSE;
}


/* These functions are for testing dbxrefs */

static ValNodePtr MakeDbxrefList (void)
{
  ValNodePtr dbxref_list = NULL;
  Int4 i;
  DbtagPtr dbtag;

  for (i = 0; legalDbXrefs [i] != NULL; i++) {
    dbtag = DbtagNew ();
    dbtag->db = StringSave (legalDbXrefs [i]);
    dbtag->tag = ObjectIdNew ();
    dbtag->tag->id = 42;
    ValNodeAddPointer (&dbxref_list, 0, dbtag);
  }

  /* legalSrcDbXrefs is contained within legalDbXrefs */

  for (i = 0; legalRefSeqDbXrefs [i] != NULL; i++) {
    dbtag = DbtagNew ();
    dbtag->db = StringSave (legalRefSeqDbXrefs [i]);
    dbtag->tag = ObjectIdNew ();
    dbtag->tag->id = 42;
    ValNodeAddPointer (&dbxref_list, 0, dbtag);
  }

  return dbxref_list;
}

static void AddDbxrefsToBioSource (BioSourcePtr biop)
{
  if (biop == NULL) return;
  if (biop->org == NULL)
  {
    biop->org = OrgRefNew();
  }

  ValNodeLink (&(biop->org->db), MakeDbxrefList());
}

static void AddDbxrefsToSeqFeat (SeqFeatPtr sfp)
{
  if (sfp == NULL) return;
  ValNodeLink (&(sfp->dbxref), MakeDbxrefList());
}

NLM_EXTERN void AddAllDbxrefsToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqFeatPtr  sfp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    AddDbxrefsToBioSource (sdp->data.ptrvalue);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  if (sfp != NULL) {
    AddDbxrefsToBioSource (sfp->data.value.ptrvalue);
    AddDbxrefsToSeqFeat (sfp);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  if (sfp != NULL) {
    AddDbxrefsToSeqFeat (sfp);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
  if (sfp != NULL) {
    AddDbxrefsToSeqFeat (sfp);
  }
}



static CharPtr organellePrefix [] = {
  NULL,
  NULL,
  "Chloroplast ",
  "Chromoplast ",
  "Kinetoplast ",
  "Mitochondrion ",
  "Plastid ",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "Cyanelle ",
  NULL,
  NULL,
  "Nucleomorph ",
  "Apicoplast ",
  "Leucoplast ",
  "Proplastid ",
  NULL,
  "Hydrogenosome ",
  NULL,
  "Chromatophore "
};

static CharPtr newOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast ",
  "chromoplast ",
  "kinetoplast ",
  "mitochondrion ",
  "plastid ",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "cyanelle ",
  NULL,
  NULL,
  "nucleomorph ",
  "apicoplast ",
  "leucoplast ",
  "proplastid ",
  NULL,
  "hydrogenosome ",
  NULL,
  "chromatophore "
};

NLM_EXTERN CharPtr FormatSourceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  CharPtr            acr = NULL;
  Boolean            addPeriod = TRUE;
  IntAsn2gbJobPtr    ajp;
  CharPtr            ana = NULL;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  CharPtr            com = NULL;
  CharPtr            common = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  CharPtr            gbacr = NULL;
  CharPtr            gbana = NULL;
  GBBlockPtr         gbp = NULL;
  GBSeqPtr           gbseq;
  CharPtr            gbsyn = NULL;
  Uint1              genome;
  CharPtr            met = NULL;
  ValNodePtr         mod = NULL;
  Int2               numacr = 0;
  Int2               numana = 0;
  Int2               numcom = 0;
  Int2               numgbacr = 0;
  Int2               numgbana = 0;
  Int2               numgbsyn = 0;
  Int2               nummet = 0;
  Int2               numsyn = 0;
  OrgModPtr          omp = NULL;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  CharPtr            prefix = " (";
  SeqDescrPtr        sdp;
  CharPtr            second = NULL;
  SeqFeatPtr         sfp;
  CharPtr            str;
  CharPtr            syn = NULL;
  CharPtr            taxname = NULL;
  StringItemPtr      ffstring, temp;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {
      if (dcontext.seqdesctype == Seq_descr_source) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      } else if (dcontext.seqdesctype == Seq_descr_genbank) {
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
      }
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (gbp != NULL) {
    common = gbp->source;
  }

  if (biop != NULL) {
    genome = biop->genome;
    if (genome <= 22) {
      if (ajp->newSourceOrg && (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT)) {
        organelle = newOrganellePrefix [genome];
      } else {
        organelle = organellePrefix [genome];
      }
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      mod = orp->mod;
      onp = orp->orgname;
      if (onp != NULL) {

        if (ajp->newSourceOrg) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            switch (omp->subtype) {
              case ORGMOD_common :
                com = omp->subname;
                numcom++;
                break;
              case ORGMOD_acronym :
                acr = omp->subname;
                numacr++;
                break;
              case ORGMOD_synonym :
                syn = omp->subname;
                numsyn++;
                break;
              case ORGMOD_anamorph :
                ana = omp->subname;
                numana++;
                break;
              case ORGMOD_gb_acronym :
                gbacr = omp->subname;
                numgbacr++;
                break;
              case ORGMOD_gb_anamorph :
                gbana = omp->subname;
                numgbana++;
                break;
              case ORGMOD_gb_synonym :
                gbsyn = omp->subname;
                numgbsyn++;
                break;
              case ORGMOD_metagenome_source :
                met = omp->subname;
                nummet++;
                break;
              default :
                break;
            }
          }

          if (numacr > 1) {
             acr = NULL;
          }
          if (numana > 1) {
             ana = NULL;
          }
          if (numcom > 1) {
             com = NULL;
          }
          if (nummet > 1) {
             met = NULL;
          }
          if (numsyn > 1) {
             syn = NULL;
          }
          if (numgbacr > 1) {
             gbacr = NULL;
          }
          if (numgbana > 1) {
             gbana = NULL;
          }
          if (numgbsyn > 1) {
             gbsyn = NULL;
          }

          if (StringHasNoText (second)) {
            second = met;
          }
          if (StringHasNoText (second)) {
            second = syn;
          }
           if (StringHasNoText (second)) {
             second = acr;
          }
          if (StringHasNoText (second)) {
            if (StringDoesHaveText (ana)) {
              second = ana;
              prefix = " (anamorph: ";
            }
          }
          if (StringHasNoText (second)) {
            second = com;
          }

          if (StringHasNoText (second)) {
            second = gbsyn;
          }
          if (StringHasNoText (second)) {
            second = gbacr;
          }
          if (StringHasNoText (second)) {
            if (StringDoesHaveText (gbana)) {
              second = gbana;
              prefix = " (anamorph: ";
            }
          }
        }
      }
      if (StringHasNoText (second)) {
        second = common;
      }
    }
  }

  /* If the organelle prefix is already on the */
  /* name, don't add it.                       */

  if (StringNCmp (organelle, taxname, StringLen (organelle)) == 0)
    organelle = "";

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    
    temp = FFGetString(ajp);

    if (ajp->newSourceOrg) {

      if (! StringHasNoText (organelle)) {
        FFAddTextToString(temp, NULL, organelle, NULL, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddTextToString(temp, NULL, taxname, NULL, FALSE, FALSE, TILDE_IGNORE);
      if (! StringHasNoText (second)) {
        FFAddTextToString(temp, prefix, second, ")", FALSE, FALSE, TILDE_IGNORE);
      }
      addPeriod = FALSE;

    } else {
      FFAddTextToString(temp, NULL, common, NULL, FALSE, FALSE, TILDE_IGNORE);
      while (mod != NULL) {
        str = (CharPtr) mod->data.ptrvalue;
        if (! StringHasNoText (str)) {
          FFAddTextToString(temp, " ", str, NULL, FALSE, FALSE, TILDE_IGNORE);
        }
        mod = mod->next;
      }
    }

    str = FFToCharPtr(temp);
    if (StringCmp (str, ".") == 0) {
      str = MemFree (str);
    }
    FFRecycleString(ajp, temp);
    /* optionally populate gbseq for XML-ized GenBank format */

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    if (gbseq != NULL) {
      gbseq->source = StringSave (str);
    }

    
    FFStartPrint(ffstring, afp->format, 0, 12, "SOURCE", 12, 5, 5, "OS", TRUE);
    if (str != NULL) {
      FFAddTextToString(ffstring, NULL, str, NULL, addPeriod, FALSE, TILDE_TO_SPACES);
    } else {
      FFAddOneChar(ffstring, '.', FALSE);
    }
    
    MemFree (str);

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    FFStartPrint(ffstring, afp->format, 0, 12, "SOURCE", 12, 5, 5, "OS", TRUE);
    FFAddTextToString(ffstring, NULL, taxname, NULL, FALSE, FALSE, TILDE_TO_SPACES);
    if ( StringICmp(taxname, common) != 0 ) {
        FFAddTextToString(ffstring, " (", common, ")", FALSE, FALSE, TILDE_TO_SPACES);
    }
  }
  
  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 0, 5, "OS");
  FFRecycleString(ajp, ffstring);
  return str;
}

NLM_EXTERN CharPtr FormatOrganismBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  Char               ch;
  CharPtr            common = NULL;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBSeqPtr           gbseq;
  Uint1              genome;
  CharPtr            lineage = NULL;
  ObjectIdPtr        oip;
  OrgModPtr          omp;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            str;
  Int4               taxid = -1;
  CharPtr            taxname = NULL;
  CharPtr            tmp;
  CharPtr            ptr;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, temp;
  Char               buf [16];

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;


  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (biop != NULL) {
    genome = biop->genome;
    if (genome <= 22) {
      organelle = organellePrefix [genome];
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      onp = orp->orgname;
      if (onp != NULL) {
        lineage = onp->lineage;
        if (StringHasNoText (lineage)) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            if (omp->subtype == ORGMOD_old_lineage) {
              lineage = omp->subname;
            }
          }
        }
      }
      for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
        dbt = (DbtagPtr) vnp->data.ptrvalue;
        if (dbt == NULL) continue;
        if (StringCmp (dbt->db, "taxon") == 0) {
          oip = dbt->tag;
          if (oip != NULL) {
            taxid = oip->id;
          }
        }
      }
    }
  }

  /* If the organelle prefix is already on the */
  /* name, don't add it.                       */

  if (StringNCmp (organelle, taxname, StringLen (organelle)) == 0)
    organelle = "";

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }
  if (StringHasNoText (lineage)) {
    lineage = "Unclassified.";
  }

  ffstring = FFGetString(ajp);
  temp = FFGetString(ajp);
  if ( ffstring == NULL || temp == NULL ) return NULL;

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    
    FFStartPrint(temp, afp->format, 2, 12, "ORGANISM", 12, 5, 5, "OC", FALSE);
    if (! ajp->newSourceOrg) {
      FFAddOneString(temp, organelle, FALSE, FALSE, TILDE_IGNORE);
    }
    if (StringNICmp (taxname, "Unknown", 7) != 0) {
      if ( GetWWW(ajp) ) { 
        if (taxid != -1) {
          FFAddOneString(temp, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
          FF_Add_NCBI_Base_URL (temp, link_tax);
          FFAddOneString(temp, "id=", FALSE, FALSE, TILDE_IGNORE);
          sprintf (buf, "%ld", (long) taxid);
          FFAddOneString(temp, buf, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, "\">", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString(temp, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
          FF_Add_NCBI_Base_URL (temp, link_tax);
          FFAddOneString(temp, "name=", FALSE, FALSE, TILDE_IGNORE);
          tmp = StringSave (taxname);
          if (tmp != NULL) {
            ptr = tmp;
            ch = *ptr;
            while (ch != '\0') {
              if (IS_WHITESP (ch)) {
                *ptr = '+';
              }
              ptr++;
              ch = *ptr;
            }
            FFAddOneString(temp, tmp, FALSE, FALSE, TILDE_IGNORE);
            MemFree (tmp);
          }
          FFAddOneString(temp, "\">", FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString(temp, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
      }
    } else {
      FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
    }
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    FFRecycleString(ajp, temp);

    temp = FFGetString(ajp);
    FFStartPrint(temp, afp->format, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    FFAddTextToString(temp, NULL, lineage, NULL, TRUE, FALSE, TILDE_TO_SPACES);
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    FFRecycleString(ajp, temp);
    /* optionally populate gbseq for XML-ized GenBank format */

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    if (gbseq != NULL) {
      temp = FFGetString(ajp);
      if (! ajp->newSourceOrg) {
        FFAddOneString(temp, organelle, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
      gbseq->organism = FFToCharPtr(temp);
      gbseq->taxonomy = StringSave (lineage);
      FFRecycleString(ajp, temp);
    }

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFStartPrint(temp, afp->format, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    FFAddTextToString(temp, NULL, lineage, NULL, TRUE, FALSE, TILDE_TO_SPACES);
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "OC");
    FFRecycleString(ajp, temp);
    if ( !StringHasNoText(organelle) ) {
      temp = FFGetString(ajp);
      if ( temp != NULL ) {
        FFStartPrint(temp, afp->format, 12, 12, NULL, 0, 5, 5, "OG", FALSE);
        FFAddTextToString(temp, NULL, organelle, NULL, TRUE, FALSE, TILDE_TO_SPACES);
        FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "OG");
        FFRecycleString(ajp, temp);
      }
    }
  }
  
  str = FFToCharPtr(ffstring);
  FFRecycleString(ajp, ffstring);
  return str;
}

/* A tilde is not an EOL if it is found in a string of the form:    */
/* /~alpahnumdot/ where alphanumdot is either alpha numeric or '.' */
/*                                                                 */
/* str points to the tilde in question.                            */
static Boolean IsTildeEOL(CharPtr str) {
  CharPtr ptr;

  if ( *(str - 1) != '/' ) return TRUE;

  ++str;

  
  for ( ptr = str; 
    IS_ALPHANUM(*ptr) || *ptr == '_' || *ptr == '-' || *ptr == '.';
    ++ptr) continue;

  return *ptr == '/' ? FALSE : TRUE;
}

/* returns a pointer to the first character past the url */
static CharPtr FindUrlEnding(CharPtr str) {
  CharPtr ptr;

  for ( ptr = str;
        !IS_WHITESP(*ptr) && *ptr != '\0' && *ptr != '(' && *ptr != '\"';
        ++ptr  ) {
    if ( *ptr == '~' ) {
      if ( IsTildeEOL(ptr) ) break;
    }
  }

  --ptr;

  /* back up over any trailing periods, commas, or parentheses */
  while ( (*ptr == '.') || (*ptr == ',') || (*ptr == ')') ) {
    --ptr;
  }

  ++ptr;

  return ptr;
}

static Boolean CommentHasSuspiciousHtml (
  IntAsn2gbJobPtr ajp,
  CharPtr searchString
)

{
  Char        ch;
  CharPtr     ptr;
  Int4        state;
  ValNodePtr  matches;

  if (StringHasNoText (searchString)) return FALSE;

  state = 0;
  ptr = searchString;
  ch = *ptr;

  while (ch != '\0') {
    matches = NULL;
    ch = TO_LOWER (ch);
    state = TextFsaNext (ajp->bad_html_fsa, state, ch, &matches);
    if (matches != NULL) {
      return TRUE;
    }
    ptr++;
    ch = *ptr;
  }

  return FALSE;
}

NLM_EXTERN void AddCommentWithURLlinks (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr prefix,
  CharPtr str,
  CharPtr suffix
)

{
  Char     ch;
  CharPtr  ptr;

  if (GetWWW (ajp) && CommentHasSuspiciousHtml (ajp, str)) {
    if (prefix != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    }
    AddCommentStringWithTildes (ffstring, str);
    if (suffix != NULL) {
      FFAddOneString(ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);
    }
    return;
  }

  while (! StringHasNoText (str)) {
    ptr = StringStr (str, "http://");
    if (ptr == NULL) {
      ptr = StringStr (str, "https://");
    }
    if (ptr == NULL) {
      if (prefix != NULL) {
        FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      }
      AddCommentStringWithTildes (ffstring, str);
      if (suffix != NULL) {
        FFAddOneString(ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);
      }
      return;
    }

    *ptr = '\0';
    AddCommentStringWithTildes (ffstring, str); 
    *ptr = 'h';

    str = ptr;
    ptr = FindUrlEnding(str);


    ch = *ptr;
    *ptr = '\0';
    if ( GetWWW(ajp) ) {
      FFAddTextToString(ffstring, "<a href=\"", str, "\">", FALSE, FALSE, TILDE_IGNORE);
      FFAddTextToString(ffstring, NULL, str, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }

    *ptr = ch;
    str = ptr;
  }
}

static void CatenateCommentInGbseq (
  GBSeqPtr gbseq,
  CharPtr str,
  Boolean compress
)

{
  Char     ch;
  CharPtr  tmp;

  if (gbseq == NULL || StringHasNoText (str)) return;

  if (StringNCmp (str, "COMMENT     ", 12) == 0) {
    str += 12;
  }
  if (gbseq->comment == NULL) {
    gbseq->comment = StringSave (str);
  } else {
    tmp = (CharPtr) MemNew (StringLen (gbseq->comment) + StringLen (str) + 4);
    StringCpy (tmp, gbseq->comment);
    StringCat (tmp, "; ");
    StringCat (tmp, str);
    gbseq->comment = MemFree (gbseq->comment);
    gbseq->comment = tmp;
  }

  tmp = gbseq->comment;
  if (tmp == NULL) return;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *tmp = ' ';
    }
    tmp++;
    ch = *tmp;
  }
  TrimSpacesAroundString (gbseq->comment);
  if (compress) {
    Asn2gnbkCompressSpaces (gbseq->comment);
  }
}


NLM_EXTERN CharPtr FormatCommentBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Boolean            add_period;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  CommentBlockPtr    cbp;
  CharPtr            db;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBSeqPtr           gbseq;
  size_t             len;
  ObjectIdPtr        oip;
  CharPtr            prefix;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  Char               sfx [32];
  CharPtr            str;
  CharPtr            suffix;
  CharPtr            title;
  StringItemPtr      ffstring;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  cbp = (CommentBlockPtr) bbp;

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  /* some comments are allocated (along with possible first COMMENT label) */

  if (! StringHasNoText (bbp->string)) {
    str = StringSave (bbp->string);
    CatenateCommentInGbseq (gbseq, str, TRUE);
    return str;
  }

  title = NULL;
  prefix = NULL;
  suffix = NULL;
  add_period = FALSE;
  sfx [0] = '\0';

  if (bbp->itemtype == OBJ_SEQDESC) {

    /* usually should reference comment, maploc, or region descriptor IDs */

    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {

      if (dcontext.seqdesctype == Seq_descr_comment) {

        title = (CharPtr) sdp->data.ptrvalue;

      } else if (dcontext.seqdesctype == Seq_descr_maploc) {

        dbt = (DbtagPtr) sdp->data.ptrvalue;
        if (dbt != NULL) {
          db = dbt->db;
          oip = dbt->tag;
          if (oip != NULL) {
            if (oip->str != NULL) {

              title = oip->str;
              prefix = ("Map location: ");

            } else if (db != NULL && oip->id != 0) {

              title = db;
              prefix = ("Map location: (Database ");
              sprintf (sfx, "; id # %ld).", (long) oip->id);
              suffix = sfx;

            }
          }
        }

      } else if (dcontext.seqdesctype == Seq_descr_region) {

        title = (CharPtr) sdp->data.ptrvalue;
        prefix = "Region: ";

      } else if (dcontext.seqdesctype == Seq_descr_name) {

        title = (CharPtr) sdp->data.ptrvalue;
        prefix = "Name: ";

      }
    }

  } else if (bbp->itemtype == OBJ_SEQFEAT) {

    /* also have to deal with comment feature across entire sequence */

    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_COMMENT) {

      title = sfp->comment;
    }
  }

  if (title == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (cbp->first) {
    FFStartPrint (ffstring, afp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
  } else {
    FFStartPrint (ffstring, afp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
  }

  str = StringSave (title);
  TrimSpacesAndJunkFromEnds (str, TRUE);
  if (! IsEllipsis (str)) {
    s_RemovePeriodFromEnd (str);
    len = StringLen (str);
    if (len > 0 && str [len - 1] != '.') {
      add_period = TRUE;
    }
  }
  AddCommentWithURLlinks(ajp, ffstring, prefix, str, suffix);
  /*
  if ( GetWWW(ajp) && prefix == NULL && suffix == NULL) {
    
    AddCommentWithURLlinks (ffstring, str);
  } else {
    FFAddTextToString (ffstring, prefix, str, suffix, FALSE, TRUE, TILDE_OLD_EXPAND);
  }
  */
  if (add_period) {
    FFAddOneChar (ffstring, '.',FALSE);
  }
  MemFree (str);

  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 5, 5, "CC");

  /*
  CatenateCommentInGbseq (gbseq, str);
  */
  CatenateCommentInGbseq (gbseq, title, FALSE);

  FFRecycleString(ajp, ffstring);
  return str;
}

/* format features section */

static Boolean is_real_id (
  SeqIdPtr sip,
  SeqIdPtr this_sip
)

{
  BioseqPtr  bsp;

  if (sip == NULL || this_sip == NULL) return FALSE;

  if (! SeqIdIn (sip, this_sip)) {
    bsp = BioseqFind (sip);
    if (bsp == NULL) return TRUE;  /* ??? */
    if (bsp->repr == Seq_repr_virtual) return FALSE;
  }

  return TRUE;
}

static Boolean FlatVirtLoc (
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  SeqPntPtr  spp;

  if (bsp == NULL || location == NULL) return FALSE;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return TRUE;
      sip = sintp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return TRUE;
      sip = spp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    default :
      break;
  }

  return FALSE;
}

static Uint1    id_order [NUM_SEQID];
static Boolean  order_initialized = FALSE;

static CharPtr lim_str [5] = { "", ">","<", ">", "<" };

NLM_EXTERN Boolean GetAccnVerFromServer (Int4 gi, CharPtr buf)

{
  AccnVerLookupFunc  func;
  SeqMgrPtr          smp;
  CharPtr            str;

  if (buf == NULL) return FALSE;
  *buf = '\0';
  smp = SeqMgrWriteLock ();
  if (smp == NULL) return FALSE;
  func = smp->accn_ver_lookup_func;
  SeqMgrUnlock ();
  if (func == NULL) return FALSE;
  str = (*func) (gi);
  if (str == NULL) return FALSE;
  if (StringLen (str) < 40) {
    StringCpy (buf, str);
  }
  MemFree (str);
  return TRUE;
}


/******************************************************************************/
/*                              FFFlatLoc functions  .                          */
/******************************************************************************/

static Boolean FF_FlatNullAhead (
  BioseqPtr bsp,
  ValNodePtr location
)

{
  SeqLocPtr  next;

  if (bsp == NULL || location == NULL) return FALSE;

  next = location->next;
  if (next == NULL) return TRUE;
  if (next->choice == SEQLOC_NULL) return TRUE;
  if (FlatVirtLoc (bsp, next)) return TRUE;

  return FALSE;
}



static void FlatLocSeqId (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip
)

{
  BioseqPtr    bsp;
  Char         buf [40];
  ObjectIdPtr  oip;
  SeqIdPtr     use_id = NULL;
  Boolean      was_lock = FALSE;

  if (ffstring == NULL || sip == NULL) return;

  buf [0] = '\0';
  bsp = BioseqFind (sip);
  if (bsp != NULL) {
    use_id = SeqIdSelect (bsp->id, id_order, NUM_SEQID);
  } else if (sip->choice == SEQID_GI) {
    if (GetAccnVerFromServer (sip->data.intvalue, buf)) {
      FFAddTextToString(ffstring, NULL, buf, ":", FALSE, FALSE, TILDE_IGNORE);
      /*AddValNodeString (head, NULL, buf, ":");*/
      return;
    }
    use_id = GetSeqIdForGI (sip->data.intvalue);
  }
  if (use_id == NULL && bsp == NULL) {
    bsp = BioseqLockById (sip);
    was_lock = TRUE;
    if (bsp != NULL) {
      use_id = SeqIdSelect (bsp->id, id_order, NUM_SEQID);
    }
  }
  if (use_id != NULL) {
    SeqIdWrite (use_id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
    if (use_id->choice == SEQID_GI) {
      ajp->relModeError = TRUE;
    }
  } else if (sip->choice == SEQID_GI) {
    SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    ajp->relModeError = TRUE;
  } else {
    SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
    if (sip->choice == SEQID_GI) {
      ajp->relModeError = TRUE;
    }
  }
  if (was_lock) {
    BioseqUnlock (bsp);
  }
  if (StringHasNoText (buf)) {
    StringCpy (buf, "?00000");
    ajp->relModeError = TRUE;
    if (use_id != NULL && use_id->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) use_id->data.ptrvalue;
      if (oip != NULL && (! StringHasNoText (oip->str))) {
        StringNCpy_0 (buf, oip->str, 13);
      }
    }
  }
  FFAddTextToString(ffstring, NULL, buf, ":", FALSE, FALSE, TILDE_IGNORE);
}



static void FlatLocCaret (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (ffstring == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (ajp, ffstring, sip);
  }

  buf [0] = '\0';
  point++; /* orginal FlatLocHalfCaret was called with point + 1 */

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)..(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) point,
                 (long) point,
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "%ld^%ld",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "%ld^%ld",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        if (fuzz->a == 3) { /* space to right */
          sprintf (buf, "%ld^%ld", (long) (point), (long) (point + 1));
        } else if (fuzz->a == 4 && point > 1) { /* space to left */
          sprintf (buf, "%ld^%ld", (long) (point - 1), (long) point);
        } else {
          index = (Uint1) fuzz->a;
          if (index > 4) {
            index = 0;
          }
          sprintf (buf, "%s%ld", lim_str [index], (long) point);
        }
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
}


static void FlatLocPoint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (ffstring == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (ajp, ffstring, sip);
  }

  buf [0] = '\0';
  point++;

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        index = (Uint1) fuzz->a;
        if (index > 4) {
          index = 0;
        }
        sprintf (buf, "%s%ld", lim_str [index], (long) point);
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
}


static void FlatLocElement (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  Boolean     minus_strand = FALSE;
  SeqBondPtr  sbp;
  SeqIntPtr   sintp;
  SeqIdPtr    sip;
  SeqPntPtr   spp;
  BioseqPtr   wholebsp;

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return;
      wholebsp = BioseqFind (sip);
      if (wholebsp == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        FlatLocPoint (ajp, ffstring, sip, bsp->id, 0, NULL);
        if (bsp->length > 0) {
          FFAddOneString(ffstring, "..", FALSE, FALSE, TILDE_IGNORE);
          FlatLocPoint (ajp, ffstring, NULL, bsp->id, bsp->length - 1, NULL);
        }
      }
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return;
      sip = sintp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (sintp->strand == Seq_strand_minus);
        if (minus_strand) {
          FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
        }
        FlatLocPoint (ajp, ffstring, sip, bsp->id, sintp->from, sintp->if_from);
        if (sintp->to > 0 &&
            (sintp->to != sintp->from ||
             sintp->if_from != NULL ||
             sintp->if_to != NULL)) {
          FFAddOneString(ffstring, "..", FALSE, FALSE, TILDE_IGNORE);
          FlatLocPoint (ajp, ffstring, NULL, bsp->id, sintp->to, sintp->if_to);
        }
        if (minus_strand) {
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
        }
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (spp->strand == Seq_strand_minus);
        if (minus_strand) {
          FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
        }
        if (spp->fuzz != NULL) {
          FlatLocCaret (ajp, ffstring, sip, bsp->id, spp->point, spp->fuzz);
        } else {
          FlatLocPoint (ajp, ffstring, sip, bsp->id, spp->point, NULL);
        }
        if (minus_strand) {
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
        }
      }
      break;
    case SEQLOC_BOND :
      sbp = (SeqBondPtr) location->data.ptrvalue;
      if (sbp == NULL) return;
      spp = sbp->a;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      FFAddOneString(ffstring, "bond(", FALSE, FALSE, TILDE_IGNORE);
      FlatLocPoint (ajp, ffstring, sip, bsp->id, spp->point, spp->fuzz);
      spp = sbp->b;
      if (spp != NULL) {
        FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
        FlatLocPoint (ajp, ffstring, NULL, bsp->id, spp->point, spp->fuzz);
      }
      FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
      break;
    default :
      /* unexpected internal complex type or unimplemented SEQLOC_FEAT */
      return;
  }
}



static void FF_FlatPackedPoint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  PackSeqPntPtr pspp,
  BioseqPtr bsp
)

{
  Uint1  dex;

  if (ffstring == NULL || pspp == NULL || bsp == NULL) return;

  for (dex = 0; dex < pspp->used; dex++) {
    FlatLocPoint (ajp, ffstring, pspp->id, bsp->id, pspp->pnts [dex], pspp->fuzz);
  }
}


static void FF_DoFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
);

static void FF_GroupFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  CharPtr prefix,
  Boolean is_flat_order
)

{
  Boolean        found_non_virt = FALSE;
  SeqIdPtr       hold_next;
  Int2           parens = 1;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;
  Boolean        special_mode = FALSE; /* join in order */

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  /* prefix will have the first parenthesis */

  FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);

  for (slp = (SeqLocPtr) location->data.ptrvalue; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) {
      if (slp != location && slp->next != NULL) {
        if (special_mode) {
          special_mode = FALSE;
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
          parens--;
        }
      }
      continue;
    }

    if (found_non_virt && slp->choice != SEQLOC_EMPTY && slp->choice != SEQLOC_NULL) {
      FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
    }

    switch (slp->choice) {
      case SEQLOC_WHOLE :
      case SEQLOC_PNT :
      case SEQLOC_BOND :
      case SEQLOC_FEAT :
        found_non_virt = TRUE;
        if (FlatVirtLoc (bsp, slp)) {
          if (slp != location && slp->next != NULL) {
            if (special_mode) {
              special_mode = FALSE;
              FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
              parens--;
            }
          }
        } else {
          FlatLocElement (ajp, ffstring, bsp, slp);
        }
        break;
      case SEQLOC_INT :
        found_non_virt = TRUE;
        if (is_flat_order && (! FF_FlatNullAhead (bsp, slp))) {
          special_mode = TRUE;
          FFAddOneString(ffstring, "join(", FALSE, FALSE, TILDE_IGNORE);
          parens++;
        }
        FlatLocElement (ajp, ffstring, bsp, slp);
        break;
      case SEQLOC_PACKED_PNT :
        found_non_virt = TRUE;
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FF_FlatPackedPoint (ajp, ffstring, pspp, bsp);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        found_non_virt = TRUE;
        hold_next = slp->next;
        slp->next = NULL;
        FF_DoFlatLoc (ajp, ffstring, bsp, slp, FALSE);
        slp->next = hold_next;
        break;
      default :
        break;
    }

  }

  while (parens > 0) {
    FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
    parens--;
  }
}




static void FF_DoFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
)

{
  Boolean        found_null;
  SeqLocPtr      next_loc;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  /* deal with complement of entire location */

  if (ok_to_complement && SeqLocStrand (location) == Seq_strand_minus) {
    slp = AsnIoMemCopy ((Pointer) location,
                        (AsnReadFunc) SeqLocAsnRead,
                        (AsnWriteFunc) SeqLocAsnWrite);
    if (slp != NULL) {
      SeqLocRevCmp (slp);
      FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
      FF_DoFlatLoc (ajp, ffstring, bsp, slp, FALSE);
      FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
    }
    SeqLocFree (slp);
    return;
  }

  /* handle each location component */

  for (slp = location; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) continue;

    /* print comma between components */

    if (slp != location) {
      FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
    }

    switch (slp->choice) {
      case SEQLOC_MIX :
      case SEQLOC_PACKED_INT :
        found_null = FALSE;
        for (next_loc = (SeqLocPtr) slp->data.ptrvalue;
         next_loc != NULL;
         next_loc = next_loc->next) {
          if (next_loc->choice == SEQLOC_NULL ||
              FlatVirtLoc (bsp, next_loc) /* ||
              LocationHasNullsBetween (slp) */ )
            found_null = TRUE;
        }
        if (found_null) {
          FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "order(", TRUE);
        } else {
          FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "join(", FALSE);
        }
        break;
      case SEQLOC_EQUIV :
        FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "one-of(", FALSE);
        break;
      case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FF_FlatPackedPoint (ajp, ffstring, pspp, bsp);
        }
        break;
      default :
        FlatLocElement (ajp, ffstring, bsp, slp);
        break;
    }

  }
}




NLM_EXTERN CharPtr FFFlatLoc (
  IntAsn2gbJobPtr ajp,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean masterStyle
)

{
  Boolean     hasNulls;
  IntFuzzPtr  fuzz = NULL;
  SeqLocPtr   loc;
  Boolean     noLeft;
  Boolean     noRight;
  Uint1       num = 1;
  SeqPntPtr   spp;
  CharPtr     str;
  SeqLocPtr   tmp;
  StringItemPtr ffstring = NULL;

  if (bsp == NULL || location == NULL) return NULL;

  ffstring = FFGetString(ajp);

  if (! order_initialized) {
    id_order [SEQID_GENBANK] = num++;
    id_order [SEQID_EMBL] = num++;
    id_order [SEQID_DDBJ] = num++;
    id_order [SEQID_OTHER] = num++;
    id_order [SEQID_TPG] = num++;
    id_order [SEQID_TPE] = num++;
    id_order [SEQID_TPD] = num++;
    id_order [SEQID_GPIPE] = num++;
    id_order [SEQID_GIBBSQ] = num++;
    id_order [SEQID_GIBBMT] = num++;
    id_order [SEQID_PRF] = num++;
    id_order [SEQID_PDB] = num++;
    id_order [SEQID_PIR] = num++;
    id_order [SEQID_SWISSPROT] = num++;
    id_order [SEQID_PATENT] = num++;
    id_order [SEQID_GI] = num++;;
    id_order [SEQID_GENERAL] = num++;
    id_order [SEQID_LOCAL] = num++;
    id_order [SEQID_GIIM] = num++;
    order_initialized = TRUE;
  }

  if (masterStyle) {

    /* map location from parts to segmented bioseq */

    if (location->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp != NULL) {
        fuzz = spp->fuzz;
      }
    }

    CheckSeqLocForPartial (location, &noLeft, &noRight);
    hasNulls = LocationHasNullsBetween (location);
    loc = SeqLocMergeEx (bsp, location, NULL, FALSE, TRUE, FALSE, hasNulls);
    if (loc == NULL) {
      tmp = TrimLocInSegment (bsp, location, &noLeft, &noRight);
      loc = SeqLocMergeEx (bsp, tmp, NULL, FALSE, TRUE, FALSE, hasNulls);
      SeqLocFree (tmp);
    }
    if (loc == NULL) {
      return StringSave ("?");
    }
    FreeAllFuzz (loc);
    SetSeqLocPartial (loc, noLeft, noRight);

    if (loc->choice == SEQLOC_PNT && fuzz != NULL) {
      spp = (SeqPntPtr) loc->data.ptrvalue;
      if (spp != NULL && spp->fuzz == NULL) {
        spp->fuzz = AsnIoMemCopy ((Pointer) fuzz,
                                  (AsnReadFunc) IntFuzzAsnRead,
                                  (AsnWriteFunc) IntFuzzAsnWrite);
      }
    }

    FF_DoFlatLoc (ajp, ffstring, bsp, loc, TRUE);

    SeqLocFree (loc);

  } else {
    FF_DoFlatLoc (ajp, ffstring, bsp, location, TRUE);
  }

  str = FFToCharPtr(ffstring);
  FFRecycleString(ajp, ffstring);
  return str;
}




static void PromoteSeqId (SeqIdPtr sip, Pointer userdata)

{
  SeqIdPtr  bestid, newid, oldid;

  bestid = (SeqIdPtr) userdata;

  newid = SeqIdDup (bestid);
  if (newid == NULL) return;

  oldid = ValNodeNew (NULL);
  if (oldid == NULL) return;

  MemCopy (oldid, sip, sizeof (ValNode));
  oldid->next = NULL;

  sip->choice = newid->choice;
  sip->data.ptrvalue = newid->data.ptrvalue;

  SeqIdFree (oldid);
  ValNodeFree (newid);

  SeqIdStripLocus (sip);
}

NLM_EXTERN SeqLocPtr SeqLocReMapEx (
  SeqIdPtr newid,
  SeqLocPtr seq_loc,
  SeqLocPtr location,
  Int4 offset,
  Boolean rev,
  Boolean masterStyle
)

{
  BioseqPtr    bsp;
  Boolean      hasNulls;
  IntFuzzPtr   fuzz = NULL;
  SeqLocPtr    loc;
  Boolean      noLeft;
  Boolean      noRight;
  SeqEntryPtr  scope;
  SeqIdPtr     sip;
  SeqLocPtr    slp = NULL;
  SeqPntPtr    spp;
  SeqLocPtr    tmp;

  if (newid == NULL || seq_loc == NULL || location == NULL) return NULL;

  if (masterStyle) {

    sip = SeqLocId (seq_loc);
    if (sip == NULL) return NULL;
    bsp = BioseqFind (sip);
    if (bsp == NULL) {
      scope = SeqEntrySetScope (NULL);
      bsp = BioseqFind (sip);
      SeqEntrySetScope (scope);
    }
    if (bsp == NULL) return NULL;
    sip = SeqIdFindBest (bsp->id, 0);

    /* map location from parts to segmented bioseq */

    if (location->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp != NULL) {
        fuzz = spp->fuzz;
      }
    }

    CheckSeqLocForPartial (location, &noLeft, &noRight);
    hasNulls = LocationHasNullsBetween (location);
    loc = SeqLocMerge (bsp, location, NULL, FALSE, TRUE, hasNulls);
    if (loc == NULL) {
      tmp = TrimLocInSegment (bsp, location, &noLeft, &noRight);
      loc = SeqLocMerge (bsp, tmp, NULL, FALSE, TRUE, hasNulls);
      SeqLocFree (tmp);
    }
    if (loc == NULL) {
      return NULL;
    }
    FreeAllFuzz (loc);
    SetSeqLocPartial (loc, noLeft, noRight);

    if (loc->choice == SEQLOC_PNT && fuzz != NULL) {
      spp = (SeqPntPtr) loc->data.ptrvalue;
      if (spp != NULL && spp->fuzz == NULL) {
        spp->fuzz = AsnIoMemCopy ((Pointer) fuzz,
                                  (AsnReadFunc) IntFuzzAsnRead,
                                  (AsnWriteFunc) IntFuzzAsnWrite);
      }
    }

    scope = SeqEntrySetScope (NULL);
    slp = SeqLocReMap (newid, seq_loc, loc, offset, rev);
    SeqEntrySetScope (scope);

    SeqLocFree (loc);

    VisitSeqIdsInSeqLoc (slp, (Pointer) sip, PromoteSeqId);
  } else {

    scope = SeqEntrySetScope (NULL);
    slp = SeqLocReMap (newid, seq_loc, location, offset, rev);
    SeqEntrySetScope (scope);
  }

  return slp;
}


/******************************************************************************/
/*                            End FFFlatLoc functions.                          */
/******************************************************************************/



static void SubSourceToQualArray (
  SubSourcePtr ssp,
  QualValPtr qvp
)

{
  SourceType  idx;
  Uint1       subtype;

  if (ssp == NULL || qvp == NULL) return;

  while (ssp != NULL) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 41;
    }
    if (subtype < 42) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].ssp == NULL) {
          qvp [idx].ssp = ssp;
        }
      }
    }
    ssp = ssp->next;
  }
}

NLM_EXTERN SourceType orgModToSourceIdx [41] = {
  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_strain,
  SCQUAL_sub_strain,
  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_variety,
  SCQUAL_serotype,
  SCQUAL_serogroup,
  SCQUAL_serovar,
  SCQUAL_cultivar,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_isolate,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,
  SCQUAL_spec_or_nat_host,
  SCQUAL_sub_species,
  SCQUAL_specimen_voucher,
  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_ecotype,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,
  SCQUAL_gb_acronym,
  SCQUAL_gb_anamorph,
  SCQUAL_gb_synonym,
  SCQUAL_culture_collection,
  SCQUAL_bio_material,
  SCQUAL_metagenome_source,
  SCQUAL_old_lineage,
  SCQUAL_old_name,
  SCQUAL_orgmod_note
};

static void OrgModToQualArray (
  OrgModPtr omp,
  QualValPtr qvp
)

{
  SourceType  idx;
  Uint1       subtype;

  if (omp == NULL || qvp == NULL) return;

  while (omp != NULL) {
    subtype = omp->subtype;
    if (subtype == 253) {
      subtype = 38;
    } else if (subtype == 254) {
      subtype = 39;
    } else if (subtype == 255) {
      subtype = 40;
    }
    if (subtype < 41) {
      idx = orgModToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].omp == NULL) {
          qvp [idx].omp = omp;
        }
      }
    }
    omp = omp->next;
  }
}

static CharPtr organelleQual [] = {
  NULL,
  NULL,
  "/organelle=\"plastid:chloroplast\"",
  "/organelle=\"plastid:chromoplast\"",
  "/organelle=\"mitochondrion:kinetoplast\"",
  "/organelle=\"mitochondrion\"",
  "/organelle=\"plastid\"",
  "/macronuclear",
  NULL,
  "/plasmid=\"\"",
  "/transposon=\"\"",
  "/insertion_seq=\"\"",
  "/organelle=\"plastid:cyanelle\"",
  "/proviral",
  NULL,
  "/organelle=\"nucleomorph\"",
  "/organelle=\"plastid:apicoplast\"",
  "/organelle=\"plastid:leucoplast\"",
  "/organelle=\"plastid:proplastid\"",
  NULL,
  "/organelle=\"hydrogenosome\"",
  NULL,
  "/organelle=\"chromatophore\""
};

NLM_EXTERN Boolean StringIsJustQuotes (
  CharPtr str
)

{
  Nlm_Uchar  ch;    /* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ' && ch != '"' && ch != '\'') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static CharPtr RemoveAllSpaces (
  CharPtr str
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch != ' ') {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

NLM_EXTERN void AddFeatureToGbseq (
  GBSeqPtr gbseq,
  GBFeaturePtr gbfeat,
  CharPtr str,
  SeqFeatPtr sfp
)

{
  Char            ch;
  CharPtr         copy;
  GBQualifierPtr  gbqual;
  GBQualifierPtr  last = NULL;
  CharPtr         ptr;
  CharPtr         qual;
  CharPtr         tmp;
  CharPtr         val;

  if (gbseq == NULL || gbfeat == NULL || StringHasNoText (str)) return;

  copy = StringSave (str);

  /* link in reverse order, to be reversed in slash block */

  gbfeat->next = gbseq->feature_table;
  gbseq->feature_table = gbfeat;

  /* now parse qualifiers */

  ptr = StringStr (copy, "                     /");
  while (ptr != NULL) {
    qual = ptr + 22;
    val = qual;
    ch = *val;
    while (ch != '=' && ch != '\n' && ch != '\0') {
      val++;
      ch = *val;
    }
    /*
    val = StringChr (qual, '=');
    if (val == NULL) {
      val = StringChr (qual, '\n');
    }
    */
    if (ch != '\0' /* val != NULL */) {
      *val = '\0';
      val++;
      if (ch == '=') {
        tmp = val;
        if (*val == '"') {
          val++;
          tmp = val;
          ch = *tmp;
          while (ch != '"' && ch != '\0') {
            tmp++;
            ch = *tmp;
          }
        }
        ptr = StringStr (tmp, "\n                     /");
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
        }
      } else {
        ptr = StringStr (val, "                     /");
        val = NULL;
      }
      gbqual = GBQualifierNew ();
      if (gbqual != NULL) {
        gbqual->name = StringSave (qual);
        if (! StringHasNoText (val)) {
          gbqual->value = StringSave (val);
          CleanQualValue (gbqual->value);
          Asn2gnbkCompressSpaces (gbqual->value);
          if (sfp != NULL) {
            if (sfp->data.choice == SEQFEAT_CDREGION &&
                StringICmp (qual, "translation") == 0) {
              RemoveAllSpaces (gbqual->value);
            } else if (sfp->data.choice == SEQFEAT_RNA &&
                       StringICmp (qual, "transcription") == 0) {
              RemoveAllSpaces (gbqual->value);
            } else if (sfp->data.choice == SEQFEAT_PROT &&
                       StringICmp (qual, "peptide") == 0) {
              RemoveAllSpaces (gbqual->value);
            }
          }
        }
      }
    } else {
      gbqual = GBQualifierNew ();
      if (gbqual != NULL) {
        gbqual->name = StringSave (qual);
      }
    }
    if (gbfeat->quals == NULL) {
      gbfeat->quals = gbqual;
    } else if (last != NULL) {
      last->next = gbqual;
    }
    last = gbqual;
  }

  MemFree (copy);
}

NLM_EXTERN CharPtr GetMolTypeQual (
  BioseqPtr bsp
)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  if (bsp == NULL) return NULL;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return NULL;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return NULL;

  switch (mip->biomol) {
    case 0 :
      switch (bsp->mol) {
        case Seq_mol_dna :
          return "unassigned DNA";
        case Seq_mol_rna :
          return "unassigned RNA";
        case Seq_mol_na :
          break;
        default :
          break;
      }
      break;
    case MOLECULE_TYPE_GENOMIC :
      switch (bsp->mol) {
        case Seq_mol_dna :
          return "genomic DNA";
        case Seq_mol_rna :
          return "genomic RNA";
        case Seq_mol_na :
          break;
        default :
          break;
      }
      break;
    case MOLECULE_TYPE_PRE_MRNA :
      return "transcribed RNA";
    case MOLECULE_TYPE_MRNA :
      return "mRNA";
    case MOLECULE_TYPE_RRNA :
      return "rRNA";
    case MOLECULE_TYPE_TRNA :
      return "tRNA";
    case MOLECULE_TYPE_SNRNA :
      return "transcribed RNA";
    case MOLECULE_TYPE_SCRNA :
      return "transcribed RNA";
    case MOLECULE_TYPE_PEPTIDE :
      break;
    case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL :
      switch (bsp->mol) {
        case Seq_mol_dna :
          return "other DNA";
        case Seq_mol_rna :
          return "other RNA";
        case Seq_mol_na :
          break;
        default :
          break;
      }
      break;
    case MOLECULE_TYPE_GENOMIC_MRNA_MIX :
      break;
    case MOLECULE_TYPE_CRNA :
      return "viral cRNA";
      break;
    case MOLECULE_TYPE_SNORNA :
      return "transcribed RNA";
      break;
    case MOLECULE_TYPE_TRANSCRIBED_RNA :
      return "transcribed RNA";
      break;
    case MOLECULE_TYPE_NCRNA :
      return "transcribed RNA";
      break;
    case MOLECULE_TYPE_TMRNA :
      return "transcribed RNA";
      break;
    case 255 :
      switch (bsp->mol) {
        case Seq_mol_dna :
          return "other DNA";
        case Seq_mol_rna :
          return "other RNA";
        case Seq_mol_na :
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }

  return NULL;
}

static ValNodePtr ParsePCRPrimerString (
  QualValPtr qvp
)

{
  CharPtr       fwd_primer_seq = NULL;
  CharPtr       rev_primer_seq = NULL;
  CharPtr       fwd_primer_name = NULL;
  CharPtr       rev_primer_name = NULL;
  SubSourcePtr  ssp;

  if (qvp == NULL) return NULL;

  ssp = qvp [SCQUAL_fwd_primer_seq].ssp;
  if (ssp != NULL) {
    fwd_primer_seq = ssp->name;
  }
  ssp = qvp [SCQUAL_rev_primer_seq].ssp;
  if (ssp != NULL) {
    rev_primer_seq = ssp->name;
  }
  ssp = qvp [SCQUAL_fwd_primer_name].ssp;
  if (ssp != NULL) {
    fwd_primer_name = ssp->name;
  }
  ssp = qvp [SCQUAL_rev_primer_name].ssp;
  if (ssp != NULL) {
    rev_primer_name = ssp->name;
  }

  return ParsePCRStrings (fwd_primer_seq, rev_primer_seq, fwd_primer_name, rev_primer_name);
}

static ValNodePtr ParseColonString (
  CharPtr strs,
  Boolean multiple
)

{
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  str = tmp;
  len = StringLen (str);
  if (len > 1 && StringChr (str, ':') != NULL /* && multiple */) {
    while (StringDoesHaveText (str)) {
      ptr = StringChr (str, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
      TrimSpacesAroundString (str);
      ValNodeCopyStr (&head, 0, str);
      str = ptr;
    }
  } else {
    ValNodeCopyStr (&head, 0, str);
  }

  MemFree (tmp);
  return head;
}

static void PrintHalfPrimer (
  ValNodePtr PNTR headp,
  CharPtr name,
  CharPtr seq,
  CharPtr nm_label,
  CharPtr sq_label,
  CharPtr prefix,
  Boolean name_only_ok,
  Boolean multiple
)

{
  ValNodePtr  name_list, seq_list, name_vnp, seq_vnp;
  CharPtr     str;

  name_list = ParseColonString (name, multiple);
  seq_list = ParseColonString (seq, multiple);

  name_vnp = name_list;
  seq_vnp = seq_list;
  if (seq_vnp != NULL) {
    while (seq_vnp != NULL) {
      if (name_vnp != NULL) {
        str = (CharPtr) name_vnp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          ValNodeCopyStr (headp, 0, prefix);
          ValNodeCopyStr (headp, 0, nm_label);
          ValNodeCopyStr (headp, 0, str);
          prefix = ", ";
        }
        name_vnp = name_vnp->next;
      }
      str = (CharPtr) seq_vnp->data.ptrvalue;
      if (StringDoesHaveText (str)) {
        ValNodeCopyStr (headp, 0, prefix);
        ValNodeCopyStr (headp, 0, sq_label);
        ValNodeCopyStr (headp, 0, str);
        prefix = ", ";
      }
      seq_vnp = seq_vnp->next;
    }
  } else if (name_only_ok) {
    while (name_vnp != NULL) {
      str = (CharPtr) name_vnp->data.ptrvalue;
      if (StringDoesHaveText (str)) {
        ValNodeCopyStr (headp, 0, prefix);
        ValNodeCopyStr (headp, 0, nm_label);
        ValNodeCopyStr (headp, 0, str);
        prefix = ", ";
      }
      name_vnp = name_vnp->next;
    }
  }

  ValNodeFreeData (name_list);
  ValNodeFreeData (seq_list);
}

static CharPtr NextPCRPrimerString (
  PcrSetPtr psp,
  Boolean isInNote,
  Boolean multiple
)

{
  ValNodePtr  head = NULL, vnp;
  CharPtr     prefix = NULL;
  CharPtr     str;

  if (psp == NULL) return NULL;

  if (StringHasNoText (psp->fwd_seq) || StringHasNoText (psp->rev_seq)) {
    if (isInNote) {
      /*
      if (StringDoesHaveText (psp->fwd_name)) {
        ValNodeCopyStr (&head, 0, prefix);
        ValNodeCopyStr (&head, 0, "fwd_name: ");
        ValNodeCopyStr (&head, 0, psp->fwd_name);
        prefix = ", ";
      }

      if (StringDoesHaveText (psp->fwd_seq)) {
        ValNodeCopyStr (&head, 0, prefix);
        ValNodeCopyStr (&head, 0, "fwd_seq: ");
        ValNodeCopyStr (&head, 0, psp->fwd_seq);
        prefix = ", ";
      }

      if (StringDoesHaveText (psp->rev_name)) {
        ValNodeCopyStr (&head, 0, prefix);
        ValNodeCopyStr (&head, 0, "rev_name: ");
        ValNodeCopyStr (&head, 0, psp->rev_name);
        prefix = ", ";
      }

      if (StringDoesHaveText (psp->rev_seq)) {
        ValNodeCopyStr (&head, 0, prefix);
        ValNodeCopyStr (&head, 0, "rev_seq: ");
        ValNodeCopyStr (&head, 0, psp->rev_seq);
        prefix = ", ";
      }
      */
      PrintHalfPrimer (&head, psp->fwd_name, psp->fwd_seq, "fwd_name: ", "fwd_seq: ", NULL, TRUE, multiple);
      if (head != NULL) {
        prefix = ", ";
      }
      PrintHalfPrimer (&head, psp->rev_name, psp->rev_seq, "rev_name: ", "rev_seq: ", prefix, TRUE, multiple);
    } else {
      return StringSave ("");
    }
  } else {
    if (isInNote) return StringSave ("");

    PrintHalfPrimer (&head, psp->fwd_name, psp->fwd_seq, "fwd_name: ", "fwd_seq: ", NULL, FALSE, multiple);
    PrintHalfPrimer (&head, psp->rev_name, psp->rev_seq, "rev_name: ", "rev_seq: ", ", ", FALSE, multiple);
  }

  if (head != NULL && isInNote) {
    vnp = ValNodeCopyStr (NULL, 0, "PCR_primers=");
    if (vnp != NULL) {
      vnp->next = head;
      head = vnp;
    }
  }

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);
  return str;
}

static void PrintHalfReaction (
  ValNodePtr PNTR headp,
  PCRPrimerPtr primers,
  CharPtr nm_label,
  CharPtr sq_label,
  CharPtr prefix,
  Boolean name_only_ok,
  Boolean multiple
)

{
  PCRPrimerPtr  ppp;

  for (ppp = primers; ppp != NULL; ppp = ppp->next) {
    if (StringDoesHaveText (ppp->seq)) {
      if (StringDoesHaveText (ppp->name)) {
        ValNodeCopyStr (headp, 0, prefix);
        ValNodeCopyStr (headp, 0, nm_label);
        ValNodeCopyStr (headp, 0, ppp->name);
        prefix = ", ";
      }
      ValNodeCopyStr (headp, 0, prefix);
      ValNodeCopyStr (headp, 0, sq_label);
      ValNodeCopyStr (headp, 0, ppp->seq);
      prefix = ", ";
    } else if (name_only_ok) {
      if (StringDoesHaveText (ppp->name)) {
        ValNodeCopyStr (headp, 0, prefix);
        ValNodeCopyStr (headp, 0, nm_label);
        ValNodeCopyStr (headp, 0, ppp->name);
        prefix = ", ";
      }
    }
  }
}

static CharPtr NextPCRReaction (
  PCRReactionPtr prp,
  Boolean isInNote,
  Boolean multiple
)

{
  Boolean       has_fwd_seq = FALSE, has_rev_seq = FALSE;
  ValNodePtr    head = NULL, vnp;
  PCRPrimerPtr  ppp;
  CharPtr       prefix = NULL, str;

  if (prp == NULL) return NULL;

  for (ppp = prp->forward; ppp != NULL; ppp = ppp->next) {
    if (StringDoesHaveText (ppp->seq)) {
      has_fwd_seq = TRUE;
    }
  }

  for (ppp = prp->reverse; ppp != NULL; ppp = ppp->next) {
    if (StringDoesHaveText (ppp->seq)) {
      has_rev_seq = TRUE;
    }
  }

  if (has_fwd_seq && has_rev_seq) {
    if (isInNote) {
      return StringSave ("");
    } else {
      PrintHalfReaction (&head, prp->forward, "fwd_name: ", "fwd_seq: ", NULL, FALSE, multiple);
      PrintHalfReaction (&head, prp->reverse, "rev_name: ", "rev_seq: ", ", ", FALSE, multiple);
    }
  } else {
    if (isInNote) {
      PrintHalfReaction (&head, prp->forward, "fwd_name: ", "fwd_seq: ", NULL, TRUE, multiple);
      if (head != NULL) {
        prefix = ", ";
      }
      PrintHalfReaction (&head, prp->reverse, "rev_name: ", "rev_seq: ", prefix, TRUE, multiple);
    } else {
      return StringSave ("");
    }
  }

  if (head != NULL && isInNote) {
    vnp = ValNodeCopyStr (NULL, 0, "PCR_primers=");
    if (vnp != NULL) {
      vnp->next = head;
      head = vnp;
    }
  }

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);
  return str;
}

/* specimen_voucher, culture_collection, bio_material hyperlinks */

#define s_atcc_base "http://www.atcc.org/SearchCatalogs/linkin?id="
#define s_bcrc_base "http://strain.bcrc.firdi.org.tw/BSAS/controller?event=SEARCH&bcrc_no="
#define s_ccmp_base "http://ccmp.bigelow.org/SD/display.php?strain=CCMP"
#define s_ccug_base "http://www.ccug.se/default.cfm?page=search_record.cfm&db=mc&s_tests=1&ccugno="
#define s_fsu_base  "http://www.prz.uni-jena.de/data.php?fsu="
#define s_ku_base   "http://collections.nhm.ku.edu/"
#define s_mvz_base  "http://mvzarctos.berkeley.edu/SpecimenDetail.cfm?guid="
#define s_pcc_base  "http://www.pasteur.fr/recherche/banques/PCC/docs/pcc"
#define s_pcmb_base "http://www2.bishopmuseum.org/HBS/PCMB/results3.asp?searchterm3="
#define s_uam_base  "http://arctos.database.museum/SpecimenDetail.cfm?GUID="

#define s_colon_pfx ":"

#define s_kui_pfx   "KU_Fish/detail.jsp?record="
#define s_kuit_pfx  "KU_Tissue/detail.jsp?record="

#define s_bcrc_sfx  "&type_id=6&keyword=;;"
#define s_pcc_sfx   ".htm"

#define s_atcc_inst  "American Type Culture Collection"
#define s_bcrc_inst  "Bioresource Collection and Research Center"
#define s_ccmp_inst  "Provasoli-Guillard National Center for Culture of Marine Phytoplankton"
#define s_ccug_inst  "Culture Collection, University of Goteborg, Department of Clinical Bacteriology"
#define s_crcm_inst  "Charles R. Conner Museum, Washington State University"
#define s_dgr_inst   "Division of Genomic Resources, University of New Mexico"
#define s_fsu_inst   "Fungal Reference Center, University of Jena"
#define s_ku_inst    "University of Kansas, Museum of Natural History"
#define s_kwp_inst   "Kenelm W. Philip Collection, University of Alaska Museum of the North"
#define s_msb_inst   "Museum of Southwestern Biology, University of New Mexico"
#define s_mvz_inst   "Museum of Vertebrate Zoology, University of California"
#define s_nbsb_inst  "National Biomonitoring Specimen Bank, U.S. Geological Survey"
#define s_pcc_inst   "Pasteur Culture Collection of Cyanobacteria"
#define s_pcmb_inst  "Pacific Center for Molecular Biodiversity"
#define s_psu_inst   "Portland State University"
#define s_uam_inst   "University of Alaska Museum of the North"
#define s_wmnu_inst  "Western New Mexico University Museum"

typedef struct vouch {
  CharPtr  sites;
  CharPtr  links;
  Boolean  prepend_institute;
  CharPtr  prefix;
  CharPtr  suffix;
  CharPtr  mouseover;
} VouchData, PNTR VouchDataPtr;

static VouchData Nlm_spec_vouchers [] = {
 { "ATCC",        s_atcc_base, FALSE, NULL,         NULL,       s_atcc_inst },
 { "BCRC",        s_bcrc_base, FALSE, NULL,         s_bcrc_sfx, s_bcrc_inst },
 { "CCMP",        s_ccmp_base, FALSE, NULL,         NULL,       s_ccmp_inst },
 { "CCUG",        s_ccug_base, FALSE, NULL,         NULL,       s_ccug_inst },
 { "CRCM:Bird",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_crcm_inst },
 { "DGR:Bird",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_dgr_inst  },
 { "DGR:Ento",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_dgr_inst  },
 { "DGR:Fish",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_dgr_inst  },
 { "DGR:Herp",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_dgr_inst  },
 { "DGR:Mamm",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_dgr_inst  },
 { "FSU<DEU>",    s_fsu_base,  FALSE, NULL,         NULL,       s_fsu_inst  },
 { "KU:I",        s_ku_base,   FALSE, s_kui_pfx,    NULL,       s_ku_inst   },
 { "KU:IT",       s_ku_base,   FALSE, s_kuit_pfx,   NULL,       s_ku_inst   },
 { "KWP:Ento",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_kwp_inst  },
 { "MSB:Bird",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_msb_inst  },
 { "MSB:Mamm",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_msb_inst  },
 { "MVZ:Bird",    s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Egg",     s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Herp",    s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Hild",    s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Img",     s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Mamm",    s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZ:Page",    s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "MVZObs:Herp", s_mvz_base,  TRUE,  s_colon_pfx,  NULL,       s_mvz_inst  },
 { "NBSB:Bird",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_nbsb_inst },
 { "PCC",         s_pcc_base,  FALSE, NULL,         s_pcc_sfx,  s_pcc_inst  },
 { "PCMB",        s_pcmb_base, FALSE, NULL,         NULL,       s_pcmb_inst },
 { "PSU:Mamm",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_psu_inst  },
 { "UAM:Bird",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Bryo",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Crus",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Ento",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Fish",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Herb",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Herp",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Mamm",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Moll",    s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAM:Paleo",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "UAMObs:Mamm", s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_uam_inst  },
 { "WNMU:Bird",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_wmnu_inst },
 { "WNMU:Fish",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_wmnu_inst },
 { "WNMU:Mamm",   s_uam_base,  TRUE,  s_colon_pfx,  NULL,       s_wmnu_inst },
 { NULL,          NULL,        FALSE, NULL,         NULL,       NULL        }
};

static Int2 VoucherNameIsValid (
  CharPtr name
)

{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return -1;
  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ' ');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (Nlm_spec_vouchers) / sizeof (Nlm_spec_vouchers [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_spec_vouchers [mid].sites, str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  /* case sensitive comparison at end enforces strictness */

  if (StringCmp (Nlm_spec_vouchers [R].sites, str) == 0) {
    return R;
  }

  return -1;
}

/* works on subname copy that it can change */

static Boolean ParseSecVoucher (
  CharPtr subname,
  CharPtr PNTR inst,
  CharPtr PNTR id
)

{
  CharPtr  ptr;
  CharPtr  tmp;

  if (StringHasNoText (subname)) return FALSE;
  if (StringLen (subname) < 5) return FALSE;
  TrimSpacesAroundString (subname);

  ptr = StringChr (subname, ':');
  if (ptr == NULL) return FALSE;

  *inst = subname;

  tmp = StringChr (ptr + 1, ':');
  if (tmp != NULL) {
    *tmp = '\0';
    tmp++;
    TrimSpacesAroundString (tmp);
    *id = tmp;
  } else {
    *ptr = '\0';
    ptr++;
    TrimSpacesAroundString (ptr);
    *id = ptr;
  }

  if (StringHasNoText (*inst) || StringHasNoText (*id)) return FALSE;

  return TRUE;
}

static void Do_www_specimen_voucher (
  StringItemPtr ffstring,
  CharPtr inst,
  CharPtr id,
  VouchDataPtr vdp
)

{
  if ( ffstring == NULL || inst == NULL || id == NULL || vdp == NULL || vdp->links == NULL ) return;

  FFAddOneString (ffstring, inst, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, vdp->links, FALSE, FALSE, TILDE_IGNORE);
  if (vdp->prepend_institute) {
    FFAddOneString (ffstring, inst, FALSE, FALSE, TILDE_IGNORE);
  }
  if (vdp->prefix != NULL) {
    FFAddOneString (ffstring, vdp->prefix, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneString (ffstring, id, FALSE, FALSE, TILDE_IGNORE);
  if (vdp->suffix != NULL) {
    FFAddOneString (ffstring, vdp->suffix, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneString(ffstring, "\"", FALSE, FALSE, TILDE_IGNORE);
  if (vdp->mouseover != NULL) {
    FFAddTextToString (ffstring, " title=\"", vdp->mouseover, "\"",
                       FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddTextToString (ffstring, ">", id, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

NLM_EXTERN void FF_www_specimen_voucher (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr subname
)

{
  Char          buf [512];
  CharPtr       inst = NULL, id = NULL;
  Int2          R;
  VouchDataPtr  vdp;

  if ( ffstring == NULL || subname == NULL ) return;
  if (! GetWWW (ajp)) { /* not in www mode */
    FFAddTextToString(ffstring, NULL, subname, NULL, FALSE, TRUE, TILDE_TO_SPACES);
    return;
  }
  StringNCpy_0 (buf, subname, sizeof (buf));
  if (! ParseSecVoucher (buf, &inst, &id)) {
    FFAddTextToString (ffstring, NULL, subname, NULL, FALSE, TRUE, TILDE_TO_SPACES);
    return;
  }
  R = VoucherNameIsValid (inst);
  if (R < 0) {
    FFAddTextToString (ffstring, NULL, subname, NULL, FALSE, TRUE, TILDE_TO_SPACES);
    return;
  }
  vdp = &(Nlm_spec_vouchers [R]);
  if (vdp == NULL || vdp->links == NULL) {
    FFAddTextToString (ffstring, NULL, subname, NULL, FALSE, TRUE, TILDE_TO_SPACES);
    return;
  }
  Do_www_specimen_voucher (ffstring, inst, id, vdp);
}

static void Do_www_lat_lon (
  StringItemPtr ffstring,
  CharPtr lat_lon
)

{
  Char     buf [128];
  Char     ch;
  CharPtr  ew = "";
  Int2     i;
  CharPtr  ns = "";
  CharPtr  ptr;
  Char     tmp [128];
  CharPtr  tokens [6];

  if ( ffstring == NULL || lat_lon == NULL ) return;

  MemSet ((Pointer) tokens, 0, sizeof (tokens));

  StringNCpy_0 (buf, lat_lon, sizeof (buf));

  i = 0;
  ptr = buf;
  ch = *ptr;
  tokens [i] = ptr;
  while (ch != '\0' && i < 5) {
    if (ch == ' ') {
      *ptr = '\0';
      ptr++;
      ch = *ptr;
      while (ch == ' ') {
        ptr++;
        ch = *ptr;
      }
      i++;
      tokens [i] = ptr;
    } else {
      ptr++;
      ch = *ptr;
    }
  }

  ptr = tokens [1];
  if (ptr != NULL && *ptr == 'S') {
    ns = "-";
  }
  ptr = tokens [3];
  if (ptr != NULL && *ptr == 'W') {
    ew = "-";
  }

  if (tokens [0] == NULL) {
    tokens [0] = "?";
  }
  if (tokens [2] == NULL) {
    tokens [2] = "?";
  }

  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
  FF_Add_NCBI_Base_URL (ffstring, link_lat_lon);
  sprintf (tmp, "lat=%s%s&lon=%s%s", ns, tokens [0], ew, tokens [2]);
  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString (ffstring, "\">", lat_lon, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

static void FF_www_lat_lon (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr lat_lon
)

{
  Boolean  format_ok = FALSE;
  FloatHi  lat = 0.0;
  FloatHi  lon = 0.0;
  Boolean  lat_in_range = FALSE;
  Boolean  lon_in_range = FALSE;

  if ( ffstring == NULL || lat_lon == NULL ) return;
  if (! GetWWW (ajp)) { /* not in www mode */
    FFAddTextToString(ffstring, NULL, lat_lon, NULL, FALSE, TRUE, TILDE_TO_SPACES);
    return;
  }
  if (StringDoesHaveText (lat_lon)) {
    IsCorrectLatLonFormat (lat_lon, &format_ok, &lat_in_range, &lon_in_range);
    if (format_ok && lat_in_range && lon_in_range) {
      if (ParseLatLon (lat_lon, &lat, &lon)) {
        Do_www_lat_lon (ffstring, lat_lon);
        return;
      }
    }
  }

  /* if any of above tests failed, default print */
  FFAddTextToString (ffstring, NULL, lat_lon, NULL, FALSE, TRUE, TILDE_TO_SPACES);
}

NLM_EXTERN CharPtr FormatSourceFeatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Boolean            add_period;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Char               buf [80];
  CharPtr            common = NULL;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBFeaturePtr       gbfeat = NULL;
  GBSeqPtr           gbseq;
  Int2               i;
  Uint1              idx;
  IntSrcBlockPtr     isp;
  Boolean            is_desc = TRUE;
  Boolean            is_gps = FALSE;
  Boolean            is_other = FALSE;
  Boolean            is_bc;
  Boolean            is_rf;
  Boolean            is_sc;
  Int2               j;
  Uint1              jdx;
  Uint1              lastomptype;
  Uint1              lastssptype;
  SeqLocPtr          location = NULL;
  CharPtr            notestr;
  SourceType PNTR    notetbl = NULL;
  Boolean            okay;
  ObjectIdPtr        oip;
  OrgModPtr          omp;
  OrgNamePtr         onp = NULL;
  OrgRefPtr          orp = NULL;
  Boolean            partial5;
  Boolean            partial3;
  CharPtr            prefix;
  PCRReactionPtr     prp;
  ValNodePtr         pset;
  PcrSetPtr          psp;
  SourceType PNTR    qualtbl = NULL;
  QualValPtr         qvp;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp = NULL;
  SeqIdPtr           sip;
  SubSourcePtr       ssp;
  CharPtr            str;
  BioseqPtr          target;
  CharPtr            taxname = NULL;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, unique;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  /* five-column feature table uses special code for formatting */

  if (ajp->format == FTABLE_FMT) {
    str = FormatFtableSourceFeatBlock (bbp, target);
    return str;
  }

  /* otherwise do regular flatfile formatting */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  isp = (IntSrcBlockPtr) bbp;

  /* could be descriptor or feature */

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
    is_desc = FALSE;
  }

  if (biop == NULL) return NULL;

  unique = FFGetString(ajp);
  if ( unique == NULL ) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  FFStartPrint (ffstring, afp->format, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
  FFAddOneString (ffstring, "source", FALSE, FALSE, TILDE_IGNORE);
  FFAddNChar(ffstring, ' ', 21 - 5 - StringLen("source"), FALSE);

  if (gbseq != NULL) {
    gbfeat = GBFeatureNew ();
    if (gbfeat != NULL) {
      gbfeat->key = StringSave ("source");
    }
  }

  location = isp->loc;

  str = FFFlatLoc (ajp, bsp, location, ajp->masterStyle);
  if ( GetWWW(ajp) ) {
    FF_www_featloc (ffstring, str);
  } else {
    FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneChar(ffstring, '\n', FALSE);

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      if (StringDoesHaveText (str)) {
        gbfeat->location = StringSave (str);
      } else {
        gbfeat->location = StringSave ("");
      }
      if (StringDoesHaveText (str)) {
        if (StringStr (str, "join") != NULL) {
          gbfeat->operator__ = StringSave ("join");
        } else if (StringStr (str, "order") != NULL) {
          gbfeat->operator__ = StringSave ("order");
        }
      }
      CheckSeqLocForPartial (location, &partial5, &partial3);
      gbfeat->partial5 = partial5;
      gbfeat->partial3 = partial3;
      if (ajp->masterStyle) {
        AddIntervalsToGbfeat (gbfeat, location, bsp);
      } else {
        AddIntervalsToGbfeat (gbfeat, location, NULL);
      }
    }
  }

  MemFree (str);

  orp = biop->org;
  if (orp != NULL) {
    taxname = orp->taxname;
    /* common = orp->common; */
  }
  if (StringHasNoText (taxname)) {
    if (ajp->flags.needOrganismQual) {
      taxname = "unknown";
      common = orp->common;
#ifdef ASN2GNBK_PRINT_UNKNOWN_ORG
    } else {
      taxname = "unknown";
      common = orp->common;
#endif
    }
  }

  sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
  if (sep == NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      is_gps = TRUE;
    }
  }

  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        is_other = TRUE;
      }
    }
  }
 
  if (ajp->refseqConventions) {
    is_other = TRUE;
  }

  /* populate qualifier table from biosource fields */

  qvp [SCQUAL_organism].str = taxname;
  qvp [SCQUAL_common_name].str = common;

  if (biop->is_focus) {
    qvp [SCQUAL_focus].ble = TRUE;
  }

  str = GetMolTypeQual (bsp);
  /*
  if (StringICmp (str, "ncRNA") == 0) {
    str = "other RNA";
  }
  */
  if (str == NULL) {
    switch (bsp->mol) {
      case Seq_mol_dna :
        str = "unassigned DNA";
        break;
      case Seq_mol_rna :
        str = "unassigned RNA";
        break;
      case Seq_mol_aa :
        break;
      default :
        str = "unassigned DNA";
        break;
    }
  }
  qvp [SCQUAL_mol_type].str = str;

  SubSourceToQualArray (biop->subtype, qvp);

  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      OrgModToQualArray (onp->mod, qvp);
    }

    if (! is_desc) {
      qvp [SCQUAL_unstructured].vnp = orp->mod;
    }
    qvp [SCQUAL_db_xref].vnp = orp->db;
  }

  if (sfp != NULL) {
    qvp [SCQUAL_org_xref].vnp = sfp->dbxref;
  }

  /* organelle currently prints /mitochondrion, /virion, etc. */

  qvp [SCQUAL_organelle].num = biop->genome;

  /* some qualifiers are flags in genome and names in subsource, print once with name */

  if (qvp [SCQUAL_ins_seq_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_insertion_seq) {
    qvp [SCQUAL_organelle].num = 0;
  }
  if (qvp [SCQUAL_plasmid_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_plasmid) {
    qvp [SCQUAL_organelle].num = 0;
  }
  /* AF095904.1
  if (qvp [SCQUAL_plastid_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_plastid) {
    qvp [SCQUAL_organelle].num = 0;
  }
  */
  if (qvp [SCQUAL_transposon_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_transposon) {
    qvp [SCQUAL_organelle].num = 0;
  }

  if (sfp != NULL) {
    qvp [SCQUAL_seqfeat_note].str = sfp->comment;
  }

  if (qvp [SCQUAL_fwd_primer_name].ssp != NULL ||
      qvp [SCQUAL_fwd_primer_seq].ssp != NULL ||
      qvp [SCQUAL_rev_primer_name].ssp != NULL ||
      qvp [SCQUAL_rev_primer_seq].ssp != NULL) {
    qvp [SCQUAL_PCR_primers].ble = TRUE;
    qvp [SCQUAL_PCR_primer_note].ble = TRUE;
  }

  if (biop->pcr_primers != NULL) {
    qvp [SCQUAL_PCR_reaction].prp = biop->pcr_primers;
  }

  if (is_other || (ajp->mode == SEQUIN_MODE || ajp->mode == DUMP_MODE)) {
    /* leave metagenome_source as a separate qualifier */
  } else {
    /* move metagenome_source to note */
    qvp [SCQUAL_metagenome_note].omp = qvp [SCQUAL_metagenome_source].omp;
    qvp [SCQUAL_metagenome_source].omp = NULL;
  }

#if 0
  if (is_other || (ajp->mode == SEQUIN_MODE || ajp->mode == DUMP_MODE)) {
    /* leave mating_type as a separate qualifier */
  } else if (qvp [SCQUAL_sex].ssp == NULL &&  qvp [SCQUAL_mating_type].ssp != NULL) {
    /* move mating_type to sex if available */
    qvp [SCQUAL_sex].ssp = qvp [SCQUAL_mating_type].ssp;
    qvp [SCQUAL_mating_type].ssp = NULL;
  }
#endif

  /* now print qualifiers from table */

  qualtbl = source_qual_order;
  if (is_desc) {
    notetbl = source_desc_note_order;
  } else {
    notetbl = source_feat_note_order;
  }

  for (i = 0, idx = qualtbl [i]; idx != 0; i++, idx = qualtbl [i]) {

    lastomptype = 0;
    lastssptype = 0;
    switch (asn2gnbk_source_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=", 
                            FALSE, FALSE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"", 
                            FALSE, FALSE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].ble) {
          FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "\n",
                            FALSE, TRUE, TILDE_IGNORE);
        }
        break;

      case Qual_class_organelle :
        j = (Int2) qvp [idx].num;
        if (j < sizeof (organelleQual) / sizeof (CharPtr)) {
          if (organelleQual [j] != NULL) {
            FFAddTextToString(ffstring, NULL, organelleQual[j], "\n",
                              FALSE, FALSE, TILDE_IGNORE);
          }
        }
        break;

      case Qual_class_orgmod :
        omp = qvp [idx].omp;
        if (lastomptype == 0 && omp != NULL) {
          lastomptype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lastomptype) {
          if (StringIsJustQuotes (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, "\"", omp->subname, "\"\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          }
          omp = omp->next;
        }
        break;

      case Qual_class_voucher :
        omp = qvp [idx].omp;
        if (lastomptype == 0 && omp != NULL) {
          lastomptype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lastomptype) {
          if (StringIsJustQuotes (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"",
                              FALSE, TRUE, TILDE_IGNORE);
            FF_www_specimen_voucher(ajp, ffstring, omp->subname);
            FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
          }
          omp = omp->next;
        }
        break;

      case Qual_class_lat_lon :
        omp = qvp [idx].omp;
        if (lastomptype == 0 && omp != NULL) {
          lastomptype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lastomptype) {
          if (StringIsJustQuotes (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"",
                              FALSE, TRUE, TILDE_IGNORE);
            FF_www_lat_lon(ajp, ffstring, omp->subname);
            FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
          }
          omp = omp->next;
        }
        break;

      case Qual_class_subsource :
        ssp = qvp [idx].ssp;
        if (lastssptype == 0 && ssp != NULL) {
          lastssptype = ssp->subtype;
        }
        while (ssp != NULL && ssp->subtype == lastssptype) {
          if (ssp->subtype == SUBSRC_germline ||
              ssp->subtype == SUBSRC_rearranged ||
              ssp->subtype == SUBSRC_transgenic ||
              ssp->subtype == SUBSRC_environmental_sample ||
              ssp->subtype == SUBSRC_metagenomic) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          } else if (StringIsJustQuotes (ssp->name)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (ssp->name)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, "\"", ssp->name, "\"\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          }
          ssp = ssp->next;
        }
        break;

      case Qual_class_pcr :
        if (qvp [idx].ble) {
          lastssptype = 0;
          pset = ParsePCRPrimerString (qvp);
          for (vnp = pset; vnp != NULL; vnp = vnp->next) {
            psp = (PcrSetPtr) vnp->data.ptrvalue;
            if (psp == NULL) continue;
            str = NextPCRPrimerString (psp, FALSE, (Boolean) (pset->next != NULL));
            if (str == NULL) continue;
            if (! StringHasNoText (str)) {
              FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                                FALSE, TRUE, TILDE_IGNORE);
              FFAddTextToString(ffstring, "\"", str, "\"\n",
                                FALSE, TRUE, TILDE_TO_SPACES);
            }
            MemFree (str);
          }
          FreePCRSet (pset);
        }
        break;

      case Qual_class_pcr_react :
        prp = qvp [idx].prp;
        while (prp != NULL) {
          str = NextPCRReaction (prp, FALSE, (Boolean) (prp->next != NULL));
          if (StringDoesHaveText (str)) {
            FFAddTextToString (ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                               FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString (ffstring, "\"", str, "\"\n",
                               FALSE, TRUE, TILDE_TO_SPACES);
          }
          MemFree (str);
          prp = prp->next;
        }
        break;

      case Qual_class_pubset :
        break;

      case Qual_class_quote :
        break;

      case Qual_class_noquote :
        break;

      case Qual_class_label :
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {

              okay = TRUE;
              if (ajp->flags.dropBadDbxref) {
                /* if RELEASE_MODE, drop unknown dbtag */

                okay = FALSE;
                if (DbxrefIsValid (dbt->db, &is_rf, &is_sc, &is_bc, NULL)) {
                  if (is_bc) {
                    /* case counts, so suppress if bad case */
                  } else if (is_rf && (is_other || is_gps)) {
                    /* allow refseq dbxrefs in source feature */
                    okay = TRUE;
                  } else if (is_sc) {
                    /* expect it to be in legalSrcDbXrefs list */
                    okay = TRUE;
                  } else {
                    /* suppress regular dbxrefs, also warn in validator */
                  }
                }

                /*
                okay = FALSE;
                for (j = 0; legalDbXrefs [j] != NULL; j++) {
                  if (StringCmp (dbt->db, legalDbXrefs [j]) == 0) {
                    okay = TRUE;
                  }
                }
                */
              }

              if (okay) {
                if (! StringHasNoText (oip->str)) {
                  if (StringLen (dbt->db) + StringLen (oip->str) < 80) {
                    sprintf (buf, "%s", oip->str);
                  }
                } else {
                  sprintf (buf, "%ld", (long) oip->id);
                }
              }
            }
          }
          if (! StringHasNoText (buf)) {
            FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
            FF_www_db_xref(ajp, ffstring, dbt->db, buf, bsp);
            FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
          }
        }
        break;

      case Qual_class_illegal :
        break;

      case Qual_class_note :
        if (! ajp->flags.srcQualsToNote) {

          /* in sequin_mode and dump_mode, all orgmods and subsources show up as separate /qualifiers */

          for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {

            lastomptype = 0;
            lastssptype = 0;
            switch (asn2gnbk_source_quals [jdx].qualclass) {

              case Qual_class_orgmod :
                if (jdx == SCQUAL_orgmod_note) break;
                omp = qvp [jdx].omp;
                if (lastomptype == 0 && omp != NULL) {
                  lastomptype = omp->subtype;
                }
                while (omp != NULL && omp->subtype == lastomptype) {
                  if (StringIsJustQuotes (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
                  } else if (! StringHasNoText (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", omp->subname, "\"\n",
                                      FALSE, TRUE, TILDE_TO_SPACES);
                  }
                  omp = omp->next;
                }
                break;

              case Qual_class_voucher :
                if (jdx == SCQUAL_orgmod_note) break;
                omp = qvp [jdx].omp;
                if (lastomptype == 0 && omp != NULL) {
                  lastomptype = omp->subtype;
                }
                while (omp != NULL && omp->subtype == lastomptype) {
                  if (StringIsJustQuotes (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
                  } else if (! StringHasNoText (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FF_www_specimen_voucher(ajp, ffstring, omp->subname);
                    FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                  }
                  omp = omp->next;
                }
                break;

              case Qual_class_subsource :
                if (jdx == SCQUAL_subsource_note) break;
                ssp = qvp [jdx].ssp;
                if (lastssptype == 0 && ssp != NULL) {
                  lastssptype = ssp->subtype;
                }
                while (ssp != NULL && ssp->subtype == lastssptype) {
                  if (ssp->subtype == SUBSRC_germline ||
                      ssp->subtype == SUBSRC_rearranged ||
                      ssp->subtype == SUBSRC_transgenic ||
                      ssp->subtype == SUBSRC_environmental_sample ||
                      ssp->subtype == SUBSRC_metagenomic) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "\n",
                                      FALSE, TRUE, TILDE_TO_SPACES);
                  } else if (StringIsJustQuotes (ssp->name)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"\"\n",
                                      FALSE, TRUE, TILDE_IGNORE);

                  } else if (! StringHasNoText (ssp->name)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", ssp->name, "\"\n",
                                      FALSE, TRUE, TILDE_TO_SPACES);
                  }
                  ssp = ssp->next;
                }
                break;

              default :
                break;
            }
          }
        }

        notestr = NULL;
        prefix = "";
        add_period = FALSE;

        if (biop->genome == 8) {
          FFAddTextToString(unique, "", "extrachromosomal", NULL, FALSE, FALSE, TILDE_IGNORE); 
          prefix = "\n";
        }

        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {

          lastomptype = 0;
          lastssptype = 0;
          switch (asn2gnbk_source_quals [jdx].qualclass) {

            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL, FALSE);
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_orgmod :
            case Qual_class_voucher :
              if ((! ajp->flags.srcQualsToNote) && jdx != SCQUAL_orgmod_note) break;
              omp = qvp [jdx].omp;
              if (lastomptype == 0 && omp != NULL) {
                lastomptype = omp->subtype;
              }
              while (omp != NULL && omp->subtype == lastomptype) {
                if (! StringHasNoText (omp->subname)) {
                  if (jdx == SCQUAL_orgmod_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }

                  str = StringSave (omp->subname);
                  add_period = s_RemovePeriodFromEnd (str);
                  if (jdx == SCQUAL_orgmod_note) {
                    FFAddString_NoRedund (unique, buf, str, NULL, FALSE);
                  } else {
                    FFAddTextToString(unique, buf, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);

                  if (jdx == SCQUAL_orgmod_note) {
                    if (add_period) {
                      prefix = ".\n";
                    } else {
                      prefix = "\n";
                    }
                  } else {
                    prefix = "; ";
                  }
                }
                omp = omp->next;
              }
              break;

            case Qual_class_subsource :
              if ((! ajp->flags.srcQualsToNote) && jdx != SCQUAL_subsource_note) break;
              ssp = qvp [jdx].ssp;
              if (lastssptype == 0 && ssp != NULL) {
                lastssptype = ssp->subtype;
              }
              while (ssp != NULL && ssp->subtype == lastssptype) {
                if (ssp->subtype == SUBSRC_germline ||
                    ssp->subtype == SUBSRC_rearranged ||
                    ssp->subtype == SUBSRC_transgenic ||
                    ssp->subtype == SUBSRC_environmental_sample ||
                    ssp->subtype == SUBSRC_metagenomic) {
                  FFAddTextToString (unique, prefix, asn2gnbk_source_quals [jdx].name, NULL, FALSE, FALSE, TILDE_IGNORE);
                  prefix = "; ";
                } else if (! StringHasNoText (ssp->name)) {
                  if (jdx == SCQUAL_subsource_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }

                  str = StringSave (ssp->name);
                  add_period = s_RemovePeriodFromEnd (str);
                  if (jdx == SCQUAL_subsource_note) {
                    FFAddString_NoRedund (unique, buf, str, NULL, FALSE);
                  } else {
                    FFAddTextToString(unique, buf, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);

                  if (jdx == SCQUAL_subsource_note) {
                    if (add_period) {
                      prefix = ".\n";
                    } else {
                      prefix = "\n";
                    }
                  } else {
                    prefix = "; ";
                 }
                }
                ssp = ssp->next;
              }
              break;

            case Qual_class_pcr :
              if (qvp [jdx].ble) {
                lastssptype = 0;
                pset = ParsePCRPrimerString (qvp);
                for (vnp = pset; vnp != NULL; vnp = vnp->next) {
                  psp = (PcrSetPtr) vnp->data.ptrvalue;
                  if (psp == NULL) continue;
                  str = NextPCRPrimerString (psp, TRUE, (Boolean) (pset->next != NULL));
                  if (str == NULL) continue;
                  if (! StringHasNoText (str)) {
                    FFAddString_NoRedund (unique, prefix, str, NULL, FALSE);
                    add_period = FALSE;
                    prefix = "; ";
                  }
                  MemFree (str);
                }
                FreePCRSet (pset);
              }
              break;

            case Qual_class_pcr_react :
              prp = qvp [jdx].prp;
              while (prp != NULL) {
                str = NextPCRReaction (prp, TRUE, (Boolean) (prp->next != NULL));
                if (StringDoesHaveText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, FALSE);
                  add_period = FALSE;
                  prefix = "; ";
                }
                MemFree (str);
                prp = prp->next;
              }
              break;

            case Qual_class_valnode :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL, FALSE);
                  add_period = FALSE;
                  prefix = "; ";
                }
              }
              break;

            default :
              break;
          }
        }
        if ( !FFEmpty(unique) ) {
          notestr = FFToCharPtr(unique);
        
          if (add_period) {
            s_AddPeriodToEnd (notestr);
          }

#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
          if (! IsEllipsis (notestr))
            s_RemovePeriodFromEnd (notestr);
#endif

          FFAddOneString (ffstring, "/note=\"", FALSE, FALSE, TILDE_IGNORE);
          if (is_desc) {
            /* AB055064.1 said TILDE_IGNORE on descriptors, but now changing policy */
            FFAddOneString (ffstring, notestr, FALSE, TRUE, /* TILDE_IGNORE */ /* TILDE_EXPAND */ TILDE_SEMICOLON);
          } else {
            /* ASZ93724.1 said TILDE_EXPAND on features, but record does not exist */
            FFAddOneString (ffstring, notestr, FALSE, TRUE, /* TILDE_EXPAND */ TILDE_SEMICOLON);
          }
          FFAddOneString (ffstring, "\"", FALSE, FALSE, TILDE_IGNORE);

          MemFree (notestr);
        }
        break;
      default :
        break;
    }
  }

  /* and then deal with the various note types separately (not in order table) */

  str = FFEndPrint(ajp, ffstring, afp->format, 21, 21, 5, 21, "FT"); 

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      AddFeatureToGbseq (gbseq, gbfeat, str, NULL);
    }
  }

  FFRecycleString(ajp, unique);
  FFRecycleString(ajp, ffstring);
  return str;
}

static void LIBCALLBACK CountBasesByStream (
  CharPtr sequence,
  Pointer userdata
)

{
  Int4Ptr  base_count;
  Char     ch;
  CharPtr  ptr;

  base_count = (Int4Ptr) userdata;

  ptr = sequence;
  ch = *ptr;
  while (ch != '\0') {
    ch = TO_UPPER (ch);
    switch (ch) {
      case 'A' :
        (base_count [0])++;
        break;
      case 'C' :
        (base_count [1])++;
        break;
      case 'G' :
        (base_count [2])++;
        break;
      case 'T' :
        (base_count [3])++;
        break;
      default :
        (base_count [4])++;
        break;
    }
    ptr++;
    ch = *ptr;
  }
}

NLM_EXTERN CharPtr FormatBasecountBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  Int4             base_count [5];
  BioseqPtr        bsp;
  Char             buf [80];
  Int2             i;
  Int4             len;
  StringItemPtr    ffstring;
  CharPtr          str;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;

  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  /* after first formatting, result is cached into bbp->string */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  for (i = 0; i < 5; i++) {
    base_count [i] = 0;
  }

  if (ajp->ajp.slp != NULL) {
    len = SeqLocLen (ajp->ajp.slp);
    SeqPortStreamLoc (ajp->ajp.slp, STREAM_EXPAND_GAPS, (Pointer) base_count, CountBasesByStream);
  } else {
    len = bsp->length;
    SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) base_count, CountBasesByStream);
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    if (base_count [4] == 0) {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t",
               (long) base_count [0], (long) base_count [1],
               (long) base_count [2], (long) base_count [3]);
    } else {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t%7ld others",
               (long) base_count [0], (long) base_count [1],
               (long) base_count [2], (long) base_count [3],
               (long) base_count [4]);
    }

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    sprintf (buf, "Sequence %ld BP; %ld A; %ld C; %ld G; %ld T; %ld other;",
             (long) len,
             (long) base_count [0], (long) base_count [1],
             (long) base_count [2], (long) base_count [3],
             (long) base_count [4]);
  }

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFAddOneString(ffstring, "XX\n", FALSE, FALSE, TILDE_IGNORE);
  }
  FFStartPrint (ffstring, afp->format, 0, 0, "BASE COUNT", 12, 5, 5, "SQ", FALSE);
  FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
  str = FFEndPrint(ajp, ffstring, afp->format, 12, 0, 5, 5, "SQ");
  FFRecycleString(ajp, ffstring);

  return str;
}

static void PrintSeqLine (
  StringItemPtr ffstring,
  FmtType format,
  CharPtr buf,
  Int4 start,
  Int4 stop
)

{
  size_t  len;
  Char    pos [16];
  Int4    pad;

  len = StringLen (buf);
  if (len > 0 && buf [len - 1] == ' ') {
    buf [len - 1] = '\0';
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {

    sprintf (pos, "%9ld", (long) (start + 1));
    FFAddOneString(ffstring, pos, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, ' ', FALSE);
    FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, '\n', FALSE);
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {

    sprintf (pos, "%8ld", (long) (stop));
    FFAddNChar(ffstring, ' ', 5, FALSE);
    FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
    pad = 72 - 5 - StringLen(buf);
    FFAddNChar(ffstring, ' ', pad, FALSE);
    FFAddOneString(ffstring, pos, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, '\n', FALSE);
  }
}

static CharPtr CompressNonBases (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_ALPHA (ch)) {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

  static Uint1 fasta_order [NUM_SEQID] = {
    33, /* 0 = not set */
    20, /* 1 = local Object-id */
    15, /* 2 = gibbsq */
    16, /* 3 = gibbmt */
    30, /* 4 = giim Giimport-id */
    10, /* 5 = genbank */
    10, /* 6 = embl */
    10, /* 7 = pir */
    10, /* 8 = swissprot */
    15, /* 9 = patent */
    20, /* 10 = other TextSeqId */
    20, /* 11 = general Dbtag */
    255, /* 12 = gi */
    10, /* 13 = ddbj */
    10, /* 14 = prf */
    12, /* 15 = pdb */
    10, /* 16 = tpg */
    10, /* 17 = tpe */
    10, /* 18 = tpd */
    10, /* 19 = gpp */
    10  /* 20 = nat */
  };

static void PrintGenome (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqLocPtr slp_head, 
  CharPtr prefix, 
  Boolean segWithParts,
  Boolean is_na
)
{
  Char         buf[40], val[166];
  Boolean      first = TRUE;
  SeqIdPtr     freeid, sid, newid;
  SeqLocPtr    slp;
  Int4         from, to, start, stop;
  BioseqPtr    bsp = NULL;
  Int2         p1=0, p2=0;

  for (slp = slp_head; slp; slp = slp->next) {
    from = to = 0;
    sid = SeqLocId(slp);
    if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
      start = from = SeqLocStart(slp);
      stop = to = SeqLocStop(slp);
    } else if (slp->choice == SEQLOC_NULL){
      sprintf(val, ",%s", "gap()");
      FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
      continue;
    } else {
      continue;
    }
    if (sid == NULL) {
      continue;
    }
    newid = NULL;
    freeid = NULL;
    buf [0] = '\0';
    if (sid->choice == SEQID_GI) {
      if (GetAccnVerFromServer (sid->data.intvalue, buf)) {
        /* no need to call GetSeqIdForGI */
      } else {
        newid = GetSeqIdForGI (sid->data.intvalue);
        if (newid != NULL) {
          freeid = newid;
        }
        if (newid != NULL && segWithParts) {
          if (newid->choice == SEQID_GIBBSQ ||
              newid->choice == SEQID_GIBBMT ||
              newid->choice == SEQID_GIIM) {
            bsp = BioseqFind (newid);
            if (bsp != NULL && bsp->repr == Seq_repr_virtual) {
              if (bsp->length > 0) {
                sprintf (val, ",gap(%ld)", (long) bsp->length);
              } else {
                sprintf(val, ",%s", "gap()");
              }
              FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
              continue;
            }
          }
        }
      }
    } else if (sid->choice == SEQID_GENERAL) {
      newid = sid;
    } else {
      newid = sid;
    }
    if (prefix != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    }
    if (first) {
      first = FALSE;
    } else {
      FFAddOneChar(ffstring, ',', FALSE);
      /*ff_AddChar(',');*/
    }
    if (! StringHasNoText (buf)) {
      /* filled in by GetAccnVerFromServer */
    } else if (newid != NULL) {
      SeqIdWrite (SeqIdSelect (newid, fasta_order, NUM_SEQID),
                 buf, PRINTID_TEXTID_ACC_VER, sizeof(buf) -1 );
    } else if (sid->choice == SEQID_GI) {
      SeqIdWrite (sid, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    }

    if (SeqLocStrand (slp) == Seq_strand_minus) {
      FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
    }
    if ( GetWWW(ajp) ) {
      if (newid == NULL) {
        newid = sid;
      }
      if (newid->choice != SEQID_GENERAL) {
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        if (is_na) {
          FF_Add_NCBI_Base_URL (ffstring, link_seqn);
        } else {
          FF_Add_NCBI_Base_URL (ffstring, link_seqp);
        }
        FFAddTextToString(ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, NULL, buf, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }
    } else {
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
    }

    if (SeqLocStrand(slp) == Seq_strand_minus) {
      sprintf (val,":%ld..%ld)", (long) start+1, (long) stop+1);
    } else {
      sprintf (val,":%ld..%ld", (long) start+1, (long) stop+1);
    }
    FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
    p1 += StringLen (val);
    p2 += StringLen (val);
    if (freeid != NULL) {
      freeid = SeqIdFree (freeid);
    }
  }
}

NLM_EXTERN CharPtr FormatContigBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BioseqPtr        bsp;
  DeltaSeqPtr      dsp;
  IntFuzzPtr       fuzz;
  GBSeqPtr         gbseq;
  Boolean          is_na;
  SeqLitPtr        litp;
  CharPtr          prefix = NULL;
  Boolean          segWithParts = FALSE;
  SeqLocPtr        slp_head = NULL;
  CharPtr          str;
  Char             tmp [16];
  Boolean          unknown;
  Char             val [32];
  StringItemPtr    ffstring;
/*  CharPtr          label;*/

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  is_na = ISA_na (bsp->mol);

  FFStartPrint (ffstring, afp->format, 0, 0, "CONTIG", 12, 5, 5, "CO", FALSE);
  /*
  if ( GetWWW(ajp) ) {
    label = "CONTIG   ";
  } else {
    label = "CONTIG";
  }
  
  FFAddOneString(ffstring, label,  FALSE, FALSE, TILDE_IGNORE);  
  FFAddNChar(ffstring, ' ', 12 - StringLen(label), FALSE);
  */

  FFAddOneString(ffstring, "join(", FALSE, FALSE, TILDE_IGNORE);

  if (bsp->seq_ext_type == 1) {

    if (bsp->repr == Seq_repr_seg && SegHasParts (bsp)) {
      segWithParts = TRUE;
    }

    slp_head = (SeqLocPtr) bsp->seq_ext;
    PrintGenome (ajp, ffstring, slp_head, prefix, segWithParts, is_na);

  } else if (bsp->seq_ext_type == 4) {

    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp=dsp->next) {
      if (dsp->choice == 1) {

        slp_head = (SeqLocPtr) dsp->data.ptrvalue;
        PrintGenome (ajp, ffstring, slp_head, prefix, FALSE, is_na);

      } else {

        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          if (litp->seq_data != NULL && litp->seq_data_type != Seq_code_gap) {
            if (litp->length == 0) {
              sprintf (val, "gap(%ld)", (long) litp->length);
              FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
            } else {
              /* don't know what to do here */
            }
          } else {
            unknown = FALSE;
            fuzz = litp->fuzz;
            if (fuzz != NULL && fuzz->choice == 4 && fuzz->a == 0) {
              unknown = TRUE;
            }
            if (unknown && litp->length > 0) {
              sprintf (tmp, "unk%ld", (long) litp->length);
            } else {
              sprintf (tmp, "%ld", (long) litp->length);
            }
            if (prefix != NULL) {
              sprintf (val, "%sgap(%s)", prefix, tmp);
            } else {
              sprintf (val, "gap(%s)", tmp);
            }
            FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
          }
        }
      }

      prefix = ",";
    }
  }

  FFAddOneChar(ffstring, ')', FALSE);

  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 5, 5, "CO");
  FFRecycleString(ajp, ffstring);

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    if (StringLen (str) > 12) {
      gbseq->contig = StringSave (str + 12);
    } else {
      gbseq->contig = StringSave (str);
    }

    CleanQualValue (gbseq->contig);
    Asn2gnbkCompressSpaces (gbseq->contig);
    StripAllSpaces (gbseq->contig);
  }

  return str;
}

static void LIBCALLBACK SaveGBSeqSequence (
  CharPtr sequence,
  Pointer userdata
)

{
  CharPtr       tmp;
  CharPtr PNTR  tmpp;

  tmpp = (CharPtr PNTR) userdata;
  tmp = *tmpp;

  tmp = StringMove (tmp, sequence);

  *tmpp = tmp;
}

static Boolean InGapBlock (
  IntAsn2gbJobPtr ajp
)

{
  return (Boolean) (ajp->seqGapCurrLen > 0);
}

static Boolean LineIsAllGaps (
  CharPtr ptr
)

{
  Char  ch;
  Int2  j;

  for (ch = *ptr, j = 0; ch != '\0' && j < 60; ptr++, ch = *ptr, j++) {
    if (ch != '-') return FALSE;
  }
  if (j == 60) return TRUE;
  return FALSE;
}

static Int2 GapAtStart (
  CharPtr ptr
)

{
  Char  ch;
  Int2  j;

  for (ch = *ptr, j = 0; ch != '\0' && j < 60; ptr++, ch = *ptr, j++) {
    if (ch != '-') return j;
  }
  return 0;
}

static void FixGapAtStart (
  CharPtr ptr,
  Char pad
)

{
  Char  ch;
  Int2  j;

  for (ch = *ptr, j = 0; ch == '-' && j < 60; ptr++, ch = *ptr, j++) {
    *ptr = pad;
  }
}

static Int2 GapAtEnd (
  CharPtr ptr
)

{
  Char  ch;
  Int2  j;
  Int2  k;

  for (ch = *ptr, j = 0, k = 0; ch != '\0' && j < 60; ptr++, ch = *ptr, j++) {
    if (ch == '-') {
      k++;
    } else {
      k = 0;
    }
  }
  return k;
}

static void FixGapAtEnd (
  CharPtr ptr,
  Char pad
)

{
  Char  ch;
  Int2  j;

  j = StringLen (ptr) - GapAtEnd (ptr);
  ptr += j;
  for (ch = *ptr; ch == '-' && j < 60; ptr++, ch = *ptr, j++) {
    *ptr = pad;
  }
}

static void FixRemainingGaps (
  CharPtr ptr,
  Char pad
)

{
  Char  ch;
  Int2  j;

  for (ch = *ptr, j = 0; ch != '\0' && j < 60; ptr++, ch = *ptr, j++) {
    if (ch == '-') {
      *ptr = pad;
    }
  }
}

static void ExpandSeqLine (
  CharPtr buf
)

{
  Char     ch;
  Int2     blk, count, lin;
  CharPtr  ptr;
  Char     seq [80];

  StringCpy (seq, buf);

  count = 0;
  blk = 0;
  lin = 0;

  ptr = seq;
  ch = *ptr;

  while (ch != '\0') {
    buf [count] = ch;
    count++;
    ptr++;
    ch = *ptr;

    blk++;
    lin++;
    if (blk >= 10 && lin < 60) {

      buf [count] = ' ';
      count++;
      blk = 0;

    }
  }

  buf [count] = '\0';
}

static Int2 ProcessGapSpecialFormat (
  Asn2gbFormatPtr afp,
  IntAsn2gbJobPtr ajp,
  BioseqPtr bsp,
  StringItemPtr ffstring,
  CharPtr buf,
  CharPtr nextchars
)

{
  Char      fmt_buf [64];
  Char      gapbuf [80];
  Int4      gi;
  Char      gi_buf [16];
  Boolean   is_na;
  Char      pad;
  SeqIdPtr  sip;
  Int2      startgapgap = 0, endgap = 0;

  is_na = ISA_na (bsp->mol);
  if (is_na) {
    pad = 'n';
  } else {
    pad = 'x';
  }

  if (LineIsAllGaps (buf)) {
    ajp->seqGapCurrLen += StringLen (buf);
    *buf = '\0';
    return 0;
  }

  startgapgap = GapAtStart (buf);
  if (InGapBlock (ajp)) {
    ajp->seqGapCurrLen += startgapgap;
    if (is_na) {
      sprintf (gapbuf, "          [gap %ld bp]", (long) ajp->seqGapCurrLen);
    } else {
      sprintf (gapbuf, "          [gap %ld aa]", (long) ajp->seqGapCurrLen);
    }
    FFAddOneString (ffstring, gapbuf, FALSE, FALSE, TILDE_TO_SPACES);
    if (GetWWW (ajp) && ajp->mode == ENTREZ_MODE && afp != NULL &&
      (ajp->format == GENBANK_FMT || ajp->format == GENPEPT_FMT)) {
      gi = 0;
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GI) {
          gi = (Int4) sip->data.intvalue;
        }
      }
      if (gi > 0) {
        sprintf(gi_buf, "%ld", (long) gi);
        sprintf(fmt_buf, "?fmt_mask=%ld", (long) EXPANDED_GAP_DISPLAY);
        if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
          StringCat (fmt_buf, "&report=gbwithparts");
        }
        FFAddOneString (ffstring, "    <a href=\"", FALSE, FALSE, TILDE_IGNORE);
        if (is_na) {
          FF_Add_NCBI_Base_URL (ffstring, link_featn);
        } else {
          FF_Add_NCBI_Base_URL (ffstring, link_featp);
        }
        FFAddOneString (ffstring, gi_buf, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, fmt_buf, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">Expand Ns", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }
    }
    FFAddOneChar (ffstring, '\n', FALSE);
    ajp->seqGapCurrLen = 0;
    FixGapAtStart (buf, ' ');
  } else if (startgapgap > 0) {
    FixGapAtStart (buf, pad);
  }

  endgap = GapAtEnd (buf);
  if (LineIsAllGaps (nextchars)) {
    FixGapAtEnd (buf, ' ');
    ajp->seqGapCurrLen += endgap;
  } else if (endgap > 0) {
    /*
    FixGapAtEnd (buf, pad);
    */
    FixGapAtEnd (buf, ' ');
    ajp->seqGapCurrLen += endgap;
  }

  FixRemainingGaps (buf, pad);

  return startgapgap;
}

/*
static void ChangeOandJtoX (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    if (ch == 'O' || ch == 'J') {
      *str = 'X';
    } else if (ch == 'o' || ch == 'j') {
      *str = 'x';
    }
    str++;
    ch = *str;
  }
}
*/

NLM_EXTERN CharPtr FormatSequenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr   ajp;
  Asn2gbSectPtr     asp;
  Int2              blk;
  BioseqPtr         bsp;
  Bioseq            bsq;
  Char              buf [80];
  Char              ch;
  Int2              count;
  Int4              extend;
  StreamFlgType     flags = STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL;
  GBSeqPtr          gbseq;
  IntAsn2gbSectPtr  iasp;
  Int2              lin;
  SeqLocPtr         loc;
  Int4              num;
  CharPtr           ptr;
  Int4              remaining;
  SeqBlockPtr       sbp;
  SeqLoc            sl;
  SeqLocPtr         slp;
  Int4              start;
  Int2              startgapgap;
  Int4              stop;
  CharPtr           str = NULL;
  CharPtr           tmp;
  StringItemPtr     ffstring;

  if (afp == NULL || bbp == NULL) return NULL;
  sbp = (SeqBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  iasp = (IntAsn2gbSectPtr) asp;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  /* if GBSeq XML, use SeqPortStream on single block */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;

    if (ajp->ajp.slp != NULL) {
      slp = ajp->ajp.slp;
      str = MemNew (sizeof (Char) * (SeqLocLen (slp) + 10));
    } else {
      str = MemNew (sizeof (Char) * (bsp->length + 10));
    }
    if (str == NULL) return NULL;

    tmp = str;
    if (ajp->ajp.slp != NULL) {
      slp = ajp->ajp.slp;
      SeqPortStreamLoc (slp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) &tmp, SaveGBSeqSequence);
    } else {
      SeqPortStream (bsp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) &tmp, SaveGBSeqSequence);
    }
    /*
    if (ISA_aa (bsp->mol) && StringDoesHaveText (str)) {
      if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
        ChangeOandJtoX (str);
      }
    }
    */
    gbseq->sequence = StringSave (str);

    tmp = gbseq->sequence;
    if (tmp == NULL) return NULL;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == '\n' || ch == '\r' || ch == '\t') {
        *tmp = ' ';
      } else if (IS_UPPER (ch)) {
        /* collab decision to present target sequence in lower case */
        *tmp = TO_LOWER (ch);
      }
      tmp++;
      ch = *tmp;
    }
    TrimSpacesAroundString (gbseq->sequence);
    CompressNonBases (gbseq->sequence);

    return str;
  }

  /* replace SeqPort with improved SeqPortStream */

  if (sbp->bases == NULL) {
    if (ajp->specialGapFormat) {
      flags = EXPAND_GAPS_TO_DASHES | STREAM_CORRECT_INVAL;
    }

    start = sbp->start;
    stop = sbp->stop;
    extend = sbp->extend;

    if (stop > start) {

      str = MemNew (sizeof (Char) * (extend - start + 3));
      if (str != NULL) {
        if (ajp->ajp.slp != NULL) {
          slp = ajp->ajp.slp;
          MemSet ((Pointer) &bsq, 0, sizeof (Bioseq));
          MemSet ((Pointer) &sl, 0, sizeof (SeqLoc));
          bsq.repr = Seq_repr_seg;
          bsq.mol = bsp->mol;
          bsq.seq_ext_type = 1;
          bsq.length = SeqLocLen (slp);
          bsq.seq_ext = &sl;
          if (slp->choice == SEQLOC_MIX || slp->choice == SEQLOC_PACKED_INT) {
            loc = (SeqLocPtr) slp->data.ptrvalue;
            if (loc != NULL) {
              sl.choice = loc->choice;
              sl.data.ptrvalue = (Pointer) loc->data.ptrvalue;
              sl.next = loc->next;
            }
          } else {
            sl.choice = slp->choice;
            sl.data.ptrvalue = (Pointer) slp->data.ptrvalue;
            sl.next = NULL;
          }
          SeqPortStreamInt (&bsq, start, extend - 1, Seq_strand_plus, flags, (Pointer) str, NULL);
        } else {
          num = SeqPortStreamInt (bsp, start, extend - 1, Seq_strand_plus, flags, (Pointer) str, NULL);
          if (num < 1) {
            /* flag possible inconsistency between bsp->length and actual sequence data length */
            ajp->relModeError = TRUE;
            return NULL;
          }
        }
        /*
        if (ISA_aa (bsp->mol) && StringDoesHaveText (str)) {
          if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
            ChangeOandJtoX (str);
          }
        }
        */
        sbp->bases = str;
      }
    }
  }

  if (sbp->bases == NULL) return NULL;

  /* format subsequence cached with SeqPortStream */

  ffstring = FFGetString (ajp);

  start = sbp->start;
  stop = sbp->stop;
  remaining = stop - start;

  count = 0;
  blk = 0;
  lin = 0;

  ptr = sbp->bases;
  ch = *ptr;

  while (ch != '\0' && remaining > 0) {
    buf [count] = (Char) (TO_LOWER (ch));
    count++;
    remaining--;
    ptr++;
    ch = *ptr;

    blk++;
    lin++;
    if (lin >= 60) {

      buf [count] = '\0';
      startgapgap = 0;
      if (ajp->specialGapFormat) {
        startgapgap = ProcessGapSpecialFormat (afp, ajp, bsp, ffstring, buf, ptr);
      }
      if (StringDoesHaveText (buf)) {
        ExpandSeqLine (buf);
        PrintSeqLine (ffstring, afp->format, buf, start + startgapgap, start + lin);
      }
      count = 0;
      blk = 0;
      lin = 0;
      start += 60;
    }
  }

  buf [count] = '\0';
  if (count > 0) {
    startgapgap = 0;
    if (ajp->specialGapFormat) {
      startgapgap = ProcessGapSpecialFormat (afp, ajp, bsp, ffstring, buf, ptr);
    }
    if (StringDoesHaveText (buf)) {
      ExpandSeqLine (buf);
      PrintSeqLine (ffstring, afp->format, buf, start + startgapgap, start + lin);
    }
  }

  str = FFToCharPtr(ffstring);

  FFRecycleString (ajp, ffstring);
  return str;
}

/*
static CharPtr insd_strd [4] = {
  NULL, "single", "double", "mixed"
};

static CharPtr insd_mol [10] = {
  "?", "DNA", "RNA", "tRNA", "rRNA", "mRNA", "uRNA", "snRNA", "snoRNA", "AA"
};

static CharPtr insd_top [3] = {
  NULL, "linear", "circular"
};
*/

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

NLM_EXTERN CharPtr FormatSlashBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  GBFeaturePtr     currf, headf, nextf;
  GBReferencePtr   currr, headr, nextr;
  GBSeqPtr         gbseq, gbtmp;
  IndxPtr          index;
  INSDSeq          is;
  /*
  Int2              moltype, strandedness, topology;
  */

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  /* sort and unique indexes */

  index = ajp->index;

  if (index != NULL) {

    MemCopy (index, &asp->index, sizeof (IndxBlock));
    MemSet (&asp->index, 0, sizeof (IndxBlock));

    index->authors = ValNodeSort (index->authors, SortVnpByString);
    index->authors = UniqueValNode (index->authors);

    index->genes = ValNodeSort (index->genes, SortVnpByString);
    index->genes = UniqueValNode (index->genes);

    index->journals = ValNodeSort (index->journals, SortVnpByString);
    index->journals = UniqueValNode (index->journals);

    index->keywords = ValNodeSort (index->keywords, SortVnpByString);
    index->keywords = UniqueValNode (index->keywords);

    index->secondaries = ValNodeSort (index->secondaries, SortVnpByString);
    index->secondaries = UniqueValNode (index->secondaries);
  }

  /* adjust XML-ized GenBank format */

  gbseq = ajp->gbseq;

  if (gbseq != NULL) {

    MemCopy (gbseq, &asp->gbseq, sizeof (GBSeq));
    MemSet (&asp->gbseq, 0, sizeof (GBSeq));

    /* reverse order of references */

    headr = NULL;
    for (currr = gbseq->references; currr != NULL; currr = nextr) {
      nextr = currr->next;
      currr->next = headr;
      headr = currr;
    }
    gbseq->references = headr;

    /* reverse order of features */

    headf = NULL;
    for (currf = gbseq->feature_table; currf != NULL; currf = nextf) {
      nextf = currf->next;
      currf->next = headf;
      headf = currf;
    }
    gbseq->feature_table = headf;
  }

  /* if generating GBSeq XML/ASN, write at each slash block */

  if (gbseq != NULL && afp->aip != NULL) {
    if (ajp->produceInsdSeq) {
      MemSet ((Pointer) &is, 0, sizeof (INSDSeq));
      is.next = (INSDSeqPtr) gbseq->next;
      is.OBbits__ = gbseq->OBbits__;
      is.locus = gbseq->locus;
      is.length = gbseq->length;
      is.strandedness = gbseq->strandedness;
      is.moltype = gbseq->moltype;
      is.topology = gbseq->topology;
      /*
      strandedness = (Int2) gbseq->strandedness;
      if (strandedness < 0 || strandedness > 3) {
        strandedness = 0;
      }
      is.strandedness = StringSave (insd_strd [strandedness]);
      moltype = (Int2) gbseq->moltype;
      if (moltype < 0 || moltype > 9) {
        moltype = 0;
      }
      is.moltype = StringSave (insd_mol [moltype]);
      topology = (Int2) gbseq->topology;
      if (topology < 0 || topology > 2) {
        topology = 0;
      }
      is.topology = StringSave (insd_top [topology]);
      */
      is.division = gbseq->division;
      is.update_date = gbseq->update_date;
      is.create_date = gbseq->create_date;
      is.update_release = gbseq->update_release;
      is.create_release = gbseq->create_release;
      is.definition = gbseq->definition;
      is.primary_accession = gbseq->primary_accession;
      is.entry_version = gbseq->entry_version;
      is.accession_version = gbseq->accession_version;
      is.other_seqids = gbseq->other_seqids;
      is.secondary_accessions = gbseq->secondary_accessions;
      is.project = gbseq->project;
      is.keywords = gbseq->keywords;
      is.segment = gbseq->segment;
      is.source = gbseq->source;
      is.organism = gbseq->organism;
      is.taxonomy = gbseq->taxonomy;
      is.references = (INSDReferencePtr) gbseq->references;
      is.comment = gbseq->comment;
      is.primary = gbseq->primary;
      is.source_db = gbseq->source_db;
      is.database_reference = gbseq->database_reference;
      is.feature_table = (INSDFeaturePtr) gbseq->feature_table;
      is.sequence = gbseq->sequence;
      is.contig = gbseq->contig;
      INSDSeqAsnWrite (&is, afp->aip, afp->atp);
    } else {
      GBSeqAsnWrite (gbseq, afp->aip, afp->atp);
    }
    if (afp->atp == NULL) {
      AsnPrintNewLine (afp->aip);
    }
    AsnIoFlush (afp->aip);

    /* clean up gbseq fields */

    gbtmp = GBSeqNew ();
    MemCopy (gbtmp, gbseq, sizeof (GBSeq));
    MemSet (gbseq, 0, sizeof (GBSeq));
    GBSeqFree (gbtmp);
  }

  /* slash always has string pre-allocated by add slash block function */

  return StringSaveNoNull (bbp->string);
}


