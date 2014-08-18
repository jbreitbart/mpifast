/*   sqnutil4.c
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
* File Name:  sqnutil4.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   12/27/2007
*
* $Revision: 1.52 $
*
* File Description: 
* This file contains functions for automatically generating definition lines.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/
#include <sqnutils.h>
#include <ncbilang.h>
#include <objfdef.h>
#include <gather.h>
#include <explore.h>
#include <edutil.h>
#include <salutil.h>
#include <tofasta.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <findrepl.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

/* This is a list of the modifiers that are of interest */
/* Note that if you modify the DefLineModifiers array, */ 
/* you should make the corresponding change to the DefLinePos enum. */

ModifierItemGlobalData DefLineModifiers[] = {
  { "Acronym"              , TRUE , ORGMOD_acronym              },
  { "Anamorph"             , TRUE , ORGMOD_anamorph             },
  { "Authority"            , TRUE , ORGMOD_authority            },
  { "Bio-material"         , TRUE,  ORGMOD_bio_material         },
  { "Biotype"              , TRUE , ORGMOD_biotype              },
  { "Biovar"               , TRUE , ORGMOD_biovar               },
  { "Breed"                , TRUE , ORGMOD_breed                },
  { "Cell-line"            , FALSE, SUBSRC_cell_line            },
  { "Cell-type"            , FALSE, SUBSRC_cell_type            },
  { "Chemovar"             , TRUE , ORGMOD_chemovar             },
  { "Chromosome"           , FALSE, SUBSRC_chromosome           },
  { "Clone"                , FALSE, SUBSRC_clone                },
  { "Clone-lib"            , FALSE, SUBSRC_clone_lib            },
  { "Collected-by"         , FALSE, SUBSRC_collected_by         },
  { "Collection-date"      , FALSE, SUBSRC_collection_date      },
  { "Common"               , TRUE , ORGMOD_common               },
  { "Country"              , FALSE, SUBSRC_country              },
  { "Cultivar"             , TRUE , ORGMOD_cultivar             },
  { "Culture-collection"   , TRUE , ORGMOD_culture_collection   },
  { "Dev-stage"            , FALSE, SUBSRC_dev_stage            },
  { "Ecotype"              , TRUE , ORGMOD_ecotype              },
  { "Endogenous-virus-name", FALSE, SUBSRC_endogenous_virus_name},
  { "Environmental-sample" , FALSE, SUBSRC_environmental_sample },
  { "Forma"                , TRUE , ORGMOD_forma                },
  { "Forma-specialis"      , TRUE , ORGMOD_forma_specialis      },
  { "Frequency"            , FALSE, SUBSRC_frequency            },
  { "Genotype"             , FALSE, SUBSRC_genotype             },
  { "Germline"             , FALSE, SUBSRC_germline             },
  { "Group"                , TRUE , ORGMOD_group                },
  { "Haplogroup"           , FALSE, SUBSRC_haplogroup           },
  { "Haplotype"            , FALSE, SUBSRC_haplotype            },
  { "Host"                 , TRUE , ORGMOD_nat_host             },
  { "Identified-by"        , FALSE, SUBSRC_identified_by        },
  { "Isolate"              , TRUE , ORGMOD_isolate              },
  { "Isolation-source"     , FALSE, SUBSRC_isolation_source     },
  { "Lab-host"             , FALSE, SUBSRC_lab_host             },
  { "Lat-lon"              , FALSE, SUBSRC_lat_lon              },
  { "Linkage-group"        , FALSE, SUBSRC_linkage_group        },
  { "Map"                  , FALSE, SUBSRC_map                  },
  { "Mating-type"          , FALSE, SUBSRC_mating_type          },
  { "Metagenomic"          , FALSE, SUBSRC_metagenomic          },
  {	"Note-OrgMod"          , TRUE,  ORGMOD_other                },
  {	"Note-SubSrc"          , FALSE, SUBSRC_other                },
  { "Pathovar"             , TRUE , ORGMOD_pathovar             },
  { "Plasmid-name"         , FALSE, SUBSRC_plasmid_name         },
  { "Plastid-name"         , FALSE, SUBSRC_plastid_name         },
  { "Pop-variant"          , FALSE, SUBSRC_pop_variant          },
  { "Rearranged"           , FALSE, SUBSRC_rearranged           },
  { "Segment"              , FALSE, SUBSRC_segment              },
  { "Serogroup"            , TRUE , ORGMOD_serogroup            },
  { "Serotype"             , TRUE , ORGMOD_serotype             },
  { "Serovar"              , TRUE , ORGMOD_serovar              },
  { "Sex"                  , FALSE, SUBSRC_sex                  },
  { "Specimen voucher"     , TRUE , ORGMOD_specimen_voucher     },
  { "Strain"               , TRUE , ORGMOD_strain               },
  { "Subclone"             , FALSE, SUBSRC_subclone             },
  { "Subgroup"             , TRUE , ORGMOD_subgroup             },
  { "Sub-species"          , TRUE , ORGMOD_sub_species          },
  { "Substrain"            , TRUE , ORGMOD_substrain            },
  { "Subtype"              , TRUE , ORGMOD_subtype              },
  { "Synonym"              , TRUE , ORGMOD_synonym              },
  { "Teleomorph"           , TRUE , ORGMOD_teleomorph           },
  { "Tissue-lib"           , FALSE, SUBSRC_tissue_lib           },
  { "Tissue-type"          , FALSE, SUBSRC_tissue_type          },
  { "Transgenic"           , FALSE, SUBSRC_transgenic           },
  { "Type"                 , TRUE , ORGMOD_type                 },
  { "Variety"              , TRUE , ORGMOD_variety              }
};

#define numDefLineModifiers (sizeof (DefLineModifiers) / sizeof (ModifierItemGlobalData))

NLM_EXTERN size_t NumDefLineModifiers (void)

{
  return numDefLineModifiers;
}

NLM_EXTERN CharPtr MergeValNodeStrings (ValNodePtr list, Boolean useReturn)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  ptr = NULL;
  if (list != NULL) {
    vnp = list;
    len = 0;
    while (vnp != NULL) {
      if (vnp->data.ptrvalue != NULL) {
        len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
      }
      vnp = vnp->next;
    }
    if (len > 0) {
      ptr = MemNew (sizeof (Char) * (len + 2));
      if (ptr != NULL) {
        vnp = list;
        tmp = NULL;
        while (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            if (tmp == NULL) {
              tmp = ptr;
            } else if (useReturn) {
              tmp = StringMove (tmp, "\n");
            } else if (IsJapanese () && (tmp - ptr > 2) &&
            		IsMBLetter (tmp - 2) && IsMBLetter (str)) {
              /* no space required between two Japanese letters. */
              tmp = tmp;
            } else if (str [0] != ',' && str [0] != ';' && str [0] != ':') {
              tmp = StringMove (tmp, " ");
            } else {
              tmp = StringMove (tmp, " ");
            }
            tmp = StringMove (tmp, str);
          }
          vnp = vnp->next;
        }
      }
    }
  }
  return ptr;
}


/* The matchFunction functions are used to identify features that meet
 * specific requirements, usually that the feature is of a particular type.
 * This function is used instead of simply using the subtype for the feature
 * because some features are identified based on the contents or presence of
 * certain modifiers.
 * Functions of this type should always return FALSE if handed a NULL argument.
 */
typedef Boolean (LIBCALLBACK *matchFunction) (
  SeqFeatPtr sfp
);

static void ListClauses (
  ValNodePtr clauselist,
  ValNodePtr PNTR strings,
  Boolean    allow_semicolons,
  Boolean    suppress_final_and
);

static void LabelClauses 
( ValNodePtr clause_list,
  Uint1      biomol,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp);

static CharPtr GetProductName 
( SeqFeatPtr cds,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp);

#define DEFLINE_FEATLIST    1
#define DEFLINE_CLAUSEPLUS  2
#define DEFLINE_REMOVEFEAT  3

typedef struct featurelabeldata {
  Boolean pluralizable;
  Boolean is_typeword_first;
  CharPtr typeword;
  CharPtr description;
  CharPtr productname;
} FeatureLabelData, PNTR FeatureLabelPtr;


typedef struct featureclause {
  ValNodePtr       featlist;
  FeatureLabelData feature_label_data;
  CharPtr          allelename;
  CharPtr          interval;
  Boolean          is_alt_spliced;
  Boolean          has_mrna;
  SeqLocPtr        slp;
  GeneRefPtr       grp;
  Boolean          clause_info_only;
  Boolean          is_unknown;
  Boolean          make_plural;
  Boolean          delete_me;
  /* this information used only for segments */
  Int2             numivals;
  Int4Ptr          ivals;
} FeatureClauseData, PNTR FeatureClausePtr;

FeatureClausePtr NewFeatureClause (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  DeflineFeatureRequestListPtr rp);

static void PluralizeConsolidatedClauseDescription (
  FeatureClausePtr fcp
);

typedef Boolean (LIBCALLBACK *ShouldRemoveFunction) (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr        bsp,
  Boolean          isLonely,
  Boolean          isRequested,
  Boolean          isSegment,
  DeflineFeatureRequestListPtr rp
);

/* This section of the code contains some functions for dealing with
 * linked lists of strings */

/* This function finds the first occurrence of "search" in one of the
 * strings in list "strings".
 * "search" could be part of the string or could be the entire string.
 */
static ValNodePtr FindStringInStrings (
  ValNodePtr strings,
  CharPtr search
)
{
  while (strings != NULL)
  {
    if (StringStr (strings->data.ptrvalue, search))
    {
      return strings;
    }
    strings = strings->next;
  }
  return NULL;
}

/* This function finds the first item in "strings" that is identical to
 * "value".
 */
extern ValNodePtr FindExactStringInStrings (
  ValNodePtr strings,
  CharPtr    value
)
{
  ValNodePtr    string_match;

  for (string_match = FindStringInStrings (strings, value);
       string_match != NULL
         && StringCmp (string_match->data.ptrvalue, value) != 0;
       string_match = FindStringInStrings (string_match->next, value))
  {}
  return string_match;
}

/* This function creates a new linked list of strings with copies of
 * contents of orig.
 */
static ValNodePtr CopyStrings (
  ValNodePtr orig
)
{
  ValNodePtr new_string_start = NULL;

  while (orig != NULL)
  {
    ValNodeAddStr (&new_string_start, 0,
      StringSave (orig->data.ptrvalue));
    orig = orig->next;
  }
  return new_string_start;
}

/*
 * This section of the code contains functions and structures for obtaining a
 * description of the organism in the record, including functions for finding
 * the combination of modifiers that will make each organism description 
 * unique.
 * The method used for determining the best combination of modifiers involves
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * I then make copies of all of the new items I added to the list and
 * add another modifier to each list, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If I run out of possible combinations, I select the best combination from
 * the list.
 * This search process occurs in FindBestCombo.  The majority of the functions
 * in this section are here to support FindBestCombo, specifically to create,
 * copy, and grow lists of combinations. 
 */

/* BioSrcDescData is used to calculate the best possible combination of
 * source and organism modifiers for uniqueness.
 * biop contains the BioSourcePtr from a sequence in the record.
 * strings contains a list of string representations of the modifiers
 * for this combination for this organism.
 */
typedef struct biosrcdescdata {
  BioSourcePtr  biop;
  ValNodePtr    strings;
  Pointer       next;
} BioSrcDescData, PNTR BioSrcDescPtr;

/* OrgGroupData is used to calculate the best possible combination of
 * source and organism modifiers for uniqueness.
 * org_list is a list of all organisms that have identical descriptions
 * using the current set of modifiers.
 * num_organisms contains the number of organisms with identical descriptions.
 */
typedef struct orggroupdata {
  BioSrcDescPtr org_list;
  Int4          num_organisms;
  Pointer       next;
} OrgGroupData, PNTR OrgGroupPtr;

/* ModifierCombinationData is used to calculate the best possible combination
 * of source and organism modifiers for uniqueness.
 * num_groups is the number of groups of organisms with identical descriptions
 *           using the modifiers specified in modifier_indices.
 * num_mods is the number of modifiers specified in modifier_indices.
 * max_orgs_in_group is the maximum number of organisms in any one group.
 * num_unique_orgs is the number of organisms that are alone in their groups
 *           i.e., their description is unique.
 * modifier_indices is the list of modifier indices for this combination.
 * group_list is the list of groups of organisms with identical descriptions
 *           using the modifiers specified in modifier_indices.
 */
typedef struct modifiercombination {
  Int4         num_groups;
  Int4         num_mods;
  Int4         max_orgs_in_group;
  Int4         num_unique_orgs;
  ValNodePtr   modifier_indices;
  OrgGroupPtr  group_list;
  Pointer      next;
} ModifierCombinationData, PNTR ModifierCombinationPtr;

static Boolean IsDeflineModifierRequiredByDefault (Boolean is_orgmod, Int2 index)
{
  if (!is_orgmod
      && (index == SUBSRC_endogenous_virus_name
          || index == SUBSRC_plasmid_name
          || index == SUBSRC_transgenic)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void AddOneSubtypeField (ValNodePtr PNTR sq_list, SourceQualDescPtr orig, CharPtr str, Uint1 subfield)
{
  SourceQualDescPtr sqdp_cpy;

  if (sq_list == NULL || orig == NULL) {
    return;
  }
  sqdp_cpy = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
  MemCpy (sqdp_cpy, orig, sizeof (SourceQualDescData));

  sqdp_cpy->name = str;
  sqdp_cpy->subfield = subfield;

  ValNodeAddPointer (sq_list, 0, sqdp_cpy);  
}


static void AddSubtypeFields (ValNodePtr PNTR sq_list, SourceQualDescPtr orig)
{
  if (sq_list == NULL || orig == NULL) return;

  if (orig->isOrgMod) {
    switch (orig->subtype) {
      case ORGMOD_specimen_voucher:
        AddOneSubtypeField (sq_list, orig, "specimen-voucher INST", 1);
        AddOneSubtypeField (sq_list, orig, "specimen-voucher COLL", 2);
        AddOneSubtypeField (sq_list, orig, "specimen-voucher SpecID", 3);
        break;  
      case ORGMOD_culture_collection:
        AddOneSubtypeField (sq_list, orig, "culture-collection INST", 1);
        AddOneSubtypeField (sq_list, orig, "culture-collection COLL", 2);
        AddOneSubtypeField (sq_list, orig, "culture-collection SpecID", 3);
        break;  
      case ORGMOD_bio_material:
        AddOneSubtypeField (sq_list, orig, "bio-material INST", 1);
        AddOneSubtypeField (sq_list, orig, "bio-material COLL", 2);
        AddOneSubtypeField (sq_list, orig, "bio-material SpecID", 3);
        break; 
    }
  } 
}


static void AddQualList (ValNodePtr PNTR list, Nlm_QualNameAssocPtr qual_list, Boolean is_orgmod, Boolean use_alternate_note_name)
{
  Int4              k;
  SourceQualDescPtr sqdp;
  
  for (k = 0; qual_list[k].name != NULL; k++) {
    if (StringHasNoText (qual_list[k].name)) {
      continue;
    } 
    sqdp = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
    if (sqdp != NULL)
    {
      if (use_alternate_note_name 
          && ((is_orgmod && qual_list[k].value == ORGMOD_other)
              || (!is_orgmod && qual_list[k].value == SUBSRC_other))) 
      {
        if (is_orgmod) {
          sqdp->name = "Note -- OrgMod";
        } else {
          sqdp->name = "Note -- SubSource";
        }
      } else {
        sqdp->name = qual_list[k].name;
      }
      sqdp->isOrgMod = is_orgmod;
      sqdp->subtype = qual_list[k].value;
      sqdp->subfield = 0;
      ValNodeAddPointer (list, 0, sqdp);
    }
    AddSubtypeFields (list, sqdp);
  }
}


static void AddNoteQual (ValNodePtr PNTR list, Boolean is_orgmod, Boolean use_alternate_note_name)
{
  SourceQualDescPtr sqdp;

  if (list == NULL) return;

  sqdp = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
  if (sqdp != NULL)
  {
    if (use_alternate_note_name) 
    {
      if (is_orgmod) {
        sqdp->name = "Note -- OrgMod";
      } else {
        sqdp->name = "Note -- SubSource";
      }
    } else {
      sqdp->name = "Note";
    }
    sqdp->isOrgMod = is_orgmod;
    if (is_orgmod) {
      sqdp->subtype = ORGMOD_other;
    } else {
      sqdp->subtype = SUBSRC_other;
    }
    sqdp->subfield = 0;
    ValNodeAddPointer (list, 0, sqdp);
  }
}


static int LIBCALLBACK SortVnpBySourceQualDesc (VoidPtr ptr1, VoidPtr ptr2)

{
  SourceQualDescPtr     str1;
  SourceQualDescPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (SourceQualDescPtr) vnp1->data.ptrvalue;
      str2 = (SourceQualDescPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL
          && str1->name != NULL && str2->name != NULL) {
        return StringICmp (str1->name, str2->name);
      }
    }
  }
  return 0;
}


extern ValNodePtr GetSourceQualDescList (Boolean get_subsrc, Boolean get_orgmod, Boolean get_discouraged, Boolean get_discontinued)
{
  ValNodePtr        source_qual_list = NULL;

  if (get_orgmod) {
    AddQualList (&source_qual_list, current_orgmod_subtype_alist, TRUE, get_subsrc);
    if (get_discouraged) {
      AddQualList (&source_qual_list, discouraged_orgmod_subtype_alist, TRUE, get_subsrc);
    }
    if (get_discontinued) {
      AddQualList (&source_qual_list, discontinued_orgmod_subtype_alist, TRUE, get_subsrc);
    }
    AddNoteQual (&source_qual_list, TRUE, get_subsrc);
  }
  if (get_subsrc) {
    AddQualList (&source_qual_list, current_subsource_subtype_alist, FALSE, get_orgmod);
    if (get_discouraged) {
      AddQualList (&source_qual_list, discouraged_subsource_subtype_alist, FALSE, get_orgmod);
    }
    if (get_discontinued) {
      AddQualList (&source_qual_list, discontinued_subsource_subtype_alist, FALSE, get_orgmod);
    }
    AddNoteQual (&source_qual_list, FALSE, get_orgmod);
  }

  source_qual_list = ValNodeSort (source_qual_list, SortVnpBySourceQualDesc);
  return source_qual_list;
}

/*
 * The CountModifiersProc is used as the callback function for
 * VisitBioSourcesInSep when we are getting a list of all the modifiers
 * that appear in the sources.  We also obtain, for each modifier class,
 * the first value seen, whether or not each value seen is unique for
 * for the modifier, and whether or not the modifier is present for all
 * sources.
 */
static void CountModifiersProc (
  BioSourcePtr biop,
  Pointer userdata
)
{
  ModifierItemLocalPtr ItemList;
  OrgModPtr     mod;
  SubSourcePtr  ssp;
  Int2 i;
  Boolean       found_this_modifier;

  if (biop == NULL) return;
  ItemList = (ModifierItemLocalPtr) userdata;

  for (i=0; i < numDefLineModifiers; i++)
  {
    found_this_modifier = FALSE;
    if (DefLineModifiers[i].isOrgMod)
    {
      if ( biop->org != NULL && biop->org->orgname != NULL) 
      {
        mod = biop->org->orgname->mod;
        while (mod != NULL
          && mod->subtype != DefLineModifiers[i].subtype)
        {
          mod = mod->next;
        }
        if (mod != NULL && mod->subname != NULL)
        {
          found_this_modifier = TRUE;
          if (ItemList[i].first_value_seen != NULL)
          {
            if (StringCmp (ItemList[i].first_value_seen, mod->subname) != 0)
            {
              ItemList[i].is_unique = FALSE;
            }
          }
          else
          {
            ItemList[i].first_value_seen = mod->subname;
          }
          if ( FindExactStringInStrings (ItemList[i].values_seen, mod->subname)
            == NULL)
          {
            ValNodeAddStr (&ItemList[i].values_seen, 0, mod->subname);
          }
          else
          {
            ItemList[i].all_unique = FALSE;
          }
        }
      }
    } else {
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != DefLineModifiers[i].subtype)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL && ssp->name != NULL)
      {
        found_this_modifier = TRUE;
        if (ItemList[i].first_value_seen != NULL)
        {
          if (StringCmp (ItemList[i].first_value_seen, ssp->name) != 0)
          {
            ItemList[i].is_unique = FALSE;
          }
        }
        else
        {
          ItemList[i].first_value_seen = ssp->name;
        }
        if ( FindExactStringInStrings (ItemList[i].values_seen, ssp->name)
          == NULL)
        {
          ValNodeAddStr (&ItemList[i].values_seen, 0, ssp->name);
        }
        else
        {
          ItemList[i].all_unique = FALSE;
        }
      }
    }
    if (found_this_modifier)
    {
      ItemList[i].any_present = TRUE;
    } else {
      ItemList[i].all_present = FALSE;
    }
  }
} 

/* The CountModifiers function visits all of the bio sources, determining
 * which modifiers are present, which modifiers have only one value,
 * which modifiers have all different values, and which modifiers are
 * present in all sources.
 * After this survey is complete, the function prepares a short summary
 * of the above information for each modifier, which is used in the 
 * definition line options dialog.
 */ 
NLM_EXTERN void CountModifiers (
  ModifierItemLocalPtr ItemList,
  SeqEntryPtr sep
)
{
  Int2 i;

  for (i=0; i < numDefLineModifiers; i++)
  {
    ItemList[i].all_present = TRUE;
    ItemList[i].is_unique = TRUE;
    ItemList[i].first_value_seen = NULL;
    ItemList[i].values_seen = NULL;
    ItemList[i].all_unique = TRUE;
  }

  VisitBioSourcesInSep (sep, ItemList, CountModifiersProc);
 
  for (i=0; i < numDefLineModifiers; i++)
  {
    if (ItemList[i].all_present && ItemList[i].all_unique)
    {
      ItemList[i].status = "All present, all unique";
    }
    else if (ItemList[i].all_present && ItemList[i].is_unique)
    {
      ItemList[i].status = "All present, one unique";
    }
    else if (ItemList[i].all_present && ! ItemList[i].is_unique)
    {
      ItemList[i].status = "All present, mixed";
    }
    else if (! ItemList[i].all_present && ItemList[i].all_unique)
    {
      ItemList[i].status = "Some missing, all unique";
    }
    else if (! ItemList[i].all_present && ItemList[i].is_unique)
    {
      ItemList[i].status = "Some missing, one unique";
    }
    else if (! ItemList[i].all_present && ! ItemList[i].is_unique)
    {
      ItemList[i].status = "Some missing, mixed";
    }
  }
}

/* The BioSrcDescData structure is used to hold a BioSourcePtr, a list
 * of strings used to describe the biosource, including the taxonomy name
 * and the values of all of the modifiers selected so far for this bio
 * source, and a pointer to the next BioSrcDescData structure in the list.
 */

/* The CopyBioSrcDescPtr function creates a copy of the linked list of
 * BioSrcDescData structures.
 */
static BioSrcDescPtr CopyBioSrcDescPtr (
  BioSrcDescPtr orig
)
{
  BioSrcDescPtr new_bsdp_start;

  if (orig == NULL) return NULL;

  new_bsdp_start = (BioSrcDescPtr) MemNew (sizeof (BioSrcDescData));
  if (new_bsdp_start == NULL) return NULL;

  new_bsdp_start->biop = orig->biop;
  new_bsdp_start->strings = CopyStrings (orig->strings);
  new_bsdp_start->next = CopyBioSrcDescPtr (orig->next);
  return new_bsdp_start;
}

/* The FreeBioSrcDescPtr function frees the memory associated with a
 * linked list of BioSrcDescData structures.
 */
static void FreeBioSrcDescPtr (
  BioSrcDescPtr bsdp
)
{
  if (bsdp == NULL) return;
  FreeBioSrcDescPtr (bsdp->next);
  bsdp->biop = NULL;
  ValNodeFreeData (bsdp->strings);
  MemFree (bsdp);
}

/* The AddQualToBioSrcDescPtr function finds the qualifier at the
 * feature_index position in the DefLineModifiers array in the
 * BioSourcePtr and adds the value for that modifier to the array
 * of strings describing the bio source.
 */
static void AddQualToBioSrcDescPtr (
  BioSrcDescPtr bsdp,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgModPtr          mod;
  SubSourcePtr       ssp;
  CharPtr            tmp;

  if (bsdp == NULL) return;
  if (bsdp->biop == NULL) return;

  if (DefLineModifiers[feature_index].isOrgMod)
  {
    if (bsdp->biop->org == NULL || bsdp->biop->org->orgname == NULL) return;
    mod = bsdp->biop->org->orgname->mod;
    while (mod != NULL
        && mod->subtype != DefLineModifiers[feature_index].subtype)
    {
      mod = mod->next;
    }
    if (mod != NULL && mod->subname != NULL)
    {
      if (mod->subtype == ORGMOD_specimen_voucher && StringNICmp (mod->subname, "personal:", 9) == 0)
      {
        tmp = mod->subname + 9;
        while (isspace (*tmp)) 
        {
          tmp++;
        }
      }
      else
      {
        tmp = mod->subname;
      }
      ValNodeCopyStr( &(bsdp->strings), 0, tmp);
    } 
  } else {
    ssp = bsdp->biop->subtype;
    while (ssp != NULL
        && ssp->subtype != DefLineModifiers[feature_index].subtype)
    {
      ssp = ssp->next;
    }
    if (ssp != NULL)
    {
      if (ssp->subtype == SUBSRC_transgenic)
      {
        ValNodeCopyStr( &(bsdp->strings), 0, "transgenic");
      }
      else if (ssp->name != NULL)
      {
        ValNodeCopyStr( &(bsdp->strings), 0, ssp->name);
      }
    }
  }
}
 
/* The CompareOrganismDescriptors function compares the contents of the
 * lists of strings for each BioSrcDesc item.
 * The function returns:
 *     -1 if org1 < org2
 *      0 if org1 = org2
 *      1 if org1 > org2
 */
static int CompareOrganismDescriptors (
  BioSrcDescPtr org1,
  BioSrcDescPtr org2
)
{
  ValNodePtr vnp1, vnp2;
  int cmpval;

  vnp1 = org1->strings;
  vnp2 = org2->strings;

  while (vnp1 != NULL && vnp2 != NULL)
  {
    cmpval = StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    if (cmpval != 0) return cmpval;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return 0;
  }
  else if (vnp1 != NULL && vnp2 == NULL)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

/* The OrgGroupData structure contains a list of BioSrcDescData items
 * for which the contents of the descriptive strings list are identical,
 * i.e., all the organisms in the group would have the same description
 * if you used the modifiers used to generate this list of strings.
 * The structure also contains the number of organisms in the list
 * so that it will be easy to tell that the OrgGroup now contains a
 * single organism with a unique description.
 */

/* The CopyOrgGroupList function creates a copy of the list of OrgGroups */
static OrgGroupPtr CopyOrgGroupList (
  OrgGroupPtr orig
)
{
  OrgGroupPtr new_ogp_start = NULL, new_ogp;
 
  if (orig == NULL) return NULL;

  new_ogp_start = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
  if (new_ogp_start == NULL) return NULL;

  new_ogp_start->num_organisms = orig->num_organisms;
  new_ogp_start->org_list = CopyBioSrcDescPtr (orig->org_list);
  new_ogp_start->next = NULL;
  orig = orig->next;
  new_ogp = new_ogp_start;
  while (orig != NULL) {
    new_ogp->next = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
    new_ogp = new_ogp->next;
    new_ogp->num_organisms = orig->num_organisms;
    new_ogp->org_list = CopyBioSrcDescPtr (orig->org_list);
    new_ogp->next = NULL;
    orig = orig->next;
  }

  return new_ogp_start;
}

/* The FreeOrgGroupPtr function frees the memory associated with a
 * list of OrgGroups */
static void FreeOrgGroupPtr (
  OrgGroupPtr ogp
)
{
  OrgGroupPtr ogp_next;

  while (ogp != NULL) {
    ogp_next = ogp->next;
    FreeBioSrcDescPtr (ogp->org_list);
    ogp = MemFree (ogp);
    ogp = ogp_next;
  }
  return;
}

/* The ReorderGroupOrgs function sorts the OrgGroup list based on the results
 * of the CompareOrganismDescriptors function.
 */
static void ReorderGroupOrgs (
  OrgGroupPtr this_group
)
{
  BioSrcDescPtr bsdp;
  BioSrcDescPtr nextBsdp;
  BioSrcDescPtr prevBsdp;
  Boolean swap_needed = TRUE;

  if (this_group->org_list == NULL) return;
  if (this_group->org_list->next == NULL) return;

  while (swap_needed)
  {
    swap_needed = FALSE;
    bsdp = this_group->org_list;
    prevBsdp = NULL;
    while (bsdp->next != NULL)
    {
      nextBsdp = bsdp->next;
      if (CompareOrganismDescriptors (bsdp, nextBsdp) > 0)
      {
        swap_needed = TRUE;
        bsdp->next = nextBsdp->next;
        nextBsdp->next = bsdp;
        if (prevBsdp == NULL)
        {
          this_group->org_list = nextBsdp;
        }
        else
        {
          prevBsdp->next = nextBsdp;
        }
        prevBsdp = nextBsdp;
      }
      else
      {
        prevBsdp = bsdp;
        bsdp = bsdp->next;
      }
    }
  }
}

/* The ReGroupOrgs function operates on a single OrgGroup item.
 * If any of the BioSrcDesc items in the group now have different
 * descriptions, the function breaks it up into smaller, homogenous OrgGroups.
 */
static void ReGroupOrgs (
  OrgGroupPtr this_group
)
{
  BioSrcDescPtr bsdp;
  OrgGroupPtr new_group;
  int num_organisms;

  if (this_group == NULL) return;
  bsdp = this_group->org_list;
  if (bsdp == NULL) return;
  num_organisms = 0;
  while (bsdp->next != NULL)
  {
    num_organisms ++;
    if (CompareOrganismDescriptors (bsdp, bsdp->next) != 0)
    {
      /* create new group to hold next set of organisms */
      new_group = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
      if (new_group == NULL) return;
      new_group->org_list = bsdp->next;
      new_group->num_organisms = this_group->num_organisms - num_organisms;
      new_group->next = this_group->next;
      this_group->next = new_group;
      this_group->num_organisms = num_organisms;
      bsdp->next = NULL;
      ReGroupOrgs (new_group);
    }
    else
    {
      bsdp = bsdp->next;
    }
  }
}

/* The AddQualToGroup function operates on a single OrgGroup item.
 * The function adds a qualifier to each BioSrcDesc item in the OrgGroup,
 * breaks the group into multiple groups if the group is no longer
 * homogenous, and sorts the new list.
 */
static void AddQualToGroup (
  OrgGroupPtr this_group,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  BioSrcDescPtr bsdp;

  if (this_group == NULL) return;

  bsdp = this_group->org_list;
  while (bsdp != NULL)
  {
    AddQualToBioSrcDescPtr (bsdp, qual, feature_index);
    bsdp= bsdp->next;
  }

  /* now reorder organisms and break up group */
  ReorderGroupOrgs (this_group);

  ReGroupOrgs (this_group);
}

/* The AddQualToGroupList function operates on a list of OrgGroup items.
 * It calls AddQualToGroup for each item in the list.
 */
static void AddQualToGroupList (
  OrgGroupPtr group_list,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgGroupPtr ogp;

  ogp = group_list;
  while (ogp != NULL)
  {
    AddQualToGroup (ogp, qual, feature_index);
    ogp = ogp->next;
  }
}

/* The CopyModifierIndices function creates a new ValNode list with the
 * same data.intvalue values for each node as the original modifier_indices
 * ValNode list.
 */
static ValNodePtr CopyModifierIndices (
  ValNodePtr modifier_indices
)
{
  ValNodePtr new_indices;

  if (modifier_indices == NULL) return NULL;
  new_indices = ValNodeNew (NULL);
  if (new_indices == NULL) return NULL;
  new_indices->choice = modifier_indices->choice;
  new_indices->data.intvalue = modifier_indices->data.intvalue;
  new_indices->next = CopyModifierIndices (modifier_indices->next);
  return new_indices;
}
 
/* The CopyModifierCombo creates a copy of a ModificationCombination item.
 * This includes creating a copy of the number and list of modifiers
 * and a copy of the number and list of OrgGroups, as well as copying the
 * maximum number of organisms in any one group and the number of unique
 * organism descriptions produced by this combination of modifiers.
 */
static ModifierCombinationPtr CopyModifierCombo (
  ModifierCombinationPtr m
)
{
  ModifierCombinationPtr newm;
  ValNodePtr  vnp;
  ValNodePtr  newval;

  newm = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  if (newm == NULL) return NULL;

  newm->next = NULL;

  /* copy list of modifier indices */
  newm->num_mods = m->num_mods;
  newm->modifier_indices = NULL;
  vnp = m->modifier_indices;
  if (vnp != NULL)
  {
    newm->modifier_indices = ValNodeNew (NULL);
    if (newm->modifier_indices == NULL) return NULL;
    newm->modifier_indices->data.intvalue = vnp->data.intvalue;
    vnp = vnp->next;
    while (vnp != NULL)
    {
      newval = ValNodeNew (newm->modifier_indices);
      if (newval == NULL) return NULL;
      newval->data.intvalue = vnp->data.intvalue;
      vnp = vnp->next;
    }
  }
  
  /* copy groups */
  newm->num_groups = m->num_groups;
  newm->group_list = CopyOrgGroupList (m->group_list);
 
  return newm; 
}

/* This function creates a new ModifierCombination item using the supplied
 * OrgGroup list.  It calculates the number of groups, maximum number of 
 * organisms in any one group, and number of unique organisms.
 * Initially there are no modifiers.
 */
static ModifierCombinationPtr NewModifierCombo (
  OrgGroupPtr group_list
)
{
  ModifierCombinationPtr newm;
  OrgGroupPtr  ogp;

  newm = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  if (newm == NULL) return NULL;

  newm->num_mods = 0;
  newm->modifier_indices = NULL;
  newm->num_unique_orgs = 0;
  
  /* copy groups */
  newm->group_list = CopyOrgGroupList (group_list);
 
  ogp = newm->group_list;
  newm->max_orgs_in_group = 0;
  newm->num_groups = 0;
  while (ogp != NULL)
  {
    if (newm->max_orgs_in_group < ogp->num_organisms)
      newm->max_orgs_in_group = ogp->num_organisms;
    if (ogp->num_organisms == 1)
      newm->num_unique_orgs ++;
    newm->num_groups ++;
    ogp = ogp->next;
  }

  newm->next = NULL;
  return newm; 
}

/* This function frees the memory associated with a list of
 * ModifierCombination items.
 */
static void FreeModifierCombo (
  ModifierCombinationPtr m
)
{
  if (m == NULL) return;
  FreeModifierCombo (m->next);
  ValNodeFree (m->modifier_indices);
  FreeOrgGroupPtr (m->group_list);
  MemFree (m);
}

/* This function adds the qualifier at the feature_index position in the
 * DefLineModifiers array to each OrgGroup in the list and recalculates
 * the maximum number of organisms in any one group and the number of
 * unique organism descriptions generated by this new combination of 
 * modifiers.
 */
static void AddQualToModifierCombo (
  ModifierCombinationPtr m,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgGroupPtr ogp;
  ValNodePtr vnp;

  if (m == NULL) return;

  /* now try adding the modifier, see if the number of groups goes up */
  /* if the number of organisms in each group is one, we can stop */
  vnp = ValNodeNew (m->modifier_indices);
  if (vnp == NULL) return;
  if (m->modifier_indices == NULL)
  {
    m->modifier_indices = vnp;
  }
  vnp->data.intvalue = feature_index;
  m->num_mods ++;
  AddQualToGroupList (m->group_list, qual, feature_index);
  ogp = m->group_list;
  m->max_orgs_in_group = 0;
  m->num_unique_orgs = 0;
  m->num_groups = 0;
  while (ogp != NULL)
  {
    if (m->max_orgs_in_group < ogp->num_organisms)
      m->max_orgs_in_group = ogp->num_organisms;
    if (ogp->num_organisms == 1)
      m->num_unique_orgs ++;
    m->num_groups ++;
    ogp = ogp->next;
  }
}

/* This function creates the initial OrgGroup list that is copied for every
 * ModifierCombination item.
 */
static void BuildTaxOrgGroupList (
  BioSourcePtr biop,
  Pointer userdata
)
{
  OrgGroupPtr   ogp;
  OrgGroupPtr   prevOgp;
  OrgGroupPtr PNTR pogp;
  BioSrcDescPtr newBsdp;
  OrgRefPtr     orp;
  int cmpval;

  pogp = (OrgGroupPtr PNTR) userdata;
  ogp = *pogp;

  newBsdp = (BioSrcDescPtr) MemNew (sizeof (BioSrcDescData));
  if (newBsdp == NULL) return;
  newBsdp->biop = biop;
  newBsdp->next = NULL;
  newBsdp->strings = NULL;

  /* add tax name as first string */
  /* later, move this into a separate function and add special handling */
  orp = biop->org;
  if (orp != NULL && orp->taxname != NULL)
  {
    ValNodeCopyStr (&(newBsdp->strings), 0, orp->taxname);
  }

  prevOgp = NULL;
  cmpval = -1;
  while (ogp != NULL && cmpval < 0)
  {
    if (ogp->org_list != NULL)
    {
      cmpval = CompareOrganismDescriptors (ogp->org_list, newBsdp);
      if (cmpval == 0)
      {
        newBsdp->next = ogp->org_list;
        ogp->org_list = newBsdp;
        ogp->num_organisms ++;
      }
    }
    if (cmpval < 0)
    {
      prevOgp = ogp;
      ogp = ogp->next;
    }
  }
  if (cmpval != 0)
  {
    /* create new group */
    ogp = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
    if (ogp == NULL) return;
    ogp->org_list = newBsdp;
    ogp->num_organisms = 1;
    ogp->next = NULL;
    if (prevOgp == NULL)
    {
      ogp->next = *pogp;
      *pogp = ogp;
    }
    else
    {
      ogp->next = prevOgp->next;
      prevOgp->next = ogp;
    }
  }
}

typedef struct bestsortdata {
  Int4    feature_index;
  Boolean all_unique;
  Boolean all_present;
  Boolean is_unique;
} BestSortData, PNTR BestSortPtr;

static Boolean Index1FoundBeforeIndex2 (
  Int4 index1,
  Int4 index2,
  ValNodePtr list
)
{
  ValNodePtr  vnp;
  BestSortPtr bsp;
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if ((bsp = vnp->data.ptrvalue) == NULL)
    {
      continue;
    }
    if (bsp->feature_index == index1) return TRUE;
    if (bsp->feature_index == index2) return FALSE;
  }
  return FALSE;
}

/* This function determines whether or not we should try adding this modifier
 * to our combination.  If we've already tried it and not added it to the list,
 * there's no reason to try adding it again.
 */
static Boolean OkToTryAddingQual (
  ModifierCombinationPtr m,
  ModifierItemLocalPtr ItemList,
  ValNodePtr           available_modifiers_list,
  Int2 feature_index
)
{
  ValNodePtr vnp;

  /* if feature_index indicates a value we don't use for best combos, skip */
  if (feature_index == DEFLINE_POS_Map)
  {
    return FALSE;
  }

  if (m == NULL) return TRUE;

  /* if feature_index is lower than anything else on list (other than */
  /* a required value, this is a repeat combination, so skip it */
  vnp = m->modifier_indices;
  while (vnp != NULL)
  {
    if (feature_index == m->modifier_indices->data.intvalue)
      return FALSE;
    if (! ItemList[m->modifier_indices->data.intvalue].required &&
      Index1FoundBeforeIndex2 (feature_index,
                               m->modifier_indices->data.intvalue, 
                               available_modifiers_list))
    {
      return FALSE;
    }
    vnp = vnp->next;
  }
  return TRUE;
}

static ValNodePtr GetListOfAvailableModifiers ( ModifierItemLocalPtr ItemList)
{
  ValNodePtr  vnp, head;
  Int2        feature_index;
  BestSortPtr bsp;

  head = NULL;
  for (feature_index = 0; feature_index < numDefLineModifiers; feature_index++)
  {
    if ( ItemList[feature_index].any_present)
    {
      bsp = (BestSortPtr) MemNew (sizeof (BestSortData));
      if (bsp == NULL) return NULL;
      bsp->feature_index = feature_index;
      bsp->all_unique = ItemList[feature_index].all_unique;
      bsp->all_present = ItemList[feature_index].all_present;
      bsp->is_unique = ItemList[feature_index].is_unique;
      vnp = ValNodeNew (head);
      if (vnp == NULL) return NULL;
      vnp->data.ptrvalue = bsp;
      if (head == NULL) head = vnp;
    }
  }
  return head;
}

static Int4 DefLineQualSortOrder [] = {
  DEFLINE_POS_Transgenic,
  DEFLINE_POS_Plasmid_name,
  DEFLINE_POS_Endogenous_virus_name,
  DEFLINE_POS_Strain,
  DEFLINE_POS_Clone,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Haplotype,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Specimen_voucher,
  DEFLINE_POS_Ecotype,
  DEFLINE_POS_Type,
  DEFLINE_POS_Serotype,
  DEFLINE_POS_Authority,
  DEFLINE_POS_Breed
};

static int LIBCALLBACK SortByImportanceAndPresence (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  BestSortPtr bsp1, bsp2; 
  Int4       num_defline_qual_sort_order, index;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;

  bsp1 = vnp1->data.ptrvalue;
  bsp2 = vnp2->data.ptrvalue;
  if (bsp1->feature_index == bsp2->feature_index) return 0;

  if (bsp1->all_present && bsp1->all_unique
    && (! bsp2->all_present || ! bsp2->all_unique))
  {
    return -1;
  }
  if (bsp2->all_present && bsp2->all_unique
    && (! bsp1->all_present || ! bsp1->all_unique))
  {
    return 1;
  }

  if ( ! bsp1->is_unique && bsp2->is_unique) return -1;
  if ( ! bsp2->is_unique && bsp1->is_unique) return 1;

  num_defline_qual_sort_order = sizeof (DefLineQualSortOrder) / sizeof (Int4);
  for (index = 0; index < num_defline_qual_sort_order; index++)
  {
    if (bsp1->feature_index == DefLineQualSortOrder [ index ]) return -1;
    if (bsp2->feature_index == DefLineQualSortOrder [ index ]) return 1;
  }

  if (bsp1->feature_index > bsp2->feature_index) return 1;
  if (bsp1->feature_index < bsp2->feature_index) return -1;
  return 0;
}


/* The function FindBestCombo tries to find the best combination of modifiers
 * to create unique organism descriptions.  This is accomplished by
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * The function then makes copies of all of the new items added to the list,
 * starting with the item pointed to by start_of_expand, and adds another
 * modifier to each combination, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If the list of possible combinations is exhausted before each organism
 * has a unique description, the function selects the combination from the
 * list with the largest number of unique organism descriptions.  If more 
 * than one combination produces the largest number of unique organisms,
 * the combination with the largest number of unique organisms and the
 * largest number of groups will be selected.
 */
static ModifierCombinationPtr FindBestCombo(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList
)
{
  OrgGroupPtr group_list;
  ModifierCombinationPtr mc_list, start_of_expand, best_found, end_of_list;
  ModifierCombinationPtr next_start_of_expand, m, newm;
  Int4 num_to_expand, next_num_to_expand;
  Int2 i;
  ValNodePtr available_modifier_list, vnp;
  BestSortPtr bsp;

  best_found = NULL;

  /* first, get list of organisms */
  group_list = NULL;
  VisitBioSourcesInSep (sep, &group_list, BuildTaxOrgGroupList);

  /* create combo with just the org groups */
  mc_list = NewModifierCombo (group_list);
  if (mc_list == NULL) return NULL;

  available_modifier_list = GetListOfAvailableModifiers (ItemList);

  /* next, add in any required qualifiers */
  for (vnp = available_modifier_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = vnp->data.ptrvalue;
    if (bsp == NULL) return NULL;
    if (ItemList[bsp->feature_index].required)
    {
      AddQualToModifierCombo (mc_list, ItemList + bsp->feature_index,
                                       bsp->feature_index);
    }
  }
  if (mc_list->max_orgs_in_group == 1)
  {
    /* we're done - they're all unique */
    best_found = mc_list;
    return best_found;
  }

  available_modifier_list = ValNodeSort (available_modifier_list,
                                         SortByImportanceAndPresence);
  start_of_expand = mc_list;
  end_of_list = mc_list;
  num_to_expand = 1;
  while (best_found == NULL && start_of_expand != NULL)
  {
    next_num_to_expand = 0;
    next_start_of_expand = NULL;
    for (i=0; i < num_to_expand && start_of_expand != NULL; i++)
    {
      /* try adding qualifiers */ 
      for (vnp = available_modifier_list;
           vnp != NULL && best_found == NULL;
           vnp = vnp->next)
      {
        bsp = vnp->data.ptrvalue;
        if (bsp == NULL) return NULL;
        if (OkToTryAddingQual (start_of_expand, ItemList,
                               available_modifier_list,
                               bsp->feature_index))
        {
          newm = CopyModifierCombo (start_of_expand);
          AddQualToModifierCombo (newm, ItemList + bsp->feature_index,
                                  bsp->feature_index);
          if (start_of_expand->num_groups >= newm->num_groups)
          {
            /* situation didn't get better, don't bother to add this one */
            FreeModifierCombo (newm);
            newm = NULL;
          }
          else if (newm->max_orgs_in_group == 1)
          {
            best_found = newm;
          }
          else
          {
            end_of_list->next = newm;
            end_of_list = end_of_list->next;
            if (next_start_of_expand == NULL)
              next_start_of_expand = newm;
            next_num_to_expand++;
          }
        }
      }
      if (start_of_expand != NULL)
      {
        start_of_expand = start_of_expand->next;
      }
    }
    num_to_expand = next_num_to_expand;
    if (start_of_expand != NULL)
    {
      start_of_expand = start_of_expand->next;
    }
  }

  if (best_found != NULL)
  {
    FreeModifierCombo (mc_list);
    return best_found;
  }
  
  /* we want to find the one with the highest number of unique organisms */
  best_found = mc_list;
  m = mc_list->next;
  while (m!= NULL)
  {
    if (m->num_unique_orgs > best_found->num_unique_orgs)
    {
      best_found = m;
    }
    else if (m->num_unique_orgs == best_found->num_unique_orgs
           && m->num_groups > best_found->num_groups)
    {
      best_found = m;
    }
    else if (m->num_unique_orgs == best_found->num_unique_orgs
           && m->num_groups == best_found->num_groups
           && m->num_mods < best_found->num_mods)
    {
      best_found = m;
    }
    m = m->next;
  }

  m = mc_list;
  while (m != NULL)
  {
    if (m != best_found)
    {
      newm = m->next;
      m->next = NULL;
      FreeModifierCombo (m);
      m = newm;
    }
    else
    {
      FreeModifierCombo (m->next);
      m->next = NULL;
      m = NULL;
    }
  }
  return best_found;
}


/* create combo with the specified modifiers */
NLM_EXTERN ValNodePtr GetModifierIndicesFromModList (
  ModifierItemLocalPtr modList
)
{
  Int4       feature_index;
  ValNodePtr modifier_indices = NULL;

  if (modList == NULL) return NULL;
  for (feature_index = 0; feature_index < numDefLineModifiers; feature_index++)
  {
    if (modList[feature_index].any_present && modList [feature_index].required)
    {
      ValNodeAddInt (&modifier_indices, 0, feature_index);
    }
  }
  return modifier_indices;
}


/* This is the callback function for sorting the modifier list.  It 
 * implements an order specified by the indexers.
 */
static Int4 DefLineQualPresentationOrder [] = {
  DEFLINE_POS_Transgenic,
  DEFLINE_POS_Strain,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Specimen_voucher,
  DEFLINE_POS_Ecotype,
  DEFLINE_POS_Type,
  DEFLINE_POS_Serotype,
  DEFLINE_POS_Authority,
  DEFLINE_POS_Breed
};

static int LIBCALLBACK SortByImportance (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr vnp1;
  ValNodePtr vnp2;
  Int4       num_defline_qual_sort_order, index;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.intvalue == vnp2->data.intvalue) return 0;

  num_defline_qual_sort_order = sizeof (DefLineQualPresentationOrder) / sizeof (Int4);
  for (index = 0; index < num_defline_qual_sort_order; index++)
  {
    if (vnp1->data.intvalue == DefLineQualPresentationOrder [ index ]) return -1;
    if (vnp2->data.intvalue == DefLineQualPresentationOrder [ index ]) return 1;
  }

  if ((vnp1->data.intvalue < 0 || vnp1->data.intvalue > numDefLineModifiers)
    && (vnp2->data.intvalue < 0 || vnp2->data.intvalue > numDefLineModifiers))
  {
    return 0;
  }
  if (vnp1->data.intvalue < 0 || vnp1->data.intvalue > numDefLineModifiers)
  {
    return 1;
  }
  if (vnp2->data.intvalue < 0 || vnp2->data.intvalue > numDefLineModifiers)
  {
    return -1;
  }

  if (DefLineModifiers [ vnp1->data.intvalue].isOrgMod
    && (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
      || vnp2->data.intvalue == DEFLINE_POS_Plasmid_name
      || vnp2->data.intvalue == DEFLINE_POS_Endogenous_virus_name))
  {
    return -1;
  }
  if (DefLineModifiers [ vnp2->data.intvalue].isOrgMod
    && (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
      || vnp1->data.intvalue == DEFLINE_POS_Plasmid_name
      || vnp1->data.intvalue == DEFLINE_POS_Endogenous_virus_name))
  {
    return 1;
  }
  
  if (vnp1->data.intvalue == DEFLINE_POS_Plasmid_name)
  {
    return -1;
  }
  if (vnp2->data.intvalue == DEFLINE_POS_Plasmid_name)
  {
    return 1;
  }
    
  if (vnp1->data.intvalue == DEFLINE_POS_Endogenous_virus_name)
  {
    return -1;
  }
  if (vnp2->data.intvalue == DEFLINE_POS_Endogenous_virus_name)
  {
    return 1;
  }

  if (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
     && vnp2->data.intvalue == DEFLINE_POS_Clone)
  {
    return 1;
  }
  if (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
     && vnp1->data.intvalue == DEFLINE_POS_Clone)
  {
    return -1;
  }

  if (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
     && vnp2->data.intvalue == DEFLINE_POS_Haplotype)
  {
    return 1;
  }
  if (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
     && vnp1->data.intvalue == DEFLINE_POS_Haplotype)
  {
    return -1;
  }

  if (vnp1->data.intvalue > vnp2->data.intvalue) return 1;
  if (vnp1->data.intvalue < vnp2->data.intvalue) return -1;
  return 0;
}

static Boolean RecordHasModifier (
  BioSourcePtr biop,
  Int4         modifier_index
)
{
  OrgModPtr     mod;
  OrgNamePtr    onp;
  SubSourcePtr  ssp;

  if (biop == NULL
    || modifier_index < 0
    || modifier_index >= numDefLineModifiers)
  {
    return FALSE;
  }
  if (DefLineModifiers[modifier_index].isOrgMod)
  {
    if (biop->org == NULL || (onp = biop->org->orgname) == NULL)
    {
      return FALSE;
    }
    mod = onp->mod;
    while (mod != NULL
        && mod->subtype != DefLineModifiers[modifier_index].subtype)
    {
      mod = mod->next;
    }
    if (mod != NULL && mod->subname != NULL)
    {
      return TRUE;
    }
  } else {
    ssp = biop->subtype;
    while (ssp != NULL && ssp->subtype != DefLineModifiers[modifier_index].subtype)
    {
      ssp = ssp->next;
    }
    if (ssp != NULL && ssp->name != NULL)
    {
      return TRUE;
    }
  }
  return FALSE;
}
  
/* This function adds in required modifiers for HIV sequences */
static void AddHIVModifierIndices (
  ValNodePtr PNTR modifier_indices,
  BioSourcePtr biop,
  ModifierItemLocalPtr modList,
  CharPtr taxName,
  Int4    clone_isolate_HIV_rule_num
)
{
  ValNodePtr  vnp;
  Boolean have_country_in_list;
  Boolean have_isolate_in_list;
  Boolean have_clone_in_list;
  Boolean have_country_mod;
  Boolean have_isolate_mod;
  Boolean have_clone_mod;

  /* special handling for HIV */
  if (StringNICmp (taxName, "HIV-1", 5) != 0
    && StringNICmp (taxName, "HIV-2", 5) != 0)
  {
    return;
  }

  have_country_in_list = FALSE;
  have_isolate_in_list = FALSE;
  have_clone_in_list = FALSE;
  have_country_mod = RecordHasModifier (biop, DEFLINE_POS_Country);
  have_isolate_mod = RecordHasModifier (biop, DEFLINE_POS_Isolate);
  have_clone_mod = RecordHasModifier (biop, DEFLINE_POS_Clone);

  if (modifier_indices != NULL)
  {
    for (vnp = *modifier_indices;
         vnp != NULL
           && (! have_country_in_list
             || ! have_isolate_in_list
             || ! have_clone_in_list);
         vnp = vnp->next)
    {
      if (vnp->data.intvalue == DEFLINE_POS_Country)
      {
        have_country_in_list = TRUE;
      }
      else if (vnp->data.intvalue == DEFLINE_POS_Isolate)
      {
        have_isolate_in_list = TRUE;
      }
      else if (vnp->data.intvalue == DEFLINE_POS_Clone)
      {
        have_clone_in_list = TRUE;
      }
    }
  }

  if ( ! have_country_in_list && have_country_mod)
  {
    vnp = ValNodeNew (*modifier_indices);
    vnp->data.intvalue = DEFLINE_POS_Country;
    if (*modifier_indices == NULL) *modifier_indices = vnp;
  }

  if ((have_clone_in_list && have_clone_mod)
      || (have_isolate_in_list && have_isolate_mod))
  {
    /* don't need HIV rule */
  }
  else
  {
    if ( ! have_isolate_in_list
        && have_isolate_mod
        && ( clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_prefer_isolate
          || clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_want_both
          || ! have_clone_mod))
    {
      vnp = ValNodeNew (*modifier_indices);
      vnp->data.intvalue = DEFLINE_POS_Isolate;
      if (*modifier_indices == NULL) *modifier_indices = vnp;
    }
    
    if ( ! have_clone_in_list
        && have_clone_mod
        && ( clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_prefer_clone
          || clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_want_both
          || ! have_isolate_mod))
    {
      vnp = ValNodeNew (*modifier_indices);
      vnp->data.intvalue = DEFLINE_POS_Clone;
      if (*modifier_indices == NULL) *modifier_indices = vnp;
    }
  }
}

/* This function looks for an OrgMod note that contains the phrase
 * "type strain of".  This function is used to determine whether
 * strain is a required modifier for the defline for this source.
 */
static Boolean HasTypeStrainComment (BioSourcePtr biop)
{
  OrgModPtr mod;
  
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL)
  {
    return FALSE;
  }
  
  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != ORGMOD_strain)
  {
    mod = mod->next;
  }
  
  if (mod == NULL)
  {
    return FALSE;
  }
  
  if (!UseOrgModifier (mod, biop->org->taxname))
  {
    return FALSE;
  }
  
  mod = biop->org->orgname->mod;
  while (mod != NULL)
  {
    if (mod->subtype == 255
        && StringISearch (mod->subname, "type strain of") != NULL)
    {
      return TRUE;
    }
    mod = mod->next;
  }
  return FALSE;
}


/* This function checks to see if there is a type strain comment on
 * the bio source.  If there is one, it checks to see whether strain
 * is already in the list of modifiers for the definition line.
 * If strain is not already in the list, it is added.
 */
static void 
AddTypeStrainModifierIndices 
(ValNodePtr PNTR modifier_indices,
 BioSourcePtr    biop)
{
  ValNodePtr vnp;
  
  if (modifier_indices == NULL || biop == NULL || ! HasTypeStrainComment (biop))
  {
    return;
  }

  for (vnp = *modifier_indices;
       vnp != NULL && vnp->data.intvalue != DEFLINE_POS_Strain;
       vnp = vnp->next)
  {
  }
  
  if (vnp == NULL)
  {
    ValNodeAddInt (modifier_indices, 0, DEFLINE_POS_Strain);
  }
}

static Boolean SpecialHandlingForSpecialTechniques (
  BioseqPtr bsp
);

/* This function checks to see if the Bioseq has a WGS technique.
 * If so, and if the strain text is not present in the taxname,
 * and strain is not already in the list of modifiers for the 
 * definition line, add strain.
 */
static void 
AddWGSModifierIndices 
(ValNodePtr PNTR modifier_indices,
 BioSourcePtr    biop,
 BioseqPtr       bsp)
{
  ValNodePtr vnp;
  OrgModPtr  omp;
  
  if (modifier_indices == NULL || biop == NULL 
      || biop->org == NULL
      || biop->org->orgname == NULL
      || biop->org->orgname->mod == NULL
      || ! SpecialHandlingForSpecialTechniques (bsp))
  {
    return;
  }

  for (vnp = *modifier_indices;
       vnp != NULL && vnp->data.intvalue != DEFLINE_POS_Strain;
       vnp = vnp->next)
  {
  }
  
  if (vnp == NULL)
  {
    omp = biop->org->orgname->mod;
    while (omp != NULL && omp->subtype != ORGMOD_strain) 
    {
      omp = omp->next;
    }
    if (omp != NULL) 
    {
      if (StringStr (biop->org->taxname, omp->subname) != NULL) 
      {
        /* don't add, present already */
      } else {
        /* add strain modifier */
        ValNodeAddInt (modifier_indices, 0, DEFLINE_POS_Strain);
      }
    }
  }
}

/* This function provides a label to be used in the definition line for
 * each modifier that requires one.  Most modifiers use a label that is 
 * similar to the name of the modifier displayed in the definition line
 * options dialog.
 */
NLM_EXTERN void AddModifierLabel (
  Boolean use_labels,
  Boolean is_orgmod,
  Uint1   subtype,
  CharPtr modifier_text
)
{
  CharPtr cp;
  if (!is_orgmod && subtype == SUBSRC_endogenous_virus_name)
  {
    StringCpy (modifier_text, "endogenous virus");
  }
  else if (is_orgmod && subtype == ORGMOD_specimen_voucher)
  {
    if (use_labels)
    {
      StringCpy (modifier_text, "voucher");
    } 
    else
    {
      modifier_text [0] = 0;
    }
  }
  else if (use_labels 
           || (!is_orgmod 
               && (subtype == SUBSRC_transgenic
                   || subtype == SUBSRC_plasmid_name)))
  {
    if (is_orgmod) 
    {
      StringCpy (modifier_text, GetOrgModQualName (subtype));
    } else {
      StringCpy (modifier_text, GetSubsourceQualName (subtype));
    }
    modifier_text[0] = tolower(modifier_text[0]);
    cp = StringStr (modifier_text, "-name");
    if (cp != NULL) *cp = 0;
  }
  else
  {
    modifier_text[0] = 0;
  }
}

typedef struct orgmodabbrevdata {
  Int2    subtype;
  CharPtr abbrev;
} OrgModAbbrevData, PNTR OrgModAbbrevPtr;

static OrgModAbbrevData orgmod_abbrevs[] = {
  { ORGMOD_variety, "var." },
  { ORGMOD_forma, "f." },
  { ORGMOD_forma_specialis, "f. sp." },
  { ORGMOD_pathovar, "pv." }
};

/* The UseOrgModifier function looks for the values of certain kinds of 
 * modifiers in the taxonomy name, so that they will not be added to the
 * definition line as modifiers if they are already present in the
 * taxonomy name.
 */
NLM_EXTERN Boolean UseOrgModifier (
  OrgModPtr mod,
  CharPtr   taxName
)
{
  CharPtr value_found;
  Int4    value_len;
  Int4    num_abbrevs, i;
  CharPtr abbrev_start;
  Boolean other_abbrev_found;

  if (mod == NULL || mod->subname == NULL) return FALSE;

  num_abbrevs = sizeof (orgmod_abbrevs) / sizeof (OrgModAbbrevData);

  /* If selected modifiers already appear in the tax Name, */
  /* don't use them in the organism description again */
  if (mod->subtype == ORGMOD_strain
    || mod->subtype == ORGMOD_variety
    || mod->subtype == ORGMOD_sub_species
    || mod->subtype == ORGMOD_forma
    || mod->subtype == ORGMOD_forma_specialis
    || mod->subtype == ORGMOD_pathovar
    || mod->subtype == ORGMOD_specimen_voucher)
  {
    value_found = StringStr (taxName, mod->subname);
    value_len = StringLen (mod->subname);
    while (value_found != NULL)
    {
      if (value_found == taxName)
      {
        value_found = StringStr (value_found + 1, mod->subname);
        continue;
      }
      if (*(value_found - 1) != ' ' && *(value_found - 1) != '(')
      {
        value_found = StringStr (value_found + 1, mod->subname);
        continue;
      }
      if (*(value_found - 1) == ')' && *(value_found + value_len) != ')')
      {
        value_found = StringStr (value_found + 1, mod->subname);
        continue;
      }
      if (*(value_found + value_len) != ' ' && *(value_found + value_len) != 0)
      {
        value_found = StringStr (value_found + 1, mod->subname);
        continue;
      }
      other_abbrev_found = FALSE;
      for (i = 0; i < num_abbrevs; i++)
      {
        abbrev_start = value_found - StringLen (orgmod_abbrevs[i].abbrev) - 1;
        if (abbrev_start > taxName
          && StringNCmp (abbrev_start,
                         orgmod_abbrevs[i].abbrev,
                         StringLen (orgmod_abbrevs[i].abbrev)) == 0)
        {
          if (mod->subtype == orgmod_abbrevs[i].subtype)
          {
            return FALSE;
          }
          else
          {
            other_abbrev_found = TRUE;
          }
        }
      }
      if ( ! other_abbrev_found 
        && ( mod->subtype == ORGMOD_strain
          || mod->subtype == ORGMOD_sub_species
          || mod->subtype == ORGMOD_specimen_voucher))
      {
        return FALSE;
      }
      value_found = StringStr (value_found + 1, mod->subname);
    }
  }
  return TRUE;
}

/* The SetRequiredModifiers function copies the default required values from
 * the global DefLineModifiers array into the local list of modifier
 * information.
 */
NLM_EXTERN void SetRequiredModifiers (
  ModifierItemLocalPtr modList
)
{
  Int4  item_index;

  for (item_index = 0; item_index < numDefLineModifiers; item_index++)
  {
    modList[item_index].required = IsDeflineModifierRequiredByDefault(DefLineModifiers[item_index].isOrgMod,
                                                                      DefLineModifiers[item_index].subtype);
  }
  
}

/* This function fixes HIV abbreviations, removes items in parentheses,
 * and trims spaces around the taxonomy name.
 */
NLM_EXTERN void CleanUpTaxName (
  CharPtr taxName,
  Boolean keep_in_paren
)
{
  CharPtr ptr;

  if (StringICmp (taxName, "Human immunodeficiency virus type 1") == 0
    || StringICmp (taxName, "Human immunodeficiency virus 1") == 0)
  {
    StringCpy (taxName, "HIV-1");
  }
  else if (StringICmp (taxName, "Human immunodeficiency virus type 2") == 0
    || StringICmp (taxName, "Human immunodeficiency virus 2") == 0)
  {
    StringCpy (taxName, "HIV-2");
  }
  else
  {
    if (! keep_in_paren)
    {
      ptr = StringStr (taxName, "(");
      if (ptr != NULL)
        *ptr = '\0';
    }
    TrimSpacesAroundString (taxName);
  }
}

/* This function gets the BioSource descriptor for the BioSeq. */
NLM_EXTERN BioSourcePtr GetBiopForBsp (
  BioseqPtr bsp
)
{
  SeqMgrDescContext  dcontext;
  SeqDescrPtr    sdp;
  BioSourcePtr    biop;

  if (bsp == NULL) return NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    return biop;
  }
  
  return NULL;
}


NLM_EXTERN Boolean IsSpName (CharPtr taxName)
{
  CharPtr cp;

  cp = StringStr (taxName, " sp.");
  /* check to make sure not "f. sp." */
  if (cp != NULL && cp[4] == ' '
      && (cp - taxName < 2 || *(cp - 2) != 'f' || *(cp - 1) != '.'))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static ValNodePtr ValNodeIntCopy (ValNodePtr orig)
{
  ValNodePtr cpy = NULL, last = NULL, vnp;

  while (orig != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = orig->choice;
    vnp->data.intvalue = orig->data.intvalue;
    if (last == NULL) {
      cpy = vnp;
    } else {
      last->next = vnp;
    }
    last = vnp;
    orig = orig->next;
  }
  return cpy;
}


NLM_EXTERN Boolean IsTSA (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mip;
  Boolean           rval = FALSE;

  if (bsp == NULL) return FALSE;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
       sdp != NULL && !rval;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_molinfo, &context)) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_tsa) {
      rval = TRUE;
    }
  }
  return rval;
}


NLM_EXTERN Boolean IsGenomeProjectIDDescriptor (SeqDescrPtr sdp) 
{
  UserObjectPtr        uop;
  ObjectIdPtr          oip;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return FALSE;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "GenomeProjectsDB") == 0) {
      return TRUE;
    }
  }
  return FALSE;
}


NLM_EXTERN SeqDescrPtr GetGenomeProjectIDDescriptor (BioseqPtr bsp)
{
  SeqDescrPtr sdp;

  if (bsp == NULL) return NULL;
  sdp = bsp->descr;
  while (sdp != NULL) {
    if (IsGenomeProjectIDDescriptor(sdp)) {
      return sdp;
    }
    sdp = sdp->next;
  }
  return NULL;
}



static Int4 GetGenomeProjectID (BioseqPtr bsp)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  UserObjectPtr     uop;
  UserFieldPtr      ufp;
  Int4              gpid = 0;

  if (bsp == NULL) return 0;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL && gpid == 0) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "GenomeProjectsDB") == 0)
    {
      ufp = uop->data;
      while (ufp != NULL && gpid == 0) {
        if (ufp->label != NULL 
            && StringCmp (ufp->label->str, "ProjectID") == 0
            && ufp->choice == 2) {
          gpid = ufp->data.intvalue;
        }
        ufp = ufp->next;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }

  return gpid;
}


/* This function generates a string describing the organism based on the
 * modifiers selected and other organism description options.
 */
static CharPtr GetOrganismDescription (
  BioseqPtr bsp, 
  ModifierItemLocalPtr modList,
  ValNodePtr   modifier_indices,
  OrganismDescriptionModifiersPtr odmp
)
{
  Char         taxName [196];
  Char         modifier_text [256];
  ValNodePtr   strings;
  BioSourcePtr biop;
  OrgModPtr    mod;
  SubSourcePtr ssp;
  ValNodePtr   vnp;
  Int2         feature_index;
  CharPtr      org_desc;
  CharPtr      cp;
  Uint4        no_semicolon_len, label_len;
  CharPtr      tmp;
  Char         id[255];
  Int4         gpid;
  SeqIdPtr     sip;
  DbtagPtr     dbtag;

  biop = NULL;
  strings = NULL;
  taxName [0] = '\0';

  biop = GetBiopForBsp (bsp);
  if (biop == NULL) return NULL;
  if (biop->org == NULL) return NULL;
  if (biop->org->taxname == NULL) return NULL;
  StringNCpy (taxName, biop->org->taxname, sizeof (taxName) - 1);
  taxName [ sizeof (taxName) - 1] = 0;

  CleanUpTaxName (taxName, odmp->keep_paren);

  if (odmp->exclude_sp && IsSpName(taxName))
  {
    return StringSave (taxName);
  }

  if (odmp->exclude_cf)
  {
    cp = StringStr (taxName, " cf.");
    if (cp != NULL)
    {
      return StringSave (taxName);
    }
  }
  
  if (odmp->exclude_aff)
  {
    cp = StringStr (taxName, " aff.");
    if (cp != NULL)
    {
      return StringSave (taxName);
    }
  }
  if (odmp->exclude_nr)
  {
    cp = StringStr (taxName, " nr.");
    if (cp != NULL)
    {
      return StringSave (taxName);
    }
  }
  
  if (biop->origin == ORG_MUT)
  {
    ValNodeAddStr (&strings, 0, StringSave ("Mutant"));
  }

  ValNodeAddStr (&strings, 0, StringSave (taxName));

  if (HasTypeStrainComment (biop))
  {
    
  } 

  /* copy modifier indices list */
  modifier_indices = ValNodeIntCopy (modifier_indices);
  AddHIVModifierIndices (&modifier_indices, biop, modList, taxName,
                         odmp->clone_isolate_HIV_rule_num);
  AddTypeStrainModifierIndices (&modifier_indices, biop);
  AddWGSModifierIndices (&modifier_indices, biop, bsp);
  
  modifier_indices = ValNodeSort (modifier_indices, SortByImportance);
  for (vnp = modifier_indices;
       vnp != NULL && (odmp->max_mods == -99 || odmp->max_mods > 0);
       vnp = vnp->next)
  {
    feature_index = vnp->data.intvalue;
    if (! odmp->use_modifiers && !IsDeflineModifierRequiredByDefault(DefLineModifiers[feature_index].isOrgMod,
                                                                     DefLineModifiers[feature_index].subtype))
    {
      /* do nothing */
    }
    else if (DefLineModifiers[feature_index].isOrgMod)
    {
      if (biop->org == NULL || biop->org->orgname == NULL) continue;
      mod = biop->org->orgname->mod;
      while (mod != NULL
        && mod->subtype != DefLineModifiers[feature_index].subtype)
      {
        mod = mod->next;
      }
      if ( UseOrgModifier (mod, taxName))
      {
        if (odmp->allow_semicolon_in_modifier) {
          no_semicolon_len = StringLen (mod->subname);
        } else {
          no_semicolon_len = StringCSpn (mod->subname, ";");
        }

        if (mod->subtype == ORGMOD_nat_host)
        {
          sprintf (modifier_text, "from ");
          if (no_semicolon_len > sizeof (modifier_text) - 6)
          {
            no_semicolon_len = sizeof (modifier_text) - 6;
          }
          StringNCpy (modifier_text + 5, mod->subname,
                      no_semicolon_len);
          modifier_text[no_semicolon_len + 5] = 0;
        }
        else
        {
          AddModifierLabel (odmp->use_labels, TRUE, mod->subtype, modifier_text);
          if (modifier_text[0] != 0)
            StringCat (modifier_text, " ");
          label_len = StringLen (modifier_text);
          if (no_semicolon_len > (Int4) sizeof (modifier_text) - label_len - 1)
          {
            no_semicolon_len = (Int4) sizeof (modifier_text) - label_len - 1;
          } 
          if (mod->subtype == ORGMOD_specimen_voucher && StringNICmp (mod->subname, "personal:", 9) == 0)
          {
            tmp = mod->subname + 9;
            while (isspace (*tmp)) 
            {
              tmp++;
            }
            if (odmp->allow_semicolon_in_modifier) {
              no_semicolon_len = StringLen (tmp);
            } else {
              no_semicolon_len = StringCSpn (tmp, ";");
            }
          }
          else
          {
            tmp = mod->subname;
          }

          StringNCat (modifier_text, tmp,
                      no_semicolon_len);
          modifier_text [ no_semicolon_len + label_len] = 0;
        }
        ValNodeCopyStr( &strings, 0, modifier_text);
        if (odmp->max_mods != -99)
          odmp->max_mods --;
      }
    } else {
      ssp = biop->subtype;
      while (ssp != NULL
          && ssp->subtype != DefLineModifiers[feature_index].subtype)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL)
      {
        if (odmp->include_country_extra || odmp->allow_semicolon_in_modifier)
        {
          no_semicolon_len = StringLen (ssp->name);
        }
        else
        {
          no_semicolon_len = StringCSpn (ssp->name, ";");
        }
        AddModifierLabel (odmp->use_labels, FALSE, ssp->subtype, modifier_text);
        if (ssp->subtype == SUBSRC_transgenic)
        {
          /* do nothing, transgenic already captured from label */
        }
        else if (ssp->subtype == SUBSRC_country)
        {
          sprintf (modifier_text, "from ");
          if (no_semicolon_len > sizeof (modifier_text) - 6)
          {
            no_semicolon_len = sizeof (modifier_text) - 6;
          }
          StringNCpy (modifier_text + 5, ssp->name, no_semicolon_len);
          modifier_text[5 + no_semicolon_len] = 0;
          if (!odmp->include_country_extra)
          {
            cp = StringChr (modifier_text, ':');
            if (cp != NULL) *cp = 0;
          }
        }
        else if (ssp->name != NULL && ssp->name[0] != 0
          && (ssp->subtype != SUBSRC_plasmid_name
            || StringCmp (ssp->name, "unnamed") != 0))
        {
          if (modifier_text[0] != 0)
            StringCat (modifier_text, " ");
          label_len = StringLen (modifier_text);
          if (no_semicolon_len > sizeof (modifier_text) - 1 - label_len)
          {
            no_semicolon_len = sizeof (modifier_text) - 1 - label_len;
          }
          StringNCat (modifier_text, ssp->name, no_semicolon_len);
          modifier_text [ no_semicolon_len + label_len ] = 0;
        }
        ValNodeCopyStr( &strings, 0, modifier_text);
        if (odmp->max_mods != -99)
          odmp->max_mods --;
      }
    }
  }

  /* add TSA project ID if necessary */
  if (IsTSA (bsp)) {
    gpid = GetGenomeProjectID (bsp);
    if (gpid > 0) {
      sprintf (id, "gpid:%d", gpid);
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GENERAL && sip->data.ptrvalue != NULL) {
          dbtag = (DbtagPtr) sip->data.ptrvalue;
          if (StringCmp (dbtag->db, id) == 0 && dbtag->tag != NULL) {
            if (dbtag->tag->str != NULL) {
              ValNodeAddPointer (&strings, 0, StringSave (dbtag->tag->str));
            } else {
              sprintf (id, "%d", dbtag->tag->id);
              ValNodeAddPointer (&strings, 0, StringSave (id));
            }
            break;
          }
        }
      }
    }
  }
  
  org_desc = MergeValNodeStrings (strings, FALSE);
  ValNodeFreeData (strings);
  modifier_indices = ValNodeFree (modifier_indices);
  return org_desc;

}

/* end of organism description section */

/* This section of code contains functions which are useful for dealing
 * with locations of features (SeqLocPtr objects).
 */

/* This function determines whether location A is on the same strand as
 * location B
 */
static Boolean AreAAndBOnSameStrand (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)
{
  Uint1 strand1;
  Uint2 strand2;

  strand1 = SeqLocStrand (slp1);
  strand2 = SeqLocStrand (slp2);
  if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus)
    return FALSE;
  else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus)
    return FALSE;
  else
    return TRUE;
}

/* This function determines whether location A is contained in or equal to
 * location B and on the same strand as location B.
 */
NLM_EXTERN Boolean IsLocAInBonSameStrand (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)
{
  if (! AreAAndBOnSameStrand ( slp1, slp2))
  {
    return FALSE;
  }
  else if ( SeqLocAinB (slp1, slp2) < 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

/* This function calculates the intersection between two locations.
 */
static SeqLocPtr SeqLocIntersection (
  SeqLocPtr slp1,
  SeqLocPtr slp2,
  BioseqPtr bsp
)
{
  SeqLocPtr diff1, diff2, result;

  diff1 = SeqLocMerge ( bsp, slp1, NULL, FALSE, TRUE, FALSE);
  diff1 = SeqLocSubtract (diff1, slp2);
  diff2 = SeqLocMerge ( bsp, slp2, NULL, FALSE, TRUE, FALSE);
  diff2 = SeqLocSubtract (diff2, slp1);
  result = SeqLocMerge ( bsp, slp1, slp2, FALSE, TRUE, FALSE);
  
  if (diff1 != NULL)
  {
    result = SeqLocSubtract (result, diff1);
    SeqLocFree (diff1);
    if (result == NULL) return NULL;
  }
  if (diff2 != NULL)
  {
    result = SeqLocSubtract (result, diff2);
    SeqLocFree (diff2);
    if (result == NULL) return NULL;
  }
  return result;
}

#define ADJACENT_TYPE_ANY        0
#define ADJACENT_TYPE_UPSTREAM   1
#define ADJACENT_TYPE_DOWNSTREAM 2

/* This function determines whether A is "next to" B and upstream or downstream
 * from B.  A cannot overlap B.  If allow_interval is TRUE, there can be
 * space between A and B.
 */
static Boolean IsAAdjacentToB (
  SeqLocPtr a,
  SeqLocPtr b,
  BioseqPtr bsp,
  Int2      adjacent_type,
  Boolean   allow_interval
)
{
  Int4      a_end, b_end;
  Uint2     strand;
 
  if (adjacent_type != ADJACENT_TYPE_ANY
    && adjacent_type != ADJACENT_TYPE_UPSTREAM
    && adjacent_type != ADJACENT_TYPE_DOWNSTREAM)
  {
    return FALSE;
  }
  
  if ( ! AreAAndBOnSameStrand (a, b))
  {
    return FALSE;
  }

  strand = SeqLocStrand (a);
  if ( adjacent_type == ADJACENT_TYPE_ANY)
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_LEFT_END);
    if ((allow_interval && b_end < a_end)
      || b_end == a_end + 1)
    {
      return TRUE;
    }
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_RIGHT_END);
    if ((allow_interval && b_end > a_end)
      || a_end == b_end + 1)
    {
      return TRUE;
    }
    return FALSE;
  }
  else if ( (strand == Seq_strand_minus
      && adjacent_type == ADJACENT_TYPE_UPSTREAM)
    || (strand != Seq_strand_minus
      && adjacent_type == ADJACENT_TYPE_DOWNSTREAM))
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_LEFT_END);
    if ((allow_interval && b_end < a_end)
      || b_end == a_end + 1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  else
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_RIGHT_END);
    if ((allow_interval && b_end > a_end)
      || a_end == b_end + 1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
}

static Boolean IsAEmptyIntervalOfB (SeqLocPtr a, SeqLocPtr b, BioseqPtr bsp)
{
  Int4 a_right, a_left, b_right, b_left, prev_right, prev_left;
  SeqLocPtr slp;
  Uint1 a_strand, b_strand;
  
  if (a == NULL || b == NULL || bsp == NULL) return FALSE;
  
  a_strand = SeqLocStrand (a);
  b_strand = SeqLocStrand (b);
  if ((a_strand == Seq_strand_minus && b_strand != Seq_strand_minus)
      || (a_strand != Seq_strand_minus && b_strand == Seq_strand_minus)) {
      return FALSE;
  }
  
  a_right = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
  a_left = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
  
  slp = SeqLocFindNext (b, NULL);
  prev_right = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
  prev_left = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
  slp = SeqLocFindNext (b, slp);
  while (slp != NULL) {
    b_right = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
    b_left = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
    if (a_left == prev_right + 1 && a_right == b_left - 1) {
      return TRUE;
    } else if (a_left == b_right + 1 && a_right == prev_left - 1) {
      return TRUE;
    } else {
      prev_right = b_right;
      prev_left = b_left;
      slp = SeqLocFindNext (b, slp);
    }
  }
  return FALSE;
}


static Boolean LocAContainsIntervalOfB (SeqLocPtr a, SeqLocPtr b)
{
  SeqLocPtr interval;
  Boolean   rval = FALSE;

  if (a == NULL || b == NULL) return FALSE;

  interval = SeqLocFindNext (b, NULL);
  while (interval != NULL && !rval) {
    if (IsLocAInBonSameStrand (interval, a)) {
      rval = TRUE;
    } else {
      interval = SeqLocFindNext (b, interval);
    }
  }
  return rval;
}


/* This section of code deals with identifying and labeling features
 * for the definition line.
 * The features currently handled are:
 *     genes
 *     exons
 *     introns
 *     LTRs
 *     3' UTRs
 *     5' UTRs
 *     CDSs
 *     rRNA
 *     mRNA
 *     misc RNA
 *     snRNA
 *     snoRNA
 *     insertion sequences
 *     integrons
 *     D-loops
 *     mRNA
 *     tRNA
 *     control regions
 *     misc feature listed as intergenic spacer in comment
 *     satellite sequences
 *     promoter regions
 *     endogenous virus source features
 *     transposons
 */

static Boolean LIBCALLBACK IsGene (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return FALSE;
  return TRUE;
}

static CharPtr GetGeneName (GeneRefPtr grp, Boolean suppress_locus_tag)
{
  ValNodePtr syn;

  if (grp == NULL) return NULL;
  if (SeqMgrGeneIsSuppressed (grp)) return NULL;
  if (StringDoesHaveText (grp->locus)) return grp->locus;
  if (! suppress_locus_tag && StringDoesHaveText (grp->locus_tag)) 
      return grp->locus_tag;
  if (StringDoesHaveText (grp->desc)) return grp->desc;
  for (syn = grp->syn; syn != NULL; syn = syn->next)
  {
    if (syn != NULL && syn->data.ptrvalue != NULL)
      return syn->data.ptrvalue;
  }
  return NULL;
}

static CharPtr GetAlleleName (GeneRefPtr grp, Boolean suppress_locus_tag)
{
  size_t  lenallele;
  size_t  lengenename;
  CharPtr  gene_name;
  CharPtr  buffer;

  if (grp == NULL) return NULL;
  if (StringHasNoText (grp->allele)) return NULL;
  gene_name = GetGeneName (grp, suppress_locus_tag);
  if (StringHasNoText (gene_name)) return NULL;
  lenallele = StringLen (grp->allele);
  lengenename = StringLen (gene_name);
  
  if (lenallele > lengenename
    && StringNICmp (gene_name, grp->allele, lengenename) == 0)
  {
    return StringSave (grp->allele);
  }
  else if (grp->allele[0] == '-')
  {
    buffer = MemNew (lenallele + lengenename + 1);
    if (buffer == NULL) return NULL;
    StringCpy (buffer, gene_name);
    StringCat (buffer, grp->allele);
  }
  else
  {
    buffer = MemNew (lenallele + lengenename + 2);
    if (buffer == NULL) return NULL;
    StringCpy (buffer, gene_name);
    StringCat (buffer, "-");
    StringCat (buffer, grp->allele);
  }

  return buffer;
}

/* This function compares the gene names and allele names of the gene
 * to see if they match.
 */
static Boolean DoGenesMatch 
(GeneRefPtr grp1,
 GeneRefPtr grp2,
 Boolean suppress_locus_tag)
{
  CharPtr name1;
  CharPtr name2;

  name1 = GetGeneName (grp1, suppress_locus_tag);
  name2 = GetGeneName (grp2, suppress_locus_tag);
  if (StringCmp (name1, name2) != 0) return FALSE;
  
  name1 = GetAlleleName (grp1, suppress_locus_tag);
  name2 = GetAlleleName (grp2, suppress_locus_tag);
  if ((name1 == NULL && name2 != NULL)
    || (name1 != NULL && name2 == NULL))
  {
    if (name1 != NULL) MemFree (name1);
    if (name2 != NULL) MemFree (name2);
    return FALSE;
  }

  if ((name1 == NULL && name2 == NULL)
           || (StringCmp (name1, name2) == 0))
  {
    if (name1 != NULL) MemFree (name1);
    if (name2 != NULL) MemFree (name2);
    return TRUE;
  }
  
  if (name1 != NULL) MemFree (name1);
  if (name2 != NULL) MemFree (name2);
  return  FALSE;
}

/* This function looks at the pseudo flag on the object itself as well as
 * the pseudo flag on the gene reference for the object (if one is present).
 */
NLM_EXTERN Boolean IsPseudo (
  SeqFeatPtr sfp
)
{
  GeneRefPtr grp;
  SeqMgrFeatContext context;

  if (sfp == NULL) return FALSE;
  if (sfp->pseudo) return TRUE;
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = sfp->data.value.ptrvalue;
  }
  else
  { 
    grp = SeqMgrGetGeneXref (sfp);
  }
  if (grp == NULL) 
  {
    if (sfp->data.choice != SEQFEAT_GENE) {
      sfp = SeqMgrGetOverlappingGene(sfp->location, &context);
      return IsPseudo(sfp);
    } else {
      return FALSE;
    }
  } else {
    return grp->pseudo;
  }
}

static Boolean LIBCALLBACK IsExon (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_exon) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsIntron (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_intron) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsExonOrIntron (SeqFeatPtr sfp)
{
  return IsExon(sfp) || IsIntron(sfp);
}

static Boolean LIBCALLBACK IsLTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_LTR) return FALSE;
  return TRUE;
}

static CharPtr GetLTRDescription (
  SeqFeatPtr sfp
)
{
  CharPtr description;
  size_t comment_len;
  if (sfp == NULL) return NULL;
  if (sfp->comment == NULL) return NULL;
  comment_len = StringLen (sfp->comment);
  if (comment_len > 3 && StringCmp (sfp->comment + comment_len - 3, "LTR") == 0)
  {
    description = (CharPtr) MemNew (comment_len - 3);
    if (description == NULL) return NULL;
    StringNCpy (description, sfp->comment, comment_len - 4);
    description[comment_len - 4] = 0;
  }
  else
  {
    description = StringSave (sfp->comment);
  }
  return description;
}

static Boolean LIBCALLBACK Is3UTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_3UTR) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK Is5UTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_5UTR) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsCDS (SeqFeatPtr sfp)
{
  if (sfp == NULL) return FALSE;
  if (sfp->data.choice == SEQFEAT_CDREGION)
    return TRUE;
  return FALSE;
}

static Boolean LIBCALLBACK IsrRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_rRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsMiscRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL 
    || (sfp->idx.subtype != FEATDEF_misc_RNA
      && sfp->idx.subtype != FEATDEF_otherRNA)) 
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK IsncRNA (
  SeqFeatPtr sfp
)
{
  RnaRefPtr rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return FALSE;
  if (sfp->idx.subtype == FEATDEF_scRNA 
      || sfp->idx.subtype == FEATDEF_snRNA 
      || sfp->idx.subtype == FEATDEF_snoRNA
      || sfp->idx.subtype == FEATDEF_ncRNA
      || sfp->idx.subtype == FEATDEF_tmRNA)
  {
    return TRUE;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type != 255 || rrp->ext.choice != 1)
  {
    return FALSE;
  }
  else if (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0
      || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
      || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0
      || IsStringInNcRNAClassList (rrp->ext.value.ptrvalue))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static CharPtr GetncRNAProduct (SeqFeatPtr sfp, Boolean use_ncrna_note)
{
  GBQualPtr gbq = NULL, q_class = NULL, q_product = NULL;
  CharPtr product = NULL;
  CharPtr tmp_class = NULL, cp;
  RnaRefPtr rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA)
  {
    return NULL;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  
  gbq = sfp->qual;
  while (gbq != NULL && q_product == NULL) {
    if (StringICmp (gbq->qual, "ncRNA_class") == 0
        && !StringHasNoText (gbq->val)) {
      q_class = gbq;
    } else if (StringICmp (gbq->qual, "product") == 0
        && !StringHasNoText (gbq->val)) {
      q_product = gbq;
    }
    gbq = gbq->next;
  }
  
  if (q_class != NULL) {
    tmp_class = StringSave (q_class->val);
    cp = tmp_class;
    while (*cp != 0) {
      if (*cp == '_') {
        *cp = ' ';
      }
      cp++;
    }
  }
  if (q_product != NULL) {
    if (tmp_class == NULL || StringStr (q_product->val, tmp_class) != NULL || StringCmp (tmp_class, "other") == 0) {
      product = StringSave (q_product->val);
    } else {
      product = (CharPtr) MemNew (sizeof (Char) * (StringLen (q_product->val) + StringLen (tmp_class) + 2));
      sprintf (product, "%s %s", q_product->val, tmp_class);
    }
  } else if (q_class != NULL) {
    if (use_ncrna_note && !StringHasNoText (sfp->comment)) {
      product = StringSave (sfp->comment);
    } else if (StringCmp (tmp_class, "other") == 0) {
      product = StringSave ("non-coding RNA");
    } else {
      product = StringSave (tmp_class);
    }
  } else if ((use_ncrna_note || 
              (rrp != NULL && rrp->type == 255
               && rrp->ext.choice == 1
               && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))
             && !StringHasNoText (sfp->comment)) {
    product = StringSave (sfp->comment);
  } else {
    product = StringSave ("non-coding RNA");
  }
  tmp_class = MemFree (tmp_class);
  cp = StringChr (product, ';');
  if (cp != NULL) {
    *cp = 0;
  }
  return product;
}


static Boolean LIBCALLBACK IsPrecursorRNA (SeqFeatPtr sfp)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_preRNA) return FALSE;
  return TRUE;
}


static CharPtr mobile_element_keywords [] = {
  "insertion sequence",
  "retrotransposon",
  "non-LTR retrotransposon",
  "transposon",
  "integron",
  "other",
  "SINE",
  "MITE",
  "LINE"
};

enum mobile_element_keyword_nums 
{
  eMobileElementInsertionSequence = 0,
  eMobileElementRetrotransposon,
  eMobileElementNonLTRRetrotransposon,
  eMobileElementTransposon,
  eMobileElementIntegron,
  eMobileElementOther,
  eMobileElementSINE,
  eMobileElementMITE,
  eMobileElementLINE
};

static Int4 StartsWithMobileElementKeyword (CharPtr txt)
{
  Int4 i, keyword_len;
    
  for (i=0; i < sizeof (mobile_element_keywords) / sizeof (CharPtr); i++) {
    keyword_len = StringLen (mobile_element_keywords[i]);
    if (StringNCmp (txt, mobile_element_keywords[i], keyword_len) == 0
        && (*(txt + keyword_len) == ':' || *(txt + keyword_len) == 0)) {      
      return i;
    }
  }
  return -1;
}

static Int4 IsMobileElementGBQual (GBQualPtr gbqual)
{
  Int4 keyword_idx;
  if (gbqual == NULL || gbqual->qual == NULL || gbqual->val == NULL) return -1;
  if (StringCmp (gbqual->qual, "mobile_element") != 0) return -1;
  keyword_idx = StartsWithMobileElementKeyword (gbqual->val);
  if (keyword_idx < 0) return -1;
  if (keyword_idx == eMobileElementOther
      && StringStr (gbqual->val, "transposable element") == NULL) {
    return -1;
  } else {
    return keyword_idx;
  }
}


static Boolean FeatureDoesNotGetPartialComplete (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  Int4 keyword_idx;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    keyword_idx = IsMobileElementGBQual(gbqual);
    if (keyword_idx == eMobileElementSINE
        || keyword_idx == eMobileElementLINE) {
      return TRUE;
    }
  }
  return FALSE;
}


NLM_EXTERN Boolean LIBCALLBACK IsMobileElement (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (IsMobileElementGBQual(gbqual) > -1) {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK IsRemovableMobileElement (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  Int4 keyword_idx;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    keyword_idx = IsMobileElementGBQual(gbqual);
    if (keyword_idx >= eMobileElementSINE
        && keyword_idx <= eMobileElementLINE) {
      return TRUE;
    }
  }
  return FALSE;
}

static CharPtr GetMobileElementTypeword (CharPtr desc_start, Int4 keyword_idx)
{
  if (keyword_idx < 0) return NULL;
  if (StringHasNoText (desc_start)) {
    return mobile_element_keywords[keyword_idx];
  }
  switch (keyword_idx) {
    case eMobileElementTransposon:
      if (StringStr (desc_start, "P-element") != NULL) {
        return "P-element";
      } else if (StringStr (desc_start, "MITE") != NULL) {
        return "MITE";
      } else {
        return mobile_element_keywords[keyword_idx];
      }
      break;
    case eMobileElementOther:
      return "transposable element";
      break;
    case eMobileElementIntegron:
      if (StringStr (desc_start, "superintegron") != NULL) {
        return "superintegron";
      } else {
        return mobile_element_keywords[keyword_idx];
      }
      break;
    default:
      return mobile_element_keywords[keyword_idx];
      break;
  }      
}


static void LIBCALLBACK GetMobileElementFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  GBQualPtr  gbqual;
  Int4    keyword_idx = -1;
  Int4    keyword_len;
  Int4    val_len;
  SeqFeatPtr      sfp;
  CharPtr         desc_start = NULL, typeword, cp;

  flp->pluralizable = TRUE;
  flp->is_typeword_first = FALSE;
  flp->typeword = NULL;
  flp->description = NULL;

  if (featlist == NULL) return;
  sfp = featlist->data.ptrvalue;
  if (sfp == NULL) return;

  gbqual = sfp->qual;
  while (gbqual != NULL 
         && (keyword_idx = IsMobileElementGBQual(gbqual)) < 0)
  {
    gbqual = gbqual->next;
  }
  if (gbqual == NULL) return;
    
  keyword_len = StringLen (mobile_element_keywords[keyword_idx]);
  desc_start = gbqual->val + keyword_len;
  while (isspace (*desc_start) || *desc_start == ':') {
    desc_start++;
  }

  /* find alternate typewords */
  typeword = GetMobileElementTypeword(desc_start, keyword_idx);
  if (typeword == NULL) return;
  keyword_len = StringLen (typeword);

  flp->typeword = StringSave (typeword);
  val_len = StringLen (desc_start);
  
  if (StringHasNoText (desc_start))
  {
    flp->is_typeword_first = FALSE;
    flp->description = NULL;
  } else if (StringCmp (desc_start, typeword) == 0) {
    /* just the keyword */
    flp->is_typeword_first = FALSE;
    flp->description = NULL;
    return;
  } else if (StringNCmp (desc_start, typeword, keyword_len) == 0) {
    /* starts with keyword */
    /* if keyword is hyphenated portion of name, no pluralization */
    if (desc_start[keyword_len] == '-') {
      flp->description = StringSave (desc_start);
      flp->typeword = MemFree (flp->typeword);
      flp->typeword = StringSave ("");
      flp->pluralizable = FALSE;
    } else {
      flp->is_typeword_first = TRUE;
      flp->description = StringSave (desc_start + keyword_len + 1);
    }
    return;
  } else if (val_len > 8 && StringCmp (desc_start + val_len - keyword_len, typeword) == 0
             && val_len - keyword_len - 1 >= 0
             && isspace (*(desc_start + val_len - keyword_len - 1))) {
    /* ends with keyword */
    flp->is_typeword_first = FALSE;
    flp->description = MemNew (val_len - keyword_len);
    if (flp->description == NULL) return;
    StringNCpy (flp->description, desc_start, val_len - keyword_len - 1);
    flp->description[val_len - keyword_len -1] = 0;
  } else if ((cp = StringStr (desc_start, typeword)) != NULL
             && cp != desc_start
             && isspace (*(cp -1))) {
    /* keyword in the middle */
    flp->description = StringSave (desc_start);
    flp->typeword = MemFree (flp->typeword);
    flp->typeword = StringSave ("");
    flp->pluralizable = FALSE;
  } else {
    /* keyword not in description */
    if (StringICmp (flp->typeword, "integron") == 0) {
      flp->is_typeword_first = FALSE;
    } else {
      flp->is_typeword_first = TRUE;
    }
    flp->description = StringSave (desc_start);
    if (StringCmp (flp->description, "") == 0) {
      flp->is_typeword_first = FALSE;
    }
  }
  if (StringCmp (flp->description, "unnamed") == 0) {
    flp->description = MemFree (flp->description);
  }
}


static Boolean LIBCALLBACK IsEndogenousVirusSequence (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (StringCmp (gbqual->qual, "endogenous_virus") == 0)
      return TRUE;
  }
  return FALSE;
}

static CharPtr GetEndogenousVirusSequenceDescription (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;

  if (sfp == NULL) return NULL;
  
  gbqual = sfp->qual;
  while (gbqual != NULL && StringCmp (gbqual->qual, "endogenous_virus") != 0)
  {
    gbqual = gbqual->next;
  }
  if (gbqual != NULL)
  {
    if (StringDoesHaveText (gbqual->val)
      && StringCmp (gbqual->val, "unnamed") != 0)
    {
      return StringSave (gbqual->val);
    }
  }
  return NULL;
}

static Boolean LIBCALLBACK IsDloop (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_D_loop) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsmRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IstRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsControlRegion (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringNCmp (sfp->comment, "control region", StringLen ("control region")) != 0)
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK IsGeneCluster (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || (StringStr (sfp->comment, "gene cluster") == NULL
        && StringStr (sfp->comment, "gene locus") == NULL))
  {
    return FALSE;
  }
  return TRUE;
}


static void LIBCALLBACK GetGeneClusterFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    cp;
  Int4       datalen;

  if (featlist == NULL || featlist->data.ptrvalue == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (StringHasNoText (main_feat->comment)) return;
  cp = StringStr (main_feat->comment, "gene cluster");
  if (cp == NULL)
  {
    cp = StringStr (main_feat->comment, "gene locus");
    if (cp == NULL) return;
    flp->typeword = StringSave ("gene locus");
  }
  else
  {
    flp->typeword = StringSave ("gene cluster");
  }
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;
  datalen = cp - main_feat->comment;
  if (datalen > 0)
  {
    flp->description = (CharPtr) MemNew ((datalen + 1) * sizeof (Char));
    StringNCpy (flp->description, main_feat->comment, datalen);
    flp->description [datalen] = 0;
    TrimSpacesAroundString (flp->description);
  }
  else
  {
    flp->description = NULL;
  }
}


static Boolean LIBCALLBACK IsIntergenicSpacer (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringStr (sfp->comment, "intergenic spacer") == NULL)
  {
    return FALSE;
  }
  return TRUE;
}

static ValNodePtr GettRNAGenesAndSpacers (CharPtr str);
static ValNodePtr FreetRNAGenesAndSpacersList (ValNodePtr list);

static Boolean LIBCALLBACK IsParsableList (
  SeqFeatPtr sfp
)
{
  ValNodePtr list;

  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL)
  {
    return FALSE;
  }
  
  list = GettRNAGenesAndSpacers (sfp->comment);
  if (list == NULL) 
  {
    return FALSE;
  } 
  else
  {
    FreetRNAGenesAndSpacersList (list);
    return TRUE;
  }
}


/* This function produces the default definition line label for a misc_feature 
 * that has the word "intergenic spacer" in the comment.  If the comment starts
 * with the word "contains", "contains" is ignored.  If "intergenic spacer"
 * appears first in the comment (or first after the word "contains", the text
 * after the words "intergenic spacer" but before the first semicolon (if any)
 * appear after the words "intergenic spacer" in the definition line.  If there
 * are words after "contains" or at the beginning of the comment before the words
 * "intergenic spacer", this text will appear in the definition line before the words
 * "intergenic spacer".
 */
static void LIBCALLBACK GetIntergenicSpacerFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    cp, buffer;
  Int4       datalen, offset;

  if (featlist == NULL || featlist->data.ptrvalue == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (StringHasNoText (main_feat->comment)) return;
  if (StringNCmp (main_feat->comment, "contains ", 9) == 0)
  {
    buffer = main_feat->comment + 9;
  }
  else if (StringNCmp (main_feat->comment, "may contain ", 12) == 0)
  {
    buffer = main_feat->comment + 12;
  }
  else
  {
    buffer = main_feat->comment;
  }
  cp = StringStr (buffer, "intergenic spacer");
  if (cp == NULL) return;
  flp->typeword = StringSave ("intergenic spacer");
  flp->pluralizable = FALSE;
  if (cp == buffer)
  {
    flp->is_typeword_first = TRUE;
    offset = StringLen ("intergenic spacer") + 1;
    if (StringNCmp (cp + offset, "and ", 4) == 0
        || *(cp + StringLen("intergenic spacer")) == 0)
    {
      flp->description = NULL;
    }
    else
    {
      flp->description = StringSave (cp + StringLen ("intergenic spacer") + 1);
      cp = StringChr (flp->description, ';');
      if (cp != NULL)
      {
        *cp = 0;
      }
    }
  }
  else
  {
    flp->is_typeword_first = FALSE;
    datalen = cp - buffer;
    flp->description = MemNew ( datalen + 1);
    if (flp->description == NULL) return;
    StringNCpy (flp->description, buffer, datalen);
    flp->description [datalen] = 0;
    TrimSpacesAroundString (flp->description);
  }
}

/* These structures are used for parsing tRNA and intergenic spacer information
 * from misc_feature comments. 
 */
typedef struct commentfeat
{
  CharPtr product_name;
  CharPtr gene_name;
} CommentFeatData, PNTR CommentFeatPtr;


static CommentFeatPtr CommentFeatFree (CommentFeatPtr cfp)
{
  if (cfp != NULL) {
    cfp->product_name = MemFree (cfp->product_name);
    cfp->gene_name = MemFree (cfp->gene_name);
    cfp = MemFree (cfp);
  }
  return cfp;
}


typedef struct intergenicspacerdef
{
  CharPtr first_gene;
  CharPtr second_gene;
} IntergenicSpacerDefData, PNTR IntergenicSpacerDefPtr;

static IntergenicSpacerDefPtr IntergenicSpacerDefFree (IntergenicSpacerDefPtr ip)
{
  if (ip != NULL) {
    ip->first_gene = MemFree (ip->first_gene);
    ip->second_gene = MemFree (ip->second_gene);
    ip = MemFree (ip);
  }
  return ip;
}


CharPtr s_tRNAGeneFromProduct (CharPtr product)
{
    CharPtr gene = NULL;

    if (StringNCmp (product, "tRNA-", 5) != 0) {
      return NULL;
    }
    product += 5;

    if (StringICmp (product, "Ala") == 0) {
        gene = "trnA";
    } else if (StringICmp (product, "Asx") == 0) {
        gene = "trnB";
    } else if (StringICmp (product, "Cys") == 0) {
        gene = "trnC";
    } else if (StringICmp (product, "Asp") == 0) {
        gene = "trnD";
    } else if (StringICmp (product, "Glu") == 0) {
        gene = "trnE";
    } else if (StringICmp (product, "Phe") == 0) {
        gene = "trnF";
    } else if (StringICmp (product, "Gly") == 0) {
        gene = "trnG";
    } else if (StringICmp (product, "His") == 0) {
        gene = "trnH";
    } else if (StringICmp (product, "Ile") == 0) {
        gene = "trnI";
    } else if (StringICmp (product, "Xle") == 0) {
        gene = "trnJ";
    } else if (StringICmp (product, "Lys") == 0) {
        gene = "trnK";
    } else if (StringICmp (product, "Leu") == 0) {
        gene = "trnL";
    } else if (StringICmp (product, "Met") == 0) {
        gene = "trnM";
    } else if (StringICmp (product, "Asn") == 0) {
        gene = "trnN";
    } else if (StringICmp (product, "Pyl") == 0) {
        gene = "trnO";
    } else if (StringICmp (product, "Pro") == 0) {
        gene = "trnP";
    } else if (StringICmp (product, "Gln") == 0) {
        gene = "trnQ";
    } else if (StringICmp (product, "Arg") == 0) {
        gene = "trnR";
    } else if (StringICmp (product, "Ser") == 0) {
        gene = "trnS";
    } else if (StringICmp (product, "Thr") == 0) {
        gene = "trnT";
    } else if (StringICmp (product, "Sec") == 0) {
        gene = "trnU";
    } else if (StringICmp (product, "Val") == 0) {
        gene = "trnV";
    } else if (StringICmp (product, "Trp") == 0) {
        gene = "trnW";
    } else if (StringICmp (product, "OTHER") == 0) {
        gene = "trnX";
    } else if (StringICmp (product, "Tyr") == 0) {
        gene = "trnY";
    } else if (StringICmp (product, "Glx") == 0) {
        gene = "trnZ";
    }
    return gene;
}


static CommentFeatPtr ParseGeneFromNoteForDefLine (CharPtr PNTR comment)
{
  CommentFeatPtr tdp;
  CharPtr    product_start, product_end, gene_start, gene_end, cp;
  Int4       product_len, gene_len;
  
  if (comment == NULL || *comment == NULL)
  {
    return NULL;
  }
  /* spacers are not genes */
  if (StringNICmp (*comment, "intergenic", 10) == 0 || StringNICmp (*comment, "spacer", 6) == 0) 
  {
    return NULL;
  }
  
  /* tRNA name must start with "tRNA-" and be followed by one uppercase letter and
   * two lowercase letters.
   */
  product_start = *comment;
  gene_start = product_start;
  while (*gene_start != 0 && !isspace (*gene_start) && *gene_start != ',') {
    gene_start++;
  }
  if (gene_start == product_start) {
    return NULL;
  }
  product_end = gene_start;

  /* if tRNA, don't require gene name, but parse if present */
  while (isspace (*gene_start)) {
    gene_start++;
  }
  if (*gene_start == '(') {
    /* parse in gene name */
    gene_start++;
    gene_end = gene_start;
    while (*gene_end != 0 && *gene_end != ')') {
      gene_end++;
    }
    if (*gene_end == 0) {
      return NULL;
    }
    cp = gene_end + 1;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  } else if (StringNICmp (gene_start, "intergenic", 10) == 0 || StringNICmp (gene_start, "spacer", 6) == 0) {
    /* spacers are not genes */
    return NULL;
  } else if (StringNCmp (product_start, "tRNA", 4) != 0) {
    /* only tRNAs can leave out gene names */
    return NULL;
  } else {
    cp = gene_start;
    gene_start = NULL;

  }

  /* skip over gene or genes if present */
  if (StringNCmp (cp, "genes", 5) == 0) {
    cp +=5;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  } else if (StringNCmp (cp, "gene", 4) == 0) {
    cp += 4;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  }

  tdp = (CommentFeatPtr) MemNew (sizeof (CommentFeatData));
  if (tdp == NULL)
  {
    return NULL;
  }
  product_len = product_end - product_start;
  tdp->product_name = (CharPtr) MemNew (sizeof (Char) * (1 + product_len));
  StringNCpy (tdp->product_name, product_start, product_len);
  tdp->product_name[product_len] = 0;
  
  if (gene_start != NULL) {
    gene_len = gene_end - gene_start;
    tdp->gene_name = (CharPtr) MemNew (sizeof (Char) * (1 + gene_len));
    StringNCpy (tdp->gene_name, gene_start, gene_len);
    tdp->gene_name[gene_len] = 0;
  }
  
  *comment = cp;
  return tdp;
}


static IntergenicSpacerDefPtr ParseIntergenicSpacerFromNoteForDef (CharPtr PNTR comment)
{
  IntergenicSpacerDefPtr idp;
  CharPtr                first_gene_start, dash, second_gene_start, second_gene_end, cp;
  Int4                   first_gene_len, second_gene_len;
  
  if (comment == NULL || *comment == NULL)
  {
    return NULL;
  }
  
  /* description must start with "trn" and be followed by one uppercase letter, followed
   * by a dash, followed by "trn", followed by one uppercase letter, followed by whitespace,
   * followed by the phrase "intergenic spacer".
   */
  first_gene_start = *comment;
  dash = first_gene_start;
  while (*dash != 0 && isalpha (*dash)) {
    dash++;
  }
  if (*dash != '-') {
    return NULL;
  }
  second_gene_start = dash + 1;
  second_gene_end = second_gene_start;
  while (*second_gene_end != 0 && isalpha (*second_gene_end)) {
    second_gene_end ++;
  }
  if (!isspace (*second_gene_end)) {
    return NULL;
  }
  cp = second_gene_end;
  while (isspace (*cp)) {
    cp++;
  }
  if (StringNCmp (cp, "intergenic spacer", 17) != 0) {
    return NULL;
  }
  
  idp = (IntergenicSpacerDefPtr) MemNew (sizeof (IntergenicSpacerDefData));
  if (idp == NULL)
  {
    return NULL;
  }
  
  first_gene_len = dash - first_gene_start;
  idp->first_gene = (CharPtr) MemNew (sizeof (Char) * (1 + first_gene_len));
  StringNCpy (idp->first_gene, first_gene_start, first_gene_len);
  idp->first_gene [first_gene_len] = 0;

  second_gene_len = second_gene_end - second_gene_start;
  idp->second_gene = (CharPtr) MemNew (sizeof (Char) * (1 + second_gene_len));
  StringNCpy (idp->second_gene, second_gene_start, second_gene_len);
  idp->second_gene [second_gene_len] = 0;

  *comment = cp + 17;
  return idp;
}

/* This creates a feature clause from a tRNADef structure. */
static FeatureClausePtr 
FeatureClauseFromParsedComment 
(CommentFeatPtr tdp,
 SeqFeatPtr misc_feat,
 Boolean    is_partial,
 BioseqPtr  bsp,
 DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  CharPtr gene_fmt = "%s (%s)";
  
  if (tdp == NULL)
  {
    return NULL;
  }
  
  fcp = NewFeatureClause ( misc_feat, bsp, rp);
  if (fcp != NULL)
  {
    fcp->feature_label_data.is_typeword_first = FALSE;
    fcp->feature_label_data.typeword = StringSave ("gene");
    if (tdp->gene_name == NULL) {
      fcp->feature_label_data.description = StringSave (tdp->product_name);
    } else {
      fcp->feature_label_data.description = (CharPtr) MemNew (sizeof (Char) * (StringLen (gene_fmt)
                                                                               + StringLen (tdp->gene_name)
                                                                               + StringLen (tdp->product_name)));
      if (fcp->feature_label_data.description != NULL)
      {
        sprintf (fcp->feature_label_data.description, gene_fmt,
                tdp->product_name, tdp->gene_name);
      }
    }
    if (is_partial)
    {
      fcp->interval = StringSave ("partial sequence");
    }
    else
    {
      fcp->interval = StringSave ("complete sequence");
    }
  }
  return fcp;
}


static CharPtr AdvancePastSeparators (CharPtr cp)
{
  if (cp == NULL || *cp == '0') return cp;

  if (*cp == ',')
  {
    cp++;
  }
  while (isspace (*cp))
  {
    cp++;
  }
  if (StringNCmp (cp, "and", 3) == 0)
  {
    cp += 3;
  }
  while (isspace (*cp))
  {
    cp++;
  }
  return cp;
}

#define MISCFEAT_TRNA_GENE 1
#define MISCFEAT_TRNA_SPACER 2

static ValNodePtr FreetRNAGenesAndSpacersList (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    if (list->choice == MISCFEAT_TRNA_GENE) {
      list->data.ptrvalue = CommentFeatFree (list->data.ptrvalue);
    } else if (list->choice == MISCFEAT_TRNA_SPACER) {
      list->data.ptrvalue = IntergenicSpacerDefFree (list->data.ptrvalue);
    }
    list->next = NULL;
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static ValNodePtr GettRNAGenesAndSpacers (CharPtr str)
{
  ValNodePtr list = NULL;
  CharPtr    cp;
  CommentFeatPtr gene, last_gene = NULL;
  IntergenicSpacerDefPtr spacer, last_spacer = NULL;
  Boolean                none_left = FALSE, names_correct = TRUE, alternating = TRUE;
  Int4                   last_item_type = 0;

  if (StringNICmp (str, "contains ", 9) == 0) {
    cp = str + 9;
  } else if (StringNICmp (str, "may contain ", 12) == 0) {
    cp = str + 12;
  } else {
    return NULL;
  }

  while (isspace (*cp))
  {
    cp ++;
  }
  
  while (!none_left && *cp != 0 && *cp != ';' && alternating && names_correct) {
    none_left = TRUE;
    gene = ParseGeneFromNoteForDefLine (&cp);
    if (gene != NULL) {
      /* if previous item was spacer, spacer names and gene names must agree */
      if (last_item_type == MISCFEAT_TRNA_SPACER && last_spacer != NULL) {
        if (gene->gene_name == NULL) {
          if (StringCmp (last_spacer->second_gene, s_tRNAGeneFromProduct (gene->product_name)) != 0) {
            names_correct = FALSE;
          }
        } else if (StringCmp (last_spacer->second_gene, gene->gene_name) != 0) {
          names_correct = FALSE;
        }
      }
      ValNodeAddPointer (&list, MISCFEAT_TRNA_GENE, gene);
      cp = AdvancePastSeparators (cp);
      none_left = FALSE;
      last_item_type = MISCFEAT_TRNA_GENE;
      last_gene = gene;
    }
 
    spacer = ParseIntergenicSpacerFromNoteForDef (&cp);
    if (spacer != NULL) {
      /* must alternate between genes and spacers */
      if (last_item_type == MISCFEAT_TRNA_SPACER) {
        alternating = FALSE;
      }
      /* spacer names and gene names must agree */
      if (last_gene != NULL) {
        if (last_gene->gene_name == NULL) {
          if (StringCmp (s_tRNAGeneFromProduct (last_gene->product_name), spacer->first_gene) != 0) {
            names_correct = FALSE;
          }
        } else if (StringCmp (last_gene->gene_name, spacer->first_gene) != 0) {
          names_correct = FALSE;
        }
      }
      ValNodeAddPointer (&list, MISCFEAT_TRNA_SPACER, spacer);
      cp = AdvancePastSeparators (cp);
      none_left = FALSE;
      last_item_type = MISCFEAT_TRNA_SPACER;
      last_spacer = spacer;
    }      
  }
  if ((*cp != 0 && *cp != ';') || !alternating || !names_correct) {
    list = FreetRNAGenesAndSpacersList (list);
  }
  return list;
}


/* This function produces a feature clause list that should replace the original
 * single clause for a misc_feat that contains a note with one or more tRNAs and
 * an intergenic spacer.
 */
static ValNodePtr 
ParsetRNAIntergenicSpacerElements 
(SeqFeatPtr misc_feat,
BioseqPtr   bsp,
DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  ValNodePtr head = NULL, vnp, list;
  Boolean partial5, partial3;
  IntergenicSpacerDefPtr spacer = NULL;
  Boolean                current_is_partial;

  if (misc_feat == NULL 
      || StringHasNoText (misc_feat->comment)) 
  {
    return NULL;
  }

  list = GettRNAGenesAndSpacers (misc_feat->comment);
  if (list != NULL) {
    if (StringNICmp (misc_feat->comment, "may contain ", 12) == 0) {
      fcp = NewFeatureClause (misc_feat, bsp, rp);
      fcp->feature_label_data.description = StringSave (misc_feat->comment + 12);
      fcp->interval = StringSave ("region");
      ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
    } else {
      CheckSeqLocForPartial (misc_feat->location, &partial5, &partial3);
      for (vnp = list; vnp != NULL; vnp = vnp->next) {
        current_is_partial = (partial5 && vnp == list) || (partial3 && vnp->next == NULL);
        if (vnp->data.ptrvalue == NULL) continue;
        if (vnp->choice == MISCFEAT_TRNA_GENE) {
          fcp = FeatureClauseFromParsedComment (vnp->data.ptrvalue, misc_feat, current_is_partial, bsp, rp);
          if (fcp != NULL) {
            ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
          }
        } else if (vnp->choice == MISCFEAT_TRNA_SPACER) {
          spacer = (IntergenicSpacerDefPtr) vnp->data.ptrvalue;
          fcp = NewFeatureClause ( misc_feat, bsp, rp);
          if (fcp != NULL)
          {
            fcp->feature_label_data.is_typeword_first = FALSE;
            fcp->feature_label_data.typeword = StringSave ("intergenic spacer");
            fcp->feature_label_data.description = (CharPtr) MemNew (10 * sizeof (Char));
            if (fcp->feature_label_data.description != NULL)
            {
              sprintf (fcp->feature_label_data.description, "%s-%s",
                      spacer->first_gene, spacer->second_gene);
            }
            if (current_is_partial)
            {
              fcp->interval = StringSave ("partial sequence");
            }
            else
            {
              fcp->interval = StringSave ("complete sequence");
            }
            ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
          }
        }
      }
    }
    list = FreetRNAGenesAndSpacersList (list);
  }
  return head;
}

static Boolean LIBCALLBACK IsSatelliteSequence (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_repeat_region)
  {
    return FALSE;
  }
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (StringCmp (gbqual->qual, "satellite") == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK IsPromoter (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_promoter) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsEndogenousVirusSourceFeature (
  SeqFeatPtr sfp
)
{
  BioSourcePtr biop;
  SubSourcePtr  ssp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return FALSE;
  if ((biop = sfp->data.value.ptrvalue) == NULL) return FALSE;
  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_endogenous_virus_name)
  {
    ssp = ssp->next;
  }
  if (ssp != NULL) return TRUE;
  return FALSE;
}

static CharPtr GetEndogenousVirusSourceFeatureDescription (
  SeqFeatPtr sfp
)
{
  BioSourcePtr biop;
  SubSourcePtr  ssp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return NULL;
  if ((biop = sfp->data.value.ptrvalue) == NULL) return NULL;
  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_endogenous_virus_name)
  {
    ssp = ssp->next;
  }
  if (ssp != NULL && ssp->name != NULL)
  {
    return StringSave (ssp->name);
  }
  return NULL;
}


static CharPtr noncoding_feature_keywords[] = {
  "similar to ",
  "contains "
};

static CharPtr find_noncoding_feature_keyword (
  CharPtr comment
)
{
  Int4 i, num_noncoding_feature_keywords, keywordlen;
  CharPtr cp, buffer;

  if (comment == NULL) return NULL;
  num_noncoding_feature_keywords = sizeof (noncoding_feature_keywords) / sizeof (CharPtr);
  for (i=0; i < num_noncoding_feature_keywords; i++)
  {
    keywordlen = StringLen (noncoding_feature_keywords [i]);
    buffer = comment;
    while ((cp = StringStr (buffer, noncoding_feature_keywords [i])) != NULL)
    {
      if ( StringNCmp (cp + keywordlen,
                       "GenBank Accession Number",
                       StringLen ("GenBank Accession Number")) != 0)
      {
        return cp + keywordlen;
      }
      else
      {
        buffer = cp + 1;
      }
    }
  }
  return NULL;
}

static Boolean LIBCALLBACK IsNoncodingProductFeat (
  SeqFeatPtr sfp
)
{
  if ( sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringStr (sfp->comment, "intergenic") != NULL
    || IsParsableList (sfp)
    || (find_noncoding_feature_keyword (sfp->comment) == NULL
      && (StringStr (sfp->comment, "nonfunctional ") == NULL
        || StringStr (sfp->comment, " due to ") == NULL)))
  {
    return FALSE;
  }

  return TRUE;
}

static CharPtr GetNoncodingProductFeatProduct (
  SeqFeatPtr sfp
)
{
  CharPtr productname;
  Int4    namelen, compare_len;
  CharPtr name_start, sep;

  if (sfp == NULL || sfp->comment == NULL) return NULL;

  if ((name_start = StringStr (sfp->comment, "nonfunctional ")) != NULL
    && (sep = StringStr (sfp->comment, " due to ")) != NULL
    && sep > name_start)
  {
    productname = StringSave (name_start);
    productname [ sep - name_start] = 0;
    return productname;
  }

  name_start = find_noncoding_feature_keyword (sfp->comment);
  if (name_start == NULL) return NULL;

  sep = StringStr (name_start, ";");
  if (sep == NULL)
  {
    namelen = StringLen (name_start);
  }
  else
  {
    namelen = sep - name_start;
  }

  productname = MemNew (namelen + 6);
  if (productname == NULL) return NULL;

  StringNCpy (productname, name_start, namelen);
  productname [namelen] = 0;
  
  /* remove sequence from end of name if present */
  compare_len = StringLen (" sequence");
  if (StringCmp (productname + namelen - compare_len, " sequence") == 0)
  {
    productname [ namelen - compare_len] = 0;
    namelen = StringLen (productname);
  }
  /* add "-like" if not present */
  compare_len = StringLen ("-like");
  if (StringCmp (productname + namelen - compare_len, "-like") != 0)
  {
    StringCat (productname, "-like");
    namelen = StringLen (productname);
  }
  return productname;
}

static Boolean LIBCALLBACK IsMiscFeat (
  SeqFeatPtr sfp
)
{
  if ( sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL)
  {
    return FALSE;
  }

  return TRUE;
}

static Boolean LIBCALLBACK IsOperon (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL 
    || sfp->idx.subtype != FEATDEF_operon)
  {
    return FALSE;
  }

  return TRUE;
}

static Boolean IsRecognizedFeature (
  SeqFeatPtr sfp
)
{
  if (IsGene (sfp)
    || IsCDS (sfp)
    || IsExon (sfp)
    || IsIntron (sfp)
    || IsLTR (sfp)
    || IsrRNA (sfp)
    || IstRNA (sfp)
    || IsmRNA (sfp)
    || IsMiscRNA (sfp)
    || IsncRNA (sfp)
    || IsPrecursorRNA (sfp)
    || Is3UTR (sfp)
    || Is5UTR (sfp)
    || IsMobileElement (sfp)
    || IsEndogenousVirusSequence (sfp)
    || IsEndogenousVirusSourceFeature (sfp)
    || IsDloop (sfp)
    || IsSatelliteSequence (sfp)
    || IsControlRegion (sfp)
    || IsIntergenicSpacer (sfp)
    || IsGeneCluster (sfp)
    || IsNoncodingProductFeat (sfp)
    || IsPromoter (sfp)
    || IsMiscFeat (sfp)
    || IsOperon (sfp))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* The following section of code contains functions for dealing with lists of
 * clauses.
 */

/* The functions for freeing the memory associated with lists of clauses 
 * are recursive.
 */
static void FreeListElement (ValNodePtr element);

/* This function simply frees the ValNodePtr, since there is no extra
 * memory associated with a DEFLINE_FEATLIST item - the sfp that is
 * pointed to by data.ptrvalue came from the sequence indexing functions
 * and should under no circumstances be freed.
 */
static void FreeFeatlist (
  ValNodePtr featlist
)
{

  if (featlist == NULL) return;
  ValNodeFree (featlist);
}
 
/* This function frees the memory associated with a FeatureClause, including
 * the memory associated with any subclauses.
 */
static void FreeClausePlusData (
  FeatureClausePtr fcp
)
{
  if (fcp->interval != NULL)
  {
    MemFree (fcp->interval);
    fcp->interval = NULL;
  }
  if (fcp->allelename != NULL)
  {
    MemFree (fcp->allelename);
    fcp->allelename = NULL;
  }
  if (fcp->feature_label_data.typeword != NULL)
  {
    MemFree (fcp->feature_label_data.typeword);
    fcp->feature_label_data.typeword = NULL;
  }
  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
    fcp->feature_label_data.description = NULL;
  }
  if (fcp->feature_label_data.productname != NULL)
  {
    MemFree (fcp->feature_label_data.productname);
    fcp->feature_label_data.productname = NULL;
  }
  if (fcp->featlist != NULL)
  {
    FreeListElement (fcp->featlist);
    fcp->featlist = NULL;
  }
  if (fcp->slp != NULL)
  {
    SeqLocFree (fcp->slp);
  }
}
    
/* This function frees the data associated with the FeatureClause
 * and then frees the ValNode.
 */
static void FreeClausePlus (
  ValNodePtr clauseplus
)
{
  FeatureClausePtr data_struct;

  if (clauseplus == NULL) return;
  data_struct = (FeatureClausePtr) clauseplus->data.ptrvalue;
  if (data_struct != NULL)
  {
    FreeClausePlusData (data_struct); 
    MemFree (data_struct);
    clauseplus->data.ptrvalue = NULL;
  }
  ValNodeFree (clauseplus);
}
 
/* This function frees a list of DEFLINE_FEATLIST, DEFLINE_REMOVEFEAT,
 * and DEFLINE_CLAUSEPLUS items, starting with the last item in the list.
 * It recursively frees memory associated with subclauses.
 */ 
static void FreeListElement (
  ValNodePtr element
)
{
  if (element == NULL) return;

  FreeListElement (element->next);
  element->next = NULL;
  if (element->choice == DEFLINE_FEATLIST 
    || element->choice == DEFLINE_REMOVEFEAT)
  {
    FreeFeatlist (element);
  }
  else if (element->choice == DEFLINE_CLAUSEPLUS)
  {
    FreeClausePlus (element);
  }
}

/* This function excises from the list pointed to by head all of the clauses
 * with the delete_me flag set to TRUE and all of the ValNodes with a choice
 * of DEFLINE_REMOVEFEAT.
 */
static void DeleteFeatureClauses (
  ValNodePtr PNTR head
)
{
  ValNodePtr vnp, prev;
  FeatureClausePtr fcp;
  Boolean          delete_this_one;

  if (head == NULL) return;

  prev = NULL;
  vnp = *head;
  while (vnp != NULL)
  {
    delete_this_one = FALSE;

    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp == NULL || fcp->delete_me || fcp->featlist == NULL) 
      {
        delete_this_one = TRUE;
      }
      else
      {
        DeleteFeatureClauses (&fcp->featlist);
        if (fcp->featlist == NULL) delete_this_one = TRUE;
      }
    }
    else if (vnp->choice == DEFLINE_REMOVEFEAT)
    {
      delete_this_one = TRUE;
    }

    if (delete_this_one)
    {
      if (prev == NULL)
      {
        *head = vnp->next;
        vnp->next = NULL;
        FreeListElement (vnp);
        if (*head == NULL) return;
        vnp = *head;
      }
      else
      {
        prev->next = vnp->next;
        vnp->next = NULL;
        FreeListElement (vnp);
        vnp = prev->next;
      }
    }
    else
    {
      prev = vnp;
      vnp = vnp->next;
    }
  }
}

/* This function counts the number of features in the feature list that
 * satisfy the itemmatch function (or all of them, if itemmatch is NULL).
 * If recurse_past_found_item, the function will not count features in
 * subclauses of features that satisfy itemmatch.
 */
static Int4 CountFeatures (
  ValNodePtr clause_list,
  matchFunction  itemmatch,
  Boolean    recurse_past_found_item
)
{
  ValNodePtr       vnp;
  Int4             num_features;
  FeatureClausePtr fcp;

  num_features = 0;
  for (vnp = clause_list;
       vnp != NULL;
       vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && (itemmatch == NULL || itemmatch (vnp->data.ptrvalue)))
    {
      num_features++;
      if (! recurse_past_found_item)
      {
        return num_features;
      }
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      num_features += CountFeatures (fcp->featlist,
                                     itemmatch,
                                     recurse_past_found_item);   
    }
  }
  return num_features;
}

/* The following section of code contains functions for grouping features. */

typedef struct matchruledata {
  matchFunction is_item;
  Int4          num_match_rules;
  matchFunction *match_rules;
} MatchRuleData, PNTR MatchRulePtr;

static void InitRuleForTopLevelClauses (MatchRulePtr mrp)
{
  if (mrp == NULL)
  {
    return;
  }
  mrp->num_match_rules = 4;
  mrp->match_rules = MemNew (mrp->num_match_rules
                                    * sizeof (matchFunction));
  if (mrp->match_rules == NULL) return;
  mrp->match_rules[0] = IsMobileElement;
  mrp->match_rules[1] = IsEndogenousVirusSourceFeature;
  mrp->match_rules[2] = IsOperon;
  mrp->match_rules[3] = IsGeneCluster;
}

static void InitRuleForBottomLevelClauses (MatchRulePtr mrp)
{
  if (mrp == NULL)
  {
    return;
  }
  mrp->num_match_rules = 6;
  mrp->match_rules = MemNew (mrp->num_match_rules
                                    * sizeof (matchFunction));
  if (mrp->match_rules == NULL) return;
  mrp->match_rules[0] = IsCDS;
  mrp->match_rules[1] = IsmRNA;
  mrp->match_rules[2] = IsGene;
  mrp->match_rules[3] = IsEndogenousVirusSourceFeature;
  mrp->match_rules[4] = IsOperon;
  mrp->match_rules[5] = IsGeneCluster;
}

/* NumGroupingRules is the number of features for which there is a list of
 * features to group under.
 * When grouping features, each feature in the list is examined sequentially.
 * If there is a rule set that applies to that feature, the entire feature 
 * list is searched for each feature type that this feature might group
 * beneath.  This preserves the biological order that was generated by the
 * original listing of features by sequence indexing.
 */
#define  NumGroupingRules 13
static MatchRulePtr InitializeGroupingRules()
{
  MatchRulePtr grouping_rules;

  grouping_rules = MemNew (NumGroupingRules * sizeof (MatchRuleData));
  if (grouping_rules == NULL) return NULL;

  grouping_rules[0].is_item = IsExon;
  grouping_rules[0].num_match_rules = 8;
  grouping_rules[0].match_rules = MemNew (grouping_rules[0].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[0].match_rules == NULL) return NULL;
  grouping_rules[0].match_rules[0] = IsCDS;
  grouping_rules[0].match_rules[1] = IsNoncodingProductFeat;
  grouping_rules[0].match_rules[2] = IsDloop;
  grouping_rules[0].match_rules[3] = IsmRNA;
  grouping_rules[0].match_rules[4] = IsGene;
  grouping_rules[0].match_rules[5] = IsEndogenousVirusSourceFeature;
  grouping_rules[0].match_rules[6] = IsOperon;
  grouping_rules[0].match_rules[7] = IsGeneCluster;

  grouping_rules[1].is_item = IsIntron;
  grouping_rules[1].num_match_rules = 8;
  grouping_rules[1].match_rules = MemNew (grouping_rules[1].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[1].match_rules == NULL) return NULL;
  grouping_rules[1].match_rules[0] = IsCDS;
  grouping_rules[1].match_rules[1] = IsNoncodingProductFeat;
  grouping_rules[1].match_rules[2] = IstRNA;
  grouping_rules[1].match_rules[3] = IsDloop;
  grouping_rules[1].match_rules[4] = IsGene;
  grouping_rules[1].match_rules[5] = IsEndogenousVirusSourceFeature;
  grouping_rules[1].match_rules[6] = IsOperon;
  grouping_rules[1].match_rules[7] = IsGeneCluster;
  
  grouping_rules[2].is_item = IsPromoter;
  InitRuleForBottomLevelClauses (grouping_rules + 2);
 
  grouping_rules[3].is_item = IsCDS;
  grouping_rules[3].num_match_rules = 5;
  grouping_rules[3].match_rules = MemNew (grouping_rules[3].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[3].match_rules == NULL) return NULL;
  grouping_rules[3].match_rules[0] = IsmRNA;
  grouping_rules[3].match_rules[1] = IsMobileElement;
  grouping_rules[3].match_rules[2] = IsEndogenousVirusSourceFeature;
  grouping_rules[3].match_rules[3] = IsOperon;
  grouping_rules[3].match_rules[4] = IsGeneCluster;
  
  grouping_rules[4].is_item = IsMobileElement;
  InitRuleForTopLevelClauses (grouping_rules + 4);
 
  grouping_rules[5].is_item = Is3UTR;
  InitRuleForBottomLevelClauses (grouping_rules + 5);
 
  grouping_rules[6].is_item = Is5UTR;
  InitRuleForBottomLevelClauses (grouping_rules + 6);
 
  grouping_rules[7].is_item = IsLTR;
  InitRuleForBottomLevelClauses (grouping_rules + 7);
 
  grouping_rules[8].is_item = IsGene;
  InitRuleForTopLevelClauses (grouping_rules + 8);
  
  grouping_rules[9].is_item = IsIntergenicSpacer;
  InitRuleForTopLevelClauses (grouping_rules + 9);

  grouping_rules[10].is_item = IsNoncodingProductFeat;
  InitRuleForTopLevelClauses (grouping_rules + 10);
 
  grouping_rules[11].is_item = IsOperon;
  InitRuleForTopLevelClauses (grouping_rules + 11);
 
  grouping_rules[12].is_item = IsGeneCluster;
  InitRuleForTopLevelClauses (grouping_rules + 12);

  return grouping_rules; 
}

static void FreeGroupingRules(
  MatchRulePtr grouping_rules
)
{
  Int4 i;

  if (grouping_rules == NULL) return;

  for (i = 0; i < NumGroupingRules; i++)
  {
    if (grouping_rules[i].match_rules != NULL)
    MemFree (grouping_rules[i].match_rules);
    grouping_rules[i].match_rules = NULL;
  }
  MemFree (grouping_rules);
}

static Boolean IsmRNASequence (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  
  if (bsp == NULL || bsp->mol != Seq_mol_rna || bsp->descr == NULL)
  {
    return FALSE;
  }
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo)
  {
    sdp = sdp->next;
  }
  if (sdp == NULL || sdp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  
  if (mip->biomol == 3)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

typedef struct matchcandidate 
{
  ValNodePtr matched_clause;
  SeqLocPtr  slp;
} MatchCandidateData, PNTR MatchCandidatePtr;

/* This function searches the search_list for features that satisfy the
 * match function and satisfy locational requirements relative to the
 * clause.
 * If more than one clause meets the match requirements, the smallest one
 * is chosen.
 */
static void FindBestMatchCandidate 
(FeatureClausePtr  clause,
 ValNodePtr        search_list,
 FeatureClausePtr  search_parent,
 matchFunction     match,
 Boolean           gene_cluster_opp_strand,
 BioseqPtr         bsp,
 MatchCandidatePtr current_candidate)
{
  ValNodePtr       search_clause;
  SeqFeatPtr       addsfp, clause_sfp;
  FeatureClausePtr searchfcp;
  SeqLocPtr        slp;

  if (clause == NULL || clause->slp == NULL || current_candidate == NULL) return;
  
  clause_sfp = (SeqFeatPtr) (clause->featlist->data.ptrvalue);

  for (search_clause = search_list;
       search_clause != NULL;
       search_clause = search_clause->next)
  {
    if (search_clause->data.ptrvalue == clause)
      continue;
    if (search_clause->choice == DEFLINE_FEATLIST
      && search_clause->data.ptrvalue != NULL)
    {
      addsfp = search_clause->data.ptrvalue;
      /* slp is the location of the feature we are trying to
       * group this feature with
       */ 
      if (search_parent != NULL)
      {
        slp = search_parent->slp;
      }
      else
      {
        slp = addsfp->location;
      }
      if (match (search_clause->data.ptrvalue))
      { 
        /* Transposons, insertion sequences, integrons, and endogenous virii
         * take subfeatures regardless of whether the subfeature is
         * on the same strand.
         * Gene Clusters can optionally take subfeatures on either
         * strand (gene_cluster_opp_strand is flag).
         * Promoters will match up to features that are adjacent.
         * Introns will match up to coding regions whose intervals
         * are adjacent to the endpoints of the intron, or to other
         * features if the intron location is inside the other feature.
         * All other feature matches must be that the feature to
         * go into the clause must fit inside the location of the
         * other clause.
         */
        if (((match == IsMobileElement
             || match == IsEndogenousVirusSourceFeature
             || (match == IsGeneCluster && gene_cluster_opp_strand))
            && SeqLocAinB (clause->slp, slp) > -1)
          || IsLocAInBonSameStrand (clause->slp, slp)
          || ( IsPromoter (clause_sfp)
            && IsAAdjacentToB (clause->slp, search_parent->slp, bsp,
                               ADJACENT_TYPE_UPSTREAM, TRUE))
          || (IsmRNASequence (bsp) 
              && match != IsMobileElement 
              && match != IsEndogenousVirusSourceFeature
              && match != IsGeneCluster)
          || (match == IsCDS 
              && IsIntron (clause_sfp)
              && IsAEmptyIntervalOfB (clause->slp, slp, bsp))
          || (match == IsCDS
              && IsExon (clause_sfp)
              && LocAContainsIntervalOfB (clause->slp, slp)))
        {
          /* if we don't already have a candidate, or if this
           * candidate's location is inside the current candidate,
           * take this candidate.
           */
          if (current_candidate->matched_clause == NULL
              || SeqLocAinB (slp, current_candidate->slp) > 0)
          {
            current_candidate->matched_clause = search_clause;
            current_candidate->slp = slp;
          }
        }
      }
    }
    else if (search_clause->choice == DEFLINE_CLAUSEPLUS 
      && search_clause->data.ptrvalue != NULL)
    {
      searchfcp = search_clause->data.ptrvalue;
      FindBestMatchCandidate (clause, searchfcp->featlist, searchfcp,
                              match, gene_cluster_opp_strand, bsp,
                              current_candidate);
    }
  }
}


/* This function iterates through the matches in the specified grouping rule.
 * If more than one match is found, the clause with the smallest location is
 * used.
 * If a match is found, the clause is added to the list of clauses for that
 * feature's parent clause.
 */
static Boolean GroupClauseByRule (
  FeatureClausePtr clause,
  ValNodePtr       search_list,
  MatchRulePtr     grouping_rule,
  Boolean          gene_cluster_opp_strand,
  BioseqPtr        bsp
)
{
  Int4               rule_index;
  MatchCandidateData mcd;  
  Boolean            rval = FALSE;
  ValNodePtr         newfeat;
  
  mcd.slp = NULL;
  mcd.matched_clause = NULL;

  for (rule_index = 0;
       rule_index < grouping_rule->num_match_rules;
       rule_index ++)
  {
  
    FindBestMatchCandidate (clause, search_list, NULL, 
                            grouping_rule->match_rules[rule_index],
                            gene_cluster_opp_strand, bsp, &mcd);
  }
  if (mcd.matched_clause != NULL)
  {
    newfeat = ValNodeNew (mcd.matched_clause);
    if (newfeat == NULL) return FALSE;
    newfeat->choice = DEFLINE_CLAUSEPLUS;
    newfeat->data.ptrvalue = clause;
    rval = TRUE;
  }
  return rval;
}


/* This function determines whether a subclause contains just a 3' UTR feature
 * and no other details.
 */
static Boolean Is3UTRClause (FeatureClausePtr clause)
{  
  if (clause == NULL
      || clause->featlist == NULL 
      || clause->featlist->choice != DEFLINE_FEATLIST
      || clause->featlist->data.ptrvalue == NULL
      || clause->featlist->next != NULL)
  {
    return FALSE;
  }
  return Is3UTR(clause->featlist->data.ptrvalue);
}


/* This function will move 3' UTRs to the end of any subfeat lists
 * so that they can be listed after the partial/complete CDS.
 */
static void Move3UTRToEndOfSubFeatList (ValNodePtr clause_list)
{
  ValNodePtr vnp, prev, last_vnp;
  FeatureClausePtr clause;
  
  if (clause_list == NULL || clause_list->next == NULL) 
  {
    return;
  }
  prev = clause_list;
  for (vnp = clause_list->next; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL && vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      clause = vnp->data.ptrvalue;
      if (Is3UTRClause (clause))
      {
        if (vnp->next != NULL)
        {
          /* move to end of clause list */
          last_vnp = vnp->next;
          while (last_vnp->next != NULL)
          {
            last_vnp = last_vnp->next;
          }
          prev->next = vnp->next;
          last_vnp->next = vnp;
          vnp->next = NULL;
        }
      }
      else
      {
        prev = vnp;
        Move3UTRToEndOfSubFeatList (clause->featlist);
      }
    }
    else
    {
      prev = vnp;
    }
  }  
}

/* This function iterates over the list of features, attempting to find and
 * apply grouping rules for each feature.
 */
static void GroupAllClauses (
  ValNodePtr PNTR clause_list,
  Boolean         gene_cluster_opp_strand,
  BioseqPtr       bsp
)
{
  MatchRulePtr     grouping_rules;
  ValNodePtr       vnp, prev;
  FeatureClausePtr clause;
  SeqFeatPtr       main_feat;
  Int4             rule_index;

  grouping_rules = InitializeGroupingRules();
  if (grouping_rules == NULL) return;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      clause = vnp->data.ptrvalue;
      if (clause->featlist != NULL
        && clause->featlist->choice == DEFLINE_FEATLIST
        && clause->featlist->data.ptrvalue != NULL)
      {
        main_feat = clause->featlist->data.ptrvalue;
        for (rule_index = 0;
             rule_index < NumGroupingRules 
               && ! grouping_rules[rule_index].is_item (main_feat);
             rule_index++)
        {
        }
        if (rule_index < NumGroupingRules)
        {
          if ( GroupClauseByRule (clause, *clause_list,
                                  grouping_rules + rule_index, 
                                  gene_cluster_opp_strand,
                                  bsp))
          {
            vnp->data.ptrvalue = NULL;
          }
        }
      }
    }
  }
  FreeGroupingRules(grouping_rules);

  vnp = *clause_list;
  prev = NULL;
  while (vnp != NULL)
  {
    if (vnp->data.ptrvalue == NULL)
    {
      if (prev == NULL)
      {
        *clause_list = vnp->next;
        vnp->next = NULL;
        ValNodeFree (vnp);
        vnp = *clause_list;
      }
      else
      {
        prev->next = vnp->next;
        vnp->next = NULL;
        ValNodeFree (vnp);
        vnp = prev->next;
      }
    }
    else
    {
      prev = vnp;
      vnp = vnp->next;
    }
  }
  
  Move3UTRToEndOfSubFeatList (*clause_list);
}

/* This function exists to handle the special case where two or more exons
 * are alternatively spliced, but there are no CDSs to represent some of the
 * alternatively spliced forms.  In order to make sure that all of the exons
 * that are alternatively spliced together appear with the CDS, they are
 * temporarily consolidated into a single clause with a location that
 * is the intersection of the exons' locations.  The clause will be
 * re-expanded after grouping by the ExpandAltSplicedExons function.
 */
static void GroupAltSplicedExons (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  Boolean         delete_now
)
{
  ValNodePtr       clause, search_clause, vnp;
  FeatureClausePtr fcp, search_fcp;
  SeqFeatPtr       sfp, search_sfp;
  SeqLocPtr        new_slp;

  if (clause_list == NULL) return;
 
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if ( ! fcp->is_alt_spliced
      || fcp->delete_me
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }
    sfp = fcp->featlist->data.ptrvalue;
    if ( ! IsExon (sfp))
    {
      continue;
    }

    for ( search_clause = clause->next; 
          search_clause != NULL
            && search_clause->choice == DEFLINE_CLAUSEPLUS
            && search_clause->data.ptrvalue != NULL
            && (search_fcp = search_clause->data.ptrvalue) != NULL
            && ! search_fcp->delete_me
            && search_fcp->is_alt_spliced
            && search_fcp->featlist != NULL
            && search_fcp->featlist->choice == DEFLINE_FEATLIST
            && (search_sfp = search_fcp->featlist->data.ptrvalue) != NULL
            && IsExon (search_sfp)
            && TestFeatOverlap (sfp, search_sfp, SIMPLE_OVERLAP) != -1;
          search_clause = search_clause->next)
    {
      vnp = ValNodeNew (fcp->featlist);
      if (vnp == NULL) return;
      vnp->choice = DEFLINE_FEATLIST;
      vnp->data.ptrvalue = search_sfp;
      search_fcp->delete_me = TRUE;
      new_slp = SeqLocIntersection (fcp->slp, search_fcp->slp, bsp);
      SeqLocFree (fcp->slp);
      fcp->slp = new_slp;
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (clause_list);
  }
}

/* This function expands a clause filled with alternatively-spliced exons
 * that was created in the GroupAltSplicedExons function.
 */
static void ExpandAltSplicedExons (
  ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr clause, rest_of_list, featlist, new_clause;
  FeatureClausePtr fcp, new_fcp;
  SeqFeatPtr sfp;

  for (clause = clause_list;
       clause != NULL;
       clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (fcp = clause->data.ptrvalue) == NULL
      || fcp->featlist == NULL)
    {
      continue;
    }
    if ( fcp->featlist->choice == DEFLINE_FEATLIST
      && (sfp = fcp->featlist->data.ptrvalue) != NULL
      && IsExon (sfp)
      && fcp->featlist->next != NULL
      && fcp->featlist->next->choice == DEFLINE_FEATLIST
      && IsExon (fcp->featlist->next->data.ptrvalue))
    {
      rest_of_list = clause->next;
      clause->next = NULL;
      for (featlist = fcp->featlist->next;
           featlist != NULL
             && featlist->choice == DEFLINE_FEATLIST
             && IsExon (featlist->data.ptrvalue);
           featlist = featlist->next)
      {
        new_clause = ValNodeNew (clause);
        if (new_clause == NULL) return;
        new_fcp = NewFeatureClause (featlist->data.ptrvalue, bsp, rp);
        if (new_fcp == NULL) return;
        new_fcp->grp = fcp->grp;
        new_fcp->is_alt_spliced = fcp->is_alt_spliced;
        new_fcp->make_plural = fcp->make_plural;
        new_clause->choice = DEFLINE_CLAUSEPLUS;
        new_clause->data.ptrvalue = new_fcp;
      }
      ValNodeFree (fcp->featlist->next);
      fcp->featlist->next = NULL;
      new_clause->next = rest_of_list;

      /* put back location for first exon - was reduced to union of 
       * all exon intervals in GroupAltSplicedExons
       */
      SeqLocFree (fcp->slp);
      sfp = fcp->featlist->data.ptrvalue;
      fcp->slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
    }
    else
    {
      ExpandAltSplicedExons (fcp->featlist, bsp, rp);
    }
  }
}



static Boolean DoFeaturesShareGene (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  Boolean share_gene = FALSE;
  SeqFeatPtr found_gene1, found_gene2;
  
  if (sfp1 != NULL && sfp2 != NULL
      && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp1))
      && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp2)))
  {
    found_gene1 = SeqMgrGetOverlappingGene (sfp1->location, NULL);
    found_gene2 = SeqMgrGetOverlappingGene (sfp2->location, NULL);
    if (found_gene1 == found_gene2 && found_gene1 != NULL)
    {
      share_gene = TRUE;
    }
  }
  return share_gene;
}

/* This function determines whether two features share the same product name */
static Boolean 
DoProductNamesMatch 
(SeqFeatPtr sfp1,
 SeqFeatPtr sfp2, 
 BioseqPtr  bsp,
 DeflineFeatureRequestListPtr rp)
{
  CharPtr productname1;
  CharPtr productname2;
  Boolean names_match = FALSE;
   
  productname1 = GetProductName (sfp1, bsp, rp);
  productname2 = GetProductName (sfp2, bsp, rp);
  if (StringHasNoText (productname1) && StringHasNoText (productname2))
  {
    names_match = TRUE;
  }
  else if (StringCmp (productname1, productname2) == 0)
  {
    names_match = TRUE;
  }
  
  productname1 = MemFree (productname1);
  productname2 = MemFree (productname2);
  
  return names_match;  
}

/* This function should combine CDSs that do not have a joined location
 * but are part of the same gene and have the same protein name.
 */
static void GroupSegmentedCDSs (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  Boolean         delete_now,
  DeflineFeatureRequestListPtr rp
)
{
  ValNodePtr       clause, search_clause, vnp;
  FeatureClausePtr fcp, search_fcp;
  SeqFeatPtr       sfp, search_sfp;
  SeqLocPtr        new_slp;

  if (clause_list == NULL) return;
 
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if (fcp->delete_me
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }
    sfp = fcp->featlist->data.ptrvalue;
    if ( ! IsCDS (sfp))
    {
      continue;
    }

    for ( search_clause = clause->next; 
          search_clause != NULL;
          search_clause = search_clause->next)
    {
      if (search_clause->choice != DEFLINE_CLAUSEPLUS
          || search_clause->data.ptrvalue == NULL
          || (search_fcp = search_clause->data.ptrvalue) == NULL
          || search_fcp->delete_me
          || search_fcp->featlist == NULL
          || search_fcp->featlist->choice != DEFLINE_FEATLIST
          || (search_sfp = search_fcp->featlist->data.ptrvalue) == NULL
          || ! IsCDS (search_sfp)
          || ! DoFeaturesShareGene (sfp, search_sfp)
          || ! DoProductNamesMatch (sfp, search_sfp, bsp, rp))
      {
        continue;
      }
      vnp = ValNodeNew (fcp->featlist);
      if (vnp == NULL) return;
      vnp->choice = DEFLINE_FEATLIST;
      vnp->data.ptrvalue = search_sfp;
      search_fcp->delete_me = TRUE;
      new_slp = SeqLocMerge (bsp, fcp->slp, search_fcp->slp,
                             FALSE, TRUE, FALSE);

      SeqLocFree (fcp->slp);
      fcp->slp = new_slp;
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (clause_list);
  }
}


/* This function searches this list for clauses to which this gene should
 * apply.  This is not taken care of by the GroupAllClauses function
 * because genes are added to clauses as a GeneRefPtr instead of as an
 * additional feature in the list, and because a gene can apply to more
 * than one clause, while other features should really only belong to
 * one clause.
 */
static Boolean AddGeneToClauses 
( SeqFeatPtr gene,
  CharPtr    gene_productname,
  ValNodePtr clause_list,
  Boolean    suppress_locus_tag)
{
  ValNodePtr    clause;
  FeatureClausePtr fcp;
  SeqFeatPtr    sfp, found_gene;
  GeneRefPtr    grp;
  Boolean    used_gene;
  
  if (gene == NULL || gene->data.value.ptrvalue == NULL) return FALSE;
  if (clause_list == NULL) return FALSE;

  used_gene = FALSE;
  grp = gene->data.value.ptrvalue;
  for (clause = clause_list; clause != NULL; clause = clause->next)
  {
    fcp = clause->data.ptrvalue;
    if (fcp == NULL || fcp->featlist == NULL) return FALSE;
    sfp = fcp->featlist->data.ptrvalue;
    if (sfp != NULL && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp))
        && (IsCDS (sfp)
            || IsrRNA (sfp)
            || IstRNA (sfp)
            || IsmRNA (sfp)
            || IsMiscRNA (sfp)
            || IsncRNA (sfp)
            || IsPrecursorRNA (sfp)
            || IsNoncodingProductFeat (sfp)))
    {
      if (fcp->grp == NULL)
      {
        found_gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
        if (found_gene != NULL)
        {
          fcp->grp = (GeneRefPtr) found_gene->data.value.ptrvalue;
        }
      }

      if (fcp->grp != NULL && DoGenesMatch (fcp->grp, grp, suppress_locus_tag))
      {
        used_gene = TRUE;
        if (gene_productname != NULL
          && fcp->feature_label_data.productname == NULL
          && IsCDS (sfp))
        {
          fcp->feature_label_data.productname =
              StringSave (gene_productname);
        }
      }
      else if (fcp->grp == NULL
        && IsLocAInBonSameStrand (sfp->location, gene->location))
      {
        fcp->grp = grp;
        used_gene = TRUE;
        if (gene_productname != NULL
          && fcp->feature_label_data.productname == NULL
          && IsCDS (sfp))
        {
          fcp->feature_label_data.productname =
              StringSave (gene_productname);
        }
      }
    }
  }
  return used_gene;
}

/* This function iterates through the list of features and calls
 * AddGeneToClauses for each gene feature it finds.
 */
static void GroupGenes (ValNodePtr PNTR clause_list, Boolean suppress_locus_tag)
{
  ValNodePtr  vnp;
  ValNodePtr  featlist;
  FeatureClausePtr fcp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS) return;
    fcp = (FeatureClausePtr) vnp->data.ptrvalue;
    if (fcp == NULL) return;

    featlist = fcp->featlist;
    if (featlist != NULL
      && featlist->choice == DEFLINE_FEATLIST
      && IsGene (featlist->data.ptrvalue))
    {
      AddGeneToClauses (featlist->data.ptrvalue,
                        fcp->feature_label_data.productname,
                        vnp->next, suppress_locus_tag);
    }
  } 
}

/* This function searches this list for clauses to which this mRNA should
 * apply.  This is not taken care of by the GroupAllClauses function
 * because when an mRNA is added to a CDS, the product for the clause is
 * replaced and the location for the clause is expanded, rather than simply
 * adding the mRNA as an additional feature in the list, and because an 
 * mRNA can apply to more than one clause, while other features should 
 * really only belong to one clause.
 */
static Boolean AddmRNAToClauses 
( SeqFeatPtr mRNA,
  ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr    clause;
  FeatureClausePtr fcp;
  SeqFeatPtr    sfp;
  Boolean    used_mRNA;
  CharPtr       productname;
  SeqLocPtr     new_slp;
  
  if (mRNA == NULL || mRNA->data.value.ptrvalue == NULL) return FALSE;
  if (clause_list == NULL) return FALSE;

  used_mRNA = FALSE;
  productname = GetProductName (mRNA, bsp, rp);
  if (productname == NULL) return TRUE;

  for (clause = clause_list; clause != NULL; clause = clause->next)
  {
    fcp = clause->data.ptrvalue;
    if (fcp == NULL || fcp->featlist == NULL) return FALSE;
    sfp = fcp->featlist->data.ptrvalue;
    if (sfp == NULL)
    {
    }
    else if (IsCDS (sfp)
      && fcp->feature_label_data.productname != NULL
      && StringCmp (fcp->feature_label_data.productname, productname) == 0)
    {
      used_mRNA = TRUE;
      fcp->has_mrna = TRUE;
      if (IsLocAInBonSameStrand (sfp->location, mRNA->location))
      {
        new_slp = SeqLocMerge (bsp, fcp->slp, mRNA->location,
                                 FALSE, TRUE, FALSE);
        if (new_slp == NULL) return FALSE;
        if (fcp->slp != NULL)
        {
          SeqLocFree (fcp->slp);
        }
        fcp->slp = new_slp;
      }
    }
    else if (fcp->feature_label_data.productname == NULL
      && (IsCDS (sfp) || IsGene (sfp))
      && (IsLocAInBonSameStrand (sfp->location, mRNA->location)
        || IsLocAInBonSameStrand (mRNA->location, sfp->location)))
    {
      fcp->feature_label_data.productname = StringSave (productname);
      used_mRNA = TRUE;
      fcp->has_mrna = TRUE;
      new_slp = SeqLocMerge (bsp, fcp->slp, mRNA->location,
                                 FALSE, TRUE, FALSE);
      if (new_slp == NULL) return FALSE;
      if (fcp->slp != NULL)
      {
        SeqLocFree (fcp->slp);
      }
      fcp->slp = new_slp;
    }
  }
  return used_mRNA;
}

/* This function iterates through the list of features and calls
 * AddmRNAToClauses for each mRNA feature it finds.
 */
static void GroupmRNAs (
  ValNodePtr PNTR clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp
)
{
  ValNodePtr  vnp;
  ValNodePtr  featlist;
  FeatureClausePtr fcp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS) return;
    fcp = (FeatureClausePtr) vnp->data.ptrvalue;
    if (fcp == NULL) return;

    featlist = fcp->featlist;
    if (featlist != NULL
      && featlist->choice == DEFLINE_FEATLIST
      && IsmRNA (featlist->data.ptrvalue))
    {
      if (AddmRNAToClauses (featlist->data.ptrvalue, *clause_list, bsp, rp))
      {
        fcp->delete_me = TRUE;
      }
    }
  } 
  DeleteFeatureClauses (clause_list);
}

/* This section of code contains functions for generating labels for
 * clauses for the definition lines.
 */

/* This function examines the specified typeword and determines whether it
 * should appear before or after the description of the feature in the
 * definition line.
 */
static Boolean IsTypeWordFirst (
  CharPtr typeword
)
{
  Int4 i;
  if (typeword == NULL) return FALSE;
  if (StringCmp (typeword, "exon") == 0
    || StringCmp (typeword, "intron") == 0
    || StringCmp (typeword, "endogenous virus") == 0)
  {
    return TRUE;
  }
  else
  {
    i = StartsWithMobileElementKeyword (typeword);
    if (i >= 0 && i != eMobileElementIntegron) {
      return TRUE;
    }
    return FALSE;
  }
}

/* This function determines the word to use to indicate what type of feature
 * is being described in the definition line.  For certain feature types,
 * the word to use in the definition line varies based on the type of
 * molecule in the record.
 */
static CharPtr GetFeatureTypeWord (
  Uint1 biomol,
  SeqFeatPtr sfp
)
{
  if (sfp == NULL) return NULL;
  if ( IsExon (sfp))
  {
    return StringSave ("exon");
  } 
  else if(IsIntron (sfp))
  {
    return StringSave ("intron");
  }
  else if (IsEndogenousVirusSequence (sfp))
  {
    return StringSave ("endogenous virus");
  }
  else if (IsControlRegion (sfp))
  {
    return StringSave ("control region");
  }
  else if (IsEndogenousVirusSourceFeature (sfp))
  {
    return StringSave ("endogenous virus");
  }
  else if (IsDloop (sfp))
  {
    return StringSave ("D-loop");
  }
  else if (IsLTR (sfp))
  {
    return StringSave ("LTR");
  }
  else if (Is3UTR (sfp))
  {
    return StringSave ("3' UTR");
  }
  else if (Is5UTR (sfp))
  {
    return StringSave ("5' UTR");
  }
  else if (IsOperon (sfp))
  {
    return StringSave ("operon");
  }
  else if (biomol == MOLECULE_TYPE_GENOMIC || biomol == MOLECULE_TYPE_CRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene");
    }
    else
    {
      return StringSave ("gene");
    }
  }
  else if ( IsrRNA (sfp) || IsncRNA (sfp))
  {
    return NULL;
  }
  else if (IsPrecursorRNA (sfp))
  {
    return StringSave ("precursor RNA");
  }
  else if (biomol == MOLECULE_TYPE_MRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene mRNA");
    }
    else
    {
      return StringSave ("mRNA");
    }
  }
  else if (biomol == MOLECULE_TYPE_PRE_MRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene precursor RNA");
    }
    else
    {
      return StringSave ("precursor RNA");
    }
  }
  else if (biomol == MOLECULE_TYPE_OTHER_GENETIC_MATERIAL)
  {
    return StringSave ("gene");
  }
  return StringSave ("");
}

/* Frequently the product associated with a feature is listed as part of the
 * description of the feature in the definition line.  This function determines
 * the name of the product associated with this specific feature.  Some
 * features will be listed with the product of a feature that is associated
 * with the feature being described - this function does not look at other
 * features to determine a product name.
 * If the feature is a misc_feat with particular keywords in the comment,
 * the product will be determined based on the contents of the comment.
 * If the feature is a CDS and is marked as pseudo, the product will be
 * determined based on the contents of the comment.
 * If the feature is a gene and has different strings in the description than
 * in the locus or locus tag, the description will be used as the product for
 * the gene.
 * If none of the above conditions apply, the sequence indexing context label
 * will be used to obtain the product name for the feature.
 */
static CharPtr GetProductName 
( SeqFeatPtr cds,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr protein_name;
  CharPtr semicolon;
  size_t len_to_copy;
  SeqMgrFeatContext  context;
  GeneRefPtr  grp;
  CharPtr gene_name;
  RnaRefPtr rrp;
  Boolean suppress_locus_tag = FALSE;

  if (cds == NULL) return NULL;
  protein_name = NULL;
  if (rp != NULL)
  {
    suppress_locus_tag = rp->suppress_locus_tags;
  }
  if (IsNoncodingProductFeat (cds))
  {
    return GetNoncodingProductFeatProduct (cds);
  }
  else if (cds->data.choice == SEQFEAT_CDREGION && cds->pseudo)
  {
    if (cds->comment != NULL)
    {
      semicolon = StringChr (cds->comment, ';');
      if (semicolon != NULL)
      {
        len_to_copy = semicolon - cds->comment;
      }
      else
      {
        len_to_copy = StringLen (cds->comment);
      }
      protein_name = MemNew (len_to_copy + 1);
      if (protein_name == NULL) return NULL;
      StringNCpy (protein_name, cds->comment, len_to_copy);
      protein_name[len_to_copy] = 0;
    }
    return protein_name;
  }
  else if (cds->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) cds->data.value.ptrvalue;
    if (grp == NULL) return NULL;
    gene_name = GetGeneName (grp, suppress_locus_tag);
    if (grp->desc != NULL
      && StringCmp (grp->desc, gene_name) != 0)
    {
      return StringSave (grp->desc);
    }
#if 0
    /* removed by request from Linda Yankie */    
    if (grp->locus_tag != NULL && ! suppress_locus_tag
      && StringCmp (grp->locus_tag, gene_name) != 0)
    {
      return StringSave (grp->locus_tag);
    }
#endif    
  }
  else if (IsncRNA (cds))
  {
    return GetncRNAProduct(cds, rp == NULL ? FALSE : rp->use_ncrna_note);
  }
  else if (IstRNA (cds) 
           && SeqMgrGetDesiredFeature (0, bsp, 0, 0, cds, &context) == cds
           && context.label != NULL)
  {
    if (StringCmp (context.label, "Xxx") == 0) {
      protein_name = StringSave ("tRNA-OTHER");
    } else {
      protein_name = MemNew ( StringLen (context.label) + 6);
      if ( protein_name == NULL) return NULL;
      sprintf (protein_name, "tRNA-%s", context.label);
    }
    return protein_name;
  }
  else if (cds->data.choice == SEQFEAT_RNA)
  {
    
    rrp = (RnaRefPtr) cds->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue))
    {
      return StringSave (rrp->ext.value.ptrvalue);
    }
  }
  else if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, cds, &context) == cds
           && context.label != NULL)
  {
    if ((IsCDS(cds) && StringCmp (context.label, "CDS") != 0)
        || (IsmRNA(cds) && StringCmp (context.label, "mRNA") != 0)
        || (! IsCDS(cds) && ! IsmRNA(cds)))
    {
      protein_name = StringSave (context.label);
      return protein_name;
    }
  }
  return NULL;
}

/* This function searches a list of features recursively for a
 * feature that satisfies the itemmatch condition and is associated with
 * the same gene as the fcp clause passed to the function.
 * This is used to obtain a product for a feature that may share a gene with
 * a product-producing feature but may not be contained in the interval of
 * the product-producing feature.
 */
static FeatureClausePtr FindProductInFeatureList (
  FeatureClausePtr fcp,
  ValNodePtr       clause_list,
  matchFunction    itemmatch,
  Boolean          suppress_locus_tag)
{
  ValNodePtr       vnp;
  FeatureClausePtr vnp_fcp;
  
  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      vnp_fcp = vnp->data.ptrvalue;
      if (DoGenesMatch (vnp_fcp->grp, fcp->grp, suppress_locus_tag)
        && vnp_fcp->featlist != NULL
        && vnp_fcp->featlist->choice == DEFLINE_FEATLIST
        && itemmatch (vnp_fcp->featlist->data.ptrvalue))
      {
        return vnp_fcp;
      }
      else
      {
        vnp_fcp = FindProductInFeatureList (fcp, vnp_fcp->featlist,
                                            itemmatch, suppress_locus_tag);
        if (vnp_fcp != NULL) return vnp_fcp;
      }
    }
  }
  return NULL;
}

/* This function uses the available information in the clause to generate
 * a description from the name of the gene (if any) and the name of the
 * product for the feature (if any).
 * If there is only a gene, the description will be the name of the gene.
 * If there is only a product, the description will be the name of the product.
 * If there is a gene and a product, the description will be the name of
 * the product followed by the name of the gene in parentheses.
 */
static CharPtr GetGeneProtDescription 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr    sfp;
  CharPtr    protein_name;
  CharPtr    gene_name;
  size_t    description_length;
  CharPtr    description;

  if (fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return NULL;
  }
  sfp = fcp->featlist->data.ptrvalue;

  description_length = 0;

  if (fcp->feature_label_data.productname != NULL)
  {
    protein_name = StringSave (fcp->feature_label_data.productname);
  }
  else
  {
    protein_name = GetProductName (sfp, bsp, rp);
    if (protein_name == NULL && IsGene (sfp))
    {
      
    }
  }
  if (protein_name != NULL)
  {
    description_length += StringLen (protein_name);
  }
     
  gene_name = GetGeneName (fcp->grp, rp == NULL ? FALSE : rp->suppress_locus_tags);
  if (gene_name != NULL)
  {
    description_length += StringLen (gene_name);
    if (protein_name != NULL)
    {
      description_length += 3;
    }
  }
  description = (CharPtr) MemNew (description_length + 1);
  if (description == NULL) return NULL;
  if (protein_name != NULL)
  {
    if (gene_name != NULL)
    {
      sprintf (description, "%s (%s)", protein_name, gene_name);
    }
    else
    {
      sprintf (description, "%s", protein_name);
    }
  }
  else
  {
    if (gene_name != NULL)
      sprintf (description, gene_name);
  }
  if (protein_name != NULL) MemFree (protein_name);
  if (StringHasNoText (description)) {
    description = MemFree (description);
  }
  return description;
}

/* This array of match functions is used to identify, in order of preference,
 * the features that might be used to generate a product for a gene-protein
 * description if the feature has not already been grouped with a product
 * feature.
 */
static matchFunction productfeatures[] = {
  IsCDS, IsmRNA, IstRNA
}; 

/* This function finds gene features without products and looks for 
 * features that might provide products for them.
 */
static void FindGeneProducts 
( ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp, productfcp;
  Int4             i, NumProductFeatureTypes;
  Boolean          suppress_locus_tag = (rp == NULL ? FALSE : rp->suppress_locus_tags);

  NumProductFeatureTypes = sizeof (productfeatures) / sizeof (matchFunction);

  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && fcp->featlist != NULL)
    {
      if (fcp->featlist->choice == DEFLINE_FEATLIST
        && IsGene (fcp->featlist->data.ptrvalue)
        && fcp->feature_label_data.productname == NULL)
      {
        productfcp = NULL;
        for (i=0; i < NumProductFeatureTypes && productfcp == NULL; i++)
        {
          productfcp = FindProductInFeatureList (fcp, clause_list,
                                                 productfeatures[i],
                                                 suppress_locus_tag);
        }
        if (productfcp != NULL)
        {
          fcp->is_alt_spliced = productfcp->is_alt_spliced;
          if (productfcp->feature_label_data.productname != NULL)
          {
            fcp->feature_label_data.productname =
                  StringSave (productfcp->feature_label_data.productname);
          }
          else
          {
            fcp->feature_label_data.productname
                  = GetProductName (productfcp->featlist->data.ptrvalue, 
                                    bsp, rp);
          }
          if (fcp->feature_label_data.description != NULL)
          {
            MemFree (fcp->feature_label_data.description);
            fcp->feature_label_data.description = NULL;
          }
          fcp->feature_label_data.description =
            GetGeneProtDescription (fcp, bsp, rp);
        }
      }
      else
      {
        FindGeneProducts (fcp->featlist, bsp, rp);
      }
    }
  }
}

static Boolean ShowInterval (
  SeqFeatPtr sfp
)
{
  if (IsSatelliteSequence (sfp) || IsExon (sfp) || IsIntron (sfp)
    || IsPromoter (sfp) || Is3UTR (sfp) || Is5UTR (sfp))
    return FALSE;
  return TRUE;
}

static CharPtr GetExonDescription (
  BioseqPtr bsp,
  SeqFeatPtr sfp
)
{
  SeqMgrFeatContext  context;
  SeqFeatPtr new_sfp;
  CharPtr    label;

  if ((new_sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, bsp, 0, 0, sfp, &context)) != sfp
      || context.label == NULL)
  {
    if ((new_sfp = SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &context)) != sfp
      || context.label == NULL)
    {
      return NULL;
    }
  }
  if ((IsExon (sfp) && StringCmp (context.label, "exon") == 0)
    || (IsIntron (sfp) && StringCmp (context.label, "intron") == 0))
  {
    return NULL;
  }
  
  label = StringSave (context.label);
  return label;
}

static CharPtr GetFeatureDescription 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr    sfp;

  if ( fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return NULL;
  }
  sfp = fcp->featlist->data.ptrvalue;
  if (sfp == NULL) return NULL;

  if (IsExon (sfp) || IsIntron (sfp))
  {
    return GetExonDescription (bsp, sfp);
  }
  else if (IsEndogenousVirusSequence (sfp))
  {
    return GetEndogenousVirusSequenceDescription (sfp);
  }
  else if (IsEndogenousVirusSourceFeature (sfp))
  {
    return GetEndogenousVirusSourceFeatureDescription (sfp);
  }
  else if (IsControlRegion (sfp))
  {
    return NULL;
  }
  else if (IsDloop (sfp))
  {
    return NULL;
  }
  else if (Is3UTR (sfp))
  {
    return NULL;
  }
  else if (Is5UTR (sfp))
  {
    return NULL;
  }
  else if (IsLTR (sfp))
  {
    return GetLTRDescription (sfp);
  }
  else
  {
    return GetGeneProtDescription (fcp, bsp, rp);
  }
}

static void LIBCALLBACK GetSatelliteFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    semicolon, colon;
  GBQualPtr  qual;
  Boolean    found = FALSE;

  flp->description = NULL; 
  flp->typeword = StringSave ("sequence");
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;
 
  if (featlist == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (main_feat == NULL) return;
  for (qual = main_feat->qual; qual != NULL && !found; qual = qual->next)
  {
    if (StringCmp (qual->qual, "satellite") == 0) 
    {
      flp->description = StringSave (qual->val);
      if ((semicolon = StringStr (flp->description, ";")) != NULL)
      {
        *semicolon = 0;
      }
      if ((colon = StringChr (flp->description, ':')) != NULL)
      {
        *colon = ' ';
      }
    }
  }
}

static void LIBCALLBACK GetPromoterFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  
  flp->description = NULL;
  flp->typeword = StringSave ("promoter region");

  if (featlist == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (main_feat == NULL) return;

  flp->description =  NULL;
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;

}

/* This function temporarily removes a 3' UTR clause from the end of
 * a clause list so that it will not be included in the list of subfeatures
 * before a CDS in the definition line.
 * The 3' UTR clause should be put back if it was not the only clause in the
 * list.
 */
static ValNodePtr Remove3UTRFromEndOfFeatList (ValNodePtr PNTR featlist)
{
  ValNodePtr vnp, prev = NULL;
  
  if (featlist == NULL || *featlist == NULL) return NULL;
  
  for (vnp = *featlist; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    prev = vnp;
  }
  if (vnp->choice == DEFLINE_CLAUSEPLUS && Is3UTRClause (vnp->data.ptrvalue))
  {
    if (prev == NULL)
    {
      *featlist = NULL;
    }
    else
    {
      prev->next = NULL;        
    }
  }
  else
  {
    vnp = NULL;
  }
  return vnp;
}

static Uint1 GetMoleculeType (BioseqPtr bsp, Uint2     entityID);
static void ConsolidateClauses (
  ValNodePtr PNTR list,
  BioseqPtr  bsp,
  Uint1      biomol,
  Boolean    delete_now,
  DeflineFeatureRequestListPtr rp);


/* This function calculates the "interval" for a clause in the definition
 * line.  The interval could be an empty string, it could indicate whether
 * the location of the feature is partial or complete and whether or not
 * the feature is a CDS, the interval could be a description of the
 * subfeatures of the clause, or the interval could be a combination of the
 * last two items if the feature is a CDS.
 */
static CharPtr GetGenericInterval 
( FeatureClausePtr fcp,
  Uint1            biomol,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr    interval;
  Boolean    partial5, partial3;
  SeqFeatPtr sfp;
  ValNodePtr featlist, strings, prev_feat;
  CharPtr    subfeatlist;
  Int4       len;
  Boolean    suppress_final_and;
  ValNodePtr utr3vnp = NULL;
  ValNodePtr last_feat;
  Uint1      molecule_type;

  if ( fcp == NULL || fcp->featlist == NULL) return NULL;
  if (fcp->is_unknown) return NULL;
  featlist = fcp->featlist;
  sfp = featlist->data.ptrvalue;
  if (sfp == NULL) return NULL;
  if ( IsExon (sfp) && fcp->is_alt_spliced)
  {
    return StringSave ("alternatively spliced");
  }
  if ( ! ShowInterval (sfp)) return NULL;

  if (IsIntergenicSpacer (sfp) && StringNCmp (sfp->comment, "may contain ", 12) == 0) {
    return StringSave ("region");
  }

  subfeatlist = NULL;
  len = 50;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  strings = NULL;
  prev_feat = NULL;
  while (featlist != NULL && featlist->choice != DEFLINE_CLAUSEPLUS)
  {
    prev_feat = featlist;
    featlist = featlist->next;
  }
  if (IsCDS (sfp))
  {
    utr3vnp = Remove3UTRFromEndOfFeatList (&featlist);
  }
  if (featlist != NULL)
  {
    suppress_final_and = FALSE;
    if (( IsCDS (sfp) && ! fcp->clause_info_only)
        || utr3vnp != NULL)
    {
      suppress_final_and = TRUE;
    }
    LabelClauses (featlist, biomol, bsp, rp);
    
    molecule_type = GetMoleculeType (bsp, bsp->idx.entityID);
    /* consolidate genes/proteins with the same names (usually hypothetical proteins) */
    ConsolidateClauses (&featlist, bsp, molecule_type, TRUE,
                        rp);

    /* make sure featlist is still intact - may have consolidated it */
    if (prev_feat == NULL)
    {
      fcp->featlist = featlist;
    }
    else
    {
      prev_feat->next = featlist;
    }

    ListClauses (featlist, &strings, FALSE, suppress_final_and);
    subfeatlist = MergeValNodeStrings (strings, FALSE);
	  ValNodeFreeData (strings);
    len += StringLen (subfeatlist) + 7;
    
    if (utr3vnp != NULL)
    {
      len += 14;
    }
  }

  interval = (CharPtr) MemNew (len * sizeof (Char));
  if (interval == NULL) return NULL;
  interval[0] = 0;

  if (StringDoesHaveText (subfeatlist))
  {
    StringCat (interval, subfeatlist);
    if ( ! IsCDS (sfp) || fcp->clause_info_only)
    {
      if (utr3vnp != NULL)
      {
        if (featlist->next != NULL)
        {
          StringCat (interval, ",");
        }
        StringCat (interval, " and 3' UTR");
        /* put 3' UTR back at end of featlist */
        if (featlist != NULL)
        {
          last_feat = featlist;
          while (last_feat != NULL && last_feat->next != NULL)
          {
            last_feat = last_feat->next;
          }
          last_feat->next = utr3vnp;
        }
      }
      if (subfeatlist != NULL) MemFree (subfeatlist);
      return interval;
    }
    if (utr3vnp == NULL)
    {
      StringCat (interval, " and ");
    }
    else
    {
      StringCat (interval, ", ");
    }
  }

  if (FeatureDoesNotGetPartialComplete (sfp)) 
  {
    /* don't add partial or complete */
  } 
  else if (partial5 || partial3)
  {
    StringCat (interval, "partial ");
  }
  else
  {
    StringCat (interval, "complete ");
  }
  if (IsCDS (sfp) && ! IsPseudo (sfp))
  {
    StringCat (interval, "cds");
    if (fcp->is_alt_spliced)
      StringCat (interval, ", alternatively spliced");
  }
  else
  {
    StringCat (interval, "sequence");
    if (IsNoncodingProductFeat (sfp) && fcp->is_alt_spliced)
    {
      StringCat (interval, ", alternatively spliced");
    }
  }
  
  if (utr3vnp != NULL)
  {
    /* tack UTR3 on at end of clause */
    if (StringDoesHaveText (subfeatlist))
    {
      StringCat (interval, ",");
    }
    StringCat (interval, " and 3' UTR");
    
    /* put 3' UTR back at end of featlist */
    if (featlist != NULL)
    {
      last_feat = featlist;
      while (last_feat != NULL && last_feat->next != NULL)
      {
        last_feat = last_feat->next;
      }
      last_feat->next = utr3vnp;
    }
  }
  
  if (subfeatlist != NULL) MemFree (subfeatlist);    
  
  return interval;
} 


/* This function is used to generate feature label information for
 * a feature clause.  It is called by the LabelFeature function if
 * a "GetFeatureLabel" function is not found for the specific feature
 * type.
 * In the future it may be advisable to create "GetFeatureLabel" functions
 * for more of the specific feature types, to reduce the number of times
 * that the feature must be identified as being a certain type.
 */ 
static void LIBCALLBACK GetGenericFeatureLabel 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  Uint1            biomol,
  FeatureLabelPtr  flp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr main_feat;
  
  if (fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return;
  }
  main_feat = fcp->featlist->data.ptrvalue;
  if (main_feat == NULL) return;

  if (flp->typeword == NULL)
  {
    flp->typeword = GetFeatureTypeWord (biomol, main_feat);
    flp->is_typeword_first = IsTypeWordFirst (flp->typeword);
    flp->pluralizable = TRUE;
  }
  if (flp->productname == NULL)
  {
    flp->productname = GetProductName (main_feat, bsp, rp);
  }
  if (flp->description == NULL
    && (! IsMiscRNA (main_feat)
      || StringStr (flp->productname, "spacer") == NULL ))
  {
    flp->description = GetFeatureDescription (fcp, bsp, rp);
  }

}

typedef void (LIBCALLBACK *GetFeatureLabelFunction) (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
);

typedef struct matchlabelfunction {
  matchFunction           itemmatch;
  GetFeatureLabelFunction labelfunction;
} MatchLabelFunctionData, PNTR MatchLabelFunctionPtr;

static MatchLabelFunctionData label_functions[] = {
 { IsSatelliteSequence, GetSatelliteFeatureLabel         },
 { IsMobileElement,     GetMobileElementFeatureLabel        },
 { IsPromoter,          GetPromoterFeatureLabel          },
 { IsIntergenicSpacer,  GetIntergenicSpacerFeatureLabel  },
 { IsGeneCluster,       GetGeneClusterFeatureLabel       }
};

typedef enum {
 DEFLINE_FEATLABEL_Satellite = 0,
 DEFLINE_FEATLABEL_Transposon,
 DEFLINE_FEATLABEL_Promoter,
 DEFLINE_FEATLABEL_IntergenicSpacer,
 DEFLINE_FEATLABEL_GeneCluster,
 NumDefLineFeatLabels
} DefLineFeatLabel;

static void LabelFeature 
( BioseqPtr        bsp,
  Uint1            biomol,
  FeatureClausePtr new_clauseplus,
  DeflineFeatureRequestListPtr rp)
{
  Int4             i;
  SeqFeatPtr       main_feat;

  if (new_clauseplus == NULL || new_clauseplus->featlist == NULL) return;

  if (new_clauseplus->featlist->choice == DEFLINE_FEATLIST)
  {
    main_feat = (SeqFeatPtr) new_clauseplus->featlist->data.ptrvalue;
    
    new_clauseplus->allelename = GetAlleleName (new_clauseplus->grp,
                                                rp == NULL ? FALSE : rp->suppress_locus_tags);
    if (new_clauseplus->interval == NULL)
    {
      new_clauseplus->interval =
                  GetGenericInterval (new_clauseplus, biomol, bsp, rp);
    }

    for (i=0; i < NumDefLineFeatLabels; i++)
    {
      if (label_functions [i].itemmatch (main_feat))
      {
        label_functions [i].labelfunction ( new_clauseplus->featlist,
                                          bsp, biomol,
                                          &new_clauseplus->feature_label_data);
        return;
      }
    }

    GetGenericFeatureLabel ( new_clauseplus, bsp, biomol, 
                           &new_clauseplus->feature_label_data, rp);
    return;
  }
}

/* This function is used to calculate the parts of a product name that
 * are "the same" for use as the name of an alternatively spliced product.
 * The common portion of the string must end at a recognized separator,
 * such as a space, comma, or dash instead of in the middle of a word.
 * The matching portions of the string could occur at the beginning or end
 * of the string, or even occasionally at the beginning and end of a
 * string, but not as the center of the string with a different beginning
 * and ending.
 */
static CharPtr FindStringIntersection (
  CharPtr str1,
  CharPtr str2,
  Boolean str1_previously_stripped
)
{
  Int4 matchleftlen = 0;
  Int4 matchlefttoken = 0;
  Int4 matchrightidx1 = 0;
  Int4 matchrightidx2 = 0;
  Int4 matchrighttoken = 0;
  CharPtr match_string;
  Int4 len1;
  Int4 len2;
  Int4 match_len;

  if (str1 == NULL || str2 == NULL) return NULL;
  if (StringCmp (str1, str2) == 0) return StringSave (str1);
  len1 = StringLen (str1);
  len2 = StringLen (str2);

  while (str1[matchleftlen] != 0 && str2[matchleftlen] != 0
         && str1[matchleftlen] == str2[matchleftlen])
  {
    if (str1 [matchleftlen] == ','
      || str1 [matchleftlen] == '-')
    {
      matchlefttoken = matchleftlen;
    }
    else if (str1 [matchleftlen] == ' '
      && matchlefttoken != matchleftlen - 1)
    {
      matchlefttoken = matchleftlen;
    }
    matchleftlen++;
  }
  if (matchleftlen == len1 && str1_previously_stripped) 
  {
    matchlefttoken = matchleftlen;
  }
  else
  {
    matchleftlen = matchlefttoken;
  }

  matchrightidx1 = len1;
  matchrightidx2 = len2;

  while (matchrightidx1 > -1 && matchrightidx2 > -1
         && str1[matchrightidx1] == str2[matchrightidx2])
  {
    if (str1 [matchrightidx1] == ' '
      || str1[matchrightidx1] == ','
      || str1[matchrightidx1] == '-')
    {
      matchrighttoken = matchrightidx1;
    }
    matchrightidx1--;
    matchrightidx2--;
  }
  if (matchrightidx1 == -1)
  {
    matchrighttoken = matchrightidx1;
  }
  else if (matchrighttoken > 0) 
  {
    matchrightidx1 = matchrighttoken;
  } 
  else if (str1_previously_stripped && matchrightidx1 < len1 - 1)
  {
    /* matchrightidx1 = matchrighttoken; */
    /* do nothing, leave right index where it is */
  }
  else
  {
    matchrightidx1 = len1;
  }

  match_len = matchleftlen;
  if (matchrightidx1 < len1 - 1)
  {
    match_len += len1 - matchrightidx1 - 1;
  }

  if (match_len <= 0) return NULL;

  match_string = MemNew (match_len + 2);
  if (match_string == NULL) return NULL;
  if (matchleftlen != 0)
  {
    StringNCpy (match_string, str1, matchleftlen);
    match_string[matchleftlen] = 0;
  }
  else
  {
    match_string[0] = 0;
  }
  if (matchrightidx1 < len1)
  {
    if (match_string[0] != 0) StringCat (match_string, " ");
    StringCat (match_string, str1 + matchrightidx1 + 1);
  }
  return match_string;
}

/* These are the words that are used to introduced the part of the protein
 * name that differs in alt-spliced products - they should not be part of
 * the alt-spliced product name.
 * Note that splice variant is listed before "variant" so that it will be
 * found first and "variant" will not be removed from "splice variant", leaving
 * splice as an orphan.
 */

static CharPtr UnwantedWords [] = {
 "splice variant",
 "splice product",
 "variant",
 "isoform"
};

static void TrimUnwantedWordsFromAltSpliceProductName (
  CharPtr productname
)
{
  Int4    num_unwanted_words, i;
  size_t  unwanted_word_len;
  CharPtr cp, tmp;
  
  num_unwanted_words = sizeof (UnwantedWords) / sizeof (CharPtr);
  for (i = 0; i < num_unwanted_words; i++)
  {
    unwanted_word_len = StringLen (UnwantedWords [i]);
    cp = StringStr (productname, UnwantedWords [i]);
    if (cp != NULL)
    {
      if (cp == productname)
      {
        /* word occurs in beginning of phrase */
        tmp = StringSave (productname + unwanted_word_len);
        StringCpy (productname, tmp);
        MemFree (tmp);
      }
      else if (cp - productname < StringLen (productname) - unwanted_word_len)
      {
        /* word occurs in middle of phrase */
        tmp = StringSave (cp + unwanted_word_len);
        StringCpy (cp - 1, tmp);
        MemFree (tmp);
      }
      else
      {
        /* word occurs at end of phrase */
        *cp = 0;
      }
    }
  }
}


static Boolean PreviouslyStripped (SeqFeatPtr cds, BioseqPtr bsp, CharPtr productname)
{
  CharPtr expected_product_name;
  Boolean rval = FALSE;
  
  if (cds == NULL || StringHasNoText (productname)) return FALSE;
  expected_product_name = GetProductName (cds, bsp, FALSE);
  if (StringCmp (productname, expected_product_name) != 0) {
    rval = TRUE;
  }
  expected_product_name = MemFree (expected_product_name);
  return rval;
}

/* This function determines whether two CDSs meet the conditions for
 * alternative splicing, and if so, it returns the name of the alternatively
 * spliced product.  In order to be alternatively spliced, the two CDSs 
 * must have the same gene, must share a complete interval, and must have
 * similarly named products.
 */
static CharPtr MeetAltSpliceRules 
( FeatureClausePtr cdsfcp1,
  FeatureClausePtr cdsfcp2,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr cds1, cds2;
  CharPtr match_string;
  Int4    res;
  
  if (cdsfcp1 == NULL || cdsfcp2 == NULL
    || cdsfcp1->featlist == NULL || cdsfcp2->featlist == NULL)
  {
    return NULL;
  }

  cds1 = cdsfcp1->featlist->data.ptrvalue;
  cds2 = cdsfcp2->featlist->data.ptrvalue;
  if (! DoGenesMatch (cdsfcp1->grp, cdsfcp2->grp, rp == NULL ? FALSE : rp->suppress_locus_tags))
    return NULL;

  if ( (res = TestFeatOverlap (cds1, cds2, COMMON_INTERVAL)) != -1)
  {
    match_string = FindStringIntersection (
                     cdsfcp1->feature_label_data.productname,
                     cdsfcp2->feature_label_data.productname,
                     PreviouslyStripped(cds1, bsp, cdsfcp1->feature_label_data.productname));
    return match_string;
  }
  return NULL;
}

/* This function is used by the FindAltSplices function to locate the
 * next CDS in a list of feature clauses.
 */
static ValNodePtr FindNextCDSClause (ValNodePtr vnp)
{
  FeatureClausePtr fcp;

  while (vnp != NULL)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp != NULL && !fcp->delete_me && fcp->featlist != NULL
        && IsCDS (fcp->featlist->data.ptrvalue))
      {
        return vnp;
      }
    }
    vnp = vnp->next;
  }
  return NULL;
}

/* This function is used by the FindAltSplices function to move the features
 * and subclauses from the second CDS in an alternatively spliced pair of 
 * CDSs to the feature clause for the first CDS, so that the subfeatures
 * can be properly listed.
 */
static void MoveSubclauses (
  FeatureClausePtr dstfcp,
  FeatureClausePtr srcfcp
)
{
  ValNodePtr dst_last_feat, dst_first_clause, dst_last_clause;
  ValNodePtr src_last_feat, src_first_clause;

  if (dstfcp == NULL || srcfcp == NULL || srcfcp->featlist == NULL) return;

  dst_first_clause = NULL;
  dst_last_clause = NULL;
  src_first_clause = NULL;

  dst_last_feat = dstfcp->featlist;
  while (dst_last_feat != NULL 
      && dst_last_feat->next != NULL
      && dst_last_feat->next->choice == DEFLINE_FEATLIST)
  {
    dst_last_feat = dst_last_feat->next;
  }
  if (dst_last_feat != NULL)
  {
    dst_first_clause = dst_last_feat->next;
  }
  dst_last_clause = dst_first_clause;
  while (dst_last_clause != NULL && dst_last_clause->next != NULL)
  {
    dst_last_clause = dst_last_clause->next;
  }

  src_last_feat = srcfcp->featlist;
  while (src_last_feat != NULL 
      && src_last_feat->next != NULL
      && src_last_feat->next->choice == DEFLINE_FEATLIST)
  {
    src_last_feat = src_last_feat->next;
  }
  if (src_last_feat != NULL)
  {
    src_first_clause = src_last_feat->next;
  }

  /* insert features before clauses */
  if (dst_last_feat == NULL)
  {
    dstfcp->featlist = srcfcp->featlist;
    dst_last_feat = src_last_feat;
  }
  else
  {
    dst_last_feat->next = srcfcp->featlist;
  }
  /* insert clauses after feats */
  if (dst_first_clause != NULL)
  {
    src_last_feat->next = dst_first_clause;
    dst_last_clause->next = src_first_clause;
  }
  srcfcp->featlist = NULL;
}
 
/* we want to look through the list for CDS features */
/* if we find two CDSs that are alternatively spliced, */
/* we replace the first alternatively spliced CDS feature */
/* with a new CDS feature that has the new protein name as */
/* a comment and a data.choice value that indicates alt splicing */
/* we remove the second alternatively spliced CDS feature from the list */

static void FindAltSplices 
( ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr  fcp1, fcp2;
  ValNodePtr cdsclause1, cdsclause2;
  ValNodePtr searchclause;
  CharPtr  combined_protein_name;
  Boolean    partial3_1, partial5_1, partial3_2, partial5_2;
  Int4       left1, left2, right1, right2;

  if (clause_list == NULL) return;
  
  cdsclause1 = FindNextCDSClause (clause_list);
  while (cdsclause1 != NULL)
  {
    fcp1 = (FeatureClausePtr) cdsclause1->data.ptrvalue;
    if (fcp1->feature_label_data.productname == NULL)
    {
      fcp1->feature_label_data.productname = 
           GetProductName (fcp1->featlist->data.ptrvalue, bsp, rp);
    }
    searchclause = cdsclause1->next;
    cdsclause2 = FindNextCDSClause (searchclause);
    while (cdsclause2 != NULL) 
    {
      fcp2 = (FeatureClausePtr) cdsclause2->data.ptrvalue;
      if (fcp2->feature_label_data.productname == NULL)
      {
        fcp2->feature_label_data.productname =
           GetProductName (fcp2->featlist->data.ptrvalue, bsp, rp);
      }
      combined_protein_name = MeetAltSpliceRules (fcp1, fcp2, bsp, rp);
      if (combined_protein_name != NULL)
      {
        /* get rid of variant, splice variant, splice product, isoform, etc.*/
        TrimUnwantedWordsFromAltSpliceProductName (combined_protein_name);

        /* get rid of trailing spaces in protein name */
        TrimSpacesAroundString (combined_protein_name);

        /* copy new protein name into first clause */
        MemFree (fcp1->feature_label_data.productname);
        fcp1->feature_label_data.productname = combined_protein_name;
        CheckSeqLocForPartial (fcp1->slp, &partial5_1, &partial3_1);
        left1 = GetOffsetInBioseq (fcp1->slp, bsp, SEQLOC_LEFT_END);
        right1 = GetOffsetInBioseq (fcp1->slp, bsp, SEQLOC_RIGHT_END);
        CheckSeqLocForPartial (fcp2->slp, &partial5_2, &partial3_2);
        left2 = GetOffsetInBioseq (fcp2->slp, bsp, SEQLOC_LEFT_END);
        right2 = GetOffsetInBioseq (fcp2->slp, bsp, SEQLOC_RIGHT_END);
        fcp1->slp = SeqLocMerge (bsp, fcp1->slp, fcp2->slp,
                                 FALSE, TRUE, FALSE);
        if (left1 == left2)
        {
          partial5_1 |= partial5_2;
        }
        else
        {
          partial5_1 = left1 < left2 ? partial5_1 : partial5_2;
        }
        if (right1 == right2)
        {
          partial3_1 |= partial3_2;
        }
        else
        {
          partial3_1 = right1 > right2 ? partial3_1 : partial3_2;
        }
        SetSeqLocPartial (fcp1->slp, partial5_1, partial3_1);
        fcp1->is_alt_spliced = TRUE;

        /* copy over fcp2 subclauses */
        MoveSubclauses (fcp1, fcp2);

        /* remove second clause */
        fcp2->delete_me = TRUE;
      }
      searchclause = cdsclause2->next;
      cdsclause2 = FindNextCDSClause (searchclause);
    }
    cdsclause1 = FindNextCDSClause (cdsclause1->next);
  } 
  DeleteFeatureClauses (&clause_list);
}

static void LabelClauses 
( ValNodePtr clause_list,
  Uint1      biomol,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr clause;
 
  clause = clause_list;
  while (clause != NULL)
  { 
    LabelFeature ( bsp, biomol, clause->data.ptrvalue, rp);
    clause = clause->next;
  }
}

static CharPtr misc_words [] = {
  "internal transcribed spacer",
  "external transcribed spacer",
  "ribosomal RNA intergenic spacer",
  "ribosomal RNA",
  "intergenic spacer"
};

typedef enum {
  MISC_RNA_WORD_INTERNAL_SPACER = 0,
  MISC_RNA_WORD_EXTERNAL_SPACER,
  MISC_RNA_WORD_RNA_INTERGENIC_SPACER,
  MISC_RNA_WORD_RNA,
  MISC_RNA_WORD_INTERGENIC_SPACER,
  NUM_MISC_RNA_WORDS
} MiscWord;

/* note - must put substrings of other separators after the longer version */
static CharPtr separators [] = {
  ", and ",
  " and ",
  ", ",
  "; "
};

#define num_separators 3


static ValNodePtr TokenListFromMiscRNAString (CharPtr str)
{
  ValNodePtr token_list = NULL;
  CharPtr cansep [num_separators];
  CharPtr token_start, next_sep, token;
  Int4    i, sep_len, datalen;
  Uint1   word_i;
  Boolean found_unparseable = FALSE;

  if ( StringStr (str, "spacer") == NULL) {
    return NULL;
  }

  token_start = str;
  for (i = 0; i < num_separators; i++) {
    cansep[i] = StringStr (token_start, separators[i]);
  }

  while (*token_start != 0 && !found_unparseable) {
    next_sep = NULL;
    sep_len = 0;
    for (i = 0; i < num_separators; i++) {
      if (cansep[i] != NULL) {
        if (cansep[i] < token_start) {
          cansep[i] = StringStr (token_start, separators[i]);
        }
      }
      if (cansep[i] != NULL && (next_sep == NULL || next_sep > cansep[i])) {
        next_sep = cansep[i];
        sep_len = StringLen (separators[i]);
      }
    }
    if (next_sep == NULL) {
      token = StringSave (token_start);
      datalen = StringLen (token);
    } else {
      datalen = next_sep - token_start;
      token = (CharPtr) MemNew (sizeof (Char) * (datalen + 1));
      StringNCpy (token, token_start, datalen);
      token[datalen] = 0;
    }
    /* determine which word is part of the token */
    for (word_i=0;
         word_i < NUM_MISC_RNA_WORDS
           && StringStr (token, misc_words [word_i]) == NULL;
         word_i++) {}
    if (word_i < NUM_MISC_RNA_WORDS) {
      ValNodeAddPointer (&token_list, word_i, token);
    } else {
      found_unparseable = TRUE;
    }      
    token_start += datalen + sep_len;
  }
  if (found_unparseable) {
    token_list = ValNodeFreeData (token_list);
  }
  return token_list;
}


static ValNodePtr 
GetFeatureClausesFromMiscRNATokens 
( ValNodePtr token_list,
  SeqFeatPtr misc_rna,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp) 
{
  ValNodePtr clause_list = NULL;
  ValNodePtr vnp;
  Boolean    partial5, partial3, unparseable = FALSE;
  CharPtr    word_loc;
  FeatureClausePtr fcp;

  if (token_list == NULL || misc_rna == NULL) {
    return NULL;
  }

  CheckSeqLocForPartial (misc_rna->location, &partial5, &partial3);
  
  for (vnp = token_list; vnp != NULL && !unparseable; vnp = vnp->next) {
    word_loc = StringStr (vnp->data.ptrvalue, misc_words [vnp->choice]);
    if (word_loc == NULL) {
      unparseable = TRUE;
    } else {
      fcp = NewFeatureClause ( misc_rna, bsp, rp);
      if (fcp == NULL) {
        unparseable = TRUE;
      } else {
        if (vnp->choice == MISC_RNA_WORD_INTERNAL_SPACER
            || vnp->choice == MISC_RNA_WORD_EXTERNAL_SPACER
            || vnp->choice == MISC_RNA_WORD_RNA_INTERGENIC_SPACER
            || vnp->choice == MISC_RNA_WORD_INTERGENIC_SPACER) {
          if (word_loc == vnp->data.ptrvalue) {
            fcp->feature_label_data.is_typeword_first = TRUE;
            fcp->feature_label_data.typeword = StringSave (misc_words [vnp->choice]);
            if (StringLen (misc_words [vnp->choice]) + 1 < StringLen (vnp->data.ptrvalue)) {
              fcp->feature_label_data.description =
                    StringSave ( ((CharPtr)vnp->data.ptrvalue) + StringLen (misc_words [vnp->choice]) + 1);
            }
          } else {
            fcp->feature_label_data.is_typeword_first = FALSE;
            fcp->feature_label_data.typeword = StringSave (misc_words [vnp->choice]);
            if (StringLen (misc_words [vnp->choice]) + 1 < StringLen (vnp->data.ptrvalue)) {
              fcp->feature_label_data.description = StringSave ( vnp->data.ptrvalue);
              fcp->feature_label_data.description [word_loc - ((CharPtr) vnp->data.ptrvalue) - 1] = 0;
            }
          }
        } else if (vnp->choice == MISC_RNA_WORD_RNA) {
          fcp->feature_label_data.description = StringSave (vnp->data.ptrvalue);
        }
        if ((vnp == token_list && partial5) || (vnp->next == NULL && partial3)) {
          fcp->interval = StringSave ("partial sequence");
        } else {
          fcp->interval = StringSave ("complete sequence");
        }
        ValNodeAddPointer (&clause_list, DEFLINE_CLAUSEPLUS, fcp);
      }
    }
  }
  if (unparseable) {
    DefLineFeatClauseListFree (clause_list);
    clause_list = NULL;
  }
  return clause_list;
}


static ValNodePtr GetMiscRNAelements 
( SeqFeatPtr misc_rna,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr buffer;
  ValNodePtr token_list, clause_list = NULL;
  FeatureClausePtr fcp;


  if (misc_rna == NULL) return NULL;
  buffer = GetProductName (misc_rna, bsp, rp);
  if (buffer == NULL) 
  {
    buffer = StringSave (misc_rna->comment);
  }
  else if (StringNCmp (buffer, misc_rna->comment, StringLen (buffer) -1) == 0
    && buffer [ StringLen (buffer) - 1] == '>')
  {
    MemFree (buffer);
    buffer = StringSave (misc_rna->comment);
  }

  if (StringNCmp (buffer, "contains ", 9) == 0) {
    token_list = TokenListFromMiscRNAString (buffer + 9);
    clause_list = GetFeatureClausesFromMiscRNATokens (token_list, misc_rna, bsp, rp);
    token_list = ValNodeFreeData (token_list);
  } else if (StringNCmp (buffer, "may contain ", 12) == 0) {
    token_list = TokenListFromMiscRNAString (buffer + 12);
    if (token_list != NULL) {
      fcp = NewFeatureClause ( misc_rna, bsp, rp);
      fcp->feature_label_data.description = StringSave (buffer + 12);
      fcp->interval = StringSave ("region");
      ValNodeAddPointer (&clause_list, DEFLINE_CLAUSEPLUS, fcp);
    }
    token_list = ValNodeFreeData (token_list);
  } else {
    token_list = TokenListFromMiscRNAString (buffer);
    clause_list = GetFeatureClausesFromMiscRNATokens (token_list, misc_rna, bsp, rp);
    token_list = ValNodeFreeData (token_list);
  }

  buffer = MemFree (buffer);
  return clause_list;
}

/* Some misc_RNA clauses have a comment that actually lists multiple
 * features.  This function creates a clause for each element in the
 * comment and inserts the list of new clauses into the feature list
 * at the point where the single previous clause was.
 */
static void ReplaceRNAClauses (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr main_feat;
  ValNodePtr clause, replacement_clauses, nextclause, vnp;

  if (clause_list == NULL || *clause_list == NULL) return;
  clause = *clause_list;
  while (clause != NULL)
  {
    nextclause = clause->next;
    fcp = (clause->data.ptrvalue);
    if (fcp == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      return;
    }
    main_feat = (SeqFeatPtr) fcp->featlist->data.ptrvalue;
  
    if (IsrRNA (main_feat) || IsMiscRNA (main_feat))
    {
      replacement_clauses = GetMiscRNAelements ( main_feat, bsp, rp );
      if (replacement_clauses != NULL)
      {
        for (vnp = replacement_clauses; vnp->next != NULL; vnp = vnp->next) {}
        vnp->next = clause->next;
        clause->next = replacement_clauses;
        fcp->delete_me = TRUE;
      }
    }
    clause = nextclause;
  }
  DeleteFeatureClauses (clause_list);
}

/* Some misc_feat clauses have a comment that lists one or more tRNAs and
 * an intergenic spacer.  This function creates a clause for each element 
 * in the comment and inserts the list of new clauses into the feature list
 * at the point where the single previous clause was.
 */
static void ReplaceIntergenicSpacerClauses (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr main_feat;
  ValNodePtr clause, replacement_clauses, nextclause, vnp;

  if (clause_list == NULL || *clause_list == NULL) return;
  clause = *clause_list;
  while (clause != NULL)
  {
    nextclause = clause->next;
    fcp = (clause->data.ptrvalue);
    if (fcp == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      return;
    }
    main_feat = (SeqFeatPtr) fcp->featlist->data.ptrvalue;
  
    if (IsParsableList (main_feat)) 
    {
      if ((replacement_clauses = ParsetRNAIntergenicSpacerElements ( main_feat, bsp, rp)) != NULL)
      {
        for (vnp = replacement_clauses; vnp->next != NULL; vnp = vnp->next) {}
        vnp->next = clause->next;
        clause->next = replacement_clauses;
        fcp->delete_me = TRUE;
      }
      else
      {
        fcp->delete_me = TRUE;
      }
    }
    clause = nextclause;
  }
  DeleteFeatureClauses (clause_list);
}

/* If we are applying a different rule for misc_feats, we need to recalculate
 * their descriptions.
 */
static void RenameMiscFeats (ValNodePtr clause_list, Uint1 biomol)
{
  ValNodePtr       vnp, featlist;
  FeatureClausePtr fcp, featlistclause;
  SeqFeatPtr       sfp;
  Int4             name_len;

  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS || vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = vnp->data.ptrvalue;
    for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
    {
      if ( featlist->data.ptrvalue == NULL)
      {
        continue;
      }
      if (featlist->choice == DEFLINE_CLAUSEPLUS)
      {
        featlistclause = featlist->data.ptrvalue;
        RenameMiscFeats (featlistclause->featlist, biomol);
        continue;
      }
      if (featlist->choice != DEFLINE_FEATLIST)
      {
        continue;
      }
      sfp = featlist->data.ptrvalue;
      if (sfp->idx.subtype != FEATDEF_misc_feature
        || sfp->comment == NULL
        || IsIntergenicSpacer (sfp)
        || IsGeneCluster (sfp)
        || IsControlRegion (sfp)) 
      {
        continue;
      }
      if (fcp->feature_label_data.description != NULL) 
      {
        fcp->feature_label_data.description
                   = MemFree (fcp->feature_label_data.description);
      }
      name_len = StringCSpn (sfp->comment, ";");
	  /* make sure we have space for terminating NULL */
      fcp->feature_label_data.description = MemNew ((name_len + 1) * sizeof (Char));
      if (fcp->feature_label_data.description == NULL) return;
      StringNCpy (fcp->feature_label_data.description, sfp->comment, name_len);
      fcp->feature_label_data.description [ name_len ] = 0;
      fcp->feature_label_data.typeword =
            MemFree (fcp->feature_label_data.typeword);
      if (biomol == MOLECULE_TYPE_GENOMIC)
      {
        fcp->feature_label_data.typeword = StringSave ("genomic sequence");
      }
      else if (biomol == MOLECULE_TYPE_MRNA)
      {
        fcp->feature_label_data.typeword = StringSave ("mRNA sequence");
      }
      else
      {
        fcp->feature_label_data.typeword = StringSave ("sequence");
      }
     
      fcp->interval = MemFree (fcp->interval);
      fcp->interval = StringSave ("");
    }
  }
}

static void RemoveUnwantedMiscFeats (
  ValNodePtr PNTR clause_list,
  Boolean delete_now
)
{
  ValNodePtr       vnp, featlist;
  FeatureClausePtr fcp, featlistclause;
  SeqFeatPtr       sfp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS || vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = vnp->data.ptrvalue;
    for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
    {
      if ( featlist->data.ptrvalue == NULL)
      {
        continue;
      }
      if (featlist->choice == DEFLINE_CLAUSEPLUS)
      {
        featlistclause = featlist->data.ptrvalue;
        RemoveUnwantedMiscFeats (&(featlistclause->featlist), FALSE);
        continue;
      }
      if (featlist->choice != DEFLINE_FEATLIST)
      {
        continue;
      }
      sfp = featlist->data.ptrvalue;
      if ( sfp->idx.subtype == FEATDEF_misc_feature
        && ! IsNoncodingProductFeat (sfp)
        && ! IsControlRegion (sfp)
        && ! IsIntergenicSpacer (sfp)
        && ! IsGeneCluster (sfp)
        && ! IsParsableList (sfp))
      {
        fcp->delete_me = TRUE;
      }
    }
  }
  DeleteFeatureClauses (clause_list);
}

/* When a feature is on the minus strand, the clauses are listed by 
 * sequence indexing in reverse biological order - we reverse the subclauses
 * for the feature in order to have them listed in the definition line
 * in biological order.
 * This is most noticeable when the main feature is a CDS with multiple
 * exons numbered sequentially.  If the exons are on the minus strand and
 * appear as 9, 8, 7, 6, we want to display them in the definition line as
 * 6, 7, 8, 9.
 */
static void ReverseClauses (
  ValNodePtr PNTR clause_list,
  matchFunction   itemmatch
)
{
  ValNodePtr vnp, last_feat, first_feat, next_item, new_list;
  FeatureClausePtr fcp;

  if (clause_list == NULL || *clause_list == NULL) return;

  last_feat = NULL;
  first_feat = NULL;
  new_list = NULL;
  vnp = *clause_list;
  while (vnp != NULL)
  {
    next_item = vnp->next;
    fcp = NULL;
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && fcp->slp != NULL
      && SeqLocStrand (fcp->slp) == Seq_strand_minus
      && fcp->featlist != NULL
      && fcp->featlist->choice == DEFLINE_FEATLIST
      && itemmatch (fcp->featlist->data.ptrvalue))
    {
      vnp->next = new_list;
      new_list = vnp;
    }
    else
    {
      if (first_feat == NULL)
      {
        first_feat = vnp;
        last_feat = vnp;
      }
      else
      {
        last_feat->next = vnp;
        last_feat = vnp;
        last_feat->next = NULL;
      }
    }
    if (fcp != NULL)
    {
      ReverseClauses (&(fcp->featlist), itemmatch);
    }
    vnp = next_item;
  }
  if (first_feat == NULL)
  {
    *clause_list = new_list;
  }
  else
  {
    last_feat->next = new_list;
    *clause_list = first_feat;
  }
}

/* This function is used to determine whether two features are both exons
 * and whether they are numerically sequential - i.e., exon 7 and exon 8
 * are a pair of consecutive exons, exon 7 and exon 9 are not, and exon 7
 * and intron 9 are not.
 */
static Boolean ClausePairIsTwoConsecutiveExons (
  ValNodePtr vnp1,
  ValNodePtr vnp2,
  BioseqPtr  bsp
)
{
  FeatureClausePtr fcp1, fcp2;
  SeqFeatPtr       exon1, exon2;
  Int4 num1, num2;
  CharPtr          exdesc1, exdesc2;

  if (vnp1 == NULL || vnp2 == NULL 
    || vnp1->choice != DEFLINE_CLAUSEPLUS
    || vnp2->choice != DEFLINE_CLAUSEPLUS
    || vnp1->data.ptrvalue == NULL
    || vnp2->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  fcp1 = vnp1->data.ptrvalue;
  fcp2 = vnp2->data.ptrvalue;
  if ( fcp1->featlist == NULL
    || fcp1->featlist->data.ptrvalue == NULL
    || fcp2->featlist == NULL
    || fcp2->featlist->data.ptrvalue == NULL
    || fcp1->featlist->choice != DEFLINE_FEATLIST
    || fcp2->featlist->choice != DEFLINE_FEATLIST
    || ! IsExon (fcp1->featlist->data.ptrvalue)
    || ! IsExon (fcp2->featlist->data.ptrvalue)
    || (fcp1->is_alt_spliced && ! fcp2->is_alt_spliced)
    || (! fcp1->is_alt_spliced && fcp2->is_alt_spliced))
  {
    return FALSE;
  }

  exon1 = (SeqFeatPtr)(fcp1->featlist->data.ptrvalue);
  exon2 = (SeqFeatPtr)(fcp2->featlist->data.ptrvalue);

  exdesc1 = GetExonDescription (bsp, exon1);
  exdesc2 = GetExonDescription (bsp, exon2);
  if (exdesc1 == NULL || exdesc2 == NULL)
  {
    if (exdesc1 != NULL) MemFree (exdesc1);
    if (exdesc2 != NULL) MemFree (exdesc2);
    return FALSE;
  }
  
  num1 = atoi (exdesc1);
  num2 = atoi (exdesc2);
  MemFree (exdesc1);
  MemFree (exdesc2);

  if (abs (num1 - num2) == 1)
  {
    return TRUE;
  }

  return FALSE; 
}

/* This function counts the number of consecutive exons in a list.
 */
static Int4 GetNumberOfConsecutiveExons (
  ValNodePtr list,
  BioseqPtr  bsp
)
{
  ValNodePtr check;
  Int4       num_exons;
 
  num_exons = 0;
  check = list->next;
  if ( ! ClausePairIsTwoConsecutiveExons (list, check, bsp)) return 0;
  
  num_exons = 2;
  while ( check != NULL
    && ClausePairIsTwoConsecutiveExons (check, check->next, bsp))
  {
    num_exons++;
    check = check->next;
  }
  return num_exons;
}

/* This function replaces a list of three or more consecutive exon clauses
 * with a single "summary" clause that gives the range of exons present -
 * i.e., if you have exons 1, 2, 3, and 4, a clause will be created that
 * contains all four of those features and has a description of "1 through 4".
 */
static void ReplaceExonClauseList (
  FeatureClausePtr fcp,
  ValNodePtr       clause,
  Int4             num_exons,
  BioseqPtr        bsp
)
{
  ValNodePtr       lastfeat, tmpclause;
  FeatureClausePtr tmpfcp;
  Int4             i;
  CharPtr          new_description;
  Int4             new_description_len;
  CharPtr          exdesc1, exdesc2;

  if (fcp == NULL || clause == NULL) return;

  lastfeat = fcp->featlist;
  while (lastfeat != NULL && lastfeat->next != NULL)
  {
    lastfeat = lastfeat->next;
  }
  tmpclause = clause->next;
  for (i=0; i < num_exons - 1 && tmpclause != NULL; i++)
  {
    tmpfcp = tmpclause->data.ptrvalue;
    tmpfcp->delete_me = TRUE;
    if (lastfeat == NULL)
    {
      fcp->featlist = tmpfcp->featlist;
    }
    else
    {
      lastfeat->next = tmpfcp->featlist;
    }
    tmpfcp->featlist = NULL;
    while (lastfeat != NULL && lastfeat->next != NULL)
    {
      lastfeat = lastfeat->next;
    }
          
    tmpclause = tmpclause->next;
  }

  exdesc1 = GetExonDescription (bsp, fcp->featlist->data.ptrvalue);
  exdesc2 = GetExonDescription (bsp, lastfeat->data.ptrvalue);
  if (exdesc1 == NULL || exdesc2 == NULL)
  {
    if (exdesc1 != NULL) MemFree (exdesc1);
    if (exdesc2 != NULL) MemFree (exdesc2);
    return;
  }
  new_description_len =
        StringLen (exdesc1)
      + StringLen (exdesc2)
      + StringLen (" through ")
      + 1;
  new_description = MemNew (new_description_len * sizeof (Char));
  if (new_description == NULL) return;
  sprintf (new_description, "%s through %s", exdesc1, exdesc2);
  MemFree (exdesc1);
  MemFree (exdesc2);
  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
  }
  fcp->feature_label_data.description = new_description;
}

/* This function recursively searches for lists of consecutive exons
 * and calls ReplaceExonClauseList to consolidate the exons into a list
 * clause.
 */
static void RenameExonSequences (
  ValNodePtr PNTR list,
  BioseqPtr       bsp,
  Boolean         delete_now
)
{
  ValNodePtr       clause;
  Int4             num_exons;
  FeatureClausePtr fcp;
 
  if (list == NULL) return; 
  clause = *list;
  while (clause != NULL)
  {
    if (clause->choice == DEFLINE_CLAUSEPLUS
      && clause->data.ptrvalue != NULL)
    {
      fcp = clause->data.ptrvalue;
      if ( ! fcp->delete_me)
      {
        num_exons = GetNumberOfConsecutiveExons (clause, bsp);
        if (num_exons > 2)
        {
          ReplaceExonClauseList (fcp, clause, num_exons, bsp);
        }
        else
        {
          RenameExonSequences (&fcp->featlist, bsp, FALSE);
        }
      }
    }
    clause = clause->next;
  }
  if (delete_now) DeleteFeatureClauses (list);
}

static CharPtr organelleByGenome [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrial",
  "plastid",
  "",
  "",
  "",
  "",
  "",
  "cyanelle",
  "",
  "",
  "",
  "apicoplast",
  "leucoplast",
  "proplastid",
  "",
  "hydrogenosome",
  NULL,
};

static CharPtr organelleByPopup [] = {
  NULL,
  "mitochondrial",
  "chloroplast",
  "kinetoplast",
  "plastid",
  "chromoplast",
  "cyanelle",
  "apicoplast",
  "leucoplast",
  "proplastid",
  NULL
};

static void 
AddProductEnding 
(CharPtr    str, 
 BioseqPtr  bsp,
 Int2       mitochloroflag,
 ValNodePtr strings)
{
  Char orgnelle [80];
  BioSourcePtr  biop;
  ValNodePtr last_string;
  Int4 num_genes;
  SubSourcePtr  ssp;

  num_genes = 0;
  biop = GetBiopForBsp (bsp);

  if (biop != NULL) {
    if (FindStringInStrings (strings, "genes"))
    {
      num_genes = 2;
    }
    else if ((last_string = FindStringInStrings (strings, "gene")) != NULL
      && last_string->next != NULL
      && (last_string = FindStringInStrings (last_string->next, "gene")) != NULL)
    {
      num_genes = 2;
    }
    else
    {
      num_genes = 1;
    }

    orgnelle [0] = '\0';
  
    switch (biop->genome) {
    case GENOME_macronuclear :
      StringCat (str, "; macronuclear");
      break;
    case GENOME_nucleomorph :
      StringCat (str, "; nucleomorph");
      break;
    case GENOME_apicoplast :
    case GENOME_chloroplast :
    case GENOME_chromoplast :
    case GENOME_kinetoplast :
    case GENOME_mitochondrion :
    case GENOME_plastid :
    case GENOME_cyanelle :
    case GENOME_leucoplast :
    case GENOME_proplastid :
    case GENOME_hydrogenosome :
      sprintf (orgnelle, "; %s", organelleByGenome [biop->genome]);
      StringCat (str, orgnelle);
      break;
    default :
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != 255)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL
        && ssp->name != NULL
        && StringStr (ssp->name, "micronuclear"))
      {
        StringCat (str, "; micronuclear");
      }
      else if (mitochloroflag > 0) {
        if (mitochloroflag > 9) {
          /* beyond list */
        }
        else {
          if (num_genes > 1)
          {
            sprintf (orgnelle, "; nuclear genes for %s products",
                     organelleByPopup [mitochloroflag]);
          }
          else 
          {
            sprintf (orgnelle, "; nuclear gene for %s product",
                     organelleByPopup [mitochloroflag]);
          }
          StringCat (str, orgnelle);
        }
      }
      break;
    }
  }  
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* AutoDef_AddEnding () -- Add an ending on to the definition line     */
/*                         after the last feature.                     */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void AutoDef_AddEnding (
  ValNodePtr   clause_list,
  ValNodePtr PNTR strings,
  BioseqPtr    bsp,
  Int2         mitochloroflag,
  Boolean      alternate_splice_flag
)
{
  Char str [200];
  ValNodePtr last_string;
  Int4 new_data_len;
  CharPtr new_data;

  str[0] = 0;
  AddProductEnding (str, bsp, mitochloroflag, *strings);
  if (alternate_splice_flag) {
    StringCat (str, ", alternatively spliced");
  }

  StringCat (str, ".");

  last_string = *strings;
  if (last_string == NULL)
  {
    ValNodeAddStr (strings, 0, StringSave ( str));
  }
  else
  {
    while (last_string->next != NULL) last_string = last_string->next;
    new_data_len = StringLen (last_string->data.ptrvalue) + StringLen (str) + 1;
    new_data = (CharPtr) MemNew (new_data_len);
    if (new_data == NULL) return;
    StringCpy (new_data, last_string->data.ptrvalue);
    StringCat (new_data, str);
    MemFree (last_string->data.ptrvalue);
    last_string->data.ptrvalue = new_data;
  }
}

static Boolean LastIntervalChangeBeforeEnd (
  FeatureClausePtr onebefore,
  FeatureClausePtr thisclause,
  ValNodePtr rest_of_list
)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp;

  if (onebefore == NULL || rest_of_list == NULL) return FALSE;
 
  if (StringCmp (onebefore->interval, thisclause->interval) == 0) return FALSE;
  
  for (vnp = rest_of_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (StringCmp (thisclause->interval, fcp->interval) != 0) return FALSE;
    }
  }
  return TRUE;
    
}

static void PluralizeClauseIntervals (
  FeatureClausePtr fcp
)
{
  CharPtr new_interval, cp;

  if (fcp->interval != NULL
    && (cp = StringStr (fcp->interval, "gene, ")) != NULL)
  {
    new_interval = MemNew (StringLen (fcp->interval) + 2);
    if (new_interval == NULL) return;
    StringCpy (new_interval, fcp->interval);
    new_interval [ cp - fcp->interval + 4] = 's';
    StringCpy (new_interval + (cp - fcp->interval) + 5,
               cp + 4);
    MemFree (fcp->interval);
    fcp->interval = new_interval;
  }
}  

static Boolean DisplayAlleleName (FeatureClausePtr thisclause)
{
  if (thisclause == NULL) return FALSE;
  if (StringCmp (thisclause->feature_label_data.typeword, "gene") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "mRNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene mRNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "precursor RNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene precursor RNA") == 0)
  {
    return TRUE;
  }
  return FALSE;
}

static void ListClauses (
  ValNodePtr clauselist,
  ValNodePtr PNTR strings,
  Boolean    allow_semicolons,
  Boolean    suppress_final_and
)
{
  FeatureClausePtr thisclause, onebefore, twobefore, oneafter, twoafter;
  Boolean print_typeword;
  Boolean print_and;
  Boolean print_comma;
  Boolean print_semicolon;
  Boolean print_comma_between_description_and_typeword;
  Boolean typeword_is_plural;
  size_t clause_len;
  CharPtr clause_string;
  Boolean oneafter_has_detail_change;
  Boolean oneafter_has_interval_change;
  Boolean oneafter_has_typeword_change;
  Boolean onebefore_has_detail_change;
  Boolean onebefore_has_interval_change;
  Boolean onebefore_has_typeword_change;
  SeqFeatPtr main_feat;
  CharPtr new_interval;
  ValNodePtr voneafter, vtwoafter;

  while (clauselist != NULL && clauselist->choice != DEFLINE_CLAUSEPLUS)
  {
    clauselist = clauselist->next;
  }
  if (clauselist == NULL) return;

  thisclause = clauselist->data.ptrvalue;
  onebefore = NULL;
  twobefore = NULL;
  
  while (thisclause != NULL)
  {
    oneafter_has_detail_change = FALSE;
    oneafter_has_interval_change = FALSE;
    oneafter_has_typeword_change = FALSE;
    onebefore_has_detail_change = FALSE;
    onebefore_has_interval_change = FALSE;
    onebefore_has_typeword_change = FALSE;
    if (onebefore != NULL)
    {
      if (StringCmp (onebefore->interval, thisclause->interval) != 0)
        onebefore_has_interval_change = TRUE;
      if (StringCmp (onebefore->feature_label_data.typeword,
                     thisclause->feature_label_data.typeword) != 0)
      {
        onebefore_has_typeword_change = TRUE;
      }
      if (onebefore_has_typeword_change || onebefore_has_interval_change
          || (DisplayAlleleName (onebefore) && StringLen (onebefore->allelename) != 0)
          || (DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) != 0))
     {
        onebefore_has_detail_change = TRUE;  
      }
    }
    voneafter = clauselist->next;
    while (voneafter != NULL && voneafter->choice != DEFLINE_CLAUSEPLUS)
    {
      voneafter = voneafter->next;
    }
    if (voneafter == NULL)
    {
      vtwoafter = NULL;
    }
    else
    {
      vtwoafter = voneafter->next;
      while (vtwoafter != NULL && vtwoafter->choice != DEFLINE_CLAUSEPLUS)
      {
        vtwoafter = vtwoafter->next;
      }
    }

    if (voneafter != NULL)
    {
      oneafter = voneafter->data.ptrvalue;
      if (StringCmp (oneafter->interval, thisclause->interval) != 0)
        oneafter_has_interval_change = TRUE;
      if (StringCmp (oneafter->feature_label_data.typeword,
                     thisclause->feature_label_data.typeword) != 0)
      {
        oneafter_has_typeword_change = TRUE;
      }
      if (oneafter_has_typeword_change  || oneafter_has_interval_change
          || (DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) != 0)
          || (DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) != 0))
      {
        oneafter_has_detail_change = TRUE;
      }
      if (vtwoafter != NULL)
      {
        twoafter = vtwoafter->data.ptrvalue;
      }
      else
      {
        twoafter = NULL;
      }
    }
    else
    {
      oneafter = NULL;
      twoafter = NULL;
    }
    print_typeword = FALSE;
    typeword_is_plural = FALSE;
    print_and = FALSE;
    print_comma = FALSE;
    print_semicolon = FALSE;

    if (thisclause->feature_label_data.is_typeword_first)
    {
      if (onebefore == NULL || onebefore_has_detail_change)
      {
        print_typeword = TRUE;
        if (oneafter != NULL && ! oneafter_has_detail_change)
        {
          typeword_is_plural = TRUE;
        }
        else if (StringStr (thisclause->feature_label_data.description, " through ") != NULL
          && StringCmp (thisclause->feature_label_data.typeword, "exon") == 0)
        {
          typeword_is_plural = TRUE;
        }
      }
    }
    else
    {
      if (oneafter == NULL || oneafter_has_detail_change)
      {
        print_typeword = TRUE;
        if (onebefore != NULL && ! onebefore_has_detail_change)
        {
          typeword_is_plural = TRUE;
        }
      }
    }

    /* when to print and before this section */
    if ( onebefore != NULL
         && ! onebefore_has_detail_change
         && (oneafter == NULL || oneafter_has_detail_change))
    {
      print_and = TRUE;
    }
    else if (oneafter == NULL && onebefore != NULL)
    {
      print_and = TRUE;
    }
    else if (onebefore != NULL
         && ! onebefore_has_interval_change
         && oneafter_has_interval_change)
    {
      print_and = TRUE;
    }
    else if ( LastIntervalChangeBeforeEnd ( onebefore, 
                                            thisclause,
                                            clauselist->next))
    {
      print_and = TRUE;
    }

    if (suppress_final_and && oneafter == NULL)
    {
      print_and = FALSE;
    }
    if (suppress_final_and && oneafter != NULL && twoafter == NULL)
    {
      print_comma = TRUE;
    }
    
    /* when to print semicolon after this section */
    /* after every interval change except when exons change "interval" */
    /* exons changing interval are going from alt-spliced to not */
    /* or vice versa, in either case we don't want a semicolon or comma */
    if (oneafter != NULL && oneafter_has_interval_change
      && (StringCmp (thisclause->feature_label_data.typeword, "exon") != 0
         || StringCmp (oneafter->feature_label_data.typeword, "exon") != 0))
    {
      print_semicolon = TRUE;
    }

    /* when to print comma after this section */
    if (onebefore != NULL && oneafter != NULL
      && ! onebefore_has_detail_change
      && ! oneafter_has_detail_change )
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && onebefore != NULL
      && ! onebefore_has_interval_change && ! oneafter_has_interval_change
      &&  onebefore_has_typeword_change &&  oneafter_has_typeword_change)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && ! oneafter_has_detail_change
      && StringCmp (twoafter->feature_label_data.typeword,
                    thisclause->feature_label_data.typeword) == 0
      && StringCmp (twoafter->interval,
                    thisclause->interval) == 0)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL  && twoafter != NULL
      && oneafter_has_typeword_change
      && StringCmp (twoafter->feature_label_data.typeword,
                    oneafter->feature_label_data.typeword) == 0
      && StringCmp (twoafter->interval,
                    oneafter->interval) == 0
      && ! print_and)
    {
      print_comma = TRUE;
    }
    else if (((oneafter_has_interval_change || oneafter == NULL)
      && StringDoesHaveText (thisclause->interval))
      || (oneafter_has_interval_change && oneafter != NULL && ! print_semicolon))
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && !oneafter_has_interval_change
      && StringCmp (thisclause->interval, twoafter->interval) == 0
      && oneafter_has_typeword_change
      && StringCmp (thisclause->feature_label_data.typeword,
                    twoafter->feature_label_data.typeword) != 0)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && onebefore != NULL && twoafter != NULL
      && ! oneafter_has_interval_change && ! onebefore_has_interval_change
      && StringCmp (thisclause->interval, twoafter->interval) == 0
      && oneafter_has_typeword_change)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && oneafter_has_typeword_change
      && StringCmp (oneafter->feature_label_data.typeword,
                    twoafter->feature_label_data.typeword) != 0
      && ! oneafter_has_interval_change
      && StringCmp (oneafter->interval, twoafter->interval) == 0)
    {
      /* spacer 1, foo RNA gene, and spacer2, complete sequence */
      /*         ^ */
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL 
      && ! oneafter_has_interval_change && StringCmp (thisclause->interval, twoafter->interval) == 0
      && ((DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) > 0)
        || (DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) > 0)))
    {
      print_comma = TRUE;  	
    }
    else if (oneafter != NULL && onebefore != NULL
      && ! oneafter_has_interval_change && ! onebefore_has_interval_change
      && ((DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) > 0)
        || (DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) > 0)))
    {
      print_comma = TRUE;  	
    }

    if (thisclause->featlist != NULL
      && thisclause->featlist->data.ptrvalue != NULL
      && StringDoesHaveText (thisclause->interval)
      && StringNCmp (thisclause->interval, "partial", 7) != 0
      && StringNCmp (thisclause->interval, "complete", 8) != 0)
    {
      main_feat = thisclause->featlist->data.ptrvalue;
      if (IsMobileElement (main_feat)
        || IsEndogenousVirusSourceFeature (main_feat) )
      {
        print_comma = FALSE;
      }
    }

    if (onebefore != NULL
      && ! onebefore_has_interval_change
      && (oneafter_has_interval_change || oneafter == NULL))
    {
      PluralizeClauseIntervals (thisclause);
    }

    if ( thisclause->make_plural )
    {
      if ((onebefore != NULL && ! onebefore_has_detail_change)
        || (oneafter != NULL && !oneafter_has_detail_change))
      {
        PluralizeConsolidatedClauseDescription (thisclause);
      }
      else
      {
        typeword_is_plural = TRUE;
      }
    }

    clause_len = StringLen (thisclause->feature_label_data.description) + 1;
    
    /* add one in case we need to add the semicolon to this clause (when
     * the interval has changed because this clause has no interval and
     * the next one does).
     */
    clause_len++;

    /* we need to place a comma between the description and the type word 
     * when the description ends with "precursor" or when the type word
     * starts with "precursor"
     */
    if ( thisclause->feature_label_data.description != NULL
      && ! thisclause->feature_label_data.is_typeword_first
      && print_typeword
      && ! StringHasNoText (thisclause->feature_label_data.typeword)
      && ((StringNCmp (thisclause->feature_label_data.typeword, "precursor", 9) == 0
            && thisclause->feature_label_data.description [StringLen (thisclause->feature_label_data.description) - 1] != ')')
          || (clause_len > StringLen ("precursor")
              && StringCmp ( thisclause->feature_label_data.description
                     + clause_len - StringLen ("precursor") - 2,
                     "precursor") == 0)))
    {
      print_comma_between_description_and_typeword = TRUE;
      clause_len += 1;
    }
    else
    {
      print_comma_between_description_and_typeword = FALSE;
    }

    if (print_typeword)
      clause_len += StringLen (thisclause->feature_label_data.typeword) + 1;
    if (typeword_is_plural)
      clause_len += 1;
    if (print_and)
      clause_len += 4;
    if (print_comma)
      clause_len += 2;
    if (DisplayAlleleName (thisclause))
    {
      clause_len += StringLen (thisclause->allelename) + 10;
      if (StringLen (thisclause->allelename) > 0) 
      {
        clause_len += StringLen (thisclause->allelename) + StringLen ("allele ");
      }
    }
    
    clause_string = (CharPtr) MemNew (clause_len);
    if (clause_string == NULL)
      return;
    clause_string[0] = 0;
    if (print_and)
      StringCat (clause_string, "and ");
    if (thisclause->feature_label_data.is_typeword_first && print_typeword
      && thisclause->feature_label_data.typeword != NULL
      && ! StringHasNoText (thisclause->feature_label_data.typeword))
    {
      StringCat (clause_string, thisclause->feature_label_data.typeword);
      if (typeword_is_plural)
        StringCat (clause_string, "s");
      if (thisclause->feature_label_data.description != NULL)
        StringCat (clause_string, " ");
    }
    if (thisclause->feature_label_data.description != NULL)
    {
      StringCat (clause_string, thisclause->feature_label_data.description);
      if (print_comma_between_description_and_typeword)
      {
        StringCat (clause_string, ",");
      }
    }
    if (! thisclause->feature_label_data.is_typeword_first && print_typeword
      && thisclause->feature_label_data.typeword != NULL
      && ! StringHasNoText (thisclause->feature_label_data.typeword))
    {
      if (!StringHasNoText (thisclause->feature_label_data.description))
        StringCat (clause_string, " ");
      StringCat (clause_string, thisclause->feature_label_data.typeword);
      if (typeword_is_plural)
        StringCat (clause_string, "s");
      if (DisplayAlleleName (thisclause)
        && thisclause->allelename != NULL)
      {
        StringCat (clause_string, ", ");
        StringCat (clause_string, thisclause->allelename);
        StringCat (clause_string, " allele");
      }
    }
    if (StringLen (clause_string) > 0 ) 
    {
      if (print_comma)
        StringCat (clause_string, ",");
      ValNodeAddStr (strings, 0, clause_string);
    }
    else 
    {
    	MemFree (clause_string);
    	clause_string = NULL;
    }
 
    if (oneafter == NULL || oneafter_has_interval_change)
    {
      if (print_semicolon) {
        if (thisclause->interval == NULL
          || StringHasNoText(thisclause->interval)) {
          if (clause_string != NULL) {
            StringCat (clause_string, ";");
          }
        } else if (thisclause->interval[StringLen (thisclause->interval)] != ';') {
          new_interval = MemNew (StringLen (thisclause->interval) + 2);
          if (new_interval == NULL) return;
          StringCpy (new_interval, thisclause->interval);
          if (allow_semicolons) 
          {
            StringCat (new_interval, ";");
          }
          else
          {
            StringCat (new_interval, ",");
          }
          MemFree (thisclause->interval);
          thisclause->interval = new_interval;
        }
      }
      if (thisclause->interval != NULL
        && !StringHasNoText (thisclause->interval))
      {
        ValNodeAddStr (strings, 0, StringSave (thisclause->interval));
      }
    }
    twobefore = onebefore;
    onebefore = thisclause;   
    thisclause = oneafter;
    clauselist = voneafter;
  }
}

static Uint1 GetMoleculeType 
(BioseqPtr bsp,
 Uint2     entityID)
{
  SeqDescPtr         sdp;
  MolInfoPtr         mip;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL) return MOLECULE_TYPE_GENOMIC;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return MOLECULE_TYPE_GENOMIC;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return MOLECULE_TYPE_GENOMIC;
  return mip->biomol;
}

static Boolean SpecialHandlingForSpecialTechniques (
  BioseqPtr bsp
)
{
  SeqDescPtr sdp;
  MolInfoPtr mip;
  ValNodePtr vnp;

  if (bsp == NULL) return MOLECULE_TYPE_GENOMIC;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, NULL);
  if (sdp == NULL)
  {
    for (sdp = bsp->descr;
         sdp != NULL && sdp->choice != Seq_descr_molinfo;
         sdp = sdp->next)
    {}
  }
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return FALSE;
  if (mip->tech == MI_TECH_htgs_0 ||
      mip->tech == MI_TECH_htgs_1 ||
      mip->tech == MI_TECH_htgs_2 ||
      mip->tech == MI_TECH_est ||
      mip->tech == MI_TECH_sts ||
      mip->tech == MI_TECH_survey ||
      mip->tech == MI_TECH_wgs) {
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
    if (vnp != NULL)
      vnp = ValNodeFreeData (vnp);
    return TRUE;
  }

  return FALSE;
}

static Boolean LIBCALLBACK ShouldRemoveExon (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  Boolean partial3, partial5;
  SeqFeatPtr main_feat;

  if (isSegment || isLonely || isRequested) return FALSE;
  if (parent_fcp == NULL
    || parent_fcp->featlist == NULL
    || parent_fcp->featlist->data.ptrvalue == NULL)
  {
    return TRUE;
  }
  
  main_feat = parent_fcp->featlist->data.ptrvalue;
  if ( IsCDS (main_feat))
  {
    CheckSeqLocForPartial (main_feat->location, &partial5, &partial3);
    if (partial5 || partial3) return FALSE;
  }
  else if (IsmRNA (main_feat) || parent_fcp->has_mrna)
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK ShouldRemoveCDS (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr description;
  Boolean retval = FALSE;

  description = GetGeneProtDescription (this_fcp, bsp, rp);
  if (StringHasNoText (description))
  {
    retval = TRUE;
  }
  if (description != NULL) MemFree (description);
  return retval;
}

static Boolean LIBCALLBACK ShouldRemoveNoncodingProductFeat (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK ShouldRemovePromoter (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  /* remove a promoter if it is in an mRNA or gene clause */
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL 
      && (parent_fcp->has_mrna 
        || (parent_fcp->featlist != NULL
           && parent_fcp->featlist->choice == DEFLINE_FEATLIST
           && parent_fcp->featlist->data.ptrvalue != NULL
           && IsmRNA (parent_fcp->featlist->data.ptrvalue))))
  {
    return TRUE;
  }
  else if (isLonely)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean LIBCALLBACK ShouldRemoveLTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL)
  {
    return TRUE;
  }
  else if (isLonely)
    return FALSE;
  else
    return TRUE;
}

static Boolean LIBCALLBACK ShouldRemove3UTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{ 
  if (isLonely || isRequested)
    return FALSE;
  else
    return TRUE;
}

static Boolean LIBCALLBACK ShouldRemove5UTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isLonely || isRequested)
    return FALSE;
  else
    return TRUE;
}

static Boolean LIBCALLBACK ShouldRemoveIntron (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL 
      && (parent_fcp->has_mrna 
        || (parent_fcp->featlist != NULL
           && parent_fcp->featlist->choice == DEFLINE_FEATLIST
           && parent_fcp->featlist->data.ptrvalue != NULL
           && IsmRNA (parent_fcp->featlist->data.ptrvalue))))
  {
    return TRUE;
  }
  else if (isLonely)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean LIBCALLBACK ShouldRemoveMobileElement
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  return (!isLonely && !isRequested);
}

static Boolean LIBCALLBACK ShouldRemoveGeneric 
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr productname;
  Boolean rval;

  rval = FALSE;
  if (IsMiscRNA (sfp) && ( productname = GetProductName (sfp, bsp, rp)) != NULL)
  {
    if (StringStr (productname, "trans-spliced leader") != NULL)
    {
      rval = TRUE;
    }
    MemFree (productname);
  }
    
  return rval;
}


static Boolean IsBioseqPrecursorRNA (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mol;
  
  if (bsp == NULL) return FALSE;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp != NULL && sdp->data.ptrvalue != NULL)
  {
		mol = (MolInfoPtr) sdp->data.ptrvalue;
    if (mol->biomol == 2)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK ShouldRemovePrecursorRNA
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  if (!isLonely && IsBioseqPrecursorRNA(bsp))
  {
    return TRUE;
  }
  else
  {
    return ShouldRemoveGeneric (sfp, parent_fcp, this_fcp, bsp, isLonely, 
                                isRequested, isSegment, rp);
  }
}


typedef struct removableitemglobal {
  matchFunction  itemmatch;
  ShouldRemoveFunction ShouldRemove;
  CharPtr  group_name;
} RemovableItemGlobalData, PNTR RemovableItemGlobalPtr;

typedef struct removableitemlocal {
/*  ButtoN  keep_request; */
  Boolean  keep;
} RemovableItemLocalData, PNTR RemovableItemLocalPtr;

static RemovableItemGlobalData remove_items[] = {
  { IsExon, ShouldRemoveExon, "Exons" },
  { IsIntron, ShouldRemoveIntron, "Introns" },
  { Is5UTR, ShouldRemove5UTR, "5' UTRs" },
  { Is3UTR, ShouldRemove3UTR, "3' UTRs" },
  { IsCDS,  ShouldRemoveCDS, "CDSs" },
  { IsPromoter, ShouldRemovePromoter, "Promoters:" },
  { IsLTR, ShouldRemoveLTR, "LTRs" },
  { IsNoncodingProductFeat,  ShouldRemoveNoncodingProductFeat, "Misc feats with comments:" },
  { IsRemovableMobileElement, ShouldRemoveMobileElement, "Optional Mobile Element" },
  { IsPrecursorRNA, ShouldRemovePrecursorRNA, "Precursor RNAs" }
};


NLM_EXTERN CharPtr GetRemovableItemName (Int4 i)
{
  if (i < 0 || i >= NumRemovableItems) {
    return NULL;
  } else {
    return remove_items[i].group_name;
  }
}

NLM_EXTERN void InitFeatureRequests (
  DeflineFeatureRequestListPtr feature_requests
)
{
  Int4 i;
  for (i=0; i < NumRemovableItems; i++)
  {
    feature_requests->keep_items[i] = FALSE;
  }
  feature_requests->add_fake_promoters = TRUE;
  feature_requests->suppress_alt_splice_phrase = FALSE;
  feature_requests->remove_subfeatures = FALSE;
  feature_requests->feature_list_type = DEFLINE_USE_FEATURES;
  feature_requests->misc_feat_parse_rule = 2;
  feature_requests->suppress_locus_tags = FALSE;
  feature_requests->suppressed_feature_list = NULL;
  feature_requests->use_ncrna_note = FALSE;
}

static Boolean RemoveCondition (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isSegment,
  DeflineFeatureRequestList *feature_requests
)
{
  Int4 i;
  if (sfp == NULL) return TRUE;
  for (i=0; i < NumRemovableItems; i++)
  {
    if (remove_items[i].itemmatch (sfp))
      return remove_items[i].ShouldRemove (sfp, parent_fcp, this_fcp, bsp,
                                           isLonely, feature_requests->keep_items[i],
                                           isSegment,
                                           feature_requests);
  }
  return ShouldRemoveGeneric(sfp, parent_fcp, this_fcp, bsp, isLonely, FALSE,
                             isSegment, feature_requests);
}

static Boolean FindOtherGeneClause 
( ValNodePtr feature_list,
  ValNodePtr me,
  GeneRefPtr grp,
  Boolean    suppress_locus_tag)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  if (grp == NULL) return FALSE;

  for (vnp = feature_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp == me) continue;
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp->delete_me) continue;
      if ( fcp->grp == grp
        || (fcp->grp != NULL && DoGenesMatch (fcp->grp, grp, suppress_locus_tag)))
      {
        return TRUE;
      }
      if ( FindOtherGeneClause (fcp->featlist, me, grp, suppress_locus_tag))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}
 
static void RemoveGenesMentionedElsewhere 
( ValNodePtr PNTR feature_list,
  ValNodePtr      search_list,
  Boolean         delete_now,
  Boolean         suppress_locus_tag)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  for (vnp = *feature_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp->featlist == NULL)
      {
        continue;
      }
      if ( IsGene (fcp->featlist->data.ptrvalue)
        && fcp->featlist->next == NULL
        && FindOtherGeneClause ( search_list, vnp, fcp->grp, suppress_locus_tag))
      {
        fcp->delete_me = TRUE;
      }
      else
      {
        RemoveGenesMentionedElsewhere ( &(fcp->featlist), search_list, FALSE, suppress_locus_tag);
      }
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void MarkUnwantedFeatureClauseForRemoval (
  ValNodePtr clause,
  BioseqPtr  bsp,
  Boolean    isLonely,
  FeatureClausePtr parent_fcp,
  Boolean    isSegment,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  FeatureClausePtr fcp;
  ValNodePtr       featlist;
  ValNodePtr       firstfeat;
  Int4             clause_count;
  SeqFeatPtr       sfp;

  if (clause == NULL
    || clause->choice != DEFLINE_CLAUSEPLUS
    || clause->data.ptrvalue == NULL)
  {
    return;
  }

  fcp = clause->data.ptrvalue;
  firstfeat = fcp->featlist;
  clause_count = 0;
  for (featlist = firstfeat;
       featlist != NULL && isLonely;
       featlist = featlist->next)
  {
    if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      clause_count ++;
      if (clause_count > 1)
      {
        isLonely = FALSE;
      }
    }
  }
    
  featlist = firstfeat;
  while (featlist != NULL)
  {  
    if (featlist->choice == DEFLINE_FEATLIST
      && featlist->data.ptrvalue != NULL)
    {
      sfp = (SeqFeatPtr) featlist->data.ptrvalue;
      if (RemoveCondition (featlist->data.ptrvalue, parent_fcp, fcp, bsp, 
                          isLonely, isSegment, feature_requests))
      {
        fcp->delete_me = TRUE;
      }
      else if (! IsGene (sfp) && ! IsmRNA (sfp))
      {
        isLonely = FALSE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS 
      && featlist->data.ptrvalue != NULL)
    {
      MarkUnwantedFeatureClauseForRemoval (featlist, bsp, isLonely, fcp,
                                           isSegment,
                                           feature_requests);
    }
    featlist = featlist->next;
  }
}
  
static void RemoveUnwantedFeatures (
  ValNodePtr PNTR list,
  BioseqPtr bsp,
  Boolean   isSegment,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  ValNodePtr       vnp;
  Boolean          isLonely;

  if (list == NULL) return;

  isLonely = TRUE;

  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->next != NULL) isLonely = FALSE;
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      MarkUnwantedFeatureClauseForRemoval (vnp, bsp, isLonely, NULL,
                                           isSegment, feature_requests);
    }
  }
  DeleteFeatureClauses (list);
}

static Boolean IsFeatureInSelectionList (SeqFeatPtr sfp, ValNodePtr feat_list)
{
  ValNodePtr       vnp;

  if (sfp == NULL || feat_list == NULL)
  {
    return FALSE;
  }
  
  for (vnp = feat_list; vnp != NULL && sfp->idx.subtype != vnp->choice; vnp = vnp->next)
  {
  }
  if (vnp == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static void MarkSuppressedFeatureClauseForRemoval (
  ValNodePtr clause,
  ValNodePtr suppressed_feature_list
)
{
  FeatureClausePtr fcp;
  ValNodePtr       featlist;
  ValNodePtr       firstfeat;
  SeqFeatPtr       sfp;

  if (clause == NULL
    || clause->choice != DEFLINE_CLAUSEPLUS
    || clause->data.ptrvalue == NULL)
  {
    return;
  }

  fcp = clause->data.ptrvalue;
  firstfeat = fcp->featlist;
    
  featlist = firstfeat;
  while (featlist != NULL)
  {  
    if (featlist->choice == DEFLINE_FEATLIST
      && featlist->data.ptrvalue != NULL)
    {
      sfp = (SeqFeatPtr) featlist->data.ptrvalue;
      if (IsFeatureInSelectionList (sfp, suppressed_feature_list))
      {
        fcp->delete_me = TRUE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS 
      && featlist->data.ptrvalue != NULL)
    {
      MarkSuppressedFeatureClauseForRemoval (featlist, suppressed_feature_list);
    }
    featlist = featlist->next;
  }
}
  
static void RemoveSuppressedFeatures (ValNodePtr PNTR list,
                                      ValNodePtr suppressed_feature_list)
{
  ValNodePtr vnp;
  
  if (list == NULL || *list == NULL || suppressed_feature_list == NULL)
  {
    return;
  }
  
  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      MarkSuppressedFeatureClauseForRemoval (vnp, suppressed_feature_list);
    }
  }
  DeleteFeatureClauses (list);  
}

static Boolean LIBCALLBACK IsMasterClause (
  SeqFeatPtr sfp
)
{
  if ( IsMobileElement (sfp)) return TRUE;
  return FALSE;
}

static void DeleteSubfeatures (
  ValNodePtr PNTR feature_list,
  Boolean         delete_now
)
{
  ValNodePtr       clause, featlist;
  FeatureClausePtr clause_fcp, fcp;

  if (feature_list == NULL) return;
  for (clause = *feature_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (clause_fcp = clause->data.ptrvalue) == NULL
      || clause_fcp->featlist == NULL)
    {
      continue;
    }
    if (clause_fcp->featlist->choice == DEFLINE_FEATLIST
      && IsMasterClause (clause_fcp->featlist->data.ptrvalue))
    {
      for (featlist = clause_fcp->featlist->next;
           featlist != NULL;
           featlist = featlist->next)
      {
        if (featlist->choice == DEFLINE_CLAUSEPLUS
          && (fcp = featlist->data.ptrvalue) != NULL)
        {
          fcp->delete_me = TRUE;
        }
      }
    }
    else
    {
      DeleteSubfeatures ( &(clause_fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void DeleteOperonAndGeneClusterSubfeatures (
  ValNodePtr PNTR feature_list,
  Boolean         delete_now
)
{
  ValNodePtr       clause, featlist;
  FeatureClausePtr clause_fcp, fcp;

  if (feature_list == NULL) return;
  for (clause = *feature_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (clause_fcp = clause->data.ptrvalue) == NULL
      || clause_fcp->featlist == NULL)
    {
      continue;
    }
    if (clause_fcp->featlist->choice == DEFLINE_FEATLIST
      && (IsOperon (clause_fcp->featlist->data.ptrvalue) 
          || IsGeneCluster (clause_fcp->featlist->data.ptrvalue)))
    {
      for (featlist = clause_fcp->featlist->next;
           featlist != NULL;
           featlist = featlist->next)
      {
        if (featlist->choice == DEFLINE_CLAUSEPLUS
          && (fcp = featlist->data.ptrvalue) != NULL)
        {
          fcp->delete_me = TRUE;
        }
      }
    }
    else
    {
      DeleteOperonAndGeneClusterSubfeatures ( &(clause_fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void RemoveFeats (
  ValNodePtr    list,
  matchFunction itemmatch
)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;
  
  if (list == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST
      && itemmatch (vnp->data.ptrvalue))
    {
      vnp->choice = DEFLINE_REMOVEFEAT;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      RemoveFeats (fcp->featlist, itemmatch);
    }
  }
}

/* A clause is "tall" if it has only one clause at any level */
static Boolean IsClauseTall (
  FeatureClausePtr fcp
)
{
  ValNodePtr featlist;
  Int4       num_clauses;
  FeatureClausePtr subclause;

  num_clauses = 0;
  if (fcp == NULL) return FALSE;
  subclause = NULL;
  if (fcp->featlist == NULL) return FALSE;
  for (featlist = fcp->featlist;
       featlist != NULL;
       featlist = featlist->next)
  {
    if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      subclause = featlist->data.ptrvalue;
      if (subclause == NULL || ! IsClauseTall (subclause))
      {
        return FALSE;
      }
      num_clauses ++;
      if (num_clauses > 1) return FALSE;
    }
  }
  if (subclause == NULL || ! subclause->feature_label_data.is_typeword_first)
  {
    return TRUE;
  }
  return FALSE;
}

static void SmashOneTallClause (
  FeatureClausePtr fcp
)
{
  FeatureClausePtr subclause;
  ValNodePtr       featlist;
  ValNodePtr       subclause_featlist;
  ValNodePtr       subclause_firstclause;
  CharPtr          new_description;
  Int4             new_description_len;
  SeqFeatPtr       main_feat;

  if (fcp == NULL || fcp->featlist == NULL) return;

  /* move features up */
  featlist = fcp->featlist;
  if (featlist->choice == DEFLINE_FEATLIST)
  {
    main_feat = fcp->featlist->data.ptrvalue;
  }
  else
  {
    main_feat = NULL;
  }
  
  while (featlist != NULL && featlist->choice != DEFLINE_CLAUSEPLUS)
  {
    featlist = featlist->next;
  }
  if (featlist == NULL) return;
  subclause = featlist->data.ptrvalue;
  if (subclause == NULL) return;
 
  /* move subclause feats to top of list */
  if (subclause->featlist != NULL
    && subclause->featlist->choice == DEFLINE_FEATLIST)
  {
    subclause_featlist = subclause->featlist;
    while (subclause->featlist != NULL
           && subclause->featlist->next != NULL
           && subclause->featlist->next->choice == DEFLINE_FEATLIST)
    {
      subclause->featlist = subclause->featlist->next;
    }
    if (subclause->featlist != NULL)
    {
      subclause_firstclause = subclause->featlist->next;
      subclause->featlist->next = fcp->featlist;
      fcp->featlist = subclause->featlist;
      subclause->featlist = subclause_firstclause;
    }
  }

  /* create new description */
  new_description_len = StringLen (subclause->feature_label_data.description)
                   + StringLen (fcp->feature_label_data.description)
                   + StringLen (fcp->feature_label_data.typeword)
                   + 4;
  new_description = (CharPtr) MemNew (new_description_len);
  if (new_description == NULL) return;
  new_description [0] = 0;
  if ( fcp->feature_label_data.is_typeword_first)
  {
    StringCat (new_description, fcp->feature_label_data.typeword);
    StringCat (new_description, " ");
  }
  StringCat (new_description, fcp->feature_label_data.description);
  if ( ! fcp->feature_label_data.is_typeword_first)
  {
    StringCat (new_description, fcp->feature_label_data.typeword);
  }

  if ( ! IsMobileElement (main_feat)
    && ! IsEndogenousVirusSourceFeature (main_feat))
  {
    StringCat (new_description, ",");
  }
  StringCat (new_description, " ");
  StringCat (new_description, subclause->feature_label_data.description);

  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
  }
  fcp->feature_label_data.description = new_description;
 
  /* move interval up */
  if (fcp->interval != NULL)
  {
    MemFree (fcp->interval);
  }
  fcp->interval = subclause->interval;
  subclause->interval = NULL;

  /* move typeword up */
  fcp->feature_label_data.typeword = subclause->feature_label_data.typeword;
  fcp->feature_label_data.is_typeword_first = 
               subclause->feature_label_data.is_typeword_first;
  subclause->feature_label_data.typeword = NULL;
  subclause->delete_me = TRUE;

}


static void SmashTallClauses (
  ValNodePtr PNTR clause_list,
  Boolean         delete_now
)
{
  ValNodePtr clause;
  FeatureClausePtr fcp;

  if (clause_list == NULL) return;
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if ( IsClauseTall (fcp))
    {
      SmashOneTallClause (fcp);
    }
    else
    {
      SmashTallClauses (& (fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (clause_list);
  }
}

static ValNodePtr RemoveAllButLastCDS (
  ValNodePtr list,
  ValNodePtr last_cds
)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  /* now remove all CDSs except the last one */
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST
      && IsCDS (vnp->data.ptrvalue))
    {
      if (last_cds != NULL)
      {
        last_cds->choice = DEFLINE_REMOVEFEAT;
      }
      last_cds = vnp;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      last_cds = RemoveAllButLastCDS (fcp->featlist, last_cds);
    }
  }
  return last_cds;
}

static Boolean OkToConsolidate (
  CharPtr last_desc,
  CharPtr new_desc,
  Boolean last_partial,
  Boolean new_partial,
  FeatureClausePtr last_fcp,
  FeatureClausePtr fcp
)
{
  if (StringCmp (last_desc, new_desc) != 0) return FALSE;
  if (new_partial != last_partial) return FALSE;
  if ( ( fcp->is_alt_spliced && ! last_fcp->is_alt_spliced)
      || (! fcp->is_alt_spliced && last_fcp->is_alt_spliced))
  {
    return FALSE;
  }
  if (fcp->featlist == NULL || last_fcp->featlist == NULL) return FALSE;
  if ( fcp->featlist->choice != DEFLINE_FEATLIST) return FALSE;
  if ( last_fcp->featlist->choice != DEFLINE_FEATLIST) return FALSE;
  if ( (IsCDS (fcp->featlist->data.ptrvalue)
        && ! IsCDS (last_fcp->featlist->data.ptrvalue)
        && ! IsGene (last_fcp->featlist->data.ptrvalue))
      || (! IsCDS (fcp->featlist->data.ptrvalue)
        && ! IsGene (fcp->featlist->data.ptrvalue)
        && IsCDS (last_fcp->featlist->data.ptrvalue)))
  {
    return FALSE;
  }
  if ((IsExon (fcp->featlist->data.ptrvalue) && !IsExon(last_fcp->featlist->data.ptrvalue))
      || (IsExon (last_fcp->featlist->data.ptrvalue) && !IsExon(fcp->featlist->data.ptrvalue))
      || (IsIntron (fcp->featlist->data.ptrvalue) && !IsIntron(last_fcp->featlist->data.ptrvalue))
      || (IsIntron (last_fcp->featlist->data.ptrvalue) && !IsIntron(fcp->featlist->data.ptrvalue)))
  {
    return FALSE;
  }
  return TRUE;
}

static void RemoveRedundantGeneFeatureFromConsolidatedClause (
  FeatureClausePtr fcp
)
{
  ValNodePtr featlist, prevfeat, tmpfeat;
  SeqFeatPtr feat1, feat2;

  prevfeat = NULL;
  featlist = fcp->featlist;
  while ( featlist != NULL
         && featlist->choice == DEFLINE_FEATLIST
         && featlist->next != NULL
         && featlist->next->choice == DEFLINE_FEATLIST)
  {
    feat1 = featlist->data.ptrvalue;
    feat2 = featlist->next->data.ptrvalue;
    if (feat1 == NULL || feat2 == NULL) return;
    if (IsGene (feat1) && ! IsGene (feat2))
    {
      if (prevfeat == NULL)
      {
        fcp->featlist = featlist->next;
        featlist->next = NULL;
        FreeListElement (featlist);
        featlist = fcp->featlist->next;
      }
      else
      {
        prevfeat->next = featlist->next;
        featlist->next = NULL;
        FreeListElement (featlist);
        featlist = prevfeat->next;
      }
    }
    else if ( !IsGene (feat1) && IsGene (feat2))
    {
      tmpfeat = featlist->next;
      featlist->next = tmpfeat->next;
      tmpfeat->next = NULL;
      FreeListElement (tmpfeat);
    }
    else
    {
      featlist = featlist->next;
    }
  }
}

static void PluralizeConsolidatedClauseDescription (
  FeatureClausePtr fcp
)
{
  CharPtr new_desc;

  /* prevent crash */
  if (fcp == NULL || fcp->feature_label_data.description == NULL) return;

  /* don't pluralize tRNA names */
  if (StringNCmp (fcp->feature_label_data.description, "tRNA-", 5) ==0) return;

  new_desc = MemNew (StringLen (fcp->feature_label_data.description) + 2);
  if (new_desc == NULL) return;

  StringCpy (new_desc, fcp->feature_label_data.description);
  StringCat (new_desc, "s");
  MemFree (fcp->feature_label_data.description);
  fcp->feature_label_data.description = new_desc;
}

static void ConsolidateClauses (
  ValNodePtr PNTR list,
  BioseqPtr  bsp,
  Uint1      biomol,
  Boolean    delete_now,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp;
  FeatureClausePtr last_cds_fcp;
  CharPtr          last_desc, new_desc;
  Boolean          last_partial, new_partial, partial3, partial5;

  if (list == NULL || *list == NULL) return;
  last_cds_fcp = NULL;
  last_desc = NULL;
  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS
      || (fcp = vnp->data.ptrvalue) == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }

    ConsolidateClauses (&(fcp->featlist), bsp, biomol, FALSE, rp);

    if (last_cds_fcp == NULL)
    {
      last_cds_fcp = fcp;
      if (fcp->feature_label_data.description == NULL)
      {
        last_desc = GetGeneProtDescription (fcp, bsp, rp);
      }
      else
      {
        last_desc = StringSave (fcp->feature_label_data.description);
      }
      CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
      if (partial5 || partial3) 
      {
        last_partial = TRUE;
      }
      else
      {
        last_partial = FALSE;
      }
    }
    else 
    {
      if (fcp->feature_label_data.description == NULL)
      {
        new_desc = GetGeneProtDescription (fcp, bsp, rp);
      }
      else
      {
        new_desc = StringSave (fcp->feature_label_data.description);
      }
      CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
      if (partial5 || partial3) 
      {
        new_partial = TRUE;
      }
      else
      {
        new_partial = FALSE;
      }
      if ( OkToConsolidate (last_desc, new_desc,
                            last_partial, new_partial,
                            last_cds_fcp, fcp))
      {
        /* two clauses have identical descriptions - combine them */
        MoveSubclauses (last_cds_fcp, fcp);
        RemoveRedundantGeneFeatureFromConsolidatedClause (last_cds_fcp);
        fcp->featlist = NULL;
        fcp->delete_me = TRUE;
        last_cds_fcp->slp = SeqLocMerge (bsp, last_cds_fcp->slp, fcp->slp,
                                         FALSE, TRUE, FALSE);
        /* if we have two clauses that are really identical instead of
         * just sharing a "prefix", make the description plural
         */
        if (StringCmp (last_cds_fcp->interval, fcp->interval) == 0)
        {
          last_cds_fcp->make_plural = TRUE;
/*          PluralizeConsolidatedClauseDescription (last_cds_fcp); */
        }

        /* Recalculate the interval */
        if (last_cds_fcp->interval != NULL)
        {
          MemFree (last_cds_fcp->interval);
        }
        last_cds_fcp->interval =
                  GetGenericInterval (last_cds_fcp, biomol, bsp, rp);
        MemFree (new_desc);
      }
      else
      {
        MemFree (last_desc);
        last_desc = new_desc;
        last_cds_fcp = fcp;
        last_partial = new_partial;
      }
    }  
  }   
  if (delete_now) 
  {
    DeleteFeatureClauses (list);
  }
} 

static void CountUnknownGenes (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp, new_fcp;
  ValNodePtr vnp, new_vnp;
  CharPtr gene_name;
  Int4 num_unknown_genes;
  
  num_unknown_genes = 0;
  vnp = *clause_list;
  new_vnp = NULL;
  new_fcp = NULL;
  while (vnp != NULL)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && ! fcp->is_unknown) 
    {
      CountUnknownGenes (&(fcp->featlist), bsp, rp);
      gene_name = GetGeneProtDescription (fcp, bsp, rp);
      if (StringCmp (gene_name, "unknown") == 0
        && fcp->featlist != NULL
        && fcp->featlist->choice == DEFLINE_FEATLIST)
      {
        if (new_fcp == NULL)
        {
          new_vnp = ValNodeNew (*clause_list);
          if (new_vnp == NULL) return;
          new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue, 
                                      bsp, rp);
          new_fcp->is_unknown = TRUE;
          new_vnp->choice = DEFLINE_CLAUSEPLUS;
          new_vnp->data.ptrvalue = new_fcp;
        }
        else
        {
          new_vnp = ValNodeNew (new_fcp->featlist);
          if (new_vnp == NULL) return;
          new_vnp->choice = DEFLINE_FEATLIST;
          new_vnp->data.ptrvalue = fcp->featlist->data.ptrvalue;
        }  
        num_unknown_genes ++;
        fcp->delete_me = TRUE;
      }
    }
    vnp = vnp->next;
  }
   
  if (num_unknown_genes > 0)
  {
    DeleteFeatureClauses (clause_list);
    if (num_unknown_genes > 1)
    {
      new_fcp->feature_label_data.typeword = StringSave ("genes");
    }
  }
}

static void ReplaceDefinitionLine (
  SeqEntryPtr sep,
  CharPtr defline
)
{
  ValNodePtr ttl;
  if (sep == NULL || defline == NULL) return;

  ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (ttl == NULL)
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
  if (ttl != NULL) {
    MemFree (ttl->data.ptrvalue);
    ttl->data.ptrvalue = defline;
    defline = NULL;
  }
  MemFree (defline);
}

FeatureClausePtr NewFeatureClause 
( SeqFeatPtr sfp,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  Boolean          partial5, partial3;

  fcp = (FeatureClausePtr) MemNew (sizeof (FeatureClauseData));
  if (fcp == NULL) return NULL;

  fcp->feature_label_data.typeword = NULL;
  fcp->feature_label_data.description = NULL;
  fcp->feature_label_data.productname = NULL;
  fcp->feature_label_data.pluralizable = FALSE;
  fcp->feature_label_data.is_typeword_first = FALSE;
  fcp->allelename = NULL;
  fcp->interval = NULL;
  fcp->featlist = NULL;
  fcp->delete_me = FALSE;
  fcp->clause_info_only = FALSE;
  fcp->make_plural = FALSE;
  fcp->is_unknown = FALSE;
  fcp->grp = NULL;
  if (sfp == NULL) return fcp;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  fcp->slp = SeqLocMerge (bsp, sfp->location, NULL,
                                 FALSE, TRUE, FALSE);
  SetSeqLocPartial (fcp->slp, partial5, partial3);
  
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    fcp->grp = sfp->data.value.ptrvalue;
  }
  else
  { 
    fcp->grp = SeqMgrGetGeneXref (sfp);
  }
  if (( IsCDS (sfp) || IsExon (sfp) || IsNoncodingProductFeat (sfp))
    && StringStr (sfp->comment, "alternatively spliced") != NULL) 
  {
    fcp->is_alt_spliced = TRUE;
  }
  else
  {
    fcp->is_alt_spliced = FALSE;
  }
  if (IsCDS (sfp))
  {
    fcp->feature_label_data.productname = GetProductName (sfp, bsp, rp);
  }
  fcp->featlist = ValNodeNew (NULL);
  if (fcp->featlist == NULL)
  {
    MemFree (fcp);
    return NULL;
  }

  fcp->featlist->data.ptrvalue = sfp;
  fcp->featlist->choice = DEFLINE_FEATLIST;
  
  return fcp;
}

static ValNodePtr GetFeatureList (BioseqPtr bsp, DeflineFeatureRequestListPtr rp)
{
  ValNodePtr        head, vnp;
  SeqFeatPtr        sfp;
  FeatureClausePtr  fcp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL) return NULL;

  /* get list of all recognized features */
  head = NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL)
  {
    if (IsRecognizedFeature (sfp))
    {
      fcp = NewFeatureClause (sfp, bsp, rp);
      if (fcp == NULL) return NULL;
      fcp->numivals = fcontext.numivals;
      fcp->ivals = fcontext.ivals;
      vnp = ValNodeNew (head);
      if (head == NULL) head = vnp;
      if (vnp == NULL) return NULL;
      vnp->data.ptrvalue = fcp;
      vnp->choice = DEFLINE_CLAUSEPLUS;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
  return head;
}

static void ExtractSegmentClauses (
  ValNodePtr segment_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR segment_clauses
);

static Boolean FeatureIsOnSegment (
  SeqFeatPtr sfp, 
  ValNodePtr segment_features
)
{
  ValNodePtr vnp, featclause;
  FeatureClausePtr fcp;

  for (vnp = segment_features; vnp != NULL; vnp = vnp->next)
  {
    fcp = vnp->data.ptrvalue;
    if (fcp != NULL)
    {
      for (featclause = fcp->featlist;
           featclause != NULL;
           featclause = featclause->next)
      {
        if (featclause->data.ptrvalue == sfp) return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean FeatureClauseIsOnSegment (
  FeatureClausePtr fcp,
  ValNodePtr segment_features
)
{
  ValNodePtr vnp;

  if (fcp == NULL || fcp->featlist == NULL) return FALSE;
  for (vnp = fcp->featlist; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && FeatureIsOnSegment (vnp->data.ptrvalue, segment_features))
    {
      return TRUE;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      if (FeatureClauseIsOnSegment (vnp->data.ptrvalue, segment_features))
      {
        return TRUE;
      }     
    }
  }
  return FALSE;
}

static FeatureClausePtr CopyMatchingClauses (
  FeatureClausePtr fcp,
  ValNodePtr segment_features
)
{
  FeatureClausePtr new_fcp, copy_clause;
  ValNodePtr       featlist, new_feat;
  Boolean          found_feat_on_segment;
  Boolean          partial5, partial3;

  new_fcp = (FeatureClausePtr) MemNew (sizeof (FeatureClauseData));
  if (new_fcp == NULL) return NULL;
  new_fcp->feature_label_data.pluralizable = 
    fcp->feature_label_data.pluralizable;
  new_fcp->feature_label_data.is_typeword_first = 
    fcp->feature_label_data.is_typeword_first;
  new_fcp->feature_label_data.typeword =
    StringSave (fcp->feature_label_data.typeword);
  new_fcp->feature_label_data.description = 
    StringSave (fcp->feature_label_data.description);
  new_fcp->feature_label_data.productname =
    StringSave (fcp->feature_label_data.productname);
  new_fcp->allelename = StringSave (fcp->allelename);
  new_fcp->interval = StringSave (fcp->interval);
  new_fcp->is_alt_spliced = fcp->is_alt_spliced;
  CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
  new_fcp->slp = (SeqLocPtr) AsnIoMemCopy (fcp->slp, (AsnReadFunc) SeqLocAsnRead,
                                           (AsnWriteFunc) SeqLocAsnWrite);
  SetSeqLocPartial (new_fcp->slp, partial5, partial3);
  new_fcp->grp = fcp->grp;
  new_fcp->has_mrna = fcp->has_mrna;
  new_fcp->delete_me = FALSE;
  new_fcp->clause_info_only = fcp->clause_info_only;
  new_fcp->featlist = NULL;
  found_feat_on_segment = FALSE;
  for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
  {
    new_feat = NULL;
    if (featlist->choice == DEFLINE_FEATLIST)
    {
      new_feat = ValNodeNew (new_fcp->featlist);
      if (new_feat == NULL) return NULL;
      new_feat->data.ptrvalue = featlist->data.ptrvalue;
      new_feat->choice = DEFLINE_FEATLIST;

      /* some portions of the clause are present for product and gene info */
      /* if they aren't actually on the segment */
      if ( segment_features == NULL
        || FeatureIsOnSegment (new_feat->data.ptrvalue, segment_features))
      {
        found_feat_on_segment = TRUE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      copy_clause = featlist->data.ptrvalue;
      if ( segment_features == NULL
        || FeatureClauseIsOnSegment ( copy_clause, segment_features))
      {
        new_feat = ValNodeNew (new_fcp->featlist);
        if (new_feat == NULL) return NULL;
        new_feat->data.ptrvalue = CopyMatchingClauses ( copy_clause,
                                                        segment_features);
        new_feat->choice = DEFLINE_CLAUSEPLUS;
      }
    }
    if (new_feat != NULL && new_fcp->featlist == NULL) 
    {
      new_fcp->featlist = new_feat;
    }
  }
  if (found_feat_on_segment)
  {
    new_fcp->clause_info_only = FALSE;
  }
  else
  {
    new_fcp->clause_info_only = TRUE;
  }
  return new_fcp; 
}

static void CopyFeatureList (
  ValNodePtr match_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR new_list
)
{
  ValNodePtr vnp, addvnp;

  for (vnp = parent_features;
       vnp != NULL;
       vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && (match_features == NULL
        || FeatureIsOnSegment ( vnp->data.ptrvalue, match_features)))
    {
      addvnp = ValNodeNew (*new_list);
      if (addvnp == NULL) return;
      addvnp->data.ptrvalue = vnp->data.ptrvalue;
      addvnp->choice = DEFLINE_FEATLIST;
      if (*new_list == NULL) *new_list = addvnp;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (match_features == NULL
        || FeatureClauseIsOnSegment ( vnp->data.ptrvalue, match_features)))
    {
      addvnp = ValNodeNew (*new_list);
      if (addvnp == NULL) return;
      addvnp->data.ptrvalue = CopyMatchingClauses ( vnp->data.ptrvalue,
                                                    match_features);
      addvnp->choice = DEFLINE_CLAUSEPLUS;
      if (*new_list == NULL) *new_list = addvnp;
    }
  }
    
}

static void ExtractSegmentClauses (
  ValNodePtr segment_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR segment_clauses
)
{
  CopyFeatureList (segment_features, parent_features, segment_clauses);
}

typedef struct segmentdeflinedata {
  BioseqPtr  parent_bsp;
  ValNodePtr parent_feature_list;
  Uint1      molecule_type;
  DeflineFeatureRequestList PNTR feature_requests;
  ModifierCombinationPtr m;
  ModifierItemLocalPtr modList;
  OrganismDescriptionModifiersPtr odmp;
  Int2 product_flag;
} SegmentDefLineData, PNTR SegmentDefLinePtr;

typedef struct segmentdeflinefeatureclausedata {
  BioseqPtr  parent_bsp;
  ValNodePtr parent_feature_list;
  Uint1      molecule_type;
  DeflineFeatureRequestList PNTR feature_requests;
  Int2            product_flag;
  Boolean         alternate_splice_flag;
  Boolean         gene_cluster_opp_strand;
  ValNodePtr PNTR list;
} SegmentDefLineFeatureClauseData, PNTR SegmentDefLineFeatureClausePtr;

typedef struct deflinefeatclause {
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  CharPtr     clauselist;
} DefLineFeatClauseData, PNTR DefLineFeatClausePtr;

NLM_EXTERN void DefLineFeatClauseListFree (ValNodePtr vnp)
{
  DefLineFeatClausePtr deflist;

  if (vnp == NULL) return;
  DefLineFeatClauseListFree (vnp->next);
  vnp->next = NULL;
  deflist = vnp->data.ptrvalue;
  if (deflist != NULL)
  {
    MemFree (deflist->clauselist);
    MemFree (deflist);
  }
  ValNodeFree (vnp);
}


static Boolean IntervalIntersectsIvals 
(Int2    numivals,
 Int4Ptr ivals,
 SeqMgrSegmentContextPtr context)
{
  Int2 idx;
  Int4 start, stop;

  if (numivals == 0 || ivals == NULL || context == NULL) return FALSE;

  for (idx = 0; idx < numivals; idx ++) {
    start = ivals [idx * 2];
    stop = ivals [idx * 2 + 1];
    if ( start <= context->cumOffset + context->to - context->from
         && stop >= context->cumOffset)
    {
      return TRUE;
    }
  }
  return FALSE;
}


/* if there are no features at all on this segment, select the genes that 
 * traverse the segment.
 */
static ValNodePtr GrabTraversingGenes 
(ValNodePtr              parent_feature_list,
 SeqMgrSegmentContextPtr context,
 BioseqPtr               parent_bsp,
 DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr  fcp, new_fcp;
  ValNodePtr        clause;
  SeqFeatPtr        sfp;
  ValNodePtr        segment_feature_list;
  ValNodePtr        vnp;

  segment_feature_list = NULL;
  for (clause = parent_feature_list;
       clause != NULL;
       clause = clause->next)
  {
    fcp = clause->data.ptrvalue;

    if (fcp != NULL  &&  fcp->featlist != NULL
        &&  (sfp = fcp->featlist->data.ptrvalue) != NULL
        &&  sfp->idx.subtype == FEATDEF_GENE
        &&  fcp->ivals != NULL && fcp->numivals > 0)
    {
      if (IntervalIntersectsIvals (fcp->numivals, fcp->ivals, context)) {
        new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue, parent_bsp,
                                    rp);
        if (new_fcp == NULL) return FALSE;
        vnp = ValNodeNew (segment_feature_list);
        if (vnp == NULL) return FALSE;
        vnp->data.ptrvalue = new_fcp;
        vnp->choice = DEFLINE_CLAUSEPLUS;
        if (segment_feature_list == NULL) segment_feature_list = vnp;
      }
    } 
  }
  return segment_feature_list;
}


static CharPtr BuildFeatureClauses (
  BioseqPtr   bsp,
  Uint1       molecule_type,
  SeqEntryPtr sep,
  ValNodePtr  PNTR feature_list,
  Boolean     isSegment,
  ValNodePtr  PNTR seg_feature_list,
  Int2        product_flag,
  Boolean     alternate_splice_flag,
  Boolean     gene_cluster_opp_strand,
  DeflineFeatureRequestList PNTR feature_requests
);

static Boolean LIBCALLBACK GetFeatureClauseForSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context)
{
  SegmentDefLineFeatureClausePtr sdlp;
  ValNodePtr        clause, tmp_parent_list;
  FeatureClausePtr  fcp, new_fcp;
  Int2              idx;
  Int4              start, stop;
  ValNodePtr        segment_feature_list, vnp;
  SeqIdPtr          sip;
  BioseqPtr         bsp;
  Uint2             entityID;
  SeqLocPtr         loc;
  DefLineFeatClausePtr deflist;

  if (slp == NULL || context == NULL) return FALSE;
  sdlp = (SegmentDefLineFeatureClausePtr) context->userdata;

  sip = SeqLocId (slp);
  
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFind (sip);

  if (bsp == NULL) return TRUE;


  segment_feature_list = NULL;
  for (clause = sdlp->parent_feature_list;
       clause != NULL;
       clause = clause->next)
  {
    fcp = clause->data.ptrvalue;

    if (fcp != NULL && fcp->ivals != NULL && fcp->numivals > 0)
    {
      idx = (fcp->numivals - 1) * 2;
      start = fcp->ivals [idx];
      stop = fcp->ivals [idx + 1];
      if ( stop <= context->cumOffset + context->to - context->from
           && stop >= context->cumOffset)
      {
        new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue,
                                    sdlp->parent_bsp, 
                                    sdlp->feature_requests);
        if (new_fcp == NULL) return FALSE;
        vnp = ValNodeNew (segment_feature_list);
        if (vnp == NULL) return FALSE;
        vnp->data.ptrvalue = new_fcp;
        vnp->choice = DEFLINE_CLAUSEPLUS;
        if (segment_feature_list == NULL) segment_feature_list = vnp;
      }
    } 
  }

  if (segment_feature_list == NULL) {
    segment_feature_list = GrabTraversingGenes (sdlp->parent_feature_list,
                                                context, sdlp->parent_bsp,
                                                sdlp->feature_requests);
  }

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  tmp_parent_list = NULL;
  CopyFeatureList (NULL, sdlp->parent_feature_list, &tmp_parent_list);

  deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
  if (deflist == NULL) return TRUE;
  deflist->sep = SeqMgrGetSeqEntryForData (bsp);
  deflist->bsp = bsp;
  deflist->clauselist = BuildFeatureClauses (sdlp->parent_bsp,
                            sdlp->molecule_type,
                            SeqMgrGetSeqEntryForData (bsp),
                            &tmp_parent_list,
                            TRUE,
                            &segment_feature_list,
                            sdlp->product_flag,
							              sdlp->alternate_splice_flag,
							              sdlp->gene_cluster_opp_strand,
                            sdlp->feature_requests);
  vnp = ValNodeNew (*(sdlp->list));
  if (vnp == NULL) return TRUE;
  if (*(sdlp->list) == NULL) *(sdlp->list) = vnp;
  vnp->data.ptrvalue = deflist;

  FreeListElement (tmp_parent_list);
  FreeListElement (segment_feature_list);
  DeleteMarkedObjects (entityID, 0, NULL);
  return TRUE;
}

static Boolean HasAnyPromoters (BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_promoter, &fcontext);
  if (sfp == NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static void AddFakePromoterClause (ValNodePtr PNTR feature_list, BioseqPtr bsp, DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr       sfp = NULL;
  
  /* create fake promoter */
  sfp = CreateNewFeature (SeqMgrGetSeqEntryForData (bsp), NULL,
                          SEQFEAT_IMP, NULL);

  sfp->location = SeqLocIntNew (0, bsp->length - 1, Seq_strand_plus, SeqIdDup (bsp->id));
  sfp->data.choice = SEQFEAT_IMP;
  sfp->idx.subtype = FEATDEF_promoter;
  /* mark promoter for deletion */
  sfp->idx.deleteme = TRUE;
  
  fcp = NewFeatureClause (sfp, bsp, rp);
  if (fcp == NULL) return;
  fcp->numivals = 1;
  fcp->ivals = (Int4Ptr) MemNew (sizeof (Int4) * 2);
  fcp->ivals[0] = 0;
  fcp->ivals[1] = bsp->length - 1;
  ValNodeAddPointer (feature_list, DEFLINE_CLAUSEPLUS, fcp);
  
}

static Boolean IsInGenProdSet (BioseqPtr bsp)
{
  BioseqSetPtr bssp;
  if (bsp == NULL || bsp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return FALSE;
  }
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  if (bssp->_class != BioseqseqSet_class_nuc_prot || bssp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return FALSE;
  }
  bssp = bssp->idx.parentptr;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/* NOTE: under some circumstances this function will create features that
 * are marked for deletion, so DeleteMarkedObjects should always be called
 * at some later point.
 */
static CharPtr BuildFeatureClauses (
  BioseqPtr   bsp,
  Uint1       molecule_type,
  SeqEntryPtr sep,
  ValNodePtr  PNTR feature_list,
  Boolean     isSegment,
  ValNodePtr  PNTR seg_feature_list,
  Int2        product_flag,
  Boolean     alternate_splice_flag,
  Boolean     gene_cluster_opp_strand,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  ValNodePtr   strings = NULL;
  ValNodePtr   clause;
  CharPtr      str = NULL;
  Char         ending_str [200];
  ValNodePtr   tmp_feat_list;
  BioSourcePtr biop;

  if ((feature_requests->feature_list_type == DEFLINE_USE_FEATURES
       || (IsmRNASequence(bsp) && IsInGenProdSet(bsp)))
      && (! isSegment || (seg_feature_list != NULL && *seg_feature_list != NULL)))
  {
    /* remove features that indexer has chosen to suppress before they are grouped
     * with other features or used to determine loneliness etc.
     */
    RemoveSuppressedFeatures (feature_list, feature_requests->suppressed_feature_list);
  
    GroupmRNAs (feature_list, bsp, feature_requests);

    /* genes are added to other clauses */
    GroupGenes (feature_list, feature_requests->suppress_locus_tags);

    if (! feature_requests->suppress_alt_splice_phrase)
    {
      /* find alt-spliced CDSs */
      FindAltSplices (*feature_list, bsp, feature_requests);
    }

    GroupAltSplicedExons (feature_list, bsp, TRUE);
    
    if (!isSegment)
    {
       /* group CDSs that have the same name and are under the same gene together */
      GroupSegmentedCDSs (feature_list, bsp, TRUE, feature_requests);
    }

    /* per Susan's request, if promoters are requested and no promoters are found, add a promoter */
    if (feature_requests->keep_items[RemovablePromoter] 
        && feature_requests->add_fake_promoters
        && !HasAnyPromoters (bsp)) {
      AddFakePromoterClause (feature_list, bsp, feature_requests);
    }

    /* now group clauses */
    GroupAllClauses ( feature_list, gene_cluster_opp_strand, bsp );

    ExpandAltSplicedExons (*feature_list, bsp, feature_requests);

    FindGeneProducts (*feature_list, bsp, feature_requests);

    if (seg_feature_list != NULL && *seg_feature_list != NULL)
    {
      tmp_feat_list = NULL; 
      ExtractSegmentClauses ( *seg_feature_list, *feature_list, &tmp_feat_list);
      FreeListElement (*feature_list);
      DeleteMarkedObjects (bsp->idx.entityID, 0, NULL);
      *feature_list = tmp_feat_list;
    }
   
    /* remove exons and other unwanted features */
    RemoveUnwantedFeatures (feature_list, bsp, isSegment, feature_requests);

    RemoveGenesMentionedElsewhere (feature_list, *feature_list, TRUE,
                                   feature_requests->suppress_locus_tags);

    if (feature_requests->remove_subfeatures)
    {
      DeleteSubfeatures (feature_list, TRUE);
    }

    DeleteOperonAndGeneClusterSubfeatures (feature_list, TRUE);

    CountUnknownGenes (feature_list, bsp, feature_requests);

    if (feature_requests->misc_feat_parse_rule == 1)
    {
      RenameMiscFeats (*feature_list, molecule_type);
    }
    else
    {
      RemoveUnwantedMiscFeats (feature_list, TRUE);
    }

    ReplaceRNAClauses (feature_list, bsp, feature_requests);

    /* take any exons on the minus strand */
    /* and reverse their order within the clause */
    ReverseClauses (feature_list, IsExonOrIntron);

    RenameExonSequences ( feature_list, bsp, TRUE);

    LabelClauses (*feature_list, molecule_type, bsp, 
                  feature_requests);

    /* parse lists of tRNA and intergenic spacer clauses in misc_feat notes */
    /* need to do this after LabelClauses, since LabelClauses labels intergenic
     * spacers with more relaxed restrictions.  The labels from LabelClauses
     * for intergenic spacers are the default values.
     */
    ReplaceIntergenicSpacerClauses (feature_list, bsp, feature_requests);

    ConsolidateClauses (feature_list, bsp, molecule_type, TRUE,
                        feature_requests);

    /* this allows genes to be listed together even if they are from */
    /* separate sequences */
/*    SmashTallClauses (feature_list, TRUE); */

    clause = *feature_list;
    ListClauses (clause, &strings, TRUE, FALSE);

    AutoDef_AddEnding (clause, &strings, bsp, 
                       product_flag, alternate_splice_flag);
    str = MergeValNodeStrings (strings, FALSE);
	  ValNodeFreeData (strings);
  }
  else if (feature_requests->feature_list_type == DEFLINE_COMPLETE_SEQUENCE)
  {
    str = StringSave (", complete sequence.");
  }
  else if (feature_requests->feature_list_type == DEFLINE_COMPLETE_GENOME)
  {
    ending_str [0] = 0;
    biop = GetBiopForBsp (bsp);
    if (biop != NULL)
    {
      switch (biop->genome) {
        case GENOME_macronuclear :
          sprintf (ending_str, "macronuclear");
          break;
        case GENOME_nucleomorph :
          sprintf (ending_str, "nucleomorph");
          break;
        case GENOME_mitochondrion :
          sprintf (ending_str, "mitochondrion");
          break;
        case GENOME_apicoplast :
        case GENOME_chloroplast :
        case GENOME_chromoplast :
        case GENOME_kinetoplast :
        case GENOME_plastid :
        case GENOME_cyanelle :
        case GENOME_leucoplast :
        case GENOME_proplastid :
        case GENOME_hydrogenosome :
          sprintf (ending_str, "%s", organelleByGenome [biop->genome]);
          break;
      }
    }
    StringCat (ending_str, ", complete genome.");
    str = StringSave (ending_str);
  }
  else
  {
    str = StringSave ("");
  }
  
  return str;
}

/* This function looks at the product names for the CDSs on the Bioseq,
 * and sets the flag for the "nuclear genes for X products" ending
 * based on the contents of the CDS products. */
static Int2 GetProductFlagFromCDSProductNames (BioseqPtr bsp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        cds = NULL;
  Int2              product_flag;
  Int2              i;
  CharPtr           found;
  Char              ch_before, ch_after;

  product_flag = 0;
  for (cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &context);
       cds != NULL && product_flag == 0;
       cds = cds->next)
  {
    for (i = 1; organelleByPopup[i] != NULL && product_flag == 0; i++)
    {
      found = StringStr (context.label, organelleByPopup[i]);
      if (found != NULL)
      {
        if (found == context.label) {
          ch_before = ' ';
        } else {
          ch_before = *(found - 1);
        }
        ch_after = *(found + StringLen (organelleByPopup[i]));
        if (ch_before == ' ' && (ch_after == 0 || ch_after == ' '))
        {
          product_flag = i;
        }
      }
    }
  }

  return product_flag;
}

NLM_EXTERN void BuildDefLineFeatClauseList (
  SeqEntryPtr sep,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list
)
{
  BioseqSetPtr    bssp;
  BioseqPtr    bsp;
  ValNodePtr    head;
  Uint1      molecule_type;
  SeqEntryPtr   nsep;
  SegmentDefLineFeatureClauseData sdld;
  DefLineFeatClausePtr deflist;
  ValNodePtr    vnp;

  if (sep == NULL || list == NULL) return;

  if ( IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if ( bssp->_class == 7 || IsPopPhyEtcSet (bssp->_class)
        || bssp->_class == BioseqseqSet_class_gen_prod_set
        || bssp->_class == BioseqseqSet_class_not_set)
    {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        BuildDefLineFeatClauseList (sep, entityID, feature_requests,
                                    product_flag, alternate_splice_flag,
                                    gene_cluster_opp_strand, list);
      }
      return;
    }
  }
    
  nsep = FindNucSeqEntry (sep);
  if (nsep != NULL)
  {
    bsp = (BioseqPtr) nsep->data.ptrvalue;
    if ( SpecialHandlingForSpecialTechniques (bsp)) 
    {
      return;
    }
    if (bsp != NULL && bsp->repr == Seq_repr_seg &&
      bsp->seq_ext != NULL && bsp->seq_ext_type == 1) 
    {
      /* get default product flag if necessary */
	  if (product_flag == -1 || product_flag == DEFAULT_ORGANELLE_CLAUSE) {
	    product_flag = GetProductFlagFromCDSProductNames (bsp);
	  }

      sdld.parent_bsp = bsp;
      sdld.molecule_type = GetMoleculeType (bsp, entityID);
      sdld.parent_feature_list = GetFeatureList (bsp, 
                                                 feature_requests);

      sdld.feature_requests =  feature_requests;
      sdld.product_flag = product_flag;
      sdld.alternate_splice_flag = alternate_splice_flag;
      sdld.gene_cluster_opp_strand = gene_cluster_opp_strand;
      sdld.list = list;
      SeqMgrExploreSegments (bsp, (Pointer) &sdld, GetFeatureClauseForSeg);
      deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
      if (deflist == NULL) return;
      deflist->sep = SeqMgrGetSeqEntryForData (bsp),
      deflist->bsp = bsp;

      deflist->clauselist = BuildFeatureClauses (bsp,
                            sdld.molecule_type,
                            SeqMgrGetSeqEntryForData (bsp),
                            &sdld.parent_feature_list,
                            FALSE,
                            NULL,
                            product_flag,
                            alternate_splice_flag,
                            gene_cluster_opp_strand,
                            feature_requests);
      vnp = ValNodeNew (*list);
      if (vnp == NULL) return;
      if (*list == NULL) *list = vnp;
      vnp->data.ptrvalue = deflist;
      FreeListElement (sdld.parent_feature_list);
      DeleteMarkedObjects (entityID, 0, NULL);
      return;
    }
  }

  if (nsep != NULL && nsep != sep)
    sep = nsep;


  if (! IS_Bioseq (sep)) return;

  /* get list of all recognized features */
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if ( SpecialHandlingForSpecialTechniques (bsp)) 
  {
    return;
  }
  molecule_type = GetMoleculeType (bsp, entityID);
  head = GetFeatureList (bsp, feature_requests);

  /* get default product flag if necessary */
  if (product_flag == -1 || product_flag == DEFAULT_ORGANELLE_CLAUSE) {
    product_flag = GetProductFlagFromCDSProductNames (bsp);
  }

  deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
  if (deflist == NULL) return;
  deflist->sep = SeqMgrGetSeqEntryForData (bsp),
  deflist->bsp = bsp;
  deflist->clauselist = BuildFeatureClauses (bsp,
                                             molecule_type,
                                             SeqMgrGetSeqEntryForData (bsp),
                                             &head,
                                             FALSE,
                                             NULL,
                                             product_flag,
                                             alternate_splice_flag,
                                             gene_cluster_opp_strand,
                                             feature_requests);
  vnp = ValNodeNew (*list);
  if (vnp == NULL) return;
  if (*list == NULL) *list = vnp;
  vnp->data.ptrvalue = deflist;
  FreeListElement (head);
  DeleteMarkedObjects (entityID, 0, NULL);
}

static Boolean IdenticalExceptForPartialComplete (CharPtr str1, CharPtr str2)
{
    CharPtr cp, word_in_first, word_in_second;
    Int4    first_len, second_len, compare_len;
    
    if (StringHasNoText (str1) && StringHasNoText (str2)) {
        return TRUE;
    } else if (StringHasNoText (str1) || StringHasNoText (str2)) {
        return FALSE;
    }
    
    word_in_first = StringISearch (str1, "partial");
    cp = StringISearch (str1, "complete");
    if (word_in_first == NULL || (cp != NULL && word_in_first > cp)) {
        word_in_first = cp;
        first_len = 8;
    } else {
        first_len = 7;
    }
    
    word_in_second = StringISearch (str2, "partial");
    cp = StringISearch (str2, "complete");
    if (word_in_second == NULL || (cp != NULL && word_in_second > cp)) {
        word_in_second = cp;
        second_len = 8;
    } else {
        second_len = 7;
    }
    
    if (word_in_first == NULL && word_in_second == NULL) {
        if (StringCmp (str1, str2) == 0) {
            return TRUE;
        } else {
            return FALSE;
        }
    } else if (word_in_first == NULL || word_in_second == NULL) {
        return FALSE;
    } else if ((compare_len = word_in_first - str1) != word_in_second - str2) {
        return FALSE;
    } else if (StringNCmp (str1, str2, compare_len) != 0) {
        return FALSE;
    } else {
        return IdenticalExceptForPartialComplete (word_in_first + first_len, word_in_second + second_len);
    }    
}


static CharPtr GetTaxnameForBsp (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;
  CharPtr           taxname = NULL;

  if (bsp != NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL && (biop = sdp->data.ptrvalue) != NULL
      && biop->org != NULL) {
      taxname = biop->org->taxname;
      }
  }
  return taxname;
}


NLM_EXTERN Boolean AreFeatureClausesUnique (ValNodePtr list)
{
  ValNodePtr vnp1, vnp2;
  DefLineFeatClausePtr deflist1, deflist2;
  CharPtr              taxname1;

  for (vnp1 = list; vnp1 != NULL && vnp1->next != NULL; vnp1 = vnp1->next)
  {
    deflist1 = vnp1->data.ptrvalue;
    if (deflist1 == NULL || deflist1->clauselist == NULL) return FALSE;
    taxname1 = GetTaxnameForBsp (deflist1->bsp);
    
    for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next)
    {
      deflist2 = vnp2->data.ptrvalue;
      if (deflist2 == NULL || deflist2->clauselist == NULL
        || (StringCmp (taxname1, GetTaxnameForBsp (deflist2->bsp)) == 0 
            && IdenticalExceptForPartialComplete (deflist1->clauselist, deflist2->clauselist)))
      {
        return FALSE;
      }
    }
  }
  return TRUE;
}


static CharPtr GetKeywordPrefix (SeqEntryPtr sep)
{
  ValNodePtr vnp;
  GBBlockPtr gbp;
  MolInfoPtr mip;
  
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_genbank, NULL);
  if (vnp == NULL) {
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL && mip->tech == MI_TECH_tsa) {
        return "TSA: ";
      }
    }
	} else {
	  gbp = (GBBlockPtr) vnp->data.ptrvalue;
	  if (gbp != NULL)
	  {
	    for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next)
	    {
	      if (StringICmp((CharPtr)vnp->data.ptrvalue, "TPA:inferential") == 0)
	      {
	        return "TPA_inf: ";
	      }
	      else if (StringICmp((CharPtr)vnp->data.ptrvalue, "TPA:experimental") == 0)
	      {
	        return "TPA_exp: ";
	      }
	    }
	  }
	}
	return "";
}

NLM_EXTERN void BuildDefinitionLinesFromFeatureClauseLists (
  ValNodePtr list,
  ModifierItemLocalPtr modList,
  ValNodePtr modifier_indices,
  OrganismDescriptionModifiersPtr odmp
)
{
  ValNodePtr vnp;
  DefLineFeatClausePtr defline_featclause;
  CharPtr    org_desc, tmp_str, keyword_prefix;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      defline_featclause = vnp->data.ptrvalue;
      
      keyword_prefix = GetKeywordPrefix (defline_featclause->sep);
      
      org_desc = GetOrganismDescription (defline_featclause->bsp,
                                         modList, modifier_indices, odmp);
      tmp_str = (CharPtr) MemNew (StringLen (keyword_prefix) 
                                  + StringLen (org_desc) 
                                  + StringLen (defline_featclause->clauselist) + 2);
      if (tmp_str == NULL) return;
      tmp_str [0] = 0;
      if (keyword_prefix != NULL)
      {
        StringCat (tmp_str, keyword_prefix);
      }
      StringCat (tmp_str, org_desc);
      if (defline_featclause->clauselist != NULL
        && defline_featclause->clauselist [0] != ','
        && defline_featclause->clauselist [0] != '.'
        && defline_featclause->clauselist [0] != 0)
      {
        StringCat (tmp_str, " ");
      }
      StringCat (tmp_str, defline_featclause->clauselist);
      tmp_str [0] = toupper (tmp_str [0]);
      ReplaceDefinitionLine (defline_featclause->sep, tmp_str);
      MemFree (org_desc);
    }
  }
}


/* This removes redundant titles on nuc-prot sets, which will not be
 * visible in the flat file if all sequences in the nuc-prot set have
 * their own title.
 */
NLM_EXTERN void RemoveNucProtSetTitles (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  this_sep;
  SeqDescrPtr  sdp, prev = NULL, sdp_next;
  
  if (sep == NULL || ! IS_Bioseq_set (sep))
  {
    return;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  for (this_sep = bssp->seq_set; this_sep != NULL; this_sep = this_sep->next)
  {
    RemoveNucProtSetTitles (this_sep);
  }
  
  if (bssp->_class != BioseqseqSet_class_nuc_prot) 
  {
    return;
  }
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp_next)
  {
    sdp_next = sdp->next;
    if (sdp->choice == Seq_descr_title)
    {
      if (prev == NULL)
      {
        bssp->descr = sdp->next;
      }
      else
      {
        prev->next = sdp->next;
      }
      sdp->next = NULL;
      SeqDescrFree (sdp);
    }
    else
    {
      prev = sdp;
    }
  }
}


NLM_EXTERN void 
AutoDefForSeqEntry 
(SeqEntryPtr sep,
 Uint2 entityID,
 OrganismDescriptionModifiersPtr odmp,
 ModifierItemLocalPtr modList,
 ValNodePtr modifier_indices,
 DeflineFeatureRequestListPtr feature_requests,
 Int2 product_flag,
 Boolean alternate_splice_flag,
 Boolean gene_cluster_opp_strand)
{
  ValNodePtr defline_clauses = NULL;

  if (sep == NULL) return;

  RemoveNucProtSetTitles (sep);
  
  SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID,
                              feature_requests,
                              product_flag, alternate_splice_flag,
                              gene_cluster_opp_strand,
                              &defline_clauses);
                              
/*  dlfp->feature_requests.suppressed_feature_list = ValNodeFree (dlfp->feature_requests.suppressed_feature_list);                               */

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, odmp);
  DefLineFeatClauseListFree (defline_clauses);
  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
}


/* functions for editing seq-locs */
NLM_EXTERN Int4 ExtendSeqLocToEnd (SeqLocPtr slp, BioseqPtr bsp, Boolean end5)
{
  Uint1          strand;
  SeqLocPtr      slp_to_change, slp_index;
  Int4           extent_to_change;
  Int4           start, stop;
  SeqIdPtr       sip;
  Int4           start_diff = 0;
  
  if (slp == NULL || bsp == NULL) return 0;

  slp_to_change = NULL;
  strand = SeqLocStrand (slp);
  switch (slp->choice)
  {
    case SEQLOC_INT:
      slp_to_change = slp;
      break;
    case SEQLOC_MIX:
  	case SEQLOC_PACKED_INT:
      sip = SeqLocId (slp);
      if (sip == NULL) return 0; /* can only process if all on one bioseq */
      slp_to_change = NULL;
      if ((strand == Seq_strand_minus && end5)
        || (strand != Seq_strand_minus && !end5))
      {
        extent_to_change = 0;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          stop = GetOffsetInBioseq (slp_index, bsp, SEQLOC_STOP);
          if (stop > extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = stop;
          }
        }
      }
      else
      {
        extent_to_change = bsp->length;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          start = GetOffsetInBioseq (slp_index, bsp, SEQLOC_START);
          if (start < extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = start;
          }
        }
      }
      break;
  }

  if (slp_to_change != NULL)
  {
    if ((strand == Seq_strand_minus && end5)
      || (strand != Seq_strand_minus && !end5))
    {
      start = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_START);
      stop = bsp->length - 1;
    }
    else
    {
      start = 0;
      stop = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_STOP);
    }
    if (end5) {
        if (strand == Seq_strand_minus) {
            start_diff = bsp->length - 1 - GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        } else {
            start_diff = GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        }
    }
    
    expand_seq_loc (start, stop, strand, slp_to_change);
  }
  return start_diff;
}


/* functions for feature conversion.  shared by sequin5 and macroapi */

NLM_EXTERN Boolean IsBioseqSetInGPS (BioseqSetPtr bssp)
{
  if (bssp == NULL) return FALSE;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) return TRUE;
  if (bssp->idx.parentptr == NULL || bssp->idx.parenttype != OBJ_BIOSEQSET) return FALSE;
  return IsBioseqSetInGPS ((BioseqSetPtr) bssp->idx.parentptr);
}


NLM_EXTERN Boolean IsBioseqInGPS (BioseqPtr bsp)
{
  if (bsp == NULL || bsp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET)
  {
    return FALSE;
  }
  else
  {
    return IsBioseqSetInGPS ((BioseqSetPtr) bsp->idx.parentptr);
  }
}


NLM_EXTERN Boolean IsFeatInGPS (SeqFeatPtr sfp)
{
  if (sfp == NULL) return FALSE;
  return IsBioseqInGPS (BioseqFindFromSeqLoc (sfp->location));
}


NLM_EXTERN RnaRefPtr RnaRefFromLabel (Uint2 featdef_to, CharPtr label, BoolPtr add_label_to_comment)
{
  RnaRefPtr rrp;
  tRNAPtr   trp = NULL;
  Boolean   just_trna_text;
  Int4      j;

  rrp = RnaRefNew ();
  if (NULL == rrp)
    return NULL;

  *add_label_to_comment = FALSE;

  switch (featdef_to) 
  {
    case FEATDEF_preRNA :
      rrp->type = 1;
      break;
    case FEATDEF_mRNA :
      rrp->type = 2;
      break;
    case FEATDEF_tRNA :
      rrp->type = 3;
      break;
    case FEATDEF_rRNA :
      rrp->type = 4;
      break;
    case FEATDEF_snRNA :
      rrp->type = 5;
      break;
    case FEATDEF_scRNA :
      rrp->type = 6;
      break;
    case FEATDEF_snoRNA :
      rrp->type = 7;
      break;
    case FEATDEF_ncRNA:
    case FEATDEF_tmRNA:
    case FEATDEF_otherRNA :
      rrp->type = 255;
      break;
    default :
      break;
  }

  if (featdef_to == FEATDEF_tRNA) 
  {
    trp = (tRNAPtr) MemNew (sizeof (tRNA));
    rrp->ext.choice = 2;
    rrp->ext.value.ptrvalue = (Pointer) trp;
    trp->aa = ParseTRnaString (label, &just_trna_text, NULL, FALSE);
    trp->aatype = 2;
    for (j = 0; j < 6; j++) {
	    trp->codon [j] = 255;
    }
    if (!just_trna_text)
    {
      *add_label_to_comment = TRUE;
    }
  } 
  else if (featdef_to == FEATDEF_ncRNA) 
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("ncRNA");
  }  
  else if (featdef_to == FEATDEF_tmRNA) 
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("tmRNA");
  } 
  else if (featdef_to == FEATDEF_otherRNA) 
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
  } 
  else if (! StringHasNoText (label))
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave (label);
  }
  return rrp;
}


NLM_EXTERN Boolean ConvertProtToProtFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  ProtRefPtr prp;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
    return FALSE;

  switch (featdef_to) {
    case FEATDEF_PROT :
      prp->processed = 0;
      break;
    case FEATDEF_preprotein :
      prp->processed = 1;
      break;
    case FEATDEF_mat_peptide_aa :
      prp->processed = 2;
      break;
    case FEATDEF_sig_peptide_aa :
      prp->processed = 3;
      break;
    case FEATDEF_transit_peptide_aa :
      prp->processed = 4;
      break;
    default :
      break;
  }
  return TRUE;
}


NLM_EXTERN void 
ApplyCDSOptionsToFeature
(SeqFeatPtr sfp,
 Boolean remove_mRNA,
 Boolean remove_gene,
 Boolean remove_transcript_id,
 Boolean keep_original)
{
  BioseqPtr         product_bsp;
  SeqFeatPtr        gene, mrna;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL) return;

  if (sfp->product != NULL) {
    product_bsp = BioseqFindFromSeqLoc (sfp->product);
    if (product_bsp != NULL && !keep_original)
    {
      product_bsp->idx.deleteme = TRUE;
    }
    if (!IsFeatInGPS (sfp) || remove_transcript_id)
    {
      sfp->product = SeqLocFree (sfp->product);
    }
  }

  if (remove_gene)
  {
    gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
    if (gene != NULL)
    {
      gene->idx.deleteme = TRUE;
    }
  }
  
  if (remove_mRNA)
  {
    mrna = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
    if (mrna != NULL)
    {
      mrna->idx.deleteme = TRUE;
    }
  }

}


NLM_EXTERN Boolean 
ConvertCDSToRNA 
(SeqFeatPtr  sfp,
 Uint2       rna_type)
{
  Char                   label [256];
  CharPtr                new_comment;
  Int4                   comment_len = 0;
  Boolean                add_label_to_comment = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;

  FeatDefLabel (sfp, label, sizeof (label), OM_LABEL_CONTENT);

  sfp->data.value.ptrvalue =
      CdRegionFree ((CdRegionPtr) sfp->data.value.ptrvalue);

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = RnaRefFromLabel (rna_type, label, &add_label_to_comment);
  
  if (add_label_to_comment && StringCmp (label, sfp->comment) != 0)
  {
    if (StringHasNoText (sfp->comment)) 
    {
      new_comment = StringSave (label);
    }
    else
    {
      comment_len = StringLen (sfp->comment) + StringLen (label) + 3;
      new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
      sprintf (new_comment, "%s; %s", sfp->comment, label);
    }
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;
  }
  /* change subtype so that feature will be reindexed */
  sfp->idx.subtype = 0;  
  
  return TRUE;
}


NLM_EXTERN Boolean ConvertGeneToRNA (SeqFeatPtr sfp, Uint2 featdef_to)
{
  Char                   label [256];
  GeneRefPtr grp;
  Boolean    add_label_to_comment = FALSE;
  CharPtr    new_comment;
  Int4       comment_len = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return FALSE;

  FeatDefLabel (sfp, label, sizeof (label), OM_LABEL_CONTENT);

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = RnaRefFromLabel (featdef_to, label, &add_label_to_comment);
  
  if (add_label_to_comment)
  {
    comment_len += StringLen (label) + 2;
  }
  if (grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && StringCmp (grp->locus, label) != 0)
    {
      comment_len += StringLen (grp->locus) + 2;
    }
    if (!StringHasNoText (grp->allele) && StringCmp (grp->allele, label) != 0)
    {
      comment_len += StringLen (grp->allele) + 2;
    }
    if (!StringHasNoText (grp->desc) && StringCmp (grp->desc, label) != 0)
    {
      comment_len += StringLen (grp->desc) + 2;
    }
    if (!StringHasNoText (grp->maploc) && StringCmp (grp->maploc, label) != 0)
    {
      comment_len += StringLen (grp->maploc) + 2;
    }
    if (!StringHasNoText (grp->locus_tag) && StringCmp (grp->locus_tag, label) != 0)
    {
      comment_len += StringLen (grp->locus_tag) + 2;
    }
  }
  if (comment_len > 0)
  {
    comment_len += StringLen (sfp->comment) + 3;
    new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
    if (!StringHasNoText (sfp->comment))
    {
      StringCat (new_comment, sfp->comment);
      StringCat (new_comment, "; ");
    }
    if (add_label_to_comment)
    {
      StringCat (new_comment, label);
      StringCat (new_comment, "; ");
    }
    /* append unused gene qualifiers */
    if (grp != NULL)
    {
      if (!StringHasNoText (grp->locus) && StringCmp (grp->locus, label) != 0)
      {
        StringCat (new_comment, grp->locus);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->allele) && StringCmp (grp->allele, label) != 0)
      {
        StringCat (new_comment, grp->allele);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->desc) && StringCmp (grp->desc, label) != 0)
      {
        StringCat (new_comment, grp->desc);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->maploc) && StringCmp (grp->maploc, label) != 0)
      {
        StringCat (new_comment, grp->maploc);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->locus_tag) && StringCmp (grp->locus_tag, label) != 0)
      {
        StringCat (new_comment, grp->locus_tag);
        StringCat (new_comment, "; ");
      }
    }
    /* remove last semicolon */
    new_comment[StringLen (new_comment) - 2] = 0;
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;

  }
  
  /* free generef */
  grp = GeneRefFree (grp);      
  
  return TRUE;
}


/* These functions are used for converting features on nucleotide sequences to
 * features on protein sequences */

/* copied from seqport.c, for the benefit of load_fuzz_to_DNA */
static Boolean add_fuzziness_to_loc (SeqLocPtr slp, Boolean less)
{
	IntFuzzPtr ifp;
	SeqIntPtr sint;
	SeqPntPtr spnt;	

	sint = NULL;
	spnt = NULL;

	if(slp->choice == SEQLOC_INT)
		sint = (SeqIntPtr) slp->data.ptrvalue;
	else
	{
		if(slp->choice == SEQLOC_PNT)
			spnt = (SeqPntPtr) slp->data.ptrvalue;
		else
			return FALSE;
	}
	ifp = IntFuzzNew();
	ifp->choice = 4;
	ifp->a = less ? 2 : 1;

	if(spnt != NULL)
		spnt->fuzz = ifp;
	else
	{
		if(less)
			sint->if_from = ifp;
		else
			sint->if_to = ifp;
	}

	return TRUE;
}


/* copied from seqport.c, for the benefit of MYdnaLoc_to_aaLoc */
static Boolean load_fuzz_to_DNA(SeqLocPtr dnaLoc, SeqLocPtr aaLoc, Boolean 
first)
{
	Uint1 strand;
	SeqPntPtr spnt;
	SeqIntPtr sint;
	IntFuzzPtr ifp;
	Boolean load, less;

	load = FALSE;
	strand = SeqLocStrand(aaLoc);
	if(aaLoc->choice == SEQLOC_INT)
	{
		sint = (SeqIntPtr) aaLoc->data.ptrvalue;
		if((first && strand != Seq_strand_minus ) || 
			(!first && strand == Seq_strand_minus))	/*the first 
Seq-loc*/
		{
			ifp = sint->if_from;
			if(ifp && ifp->choice == 4 )
				load = (ifp->a == 2);
		}
		else
		{
			ifp = sint->if_to;
			if(ifp && ifp->choice == 4)
				load = (ifp->a == 1);
		}
	}
	else if(aaLoc->choice == SEQLOC_PNT)
	{
		spnt = (SeqPntPtr) aaLoc->data.ptrvalue;
		ifp = spnt->fuzz;
		if(ifp && ifp->choice == 4)
		{
			if(first)
				load = (ifp->a == 2);
			else
				load = (ifp->a == 1);
		}
	}

	if(load)
	{
		if(SeqLocStrand(dnaLoc) == Seq_strand_minus)
			less = (first == FALSE);
		else
			less = first;
		add_fuzziness_to_loc (dnaLoc, less);
		return TRUE;
	}
	else
		return FALSE;
}	


static SeqLocPtr MYdnaLoc_to_aaLoc(SeqFeatPtr sfp, 
                                   SeqLocPtr location_loc, 
                                   Boolean merge, 
                                   Int4Ptr frame, 
                                   Boolean allowTerminator)
{
	SeqLocPtr aa_loc = NULL, loc;
	CdRegionPtr crp;
	Int4 product_len, end_pos, frame_offset;
	GatherRange gr;
	Int4 a_left = 0, a_right, last_aa = -20, aa_from, aa_to;
	SeqLocPtr slp, slp1, slp2;
	Int2 cmpval;
	SeqIdPtr aa_sip;
	BioseqPtr bsp;

	if ((sfp == NULL) || (location_loc == NULL)) return aa_loc;
	if (sfp->data.choice != 3) return aa_loc;
	if (sfp->product == NULL) return aa_loc;

	crp = (CdRegionPtr) sfp->data.value.ptrvalue;
	if(crp == NULL) return aa_loc;

    /* each interval of location_loc must be equal to or contained in
     * an interval of sfp->location
     */
    slp1 = SeqLocFindNext (sfp->location, NULL);
    slp2 = SeqLocFindNext (location_loc, NULL);
    while (slp2 != NULL && slp1 != NULL) {
      cmpval = SeqLocCompare (slp2, slp1);
      if (cmpval == SLC_A_IN_B || cmpval == SLC_A_EQ_B) {
        slp2 = SeqLocFindNext (location_loc, slp2);
      } else {
        slp1 = SeqLocFindNext (sfp->location, slp1);
      }
    }
    if (slp1 == NULL) return aa_loc;
      
	aa_sip = SeqLocId(sfp->product);
	if (aa_sip == NULL) return aa_loc;
	bsp = BioseqLockById(aa_sip);
	if (bsp == NULL) return aa_loc;
	end_pos = bsp->length - 1;
	BioseqUnlock(bsp);

	if(crp->frame == 0)
		frame_offset = 0;
	else
		frame_offset = (Int4)crp->frame-1;

	slp = NULL;
	product_len = 0;
	loc = NULL;
	while ((slp = SeqLocFindNext(sfp->location, slp))!=NULL)
	{
	   if (SeqLocOffset(location_loc, slp, &gr, 0))
	   {
			SeqLocOffset(slp, location_loc, &gr, 0);
		
			a_left = gr.left + product_len - frame_offset;
			a_right = gr.right + product_len - frame_offset;

			aa_from = a_left / 3;
			aa_to = a_right / 3;

			if (aa_from < 0)
				aa_from = 0;
			if (aa_to > end_pos)
				aa_to = end_pos;

			if (merge)
			{
				if (aa_from <= last_aa)  /* overlap due to 
codons */
					aa_from = last_aa+1;  /* set up to merge 
*/
			}

			if (aa_from <= aa_to || (allowTerminator && aa_from == aa_to + 1))
			{
				if(loc != NULL)
				{
					if(aa_loc == NULL)
						load_fuzz_to_DNA(loc, location_loc, TRUE);
					SeqLocAdd(&aa_loc, loc, merge, FALSE);
				}
				loc = SeqLocIntNew(aa_from, aa_to, 0, aa_sip);
				last_aa = aa_to;
			}
	     }

	     product_len += SeqLocLen(slp);		
	}

	if(loc != NULL)
	{
		if(aa_loc == NULL)
			load_fuzz_to_DNA(loc, location_loc, TRUE);
		load_fuzz_to_DNA(loc, location_loc, FALSE);
		SeqLocAdd(&aa_loc, loc, merge, FALSE);
	}
	if (frame != NULL)
	    *frame = a_left % 3;

	return SeqLocPackage(aa_loc);
}


NLM_EXTERN SeqLocPtr BuildProtLoc (SeqFeatPtr overlapping_cds, SeqLocPtr slp, Int4Ptr frame)
{
  SeqLocPtr tmp_loc, aa_loc, prot_loc = NULL, last_loc = NULL, next_loc;
  Boolean   partial5, partial3;
  BioseqPtr prot_bsp;
  Boolean   is_ordered = FALSE;
  Boolean   first = TRUE;
  
  prot_bsp = BioseqFindFromSeqLoc (overlapping_cds->product);
  if (prot_bsp == NULL) {
    return NULL;
  }
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  tmp_loc = SeqLocFindNext (slp, NULL);
  while (tmp_loc != NULL) {
    if (tmp_loc->choice == SEQLOC_NULL) {
      is_ordered = TRUE;
    } else {
      if (first) {
        aa_loc = MYdnaLoc_to_aaLoc (overlapping_cds, tmp_loc, FALSE, frame, FALSE);
        first = FALSE;
      } else {
        aa_loc = MYdnaLoc_to_aaLoc (overlapping_cds, tmp_loc, FALSE, NULL, FALSE);
      }
    }
    if (last_loc == NULL) {
      prot_loc = aa_loc;
    } else {
      last_loc->next = aa_loc;
    }
    last_loc = aa_loc;
    tmp_loc = SeqLocFindNext (slp, tmp_loc);
  }
  if (prot_loc != NULL && prot_loc->next != NULL) {
    tmp_loc = NULL;
    for (aa_loc = prot_loc; aa_loc != NULL; aa_loc = next_loc) {
      next_loc = aa_loc->next;
      aa_loc->next = NULL;
      
      last_loc = SeqLocMerge (prot_bsp, tmp_loc, aa_loc, FALSE, TRUE, is_ordered);
      tmp_loc = SeqLocFree (tmp_loc);
      aa_loc = SeqLocFree (aa_loc);
      tmp_loc = last_loc;
      last_loc = NULL;

      aa_loc = next_loc;
    }
    prot_loc = tmp_loc;
  }
  SetSeqLocPartial (prot_loc, partial5, partial3);
  return prot_loc;
}

NLM_EXTERN Boolean ConvertRegionToProtFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  BioseqPtr  bsp;
  ProtRefPtr prp;
  SeqFeatPtr cds;
  SeqLocPtr  location;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION)
  {
    return FALSE; 
  }
  
  /* only convert features that are on protein sequences */
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (!ISA_aa (bsp->mol))
  {
    cds = SeqMgrGetOverlappingCDS (sfp->location, NULL);
    if (cds == NULL)
    {
      return FALSE;
    }
    else
    {
      location = BuildProtLoc (cds, sfp->location, NULL);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = location;
    }
  }
  
  prp = ProtRefNew ();
  if (prp != NULL)
  {
    prp->name = ValNodeNew(NULL);
    if (prp->name != NULL)
    {
      /* use region name for protein name */
      prp->name->data.ptrvalue = sfp->data.value.ptrvalue;
      switch (featdef_to) 
      {
        case FEATDEF_PROT :
          prp->processed = 0;
          break;
        case FEATDEF_preprotein :
          prp->processed = 1;
          break;
        case FEATDEF_mat_peptide_aa :
          prp->processed = 2;
          break;
        case FEATDEF_sig_peptide_aa :
          prp->processed = 3;
          break;
        case FEATDEF_transit_peptide_aa :
          prp->processed = 4;
          break;
        default :
          break;
      }

      sfp->data.value.ptrvalue = prp;
      sfp->data.choice = SEQFEAT_PROT;
    }
  }  
  return TRUE;
}


NLM_EXTERN SeqLocPtr GetProteinLocationForNucleotideFeatureConversion (SeqLocPtr nuc_slp, BoolPtr no_cds)
{
  SeqFeatPtr cds;
  SeqMgrFeatContext cds_context;
  SeqLocPtr  prot_slp;

  cds = SeqMgrGetOverlappingCDS (nuc_slp, &cds_context);
  if (cds == NULL) {
    if (no_cds != NULL) {
      *no_cds = TRUE;
    }
    return NULL;
  } else if (no_cds != NULL) {
    *no_cds = FALSE;
  }

  prot_slp = BuildProtLoc (cds, nuc_slp, NULL);
  return prot_slp;
}



/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertImpToProt () - Convert a given import feature to a    */
/*                           protein feature.                          */
/*                                                                     */
/*    Note : Any of the Import feature's gbquals that can be converted */
/*           to protein fields are caught in the automatic cleanup     */
/*           called during reindexing, so they don't need to be        */
/*           converted here.                                           */
/*                                                                     */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertImpToProtFunc 
(SeqFeatPtr  sfp,
 Uint2       featdef_to)
{
  ImpFeatPtr ifp;
  SeqFeatPtr cds;
  SeqLocPtr  slp;
  SeqFeatPtr newSfp;
  Int4       frame;
  ProtRefPtr prp;
  SeqIdPtr   sip;
  BioseqPtr  bsp;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP)
  {
    return FALSE;
  }
  /* Get the Import Feature */

  ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
  if (NULL == ifp)
  {
    return FALSE;
  }

  /* Convert the location to a protein location */
  cds = SeqMgrGetOverlappingCDS (sfp->location, &fcontext);
  if (cds == NULL)
  {
    return FALSE;
  }

  slp = BuildProtLoc (cds, sfp->location, &frame);
  if (slp == NULL)
  {
    return FALSE;
  }

  /* Create a new generic feature */

  sip = SeqLocId (cds->product);
  if (sip == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  bsp = BioseqLockById (sip);
  if (bsp == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  newSfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_PROT, slp);
  BioseqUnlock (bsp);
  if (newSfp == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  /* Make it into a protein feature */

  prp = ProtRefNew ();
  newSfp->data.value.ptrvalue = (Pointer) prp;
  if (NULL == prp)
  {
    slp = SeqLocFree (slp);
    newSfp = SeqFeatFree (newSfp);
    return FALSE;
  }

  switch (featdef_to) {
    case FEATDEF_mat_peptide_aa :
      prp->processed = 2;
      break;
    case FEATDEF_sig_peptide_aa :
      prp->processed = 3;
      break;
    case FEATDEF_transit_peptide_aa :
      prp->processed = 4;
      break;
  }

  /* Transfer unchanged fields from old feature */

  newSfp->partial     = sfp->partial;
  newSfp->excpt       = sfp->excpt;
  newSfp->exp_ev      = sfp->exp_ev;
  newSfp->pseudo      = sfp->pseudo;
  newSfp->comment     = sfp->comment;
  newSfp->qual        = sfp->qual;
  newSfp->title       = sfp->title;
  newSfp->ext         = sfp->ext;
  newSfp->cit         = sfp->cit;
  newSfp->xref        = sfp->xref;
  newSfp->dbxref      = sfp->dbxref;
  newSfp->except_text = sfp->except_text;

  /* Null out pointers to transferred fields from old feature  */
  /* so that they don't get deleted when the feature does,     */

  sfp->comment     = NULL;
  sfp->qual        = NULL;
  sfp->title       = NULL;
  sfp->ext         = NULL;
  sfp->cit         = NULL;
  sfp->xref        = NULL;
  sfp->dbxref      = NULL;
  sfp->except_text = NULL;

  /* Mark the old feature for deletion */

  sfp->idx.deleteme = TRUE;
  return TRUE;
}


NLM_EXTERN SeqLocPtr FindNucleotideLocationForProteinFeatureConversion (SeqLocPtr slp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr cds;
  SeqLocPtr  slp_nuc;
  Boolean    partial5, partial3;

  cds = SeqMgrGetCDSgivenProduct (BioseqFindFromSeqLoc (slp), &context);
  if (NULL == cds)
  {
    return NULL;
  }

  slp_nuc = aaLoc_to_dnaLoc (cds, slp);
  if (slp_nuc == NULL) 
  {
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (partial5 && partial3) {
      slp_nuc = AsnIoMemCopy (cds->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
    }
  }
  return slp_nuc;
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertProtToImp () -                                        */
/*                                                                     */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertProtToImpFunc (SeqFeatPtr  sfp, Uint2 featdef_to)
{
  ProtRefPtr    prp;
  SeqLocPtr     slp;
  ImpFeatPtr    ifp;
  CharPtr       name;
  CharPtr       ec;
  CharPtr       activity;
  ValNodePtr    vnp;
  GBQualPtr     gbqual;
  GBQualPtr     prevGbq;
  GBQualPtr     topOfGbqList;
  DbtagPtr      dbt;
  Char          idStr[64];
  ObjectIdPtr   oip;
  Uint2         entityID;
  BioseqPtr     bsp;
  SeqAnnotPtr   old_sap, new_sap;
  SeqFeatPtr    tmp, tmp_prev;

  /* Make sure that we have a matching peptide feature */

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT)
  {
    return FALSE;
  }
  entityID = sfp->idx.entityID;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
  {
    return FALSE;
  }

  switch (sfp->idx.subtype) {
    case FEATDEF_mat_peptide_aa :
      if (2 != prp->processed)
        return FALSE;
      break;
    case FEATDEF_sig_peptide_aa :
      if (3 != prp->processed)
        return FALSE;
      break;
    case FEATDEF_transit_peptide_aa :
      if (4 != prp->processed)
        return FALSE;
      break;
  }

  /* Convert the location from the protein */
  /* to the nucleotide Bioseq.             */

  slp = FindNucleotideLocationForProteinFeatureConversion (sfp->location);
  if (NULL == slp)
    return FALSE;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  /* move feature to correct annot */
  if (sfp->idx.parenttype == OBJ_SEQANNOT
      && (old_sap = (SeqAnnotPtr) sfp->idx.parentptr) != NULL
      && old_sap->type == 1
      && (bsp = BioseqFindFromSeqLoc (sfp->location)) != NULL)
  {
    tmp_prev = NULL;
    tmp = old_sap->data;
    while (tmp != sfp) {
      tmp_prev = tmp;
      tmp = tmp->next;
    }
    if (tmp != NULL) {
      if (tmp_prev == NULL) {
        old_sap->data = tmp->next;
      } else {
        tmp_prev->next = tmp->next;
      }
      if (old_sap->data == NULL) {
        old_sap->idx.deleteme = TRUE;
      }
      new_sap = bsp->annot;
      while (new_sap != NULL && new_sap->type != 1) {
        new_sap = new_sap->next;
      }
      if (new_sap == NULL) {
        new_sap = SeqAnnotNew ();
        new_sap->type = 1;
        new_sap->next = bsp->annot;
        bsp->annot = new_sap;
      }
      sfp->next = new_sap->data;
      new_sap->data = sfp;
      sfp->idx.parentptr = new_sap;
    }
  }
        

  /* Create a new import feature and */
  /* attach it to the feature.       */

  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  /* set key */
  ifp->key = StringSave (GetFeatureNameFromFeatureType(featdef_to));

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  /* Store the protein fields as  */
  /* gbqual qualifier/value pairs */

  name = NULL;
  vnp = prp->name;
  if (vnp != NULL)
    name = vnp->data.ptrvalue;
  if (name == NULL) 
    name = prp->desc;

  if (name != NULL) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    topOfGbqList = gbqual;
    gbqual->qual = StringSave ("product");
    gbqual->val = StringSave (name);
  }

  prevGbq = gbqual;

  ec = NULL;
  vnp = prp->ec;
  if (vnp != NULL)
    ec = (CharPtr) vnp->data.ptrvalue;
  
  if (ec != NULL) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    prevGbq->next = gbqual;
    gbqual->qual = StringSave ("EC_number");
    gbqual->val = StringSave (ec);
  }

  prevGbq = gbqual;

  activity = NULL;
  vnp = prp->activity;
  if (vnp != NULL)
    activity = (CharPtr) vnp->data.ptrvalue;
  
  if (NULL != activity) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    prevGbq->next = gbqual;
    gbqual->qual = StringSave ("function");
    gbqual->val = StringSave (activity);
  }

  prevGbq = gbqual;

  for (vnp = prp->db; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (NULL == dbt ) 
      continue;
    if (! StringHasNoText (dbt->db)) {
      gbqual = GBQualNew ();
      if (NULL == gbqual)
        continue;
      prevGbq->next = gbqual;
      oip = dbt->tag;
      if (oip->str != NULL && (! StringHasNoText (oip->str))) {
        sprintf (idStr, "%s:%s", (CharPtr)dbt->tag, oip->str);
        gbqual->qual = StringSave ("db_xref");
        gbqual->val = StringSave (idStr);
      } else {
        sprintf (idStr, "%s:%ld", (CharPtr)dbt->tag, (long) oip->id);
        gbqual->qual = StringSave ("db_xref");
        gbqual->val = StringSave (idStr);
      }
      prevGbq = gbqual;
    }
  }

  /* Insert the new qualifiers in front of any existing ones */

  gbqual->next = sfp->qual;
  sfp->qual = topOfGbqList;

  /* Free the obsolete Protein reference */

  ProtRefFree (prp);
  return TRUE;
}


/* functions for converting from biosource */
NLM_EXTERN CharPtr SubSourceText (BioSourcePtr biop, Uint1 subtype, BoolPtr found)
{
  Int4 subtype_len = 0;
  SubSourcePtr ssp;
  CharPtr subtype_txt = NULL;
  
  if (biop == NULL || biop->subtype == NULL) return NULL;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == subtype) {
      if (found != NULL) *found = TRUE;
      if (!StringHasNoText (ssp->name)) {
        subtype_len += StringLen (ssp->name) + 1;
      }
    }
  }
  if (subtype_len == 0) return NULL;
  subtype_txt = (CharPtr) MemNew (sizeof (Char) * subtype_len);
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == subtype && !StringHasNoText (ssp->name)) {
      if (!StringHasNoText (subtype_txt)) {
        StringCat (subtype_txt, ";");
      }
      StringCat (subtype_txt, ssp->name);
    }
  }
  return subtype_txt;
}

NLM_EXTERN CharPtr OrgModText (BioSourcePtr biop, Uint1 subtype, BoolPtr found)
{
  Int4 subtype_len = 0;
  OrgModPtr omp;
  CharPtr subtype_txt = NULL;
  
  if (biop == NULL
     || biop->org == NULL 
     || biop->org->orgname == NULL 
     || biop->org->orgname->mod == NULL) {
    return NULL;
  }
     
  for (omp = biop->org->orgname->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == subtype) {
      if (found != NULL) *found = TRUE;
      if (!StringHasNoText (omp->subname)) {
        subtype_len += StringLen (omp->subname) + 1;
      }
    }
  }
  if (subtype_len == 0) return NULL;
  subtype_txt = (CharPtr) MemNew (sizeof (Char) * subtype_len);
  for (omp = biop->org->orgname->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == subtype && !StringHasNoText (omp->subname)) {
      if (!StringHasNoText (subtype_txt)) {
        StringCat (subtype_txt, ";");
      }
      StringCat (subtype_txt, omp->subname);
    }
  }
  return subtype_txt;
}

NLM_EXTERN CharPtr NoteText (BioSourcePtr biop, CharPtr comment)
{
  CharPtr orgmod_note, subsource_note;
  Int4    text_len = 0;
  CharPtr note_text = NULL;
  
  orgmod_note = OrgModText(biop, ORGMOD_other, NULL);
  if (!StringHasNoText (orgmod_note)) {
    text_len += StringLen (orgmod_note) + 1;
  }
  subsource_note = SubSourceText (biop, SUBSRC_other, NULL);
  if (!StringHasNoText (subsource_note)) {
    text_len += StringLen (subsource_note) + 1;
  }
  if (!StringHasNoText (comment)) {
    text_len += StringLen (comment) + 1;
  }
  
  if (text_len == 0) return NULL;  
  
  note_text = (CharPtr) MemNew (sizeof(Char) * text_len);
  if (!StringHasNoText (orgmod_note)) {
    StringCat (note_text, orgmod_note);
  }
  orgmod_note = MemFree (orgmod_note);
  if (!StringHasNoText (subsource_note)) {
    if (!StringHasNoText (note_text)) {
      StringCat (note_text, ";");
    }
    StringCat (note_text, subsource_note);
  }
  subsource_note = MemFree (subsource_note);
  
  if (!StringHasNoText (comment)) {
    if (!StringHasNoText (note_text)) {
      StringCat (note_text, ";");
    }
    StringCat (note_text, comment);
  }
  return note_text;
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertBioSrcToRepeatRegion ()                                  */
/*                                                                     */
/* 9/28/2004: Changed to convert all BioSource features with notes     */
/* instead of ones with transposon or insertion_seq qualifiers.        */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertBioSrcToRepeatRegion (SeqFeatPtr sfp, Uint2 featdef_to)
{
  BioSourcePtr  biop;
  GBQualPtr     gbqual;
  ImpFeatPtr    ifp;
  CharPtr       transposon_txt, insertion_seq_txt, note_txt;
  Boolean       is_transposon = FALSE, is_insertion_seq = FALSE;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return FALSE;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  
  transposon_txt = SubSourceText (biop, SUBSRC_transposon_name, &is_transposon);
  insertion_seq_txt = SubSourceText (biop, SUBSRC_insertion_seq_name, &is_insertion_seq);
  note_txt = NoteText (biop, sfp->comment);
  
  
  /* Create a new Import Feature */

  ifp = ImpFeatNew ();
  if (NULL == ifp)
    return FALSE;
  ifp->key = StringSave ("repeat_region");

  /* Copy relevant info from the BioSource */
  /* feature to the Import feature.        */

    
  /* Delete the old BioSource feature */

  sfp->data.value.ptrvalue = BioSourceFree (biop);

  /* Attach the new Import feature in its place */

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = ifp;
  
  if (is_transposon) {
    gbqual = GBQualNew ();
    gbqual->qual = StringSave ("mobile_element");
    gbqual->val = (CharPtr) MemNew (sizeof(Char) * (StringLen (transposon_txt) + 12));
    StringCat (gbqual->val, "transposon:");
    StringCat (gbqual->val, transposon_txt);
    gbqual->next = sfp->qual;
    sfp->qual = gbqual;
  }
  transposon_txt = MemFree (transposon_txt); 

  if (is_insertion_seq) {
    gbqual = GBQualNew ();
    gbqual->qual = StringSave ("mobile_element");
    gbqual->val = (CharPtr) MemNew (sizeof(Char) * (StringLen (insertion_seq_txt) + 19));
    StringCat (gbqual->val, "insertion sequence:");
    StringCat (gbqual->val, insertion_seq_txt);
    gbqual->next = sfp->qual;
    sfp->qual = gbqual;
  }
  insertion_seq_txt = MemFree (insertion_seq_txt); 
  
  sfp->comment = MemFree (sfp->comment);
  sfp->comment = note_txt;
  return TRUE;
}


NLM_EXTERN Boolean 
ConvertNonPseudoCDSToMiscFeat 
(SeqFeatPtr sfp,
 Boolean viral)
{
  CdRegionPtr          cdrp;
  ImpFeatPtr           ifp;
  CharPtr              noteStr;
  BioseqPtr            protBsp;
  SeqMgrFeatContext    protContext;
  CharPtr              protName = NULL;
  SeqFeatPtr           protSfp;
  ProtRefPtr           prp;
  ValNodePtr           vnp;
  Int4                 note_len = 0;
  CharPtr              viral_fmt = "nonfunctional %s due to mutation";
  CharPtr              similar_fmt = "similar to %s";

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->product == NULL)
  {
    return FALSE;
  }
      
  /* Get the CD region part of the feature, and */
  /* the associated protein bioseq.             */

  cdrp = (CdRegionPtr) sfp->data.value.ptrvalue;
  protBsp = BioseqFindFromSeqLoc (sfp->product);

  if (protBsp == NULL) return FALSE;

  /* Convert the CDS feature to a misc_feat */

  CdRegionFree (cdrp);
  sfp->data.value.ptrvalue = NULL;

  ifp = ImpFeatNew ();
  if (NULL == ifp) return FALSE;
  ifp->key = StringSave ("misc_feature");

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  /* Add a name key to the misc_feature */

  protSfp = SeqMgrGetBestProteinFeature (protBsp, &protContext);
  if (protSfp != NULL)
  {
    prp = (ProtRefPtr) protSfp->data.value.ptrvalue;

    if (prp != NULL) 
    {
      note_len = StringLen (sfp->comment) + StringLen (prp->desc) + 5;
          
      vnp = prp->name;
      if (NULL != vnp)
      {
        protName = (CharPtr) vnp->data.ptrvalue;
        if (NULL != protName) 
        {
          if (viral) {
            note_len += StringLen (viral_fmt) + StringLen (protName);
          } else {   
            note_len += StringLen (similar_fmt) + StringLen (protName);
          }
        }
      }  
      noteStr = (CharPtr) MemNew (sizeof (Char) * note_len);
        
      if (NULL != protName) {
        if (viral) {
          sprintf (noteStr, viral_fmt, protName);
        } else {
          sprintf (noteStr, similar_fmt, protName);
        }
      }
      if (!StringHasNoText (prp->desc)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, prp->desc);
      }
      if (!StringHasNoText (sfp->comment)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, sfp->comment);
      }
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = noteStr;
    }
  }

  /* set the subtype to zero so that it will be reindexed */
  sfp->idx.subtype = 0;
  return TRUE;
}


NLM_EXTERN Uint1 RnaTypeFromFeatdef (Uint2 featdef)
{
  switch (featdef) 
  {
    case FEATDEF_preRNA:
      return 1;
      break;
    case FEATDEF_mRNA:
      return 2;
      break;
    case FEATDEF_tRNA:
      return 3;
      break;
    case FEATDEF_rRNA:
      return 4;
      break;
    case FEATDEF_snRNA:
      return 8;
      break;
    case FEATDEF_scRNA:
      return 8;
      break;
    case FEATDEF_snoRNA:
      return 8;
      break;
    case FEATDEF_ncRNA:
      return 8;
      break;
    case FEATDEF_tmRNA:
      return 9; 
      break;
    case FEATDEF_otherRNA:
    default:
      return 255;
      break;
  }
}


NLM_EXTERN Boolean ConvertRegionToRNAFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  RnaRefPtr  rrp;
  CharPtr    str, new_comment;
  Boolean    add_to_comment = FALSE;
  Int4       len;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION)
  {
    return FALSE;
  }

  str = (CharPtr) sfp->data.value.ptrvalue;
  rrp = RnaRefFromLabel (featdef_to, str, &add_to_comment);

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = (Pointer) rrp;

  if (add_to_comment) {
    if (sfp->comment == NULL) {
      sfp->comment = str;
      str = NULL;
    } else {      
      len = StringLen (sfp->comment) + StringLen (str) + 3;
      new_comment = MemNew (sizeof (Char) * len);
      sprintf (new_comment, "%s; %s", sfp->comment, str);
      sfp->comment = MemFree (sfp->comment);
      str = MemFree (str);
      sfp->comment = new_comment;
    }
  } else {
    str = MemFree (str);
  }
  return TRUE;
  
  return TRUE;
}


NLM_EXTERN CharPtr GetImportFeatureName (Uint2 featdef_key)
{
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) 
  {
    if (featdef_key == key)
    {
      return curr->typelabel;
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  return NULL;
}


NLM_EXTERN Boolean ConvertRegionToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  GBQualPtr          gbqual;
  ImpFeatPtr         ifp;
  CharPtr            str;
  CharPtr            featname_to;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION) return FALSE;
  featname_to = GetImportFeatureName (featdef_to);
  ifp = ImpFeatNew ();
  if (NULL == ifp)
    return FALSE;

  str = (CharPtr) sfp->data.value.ptrvalue;
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  if (featname_to == NULL)
  {
    ifp->key = StringSave ("misc_feature");
  }
  else
  {
    ifp->key = StringSave (featname_to);
  }

  if (! StringHasNoText (str)) {
    gbqual = GBQualNew ();
    if (gbqual != NULL) {
      gbqual->qual = StringSave ("note");
      gbqual->val = str;
      gbqual->next = sfp->qual;
      sfp->qual = gbqual;
    }
  }
  return TRUE;
}


NLM_EXTERN Boolean ConvertImpToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  ImpFeatPtr         ifp;
  CharPtr            featname;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP) return FALSE;

  ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
  if (NULL == ifp)
    return FALSE;
  
  featname = GetImportFeatureName (featdef_to);
  ifp->key = MemFree (ifp->key);
  if (featname == NULL)
  {
    ifp->key = StringSave ("misc_feature");
  }
  else
  {
    ifp->key = StringSave (featname);
  }

  return TRUE;
}


NLM_EXTERN Boolean ConvertGeneToMiscFeatFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  ImpFeatPtr  ifp;
  CharPtr     new_comment;
  GeneRefPtr  grp;
  Int4        comment_len = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE)
  {
    return FALSE;
  }
  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp != NULL) 
  {
    if (!StringHasNoText (grp->locus))
    {
      comment_len += StringLen (grp->locus) + 2;
    }
    if (!StringHasNoText (grp->desc))
    {
      comment_len += StringLen (grp->desc) + 2;
    }
  }
  if (comment_len == 0) 
  {
    /* nothing to add to comment */
  }
  else
  {
    /* add one for terminating NULL */
    comment_len++;
    if (!StringHasNoText (sfp->comment))
    {
      comment_len += StringLen (sfp->comment) + 2;
    }

    new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
    /* NOTE - I don't have to check for grp == NULL because
     * comment_len would only have been > 0 if grp had existed
     * and had nonempty fields.
     */
    if (!StringHasNoText (grp->desc))
    {
      StringCat (new_comment, grp->desc);
      StringCat (new_comment, "; ");    
    }
    if (!StringHasNoText (grp->locus))
    {
      StringCat (new_comment, grp->locus);
      StringCat (new_comment, "; ");    
    }
    if (!StringHasNoText (sfp->comment))
    {
      StringCat (new_comment, sfp->comment);
      StringCat (new_comment, "; ");    
    }
    /* remove last semicolon */
    new_comment[StringLen (new_comment) - 2] = 0;
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;
  }

  sfp->data.value.ptrvalue =
    GeneRefFree ((GeneRefPtr) sfp->data.value.ptrvalue);
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");
  return TRUE;
}




/* For mat-peptide instantiation */
static SeqIdPtr MakeMatPeptideProductId (SeqLocPtr mat_peptide_loc)
{
  ObjectIdPtr oip;
  SeqIdPtr    sip;
  Char        id_str[500];
  Int4        len;
  BioseqPtr   bsp;

  if (mat_peptide_loc == NULL) return NULL;

  bsp = BioseqFindFromSeqLoc (mat_peptide_loc);
  if (bsp == NULL) return NULL;
  sip = bsp->id;
  while (sip != NULL && sip->choice != SEQID_OTHER) {
    sip = sip->next;
  }
  if (sip == NULL) {
    sip = bsp->id;
    while (sip != NULL && sip->choice != SEQID_GENBANK) {
      sip = sip->next;
    }
  }
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, 0);
  }
  SeqIdWrite (sip, id_str, PRINTID_TEXTID_ACC_ONLY, sizeof (id_str) - 41);
  len = StringLen (id_str);
  sprintf (id_str + len, ":%d-%d", SeqLocStart (mat_peptide_loc) + 1, SeqLocStop (mat_peptide_loc) + 1);
  
  oip = ObjectIdNew ();
  oip->str = StringSave (id_str);
  sip = ValNodeNew (NULL);
  sip->choice = SEQID_LOCAL;
  sip->data.ptrvalue = oip;  
  return sip;
}


static void InstantiateMatPeptideProductForProteinFeature (SeqFeatPtr sfp, Pointer data)
{
  BioseqPtr mat_bsp, prot_bsp;
  SeqEntryPtr master, sep, old;
  Int4        i, len;
  Int2        residue;
  SeqLocPtr   slp;
  ByteStorePtr src, dst;
  ProtRefPtr   prp_orig, prp_mat;
  SeqFeatPtr   sfp_mat;
  SeqDescrPtr  sdp;
  MolInfoPtr   mip;
  Boolean      partial5, partial3;
  Char         defline_buf[1024];

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mat_peptide_aa || sfp->product != NULL) {
    return;
  }

  prp_orig = sfp->data.value.ptrvalue;

  prot_bsp = BioseqFindFromSeqLoc (sfp->location);
  if (prot_bsp == NULL) {
    return;
  }
  master = GetBestTopParentForData (sfp->idx.entityID, prot_bsp);
  if (master == NULL) return;

  src = (ByteStorePtr) prot_bsp->seq_data;

  mat_bsp = BioseqNew ();
  if (mat_bsp == NULL) {
     return;
  }
  mat_bsp->mol = Seq_mol_aa;
  mat_bsp->repr = Seq_repr_raw;
  mat_bsp->seq_data_type = Seq_code_ncbieaa;
  mat_bsp->length = SeqLocLen (sfp->location);
  dst = BSNew (0);
  mat_bsp->seq_data = (SeqDataPtr) dst;
  BSSeek (dst, 0, SEEK_SET);

  for (slp = SeqLocFindNext (sfp->location, NULL);
       slp != NULL;
       slp = SeqLocFindNext (sfp->location, slp)) {
    BSSeek (src, SeqLocStart (slp), SEEK_SET);
    len = SeqLocLen (slp);
    for (i = 0; i < len; i++) {
      residue = BSGetByte (src);
      BSPutByte (dst, residue);
    }
  }

  old = SeqEntrySetScope (master);
  
  /*mat_bsp->id = MakeNewProteinSeqId (sfp->location, NULL); */
  mat_bsp->id = MakeMatPeptideProductId (sfp->location);
  SeqMgrAddToBioseqIndex (mat_bsp);
  SeqEntrySetScope (old);
  sep = SeqEntryNew ();
  if (sep != NULL) {
    sep->choice = 1;
    sep->data.ptrvalue = (Pointer) mat_bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mat_bsp, sep);
  }
  SetSeqFeatProduct (sfp, mat_bsp);
  if (sep != NULL) {
    AddSeqEntryToSeqEntry (master, sep, TRUE);
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  /* set molinfo for new protein sequence */
  sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
  if (sdp != NULL) {
    mip = MolInfoNew ();
    sdp->data.ptrvalue = (Pointer) mip;
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 13;
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      /*
      } else if (partial) {
        mip->completeness = 2;
      */
      } else {
        mip->completeness = 0;
      }
    }
  }


  /* create protein feature for mat_bsp */
  sfp_mat = CreateNewFeature (SeqMgrGetSeqEntryForData (mat_bsp), NULL,
                          SEQFEAT_PROT, NULL);

  sfp_mat->location = SeqLocIntNew (0, mat_bsp->length - 1, Seq_strand_plus, SeqIdDup (mat_bsp->id));
  SetSeqLocPartial (sfp_mat->location, partial5, partial3);
  if (partial5 || partial3) {
    sfp_mat->partial = TRUE;
  }

  prp_mat = AsnIoMemCopy (prp_orig, (AsnReadFunc) ProtRefAsnRead, (AsnWriteFunc) ProtRefAsnWrite);
  prp_mat->processed = 0;
  sfp_mat->data.value.ptrvalue = prp_mat;
  if (sfp->comment != NULL) {
    sfp_mat->comment = StringSave (sfp->comment);
  }

  /* add title */
  SeqMgrIndexFeatures (0, mat_bsp);
  if (NewCreateDefLineBuf (NULL, mat_bsp, defline_buf, sizeof (defline_buf), FALSE, FALSE)
      && !StringHasNoText (defline_buf)) {
    sdp = CreateNewDescriptor (sep, Seq_descr_title);
    sdp->data.ptrvalue = StringSave (defline_buf);
  }

}


NLM_EXTERN void InstantiateMatPeptideProducts (SeqEntryPtr sep)
{
  VisitFeaturesInSep (sep, NULL, InstantiateMatPeptideProductForProteinFeature);
}



static void ConvertLocalIdsToTSAIdsCallback (BioseqPtr bsp, Pointer data)
{
  Int4            gpid = 0;
  SeqIdPtr        sip_local = NULL;
  SeqEntryPtr     top_sep;
  SeqIdPtr        sip_new;
  DbtagPtr        dbtag;
  CharPtr         id_fmt = "gpid:%d";
  ObjectIdPtr     oip = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL || !IsTSA (bsp)) {
    return;
  }

  top_sep = (SeqEntryPtr) data;

  for (sip_local = bsp->id;
       sip_local != NULL && sip_local->choice != SEQID_LOCAL;
       sip_local = sip_local->next)
  {}
  if (sip_local == NULL) return;
  oip = sip_local->data.ptrvalue;
  if (oip == NULL) return;

  gpid = GetGenomeProjectID (bsp);

  if (gpid > 0) {
    dbtag = DbtagNew ();
    dbtag->db = MemNew (sizeof (Char) * (StringLen (id_fmt) + 15));
    sprintf (dbtag->db, id_fmt, gpid);
    dbtag->tag = ObjectIdNew ();
    if (oip->str == NULL) {
      dbtag->tag->id = oip->id;
    } else {
      dbtag->tag->str = StringSave (oip->str);
    }
    sip_new = ValNodeNew (NULL);
    sip_new->choice = SEQID_GENERAL;
    sip_new->data.ptrvalue = dbtag;
    ReplaceSeqIdWithSeqId (sip_local, sip_new, top_sep);
  }
}     
          
  
NLM_EXTERN void ConvertLocalIdsToTSAIds (SeqEntryPtr sep)
{
  VisitBioseqsInSep (sep, sep, ConvertLocalIdsToTSAIdsCallback);
}


NLM_EXTERN Int4 GetDeflinePosForFieldType (ValNodePtr field)
{
  Int4    i, rval = -1;
  CharPtr name;

  name = SummarizeFieldType (field);
  if (StringICmp (name, "specimen-voucher") == 0) {
    rval = DEFLINE_POS_Specimen_voucher;
  } else {
    for (i = 0; i < numDefLineModifiers; i++) {
      if (StringICmp (name, DefLineModifiers[i].name) == 0) {
        rval = i;
        break;
      }
    }
  }
  name = MemFree (name);
  return rval;
}


static void RemoveUnusedFieldTypes (FieldTypePtr PNTR orig_list)
{
  ValNodePtr vnp, prev = NULL, vnp_next;

  if (orig_list == NULL || *orig_list == NULL) {
    return;
  }
  for (vnp = *orig_list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (GetDeflinePosForFieldType (vnp) < 0) {
      if (prev == NULL) {
        *orig_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = FieldTypeFree (vnp);
    } else {
      prev = vnp;
    }
  }  
}


static Boolean RemoveMatchingFieldType (FieldTypePtr PNTR orig_list, FieldTypePtr match)
{
  ValNodePtr vnp, prev = NULL, vnp_next;
  Boolean    rval = FALSE;

  if (orig_list == NULL || *orig_list == NULL || match == NULL) {
    return rval;
  }

  for (vnp = *orig_list; vnp != NULL && !rval; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (CompareFieldTypes (vnp, match) == 0) {
      if (prev == NULL) {
        *orig_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = FieldTypeFree (vnp);
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean ListHasMatchingFieldType (FieldTypePtr list, FieldTypePtr match)
{
  Boolean rval = FALSE;

  if (list == NULL || match == NULL) return rval;

  while (list != NULL && !rval) {
    if (CompareFieldTypes (list, match) == 0) {
      rval = TRUE;
    } else {
      list = list->next;
    }
  }
  return rval;
}


static Int4 DefLineFieldTypeSortOrder [] = {
  Source_qual_strain,
  Source_qual_clone,
  Source_qual_isolate,
  Source_qual_haplotype,
  Source_qual_cultivar,
  Source_qual_specimen_voucher,
  Source_qual_ecotype,
  Source_qual_type,
  Source_qual_serotype,
  Source_qual_authority,
  Source_qual_breed
};


static int CompareFieldTypeByImportance (FieldTypePtr field1, FieldTypePtr field2)
{
  int rval = 0;
  Int4 index, num_defline_qual_sort_order;
  SourceQualChoicePtr scp1, scp2;

  if (field1 == NULL && field2 == NULL) {
    rval = 0;
  } else if (field1 == NULL) {
    rval = 1;
  } else if (field2 == NULL) {
    rval = -1;
  } else if (field1->choice == FieldType_source_qual && field2->choice == FieldType_source_qual) {
    scp1 = field1->data.ptrvalue;
    scp2 = field2->data.ptrvalue;
    num_defline_qual_sort_order = sizeof (DefLineFieldTypeSortOrder) / sizeof (Int4);
    for (index = 0; index < num_defline_qual_sort_order; index++)
    {
      if (scp1->data.intvalue == DefLineFieldTypeSortOrder [ index ]) return -1;
      if (scp2->data.intvalue == DefLineFieldTypeSortOrder [ index ]) return 1;
    }
    rval = CompareFieldTypes (field1, field2);
  } else {
    rval = CompareFieldTypes (field1, field2);
  }
  return rval;
}

static int LIBCALLBACK SortFieldTypeByImportance (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;

  return CompareFieldTypeByImportance (vnp1, vnp2);
}


typedef struct uniqbiosource {
  BioSourcePtr biop;
  ValNodePtr   available_fields;
  ValNodePtr   strings;
} UniqBioSourceData, PNTR UniqBioSourcePtr;

static Boolean AddQualToUniqBioSource (
  UniqBioSourcePtr u,
  FieldTypePtr     field
)
{
  CharPtr             val = NULL, tmp;
  SourceQualChoicePtr q;
  Boolean             rval = FALSE;

  if (u == NULL || field == NULL) return FALSE;
  if (field->choice != FieldType_source_qual) return FALSE;
  q = (SourceQualChoicePtr) field->data.ptrvalue;
  if (q == NULL || q->choice != SourceQualChoice_textqual) return FALSE;

  val = GetSourceQualFromBioSource (u->biop, q, NULL);
  if (StringHasNoText (val)) {
    val = MemFree (val);
  } else if (q->data.intvalue == Source_qual_specimen_voucher && StringNICmp (val, "personal:", 9) == 0) {
    tmp = StringSave (val + 9);
    val = MemFree (val);
    val = tmp;
  } else if (IsNonTextSourceQual (q->data.intvalue)) {
    val = MemFree (val);
    val = StringSave (GetSourceQualName (q->data.intvalue));
  }
  if (val != NULL) {
    ValNodeAddPointer (&(u->strings), 0, val);
    rval = TRUE;
  }
  RemoveMatchingFieldType (&(u->available_fields), field);

  u->available_fields = ValNodeSort (u->available_fields, SortFieldTypeByImportance);
  return rval;
}
 

static UniqBioSourcePtr UniqBioSourceNew (BioSourcePtr biop)
{
  UniqBioSourcePtr u;

  if (biop == NULL) return NULL;
  u = (UniqBioSourcePtr) MemNew (sizeof (UniqBioSourceData));
  u->biop = biop;
  u->available_fields = GetSourceQualFieldListFromBioSource (u->biop);
  u->strings = NULL;

  /* add tax name as first string */
  AddQualToUniqBioSource (u, u->available_fields);
  RemoveUnusedFieldTypes (&(u->available_fields));

  return u;
}

static UniqBioSourcePtr UniqBioSourceFree (UniqBioSourcePtr u)
{
  if (u != NULL) {
    u->available_fields = FieldTypeListFree (u->available_fields);
    u->strings = ValNodeFreeData (u->strings);
    u = MemFree (u);
  }
  return u;
}


static UniqBioSourcePtr UniqBioSourceCopy (UniqBioSourcePtr u)
{
  UniqBioSourcePtr u2;
  ValNodePtr       vnp;

  if (u == NULL) return NULL;
  u2 = (UniqBioSourcePtr) MemNew (sizeof (UniqBioSourceData));
  u2->biop = u->biop;
  u2->available_fields = FieldTypeListCopy (u->available_fields);
  for (vnp = u->strings; vnp != NULL; vnp = vnp->next) {
    ValNodeAddPointer (&(u2->strings), 0, StringSave (vnp->data.ptrvalue));
  }
  return u2;
}


/* The CompareOrganismDescriptors function compares the contents of the
 * lists of strings for each BioSrcDesc item.
 * The function returns:
 *     -1 if org1 < org2
 *      0 if org1 = org2
 *      1 if org1 > org2
 */
static int CompareUniqBioSource (
  UniqBioSourcePtr org1,
  UniqBioSourcePtr org2
)
{
  ValNodePtr vnp1, vnp2;
  int cmpval;

  vnp1 = org1->strings;
  vnp2 = org2->strings;

  while (vnp1 != NULL && vnp2 != NULL)
  {
    cmpval = StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    if (cmpval != 0) return cmpval;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return 0;
  }
  else if (vnp1 != NULL && vnp2 == NULL)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}


static Boolean RemoveFieldFromUniqBioSource (UniqBioSourcePtr u, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (u != NULL) {
    rval = RemoveMatchingFieldType (&(u->available_fields), field);
  }
  return rval;
}


static ValNodePtr UniqBioSourceListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = UniqBioSourceFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list= list_next;
  }
  return list;
}


static ValNodePtr UniqBioSourceListCopy (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp;
  UniqBioSourcePtr u;

  while (orig != NULL) {
    u = (UniqBioSourcePtr) orig->data.ptrvalue;
    if (u != NULL && u->biop != NULL) {
      vnp = ValNodeNew (prev);
      vnp->choice = 0;
      vnp->data.ptrvalue = UniqBioSourceCopy (u);
      if (prev == NULL) {
        list = vnp;
      }
      prev = vnp;
    }
    orig = orig->next;
  }
  return list;
}


static int LIBCALLBACK SortUniqBioSource (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return CompareUniqBioSource (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static ValNodePtr UniqBioSourceListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortUniqBioSource);
  return orig;
}


static Boolean AddQualToUniqBioSourceList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;
  while (list != NULL) {
    rval |= AddQualToUniqBioSource (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


static Boolean RemoveFieldFromUniqBioSourceList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;

  while (list != NULL) {
    rval |= RemoveFieldFromUniqBioSource (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


/* The UniqBioSrcGrp structure contains a list of UniqBioSrc items
 * for which the contents of the descriptive strings list are identical,
 * i.e., all the organisms in the group would have the same description
 * if you used the modifiers used to generate this list of strings.
 * The structure also contains the number of organisms in the list
 * so that it will be easy to tell that the UniqBioSrcGrp now contains a
 * single organism with a unique description.
 */
typedef struct uniqbiosrcgrp {
  Int4 num_biop;
  ValNodePtr biop_list;
} UniqBioSrcGrpData, PNTR UniqBioSrcGrpPtr;


static UniqBioSrcGrpPtr UniqBioSrcGrpNew (ValNodePtr biop_list)
{
  UniqBioSrcGrpPtr g;

  g = (UniqBioSrcGrpPtr) MemNew (sizeof (UniqBioSrcGrpData));
  g->num_biop = ValNodeLen (biop_list);
  g->biop_list = UniqBioSourceListCopy (biop_list);
  return g;
}


static UniqBioSrcGrpPtr UniqBioSrcGrpFree (UniqBioSrcGrpPtr g)
{
  if (g != NULL) {
    g->biop_list = UniqBioSourceListFree (g->biop_list);
    g = MemFree (g);
  }
  return g;
}


static UniqBioSrcGrpPtr UniqBioSrcGrpCopy (UniqBioSrcGrpPtr orig)
{
  UniqBioSrcGrpPtr g;

  if (orig == NULL) {
    return NULL;
  }

  g = UniqBioSrcGrpNew (orig->biop_list);
  return g;
}


static Boolean AddQualToUniqBioSrcGrp (UniqBioSrcGrpPtr g, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (g != NULL) {
    rval = AddQualToUniqBioSourceList (g->biop_list, field);
    if (rval) {
      g->biop_list = UniqBioSourceListSort (g->biop_list);
    }
  }
  return rval;
}


static Boolean RemoveFieldFromUniqBioSrcGrp (UniqBioSrcGrpPtr g, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (g != NULL) {
    rval = RemoveFieldFromUniqBioSourceList (g->biop_list, field);
  }
  return rval;
}


static FieldTypePtr GetAllPresentQualsForGroup (UniqBioSrcGrpPtr g)
{
  ValNodePtr vnp;
  FieldTypePtr match_list = NULL, ft, ft_next, ft_prev;
  UniqBioSourcePtr u;

  if (g == NULL || g->num_biop < 2) {
    return NULL;
  }

  u = g->biop_list->data.ptrvalue;
  match_list = FieldTypeListCopy (u->available_fields);
  for (vnp = g->biop_list->next; vnp != NULL && match_list != NULL; vnp = vnp->next) {
    u = vnp->data.ptrvalue;
    ft = match_list;
    ft_prev = NULL;
    while (ft != NULL) {
      ft_next = ft->next;
      if (ListHasMatchingFieldType (u->available_fields, ft)) {
        ft_prev = ft;
      } else{
        if (ft_prev == NULL) {
          match_list = ft->next;
        } else {
          ft_prev->next = ft->next;
        }
        ft->next = NULL;
        ft = FieldTypeFree (ft);
      }
      ft = ft_next;
    }
  }
  return match_list;
}


static FieldTypePtr GetAllQualsForGroup (UniqBioSrcGrpPtr g)
{
  ValNodePtr vnp;
  FieldTypePtr field_list = NULL;
  UniqBioSourcePtr u;

  if (g == NULL || g->num_biop < 2) {
    return NULL;
  }

  for (vnp = g->biop_list; vnp != NULL; vnp = vnp->next) {
    u = vnp->data.ptrvalue;
    if (u != NULL && u->available_fields != NULL) {
      ValNodeLink (&field_list, FieldTypeListCopy (u->available_fields));
    }
  }
  SortUniqueFieldTypeList (&field_list);
  return field_list;
}


static ValNodePtr UniqBioSrcGrpListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = UniqBioSrcGrpFree (list->data.ptrvalue);
    list = ValNodeFree (list); 
    list = list_next;
  }
  return list;
}


static ValNodePtr UniqBioSrcGrpListCopy (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp;
  UniqBioSrcGrpPtr g;

  while (orig != NULL) {
    g = (UniqBioSrcGrpPtr) orig->data.ptrvalue;
    if (g != NULL) {
      vnp = ValNodeNew (prev);
      vnp->choice = 0;
      vnp->data.ptrvalue = UniqBioSrcGrpCopy (g);
      if (prev == NULL) {
        list = vnp;
      }
      prev = vnp;
    }
    orig = orig->next;
  }
  return list;
}


/* NOTE - we want to sort groups from most biops to least biops */
static int LIBCALLBACK SortUniqBioSrcGrp (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  UniqBioSrcGrpPtr g1, g2;
  int              rval = 0;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL && vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      g1 = (UniqBioSrcGrpPtr) vnp1->data.ptrvalue;
      g2 = (UniqBioSrcGrpPtr) vnp2->data.ptrvalue;
      if (g1->num_biop > g2->num_biop) {
        rval = -1;
      } else if (g1->num_biop < g2->num_biop) {
        rval = 1;
      }
    }
  }
  return rval;
}


static ValNodePtr BioSrcGrpListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortUniqBioSrcGrp);
  return orig;
}


static Boolean RemoveFieldFromUniqBioSrcGrpList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;

  while (list != NULL) {
    rval |= RemoveFieldFromUniqBioSrcGrp (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


static void ReGroupUniqBioSrcGrpList (ValNodePtr list)
{
  ValNodePtr list_next, vnp;
  UniqBioSrcGrpPtr g, g2;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    vnp = g->biop_list;
    while (vnp != NULL && vnp->next != NULL 
           && CompareUniqBioSource (vnp->data.ptrvalue, vnp->next->data.ptrvalue) == 0) {
      vnp = vnp->next;
    }
    if (vnp != NULL && vnp->next != NULL) {
      g2 = UniqBioSrcGrpNew (NULL);
      g2->biop_list = vnp->next;
      g2->num_biop = ValNodeLen (vnp->next);
      vnp->next = NULL;
      g->num_biop -= g2->num_biop;
      list_next = ValNodeNew (NULL);
      list_next->data.ptrvalue = g2;
      list_next->next = list->next;
      list->next = list_next;
    }
    list = list->next;
  }
}


static Int4 FindMaxOrgsInUniqBioSrcGrpList (ValNodePtr list)
{
  Int4 max = 0;
  UniqBioSrcGrpPtr g;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    if (g != NULL && g->num_biop > max) {
      max = g->num_biop;
    }
    list = list->next;
  }
  return max;
}


static Int4 CountUniqueOrgsInUniqBioSrcGrpList (ValNodePtr list)
{
  Int4 count = 0;
  UniqBioSrcGrpPtr g;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    if (g != NULL && g->num_biop == 1) {
      count++;
    }
    list = list->next;
  }
  return count;
}    


static Boolean AddQualToUniqBioSrcGrpList (ValNodePtr list, FieldTypePtr field)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;

  vnp = list;
  while (vnp != NULL) {
    rval |= AddQualToUniqBioSrcGrp (vnp->data.ptrvalue, field);
    vnp = vnp->next;
  }
  if (rval) {
    /* regroup */
    ReGroupUniqBioSrcGrpList (list);
  }
  return rval;
}


typedef struct qualcombo {
  Int4         num_groups;
  Int4         num_mods;
  Int4         max_orgs_in_group;
  Int4         num_unique_orgs;
  FieldTypePtr field_list;
  ValNodePtr   group_list;
} QualComboData, PNTR QualComboPtr;


/* This function creates a new ModifierCombination item using the supplied
 * OrgGroup list.  It calculates the number of groups, maximum number of 
 * organisms in any one group, and number of unique organisms.
 * Initially there are no modifiers.
 */
static QualComboPtr QualComboNew (ValNodePtr grp_list)
{
  QualComboPtr newm;

  newm = (QualComboPtr) MemNew (sizeof (QualComboData));
  if (newm == NULL) return NULL;

  newm->num_mods = 0;
  newm->field_list = NULL;
  
  /* copy groups */
  newm->group_list = UniqBioSrcGrpListCopy (grp_list);
 
  newm->max_orgs_in_group = FindMaxOrgsInUniqBioSrcGrpList (newm->group_list);
  newm->num_unique_orgs = CountUniqueOrgsInUniqBioSrcGrpList (newm->group_list);
  newm->num_groups = ValNodeLen (newm->group_list);

  return newm; 
}


/* The CopyQualCombo creates a copy of a QualCombo item.
 * This includes creating a copy of the number and list of modifiers
 * and a copy of the number and list of OrgGroups, as well as copying the
 * maximum number of organisms in any one group and the number of unique
 * organism descriptions produced by this combination of modifiers.
 */
static QualComboPtr QualComboCopy (
  QualComboPtr m
)
{
  QualComboPtr newm;

  newm = QualComboNew (m->group_list);
  if (newm == NULL) return NULL;

  newm->field_list = FieldTypeListCopy (m->field_list);
  newm->num_mods = m->num_mods;
  
  return newm; 
}

/* This function frees the memory associated with a list of
 * ModifierCombination items.
 */
static QualComboPtr QualComboFree (
  QualComboPtr m
)
{
  if (m != NULL) {
    m->group_list = UniqBioSrcGrpListFree (m->group_list);
    m->field_list = FieldTypeListFree (m->field_list);
    m = MemFree (m);
  }
  return m;
}


static Boolean AddQualToQualCombo (
  QualComboPtr m,
  FieldTypePtr field
)
{
  Boolean    rval = FALSE;

  if (m == NULL || field == NULL) return rval;
  
  
  if (AddQualToUniqBioSrcGrpList (m->group_list, field)) {
    ValNodeLink (&(m->field_list), AsnIoMemCopy (field, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite));
    m->field_list = ValNodeSort (m->field_list, SortFieldTypeByImportance);
    m->num_mods ++;
    m->group_list = BioSrcGrpListSort (m->group_list);
    m->max_orgs_in_group = FindMaxOrgsInUniqBioSrcGrpList (m->group_list);
    m->num_unique_orgs = CountUniqueOrgsInUniqBioSrcGrpList (m->group_list);
    m->num_groups = ValNodeLen (m->group_list);
    rval = TRUE;
  }
  return rval;
}


static ValNodePtr LIBCALLBACK QualComboListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = QualComboFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


/* NOTE - we want to sort groups from most unique organisms to least unique organisms */
/* secondary sort - most groups to least groups */
/* tertiary sort - fewer max orgs in group to most max orgs in group */
/* fourth sort - least mods to most mods */
static int LIBCALLBACK SortQualCombo (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  QualComboPtr  g1, g2;
  FieldTypePtr  field1, field2;
  int           rval = 0;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL && vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      g1 = (QualComboPtr) vnp1->data.ptrvalue;
      g2 = (QualComboPtr) vnp2->data.ptrvalue;
      if (g1->num_unique_orgs > g2->num_unique_orgs) {
        rval = -1;
      } else if (g1->num_unique_orgs < g2->num_unique_orgs) {
        rval = 1;
      } else if (g1->num_groups > g2->num_groups) {
        rval = -1;
      } else if (g1->num_groups < g2->num_groups) {
        rval = 1;
      } else if (g1->max_orgs_in_group < g2->max_orgs_in_group) {
        rval = -1;
      } else if (g1->max_orgs_in_group > g2->max_orgs_in_group) {
        rval = 1;
      } else if (g1->num_mods < g2->num_mods) {
        rval = -1;
      } else if (g1->num_mods > g2->num_mods) {
        rval = 1;
      } else {
        /* compare modifiers */
        field1 = g1->field_list;
        field2 = g2->field_list;
        while (field1 != NULL && field2 != NULL 
               && (rval = CompareFieldTypeByImportance (field1, field2)) == 0) {
          field1 = field1->next;
          field2 = field2->next;
        }
        if (rval == 0) {
          if (field1 == NULL && field2 != NULL) {
            rval = -1;
          } else if (field1 != NULL && field2 == NULL) {
            rval = 1;
          }
        }
      }
    }
  }
  return rval;
}


static ValNodePtr QualComboListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortQualCombo);
  ValNodeUnique (&orig, SortQualCombo, QualComboListFree);
  return orig;
}



static ValNodePtr ExpandOneComboListUsingAllPresentQuals (QualComboPtr q)
{
  ValNodePtr new_list = NULL, vnp, vnp_m;
  FieldTypePtr match_list;
  UniqBioSrcGrpPtr g;
  QualComboPtr q_new;
  Boolean      found_group_improvement = FALSE;

  if (q == NULL) return NULL;
  for (vnp = q->group_list; vnp != NULL && !found_group_improvement; vnp = vnp->next) {
    g = (UniqBioSrcGrpPtr) vnp->data.ptrvalue;
    if (g->num_biop == 1) break;
    match_list = GetAllPresentQualsForGroup (g);
    for (vnp_m = match_list; vnp_m != NULL; vnp_m = vnp_m->next) {
      q_new = QualComboCopy (q);
      if (AddQualToQualCombo (q_new, vnp_m) && q_new->num_groups > q->num_groups) {
        ValNodeAddPointer (&new_list, 0, q_new);
        found_group_improvement = TRUE;
      } else {
        q_new = QualComboFree (q_new);
        RemoveFieldFromUniqBioSrcGrp (g, vnp_m);
      }
    }
    match_list = FieldTypeListFree (match_list);
  }
  return new_list;
}


static ValNodePtr ExpandOneComboListUsingAnyPresentQuals (QualComboPtr q)
{
  ValNodePtr new_list = NULL, vnp, vnp_m;
  FieldTypePtr match_list;
  UniqBioSrcGrpPtr g;
  QualComboPtr q_new;
  Boolean      found_group_improvement = FALSE;

  if (q == NULL) return NULL;
  for (vnp = q->group_list; vnp != NULL && !found_group_improvement; vnp = vnp->next) {
    g = (UniqBioSrcGrpPtr) vnp->data.ptrvalue;
    if (g->num_biop == 1) break;
    match_list = GetAllQualsForGroup (g);
    for (vnp_m = match_list; vnp_m != NULL; vnp_m = vnp_m->next) {
      q_new = QualComboCopy (q);
      if (AddQualToQualCombo (q_new, vnp_m) && q_new->num_groups > q->num_groups) {
        ValNodeAddPointer (&new_list, 0, q_new);
        found_group_improvement = TRUE;
      } else {
        q_new = QualComboFree (q_new);
        RemoveFieldFromUniqBioSrcGrp (g, vnp_m);
      }
    }
    match_list = FieldTypeListFree (match_list);
  }
  return new_list;
}


static Boolean ExpandComboList (ValNodePtr PNTR list)
{
  QualComboPtr  q;
  ValNodePtr    new_list, vnp, vnp_next, prev = NULL;
  Boolean       any_expansion = FALSE;

  if (*list == NULL) return FALSE;
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    q = (QualComboPtr) vnp->data.ptrvalue;
    new_list = ExpandOneComboListUsingAnyPresentQuals (q);

    if (new_list == NULL) {
      prev = vnp;
    } else {
      if (prev == NULL) {
        *list = new_list; 
      } else {
        prev->next = new_list;
      }
      prev = new_list;
      while (prev->next != NULL) {
        prev = prev->next;
      }
      ValNodeLink (&new_list, vnp->next);
      vnp->next = NULL;
      vnp = QualComboListFree (vnp);
      any_expansion = TRUE;
    }
    vnp = vnp_next;
  }
  return any_expansion;
}


static void BuildUniqBioSrcList (
  BioSourcePtr biop,
  Pointer userdata
)
{
  UniqBioSourcePtr u;

  u = UniqBioSourceNew (biop);
  ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, u);
}


/* The function FindBestQualCombo tries to find the best combination of modifiers
 * to create unique organism descriptions.  This is accomplished by
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * The function then makes copies of all of the new items added to the list,
 * starting with the item pointed to by start_of_expand, and adds another
 * modifier to each combination, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If the list of possible combinations is exhausted before each organism
 * has a unique description, the function selects the combination from the
 * list with the largest number of unique organism descriptions.  If more 
 * than one combination produces the largest number of unique organisms,
 * the combination with the largest number of unique organisms and the
 * largest number of groups will be selected.
 */
static QualComboPtr FindBestQualComboEx(ValNodePtr PNTR biop_list, ModifierItemLocalPtr ItemList)
{
  QualComboPtr initial_combo = NULL, best_combo = NULL;
  ValNodePtr   group_list = NULL, combo_list = NULL;
  UniqBioSrcGrpPtr g;
  SourceQualChoice scd;
  FieldType ft;
  Int4      i, qual;

  if (biop_list == NULL || *biop_list == NULL) {
    return NULL;
  }

  /* sort organisms */
  *biop_list = UniqBioSourceListSort (*biop_list);

  /* create group list */
  g = UniqBioSrcGrpNew (*biop_list);
  ValNodeAddPointer (&group_list, 0, g);  

  ReGroupUniqBioSrcGrpList (group_list);
  group_list = BioSrcGrpListSort (group_list);

  /* create combo with just the org groups */
  initial_combo = QualComboNew (group_list);
  group_list = UniqBioSrcGrpListFree (group_list);
  if (initial_combo == NULL) return NULL;

  /* add required quals */
  ft.choice = FieldType_source_qual;
  ft.data.ptrvalue = &scd;
  ft.next = NULL;
  scd.choice = SourceQualChoice_textqual;
  if (ItemList == NULL) {
    /* endogenous virus name */
    scd.data.intvalue = Source_qual_endogenous_virus_name;
    AddQualToQualCombo (initial_combo, &ft);
    /* plasmid name */
    scd.data.intvalue = Source_qual_plasmid_name;
    AddQualToQualCombo (initial_combo, &ft);
    /* transgenic */
    scd.data.intvalue = Source_qual_transgenic;
    AddQualToQualCombo (initial_combo, &ft);
  } else {
    for (i = 0; i < numDefLineModifiers; i++) {
      if (ItemList[i].required) {
        qual = GetSrcQualFromSubSrcOrOrgMod (DefLineModifiers[i].subtype, DefLineModifiers[i].isOrgMod);
        if (qual > -1) {
          scd.data.intvalue = qual;
          AddQualToQualCombo (initial_combo, &ft);
        }
      }
    }
  }  

  if (initial_combo->max_orgs_in_group == 1)
  {
    /* we're done - they're all unique */
    return initial_combo;
  }

  /* they're not unique yet.  Need to find a combination of modifiers that will make this as unique as possible */
  ValNodeAddPointer (&combo_list, 0, initial_combo);
  best_combo = initial_combo;
  while (ExpandComboList (&combo_list)) {
    /* sort after expansion */  
    combo_list = QualComboListSort (combo_list);
    best_combo = combo_list->data.ptrvalue;
    if (best_combo->max_orgs_in_group == 1) {
      break;
    } 
  }
  best_combo = QualComboCopy (best_combo);
  combo_list = QualComboListFree (combo_list);
  return best_combo;
}


static QualComboPtr FindBestQualCombo(SeqEntryPtr sep, ModifierItemLocalPtr ItemList)
{
  QualComboPtr best_combo;
  ValNodePtr   biop_list = NULL;

  /* first, get list of organisms */
  VisitBioSourcesInSep (sep, &biop_list, BuildUniqBioSrcList);

  best_combo = FindBestQualComboEx (&biop_list, ItemList);

  biop_list = UniqBioSourceListFree (biop_list);
  return best_combo;
}


static ModifierCombinationPtr ModifierCombinationFromQualCombo (QualComboPtr q)
{
  ModifierCombinationPtr m;
  FieldTypePtr           field;
  Int4                   i;

  if (q == NULL) {
    return NULL;
  }

  m = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  m->num_groups = q->num_groups;
  m->num_mods = q->num_mods;
  m->max_orgs_in_group = q->max_orgs_in_group;
  m->num_unique_orgs = q->num_unique_orgs;
  m->next = NULL;
  m->group_list = NULL;
  m->modifier_indices = NULL;
  for (field = q->field_list; field != NULL; field = field->next) {
    i = GetDeflinePosForFieldType (field);
    if (i >= 0) {
      ValNodeAddInt (&(m->modifier_indices), 0, i);
    }
  }
  return m;
}


NLM_EXTERN ValNodePtr FindBestModifiersForDeflineClauseList (
  ValNodePtr defline_clauses,
  ModifierItemLocalPtr ItemList
)

{
  QualComboPtr best_combo;
  ValNodePtr   biop_list = NULL, vnp;
  DefLineFeatClausePtr df;
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  UniqBioSourcePtr  u;
  ModifierCombinationPtr m;
  ValNodePtr modifier_indices = NULL;

  /* first, create list of organisms */
  for (vnp = defline_clauses; vnp != NULL; vnp = vnp->next) {
    df = (DefLineFeatClausePtr) vnp->data.ptrvalue;
    if (df != NULL) {
      sdp = SeqMgrGetNextDescriptor (df->bsp, NULL, Seq_descr_source, &context);
      if (sdp != NULL && sdp->data.ptrvalue != NULL) {
        u = UniqBioSourceNew (sdp->data.ptrvalue);
        ValNodeAddPointer (&(u->strings), 0, StringSave (df->clauselist));
        ValNodeAddPointer (&biop_list, 0, u);
      }
    }
  }

  best_combo = FindBestQualComboEx (&biop_list, ItemList);

  biop_list = UniqBioSourceListFree (biop_list);
  m = ModifierCombinationFromQualCombo (best_combo);
  if (m != NULL) {
    modifier_indices = CopyModifierIndices (m->modifier_indices);
    FreeModifierCombo (m);
  }
  return modifier_indices;
}


NLM_EXTERN ValNodePtr FindBestModifiersEx(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList,
  Boolean use_new
)

{
  ModifierCombinationPtr m;
  QualComboPtr q;
  ValNodePtr modifier_indices = NULL;

  if (use_new) {
    q = FindBestQualCombo (sep, ItemList);
    m = ModifierCombinationFromQualCombo (q);
  } else {
    m = FindBestCombo (sep, ItemList);
  }
  modifier_indices = CopyModifierIndices (m->modifier_indices);
  FreeModifierCombo (m);
  return modifier_indices;
}


NLM_EXTERN ValNodePtr FindBestModifiers(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList
)

{
  return FindBestModifiersEx (sep, ItemList, FALSE);
}


/* In this test function, we create a list of biosources with various combinations of modifiers,
 * and then calculate the best combination to use for the organism description.
 */
static CharPtr strings1[] = {"a", "b", "c"};
static CharPtr strings2[] = {"foo", "bar", "baz"};
static CharPtr strings3[] = {"d", "e", "f"};

static void SetBiopQual (BioSourcePtr biop, Int4 qual, CharPtr val)
{
  OrgModPtr mod;
  SubSourcePtr ssp;

  if (DefLineModifiers[qual].isOrgMod)
  {
    mod = OrgModNew ();
    mod->subtype = DefLineModifiers[qual].subtype;
    mod->subname = StringSave (val);
    mod->next = biop->org->orgname->mod;
    biop->org->orgname->mod = mod;
  } else {
    ssp = SubSourceNew ();
    ssp->subtype = DefLineModifiers[qual].subtype;
    ssp->name = StringSave (val);
    ssp->next = biop->subtype;
    biop->subtype = ssp;
  }
}


static void ClearBiopQuals (BioSourcePtr biop)
{
  biop->org->orgname->mod = OrgModSetFree (biop->org->orgname->mod);
  biop->subtype = SubSourceSetFree (biop->subtype);
}


static void PrintBiopQuals (BioSourcePtr biop, FILE *fp)
{
  OrgModPtr mod;
  SubSourcePtr ssp;

  fprintf (fp, "Taxname: %s", biop->org->taxname);
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    fprintf (fp, "\tOrgMod%d:%s", mod->subtype, mod->subname);
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    fprintf (fp, "\tSubSource%d:%s", ssp->subtype, ssp->name);
  }
  fprintf (fp, "\n");
}


static void PrintModifiers (ValNodePtr modifiers, FILE *fp)
{
  ValNodePtr vnp;

  if (modifiers == NULL) {
    fprintf (fp, "\tNo combo");
  } 
  for (vnp = modifiers; vnp != NULL; vnp = vnp->next) {
    fprintf (fp, "\t%s:%d", DefLineModifiers[vnp->data.intvalue].isOrgMod ? "OrgMod" : "SubSource", 
                            DefLineModifiers[vnp->data.intvalue].subtype);
  }
  fprintf (fp, "\n");
}


static Boolean IsNonTextDeflineQual (Int4 srcqual)
{
  if (srcqual == DEFLINE_POS_Transgenic
      || srcqual == DEFLINE_POS_Germline
      || srcqual == DEFLINE_POS_Metagenomic
      || srcqual == DEFLINE_POS_Environmental_sample
      || srcqual == DEFLINE_POS_Rearranged)
  {
    return TRUE;  
  }
  else
  {
    return FALSE;
  }
}


static void CreateOneTest (FILE *fp, Int4 i, Int4 j, Int4 k, BioSourcePtr PNTR biops, Int4 num_biops, Boolean vary1, Boolean vary2, Boolean vary3)
{
  Int4 n;
  ValNodePtr uniq_biop_list = NULL;
  QualComboPtr q;
  ModifierCombinationPtr m;

  for (n = 0; n < num_biops; n++) {
    if (i < numDefLineModifiers) {
      if (vary1) {
        SetBiopQual (biops[n], i, strings1[n]);
      } else {
        SetBiopQual (biops[n], i, strings1[0]);
      }
    }
    if (j < numDefLineModifiers) {
      if (vary2) {
        SetBiopQual (biops[n], j, strings2[n]);
      } else {
        SetBiopQual (biops[n], j, strings2[0]);
      }
    }
    if (k < numDefLineModifiers) {
      if (vary3) {
        SetBiopQual (biops[n], k, strings3[n]);
      } else {
        SetBiopQual (biops[n], k, strings3[0]);
      }
    }
    ValNodeAddPointer (&uniq_biop_list, 0, UniqBioSourceNew (biops[n]));
  }
  q = FindBestQualComboEx (&uniq_biop_list, NULL);
  m = ModifierCombinationFromQualCombo (q);
  /* print results */
  for (n = 0; n < num_biops; n++) {
    PrintBiopQuals (biops[n], fp);
  }
  PrintModifiers (m->modifier_indices, fp);
  q = QualComboFree (q);
  FreeModifierCombo (m);
  uniq_biop_list = UniqBioSourceListFree (uniq_biop_list);

  /* clear quals */
  for (n = 0; n < num_biops; n++) {
    ClearBiopQuals (biops[n]);
  }
  
}


extern void TestFindBestQualCombo (FILE *fp)
{
  BioSourcePtr biops[3];
  Int4         num_biops = 3;
  Int4         i, j, k;

  for (i = 0; i < num_biops; i++) {
    biops[i] = BioSourceNew ();
    biops[i]->org = OrgRefNew ();
    biops[i]->org->orgname = OrgNameNew ();
  }
  
  /* first try with all organism names the same */
  for (i = 0; i < num_biops; i++) {
    biops[i]->org->taxname = StringSave ("Homo sapiens");
  }

  for (i = 0; i <= numDefLineModifiers; i++) {
    if (IsNonTextDeflineQual (i)) {
      continue;
    }
    for (j = i + 1; j <= numDefLineModifiers; j++) {
      if (IsNonTextDeflineQual (j)) {
        continue;
      }
      for (k = j + 1; k <= numDefLineModifiers; k++) {
        if (IsNonTextDeflineQual (k)) {
          continue;
        }
        if (k != numDefLineModifiers) {
          /* try all the same but 1*/
          CreateOneTest (fp, i, j, k, biops, num_biops, FALSE, FALSE, TRUE);
        }
        if (j != numDefLineModifiers || k != numDefLineModifiers) {
          /* try 2 different 1 same*/
          CreateOneTest (fp, i, j, k, biops, num_biops, FALSE, TRUE, TRUE);
        }
        if (i != numDefLineModifiers || j != numDefLineModifiers || k != numDefLineModifiers) {
          /* try all different */
          CreateOneTest (fp, i, j, k, biops, num_biops, TRUE, TRUE, TRUE);
        }
      }
    }
  }
  for (i = 0; i < num_biops; i++) {
    biops[i] = BioSourceFree (biops[i]);
  }
 
}




/* collection_date has a controlled format.  
 * It is YYYY or Mmm-YYYY or DD-Mmm-YYYY where Mmm = Jan, Feb, Mar, Apr, May, 
 *                                                   Jun, Jul, Aug, Sep, Oct, 
 *                                                   Nov, Dec
 * This function will convert other formats  to this format.
 * For instance, September 12, 2004 should be converted to 12-Sep-2004
 * 12/15/2003 should be converted to 15-Dec-2003.  
 * 
 * If the date supplied is ambiguous (01/03/05), can you allow the indexer to choose which field goes in Mmm and which in DD.
 */

NLM_EXTERN Int4 ReadNumberFromToken (CharPtr token, Int4 token_len)
{
  Int4 val = 0;
  
  if (token == NULL || !isdigit (*token))
  {
    return val;
  }
  while (token_len > 0)
  {
    val *= 10;
    val += *token - '0';
    token++;
    token_len--;
  }
  
  return val;
}

static Int4 GetYearFromNumber(Int4 year)
{
	Nlm_DayTime dt;

  if (year < 1000)
  {
    GetDayTime (&dt);
    if (year + 2000 > dt.tm_year + 1901)
    {
      year += 1900;
    }
    else
    {
      year += 2000;
    }
  }
  return year;
}

NLM_EXTERN Int4 GetYearFromToken (CharPtr token, Int4 token_len)
{
  Int4        year = 0;
  
  if (token == NULL || token_len == 0 || token_len > 4)
  {
    return 0;
  }
  
  year = GetYearFromNumber(ReadNumberFromToken (token, token_len));
  
  return year;
}

static CharPtr month_abbrevs [12] =
{
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};


NLM_EXTERN CharPtr GetMonthAbbrev (Int4 n)
{
  if (n > 0 && n <= 12) {
    return month_abbrevs[n - 1];
  } else {
    return NULL;
  }
}
  

static Int4 days_in_month [12] =
{
  31, 29, 31, 30, 31, 30,
  31, 31, 30, 31, 30, 31
};


NLM_EXTERN Int4 GetDaysInMonth (Int4 n)
{
  if (n > 0 && n <= 12) {
    return days_in_month[n - 1];
  } else {
    return 0;
  }
}


NLM_EXTERN Int4 GetMonthNumFromAbbrev (CharPtr month_abbrev) 
{
  Int4 i;

  for (i = 0; i < 12; i++) {
    if (StringICmp (month_abbrev, month_abbrevs[i]) == 0) {
      return i;
    }
  }
  return -1;
}

static Int4 GetDaysInMonthByName (CharPtr month)
{
  Int4 month_num;
  
  for (month_num = 0; month_num < 12; month_num++)
  {
    if (StringCmp (month, month_abbrevs [month_num]) == 0)
    {
      return days_in_month [month_num];
    }
  }
  return 0;
}

NLM_EXTERN CharPtr GetMonthFromToken (CharPtr token, Int4 token_len)
{
  Int4    month_num;
  
  if (token == NULL || token_len == 0)
  {
    return NULL;
  }
  
  if (isdigit (*token)) 
  {
    if (token_len > 2)
    {
      return NULL;
    }
    else
    {
      month_num = ReadNumberFromToken (token, token_len);
      if (month_num == 0 || month_num > 12)
      {
        return NULL;
      }
      else
      {
        return month_abbrevs [month_num - 1];
      }
    }
  }
  else
  {
    for (month_num = 0; month_num < 12; month_num++)
    {
      if (StringNICmp (token, month_abbrevs[month_num], 3) == 0)
      {
        return month_abbrevs[month_num];
      }
    }
    return NULL;
  }
}

static Boolean
ChooseDayAndYear 
(Int4    num_1,
 Int4    num_2,
 CharPtr month,
 Boolean year_first,
 Int4Ptr day,
 Int4Ptr year)
{  
  if (day == NULL || year == NULL)
  {
    return FALSE;
  }
  
  if (num_1 == 0 && num_2 == 0)
  {
    return FALSE;
  }
  else if (num_1 == 0)
  {
    *year = 2000;
    *day = num_2;
  }
  else if (num_2 == 0)
  {
    *year = 2000;
    *day = num_1;
  }
  else if (num_1 > GetDaysInMonthByName (month))
  {
    *year = num_1;
    *day = num_2;
  }
  else if (num_2 > GetDaysInMonthByName (month))
  {
    *year = num_2;
    *day = num_1;
  }
  else if (year_first)
  {
    *year = num_1;
    *day = num_2;
  }
  else
  {
    *year = num_2;
    *day = num_1;
  }
  
  return TRUE;
}

static Boolean 
ChooseMonthAndYear
(Int4    num_1,
 Int4    num_2,
 Boolean month_first,
 CharPtr PNTR month,
 Int4Ptr year,
 BoolPtr month_ambiguous)
{
  if (year == NULL || month == NULL 
      || (num_1 == 0 && num_2 == 0)
      || (num_1 > 12 && num_2 > 12)
      || (num_1 == 0 && num_2 > 12)
      || (num_2 == 0 && num_1 > 12))
  {
    return FALSE;
  }
  
  if (num_1 == 0)
  {
    *year = 2000;
    *month = month_abbrevs[num_2 - 1];
  }
  else if (num_2 == 0)
  {
    *year = 2000;
    *month = month_abbrevs[num_1 - 1];
  }
  else if (num_1 > 12)
  {
    *year = GetYearFromNumber(num_1);
    *month = month_abbrevs [num_2 - 1];
  }
  else if (num_2 > 12)
  {
    *year = GetYearFromNumber(num_2);
    *month = month_abbrevs [num_1 - 1];
  }
  else if (month_first)
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *year = GetYearFromNumber(num_2);
    *month = month_abbrevs [num_1 - 1];
  }
  else
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *year = GetYearFromNumber(num_1);
    *month = month_abbrevs [num_2 - 1];
  }
  return TRUE;
}


static Boolean ChooseMonthAndDay 
(Int4    num_1,
 Int4    num_2,
 Boolean month_first,
 CharPtr PNTR month,
 Int4Ptr day,
 BoolPtr month_ambiguous)
{
  if (day == NULL || month == NULL || num_1 == 0 || num_2 == 0
      || (num_1 > 12 && num_2 > 12))
  {
    return FALSE;
  }
  
  if (num_1 > 12)
  {
    *day = num_1;
    *month = month_abbrevs [num_2 - 1];
  }
  else if (num_2 > 12)
  {
    *day = num_2;
    *month = month_abbrevs [num_1 - 1];
  }
  else if (month_first)
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *day = num_2;
    *month = month_abbrevs [num_1 - 1];
  }
  else
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *day = num_1;
    *month = month_abbrevs [num_2 - 1];
  }
  return TRUE;
}

NLM_EXTERN CharPtr ReformatDateStringEx (CharPtr orig_date, Boolean month_first, BoolPtr month_ambiguous)
{
  CharPtr reformatted_date = NULL, cp;
  Int4    year = 0, day = 0;
  CharPtr month = NULL;
  CharPtr token_list[3];
  Int4    token_lens[3];
  CharPtr numbers = "0123456789";
  CharPtr letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  Int4    num_tokens = 0;
  Int4    token_len;
  Int4    month_token = -1;
  Boolean is_num;
  Int4    num_1, num_2, num_3;
  
  if (StringHasNoText (orig_date))
  {
    return NULL;
  }
  
  /* divide our original date into tokens */
  /* skip over any leading spaces */
  cp = orig_date;
  while (*cp != 0 && num_tokens < 3)
  {
    is_num = FALSE;
    token_len = StringSpn (cp, numbers);  
    if (token_len == 0)
    {
      token_len = StringSpn (cp, letters);
    }
    else
    {
      is_num = TRUE;
    }
    if (token_len == 0)
    {
      cp++;
    }
    else
    {
      if (!is_num)
      {
        if (month_token == -1)
        {
          month_token = num_tokens;
        }
        else
        {
          /* already found a month string */
          return NULL;
        }
      }
      token_list [num_tokens] = cp;
      token_lens [num_tokens] = token_len;
      num_tokens ++;
      cp += token_len;
    }
  }
 
  if (num_tokens == 0 || *cp != 0)
  {
    return NULL;
  }

  if (num_tokens == 1)
  {
    if (month_token == 0)
    {
      return NULL;
    }
    year = GetYearFromToken (token_list [0], token_lens [0]);
  }
  else if (num_tokens == 2)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      year = GetYearFromToken (token_list [1], token_lens [1]);
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      year = GetYearFromToken (token_list [0], token_lens [0]);
    }
    else
    {
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      if (! ChooseMonthAndYear (num_1, num_2, month_first, &month, &year, month_ambiguous))
      {
        return NULL;
      }
    }
  }
  else if (num_tokens == 3)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      num_1 = ReadNumberFromToken (token_list [1], token_lens [1]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 2)
    {
      month = GetMonthFromToken (token_list [2], token_lens [2]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else
    {
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      num_3 = ReadNumberFromToken (token_list [2], token_lens [2]);
      
      if (num_1 > 31 || num_1 == 0)
      {
        year = num_1;
        if (! ChooseMonthAndDay (num_2, num_3, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_2 > 31 || num_2 == 0)
      {
        year = num_2;
        if (! ChooseMonthAndDay (num_1, num_3, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_3 > 31 || num_3 == 0)
      {
        year = num_3;
        if (! ChooseMonthAndDay (num_1, num_2, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_1 > 0 && num_1 < 13 && num_2 > days_in_month [num_1] && num_3 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_1 - 1];
        year = num_2;
        day = num_3;
      }
      else if (num_1 > 0 && num_1 < 13 && num_3 > days_in_month [num_1] && num_2 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_1 - 1];
        year = num_3;
        day = num_2;
      }
      else if (num_2 > 0 && num_2 < 13 && num_1 > days_in_month [num_2] && num_3 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_2 - 1];
        year = num_1;
        day = num_3;
      }
      else if (num_2 > 0 && num_2 < 13 && num_3 > days_in_month [num_2] && num_1 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_2 - 1];
        year = num_3;
        day = num_1;
      }
      else if (num_3 > 0 && num_3 < 13 && num_1 > days_in_month [num_3] && num_2 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_3 - 1];
        year = num_1;
        day = num_2;
      }
      else if (num_3 > 0 && num_3 < 13 && num_2 > days_in_month [num_3] && num_1 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_3 - 1];
        year = num_2;
        day = num_1;
      }
      else
      {
        year = num_3;
        if (! ChooseMonthAndDay (num_1, num_2, month_first, &month, &day, month_ambiguous))
        {
          year = num_1;
          if (!ChooseMonthAndDay (num_2, num_3, month_first, &month, &day, month_ambiguous))
          {
            return NULL;
          }
        }
      }
                
    }
    year = GetYearFromNumber(year);
  }
  
  if (month == NULL && day > 0)
  {
    return NULL;
  }
  
  reformatted_date = (CharPtr) MemNew (sizeof (Char) * 12);
  if (reformatted_date == NULL)
  {
    return NULL;
  }
   
  if (month == NULL)
  {
    sprintf (reformatted_date, "%d", year);
  }
  else if (day == 0)
  {
    sprintf (reformatted_date, "%s-%d", month, year);
  }
  else
  {
    sprintf (reformatted_date, "%02d-%s-%d", day, month, year);
  }
  return reformatted_date;
}

