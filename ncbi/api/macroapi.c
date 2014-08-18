/*   macro.c
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
* File Name:  macro.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/8/2007
*
* $Revision: 1.180 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <asn.h>
#include <objfeat.h>
#include <subutil.h>
#include <objmgr.h>
#include <objfdef.h>
#include <gbftdef.h>
#include <sqnutils.h>
#include <edutil.h>
#include <gather.h>
#include <ffprint.h>
#include <asn2gnbi.h>
#include <findrepl.h>
#include <utilpub.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>
#include <seqport.h>


static Boolean DoesFeatureMatchRnaType (SeqFeatPtr sfp, RnaFeatTypePtr rt);
static Int4 CompareRnaTypes (RnaFeatTypePtr rt1, RnaFeatTypePtr rt2);
static int LIBCALLBACK SortVnpByChoiceAndIntvalue (VoidPtr ptr1, VoidPtr ptr2);

/* NOTES */
/* When adding a new field type, add implementation to the following functions:
 * GetFromFieldFromFieldPair
 * GetToFieldFromFieldPair
 * BuildFieldPairFromFromField
 * FieldTypeChoiceFromFieldPairTypeChoice
 * CompareFieldTypes
 * IsObjectAppropriateForFieldValue
 * GetFieldValueForObject
 * RemoveFieldValueForObject
 * SetFieldValueForObject
 * GetObjectListForFieldType
 * GetFieldListForFieldType
 * IsFieldTypeEmpty
 * AllowFieldMulti
 * SummarizeFieldType
 * GetTargetListForRowAndColumn
 * ReportMissingTargets
 * CountObjectsForColumnFields
 */


NLM_EXTERN FeatureFieldPtr FeatureFieldCopy (FeatureFieldPtr orig)
{
  FeatureFieldPtr ff = NULL;

  if (orig != NULL) {
    ff = FeatureFieldNew();
    ff->type = orig->type;
    if (orig->field != NULL) {
      ff->field = AsnIoMemCopy (orig->field, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
    }
  }
  return ff;
}


NLM_EXTERN FieldTypePtr FieldTypeCopy (FieldTypePtr orig)
{
  FieldTypePtr ft = NULL;
  RnaQualPtr   rq, rq_orig;

  if (orig != NULL) {
    if (orig->data.ptrvalue == NULL) {
      ft = ValNodeNew (NULL);
      ft->choice = orig->choice;
    } else if (orig->choice == FieldType_feature_field) {
      ft = ValNodeNew (NULL);
      ft->choice = FieldType_feature_field;
      ft->data.ptrvalue = FeatureFieldCopy (orig->data.ptrvalue);
    } else if (orig->choice == FieldType_rna_field) {
      ft = ValNodeNew (NULL);
      ft->choice = FieldType_rna_field;
      rq_orig = (RnaQualPtr) orig->data.ptrvalue;
      rq = RnaQualNew ();
      rq->field = rq_orig->field;
      rq->type = AsnIoMemCopy (rq_orig->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      ft->data.ptrvalue = rq;
    } else {
      ft = AsnIoMemCopy (orig, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite);
    }
  }
  return ft;
}    


/* Functions for handling FieldPairs */
NLM_EXTERN FieldTypePtr GetFromFieldFromFieldPair (FieldPairTypePtr fieldpair)
{
  SourceQualChoicePtr ss = NULL;
  SourceQualPairPtr sqpp;
  FeatureFieldPairPtr fp;
  FeatureFieldPtr fs;
  RnaQualPairPtr rqp;
  RnaQualPtr rq;
  FieldTypePtr f = NULL;
  CDSGeneProtFieldPairPtr cp;
  MolinfoFieldPairPtr mp;
  StructuredCommentFieldPairPtr scfp;
  ValNodePtr vnp;

  if (fieldpair == NULL) return NULL;
  switch (fieldpair->choice) {
    case FieldPairType_source_qual:
      sqpp = (SourceQualPairPtr) fieldpair->data.ptrvalue;
      if (sqpp != NULL) {
        ss = ValNodeNew (NULL);
        ss->choice = SourceQualChoice_textqual;
        ss->data.intvalue = sqpp->field_from;
        f = ValNodeNew (NULL);
        f->choice = FieldType_source_qual;
        f->data.ptrvalue = ss;
      }
      break;
    case FieldPairType_feature_field:
      fp = (FeatureFieldPairPtr) fieldpair->data.ptrvalue;
      if (fp != NULL) {
        fs = FeatureFieldNew ();
        fs->type = fp->type;
        fs->field = (FeatQualChoicePtr) AsnIoMemCopy (fp->field_from, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
        f = ValNodeNew (NULL);
        f->choice = FieldType_feature_field;
        f->data.ptrvalue = fs;
      }
      break;
    case FieldPairType_rna_field:
      rqp = (RnaQualPairPtr) fieldpair->data.ptrvalue;
      if (rqp != NULL) {
        rq = RnaQualNew ();
        if (rqp->type != NULL) {
          rq->type = AsnIoMemCopy (rqp->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
        rq->field = rqp->field_from;
        f = ValNodeNew (NULL);
        f->choice = FieldType_rna_field;
        f->data.ptrvalue = rq;
      }
      break;
    case FieldPairType_cds_gene_prot:
      cp = (CDSGeneProtFieldPairPtr) fieldpair->data.ptrvalue;
      if (cp != NULL) {
        f = ValNodeNew (NULL);
        f->choice = FieldType_cds_gene_prot;
        f->data.intvalue = cp->field_from;
      }
      break;
    case FieldPairType_molinfo_field:
      mp = (MolinfoFieldPairPtr) fieldpair->data.ptrvalue;
      if (mp != NULL && mp->data.ptrvalue != NULL) {
        vnp = NULL;
        switch (mp->choice) {
          case MolinfoFieldPair_molecule:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_molecule;
            vnp->data.intvalue = ((MolinfoMoleculePairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_technique:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_technique;
            vnp->data.intvalue = ((MolinfoTechniquePairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_completedness:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_completedness;
            vnp->data.intvalue = ((MolinfoCompletednessPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_mol_class:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_mol_class;
            vnp->data.intvalue = ((MolinfoMolClassPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_topology:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_topology;
            vnp->data.intvalue = ((MolinfoTopologyPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_strand:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_strand;
            vnp->data.intvalue = ((MolinfoStrandPairPtr)mp->data.ptrvalue)->from;
            break;
        }
        if (vnp != NULL) {
          f = ValNodeNew (NULL);
          f->choice = FieldType_molinfo_field;
          f->data.ptrvalue = vnp;
        }
      }
      break;
    case FieldPairType_struc_comment_field:
      scfp = (StructuredCommentFieldPairPtr) fieldpair->data.ptrvalue;
      if (scfp != NULL) {        
        f = ValNodeNew (NULL);
        f->choice = FieldType_struc_comment_field;
        f->data.ptrvalue = AsnIoMemCopy (scfp->from, (AsnReadFunc) StructuredCommentFieldAsnRead, (AsnWriteFunc) StructuredCommentFieldAsnWrite);
      }
      break;
  }
  return f;
}


NLM_EXTERN FieldTypePtr GetToFieldFromFieldPair (FieldPairTypePtr fieldpair)
{
  SourceQualChoicePtr ss = NULL;
  SourceQualPairPtr sqpp;
  FeatureFieldPairPtr fp;
  FeatureFieldPtr fs;
  FieldTypePtr f = NULL;
  RnaQualPairPtr   rqp;
  RnaQualPtr       rq;
  CDSGeneProtFieldPairPtr cp;
  MolinfoFieldPairPtr     mp;
  StructuredCommentFieldPairPtr scfp;

  ValNodePtr              vnp;

  if (fieldpair == NULL) return NULL;
  switch (fieldpair->choice) {
    case FieldPairType_source_qual:
      sqpp = (SourceQualPairPtr) fieldpair->data.ptrvalue;
      if (sqpp != NULL) {
        ss = ValNodeNew (NULL);
        ss->choice = SourceQualChoice_textqual;
        ss->data.intvalue = sqpp->field_to;
        f = ValNodeNew (NULL);
        f->choice = FieldType_source_qual;
        f->data.ptrvalue = ss;
      }
      break;
    case FieldPairType_feature_field:
      fp = (FeatureFieldPairPtr) fieldpair->data.ptrvalue;
      if (fp != NULL) {
        fs = FeatureFieldNew ();
        fs->type = fp->type;
        fs->field = (FeatQualChoicePtr) AsnIoMemCopy (fp->field_to, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
        f = ValNodeNew (NULL);
        f->choice = FieldType_feature_field;
        f->data.ptrvalue = fs;
      }
      break;
    case FieldPairType_rna_field:
      rqp = (RnaQualPairPtr) fieldpair->data.ptrvalue;
      if (rqp != NULL) {
        rq = RnaQualNew ();
        if (rqp->type != NULL) {
          rq->type = AsnIoMemCopy (rqp->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
        rq->field = rqp->field_to;
        f = ValNodeNew (NULL);
        f->choice = FieldType_rna_field;
        f->data.ptrvalue = rq;
      }
      break;
    case FieldPairType_cds_gene_prot:
      cp = (CDSGeneProtFieldPairPtr) fieldpair->data.ptrvalue;
      if (cp != NULL) {
        f = ValNodeNew (NULL);
        f->choice = FieldType_cds_gene_prot;
        f->data.intvalue = cp->field_to;
      }
      break;
    case FieldPairType_molinfo_field:
      mp = (MolinfoFieldPairPtr) fieldpair->data.ptrvalue;
      if (mp != NULL && mp->data.ptrvalue != NULL) {
        vnp = NULL;
        switch (mp->choice) {
          case MolinfoFieldPair_molecule:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_molecule;
            vnp->data.intvalue = ((MolinfoMoleculePairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_technique:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_technique;
            vnp->data.intvalue = ((MolinfoTechniquePairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_completedness:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_completedness;
            vnp->data.intvalue = ((MolinfoCompletednessPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_mol_class:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_mol_class;
            vnp->data.intvalue = ((MolinfoMolClassPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_topology:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_topology;
            vnp->data.intvalue = ((MolinfoTopologyPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_strand:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_strand;
            vnp->data.intvalue = ((MolinfoStrandPairPtr)mp->data.ptrvalue)->to;
            break;
        }
        if (vnp != NULL) {
          f = ValNodeNew (NULL);
          f->choice = FieldType_molinfo_field;
          f->data.ptrvalue = vnp;
        }
      }
      break;     
    case FieldPairType_struc_comment_field:
      scfp = (StructuredCommentFieldPairPtr) fieldpair->data.ptrvalue;
      if (scfp != NULL) {        
        f = ValNodeNew (NULL);
        f->choice = FieldType_struc_comment_field;
        f->data.ptrvalue = AsnIoMemCopy (scfp->to, (AsnReadFunc) StructuredCommentFieldAsnRead, (AsnWriteFunc) StructuredCommentFieldAsnWrite);
      }
      break;
  }
  return f;
}


NLM_EXTERN FieldPairTypePtr BuildFieldPairFromFromField (FieldTypePtr field_from)
{
  SourceQualChoicePtr ss = NULL;
  SourceQualPairPtr sqpp;
  FeatureFieldPairPtr fp;
  FeatureFieldPtr fs;
  FieldTypePtr f = NULL;
  RnaQualPairPtr   rqp;
  RnaQualPtr       rq;
  CDSGeneProtFieldPairPtr cp;
  StructuredCommentFieldPairPtr scfp;
  ValNodePtr     mp;
  MolinfoMoleculePairPtr      mol_p;
  MolinfoTechniquePairPtr     tech_p;
  MolinfoCompletednessPairPtr comp_p;
  MolinfoMolClassPairPtr      class_p;
  MolinfoTopologyPairPtr      topo_p;
  MolinfoStrandPairPtr        strand_p;
  ValNodePtr              vnp;
  FieldPairTypePtr        pair = NULL;

  if (field_from == NULL) return NULL;
  switch (field_from->choice) {
    case FieldType_source_qual:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_source_qual;
      ss = (SourceQualChoicePtr) field_from->data.ptrvalue;
      if (ss != NULL && ss->choice == SourceQualChoice_textqual) {
        sqpp = SourceQualPairNew ();
        sqpp->field_from = ss->data.intvalue;
        pair->data.ptrvalue = sqpp;
      }
      break;
    case FieldType_feature_field:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_feature_field;
      fs = (FeatureFieldPtr) field_from->data.ptrvalue;
      if (fs != NULL) {
        fp = FeatureFieldPairNew ();
        fp->type = fs->type;
        fp->field_from = (FeatQualChoicePtr) AsnIoMemCopy (fs->field, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
        pair->data.ptrvalue = fp;
      }
      break;
    case FieldType_rna_field:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_rna_field;
      rq = (RnaQualPtr) field_from->data.ptrvalue;
      if (rq != NULL) {
        rqp = RnaQualPairNew ();
        if (rq->type != NULL) {
          rqp->type = AsnIoMemCopy (rq->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
        rqp->field_from = rq->field;
        pair->data.ptrvalue = rqp;
      }
      break;
    case FieldType_cds_gene_prot:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_cds_gene_prot;
      cp = CDSGeneProtFieldPairNew ();
      cp->field_from = field_from->data.intvalue;
      pair->data.ptrvalue = cp;
      break;
    case FieldType_molinfo_field:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_molinfo_field;
      vnp = field_from->data.ptrvalue;
      if (vnp != NULL) {
        switch (vnp->choice) {
          case MolinfoField_molecule:
            mol_p = MolinfoMoleculePairNew ();
            mol_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = mol_p;
            pair->data.ptrvalue = mp;
            break;
          case MolinfoField_technique:
            tech_p = MolinfoTechniquePairNew ();
            tech_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = tech_p;
            pair->data.ptrvalue = mp;
            break;
          case MolinfoField_completedness:
            comp_p = MolinfoCompletednessPairNew ();
            comp_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = comp_p;
            pair->data.ptrvalue = mp;
            break;
          case MolinfoField_mol_class:
            class_p = MolinfoMolClassPairNew ();
            class_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = class_p;
            pair->data.ptrvalue = mp;
            break;
          case MolinfoField_topology:
            topo_p = MolinfoTopologyPairNew ();
            topo_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = topo_p;
            pair->data.ptrvalue = mp;
            break;
          case MolinfoFieldPair_strand:
            strand_p = MolinfoStrandPairNew ();
            strand_p->from = vnp->data.intvalue;
            mp = ValNodeNew (NULL);
            mp->choice = MolinfoFieldPair_molecule;
            mp->data.ptrvalue = strand_p;
            pair->data.ptrvalue = mp;
            break;
        }
      }
      break;
    case FieldType_struc_comment_field:
      pair = ValNodeNew (NULL);
      pair->choice = FieldPairType_struc_comment_field;
      scfp = StructuredCommentFieldPairNew ();
      scfp->from = AsnIoMemCopy (field_from, (AsnReadFunc) StructuredCommentFieldAsnRead, (AsnWriteFunc) StructuredCommentFieldAsnWrite);
      pair->data.ptrvalue = scfp;
      break;
  }
  return pair;
}


NLM_EXTERN Uint1 FieldTypeChoiceFromFieldPairTypeChoice (Uint1 field_pair_choice)
{
  Uint1 field_type_choice = 0;

  switch (field_pair_choice) {
    case FieldPairType_source_qual:
      field_type_choice = FieldType_source_qual;
      break;
    case FieldPairType_feature_field:
      field_type_choice = FieldType_feature_field;
      break;
    case FieldPairType_rna_field:
      field_type_choice = FieldType_rna_field;
      break;
    case FieldPairType_cds_gene_prot:
      field_type_choice = FieldType_cds_gene_prot;
      break;
    case FieldPairType_molinfo_field:
      field_type_choice = FieldType_molinfo_field;
      break;
    case FieldPairType_struc_comment_field:
      field_type_choice = FieldType_struc_comment_field;
      break;
  }

  return field_type_choice;
}


/* functions for handling single fields */
NLM_EXTERN int CompareFieldTypes (FieldTypePtr vnp1, FieldTypePtr vnp2)
{
  int rval = 0;
  FeatureFieldPtr field1, field2;
  RnaQualPtr rq1, rq2;
  StructuredCommentFieldPtr scf1, scf2;
  
  if (vnp1 == NULL && vnp2 == NULL) {
    rval = 0;
  } else if (vnp1 == NULL) {
    rval = -1;
  } else if (vnp2 == NULL) {
    rval = 1;
  } else if (vnp1->choice > vnp2->choice) {
    rval = 1;
  } else if (vnp1->choice < vnp2->choice) {
    rval = -1;
  } else {
    switch (vnp1->choice) {
      case FieldType_source_qual:
      case FieldType_molinfo_field:
        vnp1 = vnp1->data.ptrvalue;
        vnp2 = vnp2->data.ptrvalue;
        rval = SortVnpByChoiceAndIntvalue (&vnp1, &vnp2);
        break;
      case FieldType_feature_field:
        field1 = (FeatureFieldPtr) vnp1->data.ptrvalue;
        field2 = (FeatureFieldPtr) vnp2->data.ptrvalue;
        if (field1 == NULL && field2 == NULL) {
          rval = 0;
        } else if (field1 == NULL) {
          rval = -1;
        } else if (field2 == NULL) {
          rval = 1;
        } else if (field1->type < field2->type) {
          rval = -1;
        } else if (field1->type > field2->type) {
          rval = 1;
        } else if (field1->field == NULL && field2->field == NULL) {
          rval = 0;
        } else if (field1->field == NULL) {
          rval = -1;
        } else if (field2->field == NULL) {
          rval = 1;
        } else if (field1->field->choice < field2->field->choice) {
          rval = -1;
        } else if (field1->field->choice > field2->field->choice) {
          rval = 1;
        } else {
          switch (field1->field->choice) {
            case FeatQualChoice_legal_qual:
              if (field1->field->data.intvalue < field2->field->data.intvalue) {
                rval = -1;
              } else if (field1->field->data.intvalue > field2->field->data.intvalue) {
                rval = 1;
              }
              break;
            case FeatQualChoice_illegal_qual:
              rval = 0;
              break;
          }
        }
        break;
      case FieldType_cds_gene_prot:
      case FieldType_pub:
      case FieldType_misc:
        if (vnp1->data.intvalue > vnp2->data.intvalue) {
          rval = 1;
        } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
          rval = -1;
        }
        break;
      case FieldType_rna_field:
        rq1 = (RnaQualPtr) vnp1->data.ptrvalue;
        rq2 = (RnaQualPtr) vnp2->data.ptrvalue;
        if (rq1 == NULL && rq2 == NULL) {
          rval = 0;
        } else if (rq1 == NULL) {
          rval = -1;
        } else if (rq2 == NULL) {
          rval = 1;
        } else if ((rval = CompareRnaTypes (rq1->type, rq2->type)) == 0) {
          if (rq1->field < rq2->field) {
            rval = -1;
          } else if (rq1->field > rq2->field) {
            rval = 1;
          } else {
            rval = 0;
          }
        }
        break;
      case FieldType_struc_comment_field:
        scf1 = (StructuredCommentFieldPtr) vnp1->data.ptrvalue;
        scf2 = (StructuredCommentFieldPtr) vnp2->data.ptrvalue;
        if (scf1 == NULL && scf2 == NULL) {
          rval = 0;
        } else if (scf1 == NULL) {
          rval = -1;
        } else if (scf2 == NULL) {
          rval = 1;
        } else if (scf1->choice < scf2->choice) {
          rval = -1;
        } else if (scf1->choice > scf2->choice) {
          rval = 1;
        } else if (scf1->choice == StructuredCommentField_named) {
          rval = StringCmp (scf1->data.ptrvalue, scf2->data.ptrvalue);
        }
        break;
    }
  }
  return rval;
}

static Boolean DoFieldTypesMatch (FieldTypePtr field1, FieldTypePtr field2)
{
  if (CompareFieldTypes (field1, field2) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Int2 FeatureTypeFromCDSGeneProtField (Uint2 cds_gene_prot_field);


NLM_EXTERN Int2 FeatureTypeFromFieldType (FieldTypePtr field)
{
  Int2 feat_type = Feature_type_any;
  FeatureFieldPtr ffp;
  RnaQualPtr      rq;

  if (field == NULL) {
    feat_type = Feature_type_any;
  } else {
    switch (field->choice) {
      case FieldType_source_qual:
        feat_type = Feature_type_biosrc;
        break;
      case FieldType_feature_field:
        ffp = (FeatureFieldPtr) field->data.ptrvalue;
        if (ffp != NULL) {
          feat_type = ffp->type;
        }
        break;
      case FieldType_rna_field:
        rq = (RnaQualPtr) field->data.ptrvalue;
        if (rq != NULL) {
          feat_type = GetFeatureTypeForRnaType (rq->type->choice);
        }
        break;
      case FieldType_cds_gene_prot:
        feat_type = FeatureTypeFromCDSGeneProtField (field->data.intvalue);
        break;
    }
  }
  return feat_type;
}


NLM_EXTERN Boolean IsFeatureFieldEmpty (FeatureFieldPtr field)
{
  if (field == NULL) return TRUE;
  if (field->field == NULL) return TRUE;
  return FALSE;
}


NLM_EXTERN Boolean IsRnaQualEmpty (RnaQualPtr rq)
{
  if (rq == NULL) return TRUE;
  return FALSE;
}


NLM_EXTERN Boolean IsFieldTypeEmpty (FieldTypePtr field)
{
  Boolean rval = TRUE;
  ValNodePtr vnp;

  if (field == NULL) return TRUE;
  switch (field->choice) {
    case FieldType_source_qual:
      if (field->data.ptrvalue != NULL) {
        rval = FALSE;
      }
      break;
    case FieldType_feature_field:
      if (!IsFeatureFieldEmpty (field->data.ptrvalue)) {
        rval = FALSE;
      }
      break;
    case FieldType_cds_gene_prot:
      rval = FALSE;
      break;
    case FieldType_pub:
      rval = FALSE;
      break;
    case FieldType_rna_field:
      rval = IsRnaQualEmpty (field->data.ptrvalue);
      break;
    case FieldType_struc_comment_field:
      vnp = field->data.ptrvalue;
      if (vnp == NULL 
          || (vnp->choice == StructuredCommentField_named && StringHasNoText (vnp->data.ptrvalue))
          || (vnp->choice != StructuredCommentField_named && vnp->choice != StructuredCommentField_database)) {
        rval = TRUE;
      } else {
        rval = FALSE;
      }
      break;
    case FieldType_misc:
      rval = FALSE;
      break;      
  }
  return rval;
}

NLM_EXTERN Boolean AllowFieldMulti (FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (field == NULL) return FALSE;
  switch (field->choice) {
    case FieldType_source_qual:
      rval = AllowSourceQualMulti (field->data.ptrvalue);
      break;
    case FieldType_feature_field:
      break;
    case FieldType_cds_gene_prot:
      break;
    case FieldType_pub:
      break;
    case FieldType_rna_field:
      break;
    case FieldType_struc_comment_field:
      break;
    case FieldType_misc:
      break;
  }
  return rval;
}


static Boolean IsUserObjectStructuredComment (UserObjectPtr uop)
{
  if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "StructuredComment") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}
      

static Boolean IsObjectAppropriateForFieldValue (Uint1 choice, Pointer data, FieldTypePtr field)
{
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  FeatureFieldPtr   fp;
  RnaQualPtr        rq;
  Boolean rval = FALSE;

  if (data == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        if (sfp->data.choice == SEQFEAT_BIOSRC) {
          rval = TRUE;
        }
      } else if (choice == OBJ_SEQDESC) {
        sdp = (SeqDescrPtr) data;
        if (sdp->choice == Seq_descr_source) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        fp = (FeatureFieldPtr) field->data.ptrvalue;
        if (fp != NULL && (fp->type == Feature_type_any || GetFeatdefFromFeatureType (fp->type) == sfp->idx.subtype)) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_rna_field :
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        rq = (RnaQualPtr) field->data.ptrvalue;
        if (rq != NULL && DoesFeatureMatchRnaType (sfp, rq->type)) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_cds_gene_prot :
      if (choice == 0) {
        rval = TRUE;
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        rval = TRUE;
      }
      break;
    case FieldType_pub:
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        if (sfp->data.choice == SEQFEAT_PUB) {
          rval = TRUE;
        }
      } else if (choice == OBJ_SEQDESC) {
        sdp = (SeqDescrPtr) data;
        if (sdp->choice == Seq_descr_pub) {
          rval = TRUE;
        }
      }
      break; 
    case FieldType_struc_comment_field:
      if (choice == OBJ_SEQDESC) {
        sdp = (SeqDescrPtr) data;
        if (sdp->choice == Seq_descr_user && IsUserObjectStructuredComment (sdp->data.ptrvalue)) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_misc:
      if (choice == OBJ_BIOSEQ && field->data.intvalue == Misc_field_genome_project_id) {
        rval = TRUE;
      } else if (choice == OBJ_SEQDESC && field->data.intvalue == Misc_field_comment_descriptor) {
        rval = TRUE;
      }
      break;      
  }
  return rval;
}


static Boolean IsObjectAppropriateForFieldPair (Uint1 choice, Pointer data, FieldPairTypePtr fieldpair)
{
  FieldTypePtr f;
  Boolean rval;

  f = GetFromFieldFromFieldPair(fieldpair);
  rval = IsObjectAppropriateForFieldValue(choice, data, f);
  f = FieldTypeFree (f);
  return rval;
}


/* structure and create/free functions for CGPSet, used for handling CDS-Gene-Prot sets */
typedef struct cgpset 
{
  ValNodePtr cds_list;
  ValNodePtr gene_list;
  ValNodePtr prot_list;
  ValNodePtr mrna_list;
} CGPSetData, PNTR CGPSetPtr;



static CGPSetPtr CGPSetNew (void)
{
  CGPSetPtr c;

  c = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  c->cds_list = NULL;
  c->gene_list = NULL;
  c->prot_list = NULL;
  c->mrna_list = NULL;
  return c;
}


static CGPSetPtr CGPSetFree (CGPSetPtr c)
{
  if (c != NULL) {
    c->cds_list = ValNodeFree (c->cds_list);
    c->gene_list = ValNodeFree (c->gene_list);
    c->prot_list = ValNodeFree (c->prot_list);
    c->mrna_list = ValNodeFree (c->mrna_list);
    c = MemFree (c);
  }
  return c;
}


static ValNodePtr FreeCGPSetList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  
  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = CGPSetFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return NULL;
}

static CGPSetPtr BuildCGPSetFromCodingRegion (SeqFeatPtr cds, BoolPtr indexing_needed);
static CGPSetPtr BuildCGPSetFromGene (SeqFeatPtr gene);
static CGPSetPtr BuildCGPSetFrommRNA (SeqFeatPtr mrna);


/* generic functions for mapping fields */

typedef struct feattypefeatdef {
  Int4 feattype;
  Int4 featdef;
  CharPtr featname;
} FeatTypeFeatDefData, PNTR FeatTypeFeatDefPtr;

static FeatTypeFeatDefData feattype_featdef[] = {
 { Feature_type_any , FEATDEF_ANY , "any" } , 
 { Feature_type_gene , FEATDEF_GENE , "gene" } , 
 { Feature_type_org , FEATDEF_ORG , "org" } , 
 { Feature_type_cds , FEATDEF_CDS , "CDS" } , 
 { Feature_type_prot , FEATDEF_PROT , "Protein" } , 
 { Feature_type_preRNA , FEATDEF_preRNA , "preRNA" } , 
 { Feature_type_mRNA , FEATDEF_mRNA , "mRNA" } , 
 { Feature_type_tRNA , FEATDEF_tRNA , "tRNA" } , 
 { Feature_type_rRNA , FEATDEF_rRNA , "rRNA" } , 
 { Feature_type_snRNA , FEATDEF_snRNA , "snRNA" } , 
 { Feature_type_scRNA , FEATDEF_scRNA , "scRNA" } , 
 { Feature_type_otherRNA , FEATDEF_otherRNA , "misc_RNA" } , 
 { Feature_type_pub , FEATDEF_PUB , "pub" } , 
 { Feature_type_seq , FEATDEF_SEQ , "seq" } , 
 { Feature_type_imp , FEATDEF_IMP , "imp" } , 
 { Feature_type_allele , FEATDEF_allele , "allele" } , 
 { Feature_type_attenuator , FEATDEF_attenuator , "attenuator" } , 
 { Feature_type_c_region , FEATDEF_C_region , "c_region" } , 
 { Feature_type_caat_signal , FEATDEF_CAAT_signal , "caat_signal" } , 
 { Feature_type_imp_CDS , FEATDEF_Imp_CDS , "imp_CDS" } , 
 { Feature_type_conflict , FEATDEF_conflict , "conflict" } , 
 { Feature_type_d_loop , FEATDEF_D_loop , "d_loop" } , 
 { Feature_type_d_segment , FEATDEF_D_segment , "d_segment" } , 
 { Feature_type_enhancer , FEATDEF_enhancer , "enhancer" } , 
 { Feature_type_exon , FEATDEF_exon , "exon" } , 
 { Feature_type_gC_signal , FEATDEF_GC_signal , "gC_signal" } , 
 { Feature_type_iDNA , FEATDEF_iDNA , "iDNA" } , 
 { Feature_type_intron , FEATDEF_intron , "intron" } , 
 { Feature_type_j_segment , FEATDEF_J_segment , "j_segment" } , 
 { Feature_type_ltr , FEATDEF_LTR , "ltr" } , 
 { Feature_type_mat_peptide , FEATDEF_mat_peptide , "mat_peptide" } , 
 { Feature_type_misc_binding , FEATDEF_misc_binding , "misc_binding" } , 
 { Feature_type_misc_difference , FEATDEF_misc_difference , "misc_difference" } , 
 { Feature_type_misc_feature , FEATDEF_misc_feature , "misc_feature" } , 
 { Feature_type_misc_recomb , FEATDEF_misc_recomb , "misc_recomb" } , 
 { Feature_type_misc_RNA , FEATDEF_otherRNA , "misc_RNA" } , 
 { Feature_type_misc_signal , FEATDEF_misc_signal , "misc_signal" } , 
 { Feature_type_misc_structure , FEATDEF_misc_structure , "misc_structure" } , 
 { Feature_type_modified_base , FEATDEF_modified_base , "modified_base" } , 
 { Feature_type_mutation , FEATDEF_mutation , "mutation" } , 
 { Feature_type_n_region , FEATDEF_N_region , "n_region" } , 
 { Feature_type_old_sequence , FEATDEF_old_sequence , "old_sequence" } , 
 { Feature_type_polyA_signal , FEATDEF_polyA_signal , "polyA_signal" } , 
 { Feature_type_polyA_site , FEATDEF_polyA_site , "polyA_site" } , 
 { Feature_type_precursor_RNA , FEATDEF_precursor_RNA , "precursor_RNA" } , 
 { Feature_type_prim_transcript , FEATDEF_prim_transcript , "prim_transcript" } , 
 { Feature_type_primer_bind , FEATDEF_primer_bind , "primer_bind" } , 
 { Feature_type_promoter , FEATDEF_promoter , "promoter" } , 
 { Feature_type_protein_bind , FEATDEF_protein_bind , "protein_bind" } , 
 { Feature_type_rbs , FEATDEF_RBS , "rbs" } , 
 { Feature_type_repeat_region , FEATDEF_repeat_region , "repeat_region" } , 
 { Feature_type_rep_origin , FEATDEF_rep_origin , "rep_origin" } , 
 { Feature_type_s_region , FEATDEF_S_region , "s_region" } , 
 { Feature_type_sig_peptide , FEATDEF_sig_peptide , "sig_peptide" } , 
 { Feature_type_source , FEATDEF_source , "source" } , 
 { Feature_type_stem_loop , FEATDEF_stem_loop , "stem_loop" } , 
 { Feature_type_sts , FEATDEF_STS , "sts" } , 
 { Feature_type_tata_signal , FEATDEF_TATA_signal , "tata_signal" } , 
 { Feature_type_terminator , FEATDEF_terminator , "terminator" } , 
 { Feature_type_transit_peptide , FEATDEF_transit_peptide , "transit_peptide" } , 
 { Feature_type_unsure , FEATDEF_unsure , "unsure" } , 
 { Feature_type_v_region , FEATDEF_V_region , "v_region" } , 
 { Feature_type_v_segment , FEATDEF_V_segment , "v_segment" } , 
 { Feature_type_variation , FEATDEF_variation , "variation" } , 
 { Feature_type_virion , FEATDEF_virion , "virion" } , 
 { Feature_type_n3clip , FEATDEF_3clip , "3clip" } , 
 { Feature_type_n3UTR , FEATDEF_3UTR , "3UTR" } , 
 { Feature_type_n5clip , FEATDEF_5clip , "5clip" } , 
 { Feature_type_n5UTR , FEATDEF_5UTR , "5UTR" } , 
 { Feature_type_n10_signal , FEATDEF_10_signal , "10_signal" } , 
 { Feature_type_n35_signal , FEATDEF_35_signal , "35_signal" } , 
 { Feature_type_site_ref , FEATDEF_site_ref , "site_ref" } , 
 { Feature_type_region , FEATDEF_REGION , "region" } , 
 { Feature_type_comment , FEATDEF_COMMENT , "comment" } , 
 { Feature_type_bond , FEATDEF_BOND , "bond" } , 
 { Feature_type_site , FEATDEF_SITE , "site" } , 
 { Feature_type_rsite , FEATDEF_RSITE , "rsite" } , 
 { Feature_type_user , FEATDEF_USER , "user" } , 
 { Feature_type_txinit , FEATDEF_TXINIT , "txinit" } , 
 { Feature_type_num , FEATDEF_NUM , "num" } , 
 { Feature_type_psec_str , FEATDEF_PSEC_STR , "psec_str" } , 
 { Feature_type_non_std_residue , FEATDEF_NON_STD_RESIDUE , "non_std_residue" } , 
 { Feature_type_het , FEATDEF_HET , "het" } , 
 { Feature_type_biosrc , FEATDEF_BIOSRC , "biosrc" } , 
 { Feature_type_preprotein , FEATDEF_preprotein , "preprotein" } , 
 { Feature_type_mat_peptide_aa , FEATDEF_mat_peptide_aa , "mat_peptide_aa" } , 
 { Feature_type_sig_peptide_aa , FEATDEF_sig_peptide_aa , "sig_peptide_aa" } , 
 { Feature_type_transit_peptide_aa , FEATDEF_transit_peptide_aa , "transit_peptide_aa" } , 
 { Feature_type_snoRNA , FEATDEF_snoRNA , "snoRNA" } , 
 { Feature_type_gap , FEATDEF_gap , "gap" } , 
 { Feature_type_operon , FEATDEF_operon , "operon" } , 
 { Feature_type_oriT , FEATDEF_oriT , "oriT" } , 
 { Feature_type_ncRNA , FEATDEF_ncRNA , "ncRNA" } , 
 { Feature_type_tmRNA , FEATDEF_tmRNA , "tmRNA" }};

#define NUM_feattype_featdef sizeof (feattype_featdef) / sizeof (FeatTypeFeatDefData)

NLM_EXTERN Int4 GetFeatdefFromFeatureType (Int4 feature_type) 
{
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef; i++) {
    if (feature_type == feattype_featdef[i].feattype) {
      return feattype_featdef[i].featdef;
    }
  }
  return FEATDEF_BAD;
}


NLM_EXTERN Int4 GetFeatureTypeFromFeatdef (Int4 featdef) 
{
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef; i++) {
    if (featdef == feattype_featdef[i].featdef) {
      return feattype_featdef[i].feattype;
    }
  }
  return FEATDEF_BAD;
}


NLM_EXTERN CharPtr GetFeatureNameFromFeatureType (Int4 feature_type)
{
  CharPtr str = NULL;
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef && str == NULL; i++) {
    if (feature_type == feattype_featdef[i].feattype) {
      str = feattype_featdef[feature_type].featname;
    }
  } 
  if (str == NULL) {
    str = "Unknown feature type";
  }
  return str;
}


static Boolean Matchnamestring (CharPtr name1, CharPtr name2)
{
  if (name1 == NULL && name2 == NULL) {
    return TRUE;
  } else if (name1 == NULL || name2 == NULL) {
    return FALSE;
  } else {
    while (*name1 != 0 && *name2 != 0) {
      while (*name1 == ' ' || *name1 == '-' || *name1 == '_') {
        name1++;
      }
      while (*name2 == ' ' || *name2 == '-' || *name2 == '_') {
        name2++;
      }
      if (tolower (*name1) != tolower(*name2)) {
        return FALSE;
      }
      name1++;
      name2++;
    }
    if (*name1 == 0 && *name2 == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  }
}


typedef struct stringalias {
  CharPtr alias;
  CharPtr canonical;
} StringAliasData, PNTR StringAliasPtr;


static CharPtr GetCanonical (CharPtr str, StringAliasPtr alias_list)
{
  Int4 i;

  if (alias_list == NULL) {
    return str;
  }
  for (i = 0; alias_list[i].alias != NULL; i++) {
    if (Matchnamestring (str, alias_list[i].alias)) {
      return alias_list[i].canonical;
    }
  }
  return str;
}


NLM_EXTERN Int4 GetFeatureTypeByName (CharPtr feat_name)
{
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef; i++) {
    if (Matchnamestring (feattype_featdef[i].featname, feat_name)) {
      return feattype_featdef[i].feattype;
    }
  }
  return -1;  
}


NLM_EXTERN void AddImportFeaturesToChoiceList (ValNodePtr PNTR feature_type_list)
{
  Int4 i, seqfeattype;
  CharPtr featname;
  ValNodePtr tmp_list = NULL;

  for (i = 1; i < NUM_feattype_featdef; i++) {
    if (feattype_featdef[i].feattype == Feature_type_gap) continue;
    seqfeattype = FindFeatFromFeatDefType (feattype_featdef[i].featdef);
    if (seqfeattype == SEQFEAT_IMP) {
      featname = GetFeatureNameFromFeatureType (feattype_featdef[i].feattype);
      if (featname != NULL) {
        ValNodeAddPointer (&tmp_list, feattype_featdef[i].feattype, StringSave (featname));
      }
    }
  }
  tmp_list = ValNodeSort (tmp_list, SortVnpByString);
  ValNodeLink (feature_type_list, tmp_list);
}



static Boolean IsMostUsedFeature (Uint1 val)
{
  if (val == Feature_type_gene
      || val == Feature_type_cds
      || val == Feature_type_prot
      || val == Feature_type_exon
      || val == Feature_type_intron
      || val == Feature_type_mRNA
      || val == Feature_type_rRNA
      || val == Feature_type_otherRNA) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static int LIBCALLBACK SortVnpByFeatureName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  Boolean     most_used1, most_used2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      most_used1 = IsMostUsedFeature (vnp1->choice);
      most_used2 = IsMostUsedFeature (vnp2->choice);
      if (most_used1 && !most_used2) {
        return -1;
      } else if (!most_used1 && most_used2) {
        return 1;
      } else {
        str1 = (CharPtr) vnp1->data.ptrvalue;
        str2 = (CharPtr) vnp2->data.ptrvalue;
        if (str1 != NULL && str2 != NULL) {
          return StringICmp (str1, str2);
        }
      }
    }
  }
  return 0;
}


NLM_EXTERN void AddAllFeaturesToChoiceList (ValNodePtr PNTR feature_type_list)
{
  Int4 i;
  CharPtr featname;
  ValNodePtr tmp_list = NULL;

  for (i = 1; i < NUM_feattype_featdef; i++) {
    if (feattype_featdef[i].feattype == Feature_type_gap) continue;
    featname = GetFeatureNameFromFeatureType (feattype_featdef[i].feattype);
    if (featname != NULL) {
      ValNodeAddPointer (&tmp_list, feattype_featdef[i].feattype, StringSave (featname));
    }
  }
  tmp_list = ValNodeSort (tmp_list, SortVnpByFeatureName);
  ValNodeLink (feature_type_list, tmp_list);
}


typedef struct featqualgbqual {
  Int4 featqual;
  Int4 gbqual;
  Int4 subfield;
  CharPtr qualname;
} FeatQualGBQualData, PNTR FeatQualGBQualPtr;

static FeatQualGBQualData featqual_gbqual[] = {
 { Feat_qual_legal_allele , GBQUAL_allele , 0,  "allele" } , 
 { Feat_qual_legal_anticodon , GBQUAL_anticodon , 0,  "anticodon" } , 
 { Feat_qual_legal_bound_moiety , GBQUAL_bound_moiety , 0,  "bound-moiety" } , 
 { Feat_qual_legal_chromosome , GBQUAL_chromosome , 0, "chromosome" } , 
 { Feat_qual_legal_citation , GBQUAL_citation , 0, "citation" } , 
 { Feat_qual_legal_codon , GBQUAL_codon , 0, "codon" } , 
 { Feat_qual_legal_codon_start , GBQUAL_codon_start , 0, "codon-start" } , 
 { Feat_qual_legal_compare , GBQUAL_compare , 0, "compare" } , 
 { Feat_qual_legal_cons_splice , GBQUAL_cons_splice , 0, "cons-splice" } , 
 { Feat_qual_legal_db_xref , GBQUAL_db_xref , 0, "db-xref" } , 
 { Feat_qual_legal_direction , GBQUAL_direction , 0, "direction" } , 
 { Feat_qual_legal_environmental_sample , 0, GBQUAL_environmental_sample , "environmental-sample" } , 
 { Feat_qual_legal_evidence , GBQUAL_evidence , 0, "evidence" } , 
 { Feat_qual_legal_exception , GBQUAL_exception , 0, "exception" } , 
 { Feat_qual_legal_experiment , GBQUAL_experiment , 0, "experiment" } , 
 { Feat_qual_legal_focus , GBQUAL_focus , 0, "focus" } , 
 { Feat_qual_legal_frequency , GBQUAL_frequency , 0, "frequency" } , 
 { Feat_qual_legal_function , GBQUAL_function , 0, "function" } , 
 { Feat_qual_legal_gene , GBQUAL_gene , 0, "locus" } , 
 { Feat_qual_legal_inference , GBQUAL_inference , 0, "inference" } , 
 { Feat_qual_legal_label , GBQUAL_label , 0, "label" } , 
 { Feat_qual_legal_locus_tag , GBQUAL_locus_tag , 0, "locus-tag" } , 
 { Feat_qual_legal_map , GBQUAL_map , 0, "map" } , 
 { Feat_qual_legal_mobile_element , GBQUAL_mobile_element , 0, "mobile-element" } , 
 { Feat_qual_legal_mobile_element_type , GBQUAL_mobile_element , 1, "mobile-element-type"} ,
 { Feat_qual_legal_mobile_element_name , GBQUAL_mobile_element , 2, "mobile-element-name"} ,
 { Feat_qual_legal_mod_base , GBQUAL_mod_base , 0, "mod-base" } , 
 { Feat_qual_legal_mol_type , GBQUAL_mol_type , 0, "mol-type" } , 
 { Feat_qual_legal_ncRNA_class , GBQUAL_ncRNA_class , 0, "ncRNA-class" } , 
 { Feat_qual_legal_note , GBQUAL_note , 0, "note" } , 
 { Feat_qual_legal_number , GBQUAL_number , 0, "number" } , 
 { Feat_qual_legal_old_locus_tag , GBQUAL_old_locus_tag , 0, "old-locus-tag" } , 
 { Feat_qual_legal_operon , GBQUAL_operon , 0, "operon" } , 
 { Feat_qual_legal_organism , GBQUAL_organism , 0, "organism" } , 
 { Feat_qual_legal_organelle , GBQUAL_organelle , 0, "organelle" } , 
 { Feat_qual_legal_partial , GBQUAL_partial , 0, "partial" } , 
 { Feat_qual_legal_phenotype , GBQUAL_phenotype , 0, "phenotype" } , 
 { Feat_qual_legal_plasmid , GBQUAL_plasmid , 0, "plasmid" } , 
 { Feat_qual_legal_product , GBQUAL_product , 0, "product" } , 
 { Feat_qual_legal_protein_id , GBQUAL_protein_id , 0, "protein-id" } , 
 { Feat_qual_legal_pseudo , GBQUAL_pseudo , 0, "pseudo" } , 
 { Feat_qual_legal_rearranged , GBQUAL_rearranged , 0, "rearranged" } , 
 { Feat_qual_legal_replace , GBQUAL_replace , 0, "replace" } , 
 { Feat_qual_legal_rpt_family , GBQUAL_rpt_family , 0, "rpt-family" } , 
 { Feat_qual_legal_rpt_type , GBQUAL_rpt_type , 0, "rpt-type" } , 
 { Feat_qual_legal_rpt_unit , GBQUAL_rpt_unit , 0, "rpt-unit" } , 
 { Feat_qual_legal_rpt_unit_seq , GBQUAL_rpt_unit_seq , 0, "rpt-unit-seq" } , 
 { Feat_qual_legal_rpt_unit_range , GBQUAL_rpt_unit_range , 0, "rpt-unit-range" } , 
 { Feat_qual_legal_satellite , GBQUAL_satellite , 0, "satellite" } , 
 { Feat_qual_legal_satellite_type , GBQUAL_satellite, 1, "satellite-type"} ,
 { Feat_qual_legal_satellite_name , GBQUAL_satellite, 2, "satellite-name"} ,
 { Feat_qual_legal_segment , GBQUAL_segment , 0, "segment" } , 
 { Feat_qual_legal_sequenced_mol , GBQUAL_sequenced_mol , 0, "sequenced-mol" } , 
 { Feat_qual_legal_standard_name , GBQUAL_standard_name , 0, "standard-name" } , 
 { Feat_qual_legal_transcript_id , GBQUAL_transcript_id , 0, "transcript-id" } , 
 { Feat_qual_legal_transgenic , GBQUAL_transgenic , 0, "transgenic" } , 
 { Feat_qual_legal_translation , GBQUAL_translation , 0, "translation" } , 
 { Feat_qual_legal_transl_except , GBQUAL_transl_except , 0, "transl-except" } , 
 { Feat_qual_legal_transl_table , GBQUAL_transl_table , 0, "transl-table" } , 
 { Feat_qual_legal_usedin , GBQUAL_usedin , 0, "usedin" } };

#define NUM_featqual_gbqual sizeof (featqual_gbqual) / sizeof (FeatQualGBQualData)


NLM_EXTERN Int4 GetNumFeatQual (void)
{
  return NUM_featqual_gbqual;
}


static Int4 GetGBQualFromFeatQual (Int4 featqual, Int4Ptr subfield) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual == featqual_gbqual[i].featqual) {
      if (subfield != NULL) {
        *subfield = featqual_gbqual[i].subfield;
      }
      return featqual_gbqual[i].gbqual;
    }
  }
  return -1;
}


static Int4 GetFeatQualByGBQualAndSubfield (Int4 gbqual, Int4 subfield)
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual_gbqual[i].gbqual == gbqual && featqual_gbqual[i].subfield == subfield) {
      return featqual_gbqual[i].featqual;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr GetFeatQualName (Int4 featqual) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual == featqual_gbqual[i].featqual) {
      return featqual_gbqual[i].qualname;
    }
  }
  return NULL;
}


NLM_EXTERN Int4 GetFeatQualByName (CharPtr qualname) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (Matchnamestring (featqual_gbqual[i].qualname, qualname)) {
      return featqual_gbqual[i].featqual;
    }
  }
  return -1;  
}


static Int4 NumGbQualSubfields (Int4 gbqual)
{
  Int4 i, num_subfields = 0;
  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual_gbqual[i].gbqual == gbqual) {
      if (featqual_gbqual[i].subfield > num_subfields) {
        num_subfields = featqual_gbqual[i].subfield;
      }
    }
  }
  return num_subfields;
}


NLM_EXTERN void AddAllFeatureFieldsToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    ValNodeAddPointer (field_list, featqual_gbqual[i].featqual, StringSave (featqual_gbqual[i].qualname));
  }
}


NLM_EXTERN CharPtr SummarizeFeatQual (ValNodePtr qual)
{
  if (qual == NULL) {
    return StringSave ("unspecified qualifier");
  } else if (qual->choice == FeatQualChoice_legal_qual) {
    return StringSave (GetFeatQualName (qual->data.intvalue));
  } else if (qual->choice == FeatQualChoice_illegal_qual) {
    return StringSave (qual->data.ptrvalue);
  } else {
    return StringSave ("unspecified qualifier");
  }
}


/* functions for RnaQual values */

/* functions for RnaType values */
typedef struct rnatypemap {
  Int4 rnatype;
  Int4 rnaval;
  Int4 featuretype;
  CharPtr rnaname;
} RnaTypeMapData, PNTR RnaTypeMapPtr;

static RnaTypeMapData rnatypemap[] = {
 { RnaFeatType_preRNA , RNA_TYPE_premsg, Feature_type_preRNA, "preRNA" } ,
 { RnaFeatType_mRNA , RNA_TYPE_mRNA, Feature_type_mRNA, "mRNA" } ,
 { RnaFeatType_tRNA , RNA_TYPE_tRNA, Feature_type_tRNA, "tRNA" } ,
 { RnaFeatType_rRNA , RNA_TYPE_rRNA, Feature_type_rRNA, "rRNA" } ,
 { RnaFeatType_ncRNA , 255 , Feature_type_ncRNA, "ncRNA" } ,
 { RnaFeatType_tmRNA , 255 , Feature_type_tmRNA, "tmRNA" } ,
 { RnaFeatType_miscRNA , 255 , Feature_type_misc_RNA, "misc_RNA" }
};

#define NUM_rnatypemap sizeof (rnatypemap) / sizeof (RnaTypeMapData)


static CharPtr GetNameForRnaType (Int4 rnatype)
{
  Int4 i;

  for (i = 0; i < NUM_rnatypemap; i++) {
    if (rnatypemap[i].rnatype == rnatype) {
      return rnatypemap[i].rnaname;
    }
  }
  return NULL;
}


static Int4 GetRnaTypeForName (CharPtr rnaname)
{
  Int4 i;

  for (i = 0; i < NUM_rnatypemap; i++) {
    if (StringCmp (rnatypemap[i].rnaname, rnaname) == 0) {
      return rnatypemap[i].rnatype;
    }
  }
  return -1;
}


static Int4 GetRnaValForRnaType (Int4 rnatype)
{
  Int4 i;

  for (i = 0; i < NUM_rnatypemap; i++) {
    if (rnatypemap[i].rnatype == rnatype) {
      return rnatypemap[i].rnaval;
    }
  }
  return -1;
}


NLM_EXTERN Int4 GetFeatureTypeForRnaType (Int4 rnatype)
{
  Int4 i;

  for (i = 0; i < NUM_rnatypemap; i++) {
    if (rnatypemap[i].rnatype == rnatype) {
      return rnatypemap[i].featuretype;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetRNATypeList (void)
{
  Int4 i;
  ValNodePtr list = NULL;

  for (i = 1; i < NUM_rnatypemap; i++) {
    ValNodeAddPointer (&list, rnatypemap[i].rnatype, StringSave (rnatypemap[i].rnaname));
  }
  return list;
}


static Boolean DoesFeatureMatchRnaType (SeqFeatPtr sfp, RnaFeatTypePtr rt)
{
  Boolean rval = FALSE;
  RnaRefPtr rrp;
  GBQualPtr gbq;
  Int4 rnaval;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) {
    return FALSE;
  }
  if (rt == NULL) return TRUE;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return FALSE;

  rnaval = GetRnaValForRnaType (rt->choice);
  if (rnaval == rrp->type) {
    switch (rt->choice) {
      case RnaFeatType_ncRNA:
        if (rrp->ext.choice == 1
            && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0) {
          if (rt->data.ptrvalue == NULL) {
            rval = TRUE;
          } else {
            gbq = sfp->qual;
            while (gbq != NULL && StringCmp (gbq->qual, "ncRNA_class") != 0) {
              gbq = gbq->next;
            }
            if (gbq != NULL && StringCmp (gbq->val, rt->data.ptrvalue) == 0) {
              rval = TRUE;
            }
          }
        }
        break;
      case RnaFeatType_tmRNA:
        if (rrp->ext.choice == 1
            && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0) {
          rval = TRUE;
        }
        break;
      case RnaFeatType_miscRNA:
        if (rrp->ext.choice == 1
            && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
            && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0) {
          rval = TRUE;
        }
        break;
      default:
        rval = TRUE;
        break;
    }
  }
  return rval;
}


static Int4 CompareRnaTypes (RnaFeatTypePtr rt1, RnaFeatTypePtr rt2)
{
  Int4 rval = 0;

  if (rt1 == NULL && rt2 == NULL) {
    rval = 0;
  } else if (rt1 == NULL) {
    rval = -1;
  } else if (rt2 == NULL) {
    rval = 1;
  } else if (rt1->choice < rt2->choice) {
    rval = -1;
  } else if (rt1->choice > rt2->choice) {
    rval = 1;
  } else if (rt1->choice == RnaFeatType_ncRNA) {
    rval = StringCmp (rt1->data.ptrvalue, rt2->data.ptrvalue);
  } else {
    rval = 0;
  }
  return rval;
}


static RnaFeatTypePtr RnaFeatTypeFromSeqFeat (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  RnaFeatTypePtr rt = NULL;
  GBQualPtr gbqual;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  switch (rrp->type) {
    case RNA_TYPE_premsg:
      rt = ValNodeNew (NULL);
      rt->choice = RnaFeatType_preRNA;
      break;
    case RNA_TYPE_mRNA:
      rt = ValNodeNew (NULL);
      rt->choice = RnaFeatType_mRNA;
      break;
    case RNA_TYPE_tRNA:
      rt = ValNodeNew (NULL);
      rt->choice = RnaFeatType_tRNA;
      break;
    case RNA_TYPE_rRNA:
      rt = ValNodeNew (NULL);
      rt->choice = RnaFeatType_rRNA;
      break;
    case 255:
      if (rrp->ext.choice == 1) {
        if (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0) {
          rt = ValNodeNew (NULL);
          rt->choice = RnaFeatType_ncRNA;
          gbqual = sfp->qual;
          while (gbqual != NULL && StringCmp (gbqual->qual, "ncRNA_class") != 0) {
            gbqual = gbqual->next;
          }
          if (gbqual != NULL) {
            rt->data.ptrvalue = StringSave (gbqual->val);
          }
        } else if (StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0) {
          rt = ValNodeNew (NULL);
          rt->choice = RnaFeatType_tmRNA;
        } else {
          rt = ValNodeNew (NULL);
          rt->choice = RnaFeatType_miscRNA;
        }
      }
      break;
  }
  return rt;
}


typedef struct rnafieldname {
  Int4 field;
  Int4    featqual;
  CharPtr fieldname;
} RnaFieldNameData, PNTR RnaFieldNamePtr;

static RnaFieldNameData rnafieldnames[] = {
 { Rna_field_product , Feat_qual_legal_product, "product" } ,
 { Rna_field_comment , Feat_qual_legal_note, "comment" } ,
 { Rna_field_codons_recognized , Feat_qual_legal_codons_recognized, "codons recognized" } ,
 { Rna_field_ncrna_class , Feat_qual_legal_ncRNA_class, "ncRNA class" } ,
 { Rna_field_anticodon , Feat_qual_legal_anticodon, "anticodon" } ,
 { Rna_field_transcript_id , Feat_qual_legal_transcript_id, "transcript ID" } ,
 { Rna_field_gene_locus , Feat_qual_legal_gene, "gene locus" } ,
 { Rna_field_gene_description , Feat_qual_legal_gene_description, "gene description" } ,
 { Rna_field_gene_maploc , Feat_qual_legal_map, "gene maploc" } ,
 { Rna_field_gene_locus_tag , Feat_qual_legal_locus_tag, "gene locus tag" } ,
 { Rna_field_gene_synonym , Feat_qual_legal_synonym, "gene synonym" } ,
 { Rna_field_gene_comment , Feat_qual_legal_gene_comment, "gene comment" }
};

#define NUM_rnafieldnames sizeof (rnafieldnames) / sizeof (RnaFieldNameData)


NLM_EXTERN CharPtr GetNameForRnaField (Int4 rnafield) 
{
  Int4 i;

  for (i = 0; i < NUM_rnafieldnames; i++) {
    if (rnafieldnames[i].field == rnafield) {
      return rnafieldnames[i].fieldname;
    }
  }
  return NULL;
}


static Int4 GetRnaFieldForName (CharPtr fieldname)
{
  Int4 i;

  for (i = 0; i < NUM_rnafieldnames; i++) {
    if (StringCmp (rnafieldnames[i].fieldname, fieldname) == 0) {
      return rnafieldnames[i].field;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetRnaFieldList (void)
{
  ValNodePtr list = NULL;
  Int4       i;

  for (i = 0; i < NUM_rnafieldnames; i++) {
    ValNodeAddPointer (&list, rnafieldnames[i].field, StringSave (rnafieldnames[i].fieldname));
  }
  return list;
}
  

static Int4 GetFeatQualForRnaField (Int4 field)
{
  Int4 i;

  for (i = 0; i < NUM_rnafieldnames; i++) {
    if (rnafieldnames[i].field == field) {
      return rnafieldnames[i].featqual;
    }
  }
  return -1;
}


NLM_EXTERN FeatureFieldPtr FeatureFieldFromRnaQual (RnaQualPtr rq)
{
  FeatureFieldPtr ffp = NULL;
  Int4 type, qual;

  if (rq == NULL || rq->type == NULL) return NULL;

  type = GetFeatureTypeForRnaType (rq->type->choice);
  qual = GetFeatQualForRnaField (rq->field);
  if (type >= 0 && qual >= 0) {
    ffp = FeatureFieldNew ();
    ffp->type = type;
    ValNodeAddInt (&(ffp->field), FeatQualChoice_legal_qual, qual);
   } 
  return ffp;
}


NLM_EXTERN RnaQualPtr RnaQualFromFeatureField (FeatureFieldPtr ffp)
{
  RnaQualPtr rq = NULL;
  Int4       i;

  if (ffp != NULL && ffp->field != NULL
      && ffp->field->choice == FeatQualChoice_legal_qual) {
    for (i = 0; i < NUM_rnafieldnames && rnafieldnames[i].featqual != ffp->field->choice; i++) {
    }
    if (i < NUM_rnafieldnames) {
      rq = RnaQualNew ();
      rq->field = rnafieldnames[i].featqual;
      rq->type = ValNodeNew (NULL);
      switch (ffp->type) {
        case Feature_type_preRNA:
        case Feature_type_precursor_RNA:
          rq->type->choice = RnaFeatType_preRNA;
          break;
        case Feature_type_mRNA:
          rq->type->choice = RnaFeatType_mRNA;
          break;
        case Feature_type_tRNA:
          rq->type->choice = RnaFeatType_tRNA;
          break;
        case Feature_type_rRNA:
          rq->type->choice = RnaFeatType_rRNA;
          break;
        case Feature_type_snRNA:
          rq->type->choice = RnaFeatType_ncRNA;
          rq->type->data.ptrvalue = StringSave ("snRNA");
          break;
        case Feature_type_scRNA:
          rq->type->choice = RnaFeatType_ncRNA;
          rq->type->data.ptrvalue = StringSave ("scRNA");
          break;
        case Feature_type_snoRNA:
          rq->type->choice = RnaFeatType_ncRNA;
          rq->type->data.ptrvalue = StringSave ("snoRNA");
          break;
        case Feature_type_otherRNA:
        case Feature_type_misc_RNA:
          rq->type->choice = RnaFeatType_miscRNA;
          break;
        case Feature_type_ncRNA:
          rq->type->choice = RnaFeatType_ncRNA;
          break;
        case Feature_type_tmRNA:
          rq->type->choice = RnaFeatType_tmRNA;
          break;
        default:
          rq = RnaQualFree (rq);
          break;
      }
    }
  }
  return rq;
}


NLM_EXTERN CharPtr SummarizeRnaType (RnaFeatTypePtr rt)
{
  CharPtr rnatypename = NULL;
  CharPtr fmt = "%s ncRNA";

  if (rt == NULL) {
    rnatypename = StringSave ("Any RNA");
  } else if (rt->choice == RnaFeatType_ncRNA) {
    if (StringHasNoText (rt->data.ptrvalue)) {
      return StringSave ("ncRNA");
    } else {
      rnatypename = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (rt->data.ptrvalue)));
      sprintf (rnatypename, fmt, rt->data.ptrvalue);
    }
  } else {
    rnatypename = StringSave (GetNameForRnaType (rt->choice));
  }
  return rnatypename;      
}


static CharPtr SummarizeRnaQual (RnaQualPtr rq)
{
  CharPtr rnatypename, qualname;
  CharPtr any_fmt = "RNA %s";
  CharPtr fmt = "%s %s";
  CharPtr s = NULL;

  if (rq == NULL) return NULL;

  qualname = GetNameForRnaField (rq->field);
  if (qualname == NULL) {
    return NULL;
  }

  rnatypename = SummarizeRnaType (rq->type);

  if (rnatypename == NULL) {
    s = (CharPtr) MemNew (sizeof (Char) * (StringLen (any_fmt) + StringLen (qualname)));
    sprintf (s, any_fmt, qualname);
  } else {
    s = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (rnatypename) + StringLen (qualname)));
    sprintf (s, fmt, rnatypename, qualname);
    rnatypename = MemFree (rnatypename);
  }
  return s; 
}


static CharPtr SummarizeStructuredCommentField (StructuredCommentFieldPtr field)
{
  CharPtr summ = NULL;

  if (field == NULL) return NULL;

  if (field->choice == StructuredCommentField_database) {
    summ = StringSave ("database");
  } else if (field->choice == StructuredCommentField_named) {
    summ = StringSave (field->data.ptrvalue);
  }
  return summ;
}


#define IS_ORGMOD 1
#define IS_SUBSRC 2
#define IS_OTHER  3

typedef struct srcqualscqual {
  Int4 srcqual;
  Int4 subtype;
  Int4 typeflag;
  Int4 subfield;
  CharPtr qualname;
} SrcQualSCQualData, PNTR SrcQualSCQualPtr;

#define kAllNotesStr "All Notes"
#define kAllQualsStr "All"

static SrcQualSCQualData srcqual_scqual[] = {
 { Source_qual_acronym , ORGMOD_acronym , IS_ORGMOD , 0 , "acronym" } , 
 { Source_qual_anamorph , ORGMOD_anamorph , IS_ORGMOD , 0 , "anamorph" } , 
 { Source_qual_authority , ORGMOD_authority , IS_ORGMOD , 0 , "authority" } , 
 { Source_qual_bio_material , ORGMOD_bio_material , IS_ORGMOD , 0 , "bio-material" } , 
 { Source_qual_bio_material_INST , ORGMOD_bio_material , IS_ORGMOD , 1 , "bio-material-inst" } , 
 { Source_qual_bio_material_COLL , ORGMOD_bio_material , IS_ORGMOD , 2 , "bio-material-coll" } , 
 { Source_qual_bio_material_SpecID , ORGMOD_bio_material , IS_ORGMOD , 3 , "bio-material-specid" } , 
 { Source_qual_biotype , ORGMOD_biotype , IS_ORGMOD , 0 , "biotype" } , 
 { Source_qual_biovar , ORGMOD_biovar , IS_ORGMOD , 0 , "biovar" } , 
 { Source_qual_breed , ORGMOD_breed , IS_ORGMOD , 0 , "breed" } , 
 { Source_qual_cell_line , SUBSRC_cell_line , IS_SUBSRC , 0 , "cell-line" } , 
 { Source_qual_cell_type , SUBSRC_cell_type , IS_SUBSRC , 0 , "cell-type" } , 
 { Source_qual_chemovar , ORGMOD_chemovar , IS_ORGMOD , 0 , "chemovar" } , 
 { Source_qual_chromosome , SUBSRC_chromosome , IS_SUBSRC , 0 , "chromosome" } , 
 { Source_qual_clone , SUBSRC_clone , IS_SUBSRC , 0 , "clone" } , 
 { Source_qual_clone_lib , SUBSRC_clone_lib , IS_SUBSRC , 0 , "clone-lib" } , 
 { Source_qual_collected_by , SUBSRC_collected_by , IS_SUBSRC , 0 , "collected-by" } , 
 { Source_qual_collection_date , SUBSRC_collection_date , IS_SUBSRC , 0 , "collection-date" } , 
 { Source_qual_common , ORGMOD_common , IS_ORGMOD , 0 , "common" } , 
 { Source_qual_common_name , 0 , IS_OTHER , 0 , "common name" } , 
 { Source_qual_country , SUBSRC_country , IS_SUBSRC , 0 , "country" } , 
 { Source_qual_cultivar , ORGMOD_cultivar , IS_ORGMOD , 0 , "cultivar" } , 
 { Source_qual_culture_collection , ORGMOD_culture_collection , IS_ORGMOD , 0 , "culture-collection" } , 
 { Source_qual_culture_collection_INST , ORGMOD_culture_collection , IS_ORGMOD , 1 , "culture-collection-inst" } , 
 { Source_qual_culture_collection_COLL , ORGMOD_culture_collection , IS_ORGMOD , 2 , "culture-collection-coll" } , 
 { Source_qual_culture_collection_SpecID , ORGMOD_culture_collection , IS_ORGMOD , 3 , "culture-collection-specid" } , 
 { Source_qual_dev_stage , SUBSRC_dev_stage , IS_SUBSRC , 0 , "dev-stage" } , 
 { Source_qual_division , 0 , IS_OTHER, 0 , "division" } ,
 { Source_qual_dosage , ORGMOD_dosage , IS_ORGMOD , 0 , "dosage" } , 
 { Source_qual_ecotype , ORGMOD_ecotype , IS_ORGMOD , 0 , "ecotype" } , 
 { Source_qual_endogenous_virus_name , SUBSRC_endogenous_virus_name , IS_SUBSRC , 0 , "endogenous-virus-name" } , 
 { Source_qual_environmental_sample , SUBSRC_environmental_sample , IS_SUBSRC , 0 , "environmental-sample" } , 
 { Source_qual_forma , ORGMOD_forma , IS_ORGMOD , 0 , "forma" } , 
 { Source_qual_forma_specialis , ORGMOD_forma_specialis , IS_ORGMOD , 0 , "forma-specialis" } , 
 { Source_qual_frequency , SUBSRC_frequency , IS_SUBSRC , 0 , "frequency" } , 
 { Source_qual_fwd_primer_name , SUBSRC_fwd_primer_name , IS_SUBSRC , 0 , "fwd-primer-name" } , 
 { Source_qual_fwd_primer_seq , SUBSRC_fwd_primer_seq , IS_SUBSRC , 0 , "fwd-primer-seq" } , 
 { Source_qual_gb_acronym , ORGMOD_gb_acronym , IS_ORGMOD , 0 , "gb-acronym" } , 
 { Source_qual_gb_anamorph , ORGMOD_gb_anamorph , IS_ORGMOD , 0 , "gb-anamorph" } , 
 { Source_qual_gb_synonym , ORGMOD_gb_synonym , IS_ORGMOD , 0 , "gb-synonym" } , 
 { Source_qual_genotype , SUBSRC_genotype , IS_SUBSRC , 0 , "genotype" } , 
 { Source_qual_germline , SUBSRC_germline , IS_SUBSRC , 0 , "germline" } , 
 { Source_qual_group , ORGMOD_group , IS_ORGMOD , 0 , "group" } , 
 { Source_qual_haplotype , SUBSRC_haplotype , IS_SUBSRC , 0 , "haplotype" } , 
 { Source_qual_identified_by , SUBSRC_identified_by , IS_SUBSRC , 0 , "identified-by" } , 
 { Source_qual_insertion_seq_name , SUBSRC_insertion_seq_name , IS_SUBSRC , 0 , "insertion-seq-name" } , 
 { Source_qual_isolate , ORGMOD_isolate , IS_ORGMOD , 0 , "isolate" } , 
 { Source_qual_isolation_source , SUBSRC_isolation_source , IS_SUBSRC , 0 , "isolation-source" } , 
 { Source_qual_lab_host , SUBSRC_lab_host , IS_SUBSRC , 0 , "lab-host" } , 
 { Source_qual_lat_lon , SUBSRC_lat_lon , IS_SUBSRC , 0 , "lat-lon" } , 
 { Source_qual_lineage , 0, IS_OTHER, 0 , "lineage" } ,
 { Source_qual_map , SUBSRC_map , IS_SUBSRC , 0 , "map" } , 
 { Source_qual_metagenome_source , ORGMOD_metagenome_source , IS_ORGMOD , 0 , "metagenome-source" } , 
 { Source_qual_metagenomic , SUBSRC_metagenomic , IS_SUBSRC , 0 , "metagenomic" } , 
 { Source_qual_old_lineage , ORGMOD_old_lineage , IS_ORGMOD , 0 , "old-lineage" } , 
 { Source_qual_old_name , ORGMOD_old_name , IS_ORGMOD , 0 , "old-name" } , 
 { Source_qual_orgmod_note , ORGMOD_other, IS_ORGMOD, 0 , "note-orgmod" } ,
 { Source_qual_pathovar , ORGMOD_pathovar , IS_ORGMOD , 0 , "pathovar" } , 
 { Source_qual_plasmid_name , SUBSRC_plasmid_name , IS_SUBSRC , 0 , "plasmid-name" } , 
 { Source_qual_plastid_name , SUBSRC_plastid_name , IS_SUBSRC , 0 , "plastid-name" } , 
 { Source_qual_pop_variant , SUBSRC_pop_variant , IS_SUBSRC , 0 , "pop-variant" } , 
 { Source_qual_rearranged , SUBSRC_rearranged , IS_SUBSRC , 0 , "rearranged" } , 
 { Source_qual_rev_primer_name , SUBSRC_rev_primer_name , IS_SUBSRC , 0 , "rev-primer-name" } , 
 { Source_qual_rev_primer_seq , SUBSRC_rev_primer_seq , IS_SUBSRC , 0 , "rev-primer-seq" } , 
 { Source_qual_segment , SUBSRC_segment , IS_SUBSRC , 0 , "segment" } , 
 { Source_qual_serogroup , ORGMOD_serogroup , IS_ORGMOD , 0 , "serogroup" } , 
 { Source_qual_serotype , ORGMOD_serotype , IS_ORGMOD , 0 , "serotype" } , 
 { Source_qual_serovar , ORGMOD_serovar , IS_ORGMOD , 0 , "serovar" } , 
 { Source_qual_sex , SUBSRC_sex , IS_SUBSRC , 0 , "sex" } , 
 { Source_qual_nat_host , ORGMOD_nat_host , IS_ORGMOD , 0 , "host" } , 
 { Source_qual_specimen_voucher , ORGMOD_specimen_voucher , IS_ORGMOD , 0 , "specimen-voucher" } , 
 { Source_qual_specimen_voucher_INST , ORGMOD_specimen_voucher , IS_ORGMOD , 1 , "specimen-voucher-inst" } , 
 { Source_qual_specimen_voucher_COLL , ORGMOD_specimen_voucher , IS_ORGMOD , 2 , "specimen-voucher-coll" } , 
 { Source_qual_specimen_voucher_SpecID , ORGMOD_specimen_voucher , IS_ORGMOD , 3 , "specimen-voucher-specid" } , 
 { Source_qual_strain , ORGMOD_strain , IS_ORGMOD , 0 , "strain" } , 
 { Source_qual_subclone , SUBSRC_subclone , IS_SUBSRC , 0 , "subclone" } , 
 { Source_qual_subgroup , ORGMOD_subgroup , IS_ORGMOD , 0 , "subgroup" } , 
 { Source_qual_subsource_note , SUBSRC_other , IS_SUBSRC , 0 , "note-subsrc" } ,
 { Source_qual_sub_species , ORGMOD_sub_species , IS_ORGMOD , 0 , "sub-species" } , 
 { Source_qual_substrain , ORGMOD_substrain , IS_ORGMOD , 0 , "substrain" } , 
 { Source_qual_subtype , ORGMOD_subtype , IS_ORGMOD , 0 , "subtype" } , 
 { Source_qual_synonym , ORGMOD_synonym , IS_ORGMOD , 0 , "synonym" } , 
 { Source_qual_taxname , 0 , IS_OTHER , 0 , "taxname" } , 
 { Source_qual_teleomorph , ORGMOD_teleomorph , IS_ORGMOD , 0 , "teleomorph" } , 
 { Source_qual_tissue_lib , SUBSRC_tissue_lib , IS_SUBSRC , 0 , "tissue-lib" } , 
 { Source_qual_tissue_type , SUBSRC_tissue_type , IS_SUBSRC , 0 , "tissue-type" } , 
 { Source_qual_transgenic , SUBSRC_transgenic , IS_SUBSRC , 0 , "transgenic" } , 
 { Source_qual_transposon_name , SUBSRC_transposon_name , IS_SUBSRC , 0 , "transposon-name" } , 
 { Source_qual_type , ORGMOD_type , IS_ORGMOD , 0 , "type" } , 
 { Source_qual_variety , ORGMOD_variety , IS_ORGMOD , 0 , "variety" } ,
 { Source_qual_all_notes , 255 , IS_OTHER , 0 , kAllNotesStr } ,
 { Source_qual_all_quals , 0 , IS_OTHER , 0, kAllQualsStr } ,
 { Source_qual_mating_type , SUBSRC_mating_type , IS_SUBSRC , 0 , "mating-type" } ,
 { Source_qual_linkage_group , SUBSRC_linkage_group , IS_SUBSRC , 0 , "linkage-group" } ,
 { Source_qual_haplogroup , SUBSRC_haplogroup, IS_SUBSRC, 0, "haplogroup"}
};

#define NUM_srcqual_scqual sizeof (srcqual_scqual) / sizeof (SrcQualSCQualData)

static StringAliasData src_qual_alias_list[] = {
  {"organism", "taxname"},
  {"organism name", "taxname"},
  {"date", "collection-date"},
  {"voucher", "specimen-voucher"},
  {"specific-host", "host"},
  { NULL, NULL}
};


static Int4 GetSubSrcQualFromSrcQual (Int4 srcqual, Int4Ptr subfield) 
{
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual == srcqual_scqual[i].srcqual) {
      if (srcqual_scqual[i].typeflag == IS_SUBSRC) {
        if (subfield != NULL) {
          *subfield = srcqual_scqual[i].subfield;
        }
        return srcqual_scqual[i].subtype;
      } else {
        return -1;
      }
    }
  }
  return -1;
}


static Int4 GetOrgModQualFromSrcQual (Int4 srcqual, Int4Ptr subfield) 
{
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual == srcqual_scqual[i].srcqual) {
      if (srcqual_scqual[i].typeflag == IS_ORGMOD) {
        if (subfield != NULL) {
          *subfield = srcqual_scqual[i].subfield;
        }
        return srcqual_scqual[i].subtype;
      } else {
        return -1;
      }
    }
  }
  return -1;
}


NLM_EXTERN Int4 GetSrcQualFromSubSrcOrOrgMod (Int4 qual, Boolean is_org_mod)
{
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (qual == srcqual_scqual[i].subtype 
        && ((is_org_mod && srcqual_scqual[i].typeflag == IS_ORGMOD)
        || (!is_org_mod && srcqual_scqual[i].typeflag == IS_SUBSRC))) {
      return srcqual_scqual[i].srcqual;
    }
  }
  return -1;
}


NLM_EXTERN Boolean IsNonTextSourceQual (Int4 srcqual)
{
  if (srcqual == Source_qual_transgenic
      || srcqual == Source_qual_germline
      || srcqual == Source_qual_metagenomic
      || srcqual == Source_qual_environmental_sample
      || srcqual == Source_qual_rearranged)
  {
    return TRUE;  
  }
  else
  {
    return FALSE;
  }
}


NLM_EXTERN Boolean IsNonTextFieldType (FieldTypePtr field)
{
  ValNodePtr vnp;

  if (field == NULL) {
    return FALSE;
  } else if (field->choice != FieldType_source_qual) {
    return FALSE;
  } else if ((vnp = field->data.ptrvalue) == NULL) {
    return FALSE;
  } else if (vnp->choice != SourceQualChoice_textqual) {
    return FALSE;
  } else {
    return IsNonTextSourceQual (vnp->data.intvalue);
  }
}


NLM_EXTERN CharPtr GetSourceQualName (Int4 srcqual)
{
  CharPtr str = NULL;
  Int4    i;

  for (i = 0; i < NUM_srcqual_scqual && str == NULL; i++) {
    if (srcqual_scqual[i].srcqual == srcqual) {
      str = srcqual_scqual[i].qualname;
    }
  }
  if (str == NULL) {
    str = "Unknown source qualifier";
  }
  return str;
}


NLM_EXTERN Int4 GetSourceQualTypeByName (CharPtr qualname)
{
  Int4    i;

  qualname = GetCanonical (qualname, src_qual_alias_list);
  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (Matchnamestring(srcqual_scqual[i].qualname, qualname)) {
      return srcqual_scqual[i].srcqual;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetSourceQualList (Boolean for_remove)
{
  ValNodePtr list = NULL, tmp = NULL, last = NULL;
  Int4 i;

  if (for_remove) {
    ValNodeAddPointer (&list, 0, StringSave (kAllQualsStr));
    last = ValNodeAddPointer (&list, 0, StringSave (kAllNotesStr));
  }
  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual_scqual[i].srcqual != Source_qual_all_notes
        && srcqual_scqual[i].srcqual != Source_qual_all_quals) {
      ValNodeAddPointer (&tmp, 0, StringSave (srcqual_scqual[i].qualname));
    }
  }
  tmp = ValNodeSort (tmp, SortVnpByString);
  if (last == NULL) {
    list = tmp;
  } else {
    last->next = tmp;
  }
  return list;
}


NLM_EXTERN ValNodePtr GetSourceQualFieldListFromBioSource (BioSourcePtr biop)
{
  SubSourcePtr ssp;
  OrgModPtr    mod;
  ValNodePtr   list = NULL, vnp;
  Int4         i;

  if (biop == NULL) {
    return NULL;
  }

  vnp = ValNodeNew (NULL);
  vnp->choice = SourceQualChoice_textqual;
  vnp->data.intvalue = Source_qual_taxname;
  ValNodeAddPointer (&list, FieldType_source_qual, vnp);

  /* add other tax values */
  if (biop->org != NULL && !StringHasNoText (biop->org->common)) {
    vnp = ValNodeNew (NULL);
    vnp->choice = SourceQualChoice_textqual;
    vnp->data.intvalue = Source_qual_common_name;
    ValNodeAddPointer (&list, FieldType_source_qual, vnp);
  }
  if (biop->org != NULL && biop->org->orgname != NULL) {
    if (!StringHasNoText (biop->org->orgname->lineage)) {
      vnp = ValNodeNew (NULL);
      vnp->choice = SourceQualChoice_textqual;
      vnp->data.intvalue = Source_qual_lineage;
      ValNodeAddPointer (&list, FieldType_source_qual, vnp);
    }
    if (!StringHasNoText (biop->org->orgname->lineage)) {
      vnp = ValNodeNew (NULL);
      vnp->choice = SourceQualChoice_textqual;
      vnp->data.intvalue = Source_qual_division;
      ValNodeAddPointer (&list, FieldType_source_qual, vnp);
    }
  }

  /* add subtypes */
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    for (i = 0;
         i < NUM_srcqual_scqual && (srcqual_scqual[i].typeflag != IS_SUBSRC || srcqual_scqual[i].subtype != ssp->subtype);
         i++) {}
    if (i < NUM_srcqual_scqual) {
      vnp = ValNodeNew (NULL);
      vnp->choice = SourceQualChoice_textqual;    
      vnp->data.intvalue = srcqual_scqual[i].srcqual;
      ValNodeAddPointer (&list, FieldType_source_qual, vnp);
    }
  }
  /* add orgmods */
  if (biop->org != NULL && biop->org->orgname != NULL) {
    for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
      for (i = 0;
          i < NUM_srcqual_scqual && (srcqual_scqual[i].typeflag != IS_ORGMOD || srcqual_scqual[i].subtype != mod->subtype);
          i++) {}
      if (i < NUM_srcqual_scqual) {
        vnp = ValNodeNew (NULL);
        vnp->choice = SourceQualChoice_textqual;    
        vnp->data.intvalue = srcqual_scqual[i].srcqual;
        ValNodeAddPointer (&list, FieldType_source_qual, vnp);
      }
    }
  }
  return list;
}


NLM_EXTERN Boolean AllowSourceQualMulti (SourceQualChoicePtr s)
{
  Boolean rval = FALSE;

  if (s == NULL || s->choice != SourceQualChoice_textqual || s->data.ptrvalue == NULL) {
    return FALSE;
  } else if (s->data.intvalue == Source_qual_culture_collection
             || s->data.intvalue == Source_qual_bio_material
             || s->data.intvalue == Source_qual_specimen_voucher
             || s->data.intvalue == Source_qual_dbxref) {
    rval = TRUE;
  }
  return rval;
}


typedef struct srclocgenome {
  Int4 srcloc;
  Int4 genome;
  CharPtr name;
} SrcLocGenomeData, PNTR SrcLocGenomePtr;

static SrcLocGenomeData srcloc_genome[] = {
 { Source_location_unknown , GENOME_unknown , " " } ,
 { Source_location_genomic , GENOME_genomic , "genomic" } ,
 { Source_location_chloroplast , GENOME_chloroplast , "chloroplast" } ,
 { Source_location_chromoplast , GENOME_chromoplast , "chromoplast" } ,
 { Source_location_kinetoplast , GENOME_kinetoplast , "kinetoplast" } ,
 { Source_location_mitochondrion , GENOME_mitochondrion , "mitochondrion" } ,
 { Source_location_plastid , GENOME_plastid , "plastid" } ,
 { Source_location_macronuclear , GENOME_macronuclear , "macronuclear" } ,
 { Source_location_extrachrom , GENOME_extrachrom , "extrachromosomal" } ,
 { Source_location_plasmid , GENOME_plasmid , "plasmid" } ,
 { Source_location_transposon , GENOME_transposon , "transposon" } ,
 { Source_location_insertion_seq , GENOME_insertion_seq , "insertion-seq" } ,
 { Source_location_cyanelle , GENOME_cyanelle , "cyanelle" } ,
 { Source_location_proviral , GENOME_proviral , "proviral" } ,
 { Source_location_virion , GENOME_virion , "virion" } ,
 { Source_location_nucleomorph , GENOME_nucleomorph , "nucleomorph" } ,
 { Source_location_apicoplast , GENOME_apicoplast , "apicoplast" } ,
 { Source_location_leucoplast , GENOME_leucoplast , "leucoplast" } ,
 { Source_location_proplastid , GENOME_proplastid , "proplastid" } ,
 { Source_location_endogenous_virus , GENOME_endogenous_virus , "endogenous-virus" } ,
 { Source_location_hydrogenosome , GENOME_hydrogenosome , "hydrogenosome" } ,
 { Source_location_chromosome , GENOME_chromosome , "chromosome" } ,
 { Source_location_chromatophore , GENOME_chromatophore , "chromatophore" } };

#define NUM_srcloc_genome sizeof (srcloc_genome) / sizeof (SrcLocGenomeData)

NLM_EXTERN Int4 GenomeFromSrcLoc (Int4 srcloc) 
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (srcloc_genome[i].srcloc == srcloc) {
      return srcloc_genome[i].genome;
    }
  }
  return -1;
}


NLM_EXTERN Int4 SrcLocFromGenome (Int4 genome) 
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (srcloc_genome[i].genome == genome) {
      return srcloc_genome[i].srcloc;
    }
  }
  return -1;
}



NLM_EXTERN CharPtr LocNameFromGenome (Int4 genome) 
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (srcloc_genome[i].genome == genome) {
      return srcloc_genome[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN Int4 GenomeFromLocName (CharPtr loc_name)
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (StringICmp (srcloc_genome[i].name, loc_name) == 0) {
      return srcloc_genome[i].genome;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetLocationList (Boolean for_remove)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (for_remove && srcloc_genome[i].srcloc == Source_location_unknown) {
      ValNodeAddPointer (&list, srcloc_genome[i].srcloc, StringSave ("any"));
    } else {
      ValNodeAddPointer (&list, srcloc_genome[i].srcloc, StringSave (srcloc_genome[i].name));
    }
  }
  list = ValNodeSort (list, SortVnpByString);
  return list;
}


typedef struct srcorigorigin {
  Int4 srcorig;
  Int4 origin;
  CharPtr name;
} SrcOrigOriginData, PNTR SrcrigOriginPtr;

static SrcOrigOriginData srcorig_origin[] = {
 { Source_origin_unknown , 0 , "unknown" } ,
 { Source_origin_natural , 1 , "natural" } ,
 { Source_origin_natmut , 2 , "natmut" } ,
 { Source_origin_mut , 3 , "mut" } ,
 { Source_origin_artificial , 4 , "artificial" } ,
 { Source_origin_synthetic , 5 , "synthetic" } ,
 { Source_origin_other , 255 , "other" } };

#define NUM_srcorig_origin sizeof (srcorig_origin) / sizeof (SrcOrigOriginData)

NLM_EXTERN Int4 OriginFromSrcOrig (Int4 srcorig) 
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (srcorig_origin[i].srcorig == srcorig) {
      return srcorig_origin[i].origin;
    }
  }
  return -1;
}


NLM_EXTERN Int4 SrcOrigFromOrigin (Int4 origin) 
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (srcorig_origin[i].origin == origin) {
      return srcorig_origin[i].srcorig;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr OriginNameFromOrigin (Int4 origin) 
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (srcorig_origin[i].origin == origin) {
      return srcorig_origin[i].name;
    }
  }
  return NULL;
}


static Int4 OriginFromOriginName (CharPtr origin_name)
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (StringCmp (srcorig_origin[i].name, origin_name) == 0) {
      return srcorig_origin[i].origin;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetOriginList (Boolean for_remove)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (for_remove && srcorig_origin[i].srcorig == Source_origin_unknown) {
      ValNodeAddPointer (&list, srcorig_origin[i].srcorig, StringSave ("any"));
    } else {
      ValNodeAddPointer (&list, srcorig_origin[i].srcorig, StringSave (srcorig_origin[i].name));
    }
  }
  return list;
}


/* special code for converting source features to source qualifier val lists */
static void SetSrcQualTextValue (ValNodePtr PNTR fields, Int4 srcqual, CharPtr val)
{
  SourceQualTextValPtr st;

  st = SourceQualTextValNew ();
  st->srcqual = srcqual;
  st->val = StringSave (val);
  ValNodeAddPointer (fields, SourceQualValChoice_textqual, st);
}


static ValNodePtr SourceQualValsFromOrgMods (OrgModPtr mod)
{
  Int4 src_qual;
  ValNodePtr fields = NULL;

  while (mod != NULL) {
    src_qual = GetSrcQualFromSubSrcOrOrgMod (mod->subtype, TRUE);
    if (src_qual > -1) {
      SetSrcQualTextValue (&fields, src_qual, mod->subname);
    }
    mod = mod->next;
  }
  return fields;
}


static ValNodePtr SourceQualValsFromSubSrcs (SubSourcePtr ssp)
{
  Int4 src_qual;
  ValNodePtr fields = NULL;

  while (ssp != NULL) {
    src_qual = GetSrcQualFromSubSrcOrOrgMod (ssp->subtype, FALSE);
    if (src_qual > -1) {
      SetSrcQualTextValue (&fields, src_qual, ssp->name);
    }
    ssp = ssp->next;
  }
  return fields;
}


static ValNodePtr SourceQualValsFromSynonyms (ValNodePtr syn)
{
  ValNodePtr fields = NULL;

  while (syn != NULL) {
    SetSrcQualTextValue (&fields, Source_qual_synonym, syn->data.ptrvalue);
    syn = syn->next;
  }
  return fields;
}


static CharPtr GetDbtagString (DbtagPtr db_tag);

static ValNodePtr SourceQualValsFromDbxrefs (ValNodePtr dbxref)
{
  ValNodePtr fields = NULL;
  CharPtr tmp;

  while (dbxref != NULL) {
    tmp = GetDbtagString (dbxref->data.ptrvalue);
    SetSrcQualTextValue (&fields, Source_qual_dbxref, tmp);
    dbxref = dbxref->next;
  }
  return fields;
}


NLM_EXTERN ValNodePtr SourceQualValsFromBioSourcePtr (BioSourcePtr biop)
{
  ValNodePtr fields = NULL;
  Int4 loc, origin;

  if (biop == NULL) {
    return NULL;
  }

  ValNodeLink (&fields, SourceQualValsFromSubSrcs (biop->subtype));

  /* genome */
  if (biop->genome != GENOME_unknown) {
    loc = SrcLocFromGenome (biop->genome);
    if (loc > -1) {
      ValNodeAddInt (&fields, SourceQualValChoice_location, loc);
    }
  }
  /* origin */
  if (origin > 0) {
    origin = SrcOrigFromOrigin (biop->origin);
    if (origin > -1) {
      ValNodeAddInt (&fields, SourceQualValChoice_origin, origin);
    }
  }
  /* TODO: need focus */


  if (biop->org != NULL) {
    if (!StringHasNoText (biop->org->taxname)) {
      SetSrcQualTextValue (&fields, Source_qual_taxname, biop->org->taxname);
    }
    /* need common */
    if (!StringHasNoText (biop->org->common)) {
      SetSrcQualTextValue (&fields, Source_qual_common, biop->org->common);
    }
    /* dbxrefs */
    ValNodeLink (&fields, SourceQualValsFromDbxrefs (biop->org->db));

    /* add synonyms */
    SourceQualValsFromSynonyms (biop->org->syn);

    if (biop->org->orgname != NULL) {
      ValNodeLink (&fields, SourceQualValsFromOrgMods (biop->org->orgname->mod));
      
      /* lineage */
      if (!StringHasNoText (biop->org->orgname->lineage)) {
        SetSrcQualTextValue (&fields, Source_qual_lineage, biop->org->orgname->lineage);
      }
      /* div */
      if (!StringHasNoText (biop->org->orgname->div)) {
        SetSrcQualTextValue (&fields, Source_qual_division, biop->org->orgname->div);
      }

      /* gcode, mgcode */
      if (biop->org->orgname->gcode > 0) {
        ValNodeAddInt (&fields, SourceQualChoice_gcode, biop->org->orgname->gcode);
      }
      if (biop->org->orgname->mgcode > 0) {
        ValNodeAddInt (&fields, SourceQualChoice_mgcode, biop->org->orgname->mgcode);
      }

    }

  }

  return fields;
}


static void SetSourceQualValOnBioSource (BioSourcePtr biop, ValNodePtr src_qual)
{
  ValNode vn;
  SourceQualTextValPtr st;

  if (biop == NULL || src_qual == NULL) {
    return;
  }

  vn.next = NULL;
  switch (src_qual->choice) {
    case SourceQualValChoice_textqual:
      st = (SourceQualTextValPtr) src_qual->data.ptrvalue;
      if (st != NULL) {
        vn.choice = SourceQualChoice_textqual;
        vn.data.intvalue = st->srcqual;
        if (AllowSourceQualMulti (src_qual)) {
          SetSourceQualInBioSource (biop, &vn, NULL, st->val, ExistingTextOption_add_qual);
        } else {
          SetSourceQualInBioSource (biop, &vn, NULL, st->val, ExistingTextOption_replace_old);
        }
      }
      break;
    case SourceQualValChoice_location:
      vn.choice = SourceQualChoice_location;
      vn.data.intvalue = src_qual->data.intvalue;
      SetSourceQualInBioSource (biop, &vn, NULL, NULL, ExistingTextOption_replace_old);
      break;
    case SourceQualValChoice_origin:
      vn.choice = SourceQualChoice_origin;
      vn.data.intvalue = src_qual->data.intvalue;
      SetSourceQualInBioSource (biop, &vn, NULL, NULL, ExistingTextOption_replace_old);
      break;
    case SourceQualValChoice_gcode:
      vn.choice = SourceQualChoice_gcode;
      vn.data.intvalue = src_qual->data.intvalue;
      SetSourceQualInBioSource (biop, &vn, NULL, NULL, ExistingTextOption_replace_old);
      break;
    case SourceQualValChoice_mgcode:
      vn.choice = SourceQualChoice_mgcode;
      vn.data.intvalue = src_qual->data.intvalue;
      SetSourceQualInBioSource (biop, &vn, NULL, NULL, ExistingTextOption_replace_old);
      break;
  }
}


NLM_EXTERN BioSourcePtr BioSourceFromSourceQualVals (ValNodePtr fields)
{
  BioSourcePtr biop = NULL;
  ValNodePtr vnp;

  if (fields != NULL) {
    biop = BioSourceNew ();

    for (vnp = fields; vnp != NULL; vnp = vnp->next) {
      SetSourceQualValOnBioSource (biop, vnp);
    }
  }
  return biop;
}





typedef struct cdsgeneprotfieldname {
  Int4 field;
  CharPtr name;
} CDSGeneProtFieldNameData, PNTR CDSGeneProtFieldNamePtr;

static CDSGeneProtFieldNameData cdsgeneprotfield_name[] = {
{ CDSGeneProt_field_cds_comment , "CDS comment" } ,
{ CDSGeneProt_field_cds_inference , "CDS inference" } ,
{ CDSGeneProt_field_codon_start , "codon-start" } ,
{ CDSGeneProt_field_gene_locus , "gene locus" } ,
{ CDSGeneProt_field_gene_description , "gene description" } ,
{ CDSGeneProt_field_gene_comment , "gene comment" } ,
{ CDSGeneProt_field_gene_inference, "gene inference" } ,
{ CDSGeneProt_field_gene_allele , "gene allele" } ,
{ CDSGeneProt_field_gene_maploc , "gene maploc" } ,
{ CDSGeneProt_field_gene_locus_tag , "gene locus tag" } ,
{ CDSGeneProt_field_gene_synonym , "gene synonym" } ,
{ CDSGeneProt_field_gene_old_locus_tag , "gene old locus tag" } ,
{ CDSGeneProt_field_mrna_product , "mRNA product" } ,
{ CDSGeneProt_field_mrna_comment , "mRNA comment" } ,
{ CDSGeneProt_field_prot_name , "protein name" } ,
{ CDSGeneProt_field_prot_description , "protein description" } ,
{ CDSGeneProt_field_prot_ec_number , "protein EC number" } ,
{ CDSGeneProt_field_prot_activity , "protein activity" } ,
{ CDSGeneProt_field_prot_comment , "protein comment" } ,
{ CDSGeneProt_field_mat_peptide_name , "mat-peptide name" } ,
{ CDSGeneProt_field_mat_peptide_description ,  "mat-peptide description" } ,
{ CDSGeneProt_field_mat_peptide_ec_number , "mat-peptide EC number" } ,
{ CDSGeneProt_field_mat_peptide_activity , "mat-peptide activity" } ,
{ CDSGeneProt_field_mat_peptide_comment , "mat-peptide comment" } };

#define NUM_cdsgeneprotfield_name sizeof (cdsgeneprotfield_name) / sizeof (CDSGeneProtFieldNameData)

NLM_EXTERN CharPtr CDSGeneProtNameFromField (Int4 field) 
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    if (cdsgeneprotfield_name[i].field == field) {
      return cdsgeneprotfield_name[i].name;
    }
  }
  return NULL;
}


static Int4 CDSGeneProtFieldFromName (CharPtr str)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    if (Matchnamestring (cdsgeneprotfield_name[i].name, str)) {
      return cdsgeneprotfield_name[i].field;
    }
  }
  return -1;
}


NLM_EXTERN void AddAllCDSGeneProtFieldsToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    ValNodeAddPointer (field_list, cdsgeneprotfield_name[i].field, StringSave (cdsgeneprotfield_name[i].name));
  }
}


static ValNodePtr MakeCDSGeneProtFieldTypeList (void)
{
  Int4 i;
  ValNodePtr field_list = NULL;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    ValNodeAddInt (&field_list, FieldType_cds_gene_prot, cdsgeneprotfield_name[i].field);
  }
  return field_list;
}


typedef struct cdsgeneprotfeatname {
  Int4 feature_type;
  CharPtr name;
} CDSGeneProtFeatNameData, PNTR CDSGeneProtFeatNamePtr;

static CDSGeneProtFeatNameData cdsgeneprotfeat_name[] = {
{ CDSGeneProt_feature_type_constraint_gene , "gene" } ,
{ CDSGeneProt_feature_type_constraint_mRNA , "mRNA" } ,
{ CDSGeneProt_feature_type_constraint_cds , "CDS" } ,
{ CDSGeneProt_feature_type_constraint_prot , "protein" } ,
{ CDSGeneProt_feature_type_constraint_mat_peptide , "mat-peptide" }};

#define NUM_cdsgeneprotfeat_name sizeof (cdsgeneprotfeat_name) / sizeof (CDSGeneProtFeatNameData)

NLM_EXTERN CharPtr CDSGeneProtFeatureNameFromFeatureType (Int4 feature_type)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfeat_name; i++) {
    if (cdsgeneprotfeat_name[i].feature_type == feature_type) {
      return cdsgeneprotfeat_name[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN void AddAllCDSGeneProtFeaturesToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfeat_name; i++) {
    ValNodeAddPointer (field_list, cdsgeneprotfeat_name[i].feature_type, StringSave (cdsgeneprotfeat_name[i].name));
  }
}


static Boolean IsCDSGeneProtFieldMatPeptideRelated (Int4 val)
{
  if (val == CDSGeneProt_field_mat_peptide_name
      || val == CDSGeneProt_field_mat_peptide_description
      || val == CDSGeneProt_field_mat_peptide_ec_number
      || val == CDSGeneProt_field_mat_peptide_activity
      || val == CDSGeneProt_field_mat_peptide_comment) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsFieldTypeMatPeptideRelated (FieldTypePtr field)
{
  Boolean rval = FALSE;
  FeatureFieldPtr ff;

  if (field == NULL) {
    rval = FALSE;
  } else if ((field->choice == FieldType_feature_field
       && (ff = field->data.ptrvalue) != NULL
       && ff->type == Feature_type_mat_peptide_aa)
      || (field->choice == FieldType_cds_gene_prot
          && IsCDSGeneProtFieldMatPeptideRelated(field->data.intvalue))) {
    rval = TRUE;
  } else {
    rval = FALSE;
  }
  return rval;
}

  
static Boolean IsConstraintChoiceMatPeptideRelated (ConstraintChoicePtr constraint)
{
  CDSGeneProtQualConstraintPtr cq;
  FieldConstraintPtr fq;
  Boolean            rval = FALSE;

  if (constraint == NULL) {
    rval = FALSE;
  } else if (constraint->choice == ConstraintChoice_cdsgeneprot_qual) {
    cq = (CDSGeneProtQualConstraintPtr) constraint->data.ptrvalue;
    if (cq != NULL && cq->field1 != NULL 
        && IsCDSGeneProtFieldMatPeptideRelated (cq->field1->data.intvalue)) {
      rval = TRUE;
    } else {
      rval = FALSE;
    }
  } else if (constraint->choice == ConstraintChoice_field) {
    fq = (FieldConstraintPtr) constraint->data.ptrvalue;
    if (fq != NULL && IsFieldTypeMatPeptideRelated (fq->field)) {
      rval = TRUE;
    } else {
      rval = FALSE;
    }
  } else {
    rval = FALSE;
  }
  return rval;
}


static Int2 FeatureTypeFromCDSGeneProtField (Uint2 cds_gene_prot_field)
{
  Int2 feat_type = Feature_type_any;

  switch (cds_gene_prot_field) {
    case CDSGeneProt_field_cds_comment:
    case CDSGeneProt_field_cds_inference:
    case CDSGeneProt_field_codon_start:
      feat_type = Feature_type_cds;
      break;
    case CDSGeneProt_field_gene_locus:
    case CDSGeneProt_field_gene_description:
    case CDSGeneProt_field_gene_comment:
    case CDSGeneProt_field_gene_allele:
    case CDSGeneProt_field_gene_maploc:
    case CDSGeneProt_field_gene_locus_tag:
    case CDSGeneProt_field_gene_synonym:
    case CDSGeneProt_field_gene_old_locus_tag:
    case CDSGeneProt_field_gene_inference:
      feat_type = Feature_type_gene;
      break;
    case CDSGeneProt_field_mrna_product:
    case CDSGeneProt_field_mrna_comment:
      feat_type = Feature_type_mRNA;
      break;
    case CDSGeneProt_field_prot_name:
    case CDSGeneProt_field_prot_description:
    case CDSGeneProt_field_prot_ec_number:
    case CDSGeneProt_field_prot_activity:
    case CDSGeneProt_field_prot_comment:
      feat_type = Feature_type_prot;
      break;
    case CDSGeneProt_field_mat_peptide_name:
    case CDSGeneProt_field_mat_peptide_description:
    case CDSGeneProt_field_mat_peptide_ec_number:
    case CDSGeneProt_field_mat_peptide_activity:
    case CDSGeneProt_field_mat_peptide_comment:
      feat_type = Feature_type_mat_peptide_aa;
      break;
  }
  return feat_type;
}


NLM_EXTERN FeatureFieldPtr FeatureFieldFromCDSGeneProtField (Uint2 cds_gene_prot_field)
{
  FeatureFieldPtr f = NULL;

  switch (cds_gene_prot_field) {
    case CDSGeneProt_field_cds_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_cds;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_cds_inference:
      f = FeatureFieldNew ();
      f->type = Feature_type_cds;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_inference;
      break;
    case CDSGeneProt_field_codon_start:
      f = FeatureFieldNew ();
      f->type = Feature_type_cds;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_codon_start;
      break;
    case CDSGeneProt_field_gene_locus:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_gene;
      break;
    case CDSGeneProt_field_gene_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_gene_description;
      break;
    case CDSGeneProt_field_gene_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_gene_allele:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_allele;
      break;
    case CDSGeneProt_field_gene_maploc:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_map;
      break;
    case CDSGeneProt_field_gene_locus_tag:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_locus_tag;
      break;
    case CDSGeneProt_field_gene_synonym:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_synonym;
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_old_locus_tag;
      break;
    case CDSGeneProt_field_gene_inference:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_inference;
      break;
    case CDSGeneProt_field_mrna_product:
      f = FeatureFieldNew ();
      f->type = Feature_type_mRNA;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_mrna_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_mRNA;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_prot_name:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_prot_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_description;
      break;
    case CDSGeneProt_field_prot_ec_number:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_ec_number;
      break;
    case CDSGeneProt_field_prot_activity:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_activity;
      break;
    case CDSGeneProt_field_prot_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_mat_peptide_name:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide_aa;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_mat_peptide_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide_aa;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_description;
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide_aa;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_ec_number;
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide_aa;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_activity;
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide_aa;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
  }
  return f;
}


static Uint2 CDSGeneProtFieldFromFeatureField (FeatureFieldPtr ffp)
{
  Uint2 cds_gene_prot_field = 0;

  if (ffp != NULL && ffp->field != NULL && ffp->field->choice == FeatQualChoice_legal_qual) {
    switch (ffp->field->data.intvalue) {
      case Feat_qual_legal_note:
        switch (ffp->type) {
          case Feature_type_cds:
            cds_gene_prot_field = CDSGeneProt_field_cds_comment;
            break;
          case Feature_type_gene:
            cds_gene_prot_field = CDSGeneProt_field_gene_comment;
            break;
          case Feature_type_mRNA:
            cds_gene_prot_field = CDSGeneProt_field_mrna_comment;
            break;
          case Feature_type_prot:
            cds_gene_prot_field = CDSGeneProt_field_prot_comment;
            break;
          case Feature_type_mat_peptide_aa:
            cds_gene_prot_field = CDSGeneProt_field_mat_peptide_comment;
            break;
        }
        break;
      case Feat_qual_legal_inference:
        switch (ffp->type) {
          case Feature_type_cds:
            cds_gene_prot_field = CDSGeneProt_field_cds_inference;
            break;
          case Feature_type_gene:
            cds_gene_prot_field = CDSGeneProt_field_gene_inference;
            break;
        }
        break;
      case Feat_qual_legal_codon_start:
        cds_gene_prot_field = CDSGeneProt_field_codon_start;
        break;
      case Feat_qual_legal_gene:
        cds_gene_prot_field = CDSGeneProt_field_gene_locus;
        break;
      case Feat_qual_legal_gene_description:
        cds_gene_prot_field = CDSGeneProt_field_gene_description;
        break;
      case Feat_qual_legal_allele:
        cds_gene_prot_field = CDSGeneProt_field_gene_allele;
        break;
      case Feat_qual_legal_map:
        cds_gene_prot_field = CDSGeneProt_field_gene_maploc;
        break;
      case Feat_qual_legal_locus_tag:
        cds_gene_prot_field = CDSGeneProt_field_gene_locus_tag;
        break;
      case Feat_qual_legal_synonym:
        cds_gene_prot_field = CDSGeneProt_field_gene_synonym;
        break;
      case Feat_qual_legal_old_locus_tag:
        cds_gene_prot_field = CDSGeneProt_field_gene_old_locus_tag;
        break;
      case Feat_qual_legal_product:
        switch (ffp->type) {
          case Feature_type_mRNA:
            cds_gene_prot_field = CDSGeneProt_field_mrna_product;
            break;
          case Feature_type_prot:
            cds_gene_prot_field = CDSGeneProt_field_prot_name;
            break;
          case Feature_type_mat_peptide_aa:
            cds_gene_prot_field = CDSGeneProt_field_mat_peptide_name;
            break;
        }
        break;
      case Feat_qual_legal_description:
        switch (ffp->type) {
          case Feature_type_gene:
            cds_gene_prot_field = CDSGeneProt_field_gene_description;
            break;
          case Feature_type_prot:
            cds_gene_prot_field = CDSGeneProt_field_prot_description;
            break;
          case Feature_type_mat_peptide_aa:
            cds_gene_prot_field = CDSGeneProt_field_mat_peptide_description;
            break;
        }
        break;
      case Feat_qual_legal_ec_number:
        switch (ffp->type) {
          case Feature_type_prot:
            cds_gene_prot_field = CDSGeneProt_field_prot_ec_number;
            break;
          case Feature_type_mat_peptide_aa:
            cds_gene_prot_field = CDSGeneProt_field_mat_peptide_ec_number;
            break;
        }
        break;
      case Feat_qual_legal_activity:
        switch (ffp->type) {
          case Feature_type_prot:
            cds_gene_prot_field = CDSGeneProt_field_prot_activity;
            break;
          case Feature_type_mat_peptide_aa:
            cds_gene_prot_field = CDSGeneProt_field_mat_peptide_activity;
            break;
        }
        break;
    }
  }
  return cds_gene_prot_field;
}


/* Molinfo fields */
typedef struct moleculetypebiomol {
  Int4 molecule_type;
  Int4 biomol;
  CharPtr name;
} MoleculeTypeBiomolData, PNTR MoleculeTypeBiomolPtr;

static MoleculeTypeBiomolData moleculetype_biomol[] = {
 { Molecule_type_unknown , 0, " " } ,
 { Molecule_type_genomic , MOLECULE_TYPE_GENOMIC , "genomic" } ,
 { Molecule_type_precursor_RNA , MOLECULE_TYPE_PRE_MRNA , "precursor RNA" } ,
 { Molecule_type_mRNA , MOLECULE_TYPE_MRNA , "mRNA" } ,
 { Molecule_type_rRNA , MOLECULE_TYPE_RRNA , "rRNA" } ,
 { Molecule_type_tRNA , MOLECULE_TYPE_TRNA , "tRNA" } ,
 { Molecule_type_genomic_mRNA , MOLECULE_TYPE_GENOMIC_MRNA_MIX , "genomic mRNA" } ,
 { Molecule_type_cRNA , MOLECULE_TYPE_CRNA , "cRNA" } ,
 { Molecule_type_transcribed_RNA, MOLECULE_TYPE_TRANSCRIBED_RNA, "transcribed RNA" } ,
 { Molecule_type_ncRNA, MOLECULE_TYPE_NCRNA, "ncRNA" } ,
 { Molecule_type_transfer_messenger_RNA, MOLECULE_TYPE_TMRNA, "tmRNA" } ,
 { Molecule_type_other, MOLECULE_TYPE_OTHER_GENETIC_MATERIAL, "other-genetic" }
};


#define NUM_moleculetype_biomol sizeof (moleculetype_biomol) / sizeof (MoleculeTypeBiomolData)

NLM_EXTERN Int4 BiomolFromMoleculeType (Int4 molecule_type) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (moleculetype_biomol[i].molecule_type == molecule_type) {
      return moleculetype_biomol[i].biomol;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr BiomolNameFromBiomol (Int4 biomol) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (moleculetype_biomol[i].biomol == biomol) {
      return moleculetype_biomol[i].name;
    }
  }
  return NULL;
}


static Int4 BiomolFromBiomolName (CharPtr biomol_name)
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (StringCmp (moleculetype_biomol[i].name, biomol_name) == 0) {
      return moleculetype_biomol[i].biomol;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetMoleculeTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    ValNodeAddPointer (&list, moleculetype_biomol[i].molecule_type, StringSave (moleculetype_biomol[i].name));
  }
  return list;
}


/* Technique fields */
typedef struct techniquetypetech {
  Int4 technique_type;
  Int4 tech;
  CharPtr name;
} TechniqueTypeTechData, PNTR TechniqueTypeTechPtr;

static TechniqueTypeTechData techniquetype_tech[] = {
 { Technique_type_unknown , MI_TECH_unknown , " " } ,
 { Technique_type_standard , MI_TECH_standard , "standard" } ,
 { Technique_type_est , MI_TECH_est , "EST" } ,
 { Technique_type_sts , MI_TECH_sts , "STS" } ,
 { Technique_type_survey , MI_TECH_survey , "survey" } ,
 { Technique_type_genetic_map , MI_TECH_genemap , "genetic map" } ,
 { Technique_type_physical_map , MI_TECH_physmap , "physical map" } ,
 { Technique_type_derived , MI_TECH_derived , "derived" } ,
 { Technique_type_concept_trans , MI_TECH_concept_trans , "concept-trans" } ,
 { Technique_type_seq_pept , MI_TECH_seq_pept , "seq-pept" } ,
 { Technique_type_both , MI_TECH_both , "both" } ,
 { Technique_type_seq_pept_overlap , MI_TECH_seq_pept_overlap , "seq-pept-overlap" } ,
 { Technique_type_seq_pept_homol , MI_TECH_seq_pept_homol, "seq-pept-homol" } ,
 { Technique_type_concept_trans_a, MI_TECH_concept_trans_a, "concept-trans-a" } ,
 { Technique_type_htgs_1, MI_TECH_htgs_1, "HTGS-1" } ,
 { Technique_type_htgs_2, MI_TECH_htgs_2, "HTGS-2" } ,
 { Technique_type_htgs_3, MI_TECH_htgs_3, "HTGS-3" } ,
 { Technique_type_fli_cDNA, MI_TECH_fli_cdna, "fli-cDNA" } ,
 { Technique_type_htgs_0, MI_TECH_htgs_0, "HTGS-0" } ,
 { Technique_type_htc, MI_TECH_htc, "HTC" } ,
 { Technique_type_wgs, MI_TECH_wgs, "WGS" } ,
 { Technique_type_barcode, MI_TECH_barcode, "BARCODE" } ,
 { Technique_type_composite_wgs_htgs, MI_TECH_composite_wgs_htgs, "composite WGS-HTGS" } ,
 { Technique_type_tsa, MI_TECH_tsa, "TSA" } ,
 { Technique_type_other, MI_TECH_other, "other" } 
};


#define NUM_techniquetype_tech sizeof (techniquetype_tech) / sizeof (TechniqueTypeTechData)

NLM_EXTERN Int4 TechFromTechniqueType (Int4 technique_type) 
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (techniquetype_tech[i].technique_type == technique_type) {
      return techniquetype_tech[i].tech;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr TechNameFromTech (Int4 tech) 
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (techniquetype_tech[i].tech == tech) {
      return techniquetype_tech[i].name;
    }
  }
  return NULL;
}


static Int4 TechFromTechName (CharPtr tech_name)
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (StringCmp (techniquetype_tech[i].name, tech_name) == 0) {
      return techniquetype_tech[i].tech;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetTechniqueTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    ValNodeAddPointer (&list, techniquetype_tech[i].technique_type, StringSave (techniquetype_tech[i].name));
  }
  return list;
}


/* Completedness fields */
typedef struct completednesstypecompleteness {
  Int4 completedness_type;
  Int4 completeness;
  CharPtr name;
} CompletednessTypeCompletenessData, PNTR CompletednessTypeCompletenessPtr;

static CompletednessTypeCompletenessData completednesstype_completeness[] = {
 { Completedness_type_unknown, 0, " " } ,
 { Completedness_type_complete, 1, "complete" } ,
 { Completedness_type_partial, 2, "partial" } ,
 { Completedness_type_no_left, 3, "no left" } ,
 { Completedness_type_no_right, 4, "no right" } ,
 { Completedness_type_no_ends, 5, "no ends" } ,
 { Completedness_type_has_left, 6, "has left" } ,
 { Completedness_type_has_right, 7, "has right" } ,
 { Completedness_type_other, 255, "other" }
};

#define NUM_completednesstype_completeness sizeof (completednesstype_completeness) / sizeof (CompletednessTypeCompletenessData)

NLM_EXTERN Int4 CompletenessFromCompletednessType (Int4 completedness_type) 
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (completednesstype_completeness[i].completedness_type == completedness_type) {
      return completednesstype_completeness[i].completeness;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr CompletenessNameFromCompleteness (Int4 completeness) 
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (completednesstype_completeness[i].completeness == completeness) {
      return completednesstype_completeness[i].name;
    }
  }
  return NULL;
}


static Int4 CompletenessFromCompletenessName (CharPtr completeness_name)
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (StringCmp (completednesstype_completeness[i].name, completeness_name) == 0) {
      return completednesstype_completeness[i].completeness;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetCompletednessTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    ValNodeAddPointer (&list, completednesstype_completeness[i].completedness_type, StringSave (completednesstype_completeness[i].name));
  }
  return list;
}


/* Molecule class fields */
typedef struct moleculeclasstypemol {
  Int4 moleculeclass_type;
  Int4 mol;
  CharPtr name;
} MoleculeClassTypeMolData, PNTR MoleculeClassTypeMolPtr;

static MoleculeClassTypeMolData moleculeclasstype_mol[] = {
 { Molecule_class_type_unknown, 0, " " } ,
 { Molecule_class_type_dna, MOLECULE_CLASS_DNA, "DNA" } ,
 { Molecule_class_type_rna, MOLECULE_CLASS_RNA, "RNA" } ,
 { Molecule_class_type_protein, MOLECULE_CLASS_PROTEIN, "protein" } ,
 { Molecule_class_type_nucleotide, MOLECULE_CLASS_NUC, "nucleotide" } ,
 { Molecule_class_type_other, 255, "other" } 
};


#define NUM_moleculeclasstype_mol sizeof (moleculeclasstype_mol) / sizeof (MoleculeClassTypeMolData)

NLM_EXTERN Int4 MolFromMoleculeClassType (Int4 moleculeclass_type) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (moleculeclasstype_mol[i].moleculeclass_type == moleculeclass_type) {
      return moleculeclasstype_mol[i].mol;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr MolNameFromMol (Int4 mol) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (moleculeclasstype_mol[i].mol == mol) {
      return moleculeclasstype_mol[i].name;
    }
  }
  return NULL;
}


static Int4 MolFromMolName (CharPtr mol_name)
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (StringCmp (moleculeclasstype_mol[i].name, mol_name) == 0) {
      return moleculeclasstype_mol[i].mol;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetMoleculeClassTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    ValNodeAddPointer (&list, moleculeclasstype_mol[i].moleculeclass_type, StringSave (moleculeclasstype_mol[i].name));
  }
  return list;
}


/* Topology fields */
typedef struct topologytypetopology {
  Int4 topology_type;
  Int4 topology;
  CharPtr name;
} TopologyTypeTopologyData, PNTR TopologyTypeTopologyPtr;

static TopologyTypeTopologyData topologytype_topology[] = {
 { Topology_type_unknown, 0, " " } ,
 { Topology_type_linear, TOPOLOGY_LINEAR, "linear" } ,
 { Topology_type_circular, TOPOLOGY_CIRCULAR, "circular" } ,
 { Topology_type_tandem, TOPOLOGY_TANDEM, "tandem" } ,
 { Topology_type_other, 255, "other" } 
};

#define NUM_topologytype_topology sizeof (topologytype_topology) / sizeof (TopologyTypeTopologyData)

NLM_EXTERN Int4 TopologyFromTopologyType (Int4 topology_type) 
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (topologytype_topology[i].topology_type == topology_type) {
      return topologytype_topology[i].topology;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr TopologyNameFromTopology (Int4 topology) 
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (topologytype_topology[i].topology == topology) {
      return topologytype_topology[i].name;
    }
  }
  return NULL;
}


static Int4 TopologyFromTopologyName (CharPtr topology_name)
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (StringCmp (topologytype_topology[i].name, topology_name) == 0) {
      return topologytype_topology[i].topology;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetTopologyTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    ValNodeAddPointer (&list, topologytype_topology[i].topology_type, StringSave (topologytype_topology[i].name));
  }
  return list;
}


/* strand fields */
typedef struct strandtypestrand {
  Int4 strand_type;
  Int4 strand;
  CharPtr name;
} StrandTypeStrandData, PNTR StrandTypeStrandPtr;

static StrandTypeStrandData strandtype_strand[] = {
 { Strand_type_unknown, 0, " " } ,
 { Strand_type_single, STRANDEDNESS_SINGLE, "single" } ,
 { Strand_type_double__, STRANDEDNESS_DOUBLE, "double" } ,
 { Strand_type_mixed, 3, "mixed" } ,
 { Strand_type_mixed_rev, 4, "mixed-rev" } ,
 { Strand_type_other, 255, "other" } 
};

#define NUM_strandtype_strand sizeof (strandtype_strand) / sizeof (StrandTypeStrandData)

NLM_EXTERN Int4 StrandFromStrandType (Int4 strand_type) 
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (strandtype_strand[i].strand_type == strand_type) {
      return strandtype_strand[i].strand;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr StrandNameFromStrand (Int4 strand) 
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (strandtype_strand[i].strand == strand) {
      return strandtype_strand[i].name;
    }
  }
  return NULL;
}


static Int4 StrandFromStrandName (CharPtr strand_name)
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (StringCmp (strandtype_strand[i].name, strand_name) == 0) {
      return strandtype_strand[i].strand;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetStrandTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    ValNodeAddPointer (&list, strandtype_strand[i].strand_type, StringSave (strandtype_strand[i].name));
  }
  return list;
}


static CharPtr GetSequenceQualValName (ValNodePtr field)
{
  CharPtr val = NULL;

  if (field == NULL) return NULL;
  switch (field->choice) {
    case MolinfoField_molecule:
      val = BiomolNameFromBiomol (BiomolFromMoleculeType (field->data.intvalue));
      break;
    case MolinfoField_technique:
      val = TechNameFromTech (TechFromTechniqueType (field->data.intvalue));
      break;
    case MolinfoField_completedness:
      val = CompletenessNameFromCompleteness (CompletenessFromCompletednessType (field->data.intvalue));
      break;
    case MolinfoField_mol_class:
      val = MolNameFromMol (MolFromMoleculeClassType (field->data.intvalue));
      break;
    case MolinfoField_topology:
      val = TopologyNameFromTopology (TopologyFromTopologyType (field->data.intvalue));
      break;
    case MolinfoField_strand:
      val = StrandNameFromStrand (StrandFromStrandType (field->data.intvalue));
      break;
  }
  return val;
}


static CharPtr GetSequenceQualName (ValNodePtr field)
{
  CharPtr str = NULL, fieldname = "invalid field", val = "invalid value";
  CharPtr fmt = "%s %s";

  if (field == NULL) return NULL;
  switch (field->choice) {
    case MolinfoField_molecule:
      fieldname = "molecule";
      val = BiomolNameFromBiomol (BiomolFromMoleculeType (field->data.intvalue));
      break;
    case MolinfoField_technique:
      fieldname = "technique";
      val = TechNameFromTech (TechFromTechniqueType (field->data.intvalue));
      break;
    case MolinfoField_completedness:
      fieldname = "completeness";
      val = CompletenessNameFromCompleteness (CompletenessFromCompletednessType (field->data.intvalue));
      break;
    case MolinfoField_mol_class:
      fieldname = "class";
      val = MolNameFromMol (MolFromMoleculeClassType (field->data.intvalue));
      break;
    case MolinfoField_topology:
      fieldname = "topology";
      val = TopologyNameFromTopology (TopologyFromTopologyType (field->data.intvalue));
      break;
    case MolinfoField_strand:
      fieldname = "strand";
      val = StrandNameFromStrand (StrandFromStrandType (field->data.intvalue));
      break;
  }
  if (val == NULL) {
    val = "Invalid value";
  }
  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (fieldname) + StringLen (val)));
  sprintf (str, fmt, fieldname, val);
  return str;
}


static ValNodePtr MakeSequenceQualFieldTypeList (void)
{
  ValNodePtr field_list = NULL;
  ValNodePtr field;

  field = ValNodeNew (NULL);
  field->choice = MolinfoField_molecule;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  field = ValNodeNew (NULL);
  field->choice = MolinfoField_technique;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  field = ValNodeNew (NULL);
  field->choice = MolinfoField_completedness;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  field = ValNodeNew (NULL);
  field->choice = MolinfoField_mol_class;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  field = ValNodeNew (NULL);
  field->choice = MolinfoField_topology;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  field = ValNodeNew (NULL);
  field->choice = MolinfoField_strand;
  field->data.ptrvalue = NULL;
  ValNodeAddPointer (&field_list, FieldType_molinfo_field, field);
  return field_list;
}


/* bond types */
typedef struct bondtype {
  Int4 macro_bond_type;
  Int4 asn1_bond_type;
  CharPtr name;
} BondTypeData, PNTR BondTypePtr;

static BondTypeData bond_type[] = {
 { Bond_type_disulfide, 1, "Disulfide" } ,
 { Bond_type_thioester, 2, "Thioester" } ,
 { Bond_type_crosslink, 3, "Crosslink" } ,
 { Bond_type_thioether, 4, "Thioether" } ,
 { Bond_type_other, 255, "Other" } 
};

#define NUM_bond_type sizeof (bond_type) / sizeof (BondTypeData)

NLM_EXTERN Int4 Asn1BondTypeFromMacroBondType (Int4 macro_bond_type) 
{
  Int4 i;

  for (i = 0; i < NUM_bond_type; i++) {
    if (bond_type[i].macro_bond_type == macro_bond_type) {
      return bond_type[i].asn1_bond_type;
    }
  }
  return -1;
}


NLM_EXTERN Int4 MacroBondTypeFromAsn1BondType (Int4 asn1_bond_type) 
{
  Int4 i;

  for (i = 0; i < NUM_bond_type; i++) {
    if (bond_type[i].asn1_bond_type == asn1_bond_type) {
      return bond_type[i].macro_bond_type;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr GetMacroBondTypeName (Int4 macro_bond_type)
{
  Int4 i;

  for (i = 0; i < NUM_bond_type; i++) {
    if (bond_type[i].macro_bond_type == macro_bond_type) {
      return bond_type[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN ValNodePtr GetBondTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_bond_type; i++) {
    ValNodeAddPointer (&list, bond_type[i].macro_bond_type, StringSave (bond_type[i].name));
  }
  return list;
}


/* site types */
typedef struct sitetype {
  Int4 macro_site_type;
  Int4 asn1_site_type;
  CharPtr name;
} SiteTypeData, PNTR SiteTypePtr;

static SiteTypeData site_type[] = {
  {Site_type_active, 1, "Active"},
  {Site_type_binding, 2, "Binding"},
  {Site_type_cleavage, 3, "Cleavage"},
  {Site_type_inhibit, 4, "Inhibit"},
  {Site_type_modified, 5, "Modified"},
  {Site_type_glycosylation, 6, "Glycosylation"},
  {Site_type_myristoylation, 7, "Myristoylation"},
  {Site_type_mutagenized, 8, "Mutagenized"},
  {Site_type_metal_binding, 9, "Metal-binding"},
  {Site_type_phosphorylation, 10, "Phosphorylation"},
  {Site_type_acetylation, 11, "Acetylation"},
  {Site_type_amidation, 12, "Amidation"},
  {Site_type_methylation, 13, "Methylation"},
  {Site_type_hydroxylation, 14, "Hydroxylation"},
  {Site_type_sulfatation, 15, "Sulfatation"},
  {Site_type_oxidative_deamination, 16, "Oxidative-deamination"},
  {Site_type_pyrrolidone_carboxylic_acid, 17, "Pyrrolidone-carboxylic-acid"},
  {Site_type_gamma_carboxyglutamic_acid, 18, "Gamma-carboxyglutamic-acid"},
  {Site_type_blocked, 19, "Blocked"},
  {Site_type_lipid_binding, 20, "Lipid-binding"},
  {Site_type_np_binding, 21, "np-binding"},
  {Site_type_dna_binding, 22, "DNA-binding"},
  {Site_type_signal_peptide, 23, "Signal-peptide"},
  {Site_type_transit_peptide, 24, "Transit-peptide"},
  {Site_type_transmembrane_region, 25, "Transmembrane-region"},
  {Site_type_nitrosylation, 26, "Nitrosylation"},
  {Site_type_other, 255, "Other"},
};


#define NUM_site_type sizeof (site_type) / sizeof (SiteTypeData)

NLM_EXTERN Int4 Asn1SiteTypeFromMacroSiteType (Int4 macro_site_type) 
{
  Int4 i;

  for (i = 0; i < NUM_site_type; i++) {
    if (site_type[i].macro_site_type == macro_site_type) {
      return site_type[i].asn1_site_type;
    }
  }
  return -1;
}


NLM_EXTERN Int4 MacroSiteTypeFromAsn1SiteType (Int4 asn1_site_type) 
{
  Int4 i;

  for (i = 0; i < NUM_site_type; i++) {
    if (site_type[i].asn1_site_type == asn1_site_type) {
      return site_type[i].macro_site_type;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr GetMacroSiteTypeName (Int4 macro_site_type)
{
  Int4 i;

  for (i = 0; i < NUM_site_type; i++) {
    if (site_type[i].macro_site_type == macro_site_type) {
      return site_type[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN ValNodePtr GetSiteTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_site_type; i++) {
    ValNodeAddPointer (&list, site_type[i].macro_site_type, StringSave (site_type[i].name));
  }
  return list;
}


/* Simple constraints */
static Boolean IsWholeWordMatch (CharPtr start, CharPtr found, Int4 match_len)
{
  Boolean rval = TRUE;
  Char    char_after;
  Char    char_before;
  
  if (match_len == 0)
  {
    rval = TRUE;
  }
  else if (start == NULL || found == NULL)
  {
    rval = FALSE;
  }
  else
  {
	  char_after = *(found + match_len);
    if (found != start)
	  {
	    char_before = *(found - 1);
	    if (isalpha ((Int4) char_before) || isdigit ((Int4) char_before))
	    {
	      rval = FALSE;
	    }
	  }
	  if (char_after != 0 && (isalpha ((Int4) char_after) || isdigit ((Int4)char_after)))
	  {
	    rval = FALSE;
	  }   
  }
  return rval;
}


NLM_EXTERN Boolean IsStringConstraintEmpty (StringConstraintPtr scp)
{
  if (scp == NULL || StringHasNoText (scp->match_text)) return TRUE;
  else return FALSE;
}


NLM_EXTERN Boolean DoesSingleStringMatchConstraint (CharPtr str, StringConstraintPtr scp)
{
  CharPtr pFound;
  Boolean rval = FALSE;
  Char    char_after = 0;
  
  if (IsStringConstraintEmpty (scp)) return TRUE;
  if (StringHasNoText (str)) return FALSE;

  switch (scp->match_location) 
  {
    case String_location_contains:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else if (scp->whole_word) 
      {
        rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
        while (!rval && pFound != NULL) 
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, scp->match_text);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, scp->match_text);
	        }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
          }
        }
      }
      else
      {
        rval = TRUE;
      }
      break;
    case String_location_starts:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      if (pFound == str)
      {
        if (scp->whole_word) 
        {
          rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
        }
        else
        {
          rval = TRUE;
        }
      }
      break;
    case String_location_ends:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      while (pFound != NULL && !rval) {
  	    char_after = *(pFound + StringLen (scp->match_text));
        if (char_after == 0)
        {
          if (scp->whole_word) 
          {
            rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
          }
          else
          {
            rval = TRUE;
          }
          /* stop the search, we're at the end of the string */
          pFound = NULL;
        }
        else
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, scp->match_text);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, scp->match_text);
	        }
        }
      }
      break;
    case String_location_equals:
      if (scp->case_sensitive) 
      {
        if (StringCmp (str, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      else
      {
        if (StringICmp (str, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      break;
    case String_location_inlist:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (scp->match_text, str);
	    }
	    else
	    {
	      pFound = StringISearch (scp->match_text, str);
	    }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else
      {
        rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (str));
        while (!rval && pFound != NULL) 
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, str);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, str);
	        }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (str));
          }
        }
      }
      if (!rval) {
        /* look for spans */
        rval = IsStringInSpanInList (str, scp->match_text);
      }
      break;
	}
	return rval;
}


NLM_EXTERN Boolean DoesStringMatchConstraint (CharPtr str, StringConstraintPtr scp)
{
  Boolean rval;

  rval = DoesSingleStringMatchConstraint (str, scp);
  if (scp != NULL && scp->not_present) {
    rval = !rval;
  }
  return rval;
}


static Boolean DoesStringListMatchConstraint (ValNodePtr list, StringConstraintPtr scp)
{
  Int4 len = 1;
  CharPtr tmp;
  Boolean rval = FALSE;
  ValNodePtr vnp;

  if (IsStringConstraintEmpty (scp)) {
    return TRUE;
  }
  if (list == NULL) return FALSE;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    len += StringLen (vnp->data.ptrvalue) + 2;
  }

  tmp = (CharPtr) MemNew (sizeof (Char) * len);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    StringCat (tmp, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (tmp, "; ");
    }
  }

  rval = DoesStringMatchConstraint (tmp, scp);
  tmp = MemFree (tmp);
  return rval;  
}



NLM_EXTERN Boolean RemoveStringConstraintPortionFromString (CharPtr PNTR str, StringConstraintPtr scp)
{
  CharPtr pFound, src, dst, cp;
  Boolean rval = FALSE;
  Int4    match_len;

  if (str == NULL || *str == NULL) return FALSE;

  if (IsStringConstraintEmpty (scp) || scp->not_present) return FALSE;

  if (scp->match_location == String_location_equals) {
    if (scp->case_sensitive) {
      if (StringCmp (*str, scp->match_text) == 0) {
        rval = TRUE;
      }
    } else {
      if (StringICmp (*str, scp->match_text) == 0) {
        rval = TRUE;
      }
    }
    if (rval == TRUE) {
      **str = 0;
    }
  } else {
    match_len = StringLen (scp->match_text);
    if (scp->case_sensitive) {
      pFound = StringSearch (*str, scp->match_text);
    } else {
      pFound = StringISearch (*str, scp->match_text);
    }
    while (pFound != NULL) {
      switch (scp->match_location) {
        case String_location_contains:
        case String_location_inlist:
          if ((!scp->whole_word && scp->match_location != String_location_inlist)
              || IsWholeWordMatch (*str, pFound, match_len)) {
            src = pFound + match_len;
            dst = pFound;
            while (*src != 0) {
              *dst = *src;
              dst++; 
              src++;
            }
            *dst = 0;
            rval = TRUE;
            cp = pFound;
          } else {
            cp = pFound + 1;
          }
          if (scp->case_sensitive) {
            pFound = StringSearch (cp, scp->match_text);
          } else {
            pFound = StringISearch (cp, scp->match_text);
          }
          break;

        case String_location_starts:
          if (pFound == *str && (!scp->whole_word || IsWholeWordMatch (*str, pFound, match_len))) {
            src = pFound + match_len;
            dst = pFound;
            while (*src != 0) {
              *dst = *src;
              dst++; 
              src++;
            }
            *dst = 0;
            rval = TRUE;
          }
          pFound = NULL;
          break;
        case String_location_ends:
          if (*(pFound + match_len) == 0 && (!scp->whole_word || IsWholeWordMatch (*str, pFound, match_len))) {
            *pFound = 0;
            rval = TRUE;
            pFound = NULL;
          } else {     
	          if (scp->case_sensitive)
	          {
	            pFound = StringSearch (pFound + 1, scp->match_text);
	          }
	          else
	          {
	            pFound = StringISearch (pFound + 1, scp->match_text);
	          }
          }
          break;
      }
    }
	}
  if (rval && StringHasNoText (*str)) {
    *str = MemFree (*str);
  }
	return rval;  
}


NLM_EXTERN Boolean IsLocationConstraintEmpty (LocationConstraintPtr lcp)
{
  Boolean rval = TRUE;

  if (lcp == NULL)
  {
    rval = TRUE;
  }
  else if (lcp->strand != Strand_constraint_any)
  {
    rval = FALSE;
  }
  else if (lcp->seq_type != Seqtype_constraint_any)
  {
    rval = FALSE;
  }
  else if (lcp->partial5 != Partial_constraint_either)
  {
    rval = FALSE;
  } 
  else if (lcp->partial3 != Partial_constraint_either)
  {
    rval = FALSE;
  }
  return rval;
}


static Boolean DoesStrandMatchConstraint (SeqLocPtr slp, LocationConstraintPtr lcp)
{
  Uint2 strand;
  Boolean rval = FALSE;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->strand == Strand_constraint_any)
  {
    rval = TRUE;
  }
  else
  {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_minus)
    {
      if (lcp->strand == Strand_constraint_minus)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
    else
    {
      if (lcp->strand == Strand_constraint_plus)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
  }
  return rval;
}


static Boolean DoesBioseqMatchSequenceType (BioseqPtr bsp, Uint2 seq_type)
{
  Boolean rval = FALSE;

  if (bsp == NULL) return FALSE;
  if (seq_type == Seqtype_constraint_any) return TRUE;

  if (ISA_na (bsp->mol) && seq_type == Seqtype_constraint_nuc)
  {
    rval = TRUE;
  }
  else if (ISA_aa (bsp->mol) && seq_type == Seqtype_constraint_prot)
  {
    rval = TRUE;
  }
  return rval;
}


static Boolean DoesSequenceTypeMatchContraint (SeqLocPtr slp, LocationConstraintPtr lcp)
{
  Boolean   rval = FALSE;
  BioseqPtr bsp;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->seq_type == Seqtype_constraint_any)
  {
    rval = TRUE;
  }
  else
  {
    bsp = BioseqFindFromSeqLoc (slp);
    rval = DoesBioseqMatchSequenceType (bsp, lcp->seq_type);
  }
  return rval;
}


static Boolean DoesLocationMatchPartialnessConstraint (SeqLocPtr slp, LocationConstraintPtr lcp)
{
  Boolean rval = FALSE;
  Boolean partial5, partial3;

  if (slp == NULL) 
  {
    rval = FALSE;
  }
  else if (lcp == NULL)
  {
    rval = TRUE;
  }
  else
  {
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (lcp->partial5 == Partial_constraint_partial && !partial5) 
    {
      rval = FALSE;
    }
    else if (lcp->partial5 == Partial_constraint_complete && partial5)
    {
      rval = FALSE;
    }
    else if (lcp->partial3 == Partial_constraint_partial && !partial3) 
    {
      rval = FALSE;
    } 
    else if (lcp->partial3 == Partial_constraint_complete && partial3) 
    {
      rval = FALSE;
    }
    else
    {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean DoesLocationMatchConstraint (SeqLocPtr slp, LocationConstraintPtr lcp)

{
  Boolean rval = FALSE;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }  
  else if (IsLocationConstraintEmpty(lcp)) 
  {
    rval = TRUE;
  }
  else if (DoesStrandMatchConstraint (slp, lcp)
           && DoesSequenceTypeMatchContraint (slp, lcp) 
           && DoesLocationMatchPartialnessConstraint (slp, lcp))
  {
    rval = TRUE;
  }
  return rval; 
}


static Boolean DoesFeatureMatchLocationConstraint (SeqFeatPtr sfp, LocationConstraintPtr constraint)
{
  BioseqPtr bsp;
  SeqFeatPtr cds;
  SeqMgrFeatContext context;
  Boolean           rval = TRUE;

  if (sfp == NULL) {
    return FALSE;
  } else if (IsLocationConstraintEmpty (constraint)) {
    return TRUE;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (constraint->strand != Strand_constraint_any) {
    if (bsp == NULL) {
      rval = FALSE;
    } else if (ISA_aa (bsp->mol)) {
      cds = SeqMgrGetCDSgivenProduct (bsp, &context);
      if (cds == NULL) {
        rval = FALSE;
      } else if (!DoesStrandMatchConstraint (cds->location, constraint)) {
        rval = FALSE;
      }
    } else {
      if (!DoesStrandMatchConstraint (sfp->location, constraint)) {
        rval = FALSE;
      }
    }
  }

  if (!DoesBioseqMatchSequenceType (bsp, constraint->seq_type)) {
    rval = FALSE;
  }

  if (!DoesLocationMatchPartialnessConstraint (sfp->location, constraint)) {
    rval = FALSE;
  }
  return rval;
}


static Boolean DoesObjectMatchLocationConstraint (Uint1 choice, Pointer data, LocationConstraintPtr constraint)
{
  SeqFeatPtr  sfp;
  SeqDescrPtr sdp;
  CGPSetPtr   cgp;
  BioseqPtr  bsp = NULL;
  BioseqSetPtr bssp;
  ValNodePtr    vnp;
  ObjValNodePtr ovp;
  SeqMgrFeatContext context;

  if (data == NULL) return FALSE;
 
  if (IsLocationConstraintEmpty(constraint)) {
    return TRUE;
  }

  if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    return DoesFeatureMatchLocationConstraint (sfp, constraint);
  } else if (choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) data;
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) ovp->idx.parentptr;
        if (bssp != NULL && bssp->seq_set != NULL && IS_Bioseq_set (bssp->seq_set)) {
          bsp = (BioseqPtr) bssp->seq_set->data.ptrvalue;
        }
      } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        bsp = (BioseqPtr) ovp->idx.parentptr;
      }
    }
    if (!DoesBioseqMatchSequenceType(bsp, constraint->seq_type)) {
      return FALSE;
    }
    if (constraint->strand != Strand_constraint_any
        || constraint->partial5 != Partial_constraint_either
        || constraint->partial3 != Partial_constraint_either) {
      if (ISA_aa (bsp->mol)) {
        sfp = SeqMgrGetCDSgivenProduct (bsp, &context);
        if (sfp == NULL) {
          return FALSE;
        } else if (!DoesLocationMatchPartialnessConstraint (sfp->location, constraint)) {
          return FALSE;
        } else if (!DoesStrandMatchConstraint (sfp->location, constraint)) {
          return FALSE;
        } else {
          return TRUE;
        }
      } else {
        return FALSE;
      }
    } else {
      return TRUE;
    }
  } else if (choice == 0) {
    if (constraint->seq_type != Seqtype_constraint_any) {
      return FALSE;
    }
    cgp = (CGPSetPtr) data;
    for (vnp = cgp->cds_list; vnp != NULL; vnp = vnp->next) {
      if (DoesFeatureMatchLocationConstraint (vnp->data.ptrvalue, constraint)) {
        return TRUE;
      }
    }
    for (vnp = cgp->gene_list; vnp != NULL; vnp = vnp->next) {
      if (DoesFeatureMatchLocationConstraint (vnp->data.ptrvalue, constraint)) {
        return TRUE;
      }
    }
    for (vnp = cgp->mrna_list; vnp != NULL; vnp = vnp->next) {
      if (DoesFeatureMatchLocationConstraint (vnp->data.ptrvalue, constraint)) {
        return TRUE;
      }
    }
    for (vnp = cgp->prot_list; vnp != NULL; vnp = vnp->next) {
      if (DoesFeatureMatchLocationConstraint (vnp->data.ptrvalue, constraint)) {
        return TRUE;
      }
    }
    return FALSE;
  } else {
    return FALSE;
  }
}


/* for parsing and editing */
NLM_EXTERN CharPtr GetTextPortionFromString (CharPtr str, TextPortionPtr text_portion)
{
  CharPtr portion = NULL;
  CharPtr found_start, found_end;
  Int4    found_len;

  if (StringHasNoText (str)) {
    return NULL;
  }
  if (text_portion == NULL) {
    return StringSave (str);
  }  
  
  if (text_portion->left_text == NULL || text_portion->left_text [0] == 0)
  {
    found_start = str;
  }
  else
  {
    if (text_portion->case_sensitive)
    {
      found_start = StringSearch (str, text_portion->left_text);
    }
    else
    {
      found_start = StringISearch (str, text_portion->left_text);
    }
    
    if (text_portion->whole_word && ! IsWholeWordMatch (str, found_start, StringLen (text_portion->left_text)))
    {
      found_start = NULL;
    }
  }
  
  if (found_start == NULL)
  {
    return NULL;
  }
  
  if (!text_portion->include_left)
  {
    found_start += StringLen (text_portion->left_text);
  }
  
  if (text_portion->right_text == NULL || text_portion->right_text [0] == 0)
  {
    found_len = StringLen (found_start);
  }
  else
  {
    if (text_portion->case_sensitive)
    {
      found_end = StringSearch (found_start, text_portion->right_text);
    }
    else
    {
      found_end = StringISearch (found_start, text_portion->right_text);
    }
    if (text_portion->whole_word && ! IsWholeWordMatch (str, found_end, StringLen (text_portion->right_text)))
    {
      found_end = NULL;
    }    
    
    if (found_end == NULL)
    {
      found_len = 0;
    }
    else if (text_portion->include_right)
    {
      found_len = (Int4)(found_end - found_start) + StringLen (text_portion->right_text);
    }
    else
    {
      found_len = found_end - found_start;
    }
  }

  if (found_len > 0)
  {
    portion = (CharPtr) MemNew (sizeof (Char) * (found_len + 1));
    StringNCpy (portion, found_start, found_len);
    portion[found_len] = 0;
  }
  return portion;
}



static CharPtr FindTextPortionLocationInString (CharPtr str, TextPortionPtr text_portion)
{
  CharPtr start, stop;

  if (str == NULL || text_portion == NULL) return FALSE;

  if (text_portion->left_text != NULL) {
    start = StringSearch (str, text_portion->left_text);
    if (start != NULL) {
      if (!text_portion->include_left) {
        start += StringLen (text_portion->left_text);
      }
    }
  } else {
    start = str;
  }
  if (start != NULL) {
    if (text_portion->right_text != NULL) { 
      stop = StringSearch (start, text_portion->right_text);
      if (stop == NULL) {
        start = NULL;
      }
    }
  }
  return start;
}


static void ReplaceStringForParse(CharPtr src_text, TextPortionPtr text_portion)
{
  CharPtr         src, dst;
  
  if (src_text == NULL || text_portion == NULL) {
    return;
  }

  dst = FindTextPortionLocationInString (src_text, text_portion);
  if (dst == NULL) return;
  if (text_portion->right_text == NULL) {
    *dst = 0;
  } else {
    src = StringSearch (src_text, text_portion->right_text);
    if (src != NULL) {
      if (text_portion->include_right) {
        src += StringLen (text_portion->right_text);
      }
      while (*src != 0) {
        *dst = *src;
        dst++;
        src++;
      }
      *dst = 0;
    }
  }
}


/* generic functions for setting field values */
NLM_EXTERN Boolean SetStringValue (CharPtr PNTR existing_val, CharPtr new_val, Uint2 existing_text)
{
  Boolean rval = FALSE;
  Int4 len;
  CharPtr tmp;

  if (existing_val == NULL) {
    return FALSE;
  }

  if (StringHasNoText (*existing_val)) {
    *existing_val = MemFree (*existing_val);
    *existing_val = StringSave (new_val);
    rval = TRUE;
  } else {
    switch (existing_text) {
      case ExistingTextOption_replace_old :
        *existing_val = MemFree (*existing_val);
        *existing_val = StringSave (new_val);
        rval = TRUE;
        break;
      case ExistingTextOption_append_semi :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s; %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_space :
        len = StringLen (new_val) + StringLen (*existing_val) + 2;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_colon :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s: %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_none :
        len = StringLen (new_val) + StringLen (*existing_val) + 1;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s%s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_semi :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s; %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_space :
        len = StringLen (new_val) + StringLen (*existing_val) + 2;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_colon :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s: %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_none :
        len = StringLen (new_val) + StringLen (*existing_val) + 1;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s%s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_leave_old :
        rval = FALSE;
    }
  }
  return rval;
}


/* NOTE: The following functions, GetTwoFieldSubfield, SetTwoFieldSubfield, and RemoveTwoFieldSubfield,
 * all assume that if only one field is present, it is subfield 1.
 */
static CharPtr GetTwoFieldSubfield (CharPtr str, Uint1 subfield)
{
  CharPtr cp;
  CharPtr new_val = NULL;
  Int4    len;

  if (StringHasNoText (str) || subfield > 2) {
    return NULL;
  }
  if (subfield == 0) {
    new_val = StringSave (str);
  } else {
    cp = StringChr (str, ':');
    if (cp == NULL) {
      if (subfield == 1) {
        new_val = StringSave (str);
      } else {
        new_val = NULL;
      }
    } else {
      if (subfield == 1) {
        len = cp - str + 1;
        new_val = (CharPtr) MemNew (sizeof (Char) * len);
        StringNCpy (new_val, str, len - 1);
        new_val[len - 1] = 0;
      } else if (!StringHasNoText (cp + 1)) {
        new_val = StringSave (cp + 1);
      }
    }
  }
  return new_val;
}


static CharPtr MakeValFromTwoFields (CharPtr PNTR fields)
{
  Boolean empty1, empty2;
  CharPtr val = NULL;

  if (fields == NULL) return NULL;

  empty1 = StringHasNoText (fields[0]);
  empty2 = StringHasNoText (fields[1]);
  if (empty1 && empty2) {
    val = NULL;
  } else if (empty1) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[1]) + 2));
    sprintf (val, ":%s", fields[1]);
  } else if (empty2) {
    val = StringSave (fields[0]);
  } else {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[0]) + StringLen (fields[1]) + 2));
    sprintf (val, "%s:%s", fields[0], fields[1]);
  }
  return val;
}

  
static Boolean RemoveTwoFieldSubfield (CharPtr PNTR existing_val, Uint1 subfield)
{
  Boolean rval = FALSE;
  CharPtr fields[2];

  if (existing_val == NULL || StringHasNoText (*existing_val) || subfield > 2) {
    return FALSE;
  }
  if (subfield == 0) {
    *existing_val = MemFree (*existing_val);
    rval = TRUE;
  } else {
    fields[0] = GetTwoFieldSubfield (*existing_val, 1);
    fields[1] = GetTwoFieldSubfield (*existing_val, 2);
    if (!StringHasNoText (fields[subfield - 1])) {
      fields[subfield - 1] = MemFree (fields[subfield - 1]);
      *existing_val = MemFree (*existing_val);
      *existing_val = MakeValFromTwoFields (fields);
      rval = TRUE;
    }
    fields[0] = MemFree (fields[0]);
    fields[1] = MemFree (fields[1]);
  }
  return rval;
}


static Boolean SetTwoFieldSubfield (CharPtr PNTR existing_val, Int4 subfield, CharPtr new_field, Uint2 existing_text)
{
  Boolean rval = FALSE;
  CharPtr fields[2];

  if (existing_val == NULL || subfield > 2 || StringHasNoText (new_field)) {
    return FALSE;
  }
  if (subfield == 0) {
    rval = SetStringValue (existing_val, new_field, existing_text);
  } else {
    fields[0] = GetTwoFieldSubfield (*existing_val, 1);
    fields[1] = GetTwoFieldSubfield (*existing_val, 2);
    if (SetStringValue (&(fields[subfield - 1]), new_field, existing_text)) {
      *existing_val = MemFree (*existing_val);
      *existing_val = MakeValFromTwoFields (fields);
      rval = TRUE;
    }
    fields[0] = MemFree (fields[0]);
    fields[1] = MemFree (fields[1]);
  }
  return rval;
}


/* NOTE: The following functions, GetThreeFieldSubfield, SetThreeFieldSubfield, and RemoveThreeFieldSubfield
 * all assume that if only one field is present, it is subfield 3.  If two fields are present, they are subfields 1 and 3.
 */
static CharPtr GetThreeFieldSubfield (CharPtr str, Uint1 subfield)
{
  CharPtr cp, cp2;
  Int4    num_colons = 0;
  CharPtr new_val = NULL;

  if (StringHasNoText (str)) {
    return NULL;
  }

  cp = StringChr (str, ':');
  while (cp != NULL) {
    num_colons ++;
    cp = StringChr (cp + 1, ':');
  }

  if (subfield == 0) {
    new_val = StringSave (str);
  } else if (subfield == 1) {
    if (num_colons == 0) {
      return NULL;
    } else {
      cp = StringChr (str, ':');
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp - str + 1));
      StringNCpy (new_val, str, cp - str);
      new_val[cp - str] = 0;
    }
  } else if (subfield == 2) {
    if (num_colons == 0 || num_colons == 1) {
      return NULL;
    } else {
      cp = StringChr (str, ':');
      cp2 = StringChr (cp + 1, ':');
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp2 - cp));
      StringNCpy (new_val, cp + 1, cp2 - cp - 1);
      new_val[cp2 - cp - 1] = 0;
    }
  } else {
    if (num_colons == 0) {
      new_val = StringSave (str);
    } else {
      cp = StringRChr (str, ':');
      new_val = StringSave (cp + 1);
    }
  }
  return new_val;
}


static CharPtr MakeValFromThreeFields (CharPtr PNTR fields)
{
  Int4 i;
  Boolean empty[3];
  CharPtr val = NULL;

  if (fields == NULL) return NULL;

  for (i = 0; i < 3; i++) {
    empty[i] = StringHasNoText (fields[i]);
  }


  if (empty[0] && empty[1] && empty[2]) {
    /* do nothing, value is now empty */
  } else if (empty[0] && empty[1]) {
    val = StringSave (fields[2]);
  } else if (empty[0] && empty[2]) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[1]) + 2));
    sprintf (val, ":%s:", fields[1]);
  } else if (empty[1] && empty[2]) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[1]) + 2));
    sprintf (val, "%s:", fields[0]);
  } else if (empty[0]) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[1]) + StringLen (fields[2]) + 3));
    sprintf (val, ":%s:%s", fields[1], fields[2]);
  } else if (empty[1]) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[0]) + StringLen (fields[2]) + 3));
    sprintf (val, "%s:%s", fields[0], fields[2]);
  } else if (empty[2]) {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[0]) + StringLen (fields[1]) + 3));
    sprintf (val, "%s:%s:", fields[0], fields[1]);
  } else {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (fields[0]) + StringLen (fields[1]) + StringLen (fields[2]) + 3));
    sprintf (val, "%s:%s:%s", fields[0], fields[1], fields[2]);
  }
  return val;
}


static Boolean RemoveThreeFieldSubfield (CharPtr PNTR existing_val, Uint1 subfield)
{
  Int4    i;
  CharPtr fields[3];
  Boolean rval = FALSE;

  if (existing_val == NULL || subfield > 3 || StringHasNoText (*existing_val)) return FALSE;

  if (subfield == 0) {
    *existing_val = MemFree (*existing_val);
    rval = TRUE;
  } else {
    for (i = 0; i < 3; i++) {
      fields[i] = GetThreeFieldSubfield (*existing_val, i + 1);
    }
    if (!StringHasNoText (fields[subfield - 1])) {
      fields[subfield - 1] = MemFree (fields[subfield - 1]);
      *existing_val = MakeValFromThreeFields (fields);
      rval = TRUE;
    }
    for (i = 0; i < 3; i++) {
      fields[i] = MemFree (fields[i]);
    }
  }
  return rval;
}


static Boolean SetThreeFieldSubfield (CharPtr PNTR existing_val, Int4 subfield, CharPtr new_field, Uint2 existing_text)
{
  Int4    i;
  CharPtr fields[3];
  Boolean rval = FALSE;

  if (existing_val == NULL || StringHasNoText (new_field) || subfield < 0 || subfield > 3) return FALSE;

  if (subfield == 0) {
    rval = SetStringValue (existing_val, new_field, existing_text);
  } else {
    for (i = 0; i < 3; i++) {
      fields[i] = GetThreeFieldSubfield (*existing_val, i + 1);
    }
    if (SetStringValue (&(fields[subfield - 1]), new_field, existing_text)) {
      *existing_val = MemFree (*existing_val);
      *existing_val = MakeValFromThreeFields (fields);
      rval = TRUE;
    }
    for (i = 0; i < 3; i++) {
      fields[i] = MemFree (fields[i]);
    }
  }
  return rval;
}


static Boolean 
SetStringsInValNodeStringList 
(ValNodePtr PNTR list,
 StringConstraintPtr scp,
 CharPtr new_val,
 Uint2   existing_text)
{
  ValNodePtr vnp;
  CharPtr    cp;
  Boolean rval = FALSE;
  
  if (list == NULL)
  {
    return FALSE;
  }

  if (*list == NULL && (scp == NULL || StringHasNoText (scp->match_text))) {
    ValNodeAddPointer (list, 0, StringSave (new_val));
    rval = TRUE;
  } else if (existing_text == ExistingTextOption_add_qual) {
      ValNodeAddPointer (list, 0, StringSave (new_val));
      rval = TRUE;
  } else if (existing_text == ExistingTextOption_replace_old) {
    if (DoesStringListMatchConstraint (*list, scp)) {
      *list = ValNodeFreeData (*list);
      vnp = ValNodeNew (NULL);
      vnp->data.ptrvalue = StringSave (new_val);
      *list = vnp;
      rval = TRUE;
    }
  } else if (existing_text == ExistingTextOption_leave_old) {
    rval = FALSE;
  } else {
    for (vnp = *list; vnp != NULL; vnp = vnp->next)
    {
      cp = (CharPtr) vnp->data.ptrvalue;
      if (DoesStringMatchConstraint (cp, scp)) {
        rval |= SetStringValue (&cp, new_val, existing_text);
        vnp->data.ptrvalue = cp;
      }
    }
  }
  return rval;
}


static Boolean SetStringInGBQualList (GBQualPtr PNTR list, ValNodePtr field, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  Boolean rval = FALSE, does_match;
  Int4 gbqual, subfield;
  CharPtr qual_name = NULL, tmp;
  GBQualPtr gbq, last_gbq = NULL;

  if (field == NULL) return FALSE;

  if (field->choice == FeatQualChoice_legal_qual) 
  {
    gbqual = GetGBQualFromFeatQual (field->data.intvalue, &subfield);
    if (gbqual > -1) {
      qual_name = ParFlat_GBQual_names [gbqual].name;
      if (existing_text == ExistingTextOption_add_qual) {
        gbq = GBQualNew ();
        gbq->qual = StringSave (qual_name);
        gbq->val = StringSave (new_val);
        if (last_gbq == NULL) {
          *list = gbq;
        } else {
          last_gbq->next = gbq;
        }
        rval = TRUE;
      } else {        
        for (gbq = *list; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, qual_name) == 0) {
            if (subfield > 0) {
              does_match = TRUE;
              if (!IsStringConstraintEmpty (scp)) {
                tmp = GetTwoFieldSubfield (gbq->val, subfield);
                does_match = DoesStringMatchConstraint (tmp, scp);
                tmp = MemFree (tmp);
              }
              if (does_match) {
                rval |= SetTwoFieldSubfield (&(gbq->val), subfield, new_val, existing_text);
              }
            } else if (DoesStringMatchConstraint (gbq->val, scp)) {
              rval |= SetStringValue (&(gbq->val), new_val, existing_text);
            }
          }
          last_gbq = gbq;
        }
        if (!rval && (scp == NULL || scp->match_text == NULL)) {
          gbq = GBQualNew ();
          gbq->qual = StringSave (qual_name);
          gbq->val = StringSave (new_val);
          if (last_gbq == NULL) {
            *list = gbq;
          } else {
            last_gbq->next = gbq;
          }
          rval = TRUE;
        }
      }
    }
  } else if (field->choice == FeatQualChoice_illegal_qual) {
    for (gbq = *list; gbq != NULL; gbq = gbq->next) {
      if (DoesStringMatchConstraint (gbq->qual, field->data.ptrvalue)
          && DoesStringMatchConstraint (gbq->val, scp)) {
        rval |= SetStringValue (&(gbq->val), new_val, existing_text);
      }
    }
  }

  return rval;
}


static Boolean IsAllDigits (CharPtr str)
{
  CharPtr cp;

  if (StringHasNoText (str)) return FALSE;

  cp = str;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (*cp == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean SetInt2ValueWithString (Int2Ptr val, CharPtr val_str, Uint2 existing_text)
{
  Char    num[15];
  CharPtr tmp = NULL;
  Boolean rval = FALSE;

  if (val == NULL) return FALSE;

  sprintf (num, "%d", *val);
  tmp = StringSave (num);
  if (SetStringValue (&tmp, val_str, existing_text)
      && IsAllDigits (tmp)) {
    *val = atoi (tmp);
    rval = TRUE;
  }
  tmp = MemFree (tmp);
  return rval;
}


static CharPtr GetInt2ValueFromString (Int2 val, StringConstraintPtr scp)
{
  Char num[15];

  sprintf (num, "%d", val);
  if (DoesStringMatchConstraint (num, scp)) {
    return StringSave (num);
  } else {
    return NULL;
  }
}


static Boolean SetObjectIdString (ObjectIdPtr oip, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  Char    num[15];
  CharPtr tmp = NULL;

  if (oip == NULL) {
    return FALSE;
  }

  if (oip->id > 0) {
    sprintf (num, "%d", oip->id);
    tmp = StringSave (num);
  } else {
    tmp = StringSaveNoNull (oip->str);
  }
  if (SetStringValue (&tmp, value, existing_text)) {
    oip->str = MemFree (oip->str);        
    oip->id = 0;
    if (IsAllDigits (tmp)) {
      oip->id = atoi (tmp);
    } else {
      oip->str = tmp;
      tmp = NULL;
    }
    rval = TRUE;
  }
  tmp = MemFree (tmp);
  return rval;
}


static CharPtr GetObjectIdString (ObjectIdPtr oip)
{
  CharPtr rval = NULL;
  Char    num[15];

  if (oip == NULL) {
    return FALSE;
  }

  if (oip->id > 0) {
    sprintf (num, "%d", oip->id);
    rval = StringSave (num);
  } else {
    rval = StringSaveNoNull (oip->str);
  }
  return rval;
}


static Boolean DoesObjectIdMatchStringConstraint (ObjectIdPtr oip, StringConstraintPtr scp)
{
  Char    num[15];
  Boolean rval = FALSE;

  if (oip == NULL) {
    return FALSE;
  } else if (IsStringConstraintEmpty (scp)) {
    return TRUE;
  } else if (oip->id > 0) {
    sprintf (num, "%d", oip->id);
    rval = DoesStringMatchConstraint (num, scp);
  } else {
    rval = DoesStringMatchConstraint (oip->str, scp);
  }
  return rval;
}


/* generic functions for getting string values */
static Int4 GetDbtagStringLen (DbtagPtr db_tag)
{
  Int4 len;
  
  if (db_tag == NULL)
  {
    return 0;
  }
  
  len = StringLen (db_tag->db) + 2;
  if (db_tag->tag != NULL)
  {
    if (db_tag->tag->str != NULL)
    {
      len += StringLen (db_tag->tag->str);
    }
    else
    {
      len += 10;
    }
  }
  return len;
}


static CharPtr GetDbtagString (DbtagPtr db_tag)
{
  Int4    len;
  CharPtr str;
  
  if (db_tag == NULL) {
    return NULL;
  }
  
  len = GetDbtagStringLen (db_tag);
  if (len == 0) {
    return NULL;
  }
  
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str != NULL) {
    StringCpy (str, db_tag->db);
    StringCat (str, ":");
    if (db_tag->tag != NULL) {
      if (db_tag->tag->str != NULL) {
        StringCat (str, db_tag->tag->str);
      } else {
        sprintf (str + StringLen (str), "%d", db_tag->tag->id);
      }
    }
  }
  return str;
}


static Boolean SetDbtagString (DbtagPtr db_tag, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  CharPtr cp;
  Int4    dbxvalid;
  CharPtr tmp;
  CharPtr twoval;
  
  if (db_tag == NULL || StringHasNoText (value)) {
    return FALSE;
  }

  cp = StringChr (value, ':');
  if (cp == NULL) {
    tmp = StringSave (db_tag->db);
    if (SetStringValue (&tmp, value, existing_text)) {
      dbxvalid = DbxrefIsValid (tmp, NULL, NULL, NULL, NULL);
      if (dbxvalid != 0) {
        db_tag->db = MemFree (db_tag->db);
        db_tag->db = tmp;
        tmp = NULL;
        rval = TRUE;
      }
    }
    if (!rval) {
      if (db_tag->tag == NULL) {
        db_tag->tag = ObjectIdNew();
      }
      rval = SetObjectIdString (db_tag->tag, value, existing_text);
    }
    tmp = MemFree (tmp);
  } else {
    twoval = StringSave (value);
    cp = StringChr (twoval, ':');
    *cp = 0;
    cp++;
    rval = SetStringValue (&(db_tag->db), twoval, existing_text);
    if (db_tag->tag == NULL) {
      db_tag->tag = ObjectIdNew ();
    }
    rval |= SetObjectIdString (db_tag->tag, cp, existing_text);
    twoval = MemFree (twoval);
  }
  return rval;
}


static Boolean SetDbxrefString (ValNodePtr PNTR list, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  ValNodePtr vnp;
  Boolean    rval = FALSE, skip;
  DbtagPtr   dbtag;
  CharPtr    cp;
  
  if (list == NULL) {
    return FALSE;
  }

  if (existing_text == ExistingTextOption_add_qual
      || (*list == NULL && (scp == NULL || StringHasNoText (scp->match_text)))) {
    dbtag = DbtagNew ();
    rval = SetDbtagString (dbtag, value, existing_text);
    if (rval) {
      ValNodeAddPointer (list, 0, dbtag);
    } else {
      dbtag = DbtagFree (dbtag);
    }
  } else {
    for (vnp = *list; vnp != NULL; vnp = vnp->next) {
      skip = FALSE;
      if (scp != NULL) {
        cp = GetDbtagString (vnp->data.ptrvalue);
        if (!DoesStringMatchConstraint (cp, scp)) {
          skip = TRUE;
        }
        cp = MemFree (cp);
      }
      if (!skip) {
        rval |= SetDbtagString (vnp->data.ptrvalue, value, existing_text);
      }
    }
  }
  return rval;
}



static CharPtr GetFirstValNodeStringMatch (ValNodePtr vnp, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (vnp != NULL && str == NULL) {
    if (!StringHasNoText (vnp->data.ptrvalue)
        && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
      str = StringSave (vnp->data.ptrvalue);
    } 
    vnp = vnp->next;
  }
  return str;
}


static Boolean RemoveValNodeStringMatch (ValNodePtr PNTR list, StringConstraintPtr scp)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp;
  Boolean    rval = FALSE;

  if (list == NULL) return FALSE;
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    if (!StringHasNoText (vnp->data.ptrvalue) 
        && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFreeData (vnp);
      rval = TRUE;
    } else {
      vnp_prev = vnp;
    }
    vnp = vnp_next;
  }
  return rval;
}


static CharPtr GetFirstGBQualMatch (GBQualPtr qual, CharPtr qual_name, Int4 subfield, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (qual != NULL && str == NULL) {
    if (StringICmp (qual->qual, qual_name) == 0) {
      str = GetTwoFieldSubfield (qual->val, subfield);
      if (StringHasNoText (str) || !DoesStringMatchConstraint (str, scp)) {
        str = MemFree (str);
      }
    } 
    qual = qual->next;
  }
  return str;
}


static CharPtr GetFirstGBQualMatchConstraintName (GBQualPtr qual, StringConstraintPtr qual_name, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (qual != NULL && str == NULL) {
    if (DoesStringMatchConstraint (qual->qual, qual_name)
        &&!StringHasNoText (qual->val)
        && DoesStringMatchConstraint (qual->val, scp)) {
      str = StringSave (qual->val);
    } 
    qual = qual->next;
  }
  return str;
}


static Boolean RemoveGBQualMatch (GBQualPtr PNTR list, CharPtr qual_name, Int4 subfield, StringConstraintPtr scp)
{
  GBQualPtr qual_prev = NULL, qual_next, qual;
  CharPtr   tmp;
  Boolean   rval = FALSE, does_match, do_remove;

  if (list == NULL) return FALSE;

  qual = *list;
  while (qual != NULL) {
    qual_next = qual->next;
    do_remove = FALSE;
    if (StringICmp (qual->qual, qual_name) == 0) {
      if (subfield > 0) {
        does_match = TRUE;
        if (!IsStringConstraintEmpty (scp)) {
          tmp = GetTwoFieldSubfield (qual->val, subfield);
          does_match = DoesStringMatchConstraint (tmp, scp);
          tmp = MemFree (tmp);
        }
        if (RemoveTwoFieldSubfield (&(qual->val), subfield)) {
          rval = TRUE;
          if (StringHasNoText (qual->val)) {
            do_remove = TRUE;
          }
        }
      } else if (DoesStringMatchConstraint (qual->val, scp)) {
        do_remove = TRUE;
      }
    } 
    if (do_remove) {
      if (qual_prev == NULL) {
        *list = qual->next;
      } else {
        qual_prev->next = qual->next;
      }
      qual->next = NULL;
      qual = GBQualFree (qual);
      rval = TRUE;
    } else {
      qual_prev = qual;
    }
    qual = qual_next;
  }
  return rval;
}


static Boolean RemoveGBQualMatchConstraintName (GBQualPtr PNTR list, StringConstraintPtr qual_name, StringConstraintPtr scp)
{
  GBQualPtr qual_prev = NULL, qual_next, qual;
  Boolean   rval = FALSE;

  if (list == NULL) return FALSE;
  qual = *list;
  while (qual != NULL) {
    qual_next = qual->next;
    if (DoesStringMatchConstraint (qual->qual, qual_name)
        && !StringHasNoText (qual->val) 
        && DoesStringMatchConstraint (qual->val, scp)) {
      if (qual_prev == NULL) {
        *list = qual->next;
      } else {
        qual_prev->next = qual->next;
      }
      qual->next = NULL;
      qual = GBQualFree (qual);
      rval = TRUE;
    } else {
      qual_prev = qual;
    }
    qual = qual_next;
  }
  return rval;
}


static CharPtr GetDbxrefString (ValNodePtr list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       len = 0;
  CharPtr    str = NULL, cp;
  
  if (list == NULL) {
    return NULL;
  }
  
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    cp = GetDbtagString (vnp->data.ptrvalue);
    if (cp != NULL && DoesStringMatchConstraint(cp, scp)) {
      len += StringLen (cp) + 1;
    }
    cp = MemFree (cp);
  }
  
  if (len == 0) {
    return NULL;
  }
  
  str = (CharPtr) MemNew ((len + 1) * sizeof (Char));
  if (str != NULL) {
    for (vnp = list; vnp != NULL; vnp = vnp->next) {
      cp = GetDbtagString (vnp->data.ptrvalue);
      if (cp != NULL && DoesStringMatchConstraint(cp, scp)) {
        StringCat (str, cp);
        StringCat (str, ";");
      }
      cp = MemFree (cp);
    }
  }
  if (StringLen (str) >1) {
    /* remove final semicolon */
    str [StringLen (str) - 2] = 0;
  }
  return str;
}


static Boolean RemoveDbxrefString (ValNodePtr PNTR list, StringConstraintPtr scp)
{
  ValNodePtr vnp, vnp_prev = NULL, vnp_next;
  CharPtr    cp;
  Boolean    rval = FALSE;
  
  if (list == NULL || *list == NULL) {
    return FALSE;
  }
  
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    cp = GetDbtagString (vnp->data.ptrvalue);
    if (DoesStringMatchConstraint(cp, scp)) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp->data.ptrvalue = DbtagFree (vnp->data.ptrvalue);
      vnp = ValNodeFree (vnp);
      rval = TRUE;
    } else {
      vnp_prev = vnp;
    }
  }
  return rval;  
}


NLM_EXTERN CharPtr GetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  RnaRefPtr  rrp;
  SeqMgrFeatContext context;
  CharPtr    str = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    str = GetFirstGBQualMatch (sfp->qual, "product", 0, scp);
  }

  if (str == NULL) {
    if (rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue)
        && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
        && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
        && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0) {
      str = StringSave (rrp->ext.value.ptrvalue);        
    } else if (rrp->ext.choice == 2 && rrp->ext.value.ptrvalue != NULL) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) != NULL
          && !StringHasNoText (context.label)
          && StringCmp (context.label, "tRNA") != 0) {
        str = (CharPtr) MemNew (sizeof (Char) + (StringLen (context.label) + 6));
        sprintf (str, "tRNA-%s", context.label);
      }
    }
    if (!DoesStringMatchConstraint(str, scp)) {
      str = MemFree (str);
    }
  }
  return str;
}


static Boolean IsParseabletRNAName (CharPtr name_string)
{
  if (StringHasNoText(name_string)) 
  {
    return TRUE;
  }
  else if (StringNICmp (name_string, "trna-", 5) != 0)
  {
    return FALSE;
  }
  else if (StringLen (name_string) != 8)
  {
    return FALSE;
  }
  else if (ParseTRnaString (name_string, NULL, NULL, TRUE) == 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static Boolean SetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  RnaRefPtr  rrp;
  Boolean rval = FALSE;
  ValNode vn;
  CharPtr cp, tmp;
  tRNAPtr trp;
  Boolean justTrnaText = FALSE;
  Uint1   codon [6];

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    vn.choice = FeatQualChoice_legal_qual;
    vn.data.intvalue = Feat_qual_legal_product;

    rval = SetStringInGBQualList (&(sfp->qual), &vn, scp, new_val, existing_text);
  }

  if (!rval) {
    if ((rrp->ext.choice == 0 || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue)))
        && (scp == NULL || scp->match_text == NULL)) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave (new_val);
      rrp->ext.choice = 1;
      rval = TRUE;
    } else if (rrp->ext.choice == 1 
                && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
                && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
                && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0
                && DoesStringMatchConstraint (rrp->ext.value.ptrvalue, scp)) {
      cp = rrp->ext.value.ptrvalue;
      rval = SetStringValue (&cp, new_val, existing_text);
      rrp->ext.value.ptrvalue = cp;
      rval = TRUE;
    } else if (rrp->ext.choice == 2) {
      tmp = GetRNAProductString (sfp, NULL);
      if (DoesStringMatchConstraint (tmp, scp)
          && SetStringValue (&tmp, new_val, existing_text)) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp == NULL) {
          trp = MemNew (sizeof (tRNA));
          trp->aatype = 0;
          MemSet (trp->codon, 255, sizeof (trp->codon));
          trp->anticodon = NULL;
          rrp->ext.value.ptrvalue = trp;
        }

        if (!IsParseabletRNAName(tmp))
        {
          if (trp->anticodon == NULL
              && trp->codon[0] == 255
              && trp->codon[1] == 255
              && trp->codon[2] == 255
              && trp->codon[3] == 255
              && trp->codon[4] == 255
              && trp->codon[5] == 255)
          {
            trp = MemFree (trp);
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = tmp;
            tmp = NULL;
            rval = TRUE;
          }
          else
          {
            vn.choice = FeatQualChoice_legal_qual;
            vn.data.intvalue = Feat_qual_legal_product;
            if (SetStringInGBQualList (&(sfp->qual), &vn, scp, new_val, existing_text)) {
              trp->aa = 0;
              rval = TRUE;
            }
          }
        }
        else
        {
          trp->aa = ParseTRnaString (tmp, &justTrnaText, codon, TRUE);
          trp->aatype = 2;
          rval = TRUE;
        }
        tmp = MemFree (tmp);
      }
    }
  }
  return rval;
}


static Boolean RemoveRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  RnaRefPtr  rrp;
  Boolean    rval = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    rval = RemoveGBQualMatch (&(sfp->qual), "product", 0, scp);
  }

  if (!rval 
      && rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue)
      && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
      && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
      && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0
      && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp)) {
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    rrp->ext.choice = 0;
    rval = TRUE;
  }
  return rval;
}


static Boolean RemovetRNACodons_Recognized (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->ext.choice != 2) {
    return FALSE;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) {
    return FALSE;
  }

  trp->codon [0] = 255;
  trp->codon [1] = 255;
  trp->codon [2] = 255;
  trp->codon [3] = 255;
  trp->codon [4] = 255;
  trp->codon [5] = 255;

  return TRUE;
}
      


static Boolean ParseCodonsRecognizedFromCommaDelimitedList (CharPtr str, Uint1Ptr codons)
{
  Int4 codon_num, k = 0, q;
  Char    ch;
  Boolean rval = TRUE;
  Uint1   codon[4];
  Uint1   code;

  if (StringHasNoText (str) || codons == NULL) {
    return FALSE;
  }

  for (codon_num = 0; codon_num < 6; codon_num++) {
    codons[codon_num] = 255;
  }
  codon_num = 0;

  while (isspace (*str)) {
    str++;
  }

  while (*str != 0 && codon_num < 6 && rval) {
    k = 0;
    q = 0;
    ch = str [k];
    while (ch != '\0' && q < 3 && rval) {
      ch = TO_UPPER (ch);
      if (StringChr ("ACGTU", ch) != NULL) {
        if (ch == 'U') {
          ch = 'T';
        }
        codon [q] = (Uint1) ch;
        q++;
      } else {
        rval = FALSE;
      }
      k++;
      ch = str [k];
    }
    if (q < 3 || isalpha (ch)) {
      rval = FALSE;
    }
    if (rval) {
      codon [q] = 0;
      if (q == 3) {
        code = IndexForCodon (codon, Seq_code_iupacna);
        if (code == INVALID_RESIDUE) {
          rval = FALSE;
        } else {
          codons [codon_num++] = code;
        }
      }
      str += 3;
      while (isspace (*str)) {
        str++;
      }
      while (*str == ',') {
        str++;
      }
      while (isspace (*str)) {
        str++;
      }
    }
  }
  if (*str != 0) {
    rval = FALSE;
  }
  return rval;
}


static Boolean SettRNACodons_Recognized (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Uint1     codon[6];
  Uint1     new_codons[6];
  Int4      codon_num, num_new, num_old = 0, i;
  Boolean   rval = FALSE, already_have;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }
  if (StringHasNoText (new_val)) {
    return FALSE;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->ext.choice != 2) {
    return FALSE;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) {
    return FALSE;
  }

  if (ParseCodonsRecognizedFromCommaDelimitedList (new_val, codon)) {
    switch (existing_text) {
      case ExistingTextOption_replace_old :
        for (codon_num = 0; codon_num < 6; codon_num++) {
          trp->codon[codon_num] = codon[codon_num];
        }
        rval = TRUE;
        break;
      case ExistingTextOption_append_semi :        
      case ExistingTextOption_append_space :
      case ExistingTextOption_append_colon :
      case ExistingTextOption_append_none :
      case ExistingTextOption_prefix_semi :
      case ExistingTextOption_prefix_space :
      case ExistingTextOption_prefix_colon :
      case ExistingTextOption_prefix_none :
      case ExistingTextOption_add_qual :
        for (num_old = 0; num_old < 6 && trp->codon[num_old] != 255; num_old++) {
          new_codons[num_old] = trp->codon[num_old];
        }
        codon_num = num_old;
        rval = TRUE;
        for (num_new = 0; num_new < 6 && codon[num_new] != 255 && rval; num_new++) {
          already_have = FALSE;
          for (i = 0; i < codon_num && !already_have; i++) {
            if (codon[num_new] == new_codons[i]) {
              already_have = TRUE;
            }
          }
          if (!already_have) {
            if (codon_num < 6) {
              new_codons[codon_num] = codon[num_new];
              codon_num++;
            } else {
              rval = FALSE;
            }
          }
        }
        if (rval) {
          for (i = 0; i < codon_num; i++) {
            trp->codon[i] = new_codons[i];
          }
          while (codon_num < 6) {
            trp->codon[codon_num++] = 255;
          }
        }
        break;
      case ExistingTextOption_leave_old :
        if (trp->codon[0] == 255) {
          for (i = 0; i < 6; i++) {
            trp->codon[i] = codon[i];
          }
          rval = TRUE;
        }
        break;
    }
  }
  return TRUE;
}
      

static CharPtr GettRNACodonsRecognized (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Int4      j;
  Char      buf[31];
  Uint1     codon [4];

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->ext.choice != 2) {
    return NULL;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) {
    return NULL;
  }

  buf[0] = 0;
  
  for (j = 0; j < 6; j++) {
    if (trp->codon [j] < 64) {
			/* Note - it is important to set the fourth character in the codon array to NULL
				* because CodonForIndex only fills in the three characters of actual codon,
				* so if you StringCpy the codon array and the NULL character is not found after
				* the three codon characters, you will write in memory you did not intend to.
				*/
			codon [3] = 0;
      if (CodonForIndex (trp->codon [j], Seq_code_iupacna, codon)) {
        if (buf[0] != 0) {
          StringCat (buf, ", ");
        }
        StringCat (buf, (CharPtr) codon);
      }
    }
  }
  if (buf[0] == 0) {
    return NULL;
  } else {
    return StringSave (buf);
  }
}


static SeqLocPtr ParseSimpleInterval (CharPtr str, BioseqPtr bsp, CharPtr PNTR end)
{
  Boolean partial_left = FALSE, partial_right = FALSE;
  Int4    left_num, right_num, swap_num;
  SeqLocPtr slp = NULL;
  Uint1     strand = Seq_strand_plus;

  if (StringHasNoText (str)) {
    return NULL;
  }

  while (isspace (*str)) {
    str++;
  }
  if (*str == '<' || *str == '>') {
    partial_left = TRUE;
    str++;
  }
  if (!isdigit (*str)) {
    return NULL;
  }
  left_num = atoi (str);
  while (isdigit (*str)) {
    str++;
  }
  while (isspace (*str) || *str == '.' || *str == '-') {
    str++;
  }
  if (*str == '<' || *str == '>') {
    partial_right = TRUE;
    str++;
  }
  if (!isdigit (*str)) {
    return NULL;
  }

  right_num = atoi (str);
  while (isdigit (*str)) {
    str++;
  }

  if (left_num > right_num) {
    swap_num = left_num;
    left_num = right_num;
    right_num = swap_num;
    strand = Seq_strand_minus;
  }

  slp = SeqLocIntNew (left_num - 1, right_num - 1, strand, SeqIdDup (SeqIdFindWorst (bsp->id)));
  SetSeqLocPartial (slp, partial_left, partial_right);

  if (end != NULL) {
    *end = str;
  }
  return slp;
} 


static void ComplementSeqLoc (SeqLocPtr slp)
{
  SeqIntPtr sip;
  Boolean   partial5 = FALSE, partial3 = FALSE;

  if (slp != NULL && slp->choice == SEQLOC_INT && slp->data.ptrvalue != NULL) {
    sip = (SeqIntPtr) slp->data.ptrvalue;
    if (sip->strand != Seq_strand_minus) {
      CheckSeqLocForPartial (slp, &partial5, &partial3);
      SetSeqLocPartial (slp, partial3, partial5);
      sip->strand = Seq_strand_minus;
    }
  }
}

  
static SeqLocPtr ParseSimpleSeqLoc (CharPtr str, BioseqPtr bsp)
{
  CharPtr cp, cp_next;
  SeqLocPtr slp = NULL, slp_first = NULL, slp_last = NULL, slp_tmp;
  Boolean is_complement = FALSE;

  if (StringHasNoText (str) || bsp == NULL) {
    return NULL;
  }

  cp = str;
  while (isspace (*cp)) {
    cp ++;
  }
  while (*cp != 0) {
    is_complement = FALSE;
    if (StringNICmp (cp, "complement", 10) == 0) {
      cp += 10;
      is_complement = TRUE;
    } else if (StringNICmp (cp, "comp", 4) == 0) {
      cp += 4;
      is_complement = TRUE;
    }
    if (*cp == '(') {
      cp++;
    }
    slp_tmp = ParseSimpleInterval (cp, bsp, &cp_next);
    if (slp_tmp == NULL) {
      slp = SeqLocFree (slp);
      return NULL;
    }
    if (is_complement) {
      ComplementSeqLoc (slp_tmp);
    }
    if (slp == NULL) {
      slp = slp_tmp;
    } else if (slp->choice == SEQLOC_INT) {
      slp_first = slp;
      slp_first->next = slp_tmp;
      slp = ValNodeNew (NULL);
      slp->choice = SEQLOC_MIX;
      slp->data.ptrvalue = slp_first;
    } else {
      ValNodeLink ((ValNodePtr PNTR) slp->data.ptrvalue, slp_tmp);
    }

    cp = cp_next;
    while (isspace (*cp)) {
      cp++;
    }
    if (*cp == ')') {
      cp++;
    }
    while (isspace (*cp)) {
      cp++;
    }
    if (*cp == ',') {
      cp++;
    }
    while (isspace (*cp)) {
      cp++;
    }
  }
  if (*cp != 0) {
    slp = SeqLocFree (slp);
  }
  return slp;
}


static Boolean SetAnticodon (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Boolean   rval = FALSE;
  SeqLocPtr slp, slp_merge;
  BioseqPtr bsp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }
  if (StringHasNoText (new_val)) {
    return FALSE;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->ext.choice != 2) {
    return FALSE;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) {
    return FALSE;
  }

  if (trp->anticodon != NULL && existing_text == ExistingTextOption_leave_old) {
    return FALSE;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) {
    return FALSE;
  }

  slp = ParseSimpleSeqLoc (new_val, bsp);
  if (slp == NULL) {
    return FALSE;
  }

  if (trp->anticodon == NULL) {
    trp->anticodon = slp;
    rval = TRUE;
  } else if (existing_text == ExistingTextOption_replace_old) {
    trp->anticodon = SeqLocFree (trp->anticodon);
    trp->anticodon = slp;
    rval = TRUE;
  } else {
    slp_merge = SeqLocMerge (bsp, trp->anticodon, slp, FALSE, FALSE, FALSE);
    slp = SeqLocFree (slp);
    trp->anticodon = SeqLocFree (trp->anticodon);
    trp->anticodon = slp_merge;
    rval = TRUE;
  }
  return rval;
}


static CharPtr GetIntervalString (SeqLocPtr slp)
{
  CharPtr fmt = "%s%d..%s%d";
  CharPtr complement_fmt = "complement(%s%d..%s%d)";
  CharPtr str = NULL;
  SeqIntPtr sip;
  Boolean   partial5 = FALSE, partial3 = FALSE;

  if (slp == NULL || slp->choice != SEQLOC_INT || slp->data.ptrvalue == NULL) {
    return NULL;
  }

  sip = (SeqIntPtr) slp->data.ptrvalue;

  CheckSeqLocForPartial (slp, &partial5, &partial3);

  if (sip->strand == Seq_strand_minus) {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (complement_fmt) + 30));
    sprintf (str, complement_fmt, partial3 ? "<" : "", 
                                  sip->from + 1,
                                  partial5 ? ">" : "",
                                  sip->to + 1);
  } else {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 30));
    sprintf (str, fmt, partial5 ? "<" : "", 
                        sip->from + 1,
                        partial3 ? ">" : "",
                        sip->to + 1);
  }
  return str;
}


static CharPtr GetAnticodonLocString (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Boolean   rval = FALSE;
  SeqLocPtr slp;
  CharPtr   str = NULL, tmp;
  ValNodePtr str_list = NULL, vnp;
  Int4       len = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->ext.choice != 2) {
    return NULL;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL || trp->anticodon == NULL) {
    return NULL;
  }

  if (trp->anticodon->choice == SEQLOC_INT) {
    str = GetIntervalString (trp->anticodon);
  } else if (trp->anticodon->choice == SEQLOC_MIX) {
    for (slp = trp->anticodon->data.ptrvalue; slp != NULL; slp = slp->next) {
      tmp = GetIntervalString (slp);
      if (tmp == NULL) {
        str_list = ValNodeFreeData (str_list);
        return StringSave ("complex location");
      } else {
        len += StringLen (tmp) + 2;
        ValNodeAddPointer (&str_list, 0, tmp);
      }
    }
    str = (CharPtr) MemNew (sizeof (Char) * len);
    str[0] = 0;
    for (vnp = str_list; vnp != NULL; vnp = vnp->next) {
      StringCat (str, vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        StringCat (str, ", ");
      }
    }
    str_list = ValNodeFreeData (str_list);
  }
  return str;
}



static SeqFeatPtr GetProtFeature (BioseqPtr protbsp)
{
  SeqMgrFeatContext fcontext;
  SeqAnnotPtr sap;
  SeqFeatPtr prot_sfp;
  ProtRefPtr prp;

  if (protbsp == NULL) return NULL;

  prot_sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &fcontext);
  if (prot_sfp == NULL) {
    sap = protbsp->annot;
    while (sap != NULL && prot_sfp == NULL) {
      if (sap->type == 1) {
        prot_sfp = sap->data;
        while (prot_sfp != NULL
               && (prot_sfp->data.choice != SEQFEAT_PROT
                   || (prp = prot_sfp->data.value.ptrvalue) == NULL
                   || prp->processed != 0)) {
          prot_sfp = prot_sfp->next;
        }
      }
      sap = sap->next;
    }
  }
  return prot_sfp;
}


static ProtRefPtr GetProtRefForFeature (SeqFeatPtr sfp)
{
  BioseqPtr  protbsp;
  SeqFeatPtr protsfp;
  ProtRefPtr prp = NULL;
  SeqFeatXrefPtr xref;

  if (sfp == NULL) return NULL;

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    xref = sfp->xref;
    while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
      xref = xref->next;
    }
    if (xref != NULL) {
      prp = xref->data.value.ptrvalue;
    }
    if (prp == NULL && sfp->product != NULL) {
      protbsp = BioseqFindFromSeqLoc (sfp->product);
      protsfp = GetProtFeature (protbsp);    
      if (protsfp != NULL) {
        prp = protsfp->data.value.ptrvalue;
      }
    }
  }
  return prp;
}


static void GetGeneInfoForFeature (SeqFeatPtr sfp, GeneRefPtr PNTR p_grp, SeqFeatPtr PNTR p_gene)
{
  GeneRefPtr grp = NULL;
  SeqFeatPtr gene = NULL;
  SeqMgrFeatContext fcontext;

  if (p_grp != NULL) {
    *p_grp = NULL;
  }
  if (p_gene != NULL) {
    *p_gene = NULL;
  }

  if (sfp == NULL) {
    return;
  }
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
    gene = sfp;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (gene != NULL) {
        grp = gene->data.value.ptrvalue;
      }
    } else if (SeqMgrGeneIsSuppressed (grp)) {
      grp = NULL;
    }
  }
  if (p_grp != NULL) {
    *p_grp = grp;
  }
  if (p_gene != NULL) {
    *p_gene = gene;
  }
}


static CharPtr GetCitationTextFromFeature (SeqFeatPtr sfp, StringConstraintPtr scp, ValNodePtr cit_list)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ValNodePtr  list = NULL, vnp;
  CharPtr     rval = NULL;
  Int4        serial_number;
  Char        buf[100];
  ValNodePtr  psp;

  if (sfp == NULL || sfp->cit == NULL) {
    return NULL;
  }

  bsp = GetSequenceForObject (OBJ_SEQFEAT, sfp);

  if (cit_list == NULL) {
    /* list not provided - must create now */
    sep = SeqMgrGetSeqEntryForData (bsp);
    list = GetCitListsForSeqEntry (sep);
    cit_list = list;
  } 

  psp = sfp->cit->data.ptrvalue;
  for (vnp = psp; vnp != NULL && rval == NULL; vnp = vnp->next) {
    
    serial_number = GetCitationNumberForMinPub (bsp, vnp, cit_list);
    if (serial_number > -1) {
      sprintf (buf, "%d", serial_number);
      if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (buf, scp)) {
        rval = StringSave (buf);
      }
    }
  }
  
  list = PubSerialNumberListFree (list);

  return rval;
}


NLM_EXTERN CharPtr GetQualFromFeatureEx (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  CharPtr   str = NULL;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  Int4      gbqual, subfield;
  SeqFeatPtr gene = NULL;
  CdRegionPtr crp;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return NULL;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return NULL;
  }

  // for gene fields
  GetGeneInfoForFeature (sfp, &grp, &gene);

  // for protein fields
  prp = GetProtRefForFeature (sfp);

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
    {
      str = StringSave (sfp->comment);
    }
  }
  /* db-xref */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue))))
  {
    str = GetDbxrefString (sfp->dbxref, scp);
  }
  /* exception */
  if (str == NULL 
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue))))
  {
    if (!StringHasNoText (sfp->except_text) && DoesStringMatchConstraint(sfp->except_text, scp))
    {
      str = StringSave (sfp->except_text);
    }
  }
  /* evidence */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue))))
  {
    if (sfp->exp_ev == 1)
    {
      str = StringSave ("experimental");
    }
    else if (sfp->exp_ev == 2)
    {
      str = StringSave ("non-experimental");
    }
    if (!DoesStringMatchConstraint(str, scp)) {
      str = MemFree (str);
    }
  }

  /* citation */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_citation)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("citation", field->field->data.ptrvalue))))
  {
    str = GetCitationTextFromFeature (sfp, scp, batch_extra == NULL ? NULL : batch_extra->cit_list);
  }

  /* fields common to some features */
  /* product */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue))))
  {
    if (prp != NULL) {
      str = GetFirstValNodeStringMatch (prp->name, scp);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      str = GetRNAProductString (sfp, scp);
    }
  }

  /* Gene fields */
  /* locus */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && DoesStringMatchConstraint(grp->locus, scp))
    {
      str = StringSave (grp->locus);
    }
  }
  /* description */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->desc) && DoesStringMatchConstraint(grp->desc, scp))
    {
      str = StringSave (grp->desc);
    }
  }
  /* maploc */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->maploc) && DoesStringMatchConstraint(grp->maploc, scp))
    {
      str = StringSave (grp->maploc);
    }
  }
  /* allele */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->allele) && DoesStringMatchConstraint(grp->allele, scp))
    {
      str = StringSave (grp->allele);
    }
  }
  /* locus_tag */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus_tag) && DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      str = StringSave (grp->locus_tag);
    }
  }
  /* synonym */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    str = GetFirstValNodeStringMatch (grp->syn, scp);
  }
  /* gene comment */
  if (str == NULL 
      && field->field->choice == FeatQualChoice_legal_qual
      && field->field->data.intvalue == Feat_qual_legal_gene_comment
      && gene != NULL
      && !StringHasNoText (gene->comment)
      && DoesStringMatchConstraint (gene->comment, scp)) {
    str = StringSave (gene->comment);
  }


  /* protein fields */
  /* note - product handled above */
  /* description */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (!StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
      str = StringSave (prp->desc);
    }
  }
  /* ec_number */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    str = GetFirstValNodeStringMatch (prp->ec, scp);
  }
  /* activity */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_activity)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    str = GetFirstValNodeStringMatch (prp->activity, scp);
  }
  
  /* tRNA qualifiers */
  /* codon-recognized */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_codons_recognized)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("codon-recognized", field->field->data.ptrvalue)))) {
    str = GettRNACodonsRecognized (sfp, scp);
  }
  /* anticodon */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_anticodon)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("anticodon", field->field->data.ptrvalue)))) {
    str = GetAnticodonLocString (sfp);
  }

  /* codon-start */
  if (field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_codon_start
      && sfp->data.choice == SEQFEAT_CDREGION) 
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp->frame == 1 || crp->frame == 0) {
      str = StringSave ("1");
    } else {
      str = (CharPtr) MemNew (sizeof (Char) * 15);
      sprintf (str, "%d", crp->frame);
    }
    if (!DoesStringMatchConstraint (str, scp)) {
      str = MemFree (str);
    }
  } 

  /* actual GenBank qualifiers */
  if (str == NULL)
  {
    if (field->field->choice == FeatQualChoice_legal_qual) 
    {
      gbqual = GetGBQualFromFeatQual (field->field->data.intvalue, &subfield);
      if (gbqual > -1) {
        str = GetFirstGBQualMatch (sfp->qual, ParFlat_GBQual_names [gbqual].name, subfield, scp);
      } else {
        /* need to do something with non-qualifier qualifiers */
      }
    } else {
      str = GetFirstGBQualMatchConstraintName (sfp->qual, field->field->data.ptrvalue, scp);
    }
  }
  return str;
}


NLM_EXTERN CharPtr GetQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp)
{
  return GetQualFromFeatureEx (sfp, field, scp, NULL);
}


static Boolean RemoveQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  RnaRefPtr  rrp;
  tRNAPtr trp;
  Int4      gbqual, subfield;
  SeqFeatPtr gene = NULL;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return FALSE;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return FALSE;
  }

  /* for gene fields */
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
    gene = sfp;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (gene != NULL) {
        grp = gene->data.value.ptrvalue;
      }
    } else if (SeqMgrGeneIsSuppressed (grp)) {
      grp = NULL;
    }
  }

  /* for protein fields */
  prp = GetProtRefForFeature (sfp);

  /* for RNA fields */
  if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  }

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->comment) && DoesStringMatchConstraint (sfp->comment, scp))
    {
      sfp->comment = MemFree (sfp->comment);
      rval = TRUE;
    }
  }
  /* db-xref */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue)))
  {
    rval = RemoveDbxrefString (&(sfp->dbxref), scp);
  }
  /* exception */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->except_text) && DoesStringMatchConstraint (sfp->except_text, scp))
    {
      sfp->except_text = MemFree (sfp->except_text);
      rval = TRUE;
    }
  }
  /* evidence */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue)))
  {
    if ((sfp->exp_ev == 1 && DoesStringMatchConstraint("experimental", scp))
        || (sfp->exp_ev == 2 && DoesStringMatchConstraint("non-experimental", scp))) {
      sfp->exp_ev = 0;
      rval = TRUE;
    }
  }

  /* citation */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_citation)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("citation", field->field->data.ptrvalue)))
  {
    if (sfp->cit != NULL) {
      sfp->cit = PubSetFree (sfp->cit);
      rval = TRUE;
    }
  }


  /* fields common to some features */
  /* product */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue)))
  {
    if (prp != NULL) {
      rval = RemoveValNodeStringMatch (&(prp->name), scp);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      rval = RemoveRNAProductString (sfp, scp);
    }
  }

  /* Gene fields */
  /* locus */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
       || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && DoesStringMatchConstraint (grp->locus, scp)) {
      grp->locus = MemFree (grp->locus);
      rval = TRUE;
    }
  }
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
       || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->desc) && DoesStringMatchConstraint(grp->desc, scp))
    {
      grp->desc = MemFree (grp->desc);
      rval = TRUE;
    }
  }
  /* maploc */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->maploc) && DoesStringMatchConstraint(grp->maploc, scp))
    {
      grp->maploc = MemFree (grp->maploc);
      rval = TRUE;
    }
  }
  /* allele */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->allele) && DoesStringMatchConstraint(grp->allele, scp))
    {
      grp->allele = MemFree (grp->allele);
      rval = TRUE;
    }
  }
  /* locus_tag */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus_tag) && DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      grp->locus_tag = MemFree (grp->locus_tag);
      rval = TRUE;
    }
  }
  /* synonym */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(grp->syn), scp);
  }
  /* gene comment */
  if (field->field->choice == FeatQualChoice_legal_qual
      && field->field->data.intvalue == Feat_qual_legal_gene_comment
      && gene != NULL
      && !StringHasNoText (gene->comment)
      && DoesStringMatchConstraint (gene->comment, scp)) {
    gene->comment = MemFree (gene->comment);
    rval = TRUE;
  }

  /* protein fields */
  /* note - product handled above */
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (!StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
      prp->desc = MemFree (prp->desc);
      rval = TRUE;
    }
  }
  /* ec_number */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(prp->ec), scp);
  }
  /* activity */
  if (((field->field->choice == FeatQualChoice_legal_qual 
        && (field->field->data.intvalue == Feat_qual_legal_activity 
            || field->field->data.intvalue == Feat_qual_legal_function))
       || (field->field->choice == FeatQualChoice_illegal_qual 
           && (DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)
               || DoesStringMatchConstraint ("function", field->field->data.ptrvalue))))
      && prp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(prp->activity), scp);
  }
  
  /* RNA fields */
  /* anticodon */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_anticodon)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("anticodon", field->field->data.ptrvalue)))
       && rrp != NULL && rrp->ext.choice == 2)
  {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL && trp->anticodon != NULL) {
      trp->anticodon = SeqLocFree (trp->anticodon);
      rval = TRUE;
    }
  }
  /* codons recognized */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_codons_recognized)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("codon-recognized", field->field->data.ptrvalue)))
       && rrp != NULL && rrp->ext.choice == 2)
  {
    rval = RemovetRNACodons_Recognized (sfp);
  }

  if (!rval) {
    /* actual GenBank qualifiers */
    if (field->field->choice == FeatQualChoice_legal_qual) 
    {
      gbqual = GetGBQualFromFeatQual (field->field->data.intvalue, &subfield);
      if (gbqual > -1) {
        rval = RemoveGBQualMatch (&(sfp->qual), ParFlat_GBQual_names [gbqual].name, subfield, scp);
      } else {
        /* need to do something with non-qualifier qualifiers */
      }
    } else {
      rval = RemoveGBQualMatchConstraintName (&(sfp->qual), field->field->data.ptrvalue, scp);
    }
  }

  return rval;
}


static Boolean ChooseBestFrame (SeqFeatPtr sfp)
{
  CdRegionPtr  crp;
  Uint1        new_frame = 0, i, orig_frame;
  ByteStorePtr bs;
  Int4         lens [3];
  Int4         max;
  Boolean      retval = TRUE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;
  
  crp = sfp->data.value.ptrvalue;
  if (crp == NULL) return FALSE;
  orig_frame = crp->frame;

  max = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    lens[i - 1] = BSLen (bs);
    BSFree (bs);
    if (lens[i - 1] > max) {
      max = lens[i - 1];
      new_frame = i;
    }
  }
  for (i = 1; i <= 3; i++) {
    if (lens [i - 1] == max && i != new_frame) {
      retval = FALSE;
    }
  }
  if (retval) {
    crp->frame = new_frame;
  } else {
    crp->frame = orig_frame;
  }
  return retval;
}


static SeqFeatPtr CreateGeneForFeature (SeqFeatPtr sfp)
{
  BioseqPtr  bsp;
  SeqFeatPtr gene = NULL;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) {
    return NULL;
  } else {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, sfp->location);
      if (gene != NULL) {
        gene->data.value.ptrvalue = GeneRefNew();
      }
    }
  }
  return gene;
}


static void AdjustProteinSequenceForReadingFrame (SeqFeatPtr cds);


static Boolean SetCitationTextOnFeature (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr value, Uint2 existing_text, ValNodePtr cit_list)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ValNodePtr  list = NULL, vnp;
  Boolean     rval = FALSE, already_present = FALSE;
  Int4        new_number, serial_number;
  ValNodePtr  min_pub, new_list;

  if (sfp == NULL) {
    return FALSE;
  }

  if (sfp->cit != NULL && existing_text == ExistingTextOption_leave_old) {
    return FALSE;
  }

  if (!IsAllDigits (value)) {
    return FALSE;
  }

  new_number = atoi (value);

  bsp = GetSequenceForObject (OBJ_SEQFEAT, sfp);

  if (cit_list == NULL) {
    /* list not provided - must create now */
    sep = SeqMgrGetSeqEntryForData (bsp);
    list = GetCitListsForSeqEntry (sep);
    cit_list = list;
  } 

  min_pub = GetMinPubForCitationNumber (bsp, new_number, cit_list);
  if (min_pub == NULL) {
    list = PubSerialNumberListFree (list);
    return FALSE;
  }

  if (existing_text == ExistingTextOption_replace_old) {
    sfp->cit = PubSetFree (sfp->cit);
    sfp->cit = ValNodeNew (NULL);
    sfp->cit->choice = 1;
    new_list = NULL;
    ValNodeLink (&new_list, AsnIoMemCopy (min_pub->data.ptrvalue, (AsnReadFunc) PubAsnRead, (AsnWriteFunc) PubAsnWrite));
    sfp->cit->data.ptrvalue = new_list;
    rval = TRUE;
  } else {    
    for (vnp = sfp->cit->data.ptrvalue; vnp != NULL && !already_present; vnp = vnp->next) {    
      serial_number = GetCitationNumberForMinPub (bsp, vnp, cit_list);
      if (serial_number == new_number) {
        already_present = TRUE;
      }
    }
    if (!already_present) {
      new_list = sfp->cit->data.ptrvalue;
      ValNodeLink (&new_list, AsnIoMemCopy (min_pub->data.ptrvalue, (AsnReadFunc) PubAsnRead, (AsnWriteFunc) PubAsnWrite));
      sfp->cit->data.ptrvalue = new_list;
      rval = TRUE;
    }
  }
  
  list = PubSerialNumberListFree (list);

  return rval;
}


static Boolean SetQualOnFeatureEx (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text, BatchExtraPtr batch_extra)
{
  Boolean rval = FALSE;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  CharPtr    tmp;
  CdRegionPtr crp;
  SeqFeatPtr  gene = NULL;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return FALSE;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return FALSE;
  }

  // for gene fields
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
    gene = sfp;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (gene != NULL) {
        grp = gene->data.value.ptrvalue;
      }
    }
  }

  // for protein fields
  prp = GetProtRefForFeature (sfp);

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (DoesStringMatchConstraint(sfp->comment, scp))
    {
      rval = SetStringValue ( &(sfp->comment), value, existing_text);
    }
  }
  /* db-xref */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue)))
  {
    rval = SetDbxrefString (&(sfp->dbxref), scp, value, existing_text);
  }
  /* exception */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue)))
  {
    if (DoesStringMatchConstraint(sfp->except_text, scp))
    {
      rval = SetStringValue ( &(sfp->except_text), value, existing_text);
    }
  }
  /* evidence */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue)))
  {
    tmp = NULL;
    if (sfp->exp_ev == 1)
    {
      tmp = StringSave ("experimental");
    }
    else if (sfp->exp_ev == 2)
    {
      tmp = StringSave ("non-experimental");
    }
    if (DoesStringMatchConstraint(tmp, scp)) {
      rval = SetStringValue (&tmp, value, existing_text);
      if (rval) {
        rval = FALSE;
        if (StringICmp (tmp, "experimental") == 0) {
          sfp->exp_ev = 1;
          rval = TRUE;
        } else if (StringICmp (tmp, "non-experimental") == 0) {
          sfp->exp_ev = 2;
          rval = TRUE;
        } else if (StringHasNoText (tmp)) {
          sfp->exp_ev = 0;
          rval = TRUE;
        }
      }
    }
    tmp = MemFree (tmp);
  }

  /* citation */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_citation)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("citation", field->field->data.ptrvalue)))
  {
    rval = SetCitationTextOnFeature (sfp, scp, value, existing_text, batch_extra == NULL ? NULL : batch_extra->cit_list);
  }

  /* fields common to some features */
  /* product */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue)))
  {
    if (prp != NULL) {
      rval = SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      rval = SetRNAProductString (sfp, scp, value, existing_text);
    }
  }

  /* Gene fields */
  /* locus */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
  {
    if (grp == NULL && IsStringConstraintEmpty (scp))
    {
      /* create new gene feature */
      gene = CreateGeneForFeature (sfp);
      if (gene != NULL)
      {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
      }
    }     
    if (grp != NULL && DoesStringMatchConstraint(grp->locus, scp))
    {
      rval = SetStringValue (&(grp->locus), value, existing_text);
    }
  }

  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->desc, scp))
    {
      rval = SetStringValue (&(grp->desc), value, existing_text);
    }
  }
  /* maploc */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->maploc, scp))
    {
      rval = SetStringValue (&(grp->maploc), value, existing_text);
    }
  }
  /* allele */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->allele, scp))
    {
      rval = SetStringValue (&(grp->allele), value, existing_text);
    }
  }
  /* locus_tag */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      rval = SetStringValue (&(grp->locus_tag), value, existing_text);
    }
  }
  /* synonym */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(grp->syn), scp, value, existing_text);
  }
  /* gene comment */
  if (field->field->choice == FeatQualChoice_legal_qual
      && field->field->data.intvalue == Feat_qual_legal_gene_comment
      && gene != NULL) {
    rval = SetStringValue (&(gene->comment), value, existing_text);
  }

  /* protein fields */
  /* note - product handled above */
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (DoesStringMatchConstraint(prp->desc, scp)) {
      rval = SetStringValue (&(prp->desc), value, existing_text);
    }
  }
  /* ec_number */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
  }
  /* activity */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_activity)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
  }
 
  /* codon start */
  /* note - if product existed before, it will be retranslated */
  if (field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_codon_start
      && sfp->data.choice == SEQFEAT_CDREGION) 
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (StringICmp (value, "best") == 0)
    {
      rval = ChooseBestFrame (sfp);
    }
    else if (StringCmp (value, "1") == 0) 
    {
      crp->frame = 1;
      rval = TRUE;
    }
    else if (StringCmp (value, "2") == 0) 
    {
      crp->frame = 2;
      rval = TRUE;
    }
    else if (StringCmp (value, "3") == 0)
    {
      crp->frame = 3;
      rval = TRUE;
    }
    if (rval && sfp->product != NULL) {
      AdjustProteinSequenceForReadingFrame (sfp);
    }
  } 

  /* tRNA fields */
  if (sfp->idx.subtype == FEATDEF_tRNA
      && ((field->field->choice == FeatQualChoice_legal_qual
           && field->field->data.intvalue == Feat_qual_legal_codons_recognized)
          || (field->field->choice == FeatQualChoice_illegal_qual
           && DoesStringMatchConstraint ("codon-recognized", field->field->data.ptrvalue))))
  {
    rval = SettRNACodons_Recognized (sfp, scp, value, existing_text);
  }

  if (sfp->idx.subtype == FEATDEF_tRNA
      && ((field->field->choice == FeatQualChoice_legal_qual
           && field->field->data.intvalue == Feat_qual_legal_anticodon)
          || (field->field->choice == FeatQualChoice_illegal_qual
           && DoesStringMatchConstraint ("anticodon", field->field->data.ptrvalue))))
  {
    rval = SetAnticodon (sfp, scp, value, existing_text);
  }


  /* actual GenBank qualifiers */
  if (!rval)
  {
    rval = SetStringInGBQualList (&(sfp->qual), field->field, scp, value, existing_text);
  }
  return rval;
}


static Boolean SetQualOnFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  return SetQualOnFeatureEx (sfp, field, scp, value, existing_text, NULL);
}

static void AddLegalFeatureField (ValNodePtr PNTR list, Uint2 featdef, Uint2 qual)
{
  FeatureFieldPtr ffield;
  Int4            gbqual, num_subfields, i, legal_qual;

  if (list == NULL) return;

  ffield = FeatureFieldNew ();
  ffield->type = GetFeatureTypeFromFeatdef (featdef);
  ValNodeAddInt (&(ffield->field), FeatQualChoice_legal_qual, qual);
  ValNodeAddPointer (list, FieldType_feature_field, ffield);

  /* also add subfields */
  gbqual = GetGBQualFromFeatQual (qual, NULL);
  num_subfields = NumGbQualSubfields (gbqual);
  for (i = 1; i <= num_subfields; i++) {
    legal_qual = GetFeatQualByGBQualAndSubfield (gbqual, i);
    if (legal_qual > -1) {
      ffield = FeatureFieldNew ();
      ffield->type = GetFeatureTypeFromFeatdef (featdef);
      ValNodeAddInt (&(ffield->field), FeatQualChoice_legal_qual, legal_qual);
      ValNodeAddPointer (list, FieldType_feature_field, ffield);
    }
  }

}


static ValNodePtr GetFieldListFromFeature (SeqFeatPtr sfp)
{
  GeneRefPtr grp = NULL;
  SeqFeatPtr gene = NULL;
  ProtRefPtr prp = NULL;
  ValNodePtr list = NULL;
  GBQualPtr  qual;
  Int4       qual_num;

  if (sfp == NULL)
  {
    return NULL;
  }

  // for gene fields
  GetGeneInfoForFeature (sfp, &grp, &gene);

  /* add gene-specific fields */
  if (grp != NULL) {
    if (!StringHasNoText (grp->locus)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_gene);
    }
    if (!StringHasNoText (grp->allele)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_allele);
    }
    if (!StringHasNoText (grp->desc)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_gene_description);
    }
    if (!StringHasNoText (grp->maploc)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_map);
    }
    if (!StringHasNoText (grp->locus_tag)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_locus_tag);
    }
    if (grp->syn != NULL) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_synonym);
    }
  }

  /* add protein-specific fields */
  prp = GetProtRefForFeature (sfp);
  if (prp != NULL) {
    /* product name */
    if (prp->name != NULL) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_product);
    }
    /* protein description */
    if (!StringHasNoText (prp->desc)) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_description);
    }
    /* ec_number */
    if (prp->ec != NULL) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_ec_number);
    }
    /* activity */
    if (prp->activity != NULL) {
      AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_activity);
    }
  }

  /* fields common to all features */
  /* note, also known as comment */
  if (!StringHasNoText (sfp->comment)) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_note);
  }
  /* db-xref */
  if (sfp->dbxref != NULL) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_db_xref);
  }
  /* exception */
  if (!StringHasNoText (sfp->except_text)) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_exception);
  }
  /* evidence */
  if (sfp->exp_ev > 0) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_evidence);
  }

  /* citation */
  if (sfp->cit != NULL) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_citation);
  }

  /* RNA specific */
  if (sfp->data.choice == SEQFEAT_RNA) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_product);
  }

  /* coding regions */
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    AddLegalFeatureField (&list, sfp->idx.subtype, Feat_qual_legal_codon_start);
  }    

  /* actual GenBank qualifiers */
  for (qual = sfp->qual; qual != NULL; qual = qual->next) 
  {
    qual_num = GetFeatQualByName (qual->qual);
    if (qual_num > -1) {
      AddLegalFeatureField (&list, sfp->idx.subtype, qual_num);
    }
  }
  return list;
}


NLM_EXTERN CharPtr GetSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint)
{
  CharPtr str = NULL;
  SubSourcePtr ssp;
  OrgModPtr mod;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1;
  Int4 subfield;
  ValNode vn;
  Char buf[15];

  if (biop == NULL || scp == NULL) return NULL;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if (biop->org != NULL && !StringHasNoText (biop->org->taxname)
            && DoesStringMatchConstraint (biop->org->taxname, constraint)) {
          str = StringSave (biop->org->taxname);
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if (biop->org != NULL && !StringHasNoText (biop->org->common)
            && DoesStringMatchConstraint (biop->org->common, constraint)) {
          str = StringSave (biop->org->common);
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if (biop->org != NULL && biop->org->orgname != NULL && !StringHasNoText (biop->org->orgname->lineage)
            && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint)) {
          str = StringSave (biop->org->orgname->lineage);
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if (biop->org != NULL && biop->org->orgname != NULL  && !StringHasNoText (biop->org->orgname->div)
            && DoesStringMatchConstraint (biop->org->orgname->div, constraint)) {
          str = StringSave (biop->org->orgname->div);
        }
      } else if (scp->data.intvalue == Source_qual_dbxref) {
        if (biop->org != NULL) {
          str = GetDbxrefString (biop->org->db, constraint);
        }
      } else if (scp->data.intvalue == Source_qual_all_notes) {
        vn.choice = SourceQualChoice_textqual;
        vn.data.intvalue = Source_qual_subsource_note;
        vn.next = NULL;
        str = GetSourceQualFromBioSource (biop, &vn, constraint);
        if (str == NULL) {
          vn.data.intvalue = Source_qual_orgmod_note;
          str = GetSourceQualFromBioSource (biop, &vn, constraint);
        }
      } else if (scp->data.intvalue == Source_qual_all_quals) {
        /* will not do */
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue, &subfield);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue, &subfield);
          for (ssp = biop->subtype; ssp != NULL && str == NULL; ssp = ssp->next) {
            if (ssp->subtype == subsrc_subtype) {
              if (StringHasNoText (ssp->name)) {
                if (IsNonTextSourceQual (scp->data.intvalue)
                    && DoesStringMatchConstraint ("TRUE", constraint)) {
                  str = StringSave ("TRUE");
                }
              } else {
                if (subfield == 0) {
                  if (DoesStringMatchConstraint (ssp->name, constraint)) {
                    str = StringSave (ssp->name);
                  }
                } else {
                  str = GetThreeFieldSubfield (ssp->name, subfield);
                  if (StringHasNoText (str) || !DoesStringMatchConstraint (str, constraint)) {
                    str = MemFree (str);
                  }
                }
              }
            }
          }
        } else {
          if (biop->org != NULL && biop->org->orgname != NULL) {
            for (mod = biop->org->orgname->mod; mod != NULL && str == NULL; mod = mod->next) {
              if (mod->subtype == orgmod_subtype) {
                if (StringHasNoText (mod->subname)) {
                  if (IsNonTextSourceQual (scp->data.intvalue)
                      && DoesStringMatchConstraint ("TRUE", constraint)) {
                    str = StringSave ("TRUE");
                  }
                } else {
                  if (subfield == 0) {
                    if (DoesStringMatchConstraint (mod->subname, constraint)) {
                      str = StringSave (mod->subname);
                    }
                  } else {
                    str = GetThreeFieldSubfield (mod->subname, subfield);
                    if (StringHasNoText (str) || !DoesStringMatchConstraint (str, constraint)) {
                      str = MemFree (str);
                    }
                  }
                }
              }
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        str = StringSave (str);
      } else {
        str = NULL;
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        str = StringSave (str);
      } else {
        str = NULL;
      }
      break;
    case SourceQualChoice_gcode:
      if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->gcode != 0) {
        sprintf (buf, "%d", biop->org->orgname->gcode);
        str = StringSave (buf);
      }
      break;
    case SourceQualChoice_mgcode:
      if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->mgcode != 0) {
        sprintf (buf, "%d", biop->org->orgname->mgcode);
        str = StringSave (buf);
      }
      break;
  }
  return str;
}


static Boolean RemoveAllSourceQualsFromBioSource (BioSourcePtr biop, StringConstraintPtr constraint)
{
  Int4 i;
  Boolean rval = FALSE;
  ValNode vn;

  vn.next = NULL;
  vn.choice = SourceQualChoice_textqual;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual_scqual[i].srcqual != Source_qual_all_quals
        && srcqual_scqual[i].srcqual != Source_qual_all_notes) {
      vn.data.intvalue = srcqual_scqual[i].srcqual;
      rval |= RemoveSourceQualFromBioSource (biop, &vn, constraint);
    }
  }
  return rval;
}

static void Lcl_RemoveOldName (OrgRefPtr orp)
{
  OrgModPtr prev = NULL, curr, next_mod;
  
  if (orp == NULL || orp->orgname == NULL) return;
  
  curr = orp->orgname->mod;
  while (curr != NULL)
  {
    next_mod = curr->next;
    if (curr->subtype == ORGMOD_old_name)
    {
      if (prev == NULL)
      {
        orp->orgname->mod = curr->next;
      }
      else
      {
        prev->next = curr->next;
      }
      curr->next = NULL;
      OrgModFree (curr);
    }
    else
    {
      prev = curr;
    }
    
    curr = next_mod;
  }
}

NLM_EXTERN Boolean RemoveSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint)
{
  SubSourcePtr ssp, ssp_prev = NULL, ssp_next;
  OrgModPtr mod, mod_prev = NULL, mod_next;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1, subfield;
  CharPtr str, tmp;
  Boolean rval = FALSE, do_remove, does_match;
  ValNode vn;

  if (biop == NULL || scp == NULL) return FALSE;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if (biop->org != NULL && !StringHasNoText (biop->org->taxname)
            && DoesStringMatchConstraint (biop->org->taxname, constraint)) {
          biop->org->taxname = MemFree (biop->org->taxname);
          RemoveTaxRef (biop->org);
          Lcl_RemoveOldName (biop->org);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if (biop->org != NULL && !StringHasNoText (biop->org->common)
            && DoesStringMatchConstraint (biop->org->common, constraint)) {
          biop->org->common = MemFree (biop->org->common);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if (biop->org != NULL && biop->org->orgname != NULL && !StringHasNoText (biop->org->orgname->lineage)
            && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint)) {
          biop->org->orgname->lineage = MemFree (biop->org->orgname->lineage);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if (biop->org != NULL && biop->org->orgname != NULL  && !StringHasNoText (biop->org->orgname->div)
            && DoesStringMatchConstraint (biop->org->orgname->div, constraint)) {
          biop->org->orgname->div = MemFree (biop->org->orgname->div);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_dbxref) {
        if (biop->org != NULL) {
          rval = RemoveDbxrefString (&(biop->org->db), constraint);
        }        
      } else if (scp->data.intvalue == Source_qual_all_notes) {
        vn.choice = SourceQualChoice_textqual;
        vn.data.intvalue = Source_qual_subsource_note;
        vn.next = NULL;
        rval |= RemoveSourceQualFromBioSource (biop, &vn, constraint);
        vn.data.intvalue = Source_qual_orgmod_note;
        rval |= RemoveSourceQualFromBioSource (biop, &vn, constraint);
      } else if (scp->data.intvalue == Source_qual_all_quals) {
        rval |= RemoveAllSourceQualsFromBioSource (biop, constraint);
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue, &subfield);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue, &subfield);
          ssp = biop->subtype;
          while (ssp != NULL) {
            ssp_next = ssp->next;
            do_remove = FALSE;
            if (ssp->subtype == subsrc_subtype) {
              if (subfield == 0) {
                if (DoesStringMatchConstraint (ssp->name, constraint)) {
                  do_remove = TRUE;
                }
              } else {
                does_match = TRUE;
                if (!IsStringConstraintEmpty (constraint)) {
                  tmp = GetThreeFieldSubfield (ssp->name, subfield);
                  does_match = DoesStringMatchConstraint (tmp, constraint);
                  tmp = MemFree (tmp);
                }
                if (does_match) {
                  rval |= RemoveThreeFieldSubfield (&(ssp->name), subfield);
                  if (StringHasNoText (ssp->name)) {
                    do_remove = TRUE;
                  }
                }
              }
            }
            if (do_remove) {
              if (ssp_prev == NULL) {
                biop->subtype = ssp->next;
              } else {
                ssp_prev->next = ssp->next;
              }
              ssp->next = NULL;
              ssp = SubSourceFree (ssp);
              rval = TRUE;
            } else {
              ssp_prev = ssp;
            }
            ssp = ssp_next;
          }
        } else {
          if (biop->org != NULL && biop->org->orgname != NULL) {
            mod = biop->org->orgname->mod;
            while (mod != NULL) {
              mod_next = mod->next;
              do_remove = FALSE;
              if (mod->subtype == orgmod_subtype) {
                if (subfield == 0) {
                  if (DoesStringMatchConstraint (mod->subname, constraint)) {
                    do_remove = TRUE;
                  }
                } else {
                  does_match = TRUE;
                  if (!IsStringConstraintEmpty (constraint)) {
                    tmp = GetThreeFieldSubfield (mod->subname, subfield);
                    does_match = DoesStringMatchConstraint (tmp, constraint);
                    tmp = MemFree (tmp);
                  }
                  if (does_match) {
                    rval |= RemoveThreeFieldSubfield (&(mod->subname), subfield);
                  }
                  if (StringHasNoText (mod->subname)) {
                    do_remove = TRUE;
                  }
                }
              }
              if (do_remove) {
                if (mod_prev == NULL) {
                  biop->org->orgname->mod = mod->next;
                } else {
                  mod_prev->next = mod->next;
                }
                mod->next = NULL;
                mod = OrgModFree (mod);
                rval = TRUE;
              } else {
                mod_prev = mod;
              }
              mod = mod_next;
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        if (scp->data.intvalue == 0 || biop->genome == GenomeFromSrcLoc (scp->data.intvalue)) {
          biop->genome = 0;
          rval = TRUE;
        }
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        if (scp->data.intvalue == 0 || biop->origin == OriginFromSrcOrig (scp->data.intvalue)) {
          biop->origin = 0;
          rval = TRUE;
        }
      }
      break; 
    case SourceQualChoice_gcode:
      if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->gcode != 0) {
        biop->org->orgname->gcode = 0;
        rval = TRUE;
      }
      break;
    case SourceQualChoice_mgcode:
      if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->mgcode != 0) {
        biop->org->orgname->mgcode = 0;
        rval = TRUE;
      }
      break;
  }
  return rval;
}


NLM_EXTERN Boolean SetSourceQualInBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint, CharPtr value, Uint2 existing_text)
{
  SubSourcePtr ssp, ssp_prev = NULL, ssp_next;
  OrgModPtr mod, mod_prev = NULL, mod_next;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1, subfield;
  CharPtr str, tmp;
  Boolean rval = FALSE, found = FALSE, does_match;
  ValNode vn;

  if (biop == NULL || scp == NULL) return FALSE;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if ((biop->org == NULL && IsStringConstraintEmpty (constraint))
            || (biop->org != NULL
                && DoesStringMatchConstraint (biop->org->taxname, constraint))) {
          if (biop->org == NULL) {
            biop->org = OrgRefNew();
          }
          rval = SetStringValue (&(biop->org->taxname), value, existing_text);
          if (rval) {
            RemoveTaxRef (biop->org);
            Lcl_RemoveOldName (biop->org);
          }
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if ((biop->org == NULL && IsStringConstraintEmpty (constraint)) 
            || (biop->org != NULL
                && DoesStringMatchConstraint (biop->org->common, constraint))) {
          if (biop->org == NULL) {
            biop->org = OrgRefNew();
          }
          rval = SetStringValue (&(biop->org->common), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if ((biop->org == NULL && IsStringConstraintEmpty (constraint)) 
            ||(biop->org != NULL && biop->org->orgname == NULL && IsStringConstraintEmpty (constraint))
            ||(biop->org != NULL && biop->org->orgname != NULL 
               && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint))) {
          if (biop->org == NULL) {
            biop->org = OrgRefNew();
          }
          if (biop->org->orgname == NULL) {
            biop->org->orgname = OrgNameNew ();
          }
          rval = SetStringValue (&(biop->org->orgname->lineage), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if ((biop->org == NULL && IsStringConstraintEmpty (constraint)) 
            || (biop->org != NULL && biop->org->orgname == NULL && IsStringConstraintEmpty (constraint))
            || (biop->org != NULL && biop->org->orgname != NULL
                && DoesStringMatchConstraint (biop->org->orgname->div, constraint))) {
          if (biop->org == NULL) {
            biop->org = OrgRefNew();
          }
          if (biop->org->orgname == NULL) {
            biop->org->orgname = OrgNameNew ();
          }
          rval = SetStringValue (&(biop->org->orgname->div), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_dbxref) {
        if (biop->org == NULL) {
          biop->org = OrgRefNew ();
        }
        rval = SetDbxrefString (&(biop->org->db), constraint, value, existing_text);
      } else if (scp->data.intvalue == Source_qual_all_notes) {
        vn.choice = SourceQualChoice_textqual;
        vn.data.intvalue = Source_qual_subsource_note;
        vn.next = NULL;
        rval |= SetSourceQualInBioSource (biop, &vn, constraint, value, existing_text);
        vn.data.intvalue = Source_qual_orgmod_note;
        rval |= SetSourceQualInBioSource (biop, &vn, constraint, value, existing_text);
      } else if (scp->data.intvalue == Source_qual_all_quals) {
        /* will not do this */
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue, &subfield);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue, &subfield);
          if (subsrc_subtype > -1) {
            if (existing_text == ExistingTextOption_add_qual) {
              /* create new subsource */
              ssp = SubSourceNew ();
              ssp->subtype = subsrc_subtype;
              rval = SetThreeFieldSubfield (&(ssp->name), subfield, value, existing_text);
              /* find last in current list */
              ssp_prev = biop->subtype;
              while (ssp_prev != NULL && ssp_prev->next != NULL) {
                ssp_prev = ssp_prev->next;
              }

              /* add to end of list */
              if (ssp_prev == NULL) {
                biop->subtype = ssp;
              } else {
                ssp_prev->next = ssp;
              }
            } else {              
              ssp = biop->subtype;
              while (ssp != NULL) {
                ssp_next = ssp->next;
                if (ssp->subtype == subsrc_subtype) {
                  if (subfield == 0) {
                    if (DoesStringMatchConstraint (ssp->name, constraint)) {
                      rval = SetStringValue (&(ssp->name), value, existing_text);
                      found = TRUE;
                    }
                  } else {
                    does_match = TRUE;
                    if (!IsStringConstraintEmpty (constraint)) {
                      tmp = GetThreeFieldSubfield (ssp->name, subfield);
                      does_match = DoesStringMatchConstraint (tmp, constraint);
                    }
                    if (does_match) {
                      rval = SetThreeFieldSubfield (&(ssp->name), subfield, value, existing_text);
                      found = TRUE;
                    }
                  }
                  if (rval && StringHasNoText (ssp->name) && !IsNonTextSourceQual(scp->data.intvalue)) {
                    if (ssp_prev == NULL) {
                      biop->subtype = ssp->next;
                    } else {
                      ssp_prev->next = ssp->next;
                    }
                    ssp->next = NULL;
                    ssp = SubSourceFree (ssp);
                  } else {
                    ssp_prev = ssp;
                  }
                } else {
                  ssp_prev = ssp;
                }
                ssp = ssp_next;
              }
              if (!found && IsStringConstraintEmpty (constraint)) {
                ssp = SubSourceNew ();
                ssp->subtype = subsrc_subtype;
                if (StringHasNoText (value) && IsNonTextSourceQual(scp->data.intvalue)) {
                  ssp->name = StringSave ("");
                } else {
                  rval = SetThreeFieldSubfield (&(ssp->name), subfield, value, existing_text);
                }
                if (ssp_prev == NULL) {
                  biop->subtype = ssp;
                } else {
                  ssp_prev->next = ssp;
                }
              }
            }
          }
        } else {
          if (existing_text == ExistingTextOption_add_qual) {
            if (biop->org == NULL) {
              biop->org = OrgRefNew();
            }
            if (biop->org->orgname == NULL) {
              biop->org->orgname = OrgNameNew();
            }
            /* create new orgmod */
            mod = OrgModNew ();
            mod->subtype = orgmod_subtype;
            rval = SetThreeFieldSubfield (&(mod->subname), subfield, value, existing_text);
            /* find last in current list */
            mod_prev = biop->org->orgname->mod;
            while (mod_prev != NULL && mod_prev->next != NULL) {
              mod_prev = mod_prev->next;
            }
            /* add to end of list */
            if (mod_prev == NULL) {
              biop->org->orgname->mod = mod;
            } else {
              mod_prev->next = mod;
            }
          } else {        
            if (biop->org != NULL && biop->org->orgname != NULL) {
              mod = biop->org->orgname->mod;
              while (mod != NULL) {
                mod_next = mod->next;
                if (mod->subtype == orgmod_subtype) {
                  if (subfield == 0) {
                    if (DoesStringMatchConstraint (mod->subname, constraint)) {
                      rval = SetStringValue (&(mod->subname), value, existing_text);
                      found = TRUE;
                    }
                  } else {
                    does_match = TRUE;
                    if (!IsStringConstraintEmpty (constraint)) {
                      tmp = GetThreeFieldSubfield (mod->subname, subfield);
                      does_match = DoesStringMatchConstraint (tmp, constraint);
                      tmp = MemFree (tmp);
                    }
                    if (does_match) {
                      rval = SetThreeFieldSubfield (&(mod->subname), subfield, value, existing_text);
                      found = TRUE;
                    }
                  }
                  if (rval && StringHasNoText (mod->subname) && !IsNonTextSourceQual(scp->data.intvalue)) {
                    if (mod_prev == NULL) {
                      biop->org->orgname->mod = mod->next;
                    } else {
                      mod_prev->next = mod->next;
                    }
                    mod->next = NULL;
                    mod = OrgModFree (mod);
                  } else {
                    mod_prev = mod;
                  }
                } else {
                  mod_prev = mod;
                }
                mod = mod_next;
              }
            }
            if (!found && IsStringConstraintEmpty (constraint)) {
              if (biop->org == NULL) {
                biop->org = OrgRefNew();
              }
              if (biop->org->orgname == NULL) {
                biop->org->orgname = OrgNameNew();
              }
              mod = OrgModNew ();
              mod->subtype = orgmod_subtype;
              rval = SetThreeFieldSubfield (&(mod->subname), subfield, value, existing_text);
              if (mod_prev == NULL) {
                biop->org->orgname->mod = mod;
              } else {
                mod_prev->next = mod;
              }
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        biop->genome = GenomeFromSrcLoc (scp->data.intvalue);
        rval = TRUE;
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        biop->origin = OriginFromSrcOrig(scp->data.intvalue);
        rval = TRUE;
      }
      break; 
    case SourceQualChoice_gcode:
      if (biop->org == NULL) {
        biop->org = OrgRefNew();
      }
      if (biop->org->orgname == NULL) {
        biop->org->orgname = OrgNameNew();
      }
      biop->org->orgname->gcode = scp->data.intvalue;
      rval = TRUE;
      break;
    case SourceQualChoice_mgcode:
      if (biop->org == NULL) {
        biop->org = OrgRefNew();
      }
      if (biop->org->orgname == NULL) {
        biop->org->orgname = OrgNameNew();
      }
      biop->org->orgname->mgcode = scp->data.intvalue;
      rval = TRUE;
      break;
  }
  return rval;
}


NLM_EXTERN BioseqPtr GetRepresentativeBioseqFromBioseqSet (BioseqSetPtr bssp)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp = NULL;

  if (bssp == NULL || (bssp->_class != BioseqseqSet_class_segset && bssp->_class != BioseqseqSet_class_nuc_prot)) {
    return NULL;
  }
  sep = bssp->seq_set;
  if (sep->data.ptrvalue == NULL) {
    bsp = NULL;
  } else if (IS_Bioseq(sep)) {
    bsp = sep->data.ptrvalue;
  } else if (IS_Bioseq_set (sep)) {
    bsp = GetRepresentativeBioseqFromBioseqSet (sep->data.ptrvalue);
  }
  return bsp;
}


NLM_EXTERN BioseqPtr GetSequenceForObject (Uint1 choice, Pointer data)
{
  BioseqPtr bsp = NULL;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  CGPSetPtr cgp;
  ValNodePtr vnp;

  if (data == NULL) return NULL;

  switch (choice) {
    case OBJ_BIOSEQ:
      bsp = (BioseqPtr) data;
      break;
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) data;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      break;
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) data;
      if (sdp->extended) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ && ovp->idx.parentptr != NULL) {
          bsp = ovp->idx.parentptr;
        } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          bsp = GetRepresentativeBioseqFromBioseqSet (ovp->idx.parentptr);
        }
      }
      break;
    case 0:
      cgp = (CGPSetPtr) data;
      for (vnp = cgp->cds_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      for (vnp = cgp->mrna_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      break;
      for (vnp = cgp->gene_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      break;
  }
  return bsp;
}


NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data)
{
  BioSourcePtr biop = NULL;
  SeqDescrPtr  sdp;
  SeqFeatPtr   sfp;
  BioseqPtr    bsp = NULL;
  SeqMgrDescContext context;

  if (data == NULL) return NULL;

  switch (choice)
  {
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) data;
      if (sdp->choice == Seq_descr_source) {
        biop = sdp->data.ptrvalue;
      }
      break;
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) data;
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = sfp->data.value.ptrvalue;
      }
      break;
  }
  if (biop == NULL) {
    bsp = GetSequenceForObject (choice, data);
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL && sdp->choice == Seq_descr_source) {
      biop = sdp->data.ptrvalue;
    }
  }
  return biop;
}


/* functions for dealing with CDS-Gene-Prot sets */
static CharPtr GetFieldValueFromCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  RnaRefPtr  rrp;
  ProtRefPtr prp;
  FeatureFieldPtr ffield;
  
  if (c == NULL) return NULL;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
    case CDSGeneProt_field_cds_inference:
    case CDSGeneProt_field_codon_start:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->cds_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        str = GetQualFromFeature (sfp, ffield, scp);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_locus:
    case CDSGeneProt_field_gene_inference:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        str = GetQualFromFeature (sfp, ffield, scp);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_description:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->desc) 
            && DoesStringMatchConstraint(grp->desc, scp))
        {
          str = StringSave (grp->desc);
        }
      }
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->allele) 
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          str = StringSave (grp->allele);
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->maploc) 
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          str = StringSave (grp->maploc);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus_tag) 
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          str = StringSave (grp->locus_tag);
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (grp->syn, scp);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          str = GetFirstGBQualMatch (sfp->qual, "old-locus-tag", 0, scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA
            && (rrp = sfp->data.value.ptrvalue) != NULL
            && rrp->ext.choice == 1
            && !StringHasNoText (rrp->ext.value.ptrvalue) 
            && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp))
        {
          str = StringSave (rrp->ext.value.ptrvalue);
        }
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->name, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          str = StringSave (prp->desc);
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->ec, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->activity, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->name, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          str = StringSave (prp->desc);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->ec, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->activity, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
  }
  return str;
}


static Boolean RemoveFieldValueFromCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  RnaRefPtr  rrp;
  ProtRefPtr prp;
  FeatureFieldPtr ffield;
  
  if (c == NULL) return FALSE;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
    case CDSGeneProt_field_cds_inference:
    case CDSGeneProt_field_codon_start:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->cds_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= RemoveQualFromFeature (sfp, ffield, scp);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_locus:
    case CDSGeneProt_field_gene_inference:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= RemoveQualFromFeature (sfp, ffield, scp);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_description:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->desc) 
            && DoesStringMatchConstraint(grp->desc, scp))
        {
          grp->desc = MemFree(grp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->allele) 
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          grp->allele = MemFree (grp->allele);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->maploc) 
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          grp->maploc = MemFree (grp->maploc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus_tag) 
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          grp->locus_tag = MemFree (grp->locus_tag);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(grp->syn), scp);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          rval |= RemoveGBQualMatch (&(sfp->qual), "old-locus-tag", 0, scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA
            && (rrp = sfp->data.value.ptrvalue) != NULL
            && rrp->ext.choice == 1
            && !StringHasNoText (rrp->ext.value.ptrvalue) 
            && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp))
        {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.choice = 0;
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->name), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          prp->desc = MemFree (prp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->ec), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->activity), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->name), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          prp->desc = MemFree (prp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->ec), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->activity), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
  }
  return rval;
}


static SeqFeatPtr CreateGeneForCGPSet (CGPSetPtr c)
{
  SeqFeatPtr gene = NULL, sfp = NULL;
  BioseqPtr  bsp;
  ValNodePtr vnp;

  if (c == NULL) return NULL;

  for (vnp = c->cds_list; vnp != NULL && sfp == NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
  }
  for (vnp = c->mrna_list; vnp != NULL && sfp == NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
  }
  if (sfp != NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, sfp->location);
      if (gene != NULL) {
        gene->data.value.ptrvalue = GeneRefNew();
      }
    }
  }
  return gene;
}


static Boolean SetFieldValueInCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  ProtRefPtr prp;
  FeatureFieldPtr ffield;
  
  if (c == NULL) return FALSE;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
    case CDSGeneProt_field_cds_inference:
    case CDSGeneProt_field_codon_start:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->cds_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= SetQualOnFeature (sfp, ffield, scp, value, existing_text);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_locus:
      if (c->gene_list == NULL && scp == NULL) {
        sfp = CreateGeneForCGPSet (c);
        if (sfp != NULL) {
          ValNodeAddPointer (&(c->gene_list), 0, sfp);
        }
      }
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->locus, scp))
        {
          rval |= SetStringValue ( &(grp->locus), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_description:
    case CDSGeneProt_field_gene_inference:
      ffield = FeatureFieldFromCDSGeneProtField (field);
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= SetQualOnFeature (sfp, ffield, scp, value, existing_text);
      }
      ffield = FeatureFieldFree (ffield);
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          rval |= SetStringValue (&(grp->allele), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          rval |= SetStringValue ( &(grp->maploc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          rval |= SetStringValue ( &(grp->locus_tag), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(grp->syn), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          rval |= RemoveGBQualMatch (&(sfp->qual), "old-locus-tag", 0, scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= SetRNAProductString (sfp, scp, value, existing_text);
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL&& DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(prp->desc, scp)) {
          rval |= SetStringValue ( &(prp->desc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(prp->desc, scp)) {
          rval |= SetStringValue ( &(prp->desc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
  }
  return rval;
}


static MolInfoPtr GetMolInfoForBioseq (BioseqPtr bsp)
{
  MolInfoPtr m = NULL;
  SeqDescrPtr sdp;

  if (bsp == NULL) return NULL;
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
    sdp = sdp->next;
  }
  if (sdp != NULL) {
    m = (MolInfoPtr) sdp->data.ptrvalue;
  }
  return m;
}
  

static CharPtr GetSequenceQualFromBioseq (BioseqPtr bsp, ValNodePtr field)
{
  CharPtr rval = NULL;
  MolInfoPtr m;

  if (bsp == NULL || field == NULL) return NULL;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = BiomolNameFromBiomol (m->biomol);
      }
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = TechNameFromTech (m->tech);
      }
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = CompletenessNameFromCompleteness (m->completeness);
      }
      break;
    case MolinfoField_mol_class:
      rval = MolNameFromMol (bsp->mol);
      break;
    case MolinfoField_topology:
      rval = TopologyNameFromTopology (bsp->topology);
      break;
    case MolinfoField_strand:
      rval = StrandNameFromStrand (bsp->strand);
      break;
  }
  if (rval != NULL) rval = StringSave (rval);
  return rval;
}


static Boolean RemoveSequenceQualFromBioseq (BioseqPtr bsp, ValNodePtr field)
{
  MolInfoPtr m;
  Boolean    rval = FALSE;

  if (bsp == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->biomol = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->tech = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->completeness = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_mol_class:
      bsp->mol = 0;
      rval = TRUE;
      break;
    case MolinfoField_topology:
      bsp->topology = 0;
      rval = TRUE;
      break;
    case MolinfoField_strand:
      bsp->strand = 0;
      rval = TRUE;
      break;
  }
  return rval;
}


static MolInfoPtr AddMolInfoToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  MolInfoPtr  m;

  sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_molinfo);
  m = MolInfoNew ();
  sdp->data.ptrvalue = m;
  return m;
}


static Boolean SetSequenceQualOnBioseq (BioseqPtr bsp, ValNodePtr field)
{
  MolInfoPtr m;
  Boolean    rval = FALSE;

  if (bsp == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->biomol = BiomolFromMoleculeType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->tech = TechFromTechniqueType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->completeness = CompletenessFromCompletednessType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_mol_class:
      bsp->mol = MolFromMoleculeClassType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_topology:
      bsp->topology = TopologyFromTopologyType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_strand:
      bsp->strand = StrandFromStrandType (field->data.intvalue);
      rval = TRUE;
      break;
  }
  return rval;
}


static CharPtr GetGenomeProjectIdFromBioseq (BioseqPtr bsp, StringConstraintPtr scp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  Char              buf[50];
  UserObjectPtr     uop;
  UserFieldPtr      ufp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "GenomeProjectsDB") == 0)
    {
      ufp = uop->data;
      while (ufp != NULL) {
        if (ufp->label != NULL 
            && StringCmp (ufp->label->str, "ProjectID") == 0
            && ufp->choice == 2) {
          sprintf (buf, "%d", ufp->data.intvalue);
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (buf, scp)) {
            return StringSave (buf);
          }
        }
        ufp = ufp->next;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }

  return NULL;
}


static Boolean RemoveGenomeProjectIdFromBioseq (BioseqPtr bsp, StringConstraintPtr scp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  Char              buf[50];
  UserObjectPtr     uop;
  UserFieldPtr      ufp;
  ObjValNodePtr     ovn;
  Boolean           rval = FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "GenomeProjectsDB") == 0)
    {
      ufp = uop->data;
      while (ufp != NULL) {
        if (ufp->label != NULL 
            && StringCmp (ufp->label->str, "ProjectID") == 0
            && ufp->choice == 2) {
          sprintf (buf, "%d", ufp->data.intvalue);
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (buf, scp)) {
            if (sdp->extended != 0) {
              ovn = (ObjValNodePtr) sdp;
              ovn->idx.deleteme = TRUE;
              rval = TRUE;
            }
          }
        }
        ufp = ufp->next;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }
  return rval;
}


static Boolean SetGenomeProjectIdOnBioseq (BioseqPtr bsp, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  Char              buf[50];
  CharPtr           tmp;
  UserObjectPtr     uop;
  UserFieldPtr      ufp;
  Boolean           rval = FALSE;

  if (bsp == NULL || !IsAllDigits (value)) {
    return FALSE;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "GenomeProjectsDB") == 0)
    {
      ufp = uop->data;
      while (ufp != NULL) {
        if (ufp->label != NULL 
            && StringCmp (ufp->label->str, "ProjectID") == 0
            && ufp->choice == 2) {
          sprintf (buf, "%d", ufp->data.intvalue);
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (buf, scp)) {
            tmp = StringSave (buf);
            if (SetStringValue (&tmp, value, existing_text) && IsAllDigits (tmp)) {
              ufp->data.intvalue = atoi (tmp);
              rval = TRUE;
            }
            tmp = MemFree (tmp);
          }
        }
        ufp = ufp->next;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }
  if (!rval && IsStringConstraintEmpty (scp)) {
    sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
    uop = CreateGenomeProjectsDBUserObject ();
    AddIDsToGenomeProjectsDBUserObject (uop, atoi (value), 0);
    sdp->data.ptrvalue = uop;
    rval = TRUE;
  }
  return rval;
}


static Boolean SetCommentDescriptor (SeqDescrPtr sdp, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  CharPtr cp;
  ObjValNodePtr ovp;
  Boolean was_empty;

  if (sdp == NULL) {
    return FALSE;
  }
  
  if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (sdp->data.ptrvalue, scp)) {
    if (StringHasNoText (sdp->data.ptrvalue)) {
      was_empty = TRUE;
    } else {
      was_empty = FALSE;
    }
    cp = sdp->data.ptrvalue;
    if (SetStringValue (&cp, value, existing_text)) {
      rval = TRUE;
    }
    sdp->data.ptrvalue = cp;
    if (was_empty) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = FALSE;
    }
  }

  return rval;
}


static CharPtr s_StringEndsWith (CharPtr str, CharPtr end)
{
  Int4 str_len, end_len;
  if (end == NULL || str == NULL) {
    return NULL;
  }
  str_len = StringLen (str);
  end_len = StringLen (end);
  if (end_len > str_len) {
    return NULL;
  }
  if (StringCmp (str + str_len - end_len, end) == 0) {
    return str + str_len - end_len;
  } else {
    return NULL;
  }
}


static CharPtr DbnameValFromPrefixOrSuffix (CharPtr val)
{
  CharPtr rval = NULL, stop;

  if (val == NULL) {
    return NULL;
  }

  if (StringNCmp (val, "##", 2) == 0) {
    val += 2;
  }
  rval = StringSave (val);
  if ((stop = s_StringEndsWith (rval, "Data-START##")) != NULL
      || (stop = s_StringEndsWith (rval, "-START##")) != NULL
      || (stop = s_StringEndsWith (rval, "-START##")) != NULL
      || (stop = s_StringEndsWith (rval, "START##")) != NULL
      || (stop = s_StringEndsWith (rval, "Data-END##")) != NULL
      || (stop = s_StringEndsWith (rval, "-END##")) != NULL
      || (stop = s_StringEndsWith (rval, "END##")) != NULL) {
    *stop = 0;
  }
  return rval;
}


static Boolean IsUserFieldStructuredCommentPrefixOrSuffix (UserFieldPtr ufp)
{
  if (ufp == NULL || ufp->label == NULL) {
    return FALSE;
  } else if (StringCmp (ufp->label->str, "StructuredCommentPrefix") == 0 
    || StringCmp (ufp->label->str, "StructuredCommentSuffix") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr GetStructuredCommentFieldFromUserObject (UserObjectPtr uop, StructuredCommentFieldPtr field, StringConstraintPtr scp)
{
  UserFieldPtr      curr;
  CharPtr           rval = NULL;

  if (!IsUserObjectStructuredComment(uop) || field == NULL) {
    return NULL;
  }

  if (field->choice == StructuredCommentField_database) {
    for (curr = uop->data; curr != NULL && rval == NULL; curr = curr->next) {
      if (IsUserFieldStructuredCommentPrefixOrSuffix(curr) && curr->choice == 1) {
        rval = DbnameValFromPrefixOrSuffix (curr->data.ptrvalue);
        if (!IsStringConstraintEmpty (scp) &&  !DoesStringMatchConstraint (rval, scp)) {
          rval = MemFree (rval);
        }
      }
    }
  } else if (field->choice == StructuredCommentField_named) {
    for (curr = uop->data; curr != NULL && rval == NULL; curr = curr->next) {
      if (curr->label != NULL && StringICmp (curr->label->str, field->data.ptrvalue) == 0) {
        if (curr->choice == 1) {
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (curr->data.ptrvalue, scp)) {
            rval = StringSave (curr->data.ptrvalue);
          }
        }
      }
    }
  } else if (field->choice == StructuredCommentField_field_name) {
    for (curr = uop->data; curr != NULL && rval == NULL; curr = curr->next) {
      if (!IsUserFieldStructuredCommentPrefixOrSuffix (curr)
          && DoesObjectIdMatchStringConstraint(curr->label, scp)) {
        rval = GetObjectIdString (curr->label);
      }
    }
  }
  return rval;
}


static Boolean RemoveStructuredCommentFieldFromUserObject (UserObjectPtr uop, ValNodePtr field, StringConstraintPtr scp)
{
  UserFieldPtr      curr, prev = NULL, ufp_next;
  Boolean           rval = FALSE, do_remove;
  CharPtr           val;

  if (!IsUserObjectStructuredComment(uop) || field == NULL) {
    return FALSE;
  }

  if (field->choice == StructuredCommentField_database) {
    for (curr = uop->data; curr != NULL; curr = ufp_next) {
      do_remove = FALSE;
      ufp_next = curr->next;
      if (IsUserFieldStructuredCommentPrefixOrSuffix (curr)
              && curr->choice == 1) {
        val = DbnameValFromPrefixOrSuffix (curr->data.ptrvalue);
        if (IsStringConstraintEmpty (scp) || !DoesStringMatchConstraint (val, scp)) {
          do_remove = TRUE;
        }
        val = MemFree (val);
      }
      if (do_remove) {
        if (prev == NULL) {
          uop->data = curr->next;
        } else {
          prev->next = curr->next;
        }
        curr->next = NULL;
        curr = UserFieldFree (curr);
        rval = TRUE;
      } else {
        prev = curr;
      }
    }
  } else if (field->choice == StructuredCommentField_named) {
    for (curr = uop->data; curr != NULL; curr = ufp_next) {
      do_remove = FALSE;
      ufp_next = curr->next;
      if (curr->label != NULL && StringICmp (curr->label->str, field->data.ptrvalue) == 0) {
        if (curr->choice == 1) {
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (curr->data.ptrvalue, scp)) {
            do_remove = TRUE;
          }
        }
      }
      if (do_remove) {
        if (prev == NULL) {
          uop->data = curr->next;
        } else {
          prev->next = curr->next;
        }
        curr->next = NULL;
        curr = UserFieldFree (curr);
        rval = TRUE;
      } else {
        prev = curr;
      }
    }
  } else if (field->choice == StructuredCommentField_field_name) {
    for (curr = uop->data; curr != NULL; curr = ufp_next) {
      do_remove = FALSE;
      ufp_next = curr->next;
      if (!IsUserFieldStructuredCommentPrefixOrSuffix (curr) && DoesObjectIdMatchStringConstraint (curr->label, scp)) {
        if (prev == NULL) {
          uop->data = curr->next;
        } else {
          prev->next = curr->next;
        }
        curr->next = NULL;
        curr = UserFieldFree (curr);
        rval = TRUE;
      } else {
        prev = curr;
      }
    }
  }
  return rval;
}


static Boolean SetStructuredCommentFieldOnUserObject (UserObjectPtr uop, StructuredCommentFieldPtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  UserFieldPtr      curr, first = NULL, last = NULL, ufp;
  Boolean           rval = FALSE;
  CharPtr           oldval, newval, fmt;
  CharPtr           prefix_fmt = "##%sData-START##";
  CharPtr           suffix_fmt = "##%sData-END##";

  if (!IsUserObjectStructuredComment(uop) || field == NULL) {
    return FALSE;
  }

  if (field->choice == StructuredCommentField_database) {
    first = uop->data;
    curr = first;
    while (curr != NULL) {
      if (IsUserFieldStructuredCommentPrefixOrSuffix (curr)
              && curr->choice == 1) {
        oldval = DbnameValFromPrefixOrSuffix (curr->data.ptrvalue);
        if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (oldval, scp)) {
          if (StringCmp (curr->label->str, "StructuredCommentPrefix") == 0) {
            fmt = prefix_fmt;
          } else {
            fmt = suffix_fmt;
          }
          SetStringValue (&oldval, value, existing_text);
          newval = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (oldval)));
          sprintf (newval, fmt, oldval);
          curr->data.ptrvalue = MemFree (curr->data.ptrvalue);
          curr->data.ptrvalue = newval;
          rval = TRUE;
        }
        oldval = MemFree (oldval);
      }
      last = curr;
      curr = curr->next;
    }
    if (!rval && IsStringConstraintEmpty (scp)) {
      /* make prefix */
      curr = UserFieldNew ();
      curr->label = ObjectIdNew ();
      curr->label->str = StringSave ("StructuredCommentPrefix");
      curr->choice = 1;
      newval = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix_fmt) + StringLen (value)));
      sprintf (newval, prefix_fmt, value);
      curr->data.ptrvalue = newval;
      curr->next = first;
      uop->data = curr;
      first = curr;

      /* make suffix */
      curr = UserFieldNew ();
      curr->label = ObjectIdNew ();
      curr->label->str = StringSave ("StructuredCommentSuffix");
      curr->choice = 1;
      newval = (CharPtr) MemNew (sizeof (Char) * (StringLen (suffix_fmt) + StringLen (value)));
      sprintf (newval, suffix_fmt, value);
      curr->data.ptrvalue = newval;
      if (last == NULL) {
        first->next = curr;
      } else {
        last->next = curr;
      }
      rval = TRUE;
    }
  } else if (field->choice == StructuredCommentField_named) {
    last = uop->data;
    for (curr = uop->data; curr != NULL; curr = curr->next) {
      if (curr->label != NULL && StringICmp (curr->label->str, field->data.ptrvalue) == 0) {
        if (curr->choice == 1) {
          if (IsStringConstraintEmpty (scp) || DoesStringMatchConstraint (curr->data.ptrvalue, scp)) {
            newval = (CharPtr) curr->data.ptrvalue;
            SetStringValue (&newval, value, existing_text);
            curr->data.ptrvalue = newval;
            rval = TRUE;
          }
        }
      }
      last = curr;
    }
    if (!rval && IsStringConstraintEmpty (scp)) {
      curr = UserFieldNew ();
      curr->label = ObjectIdNew ();
      curr->label->str = StringSave (field->data.ptrvalue);
      curr->choice = 1;
      curr->data.ptrvalue = StringSave (value);
      if (last == NULL) {
        uop->data = curr;
      } else {
        last->next = curr;
      }
      rval = TRUE;
    }
  } else if (field->choice == StructuredCommentField_field_name) {
    last = uop->data;
    for (curr = uop->data; curr != NULL; curr = curr->next) {
      if (!IsUserFieldStructuredCommentPrefixOrSuffix (curr)) {
        if (DoesObjectIdMatchStringConstraint (curr->label, scp)) {
          rval = SetObjectIdString (curr->label, value, existing_text);
        }
        last = curr;
      }
    }
    if (!rval && IsStringConstraintEmpty (scp)) {
      curr = UserFieldNew ();
      curr->label = ObjectIdNew ();
      curr->label->str = StringSave (value);
      curr->choice = 1;
      curr->data.ptrvalue = StringSave ("");
      if (last == NULL) {
        ufp = uop->data;
        curr->next = ufp->next;
        ufp->next = curr;
      } else {
        curr->next = last->next;
        last->next = curr;
      }
      rval = TRUE;
    }
  }
  return rval;
}






/* The following functions are used for getting and setting various types of data
 * in publications.
 */


static CharPtr legalMonths [] = {
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov",
  "Dec",
  NULL
};


static DatePtr ReadDateFromString (CharPtr date_str)
{
  Char      ch;
  Int2      i;
  CharPtr   ptr1, ptr2, month = NULL, day = NULL, year = NULL;
  CharPtr   str;
  Int4      day_val = 0;
  Uint1     month_num = 0;
  Int4      val;
  Int4      year_val = 0;
  DatePtr   dp = NULL;
  Boolean   critical_error = FALSE;

  if (StringHasNoText (date_str)) return NULL;

  str = StringSave (date_str);
  ptr1 = StringChr (str, '-');
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = StringChr (ptr1, '-');
    if (ptr2 != NULL) {
      *ptr2 = '\0';
      ptr2++;
      day = str;
      month = ptr1;
      year = ptr2;
    } else {
      month = str;
      year = ptr1;
    }
  } else {
    year = str;
  }

  if (day != NULL) {
    if (sscanf (day, "%ld", &day_val) != 1 || day_val < 1 || day_val > 31) {
      critical_error = TRUE;
    }
  }

  if (month != NULL) {
    for (i = 0; legalMonths [i] != NULL; i++) {
      if (StringCmp (month, legalMonths [i]) == 0) {
        month_num = i + 1;
        break;
      }
    }
    if (legalMonths [i] == NULL) critical_error = TRUE;
  }

  if (year != NULL) {
    ptr1 = year;
    ch = *ptr1;
    while (ch != '\0') {
      if (! (IS_DIGIT (ch))) critical_error = TRUE;
      ptr1++;
      ch = *ptr1;
    }
    if (sscanf (year, "%ld", &val) == 1) {
      if (val < 1700 || val > 2100) critical_error = TRUE;
      year_val = val - 1900;
    }
    else
    {
      critical_error = TRUE;
    }
  }

  str = MemFree (str);

  if (!critical_error) {
    dp = DateNew();
    dp->data[0] = 1;
    dp->data[1] = (Uint1) year_val;
    dp->data[2] = month_num;
    dp->data[3] = (Uint1) day_val;
  }
  return dp;
}


static CharPtr GetAuthorStringEx (AuthorPtr author, Boolean use_initials)
{
  CharPtr str = NULL;
  NameStdPtr n;
  Int4       len;
  Boolean    has_middle = FALSE;

  if (author == NULL || author->name == NULL) return NULL;

  switch (author->name->choice) {
    case 1: /* dbtag */
      str = GetDbtagString (author->name->data);
      break;
    case 2: /* name */ 
      n = (NameStdPtr) author->name->data;
      if (n != NULL) {
        if (use_initials) {
          len = StringLen (n->names[0]) + StringLen (n->names[4]) + 2;
          str = (CharPtr) MemNew (sizeof (Char) * (len));
          sprintf (str, "%s%s", StringHasNoText (n->names[4]) ? "" : n->names[4],
                                StringHasNoText (n->names[0]) ? "" : n->names[0]);
        } else {
          len = StringLen (n->names[1]) + StringLen (n->names[0]) + 2;
          if (StringLen (n->names[4]) > 2) {
            len += StringLen (n->names[4]) - 1;
            has_middle = TRUE;
          }
          str = (CharPtr) MemNew (sizeof (Char) * (len));
          sprintf (str, "%s%s%s%s%s", 
                  StringHasNoText (n->names[1]) ? "" : n->names[1],
                  StringHasNoText (n->names[1]) ? "" : " ",
                  has_middle ? n->names[4] + 2 : "",
                  has_middle ? " " : "",
                  StringHasNoText (n->names[0]) ? "" : n->names[0]);
        }
      }
      break;
    case 3: /* ml */
    case 4: /* str */
    case 5: /* consortium */
      str = StringSave (author->name->data);
      break;
  }
  return str;
}


static CharPtr GetAuthorString (AuthorPtr author)
{
  return GetAuthorStringEx (author, FALSE);
}


static CharPtr GetAuthorListStringEx (AuthListPtr alp, StringConstraintPtr scp, Boolean use_initials)
{
  CharPtr str = NULL, tmp;
  Int4    len = 0;
  ValNodePtr list = NULL, vnp;

  if (alp == NULL) return NULL;

  switch (alp->choice) {
    case 1:
      for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
        tmp = GetAuthorStringEx (vnp->data.ptrvalue, use_initials);
        if (tmp != NULL) {
          if (DoesStringMatchConstraint (tmp, scp)) {
            ValNodeAddPointer (&list, 0, tmp);
            len += StringLen (tmp) + 2;
          } else {
            tmp = MemFree (tmp);
          }
        }
      }
      break;
    case 2:
    case 3:
      for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
        if (vnp->data.ptrvalue != NULL && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
          ValNodeAddPointer (&list, 0, StringSave (vnp->data.ptrvalue));
          len += StringLen (vnp->data.ptrvalue) + 2;
        }
      }
      break;
  }

  if (len > 0) {
    str = (CharPtr) MemNew (sizeof (Char) * (len + 1));
    for (vnp = list; vnp != NULL; vnp = vnp->next) {
      StringCat (str, vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        StringCat (str, ", ");
      }
    }
  }
  return str;
}


static CharPtr GetAuthorListString (AuthListPtr alp, StringConstraintPtr scp)
{
  return GetAuthorListStringEx (alp, scp, FALSE);
}


static Boolean RemoveAuthorListString (AuthListPtr alp, StringConstraintPtr scp)
{
  CharPtr tmp;
  Boolean rval = FALSE;
  ValNodePtr vnp, vnp_next, vnp_prev = NULL;

  if (alp == NULL) return FALSE;

  switch (alp->choice) {
    case 1:
      for (vnp = alp->names; vnp != NULL; vnp = vnp_next) {
        vnp_next = vnp->next;
        tmp = GetAuthorString (vnp->data.ptrvalue);
        if (tmp != NULL) {
          if (DoesStringMatchConstraint (tmp, scp)) {
            if (vnp_prev == NULL) {
              alp->names = vnp->next;
            } else {
              vnp_prev->next = vnp->next;
            }
            vnp->next = NULL;
            vnp->data.ptrvalue = AuthorFree (vnp->data.ptrvalue);
            vnp = ValNodeFree (vnp);
            rval = TRUE;
          } else {
            vnp_prev = vnp;
          }
          tmp = MemFree (tmp);
        } else {
          vnp_prev = vnp;
        }
      }
      break;
    case 2:
    case 3:
      for (vnp = alp->names; vnp != NULL; vnp = vnp_next) {
        vnp_next = vnp->next;
        if (vnp->data.ptrvalue != NULL && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {        
          if (vnp_prev == NULL) {
            alp->names = vnp->next;
          } else {
            vnp_prev->next = vnp->next;
          }
          vnp->next = NULL;
          vnp = ValNodeFreeData (vnp);
          rval = TRUE;
        } else {
          vnp_prev = vnp;
        }
      }
      break;
  }

  return rval;
}


static NameStdPtr ReadNameFromString (CharPtr str, CharPtr PNTR next_name)
{
  CharPtr cp_end, cp_space;
  CharPtr p_repl1 = NULL, p_repl2 = NULL, p_repl3 = NULL;
  Char    ch_r1, ch_r2, ch_r3;
  NameStdPtr n;

  if (StringHasNoText (str)) 
  {
    if (next_name != NULL)
    {
      *next_name = NULL;
    }
    return NULL;
  }

  /* skip over any leading spaces */
  str += StringSpn (str, " \t");

  /* skip over "and" if found */
  if (StringNCmp (str, "and ", 4) == 0)
  {
    str += 4;
  }
  if (StringHasNoText (str)) {
    str = MemFree (str);
    return NULL;
  }

  cp_end = StringChr (str, ',');
  if (cp_end != NULL)
  {
    p_repl1 = cp_end;
    ch_r1 = *p_repl1;
    *cp_end = 0;
    if (next_name != NULL)
    {
      if (StringHasNoText (cp_end + 1))
      {
        *next_name = NULL;
      }
      else
      {
        *next_name = cp_end + 1;
      }
    }
  }
  else if (next_name != NULL)
  {
    *next_name = NULL;
  }

  n = NameStdNew ();  
  /* look for elements in name */
  cp_space = StringRChr (str, ' ');
  if (cp_space == NULL)
  {
    n->names[0] = StringSave (str);
  }
  else
  {
    n->names[0] = StringSave (cp_space + 1);
    while (isspace (*cp_space))
    {
      cp_space--;
    }
    p_repl2 = cp_space + 1;
    ch_r2 = *p_repl2;
    *(cp_space + 1) = 0;
    cp_space = StringChr (str, ' ');
    if (cp_space == NULL)
    {
       n->names[1] = StringSave (str);
       n->names[4] = (CharPtr) MemNew (sizeof (Char) * 3);
       sprintf (n->names[4], "%c.", *(n->names[1]));
    }
    else
    {
      p_repl3 = cp_space;
      ch_r3 = *p_repl3;
      *(cp_space) = 0;
      n->names[1] = StringSave (str);

      cp_space++;
      while (isspace (*cp_space))
      {
        cp_space++;
      }
      
      n->names[4] = (CharPtr) MemNew (sizeof (Char) * (4 + StringLen (cp_space)));
      sprintf (n->names[4], "%c.%s.", *(n->names[1]), cp_space);
    }
  }

  if (p_repl1 != NULL) {
    *p_repl1 = ch_r1;
  }
  if (p_repl2 != NULL) {
    *p_repl2 = ch_r2;
  }
  if (p_repl3 != NULL) {
    *p_repl3 = ch_r3;
  }
 
  return n;
}


static ValNodePtr ReadNameListFromString (CharPtr value)
{
  ValNodePtr names = NULL;
  AuthorPtr  ap;
  NameStdPtr n;
  CharPtr    next_cp, cp;

  cp = value;
  next_cp = NULL;
  while (cp != NULL) {
    n = ReadNameFromString (cp, &next_cp);
    if (n != NULL) {
      ap = AuthorNew ();
      ap->name = PersonIdNew ();
      ap->name->choice = 2;
      ap->name->data = n;
      ValNodeAddPointer (&names, 1, ap);
    }
    cp = next_cp;
  }
  return names;
}


static ValNodePtr FreeNameList (Uint1 choice, ValNodePtr name_list)
{
  ValNodePtr curr, next;

  curr = name_list;
  while (curr != NULL) {
    if (choice == 1)    /* std type */
        AuthorFree((AuthorPtr) curr->data.ptrvalue);
    else                      /* ml or str */
        MemFree(curr->data.ptrvalue);

    next = curr->next;
    MemFree(curr);
    curr = next;
  }
  return curr;
}




static Boolean SetAuthorListFromString (AuthListPtr alp, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  ValNodePtr name_list = NULL, vnp, vnp_prev, vnp_next, vnp_tmp;
  CharPtr    tmp;
  Boolean    rval = FALSE, found, ok_to_set = FALSE;

  if (alp == NULL || StringHasNoText (value)) return FALSE;

  /* can only combine lists if existing list is same type */
  if (alp->names == NULL || alp->choice == 1) {
    ok_to_set = TRUE;
  } else {
    switch (existing_text) {
      case ExistingTextOption_replace_old:
        if (IsStringConstraintEmpty (scp)) {
          ok_to_set = TRUE;
        }
        break;
      case ExistingTextOption_append_space:
      case ExistingTextOption_append_colon:
      case ExistingTextOption_append_none:
      case ExistingTextOption_prefix_space:
      case ExistingTextOption_prefix_colon:
      case ExistingTextOption_prefix_none:
        ok_to_set = TRUE;
        break;
    }
  }
  if (!ok_to_set) {
    return FALSE;
  }

  if (alp->names == NULL && IsStringConstraintEmpty (scp)) {
    /* no prior values - just add new list */
    name_list = ReadNameListFromString (value);
    if (name_list != NULL) {
      ValNodeLink (&alp->names, name_list);
      alp->choice = 1;
      rval = TRUE;
    }
  } else {
    switch (existing_text) {
      case ExistingTextOption_append_semi:
        name_list = ReadNameListFromString (value);
        if (IsStringConstraintEmpty (scp)) {
          /* append to list */
          ValNodeLink (&(alp->names), name_list);
          rval = TRUE;
        } else {
          /* insert in list after first match */
          vnp = alp->names;
          found = FALSE;
          while (vnp != NULL && !found) {
            tmp = GetAuthorString (vnp->data.ptrvalue);
            if (tmp != NULL && DoesStringMatchConstraint (tmp, scp)) {
              found = TRUE;
            }
            tmp = MemFree (tmp);
            if (!found) {
              vnp = vnp->next;
            }
          }
          if (found) {
            ValNodeLink (&name_list, vnp->next);
            vnp->next = name_list;
            rval = TRUE;
          }
        }
        break;
      case ExistingTextOption_prefix_semi:
        name_list = ReadNameListFromString (value);
        if (IsStringConstraintEmpty (scp)) {
          /* prepend to list */
          ValNodeLink (&name_list, alp->names);
          alp->names = name_list;
          rval = TRUE;
        } else {
          /* insert in list before first match */
          vnp = alp->names;
          vnp_prev = NULL;
          found = FALSE;
          while (vnp != NULL && !found) {
            tmp = GetAuthorString (vnp->data.ptrvalue);
            if (tmp != NULL && DoesStringMatchConstraint (tmp, scp)) {
              found = TRUE;
            }
            tmp = MemFree (tmp);
            if (!found) {
              vnp_prev = vnp;
              vnp = vnp->next;
            }
          }
          if (found) {
            if (vnp_prev == NULL) {
              ValNodeLink (&name_list, alp->names);
              alp->names = name_list;
            } else {
              ValNodeLink (&name_list, vnp_prev->next);
              vnp_prev->next = name_list;
            }
            rval = TRUE;
          }
        }
        break;
      case ExistingTextOption_replace_old:
        name_list = ReadNameListFromString (value);
        if (IsStringConstraintEmpty (scp)) {
          /* replace entire list */
          alp->names = FreeNameList (alp->choice, alp->names);
          alp->names = name_list;
          alp->choice = 1;
          rval = TRUE;
        } else {
          /* replace first author that matches with new match, remove others that match */
          vnp = alp->names;
          vnp_prev = NULL;
          found = FALSE;
          while (vnp != NULL) {
            vnp_next = vnp->next;
            tmp = GetAuthorString (vnp->data.ptrvalue);
            if (tmp != NULL && DoesStringMatchConstraint (tmp, scp)) {
              if (found) {
                if (vnp_prev == NULL) {
                  alp->names = vnp->next;
                } else {
                  vnp_prev->next = vnp->next;
                }
              } else {
                vnp_tmp = name_list;
                while (vnp_tmp->next != NULL) {
                  vnp_tmp = vnp_tmp->next;
                }
                ValNodeLink (&name_list, vnp->next);
                if (vnp_prev == NULL) {
                  alp->names = name_list;
                } else {
                  vnp_prev->next = name_list;
                }
                vnp_prev = vnp_tmp;
                found = TRUE;
                rval = TRUE;
              } 
              vnp->next = NULL;
              vnp = FreeNameList (alp->choice, vnp);
            } else {
              vnp_prev = vnp;
            }
            tmp = MemFree (tmp);
            vnp = vnp_next;
          }
        }
        break;
      case ExistingTextOption_append_space:
      case ExistingTextOption_append_colon:
      case ExistingTextOption_append_none:
      case ExistingTextOption_prefix_space:
      case ExistingTextOption_prefix_colon:
      case ExistingTextOption_prefix_none:
        vnp_prev = NULL;
        for (vnp = alp->names; vnp != NULL; vnp = vnp_next) {
          vnp_next = vnp->next;
          if (alp->choice == 1) {
            tmp = GetAuthorString (vnp->data.ptrvalue);
            if (tmp != NULL && DoesStringMatchConstraint (tmp, scp)
                && SetStringValue (&tmp, value, existing_text)) {
              name_list = ReadNameListFromString (tmp);
              if (name_list != NULL) {
                vnp_tmp = name_list;
                while (vnp_tmp->next != NULL) {
                  vnp_tmp = vnp_tmp->next;
                }
                ValNodeLink (&name_list, vnp_next);
                if (vnp_prev == NULL) {
                  alp->names = name_list;
                } else {
                  vnp_prev->next = name_list;
                }
                vnp_prev = vnp_tmp;
                vnp->next = NULL;
                vnp = FreeNameList (alp->choice, vnp);
                rval = TRUE;
                name_list = NULL;
              } else {
                vnp_prev = vnp;
              }
            } else {
              vnp_prev = vnp;
            }
            tmp = MemFree (tmp);
          } else {
            if (vnp->data.ptrvalue != NULL && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
              tmp = (CharPtr) vnp->data.ptrvalue;
              rval |= SetStringValue (&tmp, value, existing_text);
              vnp->data.ptrvalue = tmp;
            }
          }
        }
        break;
    }
  }
  if (!rval && name_list != NULL) {
    name_list = FreeNameList (1, vnp);
  }
  return rval;
}


static CharPtr GetPubFieldFromAffil (AffilPtr ap, Int4 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;

  if (ap == NULL) return NULL;

  switch (field) {
    case Publication_field_affiliation:
      if (!StringHasNoText (ap->affil) && DoesStringMatchConstraint (ap->affil, scp)) {
        str = StringSave (ap->affil);
      }
      break;
    case Publication_field_affil_div:
      if (!StringHasNoText (ap->div) && DoesStringMatchConstraint (ap->div, scp)) {
        str = StringSave (ap->div);
      }
      break;
    case Publication_field_affil_city:
      if (!StringHasNoText (ap->city) && DoesStringMatchConstraint (ap->city, scp)) {
        str = StringSave (ap->city);
      }
      break;
    case Publication_field_affil_sub:
      if (!StringHasNoText (ap->sub) && DoesStringMatchConstraint (ap->sub, scp)) {
        str = StringSave (ap->sub);
      }
      break;
    case Publication_field_affil_country:
      if (!StringHasNoText (ap->country) && DoesStringMatchConstraint (ap->country, scp)) {
        str = StringSave (ap->country);
      }
      break;
    case Publication_field_affil_street:
      if (!StringHasNoText (ap->street) && DoesStringMatchConstraint (ap->street, scp)) {
        str = StringSave (ap->street);
      }
      break;
    case Publication_field_affil_email:
      if (!StringHasNoText (ap->email) && DoesStringMatchConstraint (ap->email, scp)) {
        str = StringSave (ap->email);
      }
      break;
    case Publication_field_affil_fax:
      if (!StringHasNoText (ap->fax) && DoesStringMatchConstraint (ap->fax, scp)) {
        str = StringSave (ap->fax);
      }
      break;
    case Publication_field_affil_phone:
      if (!StringHasNoText (ap->phone) && DoesStringMatchConstraint (ap->phone, scp)) {
        str = StringSave (ap->phone);
      }
      break;
    case Publication_field_affil_zipcode:
      if (!StringHasNoText (ap->postal_code) && DoesStringMatchConstraint (ap->postal_code, scp)) {
        str = StringSave (ap->postal_code);
      }
      break;
  }
  return str;
}


static Boolean RemovePubFieldFromAffil (AffilPtr ap, Int4 field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  if (ap == NULL) return FALSE;

  switch (field) {
    case Publication_field_affiliation:
      if (!StringHasNoText (ap->affil) && DoesStringMatchConstraint (ap->affil, scp)) {
        ap->affil = MemFree (ap->affil);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_div:
      if (!StringHasNoText (ap->div) && DoesStringMatchConstraint (ap->div, scp)) {
        ap->div = MemFree (ap->div);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_city:
      if (!StringHasNoText (ap->city) && DoesStringMatchConstraint (ap->city, scp)) {
        ap->city = MemFree (ap->city);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_sub:
      if (!StringHasNoText (ap->sub) && DoesStringMatchConstraint (ap->sub, scp)) {
        ap->sub = MemFree (ap->sub);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_country:
      if (!StringHasNoText (ap->country) && DoesStringMatchConstraint (ap->country, scp)) {
        ap->country = MemFree (ap->country);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_street:
      if (!StringHasNoText (ap->street) && DoesStringMatchConstraint (ap->street, scp)) {
        ap->street = MemFree (ap->street);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_email:
      if (!StringHasNoText (ap->email) && DoesStringMatchConstraint (ap->email, scp)) {
        ap->email = MemFree (ap->email);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_fax:
      if (!StringHasNoText (ap->fax) && DoesStringMatchConstraint (ap->fax, scp)) {
        ap->fax = MemFree (ap->fax);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_phone:
      if (!StringHasNoText (ap->phone) && DoesStringMatchConstraint (ap->phone, scp)) {
        ap->phone = MemFree (ap->phone);
        rval = TRUE;
      }
      break;
    case Publication_field_affil_zipcode:
      if (!StringHasNoText (ap->postal_code) && DoesStringMatchConstraint (ap->postal_code, scp)) {
        ap->postal_code = MemFree (ap->postal_code);
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Boolean SetAffilPubField (AffilPtr ap, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  if (ap == NULL) return FALSE;

  switch (field) {
    case Publication_field_affiliation:
      if (!StringHasNoText (ap->affil) && DoesStringMatchConstraint (ap->affil, scp)) {
        rval = SetStringValue (&(ap->affil), value, existing_text);
      }
      break;
    case Publication_field_affil_div:
      if (!StringHasNoText (ap->div) && DoesStringMatchConstraint (ap->div, scp)) {
        rval = SetStringValue (&(ap->div), value, existing_text);
      }
      break;
    case Publication_field_affil_city:
      if (!StringHasNoText (ap->city) && DoesStringMatchConstraint (ap->city, scp)) {
        rval = SetStringValue (&(ap->city), value, existing_text);
      }
      break;
    case Publication_field_affil_sub:
      if (!StringHasNoText (ap->sub) && DoesStringMatchConstraint (ap->sub, scp)) {
        rval = SetStringValue (&(ap->sub), value, existing_text);
      }
      break;
    case Publication_field_affil_country:
      if (!StringHasNoText (ap->country) && DoesStringMatchConstraint (ap->country, scp)) {
        rval = SetStringValue (&(ap->country), value, existing_text);
      }
      break;
    case Publication_field_affil_street:
      if (!StringHasNoText (ap->street) && DoesStringMatchConstraint (ap->street, scp)) {
        rval = SetStringValue (&(ap->street), value, existing_text);
      }
      break;
    case Publication_field_affil_email:
      if (!StringHasNoText (ap->email) && DoesStringMatchConstraint (ap->email, scp)) {
        rval = SetStringValue (&(ap->email), value, existing_text);
      }
      break;
    case Publication_field_affil_fax:
      if (!StringHasNoText (ap->fax) && DoesStringMatchConstraint (ap->fax, scp)) {
        rval = SetStringValue (&(ap->fax), value, existing_text);
      }
      break;
    case Publication_field_affil_phone:
      if (!StringHasNoText (ap->phone) && DoesStringMatchConstraint (ap->phone, scp)) {
        rval = SetStringValue (&(ap->phone), value, existing_text);
      }
      break;
    case Publication_field_affil_zipcode:
      if (!StringHasNoText (ap->postal_code) && DoesStringMatchConstraint (ap->postal_code, scp)) {
        rval = SetStringValue (&(ap->postal_code), value, existing_text);
      }
      break;
  }
  return rval;
}


static CharPtr GetPubFieldFromImprint (ImprintPtr imprint, Int4 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  if (imprint == NULL) return NULL;

  switch (field) {
    case Publication_field_volume:
      if (!StringHasNoText (imprint->volume) && DoesStringMatchConstraint (imprint->volume, scp)) {
        str = StringSave (imprint->volume);
      }
      break;
    case Publication_field_issue:
      if (!StringHasNoText (imprint->issue) && DoesStringMatchConstraint (imprint->issue, scp)) {
        str = StringSave (imprint->issue);
      }
      break;
    case Publication_field_pages:
      if (!StringHasNoText (imprint->pages) && DoesStringMatchConstraint (imprint->pages, scp)) {
        str = StringSave (imprint->pages);
      }
      break;
    case Publication_field_date:
      str = PrintDate (imprint->date);
      if (StringHasNoText (str) || !DoesStringMatchConstraint (str, scp)) {
        str = MemFree (str);
      }
      break;
  }
  return str;
}


static Boolean RemovePubDate (DatePtr PNTR pDate, StringConstraintPtr scp)
{
  CharPtr str;
  Boolean rval = FALSE;

  if (pDate == NULL || *pDate == NULL) {
    return FALSE;
  }

  str = PrintDate (*pDate);
  if (!StringHasNoText (str) && DoesStringMatchConstraint (str, scp)) {
    *pDate = DateFree (*pDate);
    rval = TRUE;
  }
  str = MemFree (str);
  return rval;
}


static Boolean SetPubDate (DatePtr PNTR pDate, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  CharPtr tmp;
  DatePtr dp = NULL;
  Boolean rval = FALSE;

  if (pDate == NULL) {
    return FALSE;
  }
  tmp = PrintDate (*pDate);
  if (DoesStringMatchConstraint (tmp, scp)
      && SetStringValue (&tmp, value, existing_text)) {
    dp = ReadDateFromString (tmp);
    if (dp != NULL) {
      *pDate = DateFree (*pDate);
      *pDate = dp;
      rval = TRUE;
    }
  }
  tmp = MemFree (tmp);
  return rval;
}


static Boolean RemovePubFieldFromImprint (ImprintPtr imprint, Int4 field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  if (imprint == NULL) return FALSE;

  switch (field) {
    case Publication_field_volume:
      if (!StringHasNoText (imprint->volume) && DoesStringMatchConstraint (imprint->volume, scp)) {
        imprint->volume = MemFree (imprint->volume);
        rval = TRUE;
      }
      break;
    case Publication_field_issue:
      if (!StringHasNoText (imprint->issue) && DoesStringMatchConstraint (imprint->issue, scp)) {
        imprint->issue = MemFree (imprint->issue);
        rval = TRUE;
      }
      break;
    case Publication_field_pages:
      if (!StringHasNoText (imprint->pages) && DoesStringMatchConstraint (imprint->pages, scp)) {
        imprint->pages = MemFree (imprint->pages);
        rval = TRUE;
      }
      break;
    case Publication_field_date:
      rval = RemovePubDate (&(imprint->date), scp);
      break;
  }
  return rval;
}


static Boolean SetPubFieldOnImprint (ImprintPtr imprint, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;

  if (imprint == NULL) return FALSE;

  switch (field) {
    case Publication_field_volume:
      if (DoesStringMatchConstraint (imprint->volume, scp)) {
        rval = SetStringValue (&(imprint->volume), value, existing_text);
      }
      break;
    case Publication_field_issue:
      if (!StringHasNoText (imprint->issue) && DoesStringMatchConstraint (imprint->issue, scp)) {
        rval = SetStringValue (&(imprint->issue), value, existing_text);
      }
      break;
    case Publication_field_pages:
      if (!StringHasNoText (imprint->pages) && DoesStringMatchConstraint (imprint->pages, scp)) {
        rval = SetStringValue (&(imprint->pages), value, existing_text);
      }
      break;
    case Publication_field_date:
      rval = SetPubDate (&(imprint->date), scp, value, existing_text);
      break;
  }
  return rval;
}


static void SetValNodeChoices (ValNodePtr list, Uint1 new_choice)
{
  while (list != NULL) {
    list->choice = new_choice;
    list = list->next;
  }
}


static CharPtr GetPubFieldFromCitJour (CitJourPtr cjp, Int4 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  if (cjp == NULL) return NULL;

  switch (field) {
    case Publication_field_journal:
      str = GetFirstValNodeStringMatch (cjp->title, scp);
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      str = GetPubFieldFromImprint (cjp->imp, field, scp);
      break;
  }

  return str;
}


static Boolean RemovePubFieldFromCitJour (CitJourPtr cjp, Int4 field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  if (cjp == NULL) return FALSE;

  switch (field) {
    case Publication_field_journal:
      rval = RemoveValNodeStringMatch (&(cjp->title), scp);
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      rval = RemovePubFieldFromImprint (cjp->imp, field, scp);
      break;
  }

  return rval;
}


static Boolean SetPubFieldOnCitJour (CitJourPtr cjp, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  if (cjp == NULL) return FALSE;

  switch (field) {
    case Publication_field_journal:
      rval = SetStringsInValNodeStringList (&(cjp->title), scp, value, existing_text);
      SetValNodeChoices (cjp->title, 1);
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      rval = SetPubFieldOnImprint (cjp->imp, field, scp, value, existing_text);
      break;
  }

  return rval;
}


static CharPtr GetPubFieldFromCitBook (CitBookPtr cbp, Int4 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;

  if (cbp == NULL) return NULL;

  switch (field) {
    case Publication_field_title:
      str = GetFirstValNodeStringMatch (cbp->title, scp);
      break;
    case Publication_field_authors:
      str = GetAuthorListString (cbp->authors, scp);
      break;
    case Publication_field_authors_initials:
      str = GetAuthorListStringEx (cbp->authors, scp, TRUE);
      break;
    case Publication_field_affiliation:
    case Publication_field_affil_div:
    case Publication_field_affil_city:
    case Publication_field_affil_sub:
    case Publication_field_affil_country:
    case Publication_field_affil_street:
    case Publication_field_affil_email:
    case Publication_field_affil_fax:
    case Publication_field_affil_phone:
    case Publication_field_affil_zipcode:
      if (cbp->authors != NULL) {
        str = GetPubFieldFromAffil (cbp->authors->affil, field, scp);
      }
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      str = GetPubFieldFromImprint (cbp->imp, field, scp);
      break;
  }

  return str;
}


static Boolean RemovePubFieldFromCitBook (CitBookPtr cbp, Int4 field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;

  if (cbp == NULL) return FALSE;

  switch (field) {
    case Publication_field_title:
      rval = RemoveValNodeStringMatch (&(cbp->title), scp);
      break;
    case Publication_field_authors:
      rval = RemoveAuthorListString (cbp->authors, scp);
      break;
    case Publication_field_affiliation:
    case Publication_field_affil_div:
    case Publication_field_affil_city:
    case Publication_field_affil_sub:
    case Publication_field_affil_country:
    case Publication_field_affil_street:
    case Publication_field_affil_email:
    case Publication_field_affil_fax:
    case Publication_field_affil_phone:
    case Publication_field_affil_zipcode:
      if (cbp->authors != NULL) {
        rval = RemovePubFieldFromAffil(cbp->authors->affil, field, scp);
      }
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      rval = RemovePubFieldFromImprint (cbp->imp, field, scp);
      break;
  }

  return rval;
}


static Boolean SetPubFieldOnCitBook (CitBookPtr cbp, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;

  if (cbp == NULL) return FALSE;

  switch (field) {
    case Publication_field_title:
      rval = SetStringsInValNodeStringList (&(cbp->title), scp, value, existing_text);
      SetValNodeChoices (cbp->title, 1);
      break;
    case Publication_field_authors:
      rval = SetAuthorListFromString (cbp->authors, scp, value, existing_text);
      break;
    case Publication_field_affiliation:
    case Publication_field_affil_div:
    case Publication_field_affil_city:
    case Publication_field_affil_sub:
    case Publication_field_affil_country:
    case Publication_field_affil_street:
    case Publication_field_affil_email:
    case Publication_field_affil_fax:
    case Publication_field_affil_phone:
    case Publication_field_affil_zipcode:
      if (cbp->authors != NULL) {
        rval = SetAffilPubField (cbp->authors->affil, field, scp, value, existing_text);
      }
      break;
    case Publication_field_volume:
    case Publication_field_issue:
    case Publication_field_pages:
    case Publication_field_date:
      rval = SetPubFieldOnImprint (cbp->imp, field, scp, value, existing_text);
      break;
  }

  return rval;
}


NLM_EXTERN CharPtr GetPubFieldFromPub (PubPtr the_pub, Int4 field, StringConstraintPtr scp)
{
  CitGenPtr    cgp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  CitJourPtr   cjp;
  CharPtr      str = NULL;

  if (the_pub == NULL || the_pub->data.ptrvalue == NULL) return NULL;
  
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_cit:
          if (!StringHasNoText (cgp->cit) && DoesStringMatchConstraint (cgp->title, scp)) {
            str = StringSave (cgp->cit);
          }
          break;
        case Publication_field_authors:
          str = GetAuthorListString (cgp->authors, scp);
          break;
        case Publication_field_authors_initials:
          str = GetAuthorListStringEx (cgp->authors, scp, TRUE);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cgp->authors != NULL) {
            str = GetPubFieldFromAffil (cgp->authors->affil, field, scp);
          }
          break;
        case Publication_field_journal:
          str = GetFirstValNodeStringMatch (cgp->journal, scp);
          break;
        case Publication_field_volume:
          if (!StringHasNoText (cgp->volume) && DoesStringMatchConstraint (cgp->volume, scp)) {
            str = StringSave (cgp->volume);
          }
          break;
        case Publication_field_issue:
          if (!StringHasNoText (cgp->issue) && DoesStringMatchConstraint (cgp->issue, scp)) {
            str = StringSave (cgp->issue);
          }
          break;
        case Publication_field_pages:
          if (!StringHasNoText (cgp->pages) && DoesStringMatchConstraint (cgp->pages, scp)) {
            str = StringSave (cgp->pages);
          }
          break;
        case Publication_field_date:
          if (cgp->date != NULL) {
            str = PrintDate (cgp->date);
            if (StringHasNoText (str) || !DoesStringMatchConstraint (str, scp)) {
              str = MemFree (str);
            }
          }
          break;
        case Publication_field_serial_number:
          str = GetInt2ValueFromString (cgp->serial_number, scp);
          break;
        case Publication_field_title:
          if (!StringHasNoText (cgp->title) && DoesStringMatchConstraint (cgp->title, scp)) {
            str = StringSave (cgp->title);
          }
          break;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (!StringHasNoText (csp->descr) && DoesStringMatchConstraint (csp->descr, scp)) {
            str = StringSave (csp->descr);
          }
          break;
        case Publication_field_authors:
          str = GetAuthorListString (csp->authors, scp);
          break;
        case Publication_field_authors_initials:
          str = GetAuthorListStringEx (csp->authors, scp, TRUE);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (csp->authors != NULL) {
            str = GetPubFieldFromAffil (csp->authors->affil, field, scp);
          }
          break;
        case Publication_field_date:
          str = PrintDate (csp->date);
          if (StringHasNoText (str) || !DoesStringMatchConstraint (str, scp)) {
            str = MemFree (str);
          }
          break;
      }
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          str = GetFirstValNodeStringMatch (cap->title, scp);
          break;
        case Publication_field_authors:
          str = GetAuthorListString (cap->authors, scp);
          break;
        case Publication_field_authors_initials:
          str = GetAuthorListStringEx (cap->authors, scp, TRUE);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cap->authors != NULL) {
            str = GetPubFieldFromAffil (cap->authors->affil, field, scp);
          }
          break;
        default:
          if (cap->from == 1) {
            str = GetPubFieldFromCitJour (cap->fromptr, field, scp);
          } else if (cap->from == 2) {
            str = GetPubFieldFromCitBook (cap->fromptr, field, scp);
          }
          break;
      }
      break;
    case PUB_Journal:
      cjp = (CitJourPtr) the_pub->data.ptrvalue;
      str = GetPubFieldFromCitJour (cjp, field, scp);
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      str = GetPubFieldFromCitBook (cbp, field, scp);
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (!StringHasNoText (cpp->title) && DoesStringMatchConstraint (cpp->title, scp)) {
            str = StringSave (cpp->title);
          }
          break;
        case Publication_field_authors:
          str = GetAuthorListString (cpp->authors, scp);
          break;
        case Publication_field_authors_initials:
          str = GetAuthorListStringEx (cpp->authors, scp, TRUE);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cpp->authors != NULL) {
            str = GetPubFieldFromAffil (cpp->authors->affil, field, scp);
          }
          break;
      }
      break;
    default :
      break;
  }
  return str;
}


static Boolean RemovePubFieldFromPub (PubPtr the_pub, Int4 field, StringConstraintPtr scp)
{
  CitGenPtr    cgp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  Boolean      rval = FALSE;
  Char         num[15];

  if (the_pub == NULL) return FALSE;
  
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_cit:
          if (!StringHasNoText (cgp->cit) && DoesStringMatchConstraint (cgp->title, scp)) {
            cgp->cit = MemFree (cgp->cit);
            rval = TRUE;
          }
          break;
        case Publication_field_authors:
          rval = RemoveAuthorListString (cgp->authors, scp);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cgp->authors != NULL) {
            rval = RemovePubFieldFromAffil(cgp->authors->affil, field, scp);
          }
          break;
        case Publication_field_journal:
          rval = RemoveValNodeStringMatch (&(cgp->journal), scp);
          break;
        case Publication_field_volume:
          if (!StringHasNoText (cgp->volume) && DoesStringMatchConstraint (cgp->volume, scp)) {
            cgp->volume = MemFree (cgp->volume);
            rval = TRUE;
          }
          break;
        case Publication_field_issue:
          if (!StringHasNoText (cgp->issue) && DoesStringMatchConstraint (cgp->issue, scp)) {
            cgp->issue = MemFree (cgp->issue);
            rval = TRUE;
          }
          break;
        case Publication_field_pages:
          if (!StringHasNoText (cgp->pages) && DoesStringMatchConstraint (cgp->pages, scp)) {
            cgp->pages = MemFree (cgp->pages);
            rval = TRUE;
          }
          break;
        case Publication_field_date:
          rval = RemovePubDate (&(cgp->date), scp);
          break;
        case Publication_field_serial_number:
          if (cgp->serial_number > 0) {
            sprintf (num, "%d", cgp->serial_number);
            if (DoesStringMatchConstraint (num, scp)) {
              cgp->serial_number = 0;
              rval = TRUE;
            }
          }
          break;
        case Publication_field_title:
          if (!StringHasNoText (cgp->title) && DoesStringMatchConstraint (cgp->title, scp)) {
            cgp->title = MemFree (cgp->title);
            rval = TRUE;
          }
          break;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (!StringHasNoText (csp->descr) && DoesStringMatchConstraint (csp->descr, scp)) {
            csp->descr = MemFree (csp->descr);
            rval = TRUE;
          }
          break;
        case Publication_field_authors:
          rval = RemoveAuthorListString (csp->authors, scp);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (csp->authors != NULL) {
            rval = RemovePubFieldFromAffil(csp->authors->affil, field, scp);
          }
          break;
        case Publication_field_date:
          rval = RemovePubDate (&(csp->date), scp);
          break;
      }
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          rval = RemoveValNodeStringMatch (&(cap->title), scp);
          break;
        case Publication_field_authors:
          rval = RemoveAuthorListString (cap->authors, scp);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cap->authors != NULL) {
            rval = RemovePubFieldFromAffil(cap->authors->affil, field, scp);
          }
          break;
        default:
          if (cap->from == 1) {
            rval = RemovePubFieldFromCitJour (cap->fromptr, field, scp);
          } else if (cap->from == 2) {
            rval = RemovePubFieldFromCitBook (cap->fromptr, field, scp);
          }
          break;
      }
      break;
    case PUB_Journal:
      rval = RemovePubFieldFromCitJour (the_pub->data.ptrvalue, field, scp);
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      rval = RemovePubFieldFromCitBook (cbp, field, scp);
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (!StringHasNoText (cpp->title) && DoesStringMatchConstraint (cpp->title, scp)) {
            cpp->title = MemFree (cpp->title);
            rval = TRUE;
          }
          break;
        case Publication_field_authors:
          rval = RemoveAuthorListString (cpp->authors, scp);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cpp->authors != NULL) {
            rval = RemovePubFieldFromAffil(cpp->authors->affil, field, scp);
          }
          break;
      }
      break;
    default :
      break;
  }   
  return rval;
}


static Boolean SetPubFieldOnPub (PubPtr the_pub, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  CitGenPtr    cgp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  Boolean      rval = FALSE;

  if (the_pub == NULL || value == NULL) return FALSE;
  
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_cit:
          if (DoesStringMatchConstraint (cgp->cit, scp)) {
            rval = SetStringValue ( &(cgp->cit), value, existing_text);
          }
          break;
        case Publication_field_authors:
          rval = SetAuthorListFromString (cgp->authors, scp, value, existing_text);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cgp->authors != NULL) {
            rval = SetAffilPubField (cgp->authors->affil, field, scp, value, existing_text);
          }
          break;
        case Publication_field_journal:
          rval = SetStringsInValNodeStringList (&(cgp->journal), scp, value, existing_text);
          SetValNodeChoices (cgp->journal, 1);
          break;
        case Publication_field_volume:
          if (DoesStringMatchConstraint (cgp->volume, scp)) {
            rval = SetStringValue ( &(cgp->volume), value, existing_text);
          }
          break;
        case Publication_field_issue:
          if (DoesStringMatchConstraint (cgp->issue, scp)) {
            rval = SetStringValue ( &(cgp->issue), value, existing_text);
          }
          break;
        case Publication_field_pages:
          if (DoesStringMatchConstraint (cgp->pages, scp)) {
            rval = SetStringValue ( &(cgp->pages), value, existing_text);
          }
          break;
        case Publication_field_date:
          rval = SetPubDate (&(cgp->date), scp, value, existing_text);
          break;
        case Publication_field_serial_number:
          rval = SetInt2ValueWithString (&(cgp->serial_number), value, existing_text);
          break;
        case Publication_field_title:
          if (DoesStringMatchConstraint(cgp->title, scp)) {
            rval = SetStringValue ( &(cgp->title), value, existing_text);
          }
          break;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (DoesStringMatchConstraint (csp->descr, scp)) {
            rval = SetStringValue (&(csp->descr), value, existing_text);
          }
          break;
        case Publication_field_authors:
          rval = SetAuthorListFromString (csp->authors, scp, value, existing_text);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (csp->authors != NULL) {
            rval = SetAffilPubField (csp->authors->affil, field, scp, value, existing_text);
          }
          break;
        case Publication_field_date:
          rval = SetPubDate (&(csp->date), scp, value, existing_text);
          break;
      }
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          rval = SetStringsInValNodeStringList (&(cap->title), scp, value, existing_text);
          SetValNodeChoices (cap->title, 1);
          break;
        case Publication_field_authors:
          rval = SetAuthorListFromString (cap->authors, scp, value, existing_text);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cap->authors != NULL) {
            rval = SetAffilPubField (cap->authors->affil, field, scp, value, existing_text);
          }
          break;
        default:
          if (cap->from == 1) {
            rval = SetPubFieldOnCitJour (cap->fromptr, field, scp, value, existing_text);
          } else if (cap->from == 2) {
            rval = SetPubFieldOnCitBook (cap->fromptr, field, scp, value, existing_text);
          }
          break;
      }
      break;
    case PUB_Journal:
      rval = SetPubFieldOnCitJour (the_pub->data.ptrvalue, field, scp, value, existing_text);
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      rval = SetPubFieldOnCitBook (cbp, field, scp, value, existing_text);
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      switch (field) {
        case Publication_field_title:
          if (DoesStringMatchConstraint(cpp->title, scp)) {
            rval = SetStringValue ( &(cpp->title), value, existing_text);
          }
          break;
        case Publication_field_authors:
          rval = SetAuthorListFromString (cpp->authors, scp, value, existing_text);
          break;
        case Publication_field_affiliation:
        case Publication_field_affil_div:
        case Publication_field_affil_city:
        case Publication_field_affil_sub:
        case Publication_field_affil_country:
        case Publication_field_affil_street:
        case Publication_field_affil_email:
        case Publication_field_affil_fax:
        case Publication_field_affil_phone:
        case Publication_field_affil_zipcode:
          if (cpp->authors != NULL) {
            rval = SetAffilPubField (cpp->authors->affil, field, scp, value, existing_text);
          }
          break;
      }
      break;
    default :
      break;
  }
  return rval;
}



static CharPtr GetPubFieldFromObject (Uint1 choice, Pointer data, Int4 field, StringConstraintPtr scp)
{
  CharPtr rval = NULL;
  PubdescPtr pdp = NULL;
  PubPtr     pub;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;

  if (data == NULL) return NULL;
  if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = sfp->data.value.ptrvalue;
    }
  } else if (choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) data;
    if (sdp->choice == Seq_descr_pub) {
      pdp = sdp->data.ptrvalue;
    }
  }

  if (pdp == NULL) return NULL;
  for (pub = pdp->pub; pub != NULL && rval == NULL; pub = pub->next) {
    rval = GetPubFieldFromPub (pub, field, scp);
  }

  return rval;
}


static Boolean RemovePubFieldFromObject (Uint1 choice, Pointer data, Int4 field, StringConstraintPtr scp)
{
  Boolean    rval = FALSE;
  PubdescPtr pdp = NULL;
  PubPtr     pub;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;

  if (data == NULL) return FALSE;
  if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = sfp->data.value.ptrvalue;
    }
  } else if (choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) data;
    if (sdp->choice == Seq_descr_pub) {
      pdp = sdp->data.ptrvalue;
    }
  }

  if (pdp == NULL) return FALSE;

  for (pub = pdp->pub; pub != NULL; pub = pub->next) {
    rval |= RemovePubFieldFromPub (pub, field, scp);
  }
  return rval;
}


static Boolean SetPubFieldOnObject (Uint1 choice, Pointer data, Int4 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean    rval = FALSE;
  PubdescPtr pdp = NULL;
  PubPtr     pub;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;

  if (data == NULL) return FALSE;
  if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = sfp->data.value.ptrvalue;
    }
  } else if (choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) data;
    if (sdp->choice == Seq_descr_pub) {
      pdp = sdp->data.ptrvalue;
    }
  }

  if (pdp == NULL) return FALSE;

  for (pub = pdp->pub; pub != NULL; pub = pub->next) {
    rval |= SetPubFieldOnPub (pub, field, scp, value, existing_text);
  }
  return rval;
}



NLM_EXTERN Uint1 FieldTypeFromAECRAction (AECRActionPtr action)
{
  Uint1 field_type = 0;
  ApplyActionPtr a;
  EditActionPtr  e;
  ConvertActionPtr v;
  CopyActionPtr c;
  SwapActionPtr s;
  RemoveActionPtr r;
  AECRParseActionPtr p;

  if (action == NULL || action->action == NULL || action->action->data.ptrvalue == NULL) {
    return 0;
  }
  switch (action->action->choice) {
    case ActionChoice_apply:
      a = (ApplyActionPtr) action->action->data.ptrvalue;
      if (a->field != NULL) {
        field_type = a->field->choice;
      }
      break;
    case ActionChoice_edit:
      e = (EditActionPtr) action->action->data.ptrvalue;
      if (e->field != NULL) {
        field_type = e->field->choice;
      }
      break;
    case ActionChoice_convert:
      v = (ConvertActionPtr) action->action->data.ptrvalue;
      if (v->fields != NULL) {
        field_type = FieldTypeChoiceFromFieldPairTypeChoice (v->fields->choice);
      }
      break;
    case ActionChoice_copy:
      c = (CopyActionPtr) action->action->data.ptrvalue;
      if (c->fields != NULL) {
        field_type = FieldTypeChoiceFromFieldPairTypeChoice (c->fields->choice);
      }
      break;
    case ActionChoice_swap:
      s = (SwapActionPtr) action->action->data.ptrvalue;
      if (s->fields != NULL) {
        field_type = FieldTypeChoiceFromFieldPairTypeChoice (s->fields->choice);
      }
      break;
    case ActionChoice_remove:
      r = (RemoveActionPtr) action->action->data.ptrvalue;
      if (r->field != NULL) {
        field_type = r->field->choice;
      }
      break;
    case ActionChoice_parse:
      p = (AECRParseActionPtr) action->action->data.ptrvalue;
      if (p->fields != NULL) {
        field_type = FieldTypeChoiceFromFieldPairTypeChoice (p->fields->choice);
      }
      break;
  }
  return field_type;
}

typedef struct pubserialnumber {
  BioseqPtr bsp;
  Int4      serial_number;
  ValNodePtr min_pub;
} PubSerialNumberData, PNTR PubSerialNumberPtr;


static PubSerialNumberPtr PubSerialNumberNew ()
{
  PubSerialNumberPtr psn;

  psn = (PubSerialNumberPtr) MemNew (sizeof (PubSerialNumberData));
  psn->bsp = NULL;
  psn->serial_number = 0;
  psn->min_pub = NULL;

  return psn;
}


static PubSerialNumberPtr PubSerialNumberFree (PubSerialNumberPtr psn)
{
  if (psn != NULL) {
    psn->min_pub = PubSetFree (psn->min_pub);
    psn = MemFree (psn);
  }
  return psn;
}


NLM_EXTERN ValNodePtr PubSerialNumberListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = PubSerialNumberFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void CaptureRefBlockSerialNumbers
(CharPtr str,
 Pointer userdata,
 BlockType blocktype,
 Uint2 entityID,
 Uint2 itemtype,
 Uint4 itemID)
{
  CharPtr          cp;
  Int4             serial_number;
  ValNodePtr       vnp;
  BioseqPtr        bsp = NULL;
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  PubSerialNumberPtr psn;
  ValNodePtr        ppr;
  PubdescPtr        pdp = NULL;

  if (blocktype != REFERENCE_BLOCK || userdata == NULL) return;
  if (StringNICmp (str, "REFERENCE", 9) != 0) {
    return;
  }
  cp = str + 9;
  while (isspace (*cp)) {
    cp++;
  }
  if (!isdigit (*cp)) {
    return;
  }
  serial_number = atoi (cp);

  if (itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (entityID, NULL, itemID, 0, NULL, &fcontext);
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      bsp = GetSequenceForObject (OBJ_SEQFEAT, sfp);
    }
  } else if (itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, itemID, 0, NULL, &dcontext);
    if (sdp != NULL && sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      bsp = GetSequenceForObject (OBJ_SEQDESC, sdp);
    }
  }
  if (pdp != NULL && bsp != NULL) {    
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = PUB_Equiv;
      vnp->data.ptrvalue = pdp->pub;
      ppr = MinimizePub (vnp);
      ValNodeFree (vnp);
    }
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = PUB_Equiv;
      vnp->data.ptrvalue = ppr;

      psn = PubSerialNumberNew ();
      psn->bsp = bsp;
      psn->serial_number = serial_number;
      psn->min_pub = vnp;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, psn);
    }
  }
}


NLM_EXTERN ValNodePtr GetCitListsForSeqEntry (SeqEntryPtr sep)
{
  XtraBlock       xtra;
  ValNodePtr      head = NULL;
  ErrSev          level;
  Boolean         okay;
  SeqEntryPtr     oldscope;
  Uint2           entityID;

  if (sep == NULL) return NULL;
  
  MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
  xtra.ffwrite = CaptureRefBlockSerialNumbers;
  xtra.userdata = (Pointer) &head;
  level = ErrSetMessageLevel (SEV_MAX);
  oldscope = SeqEntrySetScope (sep);
  okay = SeqEntryToGnbk (sep, NULL, GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE,
                         SHOW_CONTIG_FEATURES, 0, 0, &xtra, NULL);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  SeqEntrySetScope (oldscope);
  ErrSetMessageLevel (level);
  return head;
}


NLM_EXTERN Int4 GetCitationNumberForMinPub (BioseqPtr bsp, ValNodePtr min_pub, ValNodePtr pub_list)
{
  Int4 rval = -1;
  PubSerialNumberPtr psn;
  ValNodePtr vnp, tmp;

  if (bsp == NULL || min_pub == NULL || pub_list == NULL) {
    return -1;
  }

  tmp = ValNodeNew (NULL);
  tmp->choice = PUB_Equiv;
  tmp->data.ptrvalue = min_pub;

  for (vnp = pub_list; vnp != NULL && rval == -1; vnp = vnp->next) {
    psn = (PubSerialNumberPtr) vnp->data.ptrvalue;
    if (psn->bsp == bsp) {
      if (PubLabelMatch (tmp, psn->min_pub) == 0) {
        rval = psn->serial_number;
      }
    }
  }

  tmp = ValNodeFree (tmp);

  return rval;
}


NLM_EXTERN ValNodePtr GetMinPubForCitationNumber (BioseqPtr bsp, Int4 number, ValNodePtr pub_list)
{
  ValNodePtr rval = NULL;
  PubSerialNumberPtr psn;
  ValNodePtr vnp;

  if (bsp == NULL || number < 0 || pub_list == NULL) {
    return NULL;
  }

  for (vnp = pub_list; vnp != NULL && rval == NULL; vnp = vnp->next) {
    psn = (PubSerialNumberPtr) vnp->data.ptrvalue;
    if (psn->bsp == bsp && psn->serial_number == number) {
      rval = psn->min_pub;
    }
  }

  return rval;
}


/* 
 * Some batch operations will be faster if information about the entire record is collected once
 * and reused.  The BatchExtra structure is where such data belongs.
 */
NLM_EXTERN BatchExtraPtr BatchExtraNew ()
{
  BatchExtraPtr b;

  b = (BatchExtraPtr) MemNew (sizeof (BatchExtraData));
  b->cit_list = NULL;

  return b;
}


NLM_EXTERN BatchExtraPtr BatchExtraFree (BatchExtraPtr b)
{
  if (b != NULL) {
    b->cit_list = PubSerialNumberListFree (b->cit_list);
    
    b = MemFree (b);
  }
  return b;
}


static Boolean IsCitationField (FieldTypePtr field)
{
  FeatureFieldPtr feature_field;

  if (field != NULL
      && field->choice == FieldType_feature_field
      && (feature_field = field->data.ptrvalue) != NULL
      && feature_field->field != NULL
      && ((feature_field->field->choice == FeatQualChoice_legal_qual 
           && feature_field->field->data.intvalue == Feat_qual_legal_citation)
          || (feature_field->field->choice == FeatQualChoice_illegal_qual 
          && DoesStringMatchConstraint ("citation", feature_field->field->data.ptrvalue)))) {
    return TRUE;
  } else {
    return FALSE;
  }

}


static void InitBatchExtraForField (BatchExtraPtr batch_extra, FieldTypePtr field, SeqEntryPtr sep)
{
  if (batch_extra == NULL) {
    return;
  }
  /* only need to collect citations if citation is in the list of applicable fields */
  if (IsCitationField (field)) {
    ValNodeLink (&(batch_extra->cit_list), GetCitListsForSeqEntry (sep));
  }
}


static void InitBatchExtraForAECRAction (BatchExtraPtr batch_extra, AECRActionPtr action, SeqEntryPtr sep)
{
  ValNodePtr field_list, field;

  if (batch_extra == NULL || action == NULL) {
    return;
  }

  field_list = GetFieldTypeListFromAECRAction (action);
  for (field = field_list; field != NULL; field = field->next) {
    InitBatchExtraForField (batch_extra, field, sep);
  }
  field_list = FieldTypeListFree (field_list);
}


NLM_EXTERN int LIBCALLBACK SortVnpByObject (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  CharPtr     str1, str2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice < vnp2->choice) {
        rval = -1;
      } else if (vnp1->choice > vnp2->choice) {
        rval = 1;
      } else {
        str1 = GetDiscrepancyItemText (vnp1);
        str2 = GetDiscrepancyItemText (vnp2);
        rval = StringCmp (str1, str1);
        str1 = MemFree (str1);
        str2 = MemFree (str2);
      }
    }
  }

  return rval;
}


static ValNodePtr BioseqListForObjectList (ValNodePtr object_list)
{
  ValNodePtr vnp, bsp_list = NULL;
  BioseqPtr  bsp;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    bsp = GetSequenceForObject (vnp->choice, vnp->data.ptrvalue);
    if (bsp != NULL) {
      ValNodeAddPointer (&bsp_list, OBJ_BIOSEQ, bsp);
    }
  }
  bsp_list = ValNodeSort (bsp_list, SortVnpByObject);
  ValNodeUnique (&bsp_list, SortVnpByObject, ValNodeFree);
  return bsp_list;
}


static void InitBatchExtraForAECRActionAndObjectList (BatchExtraPtr batch_extra, AECRActionPtr action, ValNodePtr object_list)
{
  ValNodePtr field_list, field;
  ValNodePtr bsp_list = NULL, vnp;
  SeqEntryPtr sep;

  if (batch_extra == NULL || action == NULL) {
    return;
  }

  field_list = GetFieldTypeListFromAECRAction (action);
  bsp_list = BioseqListForObjectList (object_list);
  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    sep = SeqMgrGetSeqEntryForData (vnp->data.ptrvalue);
    for (field = field_list; field != NULL; field = field->next) {
      InitBatchExtraForField (batch_extra, field, sep);
    }
  }
  bsp_list = ValNodeFree (bsp_list);

  field_list = FieldTypeListFree (field_list);

}





NLM_EXTERN CharPtr GetFieldValueForObjectEx (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  CharPtr str = NULL;
  FeatureFieldPtr feature_field;
  SeqDescrPtr     sdp;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      str = GetSourceQualFromBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        str = GetQualFromFeatureEx ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp, batch_extra);
      }
      break;
    case FieldType_cds_gene_prot :
      if (choice == 0) {
        str = GetFieldValueFromCGPSet ((CGPSetPtr) data, field->data.intvalue, scp);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        str = GetQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        str = GetSequenceQualFromBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
    case FieldType_pub :
      str = GetPubFieldFromObject (choice, data, field->data.intvalue, scp);
      break;
    case FieldType_rna_field :
      if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromRnaQual (field->data.ptrvalue);
        str = GetQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_struc_comment_field:
      if (choice == OBJ_SEQDESC && data != NULL) {
        sdp = (SeqDescrPtr) data;
        if (sdp != NULL && sdp->choice == Seq_descr_user) {
          str = GetStructuredCommentFieldFromUserObject (sdp->data.ptrvalue, field->data.ptrvalue, scp);
        }
      }
      break;
    case FieldType_misc:
      if (choice == OBJ_BIOSEQ && field->data.intvalue == Misc_field_genome_project_id) {
        str = GetGenomeProjectIdFromBioseq ((BioseqPtr) data, scp);
      } else if (choice == OBJ_SEQDESC && field->data.intvalue == Misc_field_comment_descriptor) {
        sdp = (SeqDescrPtr) data;
        if (sdp != NULL && sdp->choice == Seq_descr_comment && !StringHasNoText (sdp->data.ptrvalue)) {
          str = StringSave (sdp->data.ptrvalue);
        }
      }
      break;
  }
  return str;
}


NLM_EXTERN CharPtr GetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp)
{
  return GetFieldValueForObjectEx (choice, data, field, scp, NULL);
}

static Boolean RemoveFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  FeatureFieldPtr feature_field;
  SeqDescrPtr     sdp;
  ObjValNodePtr   ovp;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      rval = RemoveSourceQualFromBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        rval = RemoveQualFromFeature ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp);
      }
      break;
    case FieldType_cds_gene_prot:
      if (choice == 0) {
        rval = RemoveFieldValueFromCGPSet ((CGPSetPtr) data, field->data.intvalue, scp);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        rval = RemoveQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        rval = RemoveSequenceQualFromBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
    case FieldType_pub :
      rval = RemovePubFieldFromObject (choice, data, field->data.intvalue, scp);
      break;
    case FieldType_rna_field :
      if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromRnaQual (field->data.ptrvalue);
        rval = RemoveQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_struc_comment_field:
      if (choice == OBJ_SEQDESC && data != NULL) {
        sdp = (SeqDescrPtr) data;
        if (sdp != NULL && sdp->choice == Seq_descr_user) {
          rval = RemoveStructuredCommentFieldFromUserObject (sdp->data.ptrvalue, field->data.ptrvalue, scp);
        }
      }
      break;
    case FieldType_misc:
      if (choice == OBJ_BIOSEQ && field->data.intvalue == Misc_field_genome_project_id) {
        rval = RemoveGenomeProjectIdFromBioseq ((BioseqPtr) data, scp);
      } else if (choice == OBJ_SEQDESC && field->data.intvalue == Misc_field_comment_descriptor) {
        sdp = (SeqDescrPtr) data;
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
      break;
  }
  return rval;
}


NLM_EXTERN Boolean SetFieldValueForObjectEx (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text, BatchExtraPtr batch_extra)
{
  Boolean rval = FALSE;
  FeatureFieldPtr feature_field;
  SeqDescrPtr     sdp;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      rval = SetSourceQualInBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp, value, existing_text);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        rval = SetQualOnFeatureEx ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp, value, existing_text, batch_extra);
      }
      break;
    case FieldType_cds_gene_prot:
      if (choice == 0) {
        rval = SetFieldValueInCGPSet ((CGPSetPtr) data, field->data.intvalue, scp, value, existing_text);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        rval = SetQualOnFeatureEx ((SeqFeatPtr) data, feature_field, scp, value, existing_text, batch_extra);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field:
      if (choice == OBJ_BIOSEQ) {
        rval = SetSequenceQualOnBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
    case FieldType_pub :
      rval = SetPubFieldOnObject (choice, data, field->data.intvalue, scp, value, existing_text);
      break;
    case FieldType_rna_field :
      if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromRnaQual (field->data.ptrvalue);
        rval = SetQualOnFeatureEx ((SeqFeatPtr) data, feature_field, scp, value, existing_text, batch_extra);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_struc_comment_field:
      if (choice == OBJ_SEQDESC && data != NULL) {
        sdp = (SeqDescrPtr) data;
        if (sdp != NULL && sdp->choice == Seq_descr_user) {
          rval = SetStructuredCommentFieldOnUserObject (sdp->data.ptrvalue, field->data.ptrvalue, scp, value, existing_text);
        }
      }
      break;
    case FieldType_misc:
      if (choice == OBJ_BIOSEQ && field->data.intvalue == Misc_field_genome_project_id) {
        rval = SetGenomeProjectIdOnBioseq ((BioseqPtr) data, scp, value, existing_text);
      } else if (choice == OBJ_SEQDESC && field->data.intvalue == Misc_field_comment_descriptor) {
        sdp = (SeqDescrPtr) data;
        rval = SetCommentDescriptor (sdp, scp, value, existing_text);
      }
      break;
  }
  return rval;
}


NLM_EXTERN Boolean SetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  return SetFieldValueForObjectEx (choice, data, field, scp, value, existing_text, NULL);
}


NLM_EXTERN ValNodePtr GetFieldTypeListFromAECRAction (AECRActionPtr action)
{
  ValNodePtr field_list = NULL;
  ApplyActionPtr apply;
  EditActionPtr  edit;
  ConvertActionPtr convert;
  CopyActionPtr    copy;
  SwapActionPtr    swap;
  RemoveActionPtr  remove;
  AECRParseActionPtr parse;

  if (action == NULL) {
    return NULL;
  }

  /* todo - add fields from constraints ? */

  /* get fields from action */
  if (action->action != NULL) {
    switch (action->action->choice) {
      case ActionChoice_apply:
        apply = (ApplyActionPtr) action->action->data.ptrvalue;
        if (apply != NULL) {
          ValNodeLink (&field_list, FieldTypeCopy (apply->field));
        }
        break;
      case ActionChoice_edit:
        edit = (EditActionPtr) action->action->data.ptrvalue;
        if (edit != NULL) {
          ValNodeLink (&field_list, FieldTypeCopy (edit->field));
        }
        break;
      case ActionChoice_convert:
        convert = (ConvertActionPtr) action->action->data.ptrvalue;
        if (convert != NULL) {
          ValNodeLink (&field_list, GetFromFieldFromFieldPair (convert->fields));
          ValNodeLink (&field_list, GetToFieldFromFieldPair (convert->fields));
        }
        break;
      case ActionChoice_copy:
        copy = (CopyActionPtr) action->action->data.ptrvalue;
        if (copy != NULL) {
          ValNodeLink (&field_list, GetFromFieldFromFieldPair (copy->fields));
          ValNodeLink (&field_list, GetToFieldFromFieldPair (copy->fields));
        }
        break;
      case ActionChoice_swap:
        swap = (SwapActionPtr) action->action->data.ptrvalue;
        if (swap != NULL) {
          ValNodeLink (&field_list, GetFromFieldFromFieldPair (swap->fields));
          ValNodeLink (&field_list, GetToFieldFromFieldPair (swap->fields));
        }
        break;
      case ActionChoice_remove:
        remove = (RemoveActionPtr) action->action->data.ptrvalue;
        if (remove != NULL) {
          ValNodeLink (&field_list, FieldTypeCopy (remove->field));
        }
        break;
      case ActionChoice_parse:
        parse = (AECRParseActionPtr) action->action->data.ptrvalue;
        if (parse != NULL) {
          ValNodeLink (&field_list, GetFromFieldFromFieldPair (parse->fields));
          ValNodeLink (&field_list, GetToFieldFromFieldPair (parse->fields));
        }
        break;
    }
  }
  return field_list;
}


NLM_EXTERN Boolean AreAECRActionFieldsEqual (AECRActionPtr action1, AECRActionPtr action2)
{
  ApplyActionPtr a1, a2;
  EditActionPtr  e1, e2;
  ConvertActionPtr v1, v2;
  CopyActionPtr c1, c2;
  SwapActionPtr s1, s2;
  RemoveActionPtr r1, r2;
  AECRParseActionPtr p1, p2;
  FieldTypePtr       field1, field2;
  Boolean            rval = FALSE;

  if (action1 == NULL && action2 == NULL) {
    return TRUE;
  } else if (action1 == NULL || action2 == NULL) {
    return FALSE;
  } else if (action1->action == NULL && action2->action == NULL) {
    return TRUE;
  } else if (action1->action == NULL || action2->action == NULL) {
    return FALSE;
  } else if (action1->action->choice != action2->action->choice) {
    return FALSE;
  } else if (action1->action->data.ptrvalue == NULL && action2->action->data.ptrvalue == NULL) {
    return TRUE;
  } else if (action1->action->data.ptrvalue == NULL || action2->action->data.ptrvalue == NULL) { 
    return FALSE;
  } else {
    switch (action1->action->choice) {
      case ActionChoice_apply:
        a1 = (ApplyActionPtr) action1->action->data.ptrvalue;
        a2 = (ApplyActionPtr) action2->action->data.ptrvalue;
        rval = DoFieldTypesMatch (a1->field, a2->field);
        break;
      case ActionChoice_edit:
        e1 = (EditActionPtr) action1->action->data.ptrvalue;
        e2 = (EditActionPtr) action2->action->data.ptrvalue;
        rval = DoFieldTypesMatch (e1->field, e2->field);
        break;
      case ActionChoice_convert:
        v1 = (ConvertActionPtr) action1->action->data.ptrvalue;
        v2 = (ConvertActionPtr) action2->action->data.ptrvalue;
        field1 = GetFromFieldFromFieldPair (v1->fields);
        field2 = GetFromFieldFromFieldPair (v2->fields);
        rval = DoFieldTypesMatch (field1, field2);
        if (rval) {
          field1 = FieldTypeFree (field1);
          field2 = FieldTypeFree (field2);
          field1 = GetToFieldFromFieldPair (v1->fields);
          field2 = GetToFieldFromFieldPair (v2->fields);
          rval = DoFieldTypesMatch (field1, field2);
        }
        field1 = FieldTypeFree (field1);
        field2 = FieldTypeFree (field2);
        if (rval) {
          if ((v1->keep_original && !v2->keep_original)
            || (!v1->keep_original && v2->keep_original)) {
            rval = FALSE;
          }
        }
        break;
      case ActionChoice_copy:
        c1 = (CopyActionPtr) action1->action->data.ptrvalue;
        c2 = (CopyActionPtr) action2->action->data.ptrvalue;
        field1 = GetFromFieldFromFieldPair (c1->fields);
        field2 = GetFromFieldFromFieldPair (c2->fields);
        rval = DoFieldTypesMatch (field1, field2);
        if (rval) {
          field1 = FieldTypeFree (field1);
          field2 = FieldTypeFree (field2);
          field1 = GetToFieldFromFieldPair (c1->fields);
          field2 = GetToFieldFromFieldPair (c2->fields);
          rval = DoFieldTypesMatch (field1, field2);
        }
        field1 = FieldTypeFree (field1);
        field2 = FieldTypeFree (field2);
        break;
      case ActionChoice_swap:
        s1 = (SwapActionPtr) action1->action->data.ptrvalue;
        s2 = (SwapActionPtr) action2->action->data.ptrvalue;
        field1 = GetFromFieldFromFieldPair (s1->fields);
        field2 = GetFromFieldFromFieldPair (s2->fields);
        rval = DoFieldTypesMatch (field1, field2);
        if (rval) {
          field1 = FieldTypeFree (field1);
          field2 = FieldTypeFree (field2);
          field1 = GetToFieldFromFieldPair (s1->fields);
          field2 = GetToFieldFromFieldPair (s2->fields);
          rval = DoFieldTypesMatch (field1, field2);
        }
        field1 = FieldTypeFree (field1);
        field2 = FieldTypeFree (field2);
        break;
      case ActionChoice_remove:
        r1 = (RemoveActionPtr) action1->action->data.ptrvalue;
        r2 = (RemoveActionPtr) action2->action->data.ptrvalue;
        rval = DoFieldTypesMatch (r1->field, r2->field);
        break;
      case ActionChoice_parse:
        p1 = (AECRParseActionPtr) action1->action->data.ptrvalue;
        p2 = (AECRParseActionPtr) action2->action->data.ptrvalue;
        field1 = GetFromFieldFromFieldPair (p1->fields);
        field2 = GetFromFieldFromFieldPair (p2->fields);
        rval = DoFieldTypesMatch (field1, field2);
        if (rval) {
          field1 = FieldTypeFree (field1);
          field2 = FieldTypeFree (field2);
          field1 = GetToFieldFromFieldPair (p1->fields);
          field2 = GetToFieldFromFieldPair (p2->fields);
          rval = DoFieldTypesMatch (field1, field2);
        }
        field1 = FieldTypeFree (field1);
        field2 = FieldTypeFree (field2);
        break;
    }
  }
  return rval;
}


static Boolean IsNonTextSourceQualPresent (BioSourcePtr biop, Int4 srcqual)
{
  Int4 orgmod_subtype, subsrc_subtype, subfield;
  OrgModPtr mod;
  SubSourcePtr ssp;
  Boolean      rval = FALSE;

  if (biop == NULL) return FALSE;

  orgmod_subtype = GetOrgModQualFromSrcQual (srcqual, &subfield);
  if (orgmod_subtype == -1) {
    subsrc_subtype = GetSubSrcQualFromSrcQual (srcqual, &subfield);
    for (ssp = biop->subtype; ssp != NULL && !rval; ssp = ssp->next) {
      if (ssp->subtype == subsrc_subtype) {
        rval = TRUE;
      }
    }
  } else {
    if (biop->org != NULL && biop->org->orgname != NULL) {
      for (mod = biop->org->orgname->mod; mod != NULL && !rval; mod = mod->next) {
        if (mod->subtype == orgmod_subtype) {
          rval = TRUE;
        }
      }
    }
  }
  return rval;
}


static Boolean IsSourceQualPresent (BioSourcePtr biop, SourceQualChoicePtr scp)
{
  Boolean rval = FALSE;
  CharPtr   str;

  if (biop == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  switch (scp->choice) {
    case SourceQualChoice_textqual:
      if (IsNonTextSourceQual (scp->data.intvalue)) {
        rval = IsNonTextSourceQualPresent (biop, scp->data.intvalue);
      } else {
        str = GetSourceQualFromBioSource (biop, scp, NULL);
        if (!StringHasNoText (str)) {
          rval = TRUE;
        }
        str = MemFree (str);
      }
      break;
    case SourceQualChoice_location:
      if (biop->genome != 0) {
        rval = TRUE;
      }
      break;
    case SourceQualChoice_origin:
      if (biop->origin != 0) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


typedef struct objecthasstring
{
  StringConstraintPtr scp;
  Boolean             found;
} ObjectHasStringData, PNTR ObjectHasStringPtr;


static void LIBCALLBACK AsnWriteConstraintCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr            pchSource;
  ObjectHasStringPtr ohsp;

  ohsp = (ObjectHasStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) 
  {
	  pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	  ohsp->found |= DoesSingleStringMatchConstraint (pchSource, ohsp->scp);
  }
}


static Boolean DoesObjectMatchStringConstraint (Uint1 choice, Pointer data, StringConstraintPtr scp)

{
  ObjMgrPtr         omp;
  ObjMgrTypePtr     omtp;
  AsnIoPtr          aip;
  AsnExpOptPtr      aeop;
  ObjectHasStringData ohsd;
  SeqFeatPtr          sfp, prot;
  SeqMgrFeatContext   fcontext;
  CharPtr             search_txt;
  CGPSetPtr           c;
  ValNodePtr          vnp;
  Boolean             all_match = TRUE, any_match = FALSE, rval;
  BioseqPtr           protbsp;

  if (data == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  if (choice == 0) {
    /* CDS-Gene-Prot set */
    c = (CGPSetPtr) data;
    for (vnp = c->gene_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->cds_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->mrna_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->prot_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    if (scp->not_present) {
      rval = all_match;
    } else {
      rval = any_match;
    }        
  } else {
    omp = ObjMgrGet ();
    omtp = ObjMgrTypeFind (omp, choice, NULL, NULL);
    if (omtp == NULL) return FALSE;
    aip = AsnIoNullOpen ();
    aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteConstraintCallBack);
    ohsd.found = FALSE;
    ohsd.scp = scp;
    if (aeop != NULL) {
      aeop->user_data = (Pointer) &ohsd;
    }
    
    (omtp->asnwrite) (data, aip, NULL);
    
    if (!ohsd.found && omtp->datatype == OBJ_SEQFEAT)
    {
      sfp = (SeqFeatPtr) data;
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        protbsp = BioseqFindFromSeqLoc (sfp->product);
        prot = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &fcontext);
        if (prot != NULL) {
          (omtp->asnwrite) (prot, aip, NULL);
        }
      } else {
        if (SeqMgrFeaturesAreIndexed(sfp->idx.entityID) == 0) {
          SeqMgrIndexFeatures (sfp->idx.entityID, NULL);
        }
        if (sfp->idx.subtype == FEATDEF_tRNA) {
          sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
          ohsd.found = DoesSingleStringMatchConstraint (fcontext.label, ohsd.scp);
          if (!ohsd.found && sfp != NULL && sfp->idx.subtype == FEATDEF_tRNA)
          {
            search_txt = (CharPtr) MemNew ((StringLen (fcontext.label) + 6) * sizeof (Char));
            if (search_txt != NULL)
            {
              sprintf (search_txt, "tRNA-%s", fcontext.label);
              ohsd.found = DoesSingleStringMatchConstraint (search_txt, ohsd.scp);
              search_txt = MemFree (search_txt);
            }
          }
        }
      }
    }
    AsnIoClose (aip);
    if (scp->not_present) {
      rval = !ohsd.found;
    } else {
      rval = ohsd.found;
    }
  }
  return rval;
}


NLM_EXTERN Boolean IsSourceConstraintEmpty (SourceConstraintPtr scp)
{
  if (scp == NULL) return TRUE;

  if (scp->field1 == NULL
      && scp->field2 == NULL
      && IsStringConstraintEmpty(scp->constraint)) {
    return TRUE;
  } else {
    return FALSE;
  }
}

NLM_EXTERN Boolean DoesBiosourceMatchConstraint (BioSourcePtr biop, SourceConstraintPtr scp)
{
  Boolean rval = FALSE;
  CharPtr str1, str2;
  ValNode vn;

  if (biop == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  if (IsStringConstraintEmpty(scp->constraint)) {
    /* looking for qual present */
    if (scp->field1 != NULL && scp->field2 == NULL) {
      rval = IsSourceQualPresent (biop, scp->field1);
    } else if (scp->field2 != NULL && scp->field1 == NULL) {
      rval = IsSourceQualPresent (biop, scp->field2);
    /* looking for quals to match */
    } else if (scp->field1 != NULL && scp->field2 != NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, NULL);
      str2 = GetSourceQualFromBioSource (biop, scp->field2, NULL);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* nothing specified, automatic match */
      rval = TRUE;
    }
  } else {
    if (scp->field1 != NULL && scp->field2 == NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, scp->constraint);
      if (str1 == NULL) {
        if (scp->constraint->not_present) {
          str1 = GetSourceQualFromBioSource (biop, scp->field1, NULL);
          if (str1 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str1)) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
    } else if (scp->field2 != NULL && scp->field1 == NULL) {
      str2 = GetSourceQualFromBioSource (biop, scp->field2, scp->constraint);
      if (str2 == NULL) {
        if (scp->constraint->not_present) {
          str2 = GetSourceQualFromBioSource (biop, scp->field2, NULL);
          if (str2 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str2)) {
        rval = TRUE;
      }
      str2 = MemFree (str2);
    } else if (scp->field1 != NULL && scp->field2 != NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, scp->constraint);
      str2 = GetSourceQualFromBioSource (biop, scp->field2, scp->constraint);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* generic string constraint */
      vn.choice = Seq_descr_source;
      vn.next = NULL;
      vn.extended = 0;
      vn.data.ptrvalue = biop;
      rval = DoesObjectMatchStringConstraint (OBJ_SEQDESC, &vn, scp->constraint);
    }
  }
  return rval;
}


static Boolean DoesCGPSetMatchPseudoConstraint (CGPSetPtr c, CDSGeneProtPseudoConstraintPtr constraint)
{
  Boolean    any_pseudo = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Boolean    rval = FALSE;

  if (c == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  switch (constraint->feature) {
    case CDSGeneProt_feature_type_constraint_gene :
      for (vnp = c->gene_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_mRNA :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_cds :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_prot :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo && sfp->idx.subtype == FEATDEF_PROT) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_mat_peptide :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo && sfp->idx.subtype == FEATDEF_mat_peptide_aa) {
          any_pseudo = TRUE;
        }
      }
      break;
  }

  if ((any_pseudo && constraint->is_pseudo)
      || (!any_pseudo && !constraint->is_pseudo)) {
    rval = TRUE;
  }
  return rval;
}


NLM_EXTERN Boolean IsCDSGeneProtQualConstraintEmpty (CDSGeneProtQualConstraintPtr constraint)
{
  if (constraint == NULL) return TRUE;
  if (constraint->field1 == NULL && constraint->field2 == NULL && IsStringConstraintEmpty (constraint->constraint)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean DoesCGPSetMatchQualConstraint (CGPSetPtr c, CDSGeneProtQualConstraintPtr constraint)
{
  Boolean rval = FALSE;
  CharPtr str, str1, str2;

  if (c == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  if (IsStringConstraintEmpty (constraint->constraint)) {
    /* looking for qual present */
    if (constraint->field1 != NULL && constraint->field2 == NULL) {
      str = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, NULL);
      if (str != NULL) {
        rval = TRUE;
        str = MemFree (str);
      }
    } else if (constraint->field2 != NULL && constraint->field1 == NULL) {
      str = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, NULL);
      if (str == NULL) {
        rval = FALSE;
      } else {
        str = MemFree (str);
      }
    /* looking for quals to match */
    } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, NULL);
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, NULL);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* nothing specified, automatic match */
      rval = TRUE;
    }
  } else {
    if (constraint->field1 != NULL && constraint->field2 == NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
      if (str1 == NULL) {
        if (constraint->constraint->not_present) {
          str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
          if (str1 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str1)) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
    } else if (constraint->field2 != NULL && constraint->field1 == NULL) {
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
      if (str2 == NULL) {
        if (constraint->constraint->not_present) {
          str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
          if (str2 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str2)) {
        rval = TRUE;
      }
      str2 = MemFree (str2);
    } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* generic string constraint */
      rval = DoesObjectMatchStringConstraint (0, c, constraint->constraint);
    }
  }
  return rval;
}


static Boolean DoesSequenceHaveFeatureWithQualPresent (BioseqPtr bsp, FeatureFieldPtr feature_field, StringConstraintPtr scp)
{
  Boolean           rval = FALSE;
  SeqFeatPtr        sfp, sfp_p;
  SeqMgrFeatContext context1, context2;
  Int4              featdef;
  Uint1             seqfeattype;
  CharPtr           str;
  BioseqPtr         prot_bsp;

  if (bsp == NULL) {
    return FALSE;
  } else if (feature_field == NULL) {
    return TRUE;
  }
  featdef = GetFeatdefFromFeatureType(feature_field->type);
  seqfeattype = FindFeatFromFeatDefType (featdef);
  if (seqfeattype == SEQFEAT_PROT) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &context1);
        sfp != NULL && !rval;
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &context1)) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      for (sfp_p = SeqMgrGetNextFeature (prot_bsp, NULL, 0, featdef, &context2);
            sfp_p != NULL && !rval;
            sfp_p = SeqMgrGetNextFeature (prot_bsp, sfp_p, 0, featdef, &context2)) {
        str = GetQualFromFeature (sfp_p, feature_field, scp);
        if (str == NULL && scp != NULL) {
          if (scp->not_present) {
            str = GetQualFromFeature (sfp_p, feature_field, NULL);
            if (str == NULL) {
              rval = TRUE;
            }
          }
        } else if (!StringHasNoText (str)) {
          rval = TRUE;
        }
        str = MemFree (str);
      }
    }
  } else {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &context1);
        sfp != NULL && !rval;
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &context1)) {
      str = GetQualFromFeature (sfp, feature_field, scp);
      if (str == NULL && scp != NULL) {
        if (scp->not_present) {
          str = GetQualFromFeature (sfp, feature_field, NULL);
          if (str == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str)) {
        rval = TRUE;
      }
      str = MemFree (str);
    }
  }
  return rval;
}


static Boolean 
DoesSequenceHaveFeatureWithMatchingQuals 
(BioseqPtr bsp,
 CDSGeneProtConstraintFieldPtr f1,
 CDSGeneProtConstraintFieldPtr f2,
 StringConstraintPtr           scp)
{
  Int4              featdef;
  Uint1             seqfeattype;
  SeqFeatPtr        sfp, sfp_p;
  CharPtr           str, str2;
  SeqMgrFeatContext context1, context2;
  FeatureFieldPtr   feature_field1 = NULL, feature_field2 = NULL;
  CGPSetPtr         c;
  Boolean           b = FALSE;
  Boolean           rval = FALSE;
  BioseqPtr         prot_bsp;

  if (bsp == NULL || f1 == NULL || f2 == NULL) {
    return FALSE;
  }
  feature_field1 = FeatureFieldFromCDSGeneProtField(f1->data.intvalue);
  feature_field2 = FeatureFieldFromCDSGeneProtField(f2->data.intvalue);

  if (feature_field1 == NULL || feature_field2 == NULL) {
    feature_field1 = FeatureFieldFree (feature_field1);
    feature_field2 = FeatureFieldFree (feature_field2);
    return FALSE;
  }

  if (feature_field1->type == feature_field2->type) {
    featdef = GetFeatdefFromFeatureType(feature_field1->type);
    seqfeattype = FindFeatFromFeatDefType (featdef);
    if (seqfeattype == SEQFEAT_PROT) {
      for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &context1);
          sfp != NULL && !rval;
          sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &context1)) {
        prot_bsp = BioseqFindFromSeqLoc (sfp->product);
        for (sfp_p = SeqMgrGetNextFeature (prot_bsp, NULL, 0, featdef, &context2);
            sfp_p != NULL && !rval;
            sfp_p = SeqMgrGetNextFeature (prot_bsp, sfp_p, 0, featdef, &context2)) {
          str = GetQualFromFeature (sfp_p, feature_field1, scp);
          str2 = GetQualFromFeature (sfp_p, feature_field2, scp);
          if (str != NULL && str2 != NULL && StringCmp (str, str2) == 0) {              
            rval = TRUE;
          }
          str = MemFree (str);
          str2 = MemFree (str2);
        }
      }
    } else {
      for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &context1);
          sfp != NULL && !rval;
          sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &context1)) {
        str = GetQualFromFeature (sfp, feature_field1, scp);
        str2 = GetQualFromFeature (sfp, feature_field2, scp);
        if (str != NULL && str2 != NULL && StringCmp (str, str2) == 0) {              
          rval = TRUE;
        }
        str = MemFree (str);
        str2 = MemFree (str2);
      }
    }
  } else {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &context1);
        sfp != NULL && !rval;
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &context1)) {
      c = BuildCGPSetFromCodingRegion (sfp, &b);
      str = GetFieldValueFromCGPSet (c, f1->data.intvalue, scp);
      str2 = GetFieldValueFromCGPSet (c, f2->data.intvalue, scp);
      if (str != NULL && str2 != NULL && StringCmp (str, str2) == 0) {              
        rval = TRUE;
      }
      str = MemFree (str);
      str2 = MemFree (str2);
      c = CGPSetFree (c);
    }
  }
  return rval;
}


static Boolean DoesSequenceMatchCGPQualConstraint (BioseqPtr bsp, CDSGeneProtQualConstraintPtr constraint)
{
  FeatureFieldPtr feature_field;
  Boolean         rval = FALSE;

  if (bsp == NULL) {
    return FALSE;
  } else if (constraint == NULL) {
    return TRUE;
  }

  if (IsStringConstraintEmpty (constraint->constraint)) {
    /* looking for qual present */
    if ((constraint->field1 != NULL && constraint->field2 == NULL) 
        || (constraint->field2 != NULL && constraint->field1 == NULL)) {
      if (constraint->field1 != NULL) {
        feature_field = FeatureFieldFromCDSGeneProtField (constraint->field1->data.intvalue);
      } else {
        feature_field = FeatureFieldFromCDSGeneProtField (constraint->field2->data.intvalue);
      }
      if (feature_field != NULL) {
        rval = DoesSequenceHaveFeatureWithQualPresent (bsp, feature_field, NULL);
        feature_field = FeatureFieldFree (feature_field);
      }
    /* looking for quals to match */
    } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
      rval = DoesSequenceHaveFeatureWithMatchingQuals (bsp, constraint->field1, constraint->field2, NULL);
    } else {
      /* nothing specified, automatic match */
      rval = TRUE;
    }
  } else if ((constraint->field1 != NULL && constraint->field2 == NULL)
             || (constraint->field1 == NULL && constraint->field2 != NULL)) {
    /* one field must match constraint */
    if (constraint->field1 != NULL) {
      feature_field = FeatureFieldFromCDSGeneProtField (constraint->field1->data.intvalue);
    } else {
      feature_field = FeatureFieldFromCDSGeneProtField (constraint->field2->data.intvalue);
    }
    if (feature_field != NULL) {
      rval = DoesSequenceHaveFeatureWithQualPresent (bsp, feature_field, constraint->constraint);
      feature_field = FeatureFieldFree (feature_field);
    }
  } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
    /* two fields must match and match constraint */
    rval = DoesSequenceHaveFeatureWithMatchingQuals (bsp, constraint->field1, constraint->field2, constraint->constraint);
  } else {
    /* generic string constraint */
    rval = DoesObjectMatchStringConstraint (OBJ_BIOSEQ, bsp, constraint->constraint);
  }
  return rval;
}


static Boolean DoesSequenceInSetMatchCGPQualConstraint (BioseqSetPtr bssp, CDSGeneProtQualConstraintPtr constraint)
{
  Boolean       rval = FALSE;
  SeqEntryPtr   sep;

  if (bssp == NULL) return FALSE;
  if (constraint == NULL) return TRUE;
  
  for (sep = bssp->seq_set; sep != NULL && !rval; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      rval = DoesSequenceMatchCGPQualConstraint ((BioseqPtr) sep->data.ptrvalue, constraint);
    } else if (IS_Bioseq_set (sep)) {
      rval = DoesSequenceInSetMatchCGPQualConstraint ((BioseqSetPtr) sep->data.ptrvalue, constraint);
    }
  }
  return rval;
}


static Boolean DoesSeqDescMatchCGPQualConstraint (SeqDescrPtr sdp, CDSGeneProtQualConstraintPtr constraint)
{
  Boolean rval = FALSE;
  BioseqPtr bsp;
  ObjValNodePtr ovp;

  if (sdp == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  bsp = GetSequenceForObject (OBJ_SEQDESC, sdp);
  if (bsp == NULL) {
    if (sdp->extended) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQSET && ovp->idx.parentptr != NULL) {
        rval = DoesSequenceInSetMatchCGPQualConstraint ((BioseqSetPtr) ovp->idx.parentptr, constraint);
      }
    }
  } else {
    rval = DoesSequenceMatchCGPQualConstraint (bsp, constraint);
  }

  return rval;
}


static Boolean DoesFeatureMatchCGPQualConstraint (SeqFeatPtr sfp, CDSGeneProtQualConstraintPtr constraint)
{
  CGPSetPtr c = NULL;
  Boolean   b = FALSE;
  SeqMgrFeatContext context;
  Boolean           rval = FALSE;
  FeatureFieldPtr   ff;
  SeqFeatPtr        cds;
  CharPtr           str1 = NULL, str2 = NULL;

  if (sfp == NULL) {
    return FALSE;
  } else if (constraint == NULL) {
    return TRUE;
  }

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    c = BuildCGPSetFromCodingRegion (sfp, &b);
  } else if (sfp->data.choice == SEQFEAT_PROT) {
    cds = SeqMgrGetCDSgivenProduct (BioseqFindFromSeqLoc (sfp->location), &context);
    c = BuildCGPSetFromCodingRegion (cds, &b);
  } else if (sfp->data.choice == SEQFEAT_GENE) {
    c = BuildCGPSetFromGene (sfp);
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    c = BuildCGPSetFrommRNA (sfp);
  }

  rval = DoesCGPSetMatchQualConstraint (c, constraint);
  if (rval && sfp->idx.subtype == FEATDEF_mat_peptide_aa) {
    if (constraint->field1 != NULL) {
      if (IsCDSGeneProtFieldMatPeptideRelated (constraint->field1->data.intvalue)) {
        ff = FeatureFieldFromCDSGeneProtField (constraint->field1->data.intvalue);
        str1 = GetQualFromFeature (sfp, ff, constraint->constraint);
        ff = FeatureFieldFree (ff);
      } else {
        str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
      }
      if (str1 == NULL) {
        rval = FALSE;
      }
    }
    if (constraint->field2 != NULL) {
      if (IsCDSGeneProtFieldMatPeptideRelated (constraint->field2->data.intvalue)) {
        ff = FeatureFieldFromCDSGeneProtField (constraint->field2->data.intvalue);
        str2 = GetQualFromFeature (sfp, ff, constraint->constraint);
        ff = FeatureFieldFree (ff);
      } else {
        str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
      }
      if (str2 == NULL) {
        rval = FALSE;
      }
    }
    if (rval && constraint->field1 != NULL && constraint->field2 != NULL && StringCmp (str1, str2) != 0) {
      rval = FALSE;
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  c = CGPSetFree (c);
  return rval;
}


NLM_EXTERN Boolean IsSequenceConstraintEmpty (SequenceConstraintPtr constraint)
{
  if (constraint == NULL) return TRUE;
  if (constraint->seqtype != NULL && constraint->seqtype->choice != SequenceConstraintMolTypeConstraint_any) return FALSE;
  if (constraint->feature != Feature_type_any) return FALSE;
  if (!IsStringConstraintEmpty (constraint->id)) return FALSE;
  return TRUE;
}


static Boolean DoesSeqIDListMeetStringConstraint (SeqIdPtr sip, StringConstraintPtr string_constraint)
{
  Char       id [41];
  CharPtr    cp, cp_dst;
  SeqIdPtr   tmp;
  Boolean    match, changed;

  if (sip == NULL) 
  {
    return FALSE;
  }
  if (string_constraint == NULL)
  {
    return TRUE;
  }
  
  while (sip != NULL)
  {
    /* temporary disconnect ID from list */
    tmp = sip->next;
    sip->next = NULL;
    id [0] = '\0';
    SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
    match = DoesSingleStringMatchConstraint (id, string_constraint);
    if (!match) 
    {
      changed = FALSE;
      /* remove terminating pipe character */
      if (id[StringLen(id) - 1] == '|') 
      {
        id[StringLen(id) - 1] = 0;
        changed = TRUE;
      }
      /* remove leading pipe identifier */
      cp = StringChr (id, '|');
      if (cp != NULL)
      {
        changed = TRUE;
        cp++;
        cp_dst = id;
        while (*cp != 0) 
        {
          *cp_dst = *cp;
          cp_dst++;
          cp++;
        }
        *cp_dst = 0;
      }  
      if (changed) 
      {
        match = DoesSingleStringMatchConstraint (id, string_constraint);
      }

      /* if search text doesn't have ., try ID without version */
      if (!match && StringChr (string_constraint->match_text, '.') == NULL) 
      {
        cp = StringChr (id, '.');
        if (cp != NULL) 
        {
          *cp = 0;
          match = DoesSingleStringMatchConstraint (id, string_constraint);
        }
      }       
    }
    sip->next = tmp;

    if (match)
    {
      if (string_constraint->not_present)
      {
        return FALSE;
      }
      else
      {
        return TRUE;
      }
    }
    sip = sip->next;
  }
  if (string_constraint->not_present)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


typedef struct rnatypebiomol {
  Int4 rnatype;
  Uint1 biomol;
  CharPtr rnamolname;
} RnaTypeBiomolData, PNTR RnaTypeBiomolPtr;

static RnaTypeBiomolData rna_type_biomol[] = {
{ Sequence_constraint_rnamol_genomic , MOLECULE_TYPE_GENOMIC, "Genomic RNA" } ,
{ Sequence_constraint_rnamol_precursor_RNA , MOLECULE_TYPE_PRE_MRNA , "Precursor RNA" } ,
{ Sequence_constraint_rnamol_mRNA , MOLECULE_TYPE_MRNA , "mRNA [cDNA]" } ,
{ Sequence_constraint_rnamol_rRNA , MOLECULE_TYPE_RRNA , "Ribosomal RNA" } ,
{ Sequence_constraint_rnamol_tRNA , MOLECULE_TYPE_TRNA , "Transfer RNA" } ,
{ Sequence_constraint_rnamol_genomic_mRNA , MOLECULE_TYPE_GENOMIC_MRNA_MIX , "Genomic-mRNA" } ,
{ Sequence_constraint_rnamol_cRNA , MOLECULE_TYPE_CRNA , "cRNA" } ,
{ Sequence_constraint_rnamol_transcribed_RNA , MOLECULE_TYPE_TRANSCRIBED_RNA , "Transcribed RNA" } ,
{ Sequence_constraint_rnamol_ncRNA , MOLECULE_TYPE_NCRNA , "Non-coding  RNA" } ,
{ Sequence_constraint_rnamol_transfer_messenger_RNA , MOLECULE_TYPE_TMRNA , "Transfer-messenger RNA" } } ;

#define NUM_rna_type_biomol sizeof (rna_type_biomol) / sizeof (RnaTypeBiomolData)


NLM_EXTERN Uint1 GetBiomolForRnaType (Int4 rnatype) 
{
  Int4 i;

  for (i = 0; i <  NUM_rna_type_biomol; i++) {
    if (rna_type_biomol[i].rnatype == rnatype) {
      return rna_type_biomol[i].biomol;
    }
  }
  return 0;
}


NLM_EXTERN CharPtr GetBiomolNameForRnaType (Int4 rnatype)
{
  Int4 i;

  for (i = 0; i <  NUM_rna_type_biomol; i++) {
    if (rna_type_biomol[i].rnatype == rnatype) {
      return rna_type_biomol[i].rnamolname;
    }
  }
  return "invalid RNA type";
}

NLM_EXTERN void AddAllRNASubtypesToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  if (field_list == NULL) return;

  ValNodeAddPointer (field_list, Sequence_constraint_rnamol_any, StringSave ("Any RNA"));
  for (i = 0; i < NUM_rna_type_biomol; i++) {
    ValNodeAddPointer (field_list, rna_type_biomol[i].rnatype, StringSave (rna_type_biomol[i].rnamolname));
  }
}


static Boolean DoesSequenceMatchSequenceConstraint (BioseqPtr bsp, SequenceConstraintPtr constraint)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  MolInfoPtr mip;
  
  if (bsp == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;

  if (constraint->seqtype != NULL && constraint->seqtype->choice != SequenceConstraintMolTypeConstraint_any) {
    switch (constraint->seqtype->choice) {
      case SequenceConstraintMolTypeConstraint_nucleotide :
        if (ISA_aa (bsp->mol)) {
          return FALSE;
        }
        break;
      case SequenceConstraintMolTypeConstraint_dna :
        if (bsp->mol != Seq_mol_dna) {
          return FALSE;
        }
        break;
      case SequenceConstraintMolTypeConstraint_rna :
        if (bsp->mol != Seq_mol_rna) {
          return FALSE;
        }
        if (constraint->seqtype->data.intvalue != Sequence_constraint_rnamol_any) {
          sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
          if (sdp == NULL || sdp->data.ptrvalue == NULL || sdp->choice != Seq_descr_molinfo) {
            return FALSE;
          }
          mip = (MolInfoPtr) sdp->data.ptrvalue;
          if (GetBiomolForRnaType (constraint->seqtype->data.intvalue) != mip->biomol) {
            return FALSE;
          }
        }
        break;
      case SequenceConstraintMolTypeConstraint_protein :
        if (!ISA_aa (bsp->mol)) {
          return FALSE;
        }
        break;
    }
  }

  if (constraint->feature != Feature_type_any) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, GetFeatdefFromFeatureType (constraint->feature), &fcontext);
    if (sfp == NULL) {
      return FALSE;
    }
  }

  if (!IsStringConstraintEmpty (constraint->id) && !DoesSeqIDListMeetStringConstraint (bsp->id, constraint->id)) {
    return FALSE;
  }
  return TRUE;
}

static Boolean DoesSequenceInSetMatchSequenceConstraint (BioseqSetPtr bssp, SequenceConstraintPtr constraint)
{
  Boolean       rval = FALSE;
  SeqEntryPtr   sep;

  if (bssp == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;
  
  for (sep = bssp->seq_set; sep != NULL && !rval; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      rval = DoesSequenceMatchSequenceConstraint ((BioseqPtr) sep->data.ptrvalue, constraint);
    } else if (IS_Bioseq_set (sep)) {
      rval = DoesSequenceInSetMatchSequenceConstraint ((BioseqSetPtr) sep->data.ptrvalue, constraint);
    }
  }
  return rval;
}


static Boolean DoesObjectMatchSequenceConstraint (Uint1 choice, Pointer data, SequenceConstraintPtr constraint)
{
  BioseqPtr bsp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  Boolean       rval = FALSE;

  if (data == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;

  bsp = GetSequenceForObject (choice, data);
  if (bsp == NULL) {
    if (choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) data;
      if (sdp->extended) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET && ovp->idx.parentptr != NULL) {
          rval = DoesSequenceInSetMatchSequenceConstraint ((BioseqSetPtr) ovp->idx.parentptr, constraint);
        }
      }
    }
  } else {
    rval = DoesSequenceMatchSequenceConstraint (bsp, constraint);
  }
  return rval; 
}


NLM_EXTERN CharPtr GetPubFieldLabel (Int4 pub_field)
{
  CharPtr rval = NULL;
  switch (pub_field) {
    case Publication_field_cit:
      rval = "citation";
      break;
    case Publication_field_authors:
      rval = "authors";
      break;
    case Publication_field_journal:
      rval = "journal";
      break;
    case Publication_field_volume:
      rval = "volume";
      break;
    case Publication_field_issue:
      rval = "issue";
      break;
    case Publication_field_pages:
      rval = "pages";
      break;
    case Publication_field_date:
      rval = "date";
      break;
    case Publication_field_serial_number:
      rval = "serial number";
      break;
    case Publication_field_title:
      rval = "title";
      break;
    case Publication_field_affiliation:
      rval = "affiliation";
      break;
    case Publication_field_affil_div:
      rval = "department";
      break;
    case Publication_field_affil_city:
      rval = "city";
      break;
    case Publication_field_affil_sub:
      rval = "state";
      break;
    case Publication_field_affil_country:
      rval = "country";
      break;
    case Publication_field_affil_street:
      rval = "street";
      break;
    case Publication_field_affil_email:
      rval = "email";
      break;
    case Publication_field_affil_fax:
      rval = "fax";
      break;
    case Publication_field_affil_phone:
      rval = "phone";
      break;
    case Publication_field_affil_zipcode:
      rval = "postal code";
      break;
  }
  return rval;
}


NLM_EXTERN ValNodePtr GetPubFieldList (void)
{
  ValNodePtr         val_list = NULL;

  ValNodeAddPointer (&val_list, Publication_field_title, StringSave ("title"));
  ValNodeAddPointer (&val_list, Publication_field_authors, StringSave ("authors"));
  ValNodeAddPointer (&val_list, Publication_field_journal, StringSave ("journal"));
  ValNodeAddPointer (&val_list, Publication_field_issue, StringSave ("issue"));
  ValNodeAddPointer (&val_list, Publication_field_pages, StringSave ("pages"));
  ValNodeAddPointer (&val_list, Publication_field_serial_number, StringSave ("serial number"));
  ValNodeAddPointer (&val_list, Publication_field_date, StringSave ("date"));
  ValNodeAddPointer (&val_list, Publication_field_cit, StringSave ("citation"));
  ValNodeAddPointer (&val_list, Publication_field_affiliation, StringSave ("affiliation"));
  ValNodeAddPointer (&val_list, Publication_field_affil_div, StringSave ("department"));
  ValNodeAddPointer (&val_list, Publication_field_affil_city, StringSave ("city"));
  ValNodeAddPointer (&val_list, Publication_field_affil_sub, StringSave ("state"));
  ValNodeAddPointer (&val_list, Publication_field_affil_country, StringSave ("country"));
  ValNodeAddPointer (&val_list, Publication_field_affil_street, StringSave ("street"));
  ValNodeAddPointer (&val_list, Publication_field_affil_email, StringSave ("email"));
  ValNodeAddPointer (&val_list, Publication_field_affil_fax, StringSave ("fax"));
  ValNodeAddPointer (&val_list, Publication_field_affil_phone, StringSave ("phone"));
  ValNodeAddPointer (&val_list, Publication_field_affil_zipcode, StringSave ("postal code"));

  return val_list;
}


static ValNodePtr MakePubFieldTypeList (void)
{
  ValNodePtr field_list = NULL;

  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_title);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_authors);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_journal);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_issue);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_pages);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_serial_number);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_date);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_cit);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affiliation);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_div);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_city);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_sub);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_country);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_street);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_email);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_fax);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_phone);
  ValNodeAddInt (&field_list, FieldType_pub, Publication_field_affil_zipcode);

  return field_list;
}


NLM_EXTERN Boolean IsPublicationConstraintEmpty (PublicationConstraintPtr constraint)
{
  Boolean rval = FALSE;

  if (constraint == NULL
      || (constraint->type == Pub_type_any
          && (constraint->field == NULL 
              || IsStringConstraintEmpty (constraint->field->constraint)))) {
    rval = TRUE;
  }
  return rval;
}


NLM_EXTERN Int4 GetPubMLStatus (PubPtr the_pub)
{
  CitGenPtr  cgp;
  CitSubPtr  csp;
  CitArtPtr  cap;
  CitBookPtr cbp;
  CitJourPtr cjp;
  ImprintPtr imp = NULL;
  Int4       status = Pub_type_any;
  
  if (the_pub == NULL || the_pub->data.ptrvalue == NULL)
  {
    return Pub_type_any;
  }
  
  switch (the_pub->choice)
  {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      if (cgp->cit != NULL && StringICmp (cgp->cit, "unpublished") == 0)
      {
        status = Pub_type_unpublished;
      }
      else
      {
        status = Pub_type_published;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      status = Pub_type_submitter_block;
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      if (cap->from == 1)
      {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL)
        {
          imp = cjp->imp;
        }
      }
	  else if (cap->from == 2 || cap->from == 3) 
	  {
        cbp = (CitBookPtr) cap->fromptr;
		if (cbp != NULL) {
          imp = cbp->imp;
		}
	  }
      break;
    case PUB_Journal :
      cjp = (CitJourPtr) the_pub->data.ptrvalue;
      imp = cjp->imp;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      imp = cbp->imp;
      break;
    case PUB_Patent :
      status = Pub_type_published;
      break;
    default :
      break;
    
  }
  if (imp != NULL)
  {
    if (imp->prepub == 0)
    {
      status = Pub_type_published;
    }
    else if (imp->prepub == 2)
    {
      status = Pub_type_in_press;
    }
    else if (imp->prepub == 1 && the_pub->choice == PUB_Sub)
    {
      status = Pub_type_submitter_block;
    }
    else
    {
      status = Pub_type_unpublished;
    }
    
  }
  return status;
}


static Boolean DoesPubMatchPublicationConstraint (PubdescPtr pdp, PublicationConstraintPtr constraint)
{
  Boolean type_ok = TRUE, rval = FALSE, match_all = TRUE;
  PubPtr pub;
  CharPtr tmp;

  if (pdp == NULL) return FALSE;
  if (IsPublicationConstraintEmpty (constraint)) return TRUE;

  if (constraint->type != Pub_type_any) {
    type_ok = FALSE;
    for (pub = pdp->pub; pub != NULL && !type_ok; pub = pub->next) {
      if (GetPubMLStatus (pub) == constraint->type) {
        type_ok = TRUE;
      }
    }
  }
  if (type_ok) {
    if (constraint->field == NULL) {
      rval = TRUE;
    } else {
      if (constraint->field->constraint->not_present) {
        match_all = TRUE;
        for (pub = pdp->pub; pub != NULL && match_all; pub = pub->next) {
          tmp = GetPubFieldFromPub (pub, constraint->field->field, NULL);
          if (!DoesStringMatchConstraint (tmp, constraint->field->constraint)) {
            match_all = FALSE;
          }
          tmp = MemFree (tmp);
        }
        rval = match_all;
      } else {
        for (pub = pdp->pub; pub != NULL && !rval; pub = pub->next) {
          tmp = GetPubFieldFromPub (pub, constraint->field->field, constraint->field->constraint);
          if (tmp != NULL) {
            rval = TRUE;
          }
          tmp = MemFree (tmp);
        }
      }
    }
  }
  return rval;
}


static Boolean DoesObjectMatchPublicationConstraint (Uint1 choice, Pointer data, PublicationConstraintPtr constraint)
{
  Boolean     rval = TRUE;
  SeqFeatPtr  sfp;
  SeqDescrPtr sdp;

  if (data == NULL) return FALSE;
  if (IsPublicationConstraintEmpty (constraint)) return TRUE;

  switch (choice) {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) data;
      if (sfp->data.choice == SEQFEAT_PUB) {
        rval = DoesPubMatchPublicationConstraint (sfp->data.value.ptrvalue, constraint);
      }
      break;
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) data;
      if (sdp->choice == Seq_descr_pub) {
        rval = DoesPubMatchPublicationConstraint (sdp->data.ptrvalue, constraint);
      }
      break;
  }
  return rval;
}


NLM_EXTERN Boolean IsFieldConstraintEmpty (FieldConstraintPtr constraint)
{
  RnaQualPtr rq;
  FeatureFieldPtr ffp;

  if (constraint == NULL || constraint->field == NULL || IsStringConstraintEmpty (constraint->string_constraint)) {
    return TRUE;
  } else if (constraint->field->choice == FieldType_rna_field
             && ((rq = (RnaQualPtr)constraint->field->data.ptrvalue) == NULL
                 || rq->type == NULL)) {
    return TRUE;
  } else if (constraint->field->choice == FieldType_feature_field
    && (ffp = (FeatureFieldPtr)constraint->field->data.ptrvalue) == NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean DoesObjectMatchFeatureFieldConstraint (Uint1 choice, Pointer data, FeatureFieldPtr ffp, StringConstraintPtr string_constraint)
{
  Boolean           rval = FALSE;
  CharPtr           str;
  BioseqPtr         bsp;
  Int4              subtype;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  Boolean           not_present;
  CGPSetPtr         cgp;
  Uint2             cds_gene_prot_field;

  if (data == NULL) {
    return FALSE;
  }
  if (IsStringConstraintEmpty (string_constraint)) {
    return TRUE;
  }
  
  switch (choice) {
    case OBJ_SEQFEAT:
      str = GetQualFromFeature ((SeqFeatPtr) data, ffp, string_constraint);
      if (str != NULL) {
        rval = TRUE;
        str = MemFree (str);
      }
      break;
    case OBJ_SEQDESC:
      bsp = GetSequenceForObject (choice, data);
      if (bsp != NULL) {
        subtype = GetFeatdefFromFeatureType (ffp->type);
        not_present = string_constraint->not_present;
        string_constraint->not_present = FALSE;
        for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, subtype, &fcontext);
              !rval && sfp != NULL;
              sfp = SeqMgrGetNextFeature (bsp, sfp, 0, subtype, &fcontext)) {
          str = GetQualFromFeature (sfp, ffp, string_constraint);
          if (str != NULL) {
            rval = TRUE;
            str = MemFree (str);
          }
        }
        if (not_present) {
          rval = !rval;
          string_constraint->not_present = TRUE;
        }
      }
      break;
    case 0:
      cgp = (CGPSetPtr) data;
      cds_gene_prot_field = CDSGeneProtFieldFromFeatureField (ffp);
      if (cds_gene_prot_field > 0) {
        not_present = string_constraint->not_present;
        string_constraint->not_present = FALSE;
        str = GetFieldValueFromCGPSet (cgp, cds_gene_prot_field, string_constraint);
        if (str != NULL) {
          rval = TRUE;
          str = MemFree (str);
        }
        if (not_present) {
          rval = !rval;
          string_constraint->not_present = TRUE;
        }
      }
      break;
  }
  return rval;
}


static Boolean DoesObjectMatchFieldConstraint (Uint1 choice, Pointer data, FieldConstraintPtr constraint)
{
  Boolean rval = FALSE;
  BioSourcePtr biop;
  BioseqPtr    bsp;
  CharPtr      str;
  FeatureFieldPtr ffp;

  if (data == NULL) return FALSE;
  if (IsFieldConstraintEmpty (constraint)) {
    return TRUE;
  }

  switch (constraint->field->choice) {
    case FieldType_source_qual:
      biop = GetBioSourceFromObject (choice, data);
      if (biop != NULL) {
        str = GetSourceQualFromBioSource (biop, constraint->field->data.ptrvalue, constraint->string_constraint);
        if (str != NULL) {
          rval = TRUE;
          str = MemFree (str);
        }
      }
      break;
    case FieldType_feature_field:
      ffp = (FeatureFieldPtr) constraint->field->data.ptrvalue;
      rval = DoesObjectMatchFeatureFieldConstraint (choice, data, ffp, constraint->string_constraint);
      break;
    case FieldType_rna_field:
      ffp = FeatureFieldFromRnaQual (constraint->field->data.ptrvalue);
      rval = DoesObjectMatchFeatureFieldConstraint (choice, data, ffp, constraint->string_constraint);
      ffp = FeatureFieldFree (ffp);
      break;
    case FieldType_cds_gene_prot:
      ffp = FeatureFieldFromCDSGeneProtField (constraint->field->data.intvalue);
      rval = DoesObjectMatchFeatureFieldConstraint (choice, data, ffp, constraint->string_constraint);
      ffp = FeatureFieldFree (ffp);
      break;
    case FieldType_molinfo_field:
      bsp = GetSequenceForObject (choice, data);
      if (bsp != NULL) {
        str = GetSequenceQualFromBioseq (bsp, constraint->field->data.ptrvalue);
        if (str == NULL && constraint->string_constraint->not_present) {
          rval = TRUE;
        } else if (str != NULL && DoesStringMatchConstraint (str, constraint->string_constraint)) {
          rval = TRUE;
        }
        str = MemFree (str);
      }
      break;
/* TODO LATER */ 
    case FieldType_pub:
    case FieldType_misc:
      break;
  }
  return rval;    
}


static Boolean DoesObjectMatchConstraint (Uint1 choice, Pointer data, ConstraintChoicePtr constraint)
{
  Boolean rval = TRUE;

  if (data == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  switch (constraint->choice) {
    case ConstraintChoice_string :
      rval = DoesObjectMatchStringConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_location :
      rval = DoesObjectMatchLocationConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_field :
      rval = DoesObjectMatchFieldConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_source :
      rval = DoesBiosourceMatchConstraint (GetBioSourceFromObject (choice, data), constraint->data.ptrvalue);
      break;
    case ConstraintChoice_cdsgeneprot_qual :
      if (choice == 0) {
        rval = DoesCGPSetMatchQualConstraint (data, constraint->data.ptrvalue);
      } else if (choice == OBJ_SEQDESC) {
        rval = DoesSeqDescMatchCGPQualConstraint (data, constraint->data.ptrvalue);
      } else if (choice == OBJ_SEQFEAT) {
        rval = DoesFeatureMatchCGPQualConstraint (data, constraint->data.ptrvalue);
      } else if (choice == OBJ_BIOSEQ) {
        rval = DoesSequenceMatchCGPQualConstraint (data, constraint->data.ptrvalue);
      } else {
        rval = FALSE;
      }
      break;
    case ConstraintChoice_cdsgeneprot_pseudo :
      if (choice == 0) {
        rval = DoesCGPSetMatchPseudoConstraint (data, constraint->data.ptrvalue);
      } else {
        rval = FALSE;
      }
      break;
    case ConstraintChoice_sequence :
      rval = DoesObjectMatchSequenceConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_pub:
      rval = DoesObjectMatchPublicationConstraint (choice, data, constraint->data.ptrvalue);
      break;
  }
  return rval;
}


NLM_EXTERN Boolean DoesObjectMatchConstraintChoiceSet (Uint1 choice, Pointer data, ConstraintChoiceSetPtr csp)
{
  Boolean rval = TRUE;

  if (data == NULL) return FALSE;

  while (csp != NULL && rval) {
    rval = DoesObjectMatchConstraint (choice, data, csp);
    csp = csp->next;
  }
  return rval;
}


NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForField (FieldTypePtr field, ConstraintChoiceSetPtr csp)
{
  StringConstraintPtr scp = NULL;
  SourceConstraintPtr source_constraint;
  CDSGeneProtQualConstraintPtr cgp_constraint;
  PublicationConstraintPtr pub_constraint;
  FieldConstraintPtr       field_constraint;
  FieldType                ft;

  while (csp != NULL) {
    switch (csp->choice) {
      case ConstraintChoice_string :
        scp = csp->data.ptrvalue;
        break;
      case ConstraintChoice_source :
        source_constraint = (SourceConstraintPtr) csp->data.ptrvalue;
        if (source_constraint != NULL && source_constraint->constraint != NULL) {
          if (source_constraint->field1 != NULL) {
            ft.choice = FieldType_source_qual;
            ft.data.ptrvalue = source_constraint->field1;
            ft.next = NULL;
            if (DoFieldTypesMatch (field, &ft)) {
              scp = source_constraint->constraint;
            }
          }
          if (scp == NULL && source_constraint->field2 == NULL) {
            ft.choice = FieldType_source_qual;
            ft.data.ptrvalue = source_constraint->field2;
            ft.next = NULL;
            if (DoFieldTypesMatch (field, &ft)) {
              scp = source_constraint->constraint;
            }
          }
        } 
        break;
      case ConstraintChoice_cdsgeneprot_qual :
        cgp_constraint = (CDSGeneProtQualConstraintPtr) csp->data.ptrvalue;
        if (field->choice == FieldType_cds_gene_prot
            && cgp_constraint != NULL && cgp_constraint->constraint != NULL
            && ((cgp_constraint->field1 != NULL && cgp_constraint->field1->data.intvalue == field->data.intvalue)
                || (cgp_constraint->field2 != NULL && cgp_constraint->field2->data.intvalue == field->data.intvalue))) {
          scp = cgp_constraint->constraint;
        }
        break;
      case ConstraintChoice_pub :
        pub_constraint = csp->data.ptrvalue;
        if (pub_constraint != NULL && pub_constraint->field != NULL) {
          if (field->data.intvalue == pub_constraint->field->field
              && !IsStringConstraintEmpty (pub_constraint->field->constraint)) {
            scp = pub_constraint->field->constraint;
          }
        }
        break;
      case ConstraintChoice_field :
        field_constraint = csp->data.ptrvalue;
        if (field_constraint != NULL
            && field_constraint->field != NULL
            && DoFieldTypesMatch (field, field_constraint->field)) {
          scp = field_constraint->string_constraint;
        }
        break;
    }
    csp = csp->next;
  }
  return scp;
}


NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForFieldPair (FieldPairTypePtr fieldpair, ConstraintChoiceSetPtr csp)
{
  StringConstraintPtr scp;
  FieldTypePtr f;

  f = GetFromFieldFromFieldPair (fieldpair);
  scp = FindStringConstraintInConstraintSetForField (f, csp);
  f = FieldTypeFree (f);
  return scp;
}
 

NLM_EXTERN StringConstraintPtr StringConstraintFromFieldEdit (FieldEditPtr edit)
{
  StringConstraintPtr scp;

  if (edit == NULL || edit->find_txt == NULL) return NULL;
  scp = StringConstraintNew ();
  scp->match_text = StringSave (edit->find_txt);

  switch (edit->location) {
    case Field_edit_location_anywhere :
      scp->match_location = String_location_contains;
      break;
    case Field_edit_location_beginning :
      scp->match_location = String_location_starts;
      break;
    case Field_edit_location_end :
      scp->match_location = String_location_ends;
      break;
  }

  scp->case_sensitive = TRUE;
  scp->whole_word = FALSE;
  scp->not_present = FALSE;

  return scp;
}


static CharPtr ApplyEditToString (CharPtr str, FieldEditPtr edit)
{
  CharPtr cp_found, new_str;
  Int4 found_len, replace_len, new_len;

  if (edit == NULL) return StringSave (str);

  str = StringSave (str);
  cp_found = StringISearch (str, edit->find_txt);

  found_len = StringLen (edit->find_txt);
  replace_len = StringLen (edit->repl_txt);
  if (edit->location == Field_edit_location_beginning
      && cp_found != str) {
    cp_found = NULL;
  } 
  while (cp_found != NULL)
  {
    if (edit->location == Field_edit_location_end
        && cp_found != str + StringLen (str) - found_len) {
      cp_found = StringISearch (cp_found + found_len, edit->find_txt);
    } else {
      new_len = StringLen (str) + 1 - found_len + replace_len;
      new_str = (CharPtr) MemNew (new_len * sizeof (Char));
      if (new_str != NULL)
      {
        if (cp_found != str)
        {
          StringNCpy (new_str, str, cp_found - str);
        }
        StringCat (new_str, edit->repl_txt);
        StringCat (new_str, cp_found + found_len);
        cp_found = new_str + (cp_found - str) + replace_len;
        str = MemFree (str);
        str = new_str;
      }
      cp_found = StringISearch (cp_found, edit->find_txt);
    }
  }
  return str;
}


static void RemoveFieldNameFromString (CharPtr field_name, CharPtr str)
{
  Uint4 field_name_len;
  CharPtr src, dst;

  if (StringHasNoText (field_name) || StringHasNoText (str)) {
    return;
  }
  field_name_len = StringLen (field_name);
    
  if (!StringHasNoText (field_name) && StringNICmp(str, field_name, field_name_len) == 0
        && StringLen (str) > field_name_len 
        && str[field_name_len] == ' ')
  {
    src = str + field_name_len + 1;
    while (*src == ' ')
    {
      src++;
    }
    dst = str;
    while (*src != 0)
    {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }
}


typedef struct objectcollection {
  AECRActionPtr action;
  ValNodePtr object_list;
  BatchExtraPtr batch_extra;
} ObjectCollectionData, PNTR ObjectCollectionPtr;


static void AECRActionObjectCollectionItemCallback (Uint1 objecttype, Pointer objectdata, ObjectCollectionPtr o)
{
  ApplyActionPtr a;
  EditActionPtr e;
  ConvertActionPtr v;
  CopyActionPtr c;
  SwapActionPtr s;
  RemoveActionPtr r;
  AECRParseActionPtr p;
  CharPtr str, portion, field_name;
  StringConstraintPtr scp;
  FieldTypePtr field_from = NULL, field_to = NULL;

  if (objectdata == NULL || o == NULL) return;

  /* check to make sure object is appropriate for field and meets filter */
  switch (o->action->action->choice) {
    case ActionChoice_apply :
      a = (ApplyActionPtr) o->action->action->data.ptrvalue;
      if (a != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, a->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      break;
    case ActionChoice_edit :
      e = (EditActionPtr) o->action->action->data.ptrvalue;
      if (e != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, e->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = StringConstraintFromFieldEdit (e->edit);
        str = GetFieldValueForObjectEx (objecttype, objectdata, e->field, scp, o->batch_extra);
        if (!StringHasNoText (str)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        str = MemFree (str);
      }
      break;
    case ActionChoice_convert :
      v = (ConvertActionPtr) o->action->action->data.ptrvalue;
      if (v != NULL
          && (field_from = GetFromFieldFromFieldPair(v->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = FindStringConstraintInConstraintSetForField (field_from, o->action->constraint);
        str = GetFieldValueForObjectEx (objecttype, objectdata, field_from, scp, o->batch_extra);
        if (v->strip_name) {
          field_to = GetToFieldFromFieldPair (v->fields);
          field_name = SummarizeFieldType (field_to);
          RemoveFieldNameFromString (field_name, str);
          field_name = MemFree (field_name);
          field_to = FieldTypeFree (field_to);
        }
        if (!StringHasNoText (str)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        str = MemFree (str);
      }
      field_from = FieldTypeFree (field_from);
      break;
    case ActionChoice_copy :
      c = (CopyActionPtr) o->action->action->data.ptrvalue;
      if (c != NULL
          && (field_from = GetFromFieldFromFieldPair(c->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(c->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
    case ActionChoice_swap :
      s = (SwapActionPtr) o->action->action->data.ptrvalue;
      if (s != NULL
          && (field_from = GetFromFieldFromFieldPair(s->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(s->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
    case ActionChoice_remove :
      r = (RemoveActionPtr) o->action->action->data.ptrvalue;
      if (r != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, r->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      break;
    case ActionChoice_parse :
      p = (AECRParseActionPtr) o->action->action->data.ptrvalue;
      if (p != NULL
          && (field_from = GetFromFieldFromFieldPair(p->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(p->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = FindStringConstraintInConstraintSetForField (field_from, o->action->constraint);
        str = GetFieldValueForObjectEx (objecttype, objectdata, field_from, scp, o->batch_extra);
        portion = GetTextPortionFromString (str, p->portion);
        if (!StringHasNoText (portion)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        portion = MemFree (portion);
        str = MemFree (str);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
  }

}


static void AECRActionObjectCollectionFeatureCallback (SeqFeatPtr sfp, Pointer data)
{
  ObjectCollectionPtr o;
  if (sfp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_SEQFEAT, sfp, o);

}


static void AECRActionObjectCollectionDescriptorCallback (SeqDescrPtr sdp, Pointer data)
{
  ObjectCollectionPtr o;

  if (sdp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_SEQDESC, sdp, o);
}


static void AECRObjectCollectionBioseqCallback (BioseqPtr bsp, Pointer data)
{
  ObjectCollectionPtr o;

  if (bsp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_BIOSEQ, bsp, o);
}


NLM_EXTERN ValNodePtr GetObjectListForAECRActionEx (SeqEntryPtr sep, AECRActionPtr action, BatchExtraPtr batch_extra)
{
  ObjectCollectionData ocd;

  if (action == NULL) return NULL;

  ocd.action = action;
  ocd.object_list = NULL;
  if (batch_extra == NULL) {
    ocd.batch_extra = BatchExtraNew ();
    InitBatchExtraForAECRAction (ocd.batch_extra, action, sep);
  } else {
    ocd.batch_extra = batch_extra;
  }

  if (FieldTypeFromAECRAction (action) == FieldType_molinfo_field) {
    VisitBioseqsInSep (sep, &ocd, AECRObjectCollectionBioseqCallback);
  } else {
    VisitFeaturesInSep (sep, &ocd, AECRActionObjectCollectionFeatureCallback);
    VisitDescriptorsInSep (sep, &ocd, AECRActionObjectCollectionDescriptorCallback);
  }

  if (batch_extra != ocd.batch_extra) {
    ocd.batch_extra = BatchExtraFree (ocd.batch_extra);
  }
  return ocd.object_list;
}


NLM_EXTERN ValNodePtr GetObjectListForAECRAction (SeqEntryPtr sep, AECRActionPtr action)
{
  return GetObjectListForAECRActionEx (sep, action, NULL);
}


NLM_EXTERN ValNodePtr FreeObjectList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    if (vnp->choice == 0) {
      vnp->data.ptrvalue = CGPSetFree (vnp->data.ptrvalue);
    }
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


typedef struct buildcgpset
{
  ValNodePtr cds_list;
  ValNodePtr mrna_list;
  ValNodePtr gene_list;
} BuildCGPSetData, PNTR BuildCGPSetPtr;

static void BuildCGPSetCallback (SeqFeatPtr sfp, Pointer userdata)
{
  BuildCGPSetPtr b;

  if (sfp == NULL || sfp->idx.deleteme || userdata == NULL) return;
  b = (BuildCGPSetPtr) userdata;
  if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    ValNodeAddPointer (&(b->cds_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_GENE)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->idx.subtype == FEATDEF_mRNA)
  {
    ValNodeAddPointer (&(b->mrna_list), OBJ_SEQFEAT, sfp);
  }
  else if (SeqMgrGetGeneXref (sfp) != NULL)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
}


static CGPSetPtr BuildCGPSetFromCodingRegion (SeqFeatPtr cds, BoolPtr indexing_needed)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        gene = NULL, mrna, prot;
  BioseqPtr         protbsp;
  CGPSetPtr         cdsp;
  ProtRefPtr        prp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return NULL;

  cdsp = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  ValNodeAddPointer (&(cdsp->cds_list), 0, cds);

  gene = GetGeneForFeature (cds);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  mrna = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
  if (mrna != NULL)
  {
    ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);
    /* mark mrna, so that we'll know it's already in a set */
    mrna->idx.deleteme = TRUE;
  }

  if (cds->product != NULL)
  {
    protbsp = BioseqFindFromSeqLoc (cds->product);
    if (protbsp != NULL)
    {
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      /* if there is no full-length protein feature, make one */
      if (prot == NULL)
      {
        prp = ProtRefNew ();
        prot = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
        if (prot != NULL)
        {
          prot->data.value.ptrvalue = prp;
          if (indexing_needed != NULL)
          {
            *indexing_needed = TRUE;
          }
        }
      }
      if (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
      }
      
      /* also add in mat_peptides from protein feature */
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      while (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
        prot = SeqMgrGetNextFeature (protbsp, prot, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      }
    }
  }  
  return cdsp;
}


static CGPSetPtr BuildCGPSetFrommRNA (SeqFeatPtr mrna)
{
  SeqFeatPtr        gene;
  CGPSetPtr          cdsp;

  if (mrna == NULL || mrna->idx.deleteme || mrna->idx.subtype != FEATDEF_mRNA) return NULL;

  cdsp = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);

  gene = GetGeneForFeature (mrna);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  return cdsp;
}


static CGPSetPtr BuildCGPSetFromGene (SeqFeatPtr gene)
{
  CGPSetPtr          cdsp;

  if (gene == NULL || gene->idx.deleteme || gene->idx.subtype != FEATDEF_GENE) {
    return NULL;
  }

  cdsp = CGPSetNew ();
  ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
  return cdsp;
}


static void UnmarkFeatureList (ValNodePtr list)
{
  SeqFeatPtr sfp;

  while (list != NULL)
  {
    sfp = list->data.ptrvalue;
    if (sfp != NULL)
    {
      sfp->idx.deleteme = FALSE;
    }
    list = list->next;
  }
}


static void 
AdjustCGPObjectListForMatPeptides 
(ValNodePtr PNTR cgp_list,
 FieldTypePtr field1, 
 FieldTypePtr field2,
 ConstraintChoiceSetPtr constraints)
{
  ConstraintChoiceSetPtr mat_peptide_constraints = NULL;
  Boolean found_matpeptide_constraint = FALSE;
  ValNodePtr vnp, vnp_prev, vnp_next;
  ValNodePtr m_vnp, m_vnp_prev, m_vnp_next, mat_peptide_list;
  CGPSetPtr  cdsp;
  SeqFeatPtr sfp;

  if (cgp_list == NULL 
      || *cgp_list == NULL
      || constraints == NULL
      || (field1 == NULL && field2 == NULL) /* no fields specified */
      || (!IsFieldTypeMatPeptideRelated (field1) && !IsFieldTypeMatPeptideRelated(field2))) {
    return;
  }

  /* get list of constraints that apply to mat-peptide features */
  while (constraints != NULL) {
    if (IsConstraintChoiceMatPeptideRelated (constraints)) {        
      ValNodeLink (&mat_peptide_constraints, AsnIoMemCopy (constraints, (AsnReadFunc) ConstraintChoiceAsnRead, (AsnWriteFunc) ConstraintChoiceAsnWrite));
    }
    constraints = constraints->next;
  }
  if (mat_peptide_constraints == NULL) {
    return;
  }

  /* if both fields are mat-peptide related, or one is mat-peptide related and the other is NULL,
   * convert sets to lists of mat-peptide features 
   * otherwise just remove mat-peptide features from the prot-list that do not match the constraints.
   */
  if ((field1 != NULL && !IsFieldTypeMatPeptideRelated (field1))
      || (field2 != NULL && !IsFieldTypeMatPeptideRelated (field2))) {
    for (vnp = *cgp_list; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 0) {
        cdsp = (CGPSetPtr) vnp->data.ptrvalue;
        m_vnp_prev = NULL;
        for (m_vnp = cdsp->prot_list; m_vnp != NULL; m_vnp = m_vnp_next) {
          m_vnp_next = m_vnp->next;
          sfp = m_vnp->data.ptrvalue;
          if (sfp == NULL
              || (sfp->idx.subtype == FEATDEF_mat_peptide_aa 
                  && !DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, mat_peptide_constraints))) {
            if (m_vnp_prev == NULL) {
              cdsp->prot_list = m_vnp->next;
            } else {
              m_vnp_prev->next = m_vnp->next;
            }
            m_vnp->next = NULL;
            m_vnp = ValNodeFree (m_vnp);
          } else {
            m_vnp_prev = m_vnp;
          }
        }
      }
    }
  } else {
    vnp_prev = NULL;
    for (vnp = *cgp_list; vnp != NULL; vnp = vnp_next) {
      vnp_next = vnp->next;
      if (vnp->choice == 0) {
        mat_peptide_list = NULL;
        cdsp = (CGPSetPtr) vnp->data.ptrvalue;
        for (m_vnp = cdsp->prot_list; m_vnp != NULL; m_vnp = m_vnp->next) {
          sfp = m_vnp->data.ptrvalue;
          if (sfp->idx.subtype == FEATDEF_mat_peptide_aa
              && DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, mat_peptide_constraints)) {
            ValNodeAddPointer (&mat_peptide_list, OBJ_SEQFEAT, sfp);
          }
        }
        if (mat_peptide_list == NULL) {
          if (vnp_prev == NULL) {
            *cgp_list = vnp->next;
          } else {
            vnp_prev->next = vnp->next;
          }
          vnp->next = NULL;
          vnp = FreeObjectList (vnp);
        } else {
          m_vnp = mat_peptide_list;
          while (m_vnp->next != NULL) {
            m_vnp = m_vnp->next;
          }
          if (vnp_prev == NULL) {
            *cgp_list = mat_peptide_list;
          } else {
            vnp_prev->next = mat_peptide_list;
          }
          m_vnp->next = vnp_next;
          vnp_prev = m_vnp;
          vnp->next = NULL;
          vnp = FreeObjectList (vnp);
        }
      } else {
        vnp_prev = vnp;
      }
    }
  }
  mat_peptide_constraints = ConstraintChoiceSetFree (mat_peptide_constraints);
}


static ValNodePtr BuildCGPSetList (Uint2 entityID, AECRActionPtr act)
{
  SeqEntryPtr    sep;
  BuildCGPSetData b;
  CGPSetPtr       cdsp;
  ValNodePtr     vnp, vnp_next, vnp_prev;
  ValNodePtr     cdset_list = NULL;
  SeqFeatPtr     cds, gene, mrna;
  Boolean        need_indexing = FALSE;
  ApplyActionPtr      a;
  EditActionPtr       e;
  ConvertActionPtr    c;
  CopyActionPtr       cp;
  SwapActionPtr       s;
  AECRParseActionPtr  pa;
  RemoveActionPtr     r;
  FieldTypePtr        field_from, field_to;
  
  sep = GetTopSeqEntryForEntityID (entityID);

  b.cds_list = NULL;
  b.gene_list = NULL;
  b.mrna_list = NULL;
  
  VisitFeaturesInSep (sep, &b, BuildCGPSetCallback);

  /* build cdsets that have coding regions */
  for (vnp = b.cds_list; vnp != NULL; vnp = vnp->next)
  {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    if (cds == NULL) continue;
    cdsp = BuildCGPSetFromCodingRegion (cds, &need_indexing);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }
  if (need_indexing)
  {
    /* indexing because we have created full-length protein features */
    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* build cdsets for mrna features that don't have coding regions */
  for (vnp = b.mrna_list; vnp != NULL; vnp = vnp->next)
  {
    mrna = (SeqFeatPtr) vnp->data.ptrvalue;
    if (mrna == NULL || mrna->idx.deleteme) continue;
    cdsp = BuildCGPSetFrommRNA (mrna);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }

  /* build cdsets for lonely genes / features with gene xrefs that are not coding regions or mrnas */
  for (vnp = b.gene_list; vnp != NULL; vnp = vnp->next)
  {
    gene = (SeqFeatPtr) vnp->data.ptrvalue;
    if (gene == NULL || gene->idx.deleteme) continue;
    cdsp = BuildCGPSetFromGene (gene);
    if (cdsp != NULL) {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }

  /* now unmark features */
  UnmarkFeatureList (b.cds_list);
  UnmarkFeatureList (b.mrna_list);
  UnmarkFeatureList (b.gene_list);

  b.cds_list = ValNodeFree (b.cds_list);
  b.mrna_list = ValNodeFree (b.mrna_list);
  b.gene_list = ValNodeFree (b.gene_list);

  /* now remove sets that don't match our choice constraint */
  if (act != NULL && act->constraint != NULL) {
    vnp_prev = NULL;
    for (vnp = cdset_list; vnp != NULL; vnp = vnp_next)
    {
      vnp_next = vnp->next;
      if (!DoesObjectMatchConstraintChoiceSet (0, vnp->data.ptrvalue, act->constraint))
      {
        if (vnp_prev == NULL)
        {
          cdset_list = vnp->next;
        }
        else
        {
          vnp_prev->next = vnp->next;
        }
        vnp->next = NULL;
        FreeCGPSetList (vnp);     
      }
      else
      {
        vnp_prev = vnp;
      }
    }
  }

  /* adjust if action fields are mat-peptide specific */
  if (act != NULL && act->action != NULL && act->action->data.ptrvalue != NULL) {
    switch (act->action->choice) {
      case ActionChoice_apply:
        a = (ApplyActionPtr) act->action->data.ptrvalue;
        AdjustCGPObjectListForMatPeptides (&cdset_list, a->field, NULL, act->constraint);
        break;
      case ActionChoice_edit:
        e = (EditActionPtr) act->action->data.ptrvalue;
        AdjustCGPObjectListForMatPeptides (&cdset_list, e->field, NULL, act->constraint);
        break;
       case ActionChoice_convert:
        c = (ConvertActionPtr) act->action->data.ptrvalue;
        field_from = GetFromFieldFromFieldPair (c->fields);
        field_to = GetToFieldFromFieldPair (c->fields);
        AdjustCGPObjectListForMatPeptides (&cdset_list, field_from, field_to, act->constraint);
        field_from = FieldTypeFree (field_from);
        field_to = FieldTypeFree (field_to);
        break;
       case ActionChoice_copy:
        cp = (CopyActionPtr) act->action->data.ptrvalue;
        field_from = GetFromFieldFromFieldPair (cp->fields);
        field_to = GetToFieldFromFieldPair (cp->fields);
        AdjustCGPObjectListForMatPeptides (&cdset_list, field_from, field_to, act->constraint);
        field_from = FieldTypeFree (field_from);
        field_to = FieldTypeFree (field_to);
        break;
       case ActionChoice_swap:
        s = (SwapActionPtr) act->action->data.ptrvalue;
        field_from = GetFromFieldFromFieldPair (s->fields);
        field_to = GetToFieldFromFieldPair (s->fields);
        AdjustCGPObjectListForMatPeptides (&cdset_list, field_from, field_to, act->constraint);
        field_from = FieldTypeFree (field_from);
        field_to = FieldTypeFree (field_to);
        break;
      case ActionChoice_remove:
        r = (RemoveActionPtr) act->action->data.ptrvalue;
        AdjustCGPObjectListForMatPeptides (&cdset_list, r->field, NULL, act->constraint);
        break;
      case ActionChoice_parse:
        pa = (AECRParseActionPtr) act->action->data.ptrvalue;
        field_from = GetFromFieldFromFieldPair (pa->fields);
        field_to = GetToFieldFromFieldPair (pa->fields);
        AdjustCGPObjectListForMatPeptides (&cdset_list, field_from, field_to, act->constraint);
        field_from = FieldTypeFree (field_from);
        field_to = FieldTypeFree (field_to);
        break;
    }
  }  
  return cdset_list;
}


static void AlsoChangeMrnaForObject (Uint1 choice, Pointer data)
{
  CharPtr str;
  SeqFeatPtr sfp, mrna;
  SeqMgrFeatContext context;
  FeatureField f;

  if (choice == 0) {
    str = GetFieldValueFromCGPSet (data, CDSGeneProt_field_prot_name, NULL);
    SetFieldValueInCGPSet (data, CDSGeneProt_field_mrna_product, NULL, str, ExistingTextOption_replace_old);
    str = MemFree (str);
  } else if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      mrna = SeqMgrGetOverlappingmRNA (sfp->location, &context);
      if (mrna != NULL) {
        f.type = Feature_type_cds;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = Feat_qual_legal_product;  
        str = GetQualFromFeature (sfp, &f, NULL);
        f.type = Feature_type_mRNA;
        SetQualOnFeature (mrna, &f, NULL, str, ExistingTextOption_replace_old);
        str = MemFree (str);
        f.field = ValNodeFree (f.field);
      }
    }
  }
}


NLM_EXTERN Int4 DoApplyActionToObjectListEx (ApplyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;

  if (action == NULL || object_list == NULL) return 0;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, action->field, scp, action->value, action->existing_text, batch_extra)) {
      if (also_change_mrna) {
        AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
      }
      num_succeed ++;
    } else {
      num_fail++;
    }
  }

  return num_succeed;
}


NLM_EXTERN Int4 DoApplyActionToObjectList (ApplyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  return DoApplyActionToObjectListEx (action, object_list, also_change_mrna, scp, NULL);
}


NLM_EXTERN Int4 DoEditActionToObjectListEx (EditActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  StringConstraintPtr scp;
  CharPtr    str, new_str;

  if (action == NULL || object_list == NULL) return 0;

  scp = StringConstraintFromFieldEdit (action->edit);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, action->field, scp, batch_extra);
    new_str = ApplyEditToString (str, action->edit);
    if (StringCmp (str, new_str) != 0
        && SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, action->field, scp, new_str, ExistingTextOption_replace_old, batch_extra)) {
      if (also_change_mrna) {
        AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
      }
      num_succeed ++;
    } else {
      num_fail++;
    }
    new_str = MemFree (new_str);
    str = MemFree (str);
  }
  return num_succeed;
}


NLM_EXTERN Int4 DoEditActionToObjectList (EditActionPtr action, ValNodePtr object_list, Boolean also_change_mrna)
{
  return DoEditActionToObjectListEx (action, object_list, also_change_mrna, NULL);
}


NLM_EXTERN Int4 DoConvertActionToObjectListEx (ConvertActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str, from_val, field_name = NULL;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL || action->fields == NULL) return 0;

  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  if (action->strip_name) {
    field_name = SummarizeFieldType (field_to);
  }

  if (action->fields->choice == FieldPairType_molinfo_field) {
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
      str = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, NULL, batch_extra);
      from_val = GetSequenceQualValName (field_from->data.ptrvalue);
      if (StringCmp (str, from_val) == 0
          && SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, ExistingTextOption_replace_old, batch_extra)) {
        num_succeed ++;
      }
      str = MemFree (str);
    }
  } else {
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
      str = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, batch_extra);
      if (action->strip_name) {
        RemoveFieldNameFromString (field_name, str);
      }
      if (SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, action->existing_text, batch_extra)
          && (action->keep_original || RemoveFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp))) {
        if (also_change_mrna) {
          AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
        }
        num_succeed ++;
      } else {
        num_fail++;
      }
      str = MemFree (str);
    }
  }

  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  field_name = MemFree (field_name);

  return num_succeed;
}


NLM_EXTERN Int4 DoConvertActionToObjectList (ConvertActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  return DoConvertActionToObjectListEx (action, object_list, also_change_mrna, scp, NULL);
}


NLM_EXTERN Int4 DoCopyActionToObjectListEx (CopyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, batch_extra);
    if (SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, action->existing_text, batch_extra)) {
      if (also_change_mrna) {
        AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
      }
      num_succeed ++;
    } else {
      num_fail++;
    }
    str = MemFree (str);
  }

  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


NLM_EXTERN Int4 DoCopyActionToObjectList (CopyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  return DoCopyActionToObjectListEx (action, object_list, also_change_mrna, scp, NULL);
}

NLM_EXTERN Int4 DoSwapActionToObjectListEx (SwapActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str1, str2;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, batch_extra);
    str2 = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, batch_extra);
    if (SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str1, ExistingTextOption_replace_old, batch_extra)
        && SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, str2, ExistingTextOption_replace_old, batch_extra)) {
      if (also_change_mrna) {
        AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
      }
      num_succeed ++;
    } else {
      num_fail++;
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


NLM_EXTERN Int4 DoSwapActionToObjectList (SwapActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  return DoSwapActionToObjectListEx (action, object_list, also_change_mrna, scp, NULL);
}


NLM_EXTERN Int4 DoRemoveActionToObjectList (RemoveActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;

  if (action == NULL || object_list == NULL) return 0;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (RemoveFieldValueForObject (vnp->choice, vnp->data.ptrvalue, action->field, scp)) {
      if (also_change_mrna) {
        AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
      }
      num_succeed ++;
    } else {
      num_fail++;
    }
  }
  return num_succeed;
}


NLM_EXTERN Int4 DoParseActionToObjectListEx (AECRParseActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra)
{
  ValNodePtr vnp;
  CharPtr    str1, str2, cp;
  Int4       len, num_succeed = 0, diff;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, batch_extra);
    str2 = GetTextPortionFromString (str1, action->portion);    
    if (str2 != NULL) {
      if (action->remove_from_parsed) {
        cp = FindTextPortionLocationInString (str1, action->portion);
        if (cp != NULL) {
          len = StringLen (str2);
          if (action->remove_left &&action->portion != NULL && action->portion->left_text != NULL && !action->portion->include_left) {
            diff = StringLen (action->portion->left_text);
            cp -= diff;
            len += diff;
          }
          if (action->remove_right &&action->portion != NULL && action->portion->right_text != NULL && !action->portion->include_right) {
            diff = StringLen (action->portion->right_text);
            len += diff;
          }
          StringCpy (cp, cp + len);
          SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_from, scp, str1, ExistingTextOption_replace_old, batch_extra);
        }
      }
      if (SetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str2, action->existing_text, batch_extra)) {
        if (also_change_mrna) {
          AlsoChangeMrnaForObject (vnp->choice, vnp->data.ptrvalue);
        }
        num_succeed++;
      }
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


NLM_EXTERN Int4 DoParseActionToObjectList (AECRParseActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp)
{
  return DoParseActionToObjectListEx (action, object_list, also_change_mrna, scp, NULL);
}


static Int4 ApplyAECRActionToSeqEntry (AECRActionPtr act, SeqEntryPtr sep)
{
  StringConstraintPtr scp;
  ApplyActionPtr      a;
  ConvertActionPtr    c;
  RemoveActionPtr     r;
  ValNodePtr          object_list = NULL;
  Uint1               field_type;
  Uint2               entityID;
  Int4                num_succeed = 0;
  FieldTypePtr        field_from;
  BatchExtraPtr       batch_extra;

  if (act == NULL || act->action == NULL) return 0;

  batch_extra = BatchExtraNew ();
  InitBatchExtraForAECRAction (batch_extra, act, sep);

  field_type = FieldTypeFromAECRAction (act);
  if (field_type == FieldType_cds_gene_prot) {
    entityID = ObjMgrGetEntityIDForChoice(sep);
    object_list = BuildCGPSetList (entityID, act);
    
  } else {
    object_list = GetObjectListForAECRActionEx (sep, act, batch_extra);
  }


  switch (act->action->choice) {
    case ActionChoice_apply:
      a = (ApplyActionPtr) act->action->data.ptrvalue;
      scp = FindStringConstraintInConstraintSetForField (a->field, act->constraint);
      num_succeed = DoApplyActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, scp, batch_extra);
      break;
    case ActionChoice_edit:
      num_succeed = DoEditActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, batch_extra);
      break;
    case ActionChoice_convert:
      scp = NULL;
      if (act->constraint != NULL) {
        c = (ConvertActionPtr) act->action->data.ptrvalue;
        field_from = GetFromFieldFromFieldPair (c->fields);
        scp = FindStringConstraintInConstraintSetForField (field_from, act->constraint);
        field_from = FieldTypeFree (field_from);
      }
      num_succeed = DoConvertActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, scp, batch_extra);
      break;
    case ActionChoice_swap:
      num_succeed = DoSwapActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, NULL, batch_extra);
      break;
    case ActionChoice_copy:
      num_succeed = DoCopyActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, NULL, batch_extra);
      break;
    case ActionChoice_remove:
      r = (RemoveActionPtr) act->action->data.ptrvalue;
      scp = FindStringConstraintInConstraintSetForField (r->field, act->constraint);
      num_succeed = DoRemoveActionToObjectList (act->action->data.ptrvalue, object_list, act->also_change_mrna, scp);
      break;
    case ActionChoice_parse:
      num_succeed = DoParseActionToObjectListEx (act->action->data.ptrvalue, object_list, act->also_change_mrna, NULL, batch_extra);
      break;
  }
  object_list = FreeObjectList (object_list);
  batch_extra = BatchExtraFree (batch_extra);
  return num_succeed;
}


static AECRSamplePtr AECRSampleNew (void)
{
  AECRSamplePtr sample;

  sample = (AECRSamplePtr) MemNew (sizeof (AECRSampleData));
  MemSet (sample, 0, sizeof (AECRSampleData));
  sample->all_same = TRUE;
  return sample;
}


NLM_EXTERN AECRSamplePtr AECRSampleFree (AECRSamplePtr sample)
{
  if (sample != NULL) {
    sample->field = FieldTypeFree (sample->field);
    sample->first_value = MemFree (sample->first_value);
    sample = MemFree (sample);
  }
  return sample;
}


NLM_EXTERN ValNodePtr AECRSampleListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = AECRSampleFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static void AddTextToAECRSample (AECRSamplePtr sample, CharPtr txt)
{
  if (StringHasNoText (txt)) {
    txt = MemFree (txt);
  } else if (sample != NULL) {
    sample->num_found ++;
    if (sample->first_value == NULL) {
      sample->first_value = txt;
    } else {
      if (sample->all_same && StringCmp (sample->first_value, txt) != 0) {
        sample->all_same = FALSE;
      }
      txt = MemFree (txt);
    }
  }
}


NLM_EXTERN AECRSamplePtr GetAECRSampleFromObjectListEx (ValNodePtr object_list, FieldTypePtr field, BatchExtraPtr batch_extra)
{
  AECRSamplePtr sample;
  ValNodePtr    vnp, prot_vnp, bsp_list;
  CharPtr       txt;
  CGPSetPtr     cgp;
  SeqFeatPtr    sfp;
  BatchExtraPtr b = NULL;
  SeqEntryPtr   sep;

  if (object_list == NULL || field == NULL) return NULL;

  if (batch_extra == NULL) {
    b = BatchExtraNew ();
    batch_extra = b;
    bsp_list = BioseqListForObjectList (object_list);
    for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
      sep = SeqMgrGetSeqEntryForData (vnp->data.ptrvalue);
      InitBatchExtraForField (batch_extra, field, sep);
    }
    bsp_list = ValNodeFree (bsp_list);
  }

  sample = AECRSampleNew ();
  sample->field = FieldTypeCopy (field);
  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0 && IsFieldTypeMatPeptideRelated (field)) {
      cgp = (CGPSetPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        for (prot_vnp = cgp->prot_list; prot_vnp != NULL; prot_vnp = prot_vnp->next) {
          sfp = (SeqFeatPtr) prot_vnp->data.ptrvalue;
          if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa) {
            txt = GetFieldValueForObjectEx (OBJ_SEQFEAT, sfp, field, NULL, batch_extra);
            AddTextToAECRSample (sample, txt);
          }
        }
      }
    } else {
      txt = GetFieldValueForObjectEx (vnp->choice, vnp->data.ptrvalue, field, NULL, batch_extra);
      AddTextToAECRSample (sample, txt);
    }
  }

  b = BatchExtraFree (b);
  return sample;
}


NLM_EXTERN AECRSamplePtr GetAECRSampleFromObjectList (ValNodePtr object_list, FieldTypePtr field)
{
  return GetAECRSampleFromObjectListEx (object_list, field, NULL);
}


static void GetFieldsFromAECR (AECRActionPtr act, FieldTypePtr PNTR pField, ValNodePtr PNTR pFieldPair)
{
  ApplyActionPtr     a;
  EditActionPtr      e;
  ConvertActionPtr   c;
  SwapActionPtr      s;
  CopyActionPtr      cp;
  RemoveActionPtr    r;
  AECRParseActionPtr p;
  
  if (pField != NULL) {
    *pField = NULL;
  }
  if (pFieldPair != NULL) {
    *pFieldPair = NULL;
  }
  if (act == NULL || act->action == NULL || act->action->data.ptrvalue == NULL) {
    return;
  }

  switch (act->action->choice) {
    case ActionChoice_apply:
      if (pField != NULL) {
        a = (ApplyActionPtr) act->action->data.ptrvalue;
        *pField = a->field;
      }
      break;
    case ActionChoice_edit:
      if (pField != NULL) {
        e = (EditActionPtr) act->action->data.ptrvalue;
        *pField = e->field;
      }
      break;
    case ActionChoice_convert:
      if (pFieldPair != NULL) {
        c = (ConvertActionPtr) act->action->data.ptrvalue;
        *pFieldPair = c->fields;
      }
      break;
    case ActionChoice_swap:
      if (pFieldPair != NULL) {
        s = (SwapActionPtr) act->action->data.ptrvalue;
        *pFieldPair = s->fields;
      }
      break;
    case ActionChoice_copy:
      if (pFieldPair != NULL) {
        cp = (CopyActionPtr) act->action->data.ptrvalue;
        *pFieldPair = cp->fields;
      }
      break;
    case ActionChoice_remove:
      if (pField != NULL) {
        r = (RemoveActionPtr) act->action->data.ptrvalue;
        *pField = r->field;
      }
      break;
    case ActionChoice_parse:
      if (pFieldPair != NULL) {
        p = (AECRParseActionPtr) act->action->data.ptrvalue;
        *pFieldPair = p->fields;
      }
      break;
  }
}


NLM_EXTERN ValNodePtr LIBCALLBACK FieldTypeListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list = FieldTypeFree (list);
    list = list_next;
  }
  return list;
}


NLM_EXTERN ValNodePtr LIBCALLBACK FieldTypeListCopy (ValNodePtr orig)
{
  ValNodePtr prev = NULL, new_list = NULL, vnp;

  while (orig != NULL) {
    vnp = FieldTypeCopy (orig);
    if (prev == NULL) {
      new_list = vnp;
    } else {
      prev->next = vnp;
    }
    prev = vnp;
    orig = orig->next;
  }
  return new_list;
}


static int LIBCALLBACK SortVnpByChoiceAndIntvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 == NULL && vnp2 == NULL) {
      rval = 0;
    } else if (vnp1 == NULL) {
      rval = -1;
    } else if (vnp2 == NULL) {
      rval = 1;
    } else if (vnp1->choice > vnp2->choice) {
      rval = 1;
    } else if (vnp1->choice < vnp2->choice) {
      rval = -1;
    } else if (vnp1->data.intvalue > vnp2->data.intvalue) {
      rval = 1;
    } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
      rval = -1;
    } else {
      rval = 0;
    }
  }
  return rval;
}


/* Callback function used for sorting and uniqueing */

static int LIBCALLBACK SortVnpByFieldType (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = CompareFieldTypes (vnp1, vnp2);
  }

  return rval;
}


static void GetBioSourceFields (BioSourcePtr biop, Pointer userdata)
{

  if (biop == NULL || userdata == NULL) {
    return;
  }

  ValNodeLink ((ValNodePtr PNTR) userdata, GetSourceQualFieldListFromBioSource (biop));
}


NLM_EXTERN void SortUniqueFieldTypeList (ValNodePtr PNTR field_list)
{
  if (field_list == NULL) return;
  *field_list = ValNodeSort (*field_list, SortVnpByFieldType);
  ValNodeUnique (field_list, SortVnpByFieldType, FieldTypeListFree);
}


NLM_EXTERN ValNodePtr GetSourceQualSampleFieldList (SeqEntryPtr sep)
{
  ValNodePtr field_list = NULL;
  ValNodePtr vnp_prev = NULL, vnp, sq;
  Boolean    done = FALSE;

  VisitBioSourcesInSep (sep, &field_list, GetBioSourceFields);
  field_list = ValNodeSort (field_list, SortVnpByFieldType);
  ValNodeUnique (&field_list, SortVnpByFieldType, FieldTypeListFree);

  /* rearrange so that taxname is always first */
  for (vnp = field_list; vnp != NULL && !done; vnp = vnp->next) {
    if (vnp->choice == FieldType_source_qual
        && (sq = vnp->data.ptrvalue) != NULL
        && sq->choice == SourceQualChoice_textqual
        && sq->data.intvalue == Source_qual_taxname) {
      if (vnp_prev != NULL) {
        vnp_prev->next = vnp->next;
        vnp->next = field_list;
        field_list = vnp;
      }
      done = TRUE;
    } else {
      vnp_prev = vnp;
    }
  }

  return field_list;
}


NLM_EXTERN ValNodePtr GetSourceQualSampleFieldListForSeqEntryList (ValNodePtr list)
{
  ValNodePtr field_list = NULL;
  ValNodePtr vnp_prev = NULL, vnp, sq;
  Boolean    done = FALSE;

  if (list == NULL) {
    return NULL;
  }

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    VisitBioSourcesInSep (vnp->data.ptrvalue, &field_list, GetBioSourceFields);
  }
  field_list = ValNodeSort (field_list, SortVnpByFieldType);
  ValNodeUnique (&field_list, SortVnpByFieldType, FieldTypeListFree);

  /* rearrange so that taxname is always first */
  for (vnp = field_list; vnp != NULL && !done; vnp = vnp->next) {
    if (vnp->choice == FieldType_source_qual
        && (sq = vnp->data.ptrvalue) != NULL
        && sq->choice == SourceQualChoice_textqual
        && sq->data.intvalue == Source_qual_taxname) {
      if (vnp_prev != NULL) {
        vnp_prev->next = vnp->next;
        vnp->next = field_list;
        field_list = vnp;
      }
      done = TRUE;
    } else {
      vnp_prev = vnp;
    }
  }

  return field_list;
}


static void GetFeatureQualFieldListForAECRSampleCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr PNTR list;

  list = (ValNodePtr PNTR) data;
  if (list == NULL || sfp == NULL) return;

  ValNodeLink (list, GetFieldListFromFeature (sfp));
}

static ValNodePtr GetFeatureQualFieldList (SeqEntryPtr sep)
{
  ValNodePtr field_list = NULL;

  VisitFeaturesInSep (sep, &field_list, GetFeatureQualFieldListForAECRSampleCallback);
  field_list = ValNodeSort (field_list, SortVnpByFieldType);
  ValNodeUnique (&field_list, SortVnpByFieldType, FieldTypeListFree);
  return field_list;
}


static void GetRnaQualFieldListForAECRSampleCallback (SeqFeatPtr sfp, Pointer userdata)
{
  RnaFeatTypePtr type;
  RnaRefPtr      rrp;
  RnaQualPtr     rq;
  GBQualPtr      gbqual;
  GeneRefPtr     grp = NULL;
  SeqFeatPtr     gene = NULL;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL || userdata == NULL) {
    return;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  type = RnaFeatTypeFromSeqFeat (sfp);

  if (type == NULL) return;

  /* add product if appropriate */
  if ((type->choice == RnaFeatType_preRNA || type->choice == RnaFeatType_mRNA
        || type->choice == RnaFeatType_rRNA || type->choice == RnaFeatType_miscRNA)
       && rrp->ext.choice == 1
       && !StringHasNoText (rrp->ext.value.ptrvalue)) {
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_product;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
  } else if (type->choice == RnaFeatType_ncRNA || type->choice == RnaFeatType_tmRNA) {
    gbqual = sfp->qual;
    while (gbqual != NULL && StringCmp (gbqual->qual, "product") != 0) {
      gbqual = gbqual->next;
    }
    if (gbqual != NULL) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_product;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
  }

  /* add comment if present */
  if (!StringHasNoText (sfp->comment)) {
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_comment;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
  }

  /* add tRNA specific if appropriate */
  if (type->choice == RnaFeatType_tRNA) {
    /* codons recognized */
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_codons_recognized;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);

    /* anticodon */
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_anticodon;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
  }

  /* add ncRNA class if appropriate and present */
  if (type->choice == RnaFeatType_ncRNA) {
    gbqual = sfp->qual;
    while (gbqual != NULL && StringCmp (gbqual->qual, "ncRNA_class") != 0) {
      gbqual = gbqual->next;
    }
    if (gbqual != NULL) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_ncrna_class;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
  }

  /* add transcript ID if present */
  if (sfp->product != NULL) {
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_transcript_id;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
  }

  /* add gene fields */
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL) {
    gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
    if (gene != NULL) {
      grp = gene->data.value.ptrvalue;
    }
  }
  if (grp != NULL && !SeqMgrGeneIsSuppressed (grp)) {
    /* gene locus */
    if (!StringHasNoText (grp->locus)) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_gene_locus;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
    /* gene description */
    if (!StringHasNoText (grp->desc)) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_gene_locus;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
    /* maploc */
    if (!StringHasNoText (grp->maploc)) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_gene_maploc;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
    /* locus tag */
    if (!StringHasNoText (grp->locus_tag)) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_gene_locus_tag;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
    /* synonym */
    if (grp->syn != NULL) {
      rq = RnaQualNew ();
      rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
      rq->field = Rna_field_gene_synonym;
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
    }
  }

  /* gene comment */
  if (gene != NULL && !StringHasNoText (gene->comment)) {
    rq = RnaQualNew ();
    rq->type = AsnIoMemCopy (type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
    rq->field = Rna_field_gene_comment;
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, FieldType_rna_field, rq);
  }
}


static ValNodePtr GetRnaQualFieldList (SeqEntryPtr sep)
{
  ValNodePtr field_list = NULL;

  VisitFeaturesInSep (sep, &field_list, GetRnaQualFieldListForAECRSampleCallback);
  field_list = ValNodeSort (field_list, SortVnpByFieldType);
  ValNodeUnique (&field_list, SortVnpByFieldType, FieldTypeListFree);
  return field_list;
}


static void GetStructuredCommentFieldsCallback (SeqDescrPtr sdp, Pointer data)
{
  UserObjectPtr uop;
  UserFieldPtr  ufp;
  ValNodePtr    vnp;

  if (sdp != NULL && data != NULL && sdp->choice == Seq_descr_user
      && (uop = sdp->data.ptrvalue) != NULL
      && IsUserObjectStructuredComment (uop)) {

    ufp = uop->data;
    while (ufp != NULL) {
      if (ufp->label != NULL && ufp->label->str != NULL
          && StringCmp (ufp->label->str, "StructuredCommentPrefix") != 0
          && StringCmp (ufp->label->str, "StructuredCommentSuffix") != 0) {
        vnp = ValNodeNew (NULL);
        vnp->choice = StructuredCommentField_named;
        vnp->data.ptrvalue = StringSave (ufp->label->str);
        ValNodeAddPointer ((ValNodePtr PNTR) data, FieldType_struc_comment_field, vnp);
      }
      ufp = ufp->next;
    }
  }
}


static ValNodePtr GetStructuredCommentFieldList (SeqEntryPtr sep)
{
  ValNodePtr field_list = NULL;
  ValNodePtr dbname, field_name;

  dbname = ValNodeNew (NULL);
  dbname->choice = StructuredCommentField_database;
  ValNodeAddPointer (&field_list, FieldType_struc_comment_field, dbname);
  
  field_name = ValNodeNew (NULL);
  field_name->choice = StructuredCommentField_field_name;
  ValNodeAddPointer (&field_list, FieldType_struc_comment_field, field_name);

  VisitDescriptorsInSep (sep, &field_list, GetStructuredCommentFieldsCallback);

  field_list = ValNodeSort (field_list, SortVnpByFieldType);
  ValNodeUnique (&field_list, SortVnpByFieldType, FieldTypeListFree);
  return field_list;
}


static void CollectBioSourceDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}

static void CollectBioSourceFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void CollectFeaturesCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && data != NULL && sfp->data.choice != SEQFEAT_BIOSRC && sfp->data.choice != SEQFEAT_PUB) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void CollectPubDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_pub && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}

static void CollectPubFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void CollectBioseqCallback (BioseqPtr bsp, Pointer data)
{
  if (bsp != NULL && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void CollectNucBioseqCallback (BioseqPtr bsp, Pointer data)
{
  if (bsp != NULL && data != NULL && !ISA_aa (bsp->mol)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void AddCommentDescriptorDestinationsForBioseq (BioseqPtr bsp, ValNodePtr PNTR dest_list)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  Boolean found = FALSE;
  ObjValNodePtr ovp;

  if (bsp == NULL || dest_list == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &context)) {
    ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
    found = TRUE;
  }
  if (!found) {
    /* if no existing comment descriptor, create one, marked for delete.
     * unmark it for deletion when it gets populated.
     */
    sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_comment);
    sdp->data.ptrvalue = StringSave ("");
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.deleteme = TRUE;
    ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr CollectCommentDescriptors (SeqEntryPtr sep)
{
  ValNodePtr seq_list = NULL, vnp, desc_list = NULL;

  if (sep == NULL) {
    return NULL;
  }

  VisitBioseqsInSep (sep, &seq_list, CollectNucBioseqCallback);

  for (vnp = seq_list; vnp != NULL; vnp = vnp->next) {
    AddCommentDescriptorDestinationsForBioseq (vnp->data.ptrvalue, &desc_list);
  }
  seq_list = ValNodeFree (seq_list);
  return desc_list;
}


static void CollectStructuredCommentsCallback (SeqDescrPtr sdp, Pointer data)
{
  UserObjectPtr uop;

  if (sdp != NULL && data != NULL && sdp->choice == Seq_descr_user
      && (uop = sdp->data.ptrvalue) != NULL
      && IsUserObjectStructuredComment (uop)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


NLM_EXTERN ValNodePtr GetObjectListForFieldType (Uint1 field_type, SeqEntryPtr sep)
{
  ValNodePtr object_list = NULL;
  Uint2      entityID;

  switch (field_type) {
    case FieldType_source_qual:
      VisitDescriptorsInSep (sep, &object_list, CollectBioSourceDescCallback);
      VisitFeaturesInSep (sep, &object_list, CollectBioSourceFeatCallback);
      break;
    case FieldType_cds_gene_prot:
      entityID = ObjMgrGetEntityIDForChoice(sep);
      object_list = BuildCGPSetList (entityID, NULL);
      break;
    case FieldType_feature_field:
      VisitFeaturesInSep (sep, &object_list, CollectFeaturesCallback);
      break;
    case FieldType_molinfo_field:
      VisitBioseqsInSep (sep, &object_list, CollectBioseqCallback);
      break;
    case FieldType_pub:
      VisitDescriptorsInSep (sep, &object_list, CollectPubDescCallback);
      VisitFeaturesInSep (sep, &object_list, CollectPubFeatCallback);
      break;
    case FieldType_rna_field:
      VisitFeaturesInSep (sep, &object_list, CollectFeaturesCallback);
      break;
    case FieldType_struc_comment_field:
      VisitDescriptorsInSep (sep, &object_list, CollectStructuredCommentsCallback);
      break;
    case FieldType_misc:
      VisitBioseqsInSep (sep, &object_list, CollectNucBioseqCallback);
      ValNodeLink (&object_list, CollectCommentDescriptors (sep));
      break;
  }
  return object_list;
}


NLM_EXTERN ValNodePtr GetFieldListForFieldType (Uint1 field_type, SeqEntryPtr sep)
{
  ValNodePtr fields = NULL;

  /* get a list of the fields that are appropriate for the objects collected */
  switch (field_type) {
    case FieldType_cds_gene_prot:
      fields = MakeCDSGeneProtFieldTypeList ();
      break;
    case FieldType_source_qual:
      fields = GetSourceQualSampleFieldList (sep);
      break;
    case FieldType_feature_field:
      fields = GetFeatureQualFieldList (sep);
      break;
    case FieldType_molinfo_field:
      fields = MakeSequenceQualFieldTypeList ();
      break;
    case FieldType_pub:
      fields = MakePubFieldTypeList ();
      break;
    case FieldType_rna_field:
      fields = GetRnaQualFieldList (sep);
      break;
    case FieldType_struc_comment_field:
      fields = GetStructuredCommentFieldList (sep);
      break;
    case FieldType_misc:
      ValNodeAddInt (&fields, FieldType_misc, Misc_field_genome_project_id);
      ValNodeAddInt (&fields, FieldType_misc, Misc_field_comment_descriptor);
      break;
  }
  return fields;
}


NLM_EXTERN ValNodePtr GetAECRSampleListForSeqEntry (Uint1 field_type, SeqEntryPtr sep)
{
  ValNodePtr          object_list;
  ValNodePtr          fields = NULL, vnp;
  ValNodePtr          list = NULL;
  AECRSamplePtr       sample;
  BatchExtraPtr       batch_extra;

  object_list = GetObjectListForFieldType (field_type, sep);

  /* get a list of the fields that are appropriate for the objects collected */
  fields = GetFieldListForFieldType (field_type, sep);

  batch_extra = BatchExtraNew ();
  for (vnp = fields; vnp != NULL; vnp = vnp->next) {
    InitBatchExtraForField (batch_extra, vnp, sep);
  }
  for (vnp = fields; vnp != NULL; vnp = vnp->next) {
    sample = GetAECRSampleFromObjectListEx (object_list, vnp, batch_extra);
    if (sample != NULL && sample->num_found > 0) {
      ValNodeAddPointer (&list, 0, sample);
    } else {
      sample = AECRSampleFree (sample);
    }
  }

  batch_extra = BatchExtraFree (batch_extra);
  fields = FieldTypeListFree (fields);

  object_list = FreeObjectList (object_list);
  return list;
}


NLM_EXTERN ValNodePtr GetAECRSampleList (AECRActionPtr act, SeqEntryPtr sep)
{
  Uint1               field_type;
  Uint2               entityID;
  ValNodePtr          object_list;
  ValNodePtr          fields = NULL, vnp;
  ValNodePtr          list = NULL;
  AECRSamplePtr       sample;
  BatchExtraPtr       batch_extra;

  batch_extra = BatchExtraNew ();
  InitBatchExtraForAECRAction (batch_extra, act, sep);

  field_type = FieldTypeFromAECRAction (act);
  if (field_type == FieldType_cds_gene_prot) {
    entityID = ObjMgrGetEntityIDForChoice(sep);
    object_list = BuildCGPSetList (entityID, act);
  } else {
    object_list = GetObjectListForAECRActionEx (sep, act, batch_extra);
  }

  /* get fields used in action */
  fields = GetFieldTypeListFromAECRAction (act);

  for (vnp = fields; vnp != NULL; vnp = vnp->next) {
    sample = GetAECRSampleFromObjectListEx (object_list, vnp, batch_extra);
    if (sample != NULL && sample->num_found > 0) {
      ValNodeAddPointer (&list, 0, sample);
    } else {
      sample = AECRSampleFree (sample);
    }
  }

  fields = FieldTypeListFree (fields);

  batch_extra = BatchExtraFree (batch_extra);

  FreeObjectList (object_list);
  return list;
}


NLM_EXTERN AECRSamplePtr GetFieldSampleFromList (ValNodePtr list, FieldTypePtr field)
{
  AECRSamplePtr sample = NULL;

  while (list != NULL && sample == NULL) {
    sample = list->data.ptrvalue;
    if (sample != NULL && !DoFieldTypesMatch (sample->field, field)) {
      sample = NULL;
    }
    list = list->next;
  }
  return sample;
}


static void RemoveFieldsForWhichThereAreNoData (ValNodePtr PNTR field_list, ValNodePtr object_list)
{
  ValNodePtr vnp_prev = NULL, vnp_f, vnp_next;
  AECRSamplePtr       sample;

  if (field_list == NULL || *field_list == NULL) {
    return;
  }

  vnp_prev = NULL;
  vnp_f = *field_list;
  while (vnp_f != NULL) {
    vnp_next = vnp_f->next;
    if (vnp_f->choice == FieldType_source_qual
        || vnp_f->choice == FieldType_feature_field
        || vnp_f->choice == FieldType_rna_field) {
      vnp_prev = vnp_f;
    } else {
      sample = GetAECRSampleFromObjectList (object_list, vnp_f);
      if (sample == NULL || sample->num_found == 0) {
        if (vnp_prev == NULL) {
          *field_list = vnp_next;
        } else {
          vnp_prev->next = vnp_next;
        }
        vnp_f->next = NULL;
        vnp_f = FieldTypeFree (vnp_f);
      } else {
        vnp_prev = vnp_f;
      }
      sample = AECRSampleFree (sample);
    }
    vnp_f = vnp_next;
  }
}


NLM_EXTERN void GetAECRExistingTextList (Uint1 field_type, SeqEntryPtr sep, FILE *fp)
{
  ValNodePtr          object_list, vnp_f, vnp_o;
  ValNodePtr          fields = NULL;
  BioseqPtr           bsp;
  Char                id_buf[255];
  CharPtr             txt1 = NULL;

  object_list = GetObjectListForFieldType (field_type, sep);

  /* get a list of the fields that are appropriate for the objects collected */
  fields = GetFieldListForFieldType (field_type, sep);

  /* remove fields for which there is no data */
  RemoveFieldsForWhichThereAreNoData (&fields, object_list);

  /* add header */
  fprintf (fp, "Accession");
  for (vnp_f = fields; vnp_f != NULL; vnp_f = vnp_f->next) {
    txt1 = SummarizeFieldType (vnp_f);
    fprintf (fp, "\t%s", txt1);
    txt1 = MemFree (txt1);
  } 
  fprintf (fp, "\n");
    
  for (vnp_o = object_list; vnp_o != NULL; vnp_o = vnp_o->next) {
    bsp = GetSequenceForObject (vnp_o->choice, vnp_o->data.ptrvalue);
    if (bsp == NULL) {
      id_buf[0] = 0;
    } else {
      SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
    }
    fprintf (fp, "%s", id_buf);
    for (vnp_f = fields; vnp_f != NULL; vnp_f = vnp_f->next) {
      txt1 = GetFieldValueForObject (vnp_o->choice, vnp_o->data.ptrvalue, vnp_f, NULL);
      fprintf (fp, "\t%s", txt1 == NULL ? "" : txt1);
      txt1 = MemFree (txt1);
    }
    fprintf (fp, "\n");
  }

  fields = FieldTypeListFree (fields);

  object_list = FreeObjectList (object_list);
}


NLM_EXTERN void ExportFieldTable (Uint1 field_type, ValNodePtr src_field_list, SeqEntryPtr sep, FILE *fp)
{
  ValNodePtr          object_list, vnp_f, vnp_o;
  ValNodePtr          fields = NULL;
  BioseqPtr           bsp;
  Char                id_buf[255];
  CharPtr             txt1 = NULL;
  SeqDescrPtr         sdp;
  SeqMgrDescContext   context;

  if (field_type == 0) {
    object_list = GetObjectListForFieldType (FieldType_source_qual, sep);
  } else {
    object_list = GetObjectListForFieldType (field_type, sep);
    /* get a list of the fields that are appropriate for the objects collected */
    fields = GetFieldListForFieldType (field_type, sep);

    /* remove fields for which there is no data */
    RemoveFieldsForWhichThereAreNoData (&fields, object_list);
  }

  /* add header */
  /* accession is first column */
  fprintf (fp, "Accession");
  /* list source fields first */
  for (vnp_f = src_field_list; vnp_f != NULL; vnp_f = vnp_f->next) {
    txt1 = SummarizeFieldType (vnp_f);
    fprintf (fp, "\t%s", txt1);
    txt1 = MemFree (txt1);
  } 
  /* list fields */
  for (vnp_f = fields; vnp_f != NULL; vnp_f = vnp_f->next) {
    txt1 = SummarizeFieldType (vnp_f);
    fprintf (fp, "\t%s", txt1);
    txt1 = MemFree (txt1);
  } 
  fprintf (fp, "\n");
    
  for (vnp_o = object_list; vnp_o != NULL; vnp_o = vnp_o->next) {
    bsp = GetSequenceForObject (vnp_o->choice, vnp_o->data.ptrvalue);
    if (bsp == NULL) {
      id_buf[0] = 0;
    } else {
      SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
    }
    /* print accession */
    fprintf (fp, "%s", id_buf);
    /* print source fields */
    if (src_field_list != NULL) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
      for (vnp_f = src_field_list; vnp_f != NULL; vnp_f = vnp_f->next) {
        txt1 = GetFieldValueForObject (OBJ_SEQDESC, sdp, vnp_f, NULL);
        fprintf (fp, "\t%s", txt1 == NULL ? "" : txt1);
        txt1 = MemFree (txt1);
      }
    }
    /* get requested fields */
    for (vnp_f = fields; vnp_f != NULL; vnp_f = vnp_f->next) {
      txt1 = GetFieldValueForObject (vnp_o->choice, vnp_o->data.ptrvalue, vnp_f, NULL);
      fprintf (fp, "\t%s", txt1 == NULL ? "" : txt1);
      txt1 = MemFree (txt1);
    }
    fprintf (fp, "\n");
  }

  fields = FieldTypeListFree (fields);

  object_list = FreeObjectList (object_list);
}


/* This section handles parsing where the source field and destination field may not be on the same
 * group of objects. */
typedef struct parsesourceinfo 
{
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
  SeqDescrPtr sdp;
  SeqIdPtr    sip;
  ValNodePtr  dest_list;
  CharPtr     parse_src_txt;
} ParseSourceInfoData, PNTR ParseSourceInfoPtr;

static ParseSourceInfoPtr ParseSourceInfoNew (BioseqPtr bsp, SeqFeatPtr sfp, SeqDescrPtr sdp, SeqIdPtr sip, CharPtr parse_src_txt)
{
  ParseSourceInfoPtr psip;

  psip = (ParseSourceInfoPtr) MemNew (sizeof (ParseSourceInfoData));
  if (psip != NULL) {
    psip->bsp = bsp;
    psip->sdp = sdp;
    psip->sfp = sfp;
    psip->sip = sip;
    psip->dest_list = NULL;
    psip->parse_src_txt = parse_src_txt;
  } 
  return psip;
}


static ParseSourceInfoPtr ParseSourceInfoFree (ParseSourceInfoPtr psip)
{
  if (psip != NULL)
  {
    psip->dest_list = ValNodeFree (psip->dest_list);
    psip->parse_src_txt = MemFree (psip->parse_src_txt);
    psip = MemFree (psip);
  }
  return psip;
}

static ParseSourceInfoPtr ParseSourceInfoCopy (ParseSourceInfoPtr psip)
{
  ParseSourceInfoPtr pcopy = NULL;
  
  if (psip != NULL) 
  {
    pcopy = (ParseSourceInfoPtr) MemNew (sizeof (ParseSourceInfoData));
    if (pcopy != NULL) {
      pcopy->bsp = psip->bsp;
      pcopy->sfp = psip->sfp;
      pcopy->sdp = psip->sdp;
      pcopy->sip = psip->sip;
      pcopy->dest_list = NULL;
      pcopy->parse_src_txt = NULL;
    }
  }
  return pcopy;
}

static ValNodePtr ParseSourceListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = ParseSourceInfoFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void 
GetDeflineSourcesForBioseq 
(BioseqPtr              bsp,
 TextPortionPtr         portion,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  CharPtr            str;
  ParseSourceInfoPtr psip;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  while (sdp != NULL)
  {
    str = GetTextPortionFromString (sdp->data.ptrvalue, portion);    
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
      if (psip != NULL) {
        ValNodeAddPointer (source_list, 0, psip);
      } else {
        str = MemFree (str);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_title, &dcontext);
  }
}


static CharPtr GetIDSrc (SeqIdPtr sip, Uint1 id_type, CharPtr tag)
{
  DbtagPtr    dbt = NULL;
  ObjectIdPtr oip = NULL;
  Char        id_str[128];
  CharPtr     str_src = NULL;

  if (sip == NULL || sip->choice != id_type) return NULL;

  if (id_type == SEQID_GENERAL)
  {
    dbt = (DbtagPtr) sip->data.ptrvalue;
    if (dbt == NULL || (tag != NULL && StringCmp (dbt->db, tag) != 0)) return NULL;
    oip = dbt->tag;
  }
  else if (id_type == SEQID_LOCAL)
  {
    oip = sip->data.ptrvalue;
  }

  if (oip == NULL)
  {
    SeqIdWrite (sip, id_str, PRINTID_REPORT, sizeof (id_str));
    str_src = StringSave (id_str);
  }
  else
  {
    if (oip->str == NULL)
    {
      sprintf (id_str, "%d", oip->id);
      str_src = StringSave (id_str);
    }
    else
    {
      str_src = StringSave (oip->str);
    }
  }
  return str_src;
}


static void
GetIDSourcesForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  portion,
 Uint1           id_type,
 CharPtr         tag,
 ValNodePtr PNTR source_list)
{
  SeqIdPtr           sip;
  ParseSourceInfoPtr psip;
  CharPtr            src_str = NULL, str;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sip = bsp->id;
  while (sip != NULL)
  {
    if ((src_str = GetIDSrc (sip, id_type, tag)) != NULL) { 
      str = GetTextPortionFromString (src_str, portion); 
      if (str != NULL) {
        psip = ParseSourceInfoNew (bsp, NULL, NULL, sip, str);
        if (psip != NULL) {
          ValNodeAddPointer (source_list, 0, psip);
        } else {
          str = MemFree (str);
        }
      }
      src_str = MemFree (src_str);
    }
    sip = sip->next;
  }
}


static void
GetLocalIDSourcesForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  GetIDSourcesForBioseq (bsp, tp, SEQID_LOCAL, NULL, source_list);
}


static void GetNcbiFileSourceForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  GetIDSourcesForBioseq (bsp, tp, SEQID_GENERAL, "NCBIFILE", source_list);
}


static void StripBankitCommentForParse (SeqDescrPtr sdp, TextPortionPtr tp)
{
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  
  if (sdp == NULL || sdp->choice != Seq_descr_user || tp == NULL) {
    return;
  }
  
  /* Bankit Comments */
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
      for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
        oip = ufp->label;
        if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
          ReplaceStringForParse (ufp->data.ptrvalue, tp);
        }
      }
    }
  }
}


static void StripStructuredCommentForParse (SeqDescrPtr sdp, CharPtr comment_field, TextPortionPtr tp)
{
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;

  if (sdp == NULL || sdp->choice != Seq_descr_user || tp == NULL || StringHasNoText (comment_field)) {
    return;
  }
    
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (IsUserObjectStructuredComment (uop)) {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip != NULL && StringCmp (oip->str, comment_field) == 0) {
        ReplaceStringForParse (ufp->data.ptrvalue, tp);
      }
    }
  }
}


static void
GetBankitCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  ParseSourceInfoPtr psip;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  CharPtr            str = NULL;
  
  if (bsp == NULL || source_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    if (sdp->extended != 0) {
      /* Bankit Comments */
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
        oip = uop->type;
        if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
              str = GetTextPortionFromString (ufp->data.ptrvalue, tp);
              if (str != NULL) {
                psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
                if (psip == NULL) {
                  str = MemFree (str);
                } else {
                  ValNodeAddPointer (source_list, 0, psip);
                }
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
}


static void 
GetCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;
  ParseSourceInfoPtr psip;
  CharPtr            str;
  
  if (bsp == NULL || source_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    str = GetTextPortionFromString (sdp->data.ptrvalue, tp);
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
      if (psip == NULL) {
        str = MemFree (str);
      } else {
        ValNodeAddPointer (source_list, 0, psip);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_COMMENT, 0, &fcontext);
  while (sfp != NULL) {
    str = GetTextPortionFromString (sfp->data.value.ptrvalue, tp);
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, sfp, NULL, NULL, str);
      if (psip == NULL) {
        str = MemFree (str);
      } else {
        ValNodeAddPointer (source_list, 0, psip);
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_COMMENT, 0, &fcontext);
  }
  GetBankitCommentSourcesForBioseq (bsp, tp, source_list);
}


static void 
GetStructuredCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 CharPtr         comment_field,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  SeqMgrDescContext  dcontext;
  CharPtr            str;
  ParseSourceInfoPtr psip;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {  
    if (sdp->extended != 0
        && sdp->data.ptrvalue != NULL) {
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      if (IsUserObjectStructuredComment (uop)) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip != NULL && StringCmp (oip->str, comment_field) == 0) {
            str = GetTextPortionFromString (ufp->data.ptrvalue, tp);
            if (str != NULL) {
              psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
              if (psip == NULL) {
                str = MemFree (str);
              } else {
                ValNodeAddPointer (source_list, 0, psip);
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
}


const CharPtr nomial_keywords[] = {
"f. sp. ",
"var.",
"pv.",
"bv.",
"serovar",
"subsp." };

const Int4 num_nomial_keywords = sizeof(nomial_keywords) / sizeof (CharPtr);

static CharPtr GetTextAfterNomial (CharPtr taxname)

{
  CharPtr ptr, nomial_end;
  Int4    i;
  Boolean found_keyword = TRUE;
  
  ptr = StringChr (taxname, ' ');
  if (ptr == NULL) return NULL;
  /* skip over the first word and the spaces after it. */
  while (*ptr == ' ') {
    ptr++;
  }
  ptr = StringChr (ptr, ' ');
  /* if there are only two words, give up. */
  if (ptr == NULL) {
    return NULL;
  }
  nomial_end = ptr;
  while (*ptr == ' ') {
    ptr++;
  }
  
  while (found_keyword) {
    found_keyword = FALSE;
    /* if the next word is a nomial keyword, skip that plus the first word that follows it. */
    for (i = 0; i < num_nomial_keywords && *nomial_end != 0; i++) {
      if (StringNCmp (ptr, nomial_keywords[i], StringLen(nomial_keywords[i])) == 0) {
        ptr += StringLen(nomial_keywords[i]);
        while (*ptr == ' ' ) {
          ptr++;
        }
        nomial_end = StringChr (ptr, ' ');
        if (nomial_end == NULL) {
          nomial_end = ptr + StringLen (ptr);
        } else {          
          ptr = nomial_end;
          while (*ptr == ' ') {
            ptr++;
          }
          found_keyword = TRUE;
        }
      }
    }
  }
  return nomial_end;
}


static void 
GetOrgParseSourcesForBioSource 
(BioSourcePtr    biop,
 BioseqPtr       bsp,
 SeqDescrPtr     sdp,
 SeqFeatPtr      sfp,
 ParseSrcOrgPtr  o,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  CharPtr str = NULL, portion, tmp;
  ValNode vn;
  ParseSourceInfoPtr psip;

  if (biop == NULL || o == NULL || o->field == NULL || source_list == NULL) return;

  switch (o->field->choice) {
    case ParseSrcOrgChoice_source_qual :
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = o->field->data.intvalue;
      vn.next = NULL;
      str = GetSourceQualFromBioSource (biop, &vn, NULL);
      break;
    case ParseSrcOrgChoice_taxname_after_binomial :
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = Source_qual_taxname;
      vn.next = NULL;
      str = GetSourceQualFromBioSource (biop, &vn, NULL);
      tmp = GetTextAfterNomial (str);
      tmp = StringSave (tmp);
      str = MemFree (str);
      str = tmp;
      break;
  }
  portion = GetTextPortionFromString (str, tp);
  if (portion != NULL) {
    psip = ParseSourceInfoNew (bsp, sfp, sdp, NULL, portion);
    if (psip == NULL) {
      portion = MemFree (portion);
    } else {
      ValNodeAddPointer (source_list, 0, psip);
    }
  }
  str = MemFree (str);
}


static void GetOrgParseSourcesForBioseq (BioseqPtr bsp, ParseSrcOrgPtr o, TextPortionPtr tp, ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || o == NULL || source_list == NULL) return;

  if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor) {
    for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
         sdp != NULL;
         sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
      GetOrgParseSourcesForBioSource (sdp->data.ptrvalue, bsp, sdp, NULL, o, tp, source_list);
    }
  }

  if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext)) {
      GetOrgParseSourcesForBioSource (sfp->data.value.ptrvalue, bsp, NULL, sfp, o, tp, source_list);
    }
  }
}


typedef struct parsesrccollection {
  ParseSrcPtr src;
  TextPortionPtr portion;
  ValNodePtr src_list;
} ParseSrcCollectionData, PNTR ParseSrcCollectionPtr;


static void FindParseSourceBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  ParseSrcCollectionPtr psp;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  psp = (ParseSrcCollectionPtr) userdata;
  if (psp->src == NULL) return;

  switch (psp->src->choice)
  {
    case ParseSrc_defline:
      if (!ISA_aa (bsp->mol)) {
        GetDeflineSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
    case ParseSrc_local_id:
      if (! ISA_aa (bsp->mol) && bsp->repr != Seq_repr_seg) {
        GetLocalIDSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
    case ParseSrc_file_id:
      GetNcbiFileSourceForBioseq (bsp, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_org:
      GetOrgParseSourcesForBioseq (bsp, psp->src->data.ptrvalue, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_comment:
      GetCommentSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_structured_comment:
      GetStructuredCommentSourcesForBioseq(bsp, psp->portion, psp->src->data.ptrvalue, &(psp->src_list));
      break;
    case ParseSrc_bankit_comment:
      if (!ISA_aa (bsp->mol)) {
        GetBankitCommentSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
  }
}


static void GetOrgNamesInRecordCallback (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr PNTR org_names;
  
  if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)
      || userdata == NULL)
  {
    return;
  }
  
  org_names = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (org_names, 0, biop->org->taxname);
}


static void SetToUpper (CharPtr cp)
{
  if (cp == NULL) return;
  while (*cp != 0) {
    if (isalpha (*cp)) {
      *cp = toupper (*cp);
    }
    cp++;
  }
}


static void 
FixCapitalizationInString 
(CharPtr PNTR pTitle,
 Uint2 capitalization,
 ValNodePtr   org_names)
{
  if (pTitle == NULL || capitalization == Cap_change_none) return;

  switch (capitalization) {
    case Cap_change_tolower:
      ResetCapitalization (FALSE, *pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
    case Cap_change_toupper:
      SetToUpper (*pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
    case Cap_change_firstcap:
      ResetCapitalization (TRUE, *pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
  }
}


static void AddDeflineDestinationsForBioseq (BioseqPtr bsp, ValNodePtr PNTR dest_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || dest_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  while (sdp != NULL) {
    ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_title, &dcontext);
  }
}


static ValNodePtr GetFeatureListForNucleotideBioseq (Uint1 featdef, BioseqPtr bsp);
static ValNodePtr GetFeatureListForProteinBioseq (Uint1 featdef, BioseqPtr bsp);

static void AddFeatureDestinationsForBioseq (BioseqPtr bsp, FeatureFieldLegalPtr featfield, ValNodePtr PNTR dest_list)
{
  Int4 featdef;

  if (bsp == NULL || featfield == NULL || dest_list == NULL) return;

  featdef = GetFeatdefFromFeatureType (featfield->type);
  if (ISA_aa (bsp->mol)) {
    ValNodeLink (dest_list, GetFeatureListForProteinBioseq (featdef, bsp));
  } else {
    ValNodeLink (dest_list, GetFeatureListForNucleotideBioseq (featdef, bsp));
  }

}


static void GetBioSourceDestinationsForBioseq (BioseqPtr bsp, Uint2 object_type, ValNodePtr PNTR dest_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || dest_list == NULL)
  {
    return;
  }
  
  if (object_type == Object_type_constraint_any || object_type == Object_type_constraint_descriptor) 
  {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    while (sdp != NULL)
    {
      ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
      sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
    }
  }
  
  if (object_type == Object_type_constraint_any || object_type == Object_type_constraint_feature)
  {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    while (sfp != NULL)
    {
      ValNodeAddPointer (dest_list, OBJ_SEQFEAT, sfp);
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
    }  
  }
}


static void AddParseDestinations (ParseSourceInfoPtr psip, ParseDestPtr dst)
{
  ParseDstOrgPtr o;

  if (psip == NULL || dst == NULL) return;

  switch (dst->choice) {
    case ParseDest_defline :
      AddDeflineDestinationsForBioseq (psip->bsp, &(psip->dest_list));
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) dst->data.ptrvalue;
      if ((o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor)
          && psip->sdp != NULL && psip->sdp->choice == Seq_descr_source) {
        ValNodeAddPointer (&(psip->dest_list), OBJ_SEQDESC, psip->sdp);
      } else if ((o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature)
                 && psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_BIOSRC) {
        ValNodeAddPointer (&(psip->dest_list), OBJ_SEQFEAT, psip->sfp);
      } else {
        GetBioSourceDestinationsForBioseq (psip->bsp, o->type, &(psip->dest_list));
      }
      break;
    case ParseDest_featqual :
      AddFeatureDestinationsForBioseq (psip->bsp, dst->data.ptrvalue, &(psip->dest_list));
      break;
    case ParseDest_comment_descriptor :
      AddCommentDescriptorDestinationsForBioseq (psip->bsp, &(psip->dest_list));
      break;
    case ParseDest_dbxref :
      GetBioSourceDestinationsForBioseq (psip->bsp, Object_type_constraint_any, &(psip->dest_list));
      break;
  }
}


static Boolean SourceHasOneUndeletedDestination (ParseSourceInfoPtr source)
{
  Int4       num_seen = 0;
  ValNodePtr vnp;
  
  if (source == NULL
      || source->dest_list == NULL)
  {
    return FALSE;
  }
  
  vnp = source->dest_list;
  while (vnp != NULL && num_seen < 2)
  {
    if (vnp->choice > 1)
    {
      num_seen ++;
    }
    vnp = vnp->next;
  }
  if (num_seen == 1)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void CombineSourcesForDestinations (ValNodePtr PNTR source_list)
{
  ValNodePtr         source1_vnp, source2_vnp, dest1_vnp, dest2_vnp;
  ValNodePtr         source_new, del_vnp;
  ParseSourceInfoPtr psip1, psip2, new_psip;
  CharPtr            comb_txt;
  
  for (source1_vnp = *source_list;
       source1_vnp != NULL; 
       source1_vnp = source1_vnp->next)
  {
    psip1 = (ParseSourceInfoPtr) source1_vnp->data.ptrvalue;
    if (psip1 == NULL || psip1->dest_list == NULL)
    {
      continue;
    }
    for (source2_vnp = source1_vnp->next;
         source2_vnp != NULL; 
         source2_vnp = source2_vnp->next)
    {
      if (source2_vnp->choice > 0) 
      {
        /* already marked for deletion */
        continue;
      }
      psip2 = (ParseSourceInfoPtr) source2_vnp->data.ptrvalue;
      if (psip2 == NULL || psip2->dest_list == NULL)
      {
        continue;
      }
      for (dest1_vnp = psip1->dest_list;
           dest1_vnp != NULL; 
           dest1_vnp = dest1_vnp->next)
      {
        if (dest1_vnp->choice == 0)
        {
          /* already marked for deletion */
          continue;
        }
        for (dest2_vnp = psip2->dest_list;
             dest2_vnp != NULL;
             dest2_vnp = dest2_vnp->next)
        {
          if (dest2_vnp->choice == 0)
          {
            /* already marked for deletion */
            continue;
          }
          if (dest1_vnp->choice == dest2_vnp->choice
              && dest1_vnp->data.ptrvalue == dest2_vnp->data.ptrvalue)
          {
            comb_txt = (CharPtr) (MemNew (sizeof (Char) 
                                  * (StringLen (psip1->parse_src_txt)
                                     + StringLen (psip2->parse_src_txt)
                                     + 2)));
            StringCpy (comb_txt, psip1->parse_src_txt);
            StringCat (comb_txt, ";");
            StringCat (comb_txt, psip2->parse_src_txt);
            
            /* If the first source has a single destination, then we can 
             * add the text from the second source to the first and remove
             * the destination from the second source.
             */
            if (SourceHasOneUndeletedDestination (psip1))
            {
              
              psip1->parse_src_txt = MemFree (psip1->parse_src_txt);
              psip1->parse_src_txt = comb_txt;
              dest2_vnp->choice = 0;
            }             
            /* If the first source has more than one destination and
             * the second source has a single destination, then we can 
             * remove the repeated desination from the first source
             * and add the text from the first source to the second source.
             */
            else if (SourceHasOneUndeletedDestination (psip2))
            {
              psip2->parse_src_txt = MemFree (psip2->parse_src_txt);
              psip2->parse_src_txt = comb_txt;
              dest1_vnp->choice = 0;
            }
            /* If the first and second sources have multiple destinations,
             * we need to remove the repeated destination from both the first
             * and second source and create a new source with the combined 
             * text for just the repeated destination.
             */
            else
            {
              new_psip = ParseSourceInfoNew (NULL, NULL, NULL, NULL, comb_txt);
              ValNodeAddPointer (&(new_psip->dest_list), 
                                 dest1_vnp->choice, 
                                 dest1_vnp->data.ptrvalue);
              dest1_vnp->choice = 0;
              dest2_vnp->choice = 0;
              source_new = ValNodeNew (NULL);
              source_new->choice = 0;
              source_new->data.ptrvalue = new_psip;
              source_new->next = source1_vnp->next;
              source1_vnp->next = source_new;
            }
          }
        }
      }
      
      del_vnp = ValNodeExtractList (&(psip1->dest_list), 0);
      del_vnp = ValNodeFree (del_vnp);
      if (psip1->dest_list == NULL)
      {
        source1_vnp->choice = 1;
      }
      del_vnp = ValNodeExtractList (&(psip2->dest_list), 0);
      del_vnp = ValNodeFree (del_vnp);
      if (psip2->dest_list == NULL)
      {
        source2_vnp->choice = 1;
      }
    }
  }

  /* now remove sources deleted */
  del_vnp = ValNodeExtractList (source_list, 1);
  del_vnp = ParseSourceListFree (del_vnp); 
}


static BioseqSetPtr GetPartsForSourceDescriptorOnSegSet (SeqDescrPtr sdp)
{
  ObjValNodePtr ovp;
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;
  
  if (sdp == NULL || sdp->extended != 1) {
    return NULL;
  }
  ovp = (ObjValNodePtr) sdp;
  if (ovp->idx.parenttype != OBJ_BIOSEQSET || ovp->idx.parentptr == NULL) {
    return NULL;
  }
  bssp = (BioseqSetPtr) ovp->idx.parentptr;
  
  if (bssp->_class == BioseqseqSet_class_nuc_prot
      && IS_Bioseq_set (bssp->seq_set)
      && bssp->seq_set->data.ptrvalue != NULL) {
    bssp = (BioseqSetPtr) bssp->seq_set->data.ptrvalue;
  }
  
  if (bssp->_class == BioseqseqSet_class_segset) {
    sep = bssp->seq_set;
    while (sep != NULL) {
      if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp->_class == BioseqseqSet_class_parts) {
          return bssp;
        }
      }
      sep = sep->next;
    }
  }

  return NULL;
}


static SeqDescrPtr FindSourceDescriptorInSeqEntry (SeqEntryPtr sep)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp = NULL;
  
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sdp = bsp->descr;
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      sdp = bssp->descr;
    }
    while (sdp != NULL && sdp->choice != Seq_descr_source)
    {
      sdp = sdp->next;
    }
  }
  return sdp;
}


static SeqDescrPtr PropagateToSeqEntry (SeqEntryPtr sep, SeqDescrPtr sdp)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  new_sdp = NULL;
  
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      new_sdp = AsnIoMemCopy ((Pointer) sdp,
                              (AsnReadFunc) SeqDescrAsnRead,
                              (AsnWriteFunc) SeqDescrAsnWrite);
      ValNodeLink (&(bsp->descr), new_sdp);
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      new_sdp = AsnIoMemCopy ((Pointer) sdp,
                              (AsnReadFunc) SeqDescrAsnRead,
                              (AsnWriteFunc) SeqDescrAsnWrite);
      ValNodeLink (&(bssp->descr), new_sdp);
    }
  }
  return new_sdp;
}


static void PropagateSourceOnSegSetForParse (ValNodePtr parse_source_list)
{
  ParseSourceInfoPtr psip;
  ValNodePtr         vnp_src, vnp_dst;
  SeqDescrPtr        sdp, other_sdp;
  SeqEntryPtr        sep;
  ValNodePtr         extra_dests = NULL;
  BioseqSetPtr       parts_bssp;
  
  for (vnp_src = parse_source_list; vnp_src != NULL; vnp_src = vnp_src->next) {
    psip = (ParseSourceInfoPtr) vnp_src->data.ptrvalue;
    if (psip != NULL) {
      for (vnp_dst = psip->dest_list; vnp_dst != NULL; vnp_dst = vnp_dst->next) {
        if (vnp_dst->choice == OBJ_SEQDESC) {
          sdp = (SeqDescrPtr) vnp_dst->data.ptrvalue;
          if (sdp != NULL && sdp->choice == Seq_descr_source) {
            parts_bssp = GetPartsForSourceDescriptorOnSegSet (sdp);
            if (parts_bssp != NULL) {
              for (sep = parts_bssp->seq_set; sep != NULL; sep = sep->next) {
                if (IS_Bioseq(sep) && sep->data.ptrvalue == psip->bsp) {
                  other_sdp = FindSourceDescriptorInSeqEntry (sep);
                  if (other_sdp == NULL) {
                    other_sdp = PropagateToSeqEntry (sep, sdp);
                    ValNodeAddPointer (&extra_dests, OBJ_SEQDESC, other_sdp);
                  }
                }
              }
            
              /* set choice to 0 so master won't be a destination */
              vnp_dst->choice = 0;
            
            }
          }
        }
      }
      /* add extra destinations to list */
      ValNodeLink (&psip->dest_list, extra_dests);
      extra_dests = NULL;
    }
  }
  
}



static CharPtr GetDBxrefFromBioSource (BioSourcePtr biop, CharPtr db_name)
{
  CharPtr    rval = NULL;
  ValNodePtr vnp;
  DbtagPtr   dbtag;

  if (biop == NULL || biop->org == NULL || StringHasNoText (db_name)) {
    return NULL;
  }
  for (vnp = biop->org->db; vnp != NULL && rval != NULL; vnp = vnp->next) {
    dbtag = (DbtagPtr) vnp->data.ptrvalue;
    if (dbtag != NULL && StringCmp (db_name, dbtag->db)) {
      rval = GetObjectIdString (dbtag->tag);
    }
  }
  return rval;
}


static Boolean SetDBxrefForBioSource (BioSourcePtr biop, CharPtr db_name, CharPtr str, Uint2 existing_text)
{
  ValNodePtr    dbx;
  DbtagPtr      dbtag;
  Boolean       found = FALSE;
  Char          buf[20];
  Boolean       rval = FALSE;

  if (biop == NULL || StringHasNoText (db_name) || StringHasNoText (str)) {
    return FALSE;
  }

  if (biop->org == NULL)
  {
    biop->org = OrgRefNew();
  }
  dbx = biop->org->db;
  while (dbx != NULL && !found)
  {
    dbtag = (DbtagPtr) dbx->data.ptrvalue;
    if (dbtag != NULL && dbtag->tag != NULL
        && StringCmp (dbtag->db, db_name) == 0)
    {
      found = TRUE;
    }
    if (!found)
    {
      dbx = dbx->next;
    }
  }
  if (!found)
  {
    dbtag = DbtagNew();
    dbtag->db = StringSave (db_name);      
    ValNodeAddPointer (&(biop->org->db), 0, dbtag);
  }
  if (dbtag->tag == NULL)
  {
    dbtag->tag = ObjectIdNew();
  }
  /* if it was a number before, make it a string now */
  if (dbtag->tag->id > 0 && dbtag->tag->str == NULL)
  {
    sprintf (buf, "%d", dbtag->tag->id);
    dbtag->tag->id = 0;
    dbtag->tag->str = StringSave (buf);
  }
  rval = SetStringValue (&(dbtag->tag->str), str, existing_text);
  return rval;
}


static Int4 SetFieldForDestList (ValNodePtr dest_list, ParseDestPtr field, CharPtr str, Uint2 existing_text)
{
  ValNodePtr vnp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  CharPtr     cp;
  BioSourcePtr biop;
  ParseDstOrgPtr o;
  FeatureFieldLegalPtr fl;
  FeatureField f;
  Boolean      was_empty;
  Int4         num_succeeded = 0;

  if (dest_list == NULL || field == NULL) return 0;

  switch (field->choice) {
    case ParseDest_defline :
      for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQDESC && vnp->data.ptrvalue != NULL) {
          sdp = (SeqDescrPtr) vnp->data.ptrvalue;
          if (sdp->choice == Seq_descr_title) {
            cp = sdp->data.ptrvalue;
            if (SetStringValue (&cp, str, existing_text)) {
              num_succeeded++;
            }
            sdp->data.ptrvalue = cp;
          }
        }
      }
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) field->data.ptrvalue;
      if (o != NULL) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          if (SetSourceQualInBioSource (biop, o->field, NULL, str, existing_text)) {
            num_succeeded++;
          }
        }
      }
      break;
    case ParseDest_featqual:
      fl = (FeatureFieldLegalPtr) field->data.ptrvalue;
      if (fl != NULL) {
        f.type = fl->type;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = fl->field;        
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          if (SetQualOnFeature (vnp->data.ptrvalue, &f, NULL, str, existing_text)) {
            num_succeeded++;
          }
        }
        f.field = ValNodeFree (f.field);
      }
      break;
    case ParseDest_comment_descriptor:
      for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
        sdp = vnp->data.ptrvalue;
        if (StringHasNoText (sdp->data.ptrvalue)) {
          was_empty = TRUE;
        } else {
          was_empty = FALSE;
        }
        cp = sdp->data.ptrvalue;
        if (SetStringValue (&cp, str, existing_text)) {
          num_succeeded++;
        }
        sdp->data.ptrvalue = cp;
        if (was_empty) {
          ovp = (ObjValNodePtr) sdp;
          ovp->idx.deleteme = FALSE;
        }
      }
      break;
    case ParseDest_dbxref:
      if (!StringHasNoText (field->data.ptrvalue)) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          if (SetDBxrefForBioSource (biop, field->data.ptrvalue, str, existing_text)) {
            num_succeeded++;
          }
        }
      }
      break;
  }
  return num_succeeded;
}



static void AddToSampleForDestList (AECRSamplePtr sample, ValNodePtr dest_list, ParseDestPtr field)
{
  ValNodePtr vnp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;
  ParseDstOrgPtr o;
  FeatureFieldLegalPtr fl;
  FeatureField f;

  if (dest_list == NULL || field == NULL || sample == NULL) return;

  switch (field->choice) {
    case ParseDest_defline :
      for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQDESC && vnp->data.ptrvalue != NULL) {
          sdp = (SeqDescrPtr) vnp->data.ptrvalue;
          if (sdp->choice == Seq_descr_title) {
            AddTextToAECRSample (sample, StringSave (sdp->data.ptrvalue));
          }
        }
      }
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) field->data.ptrvalue;
      if (o != NULL) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          AddTextToAECRSample (sample, GetSourceQualFromBioSource (biop, o->field, NULL));
        }
      }
      break;
    case ParseDest_featqual:
      fl = (FeatureFieldLegalPtr) field->data.ptrvalue;
      if (fl != NULL) {
        f.type = fl->type;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = fl->field;        
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          AddTextToAECRSample (sample, GetQualFromFeature (vnp->data.ptrvalue, &f, NULL));
        }
        f.field = ValNodeFree (f.field);
      }
      break;
    case ParseDest_comment_descriptor:
      for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
        sdp = (SeqDescrPtr) vnp->data.ptrvalue;
        AddTextToAECRSample (sample, StringSave (sdp->data.ptrvalue));
      }
      break;
    case ParseDest_dbxref:
      if (!StringHasNoText (field->data.ptrvalue)) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          AddTextToAECRSample (sample, GetDBxrefFromBioSource (biop, field->data.ptrvalue));
        }
      }
      break;
  }
}


static void StripFieldForSrcList (ParseSourceInfoPtr psip, ParseSrcPtr field, TextPortionPtr text_portion)
{
  CharPtr     str;
  ParseSrcOrgPtr o;
  BioSourcePtr biop;

  if (psip == NULL || field == NULL || text_portion == NULL) return;

  switch (field->choice) {
    case ParseSrc_defline :
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_title) {
        ReplaceStringForParse (psip->sdp->data.ptrvalue, text_portion);
      }
      break;
    case ParseSrc_org :
      o = (ParseSrcOrgPtr) field->data.ptrvalue;
      if (o != NULL) {
        if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_source) {
          biop = (BioSourcePtr) psip->sdp->data.ptrvalue;
          str = GetSourceQualFromBioSource (biop, o->field, NULL);
          ReplaceStringForParse (str, text_portion);
          SetSourceQualInBioSource (biop, o->field, NULL, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        } else if (psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_BIOSRC) {
          biop = (BioSourcePtr) psip->sfp->data.value.ptrvalue;
          str = GetSourceQualFromBioSource (biop, o->field, NULL);
          ReplaceStringForParse (str, text_portion);
          SetSourceQualInBioSource (biop, o->field, NULL, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        }
      }
      break;
    case ParseSrc_comment:
      if (psip->sdp != NULL) {
        if (psip->sdp->choice == Seq_descr_user) {
          StripBankitCommentForParse (psip->sdp, text_portion);
        } else if (psip->sdp->choice == Seq_descr_comment) {
          ReplaceStringForParse (psip->sdp->data.ptrvalue, text_portion);
        }
      }
      if (psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_COMMENT) {
        ReplaceStringForParse (psip->sfp->data.value.ptrvalue, text_portion);
      }
      break;
    case ParseSrc_bankit_comment:
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_user) {
        StripBankitCommentForParse (psip->sdp, text_portion);
      }
      break;
    case ParseSrc_structured_comment:
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_user) {
        StripStructuredCommentForParse (psip->sdp, field->data.ptrvalue, text_portion);
      }
      break;
  }
}



NLM_EXTERN AECRSamplePtr GetExistingTextForParseAction (ParseActionPtr action, SeqEntryPtr sep)
{
  ParseSrcCollectionData psd;
  ParseSourceInfoPtr     psip;
  ValNodePtr             vnp;
  ValNodePtr             dest_list = NULL;
  AECRSamplePtr          sample;

  if (action == NULL || sep == NULL) return 0;

  psd.src = action->src;
  psd.portion = action->portion;
  psd.src_list = NULL;

  /* first, we need to get a list of the parse sources */  
  VisitBioseqsInSep (sep, &psd, FindParseSourceBioseqCallback);


  /* for each parse source, get a list of the destinations */
  for (vnp = psd.src_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;

    /* find destinations */
    AddParseDestinations (psip, action->dest);

    /* add destinations to list */
    ValNodeLink (&dest_list, psip->dest_list);
    psip->dest_list = NULL;
  }

  psd.src_list = ParseSourceListFree (psd.src_list);

  /* get sample for dest_list */
  sample = AECRSampleNew ();
  AddToSampleForDestList (sample, dest_list, action->dest);
  dest_list = ValNodeFree (dest_list);
  return sample;
}


static Int4 ApplyParseActionToSeqEntry (ParseActionPtr action, SeqEntryPtr sep)
{
  ParseSrcCollectionData psd;
  ParseSourceInfoPtr     psip;
  ValNodePtr             orgnames = NULL, source_list_for_removal = NULL, vnp;
  Int4                   num_succeeded = 0;

  if (action == NULL || sep == NULL) return 0;

  psd.src = action->src;
  psd.portion = action->portion;
  psd.src_list = NULL;

  /* first, we need to get a list of the parse sources */  
  VisitBioseqsInSep (sep, &psd, FindParseSourceBioseqCallback);

  if (action->capitalization != Cap_change_none) {
    /* if we will be fixing capitalization, get org names to use in fixes */
    VisitBioSourcesInSep (sep, &orgnames, GetOrgNamesInRecordCallback);
  }

  /* for each parse source, we need to get a list of the destinations */
  for (vnp = psd.src_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    if (action->remove_from_parsed) {
        ValNodeAddPointer (&source_list_for_removal, 0, ParseSourceInfoCopy (psip));
    }
    /* fix source text */
    FixCapitalizationInString (&(psip->parse_src_txt), action->capitalization, orgnames);

    /* find destinations */
    AddParseDestinations (psip, action->dest);

  }

  /* free orgname list if we created it */
  orgnames = ValNodeFree (orgnames);

  CombineSourcesForDestinations (&(psd.src_list));

  if (action->dest->choice == ParseDest_org) {
    PropagateSourceOnSegSetForParse (psd.src_list);
  }
  
  /* now do the parsing */
  for (vnp = psd.src_list; vnp != NULL; vnp = vnp->next) {
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    num_succeeded += SetFieldForDestList (psip->dest_list, action->dest, psip->parse_src_txt, action->existing_text);
  }

  /* now remove strings from sources */
  for (vnp = source_list_for_removal; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    StripFieldForSrcList (psip, action->src, action->portion);
  }

  psd.src_list = ParseSourceListFree (psd.src_list);
  return num_succeeded;
}


static void SetCdRegionGeneticCode (SeqFeatPtr cds)
{
  CdRegionPtr crp;
  SeqEntryPtr parent_sep;
  BioseqPtr   bsp;
  Int4        genCode;
  ValNodePtr  code, vnp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return;
  if (cds->data.value.ptrvalue == NULL) {
    cds->data.value.ptrvalue = CdRegionNew();
  }
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  bsp = BioseqFindFromSeqLoc (cds->location);
  if (bsp == NULL) return;
  parent_sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
  genCode = SeqEntryToGeneticCode (parent_sep, NULL, NULL, 0);

  code = ValNodeNew (NULL);
  if (code != NULL) {
    code->choice = 254;
    vnp = ValNodeNew (NULL);
    code->data.ptrvalue = vnp;
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = genCode;
    }
  }
  crp->genetic_code = code;
}

  
static void CreateDataForFeature (SeqFeatPtr sfp, Int4 feature_type)
{
  Int4 featdef, seqfeattype;
  CharPtr    label = NULL;
  RnaRefPtr  rrp;
  GBQualPtr  gbq;
  ImpFeatPtr ifp;

  featdef = GetFeatdefFromFeatureType (feature_type);
  sfp->idx.subtype = featdef;
  seqfeattype = FindFeatFromFeatDefType (featdef);
  switch (seqfeattype) {
    case SEQFEAT_GENE:
      sfp->data.value.ptrvalue = GeneRefNew();
      break;
    case SEQFEAT_CDREGION:
      sfp->data.value.ptrvalue = CdRegionNew();
      SetCdRegionGeneticCode (sfp);
      break;
    case SEQFEAT_RNA:
      rrp = RnaRefNew();
      rrp->ext.choice = 0;
      sfp->data.value.ptrvalue = rrp;
      switch (featdef) {
        case FEATDEF_preRNA:
          rrp->type = RNA_TYPE_premsg;
          break;
        case FEATDEF_mRNA:
          rrp->type = RNA_TYPE_mRNA;
          break;
        case FEATDEF_tRNA:
          rrp->type = RNA_TYPE_tRNA;
          break;
        case FEATDEF_rRNA:
          rrp->type = RNA_TYPE_rRNA;
          break;
        case FEATDEF_snRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          gbq = GBQualNew ();
          gbq->qual = StringSave ("ncRNA_class");
          gbq->val = StringSave ("snRNA");
          break;
        case FEATDEF_scRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          gbq = GBQualNew ();
          gbq->qual = StringSave ("ncRNA_class");
          gbq->val = StringSave ("scRNA");
          break;
        case FEATDEF_tmRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("tmRNA");
          break;
        case FEATDEF_ncRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          break;
      }
      break;
    case SEQFEAT_IMP:
      ifp = ImpFeatNew();
      sfp->data.value.ptrvalue = ifp;
      label = GetFeatureNameFromFeatureType (feature_type);
      ifp->key = StringSave (label);
      break;
  }
}


static void ExtraCDSCreationActions (SeqFeatPtr cds, SeqEntryPtr parent_sep)
{
  ByteStorePtr       bs;
  CharPtr            prot, ptr;
  BioseqPtr          bsp;
  Char               ch;
  Int4               i;
  SeqEntryPtr        psep, nsep;
  MolInfoPtr         mip;
  ValNodePtr         vnp, descr;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  Boolean            partial5, partial3;

  if (cds == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  /* Create corresponding protein sequence data for the CDS */

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
  if (NULL == bs)
    return;

  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (NULL == prot)
    return;

  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }

  /* Create the product protein Bioseq */
  
  bsp = BioseqNew ();
  if (NULL == bsp)
    return;
  
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  bsp->id = MakeNewProteinSeqId (cds->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  
  /* Create a new SeqEntry for the Prot Bioseq */
  
  psep = SeqEntryNew ();
  if (NULL == psep)
    return;
  
  psep->choice = 1;
  psep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
  
  /* Add a descriptor to the protein Bioseq */
  
  mip = MolInfoNew ();
  if (NULL == mip)
    return;
  
  mip->biomol = 8;
  mip->tech = 8;
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
  if (NULL == vnp)
    return;
  
  vnp->data.ptrvalue = (Pointer) mip;
  
  /**/
  
  descr = ExtractBioSourceAndPubs (parent_sep);

  AddSeqEntryToSeqEntry (parent_sep, psep, TRUE);
  nsep = FindNucSeqEntry (parent_sep);
  ReplaceBioSourceAndPubs (parent_sep, descr);
  SetSeqFeatProduct (cds, bsp);
  
  prp = ProtRefNew ();
  
  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    }
  }
}


static SeqLocPtr LocationFromApplyFeatureAction (BioseqPtr bsp, ApplyFeatureActionPtr action)
{
  LocationIntervalPtr l;
  SeqLocPtr slp = NULL;
  Uint1 strand = Seq_strand_plus;
  Int4  from, to;

  if (bsp == NULL || action == NULL || action->location == NULL) return NULL;

  if (!action->plus_strand) {
    strand = Seq_strand_minus;
  }
  if (action->location->choice == LocationChoice_interval) {
    l = (LocationIntervalPtr) action->location->data.ptrvalue;
    if (l != NULL) {
      from = MIN (l->from, l->to) - 1;
      to = MAX (l->from, l->to) - 1;
      slp = SeqLocIntNew (from, to, strand, SeqIdFindWorst (bsp->id));
    }
  } else if (action->location->choice == LocationChoice_whole_sequence) {
    slp = SeqLocIntNew (0, bsp->length - 1, strand, SeqIdFindWorst (bsp->id));
  }
  SetSeqLocPartial (slp, action->partial5, action->partial3);
  return slp;
}


static Boolean OkToApplyToBioseq (ApplyFeatureActionPtr action, BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  Int4 featdef;
  Boolean rval = TRUE;

  if (action == NULL || bsp == NULL) return FALSE;

  if (!action->add_redundant) {
    featdef = GetFeatdefFromFeatureType (action->type);
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &context);
    if (sfp != NULL) {
      rval = FALSE;
    }
  }
  return rval;
} 

static void AddParts (ApplyFeatureActionPtr action, BioseqSetPtr parts, ValNodePtr PNTR bsp_list)
{
  SeqEntryPtr sep;
  Int4         seg_num;

  if (action == NULL || !action->apply_to_parts
      || parts == NULL || parts->_class != BioseqseqSet_class_parts
      || bsp_list == NULL) {
    return;
  }

  if (action->only_seg_num > -1) {
    seg_num = 0;
    sep = parts->seq_set;
    while (seg_num < action->only_seg_num && sep != NULL) {
      sep = sep->next;
      seg_num++;
    }
    if (sep != NULL && IS_Bioseq (sep) && OkToApplyToBioseq (action, sep->data.ptrvalue)) {
      ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, sep->data.ptrvalue);
    }
  } else {
    for (sep = parts->seq_set; sep != NULL; sep = sep->next) {
      if (IS_Bioseq (sep) && OkToApplyToBioseq (action, sep->data.ptrvalue)) {
        ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, sep->data.ptrvalue);
      }
    }
  }  
}


static void AddSequenceOrParts (ApplyFeatureActionPtr action, BioseqPtr bsp, ValNodePtr PNTR bsp_list)
{
  BioseqSetPtr bssp, parts;
  SeqEntryPtr sep;

  if (action == NULL || bsp == NULL || bsp_list == NULL) return;

  if (bsp->idx.parenttype == OBJ_BIOSEQSET && bsp->idx.parentptr != NULL) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp->_class == BioseqseqSet_class_segset) {
      if (action->apply_to_parts) {
        sep = bssp->seq_set;
        while (sep != NULL && !IS_Bioseq_set (sep)) {
          sep = sep->next;
        }
        if (sep != NULL) {
          AddParts (action, sep->data.ptrvalue, bsp_list);
        }
      } else {
        if (OkToApplyToBioseq (action, bsp)) {
          ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
        }
      }       
    } else if (bssp->_class == BioseqseqSet_class_parts) {
      if (action->apply_to_parts) {
        AddParts (action, bssp, bsp_list);
      } else {
        parts = bssp;
        if (parts->idx.parenttype == OBJ_BIOSEQSET && parts->idx.parentptr != NULL) {
          bssp = (BioseqSetPtr) parts->idx.parentptr;
          if (IS_Bioseq (bssp->seq_set) && OkToApplyToBioseq (action, bssp->seq_set->data.ptrvalue)) {
            ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp_list);
          }
        }
      }
    } else {
      if (OkToApplyToBioseq (action, bsp)) {
        ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
      }
    }
  } else {
    if (OkToApplyToBioseq (action, bsp)) {
      ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
    }
  }
}

static void AddSequenceOrPartsFromSeqEntry (ApplyFeatureActionPtr action, SeqEntryPtr sep, ValNodePtr PNTR bsp_list)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  seq_set;

  if (action == NULL || sep == NULL) return;

  while (sep != NULL) {
    if (IS_Bioseq (sep)) {
      AddSequenceOrParts (action, sep->data.ptrvalue, bsp_list);
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == BioseqseqSet_class_segset) {
        /* find master segment */
        seq_set = bssp->seq_set;
        while (seq_set != NULL && !IS_Bioseq (seq_set)) {
          seq_set = seq_set->next;
        }
        if (seq_set != NULL) {
          AddSequenceOrParts (action, seq_set->data.ptrvalue, bsp_list);
        }
      } else if (bssp->_class == BioseqseqSet_class_nuc_prot) {
        /* find nucleotide sequence */
        seq_set = bssp->seq_set;
        if (seq_set != NULL) {
          if (IS_Bioseq_set (seq_set)) {
            /* nucleotide is segmented set */
            bssp = (BioseqSetPtr) seq_set->data.ptrvalue;
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset
                && bssp->seq_set != NULL && IS_Bioseq (bssp->seq_set)) {
              AddSequenceOrParts (action, bssp->seq_set->data.ptrvalue, bsp_list);
            }
          } else if (IS_Bioseq (seq_set)) {
            AddSequenceOrParts (action, seq_set->data.ptrvalue, bsp_list);
          }
        }
      } else {
        /* add from set members */
        AddSequenceOrPartsFromSeqEntry (action, bssp->seq_set, bsp_list);
      }
    }
    sep = sep->next;
  }  
}  
  

static void AdjustProteinSequenceForReadingFrame (SeqFeatPtr cds)
{
  BioseqPtr protbsp, bsp;
  ByteStorePtr bs;
  SeqFeatPtr   prot_sfp;
  Boolean      partial5, partial3;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return;

  protbsp = BioseqFindFromSeqLoc (cds->product);

  if (protbsp == NULL) {
    bsp = BioseqFindFromSeqLoc (cds->location);
    if (bsp != NULL) {
      ExtraCDSCreationActions (cds, GetBestTopParentForData (bsp->idx.entityID, bsp));
    }
  } else {
    bs = ProteinFromCdRegionExWithTrailingCodonHandling (cds,
                                              TRUE,
                                              FALSE,
                                              TRUE);
    protbsp->seq_data = (SeqDataPtr) BSFree ((ByteStorePtr)(protbsp->seq_data));
    protbsp->seq_data = (SeqDataPtr) bs;
    protbsp->length = BSLen (bs);
    prot_sfp = GetProtFeature (protbsp);
    if (prot_sfp == NULL) {
      prot_sfp = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
      prot_sfp->data.value.ptrvalue = ProtRefNew ();
      CheckSeqLocForPartial (cds->location, &partial5, &partial3);
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    } else {
      if (SeqLocLen (prot_sfp->location) != protbsp->length) {
        prot_sfp->location = SeqLocFree (prot_sfp->location);
        prot_sfp->location = SeqLocIntNew (0, protbsp->length - 1, Seq_strand_plus, SeqIdFindWorst (protbsp->id));   
        CheckSeqLocForPartial (cds->location, &partial5, &partial3);
        SetSeqLocPartial (prot_sfp->location, partial5, partial3);
        prot_sfp->partial = (partial5 || partial3);
      }
    }
  }
}


static Int4 ApplyApplyFeatureActionToSeqEntry (ApplyFeatureActionPtr action, SeqEntryPtr sep)
{
  ValNodePtr bsp_list = NULL, vnp, field_vnp;
  Int4       featdef, seqfeattype;
  BioseqPtr  bsp;
  SeqFeatPtr sfp;
  SeqLocPtr  slp;
  FeatQualLegalValPtr q;
  FeatureField f;
  SeqIdPtr   sip;
  SeqFeatPtr gene;
  Int4       num_created = 0;

  if (sep == NULL || action == NULL) return 0;

  /* first, get list of Bioseqs to apply features to */
  /* relevant values : seq_list, add_redundant, apply_to_parts, only_seg_num */
  if (action->seq_list != NULL && action->seq_list->choice == SequenceListChoice_list) {
    for (vnp = action->seq_list->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      sip = CreateSeqIdFromText (vnp->data.ptrvalue, sep);
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        AddSequenceOrParts (action, bsp, &bsp_list);
      }
    }  
  } else {
    AddSequenceOrPartsFromSeqEntry (action, sep, &bsp_list);
  }

  /* now add feature to each bioseq in list */
  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = vnp->data.ptrvalue;
    if (bsp == NULL) continue;
    featdef = GetFeatdefFromFeatureType (action->type);
    seqfeattype = FindFeatFromFeatDefType (featdef);
    slp = LocationFromApplyFeatureAction (bsp, action);
    sfp = CreateNewFeatureOnBioseq (bsp, seqfeattype, slp);
    if (sfp == NULL) continue;
    CreateDataForFeature (sfp, action->type);
    /* any extra actions */
    switch (action->type) {
      case Feature_type_cds :
        ExtraCDSCreationActions (sfp, GetBestTopParentForData (bsp->idx.entityID, bsp));
        break;
      case Feature_type_source :
        if (action->src_fields != NULL) {
          sfp->data.value.ptrvalue = ImpFeatFree (sfp->data.value.ptrvalue);
          sfp->data.choice = SEQFEAT_BIOSRC;
          sfp->data.value.ptrvalue = BioSourceFromSourceQualVals (action->src_fields);
        }
        break;
    }
    gene = NULL;
    for (field_vnp = action->fields; field_vnp != NULL; field_vnp = field_vnp->next) {
      q = (FeatQualLegalValPtr) field_vnp->data.ptrvalue;
      if (q != NULL) {
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = q->qual;        
        if (sfp->data.choice != SEQFEAT_GENE
            && (q->qual == Feat_qual_legal_gene || q->qual == Feat_qual_legal_gene_description)) {
          if (gene == NULL) {
            gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, slp);
            CreateDataForFeature (gene, Feature_type_gene);
          }
          f.type = Feature_type_gene;
          SetQualOnFeature (gene, &f, NULL, q->val, ExistingTextOption_replace_old);
        } else {
          f.type = action->type;
          SetQualOnFeature (sfp, &f, NULL, q->val, ExistingTextOption_replace_old);
        }
      }
    }
    if (action->type == Feature_type_cds) {
      /* retranslate, to account for change in reading frame */
      AdjustProteinSequenceForReadingFrame (sfp);
      /* after the feature has been created, then adjust it for gaps */
      /* Note - this step may result in multiple coding regions being created. */
      AdjustCDSLocationsForUnknownGapsCallback (sfp, NULL);
    }
    num_created++;
  }  
  return num_created;
}


typedef struct convertandremovefeaturecollection {
  Uint1 featdef;
  ValNodePtr constraint_set;
  ValNodePtr feature_list;
} ConvertAndRemoveFeatureCollectionData, PNTR ConvertAndRemoveFeatureCollectionPtr;

static void ConvertAndRemoveFeatureCollectionCallback (SeqFeatPtr sfp, Pointer data)
{
  ConvertAndRemoveFeatureCollectionPtr p;  

  if (sfp == NULL || data == NULL) return;

  p = (ConvertAndRemoveFeatureCollectionPtr) data;
  if (sfp->idx.subtype == p->featdef && DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, p->constraint_set)) {
    ValNodeAddPointer (&(p->feature_list), OBJ_SEQFEAT, sfp);
  }
}


static Int4 ApplyRemoveFeatureActionToSeqEntry (RemoveFeatureActionPtr action, SeqEntryPtr sep)
{
  ConvertAndRemoveFeatureCollectionData d;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Int4       num_deleted = 0;

  if (action == NULL) return 0;

  d.featdef = GetFeatdefFromFeatureType (action->type);
  d.constraint_set = action->constraint;
  d.feature_list = NULL;

  VisitFeaturesInSep (sep, &d, ConvertAndRemoveFeatureCollectionCallback);
  for (vnp = d.feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL) {
      sfp->idx.deleteme = TRUE;
      num_deleted ++;
    }
  }
  DeleteMarkedObjects (ObjMgrGetEntityIDForChoice(sep), 0, NULL);
  return num_deleted;
}


/* functions for converting features */

static Boolean ApplyConvertFeatureSrcOptions (SeqFeatPtr sfp, ValNodePtr src_options, Boolean keep_original)
{
  ConvertFromCDSOptionsPtr options = NULL;
  Boolean rval = FALSE;

  if (sfp == NULL) return FALSE;
  if (src_options == NULL) return TRUE;

  if (src_options->choice == ConvertFeatureSrcOptions_cds) {
    options = (ConvertFromCDSOptionsPtr) src_options->data.ptrvalue;
    if (options != NULL) {
      ApplyCDSOptionsToFeature (sfp, options->remove_mRNA, options->remove_gene, options->remove_transcript_id, keep_original);
      rval = TRUE;
    }
  }
  return rval;
} 

typedef Boolean (*ConvertFeatureFunc) PROTO ((SeqFeatPtr, Int4, ConvertFeatureDstOptionsPtr));

static void ApplyRNADestinationOptions (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  CharPtr existing_class;
  FeatureField ff;

  /* apply destination options */
  if (featdef_to == FEATDEF_ncRNA 
      && dst_options != NULL 
      && dst_options->choice == ConvertFeatureDstOptions_ncrna_class
      && !StringHasNoText (dst_options->data.ptrvalue)) {
    ff.type = Feature_type_ncRNA;
    ff.field = ValNodeNew (NULL);
    ff.field->choice = FeatQualChoice_legal_qual;
    ff.field->data.intvalue = Feat_qual_legal_ncRNA_class;
    existing_class = GetQualFromFeature (sfp, &ff, NULL);
    if (StringCmp (dst_options->data.ptrvalue, existing_class) != 0) {
      sfp->idx.subtype = FEATDEF_ncRNA;
      SetQualOnFeature (sfp, &ff, NULL, dst_options->data.ptrvalue, ExistingTextOption_append_semi);
    }
    existing_class = MemFree (existing_class);
    ff.field = ValNodeFree (ff.field);
  }      
}

static Boolean ConvertCDSToRNAFunc (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  Boolean rval;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  rval = ConvertCDSToRNA (sfp, featdef_to);
  if (rval) {
    ApplyRNADestinationOptions (sfp, featdef_to, dst_options);
  }
  return rval;
}


static Boolean ConvertGeneToRNAFunc (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  Boolean rval;

  rval = ConvertGeneToRNA (sfp, featdef_to);
  if (rval) {
    ApplyRNADestinationOptions (sfp, featdef_to, dst_options);
  }
  return rval;
}


static Boolean ConvertBioSrcToRegionFunc (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertBioSrcToRepeatRegion (sfp, featdef_to);
}


static Boolean ConvertCDSToMiscFeatFunc (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  Boolean rval = FALSE;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE; 
  }
  else if (sfp->pseudo) 
  {
    rval = ConvertOnePseudoCDSToMiscFeatEx (sfp, FALSE);
  }
  else
  {
    /* do other here */
    rval = ConvertNonPseudoCDSToMiscFeat (sfp, FALSE);
  }
  return rval;
}

static Boolean ConvertImpToProtFuncEx (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertImpToProtFunc (sfp, featdef_to);
}


static Boolean ConvertProtToImpFuncEx (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertProtToImpFunc (sfp, featdef_to);
}


static Boolean ConvertProtToProt (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertProtToProtFunc (sfp, featdef_to);
}


static Boolean ConvertImpToRNAFunc (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  RnaRefPtr          rrp;
  GBQualPtr          qual, qual_prev = NULL;
  Boolean            add_to_comment = FALSE;
  CharPtr            old_comment = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP)
  {
    return FALSE;
  }

  for (qual = sfp->qual; qual != NULL && StringCmp (qual->qual, "product") != 0; qual = qual->next) {
    qual_prev = qual;
  }
  if (qual != NULL) {
    old_comment = StringSave (qual->val);
    if (qual_prev == NULL) {
      sfp->qual = qual->next;
    } else {
      qual_prev->next = qual->next;
    }
    qual->next = NULL;
    qual = GBQualFree (qual);
  } else {
    old_comment = sfp->comment;
    sfp->comment = NULL;
  }

  rrp = RnaRefFromLabel (featdef_to, old_comment, &add_to_comment);
  
  sfp->data.value.ptrvalue = ImpFeatFree ((ImpFeatPtr) sfp->data.value.ptrvalue);
  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = (Pointer) rrp;
  SetRNAProductString (sfp, NULL, old_comment, ExistingTextOption_replace_old);
  
  if (add_to_comment) {
    SetStringValue (&(sfp->comment), old_comment, ExistingTextOption_append_semi);
  }
  old_comment = MemFree (old_comment);

  ApplyRNADestinationOptions (sfp, featdef_to, dst_options);

  return TRUE;
}


static Boolean ConvertRegionToImp (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertRegionToImpFunc (sfp, featdef_to);
}


static Boolean ConvertImpToImp (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertImpToImpFunc (sfp, featdef_to);
}


static Boolean ConvertRegionToRNA (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  Boolean rval;
  rval = ConvertRegionToRNAFunc (sfp, featdef_to);
  if (rval) {
    ApplyRNADestinationOptions (sfp, featdef_to, dst_options);
  }
  return rval;
}


static Boolean ConvertCommentToMiscFeat (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  ImpFeatPtr ifp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_COMMENT || sfp->data.value.ptrvalue != NULL)
  {
    return FALSE;
  }
  
  ifp = ImpFeatNew ();
  if (ifp != NULL) {
    ifp->key = StringSave ("misc_feature");
    sfp->data.choice = SEQFEAT_IMP;
    sfp->data.value.ptrvalue = (Pointer) ifp;
    return TRUE;
  }
  return FALSE;
}


static Boolean ConvertGeneToMiscFeat (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertGeneToMiscFeatFunc (sfp, featdef_to);
}


static Boolean ConvertRNAToImpFeat (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  CharPtr product = NULL;
  ImpFeatPtr ifp;
  Uint1      seqfeattype;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) {
    return FALSE;
  }

  seqfeattype = FindFeatFromFeatDefType (featdef_to);
  if (seqfeattype != SEQFEAT_IMP) {
    return FALSE;
  }

  product = GetRNAProductString (sfp, NULL);

  RemoveRNAProductString (sfp, NULL);

  sfp->data.value.ptrvalue = RnaRefFree (sfp->data.value.ptrvalue);

  ifp = ImpFeatNew ();
  ifp->key = StringSave (GetImportFeatureName (featdef_to));
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  SetStringValue (&(sfp->comment), product, ExistingTextOption_append_semi);
  product = MemFree (product);
  return TRUE;
}


static Boolean ConvertSiteToImpFeat (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  GBQualPtr  gbqual;
  ImpFeatPtr ifp;
  Int2       sitetype;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_SITE)
  {
    return FALSE;
  }

  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  sitetype = (Int2) sfp->data.value.intvalue;
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave (GetImportFeatureName (featdef_to));
  str = GetMacroSiteTypeName (MacroSiteTypeFromAsn1SiteType (sitetype));
  if (str != NULL) {
    gbqual = GBQualNew ();
    if (gbqual != NULL) {
      gbqual->qual = StringSave ("note");
      gbqual->val = StringSave (str);
      gbqual->next = sfp->qual;
      sfp->qual = gbqual;
    }
  }
  return TRUE;
}


static Boolean ConvertProtToRegion (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  ProtRefPtr prp;
  ValNodePtr vnp;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT)
  {
    return FALSE;
  }
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
  {
    return FALSE;
  }

  vnp = prp->name;
  if (vnp != NULL && vnp->next == NULL) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (! StringHasNoText (str)) {
      vnp->data.ptrvalue = NULL;
      sfp->data.value.ptrvalue = ProtRefFree (prp);
      sfp->data.choice = SEQFEAT_REGION;
      sfp->data.value.ptrvalue = (Pointer) str;
    }
  }
  return TRUE;
}


static Boolean ConvertRegionToProt (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  return ConvertRegionToProtFunc (sfp, featdef_to);
}


static Boolean ConvertToBond (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  SeqLocPtr   slp;
  BioseqPtr   bsp;
  SeqEntryPtr sep;
  Boolean     no_cds = FALSE;
  SeqFeatPtr  new_sfp;
  SeqIdPtr    sip;
  SeqBondPtr  sbp;
  SeqPntPtr   spp;

  if (sfp == NULL || featdef_to != FEATDEF_BOND || dst_options == NULL || dst_options->choice != ConvertFeatureDstOptions_bond) {
    return FALSE;
  }

  SeqFeatDataFree (&(sfp->data));
  sfp->data.choice = SEQFEAT_BOND;
  sfp->data.value.intvalue = Asn1BondTypeFromMacroBondType (dst_options->data.intvalue);
  
  bsp = BioseqFindFromSeqLoc (sfp->location);

  if (!ISA_aa (bsp->mol))
  {
    slp = GetProteinLocationForNucleotideFeatureConversion (sfp->location, &no_cds);
    if (no_cds || slp == NULL) {
      return FALSE;
    }
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = slp;
  }

  if (sfp->location->choice != SEQLOC_BOND) {
    sip = SeqLocId (sfp->location);
    if (sip != NULL) {
      sbp = SeqBondNew ();
      if (sbp != NULL) {
        slp = ValNodeNew (NULL);
        if (slp != NULL) {
          slp->choice = SEQLOC_BOND;
          slp->data.ptrvalue = (Pointer) sbp;
          spp = SeqPntNew ();
          if (spp != NULL) {
            spp->strand = SeqLocStrand (sfp->location);
            spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
            spp->point = SeqLocStart (sfp->location);
            sbp->a = spp;
          }
          spp = SeqPntNew ();
          if (spp != NULL) {
            spp->strand = SeqLocStrand (sfp->location);
            spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
            spp->point = SeqLocStop (sfp->location);
            sbp->b = spp;
          }
          sfp->location = SeqLocFree (sfp->location);
          sfp->location = slp;
        }
      }
    }
  }

  sfp->idx.subtype = 0;

  bsp = GetBioseqGivenSeqLoc (slp, sfp->idx.entityID);
  if (bsp == NULL) {
    return FALSE;
  } 
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return FALSE;
  }

  new_sfp = (SeqFeatPtr) AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
  sfp->idx.deleteme = TRUE;
  CreateNewFeature (sep, NULL, SEQFEAT_BOND, new_sfp);

  return TRUE;
}


static Boolean ConvertToSite (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  SeqLocPtr   slp;
  BioseqPtr   bsp;
  SeqEntryPtr sep;
  Boolean     no_cds = FALSE;
  SeqFeatPtr  new_sfp;

  if (sfp == NULL || featdef_to != FEATDEF_SITE || dst_options == NULL || dst_options->choice != ConvertFeatureDstOptions_site) {
    return FALSE;
  }

  SeqFeatDataFree (&(sfp->data));
  sfp->data.choice = SEQFEAT_SITE;
  sfp->data.value.intvalue = Asn1SiteTypeFromMacroSiteType (dst_options->data.intvalue);
  
  bsp = BioseqFindFromSeqLoc (sfp->location);

  if (!ISA_aa (bsp->mol))
  {
    slp = GetProteinLocationForNucleotideFeatureConversion (sfp->location, &no_cds);
    if (no_cds || slp == NULL) {
      return FALSE;
    }
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = slp;
  }

  sfp->idx.subtype = 0;

  bsp = GetBioseqGivenSeqLoc (slp, sfp->idx.entityID);
  if (bsp == NULL) {
    return FALSE;
  } 
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return FALSE;
  }

  new_sfp = (SeqFeatPtr) AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
  sfp->idx.deleteme = TRUE;
  CreateNewFeature (sep, NULL, SEQFEAT_SITE, new_sfp);

  return TRUE;  
}


static Boolean ConvertToRegion (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  BioseqPtr     bsp;
  RegionTypePtr r;
  Boolean       create_prot_feats, no_cds = FALSE;
  SeqLocPtr     slp;
  SeqEntryPtr   sep;
  SeqFeatPtr    new_sfp;

  if (sfp == NULL || featdef_to != FEATDEF_REGION || dst_options == NULL || dst_options->choice != ConvertFeatureDstOptions_region || dst_options->data.ptrvalue == NULL) {
    return FALSE;
  }

  r = (RegionTypePtr) dst_options->data.ptrvalue;
  create_prot_feats = !r->create_nucleotide;
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;

  if (ISA_aa (bsp->mol))
  {
    if (create_prot_feats)
    {
      slp = (SeqLocPtr) AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
    }
    else
    {
      slp = FindNucleotideLocationForProteinFeatureConversion (sfp->location);
    }
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = slp;
  }
  else if (create_prot_feats) 
  {
    slp = GetProteinLocationForNucleotideFeatureConversion (sfp->location, &no_cds);
    if (no_cds) {
      return FALSE;
    } 
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = slp;
  }

  bsp = GetBioseqGivenSeqLoc (sfp->location, sfp->idx.entityID);
  if (bsp == NULL) {
    return FALSE;
  } 

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return FALSE;
  }

  SeqFeatDataFree (&(sfp->data));
  sfp->data.choice = SEQFEAT_REGION;
  sfp->data.value.ptrvalue = sfp->comment;
  sfp->comment = NULL;

  new_sfp = (SeqFeatPtr) AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
  sfp->idx.deleteme = TRUE;
  CreateNewFeature (sep, NULL, SEQFEAT_REGION, new_sfp);
  return TRUE;  
}


static Boolean ConvertRNAToRNA (SeqFeatPtr sfp, Int4 featdef_to, ConvertFeatureDstOptionsPtr dst_options)
{
  RnaRefPtr  rrp;
  Boolean    add_to_comment = FALSE;
  CharPtr    product;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (NULL == rrp) {
    return FALSE;
  }

  product = GetRNAProductString (sfp, NULL);
  
  RemoveRNAProductString (sfp, NULL);

  sfp->data.value.ptrvalue = RnaRefFree (sfp->data.value.ptrvalue);

  sfp->data.value.ptrvalue = RnaRefFromLabel (featdef_to, product, &add_to_comment);

  SetRNAProductString (sfp, NULL, product, ExistingTextOption_replace_old);
  if (add_to_comment) {
    SetStringValue (&(sfp->comment), product, ExistingTextOption_append_semi);
  }
  product = MemFree (product);

  /* apply destination options */
  ApplyRNADestinationOptions (sfp, featdef_to, dst_options);

  sfp->idx.subtype = 0;
  return TRUE;
}


typedef struct convertfeattable {
  Uint2 seqfeat_from;
  Uint2 featdef_from;
  Uint2 seqfeat_to;
  Uint2 featdef_to;
  ConvertFeatureFunc func;
  CharPtr help_text;
} ConvertFeatTableData, PNTR ConvertFeatTablePtr;

static ConvertFeatTableData conversion_functions[] = {
  { SEQFEAT_CDREGION, FEATDEF_CDS,                SEQFEAT_RNA,    FEATDEF_ANY,
    ConvertCDSToRNAFunc,
    "Delete protein product sequence.\nClear product field if transcript ID removal was requested.\nIf converting to tRNA and anticodon value can be parsed from label, set aa value, and add any text that could not be parsed into an anticodon value to the feature note.\nIf converting to other RNA, put label in RNA product." },
  { SEQFEAT_GENE,     FEATDEF_GENE,               SEQFEAT_RNA,    FEATDEF_ANY,
    ConvertGeneToRNAFunc,
    "If converting to tRNA and anticodon value can be parsed from label, set aa value, and add any text that could not be parsed into an anticodon value to the feature note.  If converting to other RNA, put label in RNA product.  Also append gene locus, allele, description, map location, and locus tag to comment (as long as these values are not already in the label and therefore in the RNA product)." },
  { SEQFEAT_BIOSRC,   FEATDEF_BIOSRC,             SEQFEAT_IMP,    FEATDEF_repeat_region, 
    ConvertBioSrcToRegionFunc,
    "Creates a repeat_region with mobile_element qualifiers for the transposon and/or insertion sequence qualifiers on the BioSource.  All other BioSource information is discarded." },
  { SEQFEAT_CDREGION, FEATDEF_CDS,                SEQFEAT_IMP,    FEATDEF_misc_feature, 
    ConvertCDSToMiscFeatFunc,
    "Copy comment from coding region to new misc_feature and remove product field.  If not pseudo coding region, add product name from protein feature to new misc_feature comment and delete product sequence." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_PROT,   FEATDEF_ANY, 
    ConvertImpToProtFuncEx,
    "Original feature must be on nucleotide sequence and be contained in coding region location.  Coding region must have product protein sequence.  New feature is created on product protein sequence so that the translated location will be as close as possible to the original nucleotide location (may not be exact because of codon boundaries)." },
  { SEQFEAT_PROT,     FEATDEF_mat_peptide_aa,     SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertProtToImpFuncEx,
    "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
    "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
    "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
    "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_PROT,     FEATDEF_sig_peptide_aa,     SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertProtToImpFuncEx,
    "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
    "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
    "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
    "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_PROT,     FEATDEF_transit_peptide_aa, SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertProtToImpFuncEx, 
    "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
    "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
    "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
    "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_ANY,
    ConvertImpToRNAFunc,
    "Creates an RNA feature of the specified subtype.  Import feature key is discarded." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertRegionToImp,
    "Creates a misc_feature with the region name saved as a /note qualifier." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_RNA,    FEATDEF_ANY,
    ConvertRegionToRNA,
    "Creates an RNA feature with the region name as the product name." },
  { SEQFEAT_COMMENT,  FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_misc_feature,
    ConvertCommentToMiscFeat,
    "Creates a misc_feature with the same note as the original.  Note - the flatfile display for the feature is the same." },
  { SEQFEAT_GENE,     FEATDEF_GENE,               SEQFEAT_IMP,    FEATDEF_misc_feature,
    ConvertGeneToMiscFeat, 
    "Creates a misc_feature with the gene description and locus prepended to the original comment, separated by semicolons." },
  { SEQFEAT_RNA,      FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertRNAToImpFeat,
    "Creates an import feature of the specified subtype and adds the RNA product name to the comment." } ,
  { SEQFEAT_SITE,     FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertSiteToImpFeat,
    "Creates an import feature of the specified subtype with the site type name as a /note qualifier." } ,
  { SEQFEAT_PROT,     FEATDEF_mat_peptide_aa,     SEQFEAT_REGION, FEATDEF_REGION,
    NULL,
    "Creates a Region feature with the protein name as the region name." },
  { SEQFEAT_PROT,     FEATDEF_ANY,     SEQFEAT_REGION, FEATDEF_REGION,
    ConvertProtToRegion,
    "Creates a Region feature with the protein name as the region name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_PROT,   FEATDEF_ANY,
    ConvertRegionToProt, 
    "If feature is on nucleotide sequence, will create feature on protein product sequence for overlapping coding region.  Protein name will be region name." },
  { 0,                FEATDEF_ANY,                SEQFEAT_BOND,    FEATDEF_BOND,
    ConvertToBond,
    "Create Bond feature with specified bond type.  Location is a SeqLocBond with a point at the start of the original location and a point at the end of the original location.  All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { 0,                FEATDEF_ANY,                SEQFEAT_SITE,    FEATDEF_SITE,
    ConvertToSite,
    "Create Site feature with specified site type.  All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { 0,                FEATDEF_ANY,                SEQFEAT_REGION,    FEATDEF_REGION,
    ConvertToRegion,
    "Create Region feature on nucleotide sequence or protein product sequence of overlapping coding region as specified.  Use comment on feature for region name.\n"
    "All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_ANY,
    ConvertImpToImp,
    "Changes type of import feature." },
  { SEQFEAT_RNA,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_ANY,
    ConvertRNAToRNA, 
    "Changes type of RNA feature." },
  { SEQFEAT_PROT,     FEATDEF_ANY,                SEQFEAT_PROT,   FEATDEF_ANY,
    ConvertProtToProt, 
    "Changes type of protein feature." },
};


static Int4 num_convert_feature_table_lines = sizeof (conversion_functions) / sizeof (ConvertFeatTableData);

static Int4 GetConversionFunctionTableLine (Uint2 seqfeat_from, Uint2 featdef_from, Uint2 seqfeat_to, Uint2 featdef_to)
{
  Int4 i, table_line_num = -1;

  for (i = 0; i < num_convert_feature_table_lines && table_line_num == -1; i++)
  {
    if ((conversion_functions[i].seqfeat_from == 0 || conversion_functions[i].seqfeat_from == seqfeat_from)
        && (conversion_functions[i].featdef_from == FEATDEF_ANY || conversion_functions[i].featdef_from == featdef_from)
        && (conversion_functions[i].seqfeat_to == 0 || conversion_functions[i].seqfeat_to == seqfeat_to)
        && (conversion_functions[i].featdef_to == FEATDEF_ANY || conversion_functions[i].featdef_to == featdef_to))
    {
      table_line_num = i;
    }
  }
  return table_line_num;
}


NLM_EXTERN Boolean IsConversionSupported (Uint2 type_from, Uint2 type_to)
{
  Int4 line;
  Uint2 featdef_from, featdef_to, seqfeat_from, seqfeat_to;

  featdef_from = GetFeatdefFromFeatureType (type_from);
  seqfeat_from = FindFeatFromFeatDefType (featdef_from);
  featdef_to = GetFeatdefFromFeatureType (type_to);
  seqfeat_to = FindFeatFromFeatDefType (featdef_to);
  line = GetConversionFunctionTableLine (seqfeat_from, featdef_from, seqfeat_to, featdef_to);
  if (line > -1 && conversion_functions[line].func != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}



static Int4 ApplyConvertFeatureActionToSeqEntry (ConvertFeatureActionPtr action, SeqEntryPtr sep)
{
  ConvertAndRemoveFeatureCollectionData d;
  ValNodePtr vnp;
  SeqFeatPtr sfp, sfp_copy;
  Int4       num_affected = 0, table_line;
  Uint2      seqfeat_from, featdef_from, seqfeat_to, featdef_to;

  if (action == NULL) return 0;

  featdef_from = GetFeatdefFromFeatureType (action->type_from);
  seqfeat_from = FindFeatFromFeatDefType(featdef_from);
  featdef_to = GetFeatdefFromFeatureType (action->type_to);
  seqfeat_to = FindFeatFromFeatDefType (featdef_to);
  table_line = GetConversionFunctionTableLine (seqfeat_from, featdef_from, seqfeat_to, featdef_to);
  if (table_line < 0 || conversion_functions[table_line].func == NULL) {
    return 0;
  }

  d.featdef = GetFeatdefFromFeatureType (action->type_from);
  d.constraint_set = action->src_feat_constraint;
  d.feature_list = NULL;

  VisitFeaturesInSep (sep, &d, ConvertAndRemoveFeatureCollectionCallback);
  for (vnp = d.feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL) {
      sfp_copy = (SeqFeatPtr) AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
      /* add subtype value to copy */
      sfp_copy->idx.subtype = sfp->idx.subtype;
      sfp_copy->next = sfp->next;
      sfp->next = sfp_copy;
      if (conversion_functions[table_line].func (sfp_copy, featdef_to, action->dst_options)) {
        ApplyConvertFeatureSrcOptions (sfp_copy, action->src_options, action->leave_original);
        num_affected ++;
        if (!action->leave_original) {
          sfp->idx.deleteme = TRUE;
        }
      } else {
        sfp_copy->idx.deleteme = TRUE;
      }
    }
  }
  DeleteMarkedObjects (ObjMgrGetEntityIDForChoice(sep), 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);
  return num_affected;  
}


/* Functions for editing feature locations */
static Boolean DoesStrandMatch (Int4 strand_choice, Uint1 strand_val)
{
  Boolean rval = FALSE;
  
  switch (strand_choice)
  {
    case Feature_location_strand_from_any:
      rval = TRUE;
      break;
    case Feature_location_strand_from_unknown:
      if (strand_val == Seq_strand_unknown)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_plus:
      if (strand_val != Seq_strand_minus)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_minus:
      if (strand_val == Seq_strand_minus)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_both:
      if (strand_val == Seq_strand_both)
      {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Uint1 GetNewStrandValue (Int4 strand_choice, Uint1 strand_val)
{
  Uint1 rval = Seq_strand_unknown;
  
  switch (strand_choice)
  {
    case Feature_location_strand_to_reverse:
      switch (strand_val)
      {
        case Seq_strand_plus:
        case Seq_strand_unknown:
          rval = Seq_strand_minus;
          break;
        case Seq_strand_minus:
          rval = Seq_strand_plus;
          break;
        default:
          rval = strand_val;
          break;
      }
      break;
    case Feature_location_strand_to_unknown:
      rval = Seq_strand_unknown;
      break;
    case Feature_location_strand_to_plus:
      rval = Seq_strand_plus;
      break;
    case Feature_location_strand_to_minus:
      rval = Seq_strand_minus;
      break;
    case Feature_location_strand_to_both:
      rval = Seq_strand_both;
      break;
  }  
  return rval;
}


static Boolean ConvertLocationStrand (SeqLocPtr slp, Int4 fromStrand, Int4 toStrand)
{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;
  Boolean        rval = FALSE;
  Uint1          strand_orig;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL && DoesStrandMatch (fromStrand, sinp->strand)) 
        {
          strand_orig = sinp->strand;
          sinp->strand = GetNewStrandValue (toStrand, sinp->strand);
          if (strand_orig != sinp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand))
        {
          strand_orig = spp->strand;
          spp->strand = GetNewStrandValue (toStrand, spp->strand);
          if (strand_orig != spp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL && DoesStrandMatch (fromStrand, psp->strand)) 
        {
          strand_orig = psp->strand;
          psp->strand = GetNewStrandValue (toStrand, psp->strand);
          if (strand_orig != psp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          rval |= ConvertLocationStrand (loc, fromStrand, toStrand);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            strand_orig = spp->strand;
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
            if (strand_orig != spp->strand) {
              rval = TRUE;
            }
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            strand_orig = spp->strand;
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
            if (strand_orig != spp->strand) {
              rval = TRUE;
            }
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
  return rval;
}


static Boolean ApplyEditLocationStrandToSeqFeat (EditLocationStrandPtr edit, SeqFeatPtr sfp)
{
  Boolean rval = FALSE;

  if (edit == NULL || sfp == NULL) {
    return FALSE;
  }

  rval = ConvertLocationStrand (sfp->location, edit->strand_from, edit->strand_to);  
  return rval;
}


static Boolean At5EndOfSequence (SeqLocPtr slp, BioseqPtr bsp)
{
  Uint1 strand;
  Int4  start;
  Boolean at_end = FALSE;

  if (slp == NULL || bsp == NULL) return FALSE;

  strand = SeqLocStrand (slp);

  if (strand == Seq_strand_minus) {
    start = SeqLocStop (slp);
    if (start == bsp->length - 1) {
      at_end = TRUE;
    }
  } else {
    start = SeqLocStart (slp);
    if (start == 0) {
      at_end = TRUE;
    }
  }
  return at_end;
}


static Boolean HasGoodStartCodon (SeqFeatPtr sfp)
{
  ByteStorePtr bs;
  CharPtr      prot;
  Boolean     has_start = FALSE;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) {
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot != NULL && *prot == 'M') {
        has_start = TRUE;
      }
      prot = MemFree (prot);
    }
  }
  return has_start;
}


static Boolean ApplyPartial5SetActionToSeqFeat (Partial5SetActionPtr action, SeqFeatPtr sfp)
{
  Boolean      rval = FALSE;
  Boolean      make_partial = FALSE;
  Uint1        strand;
  BioseqPtr    bsp;
  CdRegionPtr  crp;
  Boolean      partial5, partial3;

  if (action == NULL || sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);

  switch (action->constraint) {
    case Partial_5_set_constraint_all:
      make_partial = TRUE;
      break;
    case Partial_5_set_constraint_at_end:
      make_partial = At5EndOfSequence (sfp->location, bsp);
      break;
    case Partial_5_set_constraint_bad_start:
      make_partial = HasGoodStartCodon (sfp);
      break;
    case Partial_5_set_constraint_frame_not_one:
      if (sfp->data.choice == SEQFEAT_CDREGION
          && (crp = sfp->data.value.ptrvalue) != NULL
          && crp->frame != 0 && crp->frame != 1) {
        make_partial = TRUE;
      }
      break;
  }

  if (make_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (!partial5) {
      SetSeqLocPartial (sfp->location, TRUE, partial3);
      if (action->extend && bsp != NULL) {
        ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
      }
      rval = TRUE; 
    }
  }
  return rval;
}


static Boolean ApplyClear5PartialToSeqFeat (Int4 action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE, clear_partial = FALSE;
  Boolean partial5, partial3;

  if (sfp == NULL) return FALSE;

  switch (action) {
    case Partial_5_clear_constraint_all:
      clear_partial = TRUE;
      break;
    case Partial_5_clear_constraint_not_at_end:
      clear_partial = !At5EndOfSequence(sfp->location, BioseqFindFromSeqLoc (sfp->location));
      break;
    case Partial_5_clear_constraint_good_start:
      clear_partial = !HasGoodStartCodon(sfp);
      break;
  }
  if (clear_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (partial5) {
      SetSeqLocPartial (sfp->location, FALSE, partial3);
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean At3EndOfSequence (SeqLocPtr slp, BioseqPtr bsp)
{
  Uint1 strand;
  Int4  stop;
  Boolean at_end = FALSE;

  if (slp == NULL || bsp == NULL) return FALSE;

  strand = SeqLocStrand (slp);

  if (strand == Seq_strand_minus) {
    stop = SeqLocStart (slp);
    if (stop == 0) {
      at_end = TRUE;
    }
  } else {
    stop = SeqLocStop (slp);
    if (stop == bsp->length - 1) {
      at_end = TRUE;
    }
  }
  return at_end;
}


static Boolean HasGoodStopCodon (SeqFeatPtr sfp)
{
  ByteStorePtr bs;
  CharPtr      prot;
  Boolean      has_stop = FALSE;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) {
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot != NULL && prot[StringLen (prot) - 1] == '*') {
        has_stop = TRUE;
      }
      prot = MemFree (prot);
    }
  }
  return has_stop;
}


static Boolean ApplyPartial3SetActionToSeqFeat (Partial3SetActionPtr action, SeqFeatPtr sfp)
{
  Boolean      rval = FALSE;
  Boolean      make_partial = FALSE;
  Uint1        strand;
  BioseqPtr    bsp;
  Boolean      partial5, partial3;

  if (action == NULL || sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);

  switch (action->constraint) {
    case Partial_3_set_constraint_all:
      make_partial = TRUE;
      break;
    case Partial_3_set_constraint_at_end:
      make_partial = At3EndOfSequence (sfp->location, bsp);
      break;
    case Partial_3_set_constraint_bad_end:
      make_partial = HasGoodStopCodon (sfp);
      break;
  }

  if (make_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (!partial3) {
      SetSeqLocPartial (sfp->location, partial5, TRUE);
      if (action->extend && bsp != NULL) {
        ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
      }
      rval = TRUE; 
    }
  }
  return rval;
}


static Boolean ApplyClear3PartialToSeqFeat (Int4 action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE, clear_partial = FALSE;
  Boolean partial5, partial3;

  if (sfp == NULL) return FALSE;

  switch (action) {
    case Partial_3_clear_constraint_all:
      clear_partial = TRUE;
      break;
    case Partial_3_clear_constraint_not_at_end:
      clear_partial = !At3EndOfSequence(sfp->location, BioseqFindFromSeqLoc (sfp->location));
      break;
    case Partial_3_clear_constraint_good_end:
      clear_partial = !HasGoodStopCodon(sfp);
      break;
  }
  if (clear_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (partial3) {
      SetSeqLocPartial (sfp->location, partial5, FALSE);
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean ApplyConvertLocationToSeqFeat (Int4 convert_location, SeqFeatPtr sfp)
{
  Boolean hasNulls, rval = FALSE;
  SeqLocPtr slp;
  BioseqPtr bsp;
  Boolean   partial5, partial3;

  if (sfp == NULL || (bsp = BioseqFindFromSeqLoc (sfp->location))== NULL) {
    return FALSE;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
	hasNulls = LocationHasNullsBetween (sfp->location);
	switch (convert_location) 
	{
	  case Convert_location_type_join :
	    if (hasNulls) 
	    {
	      slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, FALSE);
		    sfp->location = SeqLocFree (sfp->location);
		    sfp->location = slp;
		    if (bsp->repr == Seq_repr_seg) 
		    {
		      slp = SegLocToPartsEx (bsp, sfp->location, FALSE);
		      sfp->location = SeqLocFree (sfp->location);
		      sfp->location = slp;
		      hasNulls = LocationHasNullsBetween (sfp->location);
		      sfp->partial = (sfp->partial || hasNulls);
		    }
		    FreeAllFuzz (sfp->location);
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
	    }
	    break;
  	case Convert_location_type_order :
	    if (!hasNulls) 
	    {
		    slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, TRUE);
        sfp->location = SeqLocFree (sfp->location);
		    sfp->location = slp;
		    if (bsp->repr == Seq_repr_seg) 
		    {
		      slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
		      sfp->location = SeqLocFree (sfp->location);
		      sfp->location = slp;
		      hasNulls = LocationHasNullsBetween (sfp->location);
		      sfp->partial = (sfp->partial || hasNulls);
		    }
		    FreeAllFuzz (sfp->location);
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
	    }
	    break;
	  case Convert_location_type_merge :
      if (sfp->location->choice != SEQLOC_INT) {
	      slp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
	      sfp->location = SeqLocFree (sfp->location);
	      sfp->location = slp;
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
      }
	  default:
	    break;
	}
  return rval;
}


static Boolean ApplyLocationEditTypeToSeqFeat (ValNodePtr action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE;

  if (action == NULL || sfp == NULL) {
    return FALSE;
  }

  switch (action->choice) {
    case LocationEditType_strand:
      rval = ApplyEditLocationStrandToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_set_5_partial:
      rval = ApplyPartial5SetActionToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_clear_5_partial:
      rval = ApplyClear5PartialToSeqFeat (action->data.intvalue, sfp);
      break;
    case LocationEditType_set_3_partial:
      rval = ApplyPartial3SetActionToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_clear_3_partial:
      rval = ApplyClear3PartialToSeqFeat (action->data.intvalue, sfp);
      break;
    case LocationEditType_convert:
      rval = ApplyConvertLocationToSeqFeat (action->data.intvalue, sfp);
      break;
  }
  return rval;
}


static Int4 ApplyEditFeatureLocationActionToSeqEntry (EditFeatureLocationActionPtr action, SeqEntryPtr sep)
{
  ConvertAndRemoveFeatureCollectionData d;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Int4       num_affected = 0;

  if (action == NULL) return 0;

  d.featdef = GetFeatdefFromFeatureType (action->type);
  d.constraint_set = action->constraint;
  d.feature_list = NULL;

  VisitFeaturesInSep (sep, &d, ConvertAndRemoveFeatureCollectionCallback);
  for (vnp = d.feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL && ApplyLocationEditTypeToSeqFeat (action->action, sfp)) {
      num_affected++;
    }
  }
  return num_affected;
}


static void ApplyMolinfoBlockCallback (BioseqPtr bsp, Pointer data)
{
  MolinfoBlockPtr mib;
  ValNodePtr      field;
  MolInfoPtr      mip;

  if (bsp == NULL) {
    return;
  }

  mib = (MolinfoBlockPtr) data;
  if (mib == NULL) {
    return;
  }

  if (!DoesObjectMatchConstraintChoiceSet (OBJ_BIOSEQ, bsp, mib->constraint)) {
    return;
  }

  mip = GetMolInfoForBioseq (bsp);

  for (field = mib->from_list; field != NULL; field = field->next) {
    switch (field->choice) {
      case MolinfoField_molecule:
        if (mip == NULL || mip->biomol != BiomolFromMoleculeType (field->data.intvalue)) {
          return;
        }
        break;
      case MolinfoField_technique:
        if (mip == NULL || mip->tech != TechFromTechniqueType (field->data.intvalue)) {
          return;
        }
        break;
      case MolinfoField_completedness:
        if (mip == NULL || mip->completeness != CompletenessFromCompletednessType (field->data.intvalue)) {
          return;
        }
        break;
      case MolinfoField_mol_class:
        if (bsp->mol != MolFromMoleculeClassType (field->data.intvalue)) { 
          return;
        }
        break;
      case MolinfoField_topology:
        if (bsp->topology != TopologyFromTopologyType (field->data.intvalue)) {
          return;
        }
        break;
      case MolinfoField_strand:
        if (bsp->strand != StrandFromStrandType (field->data.intvalue)) {
          return;
        }
        break;
    }
  }


  for (field = mib->to_list; field != NULL; field = field->next) {
    SetSequenceQualOnBioseq (bsp, field);
  }
}


NLM_EXTERN void ApplyMolinfoBlockToSeqEntry (SeqEntryPtr sep, MolinfoBlockPtr mib)
{
  VisitBioseqsInSep (sep, mib, ApplyMolinfoBlockCallback);

}


typedef struct descriptortypename {
  Int4 descriptortype;
  Uint1 descriptor_choice;
  CharPtr descriptorname;
} DescriptorTypeNameData, PNTR DescriptorTypeNamePtr;

static DescriptorTypeNameData descriptortypename[] = {
 { Descriptor_type_all , 0 , "Any" } , 
 { Descriptor_type_title , Seq_descr_title , "Title" } , 
 { Descriptor_type_source , Seq_descr_source , "Source" } , 
 { Descriptor_type_publication , Seq_descr_pub , "Publication" } , 
 { Descriptor_type_comment , Seq_descr_comment , "Comment" } , 
 { Descriptor_type_genbank , Seq_descr_genbank , "GenBank" } , 
 { Descriptor_type_user , Seq_descr_user , "User" } , 
 { Descriptor_type_create_date , Seq_descr_create_date , "CreateDate" } , 
 { Descriptor_type_update_date , Seq_descr_update_date , "UpdateDate" } , 
 { Descriptor_type_mol_info , Seq_descr_molinfo , "MolInfo" } ,
 { Descriptor_type_structured_comment , Descriptor_type_user , "StructuredComment" }
};

#define NUM_descriptortypename sizeof (descriptortypename) / sizeof (DescriptorTypeNameData)

static Int4 GetDescriptorTypeFromDescriptorChoice (Uint1 descriptor_choice) 
{
  Int4 i;

  for (i = 0; i < NUM_descriptortypename; i++) {
    if (descriptor_choice == descriptortypename[i].descriptor_choice) {
      return descriptortypename[i].descriptortype;
    }
  }
  return -1;
}


static Uint1 GetDescriptorChoiceFromDescriptorType (Int4 descriptortype) 
{
  Int4 i;

  for (i = 0; i < NUM_descriptortypename; i++) {
    if (descriptortype == descriptortypename[i].descriptortype) {
      return descriptortypename[i].descriptor_choice;
    }
  }
  return SEQDESCR_MAX;
}


NLM_EXTERN CharPtr GetDescriptorNameFromDescriptorType (Int4 descriptortype)
{
  CharPtr str = NULL;
  Int4 i;

  for (i = 0; i < NUM_descriptortypename && str == NULL; i++) {
    if (descriptortype == descriptortypename[i].descriptortype) {
      str = descriptortypename[descriptortype].descriptorname;
    }
  } 
  if (str == NULL) {
    str = "Unknown descriptor type";
  }
  return str;
}


NLM_EXTERN void AddAllDescriptorsToChoiceList (ValNodePtr PNTR descriptor_type_list)
{
  Int4 i;
  ValNodePtr tmp_list = NULL;

  for (i = 0; i < NUM_descriptortypename; i++) {
    ValNodeAddPointer (&tmp_list, descriptortypename[i].descriptortype, StringSave (descriptortypename[i].descriptorname));
  }
  tmp_list = ValNodeSort (tmp_list, SortVnpByString);
  ValNodeLink (descriptor_type_list, tmp_list);
}



static Boolean DoesDescriptorMatchType (SeqDescrPtr sdp, Int4 descriptortype)
{
  Uint1 descriptorchoice;
  UserObjectPtr uop;

  if (sdp == NULL) {
    return FALSE;
  } else if (descriptortype == Descriptor_type_all) {
    return TRUE;
  } else if ((descriptorchoice = GetDescriptorChoiceFromDescriptorType (descriptortype)) == SEQDESCR_MAX) {
    return FALSE;
  } else if (descriptorchoice != sdp->choice) {
    return FALSE;
  } else if (descriptortype == Descriptor_type_structured_comment) {
    if ((uop = (UserObjectPtr) sdp->data.ptrvalue) == NULL
        || uop->type == NULL
        || StringCmp (uop->type->str, "StructuredComment") == 0) {
      return FALSE;
    } else {
      return TRUE;
    }
  } else {
    return TRUE;
  }
}


typedef struct removedescriptoractioncollection {
  RemoveDescriptorActionPtr action;
  ValNodePtr obj_list;
} RemoveDescriptorActionCollectionData, PNTR RemoveDescriptorActionCollectionPtr;


static void RemoveDescriptorCollectionCallback (SeqDescrPtr sdp, Pointer data)
{
  RemoveDescriptorActionCollectionPtr d;

  if (sdp == NULL || (d = (RemoveDescriptorActionCollectionPtr) data) == NULL
      || d->action == NULL) {
    return;
  }
  
  if (DoesDescriptorMatchType (sdp, d->action->type)
      && DoesObjectMatchConstraintChoiceSet (OBJ_SEQDESC, sdp, d->action->constraint)) {
    ValNodeAddPointer (&(d->obj_list), OBJ_SEQDESC, sdp);
  }
}


static Int4 ApplyRemoveDescriptorActionToSeqEntry (RemoveDescriptorActionPtr action, SeqEntryPtr sep)
{
  RemoveDescriptorActionCollectionData d;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  ValNodePtr vnp;
  Int4       num_deleted = 0;

  if (action == NULL) return 0;

  d.action = action;
  d.obj_list = NULL;

  VisitDescriptorsInSep (sep, &d, RemoveDescriptorCollectionCallback);
  for (vnp = d.obj_list; vnp != NULL; vnp = vnp->next) {
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
      num_deleted ++;
    }
  }
  DeleteMarkedObjects (ObjMgrGetEntityIDForChoice(sep), 0, NULL);
  return num_deleted;
}


static DefLineType DefLineTypeFromAutodefListType(Uint2 list_type)
{
  DefLineType deflinetype = DEFLINE_USE_FEATURES;

  switch (list_type) {
    case Autodef_list_type_feature_list:
      deflinetype = DEFLINE_USE_FEATURES;
      break;
    case Autodef_list_type_complete_sequence:
      deflinetype = DEFLINE_COMPLETE_SEQUENCE;
      break;
    case Autodef_list_type_complete_genome:
      deflinetype = DEFLINE_COMPLETE_GENOME;
      break;
  }
  return deflinetype;
}


static void ApplyAutodefActionToSeqEntry (AutodefActionPtr action, SeqEntryPtr sep)
{
  OrganismDescriptionModifiers od;
  ModifierItemLocalPtr modList;
  DeflineFeatureRequestList dfrl;
  ValNodePtr           vnp, modifier_indices = NULL;
  ValNode              field_type, source_qual_choice;
  Uint4                i;
  Int4                 defline_pos;

  od.allow_semicolon_in_modifier = FALSE;
  od.clone_isolate_HIV_rule_num = clone_isolate_HIV_rule_want_both;
  od.exclude_aff = FALSE;
  od.exclude_cf = FALSE;
  od.exclude_nr = FALSE;
  od.exclude_sp = TRUE;
  od.include_country_extra = FALSE;
  od.keep_paren = TRUE;
  od.max_mods = -99;
  od.use_labels = TRUE;
  od.use_modifiers = TRUE;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  for (i = 0; i < NumDefLineModifiers(); i++) {
    modList[i].any_present = FALSE;
    modList[i].all_present = FALSE;
    modList[i].is_unique = FALSE;
    modList[i].first_value_seen = NULL;
    modList[i].values_seen = NULL;
    modList[i].all_unique = FALSE;
    modList[i].status = NULL;
    modList[i].required = FALSE;
  }
  SetRequiredModifiers (modList);

  /* add modifiers specified in action */
  source_qual_choice.next = NULL;
  source_qual_choice.choice = SourceQualChoice_textqual;
  field_type.next = NULL;
  field_type.choice = FieldType_source_qual;
  field_type.data.ptrvalue = &source_qual_choice;

  for (vnp = action->modifiers; vnp != NULL; vnp = vnp->next) {
    source_qual_choice.data.intvalue = vnp->data.intvalue;
    defline_pos = GetDeflinePosForFieldType (&field_type);
    if (defline_pos > -1) {
      modList[defline_pos].required = TRUE;
      modList[defline_pos].any_present = TRUE;
      ValNodeAddInt (&modifier_indices, 0, defline_pos);

    }
  }

  InitFeatureRequests (&dfrl);
  dfrl.feature_list_type = DefLineTypeFromAutodefListType (action->clause_list_type);
  
  AutoDefForSeqEntry (sep, SeqMgrGetEntityIDForSeqEntry (sep), &od, modList, modifier_indices, &dfrl,
                      DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

  modList = MemFree (modList);
  modifier_indices = ValNodeFree (modifier_indices);

}


NLM_EXTERN void ApplyMacroToSeqEntry (SeqEntryPtr sep, ValNodePtr macro, Int4Ptr pNumFields, Int4Ptr pNumFeat)
{
  Int4 num_AECR = 0, num_parse = 0, num_feature = 0, num_fields = 0;

  while (macro != NULL) {
    switch (macro->choice) {
      case MacroActionChoice_aecr:
        num_AECR += ApplyAECRActionToSeqEntry ((AECRActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_parse:
        num_parse += ApplyParseActionToSeqEntry ((ParseActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_add_feature:
        num_feature += ApplyApplyFeatureActionToSeqEntry ((ApplyFeatureActionPtr) macro->data.ptrvalue, sep);
        SeqMgrIndexFeatures (ObjMgrGetEntityIDForChoice(sep), NULL);
        break;
      case MacroActionChoice_remove_feature:
        num_feature += ApplyRemoveFeatureActionToSeqEntry ((RemoveFeatureActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_edit_location:
        num_fields += ApplyEditFeatureLocationActionToSeqEntry ((EditFeatureLocationActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_convert_feature:
        num_feature += ApplyConvertFeatureActionToSeqEntry ((ConvertFeatureActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_remove_descriptor:
        num_feature += ApplyRemoveDescriptorActionToSeqEntry ((RemoveDescriptorActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_autodef:
        ApplyAutodefActionToSeqEntry ((AutodefActionPtr) macro->data.ptrvalue, sep);
        break;
    }
    macro = macro->next;
  }
  if (pNumFields != NULL) {
    *pNumFields = num_AECR + num_parse + num_fields;
  }
  if (pNumFeat != NULL) {
    *pNumFeat = num_feature;
  }
}


/* for generating text descriptions of macro objects */
NLM_EXTERN CharPtr SummarizeSourceQual (ValNodePtr field)
{
  CharPtr summ = NULL, locname, origname;
  Int4    genome, origin;
  CharPtr loc_fmt = "location %s";
  CharPtr orig_fmt = "origin %s";

  if (field == NULL) return NULL;
  switch (field->choice) {
    case SourceQualChoice_textqual:
      summ = StringSave (GetSourceQualName (field->data.intvalue));
      break;
    case SourceQualChoice_location:
      genome = GenomeFromSrcLoc (field->data.intvalue);
      locname = LocNameFromGenome (genome);
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (loc_fmt) + StringLen (locname)));
      sprintf (summ, loc_fmt, locname);
      break;
    case SourceQualChoice_origin:
      origin = OriginFromSrcOrig (field->data.intvalue);
      origname = OriginNameFromOrigin (origin);
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (orig_fmt) + StringLen (origname)));
      sprintf (summ, orig_fmt, origname);
      break;
  }
  return summ;
}


NLM_EXTERN CharPtr FeatureFieldLabel (CharPtr feature_name, ValNodePtr field)
{
  CharPtr cp;
  CharPtr label = NULL;
  CharPtr legal_fmt = "%s %s";
  CharPtr illegal_fmt = "constrained field on %s";
  
  if (feature_name == NULL) {
    feature_name = "Unknown feature";
  }

  if (field == NULL) {
    return StringSave ("missing field");
  } else if (field->choice == FeatQualChoice_legal_qual) {
    cp = GetFeatQualName (field->data.intvalue);
    if (cp == NULL) cp = "Unknown field type";
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (legal_fmt) + StringLen (feature_name) + StringLen (cp)));
    sprintf (label, legal_fmt, feature_name, cp);
  } else if (field->choice == FeatQualChoice_illegal_qual) {
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (illegal_fmt) + StringLen (feature_name)));
    sprintf (label, illegal_fmt, feature_name);
  } else {
    label = StringSave ("illegal field value");
  }
  return label;
}


NLM_EXTERN CharPtr SummarizeFieldType (ValNodePtr vnp)
{
  FeatureFieldPtr ffp;
  CharPtr str = NULL;
  CharPtr    label = NULL;
  CharPtr pub_fmt = "publication %s";

  if (vnp == NULL) {
    str = StringSave ("missing field");
  } else {
    switch (vnp->choice) {
      case FieldType_source_qual:
        str = SummarizeSourceQual (vnp->data.ptrvalue);
        break;
      case FieldType_feature_field:
        ffp = (FeatureFieldPtr) vnp->data.ptrvalue;
        if (ffp == NULL || ffp->field == NULL) {
          str = StringSave ("missing field");
        } else {
          label = GetFeatureNameFromFeatureType (ffp->type);
          str = FeatureFieldLabel (label, ffp->field);
        }
        break;
      case FieldType_cds_gene_prot:
        str = StringSaveNoNull (CDSGeneProtNameFromField (vnp->data.intvalue));
        if (str == NULL) {
          str = StringSave ("Invalid CDS-Gene-Prot Field");
        }
        break;
      case FieldType_molinfo_field:
        str = GetSequenceQualName (vnp->data.ptrvalue);
        if (str == NULL) {
          str = StringSave ("Invalid Sequence Qual Field");
        }
        break;
      case FieldType_pub:
        switch (vnp->data.intvalue) {
          case Publication_field_cit:
            str = StringSave ("publication citation");
            break;
          case Publication_field_authors:
            str = StringSave ("publication authors");
            break;
          case Publication_field_journal:
            str = StringSave ("publication journal");
            break;
          case Publication_field_volume:
            str = StringSave ("publication volume");
            break;
          case Publication_field_issue:
            str = StringSave ("publication issue");
            break;
          case Publication_field_pages:
            str = StringSave ("publication pages");
            break;
          case Publication_field_date:
            str = StringSave ("publication date");
            break;
          case Publication_field_serial_number:
            str = StringSave ("publication serial number");
            break;
          case Publication_field_title:
            str = StringSave ("publication title");
            break;
          default:
            label = GetPubFieldLabel (vnp->data.intvalue);
            if (label == NULL) {
              str = StringSave ("Invalid field type");
            } else {
              str = MemNew (sizeof (Char) * (StringLen (pub_fmt) + StringLen (label)));
              sprintf (str, pub_fmt, label);
            }
            break;
        }
        break;
      case FieldType_rna_field:
        str = SummarizeRnaQual (vnp->data.ptrvalue);
        break;
      case FieldType_struc_comment_field:
        str = SummarizeStructuredCommentField (vnp->data.ptrvalue);
        break;
      case FieldType_misc:
        if (vnp->data.intvalue == Misc_field_genome_project_id) {
          str = StringSave ("Genome Project ID");
        } else if (vnp->data.intvalue == Misc_field_comment_descriptor) {
          str = StringSave ("Comment Descriptor");
        } else {
          str = StringSave ("Invalid field type");
        }
        break;
      default:
        str = StringSave ("Invalid field type");
        break;
    }
  }
  return str;
}


NLM_EXTERN FieldTypePtr FieldTypeFromString (CharPtr str)
{
  Int4 qual_type, feat_type = -1;
  FieldTypePtr ft = NULL;
  FeatureFieldPtr ffp;
  ValNodePtr   vnp;
  CharPtr      cpy, cp;
  RnaQualPtr   rq;

  if (StringHasNoText (str)) {
    return NULL;
  }

  /* check source quals first */
  qual_type = GetSourceQualTypeByName (str);
  if (qual_type > -1) {
    vnp = ValNodeNew (NULL);
    vnp->choice = SourceQualChoice_textqual;
    vnp->data.intvalue = qual_type;
    ft = ValNodeNew (NULL);
    ft->choice = FieldType_source_qual;
    ft->data.ptrvalue = vnp;
  } else {
    /* try feature fields */
    cpy = StringSave (str);
    cp = StringChr (cpy, ' ');
    while (cp != NULL && feat_type == -1) {
      *cp = 0;
      feat_type = GetFeatureTypeByName (cpy);
      if (feat_type < 0) {
        *cp = ' ';
        cp = StringChr (cp + 1, ' ');
      }
    }
    if (feat_type > -1) {
      qual_type = GetFeatQualByName (cp + 1);
      if (qual_type > -1) {
        ffp = FeatureFieldNew ();
        ffp->type = feat_type;
        ValNodeAddInt (&ffp->field, FeatQualChoice_legal_qual, qual_type);
        ft = ValNodeNew (NULL);
        ft->choice = FieldType_feature_field;
        ft->data.ptrvalue = ffp;
      }
    }
    cpy = MemFree (cpy);
    if (ft == NULL) {
      /* try CDS-gene-prot */
      qual_type = CDSGeneProtFieldFromName (str);
      if (qual_type > -1) {
        ft = ValNodeNew (NULL);
        ft->choice = FieldType_cds_gene_prot;
        ft->data.intvalue = qual_type;
      }
    }
    if (ft == NULL) {
      /* try RNA Quals */
      cpy = StringSave (str);
      cp = StringChr (cpy, ' ');
      if (cp != NULL) {
        *cp = 0;
        feat_type = GetRnaTypeForName (cpy);
        qual_type = GetRnaFieldForName (cp + 1);
        if (feat_type > -1 && qual_type > -1) {
          rq = RnaQualNew ();
          rq->type = ValNodeNew (NULL);
          rq->type->choice = feat_type;
          rq->type->data.ptrvalue = NULL;
          rq->field = qual_type;
          ft = ValNodeNew (NULL);
          ft->choice = FieldType_rna_field;
          ft->data.ptrvalue = rq;
        }
      }
      cpy = MemFree (cpy);
    }
  }
  return ft;
}


/* for table readers that use the macro language functions */

/* MatchType is used to represent how the column should be matched.
 */

NLM_EXTERN MatchTypePtr MatchTypeNew ()
{
  MatchTypePtr match_type = MemNew (sizeof (MatchTypeData));
  match_type->data = NULL;
  match_type->match_location = String_location_equals;
  match_type->choice = eTableMatchNucID;
  return match_type;
}


NLM_EXTERN MatchTypePtr MatchTypeFree (MatchTypePtr match_type)
{
  if (match_type != NULL) {
    if (match_type->choice == eTableMatchSourceQual) {
      match_type->data = SourceQualChoiceFree (match_type->data);
    }
    match_type = MemFree (match_type);
  }
  return match_type;
}


static MatchTypePtr MatchTypeCopy (MatchTypePtr orig)
{
  MatchTypePtr match_type = NULL;
  
  if (orig != NULL) {
    match_type = MatchTypeNew();
    match_type->choice = orig->choice;
    match_type->match_location = orig->match_location;
    if (match_type->choice == eTableMatchSourceQual) {
      match_type->data = AsnIoMemCopy (orig->data, (AsnReadFunc) SourceQualChoiceAsnRead, (AsnWriteFunc) SourceQualChoiceAsnWrite);
    }
  }
  return match_type;
}


static MatchTypePtr FindMatchTypeInHeader (ValNodePtr columns)
{
  ValNodePtr col_vnp;
  MatchTypePtr match_type = NULL;
  TabColumnConfigPtr t;

  for (col_vnp = columns;
        col_vnp != NULL && match_type == NULL;
        col_vnp = col_vnp->next) {
    t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
    if (t != NULL && t->match_type != NULL) {
      match_type = MatchTypeCopy (t->match_type);
    }
  }
  return match_type;
}


NLM_EXTERN TabColumnConfigPtr TabColumnConfigNew (void)
{
  TabColumnConfigPtr t;

  t = (TabColumnConfigPtr) MemNew (sizeof (TabColumnConfigData));
  t->match_type = NULL;
  t->field = NULL;
  t->existing_text = ExistingTextOption_replace_old;
  t->constraint = NULL;
  t->skip_blank = TRUE;
  return t;
}


NLM_EXTERN TabColumnConfigPtr TabColumnConfigFree (TabColumnConfigPtr t)
{
  if (t != NULL) {
    t->field = FieldTypeFree (t->field);
    t->match_type = MatchTypeFree (t->match_type);
    t->constraint = ConstraintChoiceSetFree (t->constraint);
    t = MemFree (t);
  }
  return t;
}


NLM_EXTERN TabColumnConfigPtr TabColumnConfigCopy (TabColumnConfigPtr orig)
{
  TabColumnConfigPtr t = NULL;

  if (orig != NULL) {
    t = TabColumnConfigNew ();

    t->match_type = MatchTypeCopy (orig->match_type);
    t->existing_text = orig->existing_text;
    t->skip_blank = orig->skip_blank;
    t->match_mrna = orig->match_mrna;
    t->field = FieldTypeCopy (orig->field);
    t->constraint = AsnIoMemCopy (orig->constraint, (AsnReadFunc) ConstraintChoiceSetAsnRead, (AsnWriteFunc) ConstraintChoiceSetAsnWrite);
  }
  return t;
}


NLM_EXTERN ValNodePtr TabColumnConfigListFree (ValNodePtr columns)
{
  ValNodePtr vnp_next;

  while (columns != NULL) {
    vnp_next = columns->next;
    columns->data.ptrvalue = TabColumnConfigFree (columns->data.ptrvalue);
    columns->next = NULL;
    columns = ValNodeFree (columns);
    columns = vnp_next;
  }
  return columns;
}


NLM_EXTERN ValNodePtr TabColumnConfigListCopy (ValNodePtr orig)
{
  ValNodePtr new_list = NULL;
  TabColumnConfigPtr t;

  while (orig != NULL) {
    t = TabColumnConfigCopy (orig->data.ptrvalue);
    ValNodeAddPointer (&new_list, 0, t);
    orig = orig->next;
  }
  return new_list;
}



/* This checks the column names and returns a list of the feature fields */
NLM_EXTERN ValNodePtr ValidateFeatureFieldColumnNames (ValNodePtr header_line, ValNodePtr PNTR perr_list)
{
  ValNodePtr         header_vnp;
  ValNodePtr         err_list = NULL, col_list = NULL;
  Boolean            rval = TRUE;
  TabColumnConfigPtr t;
  FeatureFieldPtr    field;
  Int4               featqual, feat_type;
  CharPtr            first_space;
  
  if (header_line == NULL)
  {
    return FALSE;
  }
  
  header_vnp = header_line->data.ptrvalue;
  if (header_vnp == NULL || header_vnp->next == NULL)
  {
    return FALSE;
  }
  
  /* skip ID column */
  header_vnp = header_vnp->next;
  while (header_vnp != NULL && rval)
  {
    first_space = StringChr (header_vnp->data.ptrvalue, ' ');
    if (first_space != NULL) {
      *first_space = 0;
      feat_type = GetFeatureTypeByName (header_vnp->data.ptrvalue);
      featqual = GetFeatQualByName (first_space + 1);
      *first_space = ' ';
      if (feat_type < 0 || featqual < 0) {
        /* unable to recognize column name */
        ValNodeAddPointer (&err_list, 0, StringSave (header_vnp->data.ptrvalue));
        /* if we're not able to send back a list of errors, just quit now */
        if (perr_list == NULL) {
          rval = FALSE;
        }
      } else if (err_list == NULL) {
        /* if we've already found errors, don't bother collecting more fields */
        field = FeatureFieldNew ();
        field->type = feat_type;
        field->field = ValNodeNew (NULL);
        field->field->choice = FeatQualChoice_legal_qual;
        field->field->data.intvalue = featqual;
        t = TabColumnConfigNew ();
        t->field = ValNodeNew (NULL);
        t->field->choice = FieldType_feature_field;
        t->field->data.ptrvalue = field;
        ValNodeAddPointer (&col_list, 0, t);
      }
    } else {
      featqual = GetFeatQualByName (header_vnp->data.ptrvalue);
      if (featqual < 0) {
        /* unable to recognize column name */
        ValNodeAddPointer (&err_list, 0, StringSave (header_vnp->data.ptrvalue));
        /* if we're not able to send back a list of errors, just quit now */
        if (perr_list == NULL) {
          rval = FALSE;
        }
      } else if (err_list == NULL) {
        /* if we've already found errors, don't bother collecting more fields */
        field = FeatureFieldNew ();
        field->type = Feature_type_any;
        field->field = ValNodeNew (NULL);
        field->field->choice = FeatQualChoice_legal_qual;
        field->field->data.intvalue = featqual;
        t = TabColumnConfigNew ();
        t->field = ValNodeNew (NULL);
        t->field->choice = FieldType_feature_field;
        t->field->data.ptrvalue = field;
        ValNodeAddPointer (&col_list, 0, t);
      }
    }
    header_vnp = header_vnp->next;
  }
  if (err_list != NULL) {
    col_list = TabColumnConfigListFree (col_list);
    if (perr_list != NULL) {
      *perr_list = err_list;
    } else {
      err_list = ValNodeFreeData (err_list);
    }
  }
  return col_list;
}

typedef struct findgenelocustag {
  CharPtr locus_tag;
  ValNodePtr gene_list;
} FindGeneLocusTagData, PNTR FindGeneLocusTagPtr;

static void FindGeneByLocusTagBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  FindGeneLocusTagPtr p;
  SeqFeatPtr          gene;
  SeqMgrFeatContext   fcontext;

  if (bsp == NULL || userdata == NULL || !ISA_na (bsp->mol)) {
    return;
  }

  p = (FindGeneLocusTagPtr) userdata;

  gene = SeqMgrGetGeneByLocusTag (bsp, p->locus_tag, &fcontext);
  if (gene != NULL) {
    ValNodeAddPointer (&p->gene_list, OBJ_SEQFEAT, gene);
  }
}


typedef struct objbymatch {
  ValNodePtr obj_list;
  StringConstraintPtr scp;
} ObjByMatchData, PNTR ObjByMatchPtr;

static void GetFeaturesByDbxrefCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ObjByMatchPtr p;
  ValNodePtr    vnp;
  DbtagPtr      dbt;
  Char          buf[20];
  Boolean       found = FALSE;

  if (sfp == NULL || sfp->dbxref == NULL || userdata == NULL) return;
  p = (ObjByMatchPtr) userdata;

  if (IsStringConstraintEmpty (p->scp)) return;

  for (vnp = sfp->dbxref; vnp != NULL && !found; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && dbt->tag != NULL) {
      if (dbt->tag->id > 0) {
        sprintf (buf, "%d", dbt->tag->id);
        if (DoesStringMatchConstraint (buf, p->scp)) {
          found = TRUE;
        }
      } else if (DoesStringMatchConstraint (dbt->tag->str, p->scp)) {
        found = TRUE;
      }
    }
  }
  if (found) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQFEAT, sfp);
  }
}


static ValNodePtr GetFeaturesByDbxref (SeqEntryPtr sep, CharPtr dbxref, Uint1 match_location)
{
  ObjByMatchData d;

  d.scp = StringConstraintNew ();
  d.scp->match_text = StringSave (dbxref);
  d.scp->match_location = match_location;
  d.obj_list = NULL;
  VisitFeaturesInSep (sep, &d, GetFeaturesByDbxrefCallback);
  d.scp = StringConstraintFree (d.scp);
  return d.obj_list;
}


static void GetBioSourcesByTaxNameDescriptorCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjByMatchPtr p;
  BioSourcePtr  biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;
  p = (ObjByMatchPtr) userdata;

  if (IsStringConstraintEmpty (p->scp)) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop != NULL && biop->org != NULL && DoesStringMatchConstraint (biop->org->taxname, p->scp)) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQDESC, sdp);
  }

}


static void GetBioSourcesByTaxNameFeatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ObjByMatchPtr p;
  BioSourcePtr biop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;
  p = (ObjByMatchPtr) userdata;

  if (IsStringConstraintEmpty (p->scp)) return;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop != NULL && biop->org != NULL && DoesStringMatchConstraint (biop->org->taxname, p->scp)) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQFEAT, sfp);
  }

}


static ValNodePtr GetBioSourcesByTaxName (SeqEntryPtr sep, CharPtr taxname, Uint1 match_location)
{
  ObjByMatchData d;

  d.scp = StringConstraintNew ();
  d.scp->match_text = StringSave (taxname);
  d.scp->match_location = match_location;
  d.obj_list = NULL;
  VisitDescriptorsInSep (sep, &d, GetBioSourcesByTaxNameDescriptorCallback);

  VisitFeaturesInSep (sep, &d, GetBioSourcesByTaxNameFeatureCallback);
  d.scp = StringConstraintFree (d.scp);
  return d.obj_list;
}


typedef struct objbystrinfld {
  ValNodePtr obj_list;
  FieldTypePtr field;
  StringConstraintPtr scp;
} ObjByStrInFldData, PNTR ObjByStrInFldPtr;


static void GetBioSourcesBySourceQualDescriptorCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjByStrInFldPtr p;
  CharPtr      tmp;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;
  p = (ObjByStrInFldPtr) userdata;

  if (IsStringConstraintEmpty (p->scp)) return;

  tmp = GetFieldValueForObject (OBJ_SEQDESC, sdp, p->field, p->scp);
  if (tmp != NULL) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQDESC, sdp);
  }
  tmp = MemFree (tmp);
}


static void GetBioSourcesBySourceQualFeatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ObjByStrInFldPtr p;
  CharPtr          tmp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;
  p = (ObjByStrInFldPtr) userdata;

  if (IsStringConstraintEmpty (p->scp)) return;

  tmp = GetFieldValueForObject (OBJ_SEQFEAT, sfp, p->field, p->scp);
  if (tmp != NULL) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQFEAT, sfp);
  }
  tmp = MemFree (tmp);
}


static ValNodePtr GetBioSourcesBySourceQual (SeqEntryPtr sep, SourceQualChoicePtr q, CharPtr val, Uint1 match_location)
{
  ObjByStrInFldData od;

  od.scp = StringConstraintNew();
  od.scp->match_text = StringSave (val);
  od.scp->match_location = match_location;
  od.obj_list = NULL;
  od.field = ValNodeNew (NULL);
  od.field->choice = FieldType_source_qual;
  od.field->data.ptrvalue = q;

  VisitDescriptorsInSep (sep, &od, GetBioSourcesBySourceQualDescriptorCallback);

  VisitFeaturesInSep (sep, &od, GetBioSourcesBySourceQualFeatureCallback);

  od.field = ValNodeFree (od.field);
  od.scp = StringConstraintFree (od.scp);
  return od.obj_list;  
}


static void GetBioseqsByIdCallback (BioseqPtr bsp, Pointer data)
{
  ObjByMatchPtr d;
  ObjectIdPtr   oip;
  SeqIdPtr      sip;
  Boolean       found_match = FALSE;
  DbtagPtr      dbtag;
  CharPtr       cp, tmp_id;

  if (bsp == NULL || data == NULL || (d = (ObjByMatchPtr) data) == NULL) {
    return;
  }

  found_match = DoesSeqIDListMeetStringConstraint (bsp->id, d->scp);

  for (sip = bsp->id; sip != NULL && !found_match; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL && sip->data.ptrvalue != NULL) {
      dbtag = (DbtagPtr) sip->data.ptrvalue;
      if (StringCmp (dbtag->db, "NCBIFILE") == 0 && dbtag->tag != NULL) {
        if (DoesStringMatchConstraint (dbtag->tag->str, d->scp)) {
          found_match = TRUE;
        } else if ((cp = StringRChr (dbtag->tag->str, '/')) != NULL) {
          tmp_id = (CharPtr) MemNew (sizeof (Char) * (cp - dbtag->tag->str + 1));
          StringNCpy (tmp_id, dbtag->tag->str, cp - dbtag->tag->str);
          tmp_id[cp - dbtag->tag->str] = 0;
          if (DoesStringMatchConstraint (tmp_id, d->scp)) {
            found_match = TRUE;
          }
          tmp_id = MemFree (tmp_id);
        }
      }
    } else if (sip->choice == SEQID_LOCAL && (oip = sip->data.ptrvalue) != NULL
               && StringNICmp (oip->str, "bankit", 6) == 0
               && DoesStringMatchConstraint (oip->str + 6, d->scp)) {
      found_match = TRUE;
    }
  }
  if (found_match) {
    ValNodeAddPointer (&(d->obj_list), OBJ_BIOSEQ, bsp);
  }
}


static ValNodePtr FindBioseqsByMatchType (SeqEntryPtr sep, Uint1 match_location, CharPtr match_str)
{
  ObjByMatchData d;

  if (sep == NULL || StringHasNoText (match_str)) {
    return NULL;
  }
  d.scp = StringConstraintNew ();
  d.scp->match_text = StringSave (match_str);
  d.scp->match_location = match_location;
  d.obj_list = NULL;
  VisitBioseqsInSep (sep, &d, GetBioseqsByIdCallback);
  d.scp = StringConstraintFree (d.scp);
  return d.obj_list;
}


static ValNodePtr 
FindMatchForRow 
(MatchTypePtr match_type,
 CharPtr      match_str,
 Uint2        entityID,
 SeqEntryPtr  sep)
{
  ValNodePtr match_list = NULL;
  FindGeneLocusTagData fd;
  SeqFeatPtr           sfp;
  SeqMgrFeatContext    fcontext;

  if (match_type == NULL || sep == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      sfp = SeqMgrGetFeatureByFeatID (entityID, NULL, match_str, NULL, &fcontext);
      if (sfp != NULL) {
        ValNodeAddPointer (&match_list, OBJ_SEQFEAT, sfp);
      }
      break;
    case eTableMatchGeneLocusTag:
      fd.locus_tag = match_str;
      fd.gene_list = NULL;
      VisitBioseqsInSep (sep, &fd, FindGeneByLocusTagBioseqCallback);
      ValNodeLink (&match_list, fd.gene_list);
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      ValNodeLink (&match_list, FindBioseqsByMatchType (sep, match_type->match_location, match_str));
      break;
    case eTableMatchDbxref:
      match_list = GetFeaturesByDbxref (sep, match_str, match_type->match_location);
      break;
    case eTableMatchBioSource:
      match_list = GetBioSourcesByTaxName (sep, match_str, match_type->match_location);
      break;
    case eTableMatchSourceQual:
      match_list = GetBioSourcesBySourceQual (sep, match_type->data, match_str, match_type->match_location);
      break;
  }
  return match_list;
}


static ValNodePtr GetFeatureListForProteinBioseq (Uint1 featdef, BioseqPtr bsp)
{
  ValNodePtr feat_list = NULL;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;
  Int4              seqfeattype;

  if (bsp == NULL || !ISA_aa (bsp->mol)) 
  {
    return NULL;
  }

  seqfeattype = FindFeatFromFeatDefType (featdef);
  if (seqfeattype == SEQFEAT_PROT)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  else
  {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) 
    {
      if (featdef == FEATDEF_CDS)
      {
        sfp = cds;
      }
      else if (featdef == FEATDEF_GENE)
      {
        sfp = GetGeneForFeature (cds);
      }
      else if (featdef == FEATDEF_mRNA)
      {
        sfp = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
      }
      if (sfp != NULL)
      {
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
      }
    }
  }
  return feat_list;
}


static ValNodePtr GetFeatureListForNucleotideBioseq (Uint1 featdef, BioseqPtr bsp)
{
  ValNodePtr feat_list = NULL;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  Int4              seqfeattype;
  BioseqPtr         prot_bsp;

  if (bsp == NULL || ISA_aa (bsp->mol)) 
  {
    return NULL;
  }

  seqfeattype = FindFeatFromFeatDefType (featdef);
  if (seqfeattype == SEQFEAT_PROT)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &fcontext))
    {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      ValNodeLink (&feat_list, GetFeatureListForProteinBioseq (featdef, prot_bsp));
    }
  }
  else
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


static ValNodePtr GetFeaturesForGene (SeqFeatPtr gene, Uint1 featdef)
{
  BioseqPtr bsp;
  SeqFeatPtr sfp;
  ValNodePtr feat_list = NULL;
  SeqMgrFeatContext fcontext;
  Int4              start, stop, swap;

  if (gene == NULL) return NULL;

  bsp = BioseqFindFromSeqLoc (gene->location);
  start = SeqLocStart (gene->location);
  stop = SeqLocStop (gene->location);
  if (stop < start) 
  {
    swap = start;
    start = stop;
    stop = swap;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       sfp != NULL && fcontext.left < stop;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
  {
    if (fcontext.right >= start && gene == GetGeneForFeature (sfp))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


static ValNodePtr GetFeatureListForGene (Uint1 featdef, SeqFeatPtr gene)
{
  ValNodePtr feat_list = NULL, cds_list, vnp;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;
  BioseqPtr         protbsp;

  if (gene == NULL) 
  {
    return NULL;
  }

  if (featdef == FEATDEF_GENE)
  {
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, gene);
  }
  else if (FindFeatFromFeatDefType (featdef == SEQFEAT_PROT))
  {
    cds_list = GetFeaturesForGene (gene, FEATDEF_CDS);
    for (vnp = cds_list; vnp != NULL; vnp = vnp->next) 
    {
      cds = vnp->data.ptrvalue;
      if (cds != NULL)
      {
        protbsp = BioseqFindFromSeqLoc (cds->product);
        for (sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, featdef, &fcontext);
             sfp != NULL;
             sfp = SeqMgrGetNextFeature (protbsp, sfp, 0, featdef, &fcontext))
        {
          ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
        }
      }
    }
    cds_list = ValNodeFree (cds_list);
  }
  else
  {
    feat_list = GetFeaturesForGene (gene, featdef);
  }

  return feat_list;
}


static ValNodePtr AddFeaturesFromBioseqSet (BioseqSetPtr bssp, Uint1 featdef)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  Int4        seqfeattype;
  ValNodePtr  item_list = NULL;

  if (bssp == NULL) return NULL;

  seqfeattype = FindFeatFromFeatDefType (featdef);
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (sep->data.ptrvalue == NULL) continue;
    if (IS_Bioseq (sep)) {
      bsp = sep->data.ptrvalue;
      if (seqfeattype == SEQFEAT_PROT) {
        if (ISA_aa (bsp->mol)) {
          ValNodeLink (&item_list, GetFeatureListForProteinBioseq (featdef, bsp));
        }
      } else if (!ISA_aa (bsp->mol)) {
        ValNodeLink (&item_list, GetFeatureListForNucleotideBioseq (featdef, bsp));
      }
    } else if (IS_Bioseq_set (sep)) {
      ValNodeLink (&item_list, AddFeaturesFromBioseqSet (sep->data.ptrvalue, featdef));
    }
  }
  return item_list;
}


static ValNodePtr GetFeatureListForBioSourceObjects (ValNodePtr item_list, FeatureFieldPtr field)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr   bsp;
  ObjValNodePtr ovp;
  ValNodePtr  feature_list = NULL;

  if (item_list == NULL || field == NULL) return NULL;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = vnp->data.ptrvalue;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), bsp));
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = vnp->data.ptrvalue;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          ValNodeLink (&feature_list, AddFeaturesFromBioseqSet (ovp->idx.parentptr, GetFeatdefFromFeatureType(field->type)));
        } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), bsp));
        }
      }
    }
  }
  return feature_list;
}


NLM_EXTERN ValNodePtr ValNodeCopyPtr (ValNodePtr orig)
{
  ValNodePtr new_list = NULL, last_vnp = NULL, vnp;

  while (orig != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = orig->choice;
    vnp->data.ptrvalue = orig->data.ptrvalue;
    if (last_vnp == NULL) {
      new_list = vnp;
    } else {
      last_vnp->next = vnp;
    }
    last_vnp = vnp;
    orig = orig->next;
  }
  return new_list;
}


static ValNodePtr GetFeatureListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list, FeatureFieldPtr field)
{
  ValNodePtr feature_list = NULL, vnp;

  if (match_list == NULL || field == NULL || match_type == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      feature_list = ValNodeCopyPtr (match_list);
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForGene (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchProteinID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForProteinBioseq (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchDbxref:
      feature_list = ValNodeCopyPtr (match_list);
      break;
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchBioSource:
    case eTableMatchSourceQual:
      ValNodeLink (&feature_list, GetFeatureListForBioSourceObjects (match_list, field));
      break;
  }
  return feature_list;
}


static void AddBioSourcesForBioseq (BioseqPtr bsp, ValNodePtr PNTR feature_list)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  if (feature_list == NULL) return;
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
        sdp != NULL;
        sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &context)) {
    ValNodeAddPointer (feature_list, OBJ_SEQDESC, sdp);
  }
}

static void AddBioSourcesForFeature (SeqFeatPtr sfp, ValNodePtr PNTR feature_list)
{
  BioseqPtr bsp;

  if (sfp == NULL || feature_list == NULL) return;

  if (sfp->data.choice == SEQFEAT_BIOSRC) {
    ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
  } else {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    AddBioSourcesForBioseq (bsp, feature_list);
  }
}


static ValNodePtr GetBioSourceListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list, FeatureFieldPtr field)
{
  ValNodePtr feature_list = NULL, vnp;

  if (match_list == NULL || field == NULL || match_type == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_BIOSEQ) {
          AddBioSourcesForBioseq (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchDbxref:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchBioSource:
    case eTableMatchSourceQual:
      feature_list = ValNodeCopyPtr (match_list);
      break;
  }
  return feature_list;
}


static void AddPubsForBioseq (BioseqPtr bsp, ValNodePtr PNTR feature_list)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  SeqFeatPtr  sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || feature_list == NULL) return;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext)) {
    ValNodeAddPointer (feature_list, OBJ_SEQDESC, sdp);
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
       sfp != NULL; 
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext)) {
    ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
  }
}


static ValNodePtr AddPubListFromBioseqSet (BioseqSetPtr bssp)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ValNodePtr  item_list = NULL;

  if (bssp == NULL) return NULL;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (sep->data.ptrvalue == NULL) continue;
    if (IS_Bioseq (sep)) {
      bsp = sep->data.ptrvalue;
      if (!ISA_aa (bsp->mol)) {
        AddPubsForBioseq (bsp, &item_list);
      }
    } else if (IS_Bioseq_set (sep)) {
      ValNodeLink (&item_list, AddPubListFromBioseqSet (sep->data.ptrvalue));
    }
  }
  return item_list;
}


static ValNodePtr GetPubListForBioSourceObjects (ValNodePtr item_list)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr   bsp;
  ObjValNodePtr ovp;
  ValNodePtr  feature_list = NULL;

  if (item_list == NULL) return NULL;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = vnp->data.ptrvalue;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        AddPubsForBioseq (bsp, &feature_list);
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = vnp->data.ptrvalue;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          ValNodeLink (&feature_list, AddPubListFromBioseqSet (ovp->idx.parentptr));
        } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          AddPubsForBioseq (bsp, &feature_list);
        }
      }
    }
  }
  return feature_list;
}


static ValNodePtr GetPubListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list)
{
  SeqFeatPtr sfp;
  ValNodePtr vnp;
  ValNodePtr feature_list = NULL;

  if (match_type == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          AddPubsForBioseq (BioseqFindFromSeqLoc (sfp->location), &feature_list);
        }
      }
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          AddPubsForBioseq (BioseqFindFromSeqLoc (sfp->location), &feature_list);
        }
      }
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_BIOSEQ) {
          AddPubsForBioseq (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchDbxref:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          AddPubsForBioseq (BioseqFindFromSeqLoc (sfp->location), &feature_list);
        }
      }
      break;
    case eTableMatchBioSource:
    case eTableMatchSourceQual:
      feature_list = GetPubListForBioSourceObjects (match_list);
      break;
  }
  return feature_list;
}


static ValNodePtr GetSequenceListForBioSourceObjects (ValNodePtr item_list)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr   bsp;
  ObjValNodePtr ovp;
  ValNodePtr  seq_list = NULL;
  SeqEntryPtr sep;

  if (item_list == NULL) return NULL;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = vnp->data.ptrvalue;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, bsp);
        }
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = vnp->data.ptrvalue;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          sep = SeqMgrGetSeqEntryForData (ovp->idx.parentptr);
          VisitBioseqsInSep (sep, &seq_list, CollectNucBioseqCallback);
        } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          if (bsp != NULL) {
            ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, bsp);
          }
        }
      }
    }
  }
  return seq_list;
}


static ValNodePtr GetSequenceListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list)
{
  SeqFeatPtr sfp;
  ValNodePtr vnp;
  ValNodePtr seq_list = NULL;
  BioseqPtr  bsp;

  if (match_type == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          bsp = BioseqFindFromSeqLoc (sfp->location);
          if (bsp != NULL) {
            ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, bsp);
          }
        }
      }
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          bsp = BioseqFindFromSeqLoc (sfp->location);
          if (bsp != NULL) {
            ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, bsp);
          }
        }
      }
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_BIOSEQ) {
          ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, vnp->data.ptrvalue);
        }
      }
      break;
    case eTableMatchDbxref:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          sfp = (SeqFeatPtr) vnp->data.ptrvalue;
          bsp = BioseqFindFromSeqLoc (sfp->location);
          if (bsp != NULL) {
            ValNodeAddPointer (&seq_list, OBJ_BIOSEQ, bsp);
          }
        }
      }
      break;
    case eTableMatchBioSource:
    case eTableMatchSourceQual:
      seq_list = GetSequenceListForBioSourceObjects (match_list);
      break;
  }
  return seq_list;
}


static ValNodePtr GetStructuredCommentListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list)
{
  ValNodePtr seq_list, target_list = NULL, vnp;
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  seq_list = GetSequenceListForRowAndColumn (match_type, match_list);

  for (vnp = seq_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_BIOSEQ) {
      for (sdp = SeqMgrGetNextDescriptor (vnp->data.ptrvalue, NULL, Seq_descr_user, &context);
           sdp != NULL;
           sdp = sdp->next) {
        if (IsUserObjectStructuredComment (sdp->data.ptrvalue)) {
          ValNodeAddPointer (&target_list, OBJ_SEQDESC, sdp);
        }
      }
    }
  }    
  seq_list = ValNodeFree (seq_list);
  return target_list;
}


static ValNodePtr GetTargetListForRowAndColumn (MatchTypePtr match_type, ValNodePtr match_list, FieldTypePtr field, ValNodePtr constraint)
{
  ValNodePtr target_list = NULL, vnp_prev = NULL, vnp, vnp_next, tmp_list;
  FeatureFieldPtr feature_field;

  if (field == NULL || match_type == NULL) return NULL;
  switch (field->choice) {
    case FieldType_source_qual:
      target_list = GetBioSourceListForRowAndColumn (match_type, match_list, field->data.ptrvalue);
      break;
    case FieldType_feature_field:
      target_list = GetFeatureListForRowAndColumn (match_type, match_list, field->data.ptrvalue);
      break;
    case FieldType_cds_gene_prot:
      feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
      target_list = GetFeatureListForRowAndColumn (match_type, match_list, feature_field);
      feature_field = FeatureFieldFree (feature_field);
      break;
    case FieldType_pub:
      target_list = GetPubListForRowAndColumn (match_type, match_list);
      break;
    case FieldType_rna_field:
      feature_field = FeatureFieldFromRnaQual (field->data.ptrvalue);
      target_list = GetFeatureListForRowAndColumn (match_type, match_list, feature_field);
      feature_field = FeatureFieldFree (feature_field);
      break;
    case FieldType_struc_comment_field:
      target_list = GetStructuredCommentListForRowAndColumn (match_type, match_list);
      break;
    case FieldType_misc:
      if (field->data.intvalue == Misc_field_genome_project_id) {
        target_list = GetSequenceListForRowAndColumn (match_type, match_list);
      } else if (field->data.intvalue == Misc_field_comment_descriptor) {
        tmp_list = GetSequenceListForRowAndColumn (match_type, match_list);
        for (vnp = tmp_list; vnp != NULL; vnp = vnp->next) {
          AddCommentDescriptorDestinationsForBioseq (vnp->data.ptrvalue, &target_list);
        }
        tmp_list = ValNodeFree (tmp_list);
      }
      break;
  }

  /* remove targets that do not match constraint */
  vnp = target_list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    if (!DoesObjectMatchConstraintChoiceSet (vnp->choice, vnp->data.ptrvalue, constraint)) {
      if (vnp_prev == NULL) {
        target_list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFree (vnp);
    } else {
      vnp_prev = vnp;
    }
    vnp = vnp_next;
  }

  return target_list;
}


static void ReportMissingTargets (ValNodePtr PNTR perr_list, FieldTypePtr ft, CharPtr match_val, Int4 col_num, Int4 line_num)
{
  CharPtr            feat_name;
  FeatureFieldPtr    field;
  CharPtr            no_feat_fmt = "No %s feature for %s (column %d, line %d)";
  CharPtr            no_src_fmt = "No biosource for %s (column %d, line %d)";
  CharPtr            no_seq_fmt = "No sequence for %s (column %d, line %d)";
  CharPtr            no_cmt_fmt = "No structured comment for %s (column %d, line %d)";
  CharPtr            err_msg;

  if (perr_list == NULL || ft == NULL || match_val == NULL) return;

  switch (ft->choice) {
    case FieldType_source_qual:
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_src_fmt) 
                                                    + StringLen (match_val)
                                                    + 30));
      sprintf (err_msg, no_src_fmt, match_val, col_num, line_num);
      ValNodeAddPointer (perr_list, 0, err_msg);
      break;
    case FieldType_feature_field:
      field = (FeatureFieldPtr) ft->data.ptrvalue;
      if (field != NULL) {
        feat_name = GetFeatureNameFromFeatureType (field->type);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                      + StringLen (feat_name)
                                                      + StringLen (match_val)
                                                      + 30));
        sprintf (err_msg, no_feat_fmt, feat_name, match_val, col_num, line_num);
        ValNodeAddPointer (perr_list, 0, err_msg);
      }
      break;
    case FieldType_cds_gene_prot:
      field = FeatureFieldFromCDSGeneProtField (ft->data.intvalue);
      if (field != NULL) {
        feat_name = GetFeatureNameFromFeatureType (field->type);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                      + StringLen (feat_name)
                                                      + StringLen (match_val)
                                                      + 30));
        sprintf (err_msg, no_feat_fmt, feat_name, match_val, col_num, line_num);
        ValNodeAddPointer (perr_list, 0, err_msg);
      }
      field = FeatureFieldFree (field);
      break;
    case FieldType_rna_field:
      field = FeatureFieldFromRnaQual (ft->data.ptrvalue);
      if (field != NULL) {
        feat_name = GetFeatureNameFromFeatureType (field->type);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                      + StringLen (feat_name)
                                                      + StringLen (match_val)
                                                      + 30));
        sprintf (err_msg, no_feat_fmt, feat_name, match_val, col_num, line_num);
        ValNodeAddPointer (perr_list, 0, err_msg);
      }
      field = FeatureFieldFree (field);
      break;
    case FieldType_struc_comment_field:
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_cmt_fmt) + StringLen (match_val) + 30));
      sprintf (err_msg, no_cmt_fmt, match_val, col_num, line_num);
      ValNodeAddPointer (perr_list, 0, err_msg);
      break;
    case FieldType_misc:
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_seq_fmt) 
                                                    + StringLen (match_val)
                                                    + 30));
      sprintf (err_msg, no_seq_fmt, match_val, col_num, line_num);
      ValNodeAddPointer (perr_list, 0, err_msg);
      break;
  }
}


static void ReportEmptyIDColumn (ValNodePtr PNTR perr_list, Int4 line_num)
{
  CharPtr            err_msg;
  CharPtr            missing_id_fmt = "No ID for line %d";

  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_id_fmt) + 15));
  sprintf (err_msg, missing_id_fmt, line_num);
  ValNodeAddPointer (perr_list, 0, err_msg);
}

static ValNodePtr FindMatchChoiceInLine (ValNodePtr val_vnp, ValNodePtr col_vnp)
{
  TabColumnConfigPtr t;

  while (val_vnp != NULL && col_vnp != NULL) {
    t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
    if (t != NULL && t->match_type != NULL) {
      return val_vnp;
    }
    val_vnp = val_vnp->next;
    col_vnp = col_vnp->next;
  }
  return NULL;
}


NLM_EXTERN SeqFeatPtr GetmRNAForFeature (SeqFeatPtr sfp)
{
  SeqMgrFeatContext fcontext;
  BioseqPtr         pbsp;

  if (sfp == NULL) return NULL;
  if (sfp->data.choice == SEQFEAT_PROT) 
  { 
    pbsp = BioseqFindFromSeqLoc (sfp->location);
    sfp = SeqMgrGetCDSgivenProduct (pbsp, NULL);
    if (sfp == NULL) return NULL;
  }
  return SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
}


NLM_EXTERN Boolean AdjustmRNAProductToMatchProteinProduct (SeqFeatPtr sfp)
{
  SeqFeatPtr mrna;
  ProtRefPtr prp;
  RnaRefPtr  rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return FALSE;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  mrna = GetmRNAForFeature (sfp);

  if (mrna == NULL) return FALSE;

  rrp = (RnaRefPtr) mrna->data.value.ptrvalue;
  if (rrp == NULL) 
  {
    rrp = RnaRefNew();
    mrna->data.value.ptrvalue = rrp;
  }

  rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
  if (prp == NULL || prp->name == NULL || StringHasNoText (prp->name->data.ptrvalue))
  {
    rrp->ext.choice = 0;
  }
  else
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave (prp->name->data.ptrvalue);
  }
  return TRUE;
}


NLM_EXTERN Boolean IsFieldTypeCDSProduct (FieldTypePtr ft)
{
  FeatureFieldPtr field;
  Boolean         rval = FALSE;

  if (ft == NULL) return FALSE;
  if (ft->choice == FieldType_feature_field) {
    field = (FeatureFieldPtr) ft->data.ptrvalue;
    if (field != NULL && field->type == Feature_type_cds
        && field->field != NULL
        && field->field->choice == FeatQualChoice_legal_qual
        && field->field->data.intvalue == Feat_qual_legal_product) {
      rval = TRUE;
    }
  } else if (ft->choice == FieldType_cds_gene_prot) {
    if (ft->data.intvalue == CDSGeneProt_field_prot_name) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean IsFieldTypeGeneLocusTag (FieldTypePtr ft)
{
  FeatureFieldPtr field;
  RnaQualPtr      rq;
  Boolean         rval = FALSE;

  if (ft == NULL) return FALSE;
  if (ft->choice == FieldType_feature_field) {
    field = (FeatureFieldPtr) ft->data.ptrvalue;
    if (field != NULL && field->type == Feature_type_gene
        && field->field != NULL
        && field->field->choice == FeatQualChoice_legal_qual
        && field->field->data.intvalue == Feat_qual_legal_locus_tag) {
      rval = TRUE;
    }
  } else if (ft->choice == FieldType_cds_gene_prot) {
    if (ft->data.intvalue == CDSGeneProt_field_gene_locus_tag) {
      rval = TRUE;
    }
  } else if (ft->choice == FieldType_rna_field) {
    rq = (RnaQualPtr) ft->data.ptrvalue;
    if (rq != NULL && rq->field == Rna_field_gene_locus_tag) {
      rval = TRUE;
    }
  }

  return rval;
}



NLM_EXTERN ValNodePtr ValidateTabTableValues (ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, col_vnp, val_vnp;
  Int4       line_num, col_num;
  TabColumnConfigPtr t;
  ValNodePtr locus_tag_values = NULL, bad_locus_tags = NULL, vnp;
  CharPtr    bad_format_fmt = "Locus tag %s has incorrect format";
  CharPtr    dup_fmt = "Locus tag %s appears in the table more than once";
  CharPtr    inconsistent_fmt = "Locus tag prefix for %s is inconsistent";
  CharPtr    err_msg;

  if (table == NULL || columns == NULL) {
    return NULL;
  }

  for (line_vnp = table, line_num = 1;
       line_vnp != NULL;
       line_vnp = line_vnp->next, line_num++) {
    for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
         val_vnp != NULL && col_vnp != NULL;
         val_vnp = val_vnp->next, col_vnp = col_vnp->next, col_num++) {
      t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
      if (t == NULL || t->match_type != NULL || val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)) {
        continue;
      }
      if (IsFieldTypeGeneLocusTag (t->field)) {
        ValNodeAddPointer (&locus_tag_values, 0, val_vnp->data.ptrvalue);
      }
    }
  }

  bad_locus_tags = FindBadLocusTagsInList (locus_tag_values);
  for (vnp = bad_locus_tags; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case eLocusTagErrorBadFormat:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_format_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, bad_format_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
      case eLocusTagErrorDuplicate:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, dup_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
      case eLocusTagErrorInconsistentPrefix:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (inconsistent_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, inconsistent_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
    }
  }
  locus_tag_values = ValNodeFree (locus_tag_values);
  return err_list;
}


NLM_EXTERN ValNodePtr GetSequenceListsForMatchTypeInTabTable (SeqEntryPtr sep, ValNodePtr table, Int4 col, MatchTypePtr match_type, ValNodePtr PNTR p_err_list)
{
  ValNodePtr vnp_row, vnp;
  ValNodePtr sequence_lists = NULL, match_list, target_list;
  Uint2      entityID;
  Int4       num, line;
  CharPtr    no_match_fmt = "No match for %s, line %d";
  CharPtr    no_match_txt_fmt = "No match text for line %d";
  CharPtr    msg;


  if (sep == NULL || table == NULL || match_type == NULL || col < 0) {
    return NULL;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (vnp_row = table, line = 1; vnp_row != NULL; vnp_row = vnp_row->next, line++) {
    vnp = vnp_row->data.ptrvalue;
    num = 0;
    while (vnp != NULL && num < col) {
      vnp = vnp->next;
      num++;
    }
    if (vnp == NULL || StringHasNoText (vnp->data.ptrvalue)) {
      ValNodeAddPointer (&sequence_lists, 0, NULL);
      if (p_err_list != NULL) {
        msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_txt_fmt) + 15));
        sprintf (msg, no_match_txt_fmt, line);
        ValNodeAddPointer (p_err_list, 0, msg);
      }
    } else {
      match_list = FindMatchForRow (match_type, vnp->data.ptrvalue, entityID, sep);
      target_list = GetSequenceListForRowAndColumn (match_type, match_list);
      match_list = ValNodeFree (match_list);
      ValNodeAddPointer (&sequence_lists, 0, target_list);
      if (target_list == NULL && p_err_list != NULL) {
        msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (vnp->data.ptrvalue) + 15));
        sprintf (msg, no_match_fmt, vnp->data.ptrvalue, line);
        ValNodeAddPointer (p_err_list, 0, msg);
      }
    }
  }
  return sequence_lists;
}


NLM_EXTERN ValNodePtr FreeSequenceLists (ValNodePtr lists)
{
  ValNodePtr vnp;

  for (vnp = lists; vnp != NULL; vnp = vnp->next) {
    vnp->data.ptrvalue = ValNodeFree (vnp->data.ptrvalue);
  }
  lists = ValNodeFree (lists);
  return lists;
}


static ValNodePtr ReportTableSummaryLine (Int4 err_lines, Int4 total_lines, CharPtr fmt)
{
  CharPtr str;
  ValNodePtr vnp;

  str = (CharPtr) MemNew (sizeof (Char) + (StringLen (fmt) + 30));
  sprintf (str, fmt, err_lines, total_lines);
  vnp = ValNodeNew (NULL);
  vnp->data.ptrvalue = str;
  return vnp;
}


NLM_EXTERN ValNodePtr GetObjectTableForTabTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr PNTR p_err_list)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, val_vnp, col_vnp, err_vnp;
  ValNodePtr obj_table = NULL, obj_row;
  Int4       line_num = 1, col_num;
  Uint2      entityID;
  ValNodePtr match_list, match_choice, target_list;
  TabColumnConfigPtr t;
  CharPtr            err_msg;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  MatchTypePtr       match_type;
  Int4       num_empty = 0, num_missing = 0, num_no_targets = 0;


  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    if (p_err_list == NULL) {
      err_list = ValNodeFreeData (err_list);
    } else {
      *p_err_list = err_list;
    }
    return NULL;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  match_type = FindMatchTypeInHeader (columns);
  if (match_type == NULL) return NULL;

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    obj_row = NULL;
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
      num_empty++;
    } else {
      match_list = FindMatchForRow (match_type, match_choice->data.ptrvalue, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
        num_missing ++;
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          target_list = NULL;
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type != NULL 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            /* no targets */
          } else {         
            target_list = GetTargetListForRowAndColumn (match_type, match_list, t->field, t->constraint);
            if (target_list == NULL) {
              ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
              num_no_targets++;
            }
          }
          ValNodeAddPointer (&obj_row, 0, target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
    }
    ValNodeAddPointer (&obj_table, 0, obj_row);
  }

  match_type = MatchTypeFree (match_type);

  if (err_list != NULL) {    
    if (num_empty > 0) {
      err_vnp = ReportTableSummaryLine (num_empty, line_num - 1, "%d lines out of %d have no ID value");
      err_vnp->next = err_list;
      err_list = err_vnp;
    }
    if (num_no_targets > 0) {
      err_vnp = ReportTableSummaryLine (num_no_targets, line_num - 1, "%d lines out of %d have no targets");
      err_vnp->next = err_list;
      err_list = err_vnp;
    }
    if (num_missing > 0) {
      err_vnp = ReportTableSummaryLine (num_missing, line_num - 1, "%d lines out of %d have no match");
      err_vnp->next = err_list;
      err_list = err_vnp;
    }
    if (p_err_list == NULL) {
      err_list = ValNodeFreeData (err_list);
    } else {
      *p_err_list = err_list;
    }
  }  
  return obj_table;
}


NLM_EXTERN ValNodePtr FreeObjectTableForTabTable (ValNodePtr table)
{
  ValNodePtr vnp_next, vnp_row, vnp_row_next;

  while (table != NULL) {
    vnp_next = table->next;
    table->next = NULL;
    vnp_row = table->data.ptrvalue;
    while (vnp_row != NULL) {
      vnp_row_next = vnp_row->next;
      vnp_row->next = NULL;
      vnp_row->data.ptrvalue = ValNodeFree (vnp_row->data.ptrvalue);
      vnp_row = ValNodeFree (vnp_row);
      vnp_row = vnp_row_next;
    }
    table = ValNodeFree (table);
    table = vnp_next;
  }
  return table;
}


typedef struct countfeat {
  Uint1 featdef;
  Int4 num;
} CountFeatData, PNTR CountFeatPtr;


static void CountFeaturesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CountFeatPtr p;

  if (sfp == NULL || userdata == NULL) return;

  p = (CountFeatPtr) userdata;
  if (sfp->idx.subtype == p->featdef) {
    p->num++;
  }
}

static void CountBioSourceDescriptorsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  Int4Ptr p;

  p = (Int4Ptr) userdata;
  if (sdp != NULL && p != NULL && sdp->choice == Seq_descr_source) {
    (*p)++;
  }
}


static void CountPubDescriptorsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  Int4Ptr p;

  p = (Int4Ptr) userdata;
  if (sdp != NULL && p != NULL && sdp->choice == Seq_descr_pub) {
    (*p)++;
  }
}


static ValNodePtr CountObjectsForColumnFields (SeqEntryPtr sep, ValNodePtr columns)
{
  ValNodePtr count_list = NULL, vnp;
  TabColumnConfigPtr t;
  CountFeatData d;
  FeatureFieldPtr f;
  Int4 num;
  Uint1 featdef = 0;
  ValNodePtr tmp_list = NULL;

  d.featdef = 0;
  d.num = 0;
  for (vnp = columns; vnp != NULL; vnp = vnp->next) {
    num = 0;
    t = (TabColumnConfigPtr) vnp->data.ptrvalue;
    if (t != NULL && t->match_type == NULL && t->field != NULL) {
      switch (t->field->choice) {
        case FieldType_source_qual:
          if (featdef != FEATDEF_BIOSRC) {
            d.featdef = FEATDEF_BIOSRC;
            d.num = 0;
            VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            VisitDescriptorsInSep (sep, &(d.num), CountBioSourceDescriptorsCallback);
          }
          num = d.num;
          break;
        case FieldType_feature_field:
          f = (FeatureFieldPtr) t->field->data.ptrvalue;
          if (f != NULL) {
            featdef = GetFeatdefFromFeatureType(f->type);
            if (featdef != d.featdef) {
              d.featdef = featdef;
              d.num = 0;
              VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            }
            num = d.num;
          }
          break;
        case FieldType_cds_gene_prot:
          f = FeatureFieldFromCDSGeneProtField (t->field->data.intvalue);
          if (f != NULL) {
            featdef = GetFeatdefFromFeatureType(f->type);
            if (featdef != d.featdef) {
              d.featdef = featdef;
              d.num = 0;
              VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            }
            num = d.num;
          }
          f = FeatureFieldFree (f);
          break;
        case FieldType_rna_field:
          f = FeatureFieldFromRnaQual (t->field->data.ptrvalue);
          if (f != NULL) {
            featdef = GetFeatdefFromFeatureType(f->type);
            if (featdef != d.featdef) {
              d.featdef = featdef;
              d.num = 0;
              VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            }
            num = d.num;
          }
          f = FeatureFieldFree (f);
          break;
        case FieldType_pub:
          d.featdef = FEATDEF_PUB;
          d.num = 0;
          VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
          VisitDescriptorsInSep (sep, &(d.num), CountPubDescriptorsCallback);
          num = d.num;
          break;
        case FieldType_struc_comment_field:
          VisitDescriptorsInSep (sep, &tmp_list, CollectStructuredCommentsCallback);
          num = ValNodeLen (tmp_list);
          tmp_list = ValNodeFree (tmp_list);
          break;
        case FieldType_misc:
          if (t->field->data.intvalue == Misc_field_genome_project_id) {
            VisitBioseqsInSep (sep, &tmp_list, CollectNucBioseqCallback);
            num = ValNodeLen (tmp_list);
            tmp_list = ValNodeFree (tmp_list);
          } else if (t->field->data.intvalue == Misc_field_comment_descriptor) {
            tmp_list = CollectCommentDescriptors (sep);
            num = ValNodeLen (tmp_list);
            tmp_list = ValNodeFree (tmp_list);
          }
          break;
      }
    }
    ValNodeAddInt (&count_list, 0, num);
  }
  return count_list;
}


NLM_EXTERN ValNodePtr ApplyTableValuesToObjectTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table)
{
  ValNodePtr val_line_vnp, obj_line_vnp;
  ValNodePtr val_vnp, obj_vnp, col_vnp;
  ValNodePtr target_vnp;
  TabColumnConfigPtr t;
  CharPtr val, qual_name;
  ValNodePtr         err_list = NULL, count_list, count_affected_list = NULL, count_vnp, count_tot_vnp;
  CharPtr            err_msg;
  CharPtr            bad_col_val_fmt = "Did not set value for column %d, line %d";
  CharPtr            num_affected_fmt = "%d fields affected";
  CharPtr            col_num_affected_fmt = "For %s (column %d), %d items were affected out of %d total";
  Int4 num_fields_affected = 0, col_num, line_num, num_this_column;
  Boolean success;
  ValNodePtr count_msg = NULL;

  count_list = CountObjectsForColumnFields (sep, columns);

  for (val_line_vnp = table, obj_line_vnp = obj_table, line_num = 1;
       val_line_vnp != NULL && obj_line_vnp != NULL;
       val_line_vnp = val_line_vnp->next, obj_line_vnp = obj_line_vnp->next, line_num++) {
    val_vnp = val_line_vnp->data.ptrvalue;
    obj_vnp = obj_line_vnp->data.ptrvalue;
    col_vnp = columns;
    col_num = 1;
    count_vnp = count_affected_list;
    while (obj_vnp != NULL && col_vnp != NULL) {
      num_this_column = 0;
      if (obj_vnp->data.ptrvalue != NULL) {
        t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
        if (t == NULL || t->match_type != NULL 
            || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
          /* ignore column or skip blank value */
        } else {
          if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
            val = "";
          } else {
            val = val_vnp->data.ptrvalue;
          }
          for (target_vnp = obj_vnp->data.ptrvalue; target_vnp != NULL; target_vnp = target_vnp->next) {
            if (val[0] == 0) {
              success = RemoveFieldValueForObject (target_vnp->choice, target_vnp->data.ptrvalue, t->field, NULL);
            } else {
              success = SetFieldValueForObject (target_vnp->choice, target_vnp->data.ptrvalue, t->field, NULL,
                                                val_vnp->data.ptrvalue, t->existing_text);
            }
            if (success) {
              num_fields_affected++;
              num_this_column++;
              if (t->match_mrna && IsFieldTypeCDSProduct (t->field)
                  && target_vnp->choice == OBJ_SEQFEAT) {
                if (AdjustmRNAProductToMatchProteinProduct (target_vnp->data.ptrvalue)) {
                  num_fields_affected++;
                }
              }
            } else {
              err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_col_val_fmt) + 30));
              sprintf (err_msg, bad_col_val_fmt, col_num, line_num);
              ValNodeAddPointer (&err_list, 0, err_msg);
            }
          }
        }
      }
      if (val_vnp != NULL) {
        val_vnp = val_vnp->next;
      }
      if (count_vnp == NULL) {
        ValNodeAddInt (&count_affected_list, 0, num_this_column);
      } else {
        count_vnp->data.intvalue += num_this_column;
        count_vnp = count_vnp->next;
      }
      obj_vnp = obj_vnp->next;
      col_vnp = col_vnp->next;
      col_num++;
    }
  }

  /* put message at top of list for number of fields affected */
  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_affected_fmt) + 15));
  sprintf (err_msg, num_affected_fmt, num_fields_affected);
  ValNodeAddPointer (&count_msg, 0, err_msg);

  /* if any affected, list number of fields per column, and the total in the record */
  if (num_fields_affected > 0) {
    for (count_vnp = count_affected_list, count_tot_vnp = count_list, col_vnp = columns, col_num = 1;
         count_vnp != NULL && count_tot_vnp != NULL && col_vnp != NULL;
         count_vnp = count_vnp->next, count_tot_vnp = count_tot_vnp->next, col_vnp = col_vnp->next, col_num++) {
      t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
      if (t != NULL && t->match_type == NULL) {
        qual_name = SummarizeFieldType (t->field);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (col_num_affected_fmt) + StringLen (qual_name) + 45));
        sprintf (err_msg, col_num_affected_fmt, qual_name, col_num, count_vnp->data.intvalue, count_tot_vnp->data.intvalue);      
        ValNodeAddPointer (&count_msg, 0, err_msg);
        qual_name = MemFree (qual_name);
      }
    }
  }

  ValNodeLink (&count_msg, err_list);

  count_list = ValNodeFree (count_list);
  count_affected_list = ValNodeFree (count_affected_list);

  return count_msg;
}


static int LIBCALLBACK SortVnpByChoiceAndPtrvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice > vnp2->choice) {
        return 1;
      } else if (vnp1->choice < vnp2->choice) {
        return -1;
      } else if (vnp1->data.ptrvalue > vnp2->data.ptrvalue) {
        return 1;
      } else if (vnp2->data.ptrvalue < vnp2->data.ptrvalue) {
        return -1;
      } else {
        return 0;
      }
    }
  }
  return 0;
}


static ValNodePtr FindRowsForObjectInObjectTable (ValNodePtr obj_table, Int4 column, Uint1 choice, Pointer data)
{
  Int4 col_num, row_num;
  ValNodePtr line_vnp, col_vnp, obj_vnp;
  ValNodePtr match_rows = NULL;

  if (obj_table == NULL || column < 0) {
    return NULL;
  }

  for (line_vnp = obj_table, row_num = 0; line_vnp != NULL; line_vnp = line_vnp->next, row_num++) {
    col_vnp = line_vnp->data.ptrvalue;
    col_num = 0;
    while (col_num < column && col_vnp != NULL) {
      col_vnp = col_vnp->next;
      col_num++;
    }
    if (col_vnp != NULL) {
      obj_vnp = col_vnp->data.ptrvalue;
      while (obj_vnp != NULL && (obj_vnp->choice != choice || obj_vnp->data.ptrvalue != data)) {
        obj_vnp = obj_vnp->next;
      }
      if (obj_vnp != NULL) {
        ValNodeAddInt (&match_rows, 0, row_num);
      }
    }
  }
  return match_rows;
}


static CharPtr FormatMultipleDestinationErrorMessage (Int4 col_num, ValNodePtr match_rows)
{
  CharPtr multi_fmt = "Multiple rows apply to the same object for column %d.  Matching rows:";
  CharPtr err_msg;
  Char    buf[16];
  ValNodePtr vnp;

  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (multi_fmt)
                                                + 30 + 15 * ValNodeLen (match_rows)));
  sprintf (err_msg, multi_fmt, col_num);
  for (vnp = match_rows; vnp != NULL; vnp = vnp->next) {
    sprintf (buf, "%d", vnp->data.intvalue + 1);
    StringCat (err_msg, buf);
    if (vnp->next != NULL) {
      StringCat (err_msg, ",");
    }
  }
  return err_msg;
}


NLM_EXTERN ValNodePtr CheckObjTableForRowsThatApplyToTheSameDestination (ValNodePtr obj_table)
{
  Int4 col_num;
  ValNodePtr line_vnp, col_vnp, obj_vnp, vnp;
  ValNodePtr col_list = NULL, col_obj_list;
  Boolean any_column_values_left;
  ValNodePtr err_list = NULL, match_rows;
  
  /* now, for each row, get pointer to first column */
  for (line_vnp = obj_table; line_vnp != NULL; line_vnp = line_vnp->next) {
    ValNodeAddPointer (&col_list, 0, line_vnp->data.ptrvalue);
  }

  /* now for each column, make a list of all features in the column, then sort to see if there are duplicates */
  any_column_values_left = TRUE;
  col_num = 1;
  while (any_column_values_left) {
    any_column_values_left = FALSE;
    col_obj_list = NULL;
    for (vnp = col_list; vnp != NULL; vnp = vnp->next) {
      col_vnp = vnp->data.ptrvalue;
      if (col_vnp != NULL) {
        obj_vnp = col_vnp->data.ptrvalue;
        ValNodeLink (&col_obj_list, ValNodeCopyPtr (obj_vnp));
        vnp->data.ptrvalue = col_vnp->next;
        any_column_values_left = TRUE;
      }
    }
    if (col_obj_list != NULL) {
      col_obj_list = ValNodeSort (col_obj_list, SortVnpByChoiceAndPtrvalue);
      for (vnp = col_obj_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next) {
        if (vnp->choice == vnp->next->choice
            && vnp->data.ptrvalue == vnp->next->data.ptrvalue) {
          match_rows = FindRowsForObjectInObjectTable (obj_table, col_num - 1, vnp->choice, vnp->data.ptrvalue);
          /* report rows with matches */
          ValNodeAddPointer (&err_list, col_num, FormatMultipleDestinationErrorMessage (col_num, match_rows));
          match_rows = ValNodeFree (match_rows);
          /* skip over the cluster of matches */
          while (vnp->next != NULL && vnp->choice == vnp->next->choice) {
            vnp = vnp->next;
          }
        }
      }
      col_obj_list = ValNodeFree (col_obj_list);
    }
    col_num++;
  }
  col_list = ValNodeFree (col_list);
  return err_list;
}


static CharPtr GetMatchTextForLine (ValNodePtr values, ValNodePtr columns)
{
  ValNodePtr val_vnp, col_vnp;
  CharPtr    match_txt = NULL;
  TabColumnConfigPtr t;

  for (val_vnp = values, col_vnp = columns;
       val_vnp != NULL && col_vnp != NULL;
       val_vnp = val_vnp->next, col_vnp = col_vnp->next) {
    t = col_vnp->data.ptrvalue;
    if (t != NULL && t->match_type != NULL) {
      match_txt = val_vnp->data.ptrvalue;
      break;
    }
  }
  return match_txt;
}


/* Note - when creating error messages, mark summary messages with choice = 1 */
NLM_EXTERN ValNodePtr CheckObjTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table)
{
  ValNodePtr err_list = NULL, vnp;
  ValNodePtr val_line_vnp, obj_line_vnp;
  ValNodePtr val_vnp, obj_vnp, col_vnp;
  ValNodePtr col_tot = NULL, col_tot_vnp;
  Int4       line_num = 1, col_num, num_existing_text = 0;
  Uint2      entityID;
  TabColumnConfigPtr t;
  CharPtr            err_msg, str, qual_name, val;
  CharPtr            already_has_val_fmt = "%s\t%s\t%s\t%d\t%s\t%d";
  CharPtr            num_existing_text_fmt = "%d fields already have text.\nID\tOld Value\tReplacement\tColumn\tQualifier\tLine";
  CharPtr            mrna_warn_fmt = "%d coding region features have mRNAs, but %d do not.";
  CharPtr            col_tot_fmt = "For column %d, %d out of %d fields already have text.";
  ValNodePtr         target_list, feat_vnp;
  Int4               num_with_mrna = 0, num_without_mrna = 0;
  CharPtr            match_txt;
  CharPtr            new_val;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (val_line_vnp = table, obj_line_vnp = obj_table, line_num = 1;
       val_line_vnp != NULL && obj_line_vnp != NULL;
       val_line_vnp = val_line_vnp->next, obj_line_vnp = obj_line_vnp->next, line_num++) {
    val_vnp = val_line_vnp->data.ptrvalue;
    obj_vnp = obj_line_vnp->data.ptrvalue;
    col_vnp = columns;
    if (val_vnp == NULL || obj_vnp == NULL) continue;
    col_num = 1;
    col_tot_vnp = col_tot;
    if (col_tot_vnp == NULL) {
      col_tot_vnp = ValNodeAddInt (&col_tot, 0, 0);
    }
    while (obj_vnp != NULL && col_vnp != NULL) {
      if (obj_vnp->data.ptrvalue != NULL) {
        t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
        if (t == NULL || t->match_type != NULL 
            || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
          /* ignore column or skip blank value */
        } else {
          target_list = obj_vnp->data.ptrvalue;
          if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
            val = "";
          } else {
            val = val_vnp->data.ptrvalue;
          }
          for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
            /* check for existing text */
            str = GetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
            if (!StringHasNoText (str)) {
              qual_name = SummarizeFieldType (t->field);
              match_txt = GetMatchTextForLine (val_line_vnp->data.ptrvalue, columns);
              if (match_txt == NULL) {
                match_txt = "";
              }
              new_val = StringSave (str);
              SetStringValue (&new_val, val, t->existing_text);
              err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (already_has_val_fmt)
                                                           + StringLen (match_txt)
                                                           + StringLen (str)
                                                           + StringLen (new_val)
                                                           + StringLen (qual_name)
                                                           + 30));
              sprintf (err_msg, already_has_val_fmt, match_txt, str, new_val, col_num, qual_name, line_num); 
              ValNodeAddPointer (&err_list, 0, err_msg);
              num_existing_text ++;
              new_val = MemFree (new_val);
              col_tot_vnp->data.intvalue ++;
            }
            str = MemFree (str);
            /* check for mrna if changing CDS product */
            if (IsFieldTypeCDSProduct (t->field) && feat_vnp->choice == OBJ_SEQFEAT) {
              if (GetmRNAForFeature (feat_vnp->data.ptrvalue) != NULL) {
                num_with_mrna++;
              } else {
                num_without_mrna++;
              }
            }
          }
        }
      }
      if (val_vnp != NULL) {
        val_vnp = val_vnp->next;
      }
      obj_vnp = obj_vnp->next;
      col_vnp = col_vnp->next;
      col_num++;
      col_tot_vnp = col_tot_vnp->next;
      if (col_tot_vnp == NULL) {
        col_tot_vnp = ValNodeAddInt (&col_tot, 0, 0);
      }
    }
  }          
  if (num_existing_text > 0) {
    for (col_tot_vnp = col_tot, col_num = 1; col_tot_vnp != NULL; col_tot_vnp = col_tot_vnp->next, col_num++) {
      if (col_tot_vnp->data.intvalue > 0) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (col_tot_fmt) + 45));
        sprintf (err_msg, col_tot_fmt, col_num, col_tot_vnp->data.intvalue, line_num - 1);
        ValNodeAddPointer (&err_list, 1, err_msg);
      }
    }

    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_existing_text_fmt)
                                                + 15));
    sprintf (err_msg, num_existing_text_fmt, num_existing_text);
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }
  col_tot = ValNodeFree (col_tot);
  if (num_with_mrna > 0 && num_without_mrna > 0) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (mrna_warn_fmt)
                                                + 30));
    sprintf (err_msg, mrna_warn_fmt, num_with_mrna, num_without_mrna);
    vnp = ValNodeNew (NULL);
    vnp->choice = 1;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }    

  return err_list;
}


NLM_EXTERN ValNodePtr ApplyTableToFeatures (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, val_vnp, col_vnp;
  Int4       line_num = 1, col_num;
  Uint2      entityID;
  ValNodePtr match_list, match_choice, target_list, feat_vnp, vnp;
  TabColumnConfigPtr t;
  CharPtr            err_msg;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  CharPtr            bad_col_val_fmt = "Did not set value for column %d, line %d";
  CharPtr            num_affected_fmt = "%d fields affected";
  Int4               num_fields_affected = 0;
  CharPtr            val;
  Boolean            success;
  MatchTypePtr       match_type;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  match_type = FindMatchTypeInHeader (columns);

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
    } else {
      match_list = FindMatchForRow (match_type, match_choice->data.ptrvalue, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type != NULL 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            if (val_vnp != NULL) {
              val_vnp = val_vnp->next;
            }            
            continue;
          }
          
          target_list = GetTargetListForRowAndColumn (match_type, match_list, t->field, t->constraint);
          if (target_list == NULL) {
            ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
          } else {
            if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
              val = "";
            } else {
              val = val_vnp->data.ptrvalue;
            }
            for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
              if (val[0] == 0) {
                success = RemoveFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
              } else {
                success = SetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL,
                                                  val_vnp->data.ptrvalue, t->existing_text);
              }
              if (success) {
                num_fields_affected++;
                if (t->match_mrna && IsFieldTypeCDSProduct (t->field)
                    && feat_vnp->choice == OBJ_SEQFEAT) {
                  if (AdjustmRNAProductToMatchProteinProduct (feat_vnp->data.ptrvalue)) {
                    num_fields_affected++;
                  }
                }
              } else {
                err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_col_val_fmt) + 30));
                sprintf (err_msg, bad_col_val_fmt, col_num, line_num);
                ValNodeAddPointer (&err_list, 0, err_msg);
              }
            }
          }
          target_list = ValNodeFree (target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
      match_list = ValNodeFree (match_list);
    }
  }
  
  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_affected_fmt) + 15));
  sprintf (err_msg, num_affected_fmt, num_fields_affected);
  vnp = ValNodeNew (NULL);
  vnp->data.ptrvalue = err_msg;
  vnp->next = err_list;
  err_list = vnp;
  match_type = MatchTypeFree (match_type);

  return err_list;
}


NLM_EXTERN ValNodePtr CheckTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL, vnp;
  ValNodePtr line_vnp, val_vnp, col_vnp;
  Int4       line_num = 1, col_num, num_existing_text = 0;
  Uint2      entityID;
  TabColumnConfigPtr t;
  CharPtr            err_msg, str, qual_name, val;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  CharPtr            already_has_val_fmt = "%s already has value '%s' (column %d), line %d.  Replacement is '%s'";
  CharPtr            num_existing_text_fmt = "%d fields already have text.";
  ValNodePtr         match_choice, match_list;
  ValNodePtr         target_list, feat_vnp;
  MatchTypePtr       match_type;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  match_type = FindMatchTypeInHeader (columns);
  if (match_type == NULL) return NULL;

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
    } else {
      match_list = FindMatchForRow (match_type, match_choice->data.ptrvalue, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type != NULL 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            if (val_vnp != NULL) {
              val_vnp = val_vnp->next;
            }            
            continue;
          }
          target_list = GetTargetListForRowAndColumn (match_type, match_list, t->field, t->constraint);
          if (target_list == NULL) {
            ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
          } else {
            if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
              val = "";
            } else {
              val = val_vnp->data.ptrvalue;
            }
            for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
              str = GetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
              if (!StringHasNoText (str)) {
                qual_name = SummarizeFieldType (t->field);
                err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (already_has_val_fmt)
                                                            + StringLen (qual_name) + StringLen (str)  
                                                            + StringLen (val)
                                                            + 30));
                sprintf (err_msg, already_has_val_fmt, qual_name, str, col_num, line_num, val);
                ValNodeAddPointer (&err_list, col_num, err_msg);
                num_existing_text ++;
              }
              str = MemFree (str);
            }
          }
          target_list = ValNodeFree (target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
      match_list = ValNodeFree (match_list);
    }
  }          
  if (num_existing_text > 0) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_existing_text_fmt)
                                                + 15));
    sprintf (err_msg, num_existing_text_fmt, num_existing_text);
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }

  return err_list;
}


/* Reporting functions for SMART */
static void GetDescriptorPubTitles (SeqDescrPtr sdp, Pointer userdata)
{
  CharPtr title;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || userdata == NULL) {
    return;
  }

  title = GetPubFieldFromObject (OBJ_SEQDESC, sdp, Publication_field_title, NULL);
  if (title != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, title);
  }
}


static void GetFeaturePubTitles (SeqFeatPtr sfp, Pointer userdata)
{
  CharPtr title;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB || userdata == NULL) {
    return;
  }

  title = GetPubFieldFromObject (OBJ_SEQFEAT, sfp, Publication_field_title, NULL);
  if (title != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, title);
  }
}


NLM_EXTERN ValNodePtr GetPublicationTitlesInSep (SeqEntryPtr sep)
{
  ValNodePtr title_list = NULL;

  VisitDescriptorsInSep (sep, &title_list, GetDescriptorPubTitles);
  VisitFeaturesInSep (sep, &title_list, GetFeaturePubTitles);
  return title_list;
}


NLM_EXTERN ValNodePtr GetPublicationTitlesOnSep (SeqEntryPtr sep)
{
  ValNodePtr title_list = NULL;

  VisitDescriptorsOnSep (sep, &title_list, GetDescriptorPubTitles);
  VisitFeaturesOnSep (sep, &title_list, GetFeaturePubTitles);
  return title_list;
}


static void GetBankitCommentsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  UserObjectPtr uop;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;

  if (sdp == NULL || sdp->choice != Seq_descr_user || userdata == NULL) {
    return;
  }

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
      for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
        oip = ufp->label;
        if (oip != NULL 
            && StringCmp (oip->str, "AdditionalComment") == 0
            && !StringHasNoText (ufp->data.ptrvalue)) {
          ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, StringSave (ufp->data.ptrvalue));
        }
      }
    }
  }
}


NLM_EXTERN ValNodePtr GetBankitCommentsInSep (SeqEntryPtr sep)
{
  ValNodePtr comment_list = NULL;

  VisitDescriptorsInSep (sep, &comment_list, GetBankitCommentsCallback);
  return comment_list;
}


NLM_EXTERN ValNodePtr GetBankitCommentsOnSep (SeqEntryPtr sep)
{
  ValNodePtr comment_list = NULL;

  VisitDescriptorsOnSep (sep, &comment_list, GetBankitCommentsCallback);
  return comment_list;
}




