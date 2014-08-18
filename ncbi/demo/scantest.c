/*   scantest.c
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
* File Name:  scantest.c
*
* Author:  Kans
*
* Version Creation Date:   1/20/95
*
* $Revision: 6.45 $
*
* File Description: 
*       template for custom scans of ASN.1 release files
*       (was - scans through sequence records on the Entrez discs)
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <gather.h>
#include <toasn3.h>
#include <tofasta.h>
#include <asn2gnbi.h>

typedef struct appflags {
  Boolean      binary;
  Boolean      compressed;
  Boolean      log;
  Boolean      verbose;
  Boolean      normal;
  Boolean      extended;
  Boolean      flatfile;
  FILE         *fp;
  Char         id [64];
  Boolean      is_refseq;
  SeqEntryPtr  top;
} AppFlagData, PNTR AppFlagPtr;

static ByteStorePtr Se2Bs (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}

/*
static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}
*/

typedef struct chgdata {
  Boolean     rubisco;
  Boolean     rubiscoL;
  Boolean     rubiscoS;
  Boolean     rbc;
  Boolean     rbcL;
  Boolean     rbcS;
  Boolean     its;
  Boolean     sgml;
  Boolean     rnaother;
  Boolean     trnanote;
  Boolean     oldbiomol;
  Boolean     badname;
  Boolean     hasLarge;
  Boolean     hasSmall;
  Boolean     strucSpace;
  Boolean     badDbxref;
  Boolean     refDbxref;
  Boolean     srcDbxref;
  Boolean     capDbxref;
  Boolean     oldDbxref;
  Boolean     privDbxref;
  Boolean     multDbxref;
  Boolean     badOrg;
  Int4        protdesc;
  Int4        sfpnote;
  Int4        gbsource;
  Int4        cdsconf;
  Int4        cdscodon;
  AppFlagPtr  afp;
} ChangeData, PNTR ChangeDataPtr;

static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

static Boolean HasSgml (
  CharPtr str,
  AppFlagPtr afp
)

{
  Int2  ascii_len;
  Char  buf [1024];

  if (StringHasNoText (str) || afp == NULL) return FALSE;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 > sizeof (buf)) return FALSE;

  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringCmp (str, buf) != 0) {
    if (afp->verbose) {
      fprintf (afp->fp, "GML\t%s\t%s\n", afp->id, str);
      fflush (afp->fp);
    }
    return TRUE;
  }

  return FALSE;
}

static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  CharPtr        name;
  ProtRefPtr     prp;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL || afp->fp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (HasSgml (grp->locus, afp)) {
        cdp->sgml = TRUE;
      }
      if (HasSgml (grp->desc, afp)) {
        cdp->sgml = TRUE;
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (HasSgml (str, afp)) {
          cdp->sgml = TRUE;
        }
      }
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "codon") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        (cdp->cdscodon)++;
        if (afp->verbose) {
          fprintf (afp->fp, "CDN\t%s\t%s\n", afp->id, gbq->val);
          fflush (afp->fp);
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringICmp (str, "rubisco large subunit") == 0) {
          cdp->rubiscoL = TRUE;
        }
        if (StringICmp (str, "rubisco small subunit") == 0) {
          cdp->rubiscoS = TRUE;
        }
        if (StringICmp (str, "rubisco") == 0) {
          cdp->rubisco = TRUE;
        }
        if (StringICmp (str, "RbcL") == 0) {
          cdp->rbcL = TRUE;
        }
        if (StringICmp (str, "RbcS") == 0) {
          cdp->rbcS = TRUE;
        }
        if (StringICmp (str, "Rbc") == 0) {
          cdp->rbc = TRUE;
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "ncRNA") == 0 || StringCmp (name, "tmRNA") == 0) {
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
            if (afp->verbose) {
              fprintf (afp->fp, "TR3\t%s\t%s\n", afp->id, comment);
              fflush (afp->fp);
            }
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
            if (afp->verbose) {
              fprintf (afp->fp, "TR1\t%s\t%s\n", afp->id, comment);
              fflush (afp->fp);
            }
          }
        }
      }
      break;
    default:
      break;
  }
}

static void DoKeywords (
  ChangeDataPtr cdp,
  ValNodePtr keywords
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (cdp == NULL || keywords == NULL) return;

  for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringStr (str, "rbcL gene") != NULL || StringStr (str, "large subunit") != NULL) {
      cdp->hasLarge = TRUE;
    }
    if (StringStr (str, "rbcS gene") != NULL || StringStr (str, "small subunit") != NULL) {
      cdp->hasSmall = TRUE;
    }
  }
}

static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  EMBLBlockPtr   ebp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
        DoKeywords (cdp, gbp->keywords);
      }
      break;
    case Seq_descr_embl :
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      if (ebp != NULL) {
        DoKeywords (cdp, ebp->keywords);
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
}

static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
}

static void ModGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (
  BioSourcePtr biop,
  Pointer userdata
)

{
  BoolPtr         namP;
  PCRPrimerPtr    ppp;
  PCRReactionPtr  prp;

  if (biop == NULL) return;

  ModernizePCRPrimers (biop);

  namP = (BoolPtr) userdata;
  if (namP == NULL) return;

  for (prp = biop->pcr_primers; prp != NULL; prp = prp->next) {
    if (prp->forward == NULL || prp->reverse == NULL) {
      *namP = TRUE;
      return;
    }
    for (ppp = prp->forward; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
    for (ppp = prp->reverse; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
  }
}

typedef struct bsp2cds {
  BioseqPtr   bsp;
  SeqFeatPtr  cds;
} Bsp2Cds, PNTR Bsp2CdsPtr;

static void ParentCDSProc (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Bsp2CdsPtr  bcp;
  BioseqPtr   bsp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  if (sfp->product == NULL) return;
  bcp = (Bsp2CdsPtr) userdata;
  if (bcp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp != NULL) {
    if (bsp == bcp->bsp) {
      bcp->cds = sfp;
    }
  }
}

static SeqFeatPtr FindParentCDS (
  BioseqPtr bsp,
  AppFlagPtr afp
)

{
  Bsp2Cds  b2c;

  if (bsp == NULL || afp == NULL || afp->top == NULL) return NULL;

  b2c.bsp = bsp;
  b2c.cds = NULL;

  VisitFeaturesInSep (afp->top, (Pointer) &b2c, ParentCDSProc);

  return b2c.cds;
}

static void TestForRubisco (
  SeqFeatPtr sfp,
  CharPtr str,
  ChangeDataPtr cdp
)

{
  AppFlagPtr  afp;
  BioseqPtr   bsp;
  SeqFeatPtr  cds;

  if (StringHasNoText (str)) return;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL || afp->fp == NULL) return;

  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") == 0) return;
  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") == 0) return;
  if (StringStr (str, "ribulose") == NULL || StringStr (str, "bisphosphate") == NULL) return;

  if (StringStr (str, "methyltransferase") != NULL) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "METH", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "METH", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (StringStr (str, "activate") != NULL) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "ACTV", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "ACTV", afp->id);
    }
    fflush (afp->fp);
    /*
    return;
    */
  }

  if (StringStr (str, "methyltransferase") == NULL) {
    if (StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase large chain") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase-oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5 bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase, large subunit") == 0 ||
        StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxgenase") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase, large subunit") == 0 ||
        StringICmp (str, "ribulose 5-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "ribulosebisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5 bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large chain") == 0 ||
        StringICmp (str, "large subunit ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
        StringICmp (str, "ribulose-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "ribulose-1, 5-bisphosphate carboxylase/oxygenase large-subunit") == 0) {
      if (afp->verbose) {
        fprintf (afp->fp, "%s\t%s\t%s\n", "RIBBIS", afp->id, str);
      } else {
        fprintf (afp->fp, "%s %s\n", "RIBBIS", afp->id);
      }
      fflush (afp->fp);
      return;
    }
  }

  if (StringStr (str, "large") != NULL && StringStr (str, "small") == NULL) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "RIBLRG", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "RIBLRG", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (StringStr (str, "small") != NULL && StringStr (str, "large") == NULL) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "RIBSML", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "RIBSML", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (sfp != NULL && sfp->location != NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = FindParentCDS (bsp, afp);
      if (cds != NULL && StringDoesHaveText (cds->comment)) {
        if (StringStr (cds->comment, "large") != NULL && StringStr (cds->comment, "small") == NULL) {
          if (afp->verbose) {
            fprintf (afp->fp, "%s\t%s\t%s\n", "CDSLRG", afp->id, str);
          } else {
            fprintf (afp->fp, "%s %s\n", "CDSLRG", afp->id);
          }
          fflush (afp->fp);
          return;
        }
        if (StringStr (cds->comment, "small") != NULL && StringStr (cds->comment, "large") == NULL) {
          if (afp->verbose) {
            fprintf (afp->fp, "%s\t%s\t%s\n", "CDSSML", afp->id, str);
          } else {
            fprintf (afp->fp, "%s %s\n", "CDSSML", afp->id);
          }
          fflush (afp->fp);
          return;
        }
      }
    }
  }

  if (cdp->hasSmall && cdp->hasLarge) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "RIBAMB", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "RIBAMB", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (cdp->hasLarge && (! cdp->hasSmall)) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "KEYLRG", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "KEYLRG", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (cdp->hasSmall && (! cdp->hasLarge)) {
    if (afp->verbose) {
      fprintf (afp->fp, "%s\t%s\t%s\n", "KEYSML", afp->id, str);
    } else {
      fprintf (afp->fp, "%s %s\n", "KEYSML", afp->id);
    }
    fflush (afp->fp);
    return;
  }

  if (afp->verbose) {
    fprintf (afp->fp, "%s\t%s\t%s\n", "RIBREM", afp->id, str);
  } else {
    fprintf (afp->fp, "%s %s\n", "RIBREM", afp->id);
  }
  fflush (afp->fp);
}

static void TrailingCommaFix (
  CharPtr str,
  AppFlagPtr afp,
  CharPtr prefix
)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return;
  len = StringLen (str);
  if (len < 1) return;
  ch = str [len - 1];
  while (ch == ' ' && len > 2) {
    len--;
    ch = str [len - 1];
  }
  if (ch == ',') {
    if (afp != NULL && afp->verbose && afp->fp != NULL) {
      str [len] = '\0';
      if (StringHasNoText (prefix)) {
        prefix = "?";
      }
      fprintf (afp->fp, "%s\t%s\t%s\n", prefix, afp->id, str);
      fflush (afp->fp);
    }
    str [len - 1] = '_';
    str [len] = '\0';
  }
}

static void FindCommaInGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  FILE           *fp;
  GeneRefPtr     grp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;
  fp = afp->fp;
  if (fp == NULL) return;
  if (! afp->verbose) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
  str = grp->locus;
  if (StringDoesHaveText (str)) {
    if (StringChr (str, ',') != NULL) {
      fprintf (afp->fp, "LOCCOM\t%s\t%s\n", afp->id, str);
      fflush (afp->fp);
    }
    if (StringChr (str, ';') != NULL) {
      fprintf (afp->fp, "LOCSEM\t%s\t%s\n", afp->id, str);
      fflush (afp->fp);
    }
  }
  for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringChr (str, ',') != NULL) {
      fprintf (afp->fp, "SYNCOM\t%s\t%s\n", afp->id, str);
      fflush (afp->fp);
    }
    if (StringChr (str, ';') != NULL) {
      fprintf (afp->fp, "SYNSEM\t%s\t%s\n", afp->id, str);
      fflush (afp->fp);
    }
  }
}

static void RnaProtCmntTrailingCommaFix (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  ProtRefPtr     prp;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  str = sfp->comment;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, afp, "SFPCOMM");
  }

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      TrailingCommaFix (str, afp, "PRTCOMM");
      TestForRubisco (sfp, str, cdp);
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    if (rrp->ext.choice == 1) {
      str = rrp->ext.value.ptrvalue;
      if (StringDoesHaveText (str)) {
        TrailingCommaFix (str, afp, "RNACOMM");
      }
    }
  }
}

static void ReportObsoleteDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  AppFlagPtr afp
)

{
  DbtagPtr     dp;
  ObjectIdPtr  oip;
  CharPtr      str;
  ValNodePtr   vnp;

  if (list == NULL || cdp == NULL || afp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {
      str = dp->db;
      if (StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "NID") == 0 ||
          StringICmp (str, "GI") == 0) {
        cdp->privDbxref = TRUE;
        if (afp->verbose) {
          fprintf (afp->fp, "PRVDBX\t%s\t%s\n", afp->id, str);
          fflush (afp->fp);
        }
      }
      if (StringICmp (str, "SWISS-PROT") == 0 ||
          StringICmp (str, "SWISSPROT") == 0 ||
          StringICmp (str, "SPTREMBL") == 0 ||
          StringICmp (str, "SUBTILIS") == 0 ||
          StringICmp (str, "MGD") == 0 ||
          StringCmp (str, "cdd") == 0 ||
          StringICmp (str, "TrEMBL") == 0 ||
          StringICmp (str, "LocusID") == 0 ||
          StringICmp (str, "MaizeDB") == 0 ||
          StringICmp (str, "UniProt/Swiss-Prot") == 0 ||
          StringICmp (str, "UniProt/TrEMBL") == 0 ||
          StringICmp (str, "Genew") == 0 ||
          StringICmp (str, "GENEDB") == 0 ||
          StringICmp (str, "IFO") == 0 ||
          StringICmp (str, "BHB") == 0) {
        cdp->oldDbxref = TRUE;
        if (afp->verbose) {
          fprintf (afp->fp, "OLDDBX\t%s\t%s\n", afp->id, str);
          fflush (afp->fp);
        }
      }
      if (StringICmp (str, "MGD") == 0 || StringICmp (str, "MGI") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
            cdp->oldDbxref = TRUE;
            if (afp->verbose) {
              fprintf (afp->fp, "OLDDBX\t%s\t%s\n", afp->id, str);
              fflush (afp->fp);
            }
          }
        }
      } else if (StringICmp (str, "HPRD") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "HPRD_", 5) == 0) {
            cdp->oldDbxref = TRUE;
            if (afp->verbose) {
              fprintf (afp->fp, "OLDDBX\t%s\t%s\n", afp->id, str);
              fflush (afp->fp);
            }
          }
        }
      }
      oip = dp->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringChr (oip->str, ':') != NULL) {
          cdp->multDbxref = TRUE;
          if (afp->verbose) {
            fprintf (afp->fp, "MLTDBX\t%s\t%s\t%s\n", afp->id, dp->db, oip->str);
            fflush (afp->fp);
          }
        }
      }
    }
  }
}

static void ReportInvalidDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  AppFlagPtr afp,
  Boolean is_source
)

{
  Boolean     cap;
  DbtagPtr    dp;
  CharPtr     good;
  Boolean     ref;
  Boolean     src;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL || cdp == NULL || afp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {
      str = dp->db;

      if (is_source && StringCmp (str, "taxon") == 0) continue;

      if (DbxrefIsValid (str, &ref, &src, &cap, &good)) {
        if (ref && (! afp->is_refseq)) {
          cdp->refDbxref = TRUE;
          if (afp->verbose) {
            fprintf (afp->fp, "REFDBX\t%s\t%s\n", afp->id, str);
            fflush (afp->fp);
          }
        }
        if (is_source && (! src)) {
          cdp->srcDbxref = TRUE;
          if (afp->verbose) {
            fprintf (afp->fp, "SRCDBX\t%s\t%s\n", afp->id, str);
            fflush (afp->fp);
          }
        }
        if (cap) {
          cdp->capDbxref = TRUE;
          if (afp->verbose) {
            fprintf (afp->fp, "CAPDBX\t%s\t%s\n", afp->id, str);
            fflush (afp->fp);
          }
        }
      } else {
        cdp->badDbxref = TRUE;
        if (afp->verbose) {
          fprintf (afp->fp, "BADDBX\t%s\t%s\n", afp->id, str);
          fflush (afp->fp);
        }
      }
    }
  }
}

static void LookForObsoleteFeatDbxref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  GeneRefPtr     grp;
  OrgRefPtr      orp = NULL;
  ProtRefPtr     prp;

  if (sfp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        ReportObsoleteDbxref (grp->db, cdp, afp);
        ReportInvalidDbxref (grp->db, cdp, afp, FALSE);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp != NULL) {
        ReportObsoleteDbxref (prp->db, cdp, afp);
        ReportInvalidDbxref (prp->db, cdp, afp, FALSE);
      }
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      cdp->badOrg = TRUE;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
      }
    default :
      break;
  }

  if (orp != NULL) {
    ReportObsoleteDbxref (orp->db, cdp, afp);
    ReportInvalidDbxref (orp->db, cdp, afp, TRUE);
  }

  ReportObsoleteDbxref (sfp->dbxref, cdp, afp);
  ReportInvalidDbxref (sfp->dbxref, cdp, afp, FALSE);
}

static void LookForObsoleteDescDbxref (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  OrgRefPtr      orp = NULL;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      cdp->badOrg = TRUE;
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
      }
      break;
    default :
      break;
  }

  if (orp != NULL) {
    ReportObsoleteDbxref (orp->db, cdp, afp);
    ReportInvalidDbxref (orp->db, cdp, afp, TRUE);
  }
}

static void LookForSemicolonedStrains (
  BioSourcePtr biop,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  OrgModPtr      omp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  CharPtr        str;

  if (biop == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;

  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_strain) continue;
    str = omp->subname;
    if (StringHasNoText (str)) continue;
    if (StringChr (str, ';') == NULL) continue;
    if (afp->verbose) {
      fprintf (afp->fp, "STR\t%s\t%s\n", afp->id, str);
    } else {
      fprintf (afp->fp, "STR %s\n", afp->id);
    }
    fflush (afp->fp);
  }
}

static void LookForBadAuth (
  NameStdPtr nsp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  Char           ch;
  Int2           i;
  Boolean        is_bad = FALSE;
  CharPtr        prefix = "\t";
  CharPtr        str;

  if (nsp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  for (i = 0; i < 6; i++) {
    str = nsp->names [i];
    if (StringHasNoText (str)) continue;
    ch = *str;
    while (ch != '\0') {
      if (IS_DIGIT (ch)) {
        cdp->badname = TRUE;
        is_bad = TRUE;
      }
      str++;
      ch = *str;
    }
  }

  if (is_bad && afp->fp != NULL && afp->verbose) {
    fprintf (afp->fp, "%s\t%s", "AUTHOR", afp->id);
    for (i = 0; i < 6; i++) {
      str = nsp->names [i];
      if (StringHasNoText (str)) continue;
      fprintf (afp->fp, "%s%s", prefix, str);
      prefix = " | ";
    }
    fprintf (afp->fp, "\n");
    fflush (afp->fp);
  }
}

static void LookForBadPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  VisitAuthorsInPub (pdp, userdata, LookForBadAuth);
}

static void LookForStrucComment (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  ObjectIdPtr    oip;
  CharPtr        str;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "StructuredComment") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice != 1) continue;
    oip = ufp->label;
    if (oip == NULL) continue;
    str = oip->str;
    if (StringHasNoText (str)) continue;
    if (StringStr (str, " ") != NULL) {
      cdp->strucSpace = TRUE;
      if (afp->fp != NULL && afp->verbose) {
        fprintf (afp->fp, "STCCOM\t%s\t%s\n", afp->id, str);
        fflush (afp->fp);
      }
    }
  }
}

static void CommentDescrTrailingCommaFix (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  CharPtr        str;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, afp, "DSCCOMM");
  }
}

static void StripBadProtTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr            buf;
  size_t             buflen = 1001;
  ObjValNodePtr      ovp;
  SeqIdPtr           sip;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  title = (CharPtr) vnp->data.ptrvalue;
  if (StringHasNoText (title)) return;

  buf = MemNew (sizeof (Char) * (buflen + 1));
  if (buf == NULL) return;

  if (NewCreateDefLineBuf (NULL, bsp, buf, buflen, TRUE, FALSE)) {
    if (StringICmp (buf, title) != 0) {
      if (vnp->extended != 0) {
        ovp = (ObjValNodePtr) vnp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }

  MemFree (buf);
}

static void BadProtTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSep (sep, NULL, StripBadProtTitles);
}

static void DoReport (
  SeqEntryPtr sep,
  AppFlagPtr afp
)

{
  ByteStorePtr  bs = NULL, tmp = NULL;
  Boolean       bsec = FALSE, cma = FALSE, norm = FALSE, ttl = FALSE, ssec = FALSE;
  Boolean       gen = FALSE, ncr = FALSE, pcr = FALSE, nam = FALSE;
  ChangeData    cdbefore, cdafter;
  Uint2         entityID;

  if (sep == NULL || afp == NULL) return;

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  cdbefore.afp = afp;
  cdafter.afp = afp;

  CheckForChanges (sep, &cdbefore);

  if (afp->extended) {
    bs = Se2Bs (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      norm = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  VisitFeaturesInSep (sep, (Pointer) &cdbefore, RnaProtCmntTrailingCommaFix);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, CommentDescrTrailingCommaFix);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, LookForStrucComment);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, LookForObsoleteFeatDbxref);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, LookForObsoleteDescDbxref);
  VisitBioSourcesInSep (sep, (Pointer) &cdbefore, LookForSemicolonedStrains);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, FindCommaInGene);

  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    cma = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  BasicSeqEntryCleanup (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    bsec = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, (Pointer) &cdbefore, LookForBadPub);

  if (afp->extended) {
    VisitFeaturesInSep (sep, NULL, ModGenes);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      gen = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    VisitFeaturesInSep (sep, NULL, ModRNAs);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ncr = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    VisitBioSourcesInSep (sep, (Pointer) &nam, ModPCRs);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      pcr = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    entityID = ObjMgrGetEntityIDForChoice (sep);
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadProtTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ttl = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    SeriousSeqEntryCleanup (sep, NULL, NULL);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  CheckForChanges (sep, &cdafter);

  BSFree (bs);

  if (afp->extended) {
    if (ssec) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "SSEC %s\n", afp->id);
        fflush (afp->fp);
      }
    } else if (ttl) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "TITL %s\n", afp->id);
        fflush (afp->fp);
      }
    } else if (bsec) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "BSEC %s\n", afp->id);
        fflush (afp->fp);
      }
    } else if (norm) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "NORM %s\n", afp->id);
        fflush (afp->fp);
      }
    } else {
      /*
      if (afp->fp != NULL) {
        fprintf (afp->fp, "OKAY %s\n", afp->id);
        fflush (afp->fp);
      }
      */
    }
  }

  if (cma) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "CMA %s\n", afp->id);
      fflush (afp->fp);
    }
  }

  if (afp->extended) {
    if (gen) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "GEN %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (ncr) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "NCR %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (pcr) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "PCR %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (nam) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "NAM %s\n", afp->id);
        fflush (afp->fp);
      }
    }
  }

  if (cdbefore.rubiscoL) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RUL %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rubiscoS) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RUS %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rubisco) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RUB %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rbcL) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RBL %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rbcS) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RBS %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rbc) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RBC %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.its) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "ITS %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.sgml) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SGM %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rnaother) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RNA %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.trnanote) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "TRN %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.oldbiomol) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "MOL %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.badname) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "AUT %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.strucSpace) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SPC %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (! afp->verbose) {
    if (cdbefore.badDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "BDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.refDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "RDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.srcDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "SDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.capDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "CDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.privDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "PDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.oldDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "ODBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.multDbxref) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "MDBX %s\n", afp->id);
        fflush (afp->fp);
      }
    }
    if (cdbefore.cdscodon > 0) {
      if (afp->fp != NULL) {
        fprintf (afp->fp, "CDN %s\n", afp->id);
        fflush (afp->fp);
      }
    }
  }
  if (cdbefore.badOrg) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "ORG %s\n", afp->id);
      fflush (afp->fp);
    }
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "PRT %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "COM %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SRC %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "CNF %s\n", afp->id);
      fflush (afp->fp);
    }
  }
}

static void DoFlatfile (
  SeqEntryPtr sep,
  AppFlagPtr afp
)

{
#ifdef OS_UNIX
  Char    buf [256];
  Char    cmmd [256];
  size_t  ct;
  FILE    *fp;
  Char    path1 [PATH_MAX];
  Char    path2 [PATH_MAX];
  Char    path3 [PATH_MAX];

  if (sep == NULL || afp == NULL || afp->fp == NULL) return;

  fprintf (afp->fp, "%s\n", afp->id);
  fflush (afp->fp);

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  fp = FileOpen (path1, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, NORMAL_STYLE, 1, 0, 0, NULL, fp);
  }
  FileClose (fp);
  /*
  SeriousSeqEntryCleanupBulk (sep);
  */
  fp = FileOpen (path2, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, NORMAL_STYLE, 262145, 0, 0, NULL, fp);
  }
  FileClose (fp);

  sprintf (cmmd, "diff %s %s > %s", path1, path2, path3);
  system (cmmd);

  sprintf (cmmd, "cat %s", path3);
  fp = popen (cmmd, "r");
  if (fp != NULL) {
    while ((ct = fread (buf, 1, sizeof (buf), fp)) > 0) {
      fwrite (buf, 1, ct, afp->fp);
      fflush (afp->fp);
    }
    pclose (fp);
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);
#endif
}

static void DoRecord (SeqEntryPtr sep, Pointer userdata)

{
  AppFlagPtr   afp;
  BioseqPtr    fbsp;
  SeqEntryPtr  fsep;
  SeqIdPtr     sip, siphead;

  if (sep == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep == NULL) return;
  fbsp = (BioseqPtr) fsep->data.ptrvalue;
  if (fbsp == NULL) return;

  afp->is_refseq = FALSE;
  siphead = SeqIdSetDup (fbsp->id);
  for (sip = siphead; sip != NULL; sip = sip->next) {
    SeqIdStripLocus (sip);
    if (sip->choice == SEQID_OTHER) {
      afp->is_refseq = TRUE;
    }
  }
  SeqIdWrite (siphead, afp->id, PRINTID_FASTA_LONG, sizeof (afp->id));
  SeqIdSetFree (siphead);

  if (afp->log) {
    fprintf (afp->fp, "LOG %s\n", afp->id);
    fflush (afp->fp);
  }

  afp->top = sep;

  if (afp->flatfile) {
    DoFlatfile (sep, afp);
    return;
  }

  if (afp->normal || afp->extended) {
    DoReport (sep, afp);
  }
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (StringHasNoText (filename)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  if (StringStr (filename, "gbcon") != NULL ||
      StringStr (filename, "gbest") != NULL ||
      StringStr (filename, "gbgss") != NULL ||
      StringStr (filename, "gbhtg") != NULL ||
      StringStr (filename, "gbsts") != NULL) {
    printf ("Skipping %s\n", filename);
    return;
  }

  printf ("%s\n", filename);
  fflush (stdout);

  fprintf (afp->fp, "%s\n", filename);
  fflush (afp->fp);

  ScanBioseqSetRelease (filename, afp->binary, afp->compressed, (Pointer) afp, DoRecord);

  fprintf (afp->fp, "\n");
  fflush (afp->fp);
}

#define p_argInputPath    0
#define i_argInputFile    1
#define o_argOutputFile   2
#define f_argFilter       3
#define x_argSuffix       4
#define u_argRecurse      5
#define b_argBinary       6
#define c_argCompressed   7
#define l_argLog          8
#define v_argVerbose      9
#define n_argNormal      10
#define e_argExtended    11
#define q_argFlatfile    12

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Input File Name", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log SeqID", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbose Output", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Normal Tests", "F", NULL, NULL,
    TRUE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Extended Tests", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Flatfile Report", "F", NULL, NULL,
    TRUE, 'q', ARG_BOOLEAN, 0.0, 0, NULL},
};

extern Int2 Main (void)

{
  AppFlagData  afd;
  Boolean      dorecurse;
  CharPtr      filter, infile, outfile, directory, suffix;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  if (! GetArgs ("scantest", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &afd, 0, sizeof (AppFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  afd.binary = (Boolean) myargs [b_argBinary].intvalue;
  afd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  afd.log = (Boolean) myargs [l_argLog].intvalue;
  afd.verbose = (Boolean) myargs [v_argVerbose].intvalue;
  afd.normal = (Boolean) myargs [n_argNormal].intvalue;
  afd.extended = (Boolean) myargs [e_argExtended].intvalue;
  afd.flatfile = (Boolean) myargs [q_argFlatfile].intvalue;

  if (afd.flatfile) {
    if (afd.normal || afd.extended) {
      Message (MSG_FATAL, "-q cannot be used with -n or -e");
      return 1;
    }
  }

  afd.fp = FileOpen (outfile, "w");
  if (afd.fp == NULL) {
    return 0;
  }

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, dorecurse, ProcessOneRecord, (Pointer) &afd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, &afd);
  }

  FileClose (afd.fp);

  return 0;
}

