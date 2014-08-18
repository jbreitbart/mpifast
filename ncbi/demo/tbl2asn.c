/*   tbl2asn.c
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
* File Name:  tbl2asn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   5/5/00
*
* $Revision: 6.278 $
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

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <edutil.h>
#include <seqport.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <toasn3.h>
#include <valid.h>
#include <asn2gnbk.h>
#include <explore.h>
#include <tofasta.h>
#include <simple.h>
#include <suggslp.h>
#include <aliparse.h>
#include <util/creaders/alnread.h>
#include <pmfapi.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_TBL2ASN
#include <accpubseq.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <asnmacro.h>
#include <objmacro.h>
#include <macroapi.h>

#define TBL2ASN_APP_VER "13.2"

CharPtr TBL2ASN_APPLICATION = TBL2ASN_APP_VER;

typedef struct cleanupargs {
  Boolean collection_dates;
  Boolean collection_dates_month_first;
  Boolean add_notes_to_overlapping_cds_without_abc;
} CleanupArgsData, PNTR CleanupArgsPtr;

typedef struct tblargs {
  Boolean     raw2delt;
  Int2        r2dmin;
  Boolean     r2dunk100;
  Boolean     fastaset;
  Int2        whichclass;
  Boolean     deltaset;
  Boolean     alignset;
  Boolean     gapped;
  Boolean     phrapace;
  Boolean     genprodset;
  Boolean     linkbyoverlap;
  Boolean     linkbyproduct;
  Boolean     implicitgaps;
  Boolean     forcelocalid;
  Boolean     gpstonps;
  Boolean     gnltonote;
  Boolean     removeunnecxref;
  Boolean     dotaxlookup;
  Boolean     dopublookup;
  CharPtr     accn;
  CharPtr     center;
  CharPtr     organism;
  CharPtr     srcquals;
  CharPtr     comment;
  CharPtr     commentFile;
  CharPtr     tableFile;
  Boolean     findorf;
  Boolean     runonorf;
  Boolean     altstart;
  Boolean     conflict;
  Boolean     validate;
  Boolean     relaxed;
  Boolean     validate_barcode;
  Boolean     flatfile;
  Boolean     genereport;
  Boolean     seqidfromfile;
  Boolean     smartfeats;
  Boolean     smarttitle;
  Boolean     logtoterminal;
  CharPtr     aln_beginning_gap;
  CharPtr     aln_end_gap;
  CharPtr     aln_middle_gap;
  CharPtr     aln_missing;
  CharPtr     aln_match;
  Boolean     aln_is_protein;
  Boolean     save_bioseq_set;
  Boolean     apply_cmt_to_all;

  GlobalDiscrepReportPtr global_report;

  CleanupArgsData cleanup_args;
} TblArgs, PNTR TblArgsPtr;

static FILE* OpenOneFile (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix
)

{
  Char  file [FILENAME_MAX], path [PATH_MAX];

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  return FileOpen (path, "r");
}

static void WriteOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  CharPtr outfile,
  SeqEntryPtr sep,
  SubmitBlockPtr sbp,
  Boolean save_bioseq_set
)

{
  AsnIoPtr      aip;
  BioseqSetPtr  bssp;
  Char          file [FILENAME_MAX], path [PATH_MAX];
  SeqSubmit     ssb;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;

  MemSet ((Pointer) &ssb, 0, sizeof (SeqSubmit));
  ssb.sub = sbp;
  ssb.datatype = 1;
  ssb.data = (Pointer) sep;

  if (StringDoesHaveText (outfile)) {
    StringNCpy_0 (path, outfile, sizeof (path));
  } else {
    StringNCpy_0 (path, results, sizeof (path));
    sprintf (file, "%s%s", base, suffix);
    FileBuildPath (path, NULL, file);
  }

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  if (sbp != NULL) {
    SeqSubmitAsnWrite (&ssb, aip, NULL);
  } else if (save_bioseq_set && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    BioseqSetAsnWrite (bssp, aip, NULL);
  } else {
    SeqEntryAsnWrite (sep, aip, NULL);
  }

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static CharPtr compatSeverityLabel [] = {
  "NONE", "NOTE: valid", "WARNING: valid", "ERROR: valid", "REJECT: valid", "FATAL: valid", "MAX", NULL
};

static void LIBCALLBACK ValidCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  CharPtr  catname, errname;
  FILE     *fp;

  fp = (FILE *) userdata;
  if (fp == NULL) return;

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  catname = GetValidCategoryName (errcode);
  errname = GetValidErrorName (errcode, subcode);

  if (catname == NULL) {
    catname = "?";
  }
  if (errname == NULL) {
    errname = "?";
  }

  if (accession == NULL) {
    accession = "";
  }
  if (message == NULL) {
    message = "";
  }
  if (objtype == NULL) {
    objtype = "";
  }
  if (label == NULL) {
    label = "";
  }

  fprintf (fp, "%s [%s.%s] %s %s: %s",
           compatSeverityLabel [severity],
           catname, errname, message, objtype, label);
  if (location != NULL) {
    fprintf (fp, " %s", location);
  }
  if (context != NULL) {
    fprintf (fp, " %s", context);
  }
  if (product != NULL) {
    fprintf (fp, " -> %s", product);
  }
  fprintf (fp, "\n");
}


static void ValidateOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep,
  Boolean standard,
  Boolean relaxed,
  Boolean barcode
)

{
  Char            file [FILENAME_MAX], path [PATH_MAX];
  FILE            *ofp;
  ErrSev          oldErrSev;
  ValidStructPtr  vsp;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  ofp = FileOpen (path, "w");

  if (standard) {
    vsp = ValidStructNew ();
    if (vsp != NULL) {
      vsp->useSeqMgrIndexes = TRUE;
      vsp->suppressContext = TRUE;
      vsp->seqSubmitParent = TRUE;
      if (! relaxed) {
        vsp->testLatLonSubregion = TRUE;
      }
      oldErrSev = ErrSetMessageLevel (SEV_NONE);
      vsp->errfunc = ValidCallback;
      vsp->userdata = (Pointer) ofp;
      /* vsp->convertGiToAccn = FALSE; */
      ValidateSeqEntry (sep, vsp);
      ValidStructFree (vsp);
      ErrSetMessageLevel (oldErrSev);
    }
  }
  /* Barcode results if requested */
  if (barcode) {
    BarcodeValidateOneSeqEntry (ofp, sep, TRUE, FALSE, TRUE, NULL);
  }

  FileClose (ofp);
}

static void FlatfileOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char    file [FILENAME_MAX], path [PATH_MAX];
  FILE    *fp;
  ErrSev  oldErrSev;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  fp = FileOpen (path, "w");
  if (fp == NULL) return;

  oldErrSev = ErrSetMessageLevel (SEV_MAX);
  SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  ErrSetMessageLevel (oldErrSev);

  FileClose (fp);
}

/* for full-length cDNAs, allow automatic annotation of largest internal ORF */

typedef struct orfdata {
  Int4     curlen [6], bestlen [6], currstart [6], beststart [6], sublen [6];
  Boolean  inorf [6], altstart, runonorf;
  Int4     bioseq_len;
} OrfData, PNTR OrfDataPtr;

static Boolean TreatLikeStop (Int2 frame, Int4 pos, Uint1 strand, Int4 len)
{
  Int4 remainder = len % 3;
  Boolean like_stop = FALSE;

  if (strand == Seq_strand_minus) {
    if (pos < 3) {
      like_stop = TRUE;
    }
  } else {
    if (pos >= len - remainder - 3) {
      like_stop = TRUE;
    }
  }
  return like_stop;
}

static void LIBCALLBACK LookForOrfs (
  Int4 position,
  Char residue,
  Boolean atgStart,
  Boolean altStart,
  Boolean orfStop,
  Int2 frame,
  Uint1 strand,
  Pointer userdata
)

{
  Int2        idx;
  OrfDataPtr  odp;
  Boolean     start_of_seq = FALSE;

  odp = (OrfDataPtr) userdata;
  if (strand == Seq_strand_plus) {

    /* top strand */

    idx = frame;
    if (odp->inorf [idx]) {
      if (!orfStop && odp->runonorf) {
        /* treat the end of the sequence like a stop codon */
        if (TreatLikeStop(frame, position, strand, odp->bioseq_len)) {
          (odp->curlen[idx])++;
          orfStop = TRUE;
        }
      }

      if (orfStop) {
        odp->inorf [idx] = FALSE;
        if (odp->curlen [idx] > odp->bestlen [idx]) {
          odp->bestlen [idx] = odp->curlen [idx];
          odp->beststart [idx] = odp->currstart [idx];
        }
      } else {
        (odp->curlen [idx])++;
      }
    } else if (atgStart || (altStart && odp->altstart)) {
      odp->inorf [idx] = TRUE;
      odp->curlen [idx] = 1;
      odp->currstart [idx] = position - frame;
    }
  } else {

    /* bottom strand */

    idx = frame + 3;

    if (!orfStop && odp->runonorf) {
      start_of_seq = TreatLikeStop (frame, position, strand, odp->bioseq_len);
    }

    if (orfStop) {
      odp->curlen [idx] = 0;
      odp->sublen [idx] = 0;
      odp->currstart [idx] = position - frame;
    } else if (start_of_seq) {
      odp->curlen [idx] = 1;
      odp->sublen [idx] = 1;
      odp->currstart [idx] = position - frame - 3;
      if (odp->curlen [idx] > odp->bestlen [idx]) {
        odp->bestlen [idx] = odp->curlen [idx];
        odp->beststart [idx] = odp->currstart [idx];
      }
    } else if (atgStart || (altStart && odp->altstart)) {
      (odp->sublen [idx])++;
      odp->curlen [idx] = odp->sublen [idx];
      if (odp->curlen [idx] > odp->bestlen [idx]) {
        odp->bestlen [idx] = odp->curlen [idx];
        odp->beststart [idx] = odp->currstart [idx];
      }
    } else {
      (odp->sublen [idx])++;
    }
  }
}

static SeqFeatPtr AnnotateBestOrf (
  BioseqPtr bsp,
  Int2 genCode,
  Boolean altstart,
  Boolean runonorf,
  SqnTagPtr stp
)

{
  SeqFeatPtr      cds = NULL;
  CdRegionPtr     crp;
  GeneRefPtr      grp;
  Int2            i, best, idx;
  OrfData         od;
  ProtRefPtr      prp;
  SeqFeatPtr      sfp;
  SeqInt          sint;
  CharPtr         str;
  TransTablePtr   ttp;
  ValNode         vn;
  SeqFeatXrefPtr  xref;
  Boolean         partial5 = FALSE, partial3 = FALSE;

  if (bsp == NULL) return NULL;
  for (i = 0; i < 6; i++) {
    od.curlen [i] = INT4_MIN;
    od.bestlen [i] = 0;
    od.currstart [i] = 0;
    od.beststart [i] = 0;
    od.sublen [i] = INT4_MIN;
    od.inorf [i] = FALSE;
  }
  od.altstart = altstart;
  od.runonorf = runonorf;
  od.bioseq_len = bsp->length;

  /* use simultaneous 6-frame translation finite state machine */

  ttp = PersistentTransTableByGenCode (genCode);
  if (ttp != NULL) {
    TransTableProcessBioseq (ttp, LookForOrfs, (Pointer) &od, bsp);
  }
  /* TransTableFree (tbl); - now using persistent tables, free at end */
  best = -1;
  idx = -1;
  for (i = 0; i < 6; i++) {
    if (od.bestlen [i] > best) {
      best = od.bestlen [i];
      idx = i;
    }
  }
  if (idx == -1) return NULL;

  /* make feature location on largest ORF */

  if (idx < 3) {
    MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
    sint.from = od.beststart [idx] + idx;
    sint.to = sint.from + (od.bestlen [idx]) * 3 + 2;
    if (sint.to > od.bioseq_len - 1) {
      sint.to = od.bioseq_len - 1;
      partial3 = TRUE;
    }
    sint.id = SeqIdFindBest (bsp->id, 0);
    sint.strand = Seq_strand_plus;
    vn.choice = SEQLOC_INT;
    vn.extended = 0;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
  } else {
    MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
    sint.from = od.beststart [idx] + idx - 3;
    sint.to = sint.from + (od.bestlen [idx]) * 3 + 2;
    if (sint.from < 0) {
      sint.from = 0;
      partial3 = TRUE;
    }
    sint.id = SeqIdFindBest (bsp->id, 0);
    sint.strand = Seq_strand_minus;
    vn.choice = SEQLOC_INT;
    vn.extended = 0;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
  }

  SetSeqLocPartial (&vn, partial5, partial3);

  /* make CDS feature with unknown product - now check [protein=...] */

  cds = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, &vn);
  if (cds == NULL) return NULL;
  if (partial5 || partial3) {
    cds->partial = TRUE;
  }
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (crp == NULL) return NULL;
  crp->frame = 1;
  cds->data.value.ptrvalue = (Pointer) crp;

  prp = ProtRefNew ();
  if (prp == NULL) return cds;
  xref = SeqFeatXrefNew ();
  if (xref == NULL) return cds;
  xref->data.choice = SEQFEAT_PROT;
  xref->data.value.ptrvalue = (Pointer) prp;
  xref->next = cds->xref;
  cds->xref = xref;
  prp = ParseTitleIntoProtRef (stp, prp);
  if (prp->name == NULL && prp->desc == NULL) {
    prp->name = ValNodeCopyStr (NULL, 0, "unknown");
  }

  /* parse CDS comment ("note" goes to biosource) and experimental evidence */

  str = SqnTagFind (stp, "comment");
  if (StringDoesHaveText (str)) {
    cds->comment = StringSave (str);
  }

  str = SqnTagFind (stp, "evidence");
  if (StringICmp (str, "experimental") == 0) {
    cds->exp_ev = 1;
  }

  /* now check [gene=...], make gene feature if locus or synonym present */

  grp = GeneRefNew ();
  if (grp == NULL) return cds;
  grp = ParseTitleIntoGeneRef (stp, grp);
  if (grp->locus == NULL && grp->syn == NULL) {
    GeneRefFree (grp);
    return cds;
  }
  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
  if (sfp == NULL) return cds;
  sfp->data.value.ptrvalue = (Pointer) grp;

  return cds;
}

/* change all feature IDs to entered accession */

static void PromoteSeqId (
  SeqIdPtr sip,
  Pointer userdata
)

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

static void CorrectFeatureSeqIds (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqLoc (sfp->location, userdata, PromoteSeqId);
}

static void CorrectGraphSeqIds (
  SeqGraphPtr sgp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqGraph (sgp, userdata, PromoteSeqId);
}

/* source information for several common organisms sequenced by genome centers */

typedef struct orgstuff {
  CharPtr  taxname;
  CharPtr  common;
  CharPtr  lineage;
  CharPtr  division;
  Uint1    gcode;
  Uint1    mgcode;
  Int4     taxID;
} OrgStuff, PNTR OrfStuffPtr;

static OrgStuff commonOrgStuff [] = {
  {
    "Saccharomyces cerevisiae", "baker's yeast",
    "Eukaryota; Fungi; Ascomycota; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces",
    "PLN", 1, 3, 4932
  },
  {
    "Drosophila melanogaster", "fruit fly",
    "Eukaryota; Metazoa; Arthropoda; Tracheata; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Ephydroidea; Drosophilidae; Drosophila",
    "INV", 1, 5, 7227
  },
  {
    "Homo sapiens", "human",
    "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo",
    "PRI", 1, 2, 9606
  },
  {
    "Escherichia coli", "",
    "Bacteria; Proteobacteria; gamma subdivision; Enterobacteriaceae; Escherichia",
    "BCT", 11, 0, 562
  },
  {
    "Helicobacter pylori", "",
    "Bacteria; Proteobacteria; epsilon subdivision; Helicobacter group; Helicobacter",
    "BCT", 11, 0, 210
  },
  {
    "Arabidopsis thaliana", "thale cress",
    "Eukaryota; Viridiplantae; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicots; Rosidae; eurosids II; Brassicales; Brassicaceae; Arabidopsis",
    "PLN", 1, 1, 3702
  },
  {
    "Mus musculus", "house mouse",
    "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Rodentia; Sciurognathi; Muridae; Murinae; Mus",
    "ROD", 1, 2, 10090
  },
  {
    "Rattus norvegicus", "Norway rat",
    "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Rodentia; Sciurognathi; Muridae; Murinae; Rattus",
    "ROD", 1, 2, 10116
  },
  {
    "Danio rerio", "zebrafish",
    "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii; Teleostei; Euteleostei; Ostariophysi; Cypriniformes; Cyprinidae; Rasborinae; Danio",
    "VRT", 1, 2, 7955
  },
  {
    "Zea mays", "",
    "Eukaryota; Viridiplantae; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Poales; Poaceae; Zea",
    "PLN", 1, 1, 4577
  },
  {
    "Caenorhabditis elegans", "",
    "Eukaryota; Metazoa; Nematoda; Chromadorea; Rhabditida; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis",
    "INV", 1, 5, 6239
  },
  {
    "Caenorhabditis briggsae", "",
    "Eukaryota; Metazoa; Nematoda; Chromadorea; Rhabditida; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis",
    "INV", 1, 5, 6238
  },
  {
    "Anopheles gambiae", "African malaria mosquito",
    "Eukaryota; Metazoa; Arthropoda; Tracheata; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Nematocera; Culicoidea; Anopheles",
    "INV", 1, 5, 7165
  },
  {
    "Anopheles gambiae str. PEST", "African malaria mosquito",
    "Eukaryota; Metazoa; Arthropoda; Tracheata; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Nematocera; Culicoidea; Anopheles",
    "INV", 1, 5, 180454
  },
  {
    "Tetrahymena thermophila", "",
    "Eukaryota; Alveolata; Ciliophora; Oligohymenophorea; Hymenostomatida; Tetrahymenina; Tetrahymena",
    "INV", 6, 4, 5911
  },
  {
    "Pan troglodytes", "chimpanzee",
    "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Pan",
    "PRI", 1, 2, 9598
  },
  {
    "Candida albicans", "",
    "Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; mitosporic Saccharomycetales; Candida",
    "PLN", 12, 4, 5476
  },
  {
    "Candida albicans SC5314", "",
    "Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; mitosporic Saccharomycetales; Candida",
    "PLN", 12, 4, 237561
  },
  {
    "Trypanosoma brucei", "",
    "Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Trypanosoma",
    "INV", 1, 4, 5691
  },
  {
    "Trypanosoma cruzi", "",
    "Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Trypanosoma; Schizotrypanum",
    "INV", 1, 4, 5693
  },
  {
    "Oryza sativa", "",
    "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Poales; Poaceae; Ehrhartoideae; Oryzeae; Oryza",
    "PLN", 1, 1, 4530
  },
  {
    "Oryza sativa (indica cultivar-group)", "",
    "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Poales; Poaceae; Ehrhartoideae; Oryzeae; Oryza",
    "PLN", 1, 1, 39946
  },
  {
    "Oryza sativa (japonica cultivar-group)", "",
    "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; Poales; Poaceae; Ehrhartoideae; Oryzeae; Oryza",
    "PLN", 1, 1, 39947
  },
  {
    "Aspergillus nidulans FGSC A4", "",
    "Eukaryota; Fungi; Ascomycota; Pezizomycotina; Eurotiomycetes; Eurotiales; Trichocomaceae; Emericella",
    "PLN", 1, 4, 227321
  },
  {
    "environmental sequence", "",
    "unclassified; environmental samples",
    "UNA", 1, 2, 256318
  },
  {
    NULL, NULL, NULL, NULL, 0, 0, 0
  }
};

static Boolean HasTaxon (
  OrgRefPtr orp
)

{
  ValNodePtr  db;
  DbtagPtr    dbt;

  if (orp == FALSE) return FALSE;
  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL && dbt->db != NULL &&
        StringICmp (dbt->db, "taxon") == 0) return TRUE;
  }
  return FALSE;
}

static void AddMissingSourceInfo (
  BioSourcePtr biop
)

{
  ValNodePtr   db;
  DbtagPtr     dbt;
  Int2         idx;
  ObjectIdPtr  oip;
  OrgNamePtr   onp;
  OrgRefPtr    orp;
  OrfStuffPtr  osp;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;

  /* look for entry of organisms in commonOrgStuff table */

  for (idx = 0; commonOrgStuff [idx].taxname != NULL; idx++) {
    osp = &(commonOrgStuff [idx]);
    if (StringICmp (orp->taxname, osp->taxname) == 0) {
      if (StringCmp (orp->taxname, osp->taxname) != 0) {
        /* fix capitalization of supplied name if in common organism list */
        StringCpy (orp->taxname, osp->taxname);
      }
      if (StringHasNoText (orp->common) && StringDoesHaveText (osp->common)) {
        orp->common = StringSave (osp->common);
      }
      if (onp->gcode == 0) {
        onp->gcode = osp->gcode;
      }
      if (onp->mgcode == 0) {
        onp->mgcode = osp->mgcode;
      }
      if (StringHasNoText (onp->div)) {
        onp->div = StringSave (osp->division);
      }
      if (StringHasNoText (onp->lineage)) {
        onp->lineage = StringSave (osp->lineage);
      }
      if (! HasTaxon (orp)) {
        db = ValNodeNew (NULL);
        if (db != NULL) {
          dbt = DbtagNew ();
          if (dbt != NULL) {
            oip = ObjectIdNew ();
            if (oip != NULL) {
              oip->id = osp->taxID;
              dbt->db = StringSave ("taxon");
              dbt->tag = oip;
              db->data.ptrvalue = (Pointer) dbt;
              orp->db = db;
            }
          }
        }
      }
    }
  }
}

static BioseqPtr GetBioseqReferencedByAnnot (
  SeqAnnotPtr sap,
  Uint2 entityID
)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqFeatPtr    feat;
  SeqGraphPtr   graph;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (sap == NULL) return NULL;
  switch (sap->type) {
    case 1 :
      feat = (SeqFeatPtr) sap->data;
      while (feat != NULL) {
        slp = feat->location;
        if (slp != NULL) {
          bsp = BioseqFindFromSeqLoc (slp);
          if (bsp != NULL) return bsp;
        }
        feat = feat->next;
      }
      break;
    case 2 :
      align = (SeqAlignPtr) sap->data;
      while (align != NULL) {
        if (align->segtype == 1) {
          ddp = (DenseDiagPtr) align->segs;
          if (ddp != NULL) {
            for (sip = ddp->id; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 2) {
          dsp = (DenseSegPtr) align->segs;
          if (dsp != NULL) {
            for (sip = dsp->ids; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 3) {
          ssp = (StdSegPtr) align->segs;
          if (ssp != NULL && ssp->loc != NULL) {
            for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
              bsp = BioseqFindFromSeqLoc (tloc);
              if (bsp != NULL) return bsp;
            }
          }
        }
        align = align->next;
      }
      break;
    case 3 :
      graph = (SeqGraphPtr) sap->data;
      while (graph != NULL) {
        slp = graph->loc;
        if (slp != NULL) {
          bsp = BioseqFindFromSeqLoc (slp);
          if (bsp != NULL) return bsp;
        }
        graph = graph->next;
      }
      break;
    default :
      break;
  }
  return NULL;
}

static Int2 GetGenCodeForBsp (
  BioseqPtr bsp
)

{
  BioSourcePtr  biop;
  Boolean       mito;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SeqDescrPtr   sdp;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp == NULL) return 1;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return 1;
  orp = biop->org;
  if (orp == NULL) return 1;
  onp = orp->orgname;
  if (onp == NULL) return 1;
  mito = (Boolean) (biop->genome == 4 || biop->genome == 5);
  if (mito) {
    if (onp->mgcode == 0) {
      return 1;
    }
    return onp->mgcode;
  }
  if (onp->gcode == 0) {
    return 1;
  }
  return onp->gcode;
}

typedef struct gcmdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  feat;
  CharPtr     label;
} GmcData, PNTR GmcDataPtr;

static int LIBCALLBACK SortByGenePtr (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  GmcDataPtr gdp1, gdp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  gdp1 = (GmcDataPtr) vp1;
  gdp2 = (GmcDataPtr) vp2;
  if (gdp1 == NULL || gdp2 == NULL) return 0;

  if (gdp1->gene > gdp2->gene) return -1;
  if (gdp1->gene < gdp2->gene) return 1;

  if (gdp1->feat > gdp2->feat) return -1;
  if (gdp1->feat < gdp2->feat) return 1;

  return 0;
}

static void PrintOneGeneLine (
  SeqFeatPtr gene,
  SeqFeatPtr cds,
  SeqFeatPtr rna,
  CharPtr cdslabel,
  CharPtr rnalabel,
  FILE *fp
)

{
  BioseqPtr     bsp;
  ValNodePtr    db, old_locus_tag, vnp;
  DbtagPtr      dbt;
  CharPtr       desc, locus, locus_tag, cdslcl, cdsaccn, cdsgnl,
                rnaaccn, rnagnl, fbgn, gene_type, rna_type, prefix;
  GBQualPtr     gbq;
  GeneRefPtr    grp;
  ObjectIdPtr   oip;
  SeqIdPtr      sip;
  CharPtr       str;
  TextSeqIdPtr  tsip;

  if (fp == NULL) return; 

  locus = NULL;
  desc = NULL;
  locus_tag = NULL;
  old_locus_tag = NULL;

  cdslcl = NULL;
  cdsaccn = NULL;
  cdsgnl = NULL;
  rnaaccn = NULL;
  rnagnl = NULL;

  db = NULL;
  fbgn = NULL;

  gene_type = NULL;
  rna_type = NULL;

  if (gene != NULL) {
    gene_type = "gene";
    if (gene->pseudo) {
      gene_type = "pseudogene";
    }
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL) {
      if (grp->pseudo) {
        gene_type = "pseudogene";
      }
      locus = grp->locus;
      desc = grp->desc;
      locus_tag = grp->locus_tag;
      db = grp->db;
    }
    if (db == NULL) {
      db = gene->dbxref;
    }
    for (gbq = gene->qual; gbq != NULL; gbq = gbq->next) {
      if (StringICmp (gbq->qual, "old_locus_tag") != 0) continue;
      if (StringHasNoText (gbq->val)) continue;
      ValNodeCopyStr(&old_locus_tag, 0, gbq->val);
    }
    for (vnp = db; vnp != NULL; vnp = vnp->next) {
      dbt = (DbtagPtr) vnp->data.ptrvalue;
      if (dbt == NULL) continue;
      if (StringICmp (dbt->db, "FLYBASE") != 0) continue;
      oip = dbt->tag;
      if (oip == NULL) continue;
      fbgn = oip->str;
    }
  }

  if (cds != NULL) {
    if (cds->product != NULL) {
      bsp = BioseqFindFromSeqLoc (cds->product);
      if (bsp != NULL) {
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          switch (sip->choice) {
            case SEQID_LOCAL :
              oip = (ObjectIdPtr) sip->data.ptrvalue;
              if (oip == NULL) continue;
              cdslcl = oip->str;
              break;
            case SEQID_GENBANK :
            case SEQID_TPG :
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip == NULL) continue;
              cdsaccn = tsip->accession;
              break;
            case SEQID_GENERAL :
              dbt = (DbtagPtr) sip->data.ptrvalue;
              if (dbt == NULL) continue;
              if (IsSkippableDbtag (dbt)) continue;
              oip = dbt->tag;
              if (oip == NULL) continue;
              cdsgnl = oip->str;
              break;
            default :
              break;
          }
        }
      }
    }
  }

  if (rna != NULL) {
    switch (rna->idx.subtype) {
      case FEATDEF_preRNA :
        rna_type = "precursor RNA";
        break;
      case FEATDEF_mRNA :
        rna_type = "mRNA";
        break;
      case FEATDEF_tRNA :
        rna_type = "tRNA";
        break;
      case FEATDEF_rRNA :
        rna_type = "rRNA";
        break;
      case FEATDEF_otherRNA :
        rna_type = "misc RNA";
        break;
      case FEATDEF_ncRNA :
        rna_type = "ncRNA";
        for (gbq = rna->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "ncRNA_class") != 0) continue;
          if (StringDoesHaveText (gbq->val)) {
            rna_type = gbq->val;
          }
        }
        break;
      case FEATDEF_tmRNA :
        rna_type = "tmRNA";
        break;
      default :
        break;
    }
    if (rna->pseudo) {
      rna_type = "pseudo RNA";
    }
    if (rna->product != NULL) {
      bsp = BioseqFindFromSeqLoc (rna->product);
      if (bsp != NULL) {
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          switch (sip->choice) {
            case SEQID_GENBANK :
            case SEQID_TPG :
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip == NULL) continue;
              rnaaccn = tsip->accession;
              break;
            case SEQID_GENERAL :
              dbt = (DbtagPtr) sip->data.ptrvalue;
              if (dbt == NULL) continue;
              if (IsSkippableDbtag (dbt)) continue;
              oip = dbt->tag;
              if (oip == NULL) continue;
              rnagnl = oip->str;
              break;
            default :
              break;
          }
        }
      }
    }
  }

  if (StringDoesHaveText (locus_tag)) {
    fprintf (fp, "%s", locus_tag);
  } else {
    fprintf (fp, "null_gene_ltag");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (locus)) {
    fprintf (fp, "%s", locus);
  } else {
    fprintf (fp, "null_gene_locus");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (desc)) {
    fprintf (fp, "%s", desc);
  } else {
    fprintf (fp, "null_gene_desc");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (fbgn)) {
    fprintf (fp, "%s", fbgn);
  } else {
    fprintf (fp, "null_fbgn");
  }

  fprintf (fp, "\t");
  if (old_locus_tag != NULL) {
    prefix = "";
    for (vnp = old_locus_tag; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      fprintf (fp, "%s%s", prefix, str);
      prefix = ",";
    }
  } else {
    fprintf (fp, "null_old_ltag");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (cdslcl)) {
    fprintf (fp, "%s", cdslcl);
  } else {
    fprintf (fp, "null_cds_lcl");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (cdsaccn)) {
    fprintf (fp, "%s", cdsaccn);
  } else {
    fprintf (fp, "null_cds_accn");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (cdsgnl)) {
    fprintf (fp, "%s", cdsgnl);
  } else {
    fprintf (fp, "null_cds_gnl");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (rnaaccn)) {
    fprintf (fp, "%s", rnaaccn);
  } else {
    fprintf (fp, "null_rna_accn");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (rnagnl)) {
    fprintf (fp, "%s", rnagnl);
  } else {
    fprintf (fp, "null_rna_gnl");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (cdslabel)) {
    fprintf (fp, "%s", cdslabel);
  } else {
    fprintf (fp, "null_cds_product");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (rnalabel)) {
    fprintf (fp, "%s", rnalabel);
  } else {
    fprintf (fp, "null_rna_product");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (gene_type)) {
    fprintf (fp, "%s", gene_type);
  } else {
    fprintf (fp, "null_gene_type");
  }

  fprintf (fp, "\t");
  if (StringDoesHaveText (rna_type)) {
    fprintf (fp, "%s", rna_type);
  } else {
    fprintf (fp, "null_rna_type");
  }

  fprintf (fp, "\n");
}

static void GeneReportOneBsp (
  BioseqPtr bsp,
  FILE *fp
)

{
  CharPtr            cdslabel, rnalabel;
  SeqMgrFeatContext  fcontext;
  GmcDataPtr         gdp, head;
  GeneRefPtr         grp;
  Int2               i, j, k, numgene, numcds, numrna, total;
  SeqFeatPtr         matchsfp, sfp, tmp;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || fp == NULL) return;

  numgene = 0;
  numcds = 0;
  numrna = 0;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->data.choice) {
      case SEQFEAT_GENE :
        numgene++;
        break;
      case SEQFEAT_CDREGION :
        numcds++;
        break;
      case SEQFEAT_RNA :
        numrna++;
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  if (numgene == 0) return;
  total = numgene + numcds + numrna;
  if (total == 0) return;

  head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (total + 1));
  if (head == NULL) return;

  gdp = head;
  total = 0;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_CDREGION || sfp->data.choice == SEQFEAT_RNA) {
      gdp->feat = sfp;
      gdp->label = fcontext.label;
      grp = SeqMgrGetGeneXref (sfp);
      if (grp == NULL) {
        gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
      } else if (! SeqMgrGeneIsSuppressed (grp)) {
        if (StringDoesHaveText (grp->locus_tag)) {
          gdp->gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, NULL);
        } else if (StringDoesHaveText (grp->locus)) {
          gdp->gene = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, NULL);
        }
      }
      gdp++;
      total++;
    } else if (sfp->data.choice == SEQFEAT_GENE) {
      gdp->gene = sfp;
      gdp++;
      total++;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  HeapSort (head, (size_t) total, sizeof (GmcData), SortByGenePtr);

  for (i = 0; i < total; i += j) {
    sfp = head [i].gene;
    if (sfp == NULL) continue;
    numcds = 0;
    numrna = 0;
    for (j = 0; i + j < total && sfp == head [i + j].gene; j++) {
      tmp = head [i + j].feat;
      if (tmp == NULL) continue;
      if (tmp->data.choice == SEQFEAT_CDREGION) {
        numcds++;
      } else if (tmp->data.choice == SEQFEAT_RNA) {
        numrna++;
      }
    }
    cdslabel = NULL;
    rnalabel = NULL;
    if (numcds > 0) {
      for (k = 0; k < j; k++) {
        tmp = head [i + k].feat;
        if (tmp == NULL) continue;
        if (tmp->data.choice != SEQFEAT_CDREGION) continue;
        cdslabel = head [i + k].label;
        matchsfp = NULL;
        for (xref = tmp->xref; xref != NULL && matchsfp == NULL; xref = xref->next) {
          if (xref->id.choice != 0) {
            matchsfp = SeqMgrGetFeatureByFeatID (tmp->idx.entityID, NULL, NULL, xref, &fcontext);
            rnalabel = fcontext.label;
          }
        }
        PrintOneGeneLine (sfp, tmp, matchsfp, cdslabel, rnalabel, fp);
      }
    } else if (numrna > 0) {
      for (k = 0; k < j; k++) {
        tmp = head [i + k].feat;
        if (tmp == NULL) continue;
        if (tmp->data.choice != SEQFEAT_RNA) continue;
        rnalabel = head [i + k].label;
        PrintOneGeneLine (sfp, NULL, tmp, NULL, rnalabel, fp);
      }
    } else {
      PrintOneGeneLine (sfp, NULL, NULL, NULL, NULL, fp);
    }
  }

  MemFree (head);
}

static void GeneReportGenomicBsp (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  if (bsp == NULL) return;

  if (ISA_aa (bsp->mol)) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  if (mip->biomol != MOLECULE_TYPE_GENOMIC) return;

  GeneReportOneBsp (bsp, (FILE *) userdata);
}

static void GeneReportOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char    file [FILENAME_MAX], path [PATH_MAX];
  FILE    *fp;
  ErrSev  oldErrSev;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  fp = FileOpen (path, "w");
  if (fp == NULL) return;

  oldErrSev = ErrSetMessageLevel (SEV_MAX);
  VisitBioseqsInSep (sep, (Pointer) fp, GeneReportGenomicBsp);
  ErrSetMessageLevel (oldErrSev);

  FileClose (fp);
}

static void EnhanceOneCDS (
  SeqFeatPtr sfp,
  Boolean alt_splice
)

{
  DbtagPtr        dbt;
  GBQualPtr       gbq;
  Char            id [64];
  SeqIdPtr        ids, sip;
  size_t          len;
  CharPtr         name, nwstr, ptr, str;
  ObjectIdPtr     oip;
  ProtRefPtr      prp;
  Char            tmp [256];
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;

  name = NULL;
  vnp = NULL;
  prp = NULL;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) xref->data.value.ptrvalue;
    }
  }

  id [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "protein_id") == 0) {
      StringNCpy_0 (id, gbq->val, sizeof (id));
    }
  }
  if (StringDoesHaveText (id) && StringChr (id, '|') != NULL) {
    str = NULL;
    ids = SeqIdParse (id);
    for (sip = ids; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag (dbt)) continue;
      oip = dbt->tag;
      if (oip == NULL) continue;
      str = oip->str;
    }

    if (StringDoesHaveText (str)) {
      if (prp != NULL && prp->name != NULL) {
        vnp = prp->name;
        name = (CharPtr) vnp->data.ptrvalue;
      }
      if (StringDoesHaveText (name) && vnp != NULL) {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, isoform %s", str, ptr);
            len = StringLen (name) + StringLen (", ") + StringLen (tmp);
            nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
            if (nwstr != NULL) {
              StringCpy (nwstr, name);
              /*
              StringCat (nwstr, ", ");
              */
              StringCat (nwstr, " ");
              StringCat (nwstr, tmp);
              vnp->data.ptrvalue = (Pointer) nwstr;
              MemFree (name);
            }
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          len = StringLen (name) + StringLen (", ") + StringLen (str);
          nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
          if (nwstr != NULL) {
            StringCpy (nwstr, name);
            /*
            StringCat (nwstr, ", ");
            */
            StringCat (nwstr, " ");
            StringCat (nwstr, str);
            vnp->data.ptrvalue = (Pointer) nwstr;
            MemFree (name);
          }
        }
      } else {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, isoform %s", str, ptr);
            AddQualifierToFeature (sfp, "product", tmp);
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          AddQualifierToFeature (sfp, "product", str);
        }
      }
    }

    SeqIdSetFree (ids);
  }
}

static void EnhanceOneRna (
  SeqFeatPtr sfp,
  Boolean alt_splice
)

{
  DbtagPtr     dbt;
  GBQualPtr    gbq, nm_gbq;
  Char         id [64];
  SeqIdPtr     ids, sip;
  size_t       len;
  CharPtr      name, nwstr, ptr, str;
  ObjectIdPtr  oip;
  RnaRefPtr    rrp;
  Char         tmp [256];

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  name = NULL;
  nm_gbq = NULL;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp != NULL && rrp->ext.choice == 1) {
    switch (rrp->type) {
      case 1 :  /* precurrsor_RNA */
      case 2 :  /* mRNA */
      case 4 :  /* rRNA */
        name = rrp->ext.value.ptrvalue;
        break;
      case 255 :  /* misc_RNA, ncRNA, tmRNA */
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "product") == 0) {
            nm_gbq = gbq;
            name = gbq->val;
          }
        }
        break;
      case 3:  /* tRNA */
        return;
      default :
        break;
    }
  }

  id [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") == 0) {
      StringNCpy_0 (id, gbq->val, sizeof (id));
    }
  }
  if (StringDoesHaveText (id) && StringChr (id, '|') != NULL) {
    str = NULL;
    ids = SeqIdParse (id);
    for (sip = ids; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag(dbt)) continue;
      oip = dbt->tag;
      if (oip == NULL) continue;
      str = oip->str;
    }

    if (StringDoesHaveText (str)) {
      if (StringDoesHaveText (name) && StringCmp (str, name) != 0) {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, transcript variant %s", str, ptr);
            len = StringLen (name) + StringLen (", ") + StringLen (tmp);
            nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
            if (nwstr != NULL) {
              StringCpy (nwstr, name);
              /*
              StringCat (nwstr, ", ");
              */
              StringCat (nwstr, " ");
              StringCat (nwstr, tmp);
              if (nm_gbq != NULL) {
                nm_gbq->val = (Pointer) nwstr;
              } else {
                rrp->ext.value.ptrvalue = (Pointer) nwstr;
              }
              MemFree (name);
            }
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          len = StringLen (name) + StringLen (", ") + StringLen (str);
          nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
          if (nwstr != NULL) {
            StringCpy (nwstr, name);
            /*
            StringCat (nwstr, ", ");
            */
            StringCat (nwstr, " ");
            StringCat (nwstr, str);
            if (nm_gbq != NULL) {
              nm_gbq->val = (Pointer) nwstr;
            } else {
              rrp->ext.value.ptrvalue = (Pointer) nwstr;
            }
            MemFree (name);
          }
        }
      } else {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, transcript variant %s", str, ptr);
            AddQualifierToFeature (sfp, "product", tmp);
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          AddQualifierToFeature (sfp, "product", str);
        }
      }
    }

    SeqIdSetFree (ids);
  }
}

static void EnhanceFeatureAnnotation (
  SeqFeatPtr features,
  BioseqPtr bsp
)

{
  GmcDataPtr  gdp, head;
  GeneRefPtr  grp;
  Int2        i, j, k, numgene, numcds, numrna;
  SeqFeatPtr  sfp;

  if (features == NULL || bsp == NULL) return;

  numgene = 0;
  numcds = 0;
  numrna = 0;

  for (sfp = features; sfp != NULL; sfp = sfp->next) {
    switch (sfp->data.choice) {
      case SEQFEAT_GENE :
        numgene++;
        break;
      case SEQFEAT_CDREGION :
        numcds++;
        break;
      case SEQFEAT_RNA :
        numrna++;
        break;
      default :
        break;
    }
  }

  if (numgene == 0) return;

  if (numcds > 0) {
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numcds + 1));
    if (head != NULL) {
      gdp = head;
      for (sfp = features; sfp != NULL; sfp = sfp->next) {
        if (sfp->idx.subtype == FEATDEF_CDS) {
          gdp->feat = sfp;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          }
          gdp++;
        }
      }
      HeapSort (head, (size_t) numcds, sizeof (GmcData), SortByGenePtr);
      for (i = 0; i < numcds; i += j) {
        sfp = head [i].gene;
        for (j = 1; i + j < numcds && sfp == head [i + j].gene; j++) continue;
        if (j == 1) {
          /* no alt splicing */
          EnhanceOneCDS (head [i].feat, FALSE);
        } else {
          /* is alt splicing */
          for (k = 0; k < j; k++) {
            EnhanceOneCDS (head [i + k].feat, TRUE);
          }
        }
      }
    }
    MemFree (head);
  }

  if (numrna > 0) {
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numrna + 1));
    if (head != NULL) {
      gdp = head;
      for (sfp = features; sfp != NULL; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          gdp->feat = sfp;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          }
          gdp++;
        }
      }
      HeapSort (head, (size_t) numrna, sizeof (GmcData), SortByGenePtr);
      for (i = 0; i < numrna; i += j) {
        sfp = head [i].gene;
        for (j = 1; i + j < numrna && sfp == head [i + j].gene; j++) continue;
        if (j == 1) {
          /* no alt splicing */
          EnhanceOneRna (head [i].feat, FALSE);
        } else {
          /* is alt splicing */
          for (k = 0; k < j; k++) {
            EnhanceOneRna (head [i + k].feat, TRUE);
          }
        }
      }
    }
    MemFree (head);
  }
}

static BioseqPtr AttachSeqAnnotEntity (
  Uint2 entityID,
  SeqAnnotPtr sap,
  TblArgsPtr tbl
)

{
  SeqAnnotPtr  anp;
  BioseqPtr    bsp;
  Char         buf [80];
  Int2         genCode;
  SeqEntryPtr  oldscope;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp = NULL;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  if (sap == NULL || tbl == NULL) return NULL;

  bsp = GetBioseqReferencedByAnnot (sap, entityID);
  if (bsp == NULL) {
    oldscope = SeqEntrySetScope (NULL);
    if (oldscope != NULL) {
      bsp = GetBioseqReferencedByAnnot (sap, entityID);
      SeqEntrySetScope (oldscope);
    }
  }
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      genCode = GetGenCodeForBsp (bsp);
      SetEmptyGeneticCodes (sap, genCode);
    }
    if (bsp->annot == NULL) {
      bsp->annot = sap;
    } else {
      anp = bsp->annot;
      while (anp->next != NULL) {
        anp = anp->next;
      }
      anp->next = sap;
    }
    if (sfp != NULL) {
      if (tbl->smartfeats) {

        /* indexing needed to find mRNA and CDS within each gene */

        SeqMgrIndexFeatures (entityID, NULL);

        EnhanceFeatureAnnotation (sfp, bsp);
      }

      PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, tbl->genprodset, tbl->forcelocalid);
      sep = GetTopSeqEntryForEntityID (entityID);
    }
  } else {
    buf [0] = '\0';
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      if (sfp != NULL && sfp->location != NULL) {
        slp = SeqLocFindNext (sfp->location, NULL);
        if (slp != NULL) {
          sip = SeqLocId (slp);
          if (sip != NULL) {
            SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          }
        }
      }
    }
    Message (MSG_POSTERR, "Feature table identifiers %s do not match record", buf);
  }
  sep = GetTopSeqEntryForEntityID (entityID);
  return bsp;
}

static CharPtr TrimBracketsFromString (
  CharPtr str,
  SqnTagPtr stp
)

{
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  Int2     count;
  CharPtr  dst;
  CharPtr  ptr;

  if (StringHasNoText (str) || stp == NULL) return str;

  /* remove bracketed fields */

  count = 0;
  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '[') {
      if (count < stp->num_tags && (! stp->used [count])) {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
        while (ch != '\0' && ch != ']') {
          *dst = ch;
          dst++;
          ptr++;
          ch = *ptr;
        }
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
        while (ch != '\0' && ch != ']' && ch != '"') {
          ptr++;
          ch = *ptr;
        }
        if (ch == '"') {
          ptr++;
          ch = *ptr;
          while (ch != '\0' && ch != '"') {
            ptr++;
            ch = *ptr;
          }
        }
        while (ch != '\0' && ch != ']') {
          ptr++;
          ch = *ptr;
        }
        ptr++;
        ch = *ptr;
      }
      count++;
    } else {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
  }
  *dst = '\0';

  /* remove runs of whitespace characters */

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_WHITESP (ch)) {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
      while (IS_WHITESP (ch)) {
        ptr++;
        ch = *ptr;
      }
    } else {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
  }
  *dst = '\0';

  return str;
}

static Boolean HasTpaAccession (
  UserObjectPtr uop
)

{
  UserFieldPtr  curr;
  ObjectIdPtr   oip;
  CharPtr       str;
  UserFieldPtr  ufp;

  if (uop == NULL) return FALSE;
  if ((oip = uop->type) == NULL) return FALSE;
  if (StringCmp (oip->str, "TpaAssembly") != 0) return FALSE;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringDoesHaveText (str)) return TRUE;
    }
  }

  return FALSE;
}

static Boolean HasGenomeProjectDB (
  UserObjectPtr uop
)

{
  UserFieldPtr  curr;
  ObjectIdPtr   oip;
  Int4          val;

  if (uop == NULL) return FALSE;
  if ((oip = uop->type) == NULL) return FALSE;
  if (StringCmp (oip->str, "GenomeProjectsDB") != 0) return FALSE;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip == NULL || StringICmp (oip->str, "ProjectID") != 0) continue;
    if (curr->choice != 2) continue;
    val = (Int4) curr->data.intvalue;
    if (val > 0) return TRUE;
  }

  return FALSE;
}

static void GetFirstBiop (
  BioSourcePtr biop,
  Pointer userdata
)

{
  BioSourcePtr PNTR biopp;

  biopp = (BioSourcePtr PNTR) userdata;
  if (biop == NULL || biopp == NULL) return;
  if (*biopp != NULL) return;
  *biopp = biop;
}

static void ProcessOneNuc (
  Uint2 entityID,
  BioseqPtr bsp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  MolInfoPtr template_molinfo
)

{
  Boolean        addNewBiop = TRUE;
  Boolean        addNewMip = TRUE;
  BioSourcePtr   biop = NULL;
  SeqFeatPtr     cds;
  GBBlockPtr     gbp;
  Int2           genCode;
  size_t         len;
  MolInfoPtr     mip = NULL;
  Boolean        mito;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  SeqDescrPtr    sdp;
  SeqHistPtr     shp;
  SqnTagPtr      stp = NULL;
  CharPtr        str;
  CharPtr        tmp;
  CharPtr        ttl = NULL;
  UserObjectPtr  uop;
  ValNodePtr     vnp;
  SeqMgrDescContext dcontext;

  if (bsp == NULL) return;

  genCode = GetGenCodeForBsp (bsp);

  if (bsp->mol == Seq_mol_na) {
    bsp->mol = Seq_mol_dna;
  }

  if (src != NULL) {
    src = AsnIoMemCopy ((Pointer) src,
                        (AsnReadFunc) BioSourceAsnRead,
                        (AsnWriteFunc) BioSourceAsnWrite);
  } else {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (sdp != NULL) {
      src = sdp->data.ptrvalue;
      if (src != NULL) {
        addNewBiop = FALSE;
      }
    }
  }

  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
  if (vnp != NULL) {
    ttl = (CharPtr) vnp->data.ptrvalue;
  }

  if (ttl != NULL || tbl->srcquals != NULL) {
    len = StringLen (ttl) + StringLen (tbl->srcquals) + 5;
    str = (CharPtr) MemNew (len * sizeof (Char));
    if (str != NULL) {
      StringCpy (str, ttl);
      if (ttl != NULL && tbl->srcquals != NULL) {
        StringCat (str, "; ");
      }
      StringCat (str, tbl->srcquals);
      stp = SqnTagParse (str);
    }
    MemFree (str);
  }

  if (stp != NULL) {
    biop = ParseTitleIntoBioSource (stp, tbl->organism, src);
    ParseTitleIntoBioseq (stp, bsp);
    str = SqnTagFind (stp, "comment");
    if (str != NULL) {
      tmp = StringSave (str);
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_comment, (Pointer) tmp);
    }
  }
  if (biop == NULL) {
    biop = ParseTitleIntoBioSource (NULL, tbl->organism, src);
  }
  if (biop != NULL && addNewBiop) {
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_source, (Pointer) biop);
  }
  if (biop != NULL) {
    AddMissingSourceInfo (biop);
  }

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    addNewMip = FALSE;
  } else {
    mip = MolInfoNew ();
  }
  if (mip != NULL) {
    if (stp != NULL) {
      mip = ParseTitleIntoMolInfo (stp, mip);
    }
    if (mip->biomol == 0 && template_molinfo != NULL)
    {
      mip->biomol = template_molinfo->biomol;
    }
    if (mip->biomol == 0) {
      mip->biomol = MOLECULE_TYPE_GENOMIC;
    }
    if (addNewMip) {
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_molinfo, (Pointer) mip);
    }
    switch (mip->biomol) {
      case MOLECULE_TYPE_PRE_MRNA :
      case MOLECULE_TYPE_MRNA :
      case MOLECULE_TYPE_RRNA :
      case MOLECULE_TYPE_TRNA :
      case MOLECULE_TYPE_SNRNA :
      case MOLECULE_TYPE_SCRNA :
      case MOLECULE_TYPE_CRNA :
      case MOLECULE_TYPE_SNORNA :
      case MOLECULE_TYPE_TRANSCRIBED_RNA :
      case MOLECULE_TYPE_NCRNA :
      case MOLECULE_TYPE_TMRNA :
        if (bsp->mol == Seq_mol_dna) {
          str = SqnTagFind (stp, "molecule");
          if (str == NULL) {
            str = SqnTagFind (stp, "mol");
          }
          if (str != NULL) {
            if (StringICmp (str, "dna") == 0) break;
          }
          bsp->mol = Seq_mol_rna;
        }
        break;
      default :
        break;
    }
  }

  if (genCode == 0 && biop != NULL) {
    orp = biop->org;
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        mito = (Boolean) (biop->genome == 4 || biop->genome == 5);
        if (mito) {
          genCode = onp->mgcode;
        } else {
          genCode = onp->gcode;
        }
      }
    }
  }

  if (StringDoesHaveText (tbl->comment)) {
    str = StringSave (tbl->comment);
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_comment, (Pointer) str);
  }
  if (StringDoesHaveText (tbl->commentFile)) {
    str = StringSave (tbl->commentFile);
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_comment, (Pointer) str);
  }

  if (stp != NULL) {
    gbp = ParseTitleIntoGenBank (stp, NULL);
    if (gbp != NULL && (gbp->extra_accessions != NULL || gbp->keywords != NULL)) {
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_genbank, (Pointer) gbp);
    } else {
      gbp = GBBlockFree (gbp);
    }

    shp = ParseTitleIntoSeqHist (stp, NULL);
    if (shp != NULL && shp->replace_ids != NULL) {
      bsp->hist = SeqHistFree (bsp->hist);
      bsp->hist = shp;
    } else {
      shp = SeqHistFree (shp);
    }
  }

  if (stp != NULL) {
    uop = ParseTitleIntoTpaAssembly (stp, NULL);
    if (uop != NULL && HasTpaAccession (uop)) {
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_user, (Pointer) uop);
    } else {
      uop = UserObjectFree (uop);
    }
  }

  if (stp != NULL) {
    uop = ParseTitleIntoGenomeProjectsDB (stp, NULL);
    if (uop != NULL && HasGenomeProjectDB (uop)) {
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_user, (Pointer) uop);
    } else {
      uop = UserObjectFree (uop);
    }
  }

  /* look for pubmed IDs */
  if (stp != NULL) {
    AddPubsFromTitle (stp, &(bsp->descr)); 
  }

  if (tbl->findorf) {
    cds = AnnotateBestOrf (bsp, genCode, tbl->altstart, tbl->runonorf, stp);
    if (cds != NULL) {
      PromoteXrefsExEx (cds, bsp, entityID, TRUE, FALSE, FALSE, tbl->forcelocalid);
    }
  }

  TrimBracketsFromString (ttl, stp);
  if (StringDoesHaveText (ttl)) {
    str = StringSave (ttl);
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
  }

  if (stp != NULL) {
    SqnTagFree (stp);
  }

  ValNodeFreeData (vnp);
}

static void ProcessNucBioseqs (SeqEntryPtr top_sep, Uint2 entityID, BioSourcePtr src, TblArgsPtr tbl, MolInfoPtr template_molinfo)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqEntryPtr sep;

  if (top_sep == NULL || top_sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (top_sep)) {
    bsp = (BioseqPtr) top_sep->data.ptrvalue;
    if (!ISA_aa (bsp->mol)) {
      ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);
    }
  } else if (IS_Bioseq_set (top_sep)) {
    bssp = (BioseqSetPtr) top_sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      ProcessNucBioseqs (sep, entityID, src, tbl, template_molinfo);
    }
  }
}


static void ProcessOneAnnot (
  SeqAnnotPtr sap,
  Uint2 entityID,
  TblArgsPtr tbl
)

{
  BioseqPtr   bsp;
  Int2        genCode;
  SeqFeatPtr  sfp;
  SeqEntryPtr sep;

  if (sap == NULL || tbl == NULL) return;

  bsp = AttachSeqAnnotEntity (entityID, sap, tbl);
  if (bsp == NULL) return;

    sep = GetTopSeqEntryForEntityID (entityID);

  /* correct all idx parent pointers */

  AssignIDsInEntity (entityID, 0, NULL);

  genCode = GetGenCodeForBsp (bsp);

  /* coercion of SeqIds to accession moved to ProcessOneRecord->MakeAccessionID */

  /* for parsed in features or best ORF, promote CDS products to protein bioseq */

  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 1) {
      SetEmptyGeneticCodes (sap, genCode);
      sfp = (SeqFeatPtr) sap->data;
      PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, tbl->genprodset, tbl->forcelocalid);
    }
  }
  sep = GetTopSeqEntryForEntityID (entityID);
}

static void UpdateException (
  SeqFeatPtr sfp,
  CharPtr text
)

{
  size_t   len;
  CharPtr  str;

  if (sfp == NULL) return;

  sfp->excpt = TRUE;

  if (sfp->except_text == NULL) {
    sfp->except_text = StringSave (text);
  } else {
    len = StringLen (sfp->except_text) + StringLen (text) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, sfp->except_text);
    StringCat (str, ",");
    StringCat (str, text);
    sfp->except_text = MemFree (sfp->except_text);
    sfp->except_text = str;
  }
}

static void ReplaceOnePeptide (
  SimpleSeqPtr ssp,
  Boolean conflict,
  Boolean genprodset
)

{
  Uint1              aa;
  ByteStorePtr       bs;
  BioseqPtr          bsp, gen;
  SeqFeatPtr         cds;
  CdRegionPtr        crp;
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  SeqFeatPtr         prt;
  SeqDescrPtr        sdp;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  CharPtr            str, str1, str2;
  ValNodePtr         vnp;

  if (ssp == NULL || ssp->numid < 1) return;

  str = ssp->id [0];
  if (StringHasNoText (str)) {
    str = "?";
  }
  sip = MakeSeqID (str);
  bsp = BioseqFind (sip);
  SeqIdFree (sip);
  if (bsp == NULL) {
    Message (MSG_POSTERR, "Unable to find protein sequence %s", str);
  }
  if (bsp == NULL || bsp->repr != Seq_repr_raw) return;

  if (bsp->seq_data_type == Seq_code_gap) return;

  if (! ISA_aa (bsp->mol)) {
    Message (MSG_POSTERR, "Will not replace mRNA sequence %s with protein", str);
    return;
  }

  /* remove trailing X and * - now just trailing star */

  bs = ssp->seq;
  BSSeek (bs, -1, SEEK_END);
  aa = (Uint1) BSGetByte (bs);
  while (( /* aa == 'X' || */ aa == '*') && ssp->seqlen > 0) {
    BSSeek (bs, -1, SEEK_END);
    BSDelete (bs, 1);
    BSSeek (bs, -1, SEEK_END);
    aa = (Uint1) BSGetByte (bs);
  }
  ssp->seqlen = BSLen (bs);

  str1 = BSMerge (ssp->seq, NULL);
  str2 = BSMerge ((ByteStorePtr) bsp->seq_data, NULL);

  if (StringCmp (str1, str2) != 0) {

    /* swap sequence byte stores */

    bs = (ByteStorePtr) bsp->seq_data;
    bsp->seq_data = (SeqDataPtr) ssp->seq;
    ssp->seq = bs;
    bsp->length = BSLen ((ByteStorePtr) bsp->seq_data);
    bsp->seq_data_type = Seq_code_ncbieaa;

    if (genprodset) {

      /* SeqMgrGetCDSgivenProduct here would return CDS within nuc-prot set, not genomic */

      for (vnp = SeqMgrGetSfpProductList (bsp); vnp != NULL; vnp = vnp->next) {
        cds = (SeqFeatPtr) vnp->data.ptrvalue;
        if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) continue;
        gen = BioseqFindFromSeqLoc (cds->location);
        if (gen == NULL) continue;

        sdp = SeqMgrGetNextDescriptor (gen, NULL, Seq_descr_molinfo, &dcontext);
        if (sdp == NULL) continue;
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip == NULL) continue;

        if (mip->biomol == MOLECULE_TYPE_GENOMIC) {

          UpdateException (cds, "translated product replaced");

        } else if (mip->biomol == MOLECULE_TYPE_MRNA) {

          crp = (CdRegionPtr) cds->data.value.ptrvalue;
          if (crp != NULL && conflict) {

            /* mark CDS in nuc-prot set for coordinate adjustment */

            crp->conflict = TRUE;
          }
        }
      }

    } else {

      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
      if (cds != NULL) {
        UpdateException (cds, "translated product replaced");
      }
    }

    prt = SeqMgrGetBestProteinFeature (bsp, NULL);
    if (prt != NULL) {
      slp = prt->location;
      if (slp != NULL && slp->choice == SEQLOC_INT) {
        sintp = (SeqIntPtr) slp->data.ptrvalue;
        if (sintp != NULL) {
          sintp->to = bsp->length - 1;
        }
      }
    }
  }

  MemFree (str1);
  MemFree (str2);
}

static void ReplaceOneRNA (
  SimpleSeqPtr ssp,
  Boolean conflict
)

{
  ByteStorePtr       bs;
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds, mrna;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  CharPtr            str, str1, str2;

  if (ssp == NULL || ssp->numid < 1) return;

  str = ssp->id [0];
  if (StringHasNoText (str)) {
    str = "?";
  }
  sip = MakeSeqID (str);
  bsp = BioseqFind (sip);
  SeqIdFree (sip);
  if (bsp == NULL) {
    Message (MSG_POSTERR, "Unable to find mRNA sequence %s", str);
  }
  if (bsp == NULL || bsp->repr != Seq_repr_raw) return;
  if (! ISA_na (bsp->mol)) {
    Message (MSG_POSTERR, "Will not replace protein sequence %s with mRNA", str);
    return;
  }

  /* remove trailing X and * */

  bs = ssp->seq;
  ssp->seqlen = BSLen (bs);

  str1 = BSMerge (ssp->seq, NULL);
  str2 = GetSequenceByBsp (bsp);

  if (StringCmp (str1, str2) != 0) {

    /* swap sequence byte stores */

    bs = (ByteStorePtr) bsp->seq_data;
    bsp->seq_data = (SeqDataPtr) ssp->seq;
    ssp->seq = bs;
    bsp->length = BSLen ((ByteStorePtr) bsp->seq_data);
    bsp->seq_data_type = Seq_code_iupacna;
    BioseqPack (bsp);

    mrna = SeqMgrGetRNAgivenProduct (bsp, NULL);
    if (mrna != NULL) {
      UpdateException (mrna, "transcribed product replaced");

      /*
      if (conflict) {
        mrna->excpt = TRUE;
        if (StringHasNoText (mrna->except_text)) {
          mrna->except_text = StringSave ("RNA editing");
        }
      }
      */
    }

    /* make sure CDS in nuc-prot set is not longer than just-replaced RNA */

    cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
    if (cds != NULL) {
      slp = cds->location;
      if (slp != NULL && slp->choice == SEQLOC_INT) {
        sintp = (SeqIntPtr) slp->data.ptrvalue;
        if (sintp != NULL) {
          if (sintp->from == 0 && sintp->to > bsp->length - 1) {
            sintp->to = bsp->length - 1;
          }
        }
      }
    }
  }

  MemFree (str1);
  MemFree (str2);
}

static SeqLocPtr PredictOneCodingRegion (BioseqPtr nucbsp, BioseqPtr protbsp, Int2 genCode)

{
  BioseqPtr    bsp;
  SeqLocPtr    oldslp;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  slp = NULL;
  sap = SuggestCodingRegion (nucbsp, protbsp, genCode);
  if (sap != NULL && sap->type == 1) {
    sfp = (SeqFeatPtr) sap->data;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      slp = sfp->location;
      sfp->location = NULL;
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          if (bsp->repr == Seq_repr_seg) {
            oldslp = slp;
            slp = SegLocToParts (bsp, oldslp);
            FreeAllFuzz (slp);
            SeqLocFree (oldslp);
          }
        }
      }
    }
  }
  sap = SeqAnnotFree (sap);
  StripLocusFromSeqLoc (slp);
  return slp;
}

static void SuggestOnePeptide (
  BioseqPtr nucbsp,
  BioseqPtr protbsp,
  Int2 genCode
)

{
  SeqFeatPtr   cds;
  CdRegionPtr  crp;
  SeqFeatPtr   gene;
  GeneRefPtr   grp;
  Boolean      partial5;
  Boolean      partial3;
  ProtRefPtr   prp;
  SeqFeatPtr   prt;
  SeqLocPtr    slp;
  SqnTagPtr    stp;
  CharPtr      ttl;
  ValNodePtr   vnp;

  if (nucbsp == NULL || protbsp == NULL) return;
  slp = PredictOneCodingRegion (nucbsp, protbsp, genCode);
  if (slp == NULL) return;

  crp = CreateNewCdRgn (0, FALSE, genCode);
  if (crp != NULL) {
    CheckSeqLocForPartial (slp, &partial5, &partial3);

    cds = CreateNewFeatureOnBioseq (nucbsp, SEQFEAT_CDREGION, slp);
    if (cds != NULL) {
      cds->data.value.ptrvalue = (Pointer) crp;
      cds->partial |= partial5 | partial3;
      SetSeqFeatProduct (cds, protbsp);
    }

    if (protbsp->descr != NULL) {
      vnp = ValNodeExtract (&(protbsp->descr), Seq_descr_title);
      if (vnp != NULL) {
        ttl = (CharPtr) vnp->data.ptrvalue;
        if (ttl != NULL) {
          stp = SqnTagParse (ttl);
          if (stp != NULL) {

            prp = ProtRefNew ();
            prp = ParseTitleIntoProtRef (stp, prp);
            if (prp != NULL) {
              if (prp->name == NULL && prp->desc == NULL) {
                prp->name = ValNodeCopyStr (NULL, 0, "unknown");
              }
              prt = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
              if (prt != NULL) {
                prt->data.value.ptrvalue = (Pointer) prp;
                prt->partial |= partial5 | partial3;
              }
            }

            grp = GeneRefNew ();
            grp = ParseTitleIntoGeneRef (stp, grp);
            if (grp != NULL) {
              if (grp->locus == NULL && grp->syn == NULL) {
                GeneRefFree (grp);
              } else {
                gene = CreateNewFeatureOnBioseq (nucbsp, SEQFEAT_GENE, NULL);
                if (gene != NULL) {
                  gene->data.value.ptrvalue = (Pointer) grp;
                  gene->partial |= partial5 | partial3;
                  gene->location = SeqLocFree (gene->location);
                  gene->location = SeqLocMerge (nucbsp, slp, NULL, TRUE, TRUE, TRUE);
                }
              }
            }

            SqnTagFree (stp);
          }
        }

        ValNodeFreeData (vnp);
      }
    }
  }

  SeqLocFree (slp);
}

static void RnaProtTrailingCommaFix (SeqFeatPtr sfp, Pointer userdata)

{
  Char        ch;
  size_t      len;
  ProtRefPtr  prp;
  RnaRefPtr   rrp;
  CharPtr     str;
  ValNodePtr  vnp;

  if (sfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = StringLen (str);
      if (len < 1) continue;
      ch = str [len - 1];
      while (ch == ' ' && len > 2) {
        len--;
        ch = str [len - 1];
      }
      if (ch == ',') {
        str [len - 1] = '_';
        str [len] = '\0';
      }
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    if (rrp->ext.choice == 1) {
      str = rrp->ext.value.ptrvalue;
      if (StringDoesHaveText (str)) {
        len = StringLen (str);
        if (len > 0) {
          ch = str [len - 1];
          while (ch == ' ' && len > 2) {
            len--;
            ch = str [len - 1];
          }
          if (ch == ',') {
            str [len - 1] = '_';
            str [len] = '\0';
          }
        }
      }
    }
  }
}

static Uint2 ProcessOneAsn (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  CharPtr localname,
  SeqEntryPtr gsep,
  MolInfoPtr template_molinfo
)

{
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp;
  Pointer        dataptr;
  Uint2          datatype, entityID;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;

  if (fp == NULL) return 0;

  if (gsep != NULL) {
    bssp = (BioseqSetPtr) gsep->data.ptrvalue;
    if (bssp == NULL) return 0;

    SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
    GetSeqEntryParent (gsep, &parentptr, &parenttype);

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, TRUE, FALSE, TRUE, FALSE);
    if (datatype == OBJ_BIOSEQ) {
      bssp->seq_set = SeqMgrGetSeqEntryForData (dataptr);
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) dataptr, gsep);
    } else if (datatype == OBJ_BIOSEQSET) {
      bssp->seq_set = SeqMgrGetSeqEntryForData (dataptr);
      SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) dataptr, gsep);
    } else if (datatype == OBJ_SEQENTRY) {
      sep = (SeqEntryPtr) dataptr;
      bssp->seq_set = sep;
      if (IS_Bioseq (sep)) {
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) sep->data.ptrvalue, gsep);
      } else if (IS_Bioseq_set (sep)) {
        SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) sep->data.ptrvalue, gsep);
      } else return 0;
    } else return 0;

    SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
  } else {
    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, TRUE, FALSE, TRUE, FALSE);
  }
  if (dataptr == NULL) return 0;

  sep = GetTopSeqEntryForEntityID (entityID);
  bsp = FindNucBioseq (sep);
  if (bsp == NULL) {
    ObjMgrFreeByEntityID (entityID);
    return 0;
  }

  VisitFeaturesInSep (sep, NULL, RnaProtTrailingCommaFix);

  if (StringDoesHaveText (localname)) {
    sip = MakeSeqID (localname);
    if (sip != NULL) {
      bsp->id = SeqIdSetFree (bsp->id);
      bsp->id = sip;
      SeqMgrReplaceInBioseqIndex (bsp);
      VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
    }
  }

  ProcessNucBioseqs (sep, entityID, src, tbl, template_molinfo);

  return entityID;
}

static Uint2 ProcessRaw2Delt (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  CharPtr localname,
  SeqEntryPtr gsep,
  MolInfoPtr template_molinfo
)

{
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp;
  Pointer        dataptr;
  Uint2          datatype, entityID;
  Int4           gap_sizes [2];
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;

  if (fp == NULL) return 0;

  if (gsep != NULL) {
    bssp = (BioseqSetPtr) gsep->data.ptrvalue;
    if (bssp == NULL) return 0;

    SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
    GetSeqEntryParent (gsep, &parentptr, &parenttype);

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, TRUE, FALSE, TRUE, FALSE);
    if (datatype == OBJ_BIOSEQ) {
      bssp->seq_set = SeqMgrGetSeqEntryForData (dataptr);
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) dataptr, gsep);
    } else if (datatype == OBJ_BIOSEQSET) {
      bssp->seq_set = SeqMgrGetSeqEntryForData (dataptr);
      SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) dataptr, gsep);
    } else if (datatype == OBJ_SEQENTRY) {
      sep = (SeqEntryPtr) dataptr;
      bssp->seq_set = sep;
      if (IS_Bioseq (sep)) {
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) sep->data.ptrvalue, gsep);
      } else if (IS_Bioseq_set (sep)) {
        SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) sep->data.ptrvalue, gsep);
      } else return 0;
    } else return 0;

    SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
  } else {
    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, TRUE, FALSE, TRUE, FALSE);
  }
  if (dataptr == NULL) return 0;

  sep = GetTopSeqEntryForEntityID (entityID);
  bsp = FindNucBioseq (sep);
  if (bsp == NULL) {
    ObjMgrFreeByEntityID (entityID);
    return 0;
  }

  VisitFeaturesInSep (sep, NULL, RnaProtTrailingCommaFix);

  if (StringDoesHaveText (localname)) {
    sip = MakeSeqID (localname);
    if (sip != NULL) {
      bsp->id = SeqIdSetFree (bsp->id);
      bsp->id = sip;
      SeqMgrReplaceInBioseqIndex (bsp);
      VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
    }
  }

  if (bsp->repr == Seq_repr_raw) {
    if (tbl->r2dunk100) {
      gap_sizes [0] = 100;
    } else {
      gap_sizes [0] = 0;
    }
    gap_sizes [1] = -(tbl->r2dmin);

    ConvertNsToGaps (bsp, gap_sizes);
  }

  ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);

  return entityID;
}

static Uint2 ProcessGappedSet (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  SeqEntryPtr gsep,
  MolInfoPtr template_molinfo
)

{
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp;
  Uint2          entityID;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep;

  if (fp == NULL) return 0;

  if (gsep != NULL) {
    bssp = (BioseqSetPtr) gsep->data.ptrvalue;
    if (bssp == NULL) return 0;

    SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
    GetSeqEntryParent (gsep, &parentptr, &parenttype);

    bsp = ReadDeltaFasta (fp, NULL);
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bsp);
      bssp->seq_set = sep;
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, gsep);
    } else return 0;

    SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
  } else {
    bsp = ReadDeltaFasta (fp, NULL);
    if (bsp != NULL) {
      entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
    }
  }
  if (bsp == NULL) return 0;

  sep = GetTopSeqEntryForEntityID (entityID);
  bsp = FindNucBioseq (sep);
  if (bsp == NULL) {
    ObjMgrFreeByEntityID (entityID);
    return 0;
  }

  VisitFeaturesInSep (sep, NULL, RnaProtTrailingCommaFix);

  ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);

  return entityID;
}

typedef struct resqseqgph {
  Int2         index;
  SeqGraphPtr  sgp;
} ResqSeqgph, PNTR ResqSeqgphPtr;

static void RescueSeqGraphs (
  BioseqPtr bsp,
  Int2 index,
  ValNodePtr PNTR vnpp
)

{
  SeqAnnotPtr    nextsap;
  SeqGraphPtr    nextsgp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsgp;
  ResqSeqgphPtr  rsp;
  SeqAnnotPtr    sap;
  SeqGraphPtr    sgp;

  if (bsp == NULL || vnpp == NULL) return;
  sap = bsp->annot;
  prevsap = (Pointer PNTR) &(bsp->annot);
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        *(prevsgp) = sgp->next;
        sgp->next = NULL;
        rsp = (ResqSeqgphPtr) MemNew (sizeof (ResqSeqgph));
        rsp->index = index;
        rsp->sgp = sgp;
        ValNodeAddPointer (vnpp, 0, (Pointer) rsp);
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static SeqAnnotPtr NewSeqAnnotType3 (
  CharPtr name,
  SeqGraphPtr sgp
)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (StringDoesHaveText (name)) {
    SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static void OffsetAndLinkSeqGraph (
  BioseqPtr bsp,
  SeqGraphPtr sgp,
  Int2 index
)

{
  DeltaSeqPtr  dsp;
  SeqGraphPtr  lastsgp;
  Int4         len;
  SeqLitPtr    litp;
  SeqAnnotPtr  sap;
  SeqIntPtr    sintp;
  SeqLocPtr    slp;

  if (bsp == NULL || sgp == NULL || index < 1) return;
  len = 0;
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
         dsp != NULL && index > 1; dsp = dsp->next, index--) {
      if (dsp->choice == 1) {
        len += SeqLocLen ((SeqLocPtr) dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          len += litp->length;
        }
      }
    }
  }
  slp = sgp->loc;
  if (slp != NULL && slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp != NULL) {
      sintp->from += len;
      sintp->to += len;
      sintp->id = SeqIdFree (sintp->id);
      sintp->id = SeqIdDup (bsp->id);
    }
  }
  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 3) {
      for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
        continue;
      }
      lastsgp->next = sgp;
      break;
    }
  }
  if (sap == NULL) {
    if (bsp->annot != NULL) {
      for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
        continue;
      }
      sap->next = NewSeqAnnotType3 ("Phrap Graph", sgp);
    } else {
      bsp->annot = NewSeqAnnotType3 ("Phrap Graph", sgp);
    }
  }
}

static CharPtr BioseqGetLocalIdStr (
  BioseqPtr bsp
)

{
  ObjectIdPtr  oip;
  SeqIdPtr     sip;

  if (bsp == NULL) return NULL;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL && oip->str != NULL) {
        return oip->str;
      }
    }
  }
  return NULL;
}

typedef struct reqcontig {
  Int2  index;
  Char  str [41];
} ResqContig, PNTR ResqContigPtr;

#define MAX_FIELDS  8

static CharPtr ReadContigFile (
  CharPtr directory,
  CharPtr base,
  ValNodePtr PNTR fragmentgroupsp,
  CharPtr dumsp6,
  CharPtr dumt7,
  CharPtr PNTR sp6_clonep,
  CharPtr PNTR sp6_endp,
  CharPtr PNTR t7_clonep,
  CharPtr PNTR t7_endp
)

{
  Char        buf [256], instr [120];
  FileCache   fc;
  CharPtr     field [MAX_FIELDS];
  FILE        *fp;
  int         frg;
  Boolean     left_end, right_end, nonewline;
  Int4        len;
  Int2        numFields;
  CharPtr     pstring = NULL, ptr, str, sp6_end = NULL, t7_end = NULL;
  ValNodePtr  rescuedcontigs = NULL, vnp;

  fp = OpenOneFile (directory, base, ".ctg");
  if (fp == NULL) return NULL;

  FileCacheSetup (&fc, fp);

  str = FileCacheReadLine (&fc, buf, sizeof (buf), &nonewline);
  while (str != NULL) {
    MemSet ((Pointer) field, 0, sizeof (field));

   /*
   *  parse tab-delimited output line into array of fields, avoiding use of
   *  strtok so that empty columns (adjacent tabs) are properly assigned to
   *  field array
   */

    ptr = buf;
    for (numFields = 0; numFields < MAX_FIELDS && ptr != NULL; numFields++) {
      field [numFields] = ptr;
      ptr = StringChr (ptr, '\t');
      if (ptr == NULL) {
        ptr = StringChr (ptr, '\n');
      }
      if (ptr == NULL) {
        ptr = StringChr (ptr, '\r');
      }
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
    }

    if (StringDoesHaveText (field [0])) {
      StringNCpy_0 (instr, field [0], sizeof (instr) - 2);
      if (StringDoesHaveText (field [1])) {
        if (StringNICmp (field [1], "-", 1) == 0) {
          StringCat (instr, "-");
        }
      }
      ValNodeCopyStr (&rescuedcontigs, 0, instr);
      if (StringDoesHaveText (field [2])) {
        if (sscanf (field [2], "%d", &frg) == 1) {
          ValNodeCopyStr (fragmentgroupsp, (Uint1) frg, field [0]);
        }
      }
      left_end = FALSE;
      right_end = FALSE;
      if (StringDoesHaveText (field [3])) {
         if (StringDoesHaveText (field [4])) {
           if (StringNICmp (field [4], "l", 1) == 0) {
            left_end = TRUE;
           } else if (StringNICmp (field [4], "r", 1) == 0) {
             right_end = TRUE;
           }
         }
         if (StringICmp (field [3], "sp6") == 0) {
           StringCpy (dumsp6, field [0]);
           if (left_end) {
             StringCat (dumsp6, ",left");
           } else if (right_end) {
             StringCat (dumsp6, ",right");
           }
           if (sp6_clonep != NULL && *sp6_clonep == NULL) {
             *sp6_clonep = dumsp6;
           }
         } else if (StringICmp (field [3], "t7") == 0) {
           StringCpy (dumt7, field [0]);
           if (left_end) {
             StringCat (dumt7, ",left");
           } else if (right_end) {
             StringCat (dumt7, ",right");
           }
           if (t7_clonep != NULL && *t7_clonep == NULL) {
             *t7_clonep = dumt7;
           }
         }
      }
    }
    str = FileCacheReadLine (&fc, buf, sizeof (buf), &nonewline);
  }

  FileClose (fp);

  len = 0;
  for (vnp = rescuedcontigs; vnp != NULL; vnp = vnp->next) {
    len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
  }
  if (len > 1) {
    pstring = MemNew ((size_t) (len + 2));
    for (vnp = rescuedcontigs; vnp != NULL; vnp = vnp->next) {
      if (vnp != rescuedcontigs) {
        StringCat (pstring, ",");
      }
      StringCat (pstring, (CharPtr) vnp->data.ptrvalue);
    }
  }

   rescuedcontigs = ValNodeFreeData (rescuedcontigs);

  if (sp6_clonep != NULL && *sp6_clonep != NULL) {
    sp6_end = StringChr (*sp6_clonep, ',');
    if (sp6_end != NULL) {
      *sp6_end = '\0';
      sp6_end++;
      if (StringICmp (sp6_end, "left") == 0) {
        sp6_end = "left";
      } else if (StringICmp (sp6_end, "right") == 0) {
        sp6_end = "right";
      } else {
        sp6_end = NULL;
      }
    }
    if (sp6_endp != NULL) {
      *sp6_endp = sp6_end;
    }
  }
  if (t7_clonep != NULL && *t7_clonep != NULL) {
    t7_end = StringChr (*t7_clonep, ',');
    if (t7_end != NULL) {
      *t7_end = '\0';
      t7_end++;
      if (StringICmp (t7_end, "left") == 0) {
        t7_end = "left";
      } else if (StringICmp (t7_end, "right") == 0) {
        t7_end = "right";
      } else {
        t7_end = NULL;
      }
    }
    if (t7_endp != NULL) {
      *t7_endp = t7_end;
    }
  }

  return pstring;
}

static void MakeAssemblyFragments (
  BioseqPtr bsp,
  CharPtr name,
  Int2 index,
  CharPtr sp6_clone,
  CharPtr sp6_end,
  CharPtr t7_clone,
  CharPtr t7_end,
  Uint1 frag
)

{
  DeltaSeqPtr  dsp = NULL;
  Int4         from, to;
  ImpFeatPtr   ifp;
  SeqLitPtr    litp;
  SeqFeatPtr   sfp;
  SeqInt       sint;
  Char         str [128];
  Char         tmp [32];
  ValNode      vn;

  if (bsp == NULL || name == NULL || index < 1) return;
  from = 0;
  to = 0;
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
         dsp != NULL && index > 1; dsp = dsp->next, index--) {
      if (dsp->choice == 1) {
        from += SeqLocLen ((SeqLocPtr) dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          from += litp->length;
        }
      }
    }
  }
  if (dsp != NULL && dsp->choice == 2) {
    litp = (SeqLitPtr) dsp->data.ptrvalue;
    if (litp != NULL) {
      to = litp->length + from - 1;
    }
  }
  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  vn.choice = SEQLOC_INT;
  vn.data.ptrvalue = &sint;

  MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
  sint.id = SeqIdDup (SeqIdFindBest (bsp->id, 0));

  sint.from = from;
  sint.to = to;
  sint.strand = Seq_strand_plus;

  ifp = ImpFeatNew ();
  if (ifp == NULL) return;
  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, &vn);
  if (sfp == NULL) return;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");

  sprintf (str, "assembly_name:%s", name);
  if (frag > 0) {
    sprintf (tmp, "~fragment_group:%d", (int) frag);
    StringCat (str, tmp);
  }
  if (StringICmp (name, sp6_clone) == 0) {
    StringCat (str, "~clone_end:SP6");
    if (sp6_end != NULL) {
      StringCat (str, "~vector_side:");
      StringCat (str, sp6_end);
    }
  } else if (StringICmp (name, t7_clone) == 0) {
    StringCat (str, "~clone_end:T7");
    if (t7_end != NULL) {
      StringCat (str, "~vector_side:");
      StringCat (str, t7_end);
    }
  }
  sfp->comment = StringSaveNoNull (str);
}

static Uint2 ProcessPhrapAce (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  CharPtr localname,
  SeqEntryPtr gsep,
  MolInfoPtr template_molinfo,
  CharPtr directory,
  CharPtr base
)

{
  BioseqPtr      bsp, deltabsp;
  BioseqSetPtr   bssp;
  CharPtr        contigs;
  Boolean        do_contig = FALSE;
  Char           dumsp6 [64], dumt7 [64];
  Uint2          entityID;
  SeqEntryPtr    firstsep, nextsep, sep, topsep;
  Uint1          frag;
  IntFuzzPtr     ifp;
  Int2           index = 0;
  Boolean        is_unk100, lastwasraw;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  ResqContigPtr  rcp;
  ResqSeqgphPtr  rsp;
  CharPtr        seqbuf;
  SeqIdPtr       sip;
  SeqLitPtr      slp;
  CharPtr        sp6_clone = NULL, t7_clone = NULL, sp6_end = NULL, t7_end = NULL;
  ValNodePtr     rescuedcontigs = NULL, rescuedsgps = NULL, fragmentgroups = NULL, vnp, vnp2;

  if (fp == NULL) return 0;

  firstsep = ReadPhrapFile (fp);
  if (firstsep == NULL) return 0;

  dumsp6 [0] = '\0';
  dumt7 [0] = '\0';
  contigs = ReadContigFile (directory, base, &fragmentgroups, dumsp6,
                            dumt7, &sp6_clone, &sp6_end, &t7_clone, &t7_end);
  firstsep = SetPhrapContigOrder (firstsep, contigs);
  if (firstsep == NULL) return 0;
  if (contigs != NULL) {
    do_contig = TRUE;
  }

  /* always make delta, even if one component */

  bsp = FindNucBioseq (firstsep);
  if (bsp == NULL) return 0;

  sip = SeqIdSetDup (bsp->id);
  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);

  deltabsp = BioseqNew ();
  if (deltabsp == NULL) return 0;
  deltabsp->repr = Seq_repr_delta;
  deltabsp->seq_ext_type = 4;
  deltabsp->mol = Seq_mol_dna;
  deltabsp->length = 0;

  topsep = SeqEntryNew ();
  if (topsep == NULL) return 0;
  topsep->choice = 1;
  topsep->data.ptrvalue = (Pointer) deltabsp;

  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) deltabsp, topsep);

  if (gsep != NULL) {
    bssp = (BioseqSetPtr) gsep->data.ptrvalue;
    if (bssp == NULL) return 0;

    SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
    GetSeqEntryParent (gsep, &parentptr, &parenttype);

    bssp->seq_set = topsep;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) deltabsp, gsep);

    SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
  } else {
    entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) deltabsp);
  }

  lastwasraw = FALSE;
  for (sep = firstsep; sep != NULL; sep = nextsep) {
    nextsep = sep->next;
    sep->next = NULL;

    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) continue;

    if (bsp->repr == Seq_repr_raw) {

      if (lastwasraw) {
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) break;

        slp->length = 100 ;
        is_unk100 = TRUE;

        if (slp->length < 1 || is_unk100) {
          if (slp->length < 1) {
            slp->length = 0;
          }
          ifp = IntFuzzNew ();
          ifp->choice = 4;
          slp->fuzz = ifp;
        }

        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);

        deltabsp->length += slp->length;
        index++;
      }

      BioseqRawConvert (bsp, Seq_code_iupacna);
      seqbuf = BSMerge ((ByteStorePtr) bsp->seq_data, NULL);
      slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slp == NULL) continue;

      slp->length = bsp->length;
      ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
      slp->seq_data = (SeqDataPtr) BSNew (slp->length);
      slp->seq_data_type = Seq_code_iupacna;
      AddBasesToByteStore ((ByteStorePtr) slp->seq_data, seqbuf);
      MemFree (seqbuf);
      lastwasraw = TRUE;

      deltabsp->length += slp->length;
      index++;

      RescueSeqGraphs (bsp, index, &rescuedsgps);
      if (do_contig) {
        rcp = (ResqContigPtr) MemNew (sizeof (ResqContig));
        if (rcp != NULL) {
          rcp->index = index;
          StringNCpy_0 (rcp->str, BioseqGetLocalIdStr (bsp), sizeof (rcp->str));
          ValNodeAddPointer (&rescuedcontigs, 0, (Pointer) rcp);
        }
      }
    }

    SeqEntryFree (sep);
  }

  ValNodeLink (&(deltabsp->descr), vnp);
  deltabsp->id = sip;

  if (deltabsp != NULL) {
    for (vnp = rescuedsgps; vnp != NULL; vnp = vnp->next) {
      rsp = (ResqSeqgphPtr) vnp->data.ptrvalue;
      if (rsp != NULL) {
        OffsetAndLinkSeqGraph (deltabsp, rsp->sgp, (Int2) rsp->index);
      }
    }
    for (vnp = rescuedcontigs; vnp != NULL; vnp = vnp->next) {
      rcp = (ResqContigPtr) vnp->data.ptrvalue;
      if (rcp != NULL) {
        frag = 0;
        for (vnp2 = fragmentgroups; vnp2 != NULL; vnp2 = vnp2->next) {
          if (StringICmp ((CharPtr) vnp2->data.ptrvalue, rcp->str) == 0) {
            frag = (Uint1) vnp2->choice;
          }
        }
        MakeAssemblyFragments (deltabsp, rcp->str, (Int2) rcp->index,
                               sp6_clone, sp6_end, t7_clone, t7_end, frag);
      }
    }
  }
  rescuedsgps = ValNodeFreeData (rescuedsgps);
  rescuedcontigs = ValNodeFreeData (rescuedcontigs);


  if (gsep == NULL) {
    SeqMgrLinkSeqEntry (topsep, 0, NULL);
  }

  if (StringDoesHaveText (localname)) {
    sip = MakeSeqID (localname);
    if (sip != NULL) {
      bsp->id = SeqIdSetFree (bsp->id);
      bsp->id = sip;
      SeqMgrReplaceInBioseqIndex (bsp);
      VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
    }
  }

  ProcessOneNuc (entityID, deltabsp, src, tbl, template_molinfo);

  return entityID;
}

static Uint2 ProcessBulkSet (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  MolInfoPtr template_molinfo
)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Uint2         entityID;
  SeqEntryPtr   lastsep, sep, topsep;
  /*
  Pointer       dataptr;
  Uint2         datatype;
  */

  if (fp == NULL || tbl == NULL) return 0;

  bssp = BioseqSetNew ();
  if (bssp == NULL) return 0;

  switch (tbl->whichclass) {
    case 1 :
      bssp->_class = BioseqseqSet_class_pop_set;
      break;
    case 2 :
      bssp->_class = BioseqseqSet_class_phy_set;
      break;
    case 3 :
      bssp->_class = BioseqseqSet_class_mut_set;
      break;
    case 4 :
      bssp->_class = BioseqseqSet_class_eco_set;
      break;
    default :
      bssp->_class = BioseqseqSet_class_genbank;
      break;
  }

  topsep = SeqEntryNew ();
  if (topsep == NULL) return 0;
  topsep->choice = 2;
  topsep->data.ptrvalue = (Pointer) bssp;

  entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);

  lastsep = NULL;

/*
  while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, TRUE, FALSE, TRUE, FALSE)) != NULL) {
    if (datatype == OBJ_BIOSEQ) {

      sep = SeqMgrGetSeqEntryForData (dataptr);
      if (lastsep == NULL) {
        bssp->seq_set = sep;
      } else {
        lastsep->next = sep;
      }
      lastsep = sep;

      bsp = (BioseqPtr) dataptr;
      ProcessOneNuc (entityID, bsp, src, tbl);

    } else {
      ObjMgrFree (datatype, dataptr);
    }
  }
*/

  while ((bsp = ReadDeltaFasta (fp, NULL)) != NULL) {

    sep = SeqMgrGetSeqEntryForData (bsp);
    if (lastsep == NULL) {
      bssp->seq_set = sep;
    } else {
      lastsep->next = sep;
    }
    lastsep = sep;

    ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);
  }

  SeqMgrLinkSeqEntry (topsep, 0, NULL);

  return entityID;
}

static SeqEntryPtr FA2SEP (
  FILE *fp
)

{
  BioseqPtr    bsp;
  Pointer      dataptr;
  Uint2        datatype;
  SeqEntryPtr  sep;

  if (fp == NULL) return NULL;

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, TRUE, FALSE, TRUE, FALSE);
  if (datatype == OBJ_BIOSEQ) {
    sep = SeqMgrGetSeqEntryForData (dataptr);
    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        bsp = (BioseqPtr) dataptr;
        sep->choice = 1;
        sep->data.ptrvalue = bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
    }
    return sep;
  }

  return NULL;
}

static SeqEntryPtr MakeUnk100GapSep (void)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;

  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;
  bsp->repr = Seq_repr_virtual;
  bsp->mol = Seq_mol_na;
  bsp->length = 100;
  bsp->id = SeqIdParse ("lcl|unk100");
  sep->choice = 1;
  sep->data.ptrvalue = bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
  return sep;
}

static Uint2 ProcessDeltaSet (
  FILE* fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  CharPtr localname,
  SeqEntryPtr gsep,
  MolInfoPtr template_molinfo
)

{
  BioseqPtr      bsp, deltabsp;
  BioseqSetPtr   bssp;
  Uint2          entityID;
  SeqEntryPtr    firstsep, lastsep, nextsep, sep, tmp, topsep;
  IntFuzzPtr     ifp;
  Boolean        is_unk100;
  ObjectIdPtr    oip;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  CharPtr        seqbuf;
  SeqIdPtr       sip, virtid;
  SeqLitPtr      slp;
  ValNodePtr     vnp;

  if (fp == NULL) return 0;

  firstsep = NULL;
  lastsep = NULL;

  /*
  sep = FastaToSeqEntry (fp, TRUE);
  */
  sep = FA2SEP (fp);
  if (sep == NULL) return 0;

  /* loop to collect subsequent entries */

  while (sep != NULL) {
    if (firstsep == NULL) {
      firstsep = sep;
    }
    if (tbl->implicitgaps && lastsep != NULL) {
      tmp = MakeUnk100GapSep ();
      if (tmp != NULL) {
        ValNodeLink (&lastsep, tmp);
        lastsep = tmp;
      }
    }
    if (lastsep != NULL) {
      ValNodeLink (&lastsep, sep);
    }
    lastsep = sep;
    /*
    sep = FastaToSeqEntry (fp, TRUE);
    */
    sep = FA2SEP (fp);
  }

  /* if only one FASTA, treat as raw */

  if (firstsep->next == NULL) {
    bsp = FindNucBioseq (firstsep);
    if (bsp == NULL) return 0;

    if (gsep != NULL) {
      bssp = (BioseqSetPtr) gsep->data.ptrvalue;
      if (bssp == NULL) return 0;

      SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
      GetSeqEntryParent (gsep, &parentptr, &parenttype);

      bssp->seq_set = SeqMgrGetSeqEntryForData (bsp);
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, gsep);

      SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
      RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

      entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
    } else {
      entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
    }

    ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);
    return entityID;
  }

  /* now process delta */

  bsp = FindNucBioseq (firstsep);
  if (bsp == NULL) return 0;

  sip = SeqIdSetDup (bsp->id);
  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);

  deltabsp = BioseqNew ();
  if (deltabsp == NULL) return 0;
  deltabsp->repr = Seq_repr_delta;
  deltabsp->seq_ext_type = 4;
  deltabsp->mol = Seq_mol_dna;
  deltabsp->length = 0;

  topsep = SeqEntryNew ();
  if (topsep == NULL) return 0;
  topsep->choice = 1;
  topsep->data.ptrvalue = (Pointer) deltabsp;

  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) deltabsp, topsep);

  if (gsep != NULL) {
    bssp = (BioseqSetPtr) gsep->data.ptrvalue;
    if (bssp == NULL) return 0;

    SaveSeqEntryObjMgrData (gsep, &omdptop, &omdata);
    GetSeqEntryParent (gsep, &parentptr, &parenttype);

    bssp->seq_set = topsep;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) deltabsp, gsep);

    SeqMgrLinkSeqEntry (gsep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (gsep, omdptop, &omdata);

    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
  } else {
    entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) deltabsp);
  }

  for (sep = firstsep; sep != NULL; sep = nextsep) {
    nextsep = sep->next;
    sep->next = NULL;

    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) continue;

    if (bsp->repr == Seq_repr_raw) {
      BioseqRawConvert (bsp, Seq_code_iupacna);
      seqbuf = BSMerge ((ByteStorePtr) bsp->seq_data, NULL);
      slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slp == NULL) continue;

      slp->length = bsp->length;
      ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
      slp->seq_data = (SeqDataPtr) BSNew (slp->length);
      slp->seq_data_type = Seq_code_iupacna;
      AddBasesToByteStore ((ByteStorePtr) slp->seq_data, seqbuf);
      MemFree(seqbuf);

      deltabsp->length += slp->length;

    } else if (bsp->repr == Seq_repr_virtual) {
      slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slp == NULL) continue;

      slp->length = bsp->length;
      ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);

      is_unk100 = FALSE;
      virtid = bsp->id;
      if (virtid != NULL && virtid->choice == SEQID_LOCAL) {
        oip = (ObjectIdPtr) virtid->data.ptrvalue;
        if (oip != NULL) {
          if (StringCmp (oip->str, "unk100") == 0) {
            is_unk100 = TRUE;
          }
        }
      }
      if (slp->length < 1 || is_unk100) {
        if (slp->length < 1) {
          slp->length = 0;
        }
        ifp = IntFuzzNew ();
        ifp->choice = 4;
        slp->fuzz = ifp;
      }

      deltabsp->length += slp->length;
    }

    SeqEntryFree (sep);
  }

  ValNodeLink (&(deltabsp->descr), vnp);
  deltabsp->id = sip;

  if (gsep == NULL) {
    SeqMgrLinkSeqEntry (topsep, 0, NULL);
  }

  if (StringDoesHaveText (localname)) {
    sip = MakeSeqID (localname);
    if (sip != NULL) {
      bsp->id = SeqIdSetFree (bsp->id);
      bsp->id = sip;
      SeqMgrReplaceInBioseqIndex (bsp);
      VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
    }
  }

  ProcessOneNuc (entityID, deltabsp, src, tbl, template_molinfo);

  return entityID;
}

static Boolean DoSequenceLengthsMatch (
  TAlignmentFilePtr afp
)

{
  int    seq_index;
  Int4   seq_len;

  if (afp == NULL || afp->sequences == NULL || afp->num_sequences == 0) {
    return TRUE;
  }
  seq_len = StringLen (afp->sequences[0]);
  for (seq_index = 1; seq_index < afp->num_sequences; seq_index++) {
    if (StringLen (afp->sequences[seq_index]) != seq_len) {
      return FALSE;
    }
  }
  return TRUE;
}

static void ShowAlignmentNotes (
  TAlignmentFilePtr afp,
  TErrorInfoPtr error_list
)

{
  TErrorInfoPtr eip;
  Int4         index;

  for (eip = error_list; eip != NULL; eip = eip->next) {
    printf ("*****\nError category %d\n", eip->category);
    if (eip->line_num > -1) {
      printf ("Line number %d\n", eip->line_num);
    }
    if (eip->id != NULL) {
      printf ("Sequence ID %s\n", eip->id);
    }
    if (eip->message != NULL) {
      printf ("%s\n", eip->message);
    }
  }
  if (afp == NULL) {
    printf ("Catastrophic failure during reading\n");
  } else {
    printf ("Found %d sequences\n", afp->num_sequences);
    printf ("Found %d organisms\n", afp->num_organisms);
    for (index = 0; index < afp->num_sequences; index++)
    {
      printf ("\t%s\t", afp->ids [index]);
      if (index < afp->num_organisms) {
        printf ("%s\n", afp->organisms [index]);
      } else {
        printf ("No organism information\n");
      }
    }
    while (index < afp->num_organisms) {
      printf ("Unclaimed organism: %s\n", afp->organisms [index]);
      index++;
    }
  }
}

static Uint2 ProcessAlignSet (
  FILE *fp,
  BioSourcePtr src,
  TblArgsPtr tbl,
  MolInfoPtr template_molinfo
)

{
  TSequenceInfoPtr  sequence_info;
  TErrorInfoPtr     error_list;
  ReadBufferData    rbd;
  TAlignmentFilePtr afp;
  SeqEntryPtr       sep = NULL;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Char              ch;
  Uint2             entityID;
  SeqEntryPtr       tmp;
  Char              nucleotide_alphabet[] = "ABCDGHKMRSTUVWXYabcdghkmrstuvwxy";
  Char              protein_alphabet[] = "ABCDEFGHIKLMPQRSTUVWXYZabcdefghiklmpqrstuvwxyz";
  Uint1             moltype = Seq_mol_dna;

  if (fp == NULL) return 0;

  sequence_info = SequenceInfoNew ();
  if (sequence_info == NULL) return 0;

  /* format sequence options based on commandline arguments */
  /* set sequence alphabet */
  if (tbl->aln_is_protein) {
    moltype = Seq_mol_aa;
    sequence_info->alphabet = protein_alphabet;
  } else {
    moltype = Seq_mol_dna;
    sequence_info->alphabet = nucleotide_alphabet;
  }

  sequence_info->beginning_gap = MemFree (sequence_info->beginning_gap);
  if (StringHasNoText (tbl->aln_beginning_gap)) {
    sequence_info->beginning_gap = StringSave (".-?");
  } else {
    sequence_info->beginning_gap = StringSave (tbl->aln_beginning_gap);
  }
  sequence_info->middle_gap = MemFree (sequence_info->middle_gap);
  if (StringHasNoText (tbl->aln_middle_gap)) {
    sequence_info->middle_gap = StringSave ("-");
  } else {
    sequence_info->middle_gap = StringSave (tbl->aln_middle_gap);
  }
  sequence_info->end_gap = MemFree (sequence_info->end_gap);
  if (StringHasNoText (tbl->aln_end_gap)) {
    sequence_info->end_gap = StringSave (".-?");
  } else {
    sequence_info->end_gap = StringSave (tbl->aln_end_gap);
  }
  sequence_info->missing = MemFree (sequence_info->missing);
  if (StringHasNoText (tbl->aln_missing)) {
    sequence_info->missing = StringSave ("Nn?");
  } else {
    sequence_info->missing = StringSave (tbl->aln_missing);
  }
  sequence_info->match = MemFree (sequence_info->match);
  if (StringHasNoText (tbl->aln_match)) {
    sequence_info->match = StringSave (".");
  } else {
    sequence_info->match = StringSave (tbl->aln_match);
  }

  error_list = NULL;
  rbd.fp = fp;
  rbd.current_data = NULL;
  afp = ReadAlignmentFile ( AbstractReadFunction,
                            (Pointer) &rbd,
                            AbstractReportError,
                            (Pointer) &error_list,
                            sequence_info);

  ShowAlignmentNotes (afp, error_list);
  ErrorInfoFree (error_list);
  if (afp != NULL) {
    if (afp->num_organisms == 0 && src == NULL) {
      printf ("No organisms supplied!\n");
    } else if (afp->num_organisms != 0 && afp->num_organisms != afp->num_sequences) {
      printf ( "Number of organisms must match number of sequences!");
    } else {
      ch = 'y';
      if (! DoSequenceLengthsMatch (afp)) {
        printf ("Sequences are not all the same length - are you sure you want to continue?");
        ch = getchar ();
      }
      if (ch == 'y' || ch == 'Y') {
        sep = MakeSequinDataFromAlignment (afp, moltype);
      }
    }
  }
  SequenceInfoFree (sequence_info);

  AlignmentFileFree (afp);

  if (sep == NULL || sep->data.ptrvalue == NULL) return 0;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
      entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
      ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    bssp->_class = BioseqseqSet_class_phy_set;
    entityID = ObjMgrRegister (OBJ_BIOSEQSET, (Pointer) bssp);
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      if (IS_Bioseq (tmp)) {
        bsp = (BioseqPtr) tmp->data.ptrvalue;
        ProcessOneNuc (entityID, bsp, src, tbl, template_molinfo);
      }
    }
  } else return 0;

  SeqMgrLinkSeqEntry (sep, 0, NULL);

  return entityID;
}

static SeqAnnotPtr NewGraphSeqAnnot (
  CharPtr name,
  SeqGraphPtr sgp
)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (StringDoesHaveText (name)) {
    SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

typedef struct npsseqs {
  BioseqPtr  nuc;
  BioseqPtr  prot;
} NpsSeqs, PNTR NpsSeqsPtr;

static void FindNucProtSeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  NpsSeqsPtr  nsp;

  if (bsp == NULL) return;
  nsp = (NpsSeqsPtr) userdata;
  if (nsp == NULL) return;

  if (ISA_na (bsp->mol)) {
    nsp->nuc = bsp;
  } else if (ISA_aa (bsp->mol)) {
    nsp->prot = bsp;
  }
}

static Boolean InRightNps (
  CharPtr gbqval,
  SeqIdPtr protids,
  Boolean force_local_id
)

{
  Int2      adv;
  Char      id [64];
  Char      lcl [64];
  SeqIdPtr  sip = NULL;
  CharPtr   ptr;
  Boolean   rsult;
  long int  val;
  Uint4     version = 0;

  StringNCpy_0 (id, gbqval, sizeof (id));
  if (StringDoesHaveText (id)) {
    if (StringChr (id, '|') != NULL) {
      sip = SeqIdParse (id);
    } else if (force_local_id) {
      sprintf (lcl, "lcl|%s", id);
      sip = SeqIdParse (lcl);
    } else {
      adv = ValidateAccnDotVer (id);
      if (adv == 0 || adv == -5) {
        ptr = StringChr (id, '.');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          if (sscanf (ptr, "%ld", &val) == 1) {
            version = (Uint4) val;
          }
        }
        sip = SeqIdFromAccession (id, version, NULL);
      } else {
        sprintf (lcl, "lcl|%s", id);
        sip = SeqIdParse (lcl);
      }
    }
  }
  if (sip == NULL) return FALSE;
  rsult = SeqIdIn (sip, protids);
  SeqIdFree (sip);
  return rsult;
}

static void MakeNucProtCDS (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CodeBreakPtr    cbp;
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  GBQualPtr       gbq;
  Char            id [64];
  SeqFeatPtr      mrna;
  GBQualPtr       nextqual;
  NpsSeqs         ns;
  Boolean         partial5, partial3;
  GBQualPtr PNTR  prevqual;
  SeqFeatPtr      sfp;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Int4            start, stop;
  TblArgsPtr      tbl;
  SeqFeatPtr      temp;

  tbl = (TblArgsPtr) userdata;
  if (tbl == NULL) return;

  ns.nuc = NULL;
  ns.prot = NULL;
  if (VisitBioseqsInSet (bssp, (Pointer) &ns, FindNucProtSeqs) != 2) return;
  if (ns.nuc == NULL || ns.prot == NULL) return;

  cds = SeqMgrGetCDSgivenProduct (ns.prot, NULL);
  mrna = SeqMgrGetRNAgivenProduct (ns.nuc, NULL);
  if (cds == NULL || mrna == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  start = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_START);
  stop = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_STOP);

  if (start < 0 || start >= ns.nuc->length ||
      stop < 0 || stop >= ns.nuc->length) return;

  sip = SeqIdFindBest (ns.nuc->id, 0);
  if (sip == NULL) return;

  /* copy cds feature fields to paste into new cds feature */
  temp = AsnIoMemCopy (cds,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  sfp = CreateNewFeatureOnBioseq (ns.nuc, SEQFEAT_CDREGION, NULL);
  if (sfp == NULL) return;

  sfp->location = SeqLocFree (sfp->location);
  if (StringISearch (cds->except_text, "ribosomal slippage") == NULL &&
      StringISearch (cds->except_text, "ribosome slippage") == NULL &&
      StringISearch (cds->except_text, "trans splicing") == NULL &&
      StringISearch (cds->except_text, "trans-splicing") == NULL &&
      StringISearch (cds->except_text, "artificial frameshift") == NULL) {
    sfp->location = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);
  } else {
    slp = SeqLocFindNext (cds->location, NULL);
    while (slp != NULL) {
      start = GetOffsetInLoc (slp, mrna->location, SEQLOC_START);
      stop = GetOffsetInLoc (slp, mrna->location, SEQLOC_STOP);
      sfp->location = AddIntervalToLocation (sfp->location, sip, start, stop, partial5, partial3);
      slp = SeqLocFindNext (cds->location, slp);
    }
    sfp->location = SeqLocMergeEx (ns.nuc, sfp->location, NULL, FALSE, TRUE, FALSE, FALSE);
  }
  SetSeqFeatProduct (sfp, ns.prot);

  /* paste fields from temp copy of original cds */
  crp = (CdRegionPtr) temp->data.value.ptrvalue;
  sfp->data.value.ptrvalue = (Pointer) crp;

  sfp->partial = temp->partial;
  sfp->excpt = temp->excpt;
  sfp->comment = temp->comment;
  sfp->qual = temp->qual;
  sfp->title = temp->title;
  sfp->ext = temp->ext;
  sfp->cit = temp->cit;
  sfp->exp_ev = temp->exp_ev;
  sfp->xref = temp->xref;
  sfp->dbxref = temp->dbxref;
  sfp->pseudo = temp->pseudo;
  sfp->except_text = temp->except_text;

  MemFree (temp); /* do not SeqFeatFree */

  /* update code break locations */
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
    CheckSeqLocForPartial (cbp->loc, &partial5, &partial3);
    start = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_START);
    stop = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_STOP);
    if (start < 0 || start >= ns.nuc->length ||
        stop < 0 || stop >= ns.nuc->length) continue;
    cbp->loc = SeqLocFree (cbp->loc);
    cbp->loc = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);;
  }

  /* get rid of protein_id in mRNA if it matches protein Seq-id */
  gbq = mrna->qual;
  prevqual = (GBQualPtr PNTR) &(mrna->qual);
  id [0] = '\0';
  sip = NULL;
  while (gbq != NULL) {
    nextqual = gbq->next;
    if (StringICmp (gbq->qual, "protein_id") == 0 &&
        InRightNps (gbq->val, ns.prot->id, tbl->forcelocalid)) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      StringNCpy_0 (id, gbq->val, sizeof (id));
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}

/* copy gene from contig onto nuc-prot, single interval on cdna bioseq */

static void CopyGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene, copy, temp;
  GeneRefPtr         grp, xref;
  Boolean            partial5, partial3;

  /* input mrna features are multi-interval on contig */

  if (sfp->data.choice != SEQFEAT_RNA) return;

  /* find cdna product of mrna */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  /* check for gene xref */

  xref = SeqMgrGetGeneXref (sfp);
  if (xref != NULL) {
    if (SeqMgrGeneIsSuppressed (xref)) return;

    /* copy gene xref for new gene feature */

    grp = AsnIoMemCopy (xref,
                        (AsnReadFunc) GeneRefAsnRead,
                        (AsnWriteFunc) GeneRefAsnWrite);
    if (grp == NULL) return;

    /* make new gene feature on full-length of cdna */

    copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
    if (copy == NULL) return;

    copy->data.value.ptrvalue = grp;
    return;
  }

  /* overlapping gene should be single interval on contig */

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;

  CheckSeqLocForPartial (gene->location, &partial5, &partial3);

  /* copy gene feature fields to paste into new gene feature */

  temp = AsnIoMemCopy (gene,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* make new gene feature on full-length of cdna */

  copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
  if (copy == NULL) {
    SeqFeatFree (temp);
    return;
  }

  /* paste fields from temp copy of original gene */

  copy->data.value.ptrvalue = temp->data.value.ptrvalue;
  copy->partial = temp->partial;
  copy->excpt = temp->excpt;
  copy->comment = temp->comment;
  copy->qual = temp->qual;
  copy->title = temp->title;
  copy->ext = temp->ext;
  copy->cit = temp->cit;
  copy->exp_ev = temp->exp_ev;
  copy->xref = temp->xref;
  copy->dbxref = temp->dbxref;
  copy->pseudo = temp->pseudo;
  copy->except_text = temp->except_text;

  SetSeqLocPartial (copy->location, partial5, partial3);

  SeqLocFree (temp->location);
  MemFree (temp); /* do not SeqFeatFree */
}

static void CopyNcRna (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr   bsp;
  SeqFeatPtr  copy, temp;
  Boolean     partial5, partial3;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->idx.subtype != FEATDEF_ncRNA) return;

  /* find instantiated product of ncRNA */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  /* copy ncRNA feature fields to paste into new ncRNA feature */

  temp = AsnIoMemCopy (sfp,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* make new ncRNA feature on full-length of transcript */

  copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_RNA, NULL);
  if (copy == NULL) {
    SeqFeatFree (temp);
    return;
  }

  /* paste fields from temp copy of original ncRNA */

  copy->data.value.ptrvalue = temp->data.value.ptrvalue;
  copy->partial = temp->partial;
  copy->excpt = temp->excpt;
  copy->comment = temp->comment;
  copy->qual = temp->qual;
  copy->title = temp->title;
  copy->ext = temp->ext;
  copy->cit = temp->cit;
  copy->exp_ev = temp->exp_ev;
  copy->xref = temp->xref;
  copy->dbxref = temp->dbxref;
  copy->pseudo = temp->pseudo;
  copy->except_text = temp->except_text;

  SetSeqLocPartial (copy->location, partial5, partial3);

  SeqLocFree (temp->location);
  SeqLocFree (temp->product);
  MemFree (temp); /* do not SeqFeatFree */
}

static void ClearRnaProducts (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->product == NULL) return;

  sfp->product = SeqLocFree (sfp->product);
}

static void RemoveGBQualIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GBQualPtr       gbq;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;

  if (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA) return;

  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    nextqual = gbq->next;
    if (StringICmp (gbq->qual, "transcript_id") == 0 ||
        StringICmp (gbq->qual, "protein_id") == 0) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}

typedef struct dupprot {
  SeqFeatPtr  firstprot;
  SeqFeatPtr  secondprot;
} DupProt, PNTR DupProtPtr;

static void FindDupProtFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  DupProtPtr  dpp;
  ProtRefPtr  prp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  dpp = (DupProtPtr) userdata;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (dpp == NULL || prp == NULL) return;
  if (prp->processed != 0) return;
  if (dpp->firstprot == NULL) {
    dpp->firstprot = sfp;
  } else if (dpp->secondprot == NULL) {
    dpp->secondprot = sfp;
  }
}

static void ClearProtFeatStrand (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_REGION &&
      sfp->data.choice != SEQFEAT_SITE &&
      sfp->data.choice != SEQFEAT_BOND &&
      sfp->data.choice != SEQFEAT_PROT) return;

  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL) {
        if (sintp->strand != Seq_strand_unknown) {
          sintp->strand = Seq_strand_unknown;
        }
      }
    }
    slp = SeqLocFindNext (sfp->location, slp);
  }
}

static void RemoveDupProtFeats (
  BioseqPtr bsp,
  Pointer userdata
)

{
  DupProt  dp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  VisitFeaturesOnBsp (bsp, NULL, ClearProtFeatStrand);
  dp.firstprot = NULL;
  dp.secondprot = NULL;
  VisitFeaturesOnBsp (bsp, (Pointer) &dp, FindDupProtFeats);
  if (dp.firstprot == NULL || dp.secondprot == NULL) return;
  if (AsnIoMemComp ((Pointer) dp.firstprot, (Pointer) dp.secondprot, (AsnWriteFunc) SeqFeatAsnWrite)) {
    dp.firstprot->idx.deleteme = TRUE;
  }
}

/*
static void RemoveUnnecGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatXrefPtr  curr, next;
  SeqFeatXrefPtr  PNTR last;
  GeneRefPtr      grp, grpx;
  Boolean         redundantgenexref;
  SeqFeatPtr      sfpx;
  CharPtr         syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL) return;

  redundantgenexref = FALSE;
  if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if ((StringICmp (grp->locus, grpx->locus) == 0)) {
      redundantgenexref = TRUE;
    }
  } else if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if ((StringICmp (grp->locus_tag, grpx->locus_tag) == 0)) {
      redundantgenexref = TRUE;
    }
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) == 0) {
        redundantgenexref = TRUE;
      }
    }
  }

  if (redundantgenexref) {
    last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
    curr = sfp->xref;
    while (curr != NULL) {
      next = curr->next;
      if (curr->data.choice == SEQFEAT_GENE) {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      } else {
        last = &(curr->next);
      }
      curr = next;
    }
  }
}
*/

typedef struct dummysmfedata {
  Int4  max;
  Int4  num_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK T2ADummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  Int4          len;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
  }

  return TRUE;
}

static void FillInPartialGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  context;
  SeqFeatPtr         gene;
  GeneRefPtr         grp, grpx;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  if (StringDoesHaveText (grp->locus) || StringHasNoText (grp->locus_tag)) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &context);
  if (gene == NULL || gene->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) gene->data.value.ptrvalue;
  if (grpx == NULL) return;

  if (StringHasNoText (grpx->locus)) return;
  grp->locus = StringSave (grpx->locus);
}

static void RemoveUnnecGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  DummySmfeData        dsd;
  SeqMgrFeatContext    fcontext;
  SeqFeatXrefPtr PNTR  last;
  GeneRefPtr           grp, grpx;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL) return;

  if ((!StringHasNoText (grp->locus)) && (!StringHasNoText (grpx->locus))) {
    if ((StringICmp (grp->locus, grpx->locus) != 0)) return;
  } else if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grp->locus_tag)) {
    if ((StringICmp (grp->locus_tag, grpx->locus_tag) != 0)) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if ((!StringHasNoText (syn1)) && (!StringHasNoText (syn2))) {
      if ((StringICmp (syn1, syn2) != 0)) return;
    }
  }

  MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
  dsd.max = INT4_MAX;
  dsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                           LOCATION_SUBSET, (Pointer) &dsd, T2ADummySMFEProc);

  if (dsd.num_at_max < 2) {
    last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
    curr = sfp->xref;
    while (curr != NULL) {
      next = curr->next;
      if (curr->data.choice == SEQFEAT_GENE) {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      } else {
        last = &(curr->next);
      }
      curr = next;
    }
  }
}

static CharPtr RnaTypeLabel (
  SeqFeatPtr rna
)

{
  if (rna == NULL) return "RNA";
  switch (rna->idx.subtype) {
    case FEATDEF_preRNA :
      return "preRNA";
    case FEATDEF_mRNA :
      return "mRNA";
    case FEATDEF_tRNA :
      return "tRNA";
    case FEATDEF_rRNA :
      return "rRNA";
    case FEATDEF_snRNA :
      return "snRNA";
    case FEATDEF_scRNA :
      return "scRNA";
    case FEATDEF_otherRNA :
      return "otherRNA";
    case FEATDEF_snoRNA :
      return "snoRNA";
    case FEATDEF_ncRNA :
      return "ncRNA";
    case FEATDEF_tmRNA :
      return "tmRNA";
    default :
      break;
  }
  return "RNA";
}

static void AddRnaTitles (
  SeqFeatPtr rna,
  CharPtr organism
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  CharPtr            cdslabel = NULL;
  SeqMgrFeatContext  gcontext;
  CharPtr            genelabel = NULL;
  size_t             len;
  SeqFeatPtr         sfp;
  CharPtr            str;
  CharPtr            typ = NULL;

  if (rna == NULL || rna->product == NULL) return;
  bsp = BioseqFindFromSeqLoc (rna->product);
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  if (BioseqGetTitle (bsp) != NULL) return;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (sfp != NULL) {
    genelabel = gcontext.label;
    if (StringHasNoText (genelabel)) {
      genelabel = NULL;
    }
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  if (sfp != NULL) {
    cdslabel = ccontext.label;
    if (StringHasNoText (cdslabel)) {
      cdslabel = NULL;
    }
  }
  typ = RnaTypeLabel (rna); 
  len = StringLen (organism) + StringLen (genelabel) + StringLen (cdslabel) +
        StringLen (" mRNA, complete cds.") + StringLen (typ) + 10;
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str == NULL) return;
  str [0] = '\0';

  if (StringDoesHaveText (organism)) {
    StringCat (str, organism);
  }
  if (cdslabel != NULL) {
    StringCat (str, " ");
    StringCat (str, cdslabel);
  }
  if (genelabel != NULL) {
      StringCat (str, " (");
      StringCat (str, genelabel);
      StringCat (str, ")");
  }
  if (cdslabel != NULL && genelabel != NULL) {
    StringCat (str, " ");
    StringCat (str, typ);
    if (ccontext.partialL || ccontext.partialR) {
      StringCat (str, ", partial cds.");
    } else {
      StringCat (str, ", complete cds.");
    }
  } else if (genelabel != NULL) {
    StringCat (str, " ");
    StringCat (str, typ);
    StringCat (str, ".");
  }
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
}

static void MakeOneRnaTitle (
  SeqFeatPtr rna,
  SeqFeatPtr gene,
  CharPtr label,
  CharPtr organism,
  Boolean alt_splice
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  GeneRefPtr         grp;
  Char               id [64];
  CharPtr            lbl = NULL;
  size_t             len;
  CharPtr            ptr;
  CharPtr            str;
  CharPtr            typ = NULL;

  if (rna == NULL || rna->product == NULL) return;

  grp = SeqMgrGetGeneXref (rna);
  if (SeqMgrGeneIsSuppressed (grp)) return;
  if (grp == NULL && gene != NULL) {
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
  }
  if (grp == NULL) return;

  bsp = BioseqFindFromSeqLoc (rna->product);
  if (bsp == NULL) return;
  SeqIdWrite (bsp->id, id, PRINTID_TEXTID_ACC_VER, sizeof (id) - 1);

  typ = RnaTypeLabel (rna); 
  lbl = StringSaveNoNull (label);

  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);

  len = StringLen (organism) + StringLen (grp->locus_tag) + StringLen (grp->locus) +
        StringLen (id) + StringLen (" transcript variant") + StringLen (lbl) +
        StringLen (" mRNA, complete cds.") + StringLen (typ) + 20;
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str == NULL) return;
  str [0] = '\0';

  if (StringDoesHaveText (organism)) {
    StringCat (str, organism);
  }
  if (lbl != NULL) {
    StringCat (str, " ");
    ptr = StringStr (lbl, ", transcript variant ");
    if (ptr != NULL) {
      *ptr = '\0';
      ptr += 2;
      StringCat (str, lbl);
      if (StringDoesHaveText (grp->locus)) {
          StringCat (str, " (");
          StringCat (str, grp->locus);
          StringCat (str, ")");
      }
      StringCat (str, ", ");
      StringCat (str, ptr);
    } else {
      StringCat (str, lbl);
      if (StringDoesHaveText (grp->locus)) {
          StringCat (str, " (");
          StringCat (str, grp->locus);
          StringCat (str, ")");
      }
    }
  }

  StringCat (str, ", ");
  StringCat (str, typ);
  StringCat (str, ".");

  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
  MemFree (lbl);
}

static void MakeSmartRnaTitles (
  BioseqPtr bsp,
  CharPtr organism
)

{
  SeqMgrFeatContext  context;
  GmcDataPtr         gdp, head;
  GeneRefPtr         grp;
  Int2               i, j, k, numgene, numrna;
  SeqFeatPtr         sfp;

  if (bsp == NULL) return;

  numgene = 0;
  numrna = 0;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  while (sfp != NULL) {
    switch (sfp->data.choice) {
      case SEQFEAT_GENE :
        numgene++;
        break;
      case SEQFEAT_RNA :
        numrna++;
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }

  /* if (numgene == 0) return; */

  if (numrna > 0) {
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numrna + 1));
    if (head != NULL) {
      gdp = head;
      sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &context);
      while (sfp != NULL) {
        if (sfp->product != NULL) {
          gdp->feat = sfp;
          gdp->label = context.label;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          }
          gdp++;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, 0, &context);
      }
      HeapSort (head, (size_t) numrna, sizeof (GmcData), SortByGenePtr);
      for (i = 0; i < numrna; i += j) {
        sfp = head [i].gene;
        for (j = 1; i + j < numrna && sfp == head [i + j].gene; j++) continue;
        if (j == 1) {
          /* no alt splicing */
          MakeOneRnaTitle (head [i].feat, head [i].gene, head [i].label, organism, FALSE);
        } else {
          /* is alt splicing */
          for (k = 0; k < j; k++) {
            MakeOneRnaTitle (head [i + k].feat, head [i + k].gene, head [i + k].label, organism, TRUE);
          }
        }
      }
    }
    MemFree (head);
  }
}

typedef struct gosearch {
  TextFsaPtr  gotags;
  Boolean     isbad;
} GoSearch, PNTR GoSearchPtr;

static void LookForGo (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Char         ch;
  GoSearchPtr  gsp;
  CharPtr      ptr;
  Int4         state;
  ValNodePtr   matches;

  if (sfp == NULL || StringHasNoText (sfp->comment)) return;
  gsp = (GoSearchPtr) userdata;

  state = 0;
  ptr = sfp->comment;
  ch = *ptr;
  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (gsp->gotags, state, ch, &matches);
    if (matches != NULL) {
      gsp->isbad = TRUE;
    }
    ptr++;
    ch = *ptr;
  }
}

static Boolean HasGoTermsInNote (
  SeqEntryPtr sep,
  TextFsaPtr gotags
)

{
  GoSearch  gs;

  gs.gotags = gotags;
  gs.isbad = FALSE;
  VisitFeaturesInSep (sep, (Pointer) &gs, LookForGo);
  return gs.isbad;
}

static void TakeProteinsFromGPS (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqEntryPtr PNTR  lastp;
  SeqEntryPtr       sep;

  if (bsp == NULL || (! ISA_aa (bsp->mol))) return;
  lastp = (SeqEntryPtr PNTR) userdata;
  if (lastp == NULL) return;

  /* link copy after genomic sequence */

  bsp = (BioseqPtr) AsnIoMemCopy ((Pointer) bsp,
                                  (AsnReadFunc) BioseqAsnRead,
                                  (AsnWriteFunc) BioseqAsnWrite);
  sep = ValNodeAddPointer (lastp, 1, (Pointer) bsp);
  *lastp = sep;
}

static void GPStoNPS (
  SeqEntryPtr top,
  Uint2 entityID
)

{
  BioseqSetPtr  bssp;
  BioseqSetPtr  dum;
  SeqEntryPtr   last, sep;
  Uint2         parenttype;
  Pointer       parentptr;

  if (top == NULL || top->choice != 2) {
    Message (MSG_POSTERR, "GPStoNPS failed at top || top->choice");
    return;
  }
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_gen_prod_set) {
    Message (MSG_POSTERR, "GPStoNPS failed at bssp || bssp->_class");
    return;
  }

  GetSeqEntryParent (top, &parentptr, &parenttype);

  /* point to genomic Bioseq component of gps */

  sep = bssp->seq_set;
  if (sep == NULL || sep->choice != 1) {
    Message (MSG_POSTERR, "GPStoNPS failed at sep || sep->choice");
    return;
  }

  /* unlink nuc-prot sets, etc., from genomic Bioseq */

  dum = BioseqSetNew ();
  if (dum == NULL) {
    Message (MSG_POSTERR, "GPStoNPS failed at BioseqSetNew");
    return;
  }
  dum->_class = 1;
  dum->seq_set = sep->next;
  sep->next = NULL;

  last = sep;
  VisitBioseqsInSet (dum, (Pointer) &last, TakeProteinsFromGPS);

  bssp->_class = BioseqseqSet_class_nuc_prot;

  SeqMgrLinkSeqEntry (top, parenttype, parentptr);

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  VisitFeaturesInSet (bssp, NULL, ClearRnaProducts);

  move_cds (top);

  /* in case result has no proteins, demote to bioseq */

  RenormalizeNucProtSets (top, TRUE);

  /* cleanup original nuc-prot sets */

  BioseqSetFree (dum);
}

static void GeneralToNote (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr  bsp;
  Char       buf [41];
  DbtagPtr   dbt;
  size_t     len;
  SeqIdPtr   sip;
  CharPtr    str;

  if (sfp == NULL || sfp->product == NULL) return;
  if (sfp->data.choice != SEQFEAT_RNA) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice != SEQID_GENERAL) continue;
    dbt = (DbtagPtr) sip->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringICmp (dbt->db, "TMSMART") == 0 || StringICmp (dbt->db, "NCBIFILE") == 0) continue;

    SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);

    if (sfp->comment == NULL) {
      sfp->comment = StringSave (buf);
    } else {
      len = StringLen (sfp->comment) + StringLen (buf) + 5;
      str = MemNew (sizeof (Char) * len);
      StringCpy (str, sfp->comment);
      StringCat (str, "; ");
      StringCat (str, buf);
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = str;
    }
  }
}

static SeqEntryPtr PropagateDescsFromGenBankSet (
  SeqEntryPtr sep
)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqEntryPtr   firstsep = NULL;
  SeqEntryPtr   seqentry;
  ValNodePtr    sourcedescr;

  if (sep == NULL) return NULL;
  if (! IS_Bioseq_set (sep)) return sep;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return sep;
  sourcedescr = bssp->descr;
  if (sourcedescr == NULL) return sep;
  firstsep = bssp->seq_set;
  seqentry = firstsep;
  while (seqentry != NULL) {
    if (seqentry->data.ptrvalue != NULL) {
      if (seqentry->choice == 1) {
        bsp = (BioseqPtr) seqentry->data.ptrvalue;
        ValNodeLink (&(bsp->descr),
                     AsnIoMemCopy ((Pointer) sourcedescr,
                                   (AsnReadFunc) SeqDescrAsnRead,
                                   (AsnWriteFunc) SeqDescrAsnWrite));
      } else if (seqentry->choice == 2) {
        bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
        ValNodeLink (&(bssp->descr),
                     AsnIoMemCopy ((Pointer) sourcedescr,
                                   (AsnReadFunc) SeqDescrAsnRead,
                                   (AsnWriteFunc) SeqDescrAsnWrite));
      }
    }
    seqentry = seqentry->next;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  bssp->descr = SeqDescrFree (bssp->descr);
  NormalizeDescriptorOrder (sep);
  return firstsep;
}

typedef struct srcdata {
  Boolean  isSeqId;
  Boolean  isOrganism;
  Uint1    orgmodType;
  Uint1    subsourceType;
} SrcData, PNTR SrcDataPtr;

static void ParseOneOrgLabel (
  SrcDataPtr field,
  CharPtr label
)

{
  Int2  i;

  if (field == NULL || StringHasNoText (label)) return;

  if (StringICmp (label, "local_id") == 0 ||
      StringICmp (label, "local id") == 0 ||
      StringICmp (label, "SequenceID") == 0 ||
      StringICmp (label, "Sequence_ID") == 0 ||
      StringICmp (label, "Sequence ID") == 0 ||
      StringICmp (label, "SeqID") == 0 ||
      StringICmp (label, "Seq_ID") == 0 ||
      StringICmp (label, "Seq ID") == 0) {
    field->isSeqId = TRUE;
    return;
  }
  if (StringICmp (label, "organism") == 0) {
    field->isOrganism = TRUE;
    return;
  }

  i = EquivalentOrgMod (label);
  if (i != 0) {
    field->orgmodType = (Uint1) i;
    return;
  }
  i = EquivalentSubSource (label);
  if (i != 0) {
    field->subsourceType = (Uint1) i;
    return;
  }
  if (StringICmp (label, "note") == 0) {
    field->subsourceType = (Uint1) SUBSRC_other;
  }
}

static void ProcessSourceTable (
  FILE *fp
)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  CharPtr       columns [80];
  FileCache     fc;
  SrcData       fields [80];
  Int2          i, numfields;
  Char          line [4095];
  OrgModPtr     omp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  CharPtr       ptr, str;
  SeqDescrPtr   sdp;
  SeqIdPtr      sip;
  SubSourcePtr  ssp;

  if (fp == NULL) return;

  MemSet ((Pointer) fields, 0, sizeof (fields));
  numfields = 0;

  FileCacheSetup (&fc, fp);

  /* read first line with field names */

  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  if (str == NULL) return;

  TrimSpacesAroundString (str);
  while (StringDoesHaveText (str) && numfields < 78) {
    ptr = StringChr (str, '\t');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    ParseOneOrgLabel (&(fields [numfields]), str);
    numfields++;
    str = ptr;
  }

  if (! fields [0].isSeqId) return;

  /* read remaining lines with source data */

  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL) {

    MemSet ((Pointer) columns, 0, sizeof (columns));

    TrimSpacesAroundString (str);
    i = 0;
    while (StringDoesHaveText (str) && i < numfields) {
      ptr = StringChr (str, '\t');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
      TrimSpacesAroundString (str);
      columns [i] = str;
      i++;
      str = ptr;
    }

    if (StringDoesHaveText (columns [0])) {
      sip = MakeSeqID (columns [0]);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          biop = NULL;
          sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
          if (sdp != NULL) {
            biop = (BioSourcePtr) sdp->data.ptrvalue;
          }
          if (biop == NULL) {
            biop = BioSourceNew ();
            if (biop != NULL) {
              SeqDescrAddPointer (&(bsp->descr), Seq_descr_source, (Pointer) biop);
            }
          }
          if (biop != NULL) {
            for (i = 1; i < numfields; i++) {
              if (StringHasNoText (columns [i])) continue;
              if (fields [i].isOrganism) {
                if (biop->org == NULL) {
                  biop->org = OrgRefNew ();
                }
                orp = biop->org;
                if (orp != NULL) {
                  orp->taxname = MemFree (orp->taxname);
                  orp->taxname = StringSave (columns [i]);
                }
              } else if (fields [i].orgmodType > 0) {
                if (biop->org == NULL) {
                  biop->org = OrgRefNew ();
                }
                orp = biop->org;
                if (orp != NULL) {
                  if (orp->orgname == NULL) {
                    orp->orgname = OrgNameNew ();
                  }
                  onp = orp->orgname;
                  if (onp != NULL) {
                    omp = OrgModNew ();
                    if (omp != NULL) {
                      omp->subtype = (Uint1) fields [i].orgmodType;
                      omp->subname = StringSave (columns [i]);
                      omp->next = onp->mod;
                      onp->mod = omp;
                    }
                  }
                }
              } else if (fields [i].subsourceType > 0) {
                ssp = SubSourceNew ();
                if (ssp != NULL) {
                  ssp->subtype = (Uint1) fields [i].subsourceType;
                  ssp->name = StringSave (columns [i]);
                  ssp->next = biop->subtype;
                  biop->subtype = ssp;
                }
              }
            }
          }
        }
        sip = SeqIdFree (sip);
      }
    }

    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }
}

static SeqDescrPtr GetDescriptorTypeAlreadyInList (
  Uint1 descr_choice,
  SeqDescrPtr list
)

{
  while (list != NULL && list->choice != descr_choice) {
    list = list->next;
  }
  return list;
}

static void AddTemplateDescriptors (
  SeqDescrPtr PNTR current_list,
  SeqDescrPtr new_list,
  Boolean copy
)

{
  SeqDescrPtr  dsc, sdp_next, sdp;

  if (current_list == NULL || new_list == NULL) return;

  for (sdp = new_list; sdp != NULL; sdp = sdp_next) {
    sdp_next = sdp->next;
    if (sdp->choice == Seq_descr_molinfo) continue;
    if (sdp->choice == Seq_descr_source &&
        GetDescriptorTypeAlreadyInList (Seq_descr_source, *current_list) != NULL) continue;
    sdp->next = NULL;
    if (copy) {
      dsc = AsnIoMemCopy ((Pointer) sdp,
                          (AsnReadFunc) SeqDescrAsnRead,
                          (AsnWriteFunc) SeqDescrAsnWrite);
    } else {
      dsc = sdp;
    }
    ValNodeLink (current_list, (Pointer) dsc);
    sdp->next = sdp_next;
  }
}

static void GenomizeSeqId (
  SeqIdPtr sip,
  Pointer userdata
)

{
  CharPtr      accn = NULL;
  CharPtr      center;
  DbtagPtr     dbt;
  ObjectIdPtr  oip;

  if (sip == NULL || sip->choice != SEQID_LOCAL) return;
  center = (CharPtr) userdata;
  if (StringHasNoText (center)) return;

  oip = (ObjectIdPtr) sip->data.ptrvalue;
  if (oip == NULL) return;
  accn = oip->str;
  if (StringHasNoText (accn)) return;

  dbt = DbtagNew ();
  if (dbt == NULL) return;
  oip = ObjectIdNew ();
  if (oip == NULL) return;
  oip->str = StringSave (accn);
  dbt->db = StringSave (center);
  dbt->tag = oip;

  sip->data.ptrvalue = ObjectIdFree ((ObjectIdPtr) sip->data.ptrvalue);
  sip->data.ptrvalue = (Pointer) dbt;
  sip->choice = SEQID_GENERAL;
}

static void GenomizeFeatureSeqIds (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqLoc (sfp->location, userdata, GenomizeSeqId);
}

static void GenomizeGraphSeqIds (
  SeqGraphPtr sgp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqGraph (sgp, userdata, GenomizeSeqId);
}

static void MakeGenomeCenterID (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr  center;

  if (bsp == NULL) return;
  center = (CharPtr) userdata;
  if (StringHasNoText (center)) return;

  VisitSeqIdsInBioseq (bsp, userdata, GenomizeSeqId);
  SeqMgrReplaceInBioseqIndex (bsp);
  VisitFeaturesOnBsp (bsp, userdata, GenomizeFeatureSeqIds);
  VisitGraphsOnBsp (bsp, userdata, GenomizeGraphSeqIds);
}

static void MakeAccessionID (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr     accn;
  ValNodePtr  generalIDs;
  SeqIdPtr    sip;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  accn = (CharPtr) userdata;
  if (StringHasNoText (accn)) return;

  /* if existing accession, coerce all SeqIds */

  sip = SeqIdFromAccession (accn, INT2_MIN, NULL);
  if (sip == NULL) return;
  generalIDs = ValNodeExtractList (&(bsp->id), SEQID_GENERAL);
  bsp->id = SeqIdSetFree (bsp->id);
  bsp->id = sip;
  if (generalIDs != NULL) {
    ValNodeLink (&(bsp->id), generalIDs);
  }
  SeqMgrReplaceInBioseqIndex (bsp);
  VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
  VisitGraphsOnBsp (bsp, (Pointer) bsp->id, CorrectGraphSeqIds);
}

static void FindCreateDate (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BoolPtr  has_create_dateP;

  if (sdp == NULL || sdp->choice != Seq_descr_create_date || userdata == NULL) return;
  has_create_dateP = (BoolPtr) userdata;
  *has_create_dateP = TRUE;
}

static void ConvertStructuredComment (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  SeqDescrPtr    com;
  CharPtr        prefix = NULL;
  CharPtr        str;
  UserObjectPtr  uop = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;

  if (StringStr (str, "##HIVData-START##") != NULL &&
      StringStr (str, "##HIVData-END##") != NULL) {
    prefix = StringStr (str, "##HIVData-START##");
    uop = ParseStringIntoStructuredComment (NULL, str, "##HIVData-START##",
                                            "##HIVData-END##");
  } else if (StringStr (str, "##FluData-START##") != NULL &&
             StringStr (str, "##FluData-END##") != NULL) {
    prefix = StringStr (str, "##FluData-START##");
    uop = ParseStringIntoStructuredComment (NULL, str, "##FluData-START##",
                                            "##FluData-END##");
  }
  if (uop == NULL) return;

  /* if there is text before prefix, truncate existing comment and append user object */

  if (prefix != NULL) {
    *prefix = '\0';
    TrimSpacesAroundString (str);
    if (StringDoesHaveText (str)) {
      com = SeqDescrNew (NULL);
      if (com != NULL) {
        com->choice = Seq_descr_user;
        com->data.ptrvalue = uop;
        com->next = sdp->next;
        sdp->next = com;
        return;
      }
    }
  }

  /* if entire comment was structured, replace existing descriptor with user object */

  MemFree (sdp->data.ptrvalue);
  sdp->choice = Seq_descr_user;
  sdp->data.ptrvalue = uop;
}

static void CleanUpLatLonAndCountry (
  BioSourcePtr biop,
  Pointer userdata
)

{
  CharPtr       fix_lat_lon;
  Boolean       format_ok = FALSE;
  CharPtr       lat_lon = NULL;
  Boolean       lat_in_range = FALSE;
  Boolean       lon_in_range = FALSE;
  CharPtr PNTR  list;
  CharPtr       new_country;
  SubSourcePtr  ssp;

  if (biop == NULL) return;
  list = (CharPtr PNTR) userdata;
  if (list == NULL) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_country && StringDoesHaveText (ssp->name)) {
      new_country = GetCountryFix (ssp->name, list);
      if (new_country != NULL) {
        ssp->name = MemFree (ssp->name);
        ssp->name = new_country;
      }
    } else if (ssp->subtype == SUBSRC_lat_lon && StringDoesHaveText (ssp->name)) {
      lat_lon = ssp->name;
      IsCorrectLatLonFormat (lat_lon, &format_ok, &lat_in_range, &lon_in_range);
      if (! format_ok) {
        fix_lat_lon = FixLatLonFormat (lat_lon);
        if (fix_lat_lon != NULL) {
          ssp->name = MemFree (ssp->name);
          ssp->name = fix_lat_lon;
        }
      }
    }
  }
}

static void LookupPubdesc (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid = 0;
  ValNodePtr       vnp;

  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        /* ignore obsolete muids */
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        /* return on real pub */
        return;
        break;
    }
  }

  if (pmid == 0) return;

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
  }

  PubmedEntryFree (pep);
}



#ifdef INTERNAL_NCBI_ASNDISC
const PerformDiscrepancyTest taxlookup = CheckTaxNamesAgainstTaxDatabase;
#else
const PerformDiscrepancyTest taxlookup = NULL;
#endif


static void CleanupCollectionDatesMonthFirst (BioSourcePtr biop, Pointer data)
{
  SubSourcePtr ssp;
  CharPtr      reformatted_date = NULL;

  if (biop == NULL) return;

  ssp = biop->subtype;
  while (ssp != NULL)
  {
    if (ssp->subtype == SUBSRC_collection_date)
    {
      reformatted_date = ReformatDateStringEx (ssp->name, TRUE, NULL);
      if (reformatted_date != NULL)
      {
        ssp->name = MemFree (ssp->name);
        ssp->name = reformatted_date;
      }
    }
    ssp = ssp->next;
  }
}


static void CleanupCollectionDatesDayFirst (BioSourcePtr biop, Pointer data)
{
  SubSourcePtr ssp;
  CharPtr      reformatted_date = NULL;

  if (biop == NULL) return;

  ssp = biop->subtype;
  while (ssp != NULL)
  {
    if (ssp->subtype == SUBSRC_collection_date)
    {
      reformatted_date = ReformatDateStringEx (ssp->name, FALSE, NULL);
      if (reformatted_date != NULL)
      {
        ssp->name = MemFree (ssp->name);
        ssp->name = reformatted_date;
      }
    }
    ssp = ssp->next;
  }
}


static void ValNodeLinkCopy (ValNodePtr PNTR list1, ValNodePtr list2)
{
  if (list1 == NULL) return;
  while (list2 != NULL)
  {
    ValNodeAddPointer (list1, list2->choice, list2->data.ptrvalue);
    list2 = list2->next;
  }
}

static ValNodePtr FindItemListForClickableItemCategory (ValNodePtr list, CharPtr category_fmt)
{
  ClickableItemPtr cip;
  ValNodePtr       vnp;
  ValNodePtr       item_list = NULL;
  CharPtr          cp;

  if (StringLen (category_fmt) < 2) {
    return NULL;
  }
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      if (cip->description != NULL) {
        /* skip number at beginning of category title */
        cp = cip->description;
        while (isdigit (*cp)) {
          cp++;
        }
        if (StringCmp (cp, category_fmt + 2) == 0) {
          ValNodeLinkCopy (&item_list, cip->item_list);
        }
      }
      ValNodeLink (&item_list, FindItemListForClickableItemCategory (cip->subcategories, category_fmt));
    }
  }
  return item_list;
}


static void DoTbl2AsnCleanup (SeqEntryPtr sep, CleanupArgsPtr c)
{
  ValNodePtr sep_list = NULL;
  ValNodePtr discrepancy_list = NULL, item_list = NULL, vnp;
  SeqFeatPtr sfp;

  if (sep == NULL || c == NULL) {
    return;
  }

  if (c->collection_dates) {
    if (c->collection_dates_month_first) {
      VisitBioSourcesInSep (sep, NULL, CleanupCollectionDatesMonthFirst);
    } else {
      VisitBioSourcesInSep (sep, NULL, CleanupCollectionDatesDayFirst);
    }
  }
  if (c->add_notes_to_overlapping_cds_without_abc) {
    ValNodeAddPointer (&sep_list, 0, sep);
    SeqMgrIndexFeatures (ObjMgrGetEntityIDForChoice (sep), NULL);
    AddOverlappingCodingRegionDiscrepancies (&discrepancy_list, sep_list);
    sep_list = ValNodeFree (sep_list);
    item_list = FindItemListForClickableItemCategory (discrepancy_list, kOverlappingCDSNeedsNoteFmt);
    discrepancy_list = FreeClickableList (discrepancy_list);
    for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          SetStringValue (&(sfp->comment), kOverlappingCDSNoteText, ExistingTextOption_append_semi);
        }
      }
    }
    item_list = ValNodeFree (item_list);
  }
}


static void SeqEntryHasConflictingIDsCallback (BioseqPtr bsp, Pointer data)
{
  CharPtr msg, fmt = "SeqID %s is present on multiple Bioseqs in record";
  BioseqPtr bsp2;
  SeqIdPtr sip;
  DbtagPtr dbt;
  Char     buf[100];

  if (bsp == NULL || data == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL 
        && (dbt = (DbtagPtr) sip->data.ptrvalue) != NULL
        && StringICmp (dbt->db, "NCBIFILE") == 0) {
      continue;
        }
    bsp2 = BioseqFindSpecial (sip);
    if (bsp2 != NULL && bsp2 != bsp) {
      SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (buf)));
      sprintf (msg, fmt, buf);
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, msg);
    }
  }
}


static Boolean SeqEntryHasConflictingIDs (SeqEntryPtr sep)
{
  ValNodePtr errs = NULL, vnp;

  VisitBioseqsInSep (sep, &errs, SeqEntryHasConflictingIDsCallback);
  if (errs == NULL) {
    return FALSE;
  } else {
    ValNodeUnique (&errs, SortVnpByString, ValNodeFreeData);
    for (vnp = errs; vnp != NULL; vnp = vnp->next) {
      Message (MSG_POSTERR, vnp->data.ptrvalue);
    }
    errs = ValNodeFreeData (errs);
    return TRUE;
  }
}


static void ProcessOneRecord (
  SubmitBlockPtr sbp,
  PubdescPtr pdp,
  BioSourcePtr src,
  CharPtr directory,
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqDescrPtr sdphead,
  TblArgsPtr tbl,
  TextFsaPtr gotags,
  AsnIoPtr aip,
  CharPtr outfile
)

{
  AsnTypePtr         atp_bssse;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp = NULL;
  Char               buf [256];
  SeqMgrFeatContext  context;
  Pointer            dataptr;
  Uint2              datatype, entityID;
  SeqDescrPtr        descr;
  DatePtr            dp;
  BioseqSetPtr       dssp;
  Boolean            failure = FALSE;
  FileCache          fc;
  FILE               *fp;
  Int2               genCode;
  Boolean            goOn;
  SeqEntryPtr        gsep = NULL;
  Boolean            has_create_date;
  SeqGraphPtr        lastsgp;
  Int4               linenum = 0;
  CharPtr PNTR       list;
  CharPtr            localname = NULL;
  MolInfoPtr         mip;
  ErrSev             msev;
  Boolean            nonewline;
  BioseqPtr          nucbsp;
  ObjMgrDataPtr      omdp;
  CharPtr            organism;
  OrgRefPtr          orp;
  BioseqPtr          protbsp;
  SeqEntryPtr        protsep;
  CharPtr            ptr;
  SeqAnnotPtr        sap;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  CharPtr            sfx = NULL;
  SeqGraphPtr        sgp;
  SeqIdPtr           sip;
  SeqSubmitPtr       sub;
  SimpleSeqPtr       ssp;
  CharPtr            str;
  CharPtr            tblfile = NULL;
  SeqEntryPtr        tmp;
  MolInfoPtr         template_molinfo = NULL;
  ValNodePtr         cmt_errors, vnp;

  fp = OpenOneFile (directory, base, suffix);
  if (fp == NULL) return;

  if (tbl->logtoterminal) {
    Message (MSG_POSTERR, "File %s", base);
  }

  /* if genomic product set, make parent set */

  if (tbl->genprodset) {
    bssp = BioseqSetNew ();
    if (bssp == NULL) return;
    bssp->_class = BioseqseqSet_class_gen_prod_set;

    gsep = SeqEntryNew ();
    if (gsep == NULL) return;
    gsep->choice = 2;
    gsep->data.ptrvalue = (Pointer) bssp;
    SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, gsep);
  }

  if (tbl->seqidfromfile) {
    localname = base;
  }

  /* find MolInfo from template, if there is any */
  sdp = sdphead;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
    sdp = sdp->next;
  }
  if (sdp != NULL) {
    template_molinfo = (MolInfoPtr) sdp->data.ptrvalue;
  }

  /* read one or more ASN.1 or FASTA sequence files */

  if (tbl->fastaset) {
    entityID = ProcessBulkSet (fp, src, tbl, template_molinfo);
  } else if (tbl->deltaset) {
    entityID = ProcessDeltaSet (fp, src, tbl, localname, gsep, template_molinfo);
  } else if (tbl->alignset) {
    entityID = ProcessAlignSet (fp, src, tbl, template_molinfo);
  } else if (tbl->gapped) {
    entityID = ProcessGappedSet (fp, src, tbl, gsep, template_molinfo);
  } else if (tbl->phrapace) {
    entityID = ProcessPhrapAce (fp, src, tbl, localname, gsep, template_molinfo, directory, base);
  } else if (tbl->raw2delt) {
    entityID = ProcessRaw2Delt (fp, src, tbl, localname, gsep, template_molinfo);
  } else {
    entityID = ProcessOneAsn (fp, src, tbl, localname, gsep, template_molinfo);
  }
  FileClose (fp);

  if (entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (SeqEntryHasConflictingIDs (sep)) {
    return;
  }

  if (tbl->dotaxlookup) {
    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) {

      /* optionally do network taxonomy lookup - prior to instantiating mRNA and protein titles */

      Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
    }
  }

  if (tbl->dopublookup) {
    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) {

      /* optionally do network publication lookup of just PMID references */

      VisitPubdescsInSep (sep, NULL, LookupPubdesc);
    }
  }

  organism = NULL;
  if (tbl->genprodset) {
    descr = ExtractBioSourceAndPubs (bssp->seq_set);
    for (sdp = descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice != Seq_descr_source) continue;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop == NULL) continue;
      orp = biop->org;
      if (orp == NULL) continue;
      if (StringDoesHaveText (orp->taxname)) {
        organism = orp->taxname;
      }
    }
    ReplaceBioSourceAndPubs (gsep, descr);
  }

  /* read one or more feature tables from .tbl file */

  if (StringDoesHaveText (tbl->tableFile)) {
    fp = FileOpen (tbl->tableFile, "r");
    tblfile = tbl->tableFile;
  } else {
    fp = OpenOneFile (directory, base, ".tbl");
    tblfile = base;
    sfx = ".tbl";
  }
  if (fp != NULL) {

    /* indexing needed to find segmented bsp if location is on part */

    sep = GetTopSeqEntryForEntityID (entityID);

    SeqMgrIndexFeatures (entityID, NULL);

    while ((! failure) && (dataptr = ReadFeatureTableFile (fp, &datatype, NULL, &linenum, &failure)) != NULL) {
      if (datatype == OBJ_SEQANNOT) {

        sap = (SeqAnnotPtr) dataptr;
        ProcessOneAnnot (sap, entityID, tbl);

      } else {
        ObjMgrFree (datatype, dataptr);
      }
    }
    FileClose (fp);
    sep = GetTopSeqEntryForEntityID (entityID);


    if (failure) {
      if (StringHasNoText (tblfile)) {
        tblfile = "?";
      }
      ptr = StringRChr (tblfile, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        tblfile = ptr;
      }
      Message (MSG_POSTERR, "Bad feature table at line %ld of file %s%s", (long) linenum, tblfile, sfx);
    }
  }

  /* if genomic product set, copy CDS into nucprot sets */

  if (tbl->genprodset) {
    /* need to reindex to get mRNA and CDS features from cDNA and protein */
    SeqMgrIndexFeatures (entityID, NULL);
    VisitSetsInSet (bssp, (Pointer) tbl, MakeNucProtCDS);
  }

  /* read source qualifiers for set of sequences from .src file */

  fp = OpenOneFile (directory, base, ".src");
  if (fp != NULL) {

    ProcessSourceTable (fp);

    FileClose (fp);
  }

  /* read structured comments from .cmt file */
  fp = OpenOneFile (directory, base, ".cmt");
  if (fp != NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
    cmt_errors = CreateStructuredCommentsFromFile (fp, sep, tbl->apply_cmt_to_all);
    FileClose (fp);
    if (cmt_errors != NULL) {
      for (vnp = cmt_errors; vnp != NULL; vnp = vnp->next) {
        Message (MSG_POSTERR, "Error processing structured comment (.cmt) file: %s", vnp->data.ptrvalue);
      }
      cmt_errors = ValNodeFreeData (cmt_errors);
    }
  }

  /* read one or more protein sequences from .pep file */

  fp = OpenOneFile (directory, base, ".pep");
  if (fp != NULL) {

    /* indexing needed to find CDS from protein product to set conflict flag */

    SeqMgrIndexFeatures (entityID, NULL);

    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, TRUE, TRUE, TRUE)) != NULL) {
      if (datatype == OBJ_FASTA) {

        ssp = (SimpleSeqPtr) dataptr;
        ReplaceOnePeptide (ssp, tbl->conflict, tbl->genprodset);
        SimpleSeqFree (ssp);

      } else {
        ObjMgrFree (datatype, dataptr);
      }
    }
    FileClose (fp);
  }

  /* read one or more RNA sequences from .rna file */

  fp = OpenOneFile (directory, base, ".rna");
  if (fp != NULL) {

    /* indexing needed to find mRNA from transcript product to set RNA editing exception */

    SeqMgrIndexFeatures (entityID, NULL);

    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, TRUE, TRUE, TRUE)) != NULL) {
      if (datatype == OBJ_FASTA) {

        ssp = (SimpleSeqPtr) dataptr;
        ReplaceOneRNA (ssp, tbl->conflict);
        SimpleSeqFree (ssp);

      } else {
        ObjMgrFree (datatype, dataptr);
      }
    }
    FileClose (fp);
  }

  /* read one or more protein sequences from .prt file */

  fp = OpenOneFile (directory, base, ".prt");
  if (fp != NULL) {

    SeqMgrIndexFeatures (entityID, NULL);

    sep = GetTopSeqEntryForEntityID (entityID);
    nucbsp = FindNucBioseq (sep);
    if (nucbsp != NULL) {
      BioseqToGeneticCode (nucbsp, &genCode, NULL, NULL, NULL, 0, NULL);
      SetBatchSuggestNucleotide (nucbsp, genCode);

      descr = ExtractBioSourceAndPubs (sep);

      while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, TRUE, TRUE, FALSE)) != NULL) {
        if (datatype == OBJ_BIOSEQ) {

          protbsp = (BioseqPtr) dataptr;
          protsep = SeqMgrGetSeqEntryForData (protbsp);
          mip = MolInfoNew ();
          if (mip != NULL) {
            mip->biomol = 8;
            mip->tech = 13;
            sdp = CreateNewDescriptor (protsep, Seq_descr_molinfo);
            if (sdp != NULL) {
              sdp->data.ptrvalue = (Pointer) mip;
            }
          }
          AddSeqEntryToSeqEntry (sep, protsep, TRUE);
          SuggestOnePeptide (nucbsp, protbsp, genCode);

        } else {
          ObjMgrFree (datatype, dataptr);
        }
      }

      ClearBatchSuggestNucleotide ();

      ReplaceBioSourceAndPubs (sep, descr);
    }
    FileClose (fp);

    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* read one or more quality score blocks from .qvl file */

  fp = OpenOneFile (directory, base, ".qvl");
  if (fp != NULL) {

    FileCacheSetup (&fc, fp);

    goOn = TRUE;
    while (goOn) {
      str = FileCacheReadLine (&fc, buf, sizeof (buf), &nonewline);
      if (str == NULL) {
        goOn = FALSE;
      } else if (StringDoesHaveText (str)) {
        if (str [0] == '>') {
          ptr = StringChr (str, ' ');
          if (ptr == NULL) {
            ptr = StringChr (str, '\t');
          }
          if (ptr != NULL) {
            *ptr = '\0';
          }
          sip = MakeSeqID (str + 1);
          bsp = BioseqFind (sip);
          if (bsp != NULL) {
            sgp = ReadPhrapQualityFC (&fc, bsp);
            if (sgp != NULL) {
              for (sap = bsp->annot; sap != NULL; sap = sap->next) {
                if (sap->type == 3) {
                  for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
                    continue;
                  }
                  lastsgp->next = sgp;
                  break;
                }
              }
              if (sap == NULL) {
                if (bsp->annot != NULL) {
                  for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
                    continue;
                  }
                  sap->next = NewGraphSeqAnnot ("Phrap Graph", sgp);
                } else {
                  bsp->annot = NewGraphSeqAnnot ("Phrap Graph", sgp);
                }
              }
            }
          }
          SeqIdFree (sip);
        }
      }
    }
    FileClose (fp);
  }

  /* finish processing */

  if (sbp == NULL) {
    omdp = ObjMgrGetData (entityID);
    if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {

      /* if read a Seq-submit, write out a Seq-submit */

      sub = (SeqSubmitPtr) omdp->dataptr;
      if (sub != NULL && sub->datatype == 1) {
        sbp = sub->sub;
      }
    }
  }

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL) {

    if (tbl->gnltonote) {
      VisitFeaturesInSep (sep, NULL, GeneralToNote);
    }

    if (tbl->gpstonps) {
      GPStoNPS (sep, entityID);
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (! tbl->genprodset) {
      VisitFeaturesInSep (sep, NULL, RemoveGBQualIDs);
    }
    if (sdphead != NULL) {
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        AddTemplateDescriptors (&(bsp->descr), sdphead, TRUE);
      } else if (IS_Bioseq_set (sep)) {
        dssp = (BioseqSetPtr) sep->data.ptrvalue;
        AddTemplateDescriptors (&(dssp->descr), sdphead, TRUE);
      }
    }
    dp = DateCurr ();
    if (dp != NULL) {
      has_create_date = FALSE;
      VisitDescriptorsInSep (sep, (Pointer) &has_create_date, FindCreateDate);
      if (has_create_date) {
        sdp = CreateNewDescriptor (sep, Seq_descr_update_date);
      } else {
        sdp = CreateNewDescriptor (sep, Seq_descr_create_date);
      }
      if (sdp != NULL) {
        sdp->data.ptrvalue = (Pointer) dp;
      }
    }

    /* read one or more descriptors from .dsc file */

    fp = OpenOneFile (directory, base, ".dsc");
    if (fp != NULL) {

      while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, TRUE, TRUE, TRUE)) != NULL) {
        if (datatype == OBJ_SEQDESC) {

          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            AddTemplateDescriptors (&(bsp->descr), (SeqDescrPtr) dataptr, FALSE);
          } else if (IS_Bioseq_set (sep)) {
            dssp = (BioseqSetPtr) sep->data.ptrvalue;
            AddTemplateDescriptors (&(dssp->descr), (SeqDescrPtr) dataptr, FALSE);
          }

        } else {
          ObjMgrFree (datatype, dataptr);
        }
      }
      FileClose (fp);
    }

    msev = ErrSetMessageLevel (SEV_MAX);
    move_cds (sep);

    /* if reading nucleotide and protein tables, remove duplicate prot feat */
    VisitBioseqsInSep (sep, NULL, RemoveDupProtFeats);
    DeleteMarkedObjects (entityID, 0, NULL);

    /* need to reindex before extending CDS to stop codon */
    SeqMgrIndexFeatures (entityID, NULL);
    CdCheck (sep, NULL);

    /* need to reindex before copying genes, instantiating protein titles */
    SeqMgrIndexFeatures (entityID, NULL);
    EntryChangeImpFeat (sep);

    /* find locus for any gene xrefs that only have locus_tag */
    VisitFeaturesInSep (sep, NULL, FillInPartialGeneXref);

    if (tbl->removeunnecxref) {
      /* if not removed, xref will prevent locus, maploc, dbxref from being copied */
      VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
    }

    if (tbl->genprodset) {
      VisitFeaturesInSep (sep, NULL, CopyGene);
    }
    if (tbl->genprodset) {
      /* currently copying ncRNA feature onto product */
      VisitFeaturesInSep (sep, NULL, CopyNcRna);
    }
    if (! tbl->genprodset) {
    VisitFeaturesInSep (sep, NULL, ClearRnaProducts);
    }

    if (tbl->removeunnecxref) {
      /* need to reindex before removing unnecesary gene xrefs in nuc-prot sets */
      SeqMgrIndexFeatures (entityID, NULL);

      VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
    }

    if (! tbl->relaxed) {
      list = GetValidCountryList ();
      VisitBioSourcesInSep (sep, (Pointer) list, CleanUpLatLonAndCountry);
    }

    /* need to reindex so hypothetical protein titles pick up locus_tag */
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);

    if (tbl->genprodset) {
      /* need to reindex before instantiating mRNA titles */
      SeqMgrIndexFeatures (entityID, NULL);
      bsp = FindNucBioseq (sep);

      if (tbl->smarttitle) {
        MakeSmartRnaTitles (bsp, organism);
      } else {
        sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &context);
        while (sfp != NULL) {
          AddRnaTitles (sfp, organism);
          sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, 0, &context);
        }
      }
    }

    if (StringDoesHaveText (tbl->center)) {
      VisitBioseqsInSep (sep, tbl->center, MakeGenomeCenterID);
    }

    if (StringDoesHaveText (tbl->accn)) {
      bsp = FindNucBioseq (sep);
      MakeAccessionID (bsp, tbl->accn);
    }

    VisitDescriptorsInSep (sep, NULL, ConvertStructuredComment);

    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
    ErrSetMessageLevel (msev);
    /*
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    */
    ConvertFullLenSourceFeatToDesc (sep);
    ConvertFullLenPubFeatToDesc (sep);
    if (tbl->linkbyoverlap) {
      SeqMgrIndexFeatures (entityID, NULL);
      LinkCDSmRNAbyOverlap (sep);
    } else if (tbl->linkbyproduct) {
      SeqMgrIndexFeatures (entityID, NULL);
      LinkCDSmRNAbyProduct (sep);
    }

    DoTbl2AsnCleanup (sep, &(tbl->cleanup_args));
    NormalizeDescriptorOrder (sep);

    if (StringHasNoText (results)) {
      results = directory;
    }

    if (aip != NULL) {
      atp_bssse = AsnFind ("Bioseq-set.seq-set.E");
      if (atp_bssse == NULL) {
        Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
      } else if (tbl->fastaset && tbl->whichclass == 0) {
        /* already has genbank wrapper, write individual components */
        tmp = PropagateDescsFromGenBankSet (sep);
        SeqMgrClearFeatureIndexes (entityID, NULL);
        while (tmp != NULL) {
          SeqEntryAsnWrite (tmp, aip, atp_bssse);
          tmp = tmp->next;
        }
      } else {
        SeqEntryAsnWrite (sep, aip, atp_bssse);
      }
    } else {
      if (tbl->fastaset && tbl->whichclass == 0) {
        PropagateDescsFromGenBankSet (sep);
        SeqMgrClearFeatureIndexes (entityID, NULL);
      }
      WriteOneFile (results, base, ".sqn", outfile, sep, sbp, tbl->save_bioseq_set);
    }

    if (HasGoTermsInNote (sep, gotags)) {
      Message (MSG_OK, "Illegal GO term format detected in note - contact database for instructions");
    }

    if (tbl->global_report != NULL) {
      AddSeqEntryToGlobalDiscrepReport (sep, tbl->global_report, base);
    }

    if (tbl->validate || tbl->flatfile || tbl->genereport || tbl->validate_barcode) {
      if (pdp != NULL) {

        /* copy in citsub as publication for validator and flatfile */

        sdp = CreateNewDescriptor (sep, Seq_descr_pub);
        if (sdp != NULL) {
          sdp->data.ptrvalue = AsnIoMemCopy ((Pointer) pdp,
                                             (AsnReadFunc) PubdescAsnRead,
                                             (AsnWriteFunc) PubdescAsnWrite);
        }
      }
      SeqMgrIndexFeatures (entityID, 0);
      if (tbl->flatfile) {
        Message (MSG_POST, "Flatfile %s\n", base);
        FlatfileOneFile (results, base, ".gbf", sep);
      }
      if (tbl->validate || tbl->validate_barcode) {
        Message (MSG_POST, "Validating %s\n", base);
        ValidateOneFile (results, base, ".val", sep, tbl->validate, tbl->relaxed, tbl->validate_barcode);
      }
      if (tbl->genereport) {
        GeneReportOneFile (results, base, ".t2g", sep);
      }
    }
  }

  ObjMgrFreeByEntityID (entityID);
}



static CharPtr overwriteMsg = "Your template with a .sqn suffix will be overwritten.  Do you wish to continue?";

static Boolean TemplateOverwriteRisk (
  CharPtr filename,
  CharPtr single,
  CharPtr directory,
  CharPtr suffix
)

{
  Char     file [FILENAME_MAX], path [PATH_MAX];
  CharPtr  ptr;


  if (StringStr (filename, ".sqn") == NULL) return FALSE;
  if (StringDoesHaveText (single)) {
    StringNCpy_0 (file, filename, sizeof (file));
    ptr = StringStr (file, ".");
    if (ptr != NULL) {
      *ptr = '\0';
    }
    ptr = StringStr (single, ".");
    if (ptr != NULL) {
      StringCat (file, ptr);
    }
    if (StringCmp (file, single) == 0) return TRUE;
  } else if (StringDoesHaveText (directory)) {
    StringNCpy_0 (path, directory, sizeof (path));
    StringNCpy_0 (file, filename, sizeof (file));
    ptr = StringStr (file, ".");
    if (ptr != NULL) {
      *ptr = '\0';
    }
    StringCat (file, suffix);
    FileBuildPath (path, NULL, file);
    if (FileLength (path) > 0) return TRUE;
  }
  return FALSE;
}

static void FileRecurse (
  SubmitBlockPtr sbp,
  PubdescPtr pdp,
  BioSourcePtr src,
  CharPtr directory,
  CharPtr results,
  CharPtr suffix,
  Boolean recurse,
  SeqDescrPtr sdphead,
  TblArgsPtr tbl,
  TextFsaPtr gotags,
  AsnIoPtr aip,
  CharPtr outfile
)

{
  Char        path [PATH_MAX];
  CharPtr     ptr;
  CharPtr     str;
  ValNodePtr  head, vnp;

  /* get list of all files in source directory */

  head = DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringDoesHaveText (str)) {

        /* does filename have desired substring? */

        ptr = StringStr (str, suffix);

        if (ptr != NULL) {

          /* make sure detected suffix is really at end of filename */

          if (StringCmp (ptr, suffix) == 0) {
            *ptr = '\0';

            /* process file that has desired suffix (usually .fsa) */

            ProcessOneRecord (sbp, pdp, src, directory, results, str, suffix, sdphead, tbl, gotags, aip, outfile);
          }
        }
      }
    } else if (vnp->choice == 1 && recurse) {

      /* recurse into subdirectory */

      StringNCpy_0 (path, directory, sizeof (path));
      str = (CharPtr) vnp->data.ptrvalue;
      FileBuildPath (path, str, NULL);
      FileRecurse (sbp, pdp, src, path, results, suffix, recurse, sdphead, tbl, gotags, aip, outfile);
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);
}

static AsnTypePtr DoFirstPrefix (
  AsnIoPtr aip,
  SubmitBlockPtr sbp
)

{
  AsnTypePtr  atp_se, atp_ses, atp_ss, atp_ssd, atp_ssde, atp_ssdee, atp_sss, sep_atp, ssp_atp;
  DataVal     av;
  SeqEntry    se;
  SeqSubmit   ss;

  if (aip == NULL || sbp == NULL) return NULL;

  atp_ss = AsnFind ("Seq-submit");
  if (atp_ss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return NULL;
  }

  atp_sss = AsnFind ("Seq-submit.sub");
  if (atp_sss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.sub");
    return NULL;
  }

  atp_ssd = AsnFind ("Seq-submit.data");
  if (atp_ssd == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data");
    return NULL;
  }

  atp_ssde = AsnFind ("Seq-submit.data.entrys");
  if (atp_ssde == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys");
    return NULL;
  }

  atp_se = AsnFind ("Seq-entry");
  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry");
    return NULL;
  }

  atp_ses = AsnFind ("Seq-entry.set");
  if (atp_ses == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry.set");
    return NULL;
  }

  atp_ssdee = AsnFind ("Seq-submit.data.entrys.E");
  if (atp_ssdee == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys.E");
    return NULL;
  }


  ssp_atp = AsnLinkType (NULL, atp_ss);
  if (ssp_atp == NULL) return NULL;

  MemSet ((Pointer) &ss, 0, sizeof (SeqSubmit));
  MemSet ((Pointer) &se, 0, sizeof (SeqEntry));
  se.choice = 2;

  if (! AsnOpenStruct (aip, ssp_atp, (Pointer) &ss)) return NULL;

  if (! SubmitBlockAsnWrite (sbp, aip, atp_sss)) return NULL;

  av.ptrvalue = (Pointer) &se;
  if (! AsnWriteChoice (aip, atp_ssd, (Int2) 1, &av)) return NULL;

  if (! AsnOpenStruct (aip, atp_ssde, (Pointer) &se)) return NULL;

  sep_atp = AsnLinkType (atp_ssdee, atp_se);
  if (sep_atp == NULL) return NULL;

  av.ptrvalue = (Pointer) &se;
  se.choice = 2;
  if (! AsnWriteChoice (aip, sep_atp, (Int2) 2, &av)) return NULL;

  return ssp_atp;
}

static AsnTypePtr DoSecondPrefix (
  AsnIoPtr aip,
  TblArgsPtr tbl
)

{
  AsnTypePtr  atp_bsc, atp_bss, atp_bsss, atp_ses, bssp_atp;
  DataVal     av;
  BioseqSet   bs;

  atp_ses = AsnFind ("Seq-entry.set");
  if (atp_ses == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry.set");
    return NULL;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set");
    return NULL;
  }

  atp_bsc = AsnFind ("Bioseq-set.class");
  if (atp_bsc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.class");
    return NULL;
  }

  atp_bsss = AsnFind ("Bioseq-set.seq-set");
  if (atp_bsss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set");
    return NULL;
  }


  bssp_atp = AsnLinkType (atp_ses, atp_bss);
  if (bssp_atp == NULL) return NULL;

  MemSet ((Pointer) &bs, 0, sizeof (BioseqSet));

  if (! AsnOpenStruct (aip, bssp_atp, (Pointer) &bs)) return NULL;

  switch (tbl->whichclass) {
    case 1 :
      av.intvalue = BioseqseqSet_class_pop_set;
      break;
    case 2 :
      av.intvalue = BioseqseqSet_class_phy_set;
      break;
    case 3 :
      av.intvalue = BioseqseqSet_class_mut_set;
      break;
    case 4 :
      av.intvalue = BioseqseqSet_class_eco_set;
      break;
    default :
      av.intvalue = BioseqseqSet_class_genbank;
      break;
  }
  if (! AsnWrite (aip, atp_bsc, &av)) return NULL;

  if (! AsnOpenStruct (aip, atp_bsss, (Pointer) &bs.seq_set)) return NULL;

  return bssp_atp;
}

static Boolean DoFirstSuffix (
  AsnIoPtr aip,
  AsnTypePtr ssp_atp
)

{
  AsnTypePtr  atp_bsss, atp_ssde, atp_ssdee;
  BioseqSet   bs;
  SeqEntry    se;
  SeqSubmit   ss;

  if (aip == NULL || ssp_atp == NULL) return FALSE;

  atp_ssde = AsnFind ("Seq-submit.data.entrys");
  if (atp_ssde == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys");
    return FALSE;
  }

  atp_ssdee = AsnFind ("Seq-submit.data.entrys.E");
  if (atp_ssdee == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys.E");
    return FALSE;
  }

  atp_bsss = AsnFind ("Bioseq-set.seq-set");
  if (atp_bsss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set");
    return FALSE;
  }


  MemSet ((Pointer) &ss, 0, sizeof (SeqSubmit));
  MemSet ((Pointer) &se, 0, sizeof (SeqEntry));
  MemSet ((Pointer) &bs, 0, sizeof (BioseqSet));

  if (! AsnCloseStruct (aip, atp_ssde, &se)) return FALSE;

  if (! AsnCloseStruct (aip, ssp_atp, (Pointer) &ss)) return FALSE;

  AsnUnlinkType (atp_ssdee);

  return TRUE;
}

static Boolean DoSecondSuffix (
  AsnIoPtr aip,
  AsnTypePtr bssp_atp
)

{
  AsnTypePtr  atp_bsss, atp_ses;
  BioseqSet   bs;

   if (aip == NULL || bssp_atp == NULL) return FALSE;

  atp_ses = AsnFind ("Seq-entry.set");
  if (atp_ses == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry.set");
    return FALSE;
  }

  atp_bsss = AsnFind ("Bioseq-set.seq-set");
  if (atp_bsss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set");
    return FALSE;
  }


  MemSet ((Pointer) &bs, 0, sizeof (BioseqSet));

  if (! AsnCloseStruct(aip, atp_bsss, (Pointer) &bs.seq_set)) return FALSE;

  if (! AsnCloseStruct (aip, bssp_atp, (Pointer) &bs)) return FALSE;

  AsnUnlinkType (atp_ses);

  return TRUE;
}

static CharPtr ReadCommentFile (
  CharPtr filename
)

{
  FileCache   fc;
  FILE        *fp;
  ValNodePtr  head = NULL, last = NULL, vnp;
  Int4        len;
  Char        line [4096];
  Boolean     nonewline, notfirst;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (filename)) return NULL;
  fp = FileOpen (filename, "r");
  if (fp == NULL) return NULL;

  FileCacheSetup (&fc, fp);

  str = FileCacheReadLine (&fc, line, sizeof (line), &nonewline);
  while (str != NULL) {
    vnp = ValNodeCopyStr (&last, 0, str);
    if (head == NULL) {
      head = vnp;
    }
    last = vnp;

    str = FileCacheReadLine (&fc, line, sizeof (line), &nonewline);
  }

  FileClose (fp);

  if (head == NULL) return NULL;

  len = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    len += StringLen (str) + 1;
  }

  tmp = (CharPtr) MemNew (sizeof (Char) * (len + 5));
  if (tmp == NULL) return NULL;

  ptr = tmp;
  notfirst = FALSE;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str == NULL) continue;
    if (*str == '\0' || *str == ' ') {
      ptr = StringMove (ptr, "~");
    } else if (notfirst) {
      ptr = StringMove (ptr, " ");
    }
    ptr = StringMove (ptr, str);
    notfirst = TRUE;
  }

  ValNodeFreeData (head);

  return tmp;
}

static CharPtr ParseCommaField (
  CharPtr PNTR strP
)

{
  CharPtr  ptr;
  CharPtr  str;

  if (strP == NULL) return NULL;

  str = *strP;
  if (StringHasNoText (str)) {
    *strP = NULL;
    return NULL;
  }

  ptr = StringChr (str, ',');
  if (ptr == NULL) {
    *strP = NULL;
    return str;
  }

  *ptr = '\0';
  ptr++;
  if (StringHasNoText (ptr)) {
    ptr = NULL;
  }
  *strP = ptr;

  if (StringHasNoText (str)) {
    str = NULL;
  }
  return str;
}

static DatePtr DateParse (
  CharPtr str
)

{
  Int4      day = -1, month = -1, year = -1;
  DatePtr   dp;
  CharPtr   ptr;
  Char      tmp [64];
  long int  val;

  if (StringHasNoText (str)) return NULL;

  StringNCpy_0 (tmp, str, sizeof (tmp));
  ptr = StringChr (tmp, '/');
  if (ptr == NULL) {
    ptr = StringChr (tmp, '-');
  }
  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
    if (sscanf (tmp, "%ld", &val) == 1) {
      month = (Int4) val;
    }
    str = StringChr (ptr, '/');
    if (str == NULL) {
      str = StringChr (ptr, '-');
    }
    if (str != NULL) {
      *str = '\0';
      str++;
      if (sscanf (ptr, "%ld", &val) == 1) {
        day = (Int4) val;
      }
      if (sscanf (str, "%ld", &val) == 1) {
        year = (Int4) val;
     }
    }
  }

  if (month < 0 || day < 0 || year < 2000) return NULL;
  if (month > 12 || day > 31 || year > 2099) return NULL;

  dp = DateNew ();
  if (dp == NULL) return NULL;

  dp->data [0] = 1;
  dp->data [1] = (Uint1) (year - 1900);
  dp->data [2] = (Uint1) month;
  dp->data [3] = (Uint1) day;

  return dp;
}

/* Args structure contains command-line arguments */

#define p_argInputPath         0
#define r_argOutputPath        1
#define i_argInputFile         2
#define o_argOutputFile        3
#define x_argSuffix            4
#define E_argRecurse           5
#define t_argTemplate          6
#define a_argType              7
#define s_argFastaSet          8
#define g_argGenProdSet        9
#define F_argFeatIdLinks      10
#define A_argAccession        11
#define C_argCenter           12
#define n_argOrgName          13
#define j_argSrcQuals         14
#define y_argComment          15
#define Y_argCommentFile      16
#define D_argDescrsFile       17
#define f_argTableFile        18
#define k_argCdsFlags         19
#define V_argVerify           20
#define v_argValidate         21
#define b_argGenBank          22
#define q_argFileID           23
#define u_argUndoGPS          24
#define h_argGnlToNote        25
#define G_argGapFields        26
#define R_argRemote           27
#define S_argSmartFeats       28
#define Q_argSmartTitle       29
#define U_argUnnecXref        30
#define L_argLocalID          31
#define T_argTaxLookup        32
#define P_argPubLookup        33
#define W_argLogProgress      34
#define K_argBioseqSet        35
#define H_argHoldUntilPub     36
#define Z_argDiscRepFile      37
#define c_argCleanupOptions   38
#define X_argExtraFlags       39


Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Suffix", ".fsa", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'E', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Template File", NULL, NULL, NULL,
    TRUE, 't', ARG_FILE_IN, 0.0, 0, NULL},
  {"File Type\n"
   "      a Any\n"
   "      r20u Runs of 20+ Ns are gaps, 100 Ns are unknown length\n"
   "      r20k Runs of 20+ Ns are gaps, 100 Ns are known length\n"
   "      s FASTA Set (s Batch, s1 Pop, s2 Phy, s3 Mut, s4 Eco)\n"
   "      d FASTA Delta, di FASTA Delta with Implicit Gaps\n"
   "      l FASTA+Gap Alignment\n"
   "      z FASTA with Gap Lines\n"
   "      e PHRAP/ACE\n", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Read FASTAs as Set", "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Genomic Product Set", "F", NULL, NULL,
    TRUE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Feature ID Links (o by Overlap, p by Product)", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Accession", NULL, NULL, NULL,
    TRUE, 'A', ARG_STRING, 0.0, 0, NULL},
  {"Genome Center Tag", NULL, NULL, NULL,
    TRUE, 'C', ARG_STRING, 0.0, 0, NULL},
  {"Organism Name", NULL, NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Source Qualifiers", NULL, NULL, NULL,
    TRUE, 'j', ARG_STRING, 0.0, 0, NULL},
  {"Comment", NULL, NULL, NULL,
    TRUE, 'y', ARG_STRING, 0.0, 0, NULL},
  {"Comment File", NULL, NULL, NULL,
    TRUE, 'Y', ARG_FILE_IN, 0.0, 0, NULL},
  {"Descriptors File", NULL, NULL, NULL,
    TRUE, 'D', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Table File", NULL, NULL, NULL,
    TRUE, 'f', ARG_FILE_IN, 0.0, 0, NULL},
  {"CDS Flags (combine any of the following letters)\n"
   "      c Annotate Longest ORF\n"
   "      r Allow Runon ORFs\n"
   "      m Allow Alternative Starts\n"
   "      k Set Conflict on Mismatch\n", NULL, NULL, NULL,
    TRUE, 'k', ARG_STRING, 0.0, 0, NULL},
  {"Verification (combine any of the following letters)\n"
   "      v Validate with Normal Stringency\n"
   "      r Validate without Country Check\n"
   "      b Generate GenBank Flatfile\n"
   "      g Generate Gene Report\n", NULL, NULL, NULL,
    TRUE, 'V', ARG_STRING, 0.0, 0, NULL},
  {"Validate (obsolete: use -V v)", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Generate GenBank File (obsolete: use -V b)", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Seq ID from File Name", "F", NULL, NULL,
    TRUE, 'q', ARG_BOOLEAN, 0.0, 0, NULL},
  {"GenProdSet to NucProtSet", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"General ID to Note", "F", NULL, NULL,
    TRUE, 'h', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Alignment Gap Flags (comma separated fields, e.g., p,-,-,-,?,. )\n"
   "      n Nucleotide or p Protein,\n"
   "      Begin, Middle, End Gap Characters,\n"
   "      Missing Characters, Match Characters\n",  NULL, NULL, NULL,
    TRUE, 'G', ARG_STRING, 0.0, 0, NULL},
  {"Remote Sequence Record Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Smart Feature Annotation", "F", NULL, NULL,
    TRUE, 'S', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Special mRNA Titles", "F", NULL, NULL,
    TRUE, 'Q', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remove Unnecessary Gene Xref", "F", NULL, NULL,
    TRUE, 'U', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Force Local protein_id/transcript_id", "F", NULL, NULL,
    TRUE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Publication Lookup", "F", NULL, NULL,
    TRUE, 'P', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log Progress", "F", NULL, NULL,
    TRUE, 'W', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Save Bioseq-set", "F", NULL, NULL,
    TRUE, 'K', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Hold Until Publish\n"
   "      y Hold for One Year\n"
   "      mm/dd/yyyy\n", NULL, NULL, NULL,
    TRUE, 'H', ARG_STRING, 0.0, 0, NULL},
  {"Discrepancy Report Output File", NULL, NULL, NULL,
    TRUE, 'Z', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Cleanup (combine any of the following letters)\n"
   "      d Correct Collection Dates (assume month first)\n"
   "      D Correct Collection Dates (assume day first)\n"
   "      b Append note to coding regions that overlap other coding regions with similar product names and do not contain 'ABC'",
    NULL, NULL, NULL,
    TRUE, 'c', ARG_STRING, 0.0, 0, NULL},
  {"Extra Flags (combine any of the following letters)\n"
   "      C Apply comments in .cmt files to all sequences\n",  NULL, NULL, NULL,
    TRUE, 'X', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  AsnIoPtr        aip = NULL;
  Char            app [64];
  CharPtr         base;
  AsnTypePtr      bssp_atp = NULL;
  CitSubPtr       csp;
  Pointer         dataptr;
  Uint2           datatype;
  CharPtr         descrs;
  CharPtr         directory;
  DatePtr         dp;
  FILE            *fp;
  Char            gapstring [128];
  TextFsaPtr      gotags;
  CharPtr         hold;
  CharPtr         os;
  CharPtr         outfile;
  Pubdesc         pd;
  PubdescPtr      pdp = NULL;
  ValNode         pb;
  CharPtr         ptr;
  Boolean         recurse;
  Boolean         remote;
  CharPtr         results;
  SubmitBlockPtr  sbp = NULL;
  SeqDescrPtr     sdphead = NULL;
  SeqEntryPtr     sep;
  Char            sfx [32];
  BioSourcePtr    src = NULL;
  SeqSubmitPtr    ssp = NULL;
  AsnTypePtr      ssp_atp = NULL;
  Char            str [64];
  CharPtr         suffix;
  TblArgs         tbl;
  CharPtr         tmp;
  CharPtr         tmplate;
  CharPtr         disc_rep_file = NULL;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

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

  sprintf (app, "tbl2asn %s", TBL2ASN_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = NULL;
  }
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  recurse = (Boolean) myargs [E_argRecurse].intvalue;
  base = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  if (StringHasNoText (outfile)) {
    outfile = NULL;
  }
  tmplate = (CharPtr) myargs [t_argTemplate].strvalue;
  descrs = (CharPtr) myargs [D_argDescrsFile].strvalue;

  hold = (CharPtr) myargs [H_argHoldUntilPub].strvalue;

  if (StringHasNoText(directory) && StringHasNoText(base)) {
    Message (MSG_FATAL, "You must supply either an input file (-i) or an input directory (-p).\nUse -p . to specify the current directory.\n\n");
    return 1;
  }
  remote = (Boolean) myargs [R_argRemote].intvalue;

  MemSet ((Pointer) &tbl, 0, sizeof (TblArgs));

  /* -s is heavily used and will remain as an alternative to -a s */

  tbl.fastaset = (Boolean) myargs [s_argFastaSet].intvalue;

  /* process new -a type argument */

  ptr = myargs [a_argType].strvalue;
  if (StringICmp (ptr, "r20u") == 0) {
    tbl.raw2delt = TRUE;
    tbl.r2dmin = 20;
    tbl.r2dunk100 = TRUE;
  } else if (StringICmp (ptr, "r20k") == 0) {
    tbl.raw2delt = TRUE;
    tbl.r2dmin = 20;
    tbl.r2dunk100 = FALSE;
  } else if (StringICmp (ptr, "s") == 0) {
    tbl.fastaset = TRUE;
  } else if (StringICmp (ptr, "w1") == 0 || StringICmp (ptr, "s1") == 0) {
    tbl.fastaset = TRUE;
    tbl.whichclass = 1;
  } else if (StringICmp (ptr, "w2") == 0 || StringICmp (ptr, "s2") == 0) {
    tbl.fastaset = TRUE;
    tbl.whichclass = 2;
  } else if (StringICmp (ptr, "w3") == 0 || StringICmp (ptr, "s3") == 0) {
    tbl.fastaset = TRUE;
    tbl.whichclass = 3;
  } else if (StringICmp (ptr, "w4") == 0 || StringICmp (ptr, "s4") == 0) {
    tbl.fastaset = TRUE;
    tbl.whichclass = 4;
  } else if (StringICmp (ptr, "d") == 0) {
    tbl.deltaset = TRUE;
  } else if (StringICmp (ptr, "di") == 0) {
    tbl.deltaset = TRUE;
    tbl.implicitgaps = TRUE;
  } else if (StringICmp (ptr, "l") == 0) {
    tbl.alignset = TRUE;
  } else if (StringICmp (ptr, "z") == 0) {
    tbl.gapped = TRUE;
  } else if (StringICmp (ptr, "e") == 0) {
    tbl.phrapace = TRUE;
  }

  tbl.genprodset = (Boolean) myargs [g_argGenProdSet].intvalue;
  ptr = myargs [F_argFeatIdLinks].strvalue;
  if (StringICmp (ptr, "o") == 0) {
    tbl.linkbyoverlap = TRUE;
  } else if (StringICmp (ptr, "p") == 0) {
    tbl.linkbyproduct = TRUE;
  }
  tbl.forcelocalid = (Boolean) myargs [L_argLocalID].intvalue;
  tbl.gpstonps = (Boolean) myargs [u_argUndoGPS].intvalue;
  tbl.gnltonote = (Boolean) myargs [h_argGnlToNote].intvalue;
  tbl.accn = (CharPtr) myargs [A_argAccession].strvalue;
  tbl.center = (CharPtr) myargs [C_argCenter].strvalue;
  tbl.organism = (CharPtr) myargs [n_argOrgName].strvalue;
  tbl.srcquals = (CharPtr) myargs [j_argSrcQuals].strvalue;
  tbl.comment = (CharPtr) myargs [y_argComment].strvalue;
  tbl.commentFile = ReadCommentFile ((CharPtr) myargs [Y_argCommentFile].strvalue);

  ptr = myargs [k_argCdsFlags].strvalue;
  if (StringChr (ptr, 'c') != NULL) {
    tbl.findorf = TRUE;
  }
  if (StringChr (ptr, 'r') != NULL) {
    tbl.runonorf = TRUE;
    tbl.findorf = TRUE;
  }
  if (StringChr (ptr, 'm') != NULL) {
    tbl.altstart = TRUE;
  }
  if (StringChr (ptr, 'k') != NULL) {
    tbl.conflict = TRUE;
  }
  /*
  if (!tbl.findorf && tbl.runonorf) {
    Message (MSG_FATAL, "-k r cannot be used without -k c");
    return 1;
  }
  */

  /* process obsolete validate/flatfile arguments first, warn if used */

  tbl.validate = (Boolean) myargs [v_argValidate].intvalue;
  if (tbl.validate) {
    Message (MSG_POST, "-v is obsolete, use -V v instead");
  }
  tbl.flatfile = (Boolean) myargs [b_argGenBank].intvalue;
  if (tbl.flatfile) {
    Message (MSG_POST, "-b is obsolete, use -V b instead");
  }

  ptr = myargs [V_argVerify].strvalue;
  if (StringChr (ptr, 'v') != NULL) {
    tbl.validate = TRUE;
  }
  if (StringChr (ptr, 'r') != NULL) {
    tbl.validate = TRUE;
    tbl.relaxed = TRUE;
  }
  if (StringChr (ptr, 'b') != NULL) {
    tbl.flatfile = TRUE;
  }
  if (StringChr (ptr, 'g') != NULL) {
    tbl.genereport = TRUE;
  }
  if (StringChr (ptr, 'c') != NULL) {
    tbl.validate_barcode = TRUE;
  }
  

  tbl.seqidfromfile = (Boolean) myargs [q_argFileID].intvalue;
  tbl.smartfeats = (Boolean) myargs [S_argSmartFeats].intvalue;
  tbl.smarttitle = (Boolean) myargs [Q_argSmartTitle].intvalue;
  tbl.removeunnecxref = (Boolean) myargs [U_argUnnecXref].intvalue;
  tbl.dotaxlookup = (Boolean) myargs [T_argTaxLookup].intvalue;
  tbl.dopublookup = (Boolean) myargs [P_argPubLookup].intvalue;
  tbl.logtoterminal = (Boolean) myargs [W_argLogProgress].intvalue;

  tbl.save_bioseq_set = (Boolean) myargs [K_argBioseqSet].intvalue;

  /* get extra flags */
  ptr = myargs [X_argExtraFlags].strvalue;
  if (StringChr (ptr, 'C') == NULL) {
    tbl.apply_cmt_to_all = FALSE;
  } else {
    tbl.apply_cmt_to_all = TRUE;
  }

  disc_rep_file = (CharPtr) myargs [Z_argDiscRepFile].strvalue;
  if (StringHasNoText (disc_rep_file)) {
    tbl.global_report = NULL;
  } else {
    tbl.global_report = GlobalDiscrepReportNew();
    tbl.global_report->test_config = DiscrepancyConfigNew ();
    DisableTRNATests (tbl.global_report->test_config);
    ConfigureForGenomes (tbl.global_report->test_config);
    tbl.global_report->taxlookup = taxlookup;
    tbl.global_report->output_config->summary_report = FALSE;
    tbl.global_report->output_config->expand_report_categories[DISC_SUPERFLUOUS_GENE] = TRUE;
    tbl.global_report->output_config->expand_report_categories[DISC_RNA_CDS_OVERLAP] = TRUE;
    tbl.global_report->output_config->expand_report_categories[DISC_SUSPECT_PRODUCT_NAME] = TRUE;
    tbl.global_report->output_config->expand_report_categories[DISC_OVERLAPPING_CDS] = TRUE;
  }


  /* arguments for alignment reading, e.g., "p,-,-,-,?,." */

  gapstring [0] = '\0';
  ptr = (CharPtr) myargs [G_argGapFields].strvalue;
  StringNCpy_0 (gapstring, ptr, sizeof (gapstring));

  ptr = gapstring;
  tmp = ParseCommaField (&ptr);
  if (tmp != NULL) {
    if (StringChr (tmp, 'p') != NULL) {
      tbl.aln_is_protein = TRUE;
    } else if (StringChr (tmp, 'n') == NULL) {
      Message (MSG_FATAL, "-G must start with p for Protein or n for Nucleotide");
      return 1;
    }
  }
  tbl.aln_beginning_gap = ParseCommaField (&ptr);
  tbl.aln_middle_gap = ParseCommaField (&ptr);
  tbl.aln_end_gap = ParseCommaField (&ptr);
  tbl.aln_missing = ParseCommaField (&ptr);
  tbl.aln_match = ParseCommaField (&ptr);

  if (StringHasNoText (tbl.accn)) {
    tbl.accn = NULL;
  }
  if (StringHasNoText (tbl.organism)) {
    tbl.organism = NULL;
  }
  if (StringHasNoText (tbl.srcquals)) {
    tbl.srcquals = NULL;
  }
  if (StringHasNoText (tbl.comment)) {
    tbl.comment = NULL;
  }
  if (StringHasNoText (tbl.commentFile)) {
    tbl.commentFile = NULL;
  }

  if (tbl.fastaset &&
      (tbl.deltaset || tbl.phrapace || tbl.genprodset ||
       tbl.alignset || tbl.gapped)) {
    Message (MSG_FATAL, "-s cannot be used with -d, -e, -g, -l or -z");
    return 1;
  }

  if (! tbl.alignset && (StringDoesHaveText (tbl.aln_beginning_gap)
      || StringDoesHaveText (tbl.aln_end_gap)
      || StringDoesHaveText (tbl.aln_middle_gap)
      || StringDoesHaveText (tbl.aln_missing)
      || StringDoesHaveText (tbl.aln_match)
      || tbl.aln_is_protein)) {
    Message (MSG_FATAL, "-G can only be used with -a l");
    return 1;
  }

  /* arguments for cleanup */
  MemSet (&(tbl.cleanup_args), 0, sizeof (CleanupArgsData));
  ptr = (CharPtr) myargs [c_argCleanupOptions].strvalue;
  if (StringChr (ptr, 'd') != NULL) {
    if (StringChr (ptr, 'D') != NULL) {
      Message (MSG_FATAL, "Cannot use both d and D options for cleanup.  Choose one.");
      return 1;
    }
    tbl.cleanup_args.collection_dates = TRUE;
    tbl.cleanup_args.collection_dates_month_first = TRUE;
  } else if (StringChr (ptr, 'D') != NULL) {
    tbl.cleanup_args.collection_dates = TRUE;
    tbl.cleanup_args.collection_dates_month_first = FALSE;
  }

  if (StringChr (ptr, 'b') != NULL) {
    tbl.cleanup_args.add_notes_to_overlapping_cds_without_abc = TRUE;
  }
  
  if (StringHasNoText (base) && (StringDoesHaveText (tbl.accn))) {
    Message (MSG_FATAL, "Accession can be entered only for a single record");
    return 1;
  }

  /* Seq-submit or Submit-block template is optional */

  if (StringDoesHaveText (tmplate)) {
    if (TemplateOverwriteRisk (tmplate, base, directory, suffix)) {
      if (Message (MSG_YN, overwriteMsg) == ANS_NO) return 0;
    }
    fp = FileOpen (tmplate, "r");
    if (fp != NULL) {
      while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
        if (datatype == OBJ_SEQSUB) {
          ssp = (SeqSubmitPtr) dataptr;
        } else if (datatype == OBJ_SUBMIT_BLOCK) {
          sbp = (SubmitBlockPtr) dataptr;
        } else if (datatype == OBJ_SEQDESC) {
          ValNodeLink (&sdphead, (SeqDescrPtr) dataptr);
        } else {
          ObjMgrFree (datatype, dataptr);
        }
      }
      FileClose (fp);
    }

    if (ssp != NULL && sbp == NULL) {
      sbp = ssp->sub;
    }
    if (sbp == NULL) {
      Message (MSG_FATAL, "Unable to read required template file");
      return 1;
    }

    if (sbp != NULL) {
      if (ssp != NULL) {

        /* copy submit block, will free SeqSubmit before processing */

        sbp = AsnIoMemCopy ((Pointer) sbp,
                            (AsnReadFunc) SubmitBlockAsnRead,
                            (AsnWriteFunc) SubmitBlockAsnWrite);
      }
      sbp->tool = MemFree (sbp->tool);
      os = GetOpSysString ();
      if (os != NULL) {
        sprintf (str, "tbl2asn %s - %s", TBL2ASN_APPLICATION, os);
      } else {
        sprintf (str, "tbl2asn %s", TBL2ASN_APPLICATION);
      }
      sbp->tool = StringSave (str);
      MemFree (os);
      sbp->hup = FALSE;
      sbp->reldate = DateFree (sbp->reldate);
      if (StringDoesHaveText (hold)) {
        if (StringICmp (hold, "y") == 0) {
          sbp->hup = TRUE;
          dp = DateCurr ();
          sbp->reldate = dp;
          if (dp != NULL) {
            if (dp->data [0] == 1) {
              (dp->data [1])++;
            }
          }
        } else {
          dp = DateParse (hold);
          if (dp != NULL) {
            sbp->hup = TRUE;
            sbp->reldate = dp;
          }
        }
      }
      csp = sbp->cit;
      if (csp != NULL) {
        csp->date = DateFree (csp->date);
        csp->date = DateCurr ();
        MemSet ((Pointer) &pd, 0, sizeof (Pubdesc));
        MemSet ((Pointer) &pb, 0, sizeof (ValNode));
        pb.choice = PUB_Sub;
        pb.data.ptrvalue = (Pointer) csp;
        pd.pub = &pb;
        pdp = &pd;
      }
    }
    if (ssp != NULL && ssp->datatype == 1) {
      sep = (SeqEntryPtr) ssp->data;
      if (sep != NULL) {
        VisitBioSourcesInSep (sep, (Pointer) &src, GetFirstBiop);
        if (src != NULL) {

          /* copy top biosource */

          src = AsnIoMemCopy ((Pointer) src,
                              (AsnReadFunc) BioSourceAsnRead,
                              (AsnWriteFunc) BioSourceAsnWrite);
        }
      }

      /* in case template has colliding ID, free it now */

      SeqSubmitFree (ssp);
    }
  }

  if (StringDoesHaveText (descrs)) {
    if (TemplateOverwriteRisk (descrs, base, directory, suffix)) {
      if (Message (MSG_YN, overwriteMsg) == ANS_NO) return 0;
    }
    fp = FileOpen (descrs, "r");
    if (fp != NULL) {
      while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
        if (datatype == OBJ_SEQDESC) {
          ValNodeLink (&sdphead, (SeqDescrPtr) dataptr);
        } else {
          ObjMgrFree (datatype, dataptr);
        }
      }
      FileClose (fp);
    }
  }

  gotags = TextFsaNew ();
  TextFsaAdd (gotags, "go_component");
  TextFsaAdd (gotags, "go_function");
  TextFsaAdd (gotags, "go_process");

  /* register fetch functions */

  if (remote) {
#ifdef INTERNAL_NCBI_TBL2ASN
    if (! PUBSEQBioseqFetchEnable ("tbl2asn", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
#else
    PubSeqFetchEnable ();
#endif
  }

  if (remote || tbl.dopublookup) {
    PubMedFetchEnable ();
  }

  /* process one or more records */

  if (StringDoesHaveText (outfile) && StringHasNoText (base)) {
    aip = AsnIoOpen (outfile, "w");
    if (aip == NULL) {
      Message (MSG_FATAL, "Unable to open single output file");
      return 1;
    }
    ssp_atp = DoFirstPrefix (aip, sbp);
    bssp_atp = DoSecondPrefix (aip, &tbl);
  }

  if (StringDoesHaveText (base)) {
    ptr = StringRChr (base, '.');
    sfx[0] = '\0';
    if (ptr != NULL) {
      StringNCpy_0 (sfx, ptr, sizeof (sfx));
      *ptr = '\0';
    }
    tbl.tableFile = (CharPtr) myargs [f_argTableFile].strvalue;
    ProcessOneRecord (sbp, pdp, src, directory, results, base, sfx, sdphead, &tbl, gotags, aip, outfile);

  } else {

    FileRecurse (sbp, pdp, src, directory, results, suffix, recurse, sdphead, &tbl, gotags, aip, NULL);
  }

  if (aip != NULL) {
    DoSecondSuffix (aip, bssp_atp);
    DoFirstSuffix (aip, ssp_atp);
    AsnIoClose (aip);
  }

  if (tbl.global_report != NULL) {
    fp = FileOpen (disc_rep_file, "w");
    WriteGlobalDiscrepancyReport (tbl.global_report, fp);
    FileClose (fp);
    tbl.global_report = GlobalDiscrepReportFree (tbl.global_report);
  }

  if (sbp != NULL) {
    SubmitBlockFree (sbp);
  }
  if (src != NULL) {
    BioSourceFree (src);
  }

  SeqDescrFree (sdphead);

  TransTableFreeAll ();

  ECNumberFSAFreeAll ();

  TextFsaFree (gotags);

  /* close fetch functions */

  if (remote || tbl.dopublookup) {
    PubMedFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_TBL2ASN
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return 0;
}

