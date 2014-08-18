/*   cspeedtest.c
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
* File Name:  cspeedtest.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   12/17/07
*
* $Revision: 1.20 $
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
#include <objfdef.h>
#include <objsub.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <toasn3.h>
#include <pmfapi.h>
#include <tofasta.h>
#include <asn2gnbk.h>
#include <valid.h>
#include <suggslp.h>

NLM_EXTERN CharPtr NewCreateDefLine (
  ItemInfoPtr iip,
  BioseqPtr bsp,
  Boolean ignoreTitle,
  Boolean extProtTitle
);

#define CSPEEDTEST_APP_VER "1.9"

CharPtr CSPEEDTEST_APPLICATION = CSPEEDTEST_APP_VER;

typedef struct cspeedflags {
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Boolean       lock;
  Int2          type;
  Int4          maxcount;
  CharPtr       io;
  CharPtr       clean;
  CharPtr       skip;
  CharPtr       index;
  CharPtr       seq;
  CharPtr       feat;
  CharPtr       desc;
  CharPtr       verify;
  BioseqPtr     nucbsp;
  Int2          genCode;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  BioseqSet     bss;
  FILE          *ofp;
  FILE          *logfp;
} CSpeedFlagData, PNTR CSpeedFlagPtr;

static void DoVisitFeaturesTest (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  /* empty visit callback */
}

static void DoVisitCodingRegions (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr      bsp;
  CharPtr        caret5, caret3;
  CSpeedFlagPtr  cfp;
  Char           id [64];
  SeqLocPtr      loc, slp;
  Boolean        partial5, partial3;
  SeqIdPtr       sip;
  Int4           start, stop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL || cfp->ofp == NULL) return;

  loc = sfp->location;
  bsp = BioseqFindFromSeqLoc (loc);
  if (bsp == NULL) return;

  StringCpy (id, "?");
  if (sfp->product != NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL) {
      SeqIdWrite (sip, id, PRINTID_FASTA_SHORT, sizeof (id) - 1);
    }
  }

  fprintf (cfp->ofp, "%s\n", id);
  slp = SeqLocFindNext (loc, NULL);
  while (slp != NULL) {
    start = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) + 1;
    caret5 = "";
    caret3 = "";
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (partial5) {
      caret5 = "<";
    }
    if (partial3) {
      caret3 = ">";
    }
    fprintf (cfp->ofp, "%s%ld\t%s%ld\n", caret5, (long) start, caret3, (long) stop);
    slp = SeqLocFindNext (loc, slp);
  }
}

static void DoSuggestIntervals (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr        caret5, caret3;
  CSpeedFlagPtr  cfp;
  Char           id [64];
  SeqLocPtr      loc, slp;
  Boolean        partial5, partial3;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;
  Int4           start, stop;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL || cfp->ofp == NULL || cfp->nucbsp == NULL) return;

  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  SeqIdWrite (sip, id, PRINTID_FASTA_SHORT, sizeof (id) - 1);

  sap = SuggestCodingRegion (cfp->nucbsp, bsp, cfp->genCode);
  if (sap == NULL) return;
  if (sap->type == 1) {
    sfp = (SeqFeatPtr) sap->data;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      loc = sfp->location;
      if (loc != NULL) {
        fprintf (cfp->ofp, "%s\n", id);
        slp = SeqLocFindNext (loc, NULL);
        while (slp != NULL) {
          start = GetOffsetInBioseq (slp, cfp->nucbsp, SEQLOC_START) + 1;
          stop = GetOffsetInBioseq (slp, cfp->nucbsp, SEQLOC_STOP) + 1;
          caret5 = "";
          caret3 = "";
          CheckSeqLocForPartial (slp, &partial5, &partial3);
          if (partial5) {
            caret5 = "<";
          }
          if (partial3) {
            caret3 = ">";
          }
          fprintf (cfp->ofp, "%s%ld\t%s%ld\n", caret5, (long) start, caret3, (long) stop);
          slp = SeqLocFindNext (loc, slp);
        }
      }
    }
  }
  SeqAnnotFree (sap);
}

static void DoGeneOverlapPrintTest (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CSpeedFlagPtr      cfp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;
  CharPtr            str1, str2;

  if (sfp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL || cfp->ofp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE) return;
  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;

  str1 = SeqLocPrint (sfp->location);
  str2 = SeqLocPrint (gene->location);
  if (str1 != NULL && str2 != NULL) {
    fprintf (cfp->ofp, "[%s] -> [%s]\n", str1, str2);
  } else {
    fprintf (cfp->ofp, "? -> ?\n");
  }
  MemFree (str1);
  MemFree (str2);
}

static void DoGeneOverlapSpeedTest (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CSpeedFlagPtr      cfp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;

  if (sfp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL || cfp->ofp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE) return;
  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;
}

static void LIBCALLBACK EmptyStreamProc (
  CharPtr sequence,
  Pointer userdata
)

{
  /* empty stream callback */
}

static void DoFastaSeq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CSpeedFlagPtr  cfp;

  if (bsp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (cfp->ofp != NULL) {
    BioseqFastaStream (bsp, cfp->ofp, STREAM_EXPAND_GAPS, 70, 0, 0, TRUE);
  } else {
    SeqPortStream (bsp, STREAM_EXPAND_GAPS, NULL, EmptyStreamProc);
  }
}

static void DoFastaRaw (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CSpeedFlagPtr  cfp;

  if (bsp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (cfp->ofp != NULL) {
    fprintf (cfp->ofp, ">\n");
    BioseqFastaStream (bsp, cfp->ofp, STREAM_EXPAND_GAPS, 70, 0, 0, FALSE);
  } else {
    SeqPortStream (bsp, STREAM_EXPAND_GAPS, NULL, EmptyStreamProc);
  }
}

static void DoFastaDefline (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Char           buf [4096];
  CSpeedFlagPtr  cfp;
  Char           id [128];

  if (bsp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  id [0] = '\0';
  SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  buf [0] = '\0';
  CreateDefLine (NULL, bsp, buf, sizeof (buf) - 1, 0, NULL, NULL);

  if (cfp->ofp != NULL) {
    fprintf (cfp->ofp, ">%s %s\n", id, buf);
  }
}

static void DoNewFastaDefline (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqSetPtr   bssp;
  CSpeedFlagPtr  cfp;
  Char           id [128];
  CharPtr        title;

  if (bsp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (StringChr (cfp->skip, 's') != NULL) {
    if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bsp->idx.parentptr;
      if (bssp != NULL) {
        if (bssp->_class == BioseqseqSet_class_segset ||
            bssp->_class == BioseqseqSet_class_parts) return;
      }
    }
  }
  if (StringChr (cfp->skip, 'v') != NULL) {
    if (bsp->repr == Seq_repr_virtual) return;
  }

  id [0] = '\0';
  SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  title = NewCreateDefLine (NULL, bsp, FALSE, FALSE);
  if (StringHasNoText (title)) {
    title = StringSave ("?");
  }

  if (cfp->ofp != NULL) {
    fprintf (cfp->ofp, ">%s %s\n", id, title);
  }

  MemFree (title);
}

static void DoFastaComp (
  BioseqPtr bsp,
  Pointer userdata,
  Boolean ignoreExisting
)

{
  Char           buf [4096];
  CSpeedFlagPtr  cfp;
  Char           id [128];
  CharPtr        title;

  if (bsp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  id [0] = '\0';
  SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
  buf [0] = '\0';
  CreateDefLineExEx (NULL, bsp, buf, sizeof (buf) - 1, 0,
                     NULL, NULL, ignoreExisting, FALSE);
  title = NewCreateDefLine (NULL, bsp, ignoreExisting, FALSE);
  if (StringHasNoText (title)) {
    title = StringSave ("?");
  }

  if (StringCmp (buf, title) != 0) {
    if (cfp->ofp != NULL) {
      fprintf (cfp->ofp, "<  %s %s\n", id, buf);
      fprintf (cfp->ofp, ">  %s %s\n", id, title);
    }
    printf ("<  %s %s\n", id, buf);
    printf (">  %s %s\n", id, title);
    fflush (stdout);
  }

  MemFree (title);
}

static void DoFastaExist (
  BioseqPtr bsp,
  Pointer userdata
)

{
  DoFastaComp (bsp, userdata, FALSE);
}

static void DoFastaRegen (
  BioseqPtr bsp,
  Pointer userdata
)

{
  DoFastaComp (bsp, userdata, TRUE);
}

static void DoFastaFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CSpeedFlagPtr  cfp;
  CharPtr        str;

  if (sfp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (cfp->ofp != NULL) {
    str = SeqLocPrint (sfp->location);
    if (str != NULL) {
      fprintf (cfp->ofp, "> [%s]\n", str);
      MemFree (str);
    }
    SeqLocFastaStream (sfp->location, cfp->ofp, STREAM_EXPAND_GAPS, 70, 0, 0);
  } else {
    SeqPortStreamLoc (sfp->location, STREAM_EXPAND_GAPS, NULL, EmptyStreamProc);
  }
}

static void DoFastaTrans (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ByteStorePtr   bs;
  CSpeedFlagPtr  cfp;
  CharPtr        seq, str;

  if (sfp == NULL) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (sfp->data.choice != SEQFEAT_CDREGION) return;
  bs = ProteinFromCdRegion (sfp, FALSE);
  if (bs == NULL) return;

  seq = (CharPtr) BSMerge (bs, NULL);
  BSFree (bs);
  if (seq == NULL) return;

  if (cfp->ofp != NULL) {
    str = SeqLocPrint (sfp->location);
    if (str != NULL) {
      fprintf (cfp->ofp, "> (%s)\n", str);
      MemFree (str);
    }
    fprintf (cfp->ofp, "%s\n", seq);
  }

  MemFree (seq);
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

static void MarkTitles (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

static void DoProcess (
  SeqEntryPtr sep,
  Uint2 entityID,
  CSpeedFlagPtr cfp
)

{
  Char            id [64];
  ErrSev          oldErrSev;
  ValidStructPtr  vsp;

  if (sep == NULL || cfp == NULL) return;

  if (StringChr (cfp->clean, 't') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->clean, 'a') != NULL) {
    AssignIDsInEntity (entityID, 0, NULL);
  }
  if (StringChr (cfp->clean, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->clean, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
  }

  if (StringChr (cfp->index, 'f') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
  }

  if (StringChr (cfp->seq, 'c') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaExist);
  }
  if (StringChr (cfp->seq, 'C') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaRegen);
  }
  if (StringChr (cfp->seq, 's') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaSeq);
  }
  if (StringChr (cfp->seq, 'S') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaSeq);
  }
  if (StringChr (cfp->seq, 'r') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaRaw);
  }
  if (StringChr (cfp->seq, 'd') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaDefline);
  }
  if (StringChr (cfp->seq, 'D') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaDefline);
  }
  if (StringChr (cfp->seq, 'T') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, (Pointer) cfp, DoFastaDefline);
  }
  if (StringChr (cfp->seq, 'x') != NULL) {
    VisitBioseqsInSep (sep, (Pointer) cfp, DoNewFastaDefline);
  }
  if (StringChr (cfp->seq, 'X') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, (Pointer) cfp, DoNewFastaDefline);
  }
  
  if (StringChr (cfp->seq, 'f') != NULL) {
    VisitFeaturesInSep (sep, (Pointer) cfp, DoFastaFeat);
  }
  if (StringChr (cfp->seq, 't') != NULL) {
    VisitFeaturesInSep (sep, (Pointer) cfp, DoFastaTrans);
  }

  if (StringChr (cfp->feat, 'v') != NULL) {
    VisitFeaturesInSep (sep, NULL, DoVisitFeaturesTest);
  }
  if (StringChr (cfp->feat, 'g') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    VisitFeaturesInSep (sep, (Pointer) cfp, DoGeneOverlapPrintTest);
  }
  if (StringChr (cfp->feat, 'h') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    VisitFeaturesInSep (sep, (Pointer) cfp, DoGeneOverlapSpeedTest);
  }
  if (StringChr (cfp->feat, 'x') != NULL) {
  }
  if (StringChr (cfp->feat, 'o') != NULL) {
  }
  if (StringChr (cfp->feat, 'd') != NULL) {
  }
  if (StringChr (cfp->feat, 't') != NULL) {
    SeqEntryToGnbk (sep, NULL, FTABLE_FMT, SEQUIN_MODE, NORMAL_STYLE,
                    0, 0, SHOW_PROT_FTABLE, NULL, cfp->ofp);
  }
  if (StringChr (cfp->feat, 's') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    cfp->nucbsp = FindNucBioseq (sep);
    if (cfp->nucbsp != NULL) {
      BioseqToGeneticCode (cfp->nucbsp, &(cfp->genCode), NULL, NULL, NULL, 0, NULL);
      SeqIdWrite (cfp->nucbsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
      fprintf (cfp->ofp, "%s\n", id);
      VisitBioseqsInSep (sep, (Pointer) cfp, DoSuggestIntervals);
      cfp->nucbsp = NULL;
      cfp->genCode = 0;
    }
  }
  if (StringChr (cfp->feat, 'S') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    cfp->nucbsp = FindNucBioseq (sep);
    if (cfp->nucbsp != NULL) {
      BioseqToGeneticCode (cfp->nucbsp, &(cfp->genCode), NULL, NULL, NULL, 0, NULL);
      SetBatchSuggestNucleotide (cfp->nucbsp, cfp->genCode);
      SeqIdWrite (cfp->nucbsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
      fprintf (cfp->ofp, "%s\n", id);
      VisitBioseqsInSep (sep, (Pointer) cfp, DoSuggestIntervals);
      ClearBatchSuggestNucleotide ();
      cfp->nucbsp = NULL;
      cfp->genCode = 0;
    }
  }
  if (StringChr (cfp->feat, 'c') != NULL) {
    VisitFeaturesInSep (sep, (Pointer) cfp, DoVisitCodingRegions);
  }

  if (StringChr (cfp->desc, 'b') != NULL) {
  }
  if (StringChr (cfp->desc, 't') != NULL) {
  }

  if (StringChr (cfp->verify, 'v') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    vsp = ValidStructNew ();
    if (vsp != NULL) {
      vsp->useSeqMgrIndexes = TRUE;
      vsp->suppressContext = TRUE;
      vsp->seqSubmitParent = TRUE;
      vsp->testLatLonSubregion = TRUE;
      oldErrSev = ErrSetMessageLevel (SEV_NONE);
      vsp->errfunc = ValidCallback;
      vsp->userdata = (Pointer) cfp->ofp;
      /* vsp->convertGiToAccn = FALSE; */
      ValidateSeqEntry (sep, vsp);
      ValidStructFree (vsp);
      ErrSetMessageLevel (oldErrSev);
    }
  }
  if (StringChr (cfp->verify, 'b') != NULL) {
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE,
                    0, 0, 0, NULL, cfp->ofp);
  }

  if (cfp->ofp != NULL) {
    fflush (cfp->ofp);
  }
}

static void ProcessSingleRecord (
  CharPtr filename,
  CSpeedFlagPtr cfp
)

{
  AsnIoPtr      aip;
  BioseqPtr     bsp;
  ValNodePtr    bsplist = NULL;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FileCache     fc;
  FILE          *fp;
  Int1          iotype;
  Char          line [512];
  Int4          maxio = 1;
  SeqEntryPtr   sep;
  time_t        starttime, stoptime, worsttime;
  CharPtr       str;
  Int4          x;

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  if (StringChr (cfp->io, 'r') != NULL) {
    maxio = cfp->maxcount;
  }

  starttime = GetSecs ();

  for (x = 0; x < maxio; x++) {
    if (entityID != 0) {
      ObjMgrFreeByEntityID (entityID);
      entityID = 0;
      dataptr = NULL;
    }

    if (cfp->type == 1) {

      fp = FileOpen (filename, "r");
      if (fp == NULL) {
        Message (MSG_POSTERR, "Failed to open '%s'", filename);
        return;
      }

      dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

      FileClose (fp);

      entityID = ObjMgrRegister (datatype, dataptr);

    } else if (cfp->type >= 2 && cfp->type <= 5) {

      aip = AsnIoOpen (filename, cfp->binary? "rb" : "r");
      if (aip == NULL) {
        Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
        return;
      }

      switch (cfp->type) {
        case 2 :
          dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
          datatype = OBJ_SEQENTRY;
          break;
        case 3 :
          dataptr = (Pointer) BioseqAsnRead (aip, NULL);
          datatype = OBJ_BIOSEQ;
          break;
        case 4 :
          dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
          datatype = OBJ_BIOSEQSET;
          break;
        case 5 :
          dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
          datatype = OBJ_SEQSUB;
          break;
        default :
          break;
      }

      AsnIoClose (aip);

      entityID = ObjMgrRegister (datatype, dataptr);

    } else if (cfp->type == 6) {

      fp = FileOpen (filename, "r");
      if (fp == NULL) {
        Message (MSG_POSTERR, "Failed to open '%s'", filename);
        return;
      }

      dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

      FileClose (fp);

      entityID = ObjMgrRegister (datatype, dataptr);

    } else if (cfp->type == 7) {

      fp = FileOpen (filename, "r");
      if (fp == NULL) {
        Message (MSG_POSTERR, "Failed to open '%s'", filename);
        return;
      }

      FileCacheSetup (&fc, fp);

      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);

      return;

    } else {
      Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) cfp->type);
      return;
    }
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    stoptime = GetSecs ();
    worsttime = stoptime - starttime;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ASN reading time %ld seconds\n", (long) worsttime);
      fflush (cfp->logfp);
    }

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {

      if (cfp->lock) {
        starttime = GetSecs ();

        bsplist = LockFarComponents (sep);

        stoptime = GetSecs ();
        worsttime = stoptime - starttime;
        if (cfp->logfp != NULL) {
          fprintf (cfp->logfp, "Far component locking time %ld seconds\n", (long) worsttime);
          fflush (cfp->logfp);
        }
      }

      if (StringChr (cfp->io, 'w') != NULL) {
        starttime = GetSecs ();

        iotype = ASNIO_TEXT_OUT;
        if (StringChr (cfp->io, 'b') != NULL) {
          iotype = ASNIO_BIN_OUT;
        }

        for (x = 0; x < cfp->maxcount; x++) {
          aip = AsnIoNew (iotype, cfp->ofp, NULL, NULL, NULL);
          if (aip != NULL) {
            SeqEntryAsnWrite (sep, aip, NULL);
            AsnIoFree (aip, FALSE);
          }
        }

        stoptime = GetSecs ();
        worsttime = stoptime - starttime;
        if (cfp->logfp != NULL) {
          fprintf (cfp->logfp, "ASN writing time %ld seconds\n", (long) worsttime);
          fflush (cfp->logfp);
        }
      }

      starttime = GetSecs ();

      for (x = 0; x < cfp->maxcount; x++) {
        DoProcess (sep, entityID, cfp);
      }

      stoptime = GetSecs ();
      worsttime = stoptime - starttime;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "Internal processing time %ld seconds\n", (long) worsttime);
        fflush (cfp->logfp);
      }

      ObjMgrFreeByEntityID (entityID);

      bsplist = UnlockFarComponents (bsplist);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

static void ProcessMultipleRecord (
  CharPtr filename,
  CSpeedFlagPtr cfp
)

{
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  BioseqPtr    bsp;
  Char         buf [41];
  Uint2        entityID;
  FILE         *fp;
  SeqEntryPtr  fsep;
  Char         longest [41];
  Int4         numrecords, x;
  SeqEntryPtr  sep;
  time_t       starttime, stoptime, worsttime;
#ifdef OS_UNIX
  Char         cmmd [256];
  CharPtr      gzcatprog;
  int          ret;
  Boolean      usedPopen = FALSE;
#endif

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

#ifndef OS_UNIX
  if (cfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

#ifdef OS_UNIX
  if (cfp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* cfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, cfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, cfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (cfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n\n", filename);
    fflush (cfp->logfp);
  }

  longest [0] = '\0';
  worsttime = 0;
  numrecords = 0;

  atp = cfp->atp_bss;

  while ((atp = AsnReadId (aip, cfp->amp, atp)) != NULL) {
    if (atp == cfp->atp_se) {

      sep = SeqEntryAsnRead (aip, atp);
      if (sep != NULL) {

        entityID = ObjMgrGetEntityIDForChoice (sep);

        fsep = FindNthBioseq (sep, 1);
        if (fsep != NULL && fsep->choice == 1) {
          bsp = (BioseqPtr) fsep->data.ptrvalue;
          if (bsp != NULL) {
            SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
            if (cfp->logfp != NULL) {
              fprintf (cfp->logfp, "%s\n", buf);
              fflush (cfp->logfp);
            }
          }
        }

        starttime = GetSecs ();

        for (x = 0; x < cfp->maxcount; x++) {
          DoProcess (sep, entityID, cfp);
        }
        stoptime = GetSecs ();

        if (stoptime - starttime > worsttime) {
          worsttime = stoptime - starttime;
          StringCpy (longest, buf);
        }
        numrecords++;

        ObjMgrFreeByEntityID (entityID);
      }

    } else {

      AsnReadVal (aip, atp, NULL);
    }
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  if (cfp->logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (cfp->logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (cfp->logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (cfp->logfp);
  }
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  CSpeedFlagPtr  cfp;

  if (StringHasNoText (filename)) return;
  cfp = (CSpeedFlagPtr) userdata;
  if (cfp == NULL) return;

  if (cfp->batch) {
    ProcessMultipleRecord (filename, cfp);
  } else {
    ProcessSingleRecord (filename, cfp);
  }
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define i_argInputFile     1
#define o_argOutputFile    2
#define f_argFilter        3
#define x_argSuffix        4
#define a_argType          5
#define b_argBinary        6
#define c_argCompressed    7
#define l_argLockFar       8
#define L_argLogFile       9
#define R_argRemote       10
#define X_argMaxCount     11
#define O_argInOut        12
#define K_argClean        13
#define P_argSkip         14
#define I_argIndex        15
#define S_argSeq          16
#define F_argFeat         17
#define D_argDesc         18
#define V_argVerify       19

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Processing\n"
   "      f FASTA\n"
   "      l Read by Lines", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Max Repeat Count", "1", NULL, NULL,
    TRUE, 'X', ARG_INT, 0.0, 0, NULL},
  {"Input Output\n"
   "      r Read ASN.1\n"
   "      w Write Text ASN.1\n"
   "      wb Write Binary ASN.1", NULL, NULL, NULL,
    TRUE, 'O', ARG_STRING, 0.0, 0, NULL},
  {"Cleanup\n"
   "      t Remove Titles\n"
   "      a AssignIDsInEntity\n"
   "      b BasicSeqEntryCleanup\n"
   "      s SeriousSeqEntryCleanup", NULL, NULL, NULL,
    TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
  {"Skip\n"
   "      s Segmented Set Components\n"
   "      v Virtual Bioseqs", NULL, NULL, NULL,
    TRUE, 'P', ARG_STRING, 0.0, 0, NULL},
  {"Index\n"
   "      f Feature Indexing", NULL, NULL, NULL,
    TRUE, 'I', ARG_STRING, 0.0, 0, NULL},
  {"Sequence\n"
   "      c Compare FASTA Deflines\n"
   "      C Compare Regenerated FASTA Deflines\n"
   "      s FASTA of Sequence\n"
   "      S Indexed FASTA\n"
   "      r Raw FASTA no Defline\n"
   "      d Just FASTA Defline\n"
   "      D Indexed FASTA Defline\n"
   "      T Regenerate FASTA Titles\n"
   "      x New FASTA Titles\n"
   "      X Regenerate new FASTA Titles\n"
   "      f FASTA by Feature\n"
   "      t FASTA of Translation", NULL, NULL, NULL,
    TRUE, 'S', ARG_STRING, 0.0, 0, NULL},
  {"Feature\n"
   "      v Visit Features\n"
   "      g Gene Overlap Print\n"
   "      h Gene Overlap Speed\n"
   "      x Gene by Xref\n"
   "      o Operon by Overlap\n"
   "      d Feature by ID\n"
   "      t Feature Table\n"
   "      s Slow Suggest Intervals\n"
   "      S Indexed Suggest Intervals\n"
   "      c Coding Region Intervals", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Descriptor\n"
   "      b BioSource\n"
   "      t Title", NULL, NULL, NULL,
    TRUE, 'D', ARG_STRING, 0.0, 0, NULL},
  {"Verification\n"
   "      v Validate with Normal Stringency\n"
   "      b Generate GenBank Flatfile\n", NULL, NULL, NULL,
    TRUE, 'V', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char            app [64], type;
  CSpeedFlagData  cfd;
  CharPtr         directory, filter, infile, logfile, outfile, str, suffix;
  Boolean         remote;
  time_t          runtime, starttime, stoptime;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
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

  sprintf (app, "cspeedtest %s", CSPEEDTEST_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &cfd, 0, sizeof (CSpeedFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  cfd.batch = FALSE;
  cfd.binary = (Boolean) myargs [b_argBinary].intvalue;
  cfd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  cfd.lock = (Boolean) myargs [l_argLockFar].intvalue;
  cfd.type = 1;

  str = myargs [a_argType].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    type = str [0];
  } else {
    type = 'a';
  }

  type = TO_LOWER (type);
  switch (type) {
    case 'a' :
      cfd.type = 1;
      break;
    case 'e' :
      cfd.type = 2;
      break;
    case 'b' :
      cfd.type = 3;
      break;
    case 's' :
      cfd.type = 4;
      break;
    case 'm' :
      cfd.type = 5;
      break;
    case 't' :
      cfd.type = 1;
      cfd.batch = TRUE;
      break;
    case 'f' :
      cfd.type = 6;
      break;
    case 'l' :
      cfd.type = 7;
      break;
    default :
      cfd.type = 1;
      break;
  }

  remote = (Boolean) myargs [R_argRemote].intvalue;

  cfd.maxcount = myargs [X_argMaxCount].intvalue;
  if (cfd.maxcount < 1) {
    cfd.maxcount = 1;
  }

  cfd.io = myargs [O_argInOut].strvalue;
  cfd.clean = myargs [K_argClean].strvalue;
  cfd.skip = myargs [P_argSkip].strvalue;
  cfd.index = myargs [I_argIndex].strvalue;
  cfd.seq = myargs [S_argSeq].strvalue;
  cfd.feat = myargs [F_argFeat].strvalue;
  cfd.desc = myargs [D_argDesc].strvalue;
  cfd.verify = myargs [V_argVerify].strvalue;

  cfd.amp = AsnAllModPtr ();
  cfd.atp_bss = AsnFind ("Bioseq-set");
  cfd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  cfd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  cfd.atp_bsc = AsnFind ("Bioseq-set.class");
  cfd.bssp_atp = AsnLinkType (NULL, cfd.atp_bss);

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;
  if (StringDoesHaveText (logfile)) {
    cfd.logfp = FileOpen (logfile, "w");
  }

  if (remote) {
    PubSeqFetchEnable ();
  }

  if (StringDoesHaveText (outfile)) {
    cfd.ofp = FileOpen (outfile, "w");
  }

  starttime = GetSecs ();

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, FALSE, ProcessOneRecord, (Pointer) &cfd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, (Pointer) &cfd);
  }

  if (cfd.ofp != NULL) {
    FileClose (cfd.ofp);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;
  if (cfd.logfp != NULL) {
    fprintf (cfd.logfp, "Finished in %ld seconds\n", (long) runtime);
    FileClose (cfd.logfp);
  }
  printf ("Finished in %ld seconds\n", (long) runtime);

  if (remote) {
    PubSeqFetchDisable ();
  }

  return 0;
}

