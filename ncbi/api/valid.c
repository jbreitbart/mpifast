/*  valid.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  valid.c
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/94
*
* $Revision: 6.1178 $
*
* File Description:  Sequence editing utilities
*
* Modifications:  
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
*
* ==========================================================================
*/

static char    *this_module = "valid";

#define THIS_MODULE this_module

static char    *this_file = __FILE__;

#define THIS_FILE this_file

#include <ncbi.h>
#include <objfdef.h>
#include <valid.h>
#include <validerr.h>
#include <sqnutils.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <objsub.h>
#include <asn2gnbi.h>
#include <explore.h>
#include <subutil.h>
#include <tofasta.h>
#include <findrepl.h>

/*****************************************************************************
*
*   NOTE: look at all the ValidErr calls with severity=0. Some should be
*   bumped up later. Look also for string "PARSER"
*
*****************************************************************************/



#ifdef VAR_ARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif

static ValidStructPtr globalvsp;        /* for spell checker */

NLM_EXTERN void CDECL ValidErr VPROTO ((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));
static void     ValidateBioseqInst (GatherContextPtr gcp);
static void     ValidateBioseqContext (GatherContextPtr gcp);
static void     ValidateBioseqSet (GatherContextPtr gcp);
static void     ValidateGraphsOnBioseq (GatherContextPtr gcp);
static void     ValidateBioseqHist (GatherContextPtr gcp);
static void     SpellCheckSeqDescr (GatherContextPtr gcp);
NLM_EXTERN void CdTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
NLM_EXTERN void MrnaTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
NLM_EXTERN void ValidateSeqFeat (GatherContextPtr gcp);
NLM_EXTERN void ValidateSeqLoc (ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix);
NLM_EXTERN Boolean PatchBadSequence (BioseqPtr bsp);
NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf);
NLM_EXTERN void SpellCheckSeqFeat (GatherContextPtr gcp);
NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str);
NLM_EXTERN void SpliceCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     CdConflictCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     SpliceCheckEx (ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll);
static void     CdsProductIdCheck (ValidStructPtr vsp, SeqFeatPtr sfp);
static void     ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, SeqFeatPtr sfp, ValNodePtr sdp);
static void     ValidatePubdesc (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp);
static void     LookForMultiplePubs (ValidStructPtr vsp, GatherContextPtr gcp, SeqDescrPtr sdp);
static void     ValidateSfpCit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp);
static void     ValidateAffil (ValidStructPtr vsp, AffilPtr ap);
static TextFsaPtr GetSpecificECNumberFSA (void);
static TextFsaPtr GetAmbiguousECNumberFSA (void);
static TextFsaPtr GetDeletedECNumberFSA (void);
static TextFsaPtr GetReplacedECNumberFSA (void);
static Boolean ECnumberNotInList (CharPtr str);
static Boolean ECnumberWasDeleted (CharPtr str);
static Boolean ECnumberWasReplaced (CharPtr str);

/* alignment validator */
NLM_EXTERN Boolean ValidateSeqAlignWithinValidator (ValidStructPtr vsp, SeqEntryPtr sep, Boolean find_remote_bsp, Boolean do_hist_assembly);

static ValNodePtr genetic_code_name_list = NULL;

/*****************************************************************************
*
*   Perform Validation Checks on a SeqEntry
*
*****************************************************************************/

NLM_EXTERN void ValidStructClear (ValidStructPtr vsp)
{                               /* 0 out a ValidStruct */
  CharPtr         errbuf;
  Int2            cutoff;
  Boolean         patch_seq;
  SpellCheckFunc  spellfunc;
  SpellCallBackFunc spellcallback;
  Boolean         onlyspell;
  Boolean         justwarnonspell;
  Boolean         useSeqMgrIndexes;
  Boolean         suppressContext;
  Boolean         validateAlignments;
  Boolean         farIDsInAlignments;
  Boolean         alignFindRemoteBsp;
  Boolean         doSeqHistAssembly;
  Boolean         alwaysRequireIsoJTA;
  Boolean         farFetchCDSproducts;
  Boolean         farFetchMRNAproducts;
  Boolean         locusTagGeneralMatch;
  Boolean         validateIDSet;
  Boolean         seqSubmitParent;
  Boolean         justShowAccession;
  Boolean         ignoreExceptions;
  Boolean         validateExons;
  Boolean         inferenceAccnCheck;
  Boolean         testLatLonSubregion;
  Boolean         strictLatLonCountry;
  Boolean         indexerVersion;
  Int2            validationLimit;
  ValidErrorFunc  errfunc;
  Pointer         userdata;
  Boolean         convertGiToAccn;
  TextFsaPtr      sourceQualTags;
  TextFsaPtr      modifiedBases;
  Boolean         is_htg_in_sep;
  Boolean         is_barcode_sep;
  Boolean         is_refseq_in_sep;
  Boolean         is_gps_in_sep;
  Boolean         is_embl_ddbj_in_sep;
  Boolean         is_insd_in_sep;
  Boolean         only_lcl_gnl_in_sep;
  Boolean         has_gnl_prot_sep;
  Boolean         is_smupd_in_sep;
  Boolean         feat_loc_has_gi;
  Boolean         feat_prod_has_gi;
  Boolean         far_fetch_failure;

  if (vsp == NULL)
    return;

  errbuf = vsp->errbuf;
  cutoff = vsp->cutoff;
  patch_seq = vsp->patch_seq;
  spellfunc = vsp->spellfunc;
  spellcallback = vsp->spellcallback;
  onlyspell = vsp->onlyspell;
  justwarnonspell = vsp->justwarnonspell;
  useSeqMgrIndexes = vsp->useSeqMgrIndexes;
  suppressContext = vsp->suppressContext;
  validateAlignments = vsp->validateAlignments;
  farIDsInAlignments = vsp->farIDsInAlignments;
  alignFindRemoteBsp = vsp->alignFindRemoteBsp;
  doSeqHistAssembly = vsp->doSeqHistAssembly;
  alwaysRequireIsoJTA = vsp->alwaysRequireIsoJTA;
  farFetchCDSproducts = vsp->farFetchCDSproducts;
  farFetchMRNAproducts = vsp->farFetchMRNAproducts;
  locusTagGeneralMatch = vsp->locusTagGeneralMatch;
  validateIDSet = vsp->validateIDSet;
  seqSubmitParent = vsp->seqSubmitParent;
  justShowAccession = vsp->justShowAccession;
  ignoreExceptions = vsp->ignoreExceptions;
  validateExons = vsp->validateExons;
  inferenceAccnCheck = vsp->inferenceAccnCheck;
  testLatLonSubregion = vsp->testLatLonSubregion;
  strictLatLonCountry = vsp->strictLatLonCountry;
  indexerVersion = vsp->indexerVersion;
  validationLimit = vsp->validationLimit;
  errfunc = vsp->errfunc;
  userdata = vsp->userdata;
  convertGiToAccn = vsp->convertGiToAccn;
  sourceQualTags = vsp->sourceQualTags;
  modifiedBases = vsp->modifiedBases;
  is_htg_in_sep = vsp->is_htg_in_sep;
  is_barcode_sep = vsp->is_barcode_sep;
  is_refseq_in_sep = vsp->is_refseq_in_sep;
  is_gps_in_sep = vsp->is_gps_in_sep;
  is_embl_ddbj_in_sep = vsp->is_embl_ddbj_in_sep;
  is_insd_in_sep = vsp->is_insd_in_sep;
  only_lcl_gnl_in_sep = vsp->only_lcl_gnl_in_sep;
  has_gnl_prot_sep = vsp->has_gnl_prot_sep;
  is_smupd_in_sep = vsp->is_smupd_in_sep;
  feat_loc_has_gi = vsp->feat_loc_has_gi;
  feat_prod_has_gi = vsp->feat_prod_has_gi;
  far_fetch_failure = vsp->far_fetch_failure;
  MemSet ((VoidPtr) vsp, 0, sizeof (ValidStruct));
  vsp->errbuf = errbuf;
  vsp->cutoff = cutoff;
  vsp->patch_seq = patch_seq;
  vsp->spellfunc = spellfunc;
  vsp->spellcallback = spellcallback;
  vsp->onlyspell = onlyspell;
  vsp->justwarnonspell = justwarnonspell;
  vsp->useSeqMgrIndexes = useSeqMgrIndexes;
  vsp->suppressContext = suppressContext;
  vsp->validateAlignments = validateAlignments;
  vsp->farIDsInAlignments = farIDsInAlignments;
  vsp->alignFindRemoteBsp = alignFindRemoteBsp;
  vsp->doSeqHistAssembly = doSeqHistAssembly;
  vsp->alwaysRequireIsoJTA = alwaysRequireIsoJTA;
  vsp->farFetchCDSproducts = farFetchCDSproducts;
  vsp->farFetchMRNAproducts = farFetchMRNAproducts;
  vsp->locusTagGeneralMatch = locusTagGeneralMatch;
  vsp->validateIDSet = validateIDSet;
  vsp->seqSubmitParent = seqSubmitParent;
  vsp->justShowAccession = justShowAccession;
  vsp->ignoreExceptions = ignoreExceptions;
  vsp->validateExons = validateExons;
  vsp->inferenceAccnCheck = inferenceAccnCheck;
  vsp->testLatLonSubregion = testLatLonSubregion;
  vsp->strictLatLonCountry = strictLatLonCountry;
  vsp->indexerVersion = indexerVersion;
  vsp->validationLimit = validationLimit;
  vsp->errfunc = errfunc;
  vsp->userdata = userdata;
  vsp->convertGiToAccn = convertGiToAccn;
  vsp->sourceQualTags = sourceQualTags;
  vsp->modifiedBases = modifiedBases;
  vsp->is_htg_in_sep = is_htg_in_sep;
  vsp->is_barcode_sep = is_barcode_sep;
  vsp->is_refseq_in_sep = is_refseq_in_sep;
  vsp->is_gps_in_sep = is_gps_in_sep;
  vsp->is_embl_ddbj_in_sep = is_embl_ddbj_in_sep;
  vsp->is_insd_in_sep = is_insd_in_sep;
  vsp->only_lcl_gnl_in_sep = only_lcl_gnl_in_sep;
  vsp->has_gnl_prot_sep = has_gnl_prot_sep;
  vsp->is_smupd_in_sep = is_smupd_in_sep;
  vsp->feat_loc_has_gi = feat_loc_has_gi;
  vsp->feat_prod_has_gi = feat_prod_has_gi;
  vsp->far_fetch_failure = far_fetch_failure;
  return;
}

NLM_EXTERN ValidStructPtr ValidStructNew (void)
{
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) MemNew (sizeof (ValidStruct));
  return vsp;
}

NLM_EXTERN ValidStructPtr ValidStructFree (ValidStructPtr vsp)
{
  if (vsp == NULL)
    return vsp;

  MemFree (vsp->errbuf);
  TextFsaFree (vsp->sourceQualTags);
  TextFsaFree (vsp->modifiedBases);
  return (ValidStructPtr) MemFree (vsp);
}

/*****************************************************************************
*
*   ValidErr()
*
*****************************************************************************/

static void ChangeSeqIdToBestID (SeqIdPtr sip)
{
  BioseqPtr       bsp;
  SeqIdPtr        id;
  Pointer         pnt;

  if (sip == NULL)
    return;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL)
    return;
  id = SeqIdDup (SeqIdFindWorst (bsp->id));
  if (id == NULL)
    return;
  /* now remove SeqId contents to reuse SeqId valnode */
  pnt = sip->data.ptrvalue;
  switch (sip->choice) {
  case SEQID_LOCAL:            /* local */
    ObjectIdFree ((ObjectIdPtr) pnt);
    break;
  case SEQID_GIBBSQ:           /* gibbseq */
  case SEQID_GIBBMT:           /* gibbmt */
    break;
  case SEQID_GIIM:             /* giimid */
    GiimFree ((GiimPtr) pnt);
    break;
  case SEQID_GENBANK:          /* genbank */
  case SEQID_EMBL:             /* embl */
  case SEQID_PIR:              /* pir   */
  case SEQID_SWISSPROT:        /* swissprot */
  case SEQID_OTHER:            /* other */
  case SEQID_DDBJ:
  case SEQID_PRF:
  case SEQID_TPG:
  case SEQID_TPE:
  case SEQID_TPD:
  case SEQID_GPIPE:
    TextSeqIdFree ((TextSeqIdPtr) pnt);
    break;
  case SEQID_PATENT:           /* patent seq id */
    PatentSeqIdFree ((PatentSeqIdPtr) pnt);
    break;
  case SEQID_GENERAL:          /* general */
    DbtagFree ((DbtagPtr) pnt);
    break;
  case SEQID_GI:               /* gi */
    break;
  case SEQID_PDB:
    PDBSeqIdFree ((PDBSeqIdPtr) pnt);
    break;
  }
  sip->choice = id->choice;
  sip->data.ptrvalue = id->data.ptrvalue;
  SeqIdStripLocus (sip);
}

static void ChangeSeqLocToBestID (SeqLocPtr slp)
{
  SeqLocPtr       loc;
  PackSeqPntPtr   psp;
  SeqBondPtr      sbp;
  SeqIntPtr       sinp;
  SeqIdPtr        sip;
  SeqPntPtr       spp;

  while (slp != NULL) {
    switch (slp->choice) {
    case SEQLOC_NULL:
      break;
    case SEQLOC_EMPTY:
    case SEQLOC_WHOLE:
      sip = (SeqIdPtr) slp->data.ptrvalue;
      ChangeSeqIdToBestID (sip);
      break;
    case SEQLOC_INT:
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sip = sinp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        sip = spp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PACKED_PNT:
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        sip = psp->id;
        ChangeSeqIdToBestID (sip);
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:
    case SEQLOC_EQUIV:
      loc = (SeqLocPtr) slp->data.ptrvalue;
      while (loc != NULL) {
        ChangeSeqLocToBestID (loc);
        loc = loc->next;
      }
      break;
    case SEQLOC_BOND:
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        spp = (SeqPntPtr) sbp->a;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToBestID (sip);
        }
        spp = (SeqPntPtr) sbp->b;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToBestID (sip);
        }
      }
      break;
    case SEQLOC_FEAT:
      break;
    default:
      break;
    }
    slp = slp->next;
  }
}

static Int2 WorstBioseqLabel (BioseqPtr bsp, CharPtr buffer, Int2 buflen, Uint1 content)
{
  CharPtr         tmp;
  Char            label[60];
  Int2            diff, len;
  SeqIdPtr        sip;
  AsnModulePtr    amp;
  AsnTypePtr      ratp, matp;

  if ((bsp == NULL) || (buflen < 1))
    return 0;

  len = buflen;
  label[0] = '\0';

  if (content != OM_LABEL_TYPE) {
    sip = SeqIdStripLocus (SeqIdDup (SeqIdFindWorst (bsp->id)));
    SeqIdWrite (sip, label, PRINTID_FASTA_SHORT, 39);
    SeqIdFree (sip);
    if (content == OM_LABEL_CONTENT)
      return LabelCopy (buffer, label, buflen);

    diff = LabelCopyExtra (buffer, label, buflen, NULL, ": ");
    buflen -= diff;
    buffer += diff;
  }

  amp = AsnAllModPtr ();
  ratp = AsnTypeFind (amp, "Seq-inst.repr");
  matp = AsnTypeFind (amp, "Seq-inst.mol");

  label[0] = '\0';
  tmp = label;
  tmp = StringMove (tmp, AsnEnumTypeStr (ratp, (Int2) (bsp->repr)));
  tmp = StringMove (tmp, ", ");
  tmp = StringMove (tmp, AsnEnumTypeStr (matp, (Int2) (bsp->mol)));
  sprintf (tmp, " len= %ld", (long) (bsp->length));
  diff = LabelCopy (buffer, label, buflen);
  buflen -= diff;
  buffer += diff;

  if (content != OM_LABEL_SUMMARY)
    return (len - buflen);

  return (len - buflen);        /* SUMMARY not done yet */
}

static CharPtr categoryLabel [] = {
  NULL, "SEQ_INST", "SEQ_DESCR", "GENERIC", "SEQ_PKG", "SEQ_FEAT", "SEQ_ALIGN", "SEQ_GRAPH", "SEQ_ANNOT"
};

NLM_EXTERN CharPtr GetValidCategoryName (int errcode)

{
  if (errcode >= 1 && errcode < sizeof (categoryLabel)) return categoryLabel [errcode];
  return NULL;
}

static CharPtr err1Label [] = {
  NULL,
  "ExtNotAllowed",
  "ExtBadOrMissing",
  "SeqDataNotFound",
  "SeqDataNotAllowed",
  "ReprInvalid",
  "CircularProtein",
  "DSProtein",
  "MolNotSet",
  "MolOther",
  "FuzzyLen",
  "InvalidLen",
  "InvalidAlphabet",
  "SeqDataLenWrong",
  "SeqPortFail",
  "InvalidResidue",
  "StopInProtein",
  "PartialInconsistent",
  "ShortSeq",
  "NoIdOnBioseq",
  "BadDeltaSeq",
  "LongHtgsSequence",
  "LongLiteralSequence",
  "SequenceExceeds350kbp",
  "ConflictingIdsOnBioseq",
  "MolNuclAcid",
  "ConflictingBiomolTech",
  "SeqIdNameHasSpace",
  "IdOnMultipleBioseqs",
  "DuplicateSegmentReferences",
  "TrailingX",
  "BadSeqIdFormat",
  "PartsOutOfOrder",
  "BadSecondaryAccn",
  "ZeroGiNumber",
  "RnaDnaConflict",
  "HistoryGiCollision",
  "GiWithoutAccession",
  "MultipleAccessions",
  "HistAssemblyMissing",
  "TerminalNs",
  "UnexpectedIdentifierChange",
  "InternalNsInSeqLit",
  "SeqLitGapLength0",
  "TpaAssmeblyProblem",
  "SeqLocLength",
  "MissingGaps",
  "CompleteTitleProblem",
  "CompleteCircleProblem",
  "BadHTGSeq",
  "GapInProtein",
  "BadProteinStart",
  "TerminalGap",
  "OverlappingDeltaRange",
  "LeadingX",
  "InternalNsInSeqRaw",
  "InternalNsAdjacentToGap",
  "CaseDifferenceInSeqID",
  "DeltaComponentIsGi0",
  "FarFetchFailure",
  "InternalGapsInSeqRaw",
  "SelfReferentialSequence",
  "WholeComponent", 
  "TSAHistAssemblyMissing",
  "ProteinsHaveGeneralID",
  "HighNContent",
  "SeqLitDataLength0",
};

static CharPtr err2Label [] = {
  NULL,
  "BioSourceMissing",
  "InvalidForType",
  "FileOpenCollision",
  "Unknown",
  "NoPubFound",
  "NoOrgFound",
  "MultipleBioSources",
  "NoMolInfoFound",
  "BadCountryCode",
  "NoTaxonID",
  "InconsistentBioSources",
  "MissingLineage",
  "SerialInComment",
  "BioSourceNeedsFocus",
  "BadOrganelle",
  "MultipleChromosomes",
  "BadSubSource",
  "BadOrgMod",
  "InconsistentProteinTitle",
  "Inconsistent",
  "ObsoleteSourceLocation",
  "ObsoleteSourceQual",
  "StructuredSourceNote",
  "UnnecessaryBioSourceFocus",
  "RefGeneTrackingWithoutStatus",
  "UnwantedCompleteFlag",
  "CollidingPublications",
  "TransgenicProblem",
  "TaxonomyLookupProblem",
  "MultipleTitles",
  "RefGeneTrackingOnNonRefSeq",
  "BioSourceInconsistency",
  "FastaBracketTitle",
  "MissingText",
  "BadCollectionDate",
  "BadPCRPrimerSequence",
  "BadPunctuation",
  "BadPCRPrimerName",
  "BioSourceOnProtein",
  "BioSourceDbTagConflict",
  "DuplicatePCRPrimerSequence",
  "MultipleNames",
  "MultipleComments",
  "LatLonProblem",
  "LatLonFormat",
  "LatLonRange",
  "LatLonValue",
  "LatLonCountry",
  "LatLonState",
  "BadSpecificHost",
  "RefGeneTrackingIllegalStatus",
  "ReplacedCountryCode",
  "BadInstitutionCode",
  "BadCollectionCode",
  "BadVoucherID",
  "UnstructuredVoucher",
  "ChromosomeLocation",
  "MultipleSourceQualifiers",
  "UnbalancedParentheses",
  "MultipleSourceVouchers",
  "BadCountryCapitalization",
  "WrongVoucherType",
  "UserObjectProblem"
};

static CharPtr err3Label [] = {
  NULL,
  "NonAsciiAsn",
  "Spell",
  "AuthorListHasEtAl",
  "MissingPubInfo",
  "UnnecessaryPubEquiv",
  "BadPageNumbering",
  "MedlineEntryPub",
  "BadDate",
  "StructuredCitGenCit",
  "CollidingSerialNumbers",
  "EmbeddedScript",
  "PublicationInconsistency"
};

static CharPtr err4Label [] = {
  NULL,
  "NoCdRegionPtr",
  "NucProtProblem",
  "SegSetProblem",
  "EmptySet",
  "NucProtNotSegSet",
  "SegSetNotParts",
  "SegSetMixedBioseqs",
  "PartsSetMixedBioseqs",
  "PartsSetHasSets",
  "FeaturePackagingProblem",
  "GenomicProductPackagingProblem",
  "InconsistentMolInfoBiomols",
  "ArchaicFeatureLocation",
  "ArchaicFeatureProduct",
  "GraphPackagingProblem",
  "InternalGenBankSet",
  "ConSetProblem",
  "NoBioseqFound",
  "INSDRefSeqPackaging",
  "GPSnonGPSPackaging",
  "RefSeqPopSet"
};

static CharPtr err5Label [] = {
  NULL,
  "InvalidForType",
  "PartialProblem",
  "InvalidType",
  "Range",
  "MixedStrand",
  "SeqLocOrder",
  "CdTransFail",
  "StartCodon",
  "InternalStop",
  "NoProtein",
  "MisMatchAA",
  "TransLen",
  "NoStop",
  "TranslExcept",
  "NoProtRefFound",
  "NotSpliceConsensus",
  "OrfCdsHasProduct",
  "GeneRefHasNoData",
  "ExceptInconsistent",
  "ProtRefHasNoData",
  "GenCodeMismatch",
  "RNAtype0",
  "UnknownImpFeatKey",
  "UnknownImpFeatQual",
  "WrongQualOnImpFeat",
  "MissingQualOnImpFeat",
  "PseudoCdsHasProduct",
  "IllegalDbXref",
  "FarLocation",
  "DuplicateFeat",
  "UnnecessaryGeneXref",
  "TranslExceptPhase",
  "TrnaCodonWrong",
  "BothStrands",
  "CDSgeneRange",
  "CDSmRNArange",
  "OverlappingPeptideFeat",
  "SerialInComment",
  "MultipleCDSproducts",
  "FocusOnBioSourceFeature",
  "PeptideFeatOutOfFrame",
  "InvalidQualifierValue",
  "MultipleMRNAproducts",
  "mRNAgeneRange",
  "TranscriptLen",
  "TranscriptMismatches",
  "CDSproductPackagingProblem",
  "DuplicateInterval",
  "PolyAsiteNotPoint",
  "ImpFeatBadLoc",
  "LocOnSegmentedBioseq",
  "UnnecessaryCitPubEquiv",
  "ImpCDShasTranslation",
  "ImpCDSnotPseudo",
  "MissingMRNAproduct",
  "AbuttingIntervals",
  "CollidingGeneNames",
  "MultiIntervalGene",
  "FeatContentDup",
  "BadProductSeqId",
  "RnaProductMismatch",
  "MissingCDSproduct",
  "BadTrnaCodon",
  "BadTrnaAA",
  "OnlyGeneXrefs",
  "UTRdoesNotAbutCDS",
  "BadConflictFlag",
  "ConflictFlagSet",
  "LocusTagProblem",
  "CollidingLocusTags",
  "AltStartCodon",
  "PartialsInconsistent",
  "GenesInconsistent",
  "DuplicateTranslExcept",
  "TranslExceptAndRnaEditing",
  "NoNameForProtein",
  "TaxonDbxrefOnFeature",
  "UnindexedFeature",
  "CDSmRNAmismatch",
  "UnnecessaryException",
  "LocusTagProductMismatch",
  "MrnaTransFail",
  "PseudoCdsViaGeneHasProduct",
  "MissingGeneXref",
  "FeatureCitationProblem",
  "NestedSeqLocMix",
  "WrongQualOnFeature",
  "MissingQualOnFeature",
  "CodonQualifierUsed",
  "UnknownFeatureQual",
  "BadCharInAuthorName",
  "PolyATail",
  "ProteinNameEndsInBracket",
  "CDSwithMultipleMRNAs",
  "MultipleEquivBioSources",
  "MultipleEquivPublications",
  "BadFullLengthFeature",
  "RedundantFields",
  "CDSwithNoMRNAOverlap",
  "FeatureProductInconsistency",
  "ImproperBondLocation",
  "GeneXrefWithoutGene",
  "SeqFeatXrefProblem",
  "ProductFetchFailure",
  "SuspiciousGeneXref",
  "MissingTrnaAA",
  "CollidingFeatureIDs",
  "ExceptionProblem",
  "PolyAsignalNotRange",
  "OldLocusTagMismtach",
  "DuplicateGeneOntologyTerm",
  "InvalidInferenceValue",
  "HpotheticalProteinMismatch",
  "FeatureRefersToAccession",
  "SelfReferentialProduct",
  "ITSdoesNotAbutRRNA",
  "FeatureSeqIDCaseDifference",
  "FeatureLocationIsGi0",
  "GapFeatureProblem",
  "PseudoCdsHasProtXref",
  "ErroneousException",
  "SegmentedGeneProblem",
  "WholeLocation",
  "BadEcNumberFormat",
  "BadEcNumberValue",
  "EcNumberProblem",
  "VectorContamination",
  "MinusStrandProtein",
  "BadProteinName",
  "GeneXrefWithoutLocus",
  "UTRdoesNotExtendToEnd",
  "CDShasTooManyXs",
  "SuspiciousFrame",
  "TerminalXDiscrepancy",
  "UnnecessaryTranslExcept",
  "SuspiciousQualifierValue",
  "NotSpliceConsensusDonor",
  "NotSpliceConsensusAcceptor",
  "RareSpliceConsensusDonor",
  "SeqFeatXrefNotReciprocal",
  "SeqFeatXrefFeatureMissing",
  "FeatureInsideGap",
  "FeatureCrossesGap",
  "BadAuthorSuffix",
  "BadAnticodonAA",
  "BadAnticodonCodon",
  "BadAnticodonStrand",
  "UndesiredGeneSynonym",
  "UndesiredProteinName",
  "FeatureBeginsOrEndsInGap",
  "GeneOntologyTermMissingGOID",
  "PseudoRnaHasProduct",
  "PseudoRnaViaGeneHasProduct",
  "BadRRNAcomponentOrder",
  "BadRRNAcomponentOverlap",
  "MissingGeneLocusTag",
  "MultipleProtRefs",
  "BadInternalCharacter",
  "BadTrailingCharacter",
  "BadTrailingHyphen",
  "MultipleGeneOverlap",
  "BadCharInAuthorLastName",
  "PseudoCDSmRNArange",
  "ExtendablePartialProblem",
  "GeneXrefNeeded"
};

static CharPtr err6Label [] = {
  NULL,
  "SeqIdProblem",
  "StrandRev",
  "DensegLenStart",
  "StartLessthanZero",
  "StartMorethanBiolen",
  "EndLessthanZero",
  "EndMorethanBiolen",
  "LenLessthanZero",
  "LenMorethanBiolen",
  "SumLenStart",
  "AlignDimSeqIdNotMatch",
  "SegsDimSeqIdNotMatch",
  "FastaLike",
  "NullSegs",
  "SegmentGap",
  "SegsDimOne",
  "AlignDimOne",
  "Segtype",
  "BlastAligns",
  "PercentIdentity",
  "ShortAln"
};

static CharPtr err7Label [] = {
  NULL,
  "GraphMin",
  "GraphMax",
  "GraphBelow",
  "GraphAbove",
  "GraphByteLen",
  "GraphOutOfOrder",
  "GraphBioseqLen",
  "GraphSeqLitLen",
  "GraphSeqLocLen",
  "GraphStartPhase",
  "GraphStopPhase",
  "GraphDiffNumber",
  "GraphACGTScore",
  "GraphNScore",
  "GraphGapScore",
  "GraphOverlap",
  "GraphBioseqId",
  "GraphACGTScoreMany",
  "GraphNScoreMany"
};

static CharPtr err8Label [] = {
  NULL,
  "AnnotIDs",
  "AnnotLOCs"
};

NLM_EXTERN CharPtr GetValidErrorName (int errcode, int subcode)

{
  if (errcode < 1 || errcode >= sizeof (categoryLabel)) return NULL;
  switch (errcode) {
    case 1 :
      if (subcode >= 1 && subcode < sizeof (err1Label)) return err1Label [subcode];
      break;
    case 2 :
      if (subcode >= 1 && subcode < sizeof (err2Label)) return err2Label [subcode];
      break;
    case 3 :
      if (subcode >= 1 && subcode < sizeof (err3Label)) return err3Label [subcode];
      break;
    case 4 :
      if (subcode >= 1 && subcode < sizeof (err4Label)) return err4Label [subcode];
      break;
    case 5 :
      if (subcode >= 1 && subcode < sizeof (err5Label)) return err5Label [subcode];
      break;
    case 6 :
      if (subcode >= 1 && subcode < sizeof (err6Label)) return err6Label [subcode];
      break;
    case 7 :
      if (subcode >= 1 && subcode < sizeof (err7Label)) return err7Label [subcode];
      break;
    case 8 :
      if (subcode >= 1 && subcode < sizeof (err8Label)) return err8Label [subcode];
      break;
    default :
      break;
  }
  return NULL;
}

NLM_EXTERN CharPtr GetValidExplanation (int errcode, int subcode)

{
  return Nlm_GetErrLongText (THIS_MODULE, errcode, subcode);
}

static void CustValErr (ValidStructPtr vsp, ErrSev severity, int errcode, int subcode)

{
  CharPtr           accession = NULL, context = NULL, label = NULL, location = NULL,
                    message = NULL, objtype = NULL, product = NULL;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Int2              buflen, diff, wrklen;
  CharPtr           ctmp, tmp;
  Uint2             entityID = 0, itemtype = 0;
  ValidErrorFunc    errfunc;
  GatherContextPtr  gcp;
  Char              id [64];
  Uint4             itemID = 0;
  ObjValNodePtr     ovp;
  SeqDescrPtr       sdp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;
  SeqLocPtr         slp;

  if (vsp == NULL) return;
  errfunc = vsp->errfunc;
  if (errfunc == NULL) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    entityID = gcp->entityID;
    itemtype = gcp->thistype;
    itemID = gcp->itemID;
  }

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  sip = NULL;
  if (vsp->sfp != NULL) {
    sfp = vsp->sfp;
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      sip = SeqIdFindWorst (bsp->id);
    }
  } else if (vsp->descr != NULL) {
    sdp = vsp->descr;
    if (sdp != NULL && sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        bsp = (BioseqPtr) ovp->idx.parentptr;
        if (bsp != NULL) {
          sip = SeqIdFindWorst (bsp->id);
        }
      } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) ovp->idx.parentptr;
        if (bssp != NULL) {
          sep = bssp->seqentry;
          if (sep != NULL) {
            sep = FindNthBioseq (sep, 1);
            if (sep != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              if (bsp != NULL) {
                sip = SeqIdFindWorst (bsp->id);
              }
            }
          }
        }
      }
    }
  } else if (vsp->bsp != NULL) {
    bsp = vsp->bsp;
    sip = SeqIdFindWorst (bsp->id);
  } else if (vsp->bssp != NULL) {
    bssp = vsp->bssp;
    sep = bssp->seqentry;
    if (sep != NULL) {
      sep = FindNthBioseq (sep, 1);
      if (sep != NULL) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          sip = SeqIdFindWorst (bsp->id);
        }
      }
    }
  }
  if (sip != NULL) {
    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id) - 1);
    accession = id;
  }

  if (vsp->sfp != NULL) {
    objtype = "FEATURE";
  } else if (vsp->descr != NULL) {
    objtype = "DESCRIPTOR";
  } else if (vsp->bsp != NULL) {
    objtype = "BIOSEQ";
  } else if (vsp->bssp != NULL) {
    objtype = "BIOSEQ-SET";
  }

  message = vsp->errbuf;

  tmp = vsp->errbuf;
  buflen = 4000;
  while (*tmp != '\0') {
    buflen--;
    tmp++;
  }
  tmp++;
  *tmp = '\0';

  wrklen = buflen;
  if (wrklen > 2000) {
     wrklen -= 1000;
  }

  if (vsp->sfp != NULL) {
    label = tmp;
    diff = FeatDefLabel (vsp->sfp, tmp, wrklen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  } else if (vsp->descr != NULL) {
    label = tmp;
    diff = SeqDescLabel (vsp->descr, tmp, wrklen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  } else if (vsp->bsp != NULL) {
    label = tmp;
    if (vsp->convertGiToAccn) {
      diff = WorstBioseqLabel (vsp->bsp, tmp, wrklen, OM_LABEL_CONTENT);
    } else {
      diff = BioseqLabel (vsp->bsp, tmp, wrklen, OM_LABEL_BOTH);
    }
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  } else if (vsp->bssp != NULL) {
    label = tmp;
    diff = BioseqSetLabel (vsp->bssp, tmp, wrklen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
    *tmp = '\0';
    tmp++;
    *tmp = '\0';
  }

  if (vsp->sfp != NULL) {
    sfp = vsp->sfp;
  
    if (sfp->location != NULL) {
      ctmp = NULL;
      slp = NULL;
      /*
      if (vsp->suppressContext) {
        slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (slp);
        ctmp = SeqLocPrint (slp);
        SeqLocFree (slp);
      } else {
        ctmp = SeqLocPrint (sfp->location);
      }
      */
      slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (slp);
      ctmp = SeqLocPrint (slp);
      SeqLocFree (slp);
      if (ctmp != NULL) {
        if (StringLen (ctmp) > 800) {
          StringCpy (ctmp + 797, "...");
        }
        location = tmp;
        diff = LabelCopyExtra (tmp, ctmp, buflen, "[", "]");
        buflen -= diff;
        tmp += diff;
        MemFree (ctmp);
        *tmp = '\0';
        tmp++;
        *tmp = '\0';

        sip = SeqLocId (sfp->location);
        if (sip != NULL) {
          bsp = BioseqFind (sip);
          if (bsp != NULL) {
            context = tmp;
            diff = LabelCopy (tmp, "[", buflen);
            buflen -= diff;
            tmp += diff;

            diff = BioseqLabel (bsp, tmp, buflen, OM_LABEL_BOTH);
            buflen -= diff;
            tmp += diff;

            diff = LabelCopy (tmp, "]", buflen);
            buflen -= diff;
            tmp += diff;
          }
        }
        *tmp = '\0';
        tmp++;
        *tmp = '\0';
      }
    }
  
    if (sfp->product != NULL) {
      ctmp = NULL;
      slp = NULL;
      /*
      if (vsp->suppressContext) {
        slp = AsnIoMemCopy (sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (slp);
        ctmp = SeqLocPrint (slp);
        SeqLocFree (slp);
      } else {
        ctmp = SeqLocPrint (sfp->product);
      }
      */
      slp = AsnIoMemCopy (sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (slp);
      ctmp = SeqLocPrint (slp);
      SeqLocFree (slp);
      if (ctmp != NULL) {
        if (StringLen (ctmp) > 800) {
          StringCpy (ctmp + 797, "...");
        }
        product = tmp;
        diff = LabelCopyExtra (tmp, ctmp, buflen, "[", "]");
        buflen -= diff;
        tmp += diff;
        *tmp = '\0';
        tmp++;
        *tmp = '\0';
        MemFree (ctmp);
      }
    }
  } else if (vsp->descr != NULL) {
    if (vsp->bsp != NULL) {
      context = tmp;
      diff = LabelCopy (tmp, "BIOSEQ: ", buflen);
      buflen -= diff;
      tmp += diff;
      if (vsp->suppressContext || vsp->convertGiToAccn) {
        diff = WorstBioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
      *tmp = '\0';
      tmp++;
      *tmp = '\0';
    } else if (vsp->bssp != NULL) {
      context = tmp;
      diff = LabelCopy (tmp, "BIOSEQ-SET: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->suppressContext || vsp->convertGiToAccn) {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
      *tmp = '\0';
      tmp++;
      *tmp = '\0';
    }
  }

  (*errfunc) (severity, errcode, subcode, entityID, itemtype, itemID, accession,
              message, objtype, label, context, location, product, vsp->userdata);
}

#ifdef VAR_ARGS
NLM_EXTERN void CDECL ValidErr (vsp, severity, code1, code2, fmt, va_alist)
     ValidStructPtr vsp;
     int severity;
     int code1;
     int code2;
     const char     *fmt;
     va_dcl
#else
NLM_EXTERN void CDECL ValidErr (ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...)
#endif
{
  va_list           args;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Int2              buflen, diff;
  CharPtr           ctmp, tmp;
  GatherContextPtr  gcp;
  Char              id [64];
  SeqLocPtr         loc = NULL;
  ObjValNodePtr     ovp;
  SeqDescrPtr       sdp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;

  if (vsp == NULL || severity < vsp->cutoff)
    return;

  if (vsp->errbuf == NULL) {
    vsp->errbuf = MemNew (4096);
    if (vsp->errbuf == NULL)
      AbnormalExit (1);
  }
  tmp = vsp->errbuf;

  vsp->errors[severity]++;

#ifdef VAR_ARGS
  va_start (args);
#else
  va_start (args, fmt);
#endif

  gcp = vsp->gcp;
  buflen = 1023;
  vsprintf (tmp, fmt, args);
  while (*tmp != '\0') {
    buflen--;
    tmp++;
  }

  va_end (args);

  if (vsp->errfunc != NULL) {
    CustValErr (vsp, (ErrSev) (severity), code1, code2);
    vsp->errbuf[0] = '\0';
    return;
  }

  if (vsp->justShowAccession) {
    vsp->errbuf[0] = '\0';
    tmp = vsp->errbuf;
    sip = NULL;

    if (vsp->sfp != NULL) {
      sfp = vsp->sfp;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
      }
    } else if (vsp->descr != NULL) {
      sdp = vsp->descr;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
          }
        } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) ovp->idx.parentptr;
          if (bssp != NULL) {
            sep = bssp->seqentry;
            if (sep != NULL) {
              sep = FindNthBioseq (sep, 1);
              if (sep != NULL) {
                bsp = (BioseqPtr) sep->data.ptrvalue;
                if (bsp != NULL) {
                  sip = SeqIdFindWorst (bsp->id);
                }
              }
            }
          }
        }
      }
    } else if (vsp->bsp != NULL) {
      bsp = vsp->bsp;
      sip = SeqIdFindWorst (bsp->id);
    } else if (vsp->bssp != NULL) {
      bssp = vsp->bssp;
      sep = bssp->seqentry;
      if (sep != NULL) {
        sep = FindNthBioseq (sep, 1);
        if (sep != NULL) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
          }
        }
      }
    }

    if (sip != NULL) {
      SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id) - 1);
      diff = LabelCopy (tmp, id, buflen);
      buflen -= diff;
      tmp += diff;
    }

    ErrPostItem ((ErrSev) (severity), code1, code2, "%s", vsp->errbuf);
    vsp->errbuf[0] = '\0';
    return;
  }

  if (vsp->sfp != NULL) {
    diff = LabelCopy (tmp, " FEATURE: ", buflen);
    buflen -= diff;
    tmp += diff;

    diff = FeatDefLabel (vsp->sfp, tmp, buflen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;

    if (vsp->suppressContext) {
      loc = AsnIoMemCopy (vsp->sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ChangeSeqLocToBestID (loc);
      ctmp = SeqLocPrint (loc);
      SeqLocFree (loc);
    } else {
      ctmp = SeqLocPrint (vsp->sfp->location);
    }
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    if (ctmp != NULL) {
      diff = LabelCopyExtra (tmp, ctmp, buflen, " [", "]");
      buflen -= diff;
      tmp += diff;
      MemFree (ctmp);
    }

    if (!vsp->suppressContext) {
      sip = SeqLocId (vsp->sfp->location);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          diff = LabelCopy (tmp, " [", buflen);
          buflen -= diff;
          tmp += diff;

          diff = BioseqLabel (bsp, tmp, buflen, OM_LABEL_BOTH);
          buflen -= diff;
          tmp += diff;

          diff = LabelCopy (tmp, "]", buflen);
          buflen -= diff;
          tmp += diff;
        }
      }
    }
    if (vsp->sfp->product != NULL) {
      if (vsp->suppressContext) {
        loc = AsnIoMemCopy (vsp->sfp->product, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        ChangeSeqLocToBestID (loc);
        ctmp = SeqLocPrint (loc);
        SeqLocFree (loc);
      } else {
        ctmp = SeqLocPrint (vsp->sfp->product);
      }
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      if (ctmp != NULL) {
        diff = LabelCopyExtra (tmp, ctmp, buflen, " -> [", "]");
        buflen -= diff;
        tmp += diff;
        MemFree (ctmp);
      }
    }
  } else if (vsp->descr != NULL) {
    diff = LabelCopy (tmp, " DESCRIPTOR: ", buflen);
    buflen -= diff;
    tmp += diff;

    diff = SeqDescLabel (vsp->descr, tmp, buflen, OM_LABEL_BOTH);
    buflen -= diff;
    tmp += diff;
  }

  /*
     if (vsp->suppressContext)
     {
     }
     else */
  if (vsp->sfp == NULL) {       /* sfp adds its own context */
    if (vsp->bsp != NULL) {
      diff = LabelCopy (tmp, " BIOSEQ: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->bsp == NULL) {
        diff = LabelCopy (tmp, "??", buflen);
      } else if (vsp->suppressContext) {
        diff = WorstBioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqLabel (vsp->bsp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
    } else if (vsp->bssp != NULL) {
      diff = LabelCopy (tmp, " BIOSEQ-SET: ", buflen);
      buflen -= diff;
      tmp += diff;

      if (vsp->suppressContext) {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_CONTENT);
      } else {
        diff = BioseqSetLabel (vsp->bssp, tmp, buflen, OM_LABEL_BOTH);
      }
      buflen -= diff;
      tmp += diff;
    }
  }

  ErrPostItem ((ErrSev) (severity), code1, code2, "%s", vsp->errbuf);
  vsp->errbuf[0] = '\0';
}

/*****************************************************************************
*
*   Valid1GatherProc(gcp)
*     top level gather callback
*     dispatches to other levels
*
*****************************************************************************/
static Boolean Valid1GatherProc (GatherContextPtr gcp)
{
  ValidStructPtr     vsp;
  AnnotDescrPtr      desc;
  SeqAnnotPtr        sap;
  ObjectIdPtr        oip;
  Boolean            is_blast_align;
  Int2               limit;
  SeqFeatPtr         sfp;
  ValNodePtr         sdp;
  SeqGraphPtr        sgp;
  BioSourcePtr       biop;
  PubdescPtr         pdp;
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  Char               buf [64];
  Char               tmp [64];
  SeqMgrFeatContext  context;

  vsp = (ValidStructPtr) (gcp->userdata);
  vsp->gcp = gcp;               /* needed for ValidErr */

  limit = vsp->validationLimit;

  switch (gcp->thistype) {
  case OBJ_BIOSEQ:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_INST) {
        ValidateBioseqInst (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_CONTEXT) {
        ValidateBioseqContext (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_INST) {
        ValidateBioseqHist (gcp);
      }
      if (limit == VALIDATE_ALL || limit == VALIDATE_GRAPH) {
        ValidateGraphsOnBioseq (gcp);
      }
    }
    break;
  case OBJ_BIOSEQSET:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_SET) {
        ValidateBioseqSet (gcp);
      }
    }
    break;
  case OBJ_SEQANNOT:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL) {
        sap = (SeqAnnotPtr) gcp->thisitem;
        if (sap != NULL) {
          if (sap->type == 2) {
            is_blast_align = FALSE;
            desc = NULL;
            while ((desc = ValNodeFindNext (sap->desc, desc, Annot_descr_user)) != NULL) {
              if (desc->data.ptrvalue != NULL) {
                oip = ((UserObjectPtr) desc->data.ptrvalue)->type;
                if (oip != NULL && StringCmp (oip->str, "Blast Type") == 0) {
                  is_blast_align = TRUE;
                }
              }
            }
            if (is_blast_align) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_ALIGN_BlastAligns, "Record contains BLAST alignments");
            }
          }
          if (sap->type == 4) {
            vsp->bssp = NULL;
            vsp->bsp = NULL;
            vsp->descr = NULL;
            vsp->sfp = NULL;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_ANNOT_AnnotIDs, "Record contains Seq-annot.data.ids");
          }
          if (sap->type == 5) {
            vsp->bssp = NULL;
            vsp->bsp = NULL;
            vsp->descr = NULL;
            vsp->sfp = NULL;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_ANNOT_AnnotLOCs, "Record contains Seq-annot.data.locs");
          }
        }
      }
    }
    break;
  case OBJ_SEQFEAT:
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_FEAT) {
        ValidateSeqFeat (gcp);
        sfp = (SeqFeatPtr) (gcp->thisitem);
        if (sfp != NULL) {
          if (sfp->data.choice == SEQFEAT_BIOSRC) {
            biop = (BioSourcePtr) sfp->data.value.ptrvalue;
            ValidateBioSource (vsp, gcp, biop, sfp, NULL);
          }
          if (sfp->data.choice == SEQFEAT_PUB) {
            pdp = (PubdescPtr) sfp->data.value.ptrvalue;
            ValidatePubdesc (vsp, gcp, pdp);
          }
          if (sfp->cit != NULL) {
            ValidateSfpCit (vsp, gcp, sfp);
          }
          if (vsp->useSeqMgrIndexes) {
            if (SeqMgrGetDesiredFeature (gcp->entityID, NULL, 0, 0, sfp, &context) == NULL) {
              StringCpy (buf, "?");
              bsp = vsp->bsp;
              if (bsp != NULL) {
                SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
              }
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_UnindexedFeature, "Feature is not indexed on Bioseq %s", buf);
            } else {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                sip = SeqLocId (sfp->location);
                if (sip != NULL && sip->choice != SEQID_GI && sip->choice != SEQID_GIBBSQ && sip->choice != SEQID_GIBBMT) {
                  SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
                  for (sip = bsp->id; sip != NULL; sip = sip->next) {
                    if (sip->choice == SEQID_GI || sip->choice == SEQID_GIBBSQ || sip->choice == SEQID_GIBBMT) continue;
                    SeqIdWrite (sip, tmp, PRINTID_FASTA_SHORT, sizeof (tmp) - 1);
                    if (StringICmp (buf, tmp) != 0) continue;
                    if (StringCmp (buf, tmp) == 0) continue;
                    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_FeatureSeqIDCaseDifference,
                              "Sequence identifier in feature location differs in capitalization with identifier on Bioseq");
                  }
                }
              }
            }
          }
        }
      }
    }
    if (limit == VALIDATE_ALL || limit == VALIDATE_FEAT) {
      SpellCheckSeqFeat (gcp);
    }
    break;
  case OBJ_SEQGRAPH :
    if (!vsp->onlyspell) {
      if (limit == VALIDATE_ALL || limit == VALIDATE_GRAPH) {
        sgp = (SeqGraphPtr) gcp->thisitem;
        if (sgp != NULL) {
          if (StringICmp (sgp->title, "Phrap Quality") == 0 ||
              StringICmp (sgp->title, "Phred Quality") == 0 ||
              StringICmp (sgp->title, "Gap4") == 0) {
            if (sgp->flags[2] == 3) {
              sip = SeqLocId (sgp->loc);
              if (sip != NULL) {
                if (BioseqFindCore (sip) == NULL) {
                  SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBioseqId, "Bioseq not found for Graph location %s", buf);
                }
              }
            }
          }
        }
      }
    }
    break;
  case OBJ_SEQDESC:
    if (limit == VALIDATE_ALL || limit == VALIDATE_DESC) {
      SpellCheckSeqDescr (gcp);
                          /**
              ValidateSeqDescr (gcp);
              **/
      sdp = (ValNodePtr) (gcp->thisitem);
      if (sdp != NULL) {
        if (sdp->choice == Seq_descr_source) {
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          ValidateBioSource (vsp, gcp, biop, NULL, sdp);
        }
        if (sdp->choice == Seq_descr_pub) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
          ValidatePubdesc (vsp, gcp, pdp);
          LookForMultiplePubs (vsp, gcp, sdp);
        }
        if (sdp->choice == Seq_descr_mol_type) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "MolType descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_modif) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Modif descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_method) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Method descriptor is obsolete");
        }
        if (sdp->choice == Seq_descr_org) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "OrgRef descriptor is obsolete");
        }
      }
    }
    break;
  default:
    break;

  }
  return TRUE;
}


static void DiscrepanciesToValidationErrs (ValNodePtr discrepancy_list, Uint4 item_type, ValidStructPtr vsp, int severity, int code1, int code2, char *msg)
{
  ValNodePtr vnp, obj;
  ClickableItemPtr cip;

  if (discrepancy_list == NULL || vsp == NULL) {
    return;
  }
  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      if (cip->clickable_item_type == item_type) {
        if (cip->item_list == NULL) {
          DiscrepanciesToValidationErrs (cip->subcategories, item_type, vsp, severity, code1, code2, msg);
        } else {
          for (obj = cip->item_list; obj != NULL; obj = obj->next) {
            if (obj->choice == OBJ_SEQFEAT) {
              vsp->sfp = obj->data.ptrvalue;
              ValidErr (vsp, severity, code1, code2, msg);
              vsp->sfp = NULL;
            }
          }
        }
      }
    }
  }  
}


static void ValidateGeneLocusTags (SeqEntryPtr sep, ValidStructPtr vsp)
{
  ValNode vn;
  ValNodePtr discrepancy_list = NULL;

  if (sep == NULL || vsp == NULL) {
    return;
  }

  vn.choice = 0;
  vn.data.ptrvalue = sep;
  vn.next = NULL;
  
  AddDiscrepanciesForMissingOrNonUniqueGeneLocusTagsEx (&discrepancy_list, &vn, TRUE);

  DiscrepanciesToValidationErrs (discrepancy_list, DISC_GENE_MISSING_LOCUS_TAG, vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingGeneLocusTag, "Missing gene locus tag");

  discrepancy_list = FreeClickableList (discrepancy_list);
}


static void LookForAnyPubAndOrg (SeqEntryPtr sep, BoolPtr no_pub, BoolPtr no_biosrc)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqAnnotPtr     sap = NULL;
  ValNodePtr      sdp = NULL;
  SeqFeatPtr      sfp;
  SeqEntryPtr     tmp;

  if (sep == NULL || no_pub == NULL || no_biosrc == NULL)
    return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL)
      return;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      LookForAnyPubAndOrg (tmp, no_pub, no_biosrc);
    }
    sap = bssp->annot;
    sdp = bssp->descr;
  } else
    return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          *no_pub = FALSE;
        } else if (sfp->data.choice == SEQFEAT_BIOSRC) {
          *no_biosrc = FALSE;
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      *no_pub = FALSE;
    } else if (sdp->choice == Seq_descr_source) {
      *no_biosrc = FALSE;
    }
    sdp = sdp->next;
  }
}

typedef struct ftprob {
  Uint4    num_misplaced_features;
  Uint4    num_archaic_locations;
  Uint4    num_archaic_products;
  Uint4    num_misplaced_graphs;
  Uint4    num_gene_feats;
  Uint4    num_gene_xrefs;
  Uint4    num_tpa_with_hist;
  Uint4    num_tpa_without_hist;
  Boolean  has_gi;
  Boolean  loc_has_gi;
  Boolean  loc_has_just_accn;
  Boolean  loc_has_accn_ver;
  Boolean  prod_has_gi;
  Boolean  prod_has_just_accn;
  Boolean  prod_has_accn_ver;
} FeatProb, PNTR FeatProbPtr;

static void CheckFeatPacking (BioseqPtr bsp, SeqFeatPtr sfp, Uint4Ptr num_misplaced_features)
{
  SeqAnnotPtr     sap;
  BioseqSetPtr    bssp, parent;
  BioseqPtr       par;

  if (sfp->idx.parenttype == OBJ_SEQANNOT) {
    sap = (SeqAnnotPtr) sfp->idx.parentptr;
    if (sap == NULL)
      return;
    if (sap->idx.parenttype == OBJ_BIOSEQ) {
      /* if feature packaged on bioseq, must be target bioseq */
      par = (BioseqPtr) sap->idx.parentptr;
      if (par != bsp && SeqMgrGetParentOfPart (par, NULL) != bsp) {
        /* generated gap feature is an exception */
        if (par == NULL || par->id != NULL) {
          (*num_misplaced_features)++;
        }
      }
      return;
    }
    if (sap->idx.parenttype == OBJ_BIOSEQSET) {
      /* if feature packaged on set, set must contain bioseq */
      bssp = (BioseqSetPtr) sap->idx.parentptr;
      if (bssp == NULL)
        return;
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bsp->idx.parentptr;
        while (parent != NULL) {
          if (parent == bssp)
            return;
          if (parent->idx.parenttype != OBJ_BIOSEQSET) {
            (*num_misplaced_features)++;
            return;
          }
          parent = (BioseqSetPtr) parent->idx.parentptr;
        }
        (*num_misplaced_features)++;
      }
    }
  }
}

static Boolean IdIsArchaic (SeqIdPtr sip)

{
  BioseqPtr  bsp;
  DbtagPtr   dbt;
  SeqIdPtr   id;

  if (sip == NULL) return FALSE;
  if (sip->choice != SEQID_LOCAL && sip->choice != SEQID_GENERAL) return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return FALSE;
  for (id = bsp->id; id != NULL; id = id->next) {
    switch (id->choice) {
      case SEQID_GENERAL :
        if (sip->choice == SEQID_LOCAL) {
          dbt = (DbtagPtr) id->data.ptrvalue;
          if (dbt != NULL && !IsSkippableDbtag(dbt)) {
            return TRUE;
          }
        }
        break;
      case SEQID_GI :
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_PATENT :
      case SEQID_OTHER :
      case SEQID_DDBJ :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
      case SEQID_GPIPE :
        return TRUE;
      default :
        break;
    }
  }
  return FALSE;
}

static void CheckFeatLocAndProd (SeqFeatPtr sfp, FeatProbPtr fpp)

{
  SeqLocPtr  slp;

  if (sfp == NULL || fpp == NULL) return;
  if (sfp->product != NULL && IdIsArchaic (SeqLocId (sfp->product))) {
    (fpp->num_archaic_products)++;
  }
  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) {
    if (IdIsArchaic (SeqLocId (slp))) {
      (fpp->num_archaic_locations)++;
      return;
    }
    slp = SeqLocFindNext (sfp->location, slp);
  }
}

static void CheckGraphPacking (SeqGraphPtr sgp, Pointer userdata)

{
  BioseqPtr    bsp;
  FeatProbPtr  fpp;
  SeqAnnotPtr  sap;
  BioseqPtr    par;

  if (sgp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;
  bsp = BioseqFindFromSeqLoc (sgp->loc);
  if (sgp->idx.parenttype == OBJ_SEQANNOT) {
    sap = (SeqAnnotPtr) sgp->idx.parentptr;
    if (sap == NULL) return;
    if (sap->idx.parenttype == OBJ_BIOSEQ) {
      /* if graph packaged on bioseq, must be target bioseq */
      par = (BioseqPtr) sap->idx.parentptr;
      if (par != bsp && SeqMgrGetParentOfPart (par, NULL) != bsp) {
        (fpp->num_misplaced_graphs)++;
      }
      return;
    }
    (fpp->num_misplaced_graphs)++;
  }
}

static Boolean LIBCALLBACK CountMisplacedFeatures (BioseqPtr bsp, SeqMgrBioseqContextPtr bcontext)
{
  FeatProbPtr     fpp;
  SeqFeatPtr      sfp;
  SeqMgrFeatContext fcontext;

  fpp = (FeatProbPtr) bcontext->userdata;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    CheckFeatPacking (bsp, sfp, &(fpp->num_misplaced_features));
    CheckFeatLocAndProd (sfp, fpp);
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  return TRUE;
}

static void CountGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  FeatProbPtr  fpp;
  GeneRefPtr   grp;

  if (sfp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  if (sfp->data.choice == SEQFEAT_GENE) {
    (fpp->num_gene_feats)++;
  }

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  (fpp->num_gene_xrefs)++;
}

static void CountSfpLocIdTypes (SeqIdPtr sip, Pointer userdata)

{
  FeatProbPtr   fpp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  switch (sip->choice) {
    case SEQID_GI :
      fpp->loc_has_gi = TRUE;
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          if (tsip->version < 1) {
            fpp->loc_has_just_accn = TRUE;
          } else {
            fpp->loc_has_accn_ver = TRUE;
          }
        }
      }
      break;
    default :
      break;
  }
}

static void CountSfpProdIdTypes (SeqIdPtr sip, Pointer userdata)

{
  FeatProbPtr   fpp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;

  switch (sip->choice) {
    case SEQID_GI :
      fpp->prod_has_gi = TRUE;
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          if (tsip->version < 1) {
            fpp->prod_has_just_accn = TRUE;
          } else {
            fpp->prod_has_accn_ver = TRUE;
          }
        }
      }
      break;
    default :
      break;
  }
}

static void CountFeatLocIdTypes (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL || userdata == NULL) return;

  VisitSeqIdsInSeqLoc (sfp->location, userdata, CountSfpLocIdTypes);
  VisitSeqIdsInSeqLoc (sfp->product, userdata, CountSfpProdIdTypes);
}

NLM_EXTERN Boolean HasTpaUserObject (BioseqPtr bsp)

{
  SeqMgrDescContext  context;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  ValNodePtr         vnp;

  if (bsp == NULL) return FALSE;
  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (vnp != NULL) {
    uop = (UserObjectPtr) vnp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "TpaAssembly") == 0) return TRUE;
    }
    vnp = SeqMgrGetNextDescriptor (bsp, vnp, Seq_descr_user, &context);
  }
  return FALSE;
}

static void CheckTpaHist (BioseqPtr bsp, Pointer userdata)

{
  FeatProbPtr  fpp;
  SeqHistPtr   shp;
  SeqIdPtr     sip;

  if (bsp == NULL || userdata == NULL) return;
  fpp = (FeatProbPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      fpp->has_gi = TRUE;
    }
  }
  if (! HasTpaUserObject (bsp)) return;
  shp = bsp->hist;
  if (shp != NULL && shp->assembly != NULL) {
    (fpp->num_tpa_with_hist)++;
  } else {
    (fpp->num_tpa_without_hist)++;
  }
}

static Boolean IsNoncuratedRefSeq (BioseqPtr bsp, ErrSev *sev)

{
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
            StringNCmp (tsip->accession, "NP_", 3) == 0 ||
            StringNCmp (tsip->accession, "NG_", 3) == 0 ||
            StringNCmp (tsip->accession, "NR_", 3) == 0) {
          *sev = SEV_WARNING;
          return FALSE;
        }
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean IsGpipe (BioseqPtr bsp)

{
  SeqIdPtr  sip;

  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GPIPE) return TRUE;
  }
  return FALSE;
}

typedef struct vfcdata {
  ValNodePtr      uids;
  ValNodePtr      unpub;
  ValNodePtr      publshd;
  ValNodePtr      serial;
  ValidStructPtr  vsp;
} VfcData, PNTR VfcPtr;

static Boolean SkipSerialOrUIDPub (ValNodePtr vnp)

{
  CitGenPtr  cgp;

  if (vnp == NULL || vnp->next == NULL) return FALSE;
  if (vnp->choice == PUB_Muid || vnp->choice == PUB_Muid) return TRUE;
  if (vnp->choice != PUB_Gen) return FALSE;
  cgp = (CitGenPtr) vnp->data.ptrvalue;
  if (cgp == NULL) return FALSE;
  if (StringNICmp ("BackBone id_pub", cgp->cit, 15) == 0) return FALSE;
  if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) return TRUE;
  return FALSE;
}

static void MakePubTags (PubdescPtr pdp, Pointer userdata)

{
  Char        buf [1024];
  CitGenPtr   cgp;
  Int4        muid = 0, pmid = 0;
  VfcPtr      vfp;
  ValNodePtr  vnp, tmp;

  if (pdp == NULL || userdata == NULL) return;
  vfp = (VfcPtr) userdata;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid) {
      muid = vnp->data.intvalue;
    } else if (vnp->choice == PUB_PMid) {
      pmid = vnp->data.intvalue;
    } else if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL && cgp->serial_number > 0) {
        tmp = ValNodeNew (NULL);
        if (tmp != NULL) {
          tmp->data.intvalue = (Int4) cgp->serial_number;
          tmp->next = vfp->serial;
          vfp->serial = tmp;
        }
      }
    }
  }

  if (pmid != 0) {
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = 1;
      vnp->data.intvalue = pmid;
      vnp->next = vfp->uids;
      vfp->uids = vnp;
    }
  }
  if (muid != 0) {
    vnp = ValNodeNew (NULL);
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = muid;
      vnp->next = vfp->uids;
      vfp->uids = vnp;
    }
  }

  vnp = pdp->pub;
  while (vnp != NULL && SkipSerialOrUIDPub (vnp)) {
    vnp = vnp->next;
  }
  if (vnp != NULL && PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    tmp = ValNodeCopyStr (NULL, 0, buf);
    if (tmp != NULL) {
      if (pmid != 0 || muid != 0) {
        tmp->next = vfp->publshd;
        vfp->publshd = tmp;
      } else {
        tmp->next = vfp->unpub;
        vfp->unpub = tmp;
      }
    }
  }
}

static void CheckOneCit (SeqFeatPtr sfp, ValNodePtr ppr, VfcPtr vfp)

{
  Char              buf [1024];
  GatherContextPtr  gcp;
  size_t            len, lgth;
  CharPtr           str;
  Int4              uid;
  ValNodePtr        vnp;
  ValidStructPtr    vsp;

  if (sfp == NULL || ppr == NULL || vfp == NULL) return;
  vsp = vfp->vsp;
  if (vsp == NULL) return;
  gcp = vsp->gcp;

  if (gcp != NULL) {
    gcp->entityID = sfp->idx.entityID;
    gcp->itemID = sfp->idx.itemID;
    gcp->thistype = OBJ_SEQFEAT;
  }
  vsp->sfp = sfp;

  if (ppr->choice == PUB_PMid || ppr->choice == PUB_Muid) {
    uid = ppr->data.intvalue;
    for (vnp = vfp->uids; vnp != NULL; vnp = vnp->next) {
      if (uid == vnp->data.intvalue) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
              "Citation on feature refers to uid [%ld] not on a publication in the record", (long) uid);
    vsp->sfp = NULL;

  } else if (ppr->choice == PUB_Equiv) {
    return;
  
  } else {
    PubLabelUnique (ppr, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
    lgth = StringLen (buf);
    if (lgth > 0 && buf [lgth - 1] == '>') {
      buf [lgth - 1] = '\0';
     lgth--;
    }
    for (vnp = vfp->unpub; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = MIN (lgth, StringLen (str));
      if (StringNICmp (str, buf, len) == 0) return;
    }
    for (vnp = vfp->publshd; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = MIN (lgth, StringLen (str));
      if (StringNICmp (str, buf, len) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
                  "Citation on feature needs to be updated to published uid");
        vsp->sfp = NULL;
        return;
      }
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCitationProblem,
              "Citation on feature refers to a publication not in the record");
    vsp->sfp = NULL;
  }
}

static void CheckFeatCits (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr  ppr, vnp;
  VfcPtr      vfp;

  if (sfp == NULL || sfp->cit == NULL || userdata == NULL) return;
  vfp = (VfcPtr) userdata;

  vnp = sfp->cit;
  for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
    CheckOneCit (sfp, ppr, vfp);
  }
}

static void CheckForCollidingSerials (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  ValNodePtr list
)

{
  Int4        curr, last;
  Uint2       olditemtype = 0;
  Uint4       olditemid = 0;
  ValNodePtr  vnp, vnp_next;

  if (vsp == NULL || gcp == NULL || list == NULL) return;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;
  gcp->itemID = 0;
  gcp->thistype = 0;

  last = (Int4) list->data.intvalue;
  for (vnp = list->next; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    curr = (Int4) vnp->data.intvalue;
    if (last == curr) {
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_CollidingSerialNumbers,
                "Multiple publications have serial number %ld", (long) curr);
      while (vnp != NULL && vnp->data.intvalue == last) {
        vnp = vnp->next;
      }
      if (vnp == NULL) {
        vnp_next = NULL;
      } else {
        last = vnp->data.intvalue;
        vnp_next = vnp->next;
      }
    } else {
      last = curr;
    }
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
}

static void ValidateFeatCits (SeqEntryPtr sep, ValidStructPtr vsp)

{
  GatherContext  gc;
  VfcData        vfd;

  if (vsp == NULL || sep == NULL) return;
  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  MemSet ((Pointer) &vfd, 0, sizeof (VfcData));
  vfd.vsp = vsp;

  VisitPubdescsInSep (sep, (Pointer) &vfd, MakePubTags);

  VisitFeaturesInSep (sep, (Pointer) &vfd, CheckFeatCits);

  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  vfd.serial = ValNodeSort (vfd.serial, SortByIntvalue);
  CheckForCollidingSerials (vsp, vsp->gcp, vfd.serial);

  ValNodeFree (vfd.uids);
  ValNodeFreeData (vfd.unpub);
  ValNodeFreeData (vfd.publshd);
  ValNodeFree (vfd.serial);
}

static void ValidateFeatIDs (Uint2 entityID, ValidStructPtr vsp)

{
  SMFidItemPtr PNTR  array;
  BioseqExtraPtr     bspextra;
  SMFeatItemPtr      feat;
  GatherContext      gc;
  GatherContextPtr   gcp;
  SMFidItemPtr       item;
  Int4               j;
  CharPtr            last = NULL;
  Int4               num;
  ObjMgrDataPtr      omdp;
  SeqFeatPtr         sfp;

  if (entityID < 1 || vsp == NULL) return;
  omdp = ObjMgrGetData (entityID);
  if (omdp == NULL) return;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;
  array = bspextra->featsByFeatID;
  num = bspextra->numfids;
  if (array == NULL || num < 1) return;

  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));

  for (j = 0; j < num; j++) {
    item = array [j];
    if (item == NULL) continue;
    if (StringDoesHaveText (last)) {
      if (StringICmp (item->fid, last) == 0) {
        feat = item->feat;
        if (feat == NULL) continue;
        sfp = feat->sfp;
        if (sfp == NULL) continue;
        gcp = &gc;
        gcp->entityID = sfp->idx.entityID;
        gcp->itemID = sfp->idx.itemID;
        gcp->thistype = OBJ_SEQFEAT;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_CollidingFeatureIDs,
                  "Colliding feature ID %s", last);
      }
    }
    last = item->fid;
  }
}

typedef struct vsicdata {
  ValidStructPtr  vsp;
  ValNodePtr      headid;
  ValNodePtr      tailid;
} VsicData, PNTR VsicDataPtr;

static Boolean IsNCBIFileID (SeqIdPtr sip)
{
  DbtagPtr dbt;

  if (sip == NULL || sip->choice != SEQID_GENERAL) return FALSE;
  dbt = (DbtagPtr) sip->data.ptrvalue;
  if (dbt == NULL) return FALSE;
  if (StringCmp (dbt->db, "NCBIFILE") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void CaptureTextSeqIDs (BioseqPtr bsp, Pointer userdata)

{
  Char         buf [200];
  SeqIdPtr     sip;
  VsicDataPtr  vdp;
  ValNodePtr   vnp;

  if (bsp == NULL || userdata == NULL) return;
  vdp = (VsicDataPtr) userdata;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI || sip->choice == SEQID_GIBBSQ || sip->choice == SEQID_GIBBMT) continue;
    if (IsNCBIFileID (sip)) continue;
    SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    vnp = ValNodeCopyStr (&(vdp->tailid), 0, buf);
    if (vdp->headid == NULL) {
      vdp->headid = vnp;
    }
    vdp->tailid = vnp;
  }
}

static ValNodePtr UniqueValNodeCaseSensitive (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringCmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

static void ValidateSeqIdCase (SeqEntryPtr sep, ValidStructPtr vsp)

{
  CharPtr           curr;
  GatherContext     gc;
  GatherContextPtr  gcp;
  CharPtr           prev;
  VsicData          vd;
  ValNodePtr        vnp;

  if (vsp == NULL || sep == NULL) return;

  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  MemSet ((Pointer) &vd, 0, sizeof (VsicData));

  gcp = &gc;
  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  vd.vsp = vsp;

  VisitBioseqsInSep (sep, (Pointer) &vd, CaptureTextSeqIDs);
  vd.headid = ValNodeSort (vd.headid, SortVnpByString);
  vd.headid = UniqueValNodeCaseSensitive (vd.headid);

  curr = NULL;
  prev = NULL;
  for (vnp = vd.headid; vnp != NULL; vnp = vnp->next, prev = curr) {
    curr = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (curr)) continue;
    if (StringHasNoText (prev)) continue;
    if (StringICmp (curr, prev) != 0) continue;
    if (StringCmp (curr, prev) == 0) continue;
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_CaseDifferenceInSeqID,
              "Sequence identifier differs only by case - %s and %s", curr, prev);
  }

  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;

  ValNodeFreeData (vd.headid);
}

static void LookForNC (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr       is_ncp;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  if (bsp == NULL || userdata == NULL) return;
  is_ncp = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        /*
        if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
          *is_ncp = TRUE;
        }
        */
        *is_ncp = TRUE; /* any refseq now drops pubdesc message severity */
      }
    }
  }
}

static void LookForGPS (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;
  BoolPtr       is_gpsp;

  if (sep == NULL || data == NULL) return;
  is_gpsp = (BoolPtr) data;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      *is_gpsp = TRUE;
    }
  }
}

static void LookForNonGPS (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;
  BoolPtr       is_ngpsp;

  if (sep == NULL || data == NULL) return;
  is_ngpsp = (BoolPtr) data;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == BioseqseqSet_class_mut_set ||
        bssp->_class == BioseqseqSet_class_pop_set ||
        bssp->_class == BioseqseqSet_class_phy_set ||
        bssp->_class == BioseqseqSet_class_eco_set ||
        bssp->_class == BioseqseqSet_class_wgs_set) {
      *is_ngpsp = TRUE;
    }
  }
}

static void LookForEmblDdbj (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr   is_ed;
  SeqIdPtr  sip;

  if (bsp == NULL || userdata == NULL) return;
  is_ed = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {
      *is_ed = TRUE;
    }
  }
}

static void LookForGEDseqID (BioseqPtr bsp, Pointer userdata)
{
  BoolPtr         isGEDPtr;
  SeqIdPtr        sip;

  isGEDPtr = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
    case SEQID_GENBANK:
    case SEQID_EMBL:
    case SEQID_DDBJ:
    case SEQID_TPG:
    case SEQID_TPE:
    case SEQID_TPD:
      *isGEDPtr = TRUE;
      return;
    default:
      break;
    }
  }
}

static void LookForLclGnl (BioseqPtr bsp, Pointer userdata)

{
  Boolean   has_lcl_gnl = FALSE;
  Boolean   has_others = FALSE;
  BoolPtr   is_lcl_gnl_P;
  SeqIdPtr  sip;

  if (bsp == NULL || userdata == NULL) return;
  is_lcl_gnl_P = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL || sip->choice == SEQID_GENERAL) {
      has_lcl_gnl = TRUE;
    } else {
      has_others = TRUE;
    }
  }
  if (has_others) {
    *is_lcl_gnl_P = FALSE;
    return;
  }
  if (has_lcl_gnl) {
    *is_lcl_gnl_P = TRUE;
  }
}

static void LookForProteinGnl (BioseqPtr bsp, Pointer userdata)

{
  DbtagPtr  dbt;
  Boolean   has_gnl_prot = FALSE;
  BoolPtr   is_gnl_prot_P;
  SeqIdPtr  sip;

  if (bsp == NULL || userdata == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  is_gnl_prot_P = (BoolPtr) userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL) {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag (dbt)) continue;
      has_gnl_prot = TRUE;
    }
  }
  if (has_gnl_prot) {
    *is_gnl_prot_P = TRUE;
  }
}

static void LookForHTG (SeqDescrPtr sdp, Pointer userdata)

{
  BoolPtr     is_htgp;
  MolInfoPtr  mip;

  if (sdp == NULL || userdata == NULL) return;
  if (sdp->choice != Seq_descr_molinfo) return;

  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  if (mip->tech == MI_TECH_htgs_1 ||
      mip->tech == MI_TECH_htgs_2 ||
      mip->tech == MI_TECH_htgs_3 ||
      mip->tech == MI_TECH_htgs_0) {

    is_htgp = (BoolPtr) userdata;
    *is_htgp = TRUE; /* any htg now drops citsub missing affil message severity */
  }
}

static void LookForBarcode (SeqDescrPtr sdp, Pointer userdata)

{
  BoolPtr     is_barcode;
  MolInfoPtr  mip;

  if (sdp == NULL || userdata == NULL) return;
  if (sdp->choice != Seq_descr_molinfo) return;

  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  if (mip->tech == MI_TECH_barcode) {

    is_barcode = (BoolPtr) userdata;
    *is_barcode = TRUE; /* any barcode now drops bad country code message severity */
  }
}

static void LookForSMUPD (SeqDescrPtr sdp, Pointer userdata)

{
  BoolPtr        is_smupdp;
  UserObjectPtr  uop;

  if (sdp == NULL || userdata == NULL) return;
  if (sdp->choice != Seq_descr_user) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;

  if (StringICmp (uop->_class, "SMART_V1.0") == 0) {

    is_smupdp = (BoolPtr) userdata;
    *is_smupdp = TRUE;
  }
}

static void SetPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  AuthListPtr    alp;
  Char           buf [2048];
  CitGenPtr      cgp;
  CharPtr        consortium, str, tmp;
  ValNodePtr     vnp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    alp = GetAuthListPtr (pdp, NULL);
    if (alp != NULL) {
      consortium = NULL;
      str = GetAuthorsString (GENBANK_FMT, alp, &consortium, NULL, NULL);
      tmp = MemNew (StringLen (buf) + StringLen (str) + StringLen (consortium) + 10);
      if (tmp != NULL) {
        StringCpy (tmp, buf);
        if (StringDoesHaveText (str)) {
          StringCat (tmp, "; ");
          StringCat (tmp, str);
        }
        if (StringDoesHaveText (consortium)) {
          StringCat (tmp, "; ");
          StringCat (tmp, consortium);
        }
        ovp->idx.scratch = tmp;
      }
      MemFree (str);
      MemFree (consortium);
    }
  }
}

static void ClearPubScratchData (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  ovp->idx.scratch = MemFree (ovp->idx.scratch);
}

static ValNodePtr SetUpValidateGeneticCodes (void)

{
  Char            ch;
  GeneticCodePtr  codes;
  GeneticCodePtr  gcp;
  ValNodePtr      gencodelist = NULL;
  Int2            i;
  Int4            id;
  Int2            j;
  Char            name [64];
  CharPtr         ptr;
  Char            str [256];
  ValNodePtr      tmp;

  codes = GeneticCodeTableLoad ();
  if (codes != NULL) {
    for (gcp = codes; gcp != NULL; gcp = gcp->next) {
      id = 0;
      str [0] = '\0';
      for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
        switch (tmp->choice) {
          case 1 :
            if (StringLen (str) < 1) {
              StringNCpy_0 (str, (CharPtr) tmp->data.ptrvalue, sizeof (str));
              ptr = str;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '/') {
                  *ptr = '-';
                }
                ptr++;
                ch = *ptr;
              }
            }
            break;
          case 2 :
            id = tmp->data.intvalue;
            break;
          default :
            break;
        }
      }
      if (id != 7 && id != 8) {
        if (id > 0 /* && id < 30 */ ) {
          i = 0;
          if (StringLen (str + i) > 0) {
            ch = str [i];
            while (ch == ' ' || ch == ';') {
              i++;
              ch = str [i];
            }
            j = 0;
            ch = str [i + j];
            while (ch != '\0' && ch != ';') {
              name [j] = ch;
              j++;
              ch = str [i + j];
            }
            name [j] = '\0';
            i += j;
            if (ch == ';') {
              StringCat (name, ", etc.");
            }
            ValNodeCopyStr (&gencodelist, (Uint1) id, name);
          }
        }
      }
    }
  }
  return gencodelist;
}

typedef struct frd {
  ValidStructPtr    vsp;
  GatherContextPtr  gcp;
  /*
  CharPtr           string;
  */
} FindRepData, PNTR FindRepPtr;

static void FindRepValidate (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Pointer userdata)

{
  FindRepPtr        frp;
  GatherContextPtr  gcp;
  ValidStructPtr    vsp;

  frp = (FindRepPtr) userdata;
  vsp = frp->vsp;
  gcp = frp->gcp;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  ValidErr (vsp, SEV_ERROR, ERR_GENERIC_EmbeddedScript, "Script tag found in item");
}

static CharPtr findrepstrs [] = {
  "<script", "<object", "<applet", "<embed", "<form", "javascript:", "vbscript:", NULL
};

typedef struct vvmdata {
  Int2        num_mrnas;
  Boolean     accounted_for;
  Boolean     products_unique;
  Boolean     featid_matched;
  Boolean     feat_touches_gap;
  SeqFeatPtr  nearbygene;
  SeqFeatPtr  nearbycds;
  SeqFeatPtr  nearbymrna;
} VvmData, PNTR VvmDataPtr;

static void AddScratchToFeatures (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  sfp->idx.scratch = (Pointer) MemNew (sizeof (VvmData));
}

static void ClearScratchOnFeatures (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  sfp->idx.scratch = MemFree (sfp->idx.scratch);
}

static void SetupFeatureScratchData (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqFeatPtr         currcds = NULL, currmrna = NULL, currgene = NULL;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;
  VvmDataPtr         vdp;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_GENE :
        currgene = sfp;
        break;
      case FEATDEF_CDS :
        currcds = sfp;
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
          if (vdp->nearbymrna == NULL) {
            vdp->nearbymrna = currmrna;
          }
        }
        if (currgene != NULL) {
          vdp = (VvmDataPtr) currgene->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbycds == NULL) {
              vdp->nearbycds = currcds;
            }
          }
        }
        if (currmrna != NULL) {
          vdp = (VvmDataPtr) currmrna->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbycds == NULL) {
              vdp->nearbycds = currcds;
            }
          }
        }
        break;
      case FEATDEF_mRNA :
        currmrna = sfp;
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
        }
        if (currgene != NULL) {
          vdp = (VvmDataPtr) currgene->idx.scratch;
          if (vdp != NULL) {
            if (vdp->nearbymrna == NULL) {
              vdp->nearbymrna = currmrna;
            }
          }
        }
        break;
      default :
        vdp = (VvmDataPtr) sfp->idx.scratch;
        if (vdp != NULL) {
          if (vdp->nearbygene == NULL) {
            vdp->nearbygene = currgene;
          }
          if (vdp->nearbymrna == NULL) {
            vdp->nearbymrna = currmrna;
          }
          if (vdp->nearbycds == NULL) {
            vdp->nearbycds = currcds;
          }
        }
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

static void TestDeletedOrReplacedECnumbers (ValidStructPtr vsp)

{
  FileCache   fc;
  FILE        *fp = NULL;
  TextFsaPtr  fsa;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;

  /* only check first time program runs validator */

  fsa = (TextFsaPtr) GetAppProperty ("ReplacedEECNumberFSA");
  if (fsa != NULL) return;

  GetSpecificECNumberFSA ();
  GetAmbiguousECNumberFSA ();
  GetDeletedECNumberFSA ();
  GetReplacedECNumberFSA ();

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, "ecnum_replaced.txt");
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            if (! ECnumberNotInList (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replaced EC number %s still in live list", str);
            }
            if (ECnumberNotInList (ptr)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replacement EC number %s not in live list", ptr);
            }
            if (ECnumberWasDeleted (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replaced EC number %s in deleted list", str);
            }
            if (ECnumberWasDeleted (ptr)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Replacement EC number %s in deleted list", ptr);
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
    }
  }

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, "ecnum_deleted.txt");
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            if (! ECnumberNotInList (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Deleted EC number %s still in live list", str);
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
    }
  }
}

NLM_EXTERN Boolean ValidateSeqEntry (SeqEntryPtr sep, ValidStructPtr vsp)

{
  AuthListPtr     alp;
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  Uint2           entityID = 0;
  GatherScope     gs;
  BioseqSetPtr    bssp;
  SeqSubmitPtr    ssp = NULL;
  Boolean         do_many = FALSE;
  Boolean         mult_subs = FALSE;
  Boolean         farFetchProd;
  Boolean         first = TRUE;
  Int2            errors[6], i;
  Boolean         suppress_no_pubs = TRUE;
  Boolean         suppress_no_biosrc = TRUE;
  Boolean         other_sets_in_sep = FALSE;
  FeatProb        featprob;
  GatherContextPtr gcp = NULL;
  GatherContext   gc;
  SeqEntryPtr     fsep;
  BioseqPtr       fbsp = NULL;
  Int2            limit;
  SeqEntryPtr     oldsep;
  ErrSev          oldsev;
  ObjMgrDataPtr   omdp;
  SeqEntryPtr     topsep = NULL;
  SeqEntryPtr     tmp;
  ValNodePtr      bsplist;
  SubmitBlockPtr  sbp;
  ErrSev          sev;
  SeqIdPtr        sip;
  Boolean         isGPS = FALSE;
  Boolean         isPatent = FALSE;
  Boolean         isPDB = FALSE;
  FindRepData     frd;

  if (sep == NULL || vsp == NULL) return FALSE;

  genetic_code_name_list = SetUpValidateGeneticCodes ();

  vsp->useSeqMgrIndexes = TRUE; /* now always use indexing */

  for (i = 0; i < 6; i++)       /* keep errors between clears */
    errors[i] = 0;

  MemSet ((Pointer) &featprob, 0, sizeof (FeatProb));

  if (vsp->useSeqMgrIndexes) {
    entityID = ObjMgrGetEntityIDForChoice (sep);

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      oldsev = ErrSetMessageLevel (SEV_MAX);
      SeqMgrIndexFeatures (entityID, NULL);
      ErrSetMessageLevel (oldsev);
    }
    SeqMgrExploreBioseqs (entityID, NULL, (Pointer) &featprob, CountMisplacedFeatures, TRUE, TRUE, TRUE);

    topsep = GetTopSeqEntryForEntityID (entityID);
    VisitGraphsInSep (topsep, (Pointer) &featprob, CheckGraphPacking);
    VisitFeaturesInSep (topsep, (Pointer) &featprob, CountGeneXrefs);
    VisitFeaturesInSep (topsep, (Pointer) &featprob, CountFeatLocIdTypes);
    VisitBioseqsInSep (topsep, (Pointer) &featprob, CheckTpaHist);
  } else {

    /* if not using indexing, still need feature->idx.subtype now */

    entityID = ObjMgrGetEntityIDForChoice (sep);
    AssignIDsInEntity (entityID, 0, NULL);
  }

  /* Seq-submit can have multiple entries with no Bioseq-set wrapper */

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->data != NULL) {
      if (sep->next != NULL) {
        do_many = TRUE;
        mult_subs = TRUE;
      }
    }
  }

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) (sep->data.ptrvalue);
    switch (bssp->_class) {
    /* case BioseqseqSet_class_genbank: */
    case BioseqseqSet_class_pir:
    case BioseqseqSet_class_gibb:
    case BioseqseqSet_class_gi:
    case BioseqseqSet_class_swissprot:
      sep = bssp->seq_set;
      do_many = TRUE;
      break;
    case BioseqseqSet_class_gen_prod_set:
      isGPS = TRUE;
    default:
      break;
    }
  }

  /* if no pubs or biosource, only one message, not one per bioseq */

  if (mult_subs) {
    for (tmp = sep; tmp != NULL; tmp = tmp->next) {
      LookForAnyPubAndOrg (tmp, &suppress_no_pubs, &suppress_no_biosrc);
    }
  } else {
    LookForAnyPubAndOrg (sep, &suppress_no_pubs, &suppress_no_biosrc);
  }

  if (GetAppProperty ("ValidateExons") != NULL) {
    vsp->validateExons = TRUE;
  }

  vsp->is_htg_in_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_htg_in_sep), LookForHTG);
  vsp->is_barcode_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_barcode_sep), LookForBarcode);
  vsp->is_smupd_in_sep = FALSE;
  VisitDescriptorsInSep (sep, (Pointer) &(vsp->is_smupd_in_sep), LookForSMUPD);
  vsp->is_gps_in_sep = FALSE;
  SeqEntryExplore (sep, (Pointer) &(vsp->is_gps_in_sep), LookForGPS);
  other_sets_in_sep = FALSE;
  SeqEntryExplore (sep, (Pointer) &(other_sets_in_sep), LookForNonGPS);
  vsp->is_refseq_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_refseq_in_sep), LookForNC);
  vsp->is_embl_ddbj_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_embl_ddbj_in_sep), LookForEmblDdbj);
  vsp->is_insd_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->is_insd_in_sep), LookForGEDseqID);
  vsp->only_lcl_gnl_in_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->only_lcl_gnl_in_sep), LookForLclGnl);
  vsp->has_gnl_prot_sep = FALSE;
  VisitBioseqsInSep (sep, (Pointer) &(vsp->has_gnl_prot_sep), LookForProteinGnl);

  vsp->feat_loc_has_gi = featprob.loc_has_gi;
  vsp->feat_prod_has_gi = featprob.prod_has_gi;

  globalvsp = vsp;              /* for spell checker */

  while (sep != NULL) {
    vsp->far_fetch_failure = FALSE;

    /* calculate strings for LookForMultipleUnpubPubs test only once for genome product set efficiency */
    VisitDescriptorsInSep (sep, NULL, SetPubScratchData);

    MemSet (&gs, 0, sizeof (GatherScope));
    gs.scope = sep;             /* default is to scope to this set */

    ValidStructClear (vsp);
    vsp->sep = sep;

    MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
    gcp = &gc;
    gc.entityID = ObjMgrGetEntityIDForChoice (sep);
    gc.itemID = 1;
    if (IS_Bioseq (sep)) {
      gc.thistype = OBJ_BIOSEQ;
    } else {
      gc.thistype = OBJ_BIOSEQSET;
    }
    vsp->gcp = gcp;             /* above needed for ValidErr */
    vsp->suppress_no_pubs = suppress_no_pubs;
    vsp->suppress_no_biosrc = suppress_no_biosrc;

    if (vsp->is_refseq_in_sep && vsp->is_insd_in_sep) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_INSDRefSeqPackaging,
                "INSD and RefSeq records should not be present in the same set");
    }

    if (vsp->is_gps_in_sep && other_sets_in_sep) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_GPSnonGPSPackaging,
                "Genomic product set and mut/pop/phy/eco set records should not be present in the same set");
    }

    /* build seqmgr feature indices if not already done */

    bsplist = NULL;
    if (vsp->useSeqMgrIndexes) {
      entityID = ObjMgrGetEntityIDForChoice (sep);

      if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
        oldsev = ErrSetMessageLevel (SEV_MAX);
        SeqMgrIndexFeatures (entityID, NULL);
        ErrSetMessageLevel (oldsev);
      }

      /* lock all remote genome components, locations, and products in advance */

      limit = vsp->validationLimit;
      if (limit == VALIDATE_ALL || limit == VALIDATE_INST || limit == VALIDATE_HIST) {
        farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
        oldsev = ErrSetMessageLevel (SEV_WARNING);
        bsplist = LockFarComponentsEx (sep, TRUE, TRUE, farFetchProd, NULL);
        ErrSetMessageLevel (oldsev);
      }
    }

    fsep = FindNthBioseq (sep, 1);
    fbsp = NULL;
    if (fsep != NULL && IS_Bioseq (fsep)) {
      fbsp = (BioseqPtr) fsep->data.ptrvalue;
      /* report context as first bioseq */
      vsp->bsp = fbsp;
    }

    if (fbsp == NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NoBioseqFound, "No Bioseqs in this entire record.");
    } else {

      for (sip = fbsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_PATENT) {
          isPatent = TRUE;
        } else if (sip->choice == SEQID_PDB) {
          isPDB = TRUE;
        }
      }
  
      if (first) {
        TestDeletedOrReplacedECnumbers (vsp);
  
        if (suppress_no_pubs && (! vsp->seqSubmitParent)) {
          omdp = ObjMgrGetData (gc.entityID);
          if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
            sev = SEV_ERROR;
            if ((!isGPS) && (!IsNoncuratedRefSeq (fbsp, &sev)) && (! IsGpipe (fbsp))) {
              ValidErr (vsp, sev, ERR_SEQ_DESCR_NoPubFound, "No publications anywhere on this entire record.");
            }
          }
        }
        if (suppress_no_biosrc) {
          if ((!isPatent) && ((!isPDB))) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name anywhere on this entire record.");
          }
        }
  
        if (featprob.num_misplaced_features > 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_FeaturePackagingProblem, "There are %d mispackaged features in this record.", (int) featprob.num_misplaced_features);
        } else if (featprob.num_misplaced_features == 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_FeaturePackagingProblem, "There is %d mispackaged feature in this record.", (int) featprob.num_misplaced_features);
        }
  
        if (featprob.num_misplaced_graphs > 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GraphPackagingProblem, "There are %d mispackaged graphs in this record.", (int) featprob.num_misplaced_graphs);
        } else if (featprob.num_misplaced_graphs == 1) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GraphPackagingProblem, "There is %d mispackaged graph in this record.", (int) featprob.num_misplaced_graphs);
        }
  
        /*
        if (featprob.num_archaic_locations > 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureLocation, "There are %d archaic feature locations in this record.", (int) featprob.num_archaic_locations);
        } else if (featprob.num_archaic_locations == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureLocation, "There is %d archaic feature location in this record.", (int) featprob.num_archaic_locations);
        }
  
        if (featprob.num_archaic_products > 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureProduct, "There are %d archaic feature products in this record.", (int) featprob.num_archaic_products);
        } else if (featprob.num_archaic_products == 1) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_ArchaicFeatureProduct, "There is %d archaic feature product in this record.", (int) featprob.num_archaic_products);
        }
        */
  
        if (featprob.num_gene_feats == 0 && featprob.num_gene_xrefs > 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OnlyGeneXrefs, "There are %ld gene xrefs and no gene features in this record.", (long) featprob.num_gene_xrefs);
        }
  
        if (featprob.num_tpa_with_hist > 0 && featprob.num_tpa_without_hist > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_TpaAssmeblyProblem, "There are %ld TPAs with history and %ld without history in this record.",
                    (long) featprob.num_tpa_with_hist, (long) featprob.num_tpa_without_hist);
        }
  
        if (featprob.has_gi && featprob.num_tpa_without_hist > 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TpaAssmeblyProblem, "There are %ld TPAs without history in this record, but the record has a gi number assignment.",
                    (long) featprob.num_tpa_without_hist);
        }
  
        if (vsp->indexerVersion && vsp->has_gnl_prot_sep && (! vsp->is_refseq_in_sep)) {
          if (FindNucBioseq (sep) != NULL) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_ProteinsHaveGeneralID, "INDEXER_ONLY - Protein bioseqs have general seq-id.");
          }
        }
  
        first = FALSE;
      }
  
      vsp->bsp = NULL;
  
      topsep = GetTopSeqEntryForEntityID (gc.entityID);
      oldsep = SeqEntrySetScope (topsep);
  
      ValidateGeneLocusTags (topsep, vsp);
  
      VisitFeaturesInSep (sep, NULL, AddScratchToFeatures);
      VisitBioseqsInSep (sep, NULL, SetupFeatureScratchData);
  
      /* AssignIDsInEntity (gc.entityID, 0, NULL); */
  
      GatherSeqEntry (sep, (Pointer) vsp, Valid1GatherProc, &gs);
  
      if (ssp != NULL) {
        if (ssp->datatype == 1) {
          vsp->bsp = NULL;
          vsp->bssp = NULL;
          vsp->sfp = NULL;
          vsp->descr = NULL;
          vsp->gcp = NULL;
          sbp = ssp->sub;
          if (sbp != NULL) {
            csp = sbp->cit;
            if (csp != NULL) {
              alp = csp->authors;
              if (alp != NULL) {
                ValidateAffil (vsp, alp->affil);
              }
            }
            cip = sbp->contact;
            if (cip != NULL) {
              ap = cip->contact;
              if (ap != NULL) {
                ValidateAffil (vsp, ap->affil);
              }
            }
          }
        }
      }
  
      vsp->gcp = NULL;
      ValidateFeatCits (sep, vsp);
      vsp->gcp = NULL;
  
      vsp->gcp = NULL;
      ValidateFeatIDs (gc.entityID, vsp);
      vsp->gcp = NULL;
  
      vsp->gcp = NULL;
      ValidateSeqIdCase (sep, vsp);
      vsp->gcp = NULL;
  
      if (vsp->validateAlignments) {
        vsp->gcp = NULL;
        ValidateSeqAlignWithinValidator (vsp, sep, vsp->alignFindRemoteBsp, vsp->doSeqHistAssembly);
        vsp->gcp = NULL;
      }
  
      if (vsp->far_fetch_failure) {
        vsp->gcp = NULL;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_FarFetchFailure, "Far fetch failures caused some validator tests to be bypassed");
      }
  
      VisitFeaturesInSep (sep, NULL, ClearScratchOnFeatures);
  
      SeqEntrySetScope (oldsep);
  
      VisitDescriptorsInSep (sep, NULL, ClearPubScratchData);
    }

    if (vsp->useSeqMgrIndexes) {

      /* unlock all pre-locked remote genome components */

      bsplist = UnlockFarComponents (bsplist);
    }

    if (do_many) {
      for (i = 0; i < 6; i++)
        errors[i] += vsp->errors[i];
      sep = sep->next;
    } else
      sep = NULL;
  }

  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  gcp = &gc;
  gc.entityID = ObjMgrGetEntityIDForChoice (sep);
  vsp->gcp = gcp;
  frd.vsp = vsp;
  frd.gcp = gcp;

  limit = vsp->validationLimit;
  if (limit == VALIDATE_ALL) {
    /*
    frd.string = "?";
    */
    FindStringsInEntity (entityID, findrepstrs, FALSE, FALSE, FALSE, UPDATE_NEVER,
                         NULL, NULL, NULL, TRUE, FindRepValidate, (Pointer) &frd);
  }

  if (do_many) {
    for (i = 0; i < 6; i++)
      vsp->errors[i] = errors[i];
  }

  genetic_code_name_list = ValNodeFreeData (genetic_code_name_list);

  return TRUE;
}


static void ValidateSetContents (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr       bsp;
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) data;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) (sep->data.ptrvalue);
    if (ISA_aa (bsp->mol))
      vsp->protcnt++;
    else
      vsp->nuccnt++;
    if (bsp->repr == Seq_repr_seg)
      vsp->segcnt++;

  }
  return;
}


static CharPtr GetBioseqSetClass (Uint1 cl)
{
  if (cl == BioseqseqSet_class_nuc_prot)
    return ("nuc-prot");
  if (cl == BioseqseqSet_class_segset)
    return ("segset");
  if (cl == BioseqseqSet_class_conset)
    return ("conset");
  if (cl == BioseqseqSet_class_parts)
    return ("parts");
  if (cl == BioseqseqSet_class_gibb)
    return ("gibb");
  if (cl == BioseqseqSet_class_gi)
    return ("gi");
  if (cl == BioseqseqSet_class_genbank)
    return ("genbank");
  if (cl == BioseqseqSet_class_pir)
    return ("pir");
  if (cl == BioseqseqSet_class_pub_set)
    return ("pub-set");
  if (cl == BioseqseqSet_class_equiv)
    return ("equiv");
  if (cl == BioseqseqSet_class_swissprot)
    return ("swissprot");
  if (cl == BioseqseqSet_class_pdb_entry)
    return ("pdb-entry");
  if (cl == BioseqseqSet_class_mut_set)
    return ("mut-set");
  if (cl == BioseqseqSet_class_pop_set)
    return ("pop-set");
  if (cl == BioseqseqSet_class_phy_set)
    return ("phy-set");
  if (cl == BioseqseqSet_class_eco_set)
    return ("eco-set");
  if (cl == BioseqseqSet_class_gen_prod_set)
    return ("gen-prod-set");
  if (cl == BioseqseqSet_class_wgs_set)
    return ("wgs-set");
  if (cl == BioseqseqSet_class_other)
    return ("other");
  return ("not-set");
}


static BioseqSetPtr FindGenProdSetParentOfBioseqSet (BioseqSetPtr bssp)
{
  if (bssp == NULL) {
    return NULL;
  } else if (bssp->idx.parenttype != OBJ_BIOSEQSET) {
    return NULL;
  } else if ((bssp = (BioseqSetPtr)bssp->idx.parentptr) == NULL) {
    return NULL;
  } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return bssp;
  } else {
    return FindGenProdSetParentOfBioseqSet (bssp);
  }
}


static BioseqSetPtr FindGenProdSetParentOfBioseq (BioseqPtr bsp)
{
  BioseqSetPtr bssp;
  if (bsp == NULL) {
    return NULL;
  } else if (bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return NULL;
  } else if ((bssp = (BioseqSetPtr)bsp->idx.parentptr) == NULL) {
    return NULL;
  } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return bssp;
  } else {
    return FindGenProdSetParentOfBioseqSet (bssp);
  }
}


static void IfInGPSmustBeMrnaProduct (ValidStructPtr vsp, BioseqPtr bsp)

{
  /* see if in genomic product */
  if (FindGenProdSetParentOfBioseq(bsp) != NULL) {
    if (SeqMgrGetRNAgivenProduct (bsp, NULL) == NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Nucleotide bioseq should be product of mRNA feature on contig, but is not");
    }
  }
}

static void IfInGPSmustBeCDSProduct (ValidStructPtr vsp, BioseqPtr bsp)

{
  BioseqSetPtr  bssp;
  BioseqPtr     contig;
  ValNodePtr    head, vnp;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;

  /* see if in genomic product */
  if ((bssp = FindGenProdSetParentOfBioseq(bsp)) != NULL) {
    sep = bssp->seq_set;
    if (sep == NULL) return;
    if (! IS_Bioseq (sep)) return;
    contig = (BioseqPtr) sep->data.ptrvalue;
    if (contig == NULL) return;
    head = SeqMgrGetSfpProductList (bsp);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp == NULL) continue;
      if (BioseqFindFromSeqLoc (sfp->location) == contig) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Protein bioseq should be product of CDS feature on contig, but is not");
  }
}

NLM_EXTERN ValNodePtr BioseqGetSeqDescr(BioseqPtr bsp, Int2 type, ValNodePtr curr);


static void ValidateNucProtSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqDescrPtr   sdp;
  SeqEntryPtr   sep;
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp1;
  OrgRefPtr     orp;
  Int4          prot_biosource = 0;

  if (bssp->_class != BioseqseqSet_class_nuc_prot)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp == NULL) continue;
      if (ISA_na (bsp->mol)) {
        IfInGPSmustBeMrnaProduct (vsp, bsp);
      } else if (ISA_aa (bsp->mol)) {
        IfInGPSmustBeCDSProduct (vsp, bsp);
        sdp = BioseqGetSeqDescr (bsp, Seq_descr_source, NULL);
        if (sdp != NULL) {
          prot_biosource++;
        }
      }
    }

    if (!IS_Bioseq_set (sep))
      continue;

    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL)
      continue;

    if (bssp1->_class != BioseqseqSet_class_segset) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_NucProtNotSegSet,
                "Nuc-prot Bioseq-set contains wrong Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }

  if (prot_biosource > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceOnProtein,
              "Nuc-prot set has %ld proteins with a BioSource descriptor", (long) prot_biosource);
  } else if (prot_biosource > 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceOnProtein,
              "Nuc-prot set has %ld protein with a BioSource descriptor", (long) prot_biosource);
  }

  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL && StringDoesHaveText (orp->taxname)) return;
      }
    }
  }

  sep = vsp->sep;
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
      sep = bssp->seq_set;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
      }
    }
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return;
    }
  }

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceMissing,
            "Nuc-prot set does not contain expected BioSource descriptor");
}

typedef struct incons {
  Boolean     diffs;
  MolInfoPtr  mip;
} Incons, PNTR InconsPtr;

static void FindInconsistMolInfos (SeqDescrPtr sdp, Pointer userdata)

{
  InconsPtr   icp;
  MolInfoPtr  mip;

  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  icp = (InconsPtr) userdata;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (icp == NULL || mip == NULL) return;
  if (icp->mip == NULL) {
    icp->mip = mip;
  } else {
    if (icp->mip->biomol != mip->biomol) {
      icp->diffs = TRUE;
    }
  }
}

static void ValidateSegmentedSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp1;
  BioseqPtr       bsp;
  Incons          inc;
  Uint1           mol = 0;

  if (bssp->_class != BioseqseqSet_class_segset)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (mol == 0 || mol == Seq_mol_other) {
          mol = bsp->mol;
        } else if (bsp->mol != Seq_mol_other) {
          if (ISA_na (bsp->mol) != ISA_na (mol)) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_SegSetMixedBioseqs, "Segmented set contains mixture of nucleotides and proteins");
          }
        }
      }
    }

    if (!IS_Bioseq_set (sep))
      continue;

    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL)
      continue;

    if (bssp1->_class != BioseqseqSet_class_parts) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_SegSetNotParts,
                "Segmented set contains wrong Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }

  inc.diffs = FALSE;
  inc.mip = NULL;
  VisitDescriptorsInSet (bssp, (Pointer) &inc, FindInconsistMolInfos);
  if (inc.diffs) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_InconsistentMolInfoBiomols, "Segmented set contains inconsistent MolInfo biomols");
  }
}

static void ValidatePartsSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp1;
  BioseqPtr       bsp;
  Uint1           mol = 0;

  if (bssp->_class != BioseqseqSet_class_parts)
    return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (mol == 0 || mol == Seq_mol_other) {
          mol = bsp->mol;
        } else if (bsp->mol != Seq_mol_other) {
          if (ISA_na (bsp->mol) != ISA_na (mol)) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_PartsSetMixedBioseqs, "Parts set contains mixture of nucleotides and proteins");
            break;
          }
        }
      }
    }
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq_set (sep)) {
      bssp1 = sep->data.ptrvalue;
      if (bssp1 == NULL)
        continue;

      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_PartsSetHasSets,
                "Parts set contains unwanted Bioseq-set, its class is \"%s\".", GetBioseqSetClass (bssp1->_class));
      break;
    }
  }
}

static Boolean CheckForInconsistentBiosources (SeqEntryPtr sep, ValidStructPtr vsp, OrgRefPtr PNTR orpp, BioseqSetPtr top)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqEntryPtr     tmp;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioSourcePtr    biop;
  OrgRefPtr       orp;
  OrgRefPtr       firstorp;
  GatherContextPtr gcp;
  Uint2           entityID = 0, oldEntityID;
  Uint4           itemID = 0, oldItemID;
  Uint2           itemtype = 0, oldItemtype;
  size_t          len, len1, len2;
  ErrSev          sev;
  CharPtr         sp;

  if (sep == NULL || vsp == NULL || orpp == NULL)
    return FALSE;
  gcp = vsp->gcp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return FALSE;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      if (CheckForInconsistentBiosources (tmp, vsp, orpp, top))
        return TRUE;
    }
    return FALSE;
  }

  if (!IS_Bioseq (sep))
    return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return FALSE;

  biop = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    entityID = dcontext.entityID;
    itemID = dcontext.itemID;
    itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      entityID = fcontext.entityID;
      itemID = fcontext.itemID;
      itemtype = OBJ_SEQFEAT;
    }
  }
  if (biop == NULL)
    return FALSE;
  orp = biop->org;
  if (orp == NULL)
    return FALSE;

  firstorp = *orpp;
  if (firstorp == NULL) {
    *orpp = orp;
    return FALSE;
  }

  if (StringNICmp (orp->taxname, "Influenza virus ", 16) == 0 &&
      StringNICmp (firstorp->taxname, "Influenza virus ", 16) == 0 && StringNICmp (orp->taxname, firstorp->taxname, 17) == 0) {
    return FALSE;
  }

  if (StringICmp (orp->taxname, firstorp->taxname) == 0)
    return FALSE;

  sev = SEV_ERROR;
  sp = StringStr (orp->taxname, " sp. ");
  if (sp != NULL) {
    len = sp - orp->taxname + 5;
    if (StringNCmp (orp->taxname, firstorp->taxname, len) == 0) {
      sev = SEV_WARNING;
    }
  }

  if (sev == SEV_ERROR) {
    len1 = StringLen (orp->taxname);
    len2 = StringLen (firstorp->taxname);
    len = MIN (len1, len2);
    if (len > 0 && StringNCmp (orp->taxname, firstorp->taxname, len) == 0) {
      sev = SEV_WARNING;
    }
  }

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  if (top != NULL) {
    gcp->entityID = top->idx.entityID;
    gcp->itemID = top->idx.itemID;
    gcp->thistype = OBJ_BIOSEQSET;
  }

  /* only report the first one that doesn't match - but might be lower severity if not all are sp. */

  ValidErr (vsp, sev, ERR_SEQ_DESCR_InconsistentBioSources, "Population set contains inconsistent organisms.");

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  return TRUE;
}

static Boolean CheckForInconsistentMolInfos (SeqEntryPtr sep, ValidStructPtr vsp, MolInfoPtr PNTR mipp, BioseqSetPtr top)

{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqMgrDescContext  dcontext;
  Uint2              entityID = 0, oldEntityID;
  MolInfoPtr         firstmip;
  GatherContextPtr   gcp;
  Uint4              itemID = 0, oldItemID;
  Uint2              itemtype = 0, oldItemtype;
  MolInfoPtr         mip;
  ValNodePtr         sdp;
  SeqEntryPtr        tmp;

  if (sep == NULL || vsp == NULL || mipp == NULL)
    return FALSE;
  gcp = vsp->gcp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL)
      return FALSE;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      if (CheckForInconsistentMolInfos (tmp, vsp, mipp, top))
        return TRUE;
    }
    return FALSE;
  }

  if (!IS_Bioseq (sep))
    return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL || mip->biomol == MOLECULE_TYPE_PEPTIDE) return FALSE;

  firstmip = *mipp;
  if (firstmip == NULL) {
    *mipp = mip;
    return FALSE;
  }

  if (mip->biomol == firstmip->biomol) return FALSE;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  gcp->entityID = entityID;
  gcp->itemID = itemID;
  gcp->thistype = itemtype;

  if (top != NULL) {
    gcp->entityID = top->idx.entityID;
    gcp->itemID = top->idx.itemID;
    gcp->thistype = OBJ_BIOSEQSET;
  }

  /* only report the first one that doesn't match */

  ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InconsistentMolInfoBiomols, "Pop/phy/mut/eco set contains inconsistent MolInfo biomols");

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  return TRUE;
}

static void LookForMolInfoInconsistency (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  MolInfoPtr    mip = NULL;
  SeqEntryPtr   sep;

  if (bssp == NULL) return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (CheckForInconsistentMolInfos (sep, vsp, &mip, bssp))
      return;
  }
}

static void ValidatePopSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr  bssp1;
  OrgRefPtr     orp = NULL;
  SeqEntryPtr   sep;

  if (bssp->_class != BioseqseqSet_class_pop_set)
    return;

  if (vsp->is_refseq_in_sep) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_RefSeqPopSet,
              "RefSeq record should not be a Pop-set");
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  LookForMolInfoInconsistency (bssp, vsp);

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (CheckForInconsistentBiosources (sep, vsp, &orp, bssp))
      return;
  }
}

static void ValidateGenbankSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr    bssp1;
  SeqEntryPtr     sep;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }
}

static void ValidatePhyMutEcoWgsSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqSetPtr  bssp1;
  SeqEntryPtr   sep;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (!IS_Bioseq_set (sep)) continue;
    bssp1 = sep->data.ptrvalue;
    if (bssp1 == NULL) continue;

    if (bssp1->_class == BioseqseqSet_class_genbank) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_InternalGenBankSet,
                "Bioseq-set contains internal GenBank Bioseq-set");
    }
  }

  LookForMolInfoInconsistency (bssp, vsp);
}

static void ValidateGenProdSet (BioseqSetPtr bssp, ValidStructPtr vsp)

{
  BioseqPtr       bsp;
  BioseqPtr       cdna;
  SeqMgrFeatContext fcontext;
  GatherContextPtr gcp = NULL;
  CharPtr         loc = NULL;
  SeqFeatPtr      mrna;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;

  if (bssp->_class != BioseqseqSet_class_gen_prod_set)
    return;

  if (bssp->annot != NULL) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Seq-annot packaged directly on genomic product set");
  }

  sep = bssp->seq_set;
  if (!IS_Bioseq (sep))
    return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL)
    return;

  gcp = vsp->gcp;
  if (gcp == NULL)
    return;
  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;

  if (vsp->useSeqMgrIndexes) {
    mrna = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
    while (mrna != NULL) {
      cdna = BioseqFindFromSeqLoc (mrna->product);
      if (cdna == NULL) {
        gcp->itemID = mrna->idx.itemID;
        gcp->thistype = OBJ_SEQFEAT;
        loc = SeqLocPrint (mrna->product);
        if (loc == NULL) {
          loc = StringSave ("?");
        }
        sip = SeqLocId (mrna->product);
        /* okay to have far RefSeq product */
        if (sip == NULL || sip->choice != SEQID_OTHER) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Product of mRNA feature (%s) not packaged in genomic product set", loc);
        }
        MemFree (loc);
      }
      mrna = SeqMgrGetNextFeature (bsp, mrna, 0, FEATDEF_mRNA, &fcontext);
    }
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
}

static void ValidateBioseqSet (GatherContextPtr gcp)

{
  BioseqSetPtr    bssp;
  ValidStructPtr  vsp;
  SeqEntryPtr     sep;

  vsp = (ValidStructPtr) (gcp->userdata);
  bssp = (BioseqSetPtr) (gcp->thisitem);
  vsp->bssp = bssp;
  vsp->bsp = NULL;
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if (vsp->non_ascii_chars) {   /* non_ascii chars in AsnRead step */
    ValidErr (vsp, SEV_ERROR, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
    vsp->non_ascii_chars = FALSE;       /* only do once */
  }

  vsp->nuccnt = 0;
  vsp->segcnt = 0;
  vsp->protcnt = 0;

  sep = gcp->sep;

  SeqEntryExplore (sep, (Pointer) vsp, ValidateSetContents);

  switch (bssp->_class) {
  case BioseqseqSet_class_nuc_prot:
    if (vsp->nuccnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No nucleotides in nuc-prot set");
    }
    if (vsp->protcnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No proteins in nuc-prot set");
    }
    if (vsp->nuccnt > 1 && vsp->segcnt == 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_NucProtProblem, "Multiple unsegmented nucleotides in nuc-prot set");
    }
    ValidateNucProtSet (bssp, vsp);
    break;
  case BioseqseqSet_class_segset:
    if (vsp->segcnt == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_SegSetProblem, "No segmented Bioseq in segset");
    }
    ValidateSegmentedSet (bssp, vsp);
    break;
  case BioseqseqSet_class_conset:
    if (vsp->indexerVersion && (! vsp->is_refseq_in_sep)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_ConSetProblem, "Set class should not be conset");
    }
    break;
  case BioseqseqSet_class_parts:
    ValidatePartsSet (bssp, vsp);
    break;
  case BioseqseqSet_class_genbank:
    ValidateGenbankSet (bssp, vsp);
    break;
  case BioseqseqSet_class_pop_set:
    ValidatePopSet (bssp, vsp);
    break;
  case BioseqseqSet_class_mut_set:
  case BioseqseqSet_class_phy_set:
  case BioseqseqSet_class_eco_set:
  case BioseqseqSet_class_wgs_set:
    ValidatePhyMutEcoWgsSet (bssp, vsp);
    break;
  case BioseqseqSet_class_gen_prod_set:
    ValidateGenProdSet (bssp, vsp);
    break;
  /*
  case BioseqseqSet_class_other:
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_PKG_GenomicProductPackagingProblem, "Genomic product set class incorrectly set to other");
    break;
  */
  default:
    if (!((vsp->nuccnt) || (vsp->protcnt))) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_PKG_EmptySet, "No Bioseqs in this set");
    }
    break;
  }
  return;
}

static Boolean SuppressTrailingXMessage (BioseqPtr bsp)
{
  ByteStorePtr    bs;
  SeqFeatPtr      cds;
  Boolean         hasstar;
  Int4            len;
  MolInfoPtr      mip;
  SeqDescrPtr     sdp;
  CharPtr         str;

  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds != NULL) {
    bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
    if (bs != NULL) {
      str = BSMerge (bs, NULL);
      BSFree (bs);
      hasstar = FALSE;
      if (str != NULL) {
        len = StringLen (str);
        if (len > 1 && str[len - 1] == '*') {
          hasstar = TRUE;
        }
      }
      MemFree (str);
      return hasstar;
    }
  }
  sdp = BioseqGetSeqDescr (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness == 4 || mip->completeness == 5)
        return TRUE;
    }
  }
  return FALSE;
}

static void LookForSecondaryConflict (ValidStructPtr vsp, GatherContextPtr gcp, CharPtr accn, ValNodePtr extra_acc)
{
  CharPtr         str;
  ValNodePtr      vnp;

  if (vsp == NULL || gcp == NULL)
    return;
  if (StringHasNoText (accn))
    return;
  for (vnp = extra_acc; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str))
      continue;
    if (StringICmp (accn, str) == 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSecondaryAccn, "%s used for both primary and secondary accession", accn);
    }
  }
}

static void CheckSegBspAgainstParts (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)
{
  BioseqSetPtr    bssp;
  BioseqPtr       part;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  SeqLocPtr       slp;

  if (vsp == NULL || gcp == NULL || bsp == NULL)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;

  if (bsp->repr != Seq_repr_seg || bsp->seq_ext_type != 1 || bsp->seq_ext == NULL)
    return;

  sep = bsp->seqentry;
  if (sep == NULL)
    return;
  sep = sep->next;
  if (sep == NULL)
    return;
  if (!IS_Bioseq_set (sep))
    return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL)
    return;
  if (bssp->_class != BioseqseqSet_class_parts)
    return;

  sep = bssp->seq_set;
  for (slp = (ValNodePtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
    if (slp->choice == SEQLOC_NULL)
      continue;
    if (sep == NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set does not contain enough Bioseqs");
      return;
    }
    if (IS_Bioseq (sep)) {
      part = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqLocId (slp);
      if (sip != NULL && part != NULL) {
        if (!SeqIdIn (sip, part->id)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Segmented bioseq seq_ext does not correspond to parts packaging order");
          return;
        }
      }
    } else {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set component is not Bioseq");
      return;
    }
    sep = sep->next;
  }
  if (sep != NULL) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartsOutOfOrder, "Parts set contains too many Bioseqs");
  }
}

/*****************************************************************************
*
*   ValidateBioseqHist(gcp)
*      Validate one Bioseq Seq-hist
*
*****************************************************************************/
static void ValidateBioseqHist (GatherContextPtr gcp)

{
  BioseqPtr       bsp;
  Int4            gi = 0;
  SeqHistPtr      hist;
  SeqIdPtr        sip;
  ValidStructPtr  vsp;

  if (gcp == NULL) return;
  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);
  vsp->bsp_partial_val = 0;

  if (bsp == NULL) return;
  hist = bsp->hist;
  if (hist == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
    }
  }
  if (gi == 0) return;

  if (hist->replaced_by_ids != NULL && hist->replaced_by_date != NULL) {

    for (sip = hist->replaced_by_ids; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        if (gi == (Int4) sip->data.intvalue) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_HistoryGiCollision, "Replaced by gi (%ld) is same as current Bioseq", (long) gi);
        }
      }
    }
  }

  if (hist->replace_ids != NULL && hist->replace_date != NULL) {

    for (sip = hist->replace_ids; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        if (gi == (Int4) sip->data.intvalue) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_HistoryGiCollision, "Replaces gi (%ld) is same as current Bioseq", (long) gi);
        }
      }
    }
  }
}

/*****************************************************************************
*
*   ValidateBioseqInst(gcp)
*      Validate one Bioseq Seq-inst
*
*****************************************************************************/
static Boolean IsTpa (
  BioseqPtr bsp,
  Boolean has_tpa_assembly,
  BoolPtr isRefSeqP
)

{
  DbtagPtr  dbt;
  Boolean   has_bankit = FALSE;
  Boolean   has_genbank = FALSE;
  Boolean   has_gi = FALSE;
  Boolean   has_local = FALSE;
  Boolean   has_refseq = FALSE;
  Boolean   has_smart = FALSE;
  Boolean   has_tpa = FALSE;
  SeqIdPtr  sip;

  if (bsp == NULL || bsp->id == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        has_local = TRUE;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        has_genbank = TRUE;
        break;
      case SEQID_OTHER :
        has_refseq = TRUE;
        if (isRefSeqP != NULL) {
          *isRefSeqP = TRUE;
        }
        break;
      case SEQID_GI :
        has_gi = TRUE;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        has_tpa = TRUE;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "BankIt") == 0) {
            has_bankit = TRUE;
          }
          if (StringICmp (dbt->db, "TMSMART") == 0) {
            has_smart = TRUE;
          }
        }
        break;
      case SEQID_GPIPE :
        break;
      default :
        break;
    }
  }

  if (has_genbank) return FALSE;
  if (has_tpa) return TRUE;
  if (has_refseq) return FALSE;
  if (has_bankit && has_tpa_assembly) return TRUE;
  if (has_smart && has_tpa_assembly) return TRUE;
  if (has_gi) return FALSE;
  if (has_local && has_tpa_assembly) return TRUE;

  return FALSE;
}

static void ValidateIDSetAgainstDb (GatherContextPtr gcp, ValidStructPtr vsp, BioseqPtr bsp)

{
  SeqIdPtr        sip, sipset;
  SeqIdPtr        gbId = NULL;
  SeqIdPtr        dbGbId;
  DbtagPtr        generalID = NULL;
  DbtagPtr        dbGeneralID;
  Int4            gi = 0;
  Int4            dbGI;
  Char            oldGenID [128], newGenID [128];

  if (gcp != NULL && vsp != NULL && bsp != NULL && vsp->validateIDSet) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GENBANK:
          gbId = sip;
          break;
        case SEQID_GI :
          gi = (Int4) sip->data.intvalue;
          break;
        case SEQID_GENERAL :
          generalID = (DbtagPtr) sip->data.ptrvalue;
          break;
        default :
          break;
      }
    }
    if (gi == 0 && gbId != NULL) {
      gi = GetGIForSeqId (gbId);
    }
    if (gi > 0) {
      sipset = GetSeqIdSetForGI (gi);
      if (sipset != NULL) {
        dbGI = 0;
        dbGbId = NULL;
        dbGeneralID = NULL;
        oldGenID [0] = '\0';
        newGenID [0] = '\0';
        for (sip = sipset; sip != NULL; sip = sip->next) {
          switch (sip->choice) {
            case SEQID_GI :
              dbGI = (Int4) sip->data.intvalue;
              break;
            case SEQID_GENBANK:
              dbGbId = sip;
              break;
            case SEQID_GENERAL :
              dbGeneralID = (DbtagPtr) sip->data.ptrvalue;
              break;
            default :
              break;
          }
        }
        if (dbGI != gi) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_UnexpectedIdentifierChange, "New gi number (%ld) does not match one in NCBI sequence repository (%ld)", (long) gi, (long) dbGI);
        }
        if (gbId != NULL && dbGbId != NULL) {
          if (! SeqIdMatch (gbId, dbGbId)) {
            SeqIdWrite (dbGbId, oldGenID, PRINTID_FASTA_SHORT, sizeof (oldGenID));
            SeqIdWrite (gbId, newGenID, PRINTID_FASTA_SHORT, sizeof (newGenID));
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "New accession (%s) does not match one in NCBI sequence repository (%s) on gi (%ld)", newGenID, oldGenID, (long) gi);
          }
        } else if (gbId != NULL) {
          SeqIdWrite (gbId, newGenID, PRINTID_FASTA_SHORT, sizeof (newGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Gain of accession (%s) on gi (%ld) compared to the NCBI sequence repository", newGenID, (long) gi);
        } else if (dbGbId != NULL) {
          SeqIdWrite (dbGbId, oldGenID, PRINTID_FASTA_SHORT, sizeof (oldGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Loss of accession (%s) on gi (%ld) compared to the NCBI sequence repository", oldGenID, (long) gi);
        }
        if (generalID != NULL && dbGeneralID != NULL) {
          if (! DbtagMatch (generalID, dbGeneralID)) {
            DbtagLabel (dbGeneralID, oldGenID, sizeof (oldGenID));
            DbtagLabel (generalID, newGenID, sizeof (newGenID));
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "New general ID (%s) does not match one in NCBI sequence repository (%s) on gi (%ld)", newGenID, oldGenID, (long) gi);
          }
        } else if (generalID != NULL) {
          DbtagLabel (generalID, newGenID, sizeof (newGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Gain of general ID (%s) on gi (%ld) compared to the NCBI sequence repository", newGenID, (long) gi);
        } else if (dbGeneralID != NULL) {
          DbtagLabel (dbGeneralID, oldGenID, sizeof (oldGenID));
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_UnexpectedIdentifierChange, "Loss of general ID (%s) on gi (%ld) compared to the NCBI sequence repository", oldGenID, (long) gi);
        }
      }
      SeqIdSetFree (sipset);
    }
  }
}

typedef struct enrun {
  GatherContextPtr  gcp;
  ValidStructPtr    vsp;
  Int4              ncount;
  Int4              maxrun;
  Int4              seqpos;
  Int4              gapcount;
  Boolean           showAll;
  Boolean           inNrun;
  Boolean           isWGS;
} RunOfNs, PNTR RunOfNsPtr;

static void LIBCALLBACK CountAdjacentProc (CharPtr sequence, Pointer userdata)

{
  Char              ch;
  GatherContextPtr  gcp;
  RunOfNsPtr        ronp;
  CharPtr           str;
  ValidStructPtr    vsp;

  ronp = (RunOfNsPtr) userdata;
  if (sequence == NULL || ronp == NULL) return;

  str = sequence;
  ch = *str;
  while (ch != '\0') {
    (ronp->seqpos)++;
    if (ch == 'N') {
      (ronp->ncount)++;
      if (ronp->ncount > ronp->maxrun) {
        ronp->maxrun = ronp->ncount;
      }
      ronp->inNrun = TRUE;
    } else {
      if (ch == '-') {
        (ronp->gapcount)++;
      }
      if (ronp->inNrun && ronp->showAll && ronp->isWGS && ronp->ncount >= 20 && ronp->seqpos > ronp->ncount + 1) {
        vsp = ronp->vsp;
        gcp = ronp->gcp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                  (long) ronp->ncount, (long) (ronp->seqpos - ronp->ncount));
      } else if (ronp->inNrun && ronp->showAll && ronp->ncount >= 100 && ronp->seqpos > ronp->ncount + 1) {
        vsp = ronp->vsp;
        gcp = ronp->gcp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                  (long) ronp->ncount, (long) (ronp->seqpos - ronp->ncount));
      }
      ronp->ncount = 0;
      ronp->inNrun = FALSE;
    }
    str++;
    ch = *str;
  }
}

static Int4 CountAdjacentNsInSeqLit (GatherContextPtr gcp, SeqLitPtr slitp, Boolean is_na)

{
  BioseqPtr  bsp;
  RunOfNs    ron;

  if (slitp == NULL || slitp->length < 1 || slitp->seq_data == NULL) return 0;

  bsp = BioseqNew ();
  if (bsp == NULL) return 0;

  if (slitp->seq_data != NULL) {
    bsp->repr = Seq_repr_raw;
  } else {
    bsp->repr = Seq_repr_virtual;
  }
  if (is_na) {
    bsp->mol = Seq_mol_dna;
  } else {
    bsp->mol = Seq_mol_aa;
  }
  bsp->seq_data_type = slitp->seq_data_type;
  bsp->seq_data = slitp->seq_data;
  bsp->length = slitp->length;
  bsp->id = SeqIdParse ("lcl|countseqlitns");

  ron.gcp = gcp;
  ron.vsp = (ValidStructPtr) (gcp->userdata);
  ron.ncount = 0;
  ron.maxrun = 0;
  ron.seqpos = 0;
  ron.gapcount = 0;
  ron.showAll = FALSE;
  ron.inNrun = FALSE;
  ron.isWGS = FALSE;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) &ron, CountAdjacentProc);

  bsp->seq_data = NULL;

  BioseqFree (bsp);

  return ron.maxrun;
}

static Boolean HasUnparsedBrackets (CharPtr title)

{
  CharPtr  str;

  if (StringHasNoText (title)) return FALSE;

  str = StringChr (title, '[');
  if (str == NULL) return FALSE;
  str = StringChr (str, '=');
  if (str == NULL) return FALSE;
  str = StringChr (str, ']');
  if (str == NULL) return FALSE;
  return TRUE;
}

static CharPtr GetSequencePlusGapByFeature (SeqFeatPtr sfp)

{
  Int4     len;
  CharPtr  str = NULL;

  if (sfp == NULL) return NULL;
  len = SeqLocLen (sfp->location);
  if (len > 0 && len < MAXALLOC) {
    str = MemNew (sizeof (Char) * (len + 2));
    if (str != NULL) {
      SeqPortStreamLoc (sfp->location, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
    }
  }
    
  return str;
}

static Boolean IsWgsIntermediate (SeqEntryPtr sep)

{
  BioseqPtr    bsp;
  Boolean      has_gi = FALSE, is_other = FALSE, is_wgs = FALSE;
  MolInfoPtr   mip;
  SeqDescrPtr  sdp;
  SeqIdPtr     sip;

  bsp = FindNucBioseq (sep);
  if (bsp == NULL) return FALSE;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_molinfo) continue;
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL) continue;
    if (mip->tech == MI_TECH_wgs) {
      is_wgs = TRUE;
    }
  }
  if (! is_wgs) return FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      is_other = TRUE;
    } else if (sip->choice == SEQID_GI) {
      has_gi = TRUE;
    }
  }
  if (! is_other) return FALSE;
  if (has_gi) return FALSE;

  return TRUE;
}

typedef struct reusedata {
  CharPtr  seqidstr;
  Int4     from;
  Int4     to;
} ReuseData, PNTR ReuseDataPtr;

static int LIBCALLBACK SortVnpByDeltaLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  ReuseDataPtr  rdp1, rdp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  rdp1 = (ReuseDataPtr) vnp1->data.ptrvalue;
  rdp2 = (ReuseDataPtr) vnp2->data.ptrvalue;
  if (rdp1 == NULL || rdp2 == NULL) return 0;

  compare = StringICmp (rdp1->seqidstr, rdp2->seqidstr);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (rdp1->from > rdp2->from) {
    return 1;
  } else if (rdp1->from < rdp2->from) {
    return -1;
  }

  if (rdp1->to > rdp2->to) {
    return 1;
  } else if (rdp1->to < rdp2->to) {
    return -1;
  }

  return 0;
}

static void CheckDeltaForReuse (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)

{
  Char          buf [80];
  ValNodePtr    head = NULL;
  ValNodePtr    last = NULL;
  ReuseDataPtr  lastrdp = NULL;
  ReuseDataPtr  rdp;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  ValNodePtr    vnp_dsp, vnp_r;

  if (vsp == NULL || gcp == NULL || bsp == NULL) return;

  for (vnp_dsp = (ValNodePtr) bsp->seq_ext; vnp_dsp != NULL; vnp_dsp = vnp_dsp->next) {
    if (vnp_dsp->choice != 1) continue;
    slp = (SeqLocPtr) vnp_dsp->data.ptrvalue;
    if (slp == NULL) continue;
    if (slp->choice != SEQLOC_INT) continue;
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp == NULL) continue;
    sip = sintp->id;
    if (sip == NULL) continue;
    if (! SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1)) continue;
    rdp = (ReuseDataPtr) MemNew (sizeof (ReuseData));
    if (rdp == NULL) continue;
    rdp->seqidstr = StringSave (buf);
    rdp->from = sintp->from;
    rdp->to = sintp->to;
    vnp_r = ValNodeAddPointer (&last, 0, (Pointer) rdp);
    if (head == NULL) {
      head = vnp_r;
    }
    last = vnp_r;
  }

  if (head == NULL) return;

  head = ValNodeSort (head, SortVnpByDeltaLoc);

  for (vnp_r = head; vnp_r != NULL; vnp_r = vnp_r->next) {
    rdp = (ReuseDataPtr) vnp_r->data.ptrvalue;
    if (rdp == NULL) continue;
    if (lastrdp != NULL) {
      if (StringICmp (lastrdp->seqidstr, rdp->seqidstr) == 0) {
        if (lastrdp->to >= rdp->from && lastrdp->from <= rdp->to) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_OverlappingDeltaRange,
                    "Overlapping delta range %ld-%ld and %ld-%ld on a Bioseq %s",
                    (long) rdp->from + 1, (long) rdp->to + 1, (long) lastrdp->from + 1,
                    (long) lastrdp->to + 1, rdp->seqidstr);
        }
      }
    }
    lastrdp = rdp;
  }

  for (vnp_r = head; vnp_r != NULL; vnp_r = vnp_r->next) {
    rdp = (ReuseDataPtr) vnp_r->data.ptrvalue;
    if (rdp == NULL) continue;
    rdp->seqidstr = MemFree (rdp->seqidstr);
  }
  ValNodeFreeData (head);
}

static CharPtr legal_refgene_status_strings [] = {
  "Inferred",
  "Provisional",
  "Predicted",
  "Validated",
  "Reviewed",
  "Model",
  "WGS",
  "Pipeline",
  NULL
};


static void ValidateBioseqInst (GatherContextPtr gcp)
{
  Boolean         retval = TRUE;
  Int2            i, start_at, num;
  Boolean         errors[4], check_alphabet;
  static char    *repr[8] = {
    "virtual", "raw", "segmented", "constructed",
    "reference", "consensus", "map", "delta"
  };
  /*
  SeqPortPtr      spp;
  */
  Int2            residue, x, termination, gapchar;
  Boolean         gapatstart;
  Int4            len, divisor = 1, len2, len3;
  ValNode         head, vn;
  ValNodePtr      vnp, idlist;
  BioseqContextPtr bcp;
  Boolean         got_partial, is_invalid;
  int             seqtype, terminations, dashes;
  ValidStructPtr  vsp;
  BioseqPtr       bsp, bsp2;
  SeqIdPtr        sip1, sip2, sip3;
  SeqLocPtr       slp;
  SeqIntPtr       sintp;
  Char            buf1[41], buf2[41];
  SeqLitPtr       slitp;
  SeqCodeTablePtr sctp;
  MolInfoPtr      mip;
  SeqMgrDescContext context;
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  GBBlockPtr      gbp;
  GeneRefPtr      grp;
  SeqFeatPtr      gene;
  SeqMgrFeatContext genectxt;
  CharPtr         genelbl = NULL;
  SeqFeatPtr      prot;
  SeqMgrFeatContext protctxt;
  CharPtr         protlbl = NULL;
  TextSeqIdPtr    tsip;
  CharPtr         ptr, last, str, title, buf;
  Uint1           lastchoice;
  Char            ch;
  Boolean         multitoken;
  Boolean         hasGi = FALSE;
  SeqHistPtr      hist;
  Boolean         hist_asm_missing = FALSE;
  IntFuzzPtr      ifp;
  Int4            adjacent_N_gap_position;
  Boolean         adjacent_N_and_gap;
  Boolean         in_gap;
  Boolean         in_N;
  Boolean         isActiveFin = FALSE;
  Boolean         isDraft = FALSE;
  Boolean         isFullTop = FALSE;
  Boolean         isGB = FALSE;
  Boolean         isPatent = FALSE;
  Boolean         isPDB = FALSE;
  Boolean         isPreFin = FALSE;
  Boolean         isNC = FALSE;
  Boolean         isNTorNC = FALSE;
  Boolean         isNZ;
  Boolean         is_gps = FALSE;
  Boolean         isRefSeq = FALSE;
  Boolean         isSwissProt = FALSE;
  ValNodePtr      keywords;
  Boolean         last_is_gap;
  Boolean         non_interspersed_gaps;
  Int2            num_adjacent_gaps;
  Int2            num_gaps;
  Boolean         reportFastaBracket;
  SeqFeatPtr      sfp;
  SeqEntryPtr     sep;
  ErrSev          sev;
  DbtagPtr        dbt;
  SeqIdPtr        sip;
  Int2            trailingX = 0;
  Int2            numletters, numdigits, numunderscores;
  Boolean         letterAfterDigit, badIDchars;
  EMBLBlockPtr    ebp;
  SeqDescrPtr     sdp;
  SeqMgrDescContext dcontext;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;
  size_t          buflen = 1001;
  ItemInfo        ii;
  Uint1           tech;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  ObjValNodePtr   ovp;
  BioseqSetPtr    bssp;
  UserObjectPtr   uop;
  UserFieldPtr    ufp;
  ObjectIdPtr     oip;
  Boolean         hasRefGeneTracking = FALSE;
  Boolean         hasRefTrackStatus;
  Boolean         hasLegalStatus;
  Int2            accn_count = 0;
  Int2            gi_count = 0;
  Int4            runsofn;
  Int4            segnum;
  StreamCache     sc;
  RunOfNs         ron;
  Boolean         leadingX;
  Boolean         isLower;
  Boolean         isFirst;
  CharPtr         bases;
  Int4            dnalen;
  Int4            total;

  /* set up data structures */

  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);
  vsp->bsp_partial_val = 0;

  sep = vsp->sep;

  if (vsp->non_ascii_chars) {   /* non_ascii chars in AsnRead step */
    ValidErr (vsp, SEV_REJECT, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
    vsp->non_ascii_chars = FALSE;       /* only do once */
  }

  if (bsp->id == NULL) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_NoIdOnBioseq, "No ids on a Bioseq");
    return;
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    if (sip1->choice == SEQID_OTHER) {
      isRefSeq = TRUE;
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
          isNTorNC = TRUE;
        } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
          isNTorNC = TRUE;
          isNC = TRUE;
        }
      }
    } else if (sip1->choice == SEQID_GI) {
      hasGi = TRUE;
    } else if (sip1->choice == SEQID_GENBANK) {
      isGB = TRUE;
    } else if (sip1->choice == SEQID_SWISSPROT) {
      isSwissProt = TRUE;
    }

    for (sip2 = sip1->next; sip2 != NULL; sip2 = sip2->next) {
      if (SeqIdComp (sip1, sip2) != SIC_DIFF) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        SeqIdWrite (sip2, buf2, PRINTID_FASTA_SHORT, 40);
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingIdsOnBioseq, "Conflicting ids on a Bioseq: (%s - %s)", buf1, buf2);
      }
    }
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    switch (sip1->choice) {
    case SEQID_TPG:
    case SEQID_TPE:
    case SEQID_TPD:
      hist = bsp->hist;
      if (hist == NULL || hist->assembly == NULL) {
        if (ISA_na (bsp->mol) && bsp->repr != Seq_repr_seg) {
          hist_asm_missing = TRUE;
          keywords = NULL;
          vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_genbank, NULL);
          if (vnp != NULL && vnp->choice == Seq_descr_genbank) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              keywords = gbp->keywords;
            }
          }
          if (keywords == NULL) {
            vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_embl, NULL);
            if (vnp != NULL && vnp->choice == Seq_descr_embl) {
              ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
              if (ebp != NULL) {
                keywords = ebp->keywords;
              }
            }
          }
          if (keywords != NULL) {
            for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
              str = (CharPtr) vnp->data.ptrvalue;
              if (StringHasNoText (str)) continue;
              if (StringICmp (str, "TPA:reassembly") == 0) {
                hist_asm_missing = FALSE;
              }
            }
          }
          if (hist_asm_missing) {
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, 40);
            ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_HistAssemblyMissing, "TPA record %s should have Seq-hist.assembly for PRIMARY block", buf1);
          }
        }
      }
      /* continue falling through */
    case SEQID_GENBANK:
    case SEQID_EMBL:
    case SEQID_DDBJ:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        numletters = 0;
        numdigits = 0;
        letterAfterDigit = FALSE;
        badIDchars = FALSE;
        for (ptr = tsip->accession, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_UPPER (ch)) {
            numletters++;
            if (numdigits > 0) {
              letterAfterDigit = TRUE;
            }
          } else if (IS_DIGIT (ch)) {
            numdigits++;
          } else {
            badIDchars = TRUE;
          }
        }
        if (letterAfterDigit || badIDchars) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        } else if (numletters == 1 && numdigits == 5 && ISA_na (bsp->mol)) {
        } else if (numletters == 2 && numdigits == 6 && ISA_na (bsp->mol)) {
        } else if (numletters == 3 && numdigits == 5 && ISA_aa (bsp->mol)) {
        } else if (numletters == 2 && numdigits == 6 && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_seg) {
        } else if (numletters == 4 && numdigits == 8 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ)) {
        } else if (numletters == 4 && numdigits == 9 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ)) {
        } else if (numletters == 5 && numdigits == 7 && ISA_na (bsp->mol) &&
                   (sip1->choice == SEQID_GENBANK || sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ)) {
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        }
        if (vsp->useSeqMgrIndexes) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
          if (vnp != NULL) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              LookForSecondaryConflict (vsp, gcp, tsip->accession, gbp->extra_accessions);
            }
          }
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_embl, &context);
          if (vnp != NULL) {
            ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
            if (ebp != NULL) {
              LookForSecondaryConflict (vsp, gcp, tsip->accession, ebp->extra_acc);
            }
          }
        }
        if (hasGi) {
          if (tsip->version == 0) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Accession %s has 0 version", tsip->accession);
          }
        }
      }
      /* and keep going with further test */
    case SEQID_OTHER:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && tsip->name != NULL) {
        multitoken = FALSE;
        for (ptr = tsip->name, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_WHITESP (ch)) {
            multitoken = TRUE;
          }
        }
        if (multitoken) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqIdNameHasSpace, "Seq-id.name '%s' should be a single word without any spaces", tsip->name);
        }
      }
      if (tsip != NULL && tsip->accession != NULL && sip1->choice == SEQID_OTHER) {
        numletters = 0;
        numdigits = 0;
        numunderscores = 0;
        letterAfterDigit = FALSE;
        badIDchars = FALSE;
        ptr = tsip->accession;
        isNZ = (Boolean) (StringNCmp (ptr, "NZ_", 3) == 0);
        if (isNZ) {
          ptr += 3;
        }
        for (ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_UPPER (ch)) {
            numletters++;
            if (numdigits > 0 || numunderscores > 0) {
              letterAfterDigit = TRUE;
            }
          } else if (IS_DIGIT (ch)) {
            numdigits++;
          } else if (ch == '_') {
            numunderscores++;
            if (numdigits > 0 || numunderscores > 1) {
              letterAfterDigit = TRUE;
            }
          } else {
            badIDchars = TRUE;
          }
        }
        if (letterAfterDigit || badIDchars) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        } else if (isNZ && numletters == 4 && numdigits == 8 && numunderscores == 0) {
        } else if (isNZ && ValidateAccn (tsip->accession) == 0) {
        } else if (numletters == 2 && numdigits == 6 && numunderscores == 1) {
        } else if (numletters == 2 && numdigits == 8 && numunderscores == 1) {
        } else if (numletters == 2 && numdigits == 9 && numunderscores == 1) {
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Bad accession %s", tsip->accession);
        }
      }
      if (hasGi && tsip != NULL && tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
        if (sip1->choice == SEQID_DDBJ && bsp->repr == Seq_repr_seg) {
          sev = SEV_WARNING;
          /*
          ValidErr (vsp, sev, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", tsip->name);
          */
        } else {
          sev = SEV_REJECT;
          ValidErr (vsp, sev, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", tsip->name);
        }
      }
      /* and keep going with additional test */
    case SEQID_PIR:
    case SEQID_SWISSPROT:
    case SEQID_PRF:
      tsip = (TextSeqIdPtr) sip1->data.ptrvalue;
      if (tsip != NULL && StringHasNoText (tsip->accession) && ISA_na (bsp->mol)) {
        if (bsp->repr != Seq_repr_seg || hasGi) {
          if (sip1->choice != SEQID_DDBJ || bsp->repr != Seq_repr_seg) {
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_LONG, sizeof (buf1) - 1);
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadSeqIdFormat, "Missing accession for %s", buf1);
          }
        }
      }
      if (tsip != NULL && StringHasNoText (tsip->accession) &&
          StringHasNoText (tsip->name) && ISA_aa (bsp->mol)) {
        if (sip1->choice == SEQID_PIR || sip1->choice == SEQID_SWISSPROT || sip1->choice == SEQID_PRF) {
          SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_LONG, sizeof (buf1) - 1);
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_BadSeqIdFormat, "Missing identifier for %s", buf1);
        }
      }
      accn_count++;
      break;
    case SEQID_GPIPE:
      break;
    case SEQID_PATENT:
      isPatent = TRUE;
      break;
    case SEQID_PDB:
      isPDB = TRUE;
      break;
    case SEQID_GI:
      if (sip1->data.intvalue <= 0) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_ZeroGiNumber, "Invalid GI number");
      }
      gi_count++;
      break;
    case SEQID_GENERAL:
      break;
    default:
      break;
    }
  }

  if (gi_count > 0 && accn_count == 0 && (! isPDB) && bsp->repr != Seq_repr_virtual) {
    if (vsp->seqSubmitParent) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_INST_GiWithoutAccession, "No accession on sequence with gi number");
  }
  if (gi_count > 0 && accn_count > 1) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MultipleAccessions, "Multiple accessions on sequence with gi number");
  }

  /* optionally check IDs against older version in database */

  if (vsp->validateIDSet) {
    ValidateIDSetAgainstDb (gcp, vsp, bsp);
  }

  if (vsp->useSeqMgrIndexes) {
    vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
    while (vnp != NULL) {
      uop = (UserObjectPtr) vnp->data.ptrvalue;
      if (uop != NULL) {
        oip = uop->type;
        if (oip != NULL && StringICmp (oip->str, "TpaAssembly") == 0) {
          if (! IsTpa (bsp, TRUE, &isRefSeq)) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, 40);
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Non-TPA record %s should not have TpaAssembly object", buf1);
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        } else if (oip != NULL && StringICmp (oip->str, "RefGeneTracking") == 0) {
          hasRefGeneTracking = TRUE;
          hasRefTrackStatus = FALSE;
          hasLegalStatus = FALSE;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip != NULL && StringCmp (oip->str, "Status") == 0) {
              hasRefTrackStatus = TRUE;
              str = (CharPtr) ufp->data.ptrvalue;
              if (StringHasNoText (str)) {
                str = "?";
              }
              for (i = 0; legal_refgene_status_strings [i] != NULL; i++) {
                if (StringICmp (str, legal_refgene_status_strings [i]) == 0) {
                  hasLegalStatus = TRUE;
                  break;
                }
              }
              if (! hasLegalStatus) {
                olditemid = gcp->itemID;
                olditemtype = gcp->thistype;
                gcp->itemID = context.itemID;
                gcp->thistype = OBJ_SEQDESC;
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingIllegalStatus, "RefGeneTracking object has illegal Status '%s'", str);
                gcp->itemID = olditemid;
                gcp->thistype = olditemtype;
              }
            }
          }
          if (! hasRefTrackStatus) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            gcp->itemID = context.itemID;
            gcp->thistype = OBJ_SEQDESC;
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingWithoutStatus, "RefGeneTracking object needs to have Status set");
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
          if (! isRefSeq) {
            if (! IsWgsIntermediate (vsp->sep)) {
              olditemid = gcp->itemID;
              olditemtype = gcp->thistype;
              gcp->itemID = context.itemID;
              gcp->thistype = OBJ_SEQDESC;
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_RefGeneTrackingOnNonRefSeq, "RefGeneTracking object should only be in RefSeq record");
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
      vnp = SeqMgrGetNextDescriptor (bsp, vnp, Seq_descr_user, &context);
    }
  }

  for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
    bsp2 = BioseqFindSpecial (sip1);
    if (bsp2 == NULL) {
      if (!isPatent) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_IdOnMultipleBioseqs, "BioseqFind (%s) unable to find itself - possible internal error", buf1);
      }
    } else if (bsp2 != bsp) {
      if (sip1->choice == SEQID_GENERAL) {
        dbt = (DbtagPtr) sip1->data.ptrvalue;
        if (dbt != NULL && StringICmp (dbt->db, "NCBIFILE") == 0) continue;
      }
      SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_IdOnMultipleBioseqs, "SeqID %s is present on multiple Bioseqs in record", buf1);
    }
  }

  for (i = 0; i < 4; i++)
    errors[i] = FALSE;

  switch (bsp->repr) {
  case Seq_repr_virtual:
    if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
      errors[0] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_map:
    if ((bsp->seq_ext_type != 3) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_ref:
    if ((bsp->seq_ext_type != 2) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_seg:
    if ((bsp->seq_ext_type != 1) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  case Seq_repr_raw:
  case Seq_repr_const:
    if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
      errors[0] = TRUE;
    if ((bsp->seq_data_type < 1) || (bsp->seq_data_type > 11)
        || (bsp->seq_data == NULL))
      errors[2] = TRUE;
    break;
  case Seq_repr_delta:
    if ((bsp->seq_ext_type != 4) || (bsp->seq_ext == NULL))
      errors[1] = TRUE;
    if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
      errors[3] = TRUE;
    break;
  default:
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_ReprInvalid, "Invalid Bioseq->repr = %d", (int) (bsp->repr));
    return;
  }

  if (errors[0] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Bioseq-ext not allowed on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[1] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtBadOrMissing, "Missing or incorrect Bioseq-ext on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[2] == TRUE) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataNotFound, "Missing Seq-data on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (errors[3] == TRUE) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataNotAllowed, "Seq-data not allowed on %s Bioseq", repr[bsp->repr - 1]);
    retval = FALSE;
  }

  if (!retval)
    return;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  if (ISA_aa (bsp->mol)) {
    if (bsp->topology > 1) {    /* not linear */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_CircularProtein, "Non-linear topology set on protein");
    }
    if (bsp->strand > 1) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_DSProtein, "Protein not single stranded");
    }

  } else {
    if (!bsp->mol)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolNotSet, "Bioseq.mol is 0");
    else if (bsp->mol == Seq_mol_other)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolOther, "Bioseq.mol is type other");
    else if (bsp->mol == Seq_mol_na)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MolNuclAcid, "Bioseq.mol is type na");
  }

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;

  /* check sequence alphabet */
  if (((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const)) && bsp->seq_data_type != Seq_code_gap) {
    if (bsp->fuzz != NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_FuzzyLen, "Fuzzy length on %s Bioseq", repr[bsp->repr - 1]);
    }

    if (bsp->length < 1) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_InvalidLen, "Invalid Bioseq length [%ld]", (long) bsp->length);
    }

    seqtype = (int) (bsp->seq_data_type);
    switch (seqtype) {
    case Seq_code_iupacna:
    case Seq_code_ncbi2na:
    case Seq_code_ncbi4na:
    case Seq_code_ncbi8na:
    case Seq_code_ncbipna:
      if (ISA_aa (bsp->mol)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using a nucleic acid alphabet on a protein sequence");
        return;
      }
      break;
    case Seq_code_iupacaa:
    case Seq_code_ncbi8aa:
    case Seq_code_ncbieaa:
    case Seq_code_ncbipaa:
    case Seq_code_ncbistdaa:
      if (ISA_na (bsp->mol)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using a protein alphabet on a nucleic acid");
        return;
      }
      break;
    case Seq_code_gap:
      break;
    default:
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d]", (int) bsp->seq_data_type);
      return;
    }

    check_alphabet = FALSE;
    switch (seqtype) {
    case Seq_code_iupacaa:
    case Seq_code_iupacna:
    case Seq_code_ncbieaa:
    case Seq_code_ncbistdaa:
      check_alphabet = TRUE;

    case Seq_code_ncbi8na:
    case Seq_code_ncbi8aa:
      divisor = 1;
      break;

    case Seq_code_ncbi4na:
      divisor = 2;
      break;

    case Seq_code_ncbi2na:
      divisor = 4;
      break;

    case Seq_code_ncbipna:
      divisor = 5;
      break;

    case Seq_code_ncbipaa:
      divisor = 21;
      break;
    }

    len = bsp->length;
    if (len % divisor)
      len += divisor;
    len /= divisor;
    len2 = BSLen ((ByteStorePtr) bsp->seq_data);
    if (len > len2) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len2 * divisor),
                (long) bsp->length);
      return;
    } else if (len < len2) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len2 * divisor),
                (long) bsp->length);
    }

    if (check_alphabet) {       /* check 1 letter alphabets */
      switch (seqtype) {
      case Seq_code_iupacaa:
      case Seq_code_ncbieaa:
        termination = '*';
        gapchar = '-';
        break;
      case Seq_code_ncbistdaa:
        termination = 25;
        gapchar = 0;
        break;
      default:
        termination = '\0';
        gapchar = '\0';
        break;
      }
      if (! StreamCacheSetup (bsp, NULL, STREAM_EXPAND_GAPS, &sc)) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open StreamCache");
        return;
      }
      /*
      spp = SeqPortNew (bsp, 0, -1, 0, 0);
      if (spp == NULL) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open SeqPort");
        return;
      }
      */
      i = 0;
      terminations = 0;
      dashes = 0;
      gapatstart = FALSE;
      trailingX = 0;
      leadingX = FALSE;
      isLower = FALSE;
      isFirst = TRUE;
      for (len = 0; len < bsp->length; len++) {
        residue = StreamCacheGetResidue (&sc);
        /*
        residue = SeqPortGetResidue (spp);
        */
        if (!IS_residue (residue)) {
          i++;
          if (i > 10) {
            ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "More than 10 invalid residues. Checking stopped");
            /*
            SeqPortFree (spp);
            */
            if (vsp->patch_seq)
              PatchBadSequence (bsp);
            return;
          } else {
            BSSeek ((ByteStorePtr) bsp->seq_data, len, SEEK_SET);
            x = BSGetByte ((ByteStorePtr) bsp->seq_data);
            if (bsp->seq_data_type == Seq_code_ncbistdaa) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) x, (long) (len + 1));
            } else if (IS_ALPHA (x)) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue '%c' at position [%ld]", (char) x, (long) (len + 1));
            } else {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) x, (long) (len + 1));
            }
          }
        } else if (residue == termination) {
          terminations++;
          trailingX = 0;        /* suppress if followed by terminator */
        } else if (residue == gapchar) {
          dashes++;
          if (len == 0) {
            gapatstart = TRUE;
          }
        } else if (residue == 'X') {
          trailingX++;
          if (isFirst) {
            leadingX = TRUE;
          }
        } else if (ISA_na (bsp->mol) && StringChr ("EFIJLPQZ", (Char) residue) != NULL) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid nucleotide residue '%c' at position [%ld]", (char) residue, (long) (len + 1));
        } else if (! IS_ALPHA ((Char) residue)) {
          ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue '%c' at position [%ld]", (char) residue, (long) (len + 1));
        } else {
          trailingX = 0;
          if (IS_LOWER ((Char) residue)) {
            isLower = TRUE;
          }
        }
        isFirst = FALSE;
      }
      /*
      SeqPortFree (spp);
      */
      if (ISA_aa (bsp->mol) && (leadingX || trailingX > 0)) {
        /* only show leading or trailing X if product of NNN in nucleotide */
        cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
        if (cds != NULL) {
          crp = (CdRegionPtr) cds->data.value.ptrvalue;
          if (crp != NULL) {
            dnalen = SeqLocLen (cds->location);
            if (dnalen > 5) {
              bases = ReadCodingRegionBases (cds->location, dnalen, crp->frame, &total);
              len = StringLen (bases);
              if (len > 5) {
                if (StringNICmp (bases, "NNN", 3) != 0) {
                  leadingX = FALSE;
                }
                if (StringNICmp (bases + len - 3, "NNN", 3) != 0) {
                  trailingX = 0;
                }
              }
              MemFree (bases);
            }
          }
        }
      }
      if (leadingX) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LeadingX, "Sequence starts with leading X", (int) leadingX);
      }
      if (trailingX > 0 && SuppressTrailingXMessage (bsp)) {
        /* suppress if cds translation ends in '*' or 3' partial */
      } else if (trailingX > 1) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TrailingX, "Sequence ends in %d trailing Xs", (int) trailingX);
      } else if (trailingX > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_TrailingX, "Sequence ends in %d trailing X", (int) trailingX);
      }
      if (isLower) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Sequence contains lower-case characters");
      }
      if (terminations > 0 || dashes > 0) {
        cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
        grp = SeqMgrGetGeneXref (cds);
        genelbl = NULL;
        if (grp == NULL && cds != NULL) {
          gene = SeqMgrGetOverlappingGene (cds->location, &genectxt);
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
          }
        }
        if (grp != NULL && (!SeqMgrGeneIsSuppressed (grp))) {
          if (grp->locus != NULL)
            genelbl = (grp->locus);
          else if (grp->locus_tag != NULL)
            genelbl = (grp->locus_tag);
          else if (grp->desc != NULL)
            genelbl = (grp->desc);
          else if (grp->syn != NULL)
            genelbl = (CharPtr) (grp->syn->data.ptrvalue);
        }
        prot = SeqMgrGetBestProteinFeature (bsp, &protctxt);
        protlbl = protctxt.label;
      }
      if (StringHasNoText (genelbl)) {
        genelbl = "gene?";
      }
      if (StringHasNoText (protlbl)) {
        protlbl = "prot?";
      }
      if (dashes > 0) {
        if (gapatstart && dashes == 1) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadProteinStart, "gap symbol at start of protein sequence (%s - %s)", genelbl, protlbl);
        } else if (gapatstart) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadProteinStart, "gap symbol at start of protein sequence (%s - %s)", genelbl, protlbl);
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_GapInProtein, "[%d] internal gap symbols in protein sequence (%s - %s)", (dashes - 1), genelbl, protlbl);
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_GapInProtein, "[%d] internal gap symbols in protein sequence (%s - %s)", dashes, genelbl, protlbl);
        }
      }
      if (terminations) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_StopInProtein, "[%d] termination symbols in protein sequence (%s - %s)", terminations, genelbl, protlbl);
        if (!i)
          return;
      }
      if (i) {
        if (vsp->patch_seq)
          PatchBadSequence (bsp);
        return;
      }

    }
  }

  if (ISA_na (bsp->mol) && bsp->repr == Seq_repr_delta && DeltaLitOnly (bsp)) {
    if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqPortFail, "Can't open StreamCache");
      return;
    }
    in_gap = FALSE;
    in_N = FALSE;
    adjacent_N_and_gap = FALSE;
    adjacent_N_gap_position = 0;
    for (len = 0; len < bsp->length; len++) {
      residue = StreamCacheGetResidue (&sc);
      if (residue == '-') {
        if (in_N) {
          adjacent_N_and_gap = TRUE;
          if (adjacent_N_gap_position == 0) {
            adjacent_N_gap_position = len;
          }
        }
        in_N = FALSE;
        in_gap = TRUE;
      } else if (residue == 'N') {
        if (in_gap) {
          adjacent_N_and_gap = TRUE;
          if (adjacent_N_gap_position == 0) {
            adjacent_N_gap_position = len;
          }
        }
        in_gap = FALSE;
        in_N = TRUE;
      } else {
        in_gap = FALSE;
        in_N = FALSE;
      }
    }
    if (adjacent_N_and_gap) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_InternalNsAdjacentToGap,
                "Ambiguous residue N is adjacent to a gap around position %ld",
                (long) adjacent_N_gap_position);
    }
  }

  if ((bsp->repr == Seq_repr_seg) || (bsp->repr == Seq_repr_ref)) {     /* check segmented sequence */
    head.choice = SEQLOC_MIX;
    head.data.ptrvalue = bsp->seq_ext;
    head.next = NULL;
    ValidateSeqLoc (vsp, (SeqLocPtr) & head, "Segmented Bioseq");
    /* check the length */
    len = 0;
    vnp = NULL;
    while ((vnp = SeqLocFindNext (&head, vnp)) != NULL) {
      len2 = SeqLocLen (vnp);
      if (len2 > 0)
        len += len2;
    }
    if (bsp->length > len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len), (long) bsp->length);
    } else if (bsp->length < len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len), (long) bsp->length);
    }

    vnp = NULL;
    idlist = NULL;
    while ((vnp = SeqLocFindNext (&head, vnp)) != NULL) {
      sip1 = SeqLocId (vnp);
      if (sip1 != NULL) {
        SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
        ValNodeCopyStr (&idlist, vnp->choice, buf1);
      }
    }
    if (idlist != NULL) {
      idlist = ValNodeSort (idlist, SortVnpByString);
      last = (CharPtr) idlist->data.ptrvalue;
      lastchoice = (Uint1) idlist->choice;
      vnp = idlist->next;
      while (vnp != NULL) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringICmp (last, str) == 0) {
          if (vnp->choice == lastchoice && lastchoice == SEQLOC_WHOLE) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_DuplicateSegmentReferences, "Segmented sequence has multiple references to %s", str);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_DuplicateSegmentReferences,
                      "Segmented sequence has multiple references to %s that are not SEQLOC_WHOLE\n", str);
          }
        } else {
          last = (CharPtr) vnp->data.ptrvalue;
          lastchoice = (Uint1) vnp->choice;
        }
        vnp = vnp->next;
      }
      ValNodeFreeData (idlist);
    }

    vsp->bsp_partial_val = SeqLocPartialCheck ((SeqLocPtr) (&head));
    if ((vsp->bsp_partial_val) && (ISA_aa (bsp->mol))) {
      bcp = NULL;
      vnp = NULL;
      got_partial = FALSE;
      if (vsp->useSeqMgrIndexes) {
        vnp = SeqMgrGetNextDescriptor (bsp, vnp, Seq_descr_molinfo, &context);
      } else {
        bcp = BioseqContextNew (bsp);
        vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, vnp, NULL);
      }
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
        if (mip != NULL) {
          switch (mip->completeness) {
          case 2:             /* partial */
            got_partial = TRUE;
            break;
          case 3:             /* no-left */
            if (!(vsp->bsp_partial_val & SLP_START))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-left inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          case 4:             /* no-right */
            if (!(vsp->bsp_partial_val & SLP_STOP))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-right inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          case 5:             /* no-ends */
            if ((!(vsp->bsp_partial_val & SLP_STOP)) && (!(vsp->bsp_partial_val & SLP_START)))
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "No-ends inconsistent with segmented SeqLoc");
            got_partial = TRUE;
            break;
          default:
            break;
          }
        }
      }
      if (!vsp->useSeqMgrIndexes) {
        BioseqContextFree (bcp);
      }
      if (!got_partial)
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "Partial segmented sequence without MolInfo partial");
    }
  }

  mip = NULL;

  if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_raw) {

    vnp = NULL;
    if (vsp->useSeqMgrIndexes) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
    } else {
      bcp = BioseqContextNew (bsp);
      vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
      BioseqContextFree (bcp);
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
    }

    vnp = NULL;
    if (vsp->useSeqMgrIndexes) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
    } else {
      bcp = BioseqContextNew (bsp);
      vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_genbank, NULL, NULL);
      BioseqContextFree (bcp);
    }
    if (vnp != NULL) {
      gbp = (GBBlockPtr) vnp->data.ptrvalue;
      if (gbp != NULL) {
        for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringICmp (str, "HTGS_ACTIVEFIN") == 0) {
            isActiveFin = TRUE;
          } else if (StringICmp (str, "HTGS_DRAFT") == 0) {
            isDraft = TRUE;
          } else if (StringICmp (str, "HTGS_FULLTOP") == 0) {
            isFullTop = TRUE;
          } else if (StringICmp (str, "HTGS_PREFIN") == 0) {
            isPreFin = TRUE;
          }
        }
      }
    }
  }

  if (bsp->repr == Seq_repr_delta) {
    len = 0;
    for (vnp = (ValNodePtr) (bsp->seq_ext), segnum = 1; vnp != NULL; vnp = vnp->next, segnum++) {
      if (vnp->data.ptrvalue == NULL)
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "NULL pointer in delta seq_ext valnode");
      else {
        switch (vnp->choice) {
        case 1:                /* SeqLocPtr */
          slp = (SeqLocPtr) (vnp->data.ptrvalue);
          if (slp != NULL && slp->choice == SEQLOC_WHOLE) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_WholeComponent, "Delta seq component should not be of type whole");
          }
          sip3 = SeqLocId (slp);
          if (sip3 != NULL) {
            if (sip3->choice == SEQID_GI && sip3->data.intvalue <= 0) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_DeltaComponentIsGi0, "Delta component is gi|0");
            }
            for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
              if (SeqIdComp (sip1, sip3) == SIC_YES) {
                ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SelfReferentialSequence,
                          "Self-referential delta sequence");
              }
            }
          }
          len2 = SeqLocLen (slp);
          if (len2 < 0)
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "-1 length on seq-loc of delta seq_ext");
          else
            len += len2;
          sip3 = SeqLocId (slp);
          if (sip3 != NULL && slp != NULL && slp->choice == SEQLOC_INT) {
            sintp = (SeqIntPtr) slp->data.ptrvalue;
            if (sintp != NULL && (sip3->choice == SEQID_GI ||
                                 sip3->choice == SEQID_GENBANK ||
                                 sip3->choice == SEQID_EMBL ||
                                 sip3->choice == SEQID_DDBJ ||
                                 sip3->choice == SEQID_TPG ||
                                 sip3->choice == SEQID_TPE ||
                                 sip3->choice == SEQID_TPD ||
                                 sip3->choice == SEQID_OTHER)) {
              vn.choice = SEQLOC_WHOLE;
              vn.data.ptrvalue = sip3;
              vn.next = NULL;
              len3 = SeqLocLen (&vn);
              /* -1 signifies failure to lookup or not connected to lookup function */
              if (len3 != -1) {
                if (sintp->to >= len3) {
                  SeqIdWrite (sip3, buf1, PRINTID_FASTA_SHORT, 40);
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong,
                            "Seq-loc extent (%ld) greater than length of %s (%ld)",
                            (long) (sintp->to + 1), buf1, (long) len3);
                }
              }
            }
          }
          if (len2 <= 10) {
            str = SeqLocPrint ((SeqLocPtr) (vnp->data.ptrvalue));
            if (str == NULL) {
              str = StringSave ("?");
            }
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SeqLocLength, "Short length (%ld) on seq-loc (%s) of delta seq_ext", (long) len2, str);
            MemFree (str);
          }
          break;
        case 2:                /* SeqLitPtr */
          slitp = (SeqLitPtr) (vnp->data.ptrvalue);
          if (slitp->seq_data != NULL && slitp->seq_data_type != Seq_code_gap) {
            if (slitp->length == 0) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_SeqLitDataLength0, "Seq-lit of length 0 in delta chain");
            }
            sctp = SeqCodeTableFind (slitp->seq_data_type);
            if (sctp == NULL) {
              ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d] in SeqLitPtr", (int) slitp->seq_data_type);
              len += slitp->length;
              break;
            }

            start_at = (Int2) (sctp->start_at);
            num = (Int2) (sctp->num);

            switch (slitp->seq_data_type) {
            case Seq_code_iupacaa:
            case Seq_code_iupacna:
            case Seq_code_ncbieaa:
            case Seq_code_ncbistdaa:
              BSSeek ((ByteStorePtr) slitp->seq_data, 0, SEEK_SET);
              for (len2 = 1; len2 <= (slitp->length); len2++) {
                is_invalid = FALSE;
                residue = BSGetByte ((ByteStorePtr) slitp->seq_data);
                i = residue - start_at;
                if ((i < 0) || (i >= num))
                  is_invalid = TRUE;
                else if (*(sctp->names[i]) == '\0')
                  is_invalid = TRUE;
                if (is_invalid) {
                  if (slitp->seq_data_type == Seq_code_ncbistdaa)
                    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] at position [%ld]", (int) residue, (long) (len + len2));
                  else
                    ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%c] at position [%ld]", (char) residue, (long) (len + len2));
                }
              }
              break;
            default:
              break;
            }
            if (mip != NULL) {
              if (mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
                runsofn = CountAdjacentNsInSeqLit (gcp, slitp, (Boolean) ISA_na (bsp->mol));
                if (runsofn > 80) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %ld that starts at base %ld", (long) runsofn, (int) segnum, (long) (len + 1));
                }
              } else if (mip->tech == MI_TECH_wgs) {
                runsofn = CountAdjacentNsInSeqLit (gcp, slitp, (Boolean) ISA_na (bsp->mol));
                if (runsofn >= 20) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %ld that starts at base %ld", (long) runsofn, (int) segnum, (long) (len + 1));
                }
              } else if (mip->tech == MI_TECH_composite_wgs_htgs) {
                runsofn = CountAdjacentNsInSeqLit (gcp, slitp, (Boolean) ISA_na (bsp->mol));
                if (runsofn > 80) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %ld that starts at base %ld", (long) runsofn, (int) segnum, (long) (len + 1));
                }
              } else {
                runsofn = CountAdjacentNsInSeqLit (gcp, slitp, (Boolean) ISA_na (bsp->mol));
                if (runsofn > 100) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqLit, "Run of %ld Ns in delta component %ld that starts at base %ld", (long) runsofn, (int) segnum, (long) (len + 1));
                }
              }
            }
          } else if (slitp->length == 0) {
            if (isSwissProt) {
              sev = SEV_WARNING;
            } else {
              sev = SEV_ERROR;
            }
            ifp = slitp->fuzz;
            if (ifp == NULL || ifp->choice != 4 || ifp->a != 0) {
              ValidErr (vsp, sev, ERR_SEQ_INST_SeqLitGapLength0, "Gap of length 0 in delta chain");
            } else {
              ValidErr (vsp, sev, ERR_SEQ_INST_SeqLitGapLength0, "Gap of length 0 with unknown fuzz in delta chain");
            }
          }
          len += slitp->length;
          break;
        default:
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Illegal choice [%d] in delta chain", (int) (vnp->choice));
          break;
        }
      }
    }
    if (bsp->length > len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data too short [%ld] for given length [%ld]", (long) (len), (long) bsp->length);
    } else if (bsp->length < len) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_INST_SeqDataLenWrong, "Bioseq.seq_data is larger [%ld] than given length [%ld]", (long) (len), (long) bsp->length);
    }
    if (mip != NULL) {
      is_gps = FALSE;
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
          is_gps = TRUE;
        }
      }
      if ((!isNTorNC) && (! is_gps) && mip->tech != MI_TECH_htgs_0 && mip->tech != MI_TECH_htgs_1 &&
          mip->tech != MI_TECH_htgs_2 && mip->tech != MI_TECH_htgs_3 && mip->tech != MI_TECH_wgs &&
          mip->tech != MI_TECH_composite_wgs_htgs && mip->tech != MI_TECH_unknown && mip->tech != MI_TECH_standard
          && mip->tech != MI_TECH_htc && mip->tech != MI_TECH_barcode && mip->tech != MI_TECH_tsa) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadDeltaSeq, "Delta seq technique should not be [%d]", (int) (mip->tech));
      }
    }
  } else if (bsp->repr == Seq_repr_raw) {
    ron.gcp = gcp;
    ron.vsp = vsp;
    ron.ncount = 0;
    ron.maxrun = 0;
    ron.seqpos = 0;
    ron.gapcount = 0;
    ron.showAll = TRUE;
    ron.inNrun = FALSE;
    ron.isWGS = FALSE;
    if (mip == NULL) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
      }
    }
    if (mip != NULL && mip->tech == MI_TECH_wgs) {
      ron.isWGS = TRUE;
    }

    SeqPortStream (bsp, EXPAND_GAPS_TO_DASHES, (Pointer) &ron, CountAdjacentProc);

    /*
    if (ron.inNrun && ron.showAll && ron.ncount >= 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence starting at base %ld",
                (long) ron.ncount, (long) (ron.seqpos - ron.ncount + 1));
    }
    */

    if (ron.gapcount > 0 && ISA_na (bsp->mol)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalGapsInSeqRaw, "Raw nucleotide should not contain gap characters");
    }

    /*
    if (ron.maxrun >= 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_InternalNsInSeqRaw, "Run of %ld Ns in raw sequence", (long) ron.maxrun);
    }
    */
  }

  if (bsp->repr == Seq_repr_delta) {
    CheckDeltaForReuse (vsp, gcp, bsp);
  }

  sev = SEV_ERROR;
  if (mip != NULL) {
    if (mip->tech != MI_TECH_htgs_0 && mip->tech != MI_TECH_htgs_1 &&
        mip->tech != MI_TECH_htgs_2 && mip->tech != MI_TECH_htgs_3) {
      sev = SEV_WARNING;
    }
  }

  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4 && bsp->seq_ext != NULL) {
    vnp = (DeltaSeqPtr) bsp->seq_ext;
    if (vnp != NULL && vnp->choice == 2) {
      slitp = (SeqLitPtr) vnp->data.ptrvalue;
      if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
        ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "First delta seq component is a gap");
      }
    }
    last_is_gap = FALSE;
    num_adjacent_gaps = 0;
    num_gaps = 0;
    non_interspersed_gaps = FALSE;
    while (vnp->next != NULL) {
      vnp = vnp->next;
      if (vnp != NULL && vnp->choice == 2) {
        slitp = (SeqLitPtr) vnp->data.ptrvalue;
        if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
          if (last_is_gap) {
            num_adjacent_gaps++;
          }
          last_is_gap = TRUE;
          num_gaps++;
        } else {
          if (! last_is_gap) {
            non_interspersed_gaps = TRUE;
          }
          last_is_gap = FALSE;
        }
      } else {
        if (! last_is_gap) {
          non_interspersed_gaps = TRUE;
        }
        last_is_gap = FALSE;
      }
    }
    if (non_interspersed_gaps && (! hasGi) && mip != NULL &&
        (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2 /* || mip->tech == MI_TECH_htgs_3 */)) {
      if (hasRefGeneTracking) {  
        ValidErr (vsp, SEV_INFO, ERR_SEQ_INST_MissingGaps, "HTGS delta seq should have gaps between all sequence runs");
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_MissingGaps, "HTGS delta seq should have gaps between all sequence runs");
      }
    }
    if (num_adjacent_gaps > 1) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadDeltaSeq, "There are %d adjacent gaps in delta seq", (int) num_adjacent_gaps);
    } else if (num_adjacent_gaps > 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadDeltaSeq, "There is %d adjacent gap in delta seq", (int) num_adjacent_gaps);
    }
    if (vnp != NULL && vnp->choice == 2) {
      slitp = (SeqLitPtr) vnp->data.ptrvalue;
      if (slitp != NULL && (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap)) {
        ValidErr (vsp, sev, ERR_SEQ_INST_BadDeltaSeq, "Last delta seq component is a gap");
      }
    }
    if (num_gaps == 0 && mip != NULL) {
      if (/* mip->tech == MI_TECH_htgs_1 || */ mip->tech == MI_TECH_htgs_2) {
        if (VisitGraphsInSep (sep, NULL, NULL) == 0) {
          if (! isActiveFin) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadHTGSeq, "HTGS 2 delta seq has no gaps and no graphs");
          }
        }
      }
    }
  }

  if (bsp->repr == Seq_repr_raw) {
    if (mip != NULL) {
      if (/* mip->tech == MI_TECH_htgs_1 || */ mip->tech == MI_TECH_htgs_2) {
        if (VisitGraphsInSep (sep, NULL, NULL) == 0) {
          if (! isActiveFin) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_BadHTGSeq, "HTGS 2 raw seq has no gaps and no graphs");
          }
        }
      }
    }
  }

  if (mip != NULL && mip->tech == MI_TECH_htgs_3) {
    if (isDraft) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_DRAFT keyword");
    }
    if (isPreFin) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_PREFIN keyword");
    }
    if (isActiveFin) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_ACTIVEFIN keyword");
    }
    if (isFullTop) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_BadHTGSeq, "HTGS 3 sequence should not have HTGS_FULLTOP keyword");
    }
  }

  if (ISA_aa (bsp->mol)) {
    if ((bsp->length <= 3) && (bsp->length >= 0) && (!isPDB)) {
      if (mip == NULL || mip->completeness < 2 || mip->completeness > 5) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long) (bsp->length));
      }
    }

  } else {
    if ((bsp->length <= 10) && (bsp->length >= 0) && (!isPDB)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long) (bsp->length));
    }
  }

#if 0
  if (bsp->length > 350000 && (! isNTorNC)) {
    Boolean         isGenBankEMBLorDDBJ;
    Boolean         litHasData;
    if (bsp->repr == Seq_repr_delta) {
      isGenBankEMBLorDDBJ = FALSE;
      /* suppress this for data from genome annotation project */
      VisitBioseqsInSep (vsp->sep, (Pointer) &isGenBankEMBLorDDBJ, LookForGEDseqID);
      if (mip != NULL && isGenBankEMBLorDDBJ) {
        if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence, "Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_htgs_3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Phase 3 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_wgs) {
          /*
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "WGS sequence exceeds 350kbp limit");
          */
        } else {
          len = 0;
          litHasData = FALSE;
          for (vnp = (ValNodePtr) (bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
            if (vnp->choice == 2) {
              slitp = (SeqLitPtr) (vnp->data.ptrvalue);
              if (slitp != NULL) {
                if (slitp->seq_data != NULL) {
                  litHasData = TRUE;
                }
                len += slitp->length;
              }
            }
          }
          if (len > 500000 && litHasData) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_LongLiteralSequence, "Length of sequence literals exceeds 500kbp limit");
          }
        }
      }
    } else if (bsp->repr == Seq_repr_raw) {
      vnp = NULL;
      if (vsp->useSeqMgrIndexes) {
        vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
      } else {
        bcp = BioseqContextNew (bsp);
        vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
        BioseqContextFree (bcp);
      }
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
      }
      if (mip != NULL) {
        if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence, "Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_htgs_3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Phase 3 HTGS sequence exceeds 350kbp limit");
        } else if (mip->tech == MI_TECH_wgs) {
          /*
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "WGS sequence exceeds 350kbp limit");
          */
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Length of sequence exceeds 350kbp limit");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp, "Length of sequence exceeds 350kbp limit");
      }
    } else {
      /* Could be a segset header bioseq that is > 350kbp */
      /* No-op for now? Or generate a warning? */
    }
  }
#endif

  if (bsp->repr == Seq_repr_seg) {
    CheckSegBspAgainstParts (vsp, gcp, bsp);
  }

  if (ISA_na (bsp->mol) || ISA_aa (bsp->mol)) {
    vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
    if (vnp != NULL) {
      title = (CharPtr) vnp->data.ptrvalue;
      if (StringDoesHaveText (title)) {
        if (HasUnparsedBrackets (title)) {
          reportFastaBracket = TRUE;
          for (sip = bsp->id; sip != NULL; sip = sip->next) {
            if (sip->choice != SEQID_GENERAL) continue;
            dbt = (DbtagPtr) sip->data.ptrvalue;
            if (dbt == NULL) continue;
            if (StringICmp (dbt->db, "TMSMART") == 0) {
              reportFastaBracket = FALSE;
            }
            if (StringICmp (dbt->db, "BankIt") == 0) {
              reportFastaBracket = FALSE;
            }
          }
          if (reportFastaBracket) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            if (vnp->extended != 0) {
              ovp = (ObjValNodePtr) vnp;
              gcp->itemID = ovp->idx.itemID;
              gcp->thistype = OBJ_SEQDESC;
            }
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_FastaBracketTitle, "Title may have unparsed [...=...] construct");
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
      }
    }
  }

  if (ISA_aa (bsp->mol) && vsp->useSeqMgrIndexes) {
    vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
    if (vnp != NULL) {
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bsp->idx.parentptr;
        while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot) {
          if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bssp->idx.parentptr;
          } else {
            bssp = NULL;
          }
        }
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
          title = (CharPtr) vnp->data.ptrvalue;
          tech = 0;
          sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
          if (sdp != NULL) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              tech = mip->tech;
            }
          }
          buf = MemNew (sizeof (Char) * (buflen + 1));
          MemSet ((Pointer) (&ii), 0, sizeof (ItemInfo));
          /* check generated protein defline with first prp->name - new convention */
          if (buf != NULL && NewCreateDefLineBuf (&ii, bsp, buf, buflen, TRUE, FALSE)) {
            if (StringICmp (buf, title) != 0) {
              /* okay if instantiated title has single trailing period */
              len2 = StringLen (buf);
              len3 = StringLen (title);
              if (len3 == len2 + 1 && title [len3 - 1] == '.' && len3 > 3 && title [len3 - 2] != '.') {
                StringCat (buf, ".");
              }
            }
            if (StringICmp (buf, title) != 0) {
              /* also check generated protein defline with all prp->names - old convention */
              if (NewCreateDefLineBuf (&ii, bsp, buf, buflen, TRUE, TRUE)) {
                if (StringICmp (buf, title) != 0) {
                  olditemid = gcp->itemID;
                  olditemtype = gcp->thistype;
                  if (vnp->extended != 0) {
                    ovp = (ObjValNodePtr) vnp;
                    gcp->itemID = ovp->idx.itemID;
                    gcp->thistype = OBJ_SEQDESC;
                  }
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InconsistentProteinTitle,
                            "Instantiated protein title does not match automatically generated title");
                  gcp->itemID = olditemid;
                  gcp->thistype = olditemtype;
                }
              }
            }
          }
          MemFree (buf);
        }
      }
    }
  }

  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (vnp != NULL) {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness != 1 && isGB) {
        buf = MemNew (sizeof (Char) * (4097));
        if (buf != NULL && NewCreateDefLineBuf (NULL, bsp, buf, 4096, FALSE, FALSE)) {
          if (StringStr (buf, "complete genome") != NULL) {
            olditemid = gcp->itemID;
            olditemtype = gcp->thistype;
            if (vnp->extended != 0) {
              ovp = (ObjValNodePtr) vnp;
              gcp->itemID = ovp->idx.itemID;
              gcp->thistype = OBJ_SEQDESC;
            }
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteTitleProblem, "Complete genome in title without complete flag set");
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
        MemFree (buf);
      }
      if (mip->completeness != 1 && bsp->topology == 2) {
        olditemid = gcp->itemID;
        olditemtype = gcp->thistype;
        if (vnp->extended != 0) {
          ovp = (ObjValNodePtr) vnp;
          gcp->itemID = ovp->idx.itemID;
          gcp->thistype = OBJ_SEQDESC;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteCircleProblem, "Circular topology without complete flag set");
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
    }
  }

  if (ISA_na (bsp->mol) && (bsp->repr == Seq_repr_raw || (bsp->repr == Seq_repr_delta && DeltaLitOnly (bsp))) && bsp->length > 10) {
    /* check for N bases at start or stop of sequence */
    sfp = (SeqFeatPtr) MemNew (sizeof (SeqFeat));
    if (sfp == NULL) return;
    sfp->data.choice = SEQFEAT_COMMENT;

    sfp->location = AddIntervalToLocation (NULL, bsp->id, 0, 9, FALSE, FALSE);
    str = GetSequencePlusGapByFeature (sfp);
    if (str != NULL) {
      if (str [0] == 'n' || str [0] == 'N') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "NNNNNNNNNN") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalNs, "N at beginning of sequence");
      }
      if (str [0] == '-' || str [0] == '-') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "----------") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalGap, "Gap at beginning of sequence");
      }
    }
    MemFree (str);
    sfp->location = SeqLocFree (sfp->location);

    sfp->location = AddIntervalToLocation (NULL, bsp->id, bsp->length - 10, bsp->length - 1, FALSE, FALSE);
    str = GetSequencePlusGapByFeature (sfp);
    len = StringLen (str);
    if (str != NULL && len > 0) {
      if (str [len - 1] == 'n' || str [len - 1] == 'N') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "NNNNNNNNNN") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalNs, "N at end of sequence");
      }
      if (str [len - 1] == '-' || str [len - 1] == '-') {
        if (isNC || isPatent) {
          sev = SEV_WARNING;
        } else if (bsp->topology == TOPOLOGY_CIRCULAR) {
          sev = SEV_WARNING;
        } else if (StringICmp (str, "----------") == 0) {
          sev = SEV_ERROR;
        } else {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_INST_TerminalGap, "Gap at end of sequence");
      }
    }
    MemFree (str);
    sfp->location = SeqLocFree (sfp->location);

    MemFree (sfp);
  }
}

/*****************************************************************************
*
*   ValidatePubdesc(gcp)
*      Check pubdesc for missing information
*
*****************************************************************************/
static Boolean HasNoText (CharPtr str)
{
  Char            ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static Boolean HasNoName (ValNodePtr name)
{
  AuthorPtr       ap;
  NameStdPtr      nsp;
  PersonIdPtr     pid;

  if (name != NULL) {
    ap = name->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL) {
        if (pid->choice == 2) {
          nsp = pid->data;
          if (nsp != NULL) {
            if (!HasNoText (nsp->names[0])) {
              return FALSE;
            }
          }
        } else if (pid->choice == 5) {
          /* consortium */
          if (!HasNoText ((CharPtr) pid->data)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}

static void ValidateAffil (ValidStructPtr vsp, AffilPtr ap)

{
  if (ap != NULL) {
    if (ap->affil == NULL && ap->div == NULL && ap->street == NULL && ap->city == NULL &&
        ap->sub == NULL && ap->postal_code == NULL && ap->country == NULL &&
        ap->phone == NULL && ap->fax == NULL && ap->email == NULL) {
      /* no affiliation */
    } else {
      if (ap->choice == 2) {
        /*
        if (StringHasNoText (ap->city)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no city");
        }
        */
        if (StringHasNoText (ap->country)) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no country");
        }
        if (StringCmp (ap->country, "USA") == 0) {
          if (StringHasNoText (ap->sub)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no state");
          }
        }
      }
    }
  }
}

static void ValidateCitSub (ValidStructPtr vsp, CitSubPtr csp)
{
  AffilPtr        ap;
  AuthListPtr     alp;
  ValNodePtr      name;
  Boolean         hasAffil = FALSE;
  Boolean         hasName = FALSE;
  ErrSev          sev;

  if (vsp == NULL || csp == NULL) return;
  
  sev = SEV_ERROR;
  if (vsp->is_refseq_in_sep) {
    sev = SEV_WARNING;
  }
  if (vsp->is_htg_in_sep) {
    sev = SEV_WARNING;
  }

  alp = csp->authors;
  if (alp != NULL) {
    if (alp->choice == 1) {
      for (name = alp->names; name != NULL; name = name->next) {
        if (!HasNoName (name)) {
          hasName = TRUE;
        }
      }
    } else if (alp->choice == 2 || alp->choice == 3) {
      for (name = alp->names; name != NULL; name = name->next) {
        if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
          hasName = TRUE;
        }
      }
    }
    ap = alp->affil;
    if (ap != NULL) {
      if (ap->affil == NULL && ap->div == NULL && ap->street == NULL && ap->city == NULL &&
           ap->sub == NULL && ap->postal_code == NULL && ap->country == NULL &&
           ap->phone == NULL && ap->fax == NULL && ap->email == NULL) {
        /* no affiliation */
      } else {
        hasAffil = TRUE;
        if (ap->choice == 2) {
          /*
          if (StringHasNoText (ap->city)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no city");
          }
          */
          if (StringHasNoText (ap->country)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no country");
          }
          if (StringCmp (ap->country, "USA") == 0) {
            if (StringHasNoText (ap->sub)) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Submission citation affiliation has no state");
            }
          }
        }
      }
    }
  }
  if (!hasName) {
    ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Submission citation has no author names");
  }
  if (!hasAffil) {
    ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Submission citation has no affiliation");
  }
}

static void LookForMultiplePubs (ValidStructPtr vsp, GatherContextPtr gcp, SeqDescrPtr sdp)

{
  Bioseq       bs;
  Boolean      collision, otherpub;
  Int4         muid, pmid;
  SeqDescrPtr  nextpub;
  PubdescPtr   pdp;
  ValNodePtr   vnp;


  if (sdp != NULL && sdp->choice == Seq_descr_pub && sdp->extended != 0 && vsp != NULL && gcp != NULL) {
    MemSet ((Pointer) &bs, 0, sizeof (Bioseq));
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      otherpub = FALSE;
      muid = 0;
      pmid = 0;
      for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == PUB_Muid) {
          muid = vnp->data.intvalue;
        } else if (vnp->choice == PUB_PMid) {
          pmid = vnp->data.intvalue;
        } else {
          otherpub = TRUE;
        }
      }
      if (otherpub) {
        if (muid > 0 || pmid > 0) {
          collision = FALSE;
          nextpub = GetNextDescriptorUnindexed (&bs, Seq_descr_pub, sdp);
          while (nextpub != NULL) {
            pdp = (PubdescPtr) nextpub->data.ptrvalue;
            if (pdp != NULL) {
              for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == PUB_Muid) {
                  if (muid > 0 && muid == vnp->data.intvalue) {
                    collision = TRUE;
                  }
                } else if (vnp->choice == PUB_PMid) {
                  if (pmid > 0 && pmid == vnp->data.intvalue) {
                    collision = TRUE;
                  }
                }
              }
            }
            nextpub = GetNextDescriptorUnindexed (&bs, Seq_descr_pub, nextpub);
          }
          if (collision) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple publications with same identifier");
          }
        }
      }
    }
  }
}

static void LookForMultipleUnpubPubs (ValidStructPtr vsp, GatherContextPtr gcp, BioseqPtr bsp)

{
  Char               buf [2048];
  CharPtr            last, str;
  SeqMgrDescContext  dcontext;
  ValNodePtr         list = NULL, next, vnp;
  ObjValNodePtr      ovp;
  PubdescPtr         pdp;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.scratch != NULL) {
        ValNodeCopyStr (&list, 0, ovp->idx.scratch);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  if (list == NULL) return;

  list = ValNodeSort (list, SortVnpByString);
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      StringNCpy_0 (buf, str, sizeof (buf));
      StringCpy (buf + 100, "...");
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications,
                "Multiple equivalent publications annotated on this sequence [%s]", buf);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
    }
    vnp = next;
  }

  ValNodeFreeData (list);
}

static Boolean BadCharsInAuth (CharPtr str, CharPtr PNTR badauthor, Boolean allowcomma, Boolean allowperiod, Boolean last)
{
  Char     ch;
  CharPtr  ptr, stp = NULL;

  if (StringHasNoText (str)) return FALSE;
  if (last) {
      stp = StringISearch(str, "St.");
      if (stp == str) {
          stp += 2;  /* point to the period */
      }
  }

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    /* success on any of these tests are allowed values */
    if (IS_ALPHA (ch)) {
    } else if (ch == '-' || ch == '\'' || ch == ' ') {
    } else if (ch == ',' && allowcomma) {
    } else if (ch == '.' && (allowperiod || stp == ptr)) {
    } else {
      /* bad character found */
      *badauthor = str;
      return TRUE;
    }
    ptr++;
    ch = *ptr;
  }

  return FALSE;
}

static Boolean BadCharsInName (ValNodePtr name, CharPtr PNTR badauthor, BoolPtr last_name_badP)

{
  AuthorPtr    ap;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  if (name == NULL) return FALSE;
  ap = name->data.ptrvalue;
  if (ap == NULL) return FALSE;
  pid = ap->name;
  if (pid == NULL) return FALSE;

  if (pid->choice == 2) {
    nsp = pid->data;
    if (nsp == NULL) return FALSE;
    if (StringICmp (nsp->names [0], "et al.") == 0) return FALSE;
    if (BadCharsInAuth (nsp->names [0], badauthor, FALSE, FALSE, TRUE)) {
      if (last_name_badP != NULL) {
        *last_name_badP = TRUE;
      }
      return TRUE; /* last    */
    }
    if (BadCharsInAuth (nsp->names [1], badauthor, FALSE, FALSE, FALSE)) return TRUE; /* first    */
    if (BadCharsInAuth (nsp->names [4], badauthor, FALSE, TRUE, FALSE)) return TRUE;  /* initials */
    if (BadCharsInAuth (nsp->names [5], badauthor, FALSE, TRUE, FALSE)) return TRUE;  /* suffix */
  }

  return FALSE;
}

static CharPtr suffixList [] = {
  "Jr.", "Sr.", "II", "III", "IV", "V", "VI", NULL
};

static void ValidateSuffix (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp, ValNodePtr name)

{
  AuthorPtr    ap;
  Int2         i;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      suffix;

  if (name == NULL) return;
  ap = name->data.ptrvalue;
  if (ap == NULL) return;
  pid = ap->name;
  if (pid == NULL) return;

  if (pid->choice == 2) {
    nsp = pid->data;
    if (nsp == NULL) return;
    suffix = nsp->names [5];
    if (StringHasNoText (suffix)) return;
    for (i = 0; suffixList [i] != NULL; i++) {
      if (StringICmp (suffix, suffixList [i]) == 0) return;
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAuthorSuffix, "Bad author suffix %s", suffix);
  }
}

#define DATE_OKAY       0
#define EMPTY_DATE      1
#define BAD_DATE_STR    2
#define BAD_DATE_YEAR   3
#define BAD_DATE_MONTH  4
#define BAD_DATE_DAY    5
#define BAD_DATE_OTHER  6

static Int2 DateIsBad (DatePtr dp, Boolean needFullDate)

{
  if (dp == NULL) return EMPTY_DATE;

  if (dp->data [0] == 0) {
    if (dp->str == NULL) return BAD_DATE_STR;
    if (StringCmp (dp->str, "?") == 0) return BAD_DATE_STR;
    return DATE_OKAY;
  }

  if (dp->data [0] == 1) {
    if (dp->data [1] == 0) return BAD_DATE_YEAR;
    if (dp->data [2] > 12) return BAD_DATE_MONTH;
    if (dp->data [3] > 31) return BAD_DATE_DAY;
    if (needFullDate) {
      if (dp->data [2] == 0) return BAD_DATE_MONTH;
      if (dp->data [3] == 0) return BAD_DATE_DAY;
    }
    return DATE_OKAY;
  }

  return BAD_DATE_OTHER;
}

static Boolean StringAlreadyInList (
  ValNodePtr head,
  CharPtr str
)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, str) == 0) return TRUE;
  }

  return FALSE;
}

static void ValidatePubdesc (ValidStructPtr vsp, GatherContextPtr gcp, PubdescPtr pdp)
{
  AuthListPtr     alp;
  AuthorPtr       ap;
  CharPtr         badauthor;
  CitArtPtr       cap = NULL;
  CitGenPtr       cgp;
  CitJourPtr      cjp = NULL;
  ValNodePtr      conslist = NULL;
  CitSubPtr       csp;
  DatePtr         dp;
  Boolean         hasName, hasTitle, hasIsoJTA = FALSE,
                  inPress = FALSE, electronic_journal = FALSE,
                  conflicting_pmids = FALSE, redundant_pmids = FALSE,
                  conflicting_muids = FALSE, redundant_muids = FALSE,
                  unpub = FALSE;
  ImprintPtr      imp;
  Boolean         last_name_bad;
  Int4            muid = 0;
  Boolean         noVol, noPages;
  ValNodePtr      name;
  PersonIdPtr     pid;
  Int4            pmid = 0;
  CharPtr         ptr;
  ErrSev          sev;
  Int4            start;
  Int4            stop;
  CharPtr         str;
  Char            temp [64];
  ValNodePtr      title;
  Int4            uid = 0;
  long int        val;
  ValNodePtr      vnp;

  if (vsp == NULL || pdp == NULL)
    return;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
    case PUB_Gen:
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      hasName = FALSE;
      if (cgp != NULL) {
        if (!StringHasNoText (cgp->cit)) {
          if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
              StringNICmp (cgp->cit, "unpublished", 11) == 0 ||
              StringNICmp (cgp->cit, "Online Publication", 18) == 0 ||
              StringNICmp (cgp->cit, "Published Only in DataBase", 26) == 0) {
            unpub = TRUE;
          } else if (StringNICmp (cgp->cit, "(er) ", 5) == 0) {
            unpub = TRUE;
          } else {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Unpublished citation text invalid");
          }
          if (StringStr (cgp->cit, "Title=") != NULL) {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_StructuredCitGenCit, "Unpublished citation has embedded Title");
          }
          if (StringStr (cgp->cit, "Journal=") != NULL) {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_StructuredCitGenCit, "Unpublished citation has embedded Journal");
          }
        }
        /* skip if just serial number */
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) break;
        dp = cgp->date;
        if (dp == NULL) {
          if (! unpub) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date missing");
          }
        } else if (dp->str != NULL) {
          if (StringCmp (dp->str, "?") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date marked as '?'");
          }
        } else if (dp->data [1] == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date not set");
        } else if (DateIsBad (dp, FALSE) > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_GENERIC_BadDate, "Publication date has error");
        }
        alp = cgp->authors;
        if (alp != NULL) {
          if (alp->choice == 1) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoName (name)) {
                hasName = TRUE;
              }
            }
          } else if (alp->choice == 2 || alp->choice == 3) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
                hasName = TRUE;
              }
            }
          }
        }
        if (!hasName) {
          sev = SEV_ERROR;
          if (vsp->is_refseq_in_sep) {
            sev = SEV_WARNING;
          }
          ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Publication has no author names");
        }
      }
      break;
    case PUB_Muid:
      if (pmid == 0) {
        pmid = vnp->data.intvalue;
      } else if (pmid != vnp->data.intvalue) {
        conflicting_pmids = TRUE;
      } else {
        redundant_pmids = TRUE;
      }
      if (uid == 0) {
        uid = vnp->data.intvalue;
      }
      break;
    case PUB_PMid:
      if (muid == 0) {
        muid = vnp->data.intvalue;
      } else if (muid != vnp->data.intvalue) {
        conflicting_muids = TRUE;
      } else {
        redundant_muids = TRUE;
      }
      if (uid == 0) {
        uid = vnp->data.intvalue;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      if (csp != NULL) {
        ValidateCitSub (vsp, csp);
      }
      break;
    case PUB_Medline:
      ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MedlineEntryPub, "Publication is medline entry");
      break;
    case PUB_Article:
      cap = (CitArtPtr) vnp->data.ptrvalue;
      hasName = FALSE;
      hasTitle = FALSE;
      if (cap != NULL) {
        for (title = cap->title; title != NULL; title = title->next) {
          if (!HasNoText ((CharPtr) title->data.ptrvalue)) {
            hasTitle = TRUE;
          }
        }
        if (!hasTitle) {
          ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Publication has no title");
        }
        alp = cap->authors;
        if (alp != NULL) {
          if (alp->choice == 1) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoName (name)) {
                hasName = TRUE;
              }
            }
          } else if (alp->choice == 2 || alp->choice == 3) {
            for (name = alp->names; name != NULL; name = name->next) {
              if (!HasNoText ((CharPtr) name->data.ptrvalue)) {
                hasName = TRUE;
              }
            }
          }
        }
        if (!hasName) {
          ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Publication has no author names");
        }
      }

      switch (cap->from) {
      case 1:
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          hasTitle = FALSE;
          for (title = cjp->title; title != NULL; title = title->next) {
            if (title->choice == Cit_title_iso_jta) {
              hasIsoJTA = TRUE;
            }
            if (!HasNoText ((CharPtr) title->data.ptrvalue)) {
              hasTitle = TRUE;
              if (title->choice == Cit_title_name) {
                if (StringNCmp ((CharPtr) title->data.ptrvalue, "(er)", 4) == 0) {
                  electronic_journal = TRUE;
                }
              }
            }
          }
          if (!hasTitle) {
            ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo, "Journal title missing");
          }
          imp = cjp->imp;
          if (imp != NULL) {
            if (imp->pubstatus == PUBSTATUS_epublish || imp->pubstatus == PUBSTATUS_aheadofprint) {
              electronic_journal = TRUE;
            }
            if (imp->prepub == 2) {
              inPress = TRUE;
            }
            if (imp->prepub == 0 && imp->pubstatus != PUBSTATUS_aheadofprint) {
              noVol = StringHasNoText (imp->volume);
              noPages = StringHasNoText (imp->pages);
              sev = SEV_ERROR;
              if (vsp->is_refseq_in_sep) {
                sev = SEV_WARNING;
              }
              if (noVol && noPages) {
                ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal volume and pages missing");
              } else if (noVol) {
                if (! electronic_journal) {
                  ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal volume missing");
                }
              } else if (noPages) {
                ValidErr (vsp, sev, ERR_GENERIC_MissingPubInfo, "Journal pages missing");
              }
              if (! noPages) {
                sev = SEV_WARNING;
                StringNCpy_0 (temp, imp->pages, sizeof (temp));
                ptr = StringChr (temp, '-');
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                  if (sscanf (temp, "%ld", &val) == 1) {
                    start = (Int4) val;
                    if (sscanf (ptr, "%ld", &val) == 1) {
                      stop = (Int4) val;
                      if (start == 0 || stop == 0) {
                        ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering has zero value");
                      } else if (start < 0 || stop < 0) {
                        ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering has negative value");
                      } else if (start > stop) {
                        ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering out of order");
                      } else if (stop > start + 50) {
                        ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering greater than 50");
                      }
                    } else {
                      ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering stop looks strange");
                    }
                  } else if (! IS_ALPHA (temp [0])) {
                    ValidErr (vsp, sev, ERR_GENERIC_BadPageNumbering, "Page numbering start looks strange");
                  }
                }
              }
              dp = imp->date;
              if (dp == NULL) {
                ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date missing");
              } else if (dp->str != NULL) {
                if (StringCmp (dp->str, "?") == 0) {
                  ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date marked as '?'");
                }
              } else if (dp->data [1] == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "Publication date not set");
              } else if (DateIsBad (dp, FALSE) > 0) {
                ValidErr (vsp, SEV_ERROR, ERR_GENERIC_BadDate, "Publication date has error");
              }
            }
            if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub != 2) {
              if (noVol || noPages) {
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Ahead-of-print without in-press");
              }
            }
            if (imp->pubstatus == PUBSTATUS_epublish && imp->prepub == 2) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Electronic-only publication should not also be in-press");
            }
          }
        }
        break;
      default:
        break;
      }
      break;
    case PUB_Equiv:
      ValidErr (vsp, SEV_WARNING, ERR_GENERIC_UnnecessaryPubEquiv, "Publication has unexpected internal Pub-equiv");
      break;
    default:
      break;
    }
  }

  if (conflicting_pmids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple conflicting pmids in a single publication");
  } else if (redundant_pmids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple redundant pmids in a single publication");
  }
  if (conflicting_muids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple conflicting muids in a single publication");
  } else if (redundant_muids) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_CollidingPublications, "Multiple redundant muids in a single publication");
  }

  if (cap != NULL && cjp != NULL && (uid > 0 || inPress || vsp->alwaysRequireIsoJTA)) {
    if (! hasIsoJTA) {
      if (! electronic_journal) {
        ValidErr (vsp, SEV_WARNING, ERR_GENERIC_MissingPubInfo, "ISO journal title abbreviation missing");
      }
    }
  }

  alp = GetAuthListPtr (pdp, NULL);
  if (alp != NULL) {
    sev = SEV_ERROR;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_WARNING;
    }
    if (alp->choice == 1) {
      for (name = alp->names; name != NULL; name = name->next) {
        badauthor = NULL;
        last_name_bad = FALSE;
        if (BadCharsInName (name, &badauthor, &last_name_bad)) {
          if (StringHasNoText (badauthor)) {
            badauthor = "?";
          }
          if (last_name_bad) {
            ValidErr (vsp, sev, ERR_SEQ_FEAT_BadCharInAuthorLastName, "Bad characters in author %s", badauthor);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadCharInAuthorName, "Bad characters in author %s", badauthor);
          }
        }
        ValidateSuffix (vsp, gcp, pdp, name);
        ap = (AuthorPtr) name->data.ptrvalue;
        if (ap == NULL) continue;
        pid = ap->name;
        if (pid == NULL) continue;
        if (pid->choice == 5) {
          str = (CharPtr) pid->data;
          if (StringHasNoText (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Empty consortium");
            continue;
          }
          if (StringAlreadyInList (conslist, str)) {
            ValidErr (vsp, SEV_WARNING, ERR_GENERIC_PublicationInconsistency, "Duplicate consortium '%s'", str);
            continue;
          }
          ValNodeAddPointer (&conslist, 0, (Pointer) str);
        }
      }
    } else if (alp->choice == 2 || alp->choice == 3) {
      for (name = alp->names; name != NULL; name = name->next) {
        badauthor = NULL;
        if (BadCharsInAuth ((CharPtr) name->data.ptrvalue, &badauthor, TRUE, TRUE, FALSE)) {
          if (StringHasNoText (badauthor)) {
            badauthor = "?";
          }
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadCharInAuthorName, "Bad characters in author %s", badauthor);
        }
      }
    }
  }
  ValNodeFree (conslist);
}

static void ValidateSfpCit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)
{
  ValNodePtr      ppr;
  ValNodePtr      psp;

  if (vsp == NULL || sfp == NULL || sfp->cit == NULL)
    return;
  psp = sfp->cit;
  if (psp == NULL)
    return;
  for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
    if (ppr->choice == PUB_Equiv) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryCitPubEquiv, "Citation on feature has unexpected internal Pub-equiv");
      return;
    }
  }
}

typedef struct bioseqvalid
{
  ValidStructPtr  vsp;
  Boolean         is_aa;         /* bioseq is protein? */
  Boolean         is_mrna;       /* molinfo is mrna? */
  Boolean         is_prerna;     /* molinfo is precursor rna? */
  Boolean         is_artificial; /* biosource origin is artificial, synthetic, mutant? */
  Boolean         is_syn_constr; /* is organism name synthetic construct, plasmid, vector, or SYN division? */
  Boolean         got_a_pub;
  int             last_na_mol, last_na_mod, last_organelle, last_partialness, last_left_right,
                  last_biomol, last_tech, last_completeness,
                  num_full_length_src_feat, num_full_length_prot_ref,
                  num_justprot, num_preprot, num_matpep, num_sigpep, num_transpep;
  ValNodePtr      last_gb, last_embl, last_prf, last_pir, last_sp, last_pdb,
                  last_create, last_update, last_biosrc, last_orgref;
  OrgRefPtr       last_org;
  GatherContextPtr gcp;
  BioseqPtr        bsp;
}
BioseqValidStr , PNTR BioseqValidStrPtr;

static void CheckForNucProt (BioseqSetPtr bssp, Pointer userdata)
{
  BoolPtr         hasPartsP;

  if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    hasPartsP = (BoolPtr) userdata;
    *hasPartsP = TRUE;
  }
}

static void CheckForParts (BioseqSetPtr bssp, Pointer userdata)
{
  BoolPtr         hasPartsP;

  if (bssp->_class == BioseqseqSet_class_parts) {
    hasPartsP = (BoolPtr) userdata;
    *hasPartsP = TRUE;
  }
}

static Boolean DeltaOrFarSeg (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  Boolean         hasParts = FALSE;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    if (bsp->repr == Seq_repr_delta) {
      VisitSetsInSep (sep, (Pointer) &hasParts, CheckForNucProt);
      if (!hasParts)
        return TRUE;
    }
    if (bsp->repr == Seq_repr_seg) {
      VisitSetsInSep (sep, (Pointer) &hasParts, CheckForParts);
      if (!hasParts)
        return TRUE;
    }
  }
  return FALSE;
}

static void 
ValidateIntronEndsAtSpliceSiteOrGap 
(ValidStructPtr vsp, 
 SeqLocPtr slp)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  Uint1              strand;
  Int4               strt, stop, pos;
  Boolean            partial5, partial3;
  Char               buf[3];
  Char               id_buf[150];
  SeqFeatPtr         rna;
  SeqMgrFeatContext  rcontext;

  if (vsp == NULL || slp == NULL) return;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  if (partial5 && partial3) return;

  /* suppress if contained by rRNA - different consensus splice site */

  rna = SeqMgrGetOverlappingFeature (slp, 0, vsp->rrna_array, vsp->numrrna,
                                     NULL, CONTAINED_WITHIN, &rcontext);
  if (rna != NULL) return;

  /* suppress if contained by tRNA - different consensus splice site */

  rna = SeqMgrGetOverlappingFeature (slp, 0, vsp->trna_array, vsp->numtrna,
                                     NULL, CONTAINED_WITHIN, &rcontext);
  if (rna != NULL) return;


  sip = SeqLocId (slp);
  if (sip == NULL)
    return;
  
  bsp = NULL;
  if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
    bsp = BioseqLockById (sip);
  }
  if (bsp == NULL)
    return;

  BioseqLabel (bsp, id_buf, sizeof (id_buf) - 1, OM_LABEL_CONTENT);

  strt = SeqLocStart (slp);
  stop = SeqLocStop (slp);

  strand = SeqLocStrand (slp);

  if (!partial5) {
    if (strand == Seq_strand_minus) {
      SeqPortStreamInt (bsp, stop - 1, stop, Seq_strand_minus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = stop;
    } else {
      SeqPortStreamInt (bsp, strt, strt + 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = strt;
    }
    if ((buf[0] == '-' && buf[1] == '-')
        || (buf[0] == 'G' && buf[1] == 'T')
        || (buf[0] == 'G' && buf[1] == 'C')) {
      /* location is ok */
    } else if (pos == 0 || pos == bsp->length - 1) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                "Splice donor consensus (GT) not found at start of terminal intron, position %ld of %s", (long) (pos + 1), id_buf);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                "Splice donor consensus (GT) not found at start of intron, position %ld of %s", (long) (pos + 1), id_buf);
    }
  }
  if (!partial3) {
    if (strand == Seq_strand_minus) {
      SeqPortStreamInt (bsp, strt, strt + 1, Seq_strand_minus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = strt;
    } else {
      SeqPortStreamInt (bsp, stop - 1, stop, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      pos = stop;
    }
    if ((buf[0] == '-' && buf[1] == '-')
        || (buf[0] == 'A' && buf[1] == 'G')) {
      /* location is ok */
    } else if (pos == 0 || pos == bsp->length - 1) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                "Splice acceptor consensus (AG) not found at end of terminal intron, position %ld of %s, but at end of sequence", (long) (pos + 1), id_buf);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                "Splice acceptor consensus (AG) not found at end of intron, position %ld of %s", (long) (pos + 1), id_buf);
    }
  }
  BioseqUnlock (bsp);  
}


/*****************************************************************************
*
*   ValidateSeqFeatContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static Boolean ValidateSeqFeatCommon (SeqFeatPtr sfp, BioseqValidStrPtr bvsp, ValidStructPtr vsp,
                                      Int4 left, Int4 right, Int2 numivals, Uint4 featitemid, Boolean farloc, BioseqPtr bsp)
{
  BioseqSetPtr    bssp;
  Boolean         do_error;
  GatherContextPtr gcp = NULL;
  ImpFeatPtr      ifp;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  ProtRefPtr      prp;
  RnaRefPtr       rrp;
  CharPtr         str;
  SeqLocPtr       slp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  Boolean         on_seg = FALSE;
  Boolean         is_emb = FALSE;
  Boolean         is_nc = FALSE;
  Boolean         is_refseq = FALSE;
  ErrSev          sev;
  Boolean         no_nonconsensus_except;


  vsp->descr = NULL;
  vsp->sfp = sfp;

  if (featitemid > 0) {
    gcp = vsp->gcp;
    if (gcp != NULL) {
      olditemid = gcp->itemID;
      olditemtype = gcp->thistype;
      gcp->itemID = featitemid;
      gcp->thistype = OBJ_SEQFEAT;
    }
  }

  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        is_refseq = TRUE;
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            is_nc = TRUE;
          }
        }
      } else if (sip->choice == SEQID_EMBL) {
        is_emb = TRUE;
      }
    }
  }

  if (bvsp->is_aa) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      if ((left == 0) && (right == ((vsp->bsp->length) - 1))) {
        bvsp->num_full_length_prot_ref++;
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp != NULL) {
          switch (prp->processed) {
            case 0:
              bvsp->num_justprot++;
              break;
            case 1:
              bvsp->num_preprot++;
              break;
            case 2:
              bvsp->num_matpep++;
              break;
            case 3:
              bvsp->num_sigpep++;
              break;
            case 4:
              bvsp->num_transpep++;
              break;
            default:
              break;
          }
        }
      }
    }

    switch (sfp->data.choice) {
    case SEQFEAT_CDREGION:
    case SEQFEAT_RNA:
    case SEQFEAT_RSITE:
    case SEQFEAT_TXINIT:
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a protein Bioseq.");
      break;
    case SEQFEAT_GENE:
        if (bsp != NULL) {
          do_error = FALSE;
          if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
            bssp = (BioseqSetPtr) bsp->idx.parentptr;
            while (bssp != NULL) {
              switch (bssp->_class) {
              case BioseqseqSet_class_nuc_prot :
              case BioseqseqSet_class_mut_set :
              case BioseqseqSet_class_pop_set :
              case BioseqseqSet_class_phy_set :
              case BioseqseqSet_class_eco_set :
              case BioseqseqSet_class_gen_prod_set :
                do_error = TRUE;
                break;
              default :
                break;
              }
              if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
                bssp = (BioseqSetPtr) bssp->idx.parentptr;
              } else {
                bssp = NULL;
              }
            }
          }
          if (do_error) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a protein Bioseq.");
          }
        }
      break;
    default:
      break;
    }

  } else {
    switch (sfp->data.choice) {
    case SEQFEAT_PROT:
    case SEQFEAT_PSEC_STR:
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a nucleotide Bioseq.");
      break;
    default:
      break;
    }

  }

  if (bvsp->is_mrna) {
    switch (sfp->data.choice) {
    case SEQFEAT_CDREGION:
      if (numivals > 1) {
        if ((! sfp->excpt) ||
            (StringISearch (sfp->except_text, "ribosomal slippage") == NULL)) {
          sev = SEV_ERROR;
          if (is_refseq) {
            sev = SEV_WARNING;
          }
          ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Multi-interval CDS feature is invalid on an mRNA (cDNA) Bioseq.");
        }
      }
      break;
    case SEQFEAT_RNA:
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 2) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "mRNA feature is invalid on an mRNA (cDNA) Bioseq.");
      }
      break;
    case SEQFEAT_IMP:
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL && ifp->key != NULL && (!HasNoText (ifp->key))) {
        if (StringCmp (ifp->key, "intron") == 0 || StringCmp (ifp->key, "CAAT_signal") == 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an mRNA Bioseq.");
        }
      }
      break;
    default:
      break;
    }
  } else if (bvsp->is_prerna) {
    switch (sfp->data.choice) {
    case SEQFEAT_IMP:
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL && ifp->key != NULL && (!HasNoText (ifp->key))) {
        if (StringCmp (ifp->key, "CAAT_signal") == 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an pre-RNA Bioseq.");
        }
      }
      break;
    default:
      break;
    }
  }

  if (farloc && (! is_nc) && (! is_emb)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FarLocation, "Feature has 'far' location - accession not packaged in record");
  }

  if ((sfp->data.choice == SEQFEAT_PUB) || (sfp->cit != NULL))
    bvsp->got_a_pub = TRUE;

  str = (CharPtr) sfp->comment;
  if (SerialNumberInString (str)) {
    ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_SerialInComment,
              "Feature comment may refer to reference by serial number - attach reference specific comments to the reference REMARK instead.");
  }

  if (bsp != NULL && bsp->repr == Seq_repr_seg) {
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        if (SeqIdIn (sip, bsp->id)) {
          on_seg = TRUE;
        }
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (on_seg) {
      sev = SEV_ERROR;
      if (is_nc) {
        sev = SEV_WARNING;
      }
      if (! DeltaOrFarSeg (vsp->sep, sfp->location)) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_LocOnSegmentedBioseq, "Feature location on segmented bioseq, not on parts");
      }
    }
  }

  if (sfp->idx.subtype == FEATDEF_intron) {
    no_nonconsensus_except = TRUE;
    if (sfp->excpt) {
      if (StringISearch (sfp->except_text, "nonconsensus splice site") != NULL) {
        no_nonconsensus_except = FALSE;
      }
    }
    if (no_nonconsensus_except) {
      ValidateIntronEndsAtSpliceSiteOrGap (vsp, sfp->location);
    }
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  return TRUE;
}

static Boolean GeneSpansOrigin (SeqMgrFeatContextPtr context, Int4 bsplength)

{
  Int4Ptr  ivals;

  if (context == NULL || bsplength < 1) return FALSE;
  ivals = context->ivals;
  if (ivals == NULL || context->numivals != 2) return FALSE;
  if (context->strand == Seq_strand_minus) {
    if (ivals [1] == 0 && ivals [2] == bsplength - 1) return TRUE;
  } else {
    if (ivals [2] == 0 && ivals [1] == bsplength - 1) return TRUE;
  }
  return FALSE;
}

static void CheckMultiIntervalGene (SeqFeatPtr sfp, SeqMgrFeatContextPtr context, ValidStructPtr vsp, GatherContextPtr gcp)

{
  BioseqPtr     bsp;
  Int4          count;
  SeqLocPtr     mappedloc = NULL;
  Uint2         olditemtype = 0;
  Uint4         olditemid = 0;
  Boolean       segmented = FALSE;
  ErrSev        sev = SEV_ERROR;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  TextSeqIdPtr  tsip;

  if (sfp == NULL || context == NULL || vsp == NULL) return;
  if (context->numivals < 2) return;

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) return;
  }

  if (SeqLocId (sfp->location) == NULL) {
    bsp = context->bsp;
    if (bsp == NULL || bsp->repr != Seq_repr_seg) return;
    mappedloc = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
    if (mappedloc == NULL) return;
    count = 0;
    slp = SeqLocFindNext (mappedloc, NULL);
    while (slp != NULL) {
      count++;
      slp = SeqLocFindNext (mappedloc, slp);
    }
    SeqLocFree (mappedloc);
    if (count < 2) return;
    segmented = TRUE;
  }

  bsp = context->bsp;
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            sev = SEV_WARNING;
          }
        }
      } else if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {
        sev = SEV_WARNING;
      }
    }
    if (bsp->topology == 2) {
      if (context->numivals == 2 && GeneSpansOrigin (context, bsp->length)) return;
      sev = SEV_WARNING;
    }
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
    gcp->itemID = context->itemID;
    gcp->thistype = OBJ_SEQFEAT;
  }

  vsp->sfp = sfp;
  if (segmented) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_SegmentedGeneProblem,
              "Gene feature on segmented sequence should cover all bases within its extremes");
  } else {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_MultiIntervalGene,
              "Gene feature on non-segmented sequence should not have multiple intervals");
  }
  vsp->sfp = NULL;

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean LIBCALLBACK ValidateSeqFeatIndexed (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;

  bvsp = (BioseqValidStrPtr) context->userdata;
  vsp = bvsp->vsp;

  if (sfp->data.choice == SEQFEAT_GENE) {
    CheckMultiIntervalGene (sfp, context, vsp, vsp->gcp);
  }

  return ValidateSeqFeatCommon (sfp, bvsp, vsp, context->left, context->right, context->numivals, context->itemID, context->farloc, context->bsp);
}

static void ValidateSeqFeatContext (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  SeqFeatPtr      sfp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  sfp = (SeqFeatPtr) (gcp->thisitem);

  ValidateSeqFeatCommon (sfp, bvsp, vsp, gcp->extremes.left, gcp->extremes.right, 0, 0, FALSE, NULL);
}

/*****************************************************************************
*
*   CountryIsValid(name)
*      Validates subsource country against official country names
*
*****************************************************************************/

static CharPtr  Nlm_valid_country_codes [] = {
  "Afghanistan",
  "Albania",
  "Algeria",
  "American Samoa",
  "Andorra",
  "Angola",
  "Anguilla",
  "Antarctica",
  "Antigua and Barbuda",
  "Arctic Ocean",
  "Argentina",
  "Armenia",
  "Aruba",
  "Ashmore and Cartier Islands",
  "Atlantic Ocean",
  "Australia",
  "Austria",
  "Azerbaijan",
  "Bahamas",
  "Bahrain",
  "Baker Island",
  "Bangladesh",
  "Barbados",
  "Bassas da India",
  "Belarus",
  "Belgium",
  "Belize",
  "Benin",
  "Bermuda",
  "Bhutan",
  "Bolivia",
  "Bosnia and Herzegovina",
  "Botswana",
  "Bouvet Island",
  "Brazil",
  "British Virgin Islands",
  "Brunei",
  "Bulgaria",
  "Burkina Faso",
  "Burundi",
  "Cambodia",
  "Cameroon",
  "Canada",
  "Cape Verde",
  "Cayman Islands",
  "Central African Republic",
  "Chad",
  "Chile",
  "China",
  "Christmas Island",
  "Clipperton Island",
  "Cocos Islands",
  "Colombia",
  "Comoros",
  "Cook Islands",
  "Coral Sea Islands",
  "Costa Rica",
  "Cote d'Ivoire",
  "Croatia",
  "Cuba",
  "Cyprus",
  "Czech Republic",
  "Democratic Republic of the Congo",
  "Denmark",
  "Djibouti",
  "Dominica",
  "Dominican Republic",
  "East Timor",
  "Ecuador",
  "Egypt",
  "El Salvador",
  "Equatorial Guinea",
  "Eritrea",
  "Estonia",
  "Ethiopia",
  "Europa Island",
  "Falkland Islands (Islas Malvinas)",
  "Faroe Islands",
  "Fiji",
  "Finland",
  "France",
  "French Guiana",
  "French Polynesia",
  "French Southern and Antarctic Lands",
  "Gabon",
  "Gambia",
  "Gaza Strip",
  "Georgia",
  "Germany",
  "Ghana",
  "Gibraltar",
  "Glorioso Islands",
  "Greece",
  "Greenland",
  "Grenada",
  "Guadeloupe",
  "Guam",
  "Guatemala",
  "Guernsey",
  "Guinea",
  "Guinea-Bissau",
  "Guyana",
  "Haiti",
  "Heard Island and McDonald Islands",
  "Honduras",
  "Hong Kong",
  "Howland Island",
  "Hungary",
  "Iceland",
  "India",
  "Indian Ocean",
  "Indonesia",
  "Iran",
  "Iraq",
  "Ireland",
  "Isle of Man",
  "Israel",
  "Italy",
  "Jamaica",
  "Jan Mayen",
  "Japan",
  "Jarvis Island",
  "Jersey",
  "Johnston Atoll",
  "Jordan",
  "Juan de Nova Island",
  "Kazakhstan",
  "Kenya",
  "Kerguelen Archipelago",
  "Kingman Reef",
  "Kiribati",
  "Kosovo",
  "Kuwait",
  "Kyrgyzstan",
  "Laos",
  "Latvia",
  "Lebanon",
  "Lesotho",
  "Liberia",
  "Libya",
  "Liechtenstein",
  "Lithuania",
  "Luxembourg",
  "Macau",
  "Macedonia",
  "Madagascar",
  "Malawi",
  "Malaysia",
  "Maldives",
  "Mali",
  "Malta",
  "Marshall Islands",
  "Martinique",
  "Mauritania",
  "Mauritius",
  "Mayotte",
  "Mexico",
  "Micronesia",
  "Midway Islands",
  "Moldova",
  "Monaco",
  "Mongolia",
  "Montenegro",
  "Montserrat",
  "Morocco",
  "Mozambique",
  "Myanmar",
  "Namibia",
  "Nauru",
  "Navassa Island",
  "Nepal",
  "Netherlands",
  "Netherlands Antilles",
  "New Caledonia",
  "New Zealand",
  "Nicaragua",
  "Niger",
  "Nigeria",
  "Niue",
  "Norfolk Island",
  "North Korea",
  "Northern Mariana Islands",
  "Norway",
  "Oman",
  "Pacific Ocean",
  "Pakistan",
  "Palau",
  "Palmyra Atoll",
  "Panama",
  "Papua New Guinea",
  "Paracel Islands",
  "Paraguay",
  "Peru",
  "Philippines",
  "Pitcairn Islands",
  "Poland",
  "Portugal",
  "Puerto Rico",
  "Qatar",
  "Republic of the Congo",
  "Reunion",
  "Romania",
  "Russia",
  "Rwanda",
  "Saint Helena",
  "Saint Kitts and Nevis",
  "Saint Lucia",
  "Saint Pierre and Miquelon",
  "Saint Vincent and the Grenadines",
  "Samoa",
  "San Marino",
  "Sao Tome and Principe",
  "Saudi Arabia",
  "Senegal",
  "Serbia",
  "Seychelles",
  "Sierra Leone",
  "Singapore",
  "Slovakia",
  "Slovenia",
  "Solomon Islands",
  "Somalia",
  "South Africa",
  "South Georgia and the South Sandwich Islands",
  "South Korea",
  "Spain",
  "Spratly Islands",
  "Sri Lanka",
  "Sudan",
  "Suriname",
  "Svalbard",
  "Swaziland",
  "Sweden",
  "Switzerland",
  "Syria",
  "Taiwan",
  "Tajikistan",
  "Tanzania",
  "Thailand",
  "Togo",
  "Tokelau",
  "Tonga",
  "Trinidad and Tobago",
  "Tromelin Island",
  "Tunisia",
  "Turkey",
  "Turkmenistan",
  "Turks and Caicos Islands",
  "Tuvalu",
  "Uganda",
  "Ukraine",
  "United Arab Emirates",
  "United Kingdom",
  "Uruguay",
  "USA",
  "Uzbekistan",
  "Vanuatu",
  "Venezuela",
  "Viet Nam",
  "Virgin Islands",
  "Wake Island",
  "Wallis and Futuna",
  "West Bank",
  "Western Sahara",
  "Yemen",
  "Zambia",
  "Zimbabwe",
  NULL
};

static CharPtr  Nlm_formerly_valid_country_codes [] = {
  "Belgian Congo",
  "British Guiana",
  "Burma",
  "Czechoslovakia",
  "Serbia and Montenegro",
  "Siam",
  "USSR",
  "Yugoslavia",
  "Zaire",
  NULL
};

NLM_EXTERN CharPtr PNTR GetValidCountryList (void)

{
  return (CharPtr PNTR) Nlm_valid_country_codes;
}

NLM_EXTERN Boolean CountryIsValid (CharPtr name, BoolPtr old_countryP, BoolPtr bad_capP)
{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return FALSE;

  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ':');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (Nlm_valid_country_codes) / sizeof (Nlm_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_valid_country_codes [R], str) == 0) {
    if (bad_capP != NULL) {
      if (StringCmp (Nlm_valid_country_codes [R], str) != 0) {
        *bad_capP = TRUE;
      }
    }
    return TRUE;
  }

  L = 0;
  R = sizeof (Nlm_formerly_valid_country_codes) / sizeof (Nlm_formerly_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_formerly_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_formerly_valid_country_codes [R], str) == 0) {
    if (old_countryP != NULL) {
      *old_countryP = TRUE;
    }
    if (bad_capP != NULL) {
      if (StringCmp (Nlm_formerly_valid_country_codes [R], str) != 0) {
        *bad_capP = TRUE;
      }
    }
    return FALSE;
  }

  return FALSE;
}


NLM_EXTERN CharPtr GetCorrectedCountryCapitalization (CharPtr name)
{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return NULL;

  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ':');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (Nlm_valid_country_codes) / sizeof (Nlm_valid_country_codes [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_valid_country_codes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (Nlm_valid_country_codes [R], str) == 0) {
    return Nlm_valid_country_codes[R];
  }

  return NULL;
}


static CharPtr ctry_lat_lon [] = {
  "Afghanistan\tAF\t60.4\t29.3\t74.9\t38.5",
  "Albania\tAL\t19.2\t39.6\t21.1\t42.7",
  "Algeria\tAG\t-8.7\t18.9\t12.0\t37.1",
  "American Samoa\tAQ\t-171.1\t-11.1\t-171.1\t-11.0\t-170.9\t-14.4\t-169.4\t-14.2",
  "Andorra\tAN\t1.4\t42.4\t1.8\t42.7",
  "Angola\tAO\t11.6\t-18.1\t24.1\t-4.4",
  "Anguilla\tAV\t-63.2\t18.1\t-62.9\t18.3",
  "Antarctica\tAY\t",
  "Antigua and Barbuda\tAC\t-62.4\t16.9\t-62.3\t16.9\t-62.0\t16.9\t-61.7\t17.7",
  "Arctic Ocean\tXX\t",
  "Argentina\tAR\t-73.6\t-55.1\t-53.6\t-21.8",
  "Armenia\tAM\t43.4\t38.8\t46.6\t41.3",
  "Aruba\tAA\t-70.1\t12.4\t-69.8\t12.7",
  "Ashmore and Cartier Islands\tAT\t122.9\t-12.3\t123.1\t-12.1",
  "Atlantic Ocean\tXX\t",
  "Australia\tAS\t112.9\t-43.7\t153.6\t-10.0",
  "Australia: Australian Capital Territory\tXX\t148.7\t-36.0\t149.4\t-35.1",
  "Australia: Jervis Bay Territory\tXX\t150.5\t-35.2\t150.8\t-35.1",
  "Australia: New South Wales\tXX\t140.9\t-37.6\t153.6\t-28.2",
  "Australia: Northern Territory\tXX\t128.9\t-26.1\t138.0\t-10.9",
  "Australia: Queensland\tXX\t137.9\t-29.2\t153.6\t-10.0",
  "Australia: South Australia\tXX\t128.9\t-38.1\t141.0\t-26.0",
  "Australia: Tasmania\tXX\t143.8\t-43.7\t148.5\t-39.6",
  "Australia: Victoria\tXX\t140.9\t-39.6\t150.0\t-34.0",
  "Australia: Western Australia\tXX\t112.9\t-35.2\t129.0\t-13.7",
  "Austria\tAU\t9.5\t46.3\t17.2\t49.0",
  "Azerbaijan\tAJ\t45.0\t38.3\t50.6\t41.9",
  "Bahamas\tBF\t-79.7\t20.9\t-72.7\t27.2",
  "Bahrain\tBA\t50.3\t25.7\t50.7\t26.3",
  "Baker Island\tFQ\t-176.5\t0.1\t-176.5\t0.2",
  "Bangladesh\tBG\t88.0\t20.5\t92.7\t26.6",
  "Barbados\tBB\t-59.7\t13.0\t-59.4\t13.3",
  "Bassas da India\tBS\t39.6\t-21.6\t39.8\t-21.4",
  "Belarus\tBO\t23.1\t51.2\t32.8\t56.2",
  "Belgium\tBE\t2.5\t49.4\t6.4\t51.5",
  "Belize\tBH\t-89.3\t15.8\t-86.9\t18.5",
  "Benin\tBN\t0.7\t6.2\t3.9\t12.4",
  "Bermuda\tBD\t-64.9\t32.2\t-64.7\t32.4",
  "Bhutan\tBT\t88.7\t26.7\t92.1\t28.3",
  "Bolivia\tBL\t-69.7\t-22.9\t-57.5\t-9.7",
  "Bosnia and Herzegovina\tBK\t15.7\t42.5\t19.7\t45.3",
  "Botswana\tBC\t19.9\t-27.0\t29.4\t-17.8",
  "Bouvet Island\tBV\t3.3\t-54.5\t3.5\t-54.4",
  "Brazil\tBR\t-74.0\t-33.8\t-34.8\t5.0",
  "British Virgin Islands\tVI\t-64.8\t18.2\t-63.2\t18.8",
  "Brunei\tBX\t114.0\t4.0\t115.4\t5.0",
  "Bulgaria\tBU\t22.3\t41.2\t28.6\t44.2",
  "Burkina Faso\tUV\t-5.6\t9.4\t2.4\t15.1",
  "Burundi\tBY\t28.9\t-4.5\t30.8\t-2.3",
  "Cambodia\tCB\t102.3\t9.2\t107.6\t14.7",
  "Cameroon\tCM\t8.4\t1.6\t16.2\t13.1",
  "Canada\tCA\t-141.0\t41.7\t-52.6\t83.1",
  "Canada: Alberta\tXX\t-120.0\t48.9\t-110.0\t60.0",
  "Canada: British Columbia\tXX\t-139.1\t48.3\t-114.1\t60.0",
  "Canada: Manitoba\tXX\t-102.1\t48.9\t-89.0\t60.0",
  "Canada: New Brunswick\tXX\t-69.1\t44.5\t-63.8\t48.1",
  "Canada: Newfoundland and Labrador\tXX\t-67.9\t46.6\t-52.6\t60.4",
  "Canada: Northwest Territories\tXX\t-136.5\t60.0\t-102.0\t78.8",
  "Canada: Nova Scotia\tXX\t-66.4\t43.3\t-59.7\t47.0",
  "Canada: Nunavut\tXX\t-120.4\t60.0\t-61.2\t83.1",
  "Canada: Ontario\tXX\t-95.2\t41.6\t-74.3\t56.9",
  "Canada: Prince Edward Island\tXX\t-64.5\t45.9\t-62.0\t47.1",
  "Canada: Quebec\tXX\t-79.8\t45.0\t-57.1\t62.6",
  "Canada: Saskatchewan\tXX\t-110.0\t48.9\t-101.4\t60.0",
  "Canada: Yukon\tXX\t-141.0\t60.0\t-124.0\t69.6",
  "Cape Verde\tCV\t-25.4\t14.8\t-22.7\t17.2",
  "Cayman Islands\tCJ\t-81.5\t19.2\t-81.1\t19.4\t-80.2\t19.6\t-79.7\t19.8",
  "Central African Republic\tCT\t14.4\t2.2\t27.5\t11.0",
  "Chad\tCD\t13.4\t7.4\t24.0\t23.5",
  "Chile\tCI\t-75.8\t-56.0\t-66.4\t-17.5",
  "China\tCH\t73.5\t20.2\t134.8\t53.6\t108.6\t18.1\t111.1\t20.2",
  "China: Hainan\tXX\t108.6\t18.1\t111.1\t20.2",
  "Christmas Island\tKT\t105.5\t-10.6\t105.7\t-10.4",
  "Clipperton Island\tIP\t-109.3\t10.2\t-109.2\t10.3",
  "Cocos Islands\tCK\t96.8\t-12.2\t96.9\t-11.8",
  "Colombia\tCO\t-79.1\t-4.3\t-66.9\t12.5",
  "Comoros\tCN\t43.2\t-12.5\t44.5\t-11.4",
  "Cook Islands\tCW\t-159.9\t-22.0\t-157.3\t-18.8",
  "Coral Sea Islands\tCR\t",
  "Costa Rica\tCS\t-87.1\t5.4\t-87.0\t5.6\t-86.0\t8.0\t-82.6\t11.2",
  "Cote d'Ivoire\tIV\t-8.6\t4.3\t-2.5\t10.7",
  "Croatia\tHR\t13.4\t42.3\t19.4\t46.5",
  "Cuba\tCU\t-85.0\t19.8\t-74.1\t23.3",
  "Cyprus\tCY\t32.2\t34.5\t34.6\t35.7",
  "Czech Republic\tEZ\t12.0\t48.5\t18.9\t51.0",
  "Democratic Republic of the Congo\tCG\t12.2\t-13.5\t31.3\t5.4",
  "Denmark\tDA\t8.0\t54.5\t12.7\t57.7\t14.6\t54.9\t15.2\t55.3",
  "Djibouti\tDJ\t41.7\t10.9\t43.4\t12.7",
  "Dominica\tDO\t-61.5\t15.2\t-61.2\t15.6",
  "Dominican Republic\tDR\t-72.1\t17.4\t-68.3\t19.9",
  "East Timor\tTT\t124.9\t-9.5\t127.4\t-8.3",
  "Ecuador\tEC\t-92.1\t-1.5\t-89.2\t1.7\t-81.1\t-5.0\t-75.2\t1.4",
  "Ecuador: Galapagos\tXX\t-92.1\t-1.5\t-89.2\t1.7",
  "Egypt\tEG\t24.6\t21.7\t35.8\t31.7",
  "El Salvador\tES\t-90.2\t13.1\t-87.7\t14.4",
  "Equatorial Guinea\tEK\t8.4\t3.2\t8.9\t3.8\t9.2\t0.8\t11.3\t2.3",
  "Eritrea\tER\t36.4\t12.3\t43.1\t18.0",
  "Estonia\tEN\t21.7\t57.5\t28.2\t59.7",
  "Ethiopia\tET\t32.9\t3.4\t48.0\t14.9",
  "Europa Island\tEU\t40.3\t-22.4\t40.4\t-22.3",
  "Falkland Islands (Islas Malvinas)\tFK\t-61.4\t-53.0\t-57.7\t-51.0",
  "Faroe Islands\tFO\t-7.7\t61.3\t-6.3\t62.4",
  "Fiji\tFJ\t-180.0\t-20.7\t-178.2\t-15.7\t-175.7\t-19.8\t-175.0\t-15.6\t176.8\t-19.3\t180.0\t-12.5",
  "Finland\tFI\t19.3\t59.7\t31.6\t70.1",
  "France\tFR\t-5.2\t42.3\t8.2\t51.1\t8.5\t41.3\t9.6\t43.1",
  "France: Corsica\tXX\t8.5\t41.3\t9.6\t43.1",
  "French Guiana\tFG\t-54.6\t2.1\t-51.6\t5.8",
  "French Polynesia\tFP\t-154.7\t-27.7\t-134.9\t-7.8",
  "French Southern and Antarctic Lands\tFS\t68.6\t-49.8\t70.6\t-48.5",
  "Gabon\tGB\t8.6\t-4.0\t14.5\t2.3",
  "Gambia\tGA\t-16.9\t13.0\t-13.8\t13.8",
  "Gaza Strip\tGZ\t34.2\t31.2\t34.5\t31.6",
  "Georgia\tGG\t40.0\t41.0\t46.7\t43.6",
  "Germany\tGM\t5.8\t47.2\t15.0\t55.1",
  "Ghana\tGH\t-3.3\t4.7\t1.2\t11.2",
  "Gibraltar\tGI\t-5.4\t36.1\t-5.3\t36.2",
  "Glorioso Islands\tGO\t47.2\t-11.6\t47.4\t-11.5",
  "Greece\tGR\t19.3\t34.8\t28.2\t41.8",
  "Greenland\tGL\t-73.3\t59.7\t-11.3\t83.6",
  "Grenada\tGJ\t-61.8\t11.9\t-61.6\t12.3",
  "Guadeloupe\tGP\t-63.2\t17.8\t-62.8\t18.1\t-61.9\t15.8\t-61.0\t16.5",
  "Guam\tGQ\t144.6\t13.2\t145.0\t13.7",
  "Guatemala\tGT\t-92.3\t13.7\t-88.2\t17.8",
  "Guernsey\tGK\t-2.7\t49.4\t-2.4\t49.5",
  "Guinea\tGV\t-15.1\t7.1\t-7.6\t12.7",
  "Guinea-Bissau\tPU\t-16.8\t10.8\t-13.6\t12.7",
  "Guyana\tGY\t-61.4\t1.1\t-56.5\t8.6",
  "Haiti\tHA\t-74.5\t18.0\t-71.6\t20.1",
  "Heard Island and McDonald Islands\tHM\t73.2\t-53.2\t73.7\t-52.9",
  "Honduras\tHO\t-89.4\t12.9\t-83.2\t16.5",
  "Hong Kong\tHK\t113.8\t22.1\t114.4\t22.6",
  "Howland Island\tHQ\t-176.7\t0.7\t-176.6\t0.8",
  "Hungary\tHU\t16.1\t45.7\t22.9\t48.6",
  "Iceland\tIC\t-24.6\t63.2\t-13.5\t66.6",
  "India\tIN\t67.3\t8.0\t97.4\t35.5",
  "Indian Ocean\tXX\t",
  "Indonesia\tID\t95.0\t-11.1\t141.0\t5.9",
  "Iran\tIR\t44.0\t25.0\t63.3\t39.8",
  "Iraq\tIZ\t38.8\t29.1\t48.6\t37.4",
  "Ireland\tEI\t-10.7\t51.4\t-6.0\t55.4",
  "Isle of Man\tIM\t-4.9\t54.0\t-4.3\t54.4",
  "Israel\tIS\t34.2\t29.4\t35.7\t33.3",
  "Italy\tIT\t6.6\t35.4\t18.5\t47.1",
  "Jamaica\tJM\t-78.4\t17.7\t-76.2\t18.5",
  "Jan Mayen\tJN\t-9.1\t70.8\t-7.9\t71.2",
  "Japan\tJA\t122.9\t24.0\t125.5\t25.9\t126.7\t20.5\t145.8\t45.5",
  "Jarvis Island\tDQ\t-160.1\t-0.4\t-160.0\t-0.4",
  "Jersey\tJE\t-2.3\t49.1\t-2.0\t49.3",
  "Johnston Atoll\tJQ\t-169.6\t16.7\t-169.4\t16.8",
  "Jordan\tJO\t34.9\t29.1\t39.3\t33.4",
  "Juan de Nova Island\tJU\t42.6\t-17.1\t42.8\t-16.8",
  "Kazakhstan\tKZ\t46.4\t40.9\t87.3\t55.4",
  "Kenya\tKE\t33.9\t-4.7\t41.9\t4.6",
  "Kerguelen Archipelago\tXX\t",
  "Kingman Reef\tKQ\t-162.9\t6.1\t-162.4\t6.7",
  "Kiribati\tKR\t172.6\t0.1\t173.9\t3.4\t174.2\t-2.7\t176.9\t-0.5",
  "Kosovo\tKV\t20.0\t41.8\t43.3\t21.9",
  "Kuwait\tKU\t46.5\t28.5\t48.4\t30.1",
  "Kyrgyzstan\tKG\t69.2\t39.1\t80.3\t43.2",
  "Laos\tLA\t100.0\t13.9\t107.7\t22.5",
  "Latvia\tLG\t20.9\t55.6\t28.2\t58.1",
  "Lebanon\tLE\t35.1\t33.0\t36.6\t34.7",
  "Lesotho\tLT\t27.0\t-30.7\t29.5\t-28.6",
  "Liberia\tLI\t-11.5\t4.3\t-7.4\t8.6",
  "Libya\tLY\t9.3\t19.5\t25.2\t33.2",
  "Liechtenstein\tLS\t9.4\t47.0\t9.6\t47.3",
  "Lithuania\tLH\t20.9\t53.9\t26.9\t56.4",
  "Luxembourg\tLU\t5.7\t49.4\t6.5\t50.2",
  "Macau\tMC\t113.5\t22.1\t113.6\t22.2",
  "Macedonia\tMK\t20.4\t40.8\t23.0\t42.4",
  "Madagascar\tMA\t43.1\t-25.7\t50.5\t-11.9",
  "Malawi\tMI\t32.6\t-17.2\t35.9\t-9.4",
  "Malaysia\tMY\t98.9\t5.6\t98.9\t5.7\t99.6\t1.2\t104.5\t6.7\t109.5\t0.8\t119.3\t7.4",
  "Maldives\tMV\t72.6\t-0.7\t73.7\t7.1",
  "Mali\tML\t-12.3\t10.1\t4.2\t25.0",
  "Malta\tMT\t14.1\t35.8\t14.6\t36.1",
  "Marshall Islands\tRM\t160.7\t4.5\t172.0\t14.8",
  "Martinique\tMB\t-61.3\t14.3\t-60.8\t14.9",
  "Mauritania\tMR\t-17.1\t14.7\t-4.8\t27.3",
  "Mauritius\tMP\t57.3\t-20.6\t57.8\t-20.0\t59.5\t-16.9\t59.6\t-16.7",
  "Mayotte\tMF\t45.0\t-13.1\t45.3\t-12.6",
  "Mexico\tMX\t-118.5\t28.8\t-118.3\t29.2\t-117.3\t14.5\t-86.7\t32.7",
  "Micronesia\tFM\t138.0\t9.4\t138.2\t9.6\t139.6\t9.8\t139.8\t10.0\t140.5\t9.7\t140.5\t9.8\t147.0\t7.3\t147.0\t7.4\t149.3\t6.6\t149.3\t6.7\t151.5\t7.1\t152.0\t7.5\t153.5\t5.2\t153.8\t5.6\t157.1\t5.7\t160.7\t7.1\t162.9\t5.2\t163.0\t5.4",
  "Midway Islands\tMQ\t-178.4\t28.3\t-178.3\t28.4\t-177.4\t28.1\t-177.3\t28.2\t-174.0\t26.0\t-174.0\t26.1\t-171.8\t25.7\t-171.7\t25.8",
  "Moldova\tMD\t26.6\t45.4\t30.2\t48.5",
  "Monaco\tMN\t7.3\t43.7\t7.5\t43.8",
  "Mongolia\tMG\t87.7\t41.5\t119.9\t52.2",
  "Montenegro\tMJ\t18.4\t42.2\t20.4\t43.6",
  "Montserrat\tMH\t-62.3\t16.6\t-62.1\t16.8",
  "Morocco\tMO\t-13.2\t27.6\t-1.0\t35.9",
  "Mozambique\tMZ\t30.2\t-26.9\t40.8\t-10.5",
  "Myanmar\tBM\t92.1\t9.6\t101.2\t28.5",
  "Namibia\tWA\t11.7\t-29.0\t25.3\t-17.0",
  "Nauru\tNR\t166.8\t-0.6\t166.9\t-0.5",
  "Navassa Island\tBQ\t-75.1\t18.3\t-75.0\t18.4",
  "Nepal\tNP\t80.0\t26.3\t88.2\t30.4",
  "Netherlands\tNL\t3.3\t50.7\t7.2\t53.6",
  "Netherlands Antilles\tNT\t-69.2\t11.9\t-68.2\t12.4\t-63.3\t17.4\t-62.9\t18.1",
  "New Caledonia\tNC\t163.5\t-22.8\t169.0\t-19.5",
  "New Zealand\tNZ\t166.4\t-48.1\t178.6\t-34.1",
  "Nicaragua\tNU\t-87.7\t10.7\t-82.6\t15.0",
  "Niger\tNG\t0.1\t11.6\t16.0\t23.5",
  "Nigeria\tNI\t2.6\t4.2\t14.7\t13.9",
  "Niue\tNE\t-170.0\t-19.2\t-169.8\t-19.0",
  "Norfolk Island\tNF\t168.0\t-29.2\t168.1\t-29.0",
  "North Korea\tKN\t124.1\t37.5\t130.7\t43.0",
  "Northern Mariana Islands\tCQ\t144.8\t14.1\t146.1\t20.6",
  "Norway\tNO\t4.6\t57.9\t31.1\t71.2",
  "Oman\tMU\t51.8\t16.6\t59.8\t25.0",
  "Pacific Ocean\tXX\t",
  "Pakistan\tPK\t60.8\t23.6\t77.8\t37.1",
  "Palau\tPS\t132.3\t4.3\t132.3\t4.3\t134.1\t6.8\t134.7\t7.7",
  "Palmyra Atoll\tLQ\t-162.2\t5.8\t-162.0\t5.9",
  "Panama\tPM\t-83.1\t7.1\t-77.2\t9.6",
  "Papua New Guinea\tPP\t140.8\t-11.7\t156.0\t-0.9\t157.0\t-4.9\t157.1\t-4.8\t159.4\t-4.7\t159.5\t-4.5",
  "Paracel Islands\tPF\t111.1\t15.7\t111.2\t15.8",
  "Paraguay\tPA\t-62.7\t-27.7\t-54.3\t-19.3",
  "Peru\tPE\t-81.4\t-18.4\t-68.7\t0.0",
  "Philippines\tRP\t116.9\t4.9\t126.6\t21.1",
  "Pitcairn Islands\tPC\t-128.4\t-24.5\t-128.3\t-24.3",
  "Poland\tPL\t14.1\t49.0\t24.2\t54.8",
  "Portugal\tPO\t-9.5\t36.9\t-6.2\t42.1\t-31.3\t36.9\t-25.0\t39.8\t-17.3\t32.4\t-16.2\t33.2",
  "Portugal: Azores\tXX\t-31.3\t36.9\t-25.0\t39.8",
  "Portugal: Madeira\tXX\t-17.3\t32.4\t-16.2\t33.2",
  "Puerto Rico\tRQ\t-68.0\t17.8\t-65.2\t18.5",
  "Qatar\tQA\t50.7\t24.4\t52.4\t26.2",
  "Republic of the Congo\tCF\t11.2\t-5.1\t18.6\t3.7",
  "Reunion\tRE\t55.2\t-21.4\t55.8\t-20.9",
  "Romania\tRO\t20.2\t43.6\t29.7\t48.3",
  "Russia\tRS\t-180.0\t64.2\t-169.0\t71.6\t19.7\t54.3\t22.9\t55.3\t26.9\t41.1\t180.0\t81.3",
  "Rwanda\tRW\t28.8\t-2.9\t30.9\t-1.1",
  "Saint Helena\tSH\t-5.8\t-16.1\t-5.6\t-15.9",
  "Saint Kitts and Nevis\tSC\t62.9\t17.0\t62.5\t17.5",
  "Saint Lucia\tST\t-61.1\t13.7\t-60.9\t14.1",
  "Saint Pierre and Miquelon\tSB\t-56.5\t46.7\t-56.2\t47.1",
  "Saint Vincent and the Grenadines\tVC\t-61.6\t12.4\t-61.1\t13.4",
  "Samoa\tWS\t-172.8\t-14.1\t-171.4\t-13.4",
  "San Marino\tSM\t12.4\t43.8\t12.5\t44.0",
  "Sao Tome and Principe\tTP\t6.4\t0.0\t1.7\t7.5",
  "Saudi Arabia\tSA\t34.4\t15.6\t55.7\t32.2",
  "Senegal\tSG\t-17.6\t12.3\t-11.4\t16.7",
  "Serbia\tRB\t18.8\t42.2\t23.1\t46.2",
  "Seychelles\tSE\t50.7\t-9.6\t51.1\t-9.2\t52.7\t-7.2\t52.8\t-7.0\t53.0\t-6.3\t53.7\t-5.1\t55.2\t-5.9\t56.0\t-3.7\t56.2\t-7.2\t56.3\t-7.1",
  "Sierra Leone\tSL\t-13.4\t6.9\t-10.3\t10.0",
  "Singapore\tSN\t103.6\t1.1\t104.1\t1.5",
  "Slovakia\tLO\t16.8\t47.7\t22.6\t49.6",
  "Slovenia\tSI\t13.3\t45.4\t16.6\t46.9",
  "Solomon Islands\tBP\t155.5\t-11.9\t162.8\t-5.1\t165.6\t-11.8\t167.0\t-10.1\t167.1\t-10.0\t167.3\t-9.8\t168.8\t-12.3\t168.8\t-12.3",
  "Somalia\tSO\t40.9\t-1.7\t51.4\t12.0",
  "South Africa\tSF\t16.4\t-34.9\t32.9\t-22.1",
  "South Georgia and the South Sandwich Islands\tSX\t-38.3\t-54.9\t-35.7\t-53.9",
  "South Korea\tKS\t125.0\t33.1\t129.6\t38.6",
  "Spain\tSP\t-9.3\t35.1\t4.3\t43.8\t-18.2\t27.6\t-13.4\t29.5",
  "Spain: Canary Islands\tXX\t-18.2\t27.6\t-13.4\t29.5",
  "Spratly Islands\tPG\t114.0\t9.6\t115.8\t11.1",
  "Sri Lanka\tCE\t79.6\t5.9\t81.9\t9.8",
  "Sudan\tSU\t21.8\t3.4\t38.6\t23.6",
  "Suriname\tNS\t-58.1\t1.8\t-54.0\t6.0",
  "Svalbard\tSV\t10.4\t76.4\t33.5\t80.8",
  "Swaziland\tWZ\t30.7\t-27.4\t32.1\t-25.7",
  "Sweden\tSW\t10.9\t55.3\t24.2\t69.1",
  "Switzerland\tSZ\t5.9\t45.8\t10.5\t47.8",
  "Syria\tSY\t35.7\t32.3\t42.4\t37.3",
  "Taiwan\tTW\t119.3\t21.9\t122.0\t25.3",
  "Tajikistan\tTI\t67.3\t36.6\t75.1\t41.0",
  "Tanzania\tTZ\t29.3\t-11.8\t40.4\t-1.0",
  "Thailand\tTH\t97.3\t5.6\t105.6\t20.5",
  "Togo\tTO\t-0.2\t6.1\t1.8\t11.1",
  "Tokelau\tTL\t-172.6\t-9.5\t-171.1\t-8.5",
  "Tonga\tTN\t-176.3\t-22.4\t-176.2\t-22.3\t-175.5\t-21.5\t-174.5\t-20.0",
  "Trinidad and Tobago\tTD\t-62.0\t10.0\t-60.5\t11.3",
  "Tromelin Island\tTE\t54.5\t-15.9\t54.5\t-15.9",
  "Tunisia\tTS\t7.5\t30.2\t11.6\t37.5",
  "Turkey\tTU\t25.6\t35.8\t44.8\t42.1",
  "Turkmenistan\tTX\t52.4\t35.1\t66.7\t42.8",
  "Turks and Caicos Islands\tTK\t-73.8\t20.9\t-73.0\t21.3",
  "Tuvalu\tTV\t176.0\t-7.3\t177.3\t-5.6\t178.4\t-8.0\t178.7\t-7.4\t179.0\t-9.5\t179.9\t-8.5",
  "Uganda\tUG\t29.5\t-1.5\t35.0\t4.2",
  "Ukraine\tUP\t22.1\t44.3\t40.2\t52.4",
  "United Arab Emirates\tAE\t51.1\t22.4\t56.4\t26.1",
  "United Kingdom\tUK\t-8.7\t49.7\t1.8\t60.8",
  "Uruguay\tUY\t-58.5\t-35.0\t-53.1\t-30.1",
  "USA\tUS\t-124.8\t24.5\t-66.9\t49.4\t-168.2\t54.3\t-130.0\t71.4\t172.4\t52.3\t176.0\t53.0\t177.2\t51.3\t179.8\t52.1\t-179.5\t51.0\t-172.0\t52.5\t-171.5\t52.0\t-164.5\t54.5\t-164.8\t23.5\t-164.7\t23.6\t-162.0\t23.0\t-161.9\t23.1\t-160.6\t18.9\t-154.8\t22.2",
  "USA: Alabama\tXX\t-88.8\t30.1\t-84.9\t35.0",
  "USA: Alaska\tXX\t-168.2\t54.3\t-130.0\t71.4\t172.4\t52.3\t176.0\t53.0\t177.2\t51.3\t179.8\t52.1\t-179.5\t51.0\t-172.0\t52.5\t-171.5\t52.0\t-164.5\t54.5",
  "USA: Alaska, Aleutian Islands\tXX\t172.4\t52.3\t176.0\t53.0\t177.2\t51.3\t179.8\t52.1\t-179.5\t51.0\t-172.0\t52.5\t-171.5\t52.0\t-164.5\t54.5",
  "USA: Arizona\tXX\t-114.9\t31.3\t-109.0\t37.0",
  "USA: Arkansas\tXX\t-94.7\t33.0\t-89.6\t36.5",
  "USA: California\tXX\t-124.5\t32.5\t-114.1\t42.0",
  "USA: Colorado\tXX\t-109.1\t36.9\t-102.0\t41.0",
  "USA: Connecticut\tXX\t-73.8\t40.9\t-71.8\t42.1",
  "USA: Delaware\tXX\t-75.8\t38.4\t-74.9\t39.8",
  "USA: Florida\tXX\t-87.7\t24.5\t-80.0\t31.0",
  "USA: Georgia\tXX\t-85.7\t30.3\t-80.8\t35.0",
  "USA: Hawaii\tXX\t-164.8\t23.5\t-164.7\t23.6\t-162.0\t23.0\t-161.9\t23.1\t-160.6\t18.9\t-154.8\t22.2",
  "USA: Idaho\tXX\t-117.3\t41.9\t-111.0\t49.0",
  "USA: Illinois\tXX\t-91.6\t36.9\t-87.0\t42.5",
  "USA: Indiana\tXX\t-88.1\t37.7\t-84.8\t41.8",
  "USA: Iowa\tXX\t-96.7\t40.3\t-90.1\t43.5",
  "USA: Kansas\tXX\t-102.1\t36.9\t-94.6\t40.0",
  "USA: Kentucky\tXX\t-89.5\t36.5\t-82.0\t39.1",
  "USA: Louisiana\tXX\t-94.1\t28.9\t-88.8\t33.0",
  "USA: Maine\tXX\t-71.1\t43.0\t-66.9\t47.5",
  "USA: Maryland\tXX\t-79.5\t37.8\t-75.1\t39.7",
  "USA: Massachusetts\tXX\t-73.6\t41.2\t-69.9\t42.9",
  "USA: Michigan\tXX\t-90.5\t41.6\t-82.1\t48.3",
  "USA: Minnesota\tXX\t-97.3\t43.4\t-90.0\t49.4",
  "USA: Mississippi\tXX\t-91.7\t30.1\t-88.1\t35.0",
  "USA: Missouri\tXX\t-95.8\t36.0\t-89.1\t40.6",
  "USA: Montana\tXX\t-116.1\t44.3\t-104.0\t49.0",
  "USA: Nebraska\tXX\t-104.1\t40.0\t-95.3\t43.0",
  "USA: Nevada\tXX\t-120.0\t35.0\t-114.0\t42.0",
  "USA: New Hampshire\tXX\t-72.6\t42.6\t-70.7\t45.3",
  "USA: New Jersey\tXX\t-75.6\t38.9\t-73.9\t41.4",
  "USA: New Mexico\tXX\t-109.1\t31.3\t-103.0\t37.0",
  "USA: New York\tXX\t-79.8\t40.4\t-71.9\t45.0",
  "USA: North Carolina\tXX\t-84.4\t33.8\t-75.5\t36.6",
  "USA: North Dakota\tXX\t-104.1\t45.9\t-96.6\t49.0",
  "USA: Ohio\tXX\t-84.9\t38.3\t-80.5\t42.3",
  "USA: Oklahoma\tXX\t-103.1\t33.6\t-94.4\t37.0",
  "USA: Oregon\tXX\t-124.6\t41.9\t-116.5\t46.3",
  "USA: Pennsylvania\tXX\t-80.6\t39.7\t-74.7\t42.5",
  "USA: Rhode Island\tXX\t-71.9\t41.1\t-71.1\t42.0",
  "USA: South Carolina\tXX\t-83.4\t32.0\t-78.6\t35.2",
  "USA: South Dakota\tXX\t-104.1\t42.4\t-96.4\t45.9",
  "USA: Tennessee\tXX\t-90.4\t35.0\t-81.7\t36.7",
  "USA: Texas\tXX\t-106.7\t25.8\t-93.5\t36.5",
  "USA: Utah\tXX\t-114.1\t37.0\t-109.1\t42.0",
  "USA: Vermont\tXX\t-73.5\t42.7\t-71.5\t45.0",
  "USA: Virginia\tXX\t-83.7\t36.5\t-75.2\t39.5",
  "USA: Washington\tXX\t-124.8\t45.5\t-116.9\t49.0",
  "USA: West Virginia\tXX\t-82.7\t37.1\t-77.7\t40.6",
  "USA: Wisconsin\tXX\t-92.9\t42.4\t-86.3\t47.3",
  "USA: Wyoming\tXX\t-111.1\t40.9\t-104.1\t45.0",
  "Uzbekistan\tUZ\t55.9\t37.1\t73.1\t45.6",
  "Vanuatu\tNH\t166.5\t-20.3\t170.2\t-13.1",
  "Venezuela\tVE\t-73.4\t0.7\t-59.8\t12.2",
  "Viet Nam\tVM\t102.1\t8.4\t109.5\t23.4",
  "Virgin Islands\tVQ\t-65.1\t17.6\t-64.6\t18.5",
  "Wake Island\tWQ\t166.5\t19.2\t166.7\t19.3",
  "Wallis and Futuna\tWF\t-178.3\t-14.4\t-178.0\t-14.2\t-176.3\t-13.4\t-176.1\t-13.2",
  "West Bank\tWE\t34.8\t31.3\t35.6\t32.6",
  "Western Sahara\tWI\t-17.2\t20.7\t-8.7\t27.7",
  "Yemen\tYM\t41.8\t11.7\t54.5\t19.0",
  "Zambia\tZA\t21.9\t-18.1\t33.7\t-8.2",
  "Zimbabwe\tZI\t25.2\t-22.5\t33.1\t-15.6",
  NULL
};


/* one CtBlock for each discontiguous area per country */

typedef struct ctblock {
  CharPtr  country;       /* points to instance in countries list */
  FloatHi  minx;
  FloatHi  miny;
  FloatHi  maxx;
  FloatHi  maxy;
} CtBlock, PNTR CtBlockPtr;

/* one CtGrid for each 10-degree-by-10-degree area touched by a CtBlock */

typedef struct ctgrid {
  CtBlockPtr  cbp;
  Int2        xindex;
  Int2        yindex;
} CtGrid, PNTR CtGridPtr;

/* main structure for country/lat-lon lookup */

typedef struct ctset {
  ValNodePtr       countries;
  ValNodePtr       blocks;
  ValNodePtr       grids;
  CtBlockPtr PNTR  bkarray;     /* sorted by country name */
  CtGridPtr PNTR   gdarray;     /* sorted by geographic index */
  Int4             num_blocks;
  Int4             num_grids;
} CtSet, PNTR CtSetPtr;

static int LIBCALLBACK SortCbpByCountry (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  int         compare;
  CtBlockPtr  cbp1, cbp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  cbp1 = *((CtBlockPtr PNTR) ptr1);
  cbp2 = *((CtBlockPtr PNTR) ptr2);
  if (cbp1 == NULL || cbp2 == NULL) return 0;

  compare = StringICmp (cbp1->country, cbp2->country);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  return 0;
}

static int CgpGridComp (
  CtGridPtr cgp1,
  Int2 xindex,
  Int2 yindex
)

{
  if (cgp1 == NULL) return 0;

  if (cgp1->xindex > xindex) {
    return 1;
  } else if (cgp1->xindex < xindex) {
    return -1;
  }

  if (cgp1->yindex > yindex) {
    return 1;
  } else if (cgp1->yindex < yindex) {
    return -1;
  }

  return 0;
}

static int LIBCALLBACK SortCgpByGrid (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  CtBlockPtr  cbp1, cbp2;
  CtGridPtr   cgp1, cgp2;
  int         compare;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  cgp1 = *((CtGridPtr PNTR) ptr1);
  cgp2 = *((CtGridPtr PNTR) ptr2);
  if (cgp1 == NULL || cgp2 == NULL) return 0;

  compare = CgpGridComp (cgp1, cgp2->xindex, cgp2->yindex);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  cbp1 = cgp1->cbp;
  cbp2 = cgp2->cbp;
  if (cbp1 == NULL || cbp2 == NULL) return 0;

  if (cbp1->minx > cbp2->minx) {
    return 1;
  } else if (cbp1->minx < cbp2->minx) {
    return -1;
  }

  if (cbp1->maxx > cbp2->maxx) {
    return -1;
  } else if (cbp1->maxx < cbp2->maxx) {
    return 1;
  }

  if (cbp1->miny > cbp2->miny) {
    return 1;
  } else if (cbp1->miny < cbp2->miny) {
    return -1;
  }

  if (cbp1->maxy > cbp2->maxy) {
    return -1;
  } else if (cbp1->maxy < cbp2->maxy) {
    return 1;
  }

  compare = StringICmp (cbp1->country, cbp2->country);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  return 0;
}

static Int2 LatLonDegreeToIndex (
  FloatHi coord
)

{
  double  fval;
  long    ival;

  fval = coord;
  fval += 200.0;
  fval /= 10.0;
  ival = (long) fval;
  ival -= 20;

  return (Int2) ival;
}

static CtSetPtr CtSetDataFree (
  CtSetPtr csp
)

{
  if (csp == NULL) return NULL;

  ValNodeFreeData (csp->countries);
  ValNodeFreeData (csp->blocks);
  ValNodeFreeData (csp->grids);

  MemFree (csp->bkarray);
  MemFree (csp->gdarray);

  MemFree (csp);

  return NULL;
}

static Boolean ct_set_not_found = FALSE;

static CtSetPtr GetCtSetLatLonDataInt (
  CharPtr prop,
  CharPtr file,
  CharPtr PNTR local
)

{
  CtBlockPtr PNTR  bkarray;
  ValNodePtr       blocks = NULL;
  FloatHi          bounds [4];
  CtBlockPtr       cbp;
  CtGridPtr        cgp;
  ValNodePtr       countries = NULL;
  CharPtr          country;
  CtSetPtr         csp;
  FileCache        fc;
  FILE             *fp = NULL;
  CtGridPtr PNTR   gdarray;
  ValNodePtr       grids = NULL;
  Int2             hix;
  Int2             hiy;
  Int2             i;
  Int2             j = 0;
  ValNodePtr       lastblk = NULL;
  ValNodePtr       lastctry = NULL;
  ValNodePtr       lastgrd = NULL;
  Char             line [1024];
  Int2             lox;
  Int2             loy;
  Int4             num;
  Char             path [PATH_MAX];
  CharPtr          ptr;
  ErrSev           sev;
  CharPtr          str = NULL;
  double           val;
  ValNodePtr       vnp;
  CharPtr          wrk;
  Int2             x;
  Int2             y;

  csp = (CtSetPtr) GetAppProperty (prop);
  if (csp != NULL) return csp;

  if (ct_set_not_found) return NULL;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }

  if (fp != NULL) {
    FileCacheSetup (&fc, fp);
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  } else if (local != NULL) {
    str = local [j];
    if (str != NULL) {
      StringNCpy_0 (line, str, sizeof (line));
      str = line;
    }
  } else return NULL;

  while (str != NULL) {
    if (StringDoesHaveText (str)) {
      ptr = StringChr (str, '\t');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
        ptr = StringChr (ptr, '\t');
        if (ptr != NULL) {
          ptr++;
          if (StringDoesHaveText (str) && StringDoesHaveText (ptr)) {

            country = StringSave (str);

            vnp = ValNodeAddPointer (&lastctry, 0, (Pointer) country);
            if (countries == NULL) {
              countries = vnp;
            }
            lastctry = vnp;

            wrk = StringSave (ptr);
            str = wrk;
            i = 0;

            while (StringDoesHaveText (str)) {
              ptr = StringChr (str, '\t');
              if (ptr != NULL) {
                *ptr = '\0';
                ptr++;
              }

              if (sscanf (str, "%lf", &val) == 1) {
                bounds [i] = (FloatHi) val;
                i++;
                if (i > 3) {

                  cbp = (CtBlockPtr) MemNew (sizeof (CtBlock));
                  if (cbp != NULL) {
                    cbp->country = country;
                    cbp->minx = bounds [0];
                    cbp->miny = bounds [1];
                    cbp->maxx = bounds [2];
                    cbp->maxy = bounds [3];

                    vnp = ValNodeAddPointer (&lastblk, 0, (Pointer) cbp);
                    if (blocks == NULL) {
                      blocks = vnp;
                    }
                    lastblk = vnp;

                    lox = LatLonDegreeToIndex (cbp->minx);
                    loy = LatLonDegreeToIndex (cbp->miny);
                    hix = LatLonDegreeToIndex (cbp->maxx);
                    hiy = LatLonDegreeToIndex (cbp->maxy);

                    for (x = lox; x <= hix; x++) {
                      for (y = loy; y <= hiy; y++) {
                        cgp = (CtGridPtr) MemNew (sizeof (CtGrid));
                        if (cgp != NULL) {
                          cgp->cbp = cbp;
                          cgp->xindex = x;
                          cgp->yindex = y;

                          vnp = ValNodeAddPointer (&lastgrd, 0, (Pointer) cgp);
                          if (grids == NULL) {
                            grids = vnp;
                          }
                          lastgrd = vnp;
                        }
                      }
                    }
                  }

                  i = 0;
                }
              }

              str = ptr;
            }

            MemFree (wrk);
          }
        }
      }
    }

    if (fp != NULL) {
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
    } else {
      j++;
      str = local [j];
      if (str != NULL) {
        StringNCpy_0 (line, str, sizeof (line));
        str = line;
      }
    }
  }

  if (fp != NULL) {
    FileClose (fp);
  }

  if (countries == NULL || blocks == NULL || grids == NULL) {
    ct_set_not_found = TRUE;
    return NULL;
  }

  csp = (CtSetPtr) MemNew (sizeof (CtSet));
  if (csp == NULL) return NULL;

  /* now populate, heap sort arrays */

  num = ValNodeLen (blocks);

  csp->countries = countries;
  csp->blocks = blocks;
  csp->num_blocks = (Int2) num;

  bkarray = (CtBlockPtr PNTR) MemNew (sizeof (CtBlockPtr) * (num + 1));
  if (bkarray != NULL) {
    for (vnp = blocks, i = 0; vnp != NULL; vnp = vnp->next, i++) {
      cbp = (CtBlockPtr) vnp->data.ptrvalue;
      bkarray [i] = cbp;
    }

    HeapSort (bkarray, (size_t) num, sizeof (CtBlockPtr), SortCbpByCountry);
    csp->bkarray = bkarray;
  }

  num = ValNodeLen (grids);

  csp->num_grids = (Int2) num;

  gdarray = (CtGridPtr PNTR) MemNew (sizeof (CtGridPtr) * (num + 1));
  if (gdarray != NULL) {
    for (vnp = grids, i = 0; vnp != NULL; vnp = vnp->next, i++) {
      cgp = (CtGridPtr) vnp->data.ptrvalue;
      gdarray [i] = cgp;
    }

    HeapSort (gdarray, (size_t) num, sizeof (CtGridPtr), SortCgpByGrid);
    csp->gdarray = gdarray;
  }

  SetAppProperty (prop, (Pointer) csp);

  return csp;
}

static CtSetPtr GetCtSetLatLonData (
  void
)

{
  return GetCtSetLatLonDataInt ("CountryLatLonList", "country_lat_lon.txt", ctry_lat_lon);
}

NLM_EXTERN Boolean IsCountryInLatLonList (
  CharPtr country
)

{
  CtBlockPtr       cbp;
  CtBlockPtr PNTR  bkarray;
  CtSetPtr         csp;
  Int2             L, R, mid;

  if (StringHasNoText (country)) return FALSE;

  csp = GetCtSetLatLonData ();
  if (csp == NULL) return FALSE;

  bkarray = csp->bkarray;
  if (bkarray == NULL) return FALSE;

  L = 0;
  R = csp->num_blocks - 1;

  while (L < R) {
    mid = (L + R) / 2;
    cbp = bkarray [mid];
    if (cbp != NULL && StringICmp (cbp->country, country) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  cbp = bkarray [R];
  if (cbp != NULL && StringICmp (cbp->country, country) == 0) return TRUE;

  return FALSE;
}

NLM_EXTERN Boolean TestLatLonForCountry (
  CharPtr country,
  FloatHi lat,
  FloatHi lon
)

{
  CtBlockPtr       cbp;
  CtBlockPtr PNTR  bkarray;
  CtSetPtr         csp;
  Int2             L, R, mid;

  if (StringHasNoText (country)) return FALSE;

  csp = GetCtSetLatLonData ();
  if (csp == NULL) return FALSE;

  bkarray = csp->bkarray;
  if (bkarray == NULL) return FALSE;

  L = 0;
  R = csp->num_blocks - 1;

  while (L < R) {
    mid = (L + R) / 2;
    cbp = bkarray [mid];
    if (cbp != NULL && StringICmp (cbp->country, country) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  while (R < csp->num_blocks) {
    cbp = bkarray [R];
    if (cbp == NULL) return FALSE;
    if (StringICmp (cbp->country, country) != 0) return FALSE;
    if (lon >= cbp->minx && lat >= cbp->miny && lon <= cbp->maxx && lat <= cbp->maxy) return TRUE;
    R++;
  }

  return FALSE;
}

NLM_EXTERN CharPtr GuessCountryForLatLon (
  FloatHi lat,
  FloatHi lon
)

{
  CtBlockPtr      cbp;
  CtGridPtr       cgp;
  CharPtr         country = NULL;
  CtSetPtr        csp;
  CtGridPtr PNTR  gdarray;
  Int2            L, R, mid;
  Int2            x;
  Int2            y;

  csp = GetCtSetLatLonData ();
  if (csp == NULL) return NULL;

  gdarray = csp->gdarray;
  if (gdarray == NULL) return NULL;

  L = 0;
  R = csp->num_grids - 1;

  x = LatLonDegreeToIndex (lon);
  y = LatLonDegreeToIndex (lat);

  while (L < R) {
    mid = (L + R) / 2;
    cgp = gdarray [mid];
    if (cgp != NULL && CgpGridComp (cgp, x, y) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  while (R < csp->num_grids) {
    cgp = gdarray [R];
    if (cgp == NULL) return country;
    if (cgp->xindex != x || cgp->yindex != y) return country;
    cbp = cgp->cbp;
    if (cbp == NULL) return country;
    if (lon >= cbp->minx && lat >= cbp->miny && lon <= cbp->maxx && lat <= cbp->maxy) {
      country = cbp->country;
    }
    R++;
  }

  return country;
}


static CharPtr bodiesOfWater [] = {
  "Bay",
  "Canal",
  "Channel",
  "Coastal",
  "Cove",
  "Estuary",
  "Fjord",
  "Freshwater",
  "Gulf",
  "Harbor",
  "Inlet",
  "Lagoon",
  "Lake",
  "Narrows",
  "Ocean",
  "Passage",
  "River",
  "Sea",
  "Seawater",
  "Sound",
  "Strait",
  "Water",
  "Waters",
  NULL
};

static TextFsaPtr GetBodiesOfWaterFSA (void)


{
  TextFsaPtr  fsa;
  Int2        i;
  CharPtr     prop = "BodiesOfWaterFSA";

  fsa = (TextFsaPtr) GetAppProperty (prop);
  if (fsa != NULL) return fsa;

  fsa = TextFsaNew ();
  if (fsa != NULL) {
    for (i = 0; bodiesOfWater [i] != NULL; i++) {
      TextFsaAdd (fsa, bodiesOfWater [i]);
    }
  }

  SetAppProperty (prop, (Pointer) fsa);

  return fsa;
}

NLM_EXTERN Boolean StringContainsBodyOfWater (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  CharPtr     ptr;
  Int4        state;
  ValNodePtr  matches;

  if (StringHasNoText (str)) return FALSE;

  fsa = GetBodiesOfWaterFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  ptr = str;
  ch = *ptr;

  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (fsa, state, ch, &matches);
    ptr++;
    ch = *ptr;
    if (ch == '\0' || ch == ',' || ch == ':' || ch == ';' || ch == ' ') {
      if (matches != NULL) return TRUE;
      state = 0;
    }
  }

  return FALSE;
}




static CharPtr GetDash (CharPtr str)

{
  Char  ch;

  if (str == NULL) return NULL;
  ch = *str;
  while (ch != '\0') {
    if (ch == '-') return str;
    str++;
    ch = *str;
  }

  return NULL;
}

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


static Boolean CollectionDateIsValid (CharPtr name)

{
  Char      ch;
  Int2      i;
  CharPtr   ptr1, ptr2, month = NULL, day = NULL, year = NULL;
  Char      str [256];
  long int  val;

  if (StringHasNoText (name)) return FALSE;

  StringNCpy_0 (str, name, sizeof (str));
  ptr1 = GetDash (str);
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = GetDash (ptr1);
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
    if (sscanf (day, "%ld", &val) != 1 || val < 1 || val > 31) return FALSE;
    if (StringLen (day) != 2 || !isdigit(day[0]) || !isdigit(day[1])) return FALSE;
  }

  if (month != NULL) {
    for (i = 0; legalMonths [i] != NULL; i++) {
      if (StringCmp (month, legalMonths [i]) == 0) {
        break;
      }
    }
    if (legalMonths [i] == NULL) return FALSE;
  }

  if (year != NULL) {
    ptr1 = year;
    ch = *ptr1;
    while (ch != '\0') {
      if (! (IS_DIGIT (ch))) return FALSE;
      ptr1++;
      ch = *ptr1;
    }
    if (sscanf (year, "%ld", &val) == 1) {
      if (val >= 1700 && val < 2100) return TRUE;
    }
  }

  return FALSE;
}

/* This mimics a portion of the DatePtr structure,
 * but allows dates with years before 1900 because
 * the year value is Int4 instead of Uint1.
 *   data [0] : Set to 1
 *        [1] - year (- 1900)
 *        [2] - month (1-12)  optional
 *        [3] - day (1-31)     optional
 * Not bothering with time.
 */

typedef struct betterdate {
    Int4 data[8];      /* see box above */
} BetterDateData, PNTR BetterDatePtr;

static BetterDatePtr BetterDateNew()
{
  BetterDatePtr dp;

  dp = (BetterDatePtr) MemNew (sizeof (BetterDateData));
  return dp;
}

static BetterDatePtr BetterDateFree (BetterDatePtr dp)
{
  if (dp != NULL) {
    dp = MemFree (dp);
  }
  return dp;
}


static BetterDatePtr CollectionDateFromString (CharPtr name)
{
  Char      ch;
  Int2      i;
  CharPtr   ptr1, ptr2, month = NULL, day = NULL, year = NULL;
  Char      str [256];
  long int  day_val = 0;
  Int2      month_num = 0;
  long int  val, year_val = 0;
  BetterDatePtr   dp;

  if (StringHasNoText (name)) return NULL;

  StringNCpy_0 (str, name, sizeof (str));
  ptr1 = GetDash (str);
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = GetDash (ptr1);
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
    if (sscanf (day, "%ld", &day_val) != 1 || day_val < 1 || day_val > 31) return NULL;
  }

  if (month != NULL) {
    for (i = 0; legalMonths [i] != NULL; i++) {
      if (StringCmp (month, legalMonths [i]) == 0) {
        month_num = i + 1;
        break;
      }
    }
    if (legalMonths [i] == NULL) return NULL;
  }

  if (year != NULL) {
    ptr1 = year;
    ch = *ptr1;
    while (ch != '\0') {
      if (! (IS_DIGIT (ch))) return NULL;
      ptr1++;
      ch = *ptr1;
    }
    if (sscanf (year, "%ld", &val) == 1) {
      if (val < 1700 || val > 2100) return NULL;
      year_val = val - 1900;
    }
    else
    {
      return NULL;
    }
  }

  dp = BetterDateNew();
  dp->data[0] = 1;
  dp->data[1] = year_val;
  dp->data[2] = month_num;
  dp->data[3] = day_val;
  return dp;
}


static Boolean CollectionDateIsInTheFuture (CharPtr name)

{
  DatePtr   dp_now;
  BetterDatePtr dp_coll_date;
  Boolean   rval = FALSE;

  dp_coll_date = CollectionDateFromString (name);
  if (dp_coll_date == NULL) return FALSE;

  if (dp_coll_date->data[1] < 0)
  {
    /* year before 1900 */
    dp_coll_date = BetterDateFree (dp_coll_date);
    return FALSE;
  }

  dp_now = DateCurr();

  /* compare years */
  if (dp_now->data[1] < dp_coll_date->data[1])
  {
    rval = TRUE;
  }
  else if (dp_now->data[1] > dp_coll_date->data[1])
  {
    rval = FALSE;
  }
  /* years are equal - compare months */
  else if (dp_now->data[2] < dp_coll_date->data[2])
  {
    rval = TRUE;
  }
  else if (dp_now->data[2] > dp_coll_date->data[2])
  {
    rval = FALSE;
  }
  /* years and months are equal - compare days */
  else if (dp_now->data[3] < dp_coll_date->data[3])
  {
    rval = TRUE;
  }
  else
  {
    rval = FALSE;
  }

  dp_now = DateFree (dp_now);
  dp_coll_date = BetterDateFree (dp_coll_date);
  return rval;
}


static Boolean StringListIsUnique (ValNodePtr list)

{
  CharPtr     last;
  ValNodePtr  next;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL) return TRUE;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      return FALSE;
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
    }
    vnp = next;
  }

  return TRUE;
}

static CharPtr modified_base_abbrevs [] = {
  "<ac4c>",
  "<chm5u>",
  "<cm>",
  "<cmnm5s2u>",
  "<cmnm5u>",
  "<d>",
  "<fm>",
  "<gal q>",
  "<gm>",
  "<i>",
  "<i6a>",
  "<m1a>",
  "<m1f>",
  "<m1g>",
  "<m1i>",
  "<m22g>",
  "<m2a>",
  "<m2g>",
  "<m3c>",
  "<m5c>",
  "<m6a>",
  "<m7g>",
  "<mam5u>",
  "<mam5s2u>",
  "<man q>",
  "<mcm5s2u>",
  "<mcm5u>",
  "<mo5u>",
  "<ms2i6a>",
  "<ms2t6a>",
  "<mt6a>",
  "<mv>",
  "<o5u>",
  "<osyw>",
  "<p>",
  "<q>",
  "<s2c>",
  "<s2t>",
  "<s2u>",
  "<s4u>",
  "<t>",
  "<t6a>",
  "<tm>",
  "<um>",
  "<yw>",
  "<x>",
  "<OTHER>",
  NULL
};

static void InitializeModBaseFSA (ValidStructPtr vsp)

{
  Int2  i;

  vsp->modifiedBases = TextFsaNew ();
  for (i = 0; modified_base_abbrevs [i] != NULL; i++) {
    TextFsaAdd (vsp->modifiedBases, modified_base_abbrevs [i]);
  }
}

static Boolean PrimerSeqIsValid (ValidStructPtr vsp, CharPtr name, Char PNTR badch)

{
  Char        ch;
  TextFsaPtr  fsa;
  size_t      len;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  if (badch != NULL) {
    *badch = '\0';
  }

  if (vsp == NULL) return FALSE;
  if (vsp->modifiedBases == NULL) {
    InitializeModBaseFSA (vsp);
  }
  fsa = vsp->modifiedBases;
  if (fsa == NULL) return FALSE;

  if (StringHasNoText (name)) return FALSE;
  len = StringLen (name);
  if (len < 1) return FALSE;

  if (StringChr (name, ',') != NULL) {
    if (name [0] != '(' || name [len - 1] != ')') return FALSE;
  } else {
    if (StringChr (name, '(') != NULL) return FALSE;
    if (StringChr (name, ')') != NULL) return FALSE;
  }

  if (StringChr (name, ';') != NULL) return FALSE;
  if (StringChr (name, ' ') != NULL) return FALSE;

  ptr = name;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '<') {
      state = 0;
      matches = NULL;
      while (ch != '\0' && ch != '>') {
        state = TextFsaNext (fsa, state, ch, &matches);
        ptr++;
        ch = *ptr;
      }
      if (ch != '>') {
        if (badch != NULL) {
          *badch = ch;
        }
        return FALSE;
      }
      state = TextFsaNext (fsa, state, ch, &matches);
      if (matches == NULL) {
        if (badch != NULL) {
          *badch = ch;
        }
        return FALSE;
      }
    } else {
      if (ch != '(' && ch != ')' && ch != ',' && ch != ':') {
        if (! (IS_ALPHA (ch))) {
          if (badch != NULL) {
            *badch = ch;
          }
          return FALSE;
        }
        ch = TO_UPPER (ch);
        if (StringChr ("ABCDGHKMNRSTVWY", ch) == NULL) {
          if (badch != NULL) {
            ch = TO_LOWER (ch);
            *badch = ch;
          }
          return FALSE;
        }
      }
    }
    ptr++;
    ch = *ptr;
  }

  return TRUE;
}

/*
static ValNodePtr ParsePrimerSeqIntoComponents (
  CharPtr strs
)

{
  Char        ch;
  ValNodePtr  head = NULL;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  if (tmp == NULL) return NULL;

  str = tmp;
  while (StringDoesHaveText (str)) {
    ptr = str;
    ch = *ptr;

    while (ch != '\0' && ch != '(' && ch != ')' && ch != ',' && ch != ';' && ch != ':') {
      ptr++;
      ch = *ptr;
    }
    if (ch != '\0' && ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }

    TrimSpacesAroundString (str);
    if (StringDoesHaveText (str)) {
      ValNodeCopyStr (&head, 0, str);
    }

    str = ptr;
  }

  MemFree (tmp);
  return head;
}

static Boolean PrimerSeqHasDuplicates (CharPtr name)

{
  ValNodePtr  head;
  Boolean     rsult = FALSE;

  if (StringHasNoText (name)) return FALSE;

  head = ParsePrimerSeqIntoComponents (name);
  if (head == NULL) return FALSE;
  head = ValNodeSort (head, SortVnpByString);
  if (! StringListIsUnique (head)) {
    rsult = TRUE;
  }
  ValNodeFreeData (head);

  return rsult;
}
*/

static Int2 CountDigits (CharPtr str)

{
  Char  ch;
  Int2  count = 0;

  if (str == NULL) return count;
  ch = *str;
  while (IS_DIGIT (ch)) {
    count++;
    str++;
    ch = *str;
  }
  return count;
}

static Boolean LatLonIsValid (CharPtr name)

{
  Char     ch;
  Int2     count;
  CharPtr  str;

  if (StringHasNoText (name)) return FALSE;
  str = name;

  count = CountDigits (str);
  if (count < 1 || count > 2) return FALSE;
  str += count;

  ch = *str;
  if (ch == '.') {
    str++;
    count = CountDigits (str);
    if (count != 2) return FALSE;
    str += count;
  }

  ch = *str;
  if (ch != ' ') return FALSE;
  str++;
  ch = *str;
  if (ch != 'N' && ch != 'S') return FALSE;
  str++;
  ch = *str;
  if (ch != ' ') return FALSE;
  str++;

  count = CountDigits (str);
  if (count < 1 || count > 3) return FALSE;
  str += count;

  ch = *str;
  if (ch == '.') {
    str++;
    count = CountDigits (str);
    if (count != 2) return FALSE;
    str += count;
  }

  ch = *str;
  if (ch != ' ') return FALSE;
  str++;
  ch = *str;
  if (ch != 'E' && ch != 'W') return FALSE;
  str++;

  ch = *str;
  if (ch != '\0') return FALSE;

  return TRUE;
}

static CharPtr source_qual_prefixes [] = {
  "acronym:",
  "anamorph:",
  "authority:",
  "biotype:",
  "biovar:",
  "bio_material:",
  "breed:",
  "cell_line:",
  "cell_type:",
  "chemovar:",
  "chromosome:",
  "clone:",
  "clone_lib:",
  "collected_by:",
  "collection_date:",
  "common:",
  "country:",
  "cultivar:",
  "culture_collection:",
  "dev_stage:",
  "dosage:",
  "ecotype:",
  "endogenous_virus_name:",
  "environmental_sample:",
  "forma:",
  "forma_specialis:",
  "frequency:",
  "fwd_pcr_primer_name",
  "fwd_pcr_primer_seq",
  "fwd_primer_name",
  "fwd_primer_seq",
  "genotype:",
  "germline:",
  "group:",
  "haplogrop:",
  "haplotype:",
  "identified_by:",
  "insertion_seq_name:",
  "isolate:",
  "isolation_source:",
  "lab_host:",
  "lat_lon:"
  "left_primer:",
  "linkage_group:",
  "map:",
  "mating_type:",
  "metagenome_source:",
  "metagenomic:",
  "nat_host:",
  "pathovar:",
  "placement:",
  "plasmid_name:",
  "plastid_name:",
  "pop_variant:",
  "rearranged:",
  "rev_pcr_primer_name",
  "rev_pcr_primer_seq",
  "rev_primer_name",
  "rev_primer_seq",
  "right_primer:",
  "segment:",
  "serogroup:",
  "serotype:",
  "serovar:",
  "sex:",
  "specimen_voucher:",
  "strain:",
  "subclone:",
  "subgroup:",
  "substrain:",
  "subtype:",
  "sub_species:",
  "synonym:",
  "taxon:",
  "teleomorph:",
  "tissue_lib:",
  "tissue_type:",
  "transgenic:",
  "transposon_name:",
  "type:",
  "variety:",
  NULL
};

static void InitializeSourceQualTags (ValidStructPtr vsp)

{
  Int2  i;

  vsp->sourceQualTags = TextFsaNew ();
  for (i = 0; source_qual_prefixes [i] != NULL; i++) {
    TextFsaAdd (vsp->sourceQualTags, source_qual_prefixes [i]);
  }
}

static void ValidateSourceQualTags (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, CharPtr str)

{
  Char        ch;
  CharPtr     hit;
  Boolean     okay;
  CharPtr     ptr;
  CharPtr     tmp;
  Int4        state;
  ValNodePtr  matches;

  if (vsp->sourceQualTags == NULL || StringHasNoText (str)) return;
  state = 0;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (vsp->sourceQualTags, state, ch, &matches);
    if (matches != NULL) {
      hit = (CharPtr) matches->data.ptrvalue;
      if (StringHasNoText (hit)) {
        hit = "?";
      }
      okay = TRUE;
      tmp = ptr - StringLen (hit);
      if (tmp > str) {
        ch = *tmp;
        if ((! IS_WHITESP (ch)) && ch != ';') {
          okay = FALSE;
        }
      }
      if (okay) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_StructuredSourceNote,
                  "Source note has structured tag '%s'", hit);
      }
    }
    ptr++;
    ch = *ptr;
  }
}


static CharPtr GetOrgModWarning (Uint2 subtype)
{
  CharPtr warning = NULL;

  switch (subtype) {
    /*
    case ORGMOD_biovar:
      warning = "Biovar value specified is not found in taxname";
      break;
    */
    case ORGMOD_forma:
      warning = "Forma value specified is not found in taxname";
      break;
    case ORGMOD_forma_specialis:
      warning = "Forma specialis value specified is not found in taxname";
      break;
    /*
    case ORGMOD_pathovar:
      warning = "Pathovar value specified is not found in taxname";
      break;
    */
    case ORGMOD_sub_species:
      warning = "Subspecies value specified is not found in taxname";
      break;
    case ORGMOD_variety:
      warning = "Variety value specified is not found in taxname";
      break;
  }
  return warning;
}


static Boolean ValidateOrgModInTaxName (ValidStructPtr vsp, OrgModPtr mod, CharPtr taxname, Boolean varietyOK)
{
  CharPtr cp, f, warn;
  Int4    word_len, name_len;

  if (vsp == NULL || mod == NULL) return FALSE;

  name_len = StringLen (mod->subname);

  /* skip first word */
  word_len = StringCSpn (taxname, " ");
  cp = taxname + word_len;
  cp += StringSpn (cp, " ");
  /* skip second word */
  word_len = StringCSpn (cp, " ");
  cp += word_len;
  cp += StringSpn (cp, " ");  

  f = StringSearch (cp, mod->subname);
  while (f != NULL && ((f != cp && isalpha (*(f - 1))) || isalpha (*(f + name_len)))) {
    f = StringSearch (f + 1, mod->subname);
  }
  if (f == NULL) {
    warn = GetOrgModWarning (mod->subtype);
    if (warn != NULL) {
      /* variety is sorted before sub_species, so if variety was okay in taxname, can ignore missing sub_species */
      if (mod->subtype == ORGMOD_sub_species && varietyOK) return FALSE;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, warn);
      return FALSE;
    }
  }

  return TRUE;
}

/* institution:collection is now stored as a ValNode list of strings, sorted and indexed */

static Boolean inst_code_not_found = FALSE;

static ValNodePtr    ic_code_list = NULL;
static CharPtr PNTR  ic_code_data = NULL;
static Uint1 PNTR    ic_code_type = NULL;
static Int4          ic_code_len = 0;

static void SetupInstCollTable (void)

{
  Char        ch;
  FileCache   fc;
  CharPtr     file = "institution_codes.txt";
  FILE        *fp = NULL;
  Int4        i;
  ValNodePtr  last = NULL;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (ic_code_data != NULL) return;
  if (inst_code_not_found) return;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }  
  
  if (fp == NULL) {
    inst_code_not_found = TRUE;
    return;
  }
  
  FileCacheSetup (&fc, fp);
  
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL) {
    if (StringDoesHaveText (str)) {
      ch = '\0';
      ptr = StringChr (str, '\t');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
        tmp = StringChr (ptr, '\t');
        if (tmp != NULL) {
          *tmp = '\0';
          if (ptr [1] == '\0') {
            ch = *ptr;
          }
        }
      }
      TrimSpacesAroundString (str);
      vnp = ValNodeCopyStr (&last, (Uint1) ch, str);
      if (ic_code_list == NULL) {
        ic_code_list = vnp;
      }
      last = vnp;
    }
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  FileClose (fp);

  ic_code_len = ValNodeLen (ic_code_list);
  if (ic_code_len > 0) {
    ic_code_list = ValNodeSort (ic_code_list, SortVnpByString);
    ic_code_data = (CharPtr PNTR) MemNew (sizeof (CharPtr) * (ic_code_len + 1));
    if (ic_code_data != NULL) {
      for (vnp = ic_code_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        str = (CharPtr) vnp->data.ptrvalue;
        ic_code_data [i] = str;
      }
    }

    ic_code_type = (Uint1 PNTR) MemNew (sizeof (Uint1) * (ic_code_len + 1));
    if (ic_code_type != NULL) {
      for (vnp = ic_code_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        ic_code_type [i] = vnp->choice;
      }
    }
  }
}

static CharPtr CheckInstCollName (CharPtr name, Uint1Ptr typeP)

{
  Int4     L, R, mid;
  CharPtr  str;

  SetupInstCollTable ();

  if (typeP != NULL) {
    *typeP = 0;
  }

  L = 0;
  R = ic_code_len - 1;
  while (L < R) {
    mid = (L + R) / 2;
    str = ic_code_data [(int) mid];
    if (StringICmp (str, name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  if (R < 0) return NULL;

  if (typeP != NULL) {
    *typeP = ic_code_type [(int) R];
  }

  return ic_code_data [(int) R];
}

static void ValidateOrgModVoucher (ValidStructPtr vsp, OrgModPtr mod)

{
  Char     buf [512];
  CharPtr  inst = NULL, id = NULL, coll = NULL, str;
  size_t   len1, len2;
  Uint1    type;
  
  if (vsp == NULL || mod == NULL) return;
  
  StringNCpy_0 (buf, mod->subname, sizeof (buf));
  if (StringChr (buf, ':') == NULL) {
    if (mod->subtype == ORGMOD_culture_collection) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnstructuredVoucher, "Culture_collection should be structured, but is not");
    }
    return;
  }
  if (! ParseStructuredVoucher (buf, &inst, &id)) {
    if (StringHasNoText (inst) || inst [0] == ':') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Voucher is missing institution code");
    }
    if (StringHasNoText (id)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadVoucherID, "Voucher is missing specific identifier");
    }
    return;
  }
  if (inst == NULL) return;
  
  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) {
    if ((mod->subtype == ORGMOD_bio_material && type != 'b') ||
        (mod->subtype == ORGMOD_culture_collection && type != 'c') ||
        (mod->subtype == ORGMOD_specimen_voucher && type != 's')) {
      if (type == 'b') {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be bio_material", inst);
      } else if (type == 'c') {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be culture_collection", inst);
      } else if (type == 's') {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_WrongVoucherType, "Institution code %s should be specimen_voucher", inst);
      }
    }
    return;
  }

  if (StringICmp (str, inst) == 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionCode,
              "Institution code %s exists, but correct capitalization is %s", inst, str);
    return;
  }

  /* ignore personal collections */
  if (StringNICmp (inst, "personal", 8) == 0) return;

  len1 = StringLen (inst);
  len2 = StringLen (str);

  if (len1 < len2) {
    if (StringNICmp (str, inst, len1) == 0 && str [len1] == '<') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s needs to be qualified with a <COUNTRY> designation", inst);
      return;
    }
  }

  coll = StringChr (inst, ':');
  if (coll == NULL) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s is not in list", inst);
    return;
  }

  *coll = '\0';
  coll++;
  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) {
    if (StringCmp (coll, "DNA") == 0) {
      /* DNA is a valid collection for any institution (using bio_material) */
      if (type != 'b') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_WrongVoucherType, "DNA should be bio_material");
      }
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionCode,
                "Institution code %s exists, but collection %s:%s is not in list", inst, inst, coll);
    }
    return;
  }

  len1 = StringLen (inst);
  len2 = StringLen (str);
  
  if (len1 < len2) {
    if (StringNICmp (str, inst, len1) == 0 && str [len1] == '<') {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code in %s:%s needs to be qualified with a <COUNTRY> designation", inst, coll);
      return;
    }
  }
  
  ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadInstitutionCode, "Institution code %s:%s is not in list", inst, coll);
}

NLM_EXTERN Boolean VoucherInstitutionIsValid (CharPtr inst)

{
  CharPtr  str;
  Uint1    type;

  if (StringHasNoText (inst)) return FALSE;
  
  str = CheckInstCollName (inst, &type);
  if (StringCmp (str, inst) == 0) return TRUE;
  
  return FALSE;
}

/* works on subname copy that it can change */

NLM_EXTERN Boolean ParseStructuredVoucher (
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

static void ValidateLatLon (ValidStructPtr vsp, CharPtr lat_lon)

{
  Boolean format_ok = FALSE, lat_in_range = FALSE, lon_in_range = FALSE;
  CharPtr ptr;
  Char    tmp [128];

  IsCorrectLatLonFormat (lat_lon, &format_ok, &lat_in_range, &lon_in_range);

  if (! format_ok) {
    /* may have comma and then altitude, so just get lat_lon component */
    StringNCpy_0 (tmp, lat_lon, sizeof (tmp));
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      lat_lon = tmp;
      IsCorrectLatLonFormat (tmp, &format_ok, &lat_in_range, &lon_in_range);
      if (format_ok) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonFormat, "lat_lon format has extra text after correct dd.dd N|S ddd.dd E|W format");
      }
    }
  }

  if (!format_ok) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonFormat, "lat_lon format is incorrect - should be dd.dd N|S ddd.dd E|W");
  } else {
    if (!lat_in_range) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonRange, "latitude value is out of range - should be between 90.00 N and 90.00 S");
    }
    if (!lon_in_range) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonRange, "longitude value is out of range - should be between 180.00 E and 180.00 W");
    }
  }
}


static void ValidateLocationForHIV (ValidStructPtr vsp, BioSourcePtr biop, BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  Boolean location_ok = FALSE;
  MolInfoPtr mip;

  if (vsp == NULL || biop == NULL) {
    return;
  }

  if (bsp != NULL) {
    if (bsp->mol == Seq_mol_dna) {
      if (biop->genome != GENOME_proviral) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "HIV with moltype DNA should be proviral");
      }
    } else if (bsp->mol == Seq_mol_rna) {
      if (biop->genome == GENOME_unknown) {
        location_ok = TRUE;
      } else if (biop->genome == GENOME_genomic) {
        sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
        if (sdp != NULL 
            && (mip = (MolInfoPtr) sdp->data.ptrvalue) != NULL
            && mip->biomol == MOLECULE_TYPE_GENOMIC) {
          location_ok = TRUE;
            }
      }
      if (!location_ok) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "HIV with moltype RNA should have source location unset or set to genomic (on genomic RNA sequence)");
      }
    }
  }
}
/*****************************************************************************
*
*   ValidateSeqDescrContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static Boolean UnbalancedParentheses (CharPtr str)

{
  Char  ch;
  Int4  fwd_par = 0, rev_par = 0, fwd_bkt = 0, rev_bkt = 0;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch == '(') {
      fwd_par++;
    } else if (ch == ')') {
      rev_par++;
    } else if (ch == '[') {
      fwd_bkt++;
    } else if (ch == ']') {
      rev_bkt++;
    }
    str++;
    ch = *str;
  }

  if (fwd_par != rev_par) return TRUE;
  if (fwd_bkt != rev_bkt) return TRUE;

  return FALSE;
}

static CharPtr valid_sex_values [] = {
  "female",
  "male",
  "hermaphrodite",
  "unisexual",
  "bisexual",
  "asexual",
  "monoecious",
  "monecious",
  "dioecious",
  "diecious",
  NULL
};

static Boolean IsValidSexValue (CharPtr str)

{
  int  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; valid_sex_values [i] != NULL; i++) {
    if (StringICmp (str, valid_sex_values [i]) == 0) return TRUE;
  }

  return FALSE;
}

static void ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop, SeqFeatPtr sfp, ValNodePtr sdp)
{
  Char            badch;
  Boolean         bad_cap = FALSE;
  Boolean         bad_frequency;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Char            buf [256];
  Char            ch;
  Boolean         chromconf = FALSE;
  Int2            chromcount = 0;
  SubSourcePtr    chromosome = NULL;
  CharPtr         countryname = NULL;
  ValNodePtr      db;
  DbtagPtr        dbt;
  Boolean         format_ok;
  CharPtr         gb_synonym = NULL;
  Boolean         germline = FALSE;
  CharPtr         good;
  CharPtr         guess = NULL;
  Boolean         has_strain = FALSE;
  Boolean         has_fwd_pcr_seq = FALSE;
  Boolean         has_rev_pcr_seq = FALSE;
  Boolean         has_pcr_name = FALSE;
  Boolean         has_metagenome_source = FALSE;
  Int4            id;
  Boolean         is_env_sample = FALSE;
  Boolean         is_iso_source = FALSE;
  Boolean         is_mating_type = FALSE;
  Boolean         is_metagenomic = FALSE;
  Boolean         is_sex = FALSE;
  Boolean         is_specific_host = FALSE;
  Boolean         is_transgenic = FALSE;
  Boolean         isAnimal = FALSE;
  Boolean         isArchaea = FALSE;
  Boolean         isBacteria = FALSE;
  Boolean         isFungal = FALSE;
  Boolean         isPlant = FALSE;
  Boolean         isViral = FALSE;
  Boolean         is_bc;
  Boolean         is_rf;
  Boolean         is_sc;
  CharPtr         last_db = NULL;
  FloatHi         lat = 0.0;
  FloatHi         lon = 0.0;
  CharPtr         lat_lon = NULL;
  Boolean         lat_in_range;
  Boolean         lon_in_range;
  Int2            num_bio_material = 0;
  Int2            num_culture_collection = 0;
  Int2            num_specimen_voucher = 0;
  Int2            num_country = 0;
  Int2            num_lat_lon = 0;
  Int2            num_fwd_primer_seq = 0;
  Int2            num_rev_primer_seq = 0;
  Int2            num_fwd_primer_name = 0;
  Int2            num_rev_primer_name = 0;
  Boolean         old_country = FALSE;
  OrgNamePtr      onp;
  OrgModPtr       omp, nxtomp;
  OrgRefPtr       orp;
  ObjValNodePtr   ovp;
  Int4            primer_len_before;
  Int4            primer_len_after;
  ValNodePtr      pset;
  CharPtr         ptr;
  Boolean         rearranged = FALSE;
  SeqEntryPtr     sep;
  ErrSev          sev;
  SubSourcePtr    ssp;
  CharPtr         str;
  Boolean         strict = TRUE;
  CharPtr         synonym = NULL;
  Char            tmp [128];
  Boolean         varietyOK;
  CharPtr         inst1, inst2, id1, id2, coll1, coll2;
  Char            buf1 [512], buf2 [512];

  if (vsp->sourceQualTags == NULL) {
    InitializeSourceQualTags (vsp);
  }
  if (biop == NULL)
    return;
  if (biop->genome == GENOME_transposon || biop->genome == GENOME_insertion_seq) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_ObsoleteSourceLocation,
              "Transposon and insertion sequence are no longer legal locations");
  }

  if (vsp->indexerVersion && biop->genome == GENOME_chromosome) {
    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_ChromosomeLocation, "INDEXER_ONLY - BioSource location is chromosome");
  }

  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) {
        isViral = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Metazoa; ", 20) == 0) {
        isAnimal = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; ", 53) == 0 ||
                 StringNICmp (onp->lineage, "Eukaryota; Rhodophyta; ", 23) == 0 ||
                 StringNICmp (onp->lineage, "Eukaryota; stramenopiles; Phaeophyceae; ", 40) == 0) {
        isPlant = TRUE;
      } else if (StringNICmp (onp->lineage, "Bacteria; ", 10) == 0) {
        isBacteria = TRUE;
      } else if (StringNICmp (onp->lineage, "Archaea; ", 9) == 0) {
        isArchaea = TRUE;
      } else if (StringNICmp (onp->lineage, "Eukaryota; Fungi; ", 18) == 0) {
        isFungal = TRUE;
      }
    }
  }

  ssp = biop->subtype;
  while (ssp != NULL) {
    if (ssp->subtype == SUBSRC_country) {
      num_country++;
      if (countryname != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCountryCode, "Multiple country names on BioSource");
      }
      countryname = ssp->name;
      if (CountryIsValid (countryname, &old_country, &bad_cap)) {
        if (bad_cap) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCountryCapitalization, "Bad country capitalization [%s]", countryname);
        }
      } else {
        if (StringHasNoText (countryname)) {
          countryname = "?";
        }
        if (old_country) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_ReplacedCountryCode, "Replaced country name [%s]", countryname);
        } else {
          sev = SEV_ERROR;
          if (vsp->is_barcode_sep && vsp->seqSubmitParent) {
            sev = SEV_WARNING;
          }
          ValidErr (vsp, sev, ERR_SEQ_DESCR_BadCountryCode, "Bad country name [%s]", countryname);
        }
      }
    } else if (ssp->subtype == SUBSRC_chromosome) {
      chromcount++;
      if (chromosome != NULL) {
        if (StringICmp (ssp->name, chromosome->name) != 0) {
          chromconf = TRUE;
        }
      } else {
        chromosome = ssp;
      }
    } else if (ssp->subtype == SUBSRC_transposon_name || ssp->subtype == SUBSRC_insertion_seq_name) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_ObsoleteSourceQual,
                "Transposon name and insertion sequence name are no longer legal qualifiers");
    } else if (ssp->subtype == 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadSubSource, "Unknown subsource subtype %d", (int) (ssp->subtype));
    } else if (ssp->subtype == SUBSRC_other) {
      ValidateSourceQualTags (vsp, gcp, biop, ssp->name);
    } else if (ssp->subtype == SUBSRC_germline) {
      germline = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Germline qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_rearranged) {
      rearranged = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Rearranged qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_transgenic) {
      is_transgenic = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Transgenic qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_environmental_sample) {
      is_env_sample = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental_sample qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_metagenomic) {
      is_metagenomic = TRUE;
      str = ssp->name;
      if (str == NULL || str [0] != '\0') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenomic qualifier should not have descriptive text");
      }
    } else if (ssp->subtype == SUBSRC_isolation_source) {
      is_iso_source = TRUE;
    } else if (ssp->subtype == SUBSRC_sex) {
      is_sex = TRUE;
      str = ssp->name;
      if (isAnimal || isPlant) {
        /* always use /sex, do not check values at this time */
      } else if (isViral || isBacteria || isArchaea || isFungal) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /sex qualifier");
      } else if (IsValidSexValue (str)) {
        /* otherwise expect male or female, or a few others */
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /sex qualifier");
      }
    } else if (ssp->subtype == SUBSRC_mating_type) {
      is_mating_type = TRUE;
      str = ssp->name;
      if (isAnimal || isPlant || isViral) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /mating_type qualifier");
      } else if (IsValidSexValue (str)) {
        /* complain if one of the values that should go in /sex */
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Unexpected use of /mating_type qualifier");
      }
    } else if (ssp->subtype == SUBSRC_plasmid_name) {
      if (biop->genome != GENOME_plasmid) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plasmid subsource but not plasmid location");
      }
    } else if (ssp->subtype == SUBSRC_plastid_name) {
      if (StringCmp (ssp->name, "chloroplast") == 0) {
        if (biop->genome != GENOME_chloroplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource chloroplast but not chloroplast location");
        }
      } else if (StringCmp (ssp->name, "chromoplast") == 0) {
        if (biop->genome != GENOME_chromoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource chromoplast but not chromoplast location");
        }
      } else if (StringCmp (ssp->name, "kinetoplast") == 0) {
        if (biop->genome != GENOME_kinetoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource kinetoplast but not kinetoplast location");
        }
      } else if (StringCmp (ssp->name, "plastid") == 0) {    
        if (biop->genome != GENOME_plastid) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource plastid but not plastid location");
        }
      } else if (StringCmp (ssp->name, "apicoplast") == 0) {  
        if (biop->genome != GENOME_apicoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource apicoplast but not apicoplast location");
        }
      } else if (StringCmp (ssp->name, "leucoplast") == 0) {  
        if (biop->genome != GENOME_leucoplast) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource leucoplast but not leucoplast location");
        }
      } else if (StringCmp (ssp->name, "proplastid") == 0) {  
        if (biop->genome != GENOME_proplastid) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource proplastid but not proplastid location");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Plastid name subsource contains unrecognized value");
      }  
    } else if (ssp->subtype == SUBSRC_collection_date) {
      if (! CollectionDateIsValid (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionDate, "Collection_date format is not in DD-Mmm-YYYY format");
      }
      else if (CollectionDateIsInTheFuture (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCollectionDate, "Collection_date is in the future");
      }
    } else if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      num_fwd_primer_seq++;
      has_fwd_pcr_seq = TRUE;
      if (! PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR forward primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      /*
      if (PrimerSeqHasDuplicates (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                  "PCR forward primer sequence has duplicates");
      }
      */
    } else if (ssp->subtype == SUBSRC_rev_primer_seq) {
      num_rev_primer_seq++;
      has_rev_pcr_seq = TRUE;
      if (! PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence,
                  "PCR reverse primer sequence format is incorrect, first bad character is '%c'", (char) badch);
      }
      /*
      if (PrimerSeqHasDuplicates (ssp->name)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                  "PCR reverse primer sequence has duplicates");
      }
      */
    } else if (ssp->subtype == SUBSRC_fwd_primer_name) {
      num_fwd_primer_name++;
      if (StringLen (ssp->name) > 10 && PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR primer name appears to be a sequence");
      }
      has_pcr_name = TRUE;
    } else if (ssp->subtype == SUBSRC_rev_primer_name) {
      num_rev_primer_name++;
      if (StringLen (ssp->name) > 10 && PrimerSeqIsValid (vsp, ssp->name, &badch)) {
        if (badch < ' ' || badch > '~') {
          badch = '?';
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerName, "PCR primer name appears to be a sequence");
      }
      has_pcr_name = TRUE;
    } else if (ssp->subtype == SUBSRC_lat_lon) {
      num_lat_lon++;
      if (lat_lon != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonProblem, "Multiple lat_lon on BioSource");
      }
      lat_lon = ssp->name;
      ValidateLatLon (vsp, lat_lon);
    } else if (ssp->subtype == SUBSRC_frequency) {
      str = ssp->name;
      if (StringDoesHaveText (str)) {
        bad_frequency = FALSE;
        if (StringCmp (str, "0") == 0) {
          /* ignore */
        } else if (StringCmp (str, "1") == 0) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_BioSourceInconsistency, "bad frequency qualifier value %s", ssp->name);
        } else {
          ch = *str;
          if (ch == '0') {
            str++;
            ch = *str;
          }
          if (ch == '.') {
            str++;
            ch = *str;
            if (! IS_DIGIT (ch)) {
              bad_frequency = TRUE;
            } else {
              while (ch != '\0') {
                if (! IS_DIGIT (ch)) {
                  bad_frequency = TRUE;
                }
                str++;
                ch = *str;
              }
            }
          } else {
            bad_frequency = TRUE;
          }
        }
        if (bad_frequency) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "bad frequency qualifier value %s", ssp->name);
        }
      }
    } else if (ssp->subtype == SUBSRC_sex && isViral) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected sex qualifier");
    } else if (ssp->subtype == SUBSRC_cell_line && isViral) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected cell_line qualifier");
    } else if (ssp->subtype == SUBSRC_cell_type && isViral) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected cell_type qualifier");
    } else if (ssp->subtype == SUBSRC_tissue_type && isViral) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Virus has unexpected tissue_type qualifier");
    }
    ssp = ssp->next;
  }
  if (num_country > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple country qualifiers present");
  }
  if (num_lat_lon > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple lat_lon qualifiers present");
  }
  if (num_fwd_primer_seq > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple fwd_primer_seq qualifiers present");
  }
  if (num_rev_primer_seq > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple rev_primer_seq qualifiers present");
  }
  if (num_fwd_primer_name > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple fwd_primer_name qualifiers present");
  }
  if (num_rev_primer_name > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple rev_primer_name qualifiers present");
  }

  if (countryname != NULL && lat_lon != NULL) {
    IsCorrectLatLonFormat (lat_lon, &format_ok, &lat_in_range, &lon_in_range);
    if (! format_ok) {
      /* may have comma and then altitude, so just get lat_lon component */
      StringNCpy_0 (tmp, lat_lon, sizeof (tmp));
      ptr = StringChr (tmp, ',');
      if (ptr != NULL) {
        *ptr = '\0';
        lat_lon = tmp;
        IsCorrectLatLonFormat (tmp, &format_ok, &lat_in_range, &lon_in_range);
      }
    }
    if (format_ok && ParseLatLon (lat_lon, &lat, &lon)) {
      StringNCpy_0 (buf, countryname, sizeof (buf));
      ptr = StringChr (buf, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        strict = FALSE;
      }
      if (IsCountryInLatLonList (buf)) {
        if (TestLatLonForCountry (buf, lat, lon)) {
          /* match */
          if (! strict) {
            StringNCpy_0 (buf, countryname, sizeof (buf));
            ptr = StringChr (buf, ',');
            if (ptr != NULL) {
              *ptr = '\0';
            }
            ptr = StringChr (buf, ';');
            if (ptr != NULL) {
              *ptr = '\0';
            }
            if (IsCountryInLatLonList (buf)) {
              if (TestLatLonForCountry (buf, lat, lon)) {
                /* match */
              } else {
                if (vsp->strictLatLonCountry || (vsp->testLatLonSubregion && (! StringContainsBodyOfWater (countryname)))) {
                  /* passed unqualified but failed qualified country name, report at info level for now */
                  guess = GuessCountryForLatLon (lat, lon);
                  if (StringDoesHaveText (guess)) {
                    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonState,
                              "Lat_lon '%s' does not map to subregion '%s', but may be in '%s'", lat_lon, buf, guess);
                  } else {
                    ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonState,
                              "Lat_lon '%s' does not map to subregion '%s'", lat_lon, buf);
                  }
                }
              }
            }
          }
        } else if (TestLatLonForCountry (buf, -lat, lon)) {
          if (lat < 0.0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude should be set to N (northern hemisphere)");
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude should be set to S (southern hemisphere)");
          }
        } else if (TestLatLonForCountry (buf, lat, -lon)) {
          if (lon < 0.0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Longitude should be set to E (eastern hemisphere)");
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Longitude should be set to W (western hemisphere)");
          }
        /*
        } else if (TestLatLonForCountry (buf, -lat, -lon)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Both latitude and longitude appear to be in wrong hemispheres");
        */
        } else if (TestLatLonForCountry (buf, lon, lat)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_LatLonValue, "Latitude and longitude values appear to be exchanged");
        /*
        } else if (strict) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry, "Lat_lon '%s' does not map to '%s'", lat_lon, buf);
        */
        } else {
          if (vsp->strictLatLonCountry || (! StringContainsBodyOfWater (countryname))) {
            guess = GuessCountryForLatLon (lat, lon);
            if (guess != NULL) {
              ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry,
                        "Lat_lon '%s' does not map to '%s', but may be in '%s'", lat_lon, buf, guess);
            } else {
              ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_LatLonCountry,
                        "Lat_lon '%s' does not map to '%s'", lat_lon, buf);
            }
          }
        }
      }
    }
  }

  if (has_pcr_name) {
    if ((! has_fwd_pcr_seq) || (! has_rev_pcr_seq)) {
      /*
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence, "PCR primer has name but not both sequences");
      */
    }
  } else if (has_fwd_pcr_seq || has_rev_pcr_seq) {
    if (! (has_fwd_pcr_seq && has_rev_pcr_seq)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPCRPrimerSequence, "PCR primer does not have both sequences");
    }
  }

  pset = ParsePCRSet (biop);
  if (pset != NULL) {
    pset = ValNodeSort (pset, SortVnpByPCRSetSeq);
    primer_len_before = ValNodeLen (pset);
    pset = UniqueVnpByPCRSetSeq (pset);
    primer_len_after = ValNodeLen (pset);
    FreePCRSet (pset);
    if (primer_len_before != primer_len_after) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_DuplicatePCRPrimerSequence,
                "PCR primer sequence has duplicates");
    }
  }

  if (germline && rearranged) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Germline and rearranged should not both be present");
  }
  if (is_transgenic && is_env_sample) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Transgenic and environmental sample should not both be present");
  }
  if (is_metagenomic && (! is_env_sample)) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenomic should also have environmental sample annotated");
  }
  if (is_sex && is_mating_type) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Sex and mating type should not both be present");
  }

  if (biop->org != NULL 
      && biop->org->orgname != NULL
      && StringISearch (biop->org->orgname->lineage, "metagenomes") != NULL
      && !is_metagenomic) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "If metagenomes appears in lineage, BioSource should have metagenomic qualifier");
  }
  if (chromcount > 1) {
    if (chromconf) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleChromosomes, "Multiple conflicting chromosome qualifiers");
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleChromosomes, "Multiple identical chromosome qualifiers");
    }
  }
  orp = biop->org;
  if (orp != NULL) {
    /*
    if (StringICmp (orp->taxname, "Human immunodeficiency virus") == 0 ||
        StringICmp (orp->taxname, "Human immunodeficiency virus 1") == 0 ||
        StringICmp (orp->taxname, "Human immunodeficiency virus 2") == 0) {
      ValidateLocationForHIV (vsp, biop);
    } else */
    if (StringICmp (orp->taxname, "uncultured bacterium") == 0) {
      bsp = NULL;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
      } else if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
        } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) ovp->idx.parentptr;
          if (bssp != NULL) {
            sep = bssp->seqentry;
            if (sep != NULL) {
              sep = FindNthBioseq (sep, 1);
              if (sep != NULL && IS_Bioseq (sep)) {
                bsp = (BioseqPtr) sep->data.ptrvalue;
              }
            }
          }
        }
      }
      if (bsp != NULL && bsp->length >= 10000) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Uncultured bacterium sequence length is suspiciously high");
      }
    }
    if (StringNICmp (orp->taxname, "uncultured ", 11) == 0) {
      if (! is_env_sample) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Uncultured should also have /environmental_sample");
      }
    }
  }

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_germline ||
              ssp->subtype == SUBSRC_rearranged ||
              ssp->subtype == SUBSRC_transgenic ||
              ssp->subtype == SUBSRC_environmental_sample ||
              ssp->subtype == SUBSRC_metagenomic) continue;
    str = ssp->name;
    if (StringHasNoText (str)) continue;
    if (UnbalancedParentheses (str)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_UnbalancedParentheses,
                  "Unbalanced parentheses in '%s'", str);
    }
  }

  if (orp == NULL || (StringHasNoText (orp->taxname) && StringHasNoText (orp->common))) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name has been applied to this Bioseq.  Other qualifiers may exist.");
  }
  if (orp == NULL) {
    if (is_env_sample && (! is_iso_source) && (! is_specific_host)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental sample should also have isolation source or specific host annotated");
    }
    return;
  }
  onp = orp->orgname;
  if (onp == NULL || StringHasNoText (onp->lineage)) {
    if (! vsp->seqSubmitParent) { /* suppress when validator run from tbl2asn */
      sev = SEV_ERROR;
      if (vsp->is_refseq_in_sep) {
        for (db = orp->db; db != NULL; db = db->next) {
          dbt = (DbtagPtr) db->data.ptrvalue;
          if (dbt != NULL) {
            if (StringICmp (dbt->db, "taxon") == 0) {
              sev = SEV_REJECT;
            }
          }
        }
      }
      if (vsp->is_embl_ddbj_in_sep) {
        sev = SEV_WARNING;
      }
      ValidErr (vsp, sev, ERR_SEQ_DESCR_MissingLineage, "No lineage for this BioSource.");
    }
  } else {
    if (biop->genome == GENOME_kinetoplast) {
      if (StringStr (onp->lineage, "Kinetoplastida") == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrganelle, "Only Kinetoplastida have kinetoplasts");
      }
    } else if (biop->genome == GENOME_nucleomorph) {
      if (StringStr (onp->lineage, "Chlorarachniophyceae") == 0 && StringStr (onp->lineage, "Cryptophyta") == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrganelle, "Only Chlorarachniophyceae and Cryptophyta have nucleomorphs");
      }
    }

    /* warn if bacteria has organelle location */
    if (StringCmp (onp->div, "BCT") == 0 
        && biop->genome != GENOME_unknown 
        && biop->genome != GENOME_genomic
        && biop->genome != GENOME_plasmid) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Bacterial source should not have organelle location");
    }

    if (StringCmp (onp->div, "ENV") == 0 && (! is_env_sample)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "BioSource with ENV division is missing environmental sample subsource");
    }
  }
  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL) {
      if (last_db != NULL) {
        if (StringICmp (dbt->db, last_db) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceDbTagConflict, "BioSource uses db %s multiple times", last_db);
        }
      }
      last_db = dbt->db;
    }
  }
  if (onp != NULL) {
    omp = onp->mod;
    varietyOK = FALSE;
    while (omp != NULL) {
      if (omp->subtype == 0 || omp->subtype == 1) {
        ValidErr (vsp, SEV_REJECT, ERR_SEQ_DESCR_BadOrgMod, "Unknown orgmod subtype %d", (int) (omp->subtype));
      } else if (omp->subtype == ORGMOD_strain) {
        if (has_strain) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Multiple strain qualifiers on the same BioSource");
        }
        has_strain = TRUE;
      } else if (omp->subtype == ORGMOD_variety) {
        if ((StringHasNoText (onp->div) || StringICmp (onp->div, "PLN") != 0) &&
            StringStr (onp->lineage, "Cyanobacteria") == 0 &&
            StringStr (onp->lineage, "Myxogastria") == 0 &&
            StringStr (onp->lineage, "Oomycetes") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Orgmod variety should only be in plants, fungi, or cyanobacteria");
        }
        varietyOK = ValidateOrgModInTaxName (vsp, omp, orp->taxname, varietyOK);
      } else if (omp->subtype == ORGMOD_nat_host) {
        is_specific_host = TRUE;
        if (StringICmp (omp->subname, orp->taxname) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "Specific host is identical to taxname");
        }
      } else if (omp->subtype == ORGMOD_other) {
        ValidateSourceQualTags (vsp, gcp, biop, omp->subname);
      } else if (omp->subtype == ORGMOD_biovar ||
                 omp->subtype == ORGMOD_forma ||
                 omp->subtype == ORGMOD_forma_specialis ||
                 omp->subtype == ORGMOD_sub_species ||
                 omp->subtype == ORGMOD_pathovar) {         
        ValidateOrgModInTaxName (vsp, omp, orp->taxname, varietyOK);
      } else if (omp->subtype == ORGMOD_specimen_voucher) {
        num_specimen_voucher++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_culture_collection) {
        num_culture_collection++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_bio_material) {
        num_bio_material++;
        ValidateOrgModVoucher (vsp, omp);
      } else if (omp->subtype == ORGMOD_metagenome_source) {
        has_metagenome_source = TRUE;
      } else if (omp->subtype == ORGMOD_common) {
        if (StringICmp (omp->subname, orp->common) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadOrgMod, "OrgMod common is identical to Org-ref common");
        }
      } else if (omp->subtype == ORGMOD_synonym) {
        synonym = omp->subname;
      } else if (omp->subtype == ORGMOD_gb_synonym) {
        gb_synonym = omp->subname;
      }
      omp = omp->next;
    }

    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      if (omp->subtype != ORGMOD_specimen_voucher &&
          omp->subtype != ORGMOD_culture_collection &&
          omp->subtype != ORGMOD_bio_material) continue;
      nxtomp = omp->next;
      if (nxtomp == NULL) continue;
      inst1 = NULL;
      inst2 = NULL;
      id1 = NULL;
      id2 = NULL;
      coll1 = NULL;
      coll2 = NULL;
      StringNCpy_0 (buf1, omp->subname, sizeof (buf1));
      StringNCpy_0 (buf2, nxtomp->subname, sizeof (buf2));
      if (StringChr (buf1, ':') == NULL || StringChr (buf2, ':') == NULL) continue;
      if (! ParseStructuredVoucher (buf1, &inst1, &id1)) continue;
      if (! ParseStructuredVoucher (buf2, &inst2, &id2)) continue;
      if (inst1 == NULL || inst2 == NULL) continue;
      if (StringNICmp (inst1, "personal", 8) == 0) continue;
      if (StringNICmp (inst2, "personal", 8) == 0) continue;
      coll1 = StringChr (inst1, ':');
      if (coll1 != NULL) {
        *coll1 = '\0';
        coll1++;
      }
      coll2 = StringChr (inst2, ':');
      if (coll2 != NULL) {
        *coll2 = '\0';
        coll2++;
      }
      if (StringICmp (inst1, inst2) != 0) continue;
      if (StringCmp (coll1, "DNA") == 0 || StringCmp (coll2, "DNA") == 0) continue;
      if (coll1 != NULL && coll2 != NULL && StringICmp (coll1, coll2) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceVouchers, "Multiple vouchers with same institution:collection");
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceVouchers, "Multiple vouchers with same institution");
      }
    }
  }

  if (onp != NULL) {
    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      str = omp->subname;
      if (StringHasNoText (str)) continue;
      if (UnbalancedParentheses (str)) {
        if (omp->subtype == ORGMOD_old_name) continue;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_UnbalancedParentheses,
                  "Unbalanced parentheses in '%s'", str);
      }
    }
  }

  /*
  if (num_bio_material > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple bio_material qualifiers present");
  }
  if (num_culture_collection > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple culture_collection qualifiers present");
  }
  if (num_specimen_voucher > 1) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleSourceQualifiers, "Multiple specimen_voucher qualifiers present");
  }
  */
  if (is_env_sample && has_strain) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "Strain should not be present in an environmental sample");
  }
  if (is_env_sample && (! is_iso_source) && (! is_specific_host)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Environmental sample should also have isolation source or specific host annotated");
  }
  if (has_metagenome_source && (! is_metagenomic)) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceInconsistency, "Metagenome source should also have metagenomic qualifier");
  }
  if (StringDoesHaveText (synonym) && StringDoesHaveText (gb_synonym)) {
    if (StringICmp (synonym, gb_synonym) == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "OrgMod synonym is identical to OrgMod gb_synonym");
    }
  }

  for (db = orp->db; db != NULL; db = db->next) {
    id = -1;
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL && dbt->db != NULL) {

      if (DbxrefIsValid (dbt->db, &is_rf, &is_sc, &is_bc, &good)) {
        if (is_bc) {
          if (StringHasNoText (good)) {
            good = "?";
          }
          if (is_sc) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s", dbt->db, good);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s, but should not be used on an OrgRef",
                      dbt->db, good);
          }
        } else if (is_rf) {
          if (vsp->is_refseq_in_sep || vsp->is_gps_in_sep) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, 
                      "RefSeq-specific db_xref type %s should not used on an OrgRef", dbt->db);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "RefSeq-specific db_xref type %s should not used on a non-RefSeq OrgRef", dbt->db);
          }
        } else if (is_sc) {
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                    "db_xref type %s should not used on an OrgRef", dbt->db);
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", dbt->db);
      }

      /*
      dbxerr = NULL;
      dbvalid = IsDbxrefValid (dbt->db, NULL, orp, FALSE, &dbxerr);
      if (dbxerr != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, dbxerr);
        dbxerr = MemFree (dbxerr);
      }
      */
    }
  }

  if (GetAppProperty ("InternalNcbiSequin") == NULL) return;

  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL) {
      if (StringICmp (dbt->db, "taxon") == 0)
        return;
    }
  }
  if (! vsp->seqSubmitParent) { /* suppress when validator run from tbl2asn */
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_NoTaxonID, "BioSource is missing taxon ID");
  }
}

static Boolean IsXr (ValNodePtr sdp)

{
  BioseqPtr      bsp;
  ObjValNodePtr  ovp;
  SeqIdPtr       sip;
  TextSeqIdPtr   tsip;

  if (sdp->extended == 0) return FALSE;
  ovp = (ObjValNodePtr) sdp;
  if (ovp->idx.parenttype != OBJ_BIOSEQ) return FALSE;
  bsp = (BioseqPtr) ovp->idx.parentptr;
  if (bsp == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice != SEQID_OTHER) continue;
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip == NULL) continue;
    if (StringNICmp (tsip->accession, "XR_", 3) == 0) return TRUE;
  }
  return FALSE;
}

static Boolean IsSynthetic (BioseqPtr bsp)

{
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  if (biop->origin == 5) return TRUE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;
  if (StringICmp (onp->div, "SYN") == 0) return TRUE;
  return FALSE;
}

static Boolean IsMicroRNA (BioseqPtr bsp)

{
  SeqMgrFeatContext  fcontext;
  RnaRefPtr          rrp;
  SeqFeatPtr         sfp;
  CharPtr            str;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_otherRNA, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_RNA) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringStr (str, "microRNA") != NULL) return TRUE;
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_otherRNA, &fcontext);
  }
  return FALSE;
}

static Boolean IsOtherDNA (BioseqPtr bsp)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return FALSE;
  if (mip->biomol == 255) return TRUE;
  return FALSE;
}

static Boolean ValidateSeqDescrCommon (ValNodePtr sdp, BioseqValidStrPtr bvsp, ValidStructPtr vsp, Uint4 descitemid)
{
  ValNodePtr      vnp, vnp2;
  OrgRefPtr       this_org = NULL, that_org = NULL;
  int             tmpval;
  Char            buf1[20], buf2[20], ch;
  EMBLBlockPtr    ebp;
  GBBlockPtr      gbp;
  ValNodePtr      keywords = NULL;
  PubdescPtr      pdp;
  MolInfoPtr      mip;
  ObjectIdPtr     oip;
  Uint2           olditemtype = 0;
  Uint4           olditemid = 0;
  BioSourcePtr    biop;
  GatherContextPtr gcp = NULL;
  CharPtr         str, ptr;
  SeqFeatPtr      sfp;
  Boolean         tpa_exp;
  Boolean         tpa_inf;
  UserObjectPtr   uop;
  BioseqPtr       bsp;
  DatePtr         dp;
  size_t          len;
  SeqMgrFeatContext  fcontext;
  static char    *badmod = "Inconsistent GIBB-mod [%d] and [%d]";

  vsp->sfp = NULL;
  vnp = sdp;
  vsp->descr = vnp;

  if (descitemid > 0) {
    gcp = vsp->gcp;
    if (gcp != NULL) {
      olditemid = gcp->itemID;
      olditemtype = gcp->thistype;
      gcp->itemID = descitemid;
      gcp->thistype = OBJ_SEQDESC;
    }
  }

  switch (vnp->choice) {
  case Seq_descr_mol_type:
    tmpval = (int) (vnp->data.intvalue);
    switch (tmpval) {
    case 8:                    /* peptide */
      if (!bvsp->is_aa)
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with GIBB-mol = peptide");
      break;
    case 0:                    /* unknown */
    case 255:                  /* other */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "GIBB-mol unknown or other used");
      break;
    default:                   /* the rest are nucleic acid */
      if (bvsp->is_aa) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "GIBB-mol [%d] used on protein", tmpval);
      } else {
        if (bvsp->last_na_mol) {
          if (bvsp->last_na_mol != (int) vnp->data.intvalue) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent GIBB-mol [%d] and [%d]", bvsp->last_na_mol, tmpval);
          }
        } else
          bvsp->last_na_mol = tmpval;
      }
      break;
    }
    break;
  case Seq_descr_modif:
    for (vnp2 = (ValNodePtr) (vnp->data.ptrvalue); vnp2 != NULL; vnp2 = vnp2->next) {
      tmpval = (int) (vnp2->data.intvalue);
      switch (tmpval) {
      case 0:                  /* dna */
      case 1:                  /* rna */
        if (bvsp->is_aa) {      /* only temporarily on 0 */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid GIBB-mod [%d] on protein", tmpval);
        } else if (bvsp->last_na_mod) {
          if (tmpval != bvsp->last_na_mod) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_na_mod, tmpval);
          }
        } else
          bvsp->last_na_mod = tmpval;
        break;
      case 4:                  /* mitochondria */
      case 5:                  /* chloroplast */
      case 6:                  /* kinetoplast */
      case 7:                  /* cyanelle */
      case 18:                 /* macronuclear */
        if (bvsp->last_organelle) {
          if (tmpval != bvsp->last_na_mod) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_organelle, tmpval);
          }
        } else
          bvsp->last_organelle = tmpval;
        break;
      case 10:                 /* partial */
      case 11:                 /* complete */
        if (bvsp->last_partialness) {
          if (tmpval != bvsp->last_partialness) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_partialness, tmpval);
          }
        } else
          bvsp->last_partialness = tmpval;
        if ((bvsp->last_left_right) && (tmpval == 11)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_left_right, tmpval);
        }
        break;
      case 16:                 /* no left */
      case 17:                 /* no right */
        if (bvsp->last_partialness == 11) {     /* complete */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod, bvsp->last_partialness, tmpval);
        }
        bvsp->last_left_right = tmpval;
        break;
      case 255:                /* other */
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Unknown, "GIBB-mod = other used");
        break;
      default:
        break;

      }
    }
    break;
  case Seq_descr_method:
    if (!bvsp->is_aa) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
    }
    break;
  case Seq_descr_comment:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Comment descriptor needs text");
    }
    if (SerialNumberInString (str)) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_DESCR_SerialInComment,
                "Comment may refer to reference by serial number - attach reference specific comments to the reference REMARK instead.");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_comment) {
        ptr = (CharPtr) vnp2->data.ptrvalue;
        if (StringDoesHaveText (ptr) && StringICmp (str, ptr) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleComments, "Undesired multiple comment descriptors, identical text");
        }
      }
    }
    break;
  case Seq_descr_genbank:
    if (bvsp->last_gb != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple GenBank blocks");
    else
      bvsp->last_gb = vnp;
    if (vnp != NULL) {
      gbp = (GBBlockPtr) vnp->data.ptrvalue;
      if (gbp != NULL) {
        keywords = gbp->keywords;
      }
    }
    break;
  case Seq_descr_embl:
    if (bvsp->last_embl != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple EMBL blocks");
    else
      bvsp->last_embl = vnp;
    if (vnp != NULL) {
      ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
      if (ebp != NULL) {
        keywords = ebp->keywords;
      }
    }
    break;
  case Seq_descr_pir:
    if (bvsp->last_pir != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PIR blocks");
    else
      bvsp->last_pir = vnp;
    break;
  case Seq_descr_sp:
    if (bvsp->last_sp != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple SWISS-PROT blocks");
    else
      bvsp->last_sp = vnp;
    break;
  case Seq_descr_pdb:
    if (bvsp->last_pdb != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PDB blocks");
    else
      bvsp->last_pdb = vnp;
    break;
  case Seq_descr_prf:
    if (bvsp->last_prf != NULL)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PRF blocks");
    else
      bvsp->last_prf = vnp;
    break;
  case Seq_descr_create_date:
    dp = (DatePtr) vnp->data.ptrvalue;
    if (DateIsBad (dp, TRUE) > 0) {
      ValidErr (vsp, SEV_ERROR, ERR_GENERIC_BadDate, "Create date has error");
    }
    if (bvsp->last_create != NULL) {
      tmpval = (int) DateMatch ((DatePtr) vnp->data.ptrvalue, (DatePtr) (bvsp->last_create->data.ptrvalue), FALSE);
      if (tmpval) {
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (bvsp->last_create->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_dates [%s] and [%s]", buf1, buf2);
      }
    } else
      bvsp->last_create = vnp;
    if (bvsp->last_update != NULL) {
      tmpval = (int) DateMatch ((DatePtr) vnp->data.ptrvalue, (DatePtr) (bvsp->last_update->data.ptrvalue), FALSE);
      if (tmpval == 1) {
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (bvsp->last_update->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]", buf1, buf2);
      }
    }
    break;
  case Seq_descr_update_date:
    dp = (DatePtr) vnp->data.ptrvalue;
    if (DateIsBad (dp, TRUE) > 0) {
      ValidErr (vsp, SEV_ERROR, ERR_GENERIC_BadDate, "Update date has error");
    }
    if (bvsp->last_create != NULL) {
      tmpval = (int) DateMatch ((DatePtr) bvsp->last_create->data.ptrvalue, (DatePtr) (vnp->data.ptrvalue), FALSE);
      if (tmpval == 1) {
        DatePrint ((DatePtr) (bvsp->last_create->data.ptrvalue), buf1);
        DatePrint ((DatePtr) (vnp->data.ptrvalue), buf2);
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]", buf1, buf2);
      }
    }
    if (bvsp->last_update == NULL)
      bvsp->last_update = vnp;
    break;
  case Seq_descr_source:
    biop = (BioSourcePtr) vnp->data.ptrvalue;
    bsp = bvsp->bsp;
    if (biop != NULL && biop->is_focus && bsp != NULL) {
      if (ISA_aa (bsp->mol) || bsp->repr == Seq_repr_seg || SeqMgrGetParentOfPart (bsp, NULL) != NULL) {
        /* skip proteins, segmented bioseqs, or segmented parts */
      } else {
        sfp = SeqMgrGetNextFeature (bvsp->bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
        if (sfp == NULL) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_UnnecessaryBioSourceFocus, "BioSource descriptor has focus, but no BioSource feature");
        }
      }
    }
    if (biop != NULL && biop->origin == 5) {
      bsp = bvsp->bsp;
      if (! IsOtherDNA (bsp)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol other should be used if Biosource-location is synthetic");
      }
    }
    /* ValidateBioSource (vsp, gcp, biop, NULL, vnp); */
    this_org = biop->org;
    /* fall into Seq_descr_org */
  case Seq_descr_org:
    if (this_org == NULL)
      this_org = (OrgRefPtr) (vnp->data.ptrvalue);
    if (bvsp->last_org != NULL) {
      if ((this_org->taxname != NULL) && (bvsp->last_org->taxname != NULL)) {
        if (StringCmp (this_org->taxname, bvsp->last_org->taxname)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent taxnames [%s] and [%s]", this_org->taxname, bvsp->last_org->taxname);
        }
      }
    } else
      bvsp->last_org = this_org;

    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_source || vnp2->choice == Seq_descr_org) {
        that_org = NULL;
        if (vnp2->choice == Seq_descr_source) {
          that_org = ((BioSourcePtr) (vnp2->data.ptrvalue))->org;
        }
        if (that_org == NULL) {
          that_org = (OrgRefPtr) (vnp2->data.ptrvalue);
        }
        if (that_org != NULL) {
          if ((this_org->taxname != NULL) && (that_org->taxname != NULL) && StringCmp (this_org->taxname, that_org->taxname) == 0) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MultipleBioSources, "Undesired multiple source descriptors");
          }
        }
      }
    }
    break;
  case Seq_descr_title:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Title descriptor needs text");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_title) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MultipleTitles, "Undesired multiple title descriptors");
      }
    }
    len = StringLen (str);
    if (len > 4) {
      ch = str [len - 1];
      while (ch == ' ' && len > 4) {
        len--;
        ch = str [len - 1];
      }
      if (ch == '.' && len > 4) {
        len--;
        ch = str [len - 1];
      }
      if (ch == '.' || ch == ',' || ch == ';' || ch == ':') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadPunctuation, "Title descriptor ends in bad punctuation");
      }
    }
    break;
  case Seq_descr_name:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Name descriptor needs text");
    }
    for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
      if (vnp2->choice == Seq_descr_name) {
        ptr = (CharPtr) vnp2->data.ptrvalue;
        if (StringDoesHaveText (ptr) && StringICmp (str, ptr) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleNames, "Undesired multiple name descriptors, identical text");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_MultipleNames, "Undesired multiple name descriptors, different text");
        }
      }
    }
  case Seq_descr_region:
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_MissingText, "Region descriptor needs text");
    }
    break;
  case Seq_descr_user:
    uop = (UserObjectPtr) vnp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "StructuredComment") == 0) {
          if (uop->data == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_UserObjectProblem, "Structured Comment user object descriptor is empty");
          }
        }
      }
    }
    break;
  case Seq_descr_pub:
    bvsp->got_a_pub = TRUE;
    pdp = (PubdescPtr) vnp->data.ptrvalue;
    /*
       ValidatePubdesc (vsp, pdp);
     */
    break;
  case Seq_descr_molinfo:
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->biomol) {
      case MOLECULE_TYPE_PEPTIDE:      /* peptide */
        if (!bvsp->is_aa) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with Molinfo-biomol = peptide");          
        }
        break;
      case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL:
        if (! bvsp->is_artificial) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol = other genetic");
        }
        break;
      case 0:                  /* unknown */
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol unknown used");
        break;
      case 255:                /* other */
        if (! IsXr (vnp)) {
          bsp = bvsp->bsp;
          if (! IsSynthetic (bsp)) {
            if (! IsMicroRNA (bsp)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol other used");
            }
          }
        }
        break;
      default:                 /* the rest are nucleic acid */
        if (bvsp->is_aa) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol [%d] used on protein", (int) mip->biomol);
        } else {
          if (bvsp->last_biomol) {
            if (bvsp->last_biomol != (int) mip->biomol) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-biomol [%d] and [%d]", bvsp->last_biomol, (int) mip->biomol);
            }
          } else {
            bvsp->last_biomol = (int) mip->biomol;
          }
        }
        break;
      }

      if (bvsp->is_syn_constr) {
        if (mip->biomol != MOLECULE_TYPE_OTHER_GENETIC_MATERIAL && !bvsp->is_aa) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "synthetic construct should have other-genetic");
        }
        if (! bvsp->is_artificial) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "synthetic construct should have artificial origin");
        }
      } else if (bvsp->is_artificial) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "artificial origin should have other-genetic and synthetic construct");
      }
      if (bvsp->is_artificial) {
        if (mip->biomol != MOLECULE_TYPE_OTHER_GENETIC_MATERIAL && mip->biomol != MOLECULE_TYPE_PEPTIDE) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_InvalidForType, "artificial origin should have other-genetic");
        }
      }
      if (!bvsp->is_aa) {
        switch (mip->tech) {
        case MI_TECH_concept_trans:
        case MI_TECH_seq_pept:
        case MI_TECH_both:
        case MI_TECH_seq_pept_overlap:
        case MI_TECH_seq_pept_homol:
        case MI_TECH_concept_trans_a:
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
          break;
        default:
          break;
        }
      } else {
        switch (mip->tech) {
        case MI_TECH_est:
        case MI_TECH_sts:
        case MI_TECH_genemap:
        case MI_TECH_physmap:
        case MI_TECH_htgs_1:
        case MI_TECH_htgs_2:
        case MI_TECH_htgs_3:
        case MI_TECH_fli_cdna:
        case MI_TECH_htgs_0:
        case MI_TECH_htc:
        case MI_TECH_wgs:
        case MI_TECH_barcode:
        case MI_TECH_composite_wgs_htgs:
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Protein with nucleic acid sequence method");
          break;
        default:
          break;
        }
      }
      if (bvsp->last_tech) {
        if (bvsp->last_tech != (int) mip->tech) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-tech [%d] and [%d]", bvsp->last_tech, (int) mip->tech);
        }
      } else {
        bvsp->last_tech = (int) mip->tech;
      }
      if (bvsp->last_completeness) {
        if (bvsp->last_completeness != (int) mip->completeness) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-completeness [%d] and [%d]",
                    bvsp->last_completeness, (int) mip->completeness);
        }
      } else {
        bvsp->last_completeness = (int) mip->completeness;
      }
    }
    break;
  default:
    break;
  }

  if (keywords != NULL) {
    tpa_exp = FALSE;
    tpa_inf = FALSE;
    for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
      if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:experimental") == 0) {
        tpa_exp = TRUE;
      } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:inferential") == 0) {
        tpa_inf = TRUE;
      }
    }
    if (tpa_exp && tpa_inf) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "TPA:experimental and TPA:inferential should not both be in the same set of keywords");
    }
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  return TRUE;
}

static Boolean LIBCALLBACK ValidateSeqDescrIndexed (ValNodePtr sdp, SeqMgrDescContextPtr context)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;

  bvsp = (BioseqValidStrPtr) context->userdata;
  vsp = bvsp->vsp;

  return ValidateSeqDescrCommon (sdp, bvsp, vsp, context->itemID);
}

static void ValidateSeqDescrContext (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  ValNodePtr      sdp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  sdp = (ValNodePtr) (gcp->thisitem);

  ValidateSeqDescrCommon (sdp, bvsp, vsp, 0);
}

/*****************************************************************************
*
*   ValidateBioseqContextGather(gcp)
*      Gather callback for validating context on a Bioseq
*
*****************************************************************************/
static Boolean DifferentDbxrefs (ValNodePtr dbxref1, ValNodePtr dbxref2)
{
  DbtagPtr        dbt1, dbt2;
  ObjectIdPtr     oip1, oip2;

  if (dbxref1 == NULL || dbxref2 == NULL)
    return FALSE;
  dbt1 = (DbtagPtr) dbxref1->data.ptrvalue;
  dbt2 = (DbtagPtr) dbxref2->data.ptrvalue;
  if (dbt1 == NULL || dbt2 == NULL)
    return FALSE;
  if (StringICmp (dbt1->db, dbt2->db) != 0)
    return TRUE;
  oip1 = dbt1->tag;
  oip2 = dbt2->tag;
  if (oip1 == NULL || oip2 == NULL)
    return FALSE;
  if (oip1->str == NULL && oip2->str == NULL) {
    if (oip1->id != oip2->id)
      return TRUE;
  } else {
    if (StringICmp (oip1->str, oip2->str) != 0)
      return TRUE;
  }
  return FALSE;
}

static Boolean FlybaseDbxrefs (ValNodePtr vnp)

{
  DbtagPtr  dbt;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "FLYBASE") == 0 || StringCmp (dbt->db, "FlyBase") == 0) {
        return TRUE;
      }
    }
    vnp = vnp->next;
  }
  return FALSE;
}

static Boolean GPSorNTorNCorNGorNW (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return TRUE;
    }
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean IsGenBankAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GENBANK) return TRUE;
    }
  }
  return FALSE;
}

static Boolean IsEMBLAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_EMBL) return TRUE;
    }
  }
  return FALSE;
}

static Boolean IsGeneralAccn (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr  bsp;
  DbtagPtr   dbt;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag(dbt)) continue;
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean NGorNT (SeqEntryPtr sep, SeqLocPtr location, BoolPtr is_nc)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (is_nc != NULL) {
    *is_nc = FALSE;
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0 && is_nc != NULL) {
            *is_nc = TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean GPSorRefSeq (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqIdPtr      sip;

  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      return TRUE;
    }
  }
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean IsNCorNT (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean IsNCorNTorNW (SeqEntryPtr sep, SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  bsp = BioseqFindFromSeqLoc (location);
  if (bsp != NULL) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            return TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean NotPeptideException (SeqFeatPtr sfp, SeqFeatPtr last)
{
  if (sfp != NULL && sfp->excpt) {
    if (StringISearch (sfp->except_text, "alternative processing") != NULL)
      return FALSE;
  }
  if (last != NULL && last->excpt) {
    if (StringISearch (last->except_text, "alternative processing") != NULL)
      return FALSE;
  }
  return TRUE;
}

static Boolean DescsSame (AnnotDescrPtr adp1, AnnotDescrPtr adp2)

{
  if (adp1 == NULL || adp2 == NULL) return TRUE;
  if (adp1->choice != adp2->choice) return FALSE;
  if (adp1->choice == Annot_descr_name || adp1->choice == Annot_descr_title) {
    if (StringICmp ((CharPtr) adp1->data.ptrvalue, (CharPtr) adp2->data.ptrvalue) == 0) return TRUE;
  }
  return FALSE;
}

typedef struct gmcdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  feat;
} GmcData, PNTR GmcDataPtr;

static int LIBCALLBACK SortGmcByGenePtr (
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

static void ValidateLocusTagGeneral (ValidStructPtr vsp, BioseqPtr bsp)

{
  DbtagPtr           dbt;
  SeqMgrFeatContext  fcontext;
  GatherContextPtr   gcp;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  ObjectIdPtr        oip;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  BioseqPtr          prod;
  CharPtr            ptr;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  Char               tmp [64];

  if (vsp == NULL || bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.subtype == FEATDEF_CDS || sfp->idx.subtype == FEATDEF_mRNA) {
      grp = SeqMgrGetGeneXref (sfp);
      if (! SeqMgrGeneIsSuppressed (grp)) {
        if (grp == NULL) {
          gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          if (gene != NULL) {
            grp = (GeneRefPtr) sfp->data.value.ptrvalue;
          }
        }
        if (grp != NULL && StringDoesHaveText (grp->locus_tag)) {
          prod = BioseqFindFromSeqLoc (sfp->product);
          if (prod != NULL) {
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice != SEQID_GENERAL) continue;
              dbt = (DbtagPtr) sip->data.ptrvalue;
              if (dbt == NULL) continue;
              if (IsSkippableDbtag(dbt)) continue;
              oip = dbt->tag;
              if (oip == NULL) continue;
              if (StringHasNoText (oip->str)) continue;
              StringNCpy_0 (tmp, oip->str, sizeof (tmp));
              ptr = StringChr (tmp, '-');
              if (ptr != NULL) {
                *ptr = '\0';
              }
              if (StringICmp (grp->locus_tag, tmp) != 0) {
                if (gcp != NULL) {
                  gcp->itemID = sfp->idx.itemID;
                  gcp->thistype = OBJ_SEQFEAT;
                }
                vsp->descr = NULL;
                vsp->sfp = sfp;
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProductMismatch, "Gene locus_tag does not match general ID of product");
              }
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean ReplaceQualsDiffer (GBQualPtr sfpqual, GBQualPtr lastqual)

{
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  while (sfpqual != NULL && StringICmp (sfpqual->qual, "replace") != 0) {
    sfpqual = sfpqual->next;
  }
  while (lastqual != NULL && StringICmp (lastqual->qual, "replace") != 0) {
    lastqual = lastqual->next;
  }
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  if (StringICmp (sfpqual->val, lastqual->val) != 0) return TRUE;

  return FALSE;
}

static Boolean GBQualsDiffer (GBQualPtr sfpqual, GBQualPtr lastqual)

{
  if (sfpqual == NULL || lastqual == NULL) return FALSE;

  /* depends upon sorted order of gbquals imposed by BasicSeqEntryCleanup */

  while (sfpqual != NULL && lastqual != NULL) {
    if (StringICmp (sfpqual->qual, lastqual->qual) != 0) return TRUE;
    if (StringICmp (sfpqual->val, lastqual->val) != 0) return TRUE;
    sfpqual = sfpqual->next;
    lastqual = lastqual->next;
  }

  if (sfpqual != NULL || lastqual != NULL) return TRUE;

  return FALSE;
}

static CharPtr MakePubLabelString (PubdescPtr pdp)

{
  Char        buf [521];
  CitGenPtr   cgp;
  ValNodePtr  vnp;

  if (pdp == NULL) return NULL;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    return StringSaveNoNull (buf);
  }

  return NULL;
}

static CharPtr ValGetAuthorsPlusConsortium (
  AuthListPtr alp
)

{
  CharPtr  consortium;
  CharPtr  str;
  CharPtr  tmp;

  consortium = NULL;
  str = GetAuthorsString (GENBANK_FMT, alp, &consortium, NULL, NULL);
  if (str == NULL) return consortium;
  if (consortium == NULL) return str;
  tmp = MemNew (StringLen (str) + StringLen (consortium) + 5);
  if (tmp == NULL) return NULL;
  StringCpy (tmp, str);
  StringCat (tmp, "; ");
  StringCat (tmp, consortium);
  MemFree (str);
  MemFree (consortium);
  return tmp;
}

static Boolean IsIdenticalPublication (PubdescPtr pdp1, PubdescPtr pdp2)

{
  AuthListPtr  alp1, alp2;
  Boolean      rsult = TRUE;
  CharPtr      str1, str2;

  if (pdp1 == NULL || pdp2 == NULL) return FALSE;

  str1 = MakePubLabelString (pdp1);
  str2 = MakePubLabelString (pdp2);
  if (StringDoesHaveText (str1) && StringDoesHaveText (str2)) {
    if (StringICmp (str1, str2) != 0) {
       rsult = FALSE;
    }
  }
  MemFree (str1);
  MemFree (str2);
  if (! rsult) return rsult;

  alp1 = GetAuthListPtr (pdp1, NULL);
  alp2 = GetAuthListPtr (pdp2, NULL);
  if (alp1 != NULL && alp2 != NULL) {
    str1 = ValGetAuthorsPlusConsortium (alp1);
    str2 = ValGetAuthorsPlusConsortium (alp2);
    if (StringDoesHaveText (str1) && StringDoesHaveText (str2)) {
      if (StringICmp (str1, str2) != 0) {
         rsult = FALSE;
      }
    }
    MemFree (str1);
    MemFree (str2);
  }

  return rsult;
}

static Boolean IsIdenticalBioSource (BioSourcePtr biop1, BioSourcePtr biop2)

{
  DbtagPtr      dbt1, dbt2;
  ObjectIdPtr   oip1, oip2;
  OrgModPtr     omp1, omp2;
  OrgNamePtr    onp1, onp2;
  OrgRefPtr     orp1, orp2;
  SubSourcePtr  ssp1, ssp2;
  ValNodePtr    vnp1, vnp2;

  if (biop1 == NULL || biop2 == NULL) return FALSE;

  if (biop1->is_focus != biop2->is_focus) return FALSE;

  orp1 = biop1->org;
  orp2 = biop2->org;
  if (orp1 == NULL || orp2 == NULL) return FALSE;
  if (StringICmp (orp1->taxname, orp2->taxname) != 0) return FALSE;

  onp1 = orp1->orgname;
  onp2 = orp2->orgname;
  if (onp1 == NULL || onp2 == NULL) return FALSE;

  omp1 = onp1->mod;
  omp2 = onp2->mod;
  while (omp1 != NULL && omp2 != NULL) {
    if (omp1->subtype != omp2->subtype) return FALSE;
    if (StringICmp (omp1->subname, omp2->subname) != 0) return FALSE;
    omp1 = omp1->next;
    omp2 = omp2->next;
  }
  if (omp1 != NULL || omp2 != NULL) return FALSE;
  
  ssp1 = biop1->subtype;
  ssp2 = biop2->subtype;
  while (ssp1 != NULL && ssp2 != NULL) {
    if (ssp1->subtype != ssp2->subtype) return FALSE;
    if (StringICmp(ssp1->name, ssp2->name) != 0) return FALSE;
    ssp1 = ssp1->next;
    ssp2 = ssp2->next;
  }
  if (ssp1 != NULL || ssp2 != NULL) return FALSE;

  vnp1 = orp1->db;
  vnp2 = orp2->db;
  while (vnp1 != NULL && vnp2 != NULL) {
    dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
    dbt2 = (DbtagPtr) vnp2->data.ptrvalue;

    if ((dbt1 != NULL) && (dbt2 != NULL)) {
      if (StringCmp (dbt1->db, dbt2->db) != 0) return FALSE;

      oip1 = dbt1->tag;
      oip2 = dbt2->tag;
      if ((oip1 != NULL) && (oip2 != NULL)) {
        if (oip1->str != NULL) {
          if (StringICmp(oip1->str, oip2->str) != 0) return FALSE;
        } else  {
          if (oip1->id != oip2->id) return FALSE;
        }
      }
      else if (oip1 != NULL)
        return FALSE;
      else if (oip2 != NULL)
        return FALSE;
    }
    else if (dbt1 != NULL)
      return FALSE;
    else if (dbt2 != NULL)
      return FALSE;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 != NULL || vnp2 != NULL) return FALSE;

  return TRUE;
}

typedef struct lpdata {
  Int2        count;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
  Char        firstid [64];
  Boolean     products_unique;
  Boolean     featid_matched;
} LpData, PNTR LpDataPtr;

static Boolean IdXrefsAreReciprocal (
  SeqFeatPtr cds,
  SeqFeatPtr mrna
)

{
  SeqFeatXrefPtr  xref;
  Boolean         match1 = FALSE, match2 = FALSE;
  SeqFeatPtr      matchsfp;

  if (cds == NULL || mrna == NULL) return FALSE;
  if (cds->id.choice != 3 || mrna->id.choice != 3) return FALSE;

  for (xref = cds->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == mrna) {
        match1 = TRUE;
      }
    }
  }

  for (xref = mrna->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (mrna->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == cds) {
        match2 = TRUE;
      }
    }
  }

  if (match1 && match2) return TRUE;
  return FALSE;
}

static Int2 IdXrefsNotReciprocal (
  SeqFeatPtr cds,
  SeqFeatPtr mrna
)

{
  Int4            giu = 0, gip = 0;
  SeqFeatPtr      matchsfp;
  ObjectIdPtr     oip;
  SeqIdPtr        sip;
  CharPtr         tmp;
  UserFieldPtr    ufp;
  UserObjectPtr   uop;
  SeqFeatXrefPtr  xref;

  if (cds == NULL || mrna == NULL) return 0;
  if (cds->id.choice != 3 || mrna->id.choice != 3) return 0;

  for (xref = cds->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp != mrna) return 1;
    }
  }

  for (xref = mrna->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (mrna->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp != cds) return 1;
    }
  }

  if (cds->product == NULL) return 0;
  if (mrna->ext == NULL) return 0;
  uop = FindUopByTag (mrna->ext, "MrnaProteinLink");
  if (uop == NULL) return 0;
  sip = SeqLocId (cds->product);
  if (sip == NULL) return 0;
  if (sip->choice == SEQID_GI) {
    gip = (Int4) sip->data.intvalue;
  } else {
    gip = GetGIForSeqId (sip);
  }
  if (gip == 0) return 0;
  ufp = uop->data;
  if (ufp == NULL || ufp->choice != 1) return 0;
  oip = ufp->label;
  if (oip == NULL || StringICmp (oip->str, "protein seqID") != 0) return 0;
  tmp = (CharPtr) ufp->data.ptrvalue;
  if (StringHasNoText (tmp)) return 0;
  sip = MakeSeqID (tmp);
  if (sip == NULL) return 0;
  if (sip->choice == SEQID_GI) {
    giu = (Int4) sip->data.intvalue;
  } else {
    giu = GetGIForSeqId (sip);
  }
  SeqIdFree (sip);
  if (giu == 0) return 0;
  if (gip != giu) return 2;

  return 0;
}

static Boolean LIBCALLBACK FindSingleMrnaProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  Char        buf [64];
  SeqFeatPtr  cds;
  LpDataPtr   ldp;
  SeqIdPtr    sip;
  VvmDataPtr  vdp;

  ldp = (LpDataPtr) context->userdata;
  if (ldp == NULL) return TRUE;
  cds = ldp->cds;
  if (cds == NULL) return TRUE;

  if (sfp->product) {
    if (StringHasNoText (ldp->firstid)) {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, ldp->firstid, PRINTID_FASTA_LONG, sizeof (ldp->firstid) - 1);
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
      if (StringCmp (ldp->firstid, buf) == 0) {
        ldp->products_unique = FALSE;
      }
    }
  }

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->accounted_for) return TRUE;

  (ldp->count)++;
  ldp->mrna = sfp;

  if (IdXrefsAreReciprocal (cds, sfp)) {
    ldp->featid_matched = TRUE;
  }

  return TRUE;
}

/*
static Boolean LIBCALLBACK DummyCM121Proc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  return TRUE;
}
*/

static void ValidateCDSmRNAmatch (
  ValidStructPtr vsp,
  BioseqPtr bsp,
  Int2 numgene,
  Int2 numcds,
  Int2 nummrna,
  Boolean suppress_duplicate_messages
)

{
  BioSourcePtr       biop;
  ValNodePtr         cdshead = NULL;
  ValNodePtr         cdstail = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext, rcontext;
  GatherContextPtr   gcp;
  GmcDataPtr         gdp, head;
  SeqFeatPtr         gene;
  Boolean            goOn, pseudo;
  GeneRefPtr         grp;
  Int2               i, j, k, numfeats, tmpnumcds, tmpnummrna, count;
  Boolean            is_genbank = FALSE;
  LpData             ld;
  Int2               num_no_mrna = 0;
  Int4               num_repeat_regions;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  Int2               recip;
  VoidPtr            repeat_region_array;
  SeqFeatPtr         rpt_region;
  SeqDescrPtr        sdp;
  ErrSev             sev = /* SEV_INFO */ SEV_WARNING;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  VvmDataPtr         vdp;
  ValNodePtr         vnp;

  if (vsp == NULL || bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  gcp = vsp->gcp;
  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  /*
  if (GetAppProperty ("ValidateCDSmRNAoneToOne") != NULL) {
    cdsMrnaOneToOne = TRUE;
  }
  */

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      sev = SEV_WARNING;
    } else if (sip->choice == SEQID_GENBANK) {
      is_genbank = TRUE;
    }
  }

  if (is_genbank) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (sdp != NULL) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          onp = orp->orgname;
          if (onp != NULL) {
            if (StringDoesHaveText (onp->div) &&
                StringCmp (onp->div, "BCT") != 0 &&
                StringCmp (onp->div, "VRL") != 0) {
              is_genbank = FALSE;
            }
          }
        }
      }
    }
  }

  repeat_region_array = SeqMgrBuildFeatureIndex (bsp, &num_repeat_regions, 0, FEATDEF_repeat_region);

  if (numgene > 0 && numcds > 0 && nummrna > 0) {
    numfeats = numcds + nummrna;
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numfeats + 1));
    if (head != NULL) {
      gdp = head;
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (sfp->idx.subtype == FEATDEF_CDS || sfp->idx.subtype == FEATDEF_mRNA) {
          gdp->feat = sfp;
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
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }
      HeapSort (head, (size_t) numfeats, sizeof (GmcData), SortGmcByGenePtr);
      for (i = 0; i < numfeats; i += j) {
        gene = head [i].gene;
        for (j = 1; i + j < numfeats && gene == head [i + j].gene; j++) continue;
        if (j > 1 && gene != NULL) {
          /* is alt splicing */
          tmpnumcds = 0;
          tmpnummrna = 0;
          for (k = 0; k < j; k++) {
            sfp = head [i + k].feat;
            if (sfp == NULL) continue;
            if (sfp->idx.subtype == FEATDEF_CDS) {
              tmpnumcds++;
            }
            if (sfp->idx.subtype == FEATDEF_mRNA) {
              tmpnummrna++;
            }
          }
          if (tmpnumcds > 0 && tmpnummrna > 1 && tmpnumcds != tmpnummrna && (! is_genbank)) {

            if (gcp != NULL) {
              gcp->itemID = gene->idx.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = gene;
            ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNAmismatch, "mRNA count (%d) does not match CDS (%d) count for gene",
                      (int) tmpnummrna, (int) tmpnumcds);
          }
        }
      }
    }
    MemFree (head);
  }

  /* loop through CDS features, finding single unused mRNA partner */

  goOn = TRUE;
  while (goOn) {
    goOn = FALSE;
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      vdp = (VvmDataPtr) sfp->idx.scratch;
      if (vdp != NULL && (! vdp->accounted_for)) {
        vdp->num_mrnas = 0;
        ld.count = 0;
        ld.cds = sfp;
        ld.mrna = NULL;
        ld.firstid [0] = '\0';
        ld.products_unique = TRUE;
        ld.featid_matched = FALSE;
        if (sfp->excpt &&
            (StringISearch (sfp->except_text, "ribosomal slippage") != NULL ||
             StringISearch (sfp->except_text, "trans-splicing") != NULL)) {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   LOCATION_SUBSET, (Pointer) &ld, FindSingleMrnaProc);
        } else {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   CHECK_INTERVALS, (Pointer) &ld, FindSingleMrnaProc);
        }
        if (ld.count == 1 && ld.mrna != NULL) {
          vdp->accounted_for = TRUE;
          vdp->num_mrnas = ld.count;
          vdp->featid_matched = ld.featid_matched;
          vdp = (VvmDataPtr) ld.mrna->idx.scratch;
          if (vdp != NULL) {
            vdp->accounted_for = TRUE;
            goOn = TRUE;
            recip = IdXrefsNotReciprocal (sfp, ld.mrna);
            if (recip == 1) {
              if (gcp != NULL) {
                gcp->itemID = sfp->idx.itemID;
                gcp->thistype = OBJ_SEQFEAT;
              }
              vsp->descr = NULL;
              vsp->sfp = sfp;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefNotReciprocal, "CDS/mRNA unambiguous pair have erroneous cross-references");
            } else if (recip == 2) {
              if (gcp != NULL) {
                gcp->itemID = ld.mrna->idx.itemID;
                gcp->thistype = OBJ_SEQFEAT;
              }
              vsp->descr = NULL;
              vsp->sfp = ld.mrna;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "MrnaProteinLink inconsistent with feature ID cross-references");
            }
          }
        } else {
          vdp->num_mrnas = ld.count;
          vdp->products_unique = ld.products_unique;
          vdp->featid_matched = ld.featid_matched;
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL && (! is_genbank)) {
    vdp = (VvmDataPtr) sfp->idx.scratch;
    if (vdp != NULL) {
      count = vdp->num_mrnas;
      /*
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                 CHECK_INTERVALS, NULL, DummyCM121Proc);
      */
      if (count > 1) {
        if (gcp != NULL) {
          gcp->itemID = sfp->idx.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (vdp->featid_matched) {
          /* presence of reciprocal link suppresses warnings */
        } else if (vdp->products_unique) {
          /*
          if (! suppress_duplicate_messages) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_CDSwithMultipleMRNAs,
                      "CDS overlapped by %d mRNAs, but product locations are unique", (int) count);
          }
          */
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_CDSwithMultipleMRNAs,
                    "CDS overlapped by %d mRNAs, but product locations are unique", (int) count);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithMultipleMRNAs, "CDS overlapped by %d mRNAs", (int) count);
        }
      } else if (count == 0 && numgene > 0 && numcds > 0 && nummrna > 0) {
        pseudo = sfp->pseudo;
        if (! pseudo) {
          grp = SeqMgrGetGeneXref (sfp);
          if (grp != NULL) {
            pseudo = grp->pseudo;
          } else {
            gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
            if (gene != NULL) {
              pseudo = gene->pseudo;
              if (! pseudo) {
                grp = (GeneRefPtr) gene->data.value.ptrvalue;
                if (grp != NULL) {
                  pseudo = grp->pseudo;
                }
              }
            }
          }
        }
        if (! pseudo) {
          rpt_region = SeqMgrGetOverlappingFeature (sfp->location, 0, repeat_region_array, num_repeat_regions,
                                                    NULL, CONTAINED_WITHIN, &rcontext);
          if (rpt_region == NULL) {
            /*
            if (gcp != NULL) {
              gcp->itemID = sfp->idx.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = sfp;
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap, "CDS overlapped by 0 mRNAs");
            */
            vnp = ValNodeAddPointer (&cdstail, 0, (Pointer) sfp);
            if (cdshead == NULL) {
              cdshead = vnp;
            }
            cdstail = vnp;
            num_no_mrna++;
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  MemFree (repeat_region_array);

  if (num_no_mrna > 0) {
    if (num_no_mrna >= 10) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      vsp->descr = NULL;
      vsp->sfp = NULL;
      vsp->bsp = bsp;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap,
                "%d out of %d CDSs overlapped by 0 mRNAs", (int) num_no_mrna, (int) numcds);
    } else {
      for (vnp = cdshead; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp == NULL) continue;
        if (gcp != NULL) {
          gcp->itemID = sfp->idx.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDSwithNoMRNAOverlap, "CDS overlapped by 0 mRNAs");
      }
    }
  }

  ValNodeFree (cdshead);

  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
}

static Boolean HaveUniqueFeatIDXrefs (SeqFeatXrefPtr xref1, SeqFeatXrefPtr xref2)

{
  ObjectIdPtr  oip1 = NULL, oip2 = NULL;

  while (xref1 != NULL) {
    if (xref1->id.choice == 3) {
      oip1 = (ObjectIdPtr) xref1->id.value.ptrvalue;
    }
    xref1 = xref1->next;
  }

  while (xref2 != NULL) {
    if (xref2->id.choice == 3) {
      oip2 = (ObjectIdPtr) xref2->id.value.ptrvalue;
    }
    xref2 = xref2->next;
  }

  if (oip1 == NULL || oip2 == NULL) return FALSE;
  if (oip1->str == NULL && oip2->str == NULL) {
    if (oip1->id != oip2->id && oip1->id > 0 && oip2->id > 0) return TRUE;
  }

  return FALSE;
}

#define LEFT_RIBOSOMAL_SUBUNIT  1
#define INTERNAL_SPACER_1        2
#define MIDDLE_RIBOSOMAL_SUBUNIT 3
#define INTERNAL_SPACER_2        4
#define RIGHT_RIBOSOMAL_SUBUNIT  5
#define INTERNAL_SPACER_X        6
#define TRANSFER_RNA             7

static Int2 WhichRNA (SeqFeatPtr sfp)

{
  GBQualPtr  gbq;
  RnaRefPtr  rrp;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return 0;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return 0;
  if (rrp->type == 3) {
    return TRANSFER_RNA;
  }
  if (rrp->ext.choice != 1) return 0;
  str = (CharPtr) rrp->ext.value.ptrvalue;
  if (StringHasNoText (str)) return 0;
  if (rrp->type == 4) {
    if (StringNICmp (str, "small ", 6) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "18S ", 4) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "16S ", 4) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "5.8S ", 5) == 0) return MIDDLE_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "large ", 6) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "26S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "28S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "23S ", 4) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    /* variant spellings */
    if (StringNICmp (str, "18 ", 3) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "16 ", 3) == 0) return LEFT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "5.8 ", 4) == 0) return MIDDLE_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "26 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "28 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
    if (StringNICmp (str, "23 ", 3) == 0) return RIGHT_RIBOSOMAL_SUBUNIT;
  }
  if (rrp->type == 255) {
    if (StringICmp (str, "misc_RNA") == 0) {
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringICmp (gbq->qual, "product") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        str = gbq->val;
      }
    }
    if (StringICmp (str, "internal transcribed spacer 1") == 0) return INTERNAL_SPACER_1;
    if (StringICmp (str, "internal transcribed spacer 2") == 0) return INTERNAL_SPACER_2;
    /* variant spellings */
    if (StringICmp (str, "internal transcribed spacer1") == 0) return INTERNAL_SPACER_1;
    if (StringICmp (str, "internal transcribed spacer2") == 0) return INTERNAL_SPACER_2;
    if (StringICmp (str, "internal transcribed spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "ITS") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "16S-23S ribosomal RNA intergenic spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "16S-23S intergenic spacer") == 0) return INTERNAL_SPACER_X;
    if (StringICmp (str, "intergenic spacer") == 0) return INTERNAL_SPACER_X;
  }
  return 0;
}

static Boolean CDSsLinkedToDifferentMRNAs (SeqFeatPtr sfp, SeqFeatPtr last)

{
  SeqFeatPtr      mrna1 = NULL, mrna2 = NULL;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || last == NULL) return FALSE;
  if (sfp->idx.subtype != FEATDEF_CDS || last->idx.subtype != FEATDEF_CDS) return FALSE;

  for (xref = sfp->xref; xref != NULL && mrna1 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      mrna1 = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (mrna1 != NULL && mrna1->idx.subtype != FEATDEF_mRNA) {
        mrna1 = NULL;
      }
    }
  }

  for (xref = last->xref; xref != NULL && mrna2 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      mrna2 = SeqMgrGetFeatureByFeatID (last->idx.entityID, NULL, NULL, xref, NULL);
      if (mrna2 != NULL && mrna2->idx.subtype != FEATDEF_mRNA) {
        mrna2 = NULL;
      }
    }
  }

  if (mrna1 != NULL && mrna2 != NULL && mrna1 != mrna2) return TRUE;

  return FALSE;
}

static Boolean MRNAsLinkedToDifferentCDSs (SeqFeatPtr sfp, SeqFeatPtr last)

{
  SeqFeatPtr      cds1 = NULL, cds2 = NULL;
  CdRegionPtr     crp1, crp2;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || last == NULL) return FALSE;
  if (sfp->idx.subtype != FEATDEF_mRNA || last->idx.subtype != FEATDEF_mRNA) return FALSE;

  for (xref = sfp->xref; xref != NULL && cds1 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      cds1 = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (cds1 != NULL && cds1->idx.subtype != FEATDEF_CDS) {
        cds1 = NULL;
      }
    }
  }

  for (xref = last->xref; xref != NULL && cds2 == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      cds2 = SeqMgrGetFeatureByFeatID (last->idx.entityID, NULL, NULL, xref, NULL);
      if (cds2 != NULL && cds2->idx.subtype != FEATDEF_CDS) {
        cds2 = NULL;
      }
    }
  }

  if (cds1 == NULL || cds2 == NULL || cds1 == cds2) return FALSE;

  crp1 = (CdRegionPtr) cds1->data.value.ptrvalue;
  crp2 = (CdRegionPtr) cds2->data.value.ptrvalue;
  if (crp1 == NULL || crp2 == NULL) return FALSE;

  if (SeqLocCompare (cds1->location, cds2->location) != SLC_A_EQ_B) return TRUE;

  if (crp1->frame < 2 && crp2->frame < 2) return FALSE;
  if (crp1->frame != crp2->frame) return TRUE;

  return FALSE;
}

static Boolean BaseRangeIsVirtual (BioseqPtr bsp, Int4 left, Int4 right)

{
  Uint1        res;
  StreamCache  sc;

  if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) return FALSE;

  StreamCacheSetPosition (&sc, left - 1);
  res = StreamCacheGetResidue (&sc);
  if (res == '-') return FALSE;
  res = StreamCacheGetResidue (&sc);
  if (res != '-') return FALSE;

  StreamCacheSetPosition (&sc, right - 1);
  res = StreamCacheGetResidue (&sc);
  if (res != '-') return FALSE;
  res = StreamCacheGetResidue (&sc);
  if (res == '-') return FALSE;

  return TRUE;
}

static Boolean LIBCALLBACK GetFeatsInGaps (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  BioseqPtr         bsp;
  Int4              dashes;
  Int2              first = 0;
  GatherContextPtr  gcp;
  Int2              last = 0;
  Int4              len;
  SeqLocPtr         loc;
  Boolean           needToStream = TRUE;
  Int4              Ns;
  Uint2             olditemtype = 0;
  Uint4             olditemid = 0;
  Int4              plusses;
  Int2              prefix = 0;
  Int2              suffix = 0;
  Int4              realBases;
  Int2              res;
  StreamCache       sc;
  SeqIntPtr         sintp;
  ValidStructPtr    vsp;

  if (sfp == NULL || fcontext == NULL) return FALSE;
  vsp = (ValidStructPtr) fcontext->userdata;
  if (vsp == NULL) return FALSE;
  gcp = vsp->gcp;
  if (gcp == NULL) return FALSE;

  if (sfp->idx.subtype == FEATDEF_gap) return TRUE;
  loc = sfp->location;
  if (loc == NULL) return TRUE;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;

  gcp->itemID = fcontext->itemID;
  gcp->thistype = OBJ_SEQFEAT;
  vsp->sfp = sfp;


  dashes = 0;
  plusses = 0;
  Ns = 0;
  realBases = 0;

  /* special check for single interval misc_features that may exactly cover a gap */
  if (loc->choice == SEQLOC_INT && sfp->idx.subtype == FEATDEF_misc_feature) {
    sintp = (SeqIntPtr) loc->data.ptrvalue;
    if (sintp != NULL) {
      bsp = BioseqFind (sintp->id);
      if (bsp != NULL && sintp->from > 0 && sintp->to < bsp->length - 1) {
        len = SeqLocLen (loc);
        if (StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, &sc)) {
          StreamCacheSetPosition (&sc, sintp->from - 1);
          prefix = StreamCacheGetResidue (&sc);
          while ((res = StreamCacheGetResidue (&sc)) != '\0' && len > 0) {
            if (IS_LOWER (res)) {
              res = TO_UPPER (res);
            }
            if (first == 0) {
              first = res;
            }
            last = res;
            if (res == '-') {
              dashes++;
            } else if (res == '+') {
              plusses++;
            } else if (res == 'N') {
              Ns++;
            } else if (res != 0) {
              realBases++;
            }
            len--;
          }
          suffix = StreamCacheGetResidue (&sc);
          needToStream = FALSE;
        }
      }
    }
  }

  if (needToStream && StreamCacheSetup (NULL, loc, EXPAND_GAPS_TO_DASHES | KNOWN_GAP_AS_PLUS, &sc)) {
    while ((res = StreamCacheGetResidue (&sc)) != '\0') {
      if (IS_LOWER (res)) {
        res = TO_UPPER (res);
      }
      if (first == 0) {
        first = res;
      }
      last = res;
      if (res == '-') {
        dashes++;
      } else if (res == '+') {
        plusses++;
      } else if (res == 'N') {
        Ns++;
      } else if (res != 0) {
        realBases++;
      }
    }
  }

  if (dashes == 0 && plusses == 0 && Ns == 0) {
    /* ignore features that do not cover any gap characters */
  } else if (first == '-' || first == '+' || last == '-' || last == '+') {
    if (realBases > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureBeginsOrEndsInGap, "Feature begins or ends in gap");
    } else if (IS_ALPHA (prefix) && IS_ALPHA (suffix)) {
      /* ignore (misc_) features that exactly cover the gap */
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside sequence gap");
    }
  } else if (sfp->idx.subtype == FEATDEF_GENE) {
    /* ignore genes, unless they start or stop in a gap */
  } else if (dashes == 0 && plusses == 0 && Ns > 0) {
    if (realBases > 0) {
      /*
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCrossesGap, "Feature crosses gap of Ns");
      */
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside gap of Ns");
    }
  } else if (dashes > 0) {
    if (realBases > 0) {
      if (sfp->idx.subtype == FEATDEF_CDS || sfp->idx.subtype == FEATDEF_mRNA) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureCrossesGap, "Feature crosses gap of unknown length");
      }
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside gap of unknown length");
    }
  } else if (plusses > 0) {
    if (realBases > 0) {
      /* ignore feature completely containing gap of known length plus sequence on either side */
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureInsideGap, "Feature inside sequence gap");
    }
  } else {
    ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_FeatureInsideGap, "Problem with feature in gap calculation");
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;
  vsp->sfp = NULL;

  return TRUE;
}

static Boolean LIBCALLBACK MarkFeatsInGaps (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  VvmDataPtr         vdp;

  if (sfp == NULL) return TRUE;
  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp == NULL) return TRUE;

  vdp->feat_touches_gap = TRUE;

  return TRUE;
}

static void CheckBioseqForFeatsInGap (
  BioseqPtr bsp,
  ValidStructPtr vsp
)

{
  Int4               currpos = 0;
  SeqMgrFeatContext  fcontext;
  SeqLitPtr          litp;
  SeqFeatPtr         sfp;
  SeqInt             si;
  SeqIdPtr           sip;
  SeqLoc             sl;
  SeqLocPtr          slp;
  VvmDataPtr         vdp;
  ValNodePtr         vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || ISA_aa (bsp->mol)) return;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;

  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      if (slp == NULL) continue;
      currpos += SeqLocLen (slp);
    } else if (vnp->choice == 2) {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL) continue;
      if ((litp->seq_data == NULL || litp->seq_data_type == Seq_code_gap) &&
           litp->length > 0) {
        MemSet ((Pointer) &si, 0, sizeof (si));
        MemSet ((Pointer) &sl, 0, sizeof (sl));
        si.from = currpos;
        si.to = currpos + litp->length - 1;
        si.strand = Seq_strand_both;
        si.id = sip;
        sl.choice = SEQLOC_INT;
        sl.data.ptrvalue = (Pointer) &si;
        SeqMgrExploreFeatures (bsp, (Pointer) vsp, MarkFeatsInGaps, &sl, NULL, NULL);
        /* SeqMgrExploreFeatures (bsp, (Pointer) vsp, GetFeatsInGaps, &sl, NULL, NULL); */
      }
      currpos += litp->length;
    }
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    vdp = (VvmDataPtr) sfp->idx.scratch;
    if (vdp != NULL && vdp->feat_touches_gap) {
      fcontext.userdata = (Pointer) vsp;
      GetFeatsInGaps (sfp, &fcontext);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

static Boolean ReportGeneCollision (GeneRefPtr grp, GeneRefPtr lastgrp)

{
  if (grp == NULL || lastgrp == NULL) return TRUE;

  if (StringDoesHaveText (grp->locus) && StringDoesHaveText (lastgrp->locus)) {
    if (StringICmp (grp->locus, lastgrp->locus) == 0) return TRUE;
  }

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (lastgrp->locus_tag)) {
    if (StringICmp (grp->locus_tag, lastgrp->locus_tag) == 0) return TRUE;
  }

  if (StringDoesHaveText (grp->desc) && StringDoesHaveText (lastgrp->desc)) {
    if (StringICmp (grp->desc, lastgrp->desc) == 0) return FALSE;
  }

  return TRUE;
}

static Boolean ValidateBioseqContextIndexed (BioseqPtr bsp, BioseqValidStrPtr bvsp)

{
  ValidStructPtr     vsp;
  CitSubPtr          csp;
  ObjMgrDataPtr      omdp;
  SeqSubmitPtr       ssp;
  SubmitBlockPtr     sbp;
  GatherContextPtr   gcp;
  SeqFeatPtr         sfp, lastsfp;
  SeqMgrFeatContext  fcontext;
  Uint2              featdeftype = 0;
  Boolean            firstCDS;
  GeneRefPtr         grp, genomicgrp, lastgrp;
  SeqFeatPtr         last = NULL;
  Boolean            leave;
  CharPtr            label = NULL;
  CharPtr            comment = NULL;
  Int4               left = 0;
  Boolean            partialL = FALSE;
  Boolean            partialR = FALSE;
  Int4               right = 0;
  Uint1              strand = 0;
  Int2               numivals = 0;
  Int4Ptr            ivals = NULL;
  Boolean            ivalssame;
  SeqAnnotPtr        sap = NULL;
  Uint2              olditemtype = 0;
  Uint4              olditemid = 0;
  CharPtr            lastLabel;
  CharPtr            message;
  Int2               i;
  Boolean            isCuratedFlybase = FALSE;
  Boolean            isDrosophila = FALSE;
  Boolean            isGenBankAccn = FALSE;
  Boolean            isGeneralAccn = FALSE;
  Boolean            isGPSorNTorNCorNGorNW = FALSE;
  Boolean            isViral = FALSE;
  Int2               j;
  CdRegionPtr        crp;
  Uint1              frame = 0;
  Boolean            samelabel;
  int                severity;
  int                overlapPepSev;
  BioSourcePtr       biop = NULL, lastbiop;
  OrgRefPtr          orp = NULL;
  OrgNamePtr         onp = NULL;
  Int4               fiveUTRright;
  Int4               cdsRight;
  Int4               threeUTRright;
  Int4               cdscount, genecount, utr5count, utr3count;
  SeqFeatPtr         mrna, gene, cdsgene, utr5gene, utr3gene;
  PubdescPtr         pdp = NULL, lastpdp;
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  Boolean            showBadFullSource;
  Int2               numBadFullSource;
  SubSourcePtr       sbsp;
  Int2               numgene, numcds, nummrna, numcdsproducts, nummrnaproducts,
                     numcdspseudo, nummrnapseudo, lastrnatype, thisrnatype;
  Boolean            cds_products_unique = TRUE, mrna_products_unique = TRUE,
                     suppress_duplicate_messages = FALSE, pseudo;
  SeqIdPtr           sip;
  Char               buf [96];
  SeqFeatXrefPtr     xref = NULL;
  CharPtr            except_text = NULL;
  ValNodePtr         vnp, cds_prod_head = NULL, mrna_prod_head = NULL,
                     lastcdsprod = NULL, lastmrnaprod = NULL;
  StreamCache        sc;
  Int2               res;
  Int4               dashes;
  Int4               Ns;
  Int4               realBases;
  Int4               estimated_length;
  Int4               loclen;
  GBQualPtr          gbq;
  long int           val;
  SeqLocPtr          slp;
  MolInfoPtr         mip = NULL;
  SeqFeatPtr         cds;
  BioseqPtr          nbsp;

  gcp = bvsp->gcp;
  vsp = bvsp->vsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->gcp = gcp;               /* needed for ValidErr */

  vsp->rrna_array = SeqMgrBuildFeatureIndex (bsp, &(vsp->numrrna), 0, FEATDEF_rRNA);
  vsp->trna_array = SeqMgrBuildFeatureIndex (bsp, &(vsp->numtrna), 0, FEATDEF_tRNA);

  SeqMgrExploreFeatures (bsp, (Pointer) bvsp, ValidateSeqFeatIndexed, NULL, NULL, NULL);

  vsp->rrna_array = MemFree (vsp->rrna_array);
  vsp->trna_array = MemFree (vsp->trna_array);

  overlapPepSev = SEV_WARNING;
  if (GetAppProperty ("SpliceValidateAsError") != NULL) {
    overlapPepSev = SEV_ERROR;
  }

  if (gcp != NULL) {
    olditemid = gcp->itemID;
    olditemtype = gcp->thistype;
  }

  numgene = 0;
  numcds = 0;
  nummrna = 0;
  numcdsproducts = 0;
  nummrnaproducts = 0;
  numcdspseudo = 0;
  nummrnapseudo = 0;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_GENE :
        numgene++;
        break;
      case FEATDEF_CDS :
        numcds++;
        if (sfp->product != NULL) {
          numcdsproducts++;
          sip = SeqLocId (sfp->product);
          SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          if (StringDoesHaveText (buf)) {
            vnp = ValNodeCopyStr (&lastcdsprod, 0, buf);
            if (cds_prod_head == NULL) {
              cds_prod_head = vnp;
            }
            lastcdsprod = vnp;
          }
        } else {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            grp = SeqMgrGetGeneXref (sfp);
            if (grp != NULL) {
              pseudo = grp->pseudo;
            } else {
              gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
              if (gene != NULL) {
                pseudo = gene->pseudo;
                if (! pseudo) {
                  grp = (GeneRefPtr) gene->data.value.ptrvalue;
                  if (grp != NULL) {
                    pseudo = grp->pseudo;
                  }
                }
              }
            }
          }
          if (pseudo) {
            numcdspseudo++;
          }
        }
        break;
      case FEATDEF_mRNA :
        nummrna++;
        if (sfp->product != NULL) {
          nummrnaproducts++;
          sip = SeqLocId (sfp->product);
          SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          if (StringDoesHaveText (buf)) {
            vnp = ValNodeCopyStr (&lastmrnaprod, 0, buf);
            if (mrna_prod_head == NULL) {
              mrna_prod_head = vnp;
            }
            lastmrnaprod = vnp;
          }
        } else {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            grp = SeqMgrGetGeneXref (sfp);
            if (grp != NULL) {
              pseudo = grp->pseudo;
            } else {
              gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
              if (gene != NULL) {
                pseudo = gene->pseudo;
                if (! pseudo) {
                  grp = (GeneRefPtr) gene->data.value.ptrvalue;
                  if (grp != NULL) {
                    pseudo = grp->pseudo;
                  }
                }
              }
            }
          }
          if (pseudo) {
            nummrnapseudo++;
          }
        }
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = NULL;
  vsp->bsp = bsp;

  if (cds_prod_head != NULL) {
    cds_prod_head = ValNodeSort (cds_prod_head, SortVnpByString);
    cds_products_unique = StringListIsUnique (cds_prod_head);
  }
  if (mrna_prod_head != NULL) {
    mrna_prod_head = ValNodeSort (mrna_prod_head, SortVnpByString);
    mrna_products_unique = StringListIsUnique (mrna_prod_head);;
  }

  if (numcds > 0 && nummrna > 1) {
    if (numcdsproducts + numcdspseudo == numcds &&
        (nummrnaproducts + nummrnapseudo == nummrna || nummrnaproducts == 0) &&
        cds_products_unique && mrna_products_unique) {
      suppress_duplicate_messages = TRUE;
    }
    if (numcdsproducts > 0 && numcdsproducts + numcdspseudo != numcds) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d CDS features have %d product references",
                (int) numcds, (int) numcdsproducts);
    }
    if (numcdsproducts > 0 && (! cds_products_unique)) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "CDS products are not unique");
    }
    if (nummrnaproducts > 0 && nummrnaproducts + nummrnapseudo != nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d mRNA features have %d product references",
                (int) nummrna, (int) nummrnaproducts);
    }
    if (nummrnaproducts > 0 && (! mrna_products_unique)) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "mRNA products are not unique");
    }
    /*
    if (numcds > nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d CDS features and only %d mRNA features",
                (int) numcds, (int) nummrna);
    } else if (numcds < nummrna) {
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureProductInconsistency, "%d mRNA features and only %d CDS features",
                (int) nummrna, (int) numcds);
    }
    */
  }

  ValNodeFreeData (cds_prod_head);
  ValNodeFreeData (mrna_prod_head);

  /*
  SeqEntryToBioSource (vsp->sep, NULL, NULL, 0, &biop);
  */
  BioseqToGeneticCode (bsp, NULL, NULL, NULL, NULL, 0, &biop);
  if (biop != NULL) {
    orp = biop->org;
    if (orp != NULL) {
      /* curated fly source still has duplicate features */
      if (StringICmp (orp->taxname, "Drosophila melanogaster") == 0) {
        isDrosophila = TRUE;
      }
      onp = orp->orgname;
      if (onp != NULL) {
        if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0) {
          isViral = TRUE;
        }
      }
    }
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    leave = TRUE;
    if (last != NULL) {
      ivalssame = FALSE;
      if (fcontext.left == left && fcontext.right == right && fcontext.featdeftype == featdeftype) {
        if (fcontext.strand == strand || strand == Seq_strand_unknown || fcontext.strand == Seq_strand_unknown) {
          ivalssame = TRUE;
          if (fcontext.numivals != numivals || fcontext.ivals == NULL || ivals == NULL) {
            ivalssame = FALSE;
          } else {
            for (i = 0, j = 0; i < numivals; i++, j += 2) {
              if (fcontext.ivals[j] != ivals[j]) {
                ivalssame = FALSE;
              }
              if (fcontext.ivals[j + 1] != ivals[j + 1]) {
                ivalssame = FALSE;
              }
            }
          }
          if (ivalssame &&      /* StringICmp (fcontext.label, label) == 0 && */
              (fcontext.sap == sap || (fcontext.sap->desc == NULL && sap->desc == NULL) || DescsSame (fcontext.sap->desc, sap->desc))) {
            if (gcp != NULL) {
              gcp->itemID = fcontext.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = sfp;
            severity = SEV_ERROR;
            samelabel = TRUE;
            if (StringICmp (fcontext.label, label) != 0 || StringICmp (sfp->comment, comment) != 0) {
              samelabel = FALSE;
            }
            if (GBQualsDiffer (sfp->qual, last->qual)) {
              samelabel = FALSE;
            }
            if (featdeftype == FEATDEF_PUB ||
                featdeftype == FEATDEF_REGION || featdeftype == FEATDEF_misc_feature || featdeftype == FEATDEF_STS || featdeftype == FEATDEF_variation) {
              severity = SEV_WARNING;
            } else {
              if (isGPSorNTorNCorNGorNW || GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
                isGPSorNTorNCorNGorNW = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else if (isGenBankAccn || IsGenBankAccn (vsp->sep, sfp->location)) {
                isGenBankAccn = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else if (isGeneralAccn || IsGeneralAccn (vsp->sep, sfp->location)) {
                isGeneralAccn = TRUE;
                if (! isCuratedFlybase) {
                  if (isDrosophila) {
                    isCuratedFlybase = TRUE;
                  }
                }
                if (isCuratedFlybase) {
                  /* curated fly source still has duplicate features */
                  severity = SEV_WARNING;
                }
              } else {
                severity = SEV_WARNING;
              }
            }
            /* if different CDS frames, lower to warning */
            if (sfp->data.choice == SEQFEAT_CDREGION) {
              crp = (CdRegionPtr) sfp->data.value.ptrvalue;
              if (crp != NULL) {
                if (frame > 1 || crp->frame > 1) {
                  if (frame != crp->frame) {
                    severity = SEV_WARNING;
                    if (! samelabel) {
                      if (numivals == 1 && left == 0 && right == bsp->length - 1 && partialL && partialR) {
                        /* skip full length partial CDS features in different frames with different products */
                        severity = SEV_NONE;
                      }
                    }
                  }
                }
              }
            }
            if (isGPSorNTorNCorNGorNW || GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
              isGPSorNTorNCorNGorNW = TRUE;
              severity = SEV_WARNING;
            }
            if (FlybaseDbxrefs (last->dbxref) || FlybaseDbxrefs (sfp->dbxref)) {
              severity = SEV_ERROR;
            }
            if (featdeftype == FEATDEF_repeat_region) {
              severity = SEV_WARNING;
            }
            if (featdeftype == FEATDEF_SITE || featdeftype == FEATDEF_BOND) {
              severity = SEV_WARNING;
            }
            if (severity == SEV_NONE) {
              /* skip full length partial CDS features in different frames with different products */
            } else if (featdeftype == FEATDEF_REGION && DifferentDbxrefs (last->dbxref, sfp->dbxref)) {
              /* do not report if both have dbxrefs and they are different */
            } else if (featdeftype == FEATDEF_variation && ReplaceQualsDiffer (sfp->qual, last->qual)) {
              /* do not report if both have replace quals and they are different */
            } else if (CDSsLinkedToDifferentMRNAs (sfp, last)) {
              /* do not report if CDSs are linked to two different mRNAs */
            } else if (MRNAsLinkedToDifferentCDSs (sfp, last)) {
              /* do not report if mRNAs are linked to two different CDSs */
            } else if (fcontext.sap == sap) {
              if (samelabel) {
                ValidErr (vsp, severity, ERR_SEQ_FEAT_FeatContentDup, "Duplicate feature");
              } else if (featdeftype != FEATDEF_PUB) {
                if (fcontext.partialL != partialL || fcontext.partialR != partialR) {
                  /* do not report if partial flags are different */
                } else {
                  if (suppress_duplicate_messages && (featdeftype == FEATDEF_CDS || featdeftype == FEATDEF_mRNA) && HaveUniqueFeatIDXrefs (xref, sfp->xref)) {
                    /* do not report CDS or mRNA if every one has a unique product and unique featID xrefs */
                  } else if (featdeftype == FEATDEF_GENE &&
                             StringStr (sfp->except_text, "dicistronic gene") != NULL &&
                             StringStr (except_text, "dicistronic gene") != NULL &&
                             isCuratedFlybase) {
                    /* do not report genes marked dicistronic */
                  } else {
                    if (featdeftype == FEATDEF_GENE && isViral && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_CDS && isViral && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_mRNA && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_GENE && (sfp->partial || last->partial)) {
                      severity = SEV_WARNING;
                    }
                    if (featdeftype == FEATDEF_GENE && sfp->pseudo && last->pseudo) {
                      severity = SEV_WARNING;
                    }
                    if (severity == SEV_ERROR && featdeftype == FEATDEF_CDS && mip != NULL) {
                      if (mip->tech >= MI_TECH_htgs_1 && mip->tech <= MI_TECH_htgs_3) {
                        severity = SEV_WARNING;
                      }
                    }
                    if (fcontext.seqfeattype == SEQFEAT_IMP) {
                      severity = SEV_WARNING;
                    }
                    ValidErr (vsp, severity, ERR_SEQ_FEAT_DuplicateFeat, "Features have identical intervals, but labels differ");
                  }
                }
              }
            } else {
              if (samelabel) {
                ValidErr (vsp, severity, ERR_SEQ_FEAT_FeatContentDup, "Duplicate feature (packaged in different feature table)");
              } else if (featdeftype != FEATDEF_PUB) {
                if (suppress_duplicate_messages && (featdeftype == FEATDEF_CDS || featdeftype == FEATDEF_mRNA) && HaveUniqueFeatIDXrefs (xref, sfp->xref)) {
                  /* do not report CDS or mRNA if every one has a unique product and unique featID xrefs */
                } else {
                  ValidErr (vsp, /* severity */ SEV_WARNING, ERR_SEQ_FEAT_DuplicateFeat, "Features have identical intervals, but labels differ (packaged in different feature table)");
                }
              }
            }
            vsp->sfp = NULL;
            if (gcp != NULL) {
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
      if (fcontext.featdeftype == FEATDEF_mat_peptide_aa ||
          fcontext.featdeftype == FEATDEF_sig_peptide_aa || fcontext.featdeftype == FEATDEF_transit_peptide_aa) {
        if (featdeftype == FEATDEF_mat_peptide_aa || featdeftype == FEATDEF_sig_peptide_aa || featdeftype == FEATDEF_transit_peptide_aa) {
          if (fcontext.left <= right && NotPeptideException (sfp, last)) {
            if (gcp != NULL) {
              gcp->itemID = fcontext.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            buf [0] = '\0';
            cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
            if (cds != NULL) {
              nbsp = BioseqFindFromSeqLoc (cds->location);
              if (nbsp != NULL) {
                SeqIdWrite (nbsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
              }
            }
            vsp->descr = NULL;
            vsp->sfp = sfp;
            if (StringDoesHaveText (buf)) {
              ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat,
                        "Signal, Transit, or Mature peptide features overlap (parent CDS is on %s)", buf);
            } else {
              ValidErr (vsp, overlapPepSev, ERR_SEQ_FEAT_OverlappingPeptideFeat, "Signal, Transit, or Mature peptide features overlap");
            }
            vsp->sfp = NULL;
            if (gcp != NULL) {
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
    }
    if (leave) {
      last = sfp;
      left = fcontext.left;
      right = fcontext.right;
      label = fcontext.label;
      comment = sfp->comment;
      strand = fcontext.strand;
      partialL = fcontext.partialL;
      partialR = fcontext.partialR;
      featdeftype = fcontext.featdeftype;
      numivals = fcontext.numivals;
      ivals = fcontext.ivals;
      sap = fcontext.sap;
      xref = sfp->xref;
      except_text = sfp->except_text;
      frame = 0;
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        if (crp != NULL) {
          frame = crp->frame;
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  lastLabel = NULL;
  lastsfp = NULL;
  lastgrp = NULL;
  grp = NULL;
  sfp = SeqMgrGetNextFeatureByLabel (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
  while (sfp != NULL) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    label = fcontext.label;
    if (lastLabel != NULL) {
      message = NULL;
      if (StringCmp (lastLabel, label) == 0) {
        message = "Colliding names in gene features";
      } else if (StringICmp (lastLabel, label) == 0) {
        message = "Colliding names (with different capitalization) in gene features";
      }
      if (message != NULL && (ReportGeneCollision (grp, lastgrp))) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CollidingGeneNames, "%s", message);
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    lastLabel = label;
    lastsfp = sfp;
    lastgrp = grp;
    sfp = SeqMgrGetNextFeatureByLabel (bsp, sfp, SEQFEAT_GENE, 0, &fcontext);
  }

  lastLabel = NULL;
  sfp = SeqMgrGetNextGeneByLocusTag (bsp, NULL, &fcontext);
  while (sfp != NULL) {
    label = NULL;
    if (sfp->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        label = grp->locus_tag;
      }
    }
    if (lastLabel != NULL) {
      message = NULL;
      if (StringCmp (lastLabel, label) == 0) {
        message = "Colliding locus_tags in gene features";
      } else if (StringICmp (lastLabel, label) == 0) {
        message = "Colliding locus_tags (with different capitalization) in gene features";
      }
      if (message != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CollidingLocusTags, "%s", message);
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    lastLabel = label; 
    sfp = SeqMgrGetNextGeneByLocusTag (bsp, sfp, &fcontext);
  }

  /* do UTR vs. CDS check on genomic if only one CDS, still need separate minus strand logic */
  cdscount = 0;
  genecount = 0;
  utr5count = 0;
  utr3count = 0;
  cdsgene = NULL;
  utr5gene = NULL;
  utr3gene = NULL;
  strand = Seq_strand_plus;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL /* && cdscount < 2 && genecount < 2 */) {
    if (sfp->idx.subtype == FEATDEF_CDS) {
      strand = fcontext.strand;
      cdscount++;
      cdsgene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    } else if (sfp->idx.subtype == FEATDEF_GENE) {
      genecount++;
    } else if (sfp->idx.subtype == FEATDEF_5UTR) {
      utr5count++;
      utr5gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    } else if (sfp->idx.subtype == FEATDEF_3UTR) {
      utr3count++;
      utr3gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
  if (bvsp->is_mrna || cdscount == 1 && genecount < 2) {
    if (bvsp->is_mrna) {
      strand = Seq_strand_plus;
    }
    fiveUTRright = 0;
    cdsRight = 0;
    threeUTRright = 0;
    firstCDS = TRUE;

    if (strand == Seq_strand_minus) {

      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (sfp->idx.subtype == FEATDEF_3UTR && utr3count < 2) {
          if (fcontext.strand != Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "3'UTR is not on minus strand");
          } else if (threeUTRright > 0) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "Previous 3'UTR does not abut next 3'UTR");
            }
          }
          threeUTRright = fcontext.right;
        } else if (sfp->idx.subtype == FEATDEF_CDS) {
          cdsRight = fcontext.right;
          if (threeUTRright > 0 && firstCDS) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "CDS does not abut 3'UTR");
            }
          }
          firstCDS = FALSE;
        } else if (sfp->idx.subtype == FEATDEF_5UTR && utr5count < 2) {
          if (fcontext.strand != Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR is not on minus strand");
          } else if (cdsRight > 0) {
            if (cdsRight + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR does not abut CDS");
            }
          }
          threeUTRright = fcontext.right;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }

    } else {

      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        if (sfp->idx.subtype == FEATDEF_5UTR && utr5count < 2) {
          if (fcontext.strand == Seq_strand_minus) {
            if (genecount > 1 && cdsgene != NULL && utr5gene != NULL && cdsgene != utr5gene) {
              /* ignore */
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR is not on plus strand");
            }
          }
          fiveUTRright = fcontext.right;
        } else if (sfp->idx.subtype == FEATDEF_CDS) {
          cdsRight = fcontext.right;
          if (fiveUTRright > 0 && firstCDS) {
            if (fiveUTRright + 1 != fcontext.left) {
              if (genecount > 1 && cdsgene != NULL && utr5gene != NULL && cdsgene != utr5gene) {
                /* ignore */
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "5'UTR does not abut CDS");
              }
            }
          }
          firstCDS = FALSE;
        } else if (sfp->idx.subtype == FEATDEF_3UTR && utr3count < 2) {
          if (fcontext.strand == Seq_strand_minus) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "3'UTR is not on plus strand");
          } else if (threeUTRright > 0) {
            if (threeUTRright + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "Previous 3'UTR does not abut next 3'UTR");
            }
          } else if (cdsRight > 0) {
            if (cdsRight + 1 != fcontext.left) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotAbutCDS, "CDS does not abut 3'UTR");
            }
            if (bvsp->is_mrna && cdscount == 1 && utr3count == 1 && fcontext.right != bsp->length - 1) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UTRdoesNotExtendToEnd, "3'UTR does not extend to end of mRNA");
            }
          }
          threeUTRright = fcontext.right;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }
    }
  }

  if (! bvsp->is_mrna) {
    last = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &fcontext);
    if (last != NULL) {
      lastrnatype = WhichRNA (last);
      left = fcontext.left;
      right = fcontext.right;
      strand = fcontext.strand;
      partialL = fcontext.partialL;
      partialR = fcontext.partialR;
      sfp = SeqMgrGetNextFeature (bsp, last, SEQFEAT_RNA, 0, &fcontext);
      while (sfp != NULL) {
        thisrnatype = WhichRNA (sfp);
        if (fcontext.strand == strand || (strand != Seq_strand_minus && fcontext.strand != Seq_strand_minus)) {
          if (lastrnatype != 0 && thisrnatype != 0) {
            if (right + 1 < fcontext.left) {
              /* gap */
              if (BaseRangeIsVirtual (bsp, right + 1, fcontext.left)) {
                /* ignore if abuts gap */
              } else if (strand == Seq_strand_minus) {
                if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  /* okay in mitochondria */
                } else if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_2 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ITSdoesNotAbutRRNA, "ITS does not abut adjacent rRNA component");
                }
              } else {
                if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  /* okay in mitochondria */
                } else if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_1 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ITSdoesNotAbutRRNA, "ITS does not abut adjacent rRNA component");
                }
              }
            } else if (right + 1 > fcontext.left) {
              /* overlaps */
              if (strand == Seq_strand_minus) {
                if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOverlap, "tRNA overlaps adjacent rRNA component");
                } else if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_2 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOverlap, "ITS overlaps adjacent rRNA component");
                }
              } else {
                if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  /* okay in mitochondria */
                } else if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype == INTERNAL_SPACER_1 || thisrnatype == INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype == INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype == MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOverlap, "ITS overlaps adjacent rRNA component");
                }
              }
            } else {
              /* abuts */
              if (strand == Seq_strand_minus) {
                if (lastrnatype == thisrnatype && partialL && fcontext.partialR && bsp->repr == Seq_repr_seg) {
                  /* okay in segmented set */
                } else if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == LEFT_RIBOSOMAL_SUBUNIT)) {
                  /* okay in mitochondria */
                } else if ((lastrnatype == RIGHT_RIBOSOMAL_SUBUNIT && (thisrnatype != INTERNAL_SPACER_2 && thisrnatype != INTERNAL_SPACER_X)) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype != MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype != INTERNAL_SPACER_1) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype != LEFT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype != LEFT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOrder, "Problem with order of abutting rRNA components");
                }
              } else {
                if (lastrnatype == thisrnatype && partialR && fcontext.partialL && bsp->repr == Seq_repr_seg) {
                  /* okay in segmented set */
                } else if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && thisrnatype == TRANSFER_RNA) ||
                    (lastrnatype == TRANSFER_RNA && thisrnatype == RIGHT_RIBOSOMAL_SUBUNIT)) {
                  /* okay in mitochondria */
                } else if ((lastrnatype == LEFT_RIBOSOMAL_SUBUNIT && (thisrnatype != INTERNAL_SPACER_1 && thisrnatype != INTERNAL_SPACER_X)) ||
                    (lastrnatype == MIDDLE_RIBOSOMAL_SUBUNIT && thisrnatype != INTERNAL_SPACER_2) ||
                    (lastrnatype == INTERNAL_SPACER_1 && thisrnatype != MIDDLE_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_2 && thisrnatype != RIGHT_RIBOSOMAL_SUBUNIT) ||
                    (lastrnatype == INTERNAL_SPACER_X && thisrnatype != RIGHT_RIBOSOMAL_SUBUNIT)) {
                  if (gcp != NULL) {
                    gcp->itemID = fcontext.itemID;
                    gcp->thistype = OBJ_SEQFEAT;
                  }
                  vsp->descr = NULL;
                  vsp->sfp = sfp;
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadRRNAcomponentOrder, "Problem with order of abutting rRNA components");
                }
              }
            }
          }
        }
        last = sfp;
        left = fcontext.left;
        right = fcontext.right;
        strand = fcontext.strand;
        partialL = fcontext.partialL;
        partialR = fcontext.partialR;
        lastrnatype = thisrnatype;
        sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, 0, &fcontext);
      }
    }
  }

  vsp->sfp = NULL;
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }

  mrna = SeqMgrGetRNAgivenProduct (bsp, &fcontext);
  if (mrna != NULL) {
    genomicgrp = SeqMgrGetGeneXref (mrna);
    if (genomicgrp == NULL) {
      gene = SeqMgrGetOverlappingGene (mrna->location, NULL);
      if (gene != NULL) {
        genomicgrp = (GeneRefPtr) gene->data.value.ptrvalue;
      }
    }
    if (genomicgrp != NULL) {
      gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
      if (gene != NULL) {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (grp != NULL) {
          if (StringCmp (grp->locus, genomicgrp->locus) != 0 ||
              StringCmp (grp->allele, genomicgrp->allele) != 0 ||
              StringCmp (grp->desc, genomicgrp->desc) != 0 ||
              StringCmp (grp->locus_tag, genomicgrp->locus_tag) != 0) {
            if (gcp != NULL) {
              gcp->itemID = fcontext.itemID;
              gcp->thistype = OBJ_SEQFEAT;
            }
            vsp->descr = NULL;
            vsp->sfp = gene;
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenesInconsistent, "Gene on mRNA bioseq does not match gene on genomic bioseq");
            if (gcp != NULL) {
              gcp->itemID = olditemid;
              gcp->thistype = olditemtype;
            }
          }
        }
      }
    }
  }

  if (ISA_na (bsp->mol)) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_gap, &fcontext);
    while (sfp != NULL) {
      estimated_length = 0;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringICmp (gbq->qual, "estimated_length") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        if (StringICmp (gbq->val, "unknown") == 0) continue;
        if (sscanf (gbq->val, "%ld", &val) == 1) {
          estimated_length = val;
        }
      }
      if (StreamCacheSetup (NULL, sfp->location, EXPAND_GAPS_TO_DASHES, &sc)) {
        dashes = 0;
        Ns = 0;
        realBases = 0;
        while ((res = StreamCacheGetResidue (&sc)) != '\0') {
          if (IS_LOWER (res)) {
            res = TO_UPPER (res);
          }
          if (res == '-') {
            dashes++;
          } else if (res == 'N') {
            Ns++;
          } else {
            realBases++;
          }
        }
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        loclen = SeqLocLen (sfp->location);
        if (estimated_length > 0 && estimated_length != loclen) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature estimated_length %ld does not match %ld feature length",
                    (long) estimated_length, (long) loclen);
        } else if (realBases > 0 && Ns > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld real bases and %ld Ns", (long) realBases, (long) Ns);
        } else if (realBases > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld real bases", (long) realBases);
        } else if (Ns > 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature over %ld Ns", (long) Ns);
        } else if (estimated_length > 0 && dashes != estimated_length) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_GapFeatureProblem, "Gap feature estimated_length %ld does not match %ld gap characters",
                    (long) estimated_length, (long) dashes);
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_gap, &fcontext);
    }
  }
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if (ISA_aa (bsp->mol)) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (sfp != NULL) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (SeqLocStrand (slp) == Seq_strand_minus) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MinusStrandProtein, "Feature on protein indicates negative strand");
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
    }
  }
  if (gcp != NULL) {
    gcp->itemID = olditemid;
    gcp->thistype = olditemtype;
  }
  vsp->descr = NULL;
  vsp->sfp = NULL;

  lastbiop = NULL;
  lastsfp = NULL;
  numBadFullSource = 0;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  if (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      showBadFullSource = TRUE;
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        if (biop != NULL) {
          if (biop->is_focus) {
            showBadFullSource = FALSE;
          }
          for (sbsp = biop->subtype; sbsp != NULL; sbsp = sbsp->next) {
            if (sbsp->subtype == SUBSRC_transgenic) {
              showBadFullSource = FALSE;
            }
          }
        }
      }
      if (showBadFullSource) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Source feature is full length, should be descriptor");
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
  }
  /* and fall through to continue testing first and remaining source features */
  while (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      numBadFullSource++;
      if (numBadFullSource > 1) {
        if (gcp != NULL) {
          gcp->itemID = fcontext.itemID;
          gcp->thistype = OBJ_SEQFEAT;
        }
        vsp->descr = NULL;
        vsp->sfp = sfp;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Multiple full-length source features, should only be one if descriptor is transgenic");
        vsp->sfp = NULL;
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    if (biop != NULL && lastbiop != NULL) {
      if (lastsfp != NULL) {
        if (StringDoesHaveText (lastsfp->comment) && StringDoesHaveText (sfp->comment) && StringICmp (lastsfp->comment, sfp->comment) != 0) {
          /* different comments, so ignore */
        } else if (IsIdenticalBioSource (biop, lastbiop) && (! bvsp->is_artificial)) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleEquivBioSources, "Multiple equivalent source features should be combined into one multi-interval feature");
          vsp->sfp = NULL;
          if (gcp != NULL) {
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
      }
    }
    lastbiop = biop;
    lastsfp = sfp;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
  }

  lastpdp = NULL;
  lastsfp = NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
  if (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1 && fcontext.numivals == 1) {
      if (gcp != NULL) {
        gcp->itemID = fcontext.itemID;
        gcp->thistype = OBJ_SEQFEAT;
      }
      vsp->descr = NULL;
      vsp->sfp = sfp;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadFullLengthFeature, "Publication feature is full length, should be descriptor");
      vsp->sfp = NULL;
      if (gcp != NULL) {
        gcp->itemID = olditemid;
        gcp->thistype = olditemtype;
      }
    }
  }
  while (sfp != NULL) {
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    if (pdp != NULL && lastpdp != NULL) {
      if (lastsfp != NULL) {
        if (StringDoesHaveText (lastsfp->comment) && StringDoesHaveText (sfp->comment) && StringICmp (lastsfp->comment, sfp->comment) != 0) {
          /* different comments, so ignore */
        } else if (IsIdenticalPublication (pdp, lastpdp)) {
          if (gcp != NULL) {
            gcp->itemID = fcontext.itemID;
            gcp->thistype = OBJ_SEQFEAT;
          }
          vsp->descr = NULL;
          vsp->sfp = sfp;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleEquivPublications, "Multiple equivalent publication features should be combined into one multi-interval feature");
          vsp->sfp = NULL;
          if (gcp != NULL) {
            gcp->itemID = olditemid;
            gcp->thistype = olditemtype;
          }
        }
      }
    }
    lastpdp = pdp;
    lastsfp = sfp;
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext);
  }

  SeqMgrExploreDescriptors (bsp, (Pointer) bvsp, ValidateSeqDescrIndexed, NULL);

  omdp = ObjMgrGetData (gcp->entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL) {
      sbp = ssp->sub;
      if (sbp != NULL) {
        bvsp->got_a_pub = TRUE;
        csp = sbp->cit;
        /* csp = (CitSubPtr) gcp->thisitem; */
        if (gcp != NULL) {
          gcp->itemID = 1;
          gcp->thistype = OBJ_SEQSUB_CIT;
        }
        vsp->descr = NULL;
        vsp->sfp = NULL;
        vsp->bssp = NULL;
        ValidateCitSub (vsp, csp);
        if (gcp != NULL) {
          gcp->itemID = olditemid;
          gcp->thistype = olditemtype;
        }
      }
    }
  }

  ValidateCDSmRNAmatch (vsp, bsp, numgene, numcds, nummrna, suppress_duplicate_messages);

  if (vsp->locusTagGeneralMatch) {
    ValidateLocusTagGeneral (vsp, bsp);
  }

  if (ISA_na (bsp->mol) && SeqMgrGetParentOfPart (bsp, NULL) == NULL) {
    LookForMultipleUnpubPubs (vsp, gcp, bsp);
  }

  if (bsp->repr == Seq_repr_delta && ISA_na (bsp->mol)) {
    CheckBioseqForFeatsInGap (bsp, vsp);
  }

  return TRUE;
}

static Boolean ValidateBioseqContextGather (GatherContextPtr gcp)
{
  ValidStructPtr  vsp;
  BioseqValidStrPtr bvsp;
  CitSubPtr       csp;

  bvsp = (BioseqValidStrPtr) (gcp->userdata);
  vsp = bvsp->vsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->gcp = gcp;               /* needed for ValidErr */

  switch (gcp->thistype) {
  case OBJ_SEQFEAT:
    ValidateSeqFeatContext (gcp);
    break;
  case OBJ_SEQDESC:
    ValidateSeqDescrContext (gcp);
    break;
  case OBJ_SEQSUB_CIT:
    bvsp->got_a_pub = TRUE;
    csp = (CitSubPtr) gcp->thisitem;
    ValidateCitSub (vsp, csp);
    break;
  default:
    break;
  }
  return TRUE;
}


static ValNodePtr ListFeaturesContainedInLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;
  Int4              cmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    cmp = SeqLocCompare (sfp->location, slp);
    if (cmp == SLC_A_EQ_B || cmp == SLC_A_IN_B)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


/*****************************************************************************
*   FindMultiGeneOverlaps (BioseqPtr bsp, ValidStructPtr vsp)
*
*      This function reports genes that overlap two or more other genes.
*****************************************************************************/
static void FindMultiGeneOverlaps (BioseqPtr bsp, ValidStructPtr vsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  ValNodePtr overlap_list;
  Int4       num_genes_at_loc;
  GatherContextPtr  gcp;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;

  if (bsp == NULL || vsp == NULL || vsp->gcp == NULL) {
    return;
  }

  gcp = vsp->gcp;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &context);
  while (sfp != NULL) {
    overlap_list = ListFeaturesContainedInLocation (bsp, sfp->location, SEQFEAT_GENE, 0);
    num_genes_at_loc = ValNodeLen (overlap_list);
    if (num_genes_at_loc > 2) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = sfp;

      gcp->entityID = sfp->idx.entityID;
      gcp->itemID = sfp->idx.itemID;
      gcp->thistype = OBJ_SEQFEAT;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleGeneOverlap, "Gene contains %d other genes", num_genes_at_loc - 1);
    }
    overlap_list = ValNodeFree (overlap_list);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, 0, &context);
  }

  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;
}


static void ValidateTSASequenceForNs (BioseqPtr bsp, ValidStructPtr vsp)
{
  Int4 total = 0, max_stretch = 0;
  GatherContextPtr  gcp;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;
  Int4            percent_N;

  gcp = vsp->gcp;

  oldEntityID = gcp->entityID;
  oldItemID = gcp->itemID;
  oldItemtype = gcp->thistype;

  if (IsTSA (bsp)) {
    CountNsInSequence (bsp, &total, &max_stretch, FALSE);
    percent_N = (total * 100) / bsp->length;
    if (percent_N > 5) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContent, "Sequence contains %d percent Ns", percent_N);
    }
    if (max_stretch > 5) {
      vsp->bsp = bsp;
      vsp->descr = NULL;
      vsp->sfp = NULL;
      gcp->entityID = bsp->idx.entityID;
      gcp->itemID = bsp->idx.itemID;
      gcp->thistype = OBJ_BIOSEQ;

      ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_HighNContent, "Sequence has a stretch of %d Ns", max_stretch);
    }
  }
  gcp->entityID = oldEntityID;
  gcp->itemID = oldItemID;
  gcp->thistype = oldItemtype;
}


/*****************************************************************************
*
*   ValidateBioseqContext(gcp)
*      Validate one Bioseq for descriptors, features, and context
*      This is done as a second Gather, focussed on the Bioseq in
*        question.
*
*****************************************************************************/
static void ValidateBioseqContext (GatherContextPtr gcp)
{
  size_t          acclen;
  ValidStructPtr  vsp;
  BioseqPtr       bsp;
  GatherScope     gs;
  BioseqValidStr  bvs;
  SeqFeatPtr      sfp;
  ValNode         fake_whole;
  SeqIdPtr        sip;
  ValNodePtr      vnp = NULL;
  MolInfoPtr      mip = NULL;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioseqContextPtr bcp;
  Uint2           oldEntityID, oldItemtype;
  Uint4           oldItemID;
  Uint2           mipEntityID = 0, mipItemtype = 0;
  Uint4           mipItemID = 0;
  ObjMgrDataPtr   omdp;
  BioseqPtr       parent;
  PatentSeqIdPtr  psip;
  IdPatPtr        ipp;
  Boolean         isPDB = FALSE;
  Boolean         is_wgs = FALSE;
  Boolean         is_gb = FALSE;
  Boolean         is_ch_or_cm = FALSE;
  Boolean         is_nc = FALSE;
  Boolean         is_local = FALSE;
  Boolean         is_local_only = TRUE;
  Boolean         is_neg_strand_virus = FALSE;
  Boolean         is_ambisense_virus = FALSE;
  Boolean         is_synthetic = FALSE;
  Boolean         is_transgenic = FALSE;
  Boolean         has_cds = FALSE;
  ErrSev          sev;
  SubSourcePtr    ssp;
  CharPtr         str;
  TextSeqIdPtr    tsip;
  BioSourcePtr    biop;
  OrgRefPtr       orp;
  OrgNamePtr      onp;
  OrgModPtr       omp;
  /*
  Char            buf1[255];
  */

  vsp = (ValidStructPtr) (gcp->userdata);
  bsp = (BioseqPtr) (gcp->thisitem);
  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) (gcp->parentitem);

  MemSet (&gs, 0, sizeof (GatherScope));
  fake_whole.choice = SEQLOC_WHOLE;
  sip = SeqIdFindBest (bsp->id, 0);

  fake_whole.data.ptrvalue = sip;

  fake_whole.next = NULL;
  gs.target = &fake_whole;
  gs.get_feats_location = TRUE;
  gs.nointervals = TRUE;
  MemSet ((Pointer) (gs.ignore), (int) TRUE, (size_t) (sizeof (Boolean) * OBJ_MAX));
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SUBMIT_BLOCK] = FALSE;
  gs.ignore[OBJ_SEQSUB_CIT] = FALSE;

  gs.scope = vsp->sep;

  MemSet (&bvs, 0, sizeof (BioseqValidStr));
  bvs.vsp = vsp;

  /* now looking for molinfo on every bioseq (okay on segset) */
  if (bsp != NULL) {
    vnp = NULL;
    if (vsp->useSeqMgrIndexes) {
      vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (vnp != NULL) {
        mipEntityID = dcontext.entityID;
        mipItemID = dcontext.itemID;
        mipItemtype = OBJ_SEQDESC;
      }
    } else {
      bcp = BioseqContextNew (bsp);
      vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
      BioseqContextFree (bcp);
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
    }
    vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (vnp != NULL) {
      biop = (BioSourcePtr) vnp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          if (StringICmp (orp->taxname, "Human immunodeficiency virus") == 0 ||
              StringICmp (orp->taxname, "Human immunodeficiency virus 1") == 0 ||
              StringICmp (orp->taxname, "Human immunodeficiency virus 2") == 0) {
            ValidateLocationForHIV (vsp, biop, bsp);
          } else if (StringICmp (orp->taxname, "synthetic construct") == 0) {
            bvs.is_syn_constr = TRUE;
          } else if (StringISearch (orp->taxname, "vector") != NULL) {
            bvs.is_syn_constr = TRUE;
          }
          onp = orp->orgname;
          if (onp != NULL) {
            if (StringICmp (onp->div, "SYN") == 0) {
              bvs.is_syn_constr = TRUE;
              is_synthetic = TRUE;
            }
            if (StringISearch (onp->lineage, "artificial sequences") != NULL) {
              bvs.is_syn_constr = TRUE;
            }
            if (StringISearch (onp->lineage, "negative-strand viruses") != NULL) {
              is_neg_strand_virus = TRUE;
            }
            if (StringISearch (onp->lineage, "Arenavirus") != NULL ||
                StringISearch (onp->lineage, "Phlebovirus") != NULL ||
                StringISearch (onp->lineage, "Tospovirus") != NULL ||
                StringISearch (onp->lineage, "Tenuivirus") != NULL) {
              is_ambisense_virus = TRUE;
            }
            for (omp = onp->mod; omp != NULL; omp = omp->next) {
              if (omp->subtype == ORGMOD_other) {
                if (mip != NULL && (StringICmp (omp->subname, "cRNA") == 0)) {
                  oldEntityID = gcp->entityID;
                  oldItemID = gcp->itemID;
                  oldItemtype = gcp->thistype;

                  gcp->entityID = dcontext.entityID;
                  gcp->itemID = dcontext.itemID;
                  gcp->thistype = OBJ_SEQDESC;

                  if (mip->biomol == MOLECULE_TYPE_CRNA) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note redundant with molecule type");
                  } else {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note conflicts with molecule type");
                  }

                  gcp->entityID = oldEntityID;
                  gcp->itemID = oldItemID;
                  gcp->thistype = oldItemtype;
                }
              }
            }
            if (mip != NULL && mip->biomol == MOLECULE_TYPE_GENOMIC && bsp->mol == MOLECULE_CLASS_DNA) {
              if (StringNICmp (onp->lineage, "Viruses; ", 9) == 0 && StringISearch (onp->lineage, "no DNA stage") != NULL) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Genomic DNA viral lineage indicates no DNA stage");
              }
            }
          }
        }
        if (biop->origin == ORG_ARTIFICIAL || biop->origin == ORG_SYNTHETIC) {
          bvs.is_artificial = TRUE;
        }
        if (biop->origin == ORG_MUT || biop->origin == ORG_ARTIFICIAL || biop->origin == ORG_SYNTHETIC) {
          is_synthetic = TRUE;
        }
        for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
          if (ssp->subtype == SUBSRC_transgenic) {
            is_transgenic = TRUE;
          } else if (ssp->subtype == SUBSRC_other) {
            if (mip != NULL && (StringICmp (ssp->name, "cRNA") == 0)) {
              oldEntityID = gcp->entityID;
              oldItemID = gcp->itemID;
              oldItemtype = gcp->thistype;

              gcp->entityID = dcontext.entityID;
              gcp->itemID = dcontext.itemID;
              gcp->thistype = OBJ_SEQDESC;

              if (mip->biomol == MOLECULE_TYPE_CRNA) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note redundant with molecule type");
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "cRNA note conflicts with molecule type");
              }

              gcp->entityID = oldEntityID;
              gcp->itemID = oldItemID;
              gcp->thistype = oldItemtype;
            }
          }
        }
        if (is_transgenic && ISA_na (bsp->mol)) {
          if (SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext) == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TransgenicProblem, "Transgenic source descriptor requires presence of source feature");
          }
        }
      }
    }
  }

  if (is_neg_strand_virus && mip != NULL) {
    oldEntityID = gcp->entityID;
    oldItemID = gcp->itemID;
    oldItemtype = gcp->thistype;

    gcp->entityID = mipEntityID;
    gcp->itemID = mipItemID;
    gcp->thistype = mipItemtype;

    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      has_cds = TRUE;
      if (SeqLocStrand (sfp->location) == Seq_strand_minus) {
        if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with minus strand CDS should be genomic");
        }
      } else {
        if (mip->biomol != MOLECULE_TYPE_MRNA && mip->biomol != MOLECULE_TYPE_CRNA && (! is_ambisense_virus) && (! is_synthetic)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with plus strand CDS should be mRNA or cRNA");
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
    if (! has_cds) {
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
      while (sfp != NULL) {
        if (StringISearch (sfp->comment, "nonfunctional") != NULL) {
          if (SeqLocStrand (sfp->location) == Seq_strand_minus) {
            if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with nonfunctional minus strand misc_feature should be mRNA or cRNA");
            }
          } else {
            if (mip->biomol != MOLECULE_TYPE_MRNA && mip->biomol != MOLECULE_TYPE_CRNA && (! is_ambisense_virus)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BioSourceInconsistency, "Negative-strand virus with nonfunctional plus strand misc_feature should be mRNA or cRNA");
            }
          }
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext);
      }
    }

    gcp->entityID = oldEntityID;
    gcp->itemID = oldItemID;
    gcp->thistype = oldItemtype;
  }

  bvs.is_mrna = FALSE;
  bvs.is_prerna = FALSE;
  if (bsp != NULL && ISA_na (bsp->mol)) {
    if (mip != NULL) {
      if (mip->biomol == MOLECULE_TYPE_GENOMIC && mip->completeness == 1) {
        sev = SEV_ERROR;
        if (mip->tech == MI_TECH_htgs_3) {
          sev = SEV_WARNING;
        }
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_GENBANK) {
            is_gb = TRUE;
          }
        }
        if (is_gb) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
          if (vnp != NULL) {
            str = (CharPtr) vnp->data.ptrvalue;
            if (! StringHasNoText (str)) {
              if (StringISearch (str, "complete sequence") == NULL &&
                  StringISearch (str, "complete genome") == NULL) {

                oldEntityID = gcp->entityID;
                oldItemID = gcp->itemID;
                oldItemtype = gcp->thistype;

                gcp->entityID = mipEntityID;
                gcp->itemID = mipItemID;
                gcp->thistype = mipItemtype;

                if (bsp->topology == 2) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_INST_CompleteCircleProblem, "Circular topology has complete flag set, but title should say complete sequence or complete genome");
                } else {
                  ValidErr (vsp, sev, ERR_SEQ_DESCR_UnwantedCompleteFlag, "Suspicious use of complete");
                }

                gcp->entityID = oldEntityID;
                gcp->itemID = oldItemID;
                gcp->thistype = oldItemtype;
              }
            }
          }
        }
      } else if (mip->biomol == MOLECULE_TYPE_MRNA) {
        bvs.is_mrna = TRUE;
      } else if (mip->biomol == MOLECULE_TYPE_PRE_MRNA) {
        bvs.is_prerna = TRUE;
      }
      if (mip->biomol >= MOLECULE_TYPE_PRE_MRNA && mip->biomol <= MOLECULE_TYPE_SCRNA && bsp->mol == Seq_mol_dna) {
        /* - this is how we indicate an mRNA sequenced from a cDNA, so no error
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_RnaDnaConflict, "MolInfo says RNA, Bioseq says DNA");
        */
      }
    } else if (bsp->mol == Seq_mol_rna) {
      bvs.is_mrna = TRUE;       /* if no molinfo, assume rna is mrna */
    }
  }

  if (mip != NULL) {
    if (mip->tech == MI_TECH_sts ||
        mip->tech == MI_TECH_survey ||
        mip->tech == MI_TECH_wgs ||
        mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 ||
        mip->tech == MI_TECH_htgs_2 || mip->tech == MI_TECH_htgs_3) {
      if (mip->tech == MI_TECH_sts && bsp->mol == Seq_mol_rna && mip->biomol == MOLECULE_TYPE_MRNA) {
        /* there are some STS sequences derived from cDNAs, so do not report these */
      } else if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingBiomolTech, "HTGS/STS/GSS/WGS sequence should be genomic");
      } else if (bsp == NULL || (bsp->mol != Seq_mol_dna && bsp->mol != Seq_mol_na)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingBiomolTech, "HTGS/STS/GSS/WGS sequence should not be RNA");
      }
    }
  }

  if (ISA_aa (bsp->mol)) {
    bvs.is_aa = TRUE;
    /* check proteins in nuc-prot set have a CdRegion */
    if (vsp->bssp != NULL) {
      if (vsp->bssp->_class == 1) {     /* in a nuc-prot set */
        if (vsp->useSeqMgrIndexes) {
          sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (sfp == NULL) {
            sfp = SeqMgrGetPROTgivenProduct (bsp, NULL); /* now instantiating and indexing products of protein processing */
          }
        } else {
          sfp = SeqEntryGetSeqFeat (vsp->sep, 3, NULL, NULL, 1, bsp);
        }
        if (sfp == NULL)        /* no CdRegion points to this bsp */
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_PKG_NoCdRegionPtr, "No CdRegion in nuc-prot set points to this protein");
      }
    }
  }

  if (vsp->useSeqMgrIndexes) {
    bvs.gcp = gcp;
    bvs.bsp = bsp;
    ValidateBioseqContextIndexed (bsp, &bvs);
  } else {
    GatherSeqEntry (vsp->sep, &bvs, ValidateBioseqContextGather, &gs);
  }

  vsp->gcp = gcp;               /* reset the gcp pointer changed in previous gather */
  vsp->descr = NULL;
  vsp->sfp = NULL;

  if ((!bvs.got_a_pub) && (!vsp->suppress_no_pubs) && (! vsp->seqSubmitParent)) {
    omdp = NULL;
    if (gcp != NULL) {
      omdp = ObjMgrGetData (gcp->entityID);
    }
    if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
      sev = SEV_ERROR;
      if (!IsNoncuratedRefSeq (bsp, &sev)) {
        if (! IsWgsIntermediate (vsp->sep)) {
          ValidErr (vsp, sev, ERR_SEQ_DESCR_NoPubFound, "No publications refer to this Bioseq.");
        }
      }
    }
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      is_local = TRUE;
    } else {
      is_local_only = FALSE;
    }
    if (sip->choice == SEQID_PATENT) {
      psip = (PatentSeqIdPtr) sip->data.ptrvalue;
      if (psip != NULL) {
        ipp = psip->cit;
        if (ipp != NULL && StringICmp (ipp->country, "US") == 0)
          return;
      }
      return;
    } else if (sip->choice == SEQID_PDB) {
      isPDB = TRUE;
    } else if (sip->choice == SEQID_GENBANK ||
               sip->choice == SEQID_EMBL ||
               sip->choice == SEQID_DDBJ) {
      is_gb = TRUE;
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        acclen = StringLen (tsip->accession);
        if (acclen == 12) {
          is_wgs = TRUE;
        } else if (acclen == 13) {
          is_wgs = TRUE;
        /*
        } else if (StringNCmp (tsip->accession, "CH", 2) == 0 ||
                   StringNCmp (tsip->accession, "CM", 2) == 0) {
          is_ch_or_cm = TRUE;
        */
        } else if (WHICH_db_accession (tsip->accession) == ACCN_NCBI_SEGSET) {
          is_ch_or_cm = TRUE;
        }
      }
    } else if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
            StringNCmp (tsip->accession, "NP_", 3) == 0 ||
            StringNCmp (tsip->accession, "NG_", 3) == 0 ||
            StringNCmp (tsip->accession, "NR_", 3) == 0) {
          is_gb = TRUE;
        } else if (StringNCmp (tsip->accession, "NC_", 3) == 0) {
          is_nc = TRUE;
        }
      }
    }
  }
  if (! is_local) {
    is_local_only = FALSE;
  }
  if (is_wgs) {
    if (mip == NULL || mip->tech != MI_TECH_wgs) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "WGS accession should have Mol-info.tech of wgs");
    }
  } else if (mip != NULL && mip->tech == MI_TECH_wgs && is_gb) {
    if (is_ch_or_cm || is_local_only) {
      /* skip warning if CH or CM or SEQID_LOCAL only */
    } else {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Mol-info.tech of wgs should have WGS accession");
    }
  }
  if (is_nc) {
    if (mip != NULL && mip->biomol != MOLECULE_TYPE_GENOMIC && mip->biomol != MOLECULE_TYPE_CRNA) {
      if (ISA_na (bsp->mol)) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "NC nucldotide should be genomic or cRNA");
      }
    }
  }

  if ((!bvs.last_org) && (!vsp->suppress_no_biosrc))
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name has been applied to this Bioseq.  Other qualifiers may exist.");


  if ((bvs.is_aa) && (bvs.num_full_length_prot_ref == 0) && (!isPDB) && (bsp->repr != Seq_repr_virtual)) {
    parent = SeqMgrGetParentOfPart (bsp, NULL);
    if (parent == NULL || SeqMgrGetBestProteinFeature (bsp, NULL) == NULL) {

      oldEntityID = gcp->entityID;
      oldItemID = gcp->itemID;
      oldItemtype = gcp->thistype;

      if (SeqMgrGetCDSgivenProduct (bsp, &fcontext) != NULL) {
        gcp->entityID = fcontext.entityID;
        gcp->itemID = fcontext.itemID;
        gcp->thistype = OBJ_SEQFEAT;
      }

      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NoProtRefFound, "No full length Prot-ref feature applied to this Bioseq");

      gcp->entityID = oldEntityID;
      gcp->itemID = oldItemID;
      gcp->thistype = oldItemtype;
    }
  } else if (bvs.is_aa && bvs.num_full_length_prot_ref > 1) {
    if (bvs.num_justprot > 1 ||
        bvs.num_preprot > 1 ||
        bvs.num_matpep > 1 ||
        bvs.num_sigpep > 1 ||
        bvs.num_transpep > 1) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MultipleProtRefs, "%d full-length protein features present on protein",
                (int) bvs.num_full_length_prot_ref);
    }
  }

  /* now flag missing molinfo even if not in Sequin */
  if (mip == NULL && (!isPDB)) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_NoMolInfoFound, "No Mol-info applies to this Bioseq");
  }

#if 0 /* temporarily suppress */
  /* if tech is TSA, must have assembly */
  if (mip != NULL && mip->tech == MI_TECH_tsa
      && (bsp->hist == NULL || bsp->hist->assembly == NULL)) {
    SeqIdWrite (bsp->id, buf1, PRINTID_FASTA_SHORT, sizeof (buf1) - 1);
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_INST_TSAHistAssemblyMissing, "TSA record %s should have Seq-hist.assembly", buf1);    
  }
#endif

  /* look for genes that overlap two other genes */
  FindMultiGeneOverlaps (bsp, vsp);

  /* TSA checks */
  ValidateTSASequenceForNs (bsp, vsp);
}

/*****************************************************************************
*
*   ValidateSeqFeat(gcp)
*
*****************************************************************************/
static Boolean EmptyOrNullString (CharPtr str)
{
  Char            ch;

  if (str == NULL)
    return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~')
      return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

static void CheckPeptideOnCodonBoundary (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, CharPtr key)
{
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  SeqLocPtr       first = NULL, last = NULL, slp = NULL;
  Boolean         partial5, partial3;
  Int4            pos1, pos2, adjust = 0, mod1, mod2;

  cds = SeqMgrGetOverlappingCDS (sfp->location, NULL);
  if (cds == NULL)
    return;
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  if (crp == NULL)
    return;
  if (crp->frame == 2) {
    adjust = 1;
  } else if (crp->frame == 3) {
    adjust = 2;
  }

  while ((slp = SeqLocFindNext (sfp->location, slp)) != NULL) {
    last = slp;
    if (first == NULL) {
      first = slp;
    }
  }
  if (first == NULL || last == NULL)
    return;

  pos1 = GetOffsetInLoc (first, cds->location, SEQLOC_START) - adjust;
  pos2 = GetOffsetInLoc (last, cds->location, SEQLOC_STOP) - adjust;
  mod1 = pos1 % 3;
  mod2 = pos2 % 3;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5) {
    mod1 = 0;
  }
  if (partial3) {
    mod2 = 2;
  }

  if (mod1 != 0 && mod2 != 2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Start and stop of %s are out of frame with CDS codons", key);
  } else if (mod1 != 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Start of %s is out of frame with CDS codons", key);
  } else if (mod2 != 2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PeptideFeatOutOfFrame, "Stop of %s is out of frame with CDS codons", key);
  }
}

static CharPtr  legal_repeat_types[] = {
  "tandem", "inverted", "flanking", "terminal",
  "direct", "dispersed", "other", NULL
};

static CharPtr legal_cons_splice_strings [] = {
  "(5'site:YES, 3'site:YES)",
  "(5'site:YES, 3'site:NO)",
  "(5'site:YES, 3'site:ABSENT)",
  "(5'site:NO, 3'site:YES)",
  "(5'site:NO, 3'site:NO)",
  "(5'site:NO, 3'site:ABSENT)",
  "(5'site:ABSENT, 3'site:YES)",
  "(5'site:ABSENT, 3'site:NO)",
  "(5'site:ABSENT, 3'site:ABSENT)",
  NULL
};

static CharPtr legal_mobile_element_strings [] = {
  "transposon",
  "retrotransposon",
  "integron",
  "insertion sequence",
  "non-LTR retrotransposon",
  "SINE",
  "MITE",
  "LINE",
  "other",
  NULL
};

static CharPtr ecnum_ambig [] = {
  "1.-.-.-", "1.1.-.-", "1.1.1.-", "1.1.1.n", "1.1.2.-",
  "1.1.2.n", "1.1.3.-", "1.1.3.n", "1.1.4.-", "1.1.4.n",
  "1.1.5.-", "1.1.5.n", "1.1.99.-", "1.1.99.n", "1.1.n.n",
  "1.2.-.-", "1.2.1.-", "1.2.1.n", "1.2.2.-", "1.2.2.n",
  "1.2.3.-", "1.2.3.n", "1.2.4.-", "1.2.4.n", "1.2.7.-",
  "1.2.7.n", "1.2.99.-", "1.2.99.n", "1.2.n.n", "1.3.-.-",
  "1.3.1.-", "1.3.1.n", "1.3.2.-", "1.3.2.n", "1.3.3.-",
  "1.3.3.n", "1.3.5.-", "1.3.5.n", "1.3.7.-", "1.3.7.n",
  "1.3.99.-", "1.3.99.n", "1.3.n.n", "1.4.-.-", "1.4.1.-",
  "1.4.1.n", "1.4.2.-", "1.4.2.n", "1.4.3.-", "1.4.3.n",
  "1.4.4.-", "1.4.4.n", "1.4.7.-", "1.4.7.n", "1.4.99.-",
  "1.4.99.n", "1.4.n.n", "1.5.-.-", "1.5.1.-", "1.5.1.n",
  "1.5.3.-", "1.5.3.n", "1.5.4.-", "1.5.4.n", "1.5.5.-",
  "1.5.5.n", "1.5.7.-", "1.5.7.n", "1.5.8.-", "1.5.8.n",
  "1.5.99.-", "1.5.99.n", "1.5.n.n", "1.6.-.-", "1.6.1.-",
  "1.6.1.n", "1.6.2.-", "1.6.2.n", "1.6.3.-", "1.6.3.n",
  "1.6.5.-", "1.6.5.n", "1.6.6.-", "1.6.6.n", "1.6.99.-",
  "1.6.99.n", "1.6.n.n", "1.7.-.-", "1.7.1.-", "1.7.1.n",
  "1.7.2.-", "1.7.2.n", "1.7.3.-", "1.7.3.n", "1.7.7.-",
  "1.7.7.n", "1.7.99.-", "1.7.99.n", "1.7.n.n", "1.8.-.-",
  "1.8.1.-", "1.8.1.n", "1.8.2.-", "1.8.2.n", "1.8.3.-",
  "1.8.3.n", "1.8.4.-", "1.8.4.n", "1.8.5.-", "1.8.5.n",
  "1.8.7.-", "1.8.7.n", "1.8.98.-", "1.8.98.n", "1.8.99.-",
  "1.8.99.n", "1.8.n.n", "1.9.-.-", "1.9.3.-", "1.9.3.n",
  "1.9.6.-", "1.9.6.n", "1.9.99.-", "1.9.99.n", "1.9.n.n",
  "1.10.-.-", "1.10.1.-", "1.10.1.n", "1.10.2.-", "1.10.2.n",
  "1.10.3.-", "1.10.3.n", "1.10.99.-", "1.10.99.n",
  "1.10.n.n", "1.11.-.-", "1.11.1.-", "1.11.1.n", "1.11.n.n",
  "1.12.-.-", "1.12.1.-", "1.12.1.n", "1.12.2.-", "1.12.2.n",
  "1.12.5.-", "1.12.5.n", "1.12.7.-", "1.12.7.n", "1.12.98.-",
  "1.12.98.n", "1.12.99.-", "1.12.99.n", "1.12.n.n",
  "1.13.-.-", "1.13.11.-", "1.13.11.n", "1.13.12.-",
  "1.13.12.n", "1.13.99.-", "1.13.99.n", "1.13.n.n",
  "1.14.-.-", "1.14.11.-", "1.14.11.n", "1.14.12.-",
  "1.14.12.n", "1.14.13.-", "1.14.13.n", "1.14.14.-",
  "1.14.14.n", "1.14.15.-", "1.14.15.n", "1.14.16.-",
  "1.14.16.n", "1.14.17.-", "1.14.17.n", "1.14.18.-",
  "1.14.18.n", "1.14.19.-", "1.14.19.n", "1.14.20.-",
  "1.14.20.n", "1.14.21.-", "1.14.21.n", "1.14.99.-",
  "1.14.99.n", "1.14.n.n", "1.15.-.-", "1.15.1.-", "1.15.1.n",
  "1.15.n.n", "1.16.-.-", "1.16.1.-", "1.16.1.n", "1.16.3.-",
  "1.16.3.n", "1.16.8.-", "1.16.8.n", "1.16.n.n", "1.17.-.-",
  "1.17.1.-", "1.17.1.n", "1.17.3.-", "1.17.3.n", "1.17.4.-",
  "1.17.4.n", "1.17.5.-", "1.17.5.n", "1.17.99.-",
  "1.17.99.n", "1.17.n.n", "1.18.-.-", "1.18.1.-", "1.18.1.n",
  "1.18.6.-", "1.18.6.n", "1.18.n.n", "1.19.-.-", "1.19.6.-",
  "1.19.6.n", "1.19.n.n", "1.20.-.-", "1.20.1.-", "1.20.1.n",
  "1.20.4.-", "1.20.4.n", "1.20.98.-", "1.20.98.n",
  "1.20.99.-", "1.20.99.n", "1.20.n.n", "1.21.-.-",
  "1.21.3.-", "1.21.3.n", "1.21.4.-", "1.21.4.n", "1.21.99.-",
  "1.21.99.n", "1.21.n.n", "1.97.-.-", "1.97.1.-", "1.97.1.n",
  "1.97.n.n", "1.n.n.n", "2.-.-.-", "2.1.-.-", "2.1.1.-",
  "2.1.1.n", "2.1.2.-", "2.1.2.n", "2.1.3.-", "2.1.3.n",
  "2.1.4.-", "2.1.4.n", "2.1.n.n", "2.2.-.-", "2.2.1.-",
  "2.2.1.n", "2.2.n.n", "2.3.-.-", "2.3.1.-", "2.3.1.n",
  "2.3.2.-", "2.3.2.n", "2.3.3.-", "2.3.3.n", "2.3.n.n",
  "2.4.-.-", "2.4.1.-", "2.4.1.n", "2.4.2.-", "2.4.2.n",
  "2.4.99.-", "2.4.99.n", "2.4.n.n", "2.5.-.-", "2.5.1.-",
  "2.5.1.n", "2.5.n.n", "2.6.-.-", "2.6.1.-", "2.6.1.n",
  "2.6.3.-", "2.6.3.n", "2.6.99.-", "2.6.99.n", "2.6.n.n",
  "2.7.-.-", "2.7.1.-", "2.7.1.n", "2.7.2.-", "2.7.2.n",
  "2.7.3.-", "2.7.3.n", "2.7.4.-", "2.7.4.n", "2.7.6.-",
  "2.7.6.n", "2.7.7.-", "2.7.7.n", "2.7.8.-", "2.7.8.n",
  "2.7.9.-", "2.7.9.n", "2.7.10.-", "2.7.10.n", "2.7.11.-",
  "2.7.11.n", "2.7.12.-", "2.7.12.n", "2.7.13.-", "2.7.13.n",
  "2.7.99.-", "2.7.99.n", "2.7.n.n", "2.8.-.-", "2.8.1.-",
  "2.8.1.n", "2.8.2.-", "2.8.2.n", "2.8.3.-", "2.8.3.n",
  "2.8.4.-", "2.8.4.n", "2.8.n.n", "2.9.-.-", "2.9.1.-",
  "2.9.1.n", "2.9.n.n", "2.n.n.n", "3.-.-.-", "3.1.-.-",
  "3.1.1.-", "3.1.1.n", "3.1.2.-", "3.1.2.n", "3.1.3.-",
  "3.1.3.n", "3.1.4.-", "3.1.4.n", "3.1.5.-", "3.1.5.n",
  "3.1.6.-", "3.1.6.n", "3.1.7.-", "3.1.7.n", "3.1.8.-",
  "3.1.8.n", "3.1.11.-", "3.1.11.n", "3.1.13.-", "3.1.13.n",
  "3.1.14.-", "3.1.14.n", "3.1.15.-", "3.1.15.n", "3.1.16.-",
  "3.1.16.n", "3.1.21.-", "3.1.21.n", "3.1.22.-", "3.1.22.n",
  "3.1.25.-", "3.1.25.n", "3.1.26.-", "3.1.26.n", "3.1.27.-",
  "3.1.27.n", "3.1.30.-", "3.1.30.n", "3.1.31.-", "3.1.31.n",
  "3.1.n.n", "3.2.-.-", "3.2.1.-", "3.2.1.n", "3.2.2.-",
  "3.2.2.n", "3.2.n.n", "3.3.-.-", "3.3.1.-", "3.3.1.n",
  "3.3.2.-", "3.3.2.n", "3.3.n.n", "3.4.-.-", "3.4.11.-",
  "3.4.11.n", "3.4.13.-", "3.4.13.n", "3.4.14.-", "3.4.14.n",
  "3.4.15.-", "3.4.15.n", "3.4.16.-", "3.4.16.n", "3.4.17.-",
  "3.4.17.n", "3.4.18.-", "3.4.18.n", "3.4.19.-", "3.4.19.n",
  "3.4.21.-", "3.4.21.n", "3.4.22.-", "3.4.22.n", "3.4.23.-",
  "3.4.23.n", "3.4.24.-", "3.4.24.n", "3.4.25.-", "3.4.25.n",
  "3.4.n.n", "3.5.-.-", "3.5.1.-", "3.5.1.n", "3.5.2.-",
  "3.5.2.n", "3.5.3.-", "3.5.3.n", "3.5.4.-", "3.5.4.n",
  "3.5.5.-", "3.5.5.n", "3.5.99.-", "3.5.99.n", "3.5.n.n",
  "3.6.-.-", "3.6.1.-", "3.6.1.n", "3.6.2.-", "3.6.2.n",
  "3.6.3.-", "3.6.3.n", "3.6.4.-", "3.6.4.n", "3.6.5.-",
  "3.6.5.n", "3.6.n.n", "3.7.-.-", "3.7.1.-", "3.7.1.n",
  "3.7.n.n", "3.8.-.-", "3.8.1.-", "3.8.1.n", "3.8.n.n",
  "3.9.-.-", "3.9.1.-", "3.9.1.n", "3.9.n.n", "3.10.-.-",
  "3.10.1.-", "3.10.1.n", "3.10.n.n", "3.11.-.-", "3.11.1.-",
  "3.11.1.n", "3.11.n.n", "3.12.-.-", "3.12.1.-", "3.12.1.n",
  "3.12.n.n", "3.13.-.-", "3.13.1.-", "3.13.1.n", "3.13.n.n",
  "3.n.n.n", "4.-.-.-", "4.1.-.-", "4.1.1.-", "4.1.1.n",
  "4.1.2.-", "4.1.2.n", "4.1.3.-", "4.1.3.n", "4.1.99.-",
  "4.1.99.n", "4.1.n.n", "4.2.-.-", "4.2.1.-", "4.2.1.n",
  "4.2.2.-", "4.2.2.n", "4.2.3.-", "4.2.3.n", "4.2.99.-",
  "4.2.99.n", "4.2.n.n", "4.3.-.-", "4.3.1.-", "4.3.1.n",
  "4.3.2.-", "4.3.2.n", "4.3.3.-", "4.3.3.n", "4.3.n.n",
  "4.4.-.-", "4.4.1.-", "4.4.1.n", "4.4.n.n", "4.5.-.-",
  "4.5.1.-", "4.5.1.n", "4.5.n.n", "4.6.-.-", "4.6.1.-",
  "4.6.1.n", "4.6.n.n", "4.99.-.-", "4.99.1.-", "4.99.1.n",
  "4.99.n.n", "4.n.n.n", "5.-.-.-", "5.1.-.-", "5.1.1.-",
  "5.1.1.n", "5.1.2.-", "5.1.2.n", "5.1.3.-", "5.1.3.n",
  "5.1.99.-", "5.1.99.n", "5.1.n.n", "5.2.-.-", "5.2.1.-",
  "5.2.1.n", "5.2.n.n", "5.3.-.-", "5.3.1.-", "5.3.1.n",
  "5.3.2.-", "5.3.2.n", "5.3.3.-", "5.3.3.n", "5.3.4.-",
  "5.3.4.n", "5.3.99.-", "5.3.99.n", "5.3.n.n", "5.4.-.-",
  "5.4.1.-", "5.4.1.n", "5.4.2.-", "5.4.2.n", "5.4.3.-",
  "5.4.3.n", "5.4.4.-", "5.4.4.n", "5.4.99.-", "5.4.99.n",
  "5.4.n.n", "5.5.-.-", "5.5.1.-", "5.5.1.n", "5.5.n.n",
  "5.99.-.-", "5.99.1.-", "5.99.1.n", "5.99.n.n", "5.n.n.n",
  "6.-.-.-", "6.1.-.-", "6.1.1.-", "6.1.1.n", "6.1.n.n",
  "6.2.-.-", "6.2.1.-", "6.2.1.n", "6.2.n.n", "6.3.-.-",
  "6.3.1.-", "6.3.1.n", "6.3.2.-", "6.3.2.n", "6.3.3.-",
  "6.3.3.n", "6.3.4.-", "6.3.4.n", "6.3.5.-", "6.3.5.n",
  "6.3.n.n", "6.4.-.-", "6.4.1.-", "6.4.1.n", "6.4.n.n",
  "6.5.-.-", "6.5.1.-", "6.5.1.n", "6.5.n.n", "6.6.-.-",
  "6.6.1.-", "6.6.1.n", "6.6.n.n", "6.n.n.n",
  NULL
};

static CharPtr ecnum_specif [] = {
  "1.1.1.1", "1.1.1.2", "1.1.1.3", "1.1.1.4", "1.1.1.5",
  "1.1.1.6", "1.1.1.7", "1.1.1.8", "1.1.1.9", "1.1.1.10",
  "1.1.1.11", "1.1.1.12", "1.1.1.13", "1.1.1.14", "1.1.1.15",
  "1.1.1.16", "1.1.1.17", "1.1.1.18", "1.1.1.19", "1.1.1.20",
  "1.1.1.21", "1.1.1.22", "1.1.1.23", "1.1.1.24", "1.1.1.25",
  "1.1.1.26", "1.1.1.27", "1.1.1.28", "1.1.1.29", "1.1.1.30",
  "1.1.1.31", "1.1.1.32", "1.1.1.33", "1.1.1.34", "1.1.1.35",
  "1.1.1.36", "1.1.1.37", "1.1.1.38", "1.1.1.39", "1.1.1.40",
  "1.1.1.41", "1.1.1.42", "1.1.1.43", "1.1.1.44", "1.1.1.45",
  "1.1.1.46", "1.1.1.47", "1.1.1.48", "1.1.1.49", "1.1.1.50",
  "1.1.1.51", "1.1.1.52", "1.1.1.53", "1.1.1.54", "1.1.1.55",
  "1.1.1.56", "1.1.1.57", "1.1.1.58", "1.1.1.59", "1.1.1.60",
  "1.1.1.61", "1.1.1.62", "1.1.1.63", "1.1.1.64", "1.1.1.65",
  "1.1.1.66", "1.1.1.67", "1.1.1.69", "1.1.1.71", "1.1.1.72",
  "1.1.1.73", "1.1.1.75", "1.1.1.76", "1.1.1.77", "1.1.1.78",
  "1.1.1.79", "1.1.1.80", "1.1.1.81", "1.1.1.82", "1.1.1.83",
  "1.1.1.84", "1.1.1.85", "1.1.1.86", "1.1.1.87", "1.1.1.88",
  "1.1.1.90", "1.1.1.91", "1.1.1.92", "1.1.1.93", "1.1.1.94",
  "1.1.1.95", "1.1.1.96", "1.1.1.97", "1.1.1.98", "1.1.1.99",
  "1.1.1.100", "1.1.1.101", "1.1.1.102", "1.1.1.103",
  "1.1.1.104", "1.1.1.105", "1.1.1.106", "1.1.1.107",
  "1.1.1.108", "1.1.1.110", "1.1.1.111", "1.1.1.112",
  "1.1.1.113", "1.1.1.114", "1.1.1.115", "1.1.1.116",
  "1.1.1.117", "1.1.1.118", "1.1.1.119", "1.1.1.120",
  "1.1.1.121", "1.1.1.122", "1.1.1.123", "1.1.1.124",
  "1.1.1.125", "1.1.1.126", "1.1.1.127", "1.1.1.128",
  "1.1.1.129", "1.1.1.130", "1.1.1.131", "1.1.1.132",
  "1.1.1.133", "1.1.1.134", "1.1.1.135", "1.1.1.136",
  "1.1.1.137", "1.1.1.138", "1.1.1.140", "1.1.1.141",
  "1.1.1.142", "1.1.1.143", "1.1.1.144", "1.1.1.145",
  "1.1.1.146", "1.1.1.147", "1.1.1.148", "1.1.1.149",
  "1.1.1.150", "1.1.1.151", "1.1.1.152", "1.1.1.153",
  "1.1.1.154", "1.1.1.156", "1.1.1.157", "1.1.1.158",
  "1.1.1.159", "1.1.1.160", "1.1.1.161", "1.1.1.162",
  "1.1.1.163", "1.1.1.164", "1.1.1.165", "1.1.1.166",
  "1.1.1.167", "1.1.1.168", "1.1.1.169", "1.1.1.170",
  "1.1.1.172", "1.1.1.173", "1.1.1.174", "1.1.1.175",
  "1.1.1.176", "1.1.1.177", "1.1.1.178", "1.1.1.179",
  "1.1.1.181", "1.1.1.183", "1.1.1.184", "1.1.1.185",
  "1.1.1.186", "1.1.1.187", "1.1.1.188", "1.1.1.189",
  "1.1.1.190", "1.1.1.191", "1.1.1.192", "1.1.1.193",
  "1.1.1.194", "1.1.1.195", "1.1.1.196", "1.1.1.197",
  "1.1.1.198", "1.1.1.199", "1.1.1.200", "1.1.1.201",
  "1.1.1.202", "1.1.1.203", "1.1.1.205", "1.1.1.206",
  "1.1.1.207", "1.1.1.208", "1.1.1.209", "1.1.1.210",
  "1.1.1.211", "1.1.1.212", "1.1.1.213", "1.1.1.214",
  "1.1.1.215", "1.1.1.216", "1.1.1.217", "1.1.1.218",
  "1.1.1.219", "1.1.1.220", "1.1.1.221", "1.1.1.222",
  "1.1.1.223", "1.1.1.224", "1.1.1.225", "1.1.1.226",
  "1.1.1.227", "1.1.1.228", "1.1.1.229", "1.1.1.230",
  "1.1.1.231", "1.1.1.232", "1.1.1.233", "1.1.1.234",
  "1.1.1.235", "1.1.1.236", "1.1.1.237", "1.1.1.238",
  "1.1.1.239", "1.1.1.240", "1.1.1.241", "1.1.1.243",
  "1.1.1.244", "1.1.1.245", "1.1.1.246", "1.1.1.247",
  "1.1.1.248", "1.1.1.250", "1.1.1.251", "1.1.1.252",
  "1.1.1.254", "1.1.1.255", "1.1.1.256", "1.1.1.257",
  "1.1.1.258", "1.1.1.259", "1.1.1.260", "1.1.1.261",
  "1.1.1.262", "1.1.1.263", "1.1.1.264", "1.1.1.265",
  "1.1.1.266", "1.1.1.267", "1.1.1.268", "1.1.1.269",
  "1.1.1.270", "1.1.1.271", "1.1.1.272", "1.1.1.273",
  "1.1.1.274", "1.1.1.275", "1.1.1.276", "1.1.1.277",
  "1.1.1.278", "1.1.1.279", "1.1.1.280", "1.1.1.281",
  "1.1.1.282", "1.1.1.283", "1.1.1.284", "1.1.1.285",
  "1.1.1.286", "1.1.1.287", "1.1.1.288", "1.1.1.289",
  "1.1.1.290", "1.1.1.291", "1.1.1.292", "1.1.1.294",
  "1.1.1.295", "1.1.1.296", "1.1.1.297", "1.1.2.2", "1.1.2.3",
  "1.1.2.4", "1.1.2.5", "1.1.3.3", "1.1.3.4", "1.1.3.5",
  "1.1.3.6", "1.1.3.7", "1.1.3.8", "1.1.3.9", "1.1.3.10",
  "1.1.3.11", "1.1.3.12", "1.1.3.13", "1.1.3.14", "1.1.3.15",
  "1.1.3.16", "1.1.3.17", "1.1.3.18", "1.1.3.19", "1.1.3.20",
  "1.1.3.21", "1.1.3.23", "1.1.3.27", "1.1.3.28", "1.1.3.29",
  "1.1.3.30", "1.1.3.37", "1.1.3.38", "1.1.3.39", "1.1.3.40",
  "1.1.3.41", "1.1.4.1", "1.1.4.2", "1.1.5.2", "1.1.99.1",
  "1.1.99.2", "1.1.99.3", "1.1.99.4", "1.1.99.5", "1.1.99.6",
  "1.1.99.7", "1.1.99.8", "1.1.99.9", "1.1.99.10",
  "1.1.99.11", "1.1.99.12", "1.1.99.13", "1.1.99.14",
  "1.1.99.16", "1.1.99.18", "1.1.99.20", "1.1.99.21",
  "1.1.99.22", "1.1.99.23", "1.1.99.24", "1.1.99.25",
  "1.1.99.26", "1.1.99.27", "1.1.99.28", "1.1.99.29",
  "1.1.99.30", "1.1.99.31", "1.1.99.32", "1.2.1.2", "1.2.1.3",
  "1.2.1.4", "1.2.1.5", "1.2.1.7", "1.2.1.8", "1.2.1.9",
  "1.2.1.10", "1.2.1.11", "1.2.1.12", "1.2.1.13", "1.2.1.15",
  "1.2.1.16", "1.2.1.17", "1.2.1.18", "1.2.1.19", "1.2.1.20",
  "1.2.1.21", "1.2.1.22", "1.2.1.23", "1.2.1.24", "1.2.1.25",
  "1.2.1.26", "1.2.1.27", "1.2.1.28", "1.2.1.29", "1.2.1.30",
  "1.2.1.31", "1.2.1.32", "1.2.1.33", "1.2.1.36", "1.2.1.38",
  "1.2.1.39", "1.2.1.40", "1.2.1.41", "1.2.1.42", "1.2.1.43",
  "1.2.1.44", "1.2.1.45", "1.2.1.46", "1.2.1.47", "1.2.1.48",
  "1.2.1.49", "1.2.1.50", "1.2.1.51", "1.2.1.52", "1.2.1.53",
  "1.2.1.54", "1.2.1.57", "1.2.1.58", "1.2.1.59", "1.2.1.60",
  "1.2.1.61", "1.2.1.62", "1.2.1.63", "1.2.1.64", "1.2.1.65",
  "1.2.1.66", "1.2.1.67", "1.2.1.68", "1.2.1.69", "1.2.1.70",
  "1.2.1.71", "1.2.1.72", "1.2.1.73", "1.2.2.1", "1.2.2.2",
  "1.2.2.3", "1.2.2.4", "1.2.3.1", "1.2.3.3", "1.2.3.4",
  "1.2.3.5", "1.2.3.6", "1.2.3.7", "1.2.3.8", "1.2.3.9",
  "1.2.3.11", "1.2.3.13", "1.2.3.14", "1.2.4.1", "1.2.4.2",
  "1.2.4.4", "1.2.7.1", "1.2.7.2", "1.2.7.3", "1.2.7.4",
  "1.2.7.5", "1.2.7.6", "1.2.7.7", "1.2.7.8", "1.2.99.2",
  "1.2.99.3", "1.2.99.4", "1.2.99.5", "1.2.99.6", "1.2.99.7",
  "1.3.1.1", "1.3.1.2", "1.3.1.3", "1.3.1.4", "1.3.1.5",
  "1.3.1.6", "1.3.1.7", "1.3.1.8", "1.3.1.9", "1.3.1.10",
  "1.3.1.11", "1.3.1.12", "1.3.1.13", "1.3.1.14", "1.3.1.15",
  "1.3.1.16", "1.3.1.17", "1.3.1.18", "1.3.1.19", "1.3.1.20",
  "1.3.1.21", "1.3.1.22", "1.3.1.24", "1.3.1.25", "1.3.1.26",
  "1.3.1.27", "1.3.1.28", "1.3.1.29", "1.3.1.30", "1.3.1.31",
  "1.3.1.32", "1.3.1.33", "1.3.1.34", "1.3.1.35", "1.3.1.36",
  "1.3.1.37", "1.3.1.38", "1.3.1.39", "1.3.1.40", "1.3.1.41",
  "1.3.1.42", "1.3.1.43", "1.3.1.44", "1.3.1.45", "1.3.1.46",
  "1.3.1.47", "1.3.1.48", "1.3.1.49", "1.3.1.51", "1.3.1.52",
  "1.3.1.53", "1.3.1.54", "1.3.1.56", "1.3.1.57", "1.3.1.58",
  "1.3.1.60", "1.3.1.62", "1.3.1.63", "1.3.1.64", "1.3.1.65",
  "1.3.1.66", "1.3.1.67", "1.3.1.68", "1.3.1.69", "1.3.1.70",
  "1.3.1.71", "1.3.1.72", "1.3.1.73", "1.3.1.74", "1.3.1.75",
  "1.3.1.76", "1.3.1.77", "1.3.1.78", "1.3.1.79", "1.3.1.80",
  "1.3.1.81", "1.3.1.82", "1.3.2.3", "1.3.3.1", "1.3.3.3",
  "1.3.3.4", "1.3.3.5", "1.3.3.6", "1.3.3.7", "1.3.3.8",
  "1.3.3.9", "1.3.3.10", "1.3.3.11", "1.3.3.12", "1.3.5.1",
  "1.3.7.1", "1.3.7.2", "1.3.7.3", "1.3.7.4", "1.3.7.5",
  "1.3.7.6", "1.3.99.1", "1.3.99.2", "1.3.99.3", "1.3.99.4",
  "1.3.99.5", "1.3.99.6", "1.3.99.7", "1.3.99.8", "1.3.99.10",
  "1.3.99.11", "1.3.99.12", "1.3.99.13", "1.3.99.14",
  "1.3.99.15", "1.3.99.16", "1.3.99.17", "1.3.99.18",
  "1.3.99.19", "1.3.99.20", "1.3.99.21", "1.3.99.22",
  "1.3.99.23", "1.3.99.24", "1.3.99.25", "1.4.1.1", "1.4.1.2",
  "1.4.1.3", "1.4.1.4", "1.4.1.5", "1.4.1.7", "1.4.1.8",
  "1.4.1.9", "1.4.1.10", "1.4.1.11", "1.4.1.12", "1.4.1.13",
  "1.4.1.14", "1.4.1.15", "1.4.1.16", "1.4.1.17", "1.4.1.18",
  "1.4.1.19", "1.4.1.20", "1.4.1.21", "1.4.2.1", "1.4.3.1",
  "1.4.3.2", "1.4.3.3", "1.4.3.4", "1.4.3.5", "1.4.3.7",
  "1.4.3.8", "1.4.3.10", "1.4.3.11", "1.4.3.12", "1.4.3.13",
  "1.4.3.14", "1.4.3.15", "1.4.3.16", "1.4.3.19", "1.4.3.20",
  "1.4.3.21", "1.4.3.22", "1.4.4.2", "1.4.7.1", "1.4.99.1",
  "1.4.99.2", "1.4.99.3", "1.4.99.4", "1.4.99.5", "1.5.1.1",
  "1.5.1.2", "1.5.1.3", "1.5.1.5", "1.5.1.6", "1.5.1.7",
  "1.5.1.8", "1.5.1.9", "1.5.1.10", "1.5.1.11", "1.5.1.12",
  "1.5.1.15", "1.5.1.16", "1.5.1.17", "1.5.1.18", "1.5.1.19",
  "1.5.1.20", "1.5.1.21", "1.5.1.22", "1.5.1.23", "1.5.1.24",
  "1.5.1.25", "1.5.1.26", "1.5.1.27", "1.5.1.28", "1.5.1.29",
  "1.5.1.30", "1.5.1.31", "1.5.1.32", "1.5.1.33", "1.5.1.34",
  "1.5.3.1", "1.5.3.2", "1.5.3.4", "1.5.3.5", "1.5.3.6",
  "1.5.3.7", "1.5.3.10", "1.5.3.11", "1.5.3.12", "1.5.4.1",
  "1.5.5.1", "1.5.7.1", "1.5.8.1", "1.5.8.2", "1.5.99.1",
  "1.5.99.2", "1.5.99.3", "1.5.99.4", "1.5.99.5", "1.5.99.6",
  "1.5.99.8", "1.5.99.9", "1.5.99.11", "1.5.99.12", "1.6.1.1",
  "1.6.1.2", "1.6.2.2", "1.6.2.4", "1.6.2.5", "1.6.2.6",
  "1.6.3.1", "1.6.5.2", "1.6.5.3", "1.6.5.4", "1.6.5.5",
  "1.6.5.6", "1.6.5.7", "1.6.6.9", "1.6.99.1", "1.6.99.3",
  "1.6.99.5", "1.6.99.6", "1.7.1.1", "1.7.1.2", "1.7.1.3",
  "1.7.1.4", "1.7.1.5", "1.7.1.6", "1.7.1.7", "1.7.1.9",
  "1.7.1.10", "1.7.1.11", "1.7.1.12", "1.7.1.13", "1.7.2.1",
  "1.7.2.2", "1.7.2.3", "1.7.3.1", "1.7.3.2", "1.7.3.3",
  "1.7.3.4", "1.7.3.5", "1.7.7.1", "1.7.7.2", "1.7.99.1",
  "1.7.99.4", "1.7.99.6", "1.7.99.7", "1.7.99.8", "1.8.1.2",
  "1.8.1.3", "1.8.1.4", "1.8.1.5", "1.8.1.6", "1.8.1.7",
  "1.8.1.8", "1.8.1.9", "1.8.1.10", "1.8.1.11", "1.8.1.12",
  "1.8.1.13", "1.8.1.14", "1.8.1.15", "1.8.2.1", "1.8.2.2",
  "1.8.3.1", "1.8.3.2", "1.8.3.3", "1.8.3.4", "1.8.3.5",
  "1.8.4.1", "1.8.4.2", "1.8.4.3", "1.8.4.4", "1.8.4.7",
  "1.8.4.8", "1.8.4.9", "1.8.4.10", "1.8.4.11", "1.8.4.12",
  "1.8.4.13", "1.8.4.14", "1.8.5.1", "1.8.5.2", "1.8.7.1",
  "1.8.98.1", "1.8.98.2", "1.8.99.1", "1.8.99.2", "1.8.99.3",
  "1.9.3.1", "1.9.6.1", "1.9.99.1", "1.10.1.1", "1.10.2.1",
  "1.10.2.2", "1.10.3.1", "1.10.3.2", "1.10.3.3", "1.10.3.4",
  "1.10.3.5", "1.10.3.6", "1.10.99.1", "1.10.99.2",
  "1.10.99.3", "1.11.1.1", "1.11.1.2", "1.11.1.3", "1.11.1.5",
  "1.11.1.6", "1.11.1.7", "1.11.1.8", "1.11.1.9", "1.11.1.10",
  "1.11.1.11", "1.11.1.12", "1.11.1.13", "1.11.1.14",
  "1.11.1.15", "1.11.1.16", "1.12.1.2", "1.12.1.3",
  "1.12.2.1", "1.12.5.1", "1.12.7.2", "1.12.98.1",
  "1.12.98.2", "1.12.98.3", "1.12.99.6", "1.13.11.1",
  "1.13.11.2", "1.13.11.3", "1.13.11.4", "1.13.11.5",
  "1.13.11.6", "1.13.11.8", "1.13.11.9", "1.13.11.10",
  "1.13.11.11", "1.13.11.12", "1.13.11.13", "1.13.11.14",
  "1.13.11.15", "1.13.11.16", "1.13.11.17", "1.13.11.18",
  "1.13.11.19", "1.13.11.20", "1.13.11.22", "1.13.11.23",
  "1.13.11.24", "1.13.11.25", "1.13.11.26", "1.13.11.27",
  "1.13.11.28", "1.13.11.29", "1.13.11.30", "1.13.11.31",
  "1.13.11.32", "1.13.11.33", "1.13.11.34", "1.13.11.35",
  "1.13.11.36", "1.13.11.37", "1.13.11.38", "1.13.11.39",
  "1.13.11.40", "1.13.11.41", "1.13.11.43", "1.13.11.44",
  "1.13.11.45", "1.13.11.46", "1.13.11.47", "1.13.11.48",
  "1.13.11.49", "1.13.11.50", "1.13.11.51", "1.13.11.52",
  "1.13.11.53", "1.13.11.54", "1.13.11.55", "1.13.12.1",
  "1.13.12.2", "1.13.12.3", "1.13.12.4", "1.13.12.5",
  "1.13.12.6", "1.13.12.7", "1.13.12.8", "1.13.12.9",
  "1.13.12.12", "1.13.12.13", "1.13.12.14", "1.13.99.1",
  "1.13.99.3", "1.14.11.1", "1.14.11.2", "1.14.11.3",
  "1.14.11.4", "1.14.11.6", "1.14.11.7", "1.14.11.8",
  "1.14.11.9", "1.14.11.10", "1.14.11.11", "1.14.11.12",
  "1.14.11.13", "1.14.11.14", "1.14.11.15", "1.14.11.16",
  "1.14.11.17", "1.14.11.18", "1.14.11.19", "1.14.11.20",
  "1.14.11.21", "1.14.11.22", "1.14.11.23", "1.14.11.24",
  "1.14.11.25", "1.14.11.26", "1.14.11.27", "1.14.11.28",
  "1.14.12.1", "1.14.12.3", "1.14.12.4", "1.14.12.5",
  "1.14.12.7", "1.14.12.8", "1.14.12.9", "1.14.12.10",
  "1.14.12.11", "1.14.12.12", "1.14.12.13", "1.14.12.14",
  "1.14.12.15", "1.14.12.16", "1.14.12.17", "1.14.12.18",
  "1.14.12.19", "1.14.12.20", "1.14.13.1", "1.14.13.2",
  "1.14.13.3", "1.14.13.4", "1.14.13.5", "1.14.13.6",
  "1.14.13.7", "1.14.13.8", "1.14.13.9", "1.14.13.10",
  "1.14.13.11", "1.14.13.12", "1.14.13.13", "1.14.13.14",
  "1.14.13.15", "1.14.13.16", "1.14.13.17", "1.14.13.18",
  "1.14.13.19", "1.14.13.20", "1.14.13.21", "1.14.13.22",
  "1.14.13.23", "1.14.13.24", "1.14.13.25", "1.14.13.26",
  "1.14.13.27", "1.14.13.28", "1.14.13.29", "1.14.13.30",
  "1.14.13.31", "1.14.13.32", "1.14.13.33", "1.14.13.34",
  "1.14.13.35", "1.14.13.36", "1.14.13.37", "1.14.13.38",
  "1.14.13.39", "1.14.13.40", "1.14.13.41", "1.14.13.42",
  "1.14.13.43", "1.14.13.44", "1.14.13.46", "1.14.13.47",
  "1.14.13.48", "1.14.13.49", "1.14.13.50", "1.14.13.51",
  "1.14.13.52", "1.14.13.53", "1.14.13.54", "1.14.13.55",
  "1.14.13.56", "1.14.13.57", "1.14.13.58", "1.14.13.59",
  "1.14.13.60", "1.14.13.61", "1.14.13.62", "1.14.13.63",
  "1.14.13.64", "1.14.13.66", "1.14.13.67", "1.14.13.68",
  "1.14.13.69", "1.14.13.70", "1.14.13.71", "1.14.13.72",
  "1.14.13.73", "1.14.13.74", "1.14.13.75", "1.14.13.76",
  "1.14.13.77", "1.14.13.78", "1.14.13.79", "1.14.13.80",
  "1.14.13.81", "1.14.13.82", "1.14.13.83", "1.14.13.84",
  "1.14.13.85", "1.14.13.86", "1.14.13.87", "1.14.13.88",
  "1.14.13.89", "1.14.13.90", "1.14.13.91", "1.14.13.92",
  "1.14.13.93", "1.14.13.94", "1.14.13.95", "1.14.13.96",
  "1.14.13.97", "1.14.13.98", "1.14.13.99", "1.14.13.100",
  "1.14.13.101", "1.14.13.102", "1.14.13.103", "1.14.13.104",
  "1.14.13.105", "1.14.14.1", "1.14.14.3", "1.14.14.5",
  "1.14.15.1", "1.14.15.2", "1.14.15.3", "1.14.15.4",
  "1.14.15.5", "1.14.15.6", "1.14.15.7", "1.14.16.1",
  "1.14.16.2", "1.14.16.3", "1.14.16.4", "1.14.16.5",
  "1.14.16.6", "1.14.17.1", "1.14.17.3", "1.14.17.4",
  "1.14.18.1", "1.14.18.2", "1.14.19.1", "1.14.19.2",
  "1.14.19.3", "1.14.19.4", "1.14.19.5", "1.14.19.6",
  "1.14.20.1", "1.14.21.1", "1.14.21.2", "1.14.21.3",
  "1.14.21.4", "1.14.21.5", "1.14.21.6", "1.14.99.1",
  "1.14.99.2", "1.14.99.3", "1.14.99.4", "1.14.99.7",
  "1.14.99.9", "1.14.99.10", "1.14.99.11", "1.14.99.12",
  "1.14.99.14", "1.14.99.15", "1.14.99.19", "1.14.99.20",
  "1.14.99.21", "1.14.99.22", "1.14.99.23", "1.14.99.24",
  "1.14.99.26", "1.14.99.27", "1.14.99.28", "1.14.99.29",
  "1.14.99.30", "1.14.99.31", "1.14.99.32", "1.14.99.33",
  "1.14.99.34", "1.14.99.35", "1.14.99.36", "1.14.99.37",
  "1.14.99.38", "1.15.1.1", "1.15.1.2", "1.16.1.1",
  "1.16.1.2", "1.16.1.3", "1.16.1.4", "1.16.1.5", "1.16.1.6",
  "1.16.1.7", "1.16.1.8", "1.16.3.1", "1.16.8.1", "1.17.1.1",
  "1.17.1.2", "1.17.1.3", "1.17.1.4", "1.17.1.5", "1.17.3.1",
  "1.17.3.2", "1.17.3.3", "1.17.4.1", "1.17.4.2", "1.17.4.3",
  "1.17.5.1", "1.17.99.1", "1.17.99.2", "1.17.99.3",
  "1.17.99.4", "1.17.99.5", "1.18.1.1", "1.18.1.2",
  "1.18.1.3", "1.18.1.4", "1.18.6.1", "1.19.6.1", "1.20.1.1",
  "1.20.4.1", "1.20.4.2", "1.20.98.1", "1.20.99.1",
  "1.21.3.1", "1.21.3.2", "1.21.3.3", "1.21.3.4", "1.21.3.5",
  "1.21.3.6", "1.21.4.1", "1.21.4.2", "1.21.4.3", "1.21.4.4",
  "1.21.99.1", "1.97.1.1", "1.97.1.2", "1.97.1.3", "1.97.1.4",
  "1.97.1.8", "1.97.1.9", "1.97.1.10", "1.97.1.11", "2.1.1.1",
  "2.1.1.2", "2.1.1.3", "2.1.1.4", "2.1.1.5", "2.1.1.6",
  "2.1.1.7", "2.1.1.8", "2.1.1.9", "2.1.1.10", "2.1.1.11",
  "2.1.1.12", "2.1.1.13", "2.1.1.14", "2.1.1.15", "2.1.1.16",
  "2.1.1.17", "2.1.1.18", "2.1.1.19", "2.1.1.20", "2.1.1.21",
  "2.1.1.22", "2.1.1.25", "2.1.1.26", "2.1.1.27", "2.1.1.28",
  "2.1.1.29", "2.1.1.31", "2.1.1.32", "2.1.1.33", "2.1.1.34",
  "2.1.1.35", "2.1.1.36", "2.1.1.37", "2.1.1.38", "2.1.1.39",
  "2.1.1.40", "2.1.1.41", "2.1.1.42", "2.1.1.43", "2.1.1.44",
  "2.1.1.45", "2.1.1.46", "2.1.1.47", "2.1.1.48", "2.1.1.49",
  "2.1.1.50", "2.1.1.51", "2.1.1.52", "2.1.1.53", "2.1.1.54",
  "2.1.1.55", "2.1.1.56", "2.1.1.57", "2.1.1.59", "2.1.1.60",
  "2.1.1.61", "2.1.1.62", "2.1.1.63", "2.1.1.64", "2.1.1.65",
  "2.1.1.66", "2.1.1.67", "2.1.1.68", "2.1.1.69", "2.1.1.70",
  "2.1.1.71", "2.1.1.72", "2.1.1.74", "2.1.1.75", "2.1.1.76",
  "2.1.1.77", "2.1.1.78", "2.1.1.79", "2.1.1.80", "2.1.1.82",
  "2.1.1.83", "2.1.1.84", "2.1.1.85", "2.1.1.86", "2.1.1.87",
  "2.1.1.88", "2.1.1.89", "2.1.1.90", "2.1.1.91", "2.1.1.94",
  "2.1.1.95", "2.1.1.96", "2.1.1.97", "2.1.1.98", "2.1.1.99",
  "2.1.1.100", "2.1.1.101", "2.1.1.102", "2.1.1.103",
  "2.1.1.104", "2.1.1.105", "2.1.1.106", "2.1.1.107",
  "2.1.1.108", "2.1.1.109", "2.1.1.110", "2.1.1.111",
  "2.1.1.112", "2.1.1.113", "2.1.1.114", "2.1.1.115",
  "2.1.1.116", "2.1.1.117", "2.1.1.118", "2.1.1.119",
  "2.1.1.120", "2.1.1.121", "2.1.1.122", "2.1.1.123",
  "2.1.1.124", "2.1.1.125", "2.1.1.126", "2.1.1.127",
  "2.1.1.128", "2.1.1.129", "2.1.1.130", "2.1.1.131",
  "2.1.1.132", "2.1.1.133", "2.1.1.136", "2.1.1.137",
  "2.1.1.139", "2.1.1.140", "2.1.1.141", "2.1.1.142",
  "2.1.1.143", "2.1.1.144", "2.1.1.145", "2.1.1.146",
  "2.1.1.147", "2.1.1.148", "2.1.1.149", "2.1.1.150",
  "2.1.1.151", "2.1.1.152", "2.1.1.153", "2.1.1.154",
  "2.1.1.155", "2.1.1.156", "2.1.1.157", "2.1.1.158",
  "2.1.1.159", "2.1.1.160", "2.1.1.161", "2.1.1.162",
  "2.1.2.1", "2.1.2.2", "2.1.2.3", "2.1.2.4", "2.1.2.5",
  "2.1.2.7", "2.1.2.8", "2.1.2.9", "2.1.2.10", "2.1.2.11",
  "2.1.3.1", "2.1.3.2", "2.1.3.3", "2.1.3.5", "2.1.3.6",
  "2.1.3.7", "2.1.3.8", "2.1.3.9", "2.1.4.1", "2.1.4.2",
  "2.2.1.1", "2.2.1.2", "2.2.1.3", "2.2.1.4", "2.2.1.5",
  "2.2.1.6", "2.2.1.7", "2.2.1.8", "2.2.1.9", "2.3.1.1",
  "2.3.1.2", "2.3.1.3", "2.3.1.4", "2.3.1.5", "2.3.1.6",
  "2.3.1.7", "2.3.1.8", "2.3.1.9", "2.3.1.10", "2.3.1.11",
  "2.3.1.12", "2.3.1.13", "2.3.1.14", "2.3.1.15", "2.3.1.16",
  "2.3.1.17", "2.3.1.18", "2.3.1.19", "2.3.1.20", "2.3.1.21",
  "2.3.1.22", "2.3.1.23", "2.3.1.24", "2.3.1.25", "2.3.1.26",
  "2.3.1.27", "2.3.1.28", "2.3.1.29", "2.3.1.30", "2.3.1.31",
  "2.3.1.32", "2.3.1.33", "2.3.1.34", "2.3.1.35", "2.3.1.36",
  "2.3.1.37", "2.3.1.38", "2.3.1.39", "2.3.1.40", "2.3.1.41",
  "2.3.1.42", "2.3.1.43", "2.3.1.44", "2.3.1.45", "2.3.1.46",
  "2.3.1.47", "2.3.1.48", "2.3.1.49", "2.3.1.50", "2.3.1.51",
  "2.3.1.52", "2.3.1.53", "2.3.1.54", "2.3.1.56", "2.3.1.57",
  "2.3.1.58", "2.3.1.59", "2.3.1.60", "2.3.1.61", "2.3.1.62",
  "2.3.1.63", "2.3.1.64", "2.3.1.65", "2.3.1.66", "2.3.1.67",
  "2.3.1.68", "2.3.1.69", "2.3.1.70", "2.3.1.71", "2.3.1.72",
  "2.3.1.73", "2.3.1.74", "2.3.1.75", "2.3.1.76", "2.3.1.77",
  "2.3.1.78", "2.3.1.79", "2.3.1.80", "2.3.1.81", "2.3.1.82",
  "2.3.1.83", "2.3.1.84", "2.3.1.85", "2.3.1.86", "2.3.1.87",
  "2.3.1.88", "2.3.1.89", "2.3.1.90", "2.3.1.91", "2.3.1.92",
  "2.3.1.93", "2.3.1.94", "2.3.1.95", "2.3.1.96", "2.3.1.97",
  "2.3.1.98", "2.3.1.99", "2.3.1.100", "2.3.1.101",
  "2.3.1.102", "2.3.1.103", "2.3.1.104", "2.3.1.105",
  "2.3.1.106", "2.3.1.107", "2.3.1.108", "2.3.1.109",
  "2.3.1.110", "2.3.1.111", "2.3.1.112", "2.3.1.113",
  "2.3.1.114", "2.3.1.115", "2.3.1.116", "2.3.1.117",
  "2.3.1.118", "2.3.1.119", "2.3.1.121", "2.3.1.122",
  "2.3.1.123", "2.3.1.125", "2.3.1.126", "2.3.1.127",
  "2.3.1.128", "2.3.1.129", "2.3.1.130", "2.3.1.131",
  "2.3.1.132", "2.3.1.133", "2.3.1.134", "2.3.1.135",
  "2.3.1.136", "2.3.1.137", "2.3.1.138", "2.3.1.139",
  "2.3.1.140", "2.3.1.141", "2.3.1.142", "2.3.1.143",
  "2.3.1.144", "2.3.1.145", "2.3.1.146", "2.3.1.147",
  "2.3.1.148", "2.3.1.149", "2.3.1.150", "2.3.1.151",
  "2.3.1.152", "2.3.1.153", "2.3.1.154", "2.3.1.155",
  "2.3.1.156", "2.3.1.157", "2.3.1.158", "2.3.1.159",
  "2.3.1.160", "2.3.1.161", "2.3.1.162", "2.3.1.163",
  "2.3.1.164", "2.3.1.165", "2.3.1.166", "2.3.1.167",
  "2.3.1.168", "2.3.1.169", "2.3.1.170", "2.3.1.171",
  "2.3.1.172", "2.3.1.173", "2.3.1.174", "2.3.1.175",
  "2.3.1.176", "2.3.1.177", "2.3.1.178", "2.3.1.179",
  "2.3.1.180", "2.3.1.181", "2.3.1.182", "2.3.1.183",
  "2.3.1.184", "2.3.1.185", "2.3.1.186", "2.3.2.1", "2.3.2.2",
  "2.3.2.3", "2.3.2.4", "2.3.2.5", "2.3.2.6", "2.3.2.7",
  "2.3.2.8", "2.3.2.9", "2.3.2.10", "2.3.2.11", "2.3.2.12",
  "2.3.2.13", "2.3.2.14", "2.3.2.15", "2.3.3.1", "2.3.3.2",
  "2.3.3.3", "2.3.3.4", "2.3.3.5", "2.3.3.6", "2.3.3.7",
  "2.3.3.8", "2.3.3.9", "2.3.3.10", "2.3.3.11", "2.3.3.12",
  "2.3.3.13", "2.3.3.14", "2.3.3.15", "2.4.1.1", "2.4.1.2",
  "2.4.1.4", "2.4.1.5", "2.4.1.7", "2.4.1.8", "2.4.1.9",
  "2.4.1.10", "2.4.1.11", "2.4.1.12", "2.4.1.13", "2.4.1.14",
  "2.4.1.15", "2.4.1.16", "2.4.1.17", "2.4.1.18", "2.4.1.19",
  "2.4.1.20", "2.4.1.21", "2.4.1.22", "2.4.1.23", "2.4.1.24",
  "2.4.1.25", "2.4.1.26", "2.4.1.27", "2.4.1.28", "2.4.1.29",
  "2.4.1.30", "2.4.1.31", "2.4.1.32", "2.4.1.33", "2.4.1.34",
  "2.4.1.35", "2.4.1.36", "2.4.1.37", "2.4.1.38", "2.4.1.39",
  "2.4.1.40", "2.4.1.41", "2.4.1.43", "2.4.1.44", "2.4.1.45",
  "2.4.1.46", "2.4.1.47", "2.4.1.48", "2.4.1.49", "2.4.1.50",
  "2.4.1.52", "2.4.1.53", "2.4.1.54", "2.4.1.56", "2.4.1.57",
  "2.4.1.58", "2.4.1.60", "2.4.1.62", "2.4.1.63", "2.4.1.64",
  "2.4.1.65", "2.4.1.66", "2.4.1.67", "2.4.1.68", "2.4.1.69",
  "2.4.1.70", "2.4.1.71", "2.4.1.73", "2.4.1.74", "2.4.1.78",
  "2.4.1.79", "2.4.1.80", "2.4.1.81", "2.4.1.82", "2.4.1.83",
  "2.4.1.85", "2.4.1.86", "2.4.1.87", "2.4.1.88", "2.4.1.90",
  "2.4.1.91", "2.4.1.92", "2.4.1.94", "2.4.1.95", "2.4.1.96",
  "2.4.1.97", "2.4.1.99", "2.4.1.100", "2.4.1.101",
  "2.4.1.102", "2.4.1.103", "2.4.1.104", "2.4.1.105",
  "2.4.1.106", "2.4.1.109", "2.4.1.110", "2.4.1.111",
  "2.4.1.113", "2.4.1.114", "2.4.1.115", "2.4.1.116",
  "2.4.1.117", "2.4.1.118", "2.4.1.119", "2.4.1.120",
  "2.4.1.121", "2.4.1.122", "2.4.1.123", "2.4.1.125",
  "2.4.1.126", "2.4.1.127", "2.4.1.128", "2.4.1.129",
  "2.4.1.130", "2.4.1.131", "2.4.1.132", "2.4.1.133",
  "2.4.1.134", "2.4.1.135", "2.4.1.136", "2.4.1.137",
  "2.4.1.138", "2.4.1.139", "2.4.1.140", "2.4.1.141",
  "2.4.1.142", "2.4.1.143", "2.4.1.144", "2.4.1.145",
  "2.4.1.146", "2.4.1.147", "2.4.1.148", "2.4.1.149",
  "2.4.1.150", "2.4.1.152", "2.4.1.153", "2.4.1.155",
  "2.4.1.156", "2.4.1.157", "2.4.1.158", "2.4.1.159",
  "2.4.1.160", "2.4.1.161", "2.4.1.162", "2.4.1.163",
  "2.4.1.164", "2.4.1.165", "2.4.1.166", "2.4.1.167",
  "2.4.1.168", "2.4.1.170", "2.4.1.171", "2.4.1.172",
  "2.4.1.173", "2.4.1.174", "2.4.1.175", "2.4.1.176",
  "2.4.1.177", "2.4.1.178", "2.4.1.179", "2.4.1.180",
  "2.4.1.181", "2.4.1.182", "2.4.1.183", "2.4.1.184",
  "2.4.1.185", "2.4.1.186", "2.4.1.187", "2.4.1.188",
  "2.4.1.189", "2.4.1.190", "2.4.1.191", "2.4.1.192",
  "2.4.1.193", "2.4.1.194", "2.4.1.195", "2.4.1.196",
  "2.4.1.197", "2.4.1.198", "2.4.1.199", "2.4.1.201",
  "2.4.1.202", "2.4.1.203", "2.4.1.205", "2.4.1.206",
  "2.4.1.207", "2.4.1.208", "2.4.1.209", "2.4.1.210",
  "2.4.1.211", "2.4.1.212", "2.4.1.213", "2.4.1.214",
  "2.4.1.215", "2.4.1.216", "2.4.1.217", "2.4.1.218",
  "2.4.1.219", "2.4.1.220", "2.4.1.221", "2.4.1.222",
  "2.4.1.223", "2.4.1.224", "2.4.1.225", "2.4.1.226",
  "2.4.1.227", "2.4.1.228", "2.4.1.229", "2.4.1.230",
  "2.4.1.231", "2.4.1.232", "2.4.1.234", "2.4.1.236",
  "2.4.1.237", "2.4.1.238", "2.4.1.239", "2.4.1.240",
  "2.4.1.241", "2.4.1.242", "2.4.1.243", "2.4.1.244",
  "2.4.1.245", "2.4.2.1", "2.4.2.2", "2.4.2.3", "2.4.2.4",
  "2.4.2.5", "2.4.2.6", "2.4.2.7", "2.4.2.8", "2.4.2.9",
  "2.4.2.10", "2.4.2.11", "2.4.2.12", "2.4.2.14", "2.4.2.15",
  "2.4.2.16", "2.4.2.17", "2.4.2.18", "2.4.2.19", "2.4.2.20",
  "2.4.2.21", "2.4.2.22", "2.4.2.23", "2.4.2.24", "2.4.2.25",
  "2.4.2.26", "2.4.2.27", "2.4.2.28", "2.4.2.29", "2.4.2.30",
  "2.4.2.31", "2.4.2.32", "2.4.2.33", "2.4.2.34", "2.4.2.35",
  "2.4.2.36", "2.4.2.37", "2.4.2.38", "2.4.2.39", "2.4.2.40",
  "2.4.99.1", "2.4.99.2", "2.4.99.3", "2.4.99.4", "2.4.99.5",
  "2.4.99.6", "2.4.99.7", "2.4.99.8", "2.4.99.9", "2.4.99.10",
  "2.4.99.11", "2.5.1.1", "2.5.1.2", "2.5.1.3", "2.5.1.4",
  "2.5.1.5", "2.5.1.6", "2.5.1.7", "2.5.1.8", "2.5.1.9",
  "2.5.1.10", "2.5.1.11", "2.5.1.15", "2.5.1.16", "2.5.1.17",
  "2.5.1.18", "2.5.1.19", "2.5.1.20", "2.5.1.21", "2.5.1.22",
  "2.5.1.23", "2.5.1.24", "2.5.1.25", "2.5.1.26", "2.5.1.27",
  "2.5.1.28", "2.5.1.29", "2.5.1.30", "2.5.1.31", "2.5.1.32",
  "2.5.1.33", "2.5.1.34", "2.5.1.35", "2.5.1.36", "2.5.1.38",
  "2.5.1.39", "2.5.1.41", "2.5.1.42", "2.5.1.43", "2.5.1.44",
  "2.5.1.45", "2.5.1.46", "2.5.1.47", "2.5.1.48", "2.5.1.49",
  "2.5.1.50", "2.5.1.51", "2.5.1.52", "2.5.1.53", "2.5.1.54",
  "2.5.1.55", "2.5.1.56", "2.5.1.57", "2.5.1.58", "2.5.1.59",
  "2.5.1.60", "2.5.1.61", "2.5.1.62", "2.5.1.63", "2.5.1.65",
  "2.5.1.66", "2.5.1.67", "2.5.1.68", "2.5.1.69", "2.5.1.70",
  "2.5.1.71", "2.6.1.1", "2.6.1.2", "2.6.1.3", "2.6.1.4",
  "2.6.1.5", "2.6.1.6", "2.6.1.7", "2.6.1.8", "2.6.1.9",
  "2.6.1.11", "2.6.1.12", "2.6.1.13", "2.6.1.14", "2.6.1.15",
  "2.6.1.16", "2.6.1.17", "2.6.1.18", "2.6.1.19", "2.6.1.21",
  "2.6.1.22", "2.6.1.23", "2.6.1.24", "2.6.1.26", "2.6.1.27",
  "2.6.1.28", "2.6.1.29", "2.6.1.30", "2.6.1.31", "2.6.1.32",
  "2.6.1.33", "2.6.1.34", "2.6.1.35", "2.6.1.36", "2.6.1.37",
  "2.6.1.38", "2.6.1.39", "2.6.1.40", "2.6.1.41", "2.6.1.42",
  "2.6.1.43", "2.6.1.44", "2.6.1.45", "2.6.1.46", "2.6.1.47",
  "2.6.1.48", "2.6.1.49", "2.6.1.50", "2.6.1.51", "2.6.1.52",
  "2.6.1.54", "2.6.1.55", "2.6.1.56", "2.6.1.57", "2.6.1.58",
  "2.6.1.59", "2.6.1.60", "2.6.1.62", "2.6.1.63", "2.6.1.64",
  "2.6.1.65", "2.6.1.66", "2.6.1.67", "2.6.1.68", "2.6.1.70",
  "2.6.1.71", "2.6.1.72", "2.6.1.73", "2.6.1.74", "2.6.1.75",
  "2.6.1.76", "2.6.1.77", "2.6.1.78", "2.6.1.79", "2.6.1.80",
  "2.6.1.81", "2.6.1.82", "2.6.1.83", "2.6.1.84", "2.6.1.85",
  "2.6.1.86", "2.6.3.1", "2.6.99.1", "2.6.99.2", "2.7.1.1",
  "2.7.1.2", "2.7.1.3", "2.7.1.4", "2.7.1.5", "2.7.1.6",
  "2.7.1.7", "2.7.1.8", "2.7.1.10", "2.7.1.11", "2.7.1.12",
  "2.7.1.13", "2.7.1.14", "2.7.1.15", "2.7.1.16", "2.7.1.17",
  "2.7.1.18", "2.7.1.19", "2.7.1.20", "2.7.1.21", "2.7.1.22",
  "2.7.1.23", "2.7.1.24", "2.7.1.25", "2.7.1.26", "2.7.1.27",
  "2.7.1.28", "2.7.1.29", "2.7.1.30", "2.7.1.31", "2.7.1.32",
  "2.7.1.33", "2.7.1.34", "2.7.1.35", "2.7.1.36", "2.7.1.39",
  "2.7.1.40", "2.7.1.41", "2.7.1.42", "2.7.1.43", "2.7.1.44",
  "2.7.1.45", "2.7.1.46", "2.7.1.47", "2.7.1.48", "2.7.1.49",
  "2.7.1.50", "2.7.1.51", "2.7.1.52", "2.7.1.53", "2.7.1.54",
  "2.7.1.55", "2.7.1.56", "2.7.1.58", "2.7.1.59", "2.7.1.60",
  "2.7.1.61", "2.7.1.62", "2.7.1.63", "2.7.1.64", "2.7.1.65",
  "2.7.1.66", "2.7.1.67", "2.7.1.68", "2.7.1.69", "2.7.1.71",
  "2.7.1.72", "2.7.1.73", "2.7.1.74", "2.7.1.76", "2.7.1.77",
  "2.7.1.78", "2.7.1.79", "2.7.1.80", "2.7.1.81", "2.7.1.82",
  "2.7.1.83", "2.7.1.84", "2.7.1.85", "2.7.1.86", "2.7.1.87",
  "2.7.1.88", "2.7.1.89", "2.7.1.90", "2.7.1.91", "2.7.1.92",
  "2.7.1.93", "2.7.1.94", "2.7.1.95", "2.7.1.100",
  "2.7.1.101", "2.7.1.102", "2.7.1.103", "2.7.1.105",
  "2.7.1.106", "2.7.1.107", "2.7.1.108", "2.7.1.113",
  "2.7.1.114", "2.7.1.118", "2.7.1.119", "2.7.1.121",
  "2.7.1.122", "2.7.1.127", "2.7.1.130", "2.7.1.134",
  "2.7.1.136", "2.7.1.137", "2.7.1.138", "2.7.1.140",
  "2.7.1.142", "2.7.1.143", "2.7.1.144", "2.7.1.145",
  "2.7.1.146", "2.7.1.147", "2.7.1.148", "2.7.1.149",
  "2.7.1.150", "2.7.1.151", "2.7.1.153", "2.7.1.154",
  "2.7.1.156", "2.7.1.157", "2.7.1.158", "2.7.1.159",
  "2.7.1.160", "2.7.1.161", "2.7.1.162", "2.7.2.1", "2.7.2.2",
  "2.7.2.3", "2.7.2.4", "2.7.2.6", "2.7.2.7", "2.7.2.8",
  "2.7.2.10", "2.7.2.11", "2.7.2.12", "2.7.2.13", "2.7.2.14",
  "2.7.2.15", "2.7.3.1", "2.7.3.2", "2.7.3.3", "2.7.3.4",
  "2.7.3.5", "2.7.3.6", "2.7.3.7", "2.7.3.8", "2.7.3.9",
  "2.7.3.10", "2.7.4.1", "2.7.4.2", "2.7.4.3", "2.7.4.4",
  "2.7.4.6", "2.7.4.7", "2.7.4.8", "2.7.4.9", "2.7.4.10",
  "2.7.4.11", "2.7.4.12", "2.7.4.13", "2.7.4.14", "2.7.4.15",
  "2.7.4.16", "2.7.4.17", "2.7.4.18", "2.7.4.19", "2.7.4.20",
  "2.7.4.21", "2.7.4.22", "2.7.4.23", "2.7.4.24", "2.7.6.1",
  "2.7.6.2", "2.7.6.3", "2.7.6.4", "2.7.6.5", "2.7.7.1",
  "2.7.7.2", "2.7.7.3", "2.7.7.4", "2.7.7.5", "2.7.7.6",
  "2.7.7.7", "2.7.7.8", "2.7.7.9", "2.7.7.10", "2.7.7.11",
  "2.7.7.12", "2.7.7.13", "2.7.7.14", "2.7.7.15", "2.7.7.18",
  "2.7.7.19", "2.7.7.21", "2.7.7.22", "2.7.7.23", "2.7.7.24",
  "2.7.7.25", "2.7.7.27", "2.7.7.28", "2.7.7.30", "2.7.7.31",
  "2.7.7.32", "2.7.7.33", "2.7.7.34", "2.7.7.35", "2.7.7.36",
  "2.7.7.37", "2.7.7.38", "2.7.7.39", "2.7.7.40", "2.7.7.41",
  "2.7.7.42", "2.7.7.43", "2.7.7.44", "2.7.7.45", "2.7.7.46",
  "2.7.7.47", "2.7.7.48", "2.7.7.49", "2.7.7.50", "2.7.7.51",
  "2.7.7.52", "2.7.7.53", "2.7.7.54", "2.7.7.55", "2.7.7.56",
  "2.7.7.57", "2.7.7.58", "2.7.7.59", "2.7.7.60", "2.7.7.61",
  "2.7.7.62", "2.7.7.63", "2.7.7.64", "2.7.7.65", "2.7.7.66",
  "2.7.8.1", "2.7.8.2", "2.7.8.3", "2.7.8.4", "2.7.8.5",
  "2.7.8.6", "2.7.8.7", "2.7.8.8", "2.7.8.9", "2.7.8.10",
  "2.7.8.11", "2.7.8.12", "2.7.8.13", "2.7.8.14", "2.7.8.15",
  "2.7.8.17", "2.7.8.18", "2.7.8.19", "2.7.8.20", "2.7.8.21",
  "2.7.8.22", "2.7.8.23", "2.7.8.24", "2.7.8.25", "2.7.8.26",
  "2.7.8.27", "2.7.9.1", "2.7.9.2", "2.7.9.3", "2.7.9.4",
  "2.7.9.5", "2.7.10.1", "2.7.10.2", "2.7.11.1", "2.7.11.2",
  "2.7.11.3", "2.7.11.4", "2.7.11.5", "2.7.11.6", "2.7.11.7",
  "2.7.11.8", "2.7.11.9", "2.7.11.10", "2.7.11.11",
  "2.7.11.12", "2.7.11.13", "2.7.11.14", "2.7.11.15",
  "2.7.11.16", "2.7.11.17", "2.7.11.18", "2.7.11.19",
  "2.7.11.20", "2.7.11.21", "2.7.11.22", "2.7.11.23",
  "2.7.11.24", "2.7.11.25", "2.7.11.26", "2.7.11.27",
  "2.7.11.28", "2.7.11.29", "2.7.11.30", "2.7.11.31",
  "2.7.12.1", "2.7.12.2", "2.7.13.1", "2.7.13.2", "2.7.13.3",
  "2.7.99.1", "2.8.1.1", "2.8.1.2", "2.8.1.3", "2.8.1.4",
  "2.8.1.5", "2.8.1.6", "2.8.1.7", "2.8.1.8", "2.8.2.1",
  "2.8.2.2", "2.8.2.3", "2.8.2.4", "2.8.2.5", "2.8.2.6",
  "2.8.2.7", "2.8.2.8", "2.8.2.9", "2.8.2.10", "2.8.2.11",
  "2.8.2.13", "2.8.2.14", "2.8.2.15", "2.8.2.16", "2.8.2.17",
  "2.8.2.18", "2.8.2.19", "2.8.2.20", "2.8.2.21", "2.8.2.22",
  "2.8.2.23", "2.8.2.24", "2.8.2.25", "2.8.2.26", "2.8.2.27",
  "2.8.2.28", "2.8.2.29", "2.8.2.30", "2.8.2.31", "2.8.2.32",
  "2.8.2.33", "2.8.2.34", "2.8.3.1", "2.8.3.2", "2.8.3.3",
  "2.8.3.5", "2.8.3.6", "2.8.3.7", "2.8.3.8", "2.8.3.9",
  "2.8.3.10", "2.8.3.11", "2.8.3.12", "2.8.3.13", "2.8.3.14",
  "2.8.3.15", "2.8.3.16", "2.8.3.17", "2.8.4.1", "2.9.1.1",
  "3.1.1.1", "3.1.1.2", "3.1.1.3", "3.1.1.4", "3.1.1.5",
  "3.1.1.6", "3.1.1.7", "3.1.1.8", "3.1.1.10", "3.1.1.11",
  "3.1.1.13", "3.1.1.14", "3.1.1.15", "3.1.1.17", "3.1.1.19",
  "3.1.1.20", "3.1.1.21", "3.1.1.22", "3.1.1.23", "3.1.1.24",
  "3.1.1.25", "3.1.1.26", "3.1.1.27", "3.1.1.28", "3.1.1.29",
  "3.1.1.30", "3.1.1.31", "3.1.1.32", "3.1.1.33", "3.1.1.34",
  "3.1.1.35", "3.1.1.36", "3.1.1.37", "3.1.1.38", "3.1.1.39",
  "3.1.1.40", "3.1.1.41", "3.1.1.42", "3.1.1.43", "3.1.1.44",
  "3.1.1.45", "3.1.1.46", "3.1.1.47", "3.1.1.48", "3.1.1.49",
  "3.1.1.50", "3.1.1.51", "3.1.1.52", "3.1.1.53", "3.1.1.54",
  "3.1.1.55", "3.1.1.56", "3.1.1.57", "3.1.1.58", "3.1.1.59",
  "3.1.1.60", "3.1.1.61", "3.1.1.63", "3.1.1.64", "3.1.1.65",
  "3.1.1.66", "3.1.1.67", "3.1.1.68", "3.1.1.70", "3.1.1.71",
  "3.1.1.72", "3.1.1.73", "3.1.1.74", "3.1.1.75", "3.1.1.76",
  "3.1.1.77", "3.1.1.78", "3.1.1.79", "3.1.1.80", "3.1.1.81",
  "3.1.1.82", "3.1.1.83", "3.1.2.1", "3.1.2.2", "3.1.2.3",
  "3.1.2.4", "3.1.2.5", "3.1.2.6", "3.1.2.7", "3.1.2.10",
  "3.1.2.11", "3.1.2.12", "3.1.2.13", "3.1.2.14", "3.1.2.15",
  "3.1.2.16", "3.1.2.17", "3.1.2.18", "3.1.2.19", "3.1.2.20",
  "3.1.2.21", "3.1.2.22", "3.1.2.23", "3.1.2.25", "3.1.2.26",
  "3.1.2.27", "3.1.3.1", "3.1.3.2", "3.1.3.3", "3.1.3.4",
  "3.1.3.5", "3.1.3.6", "3.1.3.7", "3.1.3.8", "3.1.3.9",
  "3.1.3.10", "3.1.3.11", "3.1.3.12", "3.1.3.13", "3.1.3.14",
  "3.1.3.15", "3.1.3.16", "3.1.3.17", "3.1.3.18", "3.1.3.19",
  "3.1.3.20", "3.1.3.21", "3.1.3.22", "3.1.3.23", "3.1.3.24",
  "3.1.3.25", "3.1.3.26", "3.1.3.27", "3.1.3.28", "3.1.3.29",
  "3.1.3.31", "3.1.3.32", "3.1.3.33", "3.1.3.34", "3.1.3.35",
  "3.1.3.36", "3.1.3.37", "3.1.3.38", "3.1.3.39", "3.1.3.40",
  "3.1.3.41", "3.1.3.42", "3.1.3.43", "3.1.3.44", "3.1.3.45",
  "3.1.3.46", "3.1.3.47", "3.1.3.48", "3.1.3.49", "3.1.3.50",
  "3.1.3.51", "3.1.3.52", "3.1.3.53", "3.1.3.54", "3.1.3.55",
  "3.1.3.56", "3.1.3.57", "3.1.3.58", "3.1.3.59", "3.1.3.60",
  "3.1.3.62", "3.1.3.63", "3.1.3.64", "3.1.3.66", "3.1.3.67",
  "3.1.3.68", "3.1.3.69", "3.1.3.70", "3.1.3.71", "3.1.3.72",
  "3.1.3.73", "3.1.3.74", "3.1.3.75", "3.1.3.76", "3.1.3.77",
  "3.1.4.1", "3.1.4.2", "3.1.4.3", "3.1.4.4", "3.1.4.11",
  "3.1.4.12", "3.1.4.13", "3.1.4.14", "3.1.4.15", "3.1.4.16",
  "3.1.4.17", "3.1.4.35", "3.1.4.37", "3.1.4.38", "3.1.4.39",
  "3.1.4.40", "3.1.4.41", "3.1.4.42", "3.1.4.43", "3.1.4.44",
  "3.1.4.45", "3.1.4.46", "3.1.4.48", "3.1.4.49", "3.1.4.50",
  "3.1.4.51", "3.1.4.52", "3.1.4.53", "3.1.5.1", "3.1.6.1",
  "3.1.6.2", "3.1.6.3", "3.1.6.4", "3.1.6.6", "3.1.6.7",
  "3.1.6.8", "3.1.6.9", "3.1.6.10", "3.1.6.11", "3.1.6.12",
  "3.1.6.13", "3.1.6.14", "3.1.6.15", "3.1.6.16", "3.1.6.17",
  "3.1.6.18", "3.1.7.1", "3.1.7.2", "3.1.7.3", "3.1.7.4",
  "3.1.8.1", "3.1.8.2", "3.1.11.1", "3.1.11.2", "3.1.11.3",
  "3.1.11.4", "3.1.11.5", "3.1.11.6", "3.1.13.1", "3.1.13.2",
  "3.1.13.3", "3.1.13.4", "3.1.13.5", "3.1.14.1", "3.1.15.1",
  "3.1.16.1", "3.1.21.1", "3.1.21.2", "3.1.21.3", "3.1.21.4",
  "3.1.21.5", "3.1.21.6", "3.1.21.7", "3.1.22.1", "3.1.22.2",
  "3.1.22.4", "3.1.22.5", "3.1.25.1", "3.1.26.1", "3.1.26.2",
  "3.1.26.3", "3.1.26.4", "3.1.26.5", "3.1.26.6", "3.1.26.7",
  "3.1.26.8", "3.1.26.9", "3.1.26.10", "3.1.26.11",
  "3.1.26.12", "3.1.27.1", "3.1.27.2", "3.1.27.3", "3.1.27.4",
  "3.1.27.5", "3.1.27.6", "3.1.27.7", "3.1.27.8", "3.1.27.9",
  "3.1.27.10", "3.1.30.1", "3.1.30.2", "3.1.31.1", "3.2.1.1",
  "3.2.1.2", "3.2.1.3", "3.2.1.4", "3.2.1.6", "3.2.1.7",
  "3.2.1.8", "3.2.1.10", "3.2.1.11", "3.2.1.14", "3.2.1.15",
  "3.2.1.17", "3.2.1.18", "3.2.1.20", "3.2.1.21", "3.2.1.22",
  "3.2.1.23", "3.2.1.24", "3.2.1.25", "3.2.1.26", "3.2.1.28",
  "3.2.1.31", "3.2.1.32", "3.2.1.33", "3.2.1.35", "3.2.1.36",
  "3.2.1.37", "3.2.1.38", "3.2.1.39", "3.2.1.40", "3.2.1.41",
  "3.2.1.42", "3.2.1.43", "3.2.1.44", "3.2.1.45", "3.2.1.46",
  "3.2.1.47", "3.2.1.48", "3.2.1.49", "3.2.1.50", "3.2.1.51",
  "3.2.1.52", "3.2.1.53", "3.2.1.54", "3.2.1.55", "3.2.1.56",
  "3.2.1.57", "3.2.1.58", "3.2.1.59", "3.2.1.60", "3.2.1.61",
  "3.2.1.62", "3.2.1.63", "3.2.1.64", "3.2.1.65", "3.2.1.66",
  "3.2.1.67", "3.2.1.68", "3.2.1.70", "3.2.1.71", "3.2.1.72",
  "3.2.1.73", "3.2.1.74", "3.2.1.75", "3.2.1.76", "3.2.1.77",
  "3.2.1.78", "3.2.1.80", "3.2.1.81", "3.2.1.82", "3.2.1.83",
  "3.2.1.84", "3.2.1.85", "3.2.1.86", "3.2.1.87", "3.2.1.88",
  "3.2.1.89", "3.2.1.91", "3.2.1.92", "3.2.1.93", "3.2.1.94",
  "3.2.1.95", "3.2.1.96", "3.2.1.97", "3.2.1.98", "3.2.1.99",
  "3.2.1.100", "3.2.1.101", "3.2.1.102", "3.2.1.103",
  "3.2.1.104", "3.2.1.105", "3.2.1.106", "3.2.1.107",
  "3.2.1.108", "3.2.1.109", "3.2.1.110", "3.2.1.111",
  "3.2.1.112", "3.2.1.113", "3.2.1.114", "3.2.1.115",
  "3.2.1.116", "3.2.1.117", "3.2.1.118", "3.2.1.119",
  "3.2.1.120", "3.2.1.121", "3.2.1.122", "3.2.1.123",
  "3.2.1.124", "3.2.1.125", "3.2.1.126", "3.2.1.127",
  "3.2.1.128", "3.2.1.129", "3.2.1.130", "3.2.1.131",
  "3.2.1.132", "3.2.1.133", "3.2.1.134", "3.2.1.135",
  "3.2.1.136", "3.2.1.137", "3.2.1.139", "3.2.1.140",
  "3.2.1.141", "3.2.1.142", "3.2.1.143", "3.2.1.144",
  "3.2.1.145", "3.2.1.146", "3.2.1.147", "3.2.1.149",
  "3.2.1.150", "3.2.1.151", "3.2.1.152", "3.2.1.153",
  "3.2.1.154", "3.2.1.155", "3.2.1.156", "3.2.1.157",
  "3.2.1.158", "3.2.1.159", "3.2.1.160", "3.2.1.161",
  "3.2.1.162", "3.2.1.163", "3.2.1.164", "3.2.1.165",
  "3.2.2.1", "3.2.2.2", "3.2.2.3", "3.2.2.4", "3.2.2.5",
  "3.2.2.6", "3.2.2.7", "3.2.2.8", "3.2.2.9", "3.2.2.10",
  "3.2.2.11", "3.2.2.12", "3.2.2.13", "3.2.2.14", "3.2.2.15",
  "3.2.2.16", "3.2.2.17", "3.2.2.19", "3.2.2.20", "3.2.2.21",
  "3.2.2.22", "3.2.2.23", "3.2.2.24", "3.2.2.25", "3.3.1.1",
  "3.3.1.2", "3.3.2.1", "3.3.2.2", "3.3.2.4", "3.3.2.5",
  "3.3.2.6", "3.3.2.7", "3.3.2.8", "3.3.2.9", "3.3.2.10",
  "3.3.2.11", "3.4.11.1", "3.4.11.2", "3.4.11.3", "3.4.11.4",
  "3.4.11.5", "3.4.11.6", "3.4.11.7", "3.4.11.9", "3.4.11.10",
  "3.4.11.13", "3.4.11.14", "3.4.11.15", "3.4.11.16",
  "3.4.11.17", "3.4.11.18", "3.4.11.19", "3.4.11.20",
  "3.4.11.21", "3.4.11.22", "3.4.11.23", "3.4.11.24",
  "3.4.13.3", "3.4.13.4", "3.4.13.5", "3.4.13.7", "3.4.13.9",
  "3.4.13.12", "3.4.13.17", "3.4.13.18", "3.4.13.19",
  "3.4.13.20", "3.4.13.21", "3.4.13.22", "3.4.14.1",
  "3.4.14.2", "3.4.14.4", "3.4.14.5", "3.4.14.6", "3.4.14.9",
  "3.4.14.10", "3.4.14.11", "3.4.14.12", "3.4.15.1",
  "3.4.15.4", "3.4.15.5", "3.4.15.6", "3.4.16.2", "3.4.16.4",
  "3.4.16.5", "3.4.16.6", "3.4.17.1", "3.4.17.2", "3.4.17.3",
  "3.4.17.4", "3.4.17.6", "3.4.17.8", "3.4.17.10",
  "3.4.17.11", "3.4.17.12", "3.4.17.13", "3.4.17.14",
  "3.4.17.15", "3.4.17.16", "3.4.17.17", "3.4.17.18",
  "3.4.17.19", "3.4.17.20", "3.4.17.21", "3.4.17.22",
  "3.4.18.1", "3.4.19.1", "3.4.19.2", "3.4.19.3", "3.4.19.5",
  "3.4.19.6", "3.4.19.7", "3.4.19.9", "3.4.19.11",
  "3.4.19.12", "3.4.21.1", "3.4.21.2", "3.4.21.3", "3.4.21.4",
  "3.4.21.5", "3.4.21.6", "3.4.21.7", "3.4.21.9", "3.4.21.10",
  "3.4.21.12", "3.4.21.19", "3.4.21.20", "3.4.21.21",
  "3.4.21.22", "3.4.21.25", "3.4.21.26", "3.4.21.27",
  "3.4.21.32", "3.4.21.34", "3.4.21.35", "3.4.21.36",
  "3.4.21.37", "3.4.21.38", "3.4.21.39", "3.4.21.41",
  "3.4.21.42", "3.4.21.43", "3.4.21.45", "3.4.21.46",
  "3.4.21.47", "3.4.21.48", "3.4.21.49", "3.4.21.50",
  "3.4.21.53", "3.4.21.54", "3.4.21.55", "3.4.21.57",
  "3.4.21.59", "3.4.21.60", "3.4.21.61", "3.4.21.62",
  "3.4.21.63", "3.4.21.64", "3.4.21.65", "3.4.21.66",
  "3.4.21.67", "3.4.21.68", "3.4.21.69", "3.4.21.70",
  "3.4.21.71", "3.4.21.72", "3.4.21.73", "3.4.21.74",
  "3.4.21.75", "3.4.21.76", "3.4.21.77", "3.4.21.78",
  "3.4.21.79", "3.4.21.80", "3.4.21.81", "3.4.21.82",
  "3.4.21.83", "3.4.21.84", "3.4.21.85", "3.4.21.86",
  "3.4.21.88", "3.4.21.89", "3.4.21.90", "3.4.21.91",
  "3.4.21.92", "3.4.21.93", "3.4.21.94", "3.4.21.95",
  "3.4.21.96", "3.4.21.97", "3.4.21.98", "3.4.21.99",
  "3.4.21.100", "3.4.21.101", "3.4.21.102", "3.4.21.103",
  "3.4.21.104", "3.4.21.105", "3.4.21.106", "3.4.21.107",
  "3.4.21.108", "3.4.21.109", "3.4.21.110", "3.4.21.111",
  "3.4.21.112", "3.4.21.113", "3.4.21.114", "3.4.21.115",
  "3.4.21.116", "3.4.21.117", "3.4.21.118", "3.4.21.119",
  "3.4.21.120", "3.4.22.1", "3.4.22.2", "3.4.22.3",
  "3.4.22.6", "3.4.22.7", "3.4.22.8", "3.4.22.10",
  "3.4.22.14", "3.4.22.15", "3.4.22.16", "3.4.22.24",
  "3.4.22.25", "3.4.22.26", "3.4.22.27", "3.4.22.28",
  "3.4.22.29", "3.4.22.30", "3.4.22.31", "3.4.22.32",
  "3.4.22.33", "3.4.22.34", "3.4.22.35", "3.4.22.36",
  "3.4.22.37", "3.4.22.38", "3.4.22.39", "3.4.22.40",
  "3.4.22.41", "3.4.22.42", "3.4.22.43", "3.4.22.44",
  "3.4.22.45", "3.4.22.46", "3.4.22.47", "3.4.22.48",
  "3.4.22.49", "3.4.22.50", "3.4.22.51", "3.4.22.52",
  "3.4.22.53", "3.4.22.54", "3.4.22.55", "3.4.22.56",
  "3.4.22.57", "3.4.22.58", "3.4.22.59", "3.4.22.60",
  "3.4.22.61", "3.4.22.62", "3.4.22.63", "3.4.22.64",
  "3.4.22.65", "3.4.22.66", "3.4.22.67", "3.4.22.68",
  "3.4.23.1", "3.4.23.2", "3.4.23.3", "3.4.23.4", "3.4.23.5",
  "3.4.23.12", "3.4.23.15", "3.4.23.16", "3.4.23.17",
  "3.4.23.18", "3.4.23.19", "3.4.23.20", "3.4.23.21",
  "3.4.23.22", "3.4.23.23", "3.4.23.24", "3.4.23.25",
  "3.4.23.26", "3.4.23.28", "3.4.23.29", "3.4.23.30",
  "3.4.23.31", "3.4.23.32", "3.4.23.34", "3.4.23.35",
  "3.4.23.36", "3.4.23.38", "3.4.23.39", "3.4.23.40",
  "3.4.23.41", "3.4.23.42", "3.4.23.43", "3.4.23.44",
  "3.4.23.45", "3.4.23.46", "3.4.23.47", "3.4.23.48",
  "3.4.23.49", "3.4.24.1", "3.4.24.3", "3.4.24.6", "3.4.24.7",
  "3.4.24.11", "3.4.24.12", "3.4.24.13", "3.4.24.14",
  "3.4.24.15", "3.4.24.16", "3.4.24.17", "3.4.24.18",
  "3.4.24.19", "3.4.24.20", "3.4.24.21", "3.4.24.22",
  "3.4.24.23", "3.4.24.24", "3.4.24.25", "3.4.24.26",
  "3.4.24.27", "3.4.24.28", "3.4.24.29", "3.4.24.30",
  "3.4.24.31", "3.4.24.32", "3.4.24.33", "3.4.24.34",
  "3.4.24.35", "3.4.24.36", "3.4.24.37", "3.4.24.38",
  "3.4.24.39", "3.4.24.40", "3.4.24.41", "3.4.24.42",
  "3.4.24.43", "3.4.24.44", "3.4.24.45", "3.4.24.46",
  "3.4.24.47", "3.4.24.48", "3.4.24.49", "3.4.24.50",
  "3.4.24.51", "3.4.24.52", "3.4.24.53", "3.4.24.54",
  "3.4.24.55", "3.4.24.56", "3.4.24.57", "3.4.24.58",
  "3.4.24.59", "3.4.24.60", "3.4.24.61", "3.4.24.62",
  "3.4.24.63", "3.4.24.64", "3.4.24.65", "3.4.24.66",
  "3.4.24.67", "3.4.24.68", "3.4.24.69", "3.4.24.70",
  "3.4.24.71", "3.4.24.72", "3.4.24.73", "3.4.24.74",
  "3.4.24.75", "3.4.24.76", "3.4.24.77", "3.4.24.78",
  "3.4.24.79", "3.4.24.80", "3.4.24.81", "3.4.24.82",
  "3.4.24.83", "3.4.24.84", "3.4.24.85", "3.4.24.86",
  "3.4.25.1", "3.5.1.1", "3.5.1.2", "3.5.1.3", "3.5.1.4",
  "3.5.1.5", "3.5.1.6", "3.5.1.7", "3.5.1.8", "3.5.1.9",
  "3.5.1.10", "3.5.1.11", "3.5.1.12", "3.5.1.13", "3.5.1.14",
  "3.5.1.15", "3.5.1.16", "3.5.1.17", "3.5.1.18", "3.5.1.19",
  "3.5.1.20", "3.5.1.21", "3.5.1.22", "3.5.1.23", "3.5.1.24",
  "3.5.1.25", "3.5.1.26", "3.5.1.27", "3.5.1.28", "3.5.1.29",
  "3.5.1.30", "3.5.1.31", "3.5.1.32", "3.5.1.33", "3.5.1.35",
  "3.5.1.36", "3.5.1.38", "3.5.1.39", "3.5.1.40", "3.5.1.41",
  "3.5.1.42", "3.5.1.43", "3.5.1.44", "3.5.1.46", "3.5.1.47",
  "3.5.1.48", "3.5.1.49", "3.5.1.50", "3.5.1.51", "3.5.1.52",
  "3.5.1.53", "3.5.1.54", "3.5.1.55", "3.5.1.56", "3.5.1.57",
  "3.5.1.58", "3.5.1.59", "3.5.1.60", "3.5.1.61", "3.5.1.62",
  "3.5.1.63", "3.5.1.64", "3.5.1.65", "3.5.1.66", "3.5.1.67",
  "3.5.1.68", "3.5.1.69", "3.5.1.70", "3.5.1.71", "3.5.1.72",
  "3.5.1.73", "3.5.1.74", "3.5.1.75", "3.5.1.76", "3.5.1.77",
  "3.5.1.78", "3.5.1.79", "3.5.1.81", "3.5.1.82", "3.5.1.83",
  "3.5.1.84", "3.5.1.85", "3.5.1.86", "3.5.1.87", "3.5.1.88",
  "3.5.1.89", "3.5.1.90", "3.5.1.91", "3.5.1.92", "3.5.1.93",
  "3.5.1.94", "3.5.1.95", "3.5.1.96", "3.5.1.97", "3.5.1.98",
  "3.5.2.1", "3.5.2.2", "3.5.2.3", "3.5.2.4", "3.5.2.5",
  "3.5.2.6", "3.5.2.7", "3.5.2.9", "3.5.2.10", "3.5.2.11",
  "3.5.2.12", "3.5.2.13", "3.5.2.14", "3.5.2.15", "3.5.2.16",
  "3.5.2.17", "3.5.2.18", "3.5.3.1", "3.5.3.2", "3.5.3.3",
  "3.5.3.4", "3.5.3.5", "3.5.3.6", "3.5.3.7", "3.5.3.8",
  "3.5.3.9", "3.5.3.10", "3.5.3.11", "3.5.3.12", "3.5.3.13",
  "3.5.3.14", "3.5.3.15", "3.5.3.16", "3.5.3.17", "3.5.3.18",
  "3.5.3.19", "3.5.3.20", "3.5.3.21", "3.5.3.22", "3.5.3.23",
  "3.5.4.1", "3.5.4.2", "3.5.4.3", "3.5.4.4", "3.5.4.5",
  "3.5.4.6", "3.5.4.7", "3.5.4.8", "3.5.4.9", "3.5.4.10",
  "3.5.4.11", "3.5.4.12", "3.5.4.13", "3.5.4.14", "3.5.4.15",
  "3.5.4.16", "3.5.4.17", "3.5.4.18", "3.5.4.19", "3.5.4.20",
  "3.5.4.21", "3.5.4.22", "3.5.4.23", "3.5.4.24", "3.5.4.25",
  "3.5.4.26", "3.5.4.27", "3.5.4.28", "3.5.4.29", "3.5.4.30",
  "3.5.5.1", "3.5.5.2", "3.5.5.4", "3.5.5.5", "3.5.5.6",
  "3.5.5.7", "3.5.5.8", "3.5.99.1", "3.5.99.2", "3.5.99.3",
  "3.5.99.4", "3.5.99.5", "3.5.99.6", "3.5.99.7", "3.6.1.1",
  "3.6.1.2", "3.6.1.3", "3.6.1.5", "3.6.1.6", "3.6.1.7",
  "3.6.1.8", "3.6.1.9", "3.6.1.10", "3.6.1.11", "3.6.1.12",
  "3.6.1.13", "3.6.1.14", "3.6.1.15", "3.6.1.16", "3.6.1.17",
  "3.6.1.18", "3.6.1.19", "3.6.1.20", "3.6.1.21", "3.6.1.22",
  "3.6.1.23", "3.6.1.24", "3.6.1.25", "3.6.1.26", "3.6.1.27",
  "3.6.1.28", "3.6.1.29", "3.6.1.30", "3.6.1.31", "3.6.1.39",
  "3.6.1.40", "3.6.1.41", "3.6.1.42", "3.6.1.43", "3.6.1.44",
  "3.6.1.45", "3.6.1.52", "3.6.2.1", "3.6.2.2", "3.6.3.1",
  "3.6.3.2", "3.6.3.3", "3.6.3.4", "3.6.3.5", "3.6.3.6",
  "3.6.3.7", "3.6.3.8", "3.6.3.9", "3.6.3.10", "3.6.3.11",
  "3.6.3.12", "3.6.3.14", "3.6.3.15", "3.6.3.16", "3.6.3.17",
  "3.6.3.18", "3.6.3.19", "3.6.3.20", "3.6.3.21", "3.6.3.22",
  "3.6.3.23", "3.6.3.24", "3.6.3.25", "3.6.3.26", "3.6.3.27",
  "3.6.3.28", "3.6.3.29", "3.6.3.30", "3.6.3.31", "3.6.3.32",
  "3.6.3.33", "3.6.3.34", "3.6.3.35", "3.6.3.36", "3.6.3.37",
  "3.6.3.38", "3.6.3.39", "3.6.3.40", "3.6.3.41", "3.6.3.42",
  "3.6.3.43", "3.6.3.44", "3.6.3.46", "3.6.3.47", "3.6.3.48",
  "3.6.3.49", "3.6.3.50", "3.6.3.51", "3.6.3.52", "3.6.3.53",
  "3.6.4.1", "3.6.4.2", "3.6.4.3", "3.6.4.4", "3.6.4.5",
  "3.6.4.6", "3.6.4.7", "3.6.4.8", "3.6.4.9", "3.6.4.10",
  "3.6.4.11", "3.6.5.1", "3.6.5.2", "3.6.5.3", "3.6.5.4",
  "3.6.5.5", "3.6.5.6", "3.7.1.1", "3.7.1.2", "3.7.1.3",
  "3.7.1.4", "3.7.1.5", "3.7.1.6", "3.7.1.7", "3.7.1.8",
  "3.7.1.9", "3.7.1.10", "3.8.1.1", "3.8.1.2", "3.8.1.3",
  "3.8.1.5", "3.8.1.6", "3.8.1.7", "3.8.1.8", "3.8.1.9",
  "3.8.1.10", "3.8.1.11", "3.9.1.1", "3.10.1.1", "3.10.1.2",
  "3.11.1.1", "3.11.1.2", "3.11.1.3", "3.12.1.1", "3.13.1.1",
  "3.13.1.3", "4.1.1.1", "4.1.1.2", "4.1.1.3", "4.1.1.4",
  "4.1.1.5", "4.1.1.6", "4.1.1.7", "4.1.1.8", "4.1.1.9",
  "4.1.1.11", "4.1.1.12", "4.1.1.14", "4.1.1.15", "4.1.1.16",
  "4.1.1.17", "4.1.1.18", "4.1.1.19", "4.1.1.20", "4.1.1.21",
  "4.1.1.22", "4.1.1.23", "4.1.1.24", "4.1.1.25", "4.1.1.28",
  "4.1.1.29", "4.1.1.30", "4.1.1.31", "4.1.1.32", "4.1.1.33",
  "4.1.1.34", "4.1.1.35", "4.1.1.36", "4.1.1.37", "4.1.1.38",
  "4.1.1.39", "4.1.1.40", "4.1.1.41", "4.1.1.42", "4.1.1.43",
  "4.1.1.44", "4.1.1.45", "4.1.1.46", "4.1.1.47", "4.1.1.48",
  "4.1.1.49", "4.1.1.50", "4.1.1.51", "4.1.1.52", "4.1.1.53",
  "4.1.1.54", "4.1.1.55", "4.1.1.56", "4.1.1.57", "4.1.1.58",
  "4.1.1.59", "4.1.1.60", "4.1.1.61", "4.1.1.62", "4.1.1.63",
  "4.1.1.64", "4.1.1.65", "4.1.1.66", "4.1.1.67", "4.1.1.68",
  "4.1.1.69", "4.1.1.70", "4.1.1.71", "4.1.1.72", "4.1.1.73",
  "4.1.1.74", "4.1.1.75", "4.1.1.76", "4.1.1.77", "4.1.1.78",
  "4.1.1.79", "4.1.1.80", "4.1.1.81", "4.1.1.82", "4.1.1.83",
  "4.1.1.84", "4.1.1.85", "4.1.1.86", "4.1.2.2", "4.1.2.4",
  "4.1.2.5", "4.1.2.8", "4.1.2.9", "4.1.2.10", "4.1.2.11",
  "4.1.2.12", "4.1.2.13", "4.1.2.14", "4.1.2.17", "4.1.2.18",
  "4.1.2.19", "4.1.2.20", "4.1.2.21", "4.1.2.22", "4.1.2.23",
  "4.1.2.24", "4.1.2.25", "4.1.2.26", "4.1.2.27", "4.1.2.28",
  "4.1.2.29", "4.1.2.30", "4.1.2.32", "4.1.2.33", "4.1.2.34",
  "4.1.2.35", "4.1.2.36", "4.1.2.37", "4.1.2.38", "4.1.2.40",
  "4.1.2.41", "4.1.2.42", "4.1.3.1", "4.1.3.3", "4.1.3.4",
  "4.1.3.6", "4.1.3.13", "4.1.3.14", "4.1.3.16", "4.1.3.17",
  "4.1.3.22", "4.1.3.24", "4.1.3.25", "4.1.3.26", "4.1.3.27",
  "4.1.3.30", "4.1.3.32", "4.1.3.34", "4.1.3.35", "4.1.3.36",
  "4.1.3.38", "4.1.3.39", "4.1.3.40", "4.1.99.1", "4.1.99.2",
  "4.1.99.3", "4.1.99.5", "4.1.99.11", "4.1.99.12", "4.2.1.1",
  "4.2.1.2", "4.2.1.3", "4.2.1.4", "4.2.1.5", "4.2.1.6",
  "4.2.1.7", "4.2.1.8", "4.2.1.9", "4.2.1.10", "4.2.1.11",
  "4.2.1.12", "4.2.1.17", "4.2.1.18", "4.2.1.19", "4.2.1.20",
  "4.2.1.22", "4.2.1.24", "4.2.1.25", "4.2.1.27", "4.2.1.28",
  "4.2.1.30", "4.2.1.31", "4.2.1.32", "4.2.1.33", "4.2.1.34",
  "4.2.1.35", "4.2.1.36", "4.2.1.39", "4.2.1.40", "4.2.1.41",
  "4.2.1.42", "4.2.1.43", "4.2.1.44", "4.2.1.45", "4.2.1.46",
  "4.2.1.47", "4.2.1.48", "4.2.1.49", "4.2.1.50", "4.2.1.51",
  "4.2.1.52", "4.2.1.53", "4.2.1.54", "4.2.1.55", "4.2.1.56",
  "4.2.1.57", "4.2.1.58", "4.2.1.59", "4.2.1.60", "4.2.1.61",
  "4.2.1.62", "4.2.1.65", "4.2.1.66", "4.2.1.67", "4.2.1.68",
  "4.2.1.69", "4.2.1.70", "4.2.1.73", "4.2.1.74", "4.2.1.75",
  "4.2.1.76", "4.2.1.77", "4.2.1.78", "4.2.1.79", "4.2.1.80",
  "4.2.1.81", "4.2.1.82", "4.2.1.83", "4.2.1.84", "4.2.1.85",
  "4.2.1.87", "4.2.1.88", "4.2.1.89", "4.2.1.90", "4.2.1.91",
  "4.2.1.92", "4.2.1.93", "4.2.1.94", "4.2.1.95", "4.2.1.96",
  "4.2.1.97", "4.2.1.98", "4.2.1.99", "4.2.1.100",
  "4.2.1.101", "4.2.1.103", "4.2.1.104", "4.2.1.105",
  "4.2.1.106", "4.2.1.107", "4.2.1.108", "4.2.1.109",
  "4.2.1.110", "4.2.1.111", "4.2.1.112", "4.2.1.113",
  "4.2.2.1", "4.2.2.2", "4.2.2.3", "4.2.2.5", "4.2.2.6",
  "4.2.2.7", "4.2.2.8", "4.2.2.9", "4.2.2.10", "4.2.2.11",
  "4.2.2.12", "4.2.2.13", "4.2.2.14", "4.2.2.15", "4.2.2.16",
  "4.2.2.17", "4.2.2.18", "4.2.2.19", "4.2.2.20", "4.2.2.21",
  "4.2.2.22", "4.2.3.1", "4.2.3.2", "4.2.3.3", "4.2.3.4",
  "4.2.3.5", "4.2.3.6", "4.2.3.7", "4.2.3.8", "4.2.3.9",
  "4.2.3.10", "4.2.3.11", "4.2.3.12", "4.2.3.13", "4.2.3.14",
  "4.2.3.15", "4.2.3.16", "4.2.3.17", "4.2.3.18", "4.2.3.19",
  "4.2.3.20", "4.2.3.21", "4.2.3.22", "4.2.3.23", "4.2.3.24",
  "4.2.3.25", "4.2.3.26", "4.2.3.27", "4.2.3.28", "4.2.3.29",
  "4.2.3.30", "4.2.3.31", "4.2.3.32", "4.2.3.33", "4.2.3.34",
  "4.2.3.35", "4.2.3.36", "4.2.99.12", "4.2.99.18",
  "4.2.99.20", "4.3.1.1", "4.3.1.2", "4.3.1.3", "4.3.1.4",
  "4.3.1.6", "4.3.1.7", "4.3.1.9", "4.3.1.10", "4.3.1.12",
  "4.3.1.13", "4.3.1.14", "4.3.1.15", "4.3.1.16", "4.3.1.17",
  "4.3.1.18", "4.3.1.19", "4.3.1.20", "4.3.1.22", "4.3.1.23",
  "4.3.1.24", "4.3.1.25", "4.3.2.1", "4.3.2.2", "4.3.2.3",
  "4.3.2.4", "4.3.2.5", "4.3.3.1", "4.3.3.2", "4.3.3.3",
  "4.3.3.4", "4.4.1.1", "4.4.1.2", "4.4.1.3", "4.4.1.4",
  "4.4.1.5", "4.4.1.6", "4.4.1.8", "4.4.1.9", "4.4.1.10",
  "4.4.1.11", "4.4.1.13", "4.4.1.14", "4.4.1.15", "4.4.1.16",
  "4.4.1.17", "4.4.1.19", "4.4.1.20", "4.4.1.21", "4.4.1.22",
  "4.4.1.23", "4.4.1.24", "4.4.1.25", "4.5.1.1", "4.5.1.2",
  "4.5.1.3", "4.5.1.4", "4.5.1.5", "4.6.1.1", "4.6.1.2",
  "4.6.1.6", "4.6.1.12", "4.6.1.13", "4.6.1.14", "4.6.1.15",
  "4.99.1.1", "4.99.1.2", "4.99.1.3", "4.99.1.4", "4.99.1.5",
  "4.99.1.6", "4.99.1.7", "5.1.1.1", "5.1.1.2", "5.1.1.3",
  "5.1.1.4", "5.1.1.5", "5.1.1.6", "5.1.1.7", "5.1.1.8",
  "5.1.1.9", "5.1.1.10", "5.1.1.11", "5.1.1.12", "5.1.1.13",
  "5.1.1.14", "5.1.1.15", "5.1.1.16", "5.1.1.17", "5.1.1.18",
  "5.1.2.1", "5.1.2.2", "5.1.2.3", "5.1.2.4", "5.1.2.5",
  "5.1.2.6", "5.1.3.1", "5.1.3.2", "5.1.3.3", "5.1.3.4",
  "5.1.3.5", "5.1.3.6", "5.1.3.7", "5.1.3.8", "5.1.3.9",
  "5.1.3.10", "5.1.3.11", "5.1.3.12", "5.1.3.13", "5.1.3.14",
  "5.1.3.15", "5.1.3.16", "5.1.3.17", "5.1.3.18", "5.1.3.19",
  "5.1.3.20", "5.1.3.21", "5.1.3.22", "5.1.3.23", "5.1.99.1",
  "5.1.99.2", "5.1.99.3", "5.1.99.4", "5.2.1.1", "5.2.1.2",
  "5.2.1.3", "5.2.1.4", "5.2.1.5", "5.2.1.6", "5.2.1.7",
  "5.2.1.8", "5.2.1.9", "5.2.1.10", "5.3.1.1", "5.3.1.3",
  "5.3.1.4", "5.3.1.5", "5.3.1.6", "5.3.1.7", "5.3.1.8",
  "5.3.1.9", "5.3.1.12", "5.3.1.13", "5.3.1.14", "5.3.1.15",
  "5.3.1.16", "5.3.1.17", "5.3.1.20", "5.3.1.21", "5.3.1.22",
  "5.3.1.23", "5.3.1.24", "5.3.1.25", "5.3.1.26", "5.3.2.1",
  "5.3.2.2", "5.3.3.1", "5.3.3.2", "5.3.3.3", "5.3.3.4",
  "5.3.3.5", "5.3.3.6", "5.3.3.7", "5.3.3.8", "5.3.3.9",
  "5.3.3.10", "5.3.3.11", "5.3.3.12", "5.3.3.13", "5.3.3.14",
  "5.3.3.15", "5.3.4.1", "5.3.99.2", "5.3.99.3", "5.3.99.4",
  "5.3.99.5", "5.3.99.6", "5.3.99.7", "5.3.99.8", "5.3.99.9",
  "5.4.1.1", "5.4.1.2", "5.4.2.1", "5.4.2.2", "5.4.2.3",
  "5.4.2.4", "5.4.2.5", "5.4.2.6", "5.4.2.7", "5.4.2.8",
  "5.4.2.9", "5.4.2.10", "5.4.3.2", "5.4.3.3", "5.4.3.4",
  "5.4.3.5", "5.4.3.6", "5.4.3.7", "5.4.3.8", "5.4.4.1",
  "5.4.4.2", "5.4.4.3", "5.4.99.1", "5.4.99.2", "5.4.99.3",
  "5.4.99.4", "5.4.99.5", "5.4.99.7", "5.4.99.8", "5.4.99.9",
  "5.4.99.11", "5.4.99.12", "5.4.99.13", "5.4.99.14",
  "5.4.99.15", "5.4.99.16", "5.4.99.17", "5.4.99.18",
  "5.5.1.1", "5.5.1.2", "5.5.1.3", "5.5.1.4", "5.5.1.5",
  "5.5.1.6", "5.5.1.7", "5.5.1.8", "5.5.1.9", "5.5.1.10",
  "5.5.1.11", "5.5.1.12", "5.5.1.13", "5.5.1.14", "5.5.1.15",
  "5.5.1.16", "5.99.1.1", "5.99.1.2", "5.99.1.3", "6.1.1.1",
  "6.1.1.2", "6.1.1.3", "6.1.1.4", "6.1.1.5", "6.1.1.6",
  "6.1.1.7", "6.1.1.9", "6.1.1.10", "6.1.1.11", "6.1.1.12",
  "6.1.1.13", "6.1.1.14", "6.1.1.15", "6.1.1.16", "6.1.1.17",
  "6.1.1.18", "6.1.1.19", "6.1.1.20", "6.1.1.21", "6.1.1.22",
  "6.1.1.23", "6.1.1.24", "6.1.1.25", "6.1.1.26", "6.2.1.1",
  "6.2.1.2", "6.2.1.3", "6.2.1.4", "6.2.1.5", "6.2.1.6",
  "6.2.1.7", "6.2.1.8", "6.2.1.9", "6.2.1.10", "6.2.1.11",
  "6.2.1.12", "6.2.1.13", "6.2.1.14", "6.2.1.15", "6.2.1.16",
  "6.2.1.17", "6.2.1.18", "6.2.1.19", "6.2.1.20", "6.2.1.22",
  "6.2.1.23", "6.2.1.24", "6.2.1.25", "6.2.1.26", "6.2.1.27",
  "6.2.1.28", "6.2.1.30", "6.2.1.31", "6.2.1.32", "6.2.1.33",
  "6.2.1.34", "6.3.1.1", "6.3.1.2", "6.3.1.4", "6.3.1.5",
  "6.3.1.6", "6.3.1.7", "6.3.1.8", "6.3.1.9", "6.3.1.10",
  "6.3.1.11", "6.3.1.12", "6.3.2.1", "6.3.2.2", "6.3.2.3",
  "6.3.2.4", "6.3.2.5", "6.3.2.6", "6.3.2.7", "6.3.2.8",
  "6.3.2.9", "6.3.2.10", "6.3.2.11", "6.3.2.12", "6.3.2.13",
  "6.3.2.14", "6.3.2.16", "6.3.2.17", "6.3.2.18", "6.3.2.19",
  "6.3.2.20", "6.3.2.21", "6.3.2.22", "6.3.2.23", "6.3.2.24",
  "6.3.2.25", "6.3.2.26", "6.3.2.27", "6.3.2.28", "6.3.2.29",
  "6.3.2.30", "6.3.3.1", "6.3.3.2", "6.3.3.3", "6.3.3.4",
  "6.3.4.1", "6.3.4.2", "6.3.4.3", "6.3.4.4", "6.3.4.5",
  "6.3.4.6", "6.3.4.7", "6.3.4.8", "6.3.4.9", "6.3.4.10",
  "6.3.4.11", "6.3.4.12", "6.3.4.13", "6.3.4.14", "6.3.4.15",
  "6.3.4.16", "6.3.4.17", "6.3.4.18", "6.3.5.1", "6.3.5.2",
  "6.3.5.3", "6.3.5.4", "6.3.5.5", "6.3.5.6", "6.3.5.7",
  "6.3.5.9", "6.3.5.10", "6.4.1.1", "6.4.1.2", "6.4.1.3",
  "6.4.1.4", "6.4.1.5", "6.4.1.6", "6.4.1.7", "6.5.1.1",
  "6.5.1.2", "6.5.1.3", "6.5.1.4", "6.6.1.1", "6.6.1.2",
  NULL
};

NLM_EXTERN Boolean LookForECnumberPattern (CharPtr str)

{
  Char     ch;
  Boolean  is_ambig;
  Int2     numdashes;
  Int2     numdigits;
  Int2     numperiods;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  is_ambig = FALSE;
  numperiods = 0;
  numdigits = 0;
  numdashes = 0;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      numdigits++;
      if (is_ambig) {
        is_ambig = FALSE;
        numperiods = 0;
        numdashes = 0;
      }
      ptr++;
      ch = *ptr;
    } else if (ch == '-') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == 'n') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == '.') {
      numperiods++;
      if (numdigits > 0 && numdashes > 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      } else if (numdigits == 0 && numdashes == 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      } else if (numdashes > 1) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
      }
      numdigits = 0;
      numdashes = 0;
      ptr++;
      ch = *ptr;
    } else {
      if (numperiods == 3) {
        if (numdigits > 0 && numdashes > 0) {
        is_ambig = FALSE;
        numperiods = 0;
        numdigits = 0;
        numdashes = 0;
        } else if (numdigits > 0 || numdashes == 1) return TRUE;
      }
      ptr++;
      ch = *ptr;
      is_ambig = FALSE;
      numperiods = 0;
      numdigits = 0;
      numdashes = 0;
    }
  }

  if (numperiods == 3) {
    if (numdigits > 0 && numdashes > 0) return FALSE;
    if (numdigits > 0 || numdashes == 1) return TRUE;
  }

  return FALSE;
}

static Boolean ValidateECnumber (CharPtr str)

{
  Char     ch;
  Boolean  is_ambig;
  Int2     numdashes;
  Int2     numdigits;
  Int2     numperiods;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  is_ambig = FALSE;
  numperiods = 0;
  numdigits = 0;
  numdashes = 0;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      numdigits++;
      if (is_ambig) return FALSE;
      ptr++;
      ch = *ptr;
    } else if (ch == '-') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == 'n') {
      numdashes++;
      is_ambig = TRUE;
      ptr++;
      ch = *ptr;
    } else if (ch == '.') {
      numperiods++;
      if (numdigits > 0 && numdashes > 0) return FALSE;
      if (numdigits == 0 && numdashes == 0) return FALSE;
      if (numdashes > 1) return FALSE;
      numdigits = 0;
      numdashes = 0;
      ptr++;
      ch = *ptr;
    } else {
      ptr++;
      ch = *ptr;
    }
  }

  if (numperiods == 3) {
    if (numdigits > 0 && numdashes > 0) return FALSE;
    if (numdigits > 0 || numdashes == 1) return TRUE;
  }

  return FALSE;
}

NLM_EXTERN void ECNumberFSAFreeAll (void)

{
  CtSetPtr    csp;
  TextFsaPtr  fsa;

  fsa = (TextFsaPtr) GetAppProperty ("SpecificECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("SpecificECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("AmbiguousECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("AmbiguousECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("DeletedECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("DeletedECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("ReplacedEECNumberFSA");
  if (fsa != NULL) {
    SetAppProperty ("ReplacedEECNumberFSA", NULL);
    TextFsaFree (fsa);
  }

  fsa = (TextFsaPtr) GetAppProperty ("BodiesOfWaterFSA");
  if (fsa != NULL) {
    SetAppProperty ("BodiesOfWaterFSA", NULL);
    TextFsaFree (fsa);
  }

  csp = (CtSetPtr) GetAppProperty ("CountryLatLonList");
  if (csp != NULL) {
    SetAppProperty ("CountryLatLonList", NULL);
    CtSetDataFree (csp);
  }

  ic_code_data = MemFree (ic_code_data);
  ic_code_list = ValNodeFreeData (ic_code_list);
}

static TextFsaPtr GetECNumberFSA (CharPtr prop, CharPtr file, CharPtr PNTR local, Boolean trimAtTab)

{
  FileCache   fc;
  FILE        *fp = NULL;
  TextFsaPtr  fsa;
  Int2        i;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;

  fsa = (TextFsaPtr) GetAppProperty (prop);
  if (fsa != NULL) return fsa;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
  }

  fsa = TextFsaNew ();
  if (fsa != NULL) {
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          if (trimAtTab) {
            ptr = StringChr (str, '\t');
            if (ptr != NULL) {
              *ptr = '\0';
            }
          }
          TextFsaAdd (fsa, str);
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

    } else if (local != NULL) {
      for (i = 0; local [i] != NULL; i++) {
        str = local [i];
        if (StringDoesHaveText (str)) {
          TextFsaAdd (fsa, str);
        }
      }
    }
  }

  if (fp != NULL) {
    FileClose (fp);
  }

  SetAppProperty (prop, (Pointer) fsa);

  return fsa;
}

static TextFsaPtr GetSpecificECNumberFSA (void)

{
  return (GetECNumberFSA ("SpecificECNumberFSA", "ecnum_specific.txt", ecnum_specif, FALSE));
}

static TextFsaPtr GetAmbiguousECNumberFSA (void)

{
  return (GetECNumberFSA ("AmbiguousECNumberFSA", "ecnum_ambiguous.txt", ecnum_ambig, FALSE));
}

static TextFsaPtr GetDeletedECNumberFSA (void)

{
  return (GetECNumberFSA ("DeletedECNumberFSA", "ecnum_deleted.txt", NULL, FALSE));
}

static TextFsaPtr GetReplacedECNumberFSA (void)

{
  return (GetECNumberFSA ("ReplacedEECNumberFSA", "ecnum_replaced.txt", NULL, TRUE));
}

static Boolean ECnumberNotInList (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetSpecificECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  if (matches != NULL) return FALSE;

  fsa = GetAmbiguousECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  if (matches != NULL) return FALSE;

  return TRUE;
}

static Boolean ECnumberWasDeleted (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetDeletedECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  if (matches != NULL) return TRUE;

  return FALSE;
}

static Boolean ECnumberWasReplaced (CharPtr str)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  fsa = GetReplacedECNumberFSA ();
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  for (ptr = str, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
  }
  if (matches != NULL) return TRUE;

  return FALSE;
}

static Boolean RptUnitIsBaseRange (CharPtr str, Int4Ptr fromP, Int4Ptr toP)

{
  CharPtr   ptr;
  Char      tmp [32];
  long int  val;

  if (StringLen (str) > 25) return FALSE;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  ptr = StringStr (tmp, "..");
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  if (StringHasNoText (tmp)) return FALSE;
  if (sscanf (tmp, "%ld", &val) != 1 || val < 1) return FALSE;
  if (fromP != NULL) {
    *fromP = val - 1;
  }
  ptr += 2;
  if (StringHasNoText (ptr)) return FALSE;
  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  if (toP != NULL) {
    *toP = val - 1;
  }
  return TRUE;
}


static void ValidateRptUnit (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, GBQualPtr gbqual, Int2 qual, CharPtr key)

{
  Boolean            badchars, found, just_nuc_letters, multi_rpt_unit;
  Char               ch;
  SeqMgrFeatContext  context;
  Int4               from = -1, to = -1, ffrom, fto, ftmp;
  CharPtr            ptr, tmp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || gbqual == NULL || gbqual->val == NULL || key == NULL) return;

  found = FALSE;
  multi_rpt_unit = TRUE;
  for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    if (ch <= ' ') {
      found = TRUE;
    } else if (ch == '(' || ch == ')' || ch == ',' || ch == '.' || IS_DIGIT (ch)) {
    } else {
      multi_rpt_unit = FALSE;
    }
  }
  /*
  if (found) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
  } else if ((!multi_rpt_unit) && StringLen (gbqual->val) > 48) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
  }
  */
  if (StringICmp (key,"repeat_region") == 0 && qual == GBQUAL_rpt_unit_seq && ! multi_rpt_unit) {
    if (StringLen (gbqual->val) <= SeqLocLen (sfp->location)) {
      just_nuc_letters = TRUE;
      for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
        if (StringChr ("ACGTNacgtn", ch) == NULL) {
          just_nuc_letters = FALSE;
        }
      }
      if (just_nuc_letters) {
        tmp = GetSequenceByFeature (sfp);
        if (tmp != NULL) {
          if (StringISearch (tmp, gbqual->val) == NULL) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "repeat_region /rpt_unit and underlying sequence do not match");
          }
          MemFree (tmp);
        }
      } else {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "rpt_unit_seq qualifier contains invalid characters");
      }
    } else {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "Length of rpt_unit_seq is greater than feature length");
    }
  }

  if (qual == GBQUAL_rpt_unit_range) {
    if (RptUnitIsBaseRange (gbqual->val, &from, &to)) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) == sfp) {
        if (from < context.left || from > context.right || to < context.left || to > context.right) {
          /* could be segmented sequence */
          ffrom = SeqLocStart (sfp->location);
          fto = SeqLocStop (sfp->location);
          if (ffrom > fto) {
            ftmp = ffrom;
            ffrom = fto;
            fto = ftmp;
          }
          if (from < ffrom || from > fto || to < ffrom || to > fto) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "/rpt_unit_range is not within sequence length");
          }
        }
      }
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "/rpt_unit_range is not a base range");
    }
  }

  if (qual == GBQUAL_rpt_unit_seq) {
    badchars = FALSE;
    for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      if (ch <= ' ') {
        badchars = TRUE;
      } else if (ch == '(' || ch == ')' || IS_DIGIT (ch) || IS_ALPHA (ch)) {
      } else if (ch == ',' || ch == ';') {
      } else {
        badchars = TRUE;
      }
    }
    if (badchars) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "/rpt_unit_seq has illegal characters");
    }
  }
}

static void ValidateImpFeat (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, ImpFeatPtr ifp)

{
  Int2            adv;
  BioseqPtr       bsp;
  Char            ch;
  Boolean         failed;
  Boolean         found;
  IntFuzzPtr      fuzz;
  GBQualPtr       gbqual;
  SeqMgrFeatContext gcontext;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  Int2            i;
  Int2            index;
  Boolean         just_nuc_letters;
  Boolean         just_prt_letters;
  CharPtr         key;
  size_t          len;
  Boolean         multi_compare;
  Boolean         no_white_space;
  Boolean         only_digits;
  CharPtr         ptr;
  Int2            qual;
  Char            range[32];
  ErrSev          sev;
  SeqIntPtr       sint;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  SeqPntPtr       spp;
  CharPtr         str;
  CharPtr         tmp;
  Int2            val;

  if (vsp == NULL || gcp == NULL || sfp == NULL || ifp == NULL)
    return;
  if (StringCmp (ifp->key, "-") == 0) {
    key = StringSave ("misc_feature");
  } else {
    key = StringSaveNoNull (ifp->key);
  }
  index = GBFeatKeyNameValid (&key, FALSE);
  if (index == -1) {
    if (key != NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "Unknown feature key %s", key);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "NULL feature key");
    }
  } else if (StringICmp (key, "virion") == 0 ||
             StringICmp (key, "mutation") == 0 ||
             StringICmp (key, "allele") == 0 ||
             StringICmp (key, "Import") == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_UnknownImpFeatKey, "Feature key %s is no longer legal", key);
  } else if (StringICmp (key, "polyA_site") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    if (SeqLocStart (sfp->location) != SeqLocStop (sfp->location)) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_PolyAsiteNotPoint, "PolyA_site should be a single point");
    }
  } else if (StringICmp (key, "polyA_signal") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    if (SeqLocStart (sfp->location) == SeqLocStop (sfp->location)) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_PolyAsignalNotRange, "PolyA_signal should be a range");
    }
  } else if (StringICmp (key, "mat_peptide") == 0 ||
             StringICmp (key, "sig_peptide") == 0 ||
             StringICmp (key, "transit_peptide") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be converted to the appropriate protein feature subtype");
    CheckPeptideOnCodonBoundary (vsp, gcp, sfp, key);
  } else if (StringICmp (key, "preprotein") == 0 ||
             StringICmp (key, "proprotein") == 0) {
    sev = SEV_WARNING;
    if (vsp->is_refseq_in_sep) {
      sev = SEV_ERROR;
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be converted to the appropriate protein feature subtype");
  } else if (StringICmp (key, "mRNA") == 0 ||
             StringICmp (key, "tRNA") == 0 ||
             StringICmp (key, "rRNA") == 0 ||
             StringICmp (key, "snRNA") == 0 ||
             StringICmp (key, "scRNA") == 0 ||
             StringICmp (key, "snoRNA") == 0 ||
             StringICmp (key, "misc_RNA") == 0 ||
             StringICmp (key, "precursor_RNA") == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType,
              "RNA feature should be converted to the appropriate RNA feature subtype, location should be converted manually");
  } else if (StringICmp (key, "CDS") == 0) {
    failed = TRUE;              /* impfeat CDS must be pseudo; fail if not */
    if (sfp->pseudo) {
      failed = FALSE;
    } else {
      grp = SeqMgrGetGeneXref (sfp);
      if (grp != NULL && grp->pseudo) {
        failed = FALSE;
      } else {
        gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
        if (gene != NULL) {
          if (gene->pseudo) {
            failed = FALSE;
          } else {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (grp != NULL && grp->pseudo) {
              failed = FALSE;
            }
          }
        }
      }
    }
    for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
      if (StringCmp (gbqual->qual, "translation") == 0) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpCDShasTranslation, "ImpFeat CDS with /translation found");
      }
    }
    if (failed) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_ImpCDSnotPseudo, "ImpFeat CDS should be pseudo");
    }
  } else if (StringICmp (key, "misc_feature") == 0) {
    for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
      if (StringCmp (gbqual->qual, "standard_name") == 0) {
        if (StringCmp (gbqual->val, "Vector Contamination") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_VectorContamination, "Vector Contamination region should be trimmed from sequence");
        }
      }
    }
  }
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
      continue;
    }
    val = GBQualNameValid (gbqual->qual);
    if (val == -1) {
      if (gbqual->qual != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "Unknown qualifier %s", gbqual->qual);
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "NULL qualifier");
      }
    } else if (index != -1) {
      found = FALSE;
      for (i = 0; i < ParFlat_GBFeat[index].opt_num; i++) {
        qual = ParFlat_GBFeat[index].opt_qual[i];
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
          qual = ParFlat_GBFeat[index].mand_qual[i];
          if (qual == val) {
            found = TRUE;
            break;
          }
        }
        if (!found) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnImpFeat, "Wrong qualifier %s for feature %s", gbqual->qual, key);
        }
      }
      if (gbqual->val != NULL) {
        if (val == GBQUAL_rpt_type) {
          failed = FALSE;
          tmp = StringSave (gbqual->val);
          str = tmp;
          if (*str == '(') {
            str++;
          }
          while (!StringHasNoText (str)) {
            ptr = StringChr (str, ',');
            if (ptr == NULL) {
              ptr = StringChr (str, ')');
            }
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
            }
            found = FALSE;
            for (i = 0; legal_repeat_types[i] != NULL; i++) {
              if (StringICmp (str, legal_repeat_types[i]) == 0) {
                found = TRUE;
                break;
              }
            }
            if (!found) {
              failed = TRUE;
            }
            str = ptr;
          }
          MemFree (tmp);
          if (failed) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_rpt_unit || val == GBQUAL_rpt_unit_range || val == GBQUAL_rpt_unit_seq) {
          ValidateRptUnit (vsp, gcp, sfp, gbqual, val, key);
        } else if (val == GBQUAL_label) {
          no_white_space = TRUE;
          only_digits = TRUE;
          for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
            if (IS_WHITESP (ch)) {
              no_white_space = FALSE;
            }
            if (! IS_DIGIT (ch)) {
              only_digits = FALSE;
            }
          }
          if (only_digits || (! no_white_space)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
          }
        } else if (val == GBQUAL_replace) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
          if (bsp != NULL) {
            if (ISA_na (bsp->mol)) {
              if (StringICmp (key, "variation") == 0) {
                just_nuc_letters = TRUE;
                for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                  if (StringChr ("acgt", ch) == NULL) {
                    just_nuc_letters = FALSE;
                  }
                }
                if (!just_nuc_letters) {
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                            "%s is not a legal value for qualifier %s - should only be composed of acgt unambiguous nucleotide bases",
                            gbqual->val, gbqual->qual);
                }
              } else {
                just_nuc_letters = TRUE;
                for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                  if (StringChr ("acgtmrwsykvhdbn", ch) == NULL) {
                    just_nuc_letters = FALSE;
                  }
                }
                if (!just_nuc_letters) {
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                            "%s is not a legal value for qualifier %s - should only be composed of acgtmrwsykvhdbn nucleotide bases",
                            gbqual->val, gbqual->qual);
                }
              }
            } else if (ISA_aa (bsp->mol)) {
              just_prt_letters = TRUE;
              for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
                if (StringChr ("acdefghiklmnpqrstuvwy*", ch) == NULL) {
                  just_prt_letters = FALSE;
                }
              }
              if (!just_prt_letters) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue,
                          "%s is not a legal value for qualifier %s - should only be composed of acdefghiklmnpqrstuvwy* amino acids",
                          gbqual->val, gbqual->qual);
              }
            }
            slp = sfp->location;
            fuzz = NULL;
            if (slp != NULL && slp->choice == SEQLOC_PNT) {
              spp = (SeqPntPtr) slp->data.ptrvalue;
              if (spp != NULL) {
                fuzz = spp->fuzz;
              }
            }
            if (slp != NULL && StringLen (gbqual->val) == SeqLocLen (slp) && fuzz == NULL) {
              tmp = GetSequenceByFeature (sfp);
              if (tmp != NULL) {
                if (StringICmp (tmp, gbqual->val) == 0) {
                  ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_SuspiciousQualifierValue, "/replace already matches underlying sequence (%s)", gbqual->val);
                }
                MemFree (tmp);
              }
            }
          }
        } else if (val == GBQUAL_cons_splice) {
          found = FALSE;
          for (i = 0; legal_cons_splice_strings[i] != NULL; i++) {
            if (StringICmp (gbqual->val, legal_cons_splice_strings[i]) == 0) {
              found = TRUE;
              break;
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_mobile_element) {
          found = FALSE;
          str = NULL;
          for (i = 0; legal_mobile_element_strings[i] != NULL; i++) {
            ptr = legal_mobile_element_strings[i];
            len = StringLen (ptr);
            if (StringNICmp (gbqual->val, ptr, len) == 0) {
              found = TRUE;
              str = gbqual->val + len;
              break;
            }
          }
          if (found) {
            if (StringDoesHaveText (str) && (str [0] != ':' || str [1] == '\0')) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringNICmp (gbqual->val, "other", 5) == 0) {
              if (str [0] != ':' || str [1] == '\0') {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_frequency) {
          if (StringCmp (gbqual->val, "1") == 0 || StringCmp (gbqual->val, "1.0") == 0 || StringCmp (gbqual->val, "1.00") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is a suspicious value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_compare) {
          multi_compare = FALSE;
          ptr = gbqual->val;
          ch = *ptr;
          if (ch == '(') {
            multi_compare = TRUE;
          }
          if (! multi_compare) {
            adv = ValidateAccnDotVer (gbqual->val);
            if (adv == -5) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession missing version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv == -6) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession has bad version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv != 0) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal accession for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringChr (gbqual->val, '_') != NULL) {
              if (vsp->is_insd_in_sep) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "RefSeq accession %s cannot be used for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
        }
      }
    }
  }
  if (index != -1 && ParFlat_GBFeat[index].mand_num > 0) {
    for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
      found = FALSE;
      qual = ParFlat_GBFeat[index].mand_qual[i];
      for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
        val = GBQualNameValid (gbqual->qual);
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        if (qual == GBQUAL_citation) {
          if (sfp->cit != NULL) {
            found = TRUE;
          } else if (! StringHasNoText (sfp->comment)) {
            /* RefSeq allows conflict with accession in comment instead of sfp->cit */
            if (StringICmp (key, "conflict") == 0) {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_OTHER) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        if (StringICmp (key, "conflict") == 0 || StringICmp (key, "old_sequence") == 0) {
          /* compare qualifier can now substitute for citation qualifier for conflict and old_sequence */
          for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
            if (StringICmp (gbqual->qual, "compare") == 0 && StringDoesHaveText (gbqual->val)) {
              found = TRUE;
            }
          }
        }
      }
      if (!found) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingQualOnImpFeat, "Missing qualifier %s for feature %s", ParFlat_GBQual_names[qual].name, key);
      }
    }
  }
  if (!StringHasNoText (ifp->loc)) {
    slp = sfp->location;
    if (StringStr (ifp->loc, "one-of") != NULL) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "ImpFeat loc %s has obsolete 'one-of' text for feature %s", ifp->loc, key);
    } else if (slp != NULL && slp->choice == SEQLOC_INT) {
      sint = (SeqIntPtr) slp->data.ptrvalue;
      if (sint != NULL && sint->strand != Seq_strand_minus) {
        sprintf (range, "%ld..%ld", (long) (sint->from + 1), (long) (sint->to + 1));
        if (StringCmp (ifp->loc, range) != 0) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "ImpFeat loc %s does not equal feature location %s for feature %s", ifp->loc, range, key);
        }
      }
    }
  }
  MemFree (key);
}

static void ValidateNonImpFeat (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)
{
  Int2       adv;
  BioseqPtr  bsp;
  Char       ch;
  Boolean    failed;
  Boolean    found;
  GBQualPtr  gbqual;
  Int2       i;
  Int2       index;
  CharPtr    key;
  Boolean    multi_compare;
  Boolean    no_white_space;
  Boolean    only_digits;
  CharPtr    ptr;
  Int2       qual;
  RNAGenPtr  rgp;
  RnaRefPtr  rrp;
  ErrSev     sev;
  SeqIdPtr   sip;
  CharPtr    str;
  CharPtr    tmp;
  Int2       val;

  if (vsp == NULL || gcp == NULL || sfp == NULL)
    return;
  key = StringSaveNoNull (FeatDefTypeLabel (sfp));
  if (StringCmp (key, "Gene") == 0) {
    *key = 'g';
  }
  index = GBFeatKeyNameValid (&key, FALSE);
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
      continue;
    }
    val = GBQualNameValid (gbqual->qual);
    if (val == -1) {
      if (gbqual->qual != NULL) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          if (StringCmp (gbqual->qual, "gen_map") == 0) continue;
          if (StringCmp (gbqual->qual, "cyt_map") == 0) continue;
          if (StringCmp (gbqual->qual, "rad_map") == 0) continue;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownFeatureQual, "Unknown qualifier %s", gbqual->qual);
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownFeatureQual, "NULL qualifier");
      }
    } else if (index != -1) {
      found = FALSE;
      for (i = 0; i < ParFlat_GBFeat[index].opt_num; i++) {
        qual = ParFlat_GBFeat[index].opt_qual[i];
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
          qual = ParFlat_GBFeat[index].mand_qual[i];
          if (qual == val) {
            found = TRUE;
            break;
          }
        }
        if (!found) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Wrong qualifier %s for feature %s", gbqual->qual, key);
        }
      }
      if (gbqual->val != NULL) {
        if (val == GBQUAL_rpt_type) {
          failed = FALSE;
          tmp = StringSave (gbqual->val);
          str = tmp;
          if (*str == '(') {
            str++;
          }
          while (!StringHasNoText (str)) {
            ptr = StringChr (str, ',');
            if (ptr == NULL) {
              ptr = StringChr (str, ')');
            }
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
            }
            found = FALSE;
            for (i = 0; legal_repeat_types[i] != NULL; i++) {
              if (StringICmp (str, legal_repeat_types[i]) == 0) {
                found = TRUE;
                break;
              }
            }
            if (!found) {
              failed = TRUE;
            }
            str = ptr;
          }
          MemFree (tmp);
          if (failed) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_rpt_unit || val == GBQUAL_rpt_unit_range || val == GBQUAL_rpt_unit_seq) {
          ValidateRptUnit (vsp, gcp, sfp, gbqual, val, key);
        } else if (val == GBQUAL_label) {
          no_white_space = TRUE;
          only_digits = TRUE;
          for (ptr = gbqual->val, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
            if (IS_WHITESP (ch)) {
              no_white_space = FALSE;
            }
            if (! IS_DIGIT (ch)) {
              only_digits = FALSE;
            }
          }
          if (only_digits || (! no_white_space)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal value for qualifier %s", gbqual->qual);
          }
        } else if (val == GBQUAL_cons_splice) {
          found = FALSE;
          for (i = 0; legal_cons_splice_strings[i] != NULL; i++) {
            if (StringICmp (gbqual->val, legal_cons_splice_strings[i]) == 0) {
              found = TRUE;
              break;
            }
          }
          if (!found) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal value for qualifier %s", gbqual->val, gbqual->qual);
          }
        } else if (val == GBQUAL_compare) {
          multi_compare = FALSE;
          ptr = gbqual->val;
          ch = *ptr;
          if (ch == '(') {
            multi_compare = TRUE;
          }
          if (! multi_compare) {
            adv = ValidateAccnDotVer (gbqual->val);
            if (adv == -5) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession missing version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv == -6) {
             ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s accession has bad version for qualifier %s", gbqual->val, gbqual->qual);
            } else if (adv != 0) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "%s is not a legal accession for qualifier %s", gbqual->val, gbqual->qual);
            } else if (StringChr (gbqual->val, '_') != NULL) {
              if (vsp->is_insd_in_sep) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "RefSeq accession %s cannot be used for qualifier %s", gbqual->val, gbqual->qual);
              }
            }
          }
        }
      }
    }
  }
  if (index != -1 && ParFlat_GBFeat[index].mand_num > 0) {
    for (i = 0; i < ParFlat_GBFeat[index].mand_num; i++) {
      sev = SEV_WARNING;
      found = FALSE;
      qual = ParFlat_GBFeat[index].mand_qual[i];
      for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
        val = GBQualNameValid (gbqual->qual);
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (!found) {
        if (qual == GBQUAL_citation) {
          if (sfp->cit != NULL) {
            found = TRUE;
          } else if (! StringHasNoText (sfp->comment)) {
            /* RefSeq allows conflict with accession in comment instead of sfp->cit */
            if (StringICmp (key, "conflict") == 0) {
              bsp = BioseqFindFromSeqLoc (sfp->location);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_OTHER) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        if (StringICmp (key, "conflict") == 0 || StringICmp (key, "old_sequence") == 0) {
          /* compare qualifier can now substitute for citation qualifier for conflict and old_sequence */
          for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
            if (StringICmp (gbqual->qual, "compare") == 0 && StringDoesHaveText (gbqual->val)) {
              found = TRUE;
            }
          }
        }
      }
      if (!found) {
        if (qual == GBQUAL_ncRNA_class) {
          sev = SEV_ERROR;
          if (sfp->data.choice == SEQFEAT_RNA) {
            rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
            if (rrp != NULL) {
              if (rrp->ext.choice == 3) {
                rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
                if (rgp != NULL) {
                  if (StringDoesHaveText (rgp->_class)) {
                    found = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      if (!found) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_MissingQualOnFeature, 
                  "Missing qualifier %s for feature %s", ParFlat_GBQual_names[qual].name, key);
      }
    }
  }
  if (StringICmp (key, "mat_peptide") == 0 ||
      StringICmp (key, "sig_peptide") == 0 ||
      StringICmp (key, "transit_peptide") == 0) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      if (ISA_na (bsp->mol)) {
        sev = SEV_WARNING;
        if (vsp->is_refseq_in_sep) {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be remapped to the appropriate protein bioseq");
        CheckPeptideOnCodonBoundary (vsp, gcp, sfp, key);
      }
    }
  } else if (StringICmp (key, "preprotein") == 0 ||
      StringICmp (key, "proprotein") == 0) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      if (ISA_na (bsp->mol)) {
        sev = SEV_WARNING;
        if (vsp->is_refseq_in_sep) {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InvalidForType, "Peptide processing feature should be remapped to the appropriate protein bioseq");
      }
    }
  }
  MemFree (key);
}

/* PartialAtSpliceSiteOrGap uses code taken from SpliceCheckEx */
static Boolean PartialAtSpliceSiteOrGap (ValidStructPtr vsp, SeqLocPtr head, Uint2 slpTag, BoolPtr isgapP, BoolPtr badseqP)
{
  BioseqPtr       bsp;
  Int2            residue1, residue2;
  Boolean         rsult = FALSE;
  SeqIdPtr        sip;
  SeqLocPtr       slp = NULL, first = NULL, last = NULL;
  /*
  SeqPortPtr      spp = NULL;
  */
  Uint1           strand;
  Int4            strt, stp, donor, acceptor, len;
  StreamCache     sc;
  SeqInt          sint;
  ValNode         vn;

  if (isgapP != NULL) {
    *isgapP = FALSE;
  }
  if (badseqP != NULL) {
    *badseqP = FALSE;
  }
  if (slpTag != SLP_NOSTART && slpTag != SLP_NOSTOP)
    return FALSE;
  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    if (first == NULL) {
      first = slp;
    }
    last = slp;
  }
  if (first == NULL)
    return FALSE;

  strand = SeqLocStrand (first);
  if (SeqLocStrand (last) != strand)
    return FALSE;

  if (slpTag == SLP_NOSTART) {
    slp = first;
  } else {
    slp = last;
  }
  sip = SeqLocId (slp);
  if (sip == NULL)
    return FALSE;
  
  bsp = NULL;
  if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
    bsp = BioseqLockById (sip);
  }
  if (bsp == NULL)
    return FALSE;
  len = bsp->length;

  acceptor = SeqLocStart (slp);
  donor = SeqLocStop (slp);

  if (acceptor < 0 || acceptor >= len || donor < 0 || donor >= len) {
    /*
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range,
              "Unable to check splice consensus because feature outside range of sequence");
    */
    return FALSE;
  }

  if (strand != Seq_strand_minus) {
    if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  } else {
    sint.from = 0;
    sint.to = len - 1;
    sint.strand = strand;
    sint.id = sip;
    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
    if (! StreamCacheSetup (NULL, &vn, EXPAND_GAPS_TO_DASHES, &sc)) {
      BioseqUnlock (bsp);
      return FALSE;
    }
  }
  /* spp = SeqPortNew (bsp, 0, -1, strand, Seq_code_ncbi4na); */
  BioseqUnlock (bsp);
  /*
  if (spp == NULL)
    return FALSE;
  */

  if (strand != Seq_strand_minus) {
    strt = acceptor;
    stp = donor;
  } else {
    strt = donor;
    donor = acceptor;
    acceptor = strt;
    stp = len - donor - 1;
    strt = len - acceptor - 1;
  }

  if (slpTag == SLP_NOSTOP && stp < len - 2) {
    StreamCacheSetPosition (&sc, stp + 1);
    residue1 = StreamCacheGetResidue (&sc);
    residue2 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (stp + 1), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    */
    if (residue1 == '-' && residue2 == '-') {
      if (isgapP != NULL) {
        *isgapP = TRUE;
      }
      rsult = TRUE;
    } else if (IS_residue (residue1) && IS_residue (residue2)) {
      if ((residue1 == 'G') && (residue2  == 'T')) {
        rsult = TRUE;
      } else if ((residue1 == 'G') && (residue2 == 'C')) {
        rsult = TRUE;
      }
    } else if (badseqP != NULL) {
      *badseqP = TRUE;
    }
  } else if (slpTag == SLP_NOSTART && strt > 1) {
    StreamCacheSetPosition (&sc, strt - 2);
    residue1 = StreamCacheGetResidue (&sc);
    residue2 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (strt - 2), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    */
    if (residue1 == '-' && residue2 == '-') {
      if (isgapP != NULL) {
        *isgapP = TRUE;
      }
      rsult = TRUE;
    } else if (IS_residue (residue1) && IS_residue (residue2)) {
      if ((residue1 == 'A') && (residue2 == 'G')) {
        rsult = TRUE;
      }
    } else if (badseqP != NULL) {
      *badseqP = TRUE;
    }
  }

  /* spp = SeqPortFree (spp); */
  return rsult;
}


#if 0
static void CheckTrnaCodons (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, tRNAPtr trp)
{
  Uint1           aa = 0;
  BioseqPtr       bsp;
  Int2            code = 0;
  CharPtr         codes = NULL;
  Uint1           codon [4];
  Uint1           from;
  CharPtr         gen_code_name = NULL;
  GeneticCodePtr  gncp;
  Uint2           idx;
  Int2            j;
  Int2            k;
  ErrSev          sev = SEV_ERROR;
  SeqMapTablePtr  smtp;
  Uint1           taa;
  CharPtr         three_letter_aa = NULL;
  ValNodePtr      vnp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || trp == NULL)
    return;

  aa = 0;
  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    from = 0;
    switch (trp->aatype) {
    case 0:
      from = 0;
      break;
    case 1:
      from = Seq_code_iupacaa;
      break;
    case 2:
      from = Seq_code_ncbieaa;
      break;
    case 3:
      from = Seq_code_ncbi8aa;
      break;
    case 4:
      from = Seq_code_ncbistdaa;
      break;
    default:
      break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }

  for (j = 0; j < 6; j++) {
    if (trp->codon[j] < 64) {
      if (codes == NULL) {
        bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
        /*
        sep = GetBestTopParentForData (gcp->entityID, bsp);
        code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
        */
        BioseqToGeneticCode (bsp, &code, NULL, NULL, NULL, 0, NULL);
        gncp = GeneticCodeFind (code, NULL);
        if (gncp == NULL) {
          gncp = GeneticCodeFind (1, NULL);
          code = 1;
        }
        if (gncp == NULL)
          return;
        for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == 3) {
            codes = (CharPtr) vnp->data.ptrvalue;
          }
        }
      }
      if (codes == NULL)
        return;
      taa = codes[trp->codon[j]];
      if (aa > 0 && aa != 255) {
        if (taa != aa) {
          if (aa == 'U' || aa == 'O') {
            sev = SEV_WARNING;
          }
          if (aa == 'U' && taa == '*' && trp->codon [j] == 14) {
            /* selenocysteine normally uses TGA (14), so ignore without requiring exception in record */
          } else if (aa == 'O' && taa == '*' && trp->codon [j] == 11) {
            /* pyrrolysine normally uses TAG (11) in archaebacteria, so ignore without requiring exception in record */

            /* TAA (10) is not yet known to be used for an exceptional amino acid */
          } else if (StringISearch (sfp->except_text, "modified codon recognition") == NULL) {
            codon [0] = '\0';
            if (CodonForIndex (trp->codon [j], Seq_code_iupacna, codon)) {
              for (k = 0; k < 3; k++) {
                if (codon [k] == 'T') {
                  codon [k] = 'U';
                }
              }
              codon [3] = '\0';
            } else {
              StringCpy ((CharPtr) codon, "?");
            }
            three_letter_aa = Get3LetterSymbol (NULL, Seq_code_ncbieaa, NULL, aa);
            if (StringHasNoText (three_letter_aa)) {
              three_letter_aa = "?";
            }
            for (vnp = genetic_code_name_list; vnp != NULL; vnp = vnp->next) {
              if (vnp->choice != (Uint1) code) continue;
              gen_code_name = (CharPtr) vnp->data.ptrvalue;
              break;
            }
            if (StringHasNoText (gen_code_name)) {
              gen_code_name = "?";
            }
            ValidErr (vsp, sev, ERR_SEQ_FEAT_TrnaCodonWrong,
                      "Codon recognized by tRNA (%s) does not match amino acid (%c/%s) specified by genetic code (%d/%s)",
                      (char *) codon, (char) aa, (char *) three_letter_aa, (int) code, (char *) gen_code_name);
          }
        }
      }
    } else if (trp->codon [j] < 255) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaCodon, "tRNA codon value %d is greater than maximum 63", (int) trp->codon [j]);
    }
  }

  if (sfp->pseudo) return;

  if (aa > 0 && aa != 255) {
    /* - no gaps now that O and J are added
    if (aa <= 74) {
      shift = 0;
    } else if (aa > 79) {
      shift = 2;
    } else {
      shift = 1;
    }
    */
    if (aa != '*') {
      idx = aa - (64 /* + shift */);
    } else {
      idx = 25; /* termination */
    }
    if (idx > 0 && idx < 28) {
      /* valid trna amino acid */
    } else {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Invalid tRNA amino acid");
    }
  } else {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Missing tRNA amino acid");
  }
}
#endif

static Boolean TwoListsHaveCommonItem (
  ValNodePtr list1,
  ValNodePtr list2
)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  for (vnp1 = list1; vnp1 != NULL; vnp1 = vnp1->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    if (StringHasNoText (str1)) continue;
    for (vnp2 = list2; vnp2 != NULL; vnp2 = vnp2->next) {
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (StringHasNoText (str2)) continue;
      if (StringICmp (str1, str2) == 0) return TRUE;
    }
  }

  return FALSE;
}

static void CheckTrnaCodons (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  SeqFeatPtr sfp,
  tRNAPtr trp
)

{
  Uint1           aa = 0;
  Uint1           anticodon [4];
  BioseqPtr       bsp;
  Char            ch;
  Int2            code = 0;
  CharPtr         codes = NULL;
  Uint1           codon [4];
  CharPtr         complementBase = " TVGH  CD  M KN   YSAABW R ";
  Uint1           from;
  CharPtr         gen_code_name = NULL;
  GeneticCodePtr  gncp;
  Int2            i;
  Uint2           idx;
  Uint1           index;
  Int2            j;
  Int2            k;
  Uint1           letterToComp [256];
  Char            lttr;
  Boolean         okay;
  ValNodePtr      possibles = NULL;
  ValNodePtr      recognizes = NULL;
  StreamCache     sc;
  ErrSev          sev = SEV_ERROR;
  SeqLocPtr       slp;
  SeqMapTablePtr  smtp;
  CharPtr         str;
  Uint1           taa;
  CharPtr         three_letter_aa = NULL;
  ValNodePtr      vnp;
  CharPtr         wobble = NULL;
  Boolean         rna_editing = FALSE;

  if (vsp == NULL || gcp == NULL || sfp == NULL || trp == NULL) return;

  anticodon [0] = '\0';

  /* extract indicated amino acid */

  aa = 0;
  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    from = 0;
    switch (trp->aatype) {
    case 0:
      from = 0;
      break;
    case 1:
      from = Seq_code_iupacaa;
      break;
    case 2:
      from = Seq_code_ncbieaa;
      break;
    case 3:
      from = Seq_code_ncbi8aa;
      break;
    case 4:
      from = Seq_code_ncbistdaa;
      break;
    default:
      break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }

  three_letter_aa = Get3LetterSymbol (NULL, Seq_code_ncbieaa, NULL, aa);
  if (StringHasNoText (three_letter_aa)) {
    three_letter_aa = "?";
  }

  /* find genetic code table */

  bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  BioseqToGeneticCode (bsp, &code, NULL, NULL, NULL, 0, NULL);

  gncp = GeneticCodeFind (code, NULL);
  if (gncp == NULL) {
    gncp = GeneticCodeFind (1, NULL);
    code = 1;
  }
  if (gncp == NULL) return;

  for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != 3) continue;
    codes = (CharPtr) vnp->data.ptrvalue;
    break;
  }
  if (codes == NULL) return;

  for (vnp = genetic_code_name_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != (Uint1) code) continue;
    gen_code_name = (CharPtr) vnp->data.ptrvalue;
    break;
  }
  if (StringHasNoText (gen_code_name)) {
    gen_code_name = "?";
  }

  /* set up nucleotide complementation lookup table */

  for (i = 0; i < 256; i++) {
    letterToComp [i] = '\0';
  }
  for (ch = 'A', i = 1; ch <= 'Z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      letterToComp [(int) (Uint1) ch] = lttr;
    }
  }
  for (ch = 'a', i = 1; ch <= 'z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      letterToComp [(int) (Uint1) ch] = lttr;
    }
  }

  if (StringCmp (sfp->except_text, "RNA editing") == 0) {
    rna_editing = TRUE;
  }

  /* loop through codon_recognized array */

  for (j = 0; j < 6; j++) {
    index = (Uint1) trp->codon [j];

    if (index == 255) continue;

    if (index >= 64) {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaCodon,
                "tRNA codon value %d is greater than maximum 63",
                (int) (index));
      continue;
    }
    if (rna_editing) {
      continue;
    }

    taa = codes [index];

    codon [0] = '\0';
    if (CodonForIndex (index, Seq_code_iupacna, codon)) {
      /*
      for (k = 0; k < 3; k++) {
        if (codon [k] == 'T') {
          codon [k] = 'U';
        }
      }
      */
      codon [3] = '\0';
    } else {
      StringCpy ((CharPtr) codon, "?");
    }

    /* save codon recognized and translated amino acid for anticodon reality check */

    ValNodeCopyStr (&recognizes, taa, (CharPtr) codon);

    if (aa == 0 || aa == 255) continue;
 
    /* only report if encoded amino acid does not match indicated amino acid */

    if (taa == aa) continue;

    if (aa == 'U' || aa == 'O') {
      sev = SEV_WARNING;
    }

    /* selenocysteine normally uses TGA (14), so ignore without requiring exception in record */
    if (aa == 'U' && taa == '*' && index == 14) continue;

    /* pyrrolysine normally uses TAG (11) in archaebacteria, ignore without requiring exception */
    if (aa == 'O' && taa == '*' && index == 11) continue;

    /* TAA (10) is not yet known to be used for an exceptional amino acid, but the night is young */

    /* ignore if modified codon recognition exception is present */
    if (StringISearch (sfp->except_text, "modified codon recognition") != NULL) continue;

    for (k = 0; k < 3; k++) {
      if (codon [k] == 'T') {
        codon [k] = 'U';
      }
    }
    codon [3] = '\0';

    ValidErr (vsp, sev, ERR_SEQ_FEAT_TrnaCodonWrong,
              "Codon recognized by tRNA (%s) does not match amino acid (%c/%s) specified by genetic code (%d/%s)",
              (char *) codon, (char) aa, (char *) three_letter_aa, (int) code, (char *) gen_code_name);
  }

  /* see if anticodon is compatible with codons recognized and amino acid */

  slp = trp->anticodon;
  if (slp != NULL && SeqLocLen (slp) == 3) {

    /* read sequence under anticodon */

    if (StreamCacheSetup (NULL, slp, 0, &sc)) {
      for (i = 0; i < 3; i++) {
        ch = (Char) StreamCacheGetResidue (&sc);
        anticodon [i] = ch;
      }
      anticodon [3] = '\0';

      /* reverse complement non-wobble bases */

      codon [0] = letterToComp [(int) (Uint1) anticodon [2]];
      codon [1] = letterToComp [(int) (Uint1) anticodon [1]];
      codon [3] = '\0';

      /* expand wobble base to known binding partners */

      ch = anticodon [0];
      switch (ch) {
        case 'A' :
          wobble = "ACT";
          break;
        case 'C' :
          wobble = "G";
          break;
        case 'G' :
          wobble = "CT";
          break;
        case 'T' :
          wobble = "AG";
          break;
        default :
          break;
      }

      if (wobble != NULL) {
        for (i = 0; wobble [i] != '\0'; i++) {
          codon [2] = wobble [i];
          index = IndexForCodon (codon, Seq_code_iupacna);
          if (index < 64) {
            taa = codes [index];

            /* save possible codon recognized and translated amino acid */

            ValNodeCopyStr (&possibles, taa, (CharPtr) codon);
          }
        }
      }
    }
  }

  for (k = 0; k < 3; k++) {
    if (anticodon [k] == 'T') {
      anticodon [k] = 'U';
    }
  }
  anticodon [3] = '\0';

  if (StringHasNoText ((CharPtr) anticodon)) {
    StringCpy ((CharPtr) anticodon, "?");
  }

  /* check that codons predicted from anticodon can transfer indicated amino acid */

  if (possibles != NULL) {
    okay = FALSE;
    for (vnp = possibles; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      taa = vnp->choice;
      if (taa == aa) {
        okay = TRUE;
      }
    }
    if (! okay) {
      if (aa == 'U' && StringCmp ((CharPtr) anticodon, "UCA") == 0) {
        /* ignore TGA codon for selenocysteine */
      } else if (aa == 'O' && StringCmp ((CharPtr) anticodon, "CUA") == 0) {
        /* ignore TAG codon for pyrrolysine */
      } else if (StringISearch (sfp->except_text, "modified codon recognition") == NULL &&
                 StringISearch (sfp->except_text, "RNA editing") == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAnticodonAA,
                  "Codons predicted from anticodon (%s) cannot produce amino acid (%c/%s)",
                  (char *) anticodon, (char) aa, (char *) three_letter_aa);
      }
    }
  }

  /* check that codons recognized match codons predicted from anticodon */

  if (recognizes != NULL && possibles != NULL) {
    okay = FALSE;
    if (TwoListsHaveCommonItem (recognizes, possibles)) {
      okay = TRUE;
    }
    if (! okay) {
      if (StringISearch (sfp->except_text, "RNA editing") == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadAnticodonCodon,
                  "Codon recognized cannot be produced from anticodon (%s)",
                  (char*) anticodon);
      }
    }
  }

  ValNodeFreeData (recognizes);
  ValNodeFreeData (possibles);

  if (sfp->pseudo) return;

  if (aa == 0 || aa == 255) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Missing tRNA amino acid");
    return;
  }

  /* verify that legal amino acid is indicated */

  /* - no gaps now that O and J are added
  if (aa <= 74) {
    shift = 0;
  } else if (aa > 79) {
    shift = 2;
  } else {
    shift = 1;
  }
  */
  if (aa != '*') {
    idx = aa - (64 /* + shift */);
  } else {
    idx = 25; /* termination */
  }
  if (idx == 0 || idx >= 28) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BadTrnaAA, "Invalid tRNA amino acid");
  }
}

static void CheckRnaProductType (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, RnaRefPtr rrp)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  context;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || rrp == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  switch (rrp->type) {
    case 2 : /* mRNA */
      if (mip->biomol == MOLECULE_TYPE_MRNA) return;
      break;
    case 3 : /* tRNA */
      if (mip->biomol == MOLECULE_TYPE_TRNA) return;
      break;
    case 4 : /* rRNA */
      if (mip->biomol == MOLECULE_TYPE_RRNA) return;
      break;
    default :
      return;
  }
  ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_RnaProductMismatch, "Type of RNA does not match MolInfo of product Bioseq");
}

static BioseqSetPtr GetParentNPS (BioseqPtr bsp)
{
  BioseqSetPtr    bssp;

  if (bsp == NULL)
    return NULL;
  if (bsp->idx.parenttype != OBJ_BIOSEQSET)
    return NULL;
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot && bssp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bssp->idx.parentptr;
  }
  if (bssp->_class == BioseqseqSet_class_nuc_prot)
    return bssp;
  return NULL;
}

static Boolean NucAndProtNotInNPS (BioseqPtr nuc, BioseqPtr prot)
{
  BioseqSetPtr    bssp;

  if (nuc == NULL || prot == NULL)
    return FALSE;
  bssp = GetParentNPS (nuc);
  if (bssp == NULL)
    return TRUE;
  if (GetParentNPS (prot) != bssp)
    return TRUE;
  return FALSE;
}

static Boolean CDS5primePartialTest (
  SeqFeatPtr sfp
)

{
  BioseqPtr  nbsp;
  SeqLocPtr  slp = NULL;

  if (sfp == NULL) return FALSE;
  nbsp = BioseqFindFromSeqLoc (sfp->location);
  if (nbsp != NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      if (SeqLocStrand (slp) == Seq_strand_minus) {
        if (SeqLocStop (slp) == nbsp->length - 1) {
          return TRUE;
        }
      } else {
        if (SeqLocStart (slp) == 0) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean CDS3primePartialTest (
  SeqFeatPtr sfp
)

{
  BioseqPtr  nbsp;
  SeqLocPtr  last = NULL;
  SeqLocPtr  slp = NULL;

  if (sfp == NULL) return FALSE;
  nbsp = BioseqFindFromSeqLoc (sfp->location);
  if (nbsp != NULL) {
    last = NULL;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      last = slp;
      slp = SeqLocFindNext (sfp->location, last);
    }
    if (last != NULL) {
      if (SeqLocStrand (last) == Seq_strand_minus) {
        if (SeqLocStart (last) == 0) {
          return TRUE;
        }
      } else {
        if (SeqLocStop (last) == nbsp->length - 1) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static void CheckCDSPartial (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  context;
  MolInfoPtr         mip;
  Boolean            partial5;
  Boolean            partial3;
  SeqDescrPtr        sdp;
  ErrSev             sev;

  if (vsp == NULL || sfp == NULL) return;
  if (sfp->product == NULL) return;
  if (!vsp->useSeqMgrIndexes) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  switch (mip->completeness) {
    case 0 : /* unknown */
      break;
    case 1 : /* complete */
      if (partial5 || partial3) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is partial but protein is complete");
      }
      break;
    case 2 : /* partial */
      break;
    case 3 : /* no-left */
      if (! partial5) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' complete but protein is NH2 partial");
      }
      if (partial3) {
        sev = SEV_ERROR;
        if (CDS3primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' partial but protein is NH2 partial");
      }
      break;
    case 4 : /* no-right */
      if (! partial3) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' complete but protein is CO2 partial");
      }
      if (partial5) {
        sev = SEV_ERROR;
        if (CDS5primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' partial but protein is CO2 partial");
      }
      break;
    case 5 : /* no-ends */
      if (partial5 && partial3) {
      } else if (partial5) {
        sev = SEV_ERROR;
        if (CDS5primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 5' partial but protein has neither end");
      } else if (partial3) {
        sev = SEV_ERROR;
        if (CDS3primePartialTest (sfp)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialProblem, "CDS is 3' partial but protein has neither end");
      } else {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "CDS is complete but protein has neither end");
      }
      break;
    case 6 : /* has-left */
      break;
    case 7 : /* has-right */
      break;
    default :
      break;
  }
}

static void CheckForCommonCDSProduct (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      cds;
  CdRegionPtr     crp;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  Boolean         is_nc = FALSE;
  Boolean         is_nc_gps = FALSE;
  Boolean         is_nt = FALSE;
  Boolean         is_nw = FALSE;
  BioseqPtr       nuc;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (sfp == NULL || sfp->pseudo)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp != NULL && crp->orf)
    return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || (!SeqMgrGeneIsSuppressed (grp))) {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene != NULL) {
      if (gene->pseudo) return;
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
      if (grp != NULL && grp->pseudo) return;
    }
  }
  if (sfp->product == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) {
    sip = SeqLocId (sfp->product);
    /* okay to have far RefSeq product... */
    if (sip == NULL || sip->choice != SEQID_OTHER) {
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        /* but only if genomic product set */
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
          return;
        if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
          sep = bssp->seq_set;
          if (sep != NULL && IS_Bioseq_set (sep)) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue; 
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
              return;
          }
        }
      }
      /* or just a bioseq */
      if (sep != NULL && IS_Bioseq (sep))
        return;
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingCDSproduct, "Unable to find product Bioseq from CDS feature");
    }
    return;
  }
  nuc = BioseqFindFromSeqLoc (sfp->location);
  if (nuc != NULL) {
    for (sip = nuc->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
            is_nt = TRUE;
          } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            is_nc = TRUE;
          } else if (StringNICmp (tsip->accession, "NW_", 3) == 0) {
            is_nw = TRUE;
          }
        }
      }
    }
    if (/* (is_nc || is_nw) && */ nuc->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) nuc->idx.parentptr;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
        is_nc_gps = TRUE;
      }
    }
    if (NucAndProtNotInNPS (nuc, bsp) && (! is_nt) && (! is_nc_gps)) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CDSproductPackagingProblem, "Protein product not packaged in nuc-prot set with nucleotide");
    }
  }
  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds == NULL) return;
  if (cds != sfp) {
    /* if genomic product set, with one cds on contig and one on cdna, do not report */
    sep = vsp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
        /* feature packaging test will do final contig vs. cdna check */
        if (BioseqFindFromSeqLoc (cds->location) != BioseqFindFromSeqLoc (sfp->location))
          return;
      }
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
        sep = bssp->seq_set;
        if (sep != NULL && IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue; 
          if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set)
            if (BioseqFindFromSeqLoc (cds->location) != BioseqFindFromSeqLoc (sfp->location))
              return;
        }
      }
    }

    ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_MultipleCDSproducts, "Same product Bioseq from multiple CDS features");
  }
}

static void CheckForCommonMRNAProduct (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      gene;
  GeneRefPtr      grp;
  SeqFeatPtr      mrna;
  SeqEntryPtr     oldscope;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;

  if (sfp == NULL || sfp->pseudo)
    return;
  if (!vsp->useSeqMgrIndexes)
    return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || (!SeqMgrGeneIsSuppressed (grp))) {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene == NULL || gene->pseudo)
      return;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL && grp->pseudo)
      return;
  }
  if (sfp->product == NULL)
    return;

  oldscope = SeqEntrySetScope (vsp->sep);
  bsp = BioseqFindFromSeqLoc (sfp->product);
  SeqEntrySetScope (oldscope);
  if (bsp == NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL && sip->choice == SEQID_LOCAL) {
      sep = vsp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL) {
          if (bssp->_class == BioseqseqSet_class_gen_prod_set ||
              bssp->_class == BioseqseqSet_class_other) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MissingMRNAproduct,
            "Product Bioseq of mRNA feature is not packaged in the record");
          }
        }
      }
    }
    return;
  }

  mrna = SeqMgrGetRNAgivenProduct (bsp, NULL);
  if (mrna == NULL)
    return;
  if (mrna != sfp) {
    ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_MultipleMRNAproducts, "Same product Bioseq from multiple mRNA features");
  }
}

static void CheckForBadGeneOverlap (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr      gene, operon;
  GeneRefPtr      grp;
  ErrSev          sev = /* SEV_ERROR */ SEV_WARNING;

  if (sfp == NULL)
    return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
    return;
  gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (gene != NULL)
    return;
  gene = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, SIMPLE_OVERLAP, &fcontext);
  if (gene == NULL)
    return;
  if (IsNCorNT (vsp->sep, sfp->location)) {
    sev = SEV_WARNING;
  }
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSgeneRange, "gene overlaps CDS but does not completely contain it");
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    operon = SeqMgrGetOverlappingOperon (sfp->location, &fcontext);
    if (operon != NULL)
      return;
    ValidErr (vsp, sev, ERR_SEQ_FEAT_mRNAgeneRange, "gene overlaps mRNA but does not completely contain it");
  }
}

static void CheckForBadMRNAOverlap (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         gene = NULL;
  GeneRefPtr         grp;
  SeqFeatPtr         mrna;
  Boolean            pseudo = FALSE;
  ErrSev             sev = /* SEV_ERROR */ SEV_WARNING;

  if (sfp == NULL)
    return;

  if (sfp->pseudo) {
    pseudo = TRUE;
  }
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (SeqMgrGeneIsSuppressed (grp)) {
    } else {
      if (grp->pseudo) return;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        if (StringDoesHaveText (grp->locus_tag)) {
          gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &fcontext);
        } else if (StringDoesHaveText (grp->locus)) {
          gene = SeqMgrGetFeatureByLabel (bsp, grp->locus_tag, SEQFEAT_GENE, 0, &fcontext);
        }
        if (gene != NULL) {
          grp = (GeneRefPtr) gene->data.value.ptrvalue;
          if (grp != NULL && grp->pseudo) {
            pseudo = TRUE;
          }
        }
      }
    }
  }

  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, SIMPLE_OVERLAP, &fcontext);
  if (mrna == NULL)
    return;
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, CHECK_INTERVALS, &fcontext);
  if (mrna != NULL)
    return;
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, INTERVAL_OVERLAP, &fcontext);
  if (mrna == NULL)
    return;
  if (IsNCorNTorNW (vsp->sep, sfp->location)) {
    sev = SEV_WARNING;
  }
  if (sfp->excpt) {
    sev = SEV_WARNING;
  }
  mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, LOCATION_SUBSET, &fcontext);
  if (mrna != NULL) {
    if (StringISearch (sfp->except_text, "ribosomal slippage") == NULL && StringISearch (sfp->except_text, "trans-splicing") == NULL) {
      if (pseudo) {
        ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PseudoCDSmRNArange, "mRNA contains CDS but internal intron-exon boundaries do not match");
      } else {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNArange, "mRNA contains CDS but internal intron-exon boundaries do not match");
      }
    }
  } else {
    if (pseudo) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PseudoCDSmRNArange, "mRNA overlaps or contains CDS but does not completely contain intervals");
    } else {
      ValidErr (vsp, sev, ERR_SEQ_FEAT_CDSmRNArange, "mRNA overlaps or contains CDS but does not completely contain intervals");
    }
  }
}

/*
static void CheckForBothStrands (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  Boolean    bothstrands = FALSE, bothreverse = FALSE;
  SeqLocPtr  location, slp = NULL;
  Uint1      strand;

  if (sfp == NULL)
    return;
  location = sfp->location;
  if (location == NULL)
    return;
  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_both) {
      bothstrands = TRUE;
    } else if (strand == Seq_strand_both_rev) {
      bothreverse = TRUE;
    }
  }
  if (bothstrands && bothreverse) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (forward and reverse) strands");
  } else if (bothstrands) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (forward) strands");
  } else if (bothreverse) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BothStrands, "mRNA or CDS may not be on both (reverse) strands");
  }
}
*/

static void CheckForBothOrBothRev (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  Boolean    bothstrands = FALSE, bothreverse = FALSE, iswhole = FALSE;
  SeqLocPtr  location, slp = NULL;
  CharPtr    prefix = "Feature";
  ErrSev     sev = SEV_WARNING;
  CharPtr    suffix = "";
  Uint1      strand;

  if (sfp == NULL) return;
  location = sfp->location;
  if (location == NULL) return;

  if (sfp->idx.subtype == FEATDEF_CDS) {
    sev = SEV_ERROR;
    prefix = "CDS";
  } else if (sfp->idx.subtype == FEATDEF_mRNA) {
    sev = SEV_ERROR;
    prefix = "mRNA";
  }

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    if (slp->choice == SEQLOC_WHOLE) {
      iswhole = TRUE;
    } else {
      strand = SeqLocStrand (slp);
      if (strand == Seq_strand_both) {
        bothstrands = TRUE;
      } else if (strand == Seq_strand_both_rev) {
        bothreverse = TRUE;
      }
    }
  }
  if (bothstrands && bothreverse) {
    suffix = "(forward and reverse)";
  } else if (bothstrands) {
    suffix = "(forward)";
  } else if (bothreverse) {
    suffix = "(reverse)";
  }
  if (bothstrands || bothreverse) {
    ValidErr (vsp, sev, ERR_SEQ_FEAT_BothStrands, "%s may not be on both %s strands", prefix, suffix);
  }
  if (iswhole) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WholeLocation, "%s may not have whole location", prefix);
  }
}

static Boolean OverlappingGeneIsPseudo (SeqFeatPtr sfp)
{
  SeqFeatPtr      gene;
  GeneRefPtr      grp;

  if (sfp == NULL)
    return FALSE;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (grp->pseudo)
      return TRUE;
    return FALSE;
  }
  gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene != NULL) {
    if (gene->pseudo)
      return TRUE;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL) {
      if (grp->pseudo)
        return TRUE;
    }
  }
  return FALSE;
}

static void CheckForIllegalDbxref (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, ValNodePtr dbxref)

{
  DbtagPtr    db;
  CharPtr     good;
  Int4        id;
  Boolean     is_bc;
  Boolean     is_rf;
  Boolean     is_sc;
  ValNodePtr  vnp;

  for (vnp = dbxref; vnp != NULL; vnp = vnp->next) {
    id = -1;
    db = (DbtagPtr) vnp->data.ptrvalue;
    if (db != NULL && db->db != NULL) {

      if (DbxrefIsValid (db->db, &is_rf, &is_sc, &is_bc, &good)) {
        if (is_bc) {
          if (StringHasNoText (good)) {
            good = "?";
          }
          if (is_sc && StringICmp (db->db, "taxon") == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s, but should only be used on an OrgRef",
                      db->db, good);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "Illegal db_xref type %s, legal capitalization is %s",
                      db->db, good);
          }
        } else if (is_rf) {
          if (vsp->is_refseq_in_sep || vsp->is_gps_in_sep) {
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                      "db_xref type %s is only legal for RefSeq", db->db);
          }
        } else if (is_sc && StringICmp (db->db, "taxon") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
                    "db_xref type %s should only be used on an OrgRef", db->db);
        } else {
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
      }

      /*
      dbxerr = NULL;
      valid = IsDbxrefValid (db->db, sfp, NULL,
                             GPSorRefSeq (vsp->sep, sfp->location),
                             &dbxerr);
      if (dbxerr != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, dbxerr);
        dbxerr = MemFree (dbxerr);
      }
      */              
    }
  }
}

static CharPtr plastidtxt [] = {
  "",
  "",
  "chloroplast",
  "chromoplast",
  "",
  "",
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
  ""
};

static CharPtr legal_exception_strings [] = {
  "RNA editing",
  "reasons given in citation",
  "rearrangement required for product",
  "ribosomal slippage",
  "trans-splicing",
  "alternative processing",
  "artificial frameshift",
  "nonconsensus splice site",
  "modified codon recognition",
  "alternative start codon",
  "dicistronic gene",
  "transcribed product replaced",
  "translated product replaced",
  "transcribed pseudogene",
  NULL
};

static CharPtr refseq_exception_strings [] = {
  "unclassified transcription discrepancy",
  "unclassified translation discrepancy",
  "mismatches in transcription",
  "mismatches in translation",
  "adjusted for low-quality genome",
  NULL
};

static void ValidateExceptText (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp)

{
  Boolean  found;
  Int2     i;
  CharPtr  ptr;
  Boolean  reasons_given_except = FALSE;
  Boolean  redundant_with_comment = FALSE;
  Boolean  refseq_except = FALSE;
  ErrSev   sev = SEV_ERROR;
  CharPtr  str;
  CharPtr  tmp;

  str = StringSave (sfp->except_text);
  if (str == NULL) return;
  tmp = str;
  while (! StringHasNoText (tmp)) {
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (tmp);
    found = FALSE;
    for (i = 0; legal_exception_strings[i] != NULL; i++) {
      if (StringICmp (tmp, legal_exception_strings[i]) == 0) {
        found = TRUE;
        if (StringICmp (tmp, "reasons given in citation") == 0) {
          reasons_given_except = TRUE;
        }
        break;
      }
    }
    if (!found) {
      if (GPSorRefSeq (vsp->sep, sfp->location)) {
        for (i = 0; refseq_exception_strings[i] != NULL; i++) {
          if (StringICmp (tmp, refseq_exception_strings[i]) == 0) {
            found = TRUE;
            refseq_except = TRUE;
            break;
          }
        }
      }
      if (! found) {
        if (IsNCorNT (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_ExceptionProblem, "%s is not a legal exception explanation", tmp);
      }
    }
    if (sfp->comment != NULL && StringISearch (sfp->comment, tmp) != NULL) {
      if (StringICmp (tmp, "ribosomal slippage") != 0 && StringICmp (tmp, "trans-splicing") != 0) {
        redundant_with_comment = TRUE;
      } else if (StringICmp (sfp->comment, tmp) == 0) {
        redundant_with_comment = TRUE;
      }
    }
    tmp = ptr;
  }
  MemFree (str);
  if (redundant_with_comment) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Exception explanation text is also found in feature comment");
  }
  if (refseq_except) {
    found = FALSE;
    for (i = 0; refseq_exception_strings[i] != NULL; i++) {
      if (StringICmp (sfp->except_text, refseq_exception_strings[i]) == 0) {
        found = TRUE;
        refseq_except = TRUE;
        break;
      }
    }
    if (! found) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Genome processing exception should not be combined with other explanations");
    }
  }
  if (reasons_given_except && sfp->cit == NULL) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptionProblem, "Reasons given in citation exception does not have the required citation");
  }
}

typedef struct samecds {
  Boolean               found;
  SeqMgrFeatContextPtr  gcontext;
  Uint2                 slpTag;
  Uint1                 subtype;
  Boolean               bypassGeneTest;
} SameCds, PNTR SameCdsPtr;

static Boolean LIBCALLBACK FindSameCDS (SeqFeatPtr sfp, SeqMgrFeatContextPtr ccontext)

{
  SeqMgrFeatContextPtr  gcontext;
  Int2                  i;
  SameCdsPtr            same;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  same = (SameCdsPtr) ccontext->userdata;
  gcontext = same->gcontext;
  if (gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return TRUE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return TRUE;
  }
  if (same->subtype == FEATDEF_GENE) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right) {
      same->found = TRUE;
      return FALSE;
    }
  } else if (same->subtype == FEATDEF_mRNA) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right &&
        gcontext->numivals == ccontext->numivals) {
      if (SeqLocAinB (sfp->location, gcontext->sfp->location) >= 0) {
        if (gcontext->numivals == 1) {
          same->found = TRUE;
          return FALSE;
        } else {
          for (i = 0; i < gcontext->numivals; i++) {
            if (gcontext->ivals [2 * i] != ccontext->ivals [2 * i]) return TRUE;
            if (gcontext->ivals [2 * i + 1] != ccontext->ivals [2 * i + 1]) return TRUE;
          }
          same->found = TRUE;
          return FALSE;
        }
      }
    } else if (SeqLocAinB (sfp->location, gcontext->sfp->location) > 0) {

      if (ccontext->strand == Seq_strand_minus || gcontext->strand == Seq_strand_minus) {
        if (same->slpTag == SLP_NOSTART && gcontext->partialL) {
          if (gcontext->right == ccontext->right) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->right > ccontext->right) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        } else if (same->slpTag == SLP_NOSTOP && gcontext->partialR) {
          if (gcontext->left == ccontext->left) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->left < ccontext->left) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        }

      } else {

        if (same->slpTag == SLP_NOSTART && gcontext->partialL) {
          if (gcontext->left == ccontext->left) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->left < ccontext->left) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        } else if (same->slpTag == SLP_NOSTOP && gcontext->partialR) {
          if (gcontext->right == ccontext->right) {
            same->found = TRUE;
            return FALSE;
          }
          if (gcontext->right > ccontext->right) {
            same->bypassGeneTest = TRUE;
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}

static Boolean SameAsCDS (SeqFeatPtr sfp, Uint2 slpTag, BoolPtr bypassGeneTestP)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  Boolean            cdsFilt [SEQFEAT_MAX];
  SeqMgrFeatContext  gcontext;
  SameCds            same;
  VvmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &gcontext) != sfp) return FALSE;
  same.found = FALSE;
  same.gcontext = &gcontext;
  same.slpTag = slpTag;
  same.subtype = sfp->idx.subtype;
  same.bypassGeneTest = FALSE;

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->nearbycds != NULL) {
    cds = SeqMgrGetDesiredFeature (0, bsp, 0, 0, vdp->nearbycds, &ccontext);
    if (cds != NULL && cds->idx.subtype == FEATDEF_CDS && cds == vdp->nearbycds) {
      ccontext.userdata = (Pointer) &same;
      FindSameCDS (cds, &ccontext);
      if (same.found) {
        if (bypassGeneTestP != NULL) {
          *bypassGeneTestP = same.bypassGeneTest;
        }
        return same.found;
      }
      same.bypassGeneTest = FALSE;
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0) continue;
    cds = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, &ccontext);
    if (cds == NULL || cds->idx.subtype != FEATDEF_CDS) continue;
    ccontext.userdata = (Pointer) &same;
    FindSameCDS (cds, &ccontext);
    if (same.found) {
      if (bypassGeneTestP != NULL) {
        *bypassGeneTestP = same.bypassGeneTest;
      }
      return same.found;
    }
    same.bypassGeneTest = FALSE;
  }

  MemSet ((Pointer) &cdsFilt, 0, sizeof (cdsFilt));
  cdsFilt [SEQFEAT_CDREGION] = TRUE;
  SeqMgrExploreFeatures (bsp, (Pointer) &same, FindSameCDS, sfp->location, cdsFilt, NULL);
  if (bypassGeneTestP != NULL) {
    *bypassGeneTestP = same.bypassGeneTest;
  }
  return same.found;
}

static Boolean LIBCALLBACK FindSameMRNA (SeqFeatPtr sfp, SeqMgrFeatContextPtr ccontext)

{
  SeqMgrFeatContextPtr  gcontext;
  Int2                  i;
  SameCdsPtr            same;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA) return TRUE;
  same = (SameCdsPtr) ccontext->userdata;
  gcontext = same->gcontext;
  if (gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return TRUE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return TRUE;
  }
  if (same->subtype == FEATDEF_GENE) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right) {
      same->found = TRUE;
      return FALSE;
    }
  } else if (same->subtype == FEATDEF_CDS) {
    if (gcontext->left == ccontext->left &&
        gcontext->right == ccontext->right &&
        gcontext->numivals == ccontext->numivals) {
      if (SeqLocAinB (gcontext->sfp->location, sfp->location) >= 0) {
        if (gcontext->numivals == 1) {
          same->found = TRUE;
          return FALSE;
        } else {
          for (i = 0; i < gcontext->numivals; i++) {
            if (gcontext->ivals [2 * i] != ccontext->ivals [2 * i]) return TRUE;
            if (gcontext->ivals [2 * i + 1] != ccontext->ivals [2 * i + 1]) return TRUE;
          }
          same->found = TRUE;
          return FALSE;
        }
      }
    }
  } else if (same->subtype == FEATDEF_exon) {
    if (ccontext->strand == Seq_strand_minus || gcontext->strand == Seq_strand_minus) {
      if (same->slpTag == SLP_NOSTART && ccontext->partialL) {
        if (gcontext->right == ccontext->right) {
          same->found = TRUE;
          return FALSE;
        }
      } else if (same->slpTag == SLP_NOSTOP && ccontext->partialR) {
        if (gcontext->left == ccontext->left) {
          same->found = TRUE;
          return FALSE;
        }
      }

    } else {

      if (same->slpTag == SLP_NOSTART && ccontext->partialL) {
        if (gcontext->left == ccontext->left) {
          same->found = TRUE;
          return FALSE;
        }
      } else if (same->slpTag == SLP_NOSTOP && ccontext->partialR) {
        if (gcontext->right == ccontext->right) {
          same->found = TRUE;
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

static Boolean SameAsMRNA (SeqFeatPtr sfp, Uint2 slpTag)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  mcontext;
  SeqFeatPtr         mrna;
  Boolean            mrnaFilt [FEATDEF_MAX];
  SeqMgrFeatContext  gcontext;
  SameCds            same;
  VvmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &gcontext) != sfp) return FALSE;
  same.found = FALSE;
  same.gcontext = &gcontext;
  same.slpTag = slpTag;
  same.subtype = sfp->idx.subtype;

  vdp = (VvmDataPtr) sfp->idx.scratch;
  if (vdp != NULL && vdp->nearbymrna != NULL) {
    mrna = SeqMgrGetDesiredFeature (0, bsp, 0, 0, vdp->nearbymrna, &mcontext);
    if (mrna != NULL && mrna->idx.subtype == FEATDEF_mRNA && mrna == vdp->nearbymrna) {
      mcontext.userdata = (Pointer) &same;
      FindSameMRNA (mrna, &mcontext);
      if (same.found) {
        return same.found;
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0) continue;
    mrna = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, &mcontext);
    if (mrna == NULL || mrna->idx.subtype != FEATDEF_mRNA) continue;
    mcontext.userdata = (Pointer) &same;
    FindSameMRNA (mrna, &mcontext);
    if (same.found) {
      return same.found;
    }
  }

  MemSet ((Pointer) &mrnaFilt, 0, sizeof (mrnaFilt));
  mrnaFilt [FEATDEF_mRNA] = TRUE;
  SeqMgrExploreFeatures (bsp, (Pointer) &same, FindSameMRNA, sfp->location, NULL, mrnaFilt);
  return same.found;
}

static Boolean TestSameGene (SeqMgrFeatContextPtr ccontext, SeqMgrFeatContextPtr gcontext)

{
  if (ccontext == NULL || ccontext->sfp == NULL ||
      gcontext == NULL || gcontext->sfp == NULL ||
      gcontext->ivals == NULL || ccontext->ivals == NULL) return FALSE;
  if (gcontext->strand == ccontext->strand ||
      (ccontext->strand == Seq_strand_unknown && gcontext->strand != Seq_strand_minus) ||
      (gcontext->strand == Seq_strand_unknown && ccontext->strand != Seq_strand_minus) ||
      ccontext->strand == Seq_strand_both) {
    /* test for strands from SeqMgrGetBestOverlappingFeat, keep going if okay */
  } else {
    return FALSE;
  }
  if (gcontext->left == ccontext->left &&
      gcontext->right == ccontext->right) {
    return TRUE;
  }
  return FALSE;
}

static Boolean SameAsGene (SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         gene;
  SeqMgrFeatContext  gcontext;
  GeneRefPtr         grp;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) return FALSE;
  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &ccontext) != sfp) return FALSE;
  return TestSameGene (&gcontext, &ccontext);
}

static Boolean SplicingNotExpected (SeqFeatPtr sfp)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;

  if (StringCmp (onp->div, "BCT") == 0) return TRUE;
  if (StringCmp (onp->div, "VRL") == 0) return TRUE;
  if (StringNICmp (onp->lineage, "Bacteria; ", 10) == 0) return TRUE;
  if (StringNICmp (onp->lineage, "Archaea; ", 9) == 0) return TRUE;

  return FALSE;
}

static Boolean RareConsensusNotExpected (SeqFeatPtr sfp)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return FALSE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;
  onp = orp->orgname;
  if (onp == NULL) return FALSE;

  if (StringCmp (onp->div, "PLN") != 0) return TRUE;
  if (StringNICmp (onp->lineage, "Eukaryota; Viridiplantae; ", 26) != 0) return TRUE;

  return FALSE;
}

static Boolean HasUnderscore (CharPtr str)
{
  if (StringChr(str, '_') != NULL)
    return TRUE;
  else
    return FALSE;
}

static Boolean IsUpperCaseChar (Char ch)
{
  if (StringChr("ABCDEFGHIJKLMNOPQRSTUVWXYZ",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}

/*
static Boolean IsNumericChar (Char ch)
{
  if (StringChr("0123456789",ch) != NULL)
    return TRUE;
  else
    return FALSE;
}
*/

NLM_EXTERN Boolean IsNuclAcc (CharPtr name)

{
  if (!IsUpperCaseChar (name[0]))
    return FALSE;

  if (!HasUnderscore (name))
    return FALSE;

  return TRUE;
}

static Boolean IsCddFeat (
  SeqFeatPtr sfp
)

{
  DbtagPtr    dbt;
  ValNodePtr  vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION) return FALSE;

  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringCmp (dbt->db, "CDD") == 0 || StringCmp (dbt->db, "cdd") == 0) return TRUE;
  }

  return FALSE;
}

/* derived from ValidateSeqLoc */
static void ValidateAnticodon (ValidStructPtr vsp, SeqLocPtr slp)
{
  SeqLocPtr       tmp;
  Boolean         retval = TRUE, tmpval, mixed_strand = FALSE, unmarked_strand = FALSE,
                  ordered = TRUE, adjacent = FALSE, circular = FALSE, exception = FALSE;
  CharPtr         ctmp;
  Uint1           strand1 = Seq_strand_other, strand2 = Seq_strand_other;
  Int4            from1 = -1, from2 = -1, to1 = -1, to2 = -1;
  SeqIntPtr       sip;
  SeqPntPtr       spp;
  SeqIdPtr        id1 = NULL, id2 = NULL;
  BioseqPtr       bsp;
  SeqFeatPtr      sfp = NULL;

  if (slp == NULL)
    return;

  sfp = vsp->sfp;

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL && bsp->topology == 2) {
    circular = TRUE;
  }

  tmp = NULL;

  for (tmp = SeqLocFindNext (slp, NULL); tmp != NULL; tmp = SeqLocFindNext (slp, tmp)) {
    tmpval = TRUE;
    switch (tmp->choice) {
    case SEQLOC_INT:
      sip = (SeqIntPtr) (tmp->data.ptrvalue);
      if (sip == NULL) continue;
      strand2 = sip->strand;
      id2 = sip->id;
      from2 = sip->from;
      to2 = sip->to;
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) (tmp->data.ptrvalue);
      if (spp == NULL) continue;
      strand2 = spp->strand;
      id2 = spp->id;
      from2 = spp->point;
      to2 = spp->point;
      break;
    case SEQLOC_NULL:
      continue;
    default:
      break;
    }

    if (id1 != NULL && id2 != NULL) {
      if (SeqIdForSameBioseq (id1, id2)) {
        if ((ordered) /* && (! circular) */) {
          if (strand2 == Seq_strand_minus) {
            if (to1 < to2)
              ordered = FALSE;
            if (to2 + 1 == from1)
              adjacent = TRUE;
          } else {
            if (to1 > to2)
              ordered = FALSE;
            if (to1 + 1 == from2)
              adjacent = TRUE;
          }
        }
        if (strand1 == strand2 && from1 == from2 && to1 == to2) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_DuplicateInterval, "Duplicate anticodon exons in location");
        }
      }
    }

    if (!tmpval) {
      retval = FALSE;
      ctmp = SeqLocPrint (tmp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range, "Anticodon location [%s] out of range", ctmp);
      MemFree (ctmp);
    }

    if ((strand1 != Seq_strand_other) && (strand2 != Seq_strand_other)) {
      if (SeqIdForSameBioseq (id1, id2)) {
        if (strand1 != strand2) {
          if (strand1 == Seq_strand_plus && strand2 == Seq_strand_unknown) {
            unmarked_strand = TRUE;
          } else if (strand1 == Seq_strand_unknown && strand2 == Seq_strand_plus) {
            unmarked_strand = TRUE;
          } else {
            mixed_strand = TRUE;
          }
        }
      }
    }

    from1 = from2;
    to1 = to2;
    id1 = id2;
    strand1 = strand2;
  }

  if (sfp != NULL && sfp->excpt) {
    exception = TRUE;
  }

  if (adjacent) {
    ctmp = SeqLocPrint (slp);
    /*
    if (exception) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    */
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_AbuttingIntervals, "Adjacent intervals in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }

  if (sfp != NULL) {
    strand1 = SeqLocStrand (sfp->location);
    strand2 = SeqLocStrand (slp);
    if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadAnticodonStrand, "Anticodon should be on minus strand");
    } else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadAnticodonStrand, "Anticodon should be on plus strand");
    }
  }

  if (exception) {
    /* trans splicing exception turns off both mixed_strand and out_of_order messages */
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      return;
    }
  }

  if (mixed_strand || unmarked_strand || (!ordered)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    if (mixed_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed strands in Anticodon [%s]", ctmp);
    } else if (unmarked_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed plus and unknown strands in Anticodon [%s]", ctmp);
    }
    if (!ordered)
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqLocOrder, "Intervals out of order in Anticodon [%s]", ctmp);
    MemFree (ctmp);
    return;
  }

  /* newer check for intervals out of order on segmented bioseq */

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;

  if (SeqLocBadSortOrder (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqLocOrder, "Intervals out of order in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }

  /* newer check for mixed strand on segmented bioseq */

  if (SeqLocMixedStrands (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "Mixed strands in Anticodon [%s]", ctmp);
    MemFree (ctmp);
  }
}

static Boolean JustQuotes (CharPtr str)

{
  Char  ch;

  if (str == NULL) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch != '"' && ch != ' ') return FALSE;
    str++;
    ch = *str;
  }

  return TRUE;
}

typedef struct dummysmfedata {
  Int4        max;
  Int4        num_at_max;
  Boolean     equivalent_genes;
  GeneRefPtr  grp_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK DummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  GeneRefPtr    grp, grpx;
  Int4          len;
  Boolean       redundantgenexref;
  CharPtr       syn1, syn2;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
    dsp->equivalent_genes = FALSE;
    dsp->grp_at_max = grp;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
    grpx = dsp->grp_at_max;
    if (grpx != NULL) {
      redundantgenexref = FALSE;
      if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
        if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
          redundantgenexref = TRUE;
        }
      } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
        if (StringICmp (grp->locus, grpx->locus) == 0) {
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
    }
    if (redundantgenexref) {
      dsp->equivalent_genes = TRUE;
    }
  }

  return TRUE;
}

typedef struct govstruc {
  CharPtr  term;
  CharPtr  goid;
  CharPtr  evidence;
  Int4     pmid;
} GovStruc, PNTR GovStrucPtr;

static int LIBCALLBACK SortVnpByGvsp (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GovStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GovStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GovStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}

static void ValidateGoTerms (
  ValidStructPtr vsp,
  GatherContextPtr gcp,
  SeqFeatPtr sfp,
  UserFieldPtr entryhead,
  CharPtr qualtype
)

{
  UserFieldPtr  entry, topufp, ufp;
  CharPtr       evidence, goid, textstr;
  Char          gid [32];
  GovStrucPtr   gsp, lastgsp;
  ValNodePtr    head = NULL, vnp;
  Int2          j;
  ObjectIdPtr   oip;
  Int4          pmid;

  if (entryhead == NULL) return;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    textstr = NULL;
    evidence = NULL;
    goid = NULL;
    pmid = 0;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; goFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, goFieldType [j]) == 0) break;
      }
      if (goFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }

    if (StringDoesHaveText (textstr)) {
      gsp = (GovStrucPtr) MemNew (sizeof (GovStruc));
      if (gsp != NULL) {
        gsp->term = StringSave (textstr);
        gsp->goid = StringSave (goid);
        gsp->evidence = StringSave (evidence);
        gsp->pmid = pmid;
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL) return;
  head = ValNodeSort (head, SortVnpByGvsp);

  lastgsp = NULL;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GovStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL) continue;
    if (StringHasNoText (gsp->goid)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneOntologyTermMissingGOID, "GO term does not have GO identifier");
    }
    if (lastgsp != NULL) {
      if (StringICmp (gsp->term, lastgsp->term) == 0 || StringICmp (gsp->goid, lastgsp->goid) == 0) {
        if (gsp->pmid == lastgsp->pmid && StringICmp (gsp->evidence, lastgsp->evidence) == 0) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_DuplicateGeneOntologyTerm, "Duplicate GO term on feature");
        }
      }
    }
    lastgsp = gsp;
  }

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GovStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL) continue;
    gsp->term = MemFree (gsp->term);
    gsp->goid = MemFree (gsp->goid);
    gsp->evidence = MemFree (gsp->evidence);
  }
  ValNodeFreeData (head);
}

static void ValidateGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr      entry;
  GatherContextPtr  gcp;
  Int2              i;
  ObjectIdPtr       oip;
  ValidStructPtr    vsp;

  gcp = (GatherContextPtr) userdata;
  if (gcp == NULL) return;
  vsp = (ValidStructPtr) (gcp->userdata);
  if (vsp == NULL) return;

  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; goQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, goQualType [i]) == 0) break;
  }
  if (goQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  ValidateGoTerms (vsp, gcp, vsp->sfp, entry, goQualType [i]);
}

static void ValidateGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, ValidateGoTermsUfp);
  }
}

static void LookForAccnLocs (SeqIdPtr sip, Pointer userdata)

{
  BoolPtr       bp;
  TextSeqIdPtr  tsip;

  if (sip == NULL || userdata == NULL) return;
  bp = (BoolPtr) userdata;

  switch (sip->choice) {
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_OTHER :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringDoesHaveText (tsip->accession)) {
          *bp = TRUE;
        }
      }
      break;
    default :
      break;
  }
}

static Boolean LocationIsFar (SeqLocPtr location)

{
  BioseqPtr    bsp;
  DeltaSeqPtr  dsp;
  Boolean      is_far = FALSE;
  SeqLocPtr    loc;
  SeqEntryPtr  oldscope;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  if (location == NULL) return FALSE;

  oldscope = SeqEntrySetScope (NULL);

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice != SEQLOC_NULL) {
      sip = SeqLocId (slp);
      bsp = BioseqFind (sip);
      if (bsp == NULL) {
        is_far = TRUE;
      } else if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
        for (dsp = (DeltaSeqPtr) bsp->seq_ext;
             dsp != NULL && (! is_far);
             dsp = dsp->next) {
          if (dsp->choice != 1) continue;
          loc = (SeqLocPtr) dsp->data.ptrvalue;
          if (loc == NULL) continue;
          if (loc->choice == SEQLOC_NULL) continue;
          sip = SeqLocId (loc);
          bsp = BioseqFind (sip);
          if (bsp == NULL) {
            is_far = TRUE;
          }
        }
      } else if (bsp->repr == Seq_repr_seg && bsp->seq_ext_type == 1) {
        for (loc = (SeqLocPtr) bsp->seq_ext;
             loc != NULL && (! is_far);
             loc = loc->next) {
          if (loc == NULL) continue;
          if (loc->choice == SEQLOC_NULL) continue;
          sip = SeqLocId (loc);
          bsp = BioseqFind (sip);
          if (bsp == NULL) {
            is_far = TRUE;
          }
        }
      }
    }
    slp = SeqLocFindNext (location, slp);
  }

  SeqEntrySetScope (oldscope);

  return is_far;
}

static Boolean NoFetchFunctions (void)

{
  ObjMgrProcPtr  ompp = NULL;

  ompp = ObjMgrProcFindNext (NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_BIOSEQ, NULL);

  return (Boolean) (ompp == NULL);
}

static CharPtr infMessage [] = {
  "unknown error",
  "empty inference string",
  "bad inference prefix",
  "bad inference body",
  "single inference field",
  "spaces in inference",
  "same species misused",
  "bad inference accession",
  "bad inference accession version",
  "accession.version not public",
  "bad accession type",
  NULL
};

static CharPtr rnaNameByType [] = {
  "unknown",
  "premsg",
  "mRNA",
  "tRNA",
  "rRNA",
  "snRNA",
  "scRNA",
  "snoRNA",
  "otherRNA",
  NULL
};

static Boolean ValStrandsMatch (Uint1 featstrand, Uint1 locstrand)

{
  if (featstrand == locstrand) return TRUE;
  if (locstrand == Seq_strand_unknown && featstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_unknown && locstrand != Seq_strand_minus) return TRUE;
  if (featstrand == Seq_strand_both && locstrand != Seq_strand_minus) return TRUE;
  if (locstrand == Seq_strand_both) return TRUE;
  return FALSE;
}

static CharPtr badGeneSyn [] = {
  "alpha",
  "alternative",
  "beta",
  "cellular",
  "cytokine",
  "delta",
  "drosophila",
  "epsilon",
  "gamma",
  "HLA",
  "homolog",
  "mouse",
  "orf",
  "partial",
  "plasma",
  "precursor",
  "pseudogene",
  "putative",
  "rearranged",
  "small",
  "trna",
  "unknown function",
  "unknown protein",
  "unknown",
  "unnamed",
  NULL
};

static CharPtr badProtName [] = {
  "'hypothetical protein",
  "alpha",
  "alternative",
  "alternatively spliced",
  "bacteriophage hypothetical protein",
  "beta",
  "cellular",
  "cnserved hypothetical protein",
  "conesrved hypothetical protein",
  "conserevd hypothetical protein",
  "conserved archaeal protein",
  "conserved domain protein",
  "conserved hypohetical protein",
  "conserved hypotehtical protein",
  "conserved hypotheical protein",
  "conserved hypothertical protein",
  "conserved hypothetcial protein",
  "conserved hypothetical exported protein",
  "conserved hypothetical integral membrane protein",
  "conserved hypothetical membrane protein",
  "conserved hypothetical phage protein",
  "conserved hypothetical prophage protein",
  "conserved hypothetical protein - phage associated",
  "conserved hypothetical protein fragment 3",
  "conserved hypothetical protein, fragment",
  "conserved hypothetical protein, putative",
  "conserved hypothetical protein, truncated",
  "conserved hypothetical protein, truncation",
  "conserved hypothetical protein; possible membrane protein",
  "conserved hypothetical protein; putative membrane protein",
  "conserved hypothetical protein.",
  "conserved hypothetical protein",
  "conserved hypothetical proteins",
  "conserved hypothetical protien",
  "conserved hypothetical transmembrane protein",
  "conserved hypothetical",
  "conserved hypotheticcal protein",
  "conserved hypthetical protein",
  "conserved in bacteria",
  "conserved membrane protein",
  "conserved protein of unknown function ; putative membrane protein",
  "conserved protein of unknown function",
  "conserved protein",
  "conserved unknown protein",
  "conservedhypothetical protein",
  "conserverd hypothetical protein",
  "conservered hypothetical protein",
  "consrved hypothetical protein",
  "converved hypothetical protein",
  "cytokine",
  "delta",
  "drosophila",
  "duplicated hypothetical protein",
  "epsilon",
  "gamma",
  "HLA",
  "homeodomain protein",
  "homeodomain",
  "homolog",
  "hyopthetical protein",
  "hypotethical",
  "hypotheical protein",
  "hypothertical protein",
  "hypothetcical protein",
  "hypothetical  protein",
  "hypothetical conserved protein",
  "hypothetical exported protein",
  "hypothetical novel protein",
  "hypothetical orf",
  "hypothetical phage protein",
  "hypothetical prophage protein",
  "hypothetical protein - phage associated",
  "hypothetical protein (fragment)",
  "hypothetical protein (multi-domain)",
  "hypothetical protein (phage associated)",
  "hypothetical protein fragment ",
  "hypothetical protein fragment 1",
  "hypothetical protein predicted by genemark",
  "hypothetical protein predicted by glimmer",
  "hypothetical protein predicted by glimmer/critica",
  "hypothetical protein-putative conserved hypothetical protein",
  "hypothetical protein, conserved",
  "hypothetical protein, phage associated",
  "hypothetical protein, truncated",
  "hypothetical protein.",
  "hypothetical protein",
  "hypothetical proteins",
  "hypothetical protien",
  "hypothetical transmembrane protein",
  "hypothetical",
  "hypothetoical protein",
  "hypothteical protein",
  "identified by sequence similarity; putative; ORF located\nusing Blastx/FrameD",
  "identified by sequence similarity; putative; ORF located\nusing Blastx/Glimmer/Genemark",
  "ion channel",
  "membrane protein, putative",
  "mouse",
  "narrowly conserved hypothetical protein ",
  "novel protein",
  "orf, conserved hypothetical protein",
  "orf, hypothetical protein",
  "orf, hypothetical, fragment",
  "orf, hypothetical",
  "orf, partial conserved hypothetical protein",
  "orf; hypothetical protein",
  "orf; unknown function",
  "orf",
  "partial cds, hypothetical",
  "partial",
  "partially conserved hypothetical protein",
  "phage hypothetical protein",
  "phage-related conserved hypothetical protein",
  "phage-related protein",
  "plasma",
  "possible hypothetical protein",
  "precursor",
  "predicted coding region",
  "predicted protein (pseudogene)",
  "predicted protein family",
  "predicted protein",
  "product uncharacterised protein family",
  "protein family",
  "protein of unknown function",
  "pseudogene",
  "putative conserved protein",
  "putative exported protein",
  "putative hypothetical protein",
  "putative membrane protein",
  "putative orf; unknown function",
  "putative phage protein",
  "putative protein",
  "putative",
  "putative",
  "rearranged",
  "repeats containing protein",
  "reserved",
  "ribosomal protein",
  "similar to",
  "small hypothetical protein",
  "small",
  "transmembrane protein",
  "trna",
  "trp repeat",
  "trp-repeat protein",
  "truncated conserved hypothetical protein",
  "truncated hypothetical protein",
  "uncharacterized conserved membrane protein",
  "uncharacterized conserved protein",
  "uncharacterized conserved secreted protein",
  "uncharacterized protein conserved in archaea",
  "uncharacterized protein conserved in bacteria",
  "uncharacterized protein",
  "unique hypothetical protein",
  "unique hypothetical",
  "unknown CDS",
  "unknown function",
  "unknown gene",
  "unknown protein",
  "unknown protein",
  "unknown-related protein",
  "unknown, conserved protein",
  "unknown, hypothetical",
  "unknown; predicted coding region",
  "unknown",
  "unknown",
  "unnamed protein product",
  "unnamed",
  "very hypothetical protein",
  NULL
};

static Boolean NameInList (CharPtr name, CharPtr PNTR list, size_t numelements)

{
  Int2  L, R, mid;

  if (StringHasNoText (name) || list == NULL || numelements < 1) return FALSE;

  L = 0;
  R = numelements - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (list [mid], name) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (list [R], name) == 0) return TRUE;

  return FALSE;
}

static Boolean HasBadCharacter (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (ch == '?' || ch == '!' || ch == '~') return TRUE;
    str++;
    ch = *str;
  }

  return FALSE;
}

static Boolean EndsWithBadCharacter (CharPtr str)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return FALSE;

  len = StringLen (str);
  if (len < 1) return FALSE;

  ch = str [len - 1];
  if (ch == '_' || ch == '.' || ch == ',' || ch == ':' || ch == ';') return TRUE;

  return FALSE;
}

static Boolean EndsWithHyphen (CharPtr str)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return FALSE;

  len = StringLen (str);
  if (len < 1) return FALSE;

  ch = str [len - 1];
  if (ch == '-') return TRUE;

  return FALSE;
}


static Boolean CouldExtendPartial (SeqLocPtr slp, Boolean partial5)
{
  BioseqPtr bsp;
  Int4      pos;
  Uint1     strand;
  Char      str[4];
  Boolean   rval = FALSE;

  if (slp == NULL) {
    return FALSE;
  }

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) {
    return FALSE;
  }
  strand = SeqLocStrand (slp);

  if ((strand != Seq_strand_minus && partial5) || (strand == Seq_strand_minus && !partial5)) {
    pos = SeqLocStart (slp);
    if (pos < 2) {
      rval = TRUE;
    } else if (bsp->repr == Seq_repr_delta) {
      /* wasn't close to the sequence end, but perhaps it is close to a gap */
      SeqPortStreamInt (bsp, pos - 3, pos - 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
      if (str[0] == '-' || str[1] == '-' || str[2] == '-') {
        rval = TRUE;
      }
    }    
  } else {
    pos = SeqLocStop (slp);
    if (pos > bsp->length - 2) {
      rval = TRUE;
    } else {
      /* wasn't close to the sequence end, but perhaps it is close to a gap */
      SeqPortStreamInt (bsp, pos + 1, pos + 3, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) str, NULL);
      if (str[0] == '-' || str[1] == '-' || str[2] == '-') {
        rval = TRUE;
      }
    }
  }

  return rval;
}


NLM_EXTERN void ValidateSeqFeat (GatherContextPtr gcp)
{
  Int2            type, i, j;
  static char    *parterr[2] = { "PartialProduct", "PartialLocation" };
  static char    *parterrs[4] = {
    "Start does not include first/last residue of sequence",
    "Stop does not include first/last residue of sequence",
    "Internal partial intervals do not include first/last residue of sequence",
    "Improper use of partial (greater than or less than)"
  };
  Uint2           partials[2], errtype;
  Char            buf[80];
  CharPtr         tmp;
  ValidStructPtr  vsp;
  SeqFeatPtr      sfp;
  CdRegionPtr     crp;
  CodeBreakPtr    cbp, prevcbp;
  CharPtr         ctmp;
  RnaRefPtr       rrp;
  tRNAPtr         trp;
  GBQualPtr       gbq;
  Boolean         pseudo, excpt, conflict, codonqual,
                  anticodonqual, productqual, protidqual,
                  transidqual, ovgenepseudo;
  ImpFeatPtr      ifp;
  GeneRefPtr      grp;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  BioseqPtr       bsp;
  BioseqContextPtr bcp = NULL;
  BioSourcePtr    biop, dbiop;
  OrgNamePtr      onp;
  OrgRefPtr       orp, dorp;
  SubSourcePtr    ssp;
  Boolean         transgenic;
  Int2            biopgencode;
  Int2            cdsgencode;
  Boolean         plastid;
  GeneticCodePtr  gc;
  PubdescPtr      pdp;
  /*
  DbtagPtr        db = NULL;
  Int4            id = -1;
  */
  SeqMgrDescContext context;
  GeneRefPtr      grpx;
  SeqFeatPtr      sfpx = NULL, sfpy = NULL, prt;
  SeqFeatPtr      operon;
  Boolean         redundantgenexref;
  SeqMgrFeatContext fcontext;
  CharPtr         syn1, syn2, label = NULL;
  Uint2           oldEntityID;
  Uint4           oldItemID;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  BioseqPtr       protBsp;
  ErrSev          sev;
  Boolean         multitoken;
  Char            ch;
  CharPtr         ptr;
  Int4            anticodonlen;
  Boolean         badanticodon;
  SeqLocPtr       slp;
  Int2            count;
  DummySmfeData   dsd;
  CharPtr         str;
  size_t          len;
  Boolean         isgap;
  Boolean         badseq;
  Boolean         is_seqloc_bond;
  SeqBondPtr      sbp;
  SeqFeatXrefPtr  xref, matchxref;
  SeqFeatPtr      matchsfp, origsfp;
  Boolean         hasxref;
  CharPtr         sfp_old_locus_tag;
  CharPtr         gene_old_locus_tag;
  Boolean         bypassGeneTest;
  Boolean         dicistronic = FALSE;
  Int2            inferenceCode;
  Boolean         hasInference = FALSE;
  Boolean         hasExperiment = FALSE;
  Boolean         accn_seqid;
  SeqDescrPtr     sdp;
  SeqMgrDescContext dcontext;
  MolInfoPtr      mip;
  Boolean         farFetchProd;
  Boolean         skip;
  Boolean         is_nc = FALSE;
  Boolean         no_nonconsensus_except = TRUE;


  vsp = (ValidStructPtr) (gcp->userdata);
  sfp = (SeqFeatPtr) (gcp->thisitem);
  vsp->descr = NULL;
  vsp->sfp = sfp;
  type = (Int2) (sfp->data.choice);

  ValidateSeqLoc (vsp, sfp->location, "Location");

  ValidateSeqLoc (vsp, sfp->product, "Product");

  CheckForBothOrBothRev (vsp, sfp);

  if (vsp->feat_loc_has_gi) {
    accn_seqid = FALSE;
    VisitSeqIdsInSeqLoc (sfp->location, (Pointer) &accn_seqid, LookForAccnLocs);
    if (accn_seqid) {
      if (! vsp->is_smupd_in_sep) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureRefersToAccession, "Feature location refers to accession");
      }
    }
  }

  if (vsp->feat_prod_has_gi) {
    accn_seqid = FALSE;
    VisitSeqIdsInSeqLoc (sfp->product, (Pointer) &accn_seqid, LookForAccnLocs);
    if (accn_seqid) {
      if (! vsp->is_smupd_in_sep) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_FeatureRefersToAccession, "Feature product refers to accession");
      }
    }
  }

  farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
  partials[0] = SeqLocPartialCheckEx (sfp->product, farFetchProd);
  partials[1] = SeqLocPartialCheck (sfp->location);

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "nonconsensus splice site") != NULL) {
      no_nonconsensus_except = FALSE;
    }
  }

  if ((partials[0] != SLP_COMPLETE) || (partials[1] != SLP_COMPLETE) || (sfp->partial)) {       /* partialness */
    /* a feature on a partial sequence should be partial -- if often isn't */
    if ((!sfp->partial) && (partials[1] != SLP_COMPLETE) && (sfp->location->choice == SEQLOC_WHOLE)) {
      ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem, "On partial Bioseq, SeqFeat.partial should be TRUE");
    }
    /* a partial feature, with complete location, but partial product */
    else if ((sfp->partial) && (sfp->product != NULL) && (partials[1] == SLP_COMPLETE) && (sfp->product->choice == SEQLOC_WHOLE)
             && (partials[0] != SLP_COMPLETE)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "When SeqFeat.product is a partial Bioseq, SeqFeat.location should also be partial");
    }
    /* gene on segmented set is now 'order', should also be partial */
    else if (type == SEQFEAT_GENE && sfp->product == NULL && partials[1] == SLP_INTERNAL) {
      if (!sfp->partial) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "Gene of 'order' with otherwise complete location should have partial flag set");
      }
    }
    /* inconsistent combination of partial/complete product,location,partial flag - part 1 */
    else if (((partials[0] == SLP_COMPLETE) && (sfp->product != NULL))) {
      sev = SEV_WARNING;
      bsp = GetBioseqGivenSeqLoc (sfp->product, gcp->entityID);
      /* if not local bioseq product, lower severity */
      if (bsp == NULL) {
        sev = SEV_INFO;
      }
      tmp = StringMove (buf, "Inconsistent: ");
      if (sfp->product != NULL) {
        tmp = StringMove (tmp, "Product= ");
        if (partials[0])
          tmp = StringMove (tmp, "partial, ");
        else
          tmp = StringMove (tmp, "complete, ");
      }
      tmp = StringMove (tmp, "Location= ");
      if (partials[1])
        tmp = StringMove (tmp, "partial, ");
      else
        tmp = StringMove (tmp, "complete, ");
      tmp = StringMove (tmp, "Feature.partial= ");
      if (sfp->partial)
        tmp = StringMove (tmp, "TRUE");
      else
        tmp = StringMove (tmp, "FALSE");
      if (bsp == NULL && LocationIsFar (sfp->product) && NoFetchFunctions ()) {
        vsp->far_fetch_failure = TRUE;
      } else {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PartialsInconsistent, buf);
      }
    /* inconsistent combination of partial/complete product,location,partial flag - part 2 */
    } else if ((partials[1] == SLP_COMPLETE) || (!sfp->partial)) {
      tmp = StringMove (buf, "Inconsistent: ");
      if (sfp->product != NULL) {
        tmp = StringMove (tmp, "Product= ");
        if (partials[0])
          tmp = StringMove (tmp, "partial, ");
        else
          tmp = StringMove (tmp, "complete, ");
      }
      tmp = StringMove (tmp, "Location= ");
      if (partials[1])
        tmp = StringMove (tmp, "partial, ");
      else
        tmp = StringMove (tmp, "complete, ");
      tmp = StringMove (tmp, "Feature.partial= ");
      if (sfp->partial)
        tmp = StringMove (tmp, "TRUE");
      else
        tmp = StringMove (tmp, "FALSE");
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialsInconsistent, buf);
    }
    /* 5' or 3' partial location giving unclassified partial product */
    else if (((partials [1] & SLP_START) != 0 || ((partials [1] & SLP_STOP) != 0)) && ((partials [0] & SLP_OTHER) != 0) && sfp->partial) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "5' or 3' partial location should not have unclassified partial in product molinfo descriptor");
    }

    /* may have other error bits set as well */
    for (i = 0; i < 2; i++) {
      errtype = SLP_NOSTART;
      for (j = 0; j < 4; j++) {
        bypassGeneTest = FALSE;
        if (partials[i] & errtype) {
          if (i == 1 && j < 2 && IsCddFeat (sfp)) {
            /* suppresses  warning */
          } else if (i == 1 && j < 2 && sfp->data.choice == SEQFEAT_GENE && SameAsCDS (sfp, errtype, NULL)) {
            /*
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
              "%s: %s",
              parterr[i], parterrs[j]);
            */
          } else if (i == 1 && j < 2 && sfp->data.choice == SEQFEAT_GENE && SameAsMRNA (sfp, errtype)) {
          } else if (i == 1 && j < 2 && sfp->idx.subtype == FEATDEF_mRNA && SameAsCDS (sfp, errtype, &bypassGeneTest)) {
          } else if (i == 1 && j < 2 && sfp->idx.subtype == FEATDEF_mRNA && (! bypassGeneTest) && SameAsGene (sfp)) {

          } else if (i == 1 && j < 2 && sfp->idx.subtype == FEATDEF_exon && SameAsMRNA (sfp, errtype)) {

          } else if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
            vsp->far_fetch_failure = TRUE;
          
          } else if (i == 1 && j < 2 && sfp->data.choice == SEQFEAT_CDREGION && SameAsMRNA (sfp, errtype) &&
                     PartialAtSpliceSiteOrGap (vsp, sfp->location, errtype, &isgap, &badseq)) {
          } else if (i == 1 && j < 2 && PartialAtSpliceSiteOrGap (vsp, sfp->location, errtype, &isgap, &badseq)) {
            if (! isgap) {
              if (sfp->idx.subtype != FEATDEF_CDS || SplicingNotExpected (sfp)) {
                ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
                          "%s: %s (but is at consensus splice site)",
                          parterr[i], parterrs[j]);
              } else if (sfp->idx.subtype == FEATDEF_CDS) {
                bsp = BioseqFindFromSeqLoc (sfp->location);
                if (bsp != NULL) {
                  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
                  if (sdp != NULL) {
                    mip = (MolInfoPtr) sdp->data.ptrvalue;
                    if (mip != NULL) {
                      if (mip->biomol == MOLECULE_TYPE_MRNA) {
                        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                                  "%s: %s (but is at consensus splice site, but is on an mRNA that is already spliced)",
                                  parterr[i], parterrs[j]);
                      }
                    }
                  }
                }
              }
            }
            if (badseq) {
              ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
                "%s: %s (and is at bad sequence)",
                parterr[i], parterrs[j]);
            }
          } else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->excpt &&
                     StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
          } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 0) {
            if (no_nonconsensus_except) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                "%s: %s", parterr[i], "5' partial is not at start AND"
                " is not at consensus splice site");
            }
          } else if (sfp->data.choice == SEQFEAT_CDREGION && j == 1) {
            if (no_nonconsensus_except) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
                "%s: %s", parterr[i], "3' partial is not at stop AND"
                " is not at consensus splice site");
            }
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
              "%s: %s", parterr[i], parterrs[j]);
          }
        }
        errtype <<= 1;
      }
    }

  }

  CheckForIllegalDbxref (vsp, gcp, sfp, sfp->dbxref);
  /*
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    id = -1;
    db = vnp->data.ptrvalue;
    if (db && db->db) {
      for (i = 0; i < DBNUM; i++) {
        if (StringCmp (db->db, dbtag[i]) == 0) {
          id = i;
          break;
        }
      }
      if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
      }
    }
  }
  */

  switch (type) {
  case 1:                      /* Gene-ref */
    grp = (GeneRefPtr) (sfp->data.value.ptrvalue);
    if (grp != NULL) {
      if (EmptyOrNullString (grp->locus) &&
          EmptyOrNullString (grp->allele) && EmptyOrNullString (grp->desc) &&
          EmptyOrNullString (grp->maploc) && EmptyOrNullString (grp->locus_tag) &&
          grp->db == NULL && grp->syn == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneRefHasNoData, "There is a gene feature where all fields are empty");
      }
      if (! StringHasNoText (grp->locus_tag)) {
        multitoken = FALSE;
        for (ptr = grp->locus_tag, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
          if (IS_WHITESP (ch)) {
            multitoken = TRUE;
          }
        }
        if (multitoken) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusTagProblem, "Gene locus_tag '%s' should be a single word without any spaces", grp->locus_tag);
        }
        /* check for matching old_locus_tag */
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringCmp (grp->locus_tag, gbq->val) == 0) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_LocusTagProblem, "Gene locus_tag and old_locus_tag '%s' match", grp->locus_tag);
          }
        }        
      }
      CheckForIllegalDbxref (vsp, gcp, sfp, grp->db);
      if (StringDoesHaveText (grp->allele)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "allele") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->allele) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Redundant allele qualifier (%s) on gene", gbq->val);
            } else if (sfp->idx.subtype != FEATDEF_variation) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "Hidden allele qualifier (%s) on gene", gbq->val);
            }
          }
        }
      }
      /*
      for (vnp = grp->db; vnp != NULL; vnp = vnp->next) {
        id = -1;
        db = vnp->data.ptrvalue;
        if (db && db->db) {
          for (i = 0; i < DBNUM; i++) {
            if (StringCmp (db->db, dbtag[i]) == 0) {
              id = i;
              break;
            }
          }
          if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
          }
        }
      }
      */
      if (grp->locus != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as gene locus");
        }
      }
      if (grp->locus != NULL) {
        if (HasBadCharacter (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "Gene locus contains undesired character");
        }
        if (EndsWithBadCharacter (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "Gene locus ends with undesired character");
        }
        if (EndsWithHyphen (grp->locus)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "Gene locus ends with hyphen");
        }
      }
      if (grp->locus_tag != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus_tag, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as gene locus_tag");
        }
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->locus_tag) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "old_locus_tag has same value as gene locus_tag");
            }
          }
        }
      }
      if (grp->syn != NULL && (vsp->is_refseq_in_sep || vsp->seqSubmitParent)) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (NameInList (str, badGeneSyn, sizeof (badGeneSyn) / sizeof (badGeneSyn [0]))) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "Uninformative gene synonym '%s'", str);
          }
          if (StringDoesHaveText (grp->locus) && StringCmp (grp->locus, str) == 0) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene synonym has same value as gene locus");
          }
        }
      }
      if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grp->desc) && StringCmp (grp->locus, grp->desc) == 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene description has same value as gene locus");
      }
      if (StringHasNoText (grp->locus) && StringHasNoText (grp->desc) && grp->syn != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredGeneSynonym, "gene synonym without gene locus or description");
      }
      /* - need to ignore if curated drosophila - add to vsp internal flags for efficiency?
      if (StringDoesHaveText (grp->locus)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->locus) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "old_locus_tag has same value as gene locus");
            }
          }
        }
      }
      */
    }
    break;
  case 2:                      /* Org-ref */
    break;
  case 3:                      /* Cdregion */
    pseudo = sfp->pseudo;       /* now also uses new feature pseudo flag */
    excpt = FALSE;
    conflict = FALSE;
    codonqual = FALSE;
    crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
    if (crp != NULL) {
      conflict = crp->conflict;
    }
    protidqual = FALSE;
    transidqual = FALSE;
    ovgenepseudo = FALSE;
    gbq = sfp->qual;
    while (gbq != NULL) {
      if (StringICmp (gbq->qual, "pseudo") == 0) {
        pseudo = TRUE;
      }
      if (StringICmp (gbq->qual, "exception") == 0) {
        excpt = TRUE;
      }
      if (StringICmp (gbq->qual, "codon") == 0) {
        codonqual = TRUE;
      }
      if (StringICmp (gbq->qual, "protein_id") == 0) {
        protidqual = TRUE;
      }
      if (StringICmp (gbq->qual, "transcript_id") == 0) {
        transidqual = TRUE;
      }
      gbq = gbq->next;
    }
    if (protidqual) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "protein_id should not be a gbqual on a CDS feature");
    }
    if (transidqual) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "transcript_id should not be a gbqual on a CDS feature");
    }
    if (OverlappingGeneIsPseudo (sfp)) {
      pseudo = TRUE;
      ovgenepseudo = TRUE;
    }
    if ((!pseudo) && (!conflict)) {
      CdTransCheck (vsp, sfp);
      SpliceCheck (vsp, sfp);
    } else if (conflict) {
      CdConflictCheck (vsp, sfp);
    }
    CdsProductIdCheck (vsp, sfp);
    crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
    if (crp != NULL) {
      if (crp->code_break != NULL && StringISearch (sfp->except_text, "RNA editing") != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExceptAndRnaEditing, "CDS has both RNA editing /exception and /transl_except qualifiers");
      }
      prevcbp = NULL;
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        i = SeqLocCompare (cbp->loc, sfp->location);
        if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B)) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Code-break location not in coding region");
        }
        if (prevcbp != NULL) {
          i = SeqLocCompare (cbp->loc, prevcbp->loc);
          if (i == SLC_A_EQ_B) {
            ctmp = SeqLocPrint (cbp->loc);
            if (ctmp != NULL) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateTranslExcept, "Multiple code-breaks at same location [%s]", ctmp);
            } else {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateTranslExcept, "Multiple code-breaks at same location");
            }
            MemFree (ctmp);
          }
        }
        prevcbp = cbp;
      }
      if (excpt && (!sfp->excpt)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception flag should be set in coding region");
      }
      if (crp->orf && sfp->product != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OrfCdsHasProduct, "An ORF coding region should not have a product");
      }
      if (pseudo && sfp->product != NULL) {
        if (ovgenepseudo) {
          if (sfp->pseudo) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsHasProduct, "A pseudo coding region should not have a product");
          } else {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsViaGeneHasProduct, "A coding region overlapped by a pseudogene should not have a product");
          }
        } else {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PseudoCdsHasProduct, "A pseudo coding region should not have a product");
        }
      }
      if (pseudo && SeqMgrGetProtXref (sfp) != NULL) {
        if (NGorNT (vsp->sep, sfp->location, &is_nc) || IsEMBLAccn (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        } else {
          sev = SEV_ERROR;
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_PseudoCdsHasProtXref, "A pseudo coding region should not have a protein xref");
      }
      if (codonqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CodonQualifierUsed, "Use the proper genetic code, if available, or set transl_excepts on specific codons");
      }
      biopgencode = 0;
      cdsgencode = 0;
      bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
      if (bsp != NULL) {
        vnp = NULL;
        if (vsp->useSeqMgrIndexes) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
        } else {
          bcp = BioseqContextNew (bsp);
          vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
        }
        if (vnp != NULL && vnp->data.ptrvalue != NULL) {
          plastid = FALSE;
          biop = (BioSourcePtr) vnp->data.ptrvalue;
          orp = biop->org;
          if (orp != NULL && orp->orgname != NULL) {
            onp = orp->orgname;
            if (biop->genome == GENOME_kinetoplast ||
                biop->genome == GENOME_mitochondrion ||
                biop->genome == GENOME_hydrogenosome) {
              biopgencode = onp->mgcode;
            } else if (biop->genome == GENOME_chloroplast ||
                       biop->genome == GENOME_chromoplast ||
                       biop->genome == GENOME_plastid ||
                       biop->genome == GENOME_cyanelle ||
                       biop->genome == GENOME_apicoplast ||
                       biop->genome == GENOME_leucoplast ||
                       biop->genome == GENOME_proplastid) {
              biopgencode = 11;
              plastid = TRUE;
            } else {
              biopgencode = onp->gcode;
            }
            gc = crp->genetic_code;
            if (gc != NULL) {
              for (vnp = gc->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2) {
                  cdsgencode = (Int2) vnp->data.intvalue;
                }
              }
            }
            if (biopgencode != cdsgencode) {
              if (! vsp->seqSubmitParent) { /* suppress when validator run from tbl2asn */
                if (plastid) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenCodeMismatch,
                            "Genetic code conflict between CDS (code %d) and BioSource.genome biological context (%s) (uses code 11)", (int) cdsgencode, plastidtxt [biop->genome]);
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GenCodeMismatch,
                            "Genetic code conflict between CDS (code %d) and BioSource (code %d)", (int) cdsgencode, (int) biopgencode);
                }
              }
            }
          }
        }
        if (!vsp->useSeqMgrIndexes) {
          BioseqContextFree (bcp);
        }
      }
    }
    /* CheckForBothStrands (vsp, sfp); */
    CheckForBadGeneOverlap (vsp, sfp);
    CheckForBadMRNAOverlap (vsp, sfp);
    CheckForCommonCDSProduct (vsp, sfp);
    CheckCDSPartial (vsp, sfp);
    if (StringDoesHaveText (sfp->comment)) {
      if (LookForECnumberPattern (sfp->comment)) {
        skip = FALSE;
        bsp = BioseqFindFromSeqLoc (sfp->product);
        if (bsp != NULL && ISA_aa (bsp->mol)) {
          prt = SeqMgrGetBestProteinFeature (bsp, NULL);
          if (prt != NULL && prt->data.choice == SEQFEAT_PROT) {
            prp = (ProtRefPtr) prt->data.value.ptrvalue;
            if (prp != NULL) {
              for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (StringHasNoText (str)) continue;
                if (StringStr (sfp->comment, str) != NULL) {
                  skip = TRUE;
                }
                skip = TRUE; /* now suppress even if EC numbers are different */
              }
            }
          }
        }
        if (! skip) {
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in CDS comment");
        }
      }
    }
    break;
  case 4:                      /* Prot-ref */
    prp = (ProtRefPtr) (sfp->data.value.ptrvalue);
    if (prp != NULL) {
      if (prp->processed != 3 && prp->processed != 4) {
        vnp = prp->name;
        if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
            EmptyOrNullString (prp->desc) && prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ProtRefHasNoData, "There is a protein feature where all fields are empty");
        }
        if (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringDoesHaveText (str)) {
            len = StringLen (str);
            if (len > 1) {
              if (str [len - 1] == ']') {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ProteinNameEndsInBracket, "Protein name ends with bracket and may contain organism name");
              }
            }
            if (StringNICmp (str, "hypothetical protein XP_", 24) == 0) {
              bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
              if (bsp != NULL) {
                for (sip = bsp->id; sip != NULL; sip = sip->next) {
                  if (sip->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringICmp (tsip->accession, str + 21) != 0) {
                    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_HpotheticalProteinMismatch, "Hypothetical protein reference does not match accession");
                  }
                }
              }
            }
            if (prp->ec != NULL) {
              if (StringCmp (str, "Hypothetical protein") == 0 ||
                  StringCmp (str, "hypothetical protein") == 0 ||
                  StringCmp (str, "Unknown protein") == 0 ||
                  StringCmp (str, "unknown protein") == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProteinName, "Unknown or hypothetical protein should not have EC number");
              }
            }
            if (LookForECnumberPattern (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in protein title");
            }
          }
          if (str != NULL && sfp->comment != NULL) {
            if (StringCmp (str, sfp->comment) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as protein name");
            }
          }
          if (StringDoesHaveText (sfp->comment)) {
            if (LookForECnumberPattern (sfp->comment)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "Apparent EC number in protein comment");
            }
          }
        }
      }
      CheckForIllegalDbxref (vsp, gcp, sfp, prp->db);
      /*
      for (vnp = prp->db; vnp != NULL; vnp = vnp->next) {
        id = -1;
        db = vnp->data.ptrvalue;
        if (db && db->db) {
          for (i = 0; i < DBNUM; i++) {
            if (StringCmp (db->db, dbtag[i]) == 0) {
              id = i;
              break;
            }
          }
          if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref, "Illegal db_xref type %s", db->db);
          }
        }
      }
      */
      if (prp->name == NULL && prp->processed != 3 && prp->processed != 4) {
        if (StringDoesHaveText (prp->desc)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has description but no name");
        } else if (prp->activity != NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has function but no name");
        } else if (prp->ec != NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has EC number but no name");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NoNameForProtein, "Protein feature has no name");
        }
      }
      if (prp->desc != NULL && sfp->comment != NULL) {
        if (StringCmp (prp->desc, sfp->comment) == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RedundantFields, "Comment has same value as protein description");
        }
      }
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (! ValidateECnumber (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberFormat, "%s is not in proper EC_number format", str);
          } else if (ECnumberNotInList (str)) {
            if (ECnumberWasDeleted (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was deleted", str);
            } else if (ECnumberWasReplaced (str)) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was transferred and is no longer valid", str);
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal value for qualifier EC_number", str);
            }
          }
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "EC number should not be empty");
        }
      }
    }
    if (prp->name != NULL && (vsp->is_refseq_in_sep || vsp->seqSubmitParent)) {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
          if (NameInList (str, badProtName, sizeof (badProtName) / sizeof (badProtName [0]))) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UndesiredProteinName, "Uninformative protein name '%s'", str);
          }
        }
      }
      if (prp->name != NULL) {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "Protein name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "Protein name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "Protein name ends with hyphen");
          }
        }
      }
    break;
  case 5:                      /* RNA-ref */
    rrp = (RnaRefPtr) (sfp->data.value.ptrvalue);

    pseudo = sfp->pseudo;
    ovgenepseudo = FALSE;
    if (OverlappingGeneIsPseudo (sfp)) {
      pseudo = TRUE;
      ovgenepseudo = TRUE;
    }

    if (rrp->type == 2) {       /* mRNA */
      if (!pseudo) {
        MrnaTransCheck (vsp, sfp);      /* transcription check */
        SpliceCheck (vsp, sfp);
      }
      /* CheckForBothStrands (vsp, sfp); */
      CheckForBadGeneOverlap (vsp, sfp);
      CheckForCommonMRNAProduct (vsp, sfp);
      protidqual = FALSE;
      transidqual = FALSE;
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "protein_id") == 0) {
          protidqual = TRUE;
        }
        if (StringICmp (gbq->qual, "transcript_id") == 0) {
          transidqual = TRUE;
        }
        gbq = gbq->next;
      }
      if (protidqual) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "protein_id should not be a gbqual on an mRNA feature");
      }
      if (transidqual) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnFeature, "transcript_id should not be a gbqual on an mRNA feature");
      }
    }
    if (rrp->ext.choice == 2) { /* tRNA */
      trp = (tRNAPtr) (rrp->ext.value.ptrvalue);
      if (trp->anticodon != NULL) {
        badanticodon = FALSE;
        anticodonlen = 0;
        slp = SeqLocFindNext (trp->anticodon, NULL);
        while (slp != NULL) {
          anticodonlen += SeqLocLen (slp);
          i = SeqLocCompare (slp, sfp->location);
          if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B)) {
            badanticodon = TRUE;
          }
          slp = SeqLocFindNext (trp->anticodon, slp);
        }
        if (badanticodon) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Anticodon location not in tRNA");
        }
        if (anticodonlen != 3) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range, "Anticodon is not 3 bases in length");
        }
        ValidateAnticodon (vsp, trp->anticodon);
      }
      CheckTrnaCodons (vsp, gcp, sfp, trp);
    }
    if (rrp->type == 3) {       /* tRNA */
      anticodonqual = FALSE;
      productqual = FALSE;
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "anticodon") == 0) {
          anticodonqual = TRUE;
        } else if (StringICmp (gbq->qual, "product") == 0) {
          productqual = TRUE;
        }
        gbq = gbq->next;
      }
      if (anticodonqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed anticodon qualifier in tRNA");
      }
      if (productqual) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed product qualifier in tRNA");
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 1) { /* tRNA with string extension */
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Unparsed product qualifier in tRNA");
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) { /* tRNA with no extension */
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingTrnaAA, "Missing encoded amino acid qualifier in tRNA");
    }
    if (rrp->type == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_RNAtype0, "RNA type 0 (unknown) not supported");
    }
    if (rrp->type == 4 || rrp->type == 5 || rrp->type == 6 || rrp->type == 7) { /* rRNA, snRNA, scRNA, snoRNA */
      if (rrp->ext.choice != 1 || StringHasNoText ((CharPtr) rrp->ext.value.ptrvalue)) {
        if (! pseudo) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "%s has no name", rnaNameByType [(int) rrp->type]);
        }
      }
    }
    /*
    if (rrp->type == 255 && rrp->ext.choice == 1) {
      str = (CharPtr) rrp->ext.value.ptrvalue;
      if (StringICmp (str, "ncRNA") == 0) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "ncRNA_class") != 0) continue;
          if (StringHasNoText (gbq->val)) continue;
          if (IsStringInNcRNAClassList (gbq->val)) continue;
          ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_InvalidQualifierValue, "Illegal ncRNA_class value '%s'", gbq->val);
        }
      }
    }
    */
    if (rrp->type == 2) {
      if (rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "mRNA name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "mRNA name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "mRNA name ends with hyphen");
          }
        }
      }
    }
    if (rrp->type == 4) {
      if (rrp->ext.choice == 1) {
        str = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringDoesHaveText (str)) {
          if (HasBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadInternalCharacter, "rRNA name contains undesired character");
          }
          if (EndsWithBadCharacter (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingCharacter, "rRNA name ends with undesired character");
          }
          if (EndsWithHyphen (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadTrailingHyphen, "rRNA name ends with hyphen");
          }
        }
      }
    }
    if (sfp->product != NULL) {
      CheckRnaProductType (vsp, gcp, sfp, rrp);
    }

    if (pseudo && sfp->product != NULL && StringISearch (sfp->except_text, "transcribed pseudogene") == NULL) {
      if (ovgenepseudo) {
        if (sfp->pseudo) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaHasProduct, "A pseudo RNA should not have a product");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaViaGeneHasProduct, "An RNA overlapped by a pseudogene should not have a product");
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_PseudoRnaHasProduct, "A pseudo RNA should not have a product");
      }
    }
    break;
  case 6:                      /* Pub */
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    /*
       ValidatePubdesc (vsp, pdp);
     */
    break;
  case 7:                      /* Seq */
    break;
  case 8:                      /* Imp-feat */
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (vsp->validateExons) {

      if (ifp != NULL && StringICmp (ifp->key, "exon") == 0 && (! sfp->pseudo)) {
        SpliceCheckEx (vsp, sfp, TRUE);
      }
    }
    if (ifp != NULL) {
      ValidateImpFeat (vsp, gcp, sfp, ifp);
    }
    break;
  case 9:                      /* Region */
    break;
  case 10:                     /* Comment */
    break;
  case 11:                     /* Bond */
    break;
  case 12:                     /* Site */
    break;
  case 13:                     /* Rsite-ref */
    break;
  case 14:                     /* User-object */
    break;
  case 15:                     /* TxInit */
    break;
  case 16:                     /* Numbering */
    break;
  case 17:                     /* Secondary Structure */
    break;
  case 18:                     /* NonStdRes */
    break;
  case 19:                     /* Heterogen */
    break;
  case 20:                     /* BioSource */
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    if (biop != NULL && biop->is_focus) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_FocusOnBioSourceFeature, "Focus must be on BioSource descriptor, not BioSource feature.");
    }
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
        if (bsp != NULL) {
          vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
          if (vnp != NULL) {
            dbiop = (BioSourcePtr) vnp->data.ptrvalue;
            if (dbiop != NULL) {
              dorp = dbiop->org;
              if (dorp != NULL) {
                if (!StringHasNoText (orp->taxname)) {
                  if (StringICmp (orp->taxname, dorp->taxname) != 0) {
                    if (!dbiop->is_focus) {
                      transgenic = FALSE;
                      for (ssp = dbiop->subtype; ssp != NULL; ssp = ssp->next) {
                        if (ssp->subtype == SUBSRC_transgenic) {
                          transgenic = TRUE;
                        }
                      }
                      if (! transgenic) {
                        oldEntityID = gcp->entityID;
                        oldItemID = gcp->itemID;

                        gcp->entityID = context.entityID;
                        gcp->itemID = context.itemID;
                        gcp->thistype = OBJ_SEQDESC;

                        ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_BioSourceNeedsFocus,
                                  "BioSource descriptor must have focus or transgenic when BioSource feature with different taxname is present.");

                        gcp->entityID = oldEntityID;
                        gcp->itemID = oldItemID;
                        gcp->thistype = OBJ_SEQFEAT;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    /*
       ValidateBioSource (vsp, gcp, biop, sfp, NULL);
     */
    break;
  default:
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidType, "Invalid SeqFeat type [%d]", (int) (type));
    break;
  }
  if (type == SEQFEAT_HET) {
    /* heterogen can have mix of bonds with just "a" point specified */
    is_seqloc_bond = FALSE;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_BOND) {
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          if (sbp->a == NULL || sbp->b != NULL) {
            is_seqloc_bond = TRUE;
          }
        }
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (is_seqloc_bond) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ImproperBondLocation, "Bond location should only be on bond features");
    }
  } else if (type != SEQFEAT_BOND) {
    is_seqloc_bond = FALSE;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_BOND) {
        is_seqloc_bond = TRUE;
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (is_seqloc_bond) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ImproperBondLocation, "Bond location should only be on bond features");
    }
  }
  if (type != 8) {
    ValidateNonImpFeat (vsp, gcp, sfp);
  }
  if ((! sfp->excpt) && (! StringHasNoText (sfp->except_text))) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception text is present, but exception flag is not set");
  }
  if ((sfp->excpt) && (StringHasNoText (sfp->except_text))) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent, "Exception flag is set, but exception text is empty");
  }
  if (! StringHasNoText (sfp->except_text)) {
    ValidateExceptText (vsp, gcp, sfp);
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0 && xref->data.choice == 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "SeqFeatXref with no id or data field");
    } else if (xref->id.choice != 0) {
      matchsfp = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (matchsfp == NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefFeatureMissing, "Cross-referenced feature cannot be found");
      } else {
        hasxref = FALSE;
        for (matchxref = matchsfp->xref; matchxref != NULL; matchxref = matchxref->next) {
          if (matchxref->id.choice != 0) {
            hasxref = TRUE;
            origsfp = SeqMgrGetFeatureByFeatID (matchsfp->idx.entityID, NULL, NULL, matchxref, NULL);
            if (origsfp != sfp) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefNotReciprocal, "Cross-referenced feature does not link reciprocally");
            } else {
              if (sfp->idx.subtype == FEATDEF_CDS && matchsfp->idx.subtype == FEATDEF_mRNA) {
                /* okay */
              } else if (sfp->idx.subtype == FEATDEF_mRNA && matchsfp->idx.subtype == FEATDEF_CDS) {
                /* okay */
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "Cross-references are not between CDS and mRNA pair");
              }
            }
          }
        }
        if (! hasxref) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SeqFeatXrefProblem, "Cross-referenced feature does not have its own cross-reference");
        }
      }
    }
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    /* first check for anything other than replace */
    if (StringICmp (gbq->qual, "replace") != 0) {
      if (JustQuotes (gbq->val)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Qualifier other than replace has just quotation marks");
      }
    }
    /* now check specific gbual types */
    if (StringICmp (gbq->qual, "inference") == 0) {
      hasInference = TRUE;
      inferenceCode = ValidateInferenceQualifier (gbq->val, vsp->inferenceAccnCheck);
      if (inferenceCode != VALID_INFERENCE) {
        if (inferenceCode < VALID_INFERENCE || inferenceCode > BAD_ACCESSION_TYPE) {
          inferenceCode = VALID_INFERENCE;
        }
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidInferenceValue, "Inference qualifier problem - %s (%s)",
                  infMessage [(int) inferenceCode], (gbq->val != NULL)? gbq->val : "?");
      }
    } else if (StringICmp (gbq->qual, "experiment") == 0) {
      hasExperiment = TRUE;
    } else if (StringICmp (gbq->qual, "EC_number") == 0) {
      str = gbq->val;
      if (StringDoesHaveText (str)) {
        if (! ValidateECnumber (str)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberFormat, "%s is not in proper EC_number format", str);
        } else if (ECnumberNotInList (str)) {
          if (ECnumberWasDeleted (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was deleted", str);
          } else if (ECnumberWasReplaced (str)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "EC_number %s was replaced", str);
          } else {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadEcNumberValue, "%s is not a legal value for qualifier EC_number", str);
          }
        }
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_EcNumberProblem, "EC number should not be empty");
      }
    } else if (StringICmp (gbq->qual, "old_locus_tag") == 0) {
      if (StringChr (gbq->val, ',') != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_LocusTagProblem,
                  "old_locus_tag has comma, may contain multiple values");
      }
    }
  }
  if (sfp->exp_ev > 0 && (! hasInference) && (! hasExperiment) && (! vsp->feat_loc_has_gi)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidInferenceValue,
              "Inference or experiment qualifier missing but obsolete experimental evidence qualifier set");
  }

  if (sfp->product != NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL) {
      switch (sip->choice) {
        case SEQID_LOCAL :
      break;
        case SEQID_GENBANK :
        case SEQID_EMBL :
        case SEQID_DDBJ :
        case SEQID_OTHER :
        case SEQID_TPG :
        case SEQID_TPE :
        case SEQID_TPD :
        case SEQID_GPIPE :
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
              if (ValidateAccn (tsip->name) == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                         "Feature product should not put an accession in the Textseq-id 'name' slot");
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                         "Feature product should not use Textseq-id 'name' slot");
              }
            }
          }
          break;
        default :
          break;
      }
    }
    bsp = BioseqFindFromSeqLoc (sfp->location);
    protBsp = BioseqFindFromSeqLoc (sfp->product);
    if (bsp != NULL && protBsp != NULL) {
      if (bsp == protBsp) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_SelfReferentialProduct, "Self-referential feature product");
      }
    }
    if (protBsp != NULL && protBsp->id != NULL) {
      for (sip = protBsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_OTHER :
          case SEQID_TPG :
          case SEQID_TPE :
          case SEQID_TPD :
          case SEQID_GPIPE:
            tsip = (TextSeqIdPtr) sip->data.ptrvalue;
            if (tsip != NULL) {
              if (tsip->accession == NULL && (! StringHasNoText (tsip->name))) {
                if (ValidateAccn (tsip->name) == 0) {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                            "Protein bioseq has Textseq-id 'name' that looks"
                            " like it is derived from a nucleotide accession");
                } else {
                  ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_BadProductSeqId,
                            "Protein bioseq has Textseq-id 'name' and no accession");
                }
              }
            }
            break;
          default :
            break;
        }
      }
    }
    
  }

  if (sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, (Pointer) gcp, ValidateGoTermsSfp);
  }

  if (type != SEQFEAT_GENE) {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (sfpx != NULL) {
        grp = (GeneRefPtr) sfpx->data.value.ptrvalue;
      }
    }
    if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
      if (! StringHasNoText (grp->allele)) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringCmp (gbq->qual, "allele") == 0 && StringDoesHaveText (gbq->val)) {
            if (StringICmp (gbq->val, grp->allele) == 0) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Redundant allele qualifier (%s) on gene and feature", gbq->val);
            } else if (sfp->idx.subtype != FEATDEF_variation) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Mismatched allele qualifier on gene (%s) and feature (%s)", grp->allele, gbq->val);
            }
          }
        }
      }
    }
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) return;

    if (grp == NULL) {
      sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
      if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
      sfp_old_locus_tag = NULL;
      gene_old_locus_tag = NULL;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
          sfp_old_locus_tag = gbq->val;
        }
      }
      for (gbq = sfpx->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "old_locus_tag") == 0 && StringDoesHaveText (gbq->val)) {
          gene_old_locus_tag = gbq->val;
        }
      }
      if (StringDoesHaveText (sfp_old_locus_tag) && StringDoesHaveText (gene_old_locus_tag)) {
        if (StringICmp (sfp_old_locus_tag, gene_old_locus_tag) != 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_OldLocusTagMismtach,
                    "Old locus tag on feature (%s) does not match that on gene (%s)",
                    sfp_old_locus_tag, gene_old_locus_tag);
        }
      }
      MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
      dsd.max = INT4_MAX;
      dsd.num_at_max = 0;
      dsd.equivalent_genes = FALSE;
      dsd.grp_at_max = NULL;
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                               LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);
      if (dsd.num_at_max > 1) {
        if (dsd.equivalent_genes) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefNeeded,
                    "Feature overlapped by %d identical-length equivalent genes but has no cross-reference", (int) dsd.num_at_max);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingGeneXref,
                    "Feature overlapped by %d identical-length genes but has no cross-reference", (int) dsd.num_at_max);
        }
      }
      return;
    }

    if (StringDoesHaveText (grp->locus) /* && sfp->idx.subtype != FEATDEF_tRNA */) {
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sfpx = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
        if (sfpx != NULL) {
          sfpy = sfpx;
        }
        if (sfpx == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutGene,
                    "Feature has gene locus cross-reference but no equivalent gene feature exists");
        } else if (StringStr (sfpx->except_text, "dicistronic gene") != NULL) {
          dicistronic = TRUE;
        }
      }
    }
    if (StringDoesHaveText (grp->locus_tag)) {
      bsp = BioseqFindFromSeqLoc (sfp->location);
      if (bsp != NULL) {
        sfpx = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &fcontext);
        if (sfpx != NULL) {
          sfpy = sfpx;
        }
        if (sfpx == NULL) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutGene,
                    "Feature has gene locus_tag cross-reference but no equivalent gene feature exists");
        } else if (StringStr (sfpx->except_text, "dicistronic gene") != NULL) {
          dicistronic = TRUE;
        }
        /* look for gene xrefs with locus_tag but no locus */
        if (StringHasNoText (grp->locus)
            && sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE 
            && sfpx->data.value.ptrvalue != NULL
            && StringDoesHaveText (((GeneRefPtr)sfpx->data.value.ptrvalue)->locus)) {
            ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneXrefWithoutLocus,
                      "Feature has Gene Xref with locus_tag but no locus, gene with locus_tag and locus exists");
        }
      }
    }

    sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
    if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE)
      return;
    grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
    if (grpx == NULL)
      return;
    redundantgenexref = FALSE;
    label = fcontext.label;
    if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grp->locus_tag)) {
      if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
        redundantgenexref = TRUE;
        label = grp->locus_tag;
      }
    } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grp->locus)) {
      if (StringICmp (grp->locus, grpx->locus) == 0) {
        redundantgenexref = TRUE;
        label = grp->locus;
      }
    } else if (grp->syn != NULL && grpx->syn != NULL) {
      syn1 = (CharPtr) grp->syn->data.ptrvalue;
      syn2 = (CharPtr) grpx->syn->data.ptrvalue;
      if ((StringDoesHaveText (syn1)) && StringDoesHaveText (syn2)) {
        if (StringICmp (syn1, syn2) == 0) {
          redundantgenexref = TRUE;
          label = syn1;
        }
      }
    }
    if (redundantgenexref) {
      MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
      dsd.max = INT4_MAX;
      dsd.num_at_max = 0;
      dsd.equivalent_genes = FALSE;
      dsd.grp_at_max = NULL;
      count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                               LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);
      if (dsd.num_at_max > 1) {
        redundantgenexref = FALSE;
      }
    }
    if (redundantgenexref) {
      if (StringHasNoText (label)) {
        label = "?";
      }
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryGeneXref, "Unnecessary gene cross-reference %s", label);
    } else {
      if ((! dicistronic) && GPSorNTorNCorNGorNW (vsp->sep, sfp->location)) {
        /*
        SeqEntryToBioSource (vsp->sep, NULL, NULL, 0, &biop);
        */
        bsp = BioseqFindFromSeqLoc (sfp->location);
        BioseqToGeneticCode (bsp, NULL, NULL, NULL, NULL, 0, &biop);
        if (biop != NULL) {
          orp = biop->org;
          if (orp != NULL) {
            /* curated fly source still has duplicate features */
            if (StringICmp (orp->taxname, "Drosophila melanogaster") == 0) {
              if (StringHasNoText (label)) {
                label = "?";
              }
              if (sfpy != NULL && SeqLocAinB (sfp->location, sfpy->location) >= 0 &&
                  ValStrandsMatch (SeqLocStrand (sfp->location), SeqLocStrand (sfpy->location))) {
                /* cross-reference needed to disambiguate between multiple overlapping genes, ignore */
              } else {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_SuspiciousGeneXref, "Curated Drosophila record should not have gene cross-reference %s", label);
              }
            }
          }
        }
      }
    }
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryGeneXref, "Gene feature has gene cross-reference");
    }
    operon = SeqMgrGetOverlappingOperon (sfp->location, &fcontext);
    if (operon != NULL) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, 0, 0, 0, sfp, &fcontext) == sfp) {
        if (! StringHasNoText (fcontext.label)) {
          for (gbq = operon->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "operon") == 0) {
              if (StringICmp (gbq->val, fcontext.label) == 0) {
                ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidQualifierValue, "Operon is same as gene - %s", gbq->val);
              }
            }
          }
        }
      }
    }
  }
}

/*****************************************************************************
*
*   MrnaTransCheck (sfp, vsp)
*
*****************************************************************************/

static CharPtr bypass_mrna_trans_check [] = {
  "RNA editing",
  "reasons given in citation",
  "artificial frameshift",
  "transcribed product replaced",
  "unclassified transcription discrepancy",
  "mismatches in transcription",
  "adjusted for low-quality genome",
  NULL
};

NLM_EXTERN void MrnaTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  BioseqPtr       bsp;
  Char            ch;
  Int4            counta, countnona;
  CharPtr         farstr = "";
  ErrSev          fetchsev;
  GatherContextPtr  gcp;
  Boolean         has_errors = FALSE, unclassified_except = FALSE,
                  mismatch_except = FALSE, other_than_mismatch = FALSE;
  Int2            i;
  Boolean         is_refseq = FALSE;
  Int4            mismatch, total;
  CharPtr         mrseq, pdseq;
  Int4            mlen, plen;
  CharPtr         ptr1, ptr2;
  Boolean         report_errors = TRUE;
  ErrSev          sev;
  SeqFeat         sf;
  SeqIdPtr        sip, sip2, sip3;
  Boolean         unlockProd = FALSE;
  ValNode         vn;
  SeqDescrPtr     sdp;
  MolInfoPtr      mip;
  TextSeqIdPtr    tsip;

  if (sfp == NULL)
    return;
  if (sfp->pseudo)
    return;
  if (sfp->product == NULL)
    return;

  if (sfp->excpt && (! vsp->ignoreExceptions) && (! StringHasNoText (sfp->except_text))) {
    for (i = 0; bypass_mrna_trans_check [i] != NULL; i++) {
      if (StringISearch (sfp->except_text,  bypass_mrna_trans_check [i]) != NULL) {
        report_errors = FALSE;  /* biological exception */
      }
    }
    if (StringStr (sfp->except_text, "unclassified transcription discrepancy") != NULL) {
      unclassified_except = TRUE;
    }
    if (StringStr (sfp->except_text, "mismatches in transcription") != NULL) {
      mismatch_except = TRUE;
      report_errors = TRUE;
    }
  }

  sip = SeqLocId (sfp->product);
  if (sip == NULL)
    return;

  mrseq = GetSequenceByFeature (sfp);
  if (mrseq == NULL) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MrnaTransFail, "Unable to transcribe mRNA");
    return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL) {
    for (sip2 = bsp->id; sip2 != NULL; sip2 = sip2->next) {
      if (sip2->choice == SEQID_OTHER) {
        is_refseq = TRUE;
      }
    }
  }

  mismatch = 0;
  total = 0;

  sev = SEV_ERROR;
  gcp = vsp->gcp;
  if (gcp != NULL) {
    bsp = GetBioseqGivenSeqLoc (sfp->product, gcp->entityID);
    if (bsp == NULL) {
      /* if not local bioseq product, lower severity */
      sev = SEV_WARNING;
      if (is_refseq) {
        /* if refseq, restore higher severity */
        sev = SEV_ERROR;
      }
    }
    if (bsp == NULL && vsp->farFetchMRNAproducts) {
      if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
        bsp = BioseqLockById (sip);
      }
      if (bsp != NULL) {
        unlockProd = TRUE;
        farstr = "(far) ";
        if (sfp->partial) {
          sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
          if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              if (mip->completeness < 2 || mip->completeness > 5) {
                for (sip3 = bsp->id; sip3 != NULL; sip3 = sip3->next) {
                  if (sip3->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip3->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringNCmp (tsip->accession, "NM_", 3) == 0) {
                    /* if far NM_ record, return to lower severity */
                    sev = SEV_WARNING;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (bsp == NULL && (! vsp->farFetchMRNAproducts)) {
      goto erret;
    }
    if (bsp == NULL && sfp->product != NULL && vsp->farFetchMRNAproducts) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ProductFetchFailure, "Unable to fetch mRNA transcript");
      goto erret;
    }
  }
  if (is_refseq && unclassified_except) {
    /* if unclassified exception, drop back down to warning */
    sev = SEV_WARNING;
  }

  /* coerced feature on whole product for GetSequenceByFeature */

  MemSet ((Pointer) &sf, 0, sizeof (SeqFeat));
  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  sf.location = &vn;
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = sip;

  pdseq = GetSequenceByFeature (&sf);
  if (pdseq == NULL) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors || unclassified_except) {
      fetchsev = SEV_ERROR;
      if (sip->choice != SEQID_GI) {
        fetchsev = SEV_WARNING;
      }
      ValidErr (vsp, fetchsev, ERR_SEQ_FEAT_MrnaTransFail, "Unable to fetch mRNA transcript");
    }
  }
  if (pdseq != NULL) {
    mlen = StringLen (mrseq);
    plen = StringLen (pdseq);
    if (mlen != plen) {
      if (mlen < plen) {
        ptr1 = pdseq + mlen;
        counta = 0;
        countnona = 0;
        ch = *ptr1;
        while (ch != '\0') {
          if (ch == 'A' || ch == 'a') {
            counta++;
          } else {
            countnona++;
          }
          ptr1++;
          ch = *ptr1;
        }
        if (counta < 19 * countnona) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors) {
            ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptLen, "Transcript length [%ld] less than %sproduct length [%ld], and tail < 95%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* even if it fails polyA test, allow base-by-base comparison on common length */
        } else if (counta > 0 && countnona == 0) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PolyATail, "Transcript length [%ld] less than %sproduct length [%ld], but tail is 100%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* if it passes polyA test, allow base-by-base comparison on common length */
        } else {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors) {
            ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_PolyATail, "Transcript length [%ld] less than %sproduct length [%ld], but tail >= 95%s polyA", (long) mlen, farstr, (long) plen, "%");
          }
          plen = mlen; /* if it passes polyA test, allow base-by-base comparison on common length */
        }
      } else {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptLen, "Transcript length [%ld] greater than %sproduct length [%ld]", (long) mlen, farstr, (long) plen);
        }
      }
    }
    if (mlen == plen && mlen > 0 && StringICmp (mrseq, pdseq) != 0) {
      mismatch = 0;
      total = 0;
      ptr1 = mrseq;
      ptr2 = pdseq;
      while (total < mlen) {
        if (*ptr1 != *ptr2) {
          mismatch++;
        }
        ptr1++;
        ptr2++;
        total++;
      }
      if (mismatch > 0) {
        has_errors = TRUE;
        if (report_errors && (! mismatch_except)) {
          ValidErr (vsp, sev, ERR_SEQ_FEAT_TranscriptMismatches,
                    "There are %ld mismatches out of %ld bases between the transcript and %sproduct sequence", (long) mismatch, (long) total, farstr);
        }
      }
    }
    MemFree (pdseq);
  }

erret:

  MemFree (mrseq);

  if (unlockProd) {
    BioseqUnlock (bsp);
  }

  if (! report_errors) {
    if (! has_errors) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "mRNA has exception but passes transcription test");
    } else if (unclassified_except && (! other_than_mismatch)) {
      if (mismatch * 50 <= total) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ErroneousException,
                  "mRNA has unclassified exception but only difference is %ld mismatches out of %ld bases",
                  (long) mismatch, (long) total);
      }
    }
  }
}

/*****************************************************************************
*
*   CdTransCheck(sfp)
*       Treatment of terminal 'X'
*          If either the protein or the translation end in 'X' (usually
*          due to partial last codon) it is ignored to minimize conflicts
*          between approaches to add the X or not in this case.
*
*****************************************************************************/
static CharPtr MapToNTCoords (SeqFeatPtr sfp, SeqIdPtr protID, Int4 pos)
{
  SeqLocPtr       nslp;
  SeqLocPtr       pslp;
  CharPtr         rsult;
  SeqPntPtr       spntp;

  rsult = NULL;
  if (sfp != NULL && protID != NULL && pos >= 0) {
    spntp = SeqPntNew ();
    pslp = ValNodeNew (NULL);
    pslp->choice = SEQLOC_PNT;
    pslp->data.ptrvalue = (Pointer) spntp;
    spntp->point = pos;
    spntp->id = SeqIdDup (protID);
    nslp = aaLoc_to_dnaLoc (sfp, pslp);
    if (nslp != NULL) {
      rsult = SeqLocPrint (nslp);
    }
    SeqLocFree (pslp);
    SeqLocFree (nslp);
  }
  return rsult;
}

static Boolean Loc_is_RefSeq (SeqLocPtr location)
{
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;

  if (location == NULL)
    return FALSE;
  sip = SeqLocId (location);
  if (sip == NULL)
    return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL)
    return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean Loc_is_GEDL (SeqLocPtr location)
{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  if (location == NULL)
    return FALSE;
  sip = SeqLocId (location);
  if (sip == NULL)
    return FALSE;
  bsp = BioseqFind (sip);
  if (bsp == NULL)
    return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK) return TRUE;
    if (sip->choice == SEQID_EMBL) return TRUE;
    if (sip->choice == SEQID_DDBJ) return TRUE;
    if (sip->choice == SEQID_LOCAL) return TRUE;
  }
  return FALSE;
}

static void CdConflictCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  CharPtr       str1, str2;

  if (sfp == NULL || vsp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  str1 = GetSequenceByBsp (bsp);
  bs = TransTableTranslateCdRegion (NULL, sfp, FALSE, FALSE, TRUE);
  str2 = (CharPtr) BSMerge (bs, NULL);
  BSFree (bs);

  if (str1 != NULL && str2 != NULL && StringCmp (str1, str2) == 0) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_BadConflictFlag, "Coding region conflict flag should not be set");
  } else {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ConflictFlagSet, "Coding region conflict flag is set");
  }

  MemFree (str1);
  MemFree (str2);
}

static CharPtr bypass_cds_trans_check [] = {
  "RNA editing",
  "reasons given in citation",
  "artificial frameshift",
  "rearrangement required for product",
  "translated product replaced",
  "unclassified translation discrepancy",
  "mismatches in translation",
  "adjusted for low-quality genome",
  NULL
};

static void ValidateTranslExcept (
  ValidStructPtr vsp,
  SeqFeatPtr sfp,
  ValNodePtr codebreakhead,
  Boolean farFetchProd,
  Uint1 frame,
  ValNodePtr genetic_code
)

{
  Boolean       alt_start = FALSE;
  CdRegion      cr;
  ByteStorePtr  newprot = NULL;
  CharPtr       protseq = NULL;
  Int4          prot2len, i;
  SeqFeat       sf;
  ValNodePtr    vnp;

  MemSet ((Pointer) &sf, 0, sizeof (SeqFeat));
  MemSet ((Pointer) &cr, 0, sizeof (CdRegion));
  sf.data.choice = SEQFEAT_CDREGION;
  sf.data.value.ptrvalue = (Pointer) &cr;
  sf.location = sfp->location;
  cr.frame = frame;
  cr.genetic_code = genetic_code;

  newprot = ProteinFromCdRegionExEx (&sf, TRUE, FALSE, &alt_start, farFetchProd);
  if (newprot == NULL) return;
  protseq = BSMerge (newprot, NULL);
  BSFree (newprot);
  if (protseq == NULL) return;
  prot2len = StringLen (protseq);
  for (vnp = codebreakhead; vnp != NULL; vnp = vnp->next) {
    i = vnp->data.intvalue;
    if (i >= 0 && i < prot2len) {
      if (protseq [i] == (Char) vnp->choice) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryTranslExcept,
                  "Unnecessary transl_except %c at position %ld",
                  (char) vnp->choice, (long) (i + 1));
      }
    } else if (i == prot2len) {
      if ((Char) vnp->choice != '*') {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryTranslExcept,
                  "Unexpected transl_except %c at position %ld just past end of protein",
                  (char) vnp->choice, (long) (i + 1));
      }
    }
  }
  MemFree (protseq);
}

NLM_EXTERN void CdTransCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  ByteStorePtr    newprot = NULL;
  CharPtr         protseq = NULL;
  BioseqPtr       prot1seq = NULL, prot2seq = NULL;
  Int4            prot1len = 0, prot2len, i, len = 0, x_count = 0,
                  nonx_count = 0, xcount1 = 0, xcount2 = 0;
  CdRegionPtr     crp;
  SeqIdPtr        protid = NULL;
  Int2            residue1, residue2, stop_count = 0, mismatch = 0, ragged = 0;
  Boolean         got_stop = FALSE;
  /*
  SeqPortPtr      spp = NULL;
  */
  Uint2           part_loc = 0, part_prod = 0;
  Boolean         no_end = FALSE, no_beg = FALSE, show_stop = FALSE,
                  got_dash = FALSE, alt_start = FALSE, done;
  GBQualPtr       gb;
  ValNodePtr      vnp, code, codebreakhead = NULL;
  int             gccode = 0;
  Boolean         transl_except = FALSE, prot_ok = TRUE, is_nc = FALSE,
                  has_errors = FALSE, report_errors = TRUE,
                  unclassified_except = FALSE, mismatch_except = FALSE,
                  frameshift_except = FALSE, rearrange_except = FALSE,
                  other_than_mismatch = FALSE;
  CharPtr         nuclocstr, farstr = "";
  CodeBreakPtr    cbp;
  Int4            pos1, pos2, pos;
  SeqLocPtr       tmp;
  ErrSev          sev, trans_len_sev = SEV_ERROR;
  SeqEntryPtr     sep;
  Boolean         unlockProd = FALSE;
  StreamCache     sc;
  Boolean         isgap;
  Boolean         badseq;
  BioseqPtr       bsp;
  SeqIdPtr        sip, sip3;
  Boolean         is_ged = FALSE;
  Boolean         is_refseq = FALSE;
  Boolean         has_gi = FALSE;
  Boolean         farFetchProd;
  SeqDescrPtr     sdp;
  MolInfoPtr      mip;
  TextSeqIdPtr    tsip;

  if (sfp == NULL)
    return;

  for (gb = sfp->qual; gb != NULL; gb = gb->next) {     /* pseuogene */
    if (!StringICmp ("pseudo", gb->qual))
      return;
  }

  if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
    vsp->far_fetch_failure = TRUE;
    return;
  }

  if (sfp->excpt && (! vsp->ignoreExceptions) && (! StringHasNoText (sfp->except_text))) {
    for (i = 0; bypass_cds_trans_check [i] != NULL; i++) {
      if (StringISearch (sfp->except_text,  bypass_cds_trans_check [i]) != NULL) {
        report_errors = FALSE;  /* biological exception */
      }
    }
    if (StringStr (sfp->except_text, "unclassified translation discrepancy") != NULL) {
      unclassified_except = TRUE;
    }
    if (StringStr (sfp->except_text, "mismatches in translation") != NULL) {
      mismatch_except = TRUE;
      report_errors = TRUE;
    }
    if (StringStr (sfp->except_text, "artificial frameshift") != NULL) {
      frameshift_except = TRUE;
      report_errors = TRUE;
    }
    if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
      rearrange_except = TRUE;
    }
  }

  crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
  if (crp->code_break == NULL) {        /* check for unparsed transl_except */
    for (gb = sfp->qual; gb != NULL; gb = gb->next) {
      if (StringCmp (gb->qual, "transl_except") == 0) {
        transl_except = TRUE;
        break;
      }
    }
  } else {
    codebreakhead = MakeCodeBreakList (sfp->location, SeqLocLen (sfp->location), crp->code_break, crp->frame);
  }

  if (crp->genetic_code != NULL) {
    for (vnp = crp->genetic_code->data.ptrvalue; ((vnp != NULL) && (!gccode)); vnp = vnp->next) {
      switch (vnp->choice) {
      case 0:
        break;
      case 1:                  /* name */
        code = GeneticCodeFind (0, (CharPtr) (vnp->data.ptrvalue));
        if (code != NULL) {
          for (vnp = code->data.ptrvalue; ((vnp != NULL) && (!gccode)); vnp = vnp->next) {
            if (vnp->choice == 2)       /* id */
              gccode = (int) (vnp->data.intvalue);
          }
        }
        break;
      case 2:                  /* id */
        gccode = (int) (vnp->data.intvalue);
        break;
      default:
        gccode = 255;
        break;
      }
    }
  }

  farFetchProd = (Boolean) (vsp->farFetchCDSproducts || vsp->farFetchMRNAproducts);
  newprot = ProteinFromCdRegionExEx (sfp, TRUE, FALSE, &alt_start, farFetchProd);   /* include stop codons, do not remove trailing X/B/Z */
  if (newprot == NULL) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors || unclassified_except) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_CdTransFail, "Unable to translate");
    }
    prot_ok = FALSE;
    goto erret;
  }

  if (codebreakhead != NULL) {
    ValidateTranslExcept (vsp, sfp, codebreakhead, farFetchProd, crp->frame, crp->genetic_code);
  }

  if (alt_start && gccode == 1) {
    /* sev = SEV_WARNING; */
    sev = SEV_NONE; /* only enable for RefSeq, leave old code in for now */
    if (Loc_is_RefSeq (sfp->location)) {
      sev = SEV_ERROR;
    } else if (Loc_is_GEDL (sfp->location)) {
      sev = SEV_NONE;
    }
    if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
      if (StringStr (sfp->except_text, "alternative start codon") != NULL) {
        sev = SEV_NONE;
      }
    }
    if (sev > SEV_NONE) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_AltStartCodon, "Alternative start codon used");
      }
    }
  } else if (! alt_start) {
    if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
      if (StringStr (sfp->except_text, "alternative start codon") != NULL) {
        if (Loc_is_RefSeq (sfp->location)) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_AltStartCodon, "Unnecessary alternative start codon exception");
        }
      }
    }
  }

  part_loc = SeqLocPartialCheck (sfp->location);
  part_prod = SeqLocPartialCheckEx (sfp->product, farFetchProd);
  if ((part_loc & SLP_STOP) || (part_prod & SLP_STOP))
    no_end = TRUE;
  else {                        /* complete stop, so check for ragged end */

    len = SeqLocLen (sfp->location);
    if (crp->frame > 1)
      len -= (Int4) (crp->frame - 1);
    ragged = (Int2) (len % (Int4) (3));
    if (ragged) {
      len = SeqLocLen (sfp->location);
      cbp = crp->code_break;
      while (cbp != NULL) {
        pos1 = INT4_MAX;
        pos2 = -10;
        tmp = NULL;
        while ((tmp = SeqLocFindNext (cbp->loc, tmp)) != NULL) {
          pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_START);
          if (pos < pos1)
            pos1 = pos;
          pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_STOP);
          if (pos > pos2)
            pos2 = pos;
        }
        pos = pos2 - pos1;      /* codon length */
        if (pos >= 0 && pos <= 1 && pos2 == len - 1)
        {                       /*  a codon */
          /* allowing a partial codon at the end */
          ragged = 0;
        }

        cbp = cbp->next;
      }
    }
  }

  /* check for code break not on a codon */
  len = SeqLocLen (sfp->location);
  cbp = crp->code_break;
  while (cbp != NULL) {
    pos1 = INT4_MAX;
    pos2 = -10;
    tmp = NULL;
    while ((tmp = SeqLocFindNext (cbp->loc, tmp)) != NULL) {
      pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_START);
      if (pos < pos1)
        pos1 = pos;
      pos = GetOffsetInLoc (tmp, sfp->location, SEQLOC_STOP);
      if (pos > pos2)
        pos2 = pos;
    }
    pos = pos2 - pos1;          /* codon length */
    /* check for code break not on a codon */
    if (pos == 2 || (pos >= 0 && pos <= 1 && pos2 == len - 1)) {
      if (crp->frame == 2)
        pos = 1;
      else if (crp->frame == 3)
        pos = 2;
      else
        pos = 0;
      if ((pos1 % 3) != pos) {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExceptPhase, "transl_except qual out of frame.");
        }
      }
    }


    cbp = cbp->next;
  }

  if (crp->frame > 1) {
    if (!(part_loc & SLP_START)) {
      sev = SEV_WARNING;
      if (Loc_is_RefSeq (sfp->location)) {
        sev = SEV_ERROR;
      }
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_SuspiciousFrame, "Suspicious CDS location - frame > 1 but not 5' partial");
      }
    } else if ((part_loc & SLP_NOSTART) && (!PartialAtSpliceSiteOrGap (vsp, sfp->location, SLP_NOSTART, &isgap, &badseq))) {
      sev = SEV_INFO;
      if (Loc_is_RefSeq (sfp->location)) {
        sev = SEV_ERROR;
      }
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_SuspiciousFrame, "Suspicious CDS location - frame > 1 and not at consensus splice site");
      }
    }
  }

  if ((part_loc & SLP_START) || (part_prod & SLP_START))
    no_beg = TRUE;

  protseq = BSMerge (newprot, NULL);
  prot2len = StringLen (protseq);
  if (protseq != NULL) {
    len = prot2len;
    for (i = 0; i < len; i++) {
      residue1 = protseq [i];
      if ((i == 0) && (residue1 == '-'))
        got_dash = TRUE;
      if (residue1 == '*') {
        if (i == (len - 1))
          got_stop = TRUE;
        else
          stop_count++;
      }
      if (residue1 == 'X') {
        x_count++;
      } else {
        nonx_count++;
      }
    }
    if (x_count > nonx_count) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_CDShasTooManyXs, "CDS translation consists of more than 50%s X residues", "%");
    }
  }

  /*
  prot2len = BSLen (newprot);
  len = prot2len;
  BSSeek (newprot, 0, SEEK_SET);
  for (i = 0; i < len; i++) {
    residue1 = BSGetByte (newprot);
    if ((i == 0) && (residue1 == '-'))
      got_dash = TRUE;
    if (residue1 == '*') {
      if (i == (len - 1))
        got_stop = TRUE;
      else
        stop_count++;
    }
  }
  */

  if (stop_count > 0) {
    if (got_dash) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_StartCodon,
                  "Illegal start codon (and %ld internal stops). Probably wrong genetic code [%d]", (long) stop_count, gccode);
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InternalStop, "%ld internal stops (and illegal start codon). Genetic code [%d]", (long) stop_count, gccode);
      }
    } else {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          for (sip = bsp->id; sip != NULL; sip = sip->next) {
            switch (sip->choice) {
              case SEQID_GI :
                has_gi = TRUE;
                break;
              case SEQID_GENBANK :
              case SEQID_EMBL :
              case SEQID_DDBJ :
              case SEQID_TPG :
              case SEQID_TPE :
              case SEQID_TPD :
                is_ged = TRUE;
                break;
              case SEQID_OTHER :
                is_refseq = TRUE;
                break;
              default :
                break;
            }
          }
          if (has_gi && is_ged && (! is_refseq)) {
            sev = SEV_REJECT;
          }
        }
        ValidErr (vsp, sev, ERR_SEQ_FEAT_InternalStop, "%ld internal stops. Genetic code [%d]", (long) stop_count, gccode);
      }
    }
    prot_ok = FALSE;
    if (stop_count > 5)
      goto erret;
  } else if (got_dash) {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
    }
  }

  show_stop = TRUE;

  protid = SeqLocId (sfp->product);
  if (protid != NULL) {
    prot1seq = BioseqFind (protid);
    if (prot1seq == NULL && vsp->farFetchCDSproducts) {
      if (protid != NULL && (protid->choice != SEQID_GI || protid->data.intvalue > 0)) {
        prot1seq = BioseqLockById (protid);
      }
      if (prot1seq != NULL) {
        unlockProd = TRUE;
        farstr = "(far) ";
        if (sfp->partial) {
          sdp = GetNextDescriptorUnindexed (prot1seq, Seq_descr_molinfo, NULL);
          if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip != NULL) {
              if (mip->completeness < 2 || mip->completeness > 5) {
                for (sip3 = prot1seq->id; sip3 != NULL; sip3 = sip3->next) {
                  if (sip3->choice != SEQID_OTHER) continue;
                  tsip = (TextSeqIdPtr) sip3->data.ptrvalue;
                  if (tsip == NULL) continue;
                  if (StringNCmp (tsip->accession, "NP_", 3) == 0) {
                    /* if far NP_ record, return to lower severity */
                    trans_len_sev = SEV_WARNING;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (prot1seq == NULL && (! vsp->farFetchCDSproducts)) {
      goto erret;
    }
    if (prot1seq == NULL && sfp->product != NULL && vsp->farFetchCDSproducts) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_ProductFetchFailure, "Unable to fetch CDS product");
      goto erret;
    }
    if (prot1seq != NULL)
      prot1len = prot1seq->length;
  }

  if (prot1seq == NULL) {
    if (prot2len > 6) {
      if (! NGorNT (vsp->sep, sfp->location, &is_nc)) {
        sev = SEV_ERROR;
        if (DeltaOrFarSeg (vsp->sep, sfp->location)) {
          sev = SEV_WARNING;
        }
        if (is_nc) {
          sev = SEV_WARNING;
          sep = vsp->sep;
          if (sep != NULL && IS_Bioseq (sep)) {
            sev = SEV_NONE;
          }
        }
        if (sev != SEV_NONE) {
          has_errors = TRUE;
          other_than_mismatch = TRUE;
          if (report_errors) {
            ValidErr (vsp, sev, ERR_SEQ_FEAT_NoProtein, "No protein Bioseq given");
          }
        }
      }
    }
    goto erret;
  }

  len = prot2len;

  if ((got_stop) && (len == (prot1len + 1))) {  /* ok, got stop */
    len--;
  }

  if (! StreamCacheSetup (prot1seq, NULL, STREAM_EXPAND_GAPS, &sc)) {
    goto erret;
  }
  /*
  spp = SeqPortNew (prot1seq, 0, -1, 0, Seq_code_ncbieaa);
  if (spp == NULL)
    goto erret;
  */

  /* ignore terminal 'X' from partial last codon if present */

  done = FALSE;
  if ((!done) && (prot1len)) {
    /* prime the cache at a reasonable position near the end */
    if (prot1len > 4000) {
      StreamCacheSetPosition (&sc, prot1len - 2000);
    }
    residue1 = StreamCacheGetResidue (&sc);
  }
  while ((!done) && (prot1len)) {
    StreamCacheSetPosition (&sc, prot1len - 1);
    residue1 = StreamCacheGetResidue (&sc);
    /*
    SeqPortSeek (spp, (prot1len - 1), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    */
    if (residue1 == 'X') {        /* remove terminal X */
      prot1len--;
      xcount1++;
    }
    else
      done = TRUE;
  }
  done = FALSE;
  while ((!done) && (len)) {
    /*
    BSSeek (newprot, (len - 1), SEEK_SET);
    residue2 = BSGetByte (newprot);
    */
    residue2 = protseq [len - 1];
    if (residue2 == 'X') {
      len--;
      xcount2++;
    }
    else
      done = TRUE;
  }

  if (xcount1 != xcount2) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TerminalXDiscrepancy,
              "Terminal X count for CDS translation (%ld) and protein product sequence (%ld) are not equal",
              (long) xcount2, (long) xcount1);
  }

  if (len == prot1len) {        /* could be identical */
    StreamCacheSetPosition (&sc, 0);
    /*
    SeqPortSeek (spp, 0, SEEK_SET);
    BSSeek (newprot, 0, SEEK_SET);
    */
    for (i = 0; i < len; i++) {
      residue1 = protseq [i];
      residue2 = StreamCacheGetResidue (&sc);
      /*
      residue1 = BSGetByte (newprot);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 != residue2) {
        prot_ok = FALSE;
        if (residue2 == INVALID_RESIDUE)
          residue2 = '?';
        sev = SEV_ERROR;
        if (residue2 == 'X') {
          if (residue1 == 'B' || residue1 == 'Z' || residue1 == 'J') {
            sev = SEV_WARNING;
          }
        }
        if (mismatch == 10) {
          has_errors = TRUE;
          if (report_errors && (! mismatch_except)) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA, "More than 10 mismatches. Genetic code [%d]", gccode);
          }
        } else if (i == 0) {
          if ((sfp->partial) && (!no_beg) && (!no_end)) { /* ok, it's partial */
            has_errors = TRUE;
            other_than_mismatch = TRUE;
            if (report_errors) {
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Start of location should probably be partial");
            }
          } else if (residue1 == '-') {
            has_errors = TRUE;
            other_than_mismatch = TRUE;
            if (report_errors) {
              if (! got_dash) {
                ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
              }
            }
          } else {
            nuclocstr = MapToNTCoords (sfp, protid, i);
            if (nuclocstr != NULL) {
              has_errors = TRUE;
              if (report_errors && (! mismatch_except)) {
                ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
                          "%sResidue %ld in protein [%c] != translation [%c] at %s", farstr, (long) (i + 1), (char) residue2, (char) residue1, nuclocstr);
              }
            } else {
              has_errors = TRUE;
              if (report_errors && (! mismatch_except)) {
                ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
                          "%sResidue %ld in protein [%c] != translation [%c]", farstr, (long) (i + 1), (char) residue2, (char) residue1);
              }
            }
            MemFree (nuclocstr);
          }
        } else if (mismatch < 10) {
          nuclocstr = MapToNTCoords (sfp, protid, i);
          if (nuclocstr != NULL) {
            has_errors = TRUE;
            if (report_errors && (! mismatch_except)) {
              ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
                        "%sResidue %ld in protein [%c] != translation [%c] at %s", farstr, (long) (i + 1), (char) residue2, (char) residue1, nuclocstr);
            }
          } else {
            has_errors = TRUE;
            if (report_errors && (! mismatch_except)) {
              ValidErr (vsp, sev, ERR_SEQ_FEAT_MisMatchAA,
                        "%sResidue %ld in protein [%c] != translation [%c]", farstr, (long) (i + 1), (char) residue2, (char) residue1);
            }
          }
          MemFree (nuclocstr);
        }
        mismatch++;
      }
    }
    /*
    spp = SeqPortFree (spp);
    */
  } else {
    has_errors = TRUE;
    other_than_mismatch = TRUE;
    if (report_errors) {
      ValidErr (vsp, trans_len_sev, ERR_SEQ_FEAT_TransLen, "Given protein length [%ld] does not match %stranslation length [%ld]", prot1len, farstr, len);
    }
  }

  if ((sfp->partial) && (!mismatch)) {
    if ((!no_beg) && (!no_end)) {       /* just didn't label */
      if (!got_stop) {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "End of location should probably be partial");
        }
      } else {
        has_errors = TRUE;
        other_than_mismatch = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "This SeqFeat should not be partial");
        }
      }
      show_stop = FALSE;
    }
  }

  if (unlockProd) {
    BioseqUnlock (prot1seq);
  }


erret:
  if (show_stop) {
    if ((!got_stop) && (!no_end)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NoStop, "Missing stop codon");
      }
    } else if ((got_stop) && (no_end)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Got stop codon, but 3'end is labeled partial");
      }
    } else if ((got_stop) && (!no_end) && (ragged)) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      sev = SEV_ERROR;
      if (unclassified_except) {
        sev = SEV_WARNING;
      }
      if (report_errors || unclassified_except) {
        ValidErr (vsp, sev, ERR_SEQ_FEAT_TransLen, "Coding region extends %d base(s) past stop codon", (int) ragged);
      }
    }
  }

  if (!prot_ok) {
    if (transl_except) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExcept, "Unparsed transl_except qual. Skipped");
      }
    }
  } else {
    if (transl_except) {
      has_errors = TRUE;
      other_than_mismatch = TRUE;
      if (report_errors) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_TranslExcept, "Unparsed transl_except qual (but protein is okay). Skipped");
      }
    }
  }

  if (prot2seq != NULL)
    BioseqFree (prot2seq);
  else
    BSFree (newprot);
  /*
  SeqPortFree (spp);
  */
  MemFree (protseq);
  ValNodeFree (codebreakhead);

  if (! report_errors) {
    if (! has_errors) {
      if ((! frameshift_except) && (! rearrange_except)) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "CDS has exception but passes translation test");
      }
    } else if (unclassified_except && (! other_than_mismatch)) {
      if (mismatch * 50 <= len) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_ErroneousException,
                  "CDS has unclassified exception but only difference is %ld mismatches out of %ld residues",
                  (long) mismatch, (long) len);
      }
    }
  }
}

/*****************************************************************************
*
*   SpliceCheck(sfp)
*      checks for GT/AG rule at splice junctions
*
*****************************************************************************/
#define NOVALUE 0
#define HADGT 1
#define NOGT 2

static void SpliceCheckEx (ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll)
{
  SeqLocPtr       slp, nxt, head;
  Uint1           strand = Seq_strand_unknown;
  /*
  SeqPortPtr      spp = NULL;
  */
  SeqIdPtr        last_sip = NULL, sip, id;
  Int2            total, ctr;
  BioseqPtr       bsp = NULL;
  Int4            strt, stp, len = 0, donor, acceptor;
  Int2            residue1, residue2;
  Char            tbuf[40];
  Boolean         reportAsError, first, last, firstPartial, lastPartial, has_errors = FALSE,
                  report_errors = TRUE, checkExonDonor, checkExonAcceptor, pseudo;
  int             severity;
  Uint2           partialflag;
  Boolean         gpsOrRefSeq = FALSE;
  SeqEntryPtr     sep;
  BioseqSetPtr    bssp;
  TextSeqIdPtr    tsip;
  StreamCache     sc;
  SeqInt          sint;
  ValNode         vn;
  SeqMgrFeatContext  context;
  SeqFeatPtr      mrna, gene;
  GeneRefPtr      grp;

  if (sfp == NULL)
    return;

  if (GetAppProperty ("NcbiSubutilValidation") != NULL)
    return;                     /* suppress if NCBISubValidate */


  /* specific biological exceptions suppress check */

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "ribosomal slippage") != NULL||
        StringISearch (sfp->except_text, "artificial frameshift") != NULL ||
        StringISearch (sfp->except_text, "nonconsensus splice site") != NULL ||
        StringISearch (sfp->except_text, "adjusted for low-quality genome") != NULL) {
      report_errors = FALSE;
    }
  }

  MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
  MemSet ((Pointer) &vn, 0, sizeof (ValNode));

  head = sfp->location;
  if (head == NULL)
    return;

  if (LocationIsFar (sfp->location) && NoFetchFunctions ()) {
    vsp->far_fetch_failure = TRUE;
    return;
  }

  reportAsError = FALSE;
  if (GetAppProperty ("SpliceValidateAsError") != NULL) {
    reportAsError = TRUE;
  }

  slp = NULL;
  total = 0;
  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    total++;
    if (slp->choice == SEQLOC_EQUIV)
      return;                   /* bail on this one */
    if (total == 1)
      strand = SeqLocStrand (slp);
    else {
      if (strand != SeqLocStrand (slp)) /* bail on mixed strand */
        return;
    }
  }

  if ((!checkAll) && total < 2)
    return;
  if (total < 1)
    return;

  slp = NULL;
  ctr = 0;

  first = TRUE;
  last = FALSE;
  firstPartial = FALSE;
  lastPartial = FALSE;

  /* genomic product set or NT_ contig always relaxes to SEV_WARNING */

  sep = vsp->sep;
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      gpsOrRefSeq = TRUE;
    }
  }

  slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE);
  while (slp != NULL) {
    nxt = SeqLocFindPart (head, slp, EQUIV_IS_ONE);
    last = (Boolean) (nxt == NULL);
    partialflag = SeqLocPartialCheck (slp);
    firstPartial = (Boolean) (first && (partialflag & SLP_START));
    lastPartial = (Boolean) (last && (partialflag & SLP_STOP));
    ctr++;
    sip = SeqLocId (slp);
    if (sip == NULL)
      break;

    /* genomic product set or NT_ contig always relaxes to SEV_WARNING */
    bsp = BioseqFind (sip);
    if (bsp != NULL) {
      for (id = bsp->id; id != NULL; id = id->next) {
        if (id->choice == SEQID_OTHER) {
          tsip = (TextSeqIdPtr) id->data.ptrvalue;
          if (tsip != NULL && tsip->accession != NULL) {
            /*
            if (StringNICmp (tsip->accession, "NT_", 3) == 0) {
              gpsOrRefSeq = TRUE;
            } else if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
              gpsOrRefSeq = TRUE;
            } else if (StringNICmp (tsip->accession, "NG_", 3) == 0) {
              gpsOrRefSeq = TRUE;
            } else if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
              gpsOrRefSeq = TRUE;
            } else if (StringNICmp (tsip->accession, "NR_", 3) == 0) {
              gpsOrRefSeq = TRUE;
            }
            */
            gpsOrRefSeq = TRUE;
          }
        }
      }
    }

    if ((ctr == 1) || (!SeqIdMatch (sip, last_sip))) {
      /* spp = SeqPortFree (spp); */
      bsp = NULL;
      if (sip != NULL && (sip->choice != SEQID_GI || sip->data.intvalue > 0)) {
        bsp = BioseqLockById (sip);
      }
      if (bsp == NULL)
        break;
      len = bsp->length;
      if (strand != Seq_strand_minus) {
        if (! StreamCacheSetup (bsp, NULL, EXPAND_GAPS_TO_DASHES, &sc)) {
          BioseqUnlock (bsp);
          break;
        }
      } else {
        sint.from = 0;
        sint.to = len - 1;
        sint.strand = strand;
        sint.id = sip;
        vn.choice = SEQLOC_INT;
        vn.data.ptrvalue = (Pointer) &sint;
        vn.next = NULL;
        if (! StreamCacheSetup (NULL, &vn, EXPAND_GAPS_TO_DASHES, &sc)) {
          BioseqUnlock (bsp);
          break;
        }
      }
      /* spp = SeqPortNew (bsp, 0, -1, strand, Seq_code_ncbi4na); */
      BioseqUnlock (bsp);
      /*
      if (spp == NULL)
        break;
      */
      last_sip = sip;
    }

    acceptor = SeqLocStart (slp);
    donor = SeqLocStop (slp);

    if (acceptor < 0 || acceptor >= len || donor < 0 || donor >= len) {
      /*
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_Range,
                "Unable to check splice consensus because feature outside range of sequence");
      */
      return;
    }

    if (strand != Seq_strand_minus) {
      strt = acceptor;
      stp = donor;
    } else {
      strt = donor;
      donor = acceptor;
      acceptor = strt;
      stp = len - donor - 1;    /* orient to reverse complement seqport */
      strt = len - acceptor - 1;
    }

    checkExonDonor = FALSE;
    checkExonAcceptor = FALSE;
    if (checkAll) {
      pseudo = FALSE;
      grp = SeqMgrGetGeneXref (sfp);
      if (grp == NULL) {
        gene = SeqMgrGetOverlappingGene (sfp->location, &context);
        if (gene != NULL) {
          pseudo = gene->pseudo;
        }
      }
      if (! pseudo) {
        checkExonDonor = TRUE;
        checkExonAcceptor = TRUE;
        mrna = SeqMgrGetOverlappingmRNA (sfp->location, &context);
        if (mrna != NULL /* && (! mrna->partial) */ ) {
          if (strand != Seq_strand_minus) {
            if (donor == SeqLocStop (mrna->location) && (! context.partialR)) {
              checkExonDonor = FALSE;
            }
            if (acceptor == SeqLocStart (mrna->location) && (! context.partialL)) {
              checkExonAcceptor = FALSE;
            }
          } else {
            if (donor == SeqLocStart (mrna->location) && (! context.partialR)) {
              checkExonDonor = FALSE;
            }
            if (acceptor == SeqLocStop (mrna->location) && (! context.partialL)) {
              checkExonAcceptor = FALSE;
            }
          }
        }
      }
    }

    if (((checkExonDonor && (!lastPartial)) || ctr < total) && (stp < (len - 2))) {   /* check donor on all but last exon and on sequence */
      tbuf[0] = '\0';
      StreamCacheSetPosition (&sc, stp + 1);
      residue1 = StreamCacheGetResidue (&sc);
      residue2 = StreamCacheGetResidue (&sc);
      /*
      SeqPortSeek (spp, (stp + 1), SEEK_SET);
      residue1 = SeqPortGetResidue (spp);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 == '-' && residue2 == '-') {
        /* ignore gap, and suppress UnnecessaryException message */
        has_errors = TRUE;
      } else if (IS_residue (residue1) && IS_residue (residue2)) {
        if (residue1 != 'G' || residue2 != 'T') {        /* not T */
          if (residue1 == 'G' && residue2 == 'C') {       /* GC minor splice site */
            tbuf[0] = '\0';
            if (bsp == NULL) {
              StringCpy (tbuf, "?");
            } else if (vsp->suppressContext || vsp->convertGiToAccn) {
              WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            } else {
              BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            }
            tbuf[39] = '\0';
            if (RareConsensusNotExpected (sfp)) {
              has_errors = TRUE;
              if (report_errors) {
                ValidErr (vsp, SEV_INFO, ERR_SEQ_FEAT_RareSpliceConsensusDonor,
                          "Rare splice donor consensus (GC) found instead of (GT) after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
              }
            }
          } else {
            if (gpsOrRefSeq) {
              severity = SEV_WARNING;
            } else if (checkExonDonor) {
              severity = SEV_WARNING;
            } else if (reportAsError) {
              severity = SEV_ERROR;
            } else {
              severity = SEV_WARNING;
            }
            tbuf[0] = '\0';
            if (bsp == NULL) {
              StringCpy (tbuf, "?");
            } else if (vsp->suppressContext || vsp->convertGiToAccn) {
              WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            } else {
              BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
            }
            tbuf[39] = '\0';
            has_errors = TRUE;
            if (report_errors) {
              ValidErr (vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                        "Splice donor consensus (GT) not found after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
            }
          }
        }
      } else {
        has_errors = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusDonor,
                    "Bad sequence at splice donor after exon ending at position %ld of %s", (long) (donor + 1), tbuf);
        }
      }
    }

    if (((checkExonAcceptor && (!firstPartial)) || ctr != 1) && (strt > 1)) {
      StreamCacheSetPosition (&sc, strt - 2);
      residue1 = StreamCacheGetResidue (&sc);
      residue2 = StreamCacheGetResidue (&sc);
      /*
      SeqPortSeek (spp, (strt - 2), SEEK_SET);
      residue1 = SeqPortGetResidue (spp);
      residue2 = SeqPortGetResidue (spp);
      */
      if (residue1 == '-' && residue2 == '-') {
        /* ignore gap, and suppress UnnecessaryException message */
        has_errors = TRUE;
      } else if (IS_residue (residue1) && IS_residue (residue2)) {
        if (residue1 != 'A' || residue2 != 'G') {
          if (gpsOrRefSeq) {
            severity = SEV_WARNING;
          } else if (checkExonAcceptor) {
            severity = SEV_WARNING;
          } else if (reportAsError) {
            severity = SEV_ERROR;
          } else {
            severity = SEV_WARNING;
          }
          tbuf[0] = '\0';
          if (bsp == NULL) {
            StringCpy (tbuf, "?");
            SeqIdWrite (sip, tbuf, PRINTID_FASTA_SHORT, 39);
          } else if (vsp->suppressContext || vsp->convertGiToAccn) {
            WorstBioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
          } else {
            BioseqLabel (bsp, tbuf, 39, OM_LABEL_CONTENT);
          }
          tbuf[39] = '\0';
          has_errors = TRUE;
          if (report_errors) {
            ValidErr (vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                      "Splice acceptor consensus (AG) not found before exon starting at position %ld of %s", (long) (acceptor + 1), tbuf);
          }
        }
      } else {
        has_errors = TRUE;
        if (report_errors) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensusAcceptor,
                    "Bad sequence at splice acceptor before exon starting at position %ld of %s", (long) (acceptor + 1), tbuf);
        }
      }
    }

    first = FALSE;
    slp = nxt;
  }

  /* SeqPortFree (spp); */

  if (! report_errors) {
    if (! has_errors) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryException, "feature has exception but passes splice site test");
    }
  }
}

NLM_EXTERN void SpliceCheck (ValidStructPtr vsp, SeqFeatPtr sfp)
{
  SpliceCheckEx (vsp, sfp, FALSE);
}

/*****************************************************************************
*
*   CdsProductIdCheck (vsp, sfp)
*      code taken from asn2gnbk.c - release mode expects CDS product Bioseqs
*
*****************************************************************************/
static void CdsProductIdCheck (ValidStructPtr vsp, SeqFeatPtr sfp)

{
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  Boolean            juststop = FALSE;
  Boolean            okay = FALSE;
  SeqEntryPtr        oldscope;
  Boolean            partial5;
  Boolean            partial3;
  Boolean            pseudo = FALSE;
  SeqEntryPtr        sep;

   /* non-pseudo CDS must have /product */
   if (sfp->pseudo) {
     pseudo = TRUE;
   }
   grp = SeqMgrGetGeneXref (sfp);
   if (grp == NULL) {
     sep = GetTopSeqEntryForEntityID (sfp->idx.entityID);
     oldscope = SeqEntrySetScope (sep);
     gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
     SeqEntrySetScope (oldscope);
     if (gene != NULL) {
       grp = (GeneRefPtr) gene->data.value.ptrvalue;
       if (gene->pseudo) {
         pseudo = TRUE;
       }
     }
   }
   if (grp != NULL && grp->pseudo) {
     pseudo = TRUE;
   }
   if (sfp->location != NULL) {
     if (CheckSeqLocForPartial (sfp->location, &partial5, &partial3)) {
       if (partial5 && (! partial3)) {
         if (SeqLocLen (sfp->location) <= 5) {
           juststop = TRUE;
         }
       }
     }
   }
   if (pseudo || juststop) {
     okay = TRUE;
   } else if (sfp->product != NULL) {
     okay = TRUE;
   } else {
     if (sfp->excpt && (! StringHasNoText (sfp->except_text))) {
       if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
         okay = TRUE;
       }
     }
   }
   if (! okay) {
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MissingCDSproduct, "Expected CDS product absent");
   }
}

/*****************************************************************************
*
*   ValidateSeqLoc(vsp, slp, prefix)
*
*****************************************************************************/

static Int2 SeqLocMixCount (SeqLocPtr slp)

{
  Int2       count = 0;
  SeqLocPtr  loc;

  if (slp == NULL) return 0;

  while (slp != NULL) {
    if (slp->choice == SEQLOC_MIX) {
      count++;
      loc = (SeqLocPtr) slp->data.ptrvalue;
      count += SeqLocMixCount (loc);
    }
    slp = slp->next;
  }

  return count;
}

NLM_EXTERN void ValidateSeqLoc (ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix)
{
  SeqLocPtr       tmp, prev;
  Boolean         retval = TRUE, tmpval, mixed_strand = FALSE, unmarked_strand = FALSE,
                  ordered = TRUE, adjacent = FALSE, circular = FALSE, exception = FALSE;
  CharPtr         ctmp;
  Uint1           strand2 = 0, strand1;
  ErrSev          sev;
  SeqIntPtr       sip1, sip2, prevsip;
  SeqPntPtr       spp;
  PackSeqPntPtr   pspp;
  SeqIdPtr        id1 = NULL, id2 = NULL;
  BioseqPtr       bsp;
  SeqFeatPtr      sfp = NULL;
  Int2            zeroGi = 0;
  Char            buf [32];
  SeqIdPtr        sip;

  if (slp == NULL)
    return;

  sfp = vsp->sfp;

  tmp = NULL;
  while ((tmp = SeqLocFindNext (slp, tmp)) != NULL) {
    sip = SeqLocId (tmp);
    if (sip != NULL && sip->choice == SEQID_GI && sip->data.intvalue <= 0) {
      zeroGi++;
    }
  }
  if (zeroGi > 0) {
    StringCpy (buf, "?");
    bsp = vsp->bsp;
    if (bsp != NULL) {
      SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    }
    if (zeroGi > 1) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_FeatureLocationIsGi0, "Feature has %d gi|0 locations on Bioseq %s",
                (int) zeroGi, buf);
    } else if (zeroGi > 0) {
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_FeatureLocationIsGi0, "Feature has %d gi|0 location on Bioseq %s",
                (int) zeroGi, buf);
    }
  }

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL && bsp->topology == 2) {
    circular = TRUE;
  }

  if (SeqLocMixCount (slp) > 1) {
      retval = FALSE;
      ctmp = SeqLocPrint (slp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_NestedSeqLocMix, "%s: SeqLoc [%s] has nested SEQLOC_MIX elements", prefix, ctmp);
      MemFree (ctmp);
  }

  tmp = NULL;
  prev = NULL;
  sip1 = NULL;
  prevsip = NULL;
  strand1 = Seq_strand_other;
  while ((tmp = SeqLocFindNext (slp, tmp)) != NULL) {
    tmpval = TRUE;
    switch (tmp->choice) {
    case SEQLOC_INT:
      sip1 = prevsip;
      sip2 = (SeqIntPtr) (tmp->data.ptrvalue);
      strand2 = sip2->strand;
      id2 = sip2->id;
      tmpval = SeqIntCheck (sip2);
      if ((tmpval) && (sip1 != NULL) && (ordered) && (! circular)) {
        if (SeqIdForSameBioseq (sip1->id, sip2->id)) {
          if (strand2 == Seq_strand_minus) {
            if (sip1->to < sip2->to)
              ordered = FALSE;
            if (sip2->to + 1 == sip1->from)
              adjacent = TRUE;
          } else {
            if (sip1->to > sip2->to)
              ordered = FALSE;
            if (sip1->to + 1 == sip2->from)
              adjacent = TRUE;
          }
        }
      }
      if (prevsip != NULL) {
        if (SeqIdForSameBioseq (prevsip->id, sip2->id)) {
          if (prevsip->strand == sip2->strand && prevsip->from == sip2->from && prevsip->to == sip2->to) {
            ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_DuplicateInterval, "Duplicate exons in location");
          }
        }
      }
      prevsip = sip2;
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) (tmp->data.ptrvalue);
      strand2 = spp->strand;
      id2 = spp->id;
      tmpval = SeqPntCheck (spp);
      prevsip = NULL;
      break;
    case SEQLOC_PACKED_PNT:
      pspp = (PackSeqPntPtr) (tmp->data.ptrvalue);
      strand2 = pspp->strand;
      id2 = pspp->id;
      tmpval = PackSeqPntCheck (pspp);
      prevsip = NULL;
      break;
    case SEQLOC_NULL:
      break;
    default:
      strand2 = Seq_strand_other;
      id2 = NULL;
      prevsip = NULL;
      break;
    }
    if (!tmpval) {
      retval = FALSE;
      ctmp = SeqLocPrint (tmp);
      if (ctmp != NULL && StringLen (ctmp) > 800) {
        StringCpy (ctmp + 797, "...");
      }
      ValidErr (vsp, SEV_REJECT, ERR_SEQ_FEAT_Range, "%s: SeqLoc [%s] out of range", prefix, ctmp);
      MemFree (ctmp);

    }

    if (tmp->choice != SEQLOC_NULL) {
      if ((strand1 != Seq_strand_other) && (strand2 != Seq_strand_other)) {
        if (SeqIdForSameBioseq (id1, id2)) {
          if (strand1 != strand2) {
            if (strand1 == Seq_strand_plus && strand2 == Seq_strand_unknown) {
              unmarked_strand = TRUE;
            } else if (strand1 == Seq_strand_unknown && strand2 == Seq_strand_plus) {
              unmarked_strand = TRUE;
            } else {
              mixed_strand = TRUE;
            }
          }
        }
      }

      strand1 = strand2;
      id1 = id2;
    }
  }

  if (sfp != NULL) {

    /* Publication intervals ordering does not matter */

    if (sfp->idx.subtype == FEATDEF_PUB) {
      ordered = TRUE;
      adjacent = FALSE;
    }

    /* ignore ordering of heterogen bonds */

    if (sfp->data.choice == SEQFEAT_HET) {
      ordered = TRUE;
      adjacent = FALSE;
    }

    /* misc_recomb intervals SHOULD be in reverse order */

    if (sfp->idx.subtype == FEATDEF_misc_recomb) {
      ordered = TRUE;
    }

    /* primer_bind intervals MAY be in on opposite strands */

    if (sfp->idx.subtype == FEATDEF_primer_bind) {
      mixed_strand = FALSE;
      unmarked_strand = FALSE;
      ordered = TRUE;
    }

    if (sfp->excpt) {
      exception = TRUE;
    }
  }

  if (adjacent) {
    ctmp = SeqLocPrint (slp);
    if (exception) {
      sev = SEV_WARNING;
    } else {
      sev = SEV_ERROR;
    }
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, sev, ERR_SEQ_FEAT_AbuttingIntervals, "%s: Adjacent intervals in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }

  if (exception) {
    /* trans splicing exception turns off both mixed_strand and out_of_order messages */
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      return;
    }
  }

  if (mixed_strand || unmarked_strand || (!ordered)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    if (mixed_strand) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s]", prefix, ctmp);
    } else if (unmarked_strand) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed plus and unknown strands in SeqLoc [%s]", prefix, ctmp);
    }
    if (!ordered)
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_SeqLocOrder, "%s: Intervals out of order in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
    return;
  }

  if (sfp != NULL) {

    /* ignore special case features here as well */

    if (sfp->idx.subtype == FEATDEF_PUB ||
        sfp->data.choice == SEQFEAT_HET ||
        sfp->idx.subtype == FEATDEF_misc_recomb ||
        sfp->idx.subtype == FEATDEF_primer_bind)
      return;
  }

  /* newer check for intervals out of order on segmented bioseq */

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;

  if (SeqLocBadSortOrder (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_SeqLocOrder, "%s: Intervals out of order in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }

  /* newer check for mixed strand on segmented bioseq */

  if (SeqLocMixedStrands (bsp, slp)) {
    ctmp = SeqLocPrint (slp);
    if (ctmp != NULL && StringLen (ctmp) > 800) {
      StringCpy (ctmp + 797, "...");
    }
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s]", prefix, ctmp);
    MemFree (ctmp);
  }
}

/*****************************************************************************
*
*   SeqGraph validation section
*
*****************************************************************************/

typedef struct gphgetdata
{
  ValNodePtr      vnp;
  BioseqPtr       bsp;
}
GphGetData     , PNTR GphGetPtr;

typedef struct grphitem
{
  SeqGraphPtr     sgp;
  Int4            left;
  Int4            right;
  Int2            index;
}
GrphItem       , PNTR GrphItemPtr;

static void GetGraphsProc (SeqGraphPtr sgp, Pointer userdata)
{
  GphGetPtr       ggp;
  GrphItemPtr     gip;

  ggp = (GphGetPtr) userdata;
  if (ggp == NULL || sgp == NULL) return;
  /* only phrap or gap4 currently allowed */
  if (StringICmp (sgp->title, "Phrap Quality") == 0 || StringICmp (sgp->title, "Phred Quality") == 0 || StringICmp (sgp->title, "Gap4") == 0) {
    /* data type must be bytes */
    if (sgp->flags[2] == 3) {
      if (SeqIdIn (SeqLocId (sgp->loc), ggp->bsp->id)) {
        gip = (GrphItemPtr) MemNew (sizeof (GrphItem));
        if (gip == NULL) return;
        gip->sgp = sgp;
        gip->left = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_LEFT_END);
        gip->right = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_RIGHT_END);
        ValNodeAddPointer (&(ggp->vnp), 0, (Pointer) gip);
      }
    }
  }
  return;
}

static int LIBCALLBACK SortSeqGraphProc (VoidPtr ptr1, VoidPtr ptr2)
{
  GrphItemPtr     gip1, gip2;
  ValNodePtr      vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL)
    return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL)
    return 0;
  gip1 = (GrphItemPtr) vnp1->data.ptrvalue;
  gip2 = (GrphItemPtr) vnp2->data.ptrvalue;
  if (gip1 == NULL || gip2 == NULL)
    return 0;
  if (gip1->left > gip2->left) {
    return 1;
  } else if (gip1->left < gip2->left) {
    return -1;
  } else if (gip1->right > gip2->right) {
    return -1;
  } else if (gip2->right < gip2->right) {
    return 1;
  }
  return 0;
}

/* gets valnode list of sorted graphs in GrphItem structures */

static ValNodePtr GetSeqGraphsOnBioseq (Uint2 entityID, BioseqPtr bsp)
{
  GphGetData   ggd;
  GrphItemPtr  gip;
  Int2         index;
  ValNodePtr   vnp;

  ggd.vnp = NULL;
  ggd.bsp = bsp;
  VisitGraphsOnBsp (bsp, (Pointer) &ggd, GetGraphsProc);
  for (vnp = ggd.vnp, index = 1; vnp != NULL; vnp = vnp->next, index++) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip != NULL) {
      gip->index = index;
    }
  }
  ggd.vnp = ValNodeSort (ggd.vnp, SortSeqGraphProc);
  return ggd.vnp;
}

static Boolean NextLitLength (DeltaSeqPtr next, Int4Ptr lenp)

{
  SeqLitPtr  slp;

  if (lenp == NULL) return FALSE;
  *lenp = 0;
  if (next == NULL || next->choice != 2) return FALSE;
  slp = (SeqLitPtr) next->data.ptrvalue;
  if (slp == NULL || slp->seq_data == NULL) return FALSE;
  *lenp = slp->length;
  return TRUE;
}

static void ValidateGraphsOnBioseq (GatherContextPtr gcp)
{
  Byte            scores [400];
  ByteStorePtr    bs;
  BioseqPtr       bsp;
  Int2            k, val, index, scount;
  Int4            curroffset = 0, gphlen = 0, seqlen = 0, slplen,
                  bslen, min = INT4_MAX, max = INT4_MIN, j, lastloc = -1,
                  numBases, NsWithScore, GapsWithScore, ACGTsWithoutScore,
                  ambigWithoutScore, valsBelowMin, valsAboveMax,
                  firstN, firstACGT, firstAmbig, pos, litlen, nxtlen;
  FloatHi         pct;
  DeltaSeqPtr     dsp, next;
  Uint2           entityID, olditemtype = 0, numdsp = 0, numsgp = 0;
  Uint4           firstsgitemid = 0;
  Uint4           olditemid = 0;
  GrphItemPtr     gip;
  ValNodePtr      head, vnp;
  Boolean         outOfOrder = FALSE, fa2htgsBug = FALSE, overlaps = FALSE;
  Uint1           residue;
  SeqGraphPtr     sgp;
  SeqIntPtr       sintp;
  SeqLocPtr       slocp;
  SeqLitPtr       slp;
  StreamCache     sc;
  ValidStructPtr  vsp;
  Boolean         single_report_mode = TRUE;

  vsp = (ValidStructPtr) gcp->userdata;
  bsp = (BioseqPtr) gcp->thisitem;
  if (vsp == NULL || bsp == NULL)
    return;
  if (!ISA_na (bsp->mol))
    return;

  vsp->bsp = bsp;
  vsp->descr = NULL;
  vsp->sfp = NULL;
  vsp->bssp = (BioseqSetPtr) gcp->parentitem;

  if (SeqMgrGetParentOfPart (bsp, NULL) != NULL)
    return;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  head = GetSeqGraphsOnBioseq (entityID, bsp);
  if (head == NULL)
    return;

  olditemid = gcp->itemID;
  olditemtype = gcp->thistype;
  gcp->thistype = OBJ_SEQGRAPH;

  for (vnp = head, index = 1; vnp != NULL; vnp = vnp->next, index++) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL)
      continue;

    sgp = gip->sgp;
    if (sgp == NULL)
      continue;
    gcp->itemID = sgp->idx.itemID;
    if (firstsgitemid == 0) {
      firstsgitemid = sgp->idx.itemID;
    }

    if (gip->index != index) {
      outOfOrder = TRUE;
      if (gip->index == 129 && index == 2) {
        fa2htgsBug = TRUE;
      }
    }
    if (gip->left <= lastloc) {
      overlaps = TRUE;
    }
    lastloc = gip->right;
    min = MIN ((Int4) min, (Int4) sgp->min.intvalue);
    max = MAX ((Int4) max, (Int4) sgp->max.intvalue);

    if (sgp->min.intvalue < 0 || sgp->min.intvalue > 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphMin, "Graph min (%ld) out of range", (long) sgp->min.intvalue);
    }

    if (sgp->max.intvalue > 100) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphMax, "Graph max (%ld) out of range", (long) sgp->max.intvalue);
    }
    if (sgp->max.intvalue <= 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphMax, "Graph max (%ld) out of range", (long) sgp->max.intvalue);
    }

    gphlen += sgp->numval;
    bs = (ByteStorePtr) sgp->values;
    if (bs != NULL) {
      bslen = BSLen (bs);
      if (sgp->numval != bslen) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphByteLen, "SeqGraph (%ld) and ByteStore (%ld) length mismatch", (long) sgp->numval, (long) bslen);
      }
    }
  }
  if (outOfOrder) {
    gcp->itemID = firstsgitemid;
    if (fa2htgsBug) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphOutOfOrder, "Graph components are out of order - probably caused by old fa2htgs bug");
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphOutOfOrder, "Graph components are out of order - may be a software bug");
    }
  }
  if (overlaps) {
    gcp->itemID = firstsgitemid;
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphOverlap, "Graph components overlap, with multiple scores for a single base");
  }

  if (bsp->repr == Seq_repr_raw) {
    seqlen = bsp->length;
  } else if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
      switch (dsp->choice) {
      case 1:
        slocp = (SeqLocPtr) dsp->data.ptrvalue;
        if (slocp == NULL)
          break;
        if (slocp->choice != SEQLOC_NULL) {
          seqlen += SeqLocLen (slocp);
        }
        break;
      case 2:
        slp = (SeqLitPtr) dsp->data.ptrvalue;
        if (slp == NULL || slp->seq_data == NULL)
          break;
        seqlen += slp->length;
        break;
      default:
        break;
      }
    }
  }

  if (seqlen != gphlen && bsp->length != gphlen) {
    gcp->itemID = firstsgitemid;
    ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphBioseqLen, "SeqGraph (%ld) and Bioseq (%ld) length mismatch", (long) gphlen, (long) seqlen);
  }

  if (bsp->repr == Seq_repr_delta) {
    if (head != NULL && head->next != NULL) {
      for (dsp = (DeltaSeqPtr) (bsp->seq_ext), vnp = head; dsp != NULL && vnp != NULL; dsp = next) {
        next = dsp->next;
        gip = (GrphItemPtr) vnp->data.ptrvalue;
        if (gip == NULL)
          continue;
        sgp = gip->sgp;
        if (sgp == NULL)
          continue;
        switch (dsp->choice) {
        case 1:
          slocp = (SeqLocPtr) dsp->data.ptrvalue;
          if (slocp != NULL && slocp->choice != SEQLOC_NULL) {
            slplen = SeqLocLen (slocp);
            curroffset += slplen;
            if (sgp->numval != slplen) {
              gcp->itemID = sgp->idx.itemID;
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphSeqLocLen, "SeqGraph (%ld) and SeqLoc (%ld) length mismatch", (long) sgp->numval, (long) slplen);
            }
            numdsp++;
            if (vnp != NULL) {
              vnp = vnp->next;
              numsgp++;
            }
          }
          break;
        case 2:
          slp = (SeqLitPtr) dsp->data.ptrvalue;
          litlen = 0;
          if (slp != NULL) {
            litlen = slp->length;
          }
          if (slp != NULL && slp->seq_data != NULL) {
            while (NextLitLength (next, &nxtlen)) {
              litlen += nxtlen;
              next = next->next;
            }
            if (sgp->numval != litlen) {
              gcp->itemID = sgp->idx.itemID;
              ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphSeqLitLen, "SeqGraph (%ld) and SeqLit (%ld) length mismatch",
                        (long) sgp->numval, (long) litlen);
            }
            slocp = sgp->loc;
            if (slocp != NULL && slocp->choice == SEQLOC_INT) {
              sintp = (SeqIntPtr) slocp->data.ptrvalue;
              if (sintp != NULL) {
                if (sintp->from != curroffset) {
                  gcp->itemID = sgp->idx.itemID;
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphStartPhase, "SeqGraph (%ld) and SeqLit (%ld) start do not coincide",
                            (long) sintp->from, (long) curroffset);
                }
                if (sintp->to != litlen + curroffset - 1) {
                  gcp->itemID = sgp->idx.itemID;
                  ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphStopPhase, "SeqGraph (%ld) and SeqLit (%ld) stop do not coincide",
                            (long) sintp->to, (long) (litlen + curroffset - 1));
                }
              }
            }
            numdsp++;
            if (vnp != NULL) {
              vnp = vnp->next;
              numsgp++;
            }
          }
          if (slp != NULL) {
            curroffset += litlen;
          }
          break;
        default:
          break;
        }
      }
      for (dsp = (DeltaSeqPtr) (bsp->seq_ext), numdsp = 0; dsp != NULL; dsp = next) {
        next = dsp->next;
        switch (dsp->choice) {
        case 1:
          slocp = (SeqLocPtr) dsp->data.ptrvalue;
          if (slocp != NULL && slocp->choice != SEQLOC_NULL) {
            numdsp++;
          }
          break;
        case 2:
          slp = (SeqLitPtr) dsp->data.ptrvalue;
          if (slp != NULL && slp->seq_data != NULL) {
            while (NextLitLength (next, &nxtlen)) {
              next = next->next;
            }
            numdsp++;
          }
          break;
        default:
          break;
        }
      }
      for (vnp = head, numsgp = 0; vnp != NULL; vnp = vnp->next, numsgp++)
        continue;
      if (numdsp != numsgp) {
        gcp->itemID = firstsgitemid;
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphDiffNumber, "Different number of SeqGraph (%d) and SeqLit (%d) components", (int) numsgp, (int) numdsp);
      }
    }
  }

  numBases = 0;
  NsWithScore = 0;
  GapsWithScore = 0;
  ACGTsWithoutScore = 0;
  ambigWithoutScore = 0;
  valsBelowMin = 0;
  valsAboveMax = 0;
  firstN = -1;
  firstACGT = -1;
  firstAmbig = -1;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gip = (GrphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL)
      continue;
    sgp = gip->sgp;
    if (sgp == NULL)
      continue;
    if (! StreamCacheSetup (NULL, sgp->loc, EXPAND_GAPS_TO_DASHES, &sc)) continue;
    slplen = SeqLocLen (sgp->loc);

    bs = (ByteStorePtr) sgp->values;
    BSSeek (bs, 0, SEEK_SET);
    j = 0;
    val = 0;

    scount = (Int2) BSRead (bs, scores, sizeof (scores));
    k = 0;

    if (! single_report_mode) {
      numBases = 0;
      NsWithScore = 0;
      GapsWithScore = 0;
      ACGTsWithoutScore = 0;
      ambigWithoutScore = 0;
      valsBelowMin = 0;
      valsAboveMax = 0;
      firstN = -1;
      firstACGT = -1;
      firstAmbig = -1;
    }

    pos = gip->left;

    while ((residue = StreamCacheGetResidue (&sc)) != '\0' && j < sgp->numval) {
      if (IS_residue (residue)) {
        numBases++;
        /* val = (Int2) BSGetByte (bs); */
        if (k >= scount) {
          if (scount > 0) {
            scount = (Int2) BSRead (bs, scores, sizeof (scores));
          }
          k = 0;
        }
        if (scount > 0) {
          val = (Int2) scores [k];
          k++;
        } else {
          val = 0;
        }
        if (val < sgp->min.intvalue || val < 0) {
          valsBelowMin++;
        }
        if (val > sgp->max.intvalue || val > 100) {
          valsAboveMax++;
        }
        j++;
        switch (residue) {
        case '-': /* 0 */
          if (val > 0) {
            GapsWithScore++;
          }
          break;
        case 'A': /* 1, 2, 4, 8 */
        case 'C':
        case 'G':
        case 'T':
          if (val == 0) {
            ACGTsWithoutScore++;
            if (firstACGT == -1) {
              firstACGT = pos;
            }
          }
          break;
        case 'N': /* 15 */
          if (val > 0) {
            NsWithScore++;
            if (firstN == -1) {
              firstN = pos;
            }
          }
          break;
        default:
          if (val == 0) {
            ambigWithoutScore++;
            if (firstAmbig == -1) {
              firstAmbig = pos;
            }
          }
          break;
        }
      }
      pos++;
    }

    if (! single_report_mode) {
      gcp->itemID = sgp->idx.itemID;
      if (ACGTsWithoutScore > 0) {
        if (ACGTsWithoutScore * 10 >= numBases) {
          pct = (FloatHi) (ACGTsWithoutScore) * 100.0 / (FloatHi) numBases;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScoreMany, "%ld ACGT bases (%3.2f%s) have zero score value - first one at position %ld",
                    (long) ACGTsWithoutScore, (double) pct, "%", (long) (firstACGT + 1));
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScore, "%ld ACGT bases have zero score value - first one at position %ld",
                    (long) ACGTsWithoutScore, (long) (firstACGT + 1));
        }
      }
      if (NsWithScore > 0) {
        if (NsWithScore * 10 >= numBases) {
          pct = (FloatHi) (NsWithScore) * 100.0 / (FloatHi) numBases;
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScoreMany, "%ld N bases (%3.2f%s) have positive score value - first one at position %ld",
                    (long) NsWithScore, (double) pct, "%", (long) (firstN + 1));
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScore, "%ld N bases have positive score value - first one at position %ld",
                    (long) NsWithScore, (long) (firstN + 1));
        }
      }
      if (GapsWithScore > 0) {
        ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphGapScore, "%ld gap bases have positive score value", (long) GapsWithScore);
      }
      if (valsBelowMin > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBelow, "%ld quality scores have values below the reported minimum or 0", (long) valsBelowMin);
      }
      if (valsAboveMax > 0) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphAbove, "%ld quality scores have values above the reported maximum or 100", (long) valsAboveMax);
      }
    }
  }

  gcp->itemID = olditemid;
  gcp->thistype = olditemtype;

  if (single_report_mode) {
    if (ACGTsWithoutScore > 0) {
      if (ACGTsWithoutScore * 10 >= numBases) {
        pct = (FloatHi) (ACGTsWithoutScore) * 100.0 / (FloatHi) numBases;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScoreMany, "%ld ACGT bases (%3.2f%s) have zero score value - first one at position %ld",
                  (long) ACGTsWithoutScore, (double) pct, "%", (long) (firstACGT + 1));
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphACGTScore, "%ld ACGT bases have zero score value - first one at position %ld",
                  (long) ACGTsWithoutScore, (long) (firstACGT + 1));
      }
    }
    if (NsWithScore > 0) {
      if (NsWithScore * 10 >= numBases) {
        pct = (FloatHi) (NsWithScore) * 100.0 / (FloatHi) numBases;
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScoreMany, "%ld N bases (%3.2f%s) have positive score value - first one at position %ld",
                  (long) NsWithScore, (double) pct, "%", (long) (firstN + 1));
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphNScore, "%ld N bases have positive score value - first one at position %ld",
                  (long) NsWithScore, (long) (firstN + 1));
      }
    }
    if (GapsWithScore > 0) {
      ValidErr (vsp, SEV_ERROR, ERR_SEQ_GRAPH_GraphGapScore, "%ld gap bases have positive score value", (long) GapsWithScore);
    }
    if (valsBelowMin > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphBelow, "%ld quality scores have values below the reported minimum or 0", (long) valsBelowMin);
    }
    if (valsAboveMax > 0) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_GRAPH_GraphAbove, "%ld quality scores have values above the reported maximum or 100", (long) valsAboveMax);
    }
  }

  ValNodeFreeData (head);
}

/*****************************************************************************
*
*   PatchBadSequence(bsp)
*
*****************************************************************************/
NLM_EXTERN Boolean PatchBadSequence (BioseqPtr bsp)
{
  ByteStorePtr    newseq;
  SeqPortPtr      spp;
  Boolean         is_na;
  Uint1           seqcode;
  Int2            repchar, residue;
  Int4            i, len;

  if (bsp == NULL)
    return FALSE;
  if (!((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const)))
    return FALSE;

  is_na = ISA_na (bsp->mol);
  if (is_na) {
    seqcode = Seq_code_iupacna;
    repchar = (Int2) 'N';       /* N */
  } else {
    seqcode = Seq_code_iupacaa;
    repchar = (Int2) 'X';
  }

  spp = SeqPortNew (bsp, 0, -1, 0, seqcode);
  if (spp == NULL)
    return FALSE;

  len = bsp->length;
  newseq = BSNew (len);
  if (newseq == NULL) {
    SeqPortFree (spp);
    return FALSE;
  }

  for (i = 0; i < len; i++) {
    residue = SeqPortGetResidue (spp);
    if (residue == INVALID_RESIDUE) {
      residue = repchar;
    }
    BSPutByte (newseq, residue);
  }

  SeqPortFree (spp);
  SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data = (SeqDataPtr) newseq;
  bsp->seq_data_type = seqcode;

  BioseqRawPack (bsp);

  return TRUE;
}

static void FindABioseq (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr PNTR  bp;
  BioseqPtr       bsp;

  bp = (BioseqPtr PNTR) data;
  if (*bp != NULL)              /* already got one */
    return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) (sep->data.ptrvalue);
    *bp = bsp;
  }
  return;
}

NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf)
{
  BioseqPtr       bsp = NULL;

  if ((sep == NULL) || (buf == NULL))
    return NULL;

  *buf = '\0';
  SeqEntryExplore (sep, (Pointer) (&bsp), FindABioseq);

  if (bsp == NULL)
    return NULL;

  SeqIdPrint (bsp->id, buf, PRINTID_FASTA_LONG);
  return buf;
}

static CharPtr TrimSpacesOnEitherSide (CharPtr str)
{
  Uchar           ch;
  CharPtr         dst;
  CharPtr         ptr;

  if (str != NULL && str[0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static void CopyLetters (CharPtr dest, CharPtr source, size_t maxsize)
{
  Char            ch;
  CharPtr         tmp;

  if (dest == NULL || maxsize < 1)
    return;
  *dest = '\0';
  if (source == NULL)
    return;
  maxsize--;
  tmp = dest;
  ch = *source;
  while (maxsize > 1 && ch != '\0') {
    if (ch != '.') {
      *dest = ch;
      dest++;
      maxsize--;
    }
    source++;
    ch = *source;
  }
  *dest = '\0';
  TrimSpacesOnEitherSide (tmp);
}

static void LookForEtAl (ValidStructPtr vsp, ValNodePtr tmp)
{
  AuthorPtr       ap;
  AuthListPtr     authors = NULL;
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  CitSubPtr       csp;
  Char            first[64];
  Char            initials[16];
  Char            last[64];
  ValNodePtr      names;
  NameStdPtr      nsp;
  PersonIdPtr     pid;

  if (vsp == NULL || tmp == NULL)
    return;
  switch (tmp->choice) {
  case PUB_Article:
    cap = (CitArtPtr) (tmp->data.ptrvalue);
    authors = cap->authors;
    break;
  case PUB_Man:
  case PUB_Book:
  case PUB_Proc:
    cbp = (CitBookPtr) (tmp->data.ptrvalue);
    authors = cbp->authors;
    break;
  case PUB_Gen:
    cgp = (CitGenPtr) (tmp->data.ptrvalue);
    authors = cgp->authors;
    break;
  case PUB_Sub:
    csp = (CitSubPtr) (tmp->data.ptrvalue);
    authors = csp->authors;
    break;
  default:
    break;
  }
  if (authors == NULL || authors->choice != 1)
    return;
  for (names = authors->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL && pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL && nsp->names[0] != NULL) {
          CopyLetters (last, nsp->names[0], sizeof (last));
          CopyLetters (first, nsp->names[1], sizeof (first));
          CopyLetters (initials, nsp->names[4], sizeof (initials));
          if ((StringICmp (last, "et al") == 0) || (StringCmp (initials, "al") == 0 && StringCmp (last, "et") == 0 && first[0] == '\0')) {
            if (names->next == NULL) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl, "Author list ends in et al.");
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl, "Author list contains et al.");
            }
          }
        }
      }
    }
  }
}

static void SpellCheckPub (ValidStructPtr vsp, ValNodePtr tmp)
{
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  ValNodePtr      titles = NULL;

  if ((vsp == NULL) || (tmp == NULL))
    return;

  switch (tmp->choice) {
  case PUB_Article:
    cap = (CitArtPtr) (tmp->data.ptrvalue);
    titles = cap->title;
    break;
  case PUB_Man:
  case PUB_Book:
  case PUB_Proc:
    cbp = (CitBookPtr) (tmp->data.ptrvalue);
    titles = cbp->title;
    break;
  case PUB_Gen:
    cgp = (CitGenPtr) (tmp->data.ptrvalue);
    if (cgp->cit != NULL)
      SpellCheckString (vsp, cgp->cit);
    if (cgp->title != NULL)
      SpellCheckString (vsp, cgp->title);
    break;
  default:
    break;
  }

  if (titles != NULL) {
    for (; titles != NULL; titles = titles->next) {
      if (titles->choice == Cit_title_name)
        SpellCheckString (vsp, (CharPtr) (titles->data.ptrvalue));
    }
  }

  return;
}

static void SpellCheckSeqDescr (GatherContextPtr gcp)
{
  PubdescPtr      pdp;
  ValNodePtr      tmp, vnp;
  ValidStructPtr  vsp;

  vsp = (ValidStructPtr) (gcp->userdata);
  if (vsp == NULL)
    return;

  vnp = (ValNodePtr) (gcp->thisitem);
  if (vnp == NULL)
    return;

  vsp->descr = vnp;
  vsp->sfp = NULL;

  if (vnp->choice == Seq_descr_pub) {
    pdp = (PubdescPtr) (vnp->data.ptrvalue);
    for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next) {
      LookForEtAl (vsp, tmp);
    }
  }

  if (vsp->spellfunc == NULL)
    return;

  switch (vnp->choice) {
  case Seq_descr_title:
  case Seq_descr_region:
  case Seq_descr_comment:
    SpellCheckString (vsp, (CharPtr) (vnp->data.ptrvalue));
    break;
  case Seq_descr_pub:
    pdp = (PubdescPtr) (vnp->data.ptrvalue);
    for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next) {
      SpellCheckPub (vsp, tmp);
    }
    SpellCheckString (vsp, pdp->comment);
    break;
  default:
    break;
  }
  return;
}

NLM_EXTERN void SpellCheckSeqFeat (GatherContextPtr gcp)
{
  PubdescPtr      pdp;
  SeqFeatPtr      sfp;
  ProtRefPtr      prp;
  ValidStructPtr  vsp;
  ValNodePtr      vnp;

  vsp = (ValidStructPtr) (gcp->userdata);
  if (vsp == NULL)
    return;

  sfp = (SeqFeatPtr) (gcp->thisitem);
  if (sfp == NULL)
    return;

  vsp->descr = NULL;
  vsp->sfp = sfp;

  if (sfp->data.choice == SEQFEAT_PUB) {
    pdp = (PubdescPtr) (sfp->data.value.ptrvalue);
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      LookForEtAl (vsp, vnp);
    }
  }

  if (vsp->spellfunc == NULL)
    return;

  SpellCheckString (vsp, sfp->comment);

  switch (sfp->data.choice) {
  case 1:                      /* Gene-ref */
    break;
  case 2:                      /* Org-ref */
    break;
  case 3:                      /* Cdregion */
    break;
  case 4:                      /* Prot-ref */
    prp = (ProtRefPtr) (sfp->data.value.ptrvalue);
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      SpellCheckString (vsp, (CharPtr) (vnp->data.ptrvalue));
    SpellCheckString (vsp, prp->desc);
    break;
  case 5:                      /* RNA-ref */
    break;
  case 6:                      /* Pub */
    pdp = (PubdescPtr) (sfp->data.value.ptrvalue);
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      SpellCheckPub (vsp, vnp);
    }
    SpellCheckString (vsp, pdp->comment);
    break;
  case 7:                      /* Seq */
    break;
  case 8:                      /* Imp-feat */
    break;
  case 9:                      /* Region */
    SpellCheckString (vsp, (CharPtr) (sfp->data.value.ptrvalue));
    break;
  case 10:                     /* Comment */
    break;
  case 11:                     /* Bond */
    break;
  case 12:                     /* Site */
    break;
  case 13:                     /* Rsite-ref */
    break;
  case 14:                     /* User-object */
    break;
  case 15:                     /* TxInit */
    break;
  case 16:                     /* Numbering */
    break;
  case 17:                     /* Secondary Structure */
    break;
  case 18:                     /* NonStdRes */
    break;
  case 19:                     /* Heterogen */
    break;
  case 20:                     /* BioSource */
    break;
  default:
    break;
  }

  return;
}

NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str)
{
  if ((vsp == NULL) || (str == NULL))
    return;

  if (vsp->spellfunc == NULL)
    return;

  (*(vsp->spellfunc)) ((char *) str, (vsp->spellcallback));

  return;
}

NLM_EXTERN void SpellCallBack (char *str)
{
  ErrSev          sev;

  sev = SEV_ERROR;
  if (globalvsp != NULL && globalvsp->justwarnonspell) {
    sev = SEV_WARNING;
  }
  ValidErr (globalvsp, sev, ERR_GENERIC_Spell, "[ %s ]", (CharPtr) str);
  return;
}
