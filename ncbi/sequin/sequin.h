/*   sequin.h
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
* File Name:  sequin.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.565 $
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

#ifndef _SEQUIN_
#define _SEQUIN_

#ifdef INTERNAL_NCBI_SEQUIN
#ifndef EXTRA_SERVICES
#define EXTRA_SERVICES
#endif
#ifndef NETWORK_SAVVY_SEQUIN
#define NETWORK_SAVVY_SEQUIN
#endif
#ifndef USE_SPELL
#define USE_SPELL
#endif
#endif

#ifdef EXTRA_SERVICES
#define USE_DESKTOP
#define REPLACE_THIS
#define EDIT_LOCUS
#endif

#ifdef NETWORK_SAVVY_SEQUIN
#define USE_ENTREZ
#define USE_LOCAL
#define USE_BLAST
#define USE_MEDARCH
#define USE_TAXON
#define ALLOW_DOWNLOAD
#endif

#ifdef PUBLIC_NETWORK_SEQUIN
#define USE_DESKTOP
#define USE_ENTREZ
#define USE_LOCAL
#define USE_BLAST
#define USE_MEDARCH
#define USE_TAXON
#define ALLOW_DOWNLOAD
#endif

#include <dlogutil.h>
#include <bspview.h>
#include <objproj.h>
#include <urlquery.h>


#ifdef __cplusplus
extern "C" {
#endif


#define SEQ_PKG_SINGLE        1
#define SEQ_PKG_SEGMENTED     2
#define SEQ_PKG_GAPPED        3
#define SEQ_PKG_GENOMICCDNA   4
#define SEQ_PKG_POPULATION    5
#define SEQ_PKG_PHYLOGENETIC  6
#define SEQ_PKG_MUTATION      7
#define SEQ_PKG_ENVIRONMENT   8
#define SEQ_PKG_GENBANK       9
#define NUM_SEQ_PKG           9

#define SEQ_FMT_FASTA         1
#define SEQ_FMT_ALIGNMENT     2 
#define NUM_SEQ_FMT           2 
  
/*
#define SEQ_FMT_FASTAGAP      2
#define SEQ_FMT_PHYLIP        3
#define SEQ_FMT_NEXUS         4
#define SEQ_FMT_PAUP          5
*/

#define SEQ_ORIG_SUBMISSION   1
#define SEQ_TPA_SUBMISSION    2

typedef struct fmtblk {
  Int2         seqPackage;
  Int2         seqFormat;
  Int2         numSeqs;
  Int2         submType;
} FormatBlock, PNTR FormatBlockPtr;

typedef struct sqnblk {
  AuthorPtr    contactperson;
  AuthListPtr  citsubauthors;
  AffilPtr     citsubaffil;
  CharPtr      citsubtitle;
  DatePtr      releasedate;
  Boolean      holduntilpublished;
} SequinBlock, PNTR SequinBlockPtr;

extern CharPtr SEQUIN_APPLICATION;
extern CharPtr SEQUIN_SERVICES;
extern CharPtr SEQUIN_VERSION;

extern ForM  helpForm;

extern Boolean  useDesktop;
extern Boolean  useEntrez;
extern Boolean  useLocal;
extern Boolean  useBlast;
extern Boolean  useMedarch;
extern Boolean  useTaxon;
extern Boolean  allowDownload;
extern Boolean  extraServices;
extern Boolean  indexerVersion;
extern CharPtr  genomeCenter;

extern Boolean  leaveAsOldAsn;
extern Boolean  newAlignReader;

#ifdef WIN_MAC
extern Boolean  termListUp;
extern Boolean  docSumUp;
extern Boolean  bioseqViewUp;
#endif

extern void SwapQualifiers (IteM i);
extern void PrefixAuthorityWithOrganism (IteM i);
extern void UpdateFastaSet (IteM i);
extern void ExtendFastaSet (IteM i);
extern void ExtendAllSequencesInSet (IteM i);
extern void SeqLocAdjustByOffset (SeqLocPtr slp, Int4 offset);
extern void SplitSegmentedFeatsMenuItem (IteM i);
extern SeqFeatPtr SeqFeatCopy (SeqFeatPtr sfp);
extern SeqLocPtr SeqLocReplaceLocalID (SeqLocPtr slp,
				       SeqIdPtr  new_sip);
extern SequinBlockPtr SequinBlockFree (SequinBlockPtr sbp);

extern ForM CreateStartupForm (Int2 left, Int2 top, CharPtr title,
                               BtnActnProc startFa2htgs,
                               BtnActnProc startPhrap,
                               BtnActnProc buildContig,
                               BtnActnProc startNew,
                               BtnActnProc readExisting,
                               BtnActnProc fetchFromNet,
                               BtnActnProc showHelp,
                               BtnActnProc createSubmissionTemplate,
                               BtnActnProc quitProgram,
                               WndActnProc activateForm);

extern void DrawAbout (PaneL p);
extern Int2 AboutBoxWidth (void);
extern Int2 AboutBoxHeight (void);

extern ForM CreateFormatForm (Int2 left, Int2 top, CharPtr title,
                              BtnActnProc goToNext,
                              BtnActnProc goBack,
                              WndActnProc activateForm);

extern ForM CreateInitSubmitterForm (Int2 left, Int2 top, CharPtr title,
                                     BtnActnProc goToNext,
                                     BtnActnProc goBack,
                                     WndActnProc activateForm);

extern DialoG CreateFastaDialog (GrouP h, CharPtr title, Boolean is_na, Boolean is_mrna,
                                 CharPtr text, Boolean single, Int2Ptr seqPackagePtr);

extern ForM CreateInitOrgNucProtForm (Int2 left, Int2 top, CharPtr title,
                                      FormatBlockPtr format,
                                      BtnActnProc goToNext,
                                      BtnActnProc goBack,
                                      WndActnProc activateForm);

extern ForM CreateGenomeCenterForm (Int2 left, Int2 top, CharPtr title,
                                    BtnActnProc finish,
                                    BtnActnProc cancel,
                                    Boolean readPhrap,
                                    Boolean buildContig,
                                    WndActnProc activateForm);

extern SeqEntryPtr ImportOneGappedSequence (FILE *fp);

extern Boolean HasZeroLengthSequence (ForM newForm);
extern Boolean SequencesFormHasProteins (ForM f);
extern SeqEntryPtr GetSequencesFormProteinList (ForM f);
extern SeqEntryPtr GetSequencesFormNucleotideList (ForM f);
extern Boolean SequencesFormHasTooManyNucleotides (ForM f);

extern void AppendOrReplaceString (
  CharPtr PNTR string_loc,
  CharPtr new_value,
  Boolean PNTR asked_question,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon
);

extern void ConsolidateOrganismNotes (IteM i);
extern void ConsolidateLikeModifiersWithSemicolons (IteM i);
extern void ConsolidateLikeModifiersWithoutSemicolons (IteM i);

extern void CountryLookupWithoutCapFix (IteM i);
extern void CountryLookupWithCapFix (IteM i);
extern void ExtractProteinFeaturesFromNote (IteM i);
extern void ConvertPseudoCDSToMiscFeat (IteM i);
extern void ProcessPseudoMiscFeat (IteM i);
extern void ParseInfluenzaAVirusNames (IteM i);
extern void AddStrainAndSerotypeToInfluenzaAVirusNames (IteM i);
extern void FixupInfluenzaAVirusNames(IteM i);
extern void EditPubs (IteM i);
extern void EditPubsEx (BaseFormPtr bfp);
extern void RemovePubConsortiums (IteM i);

extern void ExtendPartialFeatures (IteM i);
extern void ExtendPartialFeaturesWithConstraint (IteM i);
extern void TrimOrganismName (IteM i);
extern void SUCSubmitterProc (IteM i);

extern CharPtr FixInfluenzaVirusName (CharPtr orig_name);

extern void ConfirmSequencesFormParsing (ForM f, FormActnFunc putItAllTogether);

extern ForM CreateHelpForm (Int2 left, Int2 top, CharPtr title,
                            CharPtr file, BtnActnProc closeForm,
                            WndActnProc activateForm);

extern void SendHelpScrollMessage (ForM f, CharPtr heading, CharPtr section);

extern void ApplyCDSFrame (IteM i);

/* The next pointer in NewObject is not used in freeing the list.  Each
block is attached individually as extra data to the appropriate menu item.
The linked list is used solely to enable and disable new feature menu items
by the target bsp->mol, or to enable and disable analysis menu items by the
ability to produce FASTA (bioseq viewer or docsum window). */

typedef struct urlparamdata {
  Uint1          type;     /* 1 = text, 2 = checkbox, 3 = popup, 4 = radio, 5 = list */
  CharPtr        param;
  CharPtr        prompt;   /* if no prompt, use param */
  CharPtr        dfault;
  CharPtr        choices;  /* choices if param is popup */
  CharPtr        group;    /* used for grouping related controls */
  CharPtr        descr;
  CharPtr        help;
} UrlParamData, PNTR UrlParamPtr;

typedef struct newobjectdata {
  Int2             kind;   /* 1 = feature creation, 2 = analysis */
  ObjMgrProcPtr    ompp;
  BaseFormPtr      bfp;
  IteM             item;
  Uint1            molgroup;
  Uint2            descsubtype;
  Boolean          bspOK;
  Boolean          dsmOK;
  Boolean          fastaNucOK;
  Boolean          fastaProtOK;
  Boolean          onlyBspTarget;
  /* the next eight fields are for the analysis menu only, for remote URLs */
  CharPtr          host_machine;
  Uint2            host_port;
  CharPtr          host_path;
  CharPtr          query;
  Uint4            timeoutsec;
  Int2             format;     /* 1 = FASTA, 2 = ASN.1 */
  Boolean          demomode;
  QueryResultProc  resultproc;
  ValNodePtr       paramlist; /* data.ptrvalue points to UrlParamData block */
  CharPtr          prefix;
  CharPtr          suffix;
  CharPtr          homepage;
  CharPtr          credits;
  CharPtr          authors;
  CharPtr          disclaimer;
  CharPtr          reference;
  Uint4            pmid;
  CharPtr          blurb;
  struct newobjectdata PNTR next;
} NewObjectData, PNTR NewObjectPtr;

#ifdef WIN_MAC
extern VoidPtr macUserDataPtr;
#endif

extern void SetupSpecialMenu (MenU m, BaseFormPtr bfp);
extern void SetupNewFeaturesMenu (MenU m, BaseFormPtr bfp);
extern void SetupNewDescriptorsMenu (MenU m, BaseFormPtr bfp);
extern void SetupNewPublicationsMenu (MenU m, BaseFormPtr bfp);
extern void SetupBatchApplyMenu (MenU s, BaseFormPtr bfp);
extern void SetupBatchEditMenu (MenU s, BaseFormPtr bfp);
extern MenU CreateAnalysisMenu (WindoW w, BaseFormPtr bfp, Boolean bspviewOK, Boolean docsumOK);
extern void SetupSequinFilters (void);
extern void SetupBioseqPageList (void);

extern Boolean LIBCALLBACK SequinOpenMimeFile (CharPtr filename);
extern Boolean LIBCALLBACK SequinOpenResultFile (CharPtr filename);
extern Boolean LIBCALLBACK SequinHandleNetResults (CharPtr filename);

extern void SequinCheckSocketsProc (void);

extern Int4 MySeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force);
extern void ValSeqEntryForm (ForM f);

extern void InitSequinExtras (void);
extern void FiniSequinExtras (void);

/* This function destroys the SequinBlockPtr */

extern Uint2 PackageFormResults (SequinBlockPtr sbp, SeqEntryPtr sep,
                                 Boolean makePubAndDefLine);

extern void EnableFeaturesPerTarget (BaseFormPtr bfp);
extern void EnableAnalysisItems (BaseFormPtr bfp, Boolean isDocSum);

extern void ExtendSeqLocToPosition (SeqLocPtr slp, Boolean end5, Int4 pos);

#define REGISTER_BIOSEQ_SEG_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit Bioseq Seg","BioseqSegEditor",OBJ_BIOSEQ_SEG,0,OBJ_BIOSEQ_SEG,0,NULL,BioseqSegEditFunc,PROC_PRIORITY_DEFAULT)
extern Int2 LIBCALLBACK BioseqSegEditFunc (Pointer data);

#define REGISTER_BIOSEQ_SET_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit Bioseq Set","BioseqSetEditor",OBJ_BIOSEQSET,0,OBJ_BIOSEQSET,0,NULL,BioseqSetEditFunc,PROC_PRIORITY_DEFAULT)
extern Int2 LIBCALLBACK BioseqSetEditFunc (Pointer data);

extern void LaunchOrfViewer (BioseqPtr bsp, Uint2 entityID, Uint4 itemID, Boolean standAlone);

extern Int2 ApplyAnnotationToAll (Int2 type, SeqEntryPtr sep,
                                  ButtoN partialLft, ButtoN partialRgt,
                                  TexT geneName, TexT protName, 
                                  TexT protDesc, TexT rnaName,
                                  TexT featcomment, TexT defline);

extern SeqFeatPtr FindBestCds (Uint2 entityID, SeqLocPtr loc, SeqLocPtr prod, SeqEntryPtr scope);

NLM_EXTERN SeqEntryPtr SequinFastaToSeqEntryEx 
  (
    FILE *fp, Boolean is_na, CharPtr PNTR errormsg,
    Boolean parseSeqId, CharPtr special_symbol
  );

NLM_EXTERN SeqEntryPtr SequinFastaToSeqEntryExEx
  (
    FILE *fp,               /* file to get sequence from */ 
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr special_symbol, /* Returns special symbol if no SeqEntry */
    BoolPtr chars_stripped  /* set to TRUE if characters other than digits
                             * were stripped from the FASTA sequence data */
  );

/* Many miscellaneous extern functions within sequin source files */

extern Boolean WriteSequinAppParam (CharPtr section, CharPtr type, CharPtr value);

extern Boolean PropagateFromGenBankBioseqSet (SeqEntryPtr sep, Boolean ask);

extern Uint2 SmartAttachSeqAnnotToSeqEntry (Uint2 entityID, SeqAnnotPtr sap, ValNodePtr PNTR err_list);
extern void HandleProjectAsn (ProjectPtr proj, Uint2 entityID);

extern CharPtr CompressSpaces (CharPtr str);
extern CharPtr SearchForString (CharPtr str, CharPtr sub, Boolean case_counts, Boolean whole_word);
extern void AddAboutAndHelpMenuItems (MenU m);
extern void NetConfigureProc (IteM i);
extern void EntrezQueryProc (IteM i);
extern void Entrez2QueryProc (IteM i);
extern void SetupEditSecondary (MenU m, BaseFormPtr bfp);
extern void SimpleCDDBlastProc (IteM i);
extern void SimpleCDDSearchFeatProc (IteM i);
extern void SimpleCDDSearchAlignProc (IteM i);
extern void ForceCleanupEntityID (Uint2 entityID);
extern void ForceTaxonFixupBtn (IteM i, ButtoN b);
extern void CommonAddOrgOrModsToDefLines (IteM i, Int2 orgmod, Int2 subsource, ButtoN b);
extern void PrefixDefLines (IteM i);
extern void MRnaFromCdsProc (Uint2 entityID);
extern void BioseqViewFormToolBar (GrouP h);
extern Boolean DoBuildContig (void);
extern void SetGenome (PopuP p);
extern PopuP ReplaceBioSourceGencodePopup (DialoG d, PopuP gencode);
extern CharPtr NameStdPtrToAuthorSpreadsheetString (NameStdPtr nsp);
extern NameStdPtr AuthorSpreadsheetStringToNameStdPtr (CharPtr txt);
extern Boolean ExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp);
extern void CommonAddOrgOrModsToDefLines (IteM i, Int2 orgmod, Int2 subsource, ButtoN b);
extern void CorrectGenCodes (SeqEntryPtr sep, Uint2 entityID);
extern void PrepareToConvertToCDS (SeqEntryPtr sep, Uint2 entityID,
                                   Uint2 subtype, CharPtr findthis);
extern void EditGenbankElements (Handle i);
extern void EditSequenceHistory (IteM i);
extern void InsertGeneLocusTagPrefix (IteM i);
extern void ReplaceRepeatRegionLocusTagWithDbxref (IteM i);
extern void FindGeneAndProtForCDS (Uint2 entityID, SeqFeatPtr cds,
                                   SeqFeatPtr PNTR gene, SeqFeatPtr PNTR prot);
extern void ExportAlignmentInterleave (IteM i);
extern void ExportAlignmentContiguous (IteM i);
extern void FixFeatureIntervals (IteM i);
extern void ConvertInnerCDSsToProteinFeatures (IteM i);
extern void CombineMultipleCDS (IteM i);

extern void NewDescriptorMenuFunc (ObjMgrProcPtr ompp, BaseFormPtr bfp, Uint2 descsubtype);
extern Boolean PropagateFromGenBankBioseqSet (SeqEntryPtr sep, Boolean ask);
extern int LIBCALLBACK SortByVnpChoice (VoidPtr ptr1, VoidPtr ptr2);
extern void PrepareToConvertToCDS (SeqEntryPtr sep, Uint2 entityID,
                                   Uint2 subtype, CharPtr findthis);
extern void ConvertToLocalProcOnlyNucs (IteM i);
extern void ConvertToLocalProcOnlyProts (IteM i);
extern void ConvertToLocalProcAll (IteM i);

extern void PromoteToBestIDProc (IteM i);
extern void PromoteToWorstIDProc (IteM i);
extern void ChangeGenBankNameToLocal (IteM i);
extern void RemoveGBIDsFromBioseqs (IteM i);
extern void RemoveGBIDsFromProteins (IteM i);
extern void RemoveGIsFromBioseqs (IteM i);

extern void CommonApplyToAllProc (BaseFormPtr bfp, Int2 type);
extern void ApplyTitle (IteM i);
extern void ApplyCDS (IteM i);
extern void ApplyRRNA (IteM i);
extern void ApplyImpFeat (IteM i);
extern void AdjustCDSLocationsForKnownAndUnknownGapsCallback (SeqFeatPtr sfp, Pointer userdata);
extern void AdjustFeaturesForGaps (IteM i);
extern void LoadTPAAccessionNumbersFromFile (IteM i);
extern void LoadSecondaryAccessionNumbersFromFile (IteM i);
extern void LoadHistoryAccessionNumbersFromFile (IteM i);
extern void LoadOrganismModifierTable (IteM i);
extern void LoadTaxConsult (IteM i);
extern void ExportOrganismTable (IteM i);
extern void LoadFeatureQualifierTable (IteM i);

extern void AddCodonListTotRNA (tRNAPtr trna, ValNodePtr codons);

extern void RemoveRedundantProproteinMiscFeats (IteM i);
extern void AddTypeStrainCommentsToAll (IteM i);
extern void AddTypeStrainCommentsWithConstraint (IteM i);
extern void RemoveSequencesFromAlignment (IteM i);
extern void RemoveSequencesFromRecord (IteM i);

extern void ParseFileToSource (IteM i);
extern void AddModToOrg (IteM i);

extern void ParseInMoreProteins (IteM i);
extern void ParseInNucUpdates (IteM i);
extern void ParseInOligoPrimers (IteM i);
extern void ParseInMoreMRNAs (IteM i);

extern void RecomputeSuggestEx (Uint2 entityID, Boolean fix_genes, Boolean recompute_all);
extern void RecomputeSuggest (IteM i);
extern void RecomputeSuggestFixGenes (IteM i);
extern void RetranslateCdRegionsEx (
  Uint2   entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds );
extern void RetranslateCdRegionsNoStop (IteM i);
extern void RetranslateCdRegionsDoStop (IteM i);
extern void RetranslateCdRegionsNoStopExceptEndCompleteCDS (IteM i);
extern void AddGlobalCodeBreak (IteM i);
extern void ParseCodonQualToCodeBreak (IteM i);
extern void CorrectCDSGenCodes (IteM i);
/* extern void CorrectCDSStartCodon (IteM i); */
/* extern Boolean RetranslateOneCDS (SeqFeatPtr sfp, Uint2 entityID, Boolean include_stop); */
extern void UpdateProteinsFromCDS (IteM i);

extern void AutoDef (IteM i);
extern void AutoDefWithOptions (IteM i);
extern void AutoDefWithoutModifiers (IteM i);
extern void AutoDefBaseFormCommon (BaseFormPtr bfp, Boolean use_form, Boolean use_modifiers);
extern void AutoDefStrain (BaseFormPtr bfp);
extern void AutoDefMiscFeat (BaseFormPtr bfp);
extern void AutoDefToolBtn (ButtoN b);
extern void AutoDefOptionsToolBtn (ButtoN b);
extern void AutoDefStrainToolBtn (ButtoN b);
extern void AutoDefMiscFeatToolBtn (ButtoN b);
extern void AutoDefEntityIDNoOptions (Uint2 entityID, Boolean use_modifiers);

extern void RemoveDefLinesToolBtn (ButtoN b);
extern void FindStringProcToolBtn (ButtoN b);
extern void FindFlatfileProcToolBtn (ButtoN b);
extern void ResolveExistingLocalIDsToolBtn (ButtoN b);
extern void GroupExplodeToolBtn (ButtoN b);

extern Int2 LIBCALLBACK MakeGroupsOf200 (Pointer data);

extern void SetBestFrame (SeqFeatPtr sfp);
extern Boolean SetBestFrameByLocation (SeqFeatPtr sfp);

extern void ViewAlignmentSummary (IteM i);

extern void SetupEditSecondary (MenU m, BaseFormPtr bfp);
extern void EditLocusProc (IteM i);

extern ValNodePtr BuildDescriptorValNodeList (void);

extern void RemoveDescriptor (IteM i);

extern void SelectDescriptor (IteM i);
extern void SelectBioseq (IteM i);
extern void SelectPubs (IteM i);

extern void FuseFeature (IteM i);

extern void MakeExonsFromCDSIntervals (IteM i);
extern void MakeExonsFromMRNAIntervals (IteM i);

extern Int2 LIBCALLBACK CreateDeleteByTextWindow (Pointer data);
extern Int2 LIBCALLBACK CreateSegregateByTextWindow (Pointer data);
extern Int2 LIBCALLBACK SegregateSetsByField (Pointer data);
extern Int2 LIBCALLBACK CreateSegregateByFeatureWindow (Pointer data);
extern Int2 LIBCALLBACK CreateSegregateByDescriptorWindow (Pointer data);
extern Int2 LIBCALLBACK CreateSegregateByMoleculeTypeWindow (Pointer data);
extern Int2 LIBCALLBACK CreateSegregateByIdWindow (Pointer data);
extern Int2 LIBCALLBACK SequesterSequences (Pointer data);
extern Int2 LIBCALLBACK RemoveExtraneousSets (Pointer data);
extern void ReverseComplementBioseqAndFeats (BioseqPtr bsp, Uint2 entityID);
extern void RemoveOrphanProteins (Uint2 entityID, SeqEntryPtr sep);
extern void RemoveTextInsideString (IteM i);
extern void RemoveTextOutsideString (IteM i);

extern void BioseqViewFormToolBar (GrouP h);

extern void FindStringProc (IteM i);
extern void FindFlatfileProc (IteM i);
extern void FindGeneProc (IteM i);
extern void FindProtProc (IteM i);
extern void FindPosProc (IteM i);

extern void SimpleUniVecScreenProc (IteM i);
extern void SimpleUniVecCoreScreenProc (IteM i);

extern Boolean MeetsStringConstraint (SeqFeatPtr sfp, CharPtr str, Boolean case_insensitive);

extern Boolean SaveSeqSubmitProc (BaseFormPtr bfp, Boolean saveAs);

extern void ExciseString (CharPtr str, CharPtr from, CharPtr to);
extern void MakeSearchStringFromAlist (CharPtr str, CharPtr name);
extern void AddToSubSource (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype);
extern void AddToOrgMod (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype);
extern Boolean AutomaticProteinProcess (SeqEntryPtr esep, SeqEntryPtr psep,
                                        Int2 code, Boolean makeMRNA, 
                                        SeqLocPtr use_this);

extern CharPtr repackageMsg;
extern BioseqPtr  updateTargetBspKludge;
extern SeqEntryPtr     globalsep;
extern Uint2           globalEntityID;
extern Char            globalPath [PATH_MAX];
extern ForM  startupForm;
extern SeqViewProcs        seqviewprocs;

extern void CommonFetchFromNet (BtnActnProc actn, BtnActnProc cancel);
extern void FetchFromNet (ButtoN b);
extern Boolean SequinEntrezInit (CharPtr appl_id, Boolean no_warnings, BoolPtr is_network);
extern void JustRegisterSeqEntry (BaseFormPtr bfp, Boolean freeit);
extern void JustRegisterSeqEntryBtn (ButtoN b);
extern void AddSubmitBlockToSeqEntry (ForM f);

extern void SqnReadAlignView (BaseFormPtr bfp, BioseqPtr target_bsp, SeqEntryPtr source_sep, Boolean do_update);
extern void DownloadAndUpdateProc (ButtoN b);
extern void DownloadAndExtendProc (ButtoN b);
extern void UpdateSeqAfterDownload (BaseFormPtr bfp, BioseqPtr oldbsp, BioseqPtr newbsp);
extern void ExtendSeqAfterDownload (BaseFormPtr bfp, BioseqPtr oldbsp, BioseqPtr newbsp);
extern void NewUpdateSequence (IteM i);
extern void NewExtendSequence (IteM i);

extern void FastaNucDirectToSeqEdProc (IteM i);

extern void ParseCodonsFromtRNAComment (IteM i);
extern void ParseAntiCodonsFromtRNAComment (IteM i);

extern void RemoveAlignment (IteM i);
extern void RemoveGraph (IteM i);
extern void RemoveSeqAnnotIDs (IteM i);
extern void RemoveSeqAnnotLOCs (IteM i);

extern void RemoveProteins (IteM i);
extern void RemoveProteinsAndRenormalize (IteM i);

extern void GlobalAddTranslExcept (IteM i);
extern void AddTranslExceptWithComment (IteM i);

extern const char *nucleotide_alphabet;
extern const char *protein_alphabet;
extern void ReadAlignment (IteM i);
extern SeqEntryPtr SeqEntryFromAlignmentFile (FILE *fp, TSequenceInfoPtr sequence_info, Uint1 moltype,
                                              CharPtr no_org_err_msg);

extern SeqAlignPtr Sqn_GlobalAlignTwoSeq (BioseqPtr bsp1, BioseqPtr bsp2, BoolPtr revcomp);

extern void SqnNewAlign (BioseqPtr bsp1, BioseqPtr bsp2, SeqAlignPtr PNTR salp);

extern void ProduceAlignmentNotes (TAlignmentFilePtr afp, TErrorInfoPtr error_list);

extern void RemoveAlignmentsWithSequence (BioseqPtr bsp, Uint2 input_entityID);
extern void FlipEntireAlignmentIfAllSequencesFlipped (SeqAnnotPtr sap, Pointer userdata);

#ifndef WIN_MAC
extern void CreateSqnInitialFormMenus (WindoW w);
#endif

#define NUM_PAGES  8

typedef struct nucprotassoc {
  Int4 position;
  SeqLocPtr loc;
  struct nucprotassoc PNTR  next;
} NucProtAssocData, PNTR NucProtAssocPtr;

typedef struct sequencesform {
  FORM_MESSAGE_BLOCK
  GrouP           pages [NUM_PAGES];
  Int2            currentPage;
  Int2            tagFromPage [NUM_PAGES];
  Int2            numPages;
  DialoG          tbs;

  Uint1           dnamolfrommolinfo;
  EnumFieldAssoc  PNTR moltypeAlist;
  ButtoN          makeAlign;
  DialoG          dnaseq;

  Int2            seqPackage;
  Int2            seqFormat;
  Int2            numSeqs;
  Int2            submType;

  ButtoN          protTechBoth;
  ButtoN          partialN;
  ButtoN          partialC;
  Boolean         makeMRNA;
  DialoG          protseq;

  DialoG          mrnaseq;
  ButtoN          partialmRNA5;
  ButtoN          partialmRNA3;

  GrouP           annotType;
  GrouP           annotGrp;
  ButtoN          partialLft;
  ButtoN          partialRgt;
  TexT            geneName;
  PrompT          protOrRnaPpt;
  TexT            protOrRnaName;
  PrompT          protDescPpt;
  TexT            protDesc;
  TexT            featcomment;
  TexT            defline;
  ButtoN          orgPrefix;

  ButtoN          nextBtn;
  ButtoN          prevBtn;
  BtnActnProc     goToNext;
  BtnActnProc     goToPrev;

  /* These are added to add modifiers on the source tab */  
  ButtoN          import_mod_btn;
  ButtoN          source_assist_btn;
  ButtoN          specify_orgs_btn;
  ButtoN          specify_locs_btn;
  ButtoN          specify_gcode_btn;
  ButtoN          specify_mgcode_btn;
  ButtoN          clear_mods_btn;
  DoC             org_doc;
  GrouP           ident_org_grp;
  DialoG          summary_dlg;
  
  /* These allow the user to specify topology and molecule */
  ButtoN          topology_btn;
  ButtoN          molecule_btn;
  
  /* This list pairs the proteins and nucleotides. */
  /* It must be freed using FreeAssociationList. */
  NucProtAssocPtr nuc_prot_assoc_list;

} SequencesForm, PNTR SequencesFormPtr;

extern ValNodePtr InsertMostUsedFeatureValNodes (ValNodePtr old_list);

extern ValNodePtr FindExactStringInStrings ( ValNodePtr strings, CharPtr value);

extern EnumFieldAssocPtr InsertMostUsedFeatureEnumFieldAssoc (
  EnumFieldAssocPtr alist
);

extern ValNodePtr BuildFeatureValNodeList (
  Boolean prefer_most_used,
  CharPtr wild_card_name,
  Int4 wild_card_value,
  Boolean skip_unusual,
  Boolean skip_import
);

extern void RemoveOldName (OrgRefPtr orp);
extern void SetTaxNameAndRemoveTaxRef (OrgRefPtr orp, CharPtr taxname);

extern void MergeToPartsJoin (IteM i);
extern void MergeToPartsOrdered (IteM i);

extern void ConvertInnerCDSsToMatPeptidesCallback (BioseqPtr bsp, Pointer userdata);
extern void MergeCDS (IteM i);

extern void InitValNodePopup (ValNodePtr list, PopuP p);
extern Int2 GetValNodePopup (PopuP p, ValNodePtr list);
extern void SetValNodePopupValue (ValNodePtr list, PopuP p, CharPtr val);

extern Uint1 FindTypeForModNameText (CharPtr cp);

typedef struct featureswithtextdata 
{
  Uint1             seqFeatChoice;
  Uint1             featDefChoice;
  CharPtr           search_text;
  Boolean           case_insensitive;
  Boolean           whole_word;
  Boolean           no_text;
  Boolean           act_when_string_not_present;
  VisitFeaturesFunc callback;
  Pointer           userdata;
} FeaturesWithTextData, PNTR FeaturesWithTextPtr;

typedef struct descriptorswithtextdata 
{
  CharPtr           search_text;
  Boolean           case_insensitive;
  Boolean           whole_word;
  Boolean           no_text;
  Boolean           act_when_string_not_present;
  VisitDescriptorsFunc callback;
  Pointer           userdata;
} DescriptorsWithTextData, PNTR DescriptorsWithTextPtr;


extern void OperateOnBioseqFeaturesWithText 
(BioseqPtr         bsp,
 Pointer           userdata);

extern void OperateOnSeqEntryFeaturesWithText (SeqEntryPtr sep, FeaturesWithTextPtr fdp);
extern void OperateOnSeqEntryDescriptorsWithText (SeqEntryPtr sep, DescriptorsWithTextPtr ddp);

extern LisT 
MakeSequenceListControl 
(GrouP g,
 SeqEntryPtr sep,
 Nlm_LstActnProc actn, 
 Pointer userdata, 
 Boolean show_nucs, 
 Boolean show_prots);
extern ValNodePtr GetSelectedSequenceList (LisT l);
extern void SelectAllSequencesInListCtrl (LisT l);
extern void UnSelectAllSequencesInListCtrl (LisT l);
extern void OffsetLocation (SeqLocPtr loc, Int4 offset, SeqIdPtr sip);

extern CharPtr kSubmitterUpdateText;
extern CharPtr kIndexerUpdateVecScreenText;
extern Boolean CreateUpdateCitSubFromBestTemplate (SeqEntryPtr top_sep, SeqEntryPtr upd_sep, CharPtr update_txt);
extern void AddCitSubToUpdatedSequence (BioseqPtr upd_bsp, Uint2 input_entityID, CharPtr update_txt);

extern Boolean AlistMessage (EnumFieldAssocPtr al, UIEnumPtr val, UIEnum dflt, CharPtr mssg);

typedef struct loginfo 
{
  FILE *fp;
  Boolean data_in_log;
  CharPtr display_title;
  Char path[PATH_MAX];  
} LogInfoData, PNTR LogInfoPtr;

extern LogInfoPtr OpenLog (CharPtr display_title);
extern void CloseLog (LogInfoPtr lip);
extern LogInfoPtr FreeLog (LogInfoPtr lip);

extern void LogCDSAmbiguousFrame (LogInfoPtr lip, SeqFeatPtr sfp);

extern void LoadGenomeProjectIDsFromFile (IteM i);
extern void RemoveEmptyGenomeProjectIDs (IteM i);

extern CharPtr SourceQualValNodeName (ValNodePtr vnp);
extern ValNodePtr SourceQualValNodeDataCopy (ValNodePtr vnp);
extern Boolean SourceQualValNodeMatch (ValNodePtr vnp1, ValNodePtr vnp2);

extern ValNodePtr GetSourceQualDescList (Boolean get_subsrc, Boolean get_orgmod, Boolean get_discouraged, Boolean get_discontinued);

extern void FeatureRemove (IteM i);
extern void ConvertFeatures (IteM i);
extern void SelectFeatures (IteM i);
extern void ReverseFeatureIntervals (IteM i);
extern void ParseDefLineToSourceQual (IteM i);
extern void ParseTaxnameToSourceQual (IteM i);
extern void ParseFlatfileToSourceQual (IteM i);
extern void ParseLocalIDToSourceQual (ButtoN b);
extern void FeatureEvidenceEditor (IteM i);
extern void FeatureExceptionEditor (IteM i);
extern void FeaturePartialEditor (IteM i);
extern void FeatureStrandEditor (IteM i);
extern void FeatureCitationEditor (IteM i);
extern void FeatureExperimentEditor (IteM i);
extern void FeatureInferenceEditor (IteM i);
extern void FeaturePseudoEditor (IteM i);
extern void ApplySourceQual (IteM i);
extern void PublicApplySourceQual (IteM i);
extern void EditSourceQual (IteM i);
extern void PublicEditSourceQual (IteM i);
extern void ConvertSourceQual (IteM i);
extern void SwapSourceQual (IteM i);
extern void RemoveSourceQual (IteM i);
extern void ApplyCDSGeneProt (IteM i);
extern void PublicApplyCDSGeneProt (IteM i);
extern void EditCDSGeneProt (IteM i);
extern void PublicEditCDSGeneProt (IteM i);
extern void ConvertCDSGeneProt (IteM i);
extern void SwapCDSGeneProt (IteM i);
extern void RemoveCDSGeneProt (IteM i);
extern void ApplyRNAQual (IteM i);
extern void PublicApplyRNAQual (IteM i);
extern void EditRNAQual (IteM i);
extern void PublicEditRNAQual (IteM i);
extern void ConvertRNAQual (IteM i);
extern void SwapRNAQual (IteM i);
extern void RemoveRNAQual (IteM i);
extern void ApplyGBQual (IteM i);
extern void PublicApplyGBQual (IteM i);
extern void EditGBQual (IteM i);
extern void PublicEditGBQual (IteM i);
extern void ConvertGBQual (IteM i);
extern void SwapGBQual (IteM i);
extern void RemoveGBQual (IteM i);
extern void ConvertLocusTagToOldLocusTag (IteM i);
extern void ExportLastLineage (IteM i);

extern void MacroApplyGBQual (IteM i);
extern void MacroApplySourceQual (IteM i);
extern void MacroApplyCDSGeneProt (IteM i);
extern void PublicMacroApplyCDSGeneProt (IteM i);
extern void MacroApplyRNAQual (IteM i);

extern void MacroRemoveGBQual (IteM i);
extern void MacroRemoveSourceQual (IteM i);
extern void MacroRemoveCDSGeneProt (IteM i);
extern void MacroRemoveRNAQual (IteM i);

extern void MacroConvertGBQual (IteM i);
extern void MacroConvertSourceQual (IteM i);
extern void MacroConvertCDSGeneProt (IteM i);
extern void MacroConvertRNAQual (IteM i);

extern void MacroSwapGBQual (IteM i);
extern void MacroSwapSourceQual (IteM i);
extern void MacroSwapCDSGeneProt (IteM i);
extern void MacroSwapRNAQual (IteM i);

extern void MacroEditGBQual (IteM i);
extern void MacroEditSourceQual (IteM i);
extern void MacroEditCDSGeneProt (IteM i);
extern void PublicMacroEditCDSGeneProt (IteM i);
extern void MacroEditRNAQual (IteM i);

extern void MacroApplyStructuredComment (IteM i);
extern void MacroEditStructuredComment (IteM i);
extern void MacroRemoveStructuredComment (IteM i);


/* constraint values */
#define LOCATION_CONSTRAINT_WHOLE_INTERVAL  1
#define LOCATION_CONSTRAINT_START_ENDPOINT  2
#define LOCATION_CONSTRAINT_STOP_ENDPOINT   3

#define LOCATION_CONSTRAINT_ANY        1
#define LOCATION_CONSTRAINT_UPSTREAM   2
#define LOCATION_CONSTRAINT_DOWNSTREAM 3
#define LOCATION_CONSTRAINT_CONTAINED  4
#define LOCATION_CONSTRAINT_NOT_IN     5
#define LOCATION_CONSTRAINT_OVERLAP    6
#define LOCATION_CONSTRAINT_EQUAL      7

#define LOCATION_CONSTRAINT_ANY_STRAND   1
#define LOCATION_CONSTRAINT_PLUS_STRAND  2
#define LOCATION_CONSTRAINT_MINUS_STRAND 3

#define LOCATION_CONSTRAINT_ANY_SEQ      1
#define LOCATION_CONSTRAINT_NUC_SEQ      2
#define LOCATION_CONSTRAINT_PROT_SEQ     3

typedef struct LocationConstraintX
{
  Int4      left;
  Int4      right;
  Int4      interval_end_choice;
  Int4      match_choice;
  Int4      strand;
  Int4      sequence_type;
} LocationConstraintXData, PNTR LocationConstraintXPtr;

typedef enum
{
  eStringConstraintContains = 1,
  eStringConstraintEquals,
  eStringConstraintStarts,
  eStringConstraintEnds,
  eStringConstraintInList 
} EStringConstraintMatchLocation;

typedef struct stringconstraint
{
  CharPtr match_text;
  Int4    match_location;
  Boolean insensitive;
  Boolean whole_word;
  Boolean not_present;
} StringConstraintData, PNTR StringConstraintXPtr;

extern StringConstraintXPtr StringConstraintXFree (StringConstraintXPtr scp);

typedef struct pseudoconstraint
{
  Boolean is_pseudo;
  Int4    featdef_type;
} PseudoConstraintData, PNTR PseudoConstraintPtr;

#define CHOICE_CONSTRAINT_ANY          1
#define CHOICE_CONSTRAINT_QUAL_PRESENT 3
#define CHOICE_CONSTRAINT_STRING       5
#define CHOICE_CONSTRAINT_MATCH        7
#define CHOICE_CONSTRAINT_PSEUDO       9

typedef struct choiceconstraint
{
  Int4                constraint_type;
  ValNodePtr          qual_choice;
  ValNodePtr          qual_choice_match;
  StringConstraintXPtr string_constraint;
  PseudoConstraintPtr pseudo_constraint;
  FreeValNodeProc     free_vn_proc;
  CopyValNodeDataProc copy_vn_proc;
} ChoiceConstraintData, PNTR ChoiceConstraintPtr;

extern ChoiceConstraintPtr ChoiceConstraintFree (ChoiceConstraintPtr scp);

typedef struct sequenceconstraint 
{
  Boolean nucs_ok;
  Boolean prots_ok;
  
  Int4                other_constraint_type;
  StringConstraintXPtr string_constraint;
  ChoiceConstraintPtr source_constraint;
  ValNodePtr          feature_list;
  
} SequenceConstraintXData, PNTR SequenceConstraintXPtr;

extern SequenceConstraintXPtr SequenceConstraintXFree (SequenceConstraintXPtr scp);
extern DialoG SequenceConstraintXDialog (GrouP g);
extern Boolean DoesSequenceMatchSequenceConstraintX (BioseqPtr bsp, SequenceConstraintXPtr scp);



typedef struct filterset 
{
  StringConstraintXPtr   scp;
  ChoiceConstraintPtr   ccp;
  LocationConstraintXPtr lcp;
  ChoiceConstraintPtr   cgp;
  StringConstraintXPtr   id_list;
} FilterSetData, PNTR FilterSetPtr;

extern void FilterSetClearText (FilterSetPtr fsp);
extern FilterSetPtr FilterSetNew (void);
extern FilterSetPtr FilterSetFree (FilterSetPtr fsp);

extern Boolean DoesStringMatchConstraintX (CharPtr pchSource, StringConstraintXPtr scp);

typedef CharPtr (*GetFeatureFieldString) PROTO ((SeqFeatPtr, ValNodePtr, FilterSetPtr));
typedef void (*SetFeatureFieldString) PROTO ((SeqFeatPtr, Pointer, FilterSetPtr));
typedef void (*RemoveFeatureFieldString) PROTO ((SeqFeatPtr, Pointer, FilterSetPtr));
typedef CharPtr (*GetDescriptorFieldString) PROTO ((SeqDescrPtr, ValNodePtr, FilterSetPtr));
typedef void (*SetDescriptorFieldString) PROTO ((SeqDescrPtr, Pointer, FilterSetPtr));
typedef void (*RemoveDescriptorFieldString) PROTO ((SeqDescrPtr, Pointer, FilterSetPtr));
typedef void (*FeatureActionProc) PROTO ((SeqFeatPtr, Pointer, FilterSetPtr));
typedef void (*DescriptorActionProc) PROTO ((SeqDescrPtr, Pointer, FilterSetPtr));
typedef Boolean (*OkToPreSample) PROTO ((Uint2 entityID));

extern void 
OperateOnSeqEntryConstrainedObjects 
(SeqEntryPtr           sep,
 FilterSetPtr          fsp,
 FeatureActionProc     feature_action,
 DescriptorActionProc  descriptor_action,
 Uint1                 seqFeatChoice,
 Uint1                 featDefChoice,
 Uint1                 descriptorChoice,
 Pointer               userdata);

extern CharPtr HandleApplyValue (CharPtr orig_text, ApplyValuePtr avp);
extern ValNodePtr 
ApplyValueToValNodeStringList 
(ValNodePtr list, Int2 choice, ApplyValuePtr avp);

typedef  Boolean (*Nlm_AcceptActnProc) PROTO((Pointer));
typedef  void  (*Nlm_CancelActnProc) PROTO ((Pointer));
typedef  void  (*Nlm_ClearActnProc) PROTO ((Pointer));
typedef  void  (*Nlm_ClearTextActnProc) PROTO ((Pointer));

extern ValNodePtr ValNodeFuncFree (ValNodePtr vnp, FreeValNodeProc free_vn_proc);

typedef struct textportion
{
  Int4    start_choice;
  CharPtr start_text;
  Int4    end_choice;
  CharPtr end_text;
  Boolean insensitive;
  Boolean whole_word;
} TextPortionXData, PNTR TextPortionXPtr;

extern TextPortionXPtr TextPortionXFree (TextPortionXPtr tp);
extern void 
FindTextPortionXInString 
(CharPtr        str, 
 TextPortionXPtr tp, 
 CharPtr PNTR   ploc, 
 Int4Ptr        plen);

extern DialoG TextPortionXDialogEx (GrouP h, Boolean inside, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
extern DialoG TextPortionXDialog (GrouP h);

#define CONVERT_TYPE_MOVE  0
#define CONVERT_TYPE_COPY  1
#define CONVERT_TYPE_SWAP  2
#define CONVERT_TYPE_PARSE 3

typedef struct convertfield
{
  ValNodePtr                  src_field_list;
  ValNodePtr                  dst_field_list;
  ExistingTextPtr             etp;
  Int2                        convert_type;
  GetFeatureFieldString       get_str_func;
  SetFeatureFieldString       set_str_func;
  RemoveFeatureFieldString    remove_str_func;
  GetDescriptorFieldString    get_d_str_func;
  SetDescriptorFieldString    set_d_str_func;
  RemoveDescriptorFieldString remove_d_str_func;
  NameFromValNodeProc         name_field_func;
  FilterSetPtr                fsp;
  TextPortionXPtr              text_portion;
  Boolean                     strip_name_from_text;
  Boolean                     remove_parsed;
} ConvertFieldData, PNTR ConvertFieldPtr;

extern DialoG StringConstraintDialogX (GrouP h, CharPtr label, Boolean clear_btn);
extern DialoG LocationConstraintXDialog (GrouP h, Boolean show_interval_controls, Boolean clear_btn);

enum pub_field_nums 
{
  PUB_FIELD_ANY = 0,
  PUB_FIELD_TITLE,
  PUB_FIELD_FIRST_NAME,
  PUB_FIELD_MIDDLE_INITIAL,
  PUB_FIELD_LAST_NAME,
  PUB_FIELD_SUFFIX,
  PUB_FIELD_CONSORTIUM,
  PUB_FIELD_INSTITUTION,
  PUB_FIELD_DEPARTMENT,
  PUB_FIELD_ADDRESS,
  PUB_FIELD_CITY,
  PUB_FIELD_STATE,
  PUB_FIELD_COUNTRY,
  PUB_FIELD_ZIP,
  PUB_FIELD_EMAIL,
  PUB_FIELD_PHONE,
  PUB_FIELD_FAX
};

enum pub_status
{
  PUB_STAT_ANY = 0,
  PUB_STAT_PUBLISHED,
  PUB_STAT_UNPUBLISHED,
  PUB_STAT_INPRESS,
  PUB_STAT_PUBLISHED_SUBMISSION
};

typedef struct pubconstraint
{
  CharPtr find_str;
  Int4    field_for_find;
  Boolean insensitive_to_case;
  Int4    pub_status;
} PubConstraintData, PNTR PubConstraintPtr;

extern PubConstraintPtr PubConstraintFree (PubConstraintPtr pcp);
extern DialoG PubConstraintDialog (GrouP h);

extern DialoG AcceptCancelDialog 
(GrouP                 parent,
 Nlm_AcceptActnProc    accept_actn,
 Nlm_CancelActnProc    cancel_actn,
 Nlm_ClearActnProc     clear_actn,
 Nlm_ClearTextActnProc clear_text_actn,
 Pointer               userdata,
 WindoW                w);
extern void EnableAcceptCancelDialogAccept (DialoG d);
extern void DisableAcceptCancelDialogAccept (DialoG d);

/* note - set sep to NULL if you don't want to limit the list to the features present */
extern ValNodePtr BuildFeatureDialogList (Boolean list_most_used_first, SeqEntryPtr sep);

extern DialoG 
FeatureSelectionDialog 
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG 
FeatureSelectionDialogEx 
(GrouP                    h,
 Boolean                  allow_multi,
 SeqEntryPtr              sep,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG
DescriptorSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG SourceQualTypeSelectionDialog 
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG 
FeatureFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Int4                     num_fields,
 CharPtr PNTR             field_names,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG 
GeneFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
extern CharPtr GetGeneFieldString (SeqFeatPtr sfp, ValNodePtr gene_field, FilterSetPtr fsp);
extern void RemoveGeneFieldString (SeqFeatPtr sfp, ValNodePtr gene_field);

extern DialoG 
MRNAFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
extern CharPtr GetmRNAFieldString (SeqFeatPtr sfp, ValNodePtr mrna_field, FilterSetPtr fsp);
extern void RemovemRNAFieldString (SeqFeatPtr sfp, ValNodePtr mrna_field);
extern CharPtr GetCDSFieldString (SeqFeatPtr sfp, ValNodePtr cds_field, FilterSetPtr fsp);
extern void RemoveCDSFieldString (SeqFeatPtr sfp, ValNodePtr cds_field);
extern CharPtr GetProteinFieldString (SeqFeatPtr sfp, ValNodePtr protein_field, FilterSetPtr fsp);

extern DialoG 
ProteinFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG 
CDSGeneProtFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);  
extern CharPtr GetCDSGeneProtField (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp);
extern Boolean 
SetCDSGeneProtField 
(SeqFeatPtr      sfp,
 ValNodePtr      vnp, 
 ApplyValuePtr   avp,
 FilterSetPtr    fsp);
extern void RemoveCDSGeneProtField (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp);
extern Uint2 FeatDefTypeFromFieldList (ValNodePtr vnp);
extern Boolean IsCDSetProteinProductChoice (ValNodePtr vnp);
extern CharPtr GetCDSGeneProtFieldName (ValNodePtr vnp);

extern DialoG
RNAAddFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
 
extern DialoG
RNARemoveFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG
RNAFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
 
extern CharPtr GetRNAFieldString (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp);

extern DialoG 
ExonFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
extern CharPtr GetExonFieldString (SeqFeatPtr sfp, ValNodePtr exon_field);
extern void RemoveExonFieldString (SeqFeatPtr sfp, ValNodePtr exon_field);

typedef  DialoG  (*FeatureFieldSelectionProc) PROTO((GrouP, Boolean, Nlm_ChangeNotifyProc, Pointer));

extern DialoG FeatureFieldChoiceDialog 
(GrouP h,
 FeatureFieldSelectionProc make_fieldlist_dlg,
 Boolean                   offer_to_remove,
 Nlm_ChangeNotifyProc      change_notify,
 Pointer                   change_userdata);

extern DialoG BioSourceStringDialog 
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG 
ConstraintChoiceDialog 
(GrouP                     h, 
 FeatureFieldSelectionProc present_func,
 FeatureFieldSelectionProc string_func,
 FreeValNodeProc           free_vn_proc,
 CopyValNodeDataProc       copy_vn_proc,
 CharPtr                   present_name,
 CharPtr                   text_name,
 Boolean                   clear_btn,
 Boolean                   use_pseudo);
extern DialoG SourceConstraintDialogX (GrouP h, Boolean clear_btn);
extern Boolean DoesOneSourceMatchConstraint (BioSourcePtr biop, ChoiceConstraintPtr scp);
extern DialoG CDSGeneProtConstraintDialog (GrouP h, Boolean clear_btn);
extern DialoG 
FilterGroup 
(GrouP h,
 Boolean has_string_constraint,
 Boolean has_source_constraint,
 Boolean has_location_constraint,
 Boolean has_cds_gene_prot_constraint, 
 Boolean has_id_list_constraint,
 CharPtr string_constraint_label);
 
typedef struct parsefield
{
  Int4        parse_field_type;
  ValNodePtr  feature_field;
  ValNodePtr  feature_subtype;
  Boolean     do_feat;
  Boolean     do_desc;
} ParseFieldData, PNTR ParseFieldPtr;

extern ParseFieldPtr ParseFieldFree (ParseFieldPtr pfp);
 
extern DialoG ParseFieldDestDialogEx 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  is_search_field,
 Boolean                  include_dbxref);
extern DialoG ParseFieldDestDialog 
(GrouP                    h, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
extern DialoG ParseFieldSourceDialog
(GrouP                    h,
 SeqEntryPtr              sep,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern DialoG SampleDialog (GrouP h);

extern NucProtAssocPtr FreeAssociationList (NucProtAssocPtr assoc_list);
extern NucProtAssocPtr 
AssignProteinsForSequenceSet 
(SeqEntryPtr nuc_list,
 SeqEntryPtr prot_list,
 Boolean always_review);

extern CharPtr 
CreateListMessage 
(CharPtr    msg_before,
 CharPtr    msg_after, 
 ValNodePtr id_list);

extern void StringToLower (CharPtr str);

extern Boolean ExportSubmitterBlockTemplate (SeqEntryPtr sep, SeqDescrPtr sdp);
extern DialoG OrganismSelectionDialog (GrouP parent, CharPtr org_name);

typedef struct getsample
{
  GetFeatureFieldString    fieldstring_func;
  GetDescriptorFieldString descrstring_func;
  
  ValNodePtr               feat_dest_list;
  ValNodePtr               descr_dest_list;
  
  ValNodePtr               requested_field;
  FreeValNodeProc          free_vn_proc;
  CopyValNodeDataProc      copy_vn_proc;

  CharPtr                  sample_text;
  Int4                     num_found;
  Boolean                  all_same;
} GetSampleData, PNTR GetSamplePtr;

extern GetSamplePtr GetSampleNew (void);
extern GetSamplePtr GetSampleFree (GetSamplePtr gsp);


extern GetSamplePtr
GetSampleForSeqEntry
(SeqEntryPtr   sep,
 Uint2         entityID, 
 ParseFieldPtr dst_field_data,
 FilterSetPtr  fsp);


extern void ApplyGDSKeyword (IteM i);
extern void ApplyTPAInferentialKeyword (IteM i);
extern void ApplyTPAExperimentalKeyword (IteM i);
extern void ApplyTPAReassemblyKeyword (IteM i);
extern void ApplyKeywordWithStringConstraint (IteM i);
extern void RemoveKeywordWithStringConstraint (IteM i);

#if defined(OS_UNIX) || defined(OS_MSWIN)
extern Int2 LIBCALLBACK CorrectRNAStrandedness (Pointer data);
extern Int2 LIBCALLBACK CorrectRNAStrandednessUseSmart (Pointer data);
#endif

extern void AddGeneFeatureFromTitle (SeqEntryPtr nucsep, CharPtr ttl,  SeqLocPtr slp);
extern SeqFeatPtr AddProteinFeatureFromDefline (SeqEntryPtr psep, CharPtr title);
extern void AddCodingRegionFieldsFromProteinTitle (CdRegionPtr crp, CharPtr title, CharPtr PNTR pcomment);
extern Boolean ReplaceImportModifierName (CharPtr PNTR orig_name, Int4 col_num);

extern SeqEntryPtr 
ImportSequencesFromFileEx
(FILE           *fp, 
 SeqEntryPtr     sep_list,
 Boolean         is_na, 
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped,
 Boolean         allow_char_stripping);

extern SeqEntryPtr 
ImportSequencesFromFile
(FILE           *fp, 
 SeqEntryPtr     sep_list,
 Boolean         is_na, 
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped);

extern void TestUpdateSequenceIndexer (IteM i);
extern void TestUpdateSequenceSubmitter (IteM i);
extern void TestExtendSequenceIndexer (IteM i);
extern void TestExtendSequenceSubmitter (IteM i);
extern void TestUpdateSequenceSetIndexer (IteM i);
extern void TestUpdateSequenceSetSubmitter (IteM i);
extern void TestExtendSequenceSetIndexer (IteM i);
extern void TestExtendSequenceSetSubmitter (IteM i);
extern void UpdateSequenceViaDownloadIndexer (IteM i);
extern void UpdateSequenceViaDownloadSubmitter (IteM i);

extern void 
ListBioseqsInSeqEntry 
(SeqEntryPtr     sep, 
 Boolean         is_na,
 Int4Ptr         seq_num, 
 ValNodePtr PNTR bioseq_list);

extern void RecomputeSuggestedIntervalsForCDS 
(Uint2          entityID,
 BioseqPtr PNTR batchbsp,
 Int4Ptr        count,
 MonitorPtr     mon,
 SeqFeatPtr     sfp);
 
 typedef struct recompdata {
  Int4        count;
  MonitorPtr  mon;
  BioseqPtr   batchbsp;
  Boolean     include_stop;
  Boolean     no_stop_at_end_of_complete_cds;
  Boolean     fix_genes;
  Uint2       entityID;
} RecompData, PNTR RecompDataPtr;

extern void RecomputeIntervalsForOneCDS (SeqFeatPtr sfp, RecompDataPtr rdp);

 
extern CharPtr ExtendProtein3 
(SeqFeatPtr sfp,
 Uint2      input_entityID,
 Boolean    force_partial);

extern SeqLocPtr 
ExpandSeqLoc 
(Int4 start,
 Int4 stop,
 Uint1 strand,
 BioseqPtr bsp,
 SeqLocPtr slp);

extern void SetSeqLocStrand (SeqLocPtr location, Uint1 strand);
extern void SetPrimerBindPairStrands (IteM i);
extern Boolean IsBioseqInAnyAlignment (BioseqPtr bsp, Uint2 input_entityID);
extern void RemoveAlignmentsWithSequence (BioseqPtr bsp, Uint2 input_entityID);
extern void FlipEntireAlignmentIfAllSequencesFlipped (SeqAnnotPtr sap, Pointer userdata);

extern void ConvertSelectedGapFeaturesToKnown (IteM i);
extern void ConvertSelectedGapFeaturesToUnknown (IteM i);
extern void CombineAdjacentGaps (IteM i);

extern void CombineAdjacentGapsOnBioseq (BioseqPtr bsp, Pointer userdata);

extern void MarkGenesWithPseudoFeaturesPseudo (IteM i);

typedef struct gaplocinfo 
{
  Int4    start_pos;
  Boolean is_known;
  Int4    length;
  Boolean replace;
} GapLocInfoData, PNTR GapLocInfoPtr;

extern void 
PrepareCodingRegionLocationsForDeltaConversionCallback
(BioseqPtr bsp, Pointer userdata);

extern void RemoveNomenclature (IteM i);

extern void ParseCollectionDateMonthFirst (IteM i);
extern void ParseCollectionDateDayFirst (IteM i);
extern void RemoveUnpublishedPublications (IteM i);
extern void RemovePublishedPublications (IteM i);
extern void RemoveUnindexedFeatures (IteM i);
extern void CopyLocusToLocusTag (IteM i);

typedef struct cdsconversionopts {
  Boolean all_are_pseudo;
  Boolean remove_mrna;
  Boolean remove_gene;
  Boolean remove_transcript_id;
  Boolean only_pseudo;
} CDSConversionOptsData, PNTR CDSConversionOptsPtr;

extern Boolean IsFeatInGPS (SeqFeatPtr sfp);

typedef struct cdstomiscfeat {
  Boolean              viral;
  Boolean              must_have_stops;
  CDSConversionOptsPtr opts;
} CDStoMiscFeatData, PNTR CDStoMiscFeatPtr;

extern void ConvertCDSToMiscFeat (SeqFeatPtr sfp, Pointer userdata);
extern void AdjustCodingRegionsEndingInGap (IteM i);

extern void ConvertCodingRegionsWithInternalKnownGapToMiscFeat (IteM i);

extern void FixOneAlignmentOverGaps (SeqAlignPtr salp, Uint2 entityID);
extern void ConsolidateSegmentsOverKnownLengthGaps (SeqAlignPtr salp);

extern void CreateReportWindow (EDiscrepancyReportType report_type);
extern void ScrollToDiscrepancyItem (ValNodePtr vnp, Pointer userdata);
extern void EditDiscrepancyItem (ValNodePtr vnp, Pointer userdata);
extern void WriteClickableListReport (FILE *fp, ValNodePtr discrepancy_list, Boolean show_all, Boolean use_feature_table_fmt);

extern void ConvertGapFeaturesToUnknown (IteM i);
extern void ChangeKnownGapLength (IteM i);
extern void AddFlankingNsToKnownLengthGaps (IteM i);

extern Int2 GetSequinAppParam (CharPtr section, CharPtr type, CharPtr dflt, CharPtr buf, Int2 buflen);
extern Boolean DoBioseqFeaturesMatchSequenceConstraintX (BioseqPtr bsp, ValNodePtr feat_list, StringConstraintXPtr scp);
extern Boolean DoesIDListMeetStringConstraint (SeqIdPtr sip, StringConstraintXPtr string_constraint);
extern Int2 LIBCALLBACK ReorderSetByAccession (Pointer data);

extern Int2 LIBCALLBACK CopyDescriptorToList (Pointer data);
 
extern void MapFeaturesToProteinSequence(IteM i);

extern void FindContig (IteM i);
extern ValNodePtr FreeSeqIdList (ValNodePtr id_list);
extern void DownloadAndDisplay (Int4 uid);
extern void LaunchDisplay (Uint2 entityID);
extern void SetBioseqViewTargetByBioseq (BaseFormPtr bfp, BioseqPtr bsp);

extern Boolean ShowDeltaReport (SeqEntryPtr sep);
extern void DeltaReport (IteM  i);

extern void ParseModifiersFromDefline (IteM i);
extern void ConvertLocalToGeneral(IteM i);

extern void FixLocusTagGeneXrefs (IteM i);
extern void FixLastExonLocNoPartial (IteM i);
extern void FixLastExonLocMakePartial (IteM i);

extern void RemoveGeneByUnderlyingFeatureType(IteM i);

extern void RemoveIntronLocationsFromCDS (IteM i);
extern void RemoveIntronLocationsFromrRNA (IteM i);
extern void RemoveIntronLocationsFromtRNA (IteM i);
extern void RemoveIntronLocationsFrommRNA (IteM i);

extern void ReverseBioseqInAlignment (SeqAlignPtr salp, Pointer userdata);

extern Int2 AddSeqAlignForSeqEntry (SeqEntryPtr sep, Uint2 entityID, Boolean choose_master, Boolean use_new_blast);

extern ValNodePtr ChooseFeaturesForConversion (ValNodePtr clickable_list, BaseFormPtr bfp, CharPtr label1, CharPtr label2);
extern void RemoveBadPubs (IteM i);
extern void VecScreenTool (IteM i);
extern void LogTrimmedLocation (LogInfoPtr lip, SeqLocPtr slp);
extern void CalculateVectorDescription (ClickableItemPtr cip);

extern void BarcodeTestTool (IteM i);

extern Boolean RelaxedSeqIdIn (SeqIdPtr sip, SeqIdPtr sip_list);

extern void NewSUC (ValNodePtr suc_list, Uint2 entityID, Boolean reverse, Boolean byblock, Boolean showsequence);
extern ValNodePtr GetSUCCommonList (SeqEntryPtr sep, Boolean reverse, Boolean byblock, Boolean showsequence, Boolean byqual);
extern ValNodePtr CategorizeSUCBlocks (ValNodePtr head);

extern void MakeGeneralIDsFromLocusTags (IteM i);

extern void ShowClickableItemList (ValNodePtr clickable_list, BaseFormPtr bfp, CharPtr win_title, CharPtr label1, CharPtr label2);

extern Int4 CountChosenDiscrepancies (ValNodePtr discrepancy_list, Boolean count_all);
extern void AddTranslExcept (SeqFeatPtr sfp, CharPtr cds_comment, Boolean use_strict, Boolean extend, Boolean adjust_gene);

enum table_data_errors 
{
  TABLE_DATA_NO_ERROR = 0,
  TABLE_DATA_ALREADY_HAS,
  TABLE_DATA_CELL_BLANK,
  TABLE_DATA_MULTIPLE_VALUES,
  TABLE_DATA_NOT_FOUND
};

extern Boolean 
GetTableOptions 
(BaseFormPtr bfp,
 ValNodePtr clickable_list,
 CharPtr win_title,
 CharPtr label1,
 CharPtr label2,
 CharPtr skip_already_txt,
 CharPtr blanks_erase_txt,
 BoolPtr skip_already_has,
 BoolPtr blanks_erase);

extern void MakeCDSmRNAPairs (IteM i);

extern void Apply16SRNA (IteM i);
extern void Apply23SRNA (IteM i);
extern void Apply18SRNA (IteM i);
extern void Apply28SRNA (IteM i);
extern void Apply26SRNA (IteM i);
extern void Apply12SRNA (IteM i);
extern void ApplySmallRNA (IteM i);
extern void ApplyLargeRNA (IteM i);

/* capitalization */
extern void GetOrgNamesInRecordCallback (BioSourcePtr biop, Pointer userdata);
extern void FixCapitalizationInTitle (CharPtr PNTR pTitle, Boolean first_is_upper, ValNodePtr org_names);
extern void FixCaseByField (IteM i);

extern void ApplyRNA_ITS (IteM i);

extern BioSourcePtr ExtractFromDeflineToBioSource (CharPtr defline, BioSourcePtr biop);
extern BioSourcePtr 
ExtractFromTitleToBioSourceOrgMod 
(CharPtr      title,
 BioSourcePtr biop, 
 CharPtr      mod_name,
 Int4         subtype);
extern BioSourcePtr 
ExtractFromTitleToBioSourceSubSource 
(CharPtr      title,
 BioSourcePtr biop, 
 CharPtr      mod_name,
 Int4         subtype);
extern BioSourcePtr 
ExtractFromTitleToBioSourceCommonName 
(CharPtr      title,
 BioSourcePtr biop);

extern void CDSmRNALinkTool (IteM i);
extern void RemoveFeatureLink (SeqFeatPtr sfp1, SeqFeatPtr sfp2);
extern void LinkTwoFeatures (SeqFeatPtr dst, SeqFeatPtr sfp);

extern void CombineToCreatePseudoGene (IteM i);

extern void BulkEditCDS (IteM i);
extern void BulkEditGene (IteM i);
extern void BulkEditRNA (IteM i);

extern void BulkEditorFeatList (Uint2 entityID, ValNodePtr feat_list);
extern void BulkEditorDescrList (Uint2 entityID, ValNodePtr descr_list);
extern void BulkEditDiscrepancy (ValNodePtr vnp, Pointer userdata);
extern void BulkEditorCheckAllDialog (DialoG dlg);

extern Uint1 GetSubtypeForBulkEdit (ValNodePtr feat_list);

extern CDSConversionOptsPtr 
GetCDSConversionOptions (Boolean all_are_pseudo, Boolean any_pseudo, Boolean any_gps, BoolPtr cancel);

extern Boolean 
ConvertOneCDSToMiscFeat 
(SeqFeatPtr sfp,
 Boolean viral,
 Boolean must_have_stops,
 CDSConversionOptsPtr opts);

extern void SuppressGenesOnFeaturesInsideMobileElements (IteM i);

extern BaseFormPtr GetBaseFormForEntityID (Uint2 entityID);

extern ValNodePtr ParseAccessionNumberListFromString (CharPtr list_str, SeqEntryPtr sep);

extern SubSourcePtr FindBadLatLon (BioSourcePtr biop);
extern DialoG LatLonTestResultsDisplay (GrouP h);
extern void LatLonTool (IteM i);
extern DialoG SpecificHostResultsDisplay (GrouP h);
extern void FixSpecificHostValues (IteM i);
extern DialoG LatLonCountryResultsDisplay (GrouP h);
extern Pointer GetLatLonCountryCorrection (Uint1 data_choice, Pointer data, Pointer metadata);
extern void LatLonCountryTool (IteM i);
extern DialoG CountryTestResultsDisplay (GrouP h, Pointer metadata);
extern void CountryFixupTool (IteM i);

extern DialoG TaxFixDisplay (GrouP h);
extern Boolean IsTaxNameBad (OrgRefPtr org);
extern void TaxFixTool (IteM i);

extern void SetTransgenicOnSourceDescWhenSourceFeatPresent (IteM i);
extern void SetFocusOnSourceDescWhenSourceFeatPresent (IteM i);

extern void MakeBadSpecificHostValueTable (IteM i);

extern const char *nucleotide_alphabet;
extern const char *protein_alphabet;

extern DialoG AlnSettingsDlg (GrouP h, Boolean allow_sequence_type);
extern TSequenceInfoPtr GetDefaultSequenceInfo (void);
extern TSequenceInfoPtr GetAlignmentOptions (Uint1Ptr moltype, TSequenceInfoPtr sequence_info);


extern void tRNAScanUpdate (IteM i);
extern void ListFailedTaxonomyLookups (IteM i);
extern SeqFeatPtr GetGeneForFeature (SeqFeatPtr sfp);

extern void LoadFeatureFieldTable (IteM i);

/* structure used for taxname options for loading values froma table or the Apply/Edit/Convert/Remove dialogs */
typedef struct taxnameoptions {
 Boolean remove_taxref;
 Boolean remove_old_name;
 Boolean remove_common;
} TaxnameOptionsData, PNTR TaxnameOptionsPtr;

extern void ApplyTaxnameOptionsToBioSource (BioSourcePtr biop, TaxnameOptionsPtr top);
extern DialoG TaxnameOptionsDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

extern void CleanTableLine (CharPtr line);

extern void ApplyTagToCodingRegionsInSourceFeatures (IteM i);
extern Boolean EntityIDAlreadyInList (Uint2 entityID, ValNodePtr entityIDList);

extern void ResolveFeatureOverlaps (IteM i);
extern void ConvertGeneralIdToLocalID (IteM i);

extern void UpdateFeatures (IteM i);

NLM_EXTERN void RemoveDupGenBankSets (SeqEntryPtr sep);
NLM_EXTERN BioseqSetPtr InsetNewSet (BioseqSetPtr bssp, Uint1 _class);

extern void AbbreviateCitSubAffilStates (IteM i);
extern void RemoveQualityScores (BioseqPtr bsp, FILE *log_fp, BoolPtr data_in_log);
extern void NewLoadFeatureQualifierTable (IteM i);
extern void NewLoadSourceQualifierTable (IteM i);
extern void ListAllSequences (BioseqPtr bsp, Pointer userdata);
extern void ChooseCategories (ValNodePtr value_list, Boolean do_choose);
extern void ChooseCategoriesByStringConstraint (ValNodePtr value_list, StringConstraintXPtr scp, Boolean do_choose);
extern void CapitalizeFirstLetterOfEveryWord (CharPtr pString);

typedef void (*BulkSetFieldFunc) PROTO ((Pointer, Pointer));
typedef Pointer (*BulkSetFieldStringFunc) PROTO ((Pointer, ApplyValuePtr));
typedef Pointer (*BulkGetFieldFunc) PROTO ((Uint1, Pointer, Pointer));
typedef CharPtr (*BulkDisplayFieldFunc) PROTO ((Pointer));
typedef void (*BulkFreeFieldFunc) PROTO ((Pointer));
typedef DialoG (*BulkCreateDlgFunc) PROTO ((GrouP, CharPtr, SeqEntryPtr));
typedef Int4 (*BulkFormatColumnFunc) PROTO ((ColPtr, CharPtr));
typedef void (*BulkDrawColumnFunc) PROTO ((Pointer, RectPtr));
typedef Pointer (*BulkReleaseCellFunc) PROTO ((Pointer));
typedef Pointer (*BulkCopyFieldFunc) PROTO ((Pointer));

typedef struct bulkedfield {
  CharPtr name;
  BulkSetFieldFunc set_func;
  BulkSetFieldStringFunc set_str_func;
  BulkGetFieldFunc get_func;
  BulkDisplayFieldFunc display_func;
  BulkFreeFieldFunc free_func;
  BulkCreateDlgFunc create_dlg_func;
  BulkFormatColumnFunc format_col_func;
  BulkDrawColumnFunc  draw_col_func;
  BulkReleaseCellFunc release_cell_func;
  BulkCopyFieldFunc   copy_func;
} BulkEdFieldData, PNTR BulkEdFieldPtr;

NLM_EXTERN Pointer BulkSetSimpleTextString (Pointer curr_val, ApplyValuePtr avp);
NLM_EXTERN CharPtr BulkDisplaySimpleText (Pointer data);
NLM_EXTERN void BulkFreeSimpleText (Pointer data);
NLM_EXTERN Int4 BulkFormatSimpleText (ColPtr col, CharPtr name);
NLM_EXTERN Pointer BulkSimpleTextCopy (Pointer data);
NLM_EXTERN DialoG BulkSimpleTextDialog (GrouP g, CharPtr name, SeqEntryPtr sep);

NLM_EXTERN void BulkEditorObjectList (Uint2 entityID, CharPtr title, ValNodePtr feat_list, BulkEdFieldPtr field_list);
NLM_EXTERN DialoG 
CreateBulkEditorDialog 
(GrouP               h,
 BulkEdFieldPtr      field_list,
 ValNodePtr          feat_list,
 SeqEntryPtr         sep, 
 Boolean             collapse_by_default,
 ClickableCallback   single_click_func,
 ClickableCallback   double_click_func);

NLM_EXTERN void ApplyBulkEditorToObjectList (DialoG d);

NLM_EXTERN void FlipSequenceIntervals (IteM i);

NLM_EXTERN void CreateRefSeqProteinIDs (IteM i);

extern void ConvertBioSourceDbxrefToFeatureDbxref (IteM i);

NLM_EXTERN void AddFluComments (IteM i);
NLM_EXTERN void CreateStructuredCommentsItem (IteM i);
NLM_EXTERN void SubmitterCreateStructuredComments (IteM i);

NLM_EXTERN SeqAlignPtr GetSeqAlignTSA (BioseqPtr bsp1, BioseqPtr bsp2);
NLM_EXTERN void EditTSAAssembly (IteM i);

NLM_EXTERN void MergeBiosources (IteM i);
NLM_EXTERN Boolean FixIDsAndTitles (SeqEntryPtr new_list, SeqEntryPtr current_list, Boolean is_nuc);

NLM_EXTERN void ExportQualifiers (IteM i);

extern Boolean IsUnwantedFeatureType (Uint1 key);
extern void ExportBankitComments (IteM i);

extern void ImportAlignmentForSeqHistInterval (IteM i);
extern void AdvancedAssemblyAlignmentEditor (IteM i);
extern void AssemblyAlignmentIntervalResolution (IteM i);

extern void ExternalSourceQualifierTableReader (IteM i);

extern void TrimPrimerSeqJunk (IteM i);

extern void SequinSeqViewFormMessage (ForM f, Int2 mssg);
extern Boolean WriteTheEntityID (Uint2 entityID, CharPtr path, Boolean binary);

extern BioseqSetPtr FindTopLevelSetForDesktopFunction (BioseqSetPtr bssp);
extern void ReorderSetByAccessionMenuItem (IteM i);
extern void DescriptorPropagateMenuItem (IteM i);
NLM_EXTERN void RepackagePartsMenuItem (IteM i);
extern void NewSegregateBioseqSetMenuItem (IteM i);
extern void SequesterSequencesMenuItem (IteM i);
extern void GetRidOfSegGapMenuItem (IteM i);
extern void GenerateSeqAlignMenuItem (IteM i);
extern void UpdateSeqAlignMenuItem (IteM i);
extern void CorrectRNAStrandednessMenuItem (IteM i);
extern void BioseqRevCompByIDMenuItem (IteM i);
extern void RemoveSetsInSetMenuItem (IteM i);

extern void LoadTaxTableReader (IteM i);

#ifdef OS_MSWIN
NLM_EXTERN Int4 RunSilent(const char *cmdline);
#endif


#ifdef __cplusplus
}
#endif

#endif /* ndef _SEQUIN_ */

