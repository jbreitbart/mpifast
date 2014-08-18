/*   sqnutils.h
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
* File Name:  sqnutils.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   9/2/97
*
* $Revision: 6.342 $
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

#ifndef _SQNUTILS_
#define _SQNUTILS_

#include <ncbi.h>
#include <sequtil.h>
#include <objpubme.h>
#include <objentgene.h>
#include <util/creaders/alnread.h>
#include <subutil.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef  void  (*Nlm_ChangeNotifyProc) PROTO ((Pointer));

NLM_EXTERN SeqEntryPtr LIBCALL GetTopSeqEntryForEntityID (Uint2 entityID);
NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForData (Uint2 entityID, BioseqPtr bsp);
NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemID (Uint2 entityID, Uint4 itemID, Uint2 itemtype);

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForDataEx (Uint2 entityID, BioseqPtr bsp, Boolean skipGenProdSet);
NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemIDEx (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Boolean skipGenProdSet);

NLM_EXTERN SeqIdPtr SeqIdFindWorst (SeqIdPtr sip);
NLM_EXTERN void ChangeSeqIdToWorstID (SeqIdPtr sip);
NLM_EXTERN void ChangeSeqLocToWorstID (SeqLocPtr slp);

NLM_EXTERN SeqIdPtr MakeSeqID (CharPtr str);
NLM_EXTERN SeqIdPtr MakeUniqueSeqID (CharPtr prefix);

NLM_EXTERN DatePtr DateAdvance (DatePtr dp, Uint1 monthsToAdd);

NLM_EXTERN SeqEntryPtr LIBCALL FindNthSeqEntry (SeqEntryPtr sep, Int2 seq);
NLM_EXTERN SeqEntryPtr LIBCALL FindNthBioseq (SeqEntryPtr sep, Int2 seq);
NLM_EXTERN SeqEntryPtr LIBCALL FindNthSequinEntry (SeqEntryPtr sep, Int2 seq);
NLM_EXTERN SeqEntryPtr LIBCALL FindNucSeqEntry (SeqEntryPtr sep);
NLM_EXTERN BioseqPtr LIBCALL FindNucBioseq (SeqEntryPtr sep);
NLM_EXTERN SeqEntryPtr LIBCALL FindBioseqSetByClass (SeqEntryPtr sep, Uint1 _class);

NLM_EXTERN Boolean LIBCALL SeqEntryHasNucs (SeqEntryPtr sep);
NLM_EXTERN Boolean LIBCALL SeqEntryHasProts (SeqEntryPtr sep);
NLM_EXTERN Boolean LIBCALL SeqEntryHasAligns (Uint2 entityID, SeqEntryPtr sep);
NLM_EXTERN Boolean LIBCALL PowerBLASTASN1Detected (SeqEntryPtr sep);

NLM_EXTERN Int2 EntityIDToGeneticCode (Uint2 entityID, BoolPtr mito, CharPtr taxname, size_t maxsize);
NLM_EXTERN Int2 SeqEntryToGeneticCode (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize);
NLM_EXTERN Int2 SeqEntryToBioSource (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize, BioSourcePtr PNTR biopp);

NLM_EXTERN Boolean BioseqToGeneticCode (
  BioseqPtr bsp,
  Int2Ptr gencodep,
  BoolPtr mitop,
  BoolPtr plastidp,
  CharPtr taxnamep,
  size_t maxsize,
  BioSourcePtr PNTR biopp
);

NLM_EXTERN SeqLocPtr   CreateWholeInterval (SeqEntryPtr sep);
NLM_EXTERN SeqFeatPtr  CreateNewFeature (SeqEntryPtr sep, SeqEntryPtr placeHere, Uint1 choice, SeqFeatPtr useThis);
NLM_EXTERN ValNodePtr  CreateNewDescriptor (SeqEntryPtr sep, Uint1 choice);

NLM_EXTERN Boolean IsPopPhyEtcSet (Uint1 _class);

/* Variants that call SeqMgrGetSeqEntryForData. The feature version allows a location
to be specified, overriding the default full-length seq-int location.  (If location is
not NULL, it copies it after deleting the existing sfp->location.)  For both functions
you still need to set the sfp->data.value.ptrvalue of the sdp->data.ptrvalue. */
NLM_EXTERN SeqFeatPtr CreateNewFeatureOnBioseq (BioseqPtr bsp, Uint1 choice, SeqLocPtr slp);
NLM_EXTERN ValNodePtr CreateNewDescriptorOnBioseq (BioseqPtr bsp, Uint1 choice);

NLM_EXTERN void        UpdateLocalId (BioseqPtr bsp, CharPtr localId);
NLM_EXTERN void        UpdateTitle (BioseqPtr bsp, CharPtr title);

NLM_EXTERN GeneRefPtr  CreateNewGeneRef (CharPtr locus, CharPtr allele,
                                     CharPtr desc, Boolean pseudo);
NLM_EXTERN ProtRefPtr  CreateNewProtRef (CharPtr name, CharPtr desc,
                                     CharPtr ec, CharPtr activity);
NLM_EXTERN CdRegionPtr CreateNewCdRgn (Uint1 frame, Boolean orf, Int2 genCode);

NLM_EXTERN void        SetSeqFeatData (SeqFeatPtr sfp, Pointer data);
NLM_EXTERN void        SetSeqFeatProduct (SeqFeatPtr sfp, BioseqPtr bsp);
NLM_EXTERN void        ResetSeqFeatInterval (SeqFeatPtr sfp);

NLM_EXTERN void        AddSeqFeatInterval (SeqFeatPtr sfp, BioseqPtr bsp, Int4 from, Int4 to,
                                       Boolean partial5, Boolean partial3);

NLM_EXTERN void        AddSeqLocPoint (SeqLocPtr PNTR old_slp, SeqIdPtr sip, Int4 location,
                                       Boolean fuzz_before, Boolean fuzz_after, Int2 strand);
NLM_EXTERN void        AddSeqFeatPoint (SeqFeatPtr sfp, BioseqPtr bsp, Int4 location, Boolean fuzz_before, Boolean fuzz_after, Int2 strand);

/* AddSeqEntryToSeqEntry and ReplaceSeqEntryWithSeqEntry leave
   the original target sep pointing to the new structure. */

NLM_EXTERN void        AddSeqEntryToSeqEntry (SeqEntryPtr target, SeqEntryPtr insert, Boolean relink);
NLM_EXTERN void        ReplaceSeqEntryWithSeqEntry (SeqEntryPtr target, SeqEntryPtr replaceWith, Boolean relink);

NLM_EXTERN void        RemoveSeqEntryFromSeqEntry (SeqEntryPtr top, SeqEntryPtr del, Boolean relink);
NLM_EXTERN void        RenormalizeNucProtSets (SeqEntryPtr sep, Boolean relink);

/* The following functions are called by the above when relink is TRUE.  Examine the
   code of ReplaceSeqEntryWithSeqEntry (in dlgutil2.c) to see how relink is treated. */

NLM_EXTERN void        GetSeqEntryParent (SeqEntryPtr target, Pointer PNTR parentptr, Uint2Ptr parenttype);

NLM_EXTERN void        SaveSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr PNTR omdptopptr, ObjMgrData PNTR omdataptr);
NLM_EXTERN void        RestoreSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr omdptop, ObjMgrData PNTR omdataptr);

/* If relink FALSE, call SeqMgrLinkSeqEntry (target, parenttype, parentptr)
   with original parent after all sequences have been added to the target. */

/* If relink FALSE, call SaveSeqEntryObjMgrData with the address of temporary
   ObjMgrDataPtr and ObjMgrData variables, and after calling SeqMgrLinkSeqEntry to
   update the link table, call RestoreSeqEntryObjMgrData with the value of the
   temporary ObjMgrDataPtr and the address of the ObjMgrData variable. */

/* ExtractBioSourceAndPubs and ReplaceBioSourceAndPubs can be called before and
   after AddSeqEntryToSeqEntry to propagate source and pub descriptors to top level. */

NLM_EXTERN ValNodePtr  ExtractBioSourceAndPubs (SeqEntryPtr sep);
NLM_EXTERN void        ReplaceBioSourceAndPubs (SeqEntryPtr sep, ValNodePtr descr);

/* SeqLocMerge combines feature intervals.  It can be used to extend the gene feature
   intervals, and (eventually) to fuse mutliple features into one. */

NLM_EXTERN SeqLocPtr SeqLocMerge (BioseqPtr target,
                                  SeqLocPtr to, SeqLocPtr from,
                                  Boolean single_interval, Boolean fuse_joints,
                                  Boolean add_null);

NLM_EXTERN SeqLocPtr SeqLocMergeEx (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                                    Boolean single_interval, Boolean fuse_joints,
                                    Boolean merge_overlaps, Boolean add_null);

NLM_EXTERN SeqLocPtr SeqLocMergeExEx (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                                      Boolean single_interval, Boolean fuse_joints,
                                      Boolean merge_overlaps, Boolean add_null,
                                      Boolean ignore_mixed);

NLM_EXTERN Boolean CheckSeqLocForPartial (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr);
NLM_EXTERN void SetSeqLocPartial (SeqLocPtr location, Boolean partial5, Boolean partial3);
NLM_EXTERN void FreeAllFuzz (SeqLocPtr location);
NLM_EXTERN Boolean LocationHasNullsBetween (SeqLocPtr location);
NLM_EXTERN Boolean SeqLocBadSortOrder (BioseqPtr bsp, SeqLocPtr slp);
NLM_EXTERN Boolean SeqLocMixedStrands (BioseqPtr bsp, SeqLocPtr slp);

/* GetBioseqGivenSeqLoc returns a segmented bioseq if the SeqLoc is to the parts */

NLM_EXTERN BioseqPtr GetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID);

NLM_EXTERN BioseqPtr GetBioseqGivenIDs (Uint2 entityID, Uint4 itemID, Uint2 itemtype);
NLM_EXTERN Uint4 GetItemIDGivenPointer (Uint2 entityID, Uint2 itemtype, Pointer lookfor);

NLM_EXTERN Uint1 FindFeatFromFeatDefType (Uint2 subtype);
NLM_EXTERN Uint1 FindFeatDefTypeFromKey (CharPtr key);
NLM_EXTERN CharPtr FindKeyFromFeatDefType (Uint1 type, Boolean forGBFF);

NLM_EXTERN Uint1 CodonToGcIndex (CharPtr codon);
NLM_EXTERN CharPtr GcIndextoCodon (Uint1 index);

/* finds bioseq from (cds) product, gets largest protein feature packaged on it */

NLM_EXTERN SeqFeatPtr LIBCALL GetBestProteinFeatureUnindexed (SeqLocPtr product);

/* resynchronizes coding regions with product protein bioseq and protein feature */

NLM_EXTERN void ResynchCodingRegionPartials (SeqEntryPtr sep);

/* resynchronizes mRNAs with product cDNA bioseq */

NLM_EXTERN void ResynchMessengerRNAPartials (SeqEntryPtr sep);

/* resynchronizes protein feature with product peptide bioseq */

NLM_EXTERN void ResynchProteinPartials (SeqEntryPtr sep);

/* individual feature callbacks for above functions */

NLM_EXTERN void ResynchMRNAPartials (SeqFeatPtr sfp, Pointer userdata);
NLM_EXTERN void ResynchCDSPartials (SeqFeatPtr sfp, Pointer userdata);
NLM_EXTERN void ResynchPeptidePartials (SeqFeatPtr sfp, Pointer userdata);

/* functions for associating CDS and parent mRNA using featureIDs */

NLM_EXTERN void ClearFeatIDs (SeqFeatPtr sfp);
NLM_EXTERN void ClearFeatIDXrefs (SeqFeatPtr sfp);

NLM_EXTERN void ClearFeatureIDs (SeqEntryPtr sep);
NLM_EXTERN Int4 FindHighestFeatureID (SeqEntryPtr sep);

NLM_EXTERN void AssignFeatureIDs (SeqEntryPtr sep);

NLM_EXTERN void OffsetFeatureIDs (SeqEntryPtr sep, Int4 offset);
NLM_EXTERN void OffsetFeatureIDXrefs (SeqEntryPtr sep, Int4 offset);

NLM_EXTERN void ReassignFeatureIDs (SeqEntryPtr sep);

NLM_EXTERN void LinkCDSmRNAbyOverlap (SeqEntryPtr sep);
NLM_EXTERN void LinkCDSmRNAbyProduct (SeqEntryPtr sep);
NLM_EXTERN void LinkCDSmRNAbyLabel (SeqEntryPtr sep);

NLM_EXTERN void StripFeatIDXrefAsnFilter (AsnIoPtr aip, AsnIoPtr aop);
NLM_EXTERN void StripSeqDataGapAsnFilter (AsnIoPtr aip, AsnIoPtr aop);
NLM_EXTERN void StripNewFeatMolInfoFieldsAsnFilter (AsnIoPtr aip, AsnIoPtr aop);
NLM_EXTERN void StripPCRPrimerAsnFilter (AsnIoPtr aip, AsnIoPtr aop);
    

/* functions to parse [org=Drosophila melanogaster] and [gene=lacZ] from titles */
/* for example, passing "gene" to SqnTagFind returns "lacZ" */

#define MAX_SQN_TAGS  32

typedef struct sqntag {
  CharPtr  query;
  Int2     num_tags;
  CharPtr  tag [MAX_SQN_TAGS];
  CharPtr  val [MAX_SQN_TAGS];
  Boolean  used [MAX_SQN_TAGS];
} SqnTag, PNTR SqnTagPtr;

NLM_EXTERN SqnTagPtr SqnTagParse (CharPtr ttl);
NLM_EXTERN SqnTagPtr SqnTagFree (SqnTagPtr stp);

NLM_EXTERN CharPtr SqnTagFind (SqnTagPtr stp, CharPtr tag);
NLM_EXTERN CharPtr SqnTagFindUnused (SqnTagPtr stp, CharPtr tag);

NLM_EXTERN void ReadTechFromString (CharPtr str, MolInfoPtr mip);
NLM_EXTERN void ReadCompletenessFromString (CharPtr str, MolInfoPtr mip);

extern Boolean StringsAreEquivalent (CharPtr str1, CharPtr str2);
NLM_EXTERN Uint1 EquivalentSubSource (CharPtr str);
NLM_EXTERN Uint1 EquivalentOrgMod (CharPtr str);
NLM_EXTERN Uint1 EquivalentSubSourceEx (CharPtr str, Boolean allow_discouraged_and_discontinued);
NLM_EXTERN Uint1 EquivalentOrgModEx (CharPtr str, Boolean allow_discouraged_and_discontinued);


/* functions to extract BioSource, MolInfo, and Bioseq information from parsed titles */

NLM_EXTERN BioSourcePtr ParseTitleIntoBioSource (
  SqnTagPtr stp,
  CharPtr organism,
  BioSourcePtr biop
);

NLM_EXTERN MolInfoPtr ParseTitleIntoMolInfo (
  SqnTagPtr stp,
  MolInfoPtr mip
);

NLM_EXTERN BioseqPtr ParseTitleIntoBioseq (
  SqnTagPtr stp,
  BioseqPtr bsp
);

NLM_EXTERN GeneRefPtr ParseTitleIntoGeneRef (
  SqnTagPtr stp,
  GeneRefPtr grp
);

NLM_EXTERN ProtRefPtr ParseTitleIntoProtRef (
  SqnTagPtr stp,
  ProtRefPtr prp
);

NLM_EXTERN GBBlockPtr ParseTitleIntoGenBank (
  SqnTagPtr stp,
  GBBlockPtr gbp
);

NLM_EXTERN SeqHistPtr ParseTitleIntoSeqHist (
  SqnTagPtr stp,
  SeqHistPtr shp
);

NLM_EXTERN SeqHistPtr ParseStringIntoSeqHist (
  SeqHistPtr shp,
  CharPtr str
);

NLM_EXTERN UserObjectPtr ParseTitleIntoTpaAssembly (
  SqnTagPtr stp,
  UserObjectPtr uop
);

NLM_EXTERN UserObjectPtr ParseTitleIntoGenomeProjectsDB (
  SqnTagPtr stp,
  UserObjectPtr uop
);

NLM_EXTERN Boolean IsGenomeProjectIDDescriptor (SeqDescrPtr sdp); 
NLM_EXTERN SeqDescrPtr GetGenomeProjectIDDescriptor (BioseqPtr bsp);


NLM_EXTERN void AddPubsFromTitle (
  SqnTagPtr stp,
  SeqDescrPtr PNTR desc_list
);

/* structured comment user object for flatfile presentation */

NLM_EXTERN UserObjectPtr ParseStringIntoStructuredComment (
  UserObjectPtr uop,
  CharPtr str,
  CharPtr prefix,
  CharPtr suffix
);


/* UseLocalAsnloadDataAndErrMsg transiently sets paths to asnload, data, and errmsg
  if they are packaged in the same directory as the executing program. */

NLM_EXTERN Boolean UseLocalAsnloadDataAndErrMsg (void);

/* GetRidOfLocusInSeqIds strips locus from all feature location and product seqIds */

NLM_EXTERN void GetRidOfLocusInSeqIds (Uint2 entityID, SeqEntryPtr sep);

NLM_EXTERN SeqLocPtr StripLocusFromSeqLoc (SeqLocPtr location);
NLM_EXTERN SeqIdPtr SeqIdStripLocus (SeqIdPtr sip);

/* LeaveBestCDD removes all but best CDD region in an area of overlapping features */

NLM_EXTERN void LeaveBestCDD (SeqEntryPtr sep);

/* ConvertPubSrcComDescsToFeats is useful when merging records */

NLM_EXTERN Boolean ConvertPubSrcComDescsToFeats (SeqEntryPtr sep, Boolean pub, Boolean src, Boolean com, Boolean toProts, Boolean PNTR asked_about_prop, Boolean PNTR propagate_descriptions, CharPtr findstring);

NLM_EXTERN void DeleteMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);

NLM_EXTERN Uint1 FindTrnaAA (CharPtr str);
NLM_EXTERN Uint1 FindTrnaAA3 (CharPtr str);
NLM_EXTERN Uint1 ParseTRnaString (CharPtr strx, BoolPtr justTrnaText, Uint1Ptr codon, Boolean noSingleLetter);
NLM_EXTERN CharPtr FindTrnaAAIndex (CharPtr str);
NLM_EXTERN Char FindResidueByName (CharPtr res_name, SeqCodeTablePtr sctp);
NLM_EXTERN ValNodePtr TokenizeTRnaString (CharPtr strx);
NLM_EXTERN Boolean ParseDegenerateCodon (tRNAPtr trp, Uint1Ptr codon);
NLM_EXTERN Boolean SerialNumberInString (CharPtr str);

/* ModernizeRNAFields uses new RNAGenPtr choice of RnaRef.ext for misc_RNA, ncRNA, tmRNA */

NLM_EXTERN void ModernizeRNAFields (
  SeqFeatPtr sfp
);

/* ModernizeGeneFields populates new GeneNomenclaturePtr field from OfficialNomenclature user object */

NLM_EXTERN void ModernizeGeneFields (
  SeqFeatPtr sfp
);

/* for sorting and uniquing valnode list by (charptr) data.ptrvalue */

NLM_EXTERN int LIBCALLBACK SortVnpByString (VoidPtr ptr1, VoidPtr ptr2);
NLM_EXTERN ValNodePtr UniqueValNode (ValNodePtr list);

/* for sorting and uniquing valnode list by data.intvalue */

NLM_EXTERN int LIBCALLBACK SortByIntvalue (VoidPtr ptr1, VoidPtr ptr2);
NLM_EXTERN ValNodePtr UniqueIntValNode (ValNodePtr list);

/* keytag sorts/uniques and then owns valnode character list */

typedef struct keytag {
  Int2               num;
  ValNodePtr         list;
  CharPtr PNTR       index; /* elements point into above valnode list */
} KeyTag;                   /* used as substructure, not allocated separately */

NLM_EXTERN void KeyTagInit (KeyTag PNTR ktp, ValNodePtr list);
NLM_EXTERN void KeyTagClear (KeyTag PNTR ktp);

NLM_EXTERN Int2 KeyFromTag (KeyTag PNTR ktp, CharPtr tag);
NLM_EXTERN CharPtr TagFromKey (KeyTag PNTR ktp, Int2 key);

/* inference qualifier utility */

#define VALID_INFERENCE            0
#define EMPTY_INFERENCE_STRING     1
#define BAD_INFERENCE_PREFIX       2
#define BAD_INFERENCE_BODY         3
#define SINGLE_INFERENCE_FIELD     4
#define SPACES_IN_INFERENCE        5
#define SAME_SPECIES_MISUSED       6
#define BAD_INFERENCE_ACCESSION    7
#define BAD_INFERENCE_ACC_VERSION  8
#define ACC_VERSION_NOT_PUBLIC     9
#define BAD_ACCESSION_TYPE        10

NLM_EXTERN Int2 ValidateInferenceQualifier (CharPtr val, Boolean fetchAccn);


/* from Colombe */
NLM_EXTERN SeqLocPtr StringSearchInBioseq (SeqIdPtr sip, CharPtr sub);

/*****************************************************************************
*
*   SequinEntryList (sep, mydata, mycallback, index, indent)
*       traverses all Seq-entry nodes beginning with sep
*       calls mycallback () at each node
*       Does enter BioseqSets of _class "parts", but ignores the
*       parts set itself
*
*****************************************************************************/

NLM_EXTERN Int4 SequinEntryList (SeqEntryPtr sep, Pointer mydata, SeqEntryFunc mycallback, Int4 index, Int2 indent);

#define SequinEntryCount( a )  SequinEntryList( a ,NULL,NULL,0,0)
#define SequinEntryExplore(a,b,c) SequinEntryList(a, b, c, 0L, 0)

/* Phrap reading function, based on sample code supplied by C. Magness, returns a SeqEntry list
of Bioseqs containing SeqGraphs, with individual reads removed and only contigs remaining */

NLM_EXTERN SeqEntryPtr ReadPhrapFile (FILE *fp);

/* Internal function to read quality scores, made available to parse separate DNA and quality score files */

NLM_EXTERN SeqGraphPtr ReadPhrapQuality (FILE *fp, BioseqPtr bsp);
NLM_EXTERN SeqGraphPtr ReadPhrapQualityFC (FileCachePtr fcp, BioseqPtr bsp);

/* SetPhrapContigOrder takes the results of ReadPhrapFile and a string indicating the order
of contigs, and returns a SeqEntryList in the desired order, with all other contigs removed */

NLM_EXTERN SeqEntryPtr SetPhrapContigOrder (SeqEntryPtr head, CharPtr contigs);

NLM_EXTERN void PrintQualityScores (BioseqPtr bsp, FILE *fp);

NLM_EXTERN void TrimSeqGraph (SeqGraphPtr sgp, Int4 num_to_trim, Boolean from_left);
NLM_EXTERN void TrimQualityScores (BioseqPtr bsp, Int4 num_to_trim, Boolean from_left);

typedef void (*QualityWriteFunc) (CharPtr buf, Uint4 buflen, Pointer userdata);

NLM_EXTERN void PrintQualityScoresToBuffer (BioseqPtr bsp, Boolean gapIsZero, Pointer userdata, QualityWriteFunc callback);

/* special function for genome contig delta sequences with far pointers */

NLM_EXTERN void PrintQualityScoresForContig (BioseqPtr bsp, Boolean gapIsZero, FILE* fp);

/* more efficient function for far genomic contig, makes separate graphs */

NLM_EXTERN SeqAnnotPtr PhrapGraphForContig (BioseqPtr bsp);

/* ReadContigList builds a far segmented bioseq from a table of accessions, starts, stops,
lengths, and (optional) strands.  Gaps of a given length (with 0 start and stop) are also
allowed. */

NLM_EXTERN SeqEntryPtr ReadContigList (FILE *fp, Boolean coordinatesOnMaster);
NLM_EXTERN SeqEntryPtr ReadContigListEx (FILE *fp, Boolean coordinatesOnMaster, CharPtr seqid, CharPtr title);

/* ReadAsnFastaOrFlatFile reads object manager-registered ASN.1, FASTA, GenBank, EMBL, GenPept,
Feature table, Restriction table, Contig table, Message response, or saved UID list, with the
option of saving FASTA results as OBJ_FASTA (SimpleSeq) to avoid ID collisions */

NLM_EXTERN Pointer ReadAsnFastaOrFlatFileEx (FILE *fp, Uint2Ptr datatypeptr, Uint2Ptr entityIDptr,
                                           Boolean forceNuc, Boolean forceProt,
                                           Boolean parseFastaSeqId, Boolean fastaAsSimpleSeq,
                                           BoolPtr chars_stripped);
NLM_EXTERN Pointer ReadAsnFastaOrFlatFile (FILE *fp, Uint2Ptr datatypeptr, Uint2Ptr entityIDptr,
                                           Boolean forceNuc, Boolean forceProt,
                                           Boolean parseFastaSeqId, Boolean fastaAsSimpleSeq);

/* ReadFeatureTableFile only handles >Feature tables */

NLM_EXTERN Pointer ReadFeatureTableFile (
  FILE *fp,
  Uint2Ptr datatypeptr,
  Uint2Ptr entityIDptr,
  Int4Ptr lineP,
  BoolPtr failP
);

/* ReadDeltaFasta reads a FASTA file, combining raw sequence and >?unk100 lines into
a delta Bioseq.  The file pointer stops at the next FASTA with a real SeqID. */

NLM_EXTERN BioseqPtr ReadDeltaFasta (FILE *fp, Uint2Ptr entityIDptr);

/* This function is identical to ReadDeltaFasta, except that the contents of 
 * chars_stripped will be set to TRUE if characters other than digits were stripped from
 * the sequence, or FALSE if not.
 */
NLM_EXTERN BioseqPtr ReadDeltaFastaEx (FILE *fp, Uint2Ptr entityIDptr, BoolPtr chars_stripped);

/* ReadDeltaFastaWithEmptyDefline reads just one delta sequence with an empty
 * definition line.
 * Calling function should make sure that fp is set to the start of the line 
 * with the empty definition line and that there is a "gap sequence ID"
 * present as the next definition line in the file.
 */
NLM_EXTERN BioseqPtr ReadDeltaFastaWithEmptyDefline (FILE *fp, Uint2Ptr entityIDptr, BoolPtr chars_stripped);

/* PromoteXrefs expands generef or protref feature cross-references (made by reading a
feature table with ReadAsnFastaOrFlatFile) to stand-alone gene features or protein features
and protein bioseqs.  It processes ALL features in the list - you give it the FIRST sfp. */

NLM_EXTERN void PromoteXrefs (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID);
NLM_EXTERN void PromoteXrefsEx (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID, Boolean include_stop,
                                Boolean remove_trailingX, Boolean gen_prod_set);
NLM_EXTERN void PromoteXrefsExEx (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID, Boolean include_stop,
                                  Boolean remove_trailingX, Boolean gen_prod_set, Boolean force_local_id);

/* SetEmptyGeneticCodes imposes genetic code on all coding regions within a feature table */

NLM_EXTERN void SetEmptyGeneticCodes (SeqAnnotPtr sap, Int2 genCode);

/* AddIntervalToLocation is a convenience function to add a single interval, and is called by
ReadAsnFastaOrFlatFile internally. */

NLM_EXTERN SeqLocPtr AddIntervalToLocation (SeqLocPtr loc, SeqIdPtr sip, Int4 start,
                                            Int4 stop, Boolean partial5, Boolean partial3);

/* AddQualifierToFeature applies cds product and gene qualifiers as protref or generef stored
as feature xrefs.  Most others (e.g., protein_id) are stored as gbquals.  PromoteXrefs can then
turn these special cases into the appropriate structures in fully expanded records. */

NLM_EXTERN void AddQualifierToFeature (SeqFeatPtr sfp, CharPtr qual, CharPtr val);

/* specialized string trimming functions */

NLM_EXTERN CharPtr TrimSpacesAndSemicolons (CharPtr str);
NLM_EXTERN CharPtr TrimSpacesAndJunkFromEnds (CharPtr str, Boolean allowEllipsis);

/* specialized cleanup for subsource and orgmod lists */
NLM_EXTERN void CleanSubSourceList (SubSourcePtr PNTR sspp, Uint1 location);
NLM_EXTERN void CleanOrgModList (OrgModPtr PNTR ompp);

/* used by original BankIt to merge multiple primer subsources */
NLM_EXTERN void CleanSubSourcePrimers (SubSourcePtr PNTR sspp);

NLM_EXTERN Boolean PubIsEffectivelyEmpty (PubdescPtr pdp);

/* extracts and reinserts descriptors in a standard order */
NLM_EXTERN void NormalizeDescriptorOrder (SeqEntryPtr sep);

/* BasicSeqEntryCleanup cleans up strings, moves gbquals to the appropriate field, and
does several other conversions, all without changing the itemID structure (which would
require reindexing) */

NLM_EXTERN void BasicSeqEntryCleanup (SeqEntryPtr sep);

/* CleanUpSeqFeat componenet of BasicSeqEntryCleanup, can be called for external features */

NLM_EXTERN void CleanUpSeqFeat (SeqFeatPtr sfp, Boolean isEmblOrDdbj, Boolean isJscan, Boolean stripSerial, ValNodePtr PNTR publist);

/* CautiousSeqEntryCleanup is a gradual consolidation and replacement of functions in SeriousSeqEntryCleanup,
which does change the itemID structure, and is intended to be safe for a retrofit of the ID database */

NLM_EXTERN void CautiousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge);

/* Convert a segmented or delta Bioseq to a raw Bioseq */

NLM_EXTERN void SegOrDeltaBioseqToRaw (BioseqPtr bsp);

NLM_EXTERN void ConvertSegSetsToDeltaSequences (SeqEntryPtr sep);

/* general purpose text finite state machine */
/* based on Practical Algorithms for Programmers by Binstock and Rex */

struct TextFsa;
typedef struct TextFsa* TextFsaPtr;

NLM_EXTERN TextFsaPtr TextFsaNew (void);
NLM_EXTERN void TextFsaAdd (TextFsaPtr tbl, CharPtr word);
NLM_EXTERN Int4 TextFsaNext (TextFsaPtr tbl, Int4 currState, Char ch, ValNodePtr PNTR matches);
NLM_EXTERN TextFsaPtr TextFsaFree (TextFsaPtr tbl);
NLM_EXTERN Boolean TextFsaGetStats (TextFsaPtr tbl, Int4Ptr highStateP, Int4Ptr numWordsP, Int4Ptr longestWordP);

/* PCR_primer manipulation functions */

typedef struct pcrset {
  CharPtr  fwd_seq;
  CharPtr  rev_seq;
  CharPtr  fwd_name;
  CharPtr  rev_name;
  Int2     orig_order;
} PcrSet, PNTR PcrSetPtr;

NLM_EXTERN ValNodePtr ParsePCRSet (BioSourcePtr biop);
NLM_EXTERN ValNodePtr ParsePCRStrings (
  CharPtr fwd_primer_seq,
  CharPtr rev_primer_seq,
  CharPtr fwd_primer_name,
  CharPtr rev_primer_name
);
NLM_EXTERN SubSourcePtr WritePCRSet (ValNodePtr pset);
NLM_EXTERN ValNodePtr FreePCRSet (ValNodePtr pset);

NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetSeq (VoidPtr ptr1, VoidPtr ptr2);
NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetOrder (VoidPtr ptr1, VoidPtr ptr2);

NLM_EXTERN ValNodePtr UniqueVnpByPCRSetSeq (ValNodePtr pset);

NLM_EXTERN void ModernizePCRPrimers (
  BioSourcePtr biop
);

/*
   very simple explore functions - VisitOn only does one chain, VisitIn goes into set components,
   they now return a count of the number of nodes visited, and the callback can be NULL if the purpose
   is simply to count nodes
*/

typedef void (*VisitDescriptorsFunc) (SeqDescrPtr sdp, Pointer userdata);
NLM_EXTERN Int4 VisitDescriptorsOnBsp (BioseqPtr bsp, Pointer userdata, VisitDescriptorsFunc callback);
NLM_EXTERN Int4 VisitDescriptorsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback);
NLM_EXTERN Int4 VisitDescriptorsInSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback);
NLM_EXTERN Int4 VisitDescriptorsOnSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback);
NLM_EXTERN Int4 VisitDescriptorsInSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback);

typedef void (*VisitFeaturesFunc) (SeqFeatPtr sfp, Pointer userdata);
NLM_EXTERN Int4 VisitFeaturesOnSap (SeqAnnotPtr sap, Pointer userdata, VisitFeaturesFunc callback);
NLM_EXTERN Int4 VisitFeaturesOnBsp (BioseqPtr bsp, Pointer userdata, VisitFeaturesFunc callback);
NLM_EXTERN Int4 VisitFeaturesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback);
NLM_EXTERN Int4 VisitFeaturesInSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback);
NLM_EXTERN Int4 VisitFeaturesOnSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback);
NLM_EXTERN Int4 VisitFeaturesInSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback);

typedef void (*VisitAlignmentsFunc) (SeqAlignPtr sap, Pointer userdata);
NLM_EXTERN Int4 VisitAlignmentsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitAlignmentsFunc callback);
NLM_EXTERN Int4 VisitAlignmentsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAlignmentsFunc callback);
NLM_EXTERN Int4 VisitAlignmentsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback);
NLM_EXTERN Int4 VisitAlignmentsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback);
NLM_EXTERN Int4 VisitAlignmentsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback);
NLM_EXTERN Int4 VisitAlignmentsInSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback);

typedef void (*VisitGraphsFunc) (SeqGraphPtr sgp, Pointer userdata);
NLM_EXTERN Int4 VisitGraphsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitGraphsFunc callback);
NLM_EXTERN Int4 VisitGraphsOnBsp (BioseqPtr bsp, Pointer userdata, VisitGraphsFunc callback);
NLM_EXTERN Int4 VisitGraphsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback);
NLM_EXTERN Int4 VisitGraphsInSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback);
NLM_EXTERN Int4 VisitGraphsOnSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback);
NLM_EXTERN Int4 VisitGraphsInSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback);

typedef void (*VisitAnnotsFunc) (SeqAnnotPtr sap, Pointer userdata);
NLM_EXTERN Int4 VisitAnnotsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAnnotsFunc callback);
NLM_EXTERN Int4 VisitAnnotsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback);
NLM_EXTERN Int4 VisitAnnotsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback);
NLM_EXTERN Int4 VisitAnnotsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback);
NLM_EXTERN Int4 VisitAnnotsInSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback);

typedef void (*VisitBioseqsFunc) (BioseqPtr bsp, Pointer userdata);
NLM_EXTERN Int4 VisitBioseqsInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioseqsFunc callback);
NLM_EXTERN Int4 VisitBioseqsInSep (SeqEntryPtr sep, Pointer userdata, VisitBioseqsFunc callback);

/* VisitSequences allows you to limit visitation to nucs or prots that aren't parts, or just to parts */

#define VISIT_MAINS 1
#define VISIT_NUCS  2
#define VISIT_PROTS 3
#define VISIT_PARTS 4

typedef void (*VisitSequencesFunc) (BioseqPtr bsp, Pointer userdata);
NLM_EXTERN Int4 VisitSequencesInSet (BioseqSetPtr bssp, Pointer userdata, Int2 filter, VisitSequencesFunc callback);
NLM_EXTERN Int4 VisitSequencesInSep (SeqEntryPtr sep, Pointer userdata, Int2 filter, VisitSequencesFunc callback);

typedef void (*VisitSetsFunc) (BioseqSetPtr bssp, Pointer userdata);
NLM_EXTERN Int4 VisitSetsInSep (SeqEntryPtr sep, Pointer userdata, VisitSetsFunc callback);
NLM_EXTERN Int4 VisitSetsInSet (BioseqSetPtr bssp, Pointer userdata, VisitSetsFunc callback);

/* visits components of pop/phy/mut/genbank sets, callback is at most nuc-prot set, can then call above functions */

typedef void (*VisitElementsFunc) (SeqEntryPtr sep, Pointer userdata);
NLM_EXTERN Int4 VisitElementsInSep (SeqEntryPtr sep, Pointer userdata, VisitElementsFunc callback);

/* visits all SeqIds within a SeqLoc, or within features, alignments, graphs, or annots */

typedef void (*VisitSeqIdFunc) (SeqIdPtr sip, Pointer userdata);
NLM_EXTERN Int4 VisitSeqIdsInSeqLoc (SeqLocPtr slp, Pointer userdata, VisitSeqIdFunc callback);

NLM_EXTERN Int4 VisitSeqIdsInBioseq (BioseqPtr bsp, Pointer userdata, VisitSeqIdFunc callback);
NLM_EXTERN Int4 VisitSeqIdsInSeqFeat (SeqFeatPtr sfp, Pointer userdata, VisitSeqIdFunc callback);
NLM_EXTERN Int4 VisitSeqIdsInSeqAlign (SeqAlignPtr sap, Pointer userdata, VisitSeqIdFunc callback);
NLM_EXTERN Int4 VisitSeqIdsInSeqGraph (SeqGraphPtr sgp, Pointer userdata, VisitSeqIdFunc callback);
NLM_EXTERN Int4 VisitSeqIdsInSeqAnnot (SeqAnnotPtr annot, Pointer userdata, VisitSeqIdFunc callback);

/* visits all sub UserFields - if the data type is 11, VisitUserFieldsInUfp recurses */

typedef void (*VisitUserFieldsFunc) (UserFieldPtr ufp, Pointer userdata);
NLM_EXTERN Int4 VisitUserFieldsInUfp (UserFieldPtr ufp, Pointer userdata, VisitUserFieldsFunc callback);
NLM_EXTERN Int4 VisitUserFieldsInUop (UserObjectPtr uop, Pointer userdata, VisitUserFieldsFunc callback);

/* visits all sub UserObjects if the data type is 12 - needed to pack multiple user objects on a single feature */

typedef void (*VisitUserObjectFunc) (UserObjectPtr uop, Pointer userdata);
NLM_EXTERN Int4 VisitUserObjectsInUop (UserObjectPtr uop, Pointer userdata, VisitUserObjectFunc callback);

/* explores sub UserObjects including "CombinedFeatureUserObjects" and finds by label  */

NLM_EXTERN UserObjectPtr FindUopByTag (UserObjectPtr top, CharPtr tag);

/* creates "CombinedFeatureUserObjects" sfp->ext to combine two user objects */

NLM_EXTERN UserObjectPtr CombineUserObjects (UserObjectPtr origuop, UserObjectPtr newuop);

/* visits all publication descriptors or features */

typedef void (*VisitPubdescsFunc) (PubdescPtr pdp, Pointer userdata);
NLM_EXTERN Int4 VisitPubdescsOnBsp (BioseqPtr bsp, Pointer userdata, VisitPubdescsFunc callback);
NLM_EXTERN Int4 VisitPubdescsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback);
NLM_EXTERN Int4 VisitPubdescsInSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback);
NLM_EXTERN Int4 VisitPubdescsOnSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback);
NLM_EXTERN Int4 VisitPubdescsInSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback);

/* visits all authors in a publication */

typedef void (*VisitAuthorFunc) (NameStdPtr nsp, Pointer userdata);
NLM_EXTERN Int4 VisitAuthorsInPub (PubdescPtr pdp, Pointer userdata, VisitAuthorFunc callback);

/* visits all biosource descriptors or features */

typedef void (*VisitBioSourcesFunc) (BioSourcePtr biop, Pointer userdata);
NLM_EXTERN Int4 VisitBioSourcesOnBsp (BioseqPtr bsp, Pointer userdata, VisitBioSourcesFunc callback);
NLM_EXTERN Int4 VisitBioSourcesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback);
NLM_EXTERN Int4 VisitBioSourcesInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback);
NLM_EXTERN Int4 VisitBioSourcesOnSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback);
NLM_EXTERN Int4 VisitBioSourcesInSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback);

/* function to scan binary ASN.1 file of entire release as Bioseq-set, simple explore from successive top seps */
/* compressed can be TRUE only on UNIX, where it does a popen on zcat to decompress on-the-fly */
/* although it now returns a count of components visited, the callback cannot be NULL for this function */

typedef void (*ScanBioseqSetFunc) (SeqEntryPtr sep, Pointer userdata);
NLM_EXTERN Int4 ScanBioseqSetRelease (CharPtr inputFile, Boolean binary, Boolean compressed, Pointer userdata, ScanBioseqSetFunc callback);

/* function to scan binary ASN.1 file of entrezgene release as Entrezgene-Set */

typedef void (*ScanEntrezgeneSetFunc) (EntrezgenePtr egp, Pointer userdata);
NLM_EXTERN Int4 ScanEntrezgeneSetRelease (CharPtr inputFile, Boolean binary, Boolean compressed, Pointer userdata, ScanEntrezgeneSetFunc callback);

/* PubMed registered fetch functionality */

NLM_EXTERN PubmedEntryPtr LIBCALL GetPubMedForUid (Int4 uid);

/* internal support type, registration function */

typedef PubmedEntryPtr (LIBCALLBACK * PubMedFetchFunc) (Int4 uid);

NLM_EXTERN void LIBCALL PubMedSetFetchFunc (PubMedFetchFunc func);

NLM_EXTERN void FirstNameToInitials (CharPtr first, CharPtr inits, size_t maxsize);

extern CharPtr MyFGetLine (FILE *fp, ValNodePtr PNTR current_data);

#if defined (WIN32)
extern char * __stdcall AbstractReadFunction (Pointer userdata);
extern void __stdcall AbstractReportError (TErrorInfoPtr err_ptr, Pointer userdata);
#else
extern char * AbstractReadFunction (Pointer userdata);
extern void AbstractReportError (TErrorInfoPtr err_ptr, Pointer userdata);
#endif

typedef struct readbuffer {
  FILE *fp;
  ValNodePtr current_data;
} ReadBufferData, PNTR ReadBufferPtr;

extern void FreeBufferedReadList (ValNodePtr vnp);

extern CharPtr AlignmentStringToSequenceString (CharPtr aln_str, Uint1 moltype);
extern SeqEntryPtr MakeSequinDataFromAlignment (TAlignmentFilePtr afp, Uint1 moltype);
extern SeqEntryPtr MakeSequinDataFromAlignmentEx (TAlignmentFilePtr afp, Uint1 moltype, Boolean check_ids);
extern SeqEntryPtr make_seqentry_for_seqentry (SeqEntryPtr sep);
extern void ProcessPseudoMiscFeatsForEntityID (Uint2 entityID);
extern Boolean ConvertOnePseudoCDSToMiscFeat (SeqFeatPtr sfp);
NLM_EXTERN Boolean ConvertOnePseudoCDSToMiscFeatEx (SeqFeatPtr sfp, Boolean remove_product);
extern void ConvertPseudoCDSToMiscFeatsForEntityID (Uint2 entityID);

extern SeqAlignPtr FindAlignmentsForBioseq (BioseqPtr bsp);
extern ValNodePtr FindAlignSeqAnnotsForBioseq (BioseqPtr bsp);
extern Boolean IsSequenceFirstInPairwise (SeqEntryPtr sep, SeqIdPtr sip);
extern Boolean RemoveSequenceFromAlignments (SeqEntryPtr sep, SeqIdPtr sip);
extern BioseqPtr ReadFastaOnly (FILE *fp,
                              Boolean forceNuc, Boolean forceProt,
                              BoolPtr chars_stripped,
                              CharPtr lastchar);
extern void MergeFeatureIntervalsToParts (SeqFeatPtr sfp, Boolean ordered);

extern void ExtendSingleGeneOnMRNA (BioseqPtr bsp, Pointer userdata);

/* structures and functions for the Discrepancy Report */
typedef void (*ClickableCallback) (ValNodePtr item_list, Pointer userdata);
typedef void (*ClickableCallbackDataFree) (Pointer userdata);

typedef struct clickableitem 
{
  Uint4                     clickable_item_type;
  CharPtr                   description;
  ValNodePtr                item_list;
  ClickableCallback         callback_func; 
  ClickableCallbackDataFree datafree_func; 
  Pointer                   callback_data;
  Boolean                   chosen;
  ValNodePtr                subcategories;
  Boolean                   expanded;
  Int4                      level;
} ClickableItemData, PNTR ClickableItemPtr;

extern ClickableItemPtr 
NewClickableItem 
(Uint4           clickable_item_type,
 CharPtr         description_fmt,
 ValNodePtr      item_list);
 
extern ClickableItemPtr ClickableItemFree (ClickableItemPtr cip);
extern ValNodePtr FreeClickableList (ValNodePtr list);

extern int LIBCALLBACK SortVnpByClickableItemDescription (VoidPtr ptr1, VoidPtr ptr2);
NLM_EXTERN int LIBCALLBACK SortVnpByClickableItemChosen (VoidPtr ptr1, VoidPtr ptr2);

extern void ExpandClickableItemList (ValNodePtr vnp);
extern void ContractClickableItemList (ValNodePtr vnp);

NLM_EXTERN void RemoveDuplicateItems (ValNodePtr PNTR item_list);

/* To add a new type of test, do ALL Of the following:
 * 1. add an item to the DiscrepancyType enum (this will fill the clickable_item_type value)
 * 2. add a collection function and declare it with the others
 * 3. add an item to discrepancy_info_list that corresponds with the position of the
 *    new enum value.  If you are combining multiple types in one collection function,
 *    be sure to list them together.
 */

typedef enum {
  DISC_GENE_MISSING = 0,
  DISC_SUPERFLUOUS_GENE,
  DISC_GENE_MISSING_LOCUS_TAG,
  DISC_GENE_DUPLICATE_LOCUS_TAG,
  DISC_GENE_LOCUS_TAG_BAD_FORMAT,
  DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX,
  DISC_NON_GENE_LOCUS_TAG,
  DISC_COUNT_NUCLEOTIDES,
  DISC_MISSING_PROTEIN_ID,
  DISC_INCONSISTENT_PROTEIN_ID_PREFIX,
  DISC_GENE_CDS_mRNA_LOCATION_CONFLICT,
  DISC_GENE_PRODUCT_CONFLICT,
  DISC_GENE_DUPLICATE_LOCUS,
  DISC_EC_NUMBER_NOTE,
  DISC_PSEUDO_MISMATCH,
  DISC_JOINED_FEATURES,
  DISC_OVERLAPPING_GENES,
  DISC_OVERLAPPING_CDS,
  DISC_CONTAINED_CDS,
  DISC_RNA_CDS_OVERLAP,
  DISC_SHORT_CONTIG,
  DISC_INCONSISTENT_BIOSRC,
  DISC_SUSPECT_PRODUCT_NAME,
  DISC_INCONSISTENT_BIOSRC_DEFLINE,
  DISC_PARTIAL_CDS_IN_COMPLETE_SEQUENCE,
  DISC_EC_NUMBER_ON_HYPOTHETICAL_PROTEIN,
  DISC_NO_TAXLOOKUP,
  DISC_BAD_TAXLOOKUP,
  DISC_SHORT_SEQUENCE,
  DISC_SUSPECT_PHRASES,
  DISC_SUSPICIOUS_NOTE_TEXT,
  DISC_COUNT_TRNA,
  DISC_DUP_TRNA,
  DISC_BADLEN_TRNA,
  DISC_STRAND_TRNA,
  DISC_COUNT_RRNA,
  DISC_DUP_RRNA,
  DISC_RNA_NO_PRODUCT,
  DISC_TRANSL_NO_NOTE,
  DISC_NOTE_NO_TRANSL,
  DISC_TRANSL_TOO_LONG,
  DISC_CDS_OVERLAP_TRNA,
  DISC_COUNT_PROTEINS,
  DISC_FEAT_OVERLAP_SRCFEAT,
  DISC_MISSING_GENPRODSET_PROTEIN,
  DISC_DUP_GENPRODSET_PROTEIN,
  DISC_MISSING_GENPRODSET_TRANSCRIPT_ID,
  DISC_DUP_GENPRODSET_TRANSCRIPT_ID,
  DISC_PERCENTN,
  DISC_N_RUNS,
  DISC_ZERO_BASECOUNT,
  DISC_ADJACENT_PSEUDOGENE,
  DISC_NO_ANNOTATION,
  DISC_INFLUENZA_DATE_MISMATCH,
  DISC_SHORT_INTRON,
  DISC_MISSING_VIRAL_QUALS,
  DISC_SRC_QUAL_PROBLEM,
  DISC_MISSING_SRC_QUAL,
  DISC_DUP_SRC_QUAL,
  DISC_DUP_SRC_QUAL_DATA,
  DISC_HAPLOTYPE_MISMATCH,
  DISC_FEATURE_MOLTYPE_MISMATCH,
  DISC_CDS_WITHOUT_MRNA,
  DISC_EXON_INTRON_CONFLICT,
  DISC_FEATURE_COUNT,
  DISC_SPECVOUCHER_TAXNAME_MISMATCH,
  DISC_GENE_PARTIAL_CONFLICT,
  DISC_FLATFILE_FIND_ONCALLER,
  DISC_CDS_PRODUCT_FIND,
  DISC_DUP_DEFLINE,
  DUP_DISC_ATCC_CULTURE_CONFLICT,
  DISC_USA_STATE,
  DISC_INCONSISTENT_MOLTYPES,
  DISC_SUBMITBLOCK_CONFLICT,
  DISC_POSSIBLE_LINKER,
  DISC_TITLE_AUTHOR_CONFLICT,
  DISC_BAD_GENE_STRAND,
  DISC_MAP_CHROMOSOME_CONFLICT,
  DISC_RBS_WITHOUT_GENE,
  DISC_CITSUBAFFIL_CONFLICT,
  DISC_REQUIRED_CLONE,
  DISC_PLASMODIUM_ISOLATE_NOT_STRAIN,
  DISC_SOURCE_QUALS_ASNDISC,
  DISC_mRNA_ON_WRONG_SEQUENCE_TYPE,
  DISC_RETROVIRIDAE_DNA,
  DISC_CHECK_AUTH_CAPS,
  DISC_CHECK_RNA_PRODUCTS_AND_COMMENTS,
  DISC_MICROSATELLITE_REPEAT_TYPE,
  DISC_MITOCHONDRION_REQUIRED,
  DISC_UNPUB_PUB_WITHOUT_TITLE,
  DISC_QUALITY_SCORES,
  DISC_INTERNAL_TRANSCRIBED_SPACER_RRNA,
  DISC_BACTERIAL_PARTIAL_PROBLEMS,
  MAX_DISC_TYPE
} DiscrepancyType;


typedef enum {
  eReportTypeDiscrepancy = 1,
  eReportTypeOnCaller
} EDiscrepancyReportType;

extern Boolean IsTestTypeAppropriateForReportType (Int4 test_type, EDiscrepancyReportType report_type);

extern void PrintDiscrepancyTestList (FILE *fp);


extern CharPtr GetDiscrepancyTestConfName (DiscrepancyType dtype);
extern CharPtr GetDiscrepancyTestSettingName (DiscrepancyType dtype);
extern DiscrepancyType GetDiscrepancyTypeFromSettingName (CharPtr setting_name);

typedef struct discrepancyconfig
{
  Boolean conf_list[MAX_DISC_TYPE];
  Boolean use_feature_table_format;
} DiscrepancyConfigData, PNTR DiscrepancyConfigPtr;

extern DiscrepancyConfigPtr DiscrepancyConfigFree (DiscrepancyConfigPtr dcp);
extern DiscrepancyConfigPtr DiscrepancyConfigNew (void);
extern DiscrepancyConfigPtr DiscrepancyConfigCopy (DiscrepancyConfigPtr dcp);
extern DiscrepancyConfigPtr ReadDiscrepancyConfig (void);
extern DiscrepancyConfigPtr ReadDiscrepancyConfigEx (CharPtr report_config_name);
extern void SaveDiscrepancyConfig (DiscrepancyConfigPtr dcp);
extern void DisableTRNATests (DiscrepancyConfigPtr dcp);
extern CharPtr SetDiscrepancyReportTestsFromString (CharPtr list, Boolean enable, DiscrepancyConfigPtr dcp);
extern void ConfigureForBigSequence (DiscrepancyConfigPtr dcp);
extern void ConfigureForGenomes (DiscrepancyConfigPtr dcp);

typedef void (*PerformDiscrepancyTest) PROTO ((ValNodePtr PNTR, ValNodePtr));

extern ValNodePtr CollectDiscrepancies (DiscrepancyConfigPtr dcp, ValNodePtr sep_list, PerformDiscrepancyTest taxlookup);
extern CharPtr GetDiscrepancyItemText (ValNodePtr vnp);
extern CharPtr GetDiscrepancyItemTextEx (ValNodePtr vnp, CharPtr filename);
extern void VisitGenProdSetFeatures (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback);
extern ValNodePtr ReplaceDiscrepancyItemWithFeatureTableStrings (ValNodePtr feat_list);
extern CharPtr GetParentLabelForDiscrepancyItem (ValNodePtr vnp);
extern void WriteDiscrepancy (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt);
extern void WriteDiscrepancyEx (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt, Boolean cmd_line, CharPtr descr_prefix, Boolean list_features_if_subcat);
extern int LIBCALLBACK SortVnpByDiscrepancyDescription (VoidPtr ptr1, VoidPtr ptr2);
extern int LIBCALLBACK SortVnpByDiscrepancyItemText (VoidPtr ptr1, VoidPtr ptr2);
extern void ValNodeReverse (ValNodePtr PNTR list);

/* Individual discrepancy test function declarations */
extern const CharPtr kOverlappingCDSNoteText;
extern const CharPtr kOverlappingCDSNeedsNoteFmt;
extern void AddOverlappingCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTagsEx (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Boolean exclude_dirsub);
extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void CheckBioSourceQuals (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);

/* structure shared by tbl2asn and discrepancy report functions */
typedef struct genprodsetdiscrepancylists {
  ValNodePtr cds_product_list;
  ValNodePtr mrna_product_list;
  ValNodePtr missing_mrna_product;
  ValNodePtr missing_protein_id;
} GenProdSetDiscrepancyListsData, PNTR GenProdSetDiscrepancyListsPtr;

extern void CheckGenProdSetsInSeqEntry (SeqEntryPtr sep, GenProdSetDiscrepancyListsPtr lists);  


typedef struct protidlists {
  ValNodePtr missing_gnl_list;
  ValNodePtr gnl_list;
} ProtIdListsData, PNTR ProtIdListsPtr;

/* structure shared by tbl2asn and discrepancy report functions */
typedef struct globaldiscrepancy {
  CharPtr str;
  Uint1   data_choice;
  Pointer data;
} GlobalDiscrepancyData, PNTR GlobalDiscrepancyPtr;

extern GlobalDiscrepancyPtr GlobalDiscrepancyNew (CharPtr str, Uint1 data_choice, Pointer data);
extern GlobalDiscrepancyPtr GlobalDiscrepancyFree (GlobalDiscrepancyPtr g);
extern ValNodePtr FreeGlobalDiscrepancyList (ValNodePtr vnp);
extern void ConvertGlobalDiscrepancyToText (GlobalDiscrepancyPtr g, Boolean use_feature_fmt, CharPtr filename);
extern void ConvertGlobalDiscrepancyListToText (ValNodePtr vnp, Boolean use_feature_fmt, CharPtr filename);
extern ValNodePtr GetGlobalDiscrepancyItem (GlobalDiscrepancyPtr g);
extern CharPtr GetGlobalDiscrepancyStr (GlobalDiscrepancyPtr g);
NLM_EXTERN int LIBCALLBACK SortVnpByGlobalDiscrepancyString (VoidPtr ptr1, VoidPtr ptr2);
extern ClickableItemPtr
ReportNonUniqueGlobalDiscrepancy 
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 CharPtr    ind_cat_fmt,
 Uint4      clickable_item_type,
 Boolean    keep_top_category);
extern ValNodePtr ReportInconsistentGlobalDiscrepancyPrefixes
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type);
extern ValNodePtr ReportInconsistentGlobalDiscrepancyStrings
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type);
extern ClickableItemPtr ReportMissingFields (ValNodePtr list, CharPtr label_fmt, Uint4 clickable_item_type);
extern ClickableItemPtr ReportBadLocusTagFormat (ValNodePtr list);
extern ClickableItemPtr FindAdjacentDuplicateLocusTagGenes (ValNodePtr locus_tag_list);
extern void FindProteinIDCallback (BioseqPtr bsp, Pointer userdata);


/* formats for global discrepancies also used by tbl2asn */
extern CharPtr discReportDuplicateLocusTagFmt;
extern CharPtr discReportOneDuplicateLocusTagFmt;
extern CharPtr discReportDuplicateProteinIDFmt;
extern CharPtr discReportOneDuplicateProteinIDFmt;
extern CharPtr discReportDuplicateTranscriptIdFmt;
extern CharPtr discReportOneDuplicateTranscriptIdFmt;
extern CharPtr discReportInconsistentLocusTagPrefixFmt;
extern CharPtr discReportMissingLocusTags;
extern CharPtr discReportInconsistentProteinIDPrefixFmt;
extern CharPtr discReportBadProteinIdFmt;
extern CharPtr discReportMissingTranscriptIDFmt;

extern CharPtr GetBioseqLabel (BioseqPtr bsp);
extern CharPtr GetBioseqSetLabel (BioseqSetPtr bssp);

NLM_EXTERN ValNodePtr ValNodeDupStringList (ValNodePtr vnp);


typedef enum {
  eLocusTagErrorBadFormat,
  eLocusTagErrorDuplicate,
  eLocusTagErrorInconsistentPrefix
} ELocusTagError;

NLM_EXTERN ValNodePtr FindBadLocusTagsInList (ValNodePtr list);


typedef struct discreportoutputconfig {
  Boolean use_feature_table_format;
  Boolean expand_report_categories[MAX_DISC_TYPE];
  Boolean summary_report;
} DiscReportOutputConfigData, PNTR DiscReportOutputConfigPtr;


typedef struct globaldiscrepreport {
  ValNodePtr  locus_tag_list;
  ValNodePtr  missing_locus_tag;
  ValNodePtr  cds_product_list;
  ValNodePtr  missing_cds_product;
  ValNodePtr  mrna_product_list;
  ValNodePtr  missing_mrna_product;
  ValNodePtr  adjacent_locus_tag_disc_list;
  ValNodePtr  missing_gnl_list;
  ValNodePtr  gnl_list;
  ValNodePtr  global_src_qual_vals;
  ValNodePtr  global_srcs;

  ValNodePtr  src_qual_repeated_list;
  ValNodePtr  src_qual_multi_list;
  ValNodePtr  feature_count_list;
  ValNodePtr  discrepancy_list;

  PerformDiscrepancyTest    taxlookup;
  DiscrepancyConfigPtr      test_config;
  DiscReportOutputConfigPtr output_config;
} GlobalDiscrepReportData, PNTR GlobalDiscrepReportPtr;

NLM_EXTERN GlobalDiscrepReportPtr GlobalDiscrepReportNew ();
NLM_EXTERN GlobalDiscrepReportPtr GlobalDiscrepReportFree (GlobalDiscrepReportPtr g);
NLM_EXTERN void AddSeqEntryToGlobalDiscrepReport (SeqEntryPtr sep, GlobalDiscrepReportPtr g, CharPtr filename);
NLM_EXTERN void WriteGlobalDiscrepancyReport (GlobalDiscrepReportPtr g, FILE *fp);


/* for the Barcode Discrepancy Test */
typedef enum {
  eBarcodeTest_Length = 0,
  eBarcodeTest_Primers,
  eBarcodeTest_Country,
  eBarcodeTest_SpecimenVoucher,
  eBarcodeTest_PercentN,
  eBarcodeTest_LAST
} EBarcodeTest;

typedef struct barcodetestconfig
{
  Boolean conf_list[eBarcodeTest_LAST];
  Int4    min_length;
  FloatLo min_n_percent;
  Boolean require_keyword;
} BarcodeTestConfigData, PNTR BarcodeTestConfigPtr;

extern BarcodeTestConfigPtr BarcodeTestConfigNew();
extern BarcodeTestConfigPtr BarcodeTestConfigFree (BarcodeTestConfigPtr cfg);

extern CharPtr GetBarcodeTestName (Int4 i);

extern Int4 GetBarcodeTestNumFromBarcodeTestName (CharPtr test_name);

typedef struct barcodetestresults
{
  Boolean failed_tests[eBarcodeTest_LAST];
  BioseqPtr bsp;
  FloatLo   n_percent;
} BarcodeTestResultsData, PNTR BarcodeTestResultsPtr;

extern BarcodeTestResultsPtr BarcodeTestResultsNew ();
extern BarcodeTestResultsPtr BarcodeTestResultsFree (BarcodeTestResultsPtr res);
extern BarcodeTestResultsPtr BarcodeTestResultsCopy (BarcodeTestResultsPtr res);
extern ValNodePtr            BarcodeTestResultsListFree (ValNodePtr res_list);

extern Boolean IsBarcodeID (SeqIdPtr sip);

extern CharPtr BarcodeTestBarcodeIdString (BioseqPtr bsp);
extern CharPtr BarcodeTestGenbankIdString (BioseqPtr bsp);

/* This one gets discrepancies by category */
extern ValNodePtr GetBarcodeDiscrepancies (ValNodePtr sep_list, BarcodeTestConfigPtr cfg);
/* This one lists accessions that fail */
extern ValNodePtr GetBarcodeFailedAccessionList (SeqEntryPtr sep, BarcodeTestConfigPtr cfg);
extern ValNodePtr GetBarcodePassFail (SeqEntryPtr sep, BarcodeTestConfigPtr cfg);
NLM_EXTERN CharPtr GetBarcodeTestFailureReasons (BarcodeTestResultsPtr res);
/* This one lists passes and failures, with reasons for failures */
extern void WriteBarcodeTestComprehensive (FILE *fp, ValNodePtr results_list);
extern void WriteBarcodeDiscrepancies (FILE *fp, ValNodePtr results_list);
extern void WriteBarcodeTestCompliance (FILE *fp, ValNodePtr results_list);
extern void WriteBarcodeTagTable (FILE *fp, ValNodePtr results_list);
NLM_EXTERN Boolean 
BarcodeValidateOneSeqEntry 
(FILE *ofp,
 SeqEntryPtr sep,
 Boolean show_all,
 Boolean use_xml,
 Boolean show_header,
 CharPtr xml_header_text);
extern void RemoveBarcodeKeywords (FILE *fp, ValNodePtr results_list);
extern void ApplyBarcodeKeywords (FILE *fp, ValNodePtr results_list);
extern Boolean PassBarcodeTests (BarcodeTestResultsPtr res);
extern Boolean HasBARCODEKeyword (BioseqPtr bsp);
extern Int4 CountPolymorphismsInBioseq (BioseqPtr bsp);


extern CharPtr ExpandDiscrepancyReportTestsFromString (CharPtr list, Boolean expand, DiscReportOutputConfigPtr dcp);
extern void CollateDiscrepancyReports (ValNodePtr PNTR discrepancy_reports);
extern void WriteAsnDiscReport (ValNodePtr discrepancy_list, FILE *ofp, DiscReportOutputConfigPtr oc, Boolean use_flag);


/* extern to allow access to subsource_subtype_alist */
typedef struct Nlm_qual_name_assoc {
   Nlm_CharPtr name; 
   Uint1       value;
} Nlm_QualNameAssoc, PNTR Nlm_QualNameAssocPtr, Nlm_QualNameAlist[];

typedef struct Nlm_name_name_assoc {
   Nlm_CharPtr name; 
   Nlm_CharPtr alias; 
   Uint1       value;
} Nlm_NameNameAssoc, PNTR Nlm_NameNameAssocPtr, Nlm_NameNameAlist[];

extern Nlm_QualNameAssoc current_orgmod_subtype_alist[];
extern Nlm_QualNameAssoc discouraged_orgmod_subtype_alist[];
extern Nlm_QualNameAssoc discontinued_orgmod_subtype_alist[];
extern Nlm_NameNameAssoc orgmod_aliases[];
extern CharPtr GetOrgModQualName (Uint1 subtype);
extern void BioSourceHasOldOrgModQualifiers (BioSourcePtr biop, BoolPtr has_discouraged, BoolPtr has_discontinued);

extern Nlm_QualNameAssoc  current_subsource_subtype_alist [];
extern Nlm_QualNameAssoc  discouraged_subsource_subtype_alist[];
extern Nlm_QualNameAssoc  discontinued_subsource_subtype_alist[];
extern Nlm_NameNameAssoc  subsource_aliases [];
extern CharPtr GetSubsourceQualName (Uint1 subtype);
extern void BioSourceHasOldSubSourceQualifiers (BioSourcePtr biop, BoolPtr has_discouraged, BoolPtr has_discontinued);
extern Boolean GeneRefMatch (GeneRefPtr grp1, GeneRefPtr grp2);
extern Boolean DbxrefsMatch (ValNodePtr vnp1, ValNodePtr vnp2, Boolean case_sensitive);
extern Boolean XrefsMatch (SeqFeatXrefPtr x1, SeqFeatXrefPtr x2);
extern Boolean ProtRefMatch (ProtRefPtr prp1, ProtRefPtr prp2);
NLM_EXTERN Boolean DoFeaturesMatch (SeqFeatPtr sfp1, SeqFeatPtr sfp2, Boolean allow_different_sequences, Boolean case_sensitive, Boolean ignore_partial);

extern void IsCorrectLatLonFormat (CharPtr lat_lon, BoolPtr format_correct, BoolPtr lat_in_range, BoolPtr lon_in_range);
extern CharPtr FixLatLonFormat (CharPtr orig_lat_lon);
extern Boolean ParseLatLon (CharPtr lat_lon, FloatHi PNTR latP, FloatHi PNTR lonP);
extern void ApplyBarcodeDbxrefsToBioseq (BioseqPtr bsp, Pointer data);

extern CharPtr GetCountryFix (CharPtr country, CharPtr PNTR country_list);

extern CharPtr ncrnaClassList[];
extern Int4 NcrnaOTHER;
extern Boolean IsStringInNcRNAClassList (CharPtr str);
extern ValNodePtr ListFeaturesInLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice);
extern ValNodePtr ListCodingRegionsContainedInSourceFeatures (SeqEntryPtr sep);
extern ValNodePtr ListFeaturesOverlappingLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice);

extern void ConvertSourceFeatDescProc (SeqFeatPtr sfp, Pointer userdata);

/* for correcting capitalization */
NLM_EXTERN void 
FixCapitalizationInElement 
(CharPtr PNTR pEl,
 Boolean      bAbbrev, 
 Boolean      bShortWords,
 Boolean      bApostrophes);

NLM_EXTERN void FixAbbreviationsInElement (CharPtr PNTR pEl);
NLM_EXTERN void FixOrgNamesInString (CharPtr str, ValNodePtr org_names);
NLM_EXTERN void ResetCapitalization (Boolean first_is_upper, CharPtr pString);

NLM_EXTERN SeqIdPtr CreateSeqIdFromText (CharPtr id_str, SeqEntryPtr sep);
NLM_EXTERN SeqLocPtr SeqLocWholeNew (BioseqPtr bsp);
NLM_EXTERN Int4 GetDeltaSeqLen (DeltaSeqPtr dsp);

typedef SeqAlignPtr (*GlobalAlignFunc) PROTO ((BioseqPtr, BioseqPtr, BoolPtr));

typedef struct adjustfeatforgap {
  ValNodePtr feature_list;
  Boolean    unknown_gaps;
  Boolean    known_gaps;
  Boolean    make_partial;
  Boolean    partial_for_pseudo;
  Boolean    trim_ends;
  Boolean    split_internal;
  GlobalAlignFunc align_func;
  ValNodePtr features_in_gap;
} AdjustFeatForGapData, PNTR AdjustFeatForGapPtr;

NLM_EXTERN AdjustFeatForGapPtr AdjustFeatForGapFree (AdjustFeatForGapPtr agp);
NLM_EXTERN Boolean FeatureOkForFeatureList (SeqFeatPtr sfp, ValNodePtr feature_list);
NLM_EXTERN void
LocationContainsGaps
(SeqLocPtr slp,
 BioseqPtr bsp,
 Boolean   unknown_gaps,
 Boolean   known_gaps,
 BoolPtr   terminal_gaps,
 BoolPtr   internal_gaps,
 BoolPtr   entirely_in_gap);

NLM_EXTERN void AdjustFeatureForGapsCallback (SeqFeatPtr sfp, Pointer data);
NLM_EXTERN void MarkFeaturesInGapsForDeletion (AdjustFeatForGapPtr afgp);
NLM_EXTERN void AdjustCDSLocationsForUnknownGapsCallback (SeqFeatPtr sfp, Pointer data);
NLM_EXTERN Boolean GapInLocation (Int4 seq_offset, Int4 length, SeqLocPtr loc);
NLM_EXTERN BioseqPtr 
AddProteinSequenceCopy 
(BioseqPtr  protbsp, 
 BioseqPtr  featbsp,
 SeqFeatPtr new_sfp,
 Uint2      entityID);
NLM_EXTERN void AdjustFrame (SeqFeatPtr sfp, BioseqPtr oldprot);
NLM_EXTERN void SetProductSequencePartials (BioseqPtr protbsp, Boolean partial5, Boolean partial3);
NLM_EXTERN void AddCDSGapComment (SeqFeatPtr sfp);


NLM_EXTERN Boolean SeqEdFixProteinFeatures (BioseqPtr oldbsp, BioseqPtr newbsp, Boolean force_fix, GlobalAlignFunc align_func);
NLM_EXTERN void SeqEdTranslateOneCDS (SeqFeatPtr sfp, BioseqPtr featbsp, Uint2 entityID, GlobalAlignFunc align_func);
NLM_EXTERN void SeqEdRemapLocation (SeqAlignPtr salp, SeqLocPtr slp, Int4 seq_len);

NLM_EXTERN CharPtr GetStateAbbreviation (CharPtr state);

typedef SeqAlignPtr (*LocalAlignFunc) PROTO ((BioseqPtr, BioseqPtr));

typedef struct transcriptomeids {
  BioseqPtr consensus_bsp;
  ValNodePtr token_list;
} TranscriptomeIdsData, PNTR TranscriptomeIdsPtr;

NLM_EXTERN TranscriptomeIdsPtr TranscriptomeIdsNew (BioseqPtr bsp, ValNodePtr token_list);
NLM_EXTERN TranscriptomeIdsPtr TranscriptomeIdsFree (TranscriptomeIdsPtr t);

NLM_EXTERN ValNodePtr TranscriptomeIdsListFree (ValNodePtr list);
NLM_EXTERN ValNodePtr GetTranscriptomeIdsList (FILE *fp, SeqEntryPtr sep, ValNodePtr PNTR err_list);
NLM_EXTERN ValNodePtr GetExistingTSATableIds (SeqEntryPtr sep);

NLM_EXTERN ValNodePtr 
ApplyTranscriptomeIdsListToSeqEntrySeqHist 
(ValNodePtr           list,
 LocalAlignFunc      aln_func,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata);
NLM_EXTERN ValNodePtr 
MakeTranscriptomeAssemblySeqHist 
(TranscriptomeIdsPtr t,
 LocalAlignFunc aln_func,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);

extern void ReverseAlignmentStrand (SeqAlignPtr salp, Int4 nth);

NLM_EXTERN Boolean HasExistingSeqHistAssembly (ValNodePtr list);
NLM_EXTERN void DeleteSeqHistAssembliesForList (ValNodePtr list);
NLM_EXTERN ValNodePtr ReportCoverageForTranscriptomeIdsListSeqHist (ValNodePtr list);

NLM_EXTERN SeqAlignPtr SortPairwiseAlignmentsByFirstSeqRange (SeqAlignPtr salp);
NLM_EXTERN ValNodePtr ReportCoverageForBioseqSeqHist (BioseqPtr bsp);
NLM_EXTERN ValNodePtr ReportConsensusMatchForBioseqSeqHist (BioseqPtr bsp);

NLM_EXTERN void ConvertLocalIdsToTSAIds (SeqEntryPtr sep);

NLM_EXTERN ValNodePtr MakeTokensFromLine (CharPtr line);

NLM_EXTERN SeqFeatPtr GetGeneForFeature (SeqFeatPtr sfp);

NLM_EXTERN Boolean IsStringInSpanInList (CharPtr str, CharPtr list);

NLM_EXTERN void ParseGoTermsFromFields (SeqEntryPtr sep);

/* for autodef */
typedef enum {
  RemovableExon = 0,
  RemovableIntron,
  Removable5UTR,
  Removable3UTR,
  RemovableCDS,
  RemovablePromoter,
  RemovableLTR,
  RemovableNoncodingProductFeat,
  RemovableMobileElement,
  RemovablePrecursorRNA,
  NumRemovableItems
} RemovableList;
NLM_EXTERN CharPtr GetRemovableItemName (Int4 i);

typedef enum {
  DEFLINE_USE_FEATURES = 1,
  DEFLINE_COMPLETE_SEQUENCE,
  DEFLINE_COMPLETE_GENOME
} DefLineType;

typedef struct deflinefeaturerequestlist {
  Boolean      keep_items[NumRemovableItems];
  Boolean      add_fake_promoters;
  Boolean      suppress_alt_splice_phrase;
  Boolean      remove_subfeatures;
  DefLineType  feature_list_type;
  Int4         misc_feat_parse_rule;
  Boolean      suppress_locus_tags;
  ValNodePtr   suppressed_feature_list;
  Boolean      use_ncrna_note;
} DeflineFeatureRequestList, PNTR DeflineFeatureRequestListPtr;

NLM_EXTERN void InitFeatureRequests (DeflineFeatureRequestListPtr feature_requests);


/* ModifierItemLocalData is used to store information about the results of
 * a search of the set of organisms in a record and the results of user
 * input to a dialog for deciding which modifiers should be used in the
 * organism description.
 */
typedef struct modifieritemlocal {
/*  ButtoN        button; */
  Boolean       any_present;
  Boolean       all_present;
  Boolean       is_unique;
  CharPtr       first_value_seen;
  ValNodePtr    values_seen;
  Boolean       all_unique;
  CharPtr       status;
  Boolean       required;
} ModifierItemLocalData, PNTR ModifierItemLocalPtr;

typedef enum {
  DEFLINE_POS_Acronym = 0,
  DEFLINE_POS_Anamorph,
  DEFLINE_POS_Authority,
  DEFLINE_POS_Bio_material,
  DEFLINE_POS_Biotype,
  DEFLINE_POS_Biovar,
  DEFLINE_POS_Breed,
  DEFLINE_POS_Cell_line,
  DEFLINE_POS_Cell_type,
  DEFLINE_POS_Chemovar,
  DEFLINE_POS_Chromosome,
  DEFLINE_POS_Clone,
  DEFLINE_POS_Clone_lib,
  DEFLINE_POS_Collected_by,
  DEFLINE_POS_Collection_date,
  DEFLINE_POS_Common,
  DEFLINE_POS_Country,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Culture_collection,
  DEFLINE_POS_Dev_stage,
  DEFLINE_POS_Ecotype,
  DEFLINE_POS_Endogenous_virus_name,
  DEFLINE_POS_Environmental_sample,
  DEFLINE_POS_Forma,
  DEFLINE_POS_Forma_specialis,
  DEFLINE_POS_Frequency,
  DEFLINE_POS_Genotype,
  DEFLINE_POS_Germline,
  DEFLINE_POS_Group,
  DEFLINE_POS_Haplogroup,
  DEFLINE_POS_Haplotype,
  DEFLINE_POS_Specific_host,
  DEFLINE_POS_Identified_by,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Isolation_source,
  DEFLINE_POS_Lab_host,
  DEFLINE_POS_Lat_lon,
  DEFLINE_POS_Linkage_group,
  DEFLINE_POS_Map,
  DEFLINE_POS_Mating_type,
  DEFLINE_POS_Metagenomic,
  DEFLINE_POS_Note_orgmod,
  DEFLINE_POS_Note_subsrc,
  DEFLINE_POS_Pathovar,
  DEFLINE_POS_Plasmid_name,
  DEFLINE_POS_Plastid_name,
  DEFLINE_POS_Pop_variant,
  DEFLINE_POS_Rearranged,
  DEFLINE_POS_Segment,
  DEFLINE_POS_Serogroup,
  DEFLINE_POS_Serotype,
  DEFLINE_POS_Serovar,
  DEFLINE_POS_Sex,
  DEFLINE_POS_Specimen_voucher,
  DEFLINE_POS_Strain,
  DEFLINE_POS_Subclone,
  DEFLINE_POS_Subgroup,
  DEFLINE_POS_Sub_species,
  DEFLINE_POS_Substrain,
  DEFLINE_POS_Subtype,
  DEFLINE_POS_Synonym,
  DEFLINE_POS_Teleomorph,
  DEFLINE_POS_Tissue_lib,
  DEFLINE_POS_Tissue_type,
  DEFLINE_POS_Transgenic,
  DEFLINE_POS_Type,
  DEFLINE_POS_Variety
} DefLinePos;

NLM_EXTERN Int4 GetDeflinePosForFieldType (ValNodePtr field);

/* ModifierItemGlobalData is used to store information about the available
 * modifiers - the name to use when displaying a list of checkboxes, whether
 * the modifier is an Organism modifier or a Source modifier, the subtype
 * to use when looking for the modifier in the organism qualifier list,
 * and whether this modifier is required by default.
 */
typedef struct modifieritemglobal {
  CharPtr       name;
  Boolean       isOrgMod;
  Uint1         subtype;
} ModifierItemGlobalData, PNTR ModifierItemGlobalPtr;

extern ModifierItemGlobalData DefLineModifiers[];
NLM_EXTERN size_t NumDefLineModifiers (void);

/* OrganismDescriptionModifiers is used to apply specific user preferences
 * for how to construct the organism descriptions - whether or not to use
 * labels for the modifiers, whether and how to limit the number of modifiers
 * used in any one organism description, whether to keep or remove modifier
 * information in parentheses in the organism taxonomy name, and whether or
 * not to apply modifiers to organisms with "sp." in the middle of the
 * taxonomy name.
 */
typedef struct organismdescriptionmodifiers {
  Boolean  use_labels;
  Int2     max_mods;
  Boolean  keep_paren;
  Boolean  exclude_sp;
  Boolean  exclude_cf;
  Boolean  exclude_aff;
  Boolean  exclude_nr;
  Boolean  include_country_extra;
  Int4     clone_isolate_HIV_rule_num;
  Boolean  use_modifiers;
  Boolean  allow_semicolon_in_modifier;
} OrganismDescriptionModifiers, PNTR OrganismDescriptionModifiersPtr;

/* These values are used for the clone_isolate_HIV_rule_num value in OrganismDescriptionModifiers */
typedef enum {
  clone_isolate_HIV_rule_prefer_clone = 1,
  clone_isolate_HIV_rule_prefer_isolate,
  clone_isolate_HIV_rule_want_both
} clone_isolate_HIV_rule_values;
  

typedef struct sourcequaldesc
{
  CharPtr       name;
  Boolean       isOrgMod;
  Uint1         subtype;
  Uint1         subfield;
} SourceQualDescData, PNTR SourceQualDescPtr;

NLM_EXTERN void SetRequiredModifiers (ModifierItemLocalPtr modList);
NLM_EXTERN void CountModifiers (ModifierItemLocalPtr ItemList, SeqEntryPtr sep);
NLM_EXTERN ValNodePtr FindBestModifiersEx(SeqEntryPtr sep, ModifierItemLocalPtr ItemList, Boolean use_new);
NLM_EXTERN ValNodePtr FindBestModifiers(SeqEntryPtr sep, ModifierItemLocalPtr ItemList);
NLM_EXTERN ValNodePtr FindBestModifiersForDeflineClauseList (ValNodePtr defline_clauses, ModifierItemLocalPtr ItemList);

NLM_EXTERN ValNodePtr GetModifierIndicesFromModList (ModifierItemLocalPtr modList);
extern void TestFindBestQualCombo (FILE *fp);


NLM_EXTERN CharPtr MergeValNodeStrings (ValNodePtr list, Boolean useReturn);

NLM_EXTERN void BuildDefLineFeatClauseList 
( SeqEntryPtr sep,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list);

NLM_EXTERN Boolean AreFeatureClausesUnique (ValNodePtr list);
NLM_EXTERN void DefLineFeatClauseListFree (ValNodePtr vnp);

NLM_EXTERN void 
BuildDefinitionLinesFromFeatureClauseLists 
(ValNodePtr list,
 ModifierItemLocalPtr modList,
 ValNodePtr modifier_indices,
 OrganismDescriptionModifiersPtr odmp);

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
 Boolean gene_cluster_opp_strand);

NLM_EXTERN Boolean IsSpName (CharPtr taxName);

#define DEFAULT_ORGANELLE_CLAUSE 10
NLM_EXTERN BioSourcePtr GetBiopForBsp (BioseqPtr bsp);
NLM_EXTERN Boolean IsLocAInBonSameStrand (SeqLocPtr slp1, SeqLocPtr slp2);
NLM_EXTERN void CleanUpTaxName (CharPtr taxName, Boolean keep_in_paren);
NLM_EXTERN Boolean UseOrgModifier (OrgModPtr mod, CharPtr   taxName);
NLM_EXTERN void AddModifierLabel 
( Boolean use_labels,
  Boolean is_orgmod,
  Uint1   subtype,
  CharPtr modifier_text);
NLM_EXTERN Boolean LIBCALLBACK IsMobileElement (SeqFeatPtr sfp);
NLM_EXTERN void RemoveNucProtSetTitles (SeqEntryPtr sep);


NLM_EXTERN ValNodePtr ReadTabTableFromFile (FILE *fp);
NLM_EXTERN ValNodePtr FreeTabTable (ValNodePtr row_list);
NLM_EXTERN ValNodePtr CountTabTableBlanks (ValNodePtr row_list);
NLM_EXTERN ValNodePtr ScanTabTableForSpecialCharacters (ValNodePtr row_list);
NLM_EXTERN void RemoveQuotesFromTabTable (ValNodePtr row_list);

NLM_EXTERN void SpecialCharFindWithContext (CharPtr PNTR strp, Pointer userdata, BoolPtr did_find, BoolPtr did_change);
NLM_EXTERN ValNodePtr FreeContextList (ValNodePtr context_list);

NLM_EXTERN Int4 ExtendSeqLocToEnd (SeqLocPtr slp, BioseqPtr bsp, Boolean end5);

/* functions for converting features */
NLM_EXTERN Boolean IsBioseqSetInGPS (BioseqSetPtr bssp);
NLM_EXTERN Boolean IsBioseqInGPS (BioseqPtr bsp);
NLM_EXTERN Boolean IsFeatInGPS (SeqFeatPtr sfp);
NLM_EXTERN void 
ApplyCDSOptionsToFeature
(SeqFeatPtr sfp,
 Boolean remove_mRNA,
 Boolean remove_gene,
 Boolean remove_transcript_id,
 Boolean keep_original);

NLM_EXTERN Boolean 
ConvertCDSToRNA 
(SeqFeatPtr  sfp,
 Uint2       rna_type);

NLM_EXTERN Boolean ConvertGeneToRNA (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertBioSrcToRepeatRegion (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN CharPtr SubSourceText (BioSourcePtr biop, Uint1 subtype, BoolPtr found);
NLM_EXTERN CharPtr OrgModText (BioSourcePtr biop, Uint1 subtype, BoolPtr found);
NLM_EXTERN CharPtr NoteText (BioSourcePtr biop, CharPtr comment);
NLM_EXTERN Boolean ConvertNonPseudoCDSToMiscFeat (SeqFeatPtr sfp, Boolean viral);

NLM_EXTERN CharPtr GetImportFeatureName (Uint2 featdef_key);
NLM_EXTERN RnaRefPtr RnaRefFromLabel (Uint2 featdef_to, CharPtr label, BoolPtr add_label_to_comment);

NLM_EXTERN SeqLocPtr GetProteinLocationForNucleotideFeatureConversion (SeqLocPtr nuc_slp, BoolPtr no_cds);
NLM_EXTERN SeqLocPtr FindNucleotideLocationForProteinFeatureConversion (SeqLocPtr slp);
NLM_EXTERN SeqLocPtr BuildProtLoc (SeqFeatPtr overlapping_cds, SeqLocPtr slp, Int4Ptr frame);
NLM_EXTERN Boolean ConvertImpToProtFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertProtToImpFunc (SeqFeatPtr  sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertRegionToProtFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertRegionToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertImpToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertRegionToRNAFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertGeneToMiscFeatFunc (SeqFeatPtr sfp, Uint2 featdef_to);
NLM_EXTERN Boolean ConvertProtToProtFunc (SeqFeatPtr sfp, Uint2 featdef_to);


NLM_EXTERN void InstantiateMatPeptideProducts (SeqEntryPtr sep);

NLM_EXTERN Boolean CodingRegionHasTranslExcept (SeqFeatPtr sfp);

NLM_EXTERN SeqEntryPtr SequenceStringToSeqEntry (CharPtr str, SeqIdPtr sip, Uint1 mol_type);

NLM_EXTERN void RevCompOneFeatForBioseq (SeqFeatPtr sfp, BioseqPtr bsp);

NLM_EXTERN ValNodePtr SplitPubListFree (ValNodePtr list);
NLM_EXTERN ValNodePtr MakeSplitPubListFromTabList (ValNodePtr PNTR tab_table, SeqEntryPtr sep, ValNodePtr PNTR err_list);
NLM_EXTERN void SplitPubsByList (ValNodePtr split_list);
/* for parsing collection dates */
NLM_EXTERN CharPtr ReformatDateStringEx (CharPtr orig_date, Boolean month_first, BoolPtr month_ambiguous);
NLM_EXTERN Int4 GetYearFromToken (CharPtr token, Int4 token_len);
NLM_EXTERN Int4 ReadNumberFromToken (CharPtr token, Int4 token_len);
NLM_EXTERN CharPtr GetMonthFromToken (CharPtr token, Int4 token_len);
NLM_EXTERN Int4 GetMonthNumFromAbbrev (CharPtr month_abbrev);
NLM_EXTERN CharPtr GetMonthAbbrev (Int4 n);
NLM_EXTERN Int4 GetDaysInMonth (Int4 n);

NLM_EXTERN ValNodePtr CreateStructuredCommentsFromFile (FILE *fp, SeqEntryPtr sep, Boolean apply_to_all);

#define ALNMGR_GAP           -2
#define ALNMGR_ROW_UNDEFINED -1

NLM_EXTERN void 
AlignmentIntervalToString 
(SeqAlignPtr salp,
 Int4        row,
 Int4        start,
 Int4        stop,
 Int4        target_row,
 Boolean     view_whole_entity,
 Uint1Ptr    seqbuf,
 Uint1Ptr    alnbuf,
 Int4 PNTR   alnbuffer_len,
 Boolean     show_substitutions);

/* Compare Functions */
NLM_EXTERN int LIBCALL ObjectIdCompare (ObjectIdPtr a, ObjectIdPtr b);
NLM_EXTERN int LIBCALL DbtagCompare (DbtagPtr a, DbtagPtr b);
NLM_EXTERN int LIBCALL OrgModSetCompare (OrgModPtr mod1, OrgModPtr mod2);
NLM_EXTERN int LIBCALL OrgNameCompare (OrgNamePtr onp1, OrgNamePtr onp2);
NLM_EXTERN int LIBCALL OrgRefCompare (OrgRefPtr orp1, OrgRefPtr orp2);

extern void CountNsInSequence (BioseqPtr bsp, Int4Ptr p_total, Int4Ptr p_max_stretch, Boolean expand_gaps);
NLM_EXTERN Boolean IsTSA (BioseqPtr bsp);

NLM_EXTERN Boolean IsPseudo (SeqFeatPtr sfp);

NLM_EXTERN Boolean ExtendPartialsToEndOrGap (SeqFeatPtr sfp);
NLM_EXTERN Boolean RetranslateOneCDS 
( SeqFeatPtr sfp,
  Uint2 entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds);
NLM_EXTERN SeqFeatPtr FindBestProtein (Uint2 entityID, SeqLocPtr product);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* ndef _SQNUTILS_ */

