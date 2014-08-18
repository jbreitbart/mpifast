static char const rcsid[] = "$Id: posit.c,v 6.87 2008/10/03 18:09:17 madden Exp $";

/* $Id: posit.c,v 6.87 2008/10/03 18:09:17 madden Exp $
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

   File name: posit.c

  Author: Alejandro Schaffer

  Contents: utilities for position-based BLAST.

  $Revision: 6.87 $ 
 *****************************************************************************

 * $Log: posit.c,v $
 * Revision 6.87  2008/10/03 18:09:17  madden
 * Change a few constants for pseudo counts (from A. Schaffer)
 *
 * Revision 6.86  2008/07/22 19:45:32  kans
 * removed unused variables
 *
 * Revision 6.85  2008/03/31 16:21:58  kans
 * CodeWarrior requred a cast to Blast_GetMatrixBackgroundFreq result
 *
 * Revision 6.84  2008/03/31 13:36:10  madden
 * Implemented a new method to compute effective observations.
 * Implemented a new entropy-based method to compute column-specific pseudocounts.
 *
 * Revision 6.83  2007/05/07 13:30:54  kans
 * added casts for Seq-data.gap (SeqDataPtr, SeqGapPtr, ByteStorePtr)
 *
 * Revision 6.82  2007/03/12 16:21:35  papadopo
 * From Kristoffer Osowski: fix a buffer overflow caused by incorrect computation of number of sequences for psiblast
 *
 * Revision 6.81  2007/01/22 19:20:55  camacho
 * From Alejandro Schaffer:
 * In posPurgeMatches, when in command-line mode, added a warning for the
 * situation in which only the query is used to construct the PSSM.
 *
 * Revision 6.80  2006/09/28 12:45:39  madden
 *        Use more accurrate, symmetric, frequency data for the PAM matrices.
 *        Make the frequency ratios for the X character in BLOSUM90 the
 *        same consistent with the other BLOSUM matrices.
 *        [from Mike Gertz]
 *
 * Revision 6.79  2006/09/25 19:27:39  madden
 *    - Include blastkar.h in posit.c to get the declarations of
 *      BlastScoreFreqNew and BlastScoreFreqDestruct
 *    - Change %lf in printf statements to %f: %lf is equivalent to %f,
 *      but is not ANSI C, so was causing warnings.
 *    - Copied the correct, 28-character frequency ratios from
 *      algo/blast/core/matrix_freq_ratios.c.
 *    [from Mike Gertz]
 *
 * Revision 6.78  2006/08/01 16:52:07  camacho
 * Fix comments
 *
 * Revision 6.77  2006/06/30 18:44:40  camacho
 * Add support for O and J to getRes, ResToInt
 *
 * Revision 6.76  2005/08/30 18:21:02  coulouri
 * From Mike Gertz:
 *    In posFreqsToMatrix while setting posSearch->posPrivateMatrix, check
 *    whether the frequency ratio is zero before taking the log.  If it is
 *    zero, set the score to BLAST_SCORE_MIN.  Do not test whether the
 *    original matrix has a score of BLAST_SCORE_MIN.
 *
 * Revision 6.75  2005/08/05 12:05:13  coulouri
 * From Mike Gertz:
 * - Changed the freq ratio of a (*,-) or (-,*) match for all matrices to
 *   be zero.  The effect of this change is that these substitutions get
 *   assigned a score of BLAST_SCORE_MIN when composition-based statistics
 *   is used in any mode.
 *
 * Revision 6.74  2005/07/28 14:57:10  coulouri
 * remove dead code
 *
 * Revision 6.73  2005/06/08 12:26:29  camacho
 * Change posReadPosFreqsScoremat to accept other encodings of Bioseq data in ASN.1 PSSM
 *
 * Revision 6.72  2004/10/12 15:06:57  papadopo
 * 1. Modify residue frequency IO to comply with new scoremat spec
 * 2. Remove check that residue frequencies read from scoremat are <= 1.0
 * 3. Pass gap open and gap extend penalties into BposComputation and
 * 	CposComputation, so that scoremats can contain them
 *
 * Revision 6.71  2004/08/23 17:09:22  papadopo
 * From Michael Gertz: move static arrays out of header and into the one file that needs them
 *
 * Revision 6.70  2004/08/05 17:30:53  camacho
 * Remove initialization of identifier as it is no longer required
 *
 * Revision 6.69  2004/07/24 18:56:01  camacho
 * Fix in posDemographics when GetSequenceWithDenseSeg cannot find sequence data
 *
 * Revision 6.68  2004/07/19 17:13:13  papadopo
 * add capability to perform input and output of residue frequencies in scoremat form; also call PSIMatrixFrequencyRatiosNew before restarting from checkpoint
 *
 * Revision 6.67  2004/07/13 13:54:15  camacho
 * Fix memory leak
 *
 * Revision 6.66  2004/06/25 21:54:51  dondosha
 * Choose ideal values for lambda and K correctly for ungapped search
 *
 * Revision 6.65  2004/06/23 14:53:29  camacho
 * Copy renamed versions of SFreqRatios and its *{New,Free} functions to avoid
 * dependency ncbitool -> blast
 *
 * Revision 6.64  2004/06/22 14:16:46  camacho
 * Changed signature of posFreqsToMatrix, added use of SFreqRatios structure from
 * algo/blast/core/ to obtain underlying matrices' frequency ratios.
 * This change results in using the frequency ratios to provide the scores
 * for the PSSM in columns where all residue frequencies are 0. Previously the
 * standard scoring matrix were used.
 *
 * Revision 6.63  2004/06/08 14:03:48  camacho
 * Alejandro Schaffer's fix to spread out gap costs in posDemographics.
 *
 * Revision 6.62  2004/05/14 12:13:09  camacho
 * Made posDemographics non-static for testing purposes.
 *
 * Revision 6.61  2003/08/04 20:43:55  dondosha
 * Test for selenocysteines when comparing checkpoint sequences with query
 *
 * Revision 6.60  2003/05/30 17:25:37  coulouri
 * add rcsid
 *
 * Revision 6.59  2003/05/13 16:02:53  coulouri
 * make ErrPostEx(SEV_FATAL, ...) exit with nonzero status
 *
 * Revision 6.58  2001/12/11 14:48:54  madden
 * Fix for ABW (reset Xcount to zero in some cases)
 *
 * Revision 6.57  2001/08/29 19:04:48  madden
 * added parameter posComputationCalled to outputPosComputation, extra printing added in revision 6.54 is suppressed if posComputationCalled is FALSE
 *
 * Revision 6.56  2001/08/06 18:09:13  madden
 * Corrected handling of X in posCancel by adding usage of Xcount
 *
 * Revision 6.55  2001/04/09 13:00:09  madden
 * Fixed error in posComputeExtents; adjustment of interval sizes when the query contained an X had been asymmetric.
 *
 * Revision 6.54  2001/04/03 19:38:24  madden
 * Changed IDENTITY_PERCENTAGE to 0.94, Added to output of -Q option in outputPosMatrix
 *
 * Revision 6.53  2001/02/16 16:11:50  dondosha
 * In WposComputation, compute posMatrix from posFreqs if seqalign argument is NULL
 *
 * Revision 6.52  2001/01/03 01:49:38  bauer
 * Changed from static to "LIBCALL":
 *  posAllocateMemory
 *  posPurgeMatches
 *  posCancel
 *  posComputeExtents
 *  posComputeSequenceWeights
 *  posCheckWeights
 *  posComputePseudoFreqs
 *  posScaling
 *
 * Revision 6.51  2000/11/24 22:07:51  shavirin
 * Fixed some memory leaks.
 *
 * Revision 6.50  2000/11/20 14:35:51  madden
 * Changed FileOpen mode for byte-encoded checkpoint files from "r" to "rb" or from "w" to "wb" to solve a problem on Windows NT.
 *
 * Revision 6.49  2000/11/13 14:00:39  madden
 * Added frequency ratios for * in all standard matrices
 *
 * Revision 6.48  2000/11/09 14:27:52  madden
 * psi-blast fixes for star character
 *
 * Revision 6.47  2000/11/01 16:25:57  madden
 * Changes from Futamura for psitblastn
 *
 * Revision 6.46  2000/10/24 16:28:29  madden
 * Changed IDENTITY_RATIO for putging near-identical matches from 0.98 to 0.95
 *
 * Revision 6.45  2000/10/10 21:46:04  shavirin
 * Added support for BLOSUM50, BLOSUM90, PAM250 with -t T
 *
 * Revision 6.44  2000/08/18 21:28:37  madden
 * support for BLOSUM62_20A and BLOSUM62_20B
 *
 * Revision 6.43  2000/07/31 16:41:01  shavirin
 * Reduced POSIT_SCALE_FACTOR from 1000 to 200 to avoid overflow
 * with BLOSUM80; moved declaration os POSIT_SCALE_FACTOR to posit.h
 *
 * Revision 6.42  2000/07/26 13:11:19  shavirin
 * Added magical "LIBCALL" to pacify Windows build function allocatePosFreqs()
 *
 * Revision 6.41  2000/07/25 18:12:04  shavirin
 * WARNING: This is no-turning-back changed related to S&W Blast from
 * Alejandro Schaffer
 *
 * Revision 6.40  2000/05/01 12:48:33  madden
 * changed rules for gaps in processing alignments with -B
 *
 * Revision 6.39  2000/03/02 21:47:07  shavirin
 * Added missing variable for POSIT_DEBUG case
 *
 * Revision 6.38  1999/12/16 19:18:00  egorov
 * Code cleanup
 *
 * Revision 6.37  1999/11/16 21:33:46  shavirin
 * Fixed bug involved posSearch->posResultSequences structure.
 *
 * Revision 6.36  1999/11/16 17:30:41  shavirin
 * Added copying use_best_align parameter in copySearchItems.
 *
 * Revision 6.35  1999/11/15 22:21:09  shavirin
 * Removed nested comments in log space.
 *
 * Revision 6.34  1999/11/15 21:48:31  shavirin
 * Added possibility to take into account best alignments of all
 * alignments (even with e-thresh value larger than threshhold)
 *
 * Revision 6.33  1999/10/21 16:15:04  shavirin
 * Removed unused array and all references to array threshSequences
 *
 * Revision 6.32  1999/09/30 14:15:29  shavirin
 * Fixed bug in the function findThreshSequences().
 *
 * Revision 6.31  1999/09/23 20:58:37  shavirin
 * Fixed some memory leaks.
 *
 * Revision 6.27  1999/09/03 17:24:38  madden
 * Eliminated use of posMaxThresh field in posSearchItems.  Recoded findThreshSequences completely and changed CposComputation so as not to use search->result_struct and thereby eliminate the hidden assumption that search->result_struct and listOfSeqAligns have the matches listed in the same order
 *
 * Revision 6.26  1999/08/04 13:27:10  madden
 * Added -B option
 *
 * Revision 6.25  1999/04/05 14:45:40  madden
 * Fixed format mismatches
 *
 * Revision 6.24  1999/03/21 19:41:51  madden
 * Added 3rd argument matrixfp to outputPosMatrix, Took some of the code in outputPosMatrix outside the #ifdef POSIT_DEBUG for use with -Q option
 *
 * Revision 6.23  1999/01/26 18:27:58  madden
 * Made functions public for AS
 *
 * Revision 6.22  1998/12/09 18:51:51  madden
 * fixed counting bug in posCancel
 *
 * Revision 6.21  1998/09/28 12:31:31  madden
 * Used BlastConstructErrorMessage
 *
 * Revision 6.20  1998/09/09 21:18:33  madden
 * AS fixed warnings
 *
 * Revision 6.19  1998/09/09 16:09:20  madden
 * Changes for PHI-BLAST
 *
 * Revision 6.18  1998/08/26 18:07:00  kans
 * fixed -v -fd warnings (AS)
 *
 * Revision 6.17  1998/06/18 18:20:22  madden
 * Fixed typo in posConvergenceTest
 *
 * Revision 6.16  1998/06/14 19:43:02  madden
 * Added function posFreqsToInformation
 *
 * Revision 6.15  1998/06/12 20:38:48  madden
 * Fix for no hits to build model situation
 *
 * Revision 6.14  1998/06/09 19:38:16  madden
 * Changes rounddown to posit_rounddown to avoid conflict
 *
 * Revision 6.13  1998/04/24 19:29:30  madden
 * Moved rescaling code to blastool.c
 *
 * Revision 6.12  1998/03/25 22:36:17  egorov
 * Change type of posRepeatSequences
 *
 * Revision 6.11  1998/03/23 18:32:30  madden
 * Fix for zero/zero problem
 *
 * Revision 6.10  1998/02/06 18:34:17  madden
 * Added check that residue was not masked in posReadCheckpoint
 *
 * Revision 6.9  1998/02/03 15:57:28  madden
 * Cpos arg in posComputePseudoFreqs set to FALSE for WposComputation call
 *
 * Revision 6.8  1998/01/02 22:19:46  madden
 * Replaced printf by ErrPostEx of SEV_WARNING
 *
 * Revision 6.7  1997/12/23 21:07:06  madden
 * Changes for checkpointing
 *
 * Revision 6.6  1997/12/12 22:14:35  kans
 * changed round to rounddown to avoid CodeWarrior 68K collision
 *
 * Revision 6.5  1997/11/19 15:29:31  madden
 * Changed OS_UNIX ifdef to POSIT_DEBUG
 *
 * Revision 6.4  1997/09/08 13:33:29  madden
 * Added posEpsilon2 to check for small numbers
 *
 * Revision 6.3  1997/09/05 22:29:13  madden
 * Check for zero denominator and replace log(2) by NCBIMATH_LN2
 *
 * Revision 6.2  1997/09/02 22:23:01  madden
 * Removed redundant calls to updateLambdaK
 *
 * Revision 6.1  1997/08/27 21:18:18  madden
 * Fixed problem with deleted matrix
 *
 * Revision 6.0  1997/08/25 18:53:48  madden
 * Revision changed to 6.0
 *
 * Revision 1.22  1997/08/20 21:35:01  madden
 * ALL_ROUNDS replaced by Boolean
 *
 * Revision 1.21  1997/08/11 15:45:24  madden
 * eliminated obsolete fields
 *
 * Revision 1.20  1997/07/28 18:35:06  madden
 * Removed and ifdefed printf
 *
 * Revision 1.19  1997/06/27 19:14:01  madden
 * Fixed two bugs in posComputeSequenceWeights for the special case where all participating sequences are identical in a block
 *
 * Revision 1.17  1997/06/23 14:42:46  madden
 *  Made posComputeSequenceWeights faster by catching the special case where the set of participating sequences does not change from one column to the next.
 *
 * Revision 1.16  1997/05/29 20:35:47  madden
 * Eliminate duplicate sequences and alignments that are 98 perc. identical and ignore columns with all identical sequence weights.
 *
 * Revision 1.15  1997/05/27 20:26:09  madden
 * Fixed problem with matrix
 *
 * Revision 1.14  1997/05/23 20:52:50  madden
 * Fixed bug in setting of matrix for psi-blast
 *
 * Revision 1.13  1997/05/22 21:25:28  madden
 * fixed memory leaks
 *
 * Revision 1.12  1997/05/16 20:56:35  madden
 * replace hard coded numbers by defines
 *
 * Revision 1.11  1997/05/16 20:09:42  madden
 * Fixes for statistical problems
 *
 * Revision 1.10  1997/05/07 21:00:03  madden
 * Call to SeqId2OrdinalId replaces call to readdb_gi2seq
 *
 * Revision 1.9  1997/05/01 15:53:25  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.8  1997/04/23  13:31:20  madden
 * Changed diagnostic output.
 *
 * Revision 1.7  1997/04/22  16:36:49  madden
 * Changes for use of psi-blast with www.
 *
 * Revision 1.6  1997/04/10  19:25:53  madden
 * Added casts, COMMAND_LINE replaced by ALL_ROUNDS.
 *
 * Revision 1.5  1997/04/09  20:01:53  madden
 * Functions CposComputation and WposComputation replace posComputations.
 *
 * Revision 1.4  1997/04/04  20:44:55  madden
 * Changed posComputation to return Int4Ptr *.
 *
 * Revision 1.3  1997/03/27  22:30:51  madden
 * Fix for Array Bounds Read.
 *
 * Revision 1.2  1997/03/11  14:38:40  madden
 * Fixes for GI's instead of ordinal numbers.
 *
 * Revision 1.1  1997/02/13  15:22:13  madden
 * Initial revision
 *
*/


#include <ncbi.h>
#include <blastpri.h>
#include <objcode.h>
#include <objseq.h>
#include <objsset.h>
#include <objscoremat.h>
#include <sequtil.h>
#include <posit.h>
#include <txalign.h>
#include <blastkar.h>

#include <algo/blast/composition_adjustment/nlm_linear_algebra.h>
#include <algo/blast/composition_adjustment/matrix_frequency_data.h>
#include <algo/blast/composition_adjustment/composition_adjustment.h>
#include <algo/blast/composition_adjustment/compo_heap.h>
#include <algo/blast/composition_adjustment/smith_waterman.h>
#include <algo/blast/composition_adjustment/redo_alignment.h>
#include <algo/blast/composition_adjustment/unified_pvalues.h>

/*small constants to test against 0*/
#define posEpsilon 0.0001
#define posEpsilon2 0.0000001
/*Representation of a gap in a motif*/
#define GAP_CHAR 0
/*Used inside a seqAlign to reprsent the presence of a gap*/
#define GAP_HERE (-1)
/*Used to check that diagnostics printing routine will work*/
#define EFFECTIVE_ALPHABET 20

#define POSIT_PERCENT 0.05
#define POSIT_NUM_ITERATIONS 10


#define POS_RESTING 0
#define POS_COUNTING 1

#define IDENTITY_RATIO 0.94


/*Allocate memory for  data structures inside posSearch used in
* position-specific caculations
* posSearch -- to be filled in 
* alphabetSize -- number of distinct characters used in the sequences
* querySize -- number of characters in the query sequence
* numSequences -- number of matching sequences potentially in the model */
void LIBCALL posAllocateMemory(posSearchItems * posSearch, 
		       Int4 alphabetSize, Int4 querySize, Int4 numSequences)
{
  Int4 i, j;  /*loop indices*/

  posSearch->posCount = (Int4 *) MemNew(querySize * sizeof(Int4));
  if (NULL == posSearch->posCount)
    exit(EXIT_FAILURE);
  for(i = 0; i < querySize; i++)
    posSearch->posCount[i] = 0;

  posSearch->posC = (Int4 **) MemNew((querySize + 1) * sizeof(Int4 *));
  if (NULL == posSearch->posC)
    exit(EXIT_FAILURE);
  for(i = 0; i <= querySize; i++) {
    posSearch->posC[i] = (Int4 *) MemNew(alphabetSize * sizeof(Int4));
    if (NULL == posSearch->posC[i])
      exit(EXIT_FAILURE);   
    for(j = 0; j < alphabetSize; j++)
      posSearch->posC[i][j] = 0;
 
  }
  posSearch->posDistinctDistrib = (Int4 **) MemNew((querySize + 1) * sizeof(Int4* ));
  for (i = 0; i <=querySize; i++)  {
    posSearch->posDistinctDistrib[i] = (Int4 *) MemNew((EFFECTIVE_ALPHABET + 1) * sizeof(Int4));
    for(j = 0; j <=EFFECTIVE_ALPHABET; j++)
      posSearch->posDistinctDistrib[i][j] = 0;
  }
  posSearch->posNumParticipating = (Int4 *) MemNew((querySize + 1) * sizeof(Int4 ));
  posSearch->posGaplessColumnWeights = (Nlm_FloatHi *) MemNew((querySize + 1) * sizeof(Nlm_FloatHi));
  if (NULL == posSearch->posGaplessColumnWeights)
    exit(EXIT_FAILURE);
  posSearch->posMatchWeights = (Nlm_FloatHi **) MemNew((querySize+1) * sizeof(Nlm_FloatHi *));
  if (NULL == posSearch->posMatchWeights)
    exit(EXIT_FAILURE);
  for (i = 0; i <= querySize ; i++) {
    posSearch->posMatchWeights[i] = (Nlm_FloatHi *) MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    if (NULL == posSearch->posMatchWeights[i])
      exit(EXIT_FAILURE);
    for(j = 0; j < alphabetSize; j++) 
      posSearch->posMatchWeights[i][j] = 0.0;
  }  

  posSearch->posMatrix = (BLAST_Score **) MemNew((querySize + 1) * sizeof(BLAST_Score *));
  posSearch->posPrivateMatrix = (BLAST_Score **) MemNew((querySize + 1) * sizeof(BLAST_Score *));
  if (NULL == posSearch->posMatrix)
    exit(EXIT_FAILURE);
  for(i = 0; i <= querySize; i++) {
    posSearch->posMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
    posSearch->posPrivateMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
    if (NULL == posSearch->posMatrix[i])
      exit(EXIT_FAILURE);   
    for(j = 0; j < alphabetSize; j++)
      posSearch->posMatrix[i][j] = 0;
 
  }

  posSearch->posSigma = (Nlm_FloatHi *) MemNew((querySize) * sizeof(Nlm_FloatHi));
  if (NULL == posSearch->posSigma)
    exit(EXIT_FAILURE);
  for(i = 0; i < querySize; i++) {
    posSearch->posSigma[i] = 0.0;
  }

  posSearch->posIntervalSizes = (Int4 *) MemNew((querySize) * sizeof(Int4));
  if (NULL == posSearch->posIntervalSizes)
    exit(EXIT_FAILURE);
  for(i=0; i < querySize; i++)
    posSearch->posIntervalSizes[i] = 0;

  posSearch->posDescMatrixLength = numSequences;
  posSearch->posDescMatrix = (posDesc **) MemNew((numSequences + 1) * sizeof(posDesc *));
  if (NULL == posSearch->posDescMatrix)
    exit(EXIT_FAILURE);
  for (i = 0; i <= numSequences; i++) {
    posSearch->posDescMatrix[i] = (posDesc *) MemNew(querySize * sizeof(posDesc));
    if (NULL == posSearch->posDescMatrix[i])
      exit(EXIT_FAILURE);
    for(j = 0; j < querySize; j++) {
      posSearch->posDescMatrix[i][j].letter = UNUSED;
      posSearch->posDescMatrix[i][j].used = FALSE;
      posSearch->posDescMatrix[i][j].e_value = 1.0;
      posSearch->posDescMatrix[i][j].leftExtent = -1;
      posSearch->posDescMatrix[i][j].rightExtent = querySize;
    }
  }  
  posSearch->posExtents = (posDesc *) MemNew(querySize * sizeof(posDesc));
  if (NULL == posSearch->posExtents)
    exit(EXIT_FAILURE);
  for(j = 0; j < querySize; j++) {
    posSearch->posExtents[j].used = FALSE;
    posSearch->posExtents[j].leftExtent = -1;
    posSearch->posExtents[j].rightExtent = querySize;
  }
   posSearch->posA = (Nlm_FloatHi *) MemNew((numSequences+ 1) * sizeof(Nlm_FloatHi));
   if (NULL == posSearch->posA)
     exit(EXIT_FAILURE);
   posSearch->posRowSigma = (Nlm_FloatHi *) MemNew((numSequences + 1) * sizeof(Nlm_FloatHi));
   if (NULL == posSearch->posRowSigma)
     exit(EXIT_FAILURE);

   /* populated in posComputePseudoFreqs or on demand */
   posSearch->stdFreqRatios = NULL; 
}

static void freePosFreqs(Nlm_FloatHi ** posFreqs, Int4 length)
{
  Int4 i;

  for (i = 0; i <= length; i++)
    MemFree(posFreqs[i]);
  MemFree(posFreqs); 
}

/*Deallocate memory allocated in posReadCheckpoint
* posSearch -- pointer to record used in building the position-specific model
* querySize -- number of characters in the query sequence
*/
void LIBCALL posCheckpointFreeMemory(posSearchItems *posSearch, Int4 querySize)
{
  Int4 i;  /*loop index*/

  freePosFreqs(posSearch->posFreqs, querySize);
  for(i = 0; i <= querySize; i++){
    MemFree(posSearch->posMatrix[i]);
    MemFree(posSearch->posPrivateMatrix[i]);
  }
  MemFree(posSearch->posMatrix);
  MemFree(posSearch->posPrivateMatrix);
}

/*Deallocate memory allocated in posAllocateMemory
* posSearch -- pointer to record used in building the position-specific model
* querySize -- number of characters in the query sequence
*/
static void posFreeMemory(posSearchItems *posSearch, Int4 querySize)
{
  Int4 i;  /*loop index*/

  MemFree(posSearch->posCount);
  MemFree(posSearch->posExtents);
  MemFree(posSearch->posSigma);

  for(i = 0; i <= querySize; i++){
    MemFree(posSearch->posC[i]);
    MemFree(posSearch->posMatrix[i]);
    MemFree(posSearch->posPrivateMatrix[i]);
    MemFree(posSearch->posMatchWeights[i]);
    MemFree(posSearch->posDistinctDistrib[i]);
  }

  MemFree(posSearch->posC);
  MemFree(posSearch->posDistinctDistrib);
  MemFree(posSearch->posNumParticipating);

  for(i = 0; i <= posSearch->posDescMatrixLength; i++)
    MemFree(posSearch->posDescMatrix[i]);

  MemFree(posSearch->posMatrix);
  MemFree(posSearch->posPrivateMatrix);
  MemFree(posSearch->posDescMatrix);
  MemFree(posSearch->posGaplessColumnWeights);
  MemFree(posSearch->posMatchWeights);
  MemFree(posSearch->posA);
  MemFree(posSearch->posRowSigma);
  MemFree(posSearch->posIntervalSizes);
  MemFree(posSearch->posUseSequences);
  posSearch->stdFreqRatios = PSIMatrixFrequencyRatiosFree(posSearch->stdFreqRatios);
  freePosFreqs(posSearch->posFreqs,querySize);
}

/*Cleanup position-specific  data structures after one pass*/
void LIBCALL posCleanup(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
  posFreeMemory(posSearch, compactSearch->qlength);
}


/*extract the e-value that applies to an entire dense
  diagonal alignment from its ScorePtr, based on similar
  code from Tom Madden*/

static Nlm_FloatHi getEvalueFromSeqAlign(SeqAlignPtr thisSeqAlign)
{
  ScorePtr thisScorePtr;

  thisScorePtr = thisSeqAlign->score;
  while ((thisScorePtr != NULL) &&
         (StringICmp(thisScorePtr->id->str, "e_value") != 0) &&
         (StringICmp(thisScorePtr->id->str, "sum_e") != 0))
    thisScorePtr = thisScorePtr->next;
  if(NULL == thisScorePtr)
    return(10.0);
  else
    return((Nlm_FloatHi) (thisScorePtr->value.realvalue));
}

/*Find the lowest e-value among all seqAligns for the sequence represented by
curSeqAlign*/
static Nlm_FloatHi minEvalueForSequence(SeqAlignPtr curSeqAlign, SeqAlignPtr listOfSeqAligns) 
{
   SeqAlignPtr testSeqAlign; /*Index into listOfSeqALigns*/
   DenseSegPtr curSegs, testSegs; /*Used to extract ids from curSeqAlign, testSeqAlign*/
   SeqIdPtr curId, testId; /*Ids of target sequences in testSeqAlign*/
   Nlm_FloatHi  returnValue; /*stores current best e-value*/
   Nlm_FloatHi  testEvalue; /*temporary e-value for one seqAlign*/
   Boolean seen;   /*have we seen a seqAlign matching the sequence yet*/

   returnValue = getEvalueFromSeqAlign(curSeqAlign);
   curSegs = (DenseSegPtr) curSeqAlign->segs;
   curId = curSegs->ids->next; 
   seen = FALSE;

   testSeqAlign = listOfSeqAligns;
   while (NULL != testSeqAlign) {
     testSegs = (DenseSegPtr) testSeqAlign->segs;

     if(testSegs->ids == NULL)
         break;

     testId = testSegs->ids->next; 
     if (SeqIdMatch(curId, testId)) {
         seen = TRUE;
	 if ((testEvalue = getEvalueFromSeqAlign(testSeqAlign)) < returnValue)
	   returnValue = testEvalue;
       }
     else
      /*if we have already seen a match and this one doesn't match,
        then stop looking*/
       if (seen)  
	 break;
     testSeqAlign = testSeqAlign->next;
   }
   return(returnValue);
}


/*Count the number of seqAligns in a list (returned) and count the number of
distinct target sequences represented (passed back in numSequences);
if useThreshold is TRUE, only those sequences with e-values below the threshold are counted.
Important assumption: All SeqAligns with  the same target sequence
are consecutive in the list*/
static Int4 countSeqAligns(SeqAlignPtr listOfSeqAligns, Int4 * numSequences, Boolean useThreshold, Nlm_FloatHi threshold)
{
    SeqAlignPtr curSeqAlign, prevSeqAlign;
    Int4 seqAlignCounter;
    DenseSegPtr curSegs;
    SeqIdPtr curId, prevId; /*Ids of target sequences in current and previous
                              SeqAlign*/
    
    seqAlignCounter = 0;
    *numSequences = 0;
    curSeqAlign = listOfSeqAligns;
    prevSeqAlign = NULL;
    while (NULL != curSeqAlign) {
        curSegs = (DenseSegPtr) curSeqAlign->segs;
        if(curSegs->ids == NULL)
            break;
        curId = curSegs->ids->next; 
        seqAlignCounter++;
        if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(curId, prevId))))
            if (!useThreshold || (threshold > minEvalueForSequence(curSeqAlign, listOfSeqAligns)))
                (*numSequences)++;
        prevSeqAlign = curSeqAlign;
        prevId = curId;
        curSeqAlign = curSeqAlign->next;
    }
    return(seqAlignCounter);
}

/*Find which sequences that match in the i-th pass have an e-value below
  the specified threshold. These sequences will be used to make the
  score matrix for the next pass*/
static void findThreshSequences(posSearchItems *posSearch, BlastSearchBlkPtr search, SeqAlignPtr listOfSeqAligns, Int4 numalign, Int4 numseq)
{

    Int4 alignIndex; /* indices for sequences and alignments*/
    SeqAlignPtr curSeqAlign, prevSeqAlign; /* pointers into list of seqAligns*/
    DenseSegPtr curSegs;  /*Item in list of seqAligns*/
    SeqIdPtr thisId, prevId; /*Ids of target sequences in current and previous
                               SeqAlign*/  
    Nlm_FloatHi thisEvalue; /*Best E-value for current sequence*/
    Int4 ordinalNumber; /*index of sequence within database*/
    
    /*Allocate boolean array to store values*/
    posSearch->posResultSequences = (Int4 *) MemNew(numseq * sizeof(Int4));
    posSearch->posResultsCounter = 0;

    curSeqAlign = listOfSeqAligns;
    prevSeqAlign = NULL;
    for(alignIndex = 0; alignIndex < numalign; alignIndex++) {
        curSegs = (DenseSegPtr) curSeqAlign->segs;
        thisId = curSegs->ids->next;
        if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(thisId, prevId)))) {
            thisEvalue = minEvalueForSequence(curSeqAlign, curSeqAlign);
            thisId = curSegs->ids->next;  /*id of target sequence is second*/
            /*Get ordinal ids of sequences in result*/
            ordinalNumber = SeqId2OrdinalId(search->rdfp, thisId);
            if(thisEvalue < (search->pbp->ethresh)) {
                posSearch->posResultSequences[posSearch->posResultsCounter] =
                    ordinalNumber;
                posSearch->posResultsCounter++;
            }
        }
        prevSeqAlign = curSeqAlign;
        prevId = thisId;
        curSeqAlign = curSeqAlign->next;
    }
}



/* Determines if  the search has converged from round to round.
*  Checks whether every new sequence found is in posSearch->posResultSequences.
*  Also sets up posSearch->posRepeatSequences, a boolean array that
*  indicates whether the sequence represented by the i-th new seqAlign is a repeat.
*  This is used in printing the sequences where they are
*  subdivided into two categories: sequences that were found previously
*  and new sequences.
*  posSearch is the data structure representing the parameters of the position-specific part
*  search represents the overall  BLAST search
*  listOfSeqAligns is one representation of the results of the current round.
*  If thissPassNum is 1, then it checks only to see that some sequence
* distinct from the query was found */
void LIBCALL posConvergenceTest(posSearchItems *posSearch, BlastSearchBlkPtr search, SeqAlignPtr listOfSeqAligns, Int4 thisPassNum)
{
  Int4 numseq;   /*Number of sequences found*/
  Int4 numalign; /* Number of items in listOfSeqAligns*/
  Int4 oldSeqIndex; /*Ordinal number of a sequence in old results (previous round)*/
  Int4 alignIndex; /*index into the list of seqAligns*/
  Boolean found;  /*Have we found the new sequence on the old list?*/
  SeqAlignPtr curSeqAlign, prevSeqAlign, startAlign; /* pointers into list of seqAligns*/
  DenseSegPtr curSegs;  /*Item in list of seqAligns*/
  SeqIdPtr thisId, prevId; /*Ids of target sequences in current and previous
			   SeqAlign*/
  Int4 ordinalNumber; /*Ordinal number of a sequence in the database*/
  Nlm_FloatHi thisEvalue; /*lowest evalue from all seqAligns for a sequence*/
  Int4 queryOffset, subjectOffset, retrievalOffset; /*offsets needed
                                                    to make a match align*/
  Int4 qplace, splace; /*index into query string and matching string*/
  Uint1Ptr q,s; /*Pointers into query and matching string*/
  Int4 queryLength; /*length of query*/
  Int4 matchLength; /* length of match*/
  Int4 subjectLength; /* length of a matching string*/
  Int4 c;  /*index into a string*/
  Int4 numsegs; /*number of segments in an alignment*/
  Int4 startQ,startS; /*Indices into array of starting positions*/


  numalign = countSeqAligns(listOfSeqAligns, &numseq, FALSE, 0.0);
  search->posConverged = TRUE;
  curSeqAlign = listOfSeqAligns;
  if (thisPassNum > 1) {
    posSearch->posRepeatSequences = (Int2Ptr) MemNew(numalign * sizeof(Int2));
    prevSeqAlign = NULL;
    for(alignIndex = 0; alignIndex < numalign; alignIndex++) {
      posSearch->posRepeatSequences[alignIndex] = (Int2) 0;
      curSegs = (DenseSegPtr) curSeqAlign->segs;
      thisId = curSegs->ids->next; 
      if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(thisId, prevId)))) {
	startAlign = curSeqAlign;
	thisEvalue = minEvalueForSequence(curSeqAlign, startAlign);
	if (thisEvalue < search->pbp->ethresh) {    
	  /*Extract the ordinal number from the SeqAlign*/
	  curSegs = (DenseSegPtr) curSeqAlign->segs;
	  thisId = curSegs->ids->next;  /*id of target sequence is second*/
	  /*Get ordinal ids of sequences in result*/
	  ordinalNumber = SeqId2OrdinalId(search->rdfp, thisId);
	  found = FALSE;
	  for(oldSeqIndex = 0; oldSeqIndex < posSearch->posResultsCounter; oldSeqIndex++)
	    if(ordinalNumber ==  posSearch->posResultSequences[oldSeqIndex]) {
	      posSearch->posRepeatSequences[alignIndex] = SEQ_ALIGN_MARK_REPEAT;
	      found = TRUE;
	      break;
	    }      
	  if (!found) 
	    search->posConverged = FALSE;
	}
      }
      else  /*both alignments come from the same sequence*/
	posSearch->posRepeatSequences[alignIndex] = posSearch->posRepeatSequences[alignIndex - 1];
      prevSeqAlign = curSeqAlign;
      prevId = thisId;
      curSeqAlign = curSeqAlign->next;
    }
    MemFree(posSearch->posResultSequences);
  }
  else {
    q = search->context[0].query->sequence;
    queryLength = search->context[0].query->length;
    prevSeqAlign = NULL;
    while (curSeqAlign != NULL) {
      curSegs = (DenseSegPtr) curSeqAlign->segs;
      s = GetSequenceWithDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
      numsegs = curSegs->numseg;
      thisId = curSegs->ids->next; 
      if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(thisId, prevId)))) {
	startAlign = curSeqAlign;
	thisEvalue = minEvalueForSequence(curSeqAlign, startAlign);
	if (thisEvalue < search->pbp->ethresh) {    
	  if (numsegs > 1) {
	    search->posConverged = FALSE;
	    return;
	  }
	  startQ = 0;
	  startS = 1;
	  queryOffset = curSegs->starts[startQ];
	  if (curSegs->starts[startS] != GAP_HERE)
	    subjectOffset = curSegs->starts[startS] - retrievalOffset;
	  else
	    subjectOffset = GAP_HERE;
	  matchLength = curSegs->lens[0];
	  if ((queryOffset != 0) || (subjectOffset != 0) ||
	      (matchLength != queryLength) || (matchLength != subjectLength)) {
	    search->posConverged = FALSE;
	    return;
	  }
	  for (c = 0, qplace = queryOffset, splace = subjectOffset;
	       c < matchLength; c++, qplace++, splace++)
	    if (s[splace] != q[qplace]) {
	      search->posConverged = FALSE;
	      return;
	    }        
	}
      }
      prevSeqAlign = curSeqAlign;
      prevId = thisId;
      curSeqAlign = curSeqAlign->next;
    }
  }
}


/*Eliminate the matches from sequence second starting at position
matchStart and extending for intervalLength characters */
void LIBCALL posCancel(posSearchItems *posSearch, compactSearchItems * compactSearch, Int4 first, Int4 second, Int4 matchStart, Int4 intervalLength)
{
  Int4 c, i;
  Boolean stillNeeded;

  for(c = matchStart, i = 0; i < intervalLength; i++, c++) {
    posSearch->posDescMatrix[second][c].used = FALSE;
    posSearch->posDescMatrix[second][c].letter = 0;
  }
  stillNeeded = FALSE;
  for (c = 0; c < compactSearch->qlength; c++)
    if (posSearch->posDescMatrix[second][c].used) {
      stillNeeded = TRUE;
      break;
    }
   if (!stillNeeded)
      posSearch->posUseSequences[second] = FALSE;
}

/*Eliminate sequences that are identical to the query and partial alignments
  that are identical in two matching sequences
  Modified by Natsuhiko Futamura to change order in which
  pairs of sequences are compared*/
void LIBCALL posPurgeMatches(posSearchItems *posSearch, compactSearchItems * compactSearch, ValNodePtr * error_return)
{
  Int4 i, j; /*index over sequences*/
  Int4 k; /*difference between pair of sequence indices*/
  Boolean matchesQuery; /*Is a matching sequence identical to the query?*/
  Int4 c; /*index over demographics of matching sequence*/
  Int4 state; /*state of checking for a match*/
  Int4 intervalLength, matchStart; /*Length and start of a matching region*/
  Int4 matchNumber; /*number of characters matching*/
  Int4 Xcount; /*number of X's in interval*/
  Int4 numSequencesInUse; /*number of sequences left in use*/

  posSearch->posUseSequences =  (Boolean *) MemNew((posSearch->posNumSequences + 1) * sizeof(Boolean));
   if (NULL == posSearch->posUseSequences)
     exit(EXIT_FAILURE);
  for(i = 0; i <= posSearch->posNumSequences; i++)
    posSearch->posUseSequences[i] = TRUE;
  for(i = 1; i <= posSearch->posNumSequences; i++) {
    matchesQuery = TRUE;
    for (c = 0; c < compactSearch->qlength; c++) {
      if ((!posSearch->posDescMatrix[i][c].used) ||
          (posSearch->posDescMatrix[i][c].letter !=
           posSearch->posDescMatrix[0][c].letter)) {
        matchesQuery = FALSE;
        break;
      }
    }
    if (matchesQuery) {
      posSearch->posUseSequences[i] = FALSE;
    }
  }
  for(j = 1; j <= posSearch->posNumSequences; j++) {
    if (!posSearch->posUseSequences[j])
      continue;
    state = POS_COUNTING;
    c = 0;
    matchStart = 0;
    intervalLength = 0;
    Xcount = 0;
    matchNumber = 0;
    while (c < compactSearch->qlength) {
      if (posSearch->posDescMatrix[j][c].used) {
	if ((posSearch->posDescMatrix[0][c].letter != Xchar) &&
	    (posSearch->posDescMatrix[j][c].letter != Xchar)) { 
	  if (state == POS_RESTING) {
	    matchStart = c;
	    intervalLength = 1;
	    state = POS_COUNTING;
	    matchNumber = 0;
	  }
	  else 
	    intervalLength++;
	  if (posSearch->posDescMatrix[j][c].used &&
	      (posSearch->posDescMatrix[0][c].letter == posSearch->posDescMatrix[j][c].letter))
	    matchNumber++;
	}
	else {
	  if (POS_COUNTING == state)
	    Xcount++;
	}
      }
      else {
	if (state == POS_COUNTING) {
	  if ((intervalLength > 0) && (matchNumber == intervalLength))
	    posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength+Xcount);
	  state = POS_RESTING;
	  Xcount = 0;
	}
      }
      c++;
    }
    if (state == POS_COUNTING) /*at end of sequence i*/
      if ((intervalLength > 0) && (matchNumber == intervalLength))
	posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength+Xcount);
  }
  
  for (k=1; k <= posSearch->posNumSequences -1; k++){
    for (i = 1; (i+k) <= posSearch->posNumSequences; i++) {
      if (!posSearch->posUseSequences[i])
	continue;
      j = i+k;
      if (!posSearch->posUseSequences[j])
	continue;

      state = POS_COUNTING;
      c = 0;
      matchStart = 0;
      intervalLength = 0;
      Xcount = 0;
      matchNumber = 0;
      while (c < compactSearch->qlength) {
	if (posSearch->posDescMatrix[i][c].used ||
	    posSearch->posDescMatrix[j][c].used) {
	  if ((posSearch->posDescMatrix[i][c].letter != Xchar) &&
	      (posSearch->posDescMatrix[j][c].letter != Xchar)) { 
	    if (state == POS_RESTING) {
	      matchStart = c;
	      intervalLength = 1;
	      state = POS_COUNTING;
	      matchNumber = 0;
	    }
	    else 
	      intervalLength++;
	    if (posSearch->posDescMatrix[i][c].used &&
		posSearch->posDescMatrix[j][c].used &&
		(posSearch->posDescMatrix[i][c].letter == posSearch->posDescMatrix[j][c].letter))
	      matchNumber++;
	  }
	  else {
	    if (POS_COUNTING == state)
	      Xcount++;
	  }
	}
	else {
	  if (state == POS_COUNTING) {
	    if ((intervalLength > 0) && ((((Nlm_FloatHi) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	      posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength+Xcount);
	    state = POS_RESTING;
	    Xcount = 0;
	  }
	}
	c++;
      }
      if (state == POS_COUNTING) /*at end of sequence i*/
	if ((intervalLength > 0) && ((((Nlm_FloatHi) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	  posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength+Xcount);
    }
  }
  if (error_return) {
    numSequencesInUse = 0;
    for(i = 0; i <= posSearch->posNumSequences; i++)
      if(posSearch->posUseSequences[i])
	numSequencesInUse++;
    if (numSequencesInUse < 2)
      BlastConstructErrorMessage("posPurgeMatches", "Due to purging near identical sequences, only the query is used to construct the position-specific score matrix\n", 1, error_return);
  }
}

static void countNumSeq(posSearchItems *posSearch,
                  compactSearchItems * compactSearch,
                  SeqAlignPtr listOfSeqAligns, Int4 *prevNumSeq)
{	
   Uint1Ptr s; /* pointer into a matching string */
   Int4  subjectLength;  /* length of subject */
   Int4  retrievalOffset;   /* retrieval offset */
   SeqAlignPtr curSeqAlign, prevSeqAlign; /* pointers into listOfSeqAligns */
   DenseSegPtr curSegs, prevSegs;  /* used to extract alignments from curSeqAlign */
   SeqIdPtr curId, prevId;  /* Used to compare sequences that come from different SeqAligns */
   Nlm_FloatHi thisEvalue;  /* evalue of current partial alignment */
   Int4 newNumSeq; /* numseq computed in another way */
   Boolean is_new_id = FALSE;	

   newNumSeq = 0;
   /*use only those sequences below e-value threshold*/
   curSeqAlign = listOfSeqAligns;
   prevSeqAlign = NULL;
   for (curSeqAlign = listOfSeqAligns; curSeqAlign != NULL;
                            curSeqAlign = curSeqAlign->next) {
      is_new_id = FALSE;
      thisEvalue = getEvalueFromSeqAlign(curSeqAlign);
      curSegs = (DenseSegPtr) curSeqAlign->segs;
      if (NULL != prevSeqAlign) {
         prevSegs = (DenseSegPtr) prevSeqAlign->segs;
         if(curSegs->ids == NULL)
            break;
         curId = curSegs->ids->next;
         prevId = prevSegs->ids->next;

         if (!(SeqIdMatch(curId, prevId)))
            is_new_id = TRUE;
      }
      if (!(compactSearch->use_best_align && is_new_id)) {
         if (thisEvalue >= compactSearch->ethresh)
            continue;
      }
      if (is_new_id == TRUE)
         newNumSeq++;
      s = GetSequenceWithDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
      SeqMgrFreeCache();
      if (s == NULL)
         continue;
      s = MemFree(s);
      prevSeqAlign = curSeqAlign;
   }
   newNumSeq++;
   /* numseq gets the highest number computed by both methods */
   if (newNumSeq > *prevNumSeq)
      *prevNumSeq = newNumSeq;
}

/*Compute general information about the sequences that matched on the
  i-th pass such as how many matched at each query position and what letter
  matched*/
void LIBCALL posDemographics(posSearchItems *posSearch, 
                             compactSearchItems * compactSearch, 
                             SeqAlignPtr listOfSeqAligns)
{
   Uint1Ptr q; /*pointers into query */
   Uint1Ptr s; /*pointer into a matching string */
   Int4 length, subjectLength;  /*length of query and subject*/
   Int4 c; /*index into a string*/
   Int4 numseq, numSeqAligns;  /*number of matching sequences and SeqAligns*/
   Int4 seqIndex;  /*index for the array of matching sequences*/
   Int4 matchLength; /*length of a match*/
   Int4  queryOffset, subjectOffset, retrievalOffset;  /*offsets needed to make a match align*/
   Int4 qplace, splace; /*index into query string and matching string*/
   SeqAlignPtr curSeqAlign, prevSeqAlign; /*pointers into listOfSeqAligns*/
   DenseSegPtr curSegs, prevSegs;  /*used to extract alignments from curSeqAlign*/
   SeqIdPtr curId, prevId;  /*Used to compare sequences that come from different SeqAligns*/
   Int4 startQ, startS; /*Indices into array of starting positions*/
   Int4 numsegs; /*Number of pieces in the gapped alignment*/
   Int4 segIndex; /*Index for which piece we are at*/
   Nlm_FloatHi thisEvalue;  /*evalue of current partial alignment*/
   Boolean is_new_id = FALSE;

   q = compactSearch->query;
   length = compactSearch->qlength;
   for(c = 0; c < length; c++) {
     posSearch->posDescMatrix[0][c].letter = (Int1) q[c];
     posSearch->posDescMatrix[0][c].used = TRUE;
     posSearch->posDescMatrix[0][c].leftExtent = 0;
     posSearch->posDescMatrix[0][c].rightExtent = length;
     posSearch->posDescMatrix[0][c].e_value = compactSearch->ethresh/2;
     posSearch->posC[c][q[c]]++;
     posSearch->posCount[c]++;
   }

   numSeqAligns = countSeqAligns(listOfSeqAligns, &numseq, 
                                 !compactSearch->use_best_align, 
                                 compactSearch->ethresh);
   
   posSearch->posNumSequences = numseq;
   /*use only those sequences below e-value threshold*/
   seqIndex = 0;
   curSeqAlign = listOfSeqAligns;
   prevSeqAlign = NULL;
   for(curSeqAlign = listOfSeqAligns; curSeqAlign != NULL; 
       curSeqAlign = curSeqAlign->next) {
       is_new_id = FALSE;

       thisEvalue = getEvalueFromSeqAlign(curSeqAlign);

       curSegs = (DenseSegPtr) curSeqAlign->segs;
       if (NULL != prevSeqAlign) {
           prevSegs = (DenseSegPtr) prevSeqAlign->segs;
           
           if(curSegs->ids == NULL) 
               break;
           
           curId = curSegs->ids->next; 
           prevId = prevSegs->ids->next;
           if (!(SeqIdMatch(curId, prevId)))
               is_new_id = TRUE;
       }
       
       if(!(compactSearch->use_best_align && is_new_id)) { 
           if (thisEvalue >= compactSearch->ethresh)
               continue;
       }
       
       if(is_new_id == TRUE)
           seqIndex++;
       
       s = GetSequenceWithDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
       if (s == NULL) {
           /* Kludge: set all of this sequence's residues to those of the query
            * so that it can be purged in posPurgeMatches */
           for (c = 0; c < length; c++) {
               posSearch->posDescMatrix[seqIndex+1][c].letter = (Int1) q[c];
               posSearch->posDescMatrix[seqIndex+1][c].used = TRUE;
               posSearch->posDescMatrix[seqIndex+1][c].e_value =
                   compactSearch->ethresh/2;
           }
           continue;
       }
       startQ = 0;
       startS = 1;
       numsegs = curSegs->numseg;
       for(segIndex = 0; segIndex < numsegs; segIndex++) {
           queryOffset = curSegs->starts[startQ];
           if (curSegs->starts[startS] != GAP_HERE) /*XX*/
               subjectOffset = curSegs->starts[startS] - retrievalOffset;
           else
               subjectOffset = GAP_HERE;
           matchLength = curSegs->lens[segIndex];
           if ((GAP_HERE ) == queryOffset) {
               ; /*do nothing, gap in query*/
           }
           else
	     if ((GAP_HERE) == subjectOffset) { /*XX*/
                   for(c = 0, qplace = queryOffset;
                       c < matchLength; c++, qplace++) {
                     /*Keep the following test if spreading out gap costs, 
                       so that in that case a lower E-value non-gap trumps
                       a higher E-value gap; if not spreading out gap costs
                       then comment out the test, so that a higher E-value
                       gap trumps a lower E-value letter*/
		     if (!posSearch->posDescMatrix[seqIndex+1][qplace].used)
		       {
			 posSearch->posDescMatrix[seqIndex + 1][qplace].used = TRUE;
			 posSearch->posDescMatrix[seqIndex + 1][qplace].letter = GAP_CHAR;
			 posSearch->posDescMatrix[seqIndex + 1][qplace].e_value = 1.0;
		       }
                   }
               }
               else {  /*no gap*/
                   for(c = 0, qplace = queryOffset, splace = subjectOffset;
                       c < matchLength; c++, qplace++, splace++) {
		     if (!posSearch->posDescMatrix[seqIndex+1][qplace].used)
		       {
			 posSearch->posDescMatrix[seqIndex+1][qplace].letter = (Int1) s[splace]; 
			 posSearch->posDescMatrix[seqIndex+1][qplace].used = TRUE; 
			 posSearch->posDescMatrix[seqIndex+1][qplace].e_value = 
			   thisEvalue;
		       }
                   }
               }
           startQ += 2;
           startS += 2;
       }
       prevSeqAlign = curSeqAlign;
       s = MemFree(s);
   } /*closes the for loop over seqAligns*/
}

void LIBCALL posComputeExtents(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Int4 seqIndex; /*index of sequence*/
   Int4 length; /*length of query*/
   Int4 qplace, qplace2; /*place in query*/
   Int4 numseq; /*number of sequences including query*/
   Uint1Ptr q; /*pointers into query */

   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   q = compactSearch->query;
   for(seqIndex = 0; seqIndex < numseq; seqIndex++) {
     if (!posSearch->posUseSequences[seqIndex+1])
       continue; /*XX*/
     if ((posSearch->posDescMatrix[seqIndex+1][0].used) 
	 && (posSearch->posDescMatrix[seqIndex+1][0].letter != GAP_CHAR))
       posSearch->posDescMatrix[seqIndex+1][0].leftExtent = 0;
     for(qplace = 1; qplace < length; qplace++)
       if(posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 if(posSearch->posDescMatrix[seqIndex+1][qplace-1].used)
	   posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent =
	     posSearch->posDescMatrix[seqIndex+1][qplace -1].leftExtent;
	 else
	   posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent = qplace;
       } 
     if ((posSearch->posDescMatrix[seqIndex+1][length-1].used)
	 && (posSearch->posDescMatrix[seqIndex+1][length-1].letter != GAP_CHAR))
       posSearch->posDescMatrix[seqIndex+1][length-1].rightExtent = length -1;
     for(qplace = length -2; qplace >= 0; qplace--)
       if(posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 if(posSearch->posDescMatrix[seqIndex+1][qplace+1].used)
	   posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent =
	     posSearch->posDescMatrix[seqIndex+1][qplace + 1].rightExtent;
	 else
	   posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent = qplace;
       }
     for(qplace = 0; qplace < length; qplace++) 
       if (posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 /* comment next if out to spread gap costs*/
	 /* if (posSearch->posDescMatrix[seqIndex+1][qplace].letter != GAP_CHAR) { */
	   posSearch->posExtents[qplace].leftExtent = MAX(posSearch->posExtents[qplace].leftExtent,
							posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent);
	   posSearch->posExtents[qplace].rightExtent = MIN(posSearch->posExtents[qplace].rightExtent,
							   posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent);
	 
	 }
     /*}*/ /*comment this out if we want to spread gap costs out */

     for(qplace = 0; qplace < length; qplace++) 
       /*used to check qplace for GAP_CHAR here*/ /*XX*/
       if (posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 posSearch->posC[qplace][posSearch->posDescMatrix[seqIndex+1][qplace].letter]++;
	 posSearch->posCount[qplace]++; /*Add to number of matches in this query position*/
       }
   }
   for(qplace = 0; qplace < length; qplace++)
     posSearch->posIntervalSizes[qplace] = posSearch->posExtents[qplace].rightExtent - 
       posSearch->posExtents[qplace].leftExtent + 1;
   for(qplace =0; qplace < length; qplace++) {
     if(Xchar == q[qplace]) {
       posSearch->posIntervalSizes[qplace] = 0;
       for(qplace2 = 0; qplace2 <qplace; qplace2++) {
	 if((Xchar != q[qplace2]) && (posSearch->posExtents[qplace2].rightExtent >= qplace))
	   posSearch->posIntervalSizes[qplace2]--;
       }
       for(qplace2 = length-1; qplace2 > qplace; qplace2--) {
	 if((Xchar != q[qplace2]) && (posSearch->posExtents[qplace2].leftExtent <= qplace))
	   posSearch->posIntervalSizes[qplace2]--;
       }
     }
   }
}
 
/*Compute weight of each sequence and letter in each position*/
void LIBCALL posComputeSequenceWeights(posSearchItems *posSearch, compactSearchItems * compactSearch, Nlm_FloatHi weightExponent)
{
   Int4 length; /*length of query*/
   Int4 numseq, seqIndex; /*number of matches, index for them*/
   Int4  i; /*index over a multi-alignment block*/
   Int4 qplace; /*index into query*/
   Nlm_FloatHi Sigma; /*Number of different characters occurring in matches within
                   a multi-alignment block, excluding identical columns*/
   Nlm_FloatHi intervalSigma; /*Same as Sigma but includes identical columns*/
   Int4 alphabetSize; /*number of characters in alphabet*/
   Int4 *participatingSequences; /*array of participating sequences at a position*/
   Int4 *oldParticipatingSequences; /*array of participating sequences at a position*/
   Int4 posLocalVariety;  /*number of different characters at a position*/
   Int4 posLocalStandardLet; /*posLocalVariety, not counting X or gap*/
   Int4 *posLocalC; /*counts of how many of each letter in this column*/
   Int4 c;
   Int4 thisSeq;
   Int4 numParticipating; /*number of sequences in this alignment block*/
   Int4 oldNumParticipating; /*number of sequences in this alignment block*/
   Boolean newSequenceSet; 
   Int4 p; /*index on sequences*/
   Nlm_FloatHi weightSum; /*Sum of intermediate sequence weights in a column 
                            used to normalize the weights, so they sum to 1*/

   alphabetSize = compactSearch->alphabetSize;
   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   participatingSequences = (Int4 *) MemNew((numseq+1) * sizeof(Int4));
   if (NULL == participatingSequences)
     exit(EXIT_FAILURE);
   oldParticipatingSequences = (Int4 *) MemNew((numseq+1) * sizeof(Int4));
   if (NULL == oldParticipatingSequences)
     exit(EXIT_FAILURE);
   posLocalC = (Int4 *) MemNew(alphabetSize * sizeof(Int4));
   if (NULL == posLocalC)
     exit(EXIT_FAILURE);
   for (qplace = 0; qplace < length; qplace++) {
     posSearch->posSigma[qplace] = 0.0;
   }
   numParticipating = 0;
   for(qplace = 0; qplace < length; qplace++) {
     posSearch->posGaplessColumnWeights[qplace] = 0.0;
     if ((posSearch->posCount[qplace] > 1) && (posSearch->posIntervalSizes[qplace] > 0)) {
       oldNumParticipating = numParticipating;
       for(p =0; p < numParticipating; p++)
         oldParticipatingSequences[p] = participatingSequences[p];
       numParticipating = 0;
       for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
         if (!posSearch->posUseSequences[seqIndex])
           continue; 
	 /* if ((posSearch->posDescMatrix[seqIndex][qplace].used) &&
	     (posSearch->posDescMatrix[seqIndex][qplace].letter != GAP_CHAR)) {
	 */
	 /*change to this if we want to spread gap costs*/
	 if (posSearch->posDescMatrix[seqIndex][qplace].used) {
	     participatingSequences[numParticipating] = seqIndex; 
	   numParticipating++;
	 }
       }
       newSequenceSet = TRUE;
       if (numParticipating == oldNumParticipating) {
         for(p = 0; p < numParticipating; p++)
           if (oldParticipatingSequences[p] != participatingSequences[p])
             break;
         if (p == numParticipating)
           newSequenceSet = FALSE;
       }
         
       if (newSequenceSet) {
	 Sigma = 0;
	 intervalSigma = 0;
	 for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
	   if (!posSearch->posUseSequences[seqIndex])
	     continue;
	   posSearch->posRowSigma[seqIndex] = 0.0;
	   posSearch->posA[seqIndex] = 0.0;
	 }
	 for (i = posSearch->posExtents[qplace].leftExtent;
	      i <= posSearch->posExtents[qplace].rightExtent; i++) {
	   posLocalVariety = 0;
	   posLocalStandardLet = 0;
	   for(c = 0; c < alphabetSize; c++)
	     posLocalC[c] = 0;
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     /*used to check for GAP here*/ /*XX*/
	     if (0 == posLocalC[posSearch->posDescMatrix[thisSeq][i].letter]) {
	       /*letter (not a gap) not seen before in this query pos.*/
	       posLocalVariety++;  
	       if ((GAP_CHAR != posSearch->posDescMatrix[thisSeq][i].letter)  &&
		 (Xchar != posSearch->posDescMatrix[thisSeq][i].letter))
		 posLocalStandardLet++;
	     }
	     posLocalC[posSearch->posDescMatrix[thisSeq][i].letter]++;
	   }
	   intervalSigma += posLocalVariety;
	   posLocalStandardLet = MIN(posLocalStandardLet,EFFECTIVE_ALPHABET);
	   posSearch->posDistinctDistrib[qplace][posLocalStandardLet]++;
	   if (posLocalVariety > 1) {
	     Sigma += posLocalVariety;
	   }
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     /*used to check for gap here*/
	     posSearch->posRowSigma[thisSeq] += 
		( 1.0 / 
  (((Nlm_FloatHi) posLocalC[posSearch->posDescMatrix[thisSeq][i].letter])
    * posLocalVariety));
	   }
	 }
       }
       else {
	 for (i = 0; i <= EFFECTIVE_ALPHABET; i++) {
	   posSearch->posDistinctDistrib[qplace][i] = posSearch->posDistinctDistrib[qplace - 1][i];	   
	 }
       }	 
       if (Sigma > 0) {
	 weightSum = 0;
	 for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = posSearch->posRowSigma[thisSeq]/
	    (posSearch->posExtents[qplace].rightExtent -
	       posSearch->posExtents[qplace].leftExtent +1);
	   /*spread gap weight here*/
           posSearch->posA[thisSeq] = pow(posSearch->posA[thisSeq],
                                          weightExponent);
	   weightSum += posSearch->posA[thisSeq];
	 }
	 for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = posSearch->posA[thisSeq]/weightSum;
	 }
       }
       else {
         for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = ((Nlm_FloatHi) 1 / (Nlm_FloatHi) numParticipating);
         }
       }
       posSearch->posSigma[qplace] = intervalSigma;
       for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 posSearch->posMatchWeights[qplace][posSearch->posDescMatrix[thisSeq][qplace].letter] += posSearch->posA[thisSeq];
         if(posSearch->posDescMatrix[thisSeq][qplace].letter)
           posSearch->posGaplessColumnWeights[qplace] += posSearch->posA[thisSeq]; 
       }
       posSearch->posNumParticipating[qplace] = numParticipating;
     }
   }
   MemFree(participatingSequences);
   MemFree(oldParticipatingSequences);
   MemFree(posLocalC);
}

static Nlm_FloatHi countsFunction(Nlm_FloatHi Sigma, Int4 intervalLength)
{
  return(Sigma / intervalLength - 1);
}

#define MAX_IND_OBSERVATIONS  400


/*initialize the expected number of observations
  use background probabilities for this matrix
  Calculate exp. # of distinct aa's as a function of independent trials   
*/ 
static void initializeExpNumObservations(double *expno, 
				    double *backgroundProbabilities)

{
int     j,k ; /*loop indices*/
double  weighted_sum; /*20 - this is how many distinct
			 amino acids are expected*/

   expno[0] = 0;
   for (j=1;j<MAX_IND_OBSERVATIONS;++j) {
     weighted_sum = 0;
     for (k=0;k<EFFECTIVE_ALPHABET;++k) 
       weighted_sum += exp(j*log(1.0-backgroundProbabilities[k]));
     expno[j] = EFFECTIVE_ALPHABET-weighted_sum;
   }
}


/*A method to estimate the effetive number of observations
  in the interval for the specified columnNumber */

static Nlm_FloatHi effectiveObservations(posSearchItems *posSearch, 
					 Int4 columnNumber, Int4 queryLength,
					 double *expno)
{
int     i,k; /*loop indices*/
double  indep; /*number of independent observations to return*/
int halfNumColumns; /*half the number of columns in the interval, rounded
                      down*/
int totalDistinctCounts; /*total number of distinct letters in columns
		     used*/
double aveDistinctAA; /*average number of distinct letters in columns used*/
int columnsAccountedFor; /*how many of the columns had their
                            distinct count totaled so far*/

  
 if (posSearch->posExtents[columnNumber].leftExtent < 0)
   return(0);
 if (posSearch->posExtents[columnNumber].rightExtent >= queryLength)
   return(0);
 
/*  Calculate the average number of distinct amino acids in the half of the
    columns within the block in question with the most distinct amino acids;
    +2 in the parentheses is for rounding up.*/

 halfNumColumns = MAX(1,(posSearch->posExtents[columnNumber].rightExtent -
			 posSearch->posExtents[columnNumber].leftExtent+2)/2);
 k = EFFECTIVE_ALPHABET;
 columnsAccountedFor = 0;
 totalDistinctCounts = 0;
 while (columnsAccountedFor < halfNumColumns) {
   totalDistinctCounts += (posSearch->posDistinctDistrib[columnNumber][k] *k);
   columnsAccountedFor += posSearch->posDistinctDistrib[columnNumber][k];
   if (columnsAccountedFor > halfNumColumns) {
     totalDistinctCounts -=
       ((columnsAccountedFor - halfNumColumns) * k);
     columnsAccountedFor = halfNumColumns;
   }
   k--;
 }
 aveDistinctAA = ((double) totalDistinctCounts)/
   ((double) columnsAccountedFor);

/*    Then use the following code to calculate the number of
        independent observations corresponding to
        aveDistinctAA.
*/

 for (i=1;i<MAX_IND_OBSERVATIONS && expno[i]<=aveDistinctAA;++i);
 indep = (i==MAX_IND_OBSERVATIONS) ? i : 
   i-(expno[i]-aveDistinctAA)/(expno[i]-expno[i-1]);
 indep = MIN(indep, posSearch->posNumParticipating[columnNumber]);	
 indep = MAX(0,indep - 1);
 return(indep);
}

static Nlm_FloatHi posit_rounddown(Nlm_FloatHi value)
{
  return (Nlm_FloatHi) Nlm_Nint(value);
}

/*check that weights add to 1 in each column */
void LIBCALL posCheckWeights(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Uint1Ptr q;  /*pointer to query*/
   Int4 length, alphabetSize; /*length of query and number of characters in alphabet*/
   Int4  a, c; /*loop indices*/
   Nlm_FloatHi runningSum; /*partial total for a column*/


   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;

   q = compactSearch->query;
   for(c = 0; c < length; c++) {
     if ((posSearch->posCount[c] > 1) && (q[c] != Xchar)) {
       runningSum = 0;
       /*       if (posSearch->posMatchWeights[c][0] > 0.0)
		printf("Stop here %d ", c); */
       for(a = 0; a < alphabetSize; a++) 
           runningSum += posSearch->posMatchWeights[c][a];
       if((runningSum < 0.99) || (runningSum > 1.01))
         ErrPostEx(SEV_ERROR, 0, 0, "\nERROR IN WEIGHTS, column %d, value %lf\n",c, runningSum);
       /* spread out gap weight here*/
       for(a = 1; a < alphabetSize; a++) 
	 if (compactSearch->standardProb[a] > posEpsilon)
	   posSearch->posMatchWeights[c][a] = posSearch->posMatchWeights[c][a] +
           (posSearch->posMatchWeights[c][0] * compactSearch->standardProb[a]);
       posSearch->posMatchWeights[c][0] = 0.0;
       runningSum = 0;
       for(a = 0; a < alphabetSize; a++) 
           runningSum += posSearch->posMatchWeights[c][a];
       if((runningSum < 0.99) || (runningSum > 1.01))
         ErrPostEx(SEV_ERROR, 0, 0, "\nERROR IN WEIGHTS, column %d, value %lf\n",c, runningSum);
     }
   }
}

/*Fill in information content per position from pseudo-count frequencies*/
static void  posFreqsToInformation(posSearchItems * posSearch, compactSearchItems * compactSearch)
{
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, alphabetSize; /*loop index and size of alphabet*/
   Nlm_FloatHi  qOverPEstimate; /*intermediate term*/
   Nlm_FloatHi  infoSum; /*information content sum for this position*/
  
   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;
   for (c = 0; c < length; c++) {
     infoSum = 0;
     for(a = 0; a < alphabetSize; a++) {
       if (compactSearch->standardProb[a] > posEpsilon) {
         qOverPEstimate = posSearch->posFreqs[c][a] / compactSearch->standardProb[a];
         if (qOverPEstimate > posEpsilon)
	   infoSum += posSearch->posFreqs[c][a] * log(qOverPEstimate)/
                    NCBIMATH_LN2;
       }
     }
     posSearch->posInformation[c] = infoSum;
   }
}

/*Convert pseudo-count frequencies to a score matrix, where standard
matrix is represented by its frequencies */
void LIBCALL posFreqsToMatrix(posSearchItems *posSearch, 
                              compactSearchItems *compactSearch)
{
   Uint1Ptr q;  /*pointer to the query*/
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, alphabetSize; /*loop index and size of alphabet*/
   Nlm_FloatHi lambda; /*Karlin-Altschul parameter*/
   Nlm_FloatHi  qOverPEstimate, value; /*intermediate terms*/
   Boolean allZeros; /*are all frequencies in a column 0?*/
   Nlm_FloatHi intermediateValue; /*intermediate value*/

   q = compactSearch->query;
   length = compactSearch->qlength;

   alphabetSize = compactSearch->alphabetSize;
   lambda = compactSearch->lambda_ideal;


   for(c = 0; c < length; c++) {
     allZeros = TRUE;
     for(a = 0; a < alphabetSize; a++) {
       /*Division compensates for multiplication in posComputePsedoFreqs*/
       if (compactSearch->standardProb[a] > posEpsilon)
	 qOverPEstimate = posSearch->posFreqs[c][a]/compactSearch->standardProb[a];
       else
        qOverPEstimate = 0.0;
       if (qOverPEstimate != 0.0)
         allZeros = FALSE;
       if (0.0 == qOverPEstimate || (compactSearch->standardProb[a] < posEpsilon))
	 posSearch->posPrivateMatrix[c][a] = BLAST_SCORE_MIN;
       else {
	 value = log(qOverPEstimate)/lambda;
	 posSearch->posPrivateMatrix[c][a] = (BLAST_Score) posit_rounddown(POSIT_SCALE_FACTOR*value);

       }
       if (((Xchar == a) || (StarChar == a)) && (compactSearch->matrix[q[c]][Xchar] != BLAST_SCORE_MIN))
	 posSearch->posPrivateMatrix[c][a] = 
	   compactSearch->matrix[q[c]][a] * POSIT_SCALE_FACTOR;
     }
     if (allZeros) {
       if ( !posSearch->stdFreqRatios ) {
         ErrPostEx(SEV_FATAL, 1, 0, "Frequency ratios for %s scoring matrix are not available\n", compactSearch->standardMatrixName);
         return;
       }

       for(a = 0; a < alphabetSize; a++) {
         posSearch->posMatrix[c][a] = compactSearch->matrix[q[c]][a];
         if (posSearch->stdFreqRatios->data[q[c]][a] == 0.0) {
           posSearch->posPrivateMatrix[c][a] = BLAST_SCORE_MIN;
         } else {
           intermediateValue = POSIT_SCALE_FACTOR *
             posSearch->stdFreqRatios->bit_scale_factor *
             log(posSearch->stdFreqRatios->data[q[c]][a])/NCBIMATH_LN2;
           posSearch->posPrivateMatrix[c][a] = Nlm_Nint(intermediateValue);
         }
       }
     }
   }
   for(a = 0; a < alphabetSize; a++) {
     posSearch->posPrivateMatrix[length][a] = posSearch->posPrivateMatrix[length][a] = BLAST_SCORE_MIN;
   }
}

/*copy position specific frequency matrix of diminesions qlength * alphabetSize*/
void LIBCALL copyPosFreqs(Nlm_FloatHi **posFreqsFrom, Nlm_FloatHi **posFreqsTo, Int4 qlength, Int4 alphabetSize)
{
  Int4 c, i; /*loop indices*/

  for (i = 0; i < qlength; i++)
    for (c = 0; c < alphabetSize; c++)
      posFreqsTo[i][c] = posFreqsFrom[i][c];
}

Nlm_FloatHi ** LIBCALL allocatePosFreqs(Int4 length, Int4 alphabetSize)
{
  Int4 c, i; /*loop indices*/
  Nlm_FloatHi ** returnArray;

  returnArray = (Nlm_FloatHi **) MemNew((length + 1) * sizeof(Nlm_FloatHi *));
  if (NULL == returnArray)
    exit(EXIT_FAILURE);
  for(i = 0; i <= length; i++) {
    returnArray[i] = (Nlm_FloatHi *) MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    if (NULL == returnArray[i])
      exit(EXIT_FAILURE);   
    for(c = 0; c < alphabetSize; c++)
      returnArray[i][c] = 0.0;
  }
  return(returnArray); 
}

/*The following constants are used in 
  posComputePseudoFreqs and columnSpecificPseudocounts */
#define PSEUDO_MULTIPLIER 500
#define PSEUDO_SMALL_INITIAL 5.5 /*small number of pseudocounts to
                                   avoid 0 probabilities in entropy-based
                                   method*/
#define PSEUDO_NUMERATOR 0.0457  /*numerator of entropy-based method*/
#define PSEUDO_EXPONENT 0.8  /*exponent of denominator*/
#define PSEUDO_MAX 1000000 /*effective infinity*/
#define ZERO_OBS_PSEUDO 30 /*arbitrary constant to use for columns with 
                             zero observations in actual data*/

static void  fillColumnProbabilities(double *probabilities, 
				     posSearchItems *posSearch, 
				     Int4 columnNumber)
{
   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/
   Int4 c; /*loop index*/

   charOrder[0] =  1;  /*A*/
   charOrder[1] =  16; /*R*/
   charOrder[2] =  13; /*N*/  
   charOrder[3] =  4;  /*D*/ 
   charOrder[4] =  3;  /*C*/
   charOrder[5] =  15; /*Q*/
   charOrder[6] =  5;  /*E*/ 
   charOrder[7] =  7;  /*G*/
   charOrder[8] =  8;  /*H*/
   charOrder[9] =  9;  /*I*/
   charOrder[10] = 11; /*L*/
   charOrder[11] = 10; /*K*/
   charOrder[12] = 12; /*M*/  
   charOrder[13] =  6; /*F*/
   charOrder[14] = 14; /*P*/
   charOrder[15] = 17; /*S*/
   charOrder[16] = 18; /*T*/
   charOrder[17] = 20; /*W*/
   charOrder[18] = 22; /*Y*/
   charOrder[19] = 19; /*V*/

   for(c = 0; c < EFFECTIVE_ALPHABET; c++) 
     probabilities[c] = posSearch->posMatchWeights[columnNumber][charOrder[c]];
}

/*adjust the probabilities by assigning observations weight
  to initialProbabilities and standardWeight to standardProbabilities*/
static void adjustColumnProbabilities(double *initialProbabilities,
				      double *probabilitiesToReturn,
				      double standardWeight, 
				      double *standardProbabilities, 
				      double observations)
{
  double intermediateSums[EFFECTIVE_ALPHABET]; /*weighted sums for each letter*/
  double overallSum; /*overall sum of weightedSums*/
  Int4 c; /*loop index*/

  overallSum = 0.0;
  for(c = 0; c < EFFECTIVE_ALPHABET; c++) {
    intermediateSums[c] =
      (initialProbabilities[c] * observations) +
      (standardProbabilities[c] * standardWeight);
    overallSum += intermediateSums[c];
  }
  for(c = 0; c < EFFECTIVE_ALPHABET; c++) 
    probabilitiesToReturn[c] = intermediateSums[c]/overallSum;
}

/*compute relative entropy of first distribution to second distribution*/

static double computeRelativeEntropy(double *newDistribution,
			      double *backgroundProbabilities)
{
   Int4 c; /*loop index*/
   double returnValue; /*value to return*/
   

   returnValue = 0;
   for(c = 0; c < EFFECTIVE_ALPHABET; c++) {
     if (newDistribution[c] > posEpsilon)
       returnValue += (newDistribution[c] * 
		       log (newDistribution[c]/backgroundProbabilities[c]));
   }
   if (returnValue < posEpsilon)
     returnValue = posEpsilon;
   return(returnValue);
}


static double columnSpecificPseudocounts(posSearchItems *posSearch, 
				  compactSearchItems *compactSearch, 
				  Int4 columnNumber, 
				  double *backgroundProbabilities,
				  double observations)
{
  double columnProbabilitiesInitial[EFFECTIVE_ALPHABET];
  double columnProbabilitiesAdjusted[EFFECTIVE_ALPHABET];
  double relativeEntropy; /*relative entropy of this column to background probs.*/
  double alpha; /*intermediate term*/
  double pseudoDenominator; /*intermediate term*/
  double returnValue;

  fillColumnProbabilities(&(columnProbabilitiesInitial[0]), posSearch, columnNumber);
  adjustColumnProbabilities(&(columnProbabilitiesInitial[0]), 
			    &(columnProbabilitiesAdjusted[0]),
			      compactSearch->standardProbWeight,
			      backgroundProbabilities, observations);
  relativeEntropy = computeRelativeEntropy(&(columnProbabilitiesAdjusted[0]),
					   backgroundProbabilities);
  pseudoDenominator = pow(relativeEntropy, compactSearch->HmethodDenominator);
  alpha = compactSearch->HmethodNumerator/pseudoDenominator;
  if (alpha < (1.0 - posEpsilon))
    returnValue = PSEUDO_MULTIPLIER * alpha/ (1- alpha);
  else
    returnValue = PSEUDO_MAX;
  /*extraOutputFile = fopen("Pseudocounts.txt","a");
  fprintf(extraOutputFile,"%s\t%5d\t%3.6lf\t%3.6lf\t%3.6lf\t%3.6lf",compactSearch->queryFileName,columnNumber+1,observations+1,relativeEntropy,alpha,returnValue);
  fprintf(extraOutputFile,"\n");
  fclose(extraOutputFile);*/

  return(returnValue);
}


Nlm_FloatHi ** LIBCALL posComputePseudoFreqs(posSearchItems *posSearch, compactSearchItems * compactSearch, Boolean Cpos)
{
   Uint1Ptr q;  /*pointer to the query*/
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, aSub, alphabetSize; /*loop indices and size of alphabet*/
   Nlm_FloatHi lambda; /*Karlin-Altschul parameter*/
   Nlm_FloatHi pseudo, numerator, denominator, qOverPEstimate; /*intermediate terms*/
   Nlm_FloatHi infoSum; /*sum used for information content*/
   Nlm_FloatHi **posFreqs; /*store frequencies*/
   double observations; /*estimated number of independent observations*/
   double  expno[MAX_IND_OBSERVATIONS+1]; /*table of expectations*/
   double pseudoWeight; /*multiplier for pseudocounts term*/
   double *backgroundProbabilities; /*background probabilities for matrix*/
   double columnCounts; /*column-specific pseudocounts*/

   q = compactSearch->query;
   length = compactSearch->qlength;

   alphabetSize = compactSearch->alphabetSize;
   lambda = compactSearch->lambda_ideal;
   posFreqs = allocatePosFreqs(length, alphabetSize);
   backgroundProbabilities = (double *) Blast_GetMatrixBackgroundFreq(compactSearch->standardMatrixName);
   if (!posSearch->stdFreqRatios) {
     posSearch->stdFreqRatios =
           PSIMatrixFrequencyRatiosNew(compactSearch->standardMatrixName);
   }

   initializeExpNumObservations(&(expno[0]),  backgroundProbabilities);
   compactSearch->standardProbWeight = PSEUDO_SMALL_INITIAL;
   compactSearch->HmethodDenominator = PSEUDO_EXPONENT;
   compactSearch->HmethodNumerator = PSEUDO_NUMERATOR;

   for(c = 0; c < length; c++) {
     if (Xchar != q[c]) {
       infoSum = 0;
       observations = effectiveObservations(posSearch,c,
                                            compactSearch->qlength,
                                            &(expno[0]));
       /* observations = countsFunction(posSearch->posSigma[c], 
	  posSearch->posIntervalSizes[c]);*/
       if (0 == compactSearch->pseudoCountConst)
	 columnCounts = columnSpecificPseudocounts(posSearch,compactSearch, c, backgroundProbabilities, observations);
	else
	 columnCounts = compactSearch->pseudoCountConst;
       if (columnCounts >= PSEUDO_MAX) {
	 pseudoWeight = ZERO_OBS_PSEUDO;
	 observations = 0;
       }
       else {
	 pseudoWeight = columnCounts;
       }
       for(a = 0; a < alphabetSize; a++) {
         if (compactSearch->standardProb[a] > posEpsilon) {
	   pseudo = 0;
           /*changed to matrix specific ratios here May 2000*/
	   for (aSub = 0; aSub < alphabetSize; aSub++)
	     if(compactSearch->matrix[a][aSub] != BLAST_SCORE_MIN) 
	       pseudo += (posSearch->posMatchWeights[c][aSub] *
			posSearch->stdFreqRatios->data[a][aSub]);
	   pseudo *= pseudoWeight;
	   numerator = pseudo + 
             (observations * posSearch->posMatchWeights[c][a]/
                compactSearch->standardProb[a]);
	   denominator = observations + pseudoWeight;
	   qOverPEstimate = numerator / denominator;
	   /*Note artificial multiplication by standard probability to
             normalize*/
           posFreqs[c][a] = qOverPEstimate * compactSearch->standardProb[a];
	 if (0.0 != qOverPEstimate && (compactSearch->standardProb[a] > posEpsilon))
	   infoSum += qOverPEstimate * compactSearch->standardProb[a] * log(qOverPEstimate)/ NCBIMATH_LN2;
          }
        else
          posFreqs[c][a] = 0.0;
       }
       if (Cpos)
	 posSearch->posInformation[c] = infoSum;
     }
     else
       for(a = 0; a < alphabetSize; a++) {
         posFreqs[c][a] = 0;
       }
   }
  return(posFreqs);
}

void LIBCALL posScaling(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
	BlastMatrixRescalePtr matrix_rescale;

	matrix_rescale = BlastMatrixRescaleNew(compactSearch->alphabetSize, 
						compactSearch->qlength,
						compactSearch->query,
						compactSearch->standardProb,
						posSearch->posMatrix,
						posSearch->posPrivateMatrix,
						compactSearch->kbp_std,
						compactSearch->kbp_psi,
						compactSearch->kbp_gap_std,
						compactSearch->kbp_gap_psi,
						compactSearch->lambda_ideal,
						compactSearch->K_ideal);

	BlastScaleMatrix(matrix_rescale, TRUE);

	matrix_rescale = BlastMatrixRescaleDestruct(matrix_rescale);

	return;
}


Int4Ptr * LIBCALL CposComputation(posSearchItems *posSearch, BlastSearchBlkPtr search, compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns, Char *ckptFileName, Boolean patternSearchStart, Int4 scorematOutput, Bioseq *query_bsp, Int4 gap_open, Int4 gap_extend, ValNodePtr * error_return, Nlm_FloatHi weightExponent)
{
    Int4 numalign, numseq; /*number of alignments and matches in previous round*/
    
    search->posConverged = FALSE;
    /*  if (patternSearchStart)
        posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, posSearch->posNumSequences);
        else {
    */
    numalign = countSeqAligns(listOfSeqAligns, &numseq, FALSE, 0.0);
    countNumSeq(posSearch, compactSearch, listOfSeqAligns, &numseq);
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numseq);
    
    if (!patternSearchStart)
        findThreshSequences(posSearch, search, listOfSeqAligns, numalign, numseq);
    posDemographics(posSearch, compactSearch, listOfSeqAligns);
    posPurgeMatches(posSearch, compactSearch, error_return);
    posComputeExtents(posSearch, compactSearch);
    posComputeSequenceWeights(posSearch, compactSearch, weightExponent);
    posCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, TRUE);
    if (NULL == search->sbp->posFreqs)
      search->sbp->posFreqs =  allocatePosFreqs(compactSearch->qlength, compactSearch->alphabetSize);
    copyPosFreqs(posSearch->posFreqs,search->sbp->posFreqs, compactSearch->qlength, compactSearch->alphabetSize);
    if (NULL != ckptFileName) {
      if (scorematOutput == NO_SCOREMAT_IO)
        posTakeCheckpoint(posSearch, compactSearch, ckptFileName, error_return);
      else
        posTakeScoremat(posSearch, compactSearch, ckptFileName,
                        scorematOutput, query_bsp, gap_open,
                        gap_extend, error_return);
    }
    posFreqsToMatrix(posSearch,compactSearch);
    posScaling(posSearch, compactSearch);
    return posSearch->posMatrix;
}

/* Top-level routine to compute position-specific matrix, when used through
the Web, one round at a time*/
Int4Ptr * LIBCALL WposComputation(compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns, Nlm_FloatHi ** posFreqs)
{
    posSearchItems *posSearch;
    Int4 i, numSeqAligns, numseq, qlength;
    Int2 alphabetSize;
    Int4Ptr *posMatrix;
    
    /* Why isn't posAllocateMemory() called? */
    posSearch = (posSearchItems *) MemNew(1 * sizeof(posSearchItems));
    qlength = compactSearch->qlength;
    alphabetSize = compactSearch->alphabetSize;

    if (listOfSeqAligns != NULL) {
       numSeqAligns = countSeqAligns(listOfSeqAligns, &numseq, FALSE, 0.0);
       countNumSeq(posSearch, compactSearch, listOfSeqAligns, &numseq);
       posAllocateMemory(posSearch, alphabetSize, 
                         qlength, numseq);
       posDemographics(posSearch, compactSearch, listOfSeqAligns);
       posPurgeMatches(posSearch, compactSearch, NULL);
       posComputeExtents(posSearch, compactSearch);
       posComputeSequenceWeights(posSearch, compactSearch, 1.0);
       posCheckWeights(posSearch, compactSearch);
       posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, 
                                                   FALSE);
       copyPosFreqs(posSearch->posFreqs,posFreqs, qlength, alphabetSize);
    } else {
       /* Assume that posFreqs are already filled, use them to calculate
          posMatrix.
          If listOfSeqAligns is NULL as a result of search, all frequencies are
          0 anyway, so there is no need to compute them. However if it is
          deliberately passed as NULL before search, this means that posFreqs
          are passed as input from PSSM.
       */
       posSearch->posFreqs = posFreqs;
       ASSERT(compactSearch->standardMatrixName);
       posSearch->stdFreqRatios =
           PSIMatrixFrequencyRatiosNew(compactSearch->standardMatrixName);
       posSearch->posMatrix = (BLAST_Score **) MemNew((qlength + 1) * sizeof(BLAST_Score *));
       posSearch->posPrivateMatrix = (BLAST_Score **) MemNew((qlength + 1) * sizeof(BLAST_Score *));
       for(i = 0; i <= qlength; i++) {
          posSearch->posMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
          posSearch->posPrivateMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
       }
    }
    posFreqsToMatrix(posSearch,compactSearch);
    posScaling(posSearch, compactSearch);
    posMatrix = posSearch->posMatrix;
    
    /* Why isn't posFreeMemory() called? */
    if (listOfSeqAligns != NULL) {
       for(i = 0; i <= qlength ; i++) {
          MemFree(posSearch->posFreqs[i]);
          MemFree(posSearch->posMatchWeights[i]);
          MemFree(posSearch->posC[i]);
       }
       
       MemFree(posSearch->posFreqs);
       MemFree(posSearch->posMatchWeights);
       MemFree(posSearch->posC);
       MemFree(posSearch->posA);
       MemFree(posSearch->posExtents);
       MemFree(posSearch->posSigma);
       MemFree(posSearch->posRowSigma);
       MemFree(posSearch->posIntervalSizes);
       MemFree(posSearch->posCount);
       
       for(i = 0; i <= posSearch->posDescMatrixLength; i++) {
          MemFree(posSearch->posDescMatrix[i]);
       }
       MemFree(posSearch->posDescMatrix);
    }

    for(i = 0; i <= qlength ; i++)
       MemFree(posSearch->posPrivateMatrix[i]);
    MemFree(posSearch->posPrivateMatrix);
    posSearch->stdFreqRatios = PSIMatrixFrequencyRatiosFree(posSearch->stdFreqRatios);
    MemFree(posSearch);

    return posMatrix;
}


static Char getRes(Char input)
{
    switch(input) 
      {
      case 0: 
	return('-');
      case 1: 
	return('A');
      case 2: 
	return('B');
      case 3: 
	return('C');
      case 4: 
	return('D');
      case 5: 
	return('E');
      case 6: 
	return('F');
      case 7: 
	return('G');
      case 8: 
	return('H');
      case 9: 
	return('I');
      case 10: 
	return('K');
      case 11: 
	return('L');
      case 12: 
	return('M');
      case 13: 
	return('N');
      case 14: 
	return('P');
      case 15: 
	return('Q');
      case 16: 
	return('R');
      case 17: 
	return('S');
      case 18: 
	return('T');
      case 19: 
	return('V');
      case 20: 
	return('W');
      case 21: 
	return('X');
      case 22: 
	return('Y');
      case 23: 
	return('Z');
      case 24: 
	return('U');
      case 25: 
	return('*');
      case 26: 
	return('O');
      case 27: 
	return('J');
      default:
        return('?');
    }
}
Uint1 LIBCALL ResToInt(Char input)
{
    switch(input) 
      {
      case '-': 
        return(0);
      case 'A': 
        return(1);
      case 'B': 
        return(2);
      case 'C': 
        return(3);
      case 'D': 
        return(4);
      case 'E': 
        return(5);
      case 'F': 
        return(6);
      case 'G': 
        return(7);
      case 'H': 
        return(8); 
      case 'I': 
        return(9);
      case 'K': 
        return(10);
      case 'L': 
        return(11);
      case 'M': 
        return(12);
      case 'N': 
        return(13);
      case 'P': 
        return(14);
      case 'Q': 
        return(15);
      case 'R': 
        return(16);
      case 'S': 
        return(17);
      case 'T': 
        return(18);
      case 'V': 
        return(19);
      case 'W': 
        return(20);
      case 'X': 
        return(21);
      case 'Y': 
        return(22);
      case 'Z': 
        return(23);
      case 'U': 
        return(24);
      case '*': 
        return(25);
      case 'O': 
        return(26);
      case 'J': 
        return(27);
      default:
        return(-1);
    }
}


/*Print out the position-specific matrix*/
void LIBCALL outputPosMatrix(posSearchItems *posSearch, compactSearchItems *compactSearch, FILE *matrixfp, Boolean posComputationCalled)
{
   Uint1Ptr q; /*query sequence*/
   Int4 i; /*loop indices*/
   Int4 c; /*index over alphabet*/
   Int4 length; /*length of query*/
   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/

   if (compactSearch->alphabetSize != PROTEIN_ALPHABET){
     ErrPostEx(SEV_ERROR, 0, 0, "\nCannot print diagnostic information because alphabet size is not %ld", (long) compactSearch->alphabetSize);
     return;
   }
   
   charOrder[0] =  1;  /*A*/
   charOrder[1] =  16; /*R*/
   charOrder[2] =  13; /*N*/  
   charOrder[3] =  4;  /*D*/ 
   charOrder[4] =  3;  /*C*/
   charOrder[5] =  15; /*Q*/
   charOrder[6] =  5;  /*E*/ 
   charOrder[7] =  7;  /*G*/
   charOrder[8] =  8;  /*H*/
   charOrder[9] =  9;  /*I*/
   charOrder[10] = 11; /*L*/
   charOrder[11] = 10; /*K*/
   charOrder[12] = 12; /*M*/  
   charOrder[13] =  6; /*F*/
   charOrder[14] = 14; /*P*/
   charOrder[15] = 17; /*S*/
   charOrder[16] = 18; /*T*/
   charOrder[17] = 20; /*W*/
   charOrder[18] = 22; /*Y*/
   charOrder[19] = 19; /*V*/

   q = compactSearch->query;
   length = compactSearch->qlength;
   
/* Used ifdef until final decision is made on output. */

#ifdef POSIT_DEBUG
   printf("\nCharacter Frequencies by positon\n");
   printf("         ");
   for (c = 0; c< EFFECTIVE_ALPHABET; c++)
      printf("  %c",getRes((Char) charOrder[c]));
   for(i=0; i < length; i++) {
     printf("\n%5d %c   ", i + 1, getRes(q[i]));
     for (c = 0; c < EFFECTIVE_ALPHABET; c++) 
       printf("%2d ", posSearch->posC[i][charOrder[c]]);
   }
   printf("\n\n");
   printf("\nposition counts used. multiplied by 10 and rounded and");
   printf("\nposition character weights used, multiplied by 10 and rounded\n");
   printf("        Counts");
   for (c = 0; c< EFFECTIVE_ALPHABET; c++)
      printf("  %c",getRes((Char) charOrder[c]));
   printf(" Extent ");
   for(i=0; i < length; i++) {
     printf("\n%5d %c   ", i + 1, getRes(q[i]));
     if ((posSearch->posCount[i] > 1) && (Xchar != q[i]))
       printf("%4d ", (Int4) posit_rounddown(10 * countsFunction
				     (posSearch->posSigma[i],posSearch->posIntervalSizes[i])));
     else
       printf("     ");
     for (c = 0; c< EFFECTIVE_ALPHABET; c++)
       if((posSearch->posMatrix[i][charOrder[c]] == BLAST_SCORE_MIN) ||
             (0.0 == posSearch->posMatchWeights[i][charOrder[c]]))
           printf(" - ");
         else
	   printf("%2d ", (Int4) posit_rounddown(10 * posSearch->posMatchWeights[i][charOrder[c]]));
     printf(" %4d",posSearch->posExtents[i].rightExtent - posSearch->posExtents[i].leftExtent +1);
   }
   printf("\n\n");
#endif
   if (NULL != matrixfp) {
     if (posComputationCalled) {
       fprintf(matrixfp,"\nLast position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n");
     }
     else {
       fprintf(matrixfp,"\nLast position-specific scoring matrix computed\n");
     }
     fprintf(matrixfp,"         ");
     for (c = 0; c< EFFECTIVE_ALPHABET; c++)
       fprintf(matrixfp,"  %c",getRes((Char) charOrder[c]));
     if (posComputationCalled) {
       for (c = 0; c< EFFECTIVE_ALPHABET; c++)
	 fprintf(matrixfp,"   %c",getRes((Char) charOrder[c]));
     }
     for(i=0; i < length; i++) {
       fprintf(matrixfp,"\n%5ld %c   ", (long) (i + 1), getRes(q[i]));
       /*fprintf(matrixfp,"\n          ");*/
       for (c = 0; c < EFFECTIVE_ALPHABET; c++) 
	 if(posSearch->posMatrix[i][charOrder[c]] == BLAST_SCORE_MIN)
	   fprintf(matrixfp,"-I ");
	 else
	   fprintf(matrixfp,"%2ld ", (long) posSearch->posMatrix[i][charOrder[c]]);
       if (posComputationCalled) {
	 for (c = 0; c < EFFECTIVE_ALPHABET; c++) 
	   if(posSearch->posMatrix[i][charOrder[c]] != BLAST_SCORE_MIN)
	     fprintf(matrixfp, "%4d", (Int4) posit_rounddown(100 * posSearch->posMatchWeights[i][charOrder[c]]));
	 fprintf(matrixfp," %5.2f", posSearch->posInformation[i]); 
	 if ((posSearch->posCount[i] > 1) && (Xchar != q[i]))
	   fprintf(matrixfp," %.2f", countsFunction(posSearch->posSigma[i],
		     posSearch->posIntervalSizes[i]) * 
                 posSearch->posGaplessColumnWeights[i]/
                 compactSearch->pseudoCountConst);
	 else
	   fprintf(matrixfp,"    0.00");
       }
     }
     fprintf(matrixfp,"\n\n");
     fprintf(matrixfp,"                      K         Lambda\n");
     fprintf(matrixfp,"Standard Ungapped    %6.4f     %6.4f\n",compactSearch->kbp_std[0]->K,compactSearch->kbp_std[0]->Lambda);
     fprintf(matrixfp,"Standard Gapped      %6.4f     %6.4f\n",compactSearch->kbp_gap_std[0]->K,compactSearch->kbp_gap_std[0]->Lambda);
     fprintf(matrixfp,"PSI Ungapped         %6.4f     %6.4f\n",compactSearch->kbp_psi[0]->K,compactSearch->kbp_psi[0]->Lambda);
     fprintf(matrixfp,"PSI Gapped           %6.4f     %6.4f\n",compactSearch->kbp_gap_psi[0]->K,compactSearch->kbp_gap_psi[0]->Lambda);
   }
}


void LIBCALL posPrintInformation(posSearchItems *posSearch, BlastSearchBlkPtr search, Int4 passNum)
{
  Int4 querySize;

  querySize = search->context[0].query->length;

/* Used ifdef until final decision is made on output. */
#ifdef POSIT_DEBUG
  {{
      Int4 c;
  
      printf("\nInformation content by position for pass %d\n", passNum);
      for(c = 0; c < querySize; c++)
          printf(" %5d", c); 
      printf("\n");
      for(c = 0; c < querySize; c++)
          printf(" %5.2lf", posSearch->posInformation[c]); 
      printf("\n");
  }}
#endif
}   
 
void LIBCALL posInitializeInformation(posSearchItems *posSearch, BlastSearchBlkPtr search)
{
  Uint1Ptr query;
  Int4 querySize;
  Int4 c, a, alphabetSize;
  BLAST_ScoreBlkPtr sbp;
  BLAST_ResFreqPtr stdrfp; /*standard frequencies*/
  Nlm_FloatHi lambda;
  Nlm_FloatHi term1, term2, term3, term4;
  Nlm_FloatHi infoSum;
 
  querySize = search->context[0].query->length;
  query = search->context[0].query->sequence;
  posSearch->posInformation = (Nlm_FloatHi *) MemNew(querySize * sizeof(Nlm_FloatHi));
  if (NULL == posSearch->posInformation)
    exit(EXIT_FAILURE);
  for(c = 0; c < querySize; c++)
    posSearch->posInformation[c] = 0.0;
  alphabetSize = search->sbp->alphabet_size;
  /*Compute standard frequencies as in BlastScoreBlkFill in blastkar.c*/
  sbp = search->sbp;
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  lambda = search->sbp->kbp[0]->Lambda;
  for(c = 0; c < querySize; c++) {
    infoSum = 0;
    for(a = 0; a < alphabetSize; a++)
      if (stdrfp->prob[a] > posEpsilon) {
        term1 = search->sbp->matrix[query[c]][a];
	term2 = term1 * lambda;
	term3 = exp(term2);
	term4 = stdrfp->prob[a] * term3;
	infoSum += term4 * log(term4/stdrfp->prob[a])/NCBIMATH_LN2;
      }
    posSearch->posInformation[c] = infoSum;
  }
  BlastResFreqFree(stdrfp);
}

/*
	Is this function used?
*/

void LIBCALL posFreeInformation(posSearchItems *posSearch)
{
  MemFree(posSearch->posInformation);
}

/*Copy a few fields from the lasrge record search into the small record
  compactSearch, so that a small amount of information
  is passed into posit.c*/
void LIBCALL copySearchItems(compactSearchItems * compactSearch, BlastSearchBlkPtr search, char * matrixName)
{
   BLAST_ResFreqPtr stdrfp; /* gets standard frequencies in prob field */
   Int4 a; /*index over characters*/

   compactSearch->query = search->context[0].query->sequence;
   compactSearch->qlength = search->context[0].query->length;
   compactSearch->gapped_calculation = search->pbp->gapped_calculation;
   compactSearch->alphabetSize = search->sbp->alphabet_size;
   compactSearch->pseudoCountConst = search->pbp->pseudoCountConst;
   compactSearch->ethresh = search->pbp->ethresh;
   compactSearch->lambda =  search->sbp->kbp[0]->Lambda;
   compactSearch->matrix = search->sbp->matrix;
   compactSearch->kbp_psi = search->sbp->kbp_psi;
   compactSearch->kbp_gap_psi = search->sbp->kbp_gap_psi;
   compactSearch->kbp_std = search->sbp->kbp_std;
   compactSearch->kbp_gap_std = search->sbp->kbp_gap_std;
   if (search->pbp->gapped_calculation) {
     compactSearch->lambda_ideal = search->sbp->kbp_ideal->Lambda;
     compactSearch->K_ideal = search->sbp->kbp_ideal->K;
   }
   else {
     compactSearch->lambda_ideal = search->sbp->kbp[0] ->Lambda;
     compactSearch->K_ideal = search->sbp->kbp[0]->K;
   }
   compactSearch->use_best_align = search->pbp->use_best_align;

   stdrfp = BlastResFreqNew(search->sbp);
   BlastResFreqStdComp(search->sbp,stdrfp); 
   compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
   if (NULL == compactSearch->standardProb)
     exit(EXIT_FAILURE);
   for(a = 0; a < compactSearch->alphabetSize; a++)
     compactSearch->standardProb[a] = stdrfp->prob[a];
   stdrfp = BlastResFreqDestruct(stdrfp);
   strcpy(compactSearch->standardMatrixName,matrixName);
}

/*allocate memory for a record of type compactSearchItems*/
compactSearchItems * LIBCALL  compactSearchNew(compactSearchItems * compactSearch)
{
   compactSearch = MemNew(1 * sizeof(compactSearchItems));
   if (NULL == compactSearch)
     exit(EXIT_FAILURE);
   return(compactSearch);
}

/*De-allocate memory for a record of type compactSearchItems*/
void LIBCALL compactSearchDestruct(compactSearchItems * compactSearch)
{

   MemFree(compactSearch->standardProb);
   MemFree(compactSearch);
}

/*Some of the following checkpointing code is taken and adapted from
code written by K. Shriram for FASTLINK.
Reference:
 A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. 
 Avoiding Recomputation in Linkage Analysis,
 Human Heredity 44(1994), pp. 225-237. */


#define  putCkptNlm_FloatHi(d, ckptFile)  (putCkptNumber(&(d),sizeof(Nlm_FloatHi),ckptFile))
#define  putCkptInt4(i, ckptFile)         (putCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  putCkptChar(c, ckptFile)         (putCkptNumber(&(c),sizeof(Char),ckptFile))
 
/* General routine for putting the internal representation of a number. */
 
static void  putCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
{
  FileWrite(numberPtr,numberSize,1,ckptFile) ;
}

/*Code to put a vector of frequencies; put only the interesting
  entries*/
static void  putFreqVector(Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
   int  vectorRef;
   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/


   charOrder[0] =  1;  /*A*/
   charOrder[1] =  16; /*R*/
   charOrder[2] =  13; /*N*/  
   charOrder[3] =  4;  /*D*/ 
   charOrder[4] =  3;  /*C*/
   charOrder[5] =  15; /*Q*/
   charOrder[6] =  5; /*E*/ 
   charOrder[7] =  7;  /*G*/
   charOrder[8] =  8;  /*H*/
   charOrder[9] =  9;  /*I*/
   charOrder[10] = 11; /*L*/
   charOrder[11] = 10; /*K*/
   charOrder[12] = 12; /*M*/  
   charOrder[13] =  6; /*F*/
   charOrder[14] = 14; /*P*/
   charOrder[15] = 17; /*S*/
   charOrder[16] = 18; /*T*/
   charOrder[17] = 20; /*W*/
   charOrder[18] = 22; /*Y*/
   charOrder[19] = 19; /*V*/

 
   for(vectorRef = 0; vectorRef < EFFECTIVE_ALPHABET; vectorRef++)
     putCkptNlm_FloatHi(theVector[charOrder[vectorRef]],ckptFile);
}

 
/* Code to put a matrix, vector-by-vector. */
static void    putCkptFreqMatrix (Nlm_FloatHi **theMatrix, Int4 length, Int4 width, FILE * ckptFile)
{
  int  matrixRef;  /*loop index*/
 
  for (matrixRef = 0; matrixRef < length ; matrixRef++ )
    putFreqVector(theMatrix[matrixRef], width, ckptFile);
}

 
/* General routine for getting the internal representation of a number. */
 
void  LIBCALL getCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
{
  FileRead(numberPtr,numberSize,1,ckptFile) ;
}

static void    getFreqVector (Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
   int  vectorRef ;

   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/


   charOrder[0] =  1;  /*A*/
   charOrder[1] =  16; /*R*/
   charOrder[2] =  13; /*N*/  
   charOrder[3] =  4;  /*D*/ 
   charOrder[4] =  3;  /*C*/
   charOrder[5] =  15; /*Q*/
   charOrder[6] =  5; /*E*/ 
   charOrder[7] =  7;  /*G*/
   charOrder[8] =  8;  /*H*/
   charOrder[9] =  9;  /*I*/
   charOrder[10] = 11; /*L*/
   charOrder[11] = 10; /*K*/
   charOrder[12] = 12; /*M*/  
   charOrder[13] =  6; /*F*/
   charOrder[14] = 14; /*P*/
   charOrder[15] = 17; /*S*/
   charOrder[16] = 18; /*T*/
   charOrder[17] = 20; /*W*/
   charOrder[18] = 22; /*Y*/
   charOrder[19] = 19; /*V*/
 
  for(vectorRef = 0; vectorRef < length; vectorRef++)
    theVector[vectorRef] = 0;
  for(vectorRef = 0; vectorRef < EFFECTIVE_ALPHABET; vectorRef++)
    getCkptNlm_FloatHi(theVector[charOrder[vectorRef]],ckptFile) ;
}

/* Code to frequency matrix, vector-by-vector. */
 
void    LIBCALL getCkptFreqMatrix (Nlm_FloatHi ** theMatrix, Int4 length, Int4 width, FILE * ckptFile)
{
  Int4  matrixRef;  /*loop index*/
 
  for (matrixRef = 0; matrixRef < length ; matrixRef++ )
    getFreqVector(theMatrix[matrixRef], width, ckptFile);
}

/*Take a checkpoint at the end of the current PSI-BLAST round, stores
query length, query, and position-specific target frequencies.
Returns TRUE if checkpoint was sucessful and FALSE otherwise. */
Boolean LIBCALL posTakeCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr *error_return)
{
  FILE * checkFile; /*file in which to take the checkpoint*/
  Int4 length; /*length of query sequence, and an index for it*/
  Int4 i; /*indices to position and alphabet */
  Char localChar; /*temporary character*/

  checkFile = FileOpen(fileName, "wb");
  if (NULL == checkFile) {
    BlastConstructErrorMessage("posTakeCheckpoint", "Could not open checkpoint file", 1, error_return);
    return(FALSE);
  }
  length = compactSearch->qlength;
  putCkptInt4(length, checkFile);
  for(i = 0; i < length; i++) {
    localChar = getRes(compactSearch->query[i]);
    putCkptChar(localChar, checkFile);
  }  
  putCkptFreqMatrix(posSearch->posFreqs,length,compactSearch->alphabetSize, checkFile);
  FileClose(checkFile);
  return(TRUE);
}

/* Like posTakeCheckpoint, posTakeScoremat will emit the position
   frequencies that have been generated. Unlike that routine, the
   file to be written is an ASN.1 encoded PssmWithParameters object. */

Boolean LIBCALL posTakeScoremat(posSearchItems *posSearch, 
                       compactSearchItems *compactSearch, 
                       CharPtr filename, Int4 scorematOutput,
                       Bioseq *query_bsp, Int4 gap_open, 
                       Int4 gap_extend, ValNodePtr *error_return)
{
  AsnIoPtr outfile = NULL;
  PssmWithParametersPtr scoremat = NULL;
  PssmIntermediateDataPtr freqs = NULL;
  PssmParametersPtr params = NULL;
  FormatRpsDbParametersPtr rpsparams = NULL;
  PssmPtr pssm = NULL;
  Int4 i, j;
  Boolean status = FALSE;

  scoremat = PssmWithParametersNew();
  if (scoremat == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not allocate PssmWithParameters", 1, error_return);
    goto bail_out;
  }

  /* Add information about the underlying score matrix.
     Note that blastpgp will ignore this information */

  params = scoremat->params = PssmParametersNew();
  if (params == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not allocate PssmParameters", 1, error_return);
    goto bail_out;
  }
  rpsparams = params->rpsdbparams = FormatRpsDbParametersNew();
  if (params == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not allocate RpsDbParameters", 1, error_return);
    goto bail_out;
  }
  rpsparams->matrixName = strdup(compactSearch->standardMatrixName);
  rpsparams->gapOpen = gap_open;
  rpsparams->gapExtend = gap_extend;

  /* Build up the objects describing the frequency ratios */

  pssm = scoremat->pssm = PssmNew();
  if (pssm == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not allocate PSSM object", 1, error_return);
    goto bail_out;
  }
  freqs = pssm->intermediateData = PssmIntermediateDataNew();
  if (freqs == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not allocate PssmIntermediateData", 1, error_return);
    goto bail_out;
  }

  pssm->isProtein = TRUE;
  pssm->numRows = compactSearch->alphabetSize;
  pssm->numColumns = compactSearch->qlength;

  for (i = 0; i < pssm->numColumns; i++) {
    for (j = 0; j < pssm->numRows; j++) {
      ValNodeAddFloat(&freqs->freqRatios, 0, posSearch->posFreqs[i][j]);
    }
  }

  /* Do not make a copy of the query bioseq; use it directly.
     The '1' below indicates a single bioseq (not a bioseq-set) */

  ValNodeAddPointer(&pssm->query, 1, query_bsp);
  if (pssm->query == NULL) {
    BlastConstructErrorMessage("posTakeScoremat", 
               "Could not attach bioseq to scoremat", 1, error_return);
    goto bail_out;
  }

  if (scorematOutput == ASCII_SCOREMAT)
     outfile = AsnIoOpen(filename, "w");
  else
     outfile = AsnIoOpen(filename, "wb");

  if (outfile == NULL) {
    ErrPostEx(SEV_FATAL, 1, 0, "Unable to open matrix output file %s\n", 
          filename);
    goto bail_out;
  }

  PssmWithParametersAsnWrite(scoremat, outfile, NULL);
  status = TRUE;

bail_out:
  AsnIoClose(outfile);

  /* explicitly free the ValNode pointing to the query bioseq.
     This will prevent the ScoreMatrix freeing routine
     from also freeing the query bioseq, which we did not
     allocate */

  pssm->query = ValNodeFree(pssm->query);

  /* free everything else */

  scoremat = PssmWithParametersFree(scoremat);
  return status;
}

static Boolean LIBCALL posReadPosFreqsScoremat(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, Int4 scorematInput, ValNodePtr * error_return)
{
  AsnIoPtr infile = NULL;
  PssmWithParametersPtr scoremat = NULL;
  PssmPtr pssm = NULL;
  PssmIntermediateDataPtr freqs = NULL;
  Int4 i, j, c;
  ValNodePtr freq_list;
  Bioseq *bsp;

  if (scorematInput == ASCII_SCOREMAT)
     infile = AsnIoOpen(fileName, "r");
  else
     infile = AsnIoOpen(fileName, "rb");

  if (infile == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0,"Could not open scoremat file\n");
    return FALSE;
  }

  scoremat = PssmWithParametersAsnRead(infile, NULL);
  AsnIoClose(infile);
  if (scoremat == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0, "Could not read scoremat from input file\n");
    return FALSE;
  }
  pssm = scoremat->pssm;
  if (pssm == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0,"Scoremat is empty\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  freqs = pssm->intermediateData;
  if (freqs == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0,"Scoremat doesn't contain intermediate data\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (freqs->freqRatios == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0,
            "Scoremat does not contain frequency ratios\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (pssm->numRows != compactSearch->alphabetSize) {
    ErrPostEx(SEV_WARNING, 0, 0, "Wrong alphabet size of %d in "
              "input scoremat\n", pssm->numRows);
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (!pssm->query || !pssm->query->data.ptrvalue) {
    ErrPostEx(SEV_WARNING, 0, 0, "Missing sequence data in input scoremat\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  bsp = (Bioseq *)(pssm->query->data.ptrvalue);
  if (pssm->numColumns != bsp->length) {
    ErrPostEx(SEV_WARNING, 0, 0, "Different sequence lengths "
              "(%d and %d) in input scoremat\n", pssm->numColumns, bsp->length);
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (pssm->numColumns != compactSearch->qlength) {
    ErrPostEx(SEV_WARNING, 0, 0, "Scoremat sequence length "
              "(%d) does not match query length (%d)\n", 
              pssm->numColumns, compactSearch->qlength);
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (!bsp->seq_data || !ISA_aa(bsp->mol)) {
    ErrPostEx(SEV_WARNING, 0, 0, 
          "Sequence within checkpoint file has no data or is not protein\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (bsp->seq_data_type == Seq_code_gap) {
    ErrPostEx(SEV_WARNING, 0, 0, 
          "Seq_code_gap passed to posReadPosFreqsScoremat\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  BSSeek((ByteStorePtr) bsp->seq_data, 0, SEEK_SET);

  /* Convert sequence data into Seq_code_ncbistdaa */
  if (bsp->seq_data_type != Seq_code_ncbistdaa) {

      ByteStore* new_byte_store = BSConvertSeq((ByteStorePtr) bsp->seq_data,
                                               Seq_code_ncbistdaa,
                                               bsp->seq_data_type,
                                               bsp->length);

      if ( !new_byte_store ) {
          ErrPostEx(SEV_FATAL, 1, 0, "Failed to convert Bioseq in ASN.1 PSSM "
                    "to Seq_code_ncbistdaa");
      }

      bsp->seq_data = (SeqDataPtr) new_byte_store;
      bsp->seq_data_type = Seq_code_ncbistdaa;
      BSSeek((ByteStorePtr) bsp->seq_data, 0, SEEK_SET);

  }

  /* verify the input query is the same as the sequence
     within the checkpoint file */

  for (i = 0; i < compactSearch->qlength; i++) {
    c = BSGetByte((ByteStorePtr) bsp->seq_data);
    if (c == EOF) {
      ErrPostEx(SEV_WARNING, 0, 0, "Premature end of sequence data\n");
      PssmWithParametersFree(scoremat);
      return FALSE;
    }
    if (c != compactSearch->query[i]) {
      if (compactSearch->query[i] == Xchar) {
        ErrPostEx(SEV_WARNING, 0, 0, 
                     "Query sequence contains '%c' at position %d; "
                     "if filtering was used, rerun the search with "
                     "filtering turned off ('-F F')\n", getRes(Xchar), i);
      }
      else {
        ErrPostEx(SEV_WARNING, 0, 0, 
                     "Query sequence contains '%c' at position %d, "
                     "while sequence withing checkpoint file contains "
                     "'%c' at this position\n", 
                     getRes(compactSearch->query[i]), i, getRes(c));
      }
      PssmWithParametersFree(scoremat);
      return FALSE;
    }
  }

  /* Read in the frequency ratios, verify they fall
     in the correct range, and verify that the linked list
     of residue frequencies is exactly as long as it should be */

  freq_list = freqs->freqRatios;
  if (pssm->byRow == FALSE) {
    for (i = 0; i < pssm->numColumns; i++) {
      for (j = 0; j < pssm->numRows; j++) {
        if (freq_list == NULL)
          break;
        posSearch->posFreqs[i][j] = freq_list->data.realvalue;

        if (posSearch->posFreqs[i][j] < 0.0) {
          ErrPostEx(SEV_WARNING, 0, 0, "position frequency (%d,%d) "
                    "out of bounds\n", i, j);
          PssmWithParametersFree(scoremat);
          return FALSE;
        }

        freq_list = freq_list->next;
      }
      if (j < pssm->numRows)
        break;
    }
  }
  else {
    for (j = 0; j < pssm->numRows; j++) {
      for (i = 0; i < pssm->numColumns; i++) {
        if (freq_list == NULL)
          break;
        posSearch->posFreqs[i][j] = freq_list->data.realvalue;

        if (posSearch->posFreqs[i][j] < 0.0) {
          ErrPostEx(SEV_WARNING, 0, 0, "position frequency (%d,%d) "
                    "out of bounds\n", i, j);
          PssmWithParametersFree(scoremat);
          return FALSE;
        }

        freq_list = freq_list->next;
      }
      if (i < pssm->numColumns)
        break;
    }
  }

  if (i < pssm->numColumns || j < pssm->numRows) {
    ErrPostEx(SEV_WARNING, 0, 0, "Not enough frequency "
              "ratios in input scoremat\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  if (freq_list != NULL) {
    ErrPostEx(SEV_WARNING, 0, 0, "Too many frequency "
              "ratios in input scoremat\n");
    PssmWithParametersFree(scoremat);
    return FALSE;
  }
  PssmWithParametersFree(scoremat);
  return TRUE;
}

static Boolean posReadPosFreqsStandard(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return)
{
  FILE * checkFile; /*file in which to take the checkpoint*/
  Int4 length1, length2, c; /*length of query sequence, and an index for it*/
  Char  nextRes; /*next residue in stored copy of the query sequence*/
  Uint1Ptr oldQuery; /*array to hold the query sequence*/

  length1 = compactSearch->qlength;

  checkFile = FileOpen(fileName, "rb");  
  if (NULL == checkFile) {
    BlastConstructErrorMessage("posReadPosFreqsStandard", "Could not open checkpoint file\n", 1, error_return);
    return(FALSE);
  }
  getCkptInt4(length2,checkFile);
  if (length1 != length2) {
    ErrPostEx(SEV_WARNING, 0, 0, "Invalid usage of checkpoint recovery; old query has length %ld, new query has length %ld", (long) length2,  (long) length1);
    BlastConstructErrorMessage("posReadPosFreqsStandard", "Failed to recover data\n", 1, error_return);
    FileClose(checkFile);
    return(FALSE);
  }
  oldQuery = (Uint1Ptr) MemNew(length1 * sizeof(Uint1));
  if (NULL == oldQuery) {
    BlastConstructErrorMessage("posReadPosFreqsStandard", "Failed to reconstruct previous query\n", 1, error_return);
    BlastConstructErrorMessage("posReadPosFreqsStandard", "Failed to recover data\n", 1, error_return);
    FileClose(checkFile);
    return(FALSE);
  }  
  for(c = 0; c < length1; c++) {
    getCkptChar(nextRes, checkFile);
    oldQuery[c] = ResToInt(nextRes);


    if ((oldQuery[c] != compactSearch->query[c]) && (oldQuery[c] != Xchar)) {
                                /* Error massage Added by Natsuhiko */
      if (compactSearch->query[c] == Xchar) {
        ErrPostEx(SEV_WARNING, 0, 0, "\nStored query has a %c at position %ld, while new query has a %c there.\n%c appears in query sequence: The query could be filtered. Run with \"-F F\" option to turn the filter off.",getRes(oldQuery[c]), (long) c, getRes(compactSearch->query[c]),
                  getRes(compactSearch->query[c]));
      }
      else{
      ErrPostEx(SEV_WARNING, 0, 0, "Stored query has a %c at position %ld, while new query has a %c there",getRes(oldQuery[c]), (long) c, getRes(compactSearch->query[c]));
      }

      BlastConstructErrorMessage("posReadPosFreqsStandard", "Failed to recover data\n", 1, error_return);
      MemFree(oldQuery);
      FileClose(checkFile);
      return(FALSE);
    }
    if ((oldQuery[c] != compactSearch->query[c]) && (Xchar==oldQuery[c])) {
      ErrPostEx(SEV_WARNING, 0, 0, "Stored query has a %c at position %ld, while new query has a %c there\n%c appears in the stored query: The stored query may be filtered. Run blastpgp with \"-F F\" option to turn the filter off",getRes(oldQuery[c]), (long) c,
                getRes(compactSearch->query[c]), getRes(oldQuery[c])); }

  }
  getCkptFreqMatrix(posSearch->posFreqs,length1,compactSearch->alphabetSize,checkFile);
  MemFree(oldQuery);
  FileClose(checkFile);
  return(TRUE);
}

/*Read a checkpoint from the end of a previous PSI-BLAST round, get
query length, query, and position-specific target frequencies.
Returns TRUE if checkpoint was read sucessfully and FALSE otherwise. */
Boolean LIBCALL  posReadCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, Int4 scorematInput, ValNodePtr * error_return)
{
  Int4 length1;    /*length of query sequence*/
  Int4 i,j;      /*indices for position and character in alphabet*/
  Boolean FreqsRead;

  BlastConstructErrorMessage("posReadCheckpoint", "Attempting to recover data from previous checkpoint\n", 1, error_return);
  length1 = compactSearch->qlength;

  /* allocate memory for the PSSMs and position frequency matrix */

  posSearch->posMatrix = (BLAST_Score **) MemNew((length1 + 1) * sizeof(BLAST_Score *));
  posSearch->posPrivateMatrix = (BLAST_Score **) MemNew((length1 + 1) * sizeof(BLAST_Score *));
  posSearch->posFreqs = (Nlm_FloatHi **) MemNew((length1 + 1) * sizeof(Nlm_FloatHi *));
  ASSERT(compactSearch->standardMatrixName);
  posSearch->stdFreqRatios =
           PSIMatrixFrequencyRatiosNew(compactSearch->standardMatrixName);
  if ((NULL == posSearch->posMatrix) || (NULL == posSearch->posPrivateMatrix) || (NULL == posSearch->posFreqs)) {

    BlastConstructErrorMessage("posReadCheckpoint", "Failed to allocate position-specific score matrix", 1, error_return);
    BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
    return(FALSE);
  }
  for(i = 0; i <= length1; i++) {
    posSearch->posMatrix[i] = (BLAST_Score *) MemNew(compactSearch->alphabetSize * sizeof(BLAST_Score));
    posSearch->posPrivateMatrix[i] = (BLAST_Score *) MemNew(compactSearch->alphabetSize * sizeof(BLAST_Score));
    posSearch->posFreqs[i] = (Nlm_FloatHi *) MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));

    if ((NULL == posSearch->posMatrix[i]) || (NULL == posSearch->posPrivateMatrix[i]) || (NULL == posSearch->posFreqs[i])) {
      BlastConstructErrorMessage("posReadCheckpoint", "Failed to allocate position-specific score matrix", 1, error_return);
      BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
      return(FALSE);
    }
    for(j = 0; j < compactSearch->alphabetSize; j++) {
      posSearch->posFreqs[i][j] = 0.0;
    }
  }

  if (scorematInput == NO_SCOREMAT_IO) {
    FreqsRead = posReadPosFreqsStandard(posSearch, compactSearch, 
                      fileName, error_return);
  }
  else {
    FreqsRead = posReadPosFreqsScoremat(posSearch, compactSearch, 
                      fileName, scorematInput, error_return);
  }
  if (FreqsRead != TRUE) {
    BlastConstructErrorMessage("posReadCheckpoint", "Data recovery failed\n", 
                                                1, error_return);
    return(FALSE);
  }
  posFreqsToInformation(posSearch,compactSearch);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  BlastConstructErrorMessage("posReadCheckpoint", "Data recovered successfully\n", 1, error_return);
  return(TRUE);
}

/* Two routines taken from */
/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */

static Int4 P_eof(FILE *f) 
{
    register Int4 ch;

    if (feof(f))
	return 1;
    if (f == stdin)
	return 0;    /* not safe to look-ahead on the keyboard! */
    ch = getc(f);
    if (ch == EOF)
	return 1;
    ungetc(ch, f);
    return 0;
}

static Boolean isBlankChar(Char thisChar)
{
  return((thisChar == ' ') || (thisChar == '\t') || (thisChar == '\n') ||
          (thisChar == '\t'));
}


/*preprocess alignment checkpoint file to find number of sequences and
number of blocks. Return number of blocks as return value
and number of sequences through a reference parameter*/
static Int4 posFindAlignmentDimensions(char * fileName, Int4 *numSeqs, ValNodePtr * error_return)
{
  FILE *checkFile;  /*checkpoint file*/
  Char nextLine[ALIGN_LINE_LENGTH];  /*line read in*/
  Boolean foundBlankLine; /*have we found a blank line yet*/
  Int4 numBlocks;  /*number of blocks to be returned*/
  Int4 testCountSeqs; /*counts number of sequences in each block
                        to ensure that each block has the same
                        number of sequences*/
   
  BlastConstructErrorMessage("posFindAlignmentDimensions", "Attempting to recover data from multiple alignment file\n", 1, error_return);
  checkFile = FileOpen(fileName, "rb");  
  if (NULL == checkFile) {
    ErrPostEx(SEV_WARNING, 0, 0, "\nCould not open alignment checkpoint file");
    BlastConstructErrorMessage("posFindAlignmentDimensions", "Could not open alignment checkpoint file\n", 1, error_return);
    return(0);
  }
  do {
    fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);
  } while (isBlankChar(nextLine[0]));
  foundBlankLine = FALSE;
  *numSeqs = 1;
  numBlocks = 0;
  while (!P_eof(checkFile) && (!foundBlankLine)) {
    fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);
    if (!(isBlankChar(nextLine[0])))
      (*numSeqs)++;
    else
      foundBlankLine = TRUE;
  }
  numBlocks = 1;
  while(!P_eof(checkFile)) {
    do {
      fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);    
    } while((!P_eof(checkFile)) && (isBlankChar(nextLine[0])));
    if (!P_eof(checkFile)) {
      numBlocks++;
      testCountSeqs = 0;
    }    
    do {
      fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);    
      testCountSeqs++;
    } while((!P_eof(checkFile)) && !(isBlankChar(nextLine[0])));
    if (!(isBlankChar(nextLine[0])))
      testCountSeqs++;
    if (testCountSeqs != (*numSeqs)) {
      ErrPostEx(SEV_WARNING, 0, 0, "\nInconsistent number of sequences across alignment blocks, first block has %d while block %d has %d sequences",(*numSeqs), numBlocks, testCountSeqs);
      BlastConstructErrorMessage("posFindAlignmentDimensions", "Could not read alignment due to different number of sequences in different blocks\n", 1, error_return);
      FileClose(checkFile);
      return(0);
    }
  }

  FileClose(checkFile);
  return(numBlocks);
}

/*Is thisChar possibly part of an alignment?*/
static Boolean isProteinChar(Char thisChar)
{

  return(((thisChar >= 'A') && (thisChar <= 'Z')) ||
         ((thisChar >= 'a') && (thisChar <= 'z')) ||
         ('-' == thisChar));
}


/*preprocess alignment checkpoint file to find the
start column and end column of each alignment block.
As a consequece the length of
the alignment can be computed and it is returned*/
static Int4 posPreprocessAlignment(char * fileName, Int4 numSeqs, Int4 numBlocks, Int4 * numCols, ValNodePtr * error_return)
{
  FILE *checkFile;  /*checkpoint file*/
  char nextLine[ALIGN_LINE_LENGTH];  /*line read in*/
  Int4 alignLength; /*number of columns in alignment, to be returned*/
  Int4 charIndex; /*index over characters in a row*/
  Int4 blockIndex; /*index for the blocks in the alignment file*/
  Int4 seqIndex;
   
  checkFile = FileOpen(fileName, "rb");  
  if (NULL == checkFile) {
    ErrPostEx(SEV_WARNING, 0, 0, "\nCould not open alignment checkpoint file");
    BlastConstructErrorMessage("posPreprocessAlignment", "Could not open alignment checkpoint file\n", 1, error_return);
    return(0);
  }

  blockIndex = 0;
  alignLength= 0;
  while (!P_eof(checkFile)) {
    do {
      fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);
    } while (isBlankChar(nextLine[0])); /*line belongs to query*/
    charIndex = 0;
    while(!(isBlankChar(nextLine[charIndex])))
      charIndex++;
    while(isBlankChar(nextLine[charIndex]))
      charIndex++;
    numCols[blockIndex] = 0;
    while (isProteinChar(nextLine[charIndex])){
      alignLength++;
      charIndex++;
      numCols[blockIndex]++;
    }
    /*skip over other sequences*/
    for (seqIndex = 0; seqIndex < numSeqs; seqIndex++) 
      fgets(nextLine, ALIGN_LINE_LENGTH,checkFile);
    blockIndex++;
  }
  FileClose(checkFile);
  return(alignLength);
}

/*Find the index of the sequence in the multiple alignment that
  matches the query sequence; if non match return -1*/
static Int4 findQuery(posDesc ** alignArray, compactSearchItems * compactSearch, Int4 numSeqs, Int4 alignLength)
{
   Uint1Ptr query; /*query sequence*/
   Int4 qlength;  /*length of query sequence*/
   Int4 seqIndex; /*index over sequences*/
   Int4 i;        /*index within a sequence*/
   Int4 queryIndex;  /*index within query*/
   Char thisRes;

   query = compactSearch->query;
   qlength = compactSearch->qlength;
   for(seqIndex = 0; seqIndex < numSeqs; seqIndex++) {
     i = 0;
     queryIndex = 0;
     while ((queryIndex < qlength) && (i < alignLength)) {
       if ('-' == alignArray[seqIndex][i].letter)
         i++;
       else {
         /*Need to keep lower-case letters*/
         thisRes = getRes(query[queryIndex]);
         /* Selenocysteines are replaced by X's in query; test for this
            possibility */
         if ((alignArray[seqIndex][i].letter == 'U' ||
             alignArray[seqIndex][i].letter == 'u') &&
             thisRes == 'X')
            thisRes = alignArray[seqIndex][i].letter;
            
         if ((thisRes != (alignArray[seqIndex][i].letter + 'A' - 'a')) &&
             (thisRes != alignArray[seqIndex][i].letter))
           /*character mismatch*/
           break;
         else {
           queryIndex++;
           i++;
         }
       }
     }
     if (queryIndex == qlength) {
       while (i < alignLength) {
         /*chew up gaps at end of alignment sequence*/
         if ('-' != alignArray[seqIndex][i].letter)
           break;
         i++;
       }
       /*found a match! */
       return(seqIndex);
     }
     else
       /*alignment string is prefix of query*/
       continue;
   }
   return (-1);

}


static posDesc** posReadAlignment(compactSearchItems *compactSearch, char * fileName, Int4 numSeqs, Int4 numBlocks, Int4 alignLength, Int4 * numCols,
ValNodePtr * error_return)
{
  FILE *checkFile; /*checkpoint file to read*/
  Char nextline[ALIGN_LINE_LENGTH];
  Int4 blockIndex;
  Int4 linePos; /* moving index for a line*/
  Int4 alignPos; /*placeholder for position alignment*/
  Int4 base; /*base for this block*/
  posDesc **returnArray; /*array of sequences to retunr*/
  Int4 i,j; /*loop indices*/
  Int4 temp; /*temporary character for swapping sequences*/
  Int4 queryIndex; /*which sequnec in the alignment is the query*/

  checkFile = FileOpen(fileName, "rb");  
  if (NULL == checkFile) {
    BlastConstructErrorMessage("posReadAlignment", "Could not open alignment checkpoint file\n", 1, error_return);
    ErrPostEx(SEV_WARNING, 0, 0, "\nCould not open alignment checkpoint file");
  }
  returnArray = (posDesc**) MemNew(numSeqs * sizeof(posDesc *));
  if (NULL == returnArray)
    exit(EXIT_FAILURE);
  for (i = 0; i < numSeqs; i++) {
    returnArray[i] = (posDesc *) MemNew(alignLength * sizeof(posDesc));
    if (NULL == returnArray[i])
      exit(EXIT_FAILURE);
    for(j = 0; j < alignLength; j++) {
      returnArray[i][j].letter = UNUSED;
      returnArray[i][j].used = FALSE;
    }
  }  
  alignPos = 0;
  base = 0;
  for(blockIndex = 0; blockIndex < numBlocks; blockIndex++){
    for(i = 0; i < numSeqs; i++) {
      do {
	fgets(nextline, ALIGN_LINE_LENGTH,checkFile);    
      } while(isBlankChar(nextline[0]));
      linePos = 0;
      while(!isBlankChar(nextline[linePos]))
        linePos++;
      while(isBlankChar(nextline[linePos]))
        linePos++;
      alignPos = base;
      while (alignPos < (base + numCols[blockIndex])) {
        if (!isProteinChar(nextline[linePos])) {
          BlastConstructErrorMessage("posReadAlignment", "Invalid character or wrong number of characters in a sequence\n", 1, error_return);
          ErrPostEx(SEV_WARNING, 0, 0, "\nInvalid character or wrong number of characters in sequence index %d\n", i+1);
        }
	returnArray[i][alignPos].letter = nextline[linePos];
	returnArray[i][alignPos].used = TRUE;
        alignPos++;
        linePos++;
      }
    }
    base += numCols[blockIndex];
  }
  FileClose(checkFile);
  queryIndex = findQuery(returnArray, compactSearch, numSeqs, alignLength);
  if (-1 == queryIndex) {
    BlastConstructErrorMessage("posReadAlignment", "None of the alignment sequences equals the query sequence\n", 1, error_return);
    ErrPostEx(SEV_WARNING, 0, 0, "\nNone of the alignment sequences equals the query sequence");
    BlastConstructErrorMessage("posReadAlignment", "Cannot recover alignment checkpoint\n", 1, error_return);
    ErrPostEx(SEV_WARNING, 0, 0, "\nCannot recover alignment checkpoint");
    exit(EXIT_FAILURE);
  }
  else {
    if (queryIndex > 0) {
      /*swap query with first sequence in alignment*/
      for (alignPos = 0; alignPos < alignLength; alignPos++) {
        temp = returnArray[0][alignPos].letter;
        returnArray[0][alignPos].letter = returnArray[queryIndex][alignPos].letter;
        returnArray[queryIndex][alignPos].letter = temp;
      }
    }
  }
  return(returnArray);
}

static void posProcessAlignment(posSearchItems *posSearch, compactSearchItems *compactSearch, char * fileName, Int4 numSeqs, Int4 numBlocks, Int4 alignLength, Int4 * numCols, ValNodePtr * error_return)
{
  Int4 queryPos, alignPos, linePos; /*placeholder for position in query and
                             alignment*/
  Int4 *queryDesc; /*position correspondence between alignment and query*/
  Int4 seqIndex; /*counter for sequences*/
  posDesc ** alignArray;
  Int4 queryLength; /*length of query sequence*/

  alignArray = posReadAlignment(compactSearch, fileName, numSeqs,  numBlocks, alignLength, numCols, error_return);
  queryDesc = (Int4 *) MemNew(alignLength * sizeof(Int4));
  if (NULL == queryDesc)
    exit(EXIT_FAILURE);
  for(alignPos = 0; alignPos < alignLength; alignPos++)
    queryDesc[alignPos] = GAP_HERE;
  alignPos = 0;
  queryPos = 0;
  for(linePos = 0; linePos < alignLength; linePos++) {
    if (alignArray[0][linePos].letter == '-') 
      queryDesc[alignPos] = GAP_HERE;
    else {
      queryDesc[alignPos] = queryPos;
      if ((alignArray[0][linePos].letter >= 'A' ) && (alignArray[0][linePos].letter <= 'Z')) {
	posSearch->posDescMatrix[0][queryPos].letter = ResToInt(alignArray[0][linePos].letter);
	posSearch->posDescMatrix[0][queryPos].used = TRUE;
	posSearch->posDescMatrix[0][queryPos].e_value = compactSearch->ethresh/2;
	posSearch->posDescMatrix[0][queryPos].leftExtent = 0;
	posSearch->posDescMatrix[0][queryPos].rightExtent = compactSearch->qlength - 1;
	posSearch->posC[queryPos][ResToInt(alignArray[0][linePos].letter)]++;
	posSearch->posCount[queryPos]++;
      }
      else {
	posSearch->posDescMatrix[0][queryPos].used = FALSE;
	posSearch->posDescMatrix[0][queryPos].letter = ResToInt(alignArray[0][linePos].letter + 'A' - 'a');
      }
      queryPos++;
    }
    alignPos++;
  }
  queryLength = queryPos;
  for(seqIndex = 1; seqIndex < numSeqs; seqIndex++) {
    for(linePos = 0; linePos < alignLength; linePos++) {
      if (queryDesc[linePos] != GAP_HERE) {
	if (!(posSearch->posDescMatrix[0][queryDesc[linePos]].used)) {
	  /*mark column as not participating*/
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].used = FALSE;
	}
	else {
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].letter = ResToInt(alignArray[seqIndex][linePos].letter);
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].used = TRUE;
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].e_value = compactSearch->ethresh/2;
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].leftExtent = 0;
	  posSearch->posDescMatrix[seqIndex][queryDesc[linePos]].rightExtent = compactSearch->qlength;
	}
      }
    }
  }
  /*make terminal gaps unused*/
  for(seqIndex = 1; seqIndex < numSeqs; seqIndex++) {
    linePos = 0;
    while((linePos < queryLength) && (posSearch->posDescMatrix[seqIndex][linePos].letter == GAP_CHAR)) {
      posSearch->posDescMatrix[seqIndex][linePos].used = FALSE;
      linePos++;
    }
    linePos = queryLength - 1;
    while((linePos >= 0) && (posSearch->posDescMatrix[seqIndex][linePos].letter == GAP_CHAR)) {
      posSearch->posDescMatrix[seqIndex][linePos].used = FALSE;
      linePos--;
    }
  }
  BlastConstructErrorMessage("posProcessAlignment", "Alignment recovered successfully\n", 1, error_return);
}


/* Top-level routine to compute position-specific matrix, when used for
one round to recover from a multiple alignment checkpoint. */
Int4Ptr * LIBCALL BposComputation(posSearchItems *posSearch, BlastSearchBlkPtr
   search, compactSearchItems * compactSearch, Char *ckptFileName, 
   Char *takeCkptFileName, Int4 scorematOutput, Bioseq *query_bsp,
   Int4 gap_open, Int4 gap_extend, ValNodePtr * error_return)
{
  Int4 numSeqs, numBlocks, alignLength; /*number of sequences, number of pieces
                        of alignment, total length of the alignment*/
  Int4 *numCols;  /*number of columns within each block*/

  search->posConverged = FALSE;

  numBlocks = posFindAlignmentDimensions(ckptFileName, &numSeqs, error_return);
  if (0 == numBlocks) {
    ErrPostEx(SEV_WARNING, 0, 0, "\nCould not recover block structure from checkpoint");
    BlastConstructErrorMessage("BposComputation", "Cannot recover alignment checkpoint\n", 1, error_return);
    return(NULL);
  }
  numCols = (Int4 *) MemNew(numBlocks * sizeof(Int4));
  if (NULL == numCols)
    exit(EXIT_FAILURE);
  alignLength = posPreprocessAlignment(ckptFileName,  numSeqs, numBlocks,  numCols, error_return);
  if (0 == alignLength) {
    ErrPostEx(SEV_WARNING, 0, 0, "\nCould not recover alignment structure from checkpoint");
    BlastConstructErrorMessage("BposComputation", "Cannot recover alignment checkpoint\n", 1, error_return);
    return(NULL);
  } 
  posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numSeqs);
  posProcessAlignment(posSearch, compactSearch, ckptFileName, numSeqs,  numBlocks, alignLength, numCols, error_return);
  MemFree(numCols);
  posSearch->posNumSequences = numSeqs;
  posPurgeMatches(posSearch, compactSearch, error_return);
  posComputeExtents(posSearch, compactSearch);
  posComputeSequenceWeights(posSearch, compactSearch, 1.0);
  posCheckWeights(posSearch, compactSearch);
  posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, TRUE);
  if (NULL == search->sbp->posFreqs)
    search->sbp->posFreqs =  allocatePosFreqs(compactSearch->qlength, compactSearch->alphabetSize);
  copyPosFreqs(posSearch->posFreqs,search->sbp->posFreqs, compactSearch->qlength, compactSearch->alphabetSize);
  if (NULL != takeCkptFileName) {
    if (scorematOutput == NO_SCOREMAT_IO)
      posTakeCheckpoint(posSearch, compactSearch, takeCkptFileName, error_return);
    else
      posTakeScoremat(posSearch, compactSearch, takeCkptFileName, scorematOutput, query_bsp, gap_open, gap_extend, error_return);
  }
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  return posSearch->posMatrix;
}

/****************************************************************************/
/* PLEASE NOTE: The following structure and the PSIMatrixFrequencyRatios*
 * functions have been copied and renamed from
 * algo/blast/core/matrix_freq_ratios.[hc] to eliminate a dependency from the
 * ncbitool library to the blast library. 
 */
#ifndef BLOSUM62_20A_SCALE_MULTIPLIER
#define BLOSUM62_20A_SCALE_MULTIPLIER 0.9666
#endif

#ifndef BLOSUM62_20B_SCALE_MULTIPLIER
#define BLOSUM62_20B_SCALE_MULTIPLIER 0.9344
#endif

/* posit.c has static frequency data for the 28 character alphabet,
 * which may be different from PROTEIN_ALPHABET */

#define POSITAA_SIZE 28

/** Underlying frequency ratios for BLOSUM45 */
static const double BLOSUM45_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 2.95043377e+00, 7.34701741e-01, 8.00397827e-01,
  6.88936672e-01, 8.25164920e-01, 5.87357723e-01, 1.08031132e+00,
  6.54086288e-01, 7.46806187e-01, 7.86209397e-01, 7.12370041e-01,
  8.21348665e-01, 7.89043130e-01, 7.08569419e-01, 8.66678731e-01,
  6.99695540e-01, 1.30031418e+00, 1.00058530e+00, 1.00992663e+00,
  5.65442334e-01, 7.50000000e-01, 6.38727873e-01, 8.41025176e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.26064996e-01},
 {0.00000000e+00, 7.34701741e-01, 3.25946658e+00, 6.00232275e-01,
  3.59386786e+00, 1.30000108e+00, 4.95670856e-01, 8.85913772e-01,
  1.08769141e+00, 4.96979579e-01, 1.02264352e+00, 4.72210752e-01,
  5.80058073e-01, 2.86239890e+00, 6.85801186e-01, 9.91862879e-01,
  8.29088221e-01, 1.05297301e+00, 9.85207016e-01, 5.20111147e-01,
  3.80841719e-01, 7.50000000e-01, 6.27051284e-01, 1.18227759e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.82061098e-01},
 {0.00000000e+00, 8.00397827e-01, 6.00232275e-01, 1.70900872e+01,
  5.32858023e-01, 5.45353890e-01, 6.01657157e-01, 5.56827599e-01,
  4.90637659e-01, 5.42801532e-01, 5.46735291e-01, 6.72663401e-01,
  6.03730225e-01, 6.80232381e-01, 4.11520916e-01, 4.86216592e-01,
  4.71814444e-01, 7.96941950e-01, 8.21766666e-01, 7.14771712e-01,
  3.33596922e-01, 7.50000000e-01, 4.88896599e-01, 5.22760621e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 6.21018471e-01},
 {0.00000000e+00, 6.88936672e-01, 3.59386786e+00, 5.32858023e-01,
  5.35580571e+00, 1.64256912e+00, 4.31018420e-01, 7.40433725e-01,
  9.76205497e-01, 4.40334375e-01, 9.41681214e-01, 4.62523949e-01,
  4.94422549e-01, 1.50174498e+00, 7.24041930e-01, 9.57667486e-01,
  7.71179788e-01, 9.29413105e-01, 8.76125496e-01, 4.93585385e-01,
  3.73386452e-01, 7.50000000e-01, 6.44927232e-01, 1.38090402e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.53699350e-01},
 {0.00000000e+00, 8.25164920e-01, 1.30000108e+00, 5.45353890e-01,
  1.64256912e+00, 3.87327599e+00, 4.97795679e-01, 5.76408577e-01,
  9.61718230e-01, 4.85270933e-01, 1.27686825e+00, 5.70784073e-01,
  6.14865688e-01, 8.93236200e-01, 9.10746443e-01, 1.53097375e+00,
  1.01074883e+00, 9.12113764e-01, 8.32873235e-01, 5.55160304e-01,
  5.19337483e-01, 7.50000000e-01, 6.16552718e-01, 2.97840479e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.36776244e-01},
 {0.00000000e+00, 5.87357723e-01, 4.95670856e-01, 6.01657157e-01,
  4.31018420e-01, 4.97795679e-01, 5.74817622e+00, 4.80346092e-01,
  6.79103477e-01, 1.06375667e+00, 5.29188869e-01, 1.30295859e+00,
  1.06255416e+00, 5.72439082e-01, 4.51176385e-01, 4.43558969e-01,
  5.89537777e-01, 6.09495403e-01, 7.16103747e-01, 9.52475503e-01,
  1.35493791e+00, 7.50000000e-01, 2.18456832e+00, 4.77074668e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.20783008e+00},
 {0.00000000e+00, 1.08031132e+00, 8.85913772e-01, 5.56827599e-01,
  7.40433725e-01, 5.76408577e-01, 4.80346092e-01, 5.07068525e+00,
  6.62087167e-01, 4.16341047e-01, 6.77879412e-01, 4.50091586e-01,
  5.84692214e-01, 1.05865660e+00, 7.02165561e-01, 6.87007231e-01,
  5.70228460e-01, 1.05800656e+00, 6.92819008e-01, 4.79214700e-01,
  5.91296285e-01, 7.50000000e-01, 5.49197180e-01, 6.18662539e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.36669292e-01},
 {0.00000000e+00, 6.54086288e-01, 1.08769141e+00, 4.90637659e-01,
  9.76205497e-01, 9.61718230e-01, 6.79103477e-01, 6.62087167e-01,
  9.51252809e+00, 4.53313059e-01, 8.90272071e-01, 6.69868446e-01,
  9.18088604e-01, 1.22006964e+00, 6.61223470e-01, 1.15049417e+00,
  9.73045615e-01, 8.54331847e-01, 7.06245757e-01, 4.56693295e-01,
  4.51816356e-01, 7.50000000e-01, 1.47204221e+00, 1.03383965e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.83746261e-01},
 {0.00000000e+00, 7.46806187e-01, 4.96979579e-01, 5.42801532e-01,
  4.40334375e-01, 4.85270933e-01, 1.06375667e+00, 4.16341047e-01,
  4.53313059e-01, 3.23256769e+00, 5.32316397e-01, 1.59618413e+00,
  1.45527106e+00, 5.64240025e-01, 6.09639867e-01, 5.77938325e-01,
  4.88387978e-01, 6.18410187e-01, 8.47505386e-01, 2.17596400e+00,
  5.64907506e-01, 7.50000000e-01, 9.06458192e-01, 5.20674297e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.24695958e+00},
 {0.00000000e+00, 7.86209397e-01, 1.02264352e+00, 5.46735291e-01,
  9.41681214e-01, 1.27686825e+00, 5.29188869e-01, 6.77879412e-01,
  8.90272071e-01, 5.32316397e-01, 3.32707189e+00, 5.53563636e-01,
  7.37955763e-01, 1.11877807e+00, 7.81202774e-01, 1.33004839e+00,
  1.94261316e+00, 8.89937552e-01, 8.84562104e-01, 5.91651856e-01,
  5.61572955e-01, 7.50000000e-01, 7.37107274e-01, 1.29718560e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.45113795e-01},
 {0.00000000e+00, 7.12370041e-01, 4.72210752e-01, 6.72663401e-01,
  4.62523949e-01, 5.70784073e-01, 1.30295859e+00, 4.50091586e-01,
  6.69868446e-01, 1.59618413e+00, 5.53563636e-01, 2.99708655e+00,
  1.73144954e+00, 4.83712850e-01, 4.77913692e-01, 6.42028706e-01,
  6.01135200e-01, 5.55659969e-01, 7.80723755e-01, 1.33363845e+00,
  6.70858407e-01, 7.50000000e-01, 9.65090110e-01, 5.98002922e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.43995990e+00},
 {0.00000000e+00, 8.21348665e-01, 5.80058073e-01, 6.03730225e-01,
  4.94422549e-01, 6.14865688e-01, 1.06255416e+00, 5.84692214e-01,
  9.18088604e-01, 1.45527106e+00, 7.37955763e-01, 1.73144954e+00,
  4.11411354e+00, 6.81741591e-01, 6.43682874e-01, 9.40467390e-01,
  7.75906233e-01, 6.60370266e-01, 8.60449567e-01, 1.23582796e+00,
  6.34345311e-01, 7.50000000e-01, 1.02316322e+00, 7.39261071e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.62161577e+00},
 {0.00000000e+00, 7.89043130e-01, 2.86239890e+00, 6.80232381e-01,
  1.50174498e+00, 8.93236200e-01, 5.72439082e-01, 1.05865660e+00,
  1.22006964e+00, 5.64240025e-01, 1.11877807e+00, 4.83712850e-01,
  6.81741591e-01, 4.47803773e+00, 6.40394172e-01, 1.03246645e+00,
  8.97848625e-01, 1.19968790e+00, 1.11473028e+00, 5.51607804e-01,
  3.89694095e-01, 7.50000000e-01, 6.05825405e-01, 9.46428796e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.15737804e-01},
 {0.00000000e+00, 7.08569419e-01, 6.85801186e-01, 4.11520916e-01,
  7.24041930e-01, 9.10746443e-01, 4.51176385e-01, 7.02165561e-01,
  6.61223470e-01, 6.09639867e-01, 7.81202774e-01, 4.77913692e-01,
  6.43682874e-01, 6.40394172e-01, 8.81911509e+00, 7.15515810e-01,
  5.81631739e-01, 7.49733904e-01, 8.56242933e-01, 5.40037335e-01,
  5.25005050e-01, 7.50000000e-01, 4.78832406e-01, 8.36159027e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.30300041e-01},
 {0.00000000e+00, 8.66678731e-01, 9.91862879e-01, 4.86216592e-01,
  9.57667486e-01, 1.53097375e+00, 4.43558969e-01, 6.87007231e-01,
  1.15049417e+00, 5.77938325e-01, 1.33004839e+00, 6.42028706e-01,
  9.40467390e-01, 1.03246645e+00, 7.15515810e-01, 4.40728842e+00,
  1.32912854e+00, 1.09183956e+00, 7.80601862e-01, 5.47266398e-01,
  6.45177884e-01, 7.50000000e-01, 8.29182983e-01, 2.62986317e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 6.16540522e-01},
 {0.00000000e+00, 6.99695540e-01, 8.29088221e-01, 4.71814444e-01,
  7.71179788e-01, 1.01074883e+00, 5.89537777e-01, 5.70228460e-01,
  9.73045615e-01, 4.88387978e-01, 1.94261316e+00, 6.01135200e-01,
  7.75906233e-01, 8.97848625e-01, 5.81631739e-01, 1.32912854e+00,
  4.74702063e+00, 7.99048209e-01, 7.15164318e-01, 5.77699501e-01,
  5.80165842e-01, 7.50000000e-01, 8.07446927e-01, 1.13238507e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.56296615e-01},
 {0.00000000e+00, 1.30031418e+00, 1.05297301e+00, 7.96941950e-01,
  9.29413105e-01, 9.12113764e-01, 6.09495403e-01, 1.05800656e+00,
  8.54331847e-01, 6.18410187e-01, 8.89937552e-01, 5.55659969e-01,
  6.60370266e-01, 1.19968790e+00, 7.49733904e-01, 1.09183956e+00,
  7.99048209e-01, 2.78188630e+00, 1.47248598e+00, 7.27836330e-01,
  4.28363793e-01, 7.50000000e-01, 7.05878947e-01, 9.80777594e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.80615182e-01},
 {0.00000000e+00, 1.00058530e+00, 9.85207016e-01, 8.21766666e-01,
  8.76125496e-01, 8.32873235e-01, 7.16103747e-01, 6.92819008e-01,
  7.06245757e-01, 8.47505386e-01, 8.84562104e-01, 7.80723755e-01,
  8.60449567e-01, 1.11473028e+00, 8.56242933e-01, 7.80601862e-01,
  7.15164318e-01, 1.47248598e+00, 3.13871529e+00, 1.04019697e+00,
  4.54128072e-01, 7.50000000e-01, 7.43457494e-01, 8.12903077e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 8.07282226e-01},
 {0.00000000e+00, 1.00992663e+00, 5.20111147e-01, 7.14771712e-01,
  4.93585385e-01, 5.55160304e-01, 9.52475503e-01, 4.79214700e-01,
  4.56693295e-01, 2.17596400e+00, 5.91651856e-01, 1.33363845e+00,
  1.23582796e+00, 5.51607804e-01, 5.40037335e-01, 5.47266398e-01,
  5.77699501e-01, 7.27836330e-01, 1.04019697e+00, 2.87075890e+00,
  4.73320057e-01, 7.50000000e-01, 8.09252575e-01, 5.52144455e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.66862396e+00},
 {0.00000000e+00, 5.65442334e-01, 3.80841719e-01, 3.33596922e-01,
  3.73386452e-01, 5.19337483e-01, 1.35493791e+00, 5.91296285e-01,
  4.51816356e-01, 5.64907506e-01, 5.61572955e-01, 6.70858407e-01,
  6.34345311e-01, 3.89694095e-01, 5.25005050e-01, 6.45177884e-01,
  5.80165842e-01, 4.28363793e-01, 4.54128072e-01, 4.73320057e-01,
  2.97023509e+01, 7.50000000e-01, 1.80096028e+00, 5.67414520e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 6.28722660e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 6.38727873e-01, 6.27051284e-01, 4.88896599e-01,
  6.44927232e-01, 6.16552718e-01, 2.18456832e+00, 5.49197180e-01,
  1.47204221e+00, 9.06458192e-01, 7.37107274e-01, 9.65090110e-01,
  1.02316322e+00, 6.05825405e-01, 4.78832406e-01, 8.29182983e-01,
  8.07446927e-01, 7.05878947e-01, 7.43457494e-01, 8.09252575e-01,
  1.80096028e+00, 7.50000000e-01, 5.75351902e+00, 6.97787623e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 9.41772708e-01},
 {0.00000000e+00, 8.41025176e-01, 1.18227759e+00, 5.22760621e-01,
  1.38090402e+00, 2.97840479e+00, 4.77074668e-01, 6.18662539e-01,
  1.03383965e+00, 5.20674297e-01, 1.29718560e+00, 5.98002922e-01,
  7.39261071e-01, 9.46428796e-01, 8.36159027e-01, 2.62986317e+00,
  1.13238507e+00, 9.80777594e-01, 8.12903077e-01, 5.52144455e-01,
  5.67414520e-01, 7.50000000e-01, 6.97787623e-01, 2.84524527e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.67250003e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 1.33300000e+00, 2.90000000e-01, 2.90000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 7.26064996e-01, 4.82061098e-01, 6.21018471e-01,
  4.53699350e-01, 5.36776244e-01, 1.20783008e+00, 4.36669292e-01,
  5.83746261e-01, 2.24695958e+00, 5.45113795e-01, 2.43995990e+00,
  1.62161577e+00, 5.15737804e-01, 5.30300041e-01, 6.16540522e-01,
  5.56296615e-01, 5.80615182e-01, 8.07282226e-01, 1.66862396e+00,
  6.28722660e-01, 7.50000000e-01, 9.41772708e-01, 5.67250003e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.36320536e+00}};


/** Underlying frequency ratios for BLOSUM50 */
static const double BLOSUM50_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 3.27354473e+00, 6.87168642e-01, 8.87624875e-01,
  6.59704230e-01, 7.97359654e-01, 5.46298170e-01, 1.10130683e+00,
  6.41220589e-01, 7.15157692e-01, 7.47622201e-01, 6.56954186e-01,
  8.53472686e-01, 7.19514908e-01, 7.14770771e-01, 8.19913708e-01,
  6.68351460e-01, 1.36359270e+00, 9.67331593e-01, 9.81625416e-01,
  4.63560176e-01, 7.50000000e-01, 5.96400452e-01, 8.05980349e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 6.80310043e-01},
 {0.00000000e+00, 6.87168642e-01, 3.67565906e+00, 5.07294595e-01,
  4.02052534e+00, 1.31892999e+00, 3.83849184e-01, 8.47577177e-01,
  1.11810068e+00, 4.07212168e-01, 9.68440788e-01, 3.91166888e-01,
  5.16058658e-01, 3.26949213e+00, 6.61247189e-01, 9.98204754e-01,
  7.59007679e-01, 1.06055276e+00, 9.55438618e-01, 4.46192804e-01,
  3.38571955e-01, 7.50000000e-01, 5.55905269e-01, 1.19634119e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 3.97605526e-01},
 {0.00000000e+00, 8.87624875e-01, 5.07294595e-01, 1.82308935e+01,
  4.27724382e-01, 4.56030174e-01, 5.64985221e-01, 5.24350848e-01,
  5.17412429e-01, 5.87186086e-01, 4.59864212e-01, 6.50074165e-01,
  6.80946766e-01, 6.01008569e-01, 4.03060607e-01, 4.81296027e-01,
  4.27834290e-01, 8.29850973e-01, 8.17890869e-01, 8.08665030e-01,
  3.12131245e-01, 7.50000000e-01, 5.38704705e-01, 4.65687383e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 6.24838486e-01},
 {0.00000000e+00, 6.59704230e-01, 4.02052534e+00, 4.27724382e-01,
  6.11147890e+00, 1.65784569e+00, 3.37276799e-01, 7.44468416e-01,
  8.83789762e-01, 3.70146565e-01, 8.90348134e-01, 3.72923686e-01,
  4.35009775e-01, 1.55790069e+00, 7.09628659e-01, 9.66997497e-01,
  6.78533894e-01, 9.21480865e-01, 8.25179634e-01, 4.31712984e-01,
  3.07050170e-01, 7.50000000e-01, 5.04473723e-01, 1.39378711e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 3.71809285e-01},
 {0.00000000e+00, 7.97359654e-01, 1.31892999e+00, 4.56030174e-01,
  1.65784569e+00, 4.43735647e+00, 4.56055858e-01, 5.39909105e-01,
  9.11487897e-01, 4.15558851e-01, 1.33108575e+00, 4.78428941e-01,
  6.05096970e-01, 9.19771370e-01, 8.38805682e-01, 1.67134451e+00,
  9.76357610e-01, 8.82100970e-01, 8.19031320e-01, 4.79075838e-01,
  5.48039829e-01, 7.50000000e-01, 6.51686488e-01, 3.38012103e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.53200481e-01},
 {0.00000000e+00, 5.46298170e-01, 3.83849184e-01, 5.64985221e-01,
  3.37276799e-01, 4.56055858e-01, 6.63625360e+00, 4.13949535e-01,
  7.54714659e-01, 9.89742646e-01, 4.44979578e-01, 1.26171801e+00,
  1.05158910e+00, 4.38699901e-01, 3.75736079e-01, 4.16967765e-01,
  4.61789222e-01, 5.52981315e-01, 5.70349305e-01, 8.59643167e-01,
  1.34807169e+00, 7.50000000e-01, 2.42442622e+00, 4.41115461e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.15257995e+00},
 {0.00000000e+00, 1.10130683e+00, 8.47577177e-01, 5.24350848e-01,
  7.44468416e-01, 5.39909105e-01, 4.13949535e-01, 5.79218671e+00,
  6.01271290e-01, 3.70370366e-01, 6.54870319e-01, 3.76903001e-01,
  5.19438610e-01, 9.69013722e-01, 6.20168128e-01, 6.41416113e-01,
  5.16985831e-01, 9.99248542e-01, 6.36269018e-01, 4.02720033e-01,
  5.04888636e-01, 7.50000000e-01, 4.67274430e-01, 5.78707493e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 3.74281590e-01},
 {0.00000000e+00, 6.41220589e-01, 1.11810068e+00, 5.17412429e-01,
  8.83789762e-01, 9.11487897e-01, 7.54714659e-01, 6.01271290e-01,
  1.04489376e+01, 4.11545408e-01, 9.45545516e-01, 5.47445792e-01,
  7.60124356e-01, 1.39406083e+00, 5.81906417e-01, 1.20911332e+00,
  9.81604707e-01, 8.22540644e-01, 6.52641826e-01, 4.13620259e-01,
  4.75002356e-01, 7.50000000e-01, 1.57000854e+00, 1.02524740e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.92911793e-01},
 {0.00000000e+00, 7.15157692e-01, 4.07212168e-01, 5.87186086e-01,
  3.70146565e-01, 4.15558851e-01, 9.89742646e-01, 3.70370366e-01,
  4.11545408e-01, 3.41093885e+00, 4.68297844e-01, 1.69677965e+00,
  1.43810563e+00, 4.50866254e-01, 5.11210841e-01, 5.02656404e-01,
  4.35318230e-01, 5.45643330e-01, 8.60722536e-01, 2.31272269e+00,
  5.22495607e-01, 7.50000000e-01, 8.26095669e-01, 4.48849604e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.38463611e+00},
 {0.00000000e+00, 7.47622201e-01, 9.68440788e-01, 4.59864212e-01,
  8.90348134e-01, 1.33108575e+00, 4.44979578e-01, 6.54870319e-01,
  9.45545516e-01, 4.68297844e-01, 3.88090096e+00, 4.79194854e-01,
  6.85231759e-01, 1.06041456e+00, 7.63679086e-01, 1.41855717e+00,
  2.06468049e+00, 8.92977250e-01, 8.44796364e-01, 5.22802406e-01,
  4.61593643e-01, 7.50000000e-01, 6.67754286e-01, 1.36451940e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.74822110e-01},
 {0.00000000e+00, 6.56954186e-01, 3.91166888e-01, 6.50074165e-01,
  3.72923686e-01, 4.78428941e-01, 1.26171801e+00, 3.76903001e-01,
  5.47445792e-01, 1.69677965e+00, 4.79194854e-01, 3.32815910e+00,
  1.78991633e+00, 4.12652855e-01, 4.43674459e-01, 5.62937315e-01,
  5.54029894e-01, 5.10605044e-01, 7.59171671e-01, 1.32427289e+00,
  6.01518374e-01, 7.50000000e-01, 8.58373419e-01, 5.10730048e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.67352044e+00},
 {0.00000000e+00, 8.53472686e-01, 5.16058658e-01, 6.80946766e-01,
  4.35009775e-01, 6.05096970e-01, 1.05158910e+00, 5.19438610e-01,
  7.60124356e-01, 1.43810563e+00, 6.85231759e-01, 1.78991633e+00,
  4.81561797e+00, 6.11514139e-01, 5.44813932e-01, 9.59915799e-01,
  6.78015272e-01, 6.82291911e-01, 8.82308943e-01, 1.21329921e+00,
  7.61522900e-01, 7.50000000e-01, 9.17831903e-01, 7.40717150e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.64874201e+00},
 {0.00000000e+00, 7.19514908e-01, 3.26949213e+00, 6.01008569e-01,
  1.55790069e+00, 9.19771370e-01, 4.38699901e-01, 9.69013722e-01,
  1.39406083e+00, 4.50866254e-01, 1.06041456e+00, 4.12652855e-01,
  6.11514139e-01, 5.28532229e+00, 6.04265819e-01, 1.03495916e+00,
  8.53785837e-01, 1.22434496e+00, 1.10885139e+00, 4.63246441e-01,
  3.75696800e-01, 7.50000000e-01, 6.16478872e-01, 9.63798879e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.27987098e-01},
 {0.00000000e+00, 7.14770771e-01, 6.61247189e-01, 4.03060607e-01,
  7.09628659e-01, 8.38805682e-01, 3.75736079e-01, 6.20168128e-01,
  5.81906417e-01, 5.11210841e-01, 7.63679086e-01, 4.43674459e-01,
  5.44813932e-01, 6.04265819e-01, 1.02035160e+01, 7.50499600e-01,
  5.18638945e-01, 7.59438841e-01, 7.48565360e-01, 5.24438149e-01,
  4.20092966e-01, 7.50000000e-01, 4.68173126e-01, 8.05053001e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.70775405e-01},
 {0.00000000e+00, 8.19913708e-01, 9.98204754e-01, 4.81296027e-01,
  9.66997497e-01, 1.67134451e+00, 4.16967765e-01, 6.41416113e-01,
  1.20911332e+00, 5.02656404e-01, 1.41855717e+00, 5.62937315e-01,
  9.59915799e-01, 1.03495916e+00, 7.50499600e-01, 4.69722165e+00,
  1.35733364e+00, 1.06872445e+00, 8.10316946e-01, 5.57384267e-01,
  7.14705591e-01, 7.50000000e-01, 7.42076535e-01, 2.82790659e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.38747839e-01},
 {0.00000000e+00, 6.68351460e-01, 7.59007679e-01, 4.27834290e-01,
  6.78533894e-01, 9.76357610e-01, 4.61789222e-01, 5.16985831e-01,
  9.81604707e-01, 4.35318230e-01, 2.06468049e+00, 5.54029894e-01,
  6.78015272e-01, 8.53785837e-01, 5.18638945e-01, 1.35733364e+00,
  5.37787401e+00, 8.04499038e-01, 7.36510915e-01, 5.12488758e-01,
  5.28823677e-01, 7.50000000e-01, 7.27125062e-01, 1.12197569e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.06393371e-01},
 {0.00000000e+00, 1.36359270e+00, 1.06055276e+00, 8.29850973e-01,
  9.21480865e-01, 8.82100970e-01, 5.52981315e-01, 9.99248542e-01,
  8.22540644e-01, 5.45643330e-01, 8.92977250e-01, 5.10605044e-01,
  6.82291911e-01, 1.22434496e+00, 7.59438841e-01, 1.06872445e+00,
  8.04499038e-01, 3.14298812e+00, 1.49727124e+00, 6.78664919e-01,
  3.92328021e-01, 7.50000000e-01, 6.51019697e-01, 9.53432893e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.24665180e-01},
 {0.00000000e+00, 9.67331593e-01, 9.55438618e-01, 8.17890869e-01,
  8.25179634e-01, 8.19031320e-01, 5.70349305e-01, 6.36269018e-01,
  6.52641826e-01, 8.60722536e-01, 8.44796364e-01, 7.59171671e-01,
  8.82308943e-01, 1.10885139e+00, 7.48565360e-01, 8.10316946e-01,
  7.36510915e-01, 1.49727124e+00, 3.55307500e+00, 1.05992520e+00,
  5.02669810e-01, 7.50000000e-01, 6.90730891e-01, 8.15700480e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.99921922e-01},
 {0.00000000e+00, 9.81625416e-01, 4.46192804e-01, 8.08665030e-01,
  4.31712984e-01, 4.79075838e-01, 8.59643167e-01, 4.02720033e-01,
  4.13620259e-01, 2.31272269e+00, 5.22802406e-01, 1.32427289e+00,
  1.21329921e+00, 4.63246441e-01, 5.24438149e-01, 5.57384267e-01,
  5.12488758e-01, 6.78664919e-01, 1.05992520e+00, 3.11745700e+00,
  4.84839541e-01, 7.50000000e-01, 7.27350506e-01, 5.09007179e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 1.72091725e+00},
 {0.00000000e+00, 4.63560176e-01, 3.38571955e-01, 3.12131245e-01,
  3.07050170e-01, 5.48039829e-01, 1.34807169e+00, 5.04888636e-01,
  4.75002356e-01, 5.22495607e-01, 4.61593643e-01, 6.01518374e-01,
  7.61522900e-01, 3.75696800e-01, 4.20092966e-01, 7.14705591e-01,
  5.28823677e-01, 3.92328021e-01, 5.02669810e-01, 4.84839541e-01,
  3.13609332e+01, 7.50000000e-01, 1.76515899e+00, 6.11743441e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 5.69808181e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 5.96400452e-01, 5.55905269e-01, 5.38704705e-01,
  5.04473723e-01, 6.51686488e-01, 2.42442622e+00, 4.67274430e-01,
  1.57000854e+00, 8.26095669e-01, 6.67754286e-01, 8.58373419e-01,
  9.17831903e-01, 6.16478872e-01, 4.68173126e-01, 7.42076535e-01,
  7.27125062e-01, 6.51019697e-01, 6.90730891e-01, 7.27350506e-01,
  1.76515899e+00, 7.50000000e-01, 6.89283261e+00, 6.86235710e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 8.45421029e-01},
 {0.00000000e+00, 8.05980349e-01, 1.19634119e+00, 4.65687383e-01,
  1.39378711e+00, 3.38012103e+00, 4.41115461e-01, 5.78707493e-01,
  1.02524740e+00, 4.48849604e-01, 1.36451940e+00, 5.10730048e-01,
  7.40717150e-01, 9.63798879e-01, 8.05053001e-01, 2.82790659e+00,
  1.12197569e+00, 9.53432893e-01, 8.15700480e-01, 5.09007179e-01,
  6.11743441e-01, 7.50000000e-01, 6.86235710e-01, 3.16905156e+00,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 4.85898712e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 2.90000000e-01, 2.90000000e-01, 2.90000000e-01,
  2.90000000e-01, 1.33300000e+00, 2.90000000e-01, 2.90000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 6.80310043e-01, 3.97605526e-01, 6.24838486e-01,
  3.71809285e-01, 4.53200481e-01, 1.15257995e+00, 3.74281590e-01,
  4.92911793e-01, 2.38463611e+00, 4.74822110e-01, 2.67352044e+00,
  1.64874201e+00, 4.27987098e-01, 4.70775405e-01, 5.38747839e-01,
  5.06393371e-01, 5.24665180e-01, 7.99921922e-01, 1.72091725e+00,
  5.69808181e-01, 7.50000000e-01, 8.45421029e-01, 4.85898712e-01,
  7.50000000e-01, 2.90000000e-01, 7.50000000e-01, 2.55759716e+00}};


/** Underlying frequency ratios for BLOSUM62 */
static const double BLOSUM62_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 3.90294070e+00, 5.64459671e-01, 8.67987664e-01,
  5.44605275e-01, 7.41264113e-01, 4.64893827e-01, 1.05686961e+00,
  5.69364849e-01, 6.32481035e-01, 7.75390239e-01, 6.01945975e-01,
  7.23150342e-01, 5.88307640e-01, 7.54121369e-01, 7.56803943e-01,
  6.12698600e-01, 1.47210399e+00, 9.84401956e-01, 9.36458396e-01,
  4.16548781e-01, 7.50000000e-01, 5.42611869e-01, 7.47274948e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.14377313e-01},
 {0.00000000e+00, 5.64459671e-01, 4.43758048e+00, 3.45226274e-01,
  4.74290926e+00, 1.33503378e+00, 3.24101420e-01, 7.38524318e-01,
  9.25449581e-01, 3.33981361e-01, 8.54849426e-01, 2.97257620e-01,
  4.04640322e-01, 4.07083696e+00, 5.53838329e-01, 9.44103648e-01,
  7.02873767e-01, 1.05798620e+00, 8.26250098e-01, 3.51280513e-01,
  2.52855433e-01, 7.50000000e-01, 4.09444638e-01, 1.18382127e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.12208474e-01},
 {0.00000000e+00, 8.67987664e-01, 3.45226274e-01, 1.95765857e+01,
  3.01454345e-01, 2.85934574e-01, 4.38990118e-01, 4.20387870e-01,
  3.55049505e-01, 6.53458801e-01, 3.49128465e-01, 6.42275633e-01,
  6.11354340e-01, 3.97802620e-01, 3.79562691e-01, 3.65781531e-01,
  3.08939296e-01, 7.38415701e-01, 7.40551692e-01, 7.55844055e-01,
  4.49983903e-01, 7.50000000e-01, 4.34203398e-01, 3.16819526e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.46828489e-01},
 {0.00000000e+00, 5.44605275e-01, 4.74290926e+00, 3.01454345e-01,
  7.39792738e+00, 1.68781075e+00, 2.98969081e-01, 6.34301019e-01,
  6.78558839e-01, 3.39015407e-01, 7.84090406e-01, 2.86613046e-01,
  3.46454634e-01, 1.55385281e+00, 5.98716826e-01, 8.97081129e-01,
  5.73200024e-01, 9.13504624e-01, 6.94789868e-01, 3.36500142e-01,
  2.32102315e-01, 7.50000000e-01, 3.45683565e-01, 1.38195506e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.07946931e-01},
 {0.00000000e+00, 7.41264113e-01, 1.33503378e+00, 2.85934574e-01,
  1.68781075e+00, 5.46952608e+00, 3.30743991e-01, 4.81267655e-01,
  9.60040718e-01, 3.30522558e-01, 1.30827885e+00, 3.72873704e-01,
  5.00342289e-01, 9.11298183e-01, 6.79202587e-01, 1.90173784e+00,
  9.60797602e-01, 9.50357185e-01, 7.41425610e-01, 4.28943130e-01,
  3.74300212e-01, 7.50000000e-01, 4.96467354e-01, 4.08949895e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.55631838e-01},
 {0.00000000e+00, 4.64893827e-01, 3.24101420e-01, 4.38990118e-01,
  2.98969081e-01, 3.30743991e-01, 8.12879702e+00, 3.40640908e-01,
  6.51990521e-01, 9.45769883e-01, 3.44043119e-01, 1.15459749e+00,
  1.00437163e+00, 3.54288952e-01, 2.87444758e-01, 3.33972402e-01,
  3.80726330e-01, 4.39973597e-01, 4.81693683e-01, 7.45089738e-01,
  1.37437942e+00, 7.50000000e-01, 2.76938063e+00, 3.31992746e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.06958025e+00},
 {0.00000000e+00, 1.05686961e+00, 7.38524318e-01, 4.20387870e-01,
  6.34301019e-01, 4.81267655e-01, 3.40640908e-01, 6.87630691e+00,
  4.92966576e-01, 2.75009722e-01, 5.88871736e-01, 2.84504012e-01,
  3.95486600e-01, 8.63711406e-01, 4.77385507e-01, 5.38649627e-01,
  4.49983999e-01, 9.03596525e-01, 5.79271582e-01, 3.36954912e-01,
  4.21690355e-01, 7.50000000e-01, 3.48714366e-01, 5.03463109e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.80638726e-01},
 {0.00000000e+00, 5.69364849e-01, 9.25449581e-01, 3.55049505e-01,
  6.78558839e-01, 9.60040718e-01, 6.51990521e-01, 4.92966576e-01,
  1.35059997e+01, 3.26288125e-01, 7.78887490e-01, 3.80675486e-01,
  5.84132623e-01, 1.22200067e+00, 4.72879831e-01, 1.16798104e+00,
  9.17048021e-01, 7.36731740e-01, 5.57503254e-01, 3.39447442e-01,
  4.44088955e-01, 7.50000000e-01, 1.79790413e+00, 1.04047242e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.58533474e-01},
 {0.00000000e+00, 6.32481035e-01, 3.33981361e-01, 6.53458801e-01,
  3.39015407e-01, 3.30522558e-01, 9.45769883e-01, 2.75009722e-01,
  3.26288125e-01, 3.99792994e+00, 3.96372934e-01, 1.69443475e+00,
  1.47774450e+00, 3.27934752e-01, 3.84662860e-01, 3.82937802e-01,
  3.54751311e-01, 4.43163582e-01, 7.79816110e-01, 2.41751209e+00,
  4.08874390e-01, 7.50000000e-01, 6.30388931e-01, 3.50796872e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.63222650e+00},
 {0.00000000e+00, 7.75390239e-01, 8.54849426e-01, 3.49128465e-01,
  7.84090406e-01, 1.30827885e+00, 3.44043119e-01, 5.88871736e-01,
  7.78887490e-01, 3.96372934e-01, 4.76433717e+00, 4.28270363e-01,
  6.25302816e-01, 9.39841129e-01, 7.03774479e-01, 1.55432308e+00,
  2.07680867e+00, 9.31919141e-01, 7.92905803e-01, 4.56542720e-01,
  3.58930071e-01, 7.50000000e-01, 5.32179333e-01, 1.40344922e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.15284382e-01},
 {0.00000000e+00, 6.01945975e-01, 2.97257620e-01, 6.42275633e-01,
  2.86613046e-01, 3.72873704e-01, 1.15459749e+00, 2.84504012e-01,
  3.80675486e-01, 1.69443475e+00, 4.28270363e-01, 3.79662137e+00,
  1.99429557e+00, 3.10043276e-01, 3.71121724e-01, 4.77325586e-01,
  4.73919278e-01, 4.28893743e-01, 6.60328975e-01, 1.31423573e+00,
  5.68037074e-01, 7.50000000e-01, 6.92059423e-01, 4.13275887e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.94078574e+00},
 {0.00000000e+00, 7.23150342e-01, 4.04640322e-01, 6.11354340e-01,
  3.46454634e-01, 5.00342289e-01, 1.00437163e+00, 3.95486600e-01,
  5.84132623e-01, 1.47774450e+00, 6.25302816e-01, 1.99429557e+00,
  6.48145121e+00, 4.74529655e-01, 4.23898024e-01, 8.64250293e-01,
  6.22623369e-01, 5.98558924e-01, 7.93801616e-01, 1.26893679e+00,
  6.10296214e-01, 7.50000000e-01, 7.08364628e-01, 6.41102583e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.78399892e+00},
 {0.00000000e+00, 5.88307640e-01, 4.07083696e+00, 3.97802620e-01,
  1.55385281e+00, 9.11298183e-01, 3.54288952e-01, 8.63711406e-01,
  1.22200067e+00, 3.27934752e-01, 9.39841129e-01, 3.10043276e-01,
  4.74529655e-01, 7.09409488e+00, 4.99932836e-01, 1.00058442e+00,
  8.58630478e-01, 1.23152924e+00, 9.84152635e-01, 3.69033853e-01,
  2.77782896e-01, 7.50000000e-01, 4.86030806e-01, 9.45834265e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.17327197e-01},
 {0.00000000e+00, 7.54121369e-01, 5.53838329e-01, 3.79562691e-01,
  5.98716826e-01, 6.79202587e-01, 2.87444758e-01, 4.77385507e-01,
  4.72879831e-01, 3.84662860e-01, 7.03774479e-01, 3.71121724e-01,
  4.23898024e-01, 4.99932836e-01, 1.28375437e+01, 6.41280589e-01,
  4.81534905e-01, 7.55503259e-01, 6.88897122e-01, 4.43082984e-01,
  2.81833164e-01, 7.50000000e-01, 3.63521119e-01, 6.64534287e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.76634549e-01},
 {0.00000000e+00, 7.56803943e-01, 9.44103648e-01, 3.65781531e-01,
  8.97081129e-01, 1.90173784e+00, 3.33972402e-01, 5.38649627e-01,
  1.16798104e+00, 3.82937802e-01, 1.55432308e+00, 4.77325586e-01,
  8.64250293e-01, 1.00058442e+00, 6.41280589e-01, 6.24442175e+00,
  1.40579606e+00, 9.65555228e-01, 7.91320741e-01, 4.66777931e-01,
  5.09360272e-01, 7.50000000e-01, 6.11094097e-01, 3.58149606e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.38898727e-01},
 {0.00000000e+00, 6.12698600e-01, 7.02873767e-01, 3.08939296e-01,
  5.73200024e-01, 9.60797602e-01, 3.80726330e-01, 4.49983999e-01,
  9.17048021e-01, 3.54751311e-01, 2.07680867e+00, 4.73919278e-01,
  6.22623369e-01, 8.58630478e-01, 4.81534905e-01, 1.40579606e+00,
  6.66557707e+00, 7.67165633e-01, 6.77754679e-01, 4.20072316e-01,
  3.95102106e-01, 7.50000000e-01, 5.55965425e-01, 1.13292384e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.25403989e-01},
 {0.00000000e+00, 1.47210399e+00, 1.05798620e+00, 7.38415701e-01,
  9.13504624e-01, 9.50357185e-01, 4.39973597e-01, 9.03596525e-01,
  7.36731740e-01, 4.43163582e-01, 9.31919141e-01, 4.28893743e-01,
  5.98558924e-01, 1.23152924e+00, 7.55503259e-01, 9.65555228e-01,
  7.67165633e-01, 3.84284741e+00, 1.61392097e+00, 5.65223766e-01,
  3.85303035e-01, 7.50000000e-01, 5.57520051e-01, 9.56235816e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.34703235e-01},
 {0.00000000e+00, 9.84401956e-01, 8.26250098e-01, 7.40551692e-01,
  6.94789868e-01, 7.41425610e-01, 4.81693683e-01, 5.79271582e-01,
  5.57503254e-01, 7.79816110e-01, 7.92905803e-01, 6.60328975e-01,
  7.93801616e-01, 9.84152635e-01, 6.88897122e-01, 7.91320741e-01,
  6.77754679e-01, 1.61392097e+00, 4.83210516e+00, 9.80943005e-01,
  4.30934144e-01, 7.50000000e-01, 5.73156574e-01, 7.60725140e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.08974203e-01},
 {0.00000000e+00, 9.36458396e-01, 3.51280513e-01, 7.55844055e-01,
  3.36500142e-01, 4.28943130e-01, 7.45089738e-01, 3.36954912e-01,
  3.39447442e-01, 2.41751209e+00, 4.56542720e-01, 1.31423573e+00,
  1.26893679e+00, 3.69033853e-01, 4.43082984e-01, 4.66777931e-01,
  4.20072316e-01, 5.65223766e-01, 9.80943005e-01, 3.69215640e+00,
  3.74456332e-01, 7.50000000e-01, 6.58038693e-01, 4.43577702e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.76339815e+00},
 {0.00000000e+00, 4.16548781e-01, 2.52855433e-01, 4.49983903e-01,
  2.32102315e-01, 3.74300212e-01, 1.37437942e+00, 4.21690355e-01,
  4.44088955e-01, 4.08874390e-01, 3.58930071e-01, 5.68037074e-01,
  6.10296214e-01, 2.77782896e-01, 2.81833164e-01, 5.09360272e-01,
  3.95102106e-01, 3.85303035e-01, 4.30934144e-01, 3.74456332e-01,
  3.81077833e+01, 7.50000000e-01, 2.10980812e+00, 4.26541694e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 5.03239261e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 5.42611869e-01, 4.09444638e-01, 4.34203398e-01,
  3.45683565e-01, 4.96467354e-01, 2.76938063e+00, 3.48714366e-01,
  1.79790413e+00, 6.30388931e-01, 5.32179333e-01, 6.92059423e-01,
  7.08364628e-01, 4.86030806e-01, 3.63521119e-01, 6.11094097e-01,
  5.55965425e-01, 5.57520051e-01, 5.73156574e-01, 6.58038693e-01,
  2.10980812e+00, 7.50000000e-01, 9.83220341e+00, 5.40805192e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.66952325e-01},
 {0.00000000e+00, 7.47274948e-01, 1.18382127e+00, 3.16819526e-01,
  1.38195506e+00, 4.08949895e+00, 3.31992746e-01, 5.03463109e-01,
  1.04047242e+00, 3.50796872e-01, 1.40344922e+00, 4.13275887e-01,
  6.41102583e-01, 9.45834265e-01, 6.64534287e-01, 3.58149606e+00,
  1.13292384e+00, 9.56235816e-01, 7.60725140e-01, 4.43577702e-01,
  4.26541694e-01, 7.50000000e-01, 5.40805192e-01, 3.89300249e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.87839626e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 1.33300000e+00, 2.50000000e-01, 2.50000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 6.14377313e-01, 3.12208474e-01, 6.46828489e-01,
  3.07946931e-01, 3.55631838e-01, 1.06958025e+00, 2.80638726e-01,
  3.58533474e-01, 2.63222650e+00, 4.15284382e-01, 2.94078574e+00,
  1.78399892e+00, 3.17327197e-01, 3.76634549e-01, 4.38898727e-01,
  4.25403989e-01, 4.34703235e-01, 7.08974203e-01, 1.76339815e+00,
  5.03239261e-01, 7.50000000e-01, 6.66952325e-01, 3.87839626e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.81516607e+00}};


/** Underlying frequency ratios for BLOSUM80 */
static const double BLOSUM80_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 4.77313697e+00, 4.77250219e-01, 7.31790935e-01,
  4.50905428e-01, 7.03141318e-01, 3.96805472e-01, 9.57108118e-01,
  5.14428496e-01, 5.43302985e-01, 7.22895459e-01, 5.04576848e-01,
  6.25509899e-01, 5.10109254e-01, 7.71095354e-01, 6.95763501e-01,
  5.55324683e-01, 1.53463504e+00, 9.79634362e-01, 8.65874830e-01,
  3.09281009e-01, 7.50000000e-01, 4.36234494e-01, 7.00327585e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 5.20181427e-01},
 {0.00000000e+00, 4.77250219e-01, 5.36249708e+00, 2.52383415e-01,
  5.75921529e+00, 1.26868736e+00, 2.51452613e-01, 6.38679494e-01,
  8.30167261e-01, 2.33516284e-01, 7.93249084e-01, 2.20539708e-01,
  3.09639510e-01, 4.86768283e+00, 4.28138134e-01, 8.50001386e-01,
  6.08526572e-01, 9.48166786e-01, 7.43109762e-01, 2.67977796e-01,
  1.78830543e-01, 7.50000000e-01, 3.07187296e-01, 1.10900996e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.25768580e-01},
 {0.00000000e+00, 7.31790935e-01, 2.52383415e-01, 2.07020176e+01,
  2.13899571e-01, 1.79625278e-01, 3.94756499e-01, 2.72055621e-01,
  2.21322329e-01, 5.81174860e-01, 2.41148280e-01, 4.92657548e-01,
  4.99414247e-01, 3.00383113e-01, 2.69170214e-01, 2.95117493e-01,
  2.32476995e-01, 5.76114817e-01, 6.01733396e-01, 6.34451252e-01,
  3.02232769e-01, 7.50000000e-01, 3.07930831e-01, 2.23671408e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 5.28325329e-01},
 {0.00000000e+00, 4.50905428e-01, 5.75921529e+00, 2.13899571e-01,
  9.10634223e+00, 1.63552603e+00, 2.34434106e-01, 5.40654444e-01,
  5.94384672e-01, 2.14047157e-01, 6.76875328e-01, 1.97166761e-01,
  2.51946380e-01, 1.58444827e+00, 4.52407652e-01, 7.63275386e-01,
  4.76639455e-01, 7.74041684e-01, 6.10882957e-01, 2.44693975e-01,
  1.44774364e-01, 7.50000000e-01, 2.45193294e-01, 1.30286928e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.03968665e-01},
 {0.00000000e+00, 7.03141318e-01, 1.26868736e+00, 1.79625278e-01,
  1.63552603e+00, 6.99475019e+00, 2.48970812e-01, 3.99359878e-01,
  9.01259830e-01, 2.64258683e-01, 1.19492933e+00, 2.75527698e-01,
  4.28576433e-01, 8.11140936e-01, 5.81414067e-01, 1.90643780e+00,
  8.31863127e-01, 8.45039731e-01, 6.85529376e-01, 3.69379335e-01,
  2.40807712e-01, 7.50000000e-01, 3.32874109e-01, 5.05418243e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.70986883e-01},
 {0.00000000e+00, 3.96805472e-01, 2.51452613e-01, 3.94756499e-01,
  2.34434106e-01, 2.48970812e-01, 9.48547379e+00, 2.48796540e-01,
  5.71854758e-01, 8.40573731e-01, 2.82718808e-01, 1.11429351e+00,
  8.93330223e-01, 2.72679266e-01, 2.37160811e-01, 2.85158291e-01,
  2.87269045e-01, 3.69100895e-01, 4.45199900e-01, 6.48888163e-01,
  1.08939386e+00, 7.50000000e-01, 2.78024963e+00, 2.62771902e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 1.00399896e+00},
 {0.00000000e+00, 9.57108118e-01, 6.38679494e-01, 2.72055621e-01,
  5.40654444e-01, 3.99359878e-01, 2.48796540e-01, 7.88244144e+00,
  3.87151686e-01, 1.84455102e-01, 4.83435068e-01, 2.10542976e-01,
  2.86067872e-01, 7.60943081e-01, 3.47300820e-01, 4.24585221e-01,
  3.77442786e-01, 7.84315177e-01, 4.92238430e-01, 2.50954947e-01,
  2.64396616e-01, 7.50000000e-01, 2.29516029e-01, 4.08980256e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.00030947e-01},
 {0.00000000e+00, 5.14428496e-01, 8.30167261e-01, 2.21322329e-01,
  5.94384672e-01, 9.01259830e-01, 5.71854758e-01, 3.87151686e-01,
  1.60694674e+01, 2.57933576e-01, 7.40110055e-01, 3.14299881e-01,
  4.32076355e-01, 1.12425153e+00, 4.19656882e-01, 1.31555625e+00,
  9.25536475e-01, 6.61407005e-01, 5.39570743e-01, 2.88933760e-01,
  3.90225420e-01, 7.50000000e-01, 1.81930455e+00, 1.05926315e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.91587251e-01},
 {0.00000000e+00, 5.43302985e-01, 2.33516284e-01, 5.81174860e-01,
  2.14047157e-01, 2.64258683e-01, 8.40573731e-01, 1.84455102e-01,
  2.57933576e-01, 4.86762150e+00, 3.13345237e-01, 1.66499837e+00,
  1.51247384e+00, 2.57799519e-01, 2.85790430e-01, 3.09071252e-01,
  2.99348100e-01, 3.78995471e-01, 7.00511896e-01, 2.49584558e+00,
  3.43150987e-01, 7.50000000e-01, 5.39308441e-01, 2.81349188e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.95548558e+00},
 {0.00000000e+00, 7.22895459e-01, 7.93249084e-01, 2.41148280e-01,
  6.76875328e-01, 1.19492933e+00, 2.82718808e-01, 4.83435068e-01,
  7.40110055e-01, 3.13345237e-01, 6.32564527e+00, 3.56851811e-01,
  5.34403407e-01, 9.38398440e-01, 5.96694133e-01, 1.52403165e+00,
  2.19214139e+00, 8.20193974e-01, 7.35729790e-01, 3.70194033e-01,
  2.41427653e-01, 7.50000000e-01, 4.08501941e-01, 1.32044155e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.39320970e-01},
 {0.00000000e+00, 5.04576848e-01, 2.20539708e-01, 4.92657548e-01,
  1.97166761e-01, 2.75527698e-01, 1.11429351e+00, 2.10542976e-01,
  3.14299881e-01, 1.66499837e+00, 3.56851811e-01, 4.46305621e+00,
  2.12274889e+00, 2.49692056e-01, 3.03099202e-01, 4.06904090e-01,
  3.62830591e-01, 3.68478686e-01, 5.60836408e-01, 1.22050154e+00,
  4.38789464e-01, 7.50000000e-01, 5.80503535e-01, 3.25631695e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.33558735e+00},
 {0.00000000e+00, 6.25509899e-01, 3.09639510e-01, 4.99414247e-01,
  2.51946380e-01, 4.28576433e-01, 8.93330223e-01, 2.86067872e-01,
  4.32076355e-01, 1.51247384e+00, 5.34403407e-01, 2.12274889e+00,
  8.88346290e+00, 3.81598352e-01, 3.61925345e-01, 8.86557630e-01,
  5.05866341e-01, 4.98438721e-01, 7.57959723e-01, 1.22414515e+00,
  5.60653516e-01, 7.50000000e-01, 5.49937808e-01, 6.03240148e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 1.87684042e+00},
 {0.00000000e+00, 5.10109254e-01, 4.86768283e+00, 3.00383113e-01,
  1.58444827e+00, 8.11140936e-01, 2.72679266e-01, 7.60943081e-01,
  1.12425153e+00, 2.57799519e-01, 9.38398440e-01, 2.49692056e-01,
  3.81598352e-01, 8.96275887e+00, 3.97867522e-01, 9.58172021e-01,
  7.73025258e-01, 1.16534759e+00, 9.08032132e-01, 2.97018980e-01,
  2.21307752e-01, 7.50000000e-01, 3.84510481e-01, 8.67215281e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.52958933e-01},
 {0.00000000e+00, 7.71095354e-01, 4.28138134e-01, 2.69170214e-01,
  4.52407652e-01, 5.81414067e-01, 2.37160811e-01, 3.47300820e-01,
  4.19656882e-01, 2.85790430e-01, 5.96694133e-01, 3.03099202e-01,
  3.61925345e-01, 3.97867522e-01, 1.51545798e+01, 5.37994539e-01,
  4.45589871e-01, 6.51704179e-01, 5.60357280e-01, 3.70233083e-01,
  1.78033069e-01, 7.50000000e-01, 2.57545983e-01, 5.64854837e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 2.96124685e-01},
 {0.00000000e+00, 6.95763501e-01, 8.50001386e-01, 2.95117493e-01,
  7.63275386e-01, 1.90643780e+00, 2.85158291e-01, 4.24585221e-01,
  1.31555625e+00, 3.09071252e-01, 1.52403165e+00, 4.06904090e-01,
  8.86557630e-01, 9.58172021e-01, 5.37994539e-01, 8.33990474e+00,
  1.39424540e+00, 8.58988713e-01, 7.24369298e-01, 4.10943414e-01,
  4.07901262e-01, 7.50000000e-01, 4.61857147e-01, 4.36001721e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.67482647e-01},
 {0.00000000e+00, 5.55324683e-01, 6.08526572e-01, 2.32476995e-01,
  4.76639455e-01, 8.31863127e-01, 2.87269045e-01, 3.77442786e-01,
  9.25536475e-01, 2.99348100e-01, 2.19214139e+00, 3.62830591e-01,
  5.05866341e-01, 7.73025258e-01, 4.45589871e-01, 1.39424540e+00,
  8.24459589e+00, 6.94509540e-01, 5.98385216e-01, 3.53719047e-01,
  2.94245493e-01, 7.50000000e-01, 4.17775411e-01, 1.04634306e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.37250515e-01},
 {0.00000000e+00, 1.53463504e+00, 9.48166786e-01, 5.76114817e-01,
  7.74041684e-01, 8.45039731e-01, 3.69100895e-01, 7.84315177e-01,
  6.61407005e-01, 3.78995471e-01, 8.20193974e-01, 3.68478686e-01,
  4.98438721e-01, 1.16534759e+00, 6.51704179e-01, 8.58988713e-01,
  6.94509540e-01, 5.10577131e+00, 1.66260189e+00, 4.93679246e-01,
  2.70773669e-01, 7.50000000e-01, 4.62005069e-01, 8.50359559e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.72716393e-01},
 {0.00000000e+00, 9.79634362e-01, 7.43109762e-01, 6.01733396e-01,
  6.10882957e-01, 6.85529376e-01, 4.45199900e-01, 4.92238430e-01,
  5.39570743e-01, 7.00511896e-01, 7.35729790e-01, 5.60836408e-01,
  7.57959723e-01, 9.08032132e-01, 5.60357280e-01, 7.24369298e-01,
  5.98385216e-01, 1.66260189e+00, 6.20547751e+00, 8.91492247e-01,
  2.85084781e-01, 7.50000000e-01, 4.74451641e-01, 7.00342047e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 6.17118219e-01},
 {0.00000000e+00, 8.65874830e-01, 2.67977796e-01, 6.34451252e-01,
  2.44693975e-01, 3.69379335e-01, 6.48888163e-01, 2.50954947e-01,
  2.88933760e-01, 2.49584558e+00, 3.70194033e-01, 1.22050154e+00,
  1.22414515e+00, 2.97018980e-01, 3.70233083e-01, 4.10943414e-01,
  3.53719047e-01, 4.93679246e-01, 8.91492247e-01, 4.58356287e+00,
  3.42175943e-01, 7.50000000e-01, 4.89223745e-01, 3.85230939e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 1.73439753e+00},
 {0.00000000e+00, 3.09281009e-01, 1.78830543e-01, 3.02232769e-01,
  1.44774364e-01, 2.40807712e-01, 1.08939386e+00, 2.64396616e-01,
  3.90225420e-01, 3.43150987e-01, 2.41427653e-01, 4.38789464e-01,
  5.60653516e-01, 2.21307752e-01, 1.78033069e-01, 4.07901262e-01,
  2.94245493e-01, 2.70773669e-01, 2.85084781e-01, 3.42175943e-01,
  4.15522183e+01, 7.50000000e-01, 2.03605072e+00, 3.04533429e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 4.00252232e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 4.36234494e-01, 3.07187296e-01, 3.07930831e-01,
  2.45193294e-01, 3.32874109e-01, 2.78024963e+00, 2.29516029e-01,
  1.81930455e+00, 5.39308441e-01, 4.08501941e-01, 5.80503535e-01,
  5.49937808e-01, 3.84510481e-01, 2.57545983e-01, 4.61857147e-01,
  4.17775411e-01, 4.62005069e-01, 4.74451641e-01, 4.89223745e-01,
  2.03605072e+00, 7.50000000e-01, 1.21940332e+01, 3.82065335e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 5.63904098e-01},
 {0.00000000e+00, 7.00327585e-01, 1.10900996e+00, 2.23671408e-01,
  1.30286928e+00, 5.05418243e+00, 2.62771902e-01, 4.08980256e-01,
  1.05926315e+00, 2.81349188e-01, 1.32044155e+00, 3.25631695e-01,
  6.03240148e-01, 8.67215281e-01, 5.64854837e-01, 4.36001721e+00,
  1.04634306e+00, 8.50359559e-01, 7.00342047e-01, 3.85230939e-01,
  3.04533429e-01, 7.50000000e-01, 3.82065335e-01, 4.78944345e+00,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.07788194e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.40000000e-01, 1.40000000e-01, 1.40000000e-01,
  1.40000000e-01, 1.33300000e+00, 1.40000000e-01, 1.40000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 5.20181427e-01, 2.25768580e-01, 5.28325329e-01,
  2.03968665e-01, 2.70986883e-01, 1.00399896e+00, 2.00030947e-01,
  2.91587251e-01, 2.95548558e+00, 3.39320970e-01, 3.33558735e+00,
  1.87684042e+00, 2.52958933e-01, 2.96124685e-01, 3.67482647e-01,
  3.37250515e-01, 3.72716393e-01, 6.17118219e-01, 1.73439753e+00,
  4.00252232e-01, 7.50000000e-01, 5.63904098e-01, 3.07788194e-01,
  7.50000000e-01, 1.40000000e-01, 7.50000000e-01, 3.18242650e+00}};


/** Underlying frequency ratios for BLOSUM90 */
static const double BLOSUM90_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 5.49812903e+00, 4.38682591e-01, 6.69393543e-01,
  4.10593787e-01, 6.67362415e-01, 3.46907686e-01, 8.55253386e-01,
  4.61993747e-01, 4.75838329e-01, 6.64089300e-01, 4.46304927e-01,
  5.48948177e-01, 4.75350260e-01, 7.12075051e-01, 6.37043670e-01,
  5.31189855e-01, 1.45279278e+00, 9.34706358e-01, 7.93712275e-01,
  2.83637106e-01, 7.50000000e-01, 3.70274337e-01, 6.55892238e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 4.58271706e-01},
 {0.00000000e+00, 4.38682591e-01, 5.85593083e+00, 2.24322103e-01,
  6.24911512e+00, 1.23816514e+00, 2.14608760e-01, 5.64122407e-01,
  7.82997018e-01, 1.96621553e-01, 7.38810616e-01, 1.94127812e-01,
  2.69937418e-01, 5.34266045e+00, 3.78549747e-01, 7.81147613e-01,
  5.45762328e-01, 8.70271567e-01, 6.66528685e-01, 2.28904046e-01,
  1.42810925e-01, 7.50000000e-01, 2.77062678e-01, 1.06526643e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 1.95138263e-01},
 {0.00000000e+00, 6.69393543e-01, 2.24322103e-01, 2.16370056e+01,
  1.76054992e-01, 1.47039970e-01, 3.92579579e-01, 2.18182122e-01,
  1.93552170e-01, 5.20896357e-01, 2.15679622e-01, 4.25126435e-01,
  4.27104871e-01, 2.87330927e-01, 2.30151009e-01, 2.42068607e-01,
  2.06242089e-01, 5.23756471e-01, 5.48342086e-01, 5.88014943e-01,
  2.50594014e-01, 7.50000000e-01, 2.92512133e-01, 1.82991170e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 4.63931905e-01},
 {0.00000000e+00, 4.10593787e-01, 6.24911512e+00, 1.76054992e-01,
  9.87094836e+00, 1.61092322e+00, 1.97472144e-01, 4.84699282e-01,
  5.46208942e-01, 1.71646922e-01, 6.20942905e-01, 1.71021465e-01,
  2.15937410e-01, 1.52110375e+00, 3.99848977e-01, 6.86343892e-01,
  4.17577647e-01, 6.96691332e-01, 5.29001708e-01, 2.06199101e-01,
  1.15609648e-01, 7.50000000e-01, 2.13739536e-01, 1.26113669e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 1.71274897e-01},
 {0.00000000e+00, 6.67362415e-01, 1.23816514e+00, 1.47039970e-01,
  1.61092322e+00, 7.91072600e+00, 2.04920324e-01, 3.47008914e-01,
  7.72534981e-01, 2.37473636e-01, 1.07737256e+00, 2.37238651e-01,
  3.61519510e-01, 7.51559518e-01, 5.28989676e-01, 1.83144568e+00,
  7.31811567e-01, 7.59912485e-01, 6.13331345e-01, 3.33942314e-01,
  2.02669731e-01, 7.50000000e-01, 2.73296982e-01, 5.61081481e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 2.37333866e-01},
 {0.00000000e+00, 3.46907686e-01, 2.14608760e-01, 3.92579579e-01,
  1.97472144e-01, 2.04920324e-01, 1.05190688e+01, 1.97046381e-01,
  5.10867985e-01, 7.73236776e-01, 2.73932585e-01, 1.05058955e+00,
  8.11125795e-01, 2.36979230e-01, 2.12448013e-01, 2.80742738e-01,
  2.62210098e-01, 3.33306129e-01, 3.90547453e-01, 5.77272331e-01,
  9.82154073e-01, 7.50000000e-01, 2.77172979e+00, 2.33605433e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 9.38207661e-01},
 {0.00000000e+00, 8.55253386e-01, 5.64122407e-01, 2.18182122e-01,
  4.84699282e-01, 3.47008914e-01, 1.97046381e-01, 8.29576319e+00,
  3.35169505e-01, 1.51428808e-01, 4.27685370e-01, 1.83265898e-01,
  2.37360208e-01, 6.67802895e-01, 2.98989408e-01, 3.80871171e-01,
  3.28608745e-01, 7.00460225e-01, 4.00474342e-01, 2.02269779e-01,
  2.34021877e-01, 7.50000000e-01, 1.88904853e-01, 3.59819671e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 1.70365676e-01},
 {0.00000000e+00, 4.61993747e-01, 7.82997018e-01, 1.93552170e-01,
  5.46208942e-01, 7.72534981e-01, 5.10867985e-01, 3.35169505e-01,
  1.85636930e+01, 2.32236241e-01, 6.80230209e-01, 2.89296785e-01,
  4.06641319e-01, 1.09210477e+00, 3.99528365e-01, 1.24778653e+00,
  8.70382542e-01, 5.94270003e-01, 4.71855284e-01, 2.49881715e-01,
  3.84396440e-01, 7.50000000e-01, 1.62620965e+00, 9.52331980e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 2.66176152e-01},
 {0.00000000e+00, 4.75838329e-01, 1.96621553e-01, 5.20896357e-01,
  1.71646922e-01, 2.37473636e-01, 7.73236776e-01, 1.51428808e-01,
  2.32236241e-01, 5.62471003e+00, 2.84042165e-01, 1.59032509e+00,
  1.47982983e+00, 2.29223921e-01, 2.60766373e-01, 2.80218761e-01,
  2.69951432e-01, 3.32350183e-01, 6.20404233e-01, 2.47466358e+00,
  2.84494454e-01, 7.50000000e-01, 4.78737209e-01, 2.53644956e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.22503669e+00},
 {0.00000000e+00, 6.64089300e-01, 7.38810616e-01, 2.15679622e-01,
  6.20942905e-01, 1.07737256e+00, 2.73932585e-01, 4.27685370e-01,
  6.80230209e-01, 2.84042165e-01, 7.40971353e+00, 3.18645229e-01,
  5.03980900e-01, 8.92677413e-01, 5.54085013e-01, 1.48494120e+00,
  2.15988586e+00, 7.31921781e-01, 6.80549543e-01, 3.28415139e-01,
  1.74946435e-01, 7.50000000e-01, 3.46172886e-01, 1.23156378e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.04624249e-01},
 {0.00000000e+00, 4.46304927e-01, 1.94127812e-01, 4.25126435e-01,
  1.71021465e-01, 2.37238651e-01, 1.05058955e+00, 1.83265898e-01,
  2.89296785e-01, 1.59032509e+00, 3.18645229e-01, 5.03145573e+00,
  2.07702392e+00, 2.24291285e-01, 2.72188366e-01, 3.74709469e-01,
  3.26170372e-01, 3.38914311e-01, 5.07385737e-01, 1.13293282e+00,
  3.72830750e-01, 7.50000000e-01, 5.08794972e-01, 2.89246563e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.63712766e+00},
 {0.00000000e+00, 5.48948177e-01, 2.69937418e-01, 4.27104871e-01,
  2.15937410e-01, 3.61519510e-01, 8.11125795e-01, 2.37360208e-01,
  4.06641319e-01, 1.47982983e+00, 5.03980900e-01, 2.07702392e+00,
  1.13394084e+01, 3.40430076e-01, 3.16108504e-01, 8.49621832e-01,
  4.63059321e-01, 4.39100157e-01, 6.79030953e-01, 1.13727027e+00,
  4.92730381e-01, 7.50000000e-01, 4.80316866e-01, 5.46178209e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 1.83504401e+00},
 {0.00000000e+00, 4.75350260e-01, 5.34266045e+00, 2.87330927e-01,
  1.52110375e+00, 7.51559518e-01, 2.36979230e-01, 6.67802895e-01,
  1.09210477e+00, 2.29223921e-01, 8.92677413e-01, 2.24291285e-01,
  3.40430076e-01, 1.03313947e+01, 3.50745318e-01, 9.04906229e-01,
  7.13097098e-01, 1.09686657e+00, 8.46059068e-01, 2.58543521e-01,
  1.78320001e-01, 7.50000000e-01, 3.59725936e-01, 8.09573592e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 2.26289963e-01},
 {0.00000000e+00, 7.12075051e-01, 3.78549747e-01, 2.30151009e-01,
  3.99848977e-01, 5.28989676e-01, 2.12448013e-01, 2.98989408e-01,
  3.99528365e-01, 2.60766373e-01, 5.54085013e-01, 2.72188366e-01,
  3.16108504e-01, 3.50745318e-01, 1.60877707e+01, 4.86643270e-01,
  4.09306303e-01, 5.93677872e-01, 4.80576271e-01, 3.23032809e-01,
  1.74477792e-01, 7.50000000e-01, 2.13074361e-01, 5.12969199e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 2.67560235e-01},
 {0.00000000e+00, 6.37043670e-01, 7.81147613e-01, 2.42068607e-01,
  6.86343892e-01, 1.83144568e+00, 2.80742738e-01, 3.80871171e-01,
  1.24778653e+00, 2.80218761e-01, 1.48494120e+00, 3.74709469e-01,
  8.49621832e-01, 9.04906229e-01, 4.86643270e-01, 9.98642121e+00,
  1.25897253e+00, 7.86784913e-01, 6.71486014e-01, 3.48766905e-01,
  3.51106744e-01, 7.50000000e-01, 4.01714365e-01, 4.91663315e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.36422330e-01},
 {0.00000000e+00, 5.31189855e-01, 5.45762328e-01, 2.06242089e-01,
  4.17577647e-01, 7.31811567e-01, 2.62210098e-01, 3.28608745e-01,
  8.70382542e-01, 2.69951432e-01, 2.15988586e+00, 3.26170372e-01,
  4.63059321e-01, 7.13097098e-01, 4.09306303e-01, 1.25897253e+00,
  9.45518393e+00, 6.27094667e-01, 5.25362150e-01, 3.16415957e-01,
  2.48889195e-01, 7.50000000e-01, 3.62347797e-01, 9.31246918e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.03390753e-01},
 {0.00000000e+00, 1.45279278e+00, 8.70271567e-01, 5.23756471e-01,
  6.96691332e-01, 7.59912485e-01, 3.33306129e-01, 7.00460225e-01,
  5.94270003e-01, 3.32350183e-01, 7.31921781e-01, 3.38914311e-01,
  4.39100157e-01, 1.09686657e+00, 5.93677872e-01, 7.86784913e-01,
  6.27094667e-01, 5.99164849e+00, 1.62868240e+00, 4.49399804e-01,
  2.30316670e-01, 7.50000000e-01, 4.11583504e-01, 7.70078852e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.36254561e-01},
 {0.00000000e+00, 9.34706358e-01, 6.66528685e-01, 5.48342086e-01,
  5.29001708e-01, 6.13331345e-01, 3.90547453e-01, 4.00474342e-01,
  4.71855284e-01, 6.20404233e-01, 6.80549543e-01, 5.07385737e-01,
  6.79030953e-01, 8.46059068e-01, 4.80576271e-01, 6.71486014e-01,
  5.25362150e-01, 1.62868240e+00, 7.25868130e+00, 8.04058132e-01,
  2.27083637e-01, 7.50000000e-01, 4.20544390e-01, 6.35332398e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 5.53180238e-01},
 {0.00000000e+00, 7.93712275e-01, 2.28904046e-01, 5.88014943e-01,
  2.06199101e-01, 3.33942314e-01, 5.77272331e-01, 2.02269779e-01,
  2.49881715e-01, 2.47466358e+00, 3.28415139e-01, 1.13293282e+00,
  1.13727027e+00, 2.58543521e-01, 3.23032809e-01, 3.48766905e-01,
  3.16415957e-01, 4.49399804e-01, 8.04058132e-01, 5.31607752e+00,
  2.99756311e-01, 7.50000000e-01, 4.04621372e-01, 3.39550748e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 1.67659508e+00},
 {0.00000000e+00, 2.83637106e-01, 1.42810925e-01, 2.50594014e-01,
  1.15609648e-01, 2.02669731e-01, 9.82154073e-01, 2.34021877e-01,
  3.84396440e-01, 2.84494454e-01, 1.74946435e-01, 3.72830750e-01,
  4.92730381e-01, 1.78320001e-01, 1.74477792e-01, 3.51106744e-01,
  2.48889195e-01, 2.30316670e-01, 2.27083637e-01, 2.99756311e-01,
  4.26750567e+01, 7.50000000e-01, 1.77446682e+00, 2.58826370e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.37037347e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 3.70274337e-01, 2.77062678e-01, 2.92512133e-01,
  2.13739536e-01, 2.73296982e-01, 2.77172979e+00, 1.88904853e-01,
  1.62620965e+00, 4.78737209e-01, 3.46172886e-01, 5.08794972e-01,
  4.80316866e-01, 3.59725936e-01, 2.13074361e-01, 4.01714365e-01,
  3.62347797e-01, 4.11583504e-01, 4.20544390e-01, 4.04621372e-01,
  1.77446682e+00, 7.50000000e-01, 1.36090374e+01, 3.21879801e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 4.96615724e-01},
 {0.00000000e+00, 6.55892238e-01, 1.06526643e+00, 1.82991170e-01,
  1.26113669e+00, 5.61081481e+00, 2.33605433e-01, 3.59819671e-01,
  9.52331980e-01, 2.53644956e-01, 1.23156378e+00, 2.89246563e-01,
  5.46178209e-01, 8.09573592e-01, 5.12969199e-01, 4.91663315e+00,
  9.31246918e-01, 7.70078852e-01, 6.35332398e-01, 3.39550748e-01,
  2.58826370e-01, 7.50000000e-01, 3.21879801e-01, 5.34819225e+00,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 2.74820979e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.20000000e-01, 1.20000000e-01, 1.20000000e-01,
  1.20000000e-01, 1.33300000e+00, 1.20000000e-01, 1.20000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 4.58271706e-01, 1.95138263e-01, 4.63931905e-01,
  1.71274897e-01, 2.37333866e-01, 9.38207661e-01, 1.70365676e-01,
  2.66176152e-01, 3.22503669e+00, 3.04624249e-01, 3.63712766e+00,
  1.83504401e+00, 2.26289963e-01, 2.67560235e-01, 3.36422330e-01,
  3.03390753e-01, 3.36254561e-01, 5.53180238e-01, 1.67659508e+00,
  3.37037347e-01, 7.50000000e-01, 4.96615724e-01, 2.74820979e-01,
  7.50000000e-01, 1.20000000e-01, 7.50000000e-01, 3.47015056e+00}};


/** Underlying frequency ratios for PAM30 */
static const double PAM30_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 7.78912912e+00, 3.02026381e-01, 1.07971742e-01,
  3.16807704e-01, 4.53014491e-01, 5.67768414e-02, 5.75744366e-01,
  8.33329143e-02, 1.98631007e-01, 9.54324770e-02, 1.15281887e-01,
  1.89338526e-01, 2.84890839e-01, 5.93402017e-01, 2.35344432e-01,
  9.10275603e-02, 8.75340369e-01, 8.27084092e-01, 4.77504036e-01,
  9.76867708e-03, 7.50000000e-01, 6.99222384e-02, 3.58158006e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.40431646e-01},
 {0.00000000e+00, 3.02026381e-01, 8.11766307e+00, 1.52757717e-02,
  8.45582330e+00, 1.47233532e+00, 2.69459344e-02, 3.37867839e-01,
  6.63755159e-01, 1.25827363e-01, 4.40914175e-01, 4.89221062e-02,
  3.39137767e-02, 7.72564409e+00, 1.02509984e-01, 3.67813095e-01,
  8.23817770e-02, 6.03773710e-01, 3.53305100e-01, 6.71234492e-02,
  3.23783774e-02, 7.50000000e-01, 1.13742904e-01, 9.91005473e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.21274943e-02},
 {0.00000000e+00, 1.07971742e-01, 1.52757717e-02, 2.75882699e+01,
  7.71413760e-03, 7.86819951e-03, 1.16056875e-02, 4.15734061e-02,
  7.76695875e-02, 1.19064797e-01, 7.78302876e-03, 6.02160155e-03,
  9.46008155e-03, 2.40417467e-02, 6.32617140e-02, 8.13944094e-03,
  6.80179741e-02, 3.79110970e-01, 6.77518036e-02, 1.32369559e-01,
  4.67559620e-03, 7.50000000e-01, 2.57579768e-01, 7.98640138e-03,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 4.01312493e-02},
 {0.00000000e+00, 3.16807704e-01, 8.45582330e+00, 7.71413760e-03,
  1.42349183e+01, 2.34246282e+00, 6.53544516e-03, 3.26331184e-01,
  2.71287596e-01, 7.98669838e-02, 2.19798242e-01, 1.38483017e-02,
  2.29035955e-02, 1.75629226e+00, 6.80928461e-02, 4.27289142e-01,
  2.95026156e-02, 2.86012439e-01, 1.99709663e-01, 6.47860871e-02,
  5.42804139e-03, 7.50000000e-01, 1.97046510e-02, 1.50786644e+00,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 3.37687752e-02},
 {0.00000000e+00, 4.53014491e-01, 1.47233532e+00, 7.86819951e-03,
  2.34246282e+00, 1.36624533e+01, 7.86638136e-03, 2.31589040e-01,
  1.81505594e-01, 1.49736909e-01, 2.26119978e-01, 4.27891596e-02,
  9.24578171e-02, 4.63622671e-01, 1.53529752e-01, 1.50713628e+00,
  4.28627549e-02, 2.25474962e-01, 1.30340902e-01, 1.09854687e-01,
  3.23568839e-03, 7.50000000e-01, 5.65446816e-02, 8.36539662e+00,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.50595678e-02},
 {0.00000000e+00, 5.67768414e-02, 2.69459344e-02, 1.16056875e-02,
  6.53544516e-03, 7.86638136e-03, 2.14068164e+01, 4.44753676e-02,
  1.31584637e-01, 4.73011553e-01, 8.55473177e-03, 4.15676862e-01,
  2.49857403e-01, 5.06072009e-02, 3.56716332e-02, 1.12390799e-02,
  4.18628794e-02, 1.13897780e-01, 5.14073961e-02, 6.66400959e-02,
  2.13219283e-01, 7.50000000e-01, 1.78617980e+00, 9.33613946e-03,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 4.32977029e-01},
 {0.00000000e+00, 5.75744366e-01, 3.37867839e-01, 4.15734061e-02,
  3.26331184e-01, 2.31589040e-01, 4.44753676e-02, 9.31368857e+00,
  4.78871557e-02, 2.41395424e-02, 8.59403011e-02, 2.78096704e-02,
  5.67940501e-02, 3.51241937e-01, 1.21951553e-01, 9.44408638e-02,
  3.97493058e-02, 5.52895766e-01, 1.35075690e-01, 1.48994667e-01,
  5.68401675e-03, 7.50000000e-01, 9.25861867e-03, 1.71822465e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 2.67022462e-02},
 {0.00000000e+00, 8.33329143e-02, 6.63755159e-01, 7.76695875e-02,
  2.71287596e-01, 1.81505594e-01, 1.31584637e-01, 4.78871557e-02,
  2.29276353e+01, 4.36863126e-02, 1.19631443e-01, 1.23251186e-01,
  2.78547216e-02, 1.11873100e+00, 2.47982049e-01, 1.35964891e+00,
  5.86952720e-01, 1.31501535e-01, 8.61653901e-02, 1.14978024e-01,
  8.12035191e-02, 7.50000000e-01, 3.23746871e-01, 6.94918113e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 9.92432853e-02},
 {0.00000000e+00, 1.98631007e-01, 1.25827363e-01, 1.19064797e-01,
  7.98669838e-02, 1.49736909e-01, 4.73011553e-01, 2.41395424e-02,
  4.36863126e-02, 1.86346064e+01, 1.27210180e-01, 6.43409762e-01,
  7.92524183e-01, 1.79107849e-01, 5.32543187e-02, 7.03383168e-02,
  1.58533228e-01, 9.68805230e-02, 4.44355269e-01, 1.93788392e+00,
  9.23811792e-03, 7.50000000e-01, 1.18573178e-01, 1.15136508e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 6.07207247e+00},
 {0.00000000e+00, 9.54324770e-02, 4.40914175e-01, 7.78302876e-03,
  2.19798242e-01, 2.26119978e-01, 8.55473177e-03, 8.59403011e-02,
  1.19631443e-01, 1.27210180e-01, 9.98780825e+00, 6.07412260e-02,
  5.62560724e-01, 6.97247226e-01, 1.06485246e-01, 3.84266618e-01,
  1.11121735e+00, 2.59857530e-01, 3.36316040e-01, 4.62631009e-02,
  1.75142268e-02, 7.50000000e-01, 4.16645387e-02, 2.95037286e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 8.07975645e-02},
 {0.00000000e+00, 1.15281887e-01, 4.89221062e-02, 6.02160155e-03,
  1.38483017e-02, 4.27891596e-02, 4.15676862e-01, 2.78096704e-02,
  1.23251186e-01, 6.43409762e-01, 6.07412260e-02, 1.00194105e+01,
  1.24218990e+00, 8.95821132e-02, 9.04970355e-02, 1.84988987e-01,
  5.31753716e-02, 5.83025512e-02, 1.02884887e-01, 4.62744892e-01,
  1.27405990e-01, 7.50000000e-01, 9.50388211e-02, 1.04757148e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.19029657e+00},
 {0.00000000e+00, 1.89338526e-01, 3.39137767e-02, 9.46008155e-03,
  2.29035955e-02, 9.24578171e-02, 2.49857403e-01, 5.67940501e-02,
  2.78547216e-02, 7.92524183e-01, 5.62560724e-01, 1.24218990e+00,
  4.65953079e+01, 4.66775487e-02, 6.48026667e-02, 2.71732304e-01,
  2.43002286e-01, 1.60899996e-01, 2.62669188e-01, 6.30878126e-01,
  1.29243603e-02, 7.50000000e-01, 2.12719960e-02, 1.70582240e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.10650779e+00},
 {0.00000000e+00, 2.84890839e-01, 7.72564409e+00, 2.40417467e-02,
  1.75629226e+00, 4.63622671e-01, 5.06072009e-02, 3.51241937e-01,
  1.11873100e+00, 1.79107849e-01, 6.97247226e-01, 8.95821132e-02,
  4.66775487e-02, 1.46457342e+01, 1.42408737e-01, 2.98864303e-01,
  1.43682999e-01, 9.72144796e-01, 5.31363673e-01, 6.98330827e-02,
  6.36210916e-02, 7.50000000e-01, 2.22758622e-01, 3.91824098e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.16595604e-01},
 {0.00000000e+00, 5.93402017e-01, 1.02509984e-01, 6.32617140e-02,
  6.80928461e-02, 1.53529752e-01, 3.56716332e-02, 1.21951553e-01,
  2.47982049e-01, 5.32543187e-02, 1.06485246e-01, 9.04970355e-02,
  6.48026667e-02, 1.42408737e-01, 1.58088136e+01, 3.76066261e-01,
  2.69813309e-01, 5.70573606e-01, 2.40194503e-01, 1.40453450e-01,
  9.02230244e-03, 7.50000000e-01, 9.80486654e-03, 2.50506944e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.92594202e-02},
 {0.00000000e+00, 2.35344432e-01, 3.67813095e-01, 8.13944094e-03,
  4.27289142e-01, 1.50713628e+00, 1.12390799e-02, 9.44408638e-02,
  1.35964891e+00, 7.03383168e-02, 3.84266618e-01, 1.84988987e-01,
  2.71732304e-01, 2.98864303e-01, 3.76066261e-01, 1.81380436e+01,
  5.84920778e-01, 1.67743101e-01, 1.50331260e-01, 1.02742037e-01,
  1.20726953e-02, 7.50000000e-01, 1.69239569e-02, 8.75457039e+00,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.50394300e-01},
 {0.00000000e+00, 9.10275603e-02, 8.23817770e-02, 6.80179741e-02,
  2.95026156e-02, 4.28627549e-02, 4.18628794e-02, 3.97493058e-02,
  5.86952720e-01, 1.58533228e-01, 1.11121735e+00, 5.31753716e-02,
  2.43002286e-01, 1.43682999e-01, 2.69813309e-01, 5.84920778e-01,
  1.89240175e+01, 3.54610075e-01, 1.06409795e-01, 7.42726505e-02,
  5.17051285e-01, 7.50000000e-01, 3.01584225e-02, 2.79081365e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 8.49660454e-02},
 {0.00000000e+00, 8.75340369e-01, 6.03773710e-01, 3.79110970e-01,
  2.86012439e-01, 2.25474962e-01, 1.13897780e-01, 5.52895766e-01,
  1.31501535e-01, 9.68805230e-02, 2.59857530e-01, 5.83025512e-02,
  1.60899996e-01, 9.72144796e-01, 5.70573606e-01, 1.67743101e-01,
  3.54610075e-01, 9.02800848e+00, 1.14463515e+00, 1.17807053e-01,
  1.79210510e-01, 7.50000000e-01, 9.65415148e-02, 2.00316511e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 6.99430663e-02},
 {0.00000000e+00, 8.27084092e-01, 3.53305100e-01, 6.77518036e-02,
  1.99709663e-01, 1.30340902e-01, 5.14073961e-02, 1.35075690e-01,
  8.61653901e-02, 4.44355269e-01, 3.36316040e-01, 1.02884887e-01,
  2.62669188e-01, 5.31363673e-01, 2.40194503e-01, 1.50331260e-01,
  1.06409795e-01, 1.14463515e+00, 1.16950751e+01, 3.90871695e-01,
  1.25703039e-02, 7.50000000e-01, 1.11995399e-01, 1.39052321e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 2.05920141e-01},
 {0.00000000e+00, 4.77504036e-01, 6.71234492e-02, 1.32369559e-01,
  6.47860871e-02, 1.09854687e-01, 6.66400959e-02, 1.48994667e-01,
  1.14978024e-01, 1.93788392e+00, 4.62631009e-02, 4.62744892e-01,
  6.30878126e-01, 6.98330827e-02, 1.40453450e-01, 1.02742037e-01,
  7.42726505e-02, 1.17807053e-01, 3.90871695e-01, 1.16090589e+01,
  5.05672115e-03, 7.50000000e-01, 8.04763342e-02, 1.06755129e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 9.07853260e-01},
 {0.00000000e+00, 9.76867708e-03, 3.23783774e-02, 4.67559620e-03,
  5.42804139e-03, 3.23568839e-03, 2.13219283e-01, 5.68401675e-03,
  8.12035191e-02, 9.23811792e-03, 1.75142268e-02, 1.27405990e-01,
  1.29243603e-02, 6.36210916e-02, 9.02230244e-03, 1.20726953e-02,
  5.17051285e-01, 1.79210510e-01, 1.25703039e-02, 5.05672115e-03,
  8.86903633e+01, 7.50000000e-01, 1.73378233e-01, 7.08668846e-03,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 9.17500230e-02},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 6.99222384e-02, 1.13742904e-01, 2.57579768e-01,
  1.97046510e-02, 5.65446816e-02, 1.78617980e+00, 9.25861867e-03,
  3.23746871e-01, 1.18573178e-01, 4.16645387e-02, 9.50388211e-02,
  2.12719960e-02, 2.22758622e-01, 9.80486654e-03, 1.69239569e-02,
  3.01584225e-02, 9.65415148e-02, 1.11995399e-01, 8.04763342e-02,
  1.73378233e-01, 7.50000000e-01, 2.84454162e+01, 3.92787209e-02,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.02140077e-01},
 {0.00000000e+00, 3.58158006e-01, 9.91005473e-01, 7.98640138e-03,
  1.50786644e+00, 8.36539662e+00, 9.33613946e-03, 1.71822465e-01,
  6.94918113e-01, 1.15136508e-01, 2.95037286e-01, 1.04757148e-01,
  1.70582240e-01, 3.91824098e-01, 2.50506944e-01, 8.75457039e+00,
  2.79081365e-01, 2.00316511e-01, 1.39052321e-01, 1.06755129e-01,
  7.08668846e-03, 7.50000000e-01, 3.92787209e-02, 8.53499117e+00,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 1.07889016e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 3.00000000e-03, 3.00000000e-03, 3.00000000e-03,
  3.00000000e-03, 1.33300000e+00, 3.00000000e-03, 3.00000000e-03},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 1.40431646e-01, 7.21274943e-02, 4.01312493e-02,
  3.37687752e-02, 7.50595678e-02, 4.32977029e-01, 2.67022462e-02,
  9.92432853e-02, 6.07207247e+00, 8.07975645e-02, 7.19029657e+00,
  1.10650779e+00, 1.16595604e-01, 7.92594202e-02, 1.50394300e-01,
  8.49660454e-02, 6.99430663e-02, 2.05920141e-01, 9.07853260e-01,
  9.17500230e-02, 7.50000000e-01, 1.02140077e-01, 1.07889016e-01,
  7.50000000e-01, 3.00000000e-03, 7.50000000e-01, 6.85288369e+00}};


/** Underlying frequency ratios for PAM70 */
static const double PAM70_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 4.89972946e+00, 6.04638890e-01, 2.42067963e-01,
  6.18433760e-01, 7.71213936e-01, 1.34833512e-01, 1.01498720e+00,
  2.21390981e-01, 4.34566961e-01, 2.49319770e-01, 2.49736666e-01,
  3.75829845e-01, 5.88646912e-01, 1.03045017e+00, 4.65602795e-01,
  2.26934456e-01, 1.34964157e+00, 1.32707297e+00, 8.08640731e-01,
  4.32618091e-02, 7.50000000e-01, 1.54222102e-01, 6.38034395e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 3.05507349e-01},
 {0.00000000e+00, 6.04638890e-01, 5.42142610e+00, 6.60540784e-02,
  5.88929014e+00, 2.24885858e+00, 7.88958522e-02, 6.48829213e-01,
  1.09845328e+00, 2.52789787e-01, 7.77906911e-01, 1.19715714e-01,
  1.39730820e-01, 4.87904539e+00, 2.73736473e-01, 7.90138593e-01,
  2.53425041e-01, 9.54703234e-01, 6.53688648e-01, 1.85063408e-01,
  7.80441148e-02, 7.50000000e-01, 2.23329473e-01, 1.61317606e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.59869478e-01},
 {0.00000000e+00, 2.42067963e-01, 6.60540784e-02, 2.48333286e+01,
  3.97130183e-02, 3.87120726e-02, 5.48382532e-02, 1.15679023e-01,
  1.64834856e-01, 2.46442788e-01, 3.83065039e-02, 3.01944614e-02,
  4.42379868e-02, 9.65904773e-02, 1.58823055e-01, 3.93861754e-02,
  1.51302811e-01, 6.65193229e-01, 1.86417892e-01, 2.74442219e-01,
  2.22251459e-02, 7.50000000e-01, 5.23320532e-01, 3.90058338e-02,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 9.54452197e-02},
 {0.00000000e+00, 6.18433760e-01, 5.88929014e+00, 3.97130183e-02,
  8.89965209e+00, 3.35659301e+00, 3.36340414e-02, 6.39387415e-01,
  5.99746161e-01, 1.93031482e-01, 4.83846284e-01, 6.35580693e-02,
  1.03502951e-01, 2.39946804e+00, 2.12686642e-01, 9.22524100e-01,
  1.38560979e-01, 5.99713481e-01, 4.42110174e-01, 1.72922381e-01,
  2.64171458e-02, 7.50000000e-01, 7.98817496e-02, 2.29587193e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.02625370e-01},
 {0.00000000e+00, 7.71213936e-01, 2.24885858e+00, 3.87120726e-02,
  3.35659301e+00, 8.65385421e+00, 3.82398528e-02, 4.91841782e-01,
  4.80150693e-01, 2.83623924e-01, 4.66181910e-01, 1.17487568e-01,
  2.12426533e-01, 9.64695360e-01, 3.47480663e-01, 2.25815498e+00,
  1.78914910e-01, 4.78627630e-01, 3.38114944e-01, 2.46972907e-01,
  1.89756038e-02, 7.50000000e-01, 1.19141830e-01, 5.86672974e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.67617543e-01},
 {0.00000000e+00, 1.34833512e-01, 7.88958522e-02, 5.48382532e-02,
  3.36340414e-02, 3.82398528e-02, 1.74546228e+01, 1.00136016e-01,
  2.76577280e-01, 8.45971180e-01, 4.12894272e-02, 8.23804865e-01,
  5.08442478e-01, 1.31366508e-01, 8.91727264e-02, 5.25046736e-02,
  9.93639027e-02, 2.15603624e-01, 1.37600314e-01, 2.07856659e-01,
  4.52242030e-01, 7.50000000e-01, 3.37849562e+00, 4.44561913e-02,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 8.30493328e-01},
 {0.00000000e+00, 1.01498720e+00, 6.48829213e-01, 1.15679023e-01,
  6.39387415e-01, 4.91841782e-01, 1.00136016e-01, 7.30864335e+00,
  1.45138751e-01, 1.07004308e-01, 2.09861286e-01, 8.09786598e-02,
  1.44276056e-01, 6.59774806e-01, 3.08399509e-01, 2.37177521e-01,
  1.22086776e-01, 9.64187863e-01, 3.70424893e-01, 3.12986176e-01,
  2.68996967e-02, 7.50000000e-01, 4.37757319e-02, 3.80863925e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 8.88316375e-02},
 {0.00000000e+00, 2.21390981e-01, 1.09845328e+00, 1.64834856e-01,
  5.99746161e-01, 4.80150693e-01, 2.76577280e-01, 1.45138751e-01,
  1.64223829e+01, 1.39752446e-01, 3.35085062e-01, 2.59072609e-01,
  1.22395235e-01, 1.67658945e+00, 4.92867985e-01, 2.18832010e+00,
  1.05779530e+00, 3.26721114e-01, 2.30266854e-01, 2.33270842e-01,
  1.86063506e-01, 7.50000000e-01, 6.11386658e-01, 1.22453853e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 2.23068949e-01},
 {0.00000000e+00, 4.34566961e-01, 2.52789787e-01, 2.46442788e-01,
  1.93031482e-01, 2.83623924e-01, 8.45971180e-01, 1.07004308e-01,
  1.39752446e-01, 1.17505537e+01, 2.66278255e-01, 1.19239531e+00,
  1.36908981e+00, 3.22065793e-01, 1.59577365e-01, 1.89299534e-01,
  3.01479599e-01, 2.53646367e-01, 7.61111058e-01, 3.00374358e+00,
  4.36969076e-02, 7.50000000e-01, 2.80489781e-01, 2.42519143e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 4.37821344e+00},
 {0.00000000e+00, 2.49319770e-01, 7.77906911e-01, 3.83065039e-02,
  4.83846284e-01, 4.66181910e-01, 4.12894272e-02, 2.09861286e-01,
  3.35085062e-01, 2.66278255e-01, 7.61026217e+00, 1.52754534e-01,
  9.47937549e-01, 1.11880255e+00, 2.56118662e-01, 7.21669881e-01,
  1.93666556e+00, 5.20459614e-01, 6.20196625e-01, 1.38845936e-01,
  7.91911207e-02, 7.50000000e-01, 1.00771644e-01, 5.77518724e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.87009175e-01},
 {0.00000000e+00, 2.49736666e-01, 1.19715714e-01, 3.01944614e-02,
  6.35580693e-02, 1.17487568e-01, 8.23804865e-01, 8.09786598e-02,
  2.59072609e-01, 1.19239531e+00, 1.52754534e-01, 8.22843955e+00,
  2.12939550e+00, 1.84817584e-01, 1.98778178e-01, 3.50784854e-01,
  1.37359422e-01, 1.48574310e-01, 2.40326213e-01, 9.08401895e-01,
  2.67168548e-01, 7.50000000e-01, 2.38980900e-01, 2.19154102e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 6.10538395e+00},
 {0.00000000e+00, 3.75829845e-01, 1.39730820e-01, 4.42379868e-02,
  1.03502951e-01, 2.12426533e-01, 5.08442478e-01, 1.44276056e-01,
  1.22395235e-01, 1.36908981e+00, 9.47937549e-01, 2.12939550e+00,
  2.85861616e+01, 1.81728699e-01, 1.75926195e-01, 4.86044449e-01,
  4.88080973e-01, 3.18976045e-01, 4.97701201e-01, 1.13284160e+00,
  5.97810639e-02, 7.50000000e-01, 9.73038916e-02, 3.31664034e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.89998090e+00},
 {0.00000000e+00, 5.88646912e-01, 4.87904539e+00, 9.65904773e-02,
  2.39946804e+00, 9.64695360e-01, 1.31366508e-01, 6.59774806e-01,
  1.67658945e+00, 3.22065793e-01, 1.11880255e+00, 1.84817584e-01,
  1.81728699e-01, 7.75354485e+00, 3.44509706e-01, 6.36668056e-01,
  3.86583495e-01, 1.36623218e+00, 8.98965212e-01, 1.99138136e-01,
  1.37893708e-01, 7.50000000e-01, 3.89624106e-01, 8.21747280e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 2.26230851e-01},
 {0.00000000e+00, 1.03045017e+00, 2.73736473e-01, 1.58823055e-01,
  2.12686642e-01, 3.47480663e-01, 8.91727264e-02, 3.08399509e-01,
  4.92867985e-01, 1.59577365e-01, 2.56118662e-01, 1.98778178e-01,
  1.75926195e-01, 3.44509706e-01, 1.18714088e+01, 6.84263791e-01,
  5.26778639e-01, 9.79291808e-01, 5.30644206e-01, 3.07279290e-01,
  4.18182164e-02, 7.50000000e-01, 4.66910129e-02, 4.94244365e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.86949727e-01},
 {0.00000000e+00, 4.65602795e-01, 7.90138593e-01, 3.93861754e-02,
  9.22524100e-01, 2.25815498e+00, 5.25046736e-02, 2.37177521e-01,
  2.18832010e+00, 1.89299534e-01, 7.21669881e-01, 3.50784854e-01,
  4.86044449e-01, 6.36668056e-01, 6.84263791e-01, 1.14706864e+01,
  1.02361055e+00, 3.73519396e-01, 3.29322200e-01, 2.32053434e-01,
  5.49721446e-02, 7.50000000e-01, 7.38645866e-02, 6.27280154e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 3.02058283e-01},
 {0.00000000e+00, 2.26934456e-01, 2.53425041e-01, 1.51302811e-01,
  1.38560979e-01, 1.78914910e-01, 9.93639027e-02, 1.22086776e-01,
  1.05779530e+00, 3.01479599e-01, 1.93666556e+00, 1.37359422e-01,
  4.88080973e-01, 3.86583495e-01, 5.26778639e-01, 1.02361055e+00,
  1.36584013e+01, 6.10010848e-01, 2.83658630e-01, 1.75829089e-01,
  9.90780192e-01, 7.50000000e-01, 8.59463892e-02, 5.47017256e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.86881035e-01},
 {0.00000000e+00, 1.34964157e+00, 9.54703234e-01, 6.65193229e-01,
  5.99713481e-01, 4.78627630e-01, 2.15603624e-01, 9.64187863e-01,
  3.26721114e-01, 2.53646367e-01, 5.20459614e-01, 1.48574310e-01,
  3.18976045e-01, 1.36623218e+00, 9.79291808e-01, 3.73519396e-01,
  6.10010848e-01, 5.20445832e+00, 1.69490837e+00, 3.03324606e-01,
  3.23868004e-01, 7.50000000e-01, 2.08076540e-01, 4.32823454e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.80278747e-01},
 {0.00000000e+00, 1.32707297e+00, 6.53688648e-01, 1.86417892e-01,
  4.42110174e-01, 3.38114944e-01, 1.37600314e-01, 3.70424893e-01,
  2.30266854e-01, 7.61111058e-01, 6.20196625e-01, 2.40326213e-01,
  4.97701201e-01, 8.98965212e-01, 5.30644206e-01, 3.29322200e-01,
  2.83658630e-01, 1.69490837e+00, 7.33707214e+00, 7.18561720e-01,
  5.45397612e-02, 7.50000000e-01, 2.27833757e-01, 3.34283233e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 3.97467804e-01},
 {0.00000000e+00, 8.08640731e-01, 1.85063408e-01, 2.74442219e-01,
  1.72922381e-01, 2.46972907e-01, 2.07856659e-01, 3.12986176e-01,
  2.33270842e-01, 3.00374358e+00, 1.38845936e-01, 9.08401895e-01,
  1.13284160e+00, 1.99138136e-01, 3.07279290e-01, 2.32053434e-01,
  1.75829089e-01, 3.03324606e-01, 7.18561720e-01, 8.21110545e+00,
  2.66588959e-02, 7.50000000e-01, 1.80843519e-01, 2.40471284e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.54065018e+00},
 {0.00000000e+00, 4.32618091e-02, 7.80441148e-02, 2.22251459e-02,
  2.64171458e-02, 1.89756038e-02, 4.52242030e-01, 2.68996967e-02,
  1.86063506e-01, 4.36969076e-02, 7.91911207e-02, 2.67168548e-01,
  5.97810639e-02, 1.37893708e-01, 4.18182164e-02, 5.49721446e-02,
  9.90780192e-01, 3.23868004e-01, 5.45397612e-02, 2.66588959e-02,
  8.06167129e+01, 7.50000000e-01, 3.76692625e-01, 3.46622138e-02,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 1.99738227e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 1.54222102e-01, 2.23329473e-01, 5.23320532e-01,
  7.98817496e-02, 1.19141830e-01, 3.37849562e+00, 4.37757319e-02,
  6.11386658e-01, 2.80489781e-01, 1.00771644e-01, 2.38980900e-01,
  9.73038916e-02, 3.89624106e-01, 4.66910129e-02, 7.38645866e-02,
  8.59463892e-02, 2.08076540e-01, 2.27833757e-01, 1.80843519e-01,
  3.76692625e-01, 7.50000000e-01, 2.31438270e+01, 9.94108657e-02,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 2.51505787e-01},
 {0.00000000e+00, 6.38034395e-01, 1.61317606e+00, 3.90058338e-02,
  2.29587193e+00, 5.86672974e+00, 4.44561913e-02, 3.80863925e-01,
  1.22453853e+00, 2.42519143e-01, 5.77518724e-01, 2.19154102e-01,
  3.31664034e-01, 8.21747280e-01, 4.94244365e-01, 6.27280154e+00,
  5.47017256e-01, 4.32823454e-01, 3.34283233e-01, 2.40471284e-01,
  3.46622138e-02, 7.50000000e-01, 9.94108657e-02, 6.04368813e+00,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 2.26204268e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 2.00000000e-02, 2.00000000e-02, 2.00000000e-02,
  2.00000000e-02, 1.33300000e+00, 2.00000000e-02, 2.00000000e-02},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 3.05507349e-01, 1.59869478e-01, 9.54452197e-02,
  1.02625370e-01, 1.67617543e-01, 8.30493328e-01, 8.88316375e-02,
  2.23068949e-01, 4.37821344e+00, 1.87009175e-01, 6.10538395e+00,
  1.89998090e+00, 2.26230851e-01, 1.86949727e-01, 3.02058283e-01,
  1.86881035e-01, 1.80278747e-01, 3.97467804e-01, 1.54065018e+00,
  1.99738227e-01, 7.50000000e-01, 2.51505787e-01, 2.26204268e-01,
  7.50000000e-01, 2.00000000e-02, 7.50000000e-01, 5.58422761e+00}};


/** Underlying frequency ratios for PAM250 */
static const double PAM250_FREQRATIOS[POSITAA_SIZE][POSITAA_SIZE] =
{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 1.51578006e+00, 1.05606466e+00, 6.26687137e-01,
  1.07002228e+00, 1.07474947e+00, 4.46290560e-01, 1.33923306e+00,
  7.32178470e-01, 8.88728816e-01, 7.66227814e-01, 6.46341735e-01,
  7.70359541e-01, 1.03988401e+00, 1.29408826e+00, 9.03456780e-01,
  7.00798809e-01, 1.28966428e+00, 1.31725380e+00, 1.04499203e+00,
  2.63609841e-01, 7.50000000e-01, 4.49686798e-01, 1.00010336e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 7.19479599e-01},
 {0.00000000e+00, 1.05606466e+00, 1.84257846e+00, 3.65113184e-01,
  2.05195422e+00, 1.82470118e+00, 3.54071966e-01, 1.11918465e+00,
  1.29683470e+00, 6.16970770e-01, 1.13010166e+00, 4.52559218e-01,
  6.04775077e-01, 1.59985542e+00, 8.46142888e-01, 1.33890747e+00,
  8.61808201e-01, 1.11706882e+00, 1.03197349e+00, 6.37811920e-01,
  2.96793416e-01, 7.50000000e-01, 4.85608019e-01, 1.61300149e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.02168751e-01},
 {0.00000000e+00, 6.26687137e-01, 3.65113184e-01, 1.56285624e+01,
  3.06772262e-01, 2.94964746e-01, 3.70976416e-01, 4.59297624e-01,
  4.51836855e-01, 5.87603719e-01, 2.86099293e-01, 2.47791243e-01,
  2.99408166e-01, 4.32746060e-01, 5.27421501e-01, 2.89534285e-01,
  4.30504374e-01, 9.89678273e-01, 6.01385571e-01, 6.38902701e-01,
  1.68818010e-01, 7.50000000e-01, 1.07830459e+00, 2.92598254e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 3.50326241e-01},
 {0.00000000e+00, 1.07002228e+00, 2.05195422e+00, 3.06772262e-01,
  2.43292288e+00, 2.19715052e+00, 2.73911818e-01, 1.14726861e+00,
  1.17343465e+00, 5.79909835e-01, 1.01987501e+00, 3.98183329e-01,
  5.48083233e-01, 1.61030872e+00, 8.05577240e-01, 1.46152937e+00,
  7.42874842e-01, 1.06886174e+00, 9.67471640e-01, 6.11601820e-01,
  2.12712296e-01, 7.50000000e-01, 3.69025897e-01, 1.87658077e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 4.53017475e-01},
 {0.00000000e+00, 1.07474947e+00, 1.82470118e+00, 2.94964746e-01,
  2.19715052e+00, 2.41404373e+00, 2.87098831e-01, 1.04389169e+00,
  1.16525001e+00, 6.26339428e-01, 9.87756703e-01, 4.63568053e-01,
  6.13603630e-01, 1.39293185e+00, 8.79669521e-01, 1.77158803e+00,
  7.82149443e-01, 9.98597151e-01, 9.13853001e-01, 6.60440693e-01,
  2.00637717e-01, 7.50000000e-01, 3.71943861e-01, 2.13407372e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.12682678e-01},
 {0.00000000e+00, 4.46290560e-01, 3.54071966e-01, 3.70976416e-01,
  2.73911818e-01, 2.87098831e-01, 8.02764197e+00, 3.32339414e-01,
  6.59545216e-01, 1.26203255e+00, 2.98137999e-01, 1.51752623e+00,
  1.04274830e+00, 4.46999215e-01, 3.50828021e-01, 3.43119584e-01,
  3.58749728e-01, 4.78353203e-01, 4.88459165e-01, 7.66309017e-01,
  1.07668873e+00, 7.50000000e-01, 4.97590308e+00, 3.11511613e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 1.44043358e+00},
 {0.00000000e+00, 1.33923306e+00, 1.11918465e+00, 4.59297624e-01,
  1.14726861e+00, 1.04389169e+00, 3.32339414e-01, 2.99056680e+00,
  6.15552047e-01, 5.57103314e-01, 6.78509416e-01, 3.94686522e-01,
  5.25369171e-01, 1.08662777e+00, 8.94077901e-01, 7.57447784e-01,
  5.54986868e-01, 1.27845679e+00, 9.98770702e-01, 7.34634729e-01,
  2.00363838e-01, 7.50000000e-01, 2.99613619e-01, 9.19064864e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 4.43694156e-01},
 {0.00000000e+00, 7.32178470e-01, 1.29683470e+00, 4.51836855e-01,
  1.17343465e+00, 1.16525001e+00, 6.59545216e-01, 6.15552047e-01,
  4.47513035e+00, 5.70614299e-01, 9.90160946e-01, 6.20055020e-01,
  6.09441413e-01, 1.43988866e+00, 9.46471632e-01, 1.96194975e+00,
  1.43036819e+00, 8.30226277e-01, 7.40733381e-01, 5.96692220e-01,
  5.44285116e-01, 7.50000000e-01, 9.78877039e-01, 1.51243666e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 6.05136779e-01},
 {0.00000000e+00, 8.88728816e-01, 6.16970770e-01, 5.87603719e-01,
  5.79909835e-01, 6.26339428e-01, 1.26203255e+00, 5.57103314e-01,
  5.70614299e-01, 2.83022396e+00, 6.40800116e-01, 1.74919442e+00,
  1.64986005e+00, 6.59934398e-01, 6.27243837e-01, 6.25316723e-01,
  6.29014923e-01, 7.22513707e-01, 1.01665478e+00, 2.33821719e+00,
  3.02670201e-01, 7.50000000e-01, 7.99863687e-01, 6.25893752e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 2.07538421e+00},
 {0.00000000e+00, 7.66227814e-01, 1.13010166e+00, 2.86099293e-01,
  1.01987501e+00, 9.87756703e-01, 2.98137999e-01, 6.78509416e-01,
  9.90160946e-01, 6.40800116e-01, 2.92620120e+00, 5.18814877e-01,
  1.09961986e+00, 1.25788411e+00, 7.70516817e-01, 1.18387275e+00,
  2.18168472e+00, 9.61645918e-01, 9.95647751e-01, 5.70096500e-01,
  4.50323287e-01, 7.50000000e-01, 3.59829368e-01, 1.07322036e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.55622697e-01},
 {0.00000000e+00, 6.46341735e-01, 4.52559218e-01, 2.47791243e-01,
  3.98183329e-01, 4.63568053e-01, 1.51752623e+00, 3.94686522e-01,
  6.20055020e-01, 1.74919442e+00, 5.18814877e-01, 3.92994485e+00,
  2.33771776e+00, 5.15595552e-01, 5.56038228e-01, 6.65507525e-01,
  5.01008003e-01, 5.24107036e-01, 6.77040993e-01, 1.53178711e+00,
  6.46886104e-01, 7.50000000e-01, 8.14904569e-01, 5.51569446e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 3.27192533e+00},
 {0.00000000e+00, 7.70359541e-01, 6.04775077e-01, 2.99408166e-01,
  5.48083233e-01, 6.13603630e-01, 1.04274830e+00, 5.25369171e-01,
  6.09441413e-01, 1.64986005e+00, 1.09961986e+00, 2.33771776e+00,
  4.40351686e+00, 6.70496227e-01, 6.21012433e-01, 7.95710012e-01,
  9.04471250e-01, 6.98581457e-01, 8.73067477e-01, 1.51039751e+00,
  3.74941590e-01, 7.50000000e-01, 5.67675865e-01, 6.92962139e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 2.13016362e+00},
 {0.00000000e+00, 1.03988401e+00, 1.59985542e+00, 4.32746060e-01,
  1.61030872e+00, 1.39293185e+00, 4.46999215e-01, 1.08662777e+00,
  1.43988866e+00, 6.59934398e-01, 1.25788411e+00, 5.15595552e-01,
  6.70496227e-01, 1.58773724e+00, 8.93169424e-01, 1.19675558e+00,
  9.99684066e-01, 1.17295383e+00, 1.10674855e+00, 6.68196499e-01,
  3.94266130e-01, 7.50000000e-01, 6.20758168e-01, 1.30744195e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.59148347e-01},
 {0.00000000e+00, 1.29408826e+00, 8.46142888e-01, 5.27421501e-01,
  8.05577240e-01, 8.79669521e-01, 3.50828021e-01, 8.94077901e-01,
  9.46471632e-01, 6.27243837e-01, 7.70516817e-01, 5.56038228e-01,
  6.21012433e-01, 8.93169424e-01, 3.84141803e+00, 1.05550534e+00,
  9.60208581e-01, 1.24349521e+00, 1.07552180e+00, 7.58279920e-01,
  2.75343404e-01, 7.50000000e-01, 3.20901839e-01, 9.56295438e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.77523805e-01},
 {0.00000000e+00, 9.03456780e-01, 1.33890747e+00, 2.89534285e-01,
  1.46152937e+00, 1.77158803e+00, 3.43119584e-01, 7.57447784e-01,
  1.96194975e+00, 6.25316723e-01, 1.18387275e+00, 6.65507525e-01,
  7.95710012e-01, 1.19675558e+00, 1.05550534e+00, 2.54088209e+00,
  1.33504131e+00, 8.89983529e-01, 8.32948914e-01, 6.49448913e-01,
  3.35176014e-01, 7.50000000e-01, 3.94549334e-01, 2.10683180e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 6.53380354e-01},
 {0.00000000e+00, 7.00798809e-01, 8.61808201e-01, 4.30504374e-01,
  7.42874842e-01, 7.82149443e-01, 3.58749728e-01, 5.54986868e-01,
  1.43036819e+00, 6.29014923e-01, 2.18168472e+00, 5.01008003e-01,
  9.04471250e-01, 9.99684066e-01, 9.60208581e-01, 1.33504131e+00,
  4.07760419e+00, 9.28286966e-01, 8.20208528e-01, 5.61449723e-01,
  1.65111328e+00, 7.50000000e-01, 3.78768704e-01, 1.02308924e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.39632804e-01},
 {0.00000000e+00, 1.28966428e+00, 1.11706882e+00, 9.89678273e-01,
  1.06886174e+00, 9.98597151e-01, 4.78353203e-01, 1.27845679e+00,
  8.30226277e-01, 7.22513707e-01, 9.61645918e-01, 5.24107036e-01,
  6.98581457e-01, 1.17295383e+00, 1.24349521e+00, 8.89983529e-01,
  9.28286966e-01, 1.44062735e+00, 1.36207698e+00, 7.99071364e-01,
  5.65535441e-01, 7.50000000e-01, 5.19882078e-01, 9.51265394e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.83974254e-01},
 {0.00000000e+00, 1.31725380e+00, 1.03197349e+00, 6.01385571e-01,
  9.67471640e-01, 9.13853001e-01, 4.88459165e-01, 9.98770702e-01,
  7.40733381e-01, 1.01665478e+00, 9.95647751e-01, 6.77040993e-01,
  8.73067477e-01, 1.10674855e+00, 1.07552180e+00, 8.32948914e-01,
  8.20208528e-01, 1.36207698e+00, 1.80640078e+00, 1.06764169e+00,
  3.05841935e-01, 7.50000000e-01, 5.30232118e-01, 8.78596534e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 7.79516037e-01},
 {0.00000000e+00, 1.04499203e+00, 6.37811920e-01, 6.38902701e-01,
  6.11601820e-01, 6.60440693e-01, 7.66309017e-01, 7.34634729e-01,
  5.96692220e-01, 2.33821719e+00, 5.70096500e-01, 1.53178711e+00,
  1.51039751e+00, 6.68196499e-01, 7.58279920e-01, 6.49448913e-01,
  5.61449723e-01, 7.99071364e-01, 1.06764169e+00, 2.69825619e+00,
  2.36794202e-01, 7.50000000e-01, 5.65836205e-01, 6.55650683e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 1.77511928e+00},
 {0.00000000e+00, 2.63609841e-01, 2.96793416e-01, 1.68818010e-01,
  2.12712296e-01, 2.00637717e-01, 1.07668873e+00, 2.00363838e-01,
  5.44285116e-01, 3.02670201e-01, 4.50323287e-01, 6.46886104e-01,
  3.74941590e-01, 3.94266130e-01, 2.75343404e-01, 3.35176014e-01,
  1.65111328e+00, 5.65535441e-01, 3.05841935e-01, 2.36794202e-01,
  5.26606678e+01, 7.50000000e-01, 9.69642510e-01, 2.59266956e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.43022417e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 4.49686798e-01, 4.85608019e-01, 1.07830459e+00,
  3.69025897e-01, 3.71943861e-01, 4.97590308e+00, 2.99613619e-01,
  9.78877039e-01, 7.99863687e-01, 3.59829368e-01, 8.14904569e-01,
  5.67675865e-01, 6.20758168e-01, 3.20901839e-01, 3.94549334e-01,
  3.78768704e-01, 5.19882078e-01, 5.30232118e-01, 5.65836205e-01,
  9.69642510e-01, 7.50000000e-01, 1.03387565e+01, 3.81794898e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 8.10366134e-01},
 {0.00000000e+00, 1.00010336e+00, 1.61300149e+00, 2.92598254e-01,
  1.87658077e+00, 2.13407372e+00, 3.11511613e-01, 9.19064864e-01,
  1.51243666e+00, 6.25893752e-01, 1.07322036e+00, 5.51569446e-01,
  6.92962139e-01, 1.30744195e+00, 9.56295438e-01, 2.10683180e+00,
  1.02308924e+00, 9.51265394e-01, 8.78596534e-01, 6.55650683e-01,
  2.59266956e-01, 7.50000000e-01, 3.81794898e-01, 2.12220220e+00,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 5.73996058e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.70000000e-01, 1.70000000e-01, 1.70000000e-01,
  1.70000000e-01, 1.33300000e+00, 1.70000000e-01, 1.70000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 7.19479599e-01, 5.02168751e-01, 3.50326241e-01,
  4.53017475e-01, 5.12682678e-01, 1.44043358e+00, 4.43694156e-01,
  6.05136779e-01, 2.07538421e+00, 5.55622697e-01, 3.27192533e+00,
  2.13016362e+00, 5.59148347e-01, 5.77523805e-01, 6.53380354e-01,
  5.39632804e-01, 5.83974254e-01, 7.79516037e-01, 1.77511928e+00,
  5.43022417e-01, 7.50000000e-01, 8.10366134e-01, 5.73996058e-01,
  7.50000000e-01, 1.70000000e-01, 7.50000000e-01, 2.91088108e+00}};


FreqRatios*
PSIMatrixFrequencyRatiosNew(const char* matrix_name)
{
    unsigned int i, j;          /* loop indices */
    FreqRatios* retval = NULL;  /* the return value */

    ASSERT(matrix_name);

    retval = (FreqRatios*) Malloc(sizeof(FreqRatios));
    if ( !retval ) {
        return NULL;
    }

    retval->data = (double**) Malloc(sizeof(double*)*PROTEIN_ALPHABET);
    if ( !retval->data ) {
        return PSIMatrixFrequencyRatiosFree(retval);
    }

    for (i = 0; i < PROTEIN_ALPHABET; i++) {
        retval->data[i] = (double*) Malloc(sizeof(double)*PROTEIN_ALPHABET);
        if ( !retval->data[i] ) {
            for (j = 0; j < i; j++) {
                retval->data[j] = MemFree(retval->data[j]);
            }
            return PSIMatrixFrequencyRatiosFree(retval);
        }
    }

    if ( !strcmp(matrix_name, "BLOSUM62") ||
         !strcmp(matrix_name, "BLOSUM62_20")) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = BLOSUM62_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "BLOSUM62_20A")) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = 
                    BLOSUM62_20A_SCALE_MULTIPLIER * BLOSUM62_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "BLOSUM62_20B")) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] =
                    BLOSUM62_20B_SCALE_MULTIPLIER * BLOSUM62_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "BLOSUM45") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = BLOSUM45_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 3;
    } else if ( !strcmp(matrix_name, "BLOSUM80") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = BLOSUM80_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "BLOSUM50") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = BLOSUM50_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "BLOSUM90") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = BLOSUM90_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "PAM30") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = PAM30_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "PAM70") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = PAM70_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else if ( !strcmp(matrix_name, "PAM250") ) {
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            for (j = 0; j < PROTEIN_ALPHABET; j++) {
                retval->data[i][j] = PAM250_FREQRATIOS[i][j];
            }
        }
        retval->bit_scale_factor = 2;
    } else {
        retval = PSIMatrixFrequencyRatiosFree(retval);
    }

    return retval;
}

FreqRatios*
PSIMatrixFrequencyRatiosFree(FreqRatios* freq_ratios)
{
    if ( !freq_ratios )
        return NULL;

    if (freq_ratios->data) {
        Uint4 i;
        for (i = 0; i < PROTEIN_ALPHABET; i++) {
            freq_ratios->data[i] = MemFree(freq_ratios->data[i]);
        }
        freq_ratios->data = MemFree(freq_ratios->data);
    }

    freq_ratios = MemFree(freq_ratios);
    return NULL;
}

/* END of copied code */
/****************************************************************************/
