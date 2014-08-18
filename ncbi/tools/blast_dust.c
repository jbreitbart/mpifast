/* $Id: blast_dust.c,v 1.40 2007/01/23 15:25:22 madden Exp $
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
 * ==========================================================================
 *
 * Authors: Richa Agarwala (based upon versions variously worked upon by Roma 
 *          Tatusov, John Kuzio, and Ilya Dondoshansky).
 *   
 * ==========================================================================
 */

/** @file blast_dust.c
 * A utility to find low complexity NA regions. This parallels functionality 
 * of dust.c from the C toolkit, but without using the structures generated 
 * from ASN.1 spec.
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_dust.c,v 1.40 2007/01/23 15:25:22 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <sequtil.h>
#include <objmgr.h>
#include <seqport.h>
#include <blastkar.h>
#include <blast_dust.h>




/* local, file scope, structures and variables */

/** localcurrents 
 * @todo expand documentation
 */
typedef struct DCURLOC { 
	Int4	cursum, curstart, curend;
	Int2	curlength;
} DCURLOC;

Uint1 *SUM_THRESHOLD;

/* local functions */

static void wo (Int4, Uint1*, Int4, DCURLOC*, Uint1*, Uint1, Int4);
static Boolean wo1 (Int4, Uint1*, Int4, DCURLOC*);
static Int4 dust_triplet_find (Uint1*, Int4, Int4, Uint1*);
static void s_GetSequence(SeqPortPtr spp, Uint1* buf, Int4 buf_length);
static SeqLocPtr s_SlpDust (SeqLocPtr slp, SeqIdPtr id, 
                     DREGION *reg, Int4 nreg, Int4 loopDustMax);


SeqLocPtr BioseqDustEx (BioseqPtr bsp, Int4 level, Int4 windowsize, Int4 linker)
{
	SeqLocPtr	slp = NULL;	/* initialize */

	SeqPortPtr	spp;

	DREGION	PNTR reg, PNTR regold;
	Int4	nreg;
	Int4 loopDustMax = 1;

        Uint1* buf = NULL;

/* error msg stuff */
/*	ErrSetOptFlags (EO_MSG_CODES | EO_SHOW_FILELINE | EO_BEEP); */
	ErrSetOptFlags (EO_MSG_CODES);

/* make sure bioseq is there */
	if (!bsp)
	{
		ErrPostEx (SEV_ERROR, 1, 1,
			  "no bioseq");
                ErrShow ();
		return slp;
	}
	if (!ISA_na (bsp->mol))
	{
		ErrPostEx (SEV_WARNING, 1, 2,
			  "not nucleic acid");
                ErrShow ();
		return slp;
	}

/* place for dusted regions */
	reg = MemNew (sizeof (DREGION));
	if (!reg)
	{
		ErrPostEx (SEV_FATAL, 1, 3,
			   "memory allocation error");
                ErrShow ();
		return slp;
	}
	reg->from = 0;
	reg->to = 0;
	reg->next = NULL;

/* do it */
	spp = SeqPortNew (bsp, 0, bsp->length, 0, Seq_code_ncbi2na);
	if (!spp)
	{
		ErrPostEx (SEV_ERROR, 1, 4,
			   "sequence port open failed");
                ErrShow ();
		MemFree (reg);
		return slp;
	}
        buf = MemNew((1+bsp->length)*sizeof(Uint1));
        s_GetSequence(spp, buf, (1+bsp->length));

	nreg = DustSegs (buf, bsp->length, 0, reg, level, windowsize, linker);
	slp = s_SlpDust (NULL, bsp->id, reg, nreg, loopDustMax);

/* clean up memory */
	SeqPortFree (spp);
        MemFree(buf);
	while (reg)
	{
		regold = reg;
		reg = reg->next;
		MemFree (regold);
	}

	return slp;
}

/* sequence location pointers */

SeqLocPtr SeqLocDustEx (SeqLocPtr this_slp,
		      Int4 level, Int4 window, Int4 linker)
{
	SeqLocPtr	next_slp, slp = NULL;

	SeqIdPtr	id;
	BioseqPtr	bsp;
	SeqPortPtr	spp;

	DREGION	PNTR reg, PNTR regold;
	Int4 nreg;
	Int4 start, end, length;
	Int2 loopDustMax = 0;

/* error msg stuff */
	ErrSetOptFlags (EO_MSG_CODES);

	if (!this_slp)
	{
		ErrPostEx (SEV_ERROR, 2, 1,
			  "no sequence location given for dusting");
                ErrShow ();
		return slp;
	}

/* place for dusted regions */
	regold = reg = MemNew (sizeof (DREGION));
	if (!reg)
	{
		ErrPostEx (SEV_FATAL, 2, 2,
			   "memory allocation error");
                ErrShow ();
		return slp;
	}
	reg->from = 0;
	reg->to = 0;
	reg->next = NULL;

/* count seqlocs */
	next_slp = NULL;
	while ((next_slp = SeqLocFindNext (this_slp, next_slp)) != NULL)
			loopDustMax++;
	if (!loopDustMax)
	{
		ErrPostEx (SEV_ERROR, 2, 3,
			   "can not find next seq loc");
                ErrShow ();
	}

/* loop for dusting as needed */
	next_slp = NULL;
	while ((next_slp = SeqLocFindNext (this_slp, next_slp)) != NULL)
	{
                Uint1* buf = NULL;
/* offsets into actual sequence */
		start = SeqLocStart (next_slp);
		end = SeqLocStop (next_slp);

/* if all goes okay should get a seqport pointer */
		id = SeqLocId (next_slp);
		if (!id)
		{
			ErrPostEx (SEV_ERROR, 2, 4,
				  "no bioseq id");
			ErrShow ();
			continue;
		}
		bsp = BioseqLockById (id);
		if (!bsp)
		{
			ErrPostEx (SEV_ERROR, 2, 5,
				  "no bioseq");
			ErrShow ();
			continue;
		}
		if (!ISA_na (bsp->mol))
		{
			ErrPostEx (SEV_WARNING, 2, 6,
				  "not nucleic acid");
			ErrShow ();
			BioseqUnlock (bsp);
			continue;
		}
		spp = SeqPortNew (bsp, start, end, 0, Seq_code_ncbi2na);
		BioseqUnlock (bsp);
		if (!spp)
		{
			ErrPostEx (SEV_ERROR, 2, 7,
				   "sequence port open failed");
			ErrShow ();
			continue;
		}
		length = spp->totlen;
                buf = MemNew((1+length)*sizeof(Uint1));
                s_GetSequence(spp, buf, (1+length));

		nreg = DustSegs (buf, length, start, reg,
				  level, window, linker);
                buf = MemFree(buf);
		SeqPortFree (spp);
		slp = s_SlpDust (slp, id, reg, nreg, loopDustMax);
/* find tail - this way avoids referencing the pointer */
		while (reg->next) reg = reg->next;

	}

/* clean up memory */
	reg = regold;
	while (reg)
	{
		regold = reg;
		reg = reg->next;
		MemFree (regold);
	}

	return slp;
}


static void
s_GetSequence(SeqPortPtr spp, Uint1* buf, Int4 buf_length)
{
        Uint1 residue;
	int index=0;

        ASSERT(spp && buf && buf_length > 0);

	while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if (residue == SEQPORT_VIRT)
			continue;
		buf[index] = residue;
		index++;
	}
        ASSERT(index < buf_length);
	buf[index] = NULLB;
        return;
}


/* Produce locations from dust results. */

static SeqLocPtr s_SlpDust (SeqLocPtr slp, SeqIdPtr id,
			  DREGION *reg, Int4 nreg, Int4 loopDustMax)
{
	SeqIntPtr	sintp;
        Int4            i;
        Boolean         flagNoPack;
        ValNodePtr      vnp = NULL;

/* point to dusted locations */
	if (nreg)
	{

/* loopDustMax == 1 forces PACKED_INT IN - PACKED_INT OUT as needed	*/
		flagNoPack = FALSE;
		if (nreg == 1 && loopDustMax == 1) flagNoPack = TRUE;

		if (!slp)
		{
			if ((slp = ValNodeNew (NULL)) == NULL)
			{
				ErrPostEx (SEV_ERROR, 6, 1,
					   "val node new failed");
				ErrShow ();
				return slp;
			}
		}

		if (flagNoPack)
		{
			slp->choice = SEQLOC_INT;
		}
		else
		{
			slp->choice = SEQLOC_PACKED_INT;
		}

		for (i = 0; i < nreg; i++)
		{
			sintp = SeqIntNew ();
			if (!sintp)
			{
				ErrPostEx (SEV_FATAL, 6, 2,
					   "memory allocation error");
				ErrShow ();
				return slp;
			}
			sintp->id = SeqIdDup (id);
			sintp->from = reg->from;
			sintp->to = reg->to;
			if (!flagNoPack) ValNodeAddPointer
					(&vnp, SEQLOC_INT, sintp);
			reg = reg->next;
		}

		if (flagNoPack)
		{
			slp->data.ptrvalue = (Pointer) sintp;
		}
		else
		{
			slp->data.ptrvalue = vnp;
		}
	}
	return slp;
}

/* entry point for dusting */

Int4 DustSegs (Uint1* sequence, Int4 length, Int4 start,
		       DREGION* reg,
		       Int4 level, Int4 windowsize, Int4 linker)
{
   Int4    len;
   Int4	i;
   Uint1* seq;
   DREGION* regold = NULL;
   DCURLOC	cloc;
   Int4	nreg;
   /* Default values. */
   const int kDustLevel = 20;
   const int kDustWindow = 64;
   const int kDustLinker = 1;
   
   /* defaults are more-or-less in keeping with original dust */
   if (level < 2 || level > 64) level = kDustLevel;
   if (windowsize < 8 || windowsize > 64) windowsize = kDustWindow;
   if (linker < 1 || linker > 32) linker = kDustLinker;
   
   nreg = 0;
   seq = (Uint1*) calloc(1, windowsize);			/* triplets */
   if (!seq) {
      return -1;
   }
   SUM_THRESHOLD = (Uint1 *) calloc(windowsize, sizeof(Uint1));  
   if (!SUM_THRESHOLD) {
      return -1;
   }
   SUM_THRESHOLD[0] = 1;
   for (i=1; i < windowsize; i++)
       SUM_THRESHOLD[i] = (level * i)/10;

   if (length < windowsize) windowsize = length;

   /* Consider smaller windows in beginning of the sequence */
   for (i = 2; i <= windowsize-1; i++) {
      len = i-1;
      wo (len, sequence, 0, &cloc, seq, 1, level);
      
      if (cloc.cursum*10 > level*cloc.curlength) {
         if (nreg &&
             regold->to + linker >= cloc.curstart+start &&
             regold->from <= cloc.curend + start + linker) {
            /* overlap windows nicely if needed */
            if (regold->to < cloc.curend +  start)
                regold->to = cloc.curend +  start;
            if (regold->from > cloc.curstart + start)
                regold->from = cloc.curstart + start;
         } else	{
            /* new window or dusted regions do not overlap */
            reg->from = cloc.curstart + start;
            reg->to = cloc.curend + start;
            regold = reg;
            reg = (DREGION*) calloc(1, sizeof(DREGION));
            if (!reg) {
               MemFree(seq);
               MemFree(SUM_THRESHOLD);
               return -1;
            }
            reg->next = NULL;
            regold->next = reg;
            nreg++;
         }
      }				/* end 'if' high score	*/
   }					/* end for */

   for (i = 1; i < length-2; i++) {
      len = (Int4) ((length > i+windowsize) ? windowsize : length-i);
      len -= 2;
      if (length >= i+windowsize)
          wo (len, sequence, i, &cloc, seq, 2, level);
      else /* remaining portion of sequence is less than windowsize */
          wo (len, sequence, i, &cloc, seq, 3, level);
      
      if (cloc.cursum*10 > level*cloc.curlength) {
         if (nreg &&
             regold->to + linker >= cloc.curstart+i+start &&
             regold->from <= cloc.curend + i + start + linker) {
            /* overlap windows nicely if needed */
            if (regold->to < cloc.curend + i + start)
                regold->to = cloc.curend + i + start;
            if (regold->from > cloc.curstart + i + start)
                regold->from = cloc.curstart + i + start;
         } else	{
            /* new window or dusted regions do not overlap */
            reg->from = cloc.curstart + i + start;
            reg->to = cloc.curend + i + start;
            regold = reg;
            reg = (DREGION*) calloc(1, sizeof(DREGION));
            if (!reg) {
               MemFree(seq);
               MemFree(SUM_THRESHOLD);
               return -1;
            }
            reg->next = NULL;
            regold->next = reg;
            nreg++;
         }
      }				/* end 'if' high score	*/
   }					/* end for */
   MemFree (seq);
   MemFree(SUM_THRESHOLD);
   return nreg;
}

static void wo (Int4 len, Uint1* seq_start, Int4 iseg, DCURLOC* cloc, 
                Uint1* seq, Uint1 FIND_TRIPLET, Int4 level)
{
	Int4 smaller_window_start, mask_window_end;
        Boolean SINGLE_TRIPLET;

	cloc->cursum = 0;
	cloc->curlength = 1;
	cloc->curstart = 0;
	cloc->curend = 0;

	if (len < 1)
		return;

        /* get the chunk of sequence in triplets */
	if (FIND_TRIPLET==1) /* Append one */
	{
		seq[len-1] = seq[len] = seq[len+1] = 0;
		dust_triplet_find (seq_start, iseg+len-1, 1, seq+len-1);
	}
	if (FIND_TRIPLET==2) /* Copy suffix as prefix and find one */
	{
		memmove(seq,seq+1,(len-1)*sizeof(Uint1));
		seq[len-1] = seq[len] = seq[len+1] = 0;
		dust_triplet_find (seq_start, iseg+len-1, 1, seq+len-1);
	}
	if (FIND_TRIPLET==3) /* Copy suffix */
		memmove(seq,seq+1,len*sizeof(Uint1));

        /* dust the chunk */
	SINGLE_TRIPLET = wo1 (len, seq, 0, cloc); /* dust at start of window */

        /* consider smaller windows only if anything interesting 
           found for starting position  and smaller windows have potential of
           being at higher level */
        if ((cloc->cursum*10 > level*cloc->curlength) && (!SINGLE_TRIPLET)) {
		mask_window_end = cloc->curend-1;
		smaller_window_start = 1;
                while ((smaller_window_start < mask_window_end) &&
                       (!SINGLE_TRIPLET)) {
			SINGLE_TRIPLET = wo1(mask_window_end-smaller_window_start,
                             seq+smaller_window_start, smaller_window_start, cloc);
                	smaller_window_start++;
	        }
	}

	cloc->curend += cloc->curstart;
}

/** returns TRUE if there is single triplet in the sequence considered */
static Boolean wo1 (Int4 len, Uint1* seq, Int4 iwo, DCURLOC* cloc)
{
   Uint4 sum;
	Int4 loop;

	Int2* countsptr;
	Int2 counts[4*4*4];
	Uint1 triplet_count = 0;

	memset (counts, 0, sizeof (counts));
/* zero everything */
	sum = 0;

/* dust loop -- specific for triplets	*/
	for (loop = 0; loop < len; loop++)
	{
		countsptr = &counts[*seq++];
		if (*countsptr)
		{
			sum += (Uint4)(*countsptr);

		    if (sum >= SUM_THRESHOLD[loop])
		    {
			if ((Uint4)cloc->cursum*loop < sum*cloc->curlength)
			{
				cloc->cursum = sum;
				cloc->curlength = loop;
				cloc->curstart = iwo;
				cloc->curend = loop + 2; /* triplets */
			}
		    }
		}
		else
			triplet_count++;
		(*countsptr)++;
	}

	if (triplet_count > 1)
		return(FALSE);
	return(TRUE);
}

/** Fill an array with 2-bit encoded triplets.
 * @param seq_start Pointer to the start of the sequence in blastna 
 *                  encoding [in]
 * @param icur Offset at which to start extracting triplets [in]
 * @param max Maximal length of the sequence segment to be processed [in]
 * @param s1 Array of triplets [out]
 * @return How far was the sequence processed?
 */
static Int4 
dust_triplet_find (Uint1* seq_start, Int4 icur, Int4 max, Uint1* s1)
{
   Int4 n;
   Uint1* s2,* s3;
   Int2 c;
   Uint1* seq = &seq_start[icur];
   Uint1 end_byte = ncbi4na_to_blastna[NULLB];
   const int k_NCBI2na_mask =  0x03;
   
   n = 0;
   
   s2 = s1 + 1;
   s3 = s1 + 2;
   
   /* set up 1 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= k_NCBI2na_mask;
   *s1 |= c;
   *s1 <<= 2;
   
   /* set up 2 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= k_NCBI2na_mask;
   *s1 |= c;
   *s2 |= c;
   
   /* triplet fill loop */
   while (n < max && (c = *seq++) != end_byte) {
      c &= k_NCBI2na_mask;
      *s1 <<= 2;
      *s2 <<= 2;
      *s1 |= c;
      *s2 |= c;
      *s3 |= c;
      s1++;
      s2++;
      s3++;
      n++;
   }
   
   return n;
}
