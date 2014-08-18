static char const rcsid[] = "$Id: salptool.c,v 6.57 2008/07/07 18:52:31 bollin Exp $";

#include <sequtil.h> /* SeqIdDupList */
#include <salpedit.h>
#include <salptool.h>
#include <salpacc.h>
#include <blast.h>
#include <simutil.h>
#include <alignval.h>
#include <blastpri.h>
#include <txalign.h>
#include <sqnutils.h>
#include <satutil.h>
#include <alignmgr.h>
#include <alignmgr2.h>
#include <spidey.h>
#include <salsap.h>

#define MIN_SEG_SIZE 10 /* used by MergeTwoDspBySIM4 , check_align_match_diagnol<-ModFilterDenseSegAlign<-ModifyAlignList<-SeqAlignConsistentDiagFilter */


/* Blast can now (2/1999) handle longer queries more efficiently.. though
 there is still a hardcoded limit to the number of HSP's */


#define PERCENT_A_RES	90 /* used by get_polyA_index<-check_polyA_tail<-SeqLocTrimPolyAtail<-compute_alignment  */


#define MAX_HANG 100 /* used in need_recompute<-compute_alignment */

#define MAX_EXTEND_OVERHANG	50 /* used by compute_alignment and SeqAlignSetFlobalFromLocal */

static Boolean is_bad_align (BioseqPtr m_bsp, BioseqPtr s_bsp, SeqAlignPtr align, Int4 min_align_len,Int4 loc_len);

static void OrderInt4(Int4Ptr x, Int4Ptr y) /* replace jzmisc: swap */
{
  Int4 temp;

	if((*x) > (*y)){
	  temp = *x;
	  *x = *y;
	  *y = temp;
	}
}

static void SeqAnnotWrite(SeqAlignPtr align)
{
        SeqAnnotPtr annot;
        AsnIoPtr aip;

        annot = SeqAnnotNew();
        annot->type = 2;
        annot->data = align;

        aip = AsnIoOpen("temp2.sat", "w");
        SeqAnnotAsnWrite(annot, aip, NULL);
        AsnIoClose(aip);

        annot->data = NULL;
        SeqAnnotFree(annot);
}

static Boolean check_align_match_diagnol(SeqAlignPtr align, ValNodePtr ddp_list, Uint1 strand)
{
	DenseSegPtr dsp;
	Int2 i;
	Int4 m_start, s_start;
	Boolean minus;
	Int4 total_len, match_len;
	ValNodePtr curr;
	DenseDiagPtr ddp;

	/*check for orientations first*/
	dsp = (DenseSegPtr) align->segs;
	minus = FALSE;
	if(dsp->strands !=NULL)
	{
		if(dsp->strands[0] != dsp->strands[1])
		{
			if(dsp->strands[1] == Seq_strand_minus  || 
				dsp->strands[0] == Seq_strand_minus)
				minus = TRUE;
		}
	}

	if(strand == Seq_strand_minus && !minus)
		return FALSE;
	if(minus && strand != Seq_strand_minus)
		return FALSE;

	total_len = 0;
	match_len = 0;
	for(i = 0; i<dsp->numseg; ++i)
	{
		m_start = dsp->starts[2*i];
		s_start = dsp->starts[2*i+1];
		if(m_start != -1 && s_start != -1 && dsp->lens[i] > MIN_SEG_SIZE)
		{
			total_len += dsp->lens[i];
			for(curr = ddp_list; curr != NULL; curr = curr->next)
			{
				ddp = (DenseDiagPtr) curr->data.ptrvalue;
				if(ddp->starts[0] == m_start 
					&& ddp->starts[1] == s_start)
				{
					match_len += dsp->lens[i];
					break;
				}
			}
		}
	}

	if(match_len < 2*MIN_SEG_SIZE || total_len < 2*MIN_SEG_SIZE)
		return FALSE;

	if(match_len == total_len)
		return TRUE;

	/*at least 40% need to be matching the fragment*/
	if(match_len > 4 * MIN_SEG_SIZE)
	{
		if((match_len *100 )/40 >= total_len)
			return TRUE;
	}
	else if(total_len < 4* MIN_SEG_SIZE)
		if((match_len *100 )/80 >= total_len)
			return TRUE;

	return FALSE;
}

				
			
/*
*
*	Filter the Dense-seg alignment. Get rid of the alignments that were not 
*	in the main diagnol
*	all the alignments MUST be for the SAME sequences
*	sip is the id of the sequence that is used for the primary ordering
*/
static Boolean ModFilterDenseSegAlign(SeqAlignPtr PNTR align, SeqIdPtr sip)
{
	SeqAlignPtr ddp_align;
	ValNodePtr head = NULL;
	SeqAlignPtr curr, prev, next;
	Uint1 strand;
        /* 	DenseSegPtr dsp; */
	BioseqPtr bsp;

	if(align == NULL || *align == NULL || sip == NULL)
		return FALSE;
	/*check if it has more than two alignments*/
	curr = *align;
	if(curr->segtype != 2)
		return FALSE;
	/* dsp = curr->segs;
	if(dsp->ids->next->choice == SEQID_LOCAL)
	{
           	SeqIdPtr t_sip;
          	ObjectIdPtr oip;
		t_sip = dsp->ids->next;
		oip = t_sip->data.ptrvalue;
		printf("%s\n", oip->str);
	} */
	

	ddp_align = SeqAlignConvertDspToDdpList(*align);
	if(ddp_align == NULL)
		return FALSE;

	bsp = BioseqLockById(sip);
	if(bsp == NULL)
		return FALSE;
	head = FilterSeqAlign(ddp_align, sip, &strand);
	BioseqUnlock(bsp);
	if(head == NULL)
	{
		SeqAlignFree(ddp_align);
		return FALSE;
	}

	curr = *align;
	prev = NULL;
	while(curr)
	{
		next = curr->next;
		if(!check_align_match_diagnol(curr, head, strand))
		{
			if(prev == NULL)
				*align = next;
			else
				prev->next = next;
			curr->next = NULL;
			SeqAlignFree(curr);
		}
		else
			prev = curr;
		curr = next;
	}

	ValNodeFree(head);
	SeqAlignFree(ddp_align);

	return TRUE;
}

static Boolean reload_seq_loc(SeqLocPtr slp, Int4 from, Int4 to, Uint1 strand)
{
	SeqIntPtr sint;

	if(slp == NULL || slp->choice != SEQLOC_INT)
		return FALSE;
	sint = (SeqIntPtr) slp->data.ptrvalue;
	sint->from = from;
	sint->to = to;
	sint->strand = strand;

	return TRUE;
}
/*
*
*	modify the locations on loc_1 and loc_2 to record the from, to for the 
*	end segment. if (head), return the first segment >= min_seg_len, else 
*	return the last segment
*/
static Boolean get_end_seg(DenseSegPtr dsp, Int4 min_seg_len, SeqLocPtr loc_1, SeqLocPtr loc_2, Boolean head)
{
	Int4 i;
	Int4 start_1, start_2;

	if(dsp == NULL || loc_1 == NULL || loc_2 == NULL)
		return FALSE;
	if(head)
		i = 0;
	else
		i = dsp->numseg -1;

	for(; head? i<dsp->numseg : i>=0; head? ++i: --i)
	{
		start_1 = dsp->starts[2*i];
		start_2 = dsp->starts[2*i+1];
		if(start_1 != -1 && start_2 != -1)
		{
			if(dsp->lens[i] >= min_seg_len)
			{
				reload_seq_loc(loc_1, start_1, 
					start_1+dsp->lens[i]-1, dsp->strands[0]);
				reload_seq_loc(loc_2, start_2, 
					start_2+dsp->lens[i]-1, dsp->strands[1]);

				return TRUE;
			}
		}
	}

	return FALSE;
}


/*
*
*	align_1 and align_2 will be  re-sorted by the ascending order in the first 
*	sequence in the alignment. The first sequence needs to aligned 
*	in the plus strand. They have to be alignment on the same 
*	diagnol
*	return NULL if any merge fails
*
*/
static SeqAlignPtr MergeTwoDspBySIM4(SeqAlignPtr align_1, SeqAlignPtr align_2,Int4 MaxGap)
{
	SeqLocPtr loc_1_1, loc_1_2, loc_2_1, loc_2_2;
	Int4 start, stop;
	Int4 match_start, match_stop;
	Boolean found = FALSE;
	BioseqPtr bsp;
	SeqAlignPtr align, t_align;
	DenseSegPtr dsp_1, dsp_2;
	BioseqPtr bsp_1, bsp_2;
	Int4 max_align_size;
	Uint1 strand;
	Int4 gap=0,gap_1, gap_2;
	Int4 start_1, stop_1, start_2, stop_2;
	SeqIdPtr sip;
	Boolean reverse = FALSE;
	Boolean has_overlap;
	Int2 order;

	if(align_1 == NULL || align_2 == NULL)
		return NULL;

	if(align_1->segtype != 2 || align_2->segtype != 2)
		return NULL;
        if(SeqAlignStart(align_1,0)>SeqAlignStart(align_2,0)) {
            VoidPtr segs;
            segs = align_1->segs;
            align_1->segs = align_2->segs;
            align_2->segs = segs;
        }

	dsp_1 = (DenseSegPtr) align_1->segs;
	if(dsp_1 && dsp_1->strands &&  dsp_1->strands[0] == Seq_strand_minus)
	{
		SeqAlignReverse(align_1, 0);
		reverse = TRUE;
	}
	dsp_2 = (DenseSegPtr) align_2->segs;
	if(dsp_2 && dsp_2->strands &&  dsp_2->strands[0] == Seq_strand_minus)
	{
		SeqAlignReverse(align_2, 0);
		reverse = TRUE;
	}

	max_align_size = MIN(SeqAlignLength(align_1), 
		SeqAlignLength(align_2));

	if(dsp_1 == NULL || dsp_2 == NULL)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, dsp_1->ids, &start_2, &stop_2, &strand); 
        /* Redundant alignment covering same region of query */
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids->next, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, dsp_1->ids->next, &start_2, &stop_2, &strand); 
        /* Redundant alignment covering same region of subject */
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids, &start, &stop, &strand); 
	loc_1_1 = SeqLocIntNew(start, stop, strand, dsp_1->ids);
	SeqAlignStartStopById(align_1, dsp_1->ids->next, &start, &stop, &strand); 
	loc_1_2 = SeqLocIntNew(start, stop, strand, dsp_1->ids->next);

	/* loc_1_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
	loc_1_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
	if(!get_end_seg(dsp_1, MIN_SEG_SIZE*2, loc_1_1, loc_1_2, FALSE))
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	} */

	loc_2_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
	loc_2_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
	if(!get_end_seg(dsp_2, MIN_SEG_SIZE*2, loc_2_1, loc_2_2, TRUE))
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		SeqLocFree(loc_2_1);
		SeqLocFree(loc_2_2);
		return NULL;
	}

	found = FALSE;
	has_overlap = FALSE;
	if(select_overlap_loc(loc_1_1, loc_2_1, loc_1_2, loc_2_2, &order))
	{
		if(order == 0)
			is_loc_overlap(loc_1_1, loc_2_1, &start, &stop, 0, max_align_size, NULL);
		else
			is_loc_overlap(loc_1_2, loc_2_2, &start, &stop, 0, max_align_size, NULL);
		found = TRUE;
		has_overlap = TRUE;
	}

	if(!found)
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		loc_1_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
		loc_1_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
		if(!get_end_seg(dsp_1, MIN_SEG_SIZE*2, loc_1_1, loc_1_2, FALSE))
		{
			SeqLocFree(loc_1_1);
			SeqLocFree(loc_1_2);
			SeqLocFree(loc_2_1);
			SeqLocFree(loc_2_2);
			return NULL;
		} 
		/*check for the overlap between the two segments*/
		if(is_loc_overlap(loc_1_1, loc_2_1, &start, &stop, 0, max_align_size, NULL))
		{
			found = TRUE;
			order = 0;
		}
		else if(is_loc_overlap(loc_1_2, loc_2_2, &start, &stop, 0, max_align_size, NULL))
		{
			found = TRUE;
			order = 1;
		}
		else
		{
			start_1 = start;
			stop_1 = stop;
			start_2 = start;
			stop_2 = stop;

			if(is_loc_overlap(loc_1_1, loc_2_1, &start_1, &stop_1, MIN_SEG_SIZE*10, max_align_size, &gap_1))
			{
				found = TRUE;
				order = 0;
                                gap = gap_1;
			}
			if(is_loc_overlap(loc_1_2, loc_2_2, &start_2, &stop_2, MIN_SEG_SIZE*10, max_align_size, &gap_2))
			{
				if(!found || gap_2 < gap_1)
				{
					order = 1;
					found = TRUE;
                                        gap = gap_2;
				}
			}
			if(found)
			{
				if(order == 0)
				{
					start = start_1;
					stop = stop_1;
				}
				else
				{
					start = start_2;
					stop = stop_2;
				}
			}
		}
	}

	SeqLocFree(loc_1_1);
	SeqLocFree(loc_1_2);
	SeqLocFree(loc_2_1);
	SeqLocFree(loc_2_2);

	if(!found)
	{
		return NULL;
	}

	/*
		===================>
			====================>
			^^^^^^^^^^^^^^^^^
			           |----------
                                    need to have some offset into the 
				    second sequence
	*/
	if(order == 0)
		sip = dsp_1->ids;
	else
		sip = dsp_1->ids->next;
	SeqAlignStartStopById(align_1, sip, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, sip, &start_2, &stop_2, &strand); 

        start = MAX(start_2, start_1);
        stop = MIN(stop_1, stop_2);
        OrderInt4(&start,&stop); /* Need reordering for gaps */

	if(dsp_1->strands[order] == Seq_strand_minus)
	{
		match_start = find_matching_position(dsp_1, stop, order);
		match_stop = find_matching_position(dsp_2, start, order);
	}
	else
	{
		match_start = find_matching_position(dsp_1, start, order);
		match_stop = find_matching_position(dsp_2, stop, order);
	}

	if(match_start == -1 || match_stop == -1)
		return NULL;
	bsp = NULL;
	if(order == 0)
		bsp = BioseqLockById(dsp_1->ids->next);
	else
		bsp = BioseqLockById(dsp_1->ids);

	if(bsp == NULL)
		return NULL;
        { 
            Int4 old_start=start,old_stop = stop;
            start = MAX(start-100,MIN(start_1,start_2));
            stop = MIN(stop+100,MAX(stop_1,stop_2));
            if(match_start<=match_stop) {
                match_start -=(old_start-start);
                match_stop +=(stop-old_stop);
                if(match_start<0) {
                    start+=match_start;
                    match_start=0;
                } 
                if(match_stop>bsp->length-1) {
                    stop-=(match_stop-(bsp->length-1));
                    match_stop = bsp->length-1;
                }
            }
            else {
                match_start +=(old_start-start);
                match_stop -=(stop-old_stop);
                if(match_stop<0) {
                    stop+=match_stop;
                    match_stop=0;
                } 
                if(match_start>bsp->length-1) {
                    start-=(match_start-(bsp->length-1));
                    match_start = bsp->length-1;
                }

            }

        }
	BioseqUnlock(bsp);


	OrderInt4(&match_start, &match_stop);
	OrderInt4(&start, &stop);
        /* HS Fix to Disallow Very Long Gaps */
        if(abs(abs(stop-start) -abs(match_stop-match_start))>MaxGap+3000)
            return NULL;
	if(order == 0)
	{
		loc_1_1 = SeqLocIntNew(start, stop, dsp_1->strands[order], dsp_1->ids);
		loc_1_2 = SeqLocIntNew(match_start, match_stop, 
			dsp_1->strands[1-order], dsp_1->ids->next);
	}
	else
	{
		loc_1_1 = SeqLocIntNew(match_start, match_stop, 
			dsp_1->strands[1-order], dsp_1->ids);
		loc_1_2 = SeqLocIntNew(start, stop, 
			dsp_1->strands[order], dsp_1->ids->next);
	}

	bsp_1 = BioseqLockById(SeqLocId(loc_1_1));
	if(bsp_1 == NULL)
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	}
	bsp_2 = BioseqLockById(SeqLocId(loc_1_2));
	if(bsp_2 == NULL)
	{
		BioseqUnlock(bsp_1);
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	}

	align = SIM4ALN_choice(loc_1_1, loc_1_2, 200, 8);
	SeqLocFree(loc_1_1);
	SeqLocFree(loc_1_2);

	if(align != NULL)
	{
		if(order == 1)
		{
			start = match_start;
			stop = match_stop;
			order = 0;
		}
		t_align = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_1, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);

		/* for debugging purpose*/
		/* t_align->next = align;
		align->next = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_2, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);
		return t_align;    */

		if(MergeTwoAlignList(t_align, &align, start, stop, order))
		{
			if(align != NULL)
				SeqAlignFree(align);
			align = t_align;
			align->next = NULL;

			t_align = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_2, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);
			t_align->next = NULL;
			MergeTwoAlignList(align, &t_align, start, stop, order);
			if(t_align != NULL)
				SeqAlignFree(t_align);
		}
		else
		{
			if(align != NULL)
				align = SeqAlignSetFree(align);
			SeqAlignFree(t_align);
		}
	}
	/*for debugging */
	/* else
	{
		align_1->next = align_2;
		align_2->next = NULL;
		SeqAnnotWrite(align_1);
		exit(1);
	} */
	BioseqUnlock(bsp_1);
	BioseqUnlock(bsp_2);

	if(reverse && align != NULL)
		SeqAlignReverse(align, 1);
	return align;
}

		

static SeqAlignPtr MergeToOneAlignment(SeqAlignPtr PNTR palign,Int4 MaxGap) 
{
	SeqAlignPtr next, prev, align;
	SeqAlignPtr merge_align;
	SeqAlignPtr curr;

	if(palign == NULL || *palign == NULL)
		return NULL;
	align = *palign;
	if(align->segtype != 2)
		return NULL;
	*palign = SeqAlignSortByRegion(align, 0);

	prev = NULL;
	align = *palign;

	while(align)
	{
		next = align->next;
		if(next != NULL)
		{
			merge_align = NULL;
			curr = *palign;
			/* if(next->next == NULL)
				printf("stop here\n"); */
			merge_align = MergeTwoDspBySIM4(align, next,MaxGap);
			if(merge_align != NULL)
			{
				if(prev == NULL)
					*palign = merge_align;
				else
					prev->next = merge_align;
				
				curr = merge_align;
				while(curr->next != NULL)
				{
					prev = curr;
					curr = curr->next;
				}
				curr->next = next->next;
				next->next = NULL;
				SeqAlignSetFree(align);
				align = curr;
			}
			else
			{
				prev = align;
				align = align->next;
			}
		}
		
		else
			align = align->next;
	}

	if(*palign != NULL && (*palign)->next!=NULL)
		*palign = SeqAlignSortByLength(*palign);
	return (*palign);
}

/*
*
*	Filter out any alignments that were not in the main 
*	diagnol
*
*/
static Boolean ModifyAlignList(SeqAlignPtr PNTR palign)
{
	SeqAlignPtr align, next;
	DenseSegPtr dsp;
	SeqAlignPtr h_align = NULL;

	if(palign == NULL || *palign == NULL)
		return FALSE;

	align = *palign;
	while(align)
	{
		next = align->next;
		if(align->segtype == 2)
		{
			dsp = (DenseSegPtr) align->segs;
			align->next = SeqAlignExtractByIds(&next,dsp->ids, dsp->ids->next);
			ModFilterDenseSegAlign(&align, dsp->ids);
			if(align != NULL)
				h_align = SeqAlignLink(align, h_align);
		}
		else
		{
			align->next = NULL;
			SeqAlignFree(align);
		}
		align = next;
	}

	*palign = h_align;
	return TRUE;
}

				

			


static Int2 get_polyA_index (Uint1Ptr res_buf, Int2 len)
{
	Int2 i, total;
	Uint1 res;
	Int2 count_A;

	count_A = 0;
	total = 0;

	for( i= len -1; i>=0; --i)
	{
		res = res_buf[i];
		if(res == 'a' || res == 'A' || res == 'n' || res == 'N')
		{
			++count_A;
			++total;
		}
		else
		{
			if(total > 5 && count_A * 100 < (total + 1) * PERCENT_A_RES)
				break;
			else
				++total;
		}
	}
	/*trim the non-A and non-Ns*/
	for(i = MAX(i, 0); i<len; )
	{
		res = res_buf[i];
		if(res != 'n' && res != 'N' && res != 'a' && res != 'A')
			++i;
		else	/*N/A residue*/
			break;
	}

	if(len-1 - i > 2)
		return i;

	return -1;
}


static SeqLocPtr check_polyA_tail (BioseqPtr bsp)
{
	Int2 len = 100;
	SeqPortPtr spp;
	Uint1 res_buf[101];
	Uint1 res;
	Int2 i;
	SeqLocPtr end_loc, begin_loc;
	Int4 start, stop;

	begin_loc = NULL;
	end_loc = NULL;

	len = MIN(bsp->length-1, 100);
	/*check for the end of the sequence for polyA tail*/
	spp = SeqPortNew(bsp, bsp->length-1-(len-1), bsp->length-1, 
		Seq_strand_plus, Seq_code_iupacna);
	i = 0;
	while((res = SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if(IS_ALPHA(res))
			res_buf[i++] = res;
	}
	SeqPortFree(spp);
	if(i > 10)
	{
		len = i;
		i = get_polyA_index (res_buf, len);
		if(i != -1)
		{
			stop = bsp->length -1;
			start = bsp->length - 1 - (len -1) + i;

			end_loc = SeqLocIntNew(start, stop, 
				Seq_strand_plus, SeqIdFindBest(bsp->id, 0));
		}

	}

	/*check for the beginning on the reverse complement of the sequence 
          for polyA tail. Usually the 3' clone of the EST*/
	spp = SeqPortNew(bsp, 0, len-1, 
		Seq_strand_minus, Seq_code_iupacna);
	i = 0;
	while((res = SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if(IS_ALPHA(res))
			res_buf[i++] = res;
	}
	SeqPortFree(spp);
	if(i > 10)
	{
		len = i;
		i = get_polyA_index (res_buf, len);
		if(i != -1)
		{
			stop = len -1 - i;
			start = 0;

			begin_loc = SeqLocIntNew(start, stop, 
				Seq_strand_minus, SeqIdFindBest(bsp->id, 0));
		}

	}
	if(begin_loc == NULL)
		return end_loc;
	else
	{
		if(end_loc != NULL)
		{
			if(SeqLocLen(end_loc) > SeqLocLen(begin_loc))
			{
				SeqLocFree(begin_loc);
				return end_loc;
			}
			else
			{
				SeqLocFree(end_loc);
				return begin_loc;
			}
		}
		return begin_loc;
	}
}

static void SeqLocTrimPolyATail(SeqLocPtr loc_2, BioseqPtr bsp_2) {
	SeqLocPtr poly_loc;
	SeqIntPtr sint;

	poly_loc = check_polyA_tail (bsp_2);
	if(poly_loc != NULL)
	{
		sint = (SeqIntPtr) loc_2->data.ptrvalue;
		if(SeqLocStart(poly_loc) > (bsp_2->length)/2)
		    sint->to = SeqLocStart(poly_loc) -1;
		else
			sint->from = SeqLocStop(poly_loc) + 1;
		SeqLocFree(poly_loc);
	}
}

static Boolean is_bad_blast_alignment(SeqAlignPtr align, SeqLocPtr slp1, SeqLocPtr slp2, Int4 min_align_len)
{
	SeqAlignPtr best_align;
	BioseqPtr bsp1, bsp2;
	Int4 max_score;

	best_align = find_best_align(align, &max_score);
	if(best_align == NULL)
		return TRUE;
	bsp1 = BioseqFind(SeqLocId(slp1));
	bsp2 = BioseqFind(SeqLocId(slp2));

	if(bsp1 == NULL || bsp2 == NULL)
		return (max_score > 50);
	return is_bad_align(bsp1, bsp2, best_align,
		min_align_len, MIN(SeqLocLen(slp1), SeqLocLen(slp2)));
}

/*
  Blindly reverse the left-right order of the seqalign... to maintain
  'display' order.
 */
NLM_EXTERN void SeqAlignReverseOrder(SeqAlignPtr align)
{
	Int2 numseg, i, j;

	if(align->segtype == SAS_DENSEG) {
	    DenseSegPtr dsp;
	    dsp = (DenseSegPtr) align->segs;
	    if(dsp) {
		numseg = dsp->numseg;
		for(i = 0, j = numseg-1; i<(numseg+1)/2 && i!=j; i++, --j)
		    {
			Int4 tmp_start,k;
			Uint1 tmp_strand;
			for(k=0;k<dsp->dim;k++) {
			    tmp_start = dsp->starts[2*i+k];
			    dsp->starts[2*i+k] = dsp->starts[2*j+k];
			    dsp->starts[2*j+k] = tmp_start;
			    tmp_strand = dsp->strands[2*i+k];
			    dsp->strands[2*i+k] = dsp->strands[2*j+k];
			    dsp->strands[2*j+k] = tmp_strand;
			}
			dsp->lens[j] = dsp->lens[i];
		    }
	    }
	} else if (align->segtype == SAS_DENDIAG) {
	    DenseDiagPtr ddp,ddp_last,ddp_next;
	    ddp = (DenseDiagPtr) align->segs;
	    ddp_last =NULL;
	    while(ddp) {
		{
		    ddp_next = ddp->next;
		    ddp->next = ddp_last;
		    ddp_last = ddp;
		    ddp = ddp_next;
		}
		align->segs = ddp_last;
	    }
	} else if (align->segtype == SAS_STD) {
	    StdSegPtr ssp,ssp_last,ssp_next;
	    ssp = (StdSegPtr) align->segs;
	    ssp_last =NULL;
	    while(ssp) {
		{
		    ssp_next = ssp->next;
		    ssp->next = ssp_last;
		    ssp_last = ssp;
		    ssp = ssp_next;
		}
		align->segs = ssp_last;
	    }	    
	} else {
	    ErrPostEx(SEV_WARNING,0,0,"SeqAlignReverse order : unsupported seqalign type %d\n",align->segtype);
	}
}

NLM_EXTERN void SeqAlignSwapSeqs(SeqAlignPtr align) {
    if(align->segtype == SAS_DENSEG) {
	DenseSegPtr dsp;
	dsp = (DenseSegPtr) align->segs;
	if(dsp) {
	    SeqIdPtr sip;
	    Int4 i;
	    sip = dsp->ids;
	    dsp->ids = dsp->ids->next;
	    sip->next = NULL;
	    dsp->ids->next = sip;
	    for(i=0;i<dsp->numseg;i++) {
		Int4 tmp_start;
		Uint1 tmp_strand;
		tmp_start = dsp->starts[2*i];
		dsp->starts[2*i] = dsp->starts[2*i+1];
		dsp->starts[2*i+1] = tmp_start;
		tmp_strand = dsp->strands[2*i];
		dsp->strands[2*i] = dsp->strands[2*i+1];
		dsp->strands[2*i+1] = tmp_strand;
	    }
	}
    } else if (align->segtype == SAS_DENDIAG) {
	DenseDiagPtr ddp;
	ddp = (DenseDiagPtr) align->segs;
	while(ddp) {
	    SeqIdPtr sip;
	    sip = ddp->id;
	    ddp->id = ddp->id->next;
	    sip->next = NULL;
	    ddp->id->next = sip;
	    {
		Int4 tmp_start;
		Uint1 tmp_strand;
		tmp_start = ddp->starts[0];
		ddp->starts[0] = ddp->starts[1];
		ddp->starts[1] = tmp_start;
		tmp_strand = ddp->strands[0];
		ddp->strands[0] = ddp->strands[1];
		ddp->strands[1] = tmp_strand;
	    }
	    ddp = ddp->next;
	}		
    } else if (align->segtype == SAS_STD) {
	StdSegPtr ssp;
	ssp = (StdSegPtr) align->segs;
	while(ssp) {
	    SeqIdPtr sip;
	    sip = ssp->ids;
	    if(sip) {
		ssp->ids = ssp->ids->next;
		sip->next = NULL;
		ssp->ids->next = sip;
	    }
	    {
		SeqLocPtr tmp_slp;
		tmp_slp = ssp->loc;
		ssp->loc = ssp->loc->next;
		ssp->loc->next = tmp_slp;
		tmp_slp->next = NULL;
	    }
	    ssp = ssp->next;
	}				
    } else {
	ErrPostEx(SEV_WARNING,0,0,"SeqAlignSwapSeqs: Do not know how to reverse SeqAlign of this type\n");
    }
    return;
}


NLM_EXTERN PSeqAlignInfoPtr SeqAlignToPSeqAlignInfo (SeqAlignPtr sap)
{
        PSeqAlignInfoPtr alip, alip_head, newalip;
        Int2 found;
        SeqAlignPtr alip_sap;
        SeqIdPtr sip;

        if (!sap)
                return NULL;
        alip = (PSeqAlignInfoPtr)MemNew(sizeof(PSeqAlignInfo));
        alip_head = alip;
        sip = SeqIdPtrFromSeqAlign(sap);
        alip->sap = sap;
        alip->sip = SeqIdDupList(sip);
        alip->next = NULL;
        sap = sap->next;
        alip->sap->next=NULL;
        while (sap)
        {
                sip = SeqIdPtrFromSeqAlign(sap);
                found = 0;
                alip = alip_head;
                while (alip && !found)
                {
                        if (SeqIdComp(sip->next, alip->sip->next))
                        {
                                alip_sap = alip->sap;
                                while(alip_sap->next !=NULL)
                                        alip_sap = alip_sap->next;
                                alip_sap->next = sap;
                                sap = sap->next;
                                alip_sap->next->next = NULL;
                                found = 1;
                        }
                        alip = alip->next;
                }
                if (found == 0)
                {
                        alip = alip_head;
                        while(alip->next != NULL)
                                alip = alip->next;
                        newalip = (PSeqAlignInfoPtr)MemNew(sizeof(PSeqAlignInfo));
                        newalip->sip = SeqIdDupList(sip);
                        newalip->sap = sap;
                        newalip->next = NULL;
                        sap = sap->next;
                        newalip->sap->next = NULL;
                        alip->next = newalip;
                }
        }
        return alip_head;
}

NLM_EXTERN SeqAlignPtr ReassembleSeqAlignFromPSeqAlignInfo(PSeqAlignInfoPtr alip)
{
        SeqAlignPtr sap, head_sap, sap_curr;
        PSeqAlignInfoPtr prevalip;

        head_sap = NULL;
        while(alip)
        {
                if (!head_sap)
                {
                        head_sap = alip->sap;
                        sap_curr = head_sap;
                } else
                {
                        sap = alip->sap;
                        sap_curr->next = sap;
                        while (sap_curr->next)
                                sap_curr = sap_curr->next;
                }
                prevalip = alip;
                alip = alip->next;
                prevalip->sip = SeqIdSetFree(prevalip->sip);
                MemFree(prevalip);
        }
        return head_sap;
}

/* Filter SeqAlign by removing bad hits and hits off the main diagonal (as defined by the
   first hit .. which is the highest scoring in blast */
static SeqAlignPtr SeqAlignConsistentDiagFilter(SeqAlignPtr align,SeqLocPtr slp1, SeqLocPtr slp2, 
					 FILE* err_fp,Int4 MaxGap) {
    Int4 old_num,curr_num,len,len1,len2;

    SeqAlignPtr prev,t_align;

    if(slp1 == NULL || slp2 == NULL || align == NULL)
	return NULL;
    len1 = SeqLocLen(slp1);
    len2 = SeqLocLen(slp2);
    if(len1 == 0 || len2 == 0)
	return NULL;

    if(align != NULL) {
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "blast.out", align);
#endif
		/*
		  Checks for harcoded conditions on the FIRST alignment. 
		  Alignment at least
		  50 bases or 90% of length if sequences are smaller than 50.
                  Also do not allow more than 10% mismatches.
		  */
	if(is_bad_blast_alignment(align, slp1, slp2, 50)) {
	    SeqAlignSetFree(align);
	    return NULL;
	}
	    /* Remove hits that are off the main diagonal of the
		   longest hit */
	ModifyAlignList(&align);
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "merge.out", align);
#endif
	if(is_bad_blast_alignment(align, slp1, slp2, 50))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
	if(align != NULL)
	    {
		old_num = SeqAlignCount(align);
		MergeToOneAlignment(&align,MaxGap);
		curr_num = SeqAlignCount(align);
		
		if(align == NULL)
		    {
			fprintf(stderr, "Fail in MergeToOneAlignment\n");
		    } 
		while(curr_num > 1 && curr_num < old_num)
		    {
			old_num = curr_num;
			MergeToOneAlignment(&align,MaxGap);
			curr_num = SeqAlignCount(align);
		    }
		
		if(curr_num > 1)
		    {
			prev = align;
			t_align = align->next;
			while(t_align)
			    {
				len = SeqAlignLength(t_align);
				if(len < 50 || SeqAlignLength(prev)/len > 5)
				    {
					prev->next = NULL;
					SeqAlignSetFree(t_align);
						break;
				    }
				else
				    {
					prev = t_align;
						t_align = t_align->next;
				    }
			    }
		    }
						
		
	    }
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "sim4.out", align);
#endif
	if(is_bad_blast_alignment(align, slp1, slp2, MIN(50, MIN(len1, len2)/5)))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
    }
    return align;
}


static Boolean is_bad_align (BioseqPtr m_bsp, BioseqPtr s_bsp, SeqAlignPtr align, Int4 min_align_len,Int4 loc_len)
{
	DenseSegPtr dsp;
	Int4 m_start, m_stop, s_start, s_stop;
	SeqPortPtr m_spp, s_spp;
	Int2  i;
	Uint1 m_res, s_res;
	Int4 mismatch = 0;
	Uint1 code;
	Int4 j;
	Uint1 strand;
	Int4 align_len;

	code = Seq_code_iupacna;
	dsp = (DenseSegPtr) align->segs;
	SeqAlignStartStopById(align, dsp->ids->next, &s_start, &s_stop, &strand);
	if(s_start<0 || s_stop>= s_bsp->length) {
	    ErrPostEx(SEV_WARNING,0,0,"find_mismatch_residue: corrupted alignment \n");
	    return TRUE;
	} else {
	    s_spp = SeqPortNew(s_bsp, s_start, s_stop, dsp->strands[1], code);
	}
	/* if(s_stop - s_start + 1 < MIN(min_align_len, loc_len-MAX(10, loc_len/10)))
		return TRUE; */

	SeqAlignStartStopById(align, dsp->ids, &m_start, &m_stop, &strand);
	if(m_start<0 || m_stop>= m_bsp->length) {
	    ErrPostEx(SEV_WARNING,0,0,"find_mismatch_residue: corrupted alignment \n");
	    return TRUE;
	} else {

	    m_spp = SeqPortNew(m_bsp, m_start, m_stop, strand, code);
	}


	align_len = 0;
	for(i = 0; i<dsp->numseg; ++i)
	{
		m_start = dsp->starts[2*i];
		s_start = dsp->starts[2*i + 1];
		if(m_start == -1)
		{
			if(dsp->lens[i] <=5)
			{
				mismatch += 1;
				++align_len;
			}
			if(s_start != -1)
				SeqPortSeek(s_spp, dsp->lens[i], SEEK_CUR);
		}
		else if(s_start == -1)
		{
			if(dsp->lens[i] <=5)
			{
				mismatch += 1;
				++align_len;
			}
			if(m_start != -1)
				SeqPortSeek(m_spp, dsp->lens[i], SEEK_CUR);
		}
		else
		{
			for(j = 0; j<dsp->lens[i]; ++j)
			{
				m_res = SeqPortGetResidue(m_spp);
				while(!IS_ALPHA(m_res) && m_res != SEQPORT_EOF)
					m_res = SeqPortGetResidue(m_spp);
				s_res = SeqPortGetResidue(s_spp);
				while(!IS_ALPHA(s_res) && s_res != SEQPORT_EOF)
					s_res = SeqPortGetResidue(s_spp);
				if(m_res != 'n' && m_res != 'N' && s_res != 'n' && s_res != 'N')
				{
					if(m_res != s_res)
						++mismatch;
				}
			}
			align_len += dsp->lens[i];
		}
		
	}

	SeqPortFree(m_spp);
	SeqPortFree(s_spp);
	if(align_len < MIN(min_align_len, loc_len-MAX(10, loc_len/10)))
		return TRUE;
	if(mismatch > 0)
		return ((mismatch*100)/align_len >10);
	else
		return FALSE;
}

static Boolean need_recompute(SeqAlignPtr align, Int4 m_len, Int4 s_len)
{
        DenseSegPtr dsp;
        Int4 m_start, m_stop;
        Int4 s_start, s_stop;
        Int2 end, i;
        Int2 small_seg_num;

        if(m_len > 10000 || s_len > 10000)
                return FALSE;
        dsp = (DenseSegPtr) align->segs;
        end = dsp->numseg-1;
        m_start = dsp->starts[0];
        m_stop = dsp->starts[2*end] + dsp->lens[end] -1;

        if(dsp->strands[1] == Seq_strand_minus)
        {
                s_start = dsp->starts[2*end +1];
                s_stop = dsp->starts[1] + dsp->lens[0] -1;
                if(s_start > MAX_HANG && (m_len - 1 - m_stop) > MAX_HANG)
                        return TRUE;
                if((s_len-1-s_stop)>MAX_HANG && m_start > MAX_HANG)
                        return TRUE;
        }
        else
        {
                s_start = dsp->starts[1];
                s_stop = dsp->starts[2*end+1] + dsp->lens[end] -1;
                if(s_start > MAX_HANG && m_start > MAX_HANG)
                        return TRUE;
                if((s_len-1-s_stop) > MAX_HANG && (m_len-1-m_stop)> MAX_HANG)
                        return TRUE;
        }

        small_seg_num = 0;
        for(i = 0; i<dsp->numseg; ++i)
        {
                if(dsp->lens[i] < 10)
                        ++small_seg_num;
        }

        return (small_seg_num >= 20);

}

static void extend_dsp_ends(Int4 first_len, Int4 second_len, SeqAlignPtr align, Int4 max_overhang)
{
	Int4 first, last;
	DenseSegPtr dsp;
	Int2 endseg;
	Int4 extend_len;
	Int4 offset;

	dsp = (DenseSegPtr) align->segs;
	endseg = dsp->numseg-1;
	if(dsp->strands[1] == Seq_strand_minus)
	{
		first = dsp->starts[2*endseg+1];
		if(first > 0 && first < max_overhang)
		{
                    if(dsp->strands[0] != Seq_strand_minus)
                        extend_len = MIN(first, first_len - dsp->starts[2*endseg] - dsp->lens[endseg]);
                    else
                        extend_len = MIN(first, dsp->starts[2*endseg]);
			if(extend_len > 0)
			{
				dsp->lens[endseg] += extend_len;
                                dsp->starts[2*endseg+1]-=extend_len;
                                if(dsp->strands[0]==Seq_strand_minus)
                                    dsp->starts[2*endseg] -= extend_len;
			}
		}

		last = dsp->starts[1] + dsp->lens[0];
		offset = second_len - last;
		if(offset > 0 && offset < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus) 
			extend_len = MIN(offset,dsp->starts[0]);
                    else
                        extend_len = MIN(offset,first_len - dsp->starts[0]-dsp->lens[0]);
			if(extend_len > 0)
			{
				dsp->lens[0] += extend_len;
                                if(dsp->strands[0]==Seq_strand_plus)
                                    dsp->starts[0]-=extend_len;
			}
		}
	}
	else
	{
		first = dsp->starts[1];
		if(first> 0 && first < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus) 
			extend_len = MIN(first, dsp->starts[0]);
                    else
                        extend_len = MIN(first, first_len -dsp->starts[0]-dsp->lens[0]);
			if(extend_len > 0)
			{
				dsp->lens[0] += extend_len;
				dsp->starts[1] -= extend_len;
                                if(dsp->strands[0]!=Seq_strand_minus)
                                    dsp->starts[0] -= extend_len;
			}
		}
		last = dsp->starts[2*endseg+ 1] + dsp->lens[endseg];
		offset = second_len  - last;

		if(offset > 0 && offset < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus)
			extend_len = MIN(offset, first_len - dsp->starts[2*endseg] - dsp->lens[endseg]);
                    else
			extend_len = MIN(offset, dsp->starts[2*endseg]);
			if(extend_len > 0)
			{
				dsp->lens[endseg] += extend_len;
                                if(dsp->strands[0]==Seq_strand_minus)
                                    dsp->starts[2*endseg]-=extend_len;
			}
		}
	}

}
	


/* Tries to Make a global alignment from a Set of local SeqAligns. 
 */
NLM_EXTERN SeqAlignPtr SeqAlignSetGlobalFromLocal(SeqAlignPtr align,SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *err_fp,Int4 MaxGap)
{
	Int4 len_1, len_2;
	Char label[101];
        if(align==NULL || loc_1 == NULL || loc_2 == NULL || bsp_1 == NULL || bsp_2 == NULL)
            return NULL;
	len_1 = SeqLocLen(loc_1);
	if(len_1 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_1), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long)len_1, label);
            }
            return NULL;
	}
	len_2 = SeqLocLen(loc_2);
	if(len_2 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_2), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long) len_1, label);
            }
            return NULL;
	}

	if(align != NULL)
	{

	    align = SeqAlignConsistentDiagFilter(align,loc_1, loc_2,err_fp,MaxGap);
	}

	if(align != NULL)
	{
		extend_dsp_ends(bsp_1->length, bsp_2->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
		extend_dsp_ends(bsp_2->length, bsp_1->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
	}
	return align;
}


 
/********************************************************
***
*** NormalizeSeqAlignId
***   Checks local seqid . if a seqid string contains "acc"
***   seqid has a correspondant sequence in db.
***   This local seqid is replaced by its gi number. 
***
***   The local sequence is compared to the sequence from 
***   the db. If the local sequence is a region of the db sequence
***   the positions in the seqalign are updaded with the offset.
***
***          Thanks to Mark for this useful function! 
**********************************************************/
static Int4 getlengthforid (SeqIdPtr sip)
{
  BioseqPtr        bsp;
  Int4             lens=0;
 
  if (sip==NULL)
     return 0;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     lens = bsp->length;
     BioseqUnlock (bsp);
  }
  return lens;
}

static void showtextalign_fromalign (SeqAlignPtr salp, CharPtr path, FILE *fp)
{
  SeqAnnotPtr sap;
  Int4        line = 80;     
  Uint4       option = 0;
  Boolean     do_close = TRUE;

  if (salp == NULL || (path==NULL && fp==NULL))
     return;
  if (path != NULL && fp == NULL) {
     fp = FileOpen (path, "a");
  }
  else
     do_close = FALSE;
  if (fp != NULL) {
/********
{{
     option = RULER_TOP;
     option|=DISPE_SHOWBLOCK;
     option|=VIEW_VARIA;
     option|=DISP_FULL_TXT;
     DDV_DisplayDefaultAlign (salp, 0, -1,-1, option, NULL, fp); 
}}
*******/
     sap = SeqAnnotNew ();
     sap->type = 2;
     sap->data = (Pointer) salp;
     option += TXALIGN_MASTER;
     option += TXALIGN_MISMATCH;
     ShowTextAlignFromAnnot (sap, line, fp, NULL, NULL, option, NULL, NULL, NULL);
     sap->data = NULL;
     SeqAnnotFree (sap);
     if (do_close) {
        FileClose(fp);
     }
  } 
}

typedef struct besthit {
   SeqAlignPtr  sap;
   SeqIdPtr     sip1;
   SeqIdPtr     sip2;
   BioseqPtr    bsp1;
   BioseqPtr    bsp2;
  CharPtr      errstr;
   Int4         nonly;
   Int4         errtype;
    struct besthit PNTR next;
} BestHit, PNTR BestHitPtr;

static int LIBCALLBACK OrderBestHits(Pointer ptr1, Pointer ptr2)
{
   BestHitPtr  hip1;
   BestHitPtr  hip2;

   hip1 = *((BestHitPtr PNTR)(ptr1));
   hip2 = *((BestHitPtr PNTR)(ptr2));
   if (hip1->errtype < hip2->errtype)
      return -1;
   else if (hip1->errtype > hip2->errtype)
      return 1;
   else
      return 0;
}

NLM_EXTERN Boolean TruncateAlignment (SeqAlignPtr salp, Int4 num_aln_pos, Boolean from_left)
{
  DenseSegPtr dsp, dsp_new;
  Int4        num_trim, seg_lens;
  Int4        segment_num;
  Int4        sequence_num;
  Int4        segment_trim = 0, orig_segment_offset;
  
  if (salp == NULL 
      || salp->segtype != 2 
      || salp->segs == NULL
      || num_aln_pos <=0) return FALSE;
      
  dsp = (DenseSegPtr) salp->segs;
  
  /* how many segments do we need to trim? */
  num_trim = 0;
  seg_lens = 0;
  if (from_left)   
  {
    while (seg_lens < num_aln_pos && num_trim < dsp->numseg) 
    {
      seg_lens += dsp->lens[num_trim];
      num_trim++;
    }
    if (seg_lens > num_aln_pos)
    {
      num_trim--;
      segment_trim = num_aln_pos - seg_lens + dsp->lens [num_trim];
    }
  }
  else
  {
    while (seg_lens < num_aln_pos && num_trim < dsp->numseg)
    {
      seg_lens += dsp->lens[dsp->numseg - 1 - num_trim];
      num_trim++;
    }
    if (seg_lens > num_aln_pos)
    {
      num_trim--;
      segment_trim = num_aln_pos - seg_lens + dsp->lens [dsp->numseg - 1 - num_trim];
    }
  }
  
  
  if (num_trim == dsp->numseg)
  {
    return FALSE;
  }
  
  dsp_new = DenseSegNew();
  dsp_new->dim = dsp->dim;
  dsp_new->ids = SeqIdDupList (dsp->ids);
  dsp_new->numseg = dsp->numseg - num_trim;
  dsp_new->starts = (Int4Ptr) MemNew (dsp_new->dim * dsp_new->numseg * sizeof (Int4));
  dsp_new->strands = (Uint1Ptr) MemNew (dsp_new->dim * dsp_new->numseg * sizeof (Uint1));  
  dsp_new->lens = (Int4Ptr) MemNew (dsp_new->numseg * sizeof (Int4));
  if (dsp->scores != NULL)
  {
    dsp_new->scores = (ScorePtr) MemNew (dsp_new->numseg * sizeof (Score));
  }
  dsp_new->strands = (Uint1Ptr) MemNew (dsp_new->dim * dsp_new->numseg * sizeof (Uint1));
  
  segment_num = 0;
  orig_segment_offset = 0;
  
  if (from_left)
  {
    /* remove segments from the left */
    orig_segment_offset = num_trim;
    /* if we're only getting rid of part of a segment, 
     * need to adjust the starts and the length
     */
    if (segment_trim > 0)
    {
      for (sequence_num = 0; sequence_num < dsp_new->dim; sequence_num++)
      {
        /* copy strands */
        dsp_new->strands[segment_num * dsp_new->dim + sequence_num]
              = dsp->strands[(orig_segment_offset + segment_num) * dsp->dim + sequence_num];
        /* adjust starts */
        if (dsp->starts[(orig_segment_offset + segment_num) * dsp->dim + sequence_num] > -1)
        {
          if (dsp_new->strands[segment_num * dsp_new->dim + sequence_num] == Seq_strand_minus)
          {
            dsp_new->starts[segment_num * dsp_new->dim + sequence_num]
                = dsp->starts[(orig_segment_offset + segment_num) * dsp->dim + sequence_num];
          }
          else
          {
            dsp_new->starts[segment_num * dsp_new->dim + sequence_num]
                = dsp->starts[(orig_segment_offset + segment_num) * dsp->dim + sequence_num] + segment_trim;
          }
        }
        else
        {
          dsp_new->starts[segment_num * dsp_new->dim + sequence_num] = -1;
        }
        if (dsp_new->scores != NULL)
        {
          /* copy scores */
          dsp_new->scores [segment_num * dsp_new->dim + sequence_num] = 
              *((ScorePtr)AsnIoMemCopy (&(dsp->scores[(num_trim + segment_num) * dsp->dim + sequence_num]),
                            (AsnReadFunc) ScoreSetAsnRead, 
                            (AsnWriteFunc)ScoreSetAsnWrite));          
        }
      }
      /* adjust len */
      dsp_new->lens [segment_num] = dsp->lens[num_trim + segment_num] - segment_trim;
      segment_num++;
    }
  }
  
  /* copy segments to the end of the alignment*/
  while (segment_num < dsp_new->numseg)
  {
    for (sequence_num = 0; sequence_num < dsp_new->dim; sequence_num++)
    {
      /* copy starts */
      dsp_new->starts[segment_num * dsp_new->dim + sequence_num]
            = dsp->starts[(orig_segment_offset + segment_num) * dsp->dim + sequence_num];
      /* copy strands */
      dsp_new->strands[segment_num * dsp_new->dim + sequence_num]
            = dsp->strands[(orig_segment_offset + segment_num) * dsp->dim + sequence_num];
      if (dsp_new->scores != NULL)
      {
        /* copy scores */
        dsp_new->scores [segment_num * dsp_new->dim + sequence_num] = 
            *((ScorePtr)AsnIoMemCopy (&(dsp->scores[(orig_segment_offset + segment_num) * dsp->dim + sequence_num]),
                          (AsnReadFunc) ScoreSetAsnRead, 
                          (AsnWriteFunc)ScoreSetAsnWrite));          
      }
    }
    /* copy len */
    dsp_new->lens [segment_num] = dsp->lens[orig_segment_offset + segment_num];
    segment_num++;
  }
  
  if (!from_left && segment_trim > 0)
  {
    /* adjust len and starts for last segment if necessary */
    segment_num--;
    /* adjust length of segment */
    dsp_new->lens [dsp_new->numseg - 1] -= segment_trim;
    /* adjust starts if not in gap and on minus strand */
    for (sequence_num = 0; sequence_num < dsp_new->dim; sequence_num++)
    {
      if (dsp_new->starts[segment_num * dsp_new->dim + sequence_num] > -1
          && dsp_new->strands[segment_num * dsp_new->dim + sequence_num] == Seq_strand_minus)
      {
        dsp_new->starts[segment_num * dsp_new->dim + sequence_num] += segment_trim;
      }
    }
  }

  /* replace old segments with new */  
  dsp = DenseSegFree (dsp);
  salp->segs = dsp_new;

  /* redo index */
  AMFreeAllIndexes (salp);
  AlnMgr2IndexSingleChildSeqAlign (salp);
      
  return TRUE;
}

static Boolean CutAlignmentAtGaps (SeqAlignPtr salp, Int4 seq_num)
{
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  DeltaSeqPtr dsp;
  SeqLitPtr   slip;
  SeqLocPtr   loc;
  Int4        seq_offset, aln_len;
  SeqAlignPtr salp_next, new_salp_mid;
  Int4        cut_pos;
  Int4        start, stop;
  
  if (salp == NULL || salp->segtype != 2 || salp->segs == NULL || seq_num < 1)
  {
    return FALSE;
  }
  
  sip = AlnMgrGetNthSeqIdPtr(salp, seq_num);

  bsp = BioseqFind (sip);
  if (bsp == NULL)
  {
    return FALSE;
  }
  
  if (bsp->repr != Seq_repr_delta)
  {
    return TRUE;
  }
  
  salp_next = salp->next;
	salp->next = NULL;
	
  AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  seq_offset = 0;
  while (dsp != NULL)
  {
    if (dsp->choice == 1)
    {
      loc = (SeqLocPtr) dsp->data.ptrvalue;
      seq_offset += SeqLocLen (loc);
    }
    else if (dsp->choice == 2)
    {
      slip = (SeqLitPtr) dsp->data.ptrvalue;
		  if (slip->seq_data == NULL && slip->fuzz != NULL && slip->fuzz->choice == 4
		      && start < seq_offset + 100 && stop > seq_offset)
		  {
        if (start < seq_offset)
        {
          aln_len = SeqAlignLength (salp);
          /* make a copy of the original alignment */
		      new_salp_mid = (SeqAlignPtr) AsnIoMemCopy (salp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);   
		      AlnMgr2IndexSingleChildSeqAlign (new_salp_mid);
          
          /* remove right end of first alignment, after start of gap */
		      cut_pos = aln_len - AlnMgr2MapBioseqToSeqAlign (salp, seq_offset, seq_num);
		      TruncateAlignment (salp, cut_pos, FALSE);
          AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
          
          /* remove left end of second alignment, before start of gap */
		      cut_pos = AlnMgr2MapBioseqToSeqAlign (new_salp_mid, seq_offset, seq_num);
          TruncateAlignment (new_salp_mid, cut_pos, TRUE);
          salp->next = new_salp_mid;
          salp = new_salp_mid;
          AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
        }
        
        if (stop > seq_offset + 100)
        {
          aln_len = SeqAlignLength (salp);
          /* make a copy of the original alignment */
		      new_salp_mid = (SeqAlignPtr) AsnIoMemCopy (salp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);   
		      AlnMgr2IndexSingleChildSeqAlign (new_salp_mid);
		      /* remove right end of first alignment, after end of gap */
		      cut_pos = aln_len - AlnMgr2MapBioseqToSeqAlign (salp, seq_offset + 100, seq_num);
		      TruncateAlignment (salp, cut_pos, FALSE);
          AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
		      
		      /* remove left end of second alignment, before end of gap */
		      cut_pos = AlnMgr2MapBioseqToSeqAlign (new_salp_mid, seq_offset + 100, seq_num);
          TruncateAlignment (new_salp_mid, cut_pos, TRUE);
          salp->next = new_salp_mid;
          salp = new_salp_mid;
          AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
        }
		  }
		  seq_offset += slip->length;
    }
    dsp = dsp->next;
  }
  salp->next = salp_next;
  
  if (salp->next != NULL)
  {
    return CutAlignmentAtGaps (salp->next, seq_num);
  }
  return TRUE;
}


static SeqAlignPtr RemoveSequencesInGapFromAlignments (SeqAlignPtr salp)
{
  SeqAlignPtr salp_next;
  SeqIdPtr    sip, sip_list, sip_remove_list = NULL, sip_tmp;
  Boolean     remove_this;
  Int4        seq_num, start, stop, seq_offset;
  DeltaSeqPtr dsp;
  SeqLitPtr   slip;
  SeqLocPtr   loc;
  Int4        num_keep = 0;
  Int4        num_all_gap = 0;
  BioseqPtr   bsp;
  
  if (salp == NULL)
  {
    return NULL;
  }
  
  salp_next = salp->next;
  salp->next = NULL;

  sip_list = SeqIdPtrFromSeqAlign(salp);
  
  for (sip = sip_list, seq_num = 1; sip != NULL; sip = sip->next, seq_num++)
  {
    AlnMgr2GetNthSeqRangeInSA(salp, seq_num, &start, &stop);
    if (start == -2 && stop == -2)
    {
      num_all_gap++;
    }
    

    bsp = BioseqFind (sip);
    if (bsp == NULL || bsp->repr != Seq_repr_delta)
    {
      num_keep ++;
      continue;
    }

    remove_this = FALSE;
    dsp = (DeltaSeqPtr) bsp->seq_ext;
    seq_offset = 0;
    while (dsp != NULL && !remove_this && seq_offset < stop)
    {
      if (dsp->choice == 1)
      {
        loc = (SeqLocPtr) dsp->data.ptrvalue;
        seq_offset += SeqLocLen (loc);
      }
      else if (dsp->choice == 2)
      {
        slip = (SeqLitPtr) dsp->data.ptrvalue;
		    if (slip->seq_data == NULL && slip->fuzz != NULL && slip->fuzz->choice == 4
		        && seq_offset <= start && seq_offset + 100 >= stop)
		    {
		      remove_this = TRUE;
		    }
		    seq_offset += slip->length;
      }
      dsp = dsp->next;
    }
    if (remove_this)
    {
      sip_tmp = SeqIdDup (sip);
      sip_tmp->next = sip_remove_list;
      sip_remove_list = sip_tmp;
    }
    else
    {
      num_keep++;
    }
  }
  
  if (num_keep < 2 || num_keep == num_all_gap)
  {
    salp = SeqAlignFree (salp);
    salp = RemoveSequencesInGapFromAlignments (salp_next);
    sip_remove_list = SeqIdFree (sip_remove_list);
  }
  else
  {
    for (sip = sip_remove_list; sip != NULL; sip = sip_tmp)
    {
      sip_tmp = sip->next;
      sip->next = NULL;
      
      salp = SeqAlignBioseqDeleteById (salp, sip);
      /* reindex after removing sequence */
      AlnMgr2IndexSingleChildSeqAlign (salp);
      
      sip = SeqIdFree (sip);
    }
    salp->next = salp_next;
    salp->next = RemoveSequencesInGapFromAlignments (salp->next);
  }
  
  return salp;
}


NLM_EXTERN SeqAlignPtr MakeDiscontiguousAlignments (SeqAlignPtr salp)
{
  Int4        seq_num;
  
  if (salp == NULL || salp->segtype != 2 || salp->segs == NULL)
  {
    return salp;
  }
  
  for (seq_num = 1; seq_num < salp->dim; seq_num++)
  {
    CutAlignmentAtGaps (salp, seq_num);
  }
  
  salp = RemoveSequencesInGapFromAlignments (salp);

  return salp;
}



