/*************************************************************************
** Copyright 2009 by Virginia Polytechnic Institute and State
** University. All rights reserved. Virginia Polytechnic Institute and
** State University (Virginia Tech) owns the mpiBLAST software and its
** associated documentation ("Software"). You should carefully read the
** following terms and conditions before using this software. Your use
** of this Software indicates your acceptance of this license agreement
** and all terms and conditions.
** 
** The parallel input/output (PIO) portion of this work was based in
** part on published research that was supported by the North Carolina
** State University and Oak Ridge National Laboratory. The actual
** implementation was completed at Virginia Tech in the summer of 2007.
** 
** Additionally, portions of this work are derivative works of the NCBI
** C Toolkit. Although, the NCBI C Toolkit is released freely without
** restriction under the terms of their license, the following files
** listed below, have been modified by Virginia Tech, and as such are
** redistributed under the terms of this license.
** 
** ncbi/api/txalign.c
** ncbi/corelib/ncbifile.c
** ncbi/object/objalign.c
** ncbi/tools/blast.c
** ncbi/tools/blastdef.h
** ncbi/tools/blastool.c
** ncbi/tools/blastutl.c
** ncbi/tools/blfmtutl.c
** ncbi/tools/ncbisam.c
** ncbi/tools/readdb.c
** ncbi/tools/readdb.h
** 
** License:
** 
** This file is part of mpiBLAST.
** 
** mpiBLAST is free software: you can redistribute it and/or modify it 
** under the terms of the GNU General Public License version 2 as published 
** by the Free Software Foundation. 
** 
** Accordingly, mpiBLAST is distributed in the hope that it will be
** useful, but WITHOUT ANY WARRANTY; without even the implied warranty
** of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
** General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with mpiBLAST. If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/
/*
 * Written by J.D. Gans and Aaron Darling
 */
#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <xmlblast.h>
#include <mblast.h>
#ifdef BLAST_CS_API
#include <objblst3.h>
#include <netblap3.h>
#endif
#include <asn.h>

#include "mpi_util.h"
#include "distributed_bioseq.h"
#include "blast_hooks.h"

BSFetchTop bsfetch_orig = NULL;
Boolean debug_bsfetch = FALSE;
Boolean debug_bslookup = FALSE;


/* a linked list of cached BioseqPtrs */
ValNodePtr bs_cache;

void setBioseqCache( ValNodePtr cache_vnp )
{
	bs_cache = cache_vnp;
}

/**
 * Retrieve a bioseq by SeqId from an mpiBLAST-specific local cache of Bioseqs
 */
BioseqPtr LIBCALLBACK mpiBLAST_BSCacheFetchFunc(SeqIdPtr sid, Uint1 ld_type)
{
	char buf[255];
	ValNodePtr cur_bs = bs_cache;
	SeqIdPrint( sid, buf, PRINTID_FASTA_LONG );
	/*if( debug_bsfetch )fprintf( stderr, "Fetching %s\n", buf );*/
	while( cur_bs != NULL ){
//		SeqIdPrint( ((BioseqPtr)cur_bs->data.ptrvalue)->id, buf2, PRINTID_FASTA_LONG );
		if( cur_bs->data.ptrvalue == NULL )
			fprintf( stderr, "Error:  Cached bioseq is NULL\n" );
		else if( ((BioseqPtr)cur_bs->data.ptrvalue)->id == NULL )
			fprintf( stderr, "Error:  Cached bioseq has NULL SeqIdPtr.\n" );
		else if( SeqIdComp( ((BioseqPtr)cur_bs->data.ptrvalue)->id, sid ) == SIC_YES ){
			/* this is the one we're looking for! */
			return (BioseqPtr)cur_bs->data.ptrvalue;
		}
		cur_bs = cur_bs->next;
	}

	/* unable to find the Bioseq in the cache */
	fprintf(stderr, "mpiBLAST_BSCacheFetchFunc: Error finding SeqId %s in cache!\n", buf);
	return NULL;
}


/* Use this function to enable distributed database searching with MPI */
Boolean Enable_MPI_Distributed_DB()
{
	SeqMgrPtr smp;
	
	smp = SeqMgrWriteLock();
	
	if(smp == NULL){
		return FALSE;
	}
	
	/* Save the original lookup function */
	bsfetch_orig = smp->bsfetch;
	smp->bsfetch = mpiBLAST_BSCacheFetchFunc;
	return SeqMgrUnlock();
}

/* Use this function to disable distributed database searching with MPI */
Boolean Disable_MPI_Distributed_DB()
{
	SeqMgrPtr smp = SeqMgrWriteLock();
	
	if(smp == NULL)
		return FALSE;
	
	if(bsfetch_orig == NULL)
		return FALSE;
	
	/* Restore the default fetch function -- see NCBI's seqmgr.c for details */
	smp->bsfetch = bsfetch_orig;

	return SeqMgrUnlock();
}




/** a linked list of pointers to buffers sent using MPI_Isend */
ValNodePtr send_buffers = NULL;

/** a linked list of MPI_Request structs corresponding to the send_buffers list */
ValNodePtr send_requests = NULL;

/** Add a send buffer and corresponding request to the linked lists */
void addSendBuffer( void* buffer, MPI_Request* request ) {
	ValNodeAddPointer( &send_buffers, 0, buffer );
	ValNodeAddPointer( &send_requests, 0, request );
}


/** Free any send buffers that have completed sending */
void freeCompletedSendBuffers() {
	ValNodePtr cur_buf = send_buffers;
	ValNodePtr cur_req = send_requests;
	ValNodePtr prev_buf = NULL;
	ValNodePtr prev_req = NULL;
	int counter = 0;
	int freed = 0;
	
	/* scan the entire list of send buffers to search for
	   requests that have completed */
	while( cur_buf != NULL ){
		int flag;
		MPI_Status status;
		ValNodePtr del_buf = NULL;
		ValNodePtr del_req = NULL;

/*		fprintf( stderr, "checking req %d\n", counter ); */

		MPI_Test( (MPI_Request*)cur_req->data.ptrvalue, &flag, &status );
		if( flag == TRUE ){
/*			fprintf( stderr, "Send %d has completed, freeing data\n", counter ); */
			/* this send has completed.  free its buffer */
			meta_MPI_Free_mem( cur_buf->data.ptrvalue );
			/* free the request record */
			free( cur_req->data.ptrvalue );

/*			fprintf( stderr, "unlinking records\n" ); */
			
			/* unlink records from lists */
			if( prev_buf )
				prev_buf->next = cur_buf->next;
			if( prev_req )
				prev_req->next = cur_req->next;
			del_buf = cur_buf;
			del_req = cur_req;
			cur_buf = cur_buf->next;
			cur_req = cur_req->next;
			/* prev pointers don't advance */
			
			/* check whether global record pointers should be updated */
			if( del_buf == send_buffers ){
				send_buffers = cur_buf;
				send_requests = cur_req;
			}
			freed++;
/*			fprintf( stderr, "freeing records\n" ); */
			/* free records */
			Nlm_Free( del_buf );
			Nlm_Free( del_req );
		}else{
			/* advance to next record */
			prev_buf = cur_buf;
			prev_req = cur_req;
			cur_buf = cur_buf->next;
			cur_req = cur_req->next;
		}
		counter++;
	}

/*	if( debug_bsfetch )
		fprintf( stderr, "%d / %d pending send operations completed\n", freed, counter );
*/
}

/** Wait for any send buffers to complete sending */
void waitForSendCompletions() {
        ValNodePtr cur_buf = send_buffers;
        ValNodePtr cur_req = send_requests;

		/* wait for all requests to complete */
        while( cur_req != NULL ){
                MPI_Status status;

                MPI_Wait( (MPI_Request*)cur_req->data.ptrvalue, &status );

                /* this send has completed.  free its buffer */
                meta_MPI_Free_mem( cur_buf->data.ptrvalue );
                /* free the request record */
                free( cur_req->data.ptrvalue );

                cur_req = cur_req->next;
                cur_buf = cur_buf->next;
        }

        /* free the ValNodeList */
        ValNodeFree( send_buffers );
        ValNodeFree( send_requests );
		send_buffers = NULL;
		send_requests = NULL;
}

/* #endif */ /* USING_MPI */
