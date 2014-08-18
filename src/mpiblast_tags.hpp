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
#ifndef __MPIBLAST_TAGS_HPP
#define __MPIBLAST_TAGS_HPP

// Tags for each type of MPI communication
// worker -> master
#define	EVENT_TYPE			1<<0	/* The type of events to the master */
//    scheduler tags
#define	DB_FRAGMENTS_TAG	1<<1	/* A list of database fragments */
#define	BLAST_RETCODE_TAG	1<<2	/* A BLAST return code */
#define	WORKER_IDLE			1<<3	/* That a worker is idle */
#define	FRAGMENT_COPY_COMPLETE	1<<4	/* That a fragment copy has completed */
#define	QUERY_TAG				1<<5	/* The query indices that should be searched */
#define DB_DATATOSEND_TAG       1<<6    /* Amount of Data to send when copying DB fragment */
#define DB_EXTENSION_DATA_TAG   1<<7      /* Actual data sent when copying DB fragment */
#define DB_FILESIZE_TAG         1<<8      /* Total amount of data to receive when copying DB fragment */
#define PROFILE_TAG				1<<9      
#define QUERY_SEGMENT_REQUEST	1<<10      
#define GROUP_FINISHED			1<<11      

#define SCHEDULER_TAGS		0|WORKER_IDLE|FRAGMENT_COPY_COMPLETE

//    writer tags
#define	BLAST_RESULTS_TAG		1<<16	/* Actual BLAST results */
#define PARTIAL_OUTPUTS_DATA 	1<<17
#define OUTPUT_OFFSETS			1<<18
#define DATA_TO_WRITE 			1<<19	
#define META_DATA_TO_WRITE 		1<<20
#define WRITER_FINALIZED		1<<21
#define QUERY_WRITE_ACK 		1<<22 /* ackknowledgement of writing a query */
#define WRITE_REQUEST			1<<23
#define WRITE_COMPLETE			1<<24
#define COLLECTIVE_READY		1<<25 /* telling writer master that I'm ready for collective write */
#define WAIT_WRITE_COMPLETE		1<<26 /* telling writer master  that I'm waiting for a write to complete*/

#define WRITER_TAGS		0|PARTIAL_OUTPUTS_DATA|META_DATA_TO_WRITE|WRITER_FINALIZED|QUERY_WRITE_ACK|WRITE_REQUEST|WRITE_COMPLETE|COLLECTIVE_READY|WAIT_WRITE_COMPLETE

// master -> worker
// assignment types
#define	ASSIGNMENT_TYPE		0	/* The type of worker assignment given */
#define	WRITER_EVENT 		1	/* writer events */
#define SEARCH_COMPLETE		2 /* Code sent to a worker when the search is complete */ 
#define COPY_FRAGMENT		3 /* Code indicating that the worker should copy a fragment */
#define SEARCH_FRAGMENT		4 /* Code telling a worker to search a fragment */
#define WORKER_QUIT		5 /* Code telling a worker to quit */
#define DO_NOTHING		6 /* place holder telling a worker to do nothing */
#define SEND_RESULTS	7 /* An assignment from a writer to send results */
#define WORKER_WRITE	8 /* Code telling a worker to print results */
#define PRINT_RESULTS	9 /* Code telling a worker to print results */
#define CHECK_RECEIVED_OFFSETS 10 /* Code telling a worker to check if the offsets of current query is processed */
#define	SHIFT_WRITE_LEADER	11	/* shift writer leader */
#define	WRITER_QUIT			12	/* tell writer worker to quit */
#define SEARCH_SEGMENT		13  // tell group master to search a segment
#define NO_MORE_SEGMENT		14  // tell group master that no more segment left
#define REQUEST_FRAGMENT_COPY	15
#define CONFIRM_COPY_FRAGMENT 	16  
#define QUERY_DATA			17

// communication tags
#define WRITER_MASTER_END 		33	/* indicate writer master has finished */

#endif
