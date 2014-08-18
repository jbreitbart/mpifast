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
#ifndef __distributed_bioseq_h__
#define __distributed_bioseq_h__

#ifdef __cplusplus
extern "C"{
#endif

#include "ncbi.h"
#include "objseq.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"

#include "mpi.h"
#include "blast_hooks.h"

/**
 * MPI Blast distributed database query functions
 */
Boolean Enable_MPI_Distributed_DB();
Boolean Disable_MPI_Distributed_DB();

extern Boolean debug_bsfetch;
extern Boolean debug_bslookup;


/** a linked list of cached BioseqPtrs */
extern ValNodePtr bs_cache;

/** a linked list of pointers to buffers sent using MPI_Isend */
extern ValNodePtr send_buffers;

/** a linked list of MPI_Request structs corresponding to the send_buffers list */
extern ValNodePtr send_requests;

/** Add a send buffer and corresponding request to the linked lists */
void addSendBuffer( void* buffer, MPI_Request* request );
/** Free any send buffers that have completed sending */
void freeCompletedSendBuffers();
/** Wait for any send buffers to complete sending */
void waitForSendCompletions();

/**
 * Sets the linked-list of BioseqPtrs that will be used during BSLookups from the cache
 * @param cache_vnp	The linked list of BioseqPtrs
 */
void setBioseqCache( ValNodePtr cache_vnp );


#ifdef __cplusplus
}
#endif

#endif /* __distributed_bioseq_h__ */
