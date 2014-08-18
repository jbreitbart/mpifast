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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi_util.h"
#include <stdlib.h>

/* a simple wrapper around MPE_Log_event that consolidates
 * all the #ifdefs in one place
 */
void meta_MPE_Log_event( int a, int b, const char* c )
{
#ifdef MPE
	MPE_Log_event( a, b, c );
#endif
}

/* a simple wrapper around meta_MPI_Alloc_mem that calls malloc() if
 * meta_MPI_Alloc_mem isn't available
 */
int meta_MPI_Alloc_mem( MPI_Aint size, MPI_Info info, void* baseptr )
{
#ifdef HAVE_MPI_ALLOC_MEM
	return MPI_Alloc_mem( size, info, baseptr );
#else
	*(void**)baseptr = malloc(size);
	return 0;
#endif
}

int meta_MPI_Free_mem( void* base )
{
#ifdef HAVE_MPI_ALLOC_MEM
	return MPI_Free_mem(base);
#else
	free(base);
	return 0;
#endif
}

double meta_MPI_Wtime() {
	return MPI_Wtime();
}
