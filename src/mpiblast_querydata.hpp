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
#ifndef __MPIBLAST_QUERYDATA_H__
#define __MPIBLAST_QUERYDATA_H__

extern "C" {
#include "ncbi.h"
}
#include "blastjob.hpp"

class QueryData {
public:
	QueryData( BlastJob* bj, ValNodePtr bs_cache ) :
		blast_job( bj ), bioseq_cache( bs_cache ),
		searched_frags( 0 ) {};
	QueryData( const QueryData& qd ) :
	blast_job( qd.blast_job ),
	bioseq_cache( qd.bioseq_cache ),
	searched_frags( qd.searched_frags ),
	bioseq_lengths( qd.bioseq_lengths ),
	bioseq_starts( qd.bioseq_starts ),
	bioseq_strands( qd.bioseq_strands )
	{}
	QueryData& operator=( const QueryData& qd ){
		blast_job = qd.blast_job;
		bioseq_cache = qd.bioseq_cache;
		searched_frags = qd.searched_frags;
		bioseq_lengths = qd.bioseq_lengths;
		bioseq_starts = qd.bioseq_starts;
		bioseq_strands = qd.bioseq_strands;
		return *this;
	}
	~QueryData() {
		delete blast_job;
	}

	/** A BlastJob tracks which segments this query has been searched against */
	BlastJob* blast_job;
	int searched_frags;	/**< counts the number of fragments that have been searched */

	// writer only data
	/** 
	 * A cache for any Bioseqs that will be needed during result output.
	 * These get sent from the worker nodes that performed the search.
	 */
	ValNodePtr bioseq_cache;
	std::vector< int > bioseq_lengths;	/**< a list of the original bioseq lengths */
	std::vector< int > bioseq_starts;	/**< a list of the original bioseq starts */
	std::vector< int > bioseq_strands;	/**< a list of the original bioseq strands */
};

#endif
