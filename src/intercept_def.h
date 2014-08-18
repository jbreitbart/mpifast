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
#ifndef __INTERCEPT_DEF_H__
#define __INTERCEPT_DEF_H__

/* tags for intercepting */

#define INT_HDR_BEGIN 		0	// begin of header
#define INT_HDR_END 		1	// end of header
#define INT_DES_BEGIN		2   // begin of one line descriptons
#define INT_DES_TAG			3   // seperate of one line descriptions
#define INT_DES_END			4   // end of one line descriptions
#define INT_ALN_BEGIN		5   // begin of alignments output
#define INT_ALN_TAG 		6   // seperate of alignments output
#define INT_ALN_END 		7   // end of alignments output
#define INT_ALN_FIRST 		8   // first alignment in the group, printing '>'
#define INT_FTR_BEGIN 		9   // begin of footer
#define INT_FTR_END 		10  // end of footer
#define INT_DBI_BEGIN 		11  // begin of db info output
#define INT_DBI_END 		12  // end of dbinfo output
#define INT_STA_BEGIN		13  // begin of statistics output
#define INT_STA_END			14  // end of statistics output

extern int pio_hijacking;

int intercept_output_tabular(int tag);
int intercept_output_pairwise(int tag);

#endif /* __INTERCEPT_DEF_H__ */
