/* $Id: blast_dust.h,v 1.16 2007/01/23 15:25:22 madden Exp $
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
 * Author: Tom Madden
 *
 */

/** @file blast_dust.h
 * DUST filtering functions. (shouldn't this be merged with blast_filter?)
 * @todo FIXME: include reference?
 */

#ifndef ALGO_BLAST_CORE__BLAST_DUST__H
#define ALGO_BLAST_CORE__BLAST_DUST__H


#ifdef __cplusplus
extern "C" {
#endif


/** endpoints
 * linked list of low-complexity offsets.
 */
typedef struct DREGION {
        struct  DREGION*        next;
        Int4    from, to;
} DREGION;


Int4 DustSegs (Uint1* sequence, Int4 length, Int4 start,
                       DREGION* reg,
                       Int4 level, Int4 windowsize, Int4 linker);

SeqLocPtr BioseqDustEx (BioseqPtr bsp, Int4 level, Int4 windowsize, Int4 linker);
SeqLocPtr SeqLocDustEx (SeqLocPtr this_slp, Int4 level, Int4 window, Int4 linker);


#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_DUST__H */
