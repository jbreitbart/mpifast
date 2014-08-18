/*
 * $Id: aceread.c,v 1.16 2008/12/22 22:40:30 bollin Exp $
 *
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
 * Authors:  Colleen Bollin
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <util/creaders/alnread.h>
#include <aceread.h>


typedef enum {
    eTrue = -1,
    eFalse = 0
} EBool;


typedef enum {
    eJustRight = 0,
    eNone,
    eTooMany,
    eTooFew, 
    eUnexpected
} EFound;


extern void PrintACEFormatErrorXMLStart (char *id, char *has_errors)
{
    if (has_errors != NULL) {
        if (*has_errors == 0) {
            printf ("<aceread>\n");
            *has_errors = 1;
        }
    }
    printf ("<message severity=\"ERROR\" seq-id=\"%s\" code=\"bad_format\">", id == NULL ? "No ID" : id);
}


extern void PrintACEFormatErrorXMLEnd (void)
{
    printf ("</message>\n");
}


extern void PrintACEFormatErrorXML (char *msg, char *id, char *has_errors)
{
    if (has_errors != NULL) {
        if (*has_errors == 0) {
            printf ("<aceread>\n");
            *has_errors = 1;
        }
    }
    printf ("<message severity=\"ERROR\" seq-id=\"%s\" code=\"bad_format\">%s</message>\n", id == NULL ? "No ID" : id, msg);
}


static void s_ReportFound (EFound val, char *label, char *id, char *has_errors)
{
    switch (val) {
        case eNone:
            PrintACEFormatErrorXMLStart (id, has_errors);
            printf ("Found no %s", label);
            PrintACEFormatErrorXMLEnd ();
            break;
        case eTooMany:
            PrintACEFormatErrorXMLStart (id, has_errors);
            printf ("Too many %s", label);
            PrintACEFormatErrorXMLEnd ();
            break;
        case eTooFew:
            PrintACEFormatErrorXMLStart (id, has_errors);
            printf ("Too few %s", label);
            PrintACEFormatErrorXMLEnd ();
            break;
        case eUnexpected:
            PrintACEFormatErrorXMLStart (id, has_errors);
            printf ("Unexpected character while reading %s", label);
            PrintACEFormatErrorXMLEnd ();
            break;
        case eJustRight:
            break;
        default:
            PrintACEFormatErrorXML ("Unknown error", id, has_errors);
            break;
    }
}


extern TGapInfoPtr GapInfoNew (void)
{
    TGapInfoPtr g;

    g = (TGapInfoPtr) malloc (sizeof (SGapInfo));
    if (g != NULL) {
        g->num_gaps = 0;
        g->gap_offsets = NULL;
    }
    return g;
}


extern void GapInfoFree (TGapInfoPtr g)
{
    if (g != NULL) {
        free (g->gap_offsets);
        free (g);
    }
}


static int s_IsGapChar (char ch, char *gap_chars)
{
    if (ch == 0 || gap_chars == NULL) {
        return 0;
    }
    while (*gap_chars != 0 && *gap_chars != ch) {
        gap_chars ++;
    }
    if (*gap_chars == ch) {
        return 1;
    } else {
        return 0;
    }
}


/* The Trace Archive Gap String is a list of the number of nucleotides to skip before adding the next gap */
extern TGapInfoPtr GapInfoFromSequenceString (char *seq_str, char *gap_chars)
{ 
    char * cp;
    int    num_gaps = 0, pos, gap_num = 0;
    TGapInfoPtr g = NULL;

    if (seq_str == NULL) return NULL;

    /* first determine number of gaps */
    cp = seq_str;
    while (*cp != 0) {
        if (s_IsGapChar(*cp, gap_chars)) {
            num_gaps++;
        }
        cp++;
    }

    g = GapInfoNew ();
    if (num_gaps > 0) {
        g->num_gaps = num_gaps;
        g->gap_offsets = malloc (g->num_gaps * sizeof (int));
        cp = seq_str;
        pos = 0;
        while (*cp != 0) {
            if (s_IsGapChar(*cp, gap_chars)) {
                g->gap_offsets[gap_num] = pos;
                gap_num++;
                pos = 0;
            } else {
                pos++;
            }
            cp++;
        }
    }
    return g;
}

extern void RemoveGapCharsFromSequenceString (char *seq_str, char *gap_chars)
{
    char *cp_src, *cp_dst;

    if (seq_str == NULL || gap_chars == NULL) {
      return;
    }

    cp_src = seq_str;
    cp_dst = seq_str;
    while (*cp_src != 0) {
        if (!s_IsGapChar(*cp_src, gap_chars)) {
            *cp_dst = *cp_src;
            cp_dst++;
        }
        cp_src++;
    }
}


/* calculate sequence position from tiling position (both values zero-based) given gap_info */
extern int SeqPosFromTilingPos (int tiling_pos, TGapInfoPtr gap_info)
{
    int pos = 0, seq_pos = 0, gap_num = 0;

    if (tiling_pos < 0 || gap_info == NULL || gap_info->num_gaps == 0) {
        return tiling_pos;
    }
    
    while (gap_num < gap_info->num_gaps && pos + gap_info->gap_offsets[gap_num] <= tiling_pos) {
        seq_pos += gap_info->gap_offsets[gap_num];
        pos += gap_info->gap_offsets[gap_num] + 1;
        gap_num++;
    }
    seq_pos += tiling_pos - pos;
    return seq_pos;
}


/* calculate sequence position from tiling position (both values zero-based) given gap_info */
extern int TilingPosFromSeqPos (int seq_pos, TGapInfoPtr gap_info)
{
    int pos = 0, tiling_pos = 0, gap_num = 0;

    if (seq_pos < 0 || gap_info == NULL || gap_info->num_gaps == 0) {
        return seq_pos;
    }

    while (gap_num < gap_info->num_gaps && pos + gap_info->gap_offsets[gap_num] <= seq_pos) {
        pos += gap_info->gap_offsets[gap_num];
        tiling_pos += gap_info->gap_offsets[gap_num] + 1;
        gap_num++;
    }
    tiling_pos += seq_pos - pos;
    return tiling_pos;
}


/* adjust gap info when sequence is trimmed */
static void AdjustGapInfoFor5Trim (TGapInfoPtr gap_info, int trim)
{
    int pos = 0;
    int num_gaps = 0;
    int i;

    if (gap_info == NULL || gap_info->num_gaps < 1 || trim < 1) {
        return;
    }

    while (num_gaps < gap_info->num_gaps && pos + gap_info->gap_offsets[num_gaps] < trim) {
        pos += gap_info->gap_offsets[num_gaps];
        num_gaps++;
    }
    if (num_gaps < gap_info->num_gaps) {
        gap_info->gap_offsets[num_gaps] -= trim - pos;
        for (i = num_gaps; i < gap_info->num_gaps; i++) {
            gap_info->gap_offsets[i - num_gaps] = gap_info->gap_offsets[i];
        }
        gap_info->num_gaps -= num_gaps;
    } else {
        free (gap_info->gap_offsets);
        gap_info->gap_offsets = NULL;
        gap_info->num_gaps = 0;
    }

}


static void AdjustGapInfoFor3Trim (TGapInfoPtr gap_info, int new_len)
{
    int pos = 0;
    int num_gaps = 0;

    if (gap_info == NULL || gap_info->num_gaps < 1) {
        return;
    }

    while (num_gaps < gap_info->num_gaps && pos + gap_info->gap_offsets[num_gaps] < new_len) {
        pos += gap_info->gap_offsets[num_gaps];
        num_gaps++;
    }
    if (num_gaps < gap_info->num_gaps) {
        gap_info->num_gaps = num_gaps;
    }
}

/* TODO: NEED TO write function for truncating on right, test function for truncating on left */

extern TContigReadPtr ContigReadNew (void)
{
    TContigReadPtr r;

    r = (TContigReadPtr) malloc (sizeof (SContigRead));
    if (r == NULL) {
        return NULL;
    }
    r->read_id = NULL;
    r->ti = 0;
    r->srr = NULL;
    r->read_seq = NULL;
    r->is_complement = 0;
    r->cons_start = 0;
    r->cons_stop = 0;
    r->gaps = NULL;
    r->local = 1;
    r->valid = 0;
    r->qual_scores = NULL;
    r->num_qual_scores = 0;
    r->tag = NULL;
    return r;
}


extern void ContigReadFree (TContigReadPtr r)
{
    if (r != NULL) {
        if (r->read_id != NULL) {
            free (r->read_id);
        }
        if (r->srr != NULL) {
            free (r->srr);
        }
        if (r->read_seq != NULL) {
            free (r->read_seq);
        }
        if (r->gaps != NULL) {
            GapInfoFree (r->gaps);
        }
        if (r->qual_scores != NULL) {
            free (r->qual_scores);
        }
        if (r->tag != NULL) {
            free (r->tag);
        }
        free (r);
    }
}


extern TBaseSegPtr BaseSegNew (void)
{
    TBaseSegPtr b;

    b = (TBaseSegPtr) malloc (sizeof (SBaseSeg));
    if (b == NULL) {
        return NULL;
    }
    b->read_id = NULL;
    b->cons_start = 0;
    b->cons_stop = 0;
    return b;
}


extern void BaseSegFree (TBaseSegPtr b)
{
    if (b != NULL) {
        if (b->read_id != NULL) {
            free (b->read_id);
        }
        free (b);
    }
}


/* reads a correctly formatted line and creates a base seg.
 */
static TBaseSegPtr s_ReadBaseSeg (char *line)
{
    TBaseSegPtr base_seg = NULL;
    char *cp;
    int   start, stop, len;

    if (line == NULL || *line != 'B' || *(line + 1) != 'S') {
        return NULL;
    }


    cp = line + 2;
    while (isspace (*cp)) {
        cp++;
    }
    if (!isdigit (*cp)) {
        return NULL;
    }
    start = atoi (cp);
    while (isdigit (*cp)) {
        cp++;
    }
    while (isspace (*cp)) {
        cp++;
    }
    if (!isdigit (*cp)) {
        return NULL;
    }
    stop = atoi (cp);
    while (isdigit (*cp)) {
        cp++;
    }
    while (isspace (*cp)) {
        cp++;
    }
    if (*cp == 0) {
        return NULL;
    }

    len = strlen (cp);

    base_seg = BaseSegNew ();
    base_seg->cons_start = start;
    base_seg->cons_stop = stop;
    base_seg->read_id = malloc (sizeof (char) * len + 1);
    strcpy (base_seg->read_id, cp);

    return base_seg;
}


extern TConsensusReadAlnPtr ConsensusReadAlnNew (int numseg)
{
    TConsensusReadAlnPtr a;
    int i;

    a = (TConsensusReadAlnPtr) malloc (sizeof (SConsensusReadAln));
    a->is_complement = 0;
    if (numseg < 1) {
        a->lens = NULL;
        a->cons_starts = NULL;
        a->read_starts = NULL;
        a->numseg = 0;
    } else {
        a->lens = (int *) malloc (sizeof (int) * numseg);
        a->cons_starts = (int *) malloc (sizeof (int) * numseg);
        a->read_starts = (int *) malloc (sizeof (int) * numseg);
        for (i = 0; i < numseg; i++) {
            a->lens[i] = 0;
            a->cons_starts[i] = 0;
            a->read_starts[0] = 0;
        }
        a->numseg = numseg;
    }
    return a;
}


extern TConsensusReadAlnPtr ConsensusReadAlnFree (TConsensusReadAlnPtr a)
{
    if (a != NULL) {
        if (a->lens != NULL) {
            free (a->lens);
            a->lens = NULL;
        }
        if (a->cons_starts != NULL) {
            free (a->cons_starts);
            a->cons_starts = NULL;
        }
        if (a->read_starts != NULL) {
            free (a->read_starts);
            a->read_starts = NULL;
        }
        free (a);
        a = NULL;
    }
    return a;
}


extern TConsensusReadAlnPtr GetConsensusReadAln (char *consensus_seq, TContigReadPtr read)
{
    TConsensusReadAlnPtr aln = NULL;
    char *c;
    char *c_start;
    char *r;
    char *r_start;
    int numseg = 0, aln_len, pos, seg, con_offset = 0, read_offset = 0;
    char con_gap_open = 0, read_gap_open = 0, gap_change;

    if (consensus_seq == NULL || read == NULL) {
        return NULL;
    }

    if (read->cons_start > 0) {
        c_start = consensus_seq + read->cons_start;
        r_start = read->read_seq + read->read_assem_start - 1;
    } else {
        c_start = consensus_seq;
        r_start = read->read_seq + read->read_assem_start - 1;
    }

    aln_len = read->cons_stop - read->cons_start + 1;
    while (*c_start == '*' && *r_start == '*') {
        c_start++;
        r_start++;
        aln_len--;
    }

    /* first, count number of segments needed */
    c = c_start;
    r = r_start;
    if (*c != '*' && *r != '*') {
      numseg++;
    }
    pos = 0;
    while (*c != 0 && *r != 0 && pos < aln_len) {
        if (*c == '*' && *r == '*') {
            /* both in gap - ignore */
        } else {
            gap_change = 0;
            if (*c == '*') {
                if (!con_gap_open) {
                    gap_change = 1;
                    con_gap_open = 1;
                }
            } else {
                if (con_gap_open) {
                    gap_change = 1;
                    con_gap_open = 0;
                }
            }
            if (*r == '*') {
                if (!read_gap_open) {
                    gap_change = 1;
                    read_gap_open = 1;
                }
            } else {
                if (read_gap_open) {
                    gap_change = 1;
                    read_gap_open = 0;
                }
            }
            if (gap_change) {
                numseg++;
            }
        }
        c++;
        r++;
        pos++;
    }

    /* create alignment */
    aln = ConsensusReadAlnNew (numseg);
    pos = 0;
    seg = 0;


    c = consensus_seq;        
    while (c < c_start) {
        if (*c != '*') {
            con_offset ++;
        }
        c++;
    }

    r = read->read_seq;
    while (r < r_start) {
        if (*r != '*') {
            read_offset ++;
        }
        r++;
    }
    
    
    if (*c_start == '*') {
        aln->cons_starts[0] = -1;
        con_gap_open = 1;
    } else {
        aln->cons_starts[0] = con_offset;
        con_gap_open = 0;
    }

    if (*r_start == '*') {
        aln->read_starts[0] = -1;
        read_gap_open = 1;
    } else {
        aln->read_starts[0] = read_offset;
        read_gap_open = 0;
    }

    c = c_start + 1;
    r = r_start + 1;
    aln->lens[0] = 1;
    pos = 1;
    
    while (*c != 0 && *r != 0 && pos < aln_len) {
        if (*c == '*' && *r == '*') {
            /* both in gap - ignore */
        } else {
            gap_change = 0;
            if (*c == '*') {
                if (!con_gap_open) {
                    gap_change = 1;
                    con_gap_open = 1;
                }
            } else {
                if (con_gap_open) {
                    gap_change = 1;
                    con_gap_open = 0;
                }
            }
            if (*r == '*') {
                if (!read_gap_open) {
                    gap_change = 1;
                    read_gap_open = 1;
                }
            } else {
                if (read_gap_open) {
                    gap_change = 1;
                    read_gap_open = 0;
                }
            }
            if (gap_change) {
                seg++;
                if (con_gap_open) {
                    aln->cons_starts[seg] = -1;
                } else if (aln->cons_starts[seg - 1] > -1) {
                    aln->cons_starts[seg] = aln->cons_starts[seg - 1] + aln->lens[seg - 1];
                } else if (seg > 1 && aln->cons_starts[seg - 2] > -1) {
                    aln->cons_starts[seg] = aln->cons_starts[seg - 2] + aln->lens[seg - 2];
                } else {
                    aln->cons_starts[seg] = con_offset;
                }
                if (read_gap_open) {
                    aln->read_starts[seg] = -1;
                } else if (aln->read_starts[seg - 1] > -1) {
                    aln->read_starts[seg] = aln->read_starts[seg - 1] + aln->lens[seg - 1];
                } else if (seg > 1 && aln->read_starts[seg - 2] > -1) {
                    aln->read_starts[seg] = aln->read_starts[seg - 2] + aln->lens[seg - 2];
                } else {
                    aln->read_starts[seg] = read_offset;
                }
            }
            aln->lens[seg]++;
        }
        c++;
        r++;
        pos++;
    }

    /* todo - adjust starts for complement */
    if (read->is_complement) {
      for (seg = 0; seg < aln->numseg; seg++) {
        if (aln->read_starts[seg] > -1) {
          aln->read_starts[seg] = read->read_len - aln->read_starts[seg] - aln->lens[seg];
        }
      }
      aln->is_complement = 1;
    }

    return aln;
}


extern TContigPtr ContigNew (void)
{
    TContigPtr c;

    c = (TContigPtr) malloc (sizeof (SContig));
    if (c == NULL) {
        return NULL;
    }
    c->consensus_id = NULL;
    c->consensus_seq = NULL;
    c->consensus_assem_len = 0;
    c->consensus_seq_len = 0;
    c->is_complement = 0;
    c->num_qual_scores = 0;
    c->qual_scores = NULL;
    c->num_reads = 0;
    c->reads = NULL;
    c->gaps = NULL;
    c->num_reads = 0;
    c->reads = NULL;
    c->num_base_segs = 0;
    c->base_segs = NULL;
    c->tag = NULL;

    return c;
}


extern void ContigFree (TContigPtr c)
{
    int i;

    if (c != NULL) {
        if (c->consensus_id != NULL) free (c->consensus_id);
        if (c->consensus_seq != NULL) free (c->consensus_seq);
        if (c->qual_scores != NULL) free (c->qual_scores);
            
        if (c->reads != NULL) {
            for (i = 0; i < c->num_reads; i++) {
                if (c->reads[i] != NULL) {
                    ContigReadFree (c->reads[i]);
                }
            }
            free (c->reads);
        }
        if (c->base_segs != NULL) {
            for (i = 0; i < c->num_base_segs; i++) {
                if (c->base_segs[i] != NULL) {
                    BaseSegFree (c->base_segs[i]);
                }
            }
            free (c->base_segs);
        }
        if (c->tag != NULL) {
            free (c->tag);
        }
        free (c);
    }
}


extern TACEFilePtr ACEFileNew ()
{
    TACEFilePtr afp;

    afp = (TACEFilePtr) malloc (sizeof (SACEFile));
    if (afp == NULL) {
        return NULL;
    }
    afp->num_contigs = 0;
    afp->contigs = NULL;

    return afp;
}


extern void ACEFileFree (TACEFilePtr afp)
{
    int i;

    if (afp != NULL) {
        for (i = 0; i < afp->num_contigs; i++) {
            ContigFree (afp->contigs[i]);
        }
        free (afp->contigs);
        free (afp);
    }      
}


static char s_IsSeqChar (char ch)
{
    if (ch == '*' || isalpha (ch)) {
        return 1;
    } else {
        return 0;
    }
}


static char s_IsEOF (char *linestring)
{
    if (linestring == NULL || linestring [0] == EOF) {
        return 1;
    } else {
        return 0;
    }
}


static char *
s_ReadSequenceFromFile
(int                  len, 
 FReadLineFunction    readfunc,
 void *               userdata,
 char *               id,
 char *               has_errors)
{
    char *seq;
    char *linestring;
    char *cp;
    int  pos = 0;

    /* copy in sequence data */
    seq = malloc (len + 1);
    linestring = readfunc (userdata);
    while (!s_IsEOF (linestring) && s_IsSeqChar (linestring [0])) {
        /* append to consensus */
        cp = linestring;
        while (s_IsSeqChar (*cp) && pos < len) {
            if (isalpha (*cp)) {
                seq [pos] = toupper (*cp);
            } else {
                seq [pos] = *cp;
            }
            pos++;
            cp++;
        }
        if (s_IsSeqChar (*cp)) {
            PrintACEFormatErrorXML ("Too many sequence characters!", id, has_errors);
            free (seq);
            return NULL;
        }
        free (linestring);
        linestring = readfunc (userdata);
    }
    free (linestring);
    if (pos < len) {
        PrintACEFormatErrorXML ("Too few sequence characters!", id, has_errors);
        free (seq);
        seq = NULL;
    } else {
        seq[pos] = 0;
    }
    return seq;
}


static char s_LineIsEmptyButNotEof (char *linestring)
{
    char *cp;
    if (s_IsEOF (linestring)) {
        return 0;
    } 

    cp = linestring;
    while (*cp != 0 && isspace (*cp)) {
        cp++;
    }
    if (*cp == 0) {
        return 1;
    } else {
        return 0;
    }
}


static void s_SkipQualScores
(FReadLineFunction    readfunc,
 void *               userdata)
{
    char * linestring;
    char * cp;
    if (readfunc == NULL) return;

    linestring = readfunc (userdata);
    while (s_LineIsEmptyButNotEof (linestring)) {
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (linestring == NULL  ||  linestring [0] == EOF || strcmp (linestring, "BQ") != 0) {
        return;
    }
    linestring = readfunc (userdata);
    while (!s_IsEOF (linestring)
           && isdigit (*(cp = linestring + strspn (linestring, " \t")))) {
        free (linestring);   
        linestring = readfunc (userdata);
    }
    free (linestring);
}


static EFound s_ReadQualScores
(TContigPtr contig,
 FReadLineFunction    readfunc,
 void *               userdata)
{
    char * linestring;
    char * cp;
    int    pos;

    if (contig == NULL || readfunc == NULL || contig->consensus_assem_len == 0) {
        return eNone;
    }

    linestring = readfunc (userdata);
    while (s_LineIsEmptyButNotEof (linestring)) {
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (linestring == NULL  ||  linestring [0] == EOF || strcmp (linestring, "BQ") != 0) {
        return eNone;
    }

    /* read quality scores */
    contig->num_qual_scores = contig->consensus_assem_len;
    /* no score for * in consensus seq */
    for (pos = 0; pos < contig->consensus_assem_len; pos++) {
      if (contig->consensus_seq[pos] == '*') {
        contig->num_qual_scores --;
      }
    }
    contig->qual_scores = malloc (sizeof (int) * contig->num_qual_scores);
    pos = 0;
    linestring = readfunc (userdata);
    while (!s_IsEOF (linestring)
           && isdigit (*(cp = linestring + strspn (linestring, " \t")))) {
        while (isdigit (*cp) && pos < contig->num_qual_scores) {
            contig->qual_scores [pos] = atoi (cp);
            pos++;
            while (isdigit (*cp)) {
                cp++;
            }
            while (isspace (*cp)) {
                cp++;
            }
        }
        if (isdigit (*cp)) {
            return eTooMany;
        }
        free (linestring);   
        linestring = readfunc (userdata);
    }
    if (pos < contig->num_qual_scores) {
        return eTooFew;
    } else {
        return eJustRight;
    }
}


static EFound s_ReadAFLines
(TContigPtr contig,
 FReadLineFunction    readfunc,
 void *               userdata,
 char **              next_line)
{
    char * linestring;
    char * cp;
    int    read_num, len;
    EFound rval = eJustRight;

    if (contig == NULL || readfunc == NULL || contig->num_reads == 0) return eNone;

    /* get AF lines */
    contig->reads = malloc (contig->num_reads * sizeof (TContigReadPtr));
    linestring = readfunc (userdata);
    while (s_LineIsEmptyButNotEof (linestring)) {
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (linestring == NULL  ||  linestring [0] == EOF || strncmp (linestring, "AF", 2) != 0) {
        *next_line = linestring;
        return eNone;
    }
    
    read_num = 0;
    while (!s_IsEOF(linestring) && read_num < contig->num_reads
           && linestring [0] == 'A' && linestring [1] == 'F' && isspace (linestring [2])) {
        contig->reads[read_num] = ContigReadNew ();
        len = strlen (linestring + 3);
        contig->reads[read_num]->read_id = malloc (len + 1);
        strcpy (contig->reads[read_num]->read_id, linestring + 3);
        cp = contig->reads[read_num]->read_id;
        while (*cp != 0 && !isspace (*cp)) {
            cp++;
        }
        if (isspace (*cp)) {
            *cp = 0;
            cp++;
        }
        if (*cp == 'C') {
            contig->reads[read_num]->is_complement = 1;
        } else if (*cp != 'U') {
            *next_line = linestring;
            return eUnexpected;
        }
        cp++;
        if (isspace (*cp)) {
            cp++;
        }
        contig->reads[read_num]->cons_start = atoi (cp) - 1;
        read_num++;
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (read_num < contig->num_reads) {
        rval = eTooFew;
    } else if (!s_IsEOF(linestring) && strncmp (linestring, "AF ", 3) == 0) {
        rval = eTooMany;
    } else {
        rval = eJustRight;
    }
    *next_line = linestring;
    return rval;
}


static EFound s_ReadBaseSegs
(TContigPtr           contig,
 int                  num_base_segs,
 char *               firstline,
 FReadLineFunction    readfunc,
 void *               userdata)
{
    char * linestring;

    if (contig == NULL || readfunc == NULL || num_base_segs == 0) return eNone;

    contig->base_segs = malloc (sizeof (TBaseSegPtr) * num_base_segs);
    contig->num_base_segs = 0;

    /* get BS lines */
    linestring = firstline;
    while (s_LineIsEmptyButNotEof (linestring)) {
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (linestring == NULL  ||  linestring [0] == EOF || strncmp (linestring, "BS", 2) != 0) {
        return eNone;
    }
    
    while (linestring != NULL  &&  linestring [0] != EOF && contig->num_base_segs < num_base_segs
           && linestring [0] == 'B' && linestring [1] == 'S' && isspace (linestring [2])) {
        contig->base_segs[contig->num_base_segs++] = s_ReadBaseSeg (linestring);
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (contig->num_base_segs < num_base_segs) {
        return eTooFew;
    } else if (linestring != NULL && linestring [0] != EOF && ! s_LineIsEmptyButNotEof (linestring)) {
        return eTooMany;
    } else {
        return eJustRight;
    }
}


static char s_IsEquivN (char ch)
{
    if (ch == 'N' || ch == 'X') {
        return 1;
    } else {
        return 0;
    }
}


/* Terminal Ns will always be trimmed in the GenBank records */
static void s_AdjustContigReadForTerminalNs (TContigReadPtr read)
{
    char * cp_src;
    char * cp_dst;
    int    len = 0;

    if (read == NULL || read->read_seq == NULL) return;
    cp_src = read->read_seq;
    while (s_IsEquivN(*cp_src)) {
        len++;
        cp_src++;
    }
    if (len > 0) {
        read->cons_start += len;
        cp_dst = read->read_seq;
        while (*cp_src != 0) {
            *cp_dst = *cp_src;
            cp_dst++;
            cp_src++;
        }
        *cp_dst = 0;
    }
    len = strlen (read->read_seq);
    cp_src = read->read_seq + len - 1;
    while (cp_src >= read->read_seq && s_IsEquivN(*cp_src)) {
        *cp_src = 0;
        cp_src--;
    }
} 


/* Terminal Ns will always be trimmed by the GenBank record */
static void s_AdjustContigForTerminalNs (TContigPtr contig)
{
    char * cp_src;
    char * cp_dst;
    int    len = 0, i;

    if (contig == NULL || contig->consensus_seq == NULL) return;
    cp_src = contig->consensus_seq;
    while (s_IsEquivN(*cp_src)) {
        len++;
        cp_src++;
    }
    if (len > 0) {
        /* adjust quality scores */
        if (contig->qual_scores != NULL) {
            contig->num_qual_scores -= len;
            for (i = 0; i < contig->num_qual_scores; i++) {
                contig->qual_scores[i] = contig->qual_scores [i + len];
            }
        }
        /* adjust reads */
        if (contig->reads != NULL) {
            for (i = 0; i < contig->num_reads; i++) {
                if (contig->reads[i] != NULL) {
                    contig->reads[i]->cons_start -= len;
                }
            }
        }
        /* adjust consensus sequence */
        cp_dst = contig->consensus_seq;
        while (*cp_src != 0) {
            *cp_dst = *cp_src;
            cp_dst++;
            cp_src++;
        }
        *cp_dst = 0;
        contig->consensus_assem_len -= len;
    }
    /* trim 3' Ns */
    len = 0;
    cp_src = contig->consensus_seq + contig->consensus_assem_len - 1;
    while (cp_src >= contig->consensus_seq && s_IsEquivN(*cp_src)) {
        *cp_src = 0;
        cp_src--;
        contig->consensus_assem_len--;
        len++;
    }
    /* truncate quality scores if 3' Ns trimmed */
    if (len > 0 && contig->qual_scores != NULL) {
        contig->num_qual_scores -= len;
    }
} 


/* Clips the sequence read in according to the QA clipping.
 * The real coordinates will be recovered when an alignment is generated between
 * the sequence in the structure and the sequence downloaded from the Trace Archive.
 */
static char ApplyQALineToRead (TContigReadPtr read, char *linestring, char *id, char *has_errors)
{
    char *cp;
    int  values[4];
    int  i = 0;

    if (read == NULL || linestring == NULL) {
        PrintACEFormatErrorXML ("File end where QA line should be", id, has_errors);
        return 0;
    }
  
    cp = linestring;
    if (*cp != 'Q') {
        PrintACEFormatErrorXMLStart (id, has_errors);
        printf ("Expected QA line, found %s", linestring);
        PrintACEFormatErrorXMLEnd ();
        return 0;
    }
    cp++;
    if (*cp != 'A') {
        PrintACEFormatErrorXMLStart (id, has_errors);
        printf ("Expected QA line, found %s", linestring);
        PrintACEFormatErrorXMLEnd ();
        return 0;
    }
    cp++;
    while (*cp != 0 && i < 4) {
        while (isspace (*cp)) {
            cp++;
        }
        if (*cp != '-' && !isdigit (*cp)) {
          PrintACEFormatErrorXML ("Found non-number on QA line", id, has_errors);
          return 0;
        }
        values[i] = atoi (cp);
        i++;
        while (*cp == '-' || isdigit (*cp)) {
            cp++;
        }
    }
    if (*cp != 0 || i < 4) {
        PrintACEFormatErrorXML ("Fewer than four numbers on line", id, has_errors);
        return 0;
    }
    if (values[0] > 0 || values[2] > 0) {
        if (values[0] > values[2]) {
            read->read_assem_start = values[0];
        } else {
            read->read_assem_start = values[2];
        }
    }

    if (values[1] > 0 && values[3] > 0) {
        if (values[1] < values[3]) {
            read->read_assem_stop = values[1];
        } else {
            read->read_assem_stop = values[3];
        }
    } else if (values[1] > 0) {
        read->read_assem_stop = values[1];
    } else if (values[3] > 0) {
        read->read_assem_stop = values[3];
    }

    /* adjust first gap position for start */
    if (read->read_assem_start > 1 && read->gaps != NULL && read->gaps->num_gaps > 0 && read->gaps->gap_offsets != NULL) {
        read->gaps->gap_offsets[0] -= read->read_assem_start - 1;
    }
        
    return 1;
}
    

/* calculate gap info for consensus sequence */
/* calculate cons_stop positions and tiling positions for each read */
static void s_CalculateContigOffsets (TContigPtr contig)
{
    int i;

    if (contig == NULL) return;

    for (i = 0; i < contig->num_reads; i++) {
        contig->reads[i]->tiling_start = contig->reads[i]->read_assem_start + contig->reads[i]->cons_start;
        contig->reads[i]->tiling_stop = contig->reads[i]->read_assem_stop + contig->reads[i]->cons_start;
        contig->reads[i]->cons_stop = SeqPosFromTilingPos (contig->reads[i]->tiling_stop - 1, contig->gaps) + 1;
        contig->reads[i]->read_start = SeqPosFromTilingPos(contig->reads[i]->read_assem_start - 1, contig->reads[i]->gaps) + 1;
        contig->reads[i]->read_stop = SeqPosFromTilingPos(contig->reads[i]->read_assem_stop - 1, contig->reads[i]->gaps) + 1;
    }

}


static int s_GetUngappedSeqLen (char *str, char *gap_chars)
{
    int len = 0;

    if (str == NULL) return 0;
    while (*str != 0) {
        if (!s_IsGapChar (*str, gap_chars)) {
            len++;
        }
        str++;
    }
    return len;
}


static char * s_AddToTagComment (char *orig, char *extra)
{
    char * tag = NULL;
    int    tag_len;

    if (orig == NULL) {
        tag = extra;
    } else {
        tag_len = strlen (orig) + strlen (extra) + 1;
        tag = malloc (sizeof (char) * (tag_len + 1));
        strcpy (tag, orig);
        strcat (tag, "\n");
        strcat (tag, extra);
        free (orig);
        free (extra);
    }
    return tag;
}


static char * s_ReadTagComment
(FReadLineFunction    readfunc,
 void *               userdata)
{
    char *linestring;
    char *tag = NULL;
    char *cp = NULL;
    char *tmp;
    int   tag_len = 0, end_len;

    linestring = readfunc (userdata);
    while (linestring != NULL  &&  linestring [0] != EOF && (cp = strchr (linestring, '}')) == NULL) {
        if (tag == NULL) {
            tag_len = strlen (linestring);
            tag = malloc (sizeof (char) * (tag_len + 1));
            strcpy (tag, linestring);
        } else {
            tag_len = tag_len + strlen (linestring) + 1;
            tmp = malloc (sizeof (char) * (tag_len + 1));
            strcpy (tmp, tag);
            strcat (tmp, "\n");
            strcat (tmp, linestring);
            free (tag);
            tag = tmp;
        }
        free (linestring);
        linestring = readfunc (userdata);
    }
    if (cp != NULL && cp > linestring) {
        end_len = cp - linestring;
        tag_len = tag_len + end_len + 1;
        tmp = malloc (sizeof (char) * (tag_len + 1));
        strcpy (tmp, tag);
        strcat (tmp, "\n");
        strncat (tmp, linestring, end_len);
        tmp[tag_len] = 0;
        free (tag);
        tag = tmp;
    }
    if (linestring != NULL) {
        free (linestring);
    }

    return tag;
}


/* Reads the portion of and ACE file for a single contig, including the reads */
static TContigPtr s_ReadContig
(char **              initline,
 FReadLineFunction    readfunc,
 void *               userdata,
 char                 make_qual_scores,
 char *               has_errors)
{
    char      *linestring;
    char      *firstline;
    char      *cp;
    int        len = 0, read_num = 0, num_base_segs = 0;
    EFound     val;
    char       found_comp_char = 0;
    TContigPtr contig = NULL;

    if (initline == NULL) return NULL;
    firstline = *initline;
    if (firstline == NULL || readfunc == NULL) return NULL;

    if (firstline [0] != 'C' || firstline [1] != 'O' || ! isspace (firstline [2])) {
        return NULL;
    }

    contig = ContigNew ();
    len = strlen (firstline + 3);
    contig->consensus_id = malloc (len + 1);
    strcpy (contig->consensus_id, firstline + 3);
 
    cp = contig->consensus_id;
    while (*cp != 0 && !isspace (*cp)) {
        cp++;
    }
    if (isspace (*cp)) {
        *cp = 0;
        cp++;
        contig->consensus_assem_len = atoi (cp);
        while (isdigit (*cp)) {
            cp++;
        }
        if (isspace (*cp)) {
            cp++;
            contig->num_reads = atoi (cp);
            while (isdigit (*cp)) {
                cp++;
            } 
            if (isspace (*cp)) {
                cp++;
                num_base_segs = atoi (cp);
                while (isdigit (*cp)) {
                    cp++;
                }
                if (isspace (*cp)) {
                    cp++;
                    found_comp_char = 1;
                    if (*cp == 'C') {
                        contig->is_complement = 1;
                    } else {
                        contig->is_complement = 0;
                    }
                } 
            }
        }
    }
    if (contig->consensus_assem_len == 0 || contig->num_reads == 0 || !found_comp_char) {
        PrintACEFormatErrorXML ("Error in consensus line", contig->consensus_id, has_errors);
        ContigFree (contig);
        return NULL;
    }
        
    /* now copy in sequence data */
    contig->consensus_seq = s_ReadSequenceFromFile (contig->consensus_assem_len, readfunc, userdata, contig->consensus_id, has_errors);
    if (contig->consensus_seq == NULL) {
        ContigFree (contig);
        return NULL;
    }

    /* record actual length of consensus seq */
    contig->consensus_seq_len = s_GetUngappedSeqLen (contig->consensus_seq, "*");

    /* calculate gap info */
    contig->gaps = GapInfoFromSequenceString (contig->consensus_seq, "*");
    
    /* read quality scores */
    if (make_qual_scores) {
        val = s_ReadQualScores (contig, readfunc, userdata);
        if (val != eNone && val != eJustRight) {
            s_ReportFound (val, "quality scores", contig->consensus_id, has_errors);
            ContigFree (contig);
            return NULL;
        }
    } else {
        s_SkipQualScores (readfunc, userdata);
    }

    /* collect reads */
    val = s_ReadAFLines (contig, readfunc, userdata, &linestring);
    if (val != eJustRight) {
        s_ReportFound (val, "AF lines", contig->consensus_id, has_errors);
        ContigFree (contig);
        if (linestring != NULL) free (linestring);
        return NULL;
    }
 
    if (num_base_segs > 0) {
        val = s_ReadBaseSegs (contig, num_base_segs, linestring, readfunc, userdata);
        if (val != eJustRight) {
            s_ReportFound (val, "base segments", contig->consensus_id, has_errors);
            ContigFree (contig);
            return NULL;
        }
    }

    
    read_num = 0;
    linestring = readfunc (userdata);
    while (linestring != NULL  &&  linestring [0] != EOF) {
        if (linestring [0] == 'R' && linestring[1] == 'D' && isspace (linestring [2])) {
            len = strlen (contig->reads[read_num]->read_id);
            if (strncmp (linestring + 3, contig->reads[read_num]->read_id, len) != 0
                || !isspace (linestring [3 + len])) {
                PrintACEFormatErrorXML ("Read IDs out of order!", contig->consensus_id, has_errors);
                ContigFree (contig);
                return NULL;
            } 
            len = atoi (linestring + 3 + len);
            contig->reads[read_num]->read_seq = s_ReadSequenceFromFile (len, readfunc, userdata, contig->reads[read_num]->read_id, has_errors);
            if (contig->reads[read_num]->read_seq == NULL) {
                ContigFree (contig);
                return NULL;
            }
            s_AdjustContigReadForTerminalNs (contig->reads[read_num]);
            contig->reads[read_num]->read_len = s_GetUngappedSeqLen (contig->reads[read_num]->read_seq, "*");
            contig->reads[read_num]->gaps = GapInfoFromSequenceString (contig->reads[read_num]->read_seq, "*");
            read_num++;
        } else if (linestring [0] == 'Q' && linestring[1] == 'A' && isspace (linestring[2])) {
            if (read_num < 1) {
                PrintACEFormatErrorXML ("Found QA line before RD!", contig->consensus_id, has_errors);
                ContigFree (contig);
                return NULL;
            } else if (!ApplyQALineToRead (contig->reads[read_num - 1], linestring, contig->reads[read_num - 1]->read_id, has_errors)) {
                PrintACEFormatErrorXML ("Error in QA line format!", contig->reads[read_num - 1]->read_id, has_errors);
                ContigFree (contig);
                return NULL;
            }
        } else if (linestring[0] == 'D' && linestring[1] == 'S' && isspace (linestring[2])) {
            /* skip DS lines */
        } else if (strncmp (linestring, "RT{", 3) == 0) {
            contig->reads[read_num - 1]->tag = s_AddToTagComment (contig->reads[read_num - 1]->tag, s_ReadTagComment (readfunc, userdata));
        } else if (strncmp (linestring, "WR{", 3) == 0) {
            contig->reads[read_num - 1]->tag = s_AddToTagComment (contig->reads[read_num - 1]->tag, s_ReadTagComment (readfunc, userdata));
        } else if (strncmp (linestring, "CT{", 3) == 0) {
            contig->tag = s_AddToTagComment (contig->tag, s_ReadTagComment (readfunc, userdata));
        } else if (strncmp (linestring, "WA{", 3) == 0) {
            contig->tag = s_AddToTagComment (contig->tag, s_ReadTagComment (readfunc, userdata));
        } else if (linestring[0] != 0) {
            /* found next line */
            *initline = linestring;
            s_AdjustContigForTerminalNs (contig);
            s_CalculateContigOffsets (contig);
            return contig;
        }
        free (linestring);
        linestring = readfunc (userdata);
    }
    *initline = NULL;
    s_AdjustContigForTerminalNs (contig);
    s_CalculateContigOffsets (contig);
    return contig;
}


/* Used to detect errors in ACE file formatting */
static char s_UnexpectedLineBetweenContigs (char *linestring)
{
    if (linestring == NULL) {
        return 0;
    } else if (linestring [0] == 'A' && linestring [1] == 'F') {
        return 1;
    } else if (linestring [0] == 'R' && linestring [1] == 'D') {
        return 1;
    } else {
        return 0;
    }
}


/* This is the main function for reading in an ACE file */
extern TACEFilePtr
ReadACEFile
(FReadLineFunction    readfunc,
 void *               userdata,
 char                 make_qual_scores,
 char *               has_errors)
{
    char *              linestring;
    TACEFilePtr         afp;
    char *              cp;
    int                 contig_num = 0, read_num = 0;
    int                 num_reads_expected = 0;

    if (readfunc == NULL) {
        return NULL;
    }

    afp = ACEFileNew ();
    if (afp == NULL) {
        return NULL;
    }
  
    linestring = readfunc (userdata);

    while (linestring != NULL  &&  linestring [0] != EOF) {
        if (linestring [0] == 'A' && linestring [1] == 'S' && isspace (linestring [2])) {
            if (num_reads_expected > 0) {
                PrintACEFormatErrorXML ("Two file header lines!", NULL, has_errors);
                ACEFileFree (afp);
                free (linestring);
                return NULL;
            }
            /* first line in file, number of contigs */
            cp = linestring + 3;
            afp->num_contigs = atoi (cp);
            afp->contigs = malloc (afp->num_contigs * sizeof (TContigPtr));
            if (afp->contigs == NULL) {
                PrintACEFormatErrorXML ("Memory allocation failed!", NULL, has_errors);
                free (linestring);
                ACEFileFree (afp);
                return NULL;
            }
            while (isdigit (*cp)) {
                cp++;
            }
            num_reads_expected = atoi (cp);
            free (linestring);
            linestring = readfunc (userdata);
        } else if (linestring [0] == 'C' && linestring [1] == 'O' && isspace (linestring [2])) {
            if (contig_num >= afp->num_contigs) {
                PrintACEFormatErrorXML ("Too many contigs!", NULL, has_errors);
                free (linestring);
                ACEFileFree (afp);
                return NULL;
            }
            afp->contigs[contig_num] = s_ReadContig (&linestring, readfunc, userdata, make_qual_scores, has_errors);
            if (afp->contigs[contig_num] == NULL) {
                PrintACEFormatErrorXMLStart (NULL, has_errors);
                printf ("Unable to read contig (%d)", contig_num);
                PrintACEFormatErrorXMLEnd ();
                ACEFileFree (afp);
                return NULL;
            }
            read_num += afp->contigs[contig_num]->num_reads;
            contig_num++;
        } else if (s_UnexpectedLineBetweenContigs (linestring)) {
            PrintACEFormatErrorXMLStart (NULL, has_errors);
            printf ("Unexpected line after contig %d:%s", read_num, linestring);
            PrintACEFormatErrorXMLEnd ();
            free (linestring);
            ACEFileFree (afp);
            return NULL;
        } else {
            free (linestring);
            linestring = readfunc (userdata);
        }
    }
    if (contig_num < afp->num_contigs) {
        PrintACEFormatErrorXML ("Not enough contigs!", NULL, has_errors);
        ACEFileFree (afp);
        afp = NULL;
    } else if (read_num < num_reads_expected) {
        PrintACEFormatErrorXML ("Not enough reads!", NULL, has_errors);
        ACEFileFree (afp);
        afp = NULL;
    }

    return afp;
}


/* This function writes out sequence characters, 60 per line. */
static void s_WriteSeq (FILE *fp, char *seq)
{
    int    i;
    char * cp;

    if (fp == NULL || seq == NULL) return;
    cp = seq;
    while (*cp != 0) {
        for (i = 0; i < 60 && *cp != 0; i++, cp++) {
            fprintf (fp, "%c", *cp);
        }
        fprintf (fp, "\n");
    }
}


/* This function writes out quality scores in the ACE file format. */
static void s_WriteQualScores (FILE *fp, TContigPtr contig)
{
    int q_pos, line_pos;

    if (fp == NULL || contig == NULL || contig->num_qual_scores == 0) return;

    fprintf (fp, "BQ\n");
    q_pos = 0;
    while (q_pos < contig->num_qual_scores) {
        line_pos = 0;
        while (line_pos < 60 && q_pos < contig->num_qual_scores) {
            if (contig->consensus_seq[q_pos] != '*') {
                fprintf (fp, "%d ", contig->qual_scores[q_pos]);
                line_pos++;
            }
            q_pos++;
        }
        fprintf (fp, "\n");
    }
}


/* NOTE - this file does not provide all of the information required for an ACE file. */
static void s_WriteContig (FILE *fp, TContigPtr contig)
{
    int i;

    if (contig == NULL) return;

    fprintf (fp, "CO %s %d %d\n\n", contig->consensus_id, contig->consensus_assem_len, contig->num_reads);
    s_WriteSeq (fp, contig->consensus_seq);
    fprintf (fp, "\n");

    s_WriteQualScores (fp, contig);
    fprintf (fp, "\n");

    for (i = 0; i < contig->num_reads; i++) {
        fprintf (fp, "AF %s %c %d\n", contig->reads[i]->read_id,
                                      contig->reads[i]->is_complement ? 'C' : 'U',
                                      contig->reads[i]->cons_start + 1);
    }
    fprintf (fp, "\n");
    for (i = 0; i < contig->num_reads; i++) {
        fprintf (fp, "RD %s %d\n", contig->reads[i]->read_id, strlen (contig->reads[i]->read_seq));
        s_WriteSeq (fp, contig->reads[i]->read_seq);
        fprintf (fp, "\n");
    }
         
}


/* NOTE - This generates an incomplete ACE file - the data structure we currently use
 * does not provide enough data to create a complete ACE file.
 */
extern void WriteACEFile (FILE *fp, TACEFilePtr afp)
{
  int i, tot_reads = 0;
  if (fp == NULL || afp == NULL) return;

  for (i = 0; i < afp->num_contigs; i++) {
    tot_reads += afp->contigs[i]->num_reads;
  }
  fprintf (fp, "AS %d %d\n\n", afp->num_contigs, tot_reads);

  for (i = 0; i < afp->num_contigs; i++) {
    s_WriteContig (fp, afp->contigs[i]);
  }   
}


/* This function generates a string that uses the FASTA+GAP method for specifying gaps
 * (dashes instead of asterisks)
 */
static char * 
s_AlignmentSeqFromContigSeq 
(char *contig_seq,
 int   cons_start,
 int   aln_len,
 int   read_start,
 int   read_stop)
{
    char * aln_seq;
    char * cp;
    int  pos = 0, i;

    aln_seq = malloc (sizeof (char) * (aln_len + 1));
    /* pad start */
    for (i = 0; i < cons_start; i++) {
        aln_seq[pos] = '-';
        pos++;
    }
    cp = contig_seq;
    if (read_start > 1) {
        i = 1;
        while (*cp != 0 && i < read_start) {
            aln_seq[pos] = '-';
            pos++;
            cp++;
            i++;
        }
    }
    while (*cp != 0 && (read_stop < 1 || i < read_stop)) {
        if (*cp == '*') {
            aln_seq[pos] = '-';
        } else {
            aln_seq[pos] = *cp;
        }
        pos++;
        cp++;
        i++;
    }
    while (pos < aln_len) {
        aln_seq[pos] = '-';
        pos++;
    }
    aln_seq[pos] = 0;
    return aln_seq;
}


static char * s_FarPointerIdFromReadId (char * read_id)
{
    char * far_id = NULL;
    far_id = malloc (sizeof (char) * (strlen (read_id) + 1));
    strcpy (far_id, read_id);
    return far_id;
}


/* This function generates an intermediate data format suitable for generating
 * a SeqEntry with an alignment.
 */
extern TAlignmentFilePtr AlignmentFileFromContig (TContigPtr contig)
{
    TAlignmentFilePtr aln;
    int               i, len;
    int               consensus_pad = 0, pad, end_pad = 0, aln_len;

    if (contig == NULL) return NULL;

    aln = AlignmentFileNew ();
    aln->num_sequences = 1 + contig->num_reads;
    aln->num_organisms = 0;
    aln->num_deflines = 0;
    aln->num_segments = 1;
    aln->ids = malloc (sizeof (char *) * aln->num_sequences);
    aln->sequences = malloc (sizeof (char *) * aln->num_sequences);
    aln->organisms = NULL;
    aln->deflines = NULL;
    aln->align_format_found = 1;
    /* calculate padding for consensus */
    for (i = 0; i < contig->num_reads; i++) {
        if (contig->reads[i]->cons_start < 0) {
            pad = 0 - contig->reads[i]->cons_start;
            if (consensus_pad < pad) {
                consensus_pad = pad;
            }
        }
        len = contig->reads[i]->cons_start + strlen (contig->reads[i]->read_seq);
        if (len > contig->consensus_assem_len) {
            pad = len - contig->consensus_assem_len;
            if (pad > end_pad) {
                end_pad = pad;
            }
        }
    }
    aln_len = consensus_pad + contig->consensus_assem_len + end_pad;
    /* seq for consensus */
    len = strlen (contig->consensus_id);
    aln->ids[0] = malloc (sizeof (char) * (len + 1));
    strcpy (aln->ids[0], contig->consensus_id);
    aln->sequences[0] = s_AlignmentSeqFromContigSeq (contig->consensus_seq,
                                                     consensus_pad,
                                                     aln_len, 0, 0);
    for (i = 1; i < aln->num_sequences; i++) {
        len = strlen (contig->reads[i - 1]->read_id);
        aln->ids[i] = s_FarPointerIdFromReadId (contig->reads[i - 1]->read_id);
        aln->sequences[i] = s_AlignmentSeqFromContigSeq (contig->reads[i - 1]->read_seq,
                                                         consensus_pad + contig->reads[i - 1]->cons_start,
                                                         aln_len,
                                                         contig->reads[i - 1]->read_assem_start,
                                                         contig->reads[i - 1]->read_assem_stop);
    }
    return aln;
}


/* The Trace Archive Gap String is a list of the number of nucleotides to skip before adding the next gap */
extern char * TraceArchiveGapStringFromACESequence (char *seq_str)
{ 
    char *cp;
    char * gap_str = NULL;
    char * print_pos;
    int    len = 0, pos;

    if (seq_str == NULL) return NULL;

    /* first determine length of gap string */
    cp = seq_str;
    while (*cp != 0) {
        if (*cp == '*' || *cp == '-') {
            len++;
        }
        cp++;
    }
    len = 15 * len + 1;
    gap_str = malloc (sizeof (char) * len);
    cp = seq_str;
    print_pos = gap_str;
    pos = 0;
    while (*cp != 0) {
        if (*cp == '*' || *cp == '-') {
            sprintf (print_pos, "%d,", pos);
            print_pos += strlen (print_pos);
            pos = 0;
        } else {
            pos++;
        }
        cp++;
    }
    /* trim final comma */
    print_pos[strlen(print_pos) - 1] = 0;
    return gap_str;
}


/* NOTE - These functions are currently incomplete */
extern void WriteTraceArchiveRead (FILE *fp, TContigReadPtr read)
{
    char *cp;
    if (fp == NULL || read == NULL) {
        return;
    }

    fprintf (fp, "<trace>\n");
    fprintf (fp, "<trace_name>%s</trace_name>\n", read->read_id);
    fprintf (fp, "<traceconsensus>");
    cp = read->read_seq;
    while (*cp != 0) {
        if (*cp != '*') {
            fprintf (fp, "%c", *cp);
        }
        cp++;
    }
    fprintf (fp, "</traceconsensus>\n"); 
    cp = TraceArchiveGapStringFromACESequence (read->read_seq);
    fprintf (fp, "<tracegaps>%s</tracegaps>\n", cp);
    free (cp);
    fprintf (fp, "</trace>\n");
}


static int s_GetTokenLen (char *str)
{
    char *cp;
    int   len = 0;

    if (str == NULL) return 0;

    cp = str;
    while (*cp != 0 && !isspace (*cp)) {
        len++;
        cp++;
    }
    return len;
}
 
   
static char * s_SkipTokens (char *str, int num_tokens)
{
    char *cp;
    int   i;

    if (str == NULL || num_tokens < 0) return NULL;

    cp = str;
    /* skip leading whitespace */
    while (isspace (*cp)) {
        cp++;
    }

    for (i = 0; i < num_tokens && *cp != 0; i++) {
        /* skip token */
        while (*cp != 0 && !isspace (*cp)) {
            cp++;
        }
        /* skip trailing whitespace */
        while (isspace (*cp)) {
            cp++;
        }
    }
    if (*cp == 0) {
        return NULL;
    } else {
        return cp;
    }
}
    

/* for reading other formats */
extern TContigReadPtr 
ReadContigFromString 
(char *str,
 char **consensus_id,
 int    id_col,
 int    seq_col, 
 int    contig_id_col,
 int    strand_col,
 int    start_col,
 int    interpret_n_col
 )
{
    TContigReadPtr read = NULL;
    char *cp;
    int len, col_num = 1, n_is_gap = 0;
    int max_col;

    if (str == NULL) {
        return NULL;
    }

    max_col = id_col;
    if (seq_col > max_col) {
      max_col = seq_col;
    }
    if (contig_id_col > max_col) {
      max_col = contig_id_col;
    }
    if (strand_col > max_col) {
      max_col = strand_col;
    }
    if (start_col > max_col) {
      max_col = start_col;
    }
    if (interpret_n_col > max_col) {
      max_col = interpret_n_col;
    }

    read = ContigReadNew ();

    cp = str;
    len = s_GetTokenLen (cp);
    while (cp != NULL && *cp != 0 && col_num <= max_col) {
        if (id_col == col_num) {
            read->read_id = malloc (len + 1);
            strncpy (read->read_id, cp, len);
            read->read_id[len] = 0;
        } else if (seq_col == col_num) {
            read->read_seq = malloc (len + 1);
            strncpy (read->read_seq, cp, len);
            read->read_seq[len] = 0;
            read->read_len = len;
        } else if (contig_id_col == col_num) {
            if (consensus_id != NULL) {
                *consensus_id = malloc (len + 1);
                strncpy (*consensus_id, cp, len);
                (*consensus_id)[len] = 0;
            }
        } else if (strand_col == col_num) {
            if (*cp == 'R' || *cp == '-') {
                read->is_complement = 1;
            }
        } else if (start_col == col_num) {
            read->cons_start = atoi (cp);
        } else if (interpret_n_col == col_num) {
            if (*cp == 'I') {
                n_is_gap = 1;
            }
        }
        /* advance to next token */
        col_num++;
        cp = s_SkipTokens (cp, 1);
        len = s_GetTokenLen (cp);
    }
            
    if (max_col > col_num) {
        ContigReadFree (read);
        read = NULL;
    } else {
        read->cons_stop = read->cons_start + read->read_len - 1;
        read->tiling_start = read->cons_start;
        read->tiling_stop = read->cons_stop;
        read->read_assem_start = 0;
        read->read_assem_stop = read->read_len - 1;
        read->read_start = 1;
        read->read_stop = read->read_len;
        if (n_is_gap) {
            /* adjust for gaps */
            read->gaps = GapInfoFromSequenceString (read->read_seq, "N");
            if (read->gaps->num_gaps > 0) {
                RemoveGapCharsFromSequenceString (read->read_seq, "N");
                read->read_stop -= read->gaps->num_gaps;
                read->read_len -= read->gaps->num_gaps;
            }
        }
    }
         
    return read;
}


extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromMAQString (char *str, char **consensus_id)
{
    TContigReadPtr read = NULL;

    read = ReadContigFromString (str, consensus_id, 1, 15, 2, 4, 3, 0);
    return read;
}


extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandMostCompressed (char *str, char **consensus_id)
{
    TContigReadPtr read = NULL;

    read = ReadContigFromString (str, consensus_id, 0, 1, 0, 5, 4, 0);
    return read;
}


extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandSanger (char *str, char **consensus_id)
{
    TContigReadPtr read = NULL;

    read = ReadContigFromString (str, consensus_id, 1, 2, 4, 6, 5, 0);
    return read;
}


extern TContigReadPtr ASSEMBLY_CALLBACK ReadFromElandStandalone (char *str, char **consensus_id)
{
    TContigReadPtr read = NULL;

    read = ReadContigFromString (str, consensus_id, 1, 2, 7, 9, 8, 10);
    return read;
}
  

#define READ_BLOCK_SIZE 50

typedef struct ReadList {
  TContigReadPtr reads[READ_BLOCK_SIZE];
  int         num_reads;
  struct ReadList * next;
} SReadList, * TReadListPtr;


static TReadListPtr ReadListNew ()
{
    TReadListPtr r;

    r = malloc (sizeof (SReadList));
    r->num_reads = 0;
    r->next = NULL;
    return r;
}


static TReadListPtr ReadListFree (TReadListPtr r)
{
    TReadListPtr r_next;
    int          i;

    while (r != NULL) {
        for (i = 0; i < r->num_reads; i++) {
            ContigReadFree (r->reads[i]);
        }
        r_next = r;
        free (r);
        r = r_next;
    }
    return r;
}


static TReadListPtr AddToReadList (TContigReadPtr read, TReadListPtr read_list)
{
    if (read_list == NULL) {
        read_list = ReadListNew();
    } else {
        while (read_list->next != NULL && read_list->num_reads == READ_BLOCK_SIZE) {
            read_list = read_list->next;
        }
        if (read_list->num_reads == READ_BLOCK_SIZE) {
            read_list->next = ReadListNew();
            read_list = read_list->next;
        }
    }
    read_list->reads[read_list->num_reads++] = read;
    return read_list;
}

   
typedef struct ConsensusReads {
  TContigPtr   contig;
  TReadListPtr read_list;
  TReadListPtr last_read;
  struct ConsensusReads * next;
} SConsensusReads, * TConsensusReadsPtr;


static TConsensusReadsPtr ConsensusReadsNew (char *consensus_id)
{
    TConsensusReadsPtr c;

    c = malloc (sizeof (SConsensusReads));
    c->contig = ContigNew ();
    if (consensus_id != NULL) {
        c->contig->consensus_id = malloc (strlen (consensus_id) + 1);
        strcpy (c->contig->consensus_id, consensus_id);
    }

    c->read_list = NULL;
    c->last_read = NULL;
    c->next = NULL;
    return c;
}


static TConsensusReadsPtr ConsensusReadsFree (TConsensusReadsPtr c)
{
    TConsensusReadsPtr c_next;

    while (c != NULL) {
        c_next = c->next;
        ContigFree (c->contig);
        c->read_list = ReadListFree (c->read_list);
        free (c);
        c = c_next;
    }
    return c;
}


static void AddReadToConsensusReads (TConsensusReadsPtr c, TContigReadPtr read)
{
    if (c != NULL && read != NULL) {
        c->last_read = AddToReadList (read, c->read_list);
        if (c->read_list == NULL) {
            c->read_list = c->last_read;
        }
    }
}

#define CONSENSUS_BLOCK_SIZE 50

typedef struct ConsensusReadsList {
  TConsensusReadsPtr contigs[CONSENSUS_BLOCK_SIZE];
  int                num_contigs;
  struct ConsensusReadsList * next;
} SConsensusReadsList, * TConsensusReadsListPtr;


static TConsensusReadsListPtr ConsensusReadsListNew ()
{
    TConsensusReadsListPtr c;

    c = malloc (sizeof (SConsensusReadsList));
    c->num_contigs = 0;
    c->next = NULL;
    return c;
}


static TConsensusReadsListPtr ConsensusReadsListFree (TConsensusReadsListPtr c)
{
    TConsensusReadsListPtr c_next;
    int                    i;

    while (c != NULL) {
        c_next = c->next;
        for (i = 0; i < c->num_contigs; i++) {
             c->contigs[i] = ConsensusReadsFree (c->contigs[i]);
        }
        free (c);
        c = c_next;
    }
    return c;
}


static TConsensusReadsPtr FindConsensusIDInConsensusReadsList (TConsensusReadsListPtr c, char *consensus_id)
{
    int i;
    TConsensusReadsPtr r = NULL;

    if (consensus_id == NULL) {
        return NULL;
    }
    while (c != NULL && r == NULL)  {
        for (i = 0; i < c->num_contigs && r == NULL; i++) {
            if (c->contigs[i] != NULL
                && c->contigs[i]->contig != NULL
                && strcmp (c->contigs[i]->contig->consensus_id, consensus_id) == 0) {
                r = c->contigs[i];
            }
        }
        c = c->next;
    }
    return r;   
}


static TConsensusReadsListPtr 
AddConsensusReadToConsensusReadsList 
(TConsensusReadsListPtr c,
 char                 * consensus_id,
 TContigReadPtr         read)
{
    TConsensusReadsPtr r = NULL;

    if (c == NULL) {
        c = ConsensusReadsListNew ();
        r = ConsensusReadsNew(consensus_id);
        c->contigs[c->num_contigs++] = r;
    } else {
        r = FindConsensusIDInConsensusReadsList (c, consensus_id);
        if (r == NULL) {
            while (c->next != NULL && c->num_contigs == CONSENSUS_BLOCK_SIZE) {
                c = c->next;
            }
            if (c->num_contigs == CONSENSUS_BLOCK_SIZE) {
                c->next = ConsensusReadsListNew ();
                c = c->next;
            }
            r = ConsensusReadsNew(consensus_id);
            c->contigs[c->num_contigs++] = r;
        }
    }
    AddReadToConsensusReads (r, read);

    return c;
}


static void MoveReadsToContigFromReadList (TContigPtr contig, TReadListPtr read_list)
{
    TReadListPtr r;
    int          n = 0, i;

    if (contig == NULL) {
        return;
    }

    for (r = read_list; r != NULL; r = r->next) {
        n += r->num_reads;
    }

    contig->num_reads = n;
    contig->reads = malloc (contig->num_reads * sizeof (TContigReadPtr));
    n = 0;

    for (r = read_list; r != NULL; r = r->next) {
        for (i = 0; i < r->num_reads; i++) {
            contig->reads[n++] = r->reads[i];
            r->reads[i] = NULL;
        }
        r->num_reads = 0;
    }
}
   

static TACEFilePtr ACEFileFromConsensusReadsList (TConsensusReadsListPtr contig_list)
{
    TACEFilePtr afp = NULL;
    TConsensusReadsListPtr c;
    int                    i, n = 0;
    
    if (contig_list == NULL || contig_list->num_contigs == 0) {
        return NULL;
    }

    afp = ACEFileNew ();
    for (c = contig_list; c != NULL; c=c->next) {
        afp->num_contigs += c->num_contigs;
    }

    afp->contigs = malloc (afp->num_contigs * sizeof (TContigPtr));
    for (c = contig_list; c != NULL; c = c->next) {
        for (i = 0; i < c->num_contigs; i++) {
            MoveReadsToContigFromReadList (c->contigs[i]->contig, c->contigs[i]->read_list);
            afp->contigs[n++] = c->contigs[i]->contig;
            c->contigs[i]->contig = NULL;
            c->contigs[i]->last_read = NULL;
        }
        c->num_contigs = 0;
    }
    return afp;
}


extern TACEFilePtr ReadAssemblyFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata,  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
 FReadFromStringFunction makeread_func) /* function to transform a string into a read */
{
    TACEFilePtr afp = NULL;
    TContigReadPtr read;
    TConsensusReadsListPtr contig_list = NULL, contig_last = NULL;
    char *linestring;
    char *consensus_id = NULL;

    if (readfunc == NULL || makeread_func == NULL) {
        return NULL;
    }
    linestring = readfunc (fileuserdata);

    while (linestring != NULL  &&  linestring [0] != EOF) {
        /* get ContigRead */
        read = makeread_func (linestring, &consensus_id);
        
        /* group with other ContigReads from the same consensus_id */
        if (read != NULL && consensus_id != NULL) {
            contig_last = AddConsensusReadToConsensusReadsList (contig_last, consensus_id, read);
            if (contig_list == NULL) {
                contig_list = contig_last;
            }
            read = NULL;
        } 
        if (consensus_id != NULL) {
            free (consensus_id);
            consensus_id = NULL;
        }
        ContigReadFree (read);
        free (linestring);
        linestring = readfunc (fileuserdata);
    }

    afp = ACEFileFromConsensusReadsList (contig_list);
    contig_list = ConsensusReadsListFree (contig_list);
    return afp;
}


extern TACEFilePtr ReadMAQFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata)  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
{
    return ReadAssemblyFile (readfunc, fileuserdata, ReadFromMAQString);
}


extern TACEFilePtr ReadElandStandaloneFile 
(FReadLineFunction    readfunc,      /* function for reading lines of 
                                       * alignment file
                                       */
 void *               fileuserdata)  /* data to be passed back each time
                                       * readfunc is invoked
                                       */
{
    return ReadAssemblyFile (readfunc, fileuserdata, ReadFromElandStandalone);
}


/* functions for writing out XML */
static void WriteTraceGapsXML (TGapInfoPtr gap_info, FILE *fp)
{
  int i;

  if (gap_info != NULL && fp != NULL) {
    fprintf (fp, "    <ntracegaps>%d</ntracegaps>\n", gap_info->num_gaps);
    if (gap_info->num_gaps > 0) {
      fprintf (fp, "    <tracegaps source=\"INLINE\">");
      for (i = 0; i < gap_info->num_gaps - 1; i++) {
        fprintf (fp, "%d ", gap_info->gap_offsets[i]);
      }
      fprintf (fp, "%d</tracegaps>\n", gap_info->gap_offsets[gap_info->num_gaps - 1]);
    } else {
      fprintf (fp, "    <tracegaps source=\"INLINE\"> </tracegaps>\n");
    }
  }
}


static void WriteTraceReadXML (TContigReadPtr read, FILE *fp)
{
  if (read != NULL && fp != NULL) {
    fprintf (fp, "<trace>\n");
    if (read->ti > 0) {
      fprintf (fp, "  <ti>%d</ti>\n", read->ti);
    }
    if (read->srr != NULL) {
      fprintf (fp, "  <srr>%s</srr>\n", read->srr);
    }
    if (read->read_id != NULL) {
      fprintf (fp, "  <trace_name>%s</trace_name>\n", read->read_id);
    }
    fprintf (fp, "  <nbasecalls>%d</nbasecalls>\n", read->read_len);
    fprintf (fp, "  <valid>\n");
    fprintf (fp, "    <start>%d</start>\n", read->read_start);
    fprintf (fp, "    <stop>%d</stop>\n", read->read_stop);
    fprintf (fp, "  </valid>\n");
    fprintf (fp, "  <tiling direction = \"%s\">\n", read->is_complement ? "REVERSE" : "FORWARD");
    fprintf (fp, "    <start>%d</start>\n", read->tiling_start);
    fprintf (fp, "    <stop>%d</stop>\n", read->tiling_stop);
    fprintf (fp, "  </tiling>\n");
    fprintf (fp, "  <traceconsensus>\n");
    fprintf (fp, "    <start>%d</start>\n", read->cons_start);
    fprintf (fp, "    <stop>%d</stop>\n", read->cons_stop);
    fprintf (fp, "  </traceconsensus>\n");
    WriteTraceGapsXML (read->gaps, fp);
    fprintf (fp, "</trace>\n");
  }
}


extern void WriteTraceAssemblyFromContig (TContigPtr contig, FILE *fp)
{
  int i;

  if (contig == NULL || fp == NULL) return;

  /* NOTE - need to add new field to TContigPtr for submitter reference, where orig ID should move to */
  fprintf (fp, "  <contig submitter_reference=\"%s\" conformation=\"LINEAR\" type=\"NEW\">\n", 
           contig->consensus_id == NULL ? "not supplied" : contig->consensus_id);

  fprintf (fp, "    <ntraces>%d</ntraces>\n", contig->num_reads);

  fprintf (fp, "    <nconbases>%d</nconbases>\n", contig->consensus_seq_len);

  /* need nbasecalls */

  if (contig->gaps == NULL) {
    fprintf (fp, "    <ncongaps>0</ncongaps>\n");
  } else {
    fprintf (fp, "    <ncongaps>%d</ncongaps>\n", contig->gaps->num_gaps);
    if (contig->gaps->num_gaps > 0) {
      fprintf (fp, "  <congaps source=\"INLINE\">");
      for (i = 0; i < contig->gaps->num_gaps - 1; i++) {
        fprintf (fp, "%d ", contig->gaps->gap_offsets[i]);
      }
      fprintf (fp, "%d</congaps>\n", contig->gaps->gap_offsets[contig->gaps->num_gaps - 1]);
    }
  }
  fprintf (fp, "    <consensus>%s</consensus>\n", 
           contig->consensus_id == NULL ? "not supplied" : contig->consensus_id);
  if (contig->num_qual_scores > 0) {
    fprintf (fp, "    <conqualities source=\"INLINE\">");
    for (i = 0; i < contig->num_qual_scores; i++) {
      fprintf (fp, "%d ", contig->qual_scores[i]);
    }
    fprintf (fp, "</conqualities>\n");
  }
  
  for (i = 0; i < contig->num_reads; i++) {
    WriteTraceReadXML (contig->reads[i], fp);
  }
  fprintf (fp, "  </contig>\n");
}


extern void
WriteTraceAssemblyHeader
(char * assembly_type,
 char * subref,
 char * center_name,
 int    taxid,
 char * description,
 char * assembly,
 int    num_contigs,
 unsigned int    num_conbases,
 int    num_reads,
 unsigned int    num_readbases,
 FILE * fp)
{
    if (fp == NULL) {
        return;
    }

    fprintf (fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

    fprintf (fp, "<assembly submitter_reference=\"%s\" type = \"%s\">\n", 
                 subref == NULL ? "Not supplied" : subref,
                 assembly_type == NULL ? "NEW" : assembly_type);
    fprintf (fp, "  <center_name>%s</center_name>\n", center_name == NULL ? "Not supplied" : center_name);\
    fprintf (fp, "  <organism descriptor=\"TAXID\">%d</organism>\n", taxid);
    fprintf (fp, "  <description>%s</description>\n", description == NULL ? "Not supplied" : description);
    fprintf (fp, "  <structure>%s</structure>\n", assembly == NULL ? "transcript assembly" : assembly);
    fprintf (fp, "  <ncontigs>%d</ncontigs>\n", num_contigs);
    fprintf (fp, "  <nconbases>%u</nconbases>\n", num_conbases);
    fprintf (fp, "  <ntraces>%d</ntraces>\n", num_reads);
    fprintf (fp, "  <nbasecalls>%u</nbasecalls>\n", num_readbases);
    fprintf (fp, "  <coverage>%f</coverage>\n", num_conbases == 0 ? 0 : (float) ((float) num_readbases/ (float) num_conbases));
}
 

extern void WriteTraceAssemblyTrailer (FILE *fp)
{
    if (fp == NULL) {
        return;
    }
    fprintf (fp, "</assembly>\n");
}


extern void 
WriteTraceAssemblyFromAceFile 
(TACEFilePtr afp,
 char *      subref,
 char *      center_name, 
 int         taxid,
 char *      description,
 FILE        *fp)
{ 
  int i, j, traces = 0;
  unsigned int conbases = 0, basecalls = 0;

  if (afp == NULL || fp == NULL) return;


  for (i = 0; i < afp->num_contigs; i++) {
    conbases += afp->contigs[i]->consensus_seq_len;
    traces += afp->contigs[i]->num_reads;
    for (j = 0; j < afp->contigs[i]->num_reads; j++) {
      basecalls += afp->contigs[i]->reads[j]->read_len;
    }
  }
  WriteTraceAssemblyHeader (NULL, subref, center_name, taxid, description, NULL, afp->num_contigs, conbases, traces, basecalls, fp);

  for (i = 0; i < afp->num_contigs; i++) {
    WriteTraceAssemblyFromContig (afp->contigs[i], fp);
  }
  WriteTraceAssemblyTrailer (fp);
}


extern void WriteFASTAFromContig
(TContigPtr contig,
 FILE       *fp)
{
    int   k;
    char *cp;

    if (contig == NULL || fp == NULL) return;
    
    fprintf (fp, ">%s\n", contig->consensus_id);
    cp = contig->consensus_seq;
    while (*cp != 0) {
        k = 0;
        while (k < 40 && *cp != 0) {
            if (*cp != '*') {
                fprintf (fp, "%c", *cp);
                k++;
            }
            cp++;
        }
        fprintf (fp, "\n");
    }
    fprintf (fp, "\n");
}


extern void
WriteFASTAFromAceFile
(TACEFilePtr afp,
 FILE        *fp)
{
  int i;

  if (afp == NULL || fp == NULL) return;
  
  for (i = 0; i < afp->num_contigs; i++) {
    WriteFASTAFromContig (afp->contigs[i], fp);
  }
}


#define kFASTASeqBufSize 100

typedef struct fastaseqbuf {
  char buf[kFASTASeqBufSize];
  int  num_used;
  struct fastaseqbuf *next;
} SFASTASeqBuf, * TFASTASeqBufPtr;


static TFASTASeqBufPtr s_FASTASeqBufNew ()
{
    TFASTASeqBufPtr s;

    s = (TFASTASeqBufPtr) malloc (sizeof (SFASTASeqBuf));
    if (s != NULL) {
        s->num_used = 0;
        s->next = NULL;
    }
    return s;
}


static void s_FASTASeqBufFree (TFASTASeqBufPtr s)
{
    TFASTASeqBufPtr s_next;

    while (s != NULL) {
        s_next = s->next;
        free (s);
        s = s_next;
    }
}


static TFASTASeqBufPtr s_AddFASTAToBuf (char *line, TFASTASeqBufPtr buf)
{
    TFASTASeqBufPtr last_buf;
    char *cp;

    if (buf == NULL) {
        buf = s_FASTASeqBufNew();
        last_buf = buf;
    } else {
        last_buf = buf;
        while (last_buf->next != NULL) {
            last_buf = last_buf->next;
        }
    }

    cp = line;
    while (*cp != 0 && *cp != '\r' && *cp != '\n') {
        while (isspace (*cp)) {
            cp++;
        }
        if (*cp != 0) {
            if (!isalpha (*cp)) {
                printf ("Found bad character in FASTA file!\n");
                s_FASTASeqBufFree (buf);
                buf = NULL;
                return buf;
            }
            if (last_buf->num_used == kFASTASeqBufSize) {
                last_buf->next = s_FASTASeqBufNew ();
                last_buf = last_buf->next;
            }
            last_buf->buf[last_buf->num_used++] = *cp;
            cp++;
        }
    }
    return buf;
}


static char * s_StripStars (char *str)
{
    char *cp_src;
    char *cp_dst;
    char *stripped;

    if (str == NULL) {
        return 0;
    }
    cp_src = str;
    stripped = (char *) malloc (sizeof (char) * (strlen (str) + 1));
    cp_dst = stripped;
    while (*cp_src != 0) {
        if (*cp_src != '*') {
            *cp_dst = *cp_src;
            cp_dst++;
        }
        cp_src++;
    }
    *cp_dst = 0;
    return stripped;
}


static int s_DoesFASTAMatchSeq (TFASTASeqBufPtr buf, char *trimmed_seq)
{
    int does_match = 1, match_len, seq_len;

    if (buf == NULL || trimmed_seq == NULL || *trimmed_seq == 0) {
        return 1;
    }

    seq_len = strlen (trimmed_seq);
    while (buf != NULL && seq_len > 0 && does_match) {
        if (seq_len < buf->num_used) {
            match_len = seq_len;
        } else {
            match_len = buf->num_used;
        }
        if (strncmp (buf->buf, trimmed_seq, match_len) == 0) {
            buf = buf->next;
            trimmed_seq += match_len;
            seq_len -= match_len;
        } else {
            does_match = 0;
        }
    }
    return does_match;
}


static char s_CompLetter (char ch)
{
    switch (ch) {
        case 'A':
            ch = 'T';
            break;
        case 'T':
            ch = 'A';
            break;
        case 'G':
            ch = 'C';
            break;
        case 'C':
            ch = 'G';
            break;
    }
    return ch;
}


static void s_RevCompSequence (char *seq)
{
    char tmp;
    int len, i;

    if (seq == NULL || *seq == 0) {
        return;
    }
    len = strlen (seq);

    for (i = 0; i < len / 2; i++) {
        tmp = seq[i];
        seq[i] = s_CompLetter (seq[len - i - 1]);
        seq[len - i - 1] = s_CompLetter (tmp);
    }
    if (len %2 > 0) {
        seq[i] = s_CompLetter (seq[i]);
    }
}


static int s_GetSequenceOffset (TFASTASeqBufPtr buf, char *trimmed_seq, char is_complement)
{
    int offset = 0, buf_offset;
    int match_found = 0;
    int match_len, seq_len;

    if (buf == NULL || trimmed_seq == NULL) {
        return 0;
    }

    trimmed_seq = s_StripStars (trimmed_seq);

    if (is_complement) {
        s_RevCompSequence (trimmed_seq);
    }

    seq_len = strlen (trimmed_seq);

    while (buf != NULL && !match_found) {
        buf_offset = 0;
        while (buf_offset < buf->num_used && !match_found) {
            if (buf->num_used - buf_offset < seq_len) {
                match_len = buf->num_used - buf_offset;
            } else {
                match_len = seq_len;
            }
            if (match_len < seq_len && buf->next == NULL) {
                /* ran out of sequence, no match */
                buf_offset = buf->num_used;
            } else if (strncmp (buf->buf + buf_offset, trimmed_seq, match_len) == 0
                && s_DoesFASTAMatchSeq (buf->next, trimmed_seq + match_len)) {
                match_found = 1;
            } else {
                buf_offset++;
                offset++;
            }
        }
        if (!match_found) {
            buf = buf->next;
        }
    }
    free (trimmed_seq);
    if (match_found) {
        return offset;
    } else {
        return -1;
    }
}


#define kSeqListBufSize 100

typedef struct seqlistbuf {
  TFASTASeqBufPtr buf[kSeqListBufSize];
  char * id_list[kSeqListBufSize];
  int  num_used;
  struct seqlistbuf *next;
} SSeqListBuf, * TSeqListBufPtr;


static TSeqListBufPtr s_SeqListBufNew ()
{
    TSeqListBufPtr s;

    s = (TSeqListBufPtr) malloc (sizeof (SSeqListBuf));
    if (s != NULL) {
        s->num_used = 0;
        s->next = NULL;
    }
    return s;
}


static void s_SeqListBufFree (TSeqListBufPtr s)
{
    TSeqListBufPtr s_next;
    int i;

    while (s != NULL) {
        s_next = s->next;
        for (i = 0; i < s->num_used; i++) {
            s_FASTASeqBufFree (s->buf[i]);
            free (s->id_list[i]);
            s->buf[i] = NULL;
        }
        free (s);
        s = s_next;
    }
}


static char * s_GetFASTAIdFromString (char * str)
{
    char * cp;
    char * id;
    int    len;

    if (str == NULL) {
        return NULL;
    }

    cp = str;
    cp += strspn (str, " >\t");
    len = strcspn (cp, " \t\r\n");
    if (len == 0) {
        return NULL;
    }
    id = (char *)malloc (len + 1);
    if (id == NULL) {
        return NULL;
    }
    strncpy (id, cp, len);
    id [ len ] = 0;
    return id;
}


static TSeqListBufPtr s_AddToSeqList (char *line, TSeqListBufPtr buf)
{
    TSeqListBufPtr last_buf;

    if (line == NULL) {
        return buf;
    }
    if (buf == NULL) {
        buf = s_SeqListBufNew();
        last_buf = buf;
    } else {
        last_buf = buf;
        while (last_buf->next != NULL) {
            last_buf = last_buf->next;
        }
    }

    if (*line == '>') {
        if (last_buf->num_used == kSeqListBufSize) {
            last_buf->next = s_SeqListBufNew();
            last_buf = last_buf->next;
        }
        last_buf->buf[last_buf->num_used] = s_FASTASeqBufNew ();
        last_buf->id_list[last_buf->num_used] = s_GetFASTAIdFromString (line);
        last_buf->num_used++;
    } else if (last_buf->num_used > 0) {
        last_buf->buf[last_buf->num_used - 1] = s_AddFASTAToBuf (line, last_buf->buf[last_buf->num_used - 1]);
    }
        
    return buf;
}


static TFASTASeqBufPtr s_GetFastaSeq (TSeqListBufPtr buf, char *id)
{
    TFASTASeqBufPtr seq = NULL;
    char *cp;
    int i, match_len;

    if (id == NULL) {
        return NULL;
    }

    cp = strchr (id, '.');
    if (cp == NULL) {
        match_len = strlen (id);
    } else {
        match_len = cp - id;
    }

    while (buf != NULL && seq == NULL) {
        for (i = 0; i < buf->num_used && seq == NULL; i++) {
            if (strncmp (id, buf->id_list[i], match_len) == 0) {
                seq = buf->buf[i];
            }
        }
        buf = buf->next;
    }
    return seq;
}


static TSeqListBufPtr s_ReadFastaFile (FReadLineFunction readfunc, void * userdata)
{
    char *linestring;
    TSeqListBufPtr fasta = NULL;

    if (readfunc == NULL) {
        return NULL;
    }
    linestring = readfunc (userdata);
    while (!s_IsEOF (linestring)) {
        fasta = s_AddToSeqList (linestring, fasta);
        free (linestring);
        linestring = readfunc (userdata);
    }
    return fasta;
}




#define kQualScoreBufSize 100

typedef struct qualscorelist {
  int scores[kQualScoreBufSize];
  int num_used;
  struct qualscorelist *next;
} SQualScoreList, * TQualScoreListPtr;


static TQualScoreListPtr s_QualScoreNew ()
{
    TQualScoreListPtr s;

    s =(TQualScoreListPtr) malloc (sizeof (SQualScoreList));
    if (s != NULL) {
        s->num_used = 0;
        s->next = NULL;
    }
    return s;
}


static void s_QualScoreFree (TQualScoreListPtr s)
{
    TQualScoreListPtr s_next;

    while (s != NULL) {
        s_next = s->next;
        free (s);
        s = s_next;
    }
}


static TQualScoreListPtr s_AddQualScores (char *line, TQualScoreListPtr scores)
{
    TQualScoreListPtr last_score;
    char *cp;

    if (scores == NULL) {
        scores = s_QualScoreNew();
        last_score = scores;
    } else {
        last_score = scores;
        while (last_score->next != NULL) {
            last_score = last_score->next;
        }
    }

    cp = line;
    while (*cp != 0 && *cp != '\r' && *cp != '\n') {
        while (isspace (*cp)) {
            cp++;
        }
        if (*cp != 0) {
            if (!isdigit (*cp)) {
                printf ("Found bad character in quality scores file!\n");
                s_QualScoreFree (scores);
                scores = NULL;
                return scores;
            }
            if (last_score->num_used == kQualScoreBufSize) {
                last_score->next = s_QualScoreNew ();
                last_score = last_score->next;
            }
            last_score->scores[last_score->num_used++] = atoi (cp);
            while (isdigit (*cp)) {
                cp++;
            }
        }
    }
    return scores;
}


static int s_AddScoresToRead (TContigReadPtr read, TQualScoreListPtr scores, TSeqListBufPtr fasta)
{
    int score_pos = 0;
    int offset = 0;
    int skip, score_len;
    char *cp;
    int *dst;
    TFASTASeqBufPtr fasta_seq;

    if (read == NULL || scores == NULL) {
        return 0;
    }

    if (fasta == NULL) {
        skip = read->read_start - 1;
    } else {
        fasta_seq = s_GetFastaSeq (fasta, read->read_id);
        if (fasta_seq == NULL) {
            printf ("Unable to locate fasta for %s\n", read->read_id);
            return 0;
        }

        skip = s_GetSequenceOffset (fasta_seq, read->read_seq, read->is_complement);
        if (skip < 0) {
            printf ("ACE read did not match FASTA read for %s\n", read->read_id);
            return 0;
        }
    }

    /* skip over scores before part used in assembly */
    while (scores != NULL && score_pos < skip) {
        if (skip - score_pos < scores->num_used) {
            offset = skip - score_pos;
            score_pos = skip;
        } else if (scores->next == NULL) {
            printf ("Not enough scores read for %s\n", read->read_id);
            return 0;
        } else {
            score_pos += kQualScoreBufSize;
            scores = scores->next;
        }
    }

    score_len = strlen (read->read_seq);
    read->qual_scores = malloc (sizeof (int) * score_len);

    if (read->is_complement) {
        /* need to read scores in reverse direction */
        cp = read->read_seq + score_len - 1;
        dst = read->qual_scores + score_len - 1;
        while (scores != NULL && read->num_qual_scores < score_len) {
            if (*cp == '*') {
                *dst = 0;
            } else {
                *dst = scores->scores[offset];
                offset++;
            }
            cp--;
            dst--;
            read->num_qual_scores++;

            if (offset == kQualScoreBufSize) {
                scores = scores->next;
                offset = 0;
            }
        }
    } else {
        cp = read->read_seq;
        dst = read->qual_scores;
        while (scores != NULL && read->num_qual_scores < score_len) {
            if (*cp == '*') {
                *dst = 0;
            } else {
                *dst = scores->scores[offset];
                offset++;
            }
            cp++;
            dst++;
            read->num_qual_scores++;

            if (offset == kQualScoreBufSize) {
                scores = scores->next;
                offset = 0;
            }
        }
    }
    if (read->num_qual_scores == score_len) {
        return 1;
    } else {
        printf ("Not enough qual scores for %s\n", read->read_id);
        return 0;
    }
}


static int s_AddQualScoresToReadsInAceFile (TACEFilePtr afp, char *id, TQualScoreListPtr scores, TSeqListBufPtr fasta)
{
    int i, j, found = 0, match_len, rval = 1;
    char *cp;

    if (afp == NULL || id == NULL) {
        return 0;
    }

    cp = strchr (id, '.');
    if (cp == NULL) {
        match_len = strlen (id);
    } else {
        match_len = cp - id;
    }

    for (i = 0; i < afp->num_contigs; i++) {
        for (j = 0; j < afp->contigs[i]->num_reads; j++) {
            if (strncmp (afp->contigs[i]->reads[j]->read_id, id, match_len) == 0) {
                found = 1;
                rval &= s_AddScoresToRead (afp->contigs[i]->reads[j], scores, fasta);              
            }
        }
    }
    if (!found) {
        printf ("Unable to locate %s in ACE file\n", id);
        rval = 0;
    }
    return rval;
}


extern int
AddReadQualScores
(TACEFilePtr          afp,
 FReadLineFunction    readfunc,
 void *               userdata,
 FReadLineFunction    fasta_readfunc,
 void *               fasta_userdata)
{
    char *linestring;
    char *score_id = NULL;
    TQualScoreListPtr scores = NULL;
    TSeqListBufPtr    fasta = NULL;
    int               rval = 1;

    if (afp == NULL || readfunc == NULL) {
        return 0;
    }

    if (fasta_readfunc != NULL) {
        fasta = s_ReadFastaFile (fasta_readfunc, fasta_userdata);
        if (fasta == NULL) {
            printf ("Unable to read FASTA file\n");
            return 0;
        }
    }

    linestring = readfunc (userdata);
    while (!s_IsEOF (linestring)) {
        if (linestring[0] == '>') {            
            if (score_id != NULL && scores != NULL) {
                /* add previously read scores to last read */
                if (s_AddQualScoresToReadsInAceFile (afp, score_id, scores, fasta) == 0) {
                    printf ("Failed to add quality scores from %s\n", score_id);
                    rval = 0;
                }
            }
            s_QualScoreFree (scores);
            scores = NULL;
            free (score_id);

            score_id = s_GetFASTAIdFromString (linestring);
        } else if (score_id != NULL) {
            scores = s_AddQualScores (linestring, scores);
        }
        free (linestring);
        linestring = readfunc (userdata);
    }

    /* handle last set of scores read */
    if (score_id != NULL && scores != NULL) {
        /* add previously read scores to last read */
        if (s_AddQualScoresToReadsInAceFile (afp, score_id, scores, fasta) == 0) {
            printf ("Failed to add quality scores from %s\n", score_id);
            rval = 0;
        }
    }
    s_SeqListBufFree (fasta);
    s_QualScoreFree (scores);
    scores = NULL;
    free (score_id);
    return rval;
}


static int s_LetterPos (char ch) {
    int rval = -1;

    switch (ch) {
        case 'A':
            rval = 0;
            break;
        case 'T':
            rval = 1;
            break;
        case 'G':
            rval = 2;
            break;
        case 'C':
            rval = 3;
            break;
        case '*':
            rval = 4;
            break;
    }
    return rval;
}


static int s_GetUngappedPosition (int gapped_pos, char *seq)
{
    int pos = 0, ungapped_pos = 0;
    int gaps_found = 0;
    char *cp;
    
    cp = seq;
    while (*cp != 0 && pos < gapped_pos ) {
        if (*cp == '*') {
            gaps_found++;
        } else {
            ungapped_pos++;
        } 
        cp++;
        pos++;
    }
    return ungapped_pos;
}


static int s_GetQualScoreForReadPos (TContigReadPtr r, int pos)
{
    if (r == NULL || pos < 0) {
        return 0;
    }

    /* note - don't need to get ungapped position because 0s are inserted when qual scores
     * are added to the reads.
     */
    /*pos = s_GetUngappedPosition (pos, r->read_seq); */
    if (pos > r->num_qual_scores) {
        return 0;
    } else {
        return r->qual_scores[pos];
    }
}


extern int ReplaceConsensusSequenceFromTraces (TContigPtr contig, char only_ns)
{
    char * consensus_buf;
    int  * new_qual_scores = NULL;
    int    num_qual_scores = 0;
    int    i, k, best, letter_pos;
    int    char_counts[5];
    char   best_ch, ch;
    int    num_best, sum_best;
    int  * consensus_qual_ptr = NULL;
    int    read_offset, len;
    int    num_change = 0;

    if (contig == NULL) {
        return 0;
    }

    consensus_buf = (char *) malloc (sizeof (char) * (contig->consensus_assem_len + 1));
    if (contig->reads[0]->num_qual_scores > 0) {
        new_qual_scores = (int *) malloc (sizeof (int) * contig->consensus_seq_len);
    }

    consensus_qual_ptr = contig->qual_scores;

    for (i = 0; i < contig->consensus_assem_len; i++) {
        if (only_ns && contig->consensus_seq[i] != 'N') {
            /* just use existing consensus character */
            consensus_buf[i] = contig->consensus_seq[i];
            /* add in qual scores */
            if (consensus_qual_ptr != NULL && new_qual_scores != NULL && contig->consensus_seq[i] != '*') {
                new_qual_scores[num_qual_scores++] = *consensus_qual_ptr;
            }
        } else {
            for (k = 0; k < 5; k++) {
                char_counts[k] = 0;
            }
            best = 0;
            best_ch = 'N';
            for (k = 0; k < contig->num_reads; k++) {
                read_offset = i - contig->reads[k]->cons_start;
                len = strlen (contig->reads[k]->read_seq);
                if (len > read_offset
                    && read_offset >= 0) {
                    ch = toupper (contig->reads[k]->read_seq[read_offset]);
                    letter_pos = s_LetterPos (ch);
                    if (letter_pos > -1) {
                      char_counts[letter_pos]++;
                      if (char_counts[letter_pos] > best
                          || (char_counts[letter_pos] == best && best_ch == '*')) {
                        best_ch = ch;
                        best = char_counts[letter_pos];
                      }
                    }
                }
            }
            if (toupper (consensus_buf[i]) != best_ch) {
                num_change++;
                consensus_buf[i] = best_ch;
                if (best_ch != '*') {
                    /* calculate quality score */
                    if (new_qual_scores != NULL) {
                        sum_best = 0;
                        num_best = 0;
                        for (k = 0; k < contig->num_reads; k++) {
                            if (contig->reads[k]->num_qual_scores > i - contig->reads[k]->cons_start
                                &&  best_ch == toupper (contig->reads[k]->read_seq[i - contig->reads[k]->cons_start])) {
                                num_best ++;
                                sum_best += s_GetQualScoreForReadPos (contig->reads[k], i - contig->reads[k]->cons_start);
                            }
                        }
                        if (num_best == 0) {
                            new_qual_scores[num_qual_scores++] = 0;
                        } else {
                            new_qual_scores[num_qual_scores++] = sum_best / num_best;
                        }
                    }
                }
            }
        }
        if (consensus_qual_ptr != NULL && contig->consensus_seq[i] != '*') {
            consensus_qual_ptr++;
        }
    }
    consensus_buf[i] = 0;

    free (contig->consensus_seq);
    contig->consensus_seq = consensus_buf;
    if (contig->qual_scores != NULL) {
        free (contig->qual_scores);
    }
    contig->qual_scores = new_qual_scores;
    contig->num_qual_scores = num_qual_scores;
    
    return num_change;
}


extern void RecalculateConsensusSequences (TACEFilePtr ace_file, char only_ns)
{
    int i;

    if (ace_file == NULL) {
        return;
    }

    for (i = 0; i < ace_file->num_contigs; i++) {
        ReplaceConsensusSequenceFromTraces(ace_file->contigs[i], only_ns);
    }

}


extern void WriteContigQualScores (TContigPtr contig, FILE *out)
{
    int i = 0, j;

    if (contig == NULL || contig->qual_scores == NULL || contig->num_qual_scores < 1 || out == NULL) {
        return;
    }
    fprintf (out, ">%s\n", contig->consensus_id);

    while (i < contig->num_qual_scores) {
        for (j = 0; j < 60 && i < contig->num_qual_scores; j++, i++) {
            fprintf (out, "%d ", contig->qual_scores[i]);
        }
        fprintf (out, "\n");
    }
    fprintf (out, "\n");
}


extern char
ProcessLargeACEFileForContigFastaAndQualScores
(FReadLineFunction    readfunc,
 void *               userdata,
 char                 make_qual_scores,
 char *               has_errors,
 ProcessContigFunc    process_func,
 void *               process_data)
{
    char *              linestring;
    char *              cp;
    int                 contig_num = 0, read_num = 0;
    int                 num_reads_expected = 0;
    int                 num_contigs = 0;
    TContigPtr          contig = NULL;
    char                rval = 1;

    if (readfunc == NULL || process_func == NULL) {
        return 0;
    }

    linestring = readfunc (userdata);

    while (linestring != NULL  &&  linestring [0] != EOF) {
        if (linestring [0] == 'A' && linestring [1] == 'S' && isspace (linestring [2])) {
            if (num_reads_expected > 0) {
                PrintACEFormatErrorXML ("Two file header lines!", NULL, has_errors);
                return 0;
            }
            /* first line in file, number of contigs */
            cp = linestring + 3;
            num_contigs = atoi (cp);
            while (isdigit (*cp)) {
                cp++;
            }
            num_reads_expected = atoi (cp);
            linestring = readfunc (userdata);
        } else if (linestring [0] == 'C' && linestring [1] == 'O' && isspace (linestring [2])) {
            if (contig_num >= num_contigs) {
                PrintACEFormatErrorXML ("Too many contigs!", NULL, has_errors);
                return 0;
            }

            contig = s_ReadContig (&linestring, readfunc, userdata, make_qual_scores, has_errors);
            if (contig == NULL) {
                PrintACEFormatErrorXMLStart (NULL, has_errors);
                printf ("Unable to read contig (%d)", contig_num);
                PrintACEFormatErrorXMLEnd ();
                return 0;
            }
            read_num += contig->num_reads;
            process_func (contig, process_data);
            ContigFree (contig);
            contig = NULL;
            contig_num++;
        } else if (s_UnexpectedLineBetweenContigs (linestring)) {
            PrintACEFormatErrorXMLStart (NULL, has_errors);
            printf ("Unexpected line after contig %d", read_num);
            PrintACEFormatErrorXMLEnd ();
            return 0;
        } else {
            linestring = readfunc (userdata);
        }
    }
    if (contig_num < num_contigs) {
        PrintACEFormatErrorXML ("Not enough contigs!", NULL, has_errors);
        rval = 0;
    } else if (read_num < num_reads_expected) {
        PrintACEFormatErrorXML ("Not enough reads!", NULL, has_errors);
        rval = 0;
    }

    return rval;
}
