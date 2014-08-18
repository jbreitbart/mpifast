/*   aceread_tst.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  aceread_tst.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   7/22/08
*
* $Revision: 1.25 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <gather.h>
#include <explore.h>
#include <lsqfetch.h>
#include <valid.h>
#include <pmfapi.h>
#ifdef INTERNAL_NCBI_ASNDISC
#include <accpubseq.h>
#include <tax3api.h>
#endif

#include "aceread.h"
#include "acerdapi.h"

typedef enum {
  i_argInputFile,
  o_argOutputFile,
  f_argFASTA,
  S_argIDSubstitutionFile,
  R_argSRRids,
  L_argSuppressIdLookup,
  Q_argMakeQualScores,
  X_argXMLFile,
  t_argTemplateFile,
  T_argTSAFields,
  C_argCenter,
  F_argFormat,
  G_argGapString,
  V_argValidateAgainstAsn1File,
  q_argReadQualScoresFile,
  r_argReadFASTAFile,
  N_argRecalculateConsensus,
  c_argChunkSize,
  n_argReadNameType,
  z_argIncludeReads,
  l_argLimitNumContigs
} EArgNum;

Args myargs [] = {
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"FASTA Output", "F", NULL, NULL,
    TRUE, 'f', ARG_BOOLEAN, 0.0, 0, NULL},
  {"ID Substitution File", "", NULL, NULL,
    TRUE, 'S', ARG_FILE_IN, 0.0, 0, NULL},
  {"Replacement IDs are SRR", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Suppress ID Lookup", "F", NULL, NULL,
    TRUE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Make Qual Scores", "T", NULL, NULL,
    TRUE, 'Q', ARG_BOOLEAN, 0.0, 0, NULL},
  {"XML Output File", "", NULL, NULL,
    TRUE, 'X', ARG_FILE_OUT, 0.0, 0, NULL },
  {"Template File", "", NULL, NULL,
    TRUE, 't', ARG_FILE_IN, 0.0, 0, NULL },
  {"TSA fields", NULL, NULL, NULL,
    TRUE, 'T', ARG_STRING, 0.0, 0, NULL },
  {"Genome Center Tag", NULL, NULL, NULL,
    TRUE, 'C', ARG_STRING, 0.0, 0, NULL},
  {"Assembly Format\n\tM MAQ\n\tE Standalone Eland\n\tA ACE", "A", NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Gap String", NULL, NULL, NULL,
    TRUE, 'G', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 File to validate against", NULL, NULL, NULL,
    TRUE, 'V', ARG_FILE_IN, 0.0, 0, NULL},
  {"Quality score file for read sequences", NULL, NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"FASTA file for read sequences (to use when trimming read quality scores)", NULL, NULL, NULL,
    TRUE, 'r', ARG_FILE_IN, 0.0, 0, NULL},
  {"Recalculate consensus sequence using read data\n\tW Whole Consensus\n\tN Ns Only", "", NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Number of contig bases per file", "50000", NULL, NULL,
    TRUE, 'c', ARG_INT, 0.0, 0, NULL},
  {"Read name type in ACE file\n\tL local trace name\n\tT TI number\n\tS SRR ID\n", "L", NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Include read sequences in ASN.1 output", "F", NULL, NULL,
    TRUE, 'z', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Limit number of contigs to read", NULL, NULL, NULL,
    TRUE, 'l', ARG_INT, 0.0, 0, NULL},
};


static FILE *OpenAceFile (CharPtr infile)
{
  FILE        *f;
  Int4        len;
#ifdef OS_UNIX
  Char            cmmd [256];
  CharPtr         gzcatprog;
  int             ret;
  Boolean         usedPopen = FALSE;
#endif

  len = StringLen (infile);
  if (StringCmp (infile + len - 3, ".gz") == 0) {
#ifdef OS_UNIX
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, infile);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", infile);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return NULL;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", infile);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return NULL;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return NULL;
        }
      }
    }
    f = popen (cmmd, "r");
    usedPopen = TRUE;
#else
    Message (MSG_POSTERR, "Unable to read gzipped files when not running in UNIX");
    return NULL;
#endif
  } else {
    f = FileOpen (infile, "r");
  }
  return f;
}


static Boolean ValidateAgainstASNFile (TACEFilePtr ace_file, CharPtr filename, char *has_errors)
{
  Pointer      dataptr;
  Uint2        datatype;
  SeqEntryPtr  sep = NULL;
  SeqSubmitPtr ssp = NULL;
  Boolean      chars_stripped = FALSE;
  FILE *fp;
  Boolean      rval = FALSE;
  

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    printf ("Unable to open %s\n", filename);
    return FALSE;
  }

  /* Read in one sequence from the file */
  dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, FALSE, FALSE,
		                   	                  TRUE, FALSE, &chars_stripped);      
  FileClose (fp);
  if (NULL == dataptr) 
  {
    printf ("Unable to read SeqEntry from %s\n", filename);
    return FALSE;
  }

  /* Convert the file data to a SeqEntry */
  
  if (datatype == OBJ_SEQENTRY)
    sep = (SeqEntryPtr) dataptr;
  else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
    sep = SeqMgrGetSeqEntryForData (dataptr);
  else if (datatype == OBJ_SEQSUB) 
  {
    ssp = (SeqSubmitPtr) dataptr;
    if (ssp != NULL && ssp->datatype == 1)
    {
      sep = (SeqEntryPtr) ssp->data;
    }
  }
  
  rval = ValidateACEFileAgainstSeqEntry (ace_file, sep, has_errors);

  if (ssp != NULL) {
    ssp = SeqSubmitFree (ssp);
  } else {
    sep = SeqEntryFree (sep);
  }
  return rval;
 
}


static Boolean StringNHasNoText (CharPtr str, Int4 n)
{
  CharPtr cp;
  Int4    i;
  if (str == NULL) return TRUE;
  cp = str;
  i = 0;
  while (i < n) {
    if (*cp == 0) return TRUE;
    if (!isspace (*cp)) return FALSE;
    cp++;
    i++;
  }
  return TRUE;
}


static Boolean BracketMatchesLabel (CharPtr cp, CharPtr cp_equal, CharPtr label) 
{
  Int4 len;

  if (cp == NULL || cp_equal == NULL || label == NULL) return FALSE;

  len = StringLen (label);
  if (StringNCmp (cp, label, len) == 0 
        && StringNHasNoText (cp + len, cp_equal - cp - len)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr GetBracketValue (CharPtr cp, CharPtr cp_end)
{
  Int4 len;
  CharPtr val = NULL;

  if (cp == NULL || cp_end == NULL || cp_end <= cp) return NULL;

  cp += StringSpn (cp, " \t");
  len = (cp_end - cp) + 1;
  val = (CharPtr) MemNew (sizeof (Char) * len);
  StringNCpy (val, cp, len - 1); 
  val [len] = 0;
  while (len > 1 && isspace (val [len-1])) {
    len--;
    val[len] = 0;
  }
  return val;
}


static Boolean
GetTSAFieldsFromString
(CharPtr str,
 CharPtr PNTR p_submitter_reference,
 CharPtr PNTR p_archive_id,
 CharPtr PNTR p_description,
 CharPtr PNTR p_assembly,
 Int4Ptr p_taxon_id)
{
  CharPtr cp, cp_next, cp_equal, cp_end;
  CharPtr subref = NULL, arch_id = NULL, desc = NULL, assembly = NULL, tmp;
  Boolean is_bad = FALSE;

  if (p_submitter_reference != NULL) {
    *p_submitter_reference = NULL;
  }
  if (p_archive_id != NULL) {
    *p_archive_id = NULL;
  }
  if (p_submitter_reference != NULL) {
    *p_description = NULL;
  }
  if (StringHasNoText (str)) {
    return TRUE;
  }

  cp = StringChr (str, '[');
  while (cp != NULL && !is_bad) {
    cp++;
    cp_next = StringChr (cp + 1, '[');
    cp_equal = StringChr (cp, '=');
    cp_end = StringChr (cp, ']');
    if (cp_equal == NULL || cp_end == NULL) {
      is_bad = TRUE;
    } else if (cp_equal > cp_end) {
      is_bad = TRUE;
    } else if (cp_next != NULL && (cp_equal > cp_next || cp_end > cp_next)) {
      is_bad = TRUE;
    } else {
      cp += StringSpn (cp, " \t");
      if (BracketMatchesLabel (cp, cp_equal, "subref")) {
        if (subref == NULL) {
          subref = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "archive_id")) {
        if (arch_id == NULL) {
          arch_id = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "desc")) {
        if (desc == NULL) {
          desc = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "assembly")) {
        if (assembly == NULL) {
          assembly = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "taxon_id")) {
        tmp = GetBracketValue (cp_equal + 1, cp_end);
        if (p_taxon_id != NULL) {
          *p_taxon_id = atoi (tmp);
        }
      } else {
        is_bad = TRUE;
      }
    }
    cp = cp_next;
  }
  if (p_submitter_reference == NULL) {
    subref = MemFree (subref);
  } else {
    *p_submitter_reference = subref;
  }
  if (p_archive_id == NULL) {
    arch_id = MemFree (arch_id);
  } else {
    *p_archive_id = arch_id;
  }
  if (p_description == NULL) {
    desc = MemFree (desc);
  } else {
    *p_description = desc;
  }
  if (p_assembly == NULL) {
    assembly = MemFree (assembly);
  } else {
    *p_assembly = assembly;
  }
  return TRUE;
}


static void PrintTraceGapsXML (TGapInfoPtr gap_info)
{
  Int4 i;

  if (gap_info != NULL) {
    printf ("    <ntracegaps>%d</ntracegaps>\n", gap_info->num_gaps);
    if (gap_info->num_gaps > 0) {
      printf ("      <tracegaps source=\"INLINE\">");
      for (i = 0; i < gap_info->num_gaps - 1; i++) {
        printf ("%d,", gap_info->gap_offsets[i]);
      }
      printf ("%d</tracegaps>\n", gap_info->gap_offsets[gap_info->num_gaps - 1]);
    }
  }
}


static void TestPosConversions (TGapInfoPtr gap_info)
{
  Int4 i, t_pos, s_pos = 0, r_pos;
  Int4 test_len = 0;

  if (gap_info != NULL && gap_info->num_gaps > 0) {
    for (i = 0; i < gap_info->num_gaps; i++) {
      test_len += gap_info->gap_offsets[i] + 1;
    }
    for (i = 0; i < test_len; i++) {
      s_pos = SeqPosFromTilingPos (i, gap_info);
      t_pos = TilingPosFromSeqPos (s_pos, gap_info);
      if (t_pos != i) {
        printf ("Failed!  %d -> SeqPosFromTilingPos -> %d -> TilingPosFromSeqPos -> %d\n",
                i, s_pos, t_pos);
      }
      r_pos = SeqPosFromTilingPos (t_pos, gap_info);
      if (r_pos != s_pos) {
        printf ("Failed!  %d -> TilingPosFromSeqPos -> %d -> SeqPosFromTilingPos -> %d\n",
                s_pos, t_pos, r_pos);
      }
      /* printf ("%d:%d:%d:%d\n", i, s_pos, t_pos, r_pos); */
    }
  }
}


static void PrintTraceReadXML (TContigReadPtr read)
{
  if (read == NULL) {
    printf ("Bad read\n");
  } else {
    printf ("<trace>\n");
    printf ("  <trace_name>%s</trace_name>\n", read->read_id == NULL ? "" : read->read_id);
    PrintTraceGapsXML (read->gaps);
    printf ("  <nbasecalls>%d</nbasecalls>\n", StringLen (read->read_seq));
    printf ("  <valid>\n");
    printf ("    <start>%d</start>\n", read->read_assem_start + 1);
    printf ("    <stop>%d</stop>\n", read->read_assem_stop + 1);
    printf ("  </valid>\n");
    printf ("  <tiling direction = \"%s\">\n", read->is_complement ? "REVERSE" : "FORWARD");
    printf ("    <start>%d</start>\n", read->cons_start + 1);
    printf ("    <start>%d</start>\n", read->cons_start + StringLen (read->read_seq) + 1);
    printf ("  </tiling>\n");
    printf ("  <consensus>\n");
    printf ("    <start>%d</start>\n", read->cons_start + 1);
    printf ("    <start>%d</start>\n", read->cons_start + StringLen (read->read_seq) + 1);
    printf ("  </consensus>\n");
    printf ("<trace>\n");
  }
}



static void TestGapInfoReading (CharPtr gap_string)
{
  TGapInfoPtr  gap_info;
  ValNodePtr   list, vnp;
  
  if (!StringHasNoText (gap_string)) {
    gap_info = GapInfoFromSequenceString(gap_string, "*");
    if (gap_info == NULL) {
      printf ("error reading");
    } else {
      PrintTraceGapsXML (gap_info);
      TestPosConversions (gap_info);
      list = GetTransitionsFromGapInfo (gap_info, 0, 0, 40);
      for (vnp = list; vnp != NULL; vnp = vnp->next) {
        printf ("%d\n", vnp->data.intvalue);
      }
    }
    GapInfoFree (gap_info);
  }
}


static void AddAlignmentToSeqEntry (DenseSegPtr dsp, SeqEntryPtr sep)
{
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;

  if (dsp == NULL || sep == NULL) return;

  sap = SeqAnnotNew ();
  sap->type = 2;

  salp = SeqAlignNew ();
  salp->type = 3;
  salp->segtype = 2;
  salp->segs = (Pointer) dsp;
  salp->dim = dsp->dim;
  sap->data = (Pointer) salp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap->next = bsp->annot;
    bsp->annot = sap;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap->next = bssp->annot;
    bssp->annot = sap;
  }
}


static void AddDescrToNucBioseqCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp, sdp_copy;

  if (bsp == NULL || !ISA_na (bsp->mol) || data == NULL) { 
    return;
  }
  sdp = (SeqDescrPtr) data;
  while (sdp != NULL) {
    if (sdp->choice != Seq_descr_pub) {
      sdp_copy = (SeqDescrPtr) AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescrAsnRead, (AsnWriteFunc) SeqDescrAsnWrite);
      sdp_copy->next = bsp->descr;
      bsp->descr = sdp_copy;
    }
    sdp = sdp->next;
  }
}

  
static SeqSubmitPtr AddSeqSubmitFromTemplate (SeqEntryPtr sep, CharPtr filename)
{
  SeqSubmitPtr   ssp = NULL;
  SubmitBlockPtr sbp;
  CitSubPtr      csp;
  FILE *fp = NULL;
  Pointer         dataptr;
  Uint2           datatype;

  if (StringHasNoText (filename)) {
    return NULL;
  }
    
  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    printf ("Unable to read template file %s\n", filename);
    return NULL;
  }

  while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
    if (datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) dataptr;
      ssp->datatype = 1;
      ssp->data = sep;
    } else if (datatype == OBJ_SUBMIT_BLOCK) {
      sbp = (SubmitBlockPtr) dataptr;
      ssp = SeqSubmitNew ();
      ssp->datatype = 1;
      ssp->data = sep;
      ssp->sub = sbp;
    } else if (datatype == OBJ_SEQDESC) {
      VisitBioseqsInSep (sep, dataptr, AddDescrToNucBioseqCallback);
      ObjMgrFree (datatype, dataptr);
    } else {
      ObjMgrFree (datatype, dataptr);
    }
  }
  FileClose (fp);
  if (ssp == NULL) {
    ssp = SeqSubmitNew ();
    ssp->datatype = 1;
    ssp->data = sep;
  }

  if (ssp->sub == NULL) {
    ssp->sub = SubmitBlockNew ();
  } 

  ssp->sub->tool = MemFree (ssp->sub->tool);
  ssp->sub->tool = StringSave ("aceread");
  ssp->sub->hup = FALSE;
  ssp->sub->reldate = DateFree (ssp->sub->reldate);
  csp = ssp->sub->cit;
  if (csp != NULL) {
    csp->date = DateFree (csp->date);
    csp->date = DateCurr ();
  }
  return ssp;
}


static Boolean AddReadQualityScores (TACEFilePtr afp, CharPtr qs_filename, CharPtr rd_filename)
{
  ReadBufferData q, r;
  Boolean use_fasta = FALSE;
  Boolean rval = FALSE;

  if (afp == NULL || StringHasNoText (qs_filename)) {
    return TRUE;
  }

  q.current_data = NULL;
  r.current_data = NULL;

  q.fp = FileOpen (qs_filename, "r");
  if (q.fp == NULL) {
    printf ("Unable to read quality score file\n");
    return FALSE;
  }

  if (!StringHasNoText (rd_filename)) {
    r.fp = FileOpen (rd_filename, "r");
    if (r.fp == NULL) {
      printf ("Unable to open read FASTA file\n");
      FileClose (q.fp);
      return FALSE;
    }
    use_fasta = TRUE;
  }

  if (AddReadQualScores (afp, AbstractReadFunction, &q, use_fasta ? AbstractReadFunction : NULL, &r) > 0) {
    rval = TRUE;
  }

  FileClose (q.fp);
  if (use_fasta) {
    FileClose (r.fp);
  }
  return rval;
}


static Boolean LIBCALL MyBioseqSetAsnWrite (BioseqSetPtr bsp, AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	Boolean retval = FALSE;

	if (aip == NULL)
		return FALSE;

	atp = AsnLinkType(orig, AsnFind ("Bioseq-set"));   /* link local tree */
	if (atp == NULL) return FALSE;

	if (bsp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

	if (! AsnOpenStruct(aip, atp, (Pointer)bsp)) goto erret;
    
  if (bsp->id != NULL)
	{
    if (! ObjectIdAsnWrite(bsp->id, aip, AsnFind ("Bioseq-set.id"))) goto erret;
	}
  if (bsp->coll != NULL)
	{
    if (! DbtagAsnWrite(bsp->coll, aip, AsnFind ("Bioseq-set.coll"))) goto erret;
	}
  if (bsp->level != INT2_MIN)
  {
    av.intvalue = bsp->level;
    if (! AsnWrite(aip, AsnFind ("Bioseq-set.level"), &av)) goto erret;
  }
  if (bsp->_class != 0)
  {
    av.intvalue = bsp->_class;
    if (! AsnWrite(aip, AsnFind ("Bioseq-set.class"), &av)) goto erret;
  }
  if (bsp->release != NULL)
  {
    av.ptrvalue = bsp->release;
    if (! AsnWrite(aip, AsnFind ("Bioseq-set.release"), &av)) goto erret;
  }
  if (bsp->date != NULL)
	{
    if (! DateAsnWrite(bsp->date, aip, AsnFind ("Bioseq-set.date"))) goto erret;
	}
  if (bsp->descr != NULL)              /* Seq-descr optional */
	{
    if (! SeqDescrAsnWrite(bsp->descr, aip, AsnFind ("Bioseq-set.descr"))) goto erret;
	}

  if (! AsnOpenStruct(aip, AsnFind ("Bioseq-set.seq-set"), (Pointer)bsp->seq_set)) goto erret;
  /* this is where we stop */
  retval = TRUE;
erret:
	AsnUnlinkType(orig);        /* unlink local tree */
	return retval;
}

static Boolean LIBCALL MySeqEntryAsnWrite (SeqEntryPtr sep, AsnIoPtr aip, AsnTypePtr orig)
{
  AsnTypePtr atp;
	DataVal av;
	Boolean retval = FALSE;

	if (aip == NULL)
		return FALSE;

	atp = AsnLinkType(orig, AsnFind ("Seq-entry"));   /* link local tree */
	if (atp == NULL) return FALSE;

	if (sep == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

	av.ptrvalue = (Pointer)sep;
  if (! AsnWriteChoice(aip, atp, (Int2)sep->choice, &av)) goto erret;
  if (sep->choice == 1)
	{
    if (! BioseqAsnWrite((BioseqPtr)sep->data.ptrvalue, aip, AsnFind ("Seq-entry.seq"))) 
    {
			goto erret;
    }
	}
  else if (sep->choice == 2)
	{
    if (! MyBioseqSetAsnWrite((BioseqSetPtr)sep->data.ptrvalue, aip, AsnFind ("Seq-entry.set")))
    {
			goto erret;
    }
	}
  /* this is where we stop */
	retval = TRUE;
erret:
  AsnUnlinkType(orig);
  return retval;
}


static Boolean MySeqSubmitAsnWrite (AsnIoPtr aip, SubmitBlockPtr sbp, SeqDescrPtr desc_list)
{
	DataVal av;
	AsnTypePtr atp;
  Boolean retval = FALSE;
	SeqEntryPtr sep;
  SeqSubmitPtr ssp = NULL;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp, sdp_copy;

	if (aip == NULL)
		return FALSE;

	atp = AsnLinkType(NULL, AsnFind ("Seq-submit"));   /* link local tree */
  if (atp == NULL)
    return FALSE;

  ssp = SeqSubmitNew ();
  ssp->sub = sbp;
  ssp->datatype = 1;
  sep = SeqEntryNew ();
  sep->choice = 2;
  bssp = BioseqSetNew ();
  bssp->_class = BioseqseqSet_class_genbank;

  if (desc_list != NULL) {
    for (sdp = desc_list; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice == Seq_descr_pub) {
        sdp_copy = (SeqDescrPtr) AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescrAsnRead, (AsnWriteFunc) SeqDescrAsnWrite);
        sdp_copy->next = bssp->descr;
        bssp->descr = sdp_copy;
      }
    }
  }

  sep->data.ptrvalue = bssp;
  ssp->data = sep;

  if (! AsnOpenStruct(aip, atp, (Pointer)ssp))
    goto erret;

	if (! SubmitBlockAsnWrite(ssp->sub, aip, AsnFind ("Seq-submit.sub"))) goto erret;

	av.ptrvalue = ssp->data;
  if (! AsnWriteChoice(aip, AsnFind ("Seq-submit.data"), (Int2)ssp->datatype, &av)) goto erret;

	if (! AsnOpenStruct(aip, AsnFind ("Seq-submit.data.entrys"), ssp->data)) goto erret;
	sep = (SeqEntryPtr) ssp->data;
	if (! MySeqEntryAsnWrite(sep, aip, AsnFind ("Seq-submit.data.entrys.E"))) goto erret;
  /* This is where we stop */
  retval = TRUE;
erret:
  ssp->sub = NULL;
  ssp = SeqSubmitFree (ssp);
	return retval;
}

static void StartSeqSubmit (AsnIoPtr aip, SubmitBlockPtr sbp, SeqDescrPtr desc_list)
{

  if (aip == NULL || aip->fp == NULL) {
    return;
  }

  if (sbp == NULL) {
    fprintf (aip->fp, "Seq-entry ::= set {\n");
    fprintf (aip->fp, "class genbank ,\n");
    fprintf (aip->fp, "seq-set {\n");
  } else {
    MySeqSubmitAsnWrite (aip, sbp, desc_list);
    AsnIoFlush (aip);
  }
}


static DenseSegPtr DenseSegFromConsensusReadAln (TConsensusReadAlnPtr aln, CharPtr contig_id, CharPtr read_id) 
{
  DenseSegPtr dsp;
  Int4        i;

  if (aln == NULL) {
    return NULL;
  }

  dsp = DenseSegNew ();
  dsp->dim = 2;
  dsp->numseg = aln->numseg;
  dsp->ids = MakeSeqID (contig_id);
  dsp->ids->next = MakeSeqID (read_id);
  dsp->starts = (Int4Ptr) MemNew (sizeof (Int4) * dsp->dim * dsp->numseg);
  dsp->lens = (Int4Ptr) MemNew (sizeof (Int4) * dsp->numseg);
  if (aln->is_complement) {
    dsp->strands = (Uint1Ptr) MemNew (sizeof (Uint1) * dsp->dim * dsp->numseg);
    for (i = 0; i < dsp->numseg; i++) {
      dsp->strands[i * 2] = Seq_strand_plus;
      dsp->strands[(i * 2) + 1] = Seq_strand_minus;
    }
  }
  for (i = 0; i < dsp->numseg; i++) {
    dsp->starts[i * 2] = aln->cons_starts[i];
    dsp->starts[(i * 2) + 1] = aln->read_starts[i];
    dsp->lens [i] = aln->lens[i];
  }
  return dsp;
}


static SeqAlignPtr SeqAlignsForConsensusAndReads (TContigPtr contig)
{
  SeqAlignPtr salp_list = NULL, salp_last = NULL, salp_tmp;
  TConsensusReadAlnPtr aln;
  DenseSegPtr dsp;
  Int4 i;

  if (contig == NULL) {
    return NULL;
  }

  for (i = 0; i < contig->num_reads; i++) {
    aln = GetConsensusReadAln (contig->consensus_seq, contig->reads[i]);
    if (aln != NULL) {
      dsp = DenseSegFromConsensusReadAln (aln, contig->consensus_id, contig->reads[i]->read_id);
      if (dsp != NULL) {
        salp_tmp = SeqAlignNew ();
        salp_tmp->type = SAT_MASTERSLAVE;
        salp_tmp->segtype = SAS_DENSEG;
        salp_tmp->segs = dsp;
        salp_tmp->dim = 2;
        if (salp_list == NULL) {
          salp_list = salp_tmp;
        } else {
          salp_last->next = salp_tmp;
        }
        salp_last = salp_tmp;
      }
    }
  }
  return salp_list;
}


static SeqEntryPtr MakeContigSeqEntryWithReads (TContigPtr contig)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  sep, sep_prev;
  Int4 i;
  SeqAlignPtr salp;

  if (contig == NULL) {
    return NULL;
  }

  bssp = BioseqSetNew ();
  bssp->_class = BioseqseqSet_class_genbank;
  bssp->seq_set = MakeSeqEntryFromContig (contig);
  salp = SeqAlignsForConsensusAndReads (contig);
  if (salp != NULL) {
    bssp->annot = SeqAnnotNew ();
    bssp->annot->type = 2;
    bssp->annot->data = salp;
  }
  sep_prev = bssp->seq_set;
  for (i = 0; i < contig->num_reads; i++) {
    sep = MakeSeqEntryFromRead (contig->reads[i]);
    sep_prev->next = sep;
    sep_prev = sep;
  }
  sep = SeqEntryNew ();
  sep->choice = 2;
  sep->data.ptrvalue = bssp;
  return sep;
}


static void WriteXMLMsgUnableToOpenFile (CharPtr has_errors, CharPtr filename)
{
  if (has_errors == NULL || filename == NULL) {
    return;
  }
  if (*has_errors == 0) {
    printf ("<aceread>\n");
    *has_errors = 1;
  }
  printf ("<message severity=\"ERROR\" seq-id=\"No ID\" code=\"bad_format\">Unable to open %s</message>\n", filename);
}


typedef struct contigcountcallback {
  Int4 num_contigs;
  Uint4 num_conbases;
  Int4 num_reads;
  Uint4 num_readbases;
  Int4  file_num;
} ContigCountCallbackData, PNTR ContigCountCallbackPtr;


static ContigCountCallbackPtr ContigCountCallbackNew ()
{
  ContigCountCallbackPtr c;

  c = (ContigCountCallbackPtr) MemNew (sizeof (ContigCountCallbackData));
  MemSet (c, 0, sizeof (ContigCountCallbackData));
  return c;
}


static ContigCountCallbackPtr SummarizeContigCountList (ValNodePtr list)
{
  ContigCountCallbackPtr summ, c;
  
  summ = ContigCountCallbackNew();
  while (list != NULL) {
    c = (ContigCountCallbackPtr) list->data.ptrvalue;
    if (c != NULL) {
      summ->num_contigs += c->num_contigs;
      summ->num_conbases += c->num_conbases;
      summ->num_reads += c->num_reads;
      summ->num_readbases += c->num_readbases;
    }
    list = list->next;
  }
  return summ;
}


typedef struct contigfilelist {
  ValNodePtr list;
  Int4 max_bases;
  ContigCountCallbackPtr current;
} ContigFileListData, PNTR ContigFileListPtr;


#define ONE_CONTIG_FOR_FIRST

static char ProcessContigCountCallback (TContigPtr contig, void *data)
{
  ContigFileListPtr list;
  Int4 i;

  list = (ContigFileListPtr) data;
  if (contig == NULL || list == NULL) {
    return 0;
  }

  if (list->current == NULL || list->current->num_conbases > list->max_bases
#ifdef ONE_CONTIG_FOR_FIRST
      || list->list->next == NULL
#endif
    ) {
    list->current = ContigCountCallbackNew();
    list->current->file_num = ValNodeLen (list->list);
    ValNodeAddPointer (&(list->list), 0, list->current);
  }

  list->current->num_contigs++;
  list->current->num_conbases += contig->consensus_seq_len;
  list->current->num_reads += contig->num_reads;

  for (i = 0; i < contig->num_reads; i++) {
    list->current->num_readbases += contig->reads[i]->read_len;
  }
  return 1;
}


typedef enum {
  eReadNameType_local = 0,
  eReadNameType_TI,
  eReadNameType_SRR } EReadNameType;

static EReadNameType ReadNameTypeFromArg (CharPtr arg)
{
  EReadNameType read_name_type = eReadNameType_local;

  if (arg != NULL) {
    if (StringNICmp (arg, "T", 1) == 0) {
      read_name_type = eReadNameType_TI;
    } else if (StringNICmp (arg, "S", 1) == 0) {
      read_name_type = eReadNameType_SRR;
    }
  }
  return read_name_type;
}


typedef struct contigcallback {
  AsnIoPtr asn1_out;
  AsnTypePtr atp;
  FILE *fasta_out;
  FILE *qual_out;
  FILE *xml_out;

  ValNodePtr file_counts_list;
  Int4 contig_count;

  CharPtr fasta_base;
  CharPtr asn_base;
  CharPtr xml_base;
  CharPtr qual_base;

  /* XML values */
  CharPtr subref;
  CharPtr center_name;
  Int4    taxid;
  CharPtr description;
  CharPtr assembly;

  Boolean recalculate_consensus;
  Boolean recalculate_only_Ns;

  Boolean no_lookup;
  Boolean is_srr;
  Boolean asn1_include_reads;

  EReadNameType read_name_type;

  SeqIdReplaceListPtr id_replacement_list;

  SubmitBlockPtr sbp;
  SeqDescrPtr desc_list;

  char *has_errors;
} ContigCallbackData, PNTR ContigCallbackPtr;


static AsnIoPtr StartAsnFile (CharPtr filename, SubmitBlockPtr sbp, SeqDescrPtr desc_list)
{
  AsnIoPtr aip;

  aip = AsnIoOpen (filename, "w");
  if (aip != NULL) {
    aip->indent_level = 1;
    aip->first[aip->indent_level] = FALSE;
    StartSeqSubmit (aip, sbp, desc_list);
  }
  return aip;
}


static AsnIoPtr EndAsnFile (AsnIoPtr aip, Boolean is_submitblock)
{
  if (aip != NULL) {
    AsnIoFlush (aip);
    if (is_submitblock) {
      fprintf (aip->fp, " } } } }\n");
    } else {
      fprintf (aip->fp, " } }\n");
    }
    AsnIoClose (aip);
    aip = NULL;
  }
  return aip;
}


static char ProcessContigCallback (TContigPtr contig, void *data)
{
  ContigCallbackPtr c;
  SeqEntryPtr       sep;
  Char              filename[300];
  ContigCountCallbackPtr count = NULL;
  ValNodePtr             tmp;
  Boolean write_out = FALSE;
  Int4    i, ti;
  char rval = 0;

  c = (ContigCallbackPtr) data;
  if (contig == NULL || c == NULL) {
    return 0;
  }

  if (c->id_replacement_list != NULL) {
    UpdateContigIds (contig, c->id_replacement_list, c->no_lookup, c->is_srr, c->has_errors);
  }

  if (c->read_name_type == eReadNameType_TI) {
    for (i = 0; i < contig->num_reads; i++) {
      if (contig->reads[i]->read_id != NULL) {
        ti = atoi (contig->reads[i]->read_id);
        if (ti < 1) {
          if (*(c->has_errors) == 0) {
            printf ("<aceread>\n");
            *(c->has_errors) = 1;
          }
          printf ("<message severity=\"ERROR\" seq-id=\"%s\" code=\"bad_format\">Non-integer value for ti</message>\n", contig->reads[i]->read_id);
        } else if (contig->reads[i]->ti == 0) {
          contig->reads[i]->ti = ti;
          free (contig->reads[i]->read_id);
          contig->reads[i]->read_id = NULL;
        } else if (ti == contig->reads[i]->ti) {
          free (contig->reads[i]->read_id);
          contig->reads[i]->read_id = NULL;
        } else {
          if (*(c->has_errors) == 0) {
            printf ("<aceread>\n");
            *(c->has_errors) = 1;
          }
          printf ("<message severity=\"ERROR\" seq-id=\"%s\" code=\"bad_format\">Conflicting values for ti</message>\n", contig->reads[i]->read_id);
        }
      }
    }
  } else if (c->read_name_type == eReadNameType_SRR) {
    for (i = 0; i < contig->num_reads; i++) {
      if (contig->reads[i]->read_id != NULL) {
        if (contig->reads[i]->srr == NULL) {
          contig->reads[i]->srr = contig->reads[i]->read_id;
          contig->reads[i]->read_id = NULL;
        } else if (StringCmp (contig->reads[i]->read_id, contig->reads[i]->srr) == 0) {
          free (contig->reads[i]->read_id);
          contig->reads[i]->read_id = NULL;
        } else {
          if (*(c->has_errors) == 0) {
            printf ("<aceread>\n");
            *(c->has_errors) = 1;
          }
          printf ("<message severity=\"ERROR\" seq-id=\"%s\" code=\"bad_format\">Conflicting values for srr</message>\n", contig->reads[i]->read_id);
        }
      }
    }
  }

  if (c->recalculate_consensus) {
    /* TODO - add read quality scores ? */

    if (ReplaceConsensusSequenceFromTraces (contig, c->recalculate_only_Ns) > 0) {
      write_out = TRUE;
    }
  } else {
    write_out = TRUE;
  }

  c->contig_count ++;

  if (write_out) {
    rval = 1;
    if (c->file_counts_list != NULL) {
      count = c->file_counts_list->data.ptrvalue;
    }

    /* write ASN.1 */
    if (c->asn1_out == NULL 
        && c->asn_base != NULL && count != NULL) {
      sprintf (filename, "%s.%d", c->asn_base, count->file_num);
      c->asn1_out = StartAsnFile (filename, c->sbp, c->desc_list);
      if (c->asn1_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c->has_errors, filename);
        rval = 0;
      }
    }
    if (c->asn1_out != NULL) {
      if (c->asn1_include_reads) {
        sep = MakeContigSeqEntryWithReads (contig);
      } else {
        sep = MakeSeqEntryFromContig (contig);
      }
      if (c->desc_list != NULL) {
        VisitBioseqsInSep (sep, c->desc_list, AddDescrToNucBioseqCallback);
      }
      SeqEntryAsnWrite(sep, c->asn1_out, c->atp);
      sep = SeqEntryFree (sep);
      if (count != NULL && c->contig_count >= count->num_contigs) {
        c->asn1_out = EndAsnFile (c->asn1_out, c->sbp != NULL);
      }
    }
    
    /* write FASTA */
    if (c->fasta_out == NULL
        && c->fasta_base != NULL && count != NULL) {
      sprintf (filename, "%s.%d", c->fasta_base, count->file_num);
      c->fasta_out = FileOpen (filename, "w");
      if (c->fasta_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c->has_errors, filename);
        rval = 0;
      }
    }
    if (c->fasta_out != NULL) {
      WriteFASTAFromContig (contig, c->fasta_out);
      if (count != NULL && c->contig_count >= count->num_contigs) {
        FileClose (c->fasta_out);
        c->fasta_out = NULL;
      }
    }

    /* write quality scores */
    if (c->qual_out == NULL
        && c->qual_base != NULL && count != NULL) {
      sprintf (filename, "%s.%d", c->qual_base, count->file_num);
      c->qual_out = FileOpen (filename, "w");
      if (c->qual_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c->has_errors, filename);
        rval = 0;
      }
    }
    if (c->qual_out != NULL) {
      WriteContigQualScores (contig, c->qual_out);
      if (count != NULL && c->contig_count >= count->num_contigs) {
        FileClose (c->qual_out);
        c->qual_out = NULL;
      }
    }

    /* write XML */
    if (c->xml_out == NULL
        && c->xml_base != NULL && count != NULL) {
      sprintf (filename, "%s.%d", c->xml_base, count->file_num);
      c->xml_out = FileOpen (filename, "w");
      WriteTraceAssemblyHeader ("UPDATE", c->subref, c->center_name, c->taxid, c->description, c->assembly,
                                count->num_contigs, count->num_conbases, count->num_reads, count->num_readbases,
                                c->xml_out);

      if (c->xml_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c->has_errors, filename);
        rval = 0;
      }
    }
    if (c->xml_out != NULL) {
      WriteTraceAssemblyFromContig (contig, c->xml_out);
      if (count != NULL && c->contig_count >= count->num_contigs) {
        WriteTraceAssemblyTrailer (c->xml_out);
        FileClose (c->xml_out);
        c->xml_out = NULL;
      }
    }
  }

  if (count != NULL && c->contig_count >= count->num_contigs) {
    tmp = c->file_counts_list;
    c->file_counts_list = tmp->next;
    tmp->next = NULL;
    tmp = ValNodeFreeData (tmp);
    c->contig_count = 0;
  }

  return 1;
}

static void ReadLargeAceFile 
(CharPtr acefile,
 CharPtr asn1_out,
 CharPtr fasta_out,
 CharPtr template_in,
 CharPtr qual_scores_out,
 CharPtr xml_out,
 CharPtr id_lookup,
 char *has_errors,
 Boolean recalculate_consensus,
 Boolean recalculate_only_Ns,
 CharPtr subref,
 CharPtr center_name,
 Int4    taxid,
 CharPtr description,
 CharPtr assembly,
 Boolean no_lookup,
 Boolean is_srr,
 Boolean make_qual_scores,
 Int4    chunk_size,
 EReadNameType read_name_type,
 Boolean include_reads)
{
  ReadBufferData    rbd;
  ContigCallbackData c;
  SeqEntryPtr old_scope;
  FILE *f;
  SeqSubmitPtr   ssp = NULL;
  CitSubPtr      csp;
  Pointer         dataptr;
  Uint2           datatype;
  SeqDescrPtr     sdp, sdp_next;
  ContigFileListData file_count_list;
  ContigCountCallbackPtr summ;

  MemSet (&c, 0, sizeof (ContigCallbackData));

  c.no_lookup = no_lookup;
  c.is_srr = is_srr;
  c.has_errors = has_errors;
  c.asn1_include_reads = include_reads;
  c.read_name_type = read_name_type;

  /* filenames */
  c.asn_base = asn1_out;
  c.asn1_out = NULL;
  c.fasta_base = fasta_out;
  c.fasta_out = NULL;
  c.qual_base = qual_scores_out;
  c.qual_out = NULL;
  c.xml_base = xml_out;
  c.xml_out = NULL;

  /* XML values */
  c.subref = subref;
  c.center_name = center_name;
  c.taxid = taxid;
  c.description = description;
  c.assembly = assembly;

  file_count_list.list = NULL;
  file_count_list.current = NULL;
  file_count_list.max_bases = chunk_size;

  rbd.fp = OpenAceFile (acefile);
  if (rbd.fp == NULL) {
    WriteXMLMsgUnableToOpenFile (c.has_errors, acefile);
    goto escape;
  }
  rbd.current_data = NULL;

  ProcessLargeACEFileForContigFastaAndQualScores ( AbstractReadFunction, &rbd, 
                                                          qual_scores_out == NULL ? make_qual_scores : TRUE,
                                                          has_errors, ProcessContigCountCallback, &file_count_list);

  FileClose (rbd.fp);
  rbd.fp = NULL;

  /* prepare XML output */
  if (c.xml_base != NULL) {
    if (chunk_size < 1) {
      summ = SummarizeContigCountList (file_count_list.list);
      c.xml_out = FileOpen (c.xml_base, "w");
      if (c.xml_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.xml_base);
        goto escape;
      }
      WriteTraceAssemblyHeader ("NEW", c.subref, c.center_name, c.taxid, c.description, c.assembly,
                                summ->num_contigs, summ->num_conbases, summ->num_reads, summ->num_readbases,
                                c.xml_out);
      summ = MemFree (summ);
      file_count_list.list = ValNodeFreeData (file_count_list.list);
    } else {
#ifdef ONE_CONTIG_FOR_FIRST
      /* temporarily, start the first file instead, which will have just one contig */
      c.xml_out = FileOpen (c.xml_base, "w");
      if (c.xml_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.xml_base);
        goto escape;
      }
      summ = (ContigCountCallbackPtr) file_count_list.list->data.ptrvalue;
      WriteTraceAssemblyHeader ("NEW", c.subref, c.center_name, c.taxid, c.description, c.assembly,
                                summ->num_contigs, summ->num_conbases, summ->num_reads, summ->num_readbases,
                                c.xml_out);
#else
      f = FileOpen (c.xml_base, "w");
      if (f == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.xml_base);
        goto escape;
      }
      WriteTraceAssemblyHeader ("NEW", c.subref, c.center_name, c.taxid, c.description, c.assembly,
                                0, 0, 0, 0,
                                f);
      WriteTraceAssemblyTrailer (f);
      FileClose (f);
#endif
    }
  } else {
    if (chunk_size < 1) {
      file_count_list.list = ValNodeFreeData (file_count_list.list);
    }
  }

  c.file_counts_list = file_count_list.list;

  /* read template file */
  c.sbp = NULL;
  c.desc_list = NULL;
  if (!StringHasNoText (template_in)) {
    f = FileOpen (template_in, "r");
    if (f == NULL) {
      WriteXMLMsgUnableToOpenFile (c.has_errors, template_in);
      goto escape;
    }
    while ((dataptr = ReadAsnFastaOrFlatFile (f, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
      if (datatype == OBJ_SEQSUB) {
        ssp = (SeqSubmitPtr) dataptr;
        c.sbp = ssp->sub;
        ssp->sub = NULL;
        ssp = SeqSubmitFree (ssp);
      } else if (datatype == OBJ_SUBMIT_BLOCK) {
        c.sbp = (SubmitBlockPtr) dataptr;
      } else if (datatype == OBJ_SEQDESC) {
        ValNodeLink (&(c.desc_list), (ValNodePtr) dataptr);
      } else {
        ObjMgrFree (datatype, dataptr);
      }
    }
    FileClose (f);
    if (c.sbp != NULL) {
      c.sbp->tool = MemFree (c.sbp->tool);
      c.sbp->tool = StringSave ("aceread");
      c.sbp->hup = FALSE;
      c.sbp->reldate = DateFree (c.sbp->reldate);
      csp = c.sbp->cit;
      if (csp != NULL) {
        csp->date = DateFree (csp->date);
        csp->date = DateCurr ();
      }
    }
  }

  c.atp = AsnFind ("Bioseq-set.seq-set.E");

  c.recalculate_consensus = recalculate_consensus;
  c.recalculate_only_Ns = recalculate_only_Ns;

  if (id_lookup != NULL) {
    f = FileOpen (id_lookup, "r");
    if (f == NULL) {
      WriteXMLMsgUnableToOpenFile (c.has_errors, id_lookup);
      goto escape;
    }
    c.id_replacement_list = ReadSeqIdPairListFromFile (f);
    SeqEntrySetScope (old_scope);
    FileClose (f);
  }

  if (chunk_size < 1) {
    if (c.asn_base != NULL) {
      c.asn1_out = StartAsnFile (c.asn_base, c.sbp, c.desc_list);
      if (c.asn1_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.asn_base);
        goto escape;
      }
    }
    if (c.fasta_base != NULL) {
      c.fasta_out = FileOpen (c.fasta_base, "w");
      if (c.fasta_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.fasta_base);
        goto escape;
      }
    }
    if (c.qual_out != NULL) {
      c.qual_out = FileOpen (c.qual_base, "w");
      if (c.qual_out == NULL) {
        WriteXMLMsgUnableToOpenFile (c.has_errors, c.qual_base);
        goto escape;
      }
    }
  }
  
  rbd.fp = OpenAceFile (acefile);
  if (rbd.fp == NULL) {
    WriteXMLMsgUnableToOpenFile (c.has_errors, acefile);
    goto escape;
  }
  rbd.current_data = NULL;

  ProcessLargeACEFileForContigFastaAndQualScores ( AbstractReadFunction, &rbd, 
                                                          qual_scores_out == NULL ? FALSE : TRUE,
                                                          has_errors, ProcessContigCallback, &c);


escape:
  FileClose (rbd.fp);
  c.id_replacement_list = SeqIdReplaceListFree (c.id_replacement_list);
  /* free c.desc_list */
  for (sdp = c.desc_list; sdp != NULL; sdp = sdp_next) {
    sdp_next = sdp->next;
    sdp->next = NULL;
    sdp = SeqDescrFree (sdp);
  }
  if (c.xml_out != NULL) {
    WriteTraceAssemblyTrailer (c.xml_out);
    FileClose (c.xml_out);
    c.xml_out = NULL;
  }
  if (c.asn1_out != NULL) {
    c.asn1_out = EndAsnFile (c.asn1_out, c.sbp != NULL);
  }
  if (c.fasta_out != NULL) {
    FileClose (c.fasta_out);
    c.fasta_out = NULL;
  }
  if (c.qual_out != NULL) {
    FileClose (c.qual_out);
    c.qual_out = NULL;
  }
  c.sbp = SubmitBlockFree (c.sbp);
}


Int2 Main (void)

{
  CharPtr      infile, outfile, xmlfile;

  ReadBufferData    rbd;
  TACEFilePtr afp;
  Int4        i, len;
  SeqEntryPtr sep;
  AsnIoPtr    aip;
  FILE *f = NULL;
  FILE *f2;
  CharPtr app = "aceread_tst";
  BioseqSetPtr bssp;
  SeqEntryPtr  last_sep = NULL;
  Uint2        entityID;
  Boolean      make_qual_scores, suppress_lookup, srr_ids, fasta_out;
  CharPtr      submitter_ref = NULL, archive_id = NULL, description = NULL, assembly = NULL;
  CharPtr      center_name = NULL;
  CharPtr      format = NULL;
  CharPtr      gap_string;
  CharPtr      asn_file = NULL;
  Int4         limit = 0;
  char         has_errors = 0;
  Boolean      recalculate_consensus = FALSE, recalculate_only_Ns = FALSE;
  CharPtr      recalculate_options;
  SeqSubmitPtr ssp;
  CharPtr      id_substitution_file = NULL;
  Int4         taxon_id = 0;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  ErrSetOpts (ERR_IGNORE, ERR_LOG_ON);

  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  PubSeqFetchEnable ();

  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  recalculate_options = (CharPtr) myargs[N_argRecalculateConsensus].strvalue;
  if (!StringHasNoText (recalculate_options)) {
    if (StringCmp (recalculate_options, "W") == 0) {
      recalculate_consensus = TRUE;
      recalculate_only_Ns = FALSE;
    } else if (StringCmp (recalculate_options, "N") == 0) {
      recalculate_consensus = TRUE;
      recalculate_only_Ns = TRUE;
    } else {
      Message (MSG_FATAL, "Invalid consensus sequence recalculation option");
      return 1;
    }
  }


  /* test gap info reading if provided */
  gap_string = (CharPtr) myargs[G_argGapString].strvalue;
  TestGapInfoReading (gap_string);

  /* limit number of contigs?  for debugging purposes */
  limit = myargs[l_argLimitNumContigs].intvalue;

  /* select format of input file */
  format = (CharPtr) myargs[F_argFormat].strvalue;
  if (StringHasNoText (format)) {
    format = "A";
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  if (StringHasNoText (infile)) {
    Message (MSG_FATAL, "Must supply input file!");
    return 1;
  }
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  xmlfile = (CharPtr) myargs[X_argXMLFile].strvalue;
  make_qual_scores = (Boolean) myargs [Q_argMakeQualScores].intvalue;
  center_name = (CharPtr) myargs[C_argCenter].strvalue;
  suppress_lookup = (Boolean) myargs [L_argSuppressIdLookup].intvalue;
  srr_ids = (Boolean) myargs[R_argSRRids].intvalue;
  fasta_out = (Boolean) myargs[f_argFASTA].intvalue;

  /* ASN.1 file to validate against */
  asn_file = (CharPtr) myargs [V_argValidateAgainstAsn1File].strvalue;

  if (!GetTSAFieldsFromString ((CharPtr) myargs [T_argTSAFields].strvalue,
                               &submitter_ref,
                               &archive_id,
                               &description,
                               &assembly,
                               &taxon_id)) {
    Message (MSG_FATAL, "Error reading TSA fields");
    return 1;
  }

  if (!StringHasNoText (xmlfile) && (StringHasNoText (center_name) || taxon_id < 1)) {
    PrintACEFormatErrorXML ("Must specify center name and taxid for XML output", NULL, &has_errors);
    printf ("</aceread>\n");
    return 1;
  }        

  len = StringLen (infile);
  if (StringHasNoText (outfile)) {
    if (len > 3 && StringCmp (infile + len - 4, ".ace") == 0) {
      outfile = StringSave (infile);
      StringCpy (outfile + len - 3, "sqn");
    } else if (len > 6 && StringCmp (infile + len - 7, ".ace.gz") == 0) {
      outfile = StringSave (infile);
      StringCpy (outfile + len - 6, "sqn");
    } else {
      outfile = (CharPtr) MemNew (sizeof (Char) * (len + 5));
      sprintf (outfile, "%s.sqn", infile);
    }
  }

  if (!StringHasNoText ((CharPtr) myargs [S_argIDSubstitutionFile].strvalue)) {
    id_substitution_file = ((CharPtr) myargs [S_argIDSubstitutionFile].strvalue);
  }

  if (StringChr (format, 'A') != NULL) {
    ReadLargeAceFile (infile, fasta_out ? NULL : outfile,
                      fasta_out ? outfile : NULL, 
                      (CharPtr) myargs[t_argTemplateFile].strvalue,
                      NULL, xmlfile, id_substitution_file, &has_errors,
                      recalculate_consensus, recalculate_only_Ns,
                      submitter_ref, center_name, taxon_id, description, assembly, 
                      suppress_lookup, srr_ids, make_qual_scores,
                      myargs[c_argChunkSize].intvalue,
                      ReadNameTypeFromArg (myargs[n_argReadNameType].strvalue),
                      (Boolean) myargs [z_argIncludeReads].intvalue);
    if (has_errors) {
      printf ("</aceread>\n");
      return 1;
    } else {
      return 0;
    }
  }

  if (id_substitution_file != NULL) {
    f = FileOpen (id_substitution_file, "r");
    if (f == NULL) {
      Message (MSG_FATAL, "Unable to open %s", id_substitution_file);
      return 1;
    }
  }

  if (StringChr (format, 'M') != NULL) {
    rbd.fp = FileOpen (infile, "r");
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }

    rbd.current_data = NULL;
    afp = ReadMAQFile (AbstractReadFunction, &rbd);
  } else if (StringChr (format, 'E') != NULL) {
    rbd.fp = FileOpen (infile, "r");
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }

    rbd.current_data = NULL;
    afp = ReadElandStandaloneFile (AbstractReadFunction, &rbd);
  } else if (StringChr (format, 'A') != NULL) { 
    rbd.fp = OpenAceFile (infile);
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }
    rbd.current_data = NULL;
    afp = ReadACEFile ( AbstractReadFunction, &rbd, make_qual_scores, &has_errors);
  } else {
    Message (MSG_FATAL, "Unrecognized format: %s\n", format);
    return 1;
  }
  FileClose (rbd.fp);
  if (afp == NULL) {
    printf ("<message severity=\"ERROR\" seq-id=\"No ID\" code=\"bad_format\">Unable to read file</message>\n");
  } else {
    if (recalculate_consensus) {
        if (!AddReadQualityScores (afp, (CharPtr) myargs [q_argReadQualScoresFile].strvalue, (CharPtr) myargs [r_argReadFASTAFile].strvalue)) {
            printf ("<message severity=\"ERROR\" seq-id=\"No ID\" code=\"bad_format\">Failed to add read quality scores</message>\n");
        } else {
            RecalculateConsensusSequences (afp, recalculate_only_Ns);
        }
    }

    if (limit > 0) {
      for (i = limit; i < afp->num_contigs; i++) {
        ContigFree (afp->contigs[i]);
        afp->contigs[i] = NULL;
      }
      afp->num_contigs = limit;
    }

    if (f != NULL) {
      UpdateAceFileIds (afp, f, suppress_lookup, srr_ids, &has_errors);
      FileClose (f);
      f = NULL;
    }
    ValidateAceFileIds (afp, &has_errors);

    if (asn_file != NULL) {
      if (ValidateAgainstASNFile (afp, asn_file, &has_errors)) {
        printf ("Validation against %s succeeded\n", asn_file);
      }
    }

    if (!StringHasNoText (xmlfile)) {
      f2 = FileOpen (xmlfile, "w");
      WriteTraceAssemblyFromAceFile (afp, submitter_ref, center_name, 0, description, f2);
      FileClose (f2);
    }

    if (fasta_out) {
      f2 = FileOpen (outfile, "w");
      WriteFASTAFromAceFile (afp, f2);
      FileClose (f2);
    } else {
      aip = AsnIoOpen (outfile, "w");
      if (aip == NULL) {
        printf ("Unable to open %s\n", outfile);
      } else {
        bssp = BioseqSetNew ();
        bssp->_class = BioseqseqSet_class_genbank;

        for (i = 0; i < afp->num_contigs; i++) {
          sep = MakeSeqEntryFromContig (afp->contigs[i]);
          if (last_sep == NULL) {
            bssp->seq_set = sep;
          } else {
            last_sep->next = sep;
          }
          last_sep = sep;          
        }
        sep = ValNodeNew (NULL);
        sep->choice = 2;
        sep->data.ptrvalue = bssp;
        bssp->seqentry = sep;
        SeqMgrLinkSeqEntry (sep, 0, NULL);
        entityID = ObjMgrGetEntityIDForChoice (sep);
        AssignIDsInEntityEx (entityID, 0, NULL, NULL);
        SeqMgrIndexFeatures (entityID, sep);
        ssp = AddSeqSubmitFromTemplate (sep, (CharPtr) myargs[t_argTemplateFile].strvalue);
        if (ssp == NULL) {
          SeqEntryAsnWrite (sep, aip, NULL);
          sep = SeqEntryFree (sep);
        } else {
          SeqSubmitAsnWrite (ssp, aip, NULL);
          ssp = SeqSubmitFree (ssp);
        }
        AsnIoClose (aip);
      }
    }
  }

  if (has_errors) {
    printf ("</aceread>\n");
  }

  return 0;

}

