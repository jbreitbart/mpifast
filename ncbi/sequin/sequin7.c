/*   sequin7.c
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
* File Name:  sequin7.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/3/98
*
* $Revision: 6.352 $
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

#ifndef CODECENTER
static char *date_of_compilation = __DATE__;
static char *time_of_compilation = __TIME__;
#else
static char *date_of_compilation = "today";
static char *time_of_compilation = "now";
#endif

#include "sequin.h"
#include <gather.h>
#include <edutil.h>
#include <cdrgn.h>
#include <subutil.h>
#include <tofasta.h>
#include <vsm.h>
#include <document.h>
#include <maputil.h>
#include <asn2gnbp.h>
#include <bspview.h>
#include <findrepl.h>
#include <toasn3.h>
#include <toporg.h>
#include <utilpub.h>
#include <salsap.h>
#include <salptool.h>
#include <salutil.h>
#include <saledit.h>
#include <explore.h>
#include <seqpanel.h>
#include <alignmgr2.h>
#include <actutils.h>
#include <tax3api.h>
#include <algo/blast/api/blast_options_api.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/blast_api.h>
#include <salstruc.h>
#include <valid.h> /* added for latloncountry conflict checking */

#define CONVERT_TO_JOIN  1
#define CONVERT_TO_ORDER 2
#define DO_NOT_CONVERT   3

NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol     /* Returns special symbol if no SeqEntry */
 );

NLM_EXTERN SeqEntryPtr SequinFastaToSeqEntryExEx
  (
    FILE *fp,               /* file to get sequence from */ 
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr special_symbol, /* Returns special symbol if no SeqEntry */
    BoolPtr chars_stripped  /* set to TRUE if characters other than digits
                             * were stripped from the FASTA sequence data */
  )
{
  BioseqPtr    bsp;
  FileCache    fc;
  Boolean      forceNuc = FALSE;
  Boolean      forceProt = FALSE;
  Pointer      dataptr;
  Uint2        datatype;
  Char         line [128];
  Int4         pos;
  SeqEntryPtr  sep = NULL;
  CharPtr      str;

  if (errormsg != NULL) {
    *errormsg = NULL;
  }
  if (special_symbol != NULL) {
    *special_symbol = NULLB;
  }
  if (is_na) {
    forceNuc = TRUE;
  } else {
    forceProt = TRUE;
  }
  dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, forceNuc, forceProt, parseSeqId, FALSE, chars_stripped);
  if (dataptr != NULL) {
    if (datatype == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) dataptr;
      sep = SeqMgrGetSeqEntryForData (bsp);
      if (sep == NULL) {
        sep = SeqEntryNew ();
        if (sep != NULL) {
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        }
      }
    }
  } else if (special_symbol != NULL) {
    /* look ahead to see what character caused inability to interpret line */
    FileCacheSetup (&fc, fp);
    /* pos = FileCacheTell (&fc); */
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
    if (str != NULL && StringDoesHaveText (str)) {
      TrimSpacesAroundString (str);
    }
    *special_symbol = line [0];
    /* seek to start of next line after one that could not be interpreted */
    pos = FileCacheTell (&fc);
    FileCacheSetup (&fc, fp);
    FileCacheSeek (&fc, pos);
    fseek (fp, pos, SEEK_SET);
  }
  /* return FastaToSeqEntryInternal((void *)fp, 2, NULL,is_na, errormsg, parseSeqId, special_symbol); */
  return sep;
}

NLM_EXTERN SeqEntryPtr SequinFastaToSeqEntryEx 
  (
    FILE *fp,               /* file to get sequence from */ 
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
  )
{
  return SequinFastaToSeqEntryExEx (fp, is_na, errormsg, parseSeqId, special_symbol, NULL);    
}

static FonT  titleFont = NULL;

#ifndef WIN_MAC
void CreateSqnInitialFormMenus (WindoW w)

{
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    AddAboutAndHelpMenuItems (m);
    if (bfp->importform != NULL || bfp->exportform != NULL) {
      if (bfp->importform != NULL) {
        FormCommandItem (m, "Import...", bfp, VIB_MSG_IMPORT);
      }
      if (bfp->exportform != NULL) {
        FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
      }
      SeparatorItem (m);
    }
    FormCommandItem (m, "Quit", bfp, VIB_MSG_QUIT);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static void DefaultMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

typedef struct startupform {
  FORM_MESSAGE_BLOCK
} StartupForm, PNTR StartupFormPtr;

static void ChangeDestination (GrouP g)

{
  Char  str [64];
  Int2  val;

  val = GetValue (g);
  switch (val) {
    case 1 :
      RemoveAppProperty ("SequinUseEMBLStyle");
      RemoveAppProperty ("SequinUseDDBJStyle");
      if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
        if (! StringHasNoText (str)) {
          if (StringICmp (str, "GenBank") != 0) {
            WriteSequinAppParam ("PREFERENCES", "DATABASE", "GenBank");
          }
        }
      }
      break;
    case 2 :
      SetAppProperty ("SequinUseEMBLStyle", (void *) 1024);
      RemoveAppProperty ("SequinUseDDBJStyle");
      WriteSequinAppParam ("PREFERENCES", "DATABASE", "EMBL");
      break;
    case 3 :
      RemoveAppProperty ("SequinUseEMBLStyle");
      SetAppProperty ("SequinUseDDBJStyle", (void *) 1024);
      WriteSequinAppParam ("PREFERENCES", "DATABASE", "DDBJ");
      break;
    default :
      break;
  }
  SetupBioseqPageList ();
}

static void CenterString (RectPtr rptr, CharPtr text, FonT fnt, Int2 inc)

{
  if (fnt != NULL) {
    SelectFont (fnt);
  }
  rptr->bottom = rptr->top + LineHeight ();
  DrawString (rptr, text, 'c', FALSE);
  rptr->top = rptr->bottom + inc;
}

extern void DrawAbout (PaneL p)

{
  RecT  r;


  if (titleFont == NULL) {
#ifdef WIN_MAC
    titleFont = GetFont ("Geneva", 18, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MSWIN
    titleFont = GetFont ("Arial", 24, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 24, TRUE, TRUE, FALSE, "");
#endif
  }

  ObjectRect (p, &r);
  InsetRect (&r, 4, 4);
  r.top += 5;
  Blue ();
  CenterString (&r, "Sequin", titleFont, 5);
  CenterString (&r, SEQUIN_VERSION, programFont, 5);
  CenterString (&r, SEQUIN_SERVICES, programFont, 10);
  CenterString (&r, "National Center for Biotechnology Information", systemFont, 5);
  CenterString (&r, "National Library of Medicine", systemFont, 5);
  CenterString (&r, "National Institutes of Health", systemFont, 10);
  CenterString (&r, "(301) 496-2475", systemFont, 5);
  CenterString (&r, "info@ncbi.nlm.nih.gov", systemFont, 0);
}

extern Int2 AboutBoxWidth (void)

{
  Int2     max;
  CharPtr  ptr;
  Char     sequinServices [60];
  Char     sequinVersion [60];
  Int2     wid;


  if (titleFont == NULL) {
#ifdef WIN_MAC
    titleFont = GetFont ("Geneva", 18, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MSWIN
    titleFont = GetFont ("Arial", 24, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 24, TRUE, TRUE, FALSE, "");
#endif
  }

  sprintf (sequinVersion, "Sequin Application Version %s", SEQUIN_APPLICATION);
  ptr = "Standard Release";
/*#ifdef USE_ENTREZ*/
  if (useEntrez || useBlast) {
    ptr = "Network Aware";
  }
/*#endif*/
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    ptr = "Indexer Services";
  }
/*#endif*/
  if (genomeCenter != NULL) {
    ptr = "Genome Center";
  }
  sprintf (sequinServices, "%s [%s]", ptr, date_of_compilation);

  SelectFont (titleFont);
  max = StringWidth ("Sequin");
  SelectFont (programFont);
  wid = StringWidth (sequinVersion);
  if (wid > max) {
    max = wid;
  }
  wid = StringWidth (sequinServices);
  if (wid > max) {
    max = wid;
  }
  SelectFont (systemFont);
  wid = StringWidth ("National Center for Biotechnology Information");
  if (wid > max) {
    max = wid;
  }
  max += 2 * stdCharWidth + 2;
  return max;
}

extern Int2 AboutBoxHeight (void)

{
  Int2  hgt;

  if (titleFont == NULL) {
#ifdef WIN_MAC
    titleFont = GetFont ("Geneva", 18, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MSWIN
    titleFont = GetFont ("Arial", 24, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 24, TRUE, TRUE, FALSE, "");
#endif
  }

  SelectFont (titleFont);
  hgt = LineHeight () + 5;
  SelectFont (programFont);
  hgt += 2 * LineHeight () + 15;
  SelectFont (systemFont);
  hgt += 5 * LineHeight () + 25;
  hgt += 18;
  return hgt;
}

extern ForM CreateStartupForm (Int2 left, Int2 top, CharPtr title,
                               BtnActnProc startFa2htgs,
                               BtnActnProc startPhrap,
                               BtnActnProc buildContig,
                               BtnActnProc startNew,
                               BtnActnProc readExisting,
                               BtnActnProc fetchFromNet,
                               BtnActnProc showHelp,
                               BtnActnProc createSubmissionTemplate,
                               BtnActnProc quitProgram,
                               WndActnProc activateForm)

{
  ButtoN          b;
  GrouP           c;
  GrouP           d;
  GrouP           k;
  PaneL           p;
  StartupFormPtr  sfp;
  Char            str [32];
  WindoW          w;
#ifndef WIN_MAC
  MenU            m;
#endif

  w = NULL;
  sfp = MemNew (sizeof (StartupForm));
  if (sfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, sfp, StdCleanupFormProc);
    sfp->form = (ForM) w;
    sfp->formmessage = DefaultMessageProc;

#ifndef WIN_MAC
    m = PulldownMenu (w, "Misc");
    CommandItem (m, "Net Configure...", NetConfigureProc);
    if (useEntrez) {
      /*
      SeparatorItem (m);
      CommandItem (m, "Entrez Query...", EntrezQueryProc);
      SeparatorItem (m);
      CommandItem (m, "Entrez2 Query...", Entrez2QueryProc);
      */
      if (extraServices) {
        SeparatorItem (m);
        CommandItem (m, "Process FASTA Nucleotide Updates", ParseInNucUpdates);
      }
    }
    if (useDesktop) {
      SeparatorItem (m);
      VSMAddToMenu (m, VSM_DESKTOP);
    }
#endif

    p = SimplePanel (w, AboutBoxWidth (), AboutBoxHeight (), DrawAbout);

    k = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (k, 3, 10);
    StaticPrompt (k, "Database for submission", 0, stdLineHeight, programFont, 'l');
    d = HiddenGroup (k, 4, 0, ChangeDestination);
    RadioButton (d, "GenBank");
    RadioButton (d, "EMBL");
    RadioButton (d, "DDBJ");
    if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
      if (StringICmp (str, "GenBank") == 0) {
        SetValue (d, 1);
      } else if (StringICmp (str, "EMBL") == 0) {
        SetValue (d, 2);
      } else if (StringICmp (str, "DDBJ") == 0) {
        SetValue (d, 3);
      } else {
        SetValue (d, 1);
      }
    } else {
      SetValue (d, 1);
    }
    ChangeDestination (d);

    c = HiddenGroup (w, 1, 0, NULL);
    SetGroupSpacing (c, 10, 5);

    if (startFa2htgs != NULL) {
      b = PushButton (c, "New FA2HTGS Submission", startFa2htgs);
      SetObjectExtra (b, sfp, NULL);
    }
    if (startPhrap != NULL) {
      b = PushButton (c, "New PHRAP Submission", startPhrap);
      SetObjectExtra (b, sfp, NULL);
    }
    if (buildContig != NULL) {
      b = PushButton (c, "Read CONTIG Instructions", buildContig);
      SetObjectExtra (b, sfp, NULL);
    }
    b = PushButton (c, "Start New Submission", startNew);
    SetObjectExtra (b, sfp, NULL);
    b = PushButton (c, "Read Existing Record", readExisting);
    SetObjectExtra (b, sfp, NULL);
    if (fetchFromNet != NULL) {
      b = PushButton (c, "Download From Entrez", fetchFromNet);
      SetObjectExtra (b, sfp, NULL);
    }
    b = PushButton (c, "Show Help", showHelp);
    SetObjectExtra (b, sfp, NULL);
    if (createSubmissionTemplate != NULL)
    {
      b = PushButton (c, "Submission Template", createSubmissionTemplate);
      SetObjectExtra (b, sfp, NULL);
    }
    b = PushButton (c, "Quit Program", quitProgram);
    SetObjectExtra (b, sfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) c, (HANDLE) k, NULL);

    RealizeWindow (w);

    if (activateForm != NULL) {
      SetActivate (w, activateForm);
    }
  }
  return (ForM) w;
}

typedef struct formatform {
  FORM_MESSAGE_BLOCK

  GrouP           package;
  GrouP           format;
  GrouP           submType;
  ButtoN          alignmentButton;
  ButtoN          originalButton;
  ButtoN          tpaButton;
  TexT            numseqs;

  Int2            restoreFormatTo;
} FormatForm, PNTR FormatFormPtr;

static Boolean allowGenomicPlusCDNA = FALSE;

static void FormatBlockPtrToFormatForm (ForM f, Pointer data)

{
  FormatBlockPtr  fbp;
  FormatFormPtr   ffp;
  Char            str [32];

  ffp = (FormatFormPtr) GetObjectExtra (f);
  fbp = (FormatBlockPtr) data;
  if (ffp == NULL) return;
  if (fbp != NULL) {
    if (fbp->seqPackage > 0 && fbp->seqPackage <= NUM_SEQ_PKG) {
      if ((! allowGenomicPlusCDNA) && fbp->seqPackage >= SEQ_PKG_GENOMICCDNA) {
        SafeSetValue (ffp->package, fbp->seqPackage - 1);
      } else {
        SafeSetValue (ffp->package, fbp->seqPackage);
      }
      if (fbp->seqPackage <= SEQ_PKG_GENOMICCDNA || fbp->seqPackage == SEQ_PKG_GENBANK) {
        SafeDisable (ffp->alignmentButton);
      } else {
        SafeEnable (ffp->alignmentButton);
      }
    } else {
      SafeSetValue (ffp->package, SEQ_PKG_SINGLE);
      SafeDisable (ffp->alignmentButton);
    }
    if (fbp->seqFormat > 0 && fbp->seqFormat <= NUM_SEQ_FMT) {
      SafeSetValue (ffp->format, fbp->seqFormat);
    } else {
      SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    }
    if (fbp->numSeqs > 0) {
      IntToStr (fbp->numSeqs, str, 0, sizeof (str));
      SafeSetTitle (ffp->numseqs, str);
    } else {
      SafeSetTitle (ffp->numseqs, "");
    }
  } else {
    SafeSetValue (ffp->package, SEQ_PKG_SINGLE);
    SafeDisable (ffp->alignmentButton);
    SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    SafeSetTitle (ffp->numseqs, "");
    ffp->restoreFormatTo = SEQ_FMT_FASTA;
  }
}

static Pointer FormatFormToFormatBlockPtr (ForM f)

{
  FormatBlockPtr  fbp;
  FormatFormPtr   ffp;
  Char            str [32];
  Int2            val;

  fbp = NULL;
  ffp = (FormatFormPtr) GetObjectExtra (f);
  if (ffp == NULL) return NULL;
  fbp = (FormatBlockPtr) MemNew (sizeof (FormatBlock));
  if (fbp == NULL) return NULL;
  fbp->seqPackage = GetValue (ffp->package);
  if ((! allowGenomicPlusCDNA) && fbp->seqPackage >= SEQ_PKG_GENOMICCDNA) {
    (fbp->seqPackage)++;
  }
  fbp->seqFormat = GetValue (ffp->format);
  fbp->submType = GetValue (ffp->submType);
  GetTitle (ffp->numseqs, str, sizeof (str));
  if (StrToInt (str, &val) && val > 0) {
    fbp->numSeqs = val;
  } else {
    fbp->numSeqs = 0;
  }
  return (Pointer) fbp;
}

static void EnableOrDisableFormats (GrouP g)

{
  FormatFormPtr  ffp;
  Int2           val;

  ffp = (FormatFormPtr) GetObjectExtra (g);
  if (ffp == NULL) return;
  val = GetValue (g);
  if ((! allowGenomicPlusCDNA) && val >= SEQ_PKG_GENOMICCDNA) {
    val++;
  }
  if (val <= SEQ_PKG_GENOMICCDNA || val == SEQ_PKG_GENBANK) {
    if (Enabled (ffp->alignmentButton)) {
      ffp->restoreFormatTo = GetValue (ffp->format);
    }
    SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    SafeDisable (ffp->alignmentButton);
  } else {
    if (! Enabled (ffp->alignmentButton)) {
      SafeSetValue (ffp->format, ffp->restoreFormatTo);
    }
    SafeEnable (ffp->alignmentButton);
  }
}

static Boolean ExportTemplateMenu (ForM f, CharPtr filename)
{
  WindoW                w;
  GrouP                 h, g1, g2, c;
  ButtoN                b;
  DialoG                org_dlg;
  TexT                  comment_txt;
  ModalAcceptCancelData acd;
  SeqEntryPtr           sep;
  BioseqSetPtr          bssp;
  BioSourcePtr          biop;
  SeqDescrPtr           sdp;
  CharPtr               org_name;
  Boolean               done;
  
  if (ANS_NO == Message (MSG_YN, "Do you want to add an organism name and comment before saving the template?"))
  {
    bssp = BioseqSetNew ();
    bssp->_class = BioseqseqSet_class_not_set;
    sep = SeqEntryNew ();
    sep->choice = 2;
    sep->data.ptrvalue = bssp;
       
    done = ExportSubmitterBlockTemplate (sep, NULL);
    if (!done)
    {
      /* if done were TRUE, sep would have been freed as part of the new SeqSubmit */
      SeqEntryFree (sep);
    }
    return done;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, "Submission Template", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g1 = NormalGroup (h, 1, 0, "Organism Name", programFont, NULL);
  org_dlg = OrganismSelectionDialog (g1, "");
  g2 = NormalGroup (h, 2, 0, "Comment", programFont, NULL);
  comment_txt = DialogText (g2, "", 30, NULL);  
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2,
                              (HANDLE) c, (HANDLE) NULL);

  Show (w);
  Select (w);
  done = FALSE;
  while (!done)
  { 
    acd.accepted = FALSE;
    acd.cancelled = FALSE;
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
  
    if (acd.cancelled)
    {
      done = TRUE;
    }
    else
    {
      bssp = BioseqSetNew ();
      sep = SeqEntryNew ();
      sep->choice = 2;
      sep->data.ptrvalue = bssp;
      bssp->_class = BioseqseqSet_class_not_set;
    
      org_name = DialogToPointer (org_dlg);
      if (!StringHasNoText (org_name))
      {
        biop = BioSourceNew ();
        biop->org = OrgRefNew ();
        biop->org->taxname = org_name;
        sdp = CreateNewDescriptor (sep, Seq_descr_source);
        sdp->data.ptrvalue = biop;
      }
      else
      {
        org_name = MemFree (org_name);
      }
    
      sdp = NULL;   
      if (!TextHasNoText (comment_txt))
      {
        sdp = CreateNewDescriptor (sep, Seq_descr_comment);
        sdp->data.ptrvalue = SaveStringFromText (comment_txt);
      }
    
      done = ExportSubmitterBlockTemplate (sep, sdp);
      if (!done)
      {
        /* if done were TRUE, sep would have been freed as part of the new SeqSubmit */
        SeqEntryFree (sep);
      }
    }
  }
  Remove (w);
  return acd.accepted;
}

static void FormatFormMessage (ForM f, Int2 mssg)

{
  switch (mssg) {
    case VIB_MSG_EXPORT :
      ExportTemplateMenu (NULL, NULL);
      break;
    case VIB_MSG_CUT :
      StdCutTextProc (NULL);
      break;
    case VIB_MSG_COPY :
      StdCopyTextProc (NULL);
      break;
    case VIB_MSG_PASTE :
      StdPasteTextProc (NULL);
      break;
    case VIB_MSG_DELETE :
      StdDeleteTextProc (NULL);
      break;
    default :
      DefaultMessageProc (f, mssg);
      break;
  }
}

static void CreateFormatFormMenus (WindoW w)
{
#ifndef WIN_MAC  
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Export Template...", bfp, VIB_MSG_EXPORT);
  }
#endif
}

static void InitFormatFormActivate (WindoW w)

{
  IteM           exportItm;
  FormatFormPtr  ffp;

  ffp = (FormatFormPtr) GetObjectExtra (w);
  if (ffp != NULL) {
    if (ffp->activate != NULL) {
      ffp->activate (w);
    }
    exportItm = FindFormMenuItem ((BaseFormPtr) ffp, VIB_MSG_EXPORT);
    SafeSetTitle (exportItm, "Export Template...");
    SafeEnable (exportItm);
  }
}

extern ForM CreateFormatForm (Int2 left, Int2 top, CharPtr title,
                              BtnActnProc goToNext,
                              BtnActnProc goBack,
                              WndActnProc activateForm)

{
  ButtoN         b;
  GrouP          c;
  FormatFormPtr  ffp;
  GrouP          g1, g2, g3;
  GrouP          h;
  PrompT         ppt;
  Char           str [32];
  WindoW         w;

  w = NULL;
  ffp = MemNew (sizeof (FormatForm));
  if (ffp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->toform = FormatBlockPtrToFormatForm;
    ffp->fromform = FormatFormToFormatBlockPtr;
    ffp->formmessage = FormatFormMessage;
    ffp->exportform = ExportTemplateMenu;

    SetGroupSpacing (w, 10, 10);

    CreateFormatFormMenus (w);

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 3, 10);

    g1 = HiddenGroup (h, 2, 0, NULL);

    allowGenomicPlusCDNA = FALSE;
    if (GetAppParam ("SEQUIN", "SETTINGS", "GENOMICPLUSTRANSCRIPTS", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        allowGenomicPlusCDNA = TRUE;
      }
    }

    ppt = StaticPrompt (g1, "Submission type", 0, 0, programFont, 'l');
    ffp->package = HiddenGroup (g1, 2, 0, EnableOrDisableFormats);
    SetObjectExtra (ffp->package, ffp, NULL);
    RadioButton (ffp->package, "Single Sequence");
    RadioButton (ffp->package, "Segmented Sequence");
    RadioButton (ffp->package, "Gapped Sequence");
    if (allowGenomicPlusCDNA) {
      RadioButton (ffp->package, "Genomic + Transcripts");
    }
    RadioButton (ffp->package, "Population Study");
    RadioButton (ffp->package, "Phylogenetic Study");
    RadioButton (ffp->package, "Mutation Study");
    RadioButton (ffp->package, "Environmental Samples");
    RadioButton (ffp->package, "Batch Submission");
    SetValue (ffp->package, SEQ_PKG_SINGLE);
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) ffp->package, NULL);

    g2 = HiddenGroup (h, 2, 0, NULL);

    ppt = StaticPrompt (g2, "Sequence data format", 0, 0, programFont, 'l');
    ffp->format = HiddenGroup (g2, -1, 0, NULL);
    SetObjectExtra (ffp->format, ffp, NULL);
    RadioButton (ffp->format, "FASTA (no alignment)");
    ffp->alignmentButton = RadioButton (ffp->format, "Alignment (FASTA+GAP, NEXUS, PHYLIP, etc.)");
    Disable (ffp->alignmentButton);
    SetValue (ffp->format, SEQ_FMT_FASTA);
    ffp->restoreFormatTo = SEQ_FMT_FASTA;
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) ffp->format, NULL);

    g3 = HiddenGroup (h, 2, 0, NULL);

    ppt = StaticPrompt (g3, "Submission category", 0, 0, programFont, 'l');
    ffp->submType = HiddenGroup (g3, -1, 0, NULL);
    SetObjectExtra (ffp->submType, ffp, NULL);
    ffp->originalButton = RadioButton (ffp->submType, "Original Submission");
    ffp->tpaButton = RadioButton (ffp->submType, "Third Party Annotation");
    SetValue (ffp->submType, SEQ_ORIG_SUBMISSION);
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) ffp->submType, NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    b = PushButton (c, " << Prev Form ", goBack);
    SetObjectExtra (b, ffp, NULL);
    b = PushButton (c, " Next Form >> ", goToNext);
    SetObjectExtra (b, ffp, NULL);

    AlignObjects (ALIGN_LEFT, (HANDLE) g1, (HANDLE) g2, (HANDLE) g3, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);

    RealizeWindow (w);

    ffp->activate = activateForm;
    SetActivate (w, InitFormatFormActivate);
  }
  return (ForM) w;
}

extern SequinBlockPtr SequinBlockFree (SequinBlockPtr sbp)

{
  if (sbp != NULL) {
    AuthorFree (sbp->contactperson);
    AuthListFree (sbp->citsubauthors);
    AffilFree (sbp->citsubaffil);
    MemFree (sbp->citsubtitle);
    DateFree (sbp->releasedate);
  }
  return NULL;
}

extern void ExciseString (CharPtr str, CharPtr from, CharPtr to)

{
  Char     ch;
  CharPtr  ptrf;
  CharPtr  ptrt;

  if (str == NULL || from == NULL || to == NULL) return;
  ptrf = StringISearch (str, from);
  if (ptrf == NULL) return;
  ptrt = StringISearch (ptrf, to);
  if (ptrt == NULL) return;
  ptrt += StringLen (to);
  ch = *ptrt;
  while (ch != '\0') {
    *ptrf = ch;
    ptrf++;
    ptrt++;
    ch = *ptrt;
  }
  *ptrf = '\0';
}

typedef struct geneextendlist {
  GeneRefPtr  grp;
  SeqLocPtr   slp;
  ObjMgrPtr   omp;
  Boolean     rsult;
  Char        label [41];
} GeneExtendList, PNTR GeneExtendPtr;

static Boolean GeneExtendFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  GeneExtendPtr  gep;
  GeneRefPtr     grp;
  Boolean        hasNulls;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  gep = (GeneExtendPtr) gcp->userdata;
  if (gep == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      omtp = ObjMgrTypeFind (gep->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        if (StringICmp (thislabel, gep->label) == 0) {
          if (SeqLocCompare (gep->slp, sfp->location) != SLC_NO_MATCH) {
            bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
            if (bsp != NULL) {
              slp = SeqLocMerge (bsp, sfp->location, gep->slp, TRUE, FALSE, FALSE);
              if (slp != NULL) {
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = slp;
                if (bsp->repr == Seq_repr_seg) {
                  slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                  sfp->location = SeqLocFree (sfp->location);
                  sfp->location = slp;
                  hasNulls = LocationHasNullsBetween (sfp->location);
                  sfp->partial = (sfp->partial || hasNulls);
                }
                FreeAllFuzz (slp);
                gep->rsult = TRUE;
              }
            }
          }
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

extern Boolean ExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp)

{
  GeneExtendList  gel;
  GatherScope     gs;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;

  if (grp == NULL || nsep == NULL || slp == NULL) return FALSE;
  gel.grp = grp;
  gel.slp = slp;
  gel.omp = ObjMgrGet ();
  gel.label [0] = '\0';
  gel.rsult = FALSE;
  omtp = ObjMgrTypeFind (gel.omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp != NULL && omtp->labelfunc != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_GENE;
      sfp->data.value.ptrvalue = (Pointer) grp;
      (*(omtp->labelfunc)) ((Pointer) sfp, gel.label, 40, OM_LABEL_CONTENT);
      sfp->data.value.ptrvalue = NULL;
      SeqFeatFree (sfp);
    }
  }
  MemSet ((Pointer)(&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherSeqEntry (nsep, (Pointer) &gel, GeneExtendFunc, &gs);
  return gel.rsult;
}

/*=====================================================================*/
/*                                                                     */
/* CreateGeneAndProtFeats() -                                          */
/*                                                                     */
/*=====================================================================*/

static void CreateGeneAndProtFeats (SeqEntryPtr nsep, SeqEntryPtr psep,
                                    SeqLocPtr slp, CdRegionPtr crp, CharPtr title,
                                    CharPtr best, size_t maxsize, CharPtr PNTR ttl)

{
  BioseqPtr   nbsp;
  BioseqPtr   pbsp;
  SeqFeatPtr  sfp;
  ProtRefPtr  prp = NULL;
  Boolean     partial5, partial3;

  if (nsep != NULL && psep != NULL && slp != NULL && crp != NULL && title != NULL) {
    if (best != NULL) {
      best [0] = '\0';
    }
    if (IS_Bioseq (nsep) && IS_Bioseq (psep)) {
      nbsp = (BioseqPtr) nsep->data.ptrvalue;
      pbsp = (BioseqPtr) psep->data.ptrvalue;
      if (nbsp != NULL && pbsp != NULL) {

        AddGeneFeatureFromTitle (nsep, title, slp);
      
        sfp = AddProteinFeatureFromDefline (psep, title);
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT) {
          prp = sfp->data.value.ptrvalue;
          CheckSeqLocForPartial (slp, &partial5, &partial3);
          SetSeqLocPartial (sfp->location, partial5, partial3);
          sfp->partial = partial5 | partial3;
        }
          
        if (prp != NULL && best != NULL)
        {
          if (prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue))
          {
	          StringNCpy_0 (best, prp->name->data.ptrvalue, maxsize);
          }
          else if (!StringHasNoText (prp->desc))
          {
	          StringNCpy_0 (best, prp->desc, maxsize);
          }
        }

        AddCodingRegionFieldsFromProteinTitle (crp, title, ttl);
      }
    }
  }
}

static Boolean  intBoxUp;
static Boolean  intBoxRsult;

static void AcceptAskProc (ButtoN b)

{
  intBoxRsult = TRUE;
  intBoxUp = FALSE;
}

static void CancelAskProc (ButtoN b)

{
  intBoxRsult = FALSE;
  intBoxUp = FALSE;
}

static SeqLocPtr AskForInterval (SeqEntryPtr sep, BioseqPtr nuc, BioseqPtr prot)

{
  GrouP       c;
  DialoG      d;
  GrouP       g;
  GrouP       m;
  SeqIdPtr    sip;
  SeqLocPtr   slp;
  Char        str [128];
  ValNodePtr  vnp;
  WindoW      w;

  slp = NULL;
  if (sep == NULL || nuc == NULL || prot == NULL) return NULL;

  if (GetAppParam ("SEQUIN", "PREFERENCES", "ASKIFSUGGESTFAILED", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      sip = SeqIdFindWorst (prot->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      Message (MSG_POSTERR, "Suggest failure for %s", str);
      return NULL;
    }
  }

  w = MovableModalWindow (-50, -33, -10, -10, "Enter coding region interval", NULL);
  g = HiddenGroup (w, -1, 0, NULL);
  m = NULL;
  SetGroupSpacing (g, 3, 10);
  if (prot->descr != NULL) {
    vnp = ValNodeFindNext (prot->descr, NULL, Seq_descr_title);
    if (vnp != NULL && vnp->data.ptrvalue != NULL) {
      m = MultiLinePrompt (g, (CharPtr) vnp->data.ptrvalue, stdCharWidth * 28, programFont);
    }
  }
  d = CreateIntervalEditorDialog (g, NULL, 4, 2, sep, TRUE, FALSE);
  c = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  DefaultButton (c, "Accept", AcceptAskProc);
  PushButton (c, "Cancel", CancelAskProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) d, (HANDLE) c, (HANDLE) m, NULL);
  Show (w);
  Select (w);
  intBoxUp = TRUE;
  intBoxRsult = FALSE;
  while (intBoxUp) {
    ProcessEventOrIdle ();
  }
  ProcessAnEvent ();
  if (intBoxRsult) {
    slp = (SeqLocPtr) DialogToPointer (d);
  }
  Remove (w);
  return slp;
}

extern Boolean AutomaticProteinProcess (SeqEntryPtr esep, SeqEntryPtr psep,
                                        Int2 code, Boolean makeMRNA, SeqLocPtr use_this)

{
  SeqFeatPtr   cds;
  CdRegionPtr  crp;
  Char         mRnaName [128];
  BioseqPtr    nbsp;
  SeqEntryPtr  nsep;
  BioseqPtr    pbsp;
  SeqFeatPtr   rna;
  RnaRefPtr    rrp;
  SeqLocPtr    slp;
  CharPtr      ttl;
  ValNodePtr   vnp;
  CharPtr      vnpstr;
  Boolean      partial5, partial3;

  if (esep == NULL || psep == NULL) return FALSE;

  nsep = FindNucSeqEntry (esep);
  if (nsep == NULL || (! IS_Bioseq (nsep)) || (! IS_Bioseq (psep))) return FALSE;

  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  pbsp = (BioseqPtr) psep->data.ptrvalue;
  if (nbsp == NULL || pbsp == NULL) return FALSE;

  cds = NULL;
  WatchCursor ();
  Update ();
  if (use_this == NULL) {
    slp = PredictCodingRegion (nbsp, pbsp, code);
    if (slp == NULL) {
      ArrowCursor ();
      Update ();
      slp = AskForInterval (nsep, nbsp, pbsp);
    }
  } else {
    slp = use_this;
  }
  if (slp == NULL) return FALSE;

  mRnaName [0] = '\0';
  ttl = NULL;
  crp = CreateNewCdRgn (0, FALSE, code);
  if (crp != NULL) {
    if (pbsp->descr != NULL) {
      vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        vnpstr = (CharPtr) vnp->data.ptrvalue;
        CreateGeneAndProtFeats (nsep, psep, slp, crp, vnpstr, mRnaName, sizeof (mRnaName), &ttl);
        TrimSpacesAroundString (vnpstr);
        if (StringHasNoText (vnpstr)) {
          ValNodeExtract (&(pbsp->descr), Seq_descr_title);
        }
      }
    }
    if (makeMRNA) {
      rrp = RnaRefNew ();
      if (rrp != NULL) {
        rrp->type = 2;
        if (! StringHasNoText (mRnaName)) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave (mRnaName);
        }
        rna = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, NULL);
        if (rna != NULL) {
          rna->data.value.ptrvalue = (Pointer) rrp;
          rna->location = SeqLocFree (rna->location);
          rna->location = AsnIoMemCopy ((Pointer) slp,
                                        (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
        }
      }
    }
    cds = CreateNewFeature (nsep, NULL, SEQFEAT_CDREGION, NULL);
    if (cds != NULL) {
      cds->data.value.ptrvalue = (Pointer) crp;
      cds->location = SeqLocFree (cds->location);
      cds->location = slp;
      slp = NULL;
      CheckSeqLocForPartial (cds->location, &partial5, &partial3);
      cds->partial |= partial5 | partial3;
      SetSeqFeatProduct (cds, pbsp);
      if (! StringHasNoText (ttl)) {
        cds->comment = ttl;
      }
    }
  }

  SeqLocFree (slp);
  return TRUE;
}

typedef struct fa2htgsform {
  FORM_MESSAGE_BLOCK

  SeqSubmitPtr       ssp;
  SeqEntryPtr        sep;

  GrouP              templateblock;
  GrouP              fastablock;
  GrouP              orderblock;
  GrouP              controlblock;
  GrouP              contigtype;

  DialoG             contigorder;
  EnumFieldAssocPtr  alists [1];
  GrouP              htgsphase;
  ButtoN             draft;
  ButtoN             fulltop;
  ButtoN             activefin;
  TexT               orgname;
  TexT               seqname;
  ButtoN             update;
  TexT               accession;
  TexT               knownlength;
  TexT               gaplength;
  TexT               remark;
  TexT               clone;
  TexT               strain;
  TexT               cultivar;
  TexT               chromosome;
  TexT               title;
  DialoG             secondaries;

  SeqEntryPtr        seplist;

  ButtoN             okBtn;
  BtnActnProc        finish;
  BtnActnProc        cancel;
  Boolean            readPhrap;
  Boolean            buildContig;

} Fa2htgsForm, PNTR Fa2htgsFormPtr;

/*------------- MakeAc2GBSeqId() -----------------------*/
/***************************************************************
*   MakeAc2GBSeqId:
*   -- return NULL if acnum == null
*                                             Hsiu-Chuan 4-18-97
****************************************************************/
static SeqIdPtr  SqnMakeAc2GBSeqId(CharPtr accession)
{
   TextSeqIdPtr tsip;
   SeqIdPtr sip;

   if (accession == NULL || *accession == '\0')
      return NULL;

   sip = ValNodeNew(NULL);
   sip->choice = SEQID_GENBANK;
   tsip = TextSeqIdNew();
   sip->data.ptrvalue = tsip;
   tsip->accession = StringSave(accession);

   return sip;

} /* MakeAc2GBSeqId */

/*----------- AddExtraAc2Entry() ----------------------------*/
/***************************************************************
*   AddExtraAc2Entry:
*                                             Hsiu-Chuan 4-11-97, modified by JK
****************************************************************/
static void SqnAddDraft2Entry (SeqEntryPtr entry, CharPtr keyword)

{
   BioseqPtr  bsp;
   ValNodePtr vnp;
   GBBlockPtr gbp;

   if (entry == NULL) return;

   bsp = (BioseqPtr)(entry->data.ptrvalue);

   for (gbp= NULL, vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
   {
       if (vnp->choice == Seq_descr_genbank)
       {
          gbp = vnp->data.ptrvalue;
          break;
       }
   }

   if (gbp == NULL)
   {
      vnp = (ValNodePtr) NewDescrOnSeqEntry (entry, Seq_descr_genbank);
      gbp = GBBlockNew();
      vnp->data.ptrvalue = (Pointer)gbp;
   }

   if (gbp != NULL) {
      ValNodeCopyStr (&gbp->keywords, 0, keyword);
   }
}

static Boolean SqnAddExtraAc2Entry (SeqEntryPtr entry , ValNodePtr extra_accs )
{
   BioseqPtr  bsp;
   ValNodePtr vnp;
   GBBlockPtr gbp;
   Char       acnum[17];
   CharPtr    p;
   Int4       i, j;
   SeqHistPtr shp;
   SeqIdPtr   sip;
   ValNodePtr tmp;

   if ((entry == NULL) || (extra_accs == NULL))
      return FALSE;

   bsp = (BioseqPtr)(entry->data.ptrvalue);

   for (gbp= NULL, vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
   {
       if (vnp->choice == Seq_descr_genbank)
       {
          gbp = vnp->data.ptrvalue;
          break;
       }
   }

   shp = bsp->hist; 

   if (gbp == NULL)
   {
      vnp = (ValNodePtr) NewDescrOnSeqEntry (entry, Seq_descr_genbank);
      gbp = GBBlockNew();
      vnp->data.ptrvalue = (Pointer)gbp;
   }
   
   for (tmp = extra_accs; tmp != NULL; tmp = tmp->next)
   {
       p = (CharPtr) tmp->data.ptrvalue;
       if (p == NULL) continue;
       for (i = 0; isalnum((Int4)(*p)) && *p != '\0'; ++p, ++i)
           acnum[i] = *p;
       acnum[i] = '\0'; 
               /* check one_letter+5digits or two_letter+6digits */
       if (i == 6 || i == 8)
       {
          if (!isalpha((Int4)(acnum[0])) || (!(isdigit((Int4)(acnum[1])) && i == 6) &&
              !(isalpha((Int4)(acnum[1])) && i == 8)))
          {
             ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
             return FALSE;
          }

          for (j = 2; j < i; ++j)
          {
              if (!(isdigit((Int4)(acnum[j]))))
              {
                 ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
                 return FALSE;
              }
          }

          ValNodeCopyStr(&gbp->extra_accessions, 0, acnum);
          sip = SqnMakeAc2GBSeqId (acnum);
          if (shp == NULL)
          {
             shp = SeqHistNew();
             bsp->hist = shp;
          }
          ValNodeLink(&shp->replace_ids, sip);
       }
       else
       {
          ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
          return FALSE;
       }

       while (!isalnum((Int4)(*p)) && *p != '\0')
           ++p;
   }

   return TRUE;

} /* AddExtraAc2Entry */

static void RescueSeqGraphs (BioseqPtr bsp, Int2 index, ValNodePtr PNTR vnpp)

{
  SeqAnnotPtr   nextsap;
  SeqGraphPtr   nextsgp;
  Pointer PNTR  prevsap;
  Pointer PNTR  prevsgp;
  SeqAnnotPtr   sap;
  SeqGraphPtr   sgp;

  if (bsp == NULL || vnpp == NULL) return;
  sap = bsp->annot;
  prevsap = (Pointer PNTR) &(bsp->annot);
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        *(prevsgp) = sgp->next;
        sgp->next = NULL;
        ValNodeAddPointer (vnpp, index, (Pointer) sgp);
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static SeqAnnotPtr NewSeqAnnotType3 (CharPtr name, SeqGraphPtr sgp)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (! StringHasNoText (name)) {
    SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static void OffsetAndLinkSeqGraph (BioseqPtr bsp, SeqGraphPtr sgp, Int2 index)

{
  DeltaSeqPtr  dsp;
  SeqGraphPtr  lastsgp;
  Int4         len;
  SeqLitPtr    litp;
  SeqAnnotPtr  sap;
  SeqIntPtr    sintp;
  SeqLocPtr    slp;

  if (bsp == NULL || sgp == NULL || index < 1) return;
  len = 0;
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
         dsp != NULL && index > 1; dsp = dsp->next, index--) {
      if (dsp->choice == 1) {
        len += SeqLocLen ((SeqLocPtr) dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          len += litp->length;
        }
      }
    }
  }
  slp = sgp->loc;
  if (slp != NULL && slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp != NULL) {
      sintp->from += len;
      sintp->to += len;
      sintp->id = SeqIdFree (sintp->id);
      sintp->id = SeqIdDup (bsp->id);
    }
  }
  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 3) {
      for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
        continue;
      }
      lastsgp->next = sgp;
      break;
    }
  }
  if (sap == NULL) {
    if (bsp->annot != NULL) {
      for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
        continue;
      }
      sap->next = NewSeqAnnotType3 ("Graphs", sgp);
    } else {
      bsp->annot = NewSeqAnnotType3 ("Graphs", sgp);
    }
  }
}

static CharPtr phrapBoilerPlate = "Sequence Quality Assessment:~ \
This entry has been annotated with sequence quality~ \
estimates computed by the Phrap assembly program.~ \
All manually edited bases have been reduced to quality zero.~ \
Quality levels above 40 are expected to have less than ~ \
1 error in 10,000 bp.~ \
Base-by-base quality values are not generally visible from the~ \
GenBank flat file format but are available as part~ \
of this entry's ASN.1 file.~----------------------~";

static void ProcessFa2htgs (Fa2htgsFormPtr ffp, SeqSubmitPtr ssp)

{
  SeqEntryPtr  sep, oldsep, the_entry, nextsep;
  NCBISubPtr nsp;
  Int2 htgs_phase = -1;
  Uint1 tech;
  CharPtr seqname = NULL, accession = NULL, orgname = NULL;
  CharPtr clone = NULL, strain = NULL, cultivar = NULL, chromosome = NULL;
  CharPtr remark = NULL, title = NULL, seqbuf = NULL;
  Int4 length = 0, cumlength = 0, gaplen;
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqLitPtr slp;
  ValNodePtr vnp, PNTR prevpnt, next, extra_accs;
  Boolean lastwasraw, draft, fulltop, activefin, usedelta = FALSE;
  Char str [64];
  long int val;
  Int2 index = 0;
  ValNodePtr rescuedsgps = NULL;
  ValNodePtr seqlitlist = NULL;
  IntFuzzPtr ifp;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;
  DatePtr dp;

  if (ffp == NULL || ssp == NULL) return;

  htgs_phase = GetValue (ffp->htgsphase) - 1;
  orgname = SaveStringFromText (ffp->orgname);
  seqname = SaveStringFromText (ffp->seqname);
  if (GetStatus (ffp->update)) {
    accession = SaveStringFromText (ffp->accession);
  }
  clone = SaveStringFromText (ffp->clone);
  strain = SaveStringFromText (ffp->strain);
  cultivar = SaveStringFromText (ffp->cultivar);
  chromosome = SaveStringFromText (ffp->chromosome);
  remark = SaveStringFromText (ffp->remark);
  title = SaveStringFromText (ffp->title);
  extra_accs = DialogToPointer (ffp->secondaries);
  draft = GetStatus (ffp->draft);
  fulltop = GetStatus (ffp->fulltop);
  activefin = GetStatus (ffp->activefin);

  length = 0;
/* may need to really calculate length */
  GetTitle (ffp->knownlength, str, sizeof (str));
  if (! StringHasNoText (str)) {
    if (sscanf (str, "%ld", &val) == 1 && val > 0) {
      length = (Int4) val;
    }
  }

  gaplen = 0;
/* now usually filling in with gaps of 100 bases */
  GetTitle (ffp->gaplength, str, sizeof (str));
  if (! StringHasNoText (str)) {
    if (sscanf (str, "%ld", &val) == 1 && val > 0) {
      gaplen = (Int4) val;
    }
  }

/* modified from fa2htgs */
   oldsep = (SeqEntryPtr)(ssp->data);  /* clear out template */
   ssp->data = NULL;
   MemFree(ssp->sub->tool);
   sprintf (str, "Sequin %s", SEQUIN_APPLICATION);
   ssp->sub->tool = StringSave (str);
   nsp = MemNew(sizeof(NCBISub));
   nsp->ssp = ssp;
   nsp->submittor_key = StringSave (genomeCenter);
   /*
   MemFree(ssp->sub->cit->descr);
   ssp->sub->cit->descr = remark;
   */

   cumlength = 0;
   index = 0;

   sep = ffp->seplist;
   if (sep != NULL && sep->next != NULL) {
     usedelta = TRUE;
   }

   if (ffp->buildContig) {
     ssp->data = (Pointer) ffp->seplist;
     the_entry = ffp->seplist;
     sep = the_entry;

     oip = ObjectIdNew ();
     oip->str = StringSave ("info");
     uop = UserObjectNew ();
     uop->type = oip;
     uop->_class = StringSave ("Genomes");

     oip = ObjectIdNew ();
     oip->id = 0;
     ufp = UserFieldNew ();
     ufp->choice = 2;
     ufp->data.intvalue = 0;
     ufp->label = oip;

     uop->data = ufp;

     if (IS_Bioseq (sep)) {
       bsp = (BioseqPtr) sep->data.ptrvalue;
       vnp = SeqDescrNew (NULL);
       vnp->choice = Seq_descr_user;
       vnp->data.ptrvalue = (Pointer) uop;
       vnp->next = bsp->descr;
       bsp->descr = vnp;
       cumlength = bsp->length;
     }

   }
   else if (htgs_phase < 3 || usedelta)
   {
      the_entry = AddDeltaSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = ffp->seplist;
      lastwasraw = FALSE;
      while (sep != NULL)
      {
         nextsep = sep->next;
         sep->next = NULL;
         bsp = (BioseqPtr)(sep->data.ptrvalue);
         if (bsp->repr == Seq_repr_raw)
         {
            if (lastwasraw) {
               slp = AddFakeGapToDeltaSeq(nsp, the_entry, gaplen);
               ValNodeAddPointer (&seqlitlist, 0, (Pointer) slp);
               index++;
               cumlength += gaplen;
            }
            BioseqRawConvert(bsp, Seq_code_iupacna);
            seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
            slp = AddLiteralToDeltaSeq(nsp, the_entry,
               bsp->length);
            AddBasesToLiteral(nsp, slp, seqbuf);
            MemFree(seqbuf);
            lastwasraw = TRUE;
            index++;
         }
         else
         {
            if (bsp->length < 0)
               bsp->length = 0;  /* -1 may be set */
            AddGapToDeltaSeq(nsp, the_entry,
               bsp->length);
            lastwasraw = FALSE;
            index++;
         }
         cumlength += bsp->length;
         RescueSeqGraphs (bsp, index, &rescuedsgps);
         SeqEntryFree(sep);
         sep = nextsep;
      }
   }
   else
   {
    the_entry = AddSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = ffp->seplist;
      nextsep = sep->next;
      sep->next = NULL;
      bsp = (BioseqPtr)(sep->data.ptrvalue);
      if (bsp->repr == Seq_repr_raw)
      {
         BioseqRawConvert(bsp, Seq_code_iupacna);
         seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
         AddBasesToBioseq(nsp, the_entry, seqbuf);
         MemFree(seqbuf);
         index++;
      }
      cumlength += bsp->length;
      RescueSeqGraphs (bsp, index, &rescuedsgps);
      SeqEntryFree(sep);
      if (nextsep != NULL) {
        ErrPostEx (SEV_ERROR ,0, 0, "Only the first contig was used for HTGS 3");
      }
      if (length > 0 && length != cumlength) {
        ErrPostEx (SEV_ERROR ,0, 0, "Length is exactly %ld, not %ld",
                   (long) cumlength, (long) length);
        length = cumlength;
      }
   }

    /* get data from template: pub, organism, and comment */
   if (IS_Bioseq(oldsep))
   {
      bsp = (BioseqPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bsp->descr);
   }
   else
   {
      bssp = (BioseqSetPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bssp->descr);
   }

   bsp = (BioseqPtr)(the_entry->data.ptrvalue);
   if (bsp != NULL) {
     bsp->length = MAX (cumlength, length);
   }

   for (vnp = *prevpnt; vnp != NULL; vnp = next)
   {
      next = vnp->next;
      if (vnp->choice == Seq_descr_pub)
      {
         *prevpnt = next;
         vnp->next = NULL;
         ValNodeLink(&(bsp->descr), vnp);
      }
      else
         prevpnt = &(vnp->next);
   }
   if (remark != NULL) {
     vnp = SeqDescrNew (NULL);
     if (vnp != NULL) {
       vnp->choice = Seq_descr_comment;
       vnp->data.ptrvalue = remark;
       ValNodeLink(&(bsp->descr), vnp);
     }
   }

   SeqEntryFree(oldsep);

   AddOrganismToEntryNew(nsp, the_entry, orgname, NULL, NULL, NULL,
                         NULL, NULL, NULL, NULL);

   AddGenomeToEntry(nsp, the_entry, 1);
   if (clone != NULL)
      AddSubSourceToEntry(nsp, the_entry, 3, clone);
   if (chromosome != NULL)
       AddSubSourceToEntry(nsp, the_entry, 1, chromosome);
   if (strain != NULL)
      AddOrgModToEntry(nsp, the_entry, 2, strain);
   if (cultivar != NULL)
      AddOrgModToEntry(nsp, the_entry, ORGMOD_cultivar, cultivar);
   if (title != NULL)
      AddTitleToEntry(nsp, the_entry, title);
   if (ffp->readPhrap) {
      AddCommentToEntry(nsp, the_entry, phrapBoilerPlate);
   }

   if (extra_accs != NULL) {
      SqnAddExtraAc2Entry(the_entry, extra_accs);
   }

   if (draft) {
      SqnAddDraft2Entry(the_entry, "HTGS_DRAFT");
   }
   if (fulltop) {
      SqnAddDraft2Entry(the_entry, "HTGS_FULLTOP");
   }
   if (activefin) {
      SqnAddDraft2Entry(the_entry, "HTGS_ACTIVEFIN");
   }

   AddBiomolToEntry(nsp, the_entry, 1);
   if (ffp->buildContig) {
   } else {
     switch (htgs_phase) {
       case 0 :
         tech = MI_TECH_htgs_0;
         break;
       case 1 :
         tech = MI_TECH_htgs_1;
         break;
       case 2 :
         tech = MI_TECH_htgs_2;
         break;
       case 3 :
         tech = MI_TECH_htgs_3;
         break;
       default :
         tech = MI_TECH_htgs_3;
         break;
     }
     AddTechToEntry(nsp, the_entry, tech);
   }
   vnp = NewDescrOnSeqEntry (the_entry, Seq_descr_create_date);
   if (vnp != NULL) {
     dp = DateCurr ();
     vnp->data.ptrvalue = (Pointer) dp;
   }

   if (bsp != NULL) {
     for (vnp = rescuedsgps; vnp != NULL; vnp = vnp->next) {
       OffsetAndLinkSeqGraph (bsp, (SeqGraphPtr) vnp->data.ptrvalue, (Int2) vnp->choice);
       vnp->data.ptrvalue = NULL;
     }
   }
   rescuedsgps = ValNodeFreeData (rescuedsgps);

   for (vnp = seqlitlist; vnp != NULL; vnp = vnp->next) {
     slp = (SeqLitPtr) vnp->data.ptrvalue;
     if (slp != NULL) {
       ifp = IntFuzzNew();
       ifp->choice = 4;    /* lim - unk*/
       slp->fuzz = ifp;
     }
   }
   seqlitlist = ValNodeFree (seqlitlist);

   MemFree (nsp);

}

static Pointer Fa2htgsToSeqSubmitPtr (ForM f)

{
  Fa2htgsFormPtr  ffp;
  SeqSubmitPtr    ssp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (f);
  if (ffp == NULL) return NULL;
  ProcessFa2htgs (ffp, ffp->ssp);
  ssp = ffp->ssp;
  ffp->ssp = NULL;
  ffp->sep = NULL;
  ffp->seplist = NULL;
  return (Pointer) ssp;
}

static void Fa2htgsFormMessage (ForM f, Int2 mssg)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (f);
  if (ffp) {
    switch (mssg) {
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

#ifndef WIN_MAC
static void CreateFa2htgsFormMenus (WindoW w)

{
  BaseFormPtr  bfp;
  MenU         m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    AddAboutAndHelpMenuItems (m);
    FormCommandItem (m, "Quit", bfp, VIB_MSG_QUIT);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static void FillInFa2htgsFields (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Fa2htgsFormPtr  ffp;
  ValNodePtr      sdp = NULL;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  ffp = (Fa2htgsFormPtr) mydata;
  if (ffp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_comment) {
      SetTitle (ffp->remark, (CharPtr) sdp->data.ptrvalue);
    }
    sdp = sdp->next;
  }
}

static void ReadTemplate (ButtoN b)

{
  AsnIoPtr        aip;
  Fa2htgsFormPtr  ffp;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep;
  SeqSubmitPtr    ssp;
  Char            str [128];

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
    ffp->ssp = SeqSubmitFree (ffp->ssp);
    aip = AsnIoOpen (path, "r");
    if (aip != NULL) {
      ffp->ssp = SeqSubmitAsnRead (aip, NULL);
    }
    AsnIoClose (aip);
  }
  if (ffp->ssp != NULL) {
    SafeShow (ffp->fastablock);
    SafeDisable (b);
    ssp = ffp->ssp;
    if (ssp->datatype == 1) {
      sep = (SeqEntryPtr) ssp->data;
      if (sep != NULL) {
        SeqEntryToGeneticCode (sep, NULL, str, sizeof (str));
        SetTitle (ffp->orgname, str);
        SeqEntryExplore (sep, (Pointer) ffp, FillInFa2htgsFields);
      }
    }
  }
}

static Boolean Fa2htgsFormOkay (Fa2htgsFormPtr ffp)

{
  if (ffp == NULL) return FALSE;
  if (ffp->ssp == NULL) return FALSE;
  if (ffp->seplist == NULL) return FALSE;
  if (GetValue (ffp->htgsphase) == 0 && (! ffp->buildContig)) return FALSE;
  if (TextHasNoText (ffp->orgname)) return FALSE;
  if (TextHasNoText (ffp->seqname)) return FALSE;
  return TRUE;
}

static void ReadFastaHtgsFile (ButtoN b)

{
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  SeqEntryPtr     head;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep;
  CharPtr         ttl;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      head = NULL;
      sep = FastaToSeqEntry (fp, TRUE);
      while (sep != NULL) {
        ValNodeLink (&head, sep);
        sep = FastaToSeqEntry (fp, TRUE);
      }
      if (head != NULL) {
        ttl = SeqEntryGetTitle (head);
        SetTitle (ffp->title, ttl);
      }
      ffp->seplist = head;
    }
    FileClose (fp);
    if (ffp->seplist != NULL) {
      SafeShow (ffp->controlblock);
      SafeDisable (b);
      if (TextHasNoText (ffp->orgname)) {
        Select (ffp->orgname);
      } else {
        Select (ffp->seqname);
      }
      if (Fa2htgsFormOkay (ffp)) {
        SafeEnable (ffp->okBtn);
      } else {
        SafeDisable (ffp->okBtn);
      }
    }
  }
}

static EnumFieldAssocPtr MakePhrapAlists (SeqEntryPtr head)

{
  EnumFieldAssocPtr  alist = NULL;
  BioseqPtr          bsp;
  Int2               count;
  Int2               j;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  Char               str [128];

  if (head == NULL) return NULL;
  count = ValNodeLen (head);
  alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) (count + 4));
  if (alist == NULL) return NULL;

  j = 0;
  alist [j].name = StringSave ("                              ");
  alist [j].value = (UIEnum) 0;
  for (j = 1, sep = head; j <= count && sep != NULL; j++, sep = sep->next) {
    if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      str [30] = '\0';
      alist [j].name = StringSave (str);
      alist [j].value = (UIEnum) j;
    }
  }
  j = count + 1;
  alist [j].name = NULL;
  alist [j].value = (UIEnum) 0;


  return alist;
}

static void ReadAPhrapFile (ButtoN b)

{
  PopuP           control;
  Int2            count;
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  SeqEntryPtr     head;
  Int2            i;
  Char            path [PATH_MAX];
  RecT            r;
  TagListPtr      tlp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    WatchCursor ();
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      head = ReadPhrapFile (fp);
      ffp->seplist = head;
    }
    FileClose (fp);
    tlp = (TagListPtr) GetObjectExtra (ffp->contigorder);
    if (tlp != NULL) {
      ffp->alists [0] = MakePhrapAlists (ffp->seplist);
      tlp->alists = ffp->alists;
      for (i = 0; i < tlp->rows; i++) {
        control = (PopuP) tlp->control [i * MAX_TAGLIST_COLS + 0];
        if (control != NULL) {
          GetPosition (control, &r);
          Reset (control);
          InitEnumPopup (control, ffp->alists [0], NULL);
          SetEnumPopup (control, ffp->alists [0], 0);
          SetPosition (control, &r);
        }
      }
      count = ValNodeLen (ffp->seplist);
      for (i = 0; i < count; i++) {
        ValNodeCopyStr (&(tlp->vnp), 0, "0");
      }
      tlp->max = MAX ((Int2) 0, (Int2) (count - tlp->rows));
      CorrectBarMax (tlp->bar, tlp->max);
      CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
    }
    SafeShow (ffp->orderblock);
    SafeDisable (b);
    ArrowCursor ();
  }
}

static void PhrapOrderChosen (ButtoN b)

{
  SeqEntryPtr     PNTR collision;
  Int2            count;
  Fa2htgsFormPtr  ffp;
  Int2            i;
  Int2            j;
  SeqEntryPtr     lastsep;
  SeqEntryPtr     nextsep;
  Boolean         okay;
  SeqEntryPtr     PNTR order;
  SeqEntryPtr     sep;
  CharPtr         str;
  TagListPtr      tlp;
  int             val;
  ValNodePtr      vnp;
  /*
  Char            contigs [256];
  Char            tmp [256];
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  Boolean         notfirst;
  */

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  tlp = (TagListPtr) GetObjectExtra (ffp->contigorder);
  if (tlp == NULL) return;
  count = ValNodeLen (tlp->vnp);
  order = MemNew (sizeof (SeqEntryPtr) * (count + 1));
  if (order == NULL) return;
  collision = MemNew (sizeof (SeqEntryPtr) * (count + 1));
  if (collision == NULL) return;

  okay = TRUE;
  for (i = 0, vnp = tlp->vnp; i < count && vnp != NULL; i++, vnp = vnp->next) {
    str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    if (str != NULL && sscanf (str, "%d", &val) == 1 && val > 0) {
      val--;
      for (j = 0, sep = ffp->seplist; j < (Int2) val && sep != NULL; j++, sep = sep->next) continue;
      if (sep != NULL) {
        if (collision [j] != NULL) {
          okay = FALSE;
        }
        collision [j] = sep;
        order [i] = sep;
      }
    }
    MemFree (str);
  }
  if (! okay) {
    Message (MSG_ERROR, "You must not select a contig more than once");
    MemFree (order);
    MemFree (collision);
    return;
  }

  okay = FALSE;
  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      okay = TRUE;
    }
  }
  /* if no contigs selected, use all in order */
  if (! okay) {
    for (j = 0, sep = ffp->seplist; j < count && sep != NULL; j++, sep = sep->next) {
      order [j] = sep;
    }
    okay = TRUE;
  }
  /*
  if (! okay) {
    Message (MSG_ERROR, "You must select at least one contig");
    MemFree (order);
    MemFree (collision);
    return;
  }
  */

  /* use spreadsheet to reorder ffp->seplist, delete unwanted items */

  /*
  contigs [0] = '\0';
  notfirst = FALSE;
  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      sep = order [i];
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      if (notfirst) {
        StringCat (contigs, " ");
      }
      StringCat (contigs, tmp);
      notfirst = TRUE;
    }
  }
  ffp->seplist = SetPhrapContigOrder (ffp->seplist, contigs);
  */

  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      sep = ffp->seplist;
      lastsep = NULL;
      while (sep != NULL && sep != order [i]) {
        lastsep = sep;
        sep = sep->next;
      }
      if (sep != NULL) {
        if (lastsep != NULL) {
          lastsep->next = sep->next;
          sep->next = NULL;
        } else {
          ffp->seplist = sep->next;
          sep->next = NULL;
        }
      }
    }
  }
  sep = ffp->seplist;
  while (sep != NULL) {
    nextsep = sep->next;
    sep->next = NULL;
    SeqEntryFree (sep);
    sep = nextsep;
  }
  ffp->seplist = NULL;

  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      ValNodeLink (&(ffp->seplist), order [i]);
    }
  }

  MemFree (order);
  MemFree (collision);
  SafeDisable (b);
  SafeHide (ffp->orderblock);
  SafeShow (ffp->controlblock);
  if (TextHasNoText (ffp->orgname)) {
    Select (ffp->orgname);
  } else {
    Select (ffp->seqname);
  }
  if (Fa2htgsFormOkay (ffp)) {
    SafeEnable (ffp->okBtn);
  } else {
    SafeDisable (ffp->okBtn);
  }
}

static void ReadContigFile (ButtoN b)

{
  BioseqPtr       bsp;
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  Boolean         onMaster;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep = NULL;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    WatchCursor ();
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      onMaster = (Boolean) (GetValue (ffp->contigtype) == 2);
      sep = ReadContigList (fp, onMaster);
      if (sep != NULL && sep->choice == 1) {
        ffp->seplist = sep;
        bsp = (BioseqPtr) sep->data.ptrvalue;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
    }
    FileClose (fp);
    ArrowCursor ();
    if (ffp->seplist != NULL) {
      SafeShow (ffp->controlblock);
      SafeDisable (b);
      SafeDisable (ffp->contigtype);
      if (TextHasNoText (ffp->orgname)) {
        Select (ffp->orgname);
      } else {
        Select (ffp->seqname);
      }
      if (Fa2htgsFormOkay (ffp)) {
        SafeEnable (ffp->okBtn);
      } else {
        SafeDisable (ffp->okBtn);
      }
    }
  }
}

static void AcceptFa2htgs (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ffp->form);
  Update ();
  if (ffp->finish != NULL) {
    ffp->finish (b);
  }
  SeqSubmitFree (ffp->ssp);
  SeqEntryFree (ffp->sep);
  Remove (ffp->form);
}

static void CancelFa2htgs (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ffp->form);
  Update ();
  if (ffp->cancel != NULL) {
    ffp->cancel (b);
  }
  SeqSubmitFree (ffp->ssp);
  SeqEntryFree (ffp->sep);
  Remove (ffp->form);
}

static void SetFa2htgsAcceptBtn (Handle control)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (control);
  if (ffp == NULL) return;
  if (Fa2htgsFormOkay (ffp)) {
    SafeEnable (ffp->okBtn);
  } else {
    SafeDisable (ffp->okBtn);
  }
}

static void SetFa2htgsUpdate (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (GetStatus (b)) {
    SafeEnable (ffp->accession);
  } else {
    SafeDisable (ffp->accession);
  }
  SetFa2htgsAcceptBtn ((Handle) b);
}

static void CleanupGenomeCenterForm (GraphiC g, VoidPtr data)

{
  Fa2htgsFormPtr  ffp;
  SeqEntryPtr     next;
  SeqEntryPtr     sep;

  ffp = (Fa2htgsFormPtr) data;
  if (ffp != NULL) {
    sep = ffp->seplist;
    while (sep != NULL) {
      next = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
      sep = next;
    }
    if (ffp->alists [0] != NULL) {
      FreeEnumFieldAlist (ffp->alists [0]);
    }
  }
  StdCleanupFormProc (g, data);
}

static CharPtr secaccstrings [] = {"Secondary", "Accessions", NULL};

static Uint2 contigorder_types [] = {
  TAGLIST_POPUP
};

static ENUM_ALIST(contigdefault_alist)
  {"                              ", 0},
END_ENUM_ALIST

static EnumFieldAssocPtr contigdefault_alists [] = {
  contigdefault_alist
};

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

extern ForM CreateGenomeCenterForm (Int2 left, Int2 top, CharPtr title,
                                    BtnActnProc finish,
                                    BtnActnProc cancel,
                                    Boolean readPhrap,
                                    Boolean buildContig,
                                    WndActnProc activateForm)

{
  ButtoN             b;
  GrouP              c;
  Fa2htgsFormPtr     ffp;
  GrouP              g;
  GrouP              h;
  PrompT             ppt;
  GrouP              q;
  GrouP              sa;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  Int2               wid;
  GrouP              x;
  GrouP              y;

  ffp = (Fa2htgsFormPtr) MemNew (sizeof (Fa2htgsForm));
  if (ffp == NULL) return NULL;
  ffp->finish = finish;
  ffp->cancel = cancel;
  ffp->readPhrap = readPhrap;
  ffp->buildContig = buildContig;

  w = FixedWindow (left, top, -10, -10, title, NULL);
  if (w == NULL) return NULL;
  SetObjectExtra (w, ffp, CleanupGenomeCenterForm);

  ffp->form = (ForM) w;
  ffp->toform = NULL;
  ffp->fromform = Fa2htgsToSeqSubmitPtr;
  ffp->formmessage = Fa2htgsFormMessage;

  ffp->ssp = NULL;
  ffp->sep = NULL;
  ffp->seplist = NULL;
  ffp->alists [0] = NULL;

#ifndef WIN_MAC
  CreateFa2htgsFormMenus (w);
#endif

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    ffp->appmessage = sepp->handleMessages;
  }

  SetGroupSpacing (w, 10, 10);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ffp->templateblock = HiddenGroup (h, -1, 0, NULL);
  b = PushButton (ffp->templateblock, "Read Seq-submit Template", ReadTemplate);
  SetObjectExtra (b, ffp, NULL);

  ffp->fastablock = HiddenGroup (h, -1, 0, NULL);
  if (readPhrap) {
    b = PushButton (ffp->fastablock, "Read PHRAP File", ReadAPhrapFile);
  } else if (buildContig) {
    SetGroupSpacing (ffp->fastablock, 10, 10);
    y = HiddenGroup (ffp->fastablock, 1, 0, NULL);
    StaticPrompt (y, "Coordinates are on", 0, stdLineHeight, programFont, 'c');
    ffp->contigtype = HiddenGroup (y, 3, 0, NULL);
    RadioButton (ffp->contigtype, "Individual Accessions");
    RadioButton (ffp->contigtype, "Master Sequence");
    SetValue (ffp->contigtype, 1);
    b = PushButton (ffp->fastablock, "Read CONTIG Instructions", ReadContigFile);
    AlignObjects (ALIGN_CENTER, (HANDLE) y, (HANDLE) b, NULL);
  } else {
    b = PushButton (ffp->fastablock, "Read FASTA File", ReadFastaHtgsFile);
  }
  SetObjectExtra (b, ffp, NULL);
  Hide (ffp->fastablock);

  x = HiddenGroup (h, 0, 0, NULL);

  ffp->orderblock = HiddenGroup (x, -1, 0, NULL);
  SetGroupSpacing (ffp->orderblock, 5, 5);
  ppt = StaticPrompt (ffp->orderblock, "Enter contigs in desired order", 0, 0, programFont, 'c');
  ffp->contigorder = CreateTagListDialogEx (ffp->orderblock, 8, 1, 2,
                                            contigorder_types, NULL,
                                            contigdefault_alists,
                                            TRUE, TRUE, NULL, NULL);
  b = PushButton (ffp->orderblock, "Proceed", PhrapOrderChosen);
  SetObjectExtra (b, ffp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) ffp->contigorder,
                (HANDLE) b, NULL);
  Hide (ffp->orderblock);

  ffp->controlblock = HiddenGroup (x, -1, 0, NULL);
  g = HiddenGroup (ffp->controlblock, 2, 0, NULL);

  ffp->htgsphase = NULL;
  if (! buildContig) {
    StaticPrompt (g, "Phase", 0, stdLineHeight, programFont, 'l');
    ffp->htgsphase = HiddenGroup (g, 4, 0, (GrpActnProc) SetFa2htgsAcceptBtn);
    SetObjectExtra (ffp->htgsphase, ffp, NULL);
    RadioButton (ffp->htgsphase, "HTGS-0");
    RadioButton (ffp->htgsphase, "HTGS-1");
    RadioButton (ffp->htgsphase, "HTGS-2");
    RadioButton (ffp->htgsphase, "HTGS-3");
    StaticPrompt (g, "HTGS_", 0, stdLineHeight, programFont, 'l');
    q = HiddenGroup (g, 3, 0, NULL);
    ffp->draft = CheckBox (q, "Draft", NULL);
    ffp->fulltop = CheckBox (q, "Fulltop", NULL);
    ffp->activefin = CheckBox (q, "Activefin", NULL);
  }

  StaticPrompt (g, "Organism", 0, dialogTextHeight, programFont, 'l');
  ffp->orgname = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->orgname, ffp, NULL);

  StaticPrompt (g, "Sequence name", 0, dialogTextHeight, programFont, 'l');
  ffp->seqname = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->seqname, ffp, NULL);

  StaticPrompt (g, "Length", 0, dialogTextHeight, programFont, 'l');
  ffp->knownlength = DialogText (g, "", 10, NULL);
  SetObjectExtra (ffp->knownlength, ffp, NULL);

  StaticPrompt (g, "Gap Length", 0, dialogTextHeight, programFont, 'l');
  ffp->gaplength = DialogText (g, "100", 10, NULL);
  SetObjectExtra (ffp->gaplength, ffp, NULL);

  ffp->update = CheckBox (g, "Update", SetFa2htgsUpdate);
  SetObjectExtra (ffp->update, ffp, NULL);
  ffp->accession = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->accession, ffp, NULL);
  Disable (ffp->accession);

  StaticPrompt (g, "Chromosome", 0, dialogTextHeight, programFont, 'l');
  ffp->chromosome = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->chromosome, ffp, NULL);

  StaticPrompt (g, "Clone", 0, dialogTextHeight, programFont, 'l');
  ffp->clone = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->clone, ffp, NULL);

  StaticPrompt (g, "Strain", 0, dialogTextHeight, programFont, 'l');
  ffp->strain = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->strain, ffp, NULL);

  StaticPrompt (g, "Cultivar", 0, dialogTextHeight, programFont, 'l');
  ffp->cultivar = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->cultivar, ffp, NULL);

  wid = MaxStringWidths (secaccstrings) + 2;
  sa = MultiLinePrompt (g, "Secondary Accessions", wid, programFont);
  ffp->secondaries = CreateVisibleStringDialog (g, 3, -1, 15);
  AlignObjects (ALIGN_MIDDLE, (HANDLE) sa, (HANDLE) ffp->secondaries, NULL);

  q = HiddenGroup (ffp->controlblock, 1, 0, NULL);

  StaticPrompt (q, "Title", 0, 0, programFont, 'c');
  ffp->title = ScrollText (q, 20, 4, programFont, TRUE, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->title, ffp, NULL);

  StaticPrompt (q, "Remark", 0, 0, programFont, 'c');
  ffp->remark = ScrollText (q, 20, 4, programFont, TRUE, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->remark, ffp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) q, NULL);
  Hide (ffp->controlblock);

  c = HiddenGroup (h, 2, 0, NULL);
  ffp->okBtn = DefaultButton (c, "Accept", AcceptFa2htgs);
  SetObjectExtra (ffp->okBtn, ffp, NULL);
  Disable (ffp->okBtn);
  b = PushButton (c, "Cancel", CancelFa2htgs);
  SetObjectExtra (b, ffp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) ffp->templateblock, (HANDLE) ffp->fastablock,
                (HANDLE) ffp->orderblock, (HANDLE) ffp->controlblock, (HANDLE) c, NULL);

  RealizeWindow (w);

  if (activateForm != NULL) {
    SetActivate (w, activateForm);
  }

  Show (w);
  Select (w);
  Select (ffp->orgname);

  return (ForM) w;
}

typedef struct convcdsdata {
  FEATURE_FORM_BLOCK

  SeqEntryPtr    sep;
  TexT           geneName;
  TexT           protName;
  TexT           featcomment;
  ButtoN         retain;
  Uint2          subtype;
  Int2           errcount;
  Char           findThis [128];
  SeqLocPtr      slp;
  SeqEntryPtr    nsep;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
} ConvCdsData, PNTR ConvCdsPtr;

static Boolean CollectCDSAncestorGatherFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  ConvCdsPtr     ccp;
  Boolean        noLeft;
  Boolean        noRight;
  ObjMgrTypePtr  omtp;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Uint2          subtype;

  if (gcp == NULL) return TRUE;

  ccp = (ConvCdsPtr) gcp->userdata;
  if (ccp == NULL ) return TRUE;

  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  omtp = ccp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return TRUE;

  sfp = (SeqFeatPtr) gcp->thisitem;
  subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
  if (subtype != ccp->subtype) return TRUE;

  bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  if (bsp == NULL) return TRUE;
  if (ISA_aa (bsp->mol)) return TRUE;
  if (bsp->repr != Seq_repr_seg) {
    sep = ccp->nsep;
    if (sep == NULL || sep->choice != 1) return TRUE;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return TRUE;
  }
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  slp = SeqLocMerge (bsp, sfp->location, ccp->slp, FALSE, TRUE, FALSE);
  if (slp == NULL) return TRUE;
  SetSeqLocPartial (slp, noLeft, noRight);

  ccp->slp = SeqLocFree (ccp->slp);
  ccp->slp = slp;

  return TRUE;
}

static void FinishConvertingToCDS (SeqEntryPtr sep, Uint2 entityID, ConvCdsPtr ccp)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  CdRegionPtr   crp;
  ValNodePtr    descr;
  Uint1         frame;
  Int2          genCode;
  Int2          i;
  Int4          len;
  Int4          lens [4];
  Int4          max;
  Boolean       noLeft;
  Boolean       noRight;
  MolInfoPtr    mip;
  SeqEntryPtr   nsep;
  SeqEntryPtr   old;
  CharPtr       prot;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  CharPtr       ptr;
  SeqFeatPtr    sfp;
  Char          str [128];
  ValNodePtr    vnp;

  if (sep == NULL || ccp == NULL) return;
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (crp == NULL) return;
  sfp = CreateNewFeature (ccp->nsep, NULL, SEQFEAT_CDREGION, NULL);
  if (sfp == NULL) {
    CdRegionFree (crp);
    return;
  }
  sfp->data.value.ptrvalue = (Pointer) crp;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = AsnIoMemCopy ((Pointer) ccp->slp,
                                (AsnReadFunc) SeqLocAsnRead,
                                (AsnWriteFunc) SeqLocAsnWrite);
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  sfp->partial = (sfp->partial || noLeft || noRight);
  if (! TextHasNoText (ccp->featcomment)) {
    sfp->comment = SaveStringFromTextAndStripNewlines (ccp->featcomment);
  }
  max = 0;
  frame = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = (Uint1) i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    len = BSLen (bs);
    BSFree (bs);
    lens [i] = len;
    if (len > max) {
      max = len;
      frame = (Uint1) i;
    }
  }
  for (i = 1; i <= 3; i++) {
    if (lens [i] == max && i != frame) {
      (ccp->errcount)++;
    }
  }
  crp->frame = frame;
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs == NULL) return;
  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot == NULL) return;
  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = (Int2) StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    /*
    if (prot [0] == '-') {
      ptr++;
    }
    */
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }
  MemFree (prot);
  if (bs == NULL) return;
  bsp = BioseqNew ();
  if (bsp == NULL) return;
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (sep);
  bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);
  psep = SeqEntryNew ();
  if (psep != NULL) {
    psep->choice = 1;
    psep->data.ptrvalue = (Pointer) bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 8;
      if (noLeft && noRight) {
        mip->completeness = 5;
      } else if (noLeft) {
        mip->completeness = 3;
      } else if (noRight) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    descr = ExtractBioSourceAndPubs (sep);
    /*
    AddSeqEntryToSeqEntry (sep, psep, FALSE);
    */
    AddSeqEntryToSeqEntry (sep, psep, TRUE);
    nsep = FindNucSeqEntry (sep);
    ReplaceBioSourceAndPubs (sep, descr);
    SetSeqFeatProduct (sfp, bsp);
    GetTitle (ccp->protName, str, sizeof (str));
    if (! StringHasNoText (str)) {
      prp = CreateNewProtRef (str, NULL, NULL, NULL);
      if (prp != NULL) {
        sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = (Pointer) prp;
          SetSeqLocPartial (sfp->location, noLeft, noRight);
          sfp->partial = (sfp->partial || noLeft || noRight);
        }
      }
    }
  }
}

static void RemoveCDSAncestorCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ConvCdsPtr     ccp;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  ObjMgrTypePtr  omtp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  Uint2          subtype;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  ccp = (ConvCdsPtr) mydata;
  if (ccp == NULL) return;
  omtp = ccp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (subtype == ccp->subtype) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          SeqFeatFree (sfp);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void ConvertToCDSCallback (SeqEntryPtr sep, Uint2 entityID, ConvCdsPtr ccp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  GeneRefPtr    grp;
  GatherScope   gs;
  SeqLocPtr     gslp;
  Boolean       hasNulls;
  Boolean       noLeft;
  Boolean       noRight;
  SeqEntryPtr   nsep;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Char          str [128];

  if (sep == NULL || ccp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ConvertToCDSCallback (sep, entityID, ccp);
      }
      return;
    }
  }
  ccp->nsep = FindNucSeqEntry (sep);
  if (ccp->nsep == NULL) return;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  ccp->slp = NULL;
  GatherEntity (entityID, (Pointer) ccp, CollectCDSAncestorGatherFunc, &gs);
  if (ccp->slp != NULL) {
    CheckSeqLocForPartial (ccp->slp, &noLeft, &noRight);
    sip = SeqLocId (ccp->slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL && ISA_na (bsp->mol)) {
        slp = SegLocToParts (bsp, ccp->slp);
        if (slp != NULL) {
          ccp->slp = SeqLocFree (ccp->slp);
          ccp->slp = slp;
          FreeAllFuzz (ccp->slp);
          SetSeqLocPartial (ccp->slp, noLeft, noRight);
        }
      }
    }
    FinishConvertingToCDS (sep, entityID, ccp);
    nsep = FindNucSeqEntry (sep);
    GetTitle (ccp->geneName, str, sizeof (str));
    if (! StringHasNoText (str)) {
      grp = CreateNewGeneRef (str, NULL, NULL, FALSE);
      if (grp != NULL) {
        if (ExtendGene (grp, nsep, ccp->slp)) {
          grp = GeneRefFree (grp);
        } else {
          sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
          if (sfp != NULL) {
            sfp->data.value.ptrvalue = (Pointer) grp;
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = AsnIoMemCopy ((Pointer) ccp->slp,
                                          (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
            sip = SeqLocId (sfp->location);
            if (sip != NULL) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) {
                gslp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
                if (gslp != NULL) {
                  sfp->location = SeqLocFree (sfp->location);
                  sfp->location = gslp;
                  if (bsp->repr == Seq_repr_seg) {
                    gslp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                    sfp->location = SeqLocFree (sfp->location);
                    sfp->location = gslp;
                    hasNulls = LocationHasNullsBetween (sfp->location);
                    sfp->partial = (sfp->partial || hasNulls);
                  }
                  FreeAllFuzz (gslp);
                  SetSeqLocPartial (sfp->location, noLeft, noRight);
                  sfp->partial = (sfp->partial || noLeft || noRight);
                }
              }
            }
          }
        }
      }
    }
  }
  ccp->slp = SeqLocFree (ccp->slp);
}

static void DoConvertToCDS (ButtoN b)

{
  ConvCdsPtr  ccp;
  CharPtr     plural;
  Char        str [128];

  ccp = GetObjectExtra (b);
  if (ccp == NULL) {
    Remove (ParentWindow (b));
    return;
  }
  GetTitle (ccp->protName, str, sizeof (str));
  if (StringHasNoText (str)) {
    Message (MSG_OK, "Protein name is required");
    return;
  }
  Hide (ccp->form);
  WatchCursor ();

  ccp->omp = ObjMgrGet ();
  if (ccp->omp != NULL) {
    ccp->omtp = ObjMgrTypeFind (ccp->omp, OBJ_SEQFEAT, NULL, NULL);
    if (ccp->omtp != NULL && ccp->omtp->subtypefunc != NULL) {
      ConvertToCDSCallback (ccp->sep, ccp->input_entityID, ccp);
      if (! GetStatus (ccp->retain)) {
        SeqEntryExplore (ccp->sep, (Pointer) ccp, RemoveCDSAncestorCallback);
      }
    }
  }

  ArrowCursor ();
  Update ();

  if (ccp->errcount > 0) {
    if (ccp->errcount > 1) {
      plural = "records";
    } else {
      plural = "record";
    }
    Message (MSG_ERROR, "Possible ambiguous frames detected in %d %s",
             (int) ccp->errcount, plural);
  }

  ObjMgrSetDirtyFlag (ccp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ccp->input_entityID, 0, 0);
  Remove (ccp->form);
  Update ();
}

extern void PrepareToConvertToCDS (SeqEntryPtr sep, Uint2 entityID,
                                   Uint2 subtype, CharPtr findthis)

{
  ButtoN             b;
  GrouP              c;
  ConvCdsPtr         ccp;
  GrouP              g;
  GrouP              h;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  if (sep == NULL || entityID == 0 || subtype == 0) return;

  ccp = (ConvCdsPtr) MemNew (sizeof (ConvCdsData));
  if (ccp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Convert to CDS", StdCloseWindowProc);
  SetObjectExtra (w, ccp, StdCleanupFormProc);
  ccp->form = (ForM) w;
  ccp->formmessage = DefaultMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ccp->appmessage = sepp->handleMessages;
  }

  ccp->input_entityID = entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ccp->sep = sep;
  ccp->subtype = subtype;
  StringNCpy_0 (ccp->findThis, findthis, sizeof (ccp->findThis));
  ccp->errcount = 0;

  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
  ccp->geneName = DialogText (g, "", 20, NULL);
  StaticPrompt (g, "Protein Name", 0, dialogTextHeight, programFont, 'l');
  ccp->protName = DialogText (g, "", 20, NULL);
  StaticPrompt (g, "Comment", 0, 4 * Nlm_stdLineHeight, programFont, 'l');
  ccp->featcomment = ScrollText (g, 20, 4, programFont, TRUE, NULL);

  ccp->retain = CheckBox (h, "Retain original features", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoConvertToCDS);
  SetObjectExtra (b, ccp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) ccp->retain, NULL);
  RealizeWindow (w);
  Show (w);
  Select (ccp->geneName);
  Update ();
}

typedef struct alignform {
  FORM_MESSAGE_BLOCK
  DoC            doc;
} AlignForm, PNTR AlignFormPtr;

static ParData bioseqParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData bioseqColFmt [] = {
  {0, 5, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
};

static ParData annotParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData annotColFmt [] = {
  {0, 10, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}
};

static ParData alignParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData alignColFmt [] = {
  {0, 15, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 5, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}
};

static Boolean AlignFormPopulateProc (GatherContextPtr gcp)

{
  AlignFormPtr      afp;
  SeqAlignPtr       align;
  Char              annotDB [32];
  Uint1             annot_type;
  BioseqContextPtr  bcp;
  BioseqPtr         bsp;
  DenseDiagPtr      ddp;
  DenseSegPtr       dsp;
  ProtRefPtr        prp;
  CharPtr           ptr;
  SeqAnnotPtr       sap;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;
  StdSegPtr         ssp;
  Char              str [256];
  Char              tmp [128];

  if (gcp == NULL) return TRUE;

  afp = (AlignFormPtr) gcp->userdata;
  if (afp == NULL ) return TRUE;

  switch (gcp->thistype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) gcp->thisitem;
      if (bsp != NULL) {
        SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
        bcp = BioseqContextNew (bsp);
        sfp = BioseqContextGetSeqFeat (bcp, SEQFEAT_PROT, NULL, NULL, 0);
        BioseqContextFree (bcp);
        if (sfp != NULL) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp != NULL) {
            if (prp->name != NULL && (! StringHasNoText (prp->name->data.ptrvalue))) {
              StringCat (str, " (");
              StringCat (str, (CharPtr) prp->name->data.ptrvalue);
              StringCat (str, ")");
            } else if (! StringHasNoText (prp->desc)) {
              StringCat (str, " (");
              StringCat (str, (CharPtr) prp->desc);
              StringCat (str, ")");
            }
          }
        }
        StringCat (str, ":\n");
        AppendText (afp->doc, str, &bioseqParFmt, bioseqColFmt, programFont);
      }
      break;
    case OBJ_SEQANNOT :
      sap = (SeqAnnotPtr) gcp->thisitem;
      if (sap != NULL && sap->type == 2) {
        get_align_annot_qual (sap, annotDB, sizeof (annotDB), &annot_type);
        if (annot_type == ANNOT_BLAST) {
          StringCpy (str, annotDB);
          StringCat (str, ":\n");
          AppendText (afp->doc, str, &annotParFmt, annotColFmt, programFont);
        }
      }
      break;
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      align = (SeqAlignPtr) gcp->thisitem;
      sip = NULL;
      if (align->segtype == 1) {
        ddp = (DenseDiagPtr) align->segs;
        if (ddp != NULL) {
          for (sip = ddp->id; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      } else if (align->segtype == 2) {
        dsp = (DenseSegPtr) align->segs;
        if (dsp != NULL) {
          for (sip = dsp->ids; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      } else if (align->segtype == 3) {
        ssp = (StdSegPtr) align->segs;
        if (ssp != NULL) {
          for (sip = ssp->ids; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      }
      if (sip != NULL) {
        bsp = BioseqLockById (sip);
        if (bsp != NULL) {
          /* SeqIdWrite (bsp->id, str, PRINTID_FASTA_LONG, sizeof (str)); */
          SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
          StringNCpy_0 (tmp, BioseqGetTitle (bsp), sizeof (tmp));
          StringCat (str, "\t");
          if (! StringHasNoText (tmp)) {
            StringCat (str, tmp);
          }
          tmp [0] = '\0';
          sep = SeqMgrGetSeqEntryForData (bsp);
          if (sep != NULL) {
            SeqEntryToBioSource (sep, NULL, tmp, sizeof (tmp), NULL);
            if (! StringHasNoText (tmp)) {
              ptr = StringStr (tmp, "(");
              if (ptr != NULL) {
                *ptr = '\0';
              }
              StringCat (str, " [");
              StringCat (str, tmp);
              StringCat (str, "]");
            }
          }
          StringCat (str, "\t");
          sprintf (tmp, "%d %d %d", (int) gcp->entityID,
                   (int) gcp->itemID, (int) gcp->thistype);
          StringCat (str, tmp);
          StringCat (str, "\n");
          AppendText (afp->doc, str, &alignParFmt, alignColFmt, programFont);
        }
        BioseqUnlock (bsp);
      }
      return TRUE;
    default :
      break;
  }
  return TRUE;
}

static void PopulateAlignForm (ForM f, SeqEntryPtr sep)

{
  AlignFormPtr  afp;
  GatherScope   gs;
  RecT          r;
  Int2          width;

  if (f == NULL || sep == NULL) return;
  afp = (AlignFormPtr) GetObjectExtra (f);
  if (afp == NULL) return;
  Reset (afp->doc);
  SetDocAutoAdjust (afp->doc, FALSE);
  ObjectRect (afp->doc, &r);
  InsetRect (&r, 4, 4);
  width = r.right - r.left;
  alignColFmt [0].pixWidth = width / 5;
  alignColFmt [1].pixWidth = width - alignColFmt [0].pixWidth;
  bioseqColFmt [0].pixWidth = width;
  annotColFmt [0].pixWidth = width;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQHIST] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.scope = sep;
  GatherEntity (afp->input_entityID, (Pointer) afp, AlignFormPopulateProc, &gs);
  AdjustDocScroll (afp->doc);
}

static void AlignMessageProc (ForM f, Int2 mssg)

{
  FILE          *fp;
  AlignFormPtr  afp;
  Char          path [PATH_MAX];

  afp = (AlignFormPtr) GetObjectExtra (f);
  if (afp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_PRINT :
        PrintDocument (afp->doc);
        break;
      case VIB_MSG_EXPORT :
        if (GetOutputFileName (path, sizeof (path), "align.txt")) {
          WatchCursor ();
#ifdef WIN_MAC
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            FileClose (fp);
          } else {
            FileCreate (path, "TEXT", "ttxt");
          }
#endif
          fp = FileOpen (path, "w");
          if (fp != NULL) {
            SaveDocument (afp->doc, fp);
            FileClose (fp);
          }
          ArrowCursor ();
        }
        break;
      default :
        if (afp->appmessage != NULL) {
          afp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void AlignFormActivate (WindoW w)

{
  AlignFormPtr  afp;
  IteM          exportItm;
  IteM          printItm;

  afp = (AlignFormPtr) GetObjectExtra (w);
#ifdef WIN_MAC
  currentFormDataPtr = (VoidPtr) afp;
#endif
  if (afp != NULL) {
    if (afp->activate != NULL) {
      afp->activate (w);
    }
    exportItm = FindFormMenuItem ((BaseFormPtr) afp, VIB_MSG_EXPORT);
    SafeSetTitle (exportItm, "Export Align Summary...");
    SafeEnable (exportItm);
    printItm = FindFormMenuItem ((BaseFormPtr) afp, VIB_MSG_PRINT);
    SafeEnable (printItm);
  }
}

static void AlignNotifyProc (DoC d, Int2 item, Int2 row, Int2 col, Boolean dblClick)

{
  unsigned int  entityID;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       txt;

  if (d == NULL || item < 1 || row < 1 || col < 1) return;
  txt = GetDocText (d, item, row, 3);
  if (! StringHasNoText (txt)) {
    if (sscanf (txt, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
      if (shftKey) {
        ObjMgrAlsoSelect (entityID, itemID, itemtype, 0, NULL);
      } else {
        ObjMgrSelect (entityID, itemID, itemtype, 0, NULL);
      }
    }
  }
  MemFree (txt);
}

static ForM CreateAlignForm (Uint2 entityID)

{
  AlignFormPtr       afp;
  GrouP              c;
  GrouP              h;
  StdEditorProcsPtr  sepp;
  WindoW             w;
#ifndef WIN_MAC
  MenU               m;
#endif

  w = NULL;
  afp = MemNew (sizeof (AlignForm));
  if (afp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Alignment Summary", StdCloseWindowProc);
    SetObjectExtra (w, afp, StdCleanupFormProc);
    afp->form = (ForM) w;
    afp->formmessage = AlignMessageProc;
    afp->input_entityID = entityID;

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      afp->appmessage = sepp->handleMessages;
    }

#ifndef WIN_MAC
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Close", (BaseFormPtr) afp, VIB_MSG_CLOSE);
    SeparatorItem (m);
    FormCommandItem (m, "Export...", (BaseFormPtr) afp, VIB_MSG_EXPORT);
    SeparatorItem (m);
    FormCommandItem (m, "Print...", (BaseFormPtr) afp, VIB_MSG_PRINT);
#endif

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    afp->doc = DocumentPanel (h, 35 * stdCharWidth, 20 * stdLineHeight);
    SetObjectExtra (afp->doc, afp, NULL);
    SetDocNotify (afp->doc, AlignNotifyProc);

    c = HiddenGroup (w, 4, 0, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) afp->doc, (HANDLE) c, NULL);

    RealizeWindow (w);
    SetActivate (w, AlignFormActivate);
  }
  return (ForM) w;
}

extern void ViewAlignmentSummary (IteM i)

{
  BaseFormPtr  bfp;
  ForM         f;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  WatchCursor ();
  Update ();
  f = CreateAlignForm (bfp->input_entityID);
  PopulateAlignForm (f, sep);
  ArrowCursor ();
  Update ();
  Show (f);
  Select (f);
}

typedef struct geneprotfind {
  SeqFeatPtr     cds;
  SeqFeatPtr     gene;
  SeqFeatPtr     prot;
  SeqLocPtr      location;
  SeqLocPtr      product;
  Int4           mingene;
  Int4           minprot;
} GeneProtFind, PNTR GeneProtPtr;

static Boolean GeneProtFindFunc (GatherContextPtr gcp)

{
  Int4         diff;
  GeneProtPtr  gpp;
  SeqFeatPtr   sfp;

  if (gcp == NULL) return TRUE;
  gpp = (GeneProtPtr) gcp->userdata;
  if (gpp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.value.ptrvalue == NULL) return TRUE;
  if (sfp->data.choice == SEQFEAT_GENE && gpp->location != NULL) {
    diff = SeqLocAinB (gpp->location, sfp->location);
    if (diff >= 0) {
      if (diff < gpp->mingene) {
        gpp->gene = sfp;
        gpp->mingene = diff;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_PROT && gpp->product != NULL) {
    diff = SeqLocAinB (gpp->product, sfp->location);
    if (diff >= 0) {
      if (diff < gpp->minprot) {
        gpp->prot = sfp;
        gpp->minprot = diff;
      }
    }
  }
  return TRUE;
}

void FindGeneAndProtForCDS (Uint2 entityID, SeqFeatPtr cds,
                            SeqFeatPtr PNTR gene, SeqFeatPtr PNTR prot)

{
  GeneProtFind  gpf;
  GatherScope   gs;

  if (gene != NULL) {
    *gene = NULL;
  }
  if (prot != NULL) {
    *prot = NULL;
  }
  if (entityID == 0 || cds == NULL) return;
  MemSet ((Pointer) (&gpf), 0, sizeof (GeneProtFind));
  gpf.cds = cds;
  gpf.gene = NULL;
  gpf.prot = NULL;
  gpf.location = cds->location;
  gpf.product = cds->product;
  gpf.mingene = INT4_MAX;
  gpf.minprot = INT4_MAX;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (entityID, (Pointer) &gpf, GeneProtFindFunc, &gs);
  if (gene != NULL) {
    *gene = gpf.gene;
  }
  if (prot != NULL) {
    *prot = gpf.prot;
  }
}


#define FIND_ASN   1
#define FIND_FLAT  2
#define FIND_GENE  3
#define FIND_PROT  4
#define FIND_POS   5

typedef struct findform {
  FORM_MESSAGE_BLOCK
  TexT            findTxt;
  TexT            replaceTxt;
  ButtoN          caseCounts;
  ButtoN          wholeWord;
  ButtoN          doSeqIdLocal;
  ButtoN          findAllBtn;
  ButtoN          replaceAllBtn;
  Int2            type;
  ButtoN          findFirstBtn;
  ButtoN          findNextBtn;
  Int4            last_paragraph_found;
} FindForm, PNTR FindFormPtr;

extern CharPtr CompressSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  Char     last;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      last = ch;
      ch = *ptr;
      if (ch != '\0' && ch < ' ') {
        *ptr = ' ';
        ch = *ptr;
      }
      while (ch != '\0' && last <= ' ' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

extern CharPtr SearchForString (CharPtr str, CharPtr sub, Boolean case_counts, Boolean whole_word)

{
  CharPtr  ptr = NULL;
  CharPtr  tmp;

  if (case_counts) {
    ptr = StringSearch (str, sub);
  } else {
    ptr = StringISearch (str, sub);
  }
  if (ptr == NULL) return NULL;
  if (whole_word) {
    if (ptr > str) {
      tmp = ptr - 1;
      if (! IS_WHITESP (*tmp)) return NULL;
    }
    tmp = ptr + StringLen (sub);
    if (*tmp != '\0' && (! IS_WHITESP (*tmp))) return NULL;
  }
  return ptr;
}


static Int4
FindInFlatFileNext 
(Asn2gbJobPtr ajp,
 CharPtr sub,
 Boolean case_counts,
 Boolean whole_word,
 Boolean stop_after_one_found,
 Int4    start_index)
{
  Int4             index;
  CharPtr          string;
  BaseBlockPtr     bbp;
  SelStructPtr     sel;
  unsigned int     entID;
  unsigned int     itmID;
  unsigned int     itmtype;
  Boolean          already = FALSE;

  if (ajp == NULL) return -1;
      
  for (index = start_index + 1; index < ajp->numParagraphs; index++) {
    string = asn2gnbk_format (ajp, (Int4) index);
    if (string != NULL && *string != '\0') {
      CompressSpaces (string);
      if (SearchForString (string, sub, case_counts, whole_word) != NULL) {
        bbp = ajp->paragraphArray [index];
        if (bbp != NULL) {
          entID = bbp->entityID;
          itmID = bbp->itemID;
          itmtype = bbp->itemtype;
          already = FALSE;
          for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
            if (sel->entityID == entID && sel->itemID == itmID && sel->itemtype == itmtype) {
              already = TRUE;
            }
          }
          if (! already) {
            ObjMgrAlsoSelect (entID, itmID, itmtype, 0, NULL);
            if (stop_after_one_found) return index;
          }
        }
      }
    }
  }
  return -1;
}


static Int4 
FindInFlatFile 
(Uint2 entityID,
 Uint4 itemID,
 Uint2 itemtype,                           
 Uint1 format,
 Uint1 mode,                           
 CharPtr sub,
 Boolean case_counts,
 Boolean whole_word,                           
 Boolean stop_after_one,
 Int4 start_index)

{
  Asn2gbJobPtr     ajp;
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  ErrSev           level;
  SeqEntryPtr      oldsep;
  SeqEntryPtr      sep = NULL;
  SeqEntryPtr      topsep;
  SeqEntryPtr      usethetop = NULL;
  Int4             retval = -1;
  FlgType          flags = SHOW_CONTIG_FEATURES | SHOW_CONTIG_SOURCES | SHOW_FAR_TRANSLATION;

  if (itemID == 0) {
    sep = GetTopSeqEntryForEntityID (entityID);
    usethetop = sep;
  } else {
    bsp = GetBioseqGivenIDs (entityID, itemID, itemtype);
    if (bsp == NULL)
      return retval;
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (bsp->repr == Seq_repr_seg) {
      sep = GetBestTopParentForData (entityID, bsp);
    }
  }

  if (sep == NULL)
    return retval;

  level = ErrSetMessageLevel (SEV_MAX);
  WatchCursor ();
  Update ();

  topsep = GetTopSeqEntryForEntityID (entityID);
  oldsep = SeqEntrySetScope (topsep);

  if (usethetop != NULL && IS_Bioseq_set (usethetop)) {
    bssp = (BioseqSetPtr) usethetop->data.ptrvalue;
    ajp = asn2gnbk_setup (NULL, bssp, NULL, (FmtType) format, (ModType) mode, NORMAL_STYLE, flags, (LckType) 0, (CstType) 0, NULL);
  } else {
    ajp = asn2gnbk_setup (bsp, NULL, NULL, (FmtType) format, (ModType) mode, NORMAL_STYLE, flags, (LckType) 0, (CstType) 0, NULL);
  }
  if (ajp != NULL) {
    retval = FindInFlatFileNext (ajp, sub, case_counts, whole_word, stop_after_one, start_index);
    asn2gnbk_cleanup (ajp);
  }

  SeqEntrySetScope (oldsep);

  ErrSetMessageLevel (level);
  ArrowCursor ();
  Update ();
  return retval;
}

static void FindFlatProc (ButtoN b)

{
  Boolean      caseCounts;
  FindFormPtr  ffp;
  CharPtr      findme;
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  findme = SaveStringFromText (ffp->findTxt);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  caseCounts = GetStatus (ffp->caseCounts);
  wholeWord = GetStatus (ffp->wholeWord);
  CompressSpaces (findme);
  ffp->last_paragraph_found = FindInFlatFile (ffp->input_entityID, ffp->input_itemID,
                  ffp->input_itemtype, GENBANK_FMT, SEQUIN_MODE,
                  findme, caseCounts, wholeWord, FALSE, -1);
  MemFree (findme);
}

static void FindFlatProcFirst (ButtoN b)

{
  Boolean      caseCounts;
  FindFormPtr  ffp;
  CharPtr      findme;
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  findme = SaveStringFromText (ffp->findTxt);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  caseCounts = GetStatus (ffp->caseCounts);
  wholeWord = GetStatus (ffp->wholeWord);
  CompressSpaces (findme);
  ffp->last_paragraph_found = FindInFlatFile (ffp->input_entityID, ffp->input_itemID,
                  ffp->input_itemtype, GENBANK_FMT, SEQUIN_MODE,
                  findme, caseCounts, wholeWord, TRUE, -1);
  MemFree (findme);
}

static void FindFlatProcNext (ButtoN b)

{
  Boolean      caseCounts;
  FindFormPtr  ffp;
  CharPtr      findme;
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  findme = SaveStringFromText (ffp->findTxt);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  caseCounts = GetStatus (ffp->caseCounts);
  wholeWord = GetStatus (ffp->wholeWord);
  CompressSpaces (findme);
  ffp->last_paragraph_found = FindInFlatFile (ffp->input_entityID, ffp->input_itemID,
                  ffp->input_itemtype, GENBANK_FMT, SEQUIN_MODE,
                  findme, caseCounts, wholeWord, TRUE, ffp->last_paragraph_found);
  MemFree (findme);
}

static void FindTextProc (TexT t)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (t);
  if (ffp != NULL) {
    if (TextLength (t) > 0) {
      SafeEnable (ffp->findAllBtn);
      SafeEnable (ffp->replaceAllBtn);
      if (ffp->type == FIND_FLAT) {
        SafeEnable (ffp->findFirstBtn);
        SafeEnable (ffp->findNextBtn);
      }
      if (ffp->type == FIND_GENE || ffp->type == FIND_PROT) {
        SafeEnable (ffp->findFirstBtn);
      }
    } else {
      SafeDisable (ffp->findAllBtn);
      SafeDisable (ffp->replaceAllBtn);
      if (ffp->type == FIND_FLAT) {
        SafeDisable (ffp->findFirstBtn);
        SafeDisable (ffp->findNextBtn);
      }
      if (ffp->type == FIND_GENE || ffp->type == FIND_PROT) {
        SafeDisable (ffp->findFirstBtn);
      }
    }
  }
}

static void CommonFindReplaceProc (ButtoN b, Boolean replace, Boolean replaceAll)

{
  Boolean      caseCounts;
  CharPtr      changeme;
  Boolean      doSeqIdLocal;
  FindFormPtr  ffp;
  CharPtr      findme;
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp != NULL) {
    findme = JustSaveStringFromText (ffp->findTxt);
    ObjMgrDeSelect (0, 0, 0, 0, NULL);
    caseCounts = GetStatus (ffp->caseCounts);
    wholeWord = GetStatus (ffp->wholeWord);
    doSeqIdLocal = GetStatus (ffp->doSeqIdLocal);
    if (replace) {
      changeme = JustSaveStringFromText (ffp->replaceTxt);
      FindReplaceInEntity (ffp->input_entityID, findme, changeme, caseCounts, wholeWord, replaceAll,
                           FALSE, UPDATE_ONCE, NULL, NULL, NULL, doSeqIdLocal, NULL, NULL);
      GetRidOfEmptyFeatsDescStrings (ffp->input_entityID, NULL);
      ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID, 0, 0);
      MemFree (changeme);
    } else {
      FindReplaceInEntity (ffp->input_entityID, findme, NULL, caseCounts, wholeWord, FALSE,
                           TRUE, UPDATE_ONCE, NULL, NULL, NULL, doSeqIdLocal, NULL, NULL);
    }
    MemFree (findme);
    Update ();
  }
}

static void FindAllProc (ButtoN b)

{
  CommonFindReplaceProc (b, FALSE, FALSE);
}

static void ReplaceAllProc (ButtoN b)

{
  CommonFindReplaceProc (b, TRUE, TRUE);
}

static SeqFeatPtr 
FindNthFeatureOnBspByLabel 
(BioseqPtr              bsp,
 CharPtr                label,
 Uint1                  seqFeatChoice,
 Uint1                  featDefChoice,
 Int4                   n,
 Int4 PNTR              last_found,
 SeqMgrFeatContext PNTR context)
{
  SMFeatItemPtr PNTR  array;
  BioseqExtraPtr      bspextra;
  Uint2               entityID;
  SMFeatItemPtr       feat;
  Int4                L;
  Int4                mid;
  Int4                num;
  ObjMgrDataPtr       omdp;
  Int4                R;
  SeqFeatPtr          sfp;
  Uint1               seqfeattype;
  Int4                index = 0;

  if (context != NULL) {
    MemSet ((Pointer) context, 0, sizeof (SeqMgrFeatContext));
  }
  if (last_found != NULL) {
    *last_found = -1;
  }

  if (bsp == NULL || StringHasNoText (label)) return NULL;

  omdp = SeqMgrGetOmdpForBioseq (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  array = bspextra->featsByLabel;
  num = bspextra->numfeats;

  if (array == NULL || num < 1) return NULL;

  if (n < 0 || n > bspextra->numfeats) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);

  /* binary search to leftmost candidate within the featsByLabel array */

  L = 0;
  R = num - 1;
  while (L < R) {
    mid = (L + R) / 2;
    feat = array [mid];
    if (feat != NULL && StringICmp (feat->label, label) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  feat = array [R];

  /* linear scan to find desired label on desired feature type */

  while (R < num && feat != NULL && StringICmp (feat->label, label) == 0) {
    sfp = feat->sfp;
    if (sfp != NULL) {
      seqfeattype = sfp->data.choice;
      if ((seqFeatChoice == 0 || seqfeattype == seqFeatChoice) &&
          (featDefChoice == 0 || feat->subtype == featDefChoice) &&
          (! feat->ignore)) {
        if (context != NULL) {
          context->entityID = entityID;
          context->itemID = feat->itemID;
          context->sfp = sfp;
          context->sap = feat->sap;
          context->bsp = feat->bsp;
          context->label = feat->label;
          context->left = feat->left;
          context->right = feat->right;
          context->dnaStop = feat->dnaStop;
          context->partialL = feat->partialL;
          context->partialR = feat->partialR;
          context->farloc = feat->farloc;
          context->strand = feat->strand;
          context->seqfeattype = seqfeattype;
          context->featdeftype = feat->subtype;
          context->numivals = feat->numivals;
          context->ivals = feat->ivals;
          context->userdata = NULL;
          context->omdp = (Pointer) omdp;
          context->index = R + 1;
        }
        if (index == n) {
          if (last_found != NULL) {
            *last_found = index;
          }
          return sfp;
        } else {
          if (last_found != NULL) {
            *last_found = index;
          }
          index++;
        }
      }
    }

    R++;
    feat = array [R];
  }

  return NULL;
}

typedef struct findnthfeatdata {
  CharPtr                findme;
  Uint1                  seqFeatChoice;
  Uint1                  featDefChoice;
  Int4                   n;
  Int4                   passed_so_far;
  SeqMgrFeatContext PNTR context;
  SeqFeatPtr             sfp;
} FindNthFeatData, PNTR FindNthFeatPtr;

static void FindNthFeatureByLabelInSeqEntryBspProc (BioseqPtr bsp, Pointer userdata)
{
  FindNthFeatPtr fnfp;
  Int4           this_search_index;
  Int4           passed_this_bsp = 0;

  if (bsp == NULL || (fnfp = (FindNthFeatPtr)userdata) == NULL
      || fnfp->sfp != NULL) {
    return;
  }

  this_search_index = fnfp->n - fnfp->passed_so_far;

  if (fnfp->seqFeatChoice == SEQFEAT_GENE || fnfp->featDefChoice == FEATDEF_GENE) {  
    fnfp->sfp = FindNthGeneOnBspByLabelOrLocusTag (bsp,
                                                   fnfp->findme,
                                                   this_search_index,
                                                   &passed_this_bsp,
                                                   fnfp->context);
  } else {
    fnfp->sfp = FindNthFeatureOnBspByLabel (bsp,
                                            fnfp->findme,
                                            fnfp->seqFeatChoice,
                                            fnfp->featDefChoice,
                                            this_search_index,
                                            &passed_this_bsp,
                                            fnfp->context);
  }
  if (fnfp->sfp == NULL) {
    fnfp->passed_so_far += passed_this_bsp;
  }
}

static SeqFeatPtr FindNthFeatureByLabelInSeqEntry
(SeqEntryPtr            sep,
 CharPtr                findme,
 Uint1                  seqFeatChoice,
 Uint1                  featDefChoice,
 Int4                   n,
 SeqMgrFeatContext PNTR context)
{
  FindNthFeatData fnf;
  
  fnf.findme = findme;
  fnf.seqFeatChoice = seqFeatChoice;
  fnf.featDefChoice = featDefChoice;
  fnf.n = n;
  fnf.context = context;
  fnf.passed_so_far = 0;
  fnf.sfp = NULL;

  VisitBioseqsInSep (sep, &fnf, FindNthFeatureByLabelInSeqEntryBspProc);

  return fnf.sfp;
}

static void FindByLabelOrPosProc (ButtoN b)

{
  Boolean            already;
  BioseqPtr          bsp;
  SeqMgrFeatContext  context;
  FindFormPtr        ffp;
  CharPtr            findme;
  SelStructPtr       sel;
  SeqEntryPtr        sep = NULL;
  SeqFeatPtr         sfp = NULL;
  Int4               val;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;

  if (ffp->input_itemID == 0) {
    sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
  } else {
    bsp = GetBioseqGivenIDs (ffp->input_entityID, ffp->input_itemID, ffp->input_itemtype);
    if (bsp == NULL) return;
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (bsp->repr == Seq_repr_seg) {
      sep = GetBestTopParentForData (ffp->input_entityID, bsp);
    }
  }
  if (sep == NULL) return;

  findme = SaveStringFromText (ffp->findTxt);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  CompressSpaces (findme);

  switch (ffp->type) {
    case FIND_GENE :
      sfp = FindNthFeatureByLabelInSeqEntry (sep, findme, SEQFEAT_GENE, 0,
                                             ffp->last_paragraph_found + 1,
                                             &context);
      if (sfp == NULL) {
        if (ffp->last_paragraph_found > -1) {
          Message (MSG_OK, "No more found");
        } else {
          Message (MSG_OK, "No matches found");
        }
        ffp->last_paragraph_found = -1;
      } else {
        ffp->last_paragraph_found ++;
      }
      break;
    case FIND_PROT :
      sfp = FindNthFeatureByLabelInSeqEntry (sep, findme, SEQFEAT_CDREGION, 0,
                                             ffp->last_paragraph_found + 1,
                                             &context);
      if (sfp == NULL) {
        ffp->last_paragraph_found = -1;
      } else {
        ffp->last_paragraph_found ++;
      }
      break;
    case FIND_POS :
      if (StrToLong (findme, &val)) {
        if (val > 0 && val <= bsp->length) {
          val--;
          sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
          while (sfp != NULL) {
            if (context.left >= val) break;
            sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
          }
        }
      }
      break;
    default :
      break;
  }

  MemFree (findme);

  if (sfp == NULL) return;

  already = FALSE;
  for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
    if (sel->entityID == context.entityID &&
        sel->itemID == context.itemID &&
        sel->itemtype == OBJ_SEQFEAT) {
      already = TRUE;
    }
  }
  if (! already) {
    ObjMgrAlsoSelect (context.entityID, context.itemID, OBJ_SEQFEAT, 0, NULL);
  }

}

static void FindByLabelOrPosProcFindFirst (ButtoN b)
{
  FindFormPtr ffp;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  ffp->last_paragraph_found = -1;
  FindByLabelOrPosProc (b);
}

static void ClearFindTextProc (ButtoN b)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  SafeSetTitle (ffp->findTxt, "");
  SafeSetTitle (ffp->replaceTxt, "");
  if (ffp->type == FIND_FLAT) {
    ObjMgrDeSelect (0, 0, 0, 0, NULL);
    ffp->last_paragraph_found = -1;
  }
  FindTextProc (ffp->findTxt);
  Select (ffp->findTxt);
}

static void FindFormMessage (ForM f, Int2 mssg)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void CopyFindReplBtn (ButtoN b)
{
  FindFormPtr ffp;
  CharPtr     str;
  
  ffp = (FindFormPtr) GetObjectExtra (b); 
  if (ffp == NULL)
  {
    return;
  }
  str = JustSaveStringFromText (ffp->findTxt);
  SetTitle (ffp->replaceTxt, str);
  str = MemFree (str);
}

static ForM CreateFindForm (Int2 left, Int2 top, CharPtr title,
                            Uint2 entityID, Uint4 itemID, Uint2 itemtype,
                            Int2 type)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  FindFormPtr        ffp;
  GrouP              j;
  GrouP              q = NULL;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  ffp = MemNew (sizeof (FindForm));
  if (ffp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->formmessage = FindFormMessage;
    ffp->input_entityID = entityID;
    ffp->input_itemID = itemID;
    ffp->input_itemtype = itemtype;
    ffp->type = type;
    ffp->last_paragraph_found = -1;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      ffp->appmessage = sepp->handleMessages;
    }

    j = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (j, 10, 10);

    g = HiddenGroup (j, 3, 0, NULL);
    StaticPrompt (g, "Find", 0, dialogTextHeight, programFont, 'l');
    ffp->findTxt = DialogText (g, "", 25, FindTextProc);
    SetObjectExtra (ffp->findTxt, ffp, NULL);
    StaticPrompt (g, "", 0, 0, programFont, 'l');
    if (type == FIND_ASN) {
      StaticPrompt (g, "Replace", 0, dialogTextHeight, programFont, 'l');
      ffp->replaceTxt = DialogText (g, "", 25, NULL);
      SetObjectExtra (ffp->replaceTxt, ffp, NULL);
      b = PushButton (g, "Copy", CopyFindReplBtn);
      SetObjectExtra (b, ffp, NULL);
    }

    if (type == FIND_ASN || type == FIND_FLAT) {
      q = HiddenGroup (w, 3, 0, NULL);
      ffp->caseCounts = CheckBox (q, "Case Sensitive", NULL);
      ffp->wholeWord = CheckBox (q, "Entire Word", NULL);
      if (indexerVersion && type == FIND_ASN) {
        ffp->doSeqIdLocal = CheckBox (q, "SeqID LOCAL", NULL);
      }
    }

    c = HiddenGroup (w, 5, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    if (type == FIND_ASN) {
      ffp->findAllBtn = DefaultButton (c, "Find All", FindAllProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
      ffp->replaceAllBtn = PushButton (c, "Replace All", ReplaceAllProc);
      SetObjectExtra (ffp->replaceAllBtn, ffp, NULL);
      Disable (ffp->replaceAllBtn);
    } else if (type == FIND_FLAT) {
      ffp->findFirstBtn = PushButton (c, "Find First", FindFlatProcFirst);
      SetObjectExtra (ffp->findFirstBtn, ffp, NULL);
      Disable (ffp->findFirstBtn);
      ffp->findNextBtn = PushButton (c, "Find Next", FindFlatProcNext);
      SetObjectExtra (ffp->findNextBtn, ffp, NULL);
      Disable (ffp->findNextBtn);
      ffp->findAllBtn = DefaultButton (c, "Find All", FindFlatProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
    } else if (type == FIND_GENE) {
      ffp->findFirstBtn = DefaultButton (c, "Find First Gene", FindByLabelOrPosProcFindFirst);
      SetObjectExtra (ffp->findFirstBtn, ffp, NULL);
      Disable (ffp->findFirstBtn);
      ffp->findAllBtn = PushButton (c, "Find Next Gene", FindByLabelOrPosProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
    } else if (type == FIND_PROT) {
     ffp->findFirstBtn = DefaultButton (c, "Find First Protein", FindByLabelOrPosProcFindFirst);
      SetObjectExtra (ffp->findFirstBtn, ffp, NULL);
      Disable (ffp->findFirstBtn);
      ffp->findAllBtn = PushButton (c, "Find Next Protein", FindByLabelOrPosProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
    } else if (type == FIND_POS) {
      ffp->findAllBtn = DefaultButton (c, "Find by Position", FindByLabelOrPosProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
    }
    b = PushButton (c, "Clear", ClearFindTextProc);
    SetObjectExtra (b, ffp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) j, (HANDLE) c, (HANDLE) q, NULL);

    RealizeWindow (w);
    Select (ffp->findTxt);
  }
  return (ForM) w;
}

extern void FindStringProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Find", bfp->input_entityID,
                        bfp->input_itemID, bfp->input_itemtype, FIND_ASN);
    Show (w);
    Select (w);
  }
}

extern void FindStringProcToolBtn (ButtoN b)

{
  BaseFormPtr  bfp;
  ForM         w;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  w = CreateFindForm (-90, -66, "Find", bfp->input_entityID,
                      bfp->input_itemID, bfp->input_itemtype, FIND_ASN);
  Show (w);
  Select (w);

}


extern void FindFlatfileProcToolBtn (ButtoN b)

{
  BaseFormPtr  bfp;
  ForM         w;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  w = CreateFindForm (-90, -66, "Flat File Find", bfp->input_entityID,
                      bfp->input_itemID, bfp->input_itemtype, FIND_FLAT);
  Show (w);
  Select (w);

}


extern void FindFlatfileProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Flat File Find", bfp->input_entityID,
                        bfp->input_itemID, bfp->input_itemtype, FIND_FLAT);
    Show (w);
    Select (w);
  }
}

extern void FindGeneProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Find", bfp->input_entityID,
                        bfp->input_itemID, bfp->input_itemtype, FIND_GENE);
    Show (w);
    Select (w);
  }
}

extern void FindProtProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Find", bfp->input_entityID,
                        bfp->input_itemID, bfp->input_itemtype, FIND_PROT);
    Show (w);
    Select (w);
  }
}

extern void FindPosProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Find", bfp->input_entityID,
                        bfp->input_itemID, bfp->input_itemtype, FIND_POS);
    Show (w);
    Select (w);
  }
}

static void RemoveRedundantSourceNote (CharPtr str, CharPtr keyword, CharPtr name)

{
  Char     ch;
  Boolean  extract;
  CharPtr  tmp1, tmp2;

  if (str == NULL || keyword == NULL || name == NULL) return;
  while (str != NULL && *str != '\0') {
    extract = TRUE;
    tmp1 = str;
    ch = *tmp1;
    while (ch == ' ' || ch == ';') {
      tmp1++;
      ch = *tmp1;
    }
    tmp2 = keyword;
    ch = TO_UPPER (*tmp1);
    while (ch != '\0' && ch == TO_UPPER (*tmp2)) {
      tmp1++;
      tmp2++;
      ch = TO_UPPER (*tmp1);
    }
    if (*tmp2 == '\0' && ch == ' ') {
      while (ch != '\0' && ch == ' ') {
        tmp1++;
        ch = *tmp1;
      }
      tmp2 = name;
      while (ch != '\0' && ch == *tmp2) {
        tmp1++;
        tmp2++;
        ch = *tmp1;
      }
      if (*tmp2 == '\0' && (ch == '\0' || ch == ' ' || ch == ';')) {
        while (ch != '\0' && ch == ' ') {
          tmp1++;
          ch = *tmp1;
        }
        /* now ready to extract */
      } else {
        extract = FALSE;
      }
    } else {
      extract = FALSE;
    }
    if (extract) {
      while (ch == ' ' || ch == ';') {
        tmp1++;
        ch = *tmp1;
      }
      tmp2 = str;
      while (ch != '\0') {
        *tmp2 = ch;
        tmp1++;
        tmp2++;
        ch = *tmp1;
      }
      *tmp2 = '\0';
    } else {
      ch = *tmp1;
      while (ch != '\0' && ch != ';') {
        tmp1++;
        ch = *tmp1;
      }
      if (ch == ';') {
        tmp1++;
        ch = *tmp1;
      }
      while (ch != '\0' && (ch == ' ' || ch == ';')) {
        tmp1++;
        ch = *tmp1;
      }
      str = tmp1;
    }
  }
}

static CharPtr TrimSpacesAndSemicolonsAroundString (CharPtr str)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch <= ' ' || ch == ';')) {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ' && ch != ';') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

extern void GetRidOfRedundantSourceNotes (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void GetRidOfRedundantSourceNotes (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  ValNodePtr         sdp;
  SubSourcePtr       ssp;
  CharPtr            str1, str2;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_source) {
      str1 = NULL;
      str2 = NULL;
      onp = NULL;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          onp = orp->orgname;
          if (onp != NULL) {
            mod = onp->mod;
            while (mod != NULL) {
              if (mod->subtype == 255) {
                str1 = mod->subname;
              }
              mod = mod->next;
            }
          }
        }
        ssp = biop->subtype;
        while (ssp != NULL) {
          if (ssp->subtype == 255) {
            str2 = ssp->name;
          }
          ssp = ssp->next;
        }
        if (str1 != NULL || str2 != NULL) {
          if (onp != NULL) {
            mod = onp->mod;
            while (mod != NULL) {
              if (mod->subtype != 255) {
                RemoveRedundantSourceNote (str1, GetOrgModQualName (mod->subtype), mod->subname);
                RemoveRedundantSourceNote (str2, GetOrgModQualName (mod->subtype), mod->subname);
              }
              mod = mod->next;
            }
          }
          ssp = biop->subtype;
          while (ssp != NULL) {
            if (ssp->subtype != 255) {
              RemoveRedundantSourceNote (str1, GetSubsourceQualName (ssp->subtype), ssp->name);
              RemoveRedundantSourceNote (str2, GetSubsourceQualName (ssp->subtype), ssp->name);
            }
            ssp = ssp->next;
          }
        }
      }
      TrimSpacesAndSemicolonsAroundString (str1);
      TrimSpacesAndSemicolonsAroundString (str2);
    }
    sdp = sdp->next;
  }
}

typedef struct convaccdata {
  CharPtr        currID;
  CharPtr        newID;
} ConvAccData, PNTR ConvAccPtr;

static void ChangeSeqIdListToAccLocalID (SeqIdPtr sip, CharPtr currID, CharPtr newID)

{
  ObjectIdPtr   oip;

  while (sip != NULL) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        oip = (ObjectIdPtr) sip->data.ptrvalue;
        if (oip != NULL) {
          if (StringCmp (oip->str, currID) == 0) {
            MemFree (oip->str);
            oip->str = StringSave (newID);
          }
        }
        break;
      default :
        break;
    }
    sip = sip->next;
  }
}

static void ChangeSeqLocListToAccLocalID (SeqLocPtr slp, CharPtr currID, CharPtr newID)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        ChangeSeqIdListToAccLocalID (sip, currID, newID);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          ChangeSeqIdListToAccLocalID (loc, currID, newID);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdListToAccLocalID (sip, currID, newID);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdListToAccLocalID (sip, currID, newID);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
}

static void ChangeAlignListToAccLocalID (SeqAlignPtr align, CharPtr currID, CharPtr newID)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  StdSegPtr     ssp;

  if (align == NULL) return;
  if (align->segtype == 1) {
    ddp = (DenseDiagPtr) align->segs;
    if (ddp != NULL) {
      ChangeSeqIdListToAccLocalID (ddp->id, currID, newID);
    }
  } else if (align->segtype == 2) {
    dsp = (DenseSegPtr) align->segs;
    if (dsp != NULL) {
      ChangeSeqIdListToAccLocalID (dsp->ids, currID, newID);
    }
  } else if (align->segtype == 3) {
    ssp = (StdSegPtr) align->segs;
    if (ssp != NULL) {
       ChangeSeqLocListToAccLocalID (ssp->loc, currID, newID);
    }
  }
}

static void CopyAccToLocalCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ConvAccPtr    cap;
  SeqAlignPtr   sal;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  cap = (ConvAccPtr) mydata;
  if (cap == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sip = bsp->id;
    ChangeSeqIdListToAccLocalID (sip, cap->currID, cap->newID);
    SeqMgrReplaceInBioseqIndex (bsp);
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        slp = sfp->location;
        ChangeSeqLocListToAccLocalID (slp, cap->currID, cap->newID);
        slp = sfp->product;
        ChangeSeqLocListToAccLocalID (slp, cap->currID, cap->newID);
        sfp = sfp->next;
      }
    } else if (sap->type == 2) {
      sal = (SeqAlignPtr) sap->data;
      while (sal != NULL) {
        ChangeAlignListToAccLocalID (sal, cap->currID, cap->newID);
        sal = sal->next;
      }
    }
    sap = sap->next;
  }
}

static void RecordAccToConvert (ValNodePtr PNTR vnpp, BioseqPtr bsp, CharPtr newID)

{
  ConvAccPtr   cap;
  ObjectIdPtr  oip;
  SeqIdPtr     sip;

  if (vnpp == NULL || bsp == NULL || StringHasNoText (newID)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL) {
        if (! StringHasNoText (oip->str)) {
          cap = (ConvAccPtr) MemNew (sizeof (ConvAccData));
          if (cap == NULL) return;
          cap->currID = StringSave (oip->str);
          cap->newID = StringSave (newID);
          ValNodeAddPointer (vnpp, 0, (Pointer) cap);
          return;
        }
      }
    }
  }
}

static Boolean IsLegalAccession (CharPtr acnum, Int2 i)

{
  Int2  j;

  if (! isalpha ((Int4)(acnum [0]))) return FALSE;
  if (!(isdigit((Int4)(acnum[1])) && i == 6) && !(isalpha((Int4)(acnum[1])) && i == 8)) return FALSE;
  for (j = 2; j < i; j++) {
    if (!(isdigit((Int4)(acnum[j])))) return FALSE;
  }
  return TRUE;
}

static void FindAccInDefCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  Char         acnum [17];
  BioseqPtr    bsp;
  Int2         i;
  CharPtr      p;
  Char         str [128];
  ValNodePtr   ttl;
  ValNodePtr   PNTR vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (! IS_Bioseq (sep)) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
  ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  while (ttl != NULL) {
    StringNCpy_0 (str, (CharPtr) ttl->data.ptrvalue, sizeof (str));
    TrimSpacesAroundString (str);
    if (! StringHasNoText (str)) {
      if (str [0] == 'a' && str [1] == 'c' && str [2] == 'c') {
        p = str + 3;
        for (i = 0; isalnum ((Int4)(*p)) && *p != '\0'; p++, i++) {
          acnum [i] = *p;
        }
        acnum [i] = '\0';
        if (i == 6 || i == 8) {
          if (IsLegalAccession (acnum, i)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            RecordAccToConvert (vnpp, bsp, str);
          }
        }
      }
    }
    ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, ttl);
  }
}

static void ConvertAccInDefToLocalIdProc (SeqEntryPtr sep)

{
  ConvAccPtr  cap;
  ValNodePtr  head;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  head = NULL;
  SeqEntryExplore (sep, (Pointer) &head, FindAccInDefCallback);
  if (head == NULL) return;
  if (Message (MSG_YN, "Convert accLNNNNN or accLLNNNNNN in title to local ID?") == ANS_NO) return;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    cap = (ConvAccPtr) vnp->data.ptrvalue;
    if (cap != NULL) {
      SeqEntryExplore (sep, (Pointer) cap, CopyAccToLocalCallback);
    }
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    cap = (ConvAccPtr) vnp->data.ptrvalue;
    if (cap != NULL) {
      MemFree (cap->currID);
      MemFree (cap->newID);
      vnp->data.ptrvalue = MemFree (cap);
    }
  }
  ValNodeFree (head);
}

static Int2  taxonCount;

static Int4 DoSeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct,
                              Boolean force, Boolean dotaxon, MonitorPtr mon)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Char          str [32];

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += DoSeqEntryToAsn3 (sep, strip, correct, force, dotaxon, mon);
      }
      return rsult;
    }
  }
/*#ifdef USE_TAXON*/
  if (dotaxon && mon != NULL) {
    taxonCount++;
    sprintf (str, "Processing Component %d", (int) taxonCount);
    MonitorStrValue (mon, str);
  }
/*#endif*/
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    if ((! force) && (! NoBiosourceOrTaxonId (sep))) return 0;
  } else {
/*#else*/
    if ((! force) && SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL) != NULL) return 0;
  }
/*#endif*/
  oldscope = SeqEntrySetScope (sep);
  if (dotaxon) {
    rsult = SeqEntryToAsn3Ex (sep, strip, correct, TRUE, NULL, Tax3MergeSourceDescr);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);
  } else {
    rsult = SeqEntryToAsn3Ex (sep, strip, correct, FALSE, NULL, NULL);
  }
  SeqEntrySetScope (oldscope);
  return rsult;
}

/* CheckSeqAlignCallback copied from salsa.c */
static void CheckSeqAlignCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAnnotPtr        sap,
                     pre;

  if (sep != NULL && sep->data.ptrvalue) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           pre=NULL;
           sap=bsp->annot;
           while (sap!= NULL)
           {
              if (sap->type == 2) {
                 if (is_dim1seqalign ((SeqAlignPtr) sap->data)) {
                    if (pre==NULL) {
                       bsp->annot=sap->next;
                       sap->next=NULL;
                       SeqAnnotFree (sap);
                       sap=bsp->annot; 
                    }
                    else {
                       pre->next=sap->next;
                       pre=sap->next;
                       sap->next=NULL; 
                       SeqAnnotFree (sap);
                       sap=pre;
                    }
                 }
                 else {
                    pre=sap;
                    sap=sap->next;
                 }
              }
              else {
                 pre=sap;
                 sap=sap->next;
              }
           }
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           pre=NULL;
           sap=bssp->annot;
           while (sap!= NULL)
           {
              if (sap->type == 2) {
                 if (is_dim1seqalign ((SeqAlignPtr) sap->data)) {
                    if (pre==NULL) {
                       bssp->annot=sap->next;
                       sap->next=NULL;
                       SeqAnnotFree (sap);
                       sap=bssp->annot; 
                    }
                    else {
                       pre=sap->next;
                       sap->next=NULL; 
                       SeqAnnotFree (sap);
                       sap=sap->next;
                    }
                 }
                 else {
                    pre=sap;
                    sap=sap->next;
                 }
              }
              else {
                 pre=sap;
                 sap=sap->next;
              }
           }
        }
     }
  }
}

/* RemoveMultipleTitles currently removes FIRST title in chain */

static void RemoveMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqDescrPtr    descr = NULL;
  SeqDescrPtr    lasttitle = NULL;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    descr = bssp->descr;
  } else return;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_title) continue;
    if (lasttitle != NULL) {
      if (lasttitle->extended != 0) {
        ovp = (ObjValNodePtr) lasttitle;
        ovp->idx.deleteme = TRUE;
      }
      lasttitle = sdp;
    } else {
      lasttitle = sdp;
    }
  }
}

extern Int4 MySeqEntryToAsn3Ex (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force, Boolean dotaxon);
extern Int4 MySeqEntryToAsn3Ex (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force, Boolean dotaxon)

{
  Uint2       entityID;
  MonitorPtr  mon;
  Boolean     needstaxfix;
  Int4        rsult;
  ErrSev      sev;

  rsult = 0;
  sev = ErrSetMessageLevel (SEV_FATAL);
  BasicSeqEntryCleanup (sep);
  EntryChangeImpFeat(sep);     /* change any CDS ImpFeat to real CdRegion */
  /* NormalizePeriodsOnInitials (sep); */ /* put periods on author initials */
  /* MoveRnaGBQualProductToName (sep); */ /* move rna gbqual product to rna-ref.ext.name */
  /* MoveProtGBQualProductToName (sep); */ /* move prot gbqual product to prot-ref.name */
  /* MoveCdsGBQualProductToName (sep); */ /* move cds gbqual product to prot-ref.name */
  /* MoveFeatGBQualsToFields (sep); */ /* move feature partial, exception to fields */
  if (indexerVersion) {
    /*
    StripTitleFromProtsInNucProts (sep);
    */
    MoveFeatsFromPartsSet (sep);
    move_cds (sep); /* move CDS features to nuc-prot set */
  }
  /* ExtendGeneFeatIfOnMRNA (0, sep); */ /* gene on mRNA is full length */
  entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  VisitBioseqsInSep (sep, NULL, ExtendSingleGeneOnMRNA);
  RemoveBioSourceOnPopSet (sep, NULL);
  /* SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback); */
  SeqEntryExplore (sep, NULL, DeleteMultipleTitles); /* do it old way in Sequin */
  /*
  SeqEntryExplore (sep, NULL, RemoveMultipleTitles);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
  */
  SeqEntryPack (sep);
  if (indexerVersion) {
    ConvertAccInDefToLocalIdProc (sep); /* if title is accXNNNNN, convert to seqID */
  }
  if (useEntrez) {
    ValidateSeqAlignandACCInSeqEntry (sep, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE); /* remove components with seqID accXNNNNN */
  }
  SeqEntryExplore (sep, NULL, CheckSeqAlignCallback); /* remove alignments with single dimension */
  needstaxfix = FALSE;
  if (! force) {
    needstaxfix = NoBiosourceOrTaxonId (sep);
  }
  if ((! force) && (! needstaxfix)) {
    ConvertFullLenSourceFeatToDesc (sep);
    ConvertFullLenPubFeatToDesc (sep);
    /* EntryStripSerialNumber(sep); */ /* strip citation serial numbers */
    EntryChangeGBSource (sep);   /* at least remove redundant information in GBBlocks */
    EntryCheckGBBlock (sep);
    /* SeqEntryMoveDbxrefs (sep); */ /* db_xref gbqual to sfp->dbxref */
    EntryMergeDupBioSources (sep);
    GetRidOfEmptyFeatsDescStrings (0, sep);
    GetRidOfLocusInSeqIds (0, sep);
    /* reindex, since CdEndCheck (from CdCheck) gets best overlapping gene */
    entityID = ObjMgrGetEntityIDForChoice (sep);
    SeqMgrIndexFeatures (entityID, NULL);
    CdCheck (sep, NULL);
    BasicSeqEntryCleanup (sep);
    ErrSetMessageLevel (sev);
    return rsult;
  }
  if (force && useTaxon) {
    dotaxon = TRUE;
  }
  mon = NULL;
  taxonCount = 0;
/*#ifdef USE_TAXON*/
  if (dotaxon) {
    WatchCursor ();
    mon = MonitorStrNewEx ("Taxonomy Lookup", 40, FALSE);
    MonitorStrValue (mon, "Processing Organism Info");
    Update ();
  }
/*#endif*/

  /* set dirty flag if no lineage or division in any biosource */
  if (dotaxon && needstaxfix) {
    entityID = ObjMgrGetEntityIDForChoice (sep);
    ObjMgrSetDirtyFlag (entityID, TRUE);
  }

  EntryMergeDupBioSources (sep); /* do before and after SE2A3 */
  if (dotaxon) {
    Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  }
  rsult = DoSeqEntryToAsn3 (sep, strip, correct, force, FALSE, mon);
/*#ifdef USE_TAXON*/
  if (dotaxon) {
    MonitorStrValue (mon, "Closing Taxon");
    Update ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
  }
/*#endif*/
  ConvertFullLenSourceFeatToDesc (sep);
  ConvertFullLenPubFeatToDesc (sep);
  /* EntryStripSerialNumber(sep); */ /* strip citation serial numbers */
  MovePopPhyMutPubs (sep);
  EntryChangeGBSource (sep);   /* remove redundant information in GBBlocks again */
  EntryCheckGBBlock (sep);
  /* SeqEntryMoveDbxrefs (sep); */ /* db_xref gbqual to sfp->dbxref */
  EntryMergeDupBioSources (sep);
  GetRidOfEmptyFeatsDescStrings (0, sep);
  GetRidOfLocusInSeqIds (0, sep);
  /* reindex, since CdEndCheck (from CdCheck) gets best overlapping gene */
  entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrClearFeatureIndexes (entityID, NULL);
  NormalizeDescriptorOrder (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  CdCheck (sep, NULL);
  BasicSeqEntryCleanup (sep);
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
  return rsult;
}

typedef struct partialformdata {
  FEATURE_FORM_BLOCK

  ValNodePtr         featlist;
  LisT               feature;
  PopuP              change5;
  PopuP              change3;
  PopuP              orderJoinState;
  GrouP              nucprotchoice;
  TexT               findThis;
  CharPtr            findThisStr;
  Int2               subtype;
  Int2               leftpolicy;
  Int2               rightpolicy;
  Int2               nucprotpolicy;
  Int2               orderjoinpolicy;
  ObjMgrPtr          omp;
  ObjMgrTypePtr      omtp;
  ButtoN             extend5_btn;
  ButtoN             extend3_btn;
  Boolean            extend5;
  Boolean            extend3;
  ButtoN             leaveDlgUp;
  Boolean            case_insensitive;
  ButtoN             case_insensitive_btn;
  Boolean            when_string_not_present;
  ButtoN             when_string_not_present_btn;
} PartialFormData, PNTR PartialFormPtr;

static Boolean CDSMeetsStringConstraint (SeqFeatPtr sfp,
				      CharPtr     findThisStr,
				      Boolean     case_insensitive)
{
  BioseqPtr		protbsp;
  SeqFeatPtr		protsfp;
  SeqMgrFeatContext	context;
  ProtRefPtr		prp;

  if (sfp == NULL) return FALSE;
  protbsp = BioseqFindFromSeqLoc (sfp->product);
  if (protbsp == NULL) return FALSE;
  protsfp = SeqMgrGetBestProteinFeature (protbsp, &context);
  if ((case_insensitive && StringISearch (context.label, findThisStr) != NULL)
    || (!case_insensitive && StringSearch (context.label, findThisStr) != NULL))
  {
    return TRUE;
  }
  if (protsfp == NULL) return FALSE;
  prp = (ProtRefPtr) protsfp->data.value.ptrvalue;
  if (prp->name != NULL)
  {
    if ((case_insensitive && StringISearch (prp->name->data.ptrvalue, findThisStr) != NULL)
      || (!case_insensitive && StringSearch (prp->name->data.ptrvalue, findThisStr) != NULL))
    {
      return TRUE;
    }
  }
  if (prp->desc != NULL)
  {
  	if ((case_insensitive && StringISearch (prp->desc, findThisStr) != NULL)
  	  || (!case_insensitive && StringSearch (prp->desc, findThisStr) != NULL))
  	{
  	  return TRUE;
  	}
  }
  return FALSE;
}

extern Boolean MeetsStringConstraint (SeqFeatPtr  sfp,
				                      CharPtr     findThisStr,
				                      Boolean     case_insensitive)
{
  GBQualPtr         gbqp;
  GeneRefPtr        grp;
  RnaRefPtr         rrp;
  SeqMgrFeatContext context;
  Boolean           have_context = FALSE;

  /* If no string constraint, then everyone matches */

  if (NULL == findThisStr)
    return TRUE;

  /* Search for the string constraint */
  /* in the feature title field */
  if (StringISearch (sfp->title, findThisStr))
    return TRUE;

  /* Search for the string constraint */
  /* in the feature comment field.    */

  if (StringISearch (sfp->comment, findThisStr))
    return TRUE;

  /* Search for the string constraint */
  /* in GB qualifiers.                */

  gbqp = sfp->qual;
  while (NULL != gbqp)
    {
      if ((NULL != gbqp->val) && StringISearch (gbqp->val, findThisStr))
	return TRUE;
      gbqp = gbqp->next;
    }

  if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) != NULL)
  {
    if (!case_insensitive && StringSearch (context.label, findThisStr) != NULL)
    {
    	return TRUE;
    }
    else if (case_insensitive && StringISearch (context.label, findThisStr) != NULL)
  	{
  	  return TRUE;
  	}
  	have_context = TRUE;
  }

  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = sfp->data.value.ptrvalue;
    if (!case_insensitive && 
        (StringSearch (grp->locus, findThisStr) != NULL
        || StringSearch (grp->desc, findThisStr) != NULL
        || StringSearch (grp->desc, findThisStr) != NULL
        || StringSearch (grp->locus_tag, findThisStr) != NULL))
    {
      return TRUE;    	
    }
    else if (case_insensitive &&
        (StringISearch (grp->locus, findThisStr) != NULL
        || StringISearch (grp->desc, findThisStr) != NULL
        || StringISearch (grp->desc, findThisStr) != NULL
        || StringISearch (grp->locus_tag, findThisStr) != NULL))
    {
      return TRUE;    	
    }
  }
  else if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    if (CDSMeetsStringConstraint (sfp, findThisStr, case_insensitive))
      return TRUE;
  }
  else if (sfp->data.choice == SEQFEAT_RNA)
  {
    rrp = sfp->data.value.ptrvalue;

    if (rrp->ext.choice == 1) {
      if ((!case_insensitive && StringSearch ((CharPtr) rrp->ext.value.ptrvalue, findThisStr) != NULL)
        || (case_insensitive && StringISearch ((CharPtr) rrp->ext.value.ptrvalue, findThisStr) != NULL))
      {
        return TRUE;
      }
    }
    else if (rrp->type == 3 && rrp->ext.choice == 2 && have_context) 
    {
      /* look for the label as it appears to the user */
      if ((!case_insensitive && StringNCmp(findThisStr, "tRNA-", 5) == 0
          && StringSearch (context.label, findThisStr + 5))
          || (case_insensitive && StringNICmp (findThisStr, "tRNA-", 5) == 0
          && StringISearch (context.label, findThisStr + 5)))
      {
      	return TRUE;
      }
    }
  }

  /* If we got to here, then the string constraint was not found */

  return FALSE;
}


extern Boolean SetBestFrameByLocation (SeqFeatPtr sfp)
{
  CdRegionPtr  crp;
  Uint1        new_frame = 0, i;
  ByteStorePtr bs;
  Int4         lens [3];
  Int4         max;
  Boolean      retval = TRUE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;
  
  crp = sfp->data.value.ptrvalue;
  if (crp == NULL) return FALSE;

  max = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    lens[i - 1] = BSLen (bs);
    BSFree (bs);
    if (lens[i - 1] > max) {
      max = lens[i - 1];
      new_frame = i;
    }
  }
  for (i = 1; i <= 3; i++) {
    if (lens [i - 1] == max && i != new_frame) {
      retval = FALSE;
    }
  }
  crp->frame = new_frame;
  return retval;
}


extern void SetBestFrame (SeqFeatPtr sfp)
{
  SeqAlignPtr  salp;
  CdRegionPtr  crp;
  BioseqPtr    old_prot, new_prot;
  ByteStorePtr bs = NULL;
  Int4         original_frame, test_frame;
  Boolean      revcomp;
  ErrSev       level;
  CharPtr      seq_str1, seq_str2;
  Int4         best_len = -1, new_aln_len;
  Int4         best_frame = -1;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_CDS) return;
  
  crp = sfp->data.value.ptrvalue;
  if (crp == NULL) return;

  old_prot = BioseqFindFromSeqLoc (sfp->product);
  if (old_prot == NULL) return;
  
  new_prot = BioseqNew ();
  new_prot->id = SeqIdParse ("lcl|CdRgnTransl");
  new_prot->repr = Seq_repr_raw;
  new_prot->mol = Seq_mol_aa;
  new_prot->seq_data_type = Seq_code_ncbieaa;
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  new_prot->seq_data = (SeqDataPtr) bs;
  new_prot->length = BSLen (bs);
  
  original_frame = crp->frame;

  /* suppress BLAST error messages when no similarity is found */
  level = ErrSetMessageLevel (SEV_MAX);
  seq_str1 = BSMerge((ByteStorePtr)(old_prot->seq_data), NULL);
  seq_str2 = BSMerge((ByteStorePtr)(new_prot->seq_data), NULL);
  
  for (test_frame = 1; test_frame <= 3; test_frame ++)
  {
    new_prot->seq_data = SeqDataFree (new_prot->seq_data, new_prot->seq_data_type);
    crp->frame = test_frame;
    new_prot->seq_data = (SeqDataPtr) ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    salp = Sequin_GlobalAlign2Seq (old_prot, new_prot, &revcomp);
    if (salp != NULL)
    {
      new_aln_len = SeqAlignLength (salp);
      if (new_aln_len > best_len)
      {
        best_len = new_aln_len;
        best_frame = test_frame;
      }
      salp = SeqAlignFree (salp);
    }
  }  	
  
  if (best_frame > -1)
  {
    crp->frame = best_frame;
  }
  else
  {
    crp->frame = original_frame;
  }

  ErrSetMessageLevel (level);
  BioseqFree (new_prot);
}


/*
static Boolean AddOrgToDefGatherFunc (GatherContextPtr gcp)

{
  CharPtr     def;
  CharPtr     ptr;
  ValNodePtr  sdp;
  CharPtr     str;
  CharPtr     text;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQDESC) return TRUE;
  text = (CharPtr) gcp->userdata;
  if (text == NULL || StringHasNoText (text)) return TRUE;
  sdp = (ValNodePtr) gcp->thisitem;
  if (sdp->choice != Seq_descr_title) return TRUE;
  def = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (def)) return TRUE;

  ptr = StringISearch (def, text);
  if (ptr != NULL && ptr == def) return TRUE;
  str = MemNew ((StringLen (text) + StringLen (def) + 4) * sizeof (Char));
  if (str != NULL) {
    StringCpy (str, text);
    StringCat (str, " ");
    StringCat (str, def);
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = str;
    ObjMgrSetDirtyFlag (gcp->entityID, TRUE);
  }
  return TRUE;
}
*/

static void AppendOrgToString (Uint2 entityID, SeqDescrPtr sdp, CharPtr text)

{
  CharPtr     def;
  CharPtr     ptr;
  CharPtr     str;

  def = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (def)) return;

  ptr = StringISearch (def, text);
  if (ptr != NULL && ptr == def) return;
  str = MemNew ((StringLen (text) + StringLen (def) + 4) * sizeof (Char));
  if (str != NULL) {
    StringCpy (str, text);
    StringCat (str, " ");
    StringCat (str, def);
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = str;
    ObjMgrSetDirtyFlag (entityID, TRUE);
  }
}

static void AddOrgToDefElement (Uint2 entityID, SeqEntryPtr sep, Int2 orgmod, Int2 subsource)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Char               ch;
  SeqMgrDescContext  dcontext;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            ptr;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  Char               str [96];
  Char               text [64];
  CharPtr            title;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 1 || bssp->_class == 2 ||
                         bssp->_class == 4)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        AddOrgToDefElement (entityID, sep, orgmod, subsource);
      }
      return;
    }
  }
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  biop = NULL;
  text [0] = '\0';
  str [0] = '\0';
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  if (biop == NULL) return;
  /* SeqEntryToBioSource (sep, NULL, str, sizeof (str) - 1, &biop); */
  if (orgmod == 0 && subsource == 0) {
    orp = biop->org;
    if (orp == NULL) return;
    StringNCpy_0 (str, orp->taxname, sizeof (str));
    /*
    ptr = StringSearch (str, "(");
    if (ptr != NULL) {
      *ptr = '\0';
    }
    */
    TrimSpacesAroundString (str);
    if ((StringICmp (str, "Human immunodeficiency virus type 1") == 0) ||
	(StringICmp (str, "Human immunodeficiency virus 1") == 0)) {
      StringCpy (str, "HIV-1");
    } else if ((StringICmp (str,"Human immunodeficiency virus type 2")==0) ||
	       (StringICmp (str,"Human immunodeficiency virus 2") == 0)) {
      StringCpy (str, "HIV-2");
    }
    str [0] = TO_UPPER (str [0]);
  } else if (biop != NULL && biop->org != NULL) {
    text [0] = '\0';
    str [0] = '\0';
    orp = biop->org;
    if (orgmod > 0) {
      onp = orp->orgname;
      if (onp != NULL) {
        mod = onp->mod;
        while (mod != NULL) {
          if (mod->subtype == orgmod) {
            StringNCpy_0 (text, mod->subname, sizeof (text));
            StringNCpy_0 (str, GetOrgModQualName (mod->subtype), sizeof (str));
          }
          mod = mod->next;
        }
      }
    } else if (subsource > 0) {
      ssp = biop->subtype;
      while (ssp != NULL) {
        if (ssp->subtype == subsource) {
          StringNCpy_0 (text, ssp->name, sizeof (text));
          StringNCpy_0 (str, GetSubsourceQualName (ssp->subtype), sizeof (str));
        }
        ssp = ssp->next;
      }
    }
    if (StringHasNoText (text)) {
      str [0] = '\0';
      text [0] = '\0';
    } else {
      StringCat (str, " ");
      ptr = str;
      while (*ptr != '\0') {
        ch = *ptr;
        *ptr = TO_LOWER (ch);
        ptr++;
      }
      StringCat (str, text);
    }
  }
  /*
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.scope = sep;
  GatherSeqEntry (sep, (Pointer) str, AddOrgToDefGatherFunc, &gs);
  */
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (sdp == NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
    if (sdp == NULL) return;
    title = (CharPtr) sdp->data.ptrvalue;
    if (title == NULL) return;
    sdp = SeqDescrAdd (&(bsp->descr));
    if (sdp == NULL) return;
    sdp->choice = Seq_descr_title;
    sdp->data.ptrvalue = StringSave (title);
  }
  if (sdp == NULL) return;
  AppendOrgToString (entityID, sdp, str);
}

static void AddOrgToDef (Uint2 entityID, SeqEntryPtr sep, Int2 orgmod, Int2 subsource)

{
  BioseqSetPtr       bssp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        AddOrgToDef (entityID, sep, orgmod, subsource);
      }
      return;
    }
  }
  AddOrgToDefElement (entityID, sep, orgmod, subsource);
}

extern void CommonAddOrgOrModsToDefLines (IteM i, Int2 orgmod, Int2 subsource, ButtoN b)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  AddOrgToDef (bfp->input_entityID, sep, orgmod, subsource);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveAlignmentCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAlignPtr    nextsalp;
  SeqAnnotPtr    nextsap;
  Pointer PNTR   prevsalp;
  Pointer PNTR   prevsap;
  SeqAlignPtr    salp;
  SeqAnnotPtr    sap;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 2) {
      salp = (SeqAlignPtr) sap->data;
      prevsalp = (Pointer PNTR) &(sap->data);
      while (salp != NULL) {
        nextsalp = salp->next;
        *(prevsalp) = salp->next;
        salp->next = NULL;
        SeqAlignFree (salp);
        salp = nextsalp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveAlignment (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveAlignmentCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void RemoveGraphCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAnnotPtr    nextsap;
  SeqGraphPtr    nextsgp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsgp;
  SeqAnnotPtr    sap;
  SeqGraphPtr    sgp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        if (sgp->flags [2] >= 1 && sgp->flags [2] <= 3) {
          *(prevsgp) = sgp->next;
          sgp->next = NULL;
          SeqGraphFree (sgp);
        } else {
          prevsgp = (Pointer PNTR) &(sgp->next);
        }
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveGraph (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveGraphCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void RemoveSeqAnnotIDsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAnnotPtr    nextsap;
  Pointer PNTR   prevsap;
  SeqAnnotPtr    sap;
  SeqIdPtr       sip;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 4) {
      sip = (SeqIdPtr) sap->data;
      SeqIdSetFree (sip);
      sap->data = NULL;
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveSeqAnnotIDs (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveSeqAnnotIDsCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void RemoveSeqAnnotLOCsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAnnotPtr    nextsap;
  Pointer PNTR   prevsap;
  SeqAnnotPtr    sap;
  SeqLocPtr      slp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 4) {
      slp = (SeqLocPtr) sap->data;
      SeqLocSetFree (slp);
      sap->data = NULL;
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveSeqAnnotLOCs (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveSeqAnnotLOCsCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void MarkProteinCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  ValNodePtr PNTR  vnpp;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_aa (bsp->mol)) {
      ValNodeAddPointer (vnpp, 0, (Pointer) bsp);
    }
  }
}

static void RemoveProteinsAndOptionallyRenormalize (IteM i, Boolean renormalize)

{
  BaseFormPtr    bfp;
  BioseqPtr      bsp;
  Uint4          itemID;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  OMProcControl  ompc;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep;
  ValNodePtr     tmp;
  ValNodePtr     vnp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  vnp = NULL;
  SeqEntryExplore (sep, (Pointer) &vnp, MarkProteinCallback);
  for (tmp = vnp; tmp != NULL; tmp = tmp->next) {
    bsp = (BioseqPtr) tmp->data.ptrvalue;
    itemID = GetItemIDGivenPointer (bfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
    if (itemID > 0) {
      MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
      ompc.do_not_reload_from_cache = TRUE;
      ompc.input_entityID = bfp->input_entityID;
      ompc.input_itemID = itemID;
      ompc.input_itemtype = OBJ_BIOSEQ;
      if (! DetachDataForProc (&ompc, FALSE)) {
        Message (MSG_POSTERR, "DetachDataForProc failed");
      }
    }
  }
  ValNodeFree (vnp);
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  if (renormalize) 
  {
    RenormalizeNucProtSets (sep, TRUE);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

extern void RemoveProteins (IteM i)
{
  RemoveProteinsAndOptionallyRenormalize (i, FALSE);
}

extern void RemoveProteinsAndRenormalize (IteM i)
{
  RemoveProteinsAndOptionallyRenormalize (i, TRUE);
}

#define EDIT_FIVE_PRIME  1
#define EDIT_THREE_PRIME 2

#define ADD_TO_END       1
#define TRIM_FROM_END    2

#define TRIM_BY_SEQUENCE 1
#define TRIM_BY_COUNT    2

typedef struct edseqendsdata {
  FEATURE_FORM_BLOCK

  TexT           seq;
  TexT           genename;
  GrouP          whichend;
  GrouP          addOrTrim;
  Int2           addOrTrimBool;
  GrouP          trimBy;
  Int2           trimByBool;
  Int4           trimCount;
  TexT           trimCountText;
  ButtoN         extendfeat;
  LisT           nuc_sequence_list_ctrl;
  ButtoN         addCitSub;
  CharPtr        seqstr;
  CharPtr        genestr;
  Int2           endval;
  Int2           frameshift;
  Boolean        adjustframe;
  Boolean        extendflag;
  BioseqPtr      extendedthis;
  SeqEntryPtr    sep;
  Boolean        rsult;
} EditSeqEnds, PNTR EditSeqPtr;

static Boolean GeneFindByNameFunc (GatherContextPtr gcp)

{
  EditSeqPtr  esp;
  GeneRefPtr  grp;
  SeqFeatPtr  sfp;

  if (gcp == NULL) return TRUE;
  esp = (EditSeqPtr) gcp->userdata;
  if (esp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return TRUE;
  if (StringICmp (grp->locus, esp->genestr) == 0) {
    esp->rsult = TRUE;
  }
  return TRUE;
}

static Boolean EditSeqEntryHasGene (BioseqPtr bsp, SeqEntryPtr sep, EditSeqPtr esp)

{
  GatherScope  gs;

  if (esp->input_entityID == 0 || esp->sep == NULL) return FALSE;
  if (StringHasNoText (esp->genestr)) return TRUE;
  esp->rsult = FALSE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  gs.scope = sep;
  MemSet ((Pointer)(gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_BIOSEQ] = FALSE;
  gs.ignore [OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore [OBJ_SEQFEAT] = FALSE;
  gs.ignore [OBJ_SEQANNOT] = FALSE;
  GatherEntity (esp->input_entityID, (Pointer) esp, GeneFindByNameFunc, &gs);
  gs.target = SeqLocFree (gs.target);
  return esp->rsult;
}

static Boolean FixACDSFunc (GatherContextPtr gcp)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  CdRegionPtr   crp;
  EditSeqPtr    esp;
  Int2          frame;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  if (gcp == NULL) return TRUE;
  esp = (EditSeqPtr) gcp->userdata;
  if (esp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return TRUE;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  sip = SeqLocId (slp);
  if (sip == NULL) return TRUE;
  bsp = BioseqFind (sip);
  if (bsp == NULL || bsp != esp->extendedthis) return TRUE;
  if (esp->adjustframe) {
    if (GetOffsetInBioseq (slp, bsp, SEQLOC_START) != 0) return TRUE;
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp == NULL) return TRUE;
    frame = crp->frame;
    if (frame == 0)
      frame = 1;
    if (esp->addOrTrimBool == ADD_TO_END)
      {
	frame--;
	frame += esp->frameshift;
	crp->frame = (frame % 3) + 1;
      }
    else if (esp->addOrTrimBool == TRIM_FROM_END)
      {
	frame = ABS(frame - esp->frameshift);
	crp->frame = 3 - (frame % 3);
      }
  } else {
    if (GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) != bsp->length - 1)
      return TRUE;
  }
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return TRUE;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return TRUE;
  if (bsp->repr != Seq_repr_raw) return TRUE;
  if (bsp->mol != Seq_mol_aa) return TRUE;
  bestprot = FindBestProtein (esp->input_entityID, sfp->product);
  bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
  if (bs == NULL) return TRUE;
  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot == NULL) return TRUE;
  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    /*
    if (prot [0] == '-') {
       ptr++;
    }
    */
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
    bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
    bsp->seq_data = (SeqDataPtr) bs;
    bsp->length = BSLen (bs);
    bsp->seq_data_type = Seq_code_ncbieaa;
  }
  if (bestprot == NULL) return TRUE;
  sep = SeqMgrGetSeqEntryForData (bsp);
  bestprot->location = SeqLocFree (bestprot->location);
  bestprot->location = CreateWholeInterval (sep);
  SetSeqLocPartial (bestprot->location, partial5, partial3);
  return TRUE;
}

static void FixAndRetranslateCDSs (BioseqPtr bsp, SeqEntryPtr sep,
                                   EditSeqPtr esp, Boolean adjustframe)

{
  GatherScope  gs;

  if (esp->input_entityID == 0 || esp->sep == NULL) return;
  esp->adjustframe = adjustframe;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  gs.scope = sep;
  MemSet ((Pointer)(gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_BIOSEQ] = FALSE;
  gs.ignore [OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore [OBJ_SEQFEAT] = FALSE;
  gs.ignore [OBJ_SEQANNOT] = FALSE;
  GatherEntity (esp->input_entityID, (Pointer) esp, FixACDSFunc, &gs);
  gs.target = SeqLocFree (gs.target);
}

static ValNodePtr CollectAndExtendSingleBaseFeatures (BioseqPtr bsp, Int2 whichend, Int4 len)

{
  SeqMgrFeatContext  context;
  ValNodePtr         head = NULL;
  SeqFeatPtr         sfp;
  SeqLocPtr          slp;
  SeqPntPtr          spp;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  while (sfp != NULL) {
    if (whichend == 1 && context.numivals == 1 && context.right == 0) {
      slp = sfp->location;
      if (slp != NULL && slp->choice == SEQLOC_PNT && slp->next == NULL) {
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL && spp->point == 0) {
          spp->point = 1;
          ValNodeAddPointer (&head, 1, (Pointer) sfp);
        }
      }
    } else if (whichend == 2 && context.numivals == 1 && context.left == bsp->length - 1) {
      slp = sfp->location;
      if (slp != NULL && slp->choice == SEQLOC_PNT && slp->next == NULL) {
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL && spp->point == bsp->length - 1) {
          spp->point = bsp->length - 2;
          ValNodeAddPointer (&head, 2, (Pointer) sfp);
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }

  return head;
}

static void ReadjustSingleBaseFeatures (ValNodePtr head, BioseqPtr bsp, Int2 whichend, Int4 len)

{
  SeqFeatPtr  sfp;
  SeqIntPtr   sintp;
  SeqLocPtr   slp;
  SeqPntPtr   spp;
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL) {
      slp = sfp->location;
      if (slp != NULL && slp->choice == SEQLOC_PNT && slp->next == NULL) {
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          if (whichend == 1) {
            sintp = SeqIntNew ();
            if (sintp != NULL) {
              sintp->from = 0;
              sintp->to = spp->point - 1;
              sintp->strand = spp->strand;
              sintp->id = spp->id;
              spp->id = NULL;
              sintp->if_from = spp->fuzz;
              spp->fuzz = NULL;
              slp->choice = SEQLOC_INT;
              slp->data.ptrvalue = (Pointer) sintp;
              SeqPntFree (spp);
            }
          } else if (whichend == 2) {
            sintp = SeqIntNew ();
            if (sintp != NULL) {
              sintp->from = spp->point + 1;
              sintp->to = spp->point + 1 + len;
              sintp->strand = spp->strand;
              sintp->id = spp->id;
              spp->id = NULL;
              sintp->if_to = spp->fuzz;
              spp->fuzz = NULL;
              slp->choice = SEQLOC_INT;
              slp->data.ptrvalue = (Pointer) sintp;
              SeqPntFree (spp);
            }
          }
        }
      }
    }
  }
}

static void TrimFromSequenceEnd (EditSeqPtr  esp,
				 SeqEntryPtr sep,
				 BioseqPtr   bsp)
{
  CharPtr currSeqStr;
  Int2    length;
  Int4    pos;
  Int4    trim_length;

  /* Get the current sequence string */

  currSeqStr = GetSequenceByBsp (bsp);

  /* Trim from either the 5' end (i.e., */
  /* the beginning of the string...     */

  if (esp->endval == EDIT_FIVE_PRIME)
  {
    /* Find end point */

    if (esp->trimByBool == TRIM_BY_SEQUENCE)
	  {
	    length = StringLen (esp->seqstr);
	    if (StringNICmp (esp->seqstr, currSeqStr, length) != 0)
	      return;
	    pos = length - 1;
      trim_length = length;
	  }
    else  if (esp->trimByBool == TRIM_BY_COUNT) 
    {
	    pos = esp->trimCount - 1;
      trim_length = esp->trimCount;
    }
    else
	    return;
      
    /* Trim from beginning of string to end point */
    
    esp->frameshift = trim_length;
    BioseqDelete (bsp->id, 0, pos, TRUE, FALSE);
    esp->extendedthis = bsp;
    FixAndRetranslateCDSs (bsp, sep, esp, TRUE);

    /* trim quality scores */
    TrimQualityScores (bsp, trim_length, TRUE);
  }
  
  /* .. or the 3' end (i.e., the */
  /* end of the string.          */
  
  else if (esp->endval == EDIT_THREE_PRIME)
  {
    /* Find trim point */

    if (esp->trimByBool == TRIM_BY_SEQUENCE)
	  {
	    length = StringLen (esp->seqstr);
	    pos = bsp->length - length;
	    if (StringICmp (esp->seqstr, &currSeqStr[pos]) != 0)
	      return;
      trim_length = length;
	  }
    else  if (esp->trimByBool == TRIM_BY_COUNT)
    {
	    pos = bsp->length - esp->trimCount;
      trim_length = esp->trimCount;
    }
    else
	    return;
      
    /* Trim from there to end of string */
    
    BioseqDelete (bsp->id, pos, bsp->length - 1, TRUE, FALSE);
    esp->extendedthis = bsp;
    FixAndRetranslateCDSs (bsp, sep, esp, FALSE);

    /* trim quality scores */
    TrimQualityScores (bsp, trim_length, FALSE);
  }
}

static void 
AddToSequenceEnd 
(EditSeqPtr  esp,
 SeqEntryPtr sep,
 BioseqPtr   bsp,
 LogInfoPtr lip)
{
  ValNodePtr    head;
  Int4          len;
  Int4          pos;
  Uint1         residue;
  SeqPortPtr    spp;
  CharPtr       str;
  Char          terminal [2];

  pos = 0;
  if (esp->endval == 2) {
    pos = bsp->length;
  }
  if (esp->extendflag) {
    esp->frameshift = 0;
    terminal [0] = '\0';
    terminal [1] = '\0';
    residue = 0;
    if (esp->endval == 2) {
      spp = SeqPortNew (bsp, bsp->length - 1, -1, 0, Seq_code_iupacna);
    } else {
      spp = SeqPortNew (bsp, 0, -1, 0, Seq_code_iupacna);
    }
    if (spp != NULL) {
      residue = SeqPortGetResidue (spp);
      if (IS_residue (residue)) {
        terminal [0] = TO_LOWER ((Char) residue);
      }
    }
    SeqPortFree (spp);
    str = MemNew ((size_t) (StringLen (esp->seqstr) + 4));
    if (str != NULL) {
      head = NULL;
      if (esp->endval == 2) {
        esp->extendedthis = bsp;
        StringCpy (str, terminal);
        StringCat (str, esp->seqstr);
        len = StringLen (esp->seqstr);
        pos = bsp->length - 1;
        head = CollectAndExtendSingleBaseFeatures (bsp, 2, len);
        insertchar (str, pos, bsp->id, bsp->mol, FALSE);
        BioseqDelete (bsp->id, bsp->length - 1, bsp->length - 1,
		      TRUE, FALSE);
        ReadjustSingleBaseFeatures (head, bsp, 2, len);
        FixAndRetranslateCDSs (bsp, sep, esp, FALSE);
      } else {
        esp->frameshift = (Int2) StringLen (esp->seqstr);
        esp->extendedthis = bsp;
        StringCpy (str, esp->seqstr);
        StringCat (str, terminal);
        len = StringLen (esp->seqstr);
        pos = 1;
        head = CollectAndExtendSingleBaseFeatures (bsp, 1, len);
        insertchar (str, pos, bsp->id, bsp->mol, FALSE);
        BioseqDelete (bsp->id, 0, 0, TRUE, FALSE);
        ReadjustSingleBaseFeatures (head, bsp, 1, len);
        FixAndRetranslateCDSs (bsp, sep, esp, TRUE);
      }
      ValNodeFree (head);
      if (lip != NULL) {
        RemoveQualityScores (bsp, lip->fp, &(lip->data_in_log));
      }
    }
    MemFree (str);
  } else {
    insertchar (esp->seqstr, pos, bsp->id, bsp->mol, FALSE);
    if (lip != NULL) {
      RemoveQualityScores (bsp, lip->fp, &(lip->data_in_log));
    }
  }
}


static void DoEditSeqEndsProc (ButtoN b)

{
  Char        ch;
  EditSeqPtr  esp;
  CharPtr     p, q;
  CharPtr     tempStr;
  ValNodePtr  sip_list, vnp;
  SeqIdPtr    sip;
  BioseqPtr   bsp;
  SeqEntryPtr sep;
  Boolean     add_cit_subs = FALSE;
  LogInfoPtr  lip;

  esp = (EditSeqPtr) GetObjectExtra (b);
  if (esp == NULL) {
    Remove (ParentWindow (b));
    return;
  }
  sip_list = GetSelectedSequenceList (esp->nuc_sequence_list_ctrl);
  if (sip_list == NULL)
  {
    Message (MSG_ERROR, "You have not specified any sequences to edit!");
    return;
  }

  Hide (esp->form);
  Update ();
  esp->seqstr = SaveStringFromText  (esp->seq);
  p = esp->seqstr;
  if (p != NULL) {
    /* remove any non-sequence characters */
    q = p;
    ch = *p;
    while (ch != '\0') {
      if (IS_ALPHA (ch)) {
        *q = ch;
        q++;
      }
      p++;
      ch = *p;
    }
    *q = '\0';
  }
  esp->genestr       = SaveStringFromText  (esp->genename);
  esp->endval        = GetValue (esp->whichend);
  esp->extendflag    = GetStatus (esp->extendfeat);
  esp->addOrTrimBool = GetValue (esp->addOrTrim);
  esp->trimByBool    = GetValue (esp->trimBy);
  tempStr            = SaveStringFromText (esp->trimCountText);
  if (tempStr != NULL)
    esp->trimCount   = atoi (tempStr);
  else
    esp->trimCount   = 0;
  
  add_cit_subs = GetStatus (esp->addCitSub);

  lip = OpenLog ("Quality Scores Affected");
  for (vnp = sip_list; vnp != NULL; vnp = vnp->next)
  {
    sip = (SeqIdPtr) vnp->data.ptrvalue;
    bsp = BioseqFind (sip);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (bsp != NULL && sep != NULL && EditSeqEntryHasGene (bsp, sep, esp))
    {
      if (esp->addOrTrimBool == 1)
        AddToSequenceEnd (esp, sep, bsp, lip);
      else
        TrimFromSequenceEnd (esp, sep, bsp);
      if (add_cit_subs)
      {
        AddCitSubToUpdatedSequence (bsp, esp->input_entityID, kSubmitterUpdateText);
      }
    }
  }
  CloseLog (lip);
  lip = FreeLog (lip);

  MemFree (esp->seqstr);
  MemFree (esp->genestr);
  ObjMgrSetDirtyFlag (esp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, esp->input_entityID, 0, 0);
  Remove (esp->form);
}

static void EditSeqMessageProc (ForM f, Int2 mssg)

{
  EditSeqPtr  esp;

  esp = (EditSeqPtr) GetObjectExtra (f);
  if (esp != NULL) {
    if (esp->appmessage != NULL) {
      esp->appmessage (f, mssg);
    }
  }
}

static void ChangeAddOrTrim_Callback (GrouP g)

{
  EditSeqPtr  esp;
  Int2        val;
  Int2        trimByState;

  val = GetValue (g);
  esp = GetObjectExtra (g);

  switch (val) {
    case ADD_TO_END :
      SafeDisable (esp->trimBy);
      SafeEnable (esp->extendfeat);
      SafeDisable (esp->trimCountText);
      SafeEnable (esp->seq);
      break;
    case TRIM_FROM_END :
      SafeDisable (esp->extendfeat);
      SafeEnable (esp->trimBy);
      trimByState = GetValue (esp->trimBy);
      switch (trimByState) {
        case TRIM_BY_SEQUENCE :
	  SafeDisable (esp->trimCountText);
	  SafeEnable (esp->seq);
	  break;
        case TRIM_BY_COUNT :
	  SafeEnable (esp->trimCountText);
	  SafeDisable (esp->seq);
	  break;
        default :
	  break;
      }
      break;
    default :
      break;
  }
}

static void ChangeTrimBy_Callback (GrouP g)

{
  EditSeqPtr  esp;
  Int2        val;

  val = GetValue (g);
  esp = GetObjectExtra (g);

  switch (val) {
    case TRIM_BY_SEQUENCE :
      SafeDisable (esp->trimCountText);
      SafeEnable (esp->seq);
      break;
    case TRIM_BY_COUNT :
      SafeEnable (esp->trimCountText);
      SafeDisable (esp->seq);
      break;
    default :
      break;
  }
}

static void SelectAllSequencesForExtend (ButtoN b)
{
  EditSeqPtr    esp;

  esp = (EditSeqPtr) GetObjectExtra (b);
  if (esp == NULL)
  {
    return;
  }
  SelectAllSequencesInListCtrl (esp->nuc_sequence_list_ctrl);  
}

static void UnSelectAllSequencesForExtend (ButtoN b)
{
  EditSeqPtr    esp;

  esp = (EditSeqPtr) GetObjectExtra (b);
  if (esp == NULL)
  {
    return;
  }
  UnSelectAllSequencesInListCtrl (esp->nuc_sequence_list_ctrl);  
}

extern void EditSeqEndsProc (IteM i);

extern void EditSeqEndsProc (IteM i)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  GrouP              c;
  EditSeqPtr         esp;
  GrouP              g;
  GrouP              h;
  GrouP              k;
  GrouP              p;
  GrouP              q;
  GrouP              r;
  GrouP              s;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  esp = (EditSeqPtr) MemNew (sizeof (EditSeqEnds));
  if (esp == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Edit Sequence Ends", StdCloseWindowProc);
  SetObjectExtra (w, esp, StdCleanupFormProc);
  esp->form = (ForM) w;
  esp->formmessage = EditSeqMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    esp->appmessage = sepp->handleMessages;
  }

  esp->input_entityID = bfp->input_entityID;
  esp->input_itemID = bfp->input_itemID;
  esp->input_itemtype = bfp->input_itemtype;

  esp->sep = sep;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);

  StaticPrompt (g, "End", 0, stdLineHeight, programFont, 'l');
  esp->whichend = HiddenGroup (g, 2, 0, NULL);
  RadioButton (esp->whichend, "5'");
  RadioButton (esp->whichend, "3'");
  SetValue (esp->whichend, 1);

  esp->addOrTrim = HiddenGroup (h, 2, 0, ChangeAddOrTrim_Callback);
  SetObjectExtra (esp->addOrTrim, esp, NULL);
  RadioButton (esp->addOrTrim, "Add to end");
  RadioButton (esp->addOrTrim, "Trim from end");
  SetValue (esp->addOrTrim, 1);

  k = HiddenGroup (h, 0, -2, NULL);
  StaticPrompt (k, "Sequence", 0, 0, programFont, 'l');
  esp->seq = ScrollText (k, 25, 5, programFont, TRUE, NULL);

  q = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (q, "Optional gene constraint", 0, dialogTextHeight,
		programFont, 'l');
  esp->genename = DialogText (q, "", 14, NULL);

  esp->extendfeat = CheckBox (h, "Extend features", NULL);

  esp->trimBy = HiddenGroup (h, 2, 0, ChangeTrimBy_Callback);
  SetObjectExtra (esp->trimBy, esp, NULL);
  RadioButton (esp->trimBy, "Trim by sequence");
  RadioButton (esp->trimBy, "Trim by count");
  SetValue (esp->trimBy, 1);
  SafeDisable (esp->trimBy);

  p = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (p, "Trim count", 0, dialogTextHeight, programFont, 'l');
  esp->trimCountText = DialogText (p, "", 5, NULL);
  SafeDisable (esp->trimCountText);
  
  s = NormalGroup (h, -1, 0, "Choose Sequences To Edit", programFont, NULL);
  esp->nuc_sequence_list_ctrl = MakeSequenceListControl (s, sep, NULL, NULL, TRUE, FALSE);
  
  r = HiddenGroup (s, 2, 0, NULL);
  b = PushButton (r, "Select All", SelectAllSequencesForExtend);
  SetObjectExtra (b, esp, NULL);
  b = PushButton (r, "Unselect All", UnSelectAllSequencesForExtend);
  SetObjectExtra (b, esp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) esp->nuc_sequence_list_ctrl,
                (HANDLE) r, NULL);

  esp->addCitSub = CheckBox (h, "Add Cit Subs to edited sequences", NULL);
  
  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoEditSeqEndsProc);
  SetObjectExtra (b, esp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) esp->extendfeat,
		(HANDLE) esp->addOrTrim, (HANDLE) esp->trimBy,
		(HANDLE) k, (HANDLE) q, (HANDLE) p, 
		(HANDLE) s, (HANDLE) esp->addCitSub, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (esp->seq);
  Update ();
}

static void SetAlignmentDim (SeqAlignPtr salp)
{
  AMAlignIndex2Ptr amaip;
  DenseSegPtr      dsp;
  
  if (salp == NULL || salp->dim > 0 || salp->saip == NULL) return;
  
  if (salp->saip->indextype == INDEX_PARENT)
  {
    amaip = (AMAlignIndex2Ptr)(salp->saip);
    salp->dim = amaip->sharedaln->dim;
  }
  else if (salp->saip->indextype == INDEX_CHILD)
  {
    dsp = (DenseSegPtr)(salp->segs);
  	salp->dim = dsp->dim;
  }
}  

static void IndexAlignmentSet (SeqAlignPtr salp)
{
  SeqAlignPtr tmp_salp, next_salp;
  
  if (salp == NULL || salp->saip != NULL) return;
  
  if (salp->next != NULL && salp->dim > 2)
  {
  	for (tmp_salp = salp; tmp_salp != NULL; tmp_salp = tmp_salp->next)
  	{
  	  next_salp = tmp_salp->next;
  	  tmp_salp->next = NULL;
      if (tmp_salp->segtype == SAS_DENSEG  &&  tmp_salp->next == NULL) {
        AlnMgr2IndexSingleChildSeqAlign(tmp_salp);
      } else {
        AlnMgr2IndexSeqAlign(tmp_salp);
      }  		
      SetAlignmentDim (tmp_salp);
      tmp_salp->next = next_salp;
  	}
  }
  else
  {
    if (salp->segtype == SAS_DENSEG  &&  salp->next == NULL) {
      AlnMgr2IndexSingleChildSeqAlign(salp);
    } else {
      AlnMgr2IndexSeqAlign(salp);
    }  
    SetAlignmentDim (salp);		
  }	
}

static void WriteSeqEntryAlignmentToFile (SeqEntryPtr sep, FILE *fp, Boolean Interleave)
{
  BioseqSetPtr bssp;
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp = NULL;

  if (sep == NULL || ! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  for (sap = bssp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 2) {
      salp = SeqAlignListDup((SeqAlignPtr) sap->data);
      IndexAlignmentSet (salp);

      if (Interleave) {
        if (salp->next != NULL)
        {
          Message (MSG_ERROR, "Unable to write segmented alignments as interleave");
          return;
        }
        WriteAlignmentInterleaveToFile (salp, fp, 40, FALSE);
      } else {
        WriteAlignmentContiguousToFile (salp, fp, 40, FALSE);
      }
      SeqAlignFree (salp);
      salp = NULL;
    }
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    WriteSeqEntryAlignmentToFile (sep, fp, Interleave);
  }
}

static void ExportAlignment (IteM i, Boolean Interleave)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;
  Char        path [PATH_MAX];
  FILE *      fp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if (! GetOutputFileName (path, sizeof (path), "")) return;
  if (! StringHasNoText (path)) {
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      WatchCursor ();
      Update();
      WriteSeqEntryAlignmentToFile (sep, fp, Interleave);
      ArrowCursor ();
      Update();
      FileClose (fp);
    } else {
      Message (MSG_ERROR, "Unable to open file");
    }
  }
}

extern void ExportAlignmentInterleave (IteM i)
{
  ExportAlignment (i, TRUE);
}

extern void ExportAlignmentContiguous (IteM i)
{
  ExportAlignment (i, FALSE);
}

static Int4 FindFeaturePos (SeqIdPtr sip, SeqAlignPtr salp, Int4 pos)
{
  Int4        aln_row;
  Int4        new_pos;
  
  if (sip == NULL || salp == NULL) return pos;
  aln_row = AlnMgr2GetFirstNForSip(salp, sip);
  if (aln_row > 0)
  {
  	new_pos = AlnMgr2MapSeqAlignToBioseq (salp, pos, aln_row);
  	if (new_pos < 0)
  	{
  	  return pos;
  	}
  	else
  	{
  	  return new_pos;
  	}
  }
  return pos;
}

static void FixFeatureIntervalSeqLoc (SeqLocPtr slp, SeqAlignPtr salp)
{
  SeqIntPtr     sintp;
  SeqBondPtr    sbp;
  SeqPntPtr     spp;
  SeqIdPtr      tmpsip;
  SeqLocPtr     this_slp;
 
  if (slp == NULL || slp->data.ptrvalue == NULL || salp == NULL) return;

  tmpsip = SeqIdPtrFromSeqAlign (salp);

  switch (slp->choice)
  {
    case SEQLOC_INT:
      sintp = slp->data.ptrvalue;
      sintp->from = FindFeaturePos (sintp->id, salp, sintp->from);
      sintp->to = FindFeaturePos (sintp->id, salp, sintp->to);
      break;
  	case SEQLOC_PNT:
  	  spp = slp->data.ptrvalue;
  	  spp->point = FindFeaturePos (spp->id, salp, spp->point);
  	  break;
    case SEQLOC_BOND:   /* bond -- 2 seqs */
      sbp = (SeqBondPtr)(slp->data.ptrvalue);
      spp = sbp->a;
  	  spp->point = FindFeaturePos (spp->id, salp, spp->point);
      spp = sbp->b;
  	  spp->point = FindFeaturePos (spp->id, salp, spp->point);
      break;
    case SEQLOC_MIX:    /* mix -- more than one seq */
    case SEQLOC_EQUIV:    /* equiv -- ditto */
    case SEQLOC_PACKED_INT:    /* packed int */
      for (this_slp = slp->data.ptrvalue; this_slp != NULL; this_slp = this_slp->next)
      {
      	FixFeatureIntervalSeqLoc (this_slp, salp);
      }
      break;
    default:
      break;	
  }
}

typedef struct featurefixdata 
{
  Uint2 entityID;
  SeqAlignPtr salp;	
} FeatureFixData, PNTR FeatureFixPtr;

static void FixFeatureIntervalCallback (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr     bsp;
  FeatureFixPtr ffp;
 
  if (sfp == NULL || userdata == NULL) return;  
  
  ffp = (FeatureFixPtr) userdata;
 
  FixFeatureIntervalSeqLoc (sfp->location, ffp->salp);
  if (sfp->idx.subtype == FEATDEF_CDS)
  {
    bsp = BioseqFindFromSeqLoc (sfp->location);
  	SeqEdTranslateOneCDS (sfp, bsp, ffp->entityID, Sequin_GlobalAlign2Seq);
  }
}

extern void FixFeatureIntervals (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;
  FeatureFixData ffd;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ffd.entityID = bfp->input_entityID;
  ffd.salp = SeqAlignListDup((SeqAlignPtr) FindSeqAlignInSeqEntry (sep, OBJ_SEQALIGN));
  if (ffd.salp == NULL)
  {
    Message (MSG_ERROR, "No alignment present - cannot remap intervals");
    return;
  }
      
  if (ffd.salp->segtype == SAS_DENSEG  &&  ffd.salp->next == NULL) 
  {
    AlnMgr2IndexSingleChildSeqAlign(ffd.salp);
  } else {
    AlnMgr2IndexSeqAlign(ffd.salp);
  }

  VisitFeaturesInSep (sep, (Pointer) &ffd, FixFeatureIntervalCallback);
  SeqAlignFree (ffd.salp);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static SeqLocPtr 
CreateProteinLoc 
(SeqLocPtr  nuc_loc,
 SeqFeatPtr top_cds,
 Int4       offset,
 Int4       cds_frame, 
 BioseqPtr  prot_bsp,
 LogInfoPtr lip)
{
  SeqLocPtr    prot_loc = NULL;
  SeqIntPtr    sintp_nuc, sintp_prot, sintp_top;
  Boolean      partial5, partial3;
  ByteStorePtr bs;
  CharPtr      prot_str;
  Boolean      ends_with_stop = FALSE;
  Int4         prot_len = 0;
  Char         ch;
  Int4         from_remainder = 0, to_remainder = 0;
  
  if (nuc_loc == NULL || top_cds == NULL || prot_bsp == NULL
      || top_cds->location == NULL
      || top_cds->location->choice != SEQLOC_INT
      || top_cds->location->data.ptrvalue == NULL)
  {
    return NULL;
  }
  if (nuc_loc->choice != SEQLOC_INT || nuc_loc->data.ptrvalue == NULL)
  {
    Message (MSG_ERROR, "Unable to translate locations that are not simple intervals");
  }
  else
  {
    sintp_top = (SeqIntPtr) (top_cds->location->data.ptrvalue);
    sintp_nuc = (SeqIntPtr)(nuc_loc->data.ptrvalue);
    
    if (sintp_nuc->from < sintp_top->from)
    {
      from_remainder = -1;
    }
    else if (cds_frame == 1)
    {
      from_remainder = (sintp_nuc->from - offset) % 3;
    }
    else if (cds_frame == 2)
    {
      from_remainder = (sintp_nuc->from - offset + 1) % 3;
    }
    else if (cds_frame == 3)
    {
      from_remainder = (sintp_nuc->from - offset + 2) % 3;
    }
    
    if (sintp_nuc->to > sintp_top->to)
    {
      to_remainder = -1;
    }
    else if (sintp_nuc->to < sintp_top->to)
    {
      to_remainder = (sintp_nuc->to + 1 - offset) % 3;
    }

    if (from_remainder != 0 || to_remainder != 0)    
    {
      if (lip != NULL && lip->fp != NULL)
      {
        fprintf (lip->fp, "Invalid coordinates for mat_peptide conversion\n");
        lip->data_in_log = TRUE;
      }
      else
      {
        Message (MSG_ERROR, "Invalid coordinates for mat_peptide conversion");
      }
    }
    else
    {  
      bs = ProteinFromCdRegionEx (top_cds, TRUE, FALSE);
      if (bs != NULL)
      {
        prot_str = BSMerge (bs, NULL);
        prot_len = StringLen (prot_str);
        ch = prot_str [prot_len - 1];
        if (ch == '*')
        {
          ends_with_stop = TRUE;
        }
        MemFree (prot_str);
        bs = BSFree (bs);
      }

      sintp_prot = SeqIntNew ();
      if (sintp_prot != NULL)
      {
        CheckSeqLocForPartial (nuc_loc, &partial5, &partial3);
        sintp_prot->id = SeqIdDup (SeqIdFindBest(prot_bsp->id, 0));
        sintp_prot->from = (sintp_nuc->from - offset) / 3;
        sintp_prot->to = (sintp_nuc->to - offset) / 3;
        if (sintp_prot->to == prot_len - 1)
        {
          if (ends_with_stop)
          {
            sintp_prot->to --;
          }
          else
          {
            partial3 = TRUE;
          }
        }
        ValNodeAddPointer (&prot_loc, SEQLOC_INT, sintp_prot);
        SetSeqLocPartial (prot_loc, partial5, partial3);
      }
    }
  }
  return prot_loc;
}

static CharPtr GetCDSProteinDesc (SeqFeatPtr cds)
{
  BioseqPtr         bsp;
  SeqFeatPtr        prot_feat;
  SeqMgrFeatContext context;
  ProtRefPtr        prp;
  
  if (cds == NULL || cds->product == NULL)
  {
    return NULL;
  }
  
  bsp = BioseqFindFromSeqLoc (cds->product);
  prot_feat = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &context);
  if (prot_feat != NULL && prot_feat->data.value.ptrvalue != NULL)
  {
    prp = (ProtRefPtr) prot_feat->data.value.ptrvalue;
    return StringSave (prp->desc);
  }
  else
  {
    return NULL;
  }
}

extern void ConvertInnerCDSsToMatPeptidesCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr        sfp = NULL;
  ValNodePtr        top_level_cds_list = NULL;
  ValNodePtr        top_level_cds_offset_list = NULL;
  SeqMgrFeatContext context, gene_context;
  ValNodePtr        vnp, offset_vnp;
  SeqFeatPtr        top_cds;
  BioseqPtr         prot_bsp;
  SeqLocPtr         prot_loc;
  Int4              offset = 0, mat_cds_frame;
  SeqFeatPtr        new_sfp;
  ProtRefPtr        prp;
  SeqFeatPtr        gene;
  CdRegionPtr       crp;
  LogInfoPtr        lip;
  
  if (bsp == NULL || ! ISA_na (bsp->mol)) return;
  lip = (LogInfoPtr) userdata;
  
  sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL)
  {
    top_cds = NULL;
    offset = 0;
    for (vnp = top_level_cds_list, offset_vnp = top_level_cds_offset_list;
         vnp != NULL && offset_vnp != NULL && top_cds == NULL; 
         vnp = vnp->next, offset_vnp = offset_vnp->next)
    {
      top_cds = (SeqFeatPtr) vnp->data.ptrvalue;
      if (top_cds != NULL)
      {
        if (SeqLocCompare (top_cds->location, sfp->location) == SLC_B_IN_A)
        {
          offset = offset_vnp->data.intvalue;
        }
        else
        {
          top_cds = NULL;
        }
      }
    }
    /* Only handling the simplest cases of SeqInt locs inside other SeqInt locs */
    if (top_cds == NULL)
    {
      if (sfp->location != NULL 
          && sfp->location->choice == SEQLOC_INT
          && sfp->location->data.ptrvalue != NULL
          && sfp->data.value.ptrvalue != NULL)
      {
        /* add to list of top level CDSs */
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        offset = context.left;
        if (crp->frame == 2)
        {
          offset += 1;
        }
        else if (crp->frame == 3)
        {
          offset += 2;
        }
        ValNodeAddPointer (&top_level_cds_list, 0, sfp);
        ValNodeAddInt (&top_level_cds_offset_list, 0, offset);
      }
    }
    else
    {
      prot_bsp = BioseqFindFromSeqLoc (top_cds->product);
      if (prot_bsp != NULL)
      {
        /* convert location by subtracting top_cds start from sfp locations and
         * divide all locations by 3
         */
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        mat_cds_frame = 1;
        if (crp->frame == 2)
        {
          mat_cds_frame = 2;
        }
        else if (crp->frame == 3)
        {
          mat_cds_frame = 3;
        }
         
        prot_loc = CreateProteinLoc (sfp->location, top_cds, 
                                     offset,
                                     mat_cds_frame, 
                                     prot_bsp, lip);
        if (prot_loc != NULL)
        {
          /* Create new feature on prot_bsp */
          new_sfp = CreateNewFeatureOnBioseq (prot_bsp, SEQFEAT_PROT, prot_loc);
          if (new_sfp != NULL)
          {
            prp = ProtRefNew ();
            if (prp != NULL)
            {
              ValNodeCopyStr (&(prp->name), 0, context.label);
              prp->processed = 2;
              prp->desc = GetCDSProteinDesc (sfp);
            }
            new_sfp->data.value.ptrvalue = prp;
            
            /* mark old product for deletion */
            prot_bsp = BioseqFindFromSeqLoc (sfp->product);
            if (prot_bsp != NULL)
            {
              prot_bsp->idx.deleteme = 1;
            }
         
            /* Mark old feature for deletion */
            sfp->idx.deleteme = 1;
            gene = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
            if (gene != NULL && SeqLocCompare (gene->location, sfp->location) == SLC_A_EQ_B)
            {
              gene->idx.deleteme = 1;
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
  ValNodeFree (top_level_cds_list);
}



extern void ConvertInnerCDSsToProteinFeatures (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;
  LogInfoPtr     lip;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  lip = OpenLog ("CDS to Mat Peptide Conversion");
  VisitBioseqsInSep (sep, lip, ConvertInnerCDSsToMatPeptidesCallback);
  CloseLog (lip);
  lip = FreeLog (lip);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();  
}

typedef struct objstringdata 
{
  StringConstraintXPtr scp;
  Boolean found;	
} ObjStringData, PNTR ObjStringPtr;

static void LIBCALLBACK AsnWriteStringConstraintCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr        pchSource;
  ObjStringPtr   osp;

  osp = (ObjStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) 
  {
	  pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
    if (DoesStringMatchConstraintX (pchSource, osp->scp))
    {
      osp->found = TRUE;
    }
  }
}

static Boolean DoesBioseqMatchStringConstraint (BioseqPtr bsp, StringConstraintXPtr scp)

{
  ObjMgrPtr         omp;
  ObjMgrTypePtr     omtp;
  AsnExpOptPtr      aeop;
  AsnIoPtr          aip;
  ObjStringData     osd;

  omp = ObjMgrGet ();
  if (omp == NULL) return FALSE;
  omtp = ObjMgrTypeFind (omp, OBJ_BIOSEQ, NULL, NULL);
  if (omtp == NULL) return FALSE;
  
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteStringConstraintCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &osd;
  }
  osd.scp = scp;

  osd.found = FALSE;
  (omtp->asnwrite) (bsp, aip, NULL);
  AsnIoClose (aip);   
  
  if (scp != NULL && scp->not_present)
  {
    osd.found = ! osd.found;
  }
  
  return osd.found;
}


extern Boolean DoBioseqFeaturesMatchSequenceConstraintX (BioseqPtr bsp, ValNodePtr feat_list, StringConstraintXPtr scp)
{
  AsnExpOptPtr            aeop;
  AsnIoPtr                aip;
  ObjStringData           osd;
  SeqFeatPtr              sfp;
  ObjMgrPtr               omp;
  ObjMgrTypePtr           omtp;
  SeqMgrFeatContext       fcontext;
  ValNodePtr              vnp;
  
  if (bsp == NULL) return FALSE;
  if (scp == NULL) return TRUE; 
  omp = ObjMgrGet ();
  if (omp == NULL) return FALSE;
  omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp == NULL) return FALSE;

  aip = AsnIoNullOpen ();
  osd.scp = scp;
  
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteStringConstraintCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &osd;
  }

  for (vnp = feat_list; vnp != NULL; vnp = vnp->next)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, vnp->choice, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, vnp->choice, &fcontext))
    {
      osd.found = FALSE;
      (omtp->asnwrite) (sfp, aip, NULL);
      if (osd.found)
      {
        if (scp->not_present)
        {
          AsnIoClose (aip);
          return FALSE;
        }
        else
        {
          AsnIoClose (aip);
          return TRUE;
        }
      }
    }
  }
  AsnIoClose (aip);
  if (scp->not_present)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}
 

typedef struct keywordform 
{
  FORM_MESSAGE_BLOCK
  DialoG string_src_dlg;
  DialoG string_constraint_dlg;
  TexT   keyword_txt;
  
  ParseFieldPtr pfp;
  FilterSetPtr  fsp;
  CharPtr       keyword;
} KeywordFormData, PNTR KeywordFormPtr;


static void ApplyKeywordCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqEntryPtr             sep;
  ValNodePtr              vnp;
  KeywordFormPtr scfp;
  GBBlockPtr              gbp;
  GetSamplePtr            gsp;
  Boolean                 ok_to_add = TRUE;
  
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL)
  {
    return;
  }
  
  scfp = (KeywordFormPtr) userdata;
  
  if (scfp->pfp != NULL && scfp->fsp != NULL)
  {
    gsp = GetSampleForSeqEntry (sep, bsp->idx.entityID, scfp->pfp, scfp->fsp);
    if (gsp == NULL || gsp->num_found == 0)
    {
      ok_to_add = FALSE;
    }
    gsp = GetSampleFree (gsp);
  }
  
  if (!ok_to_add)
  {
    return;
  }
      
	vnp = GetDescrOnSeqEntry (sep, Seq_descr_genbank);
	if (vnp == NULL) {
		vnp = NewDescrOnSeqEntry (sep, Seq_descr_genbank);
		if (vnp != NULL) {
			vnp->data.ptrvalue = (Pointer) GBBlockNew ();
		}
	}
	if (vnp == NULL) return;
	gbp = (GBBlockPtr) vnp->data.ptrvalue;
	if (gbp == NULL)
	{
	  gbp = GBBlockNew ();
	  vnp->data.ptrvalue = gbp;
	}
	if (gbp == NULL) return;
  	
	for (vnp = gbp->keywords; vnp; vnp = vnp->next) {
		if (StringCmp((CharPtr)vnp->data.ptrvalue, scfp->keyword) == 0) {
			return;
		}
	}
	ValNodeAddPointer (&(gbp->keywords), 0, StringSave (scfp->keyword));
}

static void RemoveKeywordCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqEntryPtr             sep;
  ValNodePtr              vnp, prev_keyword, next_keyword;
  KeywordFormPtr scfp;
  GBBlockPtr              gbp;
  GetSamplePtr            gsp;
  Boolean                 ok_to_remove = TRUE;
  
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL)
  {
    return;
  }
  
  scfp = (KeywordFormPtr) userdata;
  
  if (scfp->pfp != NULL && scfp->fsp != NULL)
  {
    gsp = GetSampleForSeqEntry (sep, bsp->idx.entityID, scfp->pfp, scfp->fsp);
    if (gsp == NULL || gsp->num_found == 0)
    {
      ok_to_remove = FALSE;
    }
    gsp = GetSampleFree (gsp);
  }
  
  if (!ok_to_remove)
  {
    return;
  }
      
	vnp = GetDescrOnSeqEntry (sep, Seq_descr_genbank);
  /* no GenBank descriptor, no keywords to remove */
	if (vnp == NULL) return;
	
	gbp = (GBBlockPtr) vnp->data.ptrvalue;
	/* No GBBlock, no keywords to remove */
	if (gbp == NULL) return;
  	
  prev_keyword = NULL;
	for (vnp = gbp->keywords; vnp; vnp = next_keyword) {
	  next_keyword = vnp->next;
		if (StringICmp((CharPtr)vnp->data.ptrvalue, scfp->keyword) == 0) {
			if (prev_keyword == NULL)
			{
			  gbp->keywords = next_keyword;
			}
			else
			{
			  prev_keyword->next = next_keyword;
			}
			vnp->next = NULL;
			ValNodeFreeData (vnp);
		}
		else
		{
		  prev_keyword = vnp;
		}
	}
}

static void ApplyKeyword (IteM i, CharPtr keyword)
{
  BaseFormPtr     bfp;
  SeqEntryPtr     sep;
  KeywordFormData kfd;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  kfd.pfp = NULL;
  kfd.fsp = NULL;
  kfd.keyword = keyword;
  
  VisitBioseqsInSep (sep, &kfd, ApplyKeywordCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();   
}

extern void ApplyGDSKeyword (IteM i)
{
  ApplyKeyword (i, "GDS");
}

extern void ApplyTPAInferentialKeyword (IteM i)
{
  ApplyKeyword (i, "TPA:inferential");
}

extern void ApplyTPAExperimentalKeyword (IteM i)
{
  ApplyKeyword (i, "TPA:experimental");
}

extern void ApplyTPAReassemblyKeyword (IteM i)
{
  ApplyKeyword (i, "TPA:reassembly");
}

static void DoApplyKeywords (ButtoN b)
{
  KeywordFormPtr scfp;
  SeqEntryPtr    sep;
  
  scfp = (KeywordFormPtr) GetObjectExtra (b);
  if (scfp == NULL)
  {
    return;
  }
  
  scfp->pfp = (ParseFieldPtr) DialogToPointer (scfp->string_src_dlg);
  scfp->fsp = FilterSetNew ();
  scfp->fsp->scp = (StringConstraintXPtr) DialogToPointer (scfp->string_constraint_dlg);
  scfp->keyword = SaveStringFromText (scfp->keyword_txt);
  
  sep = GetTopSeqEntryForEntityID (scfp->input_entityID);
  if (sep == NULL) return;
  
  VisitBioseqsInSep (sep, scfp, ApplyKeywordCallback);

  scfp->fsp = FilterSetFree (scfp->fsp);
  scfp->pfp = ParseFieldFree (scfp->pfp);

  ObjMgrSetDirtyFlag (scfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, scfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();   
  Remove (scfp->form);  
}

extern void ApplyKeywordWithStringConstraint (IteM i)
{
  BaseFormPtr    bfp;
  KeywordFormPtr scfp;
  WindoW         w;
  PrompT         ppt;
  GrouP          h, g, c;
  ButtoN         b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  scfp = (KeywordFormPtr) MemNew (sizeof (KeywordFormData));
  if (scfp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Apply Keywords", StdCloseWindowProc);
  SetObjectExtra (w, scfp, StdCleanupExtraProc);
  scfp->form = (ForM) w;
  scfp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Apply Keyword", 0, 0, programFont, 'l');
  scfp->keyword_txt = DialogText (g, "", 30, NULL);
  
  ppt = StaticPrompt (h, "Where", 0, 0, programFont, 'l');
  scfp->string_src_dlg = ParseFieldDestDialogEx (h, NULL, NULL, FALSE, TRUE);

  scfp->string_constraint_dlg = StringConstraintDialogX (h, NULL, FALSE);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoApplyKeywords);
  SetObjectExtra (b, scfp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g,
                              (HANDLE) ppt,
                              (HANDLE) scfp->string_src_dlg,
                              (HANDLE) scfp->string_constraint_dlg,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}

static void DoRemoveKeywords (ButtoN b)
{
  KeywordFormPtr scfp;
  SeqEntryPtr    sep;
  
  scfp = (KeywordFormPtr) GetObjectExtra (b);
  if (scfp == NULL)
  {
    return;
  }
  
  scfp->pfp = (ParseFieldPtr) DialogToPointer (scfp->string_src_dlg);
  scfp->fsp = FilterSetNew ();
  scfp->fsp->scp = (StringConstraintXPtr) DialogToPointer (scfp->string_constraint_dlg);
  scfp->keyword = SaveStringFromText (scfp->keyword_txt);
  
  sep = GetTopSeqEntryForEntityID (scfp->input_entityID);
  if (sep == NULL) return;
  
  VisitBioseqsInSep (sep, scfp, RemoveKeywordCallback);

  scfp->fsp = FilterSetFree (scfp->fsp);
  scfp->pfp = ParseFieldFree (scfp->pfp);

  ObjMgrSetDirtyFlag (scfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, scfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();   
  Remove (scfp->form);  
}

extern void RemoveKeywordWithStringConstraint (IteM i)
{
  BaseFormPtr    bfp;
  KeywordFormPtr scfp;
  WindoW         w;
  PrompT         ppt;
  GrouP          h, g, c;
  ButtoN         b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  scfp = (KeywordFormPtr) MemNew (sizeof (KeywordFormData));
  if (scfp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove Keywords", StdCloseWindowProc);
  SetObjectExtra (w, scfp, StdCleanupExtraProc);
  scfp->form = (ForM) w;
  scfp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Remove Keyword", 0, 0, programFont, 'l');
  scfp->keyword_txt = DialogText (g, "", 30, NULL);
  
  ppt = StaticPrompt (h, "Where", 0, 0, programFont, 'l');
  scfp->string_src_dlg = ParseFieldDestDialogEx (h, NULL, NULL, FALSE, TRUE);

  scfp->string_constraint_dlg = StringConstraintDialogX (h, NULL, FALSE);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoRemoveKeywords);
  SetObjectExtra (b, scfp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g,
                              (HANDLE) ppt,
                              (HANDLE) scfp->string_src_dlg,
                              (HANDLE) scfp->string_constraint_dlg,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


static CharPtr RNAstrandcmd = NULL;

typedef struct rnastrand 
{
  FORM_MESSAGE_BLOCK
  
  ParData rnaParFmt;
  ColData rnaColFmt[3];  
  
  DoC        doc;
  ButtoN     rev_feats;
  ButtoN     use_smart_btn;
  SeqEntryPtr sep;
  ValNodePtr sequence_list;
  Int4       num_sequences;
  BoolPtr    selected;
  Int2       lineheight;  
  CharPtr    database;
  Boolean    use_smart;
} RNAStrandData, PNTR RNAStrandPtr;

typedef enum 
{
  RNAstrand_PLUS = 1,
  RNAstrand_MINUS,
  RNAstrand_MIXED,
  RNAstrand_NO_HITS,
  RNAstrand_UNEXPECTED,
  RNAstrand_PARSE_ERROR,
  RNAstrand_IN_PROGRESS
} ERNAstrand_return_val;

static CharPtr RNAstrand_strings[] = 
{ "Plus", "Minus", "Mixed", "No Hits", "Unexpected", "Parse Error", "In Progress" };


typedef enum 
{
  RNASTRANDGRP_ERROR = 0,
  RNASTRANDGRP_NO_HITS,
  RNASTRANDGRP_MIXED,
  RNASTRANDGRP_MINUS,
  RNASTRANDGRP_PLUS
} ERNAstrand_group_num;

static ERNAstrand_group_num GetRNAStrandGroupNum (ValNodePtr vnp)
{
  if (vnp == NULL) return RNASTRANDGRP_ERROR;
  else if (vnp->choice == RNAstrand_NO_HITS) return RNASTRANDGRP_NO_HITS;
  else if (vnp->choice == RNAstrand_MINUS) return RNASTRANDGRP_MINUS;
  else if (vnp->choice == RNAstrand_MIXED) return RNASTRANDGRP_MIXED;
  else if (vnp->choice == RNAstrand_PLUS) return RNASTRANDGRP_PLUS;
  else return RNASTRANDGRP_ERROR;
}


static void LimitAlignmentResults (SeqAlignPtr salp, Int4 num_results)
{
  Int4        k = 0;
  SeqAlignPtr tmp_salp, last_salp = NULL;
  
  while (salp != NULL && k < num_results)
  {
    last_salp = salp;
    salp = salp->next;
    k++;
  }
  if (last_salp != NULL)
  {
    last_salp->next = NULL;
  }
  while (salp != NULL)
  {
    tmp_salp = salp->next;
    salp->next = NULL;
    salp = SeqAlignFree (salp);
    salp = tmp_salp;
  }
}


static SBlastOptions*
RNABlastOptionNew(void)

{
	SBlastOptions* options;
	Int2           rval;
	Blast_SummaryReturn *extra_returns;


  extra_returns = Blast_SummaryReturnNew();
  rval = SBlastOptionsNew("blastn", &options,
                          extra_returns);

	if (options == NULL)
		return NULL;

  /* This replaces:
   * options->expect_value = 1; 
   */
  SBlastOptionsSetEvalue(options, 1);

  /* This replaces:
   * options->filter_string = StringSave("m L"); 
   */
  SBlastOptionsSetFilterString(options, "m L");
  
  /* This replaces the following:
   * options->mb_template_length = 18;
   * options->mb_disc_type = MB_WORD_CODING; 
   * options->is_megablast_search = TRUE;
   * options->discontinuous = TRUE;
   */
  SBlastOptionsSetDiscMbParams(options, 18, MB_WORD_CODING);

  /* This replaces:
   * options->wordsize = 11; \
   */
	SBlastOptionsSetWordSize (options, 11);
	
  /* This replaces:
   * options->hitlist_size = 20; 
   */
  options->hit_options->hitlist_size = 20;
  
  /* This replaces the following:
   * options->multiple_hits_only  = TRUE;
   * options->window_size = 40; 
   */
  options->word_options->window_size = 40;
  
  /* This replaces the following:
   * options->reward = 1;
	 * options->penalty = -3;
	 * options->gap_open = 5;
	 * options->gap_extend = 2; 
	 */
  SBlastOptionsSetRewardPenaltyAndGapCosts(options, 2, -3, 5, 2, FALSE);
	
	extra_returns = Blast_SummaryReturnFree(extra_returns);
	return options;
}

#if 0
const CharPtr kRNAStrandDatabaseName = "rRNAstrand";
#else
const CharPtr kRNAStrandDatabaseName = "rRNA_blast";
#endif

static Int4
RNAScreenSequence(BioseqPtr bsp, CharPtr database, SeqAlignPtr PNTR seqalign_ptr)

{
	SBlastOptions *blast_options;
	Int2 retval=0;
	SeqAlignPtr seqalign = NULL;
	SeqLocPtr   slp;
	SBlastSeqalignArray* seqalign_arr = NULL;
	Blast_SummaryReturn *extra_returns;

	if (bsp == NULL)
		return -1;

	if (seqalign_ptr)
		*seqalign_ptr = NULL;

	blast_options = RNABlastOptionNew();
	if (blast_options == NULL)
		return -1;
	
	slp = SeqLocWholeNew(bsp);
	if (database == NULL) 
	  database = kRNAStrandDatabaseName;
	
	extra_returns = Blast_SummaryReturnNew();

  retval = Blast_DatabaseSearch(slp, NULL, database, NULL, blast_options, NULL, &seqalign_arr, NULL, extra_returns);
	extra_returns = Blast_SummaryReturnFree(extra_returns);
	blast_options = SBlastOptionsFree(blast_options);  
	
	if (retval == 0 && seqalign_arr != NULL && seqalign_arr->num_queries >0)
	{
	  seqalign = seqalign_arr->array[0];
	  seqalign_arr->array[0] = NULL;
	}
	
	seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
	
	/* limit results to first 20 alignments, as SMART does */
	LimitAlignmentResults (seqalign, 20);	
	
	if (seqalign)
	{
		if (seqalign_ptr)
			*seqalign_ptr = seqalign;
	}	


	return retval;
}


typedef struct rnastrandcollection
{
  ValNodePtr sequence_list;
  Char       path [PATH_MAX];
  CharPtr    database;
  Int4       count;
  Int4       total;
  MonitorPtr mon;
} RNAStrandCollectionData, PNTR RNAStrandCollectionPtr;

static Uint1 GetStatusForAlignmentList (SeqAlignPtr salp)
{
  Uint1 status = RNAstrand_NO_HITS;
  Uint1 strand;
  
  while (salp != NULL)
  {
    AlnMgr2IndexSingleChildSeqAlign(salp);
    strand = AlnMgr2GetNthStrand(salp, 1);
    if (status == RNAstrand_NO_HITS)
    {
      if (strand == Seq_strand_plus)
      {
        status = RNAstrand_PLUS;
      }
      else if (strand == Seq_strand_minus)
      {
        status = RNAstrand_MINUS;
      }
      else
      {
        return RNAstrand_UNEXPECTED;
      }
    }
    else if (strand == Seq_strand_plus)
    {
      if (status != RNAstrand_PLUS)
      {
        return RNAstrand_MIXED;
      }
    }
    else if (strand == Seq_strand_minus)
    {
      if (status != RNAstrand_MINUS)
      {
        return RNAstrand_MIXED;
      }
    }
    else
    {
      return RNAstrand_UNEXPECTED;
    }
    salp = salp->next;
  }
  return status;  
}

static void GetOneRNAStrandednessInfo (BioseqPtr bsp, Pointer userdata)
{
  RNAStrandCollectionPtr rscp;
  SeqAlignPtr            salp = NULL;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  rscp = (RNAStrandCollectionPtr) userdata;

  if (rscp->mon != NULL) {
    MonitorIntValue (rscp->mon, rscp->count);
  }
  (rscp->count)++;
  
  RNAScreenSequence(bsp, rscp->path, &salp);

  ValNodeAddPointer (&(rscp->sequence_list), 
                     GetStatusForAlignmentList(salp),
                     SeqIdFindBest (bsp->id, SEQID_GENBANK));
  salp = SeqAlignFree (salp);  
}

static void CountRNASequences (BioseqPtr bsp, Pointer userdata)
{
  RNAStrandCollectionPtr rscp;

  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  rscp = (RNAStrandCollectionPtr) userdata;
  
  rscp->total++;  
}

static ValNodePtr GetRNAStrandednessFromLocalDatabase (SeqEntryPtr sep, CharPtr database)
{
  RNAStrandCollectionData rscd;
  
  if (sep == NULL) return NULL;

  rscd.path [0] = '\0';
  GetAppParam ("NCBI", "NCBI", "DATA", "", rscd.path, sizeof (rscd.path));
  FileBuildPath (rscd.path, NULL, database);
  
  rscd.database = database;
  rscd.sequence_list = NULL;
  rscd.total = 0;
  rscd.count = 0;
  rscd.mon = NULL;
  
  VisitBioseqsInSep (sep, &rscd, CountRNASequences);
  if (rscd.total > 2)
  {
    rscd.mon = MonitorIntNewEx ("RNA Strand Progress", 0, rscd.total - 1, FALSE);
  }
  
  VisitBioseqsInSep (sep, &rscd, GetOneRNAStrandednessInfo);  
  
  if (rscd.mon != NULL) {
    rscd.mon = MonitorFree (rscd.mon);
    Update ();
  }
  
  return rscd.sequence_list;
}

static void GetAccessionList (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR sequence_list;
  SeqIdPtr        sip;

  if (bsp == NULL || userdata == NULL) return;

  sequence_list = (ValNodePtr PNTR) userdata;

  for (sip = bsp->id; sip != NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENBANK)
    {
      ValNodeAddPointer (sequence_list, 0, sip);
      return;
    }
  }
}


/* Looks at a portion of the list of sequences for strand correction.
 * Returns the number of sequences examined.
 */
static Int4
GetSubListForRNAStrandCorrection
(ValNodePtr start_list,
 Int4       num_seqs,
 CharPtr    RNAstrandcmd)
{
  Int4                    seq_num, cmd_len = 0, k;
  ValNodePtr              vnp;
  FILE *                  fp;
  CharPtr                 args = NULL, cp, cp2, cmmd;
  Char                    tmp_id [256];
  Char                    path [PATH_MAX];
  Char                    file_line [256];
  Boolean                 found_id;
#ifdef OS_MAC
  CharPtr                 cmd_format = "%s -a \'%s\' > %s"; /* just to allow compilation */
#endif
#ifdef OS_UNIX
  CharPtr                 cmd_format = "%s -a \'%s\' > %s";
#endif
#ifdef OS_MSWIN
  CharPtr                 cmd_format = "%s -a \"%s\" > %s";
#endif

  if (start_list == NULL || num_seqs < 1 || StringHasNoText (RNAstrandcmd))
  {
    return 0;
  }

  TmpNam (path);
 
  /* calculate length of string needed for command */
  for (vnp = start_list, seq_num = 0;
           vnp != NULL && seq_num < num_seqs;
           vnp = vnp->next, seq_num++)
  {
    SeqIdWrite (vnp->data.ptrvalue, tmp_id, PRINTID_TEXTID_ACC_ONLY, sizeof (tmp_id) - 1);
    cmd_len += StringLen (tmp_id) + 3;
  }

  args = (CharPtr) MemNew (cmd_len * sizeof (Char));
  if (args == NULL)
  {
        Message (MSG_ERROR, "Unable to allocate memory for strand script argument list");
    return 0;
  }

  cp = args;
  for (vnp = start_list, seq_num = 0;
           vnp != NULL && seq_num < num_seqs;
           vnp = vnp->next, seq_num++)
  {
    SeqIdWrite (vnp->data.ptrvalue, cp, PRINTID_TEXTID_ACC_ONLY, cmd_len - (cp - args) - 1);
    cp += StringLen (cp);
    if (vnp->next != NULL && seq_num < num_seqs - 1)
    {
#ifdef OS_UNIX
      StringCat (cp, ",");
      cp ++;
#else
      StringCat (cp, ", ");
      cp += 2;
#endif
    }
  }


  cmd_len += 3 + StringLen (cmd_format) + StringLen (RNAstrandcmd) + StringLen (path);
  cmmd = (CharPtr) MemNew (cmd_len * sizeof (Char));
  if (cmmd == NULL)
  {
    args = MemFree (args);
        Message (MSG_ERROR, "Unable to allocate memory for RNA strand script command");
    return 0;
  }

#ifdef OS_UNIX
  sprintf (cmmd, cmd_format, RNAstrandcmd, args, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, cmd_format, RNAstrandcmd, args, path);
  RunSilent (cmmd);
#endif
 
  args = MemFree (args);
  cmmd = MemFree (cmmd);
 
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return 0;
  }


  while (fgets (file_line, sizeof (file_line) - 1, fp) != NULL)
  {
    /* find SeqId that matches file line */
    cp = StringChr (file_line, '\t');
    if (cp == NULL)
    {
      continue;
    }
    *cp = 0;
    cp++;
    cp2 = StringChr (cp, '\n');
    if (cp2 != NULL)
    {
      *cp2 = 0;
    }
    found_id = FALSE;
    for (vnp = start_list, seq_num = 0;
             vnp != NULL && seq_num < num_seqs;
             vnp = vnp->next, seq_num++)
        {
      SeqIdWrite (vnp->data.ptrvalue, tmp_id, PRINTID_TEXTID_ACC_ONLY, sizeof (tmp_id) - 1);
      if (StringCmp (tmp_id, file_line) == 0)
      {
        for (k = 1; k <= RNAstrand_IN_PROGRESS; k++)
        {
          if (StringCmp (cp, RNAstrand_strings [k - 1]) == 0
              || (k == RNAstrand_MIXED
                  && StringNCmp (cp, RNAstrand_strings [k - 1],
                                 StringLen (RNAstrand_strings[k - 1])) == 0))
          {
            vnp->choice = k;
          }
        }
      }
    }
  }

  FileClose (fp);

  FileRemove (path);
  return seq_num;
}


static ValNodePtr GetStrandednessFromSMART (SeqEntryPtr sep)
{
  Char                    file_line [256];
  ValNodePtr              sequence_list = NULL, vnp;
  Int4                    num_sequences = 0, num_inspected;

  if (sep == NULL) return NULL;

  if (RNAstrandcmd == NULL) {
    if (GetAppParam ("SEQUIN", "RNACORRECT", "RNASTRAND", NULL, file_line, sizeof (file_line))) {
        RNAstrandcmd = StringSaveNoNull (file_line);
    }
  }
  if (RNAstrandcmd == NULL)
  {
    Message (MSG_ERROR, "RNASTRAND not set in config file!");
    return NULL;
  }


  VisitBioseqsInSep (sep, &sequence_list, GetAccessionList);

  if (sequence_list == NULL)
  {
    Message (MSG_ERROR, "No sequences with accession numbers found!\n");
    return NULL;
  }

  vnp = sequence_list;
  while ((num_inspected = GetSubListForRNAStrandCorrection (vnp, 25, RNAstrandcmd)) != 0)
  {
    num_sequences += num_inspected;
        while (num_inspected > 0 && vnp != NULL)
        {
          vnp = vnp->next;
          num_inspected --;
        }
  }

  return sequence_list;
}


static void DoOneReverse (ValNodePtr vnp, Boolean rev_feats, FILE *fp)
{
  BioseqPtr   bsp;
  Char        str [128];
  SeqEntryPtr sep;

  if (vnp == NULL)
  {
    return;
  }
  
  /* reverse sequence */
  bsp = BioseqFind (vnp->data.ptrvalue);
  BioseqRevComp (bsp);
  sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  VisitAlignmentsInSep (sep, (Pointer) bsp, ReverseBioseqInAlignment);
  
  if (fp != NULL && bsp != NULL && bsp->id != NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), str, PRINTID_REPORT, sizeof (str));
    fprintf (fp, "%s\n", str);
  }
  
  if (rev_feats)
  {
    ReverseBioseqFeatureStrands (bsp);
  }
}

static void RemoveAlignmentsWithSequenceCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL || salp->idx.deleteme) return;
  sip = (SeqIdPtr) userdata;
  while (sip != NULL && !sap->idx.deleteme) {
    if (FindSeqIdinSeqAlign (salp, sip)) {
  	  sap->idx.deleteme = TRUE;
  	}
  	sip = sip->next;
  }
}

extern void RemoveAlignmentsWithSequence (BioseqPtr bsp, Uint2 input_entityID)
{
  SeqEntryPtr           topsep;

  if (bsp == NULL) return;
  topsep = GetTopSeqEntryForEntityID (input_entityID);

  VisitAnnotsInSep (topsep, bsp->id, RemoveAlignmentsWithSequenceCallback);
}

extern void FlipEntireAlignmentIfAllSequencesFlipped (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  ValNodePtr  vnp;
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  Boolean     found;
  Int4 row, num_rows;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL || salp->idx.deleteme) return;
  
  
  AlnMgr2IndexSingleChildSeqAlign(salp);
  num_rows = AlnMgr2GetNumRows(salp);
  for (row = 1; row <= num_rows; row++) {
    sip = AlnMgr2GetNthSeqIdPtr(salp, row);
    found = FALSE;
    vnp = (ValNodePtr)userdata;
    while (vnp != NULL && !found) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (SeqIdOrderInBioseqIdList (sip, bsp->id) > 0) {
        found = TRUE;
      }
      vnp = vnp->next;
    }
    if (!found) return;
  }
  
  FlipAlignment(salp);      
}


static Boolean CheckForAlignmentsBeforeCorrection (RNAStrandPtr strand_info)
{
  ValNodePtr vnp, seq_in_aln = NULL, aln_bsp = NULL;
  BioseqPtr  bsp;
  Char       id_str[255];
  CharPtr    msg;
  MsgAnswer  ans;
  Boolean    retval = TRUE;
  Int4       strand_num = 0;
  Uint2      entityID = 0;


  if (strand_info == NULL || strand_info->sequence_list == NULL) return FALSE;

  for (vnp = strand_info->sequence_list, strand_num = 0; vnp != NULL; vnp = vnp->next, strand_num++) {
    if (!strand_info->selected [strand_num]) continue;
    bsp = BioseqFind (vnp->data.ptrvalue);
    if (bsp != NULL && IsBioseqInAnyAlignment (bsp, bsp->idx.entityID)) {
      entityID = bsp->idx.entityID;
      SeqIdWrite (vnp->data.ptrvalue, id_str, PRINTID_REPORT, sizeof (id_str) - 1);
      ValNodeAddPointer (&seq_in_aln, 0, StringSave (id_str));
      ValNodeAddPointer (&aln_bsp, 0, bsp);
    }
  }
  if (seq_in_aln != NULL) {
  	msg = CreateListMessage ("Sequence", 
  	                         seq_in_aln->next == NULL 
  	                         ? " is in an alignment, which will become invalid.  Do you want to remove the alignment?" 
  	                         : " are in alignments, which will become invalid unless all of the sequences in these alignments are reversed.  Do you want to remove the alignments?",
  	                         seq_in_aln);
  	seq_in_aln = ValNodeFreeData (seq_in_aln);
    ans = Message (MSG_YNC, msg);
    msg = MemFree (msg);
    if (ans == ANS_YES) {
      /* remove the alignments */
      for (vnp = aln_bsp; vnp != NULL; vnp = vnp->next) {
          bsp = vnp->data.ptrvalue;
          RemoveAlignmentsWithSequence (bsp, bsp->idx.entityID);
          DeleteMarkedObjects (bsp->idx.entityID, 0, NULL);      
      }
    } else if (ans == ANS_CANCEL) {
      retval = FALSE;
    } else {
      VisitAnnotsInSep (GetTopSeqEntryForEntityID (entityID), aln_bsp, FlipEntireAlignmentIfAllSequencesFlipped);
    }
    aln_bsp = ValNodeFree (aln_bsp);
  }
  return retval;
}

static void DoRNACorrection (ButtoN b)
{
  RNAStrandPtr strand_info;
  Int4         strand_num = 0;
  ValNodePtr   vnp;
  Boolean      rev_feats;
  FILE         *fp;
  Char         path [PATH_MAX];
  
  strand_info = (RNAStrandPtr) GetObjectExtra (b);
  if (strand_info == NULL)
  {
    return;
  }
  
  rev_feats = GetStatus (strand_info->rev_feats);

  /* check for alignments */
  if (!CheckForAlignmentsBeforeCorrection (strand_info)) {
    return;
  }

  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp != NULL) {
    for (vnp = strand_info->sequence_list, strand_num = 0;
         vnp != NULL; 
         vnp = vnp->next, strand_num++)
    {
      if (strand_info->selected [strand_num])
      {
        /* reverse sequence */
        DoOneReverse (vnp, rev_feats, fp);
      }
    }
  
    FileClose (fp);
    LaunchGeneralTextViewer (path, "Flipped Sequences");
    FileRemove (path);  
  } else {
    Message (MSG_ERROR, "Unable to open file for flip report");
  }
    
  ObjMgrSetDirtyFlag (strand_info->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, strand_info->input_entityID, 0, 0);
  Remove (strand_info->form);
  Update ();
}

static void CleanupRNAStrandFormProc (GraphiC g, Pointer data)
{
  RNAStrandPtr strand_info;
  
  if (data != NULL)  
  {
    strand_info = (RNAStrandPtr) data;
    strand_info->sequence_list = ValNodeFree (strand_info->sequence_list);
    strand_info->selected = MemFree (strand_info->selected);
    strand_info = MemFree (strand_info);
  }
}

static Int4 GetSeqNumFromListPos (Int4 list_pos, ValNodePtr sequence_list)
{
  Int4       seq_num, list_offset;
  ValNodePtr vnp;
  ERNAstrand_group_num group;
  
  if (sequence_list == NULL)
  {
    return 0;
  }

  list_offset = -1;  
  for (group = RNASTRANDGRP_ERROR; group <= RNASTRANDGRP_PLUS && list_offset != list_pos; group++) 
  {
    for (vnp = sequence_list, seq_num = -1; 
         vnp != NULL && list_offset != list_pos; 
         vnp = vnp->next, seq_num++) 
    {
      if (GetRNAStrandGroupNum(vnp) == group) 
      {
        list_offset ++;
      }
    }
  }
  
  return seq_num;
}

static void ReleaseRNAStrand (DoC d, PoinT pt)

{
  Int2            col;
  RNAStrandPtr    strand_info;
  Int2            item;
  RecT            rct;
  Int2            row;
  Int4            seq_num;

  strand_info = (RNAStrandPtr) GetObjectExtra (d);
  if (strand_info != NULL && strand_info->selected != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, &rct);
    rct.left += 1;
    rct.right = rct.left + strand_info->lineheight;
    rct.bottom = rct.top + (rct.right - rct.left);
    if (row == 1 && col == 3 && PtInRect (pt, &rct))
    {
      seq_num = GetSeqNumFromListPos (item - 1, strand_info->sequence_list);
      if (seq_num > -1 && seq_num < strand_info->num_sequences)
      {
        if (strand_info->selected [seq_num]) {
          strand_info->selected [seq_num] = FALSE;
        } else {
          strand_info->selected [seq_num] = TRUE;
        }
        InsetRect (&rct, -1, -1);
        InvalRect (&rct);
        Update ();
      }
    }
  }
}

static void DrawRNAStrand (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RNAStrandPtr strand_info;
  RecT         rct;
  RecT         doc_rect;
  Int4         seq_num;

  strand_info = (RNAStrandPtr) GetObjectExtra (d);
  
  if (strand_info == NULL || r == NULL 
      || item < 1 || item > strand_info->num_sequences 
      || firstLine != 0)
  {
    return;
  }
  
  rct = *r;
  rct.right --;
  rct.left = rct.right - strand_info->lineheight;
  rct.bottom = rct.top + (rct.right - rct.left);
  
  /* make sure we don't draw a box where we aren't drawing text */
  ObjectRect (strand_info->doc, &doc_rect);
  InsetRect (&doc_rect, 4, 4);
  if (rct.bottom > doc_rect.bottom)
  {
    return;
  }
  
  FrameRect (&rct);
  
  seq_num = GetSeqNumFromListPos (item - 1, strand_info->sequence_list);
  if (seq_num > -1 && seq_num < strand_info->num_sequences) {
    if (strand_info->selected != NULL && strand_info->selected [seq_num]) {
      MoveTo (rct.left, rct.top);
      LineTo (rct.right - 1, rct.bottom - 1);
      MoveTo (rct.left, rct.bottom - 1);
      LineTo (rct.right - 1, rct.top);
    }
  }
}

static void GetRNAStrandDesc (ValNodePtr vnp, CharPtr doc_line)
{
  if (vnp == NULL || doc_line == NULL) return;
  
  if (vnp->choice == RNAstrand_MINUS) {
    StringCat (doc_line, "\tMINUS\t\n");
  } else if (vnp->choice == RNAstrand_PLUS) {
    StringCat (doc_line, "\tPLUS\t\n");
  } else {
    StringCat (doc_line, "\t");
    if (vnp->choice == 0)
    {
      StringCat (doc_line, "Unknown error");
    }
    else
    {
      StringCat (doc_line, RNAstrand_strings [vnp->choice - 1]);
    }
    StringCat (doc_line, "\t\n");
  }   
}

static void RedrawRNAStrandDialog (RNAStrandPtr strand_info)
{
  Int4         strand_num;
  Char         doc_line [500];
  ValNodePtr   vnp;
  ERNAstrand_group_num group_num;
  
  if (strand_info == NULL)
  {
    return;
  }

  Reset (strand_info->doc);
  
  for (group_num = RNASTRANDGRP_ERROR;
       group_num <= RNASTRANDGRP_PLUS;
       group_num++) 
  {
    for (vnp = strand_info->sequence_list, strand_num = 0;
         vnp != NULL;
         vnp = vnp->next, strand_num++)
    {    
      if (GetRNAStrandGroupNum(vnp) != group_num) continue;
      SeqIdWrite (vnp->data.ptrvalue, doc_line, PRINTID_TEXTID_ACC_ONLY, 255);
      GetRNAStrandDesc (vnp, doc_line);
      AppendText (strand_info->doc, doc_line, &(strand_info->rnaParFmt), strand_info->rnaColFmt, programFont);
      if (group_num == RNASTRANDGRP_MINUS) strand_info->selected[strand_num] = TRUE;
    }
  }  

  UpdateDocument (strand_info->doc, 0, 0);  
}

static void RefreshRNAStrandDialog (ButtoN b)
{
  RNAStrandPtr strand_info;
  ValNodePtr   new_seq_list;
  
  strand_info = (RNAStrandPtr) GetObjectExtra (b);
  if (strand_info == NULL)
  {
    return;
  }

  if (strand_info->use_smart) {
    new_seq_list = GetStrandednessFromSMART (strand_info->sep);
  } else { 
    new_seq_list = GetRNAStrandednessFromLocalDatabase (strand_info->sep, strand_info->database); 
  }
  
  if (new_seq_list == NULL)
  {
    Message (MSG_ERROR, "Unable to refresh sequence list.");
    return;
  }
  strand_info->sequence_list = ValNodeFree (strand_info->sequence_list);
  strand_info->sequence_list = new_seq_list;
    
  RedrawRNAStrandDialog (strand_info);
}


static void ChangeStrandSmart (ButtoN b)
{
  RNAStrandPtr strand_info;
  
  strand_info = (RNAStrandPtr) GetObjectExtra (b);
  if (strand_info == NULL)
  {
    return;
  }
  strand_info->use_smart = GetStatus (strand_info->use_smart_btn);
  RefreshRNAStrandDialog (strand_info->use_smart_btn); 
}


static Int2 CorrectRNAStrandednessForEntityID (Uint2 entityID, Boolean use_smart)

{
  SeqEntryPtr             sep;
  ValNodePtr              sequence_list = NULL;
  RNAStrandPtr            strand_info;
  WindoW                  w;
  GrouP                   h, c;
  ButtoN                  b;
  RecT                    r;
  ButtoN                  refresh_btn;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;

  if (use_smart)
  {
    sequence_list = GetStrandednessFromSMART (sep);
  }
  else 
  {  
    sequence_list = GetRNAStrandednessFromLocalDatabase (sep, kRNAStrandDatabaseName); 
  }
  
  if (sequence_list == NULL)
  {
    return OM_MSG_RET_ERROR;
  }
  
  strand_info = (RNAStrandPtr) MemNew (sizeof (RNAStrandData));
  if (strand_info == NULL)
  {
    return OM_MSG_RET_ERROR;
  }
  
  strand_info->input_entityID = entityID;
  strand_info->sep = sep;
  strand_info->sequence_list = sequence_list;
  strand_info->num_sequences = ValNodeLen (sequence_list);
  strand_info->selected = (BoolPtr) MemNew (strand_info->num_sequences * sizeof (Boolean));
  strand_info->database = kRNAStrandDatabaseName;
  strand_info->use_smart = use_smart;

  /* initialize document paragraph format */
  strand_info->rnaParFmt.openSpace = FALSE;
  strand_info->rnaParFmt.keepWithNext = FALSE;
  strand_info->rnaParFmt.keepTogether = FALSE;
  strand_info->rnaParFmt.newPage = FALSE;
  strand_info->rnaParFmt.tabStops = FALSE;
  strand_info->rnaParFmt.minLines = 0;
  strand_info->rnaParFmt.minHeight = 0;
  
  /* initialize document column format */
  strand_info->rnaColFmt[0].pixWidth = 0;
  strand_info->rnaColFmt[0].pixInset = 0;
  strand_info->rnaColFmt[0].charWidth = 80;
  strand_info->rnaColFmt[0].charInset = 0;
  strand_info->rnaColFmt[0].font = NULL;
  strand_info->rnaColFmt[0].just = 'l';
  strand_info->rnaColFmt[0].wrap = TRUE;
  strand_info->rnaColFmt[0].bar = FALSE;
  strand_info->rnaColFmt[0].underline = FALSE;
  strand_info->rnaColFmt[0].left = FALSE;
  strand_info->rnaColFmt[0].last = FALSE;
  strand_info->rnaColFmt[1].pixWidth = 0;
  strand_info->rnaColFmt[1].pixInset = 0;
  strand_info->rnaColFmt[1].charWidth = 80;
  strand_info->rnaColFmt[1].charInset = 0;
  strand_info->rnaColFmt[1].font = NULL;
  strand_info->rnaColFmt[1].just = 'l';
  strand_info->rnaColFmt[1].wrap = TRUE;
  strand_info->rnaColFmt[1].bar = FALSE;
  strand_info->rnaColFmt[1].underline = FALSE;
  strand_info->rnaColFmt[1].left = FALSE;
  strand_info->rnaColFmt[1].last = FALSE;
  strand_info->rnaColFmt[2].pixWidth = 0;
  strand_info->rnaColFmt[2].pixInset = 0;
  strand_info->rnaColFmt[2].charWidth = 80;
  strand_info->rnaColFmt[2].charInset = 0;
  strand_info->rnaColFmt[2].font = NULL;
  strand_info->rnaColFmt[2].just = 'l';
  strand_info->rnaColFmt[2].wrap = TRUE;
  strand_info->rnaColFmt[2].bar = FALSE;
  strand_info->rnaColFmt[2].underline = FALSE;
  strand_info->rnaColFmt[2].left = FALSE;
  strand_info->rnaColFmt[2].last = TRUE;
    
  w = FixedWindow (50, 33, -10, -10, "Correct RNA Strandedness", NULL);
  SetObjectExtra (w, strand_info, CleanupRNAStrandFormProc);
  strand_info->form = (ForM) w;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  strand_info->doc = DocumentPanel (h, stdCharWidth * 27, stdLineHeight * 8);
  SetObjectExtra (strand_info->doc, strand_info, NULL);
  SetDocAutoAdjust (strand_info->doc, TRUE);
  SetDocProcs (strand_info->doc, NULL, NULL, ReleaseRNAStrand, NULL); 
  SetDocShade (strand_info->doc, DrawRNAStrand, NULL, NULL, NULL);

  SelectFont (programFont);
  strand_info->lineheight = LineHeight ();
  
  ObjectRect (strand_info->doc, &r);
  InsetRect (&r, 4, 4);
  strand_info->rnaColFmt[0].pixWidth = (r.right - r.left - strand_info->lineheight) / 2;
  strand_info->rnaColFmt[1].pixWidth = (r.right - r.left - strand_info->lineheight) / 2;
  strand_info->rnaColFmt[2].pixWidth = strand_info->lineheight;

  refresh_btn = PushButton (h, "Refresh Strand Results", RefreshRNAStrandDialog);
  SetObjectExtra (refresh_btn, strand_info, NULL);
  strand_info->rev_feats = CheckBox (h, "Also reverse features", NULL);
  SetStatus (strand_info->rev_feats, FALSE);
  strand_info->use_smart_btn = CheckBox (h, "Use SMART for strand information", ChangeStrandSmart);
  SetObjectExtra (strand_info->use_smart_btn, strand_info, NULL);
  SetStatus (strand_info->use_smart_btn, strand_info->use_smart);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Autocorrect Minus Strands", DoRNACorrection);
  SetObjectExtra (b, strand_info, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) strand_info->doc,
                              (HANDLE) refresh_btn,
                              (HANDLE) strand_info->use_smart_btn,
                              (HANDLE) strand_info->rev_feats,
                              (HANDLE) c,
                              NULL);  
  RedrawRNAStrandDialog (strand_info);
  Show (w);
  return OM_MSG_RET_OK;
}


static Int2 LIBCALLBACK CorrectRNAStrandednessEx (Pointer data, Boolean use_smart)
{
  OMProcControlPtr ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL)
  {
    Message (MSG_ERROR, "You must select something!");
    return OM_MSG_RET_ERROR;
  }

  return CorrectRNAStrandednessForEntityID (ompcp->input_entityID, use_smart);
}


extern Int2 LIBCALLBACK CorrectRNAStrandedness (Pointer data)
{
  return CorrectRNAStrandednessEx (data, FALSE);
}
extern Int2 LIBCALLBACK CorrectRNAStrandednessUseSmart (Pointer data)
{
  return CorrectRNAStrandednessEx (data, TRUE);
}

extern void CorrectRNAStrandednessMenuItem (IteM i)
{
  BaseFormPtr   bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  CorrectRNAStrandednessForEntityID (bfp->input_entityID, TRUE);
}


/* collection_date has a controlled format.  
 * It is YYYY or Mmm-YYYY or DD-Mmm-YYYY where Mmm = Jan, Feb, Mar, Apr, May, 
 *                                                   Jun, Jul, Aug, Sep, Oct, 
 *                                                   Nov, Dec
 * This function will convert other formats  to this format.
 * For instance, September 12, 2004 should be converted to 12-Sep-2004
 * 12/15/2003 should be converted to 15-Dec-2003.  
 * 
 * If the date supplied is ambiguous (01/03/05), can you allow the indexer to choose which field goes in Mmm and which in DD.
 */

typedef struct parsecollectiondate
{
  Int4    num_successful;
  Int4    num_unsuccessful;
  Boolean month_first;
} ParseCollectionDateData, PNTR ParseCollectionDatePtr;

static void ParseCollectionDateCallback (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr           ssp;
  CharPtr                reformatted_date = NULL;
  ParseCollectionDatePtr pcdp;
  
  if (biop == NULL || biop->subtype == NULL || userdata == NULL)
  {
    return;
  }
  
  pcdp = (ParseCollectionDatePtr) userdata;
  
  ssp = biop->subtype;
  while (ssp != NULL)
  {
    if (ssp->subtype == SUBSRC_collection_date)
    {
      reformatted_date = ReformatDateStringEx (ssp->name, pcdp->month_first, NULL);
      if (reformatted_date == NULL)
      {
        pcdp->num_unsuccessful ++;
      }
      else
      {
        ssp->name = MemFree (ssp->name);
        ssp->name = reformatted_date;
        pcdp->num_successful ++;
      }
    }
    ssp = ssp->next;
  }
}

static void CountParseCollectionDateCallback (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr           ssp;
  CharPtr                reformatted_date = NULL;
  ParseCollectionDatePtr pcdp;
  
  if (biop == NULL || biop->subtype == NULL || userdata == NULL)
  {
    return;
  }
  
  pcdp = (ParseCollectionDatePtr) userdata;
  
  ssp = biop->subtype;
  while (ssp != NULL)
  {
    if (ssp->subtype == SUBSRC_collection_date)
    {
      reformatted_date = ReformatDateStringEx (ssp->name, pcdp->month_first, NULL);
      if (reformatted_date == NULL)
      {
        pcdp->num_unsuccessful++;
      }
      else
      {
        MemFree (reformatted_date);
        pcdp->num_successful ++;
      }
    }
    ssp = ssp->next;
  }
}

static void ParseCollectionDate (IteM i, Boolean month_first)
{
  BaseFormPtr             bfp;
  SeqEntryPtr             sep;
  ParseCollectionDateData pcdd;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  WatchCursor ();
  Update ();

  pcdd.num_successful = 0;
  pcdd.num_unsuccessful = 0;
  pcdd.month_first = month_first;
  
  VisitBioSourcesInSep (sep, &pcdd, CountParseCollectionDateCallback);
  if (pcdd.num_unsuccessful > 0 
      && ANS_NO == Message (MSG_YN, 
                            "Sequin will be unable to reformat %d dates - do you want to continue?", pcdd.num_unsuccessful))
  {
    ArrowCursor ();
    Update ();   
    return;
  }

  pcdd.num_successful = 0;
  pcdd.num_unsuccessful = 0;  
  VisitBioSourcesInSep (sep, &pcdd, ParseCollectionDateCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();   
}

extern void ParseCollectionDateMonthFirst (IteM i)
{
  ParseCollectionDate (i, TRUE);
}

extern void ParseCollectionDateDayFirst (IteM i)
{
  ParseCollectionDate (i, FALSE);
}


static CharPtr GetIDStringForSeqEntry (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  BioseqPtr    bsp = NULL;
  SeqIdPtr     sip = NULL;
  Char         id_txt [100];
  
  if (sep == NULL || sep->data.ptrvalue == NULL)
  {
    return NULL;
  }
  
  if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_nuc_prot
        || bssp->_class == BioseqseqSet_class_segset)
    {
      sep = bssp->seq_set;
      while (sep != NULL && (IS_Bioseq_set (sep) || sep->data.ptrvalue == NULL))
      {
        sep = sep->next;
      }
      if (sep != NULL)
      {
        bsp = (BioseqPtr) sep->data.ptrvalue;
      }
    }
  }
  else if (IS_Bioseq(sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  }
  if (bsp == NULL) return NULL;
  sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  if (sip == NULL) return NULL;
  
  SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  return StringSave (id_txt);
}


static int LIBCALLBACK SortSeqEntryByIDStr (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = GetIDStringForSeqEntry (vnp1);
      str2 = GetIDStringForSeqEntry (vnp2);
      if (str1 != NULL && str2 != NULL) {
        rval = StringICmp (str1, str2);
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    }
  }
  return rval;
}

static void ReorderBioseqSetByAccession (BioseqSetPtr bssp)
{
  SeqEntryPtr       sep;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;

  sep = SeqMgrGetSeqEntryForData (bssp);
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);  

  bssp->seq_set = ValNodeSort (bssp->seq_set, SortSeqEntryByIDStr);
  
  
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  ObjMgrSetDirtyFlag (bssp->idx.entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bssp->idx.entityID, 0, 0);  
}


extern Int2 LIBCALLBACK ReorderSetByAccession (Pointer data)
{
  BioseqSetPtr      bssp;
  OMProcControlPtr  ompcp;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL)
    return OM_MSG_RET_ERROR;

  if (ompcp->input_itemtype != OBJ_BIOSEQSET || ompcp->input_data == NULL) {
    Message (MSG_ERROR, "Must select Bioseq Set!");
    return OM_MSG_RET_ERROR;
  }
  
  bssp = (BioseqSetPtr) ompcp->input_data;

  ReorderBioseqSetByAccession (bssp);
  return OM_MSG_RET_DONE;
}


extern BioseqSetPtr FindTopLevelSetForDesktopFunction (BioseqSetPtr bssp)
{
  if (bssp == NULL) {
    return NULL;
  } else if (bssp->_class != BioseqseqSet_class_genbank) {
    return bssp;
  } else if (bssp->seq_set != NULL && IS_Bioseq_set (bssp->seq_set) && bssp->seq_set->next == NULL) {
    return (BioseqSetPtr) bssp->seq_set->data.ptrvalue;
  } else {
    return bssp;
  }
}


extern void ReorderSetByAccessionMenuItem (IteM i)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-levelset!");
  } else {
    ReorderBioseqSetByAccession (FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue));
  }
}


extern void DescriptorPropagateMenuItem (IteM i) 
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-levelset!");
    return;
  }

  SetDescriptorPropagate (FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue));

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct findcontig {
  FORM_MESSAGE_BLOCK
  DoC        doc;
  Int2       clicked;
  Boolean    dblClick;
  ValNodePtr id_list;
  BioseqPtr  con_bsp;
} FindContigData, PNTR FindContigPtr;

static void ViewContigPiece (SeqIdPtr sip)
{
  SeqEntryPtr    sep = NULL;
  Int4           uid;  
  BioseqPtr      bsp;
  SeqEntryPtr    oldscope;
  Uint2          entityID = 0;
  Uint2          parenttype, last_parenttype;
  Pointer        parentptr, last_parentptr;
  OMUserDataPtr  omudp;
  BaseFormPtr    bfp;

  if (sip == NULL) {
    return;
  }
  
  oldscope = SeqEntrySetScope (NULL);
  bsp = BioseqFind (sip);
  SeqEntrySetScope (oldscope);
  if (bsp == NULL) {
    uid = GetGIForSeqId (sip);
    DownloadAndDisplay (uid);
  } else {
    if (bsp->idx.entityID == 0) {
      sep = SeqMgrGetSeqEntryForData (bsp);
      if (sep != NULL) {
        GetSeqEntryParent (sep, &parentptr, &parenttype);
        last_parentptr = parentptr;
        last_parenttype = parenttype;
        while (parentptr != NULL) {
          if (parenttype == OBJ_SEQENTRY) {
            sep = parentptr;
          } else if (parenttype == OBJ_BIOSEQ || parenttype == OBJ_BIOSEQSET) {
            sep = SeqMgrGetSeqEntryForData (parentptr);
          } else {
            sep = NULL;
          }
          if (sep != NULL) {
            last_parentptr = parentptr;
            last_parenttype = parenttype;
            GetSeqEntryParent (sep, &parentptr, &parenttype);
          }
        }
        if (last_parentptr == NULL) {
          entityID = ObjMgrRegister (OBJ_BIOSEQ, bsp);
        } else {
          entityID = ObjMgrRegister (last_parenttype, last_parentptr);
        }
      }
    } else {
      entityID = bsp->idx.entityID;
    }
    omudp = EntityAlreadyHasViewer (entityID);
    if (omudp == NULL) {
      LaunchDisplay (entityID);
      omudp = EntityAlreadyHasViewer (entityID);
    }
    if (omudp != NULL) {
      MakeViewerIndependent (entityID, omudp);
      bfp = (BaseFormPtr) omudp->userdata.ptrvalue;
      if (bfp != NULL && bfp->form != NULL) {
        Select (bfp->form);
        SeqEntrySetScope (GetTopSeqEntryForEntityID(entityID));
        SetBioseqViewTargetByBioseq (bfp, bsp);
      }
    }
  }
}

static void ClickContigPiece (DoC d, PoinT pt)
{
  Int2                      item;
  Int2                      row;
  Int2                      col;
  FindContigPtr  fcp;
  ValNodePtr     vnp;

  fcp = GetObjectExtra (d);
  if (fcp != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && fcp->clicked == item) {
      fcp->dblClick = dblClick;
    } else {
      fcp->dblClick = FALSE;
    }
    fcp->clicked = 0;
    if (item > 0 && row > 0) {
      fcp->clicked = item;
    }
    if (item > 0 && row > 0 && dblClick)
    {
      vnp = fcp->id_list;
      while (item > 1 && vnp != NULL) {
        vnp = vnp->next;
        item--;
      }
      if (vnp != NULL) {
        ViewContigPiece(vnp->data.ptrvalue);
      }
    }
  }
}

static ParData contigParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};

static Nlm_ColData contigColFmt [5] = {{0, 5, 10, 0, NULL, 'l', 1,0,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,1,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,1,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,1,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,1,0,0, TRUE}};


static void FindContigBtn (ButtoN b)
{
  FindContigPtr     fcp;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp, piece_bsp;
  DeltaSeqPtr       dsp;
  SeqLocPtr         slp;
  SeqLitPtr         slip;
  Int4              seq_len, len, feat_start, feat_offset = 0;
  SeqIdPtr          sip;
  Char              str [128];
  Char              num_str [20];
  SeqEntryPtr       oldscope;
  CharPtr           location, row_text, label;
  CharPtr           header = "Feature\t\t\tAccession\tOffset\n";
  Boolean           some_missing = FALSE;
  
  fcp = (FindContigPtr) GetObjectExtra (b);
  if (fcp == NULL) {
    return;
  }
  Reset (fcp->doc);

  fcp->id_list = FreeSeqIdList (fcp->id_list);
  
  AppendText (fcp->doc, header, &contigParFmt, contigColFmt, programFont);
  ValNodeAddPointer (&fcp->id_list, 0, NULL);
  
  for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
    if (sel->itemtype == OBJ_SEQFEAT) {
      oldscope = SeqEntrySetScope (NULL);
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      SeqEntrySetScope (oldscope);
      if (sfp != NULL) {
        oldscope = SeqEntrySetScope (NULL);
        bsp = BioseqFindFromSeqLoc(sfp->location);
        SeqEntrySetScope (oldscope);
        sip = NULL;
        piece_bsp = NULL;
        if (SeqLocStrand (sfp->location) == Seq_strand_minus) {
          feat_start = SeqLocStop (sfp->location);
        } else {
          feat_start = SeqLocStart (sfp->location);
        }
        if (bsp == NULL) {
          some_missing = TRUE;
        } else if (bsp != fcp->con_bsp) {
          piece_bsp = bsp;
          sip = SeqIdFindBest (bsp->id, 0);
          feat_offset = feat_start;
        } else if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta) {
          seq_len = 0;
          if (bsp->seq_ext_type == 4) {
            for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
                 dsp != NULL && seq_len <= feat_start && sip == NULL; 
                 dsp = dsp->next) 
            {  
              len = 0;
              if (dsp->data.ptrvalue == NULL) {
                continue;           
              } else if (dsp->choice == 1) {
                slp = (SeqLocPtr) dsp->data.ptrvalue;
                len = SeqLocLen (slp);
                if (seq_len <= feat_start && seq_len + len > feat_start) {
                  sip = SeqLocId (slp);
                  feat_offset = feat_start - seq_len;
                }
              } else if (dsp->choice == 2) {
                slip = (SeqLitPtr) dsp->data.ptrvalue;
                len = slip->length;
              }
              seq_len += len;
            }
          } else if (bsp->seq_ext_type == 1) {
            slp = (SeqLocPtr) bsp->seq_ext;
            while (slp != NULL && seq_len <= feat_start && sip == NULL) {
              len = SeqLocLen (slp);
              if (seq_len <= feat_start && seq_len + len > feat_start) {
                sip = SeqLocId (slp);
                feat_offset = feat_start - seq_len;
              } else {
                slp = slp->next;
              }
              seq_len += len;
            }
		  } else {
		    some_missing = TRUE;
		  }
        }
        
        if (piece_bsp == NULL && sip != NULL) {
          oldscope = SeqEntrySetScope (NULL);
          piece_bsp = BioseqFind (sip);
          SeqEntrySetScope (oldscope);
        }
        
        if (sip == NULL || piece_bsp == NULL) {
          some_missing = TRUE;
        } else {
          sip = SeqIdFindBest (piece_bsp->id, SEQID_GENBANK);
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
          sprintf (num_str, "%d", (int) feat_offset + 1);
          location = SeqLocPrintUseBestID (sfp->location);
          label = (CharPtr) FeatDefTypeLabel(sfp);
          row_text = (CharPtr) MemNew (sizeof (Char) * 
                                     (StringLen (label) 
                                      + StringLen (fcontext.label) 
                                      + StringLen (location) 
                                      + StringLen (str)
                                      + StringLen (num_str)
                                      + 6));
          sprintf (row_text, "%s\t%s\t%s\t%s\t%s\n", label, fcontext.label, location, str, num_str);
          location = MemFree (location);
          AppendText (fcp->doc, row_text, &contigParFmt, contigColFmt, programFont);
          row_text = MemFree (row_text);
            
          ValNodeAddPointer (&(fcp->id_list), 0, SeqIdDup (sip));
        }
      }
    }
  }
  UpdateDocument (fcp->doc, 0, 0);
  if (some_missing) {
    Message (MSG_ERROR, "Some features are on sequences that are not far records!\n");
  }
}

static void CleanupFindContigFormProc (GraphiC g, Pointer data)
{
  FindContigPtr fcp;
  
  if (data != NULL)  
  {
    fcp = (FindContigPtr) data;
    fcp->id_list = FreeSeqIdList (fcp->id_list);
  }
  StdCleanupExtraProc (g, data);  
}


static void FindConRecordCallback (BioseqPtr bsp, Pointer userdata)
{
  BioseqPtr *pbsp = (BioseqPtr *) userdata;
  if (bsp == NULL || userdata == NULL || ! ISA_na (bsp->mol)
      || (bsp->repr != Seq_repr_seg && bsp->repr != Seq_repr_delta)
      || *pbsp != NULL) {
      return;
  } else {
      *pbsp = bsp;
  }  
}


static BioseqPtr FindConRecord (Uint2 entityID)
{
  SeqEntryPtr sep;
  BioseqPtr   con_bsp = NULL;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  VisitBioseqsInSep (sep, &con_bsp, FindConRecordCallback);
  return con_bsp;
}


extern void FindContig (IteM i)
{
  BaseFormPtr   bfp;
  FindContigPtr fcp;
  WindoW        w;
  GrouP         h, c;
  ButtoN        b;
  RecT          r;
  Char          str [128];
  CharPtr       title;
  CharPtr       title_fmt = "Find Contig for Features in %s";
  BioseqPtr     con_bsp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  con_bsp = FindConRecord (bfp->input_entityID);
  if (con_bsp == NULL) {
    Message (MSG_ERROR, "No con records found!");
    return;
  }
  
  fcp = MemNew (sizeof (FindContigData));
  if (fcp != NULL) {
    fcp->input_entityID = bfp->input_entityID;
    fcp->con_bsp = con_bsp;
    SeqIdWrite (SeqIdFindBest (fcp->con_bsp->id, SEQID_GENBANK), str, PRINTID_REPORT, sizeof (str));
    title = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + StringLen (title_fmt)));
    sprintf (title, title_fmt, str);
    w = FixedWindow (-50, -33, -10, -10, title, NULL);
    title = MemFree (title);
    SetObjectExtra (w, fcp, CleanupFindContigFormProc);
    fcp->form = (ForM) w;
    fcp->formmessage = DefaultMessageProc;
    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);
    fcp->doc = DocumentPanel (h, 50 * stdCharWidth, 20 * stdLineHeight);
    SetObjectExtra (fcp->doc, fcp, NULL);
    SetDocProcs (fcp->doc, ClickContigPiece, NULL, NULL, NULL);
/*    SetDocNotify (afp->doc, AlignNotifyProc); */

    ObjectRect (fcp->doc, &r);
    InsetRect (&r, 4, 4);
  
    contigColFmt[0].pixWidth = 5 * stdCharWidth;
    contigColFmt[3].pixWidth = 7 * stdCharWidth;
    contigColFmt[4].pixWidth = 5 * stdCharWidth;
    contigColFmt[1].pixWidth = (r.right - r.left - contigColFmt[0].pixWidth - contigColFmt[3].pixWidth - contigColFmt[4].pixWidth) / 2;
    contigColFmt[2].pixWidth = (r.right - r.left - contigColFmt[0].pixWidth - contigColFmt[3].pixWidth - contigColFmt[4].pixWidth) / 2;

    c = HiddenGroup (h, 2, 0, NULL);
    b = PushButton (c, "Find Selected Contig", FindContigBtn);
    SetObjectExtra (b, fcp, NULL);
    FindContigBtn(b);
    b = PushButton (c, "Done", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) fcp->doc, (HANDLE) c, NULL);
    Show (w);
    Select (w);
  }
}


static void ListDeltaCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bsplist;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta || userdata == NULL) {
    return;
  }
  
  bsplist = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (bsplist, 0, bsp);
}


extern Boolean ShowDeltaReport (SeqEntryPtr sep)
{
  ValNodePtr bsplist = NULL, vnp;
  FILE      *fp;
  Char       path [PATH_MAX];
  BioseqPtr  bsp;
  Char       str [128];
  
  if (sep == NULL) return FALSE;
  
  VisitBioseqsInSep (sep, &bsplist, ListDeltaCallback);
  if (bsplist == NULL) return FALSE;
  
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp != NULL) {
    for (vnp = bsplist; vnp != NULL; vnp = vnp->next) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (bsp != NULL && bsp->id != NULL) {
        SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), str, PRINTID_REPORT, sizeof (str));
        fprintf (fp, "%s\n", str);
      }
    }
    FileClose (fp);
    LaunchGeneralTextViewer (path, "Delta Sequences");
    FileRemove (path);
    ValNodeFree (bsplist);
  }
  return TRUE;
}


extern void DeltaReport (IteM  i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  ShowDeltaReport(sep);
  ArrowCursor ();
  Update ();
}


static void FixLocusTagGeneXrefBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr        sfp, gene;
  SeqMgrFeatContext context;
  GeneRefPtr        grp;
  
  if (bsp == NULL) return;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (sfp->data.choice == SEQFEAT_GENE) continue;
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) continue;
    if (StringDoesHaveText (grp->locus) || StringHasNoText (grp->locus_tag)) continue;
    gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &context);
    if (gene != NULL && gene->data.choice == SEQFEAT_GENE
        && StringDoesHaveText (((GeneRefPtr)gene->data.value.ptrvalue)->locus)) {
      grp->locus = StringSave (((GeneRefPtr)gene->data.value.ptrvalue)->locus);
    }
  }
}


extern void FixLocusTagGeneXrefs (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  
  VisitBioseqsInSep (sep, NULL, FixLocusTagGeneXrefBioseqCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
  
  ArrowCursor ();
  Update ();
}


static SeqFeatPtr FindLastExonForCDS (BioseqPtr bsp, SeqFeatPtr cds)
{
  SeqFeatPtr        exon, last_exon = NULL;
  SeqMgrFeatContext cds_context, exon_context;
  Int4              cds_stop;
  
  if (bsp == NULL || cds == NULL) return NULL;
  
  cds_stop = SeqLocStop (cds->location);
  exon = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_IMP, FEATDEF_exon, &exon_context);
  while (exon != NULL && exon_context.left < cds_stop) {
    if (SeqMgrGetOverlappingCDS (exon->location, &cds_context) == cds) {
      if (SeqLocStrand (cds->location) == Seq_strand_minus) {
        return exon;
      } else {
        last_exon = exon;
      }
    }
    exon = SeqMgrGetNextFeature (bsp, exon, SEQFEAT_IMP, FEATDEF_exon, &exon_context);
  }
  return last_exon;
}


static Boolean DoesCDSEndWithStopCodon (SeqFeatPtr cds)
{
  ByteStorePtr bs;
  CharPtr      prot_str;
  Boolean      retval = FALSE;
  
  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }
  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
  if (bs == NULL) return FALSE;
  prot_str = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot_str == NULL) return FALSE;

  if (prot_str[StringLen (prot_str) - 1] == '*') {
    retval = TRUE;
  } else {
    retval = FALSE;
  }
  prot_str = MemFree (prot_str);
  return retval;
}

typedef struct lastexonfix {
  FILE *fp;
  Boolean make_partial;
} LastExonFixData, PNTR LastExonFixPtr;

static void FixLastExonLocCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr        cds, last_exon;
  SeqMgrFeatContext context;
  Int4              cds_stop, exon_stop;
  SeqIntPtr         sintp;
  LastExonFixPtr    lefp;
  Boolean           partial5, partial3;
  Boolean           was_changed;
  CharPtr           location;
  
  if (bsp == NULL || userdata == NULL) return;
  lefp = (LastExonFixPtr) userdata;
  
  for (cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       cds != NULL;
       cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, FEATDEF_CDS, &context)) {
    /* only extend the last codon if the coding region ends with a stop codon */
    if (!DoesCDSEndWithStopCodon(cds)) continue;
    last_exon = FindLastExonForCDS (bsp, cds);
    /* if last_exon location doesn't match end of CDS, extend it */
    if (last_exon != NULL && last_exon->location != NULL && last_exon->location->choice == SEQLOC_INT) {
      was_changed = FALSE;
      location = NULL;
      if (SeqLocStrand (cds->location) == Seq_strand_minus) {
        cds_stop = SeqLocStart (cds->location);
        exon_stop = SeqLocStart (last_exon->location);
        /* only extend the last exon by three base pairs to cover the last exon */
        if (cds_stop == exon_stop - 3) {
          location = SeqLocPrintUseBestID (last_exon->location);
          sintp = (SeqIntPtr) last_exon->location->data.ptrvalue;
          sintp->from = cds_stop;
          was_changed = TRUE;
        }
      } else {
        cds_stop = SeqLocStop (cds->location);
        exon_stop = SeqLocStop (last_exon->location);
        /* only extend the last exon by three base pairs to cover the last exon */
        if (cds_stop == exon_stop + 3) {
          location = SeqLocPrintUseBestID (last_exon->location);
          sintp = (SeqIntPtr) last_exon->location->data.ptrvalue;
          sintp->to = cds_stop;
          was_changed = TRUE;
        }
      }
      /* if the location changed, and the make_partial flag is set, make the last exon 3' partial */
      if (was_changed && lefp->make_partial) {
        CheckSeqLocForPartial (last_exon->location, &partial5, &partial3);
        SetSeqLocPartial (last_exon->location, partial5, TRUE);
        last_exon->partial = TRUE;
      }      
      if (location != NULL && lefp->fp != NULL) {
        fprintf (lefp->fp, "Exon at %s extended\n", location);
      }
      location = MemFree (location);
    }
  }
}

static void FixLastExonLoc (IteM i, Boolean make_partial)
{
  BaseFormPtr     bfp;
  SeqEntryPtr     sep;
  LastExonFixData lefd;
  Char            path [PATH_MAX];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  
  lefd.make_partial = make_partial;
  TmpNam (path);
  lefd.fp = FileOpen (path, "w");
  
  VisitBioseqsInSep (sep, &lefd, FixLastExonLocCallback);
  
  FileClose (lefd.fp);
  LaunchGeneralTextViewer (path, "Last Exons");
  FileRemove (path);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
  ArrowCursor ();
  Update ();
}

extern void FixLastExonLocNoPartial (IteM i)
{
  FixLastExonLoc (i, FALSE);
}


extern void FixLastExonLocMakePartial (IteM i)
{
  FixLastExonLoc (i, TRUE);
}


typedef struct removegenebyfeat
{
  FORM_MESSAGE_BLOCK
  DialoG  feature_select;
  DialoG  constraints;
  DialoG  accept_cancel;
  DoC     item_list;
  ButtoN  skip_suppressed_btn;
  ButtoN  skip_other_btn;
  PrompT  num_genes;
  
  Boolean    skip_suppressed;
  Boolean    skip_other;

  ValNodePtr feature_type_list;
  ValNodePtr feature_list;
  Int2       clicked;
  Boolean    dblClick;
  BaseFormPtr bfp;
} RemoveGeneByFeatData, PNTR RemoveGeneByFeatPtr;

static void AddGeneToList (SeqFeatPtr sfp, ValNodePtr PNTR list)
{
  ValNodePtr vnp;
  
  if (sfp == NULL || list == NULL) return;
  for (vnp = *list; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == sfp) return;
  }
  ValNodeAddPointer (list, OBJ_SEQFEAT, sfp);
}

static Boolean DoesGeneContainOtherFeatures (SeqFeatPtr gene, SeqFeatPtr one_feat)
{
  BioseqPtr         bsp;
  SeqLocPtr slp;
  Int4              tmp;
  Int4              gene_left, gene_right;
  SeqFeatPtr        other_feat;
  SeqMgrFeatContext fcontext;
  
  if (gene == NULL) return FALSE;
  
  for (slp = SeqLocFindNext (gene->location, NULL); slp != NULL; slp = slp->next) {      
    bsp = BioseqFindFromSeqLoc (gene->location);
    gene_left = SeqLocStart (slp);
    gene_right = SeqLocStop (slp);
    if (gene_left > gene_right) {
      tmp = gene_left;
      gene_left = gene_right;
      gene_right = tmp;
    }
    other_feat = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (other_feat != NULL && fcontext.left <= gene_right) {
      if (other_feat != gene && other_feat != one_feat
          && TestFeatOverlap (other_feat, gene, CONTAINED_WITHIN) > -1) {
        return TRUE;
      }
      other_feat = SeqMgrGetNextFeature (bsp, other_feat, 0, 0, &fcontext);
    }
  }
  return FALSE;
}

static void ReportGeneByFeatCallback (SeqFeatPtr sfp, Pointer userdata)
{
  RemoveGeneByFeatPtr  dlg;
  GeneRefPtr           grp;
  SeqFeatPtr           overlap_gene;
  ValNodePtr           vnp;
  SeqMgrFeatContext    fcontext;
  Boolean              is_suppressed;

  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  dlg = (RemoveGeneByFeatPtr) userdata;
  
  grp = SeqMgrGetGeneXref (sfp);
  is_suppressed = SeqMgrGeneIsSuppressed (grp);
  if (dlg->skip_suppressed && grp != NULL && is_suppressed) return;

  for (vnp = dlg->feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == sfp->idx.subtype) 
    {
      if (grp != NULL && !is_suppressed) {
        overlap_gene = SeqMgrGetGeneByLocusTag (BioseqFindFromSeqLoc(sfp->location), grp->locus_tag, &fcontext);
      } else {
        overlap_gene = SeqMgrGetOverlappingGene(sfp->location, &fcontext);
      }
      if (!dlg->skip_other || !DoesGeneContainOtherFeatures(overlap_gene, sfp)) {
        AddGeneToList(overlap_gene, &(dlg->feature_list));
      }
      return;
    }
  }
}

static Nlm_ParData removeGeneParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static Nlm_ColData removeGeneColFmt [3] = {{0, 5, 10, 0, NULL, 'l', 1,0,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,0,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,0,0,0, TRUE}};

static void PopulateRemoveGeneList (DoC doc, ValNodePtr feature_list)
{
  ValNodePtr        vnp;
  Int2              numItems;
  CharPtr           row_text;
  RecT              r;
  
  if (doc == NULL)
  {
    return;
  }
  Reset (doc);
  
  if (feature_list == NULL)
  {
    AppendText (doc, "No genes listed", NULL, NULL, programFont);
  }

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  removeGeneColFmt[0].pixWidth = 5 * stdCharWidth;
  removeGeneColFmt[1].pixWidth = (r.right - r.left - removeGeneColFmt[0].pixWidth) / 2;
  removeGeneColFmt[2].pixWidth = (r.right - r.left - removeGeneColFmt[0].pixWidth) / 2;
  
  vnp = feature_list;
  
  while (vnp != NULL)
  {
    row_text = GetDiscrepancyItemText (vnp);
    if (row_text != NULL)
    {
      if (vnp->choice == OBJ_SEQFEAT)
      {
        AppendText (doc, row_text, &removeGeneParFmt, removeGeneColFmt, programFont);
      }
      else
      {
        AppendText (doc, row_text, &removeGeneParFmt, NULL, programFont);
      }
      row_text = MemFree (row_text);
    }
    vnp = vnp->next;
  }
  GetDocParams (doc, &numItems, NULL);
  UpdateDocument (doc, 0, numItems);  
}

const CharPtr num_genes_fmt = "%d gene%s to be removed";


static void RemoveGeneByFeatChangeNotify (Pointer userdata)
{
  RemoveGeneByFeatPtr dlg;
  SeqEntryPtr         sep;
  Char                num_genes_txt[100];
  Int4                num_genes;
  
  dlg = (RemoveGeneByFeatPtr) userdata;
  if (dlg == NULL) 
  {
    return;
  }

  WatchCursor ();
  Update ();
    
  dlg->skip_suppressed = GetStatus (dlg->skip_suppressed_btn);
  dlg->skip_other = GetStatus (dlg->skip_other_btn);
  
  dlg->feature_type_list = ValNodeFree (dlg->feature_type_list);
  dlg->feature_type_list = (ValNodePtr) DialogToPointer (dlg->feature_select);
  sep = GetTopSeqEntryForEntityID(dlg->input_entityID);
  dlg->feature_list = ValNodeFree(dlg->feature_list);
  VisitFeaturesInSep (sep, dlg, ReportGeneByFeatCallback);
  PopulateRemoveGeneList (dlg->item_list, dlg->feature_list);
  
  num_genes = ValNodeLen (dlg->feature_list);
  sprintf (num_genes_txt, num_genes_fmt, num_genes, num_genes == 1 ? "" : "s");
  SetTitle (dlg->num_genes, num_genes_txt); 
  
  if (dlg->feature_list == NULL) 
  {
    DisableAcceptCancelDialogAccept (dlg->accept_cancel);
  } 
  else 
  {
    EnableAcceptCancelDialogAccept (dlg->accept_cancel);
  }
  ArrowCursor ();
  Update ();
} 

static void RemoveGeneByFeatChangeBtn (ButtoN b)
{
  RemoveGeneByFeatPtr dlg;

  dlg = (RemoveGeneByFeatPtr) GetObjectExtra (b);
  RemoveGeneByFeatChangeNotify (dlg);
}

static void RemoveGeneByFeatClear (Pointer data)
{
  RemoveGeneByFeatPtr dlg;

  dlg = (RemoveGeneByFeatPtr) data;
  if (dlg == NULL) return;
 
  PointerToDialog (dlg->feature_select, NULL);
  PointerToDialog (dlg->constraints, NULL);
}

static Boolean RemoveGeneByFeatAction (Pointer userdata)
{
  RemoveGeneByFeatPtr dlg;
  SeqEntryPtr         sep;
  ValNodePtr          vnp;
  SeqFeatPtr          sfp;
  
  if (userdata == NULL) return FALSE;
  
  dlg = (RemoveGeneByFeatPtr) userdata;
  
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  if (sep == NULL) return FALSE;

  for (vnp = dlg->feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL) {
      sfp->idx.deleteme = TRUE;
    }
  }
  DeleteMarkedObjects (dlg->input_entityID, 0, NULL);
  dlg->feature_list = ValNodeFree(dlg->feature_list);
  
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);  
  Update ();
  RemoveGeneByFeatChangeNotify(dlg);   
  return TRUE;
}

static SeqFeatPtr GetSelectedRemoveGene(ValNodePtr feature_list, Int2 item)
{
    while (feature_list != NULL && item > 1) {
        feature_list = feature_list->next;
        item--;
    }
    if (feature_list == NULL) {
      return NULL;
    } else {
      return feature_list->data.ptrvalue;
    }
}

static void ClickRemoveGene (DoC d, PoinT pt)

{
  Int2                      item, numItems;
  Int2                      row;
  Int2                col, last_selected;
  RemoveGeneByFeatPtr dlg;
  SeqFeatPtr          sfp;
  BioseqPtr           bsp;

  dlg = (RemoveGeneByFeatPtr)GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && dlg->clicked == item) {
      dlg->dblClick = dblClick;
    } else {
      dlg->dblClick = FALSE;
    }
    last_selected = dlg->clicked;
    dlg->clicked = 0;
    if (item > 0 && row > 0) {
      dlg->clicked = item;
    }
        
    if (item != last_selected)
    {
      GetDocParams (d, &numItems, NULL);
      UpdateDocument (d, 0, numItems);
    }

    if (item > 0 && row > 0)
    {
      sfp = GetSelectedRemoveGene(dlg->feature_list, item);
      if (sfp != NULL) {
        if (dblClick)
        {
          GatherProcLaunch (OMPROC_EDIT, FALSE, sfp->idx.entityID, sfp->idx.itemID,
                            OBJ_SEQFEAT, 0, 0, OBJ_SEQFEAT, 0);
        }
        else
        {
          /* need to scroll to item */
          bsp = BioseqFindFromSeqLoc (sfp->location);
          SetBioseqViewTargetByBioseq (dlg->bfp, bsp);
          ObjMgrSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
        } 
      }
    }
  }
}

static void DrawRemoveGene (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RemoveGeneByFeatPtr dlg;
  RecT                rct;

  dlg = (RemoveGeneByFeatPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    /* draw selection */
    if (item == dlg->clicked) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
  }
}

static void CleanupRemoveGeneForm (GraphiC g, VoidPtr data)

{
  RemoveGeneByFeatPtr dlg;

  dlg = (RemoveGeneByFeatPtr) data;
  if (dlg != NULL) {
    dlg->feature_list = ValNodeFree (dlg->feature_list);
    dlg->feature_type_list = ValNodeFreeData (dlg->feature_type_list);
  }
  StdCleanupFormProc (g, data);
}

extern void RemoveGeneByUnderlyingFeatureType(IteM i)
{
  BaseFormPtr       bfp;
  RemoveGeneByFeatPtr dlg;
  WindoW            w;
  GrouP             h, k, g, feat_opts;
  SeqEntryPtr       sep;
  PrompT            p;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (RemoveGeneByFeatPtr) MemNew (sizeof (RemoveGeneByFeatData));
  if (dlg == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove Genes By Underlying Feature", StdCloseWindowProc);
  SetObjectExtra (w, dlg, CleanupRemoveGeneForm);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  dlg->bfp = bfp;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  k = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (k, 10, 10);

  dlg->item_list = DocumentPanel (k, stdCharWidth * 30 + 5, stdLineHeight * 20);
  SetObjectExtra (dlg->item_list, dlg, NULL);
  SetDocAutoAdjust (dlg->item_list, FALSE);
  SetDocProcs (dlg->item_list, ClickRemoveGene, NULL, NULL, NULL);
  SetDocShade (dlg->item_list, DrawRemoveGene, NULL, NULL, NULL);
    
  g = HiddenGroup (k, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  p = StaticPrompt (g, "Remove genes that contain features of type:", 0, dialogTextHeight, systemFont, 'l');
  dlg->feature_select =  FeatureSelectionDialogEx (g, TRUE, sep,
                                                 RemoveGeneByFeatChangeNotify, 
                                                 dlg);
  feat_opts = HiddenGroup (g, -1, 0, NULL);
  dlg->skip_suppressed_btn = CheckBox (feat_opts, "Skip features with suppressed gene xrefs", RemoveGeneByFeatChangeBtn);
  SetObjectExtra (dlg->skip_suppressed_btn, dlg, NULL);
  SetStatus (dlg->skip_suppressed_btn, TRUE);
  dlg->skip_suppressed = TRUE;
  dlg->skip_other_btn = CheckBox (feat_opts, "Skip genes that contain other features", RemoveGeneByFeatChangeBtn);
  SetObjectExtra (dlg->skip_other_btn, dlg, NULL);
  SetStatus (dlg->skip_other_btn, TRUE);
  dlg->skip_other = TRUE;
  
  dlg->num_genes = StaticPrompt (g, "0 genes to be removed", 200, dialogTextHeight, systemFont, 'l');

  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) dlg->feature_select, (HANDLE) feat_opts, NULL);

  dlg->accept_cancel = AcceptCancelDialog (h, RemoveGeneByFeatAction, NULL, RemoveGeneByFeatClear, NULL, (Pointer)dlg, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) dlg->accept_cancel, NULL);
                                
  RemoveGeneByFeatChangeNotify(dlg);                                
  Show (w);
}


typedef struct intronadjust {
  Uint1 featdef;
  ValNodePtr bad_feature_list;
} IntronAdjustData, PNTR IntronAdjustPtr;


static void RemoveIntronLocationsFromFeatureCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr sfp, intron;
  SeqMgrFeatContext sfp_context, intron_context;
  IntronAdjustPtr iap;
  SeqLocPtr new_loc, tmp_loc;
  Boolean   partial5, partial3;
  
  if (bsp == NULL || userdata == NULL) return;
  
  iap = (IntronAdjustPtr)userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, iap->featdef, &sfp_context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, iap->featdef, &sfp_context)) {
    for (intron = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_IMP, FEATDEF_intron, &intron_context);
         intron != NULL && intron_context.left < sfp_context.right;
         intron = SeqMgrGetNextFeature (bsp, intron, SEQFEAT_IMP, FEATDEF_intron, &intron_context)) {
       if (intron_context.right < sfp_context.left) continue;
       CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
       tmp_loc = (SeqLocPtr) AsnIoMemCopy ((Pointer) sfp->location,
                                     (AsnReadFunc) SeqLocAsnRead,
                                     (AsnWriteFunc) SeqLocAsnWrite);
       tmp_loc = SeqLocSubtract (tmp_loc, intron->location);
       new_loc = SeqLocMerge (bsp, NULL, tmp_loc, FALSE, FALSE, FALSE);
       tmp_loc = SeqLocFree (tmp_loc);
       SetSeqLocPartial (new_loc, partial5, partial3);
       
       if (new_loc == NULL) {
         ValNodeAddPointer (&iap->bad_feature_list, OBJ_SEQFEAT, sfp);
       } else {
         sfp->location = SeqLocFree (sfp->location);
         sfp->location = new_loc;
       }
    }
  }  
}


typedef void (*ReportClickableListCallback) PROTO ((ValNodePtr));
typedef ValNodePtr (*RefreshClickableListCallback) PROTO ((Uint2));

typedef struct clickableitemlistwindow {
  FORM_MESSAGE_BLOCK
  ValNodePtr      item_list;
  DialoG          clickable_list;
  BaseFormPtr     bfp;
  ReportClickableListCallback report_func;
  RefreshClickableListCallback refresh_func;
  
} ClickableItemListWindowData, PNTR ClickableItemListWindowPtr;

static void CleanupClickableItemListWindow (GraphiC g, VoidPtr data)

{
  ClickableItemListWindowPtr drfp;

  drfp = (ClickableItemListWindowPtr) data;
  if (drfp != NULL) {
    drfp->item_list = FreeClickableList (drfp->item_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static void ClickableItemListWindowMessage (ForM f, Int2 mssg)

{
  ClickableItemListWindowPtr drfp;

  drfp = (ClickableItemListWindowPtr) GetObjectExtra (f);
  if (drfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_COPY :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        drfp->item_list = FreeClickableList (drfp->item_list);
        PointerToDialog (drfp->clickable_list, NULL);
        break;
      default :
        if (drfp->appmessage != NULL) {
          drfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Int2 LIBCALLBACK ClickableItemListWindowMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  ClickableItemListWindowPtr drfp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  drfp = (ClickableItemListWindowPtr) omudp->userdata.ptrvalue;
  if (drfp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (drfp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (drfp->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (drfp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (drfp->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}

#ifndef WIN_MAC
extern void CreateStdValidatorFormMenus (WindoW w);
#endif


static void ReportClickableList (ButtoN b)
{
  ClickableItemListWindowPtr drfp;
  ValNodePtr                 item_list;

  drfp = (ClickableItemListWindowPtr) GetObjectExtra (b);
  if (drfp == NULL || drfp->report_func == NULL || drfp->item_list == NULL) return;

  item_list = DialogToPointer (drfp->clickable_list);

  if (item_list == NULL) {
    if (ANS_YES == Message (MSG_YN, "You have not selected any items - report all?"))
    {
      item_list = drfp->item_list;
    }
    else
    {
      return;
    }

  }

  (drfp->report_func) (item_list);

  if (item_list != drfp->item_list) {
    item_list = ValNodeFree (item_list);
  }
}


static void RefreshClickableList (ButtoN b)
{
  ClickableItemListWindowPtr drfp;

  drfp = (ClickableItemListWindowPtr) GetObjectExtra (b);
  if (drfp == NULL || drfp->refresh_func == NULL) return;
  
  PointerToDialog (drfp->clickable_list, NULL);
  drfp->item_list = FreeClickableList (drfp->item_list);

  drfp->item_list = (drfp->refresh_func)(drfp->input_entityID);
  PointerToDialog (drfp->clickable_list, drfp->item_list);  
}


static void 
ShowClickableItemListEx 
(ValNodePtr  clickable_list,
 BaseFormPtr bfp,
 CharPtr     win_title,
 CharPtr     label1,
 CharPtr     label2,
 ReportClickableListCallback report_func,
 RefreshClickableListCallback refresh_func)
{
  ClickableItemListWindowPtr drfp;
  GrouP                    h;
  GrouP                    c;
  WindoW                   w;
  OMUserDataPtr            omudp;
  SeqEntryPtr              sep;
  ButtoN                   b;

  if (bfp == NULL) return;
    
  drfp = (ClickableItemListWindowPtr) MemNew (sizeof (ClickableItemListWindowData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, win_title, StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupClickableItemListWindow);
  drfp->form = (ForM) w;
  drfp->formmessage = ClickableItemListWindowMessage;

  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = ClickableItemListWindowMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif

  drfp->report_func = report_func;
  drfp->refresh_func = refresh_func;
    
  drfp->item_list = clickable_list;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = CreateClickableListDialog (h, label1, label1,
                                                    ScrollToDiscrepancyItem, EditDiscrepancyItem, bfp,
                                                    GetDiscrepancyItemText);

  PointerToDialog (drfp->clickable_list, drfp->item_list);                                                    
  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
    
  if (drfp->report_func != NULL) {
    b = PushButton (c, "Make Report", ReportClickableList);
    SetObjectExtra (b, drfp, NULL);
  }
  if (drfp->refresh_func != NULL) {
    b = PushButton (c, "Refresh", RefreshClickableList);
    SetObjectExtra (b, drfp, NULL);
  }
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);
}

extern void ShowClickableItemList (ValNodePtr clickable_list, BaseFormPtr bfp, CharPtr win_title, CharPtr label1, CharPtr label2)
{
  ShowClickableItemListEx (clickable_list, bfp, win_title, label1, label2, NULL, NULL);
}

static void RemoveIntronLocationsFromFeature (IteM i, Uint1 featdef)
{
  BaseFormPtr             bfp;
  SeqEntryPtr             sep;
  IntronAdjustData        iad;
  ValNodePtr              item_list = NULL;
  ClickableItemPtr        cip;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  WatchCursor ();
  Update ();

  iad.featdef = featdef;
  iad.bad_feature_list = NULL;
  
  VisitBioseqsInSep (sep, &iad, RemoveIntronLocationsFromFeatureCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
  
  if (iad.bad_feature_list != NULL) {
    cip = NewClickableItem (0, "%d features are completely covered by introns", iad.bad_feature_list);
    item_list = ValNodeAddPointer (&item_list, 0, cip);
    ShowClickableItemList (item_list, bfp, "Failed Intron Location Removal", "", "Unremoved Features");
  } 
}


extern void RemoveIntronLocationsFromCDS (IteM i)
{
  RemoveIntronLocationsFromFeature(i, FEATDEF_CDS);
}


extern void RemoveIntronLocationsFromrRNA (IteM i)
{
  RemoveIntronLocationsFromFeature(i, FEATDEF_rRNA);
}


extern void RemoveIntronLocationsFrommRNA (IteM i)
{
  RemoveIntronLocationsFromFeature(i, FEATDEF_mRNA);
}


extern void RemoveIntronLocationsFromtRNA (IteM i)
{
  RemoveIntronLocationsFromFeature(i, FEATDEF_tRNA);
}


extern ValNodePtr ChooseFeaturesForConversion (ValNodePtr clickable_list, BaseFormPtr bfp, CharPtr label1, CharPtr label2)
{
  ClickableItemListWindowData windata;
  GrouP                    h;
  GrouP                    c;
  ButtoN                   b;
  WindoW                   w;
  Boolean                  done;
  SeqEntryPtr              sep;
  ModalAcceptCancelData    acd;
  ValNodePtr               vnp_list, selected_feats = NULL;
  ClickableItemPtr         cip;

  if (bfp == NULL) return NULL;

  MemSet (&windata, 0, sizeof (ClickableItemListWindowData)); 
  
  windata.bfp = bfp;
  windata.input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID (windata.input_entityID);
  w = MovableModalWindow (-50, -33, -10, -10,
                          "Choose Features for Conversion",
                          StdCloseWindowProc);
  SetObjectExtra (w, &windata, NULL);
  windata.form = (ForM) w;
    
  windata.item_list = clickable_list;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
 
  windata.clickable_list = CreateClickableListDialog (h, label1, label1,
                                                    ScrollToDiscrepancyItem, EditDiscrepancyItem, bfp,
                                                    GetDiscrepancyItemText);

  PointerToDialog (windata.clickable_list, windata.item_list);                                                    
    
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Convert Selected Features", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) windata.clickable_list, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);
  Select (w);
  done = FALSE;
  while (!done)
  { 
    acd.accepted = FALSE;
    acd.cancelled = FALSE;
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
  
    if (acd.cancelled)
    {
      done = TRUE;
    }
    else
    {
      for (vnp_list = windata.item_list; vnp_list != NULL; vnp_list = vnp_list->next) {
        cip = (ClickableItemPtr) vnp_list->data.ptrvalue;
        if (cip->chosen) {
          ValNodeLink (&selected_feats, cip->item_list);
          cip->item_list = NULL;
        }
      }
      if (selected_feats == NULL) {
        Message (MSG_ERROR, "You didn't select any features to convert!  Choose features or hit Cancel");
      } else {
        done = TRUE;
      }
    }
  }
  Remove (w);
  windata.item_list = FreeClickableList (windata.item_list);
  return selected_feats;
}

static void RemoveBadPubsDescCallback (SeqDescPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovp;
  
  if (sdp == NULL || sdp->choice != Seq_descr_pub || sdp->extended == 0) {
    return;
  }
  if (IsPubContentBad (sdp->data.ptrvalue)) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.deleteme = TRUE;
  }
}

static void RemoveBadPubsFeatCallback (SeqFeatPtr sfp, Pointer userdata)
{  
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) {
    return;
  }
  if (IsPubContentBad (sfp->data.value.ptrvalue)) {
    sfp->idx.deleteme = TRUE;
  }
}

extern void RemoveBadPubs (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  VisitDescriptorsInSep (sep, NULL, RemoveBadPubsDescCallback);
  VisitFeaturesInSep (sep, NULL, RemoveBadPubsFeatCallback);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}



/* VecScreen Tool */
typedef struct vecscreentool {
  FORM_MESSAGE_BLOCK
  ValNodePtr      item_list;
  DialoG          clickable_list;
  BaseFormPtr     bfp;
  GrouP           trim_options_grp;
  ButtoN          add_citsub_btn;  

  SeqEntryPtr     top_sep;
  LogInfoPtr      lip;
  Int4            trim_option;
  ValNodePtr      trimmed_bsps;
  
} VecScreenToolData, PNTR VecScreenToolPtr;

static void CleanupVecScreenTool (GraphiC g, VoidPtr data)

{
  VecScreenToolPtr drfp;

  drfp = (VecScreenToolPtr) data;
  if (drfp != NULL) {
    drfp->item_list = FreeClickableList (drfp->item_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static void VecScreenToolMessage (ForM f, Int2 mssg)

{
  VecScreenToolPtr drfp;

  drfp = (VecScreenToolPtr) GetObjectExtra (f);
  if (drfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_COPY :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        drfp->item_list = FreeClickableList (drfp->item_list);
        PointerToDialog (drfp->clickable_list, NULL);
        break;
      default :
        if (drfp->appmessage != NULL) {
          drfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Int2 LIBCALLBACK VecScreenToolMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  VecScreenToolPtr drfp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  drfp = (VecScreenToolPtr) omudp->userdata.ptrvalue;
  if (drfp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (drfp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (drfp->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (drfp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (drfp->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}


static CharPtr GetMatchType (SeqFeatPtr sfp)
{
  GBQualPtr gbq;
  
  if (sfp == NULL) return NULL;
  gbq = sfp->qual;
  while (gbq != NULL) {
    if (StringCmp (gbq->qual, "phenotype") == 0) {
      if (StringCmp (gbq->val, "Suspect origin") == 0) {
        return "suspect";
      } else if (StringCmp (gbq->val, "Weak match") == 0) {
        return "weak";
      } else if (StringCmp (gbq->val, "Strong match") == 0) {
        return "strong";
      } else if (StringCmp (gbq->val, "Moderate match") == 0) {
        return "moderate";
      } else {
        return gbq->val;
      }
    }
    gbq = gbq->next;
  }
  return NULL;
}


static SeqLocPtr GetFeatureListLoc (ValNodePtr feat_list, VecScreenToolPtr vstp)
{
  ValNodePtr vnp;
  BioseqPtr  bsp;
  SeqFeatPtr sfp;
  Int4       loc_left = -1, loc_right = -1, tmp, start, stop;
    
  vnp = feat_list;
  while (vnp != NULL) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      start = SeqLocStart (sfp->location);
      stop = SeqLocStop (sfp->location);
      if (start > stop) {
        tmp = stop;
        stop = start;
        start = tmp;
      }
      if (loc_left < 0) {
        loc_left = start;
        loc_right = stop;
        bsp = BioseqFindFromSeqLoc (sfp->location);
      } else {
        loc_left = MIN (loc_left, start);
        loc_right = MAX (loc_right, stop);
        if (BioseqFindFromSeqLoc (sfp->location) != bsp) {
          /* quit - this is bad */
          return NULL;
        }
      }
    }
    vnp = vnp->next;
  }
  
  if (bsp == NULL) {
    /* quit - this is bad */
    return NULL;
  }

  if (vstp != NULL) {          
    if (loc_left != 0 && loc_right != bsp->length - 1) {
      /* internal - fix location */
        
      if (vstp->trim_option == 1) {
        /* trim to closest end */
        if (loc_left - 0 < bsp->length - 1 - loc_right) {
          loc_left = 0;
        } else {
          loc_right = bsp->length - 1;
        }
      } else if (vstp->trim_option == 2) {
        /* trim to 5' end */
        loc_left = 0;
      } else {
        /* trim to 3' end */
        loc_right = bsp->length - 1;
      }
    }
  }
  return SeqLocIntNew (loc_left, loc_right, Seq_strand_plus, SeqIdDup (bsp->id));        
}


extern void CalculateVectorDescription (ClickableItemPtr cip)
{
  ValNodePtr vnp;
  Int4       len, loc_left, loc_right;
  SeqLocPtr  slp;
  BioseqPtr  bsp;
  CharPtr    internal_fmt = "Internal: %d from 5' end, %d from 3' end", cp;
  Char       str[256];
  
  if (cip == NULL || cip->item_list == NULL) {
    return;
  }
  
  slp = GetFeatureListLoc (cip->item_list, NULL);
  if (slp == NULL) return;
  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  bsp = BioseqFindFromSeqLoc (slp);

  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), str, PRINTID_REPORT, sizeof (str));
  len = StringLen (str) + 2;

  if (loc_left == 0) {
    len += 8;
  } else if (loc_right == bsp->length - 1) {
    len += 8;
  } else {
    len += StringLen (internal_fmt) + 30;
  }
  for (vnp = cip->item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      cp = GetMatchType (vnp->data.ptrvalue);
      if (!StringHasNoText (cp)) {
        len += StringLen (cp) + 1;
      }
    }
  }
  cip->description = MemFree (cip->description);
  cip->description = (CharPtr) MemNew (len * sizeof (Char));
  if (loc_left == 0) {
    StringCpy (cip->description, "5' End;");
  } else if (loc_right == bsp->length - 1) {
    StringCpy (cip->description, "3' End;");
  } else {
    sprintf (cip->description, internal_fmt, loc_left, bsp->length - loc_right);
  }
  StringCat (cip->description, str);
  
  for (vnp = cip->item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      cp = GetMatchType (vnp->data.ptrvalue);
      if (!StringHasNoText (cp)) {
        StringCat (cip->description, ";");
        StringCat (cip->description, cp);
      }
    }
  }
  slp = SeqLocFree (slp);
}

static void GetVectorContaminationList (BioseqPtr bsp, Pointer userdata)
{
  Boolean          isvector;
  GBQualPtr        gbq;
  ClickableItemPtr cip = NULL;
  SeqFeatPtr       sfp;
  SeqMgrFeatContext fcontext;
  Int4              last_right;
  
  ValNodePtr PNTR feat_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || feat_list == NULL) return;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_IMP, FEATDEF_misc_feature, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_IMP, FEATDEF_misc_feature, &fcontext)) {         
    for (isvector = FALSE, gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringCmp (gbq->qual, "standard_name") == 0) {
        if (StringCmp (gbq->val, "Vector Contamination") == 0) {
          isvector = TRUE;
        }
      }
    }
    if (isvector) {
      if (cip != NULL && last_right >= fcontext.left - 1) {
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
        last_right = fcontext.right;
      } else {        
        /* calculate description for previous list */
        CalculateVectorDescription (cip);
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        if (cip != NULL) {
          cip->clickable_item_type = 0;
          cip->callback_func = NULL;
          cip->datafree_func = NULL;
          cip->callback_data = NULL;
          cip->item_list = NULL;
          cip->subcategories = NULL;
          cip->expanded = FALSE;
          cip->level = 0;
          cip->description = NULL;
  
          ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);

          ValNodeAddPointer (feat_list, 0, cip);
          last_right = fcontext.right;
        }
      }
    }
  }
  /* calculate description for last item */
  CalculateVectorDescription (cip);
}


static void VecScreenToolSelectAll (ValNodePtr vnp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      cip->chosen = TRUE;
      VecScreenToolSelectAll (cip->subcategories);
    }
    vnp = vnp->next;
  }
  
}


static void RedrawVecScreenTool (VecScreenToolPtr vstp)
{
  RecT     r;

  ObjectRect (vstp->clickable_list, &r);
  InsetRect (&r, -1, -1);
  InvalRect (&r);
}


static void VecScreenToolSelectAllBtn (ButtoN b)
{
  VecScreenToolPtr vstp;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  VecScreenToolSelectAll (vstp->item_list);  
  RedrawVecScreenTool (vstp);
}


static void VecScreenToolUnselectAll (ValNodePtr vnp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      if (cip->chosen)
      {
        cip->chosen = FALSE;
      }
      VecScreenToolUnselectAll (cip->subcategories);
    }
    vnp = vnp->next;
  }
  
}


static void VecScreenToolUnselectAllBtn (ButtoN b)
{
  VecScreenToolPtr vstp;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  VecScreenToolUnselectAll (vstp->item_list);  
  RedrawVecScreenTool (vstp);  
}


static void VecScreenToolUnselectInternal (ValNodePtr vnp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      if (cip->chosen)
      {
        if (StringNCmp (cip->description, "Internal", 8) == 0) {
          cip->chosen = FALSE;
          VecScreenToolUnselectAll (cip->subcategories);
        }
      } else {
        VecScreenToolUnselectInternal (cip->subcategories);
      }
    }
    vnp = vnp->next;
  }  
}

static void VecScreenToolUnselectInternalBtn (ButtoN b)
{
  VecScreenToolPtr vstp;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;
  VecScreenToolUnselectInternal (vstp->item_list);
  RedrawVecScreenTool (vstp);  
}


static void TrimFeatureList (ValNodePtr feat_list, VecScreenToolPtr vstp)
{
  BioseqPtr  bsp = NULL;
  SeqLocPtr  delete_loc;
  Int4       len, loc_left = -1, loc_right = -1;
  
  /* create location to trim */
  delete_loc = GetFeatureListLoc (feat_list, vstp);
  if (delete_loc == NULL) return;
  bsp = BioseqFindFromSeqLoc (delete_loc);
  loc_left = SeqLocStart (delete_loc);
  loc_right = SeqLocStop (delete_loc);
  len = SeqLocLen (delete_loc);
      
  if (loc_left == 0) {
    TrimQualityScores (bsp, len, TRUE);
  } else if (loc_right == bsp->length - 1) {
    TrimQualityScores (bsp, len, FALSE);
  }
      
  if (vstp->top_sep != NULL) {
    SeqEntryExplore (vstp->top_sep, (Pointer) delete_loc, SeqAlignDeleteByLocCallback);
  }
  SeqDeleteByLoc (delete_loc, TRUE, FALSE);  
  delete_loc = SeqLocFree (delete_loc);  
  ValNodeAddPointer (&(vstp->trimmed_bsps), OBJ_BIOSEQ, bsp);
}


static void TrimSelected (ValNodePtr vnp, Boolean do_all, VecScreenToolPtr vstp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      if (cip->chosen || do_all)
      {
        if (cip->expanded)
        {
          TrimSelected (cip->subcategories, TRUE, vstp);
        }
        else
        {
          TrimFeatureList (cip->item_list, vstp);
        }
      }
      else if (cip->expanded)
      {
        TrimSelected (cip->subcategories, FALSE, vstp);
      }
    }
    vnp = vnp->next;
  }
}

static void CollectFeatureList (ValNodePtr feat_list, VecScreenToolPtr vstp, ValNodePtr PNTR loc_list)
{
  SeqLocPtr  delete_loc;
  
  /* create location to trim */
  delete_loc = GetFeatureListLoc (feat_list, vstp);
  if (delete_loc == NULL) return;
  ValNodeAddPointer (loc_list, 0, delete_loc);
}


static void CollectSelected (ValNodePtr vnp, Boolean do_all, VecScreenToolPtr vstp, ValNodePtr PNTR loc_list)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      if (cip->chosen || do_all)
      {
        if (cip->expanded)
        {
          CollectSelected (cip->subcategories, TRUE, vstp, loc_list);
        }
        else
        {
          CollectFeatureList (cip->item_list, vstp, loc_list);
        }
      }
      else if (cip->expanded)
      {
        CollectSelected (cip->subcategories, FALSE, vstp, loc_list);
      }
    }
    vnp = vnp->next;
  }
}

static int LIBCALLBACK SortVectorLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqLocPtr   slp1, slp2;
  BioseqPtr   bsp1, bsp2;
  Int4        start1, start2;
  Char        str1[256], str2[256];
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      slp1 = vnp1->data.ptrvalue;
      slp2 = vnp2->data.ptrvalue;
      
      bsp1 = BioseqFindFromSeqLoc (slp1);
      bsp2 = BioseqFindFromSeqLoc (slp2);
      SeqIdWrite (SeqIdFindBest (bsp1->id, SEQID_GENBANK), str1, PRINTID_REPORT, sizeof (str1));
      SeqIdWrite (SeqIdFindBest (bsp2->id, SEQID_GENBANK), str2, PRINTID_REPORT, sizeof (str2));
      rval = StringCmp (str1, str2);
      if (rval == 0) {
        start1 = SeqLocStart (slp1);
        start2 = SeqLocStart (slp2);
        if (start1 > start2) {
          rval = 1;
        } else if (start1 < start2) {
          rval = -1;
        }
      }
    }
  }
  return rval;
}

static int LIBCALLBACK SortClickableItemLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  ClickableItemPtr cip1, cip2;
  SeqFeatPtr       sfp1, sfp2;
  SeqLocPtr   slp1, slp2;
  BioseqPtr   bsp1, bsp2;
  Int4        start1, start2;
  Char        str1[256], str2[256];
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {      
      cip1 = vnp1->data.ptrvalue;
      cip2 = vnp2->data.ptrvalue;
      
      if (cip1 != NULL && cip2 != NULL 
          && cip1->item_list != NULL && cip2->item_list != NULL
          && cip1->item_list->choice == OBJ_SEQFEAT
          && cip2->item_list->choice == OBJ_SEQFEAT) {
        sfp1 = (SeqFeatPtr) cip1->item_list->data.ptrvalue;
        sfp2 = (SeqFeatPtr) cip2->item_list->data.ptrvalue;
        if (sfp1 != NULL && sfp1->location != NULL
            && sfp2 != NULL && sfp2->location != NULL) {
          slp1 = sfp1->location;
          slp2 = sfp2->location;
          bsp1 = BioseqFindFromSeqLoc (slp1);
          bsp2 = BioseqFindFromSeqLoc (slp2);
          SeqIdWrite (SeqIdFindBest (bsp1->id, SEQID_GENBANK), str1, PRINTID_REPORT, sizeof (str1));
          SeqIdWrite (SeqIdFindBest (bsp2->id, SEQID_GENBANK), str2, PRINTID_REPORT, sizeof (str2));
          rval = StringCmp (str1, str2);
          if (rval == 0) {
            start1 = SeqLocStart (slp1);
            start2 = SeqLocStart (slp2);
            if (start1 > start2) {
              rval = 1;
            } else if (start1 < start2) {
              rval = -1;
            }
          }
        }
      }
    }
  }
  return rval;
}


static void LogSelected (VecScreenToolPtr vstp)
{
  ValNodePtr loc_list = NULL, vnp;
  
  if (vstp == NULL) return;
  
  CollectSelected (vstp->item_list, FALSE, vstp, &loc_list);
  loc_list = ValNodeSort (loc_list, SortVectorLoc);
    
  for (vnp = loc_list; vnp != NULL; vnp = vnp->next) {
    LogTrimmedLocation (vstp->lip, vnp->data.ptrvalue);
    vnp->data.ptrvalue = SeqLocFree (vnp->data.ptrvalue);
  }
  loc_list = ValNodeFree (loc_list);
}


static ValNodePtr SortVecScreenList (ValNodePtr item_list)
{
  ClickableItemPtr cip;
  ValNodePtr list_5 = NULL, list_3 = NULL, list_internal = NULL;
  ValNodePtr list_5_last = NULL, list_3_last = NULL, list_internal_last = NULL;
  ValNodePtr next_item, last_item;
  
  while (item_list != NULL) {
    next_item = item_list->next;
    item_list->next = NULL;
    cip = (ClickableItemPtr) item_list->data.ptrvalue;
    if (StringNCmp (cip->description, "5' End", 6) == 0) {
      if (list_5_last == NULL) {
        list_5 = item_list;
      } else {
        list_5_last->next = item_list;
      }
      list_5_last = item_list;
    } else if (StringNCmp (cip->description, "3' End", 6) == 0) {
      if (list_3_last == NULL) {
        list_3 = item_list;
      } else {
        list_3_last->next = item_list;
      }
      list_3_last = item_list;
    } else {
      if (list_internal_last == NULL) {
        list_internal = item_list;
      } else {
        list_internal_last->next = item_list;
      }
      list_internal_last = item_list;
    }
    item_list = next_item;
  }
  
  item_list = list_internal;
  last_item = list_internal_last;
  
  if (list_5 != NULL) {
    if (item_list == NULL) {
      item_list = list_5;
      last_item = list_5_last;
    } else {
      last_item->next = list_5;
      last_item = list_5_last;
    }
  }
  
  if (list_3 != NULL) {
    if (item_list == NULL) {
      item_list = list_3;
      last_item = list_3_last;
    } else {
      last_item->next = list_3;
      last_item = list_3_last;
    }
  }
  
  return item_list;  
}


static void RefreshVecScreenList (VecScreenToolPtr vstp)
{
  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = FreeClickableList (vstp->item_list);

  VisitBioseqsInSep (vstp->top_sep, &(vstp->item_list), GetVectorContaminationList);
  vstp->item_list = SortVecScreenList (vstp->item_list);
  VecScreenToolSelectAll (vstp->item_list);
  VecScreenToolUnselectInternal (vstp->item_list);
  PointerToDialog (vstp->clickable_list, vstp->item_list);  
}


static void TrimSelectedBtn (ButtoN b)
{
  VecScreenToolPtr vstp;
  ValNodePtr       vnp;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;
  WatchCursor();
  Update();
  vstp->trimmed_bsps = NULL;
  vstp->trim_option = GetValue (vstp->trim_options_grp);
  vstp->lip = OpenLog ("Trimmed Locations");
  LogSelected (vstp);
  TrimSelected (vstp->item_list, FALSE, vstp);
  CloseLog (vstp->lip);
  vstp->lip = FreeLog (vstp->lip);

  /* add cit-subs */
  if (GetStatus (vstp->add_citsub_btn)) {
    for (vnp = vstp->trimmed_bsps; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == OBJ_BIOSEQ) {
        AddCitSubToUpdatedSequence (vnp->data.ptrvalue, vstp->input_entityID, kIndexerUpdateVecScreenText);
      }
    }
  }
  vstp->trimmed_bsps = ValNodeFree (vstp->trimmed_bsps);
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = FreeClickableList (vstp->item_list);
  DeleteMarkedObjects (vstp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);

  RefreshVecScreenList(vstp);
  RedrawVecScreenTool (vstp);  
  if (vstp->item_list == NULL) {
    Remove (vstp->form);
  }
  ArrowCursor();
  Update();
}


static void RemoveSelectedFeatureList (ValNodePtr vnp)
{
  SeqFeatPtr sfp;
  
  while (vnp != NULL) {
    if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      sfp->idx.deleteme = TRUE;
    }
    vnp = vnp->next;
  }
}


static void RemoveSelected (ValNodePtr vnp, Boolean do_all, VecScreenToolPtr vstp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL)
    {
      if (cip->chosen || do_all)
      {
        if (cip->expanded)
        {
          RemoveSelected (cip->subcategories, TRUE, vstp);
        }
        else
        {
          RemoveSelectedFeatureList (cip->item_list);
        }
      }
      else if (cip->expanded)
      {
        RemoveSelected (cip->subcategories, FALSE, vstp);
      }
    }
    vnp = vnp->next;
  }
}


static void RemoveSelectedBtn (ButtoN b)
{
  VecScreenToolPtr vstp;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;
  WatchCursor();
  Update();
  
  RemoveSelected (vstp->item_list, FALSE, vstp);  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = FreeClickableList (vstp->item_list);
  DeleteMarkedObjects (vstp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);

  RefreshVecScreenList(vstp);
  RedrawVecScreenTool (vstp);  
  ArrowCursor();
  Update();
}

static void RearrangeVecScreenList (GrouP g)
{
  VecScreenToolPtr vstp;

  vstp = (VecScreenToolPtr) GetObjectExtra (g);
  if (vstp == NULL) return;
  
  if (GetValue (g) == 1) {
    vstp->item_list = SortVecScreenList (vstp->item_list);
  } else {
    vstp->item_list = ValNodeSort (vstp->item_list, SortClickableItemLoc);
  }
  PointerToDialog (vstp->clickable_list, vstp->item_list);
}


static void ReportVecScreenDiscrepancies (LogInfoPtr lip, ValNodePtr item_list)
{
  ClickableItemPtr cip;
  if (item_list == NULL || lip == NULL || lip->fp == NULL) return;

  while (item_list != NULL)
  {
    cip = (ClickableItemPtr) item_list->data.ptrvalue;
    if (cip != NULL)
    {
      if (!StringHasNoText (cip->description))
      {
        fprintf (lip->fp, "%s\n", cip->description);
        lip->data_in_log = TRUE;
      }
      ReportVecScreenDiscrepancies (lip, cip->subcategories);
    }
    item_list = item_list->next;
  }
}


static void MakeVecScreenReport (ButtoN b)
{
  VecScreenToolPtr vstp;
  LogInfoPtr       lip;

  vstp = (VecScreenToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  lip = OpenLog ("VecScreen Report");
  ReportVecScreenDiscrepancies (lip, vstp->item_list);
  CloseLog (lip);
  lip = FreeLog (lip);
}


extern void VecScreenTool (IteM i)
{
  BaseFormPtr              bfp;
  VecScreenToolPtr         drfp;
  SeqEntryPtr              sep;
  GrouP                    h, g;
  GrouP                    c, c2, sort_grp;
  ButtoN                   b, b1;
  WindoW                   w;
  OMUserDataPtr            omudp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
      
  drfp = (VecScreenToolPtr) MemNew (sizeof (VecScreenToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "VecScreen Contamination", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupVecScreenTool);
  drfp->form = (ForM) w;
  drfp->formmessage = VecScreenToolMessage;
    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = VecScreenToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = CreateClickableListDialog (h, "", "Location",
                                                    ScrollToDiscrepancyItem, EditDiscrepancyItem, bfp,
                                                    GetDiscrepancyItemText);

  RefreshVecScreenList(drfp);
  
  sort_grp = NormalGroup (h, 2, 0, "Order By", programFont, RearrangeVecScreenList);
  SetObjectExtra (sort_grp, drfp, NULL);
  SetGroupSpacing (sort_grp, 10, 10);
  RadioButton (sort_grp, "Internal, 5', 3'");
  RadioButton (sort_grp, "Accession");
  SetValue (sort_grp, 1);
  
  drfp->trim_options_grp = NormalGroup (h, 3, 0, "Internal Trim Options", programFont, NULL);
  SetGroupSpacing (drfp->trim_options_grp, 10, 10);
  RadioButton (drfp->trim_options_grp, "Trim to closest end");
  RadioButton (drfp->trim_options_grp, "Trim to 5' end");
  RadioButton (drfp->trim_options_grp, "Trim to 3' end");
  SetValue (drfp->trim_options_grp, 1);
  
  c2 = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c2, 10, 10);
  b = PushButton (c2, "Select All", VecScreenToolSelectAllBtn);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c2, "Unselect All", VecScreenToolUnselectAllBtn);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c2, "Unselect Internal", VecScreenToolUnselectInternalBtn);
  SetObjectExtra (b, drfp, NULL);

  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  b1 = PushButton (g, "Remove Selected misc_feats", RemoveSelectedBtn); 
  SetObjectExtra (b1, drfp, NULL);
  MultiLinePrompt (g, "Removes misc_feats with VecScreen hits that you do not want to trim", stdCharWidth * 15, programFont);

  drfp->add_citsub_btn = CheckBox (h, "Add CitSub update to trimmed sequences", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  
  b = PushButton (c, "Make Report", MakeVecScreenReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Trim Selected Regions", TrimSelectedBtn);
  SetObjectExtra (b, drfp, NULL);
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list,
                              (HANDLE) sort_grp,
                              (HANDLE) drfp->trim_options_grp,
                              (HANDLE) c2,
                              (HANDLE) g,
                              (HANDLE) drfp->add_citsub_btn,
                              (HANDLE) c,
                              NULL);

  RealizeWindow (w);
  
  Show (w);
}


/* Barcode Test */
typedef struct barcodeundoitem {
  ValNodePtr object_list;
	struct barcodeundoitem PNTR next;  /* next in linked list */
	struct barcodeundoitem PNTR prev;  /* prev in linked list */
} BarcodeUndoItemData, PNTR BarcodeUndoItemPtr;

typedef ValNodePtr (*BarcodeItemCopy) PROTO ((ValNodePtr));
typedef ValNodePtr (*BarcodeItemFree) PROTO ((ValNodePtr));

static BarcodeUndoItemPtr BarcodeUndoItemNew (ValNodePtr object_list, BarcodeItemCopy copy_func)
{
  BarcodeUndoItemPtr item;
  ValNodePtr         prev = NULL, next_object;

  if (object_list == NULL) return NULL;

  item = (BarcodeUndoItemPtr) MemNew (sizeof (BarcodeUndoItemData));

  item->next = NULL;
  item->prev = NULL;
  while (object_list != NULL)
  {
    if (copy_func == NULL) 
    {
      next_object = ValNodeNew(NULL);
      next_object->choice = object_list->choice;
      next_object->data.ptrvalue = object_list->data.ptrvalue;
    }
    else
    {
      next_object = copy_func (object_list);
    }
    if (next_object != NULL) 
    {
      if (prev == NULL)
      {
        item->object_list = next_object;
      }
      else
      {
        prev->next = next_object;
      }
      prev = next_object;
    }
    object_list = object_list->next;
  }
  return item;
}


static BarcodeUndoItemPtr BarcodeUndoItemFree (BarcodeUndoItemPtr item, BarcodeItemFree free_func)
{
  BarcodeUndoItemPtr item_next;
  ValNodePtr         object_next;

  while (item != NULL) 
  {
    item_next = item->next;
    if (free_func == NULL)
    {
      item->object_list = ValNodeFree (item->object_list);
    }
    else
    {
      while (item->object_list != NULL) 
      {
        object_next = item->object_list->next;
        item->object_list->next = NULL;
        item->object_list = free_func(item->object_list);
        item->object_list = object_next;
      }
    }
    item = MemFree (item);
    item = item_next;
  }
  return item;
}


typedef struct barcodeundolist {
  BarcodeUndoItemPtr list;
  BarcodeUndoItemPtr curr;
  
  BarcodeItemCopy    copy_func;
  BarcodeItemFree    free_func;
} BarcodeUndoListData, PNTR BarcodeUndoListPtr;


static BarcodeUndoListPtr BarcodeUndoListFree (BarcodeUndoListPtr undo_list)
{

  if (undo_list != NULL)
  {
    undo_list->list = BarcodeUndoItemFree (undo_list->list, undo_list->free_func);
    undo_list = MemFree (undo_list);
  }
  return undo_list;
}


static BarcodeUndoListPtr BarcodeUndoListNew (BarcodeItemCopy copy_func, BarcodeItemFree free_func)
{
  BarcodeUndoListPtr undo_list;

  undo_list = (BarcodeUndoListPtr) MemNew (sizeof (BarcodeUndoListData));
  undo_list->list = NULL;
  undo_list->curr = NULL;
  undo_list->copy_func = copy_func;
  undo_list->free_func = free_func;
  return undo_list;
}


static void AddToBarcodeUndoList (BarcodeUndoListPtr undo_list, ValNodePtr object_list)
{
  BarcodeUndoItemPtr item_new;

  if (undo_list == NULL || object_list == NULL) return;

  item_new = BarcodeUndoItemNew (object_list, undo_list->copy_func);
  if (undo_list->curr == NULL) 
  {
    /* if curr is NULL, there is nothing to undo.  list could still contain a list of 
     * actions that had been undone.
     */
    undo_list->list = BarcodeUndoItemFree (undo_list->list, undo_list->free_func);
    undo_list->list = item_new;
  }
  else
  {
    undo_list->curr->next = BarcodeUndoItemFree (undo_list->curr->next, undo_list->free_func);
    undo_list->curr->next = item_new;
    item_new->prev = undo_list->curr;
  }
  undo_list->curr = item_new;
}


static ValNodePtr GetUndoFromBarcodeUndoList (BarcodeUndoListPtr undo_list)
{
  ValNodePtr object_list = NULL;
  if (undo_list == NULL || undo_list->curr == NULL)
  {
    return NULL;
  } 
  else
  {
    object_list = undo_list->curr->object_list;
    undo_list->curr = undo_list->curr->prev;
  }
  return object_list;
}


static ValNodePtr GetRedoFromBarcodeUndoList (BarcodeUndoListPtr undo_list)
{
  ValNodePtr object_list = NULL;

  if (undo_list == NULL)
  {
    return NULL;
  } 
  else 
  {
    if (undo_list->curr == NULL) 
    {
      undo_list->curr = undo_list->list;
    }
    else if (undo_list->curr->next == NULL)
    {
      return NULL;
    } 
    else
    {
      undo_list->curr = undo_list->curr->next;
    }
    object_list = undo_list->curr->object_list;
  }
  return object_list;
}


static Boolean IsUndoAvailable (BarcodeUndoListPtr undo_list)
{
  if (undo_list == NULL || undo_list->curr == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static Boolean IsRedoAvailable (BarcodeUndoListPtr undo_list)
{
  if (undo_list == NULL || undo_list->list == NULL)
  {
    return FALSE;
  }
  
  if (undo_list->curr == NULL || undo_list->curr->next != NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static ValNodePtr BarcodeResultsVNCopyFunc (ValNodePtr vnp)
{
  BarcodeTestResultsPtr res1, res2;
  ValNodePtr            vnp_new;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  res1 = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
  res2 = BarcodeTestResultsCopy (res1);
  vnp_new = ValNodeNew(NULL);
  vnp_new->data.ptrvalue = res2;
  vnp_new->choice = vnp->choice;
  return vnp_new;
}


static ValNodePtr BarcodeResultsVNFreeFunc (ValNodePtr vnp)
{
  BarcodeTestResultsPtr res;

  if (vnp != NULL)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    res = BarcodeTestResultsFree (res);
    vnp->next = BarcodeResultsVNFreeFunc (vnp->next);
    vnp = ValNodeFree (vnp);
  }
  return vnp;
}


/* Barcode Tool */
struct barcodetool;

typedef void (*RefreshBarcodeDataFunc) PROTO ((Pointer));

#define BARCODE_TOOL_BLOCK \
  FORM_MESSAGE_BLOCK \
  DialoG          clickable_list; \
  BaseFormPtr     bfp; \
  PrompT          pass_fail_summary; \
  ButtoN          undo; \
  ButtoN          undo_all; \
  ButtoN          redo; \
  ValNodePtr      item_list; \
  SeqEntryPtr     top_sep; \
  LogInfoPtr      lip; \
  BarcodeTestConfigPtr cfg; \
  BarcodeUndoListPtr   undo_list; \
  RefreshBarcodeDataFunc refresh_func;

typedef struct barcodetool {
BARCODE_TOOL_BLOCK
} BarcodeToolData, PNTR BarcodeToolPtr;

static void CleanupBarcodeTool (GraphiC g, VoidPtr data)

{
  BarcodeToolPtr drfp;

  drfp = (BarcodeToolPtr) data;
  if (drfp != NULL) {
    drfp->item_list = BarcodeTestResultsListFree (drfp->item_list);
    drfp->undo_list = BarcodeUndoListFree (drfp->undo_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}

static void BarcodeToolMessage (ForM f, Int2 mssg)

{
  BarcodeToolPtr drfp;

  drfp = (BarcodeToolPtr) GetObjectExtra (f);
  if (drfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_COPY :
/*        CopyDiscrepancyReportToClipboard (drfp); */
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        drfp->item_list = FreeClickableList (drfp->item_list);
        PointerToDialog (drfp->clickable_list, NULL);
        break;
      default :
        if (drfp->appmessage != NULL) {
          drfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Int2 LIBCALLBACK BarcodeToolMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  BarcodeToolPtr drfp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  drfp = (BarcodeToolPtr) omudp->userdata.ptrvalue;
  if (drfp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (drfp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (drfp->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (drfp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (drfp->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}


static void RefreshBarcodeList (Pointer data)
{
  BarcodeToolPtr vstp;
  ValNodePtr pass_fail_list = NULL, vnp;
  Int4       num_fail = 0, num_pass = 0;
  Char       msg[30];

  vstp = (BarcodeToolPtr) data;
  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = BarcodeTestResultsListFree (vstp->item_list);

  vstp->item_list = GetBarcodeFailedAccessionList(vstp->top_sep, vstp->cfg);
  PointerToDialog (vstp->clickable_list, vstp->item_list);  

  pass_fail_list = GetBarcodePassFail (vstp->top_sep, vstp->cfg);
  for (vnp = pass_fail_list; vnp != NULL; vnp = vnp->next) 
  {
    if (PassBarcodeTests(vnp->data.ptrvalue)) 
    {
      num_pass ++;
    }
    else
    {
      num_fail++;
    }
  }
  pass_fail_list = BarcodeTestResultsListFree (pass_fail_list);
  sprintf (msg, "%d pass, %d fail", num_pass, num_fail);
  SetTitle (vstp->pass_fail_summary, msg);
}


static void RedrawBarcodeTool (BarcodeToolPtr vstp)
{
  RecT     r;

  ObjectRect (vstp->clickable_list, &r);
  InsetRect (&r, -1, -1);
  InvalRect (&r);
}


static void EnableUndoRedo (BarcodeToolPtr vstp)
{
  if (vstp == NULL) return;

  if (IsUndoAvailable (vstp->undo_list)) 
  {
    Enable (vstp->undo);
    Enable (vstp->undo_all);
  } 
  else 
  {
    Disable (vstp->undo);
    Disable (vstp->undo_all);
  }

  if (IsRedoAvailable (vstp->undo_list)) 
  {
    Enable (vstp->redo);
  } 
  else 
  {
    Disable (vstp->redo);
  }

}


static void BarcodeUndoButton (ButtoN b)
{
  BarcodeToolPtr vstp;
  LogInfoPtr     lip;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = GetUndoFromBarcodeUndoList (vstp->undo_list);
  EnableUndoRedo (vstp);

  lip = OpenLog ("BARCODE Keywords Put Back");

  ApplyBarcodeKeywords (lip->fp, object_list);

  RefreshBarcodeList (vstp);
  RedrawBarcodeTool (vstp);  

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
}


static void BarcodeUndoAllButton (ButtoN b)
{
  BarcodeToolPtr vstp;
  LogInfoPtr     lip;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  lip = OpenLog ("BARCODE Keywords Put Back");
  while (IsUndoAvailable (vstp->undo_list))
  {
    object_list = GetUndoFromBarcodeUndoList (vstp->undo_list);
    ApplyBarcodeKeywords (lip->fp, object_list);
  }

  EnableUndoRedo (vstp);
  RefreshBarcodeList (vstp);
  RedrawBarcodeTool (vstp);  

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
}




static void BarcodeRedoButton (ButtoN b)
{
  BarcodeToolPtr vstp;
  LogInfoPtr     lip;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = GetRedoFromBarcodeUndoList (vstp->undo_list);
  EnableUndoRedo (vstp);

  lip = OpenLog ("BARCODE Keywords Removed");

  RemoveBarcodeKeywords (lip->fp, object_list);

  RefreshBarcodeList (vstp);
  RedrawBarcodeTool (vstp);  

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);  
}


static void RemoveSelectedKeywordsBtn (ButtoN b)
{
  BarcodeToolPtr vstp;
  LogInfoPtr     lip;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = DialogToPointer (vstp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any Bioseqs - remove BARCODE keyword from all?"))
    {
      object_list = vstp->item_list;
    }
    else
    {
      return;
    }
  }
  
  lip = OpenLog ("BARCODE Keywords Removed");
  RemoveBarcodeKeywords (lip->fp, object_list);

  /* put in undo queue */
  AddToBarcodeUndoList (vstp->undo_list, object_list);

  EnableUndoRedo (vstp);  

  if (object_list != vstp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  RefreshBarcodeList (vstp);
  RedrawBarcodeTool (vstp);  

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);


}


extern void BarcodeRefreshButton (ButtoN b)
{
  BarcodeToolPtr vstp;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  if (vstp->refresh_func != NULL) 
  {
    (vstp->refresh_func) (vstp);
  }
  RedrawBarcodeTool (vstp);  
}


extern void BarcodeReportButton (ButtoN b)
{
  BarcodeToolPtr vstp;
  LogInfoPtr     lip;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = (ValNodePtr) DialogToPointer (vstp->clickable_list);
  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any Bioseqs - include all Bioseqs in report?"))
    {
      object_list = vstp->item_list;
    }
    else
    {
      return;
    }
  }

  lip = OpenLog ("BARCODE Discrepancies");
  WriteBarcodeDiscrepancies (lip->fp, object_list);
  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
  if (object_list != vstp->item_list)
  {
    object_list = ValNodeFree (object_list);
  }
}

extern void BarcodeTestComplianceReport (ButtoN b)
{
  BarcodeToolPtr vstp;
  ValNodePtr     object_list;
  LogInfoPtr     lip;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = GetBarcodePassFail (vstp->top_sep, vstp->cfg);
  lip = OpenLog ("BARCODE Compliance");
  WriteBarcodeTestCompliance (lip->fp, object_list);
  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
  object_list = BarcodeTestResultsListFree (object_list);
}


extern DialoG BarcodeTestResultsDisplay (GrouP h, BarcodeTestConfigPtr cfg);

static void BarcodeTestMakeTagTable (ButtoN b)
{
  BarcodeToolPtr vstp;
  ValNodePtr     object_list;
  LogInfoPtr     lip;
  WindoW         w;
  GrouP          h, g1, c;
  ButtoN         b2;
  TexT                  acc_list_txt;
  ModalAcceptCancelData acd;
  CharPtr               acc_list;
  ValNodePtr            id_list = NULL, vnp;
  CharPtr               genbank_id, barcode_id;
  Char                  id_str[128];
  BioseqPtr             bsp;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;
    
  w = MovableModalWindow (-20, -13, -10, -10, "Create BARCODE Tag Table", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g1 = NormalGroup (h, 3, 0, "Accessions", programFont, NULL);
  StaticPrompt (g1, "List Accessions", 0, dialogTextHeight, programFont, 'l');
  acc_list_txt = DialogText (g1, "", 30, NULL); 
  StaticPrompt (g1, "(Leave blank to use all)", 0, dialogTextHeight, programFont, 'l');
    
  c = HiddenGroup (h, 2, 0, NULL);
  b2 = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b2, &acd, NULL);
  b2 = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b2, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1,
                              (HANDLE) c, (HANDLE) NULL);

  Show (w);
  Select (w);

  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  
  if (!acd.cancelled)
  {
    lip = OpenLog ("BARCODE Tag Table");
    acc_list = SaveStringFromText (acc_list_txt);
    if (StringHasNoText (acc_list)) 
    {
      object_list = GetBarcodePassFail (vstp->top_sep, vstp->cfg);
      WriteBarcodeTagTable (lip->fp, object_list);
      object_list = BarcodeTestResultsListFree (object_list);
    }
    else
    {
      id_list = ParseAccessionNumberListFromString (acc_list, vstp->top_sep);
      for (vnp = id_list; vnp != NULL; vnp = vnp->next)
      {
        bsp = BioseqFind (vnp->data.ptrvalue);
        if (bsp == NULL)
        {
          SeqIdWrite (vnp->data.ptrvalue, id_str, PRINTID_REPORT, sizeof (id_str));
          fprintf (lip->fp, "%s\tNOT IN RECORD\t\n", id_str);
        }
        else
        {
          barcode_id = BarcodeTestBarcodeIdString (bsp);
          genbank_id = BarcodeTestGenbankIdString (bsp);
          fprintf (lip->fp, "%s\t%s\t\n", genbank_id, barcode_id);
          barcode_id = MemFree (barcode_id);
          genbank_id = MemFree (genbank_id);
        }
      }
      id_list = FreeSeqIdList (id_list);
    }      
    acc_list = MemFree (acc_list);

    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
  }
  Remove (w);
}


typedef struct tagtable {
  BioseqPtr bsp;
  CharPtr   accession;
  CharPtr   old_tag;
  CharPtr   new_tag;
} TagTableData, PNTR TagTablePtr;

static TagTablePtr TagTableFree (TagTablePtr ttp)
{
  if (ttp != NULL) 
  {
    ttp->accession = MemFree (ttp->accession);
    ttp->old_tag = MemFree (ttp->old_tag);
    ttp->new_tag = MemFree (ttp->new_tag);
    ttp = MemFree (ttp);
  }
  return ttp;
}

static ValNodePtr TagTableListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = TagTableFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static CharPtr SkipPastGnlUoguelph (CharPtr token)
{
  Int4        uoguelph_len = StringLen ("uoguelph|");

  if (token == NULL) return NULL;

  if (StringNICmp (token, "gnl|", 4) == 0)
  {
    token += 4;
  }
  if (StringNICmp (token, "uoguelph|", uoguelph_len) == 0)
  {
    token += uoguelph_len;
  }
  return token;
}


static TagTablePtr ReadReplaceTagTableLine (CharPtr line, SeqEntryPtr sep, LogInfoPtr lip)
{
  CharPtr  token, cp, delimiters = " \t,;";
  Int4     len;
  Char     ch;
  SeqIdPtr  sip;
  BioseqPtr bsp = NULL;
  CharPtr   accession = NULL, old_tag = NULL, new_tag = NULL;
  TagTablePtr ttp = NULL;
  Boolean     error = FALSE;
  Char        id_num[20];
  DbtagPtr    dbt;

  if (StringHasNoText (line) || lip == NULL || lip->fp == NULL) return NULL;

  token = line;
  while (*token != 0 && !error) 
  {
    len = StringCSpn (token, delimiters);
    if (len == 0)
    {
      token += StringSpn (token, delimiters);
    }
    else
    {
      cp = token + len;
      ch = *cp;
      *cp = 0;
      if (bsp == NULL) 
      {
        accession = StringSave (token);
        sip = CreateSeqIdFromText (token, sep);
        bsp = BioseqFind (sip);
        sip = SeqIdFree (sip);
        if (bsp == NULL)
        {
          fprintf (lip->fp, "No Bioseq found for %s\n", token);
          lip->data_in_log = TRUE;
          error = TRUE;
        }
      }
      else if (old_tag == NULL)
      {
        /* check for old tag matching what is found on sequence */
        sip = bsp->id;
        while (sip != NULL && !IsBarcodeID(sip))
        {
          sip = sip->next;
        }
        if (sip == NULL || (dbt = (DbtagPtr) sip->data.ptrvalue) == NULL || dbt->tag == NULL)
        {
          fprintf (lip->fp, "No BARCODE ID for %s\n", accession);
          error = TRUE;
        }
        else
        {
          token = SkipPastGnlUoguelph(token);
          if (dbt->tag->id > 0) 
          {
            sprintf (id_num, "%d", dbt->tag->id);
            if (StringCmp (token, id_num) != 0) 
            {
              fprintf (lip->fp, "Old Tag is incorrect for %s\n", accession);
              lip->data_in_log = TRUE;
              error = TRUE;
            }
          } 
          else if (StringCmp (token, dbt->tag->str) != 0)
          {
            fprintf (lip->fp, "Old Tag is incorrect for %s\n", accession);
            lip->data_in_log = TRUE;
            error = TRUE;
          }
          if (!error)
          {
            old_tag = StringSave (token);
          }
        }
      }
      else if (new_tag == NULL)
      {
        token = SkipPastGnlUoguelph (token);
        new_tag = StringSave (token);
      }
      else
      {
        fprintf (lip->fp, "Warning!  Too many columns in table!\n");
        lip->data_in_log = TRUE;
      }
      *cp = ch;
      token = cp;
    }
  }

  if (!error) 
  {
    if (old_tag == NULL)
    {
      fprintf (lip->fp, "No old tag for %s!\n", accession);
      lip->data_in_log = TRUE;
    }
    else if (new_tag == NULL)
    {
      fprintf (lip->fp, "No new tag for %s\n", accession);
      lip->data_in_log = TRUE;
    }
    else
    {    
      ttp = (TagTablePtr) MemNew (sizeof (TagTableData));
      ttp->bsp = bsp;
      ttp->accession = accession;
      accession = NULL;
      ttp->old_tag = old_tag;
      old_tag = NULL;
      ttp->new_tag = new_tag;
      new_tag = NULL;
    }
  }
  accession = MemFree (accession);
  old_tag = MemFree (old_tag);
  new_tag = MemFree (new_tag);
  return ttp;
}


static TagTablePtr ReadNewTagTableLine (CharPtr line, SeqEntryPtr sep, LogInfoPtr lip)
{
  CharPtr  token, cp, delimiters = " \t,;";
  Int4     len;
  Char     ch;
  SeqIdPtr  sip;
  BioseqPtr bsp = NULL;
  CharPtr   accession = NULL, new_tag = NULL;
  TagTablePtr ttp = NULL;
  Boolean     error = FALSE;

  if (StringHasNoText (line) || lip == NULL || lip->fp == NULL) return NULL;

  token = line;
  while (*token != 0 && !error) 
  {
    len = StringCSpn (token, delimiters);
    if (len == 0)
    {
      token += StringSpn (token, delimiters);
    }
    else
    {
      cp = token + len;
      ch = *cp;
      *cp = 0;
      if (bsp == NULL) 
      {
        accession = StringSave (token);
        sip = CreateSeqIdFromText (token, sep);
        bsp = BioseqFind (sip);
        sip = SeqIdFree (sip);
        if (bsp == NULL)
        {
          fprintf (lip->fp, "No Bioseq found for %s\n", token);
          lip->data_in_log = TRUE;
          error = TRUE;
        }
        else
        {
          /* check for old tag - there shouldn't be one */
          sip = bsp->id;
          while (sip != NULL && !IsBarcodeID(sip))
          {
            sip = sip->next;
          }
          if (sip != NULL)
          {
            fprintf (lip->fp, "%s already has BARCODE ID\n", accession);
            error = TRUE;
          }
        }
      }
      else if (new_tag == NULL)
      {
        token = SkipPastGnlUoguelph (token);
        new_tag = StringSave (token);
      }
      else
      {
        fprintf (lip->fp, "Warning!  Too many columns in table!\n");
        lip->data_in_log = TRUE;
      }
      *cp = ch;
      token = cp;
    }
  }

  if (!error) 
  {
    if (new_tag == NULL)
    {
      fprintf (lip->fp, "No new tag for %s\n", accession);
      lip->data_in_log = TRUE;
    }
    else
    {    
      ttp = (TagTablePtr) MemNew (sizeof (TagTableData));
      ttp->bsp = bsp;
      ttp->accession = accession;
      accession = NULL;
      ttp->old_tag = NULL;
      ttp->new_tag = new_tag;
      new_tag = NULL;
    }
  }
  accession = MemFree (accession);
  new_tag = MemFree (new_tag);
  return ttp;
}


static ValNodePtr ReadTagTable (SeqEntryPtr sep, LogInfoPtr lip, Boolean replace)
{
  Char           path [PATH_MAX];
  ValNodePtr     line_list = NULL;
  ReadBufferData rbd;
  CharPtr        line;
  TagTablePtr    ttp;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return NULL;
  
  rbd.fp = FileOpen (path, "r");
  if (rbd.fp == NULL) return NULL;
  rbd.current_data = NULL;

  line = AbstractReadFunction (&rbd);
  while (line != NULL) 
  {
    if (replace) 
    {
      ttp = ReadReplaceTagTableLine (line, sep, lip);
    }
    else
    {
      ttp = ReadNewTagTableLine (line, sep, lip);
    }
    if (ttp != NULL) 
    {
      ValNodeAddPointer (&line_list, 0, ttp);
    }
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }

  FileClose (rbd.fp);
  return line_list;
}


static void ApplyTagTable (BarcodeToolPtr vstp, Boolean replace)
{
  TagTablePtr   ttp;
  ValNodePtr    tag_list, vnp;
  LogInfoPtr     lip;
  Boolean        errors;
  SeqIdPtr       sip;
  DbtagPtr       dbt;
  
  if (vstp == NULL) return;
  
  lip = OpenLog ("Tag Table Errors");

  tag_list = ReadTagTable (vstp->top_sep, lip, replace);
  errors = lip->data_in_log;
  CloseLog (lip);
  lip = FreeLog (lip);

  if (errors) 
  {
    if (ANS_NO == Message (MSG_YN, "Errors in table.  Continue with matched rows only?"))
    {
      tag_list = TagTableListFree (tag_list);
      return;
    }
  }

  /* NOW DO SOMETHING USEFUL */
  for (vnp = tag_list; vnp != NULL; vnp = vnp->next)
  {
    ttp = (TagTablePtr) vnp->data.ptrvalue;
    if (replace)
    {
      sip = ttp->bsp->id;
      while (sip != NULL && !IsBarcodeID(sip))
      {
        sip = sip->next;
      }
    } 
    else
    {
      sip = ValNodeNew (ttp->bsp->id);
      if (ttp->bsp->id == NULL) 
      {
        ttp->bsp->id = sip;
      }
      sip->choice = SEQID_GENERAL;
      dbt = DbtagNew ();
      sip->data.ptrvalue = dbt;
      dbt->db = StringSave ("uoguelph");
      dbt->tag = ObjectIdNew ();
    }
    if (sip != NULL) 
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      dbt->tag->str = MemFree (dbt->tag->str);
      dbt->tag->id = 0;
      dbt->tag->str = StringSave (ttp->new_tag);
      SeqMgrReplaceInBioseqIndex (ttp->bsp);
    }
  }

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);

  tag_list = TagTableListFree (tag_list);

  RefreshBarcodeList(vstp);
}


static void BarcodeTestImportTagTable (ButtoN b)
{
  ApplyTagTable (GetObjectExtra (b), TRUE);
}


static void BarcodeTestApplyTagTable (ButtoN b)
{
  ApplyTagTable (GetObjectExtra (b), FALSE);
}


static void BarcodeComprehensiveReportButton (ButtoN b)
{
  BarcodeToolPtr         drfp;
  LogInfoPtr             lip;
  ValNodePtr             object_list;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  lip = OpenLog ("BARCODE Discrepancies");
  object_list = GetBarcodePassFail (drfp->top_sep, drfp->cfg);
  WriteBarcodeTestComprehensive (lip->fp, object_list);
  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
  object_list = BarcodeTestResultsListFree (object_list);
}


static void ReportPolymorphismCallback (BioseqPtr bsp, Pointer data)
{
  Int4 num_p;
  LogInfoPtr lip;
  Char       id_txt[100];

  if (bsp == NULL || data == NULL || !HasBARCODEKeyword(bsp)) {
    return;
  }

  lip = (LogInfoPtr) data;
  num_p = CountPolymorphismsInBioseq (bsp);
  if (num_p > 0) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    fprintf (lip->fp, "%s: %d polymophisms\n", id_txt, num_p);
    lip->data_in_log = TRUE;
  }

}


static void BarcodeReportPolymorphism (ButtoN b)
{
  BarcodeToolPtr         drfp;
  LogInfoPtr             lip;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  lip = OpenLog ("Barcode Sequence Polymorphism");
  VisitBioseqsInSep (GetTopSeqEntryForEntityID (drfp->input_entityID), lip, ReportPolymorphismCallback);
  CloseLog (lip);
  if (!lip->data_in_log) {
    Message (MSG_OK, "No polymorphism found");
  }
  lip = FreeLog (lip);
}


extern void BarcodeTestTool (IteM i)
{
  BaseFormPtr              bfp;
  BarcodeToolPtr         drfp;
  SeqEntryPtr              sep;
  GrouP                    h;
  GrouP                    c, c3;
  ButtoN                   b;
  WindoW                   w;
  OMUserDataPtr            omudp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
      
  drfp = (BarcodeToolPtr) MemNew (sizeof (BarcodeToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Barcode Test", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupBarcodeTool);
  drfp->form = (ForM) w;
  drfp->formmessage = BarcodeToolMessage;
  drfp->refresh_func = RefreshBarcodeList;
    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = BarcodeToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  drfp->cfg = NULL;
  drfp->undo_list = BarcodeUndoListNew (BarcodeResultsVNCopyFunc, BarcodeResultsVNFreeFunc);
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = BarcodeTestResultsDisplay (h, drfp->cfg);

  drfp->pass_fail_summary = StaticPrompt (h, "0 Pass, 0 Fail", 20 * stdCharWidth, dialogTextHeight, programFont, 'l');
  RefreshBarcodeList(drfp);

  c3 = HiddenGroup (h, 9, 0, NULL);
  SetGroupSpacing (c3, 10, 10);
  b = PushButton (c3, "Compliance Report", BarcodeTestComplianceReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c3, "Failure Report", BarcodeReportButton);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c3, "Comprehensive Report", BarcodeComprehensiveReportButton);
  SetObjectExtra (b, drfp, NULL);

  b = PushButton (c3, "Report Polymorphism", BarcodeReportPolymorphism);
  SetObjectExtra (b, drfp, NULL);

  b = PushButton (c3, "Replace Tags", BarcodeTestImportTagTable);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c3, "Add New Tags", BarcodeTestApplyTagTable);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c3, "Make Tag Table", BarcodeTestMakeTagTable);
  SetObjectExtra (b, drfp, NULL);

  b = PushButton (c3, "Refresh List", BarcodeRefreshButton);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c3, "Remove BARCODE Keyword from Selected Sequences", RemoveSelectedKeywordsBtn); 
  SetObjectExtra (b, drfp, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
    
  drfp->undo = PushButton (c, "Undo", BarcodeUndoButton);
  SetObjectExtra (drfp->undo, drfp, NULL);

  drfp->undo_all = PushButton (c, "Undo All", BarcodeUndoAllButton);
  SetObjectExtra (drfp->undo_all, drfp, NULL);
  
  drfp->redo = PushButton (c, "Redo", BarcodeRedoButton);
  SetObjectExtra (drfp->redo, drfp, NULL);

  EnableUndoRedo (drfp);

  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) drfp->pass_fail_summary, (HANDLE) c3, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);
}


typedef struct replaceid {
  SeqIdPtr orig_id;
  SeqIdPtr new_id;
} ReplaceIdData, PNTR ReplaceIdPtr;

/********************************************************************
*
* SeqLocReplaceID
*   replaces the Seq-Id in a Seq-Loc (slp) with a new Seq-Id (new_sip)
*
**********************************************************************/
static SeqLocPtr ReplaceSeqLocID (SeqLocPtr slp, SeqIdPtr orig_sip, SeqIdPtr new_sip)
{
  SeqLocPtr        curr;
  PackSeqPntPtr    pspp;
  SeqIntPtr        target_sit;
  SeqBondPtr       sbp;
  SeqPntPtr        spp;
  
  if (slp == NULL || orig_sip == NULL || new_sip == NULL) return slp;

  switch (slp->choice) {
     case SEQLOC_PACKED_INT :
     case SEQLOC_MIX :
     case SEQLOC_EQUIV :
        curr = NULL;
        while ((curr = SeqLocFindNext (slp, curr)) != NULL) {
           curr = ReplaceSeqLocID (curr, orig_sip, new_sip);
        }
        break;
     case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          if (SeqIdComp (orig_sip, pspp->id) == SIC_YES) {
            SeqIdFree (pspp->id);
            pspp->id = SeqIdDup (new_sip);
          }
        }
        break;
     case SEQLOC_EMPTY :
     case SEQLOC_WHOLE :
        if (SeqIdComp (orig_sip, slp->data.ptrvalue) == SIC_YES) {
          SeqIdFree ((SeqIdPtr) slp->data.ptrvalue);
          slp->data.ptrvalue = (Pointer) SeqIdDup (new_sip);
        }
        break;
     case SEQLOC_INT :
        target_sit = (SeqIntPtr) slp->data.ptrvalue;
        if (SeqIdComp (orig_sip, target_sit->id) == SIC_YES) {
          SeqIdFree (target_sit->id);
          target_sit->id = SeqIdDup (new_sip);
        }
        break;
     case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (SeqIdComp (orig_sip, spp->id) == SIC_YES) {
          SeqIdFree(spp->id);
          spp->id = SeqIdDup(new_sip);
        }
        break;
     case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp == NULL || sbp->a == NULL || sbp->b == NULL) break;
        if (SeqIdComp (orig_sip, sbp->a->id) == SIC_YES) {
          spp = sbp->a;
          SeqIdFree(spp->id);
          spp->id = SeqIdDup(new_sip);
        }
        if (SeqIdComp (orig_sip, sbp->b->id) == SIC_YES) {
          spp = sbp->b;
          SeqIdFree(spp->id);
          spp->id = SeqIdDup(new_sip);
        }
        break;
     default :
        break;
  }
  return slp;
}


static void ReplaceCodingRegionLocationIds (CdRegionPtr crp, SeqIdPtr orig_sip, SeqIdPtr new_sip)
{
  CodeBreakPtr cbp;
  
  if (crp == NULL) return;
  cbp = crp->code_break;
  while (cbp != NULL) {
    cbp->loc = ReplaceSeqLocID (cbp->loc, orig_sip, new_sip);
    cbp = cbp->next;
  }
}

static void ReplaceRNALocationIds (RnaRefPtr rrp, SeqIdPtr orig_sip, SeqIdPtr new_sip)
{
  tRNAPtr trp;
  
  if (rrp != NULL && rrp->ext.choice == 2) {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL) {
      trp->anticodon = ReplaceSeqLocID (trp->anticodon, orig_sip, new_sip);
    }
  }
}
  


static void ReplaceFeatureLocationIds (SeqFeatPtr sfp, Pointer userdata)
{
  ReplaceIdPtr rip;
  
  if (sfp == NULL || userdata == NULL) return;
  rip = (ReplaceIdPtr) userdata;
  if (rip->orig_id == NULL || rip->new_id == NULL) return;
  
  sfp->location = ReplaceSeqLocID (sfp->location, rip->orig_id, rip->new_id);
  sfp->product = ReplaceSeqLocID (sfp->product, rip->orig_id, rip->new_id); 
  
  if (sfp->data.value.ptrvalue == NULL) return;
  
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    ReplaceCodingRegionLocationIds (sfp->data.value.ptrvalue, rip->orig_id, rip->new_id);    
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    ReplaceRNALocationIds (sfp->data.value.ptrvalue, rip->orig_id, rip->new_id);
  }  
}


typedef struct proteinidreplacement {
  BioseqPtr bsp;
  SeqIdPtr  new_id;
} ProteinIdReplacementData, PNTR ProteinIdReplacementPtr;


typedef struct proteinidreplacementlist {
  ValNodePtr id_list;
  CharPtr    db;
} ProteinIdReplacementListData, PNTR ProteinIdReplacementListPtr;

static CharPtr GetLocusTagForProteinBioseq (BioseqPtr bsp)
{
  GeneRefPtr grp = NULL;
  SeqFeatPtr cds, gene;
  
  if (bsp == NULL || !ISA_aa (bsp->mol)) {
    return NULL;
  }
  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  if (cds == NULL) {
    return NULL;
  }
  grp = SeqMgrGetGeneXref (cds);
  if (grp == NULL) {
    gene = SeqMgrGetOverlappingGene (cds->location, NULL);
    if (gene != NULL) {
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
    }
  }

  if (grp == NULL) {
    return NULL;
  } else {
    return grp->locus_tag;
  }
}

static DbtagPtr GetDbtagId (BioseqPtr bsp) 
{
  SeqIdPtr sip;
  DbtagPtr dbt = NULL;
  
  if (bsp == NULL) return NULL;
  sip = bsp->id;
  for (sip = bsp->id; sip != NULL && dbt == NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENERAL)
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt != NULL && StringICmp (dbt->db, "TMSMART") == 0)
      {
        dbt = NULL;
      }
    }
  }
  return dbt;
}


static void BuildProteinIDList (BioseqPtr bsp, Pointer userdata)
{
  ProteinIdReplacementListPtr pid_list;
  ProteinIdReplacementPtr     pid;
  CharPtr                     locus_tag;
  SeqIdPtr                    sip;
  DbtagPtr                    dbt;
  
  if (bsp == NULL || userdata == NULL ||!ISA_aa (bsp->mol)) {
    return;
  }
  
  dbt = GetDbtagId (bsp);
  if (dbt != NULL)
  {
    /* already has general ID */
    return;
  }
  
  pid_list = (ProteinIdReplacementListPtr) userdata;
  if (StringHasNoText (pid_list->db)) {
    return;
  }
  
  locus_tag = GetLocusTagForProteinBioseq (bsp);
  
  pid = (ProteinIdReplacementPtr) MemNew (sizeof (ProteinIdReplacementData));
  pid->bsp = bsp;
  pid->new_id = NULL;
  if (!StringHasNoText (locus_tag)) {
    dbt = DbtagNew ();
    dbt->db = StringSave (pid_list->db);
    dbt->tag = ObjectIdNew ();
    dbt->tag->id = 0;
    dbt->tag->str = StringSave (locus_tag);
    sip = ValNodeNew (NULL);
    sip->choice = SEQID_GENERAL;
    sip->data.ptrvalue = dbt;
    pid->new_id = sip;
  }
  ValNodeAddPointer (&(pid_list->id_list), 0, pid);
}

static void FindLocusTagIdDb (BioseqPtr bsp, Pointer userdata)
{
  CharPtr PNTR db;
  DbtagPtr     dbt;
  
  if (bsp == NULL || !ISA_aa (bsp->mol) || userdata == NULL) return;
  db = (CharPtr PNTR) userdata;
  if (*db != NULL) return;
  
  dbt = GetDbtagId (bsp);
  if (dbt == NULL || StringHasNoText (dbt->db)) return;
  
  *db = dbt->db;
}


extern void MakeGeneralIDsFromLocusTags (IteM i)
{
  BaseFormPtr                  bfp;
  SeqEntryPtr                  sep;
  CharPtr                      db = NULL;
  ProteinIdReplacementListData pird;
  ProteinIdReplacementPtr      pid;
  ValNodePtr                   duplicated_ids = NULL;
  ValNodePtr                   missing_ids = NULL;
  ValNodePtr                   vnp;
  BioseqPtr                    bsp;
  ClickableItemPtr             cip;
  ReplaceIdData                rid;
  ValNodePtr                   clickable_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, &db, FindLocusTagIdDb);
  if (StringHasNoText (db)) {
    Message (MSG_ERROR, "Unable to locate database tag!");
    return;
  }

  pird.db = db;
  pird.id_list = NULL;
  
  VisitBioseqsInSep (sep, &pird, BuildProteinIDList);
  if (pird.id_list == NULL) {
    Message (MSG_OK, "No IDs to replace!");
    return;
  }
  
  /* build list of errors */
  for (vnp = pird.id_list; vnp != NULL; vnp = vnp->next) {
    pid = (ProteinIdReplacementPtr) vnp->data.ptrvalue;
    if (pid->new_id == NULL) {
      ValNodeAddPointer (&missing_ids, OBJ_BIOSEQ, pid->bsp);
    } else {      
      bsp = BioseqFind (pid->new_id);
      if (bsp != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (cip, 0, sizeof (ClickableItemData));
        cip->description = StringSave ("Conflicting locus_tag ID");        
        ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
        ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, pid->bsp);
        ValNodeAddPointer (&duplicated_ids, 0, cip);
      }
    }
  }
  if (duplicated_ids == NULL && missing_ids == NULL) {
    for (vnp = pird.id_list; vnp != NULL; vnp = vnp->next) {
      pid = (ProteinIdReplacementPtr) vnp->data.ptrvalue;
      rid.new_id = pid->new_id;
      rid.orig_id = pid->bsp->id;
      /* replace ID in feature locations */
      VisitFeaturesInSep (sep, &rid, ReplaceFeatureLocationIds);
      /* add ID to sequence ID list */ 
      ValNodeLink (&(pid->bsp->id), rid.new_id);
      SeqMgrReplaceInBioseqIndex (pid->bsp);
      /* add to list for display */
      ValNodeAddPointer (&clickable_list, OBJ_BIOSEQ, pid->bsp);      
    }
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    /* Show what was added */
    cip = NewClickableItem (0, "%d IDs were added", clickable_list);
    clickable_list = NULL;
    ValNodeAddPointer (&clickable_list, 0, cip);
    ShowClickableItemList (clickable_list, bfp, "New IDs", "", "Sequences");
  } else {
    if (missing_ids != NULL) {
      cip = NewClickableItem (0, "%d sequences do not have locus_tags", missing_ids);
      ValNodeAddPointer (&clickable_list, 0, cip);
    }
    if (duplicated_ids != NULL) {
      cip = NewClickableItem (0, "%d sequences have duplicated locus_tags", NULL);
      cip->subcategories = duplicated_ids;
      ValNodeAddPointer (&clickable_list, 0, cip);
    }
    Message (MSG_ERROR, "Problems generating IDs");
    ShowClickableItemList (clickable_list, bfp, "Unable to automatically generate IDs", "Error", "Sequences");
  }
  pird.id_list = ValNodeFreeData (pird.id_list);    
}


extern Boolean 
GetTableOptions 
(BaseFormPtr bfp,
 ValNodePtr clickable_list,
 CharPtr win_title,
 CharPtr label1,
 CharPtr label2,
 CharPtr skip_already_txt,
 CharPtr blanks_erase_txt,
 BoolPtr skip_already_has,
 BoolPtr blanks_erase)
{
  ClickableItemListWindowPtr drfp;
  GrouP                      h, opts = NULL;
  GrouP                      c;
  WindoW                     w;
  ButtoN                     b;
  SeqEntryPtr                sep;
  ModalAcceptCancelData      acd;
  ButtoN                     skip_already_btn = NULL, blanks_erase_btn = NULL;
  ValNodePtr                 vnp;
  ClickableItemPtr           cip;
  Boolean                    has_already = FALSE, has_blanks = FALSE;
  OMUserDataPtr              omudp;

  if (bfp == NULL) return FALSE;
    
  drfp = (ClickableItemListWindowPtr) MemNew (sizeof (ClickableItemListWindowData));
  if (drfp == NULL)
  {
    return FALSE;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = MovableModalWindow(-20, -13, -10, -10, win_title, NULL);
  SetObjectExtra (w, drfp, CleanupClickableItemListWindow);
  drfp->form = (ForM) w;
  drfp->formmessage = ClickableItemListWindowMessage;
    
  drfp->item_list = clickable_list;

  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = ClickableItemListWindowMsgFunc;
  }
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = CreateClickableListDialog (h, label1, label1,
                                                    ScrollToDiscrepancyItem, EditDiscrepancyItem, bfp,
                                                    GetDiscrepancyItemText);

  PointerToDialog (drfp->clickable_list, drfp->item_list);        

  if (skip_already_has != NULL) {
    *skip_already_has = FALSE;
  }
  if (blanks_erase != NULL) {
    *blanks_erase = FALSE;
  }

  for (vnp = clickable_list; vnp != NULL && (!has_already || !has_blanks); vnp = vnp->next) {
    cip = vnp->data.ptrvalue;
    if (cip->clickable_item_type == TABLE_DATA_ALREADY_HAS) {
      has_already = TRUE;
    } else if (cip->clickable_item_type == TABLE_DATA_CELL_BLANK) {
      has_blanks = TRUE;
    }
  }
  if (has_already || has_blanks) {
    opts = HiddenGroup (h, 0, 2, NULL);
    SetGroupSpacing (h, 10, 10);
    if (has_already) {
      if (skip_already_txt == NULL) {
        skip_already_txt = "Skip sequences that already have data";
      }
      skip_already_btn = CheckBox (opts, skip_already_txt, NULL);
    }
    if (has_blanks) {
      if (blanks_erase_txt == NULL) {
        blanks_erase_txt = "Erase data when column is blank";
      }
      blanks_erase_btn = CheckBox (opts, blanks_erase_txt, NULL);
    }
  }

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
    
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) c, (HANDLE) opts, NULL);

  Show (w);
  Select (w);
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();

  if (skip_already_btn != NULL && skip_already_has != NULL) {
    *skip_already_has = GetStatus (skip_already_btn);
  }
  if (blanks_erase_btn != NULL && blanks_erase != NULL) {
    *blanks_erase = GetStatus (blanks_erase_btn);
  }
    
  Remove (w);
  if (acd.cancelled)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }

}

static ValNodePtr GetCDSListForInterval (Int4 int_left, Int4 int_right, Uint1 strand, SeqFeatPtr first_feat)
{
  BioseqPtr         bsp;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  ValNodePtr        cds_list = NULL;

  if (first_feat == NULL) return NULL;
  bsp = BioseqFindFromSeqLoc (first_feat->location);
  if (bsp == NULL) return NULL;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
       sfp != NULL && fcontext.left <= int_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext)) {
    if (int_left <= fcontext.left
        && int_right >= fcontext.right
        && ((strand == Seq_strand_minus && fcontext.strand == Seq_strand_minus)
            || (strand != Seq_strand_minus && fcontext.strand != Seq_strand_minus))) {
      ValNodeAddPointer (&cds_list, OBJ_SEQFEAT, sfp);
    }
  }
  return cds_list;
}

static ValNodePtr GetCDSListForGene (SeqFeatPtr sfp, ValNodePtr cds_list)
{
  SeqFeatPtr        cds, gene;
  GeneRefPtr        grp;
  ValNodePtr        gene_cds_list = NULL;
  SeqMgrFeatContext fcontext;

  if (cds_list == NULL || sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return NULL;

  while (cds_list != NULL) {
    cds = cds_list->data.ptrvalue;
    if (cds != NULL) {
      grp = SeqMgrGetGeneXref (cds);
      if (grp == NULL) {
        gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
        if (gene == sfp) {
          ValNodeAddPointer (&gene_cds_list, OBJ_SEQFEAT, cds);
        }
      } else {
        if (!SeqMgrGeneIsSuppressed (grp) && GeneRefMatch (grp, sfp->data.value.ptrvalue)) {
          ValNodeAddPointer (&gene_cds_list, OBJ_SEQFEAT, cds);
        }
      }
    }
    cds_list = cds_list->next;
  }
  return gene_cds_list;
}


static Int4 GetProductLensFromProtRef (ProtRefPtr prp)
{
  Int4       lens = 0;
  ValNodePtr vnp;

  if (prp == NULL) return 0;
  for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
    lens += StringLen (vnp->data.ptrvalue) + 1;
  }
  return lens;
}

static Int4 GetProductLens (ValNodePtr cds_list)
{
  SeqFeatPtr        cds, prot_feat;
  SeqFeatXrefPtr    xref;
  Int4              lens = 0;
  BioseqPtr         prot_bsp;
  SeqMgrFeatContext fcontext;

  while (cds_list != NULL) {
    cds = cds_list->data.ptrvalue;
    xref = cds->xref;
    while (xref != NULL) {
      if (xref->data.choice == SEQFEAT_PROT) {
        lens += GetProductLensFromProtRef ((ProtRefPtr) xref->data.value.ptrvalue) + 1;
      }
      xref = xref->next;
    }
    prot_bsp = BioseqFindFromSeqLoc (cds->product);
    prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
    if (prot_feat != NULL) {
      lens += GetProductLensFromProtRef ((ProtRefPtr) prot_feat->data.value.ptrvalue) + 1;
    }
    cds_list = cds_list->next;
  }
  return lens;
}

static void AddPieceToNote (CharPtr piece, CharPtr new_note)
{
  CharPtr    cp;
  Int4       name_len;

  if (!StringHasNoText (piece)) {
    cp = StringSearch (new_note, piece);
    name_len = StringLen (piece);
    if (cp == NULL /* not found at all */
        || (cp != new_note && *(cp - 1) != ';') /* not at beginning or after semicolon */
        || (*(cp + name_len) != 0 && *(cp + name_len) != ';') /* not at end */) { 
      StringCat (new_note, piece);
      /* only add semicolon if not already at end of note */
      if (piece[name_len - 1] != ';') {
        StringCat (new_note, ";");
      }
    }
  }
}

static void AddProductNamesFromProtRef (ProtRefPtr prp, CharPtr new_note)
{
  ValNodePtr vnp;

  /* if product name other than "hypothetical protein" is available, do
   * not include "hypothetical protein" in list of product names
   */
  if (prp == NULL) return;
  for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
    if (StringCmp (vnp->data.ptrvalue, "hypothetical protein") == 0) {
      if (StringHasNoText (new_note)) {
        /* add because no other name found yet */
        AddPieceToNote (vnp->data.ptrvalue, new_note);
      }
    } else if (StringCmp (new_note, "hypothetical protein;") == 0) {
      /* replace "hypothetical protein" with new name */
      StringCpy (new_note, vnp->data.ptrvalue);
      StringCat (new_note, ";");
    } else {
      /* add as normal */
      AddPieceToNote (vnp->data.ptrvalue, new_note);
    }
  }
}

static void AddProductNames (ValNodePtr cds_list, CharPtr new_note)
{
  SeqFeatPtr        cds, prot_feat;
  SeqFeatXrefPtr    xref;
  BioseqPtr         prot_bsp;
  SeqMgrFeatContext fcontext;

  while (cds_list != NULL) {
    cds = cds_list->data.ptrvalue;
    xref = cds->xref;
    while (xref != NULL) {
      if (xref->data.choice == SEQFEAT_PROT) {
        AddProductNamesFromProtRef((ProtRefPtr) xref->data.value.ptrvalue, new_note);
      }
      xref = xref->next;
    }
    prot_bsp = BioseqFindFromSeqLoc (cds->product);
    prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
    if (prot_feat != NULL) {
      AddProductNamesFromProtRef ((ProtRefPtr) prot_feat->data.value.ptrvalue, new_note);
    }
    cds_list = cds_list->next;
  }
}

static void RemoveCDSandmRNAForGene (ValNodePtr cds_list)
{
  SeqFeatPtr        sfp, old_match;
  SeqFeatXrefPtr    xref, prev_xref, next_xref;
  Boolean           linked_by_xref;
  SeqMgrFeatContext fcontext;
  BioseqPtr         prot_bsp;

  while (cds_list != NULL) {
    sfp = cds_list->data.ptrvalue;
    sfp->idx.deleteme = TRUE;
    linked_by_xref = FALSE;
    xref = sfp->xref;
    prev_xref = NULL;
    while (xref != NULL) {
      next_xref = xref->next;
      if (xref->id.choice == 3 && xref->id.value.ptrvalue != NULL) {
        old_match = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
        if (old_match != NULL) {
          RemoveFeatureLink (sfp, old_match);
          RemoveFeatureLink (old_match, sfp);
          old_match->idx.deleteme = TRUE;
        }
        linked_by_xref = TRUE;
      } else {
        prev_xref = xref;
      }
      xref = next_xref;
    }
    if (!linked_by_xref) {
      old_match = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
      if (old_match != NULL) {
        old_match->idx.deleteme = TRUE;
      }
    }
    /* remove protein product */
    prot_bsp = BioseqFindFromSeqLoc (sfp->product);
    if (prot_bsp != NULL) {
      prot_bsp->idx.deleteme = TRUE;
    }
    cds_list = cds_list->next;
  }
}

extern void CombineToCreatePseudoGene (IteM i)
{
  BaseFormPtr                  bfp;
  SelStructPtr      sel;
  SeqFeatPtr        sfp, first_sfp;
  SeqMgrFeatContext fcontext;
  Boolean           ok_to_combine = TRUE;
  ValNodePtr        feat_list = NULL, tmp_feat;
  Int4              best_left, best_right;
  Int4              int_left = -1, int_right = -1;
  Uint1             int_strand;
  Int4              note_len, product_len;
  CharPtr           new_note = NULL, product_note = NULL;
  ValNodePtr        cds_list = NULL, gene_cds_list;
  GeneRefPtr        grp_first, grp;
  SeqIdPtr          sip;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  for (sel = ObjMgrGetSelected (); sel != NULL && ok_to_combine; sel = sel->next) {
    if (sel->entityID != bfp->input_entityID) {
      Message (MSG_ERROR, "You cannot combine genes from separate records.");
      ok_to_combine = FALSE;
    } else if (sel->itemtype != OBJ_SEQFEAT) {
      Message (MSG_ERROR, "You have selected an object that is not a feature.");
      ok_to_combine = FALSE;
    } else {
      sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp == NULL) {
        Message (MSG_ERROR, "Unable to locate selected feature.");
        ok_to_combine = FALSE;
      } else if (sfp->data.choice != SEQFEAT_GENE) {
        Message (MSG_ERROR, "You have selected a feature that is not a gene.");
        ok_to_combine = FALSE;
      } else if (feat_list == NULL
                 || best_left > fcontext.left 
                 || best_left == fcontext.left && best_right < fcontext.right) {
        best_left = fcontext.left;
        best_right = fcontext.right;
        tmp_feat = ValNodeNew (NULL);
        tmp_feat->choice = OBJ_SEQFEAT;
        tmp_feat->data.ptrvalue = sfp;
        tmp_feat->next = feat_list;
        feat_list = tmp_feat;
      } else {
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
      }

      if (int_left == -1) {
        int_left = fcontext.left;
        int_right = fcontext.right;
        int_strand = fcontext.strand;
      } else {
        if ((int_strand == Seq_strand_minus && fcontext.strand != Seq_strand_minus)
            || (int_strand != Seq_strand_minus && fcontext.strand == Seq_strand_minus)) {
          Message (MSG_ERROR, "You cannot select genes on different strands!");
          ok_to_combine = FALSE;
        } else {
          if (int_left > fcontext.left) {
            int_left = fcontext.left;
          }
          if (int_right < fcontext.right) {
            int_right = fcontext.right;
          }
        }
      }
    }
  }

  if (ok_to_combine) {
    if (feat_list == NULL || feat_list->next == NULL) {
      Message (MSG_ERROR, "You must select at least two genes to combine!");
      ok_to_combine = FALSE;
    } 
  }

  if (ok_to_combine) {
    first_sfp = feat_list->data.ptrvalue;
    grp_first = (GeneRefPtr) first_sfp->data.value.ptrvalue;

    /* get list of coding regions to use when searching for the coding regions for these genes */
    cds_list = GetCDSListForInterval (int_left, int_right, int_strand, first_sfp);
    /* get length of new gene note */
    note_len = StringLen (first_sfp->comment) + 1;

    gene_cds_list = GetCDSListForGene (first_sfp, cds_list);
    product_len = GetProductLens (gene_cds_list);
    note_len += product_len;
    gene_cds_list = ValNodeFree (gene_cds_list);
    /* get lengths of locus values, product names, and gene notes from genes after first, to use in gene note */
    for (tmp_feat = feat_list->next; tmp_feat != NULL; tmp_feat = tmp_feat->next) {
      sfp = tmp_feat->data.ptrvalue;
      grp = sfp->data.value.ptrvalue;
      if (StringCmp (grp_first->locus, grp->locus) != 0) {
        note_len += StringLen (grp->locus) + 1;
      }
      if (StringCmp (sfp->comment, first_sfp->comment) != 0) {
        note_len += StringLen (sfp->comment);
      }
      gene_cds_list = GetCDSListForGene (sfp, cds_list);
      note_len += GetProductLens (gene_cds_list) + 1;
      gene_cds_list = ValNodeFree (gene_cds_list);
    }
    new_note = (CharPtr) MemNew (sizeof (Char) * note_len);
    new_note[0] = 0;
    /* put locus tags first */
    for (tmp_feat = feat_list->next; tmp_feat != NULL; tmp_feat = tmp_feat->next) {
      sfp = tmp_feat->data.ptrvalue;
      grp = sfp->data.value.ptrvalue;
      if (StringCmp (grp_first->locus, grp->locus) != 0) {
        StringCat (new_note, grp->locus);
        StringCat (new_note, ";");
      }
    }
    /* add original gene note */\
    if (!StringHasNoText (first_sfp->comment)) {
      StringCat (new_note, first_sfp->comment);
      StringCat (new_note, ";");
    }

    /* add comments */
    for (tmp_feat = feat_list; tmp_feat != NULL; tmp_feat = tmp_feat->next) {
      sfp = tmp_feat->data.ptrvalue;
      AddPieceToNote (sfp->comment, new_note);
    }

    /* add product names */
    for (tmp_feat = feat_list; tmp_feat != NULL; tmp_feat = tmp_feat->next) {
      sfp = tmp_feat->data.ptrvalue;
      gene_cds_list = GetCDSListForGene (sfp, cds_list);
      product_note = (CharPtr) MemNew (sizeof (Char) * product_len);
      AddProductNames (gene_cds_list, product_note);
      AddPieceToNote (product_note, new_note);
      product_note = MemFree (product_note);
      gene_cds_list = ValNodeFree (gene_cds_list);
    }

    /* remove trailing semicolon, replace note */
    if (new_note[0] != 0) {
      new_note [StringLen (new_note) - 1] = 0;
      first_sfp->comment = MemFree (first_sfp->comment);
      first_sfp->comment = new_note;
    }

    /* set pseudo */
    first_sfp->pseudo = TRUE;

    /* set new location */
    sip = SeqIdDup (SeqLocId (first_sfp->location));
    first_sfp->location = SeqLocFree (first_sfp->location);
    first_sfp->location = SeqLocIntNew (int_left, int_right, int_strand, sip);
    sip = SeqIdFree (sip);

    /* remove cds and mrna for first gene */
    gene_cds_list = GetCDSListForGene (first_sfp, cds_list);
    RemoveCDSandmRNAForGene (gene_cds_list);
    gene_cds_list = MemFree (gene_cds_list);

    /* remove genes after first and related CDS and mRNA */
    for (tmp_feat = feat_list->next; tmp_feat != NULL; tmp_feat = tmp_feat->next) {
      sfp = tmp_feat->data.ptrvalue;
      gene_cds_list = GetCDSListForGene (sfp, cds_list);
      RemoveCDSandmRNAForGene (gene_cds_list);
      gene_cds_list = ValNodeFree (gene_cds_list);
      sfp->idx.deleteme = TRUE;
    }
    
    /* free interval cds list */
    cds_list = ValNodeFree (cds_list);

    /* delete objects, renormalize nuc-prot sets */
    DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    RemoveOrphanProteins (bfp->input_entityID, sep);
    RenormalizeNucProtSets (sep, TRUE);

    ObjMgrSelect (0, 0, 0, 0, NULL);

    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    Update ();
  }
  feat_list = ValNodeFree (feat_list);
}


static void FindBadLatLonSourceDesc (SeqDescrPtr sdp, Pointer userdata)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL)
  {
    return;
  }
  if (FindBadLatLon (sdp->data.ptrvalue) != NULL)
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQDESC, sdp);
  }
}


static void FindBadLatLonSourceFeat (SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL)
  {
    return;
  }
  if (FindBadLatLon (sfp->data.value.ptrvalue) != NULL)
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
  }
}



static void RefreshLatLonTool (Pointer data)
{
  BarcodeToolPtr vstp = (BarcodeToolPtr) data;
  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = ValNodeFree (vstp->item_list);

  VisitDescriptorsInSep (vstp->top_sep, &(vstp->item_list), FindBadLatLonSourceDesc);

  PointerToDialog (vstp->clickable_list, vstp->item_list);  

}


static void CleanupLatLonTool (GraphiC g, VoidPtr data)

{
  BarcodeToolPtr drfp;

  drfp = (BarcodeToolPtr) data;
  if (drfp != NULL) {
    drfp->item_list = ValNodeFree (drfp->item_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}


static void AddAltitudeToSubSourceNote (BioSourcePtr biop, CharPtr extra_text)
{
  SubSourcePtr ssp;
  CharPtr      new_note, new_note_fmt = "%s%saltitude:%s";

  if (biop == NULL || StringHasNoText (extra_text))
  {
    return;
  }

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_other)
  {
    ssp = ssp->next;
  }
  if (ssp == NULL) 
  {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_other;
    ssp->next = biop->subtype;
    biop->subtype = ssp;
  }
  new_note = (CharPtr) MemNew (sizeof (Char) * (StringLen (ssp->name)
                                                + StringLen (extra_text)
                                                + StringLen (new_note_fmt)));
  sprintf (new_note, new_note_fmt, ssp->name == NULL ? "" : ssp->name,
                                   ssp->name == NULL ? "" : "; ",
                                   extra_text);
  ssp->name = MemFree (ssp->name);
  ssp->name = new_note;
}


static void LatLonAutocorrectList (FILE *fp, ValNodePtr object_list)
{
  ValNodePtr vnp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;
  SubSourcePtr bad_ssp;
  CharPtr      fix, extra_text;

  if (fp == NULL || object_list == NULL) return;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQDESC) continue;
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_source)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      bad_ssp = FindBadLatLon (biop);
      if (bad_ssp != NULL)
      {
        fix = FixLatLonFormat (bad_ssp->name);
        if (fix != NULL) 
        {
          extra_text = StringChr (fix, ',');
          if (extra_text != NULL)
          {
            *extra_text = 0;
            extra_text++;
            while (isspace (*extra_text)) 
            {
              extra_text++;
            }
          }
          fprintf (fp, "Corrected %s to %s\n", bad_ssp->name, fix);
          bad_ssp->name = MemFree (bad_ssp->name);
          bad_ssp->name = fix;
          if (extra_text != NULL) 
          {
            AddAltitudeToSubSourceNote (biop, extra_text);
            fprintf (fp, "Moved %s to subsource note\n", extra_text);
          }
        }
        else
        {
          fprintf (fp, "Unable to correct %s\n", bad_ssp->name);
        }
      }
    }
  }
}


static void LatLonAutocorrect (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        object_list;
  LogInfoPtr        lip;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  object_list = DialogToPointer (drfp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any Lat-lon values - correct all?"))
    {
      object_list = drfp->item_list;
    }
    else
    {
      return;
    }
  }
  
  lip = OpenLog ("Lat-lon Values Corrected");
  LatLonAutocorrectList (lip->fp, object_list);

  if (object_list != drfp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  RefreshLatLonTool (drfp);
  RedrawBarcodeTool (drfp);  

  ObjMgrSetDirtyFlag (drfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, drfp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);

}


static void LatLonReport (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        vnp;
  LogInfoPtr        lip;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  SubSourcePtr      bad_ssp;
  CharPtr           fix;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  
  lip = OpenLog ("Incorrectly Formatted Lat-lon Values");
  
  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQDESC) continue;
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_source)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      bad_ssp = FindBadLatLon (biop);
      if (bad_ssp != NULL)
      {
        fix = FixLatLonFormat (bad_ssp->name);
        if (fix != NULL) 
        {
          fprintf (lip->fp, "%s (Suggested correction: %s)\n", bad_ssp->name, fix);
          fix = MemFree (fix);
        }
        else
        {
          fprintf (lip->fp, "%s (No suggested correction)\n", bad_ssp->name);
        }
        lip->data_in_log = TRUE;
      }
    }
  }

  CloseLog (lip);
  lip = FreeLog (lip);

}


static void MoveIncorrectlyFormattedLatLonToNote (ButtoN b)
{
  BarcodeToolPtr vstp;
  ValNodePtr     vnp;
  SeqDescrPtr    sdp;
  BioSourcePtr   biop;
  SubSourcePtr   note_ssp, bad_ssp, prev_note;
  CharPtr        new_txt;
  ValNodePtr     object_list;

  vstp = (BarcodeToolPtr) GetObjectExtra (b);
  if (vstp == NULL) return;

  object_list = DialogToPointer (vstp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any Lat-lon values - move all to note?"))
    {
      object_list = vstp->item_list;
    }
    else
    {
      return;
    }
  }

  for (vnp = object_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQDESC) continue;
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_source)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      bad_ssp = FindBadLatLon (biop);
      if (bad_ssp != NULL)
      {
        note_ssp = biop->subtype;
        prev_note = NULL;
        while (note_ssp != NULL && note_ssp->subtype != 255)
        {
          prev_note = note_ssp;
          note_ssp = note_ssp->next;
        }
        bad_ssp->subtype = 255;

        if (note_ssp != NULL)
        {
          new_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (note_ssp->name) + StringLen (bad_ssp->name) + 2));
          sprintf (new_txt, "%s;%s", note_ssp->name, bad_ssp->name);
          bad_ssp->name = MemFree (bad_ssp->name);
          bad_ssp->name = new_txt;
          if (prev_note == NULL)
          {
            biop->subtype = note_ssp->next;
          }
          else
          {
            prev_note->next = note_ssp->next;
          }
          note_ssp->next = NULL;
          note_ssp = SubSourceFree (note_ssp);
        }
      }
    }
  }
  if (object_list != vstp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  ObjMgrSetDirtyFlag (vstp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, vstp->input_entityID, 0, 0);
  RefreshLatLonTool (vstp);
  RedrawBarcodeTool (vstp);  
}


extern void LatLonTool (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  ValNodePtr        biosources = NULL;
  BarcodeToolPtr    drfp;
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  OMUserDataPtr            omudp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitDescriptorsInSep (sep, &biosources, FindBadLatLonSourceDesc);

  if (biosources == NULL)
  {
    Message (MSG_OK, "All lat-lon values are correctly formatted.");
    return;
  }

  biosources = ValNodeFree (biosources);

  drfp = (BarcodeToolPtr) MemNew (sizeof (BarcodeToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Lat-Lon Tool", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupLatLonTool);
  drfp->form = (ForM) w;
  drfp->formmessage = BarcodeToolMessage;

  drfp->refresh_func = RefreshLatLonTool;    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = BarcodeToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  drfp->cfg = NULL;
  drfp->undo_list = NULL;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = LatLonTestResultsDisplay (h);

  RefreshLatLonTool (drfp);
  BulkEditorCheckAllDialog (drfp->clickable_list);

  c = HiddenGroup (h, 5, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Autocorrect Values", LatLonAutocorrect);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Make Report", LatLonReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Refresh List", BarcodeRefreshButton);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Move to Note", MoveIncorrectlyFormattedLatLonToNote);
  SetObjectExtra (b, drfp, NULL);

  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);

}


static void FindLatLonCountryConflicts (SeqDescrPtr sdp, Pointer data)
{
  BioSourcePtr biop;
  SubSourcePtr ssp;
  CharPtr      orig_country = NULL, country = NULL, cp;
  Boolean      found_lat_lon = FALSE;
  FloatHi      lat, lon;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL || data == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;

  for (ssp = biop->subtype; ssp != NULL && (country == NULL || !found_lat_lon); ssp = ssp->next)
  {
    if (ssp->subtype == SUBSRC_country && !StringHasNoText (ssp->name))
    {
      orig_country = ssp->name;
      country = StringSave (orig_country);
    }
    else if (ssp->subtype == SUBSRC_lat_lon)
    {
      if (ParseLatLon (ssp->name, &lat, &lon))
      {
        found_lat_lon = TRUE;
      }
    }
  }

  cp = StringChr (country, ':');
  if (cp != NULL) 
  {
    *cp = 0;
  }

  if (found_lat_lon && IsCountryInLatLonList (country)) 
  {
    if (!TestLatLonForCountry (country, lat, lon) 
        && !TestLatLonForCountry (orig_country, lat, lon) 
        && !StringContainsBodyOfWater (orig_country))
    {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
    }
  }
  country = MemFree (country);
}


static void RefreshLatLonCountryTool (Pointer data)
{
  BarcodeToolPtr vstp = (BarcodeToolPtr) data;
  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = ValNodeFree (vstp->item_list);

  VisitDescriptorsInSep (vstp->top_sep, &(vstp->item_list), FindLatLonCountryConflicts);

  PointerToDialog (vstp->clickable_list, vstp->item_list);  

}


static void LatLonCountryReport (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        vnp;
  LogInfoPtr        lip;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  SubSourcePtr      ssp;
  CharPtr           fix, country = NULL, lat_lon = NULL;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  
  lip = OpenLog ("Incorrectly Formatted Lat-lon Values");
  
  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQDESC) continue;
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_source && sdp->data.ptrvalue != NULL)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      country = NULL;
      lat_lon = NULL;
      for (ssp = biop->subtype; ssp != NULL && (country == NULL || lat_lon == NULL); ssp = ssp->next)
      {
        if (ssp->subtype == SUBSRC_country && country == NULL)
        {
          country = ssp->name;
        }
        else if (ssp->subtype == SUBSRC_lat_lon && lat_lon == NULL)
        {
          lat_lon = ssp->name;
        }
      }
      if (country == NULL) country = "missing";
      if (lat_lon == NULL) lat_lon = "missing";
      fix = GetLatLonCountryCorrection (OBJ_SEQDESC, sdp, NULL);

      fprintf (lip->fp, "Country: %s\tLat-lon: %s\tCorrection: %s\n",
               country, lat_lon, fix == NULL ? "No suggestion" : fix);
      fix = MemFree (fix);
      lip->data_in_log = TRUE;
    }
  }

  CloseLog (lip);
  lip = FreeLog (lip);

}


static CharPtr SkipNumberInString (CharPtr str)
{
  CharPtr cp;

  cp = str;
  while (isdigit (*cp))
  {
    cp++;
  }
  if (*cp == '.')
  {
    cp++;
  }
  while (isdigit (*cp))
  {
    cp++;
  }
  return cp;
}

static Boolean GetLatLonTokens (CharPtr str, CharPtr PNTR lat, CharPtr PNTR ns, CharPtr PNTR lon, CharPtr PNTR ew)
{
  CharPtr cp;

  if (StringHasNoText (str) || lat == NULL || ns == NULL || lon == NULL || ew == NULL) return FALSE;

  cp = str;
  while (isspace (*cp)) 
  {
    cp++;
  }
  if (!isdigit (*cp)) return FALSE;
  *lat = cp;
  cp = SkipNumberInString (cp);

  while (isspace (*cp)) 
  {
    cp++;
  }

  if (*cp != 'N' && *cp != 'S') return FALSE;
  *ns = cp;
  cp++;

  while (isspace (*cp)) 
  {
    cp++;
  }

  if (!isdigit (*cp)) return FALSE;
  *lon = cp;
  cp = SkipNumberInString (cp);

  while (isspace (*cp)) 
  {
    cp++;
  }

  if (*cp != 'E' && *cp != 'W') return FALSE;
  *ew = cp;

  return TRUE;
}

static CharPtr MakeLatLonValue (CharPtr lat, Char ns, CharPtr lon, Char ew)
{
  CharPtr newval, cp;
  Int4 latlen, lonlen, len;

  if (StringHasNoText (lat) || StringHasNoText (lon))
  {
    return NULL;
  }
  cp = SkipNumberInString (lat);
  latlen = cp - lat;
  cp = SkipNumberInString (lon);
  lonlen = cp - lon;

  len = latlen + lonlen + 6;

  newval = (CharPtr) MemNew (sizeof (Char) * len);

  cp = newval;
  StringNCpy (newval, lat, latlen);
  cp += latlen;
  *(cp++) = ' ';
  *(cp++) = ns;
  *(cp++) = ' ';
  StringNCpy (cp, lon, lonlen);
  cp += lonlen;
  *(cp++) = ' ';
  *(cp++) = ew;
  *(cp++) = 0;
  
  return newval;
}


static void LatLonCountryAutocorrectList (FILE *fp, ValNodePtr object_list)
{
  ValNodePtr vnp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;
  SubSourcePtr country_ssp, lat_lon_ssp, ssp;
  FloatHi      lat, lon;
  CharPtr      country, cp;
  CharPtr      fix;

  CharPtr      pLat, pNs, pLon, pEw;

  if (fp == NULL || object_list == NULL) return;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQDESC) continue;
    sdp = vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_source)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      country_ssp = NULL;
      lat_lon_ssp = NULL;
      for (ssp = biop->subtype; ssp != NULL && (country_ssp == NULL || lat_lon_ssp == NULL); ssp = ssp->next)
      {
        if (StringHasNoText (ssp->name)) continue;
        if (ssp->subtype == SUBSRC_country && country_ssp == NULL)
        {
          country_ssp = ssp;
        }
        else if (ssp->subtype == SUBSRC_lat_lon && lat_lon_ssp == NULL && ParseLatLon (ssp->name, &lat, &lon))
        {
          lat_lon_ssp = ssp;
        }
      }
      if (country_ssp == NULL)
      {
        country = NULL;
      } else {
        country = StringSave (country_ssp->name);
        cp = StringChr (country, ':');
        if (cp != NULL) *cp = 0;
      }
      if (country != NULL && lat_lon_ssp != NULL
          && IsCountryInLatLonList (country)
          && !TestLatLonForCountry (country, lat, lon)
          && GetLatLonTokens (lat_lon_ssp->name, &pLat, &pNs, &pLon, &pEw))
      {
        fix = NULL;
        if (TestLatLonForCountry (country, -lat, lon)) 
        {
          fix = MakeLatLonValue (pLat, *pNs == 'N' ? 'S' : 'N',
                                 pLon, *pEw);
        }
        else if (TestLatLonForCountry (country, lat, -lon)) 
        {
          fix = MakeLatLonValue (pLat, *pNs,
                                 pLon, *pEw == 'E' ? 'W' : 'E');
        } else if (TestLatLonForCountry (country, lon, lat)) {
          fix = MakeLatLonValue (pLon, *pEw == 'E' ? 'N' : 'S',
                                 pLat, *pNs == 'N' ? 'E' : 'W');
        }

        if (fix != NULL) 
        {
          fprintf (fp, "Corrected %s to %s\n", lat_lon_ssp->name, fix);
          lat_lon_ssp->name = MemFree (lat_lon_ssp->name);
          lat_lon_ssp->name = fix;
        }
        else
        {
          fprintf (fp, "Unable to correct %s\n", lat_lon_ssp->name);
        }
      }
      country = MemFree (country);
    }
  }
}


static void LatLonCountryAutocorrect (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        object_list;
  LogInfoPtr        lip;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  object_list = DialogToPointer (drfp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any BioSources - correct all?"))
    {
      object_list = drfp->item_list;
    }
    else
    {
      return;
    }
  }
  
  lip = OpenLog ("Lat-lon Values Corrected");
  
  LatLonCountryAutocorrectList (lip->fp, object_list);

  if (object_list != drfp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  RefreshLatLonCountryTool (drfp);
  RedrawBarcodeTool (drfp);  

  ObjMgrSetDirtyFlag (drfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, drfp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);
}



extern void LatLonCountryTool (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  ValNodePtr        biosources = NULL;
  BarcodeToolPtr    drfp;
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  OMUserDataPtr     omudp;
  PrompT            ppt;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitDescriptorsInSep (sep, &biosources, FindLatLonCountryConflicts);

  if (biosources == NULL)
  {
    Message (MSG_OK, "No conflicts found.");
    return;
  }

  biosources = ValNodeFree (biosources);

  drfp = (BarcodeToolPtr) MemNew (sizeof (BarcodeToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Lat-Lon Country Conflict Tool", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupLatLonTool);
  drfp->form = (ForM) w;
  drfp->formmessage = BarcodeToolMessage;

  drfp->refresh_func = RefreshLatLonCountryTool;
    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = BarcodeToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  drfp->cfg = NULL;
  drfp->undo_list = NULL;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = LatLonCountryResultsDisplay (h);

  RefreshLatLonCountryTool (drfp);

  ppt = StaticPrompt (h, "Note - only hemisphere transpositions and lat-lon transpositions can be corrected",
                      0, stdLineHeight, programFont, 'c');

  c = HiddenGroup (h, 5, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Autocorrect Values", LatLonCountryAutocorrect);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Make Report", LatLonCountryReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Refresh List", BarcodeRefreshButton);
  SetObjectExtra (b, drfp, NULL);
    
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) ppt, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);

}


typedef struct specifichosttool {
  BARCODE_TOOL_BLOCK
  ButtoN caps;
  ButtoN paren;
} SpecificHostToolData, PNTR SpecificHostToolPtr;

static void CleanupSpecificHostTool (GraphiC g, VoidPtr data)

{
  SpecificHostToolPtr drfp;

  drfp = (SpecificHostToolPtr) data;
  if (drfp != NULL) {
    drfp->item_list = SpecificHostFixListFree (drfp->item_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}


static void RefreshSpecificHostTool (Pointer data)
{
  SpecificHostToolPtr vstp = (SpecificHostToolPtr) data;
  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = SpecificHostFixListFree (vstp->item_list);

  vstp->item_list = Taxon3GetSpecificHostFixesInSeqEntry (vstp->top_sep, GetStatus (vstp->caps), GetStatus (vstp->paren));

  PointerToDialog (vstp->clickable_list, vstp->item_list);  

}


static void SpecificHostAutocorrect (ButtoN b)
{
  SpecificHostToolPtr drfp;
  ValNodePtr          object_list, vnp;
  LogInfoPtr          lip;
  SpecificHostFixPtr  s;

  drfp = (SpecificHostToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  object_list = DialogToPointer (drfp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any Specific-host values - correct all?"))
    {
      object_list = drfp->item_list;
    }
    else
    {
      return;
    }
  }
  
  lip = OpenLog ("Specific-Host Values Corrected");
  for (vnp = object_list; vnp != NULL; vnp = vnp->next)
  {
    s = (SpecificHostFixPtr) vnp->data.ptrvalue;
    if (s == NULL || StringHasNoText (s->bad_specific_host)) continue;
    if (ApplyOneSpecificHostFix (vnp->data.ptrvalue))
    {
      fprintf (lip->fp, "Corrected %s (replaced %s with %s)\n",
               s->bad_specific_host, s->old_taxname, s->new_taxname);
    }
    else
    {
      fprintf (lip->fp, "Unable to correct %s\n", s->bad_specific_host);
    }
  }

  if (object_list != drfp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  RefreshSpecificHostTool (drfp);
  RedrawBarcodeTool ((BarcodeToolPtr)drfp);  

  ObjMgrSetDirtyFlag (drfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, drfp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);

}


static void SpecificHostReport (ButtoN b)
{
  SpecificHostToolPtr drfp;
  ValNodePtr          vnp;
  LogInfoPtr          lip;
  SpecificHostFixPtr  s;
  CharPtr             fix;

  drfp = (SpecificHostToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  
  lip = OpenLog ("Incorrect Specific-Host Values");
  
  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next)
  {
    s = (SpecificHostFixPtr) vnp->data.ptrvalue;
    fix = NULL;
    if (StringHasNoText (s->bad_specific_host))
    {
      continue;
    }
    else if (!StringHasNoText (s->old_taxname) && !StringHasNoText (s->new_taxname))
    {
      fix = StringSave (s->bad_specific_host);
      FindReplaceString (&fix, s->old_taxname, s->new_taxname, TRUE, TRUE);
    }
    if (fix == NULL)
    {
      fprintf (lip->fp, "%s (No suggested correction)\n", s->bad_specific_host);
    }
    else
    {
      fprintf (lip->fp, "%s (Suggestion correction: %s)\n", s->bad_specific_host, fix);
    }
    lip->data_in_log = TRUE;
    fix = MemFree (fix);
  }

  CloseLog (lip);
  lip = FreeLog (lip);
}


extern void FixSpecificHostValues (IteM i)
{
  BaseFormPtr         bfp;
  SeqEntryPtr         sep;
  SpecificHostToolPtr drfp;
  WindoW              w;
  GrouP               h, c;
  ButtoN              b;
  OMUserDataPtr       omudp;
  ValNodePtr          fix_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  fix_list = Taxon3GetSpecificHostFixesInSeqEntry (sep, FALSE, TRUE);

  if (fix_list == NULL)
  {
    Message (MSG_OK, "All specific host values look up correctly");
    return;
  }

  fix_list = SpecificHostFixListFree (fix_list);

  drfp = (SpecificHostToolPtr) MemNew (sizeof (SpecificHostToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Specific-Host Tool", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupSpecificHostTool);
  drfp->form = (ForM) w;
  drfp->formmessage = BarcodeToolMessage;

  drfp->refresh_func = RefreshSpecificHostTool;
    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = BarcodeToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  drfp->cfg = NULL;
  drfp->undo_list = NULL;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = SpecificHostResultsDisplay (h);
  drfp->caps = CheckBox (h, "Only check specific host values that start with capital letters", BarcodeRefreshButton);
  SetObjectExtra (drfp->caps, drfp, NULL);
  SetStatus (drfp->caps, TRUE);
  drfp->paren = CheckBox (h, "Check text in parentheses", BarcodeRefreshButton);
  SetObjectExtra (drfp->paren, drfp, NULL);
  SetStatus (drfp->paren, FALSE);

  RefreshSpecificHostTool (drfp);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Autocorrect Values", SpecificHostAutocorrect);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Make Report", SpecificHostReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Refresh List", BarcodeRefreshButton);
  SetObjectExtra (b, drfp, NULL);
    
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) drfp->caps, (HANDLE) drfp->paren, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);

}

static void ReportFailedTaxonomyLookups (ValNodePtr clickable_list)
{
  LogInfoPtr lip;
  ValNodePtr vnp;
  ClickableItemPtr cip;

  lip = OpenLog ("Failed Taxonomy Lookups");

  for (vnp = clickable_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL && !StringHasNoText (cip->description)) {
      fprintf (lip->fp, "%s\n", cip->description);
      lip->data_in_log = TRUE;
    }
  }
  CloseLog (lip);
  lip = FreeLog (lip);      
}


static ValNodePtr RefreshFailedTaxonomyLookups (Uint2 entityID)
{
  SeqEntryPtr       sep;
  ValNodePtr        failed_list, item_list = NULL, vnp;
  CharPtr           taxname = NULL;
  ClickableItemPtr  cip = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);

  failed_list = GetOrganismTaxLookupFailuresInSeqEntry (sep);

  /* Make ClickableItem list */ 
  for (vnp = failed_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      if (cip != NULL) {
        /* set description for previous item */
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (taxname) + 18));
        sprintf (cip->description, "%s (%d)", taxname, ValNodeLen (cip->item_list));
      }
      taxname = MemFree (taxname);
      /* get taxname to use for this list */
      taxname = vnp->data.ptrvalue;
      vnp->data.ptrvalue = NULL;
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->callback_func = BulkEditDiscrepancy;
      ValNodeAddPointer (&item_list, 0, cip);
    } else {
      /* add to clickable item list */
      if (cip != NULL) {
        ValNodeAddPointer (&(cip->item_list), vnp->choice, vnp->data.ptrvalue);
      }
    }
  }
  if (cip != NULL) {
    /* set description for last item */
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (taxname) + 18));
    sprintf (cip->description, "%s (%d)", taxname, ValNodeLen (cip->item_list));
  }
  taxname = MemFree (taxname);
  failed_list = ValNodeFree (failed_list);

  return item_list;
}


extern void ListFailedTaxonomyLookups (IteM i)
{
  BaseFormPtr       bfp;
  ValNodePtr        item_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  item_list = RefreshFailedTaxonomyLookups (bfp->input_entityID);
  if (item_list == NULL) {
    Message (MSG_OK, "No failed lookups found!");
    return;
  }

  /* display list */
  ShowClickableItemListEx (item_list, bfp, 
                           "Failed Taxonomy Lookups", "Organism Names", "BioSources",
                           ReportFailedTaxonomyLookups, RefreshFailedTaxonomyLookups);

}


static void FindTaxnames (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && data != NULL && sdp->choice == Seq_descr_source) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr FindTaxNameFixes (SeqEntryPtr sep)
{
  ValNodePtr org_list = NULL, fix_list = NULL;

  VisitDescriptorsInSep (sep, &org_list, FindTaxnames);
  fix_list = Taxon3GetTaxFixList (org_list);
  org_list = ValNodeFree (org_list);
  return fix_list;
}


static void RefreshTaxFixTool (Pointer data)
{
  BarcodeToolPtr vstp = (BarcodeToolPtr) data;

  if (vstp == NULL) return;
  
  PointerToDialog (vstp->clickable_list, NULL);
  vstp->item_list = TaxFixItemListFree (vstp->item_list);

  vstp->item_list = FindTaxNameFixes (vstp->top_sep);

  PointerToDialog (vstp->clickable_list, vstp->item_list);  

}

NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data);

static void TaxFixAutocorrectList (FILE *fp, ValNodePtr tax_fix_list)
{
  ValNodePtr vnp;
  TaxFixItemPtr t;
  BioSourcePtr biop;

  if (fp == NULL || tax_fix_list == NULL) return;

  for (vnp = tax_fix_list; vnp != NULL; vnp = vnp->next)
  {
    t = (TaxFixItemPtr) vnp->data.ptrvalue;
    if (t == NULL || t->suggested_fix == NULL) {
      continue;
    }
    biop = GetBioSourceFromObject (t->data_choice, t->data);
    if (biop != NULL) {
      if (biop->org == NULL) {
        biop->org = OrgRefNew();
      }
      fprintf (fp, "Corrected %s to %s\n", biop->org->taxname == NULL ? "missing name" : biop->org->taxname, t->suggested_fix);
      SetTaxNameAndRemoveTaxRef (biop->org, StringSave (t->suggested_fix));
    }
  }
}


static void TaxFixAutocorrect (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        object_list;
  LogInfoPtr        lip;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  ApplyBulkEditorToObjectList (drfp->clickable_list);

  object_list = DialogToPointer (drfp->clickable_list);

  if (object_list == NULL)
  {
    if (ANS_YES == Message (MSG_YN, "You have not selected any taxname values - correct all?"))
    {
      object_list = drfp->item_list;
    }
    else
    {
      return;
    }
  }
  
  lip = OpenLog ("Taxname Values Corrected");
  TaxFixAutocorrectList (lip->fp, object_list);

  if (object_list != drfp->item_list) 
  {
    object_list = ValNodeFree (object_list);
  }

  RefreshTaxFixTool (drfp);
  RedrawBarcodeTool (drfp);  

  ObjMgrSetDirtyFlag (drfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, drfp->input_entityID, 0, 0);
  Update();

  lip->data_in_log = TRUE;
  CloseLog (lip);
  lip = FreeLog (lip);

}


static void TaxFixReport (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        vnp;
  LogInfoPtr        lip;
  TaxFixItemPtr     t;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;
  
  lip = OpenLog ("Tax Name Fixes");
  
  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next)
  {
    if ((t = (TaxFixItemPtr) vnp->data.ptrvalue) == NULL) {
      continue;
    }
    if (t->suggested_fix == NULL) {
      fprintf (lip->fp, "%s (No suggested correction)\n", t->taxname);
    } else {
      fprintf (lip->fp, "%s (Suggested correction: %s)\n", t->taxname, t->suggested_fix);
      lip->data_in_log = TRUE;
    }
  }

  CloseLog (lip);
  lip = FreeLog (lip);
}


static void AddTextToFix (CharPtr text, BarcodeToolPtr drfp)
{
  ValNodePtr        vnp;
  TaxFixItemPtr     t;
  CharPtr           tmp;

  if (drfp == NULL || text == NULL) {
    return;
  }

  ApplyBulkEditorToObjectList (drfp->clickable_list);
  PointerToDialog (drfp->clickable_list, NULL);

  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next) {
    t = (TaxFixItemPtr) vnp->data.ptrvalue;
    if (t == NULL) {
      continue;
    }
    if (t->suggested_fix == NULL) {
      t->suggested_fix = (CharPtr) MemNew (sizeof (Char) * (StringLen (t->taxname) + StringLen (text) + 1));
      sprintf (t->suggested_fix, "%s%s", t->taxname, text);
    } else {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (t->suggested_fix) + StringLen (text) + 1));
      sprintf (tmp, "%s%s", t->suggested_fix, text);
      t->suggested_fix = MemFree (t->suggested_fix);
      t->suggested_fix = tmp;
    }
  }
  PointerToDialog (drfp->clickable_list, drfp->item_list);
  RedrawBarcodeTool (drfp);  

}


static void TaxFixAddSp (ButtoN b)
{
  BarcodeToolPtr    drfp;
  CharPtr           sp = " sp.";

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  AddTextToFix (sp, drfp);
}


static void TaxFixAddBacterium (ButtoN b)
{
  BarcodeToolPtr    drfp;
  CharPtr           bacterium = " bacterium";

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  AddTextToFix (bacterium, drfp);
}


static void TaxFixCopyNameToCorrection (ButtoN b)
{
  BarcodeToolPtr    drfp;
  ValNodePtr        vnp;
  TaxFixItemPtr     t;

  drfp = (BarcodeToolPtr) GetObjectExtra (b);
  if (drfp == NULL) return;

  ApplyBulkEditorToObjectList (drfp->clickable_list);
  PointerToDialog (drfp->clickable_list, NULL);

  for (vnp = drfp->item_list; vnp != NULL; vnp = vnp->next) {
    t = (TaxFixItemPtr) vnp->data.ptrvalue;
    if (t == NULL) {
      continue;
    }
    t->suggested_fix = MemFree (t->suggested_fix);
    t->suggested_fix = StringSave (t->taxname);
  }
  PointerToDialog (drfp->clickable_list, drfp->item_list);
  RedrawBarcodeTool (drfp);  
}


static void CleanupTaxFixTool (GraphiC g, VoidPtr data)

{
  BarcodeToolPtr drfp;

  drfp = (BarcodeToolPtr) data;
  if (drfp != NULL) {
    drfp->item_list = TaxFixItemListFree (drfp->item_list);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  }
  StdCleanupFormProc (g, data);
}


extern void TaxFixTool (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  ValNodePtr        biosources = NULL;
  BarcodeToolPtr    drfp;
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  OMUserDataPtr            omudp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  biosources = FindTaxNameFixes (sep);

  if (biosources == NULL)
  {
    Message (MSG_OK, "No bad taxnames.");
    return;
  }

  biosources = TaxFixItemListFree (biosources);

  drfp = (BarcodeToolPtr) MemNew (sizeof (BarcodeToolData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->top_sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Tax Fix Tool", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupTaxFixTool);
  drfp->form = (ForM) w;
  drfp->formmessage = BarcodeToolMessage;

  drfp->refresh_func = RefreshTaxFixTool;    
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = BarcodeToolMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
  
  drfp->item_list = NULL;
  drfp->cfg = NULL;
  drfp->undo_list = NULL;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  drfp->clickable_list = TaxFixDisplay (h);

  RefreshTaxFixTool (drfp);

  c = HiddenGroup (h, 7, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Apply Corrections", TaxFixAutocorrect);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Make Report", TaxFixReport);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Refresh List", BarcodeRefreshButton);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Add sp.", TaxFixAddSp);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Add bacterium", TaxFixAddBacterium);
  SetObjectExtra (b, drfp, NULL);
  b = PushButton (c, "Copy Taxname to Correction", TaxFixCopyNameToCorrection);
  SetObjectExtra (b, drfp, NULL);

  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) drfp->clickable_list, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);

}

