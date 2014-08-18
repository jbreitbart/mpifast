/*   vsmutil.c
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
* File Name:  vsmutil.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   3/3/95
*
* $Revision: 6.61 $
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

#include <vibrant.h>
#include <document.h>
#include <viewer.h>
#include <objfdef.h>
#include <prtutil.h>
#include <gather.h>
#include <vsmutil.h>
#include <dlogutil.h>
#include <bspview.h>
#include <explore.h>
#include <valid.h>
#include <subutil.h>

typedef struct errfltrdata {
  ErrSev         sev;
  int            errcode;
  int            subcode;
  CharPtr        text1;
  CharPtr        text2;
  CharPtr        text3;
} ErrFltrData, PNTR ErrFltrPtr;

typedef struct effitemdata {
  ErrSev   severity;
  int      errcode;
  int      subcode;
  Uint2    entityID;
  Uint2    itemtype;
  Uint4    itemID;
  CharPtr  sevlabel;
  CharPtr  hidden;
  CharPtr  catname;
  CharPtr  errname;
  CharPtr  accession;
  CharPtr  message;
  CharPtr  objtype;
  CharPtr  label;
  CharPtr  context;
  CharPtr  location;
  CharPtr  product;
  CharPtr  expanded;
  Pointer  userdata;
} ErrItemData, PNTR ErrItemPtr;

typedef struct validextra {
  FORM_MESSAGE_BLOCK
  ButtoN         remove;
  PopuP          verbose;
  PopuP          minlevel;
  PopuP          filter;
  DoC            doc;
  FonT           font;
  ButtoN         find;
  TexT           searchfor;
  size_t         srchtxtlen;
  PrompT         showncount;
  PrompT         summary;
  Int4           counts [6];
  Int4           totalcount;
  Int4           addedcount;
  Int4           remaining;
  Int2           clicked;
  Int2           selected;
  Boolean        dblClick;
  Boolean        shftKey;
  ErrNotifyProc  notify;
  ValNodePtr     messages;
  ValNodePtr     lastmessage;
  ValNodePtr     errorfilter;
  BaseFormPtr    bfp;
  FormActnFunc   revalProc;
  ButtoN         revalBtn;
  Boolean        okaytosetviewtarget;
  Boolean        indexerVersion;
  Int2           selected_text_start_item;
  Int2           selected_text_start_col;
  Int2           selected_text_start_row;
  Int4           selected_text_start_offset;
  Int2           selected_text_end_item;
  Int2           selected_text_end_col;
  Int2           selected_text_end_row;
  Int4           selected_text_end_offset;
  Int2           selected_text_anchor_item;
  Int2           selected_text_anchor_col;
  Int2           selected_text_anchor_row;
  Int4           selected_text_anchor_offset;
} ValidExtra, PNTR ValidExtraPtr;

static WindoW  validWindow = NULL;

static ParData valParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ParData justAccnParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};

static ColData valColFmt [] = {
  {0,  7, 15,  0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* severity */
  {0,  5, 15,  0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* error    */
  {0,  5, 45,  0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* subcode  */
  {0,  0,  0,  0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* hidden   */
  {0,  0, 60, 15, NULL, 'l', TRUE,  FALSE, FALSE, TRUE,  TRUE}   /* message  */
};

extern void FreeValidateWindow (void)

{
  if (validWindow != NULL) {
    Remove (validWindow);
    validWindow = NULL;
  }
}

static void InvalBorder (DoC d, Int2 item)

{
  Int2  bottom;
  RecT  r;
  Int2  top;

  ObjectRect (d, &r);
  InsetRect (&r, 4, 4);
  if (ItemIsVisible (d, item, &top, &bottom, NULL)) {
    r.top = top;
    r.bottom = bottom;
    r.right = r.left + 4;
    InsetRect (&r, -1, -1);
    InvalRect (&r);
  }
}


static Int2 val_columns_for_copy[] = {1, 2, 3, 5};

static void DrawTextSelection (DoC doc, ValidExtraPtr vep, Int2 item, RectPtr r)
{
  Int2 lineHeight, numRows, numCols;
  Int2 last_right;
  Int2 left_start;
  Int4 top, line_y;
  CharPtr txt;

  if (vep == NULL || r == NULL
      || item < vep->selected_text_start_item
      || item > vep->selected_text_end_item) {
    return;
  }

  if (vep->selected_text_start_item == vep->selected_text_end_item
      && vep->selected_text_start_row == vep->selected_text_end_row
      && vep->selected_text_start_offset == vep->selected_text_end_offset) {
    /* if we've only selected one char, and it's blank, don't draw it. */
    txt = GetSelectedDocText (doc, vep->selected_text_start_item, vep->selected_text_start_row,
                              vep->selected_text_start_col, vep->selected_text_start_offset,
                              vep->selected_text_end_item, vep->selected_text_end_row,
                              vep->selected_text_end_col, vep->selected_text_end_offset,
                              val_columns_for_copy, sizeof (val_columns_for_copy));
    if (StringHasNoText (txt)) {
      MemFree (txt);
      return;
    }
    MemFree (txt);
  }

  GetItemParams4 (doc, item, &top, &numRows, &numCols, &lineHeight, NULL);

  /* calculate missing rows from end first */
  if (vep->selected_text_end_item == item) {
    numRows = vep->selected_text_end_row;
    last_right = PanelOffsetFromCharOffsetEx (doc, vep->font, item, vep->selected_text_end_col, vep->selected_text_end_offset);
  } else {
    last_right = r->right;
  }

  if (vep->selected_text_start_item == item) {
    left_start = PanelOffsetFromCharOffsetEx (doc, vep->font, item, vep->selected_text_start_col, vep->selected_text_start_offset);
    line_y = r->top + (vep->selected_text_start_row) * lineHeight - 1;
    numRows -= vep->selected_text_start_row - 1;
  } else {
    left_start = PanelOffsetFromCharOffsetEx (doc, vep->font, item, 1, 0);
    line_y = r->top + lineHeight - 1;
  }

  while (numRows > 1) {
    MoveTo (left_start, line_y);
    left_start = PanelOffsetFromCharOffsetEx (doc, vep->font, item, numCols, 0);
    LineTo (r->right, line_y);
    line_y += lineHeight;
    numRows--;
  }
  MoveTo (left_start, line_y);
  LineTo (last_right, line_y);

}


static void DrawValid (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RecT           rct;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (d);
  if (vep != NULL) {
    if (item == vep->selected) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
    if (vep->indexerVersion && vep->selected_text_start_item > -1) {
      DrawTextSelection (d, vep, item, r);
    }
  }
}


static void UpdateValidTextSelection (ValidExtraPtr vep, DoC d, PoinT pt)
{
  Int2           item, row, col, offset;

  if (vep == NULL) return;

  MapDocPoint (d, pt, &item, &row, &col, NULL);

  if (item < 1 || row < 1 || col < 1) return;
  offset = GetTextSelectCharOffsetEx (pt, d, vep->font, item, row, col);
  if (item > vep->selected_text_anchor_item
      || (item == vep->selected_text_anchor_item 
          && (col > vep->selected_text_anchor_col
              || (col == vep->selected_text_anchor_col
                  && (row > vep->selected_text_anchor_row
                      || (row == vep->selected_text_anchor_row
                          && offset >= vep->selected_text_anchor_offset)))))) {
    vep->selected_text_start_item = vep->selected_text_anchor_item;
    vep->selected_text_start_row = vep->selected_text_anchor_row;
    vep->selected_text_start_col = vep->selected_text_anchor_col;
    vep->selected_text_start_offset = vep->selected_text_anchor_offset;
    vep->selected_text_end_item = item;
    vep->selected_text_end_row = row;
    vep->selected_text_end_col = col;
    vep->selected_text_end_offset = offset;
  } else {
    vep->selected_text_start_item = item;
    vep->selected_text_start_row = row;
    vep->selected_text_start_col = col;
    vep->selected_text_start_offset = offset;
    vep->selected_text_end_item = vep->selected_text_anchor_item;
    vep->selected_text_end_row = vep->selected_text_anchor_row;
    vep->selected_text_end_col = vep->selected_text_anchor_col;
    vep->selected_text_end_offset = vep->selected_text_anchor_offset;
  }
  InvalDocRows (d, 0, 0, 0);
} 


static void ClickValid (DoC d, PoinT pt)

{
  Int2           item;
  Int2           row, col;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (d);
  if (vep != NULL) {
    MapDocPointEx (d, pt, &item, &row, &col, NULL, TRUE);
    if (item > 0 && row > 0 && vep->clicked == item) {
      vep->dblClick = dblClick;
    } else {
      vep->dblClick = FALSE;
    }
    vep->clicked = 0;
    vep->shftKey = shftKey;
    if (item > 0 && row > 0) {
      vep->clicked = item;
      if (vep->indexerVersion) {
        vep->selected_text_anchor_item = item;
        vep->selected_text_anchor_col = col;
        vep->selected_text_anchor_row = row;
        vep->selected_text_anchor_offset = GetTextSelectCharOffsetEx (pt, d, vep->font, item, row, col);
        UpdateValidTextSelection (vep, d, pt);
      }
    }
  }
}


static void DragValid (DoC d, PoinT pt)
{
  ValidExtraPtr   vep;

  vep = GetObjectExtra (d);
  if (vep != NULL && vep->indexerVersion) {
    UpdateValidTextSelection (vep, d, pt);
  }
}


static Boolean FindSfpItem (GatherContextPtr gcp)

{
  SeqFeatPtr  PNTR sfpp;

  sfpp = (SeqFeatPtr PNTR) gcp->userdata;
  if (sfpp != NULL && gcp->thistype == OBJ_SEQFEAT) {
    *sfpp = (SeqFeatPtr) gcp->thisitem;
  }
  return TRUE;
}

static SeqFeatPtr GetSeqFeatGivenIDs (Uint2 entityID, Uint4 itemID, Uint2 itemtype)

{
  SeqFeatPtr  sfp;

  sfp = NULL;
  if (entityID > 0 && itemID > 0 && itemtype == OBJ_SEQFEAT) {
    GatherItem (entityID, itemID, itemtype, (Pointer) (&sfp), FindSfpItem);
  }
  return sfp;
}

static void ValDoNotify (ValidExtraPtr vep, Int2 item, Boolean select, Boolean target)

{
  BioseqPtr     bsp = NULL;
  unsigned int  entityID;
  int           errcode;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       ptr;
  Char          seqid [128];
  int           sev;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  int           subcode;
  CharPtr       str;

  if (vep != NULL && vep->doc != NULL && vep->notify != NULL) {
    str = GetDocText (vep->doc, item, 1, 4);
    if (str != NULL &&
        sscanf (str, "%d %d %d %u %u %u", &sev, &errcode,
                &subcode, &entityID, &itemID, &itemtype) == 6) {
      if (target) {
        if (itemtype == OBJ_BIOSEQ || itemtype == OBJ_SEQFEAT) {
          if (itemtype == OBJ_BIOSEQ) {
            bsp = GetBioseqGivenIDs (entityID, itemID, itemtype);
          } else {
            sfp = GetSeqFeatGivenIDs (entityID, itemID, itemtype);
            if (sfp != NULL)
            {
              bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
            }
          }
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
            SeqIdWrite (sip, seqid, PRINTID_REPORT, sizeof (seqid));
            ptr = StringChr (seqid, '|');
            if (ptr == NULL) {
              ptr = seqid;
            } else {
              ptr++;
            }
            SetBioseqViewTarget (vep->bfp, ptr);
          }
        }
      }
      if (select) {
        ObjMgrSelect (entityID, itemID, itemtype, 0, NULL);
      } else {
        ObjMgrDeSelect (entityID, itemID, itemtype, 0, NULL);
      }
      (vep->notify) ((ErrSev) sev, errcode, subcode,
                     (Uint2) entityID, itemID, (Uint2) itemtype,
                     select, vep->dblClick);
    }
    MemFree (str);
  }
}

static void ReleaseValid (DoC d, PoinT pt)

{
  Int2           item;
  Int2           old;
  Int2           row;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (d);
  if (vep != NULL) {
    ResetClip ();
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {
      if (item == vep->clicked) {
        old = vep->selected;
        vep->selected = item;
        if (old != item) {
          if (old == 0) {
            InvalBorder (d, item);
          } else {
            InvalBorder (d, MIN (item, old));
            InvalBorder (d, MAX (item, old));
            ValDoNotify (vep, old, FALSE, FALSE);
          }
          if (! Enabled (vep->remove)) {
            Enable (vep->remove);
          }
          Update ();
        }
        ValDoNotify (vep, item, TRUE, vep->okaytosetviewtarget && vep->shftKey);
      }
    } else if (vep->clicked == 0) {
      if (vep->selected != 0) {
        old = vep->selected;
        vep->selected = 0;
        InvalBorder (d, old);
        ValDoNotify (vep, old, FALSE, FALSE);
      }
      if (Enabled (vep->remove)) {
        Disable (vep->remove);
      }
      Update ();
    }
  }
}

/*
static void RemoveProc (ButtoN b)

{
  Int2           numItems;
  Int2           old;
  RecT           r;
  ValidExtraPtr  vep;
  WindoW         w;

  vep = GetObjectExtra (b);
  if (vep != NULL) {
    if (vep->selected != 0) {
      GetDocParams (vep->doc, &numItems, NULL);
      if (vep->selected <= numItems) {
        old = vep->selected;
        vep->selected = 0;
        DeleteItem (vep->doc, old);
        if (Enabled (vep->doc) && AllParentsEnabled (vep->doc) &&
            Visible (vep->doc) && AllParentsVisible (vep->doc)) {
          w = CurrentWindow ();
          Select (vep->doc);
          ObjectRect (vep->doc, &r);
          InsetRect (&r, 3, 3);
          InvalRect (&r);
          UseWindow (w);
        }
        UpdateDocument (vep->doc, 0, 0);
        Disable (b);
      }
    }
  }
}
*/

static void CloseValidWindow (WindoW w)

{
  FreeValidateWindow ();
}

static void CloseValidButton (ButtoN b)

{
  FreeValidateWindow ();
}

static ValNodePtr FreeErrItemList (ValNodePtr head)

{
  ErrItemPtr  eip;
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    eip = (ErrItemPtr) vnp->data.ptrvalue;
    if (eip == NULL) continue;
    MemFree (eip->sevlabel);
    MemFree (eip->hidden);
    MemFree (eip->catname);
    MemFree (eip->errname);
    MemFree (eip->accession);
    MemFree (eip->message);
    MemFree (eip->objtype);
    MemFree (eip->label);
    MemFree (eip->context);
    MemFree (eip->location);
    MemFree (eip->product);
    MemFree (eip->expanded);
    MemFree (eip);
  }
  return ValNodeFree (head);
}

extern void ClearValidateWindow (void)

{
  ErrFltrPtr     efp;
  int            sev;
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  if (validWindow != NULL) {
    vep = GetObjectExtra (validWindow);
    if (vep != NULL) {
      Reset (vep->doc);
      vep->selected = 0;
      vep->messages = FreeErrItemList (vep->messages);
      vep->lastmessage = NULL;
      for (vnp = vep->errorfilter; vnp != NULL; vnp = vnp->next) {
        efp = (ErrFltrPtr) vnp->data.ptrvalue;
        if (efp != NULL) {
          MemFree (efp->text1);
          MemFree (efp->text2);
          MemFree (efp->text3);
        }
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      }
      vep->errorfilter = ValNodeFree (vep->errorfilter);
      Hide (vep->filter);
      Update ();
      Reset (vep->filter);
      PopupItem (vep->filter, "ALL");
      SetValue (vep->filter, 1);
      Show (vep->filter);
      SetTitle (vep->summary, "");
      SetTitle (vep->showncount, "");
      for (sev = SEV_NONE; sev <= SEV_MAX; sev++) {
        vep->counts [sev] = 0;
      }
      vep->totalcount = 0;
      vep->addedcount = 0;
      vep->remaining = 0;
      SetTitle (vep->searchfor, "");
      vep->srchtxtlen = 0;
      if (Enabled (vep->find)) {
        Disable (vep->find);
      }
      Update ();
    }
  }
}

static Boolean doSuppressContext = TRUE;

extern Boolean ShouldSetSuppressContext (void)

{
  return doSuppressContext;
}

static Boolean doJustShowAccession = FALSE;

extern Boolean ShouldSetJustShowAccession (void)

{
  return doJustShowAccession;
}

static CharPtr StringMoveConvertTabs (CharPtr dst, CharPtr src)

{
  Char  ch;

  if (src == NULL || dst == NULL) return dst;
  ch = *src;
  while (ch != '\0') {
    if (ch < ' ') {
      ch = ' ';
    }
    *dst = ch;
    dst++;
    src++;
    ch = *src;
  }
  *dst = '\0';
  return dst;
}

static CharPtr CharPrtProc (DoC d, Int2 item, Pointer ptr)

{
  Pointer        data;
  ErrItemPtr     eip;
  size_t         len;
  CharPtr        msg_text;
  ValidExtraPtr  vep;
  Int2           verbose;

  vep = GetObjectExtra (d);
  if (vep == NULL) return NULL;
  GetItemParams4 (vep->doc, item, NULL, NULL, NULL, NULL, &data);
  if (data == NULL) return NULL;
  eip = (ErrItemPtr) data;
  verbose = GetValue (vep->verbose);

  len = StringLen (eip->sevlabel) +
        StringLen (eip->hidden) +
        StringLen (eip->catname) +
        StringLen (eip->errname) +
        StringLen (eip->accession) +
        StringLen (eip->message) +
        StringLen (eip->objtype) +
        StringLen (eip->label) +
        StringLen (eip->context) +
        StringLen (eip->location) +
        StringLen (eip->product) +
        StringLen (eip->expanded);

  msg_text = (CharPtr) MemNew (sizeof (Char) * (len + 50));
  if (msg_text == NULL) return NULL;
  ptr = msg_text;

  if (verbose == 4) {
    ptr = StringMoveConvertTabs (ptr, eip->accession);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->sevlabel);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->catname);
    ptr = StringMove (ptr, "_");
    ptr = StringMoveConvertTabs (ptr, eip->errname);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->hidden);
    ptr = StringMove (ptr, "\n");
  } else {
    ptr = StringMoveConvertTabs (ptr, eip->sevlabel);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->catname);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->errname);
    ptr = StringMove (ptr, "\t");
    ptr = StringMoveConvertTabs (ptr, eip->hidden);
    ptr = StringMove (ptr, "\n\t\t\t\t");
    ptr = StringMoveConvertTabs (ptr, eip->message);
    if (StringDoesHaveText (eip->objtype)) {
      ptr = StringMove (ptr, " ");
      ptr = StringMoveConvertTabs (ptr, eip->objtype);
    }
    if (StringDoesHaveText (eip->label)) {
      ptr = StringMove (ptr, ": ");
      ptr = StringMoveConvertTabs (ptr, eip->label);
    }
    if (StringDoesHaveText (eip->location)) {
      ptr = StringMove (ptr, " ");
      ptr = StringMoveConvertTabs (ptr, eip->location);
    }
    if (StringDoesHaveText (eip->product)) {
      ptr = StringMove (ptr, " -> ");
      ptr = StringMoveConvertTabs (ptr, eip->product);
    }
    if (verbose < 2 && StringDoesHaveText (eip->context)) {
      ptr = StringMove (ptr, " ");
      ptr = StringMoveConvertTabs (ptr, eip->context);
    } else if (verbose > 1 && StringDoesHaveText (eip->accession)) {
      ptr = StringMove (ptr, " (CONTEXT ");
      ptr = StringMoveConvertTabs (ptr, eip->accession);
      ptr = StringMove (ptr, ")");
    }
    if (verbose < 3 && StringDoesHaveText (eip->expanded)) {
      ptr = StringMove (ptr, "\n\n\t\t\t\t");
      ptr = StringMoveConvertTabs (ptr, eip->expanded);
    }
    ptr = StringMove (ptr, "\n");
  }

  return StringSave (msg_text);
}

static void RepopVal (PopuP p)

{
  ErrFltrPtr     efp;
  ErrItemPtr     eip;
  Int2           filt;
  size_t         len;
  Int2           minlev;
  Boolean        okay;
  Char           tmp [32];
  Int2           val;
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  vep = GetObjectExtra (p);
  if (vep != NULL) {
    Reset (vep->doc);
    vep->addedcount = 0;
    vep->selected = 0;
    minlev = GetValue (vep->minlevel);
    filt = GetValue (vep->filter);
    val = GetValue (vep->verbose);
    efp = NULL;
    if (vep->errorfilter != NULL && filt > 1) {
      vnp = vep->errorfilter;
      while (filt > 2 && vnp != NULL) {
        vnp = vnp->next;
        filt--;
      }
      if (vnp != NULL) {
        efp = (ErrFltrPtr) vnp->data.ptrvalue;
      }
    }
    for (vnp = vep->messages; vnp != NULL; vnp = vnp->next) {
      eip = (ErrItemPtr) vnp->data.ptrvalue;
      if (eip != NULL) {
        if (eip->severity >= minlev || eip->severity == SEV_NONE) {
          okay = FALSE;
          if (efp != NULL) {
            if (efp->errcode == eip->errcode) {
              if (efp->subcode == INT_MIN || efp->subcode == eip->subcode) {
                okay = TRUE;
              }
            }
          } else {
            okay = TRUE;
          }
          if (okay) {
            (vep->addedcount)++;
            len = StringLen (eip->accession) +
                  StringLen (eip->message) +
                  StringLen (eip->context) +
                  StringLen (eip->location) +
                  StringLen (eip->product);
            if (val == 4) {
              AppendItem (vep->doc, CharPrtProc, (Pointer) eip, FALSE,
                          (len / 50) + 1,
                          &justAccnParFmt, valColFmt, vep->font);
            } else {
              AppendItem (vep->doc, CharPrtProc, (Pointer) eip, FALSE,
                          (len / 50) + 1,
                          &valParFmt, valColFmt, vep->font);
            }
          }
        }
      }
    }
    if (vep->addedcount > 1) {
      sprintf (tmp, " %ld items shown", (long) vep->addedcount);
    } else if (vep->addedcount > 0) {
      sprintf (tmp, " %ld item shown", (long) vep->addedcount);
    } else {
      StringCpy (tmp, "");
    }
    SafeSetTitle (vep->showncount, tmp);
    UpdateDocument (vep->doc, 0, 0);
  }
}

static CharPtr FormatConsensusSpliceReport (CharPtr doc_line)
{
  CharPtr cp, cp2, feat_start, feat_end = NULL;
  CharPtr report_str = NULL;
  CharPtr msg_abbrev = NULL;
  Char    ch, ch_2;
  Int4    pos;

  cp = StringChr (doc_line, '\n');
  if (StringSearch (cp, "(AG) not found") != NULL) {
    msg_abbrev = "AG";
  } else if (StringSearch (cp, "(GT) not found") != NULL) {
    msg_abbrev = "GT";
  } else {
    /* no report */
    return NULL;
  }

  cp = StringSearch (cp, "position ");
  cp2 = StringSearch (cp, "FEATURE:");
  if (cp2 != NULL) {
    feat_start = cp2 + 8;
    if (feat_start != NULL) {
      feat_end = StringChr (feat_start, ':');
      if (feat_end != NULL) {
        ch_2 = *feat_end;
        *feat_end = 0;
      }
    }
    ch = *cp2;
    *cp2 = 0;
    if (feat_start == NULL) {
      report_str = (CharPtr) MemNew (sizeof (Char) * StringLen (cp));
    } else {
      report_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (cp) + StringLen (feat_start) + 2));
    }
    if (sscanf (cp, "position %ld of %s", &pos, report_str) == 2) {
      if (feat_start == NULL) {
        sprintf (report_str + StringLen (report_str), "\t%s at %ld", msg_abbrev, pos);
      } else {
        sprintf (report_str + StringLen (report_str), "\t%s\t%s at %ld", feat_start, msg_abbrev, pos);
      }
    } else {
      report_str = MemFree (report_str);
    }
    if (feat_end != NULL) {
      *feat_end = ch_2;
    }
    *cp2 = ch;
  }

  return report_str;
}

static CharPtr GetEcNumberReport (SeqFeatPtr sfp)
{
  SeqIdPtr  sip = NULL;
  BioseqPtr bsp;
  CharPtr   ec_number = NULL, tmp;
  GBQualPtr gbq;
  ProtRefPtr prp;
  ValNodePtr vnp;
  CharPtr    str;
  SeqFeatPtr bsp_feat = NULL, gene;
  GeneRefPtr grp;
  SeqMgrFeatContext fcontext;
  CharPtr           report_str = NULL;
  Int4              report_len;
  Char              seqid [128];

  if (sfp == NULL) return NULL;
  /* want: accession number for sequence, ec number, locus tag */

  /* find accession number */
  if (sfp->idx.subtype == FEATDEF_PROT) {
    /* if protein, want ID for coding region location */
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      bsp_feat = SeqMgrGetCDSgivenProduct (bsp, NULL);
    }
  }
  if (bsp_feat == NULL) {
    bsp_feat = sfp;
  }

  bsp = BioseqFindFromSeqLoc (bsp_feat->location);
  if (bsp == NULL) {
    bsp = BioseqFindFromSeqLoc (SeqLocFindNext (bsp_feat->location, NULL));
  }
  sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  if (sip == NULL) {
    sprintf (seqid, "Unknown location");
  } else {
    SeqIdWrite (sip, seqid, PRINTID_REPORT, sizeof (seqid));
  }

  /* find EC number in quals */
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "EC_number") == 0) {
      if (ec_number == NULL) {
        ec_number = StringSave (gbq->val);
      } else {
        tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (ec_number) + StringLen (gbq->val) + 2));
        sprintf (tmp, "%s;%s", ec_number, gbq->val);
        ec_number = MemFree (ec_number);
        ec_number = tmp;
      }
    }
  }

  /* find EC number in protref */
  if (sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (ec_number == NULL) {
        ec_number = StringSave (str);
      } else {
        tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (ec_number) + StringLen (str) + 2));
        sprintf (tmp, "%s;%s", ec_number, str);
        ec_number = MemFree (ec_number);
        ec_number = tmp;
      }
    }
  }

  if (StringHasNoText (ec_number)) {
    ec_number = MemFree (ec_number);
  }

  if (ec_number == NULL) {
    ec_number = StringSave ("Blank EC number");
  }

  /* get locus tag */
  grp = SeqMgrGetGeneXref (bsp_feat);
  if (grp == NULL) {
    gene = SeqMgrGetOverlappingGene (bsp_feat->location, &fcontext);
    if (gene != NULL && gene->idx.subtype == FEATDEF_GENE) {
      grp = gene->data.value.ptrvalue;
    }
  } else if (SeqMgrGeneIsSuppressed (grp)) {
    grp = NULL;
  }
  
  report_len = StringLen (seqid) + StringLen (ec_number) + 4;
  if (grp != NULL) {
    report_len += StringLen (grp->locus_tag);
  }
  report_str = (CharPtr) MemNew (sizeof (Char) * report_len);
  sprintf (report_str, "%s\t%s\t%s\n", seqid, ec_number, (grp != NULL && grp->locus_tag != NULL) ? grp->locus_tag : "");
  ec_number = MemFree (ec_number);
  return report_str;
}


static Boolean WriteOneBadSpecificHostSeqEntry (FILE *fp, CharPtr spec_host, SeqEntryPtr sep);

static Boolean WriteOneBadSpecificHostBioseq (FILE *fp, CharPtr spec_host, BioseqPtr bsp)
{
  Char         id_str[100];

  if (fp == NULL || spec_host == NULL || bsp == NULL || ISA_aa (bsp->mol)) return FALSE;

  SeqIdWrite (SeqIdFindBest (bsp->id, 0), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
  fprintf (fp, "%s\t%s\n", id_str, spec_host);
  return TRUE;
}

static Boolean WriteOneBadSpecificHostBioseqSet (FILE *fp, CharPtr spec_host, BioseqSetPtr bssp)
{
  SeqEntryPtr sep;
  Boolean     rval = FALSE;

  if (fp == NULL || spec_host == NULL || bssp == NULL || bssp->seq_set == NULL) return FALSE;

  if (bssp->_class == BioseqseqSet_class_nuc_prot
      || bssp->_class == BioseqseqSet_class_segset)
  {
    rval = WriteOneBadSpecificHostSeqEntry (fp, spec_host, bssp->seq_set);
  }
  else
  {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
      rval |= WriteOneBadSpecificHostSeqEntry (fp, spec_host, sep);
    }
  }
  return rval;
}

static Boolean WriteOneBadSpecificHostSeqEntry (FILE *fp, CharPtr spec_host, SeqEntryPtr sep)
{
  if (sep == NULL)
  {
    return FALSE;
  }
  else if (IS_Bioseq (sep))
  {
    return WriteOneBadSpecificHostBioseq (fp, spec_host, (BioseqPtr) sep->data.ptrvalue);
  }
  else if (IS_Bioseq_set (sep))
  {
    return WriteOneBadSpecificHostBioseqSet (fp, spec_host, (BioseqSetPtr) sep->data.ptrvalue);
  }
  else
  {
    return FALSE;
  }
}


static void GetParentForBioSource (Uint2 parenttype, Pointer parentptr, BioseqPtr PNTR bsp, BioseqSetPtr PNTR bssp)
{
  SeqAnnotPtr sap;
  if (bsp == NULL || bssp == NULL) return;

  *bsp = NULL;
  *bssp = NULL;
  if (parentptr == NULL) return;

  if (parenttype == OBJ_BIOSEQ)
  {
    *bsp = parentptr;        
  }
  else if (parenttype == OBJ_BIOSEQSET)
  {
    *bssp = parentptr;        
  }
  else if (parenttype == OBJ_SEQANNOT)
  {
    sap = (SeqAnnotPtr) parentptr;
    GetParentForBioSource (sap->idx.parenttype, sap->idx.parentptr, bsp, bssp);
  }
}


extern Boolean WriteBadSpecificHostTable (ValNodePtr bad_biop_list, FILE *fp)
{
  ValNodePtr   vnp;
  Boolean      any_in_list = FALSE;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  BioSourcePtr biop;
  ObjValNodePtr ovp;
  CharPtr       spec_host = NULL;
  OrgModPtr     mod = NULL;

  if (fp == NULL) return FALSE;

  for (vnp = bad_biop_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = NULL;
    bssp = NULL;
    if (vnp->choice == OBJ_SEQFEAT)
    {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      GetParentForBioSource (sfp->idx.parenttype, sfp->idx.parentptr, &bsp, &bssp);
    } 
    else if (vnp->choice == OBJ_SEQDESC)
    {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (sdp != NULL && sdp->extended != 0) 
      {
        ovp = (ObjValNodePtr) sdp;

        GetParentForBioSource (ovp->idx.parenttype, ovp->idx.parentptr, &bsp, &bssp);
      }
    }
    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;

    spec_host = NULL;
    mod = biop->org->orgname->mod;
    while (mod != NULL && mod->subtype != ORGMOD_nat_host)
    {
      mod = mod->next;
    }
    if (mod == NULL || StringHasNoText (mod->subname)) continue;
    
    spec_host = mod->subname;

    if (bsp != NULL )
    {
      any_in_list = WriteOneBadSpecificHostBioseq (fp, spec_host, bsp);
    }
    else if (bssp != NULL)
    {
      any_in_list = WriteOneBadSpecificHostBioseqSet (fp, spec_host, bssp);
    }    
  }

  return any_in_list;
}


typedef struct validatoreporttype {
  ButtoN PNTR btn_array;
  ValNodePtr  errorfilter;
} ValidatorReportTypeData, PNTR ValidatorReportTypePtr;


static void EnableValidatorReportTypeButtons (ButtoN b)
{
  ValidatorReportTypePtr data;
  ValNodePtr             vnp, vnp2;
  Int4                   pos, pos2;
  ErrFltrPtr             efp;

  data = (ValidatorReportTypePtr) GetObjectExtra (b);
  if (data == NULL) return;

  for (vnp = data->errorfilter, pos = 0;
       vnp != NULL;
       vnp = vnp->next, pos++)
  {
    if (data->btn_array[pos] == b)
    {
      efp = (ErrFltrPtr) vnp->data.ptrvalue;
      if (efp != NULL && efp->subcode == INT_MIN)
      {
        if (GetStatus (b)) 
        {
          /* main category checked, disable subcategories */
          for (vnp2 = vnp->next, pos2 = pos + 1;
               vnp2 != NULL && (vnp2->data.ptrvalue == NULL 
                                || ((ErrFltrPtr)vnp2->data.ptrvalue)->errcode == efp->errcode);
               vnp2 = vnp2->next, pos2++)
          {
            Disable (data->btn_array[pos2]);
          }
        } 
        else
        {
          /* main category unchecked, enable subcategories */
          for (vnp2 = vnp->next, pos2 = pos + 1;
               vnp2 != NULL && (vnp2->data.ptrvalue == NULL || ((ErrFltrPtr)vnp2->data.ptrvalue)->errcode == efp->errcode);
               vnp2 = vnp2->next, pos2++)
          {
            Enable (data->btn_array[pos2]);
          }
        } 
      }
      break;
    }
  }
}


static ValNodePtr CollectValidatorReportTypes (ValidExtraPtr vep)
{
  ErrFltrPtr  efp;
  ValNodePtr  chosen = NULL, vnp;
  WindoW      w, h, btn_grp, g1 = NULL, c;
  Int4        num_buttons, i;
  ValidatorReportTypeData data;
  ButtoN PNTR btn_array;
  ButtoN      b;
  int         last_errcode = 0;
  ModalAcceptCancelData acd;
  
  if (vep == NULL || vep->errorfilter == NULL) return NULL;

  /* if only one, just select the one */
  if (vep->errorfilter->next == NULL) 
  {
    ValNodeAddPointer (&chosen, 0, vep->errorfilter->data.ptrvalue);
    return chosen;
  }

  num_buttons = ValNodeLen (vep->errorfilter);
  btn_array = (ButtoN PNTR) MemNew (sizeof (ButtoN) * num_buttons);

  data.btn_array = btn_array;
  data.errorfilter = vep->errorfilter;

  w = MovableModalWindow (-50, -33, -10, -10, "Choose Report Items", NULL);
  SetGroupSpacing (w, 10, 10);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  btn_grp = HiddenGroup (h, 3, 0, NULL);
  for (vnp = vep->errorfilter, i = 0; vnp != NULL; vnp = vnp->next, i++) {
    efp = (ErrFltrPtr) vnp->data.ptrvalue;
    if (efp == NULL) continue;
    if (efp->subcode == INT_MIN) {
      g1 = NormalGroup (btn_grp, 0, 10, efp->text2, programFont, NULL);
      btn_array[i] = CheckBox (g1, "All", EnableValidatorReportTypeButtons);
      SetObjectExtra (btn_array[i], &data, NULL);
      last_errcode = efp->errcode;
    } else {
      if (last_errcode != efp->errcode) {
        g1 = NULL;
      }
      btn_array[i] = CheckBox (g1 == NULL ? btn_grp : g1, efp->text3 == NULL ? "" : efp->text3, NULL);
    }
  }

  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  b = DefaultButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) btn_grp, (HANDLE) c, NULL);
  RealizeWindow (w);

  Show (w);
  Select (w);
  Update ();

  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  while (!acd.accepted && ! acd.cancelled)
  {
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();

    if (acd.accepted)
    {
      for (vnp = vep->errorfilter, i = 0; vnp != NULL; vnp = vnp->next, i++)
      {
        if (vnp->data.ptrvalue != NULL && GetStatus (btn_array[i])) 
        {
          ValNodeAddPointer (&chosen, 0, vnp->data.ptrvalue);
        }
      }
    }
  }
  btn_array = MemFree (btn_array);
  Remove (w);
  return chosen;
}


static Int4 FindReportPositionForError (ErrItemPtr eip, ValNodePtr chosen)
{
  Int4       pos;
  ValNodePtr vnp;
  ErrFltrPtr efp;

  if (eip == NULL || chosen == NULL) return -1;

  for (vnp = chosen, pos = 0; vnp != NULL; vnp = vnp->next, pos++)
  {
    efp = (ErrFltrPtr) vnp->data.ptrvalue;
    if (efp != NULL && efp->errcode == eip->errcode 
        && (efp->subcode == INT_MIN || efp->subcode == eip->subcode))
    {
      return pos;
    }
  }
  return -1;
}


static Boolean MakeStandardReports (ValNodePtr chosen, ValNodePtr PNTR reports_list, FILE *fp)
{
  ValNodePtr chosen_vnp, item_vnp;
  ErrFltrPtr efp;
  Int4       i;
  CharPtr    cp;
  Boolean    found_any = FALSE;
  CharPtr    label;

  if (chosen == NULL || reports_list == NULL || fp == NULL) return FALSE;
  
  for (chosen_vnp = chosen, i = 0; chosen_vnp != NULL; chosen_vnp = chosen_vnp->next, i++)
  {
    efp = (ErrFltrPtr) chosen_vnp->data.ptrvalue;
    if (efp != NULL && reports_list[i] != NULL) 
    {
      fprintf (fp, "%s%s%s\n", efp->text2 == NULL ? "" : efp->text2,
                               efp->subcode == INT_MIN || efp->text3 != NULL ? ":" : "",
                               efp->subcode == INT_MIN ? "ALL" : efp->text3 == NULL ? "" : efp->text3);
      for (item_vnp = reports_list[i]; item_vnp != NULL; item_vnp = item_vnp->next)
      {
        cp = GetDiscrepancyItemText (item_vnp);
        if (cp != NULL)
        {
          label = GetParentLabelForDiscrepancyItem (item_vnp);
          if (label != NULL) {
            fprintf (fp, "%s:", label);
            label = MemFree (label);
          }
          fprintf (fp, "%s", cp);
          found_any = TRUE;
          cp = MemFree (cp);
        }
      }
      fprintf (fp, "\n");
    }
  }
  return found_any;
}


static void MakeValidatorReport (ButtoN b)
{
  ValidExtraPtr  vep;
  ErrItemPtr     eip;
  Int2           item = 0;
  ValNodePtr     vnp;
  CharPtr        cp;
  Char           path [PATH_MAX];
  Boolean        found_any = FALSE;
  FILE           *fp;
  ValNodePtr     consensus_splice_list = NULL;
  ValNodePtr     ecnumber_list = NULL;
  ValNodePtr     specific_host_list = NULL;
  SeqFeatPtr     sfp;
  SeqDescrPtr    sdp;
  BioseqPtr      bsp;
  CharPtr        str;
  Char           ch;
  CharPtr        tmp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  ValNodePtr        chosen = NULL;
  ValNodePtr PNTR   report_lists = NULL;
  Int4              num_reports, pos;

  vep = (ValidExtraPtr) GetObjectExtra (b);
  if (vep == NULL) return;

  chosen = CollectValidatorReportTypes (vep);
  if (chosen == NULL) return;
  num_reports = ValNodeLen (chosen);
  report_lists = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_reports);
  MemSet (report_lists, 0, sizeof (ValNodePtr) * num_reports);

  TmpNam (path); 
#ifdef WIN_MAC
  fp = FileOpen (path, "r");
  if (fp != NULL) {
    FileClose (fp);
  } else {
    FileCreate (path, "TEXT", "ttxt");
  }
#endif
  fp = FileOpen (path, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open temporary file for report");
  } else {
    for (vnp = vep->messages; vnp != NULL; vnp = vnp->next) {
      item++;
      eip = (ErrItemPtr) vnp->data.ptrvalue;
      if (eip != NULL) {
        pos = FindReportPositionForError (eip, chosen);
        if (pos == -1) continue;

        if (eip->errcode == 5 && (eip->subcode == 16 || eip->subcode == 137 || eip->subcode == 138 || eip->subcode == 139)) {
          /* ERR_SEQ_FEAT_NotSpliceConsensus */
          str = GetDocText (vep->doc, item, 0, 5);
          if (str != NULL) {
            tmp = str;
            tmp++; /* skip past first newline */
            ch = *tmp;
            while (ch != '\0') {
              if (ch < ' ') {
                *tmp = ' ';
              }
              tmp++;
              ch = *tmp;
            }
          }
          cp = FormatConsensusSpliceReport (str);
          MemFree (str);
          if (cp != NULL) {
            ValNodeAddPointer (&consensus_splice_list, 0, cp);
          }
        } else if (eip->errcode == 5 && (eip->subcode == 124 || eip->subcode == 125 || eip->subcode == 126)) {
          /* ERR_SEQ_FEAT_BadEcNumberValue */
          sfp = SeqMgrGetDesiredFeature (eip->entityID, NULL, eip->itemID, 0, NULL, &fcontext);
          if (sfp != NULL) {
            ValNodeAddPointer (&ecnumber_list, OBJ_SEQFEAT, sfp);
          }
        } else if (eip->errcode == 2 && eip->subcode == 50) {
          if (eip->itemtype == OBJ_SEQDESC) {
            sdp = SeqMgrGetDesiredDescriptor (eip->entityID, NULL, eip->itemID, 0, NULL, &dcontext);
            if (sdp != NULL) {
              ValNodeAddPointer (&specific_host_list, OBJ_SEQDESC, sdp);
            }
          } else if (eip->itemtype == OBJ_SEQFEAT) {
            sfp = SeqMgrGetDesiredFeature (eip->entityID, NULL, eip->itemID, 0, NULL, &fcontext);
            if (sfp != NULL) {
              ValNodeAddPointer (&specific_host_list, OBJ_SEQFEAT, sfp);
            }
          }
        } else if (eip->itemtype == OBJ_SEQFEAT) {
          sfp = SeqMgrGetDesiredFeature (eip->entityID, NULL, eip->itemID, 0, NULL, &fcontext);
          if (sfp != NULL) {
            ValNodeAddPointer (&(report_lists[pos]), OBJ_SEQFEAT, sfp);
          }
        } else if (eip->itemtype == OBJ_SEQDESC) {
          sdp = SeqMgrGetDesiredDescriptor (eip->entityID, NULL, eip->itemID, 0, NULL, &dcontext);
          if (sdp != NULL) {
            ValNodeAddPointer (&(report_lists[pos]), OBJ_SEQDESC, sdp);
          }
        } else if (eip->itemtype == OBJ_BIOSEQ) {
          bsp = GetBioseqGivenIDs (eip->entityID, eip->itemID, OBJ_BIOSEQ);
          if (bsp != NULL) {
            ValNodeAddPointer (&(report_lists[pos]), OBJ_BIOSEQ, bsp);
          }
        }
      }
    }
  }

  if (consensus_splice_list != NULL) {
    fprintf (fp, "Not Splice Consensus\n");
    for (vnp = consensus_splice_list; vnp != NULL; vnp = vnp->next) {
      fprintf (fp, "%s\n", vnp->data.ptrvalue);
    }
    fprintf (fp, "\n");
    found_any = TRUE;
    consensus_splice_list = ValNodeFreeData (consensus_splice_list);
  }

  if (ecnumber_list != NULL) {
    fprintf (fp, "EC Number Errors\n");
    for (vnp = ecnumber_list; vnp != NULL; vnp = vnp->next) {
      cp = GetEcNumberReport (vnp->data.ptrvalue);
      if (cp != NULL) {
        fprintf (fp, "%s", cp);
        found_any = TRUE;
        cp = MemFree (cp);
      }
    }
    fprintf (fp, "\n");
    ecnumber_list = ValNodeFree (ecnumber_list);
  }
 
  if (specific_host_list != NULL) {
    fprintf (fp, "Bad Specific-host Values\n");
    WriteBadSpecificHostTable (specific_host_list, fp);
    found_any = TRUE;
  }

  found_any |= MakeStandardReports (chosen, report_lists, fp);

  for (pos = 0; pos < num_reports; pos++) {
    report_lists[pos] = ValNodeFree (report_lists[pos]);
  }
  report_lists = MemFree (report_lists);
  chosen = ValNodeFree (chosen);      

  FileClose (fp);
  if (found_any) {
    LaunchGeneralTextViewer (path, "Validator Report");
  } else {
    Message (MSG_OK, "No reportable errors");
  }
  FileRemove (path);  
}


static Boolean ValExportProc (ForM f, CharPtr filename)

{
  Char           buf [128];
  FILE           *fp;
  Int2           i, numItems;
  Char           path [PATH_MAX];
  CharPtr        str;
  ValidExtraPtr  vep;

  vep = (ValidExtraPtr) GetObjectExtra (f);
  if (vep != NULL && vep->doc != NULL) {
    path [0] = '\0';
    StringNCpy_0 (path, filename, sizeof (path));
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
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
        if (GetValue (vep->verbose) == 4) {
           GetDocParams (vep->doc, &numItems, NULL);
          for (i = 1; i <= numItems; i++) {
            buf [0] = '\0';
            str = GetDocText (vep->doc, i, 1, 1);
            StringCat (buf, str);
            StringCat (buf, "                ");
            buf [16] = '\0';
            MemFree (str);
            str = GetDocText (vep->doc, i, 1, 2);
            StringCat (buf, str);
            StringCat (buf, "                ");
            buf [32] = '\0';
            MemFree (str);
            str = GetDocText (vep->doc, i, 1, 3);
            StringCat (buf, str);
            MemFree (str);
            fprintf (fp, "%s\n", buf);
          }
        } else {
          SaveDocument (vep->doc, fp);
        }
        FileClose (fp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void ValFindProc (ButtoN b)

{
  Boolean        found;
  Boolean        goOn;
  Int2           numItems;
  Int2           old;
  CharPtr        ptr;
  BaR            sb;
  Int2           start;
  Int4           startsAt;
  Int2           stop;
  Char           str [256];
  CharPtr        tmp;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (b);
  if (vep != NULL) {
    vep->dblClick = FALSE;
    GetDocParams (vep->doc, &numItems, NULL);
    start = vep->selected + 1;
    if (start > numItems) {
      start = 1;
    }
    stop = start;
    GetTitle (vep->searchfor, str, sizeof (str) - 1);
    tmp = str;
    while (*tmp != '\0') {
      if (*tmp < ' ') {
        *tmp = ' ';
      }
      tmp++;
    }
    found = FALSE;
    goOn = TRUE;
    while (goOn) {
      ptr = GetDocText (vep->doc, start, 0, 0);
      if (ptr != NULL) {
        tmp = ptr;
        while (*tmp != '\0') {
          if (*tmp < ' ') {
            *tmp = ' ';
          }
          tmp++;
        }
        if (StringISearch (ptr, str) != NULL) {
          found = TRUE;
          goOn = FALSE;
        }
        ptr = MemFree (ptr);
      }
      if (goOn) {
        start++;
        if (start > numItems) {
          start = 1;
        }
        if (start == stop) {
          goOn = FALSE;
        }
      }
    }
    if (found) {
      old = vep->selected;
      vep->selected = start;
      if (old != start) {
        if (old == 0) {
          InvalBorder (vep->doc, start);
        } else {
          InvalBorder (vep->doc, MIN (start, old));
          InvalBorder (vep->doc, MAX (start, old));
          ValDoNotify (vep, old, FALSE, FALSE);
        }
        if (! Enabled (vep->remove)) {
          Enable (vep->remove);
        }
        Update ();
      }
      ValDoNotify (vep, start, TRUE, FALSE);
      if (ItemIsVisible (vep->doc, start, NULL, NULL, NULL)) {
      } else {
        GetItemParams4 (vep->doc, start, &startsAt, NULL, NULL, NULL, NULL);
        sb = GetSlateVScrollBar ((SlatE) vep->doc);
        SetBarValue (sb, startsAt);
      }
    } else {
      Beep ();
    }
  }
}

static void ValTextProc (TexT t)

{
  ValidExtraPtr  vep;

  vep = GetObjectExtra (t);
  if (vep != NULL) {
    vep->srchtxtlen = TextLength (t);
    if (vep->srchtxtlen > 0) {
      if (! Enabled (vep->find)) {
        Enable (vep->find);
      }
    } else {
      if (Enabled (vep->find)) {
        Disable (vep->find);
      }
    }
  }
}

static void RevalidateProc (ButtoN b)

{
  BaseFormPtr    bfp;
  int            i;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (b);
  if (vep != NULL) {
    for (i = SEV_NONE; i <= SEV_MAX; i++) {
      vep->counts [i] = 0;
    }
    vep->totalcount = 0;
    vep->addedcount = 0;
    vep->remaining = 0;
    vep->clicked = 0;
    vep->selected = 0;
    vep->dblClick = FALSE;
    bfp = vep->bfp;
    if (bfp != NULL && vep->revalProc != NULL) {
      vep->revalProc (bfp->form);
    }
  }
}

static void SetVerbosityAndRepopulate (PopuP p)

{
  /*
  BaseFormPtr    bfp;
  */
  int            i;
  Int2           val;
  ValidExtraPtr  vep;

  vep = GetObjectExtra (p);
  if (vep != NULL) {
    val = GetValue (vep->verbose);
    if (val > 3) {
      doSuppressContext = TRUE;
      doJustShowAccession = TRUE;
    } else if (val > 1) {
      doSuppressContext = TRUE;
      doJustShowAccession = FALSE;
    } else {
      doSuppressContext = FALSE;
      doJustShowAccession = FALSE;
    }
    for (i = SEV_NONE; i <= SEV_MAX; i++) {
      vep->counts [i] = 0;
    }
    vep->totalcount = 0;
    vep->addedcount = 0;
    vep->remaining = 0;
    vep->clicked = 0;
    vep->selected = 0;
    vep->dblClick = FALSE;
    RepopVal (p);
    /*
    bfp = vep->bfp;
    if (bfp != NULL && vep->revalProc != NULL) {
      vep->revalProc (bfp->form);
    }
    */
  }
}

static void CleanupValidProc (GraphiC g, VoidPtr data)

{
  ErrFltrPtr     efp;
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  vep = (ValidExtraPtr) data;
  if (vep != NULL) {
    vep->messages = ValNodeFree (vep->messages);
    vep->messages = FreeErrItemList (vep->messages);
    vep->lastmessage = NULL;
    for (vnp = vep->errorfilter; vnp != NULL; vnp = vnp->next) {
      efp = (ErrFltrPtr) vnp->data.ptrvalue;
      if (efp != NULL) {
        MemFree (efp->text1);
        MemFree (efp->text2);
        MemFree (efp->text3);
      }
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    }
    vep->errorfilter = ValNodeFree (vep->errorfilter);
  }
  StdCleanupFormProc (g, data);
}

static void CopyValToClipboard (ValidExtraPtr vep)

{
  FILE         *fp;
  Char         path [PATH_MAX];
  CharPtr      selected_text;

  if (vep == NULL) return;
  selected_text = GetSelectedDocText (vep->doc, vep->selected_text_start_item, vep->selected_text_start_row,
                                      vep->selected_text_start_col, vep->selected_text_start_offset,
                                      vep->selected_text_end_item, vep->selected_text_end_row,
                                      vep->selected_text_end_col, vep->selected_text_end_offset,
                                      val_columns_for_copy, sizeof (val_columns_for_copy));
  if (StringHasNoText (selected_text)) {      
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      SaveDocument (vep->doc, fp);
      FileClose (fp);
      FileToClipboard (path);
    }
    FileRemove (path);
  } else {
    StringToClipboard (selected_text);
  }
  selected_text = MemFree (selected_text);
}

static void PrintValProc (ValidExtraPtr vep)

{
  if (vep == NULL) return;
  WatchCursor ();
  Update ();
  PrintDocument (vep->doc);
  ArrowCursor ();
  Update ();
}

static void ValFormMessage (ForM f, Int2 mssg)

{
  ValidExtraPtr  vep;

  vep = (ValidExtraPtr) GetObjectExtra (f);
  if (vep != NULL) {
    switch (mssg) {
      case VIB_MSG_EXPORT :
        ValExportProc (f, NULL);
        break;
      case VIB_MSG_PRINT :
        PrintValProc (vep);
        break;
      case VIB_MSG_CLOSE :
        FreeValidateWindow ();
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        CopyValToClipboard (vep);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        ClearValidateWindow ();
        break;
      default :
        if (vep->appmessage != NULL) {
          vep->appmessage (f, mssg);
        }
        break;
    }
  }
}

#ifndef WIN_MAC
extern void CreateStdValidatorFormMenus (WindoW w)

{
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    if (bfp->importform != NULL || bfp->exportform != NULL) {
      if (bfp->importform != NULL) {
        FormCommandItem (m, "Import...", bfp, VIB_MSG_IMPORT);
      }
      if (bfp->exportform != NULL) {
        FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
      }
      SeparatorItem (m);
    }
    FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
    SeparatorItem (m);
    FormCommandItem (m, "Close", bfp, VIB_MSG_CLOSE);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static CharPtr canStillSubmitText =
"Please review and correct the following errors. If you are unable to resolve " \
"an error, submit your sequences, and our annotation staff may contact you " \
"to help resolve the remaining errors when your sequences are processed.";

static CharPtr howToClickText =
"Click on an error item to select and scroll to that feature. " \
"Shift click to choose target sequence and feature if in a set of multiple sequences. " \
"Double click on an error item to launch the appropriate feature editor.";

/* CreateValidateWindowEx is hidden, allowing a revalidate button */
extern WindoW CreateValidateWindowEx (ErrNotifyProc notify, CharPtr title,
                                    FonT font, ErrSev sev, Int2 verbose,
                                    BaseFormPtr bfp, FormActnFunc revalProc,
                                    Boolean okaytosetviewtarget)

{
  ButtoN             b;
  GrouP              c;
  GrouP              f;
  GrouP              g;
  GrouP              h;
  int                i;
  GrouP              ppt;
  GrouP              q;
  StdEditorProcsPtr  sepp;
  GrouP              t;
  GrouP              btn_grp;
  ValidExtraPtr      vep;
  WindoW             w;
  Boolean indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);

  if (validWindow == NULL) {
    vep = (ValidExtraPtr) MemNew (sizeof (ValidExtra));
    if (vep != NULL) {
      vep->indexerVersion = indexerVersion;
      if (title == NULL || title [0] == '\0') {
        title = "Validation Errors";
      }
      w = FixedWindow (-50, -33, -10, -10, title, CloseValidWindow);
      SetObjectExtra (w, vep, CleanupValidProc);
      vep->form = (ForM) w;
      vep->exportform = ValExportProc;
      vep->formmessage = ValFormMessage;

#ifndef WIN_MAC
      CreateStdValidatorFormMenus (w);
#endif

      sepp = (StdEditorProcsPtr) GetAppProperty ("StdValidatorForm");
      if (sepp == NULL) {
        sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
      }
      if (sepp != NULL) {
        SetActivate (w, sepp->activateForm);
        vep->appmessage = sepp->handleMessages;
      }

      if (w != NULL) {
        vep->notify = notify;

        ppt = MultiLinePrompt (w, canStillSubmitText, stdCharWidth * 30, systemFont);

        c = HiddenGroup (w, 8, 0, NULL);
        SetGroupSpacing (c, 10, 2);
        /*
        vep->remove = PushButton (c, "Remove", RemoveProc);
        Disable (vep->remove);
        */
        StaticPrompt (c, "Message", 0, popupMenuHeight, programFont, 'l');
        vep->verbose = PopupList (c, FALSE, SetVerbosityAndRepopulate);
        SetObjectExtra (vep->verbose, vep, NULL);
        PopupItem (vep->verbose, "Verbose");
        PopupItem (vep->verbose, "Normal");
        PopupItem (vep->verbose, "Terse");
        PopupItem (vep->verbose, "Table");
        if (GetAppProperty ("InternalNcbiSequin") != NULL) {
          if (verbose < 1 || verbose > 4) {
            verbose = 1;
          }
          SetValue (vep->verbose, verbose);
          doSuppressContext = FALSE;
          doJustShowAccession = FALSE;
        } else if (verbose > 0) {
          SetValue (vep->verbose, 2);
          doSuppressContext = TRUE;
          doJustShowAccession = FALSE;
        } else {
          SetValue (vep->verbose, 3);
          doSuppressContext = TRUE;
          doJustShowAccession = FALSE;
        }
        StaticPrompt (c, "Severity", 0, popupMenuHeight, programFont, 'l');
        vep->minlevel = PopupList (c, FALSE, RepopVal); /* was SetVerbosityAndRevalidate */
        SetObjectExtra (vep->minlevel, vep, NULL);
        PopupItem (vep->minlevel, "INFO");
        PopupItem (vep->minlevel, "WARN");
        PopupItem (vep->minlevel, "ERROR");
        PopupItem (vep->minlevel, "REJECT");
        if (sev >= 1 && sev <= 4) {
          SetValue (vep->minlevel, (Int2) sev);
        } else {
          SetValue (vep->minlevel, 1);
        }
        vep->revalBtn = PushButton (c, "Revalidate", RevalidateProc);
        SetObjectExtra (vep->revalBtn, vep, NULL);
        Hide (vep->revalBtn);
        /*
        Advance (w);
        */
        Break (w);

        q = HiddenGroup (w, -2, 0, NULL);
        StaticPrompt (q, "Filter", 0, popupMenuHeight, programFont, 'l');
        vep->filter = PopupList (q, FALSE, RepopVal);
        SetObjectExtra (vep->filter, vep, NULL);
        PopupItem (vep->filter, "ALL");
        SetValue (vep->filter, 1);
        Break (w);

        f = HiddenGroup (w, 6, 0, NULL);
        vep->find = PushButton (f, "Find", ValFindProc);
        SetObjectExtra (vep->find, vep, NULL);
        Disable (vep->find);
        vep->searchfor = DialogText (f, "", 15, ValTextProc);
        SetObjectExtra (vep->searchfor, vep, NULL);
        vep->srchtxtlen = 0;
        vep->showncount = StaticPrompt (f, " 0000000 items shown", 0, dialogTextHeight, systemFont, 'l');
        SetTitle (vep->showncount, "");

        g = HiddenGroup (w, 0, -5, NULL);
        h = HiddenGroup (g, 5, 0, NULL);
        SetGroupMargins (h, 11, 0);
        SetGroupSpacing (h, 0, 2);
        /*
        StaticPrompt (h, "Severity", stdCharWidth * 5 - 2, 0, programFont, 'l');
        StaticPrompt (h, "Error", stdCharWidth * 6, 0, programFont, 'l');
        StaticPrompt (h, "Subcode", stdCharWidth * 9, 0, programFont, 'l');
        StaticPrompt (h, "Message", stdCharWidth * 20, 0, programFont, 'l');
        */
        t = MultiLinePrompt (g, howToClickText, stdCharWidth * 30, systemFont);
        vep->doc = DocumentPanel (g, stdCharWidth * 30, stdLineHeight * 20);
        SetObjectExtra (vep->doc, vep, NULL);
        SetDocAutoAdjust (vep->doc, FALSE);
        SetDocProcs (vep->doc, ClickValid, DragValid, ReleaseValid, NULL);
        SetDocShade (vep->doc, DrawValid, NULL, NULL, NULL);
        vep->font = font;

        for (i = SEV_NONE; i <= SEV_MAX; i++) {
          vep->counts [i] = 0;
        }
        vep->totalcount = 0;
        vep->addedcount = 0;
        vep->remaining = 0;
        vep->clicked = 0;
        vep->selected = 0;
        vep->dblClick = FALSE;

        valColFmt [0].pixWidth = stdCharWidth * 5;
        valColFmt [1].pixWidth = stdCharWidth * 6;
        valColFmt [2].pixWidth = stdCharWidth * 19;
        valColFmt [3].pixWidth = stdCharWidth * 20;
        valColFmt [4].pixWidth = stdCharWidth * 30;
        valColFmt [4].pixInset = stdCharWidth * 5 + 5;
        /*
        valColFmt [0].font = systemFont;
        valColFmt [1].font = systemFont;
        valColFmt [2].font = systemFont;
        */

        vep->summary = StaticPrompt (w, "", stdCharWidth * 30, 0, systemFont, 'c');

        btn_grp = HiddenGroup (w, 2, 0, NULL);
        SetGroupSpacing (btn_grp, 10, 10);
        if (indexerVersion) {
          b = PushButton (btn_grp, "Report", MakeValidatorReport);
          SetObjectExtra (b, vep, NULL);
        }

        b = PushButton (btn_grp, "Dismiss", CloseValidButton);

        AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) vep->doc, (HANDLE) t,
                      (HANDLE) vep->summary, (HANDLE) btn_grp, NULL);

        RealizeWindow (w);
        validWindow = w;
      } else {
        MemFree (vep);
      }
    }
  }
  if (validWindow != NULL) {
    vep = (ValidExtraPtr) GetObjectExtra (validWindow);
    if (vep != NULL) {
      vep->bfp = bfp;
      vep->revalProc = revalProc;
      if (vep->revalProc != NULL) {
        SafeShow (vep->revalBtn);
      } else {
        SafeHide (vep->revalBtn);
      }
      vep->okaytosetviewtarget = okaytosetviewtarget;
    }
  }
  return validWindow;
}

extern WindoW CreateValidateWindow (ErrNotifyProc notify, CharPtr title,
                                  FonT font, ErrSev sev, Int2 verbose)

{
  return CreateValidateWindowEx (notify, title, font, sev, verbose, NULL, NULL, FALSE);
}

extern void ShowValidateDoc (void)

{
  Char           str [64];
  ValidExtraPtr  vep;

  if (validWindow == NULL) return;
  vep = GetObjectExtra (validWindow);
  if (vep == NULL) return;
  RepopVal (vep->filter);
  SafeShow (vep->doc);
  SafeShow (vep->searchfor);
  sprintf (str, "INFO: %ld      WARNING: %ld      ERROR: %ld      REJECT: %ld",
           (long) vep->counts [SEV_INFO], (long) vep->counts [SEV_WARNING],
           (long) vep->counts [SEV_ERROR], (long) vep->counts [SEV_REJECT]);
  SafeSetTitle (vep->summary, str);
  if (vep->addedcount > 1) {
    sprintf (str, " %ld items shown", (long) vep->addedcount);
  } else if (vep->addedcount > 0) {
    sprintf (str, " %ld item shown", (long) vep->addedcount);
  } else {
    StringCpy (str, "");
  }
  SafeSetTitle (vep->showncount, str);
  Update ();
  ShowValidateWindow ();
  Update ();
}

extern void HideValidateDoc (void)

{
  ValidExtraPtr  vep;

  if (validWindow == NULL) return;
  vep = GetObjectExtra (validWindow);
  if (vep == NULL) return;
  SafeHide (vep->doc);
  SafeHide (vep->searchfor);
  Update ();
}

extern void ShowValidateWindow (void)

{
  ValidExtraPtr  vep;

  if (validWindow == NULL) {
    CreateValidateWindow (NULL, "Validation Errors",
                          SetSmallFont (), SEV_ERROR, FALSE);
  }
  if (validWindow != NULL) {
    vep = GetObjectExtra (validWindow);
    if (vep != NULL) {
      if (vep->srchtxtlen > 0) {
        SetTitle (vep->searchfor, "");
        vep->srchtxtlen = 0;
      }
      if (Enabled (vep->find)) {
        Disable (vep->find);
      }
    }
    if (! Visible (validWindow)) {
      if (vep->totalcount > 0) {
        Show (validWindow);
      }
    }
    if (Visible (validWindow)) {
      if (Visible (vep->doc)) {
        Select (validWindow);
      }
    }
  }
}

static CharPtr StringMoveAndConvertTabs (CharPtr dst, CharPtr src, CharPtr suffix)

{
  Char     ch;
  CharPtr  ptr;

  ptr = dst;
  if (dst != NULL) {
    if (src != NULL) {
      ptr = StringMove (dst, src);
      ch = *dst;
      while (ch != '\0') {
        if (ch < ' ') {
          *dst = ' ';
        }
        dst++;
        ch = *dst;
      }
    }
    if (suffix != NULL) {
      ptr = StringMove (ptr, suffix);
    }
  }
  return ptr;
}

static CharPtr severityLabel [] = {
  "NONE", "INFO", "WARN", "ERROR", "REJECT", "FATAL", "MAX", NULL
};

extern void RepopulateValidateFilter (void)

{
  ErrFltrPtr     curr;
  RecT           r;
  Char           str [128];
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  ShowValidateWindow ();
  vep = (ValidExtraPtr) GetObjectExtra (validWindow);
  if (vep != NULL) {
    Hide (vep->filter);
    Update ();
    Reset (vep->filter);
    vep->selected = 0;
    PopupItem (vep->filter, "ALL");
    SetValue (vep->filter, 1);
    vnp = vep->errorfilter;
    while (vnp != NULL) {
      curr = (ErrFltrPtr) vnp->data.ptrvalue;
      if (curr != NULL) {
        str [0] = '\0';
        if (! StringHasNoText (curr->text2)) {
          StringNCpy_0 (str, curr->text2, sizeof (str) - 5);
          if (curr->subcode != INT_MIN) {
            if (! StringHasNoText (curr->text3)) {
              StringCat (str, ": ");
              StringCat (str, curr->text3);
            }
          }
        } else {
          StringCpy (str, "??");
        }
        PopupItem (vep->filter, str);
      }
      vnp = vnp->next;
    }
    ObjectRect (vep->filter, &r);
    AdjustPrnt (vep->filter, &r, FALSE);
    Show (vep->filter);
    Update ();
  }
}

static void AppendFilter (ValidExtraPtr vep, CharPtr text1, CharPtr text2,
                          CharPtr text3, ErrSev sev, int errcode, int subcode)

{
  ErrFltrPtr  curr;
  ErrFltrPtr  efp;
  ValNodePtr  last;
  ValNodePtr  tmp;
  ValNodePtr  vnp;

  if (vep != NULL) {
    efp = MemNew (sizeof (ErrFltrData));
    if (efp != NULL) {
      efp->sev = sev;
      efp->errcode = errcode;
      efp->subcode = subcode;
      if (vep->errorfilter == NULL) {
        vnp = ValNodeNew (NULL);
        vep->errorfilter = vnp;
        if (vnp != NULL) {
          vnp->data.ptrvalue = efp;
          efp->text1 = StringSave (text1);
          efp->text2 = StringSave (text2);
          efp->text3 = StringSave (text3);
        }
      } else {
        last = NULL;
        for (vnp = vep->errorfilter; vnp != NULL; vnp = vnp->next) {
          curr = (ErrFltrPtr) vnp->data.ptrvalue;
          if (curr != NULL) {
            if (curr->errcode == efp->errcode &&
                curr->subcode == efp->subcode) {
              MemFree (efp);
              return;
            }
            if (curr->errcode > efp->errcode ||
                (curr->errcode == efp->errcode && curr->subcode > efp->subcode)) {
              tmp = ValNodeNew (NULL);
              if (tmp != NULL) {
                tmp->data.ptrvalue = efp;
                efp->text1 = StringSave (text1);
                efp->text2 = StringSave (text2);
                efp->text3 = StringSave (text3);
                if (last != NULL) {
                  tmp->next = vnp;
                  last->next = tmp;
                } else {
                  tmp->next = vep->errorfilter;
                  vep->errorfilter = tmp;
                }
              }
              return;
            }
          }
          last = vnp;
        }
        tmp = ValNodeNew (vep->errorfilter);
        if (tmp != NULL) {
          tmp->data.ptrvalue = efp;
          efp->text1 = StringSave (text1);
          efp->text2 = StringSave (text2);
          efp->text3 = StringSave (text3);
        }
      }
    }
  }
}

extern int LIBCALLBACK ValidErrHook (const ErrDesc *err)

{
  CharPtr        catname;
  CharPtr        errname;
  CharPtr        expanded;
  ErrItemPtr     eip;
  Int2           minlev;
  ErrSev         msg_level;
  ErrSev         severity;
  Char           str [128];
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  if (err == NULL) return 1;

  /* FileOpen report matches ERR_SEQ_DESCR_FileOpenCollision instead of ERR_SEQ_DESCR_Inconsistent */
  if (err->errcode == E_File && err->subcode == E_FOpen) return 1;
  /* suppress connector messages with 0, 0 as code and subcode - reverted after calling ErrSetLogLevel to suppress these */
  /* if (StringCmp (err->errtext, "[CONN_Read]  Cannot read data (connector \"HTTP\", error \"Closed\")") == 0) return 1; */

  severity = (ErrSev) err->severity;
  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  msg_level = ErrGetMessageLevel ();
  if (severity < msg_level) return 1;

  ShowValidateWindow ();
  if (validWindow == NULL) return 1;
  vep = GetObjectExtra (validWindow);
  if (vep == NULL) return 1;

  (vep->counts [severity])++;
  (vep->totalcount)++;
  if (vep->remaining > 0) {
    (vep->remaining)--;
  }
  if (vep->remaining < 1) {
    sprintf (str, "INFO: %ld      WARNING: %ld      ERROR: %ld      REJECT: %ld",
             (long) vep->counts [SEV_INFO], (long) vep->counts [SEV_WARNING],
             (long) vep->counts [SEV_ERROR], (long) vep->counts [SEV_REJECT]);
    if (vep->totalcount < 30) {
      vep->remaining = 1;
    } else if (vep->totalcount < 100) {
      vep->remaining = 5;
    } else if (vep->totalcount < 300) {
      vep->remaining = 20;
    } else if (vep->totalcount < 1000) {
      vep->remaining = 50;
    } else if (vep->totalcount < 3000) {
      vep->remaining = 100;
    } else if (vep->totalcount < 10000) {
      vep->remaining = 200;
    } else if (vep->totalcount < 30000) {
      vep->remaining = 500;
    } else {
      vep->remaining = 1000;
    }
    SetTitle (vep->summary, str);
    Update ();
  }

  sprintf (str, "%d %d %d %u %u %u", (int) severity, (int) err->errcode, (int) err->subcode,
           (unsigned int) err->entityID, (unsigned int) err->itemID, (unsigned int) err->itemtype);

  catname = GetValidCategoryName (err->errcode);
  errname = GetValidErrorName (err->errcode, err->subcode);
  expanded = GetValidExplanation (err->errcode, err->subcode);

  eip = (ErrItemPtr) MemNew (sizeof (ErrItemData));
  if (eip == NULL) return 1;

  eip->severity = severity;
  eip->errcode = err->errcode;
  eip->subcode = err->subcode;
  eip->entityID = err->entityID;
  eip->itemtype = err->itemtype;
  eip->itemID = err->itemID;
  eip->hidden = StringSaveNoNull (str);
  eip->sevlabel = StringSaveNoNull (severityLabel [severity]);
  eip->catname = StringSaveNoNull (catname);
  eip->errname = StringSaveNoNull (errname);
  eip->accession = NULL;
  eip->message = NULL;
  if (doJustShowAccession) {
    eip->accession = StringSaveNoNull ((CharPtr) err->errtext);
  } else {
    eip->message = StringSaveNoNull ((CharPtr) err->errtext);
  }
  eip->objtype = NULL;
  eip->label = NULL;
  eip->context = NULL;
  eip->location = NULL;
  eip->product = NULL;
  eip->expanded = StringSaveNoNull (expanded);
  eip->userdata = NULL;

    vnp = ValNodeNew (vep->lastmessage);
    if (vep->messages == NULL) {
      vep->messages = vnp;
    }
    vep->lastmessage = vnp;

  if (vnp != NULL) {
    vnp->data.ptrvalue = eip;

    minlev = GetValue (vep->minlevel);
    if (severity >= minlev || severity == SEV_NONE) {
      (vep->addedcount)++;
    }
  }

  AppendFilter (vep, severityLabel [severity], catname, NULL, severity, err->errcode, INT_MIN);
  AppendFilter (vep, severityLabel [severity], catname, errname, severity, err->errcode, err->subcode);

  return 1;
}



extern void LIBCALLBACK ValidErrCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  CharPtr        catname;
  CharPtr        errname;
  CharPtr        expanded;
  ErrItemPtr     eip;
  Int2           minlev;
  ErrSev         msg_level;
  Char           str [128];
  ValidExtraPtr  vep;
  ValNodePtr     vnp;

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  msg_level = ErrGetMessageLevel ();
  if (severity < msg_level) return;

  ShowValidateWindow ();
  if (validWindow == NULL) return;
  vep = GetObjectExtra (validWindow);
  if (vep == NULL) return;

  (vep->counts [severity])++;
  (vep->totalcount)++;
  if (vep->remaining > 0) {
    (vep->remaining)--;
  }
  if (vep->remaining < 1) {
    sprintf (str, "INFO: %ld      WARNING: %ld      ERROR: %ld      REJECT: %ld",
             (long) vep->counts [SEV_INFO], (long) vep->counts [SEV_WARNING],
             (long) vep->counts [SEV_ERROR], (long) vep->counts [SEV_REJECT]);
    if (vep->totalcount < 10) {
      vep->remaining = 1;
    } else if (vep->totalcount < 20) {
      vep->remaining = 10;
    } else if (vep->totalcount < 50) {
      vep->remaining = 25;
    } else if (vep->totalcount < 100) {
      vep->remaining = 50;
    } else if (vep->totalcount < 200) {
      vep->remaining = 100;
    } else if (vep->totalcount < 500) {
      vep->remaining = 250;
    } else if (vep->totalcount < 1000) {
      vep->remaining = 500;
    } else if (vep->totalcount < 2000) {
      vep->remaining = 1000;
    } else if (vep->totalcount < 5000) {
      vep->remaining = 2500;
    } else if (vep->totalcount < 10000) {
      vep->remaining = 5000;
    } else if (vep->totalcount < 20000) {
      vep->remaining = 10000;
    } else {
      vep->remaining = 20000;
    }
    SetTitle (vep->summary, str);
    Update ();
  }

  sprintf (str, "%d %d %d %u %u %u", (int) severity, (int) errcode, (int) subcode,
           (unsigned int) entityID, (unsigned int) itemID, (unsigned int) itemtype);

  catname = GetValidCategoryName (errcode);
  errname = GetValidErrorName (errcode, subcode);
  expanded = GetValidExplanation (errcode, subcode);

  eip = (ErrItemPtr) MemNew (sizeof (ErrItemData));
  if (eip == NULL) return;

  eip->severity = severity;
  eip->errcode = errcode;
  eip->subcode = subcode;
  eip->entityID = entityID;
  eip->itemtype = itemtype;
  eip->itemID = itemID;
  eip->hidden = StringSaveNoNull (str);
  eip->sevlabel = StringSaveNoNull (severityLabel [severity]);
  eip->catname = StringSaveNoNull (catname);
  eip->errname = StringSaveNoNull (errname);
  eip->accession = StringSaveNoNull (accession);
  eip->message = StringSaveNoNull (message);
  eip->objtype = StringSaveNoNull (objtype);
  eip->label = StringSaveNoNull (label);
  eip->context = StringSaveNoNull (context);
  eip->location = StringSaveNoNull (location);
  eip->product = StringSaveNoNull (product);
  eip->expanded = StringSaveNoNull (expanded);
  eip->userdata = userdata;

    vnp = ValNodeNew (vep->lastmessage);
    if (vep->messages == NULL) {
      vep->messages = vnp;
    }
    vep->lastmessage = vnp;

  if (vnp != NULL) {
    vnp->data.ptrvalue = eip;

    minlev = GetValue (vep->minlevel);
    if (severity >= minlev || severity == SEV_NONE) {
      (vep->addedcount)++;
    }
  }

  AppendFilter (vep, severityLabel [severity], catname, NULL, severity, errcode, INT_MIN);
  AppendFilter (vep, severityLabel [severity], catname, errname, severity, errcode, subcode);
}


typedef struct searchextra {
  ButtoN            find;
  TexT              searchfor;
  GrouP             features;
  GrouP             featurebuttons;
  ButtoN            featstouse [86];
  GrouP             descriptors;
  GrouP             descriptorbuttons;
  ButtoN            descstouse [25];
  GrouP             pubs;
  GrouP             pubbuttons;
  ButtoN            pubstouse [5];
  Uint2             entityID;
  SearchGatherProc  gather;
} SearchExtra, PNTR SearchExtraPtr;

typedef struct replaceextra {
  DoC                doc;
  FonT               font;
  ButtoN             find;
  TexT               searchfor;
  TexT               replaceWith;
  Int2               clicked;
  Int2               selected;
  Boolean            dblClick;
  Boolean            verbose;
  ReplaceNotifyProc  notify;
} ReplaceExtra, PNTR ReplaceExtraPtr;

static WindoW  searchWindow = NULL;
static WindoW  replaceWindow = NULL;

static ParData repParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};

static ColData repColFmt [] = {
  {0,  7, 75,  0, NULL, 'l', TRUE,  FALSE, FALSE, FALSE, FALSE}, /* label  */
  {0,  0,  0,  0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}   /* hidden */
};

extern void FreeSearchWindow (void)

{
  if (searchWindow != NULL) {
    Remove (searchWindow);
    searchWindow = NULL;
  }
}

static void CloseSearchWindow (WindoW w)

{
  FreeSearchWindow ();
}

static void CloseSearchButton (ButtoN b)

{
  FreeSearchWindow ();
}

static void SearchTextProc (TexT t)

{
  SearchExtraPtr  sep;

  sep = GetWindowExtra (ParentWindow (t));
  if (sep != NULL) {
    if (TextLength (t) > 0) {
      if (! Enabled (sep->find)) {
        Enable (sep->find);
      }
    } else {
      if (Enabled (sep->find)) {
        Disable (sep->find);
      }
    }
  }
}

/*
static void SearchFeatProc (GrouP g)

{
  SearchExtraPtr  sep;
  Int2            val;

  sep = GetWindowExtra (ParentWindow (g));
  if (sep != NULL) {
    val = GetValue (g);
    switch (val) {
      case 1 :
        Disable (sep->featurebuttons);
        break;
      case 2 :
        Enable (sep->featurebuttons);
        break;
      case 3 :
        Disable (sep->featurebuttons);
        break;
      default :
        break;
    }
  }
}

static void SearchDescProc (GrouP g)

{
  SearchExtraPtr  sep;
  Int2            val;

  sep = GetWindowExtra (ParentWindow (g));
  if (sep != NULL) {
    val = GetValue (g);
    switch (val) {
      case 1 :
        Disable (sep->descriptorbuttons);
        break;
      case 2 :
        Enable (sep->descriptorbuttons);
        break;
      case 3 :
        Disable (sep->descriptorbuttons);
        break;
      default :
        break;
    }
  }
}

static void SearchPubProc (GrouP g)

{
  SearchExtraPtr  sep;
  Int2            val;

  sep = GetWindowExtra (ParentWindow (g));
  if (sep != NULL) {
    val = GetValue (g);
    switch (val) {
      case 1 :
        Disable (sep->pubbuttons);
        break;
      case 2 :
        Enable (sep->pubbuttons);
        break;
      case 3 :
        Disable (sep->pubbuttons);
        break;
      default :
        break;
    }
  }
}
*/

typedef struct searchlist {
  Uint2             subtype;
  Char              thislabel [41];
  CharPtr           searchFor;
  SearchGatherProc  gather;
  ObjMgrPtr         omp;
  SearchExtraPtr    sep;
} SearchList, PNTR SearchPtr;

static StdPrintOptionsPtr searchSpop = NULL;

static Boolean SearchGatherFunc (GatherContextPtr gcp)

{
  ObjMgrTypePtr  omtp;
  SearchPtr      sp;
  Boolean        success;
  CharPtr        text;
  CharPtr        tmp;

  if (gcp == NULL) return TRUE;

  sp = (SearchPtr) gcp->userdata;
  if (sp == NULL) return TRUE;

  if (sp->gather == NULL) return TRUE;

  sp->subtype = 0;
  sp->thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQDESC ||
      gcp->thistype == OBJ_SEQFEAT) {
    omtp = ObjMgrTypeFind (sp->omp, gcp->thistype, NULL, NULL);
    if (omtp == NULL) {
      return TRUE;
    }
    if (omtp->subtypefunc != NULL) {
      sp->subtype = (*(omtp->subtypefunc)) (gcp->thisitem);
    }
    if (omtp->labelfunc != NULL) {
      (*(omtp->labelfunc)) (gcp->thisitem, sp->thislabel, 40, OM_LABEL_BOTH);
    }
  }

  /* filter here */

  if (gcp->thisitem != NULL) {
    success = FALSE;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        success = StdFormatPrint (gcp->thisitem, (AsnWriteFunc) SeqDescAsnWrite,
                                  "StdSeqDesc", searchSpop);
        break;
      case OBJ_SEQFEAT :
        success = StdFormatPrint (gcp->thisitem, (AsnWriteFunc) SeqFeatAsnWrite,
                                  "StdSeqFeat", searchSpop);
        break;
      case OBJ_SEQFEAT_CIT :
        success = StdFormatPrint (gcp->thisitem, (AsnWriteFunc) PubAsnWrite,
                                  "StdPubOnFeat", searchSpop);
        break;
      default :
        return TRUE;
    }
    if (success && searchSpop->ptr != NULL &&
        *((CharPtr) (searchSpop->ptr)) != '\0') {
      text = searchSpop->ptr;
      tmp = text;
      while (*tmp != '\0') {
        if (*tmp < ' ') {
          *tmp = ' ';
        }
        tmp++;
      }
      if (StringISearch (text, sp->searchFor) != NULL) {
        sp->gather (sp->searchFor, searchSpop->ptr, sp->thislabel,
                    gcp->entityID, gcp->itemID, gcp->thistype, sp->subtype);
      }
      searchSpop->ptr = MemFree (searchSpop->ptr);
    }
  }

  return TRUE;
}

static void SearchFindProc (ButtoN b)

{
  GatherScope     gs;
  SearchExtraPtr  sep;
  SearchList      sl;
  Char            str [256];
  CharPtr         tmp;
  WindoW          w;

  w = ParentWindow (b);
  sep = GetWindowExtra (w);
  if (sep != NULL && sep->gather != NULL && sep->entityID != 0) {
    WatchCursor ();
    if (searchSpop == NULL) {
      searchSpop = StdPrintOptionsNew (NULL);
      if (searchSpop != NULL) {
        searchSpop->newline = "\r";
        searchSpop->indent = "";
      } else {
        Message (MSG_ERROR, "StdPrintOptionsNew failed");
        return;
      }
    }
    GetTitle (sep->searchfor, str, sizeof (str) - 1);
    tmp = str;
    while (*tmp != '\0') {
      if (*tmp < ' ') {
        *tmp = ' ';
      }
      tmp++;
    }
    Hide (w);
    sl.omp = ObjMgrGet ();
    sl.searchFor = str;
    sl.gather = sep->gather;
    sl.sep = sep;
    MemSet (&gs, 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQSUB] = FALSE;
    gs.ignore[OBJ_SEQDESC] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    GatherEntity (sep->entityID, (Pointer) &sl, SearchGatherFunc, &gs);
    searchSpop = StdPrintOptionsFree (searchSpop);
    ArrowCursor ();
  }
}

static void CleanupSearchProc (WindoW w, VoidPtr data)

{
  SearchExtraPtr  sep;

  sep = (SearchExtraPtr) data;
  if (sep != NULL) {
  }
  MemFree (data);
}

/*
static CharPtr descriptorNames [] = {
  "", "Molecule Type", "Modifiers", "Method", "Name", "Title",
  "Organism", "Comment", "Numbering", "Map Location", "PIR", 
  "GenBank", "Publication", "Region", "User Object", "SWISSPROT", 
  "Neighbors", "EMBL", "Create Date", "Update Date",
  "PRF",  "PDB", "Heterogen", "BioSource", "Molecule Info", NULL
};

static CharPtr pubNames [] = {
  "", "Pub Descriptor", "Pub Feature", "Cit On Feature", NULL
};
*/

extern void CreateSearchWindow (SearchGatherProc gather, CharPtr title, Uint2 entityID)

{
  ButtoN            b;
  GrouP             f;
  GrouP             g1, g2, g3;
  SearchExtraPtr    sep;
  WindoW            w;

  if (searchWindow == NULL) {
    sep = (SearchExtraPtr) MemNew (sizeof (SearchExtra));
    if (sep != NULL) {
      if (title == NULL || title [0] == '\0') {
        title = "String Search";
      }
      w = FixedWindow (-50, -33, -10, -10, title, CloseSearchWindow);
      if (w != NULL) {
        sep->gather = gather;
        sep->entityID = entityID;

        f = HiddenGroup (w, 4, 0, NULL);
        SetGroupSpacing (f, 10, 2);
        sep->searchfor = DialogText (f, "", 20, SearchTextProc);
        sep->find = DefaultButton (f, "Find", SearchFindProc);
        Disable (sep->find);

        g1 = NULL;
        g2 = NULL;
        g3 = NULL;
        /*
        g1 = NormalGroup (w, -1, 0, "Features", systemFont, NULL);
        SetGroupSpacing (g1, 3, 10);
        sep->features = HiddenGroup (g1, 3, 0, SearchFeatProc);
        RadioButton (sep->features, "All");
        RadioButton (sep->features, "Selected");
        RadioButton (sep->features, "None");
        SetValue (sep->features, 1);
        */
        /*
        sep->featurebuttons = HiddenGroup (g1, 5, 0, NULL);
        i = 1;
        curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
        while (curr != NULL && i < 86) {
          if (key != 0) {
            sep->featstouse [i] = CheckBox (sep->featurebuttons, curr->typelabel, NULL);
            i++;
          }
          curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
        }
        */
        /*
        sep->featurebuttons = HiddenGroup (g1, 2, 0, NULL);
        i = 1;
        fdgp = DispGroupFindNext (NULL,  &key, &label);
        while (fdgp != NULL && i < 86) {
          if (key != 0) {
            sep->featstouse [i] = CheckBox (sep->featurebuttons, fdgp->groupname, NULL);
            i++;
          }
          fdgp = DispGroupFindNext (fdgp,  &key, &label);
        }
        Disable (sep->featurebuttons);
        AlignObjects (ALIGN_CENTER, (HANDLE) sep->features,
                      (HANDLE) sep->featurebuttons, NULL);
        */

        /*
        g2 = NormalGroup (w, -1, 0, "Descriptors", systemFont, NULL);
        SetGroupSpacing (g2, 3, 10);
        sep->descriptors = HiddenGroup (g2, 3, 0, SearchDescProc);
        RadioButton (sep->descriptors, "All");
        RadioButton (sep->descriptors, "Selected");
        RadioButton (sep->descriptors, "None");
        SetValue (sep->descriptors, 1);
        sep->descriptorbuttons = HiddenGroup (g2, 4, 0, NULL);
        for (i = 1; i < 25; i++) {
          sep->descstouse [i] = CheckBox (sep->descriptorbuttons, descriptorNames [i], NULL);
        }
        Disable (sep->descriptorbuttons);
        AlignObjects (ALIGN_CENTER, (HANDLE) sep->descriptors,
                      (HANDLE) sep->descriptorbuttons, NULL);

        g3 = NormalGroup (w, -1, 0, "Publications", systemFont, NULL);
        SetGroupSpacing (g3, 3, 10);
        sep->pubs = HiddenGroup (g3, 3, 0, SearchPubProc);
        RadioButton (sep->pubs, "All");
        RadioButton (sep->pubs, "Selected");
        RadioButton (sep->pubs, "None");
        SetValue (sep->pubs, 1);
        sep->pubbuttons = HiddenGroup (g3, 5, 0, NULL);
        for (i = 1; i < 4; i++) {
          sep->pubstouse [i] = CheckBox (sep->pubbuttons, pubNames [i], NULL);
        }
        Disable (sep->pubbuttons);
        AlignObjects (ALIGN_CENTER, (HANDLE) sep->pubs,
                      (HANDLE) sep->pubbuttons, NULL);
        */

        b = NULL;
#ifdef WIN_MOTIF
        b = PushButton (w, "Dismiss", CloseSearchButton);
#endif

        AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) b, (HANDLE) g1,
                      (HANDLE) g2, (HANDLE) g3, NULL);

        SetWindowExtra (w, sep, CleanupSearchProc);
        RealizeWindow (w);
        searchWindow = w;
      } else {
        MemFree (sep);
      }
    }
  }
}

extern void ShowSearchWindow (void)

{
  SearchExtraPtr  sep;

  if (searchWindow != NULL) {
    sep = GetWindowExtra (searchWindow);
    if (sep != NULL) {
      SetTitle (sep->searchfor, "");
      if (Enabled (sep->find)) {
        Disable (sep->find);
      }
    }
    if (! Visible (searchWindow)) {
      Show (searchWindow);
    }
    Select (searchWindow);
  }
}

static void DrawReplace (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RecT             rct;
  ReplaceExtraPtr  rep;

  rep = GetWindowExtra (ParentWindow (d));
  if (rep != NULL) {
    if (item == rep->selected) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
  }
}

static void ClickReplace (DoC d, PoinT pt)

{
  Int2             item;
  Int2             row;
  ReplaceExtraPtr  rep;

  rep = GetWindowExtra (ParentWindow (d));
  if (rep != NULL) {
    rep->clicked = 0;
    rep->dblClick = dblClick;
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {
      rep->clicked = item;
    }
  }
}

static void RepDoNotify (ReplaceExtraPtr rep, Int2 item, Boolean select)

{
  unsigned int  entityID;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       str;
  unsigned int  subtype;

  if (rep != NULL && rep->doc != NULL && rep->notify != NULL) {
    str = GetDocText (rep->doc, item, 1, 2);
    if (str != NULL &&
        sscanf (str, "%u %u %u %u", &entityID, &itemID, &itemtype, &subtype) == 4) {
      (rep->notify) ((Uint2) entityID, itemID, (Uint2) itemtype,
                     (Uint2) subtype, select, rep->dblClick);
    }
    MemFree (str);
  }
}

static void ReleaseReplace (DoC d, PoinT pt)

{
  Int2             item;
  Int2             old;
  ReplaceExtraPtr  rep;
  Int2             row;

  rep = GetWindowExtra (ParentWindow (d));
  if (rep != NULL) {
    ResetClip ();
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {
      if (item == rep->clicked) {
        old = rep->selected;
        rep->selected = item;
        if (old != item) {
          if (old == 0) {
            InvalBorder (d, item);
          } else {
            InvalBorder (d, MIN (item, old));
            InvalBorder (d, MAX (item, old));
            RepDoNotify (rep, old, FALSE);
          }
          /*
          if (! Enabled (rep->remove)) {
            Enable (rep->remove);
          }
          */
          Update ();
        }
        RepDoNotify (rep, item, TRUE);
      }
    } else if (rep->clicked == 0) {
      if (rep->selected != 0) {
        old = rep->selected;
        rep->selected = 0;
        InvalBorder (d, old);
        RepDoNotify (rep, old, FALSE);
      }
      /*
      if (Enabled (rep->remove)) {
        Disable (rep->remove);
      }
      */
      Update ();
    }
  }
}

extern void FreeReplaceWindow (void)

{
  if (replaceWindow != NULL) {
    Remove (replaceWindow);
    replaceWindow = NULL;
  }
}

static void CloseReplaceWindow (WindoW w)

{
  FreeReplaceWindow ();
}

static void CloseReplaceButton (ButtoN b)

{
  FreeReplaceWindow ();
}

static void CleanupReplaceProc (WindoW w, VoidPtr data)

{
  ReplaceExtraPtr  rep;

  rep = (ReplaceExtraPtr) data;
  if (rep != NULL) {
  }
  MemFree (data);
}

extern void CreateReplaceWindow (ReplaceNotifyProc notify, CharPtr title,
                                 FonT font, Boolean verbose)

{
  ButtoN           b;
  GrouP            g;
  GrouP            h;
  ReplaceExtraPtr  rep;
  WindoW           w;

  if (replaceWindow == NULL) {
    rep = (ReplaceExtraPtr) MemNew (sizeof (ReplaceExtra));
    if (rep != NULL) {
      if (title == NULL || title [0] == '\0') {
        title = "String Replace";
      }
      w = FixedWindow (-50, -33, -10, -10, title, CloseReplaceWindow);
      if (w != NULL) {
        rep->notify = notify;

        g = HiddenGroup (w, 0, -5, NULL);
        h = HiddenGroup (g, 5, 0, NULL);
        SetGroupMargins (h, 11, 0);
        SetGroupSpacing (h, 0, 2);
        rep->doc = DocumentPanel (g, stdCharWidth * 40, stdLineHeight * 20);
        SetDocAutoAdjust (rep->doc, FALSE);
        SetDocProcs (rep->doc, ClickReplace, NULL, ReleaseReplace, NULL);
        SetDocShade (rep->doc, DrawReplace, NULL, NULL, NULL);
        rep->font = font;

        rep->clicked = 0;
        rep->selected = 0;
        rep->dblClick = FALSE;
        rep->verbose = verbose;

        repColFmt [0].pixWidth = stdCharWidth * 40;
        repColFmt [1].pixWidth = stdCharWidth * 40;

        b = NULL;
#ifdef WIN_MOTIF
        b = PushButton (w, "Dismiss", CloseReplaceButton);
#endif

        AlignObjects (ALIGN_CENTER, (HANDLE) rep->doc, (HANDLE) b, NULL);

        SetWindowExtra (w, rep, CleanupReplaceProc);
        RealizeWindow (w);
        replaceWindow = w;
      } else {
        MemFree (rep);
      }
    }
  }
}

extern void LIBCALLBACK AppendReplaceMessage (CharPtr searchFor, CharPtr foundIn, CharPtr label,
                                              Uint2 entityID, Uint4 itemID,
                                              Uint2 itemtype, Uint2 subtype)

{
  CharPtr          buf;
  size_t           len;
  Int2             numItems;
  CharPtr          ptr;
  ReplaceExtraPtr  rep;
  CharPtr          str;

  ShowReplaceWindow ();
  if (replaceWindow != NULL) {
    rep = GetWindowExtra (replaceWindow);
    if (rep != NULL) {
      if (label == NULL) {
        label = "";
      }
      if (foundIn == NULL) {
        foundIn = "";
      }
      len = StringLen (foundIn) + StringLen (label) + 50;
      buf = MemGet (len, MGET_CLEAR);
      if (buf != NULL) {
        if (rep->verbose) {
          str = foundIn;
        } else {
          str = label;
        }
        ptr = StringMoveAndConvertTabs (buf, str, "\t");
        sprintf (ptr, "%u %u %u %u\n", (unsigned int) entityID,
                 (unsigned int) itemID, (unsigned int) itemtype,
                 (unsigned int) subtype);
        AppendText (rep->doc, buf, &repParFmt, repColFmt, rep->font);
        GetDocParams (rep->doc, &numItems, NULL);
        UpdateDocument (rep->doc, numItems, numItems);
        MemFree (buf);
      }
      Select (replaceWindow);
    }
  }
}

extern void LIBCALLBACK StdReplaceNotify (Uint2 entityID, Uint4 itemID, Uint2 itemtype,
                                          Uint2 subtype, Boolean select, Boolean dblClick)

{
  Int2 handled;

  if (dblClick) {
    WatchCursor ();
    handled = GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID,
                                itemtype, 0, 0, itemtype, 0);
    ArrowCursor ();
    if (handled == OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      return;
    }
  }
  ObjMgrSelect (entityID, itemID, itemtype,0,NULL);
}

extern void ShowReplaceWindow (void)

{
  if (replaceWindow == NULL) {
    CreateReplaceWindow (NULL, "String Replace", SetSmallFont (), FALSE);
  }
  if (replaceWindow != NULL) {
    if (! Visible (replaceWindow)) {
      Show (replaceWindow);
    }
    Select (replaceWindow);
  }
}

