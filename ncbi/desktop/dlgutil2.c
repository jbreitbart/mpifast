/*   dlgutil2.c
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
* File Name:  dlgutil2.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.203 $
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

#include <dlogutil.h>
#include <document.h>
#include <gather.h>
#include <subutil.h>
#include <objfdef.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <utilpub.h>
#include <objfeat.h>
#include <objseq.h>
#include <toasn3.h>
#include <explore.h>
#include <findrepl.h>
#ifdef WIN_MOTIF
#include <netscape.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

typedef struct datepage {
  DIALOG_MESSAGE_BLOCK
  TexT          year;
  PopuP         year_popup;
  Int4          start_year;
  Int4          num_years;
  PopuP         month;
  TexT          day;
} DatePage, PNTR DatePagePtr;

extern CharPtr SaveStringFromTextAndStripNewlines (TexT t)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;

  len = TextLength (t);
  if (len > 0) {
    str = MemNew (len + 1);
    if (str != NULL) {
      GetTitle (t, str, len + 1);
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch < ' ') {
          *ptr = ' ';
        }
        ptr++;
        ch = *ptr;
      }
      TrimSpacesAroundString (str);
      if (StringHasNoText (str)) {
        str = MemFree (str);
      }
      return str;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}


extern CharPtr StripNewlines (CharPtr str)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;

  if (str == NULL) return str;
  len = StringLen (str);
  if (len > 0) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch < ' ') {
        *ptr = ' ';
      }
      ptr++;
      ch = *ptr;
    }
    TrimSpacesAroundString (str);
    if (StringHasNoText (str)) {
      str = MemFree (str);
    }
  } else {
    str = MemFree (str);
  }
  return str;
}


extern void NewlinesToTildes (CharPtr str)

{
  Uchar    ch;
  CharPtr  ptr;

  if (StringHasNoText (str)) return;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch < ' ') {
      *ptr = '~';
    }
    ptr++;
    ch = *ptr;
  }
}


static void DatePtrToDatePage (DialoG d, Pointer data)

{
  DatePtr      dp;
  DatePagePtr  dpp;
  Int2         day;
  Char         str [32];
  Int2         year, val;


  dpp = (DatePagePtr) GetObjectExtra (d);
  dp = (DatePtr) data;
  if (dpp != NULL) {
    if (dp == NULL || dp->data[0] != 1) {
      SafeSetValue (dpp->month, 1);
      SafeSetTitle (dpp->day, "");
      SafeSetTitle (dpp->year, "");
      SafeSetValue (dpp->year_popup, 1);
    } else {
      /* set month */
      SetEnumPopup (dpp->month, months_alist, (UIEnum) dp->data [2]);
      /* set day */
      day = (Int2) dp->data [3];
      if (day > 0 && day <= 31) {
        sprintf (str, "%d", (int) day);
        SafeSetTitle (dpp->day, str);
      } else {
        SafeSetTitle (dpp->day, "");
      }
      /* set year */
      year = (Int2) dp->data [1];
      if (year > 0) {
        if (dpp->year_popup == NULL) {
          sprintf (str, "%d", (int) (year + 1900));
          SafeSetTitle (dpp->year, str);
        } else {
          val = year + 1900 - dpp->start_year + 1;
          if (val < 1 || val >= dpp->num_years) {
            sprintf (str, "%d", (int) (year + 1900));
            SafeSetTitle (dpp->year, str);
            dpp->year_popup = NULL;
            Show (dpp->year);
          } else {
            SetValue (dpp->year_popup, val);
          }
        }
      } else {
        if (dpp->year_popup == NULL) {
          SafeSetTitle (dpp->year, "");
        } else {
          SetValue (dpp->year_popup, 1);
        }
      }
    }
  }
}

static Pointer DatePageToDatePtr (DialoG d)

{
  DatePtr      dp;
  DatePagePtr  dpp;
  Int2         day;
  UIEnum       month;
  Char         str [32];
  Int2         year = 0;

  dp = NULL;
  dpp = (DatePagePtr) GetObjectExtra (d);
  if (dpp != NULL) {
    dp = DateNew ();
    if (dp != NULL) {
      dp->data [0] = 1;
      /* get year value */
      if (dpp->year_popup != NULL) {
        year = GetValue (dpp->year_popup);
        if (year > 0) {
          year += dpp->start_year - 1;
        } else {
          dp = DateFree (dp);
          return NULL;
        }
      } else {
        GetTitle (dpp->year, str, sizeof (str));
        if (! StringHasNoText (str)) {
          StrToInt (str, &year);
        }
      }
      if (year >= 1900) {
        dp->data [1] = (Uint1) (year - 1900);
      } else {
        dp = DateFree (dp);
        return dp;
      }
      /* get month value */
      if (GetEnumPopup (dpp->month, months_alist, &month)) {
        dp->data [2] = (Uint1) month;
      } else {
        dp = DateFree (dp);
        return dp;
      }
      /* get day value */
      GetTitle (dpp->day, str, sizeof (str));
      StrToInt (str, &day);
      dp->data [3] = (Uint1) day;
    }
  }
  return (Pointer) dp;
}

extern DialoG CreateDateDialogEx (GrouP prnt, CharPtr title, Int4 start_year, Int4 num_years)

{
  DatePagePtr  dpp;
  GrouP        f;
  GrouP        m;
  GrouP        p;
  GrouP        s;
  GrouP        year_grp;
  Char         year_buf[15];
  Int4         i;

  p = HiddenGroup (prnt, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dpp = (DatePagePtr) MemNew (sizeof (DatePage));
  if (dpp) {

    SetObjectExtra (p, dpp, StdCleanupExtraProc);
    dpp->dialog = (DialoG) p;
    dpp->todialog = DatePtrToDatePage;
    dpp->fromdialog = DatePageToDatePtr;
    dpp->testdialog = NULL;

    dpp->start_year = start_year;
    dpp->num_years = num_years;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    f = HiddenGroup (m, -6, 0, NULL);
    StaticPrompt (f, "Month", 0, popupMenuHeight, programFont, 'l');
    dpp->month = PopupList (f, TRUE, NULL);
    InitEnumPopup (dpp->month, months_alist, NULL);
    SetValue (dpp->month, 1);
    StaticPrompt (f, "Day", 0, dialogTextHeight, programFont, 'l');
    dpp->day = DialogText (f, "", 4, NULL);
    StaticPrompt (f, "Year", 0, dialogTextHeight, programFont, 'l');
    year_grp = HiddenGroup (f, 0, 0, NULL);
    if (start_year > 0 && num_years > -1) {
      dpp->year_popup = PopupList (year_grp, TRUE, NULL);
      for (i = 0; i < num_years; i++) {
        sprintf (year_buf, "%d", start_year + i);
        PopupItem (dpp->year_popup, year_buf);
      }
    }
    dpp->year = DialogText (year_grp, "", 6, NULL);
    if (dpp->year_popup != NULL) {
      SafeHide (dpp->year);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE)dpp->year, (HANDLE)dpp->year_popup, NULL);
  }

  return (DialoG) p;
}


extern DialoG CreateDateDialog (GrouP prnt, CharPtr title)

{
  return CreateDateDialogEx (prnt, title, -1, 0);
}


typedef struct featcit {
  DIALOG_MESSAGE_BLOCK
  DoC         citdoc;
  ValNodePtr  pubset;
} FeatCitPage, PNTR FeatCitPagePtr;

static ParData cofParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData cofColFmt = {0, 0, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData cofeParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData cofeColFmt = {0, 0, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static Uint1 diamondSym [] = {
  0x00, 0x10, 0x38, 0x7C, 0x38, 0x10, 0x00, 0x00
};

static void DrawCitOnFeat (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  Int2  lineHeight;
  RecT  rct;

  if (d != NULL && r != NULL && item > 0 && firstLine == 0) {
    if (item == 1) return; /* citonfeattxt or citonfeathdr explanatory text */
    GetItemParams (d, item, NULL, NULL, NULL, &lineHeight, NULL);
    rct = *r;
    rct.left += 1;
    rct.right = rct.left + 7;
    rct.top += (lineHeight - 7) / 2;
    rct.bottom = rct.top + 7;
    CopyBits (&rct, diamondSym);
    /*
    x = r->left + 1;
    y = r->top + lineHeight / 2;
    MoveTo (x, y);
    LineTo (x + 5, y);
    */
  }
}

static int LIBCALLBACK CompareStrings (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr  str1;
  CharPtr  str2;

  if (ptr1 != NULL && ptr2 != NULL) {
    str1 = *((CharPtr PNTR) ptr1);
    str2 = *((CharPtr PNTR) ptr2);
    if (str1 != NULL && str2 != NULL) {
      return StringICmp (str1, str2);
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static CharPtr citonfeattxt =
"Press 'Edit Citations' to attach publications to this feature. \
Publications must first be added to the record. For biological \
justification, and not to credit the sequencer, create publication \
with the 'Cites a feature on the sequence' scope.\n";

static CharPtr citonfeathdr =
"Press 'Edit Citations' to change publications.\n\n";

static void PubsetPtrToFeatCitPage (DialoG d, Pointer data)

{
  Int2            count;
  FeatCitPagePtr  fpp;
  Int2            i;
  Char            label [128];
  ObjMgrPtr       omp;
  ObjMgrTypePtr   omtp;
  ValNodePtr      ppr;
  ValNodePtr      psp;
  RecT            r;
  CharPtr         PNTR strs;

  fpp = (FeatCitPagePtr) GetObjectExtra (d);
  psp = (ValNodePtr) data;
  if (fpp != NULL) {
    Reset (fpp->citdoc);
    fpp->pubset = PubSetFree (fpp->pubset);
    fpp->pubset = AsnIoMemCopy (data,
                                (AsnReadFunc) PubSetAsnRead,
                                (AsnWriteFunc) PubSetAsnWrite);
    SetDocAutoAdjust (fpp->citdoc, FALSE);
    ObjectRect (fpp->citdoc, &r);
    InsetRect (&r, 4, 4);
    cofColFmt.pixWidth = r.right - r.left;
    cofColFmt.pixInset = 10;
    omp = ObjMgrGet ();
    if (omp == NULL) return;
    omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT_CIT, NULL, NULL);
    if (omtp == NULL || omtp->labelfunc == NULL) return;
    if (psp != NULL && psp->data.ptrvalue != NULL) {
      count = 0;
      for (ppr = psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
        count++;
      }
      if (count > 0) {
        strs = MemNew (sizeof (CharPtr) * (size_t) (count + 1));
        if (strs != NULL) {
          i = 0;
          for (ppr = psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
            (*(omtp->labelfunc)) (ppr, label, 127, OM_LABEL_CONTENT);
            strs [i] = StringSave (label);
            i++;
          }
          HeapSort (strs, count, sizeof (CharPtr), CompareStrings);
          AppendText (fpp->citdoc, citonfeathdr, &cofParFmt, NULL, programFont);
          for (i = 0; i < count; i++) {
            AppendText (fpp->citdoc, strs [i],
                        &cofParFmt, &cofColFmt, programFont);
          }
          for (i = 0; i < count; i++) {
            strs [i] = MemFree (strs [i]);
          }
          MemFree (strs);
        } else {
          AppendText (fpp->citdoc, citonfeattxt, &cofParFmt, NULL, programFont);
        }
      }
    } else {
      AppendText (fpp->citdoc, citonfeattxt, &cofParFmt, NULL, programFont);
    }
    SetDocShade (fpp->citdoc, DrawCitOnFeat, NULL, NULL, NULL);
    UpdateDocument (fpp->citdoc, 0, 0);
  }
}

static Pointer FeatCitPageToPubsetPtr (DialoG d)

{
  FeatCitPagePtr  fpp;
  ValNodePtr      psp;

  psp = NULL;
  fpp = (FeatCitPagePtr) GetObjectExtra (d);
  if (fpp != NULL) {
    psp = AsnIoMemCopy (fpp->pubset,
                        (AsnReadFunc) PubSetAsnRead,
                        (AsnWriteFunc) PubSetAsnWrite);
  }
  return (Pointer) psp;
}

static void CleanupCitOnFeatProc (GraphiC g, VoidPtr data)

{
  FeatCitPagePtr  fpp;

  fpp = (FeatCitPagePtr) data;
  if (fpp != NULL) {
    PubSetFree (fpp->pubset);
  }
  MemFree (data);
}

static DialoG CreateCitOnFeatDialog (GrouP h, CharPtr title)

{
  FeatCitPagePtr  fpp;
  Int2            lineHeight;
  GrouP           m;
  GrouP           p;
  GrouP           s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  fpp = (FeatCitPagePtr) MemNew (sizeof (FeatCitPage));
  if (fpp != NULL) {

    SetObjectExtra (p, fpp, CleanupCitOnFeatProc);
    fpp->dialog = (DialoG) p;
    fpp->todialog = PubsetPtrToFeatCitPage;
    fpp->fromdialog = FeatCitPageToPubsetPtr;
    fpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    SelectFont (programFont);
    lineHeight = LineHeight ();
    SelectFont (systemFont);
    cofParFmt.minHeight = lineHeight + 2;
    StaticPrompt (m, "Citations on Feature", 25 * stdCharWidth, 0, programFont, 'c');
    fpp->citdoc = DocumentPanel (m, 25 * stdCharWidth, 5 * cofParFmt.minHeight);
    SetObjectExtra (fpp->citdoc, fpp, NULL);
    SetDocAutoAdjust (fpp->citdoc, FALSE);
    fpp->pubset = NULL;
  }

  return (DialoG) p;
}

typedef struct citlist {
  Uint2    entityID;
  Uint4    itemID;
  Uint2    itemtype;
  CharPtr  label;
} CitListData, PNTR CitListDataPtr;

typedef struct featcitedit {
  DIALOG_MESSAGE_BLOCK
  DoC             allcitdoc;
  Int2            clickedItem;
  Int2            clickedRow;
  Int2            numitems;
  Int2            lineheight;
  Int2            index;
  BoolPtr         chosen;
  ValNodePtr      citlist;
  ValNodePtr      psp;
  ObjMgrPtr       omp;
  Int2            entityID;
} FeatCitEdit, PNTR FeatCitEditPtr;

static void CleanupFeatCitForm (GraphiC g, VoidPtr data)

{
  FeatCitEditPtr  fcep;

  fcep = (FeatCitEditPtr) data;
  if (fcep != NULL) {
    MemFree (fcep->chosen);
  }
  MemFree (data);
}

static void DrawFeatCit (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  FeatCitEditPtr  fcep;
  RecT            rct;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
    rct.left += 1;
    rct.right = rct.left + fcep->lineheight;
    rct.bottom = rct.top + (rct.right - rct.left);
    FrameRect (&rct);
    if (item > 0 && item <= fcep->numitems) {
      if (fcep->chosen != NULL && fcep->chosen [item - 1]) {
        MoveTo (rct.left, rct.top);
        LineTo (rct.right - 1, rct.bottom - 1);
        MoveTo (rct.left, rct.bottom - 1);
        LineTo (rct.right - 1, rct.top);
      }
    }
  }
}

static void ClickFeatCit (DoC d, PoinT pt)

{
}

static void ReleaseFeatCit (DoC d, PoinT pt)

{
  Int2            col;
  FeatCitEditPtr  fcep;
  Int2            item;
  RecT            rct;
  Int2            row;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep != NULL && fcep->chosen != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, &rct);
    rct.left += 1;
    rct.right = rct.left + fcep->lineheight;
    rct.bottom = rct.top + (rct.right - rct.left);
    if (row == 1 && col == 1 && item > 0 && item <= fcep->numitems && PtInRect (pt, &rct)) {
      if (fcep->chosen [item - 1]) {
        fcep->chosen [item - 1] = FALSE;
      } else {
        fcep->chosen [item - 1] = TRUE;
      }
      InsetRect (&rct, -1, -1);
      InvalRect (&rct);
      Update ();
    }
  }
}

static Boolean GatherAllCits (GatherContextPtr gcp)

{
  CitListDataPtr  cldp;
  FeatCitEditPtr  fcep;
  Char            label [128];
  ObjMgrTypePtr   omtp;
  PubdescPtr      pdp;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;
  fcep = (FeatCitEditPtr) gcp->userdata;
  if (fcep == NULL) return TRUE;
  label [0] = '\0';
  pdp = NULL;
  if (gcp->thistype == OBJ_SEQDESC) {
    sdp = (ValNodePtr) gcp->thisitem;
    if (sdp == NULL || sdp->choice != Seq_descr_pub) return TRUE;
    pdp = sdp->data.ptrvalue;
  } else if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return TRUE;
    pdp = sfp->data.value.ptrvalue;
  } else return TRUE;
  if (pdp == NULL) return TRUE;
  omtp = ObjMgrTypeFind (fcep->omp, gcp->thistype, NULL, NULL);
  if (omtp == NULL || omtp->labelfunc == NULL) return TRUE;
  (*(omtp->labelfunc)) (gcp->thisitem, label, 127, OM_LABEL_CONTENT);
  cldp = (CitListDataPtr) MemNew (sizeof (CitListData));
  if (cldp != NULL) {
    vnp = ValNodeNew (fcep->citlist);
    if (fcep->citlist == NULL) {
      fcep->citlist = vnp;
    }
    if (vnp != NULL) {
      vnp->data.ptrvalue = cldp;
      cldp->entityID = gcp->entityID;
      cldp->itemID = gcp->itemID;
      cldp->itemtype = gcp->thistype;
      cldp->label = StringSave (label);
      (fcep->numitems)++;
    } else {
      MemFree (cldp);
      return TRUE;
    }
  }
  return TRUE;
}

static Boolean PrintCitOnFeatItem (GatherContextPtr gcp)

{
  CharPtr  PNTR  rsultp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Boolean        success;

  rsultp = (CharPtr PNTR) gcp->userdata;
  if (rsultp != NULL) {
    success = FALSE;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          success = StdFormatPrint ((Pointer) sdp, (AsnWriteFunc) SeqDescAsnWrite,
                                    "StdSeqDesc", spop);
        } else {
          *rsultp = StringSave ("Empty Descriptor\n");
          success = TRUE;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          success = StdFormatPrint ((Pointer) sfp, (AsnWriteFunc) SeqFeatAsnWrite,
                                    "StdSeqFeat", spop);
        } else {
          *rsultp = StringSave ("Empty Feature\n");
          success = TRUE;
        }
        break;
      default :
        break;
    }
    if (success) {
      if (spop->ptr != NULL && *((CharPtr) (spop->ptr)) != '\0') {
        *rsultp = spop->ptr;
        spop->ptr = NULL;
      } else {
        *rsultp = StringSave ("Empty Data\n");
      }
    } else {
      *rsultp = StringSave ("Data Failure\n");
    }
  }
  return TRUE;
}

static CharPtr CitOnFeatPrintProc (DoC doc, Int2 item, Pointer data)

{
  unsigned int  entityID;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       rsult;
  CharPtr       str;

  rsult = NULL;
  if (data != NULL) {
    str = (CharPtr) data;
    if (sscanf (str, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
      GatherItem ((Uint2) entityID, (Uint2) itemID, (Uint2) itemtype,
                  (Pointer) &rsult, PrintCitOnFeatItem);
    }
  } else {
    rsult = StringSave ("Null Data\n");
  }
  return rsult;
}

static int LIBCALLBACK CompareCitList (VoidPtr ptr1, VoidPtr ptr2)

{
  CitListDataPtr  cldp1;
  CitListDataPtr  cldp2;
  CharPtr         str1;
  CharPtr         str2;

  if (ptr1 != NULL && ptr2 != NULL) {
    cldp1 = (CitListDataPtr) ptr1;
    cldp2 = (CitListDataPtr) ptr2;
    if (cldp1 != NULL && cldp2 != NULL) {
      str1 = cldp1->label;
      str2 = cldp2->label;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static Boolean MakeMinimalCitOnFeatItem (GatherContextPtr gcp)

{
  PubdescPtr  pdp;
  ValNodePtr  ppr;
  ValNodePtr  sdp;
  SeqFeatPtr  sfp;
  ValNodePtr  vnp;
  ValNodePtr  PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp != NULL) {
    pdp = NULL;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        }
        break;
      default :
        break;
    }
    if (pdp != NULL) {
      ppr = NULL;
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = pdp->pub;
        ppr = MinimizePub (vnp);
        ValNodeFree (vnp);
      }
      vnp = ValNodeNew (*vnpp);
      if (*vnpp == NULL) {
        *vnpp = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = ppr;
      }
    }
  }
  return TRUE;
}

static Boolean MatchMinimalCits (GatherContextPtr gcp)

{
  FeatCitEditPtr  fcep;
  PubdescPtr      pdp;
  ValNodePtr      ppr;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  fcep = (FeatCitEditPtr) gcp->userdata;
  if (fcep != NULL) {
    pdp = NULL;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        }
        break;
      default :
        break;
    }
    if (pdp != NULL && fcep->psp != NULL) {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = pdp->pub;
        for (ppr = fcep->psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
          if (PubLabelMatch (vnp, ppr) == 0) {
            fcep->chosen [fcep->index] = TRUE;
          }
        }
        ValNodeFree (vnp);
      }
    }
  }
  return TRUE;
}

static void CitListToDialog (DialoG d, Pointer userdata)
{
  FeatCitEditPtr  fcep;
  CitListDataPtr  cldp;
  CitListDataPtr  cldpp;
  Int2            count;
  GatherScope     gs;
  Int2            i;
  Int2            j;
  Char            last [128];
  CharPtr         ptr;
  RecT            r;
  Char            str [34];
  ValNodePtr      vnp;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep == NULL)
  {
    return;
  }
  Reset (fcep->allcitdoc);
  
  fcep->citlist = ValNodeFree (fcep->citlist);
  fcep->numitems = 0;
  fcep->chosen = MemFree (fcep->chosen);
  
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_SEQDESC] = FALSE;
  gs.ignore [OBJ_SEQANNOT] = FALSE;
  gs.ignore [OBJ_SEQFEAT] = FALSE;
  gs.seglevels = 1;
  GatherEntity (fcep->entityID, (Pointer) fcep, GatherAllCits, &gs);
  count = fcep->numitems;
  cldpp = (CitListDataPtr) MemNew (sizeof (CitListData) * (size_t) (count + 1));
  if (cldpp != NULL) {
    for (i = 0, vnp = fcep->citlist; i < count && vnp != NULL; i++, vnp = vnp->next) {
      cldp = vnp->data.ptrvalue;
      if (cldp != NULL) {
        cldpp [i] = *cldp;
        cldp->label = NULL;
      }
    }
    fcep->citlist = ValNodeFreeData (fcep->citlist);
    HeapSort (cldpp, count, sizeof (CitListData), CompareCitList);
    last [0] = '\0';
    ObjectRect (fcep->allcitdoc, &r);
    InsetRect (&r, 4, 4);
    cofeColFmt.pixWidth = r.right - r.left;
    cofeColFmt.pixInset = 20;
    fcep->numitems = 0;
    fcep->chosen = MemNew (sizeof (Boolean) * (size_t) (count + 1));

    fcep->psp = (ValNodePtr) userdata;
    for (i = 0, j = 0; i < count; i++) {
      cldp = &(cldpp [i]);
      if (cldp != NULL) {
        if (last == NULL || StringCmp (cldp->label, last) != 0) {
          sprintf (str, "%d %d %d", (int) cldp->entityID,
                   (int) cldp->itemID, (int) cldp->itemtype);
          ptr = StringSave (str);
          (fcep->numitems)++;
          AppendItem (fcep->allcitdoc, CitOnFeatPrintProc, (Pointer) ptr,
                      TRUE, 5, &cofeParFmt, &cofeColFmt, programFont);
          if (fcep->chosen != NULL) {
            fcep->index = j;
            GatherItem ((Uint2) cldp->entityID, (Uint2) cldp->itemID,
                        (Uint2) cldp->itemtype, (Pointer) fcep, MatchMinimalCits);
          }
          j++;
        }
        StringNCpy_0 (last, cldp->label, sizeof (last));
        cldpp [i].label = MemFree (cldpp [i].label);
      }
    }
    fcep->psp = NULL;
  }
}

static Pointer DialogToMinimizedCitList (DialoG d)
{
  FeatCitEditPtr  fcep;
  Pointer         data;
  unsigned int    entityID;
  Int2            i;
  unsigned int    itemID;
  unsigned int    itemtype;
  Int2            numItems;
  ValNodePtr      ppr = NULL, psp = NULL;
  CharPtr         str;
  
  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep == NULL)
  {
    return NULL;
  }
  
  if (fcep->chosen != NULL) {
    GetDocParams (fcep->allcitdoc, &numItems, NULL);
    for (i = 1; i <= numItems; i++) {
      if (fcep->chosen [i - 1]) {
        GetItemParams (fcep->allcitdoc, i, NULL, NULL, NULL, NULL, &data);
        if (data != NULL) {
          str = (CharPtr) data;
          if (sscanf (str, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
            GatherItem ((Uint2) entityID, (Uint2) itemID, (Uint2) itemtype,
                        (Pointer) &ppr, MakeMinimalCitOnFeatItem);
          }
        }
      }
    }
  }
 
  if (ppr != NULL)
  {
    psp = ValNodeNew (NULL);
    if (psp != NULL) {
      psp->choice = 1;
      psp->data.ptrvalue = ppr;
    }
  }
  return psp;
}

extern DialoG FeatCitEditDialog (GrouP parent, Uint2 entityID)
{
  FeatCitEditPtr  fcep;
  GrouP           p;
  
  fcep = (FeatCitEditPtr) MemNew (sizeof (FeatCitEdit));
  if (fcep == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, fcep, CleanupFeatCitForm);

  fcep->dialog = (DialoG) p;
  fcep->todialog = CitListToDialog;
  fcep->fromdialog = DialogToMinimizedCitList;
  fcep->dialogmessage = NULL;
  fcep->testdialog = NULL;

  fcep->entityID = entityID;
  fcep->omp = ObjMgrGet ();
  fcep->numitems = 0;

  fcep->allcitdoc = DocumentPanel (p, 25 * stdCharWidth, 15 * stdLineHeight);
  SetObjectExtra (fcep->allcitdoc, fcep, NULL);
  SetDocProcs (fcep->allcitdoc, ClickFeatCit, NULL, ReleaseFeatCit, NULL);
  SetDocShade (fcep->allcitdoc, DrawFeatCit, NULL, NULL, NULL);
  
  SelectFont (programFont);
  fcep->lineheight = LineHeight ();
  SelectFont (systemFont);

  return (DialoG) p;
}

typedef struct featcitationform
{
  FeatureFormPtr ffp;
  DialoG         citation_list;
   
} FeatCitationFormData, PNTR FeatCitationFormPtr;

static void AcceptFeatCit (ButtoN b)

{
  FeatCitationFormPtr  fcfp;
  FeatureFormPtr       ffp;
  ValNodePtr           psp = NULL;

  fcfp = (FeatCitationFormPtr) GetObjectExtra (b);
  if (fcfp == NULL)
  {
    return;
  }
  
  psp = DialogToPointer (fcfp->citation_list);
  ffp = (FeatureFormPtr) fcfp->ffp;
  if (ffp != NULL) {
    PointerToDialog (ffp->featcits, (Pointer) psp);
    PubSetFree (psp);
  }
  Remove (ParentWindow (b));
}

static void EditFeatCitsProc (ButtoN b)

{
  ButtoN          btn;
  GrouP           c;
  FeatCitationFormPtr  fcfp;
  FeatureFormPtr  ffp;
  WindoW          w;
  ValNodePtr      psp;

  ffp = (FeatureFormPtr) GetObjectExtra (b);
  if (ffp == NULL)
  {
    return;
  }
  fcfp = (FeatCitationFormPtr) MemNew (sizeof (FeatCitationFormData));
  if (fcfp == NULL) 
  {
    return;
  }

  WatchCursor ();
  Update ();
  w = MovableModalWindow (-50, -33, -10, -10, "Citations", NULL);
  SetObjectExtra (w, fcfp, StdCleanupExtraProc);
  
  fcfp->ffp = ffp;
  fcfp->citation_list = FeatCitEditDialog (w, ffp->input_entityID);

  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 3);
  btn = PushButton (c, "Accept", AcceptFeatCit);
  SetObjectExtra (btn, fcfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) fcfp->citation_list, (HANDLE) c, NULL);
  RealizeWindow (w);
  
  psp = DialogToPointer (ffp->featcits);
  PointerToDialog (fcfp->citation_list, psp);
  PubSetFree (psp);

  Show (w);
  Select (w);
  ArrowCursor ();
  Update ();
}

static void ChangeGenePopupOrList (Handle gene)

{
  FeatureFormPtr  ffp;
  Int2            val;

  ffp = (FeatureFormPtr) GetObjectExtra (gene);
  if (ffp != NULL) {
    val = GetValue (ffp->gene);
    if (val == 1) {
      SafeHide (ffp->newGeneGrp);
      SafeHide (ffp->editGeneBtn);
    } else if (val == 2) {
      SafeHide (ffp->editGeneBtn);
      SafeShow (ffp->newGeneGrp);
    } else {
      SafeHide (ffp->newGeneGrp);
      SafeShow (ffp->editGeneBtn);
    }
  }
}

static void ChangeSubGroup (VoidPtr data, Int2 newval, Int2 oldval)

{
  FeatureFormPtr  ffp;

  ffp = (FeatureFormPtr) data;
  if (ffp != NULL) {
    if (ffp->commonPage >= 0 && ffp->commonPage <= 6) {
      SafeHide (ffp->commonSubGrp [ffp->commonPage]);
    }
    ffp->commonPage = newval;
    if (ffp->commonPage >= 0 && ffp->commonPage <= 6) {
      SafeShow (ffp->commonSubGrp [ffp->commonPage]);
    }
  }
}

typedef struct fieldpage {
  DIALOG_MESSAGE_BLOCK
  Int2               numfields;
  CharPtr            PNTR fields;
  TexT               PNTR values;
  ButtoN             PNTR boxes;
} FieldPage, PNTR FieldPagePtr;


static Boolean ShouldSuppressGBQual(Uint1 subtype, CharPtr qual_name)
{
  if (StringHasNoText (qual_name)) {
    return FALSE;
  }

  /* always suppress experiment and inference quals */
  if (StringCmp (qual_name, "experiment") == 0 || StringCmp (qual_name, "inference") == 0) {
    return TRUE;
  }
  
  if (subtype == FEATDEF_ncRNA) {
    if (StringCmp (qual_name, "product") == 0
        || StringCmp (qual_name, "ncRNA_class") == 0) {
      return TRUE;
    }
  } else if (subtype == FEATDEF_tmRNA) {
    if (StringCmp (qual_name, "product") == 0
        || StringCmp (qual_name, "tag_peptide") == 0) {
      return TRUE;
    }
  } else if (subtype == FEATDEF_otherRNA) {
    if (StringCmp (qual_name, "product") == 0) {
      return TRUE;
    }
  }
  return FALSE;
}
    

static Boolean ShouldBeAGBQual (Uint1 subtype, Int2 qual, Boolean allowProductGBQual)

{
  if (qual < 0) return FALSE;
  if (allowProductGBQual && qual == GBQUAL_product) return TRUE;
  if (qual == GBQUAL_citation ||
      qual == GBQUAL_db_xref ||
      qual == GBQUAL_evidence ||
      qual == GBQUAL_exception ||
      qual == GBQUAL_gene ||
      qual == GBQUAL_gene_synonym ||
      qual == GBQUAL_insertion_seq ||
      qual == GBQUAL_label ||
      qual == GBQUAL_locus_tag ||
      qual == GBQUAL_note ||
      qual == GBQUAL_partial ||
      qual == GBQUAL_product ||
      qual == GBQUAL_pseudo ||
      qual == GBQUAL_rpt_unit ||
      qual == GBQUAL_transposon ||
      qual == GBQUAL_experiment ||
      qual == GBQUAL_trans_splicing ||
      qual == GBQUAL_ribosomal_slippage ||
      qual == GBQUAL_inference) {
    return FALSE;
  }
  if (qual == GBQUAL_map && subtype != FEATDEF_ANY && subtype != FEATDEF_repeat_region && subtype != FEATDEF_gap) return FALSE;
  if (qual == GBQUAL_operon && subtype != FEATDEF_ANY && subtype != FEATDEF_operon) return FALSE;
  if (Nlm_GetAppProperty ("SequinUseEMBLFeatures") == NULL) {
    if (qual == GBQUAL_usedin) {
      return FALSE;
    }
  }

  if (qual > -1 && ShouldSuppressGBQual (subtype, ParFlat_GBQual_names [qual].name)) {
    return FALSE;
  }

  return TRUE;
}

static CharPtr TrimParenthesesAndCommasAroundGBString (CharPtr str)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch < ' ' || ch == '(' || ch == ',')) {
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
      if (ch != ')' && ch != ',') {
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

static CharPtr CombineSplitGBQual (CharPtr origval, CharPtr newval)

{
  size_t   len;
  CharPtr  str = NULL;

  if (StringStr (origval, newval) != NULL) return origval;
  len = StringLen (origval) + StringLen (newval) + 5;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return origval;
  TrimParenthesesAndCommasAroundGBString (origval);
  TrimParenthesesAndCommasAroundGBString (newval);
  StringCpy (str, "(");
  StringCat (str, origval);
  StringCat (str, ",");
  StringCat (str, newval);
  StringCat (str, ")");
  /* free original string, knowing return value will replace it */
  MemFree (origval);
  return str;
}

static void CombineTitle (TexT t, CharPtr val)

{
  CharPtr  str;

  if (t == NULL || StringHasNoText (val)) return;
  str = SaveStringFromText (t);
  if (StringDoesHaveText (str)) {
    str = CombineSplitGBQual (str, val);
    SetTitle (t, str);
  } else {
    SetTitle (t, val);
  }
  MemFree (str);
}

static void GBQualPtrToFieldPage (DialoG d, Pointer data)

{
  Char          empty [4];
  FieldPagePtr  fpf;
  Int2          i;
  GBQualPtr     list;

  fpf = (FieldPagePtr) GetObjectExtra (d);
  list = (GBQualPtr) data;
  if (fpf != NULL) {
    while (list != NULL) {
      if (list->qual != NULL && list->val != NULL) {
        for (i = 0; i < fpf->numfields; i++) {
          if (StringICmp (list->qual, fpf->fields [i]) == 0) {
            if (fpf->values [i] != NULL) {
              if (*(list->val) == '\0') {
                empty [0] = '"';
                empty [1] = '"';
                empty [2] = '\0';
                SetTitle (fpf->values [i], empty);
              } else if (StringICmp (list->qual, "rpt_type") == 0 ||
                         StringICmp (list->qual, "rpt_unit") == 0 ||
                         StringICmp (list->qual, "rpt_unit_range") == 0 ||
                         StringICmp (list->qual, "rpt_unit_seq") == 0 ||
                         StringICmp (list->qual, "replace") == 0 ||
                         StringICmp (list->qual, "compare") == 0 ||
                         StringICmp (list->qual, "old_locus_tag") == 0 ||
                         StringICmp (list->qual, "usedin") == 0) {
                CombineTitle (fpf->values [i], list->val);
              } else {
                SetTitle (fpf->values [i], list->val);
              }
            } else if (fpf->boxes [i] != NULL) {
              SetStatus (fpf->boxes [i], TRUE);
            }
          }
        }
      }
      list = list->next;
    }
  }
}

extern void SeqFeatPtrToFieldPage (DialoG d, SeqFeatPtr sfp);
extern void SeqFeatPtrToFieldPage (DialoG d, SeqFeatPtr sfp)

{
  /*
  FieldPagePtr  fpf;
  Int2          i;

  fpf = (FieldPagePtr) GetObjectExtra (d);
  if (fpf != NULL && sfp != NULL) {
        for (i = 0; i < fpf->numfields; i++) {
          if (StringICmp ("exception", fpf->fields [i]) == 0) {
            if (fpf->values [i] != NULL) {
              if (sfp->except_text == NULL || *(sfp->except_text) == '\0') {
                SetTitle (fpf->values [i], "");
              } else {
                SetTitle (fpf->values [i], sfp->except_text);
              }
            }
          } else if (StringICmp ("pseudo", fpf->fields [i]) == 0) {
            if (fpf->boxes [i] != NULL && (! GetStatus (fpf->boxes [i]))) {
              SetStatus (fpf->boxes [i], sfp->pseudo);
            }
          }
        }
  }
  */

}

static Pointer FieldPageToGBQualPtr (DialoG d)

{
  FieldPagePtr  fpf;
  GBQualPtr     gbq;
  GBQualPtr     gbqlast;
  GBQualPtr     head;
  Int2          i;

  head = NULL;
  fpf = (FieldPagePtr) GetObjectExtra (d);
  if (fpf != NULL) {
    gbq = NULL;
    gbqlast = NULL;
    for (i = 0; i < fpf->numfields; i++) {
      if (fpf->fields [i] == NULL ||
          StringHasNoText (fpf->fields [i]) ||
          (fpf->values [i] == NULL && fpf->boxes [i] ==NULL) ||
          (fpf->values [i] != NULL && TextHasNoText (fpf->values [i]))) {
      } else if (fpf->boxes [i] != NULL && (! GetStatus (fpf->boxes [i]))) {
      } else {
        gbq = GBQualNew ();
        if (gbqlast == NULL) {
          head = gbq;
        } else {
          gbqlast->next = gbq;
        }
        gbqlast = gbq;
        if (gbq != NULL) {
          gbq->qual = StringSave (fpf->fields [i]);
          if (fpf->values [i] != NULL) {
            gbq->val = SaveStringFromText (fpf->values [i]);
          } else {
            gbq->val = StringSave ("");
          }
        }
      }
    }
  }
  return (Pointer) head;
}

static void CleanupFieldsPage (GraphiC g, VoidPtr data)

{
  FieldPagePtr  fpf;
  Int2          i;

  fpf = (FieldPagePtr) data;
  if (fpf != NULL) {
    if (fpf->fields != NULL) {
      for (i = 0; i < fpf->numfields; i++) {
        MemFree (fpf->fields [i]);
      }
    }
    MemFree (fpf->fields);
    MemFree (fpf->values);
    MemFree (fpf->boxes);
  }
  MemFree (data);
}

#define LEGAL_FEATURE    1
#define ILLEGAL_FEATURE  2

extern DialoG CreateImportFields (GrouP h, CharPtr name, SeqFeatPtr sfp, Boolean allowProductGBQual)

{
  FieldPagePtr    fpf;
  GrouP           g;
  GBQualPtr       gbq;
  Boolean         hasillegal;
  Int2            i;
  ImpFeatPtr      imp;
  Int2            index;
  Boolean         is_gap = FALSE;
  Int2            j;
  Int2            max;
  Int2            num;
  GrouP           p;
  PrompT          ppt;
  Int2            qual;
  Int2            seen [ParFlat_TOTAL_GBQUAL];
  SematicFeatPtr  sefp;
  Char            str [64];
  Int2            wid;

  p = HiddenGroup (h, 1, 0, NULL);
  fpf = (FieldPagePtr) MemNew (sizeof (FieldPage));
  if (fpf != NULL) {

    SetObjectExtra (p, fpf, CleanupFieldsPage);
    fpf->dialog = (DialoG) p;
    fpf->todialog = GBQualPtrToFieldPage;
    fpf->fromdialog = FieldPageToGBQualPtr;
    fpf->testdialog = NULL;

    if (sfp != NULL && sfp->data.choice == SEQFEAT_IMP) {
      imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (imp != NULL && StringICmp (imp->key, "gap") == 0) {
        is_gap = TRUE;
      }
    }

    gbq = NULL;
    if (sfp != NULL) {
      gbq = sfp->qual;
    }

    j = 0;
    hasillegal = FALSE;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      seen [i] = 0;
    }
    num = 0;

    if (name != NULL) {
      index = GBFeatKeyNameValid (&name, FALSE);
      if (index >= 0) {

        sefp = &(ParFlat_GBFeat [index]);

        for (i = 0; i < sefp->mand_num; i++) {
          qual = sefp->mand_qual [i];
          if (qual > -1 && ShouldBeAGBQual (sfp == NULL ? FEATDEF_ANY : sfp->idx.subtype, qual, allowProductGBQual)) {
            seen [qual] = LEGAL_FEATURE;
          }
        }
        if (StringCmp (name, "repeat_region") == 0) {
          seen [GBQUAL_map] = TRUE;
        }
        if (StringCmp (name, "operon") == 0) {
          seen [GBQUAL_operon] = TRUE;
        }
        for (i = 0; i < sefp->opt_num; i++) {
          qual = sefp->opt_qual [i];
          if (qual > -1 && ShouldBeAGBQual (sfp == NULL ? FEATDEF_ANY : sfp->idx.subtype, qual, allowProductGBQual)) {
            seen [qual] = LEGAL_FEATURE;
          }
        }
      }
    } else if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      /*
      seen [GBQUAL_exception] = LEGAL_FEATURE;
      seen [GBQUAL_pseudo] = LEGAL_FEATURE;
      */
    }

    while (gbq != NULL) {
      qual = GBQualNameValid (gbq->qual);
      if (qual > -1) {
        if (seen [qual] == 0 && qual != GBQUAL_experiment && qual != GBQUAL_inference) {
          seen [qual] = ILLEGAL_FEATURE;
          hasillegal = TRUE;
        }
      }
      gbq = gbq->next;
    }

    SelectFont (programFont);
    max = 0;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] > 0) {
        num++;
        if (seen [i] == LEGAL_FEATURE) {
          StringNCpy_0 (str, ParFlat_GBQual_names [i].name, sizeof (str));
          wid = StringWidth (str) + 2;
        } else {
          str [0] = '\0';
          if (name != NULL) {
            StringCpy (str, "*");
          }
          StringNCat (str, ParFlat_GBQual_names [i].name, sizeof (str) - 2);
          wid = StringWidth (str) + 2;
        }
        if (wid > max) {
          max = wid;
        }
      }
    }
    SelectFont (systemFont);

    fpf->numfields = num;
    fpf->fields = MemNew (sizeof (CharPtr) * (num + 1));
    fpf->values = MemNew (sizeof (TexT) * (num + 1));
    fpf->boxes = MemNew (sizeof (ButtoN) * (num + 1));

    g = HiddenGroup (p, 2, 0, NULL);
    j = 0;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] == LEGAL_FEATURE) {
        fpf->fields [j] = StringSave (ParFlat_GBQual_names [i].name);
        ppt = StaticPrompt (g, fpf->fields [j], max, dialogTextHeight, programFont, 'l');
        if (ParFlat_GBQual_names [i].gbclass != Class_none) {
          fpf->values [j] = DialogText (g, "", 20, NULL);
          /*
          if (i == GBQUAL_estimated_length && is_gap) {
            Disable (fpf->values [j]);
          }
          */
        } else {
          fpf->boxes [j] = CheckBox (g, "", NULL);
          AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) fpf->boxes [j], NULL);
        }
        j++;
      }
    }
    if (hasillegal && name != NULL) {
      StaticPrompt (p, "Illegal Qualifiers", 0, 0, programFont, 'c');
    }
    g = HiddenGroup (p, 2, 0, NULL);
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] == ILLEGAL_FEATURE) {
        fpf->fields [j] = StringSave (ParFlat_GBQual_names [i].name);
        str [0] = '\0';
        if (name != NULL) {
          StringCpy (str, "*");
        }
        StringNCat (str, fpf->fields [j], sizeof (str) - 2);
        ppt = StaticPrompt (g, str, max, dialogTextHeight, programFont, 'l');
        if (ParFlat_GBQual_names [i].gbclass != Class_none) {
          fpf->values [j] = DialogText (g, "", 20, NULL);
        } else {
          fpf->boxes [j] = CheckBox (g, "", NULL);
          AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) fpf->boxes [j], NULL);
        }
        j++;
      }
    }

    if (j == 0) {
      StaticPrompt (p, "See Attributes page to set legal qualifiers for this feature.",
                    0, 0, programFont, 'c');
    }
  }
  return (DialoG) p;
}

static void ChangeCannedMessage (PopuP p)

{
  FeatureFormPtr  ffp;
  Int2            val;

  ffp = (FeatureFormPtr) GetObjectExtra (p);
  if (ffp == NULL) return;
  val = GetValue (p);
  switch (val) {
    case 1 :
      if (Message (MSG_YN, "Clear the explanation field?") == ANS_YES) {
        SetTitle (ffp->exceptText, "");
        if (Message (MSG_YN, "Clear the exception flag?") == ANS_YES) {
          SetStatus (ffp->exception, FALSE);
        }
      }
      break;
    case 2 :
      SetTitle (ffp->exceptText, "RNA editing");
      SetStatus (ffp->exception, TRUE);
      break;
    case 3 :
      SetTitle (ffp->exceptText, "reasons given in citation");
      SetStatus (ffp->exception, TRUE);
      break;
    case 4 :
      SetTitle (ffp->exceptText, "ribosomal slippage");
      SetStatus (ffp->exception, TRUE);
      break;
    case 5 :
      SetTitle (ffp->exceptText, "trans-splicing");
      SetStatus (ffp->exception, TRUE);
      break;
    case 6 :
      SetTitle (ffp->exceptText, "artificial frameshift");
      SetStatus (ffp->exception, TRUE);
      break;
    case 7 :
      SetTitle (ffp->exceptText, "nonconsensus splice site");
      SetStatus (ffp->exception, TRUE);
      break;
    case 8 :
      SetTitle (ffp->exceptText, "rearrangement required for product");
      SetStatus (ffp->exception, TRUE);
      break;
    case 9 :
      SetTitle (ffp->exceptText, "alternative start codon");
      SetStatus (ffp->exception, TRUE);
      break;
    default :
      break;
  }
}

static CharPtr crossRefWarn =
"A gene mapped by cross-reference does not extend the range\n\
of the indicated gene feature.  Overlap is the usual case, and\n\
it does extend the selected gene.";

static CharPtr suppressWarn =
"This will suppress display of a gene qualifier even though\n\
there is a gene feature that overlaps this one.";

static void GeneXrefWarn (GrouP g)

{
  Int2  val;
  Boolean indexerVersion;

  indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
  if (indexerVersion) return;

  val = GetValue (g);
  if (val == 2) {
    Message (MSG_OK, "%s", crossRefWarn);
  } else if (val == 3) {
    Message (MSG_OK, "%s", suppressWarn);
  }
}

static CharPtr  commonRadioFormTabs [] = {
  "General", "Comment", "Citations", "Cross-Refs", "Evidence", "Identifiers", NULL, NULL
};

static CharPtr  commonNoCitFormTabs [] = {
  "General", "Comment", "Cross-Refs", "Evidence", "Identifiers", NULL, NULL
};

static DialoG NewCreateInferenceDialog (GrouP prnt);
extern void Nlm_LaunchGeneFeatEd (ButtoN b);


static CharPtr GetNameForFeature (SeqFeatPtr sfp)
{
  FeatDefPtr curr;
  Uint1      key;
  CharPtr    label = NULL;
  CharPtr    featname = NULL;
  CharPtr    ptr;
  Char       ch;

  if (sfp != NULL) {
    curr = FeatDefFindNext (NULL, &key, &label, sfp->idx.subtype, TRUE);
    if (curr != NULL) {
      featname = StringSave (label);
      ptr = featname;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_UPPER (ch)) {
          *ptr = TO_LOWER (ch);
        }
        ptr++;
        ch = *ptr;
      }
    }
  }
  return featname;
}


extern GrouP CreateCommonFeatureGroupEx (GrouP h, FeatureFormPtr ffp,
                                         SeqFeatPtr sfp, Boolean hasGeneControl,
                                         Boolean hasCitationTab, Boolean hasGeneSuppress)

{
  ButtoN     b;
  GrouP      c;
  PopuP      canned;
  Boolean    cdsQuals;
  GrouP      f;
  GrouP      g;
  GBQualPtr  gbq;
  Boolean    hasQuals;
  Boolean    indexerVersion;
  Char       just;
  GrouP      k;
  GrouP      m;
  GrouP      p;
  Int2       page;
  PrompT     ppt1, ppt2;
  ButtoN     pseudo = NULL;
  GrouP      q;
  GrouP      r;
  GrouP      t;
  GrouP      v;
  GrouP      x;
  GrouP      y;
  CharPtr    featname;

  c = NULL;
  if (ffp != NULL) {
    hasQuals = FALSE;
    cdsQuals = FALSE;
    if (ffp->gbquals == NULL && sfp != NULL && sfp->qual != NULL) {
      for (gbq = sfp->qual; gbq != NULL && !hasQuals; gbq = gbq->next) {
        if (!ShouldSuppressGBQual(sfp->idx.subtype, gbq->qual)) {
          hasQuals = TRUE;
        }
      }
      /*
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          cdsQuals = TRUE;
        }
      }
      */
    }
    m = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);
    ffp->commonPage = 0;
    if (cdsQuals) {
    } else if (hasQuals) {
      commonRadioFormTabs [6] = "Qualifiers";
      commonNoCitFormTabs [5] = "Qualifiers";
    }
    if (hasCitationTab) {
      ffp->commonRadio = CreateFolderTabs (m, commonRadioFormTabs, ffp->commonPage,
                                           0, 0, PROGRAM_FOLDER_TAB,
                                           ChangeSubGroup, (Pointer) ffp);
    } else {
      ffp->commonRadio = CreateFolderTabs (m, commonNoCitFormTabs, ffp->commonPage,
                                           0, 0, PROGRAM_FOLDER_TAB,
                                           ChangeSubGroup, (Pointer) ffp);
    }
    commonRadioFormTabs [6] = NULL;
    commonNoCitFormTabs [5] = NULL;

    p = HiddenGroup (m, 0, 0, NULL);

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    f = HiddenGroup (c, -1, 0, NULL);
    r = HiddenGroup (f, 7, 0, NULL);
    StaticPrompt (r, "Flags", 0, popupMenuHeight, programFont, 'l');
    ffp->partial = CheckBox (r, "Partial", NULL);
    indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
    if (! indexerVersion) {
      Disable (ffp->partial);
    }
    if (ffp->pseudo == NULL) {
      ffp->pseudo = CheckBox (r, "Pseudo", NULL);
      pseudo = ffp->pseudo; /* allows pseudo control on earlier feature-specific page */
    }
    if (indexerVersion) {
      StaticPrompt (r, "Evidence", 0, popupMenuHeight, programFont, 'l');
      ffp->evidence = PopupList (r, TRUE, NULL);
      PopupItem (ffp->evidence, " ");
      PopupItem (ffp->evidence, "Experimental");
      PopupItem (ffp->evidence, "Non-Experimental");
    }
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ffp->partial,
                  (HANDLE) pseudo, (HANDLE) ffp->evidence, NULL);
    r = HiddenGroup (f, -3, 0, NULL);
    ffp->exception = CheckBox (r, "Exception", NULL);
    StaticPrompt (r, "Explanation", 0, dialogTextHeight, programFont, 'l');
    ffp->exceptText = DialogText (r, "", 12, NULL);
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ffp->exception,
                 (HANDLE) ffp->exceptText, NULL);
    if (ffp->this_subtype == FEATDEF_CDS) {
      StaticPrompt (r, "Standard explanation", 0, popupMenuHeight, programFont, 'l');
      canned = PopupList (r, TRUE, ChangeCannedMessage);
      SetObjectExtra (canned, (Pointer) ffp, NULL);
      PopupItem (canned, " ");
      PopupItem (canned, "RNA editing");
      PopupItem (canned, "reasons given in citation");
      PopupItem (canned, "ribosomal slippage");
      PopupItem (canned, "trans-splicing");
      PopupItem (canned, "artificial frameshift");
      PopupItem (canned, "nonconsensus splice site");
      PopupItem (canned, "rearrangement required");
      PopupItem (canned, "alternative start codon");
      if (sfp != NULL && sfp->excpt) {
        if (StringICmp (sfp->except_text, "RNA editing") == 0) {
          SetValue (canned, 2);
        } else if (StringICmp (sfp->except_text, "reasons given in citation") == 0 ||
                   StringICmp (sfp->except_text, "reasons cited in publication") == 0) {
          SetValue (canned, 3);
        } else if (StringICmp (sfp->except_text, "ribosomal slippage") == 0 ||
                   StringICmp (sfp->except_text, "ribosome slippage") == 0) {
          SetValue (canned, 4);
        } else if (StringICmp (sfp->except_text, "trans-splicing") == 0 ||
                   StringICmp (sfp->except_text, "trans splicing") == 0) {
          SetValue (canned, 5);
        } else if (StringICmp (sfp->except_text, "artificial frameshift") == 0) {
          SetValue (canned, 6);
        } else if (StringICmp (sfp->except_text, "non-consensus splice site") == 0 ||
                   StringICmp (sfp->except_text, "nonconsensus splice site") == 0) {
          SetValue (canned, 7);
        } else if (StringICmp (sfp->except_text, "rearrangement required for product") == 0) {
          SetValue (canned, 8);
        } else if (StringICmp (sfp->except_text, "alternative start codon") == 0) {
          SetValue (canned, 9);
        }
      } else {
        SetValue (canned, 1);
      }
    }

    if (cdsQuals) {
      /*
      ffp->gbquals = CreateImportFields (c, NULL, sfp, FALSE);
      */
    }

    g = NULL;
    k = NULL;
    ffp->gene = NULL;
    ffp->genePopup = NULL;
    ffp->geneList = NULL;
    ffp->geneNames = NULL;
    ffp->useGeneXref = NULL;
    ffp->newGeneGrp = NULL;
    ffp->geneSymbol = NULL;
    ffp->geneAllele = NULL;
    ffp->geneDesc = NULL;
    ffp->locusTag = NULL;
    ffp->geneSynonym = NULL;
    for (page = 0; page < 8; page++) {
      ffp->commonSubGrp [page] = NULL;
    }
    page = 0;

    if (hasGeneControl) {
      g = HiddenGroup (c, -2, 0, NULL);
      StaticPrompt (g, "Gene", 0, popupMenuHeight, programFont, 'l');
      v = HiddenGroup (g, 0, 0, NULL);
      ffp->genePopup = PopupList (v, TRUE, (PupActnProc) ChangeGenePopupOrList);
      SetObjectExtra (ffp->genePopup, (Pointer) ffp, NULL);
      PopupItem (ffp->genePopup, " ");
      PopupItem (ffp->genePopup, "New:");
      SetValue (ffp->genePopup, 1);
      Hide (ffp->genePopup);
      ffp->geneList = SingleList (v, 6, 3, (LstActnProc) ChangeGenePopupOrList);
      SetObjectExtra (ffp->geneList, (Pointer) ffp, NULL);
      ListItem (ffp->geneList, " ");
      ListItem (ffp->geneList, "New:");
      SetValue (ffp->geneList, 1);
      Hide (ffp->geneList);
      k = HiddenGroup (c, 3, 0, NULL);
      StaticPrompt (k, "Map by", 0, stdLineHeight, programFont, 'l');
      ffp->useGeneXref = HiddenGroup (k, 3, 0, GeneXrefWarn);
      SetObjectExtra (ffp->useGeneXref, ffp, NULL);
      RadioButton (ffp->useGeneXref, "Overlap");
      RadioButton (ffp->useGeneXref, "Cross-reference");
      RadioButton (ffp->useGeneXref, "Suppress");
      SetValue (ffp->useGeneXref, 1);
      y = HiddenGroup (c, 0, 0, NULL);
      ffp->newGeneGrp = HiddenGroup (y, 2, 0, NULL);
      StaticPrompt (ffp->newGeneGrp, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      ffp->geneSymbol = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Allele", 0, dialogTextHeight, programFont, 'l');
      ffp->geneAllele = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Description", 0, dialogTextHeight, programFont, 'l');
      ffp->geneDesc = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Locus Tag", 0, dialogTextHeight, programFont, 'l');
      ffp->locusTag = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Synonym", 0, dialogTextHeight, programFont, 'l');
      ffp->geneSynonym = DialogText (ffp->newGeneGrp, "", 20, NULL);
      Hide (ffp->newGeneGrp);
      ffp->editGeneBtn = PushButton (y, "Edit Gene Feature", Nlm_LaunchGeneFeatEd);
      SetObjectExtra (ffp->editGeneBtn, ffp, NULL);
      Hide (ffp->editGeneBtn);
    } else if (hasGeneSuppress) {
      k = HiddenGroup (c, 3, 0, NULL);
      StaticPrompt (k, "Map by", 0, stdLineHeight, programFont, 'l');
      ffp->useGeneXref = HiddenGroup (k, 3, 0, GeneXrefWarn);
      SetObjectExtra (ffp->useGeneXref, ffp, NULL);
      RadioButton (ffp->useGeneXref, "Overlap");
      b = RadioButton (ffp->useGeneXref, "Cross-reference");
      Disable (b);
      RadioButton (ffp->useGeneXref, "Suppress");
      SetValue (ffp->useGeneXref, 1);
      y = HiddenGroup (c, 0, 0, NULL);
    }
    ffp->commonSubGrp [page] = c;
    page++;
    AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) g,
                  (HANDLE) k, (HANDLE) ffp->newGeneGrp,
                  (HANDLE) ffp->editGeneBtn, NULL);

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    q = HiddenGroup (c, 0, 2, NULL);
    StaticPrompt (q, "Comment", 25 * stdCharWidth, 0, programFont, 'c');
    if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      ffp->comment = ScrollText (q, 25, 10, programFont, TRUE, NULL);
    } else {
      ffp->comment = ScrollText (q, 25, 5, programFont, TRUE, NULL);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    if (hasCitationTab) {
      c = HiddenGroup (p, -1, 0, NULL);
      SetGroupSpacing (c, 10, 10);
      ffp->featcits = CreateCitOnFeatDialog (c, NULL);
      b = PushButton (c, "Edit Citations", EditFeatCitsProc);
      SetObjectExtra (b, ffp, NULL);
      ffp->commonSubGrp [page] = c;
      AlignObjects (ALIGN_CENTER, (HANDLE) ffp->featcits, (HANDLE) b, NULL);
      Hide (ffp->commonSubGrp [page]);
      page++;
    }

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
      just = 'c';
    } else {
      just = 'l';
      StaticPrompt (c, "This page is read-only", 15 * stdCharWidth, 0, programFont, 'c');
    }
    t = HiddenGroup (c, 2, 0, NULL);
    StaticPrompt (t, "Database", 7 * stdCharWidth, 0, programFont, just);
    StaticPrompt (t, "Object ID", 8 * stdCharWidth, 0, programFont, just);
    ffp->dbxrefs = CreateDbtagDialog (c, 3, -1, 7, 8);
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    q = HiddenGroup (c, 0, -6, NULL);
    ppt1 = StaticPrompt (q, "Experiment", 0, 0, programFont, 'c');
    ffp->experiment = CreateVisibleStringDialog (q, 3, -1, 15);
    ppt2 = StaticPrompt (q, "Inference", 0, 0, programFont, 'c');
    /*
    ffp->inference = CreateInferenceDialog (q, 3, 2, 15);
    */
    ffp->inference = NewCreateInferenceDialog (q);
    AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) ffp->experiment,
                  (HANDLE) ppt2, (HANDLE) ffp->inference, NULL);
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    t = HiddenGroup (c, 2, 0, NULL);
    SetGroupSpacing (t, 10, 30);
    StaticPrompt (t, "Feature ID for this feature", 0, dialogTextHeight, programFont, 'l');
    ffp->featid = DialogText (t, "", 8, NULL);
    StaticPrompt (t, "ID Xref to associated feature", 0, dialogTextHeight, programFont, 'l');
    ffp->fidxref = DialogText (t, "", 8, NULL);
    if (! indexerVersion) {
      Disable (ffp->featid);
      Disable (ffp->fidxref);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;


    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    x = NULL;
    if (hasQuals && ffp->gbquals == NULL) {
      x = HiddenGroup (c, -1, 0, NULL);
      featname = GetNameForFeature (sfp);
      ffp->gbquals = NewCreateImportFields (x, featname, sfp, FALSE);
      featname = MemFree (featname);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;


    AlignObjects (ALIGN_CENTER, (HANDLE) ffp->commonRadio, (HANDLE) ffp->commonSubGrp [0],
                  (HANDLE) ffp->commonSubGrp [1], (HANDLE) ffp->commonSubGrp [2],
                  (HANDLE) ffp->commonSubGrp [3], (HANDLE) ffp->commonSubGrp [4],
                  (HANDLE) ffp->commonSubGrp [5], (HANDLE) ffp->commonSubGrp [6],
                  (HANDLE) ffp->commonSubGrp [7], NULL);
  }
  return c;
}

extern GrouP CreateCommonFeatureGroup (GrouP h, FeatureFormPtr ffp,
                                       SeqFeatPtr sfp, Boolean hasGeneControl,
                                       Boolean hasCitationTab)

{
  return CreateCommonFeatureGroupEx (h, ffp, sfp, hasGeneControl, hasCitationTab, FALSE);
}

static Boolean DlgutilFindBspItem (GatherContextPtr gcp)

{
  BioseqPtr  PNTR bspp;

  bspp = (BioseqPtr PNTR) gcp->userdata;
  if (bspp != NULL && gcp->thistype == OBJ_BIOSEQ) {
    *bspp = (BioseqPtr) gcp->thisitem;
  }
  return TRUE;
}

extern void SetNewFeatureDefaultInterval (FeatureFormPtr ffp)

{
  BioseqPtr     bsp;
  SelStructPtr  sel;
  SeqIntPtr     sip;
  SeqLocPtr     slp;

  if (ffp == NULL) return;
  bsp = NULL;
  GatherItem (ffp->input_entityID, ffp->input_itemID, ffp->input_itemtype,
              (Pointer) (&bsp), DlgutilFindBspItem);
  if (bsp == NULL) return;
  slp = NULL;
  sel = ObjMgrGetSelected ();
  if (sel != NULL && sel->next == NULL && sel->entityID == ffp->input_entityID &&
      sel->itemID == ffp->input_itemID && sel->itemtype == ffp->input_itemtype) {
    if (sel->regiontype == 1 && sel->region != NULL) {
      if (GetBioseqGivenSeqLoc ((SeqLocPtr) sel->region, ffp->input_entityID) == bsp) {
        slp = AsnIoMemCopy (sel->region,
                            (AsnReadFunc) SeqLocAsnRead,
                            (AsnWriteFunc) SeqLocAsnWrite);
      }
    }
  }
  if (slp == NULL) {
    slp = ValNodeNew (NULL);
    if (slp == NULL) return;
    sip = SeqIntNew ();
    if (sip == NULL) {
      slp = SeqLocFree (slp);
      return;
    }
    slp->choice = SEQLOC_INT;
    slp->data.ptrvalue = (Pointer) sip;
    sip->from = 0;
    sip->to = bsp->length - 1;
    sip->strand = Seq_strand_plus;
    sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  }
  if (slp != NULL) {
    PointerToDialog (ffp->location, (Pointer) slp);
    SeqLocFree (slp);
  }
}

extern Boolean FileToScrollText (TexT t, CharPtr path)

{
  FILE     *fp;
  Int8     len;
  Int4     read_len;
  Int4     max;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif
#if (defined(OS_DOS) || defined (OS_NT))
  CharPtr  p;
  CharPtr  q;
#endif

  if (t != NULL && path != NULL && *path != '\0') {
    len = FileLength (path);
    max = (Int4) INT2_MAX;
#ifdef WIN_MOTIF
    max = INT4_MAX;
#endif
#ifdef WIN_MSWIN
    max = INT4_MAX;
#endif
#ifdef WIN_MAC
#ifdef OS_UNIX_DARWIN
    max = INT4_MAX;
#endif
#endif
    if (len > 0 && len < max - 4) {
      str = MemNew (sizeof (char) * ((size_t)len + 3));
      if (str != NULL) {
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          read_len = FileRead (str, sizeof (char), (size_t) len, fp);
          str [ read_len ] = 0;
#if (defined(OS_DOS) || defined (OS_NT))
          p = str;
          q = str;
          while (*p) {
            if (*p == '\r') {
              p++;
            } else {
              *q = *p;
              p++;
              q++;
            }
          }
          *q = '\0';
#endif
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\015';
            }
            p++;
          }
#endif
          FileClose (fp);
          SetTitle (t, str);
        }
        MemFree (str);
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern void ScrollTextToFile (TexT t, CharPtr path)

{
  FILE     *fp;
  size_t   len;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif

  if (t != NULL && path != NULL && *path != '\0') {
    len = TextLength (t);
    if (len > 0) {
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
        str = MemNew (sizeof (char) * (len + 3));
        if (str != NULL) {
          GetTitle (t, str, len + 1);
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\n';
            }
            p++;
          }
#endif
          FileWrite (str, sizeof (char), len, fp);
          MemFree (str);
        }
        FileClose (fp);
      }
    }
  }
}

extern void FileToClipboard (CharPtr path)

{
  FILE     *fp;
  Int8     len;
  Int4     max;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif
#if (defined(OS_DOS) || defined (OS_NT))
  CharPtr  p;
  CharPtr  q;
#endif

  if (path != NULL && *path != '\0') {
    len = FileLength (path);
#ifdef WIN_MOTIF
    max = INT4_MAX;
#else
    max = (Int4) INT2_MAX;
#endif
    if (len > 0 && len < max - 4) {
      str = MemNew (sizeof (char) * ((size_t)len + 3));
      if (str != NULL) {
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          FileRead (str, sizeof (char), (size_t) len, fp);
#if (defined(OS_DOS) || defined (OS_NT))
          p = str;
          q = str;
          while (*p) {
            if (*p == '\r') {
              p++;
            } else {
              *q = *p;
              p++;
              q++;
            }
          }
          *q = '\0';
#endif
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\015';
            }
            p++;
          }
#endif
          FileClose (fp);
          Nlm_StringToClipboard (str);
        }
        MemFree (str);
      }
    }
  }
}

typedef struct textviewform {
  FORM_MESSAGE_BLOCK

  DoC              doc;
  TexT             text;

  TexT             find;
} TextViewForm, PNTR TextViewFormPtr;

static ParData txtParFmt = {FALSE, FALSE, FALSE, FALSE, TRUE, 0, 0};
static ColData txtColFmt = {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

static void ResizeTextViewer (WindoW w)

{
  Int2             height;
  RecT             r;
  RecT             s;
  TextViewFormPtr  tfp;
  Int2             width;

  tfp = (TextViewFormPtr) GetObjectExtra (w);
  if (tfp != NULL) {
    WatchCursor ();
    ObjectRect (w, &r);
    width = r.right - r.left;
    height = r.bottom - r.top;
    if (tfp->doc != NULL) {
      GetPosition (tfp->doc, &s);
      s.right = width - s.left;
      s.bottom = height - s.left;
      SetPosition (tfp->doc, &s);
      AdjustPrnt (tfp->doc, &s, FALSE);
      txtColFmt.pixWidth = screenRect.right - screenRect.left;
      txtColFmt.pixInset = 8;
      if (Visible (tfp->doc) && AllParentsVisible (tfp->doc)) {
        UpdateDocument (tfp->doc, 0, 0);
      }
    }
    if (tfp->text != NULL) {
      GetPosition (tfp->text, &s);
      s.right = width - s.left;
      s.bottom = height - s.left;
      SetPosition (tfp->text, &s);
      AdjustPrnt (tfp->text, &s, FALSE);
    }
    ArrowCursor ();
    Update ();
  }
}

static Boolean ExportTextViewForm (ForM f, CharPtr filename)

{
  FILE             *fp;
  Char             path [PATH_MAX];
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) GetObjectExtra (f);
  if (tfp != NULL && (tfp->doc != NULL || tfp->text != NULL)) {
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
      if (tfp->doc != NULL) {
        fp = FileOpen (path, "w");
        if (fp != NULL) {
          SaveDocument (tfp->doc, fp);
          FileClose (fp);
          return TRUE;
        }
      } else if (tfp->text != NULL) {
        ScrollTextToFile (tfp->text, path);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void LIBCALL PrintTextViewForm (Pointer formDataPtr)

{
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) formDataPtr;
  if (tfp != NULL && tfp->doc != NULL) {
    PrintDocument (tfp->doc);
  }
}

static void TextViewFormMessage (ForM f, Int2 mssg)

{
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) GetObjectExtra (f);
  if (tfp != NULL) {
    switch (mssg) {
      case VIB_MSG_EXPORT :
        ExportTextViewForm (f, NULL);
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_PRINT :
        PrintTextViewForm (tfp);
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
        if (tfp->appmessage != NULL) {
          tfp->appmessage (f, mssg);
        }
        break;
    }
  }
}


static void FindInGeneralText (ButtoN b)

{
  Int2             actual;
  Char             buf [1030];
  Char             ch;
  Int2             cnt;
  Int8             cntr;
  Int2             first;
  FILE             *fp;
  Char             lastch;
  Int4             line = 0;
  ValNodePtr       matches;
  Int2             next;
  Int2             offset = 0;
  Char             path [PATH_MAX];
  CharPtr          ptr;
  Int4             state;
  CharPtr          str;
  TextFsaPtr       tbl;
  TextViewFormPtr  tfp;
  Int4             max;

  tfp = (TextViewFormPtr) GetObjectExtra (b);
  if (tfp == NULL) return;
  if (tfp->doc != NULL) {
    GetOffset (tfp->doc, NULL, &offset);
  } else if (tfp->text != NULL) {
    GetOffset (tfp->text, NULL, &offset);
  }
  first = -1;
  next = -1;
  max = INT2_MAX;

  str = SaveStringFromText (tfp->find);

  if (StringDoesHaveText (str)) {
    TmpNam (path);
    if (ExportForm (tfp->form, path)) {
      tbl = TextFsaNew ();
      if (tbl != NULL) {
        TextFsaAdd (tbl, str);
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          line = 0;
          state = 0;
          cntr = FileLength (path);
          cnt = (Int2) MIN (cntr, 1024L);
          lastch = '\0';
          while (cnt > 0 && cntr > 0 && line <= max && 
                 ((next == -1 && offset > first) || first == -1)) {
            actual = (Int2) FileRead (buf, 1, cnt, fp);
            if (actual > 0) {
              cnt = actual;
              buf [cnt] = '\0';
              ptr = buf;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '\n' || ch == '\r') {
                  if (ch == '\n' && lastch == '\r') {
                    /* do not increment line */
                  } else if (ch == '\r' && lastch == '\n') {
                    /* do not increment line */
                  } else {
                    line++;
                  }
                }
                state = TextFsaNext (tbl, state, ch, &matches);
                if (matches != NULL) {
                  if (first == -1) {
                    first = line;
                  }
                  if (next == -1 && line > offset) {
                    next = line;
                  }
                }
                lastch = ch;
                ptr++;
                ch = *ptr;
              }
              cntr -= cnt;
              cnt = (Int2) MIN (cntr, 1024L);
            } else {
              cnt = 0;
              cntr = 0;
            }
          }
        }
        FileClose (fp);
      }
      TextFsaFree (tbl);
    }
    FileRemove (path);
  }
  MemFree (str);
  if (line > max) {
    Message (MSG_ERROR, "Too many lines for search");
  }

  if (next >= 0) {
    offset = next;
  } else if (first >= 0) {
    offset = first;
  } else return;
  if (tfp->doc != NULL) {
    SetOffset (tfp->doc, 0, offset);
    Update ();
  } else if (tfp->text != NULL) {
    SetOffset (tfp->text, 0, offset);
    Update ();
  }
}

typedef struct repopulateviewer
{
  Nlm_RepopulateViewer   repopulate_func; 
  Pointer                repopulate_data;
  Nlm_RepopulateDataFree free_data_func;
  TextViewFormPtr        tfp;
  FonT                   fnt;
} RepopulateViewerData, PNTR RepopulateViewerPtr;

static void CleanupRepopulateViewer (Nlm_GraphiC g, Nlm_VoidPtr data)
{
  RepopulateViewerPtr rp;
  
  rp = (RepopulateViewerPtr) data;
  if (rp != NULL && rp->free_data_func != NULL)
  {
    (rp->free_data_func)(rp->repopulate_data);
  }
  
  StdCleanupExtraProc (g, data);
}

static void RepopulateViewer (ButtoN b)
{
  RepopulateViewerPtr rp;
  CharPtr             new_path;
  
  rp = (RepopulateViewerPtr) GetObjectExtra (b);
  
  if (rp == NULL || rp->repopulate_func == NULL)
  {
    return;
  }
  
  new_path = (rp->repopulate_func) (rp->repopulate_data);
  if (new_path == NULL)
  {
    return;
  }
  
  if (rp->tfp->text != NULL)
  {
    FileToScrollText (rp->tfp->text, new_path);
  }
  else if (rp->tfp->doc != NULL)
  {
    txtColFmt.pixWidth = screenRect.right - screenRect.left;
    txtColFmt.pixInset = 8;
    DisplayFancy (rp->tfp->doc, new_path, &txtParFmt, &txtColFmt, rp->fnt, 0);
  }
  FileRemove (new_path);
  new_path = MemFree (new_path);
}

static void 
LaunchGeneralTextViewerEx 
(CharPtr path,
 CharPtr title, 
 Boolean useScrollText,
 Nlm_RepopulateViewer   repopulate_func, 
 Pointer                repopulate_data,
 Nlm_RepopulateDataFree free_data_func)

{
  ButtoN              b;
  FonT                fnt;
  GrouP               g;
  Int2                pixheight;
  Int2                pixwidth;
  StdEditorProcsPtr   sepp;
  TextViewFormPtr     tfp;
  TextViewProcsPtr    tvpp;
  RepopulateViewerPtr rp;
  WindoW              w;
#ifndef WIN_MAC
  MenU                m;
#endif

  tfp = (TextViewFormPtr) MemNew (sizeof (TextViewForm));
  if (tfp == NULL) return;

  w = DocumentWindow (-50, -33, -10, -10, title, StdCloseWindowProc, ResizeTextViewer);
  SetObjectExtra (w, tfp, StdCleanupFormProc);
  tfp->form = (ForM) w;
  tfp->exportform = ExportTextViewForm;
  tfp->formmessage = TextViewFormMessage;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    tfp->appmessage = sepp->handleMessages;
  }

  fnt = programFont;
  pixwidth = 35 * stdCharWidth + 17;
  pixheight = 20 * stdLineHeight;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL) {
    pixwidth = MAX (pixwidth, tvpp->minPixelWidth);
    pixheight = MAX (pixheight, tvpp->minPixelHeight);
    if (tvpp->displayFont != NULL) {
      fnt = tvpp->displayFont;
    }
    if (tvpp->activateForm != NULL) {
      SetActivate (w, tvpp->activateForm);
    }
  }

#ifndef WIN_MAC
  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Export...", (BaseFormPtr) tfp, VIB_MSG_EXPORT);
  SeparatorItem (m);
  FormCommandItem (m, "Close", (BaseFormPtr) tfp, VIB_MSG_CLOSE);
  if (tvpp != NULL && tvpp->useScrollText) {
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_DELETE);
  }
#endif

#ifdef WIN_MAC
  if (useScrollText) {
    if (FileLength (path) > 32767) {
      /* text edit window has maximum length on Mac */
      useScrollText = FALSE;
    }
  }
#endif

  /* right now Find button is only in indexer Sequin */
  if (useScrollText) {
    g = HiddenGroup (w, 5, 0, NULL);
    b = PushButton (g, "Find", FindInGeneralText);
    SetObjectExtra (b, tfp, NULL);
    tfp->find = DialogText (g, "", 10, NULL);
    if (repopulate_func != NULL)
    {
      rp = (RepopulateViewerPtr) MemNew (sizeof (RepopulateViewerData));
      if (rp != NULL)
      {
        rp->repopulate_func = repopulate_func;
        rp->repopulate_data = repopulate_data;
        rp->free_data_func = free_data_func;
        rp->tfp = tfp;
        rp->fnt = fnt;
        b = PushButton (g, "Repopulate", RepopulateViewer);
        SetObjectExtra (b, rp, CleanupRepopulateViewer);
      }
    }
  }

  if (useScrollText) {
    tfp->text = ScrollText (w, (pixwidth + stdCharWidth - 1) / stdCharWidth,
                            (pixheight + stdLineHeight - 1) / stdLineHeight,
                            fnt, FALSE, NULL);
    SetObjectExtra (tfp->text, tfp, NULL);
    RealizeWindow (w);
    if (! FileToScrollText (tfp->text, path)) {
      /* SetTitle (tfp->text, "(Text is too large to be displayed in this control.)"); */
      Remove (w);
      LaunchGeneralTextViewerEx (path, title, FALSE,
                                 repopulate_func, repopulate_data, free_data_func);
      return;
    }
  } else {
    tfp->doc = DocumentPanel (w, pixwidth, pixheight);
    SetObjectExtra (tfp->doc, tfp, NULL);
    RealizeWindow (w);
    txtColFmt.pixWidth = screenRect.right - screenRect.left;
    txtColFmt.pixInset = 8;
    DisplayFancy (tfp->doc, path, &txtParFmt, &txtColFmt, fnt, 0);
    /* document.c: SaveTableItem does not strip preceeding tabs if tabCount is 0 */
  }
  Show (w);
  Select (w);
  Update ();
}

extern void LaunchGeneralTextViewerWithRepopulate 
(CharPtr                path,
 CharPtr                title, 
 Nlm_RepopulateViewer   repopulate_func,
 Pointer                repopulate_data,
 Nlm_RepopulateDataFree free_data_func)
{
  TextViewProcsPtr tvpp;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL && tvpp->useScrollText) {
    LaunchGeneralTextViewerEx (path, title, TRUE, 
                               repopulate_func, repopulate_data, free_data_func);
  } else {
    LaunchGeneralTextViewerEx (path, title, FALSE, 
                               repopulate_func, repopulate_data, free_data_func);
  }
}

extern void LaunchGeneralTextViewer (CharPtr path, CharPtr title)

{
  TextViewProcsPtr tvpp;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL && tvpp->useScrollText) {
    LaunchGeneralTextViewerEx (path, title, TRUE, NULL, NULL, NULL);
  } else {
    LaunchGeneralTextViewerEx (path, title, FALSE, NULL, NULL, NULL);
  }
}

extern void LaunchAsnTextViewer (Pointer from, AsnWriteFunc writefunc, CharPtr title)

{
  AsnIoPtr  aip;
  Char      path [PATH_MAX];

  if (from == NULL || writefunc == NULL) return;
  if (StringHasNoText (title)) {
    title = "General ASN.1 Text Viewer";
  }

  TmpNam (path);
  aip = AsnIoOpen (path, "w");
  if (aip != NULL) {
    (*writefunc) (from, aip, NULL);
    AsnIoClose (aip);
    LaunchGeneralTextViewer (path, title);
  }
  FileRemove (path);
}

#ifndef WIN_MAC
extern void CreateStdEditorFormMenus (WindoW w)

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
    FormCommandItem (m, "Close", bfp, VIB_MSG_CLOSE);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static Boolean DlgutilGetLowestStackSeqEntry (GatherContextPtr gcp)

{
  BaseFormPtr  bfp;
  Int2         i;

  if (gcp == NULL) return TRUE;
  bfp = (BaseFormPtr) gcp->userdata;
  if (bfp == NULL) return TRUE;
  if (gcp->gatherstack != NULL && gcp->numstack > 0) {
    for (i = 0; i < gcp->numstack; i++) {
      if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQ ||
          gcp->gatherstack [i].itemtype == OBJ_BIOSEQSET) {
        bfp->input_itemID = gcp->gatherstack [i].itemID;
        bfp->input_itemtype = gcp->gatherstack [i].itemtype;
      }
    }
  }
  return FALSE;
}

extern Boolean SetClosestParentIfDuplicating (BaseFormPtr bfp)

{
  Uint4              itemID;
  Uint2              itemtype;
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (bfp == NULL || sepp == NULL || (! sepp->duplicateExisting)) return FALSE;
  itemID = bfp->input_itemID;
  itemtype = bfp->input_itemtype;
  GatherItem (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype,
              (Pointer) bfp, DlgutilGetLowestStackSeqEntry);
  if (itemID == bfp->input_itemID && itemtype == bfp->input_itemtype) {
    return FALSE;
  }
  return TRUE;
}

/*****************************************************************************
*
*   Bond and Point SeqLoc dialogs
*
*****************************************************************************/

typedef struct pointpage {
  DIALOG_MESSAGE_BLOCK
  TexT               point;
  Int2               count;
  SeqEntryPtr        PNTR bsptr;
  EnumFieldAssoc     PNTR alist;
  PopuP              strand;
  PopuP              seqIdx;
  Boolean            nucsOK;
  Boolean            protsOK;
  Boolean            showIdTags;
} PointPage, PNTR PointPagePtr;

typedef struct bondpage {
  DIALOG_MESSAGE_BLOCK
  DialoG             pointA;
  DialoG             pointB;
} BondPage, PNTR BondPagePtr;

static void FillInPointProducts (SeqEntryPtr sep, Pointer mydata,
                                 Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  PointPagePtr  ppp;

  if (sep != NULL && mydata != NULL && sep->choice == 1) {
    ppp = (PointPagePtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      if ((ppp->nucsOK && ISA_na (bsp->mol)) ||
          (ppp->protsOK && ISA_aa (bsp->mol))) {
        ppp->count++;
        ppp->bsptr [ppp->count] = sep;
      }
    }
  }
}

static ENUM_ALIST(strand_alist)
{" ",             Seq_strand_unknown},  /* 0 */
{"Plus",          Seq_strand_plus},     /* 1 */
{"Minus",         Seq_strand_minus},    /* 2 */
/*
{"Both",          Seq_strand_both},
{"Reverse",       Seq_strand_both_rev},
*/
{"Other",         Seq_strand_other},    /* 255 */
END_ENUM_ALIST

static void SeqPntPtrToPointPage (DialoG d, Pointer data)


{
  BioseqPtr     bsp;
  Int2          j;
  PointPagePtr  ppp;
  SeqEntryPtr   sep;
  Int2          seq;
  SeqPntPtr     spp;
  Char          str [16];
  Uint1         strand;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp == NULL) return;
  spp = (SeqPntPtr) data;
  if (spp != NULL) {
    sprintf (str, "%ld", (long) (spp->point + 1));
    SafeSetTitle (ppp->point, str);
    seq = 0;
    strand = 0;
    bsp = BioseqFind (spp->id);
    if (bsp != NULL) {
      strand = spp->strand;
      if (strand > Seq_strand_both_rev && strand != Seq_strand_other) {
        strand = Seq_strand_unknown;
      }
      if (ppp->bsptr != NULL) {
        for (j = 1; j <= ppp->count && seq == 0; j++) {
          sep = ppp->bsptr [j];
          if (sep != NULL && sep->choice == 1) {
            if (bsp == (BioseqPtr) sep->data.ptrvalue) {
              seq = j;
            }
          }
        }
      }
    }
    SetEnumPopup (ppp->strand, strand_alist, (UIEnum) strand);
    SetEnumPopup (ppp->seqIdx, ppp->alist, (UIEnum) seq);
  } else {
    SafeSetTitle (ppp->point, "");
    SafeSetValue (ppp->strand, 0);
    SafeSetValue (ppp->seqIdx, 0);
  }
}

static Pointer PointPageToSeqPntPtr (DialoG d)

{
  BioseqPtr     bsp;
  PointPagePtr  ppp;
  SeqEntryPtr   sep;
  UIEnum        seq;
  SeqPntPtr     spp;
  Char          str [16];
  UIEnum        strand;
  Int4          val;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp == NULL) return NULL;
  spp = NULL;
  GetTitle (ppp->point, str, sizeof (str) - 1);
  if (StrToLong (str, &val) && val > 0) {
    if (GetEnumPopup (ppp->seqIdx, ppp->alist, &seq) &&
        seq > 0 && seq <= ppp->count) {
      spp = SeqPntNew ();
      if (spp != NULL) {
        spp->point = val - 1;
        if (GetEnumPopup (ppp->strand, strand_alist, &strand)) {
          spp->strand = (Uint1) strand;
        }
        sep = ppp->bsptr [(Int2) seq];
        if (sep != NULL && sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          if (bsp != NULL) {
            spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
          }
        }
      }
    }
  }
  return (Pointer) spp;
}

static void PointEditorMessage (DialoG d, Int2 mssg)

{
  PointPagePtr  ppp;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp != NULL) {
    if (mssg == VIB_MSG_INIT) {
      SeqPntPtrToPointPage (d, NULL);
    } else if (mssg == VIB_MSG_ENTER) {
      Select (ppp->point);
    } else if (mssg == VIB_MSG_RESET) {
    }
  }
}

static void CleanupPointPage (GraphiC g, VoidPtr data)

{
  Int2          j;
  PointPagePtr  ppp;

  ppp = (PointPagePtr) data;
  if (ppp != NULL) {
    MemFree (ppp->bsptr);
    if (ppp->alist != NULL) {
      for (j = 0; j <= ppp->count + 1; j++) {
        MemFree (ppp->alist [j].name);
      }
    }
    MemFree (ppp->alist);
  }
  MemFree (data);
}

static DialoG CreatePointEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep,
                                       Boolean nucsOK, Boolean protsOK)

{
  BioseqPtr     bsp;
  Int4          count;
  GrouP         f;
  Int2          j;
  GrouP         m;
  GrouP         p;
  PointPagePtr  ppp;
  CharPtr       ptr;
  GrouP         s;
  Boolean       showIdTags;
  SeqIdPtr      sip;
  Char          str [128];

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (PointPagePtr) MemNew (sizeof (PointPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, CleanupPointPage);
    ppp->dialog = (DialoG) p;
    ppp->todialog = SeqPntPtrToPointPage;
    ppp->fromdialog = PointPageToSeqPntPtr;
    ppp->dialogmessage = PointEditorMessage;
    ppp->testdialog = NULL;
    ppp->importdialog = NULL;
    ppp->exportdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    ppp->nucsOK = nucsOK;
    ppp->protsOK = protsOK;
    ppp->showIdTags = FALSE;
    ppp->count = 0;

    if (sep != NULL) {
      count = SeqEntryCount (sep);
      count += 4;
      ppp->bsptr = MemNew (sizeof (BioseqPtr) * (size_t) count);
      ppp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) count);
      ppp->count = 0;

      if (ppp->bsptr != NULL && ppp->alist != NULL) {
        SeqEntryExplore (sep, (Pointer) ppp, FillInPointProducts);
        j = 0;
        ppp->alist [j].name = StringSave ("     ");
        ppp->alist [j].value = (UIEnum) 0;
        for (j = 1; j <= ppp->count; j++) {
          sep = ppp->bsptr [j];
          if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            sip = SeqIdFindWorst (bsp->id);
            SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
            ptr = StringChr (str, '|');
            showIdTags = FALSE;
            if (ptr == NULL) {
              ptr = str;
            } else if (showIdTags) {
              ptr = str;
            } else {
              ptr++;
            }
            ppp->alist [j].name = StringSave (ptr);
            ppp->alist [j].value = (UIEnum) j;
          }
        }
        j = ppp->count + 1;
        ppp->alist [j].name = NULL;
        ppp->alist [j].value = (UIEnum) 0;
      }

    } else {
      ppp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) 4);
      if (ppp->alist != NULL) {
        j = 0;
        ppp->alist [j].name = StringSave ("     ");
        ppp->alist [j].value = (UIEnum) 0;
        j = 1;
        ppp->alist [j].name = NULL;
        ppp->alist [j].value = (UIEnum) 0;

      }
    }

    f = HiddenGroup (m, 6, 0, NULL);
    /*StaticPrompt (f, "Point", 0, dialogTextHeight, programFont, 'l');*/
    ppp->point = DialogText (f, "", 5, NULL);
    if (nucsOK) {
      /*StaticPrompt (f, "Strand", 0, popupMenuHeight, programFont, 'c');*/
      ppp->strand = PopupList (f, TRUE, NULL);
      InitEnumPopup (ppp->strand, strand_alist, NULL);
    }
    /*StaticPrompt (f, "SeqID", 0, popupMenuHeight, programFont, 'l');*/
    ppp->seqIdx = PopupList (f, TRUE, NULL);
    InitEnumPopup (ppp->seqIdx, ppp->alist, NULL);
  }

  return (DialoG) p;
}

static void SeqLocPtrToBondPage (DialoG d, Pointer data)


{
  BondPagePtr  bpp;
  SeqBondPtr   sbp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  SeqPnt       sqp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp == NULL) return;
  slp = (SeqLocPtr) data;
  if (slp != NULL) {
    if (slp->choice == SEQLOC_BOND) {
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        PointerToDialog (bpp->pointA, (Pointer) sbp->a);
        PointerToDialog (bpp->pointB, (Pointer) sbp->b);
      }
    } else {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        sqp.strand = SeqLocStrand (slp);
        sqp.id = sip;
        sqp.point = SeqLocStart (slp);
        PointerToDialog (bpp->pointA, (Pointer) &sqp);
        sqp.point = SeqLocStop (slp);
        PointerToDialog (bpp->pointB, (Pointer) &sqp);
      }
    }
  }
}

static Pointer BondPageToSeqLocPtr (DialoG d)

{
  BondPagePtr  bpp;
  SeqBondPtr   sbp;
  SeqLocPtr    slp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp == NULL) return NULL;
  slp = NULL;
  sbp = SeqBondNew ();
  if (sbp != NULL) {
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      slp->choice = SEQLOC_BOND;
      slp->data.ptrvalue = (Pointer) sbp;
      sbp->a = DialogToPointer (bpp->pointA);
      sbp->b = DialogToPointer (bpp->pointB);
    } else {
      SeqBondFree (sbp);
    }
  }
  return (Pointer) slp;
}

static void BondEditorMessage (DialoG d, Int2 mssg)

{
  BondPagePtr  bpp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp != NULL) {
    if (mssg == VIB_MSG_INIT) {
      SeqLocPtrToBondPage (d, NULL);
    } else if (mssg == VIB_MSG_ENTER) {
      SendMessageToDialog (bpp->pointA, VIB_MSG_ENTER);
    } else if (mssg == VIB_MSG_RESET) {
    }
  }
}

extern DialoG CreateBondEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep);

extern DialoG CreateBondEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep)

{
  BondPagePtr  bpp;
  GrouP        f;
  GrouP        m;
  GrouP        p;
  GrouP        s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  bpp = (BondPagePtr) MemNew (sizeof (BondPage));
  if (bpp != NULL) {

    SetObjectExtra (p, bpp, StdCleanupExtraProc);
    bpp->dialog = (DialoG) p;
    bpp->todialog = SeqLocPtrToBondPage;
    bpp->fromdialog = BondPageToSeqLocPtr;
    bpp->dialogmessage = BondEditorMessage;
    bpp->testdialog = NULL;
    bpp->importdialog = NULL;
    bpp->exportdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    f = HiddenGroup (m, 2, 0, NULL);
    StaticPrompt (f, "From", 0, popupMenuHeight, programFont, 'l');
    bpp->pointA = CreatePointEditorDialog (f, NULL, sep, FALSE, TRUE);
    StaticPrompt (f, "(To)", 0, popupMenuHeight, programFont, 'l');
    bpp->pointB = CreatePointEditorDialog (f, NULL, sep, FALSE, TRUE);
  }

  return (DialoG) p;
}

void GetRidOfEmptyFeatsDescStrings (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, GetRidOfEmptyFeatsDescCallback);
}

extern Int2 LIBCALLBACK StdVibrantEditorMsgFunc (OMMsgStructPtr ommsp)

{
  BaseFormPtr    bfp;
  OMUserDataPtr  omudp;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  bfp = (BaseFormPtr) omudp->userdata.ptrvalue;
  if (bfp == NULL) return OM_MSG_RET_ERROR;
  switch (ommsp->message) {
    case OM_MSG_DEL:
      Remove (bfp->form);
      Update ();
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}

/* launch url section */

NLM_EXTERN void LaunchEntrezURL (CharPtr database, Int4 uid, CharPtr format)

{
#ifdef WIN_MOTIF
  NS_Window  window = NULL;
#endif

  Char  url [256];

  if (uid < 1 || StringHasNoText (database) || StringHasNoText (format)) return;
  sprintf (url,
           "http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=%s&list_uids=%ld&dopt=%s",
            database, (long) uid, format);

#ifdef WIN_MAC
  Nlm_SendURLAppleEvent (url, "MOSS", NULL);
#endif
#ifdef WIN_MSWIN
  if (! Nlm_MSWin_OpenDocument (url)) {
    Message (MSG_POST, "Unable to launch browser");
  }
#endif
#ifdef WIN_MOTIF
  if (! NS_OpenURL (&window, url, NULL, TRUE)) {
    Message (MSG_POST, "Unable to launch browser");
  }
  NS_WindowFree (window);
#endif
}

extern void ModalAcceptButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->accepted = TRUE;
  }
}

extern void ModalCancelButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->cancelled = TRUE;
  }
}

extern void ModalThirdOptionButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->third_option = TRUE;
  }
}

typedef struct tabledisplay 
{
  DIALOG_MESSAGE_BLOCK
  PaneL panel;
  ValNodePtr row_list;
  Int4 frozen_header;
  Int4 frozen_left;
  Int4 table_inset;
  Int4 char_width;
  Int4 descent;
  FonT display_font;
  TableDisplayDblClick dbl_click;
  Pointer dbl_click_data;
  TableDisplayLeftInRed left_in_red;
  Pointer left_in_red_data;
} TableDisplayData, PNTR TableDisplayPtr;

extern ValNodePtr FreeTableDisplayRowList (ValNodePtr row_list)
{
  ValNodePtr row_vnp, column_list;
  
  if (row_list != NULL)
  {
    /* free table text */
    for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
    {
      column_list = (ValNodePtr) row_vnp->data.ptrvalue;
      row_vnp->data.ptrvalue = ValNodeFreeData (column_list);
    }
    row_list = ValNodeFree (row_list);
  }
  return row_list;
}

extern void PrintTableDisplayRowListToFile (ValNodePtr row_list, FILE *fp)
{
  ValNodePtr row_vnp, col_vnp, column_list;
  CharPtr    txt_val;
  
  if (row_list == NULL || fp == NULL)
  {
    return;
  }
  
  for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    column_list = (ValNodePtr) row_vnp->data.ptrvalue;
    for (col_vnp = column_list; col_vnp != NULL; col_vnp = col_vnp->next)
    {
      txt_val = (CharPtr) col_vnp->data.ptrvalue;
      if (!StringHasNoText (txt_val))
      {
        fprintf (fp, "%s", txt_val);
      }
      if (col_vnp->next == NULL)
      {
        fprintf (fp, "\n");
      }
      else
      {
        fprintf (fp, "\t");
      }
    }
  }
}

static ValNodePtr ValNodeStringListCopy (ValNodePtr orig_list)
{
  ValNodePtr new_list = NULL;
  
  if (orig_list == NULL)
  {
    return NULL;
  }
  
  new_list = ValNodeNew (NULL);
  new_list->choice = orig_list->choice;
  new_list->data.ptrvalue = StringSave (orig_list->data.ptrvalue);
  new_list->next = ValNodeStringListCopy (orig_list->next);
  return new_list;
}

extern ValNodePtr CopyTableDisplayRowList (ValNodePtr row_list)
{
  ValNodePtr new_row_list = NULL;
  
  if (row_list == NULL)
  {
    return NULL; 
  }
  
  new_row_list = ValNodeNew (NULL);
  new_row_list->choice = row_list->choice;
  new_row_list->data.ptrvalue = ValNodeStringListCopy (row_list->data.ptrvalue);
  new_row_list->next = CopyTableDisplayRowList (row_list->next);
  return new_row_list;
}

static void CleanupTableDisplayDialog (GraphiC g, VoidPtr data)
{
  TableDisplayPtr dlg;

  dlg = (TableDisplayPtr) data;
  if (dlg != NULL) {
    dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
  }
  StdCleanupExtraProc (g, data);
}

static void UpdateTableDisplayDialogScrollBars (TableDisplayPtr dlg)
{
  BaR  sb_vert;
  BaR  sb_horiz;
  Int4 start_row, start_col;
  Int4 num_rows, num_columns, visible_rows;
  Int4 new_vmax, new_hmax, old_vmax, old_hmax;
  RecT r;
  Int4 x, y;
  
  if (dlg == NULL)
  {
    return;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->panel);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->panel);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;
    
  if (dlg->row_list == NULL)
  {
    num_rows = 0;
    num_columns = 0;
  }
  else
  {
    num_rows = ValNodeLen (dlg->row_list);
    num_columns = ValNodeLen (dlg->row_list->data.ptrvalue);
  }

  ObjectRect (dlg->panel, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;
    
  visible_rows = (r.bottom - r.top - 2 * dlg->table_inset) / stdLineHeight - dlg->frozen_header;
  new_vmax = num_rows - visible_rows - 1;
  new_hmax = num_columns - dlg->frozen_left - 1;
  if (new_vmax < 0)
  {
    new_vmax = 0;
  }
  if (new_hmax < 0)
  {
    new_hmax = 0;
  }
  old_vmax = GetBarMax (sb_vert);
  old_hmax = GetBarMax (sb_horiz);
  
  if (old_vmax != new_vmax)
  {
    CorrectBarMax (sb_vert, new_vmax);
    if (start_row > new_vmax + dlg->frozen_header)
    {
      start_row = new_vmax + dlg->frozen_header;
    }
    CorrectBarValue (sb_vert, start_row - dlg->frozen_header);
    CorrectBarPage (sb_vert, 1, 1);
  }
  
  if (old_hmax != new_hmax)
  {
    CorrectBarMax (sb_horiz, new_hmax);
    if (start_col > new_hmax + dlg->frozen_left)
    {
      start_col = new_hmax + dlg->frozen_left;
    }
    CorrectBarValue (sb_horiz, start_col - dlg->frozen_left);
    CorrectBarPage (sb_horiz, 1, 1);
  }  
}

static void RowsToTableDisplayDialog (DialoG d, Pointer userdata)
{
  TableDisplayPtr dlg;
  RecT            r;
	WindoW          temport;

  dlg = (TableDisplayPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
  dlg->row_list = CopyTableDisplayRowList (userdata);
  UpdateTableDisplayDialogScrollBars (dlg);
  temport = SavePort (dlg->panel);
  Select (dlg->panel);
  ObjectRect (dlg->panel, &r);
  InvalRect (&r);  
  Update ();
  RestorePort (temport);
}

static Pointer TableDisplayDialogToRows (DialoG d)
{
  TableDisplayPtr dlg;

  dlg = (TableDisplayPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  return CopyTableDisplayRowList (dlg->row_list);
}

static Int4 
PrepareTableDisplayTextBuffer 
(CharPtr buf,
 Int4    remaining_chars_in_row,
 Int4    col_width,
 CharPtr data)
{
  Uint4 chars_to_paint = 0;
  if (buf == NULL)
  {
    return 0;
  }
  
  if (remaining_chars_in_row < col_width)
  {
    chars_to_paint = remaining_chars_in_row;
    StringNCpy (buf, data, chars_to_paint);
    buf [chars_to_paint] = 0;
  }
  else
  {
    chars_to_paint = col_width;
    StringNCpy (buf, data, chars_to_paint);
    buf [chars_to_paint] = 0;
    if (StringLen (data) > chars_to_paint && chars_to_paint > 2)
    {
      buf [chars_to_paint - 1] = '.';
      buf [chars_to_paint - 2] = '.';
      buf [chars_to_paint - 3] = '.';
    }
  }
  return chars_to_paint;
}

static void DrawTableDisplayLine (Int4 x, Int4 y, 
 ValNodePtr header_row,
 ValNodePtr data_row,
 CharPtr    buf,
 Int4       row_length,
 Int4       frozen_left,
 Int4       start_col,
 Int4       char_width,
 Int4       descent,
 Boolean    left_in_red)
{
  ValNodePtr header_vnp, data_vnp;
  Int4       x_offset, chars_to_paint, col_num;
  PoinT      pt1, pt2;
  RecT       rct;
  
  /* draw left margin */
  
  for (header_vnp = header_row, data_vnp = data_row, x_offset = 0, col_num = 0;
       header_vnp != NULL && data_vnp != NULL && x_offset < row_length && col_num < frozen_left;
       header_vnp = header_vnp->next, data_vnp = data_vnp->next, col_num++)
  {
    Gray ();
    InvertColors ();
    if (left_in_red)
    {
      Red ();
    }
    else
    {
      White ();
    }
    chars_to_paint = PrepareTableDisplayTextBuffer (buf, 
                                   (row_length - x_offset) / char_width,
                                   header_vnp->choice,
                                   data_vnp->data.ptrvalue);
        
    LoadRect (&rct, x + x_offset, y + descent,
              x + x_offset + (chars_to_paint + 2) * char_width, 
              y - stdLineHeight + descent);
    EraseRect (&rct);

    PaintStringEx ( (CharPtr)buf, x + x_offset, y);
    x_offset += (chars_to_paint + 2) * char_width;
    InvertColors ();
    Black ();
  }
  
  if (frozen_left > 0)
  {
    pt1.x = x + x_offset - 1;
    pt1.y = y;
    pt2.x = x + x_offset - 1;
    pt2.y = y - stdLineHeight;
    DrawLine (pt1, pt2);
  }

  
  while (col_num < start_col && header_vnp != NULL && data_vnp != NULL)
  {
    col_num++;
    header_vnp = header_vnp->next;
    data_vnp = data_vnp->next;
  }
  
  /* draw unfrozen columns */
  while (header_vnp != NULL && data_vnp != NULL && x_offset < row_length)
  {
    chars_to_paint = MIN (header_vnp->choice, (row_length - x_offset)/char_width);
    StringNCpy (buf, data_vnp->data.ptrvalue, chars_to_paint);
    buf [chars_to_paint] = 0;
    chars_to_paint = PrepareTableDisplayTextBuffer (buf, 
                                   (row_length - x_offset) / char_width,
                                   header_vnp->choice,
                                   data_vnp->data.ptrvalue);

    PaintStringEx ( (CharPtr)buf, x + x_offset, y);
    x_offset += (chars_to_paint + 2) * char_width;
    header_vnp = header_vnp->next;
    data_vnp = data_vnp->next;
  }
}

static void OnDrawTableDisplay (PaneL p)
{
  TableDisplayPtr dlg;
  BaR             sb_vert, sb_horiz;
  Int4            start_row, start_col;
  RecT            r;
  Int4            x, y, row, row_length;
  CharPtr         row_buffer;
  ValNodePtr      row_vnp;
  PoinT           pt1, pt2;
  Boolean         left_in_red;

  dlg = (TableDisplayPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) p);
  sb_horiz = GetSlateHScrollBar ((SlatE) p);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;

  ObjectRect (p, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;

  SelectFont (programFont); 
  
  row_length = r.right - r.left - 2;
  row_buffer = (CharPtr) MemNew (((row_length / dlg->char_width) + 1) * sizeof (Char));
  
  for (row = 0, row_vnp = dlg->row_list;
       row < dlg->frozen_header && y <= r.bottom - 2 * dlg->table_inset && row_vnp != NULL;
       row++, row_vnp = row_vnp->next)
  {
    DrawTableDisplayLine (x, y, dlg->row_list->data.ptrvalue, row_vnp->data.ptrvalue,
                          row_buffer, row_length, dlg->frozen_left, start_col, 
                          dlg->char_width, dlg->descent, FALSE);
    y += stdLineHeight;
  }
  
  while (row < start_row && row_vnp != NULL)
  {
    row++;
    row_vnp = row_vnp->next;
  }
  
  while (row_vnp != NULL && y <= r.bottom - 2 * dlg->table_inset)
  {
    left_in_red = FALSE;
    if (dlg->left_in_red != NULL)
    {
      left_in_red = (dlg->left_in_red) (row, dlg->row_list, dlg->left_in_red_data);
    }
    DrawTableDisplayLine (x, y, dlg->row_list->data.ptrvalue, row_vnp->data.ptrvalue,
                          row_buffer, row_length, dlg->frozen_left, start_col,
                          dlg->char_width, dlg->descent, left_in_red);
    row_vnp = row_vnp->next;
    y += stdLineHeight;
    row++;
  }
  
  /* draw line to separate header from remaining lines */
  if (dlg->frozen_header > 0)
  {
    Black ();
    pt1.x = x;
    pt1.y = r.top + stdLineHeight + dlg->descent;
    pt2.x = x + row_length;
    pt2.y = r.top + stdLineHeight + dlg->descent;
    DrawLine (pt1, pt2);
  }
  

}

static void OnVScrollTableDisplay (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT   r;
  WindoW temport;

  temport = SavePort (s);
  Select (s);
  ObjectRect (s, &r);
  InvalRect (&r);  
  RestorePort (temport);
  Update ();
}

static void OnHScrollTableDisplay (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT   r;
  WindoW temport;

  temport = SavePort (s);
  Select (s);
  ObjectRect (s, &r);
  InvalRect (&r);  
  RestorePort (temport);
  Update ();
}

static PoinT GetTableDisplayCell (TableDisplayPtr dlg, PoinT pt)
{
  BaR sb_horiz;
  BaR sb_vert;
  Int4 start_row, start_col;
  RecT r;
  PoinT cell_coord;
  Int4  x, y;
  ValNodePtr header_vnp;
  Int4  col_width;
  
  cell_coord.x = 0;
  cell_coord.y = 0;
  
  if (dlg == NULL || dlg->row_list == NULL)
  {
    return cell_coord;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->panel);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->panel);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;
  
  ObjectRect (dlg->panel, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = pt.x - r.left;
  y = pt.y - r.top;
  
  cell_coord.y = y / stdLineHeight;
  
  if (cell_coord.y >= dlg->frozen_header)
  {
    cell_coord.y += GetBarValue (sb_vert);
  }

  header_vnp = dlg->row_list->data.ptrvalue;
  
  col_width = 0;
  while (header_vnp != NULL && col_width + (header_vnp->choice + 2) * dlg->char_width < x
         && cell_coord.x < dlg->frozen_left)
  {
    cell_coord.x++;
    col_width += (header_vnp->choice + 2) * dlg->char_width;
    header_vnp = header_vnp->next;
  }
  
  if (cell_coord.x >= dlg->frozen_left)
  {
    /* skip over unfrozen columns not currently displayed */
    while (header_vnp != NULL && cell_coord.x < start_col)
    {
      header_vnp = header_vnp->next;
      cell_coord.x++;
    }
  
    while (header_vnp != NULL && col_width + (header_vnp->choice + 2) * dlg->char_width < x)
    {
      cell_coord.x++;
      col_width += (header_vnp->choice + 2) * dlg->char_width;
      header_vnp = header_vnp->next;
    }
  }
  return cell_coord;
}

extern CharPtr GetRowListCellText (ValNodePtr row_list, Int4 row, Int4 column)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       row_num, col_num;
  
  if (row_list == NULL || row < 0 || column < 0)
  {
    return NULL;
  }
  
  for (row_vnp = row_list, row_num = 0;
       row_vnp != NULL && row_num < row;
       row_vnp = row_vnp->next, row_num++)
  {
  }
  if (row_num != row || row_vnp == NULL)
  {
    return NULL;
  }
  for (col_vnp = row_vnp->data.ptrvalue, col_num = 0;
       col_vnp != NULL && col_num < column;
       col_vnp = col_vnp->next, col_num++)
  {
  }
  if (col_num != column || col_vnp == NULL)
  {
    return NULL;
  }
  else
  {
    return StringSave (col_vnp->data.ptrvalue);
  }  
}

static CharPtr TableDisplayGetTextForCell (TableDisplayPtr dlg, PoinT pt)
{
  if (dlg == NULL || dlg->row_list == NULL || pt.x < 0 || pt.y < 0)
  {
    return NULL;
  }
  
  return GetRowListCellText (dlg->row_list, pt.y, pt.x);
}

static void TableDisplayOnClick (PaneL p, PoinT pt)
{
  TableDisplayPtr dlg;
  Boolean         dbl_click;
  PoinT           cell_coord;
  PoinT           header_coord;
  CharPtr         cell_text;
  CharPtr         header_text;
  
  dlg = (TableDisplayPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  dbl_click = dblClick;
  if (dbl_click && dlg->dbl_click != NULL)
  {
    cell_coord = GetTableDisplayCell (dlg, pt);
    cell_text = TableDisplayGetTextForCell (dlg, cell_coord);
    header_coord.x = cell_coord.x;
    header_coord.y = 0;
    header_text = TableDisplayGetTextForCell (dlg, header_coord);
    (dlg->dbl_click) (cell_coord, header_text, cell_text, dlg->dbl_click_data);
    MemFree (cell_text);
    MemFree (header_text);
  }
}

extern FonT GetTableDisplayDefaultFont (void)
{
  FonT display_font = NULL;
  
#ifdef WIN_MAC
  display_font = ParseFont ("Monaco, 9");
#endif
#ifdef WIN_MSWIN
  display_font = ParseFont ("Courier, 9");
#endif
#ifdef WIN_MOTIF
  display_font = ParseFont ("fixed, 12");
#endif  
  return display_font;
}

extern DialoG TableDisplayDialog (GrouP parent, Int4 width, Int4 height,
                                  Int4 frozen_header, Int4 frozen_left,
                                  TableDisplayDblClick dbl_click,
                                  Pointer dbl_click_data,
                                  TableDisplayLeftInRed left_in_red,
                                  Pointer left_in_red_data)
{
  TableDisplayPtr dlg;
  GrouP           p;
  
  dlg = (TableDisplayPtr) MemNew (sizeof (TableDisplayData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupTableDisplayDialog);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = RowsToTableDisplayDialog;
  dlg->fromdialog = TableDisplayDialogToRows;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  dlg->row_list = NULL;
  dlg->frozen_header = frozen_header;
  dlg->frozen_left = frozen_left;
  dlg->table_inset = 4;
  dlg->dbl_click = dbl_click;
  dlg->dbl_click_data = dbl_click_data;
  dlg->left_in_red = left_in_red;
  dlg->left_in_red_data = left_in_red_data;
  
  dlg->display_font = GetTableDisplayDefaultFont ();

  SelectFont (dlg->display_font);
  dlg->char_width  = CharWidth ('0');
  dlg->descent = Descent ();
  
  dlg->panel = AutonomousPanel4 (p, width, height, OnDrawTableDisplay,
                               OnVScrollTableDisplay, OnHScrollTableDisplay,
                               sizeof (TableDisplayData), NULL, NULL); 
  SetObjectExtra (dlg->panel, dlg, NULL);
  SetPanelClick(dlg->panel, TableDisplayOnClick, NULL, NULL, NULL);
  
  return (DialoG) p;  
}

typedef struct multiselectdialog
{
  DIALOG_MESSAGE_BLOCK
  DoC                  doc;
  ValNodePtr           selected_list;
  ParData              listPar;
  ColData              listCol;
  Int4                 num_choices;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;    
} MultiSelectDialogData, PNTR MultiSelectDialogPtr;

static void DataToMultiSelectionDialog (DialoG d, Pointer userdata)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           vnp;
  Boolean              all_selected = FALSE;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  dlg->selected_list = ValNodeFree (dlg->selected_list);
  for (vnp = (ValNodePtr) userdata; vnp != NULL && ! all_selected; vnp = vnp->next)
  {
    if (vnp->data.intvalue == 1)
    {
      all_selected = TRUE;
    }
    else
    {
      ValNodeAddInt (&(dlg->selected_list), vnp->data.intvalue, vnp->data.intvalue);
    }
  }
  if (all_selected)
  {
    dlg->selected_list = ValNodeFree (dlg->selected_list);
    ValNodeAddInt (&(dlg->selected_list), 1, 1);
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
}

static Pointer MultiSelectionDialogToData (DialoG d)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           output_list = NULL, vnp;
 
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  
  for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
  {
    ValNodeAddInt (&output_list, vnp->choice, vnp->choice);
  }
  return output_list;
}

static void 
AddChoiceToSelection 
(Int2                 choice_num,
 Boolean              remove_other_choices,
 MultiSelectDialogPtr dlg)
{
  ValNodePtr already_sel, vnp;
  
  if (dlg == NULL || dlg->doc == NULL || choice_num < 1)
  {
    return;
  }
  
  if (choice_num == 1)
  {
    remove_other_choices = TRUE;
  }
  else 
  {
    /* if we have added a choice other than "All", remove the "All"
     * selection if it was present.
     */
    ValNodeFree (ValNodeExtractList (&(dlg->selected_list), 1));
    InvalDocRows (dlg->doc, 1, 1, 1);
  }
  
  already_sel = ValNodeExtractList (&(dlg->selected_list), choice_num);
  if (already_sel == NULL)
  {
    if (remove_other_choices)
    {
      /* delete old selections */
      for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
      {
        InvalDocRows (dlg->doc, vnp->choice, 1, 1);
      }
      dlg->selected_list = ValNodeFree (dlg->selected_list);      
      ValNodeAddInt (&dlg->selected_list, choice_num, choice_num);
      InvalDocRows (dlg->doc, choice_num, 1, 1);
    }
    else
    {
      /* add new selection */
      ValNodeAddInt (&(dlg->selected_list), choice_num, choice_num);
      InvalDocRows (dlg->doc, choice_num, 1, 1);
    }
  }
  else
  {
    already_sel = ValNodeFree (already_sel);
    InvalDocRows (dlg->doc, choice_num, 1, 1); 
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void SelectChoice (DoC d, PoinT pt)
{
  Int2      item, row;
  Boolean   remove_other_choices;
  MultiSelectDialogPtr dlg;
  
  remove_other_choices = ! ctrlKey;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  AddChoiceToSelection (item, remove_other_choices, dlg);
}

static Boolean ChoiceHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           vnp;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
  
  for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == item)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static void ChoiceOnKey (SlatE s, Char ch)
{
  MultiSelectDialogPtr dlg;
  CharPtr              str;
  Int2                 start_pos;
  Boolean              found = FALSE;
  ValNodePtr           vnp;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (s);
  if (dlg == NULL) return;

  if ( (int) ch == 0 ) return;
  
  if (isalpha (ch))
  {
    /* find the position of the last choice we added */
    for (vnp = dlg->selected_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
    {}
    
    if (vnp == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      /* start pos is one less than document row */
      start_pos ++;
      /* want to start at row after top row */
      start_pos ++;
    }
    else
    {
      /* want to start at row after currently selected row */
      start_pos = vnp->choice + 1;
    }
    
    while (!found && start_pos <= dlg->num_choices)
    {
      str = GetDocText (dlg->doc, start_pos, 1, 1);
      if (tolower (str [0]) == tolower (ch))
      {
        SetOffset (dlg->doc, 0, start_pos - 1);
        AddChoiceToSelection (start_pos, TRUE, dlg);
        found = TRUE;
      }
      str = MemFree (str);
      start_pos ++;
    }
    if (!found)
    {
      /* start searching at the top of the list */
      start_pos = 1;
      while (!found && start_pos <= dlg->num_choices)
      {
        str = GetDocText (dlg->doc, start_pos, 1, 1);
        if (tolower (str [0]) == tolower (ch))
        {
          SetOffset (dlg->doc, 0, start_pos - 1);
          AddChoiceToSelection (start_pos, TRUE, dlg);
          found = TRUE;
        }
        str = MemFree (str);
        start_pos ++;
      }
    }
  }
  else if (ch == NLM_DOWN)
  {
    /* down key */
    if (dlg->selected_list == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      start_pos ++;
    }
    else
    {
      start_pos = dlg->selected_list->choice;
    }
    start_pos ++;
    
    if (start_pos <= dlg->num_choices)
    {
      SetOffset (dlg->doc, 0, start_pos - 1);
      AddChoiceToSelection (start_pos, TRUE, dlg); 
      InvalDocRows (dlg->doc, 0, 0, 0);
    }
  }
  else if (ch == NLM_UP)
  {
    /* up key */
    if (dlg->selected_list == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      start_pos ++;
    }
    else
    {
      start_pos = dlg->selected_list->choice;
    }
    start_pos --;
    if (start_pos > 0)
    {
      SetOffset (dlg->doc, 0, start_pos - 1);
      AddChoiceToSelection (start_pos, TRUE, dlg); 
      InvalDocRows (dlg->doc, 0, 0, 0);
    }
  }
}

static DialoG 
MultiSelectDialog 
(GrouP      parent,
 ValNodePtr choice_list,
 Int4       list_height,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  MultiSelectDialogPtr dlg;
  GrouP                p;
  Int4                 height;
  Int4                 width;
  RecT                 r;
  ValNodePtr           vnp;
  
  if (choice_list == NULL)
  {
    return NULL;
  }
  dlg = (MultiSelectDialogPtr) MemNew (sizeof (MultiSelectDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToMultiSelectionDialog;
  dlg->fromdialog = MultiSelectionDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;  
  
  dlg->selected_list = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  
  SelectFont (systemFont);
  width = 0;
  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
    width = MAX (width, StringWidth (vnp->data.ptrvalue) + 2);
  }
  
  dlg->num_choices = ValNodeLen (choice_list) + 1;
  
  height = LineHeight ();
  dlg->doc = DocumentPanel (p, width, height * list_height);
  SetObjectExtra (dlg->doc, dlg, NULL);
  
  dlg->listPar.openSpace = FALSE;
  dlg->listPar.keepWithNext = FALSE;
  dlg->listPar.keepTogether = FALSE;
  dlg->listPar.newPage = FALSE;
  dlg->listPar.tabStops = FALSE;
  dlg->listPar.minLines = 0;
  dlg->listPar.minHeight = 0;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  dlg->listCol.pixWidth = r.right - r.left;
  dlg->listCol.pixInset = 0;
  dlg->listCol.charWidth = 160;
  dlg->listCol.charInset = 0;
  dlg->listCol.font = systemFont;
  dlg->listCol.just = 'l';
  dlg->listCol.wrap = FALSE;
  dlg->listCol.bar = FALSE;
  dlg->listCol.underline = FALSE;
  dlg->listCol.left = FALSE;
  dlg->listCol.last = TRUE;  
  
	AppendText (dlg->doc, "All", &(dlg->listPar), &(dlg->listCol), programFont);
  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
	  AppendText (dlg->doc, vnp->data.ptrvalue, &(dlg->listPar), &(dlg->listCol), programFont);
  }
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocProcs (dlg->doc, SelectChoice, NULL, NULL, NULL);
  SetDocShade (dlg->doc, NULL, NULL, ChoiceHighlight, NULL);
  SetSlateChar ((SlatE) dlg->doc, ChoiceOnKey);
  InvalDocument (dlg->doc);

  return (DialoG) p;
}

typedef struct selectiondialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG     multi_select_dlg;
  LisT       list_ctrl;
  PopuP      popup_ctrl;
  Int4       num_choices;
  CharPtr    err_msg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SelectionDialogData, PNTR SelectionDialogPtr;

static void ResetSelectionDialog (SelectionDialogPtr dlg)
{  
  if (dlg != NULL)
  {
    if (dlg->multi_select_dlg != NULL)
    {
      PointerToDialog (dlg->multi_select_dlg, NULL);
    }
    else if (dlg->list_ctrl != NULL)
    {
      SetValue (dlg->list_ctrl, 0);
    }
    else if (dlg->popup_ctrl != NULL)
    {
      SetValue (dlg->popup_ctrl, 0);
    }
  } 
}

static void SelectionDialogChanged (LisT l)
{
  SelectionDialogPtr dlg;

  dlg = (SelectionDialogPtr) GetObjectExtra (l);
  if (dlg == NULL)
  {
    return;
  }
    
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  } 
}

static void SelectionDialogPopupChanged (PopuP p)
{
  SelectionDialogPtr dlg;
 
  dlg = (SelectionDialogPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
    
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  } 
}

static void SelectionListToSelectionDialog (DialoG d, Pointer userdata)
{
  SelectionDialogPtr dlg;
  ValNodePtr         selected_list;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSelectionDialog (dlg);  
  selected_list = (ValNodePtr) userdata;
  if (dlg->multi_select_dlg != NULL)
  {
    PointerToDialog (dlg->multi_select_dlg, selected_list);
  }
  else if (dlg->list_ctrl != NULL)
  {
    if (selected_list == NULL)
    {
      SetValue (dlg->list_ctrl, 0);
    }
    else
    {
      SetValue (dlg->list_ctrl, selected_list->data.intvalue);
    }
  }
  else if (dlg->popup_ctrl)
  {
    if (selected_list == NULL)
    {
      SetValue (dlg->popup_ctrl, 0);
    }
    else
    {
      SetValue (dlg->popup_ctrl, selected_list->data.intvalue);
    }
  }
}


static Pointer SelectionDialogToSelectionList (DialoG d)
{
  SelectionDialogPtr dlg;
  ValNodePtr         sel_list = NULL, vnp;
  Int4               i = 0;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  if (dlg->multi_select_dlg != NULL)
  {
    sel_list = (ValNodePtr) DialogToPointer (dlg->multi_select_dlg);
    if (sel_list != NULL && sel_list->choice == 1)
    {
      sel_list = ValNodeFree (sel_list);
      for (i = 2; i <= dlg->num_choices; i++)
      {
        ValNodeAddInt (&sel_list, 0, i - 1);
      }
    }
    else
    {
      for (vnp = sel_list; vnp != NULL; vnp = vnp->next)
      {
        vnp->choice = 0;
        vnp->data.intvalue = vnp->data.intvalue - 1;
      }
    }
  }
  else
  {
    if (dlg->list_ctrl != NULL)
    {
      i = GetValue (dlg->list_ctrl);
    }
    else if (dlg->popup_ctrl != NULL)
    {
      i = GetValue (dlg->popup_ctrl);
    }
    if (i > 0)
    {
      ValNodeAddInt (&sel_list, 0, i);
    }
  }
  return (Pointer) sel_list;
}

static void SelectionDialogMessage (DialoG d, Int2 mssg)

{
  SelectionDialogPtr dlg;
  ValNode            vn;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSelectionDialog (dlg);
        break;
      case VIB_MSG_ENTER :
        if (dlg->multi_select_dlg != NULL)
        {
          Select (dlg->multi_select_dlg);
        }
        else if (dlg->list_ctrl != NULL)
        {
          Select (dlg->list_ctrl);
        }
        else if (dlg->popup_ctrl != NULL)
        {
          Select (dlg->popup_ctrl);
        }
        break;
      case NUM_VIB_MSG + 1:
        if (dlg->multi_select_dlg != NULL)
        {
          vn.next = NULL;
          vn.choice = 1;
          vn.data.intvalue = 1;
          PointerToDialog (dlg->multi_select_dlg, &vn);
        }
        else if (dlg->list_ctrl != NULL)
        {
          SetItemStatus (dlg->list_ctrl, 1, TRUE);
        }
        else if (dlg->popup_ctrl != NULL)
        {
          SetValue (dlg->popup_ctrl, 1);
        }
        SelectionDialogChanged (dlg->list_ctrl);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSelectionDialog (DialoG d)

{
  SelectionDialogPtr dlg;
  ValNodePtr         head = NULL, vnp;
  Boolean            any_selected = FALSE;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (dlg->multi_select_dlg != NULL)
    {
      vnp = DialogToPointer (dlg->multi_select_dlg);
      if (vnp != NULL)
      {
        any_selected = TRUE;
        vnp = ValNodeFree (vnp);
      }
    }
    else if (dlg->list_ctrl != NULL)
    {
      if (GetValue (dlg->list_ctrl) > 0)
      {
        any_selected = TRUE;
      }
    }
    else if (dlg->popup_ctrl != NULL)
    {
      if (GetValue (dlg->popup_ctrl) > 0)
      {
        any_selected = TRUE;
      }
    }
    if (!any_selected)
    {
      head = AddStringToValNodeChain (head, dlg->err_msg, 1);
    }
  }
  return head;
}

/* err_msg is the message to put in the results from TestDialog if nothing is selected */
/* choice_list should be a valnode list of strings to use for the names of the choices. */
/* All is automatically included as a choice if allow_multi is true. */
/* The ValNodeList returned is a list of integers indicating the position of the item
 * in the list - 1 is the first item, 2 is the second item, etc. */
extern DialoG SelectionDialogExEx 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height,
 Boolean                  force_list,
 Boolean                  force_popup)

{
  SelectionDialogPtr  dlg;
  GrouP               p;
  ValNodePtr          vnp;
  Int4                num_choices;
  Int4                list_width = 8, item_width;

  if (choice_list == NULL)
  {
    return NULL;
  }
  
  dlg = (SelectionDialogPtr) MemNew (sizeof (SelectionDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 0, 2, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SelectionListToSelectionDialog;
  dlg->fromdialog = SelectionDialogToSelectionList;
  dlg->dialogmessage = SelectionDialogMessage;
  dlg->testdialog = TestSelectionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->err_msg = err_msg;
  
  num_choices = ValNodeLen (choice_list);
  
  if (allow_multi)
  {
    dlg->multi_select_dlg = MultiSelectDialog (p, choice_list, list_height,
                                               change_notify, change_userdata);
    dlg->num_choices = num_choices + 1;                                               
  }
  else
  {
    if (force_popup || (num_choices < 20 && ! force_list) || list_height == 1)
    {
      dlg->popup_ctrl = PopupList (p, TRUE, SelectionDialogPopupChanged);
      SetObjectExtra (dlg->popup_ctrl, dlg, NULL);
      for (vnp = choice_list; vnp != NULL; vnp = vnp->next) {
        PopupItem (dlg->popup_ctrl, vnp->data.ptrvalue);
      }
    }
    else
    {
      SelectFont (systemFont);
      for (vnp = choice_list; vnp != NULL; vnp = vnp->next) {
        item_width = StringWidth (vnp->data.ptrvalue);
        list_width = MAX (list_width, item_width);
      }
      /* add padding */
      list_width += StringWidth ("W");
      list_width = list_width / Nlm_stdCharWidth;
      dlg->list_ctrl = SingleList (p, list_width, list_height, SelectionDialogChanged);
      SetObjectExtra (dlg->list_ctrl, dlg, NULL);
      for (vnp = choice_list; vnp != NULL; vnp = vnp->next) {
        ListItem (dlg->list_ctrl, vnp->data.ptrvalue);
      }      
    }
    dlg->num_choices = num_choices;
  }
  
  return (DialoG) p;
}


extern DialoG SelectionDialogEx 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height,
 Boolean                  force_list)

{
  return SelectionDialogExEx (h, change_notify, change_userdata, allow_multi, err_msg, choice_list, list_height, force_list, FALSE);
}


extern DialoG SelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height)
{
  return SelectionDialogEx (h, change_notify, change_userdata, allow_multi,
                          err_msg, choice_list, list_height, FALSE);
}

typedef struct valnodeselection
{
  DIALOG_MESSAGE_BLOCK
  DialoG           list_dlg;
  ValNodePtr       choice_list;

  Boolean             is_multi;  
  FreeValNodeProc     free_vn_proc;
  CopyValNodeDataProc copy_vn_proc;
  MatchValNodeProc    match_vn_proc;
  RemapValNodeProc    remap_vn_proc;
  
} ValNodeSelectionData, PNTR ValNodeSelectionPtr;

static void ValNodeSelectionListToDialog (DialoG d, Pointer userdata)
{
  ValNodeSelectionPtr dlg;
  ValNodePtr          item_list, vnp_list, vnp_sel, pos_list = NULL;
  Int4                i;
  Boolean             found;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  /* reset list control */
  PointerToDialog (dlg->list_dlg, NULL);  
  
  item_list = (ValNodePtr) userdata;
  for (vnp_list = item_list; vnp_list != NULL; vnp_list = vnp_list->next)
  {
    found = FALSE;
    i = 1;
    if (dlg->is_multi) {
      i++;
    }
    for (vnp_sel = dlg->choice_list;
         vnp_sel != NULL && !found;
         vnp_sel = vnp_sel->next, i++)
    {
      if ((dlg->match_vn_proc)(vnp_sel, vnp_list))
      {
        found = TRUE;
        ValNodeAddInt (&pos_list, 0, i);
      }
    }
  }
  PointerToDialog (dlg->list_dlg, pos_list);
  ValNodeFree (pos_list);  
}

static Pointer ValNodeSelectionDialogToList (DialoG d)
{
  ValNodeSelectionPtr dlg;
  ValNodePtr          item_list = NULL, vnp_list, pos_list, vnp_pos;
  ValNodePtr          vnp_copy, vnp_last = NULL, vnp_test;
  Int4                i;
  Boolean             found;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  pos_list = DialogToPointer (dlg->list_dlg);
  for (vnp_pos = pos_list; vnp_pos != NULL; vnp_pos = vnp_pos->next)
  {
    for (i = 1, vnp_list = dlg->choice_list;
         i < vnp_pos->data.intvalue && vnp_list != NULL;
         i++, vnp_list = vnp_list->next)
    {
    }
    if (i == vnp_pos->data.intvalue && vnp_list != NULL)
    {
      /* make sure we don't already have this value in the list */
      for (vnp_test = item_list, found = FALSE;
           vnp_test != NULL && !found;
           vnp_test = vnp_test->next)
      {
        found = (dlg->match_vn_proc) (vnp_list, vnp_test);
      }
      
      if (found)
      {
        continue;
      }
      vnp_copy = (dlg->copy_vn_proc) (vnp_list);
      if (vnp_last == NULL)
      {
        item_list = vnp_copy;
      }
      else
      {
        vnp_last->next = vnp_copy;
      }
      vnp_last = vnp_copy;
    }
  }
  if (dlg->remap_vn_proc != NULL)
  {
    item_list = (dlg->remap_vn_proc) (item_list);
  }
  return item_list;  
}

static void CleanupValNodeSelectionDialogForm (GraphiC g, VoidPtr data)

{
  ValNodeSelectionPtr dlg;
  ValNodePtr          vnp;

  dlg = (ValNodeSelectionPtr) data;
  if (dlg != NULL) {
    if (dlg->free_vn_proc != NULL)
    {
      for (vnp = dlg->choice_list; vnp != NULL; vnp = vnp->next)
      {
        (dlg->free_vn_proc) (vnp);
      }
    }
    dlg->choice_list = ValNodeFree (dlg->choice_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ValNodeSelectionDialogMessage (DialoG d, Int2 mssg)

{
  ValNodeSelectionPtr dlg;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        PointerToDialog (dlg->list_dlg, NULL);
        break;
      case VIB_MSG_SELECT:
        Select (dlg->list_dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->list_dlg);
        break;
      case NUM_VIB_MSG + 1:
        SendMessageToDialog (dlg->list_dlg, NUM_VIB_MSG + 1);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestValNodeSelectionDialog (DialoG d)

{
  ValNodeSelectionPtr  dlg;
  ValNodePtr           head = NULL;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    head = TestDialog (dlg->list_dlg);
  }
  return head;
}

extern DialoG ValNodeSelectionDialogExEx
(GrouP h,
 ValNodePtr               choice_list,
 Int2                     list_height,
 NameFromValNodeProc      name_proc,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 MatchValNodeProc         match_vn_proc,
 CharPtr                  err_name,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  force_list,
 Boolean                  force_popup,
 RemapValNodeProc         remap_vn_proc)
{
  ValNodeSelectionPtr  dlg;
  GrouP                p;
  ValNodePtr           choice_name_list = NULL, vnp;

  if (choice_list == NULL || name_proc == NULL
      || copy_vn_proc == NULL || match_vn_proc == NULL)
  {
    return NULL;
  }
  
  dlg = (ValNodeSelectionPtr) MemNew (sizeof (ValNodeSelectionData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupValNodeSelectionDialogForm);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = ValNodeSelectionListToDialog;
  dlg->fromdialog = ValNodeSelectionDialogToList;
  dlg->dialogmessage = ValNodeSelectionDialogMessage;
  dlg->testdialog = TestValNodeSelectionDialog;
  
  dlg->choice_list = choice_list;
  dlg->free_vn_proc = free_vn_proc;
  dlg->copy_vn_proc = copy_vn_proc;
  dlg->match_vn_proc = match_vn_proc;
  dlg->remap_vn_proc = remap_vn_proc;

  dlg->is_multi = allow_multi;

  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
    ValNodeAddPointer (&choice_name_list, 0, (name_proc) (vnp));
  }

  dlg->list_dlg = SelectionDialogExEx (p, change_notify, change_userdata,
                                   allow_multi, err_name, choice_name_list, 
                                   list_height, force_list, force_popup);
  ValNodeFreeData (choice_name_list);  
  
  return (DialoG) p;
}


extern DialoG ValNodeSelectionDialogEx
(GrouP h,
 ValNodePtr               choice_list,
 Int2                     list_height,
 NameFromValNodeProc      name_proc,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 MatchValNodeProc         match_vn_proc,
 CharPtr                  err_name,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  force_list,
 RemapValNodeProc         remap_vn_proc)
{
  return ValNodeSelectionDialogExEx (h, choice_list, list_height, name_proc, free_vn_proc, copy_vn_proc,
                                     match_vn_proc, err_name, change_notify, change_userdata, allow_multi,
                                     force_list, FALSE, remap_vn_proc);
}


extern DialoG ValNodeSelectionDialog
(GrouP h,
 ValNodePtr               choice_list,
 Int2                     list_height,
 NameFromValNodeProc      name_proc,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 MatchValNodeProc         match_vn_proc,
 CharPtr                  err_name,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi)
{
  return ValNodeSelectionDialogEx (h, choice_list, list_height,
                                   name_proc, free_vn_proc,
                                   copy_vn_proc, match_vn_proc, err_name,
                                   change_notify, change_userdata,
                                   allow_multi, FALSE, NULL);
}

extern DialoG EnumAssocSelectionDialog 
(GrouP                 h,
 Nlm_EnumFieldAssocPtr eap,
 CharPtr               err_name,
 Boolean               allow_multi,
 Nlm_ChangeNotifyProc  change_notify,
 Pointer               change_userdata)

{
  DialoG     dlg;
  ValNodePtr choice_list = NULL;
  
  if (eap == NULL)
  {
    return NULL;
  }

  while (eap->name != NULL)
  {
    if (!StringHasNoText (eap->name))
    {
      ValNodeAddPointer (&choice_list, eap->value, StringSave (eap->name));
    }
    eap++;
  }
  
  /* note - the ValNodeSelectionDialog will free the qual_choice_list when done */                                            
  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, err_name, 
                                change_notify, change_userdata, allow_multi);

  return dlg;
}

extern CharPtr ValNodeStringName (ValNodePtr vnp)
{
  if (vnp == NULL || vnp->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else
  {
    return StringSave (vnp->data.ptrvalue);
  }
}

extern void ValNodeSimpleDataFree (ValNodePtr vnp)
{
  if (vnp != NULL && vnp->data.ptrvalue != NULL)
  {
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
  }
}

extern ValNodePtr ValNodeStringCopy (ValNodePtr vnp)
{
  ValNodePtr vnp_copy = NULL;
  if (vnp != NULL)
  {
    ValNodeAddPointer (&vnp_copy, vnp->choice, StringSave (vnp->data.ptrvalue));
  }
  return vnp_copy;
}

extern Boolean ValNodeChoiceMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  if (vnp1->choice == vnp2->choice)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

extern Boolean ValNodeStringMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

typedef struct sequenceselection
{
  DIALOG_MESSAGE_BLOCK
  DialoG     sequence_list_dlg;
  ValNodePtr sequence_choice_list;
} SequenceSelectionData, PNTR SequenceSelectionPtr;

static void CleanupSequenceSelectionDialogForm (GraphiC g, VoidPtr data)

{
  SequenceSelectionPtr dlg;

  dlg = (SequenceSelectionPtr) data;
  if (dlg != NULL) {
    dlg->sequence_choice_list = ValNodeFree (dlg->sequence_choice_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ResetSequenceSelectionDialog (SequenceSelectionPtr dlg)
{  
  if (dlg != NULL)
  {
    PointerToDialog (dlg->sequence_list_dlg, NULL);
  }
}

static void SequenceSelectionListToSequenceSelectionDialog (DialoG d, Pointer userdata)
{
  SequenceSelectionPtr dlg;
  ValNodePtr           sequence_list, vnp_list, vnp_sel, pos_list = NULL;
  Int4                 i;
  SeqIdPtr             sip;
  Boolean              found;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSequenceSelectionDialog (dlg);  
  sequence_list = (ValNodePtr) userdata;
  for (vnp_list = sequence_list; vnp_list != NULL; vnp_list = vnp_list->next)
  {
    sip = (SeqIdPtr) vnp_list->data.ptrvalue;
    found = FALSE;
    while (sip != NULL && ! found)
    {
      for (vnp_sel = dlg->sequence_choice_list, i = 1;
           vnp_sel != NULL && !found;
           vnp_sel = vnp_sel->next, i++)
      {
        found = SeqIdIn (sip, vnp_sel->data.ptrvalue);
        if (found)
        {
          ValNodeAddInt (&pos_list, 0, i);
        }
      }
      sip = sip->next;
    }
  }
  PointerToDialog (dlg->sequence_list_dlg, pos_list);
  ValNodeFree (pos_list);
}

static Pointer SequenceSelectionDialogToSequenceSelectionList (DialoG d)
{
  SequenceSelectionPtr dlg;
  ValNodePtr           sequence_list = NULL, vnp_list, pos_list, vnp_pos;
  Int4                 i;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  pos_list = DialogToPointer (dlg->sequence_list_dlg);
  for (vnp_pos = pos_list; vnp_pos != NULL; vnp_pos = vnp_pos->next)
  {
    for (i = 1, vnp_list = dlg->sequence_choice_list;
         i < vnp_pos->data.intvalue && vnp_list != NULL;
         i++, vnp_list = vnp_list->next)
    {
    }
    if (i == vnp_pos->data.intvalue && vnp_list != NULL)
    {
      ValNodeAddPointer (&sequence_list, 0, vnp_list->data.ptrvalue);
    }
  }
  return sequence_list;
}

static void 
GetSequenceChoiceList 
(SeqEntryPtr sep,
 ValNodePtr PNTR list, 
 Boolean show_nucs, 
 Boolean show_prots)
{
  BioseqPtr                bsp;
  BioseqSetPtr             bssp;
  
  if (sep == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    if (!show_nucs && ISA_na (bsp->mol))
    {
      return;
    }
    if (!show_prots && ISA_aa (bsp->mol))
    {
      return;
    }
    ValNodeAddPointer (list, 0, bsp->id);
  }
  else
  {
  	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
    {
      GetSequenceChoiceList (sep, list, show_nucs, show_prots);
    }
  }
}

static void SequenceSelectionDialogMessage (DialoG d, Int2 mssg)

{
  SequenceSelectionPtr dlg;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSequenceSelectionDialog (dlg);
        break;
      case VIB_MSG_SELECT:
        Select (dlg->sequence_list_dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->sequence_list_dlg);
        break;
      case NUM_VIB_MSG + 1:
        SendMessageToDialog (dlg->sequence_list_dlg, NUM_VIB_MSG + 1);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSequenceSelectionDialog (DialoG d)

{
  SequenceSelectionPtr dlg;
  ValNodePtr           head = NULL;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    head = TestDialog (dlg->sequence_list_dlg);
  }
  return head;
}


extern DialoG SequenceSelectionDialogEx 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  allow_none,
 Boolean                  show_nucs,
 Boolean                  show_prots,
 Uint2                    entityID,
 Int4                     list_height)

{
  SequenceSelectionPtr  dlg;
  GrouP                 p;
  ValNodePtr                vnp;
  SeqEntryPtr               sep;
  SeqIdPtr                  sip;
  Char                      tmp[128];
  ValNodePtr            choice_name_list = NULL;
  
  if (!show_nucs && ! show_prots)
  {
    return NULL;
  }
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL)
  {
    return NULL;
  }

  dlg = (SequenceSelectionPtr) MemNew (sizeof (SequenceSelectionData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupSequenceSelectionDialogForm);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceSelectionListToSequenceSelectionDialog;
  dlg->fromdialog = SequenceSelectionDialogToSequenceSelectionList;
  dlg->dialogmessage = SequenceSelectionDialogMessage;
  dlg->testdialog = TestSequenceSelectionDialog;

  if (allow_none) {
    dlg->sequence_choice_list = ValNodeNew (NULL);
    dlg->sequence_choice_list->choice = 0;
    dlg->sequence_choice_list->data.ptrvalue = NULL;
  } else {
    dlg->sequence_choice_list = NULL;
  }
  GetSequenceChoiceList (sep, &dlg->sequence_choice_list, show_nucs, show_prots);
  
  
  for (vnp = dlg->sequence_choice_list; vnp != NULL; vnp = vnp->next) {
    sip = SeqIdFindWorst ((SeqIdPtr) vnp->data.ptrvalue);
    if (sip == NULL) {
      sprintf (tmp, " ");
    } else {
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
    }
    ValNodeAddPointer (&choice_name_list, 0, StringSave (tmp));
  }

  dlg->sequence_list_dlg = SelectionDialog (p, change_notify, change_userdata,
                                            allow_multi, "sequence",
                                            choice_name_list, list_height);
  ValNodeFreeData (choice_name_list); 
  return (DialoG) p;
}


extern DialoG SequenceSelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  show_nucs,
 Boolean                  show_prots,
 Uint2                    entityID)
{
  return SequenceSelectionDialogEx (h, change_notify, change_userdata, allow_multi, FALSE, show_nucs, show_prots, entityID, TALL_SELECTION_LIST);
}


extern DialoG SubSourceTypeDialog 
(GrouP                    h,
 Int2                     list_height,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  force_list,
 Boolean                  include_note)
{
  ValNodePtr subsource_list = NULL;
  Int4       i;

  for (i = 0; current_subsource_subtype_alist[i].name != NULL; i++) {
    ValNodeAddPointer (&subsource_list, current_subsource_subtype_alist[i].value, StringSave (current_subsource_subtype_alist[i].name));
  }
  if (include_note) {
    ValNodeAddPointer (&subsource_list, SUBSRC_other, StringSave ("Note"));
  }
  
  return ValNodeSelectionDialogEx (h, subsource_list, list_height, 
                                   ValNodeStringName,
                                   ValNodeSimpleDataFree, 
                                   ValNodeStringCopy,
                                   ValNodeChoiceMatch, "subsource list", 
                                   change_notify, change_userdata, allow_multi, force_list, NULL);
}


extern DialoG OrgModTypeDialog 
(GrouP                    h,
 Int2                     list_height,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  force_list,
 Boolean                  include_note)
{
  ValNodePtr orgmod_list = NULL;
  Int4       i;

  for (i = 0; current_orgmod_subtype_alist[i].name != NULL; i++) {
    ValNodeAddPointer (&orgmod_list, current_orgmod_subtype_alist[i].value, StringSave (current_orgmod_subtype_alist[i].name));
  }
  if (include_note) {
    ValNodeAddPointer (&orgmod_list, ORGMOD_other, StringSave ("Note"));
  }
  
  return ValNodeSelectionDialogEx (h, orgmod_list, list_height, 
                                   ValNodeStringName,
                                   ValNodeSimpleDataFree, 
                                   ValNodeStringCopy,
                                   ValNodeChoiceMatch, "orgmod list", 
                                   change_notify, change_userdata, allow_multi, force_list, NULL);
}


/*
static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  NULL
};

ENUM_ALIST(inference_alist)
  { " ",                     0 },
  { "similar to sequence",   1 },
  { "similar to protein",    2 },
  { "similar to DNA",        3 },
  { "similar to RNA",        4 },
  { "similar to mRNA",       5 },
  { "similar to EST",        6 },
  { "similar to other RNA",  7 },
  { "profile",               8 },
  { "nucleotide motif",      9 },
  { "protein motif",        10 },
  { "ab initio prediction", 11 },
END_ENUM_ALIST

Uint2 inference_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};

Uint2 inference_widths [] = {
  0, 0
};

static EnumFieldAssocPtr inference_popups [] = {
  inference_alist, NULL
};

extern void GBQualsToInferenceDialog (DialoG d, SeqFeatPtr sfp)

{
  Int2               best;
  Char               ch;
  GBQualPtr          gbq;
  ValNodePtr         head = NULL;
  Int2               j;
  ValNodePtr         last = NULL;
  size_t             len;
  CharPtr            rest;
  CharPtr            str;
  TagListPtr         tlp;
  Char               tmp [32];
  ValNodePtr         vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;

  if (sfp != NULL) {
    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringICmp (gbq->qual, "inference") != 0) continue;
      if (StringHasNoText (gbq->val)) continue;
      vnp = ValNodeNew (last);
      if (vnp == NULL) continue;
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;

      rest = NULL;
      best = -1;
      for (j = 0; inferencePrefix [j] != NULL; j++) {
        len = StringLen (inferencePrefix [j]);
        if (StringNICmp (gbq->val, inferencePrefix [j], len) != 0) continue;
        rest = gbq->val + len;
        best = j;
      }
      if (best >= 0 && inferencePrefix [best] != NULL) {
        if (rest != NULL) {
          ch = *rest;
          while (IS_WHITESP (ch) || ch == ':') {
            rest++;
            ch = *rest;
          }
        }
        len = StringLen (rest);
        str = MemNew (len + 16);
        if (str != NULL) {
          sprintf (tmp, "%d", (int) best);
          StringCpy (str, tmp);
          StringCat (str, "\t");
          StringCat (str, rest);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      } else {
        len + StringLen (gbq->val);
        str = MemNew (len + 8);
        if (str != NULL) {
          StringCpy (str, "0");
          StringCat (str, "\t");
          StringCat (str, gbq->val);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
    }
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) continue;
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

static void VisStringDialogToGbquals (SeqFeatPtr sfp, DialoG d, CharPtr qual)

{
  GBQualPtr   gbq, gbqlast = NULL;
  ValNodePtr  head = NULL, vnp;
  CharPtr     str;

  if (sfp == NULL || StringHasNoText (qual)) return;
  head = DialogToPointer (d);
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    gbqlast = gbq;
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    gbq = GBQualNew ();
    if (gbq == NULL) continue;
    gbq->qual = StringSave (qual);
    gbq->val = StringSave (str);
    if (gbqlast == NULL) {
      sfp->qual = gbq;
    } else {
      gbqlast->next = gbq;
    }
    gbqlast = gbq;
  }
  ValNodeFreeData (head);
}

extern void InferenceDialogToGBQuals (DialoG d, SeqFeatPtr sfp)

{
  GBQualPtr   gbq;
  GBQualPtr   gbqlast = NULL;
  Int2        j;
  size_t      len;
  CharPtr     prefix;
  CharPtr     ptr;
  CharPtr     rest;
  CharPtr     str;
  TagListPtr  tlp;
  Int2        val;
  ValNodePtr  vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL || sfp == NULL) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    gbqlast = gbq;
  }

  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    if (StringHasNoText ((CharPtr) vnp->data.ptrvalue)) continue;
    ptr = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    TrimSpacesAroundString (ptr);
    prefix = NULL;
    if (StrToInt (ptr, &val)) {
      for (j = 0; inferencePrefix [j] != NULL; j++) {
        if (j == val) {
          prefix = inferencePrefix [j];
        }
      }
    }
    MemFree (ptr);
    rest = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
    TrimSpacesAroundString (rest);
    if (StringDoesHaveText (prefix)) {
      len = StringLen (prefix) + StringLen (rest);
      str = (CharPtr) MemNew (len + 8);
      if (str != NULL) {
        if (StringDoesHaveText (prefix)) {
          StringCpy (str, prefix);
          if (StringDoesHaveText (rest)) {
            if (StringNICmp (rest, "(same species)", 14) != 0) {
              StringCat (str, ":");
            } else {
              StringCat (str, " ");
            }
          }
        }
        if (StringDoesHaveText (rest)) {
          StringCat (str, rest);
        }
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("inference");
          gbq->val = str;
          if (gbqlast == NULL) {
            sfp->qual = gbq;
          } else {
            gbqlast->next = gbq;
          }
          gbqlast = gbq;
        }
      }
    }
    MemFree (rest);
  }
}

static DialoG CreateInferenceDialog (GrouP h, Uint2 rows, Int2 spacing, Int2 width)

{
  inference_widths [1] = width;
  return CreateTagListDialog (h, rows, 2, spacing,
                              inference_types, inference_widths,
                              inference_popups, NULL, NULL);
}
*/

/* ************************ */

/* inference dialog controls, utility functions */

Uint2 accessionlist_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};

Uint2 accessionlist_widths [] = {
  0, 10
};

ENUM_ALIST(accn_type_alist)
  { " ",       0 },
  { "GenBank", 1 },
  { "EMBL",    2 },
  { "DDBJ",    3 },
  { "INSD",    4 },
  { "RefSeq",  5 },
  { "UniProt", 6 },
  { "Other",   7 },
END_ENUM_ALIST

static EnumFieldAssocPtr accessionlist_popups [] = {
  accn_type_alist, NULL
};

static CharPtr accnTypePrefix [] = {
  "",
  "GenBank",
  "EMBL",
  "DDBJ",
  "INSD",
  "RefSeq",
  "UniProt",
  "?",
  NULL
};

const Int4 numAccnTypePrefixes = sizeof (accnTypePrefix) / sizeof (CharPtr);

static Int4 GetAccnTypeNum (CharPtr str)
{
  Int4 i;

  if (StringHasNoText (str)) return 0;

  for (i = 1; i < numAccnTypePrefixes; i++)
  {
    if (StringCmp (accnTypePrefix[i], str) == 0)
    {
      return i;
    }
  }
  return 0;
}


static CharPtr ValForOneAccession (CharPtr str)
{
  CharPtr cp, val_buf = NULL;
  CharPtr val_fmt = "%d\t%s";
  Int4    db;

  if (!StringHasNoText (str))
  {
    cp = StringChr (str, '|');
    if (cp == NULL)
    {
      if ((db = GetAccnTypeNum(str)) > 0)
      {
        val_buf = MemNew (sizeof (Char) * StringLen (val_fmt));
        sprintf (val_buf, val_fmt, db, " ");
      }
      else
      {
        val_buf = MemNew (sizeof (Char) * (StringLen (val_fmt) + StringLen (str)));
        sprintf (val_buf, val_fmt, 0, str);
      }
    }
    else
    {
      *cp = 0;
      db = GetAccnTypeNum (str);
      val_buf = MemNew (sizeof (Char) * (StringLen (val_fmt) + StringLen (cp + 1)));
      sprintf (val_buf, val_fmt, db, cp + 1);
      *cp = '|';
    }
  }
  return val_buf;
}

static void AccessionListDataToDialog (DialoG d, Pointer data)
{
  TagListPtr    tlp;
  CharPtr       str, cp, val_buf;
  ValNodePtr    new_list = NULL, vnp;
  Int4          j;
  Int2          scroll_pos;

  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL) return;
  str = (CharPtr) data;
  
  cp = StringChr (str, ',');
  while (cp != NULL)
  {
    *cp = 0;
    val_buf = ValForOneAccession (str);
    if (val_buf != NULL)
    {
      ValNodeAddPointer (&new_list, 0, val_buf);
    }
    *cp = ',';
    str = cp + 1;
    cp = StringChr (str, ',');
  }
  val_buf = ValForOneAccession (str);
  if (val_buf != NULL)
  {
    ValNodeAddPointer (&new_list, 0, val_buf);
  }

  scroll_pos = 0;
  if (tlp->bar != NULL) 
  {
    scroll_pos = GetBarValue (tlp->bar);
  }
  else if (tlp->left_bar != NULL)
  {
    scroll_pos = GetBarValue (tlp->left_bar);
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = new_list;
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, (Int2) (tlp->rows-1), (Int2) (tlp->rows-1));   

  /* retain scroll position */
  if (scroll_pos > tlp->max) {
    scroll_pos = tlp->max;
  }
  
  if (tlp->bar != NULL)
  {
    CorrectBarValue (tlp->bar, scroll_pos);      
  }
  if (tlp->left_bar != NULL)
  {
    CorrectBarValue (tlp->left_bar, scroll_pos);      
  }
    
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  Update ();
 
}


static Pointer AccessionListDialogToData (DialoG d)
{
  TagListPtr    tlp;
  ValNodePtr    vnp;
  Int4          result_len = 0, db;
  CharPtr       str, acc_str, result_str = NULL;
  Boolean       first_item = TRUE;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL) return NULL;
  
  for (vnp = tlp->vnp;
       vnp != NULL;
       vnp = vnp->next)
  {
    acc_str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
    str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    if (!StringHasNoText (acc_str) || !StringHasNoText (str))
    {
      result_len += StringLen (acc_str);
      result_len += 1; /* comma */
      db = atoi (str);
      if (db >= 0 && db < numAccnTypePrefixes)
      {
        result_len += MAX (StringLen (accnTypePrefix[db]), 1);
      }
      else
      {
        result_len ++; /* will represent with ? */
      }
      result_len ++; /* for db/accession separator */
    }
    acc_str = MemFree (acc_str);
    str = MemFree (str);
  }
  
  if (result_len > 0) 
  {
    result_str = (CharPtr) MemNew (sizeof (Char) * (result_len + 1));
    for (vnp = tlp->vnp;
        vnp != NULL;
        vnp = vnp->next)
    {
      acc_str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
      str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
      db = (str == NULL ? 0 : atoi (str));
      str = MemFree (str);
      if (!StringHasNoText (acc_str) || db > 0)
      {
        if (first_item)
        {
          first_item = FALSE;
        }
        else
        {
          StringCat (result_str, ",");
        }

        if (db > 0 && db < numAccnTypePrefixes)
        {
          StringCat (result_str, accnTypePrefix[db]);
        }
        else
        {
          StringCat (result_str, "?");
        }
        StringCat (result_str, "|");
        if (!StringHasNoText (result_str)) {
          StringCat (result_str, acc_str);
        }
      }
      acc_str = MemFree (acc_str);
    }
    result_str[result_len] = 0;
  }

  return result_str;
}


static void RemoveEmptyAccessionStrings (CharPtr acc_list)
{
  CharPtr cp_prev_end = NULL, cp_src, cp_dst;

  if (acc_list == NULL) return;
  
  cp_src = acc_list;
  cp_dst = acc_list;
  cp_prev_end = acc_list;
  while (*cp_src != 0)
  {
    if (*cp_src == '|' && (*(cp_src + 1) == ',' || *(cp_src + 1) == 0))
    {
      cp_dst = cp_prev_end;
      cp_prev_end = cp_dst;
      cp_src++;
    }
    else
    {
      *cp_dst = *cp_src;      
      if (*cp_src == ',')
      {
        cp_prev_end = cp_dst;
      }
      cp_dst++;
      cp_src++;
    }
  }
  *cp_dst = 0;
}


static CharPtr insdmessage =
"GenBank, EMBL, and DDBJ records are part of the International Nucleotide " \
"Sequence Database collaboration.\nThe database prefix for the /inference " \
"qualifier in these cases is INSD by collaboration policy.";


static Boolean ReplaceDatabaseStrings (CharPtr PNTR str)
{
  Boolean changed_db = FALSE;

  if (str == NULL || *str == NULL) return FALSE;

  if (StringSearch (*str, "GenBank|") != NULL)
  {
    changed_db = TRUE;
    FindReplaceString (str, "GenBank|", "INSD|", TRUE, FALSE);
  }
  if (StringSearch (*str, "EMBL|") != NULL)
  {
    changed_db = TRUE;
    FindReplaceString (str, "EMBL|", "INSD|", TRUE, FALSE);
  }
  if (StringSearch (*str, "DDBJ|") != NULL)
  {
    changed_db = TRUE;
    FindReplaceString (str, "DDBJ|", "INSD|", TRUE, FALSE);
  }

  if (changed_db)
  {
    if (GetAppProperty ("InternalNcbiSequin") == NULL) {
      Message (MSG_OK, "%s", insdmessage);
    }
  }

  return changed_db;
}


static void ChangeInferAccessionList (Pointer data); 
static TaglistCallback AccessionListCallbacks[] = 
{ ChangeInferAccessionList, ChangeInferAccessionList };

typedef struct inferevid {
  CharPtr  prefix;     /* from inferencePrefix     */
  Boolean  species;    /* optional (same species)  */
  CharPtr  database;   /* INSD, RefSeq, etc.       */
  CharPtr  db_other;   /* other database           */
  CharPtr  accession;  /* accession.version        */
  CharPtr  program;    /* common analysis program  */
  CharPtr  pr_other;   /* other program            */
  CharPtr  version;    /* program version          */
  CharPtr  basis1;     /* profile or motif         */
  CharPtr  basis2;     /*  evidence_basis texts    */
  CharPtr  accession_list; /* accession list for alignment */
} InferEvid, PNTR InferEvidPtr;

typedef struct inferdialog {
  DIALOG_MESSAGE_BLOCK

  DoC           inferdoc;
  Int2          currItem;

  PopuP         prefix;
  ButtoN        species;
  PopuP         database;
  TexT          db_other;
  TexT          accession;
  PopuP         program;
  TexT          pr_other;
  TexT          version;
  TexT          basis1;
  TexT          basis2;
  PrompT        inf_free_program_prompt;
  PrompT        inf_free_version_prompt;
  PrompT        accession_list_program_prompt;
  PrompT        accession_list_version_prompt;
  PrompT        accession_list_prompt;
  DialoG        accession_list;

  GrouP         inf_accn_group;
  GrouP         other_db_group;
  GrouP         inf_prog_group;
  GrouP         other_pr_group;
  GrouP         inf_free_group;

  Int2          numInf;
  InferEvidPtr  evidence [128];

} InferDialog, PNTR InferDialogPtr;

static InferEvidPtr InferEvidNew (
  void
)

{
  InferEvidPtr  iep;

  iep = MemNew (sizeof (InferEvid));
  if (iep == NULL) return NULL;

  return iep;
}

static InferEvidPtr InferEvidFree (
  InferEvidPtr iep
)

{
  if (iep == NULL) return NULL;

  MemFree (iep->prefix);
  MemFree (iep->database);
  MemFree (iep->accession);
  MemFree (iep->program);
  MemFree (iep->version);
  MemFree (iep->basis1);
  MemFree (iep->basis2);
  MemFree (iep->accession_list);

  return MemFree (iep);
}

static InferEvidPtr GetInferEvid (
  InferDialogPtr idp,
  Int2 item
)

{
  InferEvidPtr  iep;

  if (idp == NULL || item < 0 || item > 127) return NULL;
  iep = idp->evidence [item];
  if (iep != NULL) return iep;

  iep = InferEvidNew ();
  if (iep != NULL) {
    /*
    iep->prefix = StringSave (" ");
    iep->database = StringSave (" ");
    iep->db_other = StringSave ("");
    iep->accession = StringSave ("");
    iep->program = StringSave (" ");
    iep->pr_other = StringSave ("");
    iep->version = StringSave ("");
    iep->basis1 = StringSave ("");
    iep->basis2 = StringSave ("");
    */
  }
  idp->evidence [item] = iep;
  return iep;
}

/* inference DoC object tables */

#define NUM_INFERENCE_LINES 3

static ParData  inferParFmt = { FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0 };

static ColData  inferColFmt [] = {
  {0, 5, 25, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* class     */
  {0, 5, 25, 2, NULL, 'l', FALSE, TRUE,  FALSE, FALSE, TRUE}   /* specifics */
};

static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  "alignment",
  NULL
};

ENUM_ALIST(inference_alist)
  { " ",                     0 },
  { "similar to sequence",   1 },
  { "similar to protein",    2 },
  { "similar to DNA",        3 },
  { "similar to RNA",        4 },
  { "similar to mRNA",       5 },
  { "similar to EST",        6 },
  { "similar to other RNA",  7 },
  { "profile",               8 },
  { "nucleotide motif",      9 },
  { "protein motif",        10 },
  { "ab initio prediction", 11 },
  { "alignment",            12 },
END_ENUM_ALIST

static CharPtr programPrefix [] = {
  "",
  "tRNAscan",
  "Genscan",
  "?",
  NULL
};

ENUM_ALIST(program_alist)
  { " ",        0 },
  { "tRNAscan", 1 },
  { "Genscan",  2 },
  { "Other",    3 },
END_ENUM_ALIST

static CharPtr PrintInferTable (
  DoC d,
  Int2 item,
  Pointer data
)

{
  CharPtr         buf;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  size_t          len;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL || item < 1 || item > 127) return NULL;
  iep = GetInferEvid (idp, item);
  if (iep == NULL) return NULL;

  len = StringLen (iep->prefix) + StringLen (iep->database) + StringLen (iep->db_other) + StringLen (iep->accession) +
        StringLen (iep->program) + StringLen (iep->pr_other) + StringLen (iep->version) + StringLen (iep->basis1) +
        StringLen (iep->basis2) + StringLen (iep->accession_list) + 50;
  buf = MemNew (len);
  if (buf == NULL) return NULL;

  if (StringHasNoText (iep->prefix)) {
    StringCat (buf, " \t \n");
    return buf;
  }

  StringCat (buf, iep->prefix);

  StringCat (buf, "\t");

  if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
    if (StringDoesHaveText (iep->accession)) {
      if (StringCmp (iep->database, "Other") == 0) {
        if (StringDoesHaveText (iep->db_other)) {
          StringCat (buf, iep->db_other);
          StringCat (buf, ":");
        }
      } else if (StringDoesHaveText (iep->database)) {
        StringCat (buf, iep->database);
        StringCat (buf, ":");
      }
      StringCat (buf, iep->accession);
    }
  } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
    if (StringCmp (iep->program, "Other") == 0) {
      if (StringDoesHaveText (iep->pr_other)) {
        StringCat (buf, iep->pr_other);
        if (StringDoesHaveText (iep->version)) {
          StringCat (buf, ":");
          StringCat (buf, iep->version);
        }
      }
    } else if (StringDoesHaveText (iep->program)) {
      StringCat (buf, iep->program);
      if (StringDoesHaveText (iep->version)) {
        StringCat (buf, ":");
        StringCat (buf, iep->version);
      }
    }
  } else if (StringDoesHaveText (iep->basis1)) {
    StringCat (buf, iep->basis1);
    if (StringDoesHaveText (iep->basis2)) {
      StringCat (buf, ":");
      StringCat (buf, iep->basis2);
      if (StringCmp (iep->prefix, "alignment") == 0 && 
          StringDoesHaveText (iep->accession_list)) {
        StringCat (buf, ":");
        StringCat (buf, iep->accession_list);
      }
    }
  } else {
    StringCat (buf, " ");
  }

  StringCat (buf, "\n");
  return buf;
}

static void ShowInferenceGroup (
  InferDialogPtr idp
)

{
  CharPtr  str;
  UIEnum   val;

  if (idp == NULL) return;
  if (GetEnumPopup (idp->prefix, inference_alist, &val)) {
    if (val >= 1 && val <= 7) {
      SafeHide (idp->inf_prog_group);
      SafeHide (idp->inf_free_group);
      SafeShow (idp->inf_accn_group);
      SafeShow (idp->species);
      str = GetEnumPopupByName (idp->database, accn_type_alist);
      if (StringCmp (str, "Other") == 0) {
        SafeShow (idp->other_db_group);
      } else {
        SafeHide (idp->other_db_group);
      }
      MemFree (str);
    } else if (val >= 8 && val <= 10) {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_prog_group);
      SafeShow (idp->inf_free_group);
      SafeHide (idp->accession_list);
      SafeShow (idp->inf_free_program_prompt);
      SafeShow (idp->inf_free_version_prompt);
      SafeHide (idp->accession_list_program_prompt);
      SafeHide (idp->accession_list_version_prompt);
      SafeHide (idp->accession_list_prompt);
    } else if (val == 11) {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_free_group);
      SafeShow (idp->inf_prog_group);
      str = GetEnumPopupByName (idp->program, program_alist);
      if (StringCmp (str, "Other") == 0) {
        SafeShow (idp->other_pr_group);
      } else {
        SafeHide (idp->other_pr_group);
      }
      MemFree (str);
    } else if (val == 12) {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_prog_group);
      SafeShow (idp->inf_free_group);
      SafeShow (idp->accession_list);
      SafeHide (idp->inf_free_program_prompt);
      SafeHide (idp->inf_free_version_prompt);
      SafeShow (idp->accession_list_program_prompt);
      SafeShow (idp->accession_list_version_prompt);
      SafeShow (idp->accession_list_prompt);
    } else {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_prog_group);
      SafeHide (idp->inf_free_group);
    }
  } else {
    SafeHide (idp->inf_accn_group);
    SafeHide (idp->species);
    SafeHide (idp->inf_prog_group);
    SafeHide (idp->inf_free_group);
  }
  Update ();
}

static void SafeSetEnumPopupByName (PopuP lst, EnumFieldAssocPtr al, CharPtr name)

{
  if (StringDoesHaveText (name)) {
    SetEnumPopupByName (lst, al, name);
  } else {
    SetEnumPopupByName (lst, al, " ");
  }
}

static void ChangeInferTableSelect (
  DoC d,
  Int2 item,
  Int2 row,
  Int2 col,
  Boolean dblClck
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            itemOld1, itemOld2;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL) return;
  if (item == 0 || row == 0 || col == 0) return;

  GetDocHighlight (d, &itemOld1, &itemOld2);
  SetDocHighlight (d, item, item);
  UpdateDocument (d, itemOld1, itemOld2);
  UpdateDocument (d, item, item);
  idp->currItem = item;

  iep = GetInferEvid (idp, item);
  if (iep != NULL) {
    ResetClip ();
    SafeSetEnumPopupByName (idp->prefix, inference_alist, iep->prefix);

    SafeSetStatus (idp->species, iep->species);
    SafeSetEnumPopupByName (idp->database, accn_type_alist, iep->database);
    SafeSetTitle (idp->db_other, iep->db_other);
    SafeSetTitle (idp->accession, iep->accession);

    SafeSetEnumPopupByName (idp->program, program_alist, iep->program);
    SafeSetTitle (idp->pr_other, iep->pr_other);
    SafeSetTitle (idp->version, iep->version);

    SafeSetTitle (idp->basis1, iep->basis1);
    SafeSetTitle (idp->basis2, iep->basis2);

    ReplaceDatabaseStrings (&(iep->accession_list));
    PointerToDialog (idp->accession_list, iep->accession_list);

    ShowInferenceGroup (idp);
  }

  Update ();
}

static void CheckExtendInferTable (
  InferDialogPtr idp
)

{
  Int2  numItems;

  if (idp == NULL) return;

  GetDocParams (idp->inferdoc, &numItems, NULL);
  if (idp->currItem == numItems) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
  }

  Update ();
}

static void ChangeInferPrefix (
  PopuP p
)

{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->prefix, inference_alist);
  iep->prefix = MemFree (iep->prefix);
  iep->prefix = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeSameSpecies (
  ButtoN b
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (b);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->species = (Boolean) (GetStatus (b));

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferDatabase (
  PopuP p
)


{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->database, accn_type_alist);
  if (StringCmp (str, "GenBank") == 0 ||
      StringCmp (str, "EMBL") == 0 ||
      StringCmp (str, "DDBJ") == 0) {
    if (GetAppProperty ("InternalNcbiSequin") == NULL) {
      Message (MSG_OK, "%s", insdmessage);
    }
    SetEnumPopupByName (idp->database, accn_type_alist, "INSD");
    str = MemFree (str);
    str = StringSave ("INSD");
  }
  iep->database = MemFree (iep->database);
  iep->database = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferDbOther (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->db_other = MemFree (iep->db_other);
  iep->db_other = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferAccession (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->accession = MemFree (iep->accession);
  iep->accession = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferProgram (
  PopuP p
)


{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->program, program_alist);
  iep->program = MemFree (iep->program);
  iep->program = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferPrOther (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->pr_other = MemFree (iep->pr_other);
  iep->pr_other = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferVersion (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->version = MemFree (iep->version);
  iep->version = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferBasis1 (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->basis1 = MemFree (iep->basis1);
  iep->basis1 = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferBasis2 (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->basis2 = MemFree (iep->basis2);
  iep->basis2 = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}


static void ChangeInferAccessionList (Pointer data)
{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) data;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->accession_list = MemFree (iep->accession_list);
  iep->accession_list = DialogToPointer (idp->accession_list);
  
  if (ReplaceDatabaseStrings (&(iep->accession_list)))
  {
    PointerToDialog (idp->accession_list, iep->accession_list);
  }

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}


static Boolean StringInList (CharPtr str, CharPtr PNTR list)

{
  Int2  i;

  if (str == NULL || list == NULL) return FALSE;

  for (i = 0; list [i] != NULL; i++) {
    if (StringICmp (str, list[i]) == 0) return TRUE;
  }

  return FALSE;
}

extern void GBQualsToInferenceDialog (DialoG d, SeqFeatPtr sfp)

{
  Int2            best;
  Char            ch;
  GBQualPtr       gbq;
  Int2            i, j, k;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  size_t          len;
  CharPtr         rest;
  CharPtr         str;
  CharPtr         tmp, tmp2;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL) return;

  if (sfp == NULL || sfp->qual == NULL) {
    Reset (idp->inferdoc);
    SetValue (idp->prefix, 0);
    SetStatus (idp->species, FALSE);
    SetValue (idp->database, 0);
    SetTitle (idp->db_other, "");
    SetTitle (idp->accession, "");
    SetValue (idp->program, 0);
    SetTitle (idp->pr_other, "");
    SetTitle (idp->version, "");
    SetTitle (idp->basis1, "");
    SetTitle (idp->basis2, "");
    SafeHide (idp->inf_accn_group);
    SafeHide (idp->inf_prog_group);
    SafeHide (idp->inf_free_group);
    idp->numInf = 0;
    idp->currItem = 1;
    for (i = 0; i < NUM_INFERENCE_LINES; i++) {
      AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                  &inferParFmt, inferColFmt, systemFont);
    }
    SetDocHighlight (idp->inferdoc, 1, 1);
    return;
  }

  idp->numInf = 0;
  idp->currItem = 1;
  Reset (idp->inferdoc);

  for (k = 0; k < 128; k++) {
    iep = idp->evidence [k];
    InferEvidFree (iep);
    idp->evidence [k] = NULL;
  }

  for (gbq = sfp->qual, k = 0; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    if (StringHasNoText (gbq->val)) continue;

    rest = NULL;
    best = -1;
    for (j = 0; inferencePrefix [j] != NULL; j++) {
      len = StringLen (inferencePrefix [j]);
      if (StringNICmp (gbq->val, inferencePrefix [j], len) != 0) continue;
      rest = gbq->val + len;
      best = j;
    }

    k++;
    iep = GetInferEvid (idp, k);
    if (iep == NULL) continue;

    str = NULL;
    if (best > 0 && inferencePrefix [best] != NULL) {
      iep->prefix = MemFree (iep->prefix);
      iep->prefix = StringSave(GetEnumName ((UIEnum) best, inference_alist));

      if (rest != NULL) {
        ch = *rest;
        while (IS_WHITESP (ch)) {
          rest++;
          ch = *rest;
        }
        if (StringNICmp (rest, "(same species)", 14) == 0) {
          iep->species = TRUE;
          rest += 14;
        } else {
          iep->species = FALSE;
        }
        ch = *rest;
        while (IS_WHITESP (ch) || ch == ':') {
          rest++;
          ch = *rest;
        }
      }
      if (StringDoesHaveText (rest)) {
        str = StringSave (rest);
      }
      tmp = StringChr (str, ':');
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        TrimSpacesAroundString (str);
        TrimSpacesAroundString (tmp);
      } else {
        TrimSpacesAroundString (str);
      }
      if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
        if (StringInList (str, accnTypePrefix)) {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull (str);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (tmp);
        } else if (tmp != NULL) {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull ("Other");
          iep->db_other = MemFree (iep->db_other);
          iep->db_other = StringSaveNoNull (str);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (tmp);
        } else {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull (" ");
          iep->db_other = MemFree (iep->db_other);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (str);
        }
      } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
        if (StringInList (str, programPrefix)) {
          iep->program = MemFree (iep->program);
          iep->program = StringSaveNoNull (str);
        } else {
          iep->program = MemFree (iep->program);
          iep->program = StringSaveNoNull ("Other");
          iep->pr_other = MemFree (iep->pr_other);
          iep->pr_other = StringSaveNoNull (str);
        }
        iep->version = MemFree (iep->version);
        iep->version = StringSaveNoNull (tmp);
      } else {
        iep->basis1 = MemFree (iep->basis1);
        iep->basis1 = StringSaveNoNull (str);
        tmp2 = NULL;
        if (StringCmp (iep->prefix, "alignment") == 0)
        {
          tmp2 = StringChr (tmp, ':');
          if (tmp2 != NULL) {
            *tmp2 = 0;
            tmp2++;
          }
        }
        iep->basis2 = MemFree (iep->basis2);
        iep->basis2 = StringSaveNoNull (tmp);
        iep->accession_list = MemFree (iep->accession_list);
        iep->accession_list = StringSaveNoNull (tmp2);
      }

    } else {
      iep->prefix = StringSave ("???");
      str = StringSave (gbq->val);
      tmp = StringChr (str, ':');
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        TrimSpacesAroundString (str);
        TrimSpacesAroundString (tmp);
      } else {
        TrimSpacesAroundString (str);
      }
      iep->basis1 = MemFree (iep->basis1);
      iep->basis1 = StringSaveNoNull (str);
      iep->basis2 = MemFree (iep->basis2);
      iep->basis2 = StringSaveNoNull (tmp);
    }

    MemFree (str);

    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);

    (idp->numInf)++;
  }

  AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
              &inferParFmt, inferColFmt, systemFont);
  k++;

  while (k < NUM_INFERENCE_LINES) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
    k++;
  }

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, 0, 0);

  ChangeInferTableSelect (idp->inferdoc, 1, 1, 1, FALSE);

  Update ();
}

extern void InferenceDialogToGBQuals (DialoG d, SeqFeatPtr sfp, Boolean convertBadToNote)

{
  CharPtr         first = NULL, second = NULL, accession_list = NULL;
  GBQualPtr       gbq, lastgbq;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            k, numItems;
  size_t          len;
  CharPtr         prefix = NULL;
  CharPtr         speciesies = NULL;
  CharPtr         str;
  UIEnum          val;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL || sfp == NULL) return;

  lastgbq = NULL;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    lastgbq = gbq;
  }

  GetDocParams (idp->inferdoc, &numItems, NULL);
  for (k = 1; k <= numItems; k++) {
    iep = GetInferEvid (idp, k);
    if (iep == NULL) continue;

    if (StringHasNoText (iep->prefix)) continue;
    gbq = GBQualNew ();
    if (gbq == NULL) continue;

    gbq->qual = StringSave ("inference");

    if (WhereInEnumPopup (inference_alist, iep->prefix, &val)) {
      if (val > 0 && val <= 12) {
        prefix = inferencePrefix [(int) val];
      }
    }
    speciesies = NULL;
    if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
      if (iep->species) {
        speciesies = " (same species)";
      }
      if (StringDoesHaveText (iep->accession)) {
        if (StringCmp (iep->database, "Other") == 0) {
          if (StringDoesHaveText (iep->db_other)) {
            first = iep->db_other;
          }
        } else if (StringDoesHaveText (iep->database)) {
          first = iep->database;
        }
        second = iep->accession;
      }
    } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
      if (StringCmp (iep->program, "Other") == 0) {
        if (StringDoesHaveText (iep->pr_other)) {
          first = iep->pr_other;
        }
      } else if (StringDoesHaveText (iep->program)) {
        first = iep->program;
      }
      second = iep->version;
    } else {
      if (StringDoesHaveText (iep->basis1)) {
        first = iep->basis1;
        second = iep->basis2;
        if (StringCmp (iep->prefix, "alignment") == 0) {
          accession_list = iep->accession_list;
          RemoveEmptyAccessionStrings (accession_list);
          ReplaceDatabaseStrings (&accession_list);
        }
      }
    }

    len = StringLen (prefix) + StringLen (speciesies) + StringLen (first) + StringLen (second) + StringLen (accession_list);
    str = MemNew (len + 6);
    if (str != NULL) {
      StringCpy (str, prefix);
      StringCat (str, speciesies);
      StringCat (str, ":");
      StringCat (str, first);
      StringCat (str, ":");
      StringCat (str, second);
      if (StringCmp (prefix, "alignment") == 0 && accession_list != NULL) {
        StringCat (str, ":");
        StringCat (str, accession_list);
      }
      gbq->val = StringSave (str);
      MemFree (str);
    } else {
      gbq->val = StringSave ("?");
    }

    /* do not allow saving of bad qualifier */
    if (convertBadToNote &&
        ValidateInferenceQualifier (gbq->val, FALSE) != VALID_INFERENCE) {
      if (StringNICmp (gbq->val, "similar to ", 11) == 0) {
        len = StringLen ("similar to ") + StringLen (first) + StringLen (second);
        str = MemNew (len + 5);
        if (str != NULL) {
          StringCpy (str, "similar to ");
          if (StringDoesHaveText (first)) {
            StringCat (str, first);
            if (StringDoesHaveText (second)) {
              StringCat (str, ":");
              StringCat (str, second);
            }
          } else if (StringDoesHaveText (second)) {
            StringCat (str, second);
          }
          gbq->val = MemFree (gbq->val);
          gbq->val = StringSave (str);
          MemFree (str);
        }
      }
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("note");
    }

    if (sfp->qual == NULL) {
      sfp->qual = gbq;
    }
    if (lastgbq != NULL) {
      lastgbq->next = gbq;
    }
    lastgbq = gbq;
    accession_list = NULL;
  }
}

static void CleanupInferProc (GraphiC g, VoidPtr data)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            k;

  idp = (InferDialogPtr) data;
  if (idp != NULL) {
    for (k = 0; k < 128; k++) {
      iep = idp->evidence [k];
      InferEvidFree (iep);
      idp->evidence [k] = NULL;
    }
  }
  StdCleanupExtraProc (g, data);
}

static DialoG NewCreateInferenceDialog (
  GrouP prnt
)

{
  GrouP           cts, tbl, g0, g1, g2, g3, g4, g5, p, prompt_grp;
  FonT            fnt;
  Int2            i, hgt, wid;
  InferDialogPtr  idp;

  idp = (InferDialogPtr) MemNew (sizeof (InferDialog));
  if (idp == NULL) return NULL;

  p = HiddenGroup (prnt, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  SetObjectExtra (p, idp, CleanupInferProc);
  idp->dialog = (DialoG) p;
  /*
  idp->todialog = GBQualToInferTable;
  idp->fromdialog = InferTableToGBQual;
  */

  SelectFont (systemFont);
  hgt = LineHeight ();
  inferColFmt [0].pixWidth = MaxAlistWidths (inference_alist) + 5;
  inferColFmt [1].pixWidth = 25 * StringWidth ("X") + 5;
  SelectFont (systemFont);

  wid = 0;
  for (i = 0; i < 2; i++) {
    wid += inferColFmt [i].pixWidth;
  }

  tbl = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (tbl, 10, 5);
  SetGroupMargins (tbl, 5, 5);

  g0 = HiddenGroup (tbl, 15, 0, NULL);
  SetGroupSpacing (g0, 0, 3);
#ifdef WIN_MSWIN
  fnt = systemFont;
#else
  fnt = programFont;
#endif
  /*
  StaticPrompt (g0, "Category", inferColFmt [0].pixWidth, 0, fnt, 'c');
  StaticPrompt (g0, "Explanation", inferColFmt [1].pixWidth, 0, fnt, 'c');
  */

  idp->inferdoc = DocumentPanel (tbl, wid + 2, NUM_INFERENCE_LINES * hgt + 2);
  SetObjectExtra (idp->inferdoc, idp, NULL);
  SetDocCache (idp->inferdoc, NULL, NULL, NULL);
  SetDocNotify (idp->inferdoc, ChangeInferTableSelect);
  idp->numInf = 0;

  for (i = 0; i < NUM_INFERENCE_LINES; i++) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
  }

  cts = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (cts, 10, 10);
  SetGroupMargins (cts, 5, 5);

  g1 = HiddenGroup (cts, -10, 0, NULL);
  SetGroupSpacing (g1, 5, 5);

  StaticPrompt (g1, "Category", 0, popupMenuHeight, programFont, 'l');
  idp->prefix = CreateEnumPopupDialog (g1, TRUE, ChangeInferPrefix, inference_alist, (UIEnum) 0, idp);

  idp->species = CheckBox (g1, "(same species)", ChangeSameSpecies);
  SetObjectExtra (idp->species, idp, NULL);
  Hide (idp->species);

  g2 = HiddenGroup (cts, 0, 0, NULL);
  SetGroupSpacing (g2, 5, 5);

  g3 = HiddenGroup (g2, -3, 0, NULL);
  SetGroupSpacing (g3, 5, 5);

  StaticPrompt (g3, "Database", 0, dialogTextHeight, programFont, 'l');
  idp->database = CreateEnumPopupDialog (g3, TRUE, ChangeInferDatabase, accn_type_alist, (UIEnum) 0, idp);
  idp->other_db_group = HiddenGroup (g3, -4, 0, NULL);
  StaticPrompt (idp->other_db_group, ":", 0, dialogTextHeight, programFont, 'l');
  idp->db_other = DialogText (idp->other_db_group, "", 8, ChangeInferDbOther);
  SetObjectExtra (idp->db_other, idp, NULL);
  Hide (idp->other_db_group);

  StaticPrompt (g3, "Accession", 0, dialogTextHeight, programFont, 'l');
  idp->accession = DialogText (g3, "", 10, ChangeInferAccession);
  SetObjectExtra (idp->accession, idp, NULL);

  idp->inf_accn_group = g3;
  Hide (idp->inf_accn_group);

  g4 = HiddenGroup (g2, -3, 0, NULL);
  SetGroupSpacing (g4, 5, 5);

  StaticPrompt (g4, "Program", 0, dialogTextHeight, programFont, 'l');
  idp->program = CreateEnumPopupDialog (g4, TRUE, ChangeInferProgram, program_alist, (UIEnum) 0, idp);
  idp->other_pr_group = HiddenGroup (g4, -4, 0, NULL);
  StaticPrompt (idp->other_pr_group, ":", 0, dialogTextHeight, programFont, 'l');
  idp->pr_other = DialogText (idp->other_pr_group, "", 8, ChangeInferPrOther);
  SetObjectExtra (idp->pr_other, idp, NULL);
  Hide (idp->other_pr_group);

  StaticPrompt (g4, "Program Version", 0, dialogTextHeight, programFont, 'l');
  idp->version = DialogText (g4, "", 3, ChangeInferVersion);
  SetObjectExtra (idp->version, idp, NULL);

  idp->inf_prog_group = g4;
  Hide (idp->inf_prog_group);

  g5 = HiddenGroup (g2, 2, 0, NULL);
  SetGroupSpacing (g5, 5, 5);

  prompt_grp = HiddenGroup (g5, 0, 0, NULL);
  idp->inf_free_program_prompt = StaticPrompt (prompt_grp, "Program or Database", 0, dialogTextHeight, programFont, 'l');
  idp->accession_list_program_prompt = StaticPrompt (prompt_grp, "Program", 0, dialogTextHeight, programFont, 'l');
  idp->basis1 = DialogText (g5, "", 10, ChangeInferBasis1);
  SetObjectExtra (idp->basis1, idp, NULL);

  prompt_grp = HiddenGroup (g5, 0, 0, NULL);
  idp->inf_free_version_prompt = StaticPrompt (prompt_grp, "Version or Accession", 0, dialogTextHeight, programFont, 'l');
  idp->accession_list_version_prompt = StaticPrompt (prompt_grp, "Version", 0, dialogTextHeight, programFont, 'l');

  idp->basis2 = DialogText (g5, "", 10, ChangeInferBasis2);
  SetObjectExtra (idp->basis2, idp, NULL);

  idp->accession_list_prompt = StaticPrompt (g5, "Accessions", 0, dialogTextHeight, programFont, 'l');
  idp->accession_list = CreateTagListDialogEx3 (g5, 3, 2, 2,
                              accessionlist_types, accessionlist_widths,
                              accessionlist_popups, TRUE, FALSE, AccessionListDataToDialog, AccessionListDialogToData,
                              AccessionListCallbacks, idp,
                              FALSE, FALSE);

  idp->inf_free_group = g5;
  Hide (idp->inf_free_group);

  AlignObjects (ALIGN_CENTER, (HANDLE) tbl, (HANDLE) cts, NULL);

  idp->numInf = 0;
  idp->currItem = 1;
  SetDocHighlight (idp->inferdoc, 1, 1);

  return (DialoG) p;
}


/* ExistingText handling dialog and structures */
typedef struct existingtextdlg 
{
  GrouP pre_app_grp;
  GrouP delim_grp;
} ExistingTextDlgData, PNTR ExistingTextDlgPtr;

static void ChangePreAppIgnoreChoice (GrouP g)
{
  ExistingTextDlgPtr etdp;
  Int4               handle_choice;
  
  etdp = (ExistingTextDlgPtr) GetObjectExtra (g);
  if (etdp == NULL)
  {
    return;
  }
  
  handle_choice = GetValue (etdp->pre_app_grp);
  if (handle_choice == 1 || handle_choice == 2)
  {
    Enable (etdp->delim_grp);
  }
  else
  {
    Disable (etdp->delim_grp);
  }
}

extern ExistingTextPtr GetExistingTextHandlerInfo (Int4 num_found, Boolean non_text)
{
  WindoW                w;
  GrouP                 h, c;
  ExistingTextDlgData   etdd;
  ButtoN                b;
  ModalAcceptCancelData acd;
  ExistingTextPtr       etp;
  Char                  txt [128];
  MsgAnswer             ans;
  PrompT                ppt;
  Int4                  handle_choice;

  if (num_found <= 0)
  {
    return NULL;
  }
  
  sprintf (txt, "%d affected fields already contain a value.  Do you wish to overwrite existing text?",
           num_found);
  ans = Message (MSG_YNC, txt, 0, dialogTextHeight, systemFont, 'l');
  if (ans == ANS_CANCEL)
  {
    etp = (ExistingTextPtr) MemNew (sizeof (ExistingTextData));
    etp->existing_text_choice = eExistingTextChoiceCancel;
    return etp;
  }
  else if (ans == ANS_YES)
  {
    etp = (ExistingTextPtr) MemNew (sizeof (ExistingTextData));
    etp->existing_text_choice = eExistingTextChoiceReplaceOld;
    return etp;
  }
    
  w = MovableModalWindow(-20, -13, -10, -10, "How to Add New Text", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  etdd.pre_app_grp = HiddenGroup (h, 0, 3, ChangePreAppIgnoreChoice);
  SetGroupSpacing (etdd.pre_app_grp, 10, 10);
  RadioButton (etdd.pre_app_grp, "Append");
  RadioButton (etdd.pre_app_grp, "Prefix");
  RadioButton (etdd.pre_app_grp, "Ignore new text");
  SetValue (etdd.pre_app_grp, 1);
  SetObjectExtra (etdd.pre_app_grp, &etdd, NULL);
  
  ppt = StaticPrompt (h, "Separate new text and old text with", 
                      0, dialogTextHeight, programFont, 'c');
  etdd.delim_grp = HiddenGroup (h, 0, 4, NULL);
  SetGroupSpacing (etdd.delim_grp, 10, 10);
  RadioButton (etdd.delim_grp, "Semicolon");
  RadioButton (etdd.delim_grp, "Space");
  RadioButton (etdd.delim_grp, "Colon");
  RadioButton (etdd.delim_grp, "Do not separate");
  SetValue (etdd.delim_grp, 1);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) etdd.pre_app_grp,
                              (HANDLE) ppt, 
                              (HANDLE) etdd.delim_grp, 
                              (HANDLE) c, 
                              NULL);
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
  etp = (ExistingTextPtr) MemNew (sizeof (ExistingTextData));
  if (acd.cancelled)
  {
    etp->existing_text_choice = eExistingTextChoiceCancel;
  }
  else
  {
    handle_choice = GetValue (etdd.pre_app_grp);
    if (handle_choice == 1)
    {
      switch (GetValue (etdd.delim_grp))
      {
        case 1:
          etp->existing_text_choice = eExistingTextChoiceAppendSemi;
          break;
        case 2:
          etp->existing_text_choice = eExistingTextChoiceAppendSpace;
          break;
        case 3:
          etp->existing_text_choice = eExistingTextChoiceAppendColon;
          break;
        case 4:
          etp->existing_text_choice = eExistingTextChoiceAppendNone;
          break;
      }
    }
    else if (handle_choice == 2)
    {
      switch (GetValue (etdd.delim_grp))
      {
        case 1:
          etp->existing_text_choice = eExistingTextChoicePrefixSemi;
          break;
        case 2:
          etp->existing_text_choice = eExistingTextChoicePrefixSpace;
          break;
        case 3:
          etp->existing_text_choice = eExistingTextChoicePrefixColon;
          break;
        case 4:
          etp->existing_text_choice = eExistingTextChoicePrefixNone;
          break;
      }
    }
    else
    {
      etp->existing_text_choice = eExistingTextChoiceLeaveOld;
    }
  }
  Remove (w);
  return etp;
}

extern CharPtr HandleExistingText (CharPtr existing_text, CharPtr new_text, ExistingTextPtr etp)
{
  CharPtr rstring = NULL;
  Int4    len;
  
  if (StringHasNoText (existing_text) || etp == NULL)
  {
    MemFree (existing_text);
    return new_text;
  }
  switch (etp->existing_text_choice)
  {
    case eExistingTextChoiceReplaceOld:
      /* replace current text with new text */
      MemFree (existing_text);
      rstring = new_text;
      break;
    case eExistingTextChoiceLeaveOld:
      /* do not change current text */
      MemFree (new_text);
      rstring = existing_text;
      break;
    case eExistingTextChoiceAppendSemi:
      /* Append new text to current text, separated by semicolon */
      len = StringLen (new_text) + StringLen (existing_text) + 4;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, existing_text);
        StringCat (rstring, "; ");
        StringCat (rstring, new_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoiceAppendSpace:
      /* Append new text to current text, separated by space */
      len = StringLen (new_text) + StringLen (existing_text) + 3;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, existing_text);
        StringCat (rstring, " ");
        StringCat (rstring, new_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoiceAppendColon:
      /* Append new text to current text, separated by colon */
      len = StringLen (new_text) + StringLen (existing_text) + 4;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, existing_text);
        StringCat (rstring, ": ");
        StringCat (rstring, new_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoiceAppendNone:
      /* Append new text to current text, no delimiter */
      len = StringLen (new_text) + StringLen (existing_text) + 1;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, existing_text);
        StringCat (rstring, new_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoicePrefixSemi:
      /* Prepend new text to current text, separated by semicolon */
      len = StringLen (new_text) + StringLen (existing_text) + 4;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, new_text);
        StringCat (rstring, "; ");
        StringCat (rstring, existing_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoicePrefixSpace:
      /* Prepend new text to current text, separated by space */
      len = StringLen (new_text) + StringLen (existing_text) + 3;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, new_text);
        StringCat (rstring, " ");
        StringCat (rstring, existing_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoicePrefixColon:
      /* Prepend new text to current text, separated by colon */
      len = StringLen (new_text) + StringLen (existing_text) + 4;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, new_text);
        StringCat (rstring, ": ");
        StringCat (rstring, existing_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;
    case eExistingTextChoicePrefixNone:
      /* prefix current text with new text */
      len = StringLen (new_text) + StringLen (existing_text) + 1;
      rstring = MemNew (len);
      if (rstring != NULL) {
        StringCpy (rstring, new_text);
        StringCat (rstring, existing_text);
        MemFree (new_text);
        MemFree (existing_text);
      }
      break;    
  }
  return rstring;
}


/* Move EditApply data and dialog here */
extern EditApplyPtr EditApplyFree (EditApplyPtr eap)
{
  if (eap != NULL)
  {
    eap->find_txt = MemFree (eap->find_txt);
    eap->repl_txt = MemFree (eap->repl_txt);
    eap->apply_txt = MemFree (eap->apply_txt);
    eap = MemFree (eap);
  }
  return eap;
}


extern EditApplyPtr EditApplyNew (void)
{
  EditApplyPtr eap;

  eap = (EditApplyPtr) MemNew (sizeof (EditApplyData));
  eap->find_location = EditApplyFindLocation_anywhere;
  return eap;
}


typedef struct EditApplydlg
{
  DIALOG_MESSAGE_BLOCK
  TexT           find_txt;
  TexT           repl_txt;
  TexT           apply_txt;

  DialoG         find_dlg;
  DialoG         repl_dlg;
  DialoG         apply_dlg;

  Int4           action_choice;
  GrouP          location_choice;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} EditApplyDlgData, PNTR EditApplyDlgPtr;

static void ResetEditApplyDlg (EditApplyDlgPtr dlg)
{
  if (dlg != NULL)
  {
    if (dlg->find_txt != NULL)
    {
      SetTitle (dlg->find_txt, "");
    }
    if (dlg->repl_txt != NULL)
    {
      SetTitle (dlg->repl_txt, "");
    }
    if (dlg->apply_txt != NULL)
    {
      SetTitle (dlg->apply_txt, "");
    }
    if (dlg->location_choice != NULL) {
      SetValue (dlg->location_choice, EditApplyFindLocation_anywhere);
    }

    PointerToDialog (dlg->find_dlg, NULL);
    PointerToDialog (dlg->repl_dlg, NULL);
    PointerToDialog (dlg->apply_dlg, NULL);

  }
}

static void EditApplyDialogChangeText (TexT t)
{
  EditApplyDlgPtr dlg;

  dlg = (EditApplyDlgPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static void EditApplyToDialog (DialoG d, Pointer userdata)
{
  EditApplyDlgPtr dlg;
  EditApplyPtr    data;
  ValNode         vn;
  
  dlg = (EditApplyDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetEditApplyDlg (dlg);
  data = (EditApplyPtr) userdata;

  vn.next = NULL;
  vn.choice = 0;

  if (data != NULL)
  {
    if (!StringHasNoText (data->find_txt))
    {
      if (dlg->find_txt != NULL)
      {
        SetTitle (dlg->find_txt, data->find_txt);
      }
      else if (dlg->find_dlg != NULL)
      {
        vn.data.ptrvalue = data->find_txt;
        PointerToDialog (dlg->find_dlg, &vn);
      }
    }

    if (!StringHasNoText (data->repl_txt))
    {
      if (dlg->repl_txt != NULL) 
      {
        SetTitle (dlg->repl_txt, data->repl_txt);
      }
      else if (dlg->repl_dlg != NULL)
      {
        vn.data.ptrvalue = data->repl_txt;
        PointerToDialog (dlg->repl_dlg, &vn);
      }
    }

    if (!StringHasNoText (data->apply_txt))
    {
      if (dlg->apply_txt != NULL)
      {
        SetTitle (dlg->apply_txt, data->apply_txt);
      }
      else if (dlg->apply_dlg != NULL)
      {
        vn.data.ptrvalue = data->apply_txt;
        PointerToDialog (dlg->apply_dlg, &vn);
      }
    }

    if (dlg->location_choice != NULL) {
      SetValue (dlg->location_choice, data->find_location);
    }
  }
}

static Pointer DialogToEditApply (DialoG d)
{
  EditApplyDlgPtr dlg;
  EditApplyPtr    data;
  ValNodePtr      vnp;
  
  dlg = (EditApplyDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  data = (EditApplyPtr) MemNew (sizeof (EditApplyData));
  if (data != NULL)
  {
    if (dlg->find_txt != NULL)
    {
      data->find_txt = JustSaveStringFromText (dlg->find_txt);
    }
    else if (dlg->find_dlg != NULL)
    {
      vnp = (ValNodePtr) DialogToPointer (dlg->find_dlg);
      if (vnp != NULL)
      {
        data->find_txt = StringSave (vnp->data.ptrvalue);
      }
      vnp = ValNodeFreeData (vnp);
    }

    if (dlg->repl_txt != NULL)
    {
      data->repl_txt = JustSaveStringFromText (dlg->repl_txt);
    }
    else if (dlg->repl_dlg != NULL)
    {
      vnp = (ValNodePtr) DialogToPointer (dlg->repl_dlg);
      if (vnp != NULL)
      {
        data->repl_txt = StringSave (vnp->data.ptrvalue);
      }
      vnp = ValNodeFreeData (vnp);
    }

    if (dlg->apply_txt != NULL)
    {
      data->apply_txt = JustSaveStringFromText (dlg->apply_txt);
    }
    else if (dlg->apply_dlg != NULL)
    {
      vnp = (ValNodePtr) DialogToPointer (dlg->apply_dlg);
      if (vnp != NULL)
      {
        data->apply_txt = StringSave (vnp->data.ptrvalue);
      }
      vnp = ValNodeFreeData (vnp);
    }

    if (dlg->location_choice != NULL) {
      data->find_location = (EditApplyFindLocation) GetValue (dlg->location_choice);
    } else {
      data->find_location = EditApplyFindLocation_anywhere;
    }
  }
  return data;
}


static void EditApplyMessage (DialoG d, Int2 mssg)

{
  EditApplyDlgPtr  dlg;

  dlg = (EditApplyDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetEditApplyDlg (dlg);
        break;
      case VIB_MSG_ENTER :
        if (dlg->find_txt != NULL)
        {
          Select (dlg->find_txt);
        }
        else if (dlg->apply_txt != NULL)
        {
          Select (dlg->apply_txt);
        }
        else if (dlg->find_dlg != NULL)
        {
          Select (dlg->find_dlg);
        }
        else if (dlg->apply_dlg != NULL)
        {
          Select (dlg->apply_dlg);
        }
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestEditApply (DialoG d)
{
  EditApplyDlgPtr dlg;
  ValNodePtr      total_err_list = NULL, err_list;
  
  dlg = (EditApplyDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }

  if (dlg->action_choice == eEditApplyChoice_Apply)
  {
    if (dlg->apply_dlg == NULL) 
    {
      if (TextHasNoText (dlg->apply_txt))
      {
        ValNodeAddPointer (&total_err_list, 0, StringSave ("apply text"));
      }
    }
    else
    {
      total_err_list = TestDialog (dlg->apply_dlg);
    }
  }
  else if (dlg->action_choice == eEditApplyChoice_Edit)
  {
    if (dlg->find_dlg == NULL)
    {
      if (TextHasNoText (dlg->find_txt))
      {
        ValNodeAddPointer (&total_err_list, 0, StringSave ("find text"));
      }
    }
    else 
    {
      total_err_list = TestDialog (dlg->find_dlg);
      err_list = TestDialog (dlg->repl_dlg);
      ValNodeLink (&total_err_list, err_list);
    }
  }
  return total_err_list;
}

static void EditApplyDialogCopy (ButtoN b)
{
  EditApplyDlgPtr dlg;
  CharPtr         str = NULL;

  dlg = (EditApplyDlgPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  str = JustSaveStringFromText (dlg->find_txt);
  SetTitle (dlg->repl_txt, str);
  str = MemFree (str);
}

extern DialoG EditApplyDialog 
(GrouP                    h,
 Int4                     action_choice, 
 CharPtr                  apply_label,
 ValNodePtr               choice_list,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  EditApplyDlgPtr dlg;
  GrouP           p, p1;
  ButtoN          b;
  ValNodePtr      cpy;
  
  dlg = (EditApplyDlgPtr) MemNew (sizeof (EditApplyDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = EditApplyToDialog;
  dlg->fromdialog = DialogToEditApply;
  dlg->dialogmessage = EditApplyMessage;
  dlg->testdialog = TestEditApply;
  dlg->action_choice = action_choice;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (choice_list == NULL)
  {
    p1 = HiddenGroup (p, 3, 0, NULL);
    SetGroupSpacing (p1, 10, 10);

    if (action_choice == eEditApplyChoice_Apply)
    {
      StaticPrompt (p1, apply_label, 0, dialogTextHeight, systemFont, 'r');
      dlg->apply_txt = DialogText (p1, "", 20, EditApplyDialogChangeText);
      SetObjectExtra (dlg->apply_txt, dlg, NULL);
      dlg->location_choice = NULL;
    }
    else if (action_choice == eEditApplyChoice_Edit)
    {
      StaticPrompt (p1, "Find", 0, dialogTextHeight, systemFont, 'r');
      dlg->find_txt = DialogText (p1, "", 18, EditApplyDialogChangeText);
      SetObjectExtra (dlg->find_txt, dlg, NULL);
      b = PushButton (p1, "Copy", EditApplyDialogCopy);
      SetObjectExtra (b, dlg, NULL);
      Hide (b);
      StaticPrompt (p1, "Replace", 0, dialogTextHeight, systemFont, 'r');
      dlg->repl_txt = DialogText (p1, "", 18, EditApplyDialogChangeText);
      SetObjectExtra (dlg->repl_txt, dlg, NULL);
      b = PushButton (p1, "Copy", EditApplyDialogCopy);
      SetObjectExtra (b, dlg, NULL);

      dlg->location_choice = HiddenGroup (p, 3, 0, NULL);
      RadioButton (dlg->location_choice, "Anywhere in field");
      RadioButton (dlg->location_choice, "At the beginning of the field");
      RadioButton (dlg->location_choice, "At the end of the field");
      SetValue (dlg->location_choice, EditApplyFindLocation_anywhere);
    }
  }
  else
  {
    if (action_choice == eEditApplyChoice_Apply)
    {
      p1 = HiddenGroup (p, 1, 0, NULL);
      SetGroupSpacing (p1, 10, 10);
      if (!StringHasNoText (apply_label))
      {
        StaticPrompt (p1, apply_label, 0, dialogTextHeight, systemFont, 'r');
      }
      cpy = ValNodeDupStringList (choice_list);
      dlg->apply_dlg = ValNodeSelectionDialog (p1, cpy, 6,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               apply_label,
                                               dlg->change_notify, dlg->change_userdata, FALSE);
    }
    else if (action_choice == eEditApplyChoice_Edit)
    {
      p1 = HiddenGroup (p, 2, 0, NULL);
      SetGroupSpacing (p1, 10, 10);
      StaticPrompt (p1, "From", 0, dialogTextHeight, systemFont, 'r');
      StaticPrompt (p1, "To", 0, dialogTextHeight, systemFont, 'r');

      cpy = ValNodeDupStringList (choice_list);
      dlg->find_dlg = ValNodeSelectionDialog (p1, cpy, 6,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "Original",
                                               dlg->change_notify, dlg->change_userdata, FALSE);
      cpy = ValNodeDupStringList (choice_list);
      dlg->repl_dlg = ValNodeSelectionDialog (p1, cpy, 6,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "New",
                                               dlg->change_notify, dlg->change_userdata, FALSE);

    }
    dlg->location_choice = NULL;
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) p1, (HANDLE) dlg->location_choice, NULL);

  return (DialoG) p;
}


/* Global Inference Editor Dialog */

extern InferenceParsePtr ParseInferenceText (CharPtr inference)
{
  CharPtr cp1, cp2;
  Int4    len;
  InferenceParsePtr ipp;

  if (StringHasNoText (inference))
  {
    return NULL;
  }

  cp1 = StringChr (inference, ':');
  if (cp1 == NULL)  
  {
    return NULL;
  }
  cp2 = StringChr (cp1 + 1, ':');
  if (cp2 == NULL)
  {
    return NULL;
  }

  ipp = (InferenceParsePtr) MemNew (sizeof (InferenceParseData));
  ipp->second_field = StringSave (cp2 + 1);

  len = cp2 - cp1;
  ipp->first_field = (CharPtr) MemNew (sizeof (Char) * len);
  StringNCpy (ipp->first_field, cp1 + 1, len - 1);
  ipp->first_field[len - 1] = 0;

  /* look for same species */
  cp2 = StringISearch (inference, "(same species)");
  if (cp2 != NULL && cp2 < cp1)
  {
    ipp->same_species = TRUE;
    cp1 = cp2;
  } 
  else
  {
    ipp->same_species = FALSE;
  }
  len = cp1 - inference + 1;
  ipp->category = (CharPtr) MemNew (sizeof (Char) * len);
  StringNCpy (ipp->category, inference, len - 1);
  ipp->category[len - 1] = 0;
  TrimSpacesAroundString (ipp->category);
  TrimSpacesAroundString (ipp->first_field);
  TrimSpacesAroundString (ipp->second_field);

  return ipp;  
}


extern CharPtr InferenceTextFromStruct (InferenceParsePtr ipp)
{
  Int4 len;
  CharPtr inference = NULL;
  CharPtr same_sp = " (same species)";

  if (ipp == NULL) return NULL;

  len = StringLen (ipp->category) + StringLen (ipp->first_field) + StringLen (ipp->second_field)
        + 3;
  if (ipp->same_species)
  {
    len += StringLen (same_sp);
  }

  inference = (CharPtr) MemNew (sizeof (Char) * len);
  sprintf (inference, "%s%s:%s:%s", ipp->category == NULL ? "" : ipp->category,
                                    ipp->same_species ? same_sp : "",
                                    ipp->first_field == NULL ? "" : ipp->first_field,
                                    ipp->second_field == NULL ? "" : ipp->second_field);
  return inference;
}


typedef enum {
  eInferenceCategorySimilarTo = 0,
  eInferenceCategoryProgram,
  eInferenceCategoryAbInitio,
  eNumInferenceCategories } EInferenceCategoryType;


static Int4 GetCategoryTypeFromNum (Int4 num)
{
  if (num > 0 && num < 8) 
  {
    return eInferenceCategorySimilarTo;
  }
  else if (num > 7 && num < 11)
  {
    return eInferenceCategoryProgram;
  }
  else if (num == 11)
  {
    return eInferenceCategoryAbInitio;
  }
  else
  {
    return -1;
  }
}


static Int4 GetCategoryNumFromName (CharPtr category)
{
  Int4 i;

  if (StringHasNoText (category)) 
  {
    return -1;
  }

  for (i = 0; inference_alist[i].name != NULL; i++)
  {
    if (StringICmp (inference_alist[i].name, category) == 0)
    {
      return i;
    }
  }
  return -1;
}


extern InferenceFieldEditPtr InferenceFieldEditFree (InferenceFieldEditPtr ifep)
{
  if (ifep != NULL)
  {
    ifep->edit_apply = EditApplyFree (ifep->edit_apply);
    ifep = MemFree (ifep);
  }
  return ifep;
}

extern InferenceEditPtr InferenceEditFree (InferenceEditPtr iep)
{
  if (iep != NULL)
  {
    iep->field_edit = InferenceFieldEditFree (iep->field_edit);
    iep = MemFree (iep);
  }
  return iep;
}

typedef struct inferencefieldeditdialog
{
  DIALOG_MESSAGE_BLOCK

  PopuP  field_category;
  PopuP  field_list[eNumInferenceCategories];
  DialoG field_editors[eNumInferenceCategories * 2];
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;

} InferenceFieldEditDialogData, PNTR InferenceFieldEditDialogPtr;


static Pointer InferenceFieldEditDataFromDialog (DialoG d)
{
  InferenceFieldEditDialogPtr dlg;
  UIEnum                      val;
  Int4                        i, j;
  InferenceFieldEditPtr       data = NULL;

  dlg = (InferenceFieldEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (GetEnumPopup (dlg->field_category, inference_alist, &val)
      && val > 0
      && (i = GetCategoryTypeFromNum (val)) > -1
      && (j = GetValue (dlg->field_list[i])) > 0
      && j < 3)
  {
    data = (InferenceFieldEditPtr) MemNew (sizeof (InferenceFieldEditData));
    data->field_category = inferencePrefix[val];
    data->field_choice = j - 1;
    data->edit_apply = DialogToPointer (dlg->field_editors[2 * i + j - 1]);
  }
  return data;
}


static void ChangeInferenceFieldChoice (PopuP p)
{
  InferenceFieldEditDialogPtr dlg;
  UIEnum                      val;
  Int4                        i, j;

  dlg = (InferenceFieldEditDialogPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  for (i = eInferenceCategorySimilarTo;
       i <= eInferenceCategoryAbInitio;
       i++)
  {
    Hide (dlg->field_list[i]);
    Hide (dlg->field_editors[2 * i]);
    Hide (dlg->field_editors[2 * i + 1]);
  }

  if (GetEnumPopup (dlg->field_category, inference_alist, &val)
      && val > 0
      && (i = GetCategoryTypeFromNum (val)) > -1)
  {
    Show (dlg->field_list[i]);
    j = GetValue (dlg->field_list[i]);
    if (j > 0 && j < 3)
    {
      Show (dlg->field_editors[2 * i + j - 1]);
    }
  }

  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static ValNodePtr MakeValNodeListFromEnum ( EnumFieldAssocPtr al)
{
  EnumFieldAssocPtr efap;
  ValNodePtr        list;

  efap = al;
  list = NULL;
  while (efap->name != NULL)
  {
    ValNodeAddStr (&list, efap->value, StringSave (efap->name));
    efap ++;
  }
  return list;
}


static DialoG 
CreateInferenceFieldEditApplyDialog 
(GrouP                h,
 Int4                 action_choice,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  InferenceFieldEditDialogPtr dlg;
  GrouP                       p, k;
  Int4                        i;
  ValNodePtr                  choice_list;

  dlg = (InferenceFieldEditDialogPtr) MemNew (sizeof (InferenceFieldEditDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (h, 3, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = InferenceFieldEditDataFromDialog;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  StaticPrompt (p, "Category", 0, popupMenuHeight, programFont, 'c');
  StaticPrompt (p, "Field", 0, popupMenuHeight, programFont, 'c');
  if (action_choice == eEditApplyChoice_Apply) 
  {
    StaticPrompt (p, "New Value", 0, popupMenuHeight, programFont, 'c');
  }
  else
  {
    StaticPrompt (p, "Convert", 0, popupMenuHeight, programFont, 'c');
  }

  dlg->field_category = PopupList (p, TRUE, ChangeInferenceFieldChoice);
  SetObjectExtra (dlg->field_category, dlg, NULL);
  InitEnumPopup (dlg->field_category, inference_alist, NULL);


  k = HiddenGroup (p, 0, 0, NULL);
  dlg->field_list[eInferenceCategorySimilarTo] = PopupList (k, TRUE, ChangeInferenceFieldChoice);
  SetObjectExtra (dlg->field_list[eInferenceCategorySimilarTo], dlg, NULL);
  PopupItem (dlg->field_list[eInferenceCategorySimilarTo], "Database");
  PopupItem (dlg->field_list[eInferenceCategorySimilarTo], "Accession");
  
  dlg->field_list[eInferenceCategoryProgram] = PopupList (k, TRUE, ChangeInferenceFieldChoice);
  SetObjectExtra (dlg->field_list[eInferenceCategoryProgram], dlg, NULL);
  PopupItem (dlg->field_list[eInferenceCategoryProgram], "Program or Database");
  PopupItem (dlg->field_list[eInferenceCategoryProgram], "Version or Accession");
  
  dlg->field_list[eInferenceCategoryAbInitio] = PopupList (k, TRUE, ChangeInferenceFieldChoice);
  SetObjectExtra (dlg->field_list[eInferenceCategoryAbInitio], dlg, NULL);
  PopupItem (dlg->field_list[eInferenceCategoryAbInitio], "Program");
  PopupItem (dlg->field_list[eInferenceCategoryAbInitio], "Program Version");

  k = HiddenGroup (p, 0, 0, NULL);
  i = 0;
  choice_list = MakeValNodeListFromEnum (accn_type_alist);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", choice_list, change_notify, change_userdata);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", NULL, change_notify, change_userdata);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", NULL, change_notify, change_userdata);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", NULL, change_notify, change_userdata);
  choice_list = MakeValNodeListFromEnum (program_alist);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", choice_list, change_notify, change_userdata);
  dlg->field_editors[i++] = EditApplyDialog (k, action_choice, "", NULL, change_notify, change_userdata);

  ChangeInferenceFieldChoice (dlg->field_category);
  return (DialoG) p;
}


typedef struct inferenceeditdialog {
  DIALOG_MESSAGE_BLOCK
  PopuP action;
  GrouP action_pages[eNumInferenceEditActions];
  PopuP category_from;
  PopuP category_to;
  
  DialoG apply_field;
  DialoG edit_field;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;

} InferenceEditDialogData, PNTR InferenceEditDialogPtr;


static Pointer InferenceEditDataFromDialog (DialoG d)
{
  InferenceEditDialogPtr dlg;
  InferenceEditPtr       data;
  Int4                   i;

  dlg = (InferenceEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  data = (InferenceEditPtr) MemNew (sizeof (InferenceEditData));
  if (data == NULL) return NULL;

  data->action = (EInferenceEditAction) (GetValue (dlg->action) - 1);

  switch (data->action)
  {
    case eInferenceRemove:
      /* no other data needed */
      break;
    case eInferenceEditCategory:
      i = GetValue (dlg->category_from);
      if (i <= 1) 
      {
        data->category_from = NULL;
      }
      else
      {
        data->category_from = inferencePrefix[i - 2];
      }
      i = GetValue (dlg->category_to);
      if (i < 1)
      {
        data->category_to = NULL;
      }
      else
      {
        data->category_to = inferencePrefix[i - 1];
      }
      break;
    case eInferenceApplyCategoryFields:
     data->field_edit = DialogToPointer (dlg->apply_field);
     break;
    case eInferenceEditCategoryFields:
     data->field_edit = DialogToPointer (dlg->edit_field);
      break;
    default:
      break;
  }

  return data;
}


static void ChangeInferenceEditAction (PopuP p)
{
  InferenceEditDialogPtr dlg;
  Int4                   i;

  dlg = (InferenceEditDialogPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  ResetClip();
  for (i = 0; i < eNumInferenceEditActions; i++)
  {
    Hide (dlg->action_pages[i]);
  }

  i = GetValue (dlg->action);
  if (i > 0 && i <= eNumInferenceEditActions)
  {
    Show (dlg->action_pages[i - 1]);
  }

  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeInferenceCategoryChoice (PopuP p)
{
  InferenceEditDialogPtr dlg;

  dlg = (InferenceEditDialogPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static ValNodePtr TestInferenceEditDialog (DialoG d)
{
  ValNodePtr             err_list = NULL;
  InferenceEditPtr       iep;

  iep = DialogToPointer (d);
  if (iep == NULL)
  {
    ValNodeAddPointer (&err_list, 0, "no values");
  }
  else
  {
    switch (iep->action)
    {
      case eInferenceRemove:
        /* nothing to check */
        break;
      case eInferenceEditCategory:
        if (StringHasNoText (iep->category_from))
        {
          ValNodeAddPointer (&err_list, 0, "missing category from");
        }
        if (StringHasNoText (iep->category_to))
        {
          ValNodeAddPointer (&err_list, 0, "missing category to");
        }
        break;
      case eInferenceApplyCategoryFields:
      case eInferenceEditCategoryFields:
        if (iep->field_edit == NULL)
        {
          ValNodeAddPointer (&err_list, 0, "missing edit data");
        }
        break;
      default:
        break;
    }
  }
  iep = InferenceEditFree (iep);

  return err_list;
}


extern DialoG CreateInferenceEditDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  InferenceEditDialogPtr dlg;
  GrouP                  p, g;
  Int4                   i;
  Nlm_EnumFieldAssocPtr  eap;

  dlg = (InferenceEditDialogPtr) MemNew (sizeof (InferenceEditDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = InferenceEditDataFromDialog;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestInferenceEditDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->action = PopupList (p, TRUE, ChangeInferenceEditAction);
  SetObjectExtra (dlg->action, dlg, NULL);
  PopupItem (dlg->action, "Remove");
  PopupItem (dlg->action, "Change Category");
  PopupItem (dlg->action, "Apply Category Fields");
  PopupItem (dlg->action, "Edit Category Fields");
  SetValue (dlg->action, eInferenceRemove + 1);

  g = HiddenGroup (p, 0, 0, NULL);
  i = 0;

  /* remove group */
  dlg->action_pages[i] = HiddenGroup (g, -1, 0, NULL);
  StaticPrompt (dlg->action_pages[i], "Hit Accept to Remove Inferences", 0, popupMenuHeight, programFont, 'c');
  i++;

  /* edit category group */
  dlg->action_pages[i] = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->action_pages[i], "Original Category", 0, popupMenuHeight, programFont, 'c');
  dlg->category_from = PopupList (dlg->action_pages[i], TRUE, ChangeInferenceCategoryChoice);
  SetObjectExtra (dlg->category_from, dlg, NULL);
  PopupItem (dlg->category_from, "Any");
  eap = inference_alist;
  while (eap->name != NULL) {
    PopupItem (dlg->category_from, eap->name);
    eap++;
  }
  
  StaticPrompt (dlg->action_pages[i], "New Category", 0, popupMenuHeight, programFont, 'c');
  dlg->category_to = PopupList (dlg->action_pages[i], TRUE, ChangeInferenceCategoryChoice);
  SetObjectExtra (dlg->category_to, dlg, NULL);
  InitEnumPopup (dlg->category_to, inference_alist, NULL);
  i++;

  dlg->action_pages[i] = HiddenGroup (g, 0, 0, NULL);
  dlg->apply_field = CreateInferenceFieldEditApplyDialog (dlg->action_pages[i], eEditApplyChoice_Apply, change_notify, change_userdata);
  i++;
  
  dlg->action_pages[i] = HiddenGroup (g, 0, 0, NULL);
  dlg->edit_field = CreateInferenceFieldEditApplyDialog (dlg->action_pages[i], eEditApplyChoice_Edit, change_notify, change_userdata);
  i++;
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->action_pages[0],
                              (HANDLE) dlg->action_pages[1],
                              (HANDLE) dlg->action_pages[2],
                              (HANDLE) dlg->action_pages[3],
                              NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->action, (HANDLE) g, NULL);

  ChangeInferenceEditAction (dlg->action);

  return (DialoG) p;
}


/* This section of code is for handling ClickableLists */

typedef struct clickableitemlist
{
  DIALOG_MESSAGE_BLOCK
  DoC doc;

  Nlm_ParData par_fmt;
  Nlm_ColData col_fmt [4];

  ClickableCallback single_click_callback;
  ClickableCallback double_click_callback;
  Pointer           click_callback_data;
  GetClickableItemText get_item_text;
  ValNodePtr           item_list;
  Int2                 selected;
} ClickableItemListDlgData, PNTR ClickableItemListDlgPtr;

static void PointerToClickableItemListDlg (DialoG d, Pointer data)
{
  ClickableItemListDlgPtr dlg;
  ValNodePtr              vnp;
  CharPtr                 row_text;
  Int2                 numItems;
  RecT                 r;

  dlg = (ClickableItemListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Reset (dlg->doc);
  
  if (dlg->get_item_text == NULL)
  {
    return;
  }

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  
  dlg->col_fmt[0].pixWidth = 5 * stdCharWidth;
  dlg->col_fmt[1].pixWidth = (r.right - r.left - dlg->col_fmt[0].pixWidth) / 3;
  dlg->col_fmt[2].pixWidth = (r.right - r.left - dlg->col_fmt[0].pixWidth) / 3;
  dlg->col_fmt[3].pixWidth = (r.right - r.left - dlg->col_fmt[0].pixWidth) / 3;

  dlg->item_list = ValNodeFree (dlg->item_list);
  dlg->item_list = (ValNodePtr) data; 
 
  if (dlg->item_list == NULL)
  {
    AppendText (dlg->doc, "No items listed", NULL, NULL, programFont);
  } else {
    for (vnp = dlg->item_list; vnp != NULL; vnp = vnp->next) {
      row_text = dlg->get_item_text (vnp);
      if (row_text != NULL)
      {
        if (vnp->choice == OBJ_SEQFEAT)
        {
          AppendText (dlg->doc, row_text, &(dlg->par_fmt), dlg->col_fmt, programFont);
        }
        else
        {
          AppendText (dlg->doc, row_text, &(dlg->par_fmt), NULL, programFont);
        }
        row_text = MemFree (row_text);
      }
    }
  }
  
  GetDocParams (dlg->doc, &numItems, NULL);
  UpdateDocument (dlg->doc, 0, numItems);  
}


static void ClickClickableItemList (DoC d, PoinT pt)

{
  Int2             item, last_selected, numItems;
  Int2             row, i;
  ClickableItemListDlgPtr dlg;
  ValNodePtr       vnp;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {  
      i = 1;
      vnp = dlg->item_list;
      while (i < item && vnp != NULL) {
        i++;
        vnp = vnp->next;
      }
      if (vnp != NULL) {
        last_selected = dlg->selected;
        dlg->selected = item;
        
        if (item != last_selected)
        {
          GetDocParams (d, &numItems, NULL);
          UpdateDocument (d, 0, numItems);
        }
    
        if (dblClick)
        {
          if (dlg->double_click_callback != NULL) {
            (dlg->double_click_callback) (vnp, dlg->click_callback_data);
          }
        } else {
          if (dlg->single_click_callback != NULL) {
            (dlg->single_click_callback) (vnp, dlg->click_callback_data);
          }
        }          
      }
    }
  }
}


static void DrawClickableItemList (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  ClickableItemListDlgPtr dlg;
  RecT             rct;

  dlg = (ClickableItemListDlgPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    /* draw selection */
    if (item == dlg->selected) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
  }
}


static void CleanupClickableItemListDlg (GraphiC g, VoidPtr data)
{
  ClickableItemListDlgPtr dlg;

  dlg = (ClickableItemListDlgPtr) data;
  if (dlg != NULL) {
    dlg->item_list = ValNodeFree (dlg->item_list);
  } 
  StdCleanupExtraProc (g, data);
}

static DialoG 
ClickableItemListDialog 
(GrouP h,
 Int4 width,
 GetClickableItemText get_item_text,
 ClickableCallback single_click_callback,
 ClickableCallback double_click_callback,
 Pointer click_callback_data)
{
  ClickableItemListDlgPtr dlg;
  GrouP                p;
  Int4 i;

  dlg = (ClickableItemListDlgPtr) MemNew (sizeof (ClickableItemListDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupClickableItemListDlg);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToClickableItemListDlg;

  dlg->get_item_text = get_item_text;
  dlg->single_click_callback = single_click_callback;
  dlg->double_click_callback = double_click_callback;
  dlg->click_callback_data = click_callback_data;
  /* initialize paragraph format */
  dlg->par_fmt.openSpace = FALSE;
  dlg->par_fmt.keepWithNext = FALSE;
  dlg->par_fmt.keepTogether = FALSE;
  dlg->par_fmt.newPage = FALSE;
  dlg->par_fmt.tabStops = FALSE;
  dlg->par_fmt.minLines = 0;
  dlg->par_fmt.minHeight = 0;

  /* initialize column format */
  for (i = 0; i < 4; i++) {
    dlg->col_fmt[i].pixWidth = 0;
    dlg->col_fmt[i].pixInset = 0;
    dlg->col_fmt[i].charWidth = 10;
    dlg->col_fmt[i].charInset = 0;
    dlg->col_fmt[i].font = NULL;
    dlg->col_fmt[i].just = 'l';
    dlg->col_fmt[i].wrap = 1;
    dlg->col_fmt[i].bar = 0;
    dlg->col_fmt[i].underline = 0;
    dlg->col_fmt[i].left = 0;
    dlg->col_fmt[i].last = FALSE;
  }
  dlg->col_fmt[0].pixInset = 5;
  dlg->col_fmt[3].last = TRUE;
  
  dlg->doc = DocumentPanel (p, width, stdLineHeight * 20);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocProcs (dlg->doc, ClickClickableItemList, NULL, NULL, NULL);
  SetDocShade (dlg->doc, DrawClickableItemList, NULL, NULL, NULL);

  return (DialoG) p;
}

typedef struct clickablelist
{
  DIALOG_MESSAGE_BLOCK
  DoC  doc;
  DialoG clickable_item_list;
  PrompT title1;
  PrompT title2;
  PrompT help1;
  PrompT help2;
  TexT find_txt;

  ClickableCallback item_single_click_callback;
  ClickableCallback item_double_click_callback;
  Pointer         item_click_callback_data;
  GetClickableItemText get_item_text;
  Boolean         dblClick;  
  Int2            clicked;
  Int2            selected;
  Int2            item_selected;  
  Int4            num_levels;
  ValNodePtr      list_list;
  Nlm_ColPtr PNTR col_fmt_array_array;

  Int2            text_select_item_start;
  Int2            text_select_row_start;
  Int2            text_select_char_start;
  Int2            text_select_item_stop;
  Int2            text_select_row_stop;
  Int2            text_select_char_stop;
  Int2            text_select_item_anchor;
  Int2            text_select_row_anchor;
  Int2            text_select_char_anchor;

  Boolean         display_chosen;
} ClickableListData, PNTR ClickableListPtr;


static Nlm_ParData clickableParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static Nlm_ColData clickableColFmt[2] = {{16, 0, 0, 0, NULL, 'l', 0,0,0,0, FALSE},
                                    {1000, 0, 0, 0, NULL, 'l', 1,0,0,0, TRUE}};


static Nlm_ColPtr PNTR FreeColumnFormatArrays (Nlm_ColPtr PNTR col_fmt_array_array, Int4 num_levels)
{
  Int4 n;
  
  if (col_fmt_array_array == NULL || num_levels < 1)
  {
    return NULL;
  }
  for (n = 0; n < num_levels; n++)
  {
    col_fmt_array_array [n] = MemFree (col_fmt_array_array [n]);
  }
  col_fmt_array_array = MemFree (col_fmt_array_array);
  return col_fmt_array_array;
}

static void CleanupClickableListDialog (GraphiC g, VoidPtr data)

{
  ClickableListPtr dlg;

  dlg = (ClickableListPtr) data;
  if (dlg != NULL) {
     dlg->col_fmt_array_array = FreeColumnFormatArrays (dlg->col_fmt_array_array, dlg->num_levels);
  }
  StdCleanupExtraProc (g, data);
}


static ClickableItemPtr GetSubItem (ValNodePtr item_list, Int2Ptr pitem)
{
  ClickableItemPtr cip = NULL;

  if (item_list == NULL || pitem == NULL)
  {
    return NULL;
  }
  while (*pitem > 0 && item_list != NULL)
  {
    (*pitem)--;
    cip = (ClickableItemPtr) item_list->data.ptrvalue;
    if (*pitem > 0)
    {
      if (cip != NULL && cip->expanded)
      {
        cip = GetSubItem (cip->subcategories, pitem);
      }
    }
    item_list = item_list->next;
  }
  if (*pitem > 0)
  {
    cip = NULL;
  }
  return cip;
}

static Int4 CountLevels (ValNodePtr item_list)
{
  Int4       num, num_sublevels = 0;
  ValNodePtr vnp;
  ClickableItemPtr cip;
  
  if (item_list == NULL) 
  {
    return 0;
  }
  
  for (vnp = item_list; vnp != NULL; vnp = vnp->next)
  {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip == NULL || cip->subcategories == NULL || !cip->expanded)
    {
      continue;
    }
    num = CountLevels (cip->subcategories);
    if (num > num_sublevels) num_sublevels = num;
  }
  
  /* one level for the top plus levels for the subcategories */

  return 1 + num_sublevels;
}


static ClickableItemPtr GetSelectedClickableList (ValNodePtr item_list, Int2 item)
{
  ClickableItemPtr cip = NULL;
  
  cip = GetSubItem (item_list, &item);

  return cip;
}


static Nlm_ColPtr PNTR GetColumnFormatArrays (Int4 num_levels, DoC doc)
{
  Int4               n, k;
  Nlm_ColPtr PNTR    col_fmt_array_array = NULL;
  RecT               r;
  Int4               doc_width;
    
  if (num_levels == 0)
  {
    return NULL;
  }
  
  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  doc_width = r.right - r.left;
  
  col_fmt_array_array = (Nlm_ColPtr PNTR) MemNew (sizeof (Nlm_ColPtr) * num_levels);
  for (n = 0; n < num_levels; n++)
  {
    col_fmt_array_array[n] = (Nlm_ColPtr) MemNew (sizeof (Nlm_ColData) * (n + 3));
    for (k = 0; k < n + 2; k++)
    {
      col_fmt_array_array[n][k].pixWidth = 16;
      col_fmt_array_array[n][k].pixInset = 0;
      col_fmt_array_array[n][k].charWidth = 0;
      col_fmt_array_array[n][k].charInset = 0;
      col_fmt_array_array[n][k].font = programFont;
      col_fmt_array_array[n][k].just = 'l';
      col_fmt_array_array[n][k].wrap = 0;
      col_fmt_array_array[n][k].bar = 0;
      col_fmt_array_array[n][k].underline = 0;
      col_fmt_array_array[n][k].left = 0;
      col_fmt_array_array[n][k].last = 0;
    }
    col_fmt_array_array[n][k].pixWidth = doc_width - ((n + 2) * 16);
    col_fmt_array_array[n][k].pixInset = 0;
    col_fmt_array_array[n][k].charWidth = 0;
    col_fmt_array_array[n][k].charInset = 0;
    col_fmt_array_array[n][k].font = programFont;
    col_fmt_array_array[n][k].just = 'l';
    col_fmt_array_array[n][k].wrap = 1;
    col_fmt_array_array[n][k].bar = 0;
    col_fmt_array_array[n][k].underline = 0;
    col_fmt_array_array[n][k].left = 0;
    col_fmt_array_array[n][k].last = 1;
  }
  return col_fmt_array_array;
}


static void AddClickableItem (ClickableListPtr dlg, ClickableItemPtr cip, Int4 level)
{
  CharPtr            item_text;
  ValNodePtr         vnp;
  Int4               n;

  if (cip == NULL)
  {
    return;
  }
  item_text = (CharPtr) MemNew (sizeof (Char) * (StringLen (cip->description) + 6 + level));
  for (n = 0; n < level; n++)
  {
    StringCat (item_text, "\t");
  }
  StringCat (item_text, " \t \t");
  StringCat (item_text, cip->description);
  StringCat (item_text, "\n");
  AppendText (dlg->doc, item_text, &clickableParFmt, dlg->col_fmt_array_array [level], programFont);
  if (cip->expanded)
  {
    for (vnp = cip->subcategories; vnp != NULL; vnp = vnp->next)
    {
      AddClickableItem (dlg, vnp->data.ptrvalue, level + 1);
    }
  }
}

static void PopulateClickableList (ClickableListPtr dlg, ValNodePtr list_list)
{
  Int2               numItems;
  Int4               num_levels;
  
  if (dlg == NULL || dlg->doc == NULL) 
  {
    return;
  }
  
  Reset (dlg->doc);
  
  num_levels = CountLevels (dlg->list_list);
  if (num_levels != dlg->num_levels)
  {
    dlg->col_fmt_array_array = FreeColumnFormatArrays (dlg->col_fmt_array_array, dlg->num_levels);
    dlg->num_levels = num_levels;
    dlg->col_fmt_array_array = GetColumnFormatArrays (dlg->num_levels, dlg->doc);
  }
  
  while (list_list != NULL)
  {
    AddClickableItem (dlg, list_list->data.ptrvalue, 0);
    list_list = list_list->next;
  }
  GetDocParams (dlg->doc, &numItems, NULL);
  UpdateDocument (dlg->doc, 0, numItems);

}


NLM_EXTERN Int2 PanelOffsetFromCharOffsetEx (DoC doc, FonT font, Int2 item, Int2 col, Int2 char_offset)
{
  Int2 numRows, numCols, lineHeight;
  Int2 left_start, width, inset, char_width;

  if (doc == NULL) return 0;
  GetItemParams4 (doc, item, NULL, &numRows, &numCols, &lineHeight, NULL);
  GetColParams (doc, item, col, &left_start, &width, &inset, NULL);

  /* need to set font so that Nlm_CharWidth works properly */
  SelectFont (font);
  char_width = Nlm_CharWidth ('A');
  return left_start + inset + char_offset * char_width;
}


static Int2 PanelOffsetFromCharOffset (ClickableListPtr dlg, Int2 item, Int2 char_offset)
{
  Int2 numCols;
  if (dlg == NULL) return 0;
  GetItemParams4 (dlg->doc, item, NULL, NULL, &numCols, NULL, NULL);

  return PanelOffsetFromCharOffsetEx (dlg->doc, dlg->col_fmt_array_array[0][0].font, item, numCols, char_offset);
}


NLM_EXTERN Int2 GetTextSelectCharOffsetEx (PoinT pt, DoC doc, FonT font, Int2 item, Int2 row, Int2 col)
{
  Int2 pixPos, pixInset;
  Int2 one_char_width;
  CharPtr txt;
  Int4    len;
  Int2    char_offset;

  if (doc == NULL) return 0;

  GetColParams (doc, item, col, &pixPos,
                NULL, &pixInset, NULL);

  SelectFont (font);
  one_char_width = Nlm_CharWidth ('A');
  if (one_char_width == 0) return 0;

  char_offset = (pt.x - pixPos - pixInset) / one_char_width;
  if (char_offset > 0) {
    txt = GetDocText (doc, item, row, col);
    len = StringLen (txt);
    if (char_offset >= len) {
      char_offset = - 1;
    }
  }

  return char_offset; 
}


static Int2 GetTextSelectCharOffset (PoinT pt, ClickableListPtr dlg, Int2 item, Int2 row, Int2 col)
{
  if (dlg == NULL) return 0;

  return GetTextSelectCharOffsetEx (pt, dlg->doc, dlg->col_fmt_array_array[0][0].font, item, row, col);
}


static void ValNodeLinkCopy (ValNodePtr PNTR list1, ValNodePtr list2)
{
  if (list1 == NULL) return;
  while (list2 != NULL)
  {
    ValNodeAddPointer (list1, list2->choice, list2->data.ptrvalue);
    list2 = list2->next;
  }
}


static ValNodePtr GetChosenItemsList (ValNodePtr clickable_list)
{
  ValNodePtr item_list = NULL;
  ClickableItemPtr cip;

  while (clickable_list != NULL) {
    cip = (ClickableItemPtr) clickable_list->data.ptrvalue;
    if (cip != NULL) {
      if (cip->chosen) {
        ValNodeLinkCopy (&item_list, cip->item_list);
      } else {
        ValNodeLink (&item_list, GetChosenItemsList (cip->subcategories));
      }
    }
    clickable_list = clickable_list->next;
  }
  return item_list;
}


static void ClickList (DoC d, PoinT pt)

{
  Int2             item, numItems;
  Int2             row;
  Int2             col;
  ClickableListPtr dlg;
  ClickableItemPtr cip;
  Int4             offset;
  Int2             first_shown;
  Boolean          rsult;
  RecT             r;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && dlg->clicked == item) {
      dlg->dblClick = dblClick;
    } else {
      dlg->dblClick = FALSE;
    }
    dlg->clicked = 0;
    if (item > 0 && row > 0) {
      dlg->clicked = item;
    }
    if (item > 0 && row > 0 && !dblClick)
    {
      cip = GetSelectedClickableList (dlg->list_list, item);
      if (cip != NULL)
      {
        if (col == cip->level + 1)
        {
          cip->chosen = !cip->chosen;
          GetDocParams (d, &numItems, NULL);
          UpdateDocument (d, 0, numItems);
          if (dlg->display_chosen) 
          {
            PointerToDialog (dlg->clickable_item_list, GetChosenItemsList (dlg->list_list));
          }
        }
        else if (col == cip->level + 2)
        {
          if (cip->subcategories != NULL) {
            cip->expanded = !cip->expanded;
            rsult = GetScrlParams4 (dlg->doc, &offset, &first_shown, NULL);
            PopulateClickableList (dlg, dlg->list_list);
            if (rsult) {
              GetItemParams4 (dlg->doc, first_shown, &offset, NULL, NULL, NULL, NULL);
              SetScrlParams4 (dlg->doc, offset);
            }
            ObjectRect (dlg->doc, &r);
            InsetRect (&r, -1, -1);
            InvalRect (&r);
          }
        } else {
          dlg->text_select_item_anchor = item;
          dlg->text_select_row_anchor = row;
          dlg->text_select_char_anchor = GetTextSelectCharOffset (pt, dlg, item, row, col);

          dlg->text_select_char_start = dlg->text_select_char_anchor;
          dlg->text_select_char_stop = dlg->text_select_char_start;
          if (dlg->text_select_char_start < 0) {
            dlg->text_select_item_start = -1;
            dlg->text_select_row_start = -1;
            dlg->text_select_item_stop = -1;
            dlg->text_select_row_stop = -1;
          } else {
            dlg->text_select_item_start = item;
            dlg->text_select_row_start = row;
            dlg->text_select_item_stop = item;
            dlg->text_select_row_stop = row;
          }
          GetDocParams (dlg->doc, &numItems, NULL);
          UpdateDocument (dlg->doc, 0, numItems);
        }
      }
    }
  }
}

static void UpdateClickableListTextSelection (ClickableListPtr dlg, Int2 item, Int2 row, Int2 col, PoinT pt)
{
  Int2 char_start, numCols;

  if (dlg == NULL || item < 1 || row < 1) return;

  /* only update for positions in the text area */
  GetItemParams4 (dlg->doc, item, NULL, NULL, &numCols, NULL, NULL);
  if (col != numCols) {
    return;
  }

  char_start = GetTextSelectCharOffset (pt, dlg, item, row, col);
  if (char_start < 0) {
    /* no updates unless mouse is in text area */
    return;
  }
    
  if (item < dlg->text_select_item_anchor) {
    /* mouse is before anchor */
    dlg->text_select_item_start = item;
    dlg->text_select_row_start = row;
    dlg->text_select_char_start = char_start;
    dlg->text_select_item_stop = dlg->text_select_item_anchor;
    dlg->text_select_row_stop = dlg->text_select_row_anchor;
    dlg->text_select_char_stop = dlg->text_select_char_anchor;
  } else if (item == dlg->text_select_item_anchor) {
    /* start and stop in the anchor item */
    dlg->text_select_item_start = item;
    dlg->text_select_item_stop = item;
    if (row < dlg->text_select_row_anchor) {
      /* mouse is before anchor */
      dlg->text_select_row_start = row;
      dlg->text_select_char_start = char_start;
      dlg->text_select_row_stop = dlg->text_select_row_anchor;
      dlg->text_select_char_stop = dlg->text_select_char_anchor;
    } else if (row == dlg->text_select_row_anchor) {
      /* start and stop in the anchor row */
      dlg->text_select_row_start = row;
      dlg->text_select_row_stop = row;
      if (char_start <= dlg->text_select_char_anchor) {
        /* mouse is before anchor */
        dlg->text_select_char_start = char_start;
        dlg->text_select_char_stop = dlg->text_select_char_anchor;
      } else {
        dlg->text_select_char_start = dlg->text_select_char_anchor;
        dlg->text_select_char_stop = char_start;
      }
    }
  } else {
    /* mouse is after anchor */
    dlg->text_select_item_start = dlg->text_select_item_anchor;
    dlg->text_select_row_start = dlg->text_select_row_anchor;
    dlg->text_select_char_start = dlg->text_select_char_anchor;
    dlg->text_select_item_stop = item;
    dlg->text_select_row_stop = row;
    dlg->text_select_char_stop = char_start;
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
}

static void DragClickableList (DoC d, PoinT pt)
{
  Int2            item, row, col, numCols;

  ClickableListPtr dlg;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    GetItemParams4 (d, item, NULL, NULL, &numCols, NULL, NULL);
    if (col == numCols) {
      UpdateClickableListTextSelection (dlg, item, row, col, pt);
    }
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

static void ActOnClickableList (ValNodePtr list_list, Int2 item)
{
  ClickableItemPtr cip;
  
  cip = GetSelectedClickableList (list_list, item);
  if (cip != NULL && cip->callback_func != NULL)
  {
    (cip->callback_func) (cip->item_list, cip->callback_data);
  }
}


static void PopulateClickableItemList (DialoG d, ClickableItemPtr cip)
{
  ValNodePtr        list = NULL;
  
  if (d == NULL) return;
  if (cip != NULL)
  {
    ValNodeLinkCopy (&list, cip->item_list);
  }
  PointerToDialog (d, list);

}


static void ReleaseClickableList (DoC d, PoinT pt)

{
  Int2           item;
  Int2           old;
  Int2           row, col;
  ClickableListPtr dlg;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    ResetClip ();
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    /* update text selection */
    UpdateClickableListTextSelection (dlg, item, row, col, pt);

    if (item > 0 && row > 0) {

      if (item == dlg->clicked) {
        old = dlg->selected;
        dlg->selected = item;
        if (old != item) {
          if (old == 0) {
            UpdateDocument (d, item, item);
          } else {
            UpdateDocument (d, old, old);
            UpdateDocument (d, item, item);
          }
          Update ();
        }
      }
    } else if (dlg->clicked == 0) {
      if (dlg->selected != 0) {
        old = dlg->selected;
        dlg->selected = 0;
        InvalBorder (d, old);
      }
      Update ();
    }
    if (dlg->selected > 0 && dlg->dblClick)
    {
      ActOnClickableList (dlg->list_list, dlg->selected);
    }
    else if (dlg->selected > 0)
    {
      dlg->item_selected = 0;
      if (!dlg->display_chosen) {
        PopulateClickableItemList (dlg->clickable_item_list, 
                                  GetSelectedClickableList (dlg->list_list,
                                                            dlg->selected));
      }
    }
  }
}


static void DrawTextSelection (ClickableListPtr dlg, Int2 item, RectPtr r)
{
  Int2 lineHeight, numRows, numCols;
  Int2 last_right;
  Int2 left_start;
  Int4 top, line_y;
  CharPtr txt;

  if (dlg == NULL || r == NULL
      || item < dlg->text_select_item_start
      || item > dlg->text_select_item_stop) {
    return;
  }

  if (dlg->text_select_item_start == dlg->text_select_item_stop
      && dlg->text_select_row_start == dlg->text_select_row_stop
      && dlg->text_select_char_start == dlg->text_select_char_stop) {
    /* if we've only selected one char, and it's blank, don't draw it. */
    txt = GetSelectedClickableListText (dlg->dialog);
    if (StringHasNoText (txt)) {
      MemFree (txt);
      return;
    }
    MemFree (txt);
  }

  GetItemParams4 (dlg->doc, item, &top, &numRows, &numCols, &lineHeight, NULL);

  /* calculate missing rows from end first */
  if (dlg->text_select_item_stop == item) {
    numRows = dlg->text_select_row_stop;
    last_right = PanelOffsetFromCharOffset (dlg, item, dlg->text_select_char_stop + 1);
  } else {
    last_right = r->right;
  }

  if (dlg->text_select_item_start == item) {
    left_start = PanelOffsetFromCharOffset (dlg, item, dlg->text_select_char_start);
    line_y = r->top + (dlg->text_select_row_start) * lineHeight - 1;
    numRows -= dlg->text_select_row_start - 1;
  } else {
    left_start = PanelOffsetFromCharOffset (dlg, item, 0);
    line_y = r->top + lineHeight - 1;
  }

  while (numRows > 1) {
    MoveTo (left_start, line_y);
    left_start = PanelOffsetFromCharOffset (dlg, item, 0);
    LineTo (r->right, line_y);
    line_y += lineHeight;
    numRows--;
  }
  MoveTo (left_start, line_y);
  LineTo (last_right, line_y);

}

static void DrawClickableList (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  ClickableListPtr dlg;
  RecT             rct;
  ClickableItemPtr cip;
  Int4             level_offset;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    cip = GetSelectedClickableList (dlg->list_list, item);
    if (cip != NULL)
    {
      level_offset = cip->level * 16;
      rct.left += level_offset;
      rct.right += level_offset;
    }

    /* draw text selection */
    DrawTextSelection (dlg, item, &rct);

    /* draw selection */
    if (item == dlg->selected) {
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }

    /* draw chosen checkboxes */
    rct.left += 5;
    rct.right = rct.left + 10;
    rct.bottom = rct.top + (rct.right - rct.left);
    FrameRect (&rct);
    
    if (cip != NULL && cip->chosen) {
      MoveTo (rct.left, rct.top);
      LineTo (rct.right - 1, rct.bottom - 1);
      MoveTo (rct.left, rct.bottom - 1);
      LineTo (rct.right - 1, rct.top);
    }
    
    /* draw open/closed checkboxes */
    if (cip!= NULL && cip->subcategories != NULL)
    {
      rct.left += 10;
      rct.right = rct.left + 10;
      rct.bottom = rct.top + (rct.right - rct.left);
      FrameRect (&rct);
      MoveTo (rct.left, (rct.top + rct.bottom) / 2);
      LineTo (rct.right - 1, (rct.top + rct.bottom) / 2);
      if (!cip->expanded)
      {
        MoveTo ((rct.left + rct.right) / 2, rct.top);
        LineTo ((rct.left + rct.right) / 2, rct.bottom - 1);
      }
    }

  }
}



static void DrawClickableListItem (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  ClickableListPtr dlg;
  RecT             rct;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    /* draw selection */
    if (item == dlg->item_selected) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
  }
}


static void ClickClickableListItem (DoC d, PoinT pt)

{
  Int2             item, last_selected, numItems;
  Int2             row;
  ClickableListPtr dlg;
  ClickableItemPtr cip;
  ValNodePtr       vnp;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {  
      cip = GetSelectedClickableList (dlg->list_list, dlg->selected);
      if (cip != NULL && cip->item_list != NULL)
      {
        vnp = cip->item_list;
        
        last_selected = dlg->item_selected;
        dlg->item_selected = item;
        
        if (item != last_selected)
        {
          GetDocParams (d, &numItems, NULL);
          UpdateDocument (d, 0, numItems);
        }
    
        /* find item in list */
        while (item > 1 && vnp != NULL)
        {
          vnp = vnp->next;
          item--;
        }
        
        if (dblClick)
        {
          if (dlg->item_double_click_callback != NULL) {
            (dlg->item_double_click_callback) (vnp, dlg->item_click_callback_data);
          }
        } else {
          if (dlg->item_single_click_callback != NULL) {
            (dlg->item_single_click_callback) (vnp, dlg->item_click_callback_data);
          }
        }          
      }
    }
  }
}


static void ClickableListToDialog (DialoG d, Pointer userdata)
{
  ClickableListPtr dlg;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->list_list = (ValNodePtr) userdata;

  PopulateClickableList (dlg, dlg->list_list);
  if (dlg->list_list != NULL) {
    dlg->selected = 1;
    dlg->item_selected = 0;
    if (dlg->display_chosen)  {
      PointerToDialog (dlg->clickable_item_list, GetChosenItemsList (dlg->list_list));
    } else {
      PopulateClickableItemList (dlg->clickable_item_list, 
                                GetSelectedClickableList (dlg->list_list,
                                                          dlg->selected));
    }
  } else {
    PointerToDialog (dlg->clickable_item_list, NULL);
  }
}

static void FindClickableText (ButtoN b) 
{
  ClickableListPtr dlg;
  CharPtr          find_txt;

  dlg = (ClickableListPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  find_txt = SaveStringFromText (dlg->find_txt);
  ScrollToNextClickableTextDescription (find_txt, dlg->dialog);
  find_txt = MemFree (find_txt);
}


static void FindPreviousClickableText (ButtoN b)
{
  ClickableListPtr dlg;
  CharPtr          find_txt;

  dlg = (ClickableListPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  find_txt = SaveStringFromText (dlg->find_txt);
  ScrollToPreviousClickableTextDescription (find_txt, dlg->dialog);
  find_txt = MemFree (find_txt);
}


static void ClickableListCopyToClipboard (DialoG d)
{
  ClickableListPtr dlg;
  CharPtr          txt;
  
  dlg = (ClickableListPtr) GetObjectExtra (d);

  if (dlg == NULL) return;

  txt = GetSelectedClickableListText (d);
  Nlm_StringToClipboard (txt);
  txt = MemFree (txt);
}

static void ClickableListMessage (DialoG d, Int2 mssg)

{
  ClickableListPtr dlg;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == VIB_MSG_COPY) {
      ClickableListCopyToClipboard(d);
    }
  }
}

static void ClickableListOnKey (SlatE s, Char ch)
{
  DoC              d;
  ClickableListPtr clp;
#ifdef WIN_MSWIN
  Int2             first_shown;
  Int4             offset;
  BaR              sb;
#endif

  if ( (int) ch == 0 ) return;

  d = (DoC) s;
  clp = (ClickableListPtr) GetObjectExtra (d);

  CaptureSlateFocus (s);
  /* later, handle control key combos */
#ifdef WIN_MSWIN
  if (ch == 3)
  {
    ClickableListCopyToClipboard (clp->dialog);
  }
  else if (ch == 11) 
  {
    /* PageUp key pressed */
    if (GetScrlParams4 (d, &offset, &first_shown, NULL) && first_shown > 0) 
    {
      sb = GetSlateVScrollBar (s);
      Nlm_Scroll (sb, SCROLL_PAGEUP);
    }
  }
  else if (ch == 12)
  {
    /* PageDown key pressed */
    if (GetScrlParams4 (d, &offset, &first_shown, NULL)) 
    {
      sb = GetSlateVScrollBar (s);
      Nlm_Scroll (sb, SCROLL_PAGEDN);
    }
  }
  else if (ch == 30)
  {
    /* Up Arrow key pressed */
    if (GetScrlParams4 (d, &offset, &first_shown, NULL) && first_shown > 0) 
    {
      GetItemParams4 (d, first_shown - 1, &offset, NULL, NULL, NULL, NULL);
      SetScrlParams4 (d, offset);
    }
  }
  else if (ch == 31)
  {
    /* Down Arrow key pressed */
    if (GetScrlParams4 (d, &offset, &first_shown, NULL)) {
      sb = GetSlateVScrollBar (s);
      if (offset < GetBarMax (sb)) 
      {
        GetItemParams4 (d, first_shown + 1, &offset, NULL, NULL, NULL, NULL);
        SetScrlParams4 (d, offset);
      }
    }
  }
#endif
}


extern DialoG 
CreateClickableListDialogExEx 
(GrouP h, 
 CharPtr label1, 
 CharPtr label2,
 CharPtr help1,
 CharPtr help2,
 ClickableCallback item_single_click_callback,
 ClickableCallback item_double_click_callback,
 Pointer         item_click_callback_data,
 GetClickableItemText get_item_text,
 Int4            left_width,
 Int4            right_width,
 Boolean         horizontal,
 Boolean         show_find,
 Boolean         display_chosen)
{
  GrouP p, pnl_grp, find_grp = NULL;
  ClickableListPtr dlg;
  RecT             r;
  ButtoN           b;

  dlg = (ClickableListPtr) MemNew (sizeof (ClickableListData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupClickableListDialog);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = ClickableListToDialog;
  dlg->fromdialog = NULL;
  dlg->dialogmessage = ClickableListMessage;
  dlg->testdialog = NULL;
  
  dlg->item_single_click_callback = item_single_click_callback;
  dlg->item_double_click_callback = item_double_click_callback;
  dlg->item_click_callback_data = item_click_callback_data;
  
  dlg->get_item_text = get_item_text;

  dlg->display_chosen = display_chosen;
  
  if (horizontal) {
    pnl_grp = HiddenGroup (p, 2, 0, NULL);
    
    if (label1 || label2) {
      dlg->title1 = StaticPrompt (pnl_grp, label1, left_width, popupMenuHeight, programFont, 'c');
      dlg->title2 = StaticPrompt (pnl_grp, label2, right_width, popupMenuHeight, programFont, 'c');
    }
    dlg->doc = DocumentPanel (pnl_grp, left_width, stdLineHeight * 20);
    dlg->clickable_item_list = ClickableItemListDialog (pnl_grp, right_width, 
                                                        get_item_text,
                                                        item_single_click_callback,
                                                        item_double_click_callback, 
                                                        item_click_callback_data);
    if (help1 || help2) {
      dlg->help1 = StaticPrompt (pnl_grp, help1, left_width, popupMenuHeight, programFont, 'c');
      dlg->help2 = StaticPrompt (pnl_grp, help2, right_width, popupMenuHeight, programFont, 'c');
    } 
  } else {
    pnl_grp = HiddenGroup (p, -1, 0, NULL);
    dlg->title1 = StaticPrompt (pnl_grp, label1, left_width, popupMenuHeight, programFont, 'c');
    dlg->doc = DocumentPanel (pnl_grp, left_width, stdLineHeight * 20);
    if (help1 || help2) {
      dlg->help1 = StaticPrompt (pnl_grp, help1, left_width, popupMenuHeight, programFont, 'c');
    }
    dlg->title2 = StaticPrompt (pnl_grp, label2, right_width, popupMenuHeight, programFont, 'c');
    dlg->clickable_item_list = ClickableItemListDialog (pnl_grp, right_width, 
                                                        get_item_text,
                                                        item_single_click_callback,
                                                        item_double_click_callback, 
                                                        item_click_callback_data);
    if (help1 || help2) {
      dlg->help2 = StaticPrompt (pnl_grp, help2, right_width, popupMenuHeight, programFont, 'c');
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->title1, (HANDLE) dlg->doc, (HANDLE) dlg->title2, (HANDLE) dlg->clickable_item_list, 
                                (HANDLE) (dlg->help1 == NULL ? dlg->help2 : dlg->help1), (HANDLE) (dlg->help1 == NULL ? NULL : dlg->help2), NULL);
  }

  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocProcs (dlg->doc, ClickList, DragClickableList, ReleaseClickableList, NULL);
  SetDocShade (dlg->doc, DrawClickableList, NULL, NULL, NULL);
  
  SetSlateChar ((SlatE) dlg->doc, ClickableListOnKey);
  
  /* adjust column width for discrepancy list */
  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  clickableColFmt[1].pixWidth = r.right - r.left - clickableColFmt[0].pixWidth;
  
  if (show_find) {
    find_grp = HiddenGroup (p, 4, 0, NULL);
    SetGroupSpacing (find_grp, 10, 10);
    StaticPrompt (find_grp, "Find Text", 0, popupMenuHeight, programFont, 'l');
    dlg->find_txt = DialogText (find_grp, "", 20, NULL);
    b = PushButton (find_grp, "<<", FindPreviousClickableText);
    SetObjectExtra (b, dlg, NULL);
    b = PushButton (find_grp, ">>", FindClickableText);
    SetObjectExtra (b, dlg, NULL);
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) pnl_grp, (HANDLE) find_grp, NULL);
  return (DialoG) p;
}


extern DialoG 
CreateClickableListDialogEx 
(GrouP h, 
 CharPtr label1, 
 CharPtr label2,
 CharPtr help1,
 CharPtr help2,
 ClickableCallback item_single_click_callback,
 ClickableCallback item_double_click_callback,
 Pointer         item_click_callback_data,
 GetClickableItemText get_item_text,
 Int4            left_width,
 Int4            right_width,
 Boolean         horizontal,
 Boolean         show_find)
{
  return CreateClickableListDialogExEx (h, label1, label2, help1, help2, 
                                        item_single_click_callback, item_double_click_callback,
                                        item_click_callback_data,
                                        get_item_text,
                                        left_width,
                                        right_width,
                                        horizontal,
                                        show_find,
                                        FALSE);
}


extern void SetClickableListDialogTitles (DialoG d, CharPtr title1, CharPtr title2, CharPtr help1, CharPtr help2)
{
  ClickableListPtr dlg;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  SafeSetTitle (dlg->title1, title1);
  SafeSetTitle (dlg->title2, title2);
  SafeSetTitle (dlg->help1, help1);
  SafeSetTitle (dlg->help2, help2);
}


extern void SetClickableListDisplayChosen (DialoG d, Boolean set)
{
  ClickableListPtr dlg;

  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  dlg->display_chosen = set;
  if (dlg->display_chosen)  {
    PointerToDialog (dlg->clickable_item_list, GetChosenItemsList (dlg->list_list));
  } else {
    PopulateClickableItemList (dlg->clickable_item_list, 
                              GetSelectedClickableList (dlg->list_list,
                                                        dlg->selected));
  }
}  


extern DialoG 
CreateClickableListDialog 
(GrouP h, 
 CharPtr label1, 
 CharPtr label2,
 ClickableCallback item_single_click_callback,
 ClickableCallback item_double_click_callback,
 Pointer         item_click_callback_data,
 GetClickableItemText get_item_text)
{
  return CreateClickableListDialogEx (h, label1, label2, NULL, NULL,
                                      item_single_click_callback,
                                      item_double_click_callback,
                                      item_click_callback_data,
                                      get_item_text, 
                                      stdCharWidth * 30,
                                      stdCharWidth * 30 + 5,
                                      TRUE, TRUE);
}


static Int4 CountExpandedSublevels (ValNodePtr subcategories)
{
  ValNodePtr vnp;
  Int4       num = 0;
  ClickableItemPtr cip;

  for (vnp = subcategories; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      num++;
      if (cip->expanded && cip->subcategories != NULL) {
        num += CountExpandedSublevels(cip->subcategories);
      }
    }
  }
  return num;
}


/* First, need to locate items that have description that contains txt.
 * Then need to make sure that they are visible.
 * Then need to figure out how to scroll to them.
 */
static Boolean FindInClickableItemDescriptions (CharPtr txt, ValNodePtr clickable_item_list, ValNodePtr PNTR row_list, Int4Ptr row_offset)
{
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  Boolean          found = FALSE;
  Int4             lower_levels;
  
  if (StringHasNoText (txt) || clickable_item_list == NULL || row_list == NULL || row_offset == NULL) return FALSE;
  
  vnp = clickable_item_list;
  while (vnp != NULL) {
    (*row_offset)++;
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (StringStr (cip->description, txt) != NULL) {
      ValNodeAddInt (row_list, 0, *row_offset);
      found = TRUE;
    }
    lower_levels = *row_offset;
    if (cip->subcategories != NULL) {
      if (FindInClickableItemDescriptions (txt, cip->subcategories, row_list, &lower_levels)) {
        cip->expanded = TRUE;
        found = TRUE;
      }
      if (cip->expanded) {
        (*row_offset) += CountExpandedSublevels (cip->subcategories);
      }
    }
    vnp = vnp->next;
  }
  
  return found;
}


extern void ScrollToNextClickableTextDescription (CharPtr txt, DialoG d)
{
  ClickableListPtr dlg;
  ValNodePtr       found_list = NULL, vnp;
  Int4             row_offset = 0, startsAt;
  Int2             numRows, numCols, lineHeight;
  
  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }

  if (!FindInClickableItemDescriptions (txt, dlg->list_list, &found_list, &row_offset)
      || found_list == NULL) {
    Message (MSG_OK, "Text not found!");
    return;
  }
 
  /* find first found txt after current selection */
  vnp = found_list;
  while (vnp != NULL && vnp->data.intvalue <= dlg->selected) {
    vnp = vnp->next;
  }
  if (vnp == NULL) {
   row_offset = found_list->data.intvalue;
  } else {
    row_offset = vnp->data.intvalue;
  }

  dlg->selected = row_offset; 
  PopulateClickableList (dlg, dlg->list_list);
  SetDocAutoAdjust (dlg->doc, TRUE);   
  UpdateDocument (dlg->doc, 0, 0);
  SetDocAutoAdjust (dlg->doc, FALSE);
  dlg->item_selected = 0;
  PopulateClickableItemList (dlg->clickable_item_list, 
                             GetSelectedClickableList (dlg->list_list,
                                                       dlg->selected));
 
  GetItemParams4 (dlg->doc, row_offset, &startsAt, &numRows,
                  &numCols, &lineHeight, NULL);
  SetScrlParams4 (dlg->doc, startsAt);
}


extern void ScrollToPreviousClickableTextDescription (CharPtr txt, DialoG d)
{
  ClickableListPtr dlg;
  ValNodePtr       found_list = NULL, vnp, vnp_prev = NULL;
  Int4             row_offset = 0, startsAt;
  Int2             numRows, numCols, lineHeight;
  
  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }

  if (!FindInClickableItemDescriptions (txt, dlg->list_list, &found_list, &row_offset)
      || found_list == NULL) {
    Message (MSG_OK, "Text not found!");
    return;
  }
 
  /* find first found txt before current selection */
  vnp = found_list;
  while (vnp != NULL && vnp->data.intvalue < dlg->selected) {
    vnp_prev = vnp;
    vnp = vnp->next;
  }
  if (vnp_prev == NULL) {
    /* use last item in list */
    vnp = found_list;
    while (vnp != NULL && vnp->next != NULL) {
      vnp = vnp->next;
    }
    row_offset = vnp->data.intvalue;
  } else {
    row_offset = vnp_prev->data.intvalue;
  }

  dlg->selected = row_offset; 
  PopulateClickableList (dlg, dlg->list_list);
  SetDocAutoAdjust (dlg->doc, TRUE);   
  UpdateDocument (dlg->doc, 0, 0);
  SetDocAutoAdjust (dlg->doc, FALSE);
  dlg->item_selected = 0;
  PopulateClickableItemList (dlg->clickable_item_list, 
                             GetSelectedClickableList (dlg->list_list,
                                                       dlg->selected));
 
  GetItemParams4 (dlg->doc, row_offset, &startsAt, &numRows,
                  &numCols, &lineHeight, NULL);
  SetScrlParams4 (dlg->doc, startsAt);
}


static CharPtr GetClickableTextFromItemList (ValNodePtr item_list, CharPtr separator) 
{
  Int4 text_len = 1;
  ValNodePtr vnp;
  CharPtr    txt;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    text_len += StringLen (vnp->data.ptrvalue) + StringLen (separator);
  }

  txt = (CharPtr) MemNew (sizeof(Char) * text_len);
  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    StringCat (txt, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (txt, separator);
    }
  }
  return txt;
}


static CharPtr
GetFragmentFromDocRow
(DoC  doc,
 Int2 item,
 Int2 col,
 Int2 row,
 Int2 char_start, Int2 char_stop)
{
  CharPtr txt = NULL, cp_txt;
  Int4    len;
 
  txt = GetDocText (doc, item, row, col);

  if (char_start == 0 && char_stop < 0) {
    /* take all of row */
    return txt;
  } else {
    len = StringLen (txt);
    if (char_stop >= len) {
      /* selection was drawn beyond length of text */
      char_stop = len - 1;
    } else if (char_stop < 0) {
      /* take all of row */
      char_stop = len - 1;
    }
    if (char_start >= len) {
      /* do nothing - selection is not within text */
      txt = MemFree (txt);
    } else {
      cp_txt = (CharPtr) MemNew (sizeof (Char) * (2 + char_stop - char_start));
      StringNCpy (cp_txt, txt + char_start, 1 + char_stop - char_start);
      cp_txt[1 + char_stop - char_start] = 0;
      txt = MemFree (txt);
      txt = cp_txt;
    }
  }  
  return txt;
}


static CharPtr 
GetFragmentsFromDocCol
(DoC  doc,
 Int2 item,
 Int2 col,
 Int2 row_start, Int2 row_stop,
 Int2 char_start, Int2 char_stop)
{
  Int2       row;
  CharPtr    txt;
  ValNodePtr fragments = NULL;

  if (row_start == row_stop) {
    txt = GetFragmentFromDocRow (doc, item, col, row_start, char_start, char_stop);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    } else {
      ValNodeAddPointer (&fragments, 0, txt);
    }
  } else {
    txt = GetFragmentFromDocRow (doc, item, col, row_start, char_start, -1);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    } else {
      ValNodeAddPointer (&fragments, 0, txt);
    }
    row = row_start + 1;
    while (row < row_stop) {
      txt = GetFragmentFromDocRow (doc, item, col, row, 0, -1);
      if (StringHasNoText (txt)) {
        txt = MemFree (txt);
      } else {
        ValNodeAddPointer (&fragments, 0, txt);
      }
      row++;
    }
    txt = GetFragmentFromDocRow (doc, item, col, row_stop, 0, char_stop);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    } else {
      ValNodeAddPointer (&fragments, 0, txt);
    }
  }
  txt = GetClickableTextFromItemList (fragments, " ");
  fragments = ValNodeFreeData (fragments);
  return txt;
}


static Boolean CollectFromThisColumn (Int2 col, Int2Ptr only_these_columns, Int2 num_col)
{
  Int2 i;
  Boolean rval = FALSE;

  if (only_these_columns == NULL) {
    rval = TRUE;
  } else {
    for (i = 0; i < num_col; i++) {
      if (only_these_columns[i] == col) {
        rval = TRUE;
      }
    }
  }
  return rval;
}

 
static CharPtr
GetFragmentsFromDocItem 
(DoC  doc,
 Int2 item, 
 Int2 col_start, Int2 col_stop,
 Int2 row_start, Int2 row_stop,
 Int2 char_start, Int2 char_stop,
 Int2Ptr only_these_columns, Int2 num_col)
{
  Int2       col;
  CharPtr    txt;
  Int2       num_rows, num_cols;
  ValNodePtr fragments = NULL;

  GetItemParams4 (doc, item, NULL, &num_rows, &num_cols, NULL, NULL);
  if (col_stop < 0) {
    col_stop = num_cols;
  }
  if (row_stop < 0) {
    row_stop = num_rows;
  }

  if (col_start == col_stop) {
    if (CollectFromThisColumn (col_start, only_these_columns, num_col)) {
      txt = GetFragmentsFromDocCol (doc, item, col_start, row_start, row_stop, char_start, char_stop);
      ValNodeAddPointer (&fragments, 0, txt);
    }
  } else {
    if (CollectFromThisColumn (col_start, only_these_columns, num_col)) {
      txt = GetFragmentsFromDocCol (doc, item, col_start, row_start, num_rows, char_start, -1);
      ValNodeAddPointer (&fragments, 0, txt);
    }
    col = col_start + 1;
    while (col < col_stop) {
      if (CollectFromThisColumn (col, only_these_columns, num_col)) {
        txt = GetFragmentsFromDocCol (doc, item, col, 1, num_rows, 0, -1);
        ValNodeAddPointer (&fragments, 0, txt);
      }
      col++;
    }
    if (CollectFromThisColumn (col_stop, only_these_columns, num_col)) {
      txt = GetFragmentsFromDocCol (doc, item, col_stop, 1, row_stop, 0, char_stop);
      ValNodeAddPointer (&fragments, 0, txt);
    }
  }
  txt = GetClickableTextFromItemList (fragments, "\t");
  fragments = ValNodeFreeData (fragments);

  return txt;
}


NLM_EXTERN CharPtr GetSelectedDocText (DoC doc, Int2 item_start, Int2 row_start, Int2 col_start, Int2 char_start,
                                       Int2 item_stop, Int2 row_stop, Int2 col_stop, Int2 char_stop,
                                       Int2Ptr only_these_columns, Int2 num_col)
{
  Int2             item;
  CharPtr          txt;
  ValNodePtr       fragment_list = NULL;
  
  if (doc == NULL || item_start < 1)
  {
    return NULL;
  }

  if (char_start < 0) char_start = 0;

  if (item_start == item_stop) {
    txt = GetFragmentsFromDocItem (doc, item_start, col_start, col_stop, row_start, row_stop, char_start, char_stop, only_these_columns, num_col); 
    ValNodeAddPointer (&fragment_list, 0, txt);
  } else {
    txt = GetFragmentsFromDocItem (doc, item_start, col_start, -1, row_start, -1, char_start, -1, only_these_columns, num_col);
    ValNodeAddPointer (&fragment_list, 0, txt);
    item = item_start + 1;
    while (item < item_stop) {
      txt = GetFragmentsFromDocItem (doc, item, 1, -1, 1, -1, 0, -1, only_these_columns, num_col);
      ValNodeAddPointer (&fragment_list, 0, txt);
      item++;
    }
    txt = GetFragmentsFromDocItem (doc, item, 1, col_stop, 1, row_stop, 0, char_stop, only_these_columns, num_col);
    ValNodeAddPointer (&fragment_list, 0, txt);
  }

  txt = GetClickableTextFromItemList (fragment_list, "\r\n");
  fragment_list = ValNodeFreeData (fragment_list);
  return txt;
}

extern CharPtr GetSelectedClickableListText (DialoG d)
{
  ClickableListPtr dlg;
  Int2             numRows;
  Int2             item, row, col, char_offset;
  Int4             text_len = 0, len = 0;
  CharPtr          txt, cp_txt;
  ValNodePtr       fragment_list = NULL, item_list = NULL;
  
  dlg = (ClickableListPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->text_select_item_start < 1)
  {
    return NULL;
  }

  item = dlg->text_select_item_start;
  row = dlg->text_select_row_start;
  char_offset = dlg->text_select_char_start;

  while (item <= dlg->text_select_item_stop) {
    if (item == dlg->text_select_item_start) {
      if (item == dlg->text_select_item_stop) {
        /* all text from one item */
        /* all text will be from last column in item */
        GetItemParams4 (dlg->doc, item, NULL, NULL, &col, NULL, NULL);
        if (row == dlg->text_select_row_stop) {
          /* all text from one row */
          txt = GetDocText (dlg->doc, item, row, col);
          len = StringLen (txt);
          if (dlg->text_select_char_stop >= len) {
            /* selection was drawn beyond length of text */
            dlg->text_select_char_stop = len - 1;
          }
          if (dlg->text_select_char_start >= len) {
            /* do nothing - selection is not within text */
          } else {
            cp_txt = (CharPtr) MemNew (sizeof (Char) * (2 + dlg->text_select_char_stop - dlg->text_select_char_start));
            StringNCpy (cp_txt, txt + dlg->text_select_char_start, 1 + dlg->text_select_char_stop - dlg->text_select_char_start);
            cp_txt[1 + dlg->text_select_char_stop - dlg->text_select_char_start] = 0;
            ValNodeAddPointer (&fragment_list, 0, cp_txt);
          }
          /* free txt read from Doc */
          txt = MemFree (txt);
        } else {
          /* take text from several rows */
          item_list = NULL;
          while (row < dlg->text_select_row_stop) {
            txt = GetDocText (dlg->doc, item, row, col);
            len = StringLen (txt);
            if (row == dlg->text_select_row_start && dlg->text_select_char_start >= len) {
              /* do nothing - selection is not within text */
              txt = MemFree (txt);
            } else {
              if (row == dlg->text_select_row_start && dlg->text_select_char_start > 0) {
                cp_txt = (CharPtr) MemNew (sizeof (Char) * (1 + StringLen (txt) - dlg->text_select_char_start));
                StringCpy (cp_txt, txt + dlg->text_select_char_start);
                txt = MemFree (txt);
                txt = cp_txt;
              }
              ValNodeAddPointer (&item_list, 0, txt);
            }
            row++;
          }
          txt = GetDocText (dlg->doc, item, row, col);
          if (dlg->text_select_char_stop < len) {
            txt[dlg->text_select_char_stop] = 0;
          }
          ValNodeAddPointer (&item_list, 0, txt);
          txt = GetClickableTextFromItemList (item_list, " ");
          ValNodeAddPointer (&fragment_list, 0, txt);
          item_list = ValNodeFreeData (item_list);
        }
      } else {
        /* take all text after first row and char offset */
        item_list = NULL;
        GetItemParams4 (dlg->doc, item, NULL, &numRows, &col, NULL, NULL);
        while (row <= numRows) {
          txt = GetDocText (dlg->doc, item, row, col);
          len = StringLen (txt);
          if (row == dlg->text_select_row_start && dlg->text_select_char_start >= len) {
            /* do nothing = selection is outside text */
            txt = MemFree (txt);
          } else {
            if (row == dlg->text_select_row_start && dlg->text_select_char_start > 0) {
              cp_txt = (CharPtr) MemNew (sizeof (Char) * (1 + StringLen (txt) - dlg->text_select_char_start));
              StringCpy (cp_txt, txt + dlg->text_select_char_start);
              txt = MemFree (txt);
              txt = cp_txt;
            }
            ValNodeAddPointer (&item_list, 0, txt);
          }
          row++;
        }
        txt = GetClickableTextFromItemList (item_list, " ");
        ValNodeAddPointer (&fragment_list, 0, txt);
        item_list = ValNodeFreeData (item_list);
      }
    } else if (item == dlg->text_select_item_stop) {
      /* take all text until last row and char offset */
      row = 1;
      item_list = NULL;
      GetItemParams4 (dlg->doc, item, NULL, NULL, &col, NULL, NULL);
      while (row < dlg->text_select_row_stop) {
        txt = GetDocText (dlg->doc, item, row, col);
        ValNodeAddPointer (&item_list, 0, txt);
      }
      txt = GetDocText (dlg->doc, item, row, col);
      if (dlg->text_select_char_stop + 1 < StringLen (txt)) {
        txt[dlg->text_select_char_stop + 1] = 0;
      }
      ValNodeAddPointer (&item_list, 0, txt);
      txt = GetClickableTextFromItemList (item_list, " ");
      ValNodeAddPointer (&fragment_list, 0, txt);
      item_list = ValNodeFreeData (item_list);
    } else {
      GetItemParams4 (dlg->doc, item, NULL, NULL, &col, NULL, NULL);
      txt = GetDocText (dlg->doc, item, 0, col);
      if (txt != NULL && txt[StringLen(txt) - 1] == '\n') {
        /* remove terminal carriage return */
        txt[StringLen(txt) - 1] = 0;
      }
      ValNodeAddPointer (&fragment_list, 0, txt);
      text_len += StringLen (txt) + 1;
    }
    item++;
  }
  txt = GetClickableTextFromItemList (fragment_list, "\r\n");
  fragment_list = ValNodeFreeData (fragment_list);
  return txt;
}

typedef struct basegbqualeditor {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer         callback_data;
} BaseGBQualEditor, PNTR BaseGBQualEditorPtr;

static void ChangeGBQualEditorPopup (PopuP p)
{
  BaseGBQualEditorPtr dlg;

  dlg = (BaseGBQualEditorPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->tlp_callback != NULL) {
    (dlg->tlp_callback)(dlg->callback_data);
  }
}

static void ChangeGBQualEditorButton (ButtoN b)
{
  BaseGBQualEditorPtr dlg;

  dlg = (BaseGBQualEditorPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->tlp_callback != NULL) {
    (dlg->tlp_callback)(dlg->callback_data);
  }
}

static void ChangeGBQualEditorText (TexT t)
{
  BaseGBQualEditorPtr dlg;

  dlg = (BaseGBQualEditorPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->tlp_callback != NULL) {
    (dlg->tlp_callback)(dlg->callback_data);
  }
}


typedef struct collectiondatedlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  PopuP month;
  PopuP day;
  PopuP year;
  Int4  start_year;
} CollectionDateDlgData, PNTR CollectionDateDlgPtr;

static void PopulateDayPopup (PopuP p, Int4 month);

static void PointerToCollectionDateDialog (DialoG d, Pointer data)
{
  CollectionDateDlgPtr dlg;
  CharPtr              val, reformatted;
  Boolean              ambiguous = FALSE;
  CharPtr              cp, cp2, month_abbrev;
  Int4                 year = 1, month = 1, day = 1;

  dlg = (CollectionDateDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  SetValue (dlg->year, 1);
  SetValue (dlg->month, 1);
  SetValue (dlg->day, 1);
  Disable (dlg->month);
  Disable (dlg->day);

  val = (CharPtr) data;

  reformatted = ReformatDateStringEx (val, TRUE, &ambiguous);
  if (StringHasNoText (reformatted) || ambiguous) {
    /* do nothing */
  } else {
    cp = StringChr (reformatted, '-');
    if (cp == NULL) {
      year = GetYearFromToken (reformatted, StringLen (reformatted));
      year = year - dlg->start_year + 2;
    } else {
      if (isdigit (*reformatted)) {
        day = ReadNumberFromToken (reformatted, cp - reformatted);
        day += 1;
        cp++;
        cp2 = StringChr (cp, '-');
        month_abbrev = GetMonthFromToken (cp, cp2 - cp);
        month = GetMonthNumFromAbbrev (month_abbrev);
        if (month > -1) {
          month += 2;
        } else {
          month = 1;
        }
        year = GetYearFromToken (cp2 + 1, StringLen (cp2 + 1));
        year = year - dlg->start_year + 2;
      }
      else
      {
        month_abbrev = GetMonthFromToken (reformatted, cp - reformatted);
        month = GetMonthNumFromAbbrev (month_abbrev);
        if (month > -1) {
          month += 2;
        } else {
          month = 1;
        }
        year = GetYearFromToken (cp + 1, StringLen (cp + 1));
        year = year - dlg->start_year + 2;
      }
    }
    SetValue (dlg->year, year);
    Enable (dlg->month);
    SetValue (dlg->month, month);
    if (month > 1) {
      PopulateDayPopup (dlg->day, month);
      Enable (dlg->day);
    }
    SetValue (dlg->day, day);
  }  
}

static Pointer CollectionDateDialogToPointer (DialoG d)
{
  CollectionDateDlgPtr dlg;
  Int4                 year, month = -1, day = 0;
  CharPtr              year_fmt = "%d";
  CharPtr              mon_year_fmt = "%s-%d";
  CharPtr              day_mon_year_fmt = "%d-%s-%d";
  Char                 date_str[100];

  dlg = (CollectionDateDlgPtr) GetObjectExtra (d);
  date_str[0] = 0;
  if (dlg != NULL) {
    year = GetValue (dlg->year);
    if (year < 2) {
      return NULL;
    } else {
      year = year + dlg->start_year - 2;
      month = GetValue (dlg->month);
      if (month < 2) {
        sprintf (date_str, year_fmt, year);
      } else {
        day = GetValue (dlg->day);
        if (day < 2) {
          sprintf (date_str, mon_year_fmt, GetMonthAbbrev (month - 1), year);
        } else {
          sprintf (date_str, day_mon_year_fmt, day - 1, GetMonthAbbrev (month - 1), year);
        }
      }
    }
  }

  return StringSave (date_str);
}

static void PopulateDayPopup (PopuP p, Int4 month)
{
  Int4 i;
  Char day[10];

  Reset (p);
  PopupItem (p, "");

  if (month > 1) {
    for (i = 1; i <= GetDaysInMonth(month - 1); i++) {
      sprintf (day, "%d", i);
      PopupItem (p, day);
    }
  }
}

static void PopulateMonthPopup (PopuP p)
{
  Int4 i;

  PopupItem (p, "");
  for (i = 0; i < 12; i++) {
    PopupItem (p, GetMonthAbbrev(i + 1));
  }
}

static void ChangeCollectionDateMonth (PopuP p)
{
  CollectionDateDlgPtr dlg;
  Int4                 month;

  dlg = (CollectionDateDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  month = GetValue (p);
  if (month < 2) {
    Disable (dlg->day);
  } else {
    PopulateDayPopup (dlg->day, month);
    Enable (dlg->day);
  }
  if (dlg->tlp_callback != NULL) {
    (dlg->tlp_callback) (dlg->callback_data);
  }
}

static Int4 PopulateYearPopup (PopuP p)
{
	Nlm_DayTime dt;
  Char        year_str[20];
  Int4        start_year, i, end_year;

  GetDayTime (&dt);
  start_year = dt.tm_year + 1901 - 10;
  end_year = start_year + 20;
  PopupItem (p, " ");
  for (i=start_year; i <= end_year; i++) {
    sprintf (year_str, "%d", i);
    PopupItem (p, year_str);
  }
  return start_year;
}

static void ChangeCollectionDateYear (PopuP p)
{
  CollectionDateDlgPtr dlg;

  dlg = (CollectionDateDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  if (GetValue (p) < 2) {
    Disable (dlg->month);
    Disable (dlg->day);
  } else {
    Enable (dlg->month);
    ChangeCollectionDateMonth (dlg->month);
  }
  if (dlg->tlp_callback != NULL) {
    (dlg->tlp_callback) (dlg->callback_data);
  }
}

extern DialoG CollectionDateDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                               TaglistCallback tlp_callback,
                               Pointer callback_data)
{
  CollectionDateDlgPtr dlg;
  GrouP           p;

  p = HiddenGroup (h, 4, 0, NULL);
  dlg = (CollectionDateDlgPtr) MemNew (sizeof(CollectionDateDlgData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToCollectionDateDialog;
  dlg->fromdialog = CollectionDateDialogToPointer;
  dlg->testdialog = NULL;
  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;

  dlg->year = PopupList (p, TRUE, ChangeCollectionDateYear);
  SetObjectExtra (dlg->year, dlg, NULL);
  dlg->start_year = PopulateYearPopup (dlg->year);
  SetValue (dlg->year, 1);

  dlg->month = PopupList (p, TRUE, ChangeCollectionDateMonth);
  SetObjectExtra (dlg->month, dlg, NULL);
  PopulateMonthPopup (dlg->month);
  SetValue (dlg->month, 1);

  dlg->day = PopupList (p, TRUE, ChangeGBQualEditorPopup);
  SetObjectExtra (dlg->day, dlg, NULL);

  return (DialoG) p;
}

extern Boolean ParseCollectionDateOk (CharPtr txt)
{
  CharPtr date_str;
  Boolean ambiguous = FALSE;
  Boolean rval = FALSE;

  if (StringHasNoText (txt)) {
    rval = TRUE;
  } else {
    date_str = ReformatDateStringEx (txt, TRUE, &ambiguous);
    if (date_str != NULL && !ambiguous) {
      rval = TRUE;
    } 
    date_str = MemFree (date_str);
  }
  return rval;
}

typedef struct rptunitrangedlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  TexT  range_start;
  TexT  range_stop;
} RptUnitRangeDlgData, PNTR RptUnitRangeDlgPtr;

static Boolean ParseRptUnitRangeOkEx (CharPtr txt, Int4Ptr pstart, Int4Ptr pstop)
{
  CharPtr cp;
  Int4    start = -1, stop = -1;
  Char    ch;
  Boolean rval = FALSE;

  if (StringHasNoText (txt)) {
    return TRUE;
  }
  
  cp = txt;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (cp != txt) {
    ch = *cp;
    *cp = 0;
    start = atoi (txt);
    *cp = ch;
    cp++;
    while (*cp == ch) {
      cp++;
    }
    txt = cp;
    while (isdigit (*cp)) {
      cp++;
    }
    if (cp != txt) {
      stop = atoi (txt);
    }
  }

  if (start > -1 && stop > -1) {
    if (pstart != NULL) {
      *pstart = start;
    }
    if (pstop != NULL) {
      *pstop = stop;
    }
    rval = TRUE;
  }
  return rval;
}

static void StringToRptUnitRangeDialog (DialoG d, Pointer data)
{
  RptUnitRangeDlgPtr dlg;
  CharPtr            val;
  Int4               start = -1, stop = -1;
  Char               num_str[100];

  dlg = (RptUnitRangeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  val = (CharPtr) data;
  
  ParseRptUnitRangeOkEx (val, &start, &stop);
  if (start > -1 && stop > -1) {
    sprintf (num_str, "%d", start);
    SetTitle (dlg->range_start, num_str);
    sprintf (num_str, "%d", stop);
    SetTitle (dlg->range_stop, num_str);
  } else {
    SetTitle (dlg->range_start, "");
    SetTitle (dlg->range_stop, "");
  }
}

static Pointer RptUnitRangeDialogToString (DialoG d)
{
  RptUnitRangeDlgPtr dlg;
  CharPtr            val = NULL;
  CharPtr            start, stop;

  dlg = (RptUnitRangeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  start = SaveStringFromText (dlg->range_start);
  stop = SaveStringFromText (dlg->range_stop);

  if (StringHasNoText (start) && StringHasNoText (stop)) {
    return StringSave ("");
  } else {
    val = (CharPtr) MemNew (sizeof (Char) * (StringLen (start) + StringLen (stop) + 3));
    sprintf (val, "%s..%s", start == NULL ? "" : start, stop == NULL ? "" : stop);
  }
  start = MemFree (start);
  stop = MemFree (stop);
  return val;
}


extern DialoG CreateRptUnitRangeDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                         TaglistCallback tlp_callback,
                                         Pointer callback_data)
{
  RptUnitRangeDlgPtr dlg;
  GrouP           p;

  p = HiddenGroup (h, 6, 0, NULL);
  dlg = (RptUnitRangeDlgPtr) MemNew (sizeof(RptUnitRangeDlgData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToRptUnitRangeDialog;
  dlg->fromdialog = RptUnitRangeDialogToString;
  dlg->testdialog = NULL;

  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;

  StaticPrompt (p, "Start", 0, dialogTextHeight, programFont, 'l');
  dlg->range_start = DialogText (p, "", 5, ChangeGBQualEditorText);
  SetObjectExtra (dlg->range_start, dlg, NULL);
  StaticPrompt (p, "Stop", 0, dialogTextHeight, programFont, 'l');
  dlg->range_stop = DialogText (p, "", 5, ChangeGBQualEditorText);
  SetObjectExtra (dlg->range_stop, dlg, NULL);


  return (DialoG) p;
}

static Boolean ParseRptUnitRangeOk (CharPtr txt)
{
  return ParseRptUnitRangeOkEx (txt, NULL, NULL);
}


typedef struct controlplusfreedlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  PopuP element_type;
  TexT  description;

  CharPtr PNTR control_words;
} ControlPlusFreeDlgData, PNTR ControlPlusFreeDlgPtr;


static Int2 GetControlNum (CharPtr txt, CharPtr PNTR control_words, CharPtr PNTR desc_start)
{
  CharPtr cp;
  Int2    i;
  Int4    keyword_len;

  if (desc_start != NULL) {
    *desc_start = NULL;
  }
  if (StringHasNoText (txt) || control_words == NULL) {
    return 0;
  }
  /* skip over any leading spaces */
  while (isspace (*txt)) {
    txt++;
  }
  cp = StringChr (txt, ':');
  if (cp == NULL) {
    for (i = 1; control_words[i] != NULL; i++) {
      if (StringICmp (txt, control_words[i]) == 0) {
        return i;
      }
    }
    return -1;
  } else {
    keyword_len = cp - txt;
    while (keyword_len > 0 && isspace (txt[keyword_len - 1])) {
      keyword_len--;
    }
    if (keyword_len == 0) {
      return 0;
    }
    for (i = 1; control_words[i] != NULL; i++) {
      if (StringNICmp (txt, control_words[i], keyword_len) == 0) {
        if (desc_start != NULL && !StringHasNoText (cp + 1)) {
          *desc_start = cp + 1;
          while (isspace(**desc_start)) {
            (*desc_start)++;
          }
        }
        return i;
      }
    }
  }
  return -1;
}

static void ChangeControlPlusFreeType (PopuP p)
{
  ControlPlusFreeDlgPtr dlg;

  dlg = (ControlPlusFreeDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  if (GetValue (dlg->element_type) <= 1) {
    Disable (dlg->description);
  } else {
    Enable (dlg->description);
  }
  ChangeGBQualEditorPopup (p);
}

static void StringToControlPlusFreeDialog (DialoG d, Pointer data)
{
  ControlPlusFreeDlgPtr dlg;
  CharPtr             val;
  CharPtr             desc_start = NULL;
  Int2                num;

  dlg = (ControlPlusFreeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  val = (CharPtr) data;

  num = GetControlNum (val, dlg->control_words, &desc_start);

  if (num >= 0) {
    SetValue (dlg->element_type, num + 1);
    SetTitle (dlg->description, desc_start == NULL ? "" : desc_start);
  }
}

static Pointer ControlPlusFreeDialogToString (DialoG d)
{
  ControlPlusFreeDlgPtr dlg;
  CharPtr             val = NULL;
  CharPtr             desc_start = NULL;
  Int2                num;

  dlg = (ControlPlusFreeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  num = GetValue (dlg->element_type);

  if (num >= 1) {
    desc_start = SaveStringFromText (dlg->description);
    if (StringHasNoText (desc_start)) {
      val = StringSave (dlg->control_words[num - 1]);
    } else {
      val = (CharPtr) MemNew (sizeof (Char) * (StringLen (dlg->control_words[num - 1]) + StringLen (desc_start) + 2));
      sprintf (val, "%s:%s", dlg->control_words[num - 1], desc_start);
    }
  }
  return (Pointer) val;
}


static DialoG CreateControlPlusFreeDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                           TaglistCallback tlp_callback,
                                           Pointer callback_data,
                                           CharPtr PNTR control_words)
{
  ControlPlusFreeDlgPtr dlg;
  GrouP           p;
  Int2            i;

  p = HiddenGroup (h, 3, 0, NULL);
  dlg = (ControlPlusFreeDlgPtr) MemNew (sizeof(ControlPlusFreeDlgData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToControlPlusFreeDialog;
  dlg->fromdialog = ControlPlusFreeDialogToString;
  dlg->testdialog = NULL;

  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;

  dlg->control_words = control_words;
  dlg->element_type = PopupList (p, TRUE, ChangeControlPlusFreeType);
  for (i = 0; dlg->control_words[i] != NULL; i++) {
    PopupItem (dlg->element_type, dlg->control_words[i]);
  }
  SetValue (dlg->element_type, 1);
  SetObjectExtra (dlg->element_type, dlg, NULL);

  dlg->description = DialogText (p, "", 15, ChangeGBQualEditorText);
  SetObjectExtra (dlg->description, dlg, NULL);

  return (DialoG) p;
}

static void CopyTextToControlPlusFreeDialog (DialoG d, CharPtr txt)
{
  ControlPlusFreeDlgPtr dlg;

  dlg = (ControlPlusFreeDlgPtr) GetObjectExtra (d);

  if (dlg == NULL || StringHasNoText (txt)) return;

  SetTitle (dlg->description, txt);
}


CharPtr mobile_element_keywords[] = 
{ " ",
  "transposon",
  "retrotransposon",
  "integron",
  "insertion sequence", 
  "non-LTR retrotransposon",              
  "SINE", 
  "MITE", 
  "LINE", 
  "other",
  NULL};


static DialoG CreateMobileElementDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                         TaglistCallback tlp_callback,
                                         Pointer callback_data)
{
  return CreateControlPlusFreeDialog (h, sep, name, tlp_callback, callback_data, mobile_element_keywords);
}

static Boolean ParseMobileElementOk (CharPtr txt)
{

  if (GetControlNum (txt, mobile_element_keywords, NULL) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}


CharPtr satellite_keywords[] = 
{ " ",
  "satellite",
  "microsatellite",
  "minisatellite",
  NULL};


static DialoG CreateSatelliteDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                         TaglistCallback tlp_callback,
                                         Pointer callback_data)
{
  return CreateControlPlusFreeDialog (h, sep, name, tlp_callback, callback_data, satellite_keywords);
}

static Boolean ParseSatelliteOk (CharPtr txt)
{

  if (GetControlNum (txt, satellite_keywords, NULL) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}


typedef struct truefalsedlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  ButtoN is_true;
} TrueFalseDlgData, PNTR TrueFalseDlgPtr;

static void PointerToTrueFalseDialog (DialoG d, Pointer data)
{
  TrueFalseDlgPtr dlg;
  CharPtr         val;

  dlg = (TrueFalseDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  val = (CharPtr) data;

  if (StringICmp (val, "TRUE") == 0) {
    SetStatus (dlg->is_true, TRUE);
  } else {
    SetStatus (dlg->is_true, FALSE);
  }
}

static Pointer TrueFalseDialogToPointer (DialoG d)
{
  TrueFalseDlgPtr dlg;

  dlg = (TrueFalseDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return StringSave ("FALSE");

  if (GetStatus (dlg->is_true)) {
    return StringSave ("TRUE");
  } else {
    return StringSave ("FALSE");
  }
}

static DialoG TrueFalseDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                               TaglistCallback tlp_callback,
                               Pointer callback_data)
{
  TrueFalseDlgPtr dlg;
  GrouP           p;

  p = HiddenGroup (h, 2, 0, NULL);
  dlg = (TrueFalseDlgPtr) MemNew (sizeof(TrueFalseDlgData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToTrueFalseDialog;
  dlg->fromdialog = TrueFalseDialogToPointer;
  dlg->testdialog = NULL;
  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;

  dlg->is_true = CheckBox (p, "", ChangeGBQualEditorButton);
  SetObjectExtra (dlg->is_true, dlg, NULL);

  return (DialoG) p;
}

static Boolean ParseTrueFalseOk (CharPtr txt)
{
  if (StringHasNoText (txt) || StringICmp (txt, "TRUE") == 0 || StringICmp (txt, "FALSE") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}

typedef struct simplelistqualdlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  PopuP list;
  CharPtr PNTR val_list;
} SimpleListQualDlgData, PNTR SimpleListQualDlgPtr;

static Int2 GetListNum (CharPtr val, CharPtr PNTR val_list)
{
  Int2 i;
  if (StringHasNoText (val)) return 0;
  for (i = 1; val_list[i] != NULL; i++) {
    if (StringICmp (val, val_list[i]) == 0) {
      return i;
    }
  }
  return -1;
}

static void StringToSimpleListQualDialog (DialoG d, Pointer data)
{
  SimpleListQualDlgPtr dlg;
  CharPtr              val;
  Int2                 num;

  dlg = (SimpleListQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  val = (CharPtr) data;
  num = GetListNum (val, dlg->val_list);
  if (num > -1) {
    SetValue (dlg->list, num + 1);
  } else {
    SetValue (dlg->list, 1);
  }
}

static Pointer SimpleListQualDialogToString (DialoG d)
{
  SimpleListQualDlgPtr dlg;
  CharPtr              val = NULL;
  Int2                 num;

  dlg = (SimpleListQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  num = GetValue (dlg->list);
  if (num > 1) {
    val = StringSave (dlg->val_list[num - 1]);
  }
  return (Pointer) val;
}

static DialoG SimpleListQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name, CharPtr PNTR list,
                                    TaglistCallback tlp_callback,
                                    Pointer callback_data)
{
  SimpleListQualDlgPtr dlg;
  GrouP           p;
  Int2            i;

  p = HiddenGroup (h, 2, 0, NULL);
  dlg = (SimpleListQualDlgPtr) MemNew (sizeof(SimpleListQualDlgData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToSimpleListQualDialog;
  dlg->fromdialog = SimpleListQualDialogToString;
  dlg->testdialog = NULL;
  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;

  dlg->list = PopupList (p, TRUE, ChangeGBQualEditorPopup);
  SetObjectExtra (dlg->list, dlg, NULL);
  dlg->val_list = list;
  if (dlg->val_list != NULL) {
    for (i = 0; dlg->val_list[i] != NULL; i++) {
      PopupItem (dlg->list, dlg->val_list[i]);
    }
    SetValue (dlg->list, 1);
  }
  return (DialoG) p;
}

typedef struct multilistqualdlg {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback; 
  Pointer callback_data;
  ButtoN  *button_list;
  TexT    last_txt;

  CharPtr PNTR val_list;
} MultiListQualDlgData, PNTR MultiListQualDlgPtr;

static void CleanupMultiListQualDlg (GraphiC g, VoidPtr data)

{
  MultiListQualDlgPtr  dlg;

  dlg = (MultiListQualDlgPtr) data;
  if (dlg != NULL) {
    dlg->button_list = MemFree (dlg->button_list);
  }
  StdCleanupExtraProc (g, data);
}

static void CheckMultiListLast (ButtoN b)
{
  MultiListQualDlgPtr dlg;

  dlg = (MultiListQualDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (GetStatus (b)) {
    Enable (dlg->last_txt);
  } else {
    Disable (dlg->last_txt);
  }
  if (dlg->tlp_callback) {
    (dlg->tlp_callback) (dlg->callback_data);
  }
}

static void ApplyOneValToMultiListQualDialog (MultiListQualDlgPtr dlg, CharPtr val, CharPtr PNTR other_text)
{
  Int2    i;
  CharPtr tmp;

  if (dlg == NULL || StringHasNoText (val) || other_text == NULL) return;

  i = GetListNum (val, dlg->val_list);
  if (i > 0) {
    SetStatus (dlg->button_list[i - 1], TRUE);
  } else {
    if (*other_text == NULL) {
      *other_text = StringSave (val);
    } else {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (*other_text) + StringLen (val) + 2));
      sprintf (tmp, "%s,%s", *other_text, val);
      *other_text = MemFree (*other_text);
      *other_text = tmp;
    }
  }
}

static void StringToMultiListQualDialog (DialoG d, Pointer data)
{
  MultiListQualDlgPtr dlg;
  CharPtr             val, cp, other_text = NULL;
  Int2                i;

  dlg = (MultiListQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  for (i = 1; dlg->val_list[i] != NULL; i++) {
    SetStatus (dlg->button_list[i - 1], FALSE);
  }
  SetStatus (dlg->button_list[i - 1], FALSE);
  SetTitle (dlg->last_txt, "");

  val = (CharPtr) data;
  if (StringHasNoText (val)) {
    return;
  }
  
  if (*val == '(' && val[StringLen(val) - 1] == ')') {
    /* parentheses list */
    val++;
    cp = StringChr (val, ',');
    while (cp != NULL) {
      *cp = 0;
      ApplyOneValToMultiListQualDialog (dlg, val, &other_text);
      *cp = ',';
      val = cp + 1;
      while (isspace(*val)) {
        val++;
      }
      cp = StringChr (val, ',');
    }
    val[StringLen(val) - 1] = 0;
    ApplyOneValToMultiListQualDialog (dlg, val, &other_text);
    val[StringLen(val) - 1] = ')';
  } else {
    ApplyOneValToMultiListQualDialog (dlg, val, &other_text);
  }
  if (dlg->last_txt != NULL)
  {
    if (!StringHasNoText (other_text)) {
      SetTitle (dlg->last_txt, other_text);
      SetStatus (dlg->button_list[i - 1], TRUE);
      Enable (dlg->last_txt);
    } else {
      Disable (dlg->last_txt);
    }
  }
  other_text = MemFree (other_text);
}

static Pointer MultiListQualDialogToString (DialoG d)
{
  MultiListQualDlgPtr dlg;
  CharPtr             val = NULL, tmp;
  Int2                i;

  dlg = (MultiListQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  for (i = 1; dlg->val_list[i] != NULL; i++) {
    if (GetStatus (dlg->button_list[i - 1])) {
      if (val == NULL) {
        val = StringSave (dlg->val_list[i]);
      } else {
        tmp = StringSave (dlg->val_list[i]);
        val = CombineSplitGBQual (val, tmp);
      }
    }
  }
  if (dlg->last_txt != NULL
      && GetStatus (dlg->button_list[i - 1]) 
      && !TextHasNoText (dlg->last_txt)) {
    if (val == NULL) {
      val = SaveStringFromText (dlg->last_txt);
    } else {
      tmp = SaveStringFromText (dlg->last_txt);
      val = CombineSplitGBQual (val, tmp);
      tmp = MemFree (tmp);
    }
  }
  return val;
}

static DialoG MultiListQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name, CharPtr PNTR list,
                                   Boolean allow_not_in_list,
                                   TaglistCallback tlp_callback,
                                   Pointer callback_data)
{
  MultiListQualDlgPtr dlg;
  GrouP           p, g, last_button = NULL;
  Int2            i, num_buttons = 1; /* start with one, for the type-in wildcard value */

  p = HiddenGroup (h, -1, 0, NULL);
  dlg = (MultiListQualDlgPtr) MemNew (sizeof(MultiListQualDlgData));

  SetObjectExtra (p, dlg, CleanupMultiListQualDlg);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToMultiListQualDialog;
  dlg->fromdialog = MultiListQualDialogToString;
  dlg->testdialog = NULL;
  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;
  dlg->val_list = list;

  /* start with 1, skip blank at start of list */
  for (i = 1; dlg->val_list[i] != NULL; i++) {
    num_buttons ++;
  }

  dlg->button_list = (ButtoN *) MemNew (num_buttons * sizeof (ButtoN));

  g = HiddenGroup (p, 4, 0, NULL);
  for (i = 1; dlg->val_list [i] != NULL; i++) {
    dlg->button_list[i - 1] = CheckBox (g, dlg->val_list[i], ChangeGBQualEditorButton);
    SetObjectExtra (dlg->button_list[i - 1], dlg, NULL);
  }
  if (allow_not_in_list)
  {
    last_button = HiddenGroup (p, 2, 0, NULL);
    dlg->button_list[i - 1] = CheckBox (last_button, "", CheckMultiListLast);
    SetObjectExtra (dlg->button_list[i - 1], dlg, NULL);
    dlg->last_txt = DialogText (last_button, "", 20, ChangeGBQualEditorText);
    SetObjectExtra (dlg->last_txt, dlg, NULL);
    Disable (dlg->last_txt);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) last_button, NULL);
    
  return (DialoG) p;
}

static Boolean ParseMultiListQualOk (CharPtr val)
{
  return TRUE;
}

CharPtr codon_start_values[] = {" ", "1", "2", "3", NULL};
static DialoG CodonStartQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                    TaglistCallback tlp_callback,
                                    Pointer callback_data)

{
  return SimpleListQualDialog (h, sep, name, codon_start_values, tlp_callback, callback_data);
}
static Boolean ParseCodonStartOk (CharPtr val)
{
  if (GetListNum (val, codon_start_values) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}

CharPtr mol_type_values[] = {" ", "genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA",
                "snoRNA", "snRNA", "scRNA", "pre-RNA", "other RNA",
                "other DNA", "viral cRNA", "unassigned DNA", "unassigned RNA", NULL};
static DialoG MolTypeQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                 TaglistCallback tlp_callback,
                                 Pointer callback_data)
{
  return SimpleListQualDialog (h, sep, name, mol_type_values, tlp_callback, callback_data);
}
static Boolean ParseMolTypeOk (CharPtr val)
{
  if (GetListNum (val, mol_type_values) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}

CharPtr organelle_values[] = {" ", "mitochondrion", "nucleomorph", "plastid",
                              "mitochondrion:kinetoplast",
                              "plastid:chloroplast", "plastid:apicoplast", 
                              "plastid:chromoplast", "plastid:cyanelle", 
                              "plastid:leucoplast", "plastid:proplastid", NULL};
static DialoG OrganelleQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                   TaglistCallback tlp_callback,
                                   Pointer callback_data)
{
  return SimpleListQualDialog (h, sep, name, organelle_values, tlp_callback, callback_data);
}
static Boolean ParseOrganelleOk (CharPtr val)
{
  if (GetListNum (val, organelle_values) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}

CharPtr rpt_type_values[] = {" ", "tandem", "inverted", "flanking", "terminal", 
                             "direct", "dispersed", "other", NULL};
static DialoG RptTypeQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                 TaglistCallback tlp_callback,
                                 Pointer callback_data)
{


  return MultiListQualDialog (h, sep, name, rpt_type_values, FALSE, tlp_callback, callback_data);
}

static Boolean ParseRptTypeOk (CharPtr val)
{
  CharPtr cpy, cp, check;
  Boolean rval = TRUE;

  if (StringHasNoText (val))
  {
    return TRUE;
  }
  else if (val[0] == '(' && val[StringLen(val) - 1] == ')')
  {
    cpy = StringSave (val + 1);
    cpy [StringLen (cpy) - 1] = 0;
    check = cpy;
    cp = StringChr (cpy, ',');
    while (cp != NULL && rval)
    {
      *cp = 0;
      TrimSpacesAroundString (check);
      if (GetListNum (check, rpt_type_values) < 0)
      {
        rval = FALSE;
      }
      check = cp + 1;
      cp = StringChr (check, ',');
    }
    if (rval)
    {
      TrimSpacesAroundString (check);
      if (GetListNum (check, rpt_type_values) < 0)
      {
        rval = FALSE;
      }
    }
    cpy = MemFree (cpy);
  }
  else if (GetListNum (val, rpt_type_values) < 0)
  {
    rval = FALSE;
  }
  return rval;
}

CharPtr direction_values[] = {" ", "LEFT", "RIGHT", "BOTH", NULL};
static DialoG DirectionQualDialog (GrouP h, SeqEntryPtr sep, CharPtr name,
                                 TaglistCallback tlp_callback,
                                 Pointer callback_data)
{
  return SimpleListQualDialog (h, sep, name, direction_values, tlp_callback, callback_data);
}
static Boolean ParseDirectionOk (CharPtr val)
{
  if (GetListNum (val, direction_values) > -1) {
    return TRUE;
  } else {
    return FALSE;
  }
}

typedef DialoG (*BuildQualEditorDialog) PROTO ((GrouP, SeqEntryPtr, CharPtr, 
                                                TaglistCallback, Pointer));
typedef void (*CopyTextToEditor) PROTO ((DialoG, CharPtr));

typedef struct gbqualeditlist {
  CharPtr               name;
  BuildQualEditorDialog build_dlg;
  ParseOK               parse_ok;
  CopyTextToEditor      copy_txt;
} GBQualEditListData, PNTR GBQualEditListPtr;

GBQualEditListData gbqual_edit_list[] = {
  {"chloroplast",          TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"chromoplast",          TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"codon_start",          CodonStartQualDialog,      ParseCodonStartOk,     NULL },
  {"collection_date",      CollectionDateDialog,      ParseCollectionDateOk, NULL },
  {"cyanelle",             TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"direction",            DirectionQualDialog,       ParseDirectionOk,      NULL },
  {"germline",             TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"kinetoplast",          TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"macronuclear",         TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"mitochondrion",        TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"mobile_element",       CreateMobileElementDialog, ParseMobileElementOk,  CopyTextToControlPlusFreeDialog },
  {"mol_type",             MolTypeQualDialog,         ParseMolTypeOk,        NULL },
  {"organelle",            OrganelleQualDialog,       ParseOrganelleOk,      NULL },
  {"partial",              TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"proviral",             TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"pseudo",               TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"rearranged",           TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"rpt_type",             RptTypeQualDialog,         ParseRptTypeOk,        NULL },
  {"rpt_unit_range",       CreateRptUnitRangeDialog,  ParseRptUnitRangeOk,   NULL } ,
  {"satellite",            CreateSatelliteDialog,     ParseSatelliteOk,      CopyTextToControlPlusFreeDialog },
  {"virion",               TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"focus",                TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"transgenic",           TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"environmental_sample", TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"ribosomal_slippage",   TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"trans_splicing",       TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {"metagenomic",          TrueFalseDialog,           ParseTrueFalseOk,      NULL },
  {NULL, NULL, NULL}};


typedef struct singlegbqualedit {
  DIALOG_MESSAGE_BLOCK
  TaglistCallback tlp_callback;
  Pointer callback_data;
  DialoG editor;
  TexT   txt;
  GrouP  choice_grp;
  ButtoN            use_editor;
  ButtoN            use_text;
  ButtoN            copy_text;
  GBQualEditListPtr edit_list;
  CharPtr           name;
  Boolean           force_string;
} SingleGBQualEditData, PNTR SingleGBQualEditPtr;

static void ChangeSingleQualEdit (GrouP g)
{
  SingleGBQualEditPtr dlg;

  dlg = (SingleGBQualEditPtr) GetObjectExtra (g);
  if (dlg != NULL) {
    if (GetValue (dlg->choice_grp) == 1) {
      SafeEnable (dlg->editor);
      Disable (dlg->txt);
    } else {
      SetValue (dlg->choice_grp, 2);
      SafeDisable (dlg->editor);
      Enable (dlg->txt);
    }
    if (dlg->tlp_callback != NULL) {
      (dlg->tlp_callback) (dlg->callback_data);
    }
  }
}

static void StringToSingleGBQualEditDialog (DialoG d, Pointer data)
{
  SingleGBQualEditPtr dlg;
  CharPtr             val;

  dlg = (SingleGBQualEditPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  val = (CharPtr) data;

  if (StringHasNoText (val) || (!dlg->force_string && dlg->edit_list != NULL && (dlg->edit_list->parse_ok(val)))) {
    SetTitle (dlg->txt, "");
    if (dlg->choice_grp != NULL) {
      SetValue (dlg->choice_grp, 1);
      Hide (dlg->use_editor);
      Hide (dlg->use_text);
      Hide (dlg->txt);
      Hide (dlg->copy_text);
      PointerToDialog (dlg->editor, val);
      Enable (dlg->editor);
      Disable (dlg->txt);
    }
  } else {
    SetTitle (dlg->txt, val);
    if (dlg->choice_grp != NULL) {
      SetValue (dlg->choice_grp, 2);
      Show (dlg->use_editor);
      Show (dlg->use_text);
      Show (dlg->txt);
      Show (dlg->copy_text);
      PointerToDialog (dlg->editor, NULL);
      Enable (dlg->txt);
      Disable (dlg->editor);
    }
  }
}

static Pointer SingleGBQualEditDialogToString (DialoG d)
{
  SingleGBQualEditPtr dlg;
  CharPtr             val;

  dlg = (SingleGBQualEditPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  if (dlg->choice_grp == NULL || GetValue (dlg->choice_grp) == 2) {
    val = SaveStringFromText (dlg->txt);
  } else {
    val = DialogToPointer (dlg->editor);
  }
  
  return val;
}

static void GBQualEditDialogMessage (DialoG d, Int2 mssg)

{
  SingleGBQualEditPtr dlg;

  dlg = (SingleGBQualEditPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == VIB_MSG_INIT 
        || mssg == VIB_MSG_ENTER) {
      if (dlg->choice_grp == NULL || GetValue (dlg->choice_grp) == 2) {
        Select (dlg->txt);
      } else {
        Select (dlg->editor);
      }
    }
  }
}

static void GBQualEditDialogCopyText (ButtoN b)
{
  SingleGBQualEditPtr dlg;
  CharPtr             txt;

  dlg = (SingleGBQualEditPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->edit_list != NULL && dlg->edit_list->copy_txt != NULL) {
    txt = SaveStringFromText (dlg->txt);
    (dlg->edit_list->copy_txt)(dlg->editor, txt);
    txt = MemFree (txt);
  }
}

static DialoG CreateSingleGBQualEditDialog (GrouP h, SeqEntryPtr sep, CharPtr name, 
                                            Boolean force_string,
                                            TaglistCallback tlp_callback, 
                                            Pointer callback_data)
{
  SingleGBQualEditPtr dlg;
  GrouP           p;
  Int4            i;

  p = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  dlg = (SingleGBQualEditPtr) MemNew (sizeof(SingleGBQualEditData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToSingleGBQualEditDialog;
  dlg->fromdialog = SingleGBQualEditDialogToString;
  dlg->testdialog = NULL;
  dlg->dialogmessage = GBQualEditDialogMessage;
  dlg->tlp_callback = tlp_callback;
  dlg->callback_data = callback_data;
  dlg->force_string = force_string;

  dlg->name = name;
  dlg->edit_list = NULL;
  for (i = 0; gbqual_edit_list[i].name != NULL && dlg->edit_list == NULL; i++) {
    if (StringICmp (gbqual_edit_list[i].name, name) == 0) {
      dlg->edit_list = gbqual_edit_list + i;
    }
  }

  dlg->copy_text = NULL;

  if (dlg->edit_list == NULL || dlg->edit_list->build_dlg == NULL || force_string) {
    dlg->txt = DialogText (p, "", 20, ChangeGBQualEditorText);
    SetObjectExtra (dlg->txt, dlg, NULL);
    dlg->editor = NULL;
    dlg->choice_grp = NULL;
  } else {
    dlg->choice_grp = HiddenGroup (p, 0, 2, ChangeSingleQualEdit);
    SetObjectExtra (dlg->choice_grp, dlg, NULL);
    dlg->use_editor = RadioButton (dlg->choice_grp, "");
    dlg->use_text = RadioButton (dlg->choice_grp, "");
    dlg->editor = (dlg->edit_list->build_dlg) (dlg->choice_grp, sep, name, tlp_callback, callback_data);
    dlg->txt = DialogText (dlg->choice_grp, "", 20, ChangeGBQualEditorText);
    SetObjectExtra (dlg->txt, dlg, NULL);
    SetValue (dlg->choice_grp, 1);
    if (dlg->edit_list->copy_txt != NULL) {
      dlg->copy_text = PushButton (dlg->choice_grp, "Copy Text", GBQualEditDialogCopyText);
      SetObjectExtra (dlg->copy_text, dlg, NULL);
    }
//    AlignObjects (ALIGN_LOWER, (HANDLE) dlg->use_text, (HANDLE) dlg->txt, NULL);
    AlignObjects (ALIGN_LOWER, (HANDLE) dlg->txt, (HANDLE) dlg->use_text, NULL);
    ChangeSingleQualEdit (dlg->choice_grp);
  }

  return (DialoG) p;
}

typedef struct newfieldpage {
  DIALOG_MESSAGE_BLOCK
  Int2               numfields;
  DialoG             PNTR editors;

  GBQualPtr          new_gbq;
  Uint1              subtype;
  Boolean            allowProductGBQual;
  Int2               last_viewed;
} NewFieldPage, PNTR NewFieldPagePtr;

static void CleanupNewFieldsPage (GraphiC g, VoidPtr data)

{
  NewFieldPagePtr  fpf;

  fpf = (NewFieldPagePtr) data;
  if (fpf != NULL) {
    MemFree (fpf->editors);
    fpf->new_gbq = GBQualFree (fpf->new_gbq);
  }
  MemFree (data);
}

static CharPtr GetDisplayQualName (CharPtr qual_name, Boolean is_legal)
{
  CharPtr display_name;

  if (is_legal) {
    display_name = StringSave (qual_name);
  } else {
    display_name = (CharPtr) MemNew (sizeof (Char) * (StringLen (qual_name) + 2));
    sprintf (display_name, "*%s", qual_name);
  }
  return display_name;
}

static Boolean QualNamesMatch (CharPtr qual_name1, CharPtr qual_name2)
{
  if (StringICmp (qual_name1, qual_name2) == 0) {
    return TRUE;
  } else if (qual_name1 != NULL 
             && *(qual_name1) == '*' 
             && StringICmp (qual_name1 + 1, qual_name2) == 0) {
    return TRUE;
  } else if (qual_name2 != NULL
             && *qual_name2 == '*'
             && StringICmp (qual_name1, qual_name2 + 1) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Boolean IsCombinableQual (Int2 qual)
{
  Boolean   combine_qual = FALSE;

  if (qual == GBQUAL_rpt_type 
      || qual == GBQUAL_rpt_unit 
      || qual == GBQUAL_rpt_unit_range
      || qual == GBQUAL_rpt_unit_seq 
      || qual == GBQUAL_replace 
      || qual == GBQUAL_compare
      || qual == GBQUAL_old_locus_tag
      || qual == GBQUAL_usedin) {
    combine_qual = TRUE;
  }
  return combine_qual;
}

static void AddTemporaryGBQual (GBQualPtr PNTR gbq_list, CharPtr qual_name, GBQualPtr gbq_feat, Boolean is_legal)
{
  Int2 qual;
  GBQualPtr gbq_last = NULL;
  GBQualPtr gbq_new = NULL;
  Boolean   found = FALSE;
  Boolean   combine_qual = FALSE;
  CharPtr   blank_val = "\"\"";

  if (gbq_list == NULL) return;

  gbq_last = *gbq_list;
  while (gbq_last != NULL && gbq_last->next != NULL) {
    if (QualNamesMatch (gbq_last->qual, qual_name)) {
      /* already added */
      return;
    }
    gbq_last = gbq_last->next;
  }
  if (gbq_last != NULL && QualNamesMatch (gbq_last->qual, qual_name)) {
    /* already added */
    return;
  }

  /* add value.  If none in feature, if true/false, add TRUE or FALSE, else use blank. */
  qual = GBQualNameValid (qual_name);

  combine_qual = IsCombinableQual (qual);

  while (gbq_feat != NULL) {
    if (StringICmp (qual_name, gbq_feat->qual) == 0) {
      if (StringHasNoText (gbq_feat->val)) {
        if (qual > -1 && ParFlat_GBQual_names [qual].gbclass == Class_none && !found) {
          /* if this is a true-false, only add one qualifier */
          gbq_new = GBQualNew ();
          gbq_new->qual = GetDisplayQualName (qual_name, is_legal);
          gbq_new->val = StringSave ("TRUE");
          if (gbq_last == NULL) {
            *gbq_list = gbq_new;
          } else {
            gbq_last->next = gbq_new;
          }
          gbq_last = gbq_new;
          found = TRUE;
        } else if (qual == GBQUAL_replace) {
          /* save blank values */
          gbq_new = NULL;
          if (found && combine_qual) {
            gbq_new = *gbq_list;
            while (gbq_new != NULL && !QualNamesMatch (gbq_feat->qual, gbq_new->qual)) {
              gbq_new = gbq_new->next;
            }
          }
          if (gbq_new == NULL) {
            /* make new qualifier */
            gbq_new = GBQualNew ();
            gbq_new->qual = GetDisplayQualName (qual_name, is_legal);
            gbq_new->val = StringSave (blank_val);
            /* add to list */
            if (gbq_last == NULL) {
              *gbq_list = gbq_new;
            } else {
              gbq_last->next = gbq_new;
            }
            gbq_last = gbq_new;
          } else {
            /* combine with previous value */
            gbq_new->val = CombineSplitGBQual (gbq_new->val, blank_val);
          }
          found = TRUE;
        } else {
          /* if not true-false or already have true-false, we'll just be adding a blank value later */    
          /* so do nothing here */
        }
      } else {
        gbq_new = NULL;
        if (found && combine_qual) {
          gbq_new = *gbq_list;
          while (gbq_new != NULL && !QualNamesMatch (gbq_feat->qual, gbq_new->qual)) {
            gbq_new = gbq_new->next;
          }
        }
        if (gbq_new == NULL) {
          /* make new qualifier */
          gbq_new = GBQualNew ();
          gbq_new->qual = GetDisplayQualName (qual_name, is_legal);
          gbq_new->val = StringSave (gbq_feat->val);
          /* add to list */
          if (gbq_last == NULL) {
            *gbq_list = gbq_new;
          } else {
            gbq_last->next = gbq_new;
          }
          gbq_last = gbq_new;
        } else {
          /* combine with previous value */
          gbq_new->val = CombineSplitGBQual (gbq_new->val, gbq_feat->val);
        }
        found = TRUE;
      }
    }
    gbq_feat = gbq_feat->next;
  }
   
  if (!found) {
    gbq_new = GBQualNew ();
    gbq_new->qual = GetDisplayQualName (qual_name, is_legal);
    if (qual > -1 && ParFlat_GBQual_names [qual].gbclass == Class_none) {
      gbq_new->val = StringSave ("FALSE");
    } else {
      gbq_new->val = StringSave ("");
    }
    if (gbq_last == NULL) {
      *gbq_list = gbq_new;
    } else {
      gbq_last->next = gbq_new;
    }
  }
}

static Boolean IsRarelyUsed (Int2 qual)
{
  if (qual == GBQUAL_allele
      || qual == GBQUAL_function
      || qual == GBQUAL_map
      || qual == GBQUAL_standard_name
      || qual == GBQUAL_old_locus_tag) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Pointer NewDialogToImportFields (DialoG d)
{
  NewFieldPagePtr fpf;
  GBQualPtr       gbq_list = NULL, gbq_last = NULL, gbq_new, gbq_it;
  Int2            qual, i;
  CharPtr         val, qual_name;

  fpf = (NewFieldPagePtr) GetObjectExtra (d);
  if (fpf == NULL) return NULL;

  for (gbq_it = fpf->new_gbq, i = 0; gbq_it != NULL; gbq_it = gbq_it->next, i++) {
    gbq_it->val = MemFree (gbq_it->val);
    gbq_it->val = DialogToPointer (fpf->editors[i]);
  }

  for (gbq_it = fpf->new_gbq; gbq_it != NULL; gbq_it = gbq_it->next) {
    if (StringHasNoText (gbq_it->val)) {
      continue;
    }
    qual_name = gbq_it->qual;
    if (qual_name != NULL && *qual_name == '*') {
      qual_name++;
    }
    qual = GBQualNameValid (qual_name);

    val = gbq_it->val;
    if (qual > -1 && ParFlat_GBQual_names [qual].gbclass == Class_none) {
      if (StringICmp (val, "TRUE") == 0) {
        /* don't put "true" in qual, just use empty string */
        val = "";
      } else if (StringICmp (val, "FALSE") == 0) {
        /* don't add FALSE qual */
        continue;
      }
      /* any other values, add as they are */
    } 

    gbq_new = GBQualNew();
    gbq_new->qual = StringSave (qual_name);
    gbq_new->val = StringSave (val);
    if (gbq_last == NULL) {
      gbq_list = gbq_new;
    } else {
      gbq_last->next = gbq_new;
    }
    gbq_last = gbq_new;
  }
  return (Pointer) gbq_list;
}

static void 
AddMandatoryAndOptionalQuals 
(CharPtr         name,
 NewFieldPagePtr fpf,
 GBQualPtr       gbq, 
 Boolean         use_rarely_used)
{
  Int2 index, i, qual;
  SematicFeatPtr  sefp;

  if (name == NULL) return;

  index = GBFeatKeyNameValid (&name, FALSE);
  if (index < 0) return;

  sefp = &(ParFlat_GBFeat [index]);
  /* add mandatory quals first */
  for (i = 0; i < sefp->mand_num; i++) {
    qual = sefp->mand_qual [i];
    if (qual > -1 
        && ShouldBeAGBQual (fpf->subtype, qual, fpf->allowProductGBQual)
        && ((use_rarely_used && IsRarelyUsed (qual))
            || (!use_rarely_used && ! IsRarelyUsed (qual)))) {
      AddTemporaryGBQual (&(fpf->new_gbq), ParFlat_GBQual_names [qual].name, gbq, TRUE);
    }
  }
  /* add optional quals next */
  for (i = 0; i < sefp->opt_num; i++) {
    qual = sefp->opt_qual [i];
    if (qual > -1 
        && ShouldBeAGBQual (fpf->subtype, qual, fpf->allowProductGBQual)
        && ((use_rarely_used && IsRarelyUsed (qual))
            || (!use_rarely_used && ! IsRarelyUsed (qual)))) {
      AddTemporaryGBQual (&(fpf->new_gbq), ParFlat_GBQual_names [qual].name, gbq, TRUE);
    }
  }
}


extern DialoG NewCreateImportFields (GrouP h, CharPtr name, SeqFeatPtr sfp, Boolean allowProductGBQual)
{
  NewFieldPagePtr fpf;
  GrouP           g, g1, g2 = NULL;
  GBQualPtr       gbq = NULL, gbq_it;
  Int2            j;
  Int2            max;
  Int2            num = 0, num_legal = 0;
  GrouP           p;
  Int2            qual;
  Int2            wid;
  PrompT          ill_q = NULL;

  p = HiddenGroup (h, -1, 0, NULL);
  fpf = (NewFieldPagePtr) MemNew (sizeof (NewFieldPage));
  if (fpf != NULL) {

    SetObjectExtra (p, fpf, CleanupNewFieldsPage);
    fpf->dialog = (DialoG) p;
    fpf->todialog = NULL;
    fpf->fromdialog = NewDialogToImportFields;
    fpf->testdialog = NULL;

    fpf->allowProductGBQual = allowProductGBQual;
    if (sfp == NULL) {
      fpf->subtype = FEATDEF_ANY;
    } else {
      fpf->subtype = sfp->idx.subtype;
      gbq = sfp->qual;
    }

    fpf->last_viewed = -1;

    /* create list of temporary GBQuals */
    /* list mandatory/optional quals that are not rarely used first */
    AddMandatoryAndOptionalQuals (name, fpf, gbq, FALSE);

    /* add rarely used quals here */
    AddMandatoryAndOptionalQuals (name, fpf, gbq, TRUE);

    /* count legal values - all others to be listed as simple strings */
    for (gbq_it = fpf->new_gbq; gbq_it != NULL; gbq_it = gbq_it->next) {
      num_legal++;
    }

    /* add legal added qualifiers next */
    gbq_it = gbq;
    while (gbq_it != NULL) {
      qual = GBQualNameValid (gbq_it->qual);
      if (qual > -1 && ShouldBeAGBQual (fpf->subtype, qual, fpf->allowProductGBQual)) {
        AddTemporaryGBQual (&(fpf->new_gbq), ParFlat_GBQual_names [qual].name, gbq, FALSE);
      } 
      gbq_it = gbq_it->next;
    }

    /* add remaining qualifiers last */
    gbq_it = gbq;
    while (gbq_it != NULL) {  
      if (!ShouldSuppressGBQual(fpf->subtype, gbq_it->qual)) {
        AddTemporaryGBQual (&(fpf->new_gbq), gbq_it->qual, gbq, FALSE);
      }
      gbq_it = gbq_it->next;
    }

    /* calculate maximum name width */
    SelectFont (systemFont);
    max = 0;
    for (gbq_it = fpf->new_gbq; gbq_it != NULL; gbq_it = gbq_it->next) {
      num++;
      wid = StringWidth (gbq_it->qual) + 2;
      if (wid > max) {
        max = wid;
      }
    }
    SelectFont (systemFont);

    fpf->numfields = num;

    if (num > 0) {
      g1 = HiddenGroup (p, 2, 0, NULL);     
      g = g1;
      fpf->editors = MemNew (sizeof (DialoG) * (num + 1));
      for (gbq_it = fpf->new_gbq, j = 0; gbq_it != NULL; gbq_it = gbq_it->next) {
        if (j == num_legal) {
          ill_q = StaticPrompt (p, "*** Illegal Qualifiers ***", 0, dialogTextHeight, programFont, 'l');
          g2 = HiddenGroup (p, 2, 0, NULL);
          g = g2;
        }
        StaticPrompt (g, gbq_it->qual, 0, dialogTextHeight, programFont, 'l');
        fpf->editors[j] = CreateSingleGBQualEditDialog (g, NULL, gbq_it->qual,
                                                        j < num_legal ? FALSE : TRUE,
                                                        NULL, NULL);
        PointerToDialog (fpf->editors[j], gbq_it->val);
        j++;
      }   
      AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) ill_q, (HANDLE) g2, NULL);  
    } else {
      fpf->editors = NULL;
      StaticPrompt (p, "See Attributes page to set legal qualifiers for this feature.",
                    0, 0, programFont, 'c');
    }
  }
  return (DialoG) p;
}


typedef struct specialcharacterdialog 
{
  TexT PNTR       text_list;
  ButtoN          accept_btn;
  Int4            num_chars;
  ValNodePtr      find_list;
} SpecialCharacterDialogData, PNTR SpecialCharacterDialogPtr;  


static void EnableSpecialCharacterAccept (TexT t)
{
  SpecialCharacterDialogPtr sd;
  Int4                      pos;
  CharPtr                   str;
  Boolean                   has_bad = FALSE;

  sd = (SpecialCharacterDialogPtr) GetObjectExtra (t);
  if (sd == NULL) return;

  for (pos = 0; pos < sd->num_chars && !has_bad; pos++)
  {
    str = SaveStringFromText (sd->text_list[pos]);
    SpecialCharFind (&str, NULL, &has_bad, NULL);
    str = MemFree (str);
  }
  if (has_bad) 
  {
    Disable (sd->accept_btn);
  }
  else
  {
    Enable (sd->accept_btn);
  }
}


static void SetWindowsSpecialCharacterDefaults (ButtoN b)
{
  SpecialCharacterDialogPtr sd;
  Int4                      pos;
  CharPtr                   str;
  ValNodePtr                vnp;

  sd = (SpecialCharacterDialogPtr) GetObjectExtra (b);
  if (sd == NULL) return;

  for (vnp = sd->find_list, pos = 0; vnp != NULL; vnp = vnp->next, pos++)
  {
    str = GetSpecialWinCharacterReplacement ((unsigned char) vnp->choice);
    SetTitle (sd->text_list[pos], str);
  }
}


static void SetMacSpecialCharacterDefaults (ButtoN b)
{
  SpecialCharacterDialogPtr sd;
  Int4                      pos;
  CharPtr                   str;
  ValNodePtr                vnp;

  sd = (SpecialCharacterDialogPtr) GetObjectExtra (b);
  if (sd == NULL) return;

  for (vnp = sd->find_list, pos = 0; vnp != NULL; vnp = vnp->next, pos++)
  {
    str = GetSpecialMacCharacterReplacement ((unsigned char) vnp->choice);
    SetTitle (sd->text_list[pos], str);
  }

}


extern Boolean FixSpecialCharactersForStringsInList (ValNodePtr find_list, CharPtr exp_text, Boolean force_fix)
{
  ValNodePtr      vnp, context_list, vnp_c;
  Int4            pos;
  CharPtr         repl;
  WindoW          w;
  GrouP           h, c, g, p1, g2;
  PrompT          p2;
  LisT            contexts;
  ButtoN          b;
  SpecialCharacterDialogData sd;
  ModalAcceptCancelData      acd;
  Char                  label[2];
  Int4                  num_contexts;
  Boolean               rval = FALSE;

  if (find_list == NULL) {
    return TRUE;
  } 

  ArrowCursor();
  Update();

  sd.find_list = find_list;
  sd.num_chars = ValNodeLen (find_list);
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  p1 = MultiLinePrompt (h, exp_text, 27 * stdCharWidth, programFont);
  p2 = StaticPrompt (h, "Choose replacement characters.", 0, 0, programFont, 'l');

  g = HiddenGroup (h, 3, 0, NULL);
  StaticPrompt (g, "Character", 0, 0, programFont, 'l');
  StaticPrompt (g, "Replacement", 0, 0, programFont, 'l');
  StaticPrompt (g, "Contexts", 0, 0, programFont, 'l');

  sd.text_list = MemNew (sizeof (TexT) * sd.num_chars);
  label[1] = 0;
  for (vnp = find_list, pos = 0; vnp != NULL; vnp = vnp->next, pos++)
  {
    label[0] = vnp->choice;
    StaticPrompt (g, label, 0, 0, programFont, 'l');
    repl = GetSpecialCharacterReplacement (vnp->choice);
    sd.text_list[pos] = DialogText (g, repl, 5, EnableSpecialCharacterAccept);
    SetObjectExtra (sd.text_list[pos], &sd, NULL);
    context_list = vnp->data.ptrvalue;
    num_contexts = ValNodeLen (context_list);
    if (num_contexts == 0)
    {
      StaticPrompt (g, "", 0, 0, programFont, 'l');
    }
    else if (num_contexts == 1)
    {
      StaticPrompt (g, *((CharPtr PNTR)context_list->data.ptrvalue), 0, 0, programFont, 'l');
    }
    else
    {
      contexts = SingleList (g, 16, MIN (num_contexts, 5), NULL);
      for (vnp_c = context_list; vnp_c != NULL; vnp_c = vnp_c->next)
      {
        ListItem (contexts, *((CharPtr PNTR)vnp_c->data.ptrvalue));
      }
      SetValue (contexts, 1);
    }
  }

  g2 = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g2, 10, 10);
  b = PushButton (g2, "Suggest Windows Replacements", SetWindowsSpecialCharacterDefaults);
  SetObjectExtra (b, &sd, NULL);
  b = PushButton (g2, "Suggest Mac Replacements", SetMacSpecialCharacterDefaults);
  SetObjectExtra (b, &sd, NULL);

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  sd.accept_btn = PushButton (c, "Replace Characters", ModalAcceptButton);
  SetObjectExtra (sd.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p1, (HANDLE) p2, (HANDLE) g, (HANDLE) g2, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.accepted)
  {
    for (vnp = find_list, pos = 0; vnp != NULL; vnp = vnp->next, pos++)
    {
      repl = JustSaveStringFromText (sd.text_list[pos]);
      label[0] = vnp->choice;
      for (vnp_c = vnp->data.ptrvalue; vnp_c != NULL; vnp_c = vnp_c->next)
      {
        FindReplaceString (vnp_c->data.ptrvalue, label, repl, TRUE, FALSE);
      }
      repl = MemFree (repl);
    }

    rval = TRUE;
  }
  else if (force_fix)
  {
    for (vnp = find_list; vnp != NULL; vnp = vnp->next)
    {
      label[0] = vnp->choice;
      for (vnp_c = vnp->data.ptrvalue; vnp_c != NULL; vnp_c = vnp_c->next)
      {
        FindReplaceString (vnp_c->data.ptrvalue, label, "#", TRUE, FALSE);
      }
    }
    rval = TRUE;
  }

  Remove (w);
  return rval;

}


NLM_EXTERN Boolean FixSpecialCharacters (Uint2 entityID)
{
  ValNodePtr      find_list = NULL;
  Boolean         rval;

  StringActionInEntity (entityID, FALSE, UPDATE_NEVER, NULL, NULL, NULL, TRUE,
                        SpecialCharFindWithContext, NULL, &find_list);

  rval = FixSpecialCharactersForStringsInList (find_list,
              "The ASN.1 contains special characters.\nIf you save it without stripping them,\nyou will be unable to load the file in Sequin.",
              FALSE);

  find_list = FreeContextList (find_list);
  if (rval) 
  {
    ObjMgrSetDirtyFlag (entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  }
  return rval;
}


NLM_EXTERN Boolean FixSpecialCharactersForObject (Uint2 datatype, Pointer objdata, CharPtr msg, Boolean force_fix, BoolPtr changed)
{
  ValNodePtr      find_list = NULL;
  Boolean         rval;

  StringActionForObject (datatype, objdata, 0, FALSE, UPDATE_NEVER,
                        SpecialCharFindWithContext, NULL, &find_list);

  rval = FixSpecialCharactersForStringsInList (find_list,
              msg == NULL ?
                     "The ASN.1 contains special characters.\nIf you save it without stripping them,\nyou will be unable to load the file in Sequin."
                     : msg,
              force_fix);

  if (changed != NULL)
  {
    if (find_list == NULL)
    {
      *changed = FALSE;
    }
    else
    {
      *changed = rval;
    }
  } 

  find_list = FreeContextList (find_list);
  return rval;
}




typedef struct featurereplace {
  SeqFeatPtr sfp;
  Boolean    chosen;
  ValNodePtr replace_list;
} FeatureReplaceData, PNTR FeatureReplacePtr;


static FeatureReplacePtr FeatureReplaceNew (SeqFeatPtr sfp)
{
  FeatureReplacePtr f;

  f = (FeatureReplacePtr) MemNew (sizeof (FeatureReplaceData));
  f->sfp = AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
  f->chosen = TRUE;
  f->replace_list = NULL;
  return f;
}


static FeatureReplacePtr FeatureReplaceCopy (FeatureReplacePtr orig)
{
  FeatureReplacePtr f = NULL;
  

  if (orig != NULL) {
    f = FeatureReplaceNew (orig->sfp);
    f->chosen = orig->chosen;
    f->replace_list = ValNodeCopyPtr (orig->replace_list);
  }
  return f;
}


static FeatureReplacePtr FeatureReplaceFree (FeatureReplacePtr f)
{
  if (f != NULL) {
    f->sfp = SeqFeatFree (f->sfp);
    f->replace_list = ValNodeFree (f->replace_list);
    f = MemFree (f);
  }
  return f;
}


NLM_EXTERN ValNodePtr FeatureReplaceListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = FeatureReplaceFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static ValNodePtr FeatureReplaceListCopy (ValNodePtr orig)
{
  ValNodePtr list = NULL;

  while (orig != NULL) {
    ValNodeAddPointer (&list, orig->choice, FeatureReplaceCopy (orig->data.ptrvalue));
    orig = orig->next;
  }
  return list;
}


NLM_EXTERN ValNodePtr FeatureReplaceListFromSeqAnnot (SeqAnnotPtr sap)
{
  SeqFeatPtr sfp;
  ValNodePtr list = NULL;

  if (sap == NULL || sap->type != 1) {
    return NULL;
  }

  sfp = sap->data;
  while (sfp != NULL) {
    ValNodeAddPointer (&list, 0, FeatureReplaceNew (sfp));
    sfp = sfp->next;
  }
  return list;
}


static Boolean DoFeaturesAutomatchForUpdate (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  GeneRefPtr grp1, grp2;
  Boolean    rval = FALSE;
  GBQualPtr  gbq;
  SeqIdPtr   sip;
  BioseqPtr  pbsp;

  if (sfp1 == NULL || sfp2 == NULL) {
    rval = FALSE;
  } else if (sfp1->data.choice != sfp2->data.choice) {
    rval = FALSE;
  } else if (sfp1->data.choice == SEQFEAT_GENE) {
    grp1 = sfp1->data.value.ptrvalue;
    grp2 = sfp2->data.value.ptrvalue;
    if (StringCmp (grp1->locus_tag, grp2->locus_tag) == 0) {
      rval = TRUE;
    }
  } else if (sfp1->data.choice == SEQFEAT_CDREGION) {
    /* match protein IDs */
    gbq = sfp1->qual;
    while (gbq != NULL && StringCmp (gbq->qual, "protein_id") != 0) {
      gbq = gbq->next;
    }
    if (gbq != NULL) {
      sip = MakeSeqID (gbq->val);
      /* look for protein ID on comparison feature */
      pbsp = BioseqFindFromSeqLoc (sfp2->product);
      if (pbsp != NULL && sip != NULL && SeqIdIn (sip, pbsp->id)) {
        rval = TRUE;
      }
      sip = SeqIdFree (sip);
    }
  }
  return rval;
}


static Boolean EntityIDAlreadyInList (Uint2 entityID, ValNodePtr entityIDList)
{
  while (entityIDList != NULL) {
    if ((Uint2)(entityIDList->data.intvalue) == entityID) {
      return TRUE;
    }
  }
  return FALSE;
}


NLM_EXTERN void ActOnFeatureReplaceList (ValNodePtr list)
{
  ValNodePtr vnp, vnp2;
  FeatureReplacePtr fr;
  BioseqPtr         bsp, product_bsp;
  SeqAnnotPtr       sap;
  SeqFeatPtr        sfp;
  ValNodePtr        entityIDList = NULL;
  SeqEntryPtr       sep;
  Int2              genCode;
  OMProcControl     ompc;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    fr = (FeatureReplacePtr) vnp->data.ptrvalue;
    if (fr != NULL && fr->chosen && fr->sfp != NULL) {
      bsp = BioseqFindFromSeqLoc (fr->sfp->location);
      if (bsp != NULL) {
        sap = SeqAnnotNew ();
        sap->type = 1;        
        sap->data = AsnIoMemCopy (fr->sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
        MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
        ompc.input_entityID = bsp->idx.entityID;
        ompc.input_itemID = GetItemIDGivenPointer (bsp->idx.entityID, OBJ_BIOSEQ, (Pointer) bsp);
        ompc.input_itemtype = OBJ_BIOSEQ;
        ompc.output_itemtype = OBJ_SEQANNOT;
        ompc.output_data = (Pointer) sap;
        if (! AttachDataForProc (&ompc, FALSE)) {
          Message (MSG_POSTERR, "Error attaching SeqAnnot");
        } else {
          for (vnp2 = fr->replace_list; vnp2 != NULL; vnp2 = vnp2->next) {
            sfp = vnp2->data.ptrvalue;
            if (sfp != NULL) {
              sfp->idx.deleteme = TRUE;
              if (sfp->product != NULL) {
                product_bsp = BioseqFindFromSeqLoc (sfp->product);
                if (product_bsp != NULL) {
                  product_bsp->idx.deleteme = TRUE;
                }
              }
            }
          }
          if (!EntityIDAlreadyInList (bsp->idx.entityID, entityIDList)) {
            ValNodeAddInt (&entityIDList, 0, bsp->idx.entityID);
          }
          sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
          genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
          SetEmptyGeneticCodes (sap, genCode);
          PromoteXrefs (sap->data, bsp, bsp->idx.entityID);
        }
      }
    }
  }
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    DeleteMarkedObjects (vnp->data.intvalue, 0, NULL);
  }
}


typedef struct featurereplacelistdlg {
  DIALOG_MESSAGE_BLOCK
  DoC doc;
  ValNodePtr list;
  Int2       selected;
  Int2       clicked;
  Boolean    dblClick;

  Nlm_ParData ParFmt;
  Nlm_ColData NewFeatColFmt[2];
  Nlm_ColData OldFeatColFmt[2];

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} FeatureReplaceListDlgData, PNTR FeatureReplaceListDlgPtr;


static void InitColAndParForFeatureReplaceList (DialoG d)
{
  FeatureReplaceListDlgPtr dlg;
  RecT r;
  Int4 doc_width;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  dlg->ParFmt.openSpace = FALSE;
  dlg->ParFmt.keepWithNext = FALSE;
  dlg->ParFmt.keepTogether = FALSE;
  dlg->ParFmt.newPage = FALSE;
  dlg->ParFmt.tabStops = FALSE;
  dlg->ParFmt.minLines = 0;
  dlg->ParFmt.minHeight = 0;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  doc_width = r.right - r.left;

  dlg->NewFeatColFmt[0].pixWidth = 16;
  dlg->NewFeatColFmt[0].pixInset = 0;
  dlg->NewFeatColFmt[0].charWidth = 0;
  dlg->NewFeatColFmt[0].charInset = 0;
  dlg->NewFeatColFmt[0].font = programFont;
  dlg->NewFeatColFmt[0].just = 'l';
  dlg->NewFeatColFmt[0].wrap = 0;
  dlg->NewFeatColFmt[0].bar = 0;
  dlg->NewFeatColFmt[0].underline = 0;
  dlg->NewFeatColFmt[0].left = 0;
  dlg->NewFeatColFmt[0].last = 0;

  dlg->NewFeatColFmt[1].pixWidth = doc_width - 16;
  dlg->NewFeatColFmt[1].pixInset = 0;
  dlg->NewFeatColFmt[1].charWidth = 0;
  dlg->NewFeatColFmt[1].charInset = 0;
  dlg->NewFeatColFmt[1].font = programFont;
  dlg->NewFeatColFmt[1].just = 'l';
  dlg->NewFeatColFmt[1].wrap = 1;
  dlg->NewFeatColFmt[1].bar = 0;
  dlg->NewFeatColFmt[1].underline = 0;
  dlg->NewFeatColFmt[1].left = 0;
  dlg->NewFeatColFmt[1].last = 1;

  
  dlg->OldFeatColFmt[0].pixWidth = 32;
  dlg->OldFeatColFmt[0].pixInset = 0;
  dlg->OldFeatColFmt[0].charWidth = 0;
  dlg->OldFeatColFmt[0].charInset = 0;
  dlg->OldFeatColFmt[0].font = programFont;
  dlg->OldFeatColFmt[0].just = 'l';
  dlg->OldFeatColFmt[0].wrap = 0;
  dlg->OldFeatColFmt[0].bar = 0;
  dlg->OldFeatColFmt[0].underline = 0;
  dlg->OldFeatColFmt[0].left = 0;
  dlg->OldFeatColFmt[0].last = 0;

  dlg->OldFeatColFmt[1].pixWidth = doc_width - 32;
  dlg->OldFeatColFmt[1].pixInset = 0;
  dlg->OldFeatColFmt[1].charWidth = 0;
  dlg->OldFeatColFmt[1].charInset = 0;
  dlg->OldFeatColFmt[1].font = programFont;
  dlg->OldFeatColFmt[1].just = 'l';
  dlg->OldFeatColFmt[1].wrap = 1;
  dlg->OldFeatColFmt[1].bar = 0;
  dlg->OldFeatColFmt[1].underline = 0;
  dlg->OldFeatColFmt[1].left = 0;
  dlg->OldFeatColFmt[1].last = 1;
}


/* calculates position of feature replace item and feature within replace item */
static Int2 GetListPosFromItem (ValNodePtr list, Int2 item, Int2Ptr feature_pos)
{
  Int2 num = 1, pos = 1, nextval;
  FeatureReplacePtr fr;
  ValNodePtr vnp;

  if (feature_pos != NULL) {
    *feature_pos = 0;
  }
  if (list == NULL || item < 1) {
    return -1;
  }

  vnp = list;

  while (num < item && vnp != NULL) {
    fr = (FeatureReplacePtr) vnp->data.ptrvalue;
    nextval = num + ValNodeLen (fr->replace_list) + 1;
    if (nextval < item ) {
      pos++;
      num = nextval;
    } else if (nextval == item) {
      pos++;
      if (feature_pos != NULL) {
        *feature_pos = 0;
      }
      return pos;
    } else {
      if (feature_pos != NULL) {
        *feature_pos = nextval - item;
      }
      return pos;
    }
    vnp = vnp->next;
  }
  return pos; 
}


static void DrawFeatureReplaceList (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  FeatureReplaceListDlgPtr dlg;
  RecT                     rct;
  FeatureReplacePtr        fr;
  ValNodePtr               vnp;
  Int2                     pos, num, feature_pos = 0;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    vnp = dlg->list;
    pos = GetListPosFromItem (dlg->list, item, &feature_pos);
    if (pos < 1) {
      /* do nothing */
    } else if (feature_pos > 0) {
      /* do nothing */
    } else {
      /* draw box for selecting feature to be imported */
      num = 1;
      while (num < pos) {
        vnp = vnp->next;
        num++;
      }
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        fr = (FeatureReplacePtr) vnp->data.ptrvalue;

        /* draw selection */
        if (item == dlg->selected) {
          rct.right = rct.left + 4;
          PaintRect (&rct);
        }

        /* draw chosen checkboxes */
        rct.left += 5;
        rct.right = rct.left + 10;
        rct.bottom = rct.top + (rct.right - rct.left);
        FrameRect (&rct);
      
        if (fr->chosen) {
          MoveTo (rct.left, rct.top);
          LineTo (rct.right - 1, rct.bottom - 1);
          MoveTo (rct.left, rct.bottom - 1);
          LineTo (rct.right - 1, rct.top);
        }
      }
    }
  }
}


static Boolean HighlightFeatureReplaceTarget (DoC doc, Int2 item, Int2 row, Int2 col)
{
  FeatureReplaceListDlgPtr dlg;
  
  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
    
  if (dlg->selected == item) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}

static CharPtr GetNewItemDescription (SeqFeatPtr sfp)
{
  CharPtr location, label, row_text;
  Char buf[129];
  
  if (sfp == NULL) {
    return StringSave ("Misc Features to Delete");
  }
  location = SeqLocPrintUseBestID (sfp->location);
  label = (CharPtr) FeatDefTypeLabel(sfp);

  FeatDefLabel (sfp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);

  row_text = (CharPtr) MemNew (sizeof (Char) * 
                              (StringLen (label) 
                              + StringLen (buf) 
                              + StringLen (location) 
                              + 6));
  sprintf (row_text, "%s:%s:%s\n", label, buf, location);
  location = MemFree (location);
  return row_text;
}


static void AddFeatureReplace (FeatureReplaceListDlgPtr dlg, FeatureReplacePtr fr)
{
  CharPtr            desc, item_text, desc_fmt = "\t%s", repl_fmt = "\tReplaces %s";
  ValNodePtr         vnp;

  if (dlg == NULL || fr == NULL)
  {
    return;
  }
  desc = GetNewItemDescription (fr->sfp);
  item_text = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + StringLen (desc)));
  sprintf (item_text, desc_fmt, desc);
  desc = MemFree (desc);
  AppendText (dlg->doc, item_text, &dlg->ParFmt, dlg->NewFeatColFmt, programFont);
  MemFree (item_text);
  for (vnp = fr->replace_list; vnp != NULL; vnp = vnp->next) {
    desc = GetNewItemDescription (vnp->data.ptrvalue);
    item_text = (CharPtr) MemNew (sizeof (Char) * (StringLen (repl_fmt) + StringLen (desc)));
    sprintf (item_text, repl_fmt, desc);
    desc = MemFree (desc);
    AppendText (dlg->doc, item_text, &dlg->ParFmt, dlg->OldFeatColFmt, programFont);
    MemFree (item_text);
  }
}


static void ClickFeatureReplaceList (DoC d, PoinT pt)

{
  Int2             item, numItems, pos, feature_pos;
  Int2             row;
  Int2             col;
  FeatureReplaceListDlgPtr dlg;
  Int4             offset;
  Int2             first_shown;
  RecT             r;
  BaR              vbar;
  Int4             scroll_pos;

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && dlg->clicked == item) {
      dlg->dblClick = dblClick;
    } else {
      dlg->dblClick = FALSE;
    }
    dlg->clicked = 0;
    if (item > 0 && row > 0) {
      dlg->clicked = item;
    }
    if (item > 0 && row > 0 && !dblClick)
    {
      vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
      scroll_pos = GetBarValue (vbar);

      pos = GetListPosFromItem (dlg->list, item, &feature_pos);
      if (pos >= 1 && col == 2 && feature_pos == 0) {
        dlg->selected = item;
      }
      GetDocParams (d, &numItems, NULL);
      UpdateDocument (d, 0, numItems);
      GetItemParams4 (dlg->doc, first_shown, &offset, NULL, NULL, NULL, NULL);
      SetScrlParams4 (dlg->doc, offset);
      ObjectRect (dlg->doc, &r);
      InsetRect (&r, -1, -1);
      InvalRect (&r);    
      SetBarValue (vbar, scroll_pos);
      if (dlg->change_notify != NULL) {
        (dlg->change_notify) (dlg->change_userdata);
      }
    }
  }
}


static void PopulateFeatureReplaceList (FeatureReplaceListDlgPtr dlg, ValNodePtr list)
{
  Int2               numItems;
  
  if (dlg == NULL || dlg->doc == NULL) 
  {
    return;
  }
  
  Reset (dlg->doc);
  
  
  while (list != NULL)
  {
    AddFeatureReplace (dlg, list->data.ptrvalue);
    list = list->next;
  }
  GetDocParams (dlg->doc, &numItems, NULL);
  UpdateDocument (dlg->doc, 0, numItems);

}


static void FeatureReplaceListToDialog (DialoG d, Pointer data)
{
  FeatureReplaceListDlgPtr dlg;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  dlg->list = FeatureReplaceListFree (dlg->list);
  dlg->list = FeatureReplaceListCopy((ValNodePtr) data);
  PopulateFeatureReplaceList (dlg, dlg->list);
}


static Pointer FeatureReplaceListFromDialog (DialoG d)
{
  FeatureReplaceListDlgPtr dlg;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  return FeatureReplaceListCopy(dlg->list);
}


static void CleanupFeatureReplaceListDialog (GraphiC g, VoidPtr data)

{
  FeatureReplaceListDlgPtr dlg;

  dlg = (FeatureReplaceListDlgPtr) data;
  if (dlg != NULL) {
    dlg->list = FeatureReplaceListFree (dlg->list);
  }
  StdCleanupExtraProc (g, data);
}


NLM_EXTERN DialoG 
FeatureReplaceListDialog 
(GrouP h, 
 Int4 width,
 Nlm_ChangeNotifyProc change_notify,
 Pointer change_userdata)
{
  FeatureReplaceListDlgPtr dlg;
  GrouP                    p;

  dlg = (FeatureReplaceListDlgPtr) MemNew (sizeof (FeatureReplaceListDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupFeatureReplaceListDialog);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = FeatureReplaceListFromDialog;
  dlg->todialog = FeatureReplaceListToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->doc = DocumentPanel (p, width, stdLineHeight * 20);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocProcs (dlg->doc, ClickFeatureReplaceList, NULL, NULL, NULL);
  SetDocShade (dlg->doc, DrawFeatureReplaceList, NULL, HighlightFeatureReplaceTarget, NULL);

  InitColAndParForFeatureReplaceList (dlg->dialog);

  return (DialoG) p;
}


NLM_EXTERN Boolean AddFeaturesToReplaceList (DialoG d, ValNodePtr feature_list)
{
  FeatureReplaceListDlgPtr dlg;
  Int2                     num, pos;
  ValNodePtr               vnp = NULL;
  FeatureReplacePtr        fr;
  Boolean                  rval = FALSE;
  BaR                      vbar;
  Int4                     scroll_pos;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || feature_list == NULL) {
    return FALSE;
  }

  vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
  scroll_pos = GetBarValue (vbar);

  if (dlg->selected >= 1) {
    pos = GetListPosFromItem (dlg->list, dlg->selected, NULL);
    num = 1;
    vnp = dlg->list;
    while (num < pos && vnp != NULL) {
      num++;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    fr = (FeatureReplacePtr) vnp->data.ptrvalue;
    if (fr != NULL) {
      ValNodeLink (&fr->replace_list, feature_list);
      PopulateFeatureReplaceList (dlg, dlg->list);
      rval = TRUE;
    }
  }
  if (!rval) {
    /* find misc item */
    vnp = dlg->list;
    while (vnp != NULL && ((fr = (FeatureReplacePtr) vnp->data.ptrvalue) == NULL
                           || fr->sfp != NULL)) {
      vnp = vnp->next;
    }
    if (vnp == NULL) {
      fr = FeatureReplaceNew (NULL);
      ValNodeAddPointer (&dlg->list, 0, fr);
    }
    ValNodeLink (&(fr->replace_list), feature_list);
    PopulateFeatureReplaceList (dlg, dlg->list);
    rval = TRUE;
  }

  SetBarValue (vbar, scroll_pos);

  return rval;
}


NLM_EXTERN SeqFeatPtr GetSelectedNewFeature (DialoG d)
{
  FeatureReplaceListDlgPtr dlg;
  FeatureReplacePtr fr;
  SeqFeatPtr sfp = NULL;
  ValNodePtr vnp;
  Int2       num, pos;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  if (dlg->selected >= 1) {
    pos = GetListPosFromItem (dlg->list, dlg->selected, NULL);
    num = 1;
    vnp = dlg->list;
    while (num < pos && vnp != NULL) {
      num++;
      vnp = vnp->next;
    }
    if (vnp != NULL) {
      fr = (FeatureReplacePtr) vnp->data.ptrvalue;
      if (fr != NULL) {
        sfp = fr->sfp;
      }
    }
  }
  return sfp;
}


NLM_EXTERN ValNodePtr RemoveFeaturesFromReplaceList (DialoG d)
{
  FeatureReplaceListDlgPtr dlg;
  Int2                     num;
  ValNodePtr               vnp, list = NULL, prev = NULL;
  FeatureReplacePtr        fr;
  BaR                      vbar;
  Int4                     scroll_pos;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->list == NULL || dlg->selected < 1) {
    return NULL;
  }

  num = 1;
  vnp = dlg->list;
  while (num < dlg->selected && vnp != NULL) {
    num++;
    prev = vnp;
    vnp = vnp->next;
  }
  if (vnp != NULL) {
    fr = (FeatureReplacePtr) vnp->data.ptrvalue;
    if (fr != NULL) {
      list = fr->replace_list;
      fr->replace_list = NULL;
      /* if we have removed features from a placeholder empty category, remove the category */
      if (fr->sfp == NULL) {
        if (prev == NULL) {
          dlg->list->next = vnp->next;
        } else {
          prev->next = vnp->next;
        }
        vnp->next = NULL;
        vnp = FeatureReplaceListFree (vnp);
      }
      vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
      scroll_pos = GetBarValue (vbar);

      PopulateFeatureReplaceList (dlg, dlg->list);
      SetBarValue (vbar, scroll_pos);
    }
  }
  return list;
}


NLM_EXTERN Boolean AutomatchFeatures (DialoG d, ValNodePtr PNTR existing_features)
{
  FeatureReplaceListDlgPtr dlg;
  ValNodePtr               vnp, vnp_f, vnp_next, list = NULL, prev = NULL;
  FeatureReplacePtr        fr;
  Boolean                  found = FALSE;

  dlg = (FeatureReplaceListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->list == NULL || existing_features == NULL || *existing_features == NULL) {
    return FALSE;
  }

  for (vnp = dlg->list; vnp != NULL; vnp = vnp->next) {
    fr = (FeatureReplacePtr) vnp->data.ptrvalue;
    prev = NULL;
    for (vnp_f = *existing_features; vnp_f != NULL; vnp_f = vnp_next) {
      vnp_next = vnp_f->next;
      if (DoFeaturesAutomatchForUpdate (fr->sfp, vnp_f->data.ptrvalue)) {
        if (prev == NULL) {
          *existing_features = vnp_f->next;
        } else {
          prev->next = vnp_f->next;
        }
        vnp_f->next = NULL;
        ValNodeLink (&(fr->replace_list), vnp_f);
        found = TRUE;
      } else {
        prev = vnp_f;
      }
    }
  }

  PopulateFeatureReplaceList (dlg, dlg->list);

  return found;
}



typedef struct featureselectlistdlg {
  DIALOG_MESSAGE_BLOCK
  DoC doc;
  ValNodePtr list;
  Int2       selected;
  Int2       clicked;
  Boolean    dblClick;

  Nlm_ParData ParFmt;
  Nlm_ColData ColFmt[2];
} FeatureSelectListDlgData, PNTR FeatureSelectListDlgPtr;


static void ClickFeatureSelectList (DoC d, PoinT pt)

{
  Int2             item, numItems;
  Int2             first_shown;
  Int2             row;
  Int2             col;
  Int4             offset;
  FeatureSelectListDlgPtr dlg;
  RecT             r;
  BaR                     vbar;
  Int4                    scroll_pos;
  

  dlg = GetObjectExtra (d);
  if (dlg != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && dlg->clicked == item) {
      dlg->dblClick = dblClick;
    } else {
      dlg->dblClick = FALSE;
    }
    dlg->clicked = 0;
    if (item > 0 && row > 0) {
      dlg->clicked = item;
      dlg->selected = item;
    }
    vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
    scroll_pos = GetBarValue (vbar);

    GetDocParams (d, &numItems, NULL);
    UpdateDocument (d, 0, numItems);
    GetItemParams4 (dlg->doc, first_shown, &offset, NULL, NULL, NULL, NULL);
    SetScrlParams4 (dlg->doc, offset);
    ObjectRect (dlg->doc, &r);
    InsetRect (&r, -1, -1);
    InvalRect (&r);
    SetBarValue (vbar, scroll_pos);
  }
}


static Boolean HighlightFeatureSelectTarget (DoC doc, Int2 item, Int2 row, Int2 col)
{
  FeatureSelectListDlgPtr dlg;
  
  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
    
  if (dlg->selected == item) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}



static void PopulateFeatureSelectList (FeatureSelectListDlgPtr dlg, ValNodePtr list)
{
  ValNodePtr              vnp;
  CharPtr                 desc;
  Int2                    numItems;
  
  Reset (dlg->doc);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    desc = GetNewItemDescription (vnp->data.ptrvalue);
    AppendText (dlg->doc, desc, &dlg->ParFmt, dlg->ColFmt, programFont);
    desc = MemFree (desc);
  }

  GetDocParams (dlg->doc, &numItems, NULL);
  UpdateDocument (dlg->doc, 0, numItems);   
}


static void FeatureSelectListToDialog (DialoG d, Pointer data)
{
  FeatureSelectListDlgPtr dlg;

  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  dlg->list = (ValNodePtr) data;
  PopulateFeatureSelectList (dlg, dlg->list);
}


static Pointer FeatureSelectListFromDialog (DialoG d)
{
  FeatureSelectListDlgPtr dlg;

  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return (Pointer) dlg->list;
}


NLM_EXTERN DialoG FeatureSelectListDialog (GrouP h, Int4 width)
{
  FeatureSelectListDlgPtr dlg;
  GrouP                   p;
  RecT                    r;
  Int4                    doc_width;

  dlg = (FeatureSelectListDlgPtr) MemNew (sizeof (FeatureSelectListDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatureSelectListToDialog;
  dlg->fromdialog = FeatureSelectListFromDialog;

  dlg->doc = DocumentPanel (p, width, stdLineHeight * 20);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocShade (dlg->doc, NULL, NULL, HighlightFeatureSelectTarget, NULL);
  SetDocProcs (dlg->doc, ClickFeatureSelectList, NULL, NULL, NULL);

  dlg->ParFmt.openSpace = FALSE;
  dlg->ParFmt.keepWithNext = FALSE;
  dlg->ParFmt.keepTogether = FALSE;
  dlg->ParFmt.newPage = FALSE;
  dlg->ParFmt.tabStops = FALSE;
  dlg->ParFmt.minLines = 0;
  dlg->ParFmt.minHeight = 0;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  doc_width = r.right - r.left;

  dlg->ColFmt[0].pixWidth = doc_width;
  dlg->ColFmt[0].pixInset = 0;
  dlg->ColFmt[0].charWidth = 0;
  dlg->ColFmt[0].charInset = 0;
  dlg->ColFmt[0].font = programFont;
  dlg->ColFmt[0].just = 'l';
  dlg->ColFmt[0].wrap = 0;
  dlg->ColFmt[0].bar = 0;
  dlg->ColFmt[0].underline = 0;
  dlg->ColFmt[0].left = 0;
  dlg->ColFmt[0].last = 1;


  return (DialoG) p;
}


NLM_EXTERN ValNodePtr RemoveSelectedFeaturesFromList (DialoG d)
{
  FeatureSelectListDlgPtr dlg;
  ValNodePtr              vnp, prev = NULL;
  Int2                    num;
  BaR                     vbar;
  Int4                    scroll_pos;

  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->selected == 0) {
    return NULL;
  }

  vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
  scroll_pos = GetBarValue (vbar);

  num = 1;
  vnp = dlg->list;
  while (num < dlg->selected && vnp != NULL) {
    prev = vnp;
    vnp = vnp->next;
    num++;
  }

  if (vnp != NULL) {
    if (prev == NULL) {
      dlg->list = vnp->next;
    } else {
      prev->next = vnp->next;
    }
    vnp->next = NULL;
  }
  PopulateFeatureSelectList (dlg, dlg->list);
  SetBarValue (vbar, scroll_pos);
  return vnp;
}


NLM_EXTERN void ScrollToMatchingFeatures (DialoG d, SeqFeatPtr sfp)
{
  FeatureSelectListDlgPtr dlg;
  ValNodePtr              vnp;
  BaR                     vbar;
  Int4                    scroll_pos;
  Int4                    start, stop, startitem;
  Int4                    range_start = -1, pos;
  RecT                    r;

  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || sfp == NULL) {
    return;
  }

  dlg->selected = 0;
  start = SeqLocStart (sfp->location);
  stop = SeqLocStop (sfp->location);

  for (vnp = dlg->list, pos = 0; vnp != NULL; vnp = vnp->next, pos++) {
    startitem = SeqLocStart (((SeqFeatPtr)vnp->data.ptrvalue)->location);
    if (startitem >= start && range_start == -1) {
      range_start = pos;
    }
    if (DoFeaturesAutomatchForUpdate (sfp, vnp->data.ptrvalue)) {
      range_start = pos;
      dlg->selected = pos + 1;
      break;
    }
    if (startitem > stop) {
      break;
    }
  }

  /* scroll to position */
  if (range_start > -1) {
    vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
    scroll_pos = GetBarMax (vbar);
    if (scroll_pos > range_start) {
      scroll_pos = range_start;
    }
    SetBarValue (vbar, scroll_pos);
  }
  ObjectRect (dlg->doc, &r);
  InsetRect (&r, -1, -1);
  InvalRect (&r);
  Update ();

}


typedef struct featsort {
  SeqFeatPtr sfp;
  Int4       start;
  Int4       stop;
} FeatSortData, PNTR FeatSortPtr;


static FeatSortPtr FeatSortNew (SeqFeatPtr sfp)
{
  FeatSortPtr f;

  if (sfp == NULL) {
    return NULL;
  }

  f = (FeatSortPtr) MemNew (sizeof (FeatSortData));
  f->sfp = sfp;
  f->start = SeqLocStart (sfp->location);
  f->stop = SeqLocStop (sfp->location);
  return f;
}


static FeatSortPtr FeatSortFree (FeatSortPtr f)
{
  if (f != NULL) {
    f = MemFree (f);
  }
  return f;
}


static ValNodePtr FeatSortListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = FeatSortFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static ValNodePtr FeatSortListFromSeqAnnot (SeqAnnotPtr sap)
{
  ValNodePtr list = NULL, prev = NULL, vnp;
  SeqFeatPtr sfp;
  FeatSortPtr f;

  if (sap == NULL || sap->type != 1) {
    return NULL;
  }

  for (sfp = sap->data; sfp != NULL; sfp = sfp->next) {
    f = FeatSortNew (sfp);
    vnp = ValNodeNew (NULL);
    vnp->data.ptrvalue = f;
    if (prev == NULL) {
      list = vnp;
    } else {
      prev->next = vnp;
    }
    prev = vnp;
  }
  return list;
}


static SeqFeatPtr SeqFeatListFromFeatSortList (ValNodePtr list)
{
  SeqFeatPtr first = NULL, last = NULL;
  FeatSortPtr f;

  while (list != NULL) {
    f = (FeatSortPtr) list->data.ptrvalue;
    if (f != NULL) {
      if (first == NULL) {
        first = f->sfp;
      } else {
        last->next = f->sfp;
      }
      last = f->sfp;
      last->next = NULL;
    }
    list = list->next;
  }
  return first;
}


static ValNodePtr FeatSortListFromFeatValNodeList (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp, vnp_orig;
  FeatSortPtr f;
  
  for (vnp_orig = orig; vnp_orig != NULL; vnp_orig = vnp_orig->next) {
    f = FeatSortNew ((SeqFeatPtr) vnp_orig->data.ptrvalue);
    vnp = ValNodeNew (NULL);
    vnp->data.ptrvalue = f;
    if (prev == NULL) {
      list = vnp;
    } else {
      prev->next = vnp;
    }
    prev = vnp;
  }
  return list;
}


static ValNodePtr FeatValNodeListFromFeatSortList (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp, vnp_orig;
  FeatSortPtr f;
  
  for (vnp_orig = orig; vnp_orig != NULL; vnp_orig = vnp_orig->next) {
    f = (FeatSortPtr) vnp_orig->data.ptrvalue;
    vnp = ValNodeNew (NULL);
    vnp->data.ptrvalue = f->sfp;
    if (prev == NULL) {
      list = vnp;
    } else {
      prev->next = vnp;
    }
    prev = vnp;
  }
  return list;
}

static int LIBCALLBACK SortVnpByFeatSort (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  FeatSortPtr f1, f2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      f1 = (FeatSortPtr) vnp1->data.ptrvalue;
      f2 = (FeatSortPtr) vnp2->data.ptrvalue;
      if (f1 != NULL && f2 != NULL) {
        if (f1->start < f2->start) {
          return -1;
        } else if (f1->start > f2->start) {
          return 1;
        } else if (f1->stop > f2->stop) {
          return -1;
        } else if (f1->stop < f2->stop) {
          return 1;
        } else if (f1->sfp->data.choice < f2->sfp->data.choice) {
          return -1;
        } else if (f1->sfp->data.choice > f2->sfp->data.choice) {
          return 1;
        } else {
          return 0;
        }
      }
    }
  }
  return 0;
}



NLM_EXTERN void SortSeqFeatInAnnot (SeqAnnotPtr sap)
{
  ValNodePtr list;

  list = FeatSortListFromSeqAnnot (sap);
  list = ValNodeSort (list, SortVnpByFeatSort);
  if (list != NULL) {
    sap->data = SeqFeatListFromFeatSortList (list);
    list = FeatSortListFree (list);
  }
}




NLM_EXTERN void AddFeaturesToList (DialoG d, ValNodePtr features)
{
  FeatureSelectListDlgPtr dlg;
  BaR                     vbar;
  Int4                    scroll_pos;
  ValNodePtr              tmp;

  dlg = (FeatureSelectListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || features == NULL) {
    return;
  }

  ValNodeLink (&dlg->list, features);
  tmp = FeatSortListFromFeatValNodeList (dlg->list);
  tmp = ValNodeSort (tmp, SortVnpByFeatSort);
  dlg->list = ValNodeFree (dlg->list);
  dlg->list = FeatValNodeListFromFeatSortList (tmp);
  tmp = FeatSortListFree (tmp);

  vbar = GetSlateVScrollBar ((SlatE) dlg->doc);
  scroll_pos = GetBarValue (vbar);
  PopulateFeatureSelectList (dlg, dlg->list);
  SetBarValue (vbar, scroll_pos);
}






