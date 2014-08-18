/*   mla2api.c
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
* File Name:  mla2api.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/30/07
*
* $Revision: 1.15 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#include <mla2api.h>
#include <urlquery.h>
#include <ncbithr.h>

/* low-level connection functions */

NLM_EXTERN CONN Mla2OpenConnection (void)

{
  return QUERY_OpenServiceQuery ("MedArch", NULL, 30);
}

#ifdef OS_MAC
#include <Events.h>
#endif

NLM_EXTERN MlaBackPtr Mla2WaitForReply (
  CONN conn
)

{
  AsnIoConnPtr  aicp;
  time_t        currtime, starttime;
  time_t        max = 0;
  MlaBackPtr    mbp = NULL;
  EIO_Status    status;
  STimeout      timeout;
#ifdef OS_MAC
  EventRecord   currEvent;
#endif

  if (conn == NULL) return NULL;

#ifdef OS_MAC
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 100;
  timeout.usec = 0;
#endif

  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) == eIO_Timeout && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
#ifdef OS_MAC
    WaitNextEvent (0, &currEvent, 0, NULL);
#endif
  }
  if (status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen ("rb", conn);
    mbp = MlaBackAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  CONN_Close (conn);

  return mbp;
}

/* high-level connection functions */

NLM_EXTERN MlaBackPtr Mla2SynchronousQuery (
  MlaRequestPtr mrp
)

{
  AsnIoConnPtr  aicp;
  CONN          conn;
  MlaBackPtr    mbp = NULL;

  if (mrp == NULL) return NULL;

  conn = Mla2OpenConnection ();

  if (conn == NULL) return NULL;

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  MlaRequestAsnWrite (mrp, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  mbp = Mla2WaitForReply (conn);

  return mbp;
}

NLM_EXTERN Boolean Mla2AsynchronousQuery (
  MlaRequestPtr mrp,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
)

{
  AsnIoConnPtr  aicp;
  CONN          conn;

  if (mrp == NULL) return FALSE;

  conn = Mla2OpenConnection ();

  if (conn == NULL) return FALSE;

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  Mla2RequestAsnWrite (mrp, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

NLM_EXTERN Int4 Mla2CheckQueue (
  QUEUE* queue
)

{
  return QUERY_CheckQueue (queue);
}

NLM_EXTERN MlaBackPtr Mla2ReadReply (
  CONN conn,
  EIO_Status status
)

{
  AsnIoConnPtr  aicp;
  MlaBackPtr    mbp = NULL;

  if (conn != NULL && status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen ("rb", conn);
    mbp = Mla2BackAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }

  return mbp;
}

/* request creation functions */

NLM_EXTERN MlaRequestPtr Mla2CreateJournalTitleRequest (
  CharPtr journal_name
)

{
  ValNodePtr     jta;
  MlaRequestPtr  mrp;
  TitleMsgPtr    tmp;

  if (StringHasNoText (journal_name)) return NULL;

  jta = ValNodeNew (NULL);
  jta->choice = Title_type_name;
  jta->data.ptrvalue = StringSave (journal_name);

  tmp = TitleMsgNew ();
  tmp->type = Title_type_all;
  tmp->title = jta;

  mrp = ValNodeNew (NULL);
  mrp->choice = MlaRequest_gettitle;
  mrp->data.ptrvalue = (Pointer) tmp;

  return mrp;
}

NLM_EXTERN MlaRequestPtr Mla2CreateCitArtJournalRequest (
  CitArtPtr cap
)

{
  CitJourPtr  cjp;
  CharPtr     str;
  ValNodePtr  ttl;

  if (cap == NULL) return NULL;
  if (cap->from != 1) return NULL;
  cjp = (CitJourPtr) cap->fromptr;
  if (cjp == NULL) return NULL;

  for (ttl = cjp->title; ttl != NULL; ttl = ttl->next) {
    if (ttl->choice == Cit_title_name ||
        ttl->choice == Cit_title_trans ||
        ttl->choice == Cit_title_jta ||
        ttl->choice == Cit_title_iso_jta ||
        ttl->choice == Cit_title_ml_jta ||
        ttl->choice == Cit_title_coden ||
        ttl->choice == Cit_title_issn) {
      str = (CharPtr) ttl->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      return Mla2CreateJournalTitleRequest (str);
    }
  }

  return NULL;
}

NLM_EXTERN MlaRequestPtr Mla2CreateCitArtMatchRequest (
  CitArtPtr cap
)

{
  MlaRequestPtr  mrp;
  ValNodePtr     pub;


  if (cap == NULL) return NULL;
  cap = AsnIoMemCopy ((Pointer) cap,
                      (AsnReadFunc) CitArtAsnRead,
                      (AsnWriteFunc) CitArtAsnWrite);

  pub = ValNodeNew (NULL);
  pub->choice = PUB_Article;
  pub->data.ptrvalue = cap;

  mrp = ValNodeNew (NULL);
  mrp->choice = MlaRequest_citmatch;
  mrp->data.ptrvalue = (Pointer) pub;

  return mrp;
}

NLM_EXTERN MlaRequestPtr Mla2CreateCitationtMatchRequest (
  CharPtr author,
  CharPtr journal,
  CharPtr volume,
  CharPtr page,
  Int4 year,
  CharPtr title
)

{
  AuthListPtr    alp;
  AuthorPtr      ap;
  CitArtPtr      cap;
  CitJourPtr     cjp;
  DatePtr        dp;
  ImprintPtr     imp;
  ValNodePtr     jta;
  MlaRequestPtr  mrp;
  ValNodePtr     names;
  NameStdPtr     nsp;
  PersonIdPtr    pid;
  ValNodePtr     pub;
  ValNodePtr     ttl;

  if (StringHasNoText (author) || StringHasNoText (journal) || year < 1900) return NULL;

  cap = CitArtNew ();
  if (cap == NULL) return NULL;

  nsp = NameStdNew ();
  nsp->names [0] = StringSave (author);

  pid = PersonIdNew ();
  pid->choice = 2;
  pid->data = nsp;

  ap = AuthorNew ();
  ap->name = pid;

  names = ValNodeNew (NULL);
  names->choice = 1;
  names->data.ptrvalue = ap;

  alp = AuthListNew ();
  alp->choice = 1;
  alp->names = names;

  cap->authors = alp;

  jta = ValNodeNew (NULL);
  jta->choice = Title_type_name;
  jta->data.ptrvalue = StringSave (journal);

  cjp = CitJourNew ();
  cjp->title = jta;

  dp = DateNew ();
  dp->data [0] = 1;
  dp->data [1] = year - 1900;

  imp = ImprintNew ();
  imp->date = dp;
  imp->volume = StringSaveNoNull (volume);
  imp->pages = StringSaveNoNull (page);

  cjp->imp = imp;

  cap->from = 1;
  cap->fromptr = cjp;

  if (StringDoesHaveText (title)) 
  {
    ttl = ValNodeNew (NULL);
    ttl->choice = Title_type_name;
    ttl->data.ptrvalue = StringSave (title);

    cap->title = ttl;
  }

  pub = ValNodeNew (NULL);
  pub->choice = PUB_Article;
  pub->data.ptrvalue = cap;

  mrp = ValNodeNew (NULL);
  mrp->choice = MlaRequest_citmatch;
  mrp->data.ptrvalue = (Pointer) pub;

  return mrp;
}

NLM_EXTERN MlaRequestPtr Mla2CreatePubFetchRequest (
  Int4 pmid
)

{
  MlaRequestPtr  mrp;

  if (pmid  < 1) return NULL;

  mrp = ValNodeNew (NULL);
  mrp->choice = MlaRequest_getpub;
  mrp->data.intvalue = pmid;

  return mrp;
}

/* reply extraction functions */

NLM_EXTERN TitleMsgListPtr Mla2ExtractJournalTitleReply (
  MlaBackPtr mbp
)

{
  TitleMsgListPtr  tlp;

  if (mbp == NULL || mbp->choice != MlaBack_gettitle) return NULL;
  tlp = (TitleMsgListPtr) mbp->data.ptrvalue;
  if (tlp == NULL) return NULL;
  mbp->data.ptrvalue = NULL;

  return tlp;
}

NLM_EXTERN Int4 Mla2ExtractCitMatchReply (
  MlaBackPtr mbp
)

{
  Int4  pmid;

  if (mbp == NULL || mbp->choice != MlaBack_citmatch) return 0;
  pmid = (Int4) mbp->data.intvalue;
  if (pmid < 1) return 0;

  return pmid;
}

NLM_EXTERN CitArtPtr Mla2ExtractPubFetchReply (
  MlaBackPtr mbp
)

{
  CitArtPtr   cap;
  ValNodePtr  pub;

  if (mbp == NULL || mbp->choice != MlaBack_getpub) return NULL;
  pub = mbp->data.ptrvalue;
  if (pub == NULL || pub->choice != PUB_Article) return NULL;
  cap = (CitArtPtr) pub->data.ptrvalue;
  if (cap == NULL) return NULL;
  pub->data.ptrvalue = NULL;

  return cap;
}

/* utility functions */

/* ml to std code modified from original in medutil.c and then in pmfapi.c */

static Boolean StrIsAllUpperCase (
  CharPtr p
)

{
  Char  ch;

  if (p == NULL) return FALSE;
  ch = *p;
  while (ch != '\0') {
    if (! IS_UPPER (ch)) return FALSE;
    p++;
    ch = *p;
  }
  return TRUE;
}

static void SplitMLAuthName (
  CharPtr name,
  CharPtr last,
  CharPtr initials,
  CharPtr suffix
)

{
  CharPtr  p, p2;
  Int2     i;
  Char     sbuf [40], ibuf [40];

  /* Clear the ibuf field and transfer the entire name to 'last',
  excluding leading and trailing spaces */

  if (name == NULL) return;

  ibuf [0] = '\0';
  sbuf [0] = '\0';
  last [0] = '\0';
  initials [0] = '\0';
  suffix [0] = '\0';
  while (*name <= ' ') {
    name++;
    if (*name == '\0') return;
  }
  StringCpy( last, name );

  for (i=StringLen (last) - 1; ((i >= 0) && (last [i] <= ' ')); i--) {
    last[i] = '\0';
  }

  /* Strip off the last token (initials or name suffix (Jr, Sr, suffix.) */

  p = StringRChr (last, (int) ' ');
  if (p != NULL) { /* more than just last name */

    /* Separate the token from the last name */

    p2 = p + 1;
    while ((p > last) && (*p == ' ')) {
      *p = '\0';
      p--;
    }

    /* If the last token is not all upper case, and there are more than
    two tokens, see if the next to the last are initials (upper case) */

    if (! StrIsAllUpperCase (p2) && (p = StringRChr (last, (int) ' ' )) != NULL) {

      /* We have at least three tokens, is the next to last initials? */

      if (StrIsAllUpperCase (p + 1)) {

        /* Yes - concatenate the last two tokens as initials */

        StringCpy (ibuf, p + 1);
        StringCpy (sbuf, p2);
        while (p > last && (*p == ' ')) {
          *p = '\0';
          p--;
        }
      }
    }
    
    if (ibuf [0] == '\0') { /* Only the last token goes in ibuf */
      StringCpy (ibuf, p2);
    }
  }

  /* now add periods to ibuf and convert suffix */

  for (p = initials, p2 = ibuf; *p2 != '\0'; p2++, p++) {
    *p = *p2;
    if (! IS_LOWER(*(p2 + 1))) { /* watch out for foreign names */
      p++;
      *p = '.';
    }
  }
  *p = '\0';

  if (sbuf [0]) {
    if (StringCmp (sbuf, "1d") == 0 || StringCmp (sbuf, "1st") == 0)
      p = StringMove (suffix, "I");
    else if (StringCmp (sbuf, "2d") == 0 || StringCmp (sbuf, "2nd") == 0)
      p = StringMove (suffix, "II");
    else if (StringCmp (sbuf, "3d") == 0 || StringCmp (sbuf, "3rd") == 0)
      p = StringMove (suffix, "III");
    else if (StringCmp (sbuf, "4th") == 0)
      p = StringMove (suffix, "IV");
    else if (StringCmp (sbuf, "5th") == 0)
      p = StringMove (suffix, "V");
    else if (StringCmp (sbuf, "6th") == 0)
      p = StringMove (suffix, "VI");
    else if (StringCmp (sbuf, "Sr") == 0)
      p = StringMove (suffix, "Sr.");
    else if (StringCmp (sbuf, "Jr") == 0)
      p = StringMove (suffix, "Jr.");
    else
      p = StringMove (suffix, sbuf);
  }
}

static ValNodePtr ChangeMLtoSTD (
  CharPtr token
)

{
  AuthorPtr   aup;
  CharPtr     eptr;
  Char        last [180], initials [40], suffix [40];
  NameStdPtr  nsp;
  PersonIdPtr pid;
  ValNodePtr  vnp;

  if (token == NULL) return NULL;
  for (eptr = token + StringLen (token) - 1;
       eptr > token && *eptr == ' ';
       eptr--) continue;

  SplitMLAuthName (token, last, initials, suffix);

  nsp = NameStdNew ();
  if (nsp == NULL) return NULL;
  nsp->names [0] = StringSave (last);
  if (initials [0] != '\0') {
    nsp->names[4] = StringSave (initials);
  }
  if (suffix[0] != '\0') {
    nsp->names[5] = StringSave (suffix);
  }
  if (nsp->names[0] != NULL) {
    pid = PersonIdNew ();
    pid->choice = 2; /* name */
    pid->data = nsp;
    aup = AuthorNew ();
    aup->name = pid;
    vnp = ValNodeNew (NULL);
    vnp->data.ptrvalue = (Pointer) aup;
    return vnp;
  }
  return NULL;
}

NLM_EXTERN void ChangeCitArtMLAuthorsToSTD (
  CitArtPtr cap
)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  ValNodePtr   curr, last = NULL, names, oldnames, tmp;
  PersonIdPtr  pid;
  CharPtr      str;

  if (cap == NULL) return;
  alp = cap->authors;
  if (alp == NULL) return;

  if (alp->choice == 1) {
    for (names = alp->names; names != NULL; names = names->next) {
      ap = names->data.ptrvalue;
      if (ap == NULL) continue;
      pid = ap->name;
      if (pid == NULL) continue;
      if (pid->choice != 3) continue;
      str = (CharPtr) pid->data;
      if (StringHasNoText (str)) continue;
      curr = ChangeMLtoSTD (str);
      if (curr == NULL) continue;
      names->data.ptrvalue = AuthorFree (ap);
      names->data.ptrvalue = curr->data.ptrvalue;
      /* do not ValNodeFreeData, since data has been switched to old author */
      ValNodeFree (curr);
    }
  }

  if (alp->choice == 2) {

    /* do not convert if too big for buffers */

    for (tmp = alp->names; tmp != NULL; tmp = tmp->next) {
      if (StringLen ((CharPtr) tmp->data.ptrvalue) > 170) return;
    }

    oldnames = alp->names;
    alp->names = NULL;
    alp->choice = 1; /* make std names */

    for (tmp = oldnames; tmp != NULL; tmp = tmp->next) {
      curr = ChangeMLtoSTD ((CharPtr) tmp->data.ptrvalue);
      if (alp->names == NULL) {
        alp->names = curr;
      }
      if (last != NULL) {
        last->next = curr;
      }
      last = curr;
    }

    ValNodeFreeData (oldnames);
  }
}

NLM_EXTERN void ChangeMlaBackMLAuthorsToSTD (
  MlaBackPtr mbp 
)

{
  CitArtPtr   cap;
  ValNodePtr  pub;

  if (mbp == NULL || mbp->choice != MlaBack_getpub) return;
  pub = mbp->data.ptrvalue;
  if (pub == NULL || pub->choice != PUB_Article) return;
  cap = (CitArtPtr) pub->data.ptrvalue;
  if (cap == NULL) return;
  ChangeCitArtMLAuthorsToSTD (cap);
}

typedef struct ejour {
  CharPtr  journal_name;
  Int2     starting_year;
} EjourData, PNTR EjourDataPtr;

static EjourData mla2_ejour_list [] = {
 {"Acta Crystallograph. Sect. F Struct. Biol. Cryst. Commun.",     0 },
 {"Acta Vet. Scand.",                                           2006 },
 {"Afr. J. Biotechnol.",                                           0 },
 {"Ambul. Surg.",                                               2005 },
 {"Ann. Clin. Microbiol. Antimicrob.",                             0 },
 {"Biol. Direct",                                                  0 },
 {"Biotecnol. Apl.",                                            2002 },
 {"BMC Biochem.",                                                  0 },
 {"BMC Bioinformatics",                                            0 },
 {"BMC Biotechnol.",                                               0 },
 {"BMC Cancer",                                                    0 },
 {"BMC Cell Biol.",                                                0 },
 {"BMC Dermatol.",                                                 0 },
 {"BMC Dev. Biol.",                                                0 },
 {"BMC Ecol.",                                                     0 },
 {"BMC Evol. Biol.",                                               0 },
 {"BMC Genet.",                                                    0 },
 {"BMC Genomics",                                                  0 },
 {"BMC Immunol.",                                                  0 },
 {"BMC Infect. Dis.",                                              0 },
 {"BMC Med. Genet.",                                               0 },
 {"BMC Microbiol.",                                                0 },
 {"BMC Mol. Biol.",                                                0 },
 {"BMC Pharmacol.",                                                0 },
 {"BMC Physiol.",                                                  0 },
 {"BMC Plant Biol.",                                               0 },
 {"BMC Struct. Biol.",                                             0 },
 {"BMC Vet. Res.",                                                 0 },
 {"Breast Cancer Res.",                                         2005 },
 {"Cancer Cell Int.",                                              0 },
 {"Cancer Immun.",                                                 0 },
 {"Cell Commun. Signal",                                           0 },
 {"Cell Struct. Funct.",                                        2005 },
 {"Cell. Mol. Biol. (Noisy-le-grand)",                          2004 },
 {"Colomb. Med.",                                               1998 },
 {"Crit. Rev. Oral Biol. Med.",                                 2002 },
 {"Dermatol. Online J.",                                           0 },
 {"Electron. J. Biotechnol.",                                      0 },
 {"Eur. J. Genet. Mol. Toxicol.",                                  0 },
 {"Evol. Bioinform. Online",                                       0 },
 {"Front. Zool.",                                                  0 },
 {"Fungal Planet",                                                 0 },
 {"Genome Biol.",                                               2005 },
 {"Geochem. Trans.",                                               0 },
 {"Hereditas",                                                  2004 },
 {"Hum. Genomics",                                              2004 },
 {"Infect. Agents Cancer",                                         0 },
 {"J. Exp. Clin. Assist. Reprod.",                                 0 },
 {"J. Plankton Res.",                                              0 },
 {"JOP",                                                           0 },
 {"Kinetoplastid Biol Dis",                                        0 },
 {"Malar. J.",                                                     0 },
 {"Microb. Cell Fact.",                                            0 },
 {"Mol. Syst. Biol.",                                              0 },
 {"Mol. Vis.",                                                     0 },
 {"Neoplasia",                                                  2005 },
 {"Nucl. Recept.",                                                 0 },
 {"Pharmacologyonline",                                            0 },
 {"Plant Methods",                                                 0 },
 {"PLoS Biol.",                                                 2006 },
 {"PLoS Comput. Biol.",                                            0 },
 {"PLoS Genet.",                                                   0 },
 {"PLoS Med.",                                                     0 },
 {"PLoS Negl Trop Dis",                                            0 },
 {"PLoS ONE",                                                      0 },
 {"PLoS Pathog.",                                                  0 },
 {"Redox Rep.",                                                 2004 },
 {"Reprod. Biol. Endocrinol.",                                     0 },
 {"Retrovirology",                                                 0 },
 {"Saline Syst.",                                                  0 },
 {"Sci. STKE",                                                     0 },
 {"ScientificWorldJournal",                                        0 },
 {"Tech. Tips Online",                                             0 },
 {"Virol. J.",                                                     0 },
 {NULL,                                                            0 }
};

NLM_EXTERN Boolean Mla2IsEPubOnlyJournal (
  CharPtr jta,
  Int2Ptr starting_yearP
)

{
  EjourDataPtr  ejp;
  Int2          L, R, mid;

  if (starting_yearP != NULL) {
    *starting_yearP = 0;
  }

  if (StringHasNoText (jta)) return FALSE;

  L = 0;
  R = sizeof (mla2_ejour_list) / sizeof (mla2_ejour_list [0]) - 1; /* -1 because now NULL terminated */

  while (L < R) {
    mid = (L + R) / 2;
    ejp = &(mla2_ejour_list [mid]);
    if (ejp != NULL && StringICmp (ejp->journal_name, jta) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  ejp = &(mla2_ejour_list [R]);
  if (ejp != NULL && StringICmp (ejp->journal_name, jta) == 0) {
    if (starting_yearP != NULL) {
      *starting_yearP = ejp->starting_year;
    }
    return TRUE;
  }

  return FALSE;
}

NLM_EXTERN Int2 Mla2CorrectCitArt (
  CitArtPtr newcap,
  CitArtPtr oldcap
)

{
  if (newcap == NULL) return 0;

  ChangeCitArtMLAuthorsToSTD (newcap);

  return 0;
}

