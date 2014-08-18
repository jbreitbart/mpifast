/*   tax3api.c
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
* File Name:  tax3api.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/8/04
*
* $Revision: 1.46 $
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
#include <objseq.h>
#include <objsset.h>
#include <tax3api.h>
#include <sqnutils.h>
#include <subutil.h>
#include <findrepl.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

/* low-level connection functions */

static Boolean text_tax_asn = FALSE;
static Boolean text_tax_set = FALSE;

#if 1
static const CharPtr tax3servicename = "TaxService3";
#else
static const CharPtr tax3servicename = "TaxService3Test";
#endif

NLM_EXTERN CONN Tax3OpenConnection (
  void
)

{
#ifdef OS_UNIX
  CharPtr  str;

  if (! text_tax_set) {
    str = (CharPtr) getenv ("TEXT_TAX_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        text_tax_asn = TRUE;
      }
    }
    text_tax_set = TRUE;
  }
#endif

  return QUERY_OpenServiceQuery (text_tax_asn ? "TaxService3Text" : tax3servicename, NULL, 30);
}

#ifdef OS_MAC
#include <Events.h>
#endif

NLM_EXTERN Taxon3ReplyPtr Tax3WaitForReply (
  CONN conn
)

{
  AsnIoConnPtr    aicp;
  time_t          currtime, starttime;
  time_t          max = 0;
  EIO_Status      status;
  STimeout        timeout;
  Taxon3ReplyPtr  t3ry = NULL;
#ifdef OS_MAC
  EventRecord     currEvent;
#endif

  if (conn == NULL) return NULL;

#ifdef OS_MAC
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 300;
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
    aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "r" : "rb", conn);
    t3ry = Taxon3ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  CONN_Close (conn);

  return t3ry;
}

/* high-level connection functions */

NLM_EXTERN Taxon3ReplyPtr Tax3SynchronousQuery (
  Taxon3RequestPtr t3rq
)

{
  AsnIoConnPtr    aicp;
  CONN            conn;
  Taxon3ReplyPtr  t3ry;

  if (t3rq == NULL) return NULL;

  conn = Tax3OpenConnection ();

  if (conn == NULL) return NULL;

  aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "w" : "wb", conn);

  Taxon3RequestAsnWrite (t3rq, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  t3ry = Tax3WaitForReply (conn);

  return t3ry;
}

NLM_EXTERN Boolean Tax3AsynchronousQuery (
  Taxon3RequestPtr t3rq,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
)

{
  AsnIoConnPtr  aicp;
  CONN          conn;

  if (t3rq == NULL) return FALSE;

  conn = Tax3OpenConnection ();

  if (conn == NULL) return FALSE;

  aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "w" : "wb", conn);

  Taxon3RequestAsnWrite (t3rq, aicp->aip, NULL);

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

NLM_EXTERN Int4 Tax3CheckQueue (
  QUEUE* queue
)

{
  return QUERY_CheckQueue (queue);
}

NLM_EXTERN Taxon3ReplyPtr Tax3ReadReply (
  CONN conn,
  EIO_Status status
)

{
  AsnIoConnPtr    aicp;
  Taxon3ReplyPtr  t3ry = NULL;

  if (conn != NULL && status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen (text_tax_asn ? "r" : "rb", conn);
    t3ry = Taxon3ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  return t3ry;
}

NLM_EXTERN Taxon3RequestPtr CreateTaxon3Request (
  Int4 taxid,
  CharPtr name,
  OrgRefPtr orp
)

{
  Taxon3RequestPtr  t2rp;

  t2rp = Taxon3RequestNew ();
  if (t2rp == NULL) return NULL;

  if (StringDoesHaveText (name)) {
    ValNodeCopyStr (&(t2rp->request), 2, name);
  } else if (taxid > 0) {
    ValNodeAddInt (&(t2rp->request), 1, taxid);
  } else if (orp != NULL) {
    orp = AsnIoMemCopy ((Pointer) orp,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
    ValNodeAddPointer (&(t2rp->request), 3, (Pointer) orp);
  }

  return t2rp;
}

NLM_EXTERN Taxon3RequestPtr CreateMultiTaxon3Request (ValNodePtr org_list)
{
  ValNodePtr vnp;
  Taxon3RequestPtr t3rp;
  OrgRefPtr orp;
  
  t3rp = Taxon3RequestNew ();
  if (t3rp == NULL) return NULL;

  for (vnp = org_list; vnp != NULL; vnp = vnp->next)
  {
    switch (vnp->choice)
    {
      case 1:
        ValNodeAddInt (&(t3rp->request), 1, vnp->data.intvalue);
        break;
      case 2:
        ValNodeCopyStr (&(t3rp->request), 2, vnp->data.ptrvalue);
        break;
      case 3:
        orp = AsnIoMemCopy (vnp->data.ptrvalue,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
        ValNodeAddPointer (&(t3rp->request), 3, (Pointer) orp);
        break;
    }
  }
  return t3rp;
}


static Boolean HasMisspellingFlag (T3DataPtr t)
{
  T3StatusFlagsPtr status;

  if (t == NULL) return FALSE;
  status = t->status;
  while (status != NULL) {
    if (StringCmp (status->property, "misspelled_name") == 0) {
      return TRUE;
    }
    status = status->next;
  }
  return FALSE;
}


static int LIBCALLBACK SortVnpByOrgRef (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return OrgRefCompare (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


NLM_EXTERN ValNodePtr Taxon3GetOrgRefList (ValNodePtr org_list)
{
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  OrgRefPtr        t3orp = NULL;
  T3ReplyPtr       trp;
  T3ErrorPtr       tep;
  ValNodePtr       uniq_list, response_list = NULL, next_org_list, last_org;
  Int4             request_num, max_requests = 2000;
  ValNodePtr PNTR  ptr_array;
  ValNodePtr       vnp, vnp_rq, vnp_rp;
  Int4             i, num_orgs;

  if (org_list == NULL) {
    return NULL;
  }

  /* make a copy of the original list - we will prepare the response list by substituting the OrgRef */
  org_list = ValNodeCopyPtr (org_list);

  /* make array to show original order of ValNodes, so that we can restore after sorting */
  num_orgs = ValNodeLen (org_list);
  ptr_array = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_orgs);
  for (vnp = org_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
    ptr_array[i] = vnp;
  }

  org_list = ValNodeSort (org_list, SortVnpByOrgRef);

  /* now make a list of just the unique requests */
  uniq_list = ValNodeCopyPtr (org_list);
  ValNodeUnique (&uniq_list, SortVnpByOrgRef, ValNodeFree);
  
  /* now break large lists into manageable chunks */
  vnp = uniq_list;
  while (vnp != NULL) {
    next_org_list = vnp->next;
    last_org = vnp; 
    request_num = 1;
    while (next_org_list != NULL && request_num < max_requests) {
      last_org = next_org_list;
      next_org_list = next_org_list->next;
      request_num++;
    }
    if (last_org != NULL) {
      last_org->next = NULL;
    }
      
    /* now create the request */
  
    t3rq = CreateMultiTaxon3Request (vnp);
    if (t3rq == NULL) return NULL;
    t3ry = Tax3SynchronousQuery (t3rq);
    Taxon3RequestFree (t3rq);
    if (t3ry != NULL) {
      for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
        switch (trp->choice) {
          case T3Reply_error :
            tep = (T3ErrorPtr) trp->data.ptrvalue;
            if (tep != NULL) {
              ErrPostEx (SEV_ERROR, 0, 0, tep->message);
            }
            if (tep != NULL && StringStr (tep->message, "ambiguous") != NULL) {
              ValNodeAddPointer (&response_list, eReturnedOrgFlag_ambiguous, NULL);
            } else {
              ValNodeAddPointer (&response_list, eReturnedOrgFlag_error, NULL);
            }
            break;
          case T3Reply_data :
            tdp = (T3DataPtr) trp->data.ptrvalue;
            if (tdp != NULL) {
              t3orp = (OrgRefPtr)(tdp->org);
              if (HasMisspellingFlag (tdp)) {
                ValNodeAddPointer (&response_list, eReturnedOrgFlag_misspelled, (Pointer) t3orp);
              } else {
                ValNodeAddPointer (&response_list, eReturnedOrgFlag_normal, (Pointer) t3orp);
              }
              tdp->org = NULL;
            }
            break;
          default :
            break;
        }
      }
      Taxon3ReplyFree (t3ry);
    }
    
    if (last_org != NULL) {
        last_org->next = next_org_list;
    }
    vnp = next_org_list;
  }  
  
  /* now put responses in list */
  vnp = uniq_list;
  vnp_rq = org_list;
  vnp_rp = response_list;

  while (vnp != NULL && vnp_rq != NULL && vnp_rp != NULL) {
    while (vnp_rq != NULL && OrgRefCompare (vnp->data.ptrvalue, vnp_rq->data.ptrvalue) == 0) {
      vnp_rq->data.ptrvalue = AsnIoMemCopy (vnp_rp->data.ptrvalue, (AsnReadFunc) OrgRefAsnRead, (AsnWriteFunc) OrgRefAsnWrite);
      vnp_rq->choice = vnp_rp->choice;
      vnp_rq = vnp_rq->next;
    }
    vnp_rp->data.ptrvalue = OrgRefFree (vnp_rp->data.ptrvalue);
    vnp_rp = vnp_rp->next;
    vnp = vnp->next;
  }
  /* if there were more requests than responses, set responses to NULL */
  while (vnp_rq != NULL) {
    vnp_rq->data.ptrvalue = NULL;
    vnp_rq = vnp_rq->next;
  }
  /* if there were more responses than requests, free extra responses */
  while (vnp_rp != NULL) {
    vnp_rp->data.ptrvalue = OrgRefFree (vnp_rp->data.ptrvalue);
    vnp_rp = vnp_rp->next;
  }
  response_list = ValNodeFree (response_list);
  uniq_list = ValNodeFree (uniq_list);

  /* now restore original order */
  for (i = 0; i < num_orgs - 1; i++) {
    ptr_array[i]->next = ptr_array[i + 1];
  }
  ptr_array[num_orgs - 1]->next = NULL;
  org_list = ptr_array[0];
  ptr_array = MemFree (ptr_array);
  
  return org_list;
}


NLM_EXTERN TaxFixItemPtr TaxFixItemNew (void)
{
  TaxFixItemPtr t;

  t = (TaxFixItemPtr) MemNew (sizeof (TaxFixItemData));
  MemSet (t, 0, sizeof (TaxFixItemData));
  return t;
}


NLM_EXTERN TaxFixItemPtr TaxFixItemCopy (TaxFixItemPtr orig)
{
  TaxFixItemPtr t = NULL;

  if (orig != NULL) {
    t = (TaxFixItemPtr) MemNew (sizeof (TaxFixItemData));
    t->data_choice = orig->data_choice;
    t->data = orig->data;
    t->response_org = AsnIoMemCopy (orig->response_org, (AsnReadFunc) OrgRefAsnRead, (AsnWriteFunc) OrgRefAsnWrite);
    if (orig->taxname != NULL) {
      t->taxname = StringSave (orig->taxname);
    }
    if (orig->suggested_fix != NULL) {
      t->suggested_fix = StringSave (orig->suggested_fix);
    }
    if (orig->rank != NULL) {
      t->rank = StringSave (orig->rank);
    }
  }
  return t;
}


NLM_EXTERN TaxFixItemPtr TaxFixItemFree (TaxFixItemPtr t)
{
  if (t != NULL) {
    t->response_org = OrgRefFree (t->response_org);
    t->taxname = MemFree (t->taxname);
    t->suggested_fix = MemFree (t->suggested_fix);
    t->rank = MemFree (t->rank);
    t = MemFree (t);
  }
  return t;
}


NLM_EXTERN ValNodePtr LIBCALLBACK TaxFixItemListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = TaxFixItemFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static Boolean LIBCALLBACK TaxFixItemOrigIsOk (ValNodePtr vnp)
{
  TaxFixItemPtr t;

  if (vnp == NULL || (t = (TaxFixItemPtr) vnp->data.ptrvalue) == NULL) {
    return TRUE;
  } else if (StringCmp (t->taxname, t->suggested_fix) == 0 && StringICmp (t->rank, "species") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr StringSum (CharPtr str1, CharPtr str2) 
{
  CharPtr sum = NULL;

  if (str1 == NULL && str2 == NULL) {
    sum = NULL;
  } else if (str1 == NULL) {
    sum = StringSave (str2);
  } else if (str2 == NULL) {
    sum = StringSave (str1);
  } else {
    sum = (CharPtr) MemNew (sizeof (Char) * (StringLen (str1) + StringLen (str2) + 1));
    sprintf (sum, "%s%s", str1, str2);
  }
  return sum;
}


static CharPtr SuggestedTaxNameFixFromOrgAndRank (CharPtr taxname, OrgRefPtr response_org, CharPtr rank)
{
  CharPtr fix = NULL, tmp;

  if (response_org == NULL) {
    return NULL;
  }

  if (StringICmp (rank, "species") == 0) {
    fix = StringSave (response_org->taxname);
  } else if (response_org->orgname != NULL) {
    if (((StringNICmp (taxname, "uncultured ", 11) == 0
           && StringICmp (taxname + 11, response_org->taxname) == 0)
         || StringICmp (taxname, response_org->taxname) == 0)
        && (StringISearch (response_org->orgname->lineage, "archaea") != NULL
            || StringISearch (response_org->orgname->lineage, "bacteria") != NULL)) {
      if (StringICmp (rank, "genus") == 0) {
        fix = StringSum (response_org->taxname, " sp.");
      } else if (StringNICmp (response_org->orgname->lineage, "bacteria", 8) == 0) {
        fix = StringSum (response_org->taxname, " bacterium");
      } else if (StringNICmp (response_org->orgname->lineage, "Archaea", 7) == 0) {
        fix = StringSum (response_org->taxname, " archaeon");
      }
      if (fix != NULL 
          && StringNICmp (fix, "uncultured ", 11) != 0) {
        tmp = fix;
        fix = StringSum ("uncultured ", tmp);
        tmp = MemFree (tmp);
      }
    }
  }
  return fix;
}


static ValNodePtr MakeTaxFixRequestList (ValNodePtr biop_list)
{
  ValNodePtr rq_list = NULL, prev = NULL;
  BioSourcePtr biop;
  OrgRefPtr  org;
  CharPtr    new_name;
  Int4       len;

  while (biop_list != NULL) {
    biop = GetBioSourceFromObject (biop_list->choice, biop_list->data.ptrvalue);
    org = AsnIoMemCopy (biop->org, (AsnReadFunc) OrgRefAsnRead, (AsnWriteFunc) OrgRefAsnWrite);
    if ((len = StringLen (org->taxname)) > 3 && StringCmp (org->taxname + len - 3, " sp") == 0) {
      new_name = StringSum (org->taxname, ".");
      org->taxname = MemFree (org->taxname);
      org->taxname = new_name;
    }

    ValNodeAddPointer (&prev, 3, org);
    if (rq_list == NULL) {
      rq_list = prev;
    }
    biop_list = biop_list->next;
  }
  return rq_list;
}


static void CheckSuggestedFixes (ValNodePtr tax_fix_list) 
{
  ValNodePtr rq_list = NULL, rp_list = NULL, prev, vnp_rq, vnp_rp;
  ValNodePtr vnp, next_org_list, last_org;
  Int4             request_num, max_requests = 2000;
  TaxFixItemPtr t;
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  T3ReplyPtr       trp;
  T3ErrorPtr       tep;
  T3StatusFlagsPtr tfp;
  OrgRefPtr        org;
  Boolean          is_species;

  prev = NULL;
  for (vnp = tax_fix_list; vnp != NULL; vnp = vnp->next) {
    t = (TaxFixItemPtr) vnp->data.ptrvalue;
    if (t != NULL && t->suggested_fix != NULL) {
      ValNodeAddPointer (&prev, 2, StringSave (t->suggested_fix));
    }
    if (rq_list == NULL) {
      rq_list = prev;
    }
  }

  /* now break large lists into manageable chunks */
  vnp = rq_list;
  while (vnp != NULL) {
    next_org_list = vnp->next;
    last_org = vnp; 
    request_num = 1;
    while (next_org_list != NULL && request_num < max_requests) {
      last_org = next_org_list;
      next_org_list = next_org_list->next;
      request_num++;
    }
    if (last_org != NULL) {
      last_org->next = NULL;
    }
      
    /* now create the request */
  
    t3rq = CreateMultiTaxon3Request (vnp);
    if (t3rq == NULL) return;
    t3ry = Tax3SynchronousQuery (t3rq);
    Taxon3RequestFree (t3rq);
    if (t3ry != NULL) {
      for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
        switch (trp->choice) {
          case T3Reply_error :
            tep = (T3ErrorPtr) trp->data.ptrvalue;
            ValNodeAddPointer (&rp_list, 0, NULL);
            break;
          case T3Reply_data :
            tdp = (T3DataPtr) trp->data.ptrvalue;
            is_species = FALSE;
            if (tdp != NULL) {
              for (tfp = tdp->status; tfp != NULL; tfp = tfp->next) {
                if (StringICmp (tfp->property, "rank") == 0
                    && tfp->Value_value != NULL
                    && tfp->Value_value->choice == Value_value_str
                    && StringICmp (tfp->Value_value->data.ptrvalue, "species") == 0) {
                  is_species = TRUE;
                }
              }
            }
            if (is_species) {
              org = (OrgRefPtr) tdp->org;
              ValNodeAddPointer (&rp_list, 0, StringSave (org->taxname));
            } else {
              ValNodeAddPointer (&rp_list, 0, NULL);
            }
            break;
          default :
            break;
        }
      }
      Taxon3ReplyFree (t3ry);
    }
    
    if (last_org != NULL) {
        last_org->next = next_org_list;
    }
    vnp = next_org_list;
  }  
  rq_list = ValNodeFreeData (rq_list);

  /* adjust suggested fixes */
  vnp_rq = tax_fix_list;
  vnp_rp = rp_list;

  while (vnp_rq != NULL && vnp_rp != NULL) {
    while (vnp_rq != NULL && ((t = (TaxFixItemPtr) vnp_rq->data.ptrvalue) == NULL || t->suggested_fix == NULL)) {
      vnp_rq = vnp_rq->next;
    }
    if (t != NULL) {
      t->suggested_fix = MemFree (t->suggested_fix);
      t->suggested_fix = vnp_rp->data.ptrvalue;
      vnp_rp->data.ptrvalue = NULL;
      vnp_rq = vnp_rq->next;
      vnp_rp = vnp_rp->next;
    }
  }
  rp_list = ValNodeFreeData (rp_list);
}


NLM_EXTERN ValNodePtr Taxon3GetTaxFixList (ValNodePtr biop_list)
{
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  T3ReplyPtr       trp;
  T3ErrorPtr       tep;
  ValNodePtr       uniq_list, response_list = NULL, next_org_list, last_org, request_list;
  Int4             request_num, max_requests = 2000;
  ValNodePtr PNTR  ptr_array;
  ValNodePtr       vnp, vnp_rq, vnp_rp, vnp_b;
  T3StatusFlagsPtr tfp;
  Int4             i, num_orgs;
  TaxFixItemPtr    t;
  BioSourcePtr     biop;

  if (biop_list == NULL) {
    return NULL;
  }

  /* make a copy of the original list, removing uncultured */
  request_list = MakeTaxFixRequestList (biop_list);

  /* make array to show original order of ValNodes, so that we can restore after sorting */
  num_orgs = ValNodeLen (request_list);
  ptr_array = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_orgs);
  for (vnp = request_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
    ptr_array[i] = vnp;
  }

  request_list = ValNodeSort (request_list, SortVnpByOrgRef);

  /* now make a list of just the unique requests */
  uniq_list = ValNodeCopyPtr (request_list);
  ValNodeUnique (&uniq_list, SortVnpByOrgRef, ValNodeFree);
  
  /* now break large lists into manageable chunks */
  vnp = uniq_list;
  while (vnp != NULL) {
    next_org_list = vnp->next;
    last_org = vnp; 
    request_num = 1;
    while (next_org_list != NULL && request_num < max_requests) {
      last_org = next_org_list;
      next_org_list = next_org_list->next;
      request_num++;
    }
    if (last_org != NULL) {
      last_org->next = NULL;
    }
      
    /* now create the request */
  
    t3rq = CreateMultiTaxon3Request (vnp);
    if (t3rq == NULL) return NULL;
    t3ry = Tax3SynchronousQuery (t3rq);
    Taxon3RequestFree (t3rq);
    if (t3ry != NULL) {
      for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
        switch (trp->choice) {
          case T3Reply_error :
            tep = (T3ErrorPtr) trp->data.ptrvalue;
            t = TaxFixItemNew ();
            ValNodeAddPointer (&response_list, 0, t);
            break;
          case T3Reply_data :
            tdp = (T3DataPtr) trp->data.ptrvalue;
            if (tdp != NULL) {
              t = TaxFixItemNew ();
              t->response_org = (OrgRefPtr)(tdp->org);
              tdp->org = NULL;
              for (tfp = tdp->status; tfp != NULL; tfp = tfp->next) {
                if (StringICmp (tfp->property, "rank") == 0
                    && tfp->Value_value != NULL
                    && tfp->Value_value->choice == Value_value_str) {
                  t->rank = StringSave (tfp->Value_value->data.ptrvalue);
                }
              }
              t->taxname = StringSave (t->response_org->taxname);
              t->suggested_fix = SuggestedTaxNameFixFromOrgAndRank (t->taxname, t->response_org, t->rank);
              ValNodeAddPointer (&response_list, 0, t);
            }
            break;
          default :
            break;
        }
      }
      Taxon3ReplyFree (t3ry);
    }
    
    if (last_org != NULL) {
        last_org->next = next_org_list;
    }
    vnp = next_org_list;
  }  

  CheckSuggestedFixes (response_list);
  
  /* now put responses in list */
  vnp = uniq_list;
  vnp_rq = request_list;
  vnp_rp = response_list;

  while (vnp != NULL && vnp_rq != NULL && vnp_rp != NULL) {
    while (vnp_rq != NULL && OrgRefCompare (vnp->data.ptrvalue, vnp_rq->data.ptrvalue) == 0) {
      t = TaxFixItemCopy (vnp_rp->data.ptrvalue);
      vnp_rq->data.ptrvalue = t;
      vnp_rq = vnp_rq->next;
    }
    vnp_rp = vnp_rp->next;
    vnp = vnp->next;
  }
  /* if there were more requests than responses, set responses to NULL */
  while (vnp_rq != NULL) {
    vnp_rq->data.ptrvalue = NULL;
    vnp_rq = vnp_rq->next;
  }

  /* free response list */
  response_list = TaxFixItemListFree (response_list);

  uniq_list = ValNodeFree (uniq_list);

  /* now restore original order */
  for (i = 0; i < num_orgs - 1; i++) {
    ptr_array[i]->next = ptr_array[i + 1];
  }
  ptr_array[num_orgs - 1]->next = NULL;
  request_list = ptr_array[0];
  ptr_array = MemFree (ptr_array);

  /* now reassociate with original objects */
  for (vnp_b = biop_list, vnp_rp = request_list; vnp_b != NULL && vnp_rp != NULL; vnp_b = vnp_b->next, vnp_rp = vnp_rp->next) {
    t = vnp_rp->data.ptrvalue;
    t->data_choice = vnp_b->choice;
    t->data = vnp_b->data.ptrvalue;
    t->taxname = MemFree (t->taxname);
    biop = GetBioSourceFromObject (t->data_choice, t->data);
    if (biop != NULL && biop->org != NULL) {
      t->taxname = StringSave (biop->org->taxname);
    }
  }

  /* now remove items for which the original and suggested taxnames are the same */
  
  ValNodePurge (&request_list, TaxFixItemOrigIsOk, TaxFixItemListFree);
  return request_list;
}


NLM_EXTERN OrgRefPtr Taxon3GetOrg (OrgRefPtr orp)

{
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3DataPtr        tdp;
  OrgRefPtr        t3orp = NULL;
  T3ReplyPtr        trp;
  T3ErrorPtr        tep;
	
  if (orp == NULL) return NULL;
  
  t3rq = CreateTaxon3Request (0, NULL, orp);
  if (t3rq == NULL) return NULL;
  t3ry = Tax3SynchronousQuery (t3rq);
  Taxon3RequestFree (t3rq);
  if (t3ry != NULL) {
    for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
      switch (trp->choice) {
        case T3Reply_error :
          tep = (T3ErrorPtr) trp->data.ptrvalue;
          if (tep != NULL) {
            ErrPostEx (SEV_ERROR, 0, 0, tep->message);
          }
          break;
        case T3Reply_data :
          tdp = (T3DataPtr) trp->data.ptrvalue;
          if (tdp != NULL) {
            t3orp = (OrgRefPtr)(tdp->org);
            tdp->org = NULL;
          }
          break;
        default :
          break;
      }
    }
    Taxon3ReplyFree (t3ry);
  }
  
  return t3orp;
}

static Boolean DoOrgIdsMatch(BioSourcePtr b1, BioSourcePtr b2)
{
  DbtagPtr d1 = NULL, d2 = NULL;
  ValNodePtr vnp;
	
  if (b1 == NULL || b2 == NULL) 
  {
    return FALSE;
  }
  if (b1->org ==  NULL || b2->org == NULL) 
  {
    return FALSE;
  }
  for (vnp = b1->org->db; vnp; vnp = vnp->next) 
  {
    d1 = (DbtagPtr) vnp->data.ptrvalue;
    if (StringCmp(d1->db, "taxon") == 0) 
    {
      break;
    }
  }
  for (vnp = b2->org->db; vnp; vnp = vnp->next) 
  {
    d2 = (DbtagPtr) vnp->data.ptrvalue;
	if (StringCmp(d2->db, "taxon") == 0) 
	{
      break;
	}
  }
  if (d1 && d2) 
  {
	if (d1->tag->id == d2->tag->id) 
	{
      return TRUE;
	}
  }
  else if (StringICmp(b1->org->taxname, b2->org->taxname) == 0) 
  {
	return TRUE;
  }
  return FALSE;
}

static BioSourcePtr Tax3BioSourceMerge(BioSourcePtr host, BioSourcePtr guest)
{
  SubSourcePtr ssp, sp, last_ssp;
  OrgModPtr omp, homp, last_omp;
  OrgNamePtr	onp;
	
  if (host == NULL && guest == NULL) 
  {
    return NULL;
  }
  if (host == NULL && guest != NULL) 
  {
	host = AsnIoMemCopy(guest, (AsnReadFunc) BioSourceAsnRead, 
		   						(AsnWriteFunc) BioSourceAsnWrite);
	return host;
  }
  if (host != NULL && guest == NULL) 
  {
    return host;
  }
  if (host->genome == 0 && guest->genome != 0) 
  {
    host->genome = guest->genome;
  }
  if (host->origin == 0 && guest->origin != 0) 
  {
    host->origin = guest->origin;
  }
  last_ssp = host->subtype;
  while (last_ssp != NULL && last_ssp->next != NULL)
  {
  	last_ssp = last_ssp->next;
  }
  for (ssp = guest->subtype; ssp; ssp = ssp->next) 
  {
    sp = AsnIoMemCopy(ssp, (AsnReadFunc) SubSourceAsnRead, 
		   						(AsnWriteFunc) SubSourceAsnWrite);
    if (last_ssp == NULL)
    {
      host->subtype = sp;
    }
    else
    {
      last_ssp->next = sp;
      last_ssp = sp;
    }
  }
  if (guest->org->orgname) 
  {
   	if ((onp = host->org->orgname)	== NULL) 
   	{
   	  onp = OrgNameNew();
   	  host->org->orgname = onp;
    }	
    last_omp = onp->mod;		
    while (last_omp != NULL && last_omp->next != NULL)
    {
      last_omp = last_omp->next;
    }
    for (omp = guest->org->orgname->mod; omp; omp = omp->next) 
    {
      homp = AsnIoMemCopy(omp, (AsnReadFunc) OrgModAsnRead, 
		   						(AsnWriteFunc) OrgModAsnWrite);
      if (last_omp == NULL)
      {
      	onp->mod = homp;
      }
      else
      {
      	last_omp->next = homp;
      	last_omp = homp;
      }
    }
  }
  return host;
}


/**************************************************************************
*	Compare BioSources in one bioseq->descr using Taxonomy to find
*	their join parent
*	merge if organisms are the same or create a feature if different
*
**************************************************************************/
NLM_EXTERN void Tax3MergeSourceDescr (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr    bsp = NULL;
	ValNodePtr   vnp, newlist;
	SeqFeatPtr   sfp;
	BioSourcePtr first_biop = NULL;
	BioSourcePtr other_biop;
	BioSourcePtr tmp_biop;
	ObjValNodePtr ovp;

	if (!IS_Bioseq(sep)) {
		return;
	}
	newlist = (ValNodePtr) data;
	bsp = (BioseqPtr) sep->data.ptrvalue;
	if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const) 
			&& (bsp->repr != Seq_repr_delta))
		return;

	if (! ISA_na(bsp->mol))
		return;
	
	/* add the descriptors in newlist to the end of the list in bsp->descr*/
	if (bsp->descr == NULL)
	{
	  bsp->descr = newlist;
	}
	else
	{
	  for (vnp = bsp->descr; vnp->next != NULL; vnp = vnp->next)
	  {	
	  }
	  vnp->next = newlist;
	}
	
	/* now find the first source descriptor in bsp->descr that has an org*/
    /* note - we can't use SeqMgrGetNextDescriptor here because we have just
     * added to the descriptors, so they are not indexed. */
	for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
	{
	  if (vnp->choice != Seq_descr_source) continue;
	  if (vnp->data.ptrvalue == NULL)
	  {
	  	ErrPostStr(SEV_WARNING, 0, 0, "Source descriptor missing data");
	  	if (vnp->extended)
	  	{
	  	  ovp = (ObjValNodePtr) vnp;
	  	  ovp->idx.deleteme = TRUE;
	  	}
	  }
	  if (first_biop == NULL)
	  {
	  	first_biop = vnp->data.ptrvalue;
	  }
	  else
	  {
		other_biop = vnp->data.ptrvalue;
		/* detach biosource pointer from descr, so that it will not be freed
		 * when the descriptor is deleted.
		 */
		vnp->data.ptrvalue = NULL;
        if (vnp->extended)
        {
          ovp = (ObjValNodePtr) vnp;
	  	  ovp->idx.deleteme = TRUE;
        }
        if (DoOrgIdsMatch(first_biop, other_biop)) 
		{
		  /* merge the two sources */
		  tmp_biop = Tax3BioSourceMerge(first_biop, other_biop);
		  if (tmp_biop == NULL)
		  {
		  	ErrPostStr (SEV_WARNING, 0, 0, "Failed to merge biosources");
		  }
		  else
		  {
		  	first_biop = tmp_biop;
		  }
		  other_biop = BioSourceFree (other_biop);
		} else {
		  /* create a source feature */
		  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_BIOSRC, NULL);
		  if (sfp != NULL)
		  {
            sfp->data.value.ptrvalue = other_biop;
		  }
        }
	  }
	}
	return;
}

static Int4 GetTaxIdFromOrgRef (OrgRefPtr orp)
{
  Int4       tax_id = -1;
  ValNodePtr vnp;
  DbtagPtr   d;

  if (orp != NULL)
  {
    for (vnp = orp->db; vnp != NULL; vnp = vnp->next) 
    {
      d = (DbtagPtr) vnp->data.ptrvalue;
      if (StringCmp(d->db, "taxon") == 0) 
      {
        tax_id = d->tag->id;
        break;
      }
    }
  }
  return tax_id;
}

NLM_EXTERN Int4 Taxon3GetTaxIdByOrgRef (OrgRefPtr orp)
{
  OrgRefPtr  orp_repl;
  Int4       tax_id = -1;
  
  if (orp == NULL) return -1;
  
  orp_repl = Taxon3GetOrg (orp);
  tax_id = GetTaxIdFromOrgRef (orp_repl);
  OrgRefFree (orp_repl);
  
  return tax_id;
}

NLM_EXTERN OrgRefPtr Taxon3GetOrgRefByName (CharPtr orgname)
{
  OrgRefPtr request, org;
  
  request = OrgRefNew ();
  if (request == NULL) return NULL;
  request->taxname = orgname;
  org = Taxon3GetOrg (request);
  request->taxname = NULL;
  OrgRefFree (request);
  return org;
}

NLM_EXTERN Int4 Taxon3GetTaxIdByName (CharPtr orgname)
{
  OrgRefPtr orp;
  Int4      tax_id;
  
  orp = Taxon3GetOrgRefByName (orgname);
  tax_id = GetTaxIdFromOrgRef (orp);

  OrgRefFree(orp);
  return tax_id;
}

static void AddBioSourceToList (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr PNTR list;
  
  if (biop == NULL || userdata == NULL) return;
  list = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (list, 4, (Pointer) biop);
}

NLM_EXTERN void Taxon3ReplaceOrgInSeqEntry (SeqEntryPtr sep, Boolean keep_syn)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   biop_vnp, response_vnp;
  BioSourcePtr biop;
  OrgRefPtr    swap_org, response_org;
  
  VisitBioSourcesInSep (sep, &biop_list, AddBioSourceToList);

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    biop = (BioSourcePtr) biop_vnp->data.ptrvalue;
    ValNodeAddPointer (&request_list, 3, biop->org);
  }
  response_list = Taxon3GetOrgRefList (request_list);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
    return;
  }

  for (biop_vnp = biop_list, response_vnp = response_list;
       biop_vnp != NULL && response_vnp != NULL;
       biop_vnp = biop_vnp->next, response_vnp = response_vnp->next)
  {
    biop = (BioSourcePtr) biop_vnp->data.ptrvalue;
    swap_org = biop->org;
    response_org = response_vnp->data.ptrvalue;
    if (response_org != NULL)
    {
      biop->org = response_org;
      response_vnp->data.ptrvalue = NULL;
      OrgRefFree (swap_org);
      if (! keep_syn)
      {
        biop->org->syn = ValNodeFreeData(biop->org->syn);
      }
    }
  }
  ValNodeFree (request_list);
  ValNodeFree (response_list);
  ValNodeFree (biop_list);   
}


static void GetBioSourceFeaturesForCheck (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR list = (ValNodePtr PNTR) userdata;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || list == NULL
      || sfp->data.value.ptrvalue == NULL) {
    return;
  }
  ValNodeAddPointer (list, OBJ_SEQFEAT, sfp);
}


static void GetBioSourceDescriptorsForCheck (SeqDescrPtr sdp, Pointer userdata)
{
  ValNodePtr PNTR list = (ValNodePtr PNTR) userdata;
  if (sdp == NULL || sdp->choice != Seq_descr_source || list == NULL
      || sdp->data.ptrvalue == NULL) {
    return;
  }
  ValNodeAddPointer (list, OBJ_SEQDESC, sdp);
}


static DbtagPtr GetTaxonXref (OrgRefPtr org)
{
  ValNodePtr vnp;
  DbtagPtr   dbt = NULL;
  
  if (org == NULL) return NULL;
  vnp = org->db;
  while (vnp != NULL && dbt == NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringICmp ((CharPtr) dbt->db, "taxon") != 0) {
      dbt = NULL;
    }
    vnp = vnp->next;
  }
  return dbt;
}
  
static Boolean DoTaxonIdsMatch (OrgRefPtr org1, OrgRefPtr org2)
{
  DbtagPtr   dbt1 = NULL, dbt2 = NULL;
  
  if (org1 == NULL || org2 == NULL) return FALSE;
  
  dbt1 = GetTaxonXref (org1);
  if (dbt1 == NULL) return FALSE;
  dbt2 = GetTaxonXref (org2);
  if (dbt2 == NULL) return FALSE;
  
  return DbtagMatch(dbt1, dbt2);
}


NLM_EXTERN void Taxon3CheckOrgInSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR not_found, ValNodePtr PNTR bad_match)
{
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   biop_vnp, response_vnp;
  BioSourcePtr biop;
  OrgRefPtr    orig_org, response_org;
  ValNodePtr   item_list = NULL;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  
  VisitFeaturesInSep (sep, &item_list, GetBioSourceFeaturesForCheck);
  VisitDescriptorsInSep (sep, &item_list, GetBioSourceDescriptorsForCheck);
  
  for (biop_vnp = item_list; biop_vnp != NULL; biop_vnp = biop_vnp->next) {
    biop = NULL;
    if (biop_vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) biop_vnp->data.ptrvalue;  
      if (sfp != NULL) {  
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;      
      }
    } else if (biop_vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) biop_vnp->data.ptrvalue;
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop != NULL) {
      ValNodeAddPointer (&request_list, 3, biop->org);
    }
  }

  response_list = Taxon3GetOrgRefList (request_list);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
    ValNodeFree (request_list);
    ValNodeFree (item_list);
    return;
  }

  for (biop_vnp = item_list, response_vnp = response_list;
       biop_vnp != NULL && response_vnp != NULL;
       biop_vnp = biop_vnp->next, response_vnp = response_vnp->next)
  {
    response_org = response_vnp->data.ptrvalue;  
    biop = NULL;
    orig_org = NULL;
    if (biop_vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) biop_vnp->data.ptrvalue;    
      if (sfp != NULL) {  
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      }
    } else if (biop_vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) biop_vnp->data.ptrvalue;
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop == NULL) {
      Message (MSG_POST, "Error collecting data");
      ValNodeFree (request_list);
      ValNodeFree (item_list);
      return;
    } else {
      orig_org = biop->org;
      if (orig_org != NULL) {
        if (response_org == NULL) {
          ValNodeAddPointer (not_found, biop_vnp->choice, biop_vnp->data.ptrvalue);          
        } else if (StringCmp (orig_org->taxname, response_org->taxname) != 0) {
          ValNodeAddPointer (bad_match, biop_vnp->choice, biop_vnp->data.ptrvalue);
        } else if (!DoTaxonIdsMatch(orig_org, response_org)) {
          ValNodeAddPointer (bad_match, biop_vnp->choice, biop_vnp->data.ptrvalue);
        }        
      }
    }
    OrgRefFree (response_org);
  }
  ValNodeFree (request_list);
  ValNodeFree (response_list);
  ValNodeFree (item_list);   
}


NLM_EXTERN void CheckTaxNamesAgainstTaxDatabase (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp;
  SeqEntryPtr sep;
  SeqEntryPtr orig_scope;
  ValNodePtr  not_found = NULL, bad_match = NULL;
  CharPtr     bad_match_fmt = "%d tax names do not match taxonomy lookup.";
  CharPtr     no_match_fmt = "%d organisms are not found in taxonomy lookup.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  
  orig_scope = SeqEntryGetScope ();
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    SeqEntrySetScope (sep);
    Taxon3CheckOrgInSeqEntry (sep, &not_found, &bad_match);
  }
  SeqEntrySetScope (orig_scope);
  if (not_found != NULL) {
    dip = NewClickableItem (DISC_NO_TAXLOOKUP, no_match_fmt, not_found);
    dip->subcategories = NULL;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
  if (bad_match != NULL) {
    dip = NewClickableItem (DISC_BAD_TAXLOOKUP, bad_match_fmt, bad_match);
    dip->subcategories = NULL;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static ValNodePtr FreeOrgRefValNodeList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  OrgRefPtr  org;

  while (vnp != NULL)
  { 
    vnp_next = vnp->next;
    vnp->next = NULL;
    org = (OrgRefPtr) vnp->data.ptrvalue;
    vnp->data.ptrvalue = OrgRefFree (org);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static Boolean EndsWithSp (CharPtr str)
{
  Int4 len;

  if (StringHasNoText (str)) return FALSE;
  len = StringLen (str);
  if (len < 4) return FALSE;
  if (StringCmp (str + len - 4, " sp.") == 0) return TRUE;
  return FALSE;
}


static CharPtr RemoveSp (CharPtr orig)
{
  CharPtr cpy = NULL;
  Int4    len;

  len = StringLen (orig);
  if (len >= 4 && StringCmp (orig + len - 4, " sp.") == 0) {
    cpy = (CharPtr) MemNew (sizeof (Char) * len - 3);
    StringNCpy (cpy, orig, len - 4);
    cpy[len - 4] = 0;
  }
  return cpy;
}

  
static void AddRequestOrgForString (CharPtr str, CharPtr host, ValNodePtr PNTR request_list, ValNodePtr PNTR req_host_list)
{
  OrgRefPtr    request_org;
  CharPtr      cp, cpy;

  if (StringHasNoText (str) || host == NULL || request_list == NULL || req_host_list == NULL)
  {
    return;
  }

  /* if ends with " sp.", remove " sp." */
  cpy = RemoveSp (host);
  if (cpy != NULL) {
    request_org = OrgRefNew();
    request_org->taxname = StringSave (cpy);
    ValNodeAddPointer (request_list, 3, request_org);
    ValNodeAddPointer (req_host_list, 0, StringSave (host));
  } else {
    request_org = OrgRefNew();
    request_org->taxname = StringSave (str);
    ValNodeAddPointer (request_list, 3, request_org);
    ValNodeAddPointer (req_host_list, 0, StringSave (host));

     
    /* if more than one word, try chopping off last to see if abbreviated name looks up */
    cp = StringRChr (str, ' ');
    if (cp != NULL)
    {
      cpy = StringSave (str);    
      cp = StringRChr (cpy, ' ');
      if (cp != NULL)
      {
        *cp = 0;
        AddRequestOrgForString (cpy, host, request_list, req_host_list);
      }
      cpy = MemFree (cpy);
    }
  }
}

typedef struct specifichostcheck {
  CharPtr      spec_host;
  ValNodePtr   request_list;  /* ValNodeList of orgs */
  ValNodePtr   response_list; /* ValNodeList of orgs */
  ValNodePtr   biop_list;     /* ValNodeList of sources with this spec_host value */
} SpecificHostCheckData, PNTR SpecificHostCheckPtr;


static ValNodePtr SpecificHostCheckListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  SpecificHostCheckPtr p;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    p = (SpecificHostCheckPtr) vnp->data.ptrvalue;
    if (p != NULL)
    {
      p->request_list = FreeOrgRefValNodeList (p->request_list);
      p->response_list = FreeOrgRefValNodeList (p->response_list);
      p->spec_host = MemFree (p->spec_host);
      p->biop_list = ValNodeFree (p->biop_list);
    }
    vnp = ValNodeFreeData (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static ValNodePtr SortSpecificHostOrgs (ValNodePtr host_list, ValNodePtr request_list, ValNodePtr response_list)
{
  ValNodePtr           check_list = NULL;
  SpecificHostCheckPtr p = NULL;
  CharPtr              host, prev_host = NULL;

  while (host_list != NULL
         && request_list != NULL
         && response_list != NULL)
  {
    host = (CharPtr) host_list->data.ptrvalue;
    if (StringCmp (host, prev_host) != 0)
    {
      p = (SpecificHostCheckPtr) MemNew (sizeof (SpecificHostCheckData));
      p->spec_host = StringSave (host);
      ValNodeAddPointer (&check_list, 0, p);
      prev_host = host;
    }
    ValNodeAddPointer (&(p->request_list), request_list->choice, request_list->data.ptrvalue);
    ValNodeAddPointer (&(p->response_list), response_list->choice, response_list->data.ptrvalue);
    request_list->data.ptrvalue = NULL;
    response_list->data.ptrvalue = NULL;
    host_list = host_list->next;
    request_list = request_list->next;
    response_list = response_list->next;
  }
  return check_list;        
}


static Boolean StringAlreadyInValNodeList (CharPtr str, ValNodePtr list) 
{
  if (StringHasNoText (str))
  {
    return TRUE;
  }
  
  while (list != NULL)
  {
    if (StringCmp (str, list->data.ptrvalue) == 0)
    {
      return TRUE;
    }
    list = list->next;
  }
  return FALSE;
}


static BioSourcePtr GetBioSourceFromValNode (ValNodePtr vnp)
{
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioSourcePtr biop = NULL;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  } 
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  return biop;
}


static CharPtr extract_list[] = {
  "cf.",
  "cf ",
  "aff ",
  "aff.",
  "near",
  "nr.",
  "nr ",
  NULL};

static void AdjustSpecificHostForTaxServer (CharPtr spec_host)
{
  CharPtr cp, src, dst;
  Int4 i;

  /* ignore separator words */
  for (i = 0; extract_list[i] != NULL; i++) {
    if ((cp = StringSearch (spec_host, extract_list[i])) != NULL && cp > spec_host && isspace (*(cp - 1))) {
      src = cp + StringLen (extract_list[i]);
      dst = cp;
      while (isspace (*src)) {
        src++;
      }
      while (*src != 0) {
        *dst = *src;
        dst++;
        src++;
      }
      *dst = 0;
    }
  }
}


static void AddBioSourcesToSpecificHostChecklist (ValNodePtr biop_list, ValNodePtr check_list)
{
  ValNodePtr biop_vnp, last_vnp = NULL, stop_search;
  BioSourcePtr biop;
  OrgModPtr    mod;
  SpecificHostCheckPtr p;
  CharPtr tmp;

  if (biop_list == NULL || check_list == NULL) return;

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {

    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL) continue;

    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    mod = biop->org->orgname->mod;
    while (mod != NULL)
    {
      if (mod->subtype == ORGMOD_nat_host
          && !StringHasNoText (mod->subname))
      {
        if (last_vnp == NULL)
        {
          last_vnp = check_list;
          stop_search = NULL;
        }
        else
        {
          stop_search = last_vnp;
        }
        tmp = StringSave (mod->subname);
        AdjustSpecificHostForTaxServer (tmp);
        p = NULL;
        while (last_vnp != NULL 
               && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
               && StringCmp (p->spec_host, tmp) != 0)
        {
          p = NULL;
          last_vnp = last_vnp->next;
        }
        if (p == NULL && stop_search != NULL)
        {
          last_vnp = check_list;
          while (last_vnp != stop_search 
                 && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
                 && StringCmp (p->spec_host, tmp) != 0)
          {
            p = NULL;
            last_vnp = last_vnp->next;
          }
        }
        tmp = MemFree (tmp);
        if (p != NULL)
        {
          ValNodeAddPointer (&(p->biop_list), biop_vnp->choice, biop_vnp->data.ptrvalue);
        }
      }
      mod = mod->next;
    }
  }
}


static Boolean ShouldCheckSpecificHostValueForValidator (CharPtr spec_host)
{
  if (StringHasNoText (spec_host) || !isupper (*spec_host)) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static CharPtr GetSpecificHostValueToCheckForValidator (CharPtr spec_host)
{
  CharPtr cp, check_val = NULL;
  Int4    len = 0;

  if (ShouldCheckSpecificHostValueForValidator(spec_host)) {
    cp = spec_host;
    /* skip first word */
    while (*cp != 0 && !isspace (*cp)) {
      cp++;
      len++;
    }
    while (isspace (*cp)) {
      cp++;
      len++;
    }
    if (*cp != '(' && StringNCmp (cp, "sp.", 3) != 0 && *cp != 0) {
      /* collect second word */
      while (*cp != 0 && !isspace (*cp)) {
        cp++;
        len++;
      }
    }
    check_val = (CharPtr) MemNew (sizeof (Char) * (len + 1));
    StringNCpy (check_val, spec_host, len);
    check_val[len] = 0;
    TrimSpacesAroundString (check_val);
  }
  return check_val;
}

static Boolean ShouldCheckSpecificHostInBioSource (BioSourcePtr biop)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }
  for (mod = biop->org->orgname->mod; mod != NULL && !rval; mod = mod->next) {
    if (mod->subtype == ORGMOD_nat_host) {
      rval = ShouldCheckSpecificHostValueForValidator (mod->subname);
    }
  }
  return rval;
}



static void AddValidatorSpecificHostBioSourceFeatToList (SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;

  if (ShouldCheckSpecificHostInBioSource (sfp->data.value.ptrvalue))
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
  }
}


static void AddValidatorSpecificHostBioSourceDescToList (SeqDescrPtr sdp, Pointer userdata)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  if (ShouldCheckSpecificHostInBioSource (sdp->data.ptrvalue))
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr GetValidatorSpecificHostBioSourceList (SeqEntryPtr sep)
{
  ValNodePtr list = NULL;

  VisitFeaturesInSep (sep, &list, AddValidatorSpecificHostBioSourceFeatToList);
  VisitDescriptorsInSep (sep, &list, AddValidatorSpecificHostBioSourceDescToList);
  return list;
}


static void 
FormatValidatorSpecificHostRequests 
(ValNodePtr spec_host_list,
 ValNodePtr PNTR request_list,
 ValNodePtr PNTR req_host_list)
{
  ValNodePtr vnp;
  CharPtr    orig;
  OrgRefPtr  request_org;
  
  /* now format requests for unique specific_host values */
  for (vnp = spec_host_list; vnp != NULL; vnp = vnp->next)
  {
    orig = (CharPtr) vnp->data.ptrvalue;
    request_org = OrgRefNew();
    request_org->taxname = GetSpecificHostValueToCheckForValidator (orig);
    ValNodeAddPointer (request_list, 3, request_org);
    ValNodeAddPointer (req_host_list, 0, StringSave (orig));    
  }
}

static Boolean MatchesSynonym (CharPtr txt, OrgRefPtr response_org)
{
  ValNodePtr syn;
  Boolean    rval = FALSE;
  if (StringHasNoText (txt) || response_org == NULL) return FALSE;

  for (syn = response_org->syn; syn != NULL && !rval; syn = syn->next)
  {
    if (StringCmp (txt, syn->data.ptrvalue) == 0)
    {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean MatchesGenBankSynonym (CharPtr txt, OrgRefPtr response_org)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (StringHasNoText (txt) || response_org == NULL || response_org->orgname == NULL) return FALSE;
  mod = response_org->orgname->mod;
  while (mod != NULL) 
  {
    if ((mod->subtype == ORGMOD_gb_synonym || mod->subtype == ORGMOD_old_name) && StringCmp (txt, mod->subname) == 0)
    {
      rval = TRUE;
    }
    mod = mod->next;
  }
  return rval;
}


static ValNodePtr GetListOfUniqueSpecificHostValues (ValNodePtr biop_list)
{
  ValNodePtr   biop_vnp;
  BioSourcePtr biop;
  OrgModPtr    mod;
  ValNodePtr   spec_host_list = NULL;
  CharPtr      tmp;
  
  /* get a list of unique specific_host values */
  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    if (biop_vnp->data.ptrvalue == NULL) continue;
    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    mod = biop->org->orgname->mod;
    while (mod != NULL)
    {
      if (mod->subtype == ORGMOD_nat_host
          && !StringHasNoText (mod->subname))
      {
        tmp = StringSave (mod->subname);
        AdjustSpecificHostForTaxServer (tmp);
        ValNodeAddPointer (&spec_host_list, 0, tmp);
      }
      mod = mod->next;
    }
  }
  spec_host_list = ValNodeSort (spec_host_list, SortVnpByString);
  ValNodeUnique (&spec_host_list, SortVnpByString, ValNodeFreeData);
  return spec_host_list;
}


static Boolean StringIsExactMatchForOrgRef (CharPtr str, OrgRefPtr org)
{
  if (StringHasNoText (str) || org == NULL) {
    return FALSE;
  } else if (StringCmp (org->taxname, str) == 0
             || StringCmp (org->common, str) == 0 
             || MatchesSynonym (str, org) 
             || MatchesGenBankSynonym (str, org)) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static CharPtr FindMatchInOrgRef (CharPtr str, OrgRefPtr org)
{
  ValNodePtr syn;
  OrgModPtr  mod;
  CharPtr    rval = NULL;

  if (StringHasNoText (str) || org == NULL) {
    rval = NULL;
  } else if (StringICmp (org->taxname, str) == 0) {
    rval = org->taxname;
  } else if (StringICmp (org->common, str) == 0) {
    rval = org->common;
  } else {
    for (syn = org->syn; syn != NULL && rval == NULL; syn = syn->next) {
      if (StringICmp (str, syn->data.ptrvalue) == 0) {
        rval = syn->data.ptrvalue;
      }
    }
    if (org->orgname != NULL) {
      for (mod = org->orgname->mod; mod != NULL && rval == NULL; mod = mod->next) {
        if ((mod->subtype == ORGMOD_gb_synonym || mod->subtype == ORGMOD_old_name)
            && StringICmp (str, mod->subname) == 0) {
          rval = mod->subname;
        }
      }
    }
  }
  return rval;
}


/* Want to check that specific host names are valid */
NLM_EXTERN void 
Taxon3ValidateSpecificHostsInSeqEntry 
(SeqEntryPtr sep,
 ValNodePtr PNTR misspelled_list,
 ValNodePtr PNTR bad_caps_list,
 ValNodePtr PNTR ambiguous_list,
 ValNodePtr PNTR unrecognized_list)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   req_host_list = NULL, spec_host_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   response_vnp, request_vnp;
  ValNodePtr   check_list, check_vnp;
  OrgRefPtr    request_org, response_org;
  SpecificHostCheckPtr p;
  Boolean              has_match;
  ErrSev               level;
  Boolean              misspelled_flag;
  Boolean              bad_caps_flag;
  Boolean              ambiguous_flag;
  CharPtr              match;
    
  biop_list = GetValidatorSpecificHostBioSourceList (sep);

  /* get a list of unique specific_host values */
  spec_host_list = GetListOfUniqueSpecificHostValues (biop_list);

  /* now format requests for unique specific_host values */
  FormatValidatorSpecificHostRequests (spec_host_list, &request_list, &req_host_list);

  spec_host_list = ValNodeFreeData (spec_host_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    /* resort requests so that we can check all responses for the same BioSource together */
    check_list = SortSpecificHostOrgs (req_host_list, request_list, response_list);
    AddBioSourcesToSpecificHostChecklist (biop_list, check_list);  

    /* now look at responses */
    check_vnp = check_list;
    while (check_vnp != NULL)
    {
      p = (SpecificHostCheckPtr) check_vnp->data.ptrvalue;
      if (p != NULL)
      {
        has_match = FALSE;
        misspelled_flag = FALSE;
        bad_caps_flag = FALSE;
        ambiguous_flag = FALSE;

        request_vnp = p->request_list;
        response_vnp = p->response_list;
        while (!has_match && request_vnp != NULL && response_vnp != NULL)
        {
          request_org = (OrgRefPtr) request_vnp->data.ptrvalue;
          response_org = (OrgRefPtr) response_vnp->data.ptrvalue;
          if (response_vnp->choice == eReturnedOrgFlag_misspelled)
          {
            misspelled_flag = TRUE;
          }
          else if (response_vnp->choice == eReturnedOrgFlag_ambiguous)
          {
            ambiguous_flag = TRUE;
          }
          else
          {
            match = FindMatchInOrgRef (request_org->taxname, response_org);
            if (StringCmp (match, request_org->taxname) == 0)
            {
              has_match = TRUE;
            }
            else if (StringICmp (match, request_org->taxname) == 0)
            {
              bad_caps_flag = TRUE;
            }
          }  
          request_vnp = request_vnp->next;
          response_vnp = response_vnp->next;
        }     
        if (!has_match)
        {
          /* add to the list of bad */
          if (misspelled_flag) {
            if (misspelled_list != NULL) {
              ValNodeLink (misspelled_list, p->biop_list);
              p->biop_list = NULL;
            }
          } else if (bad_caps_flag) {
            if (bad_caps_list != NULL) {
              ValNodeLink (bad_caps_list, p->biop_list);
              p->biop_list = NULL;
            }
          } else if (ambiguous_flag) {
            if (ambiguous_list != NULL) {
              ValNodeLink (ambiguous_list, p->biop_list);
              p->biop_list = NULL;
            }
          } else {
            if (unrecognized_list != NULL) {
              ValNodeLink (unrecognized_list, p->biop_list);
              p->biop_list = NULL;
            }
          }
        }
      }
      check_vnp = check_vnp->next;
    }
    check_list = SpecificHostCheckListFree (check_list);
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);
  req_host_list = ValNodeFreeData (req_host_list);
}


typedef struct spechostgather {
  ValNodePtr list;
  Boolean    caps; /* if true, check only when first letter of first word is capitalized */
  Boolean    paren; /* if true, check portion inside parentheses as separate string */
} SpecHostGatherData, PNTR SpecHostGatherPtr;


static Boolean ShouldCheckSpecificHostString (CharPtr str, SpecHostGatherPtr p)
{
  CharPtr cp_start;
  Boolean rval = FALSE;

  if (StringHasNoText (str) || p == NULL) {
    return FALSE;
  }

  if (!p->caps) {
    rval = TRUE;
  } else if (isupper (*str)) {
    rval = TRUE;
  } else if (p->paren) {
    cp_start = StringChr (str, '(');
    if (cp_start != NULL && ShouldCheckSpecificHostString (cp_start + 1, p)) {
      rval = TRUE;
    } 
  }
  return rval;  
}
  
    
static Boolean HasSpecificHostToBeChecked (BioSourcePtr biop, SpecHostGatherPtr p)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || p == NULL) return FALSE;
  
  for (mod = biop->org->orgname->mod; mod != NULL && !rval; mod = mod->next) {
    if (mod->subtype == ORGMOD_nat_host && ShouldCheckSpecificHostString (mod->subname, p)) {
      rval = TRUE;
    }
  }
  return rval;
}


static void AddSpecificHostBioSourceFeatToList (SeqFeatPtr sfp, Pointer userdata)
{
  SpecHostGatherPtr p;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;

  p = (SpecHostGatherPtr) userdata;
  if (HasSpecificHostToBeChecked (sfp->data.value.ptrvalue, p))
  {
    ValNodeAddPointer (&(p->list), OBJ_SEQFEAT, sfp);
  }
}


static void AddSpecificHostBioSourceDescToList (SeqDescrPtr sdp, Pointer userdata)
{
  SpecHostGatherPtr p;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  p = (SpecHostGatherPtr) userdata;
  if (HasSpecificHostToBeChecked (sdp->data.ptrvalue, p))
  {
    ValNodeAddPointer (&(p->list), OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr GetSpecificHostBioSourceList (SeqEntryPtr sep, Boolean caps, Boolean paren)
{
  SpecHostGatherData   d;
  
  d.caps = caps;
  d.paren = paren;
  d.list = NULL;
  VisitFeaturesInSep (sep, &d, AddSpecificHostBioSourceFeatToList);
  VisitDescriptorsInSep (sep, &d, AddSpecificHostBioSourceDescToList);
  return d.list;
}


static void 
FormatSpecificHostRequests 
(ValNodePtr spec_host_list,
 ValNodePtr PNTR request_list,
 ValNodePtr PNTR req_host_list,
 Boolean caps,
 Boolean paren)
{
  ValNodePtr vnp;
  CharPtr    orig, cp, str, cp2 = NULL;
  
  /* now format requests for unique specific_host values */
  for (vnp = spec_host_list; vnp != NULL; vnp = vnp->next)
  {
    orig = (CharPtr) vnp->data.ptrvalue;
    /* if we have a value in parentheses, submit it separately */
    cp = StringChr (orig, '(');
    if (cp != NULL)
    {
      cp2 = StringChr (cp, ')');
    }
    if (cp != NULL && cp2 != NULL 
        && ((cp > orig && orig[StringLen (orig) - 1] == ')') /* ends with paren */
            || (cp == orig))) /* starts with paren */
    {
      if (cp > orig && orig[StringLen (orig) - 1] == ')')
      {
        str = StringSave (orig);
        /* remove trailing parenthesis */
        str [StringLen(str) - 1] = 0;

        cp = str + (cp - orig);

        /* remove opening parenthesis */
        *cp = 0;
        cp++;
      }
      else
      {
        str = StringSave (orig);
        /* remove leading parenthesis */
        str[0] = ' ';
        cp = str + (cp2 - orig);
        /* remove trailing parenthesis */
        *cp = 0; 
        cp++;
      }
      TrimSpacesAroundString (cp);
      TrimSpacesAroundString (str);
      if (paren && (!caps || isupper (*cp))) {
        AddRequestOrgForString (cp, orig, request_list, req_host_list);
      }
      if (!caps || isupper (*str)) {
        AddRequestOrgForString (str, orig, request_list, req_host_list);
      }
    }
    else
    {
      if (!caps || isupper (*orig)) {
        AddRequestOrgForString (orig, orig, request_list, req_host_list);
      }
    }
  }
}


typedef struct replacementpair {
  CharPtr find;
  CharPtr repl;
} ReplacementPairData, PNTR ReplacementPairPtr;

static ReplacementPairPtr ReplacementPairNew (CharPtr find, CharPtr repl)
{
  ReplacementPairPtr r;

  r = (ReplacementPairPtr) MemNew (sizeof (ReplacementPairData));
  r->find = StringSave (find);
  r->repl = StringSave (repl);
  return r;
}

static ReplacementPairPtr ReplacementPairFree (ReplacementPairPtr r)
{
  if (r != NULL) {
    r->find = MemFree (r->find);
    r->repl = MemFree (r->repl);
    r = MemFree (r);
  }
  return r;
}

static ValNodePtr ReplacementPairListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = ReplacementPairFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static SpecificHostFixPtr 
SpecificHostFixNew 
(ValNodePtr feat_or_desc,
 CharPtr    bad_host,
 CharPtr    old_taxname,
 CharPtr    new_taxname,
 Uint1      fix_type)
{
  SpecificHostFixPtr s;

  s = (SpecificHostFixPtr) MemNew (sizeof (SpecificHostFixData));
  if (feat_or_desc != NULL) 
  {
    s->feat_or_desc = ValNodeNew(NULL);
    s->feat_or_desc->choice = feat_or_desc->choice;
    s->feat_or_desc->data.ptrvalue = feat_or_desc->data.ptrvalue;
  }
  s->bad_specific_host = StringSave (bad_host);
  s->old_taxname = StringSave (old_taxname);
  s->new_taxname = StringSave (new_taxname);
  s->fix_type = fix_type;
  return s;
}


static SpecificHostFixPtr SpecificHostFixFree (SpecificHostFixPtr s)
{
  if (s != NULL)
  {
    s->feat_or_desc = ValNodeFree (s->feat_or_desc);
    s->bad_specific_host = MemFree (s->bad_specific_host);
    s->old_taxname = MemFree (s->old_taxname);
    s->new_taxname = MemFree (s->new_taxname);
    s = MemFree (s);
  }
  return s;
}


extern ValNodePtr SpecificHostFixListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = SpecificHostFixFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static ValNodePtr GetFixesForOneSpecificHostValue (SpecificHostCheckPtr p)
{
  CharPtr      prev_success = NULL, new_val, prev_fail = NULL;
  Boolean      fix_needed = FALSE;
  ValNodePtr   suggested_fixes = NULL;
  OrgRefPtr    request_org, response_org;
  ValNodePtr   biop_vnp, response_vnp, request_vnp, vnp;
  SpecificHostFixPtr s;
  ValNodePtr         fix_list = NULL;
  ReplacementPairPtr r;
  Uint1              fix_type;
  Boolean            add_nontrunc_fix;
  Boolean            ambiguous = FALSE;

  if (p == NULL) return NULL;

  request_vnp = p->request_list;
  response_vnp = p->response_list;
  
  while (request_vnp != NULL && response_vnp != NULL)
  {
    request_org = (OrgRefPtr) request_vnp->data.ptrvalue;
    response_org = (OrgRefPtr) response_vnp->data.ptrvalue;
    if (prev_success != NULL 
        && StringNCmp (request_org->taxname, prev_success, StringLen (request_org->taxname)) == 0) {
      /* we don't need to check this one */
    } else if (response_org == NULL) {
      fix_needed = TRUE;
      if (response_vnp->choice == eReturnedOrgFlag_ambiguous) {
        ambiguous = TRUE;
      }
      if (prev_fail == NULL) {
        prev_fail = request_org->taxname;
      } else if (StringNCmp (prev_fail, request_org->taxname, StringLen (request_org->taxname)) != 0) {
        if (response_vnp->choice == eReturnedOrgFlag_ambiguous) {
          ValNodeAddPointer (&suggested_fixes, eSpecificHostFix_ambiguous, ReplacementPairNew (request_org->taxname, NULL));
        } else {
          ValNodeAddPointer (&suggested_fixes, eSpecificHostFix_unrecognized, ReplacementPairNew (request_org->taxname, NULL));
        }
        prev_fail = request_org->taxname;      
      }
    } else {
      prev_success = request_org->taxname;
      add_nontrunc_fix = FALSE;
      if (response_vnp->choice == eReturnedOrgFlag_misspelled) {
        fix_needed = TRUE;
        fix_type = eSpecificHostFix_spelling;
        new_val = response_org->taxname;
        add_nontrunc_fix = TRUE;
      } else {
        new_val = FindMatchInOrgRef (request_org->taxname, response_org);
        if (new_val == NULL) {
          fix_needed = TRUE;
          fix_type = eSpecificHostFix_replacement;
          new_val = response_org->taxname;
          add_nontrunc_fix = TRUE;
        } else if (StringCmp (new_val, request_org->taxname) != 0) {
          fix_needed = TRUE;
          fix_type = eSpecificHostFix_capitalization;
          add_nontrunc_fix = TRUE;
        }
      }

      /* add fix to truncate and correct spelling and capitalization first */
      /* this way the truncation won't fail when it looks for the old version that's already been corrected */
      if (prev_fail != NULL) {
        if (StringNCmp (prev_fail, request_org->taxname, StringLen (request_org->taxname)) == 0) {
          if (new_val != NULL) {
            ValNodeAddPointer (&suggested_fixes, eSpecificHostFix_truncation, ReplacementPairNew (prev_fail, new_val));            
            fix_needed = TRUE;
          }
        } else {
          ValNodeAddPointer (&suggested_fixes, eSpecificHostFix_unrecognized, ReplacementPairNew (prev_fail, NULL));
        }
      }
      /* add fix for just spelling and capitalization after */
      if (add_nontrunc_fix) {
        ValNodeAddPointer (&suggested_fixes, fix_type, ReplacementPairNew (request_org->taxname, new_val));
      }

      prev_fail = NULL;
    }
    request_vnp = request_vnp->next;
    response_vnp = response_vnp->next;
  }

  if (fix_needed) {
    for (biop_vnp = p->biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next) {
      if (suggested_fixes == NULL) {
        s = SpecificHostFixNew (biop_vnp, p->spec_host, p->spec_host, NULL, ambiguous ? eSpecificHostFix_ambiguous : eSpecificHostFix_unrecognized);
        ValNodeAddPointer (&fix_list, 0, s);
      } else {
        for (vnp = suggested_fixes; vnp != NULL; vnp = vnp->next) {
          r = (ReplacementPairPtr) vnp->data.ptrvalue;
          s = SpecificHostFixNew (biop_vnp, p->spec_host, r->find, r->repl, vnp->choice);
          ValNodeAddPointer (&fix_list, 0, s);
        }
      }
    }
  }
  suggested_fixes = ReplacementPairListFree (suggested_fixes);
  return fix_list;
}


NLM_EXTERN ValNodePtr Taxon3GetSpecificHostFixesInSeqEntry (SeqEntryPtr sep, Boolean caps, Boolean paren)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   req_host_list = NULL, spec_host_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   check_list, check_vnp;
  SpecificHostCheckPtr p;
  ErrSev               level;
  ValNodePtr           fix_list = NULL;
  
  biop_list = GetSpecificHostBioSourceList (sep, caps, paren);

  /* get a list of unique specific_host values */
  spec_host_list = GetListOfUniqueSpecificHostValues (biop_list);

  /* now format requests for unique specific_host values */
  FormatSpecificHostRequests (spec_host_list, &request_list, &req_host_list, caps, paren);

  spec_host_list = ValNodeFreeData (spec_host_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    /* resort requests so that we can check all responses for the same BioSource together */
    check_list = SortSpecificHostOrgs (req_host_list, request_list, response_list);
    AddBioSourcesToSpecificHostChecklist (biop_list, check_list);  

    /* now look at responses */
    check_vnp = check_list;
    while (check_vnp != NULL)
    {
      p = (SpecificHostCheckPtr) check_vnp->data.ptrvalue;
      ValNodeLink (&fix_list, GetFixesForOneSpecificHostValue (p));
      check_vnp = check_vnp->next;
    }
    check_list = SpecificHostCheckListFree (check_list);
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);
  req_host_list = ValNodeFreeData (req_host_list);

  return fix_list;
}


extern Boolean ApplyOneSpecificHostFix (SpecificHostFixPtr s)
{
  BioSourcePtr biop = NULL;
  Boolean      rval = FALSE;
  CharPtr      new_spec_host;
  ValNode      vn;

  if (s == NULL || s->feat_or_desc == NULL 
      || StringHasNoText (s->bad_specific_host)
      || StringHasNoText (s->old_taxname)
      || StringHasNoText (s->new_taxname)) {
    return rval;
  }
  biop = GetBioSourceFromValNode (s->feat_or_desc);
  if (biop == NULL) return rval;

  vn.choice = SourceQualChoice_textqual;
  vn.data.intvalue = Source_qual_nat_host;
  vn.next = NULL;

  new_spec_host = GetSourceQualFromBioSource (biop, &vn, NULL);  
  FindReplaceString (&new_spec_host, s->old_taxname, s->new_taxname, TRUE, TRUE);
  if (StringCmp (new_spec_host, s->bad_specific_host) != 0)
  {
    rval = SetSourceQualInBioSource (biop, &vn, NULL, new_spec_host, ExistingTextOption_replace_old);
  }
  new_spec_host = MemFree (new_spec_host);
  return rval;
}

static void AddBioSourceFeatToList (SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
}


static void AddBioSourceDescToList (SeqDescrPtr sdp, Pointer userdata)
{

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQDESC, sdp);
}


static ValNodePtr GetBioSourceList (SeqEntryPtr sep)
{
  ValNodePtr list = NULL;
  
  VisitFeaturesInSep (sep, &list, AddBioSourceFeatToList);
  VisitDescriptorsInSep (sep, &list, AddBioSourceDescToList);
  return list;
}


static ValNodePtr GetListOfOrganismNames (ValNodePtr biop_list)
{
  ValNodePtr   biop_vnp;
  BioSourcePtr biop;
  ValNodePtr   list = NULL;
  
  /* get a list of unique specific_host values */
  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {
    if (biop_vnp->data.ptrvalue == NULL) continue;
    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)) continue;
    if (!StringAlreadyInValNodeList (biop->org->taxname, list))
    {
      ValNodeAddPointer (&list, 0, biop->org->taxname);
    }
  }
  return list;
}


static void AddBioSourcesToChecklist (ValNodePtr biop_list, ValNodePtr check_list)
{
  ValNodePtr biop_vnp, last_vnp = NULL, stop_search;
  BioSourcePtr biop;
  SpecificHostCheckPtr p;

  if (biop_list == NULL || check_list == NULL) return;

  for (biop_vnp = biop_list; biop_vnp != NULL; biop_vnp = biop_vnp->next)
  {

    biop = GetBioSourceFromValNode (biop_vnp);
    if (biop == NULL) continue;

    if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) continue;
    if (last_vnp == NULL)
    {
      last_vnp = check_list;
      stop_search = NULL;
    }
    else
    {
      stop_search = last_vnp;
    }
    p = NULL;
    while (last_vnp != NULL 
           && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
           && StringCmp (p->spec_host, biop->org->taxname) != 0)
    {
      p = NULL;
      last_vnp = last_vnp->next;
    }
    if (p == NULL && stop_search != NULL)
    {
      last_vnp = check_list;
      while (last_vnp != stop_search 
              && (p = (SpecificHostCheckPtr) last_vnp->data.ptrvalue) != NULL
              && StringCmp (p->spec_host, biop->org->taxname) != 0)
      {
        p = NULL;
        last_vnp = last_vnp->next;
      }
    }

    if (p != NULL)
    {
      ValNodeAddPointer (&(p->biop_list), biop_vnp->choice, biop_vnp->data.ptrvalue);
    }
  }
}


static ValNodePtr GetBioSourcesWithTaxName (CharPtr taxname, ValNodePtr biop_list)
{
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;
  ValNodePtr match_list = NULL, vnp;

  if (StringHasNoText (taxname) || biop_list == NULL) return NULL;

  for (vnp = biop_list; vnp != NULL; vnp = vnp->next) {
    biop = NULL;
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp != NULL && sdp->choice == Seq_descr_source) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      }
    }
    if (biop != NULL && biop->org != NULL && StringCmp (taxname, biop->org->taxname) == 0) {
      ValNodeAddPointer (&match_list, vnp->choice, vnp->data.ptrvalue);
    }
  }
  return match_list;
}


NLM_EXTERN ValNodePtr GetOrganismTaxLookupFailuresInSeqEntry (SeqEntryPtr sep)
{
  ValNodePtr   biop_list = NULL;
  ValNodePtr   unique_list = NULL;
  ValNodePtr   request_list = NULL;
  ValNodePtr   response_list = NULL;
  ValNodePtr   req_vnp, resp_vnp;
  ErrSev               level;
  ValNodePtr           failed_list = NULL, vnp;
  OrgRefPtr            request_org;
  
  biop_list = GetBioSourceList (sep);

  /* get a list of unique specific_host values */
  unique_list = GetListOfOrganismNames (biop_list);

  /* now format requests for unique taxname values */
  for (vnp = unique_list; vnp != NULL; vnp = vnp->next) 
  {
    request_org = OrgRefNew();
    request_org->taxname = StringSave (vnp->data.ptrvalue);
    ValNodeAddPointer (&request_list, 3, request_org);
  }

  unique_list = ValNodeFree (unique_list);

  level = ErrSetMessageLevel (SEV_MAX);
  response_list = Taxon3GetOrgRefList (request_list);
  ErrSetMessageLevel (level);
 
  if (ValNodeLen (response_list) != ValNodeLen (request_list))
  {
    Message (MSG_POST, "Unable to retrieve information from tax server");
  }
  else
  {
    for (req_vnp = request_list, resp_vnp = response_list;
         req_vnp != NULL && resp_vnp != NULL;
         req_vnp = req_vnp->next, resp_vnp = resp_vnp->next)
    {
      if (resp_vnp->data.ptrvalue == NULL)
      {        
        request_org = (OrgRefPtr) req_vnp->data.ptrvalue;
        vnp = GetBioSourcesWithTaxName (request_org->taxname, biop_list);
        if (vnp != NULL) {
          ValNodeAddPointer (&failed_list, 0, StringSave (request_org->taxname));
          ValNodeLink (&failed_list, vnp);
        }
      }
    }
  }

  biop_list = ValNodeFree (biop_list);
  request_list = FreeOrgRefValNodeList (request_list);
  response_list = FreeOrgRefValNodeList (response_list);

  return failed_list;  
}

