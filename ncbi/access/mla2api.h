/*   mla2api.h
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
* File Name:  mla2api.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/30/07
*
* $Revision: 1.7 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#ifndef _MLA2API_
#define _MLA2API_

#include <ncbi.h>
#include <asn.h>
#include <mlkludge.h>
#include <objmla2.h>
#include <urlquery.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* low-level connection functions */

NLM_EXTERN CONN Mla2OpenConnection (void);

NLM_EXTERN MlaBackPtr Mla2WaitForReply (
  CONN conn
);

/*
 Mla2SynchronousQuery opens connection, sends
 MlaRequestPtr ASN.1 query, and waits for reply,
 cleaning up connection afterwards.
*/

NLM_EXTERN MlaBackPtr Mla2SynchronousQuery (
  MlaRequestPtr mrp
);

/*
 Mla2AsynchronousQuery opens connection, sends
 request, and queues completion routine using urlquery
 queueing mechanism.

 Mla2CheckQueue should be called several times a
 second with a timer.  It calls QUERY_CheckQueue to
 poll connection, which calls completion routine when
 result is available, cleaning up connection afterwards.

 Mla2ReadReply takes conn and status parameters from
 completion routine and reads MlaBackPtr.
*/

NLM_EXTERN Boolean Mla2AsynchronousQuery (
  MlaRequestPtr mrp,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
);

NLM_EXTERN Int4 Mla2CheckQueue (
  QUEUE* queue
);

NLM_EXTERN MlaBackPtr Mla2ReadReply (
  CONN conn,
  EIO_Status status
);

/* request creation functions */

NLM_EXTERN MlaRequestPtr Mla2CreateJournalTitleRequest (
  CharPtr journal_name
);

NLM_EXTERN MlaRequestPtr Mla2CreateCitArtJournalRequest (
  CitArtPtr cap
);

NLM_EXTERN MlaRequestPtr Mla2CreateCitArtMatchRequest (
  CitArtPtr cap
);

NLM_EXTERN MlaRequestPtr Mla2CreateCitationtMatchRequest (
  CharPtr author,
  CharPtr journal,
  CharPtr volume,
  CharPtr page,
  Int4 year,
  CharPtr title
);

NLM_EXTERN MlaRequestPtr Mla2CreatePubFetchRequest (
  Int4 pmid
);

/*
   reply extraction functions - these do not free the MlaBackPtr,
   but unlink the internal structure so MlaBackFree can be called
*/

NLM_EXTERN TitleMsgListPtr Mla2ExtractJournalTitleReply (
  MlaBackPtr mbp
);

NLM_EXTERN Int4 Mla2ExtractCitMatchReply (
  MlaBackPtr mbp
);

NLM_EXTERN CitArtPtr Mla2ExtractPubFetchReply (
  MlaBackPtr mbp
);

/* utility functions */

NLM_EXTERN void ChangeCitArtMLAuthorsToSTD (
  CitArtPtr cap
);

NLM_EXTERN void ChangeMlaBackMLAuthorsToSTD (
  MlaBackPtr mbp 
);

NLM_EXTERN Boolean Mla2IsEPubOnlyJournal (
  CharPtr jta,
  Int2Ptr starting_yearP
);

/*
   Mla2CorrectCitArt takes the original and fetched CitArt and
   fixes author types, compares author lists for differences,
   restores consortium authors, and other automatic cleanup
   functions to be determined with experience
*/

NLM_EXTERN Int2 Mla2CorrectCitArt (
  CitArtPtr newcap,
  CitArtPtr oldcap
);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _MLA2API_ */

