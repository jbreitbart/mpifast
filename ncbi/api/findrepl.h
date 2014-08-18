/*
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
* File Name: findrepl.h
*
* Author:  Jonathan Kans, Tim Ford
*
* Version Creation Date:   10/17/00
*
* File Description:
*   Complete redesign of find/replace from original of Yuri Sadykov
*
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
*/

#ifndef __FINDREPL_H__
#define __FINDREPL_H__

#include <ncbi.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */


/*** defines for send_update parameter below ***/

#define UPDATE_NEVER 0  /* never send ObjMgrMessage (OM_MSG_UPDATE... */
#define UPDATE_EACH  1  /* send it on each replace */
#define UPDATE_ONCE  2  /* send once for whole entityID, if any replacements occur */

typedef void (*FindReplProc) (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Pointer userdata);

typedef void (*StringActionFunc) (CharPtr PNTR strp, Pointer userdata, BoolPtr did_find, BoolPtr did_change);
/* EXTERNAL FIND-REPLACE FUNCTIONS */
NLM_EXTERN void StringActionInEntity (
  Uint2 entityID,
  Boolean select_item,
  Int2 send_update,
  BoolPtr descFilter,
  BoolPtr featFilter,
  BoolPtr seqidFilter,
  Boolean do_seqid_local,
  StringActionFunc action_func,
  FindReplProc callback,
  Pointer userdata
);

NLM_EXTERN void StringActionForObject (
  Uint2   datatype,
  Pointer objdata,
  Uint2 entityID,
  Boolean select_item,
  Int2 send_update,
  StringActionFunc action_func,
  FindReplProc callback,
  Pointer userdata
);

NLM_EXTERN void SpecialCharReplace (CharPtr PNTR strp, Pointer userdata, BoolPtr did_find, BoolPtr did_change);
NLM_EXTERN void SpecialCharFind (CharPtr PNTR strp, Pointer userdata, BoolPtr did_find, BoolPtr did_change);
NLM_EXTERN CharPtr GetSpecialCharacterReplacement (unsigned char ch);
NLM_EXTERN CharPtr GetSpecialWinCharacterReplacement (unsigned char ch);
NLM_EXTERN CharPtr GetSpecialMacCharacterReplacement (unsigned char ch);

NLM_EXTERN void FindReplaceInEntity (
  Uint2 entityID,
  CharPtr find_string,
  CharPtr replace_string,
  Boolean case_counts,
  Boolean whole_word,
  Boolean do_replace,
  Boolean select_item,
  Int2 send_update,
  BoolPtr descFilter,
  BoolPtr featFilter,
  BoolPtr seqidFilter,
  Boolean do_seqid_local,
  FindReplProc callback,
  Pointer userdata
);

NLM_EXTERN void FindReplaceString (
  CharPtr PNTR strp,
  CharPtr find_string,
  CharPtr replace_string,
  Boolean case_counts,
  Boolean whole_word
);

NLM_EXTERN void FindStringsInEntity (
  Uint2 entityID,
  CharPtr PNTR find_strings,
  Boolean case_counts,
  Boolean whole_word,
  Boolean select_item,
  Int2 send_update,
  BoolPtr descFilter,
  BoolPtr featFilter,
  BoolPtr seqidFilter,
  Boolean do_seqid_local,
  FindReplProc callback,
  Pointer userdata
);


extern void RemoveTaxRef (OrgRefPtr orp);


#ifdef __cplusplus
extern "C" }
#endif	/* __cplusplus */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif	/* __FINDREPL_H__ */
