/*   ncbifile.c
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
* File Name:  ncbifile.c
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   3/4/91
*
* $Revision: 6.41 $
*
* File Description: 
*     portable file routines
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#define THIS_MODULE g_corelib
#define THIS_FILE  _this_file

#include <ncbilcl.h>
#include <assert.h>

#include "corepriv.h"
#ifdef MPIBLAST_PIO
#include "virtual_files.hpp"
#endif

#ifdef OS_MAC
#ifdef OS_UNIX_DARWIN
#include <Carbon.h>
#endif
#include <Navigation.h>
#include <Script.h>
#endif

#ifdef OS_UNIX_SUN
#define DEFAULT_CDROM "/dev/sr0"
#define DEFAULT_RAW_CDROM "/dev/rsr0"
#endif

#if defined(PROC_MIPS) && !defined(OS_UNIX_LINUX)
#define DEFAULT_CDROM "/dev/scsi/sc0d4l0"
#endif

#ifdef OS_UNIX
#ifndef DEFAULT_CDROM
#define DEFAULT_CDROM "/dev/cdrom"
#endif
#endif

#if defined(OS_MSWIN) || defined (OS_NT)
#ifdef COMP_MSC
#ifndef mkdir
#define mkdir _mkdir
#endif
#ifndef stat
#define stat _stat
#endif
#endif
#endif

#ifdef OS_VMS
#ifndef DEFAULT_CDROM
#define DEFAULT_CDROM "cdrom:"
#endif
#endif

extern char *g_corelib;
static char * _this_file = __FILE__;

#ifdef OS_MAC
#define INLINE_MOREFILES
/* #include "FullPath.h" */
/* #include "MoreFilesExtra.h" */
#ifdef INLINE_MOREFILES

/* MoreFilesExtras.c */
/* ----------------- */

/*
**	Apple Macintosh Developer Technical Support
**
**	A collection of useful high-level File Manager routines.
**
**	by Jim Luther, Apple Developer Technical Support Emeritus
**
**	File:		MoreFilesExtras.c
**
**	Copyright © 1992-1999 Apple Computer, Inc.
**	All rights reserved.
**
**	You may incorporate this sample code into your applications without
**	restriction, though the sample code has been provided "AS IS" and the
**	responsibility for its operation is 100% yours.  However, what you are
**	not permitted to do is to redistribute the source as "DSC Sample Code"
**	after having made changes. If you're going to re-distribute the source,
**	we require that you make it clear in the source that the code was
**	descended from Apple Sample Code, but that you've made changes.
*/

#include <Types.h>
#include <Traps.h>
#include <OSUtils.h>
#include <Errors.h>
#include <Files.h>
#include <Devices.h>
#include <Finder.h>
#include <Folders.h>
#include <FSM.h>
#include <Disks.h>
#include <Gestalt.h>
#include <TextUtils.h>
#include <Script.h>
#include <Math64.h>
#include <CodeFragments.h>
#include <stddef.h>

#define	__COMPILINGMOREFILES

#if 0
#include "MoreFiles.h"
#include "MoreFilesExtras.h"
#include "MoreDesktopMgr.h"
#include "FSpCompat.h"
#endif

/*
**	GetVolumeInfoNoName uses pathname and vRefNum to call PBHGetVInfoSync
**	in cases where the returned volume name is not needed by the caller.
**	The pathname and vRefNum parameters are not touched, and the pb
**	parameter is initialized by PBHGetVInfoSync except that ioNamePtr in
**	the parameter block is always returned as NULL (since it might point
**	to the local tempPathname).
**
**	I noticed using this code in several places, so here it is once.
**	This reduces the code size of MoreFiles.
*/
static pascal	OSErr	GetVolumeInfoNoName(ConstStr255Param pathname,
									short vRefNum,
									HParmBlkPtr pb)
{
	Str255 tempPathname;
	OSErr error;
	
	/* Make sure pb parameter is not NULL */ 
	if ( pb != NULL )
	{
		pb->volumeParam.ioVRefNum = vRefNum;
		if ( pathname == NULL )
		{
			pb->volumeParam.ioNamePtr = NULL;
			pb->volumeParam.ioVolIndex = 0;		/* use ioVRefNum only */
		}
		else
		{
			BlockMoveData(pathname, tempPathname, pathname[0] + 1);	/* make a copy of the string and */
			pb->volumeParam.ioNamePtr = (StringPtr)tempPathname;	/* use the copy so original isn't trashed */
			pb->volumeParam.ioVolIndex = -1;	/* use ioNamePtr/ioVRefNum combination */
		}
		error = PBHGetVInfoSync(pb);
		pb->volumeParam.ioNamePtr = NULL;	/* ioNamePtr may point to local	tempPathname, so don't return it */
	}
	else
	{
		error = paramErr;
	}
	return ( error );
}

/*****************************************************************************/

static pascal	OSErr	DetermineVRefNum(ConstStr255Param pathname,
								 short vRefNum,
								 short *realVRefNum)
{
	HParamBlockRec pb;
	OSErr error;

	error = GetVolumeInfoNoName(pathname,vRefNum, &pb);
	if ( error == noErr )
	{
		*realVRefNum = pb.volumeParam.ioVRefNum;
	}
	return ( error );
}

/*****************************************************************************/

static pascal	OSErr GetCatInfoNoName(short vRefNum,
							   long dirID,
							   ConstStr255Param name,
							   CInfoPBPtr pb)
{
	Str31 tempName;
	OSErr error;
	
	/* Protection against File Sharing problem */
	if ( (name == NULL) || (name[0] == 0) )
	{
		tempName[0] = 0;
		pb->dirInfo.ioNamePtr = tempName;
		pb->dirInfo.ioFDirIndex = -1;	/* use ioDirID */
	}
	else
	{
		pb->dirInfo.ioNamePtr = (StringPtr)name;
		pb->dirInfo.ioFDirIndex = 0;	/* use ioNamePtr and ioDirID */
	}
	pb->dirInfo.ioVRefNum = vRefNum;
	pb->dirInfo.ioDrDirID = dirID;
	error = PBGetCatInfoSync(pb);
	pb->dirInfo.ioNamePtr = NULL;
	return ( error );
}

/*****************************************************************************/

static pascal	OSErr	GetDirectoryID(short vRefNum,
							   long dirID,
							   ConstStr255Param name,
							   long *theDirID,
							   Boolean *isDirectory)
{
	CInfoPBRec pb;
	OSErr error;

	error = GetCatInfoNoName(vRefNum, dirID, name, &pb);
	if ( error == noErr )
	{
		*isDirectory = (pb.hFileInfo.ioFlAttrib & kioFlAttribDirMask) != 0;
		if ( *isDirectory )
		{
			*theDirID = pb.dirInfo.ioDrDirID;
		}
		else
		{
			*theDirID = pb.hFileInfo.ioFlParID;
		}
	}
	
	return ( error );
}

/*****************************************************************************/

static pascal	OSErr	FSpGetDirectoryID(const FSSpec *spec,
								  long *theDirID,
								  Boolean *isDirectory)
{
	return ( GetDirectoryID(spec->vRefNum, spec->parID, spec->name,
			 theDirID, isDirectory) );
}

#endif  /* INLINE_MOREFILES */
#endif  /* OS_MAC */

/*****************************************************************************
*
*   Macintosh file utilities
*
*****************************************************************************/

#ifdef OS_MAC

static OSErr MacPathname2FSSpec(const char *inPathname, FSSpec *outFSS)
{
	OSErr err;
	size_t len;
	char *p;
	short vRefNum;
	long dirID;
	FSSpec fss;
	
	if (inPathname == NULL || outFSS == NULL) {
		return paramErr;
	}
	
	err = HGetVol(NULL, &vRefNum, &dirID);  /* default volume and directory */
	if (err != noErr) return err;
	
	len = strlen(inPathname);
	
	p = strchr(inPathname, ':');
	if (p == NULL) {
		/* Partial pathname -- filename only */
		Str31 filename;
		assert(len <= 31);
		c2pstrcpy(filename, inPathname);
		err = FSMakeFSSpec(vRefNum, dirID, filename, outFSS);
	} else {
		Str31 name;
		int nameLen;
		if (inPathname[0] == ':') {
			/* Relative pathname including directory path */
			
		} else {
			/* Absolute pathname */
			/* Str31 volName;  We would use Str28 if it was defined -- 27, plus 1 for ':'. */
			nameLen = p - inPathname;
			assert(nameLen <= 27);
			name[0] = nameLen + 1;
			memcpy(name + 1, inPathname, nameLen + 1);  /* Copy the volume name and the colon. */
			err = DetermineVRefNum(name, 0, &vRefNum);
			if (err != noErr) return err;
			dirID = 2;
		}
		/* vRefNum and dirID now specify the directory in which we should descend
		   the path pointed to by p (pointing to the first colon). */
		p++;
		while (p != NULL && *p != '\0') {
			char *q = strchr(p, ':');
			if (q != NULL) {
				Boolean isDir;
				nameLen = q - p;
				assert(nameLen <= 31);
				name[0] = nameLen;
				memcpy(name + 1, p, nameLen);
				err = FSMakeFSSpec(vRefNum, dirID, name, &fss);
				if (err != noErr) return err;
				if (q[1] == '\0') {
					p = NULL;
					*outFSS = fss;
				} else {
					err = FSpGetDirectoryID(&fss, &dirID, &isDir);
					assert(isDir == true);
					if (err != noErr) return err;
					p = q + 1;
				}
			} else {
				q = strchr(p, '\0');  /* go to end of string */
				nameLen = q - p;
				assert(nameLen > 0);
				assert(nameLen <= 31);
				c2pstrcpy(name, p);
				p = NULL;
				err = FSMakeFSSpec(vRefNum, dirID, name, outFSS);
			}
		}
	}
	return err;
}

#if 0
static OSErr MacFSSpec2FullPathname(const FSSpec *inFSS, char **outPathname)
{
	OSErr err;
	Handle h;
	short fullPathLength;
	static char *fullPath = NULL;
	
	if (fullPath != NULL) {
		Nlm_Free(fullPath);
		fullPath = NULL;
	}
	err = FSpGetFullPath(inFSS, &fullPathLength, &h);
	if (err != noErr) return err;
	
	assert(fullPathLength >= 2);  /* An absolute pathname must be at least two chars long */
	fullPath = (char *)Nlm_Malloc(fullPathLength + 1);
	if (fullPath == NULL) {
		err = memFullErr;
	} else {
		strncpy(fullPath, *h, fullPathLength);
	}
	
	DisposeHandle(h);
	
	*outPathname = fullPath;
	return err;
}
#endif

static OSErr MacCreateDirectory(const char *inPathname)
{
	OSErr err;
	FSSpec fss;
	ScriptCode scriptTag = 0;
	long createdDirID;
	
	err = MacPathname2FSSpec(inPathname, &fss);
	if (err != noErr) return err;
	
	err = FSpDirCreate(&fss, scriptTag, &createdDirID);
	return err;
}

#endif

/*****************************************************************************
*
*   FileOpen(filename, mode)
*     if (filename == "stdin" or "stdout" or "stderr"
*           returns those predefined
*           streams on non-windowing systems)
*
*****************************************************************************/

/* p_churchill 12/99 removed MPW conditional compilation
 */

static Nlm_FileOpenHook _hookFile = NULL;


NLM_EXTERN FILE * LIBCALL Nlm_FileOpen(const char *filename, const char *mode)
{
  FILE *f = NULL;

  if ( _hookFile )
    return _hookFile(filename, mode);

#if defined(WIN_DUMB)
    if ( Nlm_HasConsole )
      {
        if      ( !StringCmp("stdin",  filename))
          f = stdin;
        else if ( !StringCmp("stdout", filename) )
          f = stdout;
        else if ( !StringCmp("stderr", filename))
          f = stderr;
        else
          f = fopen(filename, mode);

#if defined(WIN32)  &&  ! defined(COMP_METRO) 
        if (strchr(mode, 'b')  &&
            (f == stdin  ||  f == stdout  ||  f == stderr))
          setmode(fileno(f), O_BINARY);
#endif
      }
    else
      f = fopen(filename, mode);

#elif defined(OS_MAC)
  {
    OSType    fCreator;
    Nlm_Int2  fError;
    FInfo     fInfo;
    OSType    fType;
    Nlm_Char  temp [256];

    Nlm_StringNCpy_0(temp, filename, sizeof(temp));
    Nlm_CtoPstr ((Nlm_CharPtr) temp);
    fError = HGetFInfo( 0, 0, (StringPtr) temp, &fInfo);
    if (fError == noErr) {
      fCreator = fInfo.fdCreator;
      fType = fInfo.fdType;
    } else {
        fCreator = '    ';
        if (strchr(mode, 'b') != NULL)
            fType = '    ';
        else
            fType = 'TEXT';
    }
    f = fopen( filename, mode);

    fError = HGetFInfo( 0, 0, (StringPtr) temp, &fInfo);
    if (fError == noErr) {
      fInfo.fdCreator = fCreator;
      fInfo.fdType = fType;
      fError = HSetFInfo ( 0, 0, (StringPtr) temp,&fInfo);
    }
  } /* def OS_MAC */

#elif defined(OS_VMS) && defined(DCC4DW12)
  /* never used */ 
  f = fopen (filename, mode);
  if (f  &&
      fstat(fileno(f), &statbuf) == 0  &&
      statbuf.st_fab_rfm == FAB$C_UDF)
    {
      fclose(f);
      f = fopen(filename,mode,"ctx=stm");
    }

#else
    f = fopen(filename, mode);  
#endif

  if (f == NULL)
    ErrPostEx(SEV_INFO, E_File, E_FOpen, "FileOpen(\"%s\",\"%s\") failed",
              filename, mode);
		
  return f;
}

/*****************************************************************************
*
*   SetFileOpenHook(hook)
*
*****************************************************************************/

NLM_EXTERN void LIBCALL Nlm_SetFileOpenHook (Nlm_FileOpenHook hook)
{
	_hookFile = hook;
}

/*****************************************************************************
*
*   FileClose(fp)
*
*****************************************************************************/

NLM_EXTERN void LIBCALL  Nlm_FileClose (FILE *stream)
{
  if (stream == NULL)
    return;

#ifdef WIN_DUMB    
  if (stream == stdin  ||  stream == stdout  ||  stream == stderr)
    {
#if defined(WIN32)  &&  ! defined(COMP_METRO) 
      setmode(fileno(stream), O_TEXT);
#endif
      return;
    }
#endif

  fclose(stream);
}

/*****************************************************************************
*   FileRead(buf, size, fp)
*****************************************************************************/
NLM_EXTERN size_t LIBCALL Nlm_FileRead
(void *ptr, size_t size, size_t n, FILE *stream)
{
  if (n  &&  (SIZE_MAX / n) < size) {
    ErrPostEx(SEV_WARNING,E_Programmer,0,"FileRead: size > SIZE_MAX");
    return 0;
  }
  if (!ptr  ||  !stream)
    return 0;

  return fread(ptr,size,n,stream);
}

/*****************************************************************************
*   FileWrite(buf, size, fp)
*****************************************************************************/
NLM_EXTERN size_t LIBCALL Nlm_FileWrite
(const void *ptr, size_t size, size_t n, FILE *stream)
{
    size_t cnt;
    if (n   &&  (SIZE_MAX / n)  <  size) {
        ErrPostEx(SEV_WARNING,E_Programmer,0,"FileWrite:  size > SIZE_MAX");
        return 0;
    }
    if (!ptr  ||  !stream || !size)
        return 0;
    
    cnt = fwrite(ptr,size,n,stream);
    if (cnt != n)
        ErrPostEx(SEV_FATAL,E_File,E_FWrite,"File write error");
    
    return cnt;
}

/*****************************************************************************
*
*   FilePuts(ptr, fp)
*
*****************************************************************************/
NLM_EXTERN int LIBCALL  Nlm_FilePuts (const char *ptr, FILE *fp)
{
	int retval;

	if ((ptr == NULL) || (fp == NULL))
    	return EOF;
	if ((retval = fputs(ptr,fp)) ==EOF)
		ErrPostEx(SEV_FATAL,E_File,E_FWrite,"File write error");
	return retval;
}

/*****************************************************************************
*
*   FileGets()
*
*****************************************************************************/
NLM_EXTERN char * LIBCALL  Nlm_FileGets (Nlm_CharPtr ptr, size_t size, FILE *fp)
{
#if defined(OS_MAC) || defined (OS_UNIX_DARWIN)
    int         ch;
	int         count;
	Nlm_CharPtr tmp;
#endif

	if ((ptr == NULL) || (size < 1) || (fp == NULL))
		return NULL;
#if defined(OS_MAC) || defined (OS_UNIX_DARWIN)
	ch = fgetc (fp);
	count = 0;
	tmp = ptr;
	while (ch != EOF && ch != '\0' && ch != '\n' && ch != '\r' && count < size - 2) {
	  *tmp = ch;
	  tmp++;
	  count++;
	  ch = fgetc (fp);
	}
	if (ch == '\n' || ch == '\r') {
	  *tmp = '\n';
	  tmp++;
	  count++;
	} else if (ch != EOF && ch != '\0') {
	  *tmp = ch;
	  tmp++;
	  count++;
	}
	*tmp = '\0';
	if (count < 1)
		return NULL;
	return ptr;
#else
	return fgets(ptr,size,fp);
#endif
}


/*****************************************************************************
*
*   FileBuildPath()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL  Nlm_FileBuildPath (Nlm_CharPtr root, Nlm_CharPtr sub_path, Nlm_CharPtr filename)

{
    Nlm_CharPtr tmp;
    Nlm_Boolean dir_start = FALSE;
#ifdef OS_VMS
  Nlm_Boolean had_root = FALSE;
#endif

    if (root == NULL)              /* no place to put it */
        return NULL;

    tmp = root;
    if (*tmp != '\0')                /* if not empty */
    {
#ifndef OS_VMS
        dir_start = TRUE;
#else
        had_root = TRUE;
#endif
        while (*tmp != '\0')
        {
#ifdef OS_VMS
            if (*tmp == '[')
                dir_start = TRUE;
#endif
            tmp++;
        }

        if ((*(tmp - 1) != DIRDELIMCHR) && (dir_start))
        {
            *tmp = DIRDELIMCHR;
            tmp++; *tmp = '\0';
        }
    }

    if (sub_path != NULL)
    {
#ifdef OS_VMS
        if (dir_start)
        {
            *(tmp-1) = '.';
            if (*sub_path == '[')
                sub_path++;
        }
        else if ((had_root) && (*sub_path != '['))
        {
            *tmp = '[';
            tmp++; *tmp = '\0';
        }
#else
        if ((dir_start) && (*sub_path == DIRDELIMCHR))
            sub_path++;
#endif
        tmp = StringMove(tmp, sub_path);
        if (*(tmp-1) != DIRDELIMCHR)
        {
            *tmp = DIRDELIMCHR;
            tmp++; *tmp = '\0';
        }
    }

    if (filename != NULL)
        StringMove(tmp, filename);

    return root;
}

/*****************************************************************************
*
*   FileNameFind()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FileNameFind (Nlm_CharPtr pathname)

{
  Nlm_CharPtr  filename;
  Nlm_Int2     len;

  if (pathname != NULL) {
    len = Nlm_StringLen (pathname);
    filename = &(pathname [len]);
    while (len > 0 && pathname [len - 1] != DIRDELIMCHR) {
      len--;
      filename--;
    }
    return filename;
  } else {
    return NULL;
  }
}


/*****************************************************************************
*
*   FilePathFind()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FilePathFind(const Nlm_Char* fullname)
{
  Nlm_CharPtr str;
  size_t	     len = Nlm_StringLen(fullname);
  if ( !len )
    return 0;

  while (len  &&  fullname[len] != DIRDELIMCHR)
    len--;

  str = (Nlm_Char*)Nlm_MemGet(len + 1, MGET_ERRPOST);
  Nlm_MemCpy(str, fullname, len);
  str[len] = '\0';
  return str;
}


/*****************************************************************************
*
*   FileLength()
*
*****************************************************************************/
NLM_EXTERN Nlm_Int8 LIBCALL Nlm_FileLength(Nlm_CharPtr fileName)
{
  Nlm_Int8 file_len = Nlm_FileLengthEx(fileName);
#ifdef MPIBLAST_PIO
  Nlm_Int8 vfile_len = 0;
  if(use_virtual_frags) {
	  vfile_len = VFMGetFileLength(fileName);
	  if(vfile_len > 0) {
		return vfile_len;
	  }
  }
#endif
  return (file_len > 0) ? file_len : 0;
}


/*****************************************************************************
*
*   FileLengthEx()
*
*****************************************************************************/
NLM_EXTERN Nlm_Int8 LIBCALL Nlm_FileLengthEx(const Nlm_Char* fileName)
{
  if (!fileName  ||  !*fileName)
    return -1;

#ifdef OS_MAC
  {{
    OSErr           err;
    HParamBlockRec  params;
    Nlm_Char        path[256];

	Nlm_MemSet ((Nlm_VoidPtr) &params, 0, sizeof (HParamBlockRec));
    Nlm_StringNCpy_0(path, fileName, sizeof(path));
    Nlm_CtoPstr((Nlm_CharPtr) path);
    params.fileParam.ioNamePtr = (StringPtr)path;
    params.fileParam.ioVRefNum = 0;
    params.fileParam.ioFDirIndex = 0;
    err = PBHGetFInfo(&params, FALSE);
    return (err == noErr) ?
            params.fileParam.ioFlLgLen : -1;
  }}
#else
  {{
    struct stat sbuf;
    return (stat(fileName, &sbuf) == 0) ? sbuf.st_size : -1;
  }}
#endif
}


/*****************************************************************************
*
*   FileDelete()
*
*****************************************************************************/
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRemove (Nlm_CharPtr fileName)

{
  Nlm_Char  local [256];

  if (fileName != NULL && fileName [0] != '\0') {
    Nlm_StringNCpy_0(local, fileName, sizeof(local));
    return (Nlm_Boolean) (remove (local) == 0);
  } else {
    return FALSE;
  }
}

/*****************************************************************************
*
*   FileRename()
*
*****************************************************************************/
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRename (Nlm_CharPtr oldFileName, Nlm_CharPtr newFileName)

{
  Nlm_Char  localnew [256];
  Nlm_Char  localold [256];

  if (oldFileName != NULL && oldFileName [0] != '\0'
    && newFileName != NULL && newFileName [0] != '\0') {
    Nlm_StringNCpy_0(localold, oldFileName, sizeof(localold));
    Nlm_StringNCpy_0(localnew, newFileName, sizeof(localnew));
    return (Nlm_Boolean) (rename (localold, localnew) == 0);
  } else {
    return FALSE;
  }
}

/*****************************************************************************
*
*   FileCreate()
*
*****************************************************************************/
#ifdef OS_MAC
static OSType Nlm_GetOSType (Nlm_CharPtr str, OSType dfault)

{
  OSType  rsult;

  rsult = dfault;
  if (str != NULL && str [0] != '\0') {
    rsult = *(OSType*) str;
  }
  return rsult;
}
#endif

NLM_EXTERN void LIBCALL Nlm_FileCreate (Nlm_CharPtr fileName, Nlm_CharPtr type, Nlm_CharPtr creator)

{
#ifdef OS_MAC
  Nlm_Int2  fError;
  Nlm_Char  temp [256];
  OSType    fType;
  OSType    fCreator;
  FSSpec    spec;
#else
  FILE      *fp;
#endif

  if (fileName != NULL && fileName [0] != '\0') {

#ifdef OS_MAC
    /* note: the following assumes either full pathname or that the current
       directory is the proper location to find/create the file */

    Nlm_StringNCpy_0(temp, fileName, sizeof(temp));
    Nlm_CtoPstr ( temp);
    fError = FSMakeFSSpec( 0, 0, (StringPtr)temp, &spec);
    
    /* file not found, so create it... */
    if( fError == fnfErr){
        fType = Nlm_GetOSType (type, 'TEXT');
        fCreator = Nlm_GetOSType (creator, '    ');
        FSpCreate( &spec, fCreator, fType, smSystemScript);
    }
#else
    fp = Nlm_FileOpen (fileName, "w");
    if (fp != NULL) {
      Nlm_FileClose (fp);
    }
#endif
  }
}

/*****************************************************************************
*
*   CreateDir(pathname)
*
*****************************************************************************/

#ifdef OS_MAC
NLM_EXTERN Nlm_Boolean LIBCALL  Nlm_CreateDir (Nlm_CharPtr pathname)
{
  if (pathname != NULL && pathname [0] != '\0') {
    OSErr err;
    err = MacCreateDirectory(pathname);
    return (Nlm_Boolean) (err == noErr || err == dupFNErr);
  }
  return FALSE;
}
#endif

#ifndef OS_MAC
NLM_EXTERN Nlm_Boolean LIBCALL  Nlm_CreateDir (Nlm_CharPtr pathname)
{
#ifndef OS_VMS
  size_t          len;
  Nlm_Char        path[PATH_MAX];
#endif
#ifdef OS_UNIX
  mode_t          oldmask;
#endif
  Nlm_Boolean     rsult = FALSE;

  if (pathname != NULL && pathname [0] != '\0') {
#if defined(OS_MSWIN) || defined(OS_NT)
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    len = Nlm_StringLen (path);
    if (len > 0 && path [len - 1] == DIRDELIMCHR) {
        path [len - 1] = '\0';
    }
    rsult = (Nlm_Boolean) (mkdir ((char *) path) == 0);
    if (errno == EACCES) { /* it's O.K. if it was already there */
	rsult = TRUE;
    }
#endif
#ifdef OS_UNIX
    oldmask = umask (0000);
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    len = Nlm_StringLen (path);
    if (len > 0 && path [len - 1] == DIRDELIMCHR) {
        path [len - 1] = '\0';
    }
    rsult = (Nlm_Boolean) (mkdir ((char *) path, 0755) == 0);
    if (errno == EEXIST) { /* it's O.K. if it was already there */
	rsult = TRUE;
    }
    umask (oldmask);
#endif
#ifdef OS_VMS
    rsult = (Nlm_Boolean) (mkdir ((char *) pathname, 0755) == 0);
#endif
  }
  return rsult;
}
#endif

/*****************************************************************************
*
*   DirectoryContents()
*
*****************************************************************************/

#ifdef OS_UNIX_DARWIN
#include <dirent.h>
#endif

NLM_EXTERN ValNodePtr LIBCALL Nlm_DirCatalog (Nlm_CharPtr pathname)

{
#ifdef OS_MAC
  long            dirID;
  OSErr           err;
  short           index;
  unsigned short  num;
  Nlm_Char        path[PATH_MAX];
  CInfoPBRec      pbc;
  HParamBlockRec  pbh;
  short           vRefNum;
#endif
#ifdef OS_UNIX
  Nlm_Uint1       choice;
#ifdef OS_UNIX_DARWIN
  DIR             *dirp;
  struct dirent   *dep;
#else
  Nlm_Char        buf [256];
  Nlm_Char        ch;
  Nlm_Char        cmmd [PATH_MAX + 20];
  FILE            *fp;
  Nlm_CharPtr     ptr;
#endif
#endif
  ValNodePtr      vnp = NULL;

  if (pathname != NULL && pathname [0] != '\0') {
#ifdef OS_MAC
    Nlm_StringNCpy_0 (path, pathname, sizeof (path));
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    Nlm_MemSet ((Nlm_VoidPtr) (&pbh), 0, sizeof (HParamBlockRec));
    pbh.volumeParam.ioNamePtr = (StringPtr) path;
    pbh.volumeParam.ioVolIndex = -1;
    err = PBHGetVInfo (&pbh, FALSE);
    if (err != noErr) return NULL;
    vRefNum = pbh.volumeParam.ioVRefNum;
    Nlm_StringNCpy_0 (path, pathname, sizeof (path));
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    Nlm_MemSet ((Nlm_VoidPtr) (&pbc), 0, sizeof (CInfoPBRec));
    pbc.dirInfo.ioNamePtr = (StringPtr) path;
    pbc.dirInfo.ioVRefNum = vRefNum;
    err = PBGetCatInfo (&pbc, FALSE);
    if (err != noErr) return NULL;
    if (pbc.dirInfo.ioFlAttrib & 16) {
      num = pbc.dirInfo.ioDrNmFls;
      dirID = pbc.dirInfo.ioDrDirID;
      for (index = 1; index <= num; index++) {
        Nlm_MemSet ((Nlm_VoidPtr) (&pbc), 0, sizeof (CInfoPBRec));
        pbc.dirInfo.ioNamePtr = (StringPtr) path;
        pbc.dirInfo.ioVRefNum = vRefNum;
        pbc.dirInfo.ioFDirIndex = index;
        pbc.dirInfo.ioDrDirID = dirID;
        pbc.dirInfo.ioACUser = 0;
        err = PBGetCatInfo (&pbc, FALSE);
        if (err == noErr) {
          Nlm_PtoCstr ((Nlm_CharPtr) path);
          if (pbc.dirInfo.ioFlAttrib & 16) {
            ValNodeCopyStr (&vnp, 1, path);
          } else {
            ValNodeCopyStr (&vnp, 0, path);
          }
        }
      }
    }
#endif
#if defined(WIN32)
    {{
      Nlm_Char x_path[PATH_MAX];
      WIN32_FIND_DATA fData;
      HANDLE hFindFile;
      Nlm_StringNCpy_0(x_path, pathname, sizeof(x_path) - 5);
      Nlm_StringCat(x_path, "\\*.*");
      hFindFile = FindFirstFile(x_path, &fData);
      if (hFindFile == INVALID_HANDLE_VALUE)
        return 0;
      do {
        if (fData.cFileName[0] != '.'  ||
            (fData.cFileName[1] != '.'  &&  fData.cFileName[1] != '\0'))
          ValNodeCopyStr
            (&vnp, (fData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ? 1 : 0,
             fData.cFileName);
      } while ( FindNextFile(hFindFile, &fData) );
      FindClose(hFindFile);
    }}
#endif
#ifdef OS_UNIX
#ifdef OS_UNIX_DARWIN
    dirp = opendir(pathname);
    if (dirp == NULL) return NULL;
    while ((dep = readdir(dirp)) != NULL) {
      /* ignore 'invisible' files. */
      if (dep->d_namlen < 1 || dep->d_name[0] == '.')
        continue;
      if (dep->d_type == DT_DIR) /* directory */
        choice = 1;
      else          /* all other file types. */
        choice = 0;
      ValNodeCopyStr (&vnp, choice, dep->d_name);
    }
    closedir(dirp);
#else
    sprintf (cmmd, "ls -1p %s 2>/dev/null", pathname);
    fp = popen (cmmd, "r");
    if (fp == NULL) return NULL;
    while (Nlm_FileGets (buf, sizeof (buf), fp) != NULL) {
      ptr = buf;
      ch = *ptr;
      while (ch != '\0' && ch != '\n' && ch != '\r') {
        ptr++;
        ch = *ptr;
      }
      *ptr = '\0';
      choice = 0;
      ptr = Nlm_StringChr (buf, '/');
      if (ptr != NULL) {
        *ptr = '\0';
        choice = 1;
      }
      ValNodeCopyStr (&vnp, choice, buf);
    }
    pclose (fp);
#endif
#endif
#ifdef OS_VMS
#endif
  }
  return vnp;
}

/*****************************************************************************
*
*   general file recursion functio
*
*****************************************************************************/

NLM_EXTERN Nlm_Int4 Nlm_DirExplore (
  Nlm_CharPtr directory,
  Nlm_CharPtr filter,
  Nlm_CharPtr suffix,
  Nlm_Boolean recurse,
  Nlm_DirExpProc proc,
  Nlm_VoidPtr userdata
)

{
  Nlm_Int4     count = 0;
  Nlm_Char     file [FILENAME_MAX], path [PATH_MAX];
  Nlm_CharPtr  ptr, str;
  ValNodePtr   head, vnp;

  if (proc == NULL) return 0;
  if (Nlm_StringHasNoText (directory) /* || Nlm_StringHasNoText (suffix) */ ) return 0;

  /* get list of all files in source directory */

  head = Nlm_DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      str = (Nlm_CharPtr) vnp->data.ptrvalue;
      if (! Nlm_StringHasNoText (str)) {

        /* check filename for indicated suffix */

        ptr = Nlm_StringStr (str, suffix);
        if (ptr != NULL) {

          /* make sure detected suffix is really at end of filename */

          if (Nlm_StringCmp (ptr, suffix) == 0) {
            *ptr = '\0';
          } else {
            ptr = NULL;
          }
        }

        if (Nlm_StringHasNoText (suffix) || ptr != NULL) {

          Nlm_StringNCpy_0 (path, directory, sizeof (path));
          sprintf (file, "%s%s", str, suffix);
          Nlm_FileBuildPath (path, NULL, file);

          /* check full path/file name for desired filter */

          if (Nlm_StringHasNoText (filter) || Nlm_StringStr (path, filter) != NULL) {

            /* process file that satisfies optional filter and suffix constraints */

            proc (path, userdata);
            count++;
          }
        }
      }
    } else if (vnp->choice == 1 && recurse) {

      /* recurse into subdirectory */

      Nlm_StringNCpy_0 (path, directory, sizeof (path));
      str = (Nlm_CharPtr) vnp->data.ptrvalue;
      Nlm_FileBuildPath (path, str, NULL);

      count += Nlm_DirExplore (path, filter, suffix, recurse, proc, userdata);
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);

  return count;
}

/*****************************************************************************
*
*   TmpNam()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_TmpNam (Nlm_CharPtr s)

{
#ifdef TEMPNAM_AVAIL
    char *filename;
    static Nlm_Char save_filename[PATH_MAX];

    /* emulate tmpnam(), except get the benefits of tempnam()'s ability to */
    /* place the files in another directory specified by the environment   */
    /* variable TMPDIR                                                     */
#ifdef OS_UNIX_DARWIN
    filename = tempnam("/tmp", "ncbi.");
#else
    filename = tempnam(NULL, NULL);
#endif
    if (s == NULL)
    { /* return pointer to static string */
        if (filename != NULL) {
          strcpy ((char *) save_filename, (char *) filename);
          free ((void *) filename);
        } else {
          save_filename [0] = '\0';
        }
        return save_filename;
    } else {
        if (filename != NULL) {
          strcpy ((char *) save_filename, (char *) filename);
          Nlm_StrCpy (s, save_filename);
          free ((void *) filename);
        } else {
          *s = '\0';
        }
        return s;
    }
#else
#ifdef OS_MAC
    static Nlm_Char  directory [PATH_MAX];
    OSErr        err;
    long         gesResponse;
    long         newDirID;
    short        newVRefNum;
    CInfoPBRec   params;
    Nlm_Char     temp [PATH_MAX];
    Nlm_CharPtr  tmp;
    Nlm_Boolean  useTempFolder;
    char * filename;

    useTempFolder = FALSE;
    if (! Gestalt (gestaltFindFolderAttr, &gesResponse) &&
        (gesResponse & (1 << gestaltFindFolderPresent))) {
      err = FindFolder(kOnSystemDisk, kTemporaryFolderType,
                       kCreateFolder, &newVRefNum, &newDirID);
      if (err == noErr) {
        useTempFolder = TRUE;
      }
    }
    filename = tmpnam (NULL);
    if (useTempFolder) {
      temp [0] = '\0';
      params.dirInfo.ioNamePtr = (StringPtr) directory;
      params.dirInfo.ioDrParID = newDirID;
      do {
        params.dirInfo.ioVRefNum = newVRefNum;
        params.dirInfo.ioFDirIndex = -1;
        params.dirInfo.ioDrDirID = params.dirInfo.ioDrParID;
        err = PBGetCatInfo (&params, FALSE);
        Nlm_PtoCstr ((Nlm_CharPtr) directory);
        Nlm_StringCat (directory, DIRDELIMSTR);
        Nlm_StringCat (directory, temp);
        Nlm_StringCpy (temp, directory);
      } while (params.dirInfo.ioDrDirID != fsRtDirID);
      tmp = Nlm_StringMove (directory, temp);
      tmp = Nlm_StringMove (tmp, (Nlm_CharPtr) filename);
      if (s == NULL) {
          return (Nlm_CharPtr) directory;
      } else {
          s [0] = '\0';
          Nlm_StringCpy (s, directory);
          return s;
      }
    } else {
      if (s == NULL) {
          return (Nlm_CharPtr) filename;
      } else {
          s [0] = '\0';
          Nlm_StringCpy (s, filename);
          return s;
      }
    }
#else
    char * filename;

    filename = tmpnam (NULL);
    if (s == NULL) {
        return (Nlm_CharPtr) filename;
    } else {
        s [0] = '\0';
        Nlm_StringCpy (s, filename);
        return s;
    }
#endif
#endif
}

/* FileCache provides buffered text read for handling Unix, Mac, and DOS line endings gracefully */

/* attach file pointer (text read mode expected) to cache object (usually on stack) */

#ifdef OS_MSWIN
#include <fcntl.h>
#include <io.h>
#endif

NLM_EXTERN Nlm_Boolean Nlm_FileCacheSetup (
  Nlm_FileCache PNTR fcp,
  FILE *fp
)

{
  if (fp == NULL || fcp == NULL) return FALSE;

#ifdef OS_MSWIN
  _setmode (_fileno (fp), _O_BINARY);
#endif

  MemSet ((Nlm_VoidPtr) fcp, 0, sizeof (Nlm_FileCache));

  fcp->fp = fp;
  fcp->offset = ftell (fp);

  return TRUE;
}

static void Nlm_FileCacheReadBlock (
  Nlm_FileCache PNTR fcp
)

{
  int  total;

  if (fcp == NULL || fcp->fp == NULL) return;

  if (fcp->ctr >= fcp->total) {
    fcp->offset += (Nlm_Int4) fcp->total;
    fcp->ctr = 0;
    fcp->total = 0;

    fcp->buf [0] = '\0';
    total = (int) Nlm_FileRead ((Nlm_VoidPtr) fcp->buf, sizeof (Nlm_Char), (size_t) 512, fcp->fp);
    if (total < 0 || total > 512) {
      total = 512;
    }

    fcp->buf [total] = '\0';
    fcp->total = Nlm_StringLen (fcp->buf);
  }
}

/* equivalent of getc */

static Nlm_Char Nlm_FileCacheGetChar (
  Nlm_FileCache PNTR fcp
)

{
  Nlm_Char  ch = '\0', nxt;

  if (fcp == NULL || fcp->fp == NULL) return ch;

  /* read a fresh block if buffer is empty */

  if (fcp->ctr >= fcp->total) {
    Nlm_FileCacheReadBlock (fcp);
  }

  /* get next character in buffer */

  if (fcp->ctr < fcp->total) {
    ch = fcp->buf [(int) fcp->ctr];
    (fcp->ctr)++;
  }

  if (ch == '\n' || ch == '\r') {
    if (fcp->ctr >= fcp->total) {
      Nlm_FileCacheReadBlock (fcp);
    }
    if (fcp->ctr < fcp->total) {

      /* look for carriage return / linefeed pair - DOS file read on Mac or Unix platform */

      nxt = fcp->buf [(int) fcp->ctr];

      /* advance past second character in cr/lf pair */

      if (ch == '\n' && nxt == '\r') {
        (fcp->ctr)++;
      } else if (ch == '\r' && nxt == '\n') {
        (fcp->ctr)++;
      }
    }

    /* cr or lf returned as newline */

    ch = '\n';
  }

  return ch;
}

/* equivalent of fgets, leaves /n at end of string */

NLM_EXTERN Nlm_CharPtr Nlm_FileCacheGetString (
  Nlm_FileCache PNTR fcp,
  Nlm_CharPtr str,
  size_t size
)

{
  Nlm_Char     ch;
  Nlm_Uint2     count;
  Nlm_CharPtr  ptr;

  if (fcp == NULL || fcp->fp == NULL || str == NULL || size < 1) return NULL;

  ch = Nlm_FileCacheGetChar (fcp);
  count = 0;
  ptr = str;

  while (ch != '\0' && ch != '\n' && ch != '\r' && count < size - 2) {
    *ptr = ch;
    ptr++;
    count++;
    ch = Nlm_FileCacheGetChar (fcp);
  }

  if (ch == '\n' || ch == '\r') {
    *ptr = '\n';
    ptr++;
    count++;
  } else if (ch != '\0') {
    *ptr = ch;
    ptr++;
    count++;
  }
  *ptr = '\0';

  if (count < 1) return NULL;

  return str;
}

/* smarter fgets removes newline from end of string */

NLM_EXTERN Nlm_CharPtr Nlm_FileCacheReadLine (
  Nlm_FileCache PNTR fcp,
  Nlm_CharPtr str,
  size_t size,
  Nlm_BoolPtr nonewline
)

{
  Nlm_Char     ch;
  Nlm_CharPtr  ptr;
  Nlm_CharPtr  tmp;

  if (fcp == NULL || fcp->fp == NULL || str == NULL || size < 1) return NULL;
  *str = '\0';
  tmp = Nlm_FileCacheGetString (fcp, str, size);
  if (tmp != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch != '\n' && ch != '\r') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
    if (nonewline != NULL) {
      if (ch != '\n' && ch != '\r') {
        *nonewline = TRUE;
      } else {
        *nonewline = FALSE;
      }
    }
  }
  return tmp;
}

NLM_EXTERN void Nlm_FileCacheSeek (
  Nlm_FileCache PNTR fcp,
  Nlm_Int4 pos
)

{
  if (fcp == NULL || fcp->fp == NULL) return;

  /*
  if (fcp->offset <= pos && fcp->offset + (Nlm_Int4) fcp->total >= pos) {
    fcp->ctr = (Nlm_Int2) (pos - fcp->offset);
    return;
  }
  */

  fcp->ctr = 0;
  fcp->total = 0;
  fcp->offset = pos;

  fseek (fcp->fp, pos, SEEK_SET);
}

NLM_EXTERN Nlm_Int4 Nlm_FileCacheTell (
  Nlm_FileCache PNTR fcp
)

{
  Nlm_Int4  bytes;
  Nlm_Int4  offset;

  if (fcp == NULL || fcp->fp == NULL) return 0L;

  offset = ftell (fcp->fp);
  bytes = (Nlm_Int4) (fcp->total - fcp->ctr);
  offset -= bytes;

  return offset;
}

NLM_EXTERN Nlm_Boolean Nlm_FileCacheFree (
  Nlm_FileCache PNTR fcp,
  Nlm_Boolean restoreFilePos
)

{
  Nlm_Int4  pos;

  if (fcp == NULL || fcp->fp == NULL) return FALSE;

  if (restoreFilePos) {

    /* correct position of file pointer */

    pos = Nlm_FileCacheTell (fcp);
    fseek (fcp->fp, pos, SEEK_SET); 
  }

  MemSet ((Nlm_VoidPtr) fcp, 0, sizeof (Nlm_FileCache));

  return TRUE;
}

