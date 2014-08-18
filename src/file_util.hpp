/*************************************************************************
** Copyright 2009 by Virginia Polytechnic Institute and State
** University. All rights reserved. Virginia Polytechnic Institute and
** State University (Virginia Tech) owns the mpiBLAST software and its
** associated documentation ("Software"). You should carefully read the
** following terms and conditions before using this software. Your use
** of this Software indicates your acceptance of this license agreement
** and all terms and conditions.
** 
** The parallel input/output (PIO) portion of this work was based in
** part on published research that was supported by the North Carolina
** State University and Oak Ridge National Laboratory. The actual
** implementation was completed at Virginia Tech in the summer of 2007.
** 
** Additionally, portions of this work are derivative works of the NCBI
** C Toolkit. Although, the NCBI C Toolkit is released freely without
** restriction under the terms of their license, the following files
** listed below, have been modified by Virginia Tech, and as such are
** redistributed under the terms of this license.
** 
** ncbi/api/txalign.c
** ncbi/corelib/ncbifile.c
** ncbi/object/objalign.c
** ncbi/tools/blast.c
** ncbi/tools/blastdef.h
** ncbi/tools/blastool.c
** ncbi/tools/blastutl.c
** ncbi/tools/blfmtutl.c
** ncbi/tools/ncbisam.c
** ncbi/tools/readdb.c
** ncbi/tools/readdb.h
** 
** License:
** 
** This file is part of mpiBLAST.
** 
** mpiBLAST is free software: you can redistribute it and/or modify it 
** under the terms of the GNU General Public License version 2 as published 
** by the Free Software Foundation. 
** 
** Accordingly, mpiBLAST is distributed in the hope that it will be
** useful, but WITHOUT ANY WARRANTY; without even the implied warranty
** of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
** General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with mpiBLAST. If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/
#ifndef __MPIBLAST_FILE_UTIL
#define __MPIBLAST_FILE_UTIL

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>
#include <vector>
#include "mpiblast_types.h"

// copy styles 
#define COPY_VIA_NONE 0 /* don't copy db fragments. useful for PVFS,GPFS,PFS */
#define COPY_VIA_CP   1 /* copy db fragments using cp */
#define COPY_VIA_RCP  2 /* copy db fragments using rcp */
#define COPY_VIA_SCP  3 /* copy db fragments using scp */
#define COPY_VIA_MPI  4 /* copy db fragments using MPI */

#define COPY_VIA_XCP  1 /* copyFile calls that use this need to be changed. this is a hack for right now */



#ifdef WIN32
const std::string PATH_SEPARATOR = "\\";
#else
#include "unistd.h"
#include <sys/types.h>
#include <sys/param.h>
const std::string PATH_SEPARATOR = "/";
#if defined(__MACH__)
# include <sys/mount.h>
# define STATFS statfs
#else
# include <sys/statvfs.h>
# define STATFS statvfs
#endif

#endif

/** Extension for a database's local fragment list file */
#define	FRAG_LIST_EXTENSION 	".mbf"

bool doesFileExist(const std::string& file_name );
std::string getPath(const std::string& file_path );
std::string getFilename(const std::string& file_path);
void getTempFileName(std::string& tempname);
int moveFile(const std::string& src, const std::string& dest);
int copyFile(const std::string& src, const std::string& dest, const int fileType);
int removeFile(const std::string& filename, bool verbose = false);

bool removeDirectory( const std::string& dir );

/** get the size in bytes of a particular file */
uint64 statFileSize( const char* path );

/** get the modification time of a file */
time_t statFileMTime( const char* path );

/** 
 * Check that a directory exists and if not try to create it
 * @return true if the directory exists or was created, false otherwise
 */
bool checkCreateDir( const std::string& path );
bool checkCreateMultiDir( const std::string& path );

/** 
 * Check that we have write permissions to a directory
 * @return true if we have write permission, false otherwise
 */
bool checkDirWritePerms( const std::string& path );

/**
 * Set an environment variable in a cross platform manner
 */
void setEnvironmentVariable( const char* var );

/**
 * Check whether a path specifies an rcp or scp host
 * If a path contains a ":" character and the ":" is not part of a
 * Windows drive identifier then return true.  fixme: how can
 * we distinguish Mac OS X paths with : from remote hosts?
 */
bool isRemotePath( const std::string& path );

/**
 * Register a file name to be deleted before the process exits
 * When passed an empty string, it does not add to the list of files to delete
 * @param fname	The name of a file to delete, empty strings are ignored
 * @return A vector of file names registered for deletion
 */
std::vector< std::string >& registerFileToDelete( std::string fname = "" );

/* Returns free space in KB for a filesystem in a given path */
long freeSpaceKB(const std::string& path);

/**
 * Deletes files associated with mpiBLAST from a given directory
 * Recursively deletes temp directories associated with mpiBLAST
 * this function needs to execute the equivalent of
 * rm -fr /$TMPDIR/mpiblast_tmpXXXXXX
 * rm -fr /$TMPDIR/*.[np]?? /$TMPDIR/*.mbf /$TMPDIR/*.dbs
 * rm -fr /$TMPDIR/*.[np]??XXXXXX /$TMPDIR/*.mbfXXXXXX /$TMPDIR/*.dbsXXXXXX
 */
void deleteMpiBlastFiles( const std::string& directory );
int Posix_Set_lock(int fd, int cmd, int type, off_t offset, int whence, off_t len);

#endif // __MPIBLAST_FILE_UTIL
