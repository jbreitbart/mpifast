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

#include "mpi.h"
#include "file_util.hpp"
#include "mpiblast.hpp"
#include "mpiblast_types.h"
using namespace std;

#ifdef WIN32
#include <windows.h>
#else
#include <dirent.h>
#include "sys/types.h"
#include "sys/stat.h"
#endif

bool doesFileExist( const string& file_name )
{
	// this is a platform agnostic method of checking for existence,
	// but opening the file is somewhat undesirable.
	ifstream check_file( file_name.c_str(), ios::in );
	if( check_file.is_open() ){
		check_file.close();
		return true;
	}
	return false;
}

/** Utility function to get the path component of a file path */
string getPath( const string& file_path )
{
  string::size_type slashloc;
  string::size_type search_loc = string::npos;
  for(;;){
    slashloc = file_path.rfind( PATH_SEPARATOR, search_loc );
    if( slashloc > 0 && slashloc != string::npos )
      if( file_path[ slashloc - 1 ] == '\\' ){
	search_loc = slashloc - 1;
	continue;
      }
    break;
  }
  string path_str;
  if( slashloc != string::npos )
    path_str = file_path.substr( 0, slashloc + 1 );
  else
    path_str = "";

  return path_str;
}

/** Utility function to get the file name component of a file path */
string getFilename( const string& file_path )
{
  string::size_type slashloc;
  string::size_type search_loc = string::npos;
  for(;;){
    slashloc = file_path.rfind( PATH_SEPARATOR, search_loc );
    if( slashloc > 0 && slashloc != string::npos )
      if( file_path[ slashloc - 1 ] == '\\' ){
	search_loc = slashloc - 1;
	continue;
      }
    break;
  }
  string file_str;
  if( slashloc != string::npos )
    file_str = file_path.substr( slashloc + 1 );
  else
    file_str = file_path;

  return file_str;
}

/** Utility function to get a name for a temp file */
void getTempFileName( string& tempname )
{

	//char* temp_char_name = new char[ tempname.size() + 1 ];
	char* temp_char_name = (char*)malloc((tempname.size() + 1)*sizeof(char));
	
	if( temp_char_name == 0 ){
		throw __FILE__ "(getTempFileName): Unable to allocate memory for temp file" ;
	}
	
	strcpy( temp_char_name, tempname.c_str() );
	
	if( debug_msg ){
		LOG_MSG << "Temp name base: " << tempname << endl;
	}
	
#ifdef WIN32
	char* win_name = _mktemp( temp_char_name );
	tempname = win_name;
	// mkstemp actually creates the file, windows mktemp() doesn't.
	ofstream tmp_file( win_name );
#else
	int fd = mkstemp( temp_char_name );
	
	if(fd < 0){
		cerr << "Unable to make temp file name " << temp_char_name 
		     << " from " << tempname << endl;
		cerr << "fd = " << fd << endl;		
		cerr << "Error Msg = \"" << strerror(errno) << "\"" << endl;
		
		throw __FILE__ "(getTempFileName): Unable to make temp file" ;
	}
	
	// Close the file descriptor (we only needed to create an empty file)
	close(fd);
	
	tempname = temp_char_name;
#endif
		
	if( debug_msg ){
		LOG_MSG << "Got temp name: " << tempname << endl;
	}
	
	//delete[] temp_char_name;
	free(temp_char_name);
}

/** Utility function to move a file */
int moveFile( const string& src, const string& dest )
{
#ifdef WIN32
	return !MoveFileEx( src.c_str(), dest.c_str(), MOVEFILE_COPY_ALLOWED | MOVEFILE_REPLACE_EXISTING );
#else
	
	// Updated version that will use rcp if requested.
	// Note that we are expecting only the destination string to 
	// possibly have a hostname prepended to it.
	string::size_type dest_pos = dest.find(":", 0);
	
	// Do we need to use rcp?
	if( dest_pos == string::npos ){
		// Don't use rcp, use mv
		// This is the original behavior
		
		// First, check to make sure that src != dest
		if(src == dest){
			return EXIT_SUCCESS;
		}
		
		string mv_cmd;

		mv_cmd = "mv " + src + " " + dest;
		
		//return system( mv_cmd.c_str() );
		return rename( src.c_str(), dest.c_str() );
	}
	else{
		// Don't use cp, use rcp
		
		// First, check to make sure that src != dest
		// Note that we trust the HOSTNAME environment variable
		// Build the fully qualified src path, i.e.:
		//	hostname:/local/path/to/src
		string full_src_name = getenv( "HOSTNAME" );
		
		// Do we need to append the full path?
		if(src.find(PATH_SEPARATOR, 0) == string::npos){
			// Yes, append the path to src
			char dir_buf[1024];
			
			string curr_dir = getcwd(dir_buf, 1024);
			
			full_src_name = full_src_name + ":" + curr_dir + PATH_SEPARATOR + src;
			
		}
		else{
			// src already contains a path. Lets hope its an absolute path!
			full_src_name = full_src_name + ":" + src;
		}
		
		if(full_src_name == dest){
			return EXIT_SUCCESS;
		}
		
		string rcp_cmd;
		int ret_value;
		
		rcp_cmd = "rcp " + src + " " + dest;
		
		ret_value = system( rcp_cmd.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			// DEBUG
			cerr << "rcp command failed!" << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			
			return ret_value;
		}
		
		// Unlink the src file to complete the move
		ret_value = unlink( src.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			// DEBUG
			cerr << "unlink command failed!" << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			
			return ret_value;
		}
		
		return ret_value;
	}
	
#endif
}

/** Utility function to copy a file */
int copyFile( const string& src, const string& dest, const int copyType )
{
#ifdef WIN32
	return !CopyFile( src.c_str(), dest.c_str(), false );
#else
	int ret_value;
	string xcp_cmd;

	// What type of "xcp" to use
	if( copyType == COPY_VIA_CP ){
		xcp_cmd = "cp " + src + " " + dest;
	} else if (copyType == COPY_VIA_RCP) {
		xcp_cmd = "rcp " + src + " " + dest;
	} else if (copyType == COPY_VIA_SCP) {
		xcp_cmd = "scp " + src + " " + dest;
	} else if (copyType == COPY_VIA_NONE) {
		return EXIT_SUCCESS;
	} else {
		//UNKNOWN CP command what should we do?
		xcp_cmd = "";
	}
	
	ret_value = system( xcp_cmd.c_str() );
		
	if(ret_value != EXIT_SUCCESS){
		cerr << "cp command failed!" << endl;
		cerr << "command: " << xcp_cmd << endl;
		cerr << "source = " << src << endl;
		cerr << "dest = " << dest << endl;
		cerr << "ret_value = " << ret_value << endl;
	}
		
	return ret_value;
#endif
}

/** Utility function to delete a file */
int removeFile( const string& filename, bool verbose )
{
#ifdef WIN32
	int success = DeleteFile( filename.c_str() );
	if( !success && verbose ){
		char msg[5001];
		memset( msg, 0, 5001 );
		int errcode = GetLastError();
		FormatMessage( FORMAT_MESSAGE_FROM_SYSTEM, NULL, errcode, 0, msg, 5000, NULL );
		cerr << "DeleteFile: " << msg << endl;
	}
	return success == 0;
#else
	errno = 0;
	if ( unlink(filename.c_str()) == -1){
		string err_str = "Could not remove " + filename +":";
		perror(err_str.c_str());
		return errno;
	}
	if (verbose)
		cerr << "Removed " << filename << endl;
	return 0;
#endif
}

/** Get the size of a file */
uint64 statFileSize( const char * path ) {
#ifdef WIN32
	WIN32_FILE_ATTRIBUTE_DATA file_data;
	uint64 f_size;
	GetFileAttributesEx( path, GetFileExInfoStandard, (void*)&file_data );
	f_size = file_data.nFileSizeHigh;
	f_size <<= 32;
	f_size += file_data.nFileSizeLow;
	return f_size;
#else
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
	}
	return stat_data.st_size;
#endif
}

/** Get the last modified time of a file */
time_t statFileMTime( const char * path ) {
#ifdef WIN32
	struct __stat64 stat_data;
	if( _stat64( path, &stat_data ) ){
		perror(path);
	}
	return stat_data.st_mtime;
#else
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
	}
	return stat_data.st_mtime;
#endif
}

/* Check that a directory exists and if not try to create it, support multi-layers */
bool checkCreateMultiDir( const string& path ){

	string last_char = path.substr(path.length() - 1);
	string full_path;

	if(last_char == "/") {
		full_path = path;
	} else {
		full_path = path + "/";
	}

	string::size_type start = 1;
	while(true) {
		string::size_type loc = full_path.find("/", start);
		if(loc == string::npos) {
			break;
		}

		string curr_path = full_path.substr(0, loc);

		if(!checkCreateDir(curr_path)) {
			return 0;
		}
		
		start = loc + 1;
	}

	return 1;
}

/** 
 * Check that a directory exists and if not try to create it
 * @return true if the directory exists or was created, false otherwise
 */
bool checkCreateDir( const string& path ){
	if( isRemotePath(path) )
		return true;	// this is a remote path for scp/rcp.  just assume it exists
		
#ifdef WIN32
	SECURITY_ATTRIBUTES secattr;
	secattr.nLength = sizeof( SECURITY_ATTRIBUTES );
	secattr.lpSecurityDescriptor = NULL;
	secattr.bInheritHandle = FALSE;
	int rval = CreateDirectory( path.c_str(), &secattr );
	if( rval != 0 ){
		return GetLastError() == ERROR_ALREADY_EXISTS;
	}
	return TRUE;
#else

	int rval;
	struct stat stat_buf;
	mode_t mode = 0x1fd; // mode 775: set the perms lax and let umask control
	if ((rval = mkdir( path.c_str(), mode )) == -1 && errno == EEXIST) {
		stat( path.c_str(), &stat_buf );
		return S_ISDIR( stat_buf.st_mode );
	}
	// if the directory was successfully created then return true
	return rval == 0;
#endif
}

/** 
 * Check that we have write permissions to a directory
 * @return true if we have write permission, false otherwise
 */
bool checkDirWritePerms( const string& path ){
#ifdef WIN32
	// this is a platform agnostic method of checking for write perms,
	// but creating a file is somewhat undesirable.
	// Note for all you security buffs:  the file removal here introduces 
	// a potentially dangerous race condition.  I'm the lazy programmer
	// that you guys hate. 
	string check_fname = path + PATH_SEPARATOR + "please_delete_me.now";
	bool already_exists = doesFileExist( check_fname );
	ofstream check_file( check_fname.c_str(), ios::out | ios::ate );
	if( check_file.is_open() ){
		check_file.close();
		removeFile( check_fname );
		return true;
	}
	removeFile( check_fname );
	return false;
#else
	return access( path.c_str(), R_OK|W_OK|X_OK ) == 0;
#endif
}

/**
 * Set an environment variable in a cross platform manner
 */
void setEnvironmentVariable( const char* var ){
#ifdef WIN32
	//Silly Windows, silly underscores
	_putenv( var );
#else
	putenv( const_cast<char*>(var) );
#endif
}

bool isRemotePath( const std::string& path ){
	string::size_type colon_loc = path.find(":");
	if( colon_loc == string::npos ||
		(colon_loc == 1 && path.find("\\") == 2) ||	// watch out for windows drive specifier
		path.find("/") < colon_loc)
	{
		return false;
	}
	return true;
}

std::vector< std::string >& registerFileToDelete( std::string fname ) {
	// since this vector is needed when atexit() is called we allocate it
	// on the heap so its destructor won't get called
	static vector< string >* files = new vector< string >();
	if( fname != "" )
		files->push_back( fname );
	return *files;
}

/**
 * determine how much free disk space is at the given path
 */
long freeSpaceKB(const std::string& path){
#ifndef WIN32
	struct STATFS sfs;
	int rc = STATFS(path.c_str(), &sfs);
	size_t avail = sfs.f_bavail;
	return ((sfs.f_bsize*avail)/1024);
#else
	return LONG_MAX;
#endif
}

bool removeDirectory( const string& dir ){
	bool success = false;
#ifdef WIN32
	success = RemoveDirectory( dir.c_str() );
#else
	success = remove( dir.c_str() ) == 0;
#endif
	return success;
}

void deleteMpiBlastFiles( const std::string& directory ){
	if( directory.length() == 0 )
		return;	// something must be wrong-- at least "." should be given

	static vector< string > extensions;
	// if this is the first time through this function
	// build the list of file extensions that should be deleted
	if( extensions.size() == 0 ){
		const vector< string >& Nextensions = FragmentExtensions( string("n") );
		const vector< string >& Pextensions = FragmentExtensions( string("p") );
		extensions.insert( extensions.end(), Nextensions.begin(), Nextensions.end() );
		extensions.insert( extensions.end(), Pextensions.begin(), Pextensions.end() );
		extensions.push_back( ".nal" );
		extensions.push_back( ".pal" );
		extensions.push_back( ".mbf" );
		extensions.push_back( ".dbs" );
	}
	const string& tempdir_base = MpiBlastConfig::TempDirBaseName();


	// Search each directory for mpiBLAST-related files
	// delete files as they are found.
	// Recursively search mpiBLAST-related temp directories

// platform-specific looping over all files in a directory 
// initiate loop here
#ifdef WIN32
	WIN32_FIND_DATA FindFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;
	string m_dir = directory + PATH_SEPARATOR + "*";
	hFind = FindFirstFile(m_dir.c_str(), &FindFileData);
	if( hFind == INVALID_HANDLE_VALUE )
		return;
	do{
		string filename = FindFileData.cFileName;
		bool is_directory = (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
#else
	DIR* dirp = opendir(directory.c_str());
	dirent* dp;

	while((dp = readdir(dirp)) != NULL) {
		string filename = dp->d_name;
		string m_fullpath = directory + PATH_SEPARATOR + filename;
		struct stat stat_data;
		if( stat( m_fullpath.c_str(), &stat_data) )
			perror( m_fullpath.c_str() );
		bool is_directory = S_ISDIR( stat_data.st_mode );
#endif

		string fullpath = directory + PATH_SEPARATOR + filename;

		// check whether we're looking at a directory
		if( is_directory ){
			// recurse if it's an mpiblast directory.
			if( filename.find( tempdir_base ) == 0 ){
				// this matches, recurse
				deleteMpiBlastFiles( fullpath );
				removeDirectory( fullpath );
			}
			continue;
		}

		// check whether this file has an extension we're interested in
		for( uint extI = 0; extI < extensions.size(); extI++ ){
			string::size_type pos = string::npos;
			if( (pos = filename.find( extensions[ extI ] )) == string::npos )
				continue;	// no match
			// check that filename conforms to either .ext or .extXXXXXX
			// since some temp files have that name structure
			if( pos != filename.size() - 4 && pos != filename.size() - 10 )
				continue;	// this is not an mpiblast file but it somehow
							// managed to get created with an mpiblast extension
							// somewhere in its name.
			
			// this file matches the extension
			removeFile( fullpath, false );
		}

// platform-specific loop termination
#ifdef WIN32
	}while( FindNextFile(hFind, &FindFileData) );
	FindClose( hFind );
#else
	}
	closedir(dirp);
#endif

}

int Posix_Set_lock(int fd, int cmd, int type, off_t offset, int whence,
         off_t len)
{
    int err, error_code;
    struct flock lock;

    lock.l_type = type;
    lock.l_start = offset;
    lock.l_whence = whence;
    lock.l_len = len;

    do {
    err = fcntl(fd, cmd, &lock);
    } while (err && (errno == EINTR));

    if (err && (errno != EBADF)) {
    fprintf(stderr, "File locking failed in ADIOI_Set_lock. If the file system is NFS, you need to use NFS version 3 and mount the directory with the 'noac' option (no attribute caching).\n");     MPI_Abort(MPI_COMM_WORLD, 1);
    }

    error_code = (err == 0) ? MPI_SUCCESS : MPI_ERR_UNKNOWN;
    return error_code;
}
