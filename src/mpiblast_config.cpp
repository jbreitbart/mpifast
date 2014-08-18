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

#include "mpiblast_config.hpp"
#include "file_util.hpp"
#include <iostream>

#ifdef WIN32
// for windows mktemp
#include <io.h>
#endif

using namespace std;

extern "C"{
#include "ncbi.h"
//#include "ncbiwin.h"
}

MpiBlastConfig::MpiBlastConfig( const string& filename ) {
	
	if( filename == ".ncbirc" ){
		//
		// try reading the .ncbirc config file for mpiBLAST values
		//
		char path_buffer[ PATH_MAX ] = "";
		char* tmp_path;
	
		// read shared path
		FindPath( "ncbi", "mpiBLAST", "Shared", path_buffer, PATH_MAX );
		shared_db_path = path_buffer;
		// try environment if .ncbirc failed
		if( shared_db_path.length() == 0 && 
			(tmp_path = getenv( "MPIBLAST_SHARED" )) != NULL ){
			shared_db_path = tmp_path;
		}
		// if it still wasn't found then throw an exception!
		if( shared_db_path.length() == 0 )
			throw "Unable to read mpiBLAST shared storage path from either .ncbirc or the MPIBLAST_SHARED environment variable.\nThe mpiBLAST configuration in .ncbirc should look like:\n[mpiBLAST]\nShared=/path/to/shared/storage\nLocal=/path/to/local/storage";
			
		// read local path
		FindPath( "ncbi", "mpiBLAST", "Local", path_buffer, PATH_MAX );
		local_db_path = path_buffer;

		if( local_db_path.length() != 0 && 
			!checkCreateMultiDir( local_db_path ) && 
			(tmp_path = getenv( "MPIBLAST_LOCAL" )) != NULL ){
			cerr << "Warning:  local db path given as " << local_db_path << " in .ncbirc is unreadable.  Trying MPIBLAST_LOCAL environment variable.\n";
			local_db_path = tmp_path;
		}
		if( local_db_path.length() == 0 && 
			(tmp_path = getenv( "MPIBLAST_LOCAL" )) != NULL ){
			local_db_path = tmp_path;
		}

		if( local_db_path.length() != 0 && 
			!checkCreateMultiDir( local_db_path ) && 
			(tmp_path = getenv( "TMPDIR" )) != NULL ){
			cerr << "Warning:  local db path given as " << local_db_path << " in .ncbirc is unreadable.  Trying $TMPDIR/mpiBLAST_local_db.\n";
			local_db_path = tmp_path + PATH_SEPARATOR + "mpiBLAST_local_db";
		}

		// try setting it to $TMPDIR/mpiBLAST_local_db and print a warning
		if( local_db_path.length() == 0 && 
			(tmp_path = getenv( "TMPDIR" )) != NULL ){
			local_db_path = tmp_path + PATH_SEPARATOR + "mpiBLAST_local_db";
			cerr << "Warning: mpiBLAST local storage path unspecified in either .ncbirc or the environment, defaulting to " << local_db_path << endl;
		}
#ifdef WIN32
		// try setting it to $TEMP/mpiBLAST_local_db under windows
		if( local_db_path.length() == 0 && 
			(tmp_path = getenv( "TEMP" )) != NULL ){
			local_db_path = tmp_path + PATH_SEPARATOR + "mpiBLAST_local_db";
			cerr << "Warning: mpiBLAST local storage path unspecified in either .ncbirc or the environment, defaulting to " << local_db_path << endl;
		}
#else
		// finally default to /tmp/mpiBLAST_local_db on unix
		if( local_db_path.length() == 0  ){
			local_db_path = "/tmp/mpiBLAST_local_db";
			cerr << "Warning: mpiBLAST local storage path unspecified in either .ncbirc or the environment, defaulting to " << local_db_path << endl;
		}
#endif		
		// if it still wasn't found then throw an exception!
		if( local_db_path.length() == 0 )
			throw "Unable to read mpiBLAST local storage path from either .ncbirc or the MPIBLAST_LOCAL environment variable.\nThe mpiBLAST configuration in .ncbirc should look like:\n[mpiBLAST]\nShared=/path/to/shared/storage\nLocal=/path/to/local/storage";
		
		/*
		  This first check is no longer needed, since shared storage isn't necessary when using --copy-via=mpi

		  It could be brought back in, for those cases where shared storage is necessary.
		// make sure shared and local paths exist
		if( !checkCreateDir( shared_db_path ) ){
			string* error_msg = new string("The shared storage directory " + shared_db_path + " does not exist\n");
			throw error_msg->c_str();
		}
		*/

		if( !checkCreateMultiDir( local_db_path ) ){
			string* error_msg = new string("The local storage directory " + local_db_path + " does not exist\n");
			throw error_msg->c_str();
		}
		
		if( !checkDirWritePerms( local_db_path ) ){
			string* error_msg = new string("The local storage directory " + local_db_path + " appears to be write-protected\n");
			throw error_msg->c_str();
		}

		// ensure that the last character is a path separator
		if( shared_db_path.substr( shared_db_path.size() - 1 ) != PATH_SEPARATOR )
			shared_db_path += PATH_SEPARATOR;
		if( local_db_path.substr( local_db_path.size() - 1 ) != PATH_SEPARATOR )
			local_db_path += PATH_SEPARATOR;

		return;
	}

	// Do we need to make a local copy of the config file?
  
  	if( !isRemotePath( filename ) ){
  		// The file is local or accessible via nfs
		config_filename = filename;
  	}
  	else{
  		// Make a copy of the file in the current directory
		// Define the local name by (a) stripping off the hostname and
		// (b) appending a random six character string to the end.
	  	string::size_type pos = filename.find(":", 0);
		pos ++;
		config_filename = filename.substr(pos, filename.length() - pos) + "XXXXXX";
  
  		try{
			getTempFileName(config_filename);
		}
		catch(const char*){
			cerr << "Error creating temp config file:" << endl;
			throw "Error creating temp config file";
		}
	
		copyFile(filename, config_filename, COPY_VIA_XCP);
  	}
  
  
  	// read config file
  	ifstream config_file( config_filename.c_str() );
	
  	if( !config_file.is_open() ){
    		cerr << "Error opening configuration file: " << config_filename << endl;
		cerr << "Try specifying the configuration file location with the --config-file option\n";
    		throw __FILE__ "(MpiBlastConfig): Unable to open config file";
  	}
  
  	// read shared path name
  	getline( config_file, shared_db_path );
	
  	// read local path name
  	getline( config_file, local_db_path );
	
  	// read BLAST binaries path name
  	getline( config_file, blast_path );
	
  	config_file.close();
  
  	// Delete the temporary local file name if we created one
  	if( isRemotePath( filename ) ){
  		unlink( config_filename.c_str() );
  	}

	// ensure that the last character is a path separator
	if( shared_db_path.substr( shared_db_path.size() - 1 ) != PATH_SEPARATOR )
		shared_db_path += PATH_SEPARATOR;
	if( local_db_path.substr( local_db_path.size() - 1 ) != PATH_SEPARATOR )
		local_db_path += PATH_SEPARATOR;
}

MpiBlastConfig::MpiBlastConfig( const MpiBlastConfig& mbc )
{
  *this = mbc;
}

MpiBlastConfig& MpiBlastConfig::operator=( const MpiBlastConfig& mbc )
{
  config_filename = mbc.config_filename;
  local_db_path = mbc.local_db_path;
  shared_db_path = mbc.shared_db_path;
  blast_path = mbc.blast_path;
  return *this;
}

MpiBlastConfig::MpiBlastConfig( const string& local, const string& shared )
{
	local_db_path = local;
	shared_db_path = shared;

	char* tmp_path;
	if( local_db_path.length() != 0 && 
		!checkCreateMultiDir( local_db_path ) && 
		(tmp_path = getenv( "MPIBLAST_LOCAL" )) != NULL ){
		cerr << "Warning:  local db path given as " << local_db_path << " in .ncbirc is unreadable.  Trying MPIBLAST_LOCAL environment variable.\n";
		local_db_path = tmp_path;
	}
	if( local_db_path.length() == 0 && 
		(tmp_path = getenv( "MPIBLAST_LOCAL" )) != NULL ){
		local_db_path = tmp_path;
	}

	if( local_db_path.length() != 0 && 
		!checkCreateMultiDir( local_db_path ) && 
		(tmp_path = getenv( "TMPDIR" )) != NULL ){
		cerr << "Warning:  local db path given as " << local_db_path << " in .ncbirc is unreadable.  Trying $TMPDIR/mpiBLAST_local_db.\n";
		local_db_path = tmp_path + PATH_SEPARATOR + "mpiBLAST_local_db";
	}

	// try setting it to $TMPDIR/mpiBLAST_local_db and print a warning
	if( local_db_path.length() == 0 && 
		(tmp_path = getenv( "TMPDIR" )) != NULL ){
		local_db_path = tmp_path + PATH_SEPARATOR + "mpiBLAST_local_db";
		cerr << "Warning: mpiBLAST local storage path unspecified in either .ncbirc or the environment, defaulting to " << local_db_path << endl;
	}
	// finally default to /tmp/mpiBLAST_local_db on unix
	if( local_db_path.length() == 0  ){
		local_db_path = "/tmp/mpiBLAST_local_db";
		cerr << "Warning: mpiBLAST local storage path unspecified in either .ncbirc or the environment, defaulting to " << local_db_path << endl;
	}
	// if it still wasn't found then throw an exception!
	if( local_db_path.length() == 0 )
		throw "Unable to read mpiBLAST local storage path from either .ncbirc or the MPIBLAST_LOCAL environment variable.\nThe mpiBLAST configuration in .ncbirc should look like:\n[mpiBLAST]\nShared=/path/to/shared/storage\nLocal=/path/to/local/storage";
	
	if( !checkCreateMultiDir( local_db_path ) ){
		string* error_msg = new string("The local storage directory " + local_db_path + " does not exist\n");
		throw error_msg->c_str();
	}
	
	if( !checkDirWritePerms( local_db_path ) ){
		string* error_msg = new string("The local storage directory " + local_db_path + " appears to be write-protected\n");
		throw error_msg->c_str();
	}

	// ensure that the last character is a path separator
	if( shared_db_path.substr( shared_db_path.size() - 1 ) != PATH_SEPARATOR )
		shared_db_path += PATH_SEPARATOR;
	if( local_db_path.substr( local_db_path.size() - 1 ) != PATH_SEPARATOR )
		local_db_path += PATH_SEPARATOR;

	return;
}

string MpiBlastConfig::defaultConfigFileName()
{
	static string file_name;

// check for definition of the default config file location
#ifdef DEFAULT_CONFIG_FILENAME
	file_name = DEFAULT_CONFIG_FILENAME;
#endif

#ifdef WIN32
	if( file_name == "" ){
		const char* home_path = getenv( "USERPROFILE" );
		if( home_path != NULL ){
			file_name = home_path;
			file_name += PATH_SEPARATOR;
			file_name += ".mpiblastrc";
			if( doesFileExist( file_name ) )
				return file_name;
		}
		// try the windows system directory 
		const char* winnt_path = getenv( "windir" );
		if( winnt_path != NULL ){
			file_name = winnt_path;
			file_name += PATH_SEPARATOR;
			file_name += "mpiblast.ini";
		}else
			file_name = "mpiblast.conf";
	}
#else
	if( file_name == "" ){
		const char* home_path = getenv( "HOME" );
		if( home_path != NULL ){
			file_name = home_path;
			file_name += PATH_SEPARATOR;
			file_name += ".mpiblastrc";
			if( doesFileExist( file_name ) )
				return file_name;
		}
		// use the mpiblast.conf in the installation path
		file_name = INSTALL_PREFIX;
		file_name += PATH_SEPARATOR + "etc" + PATH_SEPARATOR + "mpiblast.conf";
	}
#endif
	return file_name;
}

const string& MpiBlastConfig::TempDirBaseName() {
	static string base_name = "mpiblast_temp";
	return base_name;
}

void MpiBlastConfig::createLocalTempDir() 
{
	string tmpdir_basename = local_db_path + TempDirBaseName() + "XXXXXX";
	char* tmpdirname = new char[ tmpdir_basename.length() + 1 ];
	strncpy( tmpdirname, tmpdir_basename.c_str(), tmpdir_basename.length() + 1 );
#ifndef WIN32
	char* success = mkdtemp( tmpdirname );
#else
	char* success = _mktemp( tmpdirname );
	// windows doesn't actually create the directory...
	if( success != NULL )
		checkCreateDir( tmpdirname );
#endif
	if( success == NULL ){
		perror( "Error creating temp directory" );
		throw tmpdirname;
	}
	local_db_path = tmpdirname + PATH_SEPARATOR;
	delete[] tmpdirname;
}


