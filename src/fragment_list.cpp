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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#include "virtual_files.hpp"
#include "fragment_list.hpp"
#include "mpiblast.hpp"
#include "mpiblast_util.hpp"

extern "C" {
#include "blast_hooks.h"
#include "readdb.h"
}

using namespace std;

bool lock_mbf = false;

/*
 * This function uses platform-specific I/O since accesses to the
 * fragment list file must be synchronous (we lock the file)
 */
FragmentListFile::FragmentListFile( const string& filename, const MpiBlastConfig& config, const std::string& database_name, const std::string& db_type, bool in_mem ) : 
config( config ), 
database_name( database_name ),
db_type( db_type )
{
  frag_filename = filename;

  is_in_mem = in_mem;

  if(is_in_mem) {
	  VFM *vfm=VFM::Instance();
	  set <int>& virtual_frags = vfm->GetFragsList();
	  set <int>::iterator sit;
	  for(sit=virtual_frags.begin(); sit!=virtual_frags.end(); sit++) {
		  fragment_list.push_back(*sit);
		  fragment_set.insert(*sit);
	  }

	  return;
  }
  
#ifdef WIN32
	HANDLE fragfile = openAndReadWin32();
	if( fragfile == INVALID_HANDLE_VALUE )
		throw "Unable to access fragment list file";
	CloseHandle( fragfile );
#else
	int fragfile = openAndReadPOSIX();
	if( fragfile == -1 )
		throw "Unable to access fragment list file";
	close( fragfile );
#endif
}

FragmentListFile& FragmentListFile::operator=( const FragmentListFile& mbc ){
  fragment_list = mbc.fragment_list;
  fragment_set = mbc.fragment_set;
  frag_filename = mbc.frag_filename;

	config = mbc.config;
	database_name = mbc.database_name;
	db_type = mbc.db_type ;
	is_in_mem = mbc.is_in_mem;

	return *this;
}

/**
 * Populate the fragment list from an input stream containing a fragment list
 */
void FragmentListFile::parseStream( istream& frag_stream ){
	int fragmentI;
	fragment_list.clear();
	fragment_set.clear();
	while( frag_stream >> fragmentI ){
		fragment_list.push_back( fragmentI );
		fragment_set.insert( fragmentI );
	}
}

#ifdef WIN32

/**
 * Open a fragment list file (r/w) and read its contents,
 * return the file handle.  Leaves the file pointer at the end
 * of the file.  Returns INVALID_HANDLE_VALUE on error.
 */
HANDLE FragmentListFile::openAndReadWin32(){
	HANDLE fragfile;
	// this is really lame:  lock the file when opening
	// and if it's already open by another process, loop until
	// it can be opened
	while(true){
		fragfile = CreateFile( frag_filename.c_str(), 
			GENERIC_READ | GENERIC_WRITE,
			0, 0, OPEN_ALWAYS, 0, 0 );
		if( fragfile == INVALID_HANDLE_VALUE &&
			GetLastError() == ERROR_SHARING_VIOLATION )
						// the file is locked by another process
			Sleep(25);	// sleep 25 milliseconds and try again
		else
			break;	// we either opened it successfully or failed fatally
	}

	if( fragfile == INVALID_HANDLE_VALUE ){
		cerr << "Error opening " << frag_filename << endl;
		return INVALID_HANDLE_VALUE;	// failed
	}
	
	// determine the file size
	// read the whole file into a stringstream
	// parse file using c++ i/o functionality
	DWORD fragfile_size = GetFileSize( fragfile, 0 );
	char* buf = new char[ fragfile_size + 1 ];
	DWORD bytes_read = 0;
	if( !ReadFile( fragfile, buf, fragfile_size, &bytes_read, 0 ) ){
		// an error occurred... do something
		cerr << "Error reading " << frag_filename << endl;
		CloseHandle( fragfile );
		return INVALID_HANDLE_VALUE;
	}
	buf[ fragfile_size ] = 0;

	stringstream ff_sstream( buf );
	delete[] buf;
	parseStream( ff_sstream );
	return fragfile;
}

/**
 * Open a fragment list file and add a fragment if it doesn't
 * already exist in the file
 */
void FragmentListFile::addFragment( int fragmentI ){
	HANDLE fragfile = openAndReadWin32();

	if( fragfile == INVALID_HANDLE_VALUE ){
		// errored out
		return;
	}

	// add the fragment if its not already in the set
	if( fragment_set.find( fragmentI ) == fragment_set.end() ){
		// add it
		char fragno[8]; //space for the frag number
		int nchars = sprintf(fragno, "%03d\n", fragmentI);
		uint64 written;
		WriteFile( fragfile, fragno, nchars, &written, 0 );
	}

	CloseHandle( fragfile );
	fragment_set.insert( fragmentI );
	fragment_list.push_back( fragmentI );
}

#else

/**
 * Open a fragment list file (r/w) and read its contents,
 * return the file descriptor.  Leaves the file pointer at the end
 * of the file.  Returns -1 on error.
 */
int FragmentListFile::openAndReadPOSIX(){
	int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;
	int fragfile = open( frag_filename.c_str(), O_RDWR | O_CREAT, mode );
	if( fragfile < 0 ){
		cerr << "Error opening " << frag_filename << endl;
		perror( "open" );
		return -1;
	}
	if (debug_msg)
		LOG_MSG << "Locking fragment list" << endl;
	int lockval = getlock( fragfile, F_SETLKW, F_WRLCK );
	if ( lockval == -1 ){
		// error locking the file
		cerr << "Error locking file " << frag_filename << endl;
		perror( "getlock" );
		return -1;
	}
	if (debug_msg)
		LOG_MSG << "Locked fragment list" << endl;

	int fragfile_size = statFileSize( frag_filename.c_str() );
	
	char* buf = new char[ fragfile_size + 1 ];
	int bytes_read = read( fragfile, buf, fragfile_size );
	buf[ fragfile_size ] = 0;
	if( bytes_read != fragfile_size ){
		// an error occurred... do something
		cerr << "Error reading from " << frag_filename << endl;
		close( fragfile );
		return -1;
	}

	stringstream ff_sstream( buf );
	delete[] buf;
	parseStream( ff_sstream );
	return fragfile;
}

void FragmentListFile::addFragment( int fragmentI ){

	if(use_virtual_frags) {
		VFM::Instance()->InsertDBFragment( fragmentI );
	} else {
		int fragfile = openAndReadPOSIX();

		if( fragfile == -1 ){
			// errored out
			return;
		}

		// add the fragment if its not already in the set
		if( fragment_set.find( fragmentI ) == fragment_set.end() ){
			// add it
			char fragno[5]; //space for the frag number
			int nchars = sprintf( fragno, "%03d\n", fragmentI );
			write( fragfile, fragno, nchars );
		}

		close( fragfile );
	}	

	fragment_set.insert( fragmentI );
	fragment_list.push_back( fragmentI );
}

#endif

bool FragmentListFile::contains( int fragmentI ){
	return fragment_set.find( fragmentI ) != fragment_set.end();
}

void FragmentListFile::setTimestamps( char* all_dates, int total_datelen ){
	char* tmp_date = all_dates;
	while( tmp_date != all_dates + total_datelen ){
		timestamps.push_back( tmp_date );
		tmp_date += strlen( tmp_date ) + 1;
	}
}

void FragmentListFile::checkTimestamps(){
	// if the remote frag has a different modification time than the local one
	// we pretend the local one doesn't exist
	vector< int > up_to_date_fragments;
	char frag_num_buff[4];
	for(int fragI = 0 ; fragI < (int)fragment_list.size(); fragI++){
		if( fragI >= timestamps.size() )
			break;	// the new database has fewer frags than the old...

		// always use 3 digit fragment ids
		sprintf(frag_num_buff,"%03d",fragment_list[ fragI ]);
		string local_fragname = config.localPath() + database_name + "." + frag_num_buff;
		char* m_dbname = (char*)malloc( local_fragname.length() + 1 );
		strcpy( m_dbname, local_fragname.c_str() );
		ReadDBFILEPtr rdfp = readdb_new_ex2( m_dbname, db_type == "p", READDB_NEW_DO_REPORT, NULL, NULL );
		free( m_dbname );
		if( rdfp == NULL ){
			cerr << "readdb unable to open " << local_fragname << endl;
			continue;
		}
		char* date_str = readdb_get_date( rdfp );

		if( date_str == timestamps[ fragment_list[ fragI ] ] ){
			up_to_date_fragments.push_back( fragment_list[ fragI ] );
		}else if( debug_msg ){
			LOG_MSG << "Out of date fragment:\n";
			LOG_MSG << "Local frag " << local_fragname << " created on " << date_str << endl ;
			LOG_MSG << "Remote frag " << frag_num_buff << " created on " << timestamps[ fragment_list[ fragI ] ] << endl ;
		}	

		rdfp = readdb_destruct( rdfp );
	}
	fragment_list = up_to_date_fragments;
	fragment_set = set< int >( fragment_list.begin(), fragment_list.end() );
}

void FragmentListFile::SendList( int dest )
{

	SendIntVec(group_comm, fragment_list, dest, DB_FRAGMENTS_TAG);
}
	
vector<string> makeFragmentExtensions( string db_type);
const vector<string>& FragmentExtensions( const string& db_type)
{
  static vector<string> n_extensions = makeFragmentExtensions( "n" );
  static vector<string> p_extensions = makeFragmentExtensions( "p" );
  if( db_type == "n" ){
  	return n_extensions;
  }
  return p_extensions;
}

/**
 * Because DB "indices" are required, all 7 extensions have become
 * mandatory.
 */
vector<string> makeFragmentExtensions(string db_type)
{
	vector< string > extensions;
	
	extensions.push_back( "." + db_type + "hr" );
	extensions.push_back( "." + db_type + "in" );
	extensions.push_back( "." + db_type + "sq" );
	extensions.push_back( "." + db_type + "nd" );
	extensions.push_back( "." + db_type + "ni" );
	extensions.push_back( "." + db_type + "sd" );
	extensions.push_back( "." + db_type + "si" );
	
	return extensions;
}

