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
#ifndef __MPI_BLAST_FRAGMENT_LIST
#define __MPI_BLAST_FRAGMENT_LIST

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include "mpiblast_config.hpp"

/**
 * API for reading and modifying a database fragment list file.
 * The fragment list file tracks the database fragments currently on
 * local storage.
 */
class FragmentListFile {
public:
	FragmentListFile(){};
	/**
	 * Opens and reads a fragment list file, throws a (const char *) exception
	 * if something fails.
	 */
  	FragmentListFile( const std::string& filename, const MpiBlastConfig& config, const std::string& database_name, const std::string& db_type, bool in_mem = false );
 	FragmentListFile( const FragmentListFile& dbsf )
	{
		*this = dbsf;
	};
	
  	FragmentListFile& operator=( const FragmentListFile& mbc );
  	~FragmentListFile()
	{
		// Do nothing
	};
  
	/**
	* Populate the fragment list from an input stream containing a fragment list
	*/
	void parseStream( std::istream& frag_stream );

#ifdef WIN32
	/**
	* Open a fragment list file (r/w) and read its contents,
	* return the file handle.  Leaves the file pointer at the end
	* of the file.  Returns INVALID_HANDLE_VALUE on error.
	*/
	void* openAndReadWin32();
#else
	/**
	* Open a fragment list file (r/w) and read its contents,
	* return the open file descriptor.  Leaves the file pointer at the end
	* of the file.  Returns -1 on error.
	*/
	int openAndReadPOSIX();
#endif

	/** 
  	 * Adds the fragment ID given in <code>fragmentI</code>
  	 * to the fragment list file
  	 */
  	void addFragment( int fragmentI );
	
	/**
	 * Set the database timestamps for the master copy of the database
	 * @param	all_dates	a series of null-terminated datestamp strings, 
	 *						one for every fragment in the FULL database
	 * @param	total_datelen	The total length of space allocated for all_dates
	 */
	void setTimestamps( char* all_dates, int total_datelen );

	/**
	 * Remove any out-of-date fragments from fragment_list and fragment_set
	 */
	void checkTimestamps();

	/**
	 * Send the fragment list to another node.  Uses sendIntVec()
	 */
	void SendList( int dest );
	
	/**
	 * Returns the number of fragments currently on local storage
	 */
	inline int fragmentCount()
	{
		return (int)(fragment_list.size());
	};
	
	/**
	 * Returns true if the fragment list contains fragmentI
	 */
	bool contains( int fragmentI );
private:
  	std::vector< int > fragment_list;
 	std::set< int > fragment_set;
  	std::string frag_filename;
	
	// these variables are for db updates,
	// ideally they should all be tucked inside
	// the MpiBlastConfig and there should be a
	// single per-process instance of MpiBlastConfig
	MpiBlastConfig config;
	std::string database_name;
	std::string db_type;
	std::vector< std::string > timestamps;
	bool is_in_mem;
};


/**
 * File name extensions for sequence database files
 */
const std::vector< std::string >& FragmentExtensions( const std::string& db_type);

extern bool lock_mbf;	/**< When true use file-locking whenever accessing a fragment list */

#endif
