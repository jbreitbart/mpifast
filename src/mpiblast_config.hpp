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
#ifndef __MPI_BLAST_CONFIG
#define __MPI_BLAST_CONFIG

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

/**
 * API for reading the mpiBLAST configuration file.
 */
class MpiBlastConfig {
public:
	MpiBlastConfig(){};

	/**
	 * If no configuration file name is specified (default), this
	 * first attempts to read configuration values from an [mpiBLAST] section
	 * of a .ncbirc configuration file
	 * If a given a configuration file name, it reads the old mpiBLAST format
	 * configuration file which specifies:
	 * /path/to/shared/storage
	 * /path/to/local/storage
	 *
	 * Throws a (const char*) exception if something fails.
	 */
	MpiBlastConfig( const std::string& filename );
	MpiBlastConfig( const MpiBlastConfig& mbc );
	MpiBlastConfig& operator=( const MpiBlastConfig& mbc );
	MpiBlastConfig( const std::string& local, const std::string& shared);

	/** 
	 * Returns the path to shared storage
	 * The returned path always ends with a path separator character 
	 */
	const std::string& sharedPath() const { return shared_db_path; }

	/**
	 * Returns the path to local storage 
	 * The returned path always ends with a path separator character 
	 */
	const std::string& localPath() const { return local_db_path; }

	/** Returns the path to the NCBI BLAST binaries */
	const std::string& blastPath() const { return blast_path; }

	/** Returns the default path to the configuration file */
	static std::string defaultConfigFileName();

	/** Returns the base name for temporary directories created by mpiBLAST 
	  * add XXXXXX to this before getting a directory name with mkdtemp() 
	  */
	static const std::string& TempDirBaseName();

	/** creates a temp directory under the current local storage directory
	    and sets the local storage to the new directory */
	void createLocalTempDir();

private:
  std::string config_filename;	/**< The path to the config file */
  std::string local_db_path;	/**< The path to the local database storage */
  std::string shared_db_path;	/**< The path to the shared database storage */
  std::string blast_path;		/**< The path to the NCBI BLAST binaries (obsolete) */
};

#endif
