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
#include "db_spec.hpp"
#include "mpiblast.hpp"
using namespace std;

DbSpecFile::DbSpecFile( const string& filename )
{
  ifstream spec_file( filename.c_str() );

  if( !spec_file.is_open() ){
    cerr << "(" << my_rank << ") Error opening: " << filename << endl;
    throw __FILE__ "(DbSpecFile): Unable to open file";
  }

  string db_name;
  string nc;

  getline( spec_file, db_name );
  string db_size_line;
  getline( spec_file, db_size_line );
  stringstream db_sss( db_size_line );
  db_sss >> db_size;
  spec_file >> fragment_count;
  spec_file.close();
}

DbSpecFile::DbSpecFile( const DbSpecFile& dbsf )
{
  *this = dbsf;
}

DbSpecFile& DbSpecFile::operator=( const DbSpecFile& mbc )
{

  db_size = mbc.db_size;
  fragment_count = mbc.fragment_count;
  return *this;
}

void DbSpecFile::write( ostream& os, const string& db_name, uint64 db_size, uint fragment_count ) {
	os << db_name << endl;
	os << db_size << endl;
	os << fragment_count << endl;
}
