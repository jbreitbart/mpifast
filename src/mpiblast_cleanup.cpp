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
#include <mpi.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "mpiblast_config.hpp"
#include "mpiblast_util.hpp"
#include "fragment_list.hpp"
#include "file_util.hpp"

using namespace std;

int main( int argc, char* argv[]){
	//Default to .ncbirc, else 1st arg is mpiblast config-file
	//Does NOT work with 1st arg == /path/to/.ncbirc
	string config_filename = ".ncbirc";
	if( argc != 1 && argc != 2 ){
		cerr << "Usage: " << argv[0] << " [/path/to/mpiblast.conf]" << endl;
		return -1;
	}
	if( argc == 2 )
		config_filename = argv[1];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MpiBlastConfig config( config_filename );

		// don't delete shared storage!!
		if (config.localPath() == config.sharedPath()){
			if( my_rank == 0 ){
				cerr << "Shared and local storage are identical.\n";
				cerr << "This program will not delete files from shared storage\n";
			}
			MPI_Finalize();
			return -1;
		}


	try {
		vector< string > arg_vector;
		arg_vector.push_back( argv[0] );
		initNCBI( arg_vector );

		deleteMpiBlastFiles( config.localPath() );
	}
	catch (char const* message){
		cerr << message << endl;
	}
	MPI_Finalize();
	return 0;
}
