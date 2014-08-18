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
#ifndef __VIRTUAL_FILES_HPP__
#define __VIRTUAL_FILES_HPP__

#ifdef __cplusplus
#include "ncbi.h"
#include "readdb.h"

#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;

enum VFTYPE {	MAPREAD,				//mapped read file
				MAPWRITE,				//mapped write file
				MAPVIRTUAL				//mapped virtual file
};

class VFile {	// File in memory
public:
	VFile(const std::string& afile_name, void* afp, VFTYPE atype, Int4 alength) :
		vf_name(afile_name), vf_fp(afp), vf_type(atype), vf_len(alength) {}
	~VFile();
		
	void *GetPointer() { return vf_fp; }
	Int4 GetLen() { return vf_len; }
private:
	string vf_name;	// file name to be simulated
	void* vf_fp;	// file pointer, could be NlmMFILEPtr or pointer to a mem block
	Int4 vf_len;	// length of the file, only for VF_MEM type
	VFTYPE vf_type;
};

typedef std::map<std::string, VFile *> VFMAP; 
typedef VFMAP::iterator VFMAPIT;

class VFM {	// Virtual files manager
public:
	static VFM* Instance();
	NlmMFILEPtr InsertMapFile(const string& afilename, int mode, Int4 aSize=0);
	bool Erase(const std::string& afilename);
	void* GetFilePtr(const std::string& afilename);
	Int8 GetFileLength(const std::string& afilename);
	void PrintMapFiles();

	int InsertDBFragment(int frag_id);
	int LoadVirtualFile(const string& src, const string& dest);
	int DumpVirtualFile(const string& dest);
	set <int>& GetFragsList(); // only support single database for now
	void Destroy();
private:
	VFMAP vfm_map;
	set <int> _local_frags;
	static VFM* _instance;
	VFM() {};
	VFM(const VFM&);
	VFM& operator=(VFM &);
};

//Wrapper function for calling from C program
extern "C" {
#endif

Int8 VFMGetFileLength(char* afilename);
extern int use_virtual_frags; /* 0-false; 1-true */

#ifdef __cplusplus
}
#endif

#endif //__VIRTUAL_FILES_H__
