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
#include "virtual_files.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_writer.hpp"
#include "file_util.hpp"

extern "C" {
extern NlmMFILEPtr LIBCALL NlmOpenMFILEForWrite (CharPtr name, Nlm_Int8 length);
extern Nlm_MemMapPtr Nlm_MemMapInitForWrite(const Nlm_Char PNTR name, Nlm_Int8 length);
extern NlmMFILEPtr LIBCALL NlmOpenMFILEVirtual (CharPtr name, Nlm_Int8 length);
extern NlmMFILEPtr LIBCALL NlmCloseMFILEVirtual (NlmMFILEPtr mfp);
NlmMFILEPtr VFMFileOpen(char* afilename);
}

#if defined(ARCH_BGL) || defined(ARCH_BGP)
int use_virtual_frags = 1;
#else
int use_virtual_frags = 0;
#endif

VFile::~VFile() {
	if(vf_fp!=NULL) {
		switch(vf_type) {
		case MAPREAD:
			NlmCloseMFILE((NlmMFILEPtr)vf_fp);
			break;
		//case MAPWRITE:
		//	NlmCloseMFILE((NlmMFILEPtr)vf_fp);
		//	break;
		case MAPVIRTUAL:
			NlmCloseMFILEVirtual((NlmMFILEPtr)vf_fp);
			break;
		}
	}
}

VFM* VFM::_instance=NULL;

VFM* VFM::Instance() {
	if(_instance==NULL) {
		_instance=new VFM();
	}
	return _instance;
}

NlmMFILEPtr VFM::InsertMapFile(const string& afilename, int amode, Int4 asize) {
	// check duplicate file
	VFMAPIT it=vfm_map.find(afilename);
	if(it!=vfm_map.end()) {
		throw __FILE__ "Cannot insert duplicate file";
	}

	if(debug_msg) {
		LOG_MSG << "Inserting virtul file: " << afilename <<  endl;
	}
	
	NlmMFILEPtr mfp=NULL;
	VFTYPE type;
	switch(amode) {
	case 0: //read only mapped file
		mfp=NlmOpenMFILE((char *)afilename.c_str());
		type=MAPREAD;
		break;
	case 1: //mapped file for write
		//mfp=NlmOpenMFILEForWrite((char *)afilename.c_str(), asize);
		//type=MAPWRITE;
		//break;
	case 2: //virtual mapped file, no actual disk file
		mfp=NlmOpenMFILEVirtual((char *)afilename.c_str(), asize);
		type=MAPVIRTUAL;
		break;
	default:
		break;
	}
	
	if(mfp!=NULL) {	
		vfm_map.insert(VFMAP::value_type(afilename, new VFile(afilename, mfp, type, asize)));
	} else {
		throw __FILE__ " cannot create virtual file";
	}
	return mfp;
}

void VFM::Destroy() {
	VFMAPIT vit;
	for(vit=vfm_map.begin(); vit!=vfm_map.end(); vit++) {
		delete vit->second;
	}

	vfm_map.clear();
}

set <int>& VFM::GetFragsList() {
	return _local_frags;
}

bool VFM::Erase(const std::string& afilename) {
	VFMAPIT it=vfm_map.find(afilename);	

	if(it!=vfm_map.end()) {
		delete it->second;
	}

	vfm_map.erase(afilename);
	return true;
}

void* VFM::GetFilePtr(const std::string& afilename) {
	VFMAPIT it=vfm_map.find(afilename);

	if(it!=vfm_map.end()) {
		return (it->second)->GetPointer();
	} else {
		return NULL;
	}
}

Int8 VFM::GetFileLength(const std::string& afilename) {
	NlmMFILEPtr mfp=(NlmMFILEPtr)GetFilePtr(afilename);
	if(mfp==NULL) {
		return 0;
	}
	
	return mfp->mem_mapp->file_size;
}

void VFM::PrintMapFiles() {
	int count=0;
	for(VFMAPIT it=vfm_map.begin(); it!=vfm_map.end(); it++) {
		//cout<<"virtual file "<<count++<<" "<<it->first<<endl;
	}	
}

int VFM::InsertDBFragment(int frag_id) {
	_local_frags.insert(frag_id);
}

int VFM::LoadVirtualFile(const string& src, const string& dest) {
	
	uint64 file_len = statFileSize(src.c_str());

	if(file_len == 0) {
		// throw __FILE__ "VFM::LoadVirtualFile - Cannot load empty file";
		return 0;
	}

	NlmMFILEPtr mfp = InsertMapFile(dest, MAPVIRTUAL, file_len);
	
	// copy file content into the virtual file buffer 
	if(io_function == MPI_IO_FUNC) {
		MPI_Status status;
		MPI_File mpifh;
		MPI_File_open(MPI_COMM_SELF, (char *)(src.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifh);
		MPI_File_read_at(mpifh, 0, mfp->mmp_begin, file_len, MPI_BYTE, &status);
		MPI_File_close(&mpifh);
	} else if (io_function == POSIX_IO_FUNC) {
		int fd;
#ifdef O_LARGEFILE
		fd = open(src.c_str(), O_RDONLY|O_LARGEFILE, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#else
		fd = open(src.c_str(), O_RDONLY, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
#endif
		if(fd < 0) {
			return -1;
		}
		
		int read_len = read(fd, mfp->mmp_begin, file_len);
		if(read_len != file_len) {
			return -1;
		}

	} else {
		return -1;
	}

	return 0;
}

int VFM::DumpVirtualFile(const string& filename) {
	VFMAPIT it = vfm_map.find(filename);
	if(it == vfm_map.end()) {
		throw __FILE__ "VFM::DumpVirtualFile -- invalid filename";
	}

	NlmMFILEPtr mfp = (NlmMFILEPtr)(it->second->GetPointer());
	Nlm_Int8 file_length = mfp->mem_mapp->file_size;
	ofstream ofs(filename.c_str()); 
	if(!ofs.is_open()) {
		throw __FILE__ "VFM::DumpVirtualFile -- cannot open file for write";
	}
	
	ofs.write((char*)mfp->mmp_begin, file_length);
	return 0;
}

NlmMFILEPtr VFMFileOpen(char* afilename)
{
	VFM *vfm=VFM::Instance();
	NlmMFILEPtr mfp=(NlmMFILEPtr)vfm->GetFilePtr(string(afilename));
	if(mfp!=NULL) {
		mfp->mmp=mfp->mmp_begin;
	}
	return mfp;
}

Int8 VFMGetFileLength(char* afilename) {
	VFM *vfm=VFM::Instance();
	return vfm->GetFileLength(string(afilename));
}

extern "C" {
NlmMFILEPtr LIBCALL NlmOpenMFILEVirtual (CharPtr name, Nlm_Int8 length) {

    NlmMFILEPtr mfp;
    Nlm_MemMapPtr mem_mapp;

    if (!name || name[0] == '\0')
        return NULL;

    if ((mfp=(NlmMFILEPtr) MemNew(sizeof(NlmMFILE))) == NULL)
        return NULL;

    /* Default is FALSE. */
    mfp->mfile_true = FALSE;

    mfp->mmp_begin = NULL;

    if((mem_mapp = (Nlm_MemMapPtr)Nlm_MemNew(sizeof(Nlm_MemMap))) == NULL)
		return NULL;

	mem_mapp->file_size=length;
	if((mem_mapp->mmp_begin=(Nlm_CharPtr)Nlm_MemNew(length))==NULL) {
		return NULL;
	}

	mfp->mem_mapp=mem_mapp;
	mfp->mmp_madvise_end = mfp->mmp_begin = mfp->mmp = (Uint1Ptr) mfp->mem_mapp->mmp_begin;
	if (mfp->mmp_begin != NULL)
	{
		mfp->mfile_true = TRUE;
		mfp->mmp_end = mfp->mmp_begin + mfp->mem_mapp->file_size;
	}
	
    /* contents have been allocated. */
    mfp->contents_allocated = TRUE;

    return mfp;
}

NlmMFILEPtr LIBCALL NlmCloseMFILEVirtual (NlmMFILEPtr mfp)
{
    if (mfp == NULL)
        return NULL;

    /* Have the contents been allocated, or is this just an attachemnt? */
    if (mfp->contents_allocated)
    {

        if (mfp->mfile_true == TRUE)
        {
            MemFree(mfp->mem_mapp);
        }
    }

    mfp = (NlmMFILEPtr) MemFree(mfp);
    return mfp;

}

//Nlm_MemMapPtr Nlm_MemMapInitForWrite(const Nlm_Char PNTR name, Nlm_Int8 length) {
//  Nlm_MemMapPtr mem_mapp;
//  if (!Nlm_MemMapAvailable()  ||  !name  ||  !*name  ||
//      (mem_mapp = (Nlm_MemMapPtr)Nlm_MemNew(sizeof(Nlm_MemMap))) == NULL)
//    return NULL;
//
//  for (;;) {{ /* (quasi-TRY block) */
//	  mem_mapp->file_size=length;
//
//#ifdef WIN32
//    {{
//      char x_name[MAX_PATH], *str;
//      Nlm_StringNCpy_0(x_name, name, sizeof(x_name));
//      for (str = x_name;  *str;  str++)
//        if (*str == '\\')
//          *str = '/';  /* name of a file-mapping object cannot contain '\' */
//
//      if ( !(mem_mapp->hMap =
//             OpenFileMapping(FILE_MAP_READ, FALSE, x_name)) )
//        { /* If failed to attach to an existing file-mapping object then
//           * create a new one(based on the specified file) */
//          HANDLE hFile= CreateFile(name, GENERIC_READ, FILE_SHARE_READ, NULL,
//                                   OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
//          if (hFile == INVALID_HANDLE_VALUE)
//            break;
//
//          mem_mapp->hMap = CreateFileMapping(hFile, NULL, PAGE_READONLY,
//                                             0, 0, x_name);
//          CloseHandle( hFile );
//          if ( !mem_mapp->hMap )
//            break;
//        }
//
//      if ( !(mem_mapp->mmp_begin =
//             MapViewOfFile(mem_mapp->hMap, FILE_MAP_READ,
//                           0, 0, mem_mapp->file_size)) ) {
//        CloseHandle( mem_mapp->hMap );
//        break;
//      }
//    }}
//
//#elif defined(MMAP_AVAIL)
//    {{  /* UNIX memory mapping. */
//      int fd = open(name, O_RDWR|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
//      if (fd < 0)
//        break;
//      // write last bit;
//      if(lseek(fd, length-1, SEEK_SET)!=-1) {
//		  if(write(fd, "", 1)==1) {
//			  mem_mapp->mmp_begin = mmap(NULL, mem_mapp->file_size, PROT_READ|PROT_WRITE,
//                                 MAP_SHARED, fd, 0);
//          }
//      }
//      close(fd);
//      if ((void*) mem_mapp->mmp_begin == (void*) MAP_FAILED)
//        break;
//    }}
//#endif
//
//    /* Success */
//    return mem_mapp;
//  }}
//
//  /* Error;  cleanup */
//  Nlm_MemFree(mem_mapp);
//  return NULL;
//	
//}

//NlmMFILEPtr LIBCALL NlmOpenMFILEForWrite (CharPtr name, Nlm_Int8 length) {
//	
//    NlmMFILEPtr mfp;
//
//    if (!name || name[0] == '\0')
//        return NULL;
//
//    if ((mfp=(NlmMFILEPtr) MemNew(sizeof(NlmMFILE))) == NULL)
//        return NULL;
//
//    /* Default is FALSE. */
//    mfp->mfile_true = FALSE;
//
//    mfp->mmp_begin = NULL;
//
//    if (Nlm_MemMapAvailable() == TRUE)
//    {     /* IF mem-map fails, open as a regular file. */
//                if((mfp->mem_mapp = Nlm_MemMapInitForWrite(name, length)) != NULL)
//        { /* copy this pointer to where it's convenient. */
//            mfp->mmp_madvise_end = mfp->mmp_begin = mfp->mmp = (Uint1Ptr) mfp->mem_mapp->mmp_begin;
//            if (mfp->mmp_begin != NULL)
//            {
//                mfp->mfile_true = TRUE;
//                mfp->mmp_end = mfp->mmp_begin + mfp->mem_mapp->file_size;
//            }
//        }
//    }
//
//    if (mfp->mmp_begin == NULL)
//    {
//        mfp->fp = FileOpen(name, "wb");
//        if (mfp->fp == NULL)
//        {
//            mfp = (NlmMFILEPtr) MemFree(mfp);
//            return NULL;
//        }
//    }
//
//    /* contents have been allocated. */
//    mfp->contents_allocated = TRUE;
//
//    return mfp;
//}

}
