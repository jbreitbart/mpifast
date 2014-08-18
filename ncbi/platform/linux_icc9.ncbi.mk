#
# $Id: linux_icc9.ncbi.mk,v 1.9 2009/02/04 14:24:03 lavr Exp $
#
# ICC 9.0 with optimization options for Pentium 4 processor

NCBI_DEFAULT_LCL = lnx
NCBI_MAKE_SHELL = /bin/sh
NCBI_AR=xiar
# Defining _GCC_NEXT_LIMITS_H ensures that <limits.h> chaining doesn't
# stop short, as can otherwise happen (at least with ICC 10).
NCBI_CC = icc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_GCC_NEXT_LIMITS_H -D_GNU_SOURCE
NCBI_CFLAGS1 = -c
NCBI_LDFLAGS1 = -fast
NCBI_OPTFLAG = -fast
#
#for Pentium4 you can try:
#NCBI_LDFLAGS1 = -i-static -O3 -tpp7 -mcpu=pentium4 -xW -march=pentium4
#NCBI_OPTFLAG = -O3 -tpp7 -mcpu=pentium4 -xW -march=pentium4
#
NCBI_BIN_MASTER = /home/coremake/ncbi/bin
NCBI_BIN_COPY = /home/coremake/ncbi/bin
NCBI_INCDIR = /home/coremake/ncbi/include
NCBI_LIBDIR = /home/coremake/ncbi/lib
NCBI_ALTLIB = /home/coremake/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/X11R6/include -L/usr/X11R6/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXt -lX11 -lXext
#warning! If you have only dynamic version of Motif or Lesstif
#you should delete -Wl,-Bstatic sentence from the next line:
NCBI_DISTVIBLIBS = -L/usr/X11R6/lib -Wl,-Bstatic -lXm -Wl,-Bdynamic -lXmu -lXt -lX11 -lXext -lXp
NCBI_OTHERLIBS = -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = -lpthread
NCBI_OTHERLIBS_MT = $(NCBI_MT_OTHERLIBS) -lm
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

# uncomment OPENGL_TARGETS to build OpenGL apps; do not change
# OPENGL_NCBI_LIBS! However, may need to set
# OPENGL_INCLUDE and OPENGL_LIBS to suit local environment
# OPENGL_TARGETS = Cn3D
OPENGL_NCBI_LIBS = LIB400=libvibrantOGL.a LIB3000=libncbicn3dOGL.a
OPENGL_INCLUDE = -I/usr/X11R6/include
OPENGL_LIBS = -L/usr/X11R6/lib -lGL -lGLU
NCBI_OGLLIBS = -L/usr/X11R6/lib -lGL -lGLU

# uncomment (and change appropriately) these lines to build PNG
# output support into Cn3D (OpenGL version only)
#LIBPNG_DIR = /home/paul/Programs/libpng
#ZLIB_DIR = /home/paul/Programs/zlib

NCBI_LBSM_SRC = ncbi_lbsmd_stub.c
NCBI_LBSM_OBJ = ncbi_lbsmd_stub.o
