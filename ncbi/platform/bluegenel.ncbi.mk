NCBI_DEFAULT_LCL = plx
NCBI_MAKE_SHELL = /bin/sh
NCBI_AR= /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-ar
NCBI_CC = /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-gcc -I/bgl/BlueLight/ppcfloor/blrts-gnu/include
NCBI_CFLAGS1 = -c -DARCH_BGL -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
NCBI_LDFLAGS1 = -O3
NCBI_OPTFLAG =  -O3
NCBI_BIN_MASTER = 
NCBI_BIN_COPY = 
NCBI_INCDIR = 
NCBI_LIBDIR = 
NCBI_ALTLIB = 
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/X11R6/include -L/usr/X11R6/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXt -lX11 -lXext
#warning! If you have only dynamic version of Motif or Lesstif
#you should delete -Wl,-Bstatic sentence from the next line:
NCBI_DISTVIBLIBS = -L/usr/X11R6/lib -Wl,-Bstatic -lXm -Wl,-Bdynamic -lXmu -lXt -lX11 -lXext -lXp
NCBI_OTHERLIBS = -lm -L/bgl/BlueLight/ppcfloor/bglsys/lib -lrts.rts -ldevices.rts -lc -lnss_files -lnss_dns -lresolv 
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = 
#NCBI_OTHERLIBS_MT = $(NCBI_MT_OTHERLIBS) -lm
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

# uncomment OPENGL_TARGETS to build OpenGL apps; do not change
# OPENGL_NCBI_LIBS! However, may need to set
# OPENGL_INCLUDE and OPENGL_LIBS to suit local environment
# OPENGL_TARGETS = Cn3D
OPENGL_NCBI_LIBS = LIB400=libvibrantOGL.a LIB3000=libncbicn3dOGL.a
OPENGL_INCLUDE =
OPENGL_LIBS = -L/usr/X11R6/lib -lGL -lGLU
NCBI_OGLLIBS = -L/usr/X11R6/lib -lGL -lGLU

# uncomment (and change appropriately) these lines to build PNG
# output support into Cn3D (OpenGL version only)
#LIBPNG_DIR = /usr/lib
#ZLIB_DIR = /usr/lib

NCBI_LBSM_SRC = ncbi_lbsmd_stub.c
NCBI_LBSM_OBJ = ncbi_lbsmd_stub.o
