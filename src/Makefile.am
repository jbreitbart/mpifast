AM_CPPFLAGS = -DINSTALL_PREFIX="\"@prefix@\"" @COMPILE_FLAGS@ -D_FILE_OFFSET_BITS=64 \
	-D_LARGEFILE64_SOURCE -DMPIBLAST_PIO \
	-I$(top_srcdir)/include/ -I@CONFIG_NCBIDIR@/include/ $(MPE_INCS) \
	-I@CONFIG_NCBIDIR@/corelib/ 
AM_CXXFLAGS = -DMPICH_IGNORE_CXX_SEEK -DINSTALL_PREFIX="\"@prefix@\"" \
	@COMPILE_FLAGS@ -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE \
	-I$(top_srcdir)/include/ -I@CONFIG_NCBIDIR@/include/ $(MPE_INCS) \
	-I@CONFIG_NCBIDIR@/corelib/ 

if MPEXXX
MPE_INCS = -I@CONFIG_MPEDIR@/include
MPE_LIBS = -L@CONFIG_MPEDIR@/lib -llmpe -lmpe
else
MPE_INCS =
MPE_LIBS =
endif

bin_PROGRAMS = mpiblast mpiformatdb mpiblast_cleanup

mpiblastsrc = blastjob.cpp db_spec.cpp db_spec.hpp file_util.cpp file_util.hpp \
	fragment_list.cpp fragment_list.hpp mpiblast_util.cpp mpiblast_util.hpp \
	mpiblast.cpp blast_hooks.h blastjob.hpp mpiblast.hpp mpiblast_config.cpp \
	mpiblast_config.hpp blast_wrapper.cpp blast_wrapper.hpp mpiblast_querydata.cpp \
	mpiblast_querydata.hpp mpiblast_scheduler.cpp mpiblast_scheduler.hpp \
	mpiblast_writer.cpp mpiblast_writer.hpp mpiblast_types.h distributed_bioseq.h \
	distributed_bioseq.cpp pio_intercept.cpp pio_intercept.h virtual_files.hpp \
	virtual_files.cpp query_manager.cpp query_manager.hpp mpi_util.cpp mpi_util.h
mpiblast_SOURCES = $(mpiblastsrc)

mpiblast_LDADD = @GETOPT_LONG_SYSTEM@ @CONFIG_NCBIDIR@/build/ncbithr.o \
	@CONFIG_NCBIDIR@/build/blast_hooks.o -lmyncbitool -lblastcompadj -lmyncbiobj \
	-lmyncbi -L. -L@CONFIG_NCBIDIR@/lib/ $(MPE_LIBS) -lm @LINK_FLAGS@
mpiblast_DEPENDENCIES = @GETOPT_LONG_SYSTEM@
mpiblast_CPPFLAGS = -DUSING_MPI -DMPICH_IGNORE_CXX_SEEK

mpiformatdb_SOURCES = mpiformatdb.cpp mpiblast_config.cpp mpiblast_config.hpp \
	file_util.cpp file_util.hpp mpiblast_util.cpp mpiblast_util.hpp \
	fragment_list.cpp fragment_list.hpp db_spec.hpp db_spec.cpp \
	mpi_util.cpp mpi_util.h virtual_files.cpp virtual_files.hpp

mpiformatdb_LDADD = @GETOPT_LONG_SYSTEM@ @CONFIG_NCBIDIR@/build/ncbithr.o \
	@CONFIG_NCBIDIR@/build/mpiblast_formatdb.o -lncbitool -lblastcompadj -lncbiobj \
	-lncbi -L. -L@CONFIG_NCBIDIR@/lib/ $(MPE_LIBS) -lm @LINK_FLAGS@
mpiformatdb_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mpiblast_cleanupsrc = mpiblast_cleanup.cpp mpiblast_config.cpp file_util.cpp \
	mpiblast_util.cpp fragment_list.cpp mpi_util.cpp mpi_util.h virtual_files.cpp \
	virtual_files.hpp
mpiblast_cleanup_SOURCES = $(mpiblast_cleanupsrc)
mpiblast_cleanup_LDADD = @CONFIG_NCBIDIR@/build/ncbithr.o  -lncbitool -lblastcompadj \
	-lncbiobj -lncbi -L. -L@CONFIG_NCBIDIR@/lib/ $(MPE_LIBS) -lm @LINK_FLAGS@

mpiformatdb_INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_NCBIDIR@/include/ $(MPE_INCS) 
