########################################################################
#       Generic command-line makefile for NCBI applications
#
#   This assumes the following variables are set in a file called "ncbi.mk"
#     which is included by this one. In this case it resides in a directory
#     set by the environment variable "NCBI"
#
#    NCBI_CC = default compiler to use
#    NCBI_CFLAGS1 = default compiler switches (not optimization)
#    NCBI_LDFLAGS1 = default compiler switches when linking (not optimization)
#    NCBI_OPTFLAG = default optimization (-O)
#    NCBI_INCDIR = default toolkit include paths
#    NCBI_LIBDIR = default toolkit library paths
#    NCBI_ALTLIB  = debugging toolkit library paths
#    NCBI_VIBFLAG = additional include paths/defines for Vibrant apps
#    NCBI_VIBLIBS  = additional libraries needed to link Vibrant apps
#    NCBI_OTHERLIBS = additional libraries needed to link the system
#
#   This enables platform independent makefiles for ncbi toolkit apps. You
#   can change any default behavior by setting the variable on the command
#   line.
#
#   Use gcc complier.   "make NCBI_CC=gcc" or "make CC=gcc"
#   Change optimization to debug
#                "make NCBI_OPTFLAG="-g" NCBI_LIBPATH=$NCBI_ALTLIB"
#                or
#                "make OPTFLAG="-g" LIBPATH=$NCBI_ALTLIB"
#
#   You can also change the behavior permanently in the makefile by changing
#    the assignment of the internal variable to the one from the include
#    file:
#
#   Use gcc compiler.
#    Change   CC=$(NCBI_CC)
#    To       CC=gcc
#    in the makefile itself
#
#   Make debugging version
#    OPTFLAG=-g
#    LIBPATH=$(NCBI_ALTDIR)
#
#   You need to specify the EXE and SRC symbols
#
#######################################################################
#
# default flags for compiling and loading
#
#######################################################################

include $(NCBI)/ncbi.mk

CC = $(NCBI_CC)
CFLAGS1 = $(NCBI_CFLAGS1)
LDFLAGS1 = $(NCBI_LDFLAGS1)
# OPTFLAG = $(NCBI_OPTFLAG)
OPTFLAG = -g
OTHERLIBS = $(NCBI_OTHERLIBS)
VIBLIBS= $(NCBI_VIBLIBS)
VIBFLAG= $(NCBI_VIBFLAG)
#NCBI_INCDIR = ../include
#NCBI_LIBDIR = ../lib
INCPATH = $(NCBI_INCDIR)
LIBPATH = $(NCBI_LIBDIR)
DEF =
##### put your source filename here:
SRC = 
##### put executable name here:
EXE = mmdborgcpy pdbeast tax_index orgmatch


CFLAGS = $(CFLAGS1) $(OPTFLAG) -I$(INCPATH) $(VIBFLAG) $(DEF)
LDFLAGS = -I$(INCPATH) $(OPTFLAG) -L$(LIBPATH) $(VIBFLAG) $(DEF)

LIB1 = -lncbi
LIB2 = -lncbiobj
LIB3 = -lncbicdr
LIB4 = -lvibrant
LIB5 = -lncbiacc
LIB6 = -lnetcli
LIB7 = -lnetentr
LIB8 = -lncbiNacc
LIB9 = -lncbiCacc
# LIB10 is reserved for NCSA socket library
LIB10 =
LIB11 = -lncbimla
LIB12 = -lncbitax
LIB13 = -lncbiid1
#LIB14 = -lncbibls0
LIB15 = -lnetcliE
LIB16 = -lnetcliS
LIB17 = -lnetcliES
LIB18 = -lncbibls1
LIB19 = -lncbispel
# LIB20 is for the NCBI desktop utilities
LIB20 = -lncbidesk
LIB21 = -lncbibls2
LIB22 = -lncbimmdb
LIB23 = -lncbitool
LIB24 = -lncbisugg
LIB25 = -lncbiwww
LIB26 = -lncbitax1
LIB27 = -lncbimsc1
LIB29 = -lncbitxc2
LIB30 = -lncbicn3d
LIB31 = -lvibnet
LIB32 = -lncbiwww2
LIB33 = -lblastapi
LIB34 = -lblast
LIB35 = -lblastcompadj

## To clean out the directory without removing make
##

## Implicit actions
##
## if need a .o, compile the .c
##

.c.o :
	$(CC) $(CFLAGS) $<

###############################################################################
##
## This is the makefile for the following products, used for MMDB processing:
##
##	getmuid2 - retrieves Medline ids and inserts them into (ASN1) biostrucs
##
##	getgid2 - retrieves sequence ids and inserts them into (ASN1) biostrucs
##
##	struc2seq - creates sequence entries from (ASN1) biostrucs already
##		containing Medline ids
##
##	makeset - creates a single binary file containing a list of all
##		sequence entries
##
##	pdbset - creates a singel binary file containing a set of biostrucs
##
###############################################################################

all : $(EXE)

pdbeast: interface.c pdbeast.c
	$(CC) -o pdbeast -DX11 -I./ -I/~/ncbi/build $(LDFLAGS) \
	interface.c pdbeast.c  \
	$(LIB29) $(LIB20) $(LIB33) $(LIB34) $(LIB23) $(LIB4) $(LIB22) $(LIB8) $(LIB7) $(LIB6) $(LIB3) \
 	$(LIB2) $(LIB1) $(LIB35) $(VIBLIBS) $(OTHERLIBS)

mmdborgcpy: mmdborgcpy.c
	$(CC) -o mmdborgcpy -DX11 -I./ -I/~/ncbi/build $(LDFLAGS) \
	mmdborgcpy.c \
	$(LIB29) $(LIB22) $(LIB8) $(LIB7) $(LIB6) $(LIB3) \
 	$(LIB2) $(LIB1) $(LIB35)  $(OTHERLIBS)

tax_index: tax_index.c
	$(CC) -o tax_index -DX11 -I./ -I/~/ncbi/build $(LDFLAGS) \
	tax_index.c  \
	$(LIB29) $(LIB22) $(LIB8) $(LIB7) $(LIB6) $(LIB3) \
		  $(LIB2) $(LIB1) $(LIB35) $(OTHERLIBS)

orgmatch: orgmatch.c
	$(CC) -o orgmatch -DX11 -I./ -I/~/ncbi/build $(LDFLAGS) orgmatch.c  \
	$(LIB29) $(LIB22) $(LIB8) $(LIB7) $(LIB6) $(LIB3) \
	$(LIB2) $(LIB1) $(LIB35)  $(OTHERLIBS)

