
.PHONY: ncbi

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = LICENSE blast.cgi WWWBlastWrap.pl msvc_projects/mpiblast.sln \
	msvc_projects/mpiblast.vcproj msvc_projects/mpiformatdb.vcproj \
	TODO mpiblast.nsi mpiblast.dox ncbi_Oct2004_evalue.patch ncbi_Jun2005_evalue.patch \
	ech_10k.fas
SUBDIRS = src 

ncbi:
	if test "@ARCH@" == "normal"; then \
		./ncbi/make/makedis.csh ; \
	else \
		./ncbi/make/makedis.csh @ARCH@ ; \
	fi
