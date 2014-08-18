/*   wprint.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  wprint.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov
*
* Version Creation Date:   7/15/95
*
* $Revision: 6.74 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/


#include <ncbi.h>
#include <objsset.h>
#include <prtutil.h>
#include <seqport.h>
#include <sequtil.h>
#include <asn2ffp.h>
#include <ffprint.h>

#ifdef ENABLE_ENTREZ
#include <accentr.h>
#endif

#define MAX_WWWBUF 328
#define MAX_WWWLOC 2

static CharPtr lim_str [MAX_WWWLOC] = {">","<"};
static CharPtr www_lim_str [MAX_WWWLOC] = {"&gt;","&lt;"};

static Char link_ff[MAX_WWWBUF];
static Char link_muid[MAX_WWWBUF];
static Char link_seq[MAX_WWWBUF];
static Char link_ace[MAX_WWWBUF];
static Char link_tax[MAX_WWWBUF];
static Char link_code[MAX_WWWBUF];
static Char link_fly[MAX_WWWBUF];
static Char link_cog[MAX_WWWBUF];
static Char link_sgd[MAX_WWWBUF];
static Char link_gdb[MAX_WWWBUF];
static Char link_ck[MAX_WWWBUF];
static Char link_rice[MAX_WWWBUF];
static Char link_sp[MAX_WWWBUF];
static Char link_pir[MAX_WWWBUF];
static Char link_pdb[MAX_WWWBUF];
static Char link_gdb_map[MAX_WWWBUF];
static Char link_UniSTS[MAX_WWWBUF];
static Char link_dbSTS[MAX_WWWBUF];
static Char link_dbEST[MAX_WWWBUF];
static Char link_omim[MAX_WWWBUF];
static Char link_locus[MAX_WWWBUF];
static Char link_snp[MAX_WWWBUF];
static Char link_ratmap[MAX_WWWBUF];
static Char link_rgd[MAX_WWWBUF];
static Char link_mgd[MAX_WWWBUF];
static Char link_fly_fban[MAX_WWWBUF];
static Char link_fly_fbgn[MAX_WWWBUF];
static Char link_cdd[MAX_WWWBUF];
static Char link_niaest[MAX_WWWBUF];
static Char link_worm_sequence[MAX_WWWBUF];
static Char link_worm_locus[MAX_WWWBUF];
static Char link_imgt[MAX_WWWBUF];
static Char link_ifo[MAX_WWWBUF];
static Char link_jcm[MAX_WWWBUF];
static Char link_isfinder[MAX_WWWBUF];
static Char link_gabi[MAX_WWWBUF];
static Char link_fantom[MAX_WWWBUF];

#define DEF_LINK_FF  "/cgi-bin/Entrez/getfeat?"

#define DEF_LINK_ACE  "http://www.ncbi.nlm.nih.gov/AceView/hs.cgi?"
#define DEF_LINK_SEQ  "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"
#define DEF_LINK_MUID  "/entrez/utils/qmap.cgi?"

#define DEF_LINK_TAX  "/Taxonomy/Browser/wwwtax.cgi?"
#define DEF_LINK_CODE "/Taxonomy/Utils/wprintgc.cgi?"

#define DEF_LINK_COG "http://www.ncbi.nlm.nih.gov/cgi-bin/COG/palox?"
#define DEF_LINK_MGD "http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:"
#define DEF_LINK_FBGN "http://flybase.bio.indiana.edu/.bin/fbidq.html?"
#define DEF_LINK_FBAN "http://www.fruitfly.org/cgi-bin/annot/fban?"

/*
#define DEF_LINK_FLY "/cgi-bin/Entrez/referer?http://morgan.harvard.edu/htbin-post/gene.script%3f"
*/
#define DEF_LINK_FLY "http://flybase.bio.indiana.edu/.bin/fbidq.html?"
/*
#define DEF_LINK_SGD "/cgi-bin/Entrez/referer?http://genome-www.stanford.edu/cgi-bin/dbrun/SacchDB%3ffind+SGDID+"
*/

#define DEF_LINK_SGD "/cgi-bin/Entrez/referer?http://genome-www4.stanford.edu/cgi-bin/SGD/locus.pl?locus="

#define DEF_LINK_GDB "http://gdbwww.gdb.org/gdb-bin/genera/genera/hgd/DBObject/GDB:"
#define DEF_LINK_CK "http://flybane.berkeley.edu/cgi-bin/cDNA/CK_clone.pl?db=CK&dbid="
/*
#define DEF_LINK_RICE "http://genome.cornell.edu:80/cgi-bin/webace?db=ricegenes&class=Probe&object="
*/
#define DEF_LINK_RICE "http://ars-genome.cornell.edu/cgi-bin/WebAce/webace?db=ricegenes&class=Marker&object="
#define DEF_LINK_SP "/cgi-bin/Entrez/referer?http://expasy.hcuge.ch/cgi-bin/sprot-search-ac%3f"
#define DEF_LINK_PDB "/cgi-bin/Entrez/referer?http://expasy.hcuge.ch/cgi-bin/get-pdb-entry%3f"
#define DEF_LINK_GDB_PREFIX "http://gdbwww.gdb.org/gdb-bin/gdb/browser/bin/locq?ACTION=query&cyto="

#define DEF_LINK_GDB_SUFFIX "&match=Inclusive&order=Locus+Location"

#define DEF_LINK_UniSTS "http://www.ncbi.nlm.nih.gov/genome/sts/sts.cgi?uid="
#define DEF_LINK_dbSTS "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"
#define DEF_LINK_dbEST "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"
#define DEF_LINK_OMIM "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id="
#define DEF_LINK_LOCUS "http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l="
#define DEF_LINK_SNP "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs="
#define DEF_LINK_RATMAP "http://ratmap.gen.gu.se/action.lasso?-database=RATMAPfmPro&-layout=Detail&-response=/RM/Detail+Format.html&-search&-recid="
#define DEF_LINK_RGD "http://rgd.mcw.edu/query/query.cgi?id="

#define DEF_LINK_CDD "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="
#define DEF_LINK_NIAEST "http://lgsun.grc.nia.nih.gov/cgi-bin/pro3?sname1="
#define DEF_LINK_WORM_SEQUENCE "http://www.wormbase.org/db/seq/sequence?name="
#define DEF_LINK_WORM_LOCUS "http://www.wormbase.org/db/gene/locus?name="
#define DEF_LINK_IMGT "http://imgt.cines.fr:8104/cgi-bin/IMGTlect.jv?query=202+"
#define DEF_LINK_IFO "http://www.ifo.or.jp/index_e.html"
#define DEF_LINK_JCM "http://www.jcm.riken.go.jp/cgi-bin/jcm/jcm_number?JCM="
#define DEF_LINK_ISFINDER "http://www-is.biotoul.fr/scripts/is/is_spec.idc?name="
#define DEF_LINK_GABI "https://gabi.rzpd.de/cgi-bin-protected/GreenCards.pl.cgi?Mode=ShowBioObject&BioObjectName="
#define DEF_LINK_FANTOM "http://fantom.gsc.riken.go.jp/db/view/main.cgi?masterid="

/* now other data bases are linked to Entrez. may be changed later 
static char *link_epd = 
"/htbin-post/Entrez/query";
static char *link_tfd = 
"/htbin-post/Entrez/query";
*/
static Boolean www = FALSE;
static Pointer data = NULL;
static HeadTailProc head = NULL;
static HeadTailProc tail = NULL;

NLM_EXTERN CharPtr PrintDate(NCBI_DatePtr date)
{
	CharPtr retval = NULL;
	char month[4], year[5], day[3];
	char result[15];

	if ( date -> data[0] == 0){
/*---string---*/
		if (StringICmp(date -> str,"Not given") != 0){
			retval = StringSave(date -> str);
		}
	} else {
/*---standard---*/
		if (date->data[1] && date->data[2] && date->data[3]) {
			if ((int) (date -> data[1]) < 30) {
				sprintf(year, "%4d", (int) (date -> data[1] + 2000));
			} else {
				sprintf(year, "%4d", (int) (date -> data[1] + 1900));
			}
			sprintf(day, "%d", (int) (date -> data[3] ));
			StringCpy(month, NCBI_months[date->data[2] -1 ]);
			sprintf(result,"%s %s, %s.\n", month, day, year);
			retval = StringSave(result);
		}
	}
	return retval;
}

void PrintXrefButton PROTO((FILE *fp, SeqEntryPtr sep));

NLM_EXTERN void LIBCALL init_www(void)
{
	www = TRUE;
	GetAppParam("NCBI", "WWWENTREZ", "LINK_FF", DEF_LINK_FF, 
		link_ff, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_MUID", DEF_LINK_MUID, 
		link_muid, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_ACE", DEF_LINK_ACE, 
		link_ace, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, 
		link_seq, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_TAX", DEF_LINK_TAX, 
		link_tax, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_CODE", DEF_LINK_CODE, 
		link_code, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_FLY", DEF_LINK_FLY, 
		link_fly, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_COG", DEF_LINK_COG, 
		link_cog, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_SGD", DEF_LINK_SGD, 
		link_sgd, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_SGD", DEF_LINK_GDB, 
		link_gdb, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_CK", DEF_LINK_CK, 
		link_ck, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_RICE", DEF_LINK_RICE, 
		link_rice, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_SP", DEF_LINK_SP, 
		link_sp, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_PDB", DEF_LINK_PDB, 
		link_pdb, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_OMIM", DEF_LINK_OMIM, 
		link_omim, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_UniSTS", DEF_LINK_UniSTS, 
		link_UniSTS, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_dbSTS", DEF_LINK_dbSTS, 
		link_dbSTS, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_dbEST", DEF_LINK_dbEST, 
		link_dbEST, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_LOCUS", DEF_LINK_LOCUS, 
		link_locus, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_SNP", DEF_LINK_SNP, 
		link_snp, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_RATMAP", DEF_LINK_RATMAP, 
		link_ratmap, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_RGD", DEF_LINK_RGD, 
		link_rgd, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_MGD", DEF_LINK_MGD, 
		link_mgd, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_FBGN", DEF_LINK_FBGN, 
		link_fly_fbgn, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_FBAN", DEF_LINK_FBAN, 
		link_fly_fban, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_CDD", DEF_LINK_CDD, 
		link_cdd, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_NIAEST", DEF_LINK_NIAEST, 
		link_niaest, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_WORM_SEQUENCE", DEF_LINK_WORM_SEQUENCE, 
		link_worm_sequence, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_WORM_LOCUS", DEF_LINK_WORM_LOCUS, 
		link_worm_locus, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_IMGT", DEF_LINK_IMGT, 
		link_imgt, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_IFO", DEF_LINK_IFO, 
		link_ifo, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_JCM", DEF_LINK_JCM, 
		link_jcm, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_ISFINDER", DEF_LINK_ISFINDER,
		link_isfinder, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_GABI", DEF_LINK_GABI,
		link_gabi, MAX_WWWBUF);
	GetAppParam("NCBI", "WWWENTREZ", "LINK_FANTOM", DEF_LINK_FANTOM,
		link_fantom, MAX_WWWBUF);

}

NLM_EXTERN void LIBCALL fini_www(void)
{
	www = FALSE;
}

NLM_EXTERN void LIBCALL head_tail_ff(Pointer mydata, HeadTailProc headfun, HeadTailProc tailfun)
{
	data = mydata;
	head = headfun;
	tail = tailfun;
}

NLM_EXTERN void LIBCALL head_www(FILE *fp, SeqEntryPtr sep)
{

	if (head != NULL) {
		head(data, fp);
		return;
	}
	if (!www) {
		return;
	}
	fprintf(fp, "Content-type: text/html\n\n");
	fprintf(fp, "<HTML>\n");
	fprintf(fp, "<HEAD><TITLE>%s", "GenBank entry");
	fprintf(fp, "</TITLE></HEAD>\n");
	fprintf(fp, "<BODY>\n");
	fprintf(fp, "<hr>\n");
	fprintf(fp, "<pre>");
} 

NLM_EXTERN void LIBCALL tail_www(FILE *fp)
{
	if (tail != NULL) {
		tail(data, fp);
		return;
	}
	if (!www) {
		return;
	}
	fprintf(fp, "</pre>\n");
	fprintf(fp, "<hr>\n");
	fprintf(fp, "</BODY>\n");
	fprintf(fp, "</HTML>\n");
}

NLM_EXTERN Boolean LIBCALL get_www(void)
{
	return www;
}

/**************************************************************************
*	source is coming from GBBlock.source so we have to find real 
*	organism name for the hot link 
***************************************************************************/
NLM_EXTERN Boolean LIBCALL www_source(CharPtr orgname, OrgRefPtr orp)
{
	Int2	l, ll, lll;
	CharPtr	s;
	
	if (www && orp) {
		l = StringLen(link_tax);
		lll = StringLen(orgname);
		ll = StringLen("<a href=%sname=%s>");
		s = (CharPtr)MemNew(l+ ll + lll);
		sprintf(s, "<a href=%sname=%s>", link_tax, 
					orp->taxname?orp->taxname:orp->common);
		AddLink(s);
		MemFree(s);
		ff_AddString(orgname);
		AddLink("</a>");
	} else {
		ff_AddString(orgname);
	}
		return TRUE;

}			
NLM_EXTERN Boolean LIBCALL www_organism(CharPtr orgname, Int4 id)
{
	Int2	l, ll, lll;
	CharPtr	s, linkname, ss;
	
	if (www) {
		l = StringLen(link_tax);
		lll = StringLen(orgname);
		linkname = StringSave(orgname);
		for (ss = linkname; *ss != '\0'; ss++) {
			if (*ss == ' ') {
				*ss = '+';
			}
		}
		if (id != -1) {
			ll = StringLen("<a href=%sid=%d>");
		} else {
			ll = StringLen("<a href=%sname=%s>");
		}
		s = (CharPtr)MemNew(l+ ll + lll);
		if (id != -1) {
			sprintf(s, "<a href=%sid=%d>", link_tax, id);
		} else {
			sprintf(s, "<a href=%sname=%s>", link_tax, linkname);
		}
		AddLink(s);
		MemFree(s);
		ff_AddString(orgname);
		AddLink("</a>");
		MemFree(linkname);
	} else {
		ff_AddString(orgname);
	}
		return TRUE;
}			

NLM_EXTERN Boolean LIBCALL www_taxid(CharPtr orgname, Int4 id)
{
	Int2	l, ll, lll;
	CharPtr	s, linkname, ss;
	
	if (www) {
		l = StringLen(link_tax);
		lll = StringLen(orgname);
		linkname = StringSave(orgname);
		for (ss = linkname; *ss != '\0'; ss++) {
			if (*ss == ' ') {
				*ss = '+';
			}
		}
		if (id != -1) {
			ll = StringLen("<a href=%sid=%d>");
		} else {
			ll = StringLen("<a href=%sname=%s>");
		}
		s = (CharPtr)MemNew(l+ ll + lll);
		if (id != -1) {
			sprintf(s, "<a href=%sid=%d>", link_tax, id);
		} else {
			sprintf(s, "<a href=%sname=%s>", link_tax, linkname);
		}
		AddLink(s);
		MemFree(s);
		ff_AddString(orgname);
	/*	ff_AddInteger(" (%d)", id);*/
		AddLink("</a>");
		MemFree(linkname);
	} else {
		ff_AddString(orgname);
	}
		return TRUE;
}			

NLM_EXTERN Boolean LIBCALL www_featkey(CharPtr key, Int4 gi, Int2 entityID, Uint4 itemID)
{
	Int2	l, ll;
	CharPtr	s;
	
	if (www) {
		l = StringLen(link_ff);
		ll = StringLen("<a href=%sgi=%ld&id=%d&entity=%d>");
		s = (CharPtr)MemNew(l+ ll + 3*7);
		sprintf(s, "<a href=%sgi=%ld&id=%d&entity=%d>", 
							link_ff, (Int4) gi, itemID, entityID);
		AddLink(s);
		MemFree(s);
		ff_AddString(key);
		AddLink("</a>");
	} else {
		ff_AddString(key);
	}
		return TRUE;
}			

NLM_EXTERN Boolean LIBCALL www_gcode(CharPtr gcode)
{
	Int2	l, ll, gc, lll = 1;
	CharPtr	s;
	
	if (www) {
		gc = atoi(gcode);
		if (gc >= 10)
			lll = 2;
		l = StringLen(link_code);
			ll = StringLen("<a href=%smode=c#SG%d>");
		s = (CharPtr)MemNew(l+ ll + lll);
			sprintf(s, "<a href=%smode=c#SG%d>", link_code, gc);
		AddLink(s);
		MemFree(s);
		ff_AddInteger("%d", gc);
		AddLink("</a>");
	} else {
		ff_AddString(gcode);
	}
		return TRUE;
}			

NLM_EXTERN Boolean LIBCALL www_muid(Int4 muid)
{
	Int2	l, ll;
	CharPtr	s;
	
	if (www) {
		l = StringLen(link_muid);
		ll = StringLen("<a href=%suid=%ld&form=6&db=m&Dopt=r>");
		s = (CharPtr)MemNew(l+ ll + 10);
		sprintf(s, "<a href=%suid=%ld&form=6&db=m&Dopt=r>", 
												link_muid, (Int4) muid);
		AddLink(s);
		MemFree(s);
		ff_AddInteger("%ld", (long) muid);
		AddLink("</a>");
	} else {
		ff_AddInteger("%ld", (long) muid);
	}
		return TRUE;
}			

NLM_EXTERN Boolean LIBCALL www_extra_acc(CharPtr acc, Boolean ncbi)
{
	Int2	l, ll;
	CharPtr	s, prefix;
	
	if (www && !ncbi) {
		l = StringLen(link_seq);
		prefix = "<a href=%sval=%s>"; 
		ll = StringLen(prefix); 
		s = (CharPtr)MemNew(l+ ll + 10);
		sprintf(s, prefix, link_seq, acc);
		AddLink(s);
		MemFree(s);
		ff_AddString( acc);
		AddLink("</a>");
	} else {
		ff_AddString( acc);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_genpept_gi(CharPtr str)
{
	Int2	l, ll;
	Int4 	gi;
	CharPtr	s, prefix;
	
	if(www) {
			l = StringLen(link_seq);
			prefix = "<a href=%sval=%ld>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 10);
		if (StringNCmp(str, "gi|", 3) == 0) {
			str += 3;
			gi = atoi(str);
			ff_AddString("gi|");
			sprintf(s, prefix, link_seq, gi);
			AddLink(s);
			MemFree(s);
			ff_AddInteger("%ld", (long) gi);
			AddLink("</a>");
			for(; IS_DIGIT(*str); str++);
			ff_AddString(str);
		}
	} else {
		ff_AddString(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_dbsource(CharPtr str, Boolean first, Uint1 choice)
{
	Int2	l, ll;
	CharPtr	s, prefix=NULL, p, text, link=NULL;
	
	if(www) {
		if (choice == SEQID_PIR /*|| choice == SEQID_SWISSPROT*/) {
			link = link_seq;
			prefix = "<a href=%sval=%s>";
		} else if (choice == SEQID_PDB || choice == SEQID_PRF) {
			link = link_seq;
			prefix = "<a href=%sval=%s>";
		} else if (choice == SEQID_EMBL || choice == SEQID_GENBANK || 
				choice == SEQID_DDBJ || choice == SEQID_GIBBSQ || 
				choice == SEQID_GIBBMT || choice == SEQID_GI || 
				choice == SEQID_GIIM || choice == SEQID_OTHER ||
				choice == SEQID_TPG || choice == SEQID_TPE || choice == SEQID_TPD)  {
			link = link_seq;
			prefix = "<a href=%sval=%s>";
		} else {
			ff_AddStringWithTildes(str);
			return TRUE;
		}
		l = StringLen(link);
		ll = StringLen(prefix); 
		if ((p = StringStr(str, "accession")) != NULL) {
			p += 9;
			while (*p == ' ') {
				p++;
			}
			text = TextSave(str, p-str);
			if (first == FALSE) {
				ff_AddString(", ");
			}
			ff_AddString(text);
			MemFree(text);
			text = StringSave(p);
			if (text[StringLen(text)-1] == ';') {
				text[StringLen(text)-1] = '\0';
			}
			s = (CharPtr)MemNew(l+ ll + StringLen(text) + 1);
			sprintf(s, prefix, link, text);
			AddLink(s);
			MemFree(s);
			MemFree(text);
			ff_AddString(p);
			AddLink("</a>");
		} else {
			if (first == FALSE) {
				ff_AddString(", ");
			}
			ff_AddString(str);
		}
	} else {
		ff_AddStringWithTildes(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean www_coded_by(CharPtr str)
{
	Int2	l, ll;
	CharPtr	s, ss, prefix, p, text, link;
	
	if(www) {
		link = link_seq;
		prefix = "<a href=%sval=%s>";
		p = StringChr(str, ':');
		while (p != NULL) {
			for (ss = p-1; ss > str && IS_ALPHANUM(*ss); ss--);
			if (ss > str) {
				ss++;
				text = TextSave(str, ss-str);
				ff_AddString(text);
				MemFree(text);
			}
			text = TextSave(ss, p-ss);
			l = StringLen(link);
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + StringLen(text));
			sprintf(s, prefix, link, text);
			AddLink(s);
			MemFree(s);
			ff_AddString(text);
			AddLink("</a>");
			str = p;
			p = StringChr(str+1, ':');
		}
		ff_AddString(str);
	} else {
		ff_AddString(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_map(CharPtr str)
{
	Int2	l, ll, lll, llll;
	CharPtr	s, prefix;

	if (www) {
			l = StringLen(DEF_LINK_GDB_PREFIX);
			prefix = "<a href=%s%s%s>"; 
			ll = StringLen(prefix); 
			lll = StringLen(str); 
			llll = StringLen(DEF_LINK_GDB_SUFFIX); 
			s = (CharPtr)MemNew(l+ ll + lll + llll);
			sprintf(s, prefix, DEF_LINK_GDB_PREFIX, str, DEF_LINK_GDB_SUFFIX);
			AddLink(s);
			MemFree(s);
			ff_AddString(str);
			AddLink("</a>");
	} else {
		ff_AddString(str);
	}
	return TRUE;

}
NLM_EXTERN Boolean LIBCALL www_protein_id(CharPtr str)
{
	Int2	l, ll;
	CharPtr	s, prefix;
	
	if (www) {
		l = StringLen(link_seq);
		prefix = "<a href=%sval=%s>"; 
		ll = StringLen(prefix) + StringLen(str); 
		s = (CharPtr)MemNew(l + ll);
		sprintf(s, prefix, link_seq, str);
		AddLink(s);
		MemFree(s);
		ff_AddString(str);
		AddLink("</a>");
	} else {
		ff_AddString(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_db_xref(CharPtr str)
{
	Int2	l, ll;
	Int4 	gi;
	CharPtr	s, prefix, ss, p, pp;
	Boolean nothing = TRUE;
	Char id[10];
	Int2 id1, id2;
	
	if (www) {
		while ((p= StringStr(str, "FLYBASE:")) != NULL) {
			nothing = FALSE;
			p += StringLen("FlyBase:");
			l = StringLen(link_fly);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 20);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			for (pp = p; *pp != '\0' && *pp != ';'; pp++);
			ss = (CharPtr)MemNew(pp-p+1);
			StringNCpy(ss, p, pp-p);
			sprintf(s, prefix, link_fly, ss);
			AddLink(s);
			MemFree(s);
			ff_AddString(ss);			
			AddLink("</a>");
			ff_AddChar(';');
			str = pp+1;
			
		}
		if (( p = StringStr(str, "COG:")) != NULL) {
			nothing = FALSE;
			p += StringLen("COG:");
			l = StringLen(link_cog) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_cog, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "FLYBASE:FBa")) != NULL) {
			nothing = FALSE;
			p += StringLen("FLYBASE:");
			l = StringLen(link_fly_fban) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_fly_fban, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "FLYBASE:FBgn")) != NULL) {
			nothing = FALSE;
			p += StringLen("FLYBASE:");
			l = StringLen(link_fly_fbgn) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_fly_fbgn, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "PID:g")) != NULL) {
			nothing = FALSE;
			p += StringLen("PID:g");
			l = StringLen(link_seq);
			prefix = "<a href=%sval=%ld>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 10);
			gi = atoi(p);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			sprintf(s, prefix, link_seq, gi);
			AddLink(s);
			MemFree(s);
			ff_AddInteger("%ld", (long) gi);
			AddLink("</a>");
		} 
/*
		if (( p = StringStr(str, "SWISS-PROT:")) != NULL) {
			nothing = FALSE;
			p += StringLen("SWISS-PROT:");
			l = StringLen(link_seq) + StringLen(p);
			prefix = "<a href=%sval=%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_seq, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		}
*/		
		if (( p = StringStr(str, "UniSTS:")) != NULL) {
			nothing = FALSE;
			p += StringLen("UniSTS:");
			l = StringLen(link_UniSTS) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_UniSTS, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "dbEST:")) != NULL) {
			nothing = FALSE;
			p += StringLen("dbEST:");
			l = StringLen(link_dbEST) + StringLen(p);
			prefix = "<a href=%sval=gnl|dbest|%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_dbEST, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "dbSTS:")) != NULL) {
			nothing = FALSE;
			p += StringLen("dbSTS:");
			l = StringLen(link_dbSTS) + StringLen(p);
			prefix = "<a href=%sval=gnl|dbsts|%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_dbSTS, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "LocusID:")) != NULL) {
			nothing = FALSE;
			p += StringLen("LocusID:");
			l = StringLen(link_locus) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_locus, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "InterimID:")) != NULL) {
			nothing = FALSE;
			p += StringLen("InterimID:");
			l = StringLen(link_locus) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_locus, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "MIM:")) != NULL) {
			nothing = FALSE;
			p += StringLen("MIM:");
			l = StringLen(link_omim) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_omim, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "SGD:")) != NULL) {
			nothing = FALSE;
			p += StringLen("SGD:");
			l = StringLen(link_sgd) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_sgd, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "IMGT/LIGM:")) != NULL) {
			nothing = FALSE;
			p += StringLen("IMGT/LIGM:");
			l = StringLen(link_imgt) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_imgt, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "GDB:")) != NULL) {
			p += StringLen("GDB:");
			if (( pp = StringStr(p, "G00-")) != NULL) {
				pp += StringLen("G00-");
				id1 = atoi(pp);
				for (; *pp != '\0' && *pp != '-'; pp++);
				id2 = atoi(pp+1);
				sprintf(id, "%d%d", id1, id2);
				l = StringLen(link_gdb) + StringLen(p);
				prefix = "<a href=%s%s>"; 
				ll = StringLen(prefix); 
				s = (CharPtr)MemNew(l + ll);
				ss = (CharPtr)MemNew(p-str+1);
				StringNCpy(ss, str, p-str);
				ff_AddString(ss);
				MemFree(ss);
				sprintf(s, prefix, link_gdb, id);
				AddLink(s);
				MemFree(s);
				ff_AddString(p);
				AddLink("</a>");
				nothing = FALSE;
			} else if (IS_DIGIT(*p)) {
				l = StringLen(link_gdb) + StringLen(p);
				prefix = "<a href=%s%s>"; 
				ll = StringLen(prefix); 
				s = (CharPtr)MemNew(l + ll);
				ss = (CharPtr)MemNew(p-str+1);
				StringNCpy(ss, str, p-str);
				ff_AddString(ss);
				MemFree(ss);
				sprintf(s, prefix, link_gdb, p);
				AddLink(s);
				MemFree(s);
				ff_AddString(p);
				AddLink("</a>");
				nothing = FALSE;
			}
		} 
		if (( p = StringStr(str, "CK:")) != NULL) {
			nothing = FALSE;
			p += StringLen("CK:");
			l = StringLen(link_ck) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			sprintf(s, prefix, link_ck, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "RiceGenes:")) != NULL) {
			nothing = FALSE;
			p += StringLen("RiceGenes:");
			l = StringLen(link_rice) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			sprintf(s, prefix, link_rice, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "dbSNP:")) != NULL) {
			nothing = FALSE;
			p += StringLen("dbSNP:");
			l = StringLen(link_snp) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_snp, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "RATMAP:")) != NULL) {
			nothing = FALSE;
			p += StringLen("RATMAP:");
			l = StringLen(link_ratmap) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_ratmap, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "RGD:")) != NULL) {
			nothing = FALSE;
			p += StringLen("RGD:");
			l = StringLen(link_rgd) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_rgd, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "MGD:")) != NULL) {
			nothing = FALSE;
			p += StringLen("MGD:");
			l = StringLen(link_mgd) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			if (StringNICmp (p, "MGI:", 4) == 0) {
				p += 4;
			}
			sprintf(s, prefix, link_mgd, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "CDD:")) != NULL) {
			nothing = FALSE;
			p += StringLen("CDD:");
			l = StringLen(link_cdd) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_cdd, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "niaEST:")) != NULL) {
			nothing = FALSE;
			p += StringLen("niaEST:");
			l = StringLen(link_niaest) + StringLen(p);
			prefix = "<a href=%s%s&val=1>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_niaest, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "WormBase:")) != NULL) {
			CharPtr worm_link = NULL;
			nothing = FALSE;
			p += StringLen("WormBase:");
			if (StringStr (str, "unc") != NULL) {
				worm_link = link_worm_locus;
				l = StringLen(worm_link) + StringLen(p);
				prefix = "<a href=%s%s;class=Locus>"; 
			} else {
				worm_link = link_worm_sequence;
				l = StringLen(worm_link) + StringLen(p);
				prefix = "<a href=%s%s;class=Sequence>"; 
			}
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, worm_link, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "IFO:")) != NULL) {
			nothing = FALSE;
			p += StringLen("IFO:");
			l = StringLen(link_ifo) + StringLen(p);
			prefix = "<a href=%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_ifo);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "JCM:")) != NULL) {
			nothing = FALSE;
			p += StringLen("JCM:");
			l = StringLen(link_jcm) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_jcm, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		}
		if (( p = StringStr(str, "ISFinder:")) != NULL) {
			nothing = FALSE;
			p += StringLen("ISFinder:");
			l = StringLen(link_isfinder) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_isfinder, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		}
		if (( p = StringStr(str, "GABI:")) != NULL) {
			nothing = FALSE;
			p += StringLen("GABI:");
			l = StringLen(link_gabi) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_gabi, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		}
		if (( p = StringStr(str, "FANTOM_DB:")) != NULL) {
			nothing = FALSE;
			p += StringLen("FANTOM_DB:");
			l = StringLen(link_fantom) + StringLen(p);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			sprintf(s, prefix, link_fantom, p);
			AddLink(s);
			MemFree(s);
			ff_AddString(p);
			AddLink("</a>");
		}
		if (nothing) {
			ff_AddString(str);
		}
	} else {
		ff_AddString(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_note_gi(CharPtr str)
{
	Int2	l, ll;
	Int4 	gi;
	CharPtr	s, prefix, ss, p, pp;
	Boolean nothing = TRUE;
	
	if (www) {
		while ((p= StringStr(str, "FlyBase: ")) != NULL) {
			nothing = FALSE;
			p += StringLen("FlyBase: ");
			l = StringLen(link_fly);
			prefix = "<a href=%s%s>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 20);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			while (*p == ' ')
				p++;
			for (pp = p; *pp != '\0' && *pp != ';'; pp++);
			ss = (CharPtr)MemNew(pp-p+1);
			StringNCpy(ss, p, pp-p);
			sprintf(s, prefix, link_fly, ss);
			AddLink(s);
			MemFree(s);
			ff_AddString(ss);			
			AddLink("</a>");
			ff_AddChar(';');
			str = pp+1;
			
		}
		if (( p = StringStr(str, "NCBI gi: ")) != NULL) {
			nothing = FALSE;
			p += StringLen("NCBI gi: ");
			l = StringLen(link_seq);
			prefix = "<a href=%sval=%ld>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 10);
			gi = atoi(p);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			sprintf(s, prefix, link_seq, gi);
			AddLink(s);
			MemFree(s);
			ff_AddInteger("%ld", (long) gi);
			AddLink("</a>");
		} 
		if (( p = StringStr(str, "AceView:")) != NULL) {
			nothing = FALSE;
		/*	p += StringLen("AceView:");*/
			gi = atoi(p + StringLen("AceView:"));
			l = StringLen(link_ace) + StringLen(p);
			prefix = "<a href=%sl=%ld>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-str+1);
			StringNCpy(ss, str, p-str);
			ff_AddString(ss);
			MemFree(ss);
			sprintf(s, prefix, link_ace, gi);
			AddLink(s);
			MemFree(s);
			ff_AddString("AceView");
		/*	ff_AddInteger("%ld", (long) gi);*/
			AddLink("</a>");
			
		} 
		if (nothing) {
			ff_AddString(str);
		}
	} else {
		ff_AddString(str);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL www_xref(CharPtr str, Uint1 xref_class)
{
	Int2	l, ll;
	CharPtr	s, link, prefix, ss;
	
	if (www) {
		ss = (CharPtr)MemNew(StringLen(str) + 1);
		StringCpy(ss, str);
		if (xref_class == 5) {
			link = link_seq;    /*link_sp*/
			prefix = "<a href=%sval=%s>"; 
		} else if (xref_class == 8) {
			link = link_seq;  /*link_epd*/
			prefix = "";
		} else if (xref_class == 10) {
			link = link_seq;  /* link_tfd */
			prefix = "";
		} else if (xref_class == 11) {
			link = link_fly;
			prefix = "<a href=%s%s>";
		}
		if (link && prefix) {
			l = StringLen(link);
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l+ ll + 10);
			if (*(ss + StringLen(ss)  - 1) == '.') 
				*(ss + StringLen(ss)  - 1) = '\0';
				
			sprintf(s, prefix, link, ss);
			AddLink(s);
			MemFree(s);
		}
		ff_AddString( str);
		AddLink("</a>");
		MemFree(ss);
	} else {
		ff_AddString( str);
	}
	return TRUE;
}
NLM_EXTERN Boolean LIBCALL www_xref_button(FILE *fp, CharPtr str, Uint1 xref_class, Uint1 db)
{
	Int2	l, ll;
	CharPtr	s, link, prefix, ss;
	
	if (www) {
		ss = (CharPtr)MemNew(StringLen(str) + 1);
		StringCpy(ss, str);
		if (xref_class == 5) {
			link = link_seq;  /* link_sp */
			prefix = "<a href=%sval=%s>"; 
		} else if (xref_class == 8) {
			link = link_seq;  /* link_epd */
			prefix = "";
		} else if (xref_class == 10) {
			link = link_seq;   /*link_tfd*/
			prefix = "";
		} else if (xref_class == 11) {
			link = link_fly;
			prefix = "<a href=%s%s>";
		} if (xref_class == 255) {    /*for PIR and S-P  use db */
			link = link_seq;
			prefix = "<a href=%sval=%s>"; 
		}
		if (link && prefix) {
			l = StringLen(link);
			ll = StringLen(prefix); 
			if (db == SEQID_PIR) {
				s = (CharPtr)MemNew(l+ ll + 20);
				sprintf(s, prefix, link, str, str);
			} else {
				s = (CharPtr)MemNew(l+ ll + 10);
				if (*(ss + StringLen(ss)  - 1) == '.') {
					*(ss + StringLen(ss)  - 1) = '\0';
				}
				sprintf(s, prefix, link, ss);
			}
			fprintf(fp, "%s<img src=http://www.ncbi.nlm.nih.gov/cgi-bin/gorf/butt?%s align=middle border=0></a>", s, str);
			MemFree(s);
		}
		MemFree(ss);
	}
	return TRUE;
}

NLM_EXTERN Boolean LIBCALL PrintSPBlock (Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	Int2	l, ll;
	CharPtr	link, prefix, s;
	ValNodePtr vnp=NULL, v;
	CharPtr string, acc = NULL;
	ValNodePtr ds_vnp, tvnp;
	DescrStructPtr dsp;
	SPBlockPtr spb;
	SeqIdPtr sid;
	TextSeqIdPtr tid;
	Boolean has_link = FALSE, first;
	DbtagPtr	db;

	tvnp = GatherDescrListByChoice(ajp, gbp, Seq_descr_sp); 
	for (ds_vnp= tvnp; ds_vnp; ds_vnp=ds_vnp->next) {
		dsp = (DescrStructPtr) ds_vnp->data.ptrvalue;
		vnp = dsp->vnp;
		gbp->descr = dsp;
		spb = (SPBlockPtr)vnp->data.ptrvalue;
		MemFree(vnp);
		if (spb->_class == 1) {
			ff_AddString("class: standard.");
			NewContLine();
		} else if (spb->_class == 2) {
			ff_AddString("class: preliminary.");
			NewContLine();
		}
		if (spb->extra_acc) {
			ff_AddString("extra accessions:");
			for (v = spb->extra_acc; v; v= v->next) {
				ff_AddString((CharPtr)v->data.ptrvalue);
				ff_AddString(",");
			}
		}
		if (spb->imeth) {
			ff_AddString("seq starts with Met");
		}
		if (spb->plasnm) {
			ff_AddString("plasmid:");
			for (v = spb->plasnm; v; v= v->next) {
				ff_AddString((CharPtr)v->data.ptrvalue);
				ff_AddString(",");
			}
		}
		if (spb->created) {
			string = PrintDate(spb->created);
			ff_AddString("created: ");
			ff_AddString(string);
		}
		if (spb->sequpd) {
			string = PrintDate(spb->sequpd);
			ff_AddString("sequence updated: ");
			ff_AddString(string);
		}
		if (spb->annotupd) {
			string = PrintDate(spb->annotupd);
			ff_AddString("annotation updated: ");
			ff_AddString(string);
		}
		if (spb->seqref) {
			ff_AddString("xrefs: ");
			first = TRUE;
			for (sid = spb->seqref; sid; sid= sid->next) {
				link = link_seq;
				if (first == FALSE) {
					ff_AddString(", ");
				}
				has_link = TRUE;
				first = FALSE;
				if(sid->choice == SEQID_GI) {
					prefix = "<a href=%sval=%ld>"; 
					ff_AddString("gi: ");
					acc = (CharPtr)MemNew(10);
					sprintf(acc, "%ld", (long) sid->data.intvalue);
				} else {
					prefix = "<a href=%sval=%s>";
					tid = (TextSeqIdPtr)sid->data.ptrvalue;
					acc = tid->accession;
				}
				switch(sid->choice) {
				case SEQID_GENBANK:
					ff_AddString("genbank accession ");
					break; 
				case SEQID_EMBL:
					ff_AddString("embl accession ");
					break; 
				case SEQID_PIR:
					ff_AddString("pir locus ");
					break; 
				case SEQID_SWISSPROT:
					ff_AddString("swissprot accession ");
					break; 
				case SEQID_DDBJ:
					ff_AddString("ddbj accession ");
					break; 
				case SEQID_PRF:
					ff_AddString("prf accession ");
					break; 
				case SEQID_GI:
					ff_AddString("gi: ");
					break; 
				case SEQID_TPG:
					ff_AddString("genbank third party accession ");
					break; 
				case SEQID_TPE:
					ff_AddString("embl third party accession ");
					break; 
				case SEQID_TPD:
					ff_AddString("ddbj third party accession ");
					break; 
				default:
					acc = NULL;
					break; 
				}
				if (www && has_link) {
					if (acc && link && prefix) {
						l = StringLen(link);
						ll = StringLen(prefix); 
						if (sid->choice == SEQID_PIR) {
							s = (CharPtr)MemNew(l+ ll + 40);
							sprintf(s, prefix, link, acc, acc);
						} else if (sid->choice == SEQID_GI) {
							s = (CharPtr)MemNew(l+ ll + 10);
							sprintf(s, prefix, link, sid->data.intvalue);
						} else {
							s = (CharPtr)MemNew(l+ ll + 10);
							sprintf(s, prefix, link, acc);
						}
					}
					AddLink(s);
				/*	MemFree(s); */
					ff_AddString( acc);
					AddLink("</a>");
				} else if (acc) {
					ff_AddString( acc);
				}
			}
		}
		first = TRUE;
		for (vnp = spb->dbref; vnp; vnp=vnp->next) {
			db = (DbtagPtr)vnp->data.ptrvalue;
			has_link = FALSE;
			if (first) {
				NewContLine();
				ff_AddString("xrefs (non-sequence databases): ");
				first = FALSE;
			} else {
				ff_AddString(", ");
			}
			ff_AddString(db->db);
			if (StringCmp(db->db, "MIM") == 0) {
				prefix = "<a href=%s%s>";
				l = StringLen(link_omim);
				ll = StringLen(prefix); 
				s = (CharPtr)MemNew(l+ ll + 10);
				has_link = TRUE;
			} 
			if (db->tag && db->tag->str) {
				ff_AddString(" ");				
				if (www && has_link) {
					sprintf(s, prefix, link_omim, db->tag->str);
					AddLink(s);
					MemFree(s);
					ff_AddString(db->tag->str);
					AddLink("</a>");
				} else {				
					ff_AddString(db->tag->str);
				}
			} else if (db->tag && db->tag->id) {
				ff_AddString(" ");
				if (www && has_link) {
					sprintf(s, "<a href=%s%d>", link_omim, db->tag->id);
					AddLink(s);
					MemFree(s);
					ff_AddInteger("%d", db->tag->id);
					AddLink("</a>");
				} else {				
					ff_AddInteger("%d", db->tag->id);
				}
			}
		}
	}
	return TRUE;
}

NLM_EXTERN CharPtr LIBCALL www_featloc(CharPtr loc)
{
	Int2 i, n, k, l1, l2;
	CharPtr buf, ptr, s, eptr;
	
	if (loc == NULL) {
		return NULL;
	}
	n = StringLen(loc);
	eptr = loc + n - 1;
	for (i = 0, k = 0; i < MAX_WWWLOC; i++) {
		s = loc;
		while (s < eptr) {
			if ((ptr = StringStr(s, lim_str[i])) != NULL) {
				n += StringLen(www_lim_str[i]) - StringLen(lim_str[i]);
				k++;
				s = ptr + StringLen(lim_str[i]);
			} else {
				break;
			}
		}
	}
	if (k == 0) {
		return StringSave(loc);
	}
	buf = (CharPtr)MemNew(n + 1);
	StringCpy(buf, loc);
	for (i = 0; i < MAX_WWWLOC; i++) {
		s = buf;
		while ((ptr = StringStr(s, lim_str[i])) != NULL) {
			l1 = StringLen(www_lim_str[i]);
			l2 = StringLen(lim_str[i]);
			if (l1 != l2) {
				MemMove(ptr+l1, ptr+l2, StringLen(ptr+l2));
			}
			MemCpy(ptr, www_lim_str[i], l1);
			s = ptr + l1;
		}
	}
	
	return buf;
}

static void LitPrintGenome(SeqLitPtr slp)
{
	static Char		val[166];

	if (slp->seq_data != NULL)         /* not a gap */
	{
		if (slp->length == 0)  /* unknown length */
		{
			sprintf(val, "gap(%d)", slp->length);
			ff_AddString(val);
		} else {
/* don't know what to do here */
		}
	} else {                  /* gap length was set */
			sprintf(val, ",gap(%d)", slp->length);
			ff_AddString(val);
	}
}

static void LocPrintGenome(Asn2ffJobPtr ajp, GBEntryPtr gbp, SeqLocPtr slp_head)
{
	SeqLocPtr	slp;
	Boolean		first = TRUE;
	static Char		buf[14], val[166], temp[166];
	SeqIdPtr	sid, newid;
	Int4 		from, to, start, stop;
	BioseqPtr 	bsp = NULL, b = NULL;
	SeqEntryPtr sep = NULL;
	Int4		uid;
	Boolean		is_link = FALSE;
	Int2 p1=0, p2=0;
	static Uint1 fasta_order[NUM_SEQID] = { 
 	33, /* 0 = not set */
	20, /* 1 = local Object-id */
	15,  /* 2 = gibbsq */
	16,  /* 3 = gibbmt */
	30, /* 4 = giim Giimport-id */
	10, /* 5 = genbank */
	10, /* 6 = embl */
	10, /* 7 = pir */
	10, /* 8 = swissprot */
	15,  /* 9 = patent */
	20, /* 10 = other TextSeqId */
	20, /* 11 = general Dbtag */
	255,  /* 12 = gi */
	10, /* 13 = ddbj */
	10, /* 14 = prf */
	12, /* 15 = pdb */
    10,  /* 16 = tpg */
    10,  /* 17 = tpe */
    10,  /* 18 = tpd */
    10,  /* 19 = gpp */
    10   /* 20 = nat */
    };
	
	Int2	l, ll;
	CharPtr	s, prefix;
	
#ifdef ENABLE_ENTREZ
	if ( !EntrezInit("asn2ff", FALSE, &is_network) ) {
		return;
	}
#endif
	for (slp = slp_head; slp; slp = slp->next) {
		is_link = FALSE;
		from = to = 0;
		sid = SeqLocId(slp);
		if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
			start = from = SeqLocStart(slp);
			stop = to = SeqLocStop(slp);
		} else if (slp->choice == SEQLOC_NULL){
			sprintf(val, ",%s", "gap()");
			ff_AddString(val);
			continue;
		} else {
			continue;
		}
		if (sid == NULL) {
			continue;
		}
		if (sid->choice == SEQID_GI) {
			newid = GetSeqIdForGI(sid->data.intvalue);
		} else if (sid->choice == SEQID_GENERAL) {
#ifdef ENABLE_ENTREZ
			if ((uid = GetUniGeneIDForSeqId(sid)) != 0)
				sep = EntrezSeqEntryGet(uid, -1);
			if (sep && IS_Bioseq(sep)){
			 bsp = (BioseqPtr) sep -> data.ptrvalue;
			}
			if (bsp && (bsp->seq_ext_type == 1 || bsp->seq_ext_type == 4)) {
				lcur = lprev = 0;
				for (sslp = slp_head; sslp; sslp = sslp->next) {
					is_link = FALSE;
					lprev = lcur;
					lcur += SeqLocLen(sslp);
					beg = SeqLocStart(sslp);
					end = SeqLocStop(sslp);
					sid = SeqLocId(sslp);
					if (from > lcur || to < lprev) {
						continue;
					}
					if (from > lprev) { /* first */
						start = beg + from - lprev;
					} else {
						start = beg;  /* middle */
					}
					if (to > lcur) {  /* middle */
						stop = end;
					} else {
						stop = beg + to - lprev;
					}
					if (first) {
						first = FALSE;
					} else {
						ff_AddChar(',');
					}
					if (sid->choice == SEQID_GI) {
						newid = GetSeqIdForGI(sid->data.intvalue);
						if (newid == NULL) {
							newid = sid;
						}
					} else {
						newid = sid;
					}
					if (ajp->show_version) {
					SeqIdWrite(SeqIdSelect(newid, fasta_order, NUM_SEQID), buf, PRINTID_TEXTID_ACC_VER, 13);
					} else {
					SeqIdWrite(SeqIdSelect(newid, fasta_order, NUM_SEQID), buf, PRINTID_TEXTID_ACCESSION, 10);
					}
					if (www && newid) {
						l = StringLen(link_seq);
						if (newid->choice == SEQID_GENERAL) {
							uid = GetUniGeneIDForSeqId(newid);
							if (uid != 0) {
								prefix =
								"<a href=%sval=%d>"; 
								ll = StringLen(prefix); 
								s = (CharPtr)MemNew(l+ ll + 10);
								sprintf(s, prefix, link_seq, uid);
								if (SeqLocStrand(slp) == Seq_strand_minus) {
									ff_AddString("complement(");
								}
								ff_AddString(buf);
								p1 = StringLen(buf);
								p2 = 0;
								is_link = TRUE;
							} else {
								if (SeqLocStrand(slp) == Seq_strand_minus) {
									ff_AddString("complement(");
								}
								ff_AddString(buf);
							}
						} else {
							prefix = 
							"<a href=%sval=%s>";
							ll = StringLen(prefix); 
							s = (CharPtr)MemNew(l+ ll + 10);
							sprintf(s, prefix, link_seq, buf);
							if (SeqLocStrand(slp) == Seq_strand_minus) {
								ff_AddString("complement(");
							}
							ff_AddString(buf);
							p1 = StringLen(buf);
							p2 = 0;
							is_link = TRUE;
						}
					} else {
						if (SeqLocStrand(slp) == Seq_strand_minus) {
							ff_AddString("complement(");
						}
						ff_AddString(buf);
					}
					if (SeqLocStrand(slp) == Seq_strand_minus) {
					sprintf(val,":%ld..%ld)",(long) start+1, (long) stop+1 );
					} else {
					sprintf(val,":%ld..%ld",(long) start+1, (long) stop+1 );
					}
					ff_AddString(val);
					p1 += StringLen(val);
					p2 += StringLen(val);
					if (is_link) {
						AddLinkLater(s, p1);
						MemFree(s);
						AddLinkLater("</a>", p2);
					}
				}
				continue;
			}
			newid = sid;
#else
			newid = sid;
#endif
		} else {
			newid = sid;
		}
		if (first) {
			first = FALSE;
		} else {
			ff_AddChar(',');
		}
		if (ajp->show_version) {
			SeqIdWrite(SeqIdSelect(newid, fasta_order, NUM_SEQID),
				buf, PRINTID_TEXTID_ACC_VER, 13);
		} else {
			SeqIdWrite(SeqIdSelect(newid, fasta_order, NUM_SEQID),
				buf, PRINTID_TEXTID_ACCESSION, 10);
		}
		if (www && newid) {
			l = StringLen(link_seq);
			if (newid->choice == SEQID_GENERAL) {
				uid = GetUniGeneIDForSeqId(newid);
				if (uid != 0) {
					prefix = "<a href=%sval=%ld>";
					ll = StringLen(prefix); 
					s = (CharPtr)MemNew(l+ ll + 10);
					sprintf(s, prefix, link_seq, uid);
					if (SeqLocStrand(slp) == Seq_strand_minus) {
						ff_AddString("complement(");
					}
					ff_AddString( buf);
					p1 = StringLen(buf);
					p2 = 0;
					is_link = TRUE;
				} else {
					ff_AddString( buf);
				}
			} else {
				prefix = "<a href=%sval=%s>";
				ll = StringLen(prefix); 
				s = (CharPtr)MemNew(l+ ll + 10);
				sprintf(s, prefix, link_seq, buf);
				if (SeqLocStrand(slp) == Seq_strand_minus) {
					ff_AddString("complement(");
				}
				ff_AddString( buf);
				p1 = StringLen(buf);
				p2 = 0;
				is_link = TRUE;
			} 
		} else {
			if (SeqLocStrand(slp) == Seq_strand_minus) {
				ff_AddString("complement(");
			}
			ff_AddString( buf);
		}
		if (SeqLocStrand(slp) == Seq_strand_minus) {
			sprintf(val,":%ld..%ld)",(long) start+1, (long) stop+1 );
		} else {
			sprintf(val,":%ld..%ld",(long) start+1, (long) stop+1 );
		}
		ff_AddString(val);
		p1 += StringLen(val);
		p2 += StringLen(val);
		if (is_link) {
			AddLinkLater(s, p1);
			MemFree(s);
			AddLinkLater("</a>", p2);
		}
	}
}

void PrintGenome(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	SeqLocPtr	slp_head=NULL;
	Boolean		first = TRUE;
	DeltaSeqPtr dsp;
	SeqLitPtr 	litp;
	
	
#ifdef ENABLE_ENTREZ
	if ( !EntrezInit("asn2ff", FALSE, &is_network) ) {
		return;
	}
#endif
	ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
	if (www && ajp->format == GRAPHIK_FMT) {
		if (gbp->bsp && gbp->bsp->seq_ext == NULL) {
			return;
		}
		ff_AddString("<b>CONTIG</b>&nbsp;&nbsp;&nbsp;");
	} else {
		ff_AddString("CONTIG");
	}
	TabToColumn(13);
	ff_AddString("join(");
	if (gbp->bsp->seq_ext_type == 1) {
		slp_head = (SeqLocPtr) gbp->bsp->seq_ext;
		LocPrintGenome(ajp, gbp, slp_head);
	} else if (gbp->bsp->seq_ext_type == 4) {
		for (dsp = (DeltaSeqPtr) gbp->bsp->seq_ext; dsp; dsp=dsp->next) {
			if (dsp->choice == 1) {  /* SeqLoc */
				slp_head = (SeqLocPtr)(dsp->data.ptrvalue);
				if (first == FALSE) {
					ff_AddString(",");
				}
				LocPrintGenome(ajp, gbp, slp_head);
			} else {
				litp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (litp == NULL) continue;
				LitPrintGenome(litp);
			}
			if (first)
				first = FALSE;
		}
	}
	ff_AddChar(')');
	ff_EndPrint();
	/*PrintTerminator();*/

}

NLM_EXTERN void LIBCALL www_accession (CharPtr string)
{
	CharPtr	s, prefix=NULL, link=NULL;

	if (string == NULL) {
		return;
	}
	if (!www) {
		ff_AddString(string);
	} else {
			link = link_seq;
			prefix = "<a href=%sval=%s>";
			s = (CharPtr)MemNew(StringLen(link_seq)+ StringLen(prefix) + StringLen(string) + 10);
			sprintf(s, prefix, link, string);
			AddLink(s);
			MemFree(s);
			ff_AddString(string);
			AddLink("</a>");
	}
	return;
}
static Boolean iscospa(Char c)
{
	return (isspace(c) || c == ','  || c == ' ') ;
}


NLM_EXTERN void LIBCALL www_PrintComment (CharPtr string, Boolean identifier, Uint1 format)
{
	Int2	lpref, l, ll;
	Int4	gi;
	CharPtr	s, prefix=NULL, p, pp, link=NULL, www_str, acc, ss;
	Boolean isfirst = TRUE;
	
	if (string == NULL) {
		return;
	}
	if (format == EMBL_FMT || format == PSEUDOEMBL_FMT ||
		format == EMBLPEPT_FMT) {
		if (identifier == TRUE)
			PrintXX();
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "CC");
	} else {
		ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
		if (identifier == TRUE) {
			if (format == GRAPHIK_FMT && www) {
				ff_AddString("<BR><b>COMMENT</b>&nbsp;&nbsp;&nbsp;");
			} else {
				ff_AddString("COMMENT");
			}
		}
		TabToColumn(13);
	}
	if (!www) {
		ff_AddStringWithTildes(string);
		ff_EndPrint();
		return;
	} 
	if ((p = StringStr(string, "REFSEQ")) != NULL) {
		acc = TextSave(string, p-string);
		ff_AddString(acc);
		MemFree(acc);
		p += 7;
		AddLink("<a href=http://www.ncbi.nlm.nih.gov/LocusLink/refseq.html>");
		ff_AddString("REFSEQ:");
		AddLink("</a>");
		string = p;
		if ((p = StringStr(string, "was derived from ")) != NULL) {
			p += StringLen("was derived from ");
			s = TextSave(string, p-string);
			ff_AddString(s);
			MemFree(s);
			prefix = "<a href=%sval=%s>";
			lpref = StringLen(link_seq)+ StringLen(prefix);
			for (; isspace(*p); p++) ;
			isfirst = TRUE;
			do {
				if (isfirst) {
					isfirst = FALSE;
				} else {
					ff_AddString(", ");
				}
				for (pp=p; *pp != '\0' && !iscospa(*pp); pp++) ;
				acc = TextSave(p, pp-p);
				if (acc[StringLen(acc)-1] == '.') {
					acc[StringLen(acc)-1] = '\0';
				}
				s = (CharPtr)MemNew(lpref + StringLen(acc)+1);
				sprintf(s, prefix, link_seq, acc);
				AddLink(s);
				MemFree(s);
				ff_AddString(acc);
				MemFree(acc);
				AddLink("</a>");
				for (p = pp; iscospa(*p); p++) ;
			} while (*p != '\0');
			ff_AddString(".");
		} else {
			ff_AddString(string);
		}
		ff_EndPrint();
		return;
	}
	if (( p = StringStr(string, "AceView:")) != NULL) {
		/*	p += StringLen("AceView:");*/
			gi = atoi(p + StringLen("AceView:"));
			l = StringLen(link_ace) + StringLen(p);
			prefix = "<a href=%sl=%ld>"; 
			ll = StringLen(prefix); 
			s = (CharPtr)MemNew(l + ll);
			ss = (CharPtr)MemNew(p-string+1);
			StringNCpy(ss, string, p-string);
			ff_AddString(ss);
			MemFree(ss);
			sprintf(s, prefix, link_ace, gi);
			AddLink(s);
			MemFree(s);
			ff_AddString("AceView");
		/*	ff_AddInteger("%ld", (long) gi);*/
			AddLink("</a>");
		ff_EndPrint();
		return;
	}
	while ((p = StringStr(string, "gi:")) != NULL) {
		p += 3;
		s = TextSave(string, p-string);
		ff_AddString(s);
		gi = atoi(p);
		link = link_seq;
		prefix = "<a href=%sval=%d>";
		s = (CharPtr)MemNew(StringLen(link_seq)+ StringLen(prefix) + 12);
		sprintf(s, prefix, link, gi);
		AddLink(s);
		MemFree(s);
		ff_AddInteger("%d", gi);
		AddLink("</a>");
		while (IS_DIGIT(*p)) 
			p++;
		string = p;
	}
	www_str = www_featloc(string);
	ff_AddStringWithTildes(www_str);
	MemFree(www_str);
	ff_EndPrint();
	return;
}	/* PrintComment */


