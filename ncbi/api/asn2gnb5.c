/*   asn2gnb5.c
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
* File Name:  asn2gnb5.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* Version Creation Date:   10/21/98
*
* $Revision: 1.149 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objpubme.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbi.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

/* URLs */


static CharPtr link_muid = "http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=Retrieve&db=pubmed&list_uids=";

static CharPtr link_uspto = "http://patft.uspto.gov/netacgi/nph-Parser?patentnumber=";

static CharPtr link_cambia = "http://www.patentlens.net/patentlens/simple.cgi?patnum=";



/* www utility functions */

NLM_EXTERN Boolean GetWWW (IntAsn2gbJobPtr ajp) {
    return ajp->www;
}

NLM_EXTERN void FiniWWW (IntAsn2gbJobPtr ajp) {
    ajp->www = FALSE;
}

NLM_EXTERN void InitWWW (IntAsn2gbJobPtr ajp)
{
  ajp->www = TRUE;
}

NLM_EXTERN void FF_www_featloc(StringItemPtr ffstring, CharPtr loc)
{
  CharPtr ptr;

  if (loc == NULL) return;

  for ( ptr = loc; *ptr != '\0'; ++ptr ) {
    switch (*ptr) {
    case '<' :
      /*FFAddOneString (ffstring, "<", FALSE, FALSE, TILDE_IGNORE);*/
      FFAddOneString (ffstring, "&lt;", FALSE, FALSE, TILDE_IGNORE);
      break;
    case '>' :
      /*FFAddOneString (ffstring, ">", FALSE, FALSE, TILDE_IGNORE);*/
      FFAddOneString (ffstring, "&gt;", FALSE, FALSE, TILDE_IGNORE);
      break;
    default:
      FFAddOneChar(ffstring, *ptr, FALSE);
      break;
    }
  }
}


/* ************** */

typedef struct dbxrefurldata {
  CharPtr  tag;
  CharPtr  url;
} UrlData, PNTR UrlDataPtr;

static UrlData Nlm_url_base [] = {
  {"AceView/WormGenes",     "http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=worm&c=gene&q="},
  {"AFTOL",                 "http://aftol.biology.duke.edu/pub/displayTaxonInfo?aftol_id="},
  {"ApiDB",                 "http://www.apidb.org/apidb/showRecord.do?name=GeneRecordClasses.ApiDBGeneRecordClass&primary_key="},
  {"ApiDB_CryptoDB",        "http://cryptodb.org/cryptodb/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="},
  {"ApiDB_PlasmoDB",        "http://www.plasmodb.org/plasmo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="},
  {"ApiDB_ToxoDB",          "http://www.toxodb.org/toxo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&project_id=&primary_key="},
  {"ASAP",                  "https://asap.ahabs.wisc.edu/annotation/php/feature_info.php?FeatureID="},
  {"ATCC",                  "http://www.atcc.org/SearchCatalogs/linkin?id="},
  {"Axeldb",                "http://www.dkfz-heidelberg.de/tbi/services/axeldb/clone/xenopus?name="},
  {"BEETLEBASE",            "http://www.beetlebase.org/cgi-bin/report.cgi?name="},
  {"BioHealthBase",         "http://www.biohealthbase.org/GSearch/fluSegmentDetails.do?bhbSubmissionId="},
  {"BOLD",                  "http://www.boldsystems.org/connectivity/specimenlookup.php?processid="},
  {"CCDS",                  "http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA="},
  {"CDD",                   "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="},
  {"CK",                    "http://flybane.berkeley.edu/cgi-bin/cDNA/CK_clone.pl?db=CK&dbid="},
  {"COG",                   "http://www.ncbi.nlm.nih.gov/COG/new/release/cow.cgi?cog="},
  {"dbClone",               "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=clone&cmd=Retrieve&list_uids="},
  {"dbCloneLib",            "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=clonelib&cmd=Retrieve&list_uids="},
  {"dbEST",                 "http://www.ncbi.nlm.nih.gov/nucest/"},
  {"dbProbe",               "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=probe&cmd=Retrieve&list_uids="},
  {"dbSNP",                 "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs="},
  {"dbSTS",                 "http://www.ncbi.nlm.nih.gov/nuccore/"},
  {"dictyBase",             "http://dictybase.org/db/cgi-bin/gene_page.pl?dictybaseid="},
  {"ECOCYC",                "http://biocyc.org/ECOLI/new-image?type=GENE&object="},
  {"EcoGene",               "http://ecogene.org/geneInfo.php?eg_id="},
  {"ERIC",                  "http://www.ericbrc.org/genbank/dbxref/"},
  {"FANTOM_DB",             "http://fantom.gsc.riken.jp/db/annotate/main.cgi?masterid="},
  {"FLYBASE",               "http://flybase.bio.indiana.edu/.bin/fbidq.html?"},
  {"GABI",                  "http://www.gabipd.org/database/cgi-bin/GreenCards.pl.cgi?Mode=ShowSequence&App=ncbi&SequenceId="},
  {"GeneDB",                "http://www.genedb.org/genedb/Dispatcher?formType=navBar&submit=Search+for&organism=All%3Apombe%3Acerevisiae%3Adicty%3Aasp%3Atryp%3Aleish%3Amalaria%3Astyphi%3Aglossina&desc=yes&ohmr=%2F&name="},
  {"GeneID",                "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids="},
  {"GO",                    "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&depth=1&query=GO:"},
  {"GRIN",                  "http://www.ars-grin.gov/cgi-bin/npgs/acc/display.pl?"},
  {"H-InvDB",               "http://www.h-invitational.jp"},
  {"HGNC",                  "http://www.genenames.org/data/hgnc_data.php?hgnc_id="},
  {"HPRD",                  "http://www.hprd.org/protein/"},
  {"HSSP",                  "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-newId+-e+hssp-ID:"},
  {"IMGT/GENE-DB",          "http://imgt.cines.fr/cgi-bin/GENElect.jv?species=Homo+sapiens&query=2+"},
  {"IMGT/LIGM",             "http://imgt.cines.fr:8104/cgi-bin/IMGTlect.jv?query=202+"},
  {"InterimID",             "http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l="},
  {"InterPro",              "http://www.ebi.ac.uk/interpro/ISearch?mode=ipr&query="},
  {"ISD",                   "http://www.flu.lanl.gov/search/view_record.html?accession="},
  {"ISFinder",              "http://www-is.biotoul.fr/scripts/is/is_spec.idc?name="},
  {"JCM",                   "http://www.jcm.riken.go.jp/cgi-bin/jcm/jcm_number?JCM="},
  {"JGIDB",                 "http://genome.jgi-psf.org/cgi-bin/jgrs?id="},
  {"LocusID",               "http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l="},
  {"MaizeGDB",              "http://www.maizegdb.org/cgi-bin/displaylocusrecord.cgi?id="},
  {"MGI",                   "http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:"},
  {"MIM",                   "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id="},
  {"miRBase",               "http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?acc="},
  {"NBRC",                  "http://www.nbrc.nite.go.jp/NBRC2/NBRCCatalogueDetailServlet?ID=NBRC&CAT="},
  {"NextDB",                "http://nematode.lab.nig.ac.jp/cgi-bin/db/ShowGeneInfo.sh?celk="},
  {"niaEST",                "http://lgsun.grc.nia.nih.gov/cgi-bin/pro3?sname1="},
  {"NMPDR",                 "http://www.nmpdr.org/linkin.cgi?id="},
  {"NRESTdb",               "http://genome.ukm.my/nrestdb/db/single_view_est.php?id="},
  {"Osa1",                  "http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice/?name="},
  {"Pathema",               "http://pathema.jcvi.org/cgi-bin/Burkholderia/shared/GenePage.cgi?all=1&locus="},
  {"PBmice",                "http://www.idmshanghai.cn/PBmice/DetailedSearch.do?type=insert&id="},
  {"PBR",                   "http://www.poxvirus.org/query.asp?web_id="},
  {"PDB",                   "http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId="},
  {"PFAM",                  "http://pfam.sanger.ac.uk/family?acc="},
  {"PGN",                   "http://pgn.cornell.edu/cgi-bin/search/seq_search_result.pl?identifier="},
  {"PseudoCap",             "http://www.pseudomonas.com/getAnnotation.do?locusID="},
  {"RAP-DB",                "http://rapdb.dna.affrc.go.jp/cgi-bin/gbrowse_details/latest?name="},
  {"RATMAP",                "http://ratmap.gen.gu.se/ShowSingleLocus.htm?accno="},
  {"REBASE",                "http://rebase.neb.com/rebase/enz/"},
  {"RFAM",                  "http://www.sanger.ac.uk/cgi-bin/Rfam/getacc?"},
  {"RGD",                   "http://rgd.mcw.edu/query/query.cgi?id="},
  {"RiceGenes",             "http://ars-genome.cornell.edu/cgi-bin/WebAce/webace?db=ricegenes&class=Marker&object="},
  {"SEED",                  "http://www.theseed.org/linkin.cgi?id="},
  {"SGD",                   "http://db.yeastgenome.org/cgi-bin/SGD/locus.pl?locus="},
  {"SubtiList",             "http://genolist.pasteur.fr/SubtiList/genome.cgi?external_query+"},
  {"TAIR",                  "http://www.arabidopsis.org/servlets/TairObject?type=locus&name="},
  {"taxon",                 "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?"},
  {"TIGRFAM",               "http://cmr.tigr.org/tigr-scripts/CMR/HmmReport.cgi?hmm_acc="},
  {"UniGene",               "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=unigene&Cmd=Search&Term="},
  {"UniProtKB/Swiss-Prot",  "http://www.uniprot.org/uniprot/"},
  {"UniProtKB/TrEMBL",      "http://www.uniprot.org/uniprot/"},
  {"UniSTS",                "http://www.ncbi.nlm.nih.gov/genome/sts/sts.cgi?uid="},
  {"UNITE",                 "http://unite.ut.ee/bl_forw.php?nimi="},
  {"VBASE2",                "http://www.dnaplot.de/vbase2/vgene.php?id="},
  {"VBRC",                  "http://vbrc.org/query.asp?web_view=curation&web_id="},
  {"VectorBase",            "http://www.vectorbase.org/Genome/BRCGene/?feature="},
  {"WorfDB",                "http://worfdb.dfci.harvard.edu/search.pl?form=1&search="},
  {"WormBase",              "http://www.wormbase.org/db/gene/gene?class=CDS;name="},
  {"Xenbase",               "http://www.xenbase.org/gene/showgene.do?method=display&geneId="},
  {"ZFIN",                  "http://zfin.org/cgi-bin/webdriver?MIval=aa-markerview.apg&OID="},
};

static Int2 DbNameIsValid (
  CharPtr db
)

{
  Int2  L, R, mid;

  if (StringHasNoText (db)) return -1;

  L = 0;
  R = sizeof (Nlm_url_base) / sizeof (Nlm_url_base [0]);

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (Nlm_url_base [mid].tag, db) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  /* case sensitive comparison at end enforces strictness */

  if (StringCmp (Nlm_url_base [R].tag, db) == 0) {
    return R;
  }

  return -1;
}

static void FF_www_get_url (
  StringItemPtr ffstring,
  CharPtr db,
  CharPtr identifier,
  BioseqPtr bsp
)

{
  CharPtr       base = NULL, prefix = NULL, profix = NULL, ident = NULL, suffix = NULL, url = NULL, ptr, str;
  Char          ch, buf [128], id [20], taxname [128];
  Boolean       is_accession;
  /*
  Boolean       is_numeric;
  */
  Int2          R;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  if (ffstring == NULL || StringHasNoText (db) || StringHasNoText (identifier)) return;

  while (*identifier == ' ') {
    identifier++;
  }
  ident = identifier;

  R = DbNameIsValid (db);
  if (R < 0) {
    FFAddOneString (ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
    return;
  }

  url = Nlm_url_base [R].url;

  /* NCBI URL can be overridden by configuration file */

  if (StringNICmp (url, "http://www.ncbi.nlm.nih.gov/", 28) == 0) {
    if (GetAppParam ("NCBI", "WWWENTREZ", "NCBI_URL_BASE", NULL, buf, sizeof (buf))) {
      if (StringDoesHaveText (buf)) {
        base = buf;
        url += 28;
      }
    }
  }

  /* special cases */


  if (StringCmp (db, "BioHealthBase") == 0) {

    is_accession = FALSE;

    if (bsp != NULL && ValidateAccn (identifier) == 0) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_TPG :
          case SEQID_TPE :
          case SEQID_TPD :
          case SEQID_OTHER :
            tsip = (TextSeqIdPtr) sip->data.ptrvalue;
            if (tsip != NULL) {
              if (StringICmp (tsip->accession, identifier) == 0) {
                is_accession = TRUE;
              }
            }
            break;
          default :
            break;
        }
      }
    }

    if (is_accession && bsp != NULL &&
        BioseqToGeneticCode (bsp, NULL, NULL, NULL, taxname, sizeof (taxname), NULL)) {
      ptr = StringChr (taxname, ' ');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      ch = taxname [0];
      if (IS_UPPER (ch)) {
        taxname [0] = TO_LOWER (ch);
      }
      if (StringNICmp (taxname, "influenza", 9) == 0) {
        url = "http://www.biohealthbase.org/GSearch/fluSegmentDetails.do?decorator=";
        prefix = taxname;
        profix = "&ncbiGenomicAccession=";
      } else {
        url = "http://www.biohealthbase.org/GSearch/details.do?decorator=";
        prefix = taxname;
        profix = "&locus=";
      }
    }

  } else if (StringCmp (db, "dbSTS") == 0) {

    /*
    is_numeric = TRUE;
    str = identifier;
    ch = *str;
    while (ch != '\0') {
      if (! IS_DIGIT (ch)) {
        is_numeric = FALSE;
      }
      str++;
      ch = *str;
    }

    if (is_numeric) {
      prefix = "val=gnl|dbsts|";
    } else if (ValidateAccn (identifier) == 0) {
      prefix = "val=";
    } else {
      FFAddOneString (ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
      return;
    }
    */

  } else if (StringCmp (db, "FLYBASE") == 0) {

    if (StringStr (identifier, "FBa") != NULL ) {
      url = "http://www.fruitfly.org/cgi-bin/annot/fban?";
    }

  } else if (StringCmp (db, "dictyBase") == 0) {

    /*
    if (StringChr (identifier, '_') == NULL) {
      url = "http://dictybase.org/db/cgi-bin/feature_page.pl?primary_id=";
    }
    */

  } else if (StringCmp (db, "GDB") == 0) {

    str = StringStr (identifier, "G00-");
    if (str != NULL) {
      ptr = id;
      str += 4;
      ch = *str;
      while (ch != '\0') {
        if (ch != '-') {
          *ptr = ch;
          ptr++;
        }
        str++;
        ch = *str;
      }
      *ptr = '\0';
      ident = id;
    } else {
      ch = *identifier;
      if (! IS_DIGIT (ch)) {
        FFAddOneString (ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
        return;
      }
    }

  } else if (StringCmp (db, "H_InvDB") == 0) {

    if (StringStr (identifier, "HIT") != NULL) {
      url = "http://www.jbirc.aist.go.jp/hinv/hinvsys/servlet/ExecServlet?KEN_INDEX=0&KEN_TYPE=30&KEN_STR=";
    } else if (StringStr (identifier, "HIX") != NULL) {
      url = "http://www.jbirc.aist.go.jp/hinv/hinvsys/servlet/ExecServlet?KEN_INDEX=0&KEN_TYPE=31&KEN_STR=";
    }

  } else if (StringCmp (db, "IMGT/GENE-DB") == 0) {

    if (bsp != NULL && BioseqToGeneticCode (bsp, NULL, NULL, NULL, taxname, sizeof (taxname), NULL)) {
      if (StringCmp (taxname, "Homo sapiens") == 0) {
        url = "http://imgt.cines.fr/cgi-bin/GENElect.jv?species=Homo+sapiens&query=2+";
      }
      if (StringCmp (taxname, "Mus musculus") == 0) {
        url = "http://imgt.cines.fr/cgi-bin/GENElect.jv?species=Mus+musculus&query=2+";
      }
    }

  } else if (StringCmp (db, "niaEST") == 0) {

    suffix = "&val=1";

  } else if (StringCmp (db, "RAP-DB") == 0) {

    suffix = ";class=locus_id";

  } else if (StringCmp (db, "REBASE") == 0) {

    suffix = ".html";

  } else if (StringCmp (db, "taxon") == 0) {

    ch = *identifier;
    if (IS_DIGIT (ch)) {
      prefix = "id=";
    } else {
      prefix = "name=";
    }

  }

  /* now generate URL */

  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
  if (StringDoesHaveText (base)) {
    FFAddOneString (ffstring, base, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneString (ffstring, url, FALSE, FALSE, TILDE_IGNORE);
  if (StringDoesHaveText (prefix)) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
  }
  if (StringDoesHaveText (profix)) {
    FFAddOneString (ffstring, profix, FALSE, FALSE, TILDE_IGNORE);
  }
  if (StringDoesHaveText (ident)) {
    FFAddOneString (ffstring, ident, FALSE, FALSE, TILDE_IGNORE);
  }
  if (StringDoesHaveText (suffix)) {
    FFAddOneString (ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddTextToString (ffstring, "\">", identifier, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

NLM_EXTERN void FF_www_db_xref (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr db,
  CharPtr identifier,
  BioseqPtr bsp
)
{
  if (ffstring == NULL || StringHasNoText (db) || StringHasNoText (identifier)) return;

  if (GetWWW (ajp)) {
    FFAddTextToString (ffstring, db, ":", NULL, FALSE, FALSE, TILDE_IGNORE);
    FF_www_get_url (ffstring, db, identifier, bsp);
  } else {
    FFAddTextToString (ffstring, db, ":", identifier, FALSE, FALSE, TILDE_IGNORE);
  }
}

NLM_EXTERN void FF_Add_NCBI_Base_URL (
  StringItemPtr ffstring,
  CharPtr url
)

{
  CharPtr  base = NULL;
  Char     buf [128];

  if (ffstring == NULL || StringHasNoText (url)) return;

  /* NCBI URL can be overridden by configuration file */

  if (StringNICmp (url, "http://www.ncbi.nlm.nih.gov/", 28) == 0) {
    if (GetAppParam ("NCBI", "WWWENTREZ", "NCBI_URL_BASE", NULL, buf, sizeof (buf))) {
      if (StringDoesHaveText (buf)) {
        base = buf;
        url += 28;
      }
    }
  }

  if (StringDoesHaveText (base)) {
    FFAddOneString (ffstring, base, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneString (ffstring, url, FALSE, FALSE, TILDE_IGNORE);
}


/* ************** */


/* public function to get URLs for collaboration-approved db_xrefs */

static Boolean links_loaded = FALSE;

NLM_EXTERN CharPtr asn2gnbk_dbxref (
  DbtagPtr dbt
)

{
  IntAsn2gbJobPtr  ajp;
  Char             buf [80];
  StringItemPtr    ffstring;
  ObjectIdPtr      oip;
  CharPtr          ptr;
  CharPtr          str;
  CharPtr          tmp;

  if (dbt == NULL) return NULL;
  if (StringHasNoText (dbt->db)) return NULL;
  oip = dbt->tag;
  if (oip == NULL) return NULL;

  if (! StringHasNoText (oip->str)) {
    if (StringLen (dbt->db) + StringLen (oip->str) < 80) {
      sprintf (buf, "%s", oip->str);
    }
  } else {
    sprintf (buf, "%ld", (long) oip->id);
  }

  ajp = (IntAsn2gbJobPtr) MemNew (sizeof (IntAsn2gbJob));
  if (ajp == NULL) return NULL;
  ffstring = FFGetString (ajp);
  if ( ffstring == NULL ) return NULL;

  if (! links_loaded) {
    InitWWW (ajp);
    links_loaded = TRUE;
  }
  ajp->www = TRUE;

  FF_www_db_xref (ajp, ffstring, dbt->db, buf, NULL);

  ajp->www = FALSE;

  str = FFToCharPtr (ffstring);

  FFRecycleString (ajp, ffstring);
  /*
  MemFree (ajp);
  */
  asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  tmp = StringChr (str, '<');
  if (tmp != NULL) {
    ptr = StringSave (tmp);
    tmp = StringChr (ptr, '>');
    if (tmp != NULL) {
      tmp++;
      *tmp = '\0';
    }
    MemFree (str);
    str = ptr;
  } else {
    str = MemFree (str);
  }

  return str;
}

/* format references section */

NLM_EXTERN AuthListPtr GetAuthListPtr (
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  AuthListPtr  alp = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  ValNodePtr   vnp;

  if (csp != NULL) {
    alp = csp->authors;
    if (alp != NULL) return alp;
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          alp = cgp->authors;
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          alp = csp->authors;
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          alp = cap->authors;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          alp = cbp->authors;
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          alp = cpp->authors;
        }
        break;
      default :
        break;
    }

    if (alp != NULL) return alp;
  }

  return NULL;
}

static CharPtr MakeSingleAuthorString (
  FmtType format,
  CharPtr prefix,
  CharPtr name,
  CharPtr initials,
  CharPtr suffix,
  IndxPtr index,
  GBReferencePtr gbref
)

{
  Char     ch;
  Char     dummy [10];
  size_t   len;
  CharPtr  nametoindex;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (name == NULL) return NULL;

  /* !!! clean up 'et al' as (presumably) last author !!! */

  /* !!! temporary to suppress diff !!! */
  {
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    if (StringCmp (prefix, " and ") == 0) {
      prefix = NULL;
      dummy [0] = ' ';
      StringNCpy_0 (dummy + 1, name, sizeof (dummy) - 1);
      name = dummy;
    }
  }
  }
  /*
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    name = "et al.";
    if (StringCmp (prefix, " and ") == 0) {
      prefix = ", ";
    }
  }
  */

  len = StringLen (name) + StringLen (initials) + StringLen (suffix) + StringLen (prefix);
  str = MemNew (sizeof (Char) * (len + 4));
  if (str == NULL) return NULL;

  ptr = str;
  if (! StringHasNoText (prefix)) {
    ptr = StringMove (ptr, prefix);
  }
  nametoindex = ptr;

  /* initials and suffix to support structured name fields */

  tmp = StringMove (ptr, name);
  if (! StringHasNoText (initials)) {
    tmp = StringMove (tmp, ",");
    tmp = StringMove (tmp, initials);
  }
  if (! StringHasNoText (suffix)) {
    tmp = StringMove (tmp, " ");
    tmp = StringMove (tmp, suffix);
  }

  /* optionally populate indexes for NCBI internal database */

  if (index != NULL) {
    ValNodeCopyStrToHead (&(index->authors), 0, nametoindex);
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbref != NULL) {
    ValNodeCopyStr (&(gbref->authors), 0, nametoindex);
  }

  /* if embl, remove commas in individual names, starting after prefix */

  if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    tmp = ptr;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == ',') {
        *tmp = ' ';
      }
      tmp++;
      ch = *tmp;
    }
  }

  return str;
}

NLM_EXTERN CharPtr GetAuthorsString (
  FmtType format,
  AuthListPtr alp,
  CharPtr PNTR consortP,
  IndxPtr index,
  GBReferencePtr gbref
)

{
  AuthorPtr    ap;
  ValNodePtr   clist;
  ValNodePtr   conslist;
  Int2         count;
  ValNodePtr   head = NULL;
  ValNodePtr   names;
  ValNodePtr   next;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  ValNodePtr   pidlist;
  CharPtr      prefix = NULL;
  CharPtr      str;
  ValNodePtr   vnp;

  if (alp == NULL) return NULL;

  /*
  alp = AsnIoMemCopy ((Pointer) alp,
                      (AsnReadFunc) AuthListAsnRead,
                      (AsnWriteFunc) AuthListAsnWrite);
  if (alp == NULL) return NULL;
  */

  count = 0;
  if (alp->choice == 1) {

    pidlist = NULL;
    conslist = NULL;

    for (names = alp->names; names != NULL; names = names->next) {
      ap = (AuthorPtr) names->data.ptrvalue;
      if (ap == NULL) continue;
      pid = ap->name;
      if (pid == NULL) continue;
      if (pid->choice == 2 || pid->choice == 3 || pid->choice == 4) {
        ValNodeAddPointer (&pidlist, 0, (Pointer) pid);
      } else if (pid->choice == 5) {
        ValNodeAddPointer (&conslist, 0, (Pointer) pid);
      }
    }

    for (vnp = pidlist; vnp != NULL; vnp = vnp->next) {
      next = vnp->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = NULL;
      pid = (PersonIdPtr) vnp->data.ptrvalue;
      if (pid->choice == 2) {
        nsp = (NameStdPtr) pid->data;
        if (nsp != NULL) {
          if (! StringHasNoText (nsp->names [0])) {
            str = MakeSingleAuthorString (format, prefix, nsp->names [0], nsp->names [4], nsp->names [5], index, gbref);
          } else if (! StringHasNoText (nsp->names [3])) {
            str = MakeSingleAuthorString (format, prefix, nsp->names [3], NULL, NULL, index, gbref);
          }
        }
      } else if (pid->choice == 3 || pid->choice == 4) {
        str = MakeSingleAuthorString (format, prefix, (CharPtr) pid->data, NULL, NULL, index, gbref);
      }
      if (str != NULL) {
        ValNodeAddStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }

    prefix = NULL;
    clist = NULL;
    for (vnp = conslist; vnp != NULL; vnp = vnp->next) {
      str = NULL;
      pid = (PersonIdPtr) vnp->data.ptrvalue;
      if (pid->choice == 5) {
        str = MakeSingleAuthorString (format, prefix, (CharPtr) pid->data, NULL, NULL, index, NULL);
        if (str != NULL) {
          ValNodeAddStr (&clist, 0, str);
        }
        prefix = "; ";
      }
    }
    if (clist != NULL) {
      str = MergeFFValNodeStrs (clist);
      if ((! StringHasNoText (str)) && consortP != NULL && *consortP == NULL) {
        *consortP = StringSave (str);
      }

      /* optionally populate gbseq for XML-ized GenBank format */

      if (gbref != NULL) {
        gbref->consortium = StringSave (str);
      }

      str = MemFree (str);
      ValNodeFreeData (clist);
    }

    ValNodeFree (pidlist);
    ValNodeFree (conslist);

  } else if (alp->choice == 2 || alp->choice == 3) {
    for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
      next = vnp->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = MakeSingleAuthorString (format, prefix, (CharPtr) vnp->data.ptrvalue, NULL, NULL, index, gbref);
      if (str != NULL) {
        ValNodeAddStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }
  }

  str = MergeFFValNodeStrs (head);

  ValNodeFreeData (head);

  /*
  AuthListFree (alp);
  */

  return str;
}

/*
Strips all spaces in string in following manner. If the function
meet several spaces (spaces and tabs) in succession it replaces them
with one space. Strips all spaces after '(' and before ')'
*/

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

static Boolean AllCaps (
  CharPtr p
)

{
  if (p == NULL) return FALSE;

  for (p++; p != NULL && *p != '\0'; p++) {
    if (IS_LOWER (*p)) return FALSE;
  }
  return TRUE;
}

static void CleanEquals (
  CharPtr p
)

{
  if (p == NULL) return;

  for (; *p != '\0'; p++) {
    if (*p == '\"') {
      *p = '\'';
    }
  }
}

static CharPtr GetPubTitle (
  FmtType format,
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  Char             ch;
  CitPatPtr        cpp;
  MedlineEntryPtr  mep;
  CharPtr          ptr;
  CharPtr          title = NULL;
  ValNodePtr       ttl = NULL;
  ValNodePtr       vnp;

  if (csp != NULL) {
    if (format == GENBANK_FMT || format == GENPEPT_FMT) {
      title = "Direct Submission";
      return StringSave (title);
    } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
      return NULL;
    }
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (! StringHasNoText (cgp->title)) return StringSave (cgp->title);
          if (! StringHasNoText (cgp->cit)) {
            ptr = StringStr (cgp->cit, "Title=\"");
            if (ptr != NULL) {
              title = StringSave (ptr + 7);
              for (ptr = title; *ptr != '\0'; ptr++) {
                if (*ptr == '"') {
                  *ptr = '\0';
                  break;
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (format == GENBANK_FMT || format == GENPEPT_FMT) {
            title = "Direct Submission";
            return StringSave (title);
          } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
            return NULL;
          }
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            ttl = cap->title;
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          ttl = cap->title;
        }
        break;
      /* case PUB_Book : */
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          ttl = cbp->title;
          if (ttl != NULL) {
            title = (CharPtr) ttl->data.ptrvalue;
            if (! StringHasNoText (title)) {
              title = StringSave (title);
              if (StringLen (title) > 3) {
                ch = *title;
                if (IS_LOWER (ch)) {
                  *title = TO_UPPER (ch);
                }
                ptr = title;
                if (AllCaps (ptr)) {
                  for (ptr++; ptr != NULL && *ptr != '\0'; ptr++) {
                    ch = *ptr;
                    *ptr = TO_LOWER (ch);
                  }
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          title = cpp->title;
          if (! StringHasNoText (title)) {
            return StringSave (title);
          }
        }
        break;
      default :
        break;
    }

    if (ttl != NULL) {
      title = (CharPtr) ttl->data.ptrvalue;
      if (! StringHasNoText (title)) {
        return StringSave (title);
      }
    }
  }

  return NULL;
}

static void CleanPubTitle (
  CharPtr title
)

{
  CharPtr  p;
  Boolean  remove_it;

  if (title == NULL) return;

  CleanEquals (title);

  for (p = title + StringLen (title) - 1; p > title + 2; p--) {
    if (*p == ' ') {
      *p = '\0';
    } else if (*p == '.') {
      remove_it = FALSE;
      if (p > title + 5) {
        if (*(p - 1) != '.' || *(p - 2) != '.') {
          remove_it = TRUE;
        }
      }
      if (remove_it) {
        *p = '\0';
      }
      break;
    } else {
      break;
    }
  }
}

/*
medline type page numbering is expanded (e.g., 125-35 -> 125-135,
F124-34 -> F124-F134, 12a-c -> 12a-12c).
If only one page is given, this is output without a dash.
Expanded numbering is validated to ensure that the
first number is smaller than or equal to the second and
that the first letter is less than or identical to the second
(i.e., a < c).  If the input is all letters (i.e., roman numerals)
this is not validated.

Return values:
 0 : valid page numbering.
-1 : invalid page numbering.
*/

#define MAX_PAGE_DIGITS 12

static Int2 FixPages (
  CharPtr out_pages,
  CharPtr in_pages
)

{
  Boolean dash=TRUE, first_alpha;
  Char firstbegin[MAX_PAGE_DIGITS];
  Char secondbegin[MAX_PAGE_DIGITS];
  Char firstend[MAX_PAGE_DIGITS];
  Char secondend[MAX_PAGE_DIGITS];
  Char temp[MAX_PAGE_DIGITS];
  CharPtr alphabegin, numbegin, alphaend, numend, ptr, in=in_pages;
  Int2 diff, index, retval=0;
  Int2 length_nb, length_ab, length_ne, length_ae;
  Int4 num1=0, num2=0;

  if (in_pages == NULL) return retval;

  while (*in != '\0')
  {      /* Check for digits in input*/
    if (IS_DIGIT(*in))
      break;
    in++;
  }

  if (*in == '\0' || (in != in_pages && *(in-1) == ' '))
  {    /* if all letters (i.e. roman numerals), put out. */
    out_pages = StringCpy(out_pages, in_pages);
    return retval;
  }

  in = in_pages;
  if (IS_DIGIT(*in))
  {      /* Do digits come first? */
    first_alpha = FALSE;
    index=0;
    while (IS_DIGIT(*in) || *in == ' ')
    {
      firstbegin[index] = *in;
      if (*in != ' ')
        index++;
      in++;
      if (*in == '-')
        break;

    }
    firstbegin[index] = '\0';
    index=0;
    if (*in != '-')
    {    /* After digits look for letters. */
      while (IS_ALPHA(*in)  || *in == ' ')
      {
        secondbegin[index] = *in;
        index++;
        in++;
        if (*in == '-')
          break;
      }
    }
    secondbegin[index] = '\0';
    if (*in == '-')    /* if dash is not present, note */
      in++;
    else
      dash=FALSE;
    index=0;
    while (IS_DIGIT(*in) || *in == ' ')
    {      /* Look for digits.  */
      firstend[index] = *in;
      if (*in != ' ')
        index++;
      in++;
    }
    firstend[index] = '\0';
    index=0;
    if (*in != '\0')
    {      /* Look for letters again. */
      while (IS_ALPHA(*in)  || *in == ' ')
      {
        secondend[index] = *in;
        index++;
        in++;
      }
    }
    secondend[index] = '\0';
  }
  else
  {      /* Do letters come first? */
    first_alpha = TRUE;
    index=0;
    while (IS_ALPHA(*in) || *in == ' ')
    {
      firstbegin[index] = *in;
      index++;
      in++;
      if (*in == '-')
        break;
    }
    firstbegin[index] = '\0';
    index=0;
    if (*in != '-')
    {    /* After letters look for digits.   */
      while (IS_DIGIT(*in)  || *in == ' ')
      {
        secondbegin[index] = *in;
        if (*in != ' ')
          index++;
        in++;
        if (*in == '-')
          break;
      }
    }
    secondbegin[index] = '\0';
    if (*in == '-')    /* Note if dash is missing. */
      in++;
    else
      dash=FALSE;
    index=0;
    while (IS_ALPHA(*in) || *in == ' ')
    {    /* Look for letters again. */
      firstend[index] = *in;
      index++;
      in++;
    }
    firstend[index] = '\0';
    index=0;
    if (*in != '\0')
    {    /* Any digits here? */
      while (IS_DIGIT(*in)  || *in == ' ')
      {
        secondend[index] = *in;
        if (*in != ' ')
          index++;
        in++;
      }
    }
    secondend[index] = '\0';
  }

  if (first_alpha)
  {
    alphabegin = firstbegin;
    numbegin = secondbegin;
    alphaend = firstend;
    numend = secondend;
  }
  else
  {
    numbegin = firstbegin;
    alphabegin = secondbegin;
    numend = firstend;
    alphaend = secondend;
  }

  length_nb = StringLen(numbegin);
  length_ab = StringLen(alphabegin);
  length_ne = StringLen(numend);
  length_ae = StringLen(alphaend);

  /* If no dash, but second letters or numbers present, reject. */
  if (dash == FALSE)
  {
    if (length_ne != 0 || length_ae != 0)
      retval = -1;
  }
  /* Check for situations like "AAA-123" or "222-ABC". */
  if (dash == TRUE)
  {
    if (length_ne == 0 && length_ab == 0)
      retval = -1;
    else if (length_ae == 0 && length_nb == 0)
      retval = -1;
  }

  /* The following expands "F502-512" into "F502-F512" and
  checks, for entries like "12a-12c" that a > c.  "12aa-12ab",
  "125G-137A", "125-G137" would be rejected. */
  if (retval == 0)
  {
    if (length_ab > 0)
    {
      if (length_ae > 0)
      {
        if (StringCmp(alphabegin, alphaend) != 0)
        {
          if (length_ab != 1 || length_ae != 1)
            retval = -1;
          else if (*alphabegin > *alphaend)
            retval = -1;
        }
      }
      else
      {
        alphaend = alphabegin;
        length_ae = length_ab;
      }
    }
    else if (length_ae > 0)
      retval = -1;
  }

/* The following expands "125-37" into "125-137".  */
  if (retval == 0)
  {
    if (length_nb > 0)
    {
      if (length_ne > 0)
      {
        diff = length_nb - length_ne;
        if (diff > 0)
        {
          index=0;
          while (numend[index] != '\0')
          {
            temp[index+diff] = numend[index];
            index++;
          }
          temp[index+diff] = numend[index];
          for (index=0; index<diff; index++)
            temp[index] = numbegin[index];
          index=0;
          while (temp[index] != '\0')
          {
            numend[index] = temp[index];
            index++;
          }
          numend[index] = temp[index];
        }
      }
      else
      {
        numend = numbegin;
        length_ne = length_nb;
      }

    }
    else if (length_ne > 0)
      retval = -1;
  /* Check that the first number is <= the second (expanded) number. */
    if (retval == 0)
    {
  /*    sscanf(numbegin, "%ld", &num_type);
      num1 = (Int4) num_type;
      sscanf(  numend, "%ld", &num_type);
      num2 = (Int4) num_type;
  */
      num1 = (Int4) atol(numbegin);
      num2 = (Int4) atol(numend);
      if (num2 < num1)
        retval = -1;
    }
  }

  if (retval == -1)
  {
    out_pages = StringCpy(out_pages, in_pages);
  }
  else
  {
    ptr = out_pages;
  /* Place expanded and validated page numbers into "out_pages". */
    if (first_alpha)
    {
      while (*alphabegin != '\0')
      {
        *ptr = *alphabegin;
        alphabegin++;
        ptr++;
      }
      while (*numbegin != '\0')
      {
        *ptr = *numbegin;
        numbegin++;
        ptr++;
      }
      if (dash == TRUE)
      {
        *ptr = '-';
        ptr++;
        while (*alphaend != '\0')
        {
          *ptr = *alphaend;
          alphaend++;
          ptr++;
        }
        while (*numend != '\0')
        {
          *ptr = *numend;
          numend++;
          ptr++;
        }
      }
      *ptr = '\0';
    }
    else
    {
      while (*numbegin != '\0')
      {
        *ptr = *numbegin;
        numbegin++;
        ptr++;
      }
      while (*alphabegin != '\0')
      {
        *ptr = *alphabegin;
        alphabegin++;
        ptr++;
      }
      if (dash == TRUE)
      {
        *ptr = '-';
        ptr++;
        while (*numend != '\0')
        {
          *ptr = *numend;
          numend++;
          ptr++;
        }
        while (*alphaend != '\0')
        {
          *ptr = *alphaend;
          alphaend++;
          ptr++;
        }
      }
      *ptr = '\0';
    }
  }
  return retval;
}

/* !!! still need to add StripParanthesis equivalent !!! */

static void DoSup (
  ValNodePtr PNTR head,
  CharPtr issue,
  CharPtr part_sup,
  CharPtr part_supi
)

{
  size_t   len;
  CharPtr  str;
  CharPtr  temp;

  len = StringLen (issue) + StringLen (part_sup) + StringLen (part_supi) + 25;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return;
  temp = str;

  if (! StringHasNoText (part_sup)) {
    *temp = ' ';
    temp++;
    temp = StringMove (temp, part_sup);
  }
  if (StringHasNoText (issue) && StringHasNoText (part_supi)) {
    ValNodeCopyStr (head, 0, str);
    MemFree (str);
    return;
  }
  *temp = ' ';
  temp++;
  *temp = '(';
  temp++;
  if (! StringHasNoText (issue)) {
    temp = StringMove (temp, issue);
  }
  if (! StringHasNoText (part_supi)) {
    *temp = ' ';
    temp++;
    temp = StringMove (temp, part_supi);
  }
  *temp = ')';
  temp++;
  ValNodeCopyStr (head, 0, str);
  MemFree (str);
}

static CharPtr FormatCitJour (
  FmtType format,
  Boolean citArtIsoJta,
  CitJourPtr cjp
)

{
  Char        buf [256];
  DatePtr     dp;
  Boolean     electronic_journal = FALSE;
  ValNodePtr  head = NULL;
  ImprintPtr  imp;
  CharPtr     issue = NULL;
  Char        pages [128];
  CharPtr     part_sup = NULL;
  CharPtr     part_supi = NULL;
  CharPtr     rsult = NULL;
  CharPtr     title;
  ValNodePtr  ttl;
  CharPtr     volume;
  Char        year [8];

  if (cjp == NULL) return NULL;

  ttl = cjp->title;
  if (ttl == NULL) return NULL;

  /* always use iso_jta title if present */

  while (ttl != NULL && ttl->choice != Cit_title_iso_jta) {
    ttl = ttl->next;
  }

  imp = cjp->imp;
  if (imp == NULL) return NULL;

  /* release mode requires iso_jta title */

  if (imp->pubstatus == 3 || imp->pubstatus == 10) {
    electronic_journal = TRUE;
  }

  if (ttl == NULL) {
    ttl = cjp->title;
    if (ttl != NULL && ttl->choice == Cit_title_name) {
      title = (CharPtr) ttl->data.ptrvalue;
      if (title != NULL && StringNCmp (title, "(er)", 4) == 0) {
        electronic_journal = TRUE;
      }
    }
    if (citArtIsoJta && (! electronic_journal)) return NULL;
  }

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished %s", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  /*
  if (imp->pubstatus == 3 || imp->pubstatus == 10) {
    ValNodeCopyStr (&head, 0, "(er) ");
  }
  */

  ValNodeCopyStr (&head, 0, title);

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  if (electronic_journal) {
    StringNCpy_0 (pages, imp->pages, sizeof (pages));
  } else {
    FixPages (pages, imp->pages);
  }

  if (! StringHasNoText (volume)) {
    AddValNodeString (&head, " ", volume, NULL);
  }

  if ((! StringHasNoText (volume)) || (! StringHasNoText (pages))) {
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ", ", pages, NULL);
    }
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ":", pages, NULL);
    } else if (imp->prepub == 2 || (StringHasNoText (volume))) {
      ValNodeCopyStr (&head, 0, " 0:0-0");
    }
  }

  ValNodeCopyStr (&head, 0, year);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    } else if (imp->pubstatus == 10 && StringHasNoText (pages)) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr MakeAffilStr (
  AffilPtr afp
)

{
  ValNodePtr  head = NULL;
  CharPtr     prefix = "";
  CharPtr     rsult = NULL;

  if (afp == NULL) return NULL;

  if (! StringHasNoText (afp->affil)) {
    ValNodeCopyStr (&head, 0, afp->affil);
    prefix = ", ";
  }

  if (afp->choice == 2) {
    if (! StringHasNoText (afp->div)) {
      AddValNodeString (&head, prefix, afp->div, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->street)) {
      AddValNodeString (&head, prefix, afp->street, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->city)) {
      AddValNodeString (&head, prefix, afp->city, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->sub)) {
      AddValNodeString (&head, prefix, afp->sub, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->country)) {
      AddValNodeString (&head, prefix, afp->country, NULL);
      prefix = ", ";
    }
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetAffil (
  AffilPtr afp
)

{
  Boolean need_comma=FALSE;
  CharPtr string=NULL, temp, ptr;
  Char ch;
  Int2 aflen=15;

  if (afp == NULL) return NULL;
  if (afp) {
    if (afp -> choice == 1){
      if (afp -> affil){
        aflen += StringLen(afp -> affil);
      }
    }else if (afp -> choice == 2){
      aflen += StringLen (afp -> affil) +
      StringLen (afp -> div) +
      StringLen (afp -> city) +
      StringLen (afp -> sub) +
      StringLen (afp -> street) +
      StringLen (afp -> country) + StringLen(afp->postal_code);
    }

    temp = string = MemNew(aflen);

    if ( afp -> choice == 1){
       if (afp -> affil){
        ptr = afp->affil;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
       }
    }else if (afp -> choice == 2){

      if( afp -> div) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->div;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if(afp -> affil) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->affil;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if(afp -> street) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->street;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> city) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->city;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> sub) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->sub;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> postal_code){
        *temp = ' ';
        temp++;
        ptr = afp->postal_code;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
      }

      if( afp -> country){
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->country;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }
    }
    temp++;
    *temp = '\0';
  }

    /* convert double quotes to single quotes */

    ptr = string;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }

  return string;
}

NLM_EXTERN CharPtr GetFlatFileAffilString (AffilPtr afp)
{
  return GetAffil (afp);
}


static CharPtr FormatCitBookArt (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  AuthListPtr  alp;
  CharPtr      book_title = NULL;
  Char         buf [256];
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      issue = NULL;
  ValNodePtr   names = NULL;
  Char         pages [128];
  CharPtr      part_sup = NULL;
  CharPtr      part_supi = NULL;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      title;
  ValNodePtr   ttl;
  ValNodePtr   vnp;
  CharPtr      volume;
  Char         year [8];

  if (cbp == NULL) return NULL;

  ttl = cbp->title;
  if (ttl == NULL) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "(%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, "(");
      StringNCat (year, dp->str, 4);
      StringNCat (year, ")", 1);
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished %s", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  ValNodeCopyStr (&head, 0, "(in) ");

  alp = cbp->authors;
  if (alp != NULL) {
    str = GetAuthorsString (format, alp, NULL, NULL, NULL);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
      names = alp->names;
      if (names != NULL) {
        if (names->next != NULL) {
          ValNodeCopyStr (&head, 0, " (Eds.);");
        } else {
          ValNodeCopyStr (&head, 0, " (Ed.);");
        }
      }
      ValNodeCopyStr (&head, 0, "\n");
    }
    MemFree (str);
  }

  book_title = StringSaveNoNull (title);
  vnp = ValNodeAddStr (&head, 0, book_title);
  if (book_title != NULL) {

    /* make book title all caps */

    title = book_title;
    ch = *title;
    while (ch != '\0') {
      *title = TO_UPPER (ch);
      title++;
      ch = *title;
    }
  }

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  FixPages (pages, imp->pages);

  if ((! StringHasNoText (volume)) && (StringCmp (volume, "0") != 0)) {
    AddValNodeString (&head, ", Vol. ", volume, NULL);
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (! StringHasNoText (pages)) {
    AddValNodeString (&head, ": ", pages, NULL);
  }

  if (book_title != NULL) {
    ValNodeCopyStr (&head, 0, ";\n");
  }

  afp = imp->pub;
  if (afp != NULL) {
    str = MakeAffilStr (afp);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
      ValNodeCopyStr (&head, 0, " ");
      MemFree (str);
    }
  }

  AddValNodeString (&head, NULL, year, NULL);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitBook (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr   afp;
  char       year[5];
  CharPtr    bookTitle=NULL;
  CharPtr    retval = NULL;
  CharPtr    temp;
  DatePtr    dp;
  ImprintPtr ip;
  int        aflen = 0;
  CharPtr    p;
  CharPtr    affilStr = NULL;

  /* Check parameters */

  if (cbp == NULL)
    return NULL;

  if ( cbp -> othertype != 0)
    return NULL;

  ip = cbp -> imp;

  /* Format the year */

  dp = ip -> date;
  year[0] = '\0';

  if ( dp -> data[0] == 1)
    sprintf(year,"%ld",(long) ( 1900+dp -> data[1]));
  else
    {
      StringNCpy( (CharPtr) year, (CharPtr) dp -> str, (size_t) 4);
      year[4] = '\0';
    }

  /* Get the book title */

  if (cbp->title)
    bookTitle = StringSave(cbp -> title -> data.ptrvalue);

  /* Get the affiliation length */

  if ( ip -> pub){
    afp = ip -> pub;
    aflen = StringLen(afp -> affil)+ 5;
    if ( afp -> choice == 2){
      aflen += 3 + StringLen(afp -> div);
      aflen += 3 + StringLen(afp -> street);
      aflen += 3 + StringLen(afp -> city);
      aflen += 3 + StringLen(afp -> sub);
      aflen += 3 + StringLen(afp -> country);
    }
  } else{
    aflen = 22;
  }
  if (ip->prepub == 2)
    aflen += 10;

  /* Create a Char String big enough to hold */
  /* the title, year, and affiliation.       */

  temp = retval = MemNew( (size_t) (30+StringLen( bookTitle)+StringLen( year) + aflen) );

  /* Convert the title to upper case and */
  /* add it to the string.               */

  for ( p = bookTitle; *p; p++)
    *p = TO_UPPER(*p);

  /* temp = StringMove(temp, "Book: "); */
  temp = StringMove(temp, "(in) ");
  temp = StringMove(temp, bookTitle);
  temp = StringMove(temp, ".");

  /* Add the affiliation to the string */

  if ( ip -> pub)
    {
      afp = ip -> pub;
      *temp = ' ';
      temp++;
      affilStr = MakeAffilStr(afp);
      temp = StringMove(temp,affilStr);
    }

  /* Add the year to the string */

  if (year[0] != '\0')
    {
      if (affilStr != NULL)
        temp = StringMove(temp," (");
      else
        temp = StringMove(temp, "(");
      temp = StringMove(temp, year);
      temp = StringMove(temp, ")");
    }

  /* If in press, add note */

  if (ip->prepub == 2)
    temp = StringMove(temp, ", In press");

  /* Clean up and return */

  if (bookTitle)
    MemFree(bookTitle);

  return retval;

}

static CharPtr FormatThesis (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      ptr;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      suffix = NULL;
  Char         year [8];

  if (cbp == NULL) return NULL;
  if (cbp->othertype != 2 || cbp->let_type != 3) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "%ld", (long) (1900 + dp->data [1]));
      }
    } else {
      StringNCpy (year, dp->str, (size_t) 4);
      year [4] = '\0';
    }
  }

  AddValNodeString (&head, "Thesis (", year, ")");

  if (imp->prepub == 2) {
    suffix = ", In press";
  }

  str = NULL;
  afp = imp->pub;
  if (afp != NULL) {
    if (afp->choice == 1) {
      str = StringSave (afp->affil);
    } else if (afp->choice == 2) {
      str = MakeAffilStr (afp);
    }
  }

  if (str != NULL) {

    /* convert double quotes to single quotes */

    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }
    AddValNodeString (&head, " ", str, suffix);
    MemFree (str);
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitArt (
  FmtType format,
  Boolean citArtIsoJta,
  CitArtPtr cap
)

{
  CitBookPtr  cbp;
  CitJourPtr  cjp;
  CharPtr     rsult = NULL;

  if (cap == NULL) return NULL;

  switch (cap->from) {
    case 1 :
      cjp = (CitJourPtr) cap->fromptr;
      if (cjp != NULL) {
        rsult = FormatCitJour (format, citArtIsoJta, cjp);
      }
      break;
    case 2 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBookArt (format, cbp);
      }
      break;
    case 3 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBookArt (format, cbp);
      }
      break;
    default :
      break;
  }

  return rsult;
}

static CharPtr FormatCitPat (
  FmtType format,
  ModType mode,
  CitPatPtr cpp,
  SeqIdPtr seqidp,
  IntAsn2gbJobPtr ajp
)

{
  AffilPtr       afp;
  AuthListPtr    alp;
  IdPatPtr       cit;
  CharPtr        consortium = NULL;
  Char           date [40];
  ValNodePtr     head = NULL;
  Boolean        is_us_pre_grant = FALSE;
  CharPtr        prefix = NULL;
  CharPtr        rsult = NULL;
  SeqIdPtr       sip;
  CharPtr        str;
  CharPtr        suffix = NULL;
  PatentSeqIdPtr psip;
  Int4           pat_seqid = 0;
  Char           buf[10];

  if (cpp == NULL) return NULL;

  if (StringHasNoText (cpp->number) &&
      StringDoesHaveText (cpp->app_number) &&
      StringCmp (cpp->country, "US") == 0 &&
      mode != RELEASE_MODE) {
    for (sip = seqidp; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_PATENT) continue;
      psip = (PatentSeqIdPtr) sip->data.ptrvalue;
      if (psip == NULL) continue;
      cit = psip->cit;
      if (cit == NULL) continue;
      if (StringDoesHaveText (cit->app_number)) {
        is_us_pre_grant = TRUE;
      }
    }
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (is_us_pre_grant) {
      ValNodeCopyStr (&head, 0, "Pre-Grant Patent: ");
      suffix = " ";
    } else {
      ValNodeCopyStr (&head, 0, "Patent: ");
      suffix = " ";
    }
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, "Patent number ");
  }

  if (! StringHasNoText (cpp->country)) {
    AddValNodeString (&head, NULL, cpp->country, suffix);
  }

  if (! StringHasNoText (cpp->number)) {
    if (ajp != NULL && GetWWW (ajp) && StringCmp (cpp->country, "US") == 0) {
      ValNodeCopyStr (&head, 0, "<a href=\"");
      ValNodeCopyStr (&head, 0, link_uspto);
      ValNodeCopyStr (&head, 0, cpp->number);
      ValNodeCopyStr (&head, 0, "\">");
      ValNodeCopyStr (&head, 0, cpp->number);
      ValNodeCopyStr (&head, 0, "</a>");
    } else {
      ValNodeCopyStr (&head, 0, cpp->number);
    }
  } else if (! StringHasNoText (cpp->app_number)) {
    if (is_us_pre_grant) {
      AddValNodeString (&head, NULL, cpp->app_number, NULL);
    } else {
      AddValNodeString (&head, "(", cpp->app_number, ")");
    }
  }

  if (! StringHasNoText (cpp->doc_type)) {
    AddValNodeString (&head, "-", cpp->doc_type, NULL);
  }

  /* pat_seqid test */

  for (sip = seqidp; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_PATENT) {
      psip = (PatentSeqIdPtr) sip -> data.ptrvalue;
      if (psip != NULL) {
        pat_seqid = psip->seqid;
      }
    }
  }
  if (pat_seqid > 0) {
    if (format == EMBL_FMT) {
      sprintf(buf,"%s%ld%s", "/", (long) pat_seqid, ", ");
      ValNodeCopyStr (&head, 0, buf);
    } else {
      sprintf(buf,"%s%ld ", " ", (long) pat_seqid);
      ValNodeCopyStr (&head, 0, buf);
    }
  } else {
    ValNodeCopyStr (&head, 0, " ");
  }

  /* Date */

  date [0] = '\0';
  if (cpp->date_issue != NULL) {
    DateToFF (date, cpp->date_issue, FALSE);
  } else if (cpp->app_date != NULL) {
    DateToFF (date, cpp->app_date, FALSE);
  }
  if (! StringHasNoText (date)) {
    ValNodeCopyStr (&head, 0, date);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ";");
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ".");
  }

  alp = cpp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      suffix = NULL;
      if (afp->choice == 2) {
        suffix = ";";
      }

      /* If any of the affiliation fields are */
      /* non-blank, put them on a new line.   */

      if ((! StringHasNoText (afp->affil)) ||
          (! StringHasNoText (afp->street)) ||
          (! StringHasNoText (afp->div)) ||
          (! StringHasNoText (afp->city)) ||
          (! StringHasNoText (afp->sub)) ||
          (! StringHasNoText (afp->country)))
        ValNodeCopyStr (&head, 0, "\n");

      /* Write out the affiliation fields */

      if (! StringHasNoText (afp->affil)) {
        AddValNodeString (&head, NULL, afp->affil, suffix);
        prefix = " ";
      }
      if (! StringHasNoText (afp->street)) {
        AddValNodeString (&head, prefix, afp->street, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->div)) {
        AddValNodeString (&head, prefix, afp->div, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->city)) {
        AddValNodeString (&head, prefix, afp->city, NULL);
        prefix = ", ";
      }
      if (! StringHasNoText (afp->sub)) {
        AddValNodeString (&head, prefix, afp->sub, NULL);
      }
      if (! StringHasNoText (afp->country)) {
        AddValNodeString (&head, ";\n", afp->country, ";");
      }
    }
  }

  alp = cpp->assignees;
  if (alp != NULL) {
    str = GetAuthorsString (format, alp, &consortium, NULL, NULL);
    afp = alp->affil;
    if (afp != NULL) {
      suffix = NULL;
      if (afp->choice == 2) {
        suffix = ";";
      }

      /* If any of the affiliation fields are */
      /* non-blank, put them on a new line.   */

      if ((! StringHasNoText (str)) ||
          (! StringHasNoText (consortium)) ||
          (! StringHasNoText (afp->affil)) ||
          (! StringHasNoText (afp->street)) ||
          (! StringHasNoText (afp->div)) ||
          (! StringHasNoText (afp->city)) ||
          (! StringHasNoText (afp->sub)) ||
          (! StringHasNoText (afp->country)))
        ValNodeCopyStr (&head, 0, "\n");

      if (! StringHasNoText (str)) {
        AddValNodeString (&head, NULL, str, ";");
        prefix = " ";
      }
      if (! StringHasNoText (consortium)) {
        AddValNodeString (&head, NULL, consortium, ";");
        prefix = " ";
      }

      /* Write out the affiliation fields */

      if (! StringHasNoText (afp->affil)) {
        AddValNodeString (&head, NULL, afp->affil, suffix);
        prefix = " ";
      }
      if (! StringHasNoText (afp->street)) {
        AddValNodeString (&head, prefix, afp->street, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->div)) {
        AddValNodeString (&head, prefix, afp->div, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->city)) {
        AddValNodeString (&head, prefix, afp->city, NULL);
        prefix = ", ";
      }
      if (! StringHasNoText (afp->sub)) {
        AddValNodeString (&head, prefix, afp->sub, NULL);
      }
      if (! StringHasNoText (afp->country)) {
        AddValNodeString (&head, ";\n", afp->country, ";");
      }
    }
    MemFree (consortium);
    MemFree (str);
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  /*
  s_StringCleanup(rsult);
  */

  return rsult;
}

static CharPtr FormatCitGen (
  FmtType format,
  Boolean dropBadCitGens,
  Boolean noAffilOnUnpub,
  CitGenPtr cgp
)

{
  CharPtr      affil = NULL;
  AuthListPtr  alp = NULL;
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  CharPtr      inpress = NULL;
  CharPtr      journal = NULL;
  Char         pages [128];
  CharPtr      prefix = NULL;
  CharPtr      ptr;
  CharPtr      rsult = NULL;
  Char         year [8];

  if (cgp == NULL) return NULL;

  if (cgp->journal == NULL && StringNICmp (cgp->cit, "unpublished", 11) == 0) {
    if (noAffilOnUnpub) {

      /* !!! temporarily put date in unpublished citation for QA !!! */

      if (dropBadCitGens) {
        year [0] = '\0';
        dp = cgp->date;
        if (dp != NULL) {
          if (dp->data [0] == 1) {
            if (dp->data [1] != 0) {
              sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
            }
          } else {
            StringCpy (year, " (");
            StringNCat (year, dp->str, 4);
            StringCat (year, ")");
          }
        }
        AddValNodeString (&head, NULL, "Unpublished", NULL);
        AddValNodeString (&head, NULL, year, NULL);
        rsult = MergeFFValNodeStrs (head);
        ValNodeFreeData (head);
        return rsult;
      }

      /* !!! remove above section once QA against asn2ff is done !!! */

      return StringSave ("Unpublished");
    }

    alp = cgp->authors;
    if (alp != NULL) {
      affil = GetAffil (alp->affil);
      if (! StringHasNoText (affil)) {
        rsult = MemNew ((size_t) StringLen (affil) + (size_t) StringLen (cgp->cit) + 15);
        StringCpy (rsult, "Unpublished ");
        StringCat (rsult, affil);
        TrimSpacesAroundString (rsult);
        return rsult;
      }
    }

    rsult = StringSave (cgp->cit);
    TrimSpacesAroundString (rsult);
    return rsult;
  }

  year [0] = '\0';
  dp = cgp->date;
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  pages [0] = '\0';
  if (cgp->pages != NULL) {
    FixPages (pages, cgp->pages);
  }

  if (cgp->journal != NULL) {
    journal = (CharPtr) cgp->journal->data.ptrvalue;
  }
  if (cgp->cit != NULL) {
    ptr = StringStr (cgp->cit, "Journal=\"");
    if (ptr != NULL) {
      journal = ptr + 9;
    } else if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
               StringNICmp (cgp->cit, "unpublished", 11) == 0) {

      if ((! dropBadCitGens) || journal != NULL) {
        inpress = cgp->cit;
      } else {
        inpress = "Unpublished";
      }
    } else if (StringNICmp (cgp->cit, "Online Publication", 18) == 0 ||
               StringNICmp (cgp->cit, "Published Only in DataBase", 26) == 0 ||
               StringNICmp (cgp->cit, "In press", 8) == 0 ) {
      inpress = cgp->cit;
    } else if (StringNICmp (cgp->cit, "(er) ", 5) == 0) {
      journal = cgp->cit;
    } else if ((! dropBadCitGens) && journal == NULL) {
      journal = cgp->cit;
    }
  }
  if (journal != NULL) {
    journal = StringSave (journal);
    for (ptr = journal, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      if (ch == '=' || ch == '\"') {
        *ptr = '\0';
      }
    }
    ValNodeAddStr (&head, 0, journal);
    prefix = " ";
  }

  if (! StringHasNoText (inpress)) {
    AddValNodeString (&head, prefix, inpress, NULL);
    prefix = " ";
  }

  if (! StringHasNoText (cgp->volume)) {
    AddValNodeString (&head, prefix, cgp->volume, NULL);
  }

  if (! StringHasNoText (pages)) {
    if (format == GENBANK_FMT) {
      AddValNodeString (&head, ", ", pages, NULL);
    } else if (format == EMBL_FMT) {
      AddValNodeString (&head, ":", pages, NULL);
    }
  }

  if (! StringHasNoText (year)) {
    AddValNodeString (&head, NULL, year, NULL);
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitSub (
  FmtType format,
  CitSubPtr csp
)

{
  CharPtr      affil;
  AffilPtr     afp;
  AuthListPtr  alp;
  Char         buf [256];
  Char         date [40];
  ValNodePtr   head = NULL;
  CharPtr      rsult = NULL;

  if (csp == NULL) return NULL;

  date [0] = '\0';
  if (csp->date != NULL) {
    DateToFF (date, csp->date, TRUE);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "\?\?-\?\?\?-\?\?\?\?");
  }

  sprintf (buf, "Submitted (%s)", date);
  ValNodeCopyStr (&head, 0, buf);

  alp = csp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      affil = GetAffil (afp);
      if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
        if (StringNCmp(affil, " to the EMBL/GenBank/DDBJ databases.", 36) != 0) {
          ValNodeCopyStr (&head, 0, " to the EMBL/GenBank/DDBJ databases.\n");
        } else {
          ValNodeCopyStr (&head, 0, " ");
        }
      } else {
        ValNodeCopyStr (&head, 0, " ");
      }
      ValNodeCopyStr (&head, 0, affil);
      MemFree (affil);
    } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
      ValNodeCopyStr (&head, 0, " to the EMBL/GenBank/DDBJ databases.\n");
    }
  }

  rsult = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetPubJournal (
  FmtType format,
  ModType mode,
  Boolean dropBadCitGens,
  Boolean noAffilOnUnpub,
  Boolean citArtIsoJta,
  PubdescPtr pdp,
  CitSubPtr csp,
  SeqIdPtr seqidp,
  IndxPtr index,
  IntAsn2gbJobPtr ajp
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  CitPatPtr        cpp;
  CharPtr          journal = NULL;
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (csp != NULL) {
    return FormatCitSub (format, csp);
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
            if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
              break; /* skip just serial number */
            }
          }
          journal = FormatCitGen (format, dropBadCitGens, noAffilOnUnpub, cgp);
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          journal = FormatCitSub (format, csp);
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            journal = FormatCitArt (format, citArtIsoJta, cap);
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          journal = FormatCitArt (format, citArtIsoJta, cap);
        }
        break;
      case PUB_Book :
      case PUB_Proc :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatCitBook (format, cbp);
        }
        break;
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatThesis (format, cbp);
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          journal = FormatCitPat (format, mode, cpp, seqidp, ajp);
        }
        break;
      default :
        break;
    }

    /* optionally populate indexes for NCBI internal database */

    if (index != NULL && journal != NULL) {

      /* skip non-informative cit-gens */

      if (StringNICmp (journal, "submitted", 8) == 0 ||
          StringNICmp (journal, "unpublished", 11) == 0 ||
          StringNICmp (journal, "Online Publication", 18) == 0 ||
          StringNICmp (journal, "Published Only in DataBase", 26) == 0) {
      } else {
        ValNodeCopyStrToHead (&(index->journals), 0, journal);
      }
    }

    if (journal != NULL) return journal;
  }

  return NULL;
}

static Int4 GetMuid (
  PubdescPtr pdp
)

{
  ArticleIdPtr     aip;
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->uid;
        }
        break;
      case PUB_Muid :
        return vnp->data.intvalue;
      case PUB_Article:
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap!= NULL && cap->ids != NULL) {
          for (aip = cap->ids; aip != NULL; aip = aip->next) {
            if (aip->choice == ARTICLEID_MEDLINE) {
              return aip->data.intvalue;
            }
          }
        }
      default :
        break;
    }
  }

  return 0;
}

static Int4 GetPmid (
  PubdescPtr pdp
)

{
  ArticleIdPtr       aip;
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->pmid;
        }
        break;
      case PUB_PMid :
        return vnp->data.intvalue;
      case PUB_Article:
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap!= NULL && cap->ids != NULL) {
          for (aip = cap->ids; aip != NULL; aip = aip->next) {
            if (aip->choice == ARTICLEID_PUBMED) {
              return aip->data.intvalue;
            }
          }
        }
      default :
        break;
    }
  }

  return 0;
}

NLM_EXTERN CharPtr CleanQualValue (
  CharPtr str
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t' || ch == '"') {
      *dst = ' ';
      dst++;
    } else {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

NLM_EXTERN CharPtr Asn2gnbkCompressSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  Char     last;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      last = ch;
      ch = *ptr;
      if (ch != '\0' && ch < ' ') {
        *ptr = ' ';
        ch = *ptr;
      }
      while (ch != '\0' && last <= ' ' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

NLM_EXTERN CharPtr StripAllSpaces (
  CharPtr str
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == ' ' || ch == '\t') {
    } else {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static CharPtr remarksText [] = {
  "full automatic", "full staff_review", "full staff_entry",
  "simple staff_review", "simple staff_entry", "simple automatic",
  "unannotated automatic", "unannotated staff_review", "unannotated staff_entry",
  NULL
};

static void AddReferenceToGbseq (
  GBSeqPtr gbseq,
  GBReferencePtr gbref,
  CharPtr str,
  RefBlockPtr rbp,
  BioseqPtr bsp
)

{
  Char            buf [32];
  CharPtr         copy;
  ValNodePtr      head = NULL;
  IntRefBlockPtr  irp;
  SeqLocPtr       loc;
  CharPtr         ptr;
  CharPtr         ref;
  SeqLocPtr       slp;
  Int4            start;
  Int4            stop;
  CharPtr         tmp;

  if (gbseq == NULL || gbref == NULL || StringHasNoText (str) || rbp == NULL || bsp == NULL) return;

  copy = StringSave (str);

  /* link in reverse order, to be reversed in slash block */

  gbref->next = gbseq->references;
  gbseq->references = gbref;

  /* now parse or make ASN required default values for remaining fields */

  if (StringNCmp (copy, "REFERENCE   ", 12) == 0) {
    ref = copy + 12;
    ptr = StringStr (ref, "\n  AUTHORS");
    if (ptr == NULL) {
      ptr = StringStr (ref, "\n  CONSRTM");
    }
    if (ptr == NULL) {
      ptr = StringStr (ref, ")\n");
      if (ptr != NULL) {
        ptr++;
      }
    }
    if (ptr != NULL) {
      *ptr = '\0';
      /* gbref->reference = StringSave (ref); */
      sprintf (buf, "%d", (int) rbp->serial);
      gbref->reference = StringSave (buf);
    }
  }

  if (gbref->reference == NULL) {
    gbref->reference = StringSave ("?");
  }

  CleanQualValue (gbref->reference);
  Asn2gnbkCompressSpaces (gbref->reference);

  if (gbref->journal == NULL) {
    gbref->journal = StringSave ("?");
  }

  CleanQualValue (gbref->journal);
  Asn2gnbkCompressSpaces (gbref->journal);

  MemFree (copy);

  if (rbp->sites == 1 || rbp->sites == 2) {
    gbref->position = StringSave ("sites");
  } else if (rbp->sites == 3) {
  } else {
    irp = (IntRefBlockPtr) rbp;
    loc = irp->loc;
    if (loc != NULL) {
      slp = SeqLocFindNext (loc, NULL);
      while (slp != NULL) {
        start = SeqLocStart (slp) + 1;
        stop = SeqLocStop (slp) + 1;
        if (head == NULL) {
          sprintf (buf, "%ld..%ld", (long) start, (long) stop);
        } else {
          sprintf (buf, "; %ld..%ld", (long) start, (long) stop);
        }
        ValNodeCopyStr (&head, 0, buf);
        slp = SeqLocFindNext (loc, slp);
      }
      tmp = MergeFFValNodeStrs (head);
      ValNodeFreeData (head);
      gbref->position = tmp;
    } else {
      start = 1;
      stop = bsp->length;
      sprintf (buf, "%ld..%ld", (long) start, (long) stop);
      gbref->position = StringSave (buf);
    }
  }
}

static Boolean IsCitSub (
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  ValNodePtr  vnp;

  if (csp != NULL) return TRUE;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Sub) return TRUE;
  }
  return FALSE;
}

static void  FF_www_muid(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  Int4 muid
)
{
  Char numbuf[40];
  
  if ( GetWWW(ajp) ) {
    FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    FF_Add_NCBI_Base_URL (ffstring, link_muid);
    sprintf (numbuf, "%ld", (long)muid);
    FFAddTextToString (ffstring, NULL, numbuf, "\">", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    sprintf(numbuf, "%ld", (long)muid);
    FFAddOneString (ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
  }
}

static Uint1 GetJournalPubStatus (PubdescPtr pdp)

{
  CitArtPtr   cap;
  CitJourPtr  cjp;
  ImprintPtr  imp;
  ValNodePtr  vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Article) continue;
    cap = (CitArtPtr) vnp->data.ptrvalue;
    if (cap == NULL) continue;
    if (cap->from != 1) continue;
    cjp = (CitJourPtr) cap->fromptr;
    if (cjp == NULL) continue;
    imp = cjp->imp;
    if (imp == NULL) continue;
    return imp->pubstatus;
  }

  return 0;
}

NLM_EXTERN CharPtr FormatReferenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  SeqMgrAndContext   acontext;
  AnnotDescPtr       adp;
  IntAsn2gbJobPtr    ajp;
  AuthListPtr        alp;
  Asn2gbSectPtr      asp;
  BioseqPtr          bsp;
  Char               buf [150];
  CitArtPtr          cap;
  Char               ch;
  CitJourPtr         cjp;
  Boolean            citArtIsoJta;
  CharPtr            consortium;
  CitPatPtr          cpp;
  CitRetractPtr      crp;
  CitSubPtr          csp = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int4               gibbsq;
  GBReferencePtr     gbref = NULL;
  GBSeqPtr           gbseq;
  GBXrefPtr          gxp;
  ValNodePtr         head;
  Int2               i;
  ArticleIdPtr       ids;
  ImprintPtr         imp;
  IndxPtr            index;
  IntRefBlockPtr     irp;
  size_t             len;
  SeqLocPtr          loc = NULL;
  MedlineEntryPtr    mep;
  Int4               muid = 0;
  Boolean            needsPeriod = FALSE;
  SeqLocPtr          nextslp;
  Boolean            notFound;
  ObjMgrDataPtr      omdp;
  PubdescPtr         pdp = NULL;
  PubdescPtr         pdpcopy = NULL;
  PubmedEntryPtr     pep = NULL;
  Int4               pmid = 0;
  CharPtr            prefix = NULL;
  Uint1              pubstatus;
  CharPtr            pubstatnote;
  RefBlockPtr        rbp;
  ValNodePtr         remarks = NULL;
  CharPtr            remprefix = NULL;
  SubmitBlockPtr     sbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp = NULL;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  SeqSubmitPtr       ssp;
  Int4               start;
  Int4               stop;
  CharPtr            str = NULL;
  Boolean            strict_isojta;
  CharPtr            suffix = NULL;
  BioseqPtr          target;
  CharPtr            tmp;
  Boolean            trailingPeriod = TRUE;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, temp;

  if (afp == NULL || bbp == NULL) return NULL;
  rbp = (RefBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;

  /* five-column feature table uses special code for formatting */

  if (ajp->format == FTABLE_FMT) {
    irp = (IntRefBlockPtr) bbp;
    if (irp->loc != NULL) {
      if (irp->rb.pmid != 0 || irp->rb.muid != 0) {
        head = NULL;
        PrintFtableIntervals (&head, target, irp->loc, "REFERENCE");
        if (irp->rb.pmid != 0) {
          sprintf (buf, "\t\t\tpmid\t%ld\n", (long) irp->rb.pmid);
          ValNodeCopyStr (&head, 0, buf);
        } else if (irp->rb.muid != 0) {
          sprintf (buf, "\t\t\tmuid\t%ld\n", (long) irp->rb.muid);
          ValNodeCopyStr (&head, 0, buf);
        }
        str = MergeFFValNodeStrs (head);
        ValNodeFreeData (head);
      }
    }
    return str;
  }

  /* otherwise do regular flatfile formatting */

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (! StringHasNoText (rbp->string)) return StringSave (rbp->string);

  /* could be descriptor, feature, annotdesc, or submit block citation */

  if (rbp->itemtype == OBJ_SEQDESC) {

    sdp = SeqMgrGetDesiredDescriptor (rbp->entityID, NULL, rbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQFEAT) {

    sfp = SeqMgrGetDesiredFeature (rbp->entityID, NULL, rbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_ANNOTDESC) {

    adp = SeqMgrGetDesiredAnnotDesc (rbp->entityID, NULL, rbp->itemID, &acontext);
    if (adp != NULL && acontext.annotdesctype == Annot_descr_pub) {
      pdp = (PubdescPtr) adp->data.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQSUB_CIT) {

    omdp = ObjMgrGetData (rbp->entityID);
    if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) omdp->dataptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sbp = ssp->sub;
        if (sbp != NULL) {
          csp = sbp->cit;
        }
      }
    }
  }

  if (pdp == NULL && csp == NULL) return NULL;

  temp = FFGetString(ajp);
  if ( temp == NULL ) {
    FFRecycleString(ajp, ffstring);
    return NULL;
  }

  /* any justuids left at this point is RefSeq protein, and should be fetched */

  irp = (IntRefBlockPtr) rbp;
  if (irp->justuids) {
    if (rbp->pmid != 0) {
      pep = GetPubMedForUid (rbp->pmid);
    } else if (rbp->muid != 0) {
      pep = GetPubMedForUid (rbp->muid);
    }
    if (pep != NULL) {
      mep = (MedlineEntryPtr) pep->medent;
      if (mep != NULL && mep->cit != NULL) {
        pdpcopy = AsnIoMemCopy ((Pointer) pdp,
                                 (AsnReadFunc) PubdescAsnRead,
                                 (AsnWriteFunc) PubdescAsnWrite);
        cap = AsnIoMemCopy ((Pointer) mep->cit,
                            (AsnReadFunc) CitArtAsnRead,
                            (AsnWriteFunc) CitArtAsnWrite);
        vnp = ValNodeAddPointer (&(pdpcopy->pub), PUB_Article, (Pointer) cap);
        pdp = pdpcopy;
      }
    }
  }

  /* print serial number */
  FFStartPrint(temp, afp->format, 0, 12, "REFERENCE", 12, 5, 5, "RN", TRUE);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    if (rbp->serial > 99) {
      sprintf (buf, "%d ", (int) rbp->serial);
    } else {
      sprintf (buf, "%d", (int) rbp->serial);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    sprintf (buf, "[%d]", (int) rbp->serial);
  }

  FFAddOneString (temp, buf, FALSE, FALSE, TILDE_TO_SPACES);

  /* print base range */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    if (rbp->sites != 3) {
      FFAddNChar(temp, ' ', 15 - temp->pos, FALSE);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    if (rbp->sites == 0) {
      FFLineWrap(ffstring, temp, 0, 5, ASN2FF_EMBL_MAX, "RN");   
      FFRecycleString(ajp, temp);
      temp = FFGetString(ajp);
      FFStartPrint(temp, afp->format, 0, 0, NULL, 0, 5, 5, "RP", FALSE);
    }
  }

  if (rbp->sites == 1 || rbp->sites == 2) {

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddOneString (temp, "(sites)", FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  } else if (rbp->sites == 3) {
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  } else {
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddNChar(temp, ' ', 15 - temp->pos, FALSE);
      if (afp->format == GENBANK_FMT) {
        FFAddOneString (temp, "(bases ", FALSE, FALSE, TILDE_TO_SPACES);
      } else {
        FFAddOneString (temp, "(residues ", FALSE, FALSE, TILDE_TO_SPACES);
      }
    }

    irp = (IntRefBlockPtr) rbp;
    loc = irp->loc;

    if (loc != NULL) {
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        suffix = "; ";
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        suffix = ", ";
      }

      slp = SeqLocFindNext (loc, NULL);
      while (slp != NULL) {
        nextslp = SeqLocFindNext (loc, slp);
        start = SeqLocStart (slp) + 1;
        stop = SeqLocStop (slp) + 1;
        if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
          sprintf (buf, "%ld to %ld", (long) start, (long) stop);
        } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
          sprintf (buf, "%ld-%ld", (long) start, (long) stop);
        }
        if (nextslp == NULL) {
          suffix = NULL;
        }
        FFAddTextToString (temp, NULL, buf, suffix, FALSE, FALSE, TILDE_TO_SPACES);
        slp = nextslp;
      }

    } else {

      /* code still used for ssp->cit */

      start = 1;
      stop = bsp->length;
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        sprintf (buf, "%ld to %ld", (long) start, (long) stop);
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        sprintf (buf, "%ld-%ld", (long) start, (long) stop);
      }
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_TO_SPACES);
    }

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddOneString (temp, ")", FALSE, FALSE, TILDE_TO_SPACES);
    }
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  }

  if (gbseq != NULL) {
    gbref = GBReferenceNew ();
  }

  /* print author list */

  str = NULL;
  consortium = NULL;

  alp = GetAuthListPtr (pdp, csp);
  if (alp != NULL) {
    str = GetAuthorsString (afp->format, alp, &consortium, index, gbref);
    TrimSpacesAroundString (str);
  }

  if (str != NULL || StringHasNoText (consortium)) {
    FFRecycleString(ajp, temp);
    temp = FFGetString(ajp);
    FFStartPrint(temp, afp->format, 2, 12, "AUTHORS", 12, 5, 5, "RA", FALSE);

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      suffix = NULL;
      trailingPeriod = TRUE;
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      trailingPeriod = FALSE;
      len = StringLen (str);
      if (len > 0 && str [len - 1] != '.') {
        suffix = ".;";
      } else {
        suffix = ";";
      }
    }

    /* if no authors were found, period will still be added by this call */
    if (str != NULL) {
      FFAddTextToString (temp, NULL, str, suffix, trailingPeriod, FALSE, TILDE_TO_SPACES);
    } else if (StringHasNoText (consortium)) {
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        FFAddOneChar(temp, '.', FALSE);
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        FFAddOneChar(temp, ';', FALSE);
      }    
    }

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RA");
    }
  }
  MemFree (str);

  /* print consortium */

  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  if (! StringHasNoText (consortium)) {
    FFStartPrint (temp, afp->format, 2, 12, "CONSRTM", 12, 5, 5, "RG", FALSE);
    FFAddTextToString (temp, NULL, consortium, suffix, FALSE, FALSE, TILDE_TO_SPACES);
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RG");
    }
  }
  MemFree (consortium);

  /* print title */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);

  str = GetPubTitle (afp->format, pdp, csp);
  CleanPubTitle (str);
  StrStripSpaces (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    prefix = NULL;
    suffix = NULL;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    if (str != NULL) {
      prefix = "\"";
      suffix = "\";";
    } else {
      prefix = NULL;
      suffix = ";";
    }
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    if (! StringHasNoText (str)) {
      FFStartPrint (temp, afp->format, 2, 12, "TITLE", 12, 5, 5, "RT", FALSE);

      FFAddTextToString (temp, prefix, str, suffix, FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFStartPrint (temp, afp->format, 2, 12, "TITLE", 12, 5, 5, "RT", FALSE);
    if (! StringHasNoText (str)) {

      FFAddTextToString (temp, prefix, str, suffix, FALSE, FALSE, TILDE_TO_SPACES);

    } else {
      FFAddOneChar (temp, ';', FALSE);
    }
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RT");
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->title = StringSaveNoNull (str);
    }
  }

  MemFree (str);

  /* print journal */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);

  FFStartPrint (temp, afp->format, 2, 12, "JOURNAL", 12, 5, 5, "RL", FALSE);

  /* Only GenBank/EMBL/DDBJ require ISO JTA in ENTREZ/RELEASE modes (RefSeq should later) */

  citArtIsoJta = ajp->flags.citArtIsoJta;
  strict_isojta = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ ||
        /* sip->choice == SEQID_OTHER || */
        sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) {
      strict_isojta = TRUE;
    }
  }
  if (! strict_isojta) {
    citArtIsoJta = FALSE;
  }

  str = GetPubJournal (afp->format, ajp->mode, ajp->flags.dropBadCitGens,
                       ajp->flags.noAffilOnUnpub, citArtIsoJta, pdp, csp,
                       bsp->id, index, ajp);
  if (str == NULL) {
    str = StringSave ("Unpublished");
  }
  StrStripSpaces (str);
  TrimSpacesAroundString (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    needsPeriod = FALSE;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    if (! IsCitSub (pdp, csp)) {
      needsPeriod = TRUE;
    }
  }

  FFAddOneString (temp, str, FALSE, FALSE, TILDE_IGNORE);
  if (needsPeriod) {
    FFAddOneChar(temp, '.', FALSE);
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->journal = StringSaveNoNull (str);
      tmp = gbref->journal;
      if (tmp != NULL) {
        ch = *tmp;
        while (ch != '\0') {
          if (ch == '\n' || ch == '\r' || ch == '\t') {
            *tmp = ' ';
          }
          tmp++;
          ch = *tmp;
        }
        TrimSpacesAroundString (gbref->journal);
      }
    }
  }

  MemFree (str);
  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
  } else {
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RL");
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      if (pdp != NULL) {
        for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == PUB_Article) {
            cap = (CitArtPtr) vnp->data.ptrvalue;
            if (cap != NULL) {
              for (ids = cap->ids; ids != NULL; ids = ids->next) {
                if (ids->choice == ARTICLEID_DOI) {
                  tmp = (CharPtr) ids->data.ptrvalue;
                  if (StringDoesHaveText (tmp)) {
                    gxp = GBXrefNew ();
                    if (gxp != NULL) {
                      gxp->dbname = StringSave ("doi");
                      gxp->id = StringSave (tmp);
                      gxp->next = gbref->xref;
                      gbref->xref = gxp;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  /* print muid */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  
  pmid = GetPmid (pdp);
  muid = GetMuid (pdp);

  if (pmid == 0 && muid > 0) {
    FFStartPrint (temp, afp->format, 2, 12, "MEDLINE", 12, 5, 5, "RX", FALSE);

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FF_www_muid (ajp, temp, muid);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      sprintf (buf, "MEDLINE; %ld.", (long) muid);
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RX");
    }
  }

  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  
  if (pmid > 0) {
    FFStartPrint (temp, afp->format, 3, 12, "PUBMED", 12, 5, 5, "RX", FALSE);
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FF_www_muid (ajp, temp, pmid);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      sprintf (buf, "PUBMED; %ld.", (long) pmid);
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RX");
    }
  }
  FFRecycleString(ajp, temp);

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->pubmed = pmid;
    }
  }

  if (pdp == NULL) {

    if (csp != NULL) {
      if (! StringHasNoText (csp->descr)) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);
 
        ValNodeCopyStr (&remarks, 0, csp->descr);
        FFStartPrint (temp, afp->format, 2, 12, "REMARK", 12, 5, 5, NULL, FALSE);
        /* FFAddOneString (temp, csp->descr, FALSE, TRUE, TILDE_EXPAND); */
        AddCommentWithURLlinks(ajp, temp, NULL, csp->descr, NULL);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
      }
    }

    str = FFToCharPtr(ffstring);

    if (gbseq != NULL) {
      if (gbref != NULL) {
        AddReferenceToGbseq (gbseq, gbref, str, rbp, bsp);
      }
    }

    FFRecycleString(ajp, ffstring);
    FFRecycleString(ajp, temp);
    if (pep != NULL) {
      PubmedEntryFree (pep);
    }
    if (pdpcopy != NULL) {
      PubdescFree (pdpcopy);
    }

    return str;
  }


  /* !!! remainder of fields are only for GenBank !!! */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    prefix = "REMARK";

    cpp = NULL;
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Patent) {
        cpp = (CitPatPtr) vnp->data.ptrvalue;
      }
    }
    if (cpp != NULL && ajp != NULL && ajp->mode == ENTREZ_MODE) {
      if (StringCmp (cpp->country, "US") == 0) {
        if (StringDoesHaveText (cpp->number)) {
          FFRecycleString(ajp, temp);
          temp = FFGetString(ajp);

          sprintf (buf, "CAMBIA Patent Lens: %s %s", cpp->country, cpp->number);
          if (remprefix != NULL) {
            ValNodeCopyStr (&remarks, 0, remprefix);
          }
          ValNodeCopyStr (&remarks, 0, buf);
          remprefix = "; ";
          FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
          if (GetWWW (ajp)) {
            sprintf (buf, "CAMBIA Patent Lens: %s ", cpp->country);
            FFAddOneString (temp, buf, FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, "<a href=\"", FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, link_cambia, FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, cpp->country, FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, cpp->number, FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, "#list\">", FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, cpp->number, FALSE, FALSE, TILDE_EXPAND);
            FFAddOneString (temp, "</a>", FALSE, FALSE, TILDE_EXPAND);
          } else {
            FFAddOneString (temp, buf, FALSE, FALSE, TILDE_EXPAND);
          }
          FFLineWrap (ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
          prefix = NULL;
        }
      }
    }

    if (pdp->comment != NULL) {
      for (i = 0, notFound = TRUE; notFound && remarksText [i] != NULL; i++) {
        if (StringCmp (pdp->comment, remarksText [i]) == 0) {
          notFound = FALSE;
        }
      }
      if (notFound) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        if (remprefix != NULL) {
          ValNodeCopyStr (&remarks, 0, remprefix);
        }
        ValNodeCopyStr (&remarks, 0, pdp->comment);
        remprefix = "; ";
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, pdp->comment, FALSE, TRUE, TILDE_EXPAND);
        /* AddCommentWithURLlinks(ajp, temp, NULL, pdp->comment, NULL); */
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;

        if (gbseq != NULL) {
          if (gbref != NULL) {
            /*
            gbref->remark = StringSave (pdp->comment);
            */
          }
        }

      }
    }

    gibbsq = 0;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GIBBSQ) {
        gibbsq = sip->data.intvalue;
      }
    }
    csp = NULL;
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Sub) {
        csp = (CitSubPtr) vnp->data.ptrvalue;
      }
    }
    if (gibbsq > 0 /* && csp == NULL */) {
      FFRecycleString(ajp, temp);
      temp = FFGetString(ajp);

      sprintf (buf, "GenBank staff at the National Library of Medicine created this entry [NCBI gibbsq %ld] from the original journal article.", (long) gibbsq);
      if (remprefix != NULL) {
        ValNodeCopyStr (&remarks, 0, remprefix);
      }
      ValNodeCopyStr (&remarks, 0, buf);
      remprefix = "; ";
      FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_EXPAND);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
      prefix = NULL;

      /* gibbsq comment section (fields may be copied from degenerate pubdesc) */

      str = pdp->fig;
      if (StringHasNoText (str)) {
        str = irp->fig;
      }
      if (! StringHasNoText (str)) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        sprintf (buf, "This sequence comes from %s", str);
        if (remprefix != NULL) {
          ValNodeCopyStr (&remarks, 0, remprefix);
        }
        ValNodeCopyStr (&remarks, 0, buf);
        remprefix = "; ";
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, buf, TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

      if (pdp->poly_a || irp->poly_a) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        if (remprefix != NULL) {
          ValNodeCopyStr (&remarks, 0, remprefix);
        }
        ValNodeCopyStr (&remarks, 0, "Polyadenylate residues occurring in the figure were omitted from the sequence.");
        remprefix = "; ";
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, "Polyadenylate residues occurring in the figure were omitted from the sequence.", TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

      str = pdp->maploc;
      if (StringHasNoText (str)) {
        str = irp->maploc;
      }
      if (! StringHasNoText (str)) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        sprintf (buf, "Map location: %s", str);
        if (remprefix != NULL) {
          ValNodeCopyStr (&remarks, 0, remprefix);
        }
        ValNodeCopyStr (&remarks, 0, buf);
        remprefix = "; ";
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, buf, TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

    }

    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Article) {
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL && cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
            if (imp != NULL) {
              crp = imp->retract;
              if (crp != NULL) {
                if (crp->type == 1) {
                  FFRecycleString(ajp, temp);
                  temp = FFGetString(ajp);

                  len = StringLen (crp->exp) + 30;
                  str = MemNew (sizeof (Char) * len);
                  if (str != NULL) {
                    StringCpy (str, "Retracted");
                    if (StringDoesHaveText (crp->exp)) {
                      StringCat (str, ":[");
                      StringCat (str, crp->exp);
                      StringCat (str, "]");
                    }
                    if (remprefix != NULL) {
                      ValNodeCopyStr (&remarks, 0, remprefix);
                    }
                    ValNodeCopyStr (&remarks, 0, str);
                    remprefix = "; ";
                    str = MemFree (str);
                  }
                  FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
                  FFAddOneString (temp, "Retracted", FALSE, FALSE, TILDE_TO_SPACES);
                  if (StringDoesHaveText (crp->exp)) {
                    FFAddTextToString (temp, ":[", crp->exp, "]", FALSE, TRUE, TILDE_EXPAND);
                  }
                  FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
                  prefix = NULL;
                } else if (crp->type == 3) {
                  FFRecycleString(ajp, temp);
                  temp = FFGetString(ajp);

                  len = StringLen (crp->exp) + 30;
                  str = MemNew (sizeof (Char) * len);
                  if (str != NULL) {
                    StringCpy (str, "Erratum");
                    if (StringDoesHaveText (crp->exp)) {
                      StringCat (str, ":[");
                      StringCat (str, crp->exp);
                      StringCat (str, "]");
                    }
                    if (remprefix != NULL) {
                      ValNodeCopyStr (&remarks, 0, remprefix);
                    }
                    ValNodeCopyStr (&remarks, 0, str);
                    remprefix = "; ";
                    str = MemFree (str);
                  }
                  FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
                  FFAddOneString (temp, "Erratum", FALSE, FALSE, TILDE_TO_SPACES);
                  if (StringDoesHaveText (crp->exp)) {
                    FFAddTextToString (temp, ":[", crp->exp, "]", FALSE, TRUE, TILDE_EXPAND);
                  }
                  FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
                  prefix = NULL;
                } else if (crp->type == 4) {
                  FFRecycleString(ajp, temp);
                  temp = FFGetString(ajp);

                  len = StringLen (crp->exp) + 30;
                  str = MemNew (sizeof (Char) * len);
                  if (str != NULL) {
                    StringCpy (str, "Correction");
                    if (StringDoesHaveText (crp->exp)) {
                      StringCat (str, " to:[");
                      StringCat (str, crp->exp);
                      StringCat (str, "]");
                    }
                    if (remprefix != NULL) {
                      ValNodeCopyStr (&remarks, 0, remprefix);
                    }
                    ValNodeCopyStr (&remarks, 0, str);
                    remprefix = "; ";
                    str = MemFree (str);
                  }
                  FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
                  FFAddOneString (temp, "Correction", FALSE, FALSE, TILDE_TO_SPACES);
                  if (StringDoesHaveText (crp->exp)) {
                    FFAddTextToString (temp, " to:[", crp->exp, "]", FALSE, TRUE, TILDE_EXPAND);
                  }
                  FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
                  prefix = NULL;
                }
              }
            }
          }
        }
      } else if (vnp->choice == PUB_Sub) {
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (! StringHasNoText (csp->descr)) {
            FFRecycleString(ajp, temp);
            temp = FFGetString(ajp);

            if (remprefix != NULL) {
              ValNodeCopyStr (&remarks, 0, remprefix);
            }
            ValNodeCopyStr (&remarks, 0, csp->descr);
            remprefix = "; ";
            FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
            /* FFAddOneString (temp, csp->descr, FALSE, TRUE, TILDE_EXPAND); */
            AddCommentWithURLlinks(ajp, temp, NULL, csp->descr, NULL);
            FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
            prefix = NULL;
          }
        }
      }
    }

    pubstatnote = NULL;
    pubstatus = GetJournalPubStatus (pdp);
    if (pubstatus == 3) {
      pubstatnote = "Publication Status: Online-Only";
    } else if (pubstatus == 10) {
      pubstatnote = "Publication Status: Available-Online prior to print";
    }
    if (StringDoesHaveText (pubstatnote)) {
      FFRecycleString(ajp, temp);
      temp = FFGetString(ajp);

      if (remprefix != NULL) {
        ValNodeCopyStr (&remarks, 0, remprefix);
      }
      ValNodeCopyStr (&remarks, 0, pubstatnote);
      remprefix = "; ";
      FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
      FFAddOneString (temp, pubstatnote, FALSE, FALSE, TILDE_EXPAND);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
      prefix = NULL;
    }

  }

  str = FFToCharPtr(ffstring);

  if (gbseq != NULL) {
    if (gbref != NULL) {
      if (remarks != NULL) {
        gbref->remark = MergeFFValNodeStrs (remarks);
      }

      AddReferenceToGbseq (gbseq, gbref, str, rbp, bsp);
    }
  }
  ValNodeFreeData (remarks);

  FFRecycleString(ajp, ffstring);
  FFRecycleString(ajp, temp);
  if (pep != NULL) {
    PubmedEntryFree (pep);
  }
  if (pdpcopy != NULL) {
    PubdescFree (pdpcopy);
  }

  return str;
}


