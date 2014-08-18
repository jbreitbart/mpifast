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
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define LIBC	"libc.so.6"
#endif

#include "pio_intercept.h"
#include "mpiblast_writer.hpp"
#include <stdarg.h>
#include <unistd.h>
#ifdef HAVE_DLL
#include <dlfcn.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <assert.h>

#include <codon.h>
#include <ncbimisc.h>
#include <salpacc.h>
#include <salpstat.h>
#include <fdlKludge.h>
#include <intercept_def.h>

#include <iostream>
#include <mpiblast_util.hpp>

#define PARSE_BUF_SIZE 3*1024*1024 
#define BUFFER_LENGTH  2048

extern FILE* dbgfp;

int pio_hijacking = 0;
int hijack_default = 0;

static void *libc_handle;
static char *print_buf; /* parsing buffer */
static int print_buf_size; /* size of parsing buffer */
static char *parse_ptr; 		/* current buffer position */
static ByteStorePtr int_bs; /* store data to be parsed */
static int intercepting; /* indicate whether output is intercepting*/
static size_t (*__fwrite)(const void *ptr, size_t size, size_t nmemb, FILE *stream);

static FILE* int_stream;  /* the stream to be intercepted */
//static int int_query; /* the query to be intercepted */
static ValNodePtr int_orplist; /* store intercepted data */
static QueryOutputInfoPtr int_qoip;  /* query output info */
static int int_format; /* output format */
static int int_HTML;   /* HTML output */
static int int_mode;   /* intercept mode */
static int state=0;    /* processing state */

/* initialize output intercepting setting */
int init_intercept_lib() {
	char *errmsg;
	
	pio_hijacking = 1;
	print_buf_size = PARSE_BUF_SIZE;

	print_buf = (char*) malloc (print_buf_size);
	if(print_buf == NULL) {
		return 1;
	}
	memset(print_buf, 0, print_buf_size);
	parse_ptr = print_buf;

#ifdef HAVE_DLL
	/* get original fwrite address */
	libc_handle = dlopen(LIBC, RTLD_LAZY);
	if(libc_handle==NULL) {
		fprintf(stderr, "Error opening C runtime library: %s. Use CONFIG_LIBC in the configure script to specify the correct libc.\n", dlerror());
		return 1;
	}
	__fwrite = (size_t (*)(const void*, size_t, size_t, FILE*))dlsym(libc_handle, "fwrite");
	errmsg = dlerror();
	if(errmsg!=NULL) {
		fprintf(stderr, "dlsym error %s\n", errmsg);
		return 1;
	}
#endif

	return 0;
}

/* finalize */
int fin_intercept_lib() {
	free (print_buf);
#ifdef HAVE_DLL
	dlclose(libc_handle);
#endif
	return 0;
}

/* hijack fprintf call, redirect it to the memory buffers */
/* in order to prevent print_buf overflow, we can intercept write() call instead */
int fprintf(FILE *stream, const char *format, ...) {
	va_list arg, arg1, arg2;
	int done;
	int print_len;
	char errmsg[100];
	
	va_start(arg, format);
	
	/* make a copy of arg, because in some system arg will be changed after calling ?printf function*/
	va_copy(arg1, arg);

	if(intercepting && (stream == int_stream) ) { /* we want to intercept this stream */
		done = vsprintf(print_buf, format, arg);
		
		if(print_buf[print_buf_size - 1] != 0) { /* print_buf overflow */
			memset(errmsg, 0, 100);
			sprintf(errmsg, "fprintf(): print buffer overflow\n");
#ifdef HAVE_DLL
			__fwrite(errmsg, 100, sizeof(char), stderr);
#else
			write(fileno(stderr), errmsg, 100 * sizeof(char));
#endif
			exit(-1);
		}
		
		print_len = strlen(print_buf);
		/* write intercepted data to the buffer */
		BSWrite(int_bs, print_buf, print_len);
	}

	va_end(arg);

	/* continue original fprintf */
	if(stream != int_stream) {
		done = vfprintf(stream, format, arg1);
	}
	va_end(arg1);

	return done;
}

/* hijack fwrite call, redirect it to the memory buffers */
size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream) {
	
	if(intercepting && stream == int_stream) {
		/* write intercepted data to the buffer */
		BSWrite(int_bs, (void*)ptr, size*nmemb);
		return nmemb;
	} else {
#ifdef HAVE_DLL
		return __fwrite(ptr, size, nmemb, stream);
#else
		return write(fileno(stream), ptr, size * nmemb);
#endif
	}
}

/* Create a list of output records from a list of sequence alginments. This list will be
 * used to store intercepted outputs.
 */
int create_output_records(SeqAlignPtr seqalign, ValNodePtr PNTR orplist, int maxnum) {
	char buffer[BUFFER_LENGTH+1];
    BioseqPtr bsp;
    SeqIdPtr bestid, gi_list, subject_id, sip_list=NULL, last_id;
	Boolean same_id, found_score;
    Boolean retval = FALSE;
	int numalign = 0;
    Nlm_FloatHi bit_score, evalue;
    Int4 gi = 0, number, score, copy_len;
	OutputRecordPtr curr_orp;
    ValNodePtr last_node=NULL, curr_node=NULL;
	Boolean is_first=TRUE;
	int count=0;

	last_node = NULL;
	*orplist = NULL;
	curr_orp = NULL;
	last_id = NULL;
    while (seqalign) {
            
		subject_id = SeqIdDup(TxGetSubjectIdFromSeqAlign(seqalign));

		same_id = FALSE;
		if(last_id && SeqIdComp(subject_id, last_id) == SIC_YES) {
			same_id = TRUE;
		}
		
		last_id = SeqIdFree(last_id);
		last_id = SeqIdDup(subject_id);
		subject_id = SeqIdFree(subject_id);
		
		found_score = GetScoreAndEvalue(seqalign, &score, &bit_score, &evalue, &number);
		/* if the ID has been seen before, check that proper values are saved. */
		if (same_id == TRUE) {
			if (bit_score > curr_orp->bit_score)
				curr_orp->bit_score = bit_score;
			if (evalue < curr_orp->evalue)
				curr_orp->evalue = evalue;
		} else {

			if(++count > maxnum) {
				break;
			}
			
			/* add new output record in the list */
			if(curr_orp) {
				curr_node=ValNodeAddPointer(&last_node, 0, curr_orp);
				if(is_first) {
					*orplist = last_node;
					is_first = FALSE;
				}
				last_node = curr_node;
			}
			
			curr_orp = (OutputRecordPtr)malloc(sizeof(OutputRecord));
			CHECK_NULPTR(curr_orp);

			curr_orp->id = NULL;
			curr_orp->des_size = 0;
			curr_orp->des_data = NULL;
			curr_orp->aln_size = 0;
			curr_orp->aln_data = NULL;
			curr_orp->des_offset = -1;
			curr_orp->aln_offset = -1;

			curr_orp->bit_score = bit_score;
			curr_orp->evalue = evalue;
			
			retval = TRUE;
		}
		seqalign = seqalign->next;
		numalign++;
	}
	if(curr_orp) {
		curr_node=ValNodeAddPointer(&last_node, 0, curr_orp);
		if(is_first) {
			*orplist = last_node;
			is_first = FALSE;
		}
	}
    last_id = SeqIdFree(last_id);

	return retval;
}

/* debug function, print out a list of output records */
void print_output_records(FILE* fp, ValNodePtr orplist) {
	OutputRecordPtr curr_orp;
	int count;
	
	count=-1;
	while(orplist) {
		curr_orp = (OutputRecordPtr)(orplist->data.ptrvalue);
		fprintf(fp, "\nrecord %d:\n", ++count);
		fflush(fp);
		print_output_record(fp, curr_orp);
		orplist = orplist->next;
	}
	return;
}

/* debug function, print out a output record */
void print_output_record(FILE *fp, OutputRecordPtr curr_orp) {
	fprintf(fp, "\nrank = %d, fragment=%d\n", curr_orp->rank, curr_orp->frag_id);
	fprintf(fp, "bit_score = %lf, evalue=%lf\n", curr_orp->bit_score, curr_orp->evalue);
	fprintf(fp, "description size = %d\n", curr_orp->des_size);
	fflush(fp);
	if(curr_orp->des_size>0 && curr_orp->des_data!=NULL) {
		fwrite(curr_orp->des_data, curr_orp->des_size, 1, fp);
	}
	fprintf(fp, "alignment size = %d\n", curr_orp->aln_size);
	fflush(fp);
	if(curr_orp->aln_size>0 && curr_orp->aln_data!=NULL) {
		fwrite(curr_orp->aln_data, curr_orp->aln_size, 1, fp);
	}
	fprintf(fp, "des_offset = %ld\n", curr_orp->des_offset);
	fprintf(fp, "aln_offset=%ld\n", curr_orp->aln_offset);
	
	fflush(fp);
}

/* start intercepting output for a query */
int start_intercept(FILE *stream, int format, Boolean HTML, int mode) {
	int_stream = stream;
	int_mode = mode;
	int_format = format;
	int_HTML = HTML;
	intercepting = 1;
	parse_ptr = print_buf;
	
	int_bs = BSNew(PARSE_BUF_SIZE);
	int_qoip = new_query_output_infop();

	return 0;
}

/* end intercepting output for a query */
int end_intercept(int query_id, int frag_id, ByteStorePtr bs_ptr) {
	int ret = 0;

	intercepting = 0;
	state = 0;

	if(int_format == 7) {
		ret = parse_output_xml();
	}
	
	/* free int_bs */
	if(int_bs!=NULL) {
		BSFree(int_bs);
	}
	
	if(int_mode == INT_INFO_RECORDS) {
		/* add current intercepted outputs to local list */
		cProcessLocalOutputs(query_id, frag_id, int_orplist, int_qoip, bs_ptr);
	} else if (int_mode == INT_INFO) {
		if(hijack_default) {
			cUpdateDefaultReportInfo(int_qoip);
		} else {
			cUpdateQueryOutputInfo(query_id, int_qoip);
		}
	}
	
	/* free int_orplist */
	int_orplist = ValNodeFree(int_orplist);
	
	/* free int_qoip */
	free_query_output_infop(int_qoip);

	if(ret == -1) {
		cerr << "WARNING: results parsing errors" << endl;
	}

	return ret;
}

/* start intercepting the output stream */
int start_comm_intercept(FILE *stream) {
	int_stream = stream;
	intercepting = 1;
	
	int_bs = BSNew(PARSE_BUF_SIZE);
}

/* end intercepting the output stream and return content in buffer */
int end_comm_intercept(char** buffer) {
	int copy_len = 0;

	intercepting = 0;

	copy_len = BSTell(int_bs);
	/* reset buffer */
	BSSeek(int_bs, 0, SEEK_SET);

	if(copy_len > 0) {
		*buffer = (char*) malloc (copy_len);
		CHECK_NULPTR(*buffer);
		BSRead(int_bs, *buffer, copy_len);
		BSSeek(int_bs, 0, SEEK_SET);
	} 

	if(int_bs != NULL) {
		BSFree(int_bs);
	}
			
	return copy_len;
}

QueryOutputInfoPtr new_query_output_infop(void) {
	QueryOutputInfoPtr qoip;
	
	qoip = (QueryOutputInfoPtr)malloc(sizeof(QueryOutputInfo));
	CHECK_NULPTR(qoip);
	qoip->header = NULL;
	qoip->header_size = 0;
	qoip->des_hdr = NULL;
	qoip->des_hdr_size = 0;
	qoip->des_ftr = NULL;
	qoip->des_ftr_size = 0;
	qoip->no_hits = NULL;
	qoip->no_hits_size = 0;
	qoip->stat = NULL;
	qoip->stat_size = 0;
	qoip->footer = NULL;
	qoip->footer_size = 0;
	
	return qoip;
}

int free_query_output_infop(QueryOutputInfoPtr qoip) {
	if(qoip->header != NULL)	free(qoip->header);
	if(qoip->des_hdr != NULL)	free(qoip->des_hdr);
	if(qoip->des_ftr != NULL)	free(qoip->des_ftr);
	if(qoip->no_hits != NULL)	free(qoip->no_hits);
	if(qoip->stat != NULL)		free(qoip->stat);
	if(qoip->footer != NULL)	free(qoip->footer);
	free(qoip);
	qoip = NULL;
	
	return 0;
}

DefaultReportPtr new_default_reportp(void) {
	DefaultReportPtr drp;

	drp = (DefaultReportPtr)malloc(sizeof(DefaultReport));
	CHECK_NULPTR(drp);
	drp->ref = NULL;
	drp->ref_size = 0;
	drp->db_info = NULL;
	drp->db_info_size = 0;
	drp->no_hits = NULL;
	drp->no_hits_size = 0;
	drp->db_report = NULL;
	drp->db_report_size = 0;
	drp->footer = NULL;
	drp->db_report_size = 0;

	return drp;
}

int free_default_reportp(DefaultReportPtr drp) {
	if(drp->ref != NULL) 		free(drp->ref);
	if(drp->db_info != NULL) 	free(drp->db_info);
	if(drp->no_hits != NULL) 	free(drp->no_hits);
	if(drp->db_report != NULL) 	free(drp->db_report);
	if(drp->footer != NULL) 	free(drp->footer);
	free(drp);

	return 0;
}

/* must be called before intercepting output records */
int set_results_info(ValNodePtr orplist, int effective_len) {
	int_orplist = orplist;
	return 0;
}

/* intercept output for pairwise format (m=0) */
int intercept_output_pairwise(int tag) {
	static ValNodePtr curr_node;
	static long last_pos;
	long adjust_len;

	OutputRecordPtr curr_orp;
	int print_len, copy_len;
	int copying = 0;

	/* processing state machine */
	switch(state) {
	case 0: /* just started */
		if(tag == INT_HDR_BEGIN) {
			state=1;
		}
		break;
	case 1: /* record search info */
		if(tag == INT_HDR_END) {
			copy_len = BSTell(int_bs);
			/* reset buffer */
			BSSeek(int_bs, 0, SEEK_SET);
			
			if(copy_len > 0) {
				/* copy header output here */
				int_qoip->header_size = copy_len;
				int_qoip->header = (char*) malloc (copy_len);
				CHECK_NULPTR(int_qoip->header);
				BSRead(int_bs, int_qoip->header, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			}
			
			/* skip output records if only collect info */
			if(int_mode == INT_INFO) {
				state = 6;
			}
		} else if (tag == INT_DES_BEGIN) { 
			state = 2;
			BSSeek(int_bs, 0, SEEK_SET);
		} 
		break;
	case 2: /* record description header*/
		if(tag == INT_DES_TAG) {
			copy_len = BSTell(int_bs);
			/* reset int_bs */
			BSSeek(int_bs, 0, SEEK_SET);

			if(copy_len > 0) {
				int_qoip->des_hdr_size = copy_len;
				int_qoip->des_hdr = (char*)malloc(copy_len);
				CHECK_NULPTR(int_qoip->des_hdr);
				BSRead(int_bs, int_qoip->des_hdr, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			}
		
			curr_node = int_orplist;
			state = 3;
		}
		break;
	case 3: /* record description data */
		if(tag == INT_DES_TAG) {
			assert(curr_node);
			curr_orp = (OutputRecordPtr)(curr_node->data.ptrvalue);
			copy_len = BSTell(int_bs);
			curr_orp->des_size = copy_len;
			curr_orp->des_data = (char*)malloc(copy_len);
			CHECK_NULPTR(curr_orp->des_data);
			BSSeek(int_bs, 0, SEEK_SET);
			BSRead(int_bs, curr_orp->des_data, copy_len);
			BSSeek(int_bs, 0, SEEK_SET);
			
			curr_node = curr_node->next;
		} else if (tag == INT_DES_END) {
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);
			/* intercept description tail here*/
			if(copy_len > 0) {
				int_qoip->des_ftr_size = copy_len;
				int_qoip->des_ftr = (char*)malloc(copy_len);
				CHECK_NULPTR(int_qoip->des_ftr);
				BSRead(int_bs, int_qoip->des_ftr, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			} 

			state = 4;
		}
		break;
	case 4:
		if(tag == INT_ALN_TAG) {
			curr_node = int_orplist;
			BSSeek(int_bs, 0, SEEK_SET);
		} else if (tag == INT_ALN_FIRST) {
			state = 5;
		} 
		break;
	case 5: /* record alignment data */
		if(tag == INT_ALN_TAG) {
			last_pos = BSTell(int_bs);
		} else if (tag == INT_ALN_FIRST) {
			copying = 1;   
		} else if(tag == INT_ALN_END) {
			copying = 1;
			last_pos =  BSTell(int_bs);
			state = 6;
		}

		if(copying) {
			assert(curr_node);
			curr_orp = (OutputRecordPtr)(curr_node->data.ptrvalue);
			copy_len = last_pos;
			adjust_len = BSTell(int_bs) - last_pos;
			curr_orp->aln_size = copy_len;
			curr_orp->aln_data = (char*)malloc(copy_len);
			CHECK_NULPTR(curr_orp->aln_data);
			BSSeek(int_bs, 0, SEEK_SET);
			BSRead(int_bs, curr_orp->aln_data, copy_len);
			BSSeek(int_bs, 0, SEEK_SET);
			
			curr_node = curr_node->next;

			/* move data since last segment to the head of buffer */
			
			BSDelete(int_bs, last_pos);
			if(adjust_len > 0) {
				BSSeek(int_bs, adjust_len, SEEK_SET);
			}
		}

		break;
	case 6: // collect statistics
		if(tag == INT_STA_BEGIN) {
			if(int_mode == INT_INFO) {
				/* copy no hits promt */
				copy_len = BSTell(int_bs);
				BSSeek(int_bs, 0, SEEK_SET);
				/* intercept description tail here*/
				if(copy_len > 0) {
					int_qoip->no_hits_size = copy_len;
					int_qoip->no_hits = (char*)malloc(copy_len);
					CHECK_NULPTR(int_qoip->no_hits);
					BSRead(int_bs, int_qoip->no_hits, copy_len);
					BSSeek(int_bs, 0, SEEK_SET);
				} 
			} else {
				BSSeek(int_bs, 0, SEEK_SET);
			}
		} else if (tag == INT_STA_END) {
			/* copy statistics here */
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);
			/* intercept description tail here*/
			if(copy_len > 0) {
				int_qoip->stat_size = copy_len;
				int_qoip->stat = (char*)malloc(copy_len);
				CHECK_NULPTR(int_qoip->stat);
				BSRead(int_bs, int_qoip->stat, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			} 
			
			state = 7;
		}
	
		break;
	case 7: /* footer */
		if(tag == INT_FTR_BEGIN) {
			BSSeek(int_bs, 0, SEEK_SET);
		} else if (tag == INT_FTR_END) {
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);
			/* intercept description tail here*/
			if(copy_len > 0) {
				int_qoip->footer_size = copy_len+1;
				int_qoip->footer = (char*)malloc(copy_len+1);
				CHECK_NULPTR(int_qoip->footer);
				BSRead(int_bs, int_qoip->footer, copy_len);
				int_qoip->footer[copy_len] = '\n';
				BSSeek(int_bs, 0, SEEK_SET);
			} 
			
			state = 8;
		}
	
		break;
	default:
		break;
	}

	return 0;
}

/* intercept output for tabular format (m=8,9))*/
int intercept_output_tabular(int tag) {
	static ValNodePtr curr_node;
	static long last_pos;

	OutputRecordPtr curr_orp;
	int print_len, copy_len;
	int copying = 0;

	/* processing state machine */
	switch(state) {
	case 0: /* just started */
		if(tag == INT_HDR_BEGIN) {
			state=1;
		}
		break;
	case 1: /* record search info */
		if(tag == INT_HDR_END) {
			copy_len = BSTell(int_bs);
			/* capture header here */
			BSSeek(int_bs, 0, SEEK_SET);
			
			if(copy_len > 0) {
				/* copy header output here */
				int_qoip->header_size = copy_len;
				int_qoip->header = (char*) malloc (copy_len);
				CHECK_NULPTR(int_qoip->header);
				BSRead(int_bs, int_qoip->header, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			}
			
			if(int_mode == INT_INFO) 
				state = 4;
				
		} else if ( tag == INT_ALN_BEGIN ) { 
			state = 2;
			BSSeek(int_bs, 0, SEEK_SET);
		} 
		
		break;
	case 2:
		if(tag == INT_ALN_TAG) {
			curr_node = int_orplist;
			BSSeek(int_bs, 0, SEEK_SET);
			state = 3;
		}
		break;
	case 3: /* record segment data */
		if(tag == INT_ALN_TAG || tag == INT_ALN_END) {
			assert(curr_node);
			curr_orp = (OutputRecordPtr)(curr_node->data.ptrvalue);
			copy_len = BSTell(int_bs);
			
			curr_orp->aln_size = copy_len;
			curr_orp->aln_data = (char*)malloc(copy_len);
			CHECK_NULPTR(curr_orp->aln_data);
			BSSeek(int_bs, 0, SEEK_SET);
			BSRead(int_bs, curr_orp->aln_data, copy_len);
			
			curr_node = curr_node->next;
			BSSeek(int_bs, 0, SEEK_SET);
		} 
		
		if(tag == INT_ALN_END) {
			state = 4;
		}

		break;
	case 4: /* footer */
		if(tag == INT_FTR_BEGIN) {
			BSSeek(int_bs, 0, SEEK_SET);
		} else if (tag == INT_FTR_END) {
			/* copy statistics here */
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);
			/* intercept description tail here*/
			if(copy_len > 0) {
				int_qoip->footer_size = copy_len + 1;
				int_qoip->footer = (char*)malloc(copy_len + 1);
				CHECK_NULPTR(int_qoip->footer);
				BSRead(int_bs, int_qoip->footer, copy_len);
				int_qoip->footer[copy_len] = '\n';
				BSSeek(int_bs, 0, SEEK_SET);
			} 
			
			state = 5;
		}
	
		break;
	default:
		break;
	}
	return 0;
}

/* intercept ouput for asn format (m=10,11)*/
int intercept_output_asn(int tag) {
	int copy_len, i;
	char byte;
	static ValNodePtr curr_node;
	OutputRecordPtr curr_orp;
	
	switch(state) {
	case 0:
		if( tag == INT_HDR_BEGIN ) {
			state = 1;
		}
		break;
	case 1:
		/* copy header here */
		if(tag == INT_HDR_END) {
			copy_len = BSTell(int_bs);
			if(copy_len > 0) {
				BSSeek(int_bs, 0, SEEK_SET);
				int_qoip->header_size = copy_len;
				int_qoip->header = (char*)malloc(copy_len);
				CHECK_NULPTR(int_qoip->header);
				BSRead(int_bs, int_qoip->header, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			}
				
			curr_node = int_orplist;
			state = 2;
		}
		break;
	case 2:
		/* copy alignment data */
		if(tag == INT_ALN_END || tag == INT_FTR_BEGIN) {
			assert(curr_node);
			curr_orp = (OutputRecordPtr)(curr_node->data.ptrvalue);
			
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);

			if(int_format == 10) { /* text asn output, take care of open struct */
				for(i=0; i<10; i++) { /* use 10 as scan upper bound to avoid dead loop incase something wrong happed */
					byte = BSGetByte(int_bs);
					if(byte == '{') {
						break;
					}
				}
				/* adjust copy_len here */
				copy_len -= (i+1);
			}

			curr_orp->aln_size = copy_len;
			curr_orp->aln_data = (char*)malloc(copy_len);
			CHECK_NULPTR(curr_orp->aln_data);
			BSRead(int_bs, curr_orp->aln_data, copy_len);
			BSSeek(int_bs, 0, SEEK_SET);
			
			curr_node = curr_node->next;
		} 
		
		if(tag == INT_FTR_BEGIN) { /* the last alignment */
			state = 3;	
		}
		break;
	case 3:
		if(tag == INT_FTR_END) {
			/* copy footer */
			copy_len = BSTell(int_bs);
			BSSeek(int_bs, 0, SEEK_SET);
			/* intercept description tail here*/
			if(copy_len > 0) {
				int_qoip->footer_size = copy_len;
				int_qoip->footer = (char*)malloc(copy_len);
				CHECK_NULPTR(int_qoip->footer);
				BSRead(int_bs, int_qoip->footer, copy_len);
				BSSeek(int_bs, 0, SEEK_SET);
			} 
			state = 4;
		}
		
		break;
	}
	
	return 0;
}

/* parsing outupt for xml format, m=11 */
int parse_output_xml() {
	ValNodePtr curr_node;
	OutputRecordPtr curr_orp;
	char* find_ptr;
	int copy_len;

	char* parse_buf;
	int parse_buf_len;
	
	parse_buf_len = BSLen(int_bs);
	parse_buf = (char*)malloc(parse_buf_len + 1);
	parse_buf[parse_buf_len] = 0;
	CHECK_NULPTR(parse_buf);

	BSMerge(int_bs, parse_buf);
	BSFree(int_bs);
	int_bs = NULL;
	
	curr_node = int_orplist;
	
	parse_ptr = parse_buf;

	/* header */
	find_ptr = strstr(parse_ptr, "</Iteration_iter-num>");
	assert(find_ptr);
	copy_len = find_ptr - parse_ptr + strlen("</Iteration_iter-num>\n");
	int_qoip->header_size = copy_len;
	int_qoip->header = (char*)malloc(copy_len);
	CHECK_NULPTR(int_qoip->header);
	memcpy(int_qoip->header, parse_ptr, copy_len);
	parse_ptr += copy_len;
	
	if(int_mode == INT_INFO_RECORDS) {
		/* description header */
		find_ptr = strstr(parse_ptr, "<Iteration_hits>");
		assert(find_ptr);
		copy_len = find_ptr - parse_ptr + strlen("<Iteration_hits>\n");
		int_qoip->des_hdr_size = copy_len;
		int_qoip->des_hdr = (char*)malloc(copy_len);
		CHECK_NULPTR(int_qoip->des_hdr);
		memcpy(int_qoip->des_hdr, parse_ptr, copy_len);
		
		while((find_ptr=strstr(parse_ptr, "</Hit_num>"))!=NULL) {
			find_ptr += strlen("</Hit_num>\n");
			parse_ptr=strstr(find_ptr, "</Hit>");
			if(parse_ptr==NULL) {
				// fprintf(stderr, "The xml output is corrupted.\n");
				return -1;
			}
			
			if(curr_node == NULL) {
				if(debug_msg) {
					LOG_MSG << "WARNING: parse_output_xml(): mismatched alignments and output" << endl;
				}
				break;
			}

			copy_len = parse_ptr + 7 - find_ptr; /* 7=strlen("</Hit>\n") */
			curr_orp = (OutputRecordPtr)(curr_node->data.ptrvalue);
			assert(curr_orp);
			curr_orp->aln_size = copy_len;
			curr_orp->aln_data = (char*)malloc(copy_len);
			CHECK_NULPTR(curr_orp->aln_data);
			memcpy(curr_orp->aln_data, find_ptr, copy_len);
			
			curr_node = curr_node->next;
		}
		
		/* description footer */
		parse_ptr += 7; /* 7=strlen("</Hit>\n") */
		find_ptr = strstr(parse_ptr, "</Iteration_hits>");
		copy_len = find_ptr - parse_ptr + strlen("</Iteration_hits>\n");
		int_qoip->des_ftr_size = copy_len;
		int_qoip->des_ftr = (char*)malloc(copy_len);
		CHECK_NULPTR(int_qoip->des_ftr);
		memcpy(int_qoip->des_ftr, parse_ptr, copy_len);
	} 
	
	/* stat */
	find_ptr = strstr(parse_ptr, "</Iteration_stat>");
	if(find_ptr != NULL) {
		copy_len = find_ptr - parse_ptr + strlen("</Iteration_stat>\n");
		int_qoip->stat_size = copy_len;
		int_qoip->stat = (char*)malloc(copy_len);
		CHECK_NULPTR(int_qoip->stat);
		memcpy(int_qoip->stat, parse_ptr, copy_len);
		parse_ptr += copy_len;
	}
	
	if(int_mode == INT_INFO) {
		find_ptr = strstr(parse_ptr, "</Iteration_message>");
		copy_len = find_ptr - parse_ptr + strlen("</Iteration_message>\n");
		int_qoip->no_hits_size = copy_len;
		int_qoip->no_hits = (char*)malloc(copy_len);
		CHECK_NULPTR(int_qoip->no_hits);
		memcpy(int_qoip->no_hits, parse_ptr, copy_len);
		parse_ptr += copy_len;
	}
	
	/* footer */
	find_ptr = strstr(parse_ptr, "</BlastOutput>");
	copy_len = find_ptr - parse_ptr + strlen("</BlastOutput>\n");
	int_qoip->footer_size = copy_len;
	int_qoip->footer = (char*)malloc(copy_len);
	CHECK_NULPTR(int_qoip->footer);
	memcpy(int_qoip->footer, parse_ptr, copy_len);
	
	free(parse_buf);
	
	return 0;
}

/* compare scores of two output records */
int OutputRecordCmp(OutputRecordPtr lorp, OutputRecordPtr rorp) {
	int ret;
	
	if(lorp->evalue < rorp->evalue) 
		return -1;
	if(lorp->evalue > rorp->evalue) 
		return 1;
	if(lorp->bit_score > rorp->bit_score)
		return -1;
	if(lorp->bit_score < rorp->bit_score)
		return 1;

	if(lorp->frag_id > rorp->frag_id)
		return -1;
	if(lorp->frag_id < rorp->frag_id)
		return 1;

	return 0;
}

