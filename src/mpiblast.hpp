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
#ifndef __MPIBLAST_HPP__
#define __MPIBLAST_HPP__

#include "mpi.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <queue>
#include <list>

#ifdef __GNUG__
//#include "unistd.h"
#endif
//#include "stdio.h"

#include "file_util.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_config.hpp"
#include "mpiblast_types.h"
#include "fragment_list.hpp"
#include "blastjob.hpp"
#include "blast_wrapper.hpp"
#include "mpiblast_querydata.hpp"
#include "mpiblast_scheduler.hpp"
#include "mpiblast_writer.hpp"

extern "C" {
#include "blastdef.h"
#include "ncbi.h"
#include "objseq.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"
}

#include "mpiblast_tags.hpp"

/**
 * This class stores information related to each query--
 * A BlastJob to track how much of the query has completed
 * and the bioseq_cache for bioseqs that go along with the query results
 */
/*
//class QueryData {
//public:
//	QueryData( BlastJob* bj, ValNodePtr bs_cache ) :
//		blast_job( bj ), bioseq_cache( bs_cache ),
//		searched_frags( 0 ) {};
//	QueryData( const QueryData& qd ) :
//	blast_job( qd.blast_job ),
//	bioseq_cache( qd.bioseq_cache ),
//	searched_frags( qd.searched_frags ),
//	bioseq_lengths( qd.bioseq_lengths ),
//	bioseq_starts( qd.bioseq_starts ),
//	bioseq_strands( qd.bioseq_strands )
//	{}
//	QueryData& operator=( const QueryData& qd ){
//		blast_job = qd.blast_job;
//		bioseq_cache = qd.bioseq_cache;
//		searched_frags = qd.searched_frags;
//		bioseq_lengths = qd.bioseq_lengths;
//		bioseq_starts = qd.bioseq_starts;
//		bioseq_strands = qd.bioseq_strands;
//		return *this;
//	}
//
//	/** A BlastJob tracks which segments this query has been searched against */
//	BlastJob* blast_job;
//	int searched_frags;	/**< counts the number of fragments that have been searched */
//
//	// writer only data
//	/** 
//	 * A cache for any Bioseqs that will be needed during result output.
//	 * These get sent from the worker nodes that performed the search.
//	 */
//	ValNodePtr bioseq_cache;
//	std::vector< int > bioseq_lengths;	/**< a list of the original bioseq lengths */
//	std::vector< int > bioseq_starts;	/**< a list of the original bioseq starts */
//	std::vector< int > bioseq_strands;	/**< a list of the original bioseq strands */
//};

typedef struct ResultData 
{
	/** the search complete notification message sent by the worker 
	 *  The message is an array of integers of the form:
	 *  result_id, frag_id, searched count, searched list, success count, success list, data size
	 *  result_id:	a unique identifier the worker has assigned to this result set.
	 *				This ID gets used to request the results from the worker
	 *  frag_id:	The id of the database fragment that was searched
	 *  searched count:	The number of queries searched
	 *  searched list:	The list of searched queries
	 *  success count:	The number of queries that had hits
	 *  success list:	The list of queries with hits
	 *  data size:	The size of the result buffer in bytes
	 */
	int* notification;
	/** A buffer to store result data */
	unsigned char* buf;
	/** The MPI_Request object tracking the progress of the MPI_Irecv() operation to
	 *  receive the results */
	MPI_Request* req;
} ResultData_t;

/**
 * The MpiBlast class provides the functionality to perform BLAST
 * searches on a cluster of nodes connected with MPI.
 */
class MpiBlast{
	public:
		MpiBlast()
		{ 
			log_stream = &std::cerr;  
			remove_db = false;
			output_format = -1;
			copy_via_set = false; 
			copy_via = COPY_VIA_CP;	// default to standard filesystem copy
			disable_mpi_db = false;	// default to workers sending bioseqs
			concurrent_accesses = 1;
			query_batch = 1;
			result_id = 0;
			
			ncbiblast = NULL;
			myscheduler = NULL;
			writer_m = NULL;
			writer_w = NULL;
			
			ncbiblast = new BLAST();
			_irm = new IrecvMgr();
		}
		
		~MpiBlast();
//
// begin Writer-specific functions
//
		/**
		 * Writer - Reads results in ASN.1 format from a worker
		 */
		void receiveResults( ResultData_t& rd, std::vector<SeqAlignPtr>& results );
		
		/**
		 * Writer - Sends all DB fragment files for "fragment_id" to worker 
		 * "worker_id" using MPI
		 */
		bool sendFragsViaMPI(int fragment_id, int worker_id);
                
		/**
		 * Writer - process received BLAST results
		 * scan through the list of pending result receive operations
		 * when an operation has completed, merge the results with any
		 * previously received results for those queries
		 */
		void processResultReceiveCompletions();
		void readPartialBioseqs( SeqAlignPtr sap_p, AsnIoMemPtr aimp, int queryI );

		/**
		 * Writer - check whether a query has been completely searched, and if so,
		 *          write it using the NCBI toolbox code
		 */
		void outputQueryResults();

//
// begin Worker-specific functions
//
		/**
		 * Worker -  Receives all DB fragment files for "fragment_id"
		 * from rank 0 using MPI
		 */
		bool recvFragsViaMPI(int fragment_id);

		/**
		 * Worker - fill alias file with the names of database fragments
		 */		
		void WriteAliasFile( const std::string &alias_filename, const std::vector< int > &fragment_list);
		void DistributeDB();
		void ReadDBSpecs(const string& spec_file);
		void CreateGroupComm();

		/**
		 * Worker - Tell the writer the search has completed and how 
		 *          large the generated results are
		 * @param	dest	The writer process rank
		 * @param	frag_id	The fragment id that was searched
		 * @param	queries	A vector of query indexes that were searched
		 */
		void master();
		void super_master();
		void scheduler();
		void writer();

		/**
		 * A 'worker' process that actually searches the database
		 */
		void worker();

		/**
		 * Call this with the program arguments to start mpiBLAST
		 */
		int main( int argc, char* argv[] );
		/**
		 * parse input arguments 
		 */
		void ParseArguments(int argc, char* argv[]);
		int CalculateFirstAssignment();
		
		/** True if the user requested the database be removed */
		bool remove_db;	

	protected:
		MpiBlastConfig config;	/**< Class to read in configuration file */
		FragmentListFile frag_list_file;	/**< Class to read and modify local fragment list file */
		std::string exec_path;		/**< Directory that mpiblast was executed from */
		std::ofstream log_file, pro_phile;	/**< output file streams for logging and profiling information */
		std::string query_filename;	/**< Full path to query file */
		std::string output_filename;	/**< Full path to output file */
		std::string database_name;	/**< Name of database to query */
		std::string blast_cl;		/**< Base BLAST command line including user specified options */
		std::string blast_type;		/**< Type of blast search ( blastn, blastp, etc. ) */
		std::string db_type;			/**< n for nucleotide, p for protein */
		std::string query_file;		/**< Just the filename of the query file (pathname stripped) */
		std::string local_query_filename;	/**< path to a local copy of the query for workers */
		std::string filter_filename;	/**< Full path to gi filter file */
		std::string filter_file;		/**< Just the filename of the gi filter file (pathname stripped) */
		std::string local_filter_filename;	/**< path to a local copy of the gi filter for workers */
		std::string log_file_name;   /**< Log file name */
		std::string config_f_name;	/**< path to configuration file */
		bool disable_mpi_db;	/**< when true, workers don't send bioseqs to writer, instead
									     writer looks them up on disk */
		bool copy_via_set;      /**< whether copy method has been set */

		std::vector< std::string > scheduler_opts;	/**< Command line options for the master node to generate output */
		std::vector< std::string > writer_opts;	/**< Command line options for the master node to generate output */
		std::vector< std::string > worker_opts;	/**< Command line opts for each worker node to search correctly */
		
		int concurrent_accesses;	/**< number of concurrent accesses to shared storage */
		int query_batch;		/** number of queries assigned to the workers as a batch */
		int copy_via;			/**< how to copy the fragments */
		unsigned int total_frags;
		
		BLAST* ncbiblast;
		Scheduler* myscheduler;
		WriterMaster* writer_m;
		WriterWorker* writer_w;
		IrecvMgr* _irm;
		
		// master-specific data structures
		/** Classes to track progress of queries and distribution of database */
		// std::map< int, std::vector< std::pair< int, SeqAlignPtr > > > results_map;
		int output_queryI;	/**< Stores the index of the next query to output */
		// std::list< ResultData > result_data;	/**< A list of results that are in the process of being recv'd */
		
		int query_count;	/**< The number of queries in the query file */
		int output_format;	/**< The user-selected output format, defaults to -1 */

		// worker-specific data structures
		// std::map< int, std::pair< int, unsigned char* > > result_id_data_map; /**< maps result ID's to a pair of buffer size and result buffer */
		int result_id;

private:
		/**
		 * This internal function gets called by the other sendResults()
		 */
		// void sendResults( int dest, int queryI);

};

unsigned int determineLastQueryNum(const std::string outputfile);
void makeNewQuery(std::string& queryfile, std::string& query_filename, unsigned int lastquerynum);
unsigned int fragmentCount(const std::string sharedpath, const std::string dbname, const std::string dbtype);
#endif  // __mpiblast_hpp__
