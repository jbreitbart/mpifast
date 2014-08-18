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

/*
 * mpiformatdb - a program to format BLAST databases into several fragments 
 * for mpiblast.
 * 
 * Overall program flow is as follows:
 * 1)  parse and sanity check command line options
 * 2)  read input sequence database in FastA format
 *     - index the location and length of each sequence entry 
 *       in the input file
 * 3)  write a temp file with the sequences reordered to be 
 *     largest first
 * 4)  call the hacked version of formatdb's Main() function that 
 *     always creates the desired number of database fragments.
 *     Since the hacked formatdb always adds the next sequence to
 *     the smallest fragment, fragment sizes are guaranteed to be
 *     equal.
 * 5)  create and update supporting files like the .mbf, .dbs, and
 *     the .{n,p}al alias file
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// #include "mpi.h"

#include "mpiblast_config.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_types.h"
#include "fragment_list.hpp"
#include "file_util.hpp"

#include <signal.h>
#include <map>
#include <cstdio>
#include <sstream>
using namespace std;

extern "C"{
#include "ncbi.h"
//#include "ncbiwin.h"
}

#include "getopt.h"

static string tmp_reorder_fname;
static bool print_status = true;	/**< Set this to false to not print %complete counter */
static bool use_basename = false;

/* for creating db specification file */
Int8 global_db_len = 0;
Int4 global_dbseq_num = 0;
bool create_dbs = true;

void print_usage( const char* pname ){
	cout << "Usage: " << pname << " <--nfrags=<NumberOfFrags>|--frag-size=<SizeOfFrags> (in MB)> <formatdb arguments>\n";
	cout << "Options:\n";
	cout << "\t-f  --config-file=<CONFIGFILE>  Specify the location of the mpiBLAST configuration file\n";
	cout << "\t    --nfrags=<NumberOfFrags>    Specifies the number of fragments to split the database into\n";
	cout << "\t    --frag-size=<SizeOfFrags>   Specifies the size for database fragments, can not be used with --nfrags\n";
	cout << "\t    --version                   Prints the version number\n";
	cout << "\t    --debug                     Prints out debugging information\n";
	cout << "\t    --enable-reorder            Enable the sequences reordering..\n";
	cout << "\t    --quiet                     Don't show % complete counter.\n";
	cout << "\t    --update_db=<DBtoUpdate>    Update database by adding additional fragments\n";
	cout << "\nExample:\n";
	cout << "mpiformatdb -i ecoli.nt --nfrags=5\n\n";
}

// silly Windows defines _snprintf instead of snprintf
// #if defined(WIN32) && !defined(snprintf)
// #define snprintf _snprintf
// #endif // WIN32

#define BUFFER_SIZE 10000

/**
 * Records location of a FastA sequence within a file
 */
typedef struct FastAEntry {
	FILE* file;	/**< the file containing this sequence */
	size_t start;	/**< file offset of the first byte of the header */
	size_t end;	/**< file offset of the first byte of the next header - 1 */
} FastAEntry_t;

/**
 * count the number of characters in each sequence entry of seq_file (assuming it is FastA format)
 * @param seq_file	A C File object pointing to the open sequence file
 * @param seq_lengths	A map that will be filled with the lengths of each sequence entry 
 *			and the corresponding file offset of each entry
 * @param nucleotide	This variable is set to true if > 75% of the sequence looks like nucleotides
 * 			e.g. A, C, G, T
 * @return	The total number of sequence characters in the file
 */
uint64 indexSeqLengths( string& input_filename , multimap< uint64, FastAEntry_t >& seq_lengths, bool& nucleotide ){

	// open the FastA sequence file
 	FILE *seq_file;
	if (!(seq_file = fopen( input_filename.c_str(), "r" ))) {
 		cerr << "Error: Unable to open " << input_filename << endl;
		exit(-1);
 	}

	int i;
	FastAEntry_t fat;
	fat.file = seq_file;
	fat.start = 0;
	fat.end = 0;
	seq_lengths.clear();
	uint64 total_chars = 0, counter = 0 ;
	uint64 file_size = statFileSize(input_filename.c_str());
	uint64 nucleotide_count = 0;
	uint64 cur_seq_chars = 0;
	uint64 cur_file_char = 0;
	
	char buffer[BUFFER_SIZE];
	char defline = 0; //defline is used as bool
	nucleotide = false;
	cout << "Reading input file\n";
 	while ((fgets(buffer, BUFFER_SIZE, seq_file)) != 0) {
		for (i=0; buffer[i] != '\0'; i++) {
			cur_file_char++;
			if (defline == 0 && buffer[i] == '>') {
				defline = 1;
				
				global_dbseq_num++;
				
				if( total_chars > 0 ){
					fat.end = cur_file_char - 1;
					seq_lengths.insert( multimap< uint64, FastAEntry_t >::value_type( cur_seq_chars, fat ) );
				}
				cur_seq_chars = 0;
				fat.start = fat.end;
				
			} else if ((buffer[i] == '\n') || (buffer[i] == '\r')) {
				defline = 0;
			} else if (defline == 0) {
				if( buffer[i] == 'a' || buffer[i] == 'A' ||
					buffer[i] == 'c' || buffer[i] == 'C' ||
					buffer[i] == 'g' || buffer[i] == 'G' ||
					buffer[i] == 't' || buffer[i] == 'T' ){
					nucleotide_count++;
				}
				cur_seq_chars++;
				total_chars++;
			}
		}
		if ( print_status && (((cur_file_char % 100000) == 0)  || cur_file_char == file_size) )
			printProgress(cur_file_char, file_size);

		counter++ ;
	}	
	// add the final length/offset combination
	
	fat.end = cur_file_char;
	seq_lengths.insert( multimap< uint64, FastAEntry_t >::value_type( cur_seq_chars, fat ) );
	
	// do nucleotide calculation
	uint64 tmp_chars = total_chars;
	tmp_chars *= 3;
	tmp_chars /= 4;
	nucleotide = nucleotide_count > tmp_chars;
	
	global_db_len = total_chars;
	
	cout << "Done, read " << counter << " lines\n" ;	
	
	return total_chars;
}

/**
 * Rewrite the input database entries to be in descending order based on entry length
 */
void reorderInputSequences( FILE* out_file,  multimap< uint64, FastAEntry_t >& seq_lengths ){
	multimap< uint64, FastAEntry_t >::iterator seqlen_iter = seq_lengths.end();
	float64 total_seqs = seq_lengths.size() ;
	uint64 output_progress = 0;
	// write out each sequence
	cout << "Reordering " << seq_lengths.size() << " sequence entries\n";
	while( seqlen_iter != seq_lengths.begin() ){
		seqlen_iter--;
		fseek( seqlen_iter->second.file, seqlen_iter->second.start, SEEK_SET );
		char buffer[BUFFER_SIZE];
		char defline = 0; //defline is used as bool
		size_t rdata;
		size_t data_len = seqlen_iter->second.end - seqlen_iter->second.start; 
		size_t rsize = BUFFER_SIZE < data_len ? BUFFER_SIZE : data_len;
	 	while ((rdata = fread(buffer, 1, rsize, seqlen_iter->second.file)) != 0) {
			data_len -= rdata;
			rsize = BUFFER_SIZE < data_len ? BUFFER_SIZE : data_len;
			fwrite( buffer, 1, rdata, out_file );
		}
		fprintf( out_file, "\n" );
		output_progress++ ;
		if ( print_status && (((output_progress % 10) == 0)  || output_progress == total_seqs)) 
			printProgress(output_progress, total_seqs);
	}
	fflush( out_file );
}


/**
 * add db frags to existing alias file for update_db
 */
void modifyAliasFile( const string& alias_orig, int high_frag_old, int high_frag_new ) {
	string alias_temp = alias_orig + "XXXXXX";
	getTempFileName( alias_temp );
	ifstream alias_in( alias_orig.c_str() );
	ofstream alias_out( alias_temp.c_str() );
	string line;
	string frag_addition = alias_orig.substr( 0, alias_orig.length() - 3 );
	char fragnum[4];
	try{
	while( getline( alias_in, line ) ){
		alias_out << line;
		// look for DBLIST token
		if( line.size() > 6 && line.substr(0,7) == "DBLIST " ) {
			for( int n = high_frag_old; n < high_frag_new; n++ ){
				// we always use 3 digit fragment identifiers
				sprintf(fragnum,"%03d",n);
				alias_out << " " << frag_addition << fragnum;
			}
		}
		alias_out << endl;
	}
	rename( alias_temp.c_str(), alias_orig.c_str() );
	}catch( exception& e ){
		cerr << e.what() << endl;
	}
	return;
}

void terminateProgram(int sig){
	cerr << "mpiformatdb bailing out with signal " << sig << endl;
	cerr << "removing " << tmp_reorder_fname << endl;
	removeFile( tmp_reorder_fname );
	exit(sig);
}
/**
 * Write an alias file (included to work around a windows-specific toolbox bug)
 * This is a rough copy of MpiBlast::WriteAliasFile
 */
void WriteAliasFile( const MpiBlastConfig& config, const string& database_name, 
					const string& alias_filename, const vector< int >& fragment_list)
{
	if(fragment_list.size() == 0){
		return;	// don't write a file
		// this isn't an error condition when # workers > # db frags
//		throw __FILE__ "(WriteAliasFile) Empty fragment_list";
	}
	
	ofstream alias_file( alias_filename.c_str() );

	if( !alias_file.is_open() ){
		cerr << "Error opening " << alias_filename << endl;
		throw __FILE__ "(MpiBlast::worker): Unable to open alias file";
	}

	alias_file << "TITLE " << config.localPath() << database_name << endl;
	alias_file << "DBLIST";
	
	for( uint iter = 0; iter != fragment_list.size(); iter++ ){
		alias_file << " " << database_name << ".";
		char fragid[8];
		
		memset(fragid, 0, 8);

		// always use 3 digit fragment identifiers
		sprintf( fragid, "%03d", fragment_list[ iter ] );

		alias_file << fragid;
		if( debug_msg ){
			LOG_MSG << "node " << rank << " fragid: " << fragid << endl;
		}
	}
	alias_file << endl;
	alias_file.close();
		
}

	
int main( int argc, char* argv[] ){
	MpiBlastConfig config;	/**< configuration file parser */
	string localPath, blast_cl, blast_type, input_filename;
	multimap< uint64, FastAEntry_t > seq_lengths;	/**< sequence entry file offsets and lengths */
	
	log_stream = &cout;

	if( argc < 2 ){
		if( argc == 0 )
			print_usage( "mpiformatdb" );
		else
			print_usage( argv[0] );
	}
	
	signal( SIGINT, terminateProgram );
	signal( SIGTERM, terminateProgram );
	signal( SIGSEGV, terminateProgram );


	//
	// must read -i and -p options
	//
	const char* options = "-i:p:n:N:o:t:";
	int ac = argc;
	char** av = argv;
	int opt_code, opt_type = 0;
	opterr = 0;
	int config_opt = 0;
	int long_index = 0;
	vector< string > formatdb_opts;
	formatdb_opts.push_back( "formatdb" );
	
	string config_f_name = ".ncbirc";
	bool print_version = false;
	int nfrags = -1;
	int frag_size = -1;
	bool nucleotide = false;	/**< this gets set to true later if the db appears to contain DNA */
	bool reorder = false;	/**< Set this to false to skip the reordering by size step */
	bool incompatible_option = false;	/**< Set to true if a command line option incomatible with db segmentation is given */
	string base_db_name;	/**< Stores the user specified -n base db name option */
	bool user_indexing_override = false;	/**< Set to true when the user explicitly gives -o F on the command line */
	// the full path name of the db to update
	string update_db;
	string db_title;
	
	struct option long_opts[] = {
		{"config-file", required_argument, &config_opt, 1 },	// this option is obsolete and only remains
									// to notify users of obsolescence
		{"nfrags", required_argument, &config_opt, 2 },
		{"frag-size", required_argument, &config_opt, 3 },
		{"version", no_argument, &config_opt, 6 },
		{"debug", no_argument, &config_opt, 8 },
		{"update_db", required_argument, &config_opt, 9 },
		{"enable-reorder", no_argument, &config_opt, 10 },
		{"quiet", no_argument, &config_opt, 11 },

		{0,0,0,0}	// signifies termination of option list
	};

	// get the mpiformatdb and interesting formatdb options with getopt_long
	while((opt_code = getopt_long( ac, av, options, long_opts, &long_index ) ) != EOF ){
		switch( opt_code ){
			case 0:
				if( config_opt == 1 ){
//					cerr << "The mpiBLAST config file is no longer used.  Instead, ";
//					cerr << "configuration options should be specified directly in the ";
//					cerr << ".ncbirc configuration file\n";
					config_f_name = optarg;
				}
				if( config_opt == 2 ){
					nfrags = atoi( optarg );
				}
				if( config_opt == 3 ){
					frag_size = atoi( optarg );
				}
				if( config_opt == 6 ){
					print_version = true;
				}
				
				if( config_opt == 8 ){
					debug_msg = true;
				}
	  			if( config_opt == 9 )
	   				update_db = optarg;
	  			if( config_opt == 10 )
	   				reorder = true;
	  			if( config_opt == 11 )
	   				print_status = false;
				break;
/*  mpiBLAST.formatdb no longer recognizes these options to avoid argument overlap with NCBI
			case 'f':
				config_f_name = optarg;
				break;
			case 'S':
				frag_size = atoi( optarg );
				break;
*/
			case 'o':
				if( strcmp( optarg, "F" ) == 0 ){
					addOpt( formatdb_opts, 'o', optarg );
//					cerr << "Warning, mpiBLAST requires database indices.  ";
//					cerr << "Allowing -o F anyways...\n";
					user_indexing_override = true;
				}
				break;	// we force -o T anyways
			case 'N':
				nfrags = atoi( optarg );
				break;
			case 'i':
				input_filename = optarg;
				break;
			case 'p':
				addOpt( formatdb_opts, 'p', optarg );
				if( strcmp( optarg, "F" ) == 0 )
					blast_type = "n";
				else
					blast_type = "p";
				break;
			case 'n':
				base_db_name = optarg;
				use_basename = true;
				break;	// override the user specified -n
			// scan for options incompatible with db formatting
			case 'a':
			case 'b':
			case 'F':
			case 'B':
				addOpt( formatdb_opts, opt_code, optarg );
				incompatible_option = true;
				break;
			case 't':
				db_title = optarg;
				addOpt( formatdb_opts, opt_code, optarg );
				break;
			case 1:
				formatdb_opts.push_back( optarg );
				break;
			case '?':
				addOpt( formatdb_opts, optopt, optarg );
				opt_type = optopt;
				break;
			default:
				addOpt( formatdb_opts, opt_code, optarg );
//				LOG_MSG << "Default " << (char)optopt << endl;
				opt_type = optopt;
				break;
		}
	}
	if( print_version ){
		cout << PACKAGE << " version " << VERSION << endl;
		cout << "Built " << __DATE__ << " at " << __TIME__ << endl;
		return 0;
	}

	// just in case an argument was left behind -- go get it
	for( int optI = optind; optI < ac; optI++ ){
		formatdb_opts.push_back( av[ optI ] );
	}
	
	//If we didn't specify either of -N or -S just run standard formatdb
	if( nfrags <= 0 && frag_size <= 0 ){
		if( input_filename != "" )
			addOpt( formatdb_opts, 'i', input_filename.c_str() );
		if( base_db_name != "" )
			addOpt( formatdb_opts, 'n', base_db_name.c_str() );
		string formatdb_cl;
		for( uint optI = 0; optI < formatdb_opts.size(); optI++ )
			formatdb_cl += formatdb_opts[ optI ] + " ";
			
		LOG_MSG << "Executing: " << formatdb_cl << endl;
		// execute formatdb
		initNCBI( formatdb_opts );
	
		int retval = Main();

		cleanupNCBI();
		return retval;
	}
	
	if( input_filename != "" && incompatible_option ){
		cerr << "Error!  The options -a -b -F and -B can't be used when ";
		cerr << "attempting to format a database with mpiformatdb\n";
		return -1;
	}
	if( input_filename == "stdin" && reorder ){
		cerr << "Error!  -i stdin is incompatible with input sequence reordering.\n";
		cerr << "Add --skip-reorder to disable sequence reordering\n";
		return -1;
	}
	if( input_filename == "stdin" && db_title == "" ){
		cerr << "Error!  -i stdin requires that a DB title is specified using -t <title>.\n";
		return -1;
	}
	if( input_filename == "stdin" && blast_type == "" ){
		cerr << "Error!  -i stdin requires that the sequence type is specified using -p <T|F>.\n";
		return -1;
	}

	// read the configuration file
	try{
		config = MpiBlastConfig( config_f_name );
	}catch(const char *error){
		cerr << "Fatal Error:" << endl;
		cerr << error << endl;
		return -1;
	}
	localPath = config.localPath();
	
	// check for valid option combinations
	if( input_filename == "" ){
		print_usage( argv[0] );
		cerr << "You must specify an input file in FastA format\n";
		return -1;
	}

	//If update_db, check to make sure that we have a pre-formatted db
	string formatted_db = getFilename( input_filename );
	if( input_filename == "stdin" )
		formatted_db = db_title;
	int high_frag_old = 0;
	
	
	// make sure the output directory is writable
	if( !checkDirWritePerms( config.sharedPath() ) ){
		cerr << "Error: unable to write to shared storage path: " << config.sharedPath() << endl;
		return -1;
	}
	
	// read in the input file to ascertain the number and lengths of 
	// sequence entries it contains
	uint64 bases = 0;
	if( input_filename != "stdin" )
		bases = indexSeqLengths( input_filename, seq_lengths, nucleotide );
	bases /= 1000000;	// roughly convert to MB 
	
	// if we get here, we will be fragmenting the database.
	// check that the user specified a fragmentation	
	if( nfrags > 0 && frag_size > 0 ){
		print_usage( argv[0] );
		cerr << "You may only specify one of --nfrags or --frag-size\n";
		return -1;
	}
	
	// if the number of frags wasn't specified, figure it out using frag size
	if( nfrags <= 0 ){
		nfrags = bases / frag_size;
		if( nfrags == 0 || bases % frag_size != 0 )
			nfrags++;
	}	
		
	if( input_filename != "stdin" && nfrags > seq_lengths.size() ){
		cerr << "WARNING: The database contains only " << seq_lengths.size() << " sequence entries ";
		cerr << "and can't be broken into " << nfrags << " fragments.  Some fragments will ";
		cerr << "remain empty\n";
		nfrags = seq_lengths.size();
	}

	string db_input_fname = input_filename;	/**< contains the db file name ultimately given to formatdb */
	if( reorder ){
		// write out a file of sequence entries re-ordered on length
 		FILE *tmp_reorder_fp;
		char* tmp_path;

		// check $TMPDIR for a temp directory
		if( (tmp_path = getenv( "TMPDIR" )) != NULL ){
			tmp_reorder_fname = tmp_path + PATH_SEPARATOR + "reorderXXXXXX";
		}else{
#ifdef WIN32
			tmp_reorder_fname = input_filename + ".XXXXXX";
#else
			tmp_reorder_fname = "/tmp/reorderXXXXXX";
#endif
		}
		getTempFileName( tmp_reorder_fname );
		if (!(tmp_reorder_fp = fopen( tmp_reorder_fname.c_str(), "w" ))) {
 			cerr << "Error: Unable to open " << tmp_reorder_fname << " for writing\n";
			return -1;
 		}
		reorderInputSequences( tmp_reorder_fp, seq_lengths );
		fclose( tmp_reorder_fp );	// done with the file, close it
		db_input_fname = tmp_reorder_fname;
	}
	// add the input db file to the formatdb command line
	addOpt( formatdb_opts, 'i', db_input_fname.c_str() );
	
	// check whether the user specified a database type
	if( blast_type == "" ){
		if( nucleotide ){
			blast_type = "n";
			cerr << "Database type unspecified, assuming nucleotide\n";
			addOpt( formatdb_opts, 'p', "F" );
		}else{
			blast_type = "p";
			cerr << "Database type unspecified, assuming protein\n";
			addOpt( formatdb_opts, 'p', "T" );
		}
	}else if( input_filename != "stdin" ){
		if( blast_type == "p" && nucleotide )
			cerr << "WARNING: Protein database type specified, but the database appears to contain nucleotides!\n";
		else if( blast_type == "n" && !nucleotide )
			cerr << "WARNING: Nucleotide database type specified, but the database appears to contain amino acids!\n";
	}
	
/*
 * Here we do checks for updating a database, then read in the necessary info regarding the old one.
 * Next we make sure we're working with the same type of db (nucleotide or protein)
 */
	if ( update_db.size() > 0 ) {
		string al_orig = update_db + "." + blast_type + "al";
		if ( !doesFileExist(al_orig) ){
			cerr << al_orig << " does not exist.  Make sure the update and existing database are both the same type (NT or AA)\n";
			return -1;
		}
	}

	// specify the number of fragments using the -N option to the mpiblast_formatdb command line
	formatdb_opts.push_back( "-N" );
	ostringstream oss;
	oss << nfrags;
	formatdb_opts.push_back( oss.str().c_str() );
	
	// specify the output location for the database with -n
	string base_name;
	if(use_basename) {
	    base_name = config.sharedPath() + formatted_db;
		formatdb_opts.push_back( "-n" );
		formatdb_opts.push_back( base_name );
	} else {
        base_name = input_filename;
    }
	
	// force index construction with -o T
	if( !user_indexing_override ){
		formatdb_opts.push_back( "-o" );
		formatdb_opts.push_back( "T" );
	}

	cout << "Breaking " << formatted_db << " into " << nfrags << 
		" fragments\n";
	
	string formatdb_cl;
	for( uint optI = 0; optI < formatdb_opts.size(); optI++ )
		formatdb_cl += formatdb_opts[ optI ] + " ";
	LOG_MSG << "Executing: " << formatdb_cl << endl;
	
	//
	// remove previously created databases with the same name
	//
	const vector< string >& extensions = FragmentExtensions( blast_type );
	if ( !update_db.size()) {
		int fragI = 0;
		while(true){
			char number[4];
			// we always use 3 digit fragment ids
			snprintf( number, 4, "%03d", fragI );
			uint extI = 0;
			for( ; extI < extensions.size(); extI++ ){
				string old_frag = base_name + "." + number + extensions[extI];
				if( extensions[extI].find("in") != string::npos && !doesFileExist( old_frag ) )
					break;
				if( !doesFileExist( old_frag ) )
					continue;
				removeFile( old_frag, false );
			}
			// break out if the .[p|n]in file didn't exist
			if( extI < extensions.size() )
				break;
			// move on to the next fragment
			fragI++;
		}
	}

	// execute formatdb
	initNCBI( formatdb_opts );
	
	int retval = Main();

	cleanupNCBI();

	// delete the temporary file used for sequence reordering
	if( reorder ){
		if( removeFile( db_input_fname, true ) ){
			cerr << "Failed to delete temp file: " << db_input_fname << endl;
			cerr << "Please delete it manually.\n";
		}
	}

	if( retval != 0 ){
		cerr << "There was an error executing formatdb.  Check formatdb.log\n";
		
		// Only return error values <= 0 (since this program now returns the number of
		// fragments actually created).
		return retval > 0 ? -retval : retval;
	}

	// if formatdb exited cleanly then the requested number of fragments were created
	cout << "Created " << nfrags << " fragments.\n";
	
	// ensure that files with all 7 fragment extensions were created
	for( int fragI = 0; fragI < nfrags; fragI++ ){
		char number[4];
		snprintf( number, 4, "%03d", fragI );
		for( uint extI = 0; extI < extensions.size(); extI++ ){
			string frag_fname = base_name + "." + number + extensions[extI];
			ifstream frag_fstream( frag_fname.c_str() );
			if( !frag_fstream.is_open() )
				ofstream frag_creator( frag_fname.c_str() );
		}
	}

	// if there's only a single fragment make an alias file for it anyways
	string al_dest = config.sharedPath() + formatted_db + "." + blast_type + "al";
	if( nfrags == 1 ){
		ofstream alias_out( al_dest.c_str() );
		alias_out << "DBLIST " << base_name << ".000" << endl;
		alias_out.close();
	}

#ifdef WIN32
	// 11/2004:
	// win32 ncbi toolbox has a bug which causes the db path to be concatenated
	// twice into the fragment file path when reading the database.  to workaround,
	// create an alias file without the full path given in each fragment reference
	vector< int > frag_int_list;
	for( int fragI = 0; fragI < nfrags; fragI++ )
		frag_int_list.push_back( fragI );
	WriteAliasFile( config, formatted_db, al_dest, frag_int_list );
#endif

	// If a database is being updated rename the new fragments and add them
	// to the existing alias file.
	if (update_db.size() > 0) {
		char new_number[4];
		char old_number[4];

		for (int n = 0; n < nfrags ; n++){
			// we always use 3 digit fragment ids
			snprintf( new_number, 4, "%03d", (n+high_frag_old) );
			snprintf( old_number, 4, "%03d", n );
			for( uint extI = 0; extI < extensions.size(); extI++ ){
				string old_frag = base_name + "." + old_number + extensions[extI];
				string new_frag = update_db + "." + new_number + extensions[extI];
				rename( old_frag.c_str(), new_frag.c_str() );
			}
		}
		//Remove created [.nal|.pal] file and modify original one
		remove( al_dest.c_str() );
		string al_orig = update_db + "." + blast_type + "al";
		modifyAliasFile( al_orig, high_frag_old, (nfrags + high_frag_old) );

	}
	
	// Create .mbf file in SharedStorage for --copy-via=none
	string mbf_file;
    if(use_basename) {
        mbf_file = config.sharedPath() + formatted_db + ".mbf";
    } else {
        mbf_file = base_name + ".mbf";
    }
	if (update_db.size() > 0)  //If we're updating, use the old db name. else remove the file
		mbf_file = update_db + ".mbf";
	else 
		remove(mbf_file.c_str());
	//Now open up the output stream and write to it. high_frag_old == 0 if ! updating
	ofstream mbf_out(mbf_file.c_str(),ios::app);
	char fragid[5];
	for (int n = high_frag_old; n < (nfrags+high_frag_old); n++){
		snprintf(fragid,5, "%03d\n",n);
		mbf_out << fragid ;
	}
	mbf_out.close();
	
	if(create_dbs) {
		string dbs_file;
        if(use_basename) {
            dbs_file = config.sharedPath() + formatted_db + ".dbs";
        } else {
            dbs_file = base_name + ".dbs";
        }
		ofstream dbs_out(dbs_file.c_str());
		if(!dbs_out.is_open()) {
			cerr << "Error opening database specification file!" << endl;
			return -1;
		}
		dbs_out << global_db_len << endl;
		dbs_out << global_dbseq_num << endl;
		dbs_out.close();
	}

	cout << "<<< Please make sure the formatted database fragments are placed in " << config.sharedPath() << " before executing mpiblast. >>> \n" << endl;
			
	return 0;
}

