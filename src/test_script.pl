#!/usr/bin/perl -w
#
# Script to test mpiBLAST execution under various conditions
# (c) 2k3 aaron darling
#
use strict;

# configure these parameters according to your setup
my $nt_database = "yeast.nt";
my $aa_database = "nr";
my $nt_query = "~/development/mpiblast2/src/ech_query.fas";
my $aa_query = "~/development/testing/mpiblast/test3.aa";
my $blast_executable = "mpirun -np 6 ./mpiblast";
#my $blast_executable = "/home/koadman/software/ncbi/build/blastall";
#my $blast_executable = "mpirun -np 6 /mnt/isis/Software/mpiblast/mpiblast";
my $output_dir = "test_results/";
#my $output_dir = "ncbi_test_results/";
my $output_types = 11;

# actually do the work
# run [t]blast[npx] tests for each output type
my $test_cl;
$test_cl = "$blast_executable -p blastn -i $nt_query -d $nt_database -o $output_dir/blastn_test.typ";
run_output_tests();
$test_cl = "$blast_executable -p blastx -i $nt_query -d $aa_database -o $output_dir/blastx_test.typ";
run_output_tests();
$test_cl = "$blast_executable -p tblastn -i $aa_query -d $nt_database -o $output_dir/tblastn_test.typ";
run_output_tests();
$test_cl = "$blast_executable -p blastp -i $aa_query -d $aa_database -o $output_dir/blastp_test.typ";
run_output_tests();
$test_cl = "$blast_executable -p tblastx -i $nt_query -d $nt_database -o $output_dir/tblastx_test.typ";
run_output_tests();

sub run_output_tests {
	for( my $typeI = -1; $typeI <= $output_types; $typeI++ ){
		my $output_type_cl = $test_cl;
		$output_type_cl .= "$typeI -m $typeI" if $typeI >= 0;
		$output_type_cl .= " -J " if $typeI >= 10;
		
		# tack on debugging output
#		$output_type_cl .= " --debug 2> last_run.debug";
		
		print $output_type_cl."\n";
		`$output_type_cl`;
	}

}
