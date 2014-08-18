#!/usr/bin/perl -w
#parse clog_print output
use strict;
no strict "refs";
my @nodes;
my %nodehash;
my $debug;
my %cats;
if ($ARGV[0] =~ m/-d/) {$debug=shift;}

foreach my $file (@ARGV){
#read to collect timming data
	open(CLOG,$file);
	while(<CLOG>){	

#map numbers to descriptions
		if (/start=(\d+).*desc=(\w{3}\w+)/){
			$cats{$1} = $2;
		}

		next unless (/^ts=(\d+\.\d+) type=raw.* len=(\d+), pid=(\d+) id=(\d+) data=(\d+) srcid=(\d+) desc=(\w*).*\n$/);
		my ($time , $len, $node, $id, $data, $srcid) = ($1,$2,$3,$4,$5,$6);
		push(@{$nodehash{$node}}, [$id,$time]);
	}
	foreach my $node (sort {$a <=> $b} keys %nodehash){
		my %times;
		foreach my $data (@{$nodehash{$node}}){
			my @dat = @{$data};
			push(@{$times{$dat[0]}},$dat[1]);
		}
		my $count=-1;
		foreach my $id (sort {$a <=> $b} keys %times){
			$count++;
			next unless ($count % 2 == 0);
			my $endid = $id+1;
			my $totaltime=0;
			for(my $i=0;$i<@{$times{$id}};$i++){
				$totaltime+=(${$times{$endid}}[$i]-${$times{$id}}[$i]);
				print "$id - $endid , # $i = ",${$times{$endid}}[$i]-${$times{$id}}[$i],"\n" if ($debug);
			}
			foreach (keys %cats) {
				if($_ eq $id ){
					$id = $cats{$_};
				}
			}
			printf ("%2d %15s %12f %2d\n",$node,$id,$totaltime,($#{$times{$id}}+1));
		}
	}
}

