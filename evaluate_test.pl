#!/usr/bin/perl -w
# Badri Adhikari, 8/5/2014

# Build models for all 7 test cases. 
# Please refer the 'readme.txt' attached for more information about the 7 test cases.

use strict;
use warnings;
use Carp;
use File::Basename;

my $dir_in  = "./test/input";
my $dir_out = "./test/output";
my $tmscore = "/home/bap54/bin/TMscore";

my %test_case = qw/1eaz 1 1g7r 1 1guu 1 1qjp 1 1smx 1 1vjk 1 helix 1/;
print "PDB\tTM-score\tRMSD\tMODEL\n";
foreach my $pdb (sort keys %test_case){
	my @pdb_list = load_pdb("$dir_out/$pdb/stage1");
	@pdb_list = load_pdb("$dir_out/$pdb/stage2") if -d "$dir_out/$pdb/stage2";
	my $native = "$dir_in/$pdb.pdb";
	my $best_pdb;
	my %best_tm = ();
	$best_tm{"tm-score"} = 0;
	foreach my $p (@pdb_list){
		next if $p =~ /sub_embed/;
		my %tm = evaluate_using_TMscore($p, $native);
		next if $tm{"tm-score"} < $best_tm{"tm-score"};
		$best_pdb = $p;
		%best_tm = %tm;
	}
	confess "\nERROR! could not find pdb for $pdb!" if not defined $best_pdb;
	printf $pdb."\t%.2f\t%.2f\t$best_pdb\n", $best_tm{"tm-score"}, $best_tm{"rmsd"};
}

sub load_pdb{
	my $dir_chains = shift;
	confess "\n:( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdb_list = <$dir_chains/*.pdb>;
	if(not (@pdb_list)){
		@pdb_list = <$dir_chains/*.ent>;
	}
	confess "\nERROR! Directory $dir_chains has no pdb files!\n" unless(@pdb_list);
	return @pdb_list;
}

sub evaluate_using_TMscore{
	my $predicted = shift;
	my $native = shift;
	confess "\nERROR! Predicted pdb $predicted does not exit!" if not -f $predicted;
	confess "\nERROR! Native pdb $native does not exist!" if not -f $native;
	my @results = `$tmscore $predicted $native | grep -e RMSD\\ of -e TM-score\\ \\ \\  -e MaxSub-score -e GDT-TS-score -e GDT-HA-score`;
	if (not defined $results[0]){
		print "Executing: [$tmscore $predicted $native]\n";
		system("$tmscore $predicted $native");
		confess "\nError! TM-score failed!";
	}
	$results[0] =~ s/RMSD of  the common residues//g;
	$results[1] =~ s/TM-score//g;
	$results[2] =~ s/MaxSub-score//g;
	$results[3] =~ s/GDT-TS-score//g;
	$results[4] =~ s/GDT-HA-score//g;
	for(my $i = 0; $i <= 4; $i++ ){
        $results[$i] =~ s/\s+//g;
        $results[$i] =~ s/=//;
	}
	my %result = ();
	$result{"rmsd"}		= substr $results[0], 0, 5;
	$result{"tm-score"}	= substr $results[1], 0, 6;
	$result{"max-sub"}	= substr $results[2], 0, 6;
	$result{"gdt-ts"}	= substr $results[3], 0, 6;
	$result{"gdt-ha"}	= substr $results[4], 0, 6;
	return %result;
}

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "\nEXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "\nERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}
