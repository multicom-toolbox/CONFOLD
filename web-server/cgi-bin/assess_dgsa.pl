#!/usr/bin/perl -w
# Badri Adhikari, 1/17/2014

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use lib dirname(abs_path($0));
use confold;

my $dircns = shift;
die "dircns? Why! Why!" if not (defined $dircns && -d $dircns);
$dircns = abs_path($dircns);

print "[$0]: ".(localtime)."\n";
my @seq_files = <$dircns/*.fasta>;
die ":( fasta not found!" if not defined $seq_files[0];
my $seqId = basename($seq_files[0], ".fasta");
my $seq = seq_fasta("$dircns/${seqId}.fasta");
my @ss_files = <$dircns/*.ss>;
my $fileSS = $ss_files[0];
chdir $dircns or die $!;
my @pdbList = load_pdb($dircns);

my %listTbl = ();
$listTbl{"$dircns/contact.tbl"}  = 1 if (-f "$dircns/contact.tbl"  and count_lines("$dircns/contact.tbl") > 0);
$listTbl{"$dircns/ssnoe.tbl"}    = 1 if (-f "$dircns/ssnoe.tbl"    and count_lines("$dircns/ssnoe.tbl") > 0);
$listTbl{"$dircns/hbond.tbl"}    = 1 if (-f "$dircns/hbond.tbl"    and count_lines("$dircns/hbond.tbl") > 0);
$listTbl{"$dircns/dihedral.tbl"} = 1 if (-f "$dircns/dihedral.tbl" and count_lines("$dircns/dihedral.tbl") > 0);

print "\n";
printf "%9s : $dircns\n", "dir";                                   
printf "%9s : ${seqId}.fasta\n", "fasta";                             
printf "%9s : ".length($seq)."\n", "length";
foreach my $fileTbl (sort keys %listTbl){
	printf "%9s : ".count_lines($fileTbl)." lines\n", basename($fileTbl,".tbl");      
	# also verify if CNS actually accepted all the distance restraints provided
	next if $fileTbl =~ m/dihedral.tbl/;
	my $search_string = "N1";
	$search_string = "N2" if $fileTbl eq "$dircns/ssnoe.tbl";
	$search_string = "HBND" if $fileTbl eq "$dircns/hbond.tbl";
	confess ":( log file not found ($dircns/dg_sa.log)" if not -f "$dircns/dg_sa.log";
	my $result = `grep -m 1 -e NOEPRI.*$search_string $dircns/dg_sa.log`;
	my @C = split /\s+/, $result;
	my $count = $C[($#C)-1];
	confess ":( CNS did not accept all restraints of $fileTbl! Something wrong somewhere! Only $count accepted" if ($count != count_lines($fileTbl));
}

# remove "trial" structure of corresponding "accepted" structure because they are same
for(my $i = 0; $i <= 1000; $i++){
	next if not -f "${seqId}a_$i.pdb";
	print "\ndeleting ${seqId}_$i.pdb because ${seqId}a_$i.pdb exists!";
	system_cmd("rm $dircns/${seqId}_$i.pdb");
}

my %energyNoe = ();
foreach my $pdb (@pdbList) {
	next if $pdb =~ m/sub_embed/;
	next if $pdb =~ m/extended/;
	$energyNoe{$pdb} = get_cns_energy($pdb, "noe");
}
my $bestPdb = (sort {$energyNoe{$a} <=> $energyNoe{$b}} keys %energyNoe)[0];

print "\n";
print "$seq [input]\n";
print "".ss_line_ss($fileSS)." [input]\n" if (defined $fileSS && -f $fileSS);
foreach my $fileTbl (sort keys %listTbl){
	my $flagDihedral = 0;
	$flagDihedral = 1 if $fileTbl =~ m/dihedral.tbl/;
	print "".coverage_tbl("${seqId}.fasta", $fileTbl, $flagDihedral)."\n";
}
print "".dssp_ss_pdb($bestPdb)." [best prediction, ".basename($bestPdb, ".pdb")."]\n";

print "\n";
print "           ENERGY            CLASH     SS              NOE SATISFIED(±0.2A)            SUM OF DEVIATIONS >= 0.2     PDB\n";
print "--------------------------  -------  -------  ---------------------------------------  -------------------------  --------\n";
print "TOTAL  VDW    BOND   NOE    2.5 3.5  H   E    CONTACTS  SS-NOE    HBONDS    DIHEDRAL   CONTACTS SS-NOE   HBONDS\n"; 
foreach my $pdb (sort {$energyNoe{$b} <=> $energyNoe{$a}} keys %energyNoe){
	my ($e1, $e2, $e3, $e4);
	my ($c1, $c2);
	my ($h, $e) = ("-", "-");
	my ($n1, $n2, $n3, $n4);
	my ($s1, $s2, $s3);
	$e1 = $e2 = $e3 = $e4 = $c1 = $c2 = $h = $e = $n1 = $n2 = $n3 = $n4 = $s1 = $s2 = $s3 = "-";
	$e1 = get_cns_energy($pdb, "overall");
	$e2 = get_cns_energy($pdb, "vdw");
	$e3 = get_cns_energy($pdb, "bon");
	$e4 = $energyNoe{$pdb};
	$c1 = clash_count($pdb, 2.5);
	$c2 = clash_count($pdb, 3.5);
	$h  = count_ss_match($pdb, $fileSS, "H") if (defined $fileSS && -f $fileSS);
	$e  = count_ss_match($pdb, $fileSS, "E") if (defined $fileSS && -f $fileSS);
	$n1 = count_satisfied_tbl_rows($pdb, "$dircns/contact.tbl", "noe") if defined $listTbl{"$dircns/contact.tbl"};
	$n2 = count_satisfied_tbl_rows($pdb, "$dircns/ssnoe.tbl", "noe")   if defined $listTbl{"$dircns/ssnoe.tbl"};
	$n3 = count_satisfied_tbl_rows($pdb, "$dircns/hbond.tbl", "noe")   if defined $listTbl{"$dircns/hbond.tbl"};
	$n4 = count_satisfied_tbl_rows($pdb, "$dircns/dihedral.tbl", "dihedral")   if defined $listTbl{"$dircns/dihedral.tbl"};
	$s1 = sum_noe_dev($pdb, "$dircns/contact.tbl") if defined $listTbl{"$dircns/contact.tbl"};
	$s2 = sum_noe_dev($pdb, "$dircns/ssnoe.tbl")   if defined $listTbl{"$dircns/ssnoe.tbl"};
	$s3 = sum_noe_dev($pdb, "$dircns/hbond.tbl")   if defined $listTbl{"$dircns/hbond.tbl"};
	printf "%-6s %-6s %-6s %-6s %-3s %-3s  %-3s %-3s  ",  $e1, $e2, $e3, $e4, $c1, $c2, $h, $e;
	printf "%-9s %-9s %-9s %-9s  %-8s %-8s %-8s %-25s\n", $n1, $n2, $n3, $n4, $s1, $s2, $s3, basename($pdb, ".pdb"); 
}

print "\n";
foreach my $pdb (sort {$energyNoe{$b} <=> $energyNoe{$a}} keys %energyNoe){
	my $ss = dssp_ss_pdb($pdb);
	$ss =~ s/C/-/g;
	printf $ss." [".basename($pdb, ".pdb")."]\n";
}

foreach my $fileTbl (sort keys %listTbl){
	next if $fileTbl =~ m/dihedral/;
	print "\n";
	foreach my $pdb (sort {$energyNoe{$b} <=> $energyNoe{$a}} keys %energyNoe){
		print "".noe_tbl_violation_coverage($pdb, $fileTbl)." [ violation of ".basename($fileTbl)." in ".basename($pdb, ".pdb")." ]\n";
	}
}

print "\n";
my $i = 1;
foreach( sort {$energyNoe{$a} <=> $energyNoe{$b}} keys %energyNoe){
	print "model$i.pdb <= $_\n";
	system_cmd("mv $_ ${seqId}_model$i.pdb");
	$i++;
	last if $i > 5;
}
system_cmd("rm $dircns/dg_sa.log");

print "[$0]: ".(localtime)."\n";
print "################################################################################\n";
