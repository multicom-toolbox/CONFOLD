#!/usr/bin/perl -w
# Badri Adhikari, 1/23/2014

package confold;

use strict;
use warnings;
use Carp;
use Cwd;
use File::Basename;

BEGIN {
	require Exporter;
	our $VERSION = 1.0;
	our @ISA     = qw(Exporter);
	our @EXPORT  = qw(%AA3TO1 %AA1TO3 @AA3 @AA1);
	push @EXPORT, 'all_xyz_pdb';
	push @EXPORT, 'ca_xyz_pdb';
	push @EXPORT, 'calc_dist';
	push @EXPORT, 'carr2tbl';
	push @EXPORT, 'cb_xyz_pdb';
	push @EXPORT, 'cbrr2cb_dist';
	push @EXPORT, 'cbrr2r1a1r2a2';
	push @EXPORT, 'cbrr2tbl';
	push @EXPORT, 'cbrr_precision';
	push @EXPORT, 'ccmpred2rr';
	push @EXPORT, 'clash_count';
	push @EXPORT, 'contacts_with_min_all_atom_dist';
	push @EXPORT, 'count_lines';
	push @EXPORT, 'count_satisfied_tbl_rows';
	push @EXPORT, 'count_ss_element';
	push @EXPORT, 'count_ss_match';
	push @EXPORT, 'count_ss';
	push @EXPORT, 'coverage_tbl';
	push @EXPORT, 'dssp_result';
	push @EXPORT, 'dssp_ss_pdb';
	push @EXPORT, 'each_residue_sec_pdb';
	push @EXPORT, 'each_residue_sec_ss';
	push @EXPORT, 'evaluate_using_TMscore';
	push @EXPORT, 'extract_chain';
	push @EXPORT, 'fasta2residues_hash';
	push @EXPORT, 'flatten_fasta';
	push @EXPORT, 'get_cns_energy';
	push @EXPORT, 'load_pdb';
	push @EXPORT, 'min_all_atom_dist_in_top_model';
	push @EXPORT, 'noe_tbl_violation_coverage';
	push @EXPORT, 'parse_pdb_row';
	push @EXPORT, 'pdb2rr';
	push @EXPORT, 'pdb_each_residue_atoms';
	push @EXPORT, 'print2file';
	push @EXPORT, 'print2line';
	push @EXPORT, 'print_log';
	push @EXPORT, 'psipred2ss';
	push @EXPORT, 'reindex_chain_with_gaps';
	push @EXPORT, 'reindex_chain';
	push @EXPORT, 'reindex_rr';
	push @EXPORT, 'res_num_res_name';
	push @EXPORT, 'rr2chimera_ca';
	push @EXPORT, 'rr2chimera_cb';
	push @EXPORT, 'rr2contacts_hash';
	push @EXPORT, 'rr_coverage';
	push @EXPORT, 'rr_remove_unsatisfied_top_model';
	push @EXPORT, 'rr_rows_thresholds';
	push @EXPORT, 'seq_chain_with_gaps';
	push @EXPORT, 'seq_chain';
	push @EXPORT, 'seq_fasta';
	push @EXPORT, 'seq_rr';
	push @EXPORT, 'seq_ss';
	push @EXPORT, 'ss2residues_hash';
	push @EXPORT, 'ss_line_ss';
	push @EXPORT, 'ssnoe_tbl_min_pdb_dist';
	push @EXPORT, 'strand_count';
	push @EXPORT, 'sum_noe_dev';
	push @EXPORT, 'system_cmd';
	push @EXPORT, 'tbl2rows_hash';
	push @EXPORT, 'tmscore_using_TMscore';
	push @EXPORT, 'trim_fasta';
	push @EXPORT, 'trim_pdb';
	push @EXPORT, 'write_cns_seq';
}

# http://proteopedia.org/wiki/index.php/Standard_Residues and http://prody.csb.pitt.edu/manual/reference/atomic/flags.html
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;
our @AA3;
our @AA1;

foreach my $residue (keys %AA3TO1) {
    push(@AA3,$residue);
}
foreach my $residue (keys %AA1TO3) {
    push(@AA1,$residue);
}

my $program_tmscore = "/var/www/cgi-bin/confold/bin/TMscore";
my $program_dssp = "/var/www/cgi-bin/confold/bin/dssp-2.0.4-linux-amd64";

our $ERROR_ANGLE = 2.0;
our $ERROR_DISTANCE = 0.2;

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "EXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		print "\nCOMMAND   : [$command]\n";
		print "PWD       : ".getcwd."\n";
		print "EXIT-CODE : $exit_code\n";
		print "ERROR MSG : $!\n";
		print "TIME      : ".(localtime)."\n";
		if(defined $log){
			print "LOG       : $log\n";
			system("tail -n 100 $log");
		}
		confess "EXITING!!\n";
	}
}
sub count_lines{
	my $file = shift;
	my $lines = 0;
	return "N/A" if not -f $file;
	open FILE, $file or confess "ERROR! Could not open $file! $!";
	while (<FILE>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		next if not defined $_;
		next if length($_) < 1;
		$lines ++;
	}	
	close FILE;
	return $lines;
}
sub print_log{
	my $file_log = shift;
	my $message = shift;
	if (-f $file_log){
		open  LOG, ">>$file_log" or confess $!;
		print LOG "[".(localtime)."] ".$message."\n";
		close LOG;
	}
	else{
		open  LOG, ">$file_log" or confess $!;
		print LOG "[".(localtime)."] ".$message."\n";
		close LOG;
	}
}
sub print2file{
	my $file = shift;
	my $message = shift;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message."\n";
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message."\n";
		close FILE;
	}
}
sub print2line{
	my $file = shift;
	my $message = shift;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message;
		close FILE;
	}
}
sub tbl2rows_hash{
	my $file_tbl = shift;
	confess ":(" if not -f $file_tbl;
	my %noe = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		$_ =~ s/^\s+//;
		next if $_ !~ m/^assign/;
		#next if $_ =~ m/or/;
		$_ =~ s/^\s+//;
		$_ =~ s/\(/ /g;
		$_ =~ s/\)/ /g;
		$noe{$_} = 1;
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe;
	return %noe;
}
sub ssnoe_tbl_min_pdb_dist{
	my $file_tbl = shift;
	my $file_pdb = shift;
	confess ":(" if not -f $file_tbl;
	confess ":(" if not -f $file_pdb;
	my %noe_hash = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		chomp $_;
		$_ =~ s/^\s+//;
		$_ =~ s/\)/ /g;
		$_ =~ s/\(/ /g;
		confess "$_" if $_ !~ m/^assign/;
		my @C = split /\s+/, $_;
		if ($C[6] eq "or" and $C[17] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16]." ".$C[19]." ".$C[22];
			$noe_hash{$_}{"distance"} = $C[23]." ".$C[24]." ".$C[25];
		}
		elsif ($C[6] eq "or" and $C[17] ne "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		elsif ($C[6] ne "or" and $C[11] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10]." ".$C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		else{
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10];
			$noe_hash{$_}{"distance"} = $C[11]." ".$C[12]." ".$C[13];
		}
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe_hash;
	my %xyzPDB = all_xyz_pdb($file_pdb);
	foreach (sort keys %noe_hash){
		my $left = $noe_hash{$_}{"left"};
		my $right = $noe_hash{$_}{"right"};
		my $distance = $noe_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my %left_list = ();
		my %right_list = ();
		for(my $i = 0; $i <= $#L; $i = $i+2){
			$left_list{$L[$i]." ".$L[$i+1]} = 1;
		}
		for(my $i = 0; $i <= $#R; $i = $i+2){
			$right_list{$R[$i]." ".$R[$i+1]} = 1;
		}
		my $distance_pdb = 1000.0;
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $xyzPDB{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $xyzPDB{$R[0]." ".uc($R[1])};
				my $d = calc_dist($xyzPDB{$L[0]." ".uc($L[1])}, $xyzPDB{$R[0]." ".uc($R[1])});
				$distance_pdb = $d if $d < $distance_pdb;
			}
		}
		$noe_hash{$_}{"pdb_distance"} = $distance_pdb;
	}
	return %noe_hash;
}
sub coverage_tbl{
	my $file_fasta = shift;
	my $file_tbl = shift;
	my $flag_dihe = shift;
	confess ":(" if not -f $file_fasta;
	confess ":(" if not -f $file_tbl;
	$flag_dihe = 0 if not defined $flag_dihe;
	my $seq = seq_fasta($file_fasta);
	my $L = length $seq;
	my $cov = $seq;
	$cov =~ s/[A-Z]/-/g;
	my %noe = tbl2rows_hash($file_tbl);
	foreach (keys %noe){
		#assign (resid 123 and name CA) (resid  58 and name CA) 3.60 0.10 3.40
		my @C = split /\s+/, $_;
		my $r1 = $C[2]; my $r2 = $C[7];
		# the case when "ca or cb" restraints are provided
		#assign ((resid 123 and name ca) or (resid 123 and name cb)) ((resid 58 and name ca) or (resid 58 and name cb)) 3.60 0.10 3.40
		$r2 = $C[13] if $C[6] eq "or";
		$r2 = $C[17] if $flag_dihe;
		my $c1 = substr $cov, ($r1 - 1), 1;
		my $c2 = substr $cov, ($r2 - 1), 1;
		if ($c1 eq "-" ){
			$c1 = 1;
		}
		elsif ($c1 eq "*" ){
			$c1 = "*";
		}
		else{
			$c1++;
			$c1 = "*" if ($c1 > 9);
		}
		if ($c2 eq "-" ){
			$c2 = 1;
		}
		elsif ($c2 eq "*" ){
			$c2 = "*";
		}
		else{
			$c2++;
			$c2 = "*" if ($c2 > 9);
		}
		substr $cov, ($r1 - 1), 1, $c1;
		substr $cov, ($r2 - 1), 1, $c2;
	}
	my $cov2 = $cov;
	$cov2 =~ s/-//g;
	return sprintf "$cov [%12s : %3s restraints touching %s residues]", basename($file_tbl), (scalar keys %noe), length($cov2);
}
sub count_ss_match{
	my $file_pdb = shift;
	my $file_sec = shift;
	my $se = shift; # secondary structure element
	confess ":(" if not -f $file_pdb;
	confess ":(" if not -f $file_sec;
	confess ":(" if not defined $se;
	my %residue_ss1 = each_residue_sec_ss($file_sec);
	my %residue_ss2 = each_residue_sec_pdb($file_pdb);
	my $count = 0;
	foreach my $r (keys %residue_ss1){
		next if $residue_ss1{$r} ne $se;
		$count++ if $residue_ss1{$r} eq $residue_ss2{$r};
	}
	return $count;
}
sub count_ss_element{
	my $file_pdb_or_ss = shift;
	confess ":(" if not -f $file_pdb_or_ss;
	my $se = shift; # secondary structure element
	my %residue_ss = each_residue_ss($file_pdb_or_ss);
	foreach my $r (keys %residue_ss){
		delete $residue_ss{$r} if $residue_ss{$r} ne $se;
	}
	return (scalar keys %residue_ss);
}
sub count_satisfied_tbl_rows{
	my $file_pdb = shift;
	my $file_tbl = shift;
	my $type    = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not -f $file_tbl;
	confess ":(" if not ($type eq "noe" or $type eq "dihedral");
	my $count = 0;
	my $total = 0;
	my %log_rows = ();
	if ($type eq "dihedral"){
		my %noe = tbl2rows_hash($file_tbl);;
		my %phi = dssp_result($file_pdb, "phi");
		my %psi = dssp_result($file_pdb, "psi");
		foreach (sort keys %noe){
			chomp $_;
			#assign resid   1 and name ca resid   2 and name ca resid   3 and name ca resid   4 and name ca  5.0 -164.6 5.0 2 ! alpha for 2
			my @C = split /\s+/, $_;
			my $angle_true = 0.0;
			if (uc($C[5]) eq "C" and uc($C[10]) eq "N" and uc($C[15]) eq "CA" and uc($C[20]) eq "C"){
				$angle_true = $phi{$C[2]};
			}
			elsif (uc($C[5]) eq "N" and uc($C[10]) eq "CA" and uc($C[15]) eq "C" and uc($C[20]) eq "N"){
				$angle_true = $psi{$C[2]};
			}
			else{
				confess "Undefined dihedral angle $_";
			}
			my $viol_flag = 1;
			my $d = abs($angle_true - $C[22]);
			$d = abs(360.0 - $d) if $d > 180.0;
			if ($d < ($C[23] + 2.0)){
				$count ++;
				$viol_flag = 0;
			}
			$total++;
			$log_rows{(sprintf "%3s\t%.2f\t%.2f # $_", $viol_flag, $d, $angle_true)} =  $viol_flag;
		}
	}
	else{
		my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
		foreach (sort keys %tbl_hash){
			my $viol_flag = 1;
			my $distance = $tbl_hash{$_}{"distance"};
			my @D = split /\s+/, $distance;
			my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
			my $deviation = $pdb_distance - ( $D[0] + $D[2] );
			if ($pdb_distance < ( $D[0] + $D[2] + 0.2) ){
				$count ++;
				$viol_flag = 0;
				$deviation = 0.0;
			}
			if ($pdb_distance < ( $D[0] - $D[1] - 0.2) ){
				$count--;
				$viol_flag = 1;
				$deviation = -($D[0] - $D[1] - $pdb_distance);
			}
			
			$log_rows{(sprintf "%3s\t%.2f\t%.2f # $_", $viol_flag, $deviation, $pdb_distance)} = $viol_flag;
			$total++;
		}
	}
	open NOEDEV, ">".basename($file_tbl,".tbl")."_violation.txt" or confess $!;
	print NOEDEV "#NOE violation check; $file_pdb against $file_tbl\n";
	print NOEDEV "#violation-flag, deviation, actual-measurement, Input-NOE-restraint\n";
	foreach (sort {$log_rows{$b} <=> $log_rows{$a}} keys %log_rows){
		print NOEDEV $_."\n";
	}
	close NOEDEV;
	return $count."/".$total;
}
sub noe_tbl_violation_coverage{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not -f $file_tbl;
	my $cov = seq_chain($file_pdb);
	$cov =~ s/[A-Z]/-/g;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	my %xyz = all_xyz_pdb($file_pdb);
	foreach (sort keys %tbl_hash){
		my $left = $tbl_hash{$_}{"left"};
		my $right = $tbl_hash{$_}{"right"};
		my $distance = $tbl_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		substr $cov, $L[0] - 1, 1, "f" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $R[0] - 1, 1, "f" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $L[0] - 1, 1, "f" if $pdb_distance < $D[0] - $D[1] - 0.2;
		substr $cov, $R[0] - 1, 1, "f" if $pdb_distance < $D[0] - $D[1] - 0.2;
	}
	return $cov;
}
sub sum_noe_dev{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not -f $file_tbl;
	my $sum_dev = 0.0;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	foreach (sort keys %tbl_hash){
		my $viol_flag = 1;
		my @D = split /\s+/, $tbl_hash{$_}{"distance"};
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		if ($pdb_distance > ( $D[0] + $D[2] + 0.2) ){
			$sum_dev += $pdb_distance - ( $D[0] + $D[2] );
		}
		if ($pdb_distance < ( $D[0] - $D[1] - 0.2) ){
			$sum_dev += ( $D[0] - $D[1] ) - $pdb_distance;
		}
	}
	return sprintf "%.2f", $sum_dev;
}
sub get_cns_energy{
	my $cns_pdb = shift;
	my $energy_term = shift; # "overall", "vdw", "bon", "noe"
	confess ":(" if not -f $cns_pdb;
	my $value = "X";
	open CHAIN, $cns_pdb or confess $!;
	while(<CHAIN>){
		chomp $_;
		next if $_ !~ m/^REMARK\ $energy_term/;
		$_ =~ s/\s+//g;
		my @C = split /=/, $_;
		$value = $C[1];
	}
	close CHAIN;
	return int($value);
}
sub pdb_each_residue_atoms{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_atoms = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_atoms{parse_pdb_row($_,"rnum")} = "";
	}
	close CHAIN;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_atoms{parse_pdb_row($_,"rnum")} .= parse_pdb_row($_,"aname")." ";
	}
	close CHAIN;
	confess "ERROR! no atoms in file $chain!" if not scalar keys %rnum_atoms;
	return %rnum_atoms;
}
sub load_pdb{
	my $dir_chains = shift;
	confess ":( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdbList = <$dir_chains/*.pdb>;
	if(not (@pdbList)){
		@pdbList = <$dir_chains/*.ent>;
	}
	confess "ERROR! Directory $dir_chains has no pdb files!\n" unless(@pdbList);
	return @pdbList;
}
sub res_num_res_name{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}
sub ca_xyz_pdb{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %xyz_pdb = all_xyz_pdb($chain);
	my %rnum_rname = res_num_res_name($chain);
	my %ca_xyz = ();
	foreach (keys %xyz_pdb){
		my @C = split /\s+/, $_;
		next if $C[1] ne "CA";
		$ca_xyz{$C[0]} = $xyz_pdb{$_};
	}	
	confess ":(" if not scalar keys %ca_xyz;
	return %ca_xyz;
}
sub cb_xyz_pdb{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %xyz_pdb = all_xyz_pdb($chain);
	my %rnum_rname = res_num_res_name($chain);
	my %cb_xyz = ();
	foreach (keys %xyz_pdb){
		my @C = split /\s+/, $_;
		$cb_xyz{$C[0]} = $xyz_pdb{$_} if $C[1] eq "CB";
		$cb_xyz{$C[0]} = $xyz_pdb{$_} if ($rnum_rname{$C[0]} eq "GLY" and $C[1] eq "CA");
	}
	# if CB is not available, get the co-ordinates of CA
	# since, it is a chain, no residue number can have missing co-ordinate
	foreach (keys %rnum_rname){
		$cb_xyz{$_} = $xyz_pdb{$_." CA"} if not defined $cb_xyz{$_};
		confess ":( [$chain has no CA/CB for residue $_]" if not defined $cb_xyz{$_};
	}
	confess "ERROR! file $chain has no Cb atoms!" if not scalar keys %cb_xyz;
	return %cb_xyz;
}
sub all_xyz_pdb{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "ERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	return %xyz_pdb;
}
sub pdb2rr{
	my $chain = shift;
	my $rr = shift;
	my $seq_separation = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	$seq_separation = 6 if not defined $seq_separation;
	my $dist_threshold = 8.0;
	my %contacts_pdb = ();
	my %cb_xyz1 = cb_xyz_pdb($chain);
	my %cb_xyz2 = %cb_xyz1;
	foreach my $r1 (sort keys %cb_xyz1){
		foreach my $r2 (sort keys %cb_xyz2){
			next if ($r1 >= $r2);
			next if (($r2 - $r1) < $seq_separation);
			my $d = calc_dist($cb_xyz1{$r1}, $cb_xyz2{$r2});
			$contacts_pdb{"$r1 $r2"} = sprintf "%.3f", $d if ($d <= $dist_threshold);
		}
	}
	open RR, ">$rr" or confess $!;
	print RR "".seq_chain_with_gaps($chain)."\n";
	foreach (sort keys %contacts_pdb){
		print RR "$_ 0 8 ".$contacts_pdb{$_}."\n";
	}
	close RR;
}
sub seq_chain_with_gaps{
	my $chain = shift;
	my $flag = shift; # flag 1 if the left side dashes of the sequence are not wanted
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $start = 1;
	# if flagged, keep start for trimming
	if (defined $flag){
		open CHAIN, $chain or confess $!;
		while(<CHAIN>){
			next if $_ !~ m/^ATOM/;
			next unless (parse_pdb_row($_,"aname") eq "CA");
			confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
			$start = parse_pdb_row($_,"rnum");
			last;
		}
		close CHAIN;
	}
	# 1.find end residue number
	my $end;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		$end = parse_pdb_row($_,"rnum");
	}
	close CHAIN;
	# 2.initialize
	my $seq = "";
	for (my $i = 1; $i <= $end; $i++){
		$seq .= "-";
	}
	# 3.replace with residues
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $rnum = parse_pdb_row($_,"rnum");
		$rnum =~ s/[A-G]//g;
		substr $seq, ($rnum - 1), 1, $AA3TO1{parse_pdb_row($_,"rname")}; 
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return (substr $seq, $start - 1);
}
sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}
sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}
sub tmscore_using_TMscore{
	my $predicted = shift;
	my $native = shift;
	confess ":(" if not -f $predicted;
	confess ":(" if not -f $native;
	my %results = evaluate_using_TMscore($predicted, $native);
	return $results{"tm-score"};
}
sub evaluate_using_TMscore{
	my $predicted = shift;
	my $native = shift;
	confess ":( predicted does not exit!" if not -f $predicted;
	confess ":( native does not exist!" if not -f $native;
	my @results = `$program_tmscore $predicted $native | grep -e RMSD\\ of -e TM-score\\ \\ \\  -e MaxSub-score -e GDT-TS-score -e GDT-HA-score`;
	if (not defined $results[0]){
		print "Executing: [$program_tmscore $predicted $native]\n";
		system("$program_tmscore $predicted $native");
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
sub clash_count{
	my $file_pdb = shift;
	my $threshold = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not defined $threshold;
	my $count = 0;
	my %ca_xyz = ca_xyz_pdb($file_pdb);
	foreach my $r1 (sort keys %ca_xyz){
		my @R1 = split /\s+/, $ca_xyz{$r1};
		my $x1 = $R1[0]; my $y1 = $R1[1]; my $z1 = $R1[2];
		foreach my $r2 (sort keys %ca_xyz){
			next if $r1 >= $r2;
			my @R2 = split /\s+/, $ca_xyz{$r2};
			my $x2 = $R2[0]; my $y2 = $R2[1]; my $z2 = $R2[2];
			my $d = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			if ($d <= $threshold){
				$count++;
			}
		}
	}
	return $count;
}
sub reindex_chain{
	my $file_pdb = shift;
	my $index = shift;
	my $out_pdb = shift;
	confess "ERROR! file $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! index $index is invalied!" if not defined $index;
	open PDBFILE, $file_pdb or confess $!;
	my @pdb_lines = <PDBFILE>;
	close PDBFILE;
	# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
	my $res_counter = $index - 1;
	my $atom_counter = 0;
	my $prev_res_num = "XX";
	open OUTPDB, ">$out_pdb" or confess $!;
	foreach (@pdb_lines) {
		next if $_ !~ m/^ATOM/;
		next if not ((parse_pdb_row($_,"altloc") eq "") or (parse_pdb_row($_,"altloc") eq "A"));
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		my $this_rnum = parse_pdb_row($_,"rnum");
		if ($prev_res_num ne $this_rnum) {
			$prev_res_num = $this_rnum;
			$res_counter++;
		}
		$atom_counter++;
		my $rnum_string = sprintf("%4s", $res_counter);
		my $anum_string = sprintf("%5s", $atom_counter);
		my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
		print OUTPDB $row;
	}
	print OUTPDB "END\n";
	close OUTPDB;
}
sub reindex_chain_with_gaps{
	my $file_pdb = shift;
	my $index = shift;
	my $out_pdb = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not defined $index;
	open PDBFILE, $file_pdb or confess $!;
	my @pdb_lines = <PDBFILE>;
	close PDBFILE;
	# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
	my $first_residue;
	foreach (@pdb_lines) {
		next if $_ !~ m/^ATOM/;
		next if not ((parse_pdb_row($_,"altloc") eq "") or (parse_pdb_row($_,"altloc") eq "A"));
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		$first_residue = parse_pdb_row($_,"rnum");
		last;
	}
	my $atom_counter = 0;
	open OUTPDB, ">$out_pdb" or confess $!;
	foreach (@pdb_lines) {
		next if $_ !~ m/^ATOM/;
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		my $alt_loc = parse_pdb_row($_,"altloc");
		my $rnum = parse_pdb_row($_,"rnum");
		my $rnum_digits = $rnum;
		my $rnum_alphabets = $rnum;
		$rnum_alphabets =~ s/[0-9]//g;
		$rnum_digits =~ s/[A-Z]//g;
		$atom_counter++;
		my $rnum_string = sprintf "%4s", ($index + ($rnum_digits - $first_residue)).$rnum_alphabets;
		my $anum_string = sprintf "%5s", $atom_counter;
		my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
		print OUTPDB $row;
	}
	print OUTPDB "END\n";
	close OUTPDB;
}
sub extract_chain{
	my $file_pdb = shift;
	my $out_pdb = shift;
	my $chain_id = shift;
	confess ":(" if not -f $file_pdb;
	open PDBFILE, $file_pdb or confess $!;
	open OUT, ">$out_pdb" or confess $!;
	while(<PDBFILE>) {
		last if $_ =~ m/ENDMDL/;
		next if $_ !~ m/^ATOM/;
		next if $chain_id ne "".parse_pdb_row($_,"chain");
		print OUT $_;
	}
	close PDBFILE;
	close OUT;
}
sub trim_pdb{
	my $file_pdb  = shift;
	my $start  = shift;
	my $end    = shift;
	my $out_pdb = shift;
	open PDBFILE, $file_pdb or confess $!;
	my @pdb_lines = <PDBFILE>;
	close PDBFILE;
	open OUTPDB, ">$out_pdb" or confess $!;
	foreach (@pdb_lines) {
		next if $_ !~ m/^ATOM/;
		next if not defined $AA3TO1{parse_pdb_row($_,"rname")};
		my $this_rnum = parse_pdb_row($_,"rnum");
		next if $this_rnum < $start;
		next if $this_rnum > $end;
		print OUTPDB $_;
	}
	close OUTPDB;
}
sub each_residue_sec_pdb{
	my $file_pdb = shift;
	my %res = ();
	my $ss = dssp_ss_pdb($file_pdb);
	foreach (my $i = 0; $i <= length($ss); $i++){
		$res{$i+1} = substr $ss, $i, 1;
	}
	return %res;
}
sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}
sub rr2contacts_hash{
	my $file_rr = shift;
	my $seq_sep = shift;
	my $count = shift;
	$seq_sep = 1 if not defined $seq_sep;
	$count = 100000 if not defined $count;
	confess ":(" if not -f $file_rr;
	my %contacts = ();
	open RR, $file_rr or confess $!;
	while(<RR>){
		next unless ($_ =~ /[0-9]/);
		$_ =~ s/^\s+//;
		next unless ($_ =~ /^[0-9]/);
		my @C = split(/\s+/, $_);
		confess "ERROR! Expecting a pair in row [".$_."]!\n" if (not defined $C[0] || not defined $C[1]);
		next if (abs($C[1] - $C[0]) < $seq_sep);
		if(defined $C[4]){
			$contacts{$C[0]." ".$C[1]} = $C[4];
		}
		elsif(defined $C[2] && $C[2] != 0){
			$contacts{$C[0]." ".$C[1]} = $C[2];
		}
		else{
			confess "ERROR! Confidence column not defined in row [".$_."] in file $file_rr!\n";
		}
		last if (scalar keys %contacts) == $count;
	}
	close RR;
	return %contacts;
}
sub reindex_rr{
	my $file_rr = shift;
	my $map_x_to_1 = shift;
	my $seq_len = shift;
	my $file_out_rr = shift;
	confess ":(" if not -f $file_rr;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	open OUT, ">$file_out_rr" or confess $!;
	while(<RR>){
		chomp $_;
		if ($_ =~ m/^[A-Z]/){
			print OUT (substr $_, $map_x_to_1-1, $seq_len)."\n";
		}
		else{
			my @C = split /\s+/, $_;
			print OUT ($C[0]-$map_x_to_1+1)." ".($C[1]-$map_x_to_1+1)." ".$C[2]." ".$C[3]." ".$C[4]."\n";
			# I assume that contacts not in the region interest are absent in the input already!
			confess ":(" if $C[0]-$map_x_to_1+1 > $seq_len;
			confess ":(" if $C[1]-$map_x_to_1+1 > $seq_len;
		}
	}
	close OUT;
	close RR;
}
sub seq_rr{
	my $file_rr = shift;
	confess ":(" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$_ =~ s/^\s+//;
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	confess ":( no sequence header in $file_rr" if not defined $seq;
	return $seq;
}
sub carr2tbl{
	my $file_rr  = shift;
	my $file_tbl = shift;
	confess ":(" if not -f $file_rr;
	my %contacts = rr2contacts_hash($file_rr);
	confess "ERROR! Empty contacts file!" if not scalar keys %contacts;
	my %rr_thres = rr_rows_thresholds($file_rr);
	open NOE, ">$file_tbl" or confess $!;
	foreach (keys %contacts){
		my @C = split /\s+/, $_;
		my $distance = sprintf("%.2f", 3.6);
		my $negdev	 = sprintf("%.2f", 0.1);
		confess "ERROR! threshold not defined for ".$C[0]." ".$C[1] if not defined $rr_thres{$C[0]." ".$C[1]};
		my $posdev	 = sprintf("%.2f", ($rr_thres{$C[0]." ".$C[1]} - 3.6));
		printf NOE "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f\n", $C[0], "ca", $C[1], "ca", $distance, $negdev, $posdev;
	}
	close NOE;
}
sub cbrr2tbl{
	my $file_rr  = shift;
	my $file_tbl = shift;
	confess ":(" if not -f $file_rr;
	my %r1a1r2a2 = cbrr2r1a1r2a2($file_rr, 1);
	my %rr_thres = rr_rows_thresholds($file_rr);
	my $seq = seq_rr($file_rr);
	open NOE, ">$file_tbl" or confess $!;
	foreach (sort {$r1a1r2a2{$b} <=> $r1a1r2a2{$a}} keys %r1a1r2a2){
		my @C = split /\s+/, $_;
		my $distance = sprintf("%.2f", 3.6);
		my $negdev	 = sprintf("%.2f", 0.1);
		confess "ERROR! threshold not defined for ".$C[0]." ".$C[2] if not defined $rr_thres{$C[0]." ".$C[2]};
		my $posdev	 = sprintf("%.2f", ($rr_thres{$C[0]." ".$C[2]} - 3.6));
		printf NOE "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f\n", $C[0], $C[1], $C[2], $C[3], $distance, $negdev, $posdev;
	}
	close NOE;
}
sub rr_rows_thresholds{
	my $file_rr = shift;
	confess ":(" if not -f $file_rr;
	my %rr_thres = ();
	open RR, $file_rr or confess $!;
	while(<RR>){
		next unless ($_ =~ /[0-9]/);
		$_ =~ s/^\s+//;
		next unless ($_ =~ /^[0-9]/);
		my @C = split(/\s+/, $_);
		confess "ERROR! Expecting a pair in row [".$_."]!\n" if (not defined $C[0] || not defined $C[1]);
		confess "ERROR! Confidence column not defined in row [".$_."] in file $file_rr!\n" if not defined $C[4];
		$rr_thres{$C[0]." ".$C[1]} = $C[3];
	}
	close RR;
	return %rr_thres;
}
sub cbrr2r1a1r2a2{
	my $file_rr  = shift;
	my $flag_confidence = shift;
	confess ":(" if not -f $file_rr;
	my %contacts = rr2contacts_hash($file_rr);
	my %r1a1r2a2 = ();
	my $seq = seq_rr($file_rr);
	foreach (keys %contacts) {
		my @pair = split(/\s+/,$_);
		my $r1 = substr($seq, ($pair[0] -1), 1);
		my $r2 = substr($seq, ($pair[1] -1), 1);
		my $ca1 = "cb";
		my $ca2 = "cb";
		$ca1 = "ca" if($r1 eq "G");
		$ca2 = "ca" if($r2 eq "G");
		$r1a1r2a2{$pair[0]." ".$ca1." ".$pair[1]." ".$ca2} = $pair[0]+$pair[1];
		$r1a1r2a2{$pair[0]." ".$ca1." ".$pair[1]." ".$ca2} = $contacts{$_} if defined $flag_confidence;
	}
	return %r1a1r2a2;
}
sub rr2chimera_ca{
	my $file_rr = shift;
	my $file_chimera = shift;
	confess ":(" if not -f $file_rr;
	my %r1a1r2a2 = rr2contacts_hash($file_rr);
	open CHIMERA, ">$file_chimera" or confess $!;
	foreach (sort {$r1a1r2a2{$b} <=> $r1a1r2a2{$a}} (keys %r1a1r2a2)){
		my @C = split /\s+/, $_;
		print CHIMERA "distance :".$C[0]."\@ca :".$C[2]."\@ca\n" ;
	}
	close CHIMERA;
}
sub rr2chimera_cb{
	my $file_rr = shift;
	my $file_chimera = shift;
	confess ":(" if not -f $file_rr;
	my %r1a1r2a2 = cbrr2r1a1r2a2($file_rr);
	open CHIMERA, ">$file_chimera" or confess $!;
	foreach (sort {$r1a1r2a2{$a} <=> $r1a1r2a2{$b}} (keys %r1a1r2a2)){
		my @C = split /\s+/, $_;
		print CHIMERA "distance :".$C[0]."\@".$C[1]." :".$C[2]."\@".$C[3]."\n" ;
	}
	close CHIMERA;
}
sub rr_coverage{
	my $file_rr = shift;
	my $seq = seq_rr($file_rr);
	confess ":(" if not -f $file_rr;
	my $cov = $seq;
	$cov =~ s/[A-Z]/-/g;
	my $lr = 0; my $mr = 0; my $sr = 0; my $nr = 0;
	my %contacts = rr2contacts_hash($file_rr);
	foreach my $pair (keys %contacts){
		my @C = split(/\s+/, $pair);
		my $r1 = $C[0]; my $r2 = $C[1];
		my $c1 = substr $cov, ($r1 - 1), 1;
		my $c2 = substr $cov, ($r2 - 1), 1;
		if ($c1 eq "-" ){
			$c1 = 1;
		}
		elsif ($c1 eq "*" ){
			$c1 = "*";
		}
		else{
			$c1++;
			$c1 = "*" if ($c1 > 9);
		}
		if ($c2 eq "-" ){
			$c2 = 1;
		}
		elsif ($c2 eq "*" ){
			$c2 = "*";
		}
		else{
			$c2++;
			$c2 = "*" if ($c2 > 9);
		}
		if(abs($C[0] - $C[1]) >= 24){
			$lr++;
		}
		elsif(abs($C[0] - $C[1]) >= 12){
			$mr++;
		}
		elsif(abs($C[0] - $C[1]) >= 6){
			$sr++;
		}
		else{
			$nr++;
		}
		substr $cov, ($r1 - 1), 1, $c1;
		substr $cov, ($r2 - 1), 1, $c2;
	}
	my $cov2 = $cov;
	$cov2 =~ s/-//g;
	return $cov." [".(scalar keys %contacts)." contacts, touching ".length($cov2)." residues, $lr LR, $mr MR, $sr SR, and $nr others]";
}
sub contacts_with_min_all_atom_dist{
	my $file_rr = shift;
	my $file_pdb = shift;
	confess ":(" if not -f $file_rr;
	confess ":(" if not -f $file_pdb;
	my %contacts = rr2contacts_hash($file_rr);
	confess ":(" if not %contacts;
	my %xyzPDB = all_xyz_pdb($file_pdb);
	my %pdb_atom_list = pdb_each_residue_atoms($file_pdb);
	foreach (keys %contacts){
		my @C = split /\s+/, $_;
		my %left_list = ();
		my %right_list = ();
		my @A1 = split /\s+/, $pdb_atom_list{$C[0]};
		my @A2 = split /\s+/, $pdb_atom_list{$C[1]};
		foreach (@A1){
			$left_list{$C[0]." ".$_} = 1;
		}
		foreach (@A2){
			$right_list{$C[1]." ".$_} = 1;
		}
		my $distance_pdb = 1000.0;
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $xyzPDB{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $xyzPDB{$R[0]." ".uc($R[1])};
				my @c1 = split(/\s+/,$xyzPDB{uc($L[0])." ".uc($L[1])});
				my $x1 = $c1[0]; my $y1 = $c1[1]; my $z1 = $c1[2];
				my @c2 = split(/\s+/,$xyzPDB{uc($R[0])." ".uc($R[1])});
				my $x2 = $c2[0]; my $y2 = $c2[1]; my $z2 = $c2[2];
				my $d = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
				$distance_pdb = $d if $d < $distance_pdb;
			}
		}
		$contacts{$_} = $distance_pdb;
	}
	return %contacts;
}
sub cbrr2cb_dist{
	my $file_rr = shift;
	my $file_pdb = shift;
	confess ":( $file_rr does not exist!" if not -f $file_rr;
	confess ":( $file_pdb does not exist!" if not -f $file_pdb;
	my %contacts = rr2contacts_hash($file_rr);
	confess ":(" if not %contacts;
	my %cb_xyz = cb_xyz_pdb($file_pdb);
	foreach (keys %contacts){
		my @C = split /\s+/, $_;
		next if not defined $cb_xyz{$C[0]};
		next if not defined $cb_xyz{$C[1]};
		confess "$file_pdb does not have ".$C[0]."\n" if not defined $cb_xyz{$C[0]};
		confess "$file_pdb does not have ".$C[1]."\n" if not defined $cb_xyz{$C[1]};
		my @c1 = split /\s+/, $cb_xyz{$C[0]};
		my $x1 = $c1[0]; my $y1 = $c1[1]; my $z1 = $c1[2];
		my @c2 = split /\s+/, $cb_xyz{$C[1]};
		my $x2 = $c2[0]; my $y2 = $c2[1]; my $z2 = $c2[2];
		my $d = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
		$contacts{$_} = sprintf "%.2f", $d;
	}
	confess ":(" if not scalar %contacts;
	return %contacts;
}
sub cbrr_precision{
	my $file_rr = shift;
	my $file_pdb = shift;
	my $file_log = shift;
	my %cb_dist = cbrr2cb_dist($file_rr, $file_pdb);
	my %rr_thres = rr_rows_thresholds($file_rr);
	my $satisfied = 0;
	system_cmd("rm -f $file_log");
	print2file($file_log, "DIST R1 R2");
	foreach (keys %cb_dist){
		$satisfied++ if $cb_dist{$_} <= $rr_thres{$_};
		print2file($file_log, sprintf "%-.3f $_", $cb_dist{$_});
	}
	return sprintf "%.3f", ($satisfied/(scalar keys %cb_dist));
}
sub rr_remove_unsatisfied_top_model{
	my $file_rr  = shift;
	my $file_pdb = shift;
	my $file_new_rr  = shift;
	my $file_log = shift;
	confess ":(" if not -f $file_rr;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not defined $file_new_rr;
	my %contacts = min_all_atom_dist_in_top_model($file_rr, $file_pdb, "log_min_all_atom_dist_in_top_model.txt");
	my %confidence = rr2contacts_hash($file_rr);
	my %rr_thres = rr_rows_thresholds($file_rr);
	system_cmd("rm -f $file_new_rr");
	system_cmd("rm -f $file_log");
	print2file($file_new_rr, seq_rr($file_rr));
	foreach (keys %contacts){
		print2file($file_new_rr, $_." 0 ".$rr_thres{$_}." ".$confidence{$_}) if $contacts{$_} <= $rr_thres{$_};
		print2file($file_log, "$_ ignored") if $contacts{$_} > $rr_thres{$_};
	}
}
sub min_all_atom_dist_in_top_model{
	# find the minimum of all atom to all atom distance between two residues, in model1.pdb
	my $file_rr = shift;
	my $file_pdb = shift;
	my $file_log  = shift;
	confess ":(" if not -f $file_rr;
	confess ":(" if not -f $file_pdb;
	my %contacts = rr2contacts_hash($file_rr);
	my %contacts_atom_list = ();
	confess ":(" if not %contacts;
	my %pdb_atom_list = pdb_each_residue_atoms($file_pdb);
	foreach (keys %contacts){
		my @C = split /\s+/, $_;
		my %left_list = ();
		my %right_list = ();
		my @A1 = split /\s+/, $pdb_atom_list{$C[0]};
		my @A2 = split /\s+/, $pdb_atom_list{$C[1]};
		foreach (@A1){
			$left_list{$C[0]." ".$_} = 1;
		}
		foreach (@A2){
			$right_list{$C[1]." ".$_} = 1;
		}
		$contacts_atom_list{$_." left"} = \%left_list;
		$contacts_atom_list{$_." right"} = \%right_list;
	}
	my %pdb_xyz = all_xyz_pdb($file_pdb);
	open LOG, ">$file_log" or confess $!;
	printf LOG "Minimum all atom distance in $file_pdb:";
	printf LOG "\n%-10s %-10s %-20s", "contact", "min_dist", "atom-pair";
	foreach (keys %contacts){
		my %left_list = %{$contacts_atom_list{$_." left"}};
		my %right_list = %{$contacts_atom_list{$_." right"}};
		my $min_dist = 1000;
		my $atom_pair = "";
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $pdb_xyz{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $pdb_xyz{$R[0]." ".uc($R[1])};
				my @c1 = split(/\s+/,$pdb_xyz{uc($L[0])." ".uc($L[1])});
				my $x1 = $c1[0]; my $y1 = $c1[1]; my $z1 = $c1[2];
				my @c2 = split(/\s+/,$pdb_xyz{uc($R[0])." ".uc($R[1])});
				my $x2 = $c2[0]; my $y2 = $c2[1]; my $z2 = $c2[2];
				my $dist = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
				if($min_dist > $dist){
					$min_dist = sprintf "%.2f", $dist;
					$atom_pair = $le." - ".$ri;
				}
			}
		}
		$contacts{$_} = $min_dist;
		printf LOG "\n%-10s %-10s %-20s", $_, $min_dist, $atom_pair;
	}
	close LOG;
	return %contacts;
}
sub ccmpred2rr{
	my $file_fasta = shift;
	my $file_ccmpred = shift;
	my $file_rr = shift;
	my %conf = ();
	open CCM, $file_ccmpred or confess $!;
	my $i = 1;
	while(<CCM>){
		my @C = split /\s+/, $_;
		for(my $j = 0; $j <= $#C; $j++){
			my $pair = $i." ".($j+1);
			$pair = ($j+1)." ".$i if ($j+1) < $i;
			my $confidence = $C[$j];
			$confidence = $conf{$pair} if (defined $conf{$pair} && $conf{$pair} > $confidence);
			$conf{$pair} = $confidence;
		}
		$i++;
	}
	close CCM;
	open RR, ">$file_rr" or confess $!;
	print RR "".seq_fasta($file_fasta)."\n";
	foreach (sort {$conf{$b} <=> $conf{$a}} keys %conf){
		my @C = split /\s+/, $_;
		next if abs($C[0] - $C[1]) < 6;
		print RR $_." 0 8 ".$conf{$_}."\n";
	}
	close RR;
}
sub seq_ss{
	my $file_ss = shift;
	confess ":(" if not -f $file_ss;
	my $header = `head -n 1 $file_ss`;
	confess "$file_ss does not have header!" if $header !~ m/\>/;
	my $seq_ss = `head -n 2 $file_ss | tail -n 1`;
	chomp $seq_ss;
	$seq_ss =~ tr/\r//d; # chomp does not remove \r
	return $seq_ss;
}
sub ss_line_ss{
	my $file_ss = shift;
	confess ":(" if not -f $file_ss;
	my $header = `head -n 1 $file_ss`;
	confess "$file_ss does not have header!" if $header !~ m/\>/;
	my $ss = `head -n 3 $file_ss | tail -n 1`;
	chomp $ss;
	return $ss;
}
sub psipred2ss{
	my $file_psipred = shift;
	my $out_file_ss = shift;
	confess ":(" if not -f $file_psipred;
	my $seq = ""; 
	my $ss  = "";
	open INPUT, $file_psipred or confess $!;
	while(<INPUT>){
		$_ =~ s/^\s+//;
		next if $_ !~ m/^[0-9]/;
		my @columns = split(/\s+/, $_);
		$seq .= $columns[1];
		$ss .= $columns[2];
	}
	close INPUT;
	open OUTPUT, ">$out_file_ss" or confess $!;
	print OUTPUT ">".basename($file_psipred)."\n";
	print OUTPUT "$seq\n";
	print OUTPUT "$ss\n";
	close OUTPUT;
}
sub strand_count{
	my $file_ss = shift;
	confess ":(" if not -f $file_ss;
	open INPUT, $file_ss or confess $!;
	my @rows = <INPUT>;
	chomp @rows;
	close INPUT;
	my $ss = $rows[2];
	my $count = 0;
	for(my $i = 0; $i <= length($ss); $i++){
		next if (substr $ss, $i, 1) ne "E";
		next if not defined (substr $ss, $i-1, 1);
		next if not defined (substr $ss, $i, 1);
		next if not defined (substr $ss, $i+1, 1);
		$count++ if ((substr $ss, $i-1, 1) eq "E" and (substr $ss, $i, 1) eq "E" and (substr $ss, $i+1, 1) ne "E");
	}
	return $count;
}
sub dssp_result{
	my $file_pdb = shift;
	my $selection = shift;
	confess ":(" if not -f $file_pdb;
	confess ":(" if not defined $selection;
	my %RESIDUE = ();
	my %SS = ();
	my %PHI = ();
	my %PSI = ();
	my $command = "$program_dssp $file_pdb | grep -C 0 -A 1000 \'  #  RESIDUE\' | tail -n +2";
	my @dssp_rows = `$command`;
	foreach(@dssp_rows){
		my $rnum = substr($_,  5, 6);
		my $res  = substr($_, 13, 1);
		my $sstr = substr($_, 16, 1);
		my $phia = substr($_,103, 6);
		my $psia = substr($_,109, 6);
		$rnum =~ s/\s+//g;
		$res  =~ s/\s+//g;
		$sstr =~ s/\s+//g;
		$phia =~ s/\s+//g;
		$psia =~ s/\s+//g;
		# alternate residue locations may have alphabets
		#$rnum =~ s/[^0-9]//g;
		$res  =~ s/[^A-Z]//g;
		$sstr =~ s/[^A-Z]//g;
		next if length($rnum) < 1;
		confess ":( residue not defined for $rnum" if length($res) < 1;
		confess ":( phi not defined for $rnum" if length($phia) < 1;
		confess ":( psi not defined for $rnum" if length($psia) < 1;
		$sstr = "C" if length($sstr) < 1;
		$RESIDUE{$rnum} = $res;
		$SS{$rnum}      = $sstr;
		$PHI{$rnum}     = $phia;
		$PSI{$rnum}     = $psia;
	}

	# special temporary change
	return if not scalar keys %PHI;

	confess "$file_pdb seems empty!" if not scalar keys %PHI;
	return %SS if ($selection eq "ss");
	return %PHI if ($selection eq "phi");
	return %PSI if ($selection eq "psi");
	confess "ERROR! Invalid selection string!";
}
sub dssp_ss_pdb{
	my $file_pdb = shift;
	confess ":(" if not -f $file_pdb;
	my %ss = dssp_result($file_pdb, "ss");
	my $ssrow = "";
	foreach(sort {$a <=> $b} keys %ss){
		$ssrow .= $ss{$_};
	}
	$ssrow =~ s/\./C/g;
	$ssrow =~ s/S/C/g;
	$ssrow =~ s/T/C/g;
	$ssrow =~ s/B/C/g;
	$ssrow =~ s/G/C/g;
	return $ssrow;
}
sub each_residue_sec_ss{
	my $file_ss = shift;
	my %res = ();
	my $ss = ss_line_ss($file_ss);
	foreach (my $i = 0; $i < length($ss); $i++){
		$res{$i+1} = substr $ss, $i, 1;
	}
	return %res;
}
sub count_ss{
	my $file_pdb = shift;
	my $char = shift;
	confess ":(" if not -f $file_pdb;
	my $count = 0;
	my %ssHash = dssp_result($file_pdb, "ss");
	foreach (keys %ssHash){
		$count++ if $ssHash{$_} eq $char;
	}
	return $count;
}
sub ss2residues_hash{
	my $file_ss = shift;
	confess ":(" if not -f $file_ss;
	my $seq = seq_ss($file_ss);
	my %res = ();
	foreach (my $i = 0; $i <= length($seq); $i++){
		$res{$i+1} = substr $seq, $i, 1;
	}
	return %res;
}
sub seq_fasta{
	my $file_fasta = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = "";
	open FASTA, $file_fasta or confess $!;
	while (<FASTA>){
		next if (substr($_,0,1) eq ">"); 
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$seq .= $_;
	}
	close FASTA;
	return $seq;
}
sub write_cns_seq{
	my $file_fasta = shift;
	my $file_cns_seq = shift;
	confess ":(" if not -f $file_fasta;
	my @seq = split //, seq_fasta($file_fasta);
	my $three_letter_seq;
	foreach (@seq) {
		$three_letter_seq .= $AA1TO3{$_}." ";
	}
	open SEQUENCE, ">$file_cns_seq" or confess $!;
	while($three_letter_seq){
		if(length ($three_letter_seq) <= 64 ){
			print SEQUENCE $three_letter_seq;
			print SEQUENCE "\n";
			$three_letter_seq = ();
		}
		else{
			print SEQUENCE substr($three_letter_seq, 0, 64);
			print SEQUENCE "\n";
			$three_letter_seq = substr($three_letter_seq, 64);
		}
	}
	close SEQUENCE;
}
sub flatten_fasta{
	my $file_fasta = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	open FASTA, $file_fasta or confess $!;
	my $header = <FASTA>;
	chomp $header;
	close FASTA;
	open FASTA, ">$file_fasta" or confess $!;
	print FASTA "$header\n";
	print FASTA "$seq\n";
	close FASTA;
}
sub trim_fasta{
	my $file_fasta = shift;
	my $start = shift;
	my $end = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	open FASTA, $file_fasta or confess $!;
	my $header = <FASTA>;
	chomp $header;
	close FASTA;
	my $trimmmed_seq = substr $seq, $start - 1, $end - $start + 1;
	open FASTA, ">$file_fasta" or confess $!;
	print FASTA "$header\n";
	print FASTA "$trimmmed_seq\n";
	close FASTA;
}
sub fasta2residues_hash{
	my $file_fasta = shift;
	confess ":(" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	my %res = ();
	foreach (my $i = 0; $i < length($seq); $i++){
		$res{$i+1} = substr $seq, $i, 1;
	}
	return %res;
}

1;
