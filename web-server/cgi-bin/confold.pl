#!/usr/bin/perl -w
# Badri Adhikari, 1/23/2015

# v0.6, 1/26/2015
# setting rep2 radius to 0.8 only

# v0.5, 1/17/2015
# all subroutines merged into confold.pm
# secondary structure restraints hard-coded

# v0.4, 1/14/2015
# if all contacts are provided rep2 is 0.8 and contact detection threshold is 6.5A
# setting atom select scheme 6 as the default scheme

# v0.3, 1/10/2015
# experiments with appropriate atom selection for distance geometry
# setting default type as 6 (backbone and cb and h)

# v0.2, 1/5/2015
# β-sheet detection code updated, up to 2 residues shifting considered, 10 types of configurations
# for β-sheet detection, hydrogen bond distance calculation added besides distance calculations

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use lib dirname(abs_path($0));
use Getopt::Long;
use confold;

my $param_info = <<EOF;
Build models using CONFOLD
Parameter      Type    Default Description 
dir_job      : string : ?    : Directory where three files are expected 
allrr        : flag   : 0    : If supplied uses all contacts in the rr file. 
                               If not, 10 sets of contacts are used, 0.4L, 0.6L, etc. up to 2.2L
norr         : flag   : 0    : If provided, contacts are not used to build models. 
                               If not, and contacts are absent, error is thrown   
selectrr     : flag   : 0    : If provided, 20 sets of contact contacts are not used to build models, top-0.4L to top-2.2L 
rr_type      : string : cb   : ca for Cα-Cα contacts and cb for Cβ-Cβ contacts
lambda       : float  : 0.5  : Choose set of secondary structure restraints, 0.4 for predictions and 1.0 for reconstructions
stage2       : flag   : 0    : Build stage2 models using contact filtering and beta sheet detection
sheet_detect : flag   : 0    : Build additional set of stage2 models using contact filtering only
con_filter   : flag   : 0    : Build additional set of stage2 models using beta sheet detection only
cont_wt      : float  : 10   : Weight for contact restraints
ss_wt1       : float  : 5    : First scheme of weighting of secondary structure restraints relative to contact restraints. 
                               1 means contact and secondary structure restraints have equal weights
ss_wt2       : float  : 0.5  : Second scheme of weighting of secondary structure restraints relative to contact restraints. 
                               A different weight is required.
no_ss_wt2    : flag   : 0    : If provided, the second weighting scheme is not used
cpu          : int    : 20   : Number of processors to use
run_betapro  : flag   : 0    : 
betapro      : string : -    : If betapro is already run, supply the str you have
pairing      : string : -    : Supply the pairing file you have
pair_thres   : string : -    : pairing threshold
flush_dirs   : flag   : 0    : If provided, deletes all existing directories and files
rep2         : string : -    : 0.8 OR 0.85
atom_scheme  : int    : 2    : Atom selection scheme for distance geometry
                               (1) as-is / ca + ha + n + hn + c + cb + cg (2) as-is + o (3) as-is + o + h 
                               (2) as-is + o (3) as-is + o + h (4) c + cα + n + o (backbone atoms)
                               (5) bkbone + cβ (6) bkbone + cβ + h (7)  bkbone + cβ + h + cg                              
                               Notes: 
                               (a) with option 6, error "%MATROT error encountered: Error in internal consistency check" was observed for pdbs 1QF9 and 3TGI
                               (b) with option 7, chirality issues were faced on target 1QJP
                               (c) with option 2 (best), maximum reconstruction accuracy and sheet_detection results for 150 fragfold proteins.

Pairing File Format:
 - 6 column file with columns a, b, c, d, t, and f
 - a-b and c-d are residue strands, for example 2-7 and 20-25
 - t is the pairing type (A or P), and f is the confidence of pairing
 - a must always be less than b, c must be less than d if parallel, and c must be greater than d if anti-parallel
EOF

my ($dir_job, $flag_all_rr, $flag_no_rr, $flag_select_rr, $rr_type, $lambda);
my ($flag_stage2, $flag_sheet_detect, $flag_con_filter, $con_wt, $ss_wt1);
my ($ss_wt2, $flag_no_ss_wt2, $cpu, $flag_run_betapro, $file_betapro); 
my ($file_pairing, $pair_caca_thres, $flag_flush, $atom_select_scheme, $rep2);

GetOptions(
	"dir_job=s"		=> \$dir_job,
	"allrr"			=> \$flag_all_rr,
	"norr"			=> \$flag_no_rr,
	"selectrr"		=> \$flag_select_rr,
	"rr_type=s"		=> \$rr_type,
	"lambda=s"		=> \$lambda,
	"stage2"		=> \$flag_stage2,
	"sheet_detect"	=> \$flag_sheet_detect, 
	"con_filter"	=> \$flag_con_filter, 
	"cont_wt=s" 	=> \$con_wt,
	"ss_wt1=s"		=> \$ss_wt1,
	"ss_wt2=s"		=> \$ss_wt2,
	"no_ss_wt2"		=> \$flag_no_ss_wt2,
	"cpu=i"			=> \$cpu,
	"runbetapro"	=> \$flag_run_betapro,
	"betapro=s"		=> \$file_betapro,
	"pairing=s"		=> \$file_pairing,
	"pair_thres=s"	=> \$pair_caca_thres,
	"flush_dirs"	=> \$flag_flush,
	"atom_scheme=i"	=> \$atom_select_scheme,
	"rep2=s"		=> \$rep2)
or confess "Error in command line arguments\n";

# apply_defaults
$rr_type = "cb" if !$rr_type and !$flag_no_rr; 
$lambda = 0.5 if !$lambda; # this actually does not matter very much
$con_wt = 10 if !$con_wt;
$ss_wt1 = 5 if !$ss_wt1;
$ss_wt2 = 0.5 if !$ss_wt2 and !$flag_no_ss_wt2;
$cpu = 20 if !$cpu;
$atom_select_scheme = 6 if !$atom_select_scheme;

# Fixed Global parameters, common to all stages, xL, and sec_wt
# Parameter     Description 
# model_count   if model_count is not provided, extended structure is only generated
# dihed_wt1     dihedral angle weight parameter 1 (default is noe_wt * class_n2_wt)
# dihed_wt2     dihedral angle weight parameter 2 (default is noe_wt * class_n2_wt)
# rep1          initial van der Waals repel radius (md.cool.init.rad)
# rep2          final van der Waals repel radius(md.cool.fina.rad)
# mini          number of minimization steps
# bprothres     betapro confidence threshold
# pairthres     distance threshold for detecting beta pairs

my $model_count = 20;
my $rep1 = 1.0;
#my $rep2 = 0.8; # 0.8 is good when using true contacts, but 0.85 for predicted contacts. For EVFOLD proteins, average accuracy is lower with 0.8.
my $mini = 15000;
my $bpro_thres = 0.1;
# while experimenting with true contacts, I found that in 8A, a strand next to pairing strand is also considered. So, setting it around 7A.
#my $pair_caca_thres = 7.0; 
my $pair_hbond_thres = $pair_caca_thres + 3.0; # while experimenting with true contacts, I found that in 8A, a strand next to pairing strand is also considered. So, setting it around 7A.
my $pair_policy = "single_model";

print "[$0]: ".(localtime)."\n";
validate_params();
print_params();

# Global variables and all job options
my ($file_seq, $file_sec, $file_rr, $seq_id, $seq_len);
my ($dir_cns_install, $dir_scripts);
my %stage_list = ();
my %sec_wt_list = ();
my %xL = ();
my %residues = ();
my ($dir_log, $dir_extended, $dir_betapro);
my ($program_betapro1, $program_betapro2);
load_global_variables();
print_global_variables();

# Global variables for generating secondary structure restraints
my %ATOMTYPE   = qw( CA 1 N 1 C 1 O 1 );
my %SHIFT      = qw( 0 1 +1 1 -1 1 );
my %restraints_dihedral = ();
my %restraints_strandOO = ();
my %restraints_distance = ();
my %restraints_hbond    = ();
load_ss_restraints($lambda, "$dir_log/sec_restraints.txt");

if ($flag_run_betapro){
	print "\nrun BETApro";
	if (not -f $file_betapro){
		system_cmd("cp $file_sec ./$seq_id.ss");
		system_cmd("$program_betapro1 $seq_id.ss betapro.$seq_id.matrix &> LOG_betapro.txt");
		system_cmd("$program_betapro2 betapro.$seq_id.matrix $seq_id &> LOG_betapro.txt");
		$file_betapro = "betapro/$seq_id.str";
		confess "ERROR! Expected BETApro str file not found!" if not -f $file_betapro;
	}
}
if($file_betapro){
	print "\nload BETApro pairing file";
	$file_pairing = "$dir_job/pairing.txt";
	str2pairing_info($file_betapro, $file_pairing);
	confess "ERROR! BETApro file conversion failed!" if not -f $file_pairing;
}

print "\nbuild extended mtf and pdb";
mkdir $dir_extended or die $! if not -d $dir_extended;
chdir $dir_extended;
system_cmd("cp $dir_job/*.fasta ./");
build_extended() if not -f "extended.pdb";

# build models, for all job options
# Local variables
my ($current_stage, $current_sec_wt, $current_topx, $current_file_sec, $current_job_dir, $current_file_pairing);
foreach (sort {$stage_list{$a} <=> $stage_list{$b}} keys %stage_list){
	$current_stage = $_;
	my $dir_stage = "$dir_job/$current_stage";
	mkdir $dir_stage or die $! if not -d $dir_stage;
	foreach (keys %sec_wt_list){
		$current_sec_wt = $_;
		foreach (sort keys %xL){
			$current_topx = $_;
			print "\n\nstarting $current_stage for sec-wt $current_sec_wt xL $current_topx ..";
			$| = 1;
			if ($current_stage ne "stage1"){
				my $stage1_top_model = "$dir_job/stage1/sec_${current_sec_wt}_${current_topx}L/${seq_id}_model1.pdb";
				my $stage1_models = "$dir_job/stage1/sec_${current_sec_wt}_${current_topx}L";
				print "\nwait for corresponding stage 1 to finish..";
				system_cmd("touch $dir_job/check.running");
				while(1){
					sleep 3;
					confess "ERROR! There is nothing running!" if not -f "$dir_job/check.running";
					if(-f $stage1_top_model){
						system_cmd("rm -f $dir_job/check.running");
						last;
					}
					confess "ERROR! There is no extended.mtf in stage1. Is the stage1 running at $stage1_models" if not -f "$stage1_models/extended.mtf";
					confess "ERROR! stage1 failed! Check at $stage1_models" if -f "$stage1_models/iam.failed";
					open JOBLOG, ">$dir_job/check.running" or confess $!;
					print JOBLOG "".(localtime)."\n";
					print JOBLOG " checking $stage1_models models.. looks like stage1 is still running..\n";
					close JOBLOG;
				}
			}
			print "\nwait for processor to be free..";
			while(1){
				sleep 3;
				confess "ERROR! dir_job deleted? $dir_job !" if not -d $dir_job;
				my @running = `ps -u bap54 -f | grep cns_solve`;
				last if $#running < $cpu;
			}
			$current_file_pairing = $file_pairing if $file_pairing;
			$current_job_dir = "$dir_stage/sec_${current_sec_wt}_${current_topx}L";
			mkdir $current_job_dir or die $! if not -d $current_job_dir;
			chdir $current_job_dir or die $!;
			if (-f "./${seq_id}_model1.pdb"){
				print "\nlast execution results already exist!";
			}
			else{
				system_cmd("cp $dir_extended/extended.* ./");
				system_cmd("cp $dir_extended/*.fasta ./");
				system_cmd("cp $dir_job/*.ss    ./");
				$current_file_sec = "$current_job_dir/$seq_id.ss";
				contact_restraints();
				sec_restraints();
				build_models_in_background();
			}
			# flush all temporary variables, just to be safe
			$current_topx     = undef;
			$current_file_sec = undef;
			$current_job_dir  = undef;
			$current_file_pairing = undef;;
		}
		$current_sec_wt = undef;
	}
	$current_stage = undef;
}

print "\nwait for all jobs to finish..";
while(1){
	sleep 5;
	my @count = `find $dir_job/ -name iam.running`;
	chomp @count;
	my $count = 0;
	foreach (@count){
		next if length($_) < 10;
		$count++ if $_ =~ m/iam.running/;
		confess "ERROR! failed at $_!" if $_ =~ m/failed/;
	}
	last if $count == 0;
}
print "\n\n[$0]: ".(localtime)."\n";
print "################################################################################\n";

sub validate_params{
	print_usage("job directory not supplied") if !$dir_job;
	print_usage("one of the options allrr, norr, selectrr must be provided") if (!$flag_all_rr and !$flag_no_rr and !$flag_select_rr);
	print_usage("dgsa atom selection needed") if !$atom_select_scheme;
	print_usage("job directory $dir_job does not exist") if not -d $dir_job;
	if (not $flag_no_rr){
		print_usage("contact type ca/cb not defined") if not $rr_type;
		print_usage("rr type must be ca or cb") if not ($rr_type eq "ca" or $rr_type eq "cb");
	}
	$dir_job = abs_path($dir_job);
	print_usage("lambda must be between 0.1 and 10.0") if ($lambda > 10.0 or $lambda < 0.1);
	print_usage("please use different secondary structure weights") if ($ss_wt2 and $ss_wt1 eq $ss_wt2);
	print_usage("BETApro file $file_betapro does not exist") if ($file_betapro and not -f $file_betapro);
	print_usage("pairing file $file_pairing does not exist") if ($file_pairing and not -f $file_pairing);
}

sub print_params{
	print "\nParameters:";
	print "\n dir_job           $dir_job          ";
	print "\n flag_all_rr       $flag_all_rr      " if $flag_all_rr;
	print "\n flag_no_rr        $flag_no_rr       " if $flag_no_rr;
	print "\n flag_select_rr    $flag_select_rr   " if $flag_select_rr;
	print "\n rr_type           $rr_type          " if $rr_type;
	print "\n lambda            $lambda            " if $lambda;
	print "\n flag_stage2       $flag_stage2      " if $flag_stage2;
	print "\n flag_sheet_detect $flag_sheet_detect" if $flag_sheet_detect;
	print "\n flag_con_filter   $flag_con_filter  " if $flag_con_filter;
	print "\n con_wt            $con_wt           " if $con_wt;
	print "\n ss_wt1            $ss_wt1           " if $ss_wt1;
	print "\n ss_wt2            $ss_wt2           " if $ss_wt2;
	print "\n flag_no_ss_wt2    $flag_no_ss_wt2   " if $flag_no_ss_wt2;
	print "\n cpu               $cpu              " if $cpu;
	print "\n flag_run_betapro  $flag_run_betapro " if $flag_run_betapro;
	print "\n file_betapro      $file_betapro     " if $file_betapro;
	print "\n file_pairing      $file_pairing     " if $file_pairing;
	print "\n atom_scheme       $atom_select_scheme";
	print "\n";
	print "\nFixed Parameters:";
	print "\n model_count       $model_count      " if $model_count;
	print "\n rep1              $rep1             " if $rep1;
	print "\n rep2              $rep2             " if $rep2;
	print "\n mini              $mini             " if $mini;
	print "\n bpro_thres        $bpro_thres       " if $bpro_thres;
	print "\n pair_caca_thres   $pair_caca_thres  " if $pair_caca_thres;
	print "\n pair_hbond_thres  $pair_hbond_thres " if $pair_hbond_thres;
	print "\n pair_policy       $pair_policy      " if $pair_policy; 
	print "\n";
}

sub load_global_variables{
	if($flag_flush){
		system_cmd("rm -rf $dir_job/*.log");
		system_cmd("rm -rf $dir_job/logs");
		system_cmd("rm -rf $dir_job/stage*");
		system_cmd("rm -rf $dir_job/extended");
		system_cmd("rm -rf $dir_job/sheet_detect");
		system_cmd("rm -rf $dir_job/con_filter");
		system_cmd("rm -rf $dir_job/betapro");
	}
	# load directories names
	$dir_cns_install = "/var/www/cgi-bin/confold/bin/cns_solve_1.3";
	$dir_extended = "$dir_job/extended";
	$dir_log = "$dir_job/logs";
	mkdir $dir_log or die $! if not -d $dir_log;
	$dir_scripts = "".dirname(abs_path($0));
	$program_betapro1 = "/rose/space1/tools_badri/tools/betapro/betapro-1.0/bin/predict_beta_ss.sh";
	$program_betapro2 = "/rose/space1/tools_badri/tools/betapro/betapro-1.0/bin/predict_strand.sh";
	$dir_betapro = "$dir_job/betapro" if $flag_run_betapro;
	mkdir $dir_betapro or die $! if $flag_run_betapro and not -d $dir_betapro;
	# load required files
	my @seq_files = <$dir_job/*.fasta>;
	confess "ERROR! fasta file not found in $dir_job!" if not defined $seq_files[0];
	confess "ERROR! multiple fasta files in $dir_job!" if defined $seq_files[1];
	$file_seq = $seq_files[0];
	$seq_id = basename($file_seq, ".fasta");
	flatten_fasta($file_seq);
	$seq_len = length(seq_fasta($file_seq));
	my @sec_files = <$dir_job/*.ss>;
	confess "ERROR! SS file not found in $dir_job!" if not defined $sec_files[0];
	confess "ERROR! multiple sec files in $dir_job!" if defined $sec_files[1];
	$file_sec = $sec_files[0];
	confess "ERROR! fasta and ss file name must match!" if (basename($file_seq, ".fasta") ne basename($file_sec, ".ss"));
	confess "ERROR! fasta sequence and and ss sequence do not match!" if (seq_fasta($file_seq) ne seq_ss($file_sec));
	if (!$flag_no_rr){
		my @rr_files = <$dir_job/*.rr>;
		confess "ERROR! RR file not found in $dir_job!" if not defined $rr_files[0];
		confess "ERROR! multiple rr files in $dir_job!" if defined $rr_files[1];
		$file_rr = $rr_files[0];
		confess "ERROR! fasta sequence and and rr sequence do not match!" if (seq_fasta($file_seq) ne seq_rr($file_rr));
	}
	%residues = fasta2residues_hash($file_seq);
	# load stage information
	$stage_list{"stage1"} = 1;
	$stage_list{"stage2"} = 2 if $flag_stage2;
	$stage_list{"sheet_detect"} = 2 if $flag_sheet_detect;
	$stage_list{"con_filter"} = 2 if $flag_con_filter;
	# load global variables for no contacts
	if ($flag_no_rr){
		$sec_wt_list{$ss_wt1} = 1;
		$xL{"none"} = 1;
	}
	# for true contacts
	elsif ($flag_all_rr){
		$sec_wt_list{$ss_wt1} = 1;
		$sec_wt_list{$ss_wt2} = 2 if $ss_wt2;
		$xL{"all"} = 1;
	}
	# for predicted contacts
	else{
		$sec_wt_list{$ss_wt1} = 1;
		$sec_wt_list{$ss_wt2} = 2 if $ss_wt2;
		%xL = qw( 0.4 1 0.6 1 0.8 1 1.0 1 1.2 1 1.4 1 1.6 1 1.8 1 2.0 1 2.2 1 );
	}
}

sub print_global_variables{
	print "\nGlobal Variables:";
	print "\n file_seq           $file_seq         ";
	print "\n file_sec           $file_sec         ";
	print "\n file_rr            $file_rr          " if $file_rr;
	print "\n seq_id             $seq_id           ";
	print "\n seq_len            $seq_len          ";
	print "\n dir_cns_install    $dir_cns_install  ";
	print "\n dir_log            $dir_log          ";
	print "\n dir_extended       $dir_extended     ";
	print "\n stages             ";
	foreach (keys %stage_list){
		print $_." ";
	}
	print "\n secondary weights  ";
	foreach (keys %sec_wt_list){
		print $_." ";
	}
	if ($flag_no_rr){
		print "\n";
		return;
	}
	print "\n contact selection  ";
	foreach (sort keys %xL){
		print $_." ";
	}
	my %all_contacts = rr2contacts_hash($file_rr);
	my %standard_contacts = rr2contacts_hash($file_rr, 1);
	print "\n total input rr     ".(scalar keys %all_contacts);
	if ((scalar keys %all_contacts) ne (scalar keys %standard_contacts)){
		print "\n standard contacts  ".(scalar keys %standard_contacts);
		warn "\nWARNING! some input contacts are below 6 sequence separation!"
	}
	print "\n";
}

sub build_extended{
	my $file_seq = "$seq_id.fasta";
	confess "ERROR! expected fasta file not found in the extended directory!" if not -f $file_seq;
	flatten_fasta($file_seq);
	write_cns_seq($file_seq, "input.seq");
	system_cmd("cp $dir_cns_install/inputs/general/generate_seq.inp ./generate_seq.inp");
	system_cmd("chmod +rw ./generate_seq.inp");
	system_cmd("sed -i s/amy_start.seq/input.seq/g ./generate_seq.inp");
	system_cmd("sed -i s/generate_seq.mtf/extended.mtf/ ./generate_seq.inp");
	# Changing hydrogen_flag to true/false does not work. A bug?
	system_cmd("sed -i s/hydrogen_flag=false/hydrogen_flag=true/g generate_seq.inp");
	# So, I remove it completely.
	system_cmd("sed -i \'s/delete selection=( hydrogen ) end//g\' generate_seq.inp");
	system_cmd("sed -i s/CNS_TOPPAR:protein_rep.param/CNS_TOPPAR:protein.param/g generate_seq.inp");
	system_cmd("cp --no-preserve=mode $dir_cns_install/inputs/nmr_calc/generate_extended.inp ./generate_extended.inp");
	system_cmd("sed -i s/generate_extended.mtf/extended.mtf/g ./generate_extended.inp");
	system_cmd("sed -i s/amy_extended.pdb/extended.pdb/g ./generate_extended.inp");
	system_cmd("sed -i s/CNS_TOPPAR:protein-allhdg5-4.param/CNS_TOPPAR:protein.param/g ./generate_extended.inp");
	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                     \n";
	print JOB "# CNS-CONFIGURATION                             \n";
	print JOB "source $dir_cns_install/cns_solve_env.sh        \n";
	print JOB "export KMP_AFFINITY=none                        \n";
	print JOB "$dir_cns_install/intel-x86_64bit-linux/bin/cns_solve < generate_seq.inp > generate_seq.log \n";
	print JOB "$dir_cns_install/intel-x86_64bit-linux/bin/cns_solve < generate_extended.inp > generate_extended.log \n";
	close JOB;
	system_cmd("chmod +x job.sh");
	system_cmd("./job.sh", "job.log");
	confess "FAILED! extended.mtf not found!" if not -f "extended.pdb";
	system_cmd("rm -f generate_seq.log");
}

sub contact_restraints{
	return if $flag_no_rr;
	my $stage1_rr = "$dir_job/stage1/sec_${current_sec_wt}_${current_topx}L/contact.rr";
	my $stage1_top_model = "$dir_job/stage1/sec_${current_sec_wt}_${current_topx}L/${seq_id}_model1.pdb";
	if($current_stage eq "stage1"){
		if($flag_all_rr){
			system_cmd("cp $file_rr contact.rr");
		}
		elsif($flag_select_rr){
			my $top_xl = int(($current_topx * $seq_len) + 0.5) + 1; # +1 to account for header line
			system_cmd("head -n $top_xl $file_rr > contact.rr");
		}
		else{
			confess "ERROR! neither all rr nor select rr!";
		}
	}
	elsif($current_stage eq "sheet_detect"){
		system_cmd("cp $stage1_rr ./");
	}
	elsif($current_stage eq "con_filter" or $current_stage eq "stage2"){
		rr_remove_unsatisfied_top_model($stage1_rr, $stage1_top_model, "contact.rr", "rr_filter.log");
	}
	else{
		confess ":( Unexpected";
	}
	cbrr2tbl("contact.rr", "contact.tbl") if $rr_type eq "cb";
	carr2tbl("contact.rr", "contact.tbl") if $rr_type eq "ca";
}

sub build_models_in_background{
	my $dihed_wt1 = $con_wt * $current_sec_wt;
	my $dihed_wt2 = $con_wt * $current_sec_wt;
	my ($tbl_contact, $ss_horiz, $tbl_ssnoe, $tbl_hbond, $tbl_dihed);
	confess "ERROR! expected fasta file not found!" if not -f $file_seq;
	my $secondary = "$seq_id.ss";
	confess "ERROR! expected secondary structure file not found!" if not -f $secondary;
	# load restraints and print to logs
	my %tbl_list = ();
	$tbl_list{"contact"}  = "contact.tbl"  if (-f "contact.tbl"  and count_lines("contact.tbl") > 0);
	$tbl_list{"ssnoe"}    = "ssnoe.tbl"    if (-f "ssnoe.tbl"    and count_lines("ssnoe.tbl") > 0);
	$tbl_list{"hbond"}    = "hbond.tbl"    if (-f "hbond.tbl"    and count_lines("hbond.tbl") > 0);
	$tbl_list{"dihedral"} = "dihedral.tbl" if (-f "dihedral.tbl" and count_lines("dihedral.tbl") > 0);
	warn "\n\nWarning! contact.tbl not found! building secondary structure models!" if not defined $tbl_list{"contact"};
	warn "\nWarning! dihedral.tbl not found! build models with no secondary structures!\n" if not defined $tbl_list{"dihedral"};
	rr2chimera_ca($tbl_list{"contact"}, "contact.chimera", "ca") if defined $tbl_list{"contact"};
	print "\n".seq_fasta("${seq_id}.fasta");
	print "\n".ss_line_ss($current_file_sec);
	foreach my $tbl (sort keys %tbl_list){
		my $flag_dihe = 0;
		$flag_dihe = 1 if $tbl eq "dihedral";
		print "\n".coverage_tbl("${seq_id}.fasta", $tbl_list{$tbl}, $flag_dihe);
	}
	# prepare CNS task files
	system_cmd("rm -f iam.*");
	system_cmd("cp $dir_cns_install/inputs/nmr_calc/dg_sa.inp ./");
	system_cmd("chmod +rw ./dg_sa.inp");
	system_cmd("cp $dir_scripts/scalecoolsetupedited ./");
	system_cmd("cp $dir_scripts/scalehotedited ./");
	system_cmd("sed -i \"s/input.noe.scale \\* 10/input.noe.scale \\* $current_sec_wt/\" scalehotedited");
	system_cmd("sed -i \"s/input.noe.scale \\* 10/input.noe.scale \\* $current_sec_wt/\" scalecoolsetupedited");
	system_cmd("sed -i s/CNS_TOPPAR:protein-allhdg5-4.param/CNS_TOPPAR:protein.param/g dg_sa.inp");
	system_cmd("sed -i s/CNS_NMRMODULE:scalehot/CNS_CUSTOMMODULE:scalehotedited/g dg_sa.inp");
	system_cmd("sed -i s/CNS_NMRMODULE:scalecoolsetup/CNS_CUSTOMMODULE:scalecoolsetupedited/g dg_sa.inp");
	system_cmd("sed -i s/il8.mtf/extended.mtf/ dg_sa.inp");
	system_cmd("sed -i s/il8.pdb/extended.pdb/ dg_sa.inp");
	system_cmd("sed -i s/pdb.dg.count=10/pdb.dg.count=${model_count}/ dg_sa.inp");
	system_cmd("sed -i s/pdb.end.count=10/pdb.end.count=${model_count}/ dg_sa.inp");
	system_cmd("sed -i s/flg.calc.ave.struct=true/flg.calc.ave.struct=false/ dg_sa.inp");
	system_cmd("sed -i s/flg.calc.ave.accpt=true/flg.calc.ave.accpt=false/ dg_sa.inp");
	if ($atom_select_scheme == 1){
		# leave default
	}
	elsif ($atom_select_scheme == 2){
		# just add oxygen
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name c or name cb\\* or name cg\\* or name o/\" dg_sa.inp");
	}
	elsif ($atom_select_scheme == 3){
		# just add oxygen and hydrogen
		system_cmd("sed -i \"s/name ca or name ha or name n or name hn/name ca or name ha or name n or name hn or name h/\" dg_sa.inp");
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name c or name cb\\* or name cg\\* or name o/\" dg_sa.inp");
	}
	elsif ($atom_select_scheme == 4){
		# backbone atoms only
		system_cmd("sed -i \"s/name ca or name ha or name n or name hn/name c or name ca/\" dg_sa.inp");
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name n or name o/\" dg_sa.inp");
	}
	elsif ($atom_select_scheme == 5){
		# backbone atoms with cb
		system_cmd("sed -i \"s/name ca or name ha or name n or name hn/name c or name ca or name cb/\" dg_sa.inp");
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name n or name o/\" dg_sa.inp");
	}
	elsif ($atom_select_scheme == 6){
		# backbone atoms with cb and hydrogen
		# atom selection according to instructions in the NIH-XPLORE manual
		system_cmd("sed -i \"s/name ca or name ha or name n or name hn/name c or name ca or name cb/\" dg_sa.inp");
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name n or name o or name h/\" dg_sa.inp");
	}
	elsif ($atom_select_scheme == 7){
		# backbone atoms with cb, hydrogen and cg
		system_cmd("sed -i \"s/name ca or name ha or name n or name hn/name c or name ca or name cb/\" dg_sa.inp");
		system_cmd("sed -i \"s/or name c or name cb\\* or name cg\\*/or name n or name o or name h or name cg\\*/\" dg_sa.inp");
	}
	else{
		confess "ERROR! Invalid atom_select_scheme $atom_select_scheme!";
	}
	system_cmd("sed -i s/md.cool.init.rad=0.9/md.cool.init.rad=$rep1/ dg_sa.inp");
	system_cmd("sed -i s/md.cool.fina.rad=0.8/md.cool.fina.rad=$rep2/ dg_sa.inp");
	system_cmd("sed -i s/md.cool.noe=50/md.cool.noe=$con_wt/ dg_sa.inp");
	system_cmd("sed -i s/md.pow.noe=50/md.pow.noe=$con_wt/ dg_sa.inp");
	system_cmd("sed -i s/md.cool.cdih=200/md.cool.cdih=$dihed_wt1/ dg_sa.inp");
	system_cmd("sed -i s/md.pow.cdih=400/md.pow.cdih=$dihed_wt2/ dg_sa.inp");
	system_cmd("sed -i s/md.pow.step=200/md.pow.step=$mini/ dg_sa.inp");
	system_cmd("sed -i s/il8_noe.tbl//g dg_sa.inp")                              if not defined $tbl_list{"contact"};
	system_cmd("sed -i s/il8_hbonds.tbl//g dg_sa.inp")                           if not defined $tbl_list{"hbond"};
	system_cmd("sed -i s/il8_dihe.tbl//g dg_sa.inp")                             if not defined $tbl_list{"dihedral"};
	system_cmd("sed -i s/il8_noe.tbl/contact.tbl/g dg_sa.inp")                   if defined $tbl_list{"contact"};
	system_cmd("sed -i s/nmr.noe.file.2=\\\"\\\"/nmr.noe.file.2=\\\"ssnoe.tbl\\\"/g dg_sa.inp") if defined $tbl_list{"ssnoe"};
	system_cmd("sed -i s/il8_hbonds.tbl/hbond.tbl/g dg_sa.inp")                  if defined $tbl_list{"hbond"};
	system_cmd("sed -i s/il8_dihe.tbl/dihedral.tbl/g dg_sa.inp")                 if defined $tbl_list{"dihedral"};
	system_cmd("sed -i s/pdb.out.name=\\\"dg\\\"/pdb.out.name=\\\"${seq_id}\\\"/g dg_sa.inp");

	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                       \n";
	print JOB "echo \"starting cns..\"                           \n";
	print JOB "touch iam.running                                 \n";
	print JOB "# CNS-CONFIGURATION                               \n";
	print JOB "source $dir_cns_install/cns_solve_env.sh          \n";
	print JOB "export KMP_AFFINITY=none                          \n";
	print JOB "export CNS_CUSTOMMODULE=$current_job_dir          \n";
	print JOB "$dir_cns_install/intel-x86_64bit-linux/bin/cns_solve < dg_sa.inp > dg_sa.log \n";
	print JOB "if [ -f \"${seq_id}_${model_count}.pdb\" ]; then  \n";
	print JOB "   echo \"trial structures written.\"             \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   echo \"running assess.pl to rank models..\"    \n";
	print JOB "   $dir_scripts/assess_dgsa.pl ./ $seq_id &> ./assess.log \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	print JOB "if [ -f \"${seq_id}a_${model_count}.pdb\" ]; then \n";
	print JOB "   echo \"accepted structures written.\"          \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   echo \"running assess.pl to rank models..\"    \n";
	print JOB "   $dir_scripts/assess_dgsa.pl ./ $seq_id &> ./assess.log \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	print JOB "tail -n 30 dg_sa.log                              \n";
	print JOB "echo \"ERROR! Final structures not found!\"       \n";
	print JOB "echo \"CNS FAILED!\"                              \n";
	print JOB "mv iam.running iam.failed                         \n";
	close JOB;
	print "\n"."starting $current_job_dir/job.sh";
	system_cmd("chmod +x $current_job_dir/job.sh");
	system_cmd("$current_job_dir/job.sh &> log_job.txt &");
}

sub sec_restraints{
	chdir $current_job_dir or die $!;
	my $stage1_top_model = "$dir_job/stage1/sec_${current_sec_wt}_${current_topx}L/${seq_id}_model1.pdb" if $current_stage ne "stage1";
	system_cmd("rm -f secondary_restraints.log");
	print_log("secondary_restraints.log", seq_fasta($file_seq)." [$file_seq]");
	print_log("secondary_restraints.log", ss_line_ss($current_file_sec)." [$current_file_sec]");
	confess "\nSequence mismatch! ".seq_ss($current_file_sec)." [$current_file_sec]" if (seq_fasta($file_seq) ne seq_ss($current_file_sec));
	# helix restraints; start writing to the files hbond.tbl, ssnoe.tbl and dihedral.tbl
	print_helix_restraints();
	# strand and/or sheet restraints
	if (strand_count($current_file_sec) == 0){
		print_log("secondary_restraints.log", "no strand restraints");
		return;
	}
	if (strand_count($current_file_sec) == 1){
		warn "\nWARNING! only 1 strand! please check your secondary structure file!";
		print_log("secondary_restraints.log", "WARNING! only 1 strand! please check your secondary structure file!");
		return;
	}
	if($current_file_pairing){
		validate_pairing_file($current_file_pairing);
		write_pairing_hbonds("pairing.hbonds");
	}
	elsif($current_stage eq "stage2" or $current_stage eq "sheet_detect"){
		detect_hbonds_from_model($stage1_top_model, "pairing.hbonds");
		$current_file_pairing = "pairing.txt" if -s "pairing.txt";
		print "\nWARNING! pair detection did not write any hbonds file!" if (not -s "pairing.hbonds") or (not -s $current_file_pairing);
	}
	if($current_file_pairing){
		update_sec_using_pairing_info($current_file_sec, $current_file_pairing, "updated.ss");
		system_cmd("mv $current_file_sec $current_file_sec.old");
		system_cmd("mv updated.ss $current_file_sec");
	}
	# continue to write to hbond.tbl, ssnoe.tbl and dihedral.tbl
	strand_and_sheet_tbl("pairing.hbonds") if $current_file_pairing;
	strand_and_sheet_tbl() if !$current_file_pairing;
}

sub print_helix_restraints{
	my %res_sec = each_residue_sec_ss($current_file_sec);
	my %all_residues = %res_sec;
	foreach (keys %res_sec){
		delete $res_sec{$_} if $res_sec{$_} ne "H";
	}
	if (scalar keys %res_sec){
		print_log("secondary_restraints.log", "write helix tbl restraints");
		open DIH, ">dihedral.tbl" or confess $!;
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			my @PHI = split /\s+/, $restraints_dihedral{"H PHI"};
			my @PSI = split /\s+/, $restraints_dihedral{"H PSI"};
			# 9/1/2014: I found that if I remove phi angle from the pre-first residue, the first residue cannot form helix most of the times. I checked 4/5 proteins and 2 dozen helices.
			# So, adding the pre-first phi angle as well (by removing the condition for  res_sec{$i-1} to exist)
			printf DIH "assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0 %7s %7s 2 !helix phi\n", $i-1, $i, $i, $i, $PHI[0], $PHI[1] if defined $all_residues{$i+1};
			# 9/1/2014: Even if we don't have PSI for the last residue, I can see that last residue form helix, almost in all cases.
			# 9/2/2014: I can see that in some cases like target Tc767, only 106 helix residues were formed while the input is 111. And all of it is because of the psi angle at the end.
			# So, adding the post-last psi angle as well (by removing the condition for res_sec{$i+1} to exist)
			printf DIH "assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0 %7s %7s 2 !helix psi\n", $i, $i, $i, $i+1, $PSI[0], $PSI[1] if defined $all_residues{$i+1};
		}
		close DIH;
		open HBO, ">hbond.tbl" or confess $!;
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			next if not defined $res_sec{$i+1};
			next if not defined $res_sec{$i+2};
			next if not defined $res_sec{$i+3};
			next if not defined $res_sec{$i+4};
			my @HR = split /\s+/, $restraints_hbond{"H"};
			# In case of Proline, use N atom instead of H
			my $HATOM = "H";
			$HR[0] -= 1.0;
			if($residues{$i+4} eq "P"){
				$HR[0] += 1.0;
				$HATOM  = "N";
			}
			printf HBO "assign (resid %3d and name O) (resid %3d and name $HATOM) %4s %4s %4s !helix\n", $i, $i+4, $HR[0], $HR[1], $HR[2];
		}
		close HBO;
		open NOE, ">ssnoe.tbl" or confess $!;
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			next if not defined $res_sec{$i+1};
			next if not defined $res_sec{$i+2};
			next if not defined $res_sec{$i+3};
			next if not defined $res_sec{$i+4};
			foreach my $A (sort keys %ATOMTYPE){
				foreach my $SH (sort keys %SHIFT){
					my @AR = split /\s+/, $restraints_distance{"H $A-$A O $SH"};
					# found restraints even for non-helix SS
					next if not defined $res_sec{$i+4+$SH};
					confess ":(" if not defined ($i and $A and ($i+4+$SH) and $AR[0] and $AR[1] and $AR[2]);
					printf NOE "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f !helix\n", $i, $A, ($i+4+$SH), $A, $AR[0], $AR[1], $AR[2];
				}
			}
		}
		close NOE;
	}
	else{
		print_log("secondary_restraints.log", "no helix predictions!");
	}
}

sub write_pairing_hbonds{
	my $file_hbonds = shift;
	my %res_ss = each_residue_sec_ss($current_file_sec);
	my %pairs = pairing_info2hash($current_file_pairing);
	print_log("secondary_restraints.log", "initial pairing");
	print_log("secondary_restraints.log", "A/P  R1-R2   R1-R2  Confidence");
	foreach (sort {$pairs{$a} <=> $pairs{$b}} keys %pairs){
		my @A = split /\s+/, $_;
		print_log("secondary_restraints.log", sprintf "%3s %3d-%-3d %3d-%-3d %.2f", $A[4], $A[0], $A[1], $A[2], $A[3], $A[5]);
	}
	#build hydrogen bonding pattern, and shift strands if necessary, for pairing
	my %hbond_connectors = ();
	my %hbonds = ();
	my %slided_pairs = ();
	my %possible_slides = qw( 0 1 L1 2 R1 3 L2 4 R2 5 ); # total 10 possible configurations
	my %possible_configs = qw( a 1 b 2 c 3 d 4 ); # parallel can have c and d hydrogen bonding pattern as well
	print_log("secondary_restraints.log", "Selected Configurations:");
	foreach my $pair (sort {$pairs{$a} <=> $pairs{$b}} keys %pairs){
		my @A = split /\s+/, $pair;
		if ((abs($A[0] - $A[1]) < 2) or abs($A[2] - $A[3]) < 2){
			delete $pairs{$pair};
			next;
		}
		my $selected_shift_config = undef;
		foreach my $slide (sort {$possible_slides{$a} <=> $possible_slides{$b}} keys %possible_slides){
			last if defined $selected_shift_config;
			foreach my $config (sort {$possible_configs{$a} <=> $possible_configs{$b}} keys %possible_configs){
				last if defined $selected_shift_config;
				my @A = split /\s+/, $pair;
				next if ($A[4] eq "A" and not ($config eq "a" or $config eq "b")); # anti-parallel do not have hbond-patterns c and d
				my $no_of_residues = $slide;
				$no_of_residues =~ s/[A-Z]//g;
				# first strand is sliding forward 
				if ($slide =~ /^R/){
					$A[1] = $A[1] - $no_of_residues if $A[4] eq "A";
					$A[2] = $A[2] - $no_of_residues if $A[4] eq "A";
					$A[1] = $A[1] - $no_of_residues if $A[4] eq "P";
					$A[2] = $A[2] + $no_of_residues if $A[4] eq "P";
				}
				# first strand is sliding backward
				elsif ($slide =~ /^L/){
					$A[0] = $A[0] + $no_of_residues if $A[4] eq "A";
					$A[3] = $A[3] + $no_of_residues if $A[4] eq "A";
					$A[0] = $A[0] + $no_of_residues if $A[4] eq "P";
					$A[3] = $A[3] - $no_of_residues if $A[4] eq "P";
				}
				next if ((abs($A[0] - $A[1]) < 2) or abs($A[2] - $A[3]) < 2);
				warn "\nWARNING! 2 residue shifting is needed for pair $pair slide $slide config $config!" if $no_of_residues == 2;
				my %ideal_hbonds = make_ideal_hbonds($A[0], $A[1], $A[2], $A[3], $A[4], $config);
				my %connectors = hbonds2connectors(\%ideal_hbonds);
				my $flag_reject = 0;
				foreach (keys %connectors){
					$flag_reject = 1 if (defined $hbond_connectors{$_});
				}
				if (not $flag_reject){
					$selected_shift_config = $A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." ".$config;
					print_log("secondary_restraints.log", "pair:$pair slide:$slide selected:".$A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." ".$config);
				}
			}
		}
		confess ":( $_" if not defined $selected_shift_config;
		my @S = split /\s+/, $selected_shift_config;
		my %ideal_hbonds = make_ideal_hbonds($S[0], $S[1], $S[2], $S[3], $S[4], $S[5]);
		my %connectors = hbonds2connectors(\%ideal_hbonds);
		foreach (keys %ideal_hbonds){
			$hbonds{$_} = $S[4];
		}
		foreach (keys %connectors){
			$hbond_connectors{$_} = 1;
		}
		$slided_pairs{$S[0]." ".$S[1]." ".$S[2]." ".$S[3]." ".$S[4]." ".$S[5]} = $S[0];
	}
	open HBOND, ">$file_hbonds" or confess $!;
	foreach (sort keys %hbonds){
		print HBOND $_." ".$hbonds{$_}."\n";
	}
	close HBOND;
}

sub strand_and_sheet_tbl{
	my $file_pairing_hbond = shift;
	my %res_ss = each_residue_sec_ss($current_file_sec);
	my %residues = ss2residues_hash($current_file_sec);
	my %unpaired_residues = %res_ss;
	my %paired_residues = ();
	my %res_ssE = %res_ss;
	foreach (sort keys %res_ssE){
		delete $res_ssE{$_} if $res_ssE{$_} ne "E";
	}
	my %hbonds = ();
	if ($current_file_pairing){
		%paired_residues = paired_residues($current_file_pairing);
		foreach (sort keys %unpaired_residues){
			delete $unpaired_residues{$_} if $unpaired_residues{$_} ne "E";
		}
		foreach (sort keys %unpaired_residues){
			delete $unpaired_residues{$_} if defined $paired_residues{$_};
		}
		%hbonds = load_hbonds($file_pairing_hbond);
		confess ":(" if not scalar %paired_residues;
		confess ":(" if not scalar %hbonds;

		open HBOND, ">>hbond.tbl" or confess $!;
		foreach (sort keys %hbonds){
			my @H = split /\s+/, $_;
			my @HR = split /\s+/, $restraints_hbond{$hbonds{$_}};
			confess ":( distance not defined ".$hbonds{$_} if not (defined $HR[0] and defined $HR[1] and defined $HR[2]);
			# In case of Proline, use N atom instead of H
			$HR[0] -= 1.0;
			if($residues{$H[0]} eq "P" and $H[1] eq "H"){
				$HR[0] += 1.0;
				$H[1]  = "N";
			}
			if($residues{$H[2]} eq "P" and $H[3] eq "H"){
				$HR[0] += 1.0;
				$H[3]  = "N";
			}
			printf HBOND "assign (resid %3d and name %1s) (resid %3d and name %1s) %4s %4s %4s !sheet\n", $H[0], $H[1], $H[2], $H[3], $HR[0], $HR[1], $HR[2];
		}
		close HBOND;
		open SSNOE, ">>ssnoe.tbl" or confess $!;
		foreach my $hb (sort keys %hbonds){
			my @H = split /\s+/, $hb;
			my @hbondConnectors = ( $H[0]." ".$H[1]." ".$H[2], $H[2]." ".$H[3]." ".$H[0]); 
			foreach (@hbondConnectors){
				my @HBC = split /\s+/, $_;
				foreach my $A (sort keys %ATOMTYPE){
					foreach my $S (sort keys %SHIFT){
						# A   O-O   O -1 
						my $distances = $restraints_distance{$hbonds{$hb}." ".$A."-".$A." ".$HBC[1]." ".$S};
						my @DIST = split /\s+/, $distances;
						confess ":( distance not defined ".$hbonds{$hb}." ".$A."-".$A." ".$HBC[1]." ".$S if not defined $DIST[2];
						next if not defined $paired_residues{($HBC[2]+$S)};
						# observing decrease in quality of many sheet detected proteins, I am adding this
						# next if not defined $paired_atoms{($HBC[2]+$S)." ".$A};
						printf SSNOE "assign (resid %3d and name %2s) (resid %3d and name %2s) %2.2f %2.2f %2.2f !sheet\n", $HBC[0], $A, ($HBC[2]+$S), $A, $DIST[0], $DIST[1], $DIST[2];
					}
				}
			}
		}
		close SSNOE;
	}

	# Identify the strands that are not used for pairing and generate generic dihedral restraints for them
	open DIH, ">>dihedral.tbl" or confess $!;
	foreach my $i (sort {$a <=> $b} keys %res_ssE){
		my @SPHI = ();
		my @SPSI = ();
		my $strand_type = "unpaired E residue";
		if (defined $paired_residues{$i}){
			@SPHI = split /\s+/, $restraints_dihedral{$paired_residues{$i}." PHI"};
			@SPSI = split /\s+/, $restraints_dihedral{$paired_residues{$i}." PSI"};
			$strand_type = "paired E residue";
		}
		else{
			@SPHI = split /\s+/, $restraints_dihedral{"U PHI"};
			@SPSI = split /\s+/, $restraints_dihedral{"U PSI"};
		}
		if (defined $res_ss{$i-1} and $res_ss{$i-1} eq "E"){
			printf DIH "assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0 %7s %7s 2 !$strand_type phi\n", $i-1, $i, $i, $i, $SPHI[0], $SPHI[1];
		}
		if (defined $res_ss{$i+1} and $res_ss{$i+1} eq "E"){
			printf DIH "assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0 %7s %7s 2 !$strand_type psi\n", $i, $i, $i, $i+1, $SPSI[0], $SPSI[1];
		}
	}
	close DIH;

	# Identify the strands that are not used for pairing and generate generic restraints for them
	open NOE, ">>ssnoe.tbl" or confess $!;
	foreach my $i (sort {$a <=> $b} keys %res_ssE){
		my @SD = ();
		my $strand_type = "unpaired E residue";
		if (defined $paired_residues{$i}){
			@SD = split /\s+/, $restraints_strandOO{$paired_residues{$i}};
			confess ":(" if (!$SD[0] or !$SD[1] or !$SD[2]);
			$strand_type = "paired E residue";
		}
		else{
			@SD = split /\s+/, $restraints_strandOO{"U"};
			confess ":(" if (!$SD[0] or !$SD[1] or !$SD[2]);
		}
		next if not defined $res_ssE{$i+1};
		next if $res_ssE{$i+1} ne "E";
		printf NOE "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f !$strand_type\n", $i, "O", $i+1, "O", $SD[0], $SD[1], $SD[2];
	}
	close NOE;
}

sub load_ss_restraints{
	my $lambda = shift;
	my $log_reference = shift;
	
	# T      Helix or Parallel or anti-parallel or Unknown Strand Type
	# A1_A2  Atom1-Atom2 pair
	# Ref    Hydrogen bonding connector atom (reference hbond connector)
	# N      Neighborhood residue shifting on the hbond connector of R2 side. For example, If R1:N and R2:O have hbond and S = +1, the restraint A1-A2 are for R1 and (R2+1)
	# Note: hbond distances are the distances between Nitrogen and Oxygen
	# Places to verify:
	# http://www.beta-sheet.org/page29/page51/page53/index.html
	# In this model, the HP sheet is composed of identical straight helical chains with phi = -122 degrees, psi = 135 degrees, and a slightly non-linear interchain H-bond angle delta of 170 degrees.
	# http://en.wikipedia.org/wiki/Alpha_helix
	# Residues in α-helices typically adopt backbone (φ, ψ) dihedral angles around (-60°, -45°), as shown in the image at right
	# the H to O distance is about 2 Å (0.20 nm

	# T A Mean Standard_deviation
	$restraints_dihedral{"A PSI"} = "136.91 ".($lambda * 17.39);
	$restraints_dihedral{"A PHI"} = "-120.89 ".($lambda * 21.98);
	$restraints_dihedral{"P PSI"} = "130.96 ".($lambda * 16.66);
	$restraints_dihedral{"P PHI"} = "-115.00 ".($lambda * 20.31);
	$restraints_dihedral{"U PSI"} = "134.95 ".($lambda * 17.65);
	$restraints_dihedral{"U PHI"} = "-118.91 ".($lambda * 21.73);
	$restraints_dihedral{"H PSI"} = "-41.51 ".($lambda * 9.84);
	$restraints_dihedral{"H PHI"} = "-63.47 ".($lambda * 9.20);

	#T A1_A2 Ref N Mean Standard_deviation
	$restraints_distance{"A O-O O +1"} = get_dist_neg_pos(7.73, 0.59, $lambda);
	$restraints_distance{"A O-O O -1"} = get_dist_neg_pos(4.84, 0.16, $lambda);
	$restraints_distance{"A O-O O 0"} = get_dist_neg_pos(3.57, 0.28, $lambda);
	$restraints_distance{"A O-O H +1"} = get_dist_neg_pos(7.76, 0.60, $lambda);
	$restraints_distance{"A O-O H -1"} = get_dist_neg_pos(4.90, 0.45, $lambda);
	$restraints_distance{"A O-O H 0"} = get_dist_neg_pos(3.58, 0.31, $lambda);
	$restraints_distance{"A C-C O +1"} = get_dist_neg_pos(7.66, 0.52, $lambda);
	$restraints_distance{"A C-C O -1"} = get_dist_neg_pos(4.80, 0.17, $lambda);
	$restraints_distance{"A C-C O 0"} = get_dist_neg_pos(4.96, 0.21, $lambda);
	$restraints_distance{"A C-C H +1"} = get_dist_neg_pos(7.65, 0.51, $lambda);
	$restraints_distance{"A C-C H -1"} = get_dist_neg_pos(4.85, 0.34, $lambda);
	$restraints_distance{"A C-C H 0"} = get_dist_neg_pos(4.96, 0.21, $lambda);
	$restraints_distance{"A N-N O +1"} = get_dist_neg_pos(5.09, 0.34, $lambda);
	$restraints_distance{"A N-N O -1"} = get_dist_neg_pos(6.86, 0.40, $lambda);
	$restraints_distance{"A N-N O 0"} = get_dist_neg_pos(4.42, 0.24, $lambda);
	$restraints_distance{"A N-N H +1"} = get_dist_neg_pos(5.04, 0.21, $lambda);
	$restraints_distance{"A N-N H -1"} = get_dist_neg_pos(6.85, 0.45, $lambda);
	$restraints_distance{"A N-N H 0"} = get_dist_neg_pos(4.43, 0.25, $lambda);
	$restraints_distance{"A CA-CA O +1"} = get_dist_neg_pos(6.43, 0.41, $lambda);
	$restraints_distance{"A CA-CA O -1"} = get_dist_neg_pos(5.67, 0.28, $lambda);
	$restraints_distance{"A CA-CA O 0"} = get_dist_neg_pos(5.26, 0.24, $lambda);
	$restraints_distance{"A CA-CA H +1"} = get_dist_neg_pos(6.38, 0.36, $lambda);
	$restraints_distance{"A CA-CA H -1"} = get_dist_neg_pos(5.71, 0.40, $lambda);
	$restraints_distance{"A CA-CA H 0"} = get_dist_neg_pos(5.27, 0.25, $lambda);
	$restraints_distance{"P O-O O +1"} = get_dist_neg_pos(7.90, 0.61, $lambda);
	$restraints_distance{"P O-O O -1"} = get_dist_neg_pos(4.86, 0.16, $lambda);
	$restraints_distance{"P O-O O 0"} = get_dist_neg_pos(3.78, 0.34, $lambda);
	$restraints_distance{"P O-O H +1"} = get_dist_neg_pos(4.92, 0.40, $lambda);
	$restraints_distance{"P O-O H -1"} = get_dist_neg_pos(8.02, 0.60, $lambda);
	$restraints_distance{"P O-O H 0"} = get_dist_neg_pos(3.78, 0.32, $lambda);
	$restraints_distance{"P C-C O +1"} = get_dist_neg_pos(8.03, 0.51, $lambda);
	$restraints_distance{"P C-C O -1"} = get_dist_neg_pos(4.82, 0.17, $lambda);
	$restraints_distance{"P C-C O 0"} = get_dist_neg_pos(5.21, 0.25, $lambda);
	$restraints_distance{"P C-C H +1"} = get_dist_neg_pos(4.88, 0.34, $lambda);
	$restraints_distance{"P C-C H -1"} = get_dist_neg_pos(7.87, 0.44, $lambda);
	$restraints_distance{"P C-C H 0"} = get_dist_neg_pos(5.22, 0.22, $lambda);
	$restraints_distance{"P N-N O +1"} = get_dist_neg_pos(8.14, 0.35, $lambda);
	$restraints_distance{"P N-N O -1"} = get_dist_neg_pos(4.86, 0.40, $lambda);
	$restraints_distance{"P N-N O 0"} = get_dist_neg_pos(5.13, 0.32, $lambda);
	$restraints_distance{"P N-N H +1"} = get_dist_neg_pos(4.80, 0.18, $lambda);
	$restraints_distance{"P N-N H -1"} = get_dist_neg_pos(7.54, 0.69, $lambda);
	$restraints_distance{"P N-N H 0"} = get_dist_neg_pos(5.10, 0.28, $lambda);
	$restraints_distance{"P CA-CA O +1"} = get_dist_neg_pos(8.55, 0.37, $lambda);
	$restraints_distance{"P CA-CA O -1"} = get_dist_neg_pos(4.90, 0.29, $lambda);
	$restraints_distance{"P CA-CA O 0"} = get_dist_neg_pos(6.21, 0.26, $lambda);
	$restraints_distance{"P CA-CA H +1"} = get_dist_neg_pos(4.90, 0.28, $lambda);
	$restraints_distance{"P CA-CA H -1"} = get_dist_neg_pos(7.49, 0.60, $lambda);
	$restraints_distance{"P CA-CA H 0"} = get_dist_neg_pos(6.24, 0.24, $lambda);
	$restraints_distance{"H O-O O +1"} = get_dist_neg_pos(8.40, 0.27, $lambda);
	$restraints_distance{"H O-O O -1"} = get_dist_neg_pos(4.99, 0.16, $lambda);
	$restraints_distance{"H O-O O 0"} = get_dist_neg_pos(6.12, 0.26, $lambda);
	$restraints_distance{"H O-O H +1"} = get_dist_neg_pos(5.03, 0.31, $lambda);
	$restraints_distance{"H O-O H -1"} = get_dist_neg_pos(8.43, 0.32, $lambda);
	$restraints_distance{"H O-O H 0"} = get_dist_neg_pos(6.12, 0.26, $lambda);
	$restraints_distance{"H C-C O +1"} = get_dist_neg_pos(8.16, 0.24, $lambda);
	$restraints_distance{"H C-C O -1"} = get_dist_neg_pos(4.87, 0.13, $lambda);
	$restraints_distance{"H C-C O 0"} = get_dist_neg_pos(6.09, 0.23, $lambda);
	$restraints_distance{"H C-C H +1"} = get_dist_neg_pos(4.89, 0.23, $lambda);
	$restraints_distance{"H C-C H -1"} = get_dist_neg_pos(8.17, 0.25, $lambda);
	$restraints_distance{"H C-C H 0"} = get_dist_neg_pos(6.09, 0.23, $lambda);
	$restraints_distance{"H N-N O +1"} = get_dist_neg_pos(8.07, 0.23, $lambda);
	$restraints_distance{"H N-N O -1"} = get_dist_neg_pos(4.84, 0.19, $lambda);
	$restraints_distance{"H N-N O 0"} = get_dist_neg_pos(6.10, 0.20, $lambda);
	$restraints_distance{"H N-N H +1"} = get_dist_neg_pos(4.81, 0.13, $lambda);
	$restraints_distance{"H N-N H -1"} = get_dist_neg_pos(8.08, 0.21, $lambda);
	$restraints_distance{"H N-N H 0"} = get_dist_neg_pos(6.10, 0.20, $lambda);
	$restraints_distance{"H CA-CA O +1"} = get_dist_neg_pos(8.63, 0.28, $lambda);
	$restraints_distance{"H CA-CA O -1"} = get_dist_neg_pos(5.13, 0.20, $lambda);
	$restraints_distance{"H CA-CA O 0"} = get_dist_neg_pos(6.16, 0.26, $lambda);
	$restraints_distance{"H CA-CA H +1"} = get_dist_neg_pos(5.14, 0.21, $lambda);
	$restraints_distance{"H CA-CA H -1"} = get_dist_neg_pos(8.64, 0.26, $lambda);
	$restraints_distance{"H CA-CA H 0"} = get_dist_neg_pos(6.16, 0.26, $lambda);

	# T Mean Standard_deviation
	$restraints_strandOO{"A"} = get_dist_neg_pos(4.57, 0.30, $lambda);
	$restraints_strandOO{"P"} = get_dist_neg_pos(4.57, 0.29, $lambda);
	$restraints_strandOO{"U"} = get_dist_neg_pos(4.57, 0.30, $lambda);

	# T Mean Standard_deviation
	$restraints_hbond{"A"} = get_dist_neg_pos(2.92, 0.16, $lambda);
	$restraints_hbond{"P"} = get_dist_neg_pos(2.93, 0.16, $lambda);
	$restraints_hbond{"H"} = get_dist_neg_pos(2.99, 0.17, $lambda);

	confess ":( dihe restraints could not be loaded" if not scalar keys %restraints_dihedral;
	confess ":( dstr restraints could not be loaded" if not scalar keys %restraints_strandOO;
	confess ":( dist restraints could not be loaded" if not scalar keys %restraints_distance;
	confess ":( hbnd restraints could not be loaded" if not scalar keys %restraints_hbond;

	system_cmd("rm -f $log_reference");
	foreach (keys %restraints_dihedral){ print2file($log_reference, $_."\t".$restraints_dihedral{$_});}
	foreach (keys %restraints_strandOO){ print2file($log_reference, $_."\t".$restraints_strandOO{$_});}
	foreach (keys %restraints_distance){ print2file($log_reference, $_."\t".$restraints_distance{$_});}
	foreach (keys %restraints_hbond){ print2file($log_reference, $_."\t".$restraints_hbond{$_});}
}

sub get_dist_neg_pos{
	my $mean = shift;
	my $devi = shift;
	my $lambda = shift;
	my $lambda_devi = sprintf "%.1f", $lambda * $devi;
	return $mean." ".$lambda_devi." ".$lambda_devi;
}

sub print_usage{
	my $error_message = shift;
	print "ERROR! $error_message!\n" if $error_message;
	print $param_info;
	print "\n";
	exit 1;
}

sub detect_hbonds_from_model{
	my $file_model = shift;
	my $file_hbond = shift;
	my %strands = strand_list($file_sec);
	my %final_pairs = ();
	my %final_aligned_pairs = ();
	confess ":( model does not exist $file_model " if not -f $file_model;
	my %pair_dist = ();
	my %all_pdb_all_xyz = ();
	my %pdb_xyz = all_xyz_pdb($file_model);
	# list all the possible pairings
	my %pairs_subsets_info = ();
	my %pairs_subsets_dist = ();
	foreach my $s1 (keys %strands){
		foreach my $s2 (keys %strands){
			next if $s2 eq $s1;
			my @A = split /\s+/, $s1;
			my @B = split /\s+/, $s2;
			next if $A[0] > $B[0];
			next if abs($A[0]-$A[1]) < 2; # implies at least 3 residues
			next if abs($B[0]-$B[1]) < 2; # implies at least 3 residues
			my %strand1_trims = ();
			my %strand2_trims = ();
			# (a) they are already of equal length
			if ($strands{$s1} == $strands{$s2}){  
				$strand1_trims{$s1} = 1;
				$strand2_trims{$s2} = 1;
			}
			# (b) find subsets of s1 to match the length as s2
			elsif ($strands{$s1} >= $strands{$s2}){
				$strand2_trims{$s2} = 1;
				for (my $i = 0; $i <= ($strands{$s1} - $strands{$s2}); $i++){
					confess ":(" if ($A[0]+$i+$strands{$s2}-1) > $A[1];
					$strand1_trims{($A[0]+$i)." ".($A[0]+$i+$strands{$s2}-1)} = 1;
				}
			}
			# (c) find subsets of s2 of same length as s1
			else{
				$strand1_trims{$s1} = 1;
				for (my $i = 0; $i <= ($strands{$s2} - $strands{$s1}); $i++){
					confess ":(" if ($B[0]+$i+$strands{$s1}-1) > $B[1];
					$strand2_trims{($B[0]+$i)." ".($B[0]+$i+$strands{$s1}-1)} = 1;
				}
			}
			my %subsets = ();
			my %possible_slides = qw( 0 1 L1 2 R1 3 L2 4 R2 5 ); # total 10 possible configurations
			my %possible_configs = qw( a 1 b 2 c 3 d 4);
			my %possible_types = qw( P 1 A 2 );
			foreach my $subset1 (keys %strand1_trims){
				foreach my $subset2 (keys %strand2_trims){
					my @X = split /\s+/, $subset1;
					my @Y = split /\s+/, $subset2;
					confess ":( strand size must be equal" if abs($X[0]-$X[1]) != abs($Y[0]-$Y[1]);
					my $i = 1;
					foreach my $slide (sort {$possible_slides{$a} <=> $possible_slides{$b}} keys %possible_slides){
						foreach my $config (sort {$possible_configs{$a} <=> $possible_configs{$b}} keys %possible_configs){
							foreach my $type (keys %possible_types){
								my @A = ($X[0], $X[1], $Y[0], $Y[1]);
								next if ($type eq "A" and not ($config eq "a" or $config eq "b")); # anti-parallel do not have hbond-patterns c and d
								my $no_of_residues = $slide;
								$no_of_residues =~ s/[A-Z]//g;
								# first strand is sliding forward
								if ($slide =~ /^R/){
									$A[1] = $A[1] - $no_of_residues;
									$A[3] = $A[3] - $no_of_residues if $type eq "A";
									$A[2] = $A[2] + $no_of_residues if $type eq "P";
								}
								# first strand is sliding backward
								elsif ($slide =~ /^L/){
									$A[0] = $A[0] + $no_of_residues;
									$A[2] = $A[2] + $no_of_residues if $type eq "A";
									$A[3] = $A[3] - $no_of_residues if $type eq "P";
								}
								confess "ERROR! strand size must be equal!" if abs($A[0]-$A[1]) != abs($A[2]-$A[3]);
								next if abs($A[0]-$A[1]) < 2; # must have at least 3 residues
								next if abs($A[2]-$A[3]) < 2; # must have at least 3 residues
								$subsets{$A[0]." ".$A[1]." ".$A[2]." ".$A[3]." $type $config"} = "[pattern $i] slide $slide type $type config $config " if $type eq "P";
								$subsets{$A[0]." ".$A[1]." ".$A[3]." ".$A[2]." $type $config"} = "[pattern $i] slide $slide type $type config $config " if $type eq "A";
								$i++;
							}
						}
					}
				}
			}
			$pairs_subsets_info{$s1." ".$s2} = \%subsets;
		}
	}
	system_cmd("rm -f pairs_and_distances.log");
	print2file("pairs_and_distances.log", "original_pair subset avg_ca_ca_dist avg_hbond_dist overall_farness pair_information");
	# for each strand pair, for each subset pair, compute average of distances in pdb files, for each A and P
	foreach my $pairing (keys %pairs_subsets_info){
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = ();
		foreach my $subset (keys %subsets_info){
			my @A = split /\s+/, $subset;
			# Step (1) measure the average backbone Ca-Ca distance
			my $avg_ca_ca = avg_ca_ca_dist(\%pdb_xyz, $A[0], $A[1], $A[2], $A[3], $A[4]);
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f [%s]", $pairing, $subset, $avg_ca_ca, $subsets_info{$subset});
			next if $avg_ca_ca > $pair_caca_thres;
			# Step (2) measure the average hbond distance for this configuration
			my $avg_hb_dist = avg_hbond_dist(\%pdb_xyz, $A[0], $A[1], $A[2], $A[3], $A[4], $A[5]);
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f : %.3f [%s]", $pairing, $subset, $avg_ca_ca, $avg_hb_dist, $subsets_info{$subset});
			next if $avg_hb_dist > $pair_hbond_thres;
			# Step (3) average them.
			my $farness = ($avg_ca_ca + $avg_hb_dist)/2;
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f : %.3f : %.3f : %s (selected)", $pairing, $subset, $avg_ca_ca, $avg_hb_dist, $farness, $subsets_info{$subset});
			$subsets_dist{$subset} = $farness if not defined $subsets_dist{$subset};
			$subsets_dist{$subset} = $farness if $subsets_dist{$subset} > $farness;
		}
		$pairs_subsets_dist{$pairing} = \%subsets_dist if scalar keys %subsets_dist;
	}
	if (not scalar keys %pairs_subsets_dist){
		print2file("pairs_and_distances.log", ":( not any eligible pairs!");
		return;
	}
	# best looking pairs
	print2file("pairs_and_distances.log", "\nBest looking pairs:");
	foreach my $pairing (sort keys %pairs_subsets_dist){
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		foreach (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist){
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %s : %s", $pairing, $_, $subsets_dist{$_}, $subsets_info{$_});
			last;
		}
	}
	# rank pairs
	my %pairs_rank = ();
	foreach my $pairing (keys %pairs_subsets_dist){
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		my $best_subset = (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist)[0];
		$pairs_rank{$pairing} = $subsets_dist{$best_subset};
	}
	# keep the best subset having minimum distance AND not violating the hbond pattern
	my %hbonds = ();
	my %final_pairing_info = ();
	my %final_pairing_dist = ();
	print2file("pairs_and_distances.log", "\nSelecting pairs:");
	foreach my $pairing (sort {$pairs_rank{$a} <=> $pairs_rank{$b}} keys %pairs_rank){
		print2file("pairs_and_distances.log", "\t$pairing");
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		# find the direction of top 1
		my $top1 = (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist)[0];
		my @T = split /\s+/, $top1;
		my $direction_of_top1 = $T[4];
		foreach my $best_subset (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist){
			print2line("pairs_and_distances.log", "\t\t$best_subset");
			my @T = split /\s+/, $best_subset;
			if($T[4] ne $direction_of_top1){
				print2line("pairs_and_distances.log", " -> direction is different.. ignoring rest!");
				last;
			}
			my $addable = 0;
			$addable = is_addable_to_hbond_list($best_subset, %hbonds);
			if ($addable){
				my @S = split /\s+/, $best_subset;
				my %this_hbonds = make_ideal_hbonds($S[0], $S[1], $S[2], $S[3], $S[4], $S[5]);
				foreach (keys %this_hbonds){
					$hbonds{$_} = $S[4];
				}
				$final_pairing_info{$best_subset} = $subsets_info{$best_subset};
				$final_pairing_dist{$best_subset} = $subsets_dist{$best_subset};
				last;
			}
			else{
				print2line("pairs_and_distances.log", " -> could not add because of hbond issues!");
			}
			print2file("pairs_and_distances.log", "");
		}
		print2file("pairs_and_distances.log", "");
	}
	if (not scalar keys %hbonds){
		print2file("pairs_and_distances.log", "Quitting because no strands could be paired!");
		return;
	}
	print2file("pairs_and_distances.log", "\nHbonds:");
	foreach (keys %hbonds){
		print2file("pairs_and_distances.log", $_." ".$hbonds{$_});
	}
	# log final pairs
	print2file("pairs_and_distances.log", "\nFinal pairs:");
	foreach (sort {$final_pairing_dist{$a} <=> $final_pairing_dist{$b}} keys %final_pairing_dist){
		print2file("pairs_and_distances.log", sprintf "%15s : %.3f : %s", $_, $final_pairing_dist{$_}, $final_pairing_info{$_});
	}
	# print to the output files, finally
	system_cmd("rm -f pairing.txt");
	foreach (sort {$final_pairing_dist{$a} <=> $final_pairing_dist{$b}} keys %final_pairing_dist){
		my @A = split /\s+/, $_;
		# converting distance to confidence
		print2file("pairing.txt", sprintf $A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." %.2f ".$A[5]."  # average distance ".$final_pairing_dist{$_}." [".$final_pairing_info{$_}."]", (1/$final_pairing_dist{$_}));
	}
	system_cmd("rm -f $file_hbond");
	foreach (sort keys %hbonds){
		print2file($file_hbond, $_." ".$hbonds{$_});
	}
}

sub avg_ca_ca_dist{
	my %pdb_xyz = %{shift @_};
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	my $sum_dist = 0.0;
	my $count_dist = 0;
	confess "ERROR! type must be A or P in $a $b $c $d $T!" if not ($T eq "A" or $T eq "P");
	if ($T eq "P"){
		foreach (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i++, $j++){
			my @row1 = split /\s+/, $pdb_xyz{$i." CA"};
			my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
			my @row2 = split /\s+/, $pdb_xyz{$j." CA"};
			my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
			$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			$count_dist++;
		}
	}
	#foreach (my $i = $A[0], my $j = $A[3]; $i <= $A[1] and $j >= $A[2]; $i++, $j--){
	else{
		confess "ERROR! anti-parallel config must have c > d in $a $b $c $d $T!" if $c < $d;
		foreach (my $i = $a, my $j = $c; $i <= $b and $j >= $d; $i++, $j--){
			my @row1 = split /\s+/, $pdb_xyz{$i." CA"};
			my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
			my @row2 = split /\s+/, $pdb_xyz{$j." CA"};
			my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
			$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			$count_dist++;
		}
	}
	return sprintf "%.3f", ($sum_dist/$count_dist);
}

sub avg_hbond_dist{
	my %pdb_xyz = %{shift @_};
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	my $C = shift;
	my $sum_dist = 0.0;
	my $count_dist = 0;
	my %ideal_hbonds = make_ideal_hbonds($a, $b, $c, $d, $T, $C);
	foreach (keys %ideal_hbonds){
		# some complications here because Proline does not have H atom.
		my @H = split /\s+/, $_;
		my @row1;
		@row1 = split /\s+/, $pdb_xyz{$H[0]." N"} if ($residues{$H[0]} eq "P" and $H[1] eq "H");
		@row1 = split /\s+/, $pdb_xyz{$H[0]." ".$H[1]} if not ($residues{$H[0]} eq "P" and $H[1] eq "H");
		my @row2;
		@row2 = split /\s+/, $pdb_xyz{$H[2]." N"} if ($residues{$H[2]} eq "P" and $H[3] eq "H");
		@row2 = split /\s+/, $pdb_xyz{$H[2]." ".$H[3]} if not ($residues{$H[2]} eq "P" and $H[3] eq "H");
		my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
		my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
		$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
		$count_dist++;
	}
	return sprintf "%.3f", ($sum_dist/$count_dist);
}

sub validate_pairing_file{
	my $file_pairing_info = shift;
	open PAIR, $file_pairing_info or confess $!;
	while (<PAIR>){
		chomp $_;
		$_ =~ s/^\s+//;
		next unless $_;
		confess "Error! Numbers expected at beginning in $_" if $_ !~ /^[0-9]/; 
		my @C = split /\s+/, $_;
		for( my $i=0; $i <= 5; $i++){
			confess "Error! Column $i is missing in $_" if not defined $C[$i];
		}
		#confess "Unequal length pair in $_" if (abs($C[0]-$C[1]) != abs($C[2]-$C[3]));
		confess "Error! a must be less than b in $_" if ($C[0] >= $C[1]);
		confess "Error! Pairing must be A or P in $_" if ($C[4] ne "A" and $C[4] ne "P");
		confess "Error! c must be less than d for parallel pairing in $_" if ($C[2] >= $C[3] and $C[4] eq "P");
		confess "Error! c must be greater than d for anti-parallel pairing in $_" if ($C[2] <= $C[3] and $C[4] eq "A");
	}
}

sub str2pairing_info{
	my $file_bpro_str = shift;
	my $file_pairing_info = shift;
	my %pairs = bpro_pairs_filtered($file_bpro_str);
	open PAIR, ">$file_pairing_info" or confess $!;
	foreach (sort {$pairs{$a} <=> $pairs{$b}} keys %pairs){
		my @I = split(/\s+/, $_);
		printf PAIR "%d %d %d %d %s %.2f\n", $I[3], $I[4], $I[5], $I[6], $I[2], $I[$#I];
		confess ":( expecting A or P!" if (not ($I[2] eq "P" or $I[2] eq "A"));
	}
	close PAIR;
}

# If the pairing information is out of the strand residues, extend the secondary structure information
sub update_sec_using_pairing_info{
	my $file_sec = shift;
	my $file_pairing_info = shift;
	my $file_updated_ss = shift;
	my %pairing_info = pairing_info2hash($file_pairing_info);
	my %strands = strand_list($file_sec);
	my %paired_residues = ();
	foreach my $pair (keys %pairing_info){
		my @A = split /\s+/, $pair;
		confess "ERROR! $pair pair not defined!" if not (defined $A[0] and defined $A[1]);
		for(my $i = $A[0]; $i <= $A[1]; $i++){
			$paired_residues{$i} = 1;
		}
		for(my $i = $A[2]; $i <= $A[3]; $i++){
			$paired_residues{$i} = 1;
		}
		for(my $i = $A[3]; $i <= $A[2]; $i++){
			$paired_residues{$i} = 1;
		}
	}
	my %res_ss = each_residue_sec_ss($file_sec);
	foreach (sort {$a <=> $b} keys %res_ss){
		$res_ss{$_} = "E" if defined $paired_residues{$_};
	}
	my $line_sec = "";
	foreach (sort {$a <=> $b} keys %res_ss){
		$line_sec .= $res_ss{$_};
	}
	my $line_updated_E = "";
	foreach (sort {$a <=> $b} keys %res_ss){
		if (defined $paired_residues{$_}){
			$line_updated_E .= "E";
			next;
		}
		$line_updated_E .= "-";
	}
	confess " $line_sec :( " if length(ss_line_ss($file_sec)) ne length($line_sec);
	print "\n".ss_line_ss($file_sec)." [".basename($file_sec)."]";
	print "\n".$line_sec." [updated]";
	print "\n".$line_updated_E." [paired strands]";
	open SSNEW, ">$file_updated_ss" or confess $!;
	print SSNEW ">$seq_id\n";
	print SSNEW "".seq_ss($file_sec)."\n";
	print SSNEW "$line_sec\n";
	close SSNEW;
}

sub pairing_info2hash{
	my $file_pairing_info = shift;
	my %pairing_rows = ();
	open PAIR, $file_pairing_info or confess $!;
	while (<PAIR>){
		chomp $_;
		next if not defined $_;
		$_ =~ s/^\s+//;
		my @A = split /\s+/, $_;
		confess "ERROR! all columns not defined in $_" if not (defined $A[0] and defined $A[1] and defined $A[2] and defined $A[3] and defined $A[4]); 
		$pairing_rows{$_} = $A[0];
	}
	close PAIR;
	return %pairing_rows;
}

sub paired_residues{
	my $file_pairing_info = shift;
	my %pairing_info = pairing_info2hash($file_pairing_info);
	my %paired_residues = ();
	foreach my $pair (keys %pairing_info){
		my @A = split /\s+/, $pair;
		for(my $i = $A[0]; $i <= $A[1]; $i++){
			$paired_residues{$i} = $A[4];
		}
		for(my $i = $A[2]; $i <= $A[3]; $i++){
			$paired_residues{$i} = $A[4];
		}
		for(my $i = $A[3]; $i <= $A[2]; $i++){
			$paired_residues{$i} = $A[4];
		}
	}
	return %paired_residues;
}

sub load_hbonds{
	my $file_hbonds = shift;
	my %hbond_rows = ();
	open HBOND, $file_hbonds or confess $!;
	while (<HBOND>){
		my @A = split /\s+/, $_;
		$hbond_rows{$A[0]." ".$A[1]." ".$A[2]." ".$A[3]} = $A[4];
	}
	close HBOND;
	return %hbond_rows;
}

sub hbonds2connectors{
	my $ref_hbonds = shift;
	my %hbonds = %{$ref_hbonds};
	my %connectors = ();
	foreach (sort keys %hbonds){
		my @H = split /\s+/, $_;
		$connectors{$H[0]." ".$H[1]} = 1;
		$connectors{$H[2]." ".$H[3]} = 1;
	}
	return %connectors;
}

sub strand_list{
	my $file_sec = shift;
	my %each_residue_ss = each_residue_sec_ss($file_sec);
	my %strand_list = ();
	my $start = 1;
	foreach my $r (sort {$a <=> $b} keys %each_residue_ss){
		next if $each_residue_ss{$r} ne "E";
		$start = $r if (defined $each_residue_ss{$r-1} && $each_residue_ss{$r-1} ne "E");
		$strand_list{$start." ".$r} = $r - $start + 1 if (defined $each_residue_ss{$r+1} && $each_residue_ss{$r+1} ne "E");
	}
	return %strand_list;
}

sub make_ideal_hbonds{
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	my $C = shift; # configuration a or b (look at the pictures folder) 
	confess "ERROR! unknown type! type must be P or A" if not ($T eq "P" or $T eq "A"); 
	confess "ERROR! unknown configuration! configuration must be a or b" if not ($C eq "a" or $C eq "b" or $C eq "c" or $C eq "d"); 
	my %hbond_list = ();
	if ($T eq "A"){
		confess "ERROR! in anti-parallel strand pairing in [$T $a-$b:$c-$d]. In [a-b:c-d] c must be greater than d" if $c < $d;
		if ($C eq "a"){
			for (my $i = $a, my $j = $c; $i <= $b and $j >= $d; $i = $i+2,  $j = $j-2){
				$hbond_list{"$i H $j O"} = $T;
				$hbond_list{"$i O $j H"} = $T;
			}
		}
		else{
			for (my $i = $a+1, my $j = $c-1; $i <= $b and $j >= $d; $i = $i+2,  $j = $j-2){
				$hbond_list{"$i H $j O"} = $T;
				$hbond_list{"$i O $j H"} = $T;
			}
		}
	}
	else{
		confess "ERROR! in parallel strand pairing in [$T $a-$b:$c-$d]. In [a-b:c-d] c must be less than d" if $d < $c;
		# 9/2/2014:
		# I noticed that my output number of EEE in parallel pairs is less than input, in my "test1" results.
		# Looking at how Chimera shows beta sheets and how DSSP calculates E for a PDB, it seems this:
		# In a native PDB (I used 1a05a's first parallel), if there are two strands 3-6 and 36-39, forming a parallel beta sheet, then the actual residues involved in hydrogen bonding are 2-6 and 36-40.
		# But, somehow, DSSP and Chimera both ignore the first hydrogen-bonded residue and last hydrogen-bonded residue.
		# This makes me start hydrogen bonding from $a-1 and go up to $d+1
		# 12/24/2014:
		# Seeing 1yipA's hydrogen bonding 122N-158O, in reindexed pdb, and looking at DSSP's interpretations, it 
		# is clear that pairs starting with N-O hydrogen bonds are ignored, and so my assumption of ideal h-bond starting with O-H bonding may
		# be all right.
		# 12/26/2014:
		# I know that this is not appropriate in terms of reconstruction.
		# I must go up to 1 more and so should j. But, that way, it will work great for reconstruction but not for prediction (my intuition)
		if ($C eq "a"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{"$i O ".($j+1)." H"} = $T if $j+1 <= $d;
				$hbond_list{"".($i+2)." H ".($j+1)." O"} = $T if ($i+2 <= $b and $j+1 <= $d);
			}
		}
		elsif ($C eq "b"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{($i+1)." H $j O"} = $T if $i+1 <= $b;
				$hbond_list{($i+1)." O ".($j+2)." H"} = $T if ($i+1 <= $b and $j+2 <= $d);
			}
		}
		elsif ($C eq "c"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{"$i H ".($j+1)." O"} = $T if $j+1 <= $d;
				$hbond_list{"".($i+2)." O ".($j+1)." H"} = $T if ($i+2 <= $b and $j+1 <= $d);
			}
		}
		else{
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{($i+1)." O $j H"} = $T if $i+1 <= $b;
				$hbond_list{($i+1)." H ".($j+2)." O"} = $T if ($i+1 <= $b and $j+2 <= $d);
			}
		}
	}
	return %hbond_list;
}

sub bpro_pairs_filtered{
	my $file_bpro = shift;
	confess "Betacon file $file_bpro?" if not -f $file_bpro;
	my %pair_confid = ();
	open BETACON, $file_bpro or confess $!;
	while (<BETACON>){
		my $line = $_;
		next if $line !~ m/--/;
		next if $line =~ m/---/;
		next if $line =~ m/x--y/;
		next if $line =~ m/strand/;
		chomp $line;
		$line =~ s/\s+//g; # 14--15:A:[234-240:252-246]:4.07
		$line =~ s/\[//g;
		$line =~ s/\]//g;
		$line =~ s/--/ /g;
		$line =~ s/:/ /g;
		$line =~ s/-/ /g;  # 14 15 A 234 240 252 246 4.07
		my @C = split(/\s+/, $line);
		for(my $i = 0; $i <= $#C; $i++){
			confess "column $i+1 not defined in $_ in $file_bpro" if not defined $C[$i];
		}
		# Don't keep low confidence pairs
		next if $C[$#C] < $bpro_thres;
		$pair_confid{$line} = $C[$#C];
	}
	close BETACON;
	# Get rid of very short strands
	foreach (sort {$pair_confid{$b} <=> $pair_confid{$a}} (keys %pair_confid)){
		my @I = split(/\s+/, $_);
		if (abs($I[3]-$I[4]) < 2){
			delete $pair_confid{$_};
			next;
		}
		if (abs($I[6]-$I[5]) < 2){
			delete $pair_confid{$_};
			next;
		}
	}
	# Get rid of cases when one strand is involved in three pairings. 
	# Based on confidence, remove the third pair, in such a case.
	my %strandPresence = ();
	foreach my $pair (sort {$pair_confid{$b} <=> $pair_confid{$a}} (keys %pair_confid)){
		my @I = split(/\s+/, $pair);
		my @strands = ( $I[0], $I[1]); 
		foreach my $st (@strands){
			if(not defined $strandPresence{$st}){
				$strandPresence{$st} = 1;
			}
			else{
				$strandPresence{$st}++;
				if ($strandPresence{$st} eq 3){
					delete $pair_confid{$pair};
					$strandPresence{$st}--;
				}
			}
		}
	}
	foreach (keys %pair_confid){
		my @I = split(/\s+/, $_);
		$pair_confid{$_} = $I[3];
	}
	print " loaded beta pro predictions with confidence >= $bpro_thres\n";
	return %pair_confid;
}

sub is_addable_to_hbond_list{
	my $pairing = shift;
	my %hbonds = @_;
	my %connectors = ();
	my @P = split /\s+/, $pairing;
	my %ideal_hbonds = make_ideal_hbonds($P[0], $P[1], $P[2], $P[3], $P[4], $P[5]);
	foreach (keys %hbonds){
		my @H = split /\s+/, $_;
		$connectors{$H[0]." ".$H[1]} = 1;
		$connectors{$H[2]." ".$H[3]} = 1;
	}
	my $flag_problem = 0;
	foreach (keys %ideal_hbonds){
		my @H = split /\s+/, $_;
		$flag_problem = 1 if defined $connectors{$H[0]." ".$H[1]};
		$flag_problem = 1 if defined $connectors{$H[2]." ".$H[3]};
	}
	return 1 if not $flag_problem;
	return 0;
}
