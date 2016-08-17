#!/usr/bin/perl
use strict;
use warnings;
use CGI qw(param);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser); 
use File::Basename;
use MIME::Lite;
use Cwd 'abs_path';
use lib dirname(abs_path($0));
use confold;

$CGI::POST_MAX = 10000;
$CGI::DISABLE_UPLOADS = 1;

print "Content-type: text/html\n\n";

# input parameters
my $param_email       = param("email");
my $param_id          = param("id");
my $param_seq         = param("protein_sequence");
my $param_sec         = param("protein_sec");
my $param_rr          = param("rr");
my $param_rr_subset   = param("rr_subset");
my $param_rr_type     = param("rr_type");
my $param_lambda      = param("lambda");
my $param_stage2      = param("stage2");
my $param_sheet_thr   = param("sheet_threshold");
my $param_con_wt      = param("con_wt");
my $param_sec_wt      = param("sec_wt");
my $param_atom_scheme = param("atom_scheme");
my $param_rep2        = param("rep2");
my $param_pairing     = param("pairing");
my $remote_addr       = $ENV{'REMOTE_ADDR'};
my $now = (localtime);
$now    =~ s/\s+/_/g;
$now    =~ s/:/_/g;
my $job_id = $param_id."_".$now;

# global parameters
my $dir_all_jobs = "/var/www/html/confold/jobs";
my $file_history = "/var/www/html/confold/logs/history.log";
my $file_error   = "/var/www/html/confold/logs/errors.log";
my $file_log     = "$dir_all_jobs/confold_$job_id.log";
my $web_addr     = "http://protein.rnet.missouri.edu/confold";
my $stg2_folder;
if($param_stage2 ne "1"){
	$stg2_folder = "sheet_detect" if $param_stage2 eq "2";
	$stg2_folder = "con_filter" if $param_stage2 eq "3";
	$stg2_folder = "stage2" if $param_stage2 eq "4";
}

validate_parameters();
sleep(1);

my $pid=fork;
if (!$pid){
	close(STDIN); close(STDOUT);close(STDERR);
	# make platform
	my $dir_job = "$dir_all_jobs/$job_id";
	mkdir $dir_job or log_n_exit($file_error, "could not make $dir_job"); ;
	chdir $dir_job or log_n_exit($file_error, "could not change to $dir_job");
	# record and start job
	log_parameters();
	create_files();
	my $options = "-dir_job $dir_job";
	$options .= " --allrr";
	$options .= " --no_ss_wt2";
	$options .= " -rr_type $param_rr_type";
	$options .= " -lambda $param_lambda";
	$options .= " -cpu 5";
	$options .= " -atom_scheme $param_atom_scheme";
	$options .= " --sheet_detect" if $param_stage2 eq "2";
	$options .= " --con_filter"   if $param_stage2 eq "3";
	$options .= " --stage2"       if $param_stage2 eq "4";
	$options .= " -rep2 $param_rep2";
	$options .= " -pair_thres $param_sheet_thr";
	$options .= " -pairing $dir_job/pairing.txt" if defined $param_pairing;
	$options .= " -cont_wt $param_con_wt";
	$options .= " -ss_wt1 $param_sec_wt";
	exec_system(dirname(abs_path($0))."/confold.pl $options &> $file_log", $file_error); 
	# zip stage 1 and stage 2 top 5 models
	my $results_location = "$dir_job/to_email";
	mkdir $results_location or log_n_exit($file_error, "could not make $results_location");
	chdir $results_location or log_n_exit($file_error, "could not change to $results_location");
	for (my $i = 1; $i <= 5; $i++){
		exec_system("cp $dir_job/stage1/sec_${param_sec_wt}_allL/${param_id}_model$i.pdb ./${param_id}_stage1_$i.pdb", $file_error);
		exec_system("cp $dir_job/$stg2_folder/sec_${param_sec_wt}_allL/${param_id}_model$i.pdb ./${param_id}_${stg2_folder}_$i.pdb", $file_error) if $param_stage2 ne "1";
	}
	exec_system("zip ${param_id}_models.zip *.pdb", $file_error);
	send_email();
}
else{
	open HEADER, "/var/www/html/confold/header.html"; 
	while(<HEADER>){
		chomp $_;
		print $_;
	}
	close HEADER;
	print "<tr bgcolor=\"#EEEEEE\">";
	print "<td align=\"left\">";
	print "<p>";
	print "<b>Thank you! Your job is submitted!</b></br></br>";
	print "All files related to this job (id = $job_id) can be browsed <a href=\"$web_addr/jobs/$job_id/\">here</a>.</br></br>";
	print "Results will be available for download at <a href=\"$web_addr/jobs/$job_id/to_email/${param_id}_models.zip\">${param_id}_models.zip</a> when the job is finished.</br>";
	print "This download URL will also be emailed to $param_email.</br>";
	print "Usually it takes less than an hour (depends on the sequence length and our server load).</br></br>";
	print "Log files:</br>";
	print " (a) <a href=\"$web_addr/jobs/confold_$job_id.log\">Main Log</a><br>";
	print " (b) <a href=\"$web_addr/jobs/$job_id/stage1/sec_${param_sec_wt}_allL/assess.log\">Stage 1 model assessment log</a> (available after stage 1 is finished)<br>";
	print " (c) <a href=\"$web_addr/jobs/$job_id/$stg2_folder/sec_${param_sec_wt}_allL/assess.log\">Stage 2 model assessment log</a> (available after stage 2 is finished)<br>" if defined $stg2_folder;
	print "</p>";
	print "</td></tr>";	
	open FOOTER, "/var/www/html/confold/footer.html"; 
	while(<FOOTER>){
		chomp $_;
		print $_;
	}
	close FOOTER;
	exit;
}

sub validate_parameters{
	# strict sequence check
	chomp $param_seq;
	if (length($param_seq) < 5){
		error_exit("sequence $param_seq is undefined or too short");
	}
	if (length($param_seq) > 500){
		error_exit("sequence $param_seq is too long. 500 is the limit.");
	}
	for(my $i = 1; $i <= length($param_seq); $i++){
		my $char = substr $param_seq, $i-1, 1;
		if (not defined $AA1TO3{$char}){
			error_exit("undefined amino acid $char in $param_seq");
		}
	}
	if (length($param_email) < 5){
		error_exit("check email $param_email");
	}
	# string secondary structure check
	$param_sec = undef if (defined $param_sec and length($param_sec) < 5);
	if (defined $param_sec){
		chomp $param_sec;
		for(my $i = 1; $i <= length($param_sec); $i++){
			my $char = substr $param_sec, $i-1, 1;
			if (not ($char eq "H" or $char eq "C" or $char eq "E")){
				error_exit("undefined secondary structure unit $char in $param_sec");
			}
		}
		if(length($param_sec) ne length($param_seq)){
			error_exit("length of protein sequence and secondary structure do not match! $param_seq $param_sec");
		}
	}
	# id cannot be very short
	if (length($param_id) < 3){
		error_exit("too short id $param_id");
	}
	# rr input must be in CASP format
	$param_rr = undef if (defined $param_rr and length($param_rr) < 5);
	if(defined $param_rr){
		$param_rr =~ s/\r//g;
		my @RR = split /\n/, $param_rr;
		foreach my $line (@RR){
			next if length($line) < 2;
			my @C = split /\s+/, $line;
			error_exit("column 1 not defined in $line in rr input") if not defined $C[0];
			error_exit("column 2 not defined in $line in rr input") if not defined $C[1];
			error_exit("column 3 not defined in $line in rr input") if not defined $C[2];
			error_exit("column 4 not defined in $line in rr input") if not defined $C[3];
			error_exit("column 5 not defined in $line in rr input") if not defined $C[4];
		}
	}
	# check pairing input format
	$param_pairing = undef if (defined $param_pairing and length($param_pairing) < 5);
	if(defined $param_pairing){
		$param_pairing =~ s/\r//g;
		my @RR = split /\n/, $param_pairing;
		foreach my $line (@RR){
			next if length($line) < 2;
			my @C = split /\s+/, $line;
			error_exit("column 1 not defined in $line in rr input") if not defined $C[0];
			error_exit("column 2 not defined in $line in rr input") if not defined $C[1];
			error_exit("column 3 not defined in $line in rr input") if not defined $C[2];
			error_exit("column 4 not defined in $line in rr input") if not defined $C[3];
			error_exit("column 5 not defined in $line in rr input") if not defined $C[4];
			error_exit("column 6 not defined in $line in rr input") if not defined $C[5];
		}
	}
}

sub log_parameters{
	print2file($file_history, "email: ".$param_email);
	print2file($file_history, "id   : ".$param_id);
	print2file($file_history, "seq  : ".$param_seq);
	print2file($file_history, "sec  : ".$param_sec) if defined $param_sec;
	print2file($file_history, "rr   : ".substr($param_rr, 1, 500)) if defined $param_rr;
	print2file($file_history, "sub  : ".$param_rr_subset);
	print2file($file_history, "type : ".$param_rr_type);
	print2file($file_history, "lam  : ".$param_lambda);
	print2file($file_history, "stg2 : ".$param_stage2);
	print2file($file_history, "bdetc: ".$param_sheet_thr);
	print2file($file_history, "conwt: ".$param_con_wt);
	print2file($file_history, "secwt: ".$param_sec_wt);
	print2file($file_history, "atoms: ".$param_atom_scheme);
	print2file($file_history, "pairs: ".$param_pairing) if defined $param_pairing;
	print2file($file_history, "time : ".$now);
	print2file($file_history, "jobid: ".$job_id);
	print2file($file_history, "ip   : ".$remote_addr."\n");
}

sub create_files{
	print2file("$param_id.fasta", ">$param_id");
	print2file("$param_id.fasta", $param_seq);
	print2file("$param_id.ss", ">$param_id");
	print2file("$param_id.ss", $param_seq);
	print2file("$param_id.ss", $param_sec);
	if ($param_rr_subset eq "all"){
		print2file("$param_id.rr", $param_seq);
		print2file("$param_id.rr", $param_rr);
	}
	else{
		print2file("$param_id.rr.original", $param_seq);
		print2file("$param_id.rr.original", $param_rr);
		my $seq_len = length($param_seq); 
		my $top_xl = int(($param_rr_subset * $seq_len) + 0.5) + 1; # +1 to account for header line
		system_cmd("head -n $top_xl $param_id.rr.original > contact.rr");
	}
	if (defined $param_pairing){
		print2file("pairing.txt", $param_pairing);
	}
}

sub log_n_exit{
	my $file_log = shift;
	my $message = shift;
	if (-f $file_log){
		open  LOG, ">>$file_log" or die $!;
		print LOG "[".(localtime)."] ".$message."\n";
		close LOG;
	}
	else{
		open  LOG, ">$file_log" or die $!;
		print LOG "[".(localtime)."] ".$message."\n";
		close LOG;
	}
	exit;
}

sub exec_system{
	my $command = shift;
	my $log = shift;
	my $code = system($command);
	my $error = $!;
	if($code != 0){
		print_log($log, "COMMAND : [$command]");
		print_log($log, "ERROR   : $error");
		print2file($log, "");
		exit;
	}
}

sub error_exit{
	my $message = shift;
	open HEADER, "/var/www/html/confold/header.html"; 
	while(<HEADER>){
		chomp $_;
		print $_;
	}
	close HEADER;
	print "<tr bgcolor=\"#EEEEEE\">";
	print "<td align=\"middle\">";
	print "<p>";
	print "ERROR! $message";
	print "</p>";
	print "</td></tr>";
	open FOOTER, "/var/www/html/confold/footer.html"; 
	while(<FOOTER>){
		chomp $_;
		print $_;
	}
	close FOOTER;
	exit;
}

sub send_email{
	my $body = "Your job $job_id has completed!";
	$body .= "\n\nDownload top models at $web_addr/jobs/$job_id/to_email/${param_id}_models.zip";
	$body .= "\n\nAll files related to this job can be browesed at $web_addr/jobs/$job_id/";
	$body .= "\n\nMain log file $web_addr/jobs/confold_$job_id.log";
	$body .= "\nStage 1 models' assessment at $web_addr/jobs/$job_id/stage1/sec_${param_sec_wt}_allL/assess.log";
	$body .= "\nStage 2 models' assessment at $web_addr/jobs/$job_id/$stg2_folder/sec_${param_sec_wt}_allL/assess.log" if defined $stg2_folder;
	print2file($file_history, "attempting to send email for job $job_id ..");
	my $msg = MIME::Lite->new(
			From => "CONFOLD SERVER<chengji\@missouri.edu>",
			To => "$param_email",
			Subject => "CONFOLD job finished for $param_id",
			Bcc => "bap54\@mail.missouri.edu",
			Data => $body);
	$msg->send;
	print2file($file_history, "finished sending email for job $job_id.\n");
}