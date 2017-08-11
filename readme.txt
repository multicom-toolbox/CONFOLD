-----------------------------------------------------------
Dependencies
-----------------------------------------------------------
1. CONFOLD is implemented to run in Linux environment.
   It was tested in "x86_64 GNU/Linux" OS.
2. Perl v5.10.1 was used during development and testing;
   But it should run with other versions of perl as well.
3. CNS suite
4. DSSP
5. TM-score (only for executing the comprehensive test suite)

-----------------------------------------------------------
Installation
-----------------------------------------------------------
1. Download DSSP
   1.1 Download DSSP
       $ wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64
   1.2 Make it executable
       $ chmod +x dssp-2.0.4-linux-amd64
   1.3 Test it
       $ ./dssp-2.0.4-linux-amd64

2. Install CNS suite
   2.1. To download CNS suite, provide your academic profile related 
        information at http://cns-online.org/cns_request/. An email
        with (a) link to download, (b) login, and (c) password
        will be sent to you. Follow the link, possibly
        http://cns-online.org/download/, and download 
        CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
   2.2. Unzip
        $ tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
   2.3. Change directory to cns_solve
        $ cd cns_solve_1.3
   2.4. Unhide the file '.cns_solve_env_sh'
        $ mv .cns_solve_env_sh cns_solve_env.sh
   2.5. Edit 'cns_solve_env.sh' and 'cns_solve_env' to replace
        '_CNSsolve_location_' with CNS installation directory.
        For instance, if your CNS installation path is
        '/home/user/programs/cns_solve_1.3' replace
        '_CNSsolve_location_' with this path
   2.6. Test CNS installation
        $ source cns_solve_env.sh
        $ cd test 
        $ ../bin/run_tests -tidy *.inp
 
3. Download confold_v1.0.tar.gz
   3.1 Download confold_v1.0.tar.gz if you don't have it.
   3.2 Untar
       $ tar zxvf confold_v1.0.tar.gz
   3.3 Change directory to confold_v1.0
       $ cd confold_v1.0

4. Change variable values in the confold.pl file
   4.1 Change the path of the variable $dssp to DSSP executable
   4.2 Change the path of the variable $cns_suite 
       to CNS installation directory
   4.3 Make it executable
       $chmod +x confold.pl
  
5. Test CONFOLD
   5.1 Execute "perl ./confold.pl" or "./confold.pl"
       It should print the usage information.
   5.2 Test using a short example
       $ ./confold.pl -rrtype ca -stage2 1 -mcount 5 -seq ./test/input/short.fasta -ss ./test/input/short.ss -rr ./test/input/short.rr -o ./test/output/short
   5.3 (Optional) Visualize the top model 'short_model1.pdb' 
       in ./test/output/stage2/ folder using a pdb visualization tool
       like USEF Chimera or PyMol or JMol.   
   5.4 (Optional) For a more comprehensive testing see the section below.

To learn about execution time of CONFOLD please visit
http://protein.rnet.missouri.edu/confold/tool.php

-----------------------------------------------------------
File Formats
-----------------------------------------------------------
Below is the description of the file formats of the input files.
The best way to learn about these files, however, is to see
the examples in ./test/input/ folder.

Fasta: This file format is described at 
https://en.wikipedia.org/wiki/FASTA_format. Although not recommended
for readability, the lines may be longer than 80 characters as well.

Contacts: CONFOLD accepts CASP's RR format files as input. See, 
http://predictioncenter.org/casprol/index.cgi?page=format#RR.
For simpliciy, CONFOLD also accepts files having sequence in the 
first line followed by contact rows.

Secondary Structure: The format is same as FASTA file format with 
the residue names replaced by their 3-stage secondary structure 
(i.e. H, E, or C).

Beta-sheet Pairing File:
- 5 columns a, b, c, d, and t in each row
- a-b and c-d are residue strands, for example 2-7 and 20-25
- t is the pairing type (A or P)
- a must always be less than b
- c must be less than d if parallel and greater than d if anti-parallel

-----------------------------------------------------------
Comprehensive Testing
-----------------------------------------------------------
We have collected 7 test cases to test the CONFOLD perl script, 
CNS suite installation, and other configurations needed to run CONFOLD. 
An image of all these proteins, input.png, is in 
the ./test/input/ directory. Each test case uses contacts and 
secondary structure as input, either true or predicted. 
However, some test cases use only secondary structure information as input. 
The test cases are listed below. 

1. Reconstruction of 1VJK, a helix and anti-parallel 
   mixed protein, using true contacts only
2. Reconstruct a 50 residue long helix
3. Fold 1EAZ, a helix and anti-parallel mixed protein, 
   using predicted contacts and secondary structures
4. Fold 1GUU, a helical protein, using predicted 
   contacts and secondary structures
5. Fold 1SMX, a parallel and anti-parallel beta-sheet 
   protein, using true contacts and secondary structures
6. Reconstruct 1QJP, an anti-parallel beta barrel protein, using 
   true secondary structure and pairing information
7. Reconstruct 1G7R, a helix and parallel beta-sheet 
   protein, using true contacts, secondary structures and pairing information

For a comprehensive test of CONFOLD, please run 
the 7 CONFOLD jobs below and evaluate the models against
the native models provided in the input using TM-score.

./confold.pl -seq test/input/helix.fasta -ss test/input/helix.ss -o test/output/helix -sswt 10 -mcount 20 &> test/output/helix.log &
./confold.pl -seq test/input/1vjk.fasta  -rr test/input/1vjk.rr  -o test/output/1vjk  -contwt 50 -pthres 6.5 -rep2 0.8 -lambda 1.0  -mcount 20 &> test/output/1vjk.log &
./confold.pl -seq test/input/1eaz.fasta  -rr test/input/1eaz.rr -ss test/input/1eaz.ss -o test/output/1eaz -selectrr 1.0L -stage2 1 -mcount 20 &> test/output/1eaz.log &
./confold.pl -seq test/input/1guu.fasta  -rr test/input/1guu.rr -ss test/input/1guu.ss -o test/output/1guu -selectrr 0.8L -stage2 1 -mcount 20 &> test/output/1guu.log &
./confold.pl -seq test/input/1smx.fasta  -rr test/input/1smx.rr -ss test/input/1smx.ss -o test/output/1smx -stage2 3 -rrtype ca -mcount 20 &> test/output/1smx.log &
./confold.pl -seq test/input/1qjp.fasta -pair test/input/1qjp.pair -ss test/input/1qjp.ss -o test/output/1qjp -mcount 20 &> test/output/1qjp.log &
./confold.pl -seq test/input/1g7r.fasta -pair test/input/1g7r.pair -ss test/input/1g7r.ss -rr test/input/1g7r.rr -o test/output/1g7r -mcount 20 -rrtype ca &> test/output/1g7r.log &

Expected Results:
PDB     TM-score        RMSD    MODEL
1eaz    0.74    2.68    ./test/output/1eaz/stage2/1eaz_12.pdb
1g7r    0.67    2.98    ./test/output/1g7r/stage1/1g7r_5.pdb
1guu    0.59    4.37    ./test/output/1guu/stage2/1guu_11.pdb
1qjp    0.64    4.30    ./test/output/1qjp/stage1/1qjp_1.pdb
1smx    0.63    3.35    ./test/output/1smx/stage2/1smx_20.pdb
1vjk    0.79    1.90    ./test/output/1vjk/stage1/1vjk_model2.pdb
helix   0.98    0.34    ./test/output/helix/stage1/helix_20.pdb

For a more rigorous testing, please download the script and inputs
for 150 proteins in fragfold benchmark set at
http://protein.rnet.missouri.edu/confold/tool.html
 
-----------------------------------------------------------
Please cite:
"CONFOLD: Residue-Residue Contact-guided ab initio Protein Folding",
Proteins: Structure, Function, and Bioinformatics, 2015.
B. Adhikari, D. Bhattacharya, R. Cao, J. Cheng. 

-----------------------------------------------------------
bap54@mail.missouri.edu (developer)
chengji@missouri.edu (PI)
