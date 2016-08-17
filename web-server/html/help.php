<?php include("header.html");?>

<tr bgcolor="#EEEEEE">
<td>
<p ><b>Job Submission</b></p>
<p >We have setup our web application such that you can check if the server is running at every step.</p>
<ul>
<li >After you submit your job, to make sure that your job was submitted, open the job log history <b><a href="http://protein.rnet.missouri.edu/confold/logs/history.log">here</a></b> and find your input at the end of the file. If you cannot open this file, something already has gone wrong. We may or may not be working on fixing the server, depending on our priorities.</li>
<li >The next step is to check the folder where you jobs are running. Your jobs folder can be browsed online at the location is shown in the page that is opened after the jobs you submit your job. The URL is also in the email that is sent to you. Refresh and check if files are being created and the folder contents seem to be updating.</li>
</ul>
</td>
</tr>
<tr>
<td>
<p ><b>Secondary Structure Related Parameters</b></p>
<p >
Lambda
<ul>
<li >Lambda is the parameter to loosen the secondary structure restraints.</li>
<li >The upper bounds and lower bounds of all secondary structure restraints are multiplied by lambda.</li>
</ul>
</p>

<p >
Sheet-detection threshold
<ul>
<li >How far to go to detect beta strand pairs in stage 1 model? (distance in Angstroms)</li>
<li >Selecting a lower value like 6.5 will pair really close strands.</li>
<li >Selecting a high value like 8.5 will pair many strands, possibly false positives.</li>
<li >Default value of 7.0 works best for true as well as predicted contacts.</li>
</ul>
</p>

<p >
Restraints Weight
<ul>
<li >This value controls the ratio of contact restraints and secondary structure restraints weight.</li>
<li >For example, selecting 0.5 sets half of contact restraints weight to secondary structure restraints.</li>
</ul>
</p>

<p >
Pairing Information :: Pairing File Format
</p>
<ul>
<li >6 columns a, b, c, d, t, and f in each row</li>
<li >a-b and c-d are residue strands, for example 2-7 and 20-25</li>
<li >t is the pairing type (A or P), and f is the confidence of pairing</li>
<li >a must always be less than b</li>
<li >c must be less than d if parallel and greater than d if anti-parallel</li>
</ul>
</td>
</tr>

<tr bgcolor="#EEEEEE">
<td>
<p ><b>Contact Related Parameters</b></p>
<p >Select top-xL Contacts</p>
<ul>
<li >How many contacts should be used as restraints?</li>
<li >For example, selecting top-0.4L will use only top-40 contacts if the sequence is 100 residues long.</li>
</ul>

<p >Contact Type</p>
<ul>
<li >Select cb is the input contacts are between Carbon-beta atoms of the residues in contact, otherwise select ca.</li>
</ul>

<p >Contact Restraints Weight</p>
<ul>
<li >contact restraints weight</li>
</ul>
</td>
</tr>
<tr><td>
<p ><b>CONFOLD and CNS Parameters</b></p>
<p >Run Stage2 Flag</p>
<ul>
<li >How to build models in the second stage?</li>
<ul>
<li >no stage 2 - don't run stage 2 at all (reduces the execution time by half)</li>
<li >sheet detect - detect beta-strand pairings but don't filter any unsatisfied contacts</li>
<li >contact filter - filter any unsatisfied contacts but don't detect strand pairings</li>
<li >both - filter any unsatisfied contacts and detect beta-pairings</li>
</ul>
</li>
</ul>

<p >Atom Selection Scheme</p>
<ul>
<li >What atoms to use for distance geometry?</li>
<li >'existing' refers to the existing list of atoms in the CNS suite.</li>
</ul>

<p >Second Repel Radius</p>
<ul>
<li >Second repel radius in the CNS DGSA script. Default value of 0.8 usually gives best results.</li>
</ul>
</td>
</tr>
<?php include("footer.html");?>



