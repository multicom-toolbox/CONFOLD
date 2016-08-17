<?php include("header.html");?>
<form method="post" action="http://protein.rnet.missouri.edu/cgi-bin/confold/main.cgi">
		<tr bgcolor="#EEEEEE">
			<td align="right" width="30%"><b>E-mail Address </b></td>
			<td><input id="email" name="email" size="30" value="" required/></td>
		</tr>
		<tr bgcolor="#EEEEEE">
			<td align="right"><b>Job Id </b></td>
			<td><input id="id" name="id" size="30" value=""/></td>
		</tr>
		<tr bgcolor="#EEEEEE">
			<td valign="top" align="right"><b>Sequence </b></td>
			<td><textarea id="protein_sequence" name="protein_sequence" rows="5" cols="80" required></textarea></td>
			
		</tr >
		<tr bgcolor="#EEFFFF"><td></br></td><td></td></tr>
		<tr bgcolor="#EEFFFF">
			<td valign="top" align="right"><b>Secondary Structure</b></td>
			<td><textarea id="protein_sec" name="protein_sec" rows="5" cols="80"></textarea></td>
		</tr>
		<tr bgcolor="#EEFFFF"><td></td><td><a href="#" id="toggle" onClick="toggle_it('lambda');toggle_it('sheet');toggle_it('sswt');toggle_it('pair')">more SS options</a></td></tr>
		
		<tr bgcolor="#EEFFFF" id="lambda" style="display:none">
			<td align="right" title="Lambda is the parameter to loosen the secondary structure restraints. 
The upper bounds and lower bounds of all secondary structure restraints are multiplied by lambda."><b>lambda</b></td>
			<td><select id="lambda" name="lambda">
				<option value="0.5">0.5</option>
				<option value="0.6">0.6</option>
				<option value="0.7">0.7</option>
				<option value="0.8">0.8</option>
				<option value="0.9">0.9</option>
				<option value="1.0">1.0</option>
			</select></td>
		</tr>
		<tr bgcolor="#EEFFFF" id="sheet" style="display:none">
			<td align="right" title="How far to go to detect beta strand pairs in stage 1 model? (distance in Angstroms)
Selecting a lower value like 6.5 will pair really close strands.
Selecting a high value like 8.5 will pair many strands, possibly false positives.
Default value of 7.0 works best for true as well as predicted contacts."><b>sheet detection threshold</b></td>
			<td><select id="sheet_threshold" name="sheet_threshold">
				<option value="7.0">7.0</option>
				<option value="6.5">6.5</option>
				<option value="7.5">7.5</option>
				<option value="8.0">8.0</option>
				<option value="8.5">8.5</option>
			</select></td>
		</tr>
		<tr bgcolor="#EEFFFF" id="sswt" style="display:none">
			<td align="right" title="This value controls the ratio of contact restraints and secondary structure restraints weight.
For example, selecting 0.5 sets half of contact restraints weight to secondary structure restraints."><b>restraints weight</b></td>
			<td><select id="sec_wt" name="sec_wt">
				<option value="1">1</option>
				<option value="0.5">0.5</option>
				<option value="0.1">0.1</option>
				<option value="5">5</option>
			</select></td>
		</tr>
		<tr bgcolor="#EEFFFF" id="pair" style="display:none">
			<td valign="top" title="Pairing File Format:
- 6 column file with columns a, b, c, d, t, and f
- a-b and c-d are residue strands, for example 2-7 and 20-25
- t is the pairing type (A or P), and f is the confidence of pairing
- a must always be less than b
- c must be less than d if parallel and greater than d if anti-parallel" align="right"><b>pairing information</b></td>
			<td><textarea id="pairing" name="pairing" rows="8" cols="80" placeholder="Pairing File Format:
- 6 column file with columns a, b, c, d, t, and f
- a-b and c-d are residue strands, for example 2-7 and 20-25
- t is the pairing type (A or P), and f is the confidence of pairing
- a must always be less than b
- c must be less than d if parallel and greater than d if anti-parallel"></textarea></td>
		</tr>
		<tr>
		<tr bgcolor="#EEDDDD"><td></br></td><td></td></tr>
		<tr bgcolor="#EEDDDD">
			<td valign="top" title="Contacts in CASP RR format. Must be sorted by confidence (highest conf on top)" align="right"><b>Contacts</b> </td>
			<td><textarea id="rr" name="rr" rows="5" cols="80" required></textarea></td>
		</tr>
		<tr bgcolor="#EEDDDD"><td></td><td><a href="#" id="toggle" onClick="toggle_it('top');toggle_it('contype');toggle_it('conwt')">more RR options</a></td></tr>
		<tr  bgcolor="#EEDDDD"  id="top" style="display:none">
			<td align="right" title="How many contacts should be used as restraints?
For example, selecting top-0.4L will use only top-40 contacts if the sequence is 100 residues long." ><b>select top-xL contacts</b></td>
			<td><select id="rr_subset" name="rr_subset">
				<option value="all">all</option>
				<option value="0.4">top-0.4L</option>
				<option value="0.6">top-0.6L</option>
				<option value="0.8">top-0.8L</option>
				<option value="1.0">top-1.0L</option>
				<option value="1.2">top-1.2L</option>
				<option value="1.4">top-1.4L</option>
				<option value="1.6">top-1.6L</option>
				<option value="1.8">top-1.8L</option>
				<option value="2.0">top-2.0L</option>
				<option value="2.2">top-2.2L</option>
			</select></td>
		</tr>
		<tr  bgcolor="#EEDDDD"  id="contype" style="display:none">
			<td align="right" title="select cb is the input contacts are between Carbon-beta atoms of the residues in contact, otherwise select ca"><b>contact type</b></td>
			<td><select id="rr_type" name="rr_type">
				<option value="cb">cb</option>
				<option value="ca">ca</option>
			</select></td>
		</tr>
		<tr bgcolor="#EEDDDD" id="conwt" style="display:none">
			<td align="right"><b>contact restraints weight</b></td>
			<td><select id="con_wt" name="con_wt">
				<option value="10">10</option>
				<option value="1">1</option>
				<option value="5">5</option>
				<option value="50">50</option>
			</select></td>
		</tr>
		<tr bgcolor="#CCDDDD"><td></br></td><td></td></tr>
		<tr bgcolor="#CCDDDD"><td></td><td><a href="#" id="toggle3" onClick="toggle_it('stage2');toggle_it('scheme');toggle_it('repel')">more options</a></td></tr>
		<tr bgcolor="#CCDDDD" id="stage2" style="display:none">
			<td align="right" title="How to build models in the second stage?
no stage 2 - don't run stage 2 at all (reduces the execution time by half)
sheet detect - detect beta-strand pairings but don't filter any unsatisfied contacts
contact filter - filter any unsatisfied contacts but don't detect strand pairings
both - filter any unsatisfied contacts and detect beta-pairings
"><b>run stage2</b></td>
			<td><select id="stage2" name="stage2">
				<option value="1">no stage 2</option>
				<option value="2">sheet detect only</option>
				<option value="3">contact filter only</option>
				<option value="4">sheet detect and contact filter</option>
			</select></td>
		</tr>
		<tr bgcolor="#CCDDDD" id="scheme" style="display:none">
			<td align="right" title="What atoms to use for distance geometry?
'existing' refers to the existing list of atoms in the CNS suite."><b>atom selection scheme</b></td>
			<td><select id="atom_scheme" name="atom_scheme">
				<option value="2">existing, o</option>
				<option value="1">existing (ca, ha, n, hn, c, cb, cg)</option>
				<option value="3">existing, o, h</option>
				<option value="4">backbone atoms (c, ca, n, o)</option>
				<option value="5">backbone atoms, cb</option>
				<option value="6">backbone atoms, cb, h</option>
				<option value="7">backbone atoms, cb, h, cg</option>
			</select></td>
		</tr>
		<tr bgcolor="#CCDDDD" id="repel" style="display:none">
			<td align="right" title="Second repel radius in the CNS DGSA script. Default value of 0.8 usually gives best results."><b>second repel radius</b></td>
			<td><select id="rep2" name="rep2">
				<option value="0.8">0.8</option>
				<option value="0.85">0.85</option>
			</select></td>
		</tr>

			
		</tr>
		<tr>
			<td></td><td style="text-align: right; color: blue; font-size : 90%;"><a onclick="clearForm()">reset-fields</a> &nbsp; <a onclick="fillExample1()">example1(test)</a> &nbsp; <a onclick="fillExample2()">example2(true-rr)</a> &nbsp;<a onclick="fillExample3()">example3(predicted-rr)</a>&nbsp;</td>
		</tr>

</table>
<p align="center">
	<input name="submit"  type="submit" value="Run Job"/>
</p>
</form>
<?php include("footer.html");?>
