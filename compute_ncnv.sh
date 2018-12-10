#!/bin/sh

state=$1 # protomer or oligomer
filename=$2
bcl="/hd0/lib14/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe protein:PDBConvert"
database=/dors/meilerlab/apps/PISCES/opm/2014

if [ $state == "protomer" ]; then
	echo "You requested to compute protomeric neighbor counts and neighbor vectors"
	while read -r line
	do
		echo "PDB ID read from file: $line"
		pdbid="$line"
		midid=${line:1:2}
		chainid=${line:4:1}
		# check to see if contact numbers have been computed
		if [ -s "$database/$midid/$pdbid".proto_ncnv ]; then
			echo "$database/$midid/$pdbid.proto_ncnv already exists, skip to the next protein"
		else    
			echo "Compute protomeric neighbor counts and neighbor vectors for $pdbid"        
			$bcl $database/$midid/$pdbid.pdb -sasa 3 -output_prefix $database/$midid/$pdbid.proto.
			# rename the ncnv file
			echo "rename ncnv file: mv $database/$midid/$pdbid.proto.$chainid.ncnv $database/$midid/$pdbid.proto_ncnv"
			mv $database/$midid/$pdbid.proto.$chainid.ncnv $database/$midid/$pdbid.proto_ncnv 
		fi
	done < $filename
else
	echo "You requested to compute oligomeric neighbor counts and neighbor vectors"
	while read -r line
	do
		echo "PDB ID read from file: $line"
		# use four-letter id instead
		pdbid=${line:0:4}
		midid=${line:1:2}
		chainid=${line:4:1}
		# check to see if contact numbers have been computed
		if [ -s "$database/$midid/$pdbid$chainid".oligo_ncnv ]; then
			echo "$database/$midid/$pdbid$chainid.oligo_ncnv already exists, skip to the next protein"
		else    
			echo "Compute oligomeric neighbor counts and neighbor vectors for $pdbid"        
			$bcl $database/$midid/$pdbid.bcl.pdb -sasa 3 -output_prefix $database/$midid/$pdbid.oligo.
			# rename the ncnv file
			echo "rename ncnv file: mv $database/$midid/$pdbid.oligo.$chainid.ncnv $database/$midid/$pdbid$chainid.oligo_ncnv"
			mv $database/$midid/$pdbid.oligo.$chainid.ncnv $database/$midid/$pdbid$chainid.oligo_ncnv 
		fi
	done < $filename
fi
