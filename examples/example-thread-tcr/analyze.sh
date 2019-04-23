#!/bin/bash

echo extracting relaxed structures
# make a directory for relaxed structures
mkdir relaxed_structures

# go into the directory 
cd relaxed_structures

# extract relaxed pdb from the files with the name .pdb_decoys_relaxed.out 
extract_pdbs.linuxgccrelease -in:file:silent ../*.pdb_decoys_relaxed.out -out:pdb

# go up a directory
cd ..

echo calculating interface energy
# make a directory for interface energy
mkdir interface_energies

# go into the directory
cd interface_energies

# analyze interface energies for peptide-MHC structures
# set cutpoint to the last mhc residue before the peptide it will always be 375
/home/snerli/rosetta/Rosetta_ref2015/main/source/bin/InterfaceAnalyzer.linuxgccrelease -database /home/snerli/rosetta/Rosetta_ref2015/main/database -in:file:silent ../*.out -cutpoint 375 -out:file:score_only interface_score.sc

# extract interface scores from the  dG_separated score column 
awk 'FNR > 2 {print $7"\t"$NF}' interface_score.sc > results_interface_score.txt

# go up a directoy
cd ..

echo done

