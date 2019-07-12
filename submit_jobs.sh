#!/bin/bash

native="/home/snerli/test_rmhc/mhc-pep-threader/examples/example-thread-peptide/native.pdb"
cutpoint=180

for filename in /home/snerli/test_rmhc/mhc-pep-threader/examples/example-thread-peptide/*relaxed*; do
	name=`basename $filename .pdb`
	movemap_name=`echo $filename | sed 's|\(.*\)_.*_.*|\1|'`
	ext=".movemap"
	move=$movemap_name$ext
	qsub -pe orte 1 production.RosettaMHC.job $filename $move $native $cutpoint
done
