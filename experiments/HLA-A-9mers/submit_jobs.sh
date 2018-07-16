#!/bin/bash

for filename in *.pdb; do
	name=`basename $filename .pdb`
	ext=".movemap"
	move=$name$ext
	qsub -pe orte 32 production.relax.job $filename $move
done
