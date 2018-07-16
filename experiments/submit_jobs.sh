#!/bin/bash

for filename in *_with_*.pdb; do
	name=`basename $filename .pdb`
	ext=".movemap"
	move=$name$ext
	qsub -pe orte 5 production.relax.job $filename $move
done
