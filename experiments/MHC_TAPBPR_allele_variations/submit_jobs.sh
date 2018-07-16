#!/bin/bash

for filename in *.pdb; do
	qsub -pe orte 10 production.relax.job $filename
done
