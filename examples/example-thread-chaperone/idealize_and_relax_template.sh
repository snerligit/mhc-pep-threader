#!/bin/bash

template=5opi

# add protons
score_jd2.linuxgccrelease -ignore_unrecognized_res -in:file:s "$template".pdb -out:pdb

# idealize
idealize_jd2.linuxgccrelease -ignore_unrecognized_res -in:file:s "$template"_0001.pdb -out:pdb

# relax
relax.linuxgccrelease -in:file:s "$template"_0001_0001.pdb -out:pdb

wait $!

#remove unneeded files
rm "$template".pdb
rm "$template"_0001.pdb
rm "$template"_0001_0001.pdb
mv "$template"_0001_0001_0001.pdb "$template"_idealized-relaxed.pdb
