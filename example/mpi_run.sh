mpiexec -np 4 python3 ../main.py -nstruct 2 -relax_after_threading -template_pdb 3mrm.pdb -mhcs mhc_list -peptides pep_list -mhc_chain A -peptide_chain P -pep_start_index 181 -interface_cutpoint 180
