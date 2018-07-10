# MHC peptide threading protocol

This protocol requires the following libraries

    * Python3.x

    * PyRosetta4

    * Biopython

    * scikit-bio

    * clustal omega (command line interface)
    
    * Rosetta

Command to run this protocol for the example input from within the example folder
```
python3 ../main.py -template_pdb 3mrm.pdb -peptide pep_list -mhcs mhc_list -mhc_chain A -peptide_chain P
```

This method will output threaded structure in the pdb format and a movemap file for each target of interest. These files can be used directly within Rosetta to carry out energy minimization using the following command:
```
<path_to_rosetta_binaries>/relax.linuxgccrelease -database <path_to_rosetta>/main/database -s threaded_structure.pdb -in:file:movemap threaded_structure.movemap -nstruct 10 
```

The above command will create 10 relaxed structures whose binding energies can be extracted using Rosetta's InterfaceAnalyzer protocol.  
```
<path_to_rosetta_binaries>/InterfaceAnalyzer.linuxgccrelease -database <path_to_rosetta>/main/database -s threaded_structure_relaxed.pdb -out:file:score_only threaded_structure.ia 
```
threaded_structure.ia file will consist of the binding energy value. Look for column 6 or 7 which says dG_separated.

CAVEAT: The relaxed structures will be in a single chain format. However, InterfaceAnalyzer expects two chains to analyze interface. This can be created using PyMOL. There are other alternatives. Please contact me if you are seeking other alternatives

If you use this work, please cite:

Toor JS, Rao AA, McShan AC, Yarmarkovich M, Nerli S, Yamaguchi K, Madejska AA, Nguyen S, Tripathi S, Maris JM, Salama SR, Haussler D and Sgourakis NG (2018) A Recurrent Mutation in Anaplastic Lymphoma Kinase with Distinct Neoepitope Conformations. Front. Immunol. 9:99. doi: 10.3389/fimmu.2018.00099
