# MHC peptide threading protocol

This protocol requires the following libraries

    * Python3.x

    * PyRosetta4

    * Biopython

    * scikit-bio

Command to run this protocol for the example input from within the example folder
```
python3 ../main.py -template_pdb 3mrm.pdb -peptide pep_list -mhcs mhc_list -mhc_chain A -peptide_chain P
```

If you use this work, please cite:

Toor JS, Rao AA, McShan AC, Yarmarkovich M, Nerli S, Yamaguchi K, Madejska AA, Nguyen S, Tripathi S, Maris JM, Salama SR, Haussler D and Sgourakis NG (2018) A Recurrent Mutation in Anaplastic Lymphoma Kinase with Distinct Neoepitope Conformations. Front. Immunol. 9:99. doi: 10.3389/fimmu.2018.00099
