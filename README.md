# MHC peptide threading protocol

This protocol requires the following libraries

    * Python3.x

    * PyRosetta4

    * Biopython

    * scikit-bio

Command to run this protocol for the example input
```
python3 main.py -template_pdb example/template_01.pdb -fasta example/target.fasta -pep_start_index 180 -nstruct 1
```

If you use this work, please cite:

Toor JS, Rao AA, McShan AC, Yarmarkovich M, Nerli S, Yamaguchi K, Madejska AA, Nguyen S, Tripathi S, Maris JM, Salama SR, Haussler D and Sgourakis NG (2018) A Recurrent Mutation in Anaplastic Lymphoma Kinase with Distinct Neoepitope Conformations. Front. Immunol. 9:99. doi: 10.3389/fimmu.2018.00099
