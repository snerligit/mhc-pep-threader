# MHC-peptide threading protocol

## References

* Toor JS, Rao AA, McShan AC, Yarmarkovich M, Nerli S, Yamaguchi K, Madejska AA, Nguyen S, Tripathi S, Maris JM, Salama SR, Haussler D and Sgourakis NG (2018) A Recurrent Mutation in Anaplastic Lymphoma Kinase with Distinct Neoepitope Conformations. Front. Immunol. 9:99. doi: 10.3389/fimmu.2018.00099

## Introduction


## Algorithm

RosettaMHC has following stages:

* Treat template structure: This stage takes a given template structure in PDB format, cleans the PDB file, idealizes and refines it using respective Rosetta methods. This stage is optional. It can be performed using RosettaMHC protocol pipeline or performed individually using Rosetta binaries.

* Align template and target: This phase aligns sequences of given targets with template structure using Clustal Omega command line interface (See Software Requirements on how to obtain Clustal Omega). The alignment is stored in a Rosetta compatible alignment file called [grishin](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Grishan-format-alignment) file.

* Threading: Here, target sequence is threaded onto the template structure using partial_thread protocol in Rosetta.

* Refinement: The threaded structure is further refined using Rosetta's refinement protocol (called relax).

* Calculating binding affinities: Binding affinities can be calculated between (i) peptides and MHC molecules, (ii) peptide-MHC and TCR molecules, and (iii) MHC molecules and chaperones using InterfaceAnalyzer protocol in Rosetta.

We currently provide a database of MHC, beta2m and chaperone molecules. Each database file consists of following number of sequences. These sequences are obtained from the following sources:

MHC - 2895 sequences of length 180 amino acids obtained from EBI.
MHC - 631 full-length sequences.
beta2m - 5 full-length sequences.
Chaperone - 2 full-length sequences.

## How to run this protocol?

### Software requirements
This protocol requires the following libraries/binaries:

    * Python3.x - Python3 is available for download [here](https://www.python.org/downloads/). Detailed instructions to install python for Linux and MacOS machines are available at [Guide for Linux machines](https://docs.python-guide.org/starting/install3/linux/) and [Guide for MacOS machines](https://docs.python-guide.org/starting/install3/osx/).

    * [PyRosetta4](http://www.pyrosetta.org/dow)

    * [Biopython](https://biopython.org/wiki/Download)

    * Clustal omega (command line interface), which is available for download [here](http://www.clustal.org/omega/).

    * [Rosetta](https://www.rosettacommons.org/software) (optional)

### Input Files

* A template PDB (required). This file is required to thread a target sequence onto it.

* MHC list file. This is the file that contains list of MHC molecules named in standard format. For example: A*02:01. The list of supported HLA sequences can be viewed using a command line shown under Examples section.

* Peptide list file. This is a file containing list of peptide sequences in fasta format.

* TCR sequence file. This is a file containing list of TCR sequences in fasta format.

* Chaperone list file. This is a file similar to MHC list file. Currently, we only have tapasin and tapbpr sequences in our database. If you don't want to provide chaperones, you need not provide the file or simply add "none" in the chaperone list file.

* beta2m list file. This is a file similar to MHC list file. Currently, we have sequences for human, mouse, chicken, cow and rat beta2m in our database. If you don't want to provide beta2m sequence, you need not provide the file or simply add "none" in the beta2m list file.

### Output Files


## Basic options

## Examples

Below are the examples that showcase how RosettaMHC can be utilized to understand MHC-peptide/MHC-peptide-TCR/MHC-Chaperone systems.

### Modeling MHC-peptide complexes

### Modeling MHC-peptide-TCR complexes

### Modeling MHC-Chaperone complexes
