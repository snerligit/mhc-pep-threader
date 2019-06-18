# MHC-peptide threading protocol

## References

* "RosettaMHC; Structure-based modeling of neoantigen/HLA complexes” Nerli et al., manuscript in preparation

* Toor JS, Rao AA, McShan AC, Yarmarkovich M, Nerli S, Yamaguchi K, Madejska AA, Nguyen S, Tripathi S, Maris JM, Salama SR, Haussler D and Sgourakis NG (2018) A Recurrent Mutation in Anaplastic Lymphoma Kinase with Distinct Neoepitope Conformations. Front. Immunol. 9:99. doi: 10.3389/fimmu.2018.00099

## Introduction

RosettaMHC is an automated structure-based method used to study peptide/MHC-I biomolecular systems in a high-throughput manner. We can use this method to

* Model peptide/MHC-I structures and predict peptides that bind to MHC-I molecules,
* Identify altered peptide ligands, and
* Filter molecular chaperones that load or edit peptides in the MHC-I molecules.

RosettaMHC uses homology modeling and energy minimization to predict structures of peptide/MHC-I/(chaperone or TCR) molecules. From the modeled complexes, interface analyses is carried out which reports on the affinities between (i) peptide and MHC-I, (ii) peptide/MHC-I and TCR, or (iii) peptide/MHC-I and chaperone. To run this method, a user has to provide X-ray antigen/MHC-I/(chaperone or TCR) (or template) structure alongside a list of names of MHC-I molecules, antigen and chaperone (or TCR) sequences (See Input files section below). RosettaMHC automatically prepares the template, aligns template and target sequences, performs homology modeling and energy minimization, and predicts binding affinities (See Algorithm section below). RosettaMHC is implemented using PyRosetta software libraries and supports MPI-based parallelization.

We carried out benchmark calculations using RosettaMHC and found that (i) the models of complexes are typically within 2 Å from the reference structures (tested using X-ray structures in the PDB) (See Reference) and (ii) predicted non-self-antigen with HLA-A*01 (a sub-group of MHC-I) model was very close to the solved X-ray structure (See Reference).

Limitation: This method works when the sequences of antigens in both template and target are of same length.

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

* Clustal omega (command line interface), which is available for download [here](http://www.clustal.org/omega/). The binary as per the operating system must be downloaded and renamed to clustalo. The PATH for clustalo has to be added to bashrc or pash_profile by using export PATH="<Path to clsutalo>":$PATH. 

* [Rosetta](https://www.rosettacommons.org/software) (optional)

### Input Files

* A template PDB (required). This file is required to thread a target sequence onto it.

* MHC list file. This is the file that contains list of MHC molecules named in standard format. For example: A*02:01. The list of supported HLA sequences can be viewed using a command line shown under Examples section.

* Peptide list file. This is a file containing list of peptide sequences in fasta format.

* TCR sequence file. This is a file containing list of TCR sequences in fasta format.

* Chaperone list file. This is a file similar to MHC list file. Currently, we only have tapasin and tapbpr sequences in our database. If you don't want to provide chaperones, you need not provide the file or simply add "none" in the chaperone list file.

* beta2m list file. This is a file similar to MHC list file. Currently, we have sequences for human, mouse, chicken, cow and rat beta2m in our database. If you don't want to provide beta2m sequence, you need not provide the file or simply add "none" in the beta2m list file.

### Output Files

* Cleaned template PDB file. Subsequently, if a template is idealized and relaxed in the pipeline, the idealized and refined structure is also produced.

* Clustal omega input fasta file and corresponding alignment file.

* Target sequence fasta file.

* Rosetta format alignment file (grishin).

* [Movemap](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/movemap-file) file .

* Corresponding threaded PDB and refined files.

* A csv file containing binding affinity measure (as per Rosetta energy function).

## Basic options

    * -list_mhcs : Lists all the MHCs for which sequences are available in the database.

    * -template_pdb : Provide template structure in PDB to perform threading.

    * -mhcs : Provide the list of names of MHCs in the file, if you want to include all, just type "all" in the file.

    * -beta2m :Provide the file containing names of beta2m, Choices include: [humanbeta2m, mousebeta2m, chickenbeta2m, bovinebeta2m, ratbeta2m, none].

    * -peptides : Provide fasta file with peptide sequences that need to be threaded.

    * -tcr : Provide the file containing tcr sequence in fasta format.

    * -chaperone : Provide the file containing names of chaperones, Choices include: [tapasin, tapbpr, none].

    * -mhc_trim_length : Provide the number of residue from which the mhc should be trimmed, default value=181)

    * -no_trim_mhc : Should we model the whole complex.

    * -idealize_relax : idealize and relax template structure before threading.

    * -relax_after_threading : idealize and relax template structure before threading.

    * -mhc_chain : Provide mhc chain id in the template.

    * -peptide_chain : Provide peptide chain id in the template.

    * -pep_start_index : Provide peptide start index.

    * -groove_distance : Provide distance to select nearest groove residues from the peptide, default value=3.5 Å

    * -interface_cutpoint : Last residue index that separates the interfaces for which you are calculating binding energies, default value = 0.

    * -out_file : Output file name in csv format to write the binding energies", default name="binding_energies.csv".

    * -nstruct : number of times a threaded structure should be relaxed", type=int, default value=1.

## Examples

Below are the examples that showcase how RosettaMHC can be utilized to understand MHC-peptide/MHC-peptide-TCR/MHC-Chaperone systems.

### Modeling MHC-peptide complexes

Example scripts in run.sh or mpi_run.sh used to model and extract binding energies of peptides to MHCs can be found under examples/example-thread-peptide.

### Modeling MHC-peptide-TCR complexes

Example scripts in run.sh or mpi_run.sh used to model and extract binding energies of peptides-MHCs to TCR molecules can be found under examples/example-thread-tcr.

### Modeling MHC-Chaperone complexes

Example scripts in run.sh or mpi_run.sh used to model and extract binding energies of MHCs to chaperone molecules can be found under examples/example-thread-chaperone.
