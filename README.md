# Introduction
A repository to share our useful scripts in, between our group and others.

## Structure
Each member who would like to add some useful scripts here can start a new directory with their name and add inside their scripts and a useful README explaining their scripts. To keep this tidy, also update this README with a section name (your name), and under it a very quick and short explanation if what you have in your directory. Therefore, one can access the repository, see the general description and then go to that specific directory to get more info. Also make sure inputs and outputs are clear for each script you add, just to make it easy for others to share and use

## Topics

#### Alphafold
* run alphafold in parallel for multiple sequences\
  [Code](Roman/run_multifold.py) | [Details](Roman/README.md#Scripts)

#### PDB handling
* extract a specific chain from a pdb-file \
  [Code](Roman/chain.py) | [Details](Roman/README.md#Scripts)

### Fawaz
The following scripts I have added:
* [extract clusters](Fawaz/extract_clusters.py) that extract clusters produced from `mmseqs2` clustering into separate fasta file, i.e. each cluster is in a separate fasta file
* [separate proteins](Fawaz/separate_proteins.py) that extract all proteins in annotation files given into one big fasta file, the proteins separated will be named according to the original annotation file they came from, the gene/protein name/id and the coordinates.
* [fasta fastq stats](Fawaz/fasta_fastq_stats.sh) this is a very simple (kinda hacky) bash script that counts the number of reads in a fasta or fastq file (can also be gzipped), and the average length of the reads and the total number of sequences in that file
* [extract protein from patric data](Fawaz/extract_patric_protein.py) this script takes a patric assembly and patric tab annotation and a FIGfamily id (e.g. FIG00000080) and extracts the amino acid sequence corresponding to that gene if it exists in the .tab file given for that specific strain

### Ilya
* [PDB to graph](Ilya/pdb_to_pyg.py) that takes either single pdb file or a batch and converts them to torch_geometric.data.Data-like dictionary. Output will always be a pickle file containing either a single dict or a pandas DataFrame of them.
* [PDB to fasta](Ilya/pdb_to_fasta.py) that takes a directory of pdb files and puts all their sequences into a single fasta file.

### Roman

* [MultiFold](Roman/run_multifold.py) run alphafold in parallel for multiple sequences.
* [Chain from PDB](Roman/chain.py) extract a specific chain from a PDB file.
* [Crawl PubChem](Roman/pubchem_crawler.py) Download and scan all compounds in the PubChem database

See the [README](Roman/README.md) for detailed explanations

### Alper

* [GBK Parser](Alper/gbk_parser.py) parse gbk files and store them in list of records class.

### Aldo

* [Dataset_generator](Aldo/data_generator/dataset_generator.ipynb) This notebook is designed to format the information
from binding DB into pandas dataframe. In order to do it employs other useful scripts by its own.
Those other scripts are:
  * Computes molecular and chemical descriptor for all ligands in dataset
  * Download pdb structures from pandas df column within ids.
  * Clean non-protein molecules from pdb structures
  * Pdb2fasta (Ilya/pdb_to_fasta.py) is adapted to add fasta sequences innto pandas dataframe
  * Class balance of dataset within user-tune hyperparameters
  * Filtering of dataset by enzymatic activity (EC number)
  * Computes the RMSD and sequences identify of all pdb 3D structures pool.
* [Commands for RINDTI installation | june 2023](Aldo/rindti_installation/rindti_commands_for_installation.txt) This is a txt file that contains the commands that worked for
me to succesfully install RINDTI model at June 2023. 

See the [Requirements](Alper/requirements.md) for required installations
