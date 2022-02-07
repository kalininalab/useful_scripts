# Introduction
A repository to share our useful scripts in, between our group and others.

## Structure
Each member who would like to add some useful scripts here can start a new directory with their name and add inside their scripts and a useful README explaining their scripts. To keep this tidy, also update this README with a section name (your name), and under it a very quick and short explanation if what you have in your directory. Therefore, one can access the repository, see the general description and then go to that specific directory to get more info. Also make sure inputs and outputs are clear for each script you add, just to make it easy for others to share and use

### Fawaz
The following scripts I have added:
* [extract clusters](Fawaz/extract_clusters.py) that extract clusters produced from `mmseqs2` clustering into separate fasta file, i.e. each cluster is in a separate fasta file
* [separate proteins](Fawaz/separate_proteins.py) that extract all proteins in annotation files given into one big fasta file, the proteins separated will be named according to the original annotation file they came from, the gene/protein name/id and the coordinates.
* [fasta fastq stats](Fawaz/fasta_fastq_stats.sh) this is a very simple (kinda hacky) bash script that counts the number of reads in a fasta or fastq file (can also be gzipped), and the average length of the reads and the total number of sequences in that file

### Ilya
* [PDB to graph](Ilya/pdb_to_pyg.py) that takes either single pdb file or a batch and converts them to torch_geometric.data.Data-like dictionary. Output will always be a pickle file containing either a single dict or a pandas DataFrame of them.
