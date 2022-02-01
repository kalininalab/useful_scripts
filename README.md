# Introduction
A repository to share our useful scripts in, between our group and others.

## Structure
Each member who would like to add some useful scripts here can start a new directory with their name and add inside their scripts and a useful README explaining their scripts. To keep this tidy, also update this README with a section name (your name), and under it a very quick and short explanation if what you have in your directory. Therefore, one can access the repository, see the general description and then go to that specific directory to get more info. Also make sure inputs and outputs are clear for each script you add, just to make it easy for others to share and use

### Fawaz
The following scripts I have added:
* [extract clusters](Fawaz/extract_clusters.py) that extract clusters produced from `mmseqs2` clustering into separate fasta file, i.e. each cluster is in a separate fasta file
* [separate proteins](Fawaz/separate_proteins.py) that extract all proteins in annotation files given into one big fasta file, the proteins separated will be named according to the original annotation file they came from, the gene/protein name/id and the coordinates.
