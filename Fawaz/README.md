## Separate Proteins
This script takes a txt file with a list of annotation files and an output file name/path, the annotation files should be in gbf, gbff, gbk and as is or gzipped version of them. Then it separates all proteins in these annotation files and produces one big files with all the proteins, the protein sequences produced will be names with the annotation file they originated with then the protein id and the coordinates of that protein in its genome. Can be used if you have a lot of annotations and you want to extract a lot of proteins for clustering later or any other use.
Internally, it's simply a Protein class that stores important information about the protein, and it loops through all the files given, open then with `Bio.SeqIO` parsing function, and checks if it's a CDS, then it extracts the translated sequences (amino acids) and if it has a name or just an id and writes these results to the output file.

## Extract Clusters
This scripts is dedicated to separate clusters into individual FASTA files from `mmseqs2` output. It takes as input a list of cluster names txt file, where each cluster name is in one line, so just the unique cluster names that one can get from the TSV file produces with cluster names, then it takes the `mmseqs2` fasta output with all sequences and output directory. The script will only output the clusters in the txt file, so if the user only wants to output certain clusters as a separate FASTA file, can only include those cluseter names in the txt file.

## Fasta and Fastq Stats
Simple bash script that takes one arguemnt with the following allowed extentions (.fasta, .fa, .fastq, .fq) with or without .gz at the end, and number of reads, total length (sum of all reads), and average length of reads in the file
