import os
import sys
import gzip
import pdb
from Bio import SeqIO


class Protein:
    def __init__(self, name, p_id, product, start, end, translation, parent_dir):
        self.name = name
        self.id = p_id
        self.product = product
        self.start = start
        self.end = end
        self.translation = translation
        self.parent_dir = parent_dir

    def return_fasta(self):
        """
        generate FASTA sequence with a sequence name
        """
        seq_name = [self.parent_dir]
        # I can also add 'product' to this list but product was causing problems
        # because it has special characters sometimes
        for k in ['name', 'id']:
            if getattr(self, k):
                seq_name.append(getattr(self, k))
        seq_name.append("coord_" + str(self.start) + "-" + str(self.end))

        return ">" + "_".join(seq_name) + "\n" + self.translation + "\n"


def write_sequences(proteins, prot_name, out_file):
    """
    takes the proteins dict and the protein name
    outputs the sequences for that protein in the out_file in fasta format
    """
    with open(out_file, "w") as out_file:
        if prot_name in proteins:
            for prot in proteins[prot_name]:
                out_file.write(prot.return_fasta())


# check all gbf and gbff files in directory given
if len(sys.argv) < 3:
    print("Takes a txt file with a list of gbk files and an output file name/path")
    print("accepted file types (gbf, gbff, gbk) or gzipped version of these")
    sys.exit()

in_list = sys.argv[1]
out_file = sys.argv[2]
gene_files = []

with open(in_list, "r") as infile:
    for f in infile:
        f = f.strip()
        gene_files.append(f)

"""
I will loop through each file store the features and check if they are CDS
if so, then I take the qualifiers and take the "gene" key for the name
and the "translation" key for the amino acid sequence
I should keep this information for the name of the sequences
file name the sequence came from, gene name, location (start, end)

to get start and end, I take the feature.location.nofuzzy_start and nofuzzy_end
"""

# proteins = dict()
outfile = open(out_file, "w")

for f in gene_files:
    # if gzipped then opening in text mode then giving the handle to seqio
    if f.endswith("gz"):
        handle = gzip.open(f, "rt")
    else:
        handle = open(f, "r")

    for rec in SeqIO.parse(handle, 'genbank'):
        for feature in rec.features:
            if feature.type == "CDS":

                # If the script was running from somewhere else, I take the name of the directory
                # where the genbank files are and use it as parent_dir
                # if the scrip was ran in the same directory the I get the current working directoy
                # and take the last part
                splits = f.split(os.sep)
                if (len(f.split(os.sep)) == 2) and (splits[0] == "."):
                    parent_dir = os.getcwd().split(os.sep)[-1]
                else:
                    parent_dir = f.split(os.sep)[-2]

                # if "gene" not in feature.qualifiers:
                if "translation" in feature.qualifiers:
                    try:
                        name = feature.qualifiers['gene'][0]
                    except KeyError:
                        name = ""

                    try:
                        p_id = feature.qualifiers['protein_id'][0]
                    except KeyError:
                        p_id = ""

                    try:
                        product_name = feature.qualifiers['product'][0]
                        product_name = product_name.replace(" ", "-")
                    except KeyError:
                        product_name = ""
                    # pdb.set_trace()
                    prot = Protein(name=name, p_id=p_id, product=product_name,
                                   start=feature.location.nofuzzy_start, end=feature.location.nofuzzy_end,
                                   translation=feature.qualifiers["translation"][0], parent_dir=parent_dir)
                    outfile.write(prot.return_fasta())
                    # if prot.id not in proteins:
                    #     proteins[prot.id] = [prot]
                    # else:
                    #     proteins[prot.id].append(prot)
                # else:
                #     if "translation" in feature.qualifiers:
                #
                #         prot = Protein(name=feature.qualifiers["gene"][0], p_id=feature.qualifiers["protein_id"][0],
                #                        start=feature.location.nofuzzy_start, end=feature.location.nofuzzy_end,
                #                        translation=feature.qualifiers["translation"][0], parent_dir=parent_dir)
                #         outfile.write(prot.return_fasta())
                #         # if prot.name not in proteins:
                #         #     proteins[prot.name] = [prot]
                #         # else:
                #         #     proteins[prot.name].append(prot)

outfile.close()
## TESTING, need to delete later
# random_prots = ['fldA', 'rpoA', 'yoaE', 'ymgG', 'ydiK', 'ybeL', 'cdsA', 'yqhC', 'rsxG', 'arcB', 'alkA', 'ycaL']
# for prot in random_prots:
#     out_f_name = "e_coli_" + prot + ".fasta"
#     with open(out_f_name, "w") as out_file:
#         for protein in proteins[prot]:
#             out_file.write(protein.return_fasta())

# I might add a module that calls clustalo on each outputted protein to generate an MSA and maybe remove the original
# to save space
# $clustalo -i all_proteins.fasta -o all_proteins.aln --threads 10 --distmat-out=all_proteins.mat --force --full --percent-id
