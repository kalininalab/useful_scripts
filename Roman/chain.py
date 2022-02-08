from Bio.PDB import PDBParser, PDBIO
import sys


def extract_chain(pdb_name, pdb_file, search_chain, output_file):
	io = PDBIO()
	pdb = PDBParser().get_structure(pdb_name, pdb_file)

	for chain in pdb.get_chains():
		if chain.get_id() == search_chain:
			io.set_structure(chain)
			io.save(output_file)


if __name__ == "__main__":
	extract_chain(*sys.argv[1:5])
