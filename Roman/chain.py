from Bio.PDB import PDBParser, PDBIO
import os
from urllib.request import urlretrieve


def extract_chain(pdb_name, pdb_file, search_chain, output_file):
	io = PDBIO()
	pdb = PDBParser().get_structure(pdb_name, pdb_file)

	for chain in pdb.get_chains():
		if chain.get_id() == search_chain:
			io.set_structure(chain)
			io.save(output_file)


def download_pdb(pdb_code, output_file):
	if "_" in pdb_code:
		chain = pdb_code.split("_")[1]
		pdb_code = pdb_code.split("_")[0]
	else:
		return
	url = "https://files.rcsb.org/download/" + pdb_code + ".pdb"
	try:
		urlretrieve(url, output_file)
		extract_chain(pdb_code, output_file, chain, output_file)
	except Exception as e:
		print(pdb_code, "->", e)


def main():
	root_dir = "./outdated/Lectins/datasets/resources/oracle/structures/tmp/"
	target_dir = "./structures"
	for j, folder in enumerate(sorted(os.listdir(root_dir))):
		print(folder)
		if folder + ".pdb" in os.listdir(target_dir):
			continue
		content = os.listdir(os.path.join(root_dir, folder))
		if "msas" in content:
			if os.path.exists(os.path.join(root_dir, folder, "msas", "pdb_hits.hhr")):
				pdb_hits = open(os.path.join(root_dir, folder, "msas", "pdb_hits.hhr"))
				for i, line in enumerate(pdb_hits.readlines()):
					if i == 9:
						pdb = [x for x in line.split(" ") if len(x) > 0][1]
						download_pdb(pdb, os.path.join("./structures", folder + ".pdb"))
					elif i > 9:
						break
			elif os.path.exists(os.path.join(root_dir, folder, "ranked_0.pdb")):
				copyfile(os.path.join(root_dir, folder, "ranked_0.pdb"), os.path.join(target_dir, folder + ".pdb"))


if __name__ == "__main__":
	main()

