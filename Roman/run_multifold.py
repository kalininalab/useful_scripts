import sys
from multiprocessing.pool import ThreadPool
import subprocess
from shutil import copyfile
import os
from urllib.request import urlretrieve
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser


class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accecpt_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        return 0


def download_pdb(pdb_code, output_file):
    if "_" in pdb_code:
        chain = pdb_code.split("_")[1]
        pdb_code = pdb_code.split("_")[0]
    else:
        chain = None
    url = "https://files.rcsb.org/download/" + pdb_code + ".pdb"
    try:
        urlretrieve(url, output_file)
        if chain is not None:
            p = PDBParser(PERMISSIVE=1)
            structure = p.get_structure(output_file, output_file)
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_file, ChainSelect(chain))
    except Exception as e:
        print(pdb_code, "->", e)


def postprocess(output_dir, name):
    if os.path.exists(os.path.join(output_dir, "tmp", name, "ranked_0.pdb")):
        copyfile(os.path.join(output_dir, "tmp", name, "ranked_0.pdb"), os.path.join(output_dir, name + ".pdb"))
    elif os.path.exists(os.path.join(output_dir, "tmp", name, "msas", "pdb_hits.hhr")):
        pdb_hits = open(os.path.join(output_dir, "tmp", name, "msas", "pdb_hits.hhr"))
        for i, line in enumerate(pdb_hits.readlines()):
            if i == 9:
                pdb = [x for x in line.split(" ") if len(x) > 0][1]
                download_pdb(pdb, os.path.join(output_dir, name + ".pdb"))
            elif i > 9:
                break


def run_single(data, file_name, output_dir, time, only_post=False):
    try:
        name = os.path.basename(file_name).split(".")[0]
        if not only_post:
            subprocess.run(
                ["run_alphafold.sh", "-d", data, "-f", file_name, "-o", os.path.join(output_dir, "tmp"), "-t", time,
                 "-g", "False"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        postprocess(output_dir, name)
    except Exception as e:
        print(file_name, "->", e)


def main():
    data, input_dir, output_dir, time = sys.argv[1:5]
    kernels = int(sys.argv[5])
    post_only = False if len(sys.argv) < 7 else sys.argv[6] == "-p"
    input_files = [filename for filename in os.listdir(input_dir) if filename.endswith(".fasta")]

    pool = ThreadPool(processes=kernels)
    tasks = [None for _ in range(len(input_files))]
    for i, file_name in enumerate(input_files):
        print("Start process", i)
        tasks[i] = pool.apply_async(run_single, (data, os.path.join(input_dir, file_name), output_dir, time, post_only))
    for i in range(len(tasks)):
        tasks[i].get()
        print("Collect process", i)
    print("Finished")


if __name__ == "__main__":
    main()
