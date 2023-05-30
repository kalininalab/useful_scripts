import os
import sys

from joblib import Parallel, delayed
import pymol
from pymol import cmd
from tqdm import tqdm


def convert_single(mol2_file, pdb_folder):
    if not mol2_file.endswith('.mol2'):
        return
    pymol.finish_launching(['pymol', '-cq'])
    cmd.set('pdb_use_ter_records', 1)  # Use TER records for chain breaks
    cmd.load(mol2_file, 'mol2_object')
    # Check if the cmd object has the h_add attribute
    add_hydrogen = hasattr(cmd, 'h_add')
    if not add_hydrogen:
        cmd.h_add('mol2_object')  # Add hydrogen atoms
    cmd.sort()  # Sort atoms
    cmd.save(os.path.join(pdb_folder, mol2_file.split("/")[-1][:-4] + 'pdb'), 'mol2_object')
    cmd.delete("all")


def convert(mol_folder, pdb_folder, num_threads):
    os.makedirs(pdb_folder, exist_ok=True)
    Parallel(n_jobs=int(num_threads))(
        delayed(convert_single)(os.path.join(mol_folder, file), pdb_folder) for file in tqdm(os.listdir(mol_folder))
    )

    # Quit PyMOL
    cmd.quit()


if __name__ == '__main__':
    convert(*sys.argv[1:4])
