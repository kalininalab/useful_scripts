import os
import sys

import pymol
from pymol import cmd


def convert(mol_folder, pdb_folder):
    os.makedirs(pdb_folder, exist_ok=True)
    pymol.finish_launching(['pymol', '-cq'])

    count = len(os.listdir(mol_folder))
    for i, file in enumerate(os.listdir(mol_folder)):
        if file.endswith('.mol2'):
            print(f"\r{i} / {count} | {file}", end="")
            cmd.set('pdb_use_ter_records', 1)  # Use TER records for chain breaks
            cmd.load(os.path.join(mol_folder, file), 'mol2_object')
            # Check if the cmd object has the h_add attribute
            add_hydrogen = hasattr(cmd, 'h_add')
            if not add_hydrogen:
                cmd.h_add('mol2_object')  # Add hydrogen atoms
            cmd.sort()  # Sort atoms
            cmd.save(os.path.join(pdb_folder, os.path.splitext(file)[0] + '.pdb'), 'mol2_object')

    # Quit PyMOL
    cmd.quit()


if __name__ == '__main__':
    convert(*sys.argv[1:3])
