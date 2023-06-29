""" What does this code does?
This notebook generates 'lig', 'inter', and 'prot' files
from TSV files retrieved from the Binding DB dataset. Those
files are properly formatted to be used at RINDTI model
https://rindti.readthedocs.io/en/master/index.html

The files (TSV, TXT, or images) will contain information 
about:
    Class balance: This code will downsample the dataset to easily-tune 
    a max radio of 1:3. Minimum samples per class is also easily-tune parameter.
    results/active-inactives_info.txt file generated will provide the class balance per target.
    results/class_balance.txt file will contain global information of dataset balance.
    This code will balance each target active-inactives radio, NOT just whole set.

    Molecular descriptors: This code will generate a image of molecular descriptor of 
    dataset distribution of all dataset after lipinski rules filtering.
    Also complete dataset information can be generated from this code as those molecular
    descriptors are generated from rdkit in this notebook.
    
    pbd2fasta: This code also contain codes that are usuful to get fasta sequences from pdbs
    and place that information in new cells in dataframe. This code could be reused in 
    future as its own as it based on: https://github.com/kalininalab/useful_scripts/tree/main/Ilya

    pdb-cleaning: This code also contains a cell code that is useful to clean pdbs from other
    non-protein molecules, as ligands, water or metals and get only protein molecules.

    pdb folder: This code also generated a folder with pdb structures and renames it as pdb id.pdb,
    same name that would be used in prot.tsv id.

    Active Threshold: This code also create the active-inactive classification based
    on a cutoff of 10 mM. Whereas less than 10 mM indicates positive activity.

    Filter by EC number: In case to be interest in specific class of enzymatic activity,
    it could be filter by EC number. However, this highly depends on pdb files, since
    that number is retrieved from the information contained in that pdb.ent file. If non-existent
    will not be parse to dataframe.

    Formatted dataset: The dataset is properly formatted to be used at RINDTI model
    https://rindti.readthedocs.io/en/master/index.html 
"""
