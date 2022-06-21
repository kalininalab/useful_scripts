import os.path
import sys
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.rdmolops import AddHs, GetAdjacencyMatrix
import pandas as pd
import numpy as np
import urllib.request
import gzip
import shutil


def download(num, log_file):
    """
    Download the SDF file of the compounds in the specified interval into a SCRATCH_SAS folder.

    Args:
        num (int): Lower number of the interval to be downloaded
        log_file (file-pointer): output stream to write the log to

    Returns:
        absolute filepath of the stored xyz.sdf.gz file
    """
    c_name = f"Compound_{num + 1:09d}_{num + 500_000:09d}.sdf"
    filename = f"/scratch/SCRATCH_SAS/roman/pubchem/" + c_name
    print(filename, file=log_file)
    if os.path.exists(filename):
        return filename
    urllib.request.urlretrieve(
        f"https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/" + c_name + ".gz",
        filename + ".gz"
    )
    # unzip it, this might fail for some files, then use zcat [filename.sdf.gz] > [filename.sdf]
    with gzip.open(filename + ".gz", 'rb') as f_in:
        with open(filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return filename


def convert(datafile):
    """
    Convert an SDF file into a csv file for better readability. The csv file will have the exact same name and location
    except for being a csv file with .csv ending.

    Args:
        datafile (str): filepath to the sdf file
    """
    if os.path.exists(datafile[:-3] + "csv"):
        return
    try:
        # just let RDKit do the magic and hope it doesn't crash
        tmp = PandasTools.LoadSDF(datafile)
        tmp.to_csv(datafile[:-3] + "csv")
    except Exception:
        return


def process_file(num, log_file=None):
    """
    Process the next file, i.e. download, extract, convert, and filter it. The filtered compounds will be stored in a
    results.tsv file in the working directory. By opening the file in the "appending"-mode, I ensure, all results will
    go there (also for all following and previous files, they will be in the same file).

    Args:
        num (int): Next interval-lower-bound to be processed.
        log_file (file-pointer): output stream to write the log to.
    """
    filename = download(num, log_file)
    if not os.path.exists(f"{filename[:-3]}.csv"):
        print("Convert SDF to CSV ...", file=log_file)
        convert(filename)
        print("\tConversion completed", file=log_file)
    filtering(read(f"{filename[:-3]}csv"), output=open("./result.tsv", "a"))


def read(datafile):
    """
    Read in the datafile and return a generator for every line's CID and ISO-SMILES as tuple.
    Args:
        datafile (str): filepath to the csv file to be read in line by line
    Yields:
        Every line's PubChem CID (Compound ID) and the Isomeric SMILES string of every compound
    """
    df = pd.read_csv(datafile)
    for index, row in df.iterrows():
        if index == 0:
            continue
        yield row["PUBCHEM_COMPOUND_CID"], row["PUBCHEM_OPENEYE_ISO_SMILES"]


def filtering(reads, output=None):
    """
    Filtering the compounds based on what read yields.
    The filtering is actually not commented as this is custom and everyone might apply different filters.

    Args:
        reads (Tuple[str, str]): Tuple of CID and ISO-SMILES string from Compound file.
        output (file-pointer): output stream to a file where to store the compounds in that passed this filter.
    """
    count = 0
    passed = 0
    for i, (cid, smiles) in enumerate(reads):
        print(f"\r{i:6}/500000", end="")
        try:
            count += 1

            if "+" in smiles or "-" in smiles or "." in smiles:
                continue

            mol = Chem.MolFromSmiles(smiles)
            ccount = sum([a.GetAtomicNum() == 6 for a in mol.GetAtoms()])
            ocount = sum([a.GetAtomicNum() == 8 for a in mol.GetAtoms()])
            rcount = sum([a.GetAtomicNum() not in [6, 8] for a in mol.GetAtoms()])
            rings = mol.GetRingInfo().AtomRings()
            tmp = [l for sublist in rings for l in sublist]

            if not (
                ccount * 0.5 < ocount < ccount * 2
                and (ccount + ocount) / 4 > rcount
                and len(rings) > 0
                and all([len(r) in [5, 6] for r in rings])
                and all([sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in r) == 1 for r in rings])
                and all([sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in r) == (len(r) - 1) for r in rings])
                and len(set(tmp)) == len(tmp)
            ):
                continue

            adj = np.array(GetAdjacencyMatrix(mol))
            adj_5 = np.linalg.matrix_power(adj, 5)
            tmp = set(tmp)

            if not (
                all([len(set(np.nonzero(adj_5[:, a.GetIdx()])[0].flatten()).intersection(tmp)) != 0 for a in mol.GetAtoms()])
            ):
                continue

            # After passing all filters, save the compound and increase the statistics
            passed += 1
            print(cid, "\t", smiles, file=output)
        except Exception as e:
            raise e
    print(f"{passed}/{count}: {passed / count:.5}")


def main():
    """
    Iterate over all Compound file of the PubChem database and filter them according to some filters.
    These filters might be implemented in the filtering-method according to the given comments/description.

    !!! ATTENTION: THE SCRIPT WILL DOWNLOAD AND STORE THE EXTRACTED PUBCHEM DATABASE ON YOUR HARDDRIVE
    (e.g. scratch/SCRATCH_SAS/.../pubchem). THIS WILL NEED APPROXIMATELY 2 TB!!! (21.06.2022)
    MAKE SURE, YOUR TARGET DIRECTORY HAS ENOUGH SPACE!!!
    """
    log_file = open("log.txt", "w")
    i = 0
    while i < 163_000_001:
        try:
            process_file(i)
        except Exception as e:
            print("Failure with", i, file=log_file)
            print(e, file=log_file)
        i += 500_000


if __name__ == '__main__':
    main()
