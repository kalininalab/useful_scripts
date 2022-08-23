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
    filename = f"/scratch/SCRATCH_SAS/roman/pubchem/Compound_{num + 1}_{num + 500_000}.sdf"
    print(filename, file=log_file)
    if os.path.exists(filename):
        return filename
    urllib.request.urlretrieve(
        f"https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_{num + 1:09d}_{num + 500_000:09d}.sdf.gz",
        filename + ".gz"
    )
    with gzip.open(filename + ".gz", 'rb') as f_in:
        with open(filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return filename


def convert(datafile):
    if os.path.exists(datafile[:-3] + "csv"):
        return
    try:
        tmp = PandasTools.LoadSDF(datafile)
        tmp.to_csv(datafile[:-3] + "csv")
    except Exception:
        return


def read(datafile):
    df = pd.read_csv(datafile)
    for index, row in df.iterrows():
        if index == 0:
            continue
        yield row["PUBCHEM_COMPOUND_CID"], row["PUBCHEM_OPENEYE_ISO_SMILES"]


def filtering(reads, output=None):
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
            adj_2 = np.linalg.matrix_power(adj, 2)
            adj_5 = np.linalg.matrix_power(adj, 5)
            tmp = set(tmp)

            if not (
                all([len(set(np.nonzero(adj_5[:, a.GetIdx()])[0].flatten()).intersection(tmp)) != 0 for a in mol.GetAtoms()])
            ):
                continue

            passed += 1
            print(cid, "\t", smiles, file=output)
        except Exception as e:
            raise e
    print(f"{passed}/{count}: {passed / count:.5}", file=output)


def next_file(num, log_file=None):
    filename = download(num, log_file)
    if not os.path.exists(f"{filename[:-3]}.csv"):
        print("Convert SDF to CSV ...", file=log_file)
        convert(filename)
        print("\tConversion completed", file=log_file)
    filtering(read(f"{filename[:-3]}csv"), output=open("./result_poly.tsv", "a"))


def main():
    log_file = open("log.txt", "w")
    i = 0
    while i < 163_000_001:
        try:
            next_file(i)
        except Exception as e:
            print("Failure with", i, file=log_file)
            print(e, file=log_file)
        i += 500_000


if __name__ == '__main__':
    # convert("PubChem_compound_polymer.sdf")
    # filtering(read("PubChem_compound_polymer.csv"), open("./result.tsv", "w"))
    main()

    # filename = "PubChem_compound_text_aspirin_records"
    # filename = "Compound_000000001_000500000"
    # if not os.path.exists(f"{filename}.csv"):
    #     print("Convert SDF to CSV ...")
    #     convert(f"{filename}.sdf")
    #     print("\tConversion completed")
    # filtering(read(f"{filename}.csv"), output=open("./result.tsv", "w"))
