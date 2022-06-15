import os.path

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.rdmolops import AddHs
import pandas as pd
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
            stats = {
                "ccount":sum([a.GetAtomicNum() == 6 for a in mol.GetAtoms()]),
                "ocount":sum([a.GetAtomicNum() == 8 for a in mol.GetAtoms()]),
                "rcount":sum([a.GetAtomicNum() not in [6, 8] for a in mol.GetAtoms()]),
                "rings": mol.GetRingInfo().AtomRings(),
            }
            if not (stats["ccount"] <
            if len(stats["rings"]) > 0:
                stats["ring_atoms"] = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in stats["rings"][0]]
                stats["ring_bonds"] = sum([
                    b.GetBondType() == Chem.rdchem.BondType.SINGLE for b in mol.GetBonds()
                    if b.GetBeginAtomIdx() in stats["rings"][0] and b.GetEndAtomIdx() in stats["rings"][0]
                ]),
            if all(5 <= len(r) <= 6 for r in stats["ring_atoms"]) \
                    and all(stats["ring_atoms"]) == stats["ring_bonds"][0] \
                    and sum([x == 8 for x in stats["ring_atoms"]]) == 1 \
                    and sum([x == 6 for x in stats["ring_atoms"]]) == len(stats["ring_atoms"]) - 1:
                passed += 1
                print(cid, "\t", smiles, file=output)  # , "|", Chem.rdmolops.GetFormalCharge(mol))
        except Exception:
            pass
    print(f"{passed}/{count}: {passed / count:.5}", file=output)


def next_file(num, log_file=None):
    filename = download(num, log_file)
    if not os.path.exists(f"{filename[:-3]}.csv"):
        print("Convert SDF to CSV ...", file=log_file)
        convert(filename)
        print("\tConversion completed", file=log_file)
    filtering(read(f"{filename[:-3]}csv"), output=open("./result.tsv", "a"))


def main():
    log_file = open("log.txt", "w")
    i = 1_000_000
    while i < 163000001:
        try:
            next_file(i)
        except Exception as e:
            print("Failure with", i, file=log_file)
            print(e, file=log_file)
        i += 500_000


if __name__ == '__main__':
    main()

    # filename = "PubChem_compound_text_aspirin_records"
    # filename = "Compound_000000001_000500000"
    # if not os.path.exists(f"{filename}.csv"):
    #     print("Convert SDF to CSV ...")
    #     convert(f"{filename}.sdf")
    #     print("\tConversion completed")
    # filtering(read(f"{filename}.csv"), output=open("./result.tsv", "w"))
