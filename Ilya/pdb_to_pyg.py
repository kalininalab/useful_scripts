""""""

import torch

def list_to_dict(l):
    """Convert list to dict"""
    return {val: i for i, val in enumerate(l)}


def onehot_encode(position: int, count: int) -> list:
    """One-hot encode position
    Args:
        position (int): Which entry to set to 1
        count (int): Max number of entries.
    Returns:
        list: list with zeroes and 1 in <position>
    """
    t = [0] * (count)
    t[position - 1] = 1
    return t


node_encoding = list_to_dict(
    [
        "ala",
        "arg",
        "asn",
        "asp",
        "cys",
        "gln",
        "glu",
        "gly",
        "his",
        "ile",
        "leu",
        "lys",
        "met",
        "phe",
        "pro",
        "ser",
        "thr",
        "trp",
        "tyr",
        "val",
    ]
)


def encode_residue(residue: str, node_feats: str):
    """Encode a residue"""
    residue = residue.lower()
    if node_feats == "label":
        return node_encoding[residue] + 1
    elif node_feats == "onehot":
        return onehot_encode(node_encoding[residue], len(node_encoding))
    else:
        raise ValueError("Unknown node_feats type!")


class Residue:
    """Residue class"""

    def __init__(self, line: str) -> None:
        self.name = line[17:20].strip()
        self.num = int(line[22:26].strip())
        self.chainID = line[21].strip()
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())


class Structure:
    """Structure class"""

    def __init__(self, filename: str, node_feats: str) -> None:
        self.residues = {}
        self.parse_file(filename)
        self.node_feats = node_feats

    def parse_file(self, filename: str) -> None:
        """Parse PDB file"""
        for line in open(filename, "r"):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res = Residue(line)
                self.residues[res.num] = res

    def get_coords(self) -> torch.Tensor:
        """Get coordinates of all atoms"""
        coords = [[res.x, res.y, res.z] for res in self.residues.values()]
        return torch.tensor(coords)

    def get_nodes(self) -> torch.Tensor:
        """Get features of all nodes of a graph"""
        return torch.tensor([encode_residue(res.name, self.node_feats) for res in self.residues.values()])

    def get_edges(self, threshold: float) -> torch.Tensor:
        """Get edges of a graph using threshold as a cutoff"""
        coords = self.get_coords()
        dist = torch.cdist(coords, coords)
        edges = torch.where(dist < threshold)
        edges = torch.cat([arr.view(-1, 1) for arr in edges], axis=1)
        edges = edges[edges[:, 0] != edges[:, 1]]
        return edges.t()

    def get_graph(self, threshold: float) -> dict:
        """Get a graph using threshold as a cutoff. 
        Resulting dict can be converted to torch_geometric.data.Data with Data(**d) where d is the output of this function."""
        nodes = self.get_nodes()
        edges = self.get_edges(threshold)
        return dict(x=nodes, edge_index=edges)


if __name__ == "__main__":
    import os
    import os.path as osp
    import pandas as pd
    from jsonargparse import CLI
    from joblib import Parallel, delayed
    from tqdm import tqdm
    import pickle

    def run(filename:str, output:str, threshold:float=5, node_feats:str="label") -> None:
        """Calculate the graph form for a single pdb file.
        Result is a pickle file containing a dictionary, keys are "x" and "edge_index" (consistent with torch_geometric.data.Data).

        Args:
            filename (str): path to original PDB file
            output (str): path to the output pickle file
            threshold (float, optional): Distance threshold in Angstroms. Defaults to 5.
            node_feats (str, optional): How to encode node features. Can be "onehot" or "label". Defaults to "label".
        """        
        struct = Structure(filename, node_feats).get_graph(threshold)
        with open(output, "wb") as file:
            file.write(pickle.dump(struct))


    def run_batch(pdb_dir: str, output: str, threads: int = 1, threshold: float = 5, node_feats: str = "label"):
        """Run the pipeline for a directory full of pdb files.
        Result is a pickle file containing a DataFrame of dictionaries, in each 
        keys are "x" and "edge_index" (consistent with torch_geometric.data.Data).

        Args:
            pdb_dir (str): path to directory of pdb files
            output (str): path to the output pickle file
            threads (int): number of threads to use
            threshold (float, optional): Distance threshold in Angstroms. Defaults to 5.
            node_feats (str, optional): How to encode node features. Can be "onehot" or "label". Defaults to "label".
        """     

        def get_graph(filename: str) -> dict:
            """Calculate a single graph from a file"""
            return Structure(filename, node_feats).get_graph(threshold)

        pdbs = [osp.join(pdb_dir, x) for x in os.listdir(pdb_dir)]
        data = Parallel(n_jobs=threads)(delayed(get_graph)(i) for i in tqdm(pdbs))
        df = pd.DataFrame(pd.Series(data, name="data"))
        df["filename"] = pdbs
        df["ID"] = df["filename"].apply(lambda x: x.split("/")[-1].split(".")[0])
        df.set_index("ID", inplace=True)
        df.drop("filename", axis=1).to_pickle(output)

    cli = CLI([run, run_batch])