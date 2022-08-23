import os.path
import shutil
import sys
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "fasta",
    type=str,
    help="Path to fasta file with sequences to query. "
         "This file must contain all sequences you want to get structures for.",
)
parser.add_argument(
    "output",
    type=str,
    help="Folder to store the output in. Every entry in the fasta file will get it's own folder.\n"
         "This folder must be empty if it already exists.",
)
parser.add_argument(
    "--openfold-script",
    type=str,
    dest="openfold_script",
    default="",
    help="Path to the pythonscript to execute OpenFold."
)
parser.add_argument(
    "--omegafold-script",
    type=str,
    dest="omegafold_script",
    default="",
    help="Path to the pythonscript to execute OmegaFold."
)
parser.add_argument(
    "-a",
    "--alphafold",
    type=str,
    default="",
    help="Path to the AlphaFold weights to use.",
)
parser.add_argument(
    "-o",
    "--openfold",
    type=str,
    default="",
    help="Path to the OpenFold weights to use.",
)
parser.add_argument(
    "-m",
    "--omegafold",
    type=str,
    default="",
    help="Path to the OmegaFold weights to use.",
)
parser.add_argument(
    "-d",
    "--data",
    type=str,
    default="",
    help="Path to directory containing all databases necessary to run AlphaFold or OpenFold.\n"
         "! This argument is required when using AlphaFold or OpenFold !\n"
         "! This argument mus be passed as absolute path !",
)

args = parser.parse_args(sys.argv[1:])

if args.alphafold != "" and args.openfold_script == "":
    print("In order to execute AlphaFold, both the AlphaFold weights and the script of OpenFold must be provided.")
    exit(1)

if args.openfold != "" and args.openfold_script == "":
    print("In order to execute OpenFold, both the OpenFold weights and the script of OpenFold must be provided.")
    exit(1)

if args.openfold_script != "" and args.alphafold == "" and args.openfold == "":
    print("In order to execute AlphaFold or OpenFold, either the AlphaFold weights or the OpenFold weights must be "
          "provided.")
    exit(1)

if (args.omegafold != "") != (args.omegafold_script != ""):
    print("In order to execute OmegaFold, both the OmegaFold weights and the script of OmegaFold must be provided.")
    exit(1)

run_alphafold = args.alphafold != "" and args.openfold_script != ""
run_openfold = args.openfold != "" and args.openfold_script != ""
run_omegafold = args.omegafold != "" and args.omegafold_script != ""

if not any([run_alphafold, run_openfold, run_omegafold]):
    print("Either weights for AlphaFold, OpenFold, or OmegaFold have to be passed.")
    exit(1)

if not os.path.exists(args.fasta):
    print("Provided FASTA file is invalid.")
    exit(1)

if (run_alphafold or run_openfold) and args.data == "":
    print("In order to use AlphaFold or OpenFold, you have to provide the path to the datasets.")
    exit(1)

os.makedirs(args.output, exist_ok=True)
alphafold_output = os.path.join(args.output, "alphafold")
openfold_output = os.path.join(args.output, "openfold")
fasta_input = os.path.join(args.output, "fasta")
omegafold_output = os.path.join(args.output, "omegafold")

if run_alphafold:
    os.makedirs(alphafold_output, exist_ok=True)
    os.makedirs(fasta_input, exist_ok=True)
if run_openfold:
    os.makedirs(openfold_output, exist_ok=True)
    os.makedirs(fasta_input, exist_ok=True)
if run_omegafold:
    os.makedirs(omegafold_output, exist_ok=True)

if len(os.listdir(args.output)) != 0:
    print("Output path is not empty.")
    exit(1)

with open(args.fasta, "r") as fasta:
    seqs = []
    head = ""
    for line in fasta.readlines():
        if head != "":
            seqs.append((head, line.strip()))
            head = ""
        else:
            head = line.strip()

# TODO: make this clean in own repo with copied/forked repos
# Repo.clone_from("https://github.com/aqlaboratory/openfold", os.path.join(args.output, "openfold"))
# Repo.clone_from("https://github.com/HeliXonProtein/OmegaFold.git", os.path.join(args.output, "omegafold"))

if run_alphafold or run_openfold:
    for head, seq in seqs:

        name = head
        if name.startswith(">"):
            name = name[1:]
        if "|" in name:
            name = head.split("|")[1]

        prot_filename = os.path.join(fasta_input, name + ".fasta")
        with open(prot_filename, "w") as fasta:
            print(head, file=fasta)
            print(seq, file=fasta)

    specifics = []
    if run_alphafold:
        specifics.append(f"--output_dir {os.path.abspath(alphafold_output)} "
                         f"--jax_param_path {os.path.abspath(args.alphafold)} ")
    if run_openfold:
        specifics.append(f"--output_dir {os.path.abspath(openfold_output)} "
                         f"--openfold_checkpoint_path {os.path.abspath(args.openfold)} ")

    for output_and_model_parameters in specifics:
        os.system(
            f"python "
            f"{args.openfold_script} "
            f"{os.path.abspath(fasta_input)} "
            f"{os.path.join(args.data, 'pdb_mmcif', 'mmcif_files')} "
            f"{output_and_model_parameters} "
            f"--uniref90_database_path {os.path.join(args.data, 'uniref90', 'uniref90.fasta')} "
            f"--mgnify_database_path {os.path.join(args.data, 'mgnify', 'mgy_clusters_2018_12.fa')} "
            f"--pdb70_database_path {os.path.join(args.data, 'pdb70', 'pdb70')} "
            f"--uniclust30_database_path {os.path.join(args.data, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')} "
            f"--bfd_database_path {os.path.join(args.data, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')} "
            f"--jackhmmer_binary_path {shutil.which('jackhmmer')} "
            f"--hhblits_binary_path {shutil.which('hhblits')} "
            f"--hhsearch_binary_path {shutil.which('hhsearch')} "
            f"--kalign_binary_path {shutil.which('kalign')} "
            f"--model_device cuda:0 "
        )

if run_omegafold:
    os.system(f"python {args.omegafold_script} {args.fasta} {omegafold_output} --weights_file {args.omegafold}")
