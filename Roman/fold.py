import os
from os.path import join as pj, abspath as ap
import shutil
import sys
from argparse import ArgumentParser, Namespace
from typing import Tuple, List

from git import Repo


def error(msg: str) -> None:
    """
    Print an error message and abort the program.

    Args:
        msg: The error message to print.
    """
    print(msg)
    exit(1)


def read_fasta(r_filename: str) -> List[Tuple[str, str]]:
    """
    Read a fasta file into a list of pairs of sequence head and their actual amino acid sequence.

    Args:
        r_filename: Filename of the fasta file to be read in

    Return:
        List of tuples of sequence heads and amino acid sequences
    """
    with open(r_filename, "r") as r_fasta:
        r_seqs = []
        r_head = ""
        for r_line in r_fasta.readlines():
            if r_head != "":
                r_seqs.append((r_head, r_line.strip()))
                r_head = ""
            else:
                r_head = r_line.strip()
    return r_seqs


def parse_args() -> Namespace:
    """Parse the arguments"""
    parser = ArgumentParser(
        description="With this script you can easily install and run AlphaFold2, OpenFold (fast, PyTorch-based AF2), "
                    "and OmegaFold (db-free AF2) on input aminoacid sequences to compute their folding. So far, only "
                    "monomer-folding is supported by this script, multimer folds will follow (on request)."
    )
    parser.add_argument(
        "fasta",
        type=str,
        help="Path to fasta file or directory with sequences to query.\n"
             "If it's a file, every sequence will be treated a individual proteins.\n"
             "If it's a directory, you can also provide multiple chains for one protein to perform multimer folding",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Folder to store the output in. Every entry in the fasta file will get it's own folder.\n"
             "This folder must be empty if it already exists.",
    )
    parser.add_argument(
        "--openfold_script",
        type=str,
        default="",
        help="Path to the pythonscript to execute OpenFold."
    )
    parser.add_argument(
        "--omegafold_script",
        type=str,
        default="",
        help="Path to the pythonscript to execute OmegaFold and AlphaFold."
    )
    parser.add_argument(
        "--esmfold_script",
        type=str,
        default="",
        help="Path to the pythonscript to execute ESMFold."
    )
    parser.add_argument(
        "-a",
        "--alphafold_weights",
        type=str,
        default="",
        help="Path to the AlphaFold weights to use.",
    )
    parser.add_argument(
        "-o",
        "--openfold_weights",
        type=str,
        default="",
        help="Path to the OpenFold weights to use.",
    )
    parser.add_argument(
        "-m",
        "--omegafold_weights",
        type=str,
        default="",
        help="Path to the OmegaFold weights to use.",
    )
    parser.add_argument(
        "-e",
        "--esmfold",
        action='store_true',
        default=False,
        help="Flag indicating to run ESMFold."
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
    parser.add_argument(
        "-i",
        "--init",
        type=str,
        default="",
        help="Path to empty directory where to store the model scripts and their weights. This command will download "
             "the OmegaFold github repository, the OpenFold github repository, and the weights for all four models: "
             "AlphaFold2, OpenFold, OmegaFold, and ESMFold. This command will automatically set/overwrite all other "
             "optional parameters (except for -d and the weights if set) and execute a test run with the specified "
             "fasta file.\n! Please make sure to have enough memory on your disk for this step !\n"
             "This can also be used to easily set the optional parameters if the folder structure has been generated "
             "by a previous use of the -i parameter."
    )
    return parser.parse_args(sys.argv[1:])


# parse the arguments
args = parse_args()

# If an initialization folder is provided,
if args.init != "":
    af_weights_dir = pj(args.init, "weights", "alphafold")
    of_weights_dir = pj(args.init, "weights", "openfold")
    mf_weights_dir = pj(args.init, "weights", "omegafold")

    if len(os.listdir(args.init)) == 0:
        # create the folder if not exists
        os.makedirs(args.init, exist_ok=True)

        # Clone the repositories for OpenFold and OmegaFold
        # TODO: make this clean in own repo with copied/forked repos to keep everything working in case of later commits
        Repo.clone_from("https://github.com/aqlaboratory/openfold", pj(args.init, "openfold"))
        os.system(f"pip install {pj(args.init, 'openfold')}/.")
        os.system(
            f"wget "
            f"-P {pj(args.init, 'openfold', 'openfold', 'resources')} "
            f"\"https://git.scicore.unibas.ch/schwede/openstructure/-/raw/"
            f"7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt\""
        )
        Repo.clone_from("https://github.com/HeliXonProtein/OmegaFold.git", pj(args.init, "omegafold"))
        os.system(
            f"wget "
            f"-P {pj(args.init, 'esmfold')} "
            f"https://raw.githubusercontent.com/facebookresearch/esm/main/scripts/esmfold_inference.py"
        )

        # create folder to store the weights
        os.makedirs(af_weights_dir, exist_ok=True)
        os.makedirs(of_weights_dir, exist_ok=True)
        os.makedirs(mf_weights_dir, exist_ok=True)

        # OpenFold will be downloaded automatically

        # Download the AlphaFold weights
        os.system(f"wget -P {af_weights_dir} https://storage.googleapis.com/alphafold/alphafold_params_2022-01-19.tar")
        os.system(
            f"tar --extract --verbose --file={pj(af_weights_dir, 'alphafold_params_2022-01-19.tar')} "
            f"--directory={af_weights_dir} --preserve-permissions"
        )
        os.system(f"rm {pj(af_weights_dir, 'alphafold_params_2022-01-19.tar')}")

        # Download the OpenFold weights
        os.system(f"git clone https://huggingface.co/nz/OpenFold {of_weights_dir}")
        os.system(f"rm -rf {pj(of_weights_dir, '.git')}")

    # set the arguments for the following
    args.omegafold_script = pj(args.init, "omegafold", "main.py")
    if args.omegafold_weights == "":
        args.omegafold_weights = pj(mf_weights_dir, "omegafold_weights.pt")
    args.openfold_script = pj(args.init, "openfold", "run_pretrained_openfold.py")
    if args.openfold_weights == "":
        args.openfold_weights = pj(of_weights_dir, "finetuning_*.pt")
    if args.alphafold_weights == "":
        args.alphafold_weights = pj(af_weights_dir, "params_model_*.npz")
    args.esmfold_script = pj(args.init, "esmfold", "esmfold_inference.py")

# Check if the script to use the AlphaFold weights is provided
if args.alphafold_weights != "" and args.openfold_script == "":
    error("In order to execute AlphaFold, both the AlphaFold weights and the script of OpenFold must be provided.")

# Check if the script to use the OpenFold weights is provided
if args.openfold_weights != "" and args.openfold_script == "":
    error("In order to execute OpenFold, both the OpenFold weights and the script of OpenFold must be provided.")

# Check if the OpenFold script is provided, at least one of the weights are provided
if args.openfold_script != "" and args.alphafold_weights == "" and args.openfold_weights == "":
    error("In order to execute AlphaFold or OpenFold, either the AlphaFold weights or the OpenFold weights must be "
          "provided.")

# Check that either OmegaFold script and weights are provided or both are not provided
if (args.omegafold_weights != "") != (args.omegafold_script != ""):
    error("In order to execute OmegaFold, both the OmegaFold weights and the script of OmegaFold must be provided.")

# Determine what to execute
run_alphafold = args.alphafold_weights != "" and args.openfold_script != ""
run_openfold = args.openfold_weights != "" and args.openfold_script != ""
run_omegafold = args.omegafold_weights != "" and args.omegafold_script != ""
run_esmfold = args.esmfold_script != "" and args.esmfold_script != ""

# Check if anything is provided to execute
if not any([run_alphafold, run_openfold, run_omegafold, run_esmfold]):
    error("Either weights for AlphaFold, OpenFold, or OmegaFold have to be passed.")

# Check if the fasta path is valid
if not os.path.exists(args.fasta):
    error("Provided FASTA path is neither a file nor a directory.")

# Check if data is provided if necessary
if (run_alphafold or run_openfold) and args.data == "":
    error("In order to use AlphaFold or OpenFold, you have to provide the path to the datasets.")

os.makedirs(args.output, exist_ok=True)
if len(os.listdir(args.output)) != 0:
    print("WARNING: Output directory is not empty. Existing data might be overwritten. Continue? [y]/n")
    char = input()
    if char not in ["", "y"]:
        exit(0)

# set up directory structure for folding
fasta_input = pj(args.output, "fasta")
alphafold_output = pj(args.output, "alphafold")
openfold_output = pj(args.output, "openfold")
omegafold_output = pj(args.output, "omegafold")
esmfold_output = pj(args.output, "esmfold")
fold_output = pj(args.output, "folds")

if run_alphafold:
    os.makedirs(alphafold_output, exist_ok=True)
    os.makedirs(fasta_input, exist_ok=True)
if run_openfold:
    os.makedirs(openfold_output, exist_ok=True)
    os.makedirs(fasta_input, exist_ok=True)
if run_omegafold:
    os.makedirs(omegafold_output, exist_ok=True)
if run_esmfold:
    os.makedirs(esmfold_output, exist_ok=True)
os.makedirs(fold_output, exist_ok=True)

# copy the fasta file into single sequence fasta files for AlphaFold and OpenFold
if os.path.isfile(args.fasta) and (run_alphafold or run_openfold):
    seqs = read_fasta(args.fasta)

    for head, seq in seqs:
        name = head
        if name.startswith(">"):
            name = name[1:]
        if "|" in name:
            name = head.split("|")[1]

        prot_filename = pj(fasta_input, name + ".fasta")
        with open(prot_filename, "w") as fasta:
            print(head, file=fasta)
            print(seq, file=fasta)

# For multimer folding setup fasta files accordingly
elif os.path.isdir(args.fasta):
    error("Multimer folding not supported so far.")
    fasta_input = args.fasta
    with open(pj(args.output, "all_sequences.fasta"), "w") as fasta:
        for filename in os.listdir(fasta_input):
            seqs = read_fasta(pj(fasta_input, filename))
            print(f">{filename.split('.')[0]}", file=fasta)
            print(":".join(value for _, value in seqs), file=fasta)

if run_alphafold or run_openfold:
    specifics = []
    if run_openfold:
        # Setup OpenFold specific parameters
        if '*' in args.openfold_weights:
            # Run multiple weight initializations for OpenFold
            for i in range(2, 6):
                specifics.append(
                    f"--output_dir {ap(openfold_output)} "
                    f"--openfold_checkpoint_path {ap(args.openfold_weights.replace('*', str(i)))} "
                    f"--output_postfix openfold_model_{i} "
                )
        else:
            # Run only a single initialization for AlphaFold
            specifics.append(
                f"--output_dir {ap(openfold_output)} "
                f"--openfold_checkpoint_path {ap(args.openfold_weights)} "
                f"--output_postfix openfold "
            )
    if run_alphafold:
        # Setup OpenFold specific parameters
        if '*' in args.openfold_weights:
            # Run multiple weight initializations for AlphaFold
            for i in range(1, 6):
                specifics.append(
                    f"--output_dir {ap(alphafold_output)} "
                    f"--jax_param_path {ap(args.alphafold_weights.replace('*', str(i)))} "
                    f"--output_postfix alphafold_model_{i} "
                )
        else:
            # Run only a single initialization for AlphaFold
            specifics.append(
                f"--output_dir {ap(alphafold_output)} "
                f"--jax_param_path {ap(args.alphafold_weights)} "
                f"--output_postfix alphafold "
            )

    for output_and_model_parameters in specifics:
        # Run OpenFold script with folding specific parameters
        cmd = \
            f"python " \
            f"{args.openfold_script} " \
            f"{ap(fasta_input)} " \
            f"{pj(args.data, 'pdb_mmcif', 'mmcif_files')} " \
            f"{output_and_model_parameters} " \
            f"--uniref90_database_path {pj(args.data, 'uniref90', 'uniref90.fasta')} " \
            f"--mgnify_database_path {pj(args.data, 'mgnify', 'mgy_clusters_2018_12.fa')} " \
            f"--pdb70_database_path {pj(args.data, 'pdb70', 'pdb70')} " \
            f"--uniclust30_database_path {pj(args.data, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')} " \
            f"--bfd_database_path {pj(args.data, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')} " \
            f"--jackhmmer_binary_path {shutil.which('jackhmmer')} " \
            f"--hhblits_binary_path {shutil.which('hhblits')} " \
            f"--hhsearch_binary_path {shutil.which('hhsearch')} " \
            f"--kalign_binary_path {shutil.which('kalign')} " \
            f"--model_device cuda:0 "
        if os.path.exists(pj(ap(openfold_output), "alignments")):
            cmd += f"--use_precomputed_alignments {pj(ap(openfold_output), 'alignments')} "
        # print(cmd)
        os.system(cmd)

        # copy results into one central folder
        os.system(f"cp {pj(output_and_model_parameters.split(' ')[1], 'predictions', '*')} {fold_output}")

if run_omegafold:
    # Run OmegaFold
    cmd = f"python " \
          f"{args.omegafold_script} " \
          f"{args.fasta if os.path.isfile(args.fasta) else pj(ap(args.output), 'all_sequences.fasta')} " \
          f"{omegafold_output} " \
          f"--weights_file {args.omegafold_weights} " \
          f"--device cuda:0 " \
          f"--num_cycle 2 "
    # print(cmd)
    os.system(cmd)

    # add postfix to result files
    for file in os.listdir(omegafold_output):
        os.system(f"mv {pj(omegafold_output, file)} {pj(omegafold_output, file[:-4] + '_omegafold.pdb')}")

    # copy results into one central folder
    os.system(f"cp {pj(omegafold_output, '*')} {fold_output}")

if run_esmfold:
    cmd = f"python " \
          f"{args.esmfold_script} " \
          f"-i {args.fasta if os.path.isfile(args.fasta) else pj(ap(args.output), 'all_sequences.fasta')} " \
          f"-o {esmfold_output} "
    print(cmd)
    os.system(cmd)

    # add postfix to result files
    for file in os.listdir(omegafold_output):
        os.system(f"mv {pj(omegafold_output, file)} {pj(omegafold_output, file[:-4] + '_esmfold.pdb')}")

    # copy results into one central folder
    os.system(f"cp {pj(omegafold_output, '*')} {fold_output}")
