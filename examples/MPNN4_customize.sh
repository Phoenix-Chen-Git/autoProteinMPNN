#!/bin/bash
#SBATCH -J ProteinMPNN_non_fixed
#SBATCH -p gpu_l40
#SBATCH -N 1
#SBATCH -o ProteinMPNN_non_fixed_%j.out
#SBATCH -e ProteinMPNN_non_fixed_%j.err
#SBATCH --no-requeue
#SBATCH -A yulongli_g1
#SBATCH --qos=yulonglil40
#SBATCH --gres=gpu:1
#SBATCH --overcommit
#SBATCH --mincpus=8

usage() {
    cat <<'EOF'
Usage: sbatch MPNN4_customize.sh [--pdb_dir DIR] [--output_dir DIR] [--chains "A B"] [--design_positions "..."]

Options:
  --pdb_dir DIR      Directory containing input PDB complexes (default ../inputs/PDB_complexes/pdbs/)
  --output_dir DIR   Directory where outputs will be written (default ../outputs/example_4_non_fixed_outputs)
  --chains STR       Space-separated chains to design (default "A C")
  --design_positions STR  Residue ranges for --specify_non_fixed (default "1 2 3 4 5 6 7 8 9 10, 3 4 5 6 7 8")
  -h, --help         Show this help message and exit
EOF
    exit 0
}

folder_with_pdbs="../inputs/PDB_complexes/pdbs/"
output_dir="../outputs/example_4_non_fixed_outputs"
chains_to_design="A C"
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
design_only_positions="1 2 3 4 5 6 7 8 9 10, 3 4 5 6 7 8" #design only these residues; use flag --specify_non_fixed

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pdb_dir)
            folder_with_pdbs="$2"
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --chains)
            chains_to_design="$2"
            shift 2
            ;;
        --design_positions)
            design_only_positions="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            ;;
    esac
done

if [ ! -d "$output_dir" ]
then
    mkdir -p "$output_dir"
fi


path_for_parsed_chains="$output_dir/parsed_pdbs.jsonl"
path_for_assigned_chains="$output_dir/assigned_pdbs.jsonl"
path_for_fixed_positions="$output_dir/fixed_pdbs.jsonl"

python ../helper_scripts/parse_multiple_chains.py --input_path="$folder_with_pdbs" --output_path="$path_for_parsed_chains"

python ../helper_scripts/assign_fixed_chains.py --input_path="$path_for_parsed_chains" --output_path="$path_for_assigned_chains" --chain_list "$chains_to_design"

python ../helper_scripts/make_fixed_positions_dict.py --input_path="$path_for_parsed_chains" --output_path="$path_for_fixed_positions" --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python ../protein_mpnn_run.py \
    --jsonl_path "$path_for_parsed_chains" \
    --chain_id_jsonl "$path_for_assigned_chains" \
    --fixed_positions_jsonl "$path_for_fixed_positions" \
    --out_folder "$output_dir" \
        --num_seq_per_target 10 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 10
