import os
import shutil
import argparse
import subprocess
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, MMCIFParser, Polypeptide
from Bio import Align
from Bio.Align import substitution_matrices

def get_sequence_from_pdb(pdb_file_path):
    parser = PDBParser()
    try:
        structure = parser.get_structure("protein", pdb_file_path)
    except Exception:
        parser = MMCIFParser()
        structure = parser.get_structure("protein", pdb_file_path)

    sequences = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            chain_sequence = []
            for residue in chain:
                resname = residue.get_resname()
                if resname in Polypeptide.protein_letters_3to1:
                    chain_sequence.append(Polypeptide.protein_letters_3to1[resname])
            if chain_sequence:
                sequences[chain_id] = "".join(chain_sequence)
    return sequences

def group_pdbs(generated_pdbs_dir, output_base_dir):
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    print(f"Organizing generated PDBs from {generated_pdbs_dir} into {output_base_dir}...")
    for pdb_file_name in os.listdir(generated_pdbs_dir):
        if pdb_file_name.endswith(".pdb"):
            pdb_path = os.path.join(generated_pdbs_dir, pdb_file_name)
            sequences = get_sequence_from_pdb(pdb_path)
            
            if sequences:
                # Take the length of the longest chain for grouping
                seq_length = max([len(s) for s in sequences.values()])
                length_dir = os.path.join(output_base_dir, f"length_{seq_length}")
                if not os.path.exists(length_dir):
                    os.makedirs(length_dir)
                shutil.copy(pdb_path, os.path.join(length_dir, pdb_file_name))
            else:
                print(f"Warning: No sequences found in {pdb_file_name}. Skipping.")
    print("Organization complete.\n")

def find_unaligned_ggg_indices(original_sequences, generated_sequences):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    results = [] # List of tuples (chain_id, [list of 1-based indices])

    for gen_chain_id, gen_seq in generated_sequences.items():
        gen_seq_mask = [False] * len(gen_seq)
        for orig_chain_id, orig_seq in original_sequences.items():
            alignments = aligner.align(gen_seq, orig_seq)
            if alignments:
                best_alignment = max(alignments, key=lambda x: x.score)
                aligned_gen_seq = best_alignment[0]
                aligned_orig_seq = best_alignment[1]
                
                gen_idx = 0
                for i in range(len(aligned_gen_seq)):
                    gen_char = aligned_gen_seq[i]
                    orig_char = aligned_orig_seq[i]
                    if gen_char != '-':
                        if orig_char != '-':
                             gen_seq_mask[gen_idx] = True
                        gen_idx += 1

        unaligned_indices = []
        current_ggg_start = -1
        
        def process_range(start, end):
             length = end - start
             if length >= 3:
                 subseq = gen_seq[start:end]
                 if subseq.strip('G') == '':
                     for k in range(start, end):
                         unaligned_indices.append(k + 1)

        for i in range(len(gen_seq)):
            is_g = (gen_seq[i] == 'G')
            is_unaligned = (not gen_seq_mask[i])

            if is_g and is_unaligned:
                if current_ggg_start == -1:
                    current_ggg_start = i
            else:
                if current_ggg_start != -1:
                    process_range(current_ggg_start, i)
                    current_ggg_start = -1
        
        if current_ggg_start != -1:
            process_range(current_ggg_start, len(gen_seq))
            
        if unaligned_indices:
            results.append((gen_chain_id, unaligned_indices))

    return results

def main():
    parser = argparse.ArgumentParser(description="Auto run MPNN design on generated PDBs with unaligned GGG regions.")
    parser.add_argument("--original_pdb", type=str, required=True, help="Path to the original PDB file.")
    parser.add_argument("--generated_pdbs_folder", type=str, required=True, help="Path to the directory containing generated PDB files.")
    parser.add_argument("--output_dir", type=str, default="outputs", help="Directory where MPNN outputs will be stored.")
    parser.add_argument("--grouped_pdbs_dir", type=str, default="grouped_pdbs", help="Directory to store grouped PDBs.")

    args = parser.parse_args()

    # 1. Group PDBs
    group_pdbs(args.generated_pdbs_folder, args.grouped_pdbs_dir)

    # 2. Analyze and Run MPNN
    original_sequences = get_sequence_from_pdb(args.original_pdb)

    groups = [d for d in os.listdir(args.grouped_pdbs_dir) if os.path.isdir(os.path.join(args.grouped_pdbs_dir, d)) and d.startswith("length_")]
    groups.sort() 
    
    for group_name in groups:
        group_path = os.path.join(args.grouped_pdbs_dir, group_name)
        pdb_files = [f for f in os.listdir(group_path) if f.endswith(".pdb")]
        if not pdb_files:
            continue
            
        rep_pdb = pdb_files[0]
        rep_pdb_path = os.path.join(group_path, rep_pdb)
        
        print(f"Processing group: {group_name} (Representative: {rep_pdb})")
        
        gen_sequences = get_sequence_from_pdb(rep_pdb_path)
        
        targets = find_unaligned_ggg_indices(original_sequences, gen_sequences)
        
        if not targets:
            print(f"  No unaligned GGG regions found for {group_name}. Skipping.")
            continue
            
        chains_arg = []
        pos_arg = []
        
        for chain_id, indices in targets:
            chains_arg.append(chain_id)
            pos_str = " ".join(map(str, indices))
            pos_arg.append(pos_str)
            
        chains_str = " ".join(chains_arg)
        pos_str_full = ", ".join(pos_arg)
        
        abs_pdb_dir = os.path.abspath(group_path)
        abs_output_dir = os.path.abspath(os.path.join(args.output_dir, f"MPNN_design_{group_name}"))
        
        print(f"  Chains: {chains_str}")
        print(f"  Positions (truncated): {pos_str_full[:50]}...")
        
        subprocess_args = [
            "bash",
            "MPNN4_customize.sh",
            "--pdb_dir", abs_pdb_dir,
            "--output_dir", abs_output_dir,
            "--chains", chains_str,
            "--design_positions", pos_str_full
        ]
        
        cwd = "examples"
        
        try:
             subprocess.run(subprocess_args, cwd=cwd, check=True)
             print(f"  Completed {group_name}\n")
        except subprocess.CalledProcessError as e:
             print(f"  Error running MPNN for {group_name}: {e}\n")

if __name__ == "__main__":
    main()