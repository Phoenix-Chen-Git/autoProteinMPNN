import os
import shutil
import argparse
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, MMCIFParser, Polypeptide
from Bio import Align
from Bio.Align import substitution_matrices

def get_sequence_from_pdb(pdb_file_path):
    """
    Parses a PDB file and returns a dictionary of sequences with chain IDs as keys,
    and a dictionary of residue lists with chain IDs as keys.
    """
    parser = PDBParser()
    try:
        structure = parser.get_structure("protein", pdb_file_path)
    except Exception:
        parser = MMCIFParser() # Try CIF if PDB fails
        structure = parser.get_structure("protein", pdb_file_path)

    sequences = {}
    residue_lists = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            chain_sequence = []
            chain_residues = []
            for residue in chain:
                resname = residue.get_resname()
                if resname in Polypeptide.protein_letters_3to1: # Check for standard amino acids
                    chain_sequence.append(Polypeptide.protein_letters_3to1[resname])
                    chain_residues.append(residue)
            if chain_sequence: # Only add if a sequence was actually built
                sequences[chain_id] = "".join(chain_sequence)
                residue_lists[chain_id] = chain_residues
    return sequences, residue_lists

def find_unaligned_ggg(original_sequences, generated_sequences, generated_residue_lists):
    """
    Performs alignment and finds continuous 'G' stretches in generated sequences
    that do not align with ANY original chain.
    """
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    results = []

    for gen_chain_id, gen_seq in generated_sequences.items():
        # Initialize a mask for the generated sequence: False means unaligned
        gen_seq_mask = [False] * len(gen_seq)

        # Iteratively align against ALL original chains
        for orig_chain_id, orig_seq in original_sequences.items():
            alignments = aligner.align(gen_seq, orig_seq)
            if alignments:
                # Take the best alignment
                best_alignment = max(alignments, key=lambda x: x.score)
                aligned_gen_seq = best_alignment[0]
                aligned_orig_seq = best_alignment[1]
                
                # Update mask: Mark positions as True if they align to a residue (not a gap) in the original
                gen_idx = 0
                for i in range(len(aligned_gen_seq)):
                    gen_char = aligned_gen_seq[i]
                    orig_char = aligned_orig_seq[i]
                    
                    if gen_char != '-':
                        # If this position in gen_seq aligns to a residue in orig_seq (not a gap)
                        if orig_char != '-':
                             gen_seq_mask[gen_idx] = True
                        gen_idx += 1

        # Now find continuous 'G' stretches that are still unmasked (False)
        current_ggg_start = -1
        for i in range(len(gen_seq)):
            is_g = (gen_seq[i] == 'G')
            is_unaligned = (not gen_seq_mask[i])

            if is_g and is_unaligned:
                if current_ggg_start == -1:
                    current_ggg_start = i
            else:
                if current_ggg_start != -1:
                    # End of a potential GGG stretch
                    # Check length filter (e.g., >= 3) to avoid single G noise
                    length = i - current_ggg_start
                    if length >= 3:
                        # Verify strictly 'G' (redundant here but good for safety)
                         if gen_seq[current_ggg_start:i].strip('G') == '':
                            start_res = generated_residue_lists[gen_chain_id][current_ggg_start]
                            end_res = generated_residue_lists[gen_chain_id][i - 1]
                            results.append({
                                "chain_id": gen_chain_id,
                                "start_res": start_res.get_id()[1],
                                "end_res": end_res.get_id()[1],
                                "sequence": gen_seq[current_ggg_start:i]
                            })
                    current_ggg_start = -1
        
        # Check if sequence ended with a GGG stretch
        if current_ggg_start != -1:
            length = len(gen_seq) - current_ggg_start
            if length >= 3:
                 if gen_seq[current_ggg_start:].strip('G') == '':
                    start_res = generated_residue_lists[gen_chain_id][current_ggg_start]
                    end_res = generated_residue_lists[gen_chain_id][len(gen_seq) - 1]
                    results.append({
                        "chain_id": gen_chain_id,
                        "start_res": start_res.get_id()[1],
                        "end_res": end_res.get_id()[1],
                        "sequence": gen_seq[current_ggg_start:]
                    })

    return results

def main(original_pdb_path, generated_pdbs_dir, output_base_dir="grouped_pdbs"):
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    # 1. Organize Generated PDBs by Length
    print(f"Organizing generated PDBs from {generated_pdbs_dir}...")
    for pdb_file_name in os.listdir(generated_pdbs_dir):
        if pdb_file_name.endswith(".pdb"):
            pdb_path = os.path.join(generated_pdbs_dir, pdb_file_name)
            sequences, _ = get_sequence_from_pdb(pdb_path)
            
            # Assuming the first chain's length for grouping, or taking the longest
            if sequences:
                # Take the length of the longest chain for grouping
                seq_length = max([len(s) for s in sequences.values()])
                length_dir = os.path.join(output_base_dir, f"length_{seq_length}")
                if not os.path.exists(length_dir):
                    os.makedirs(length_dir)
                shutil.copy(pdb_path, os.path.join(length_dir, pdb_file_name))
            else:
                print(f"Warning: No sequences found in {pdb_file_name}. Skipping.")
    print("Organization complete.")

    # Get original PDB sequences
    original_sequences, _ = get_sequence_from_pdb(original_pdb_path)
    if not original_sequences:
        print(f"Error: No sequences found in original PDB: {original_pdb_path}")
        return

    # 2. Analyze Each Length Group
    print("Analyzing each length group...")
    for length_dir_name in os.listdir(output_base_dir):
        length_dir_path = os.path.join(output_base_dir, length_dir_name)
        if os.path.isdir(length_dir_path) and length_dir_name.startswith("length_"):
            pdb_files_in_group = [f for f in os.listdir(length_dir_path) if f.endswith(".pdb")]
            if pdb_files_in_group:
                # 2.a. Select one representative PDB
                representative_pdb_name = pdb_files_in_group[0]
                representative_pdb_path = os.path.join(length_dir_path, representative_pdb_name)
                
                print(f"  Analyzing group {length_dir_name} using {representative_pdb_name}")
                
                # 2.b. Extract sequences and residue lists for the representative PDB
                generated_sequences, generated_residue_lists = get_sequence_from_pdb(representative_pdb_path)
                if not generated_sequences:
                    print(f"    Warning: No sequences found in representative PDB {representative_pdb_name}. Skipping group.")
                    continue

                # 3. & 4. & 5. Align, find GGG, and output
                unaligned_ggg_regions = find_unaligned_ggg(original_sequences, generated_sequences, generated_residue_lists)
                
                if unaligned_ggg_regions:
                    print(f"    Unaligned 'GGG...GGG' regions found in {representative_pdb_name}:")
                    for region in unaligned_ggg_regions:
                        position_list = ' '.join(map(str, range(region['start_res'], region['end_res'] + 1)))
                        print(f"      '{region['chain_id']}' '{position_list}'")
                else:
                    print(f"    No unaligned 'GGG...GGG' regions found in {representative_pdb_name}.")
            else:
                print(f"  No PDB files found in group {length_dir_name}.")
    print("Analysis complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze generated PDBs for unaligned GGG stretches.")
    parser.add_argument("--original_pdb", type=str, required=True,
                        help="Path to the original PDB file.")
    parser.add_argument("--generated_pdbs_folder", type=str, required=True,
                        help="Path to the directory containing generated PDB files.")
    parser.add_argument("--output_base_dir", type=str, default="grouped_pdbs",
                        help="Base directory to store organized PDBs and analysis results. Defaults to 'grouped_pdbs'.")

    args = parser.parse_args()
    
    main(args.original_pdb, args.generated_pdbs_folder, args.output_base_dir)
