import os
import glob
import re
import argparse

def split_fasta(root_dir):
    """
    Searches for .fa files recursively in root_dir, 
    splits multi-sequence fasta files into individual single-sequence fasta files,
    and saves them directly in the root_dir.
    Ignores sequences without 'sample=' in the header.
    """
    # Find all .fa files recursively
    # We assume the user wants to process files ending in .fa
    fa_files = glob.glob(os.path.join(root_dir, '**', '*.fa'), recursive=True)
    
    # Filter out files that look like the output files (ending in _sample_X.fa) to avoid re-processing
    # This prevents infinite loops if the script is run multiple times
    input_files = []
    for f in fa_files:
        if re.search(r'_sample_\d+\.fa$', f):
            continue
        input_files.append(f)
        
    if not input_files:
        print(f"No suitable input .fa files found in {root_dir}")
        return

    print(f"Found {len(input_files)} input .fa files. Processing...")

    for fa_file in input_files:
        dir_name = os.path.dirname(fa_file)
        base_name = os.path.basename(fa_file)
        file_stem = os.path.splitext(base_name)[0]
        
        # Output directly to root_dir
        output_dir = root_dir
        
        try:
            with open(fa_file, 'r') as f:
                content = f.read()
        except Exception as e:
            print(f"Error reading {fa_file}: {e}")
            continue
        
        # Split by '>'
        entries = content.split('>')
        
        count = 0
        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue
            
            lines = entry.split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:]) # Join remaining lines as sequence
            
            # Determine a filename based on the header
            # ProteinMPNN headers typically look like: 
            # >design4_1, score=... (Original)
            # >T=0.1, sample=1, score=... (Generated)
            
            sample_match = re.search(r'sample=(\d+)', header)
            
            if sample_match:
                sample_num = sample_match.group(1)
                new_filename = f"{file_stem}_sample_{sample_num}.fa"
            else:
                # Skip original sequences or anything without a sample number
                continue

            output_path = os.path.join(output_dir, new_filename)
            
            try:
                with open(output_path, 'w') as out_f:
                    out_f.write(f">{header}\n{sequence}\n")
                count += 1
            except Exception as e:
                print(f"Error writing to {output_path}: {e}")

        print(f"Processed {base_name}: {count} sequences extracted to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split ProteinMPNN multi-sequence FASTA files into individual sample files.")
    parser.add_argument("directory", help="The directory to search for .fa files recursively.")
    args = parser.parse_args()

    if os.path.exists(args.directory):
        split_fasta(args.directory)
    else:
        print(f"Error: Directory '{args.directory}' not found.")