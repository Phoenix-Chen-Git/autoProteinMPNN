import os
import argparse
import sys

def trim_sequence(seq, start_seq, end_seq):
    """
    Trims the sequence to keep [start_seq...middle...end_seq].
    Returns the trimmed sequence or None if pattern not found.
    """
    start_idx = seq.find(start_seq)
    if start_idx == -1:
        return None
    
    # Find end_seq after start_seq
    # We search starting from where start_seq begins to ensure we find an end after the start.
    # If the user implies [start]...[end], end must be after start.
    end_idx = seq.find(end_seq, start_idx + len(start_seq))
    
    if end_idx == -1:
        return None
    
    # Return the slice including start_seq and end_seq
    return seq[start_idx : end_idx + len(end_seq)]

def process_file(file_path, start_seq, end_seq):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    new_content = []
    header = None
    seq_lines = []
    file_modified = False
    
    def process_entry(h, s):
        trimmed = trim_sequence(s, start_seq, end_seq)
        if trimmed:
            return h + "\n" + trimmed + "\n"
        else:
            # If pattern not found, keep original but warn
            print(f"Warning: Pattern [start]...[end] not found in {file_path} for {h.strip()}. Keeping original.")
            return h + "\n" + s + "\n"

    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith('>') :
            if header:
                full_seq = "".join(seq_lines)
                new_content.append(process_entry(header, full_seq))
            header = line
            seq_lines = []
        else:
            seq_lines.append(line)
            
    if header:
        full_seq = "".join(seq_lines)
        new_content.append(process_entry(header, full_seq))
    
    # Write back to file
    try:
        with open(file_path, 'w') as f:
            f.writelines(new_content)
        # print(f"Processed {file_path}") 
    except Exception as e:
        print(f"Error writing to {file_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Trim FASTA sequences to contain only [start_seq]...[end_seq].")
    parser.add_argument("start_seq", help="The starting DNA/Protein sequence (inclusive).")
    parser.add_argument("end_seq", help="The ending DNA/Protein sequence (inclusive).")
    parser.add_argument("fasta_dir", help="Directory to search for .fa/.fasta files.")
    
    args = parser.parse_args()
    
    start_seq = args.start_seq
    end_seq = args.end_seq
    target_dir = args.fasta_dir

    if not os.path.isdir(target_dir):
        print(f"Error: Directory \'{target_dir}\' does not exist.")
        sys.exit(1)

    print(f"Trimming sequences in \'{target_dir}\'...")
    print(f"Start sequence: {start_seq}")
    print(f"End sequence:   {end_seq}")
    
    count = 0
    for file_name in os.listdir(target_dir):
        if file_name.endswith('.fa') or file_name.endswith('.fasta'):
            file_path = os.path.join(target_dir, file_name)
            if os.path.isfile(file_path):
                process_file(file_path, start_seq, end_seq)
                count += 1
                
    print(f"Finished processing {count} files.")

if __name__ == "__main__":
    main()
