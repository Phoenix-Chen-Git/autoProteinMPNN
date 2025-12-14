# Auto Run ProteinMPNN (autoProteinMPNN)

This repo includes a small “auto-run” workflow on top of ProteinMPNN that:

1. Groups generated PDBs by length.
2. Detects which residues should be redesigned (two modes).
3. Runs ProteinMPNN only on those residues.
4. Splits the designed FASTA output into per-target FASTA files.

## Files

### `auto_run_mpnn.py`
Python entrypoint that performs analysis + launches ProteinMPNN.

**What it does**
- Reads sequences from PDB/mmCIF (per chain).
- Groups all PDBs under `--generated_pdbs_folder` into `--grouped_pdbs_dir/length_<N>/` based on the longest chain length.
- For each length group, picks one representative PDB to determine design positions.
- Runs `examples/MPNN4_customize.sh` once per group with:
  - the PDB directory for that group
  - which chain(s) to design
  - which residue positions to design
  - temperature and number of sequences per PDB

**Two modes for choosing designable positions**
- **Auto-alignment mode** (provide `--original_pdb`):
  - Aligns each generated chain sequence to the original chain sequences.
  - Finds *unaligned* stretches of glycine ("GGG…") of length ≥ 3.
  - Only those unaligned GGG residues become designable.
- **Fixed-sequences mode** (provide `--fixed_sequences`):
  - Treats the given amino-acid sequences as “must keep fixed”.
  - Marks all occurrences of each fixed sequence in each chain as fixed.
  - Everything not part of those fixed sequences becomes designable.

**Key CLI options**
- `--generated_pdbs_folder DIR` (required)
- `--output_dir DIR` (default: `outputs`)
- `--grouped_pdbs_dir DIR` (default: `grouped_pdbs`)
- `--original_pdb PATH` (required unless using `--fixed_sequences`)
- `--fixed_sequences SEQ [SEQ ...]` (alternative to `--original_pdb`)
- `--sampling_temp STR` (default: `0.1`)
- `--num_seq_per_target INT` (default: `10`)  ← *this is “sequences per PDB”*

**Examples**
Auto-alignment mode:
```bash
python auto_run_mpnn.py \
  --original_pdb design4/original/design4.pdb \
  --generated_pdbs_folder design4/pdbs_test \
  --output_dir test_output \
  --grouped_pdbs_dir test_output_grouped \
  --sampling_temp 0.2 \
  --num_seq_per_target 25
```

Fixed-sequences mode:
```bash
python auto_run_mpnn.py \
  --generated_pdbs_folder design4/pdbs_test \
  --output_dir test_output \
  --grouped_pdbs_dir test_output_grouped \
  --fixed_sequences "SEQUENCE_1" "SEQUENCE_2" \
  --sampling_temp 0.15 \
  --num_seq_per_target 10
```

---

### `run_autoProteinMPNN.srp`
A SLURM wrapper script that:
- Activates the `mlfold` conda environment.
- Calls `auto_run_mpnn.py` with the right arguments.
- Runs `split_mpnn_seqs.py` to split MPNN FASTA outputs.

**Options**
- `-t <temp>`: sampling temperature passed through to ProteinMPNN
- `-n <num>`: number of sequences generated per PDB/target

**Usage**
```bash
sbatch run_autoProteinMPNN.srp [-t SAMPLING_TEMP] [-n NUM_SEQ_PER_TARGET] \
  <generated_pdbs_dir> <output_dir> [original_pdb_or_fixed_seq1] [fixed_seq2] ...
```

Examples:
```bash
sbatch run_autoProteinMPNN.srp -t 0.2 -n 25 design4/pdbs_test test_output design4/original/design4.pdb
```

```bash
sbatch run_autoProteinMPNN.srp -t 0.15 -n 10 design4/pdbs_test test_output "SEQUENCE_1" "SEQUENCE_2"
```

## Notes
- Temperature (`--sampling_temp`) can be a single value like `0.1` (and ProteinMPNN also supports multiple temperatures as a space-separated string, if you want diversity).
- The wrapper ensures `--num_seq_per_target` and `--sampling_temp` are forwarded end-to-end (they used to be hard-coded in `examples/MPNN4_customize.sh`).
