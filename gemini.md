use conda env: mlfold.
This device is using M5 silicon of apple.
test data: is under test/ folder
test input: test/origin_pdb/design4.pdb
test output: test/pdbs

Auto Run MPNN Usage:
python auto_run_mpnn.py --original_pdb design4/original/design4.pdb --generated_pdbs_folder design4/pdbs_test --output_dir test_output --grouped_pdbs_dir test_output_grouped

SLURM Wrapper Usage:
sbatch run_autoProteinMPNN.srp <original_pdb> <generated_pdbs_dir> <output_dir>