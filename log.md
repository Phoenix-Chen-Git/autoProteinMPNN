# 2025.12.8
We have original pdb, generated pbds.
(two input --original_pdb --generated_pdbs)


We want, group the generated pdbs, get their chain and positions to edit.
(analyze_pdbs.py --original_pdb --generated_pdbs to --output_dir)


design for all of them.
(--original_pdb --generated_pdbs --output_dir)


separate the fa into fas.


trimm using start and end.


## analyze_pdbs.py
usage:
python analyze_pdbs.py --original_pdb test/origin_pdb/design4.pdb --generated_pdbs_folder test/pdbs/ --output_dir *


## MPNN4_customize.sh
usage:
sbatch MPNN4_customize.sh [--pdb_dir DIR] [--output_dir DIR] [--chains "A B"] [--design_positions "..."]


## auto_run_mpnn.py
usage:
python auto_run_mpnn.py --original_pbd test/origin_pdb/design4.pdb --generated_pdbs_folder test/pdbs/ --output_dir test_output
this script's function can cover the use of analyze_pdbs.py