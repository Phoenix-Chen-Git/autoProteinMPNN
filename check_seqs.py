
from Bio import PDB
from Bio.SeqUtils import seq1

def get_seq(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)
    seqs = []
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                if PDB.is_aa(residue):
                    seq.append(seq1(residue.get_resname()))
            seqs.append("".join(seq))
    return seqs

print("Original:", get_seq("test/origin_pdb/design4.pdb"))
print("Designed:", get_seq("test/pdbs/design4_0.pdb"))
