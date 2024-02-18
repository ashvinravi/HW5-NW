# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    a = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    gg_score = a.align(gg_seq, hs_seq)[0]
    mm_score = a.align(mm_seq, hs_seq)[0]
    br_score = a.align(br_seq, hs_seq)[0]
    tt_score = a.align(tt_seq, hs_seq)[0]
    species = ["Gallus gallus", "Mus musculus", "Balaeniceps rex", "Tursiops truncatus"]
    scores = [gg_score, mm_score, br_score, tt_score]
    indices = np.argsort(scores)[::-1]
    for i in indices:
        print(species[i] + ": " + str(scores[i]))

    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass
    

if __name__ == "__main__":
    main()
