# Importing Dependencies
# import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    a = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    alignment = a.align(seq1, seq2)
    alignment_score = 4.0
    backtrace_matrix = np.array([[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [-np.inf, 0.0, 2.0, 2, 2], [-np.inf, 1, 0, 0, 0], [-np.inf, 1, 1, 0, 0]])
    print(a._back)
    assert (np.array_equal(a._back, backtrace_matrix) == True)
    assert (alignment_score == alignment[0])

    
def test_nw_backtrace(): 
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    a = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    assert ( a.align(seq1, seq2) == (4.0, 'MYQR', 'M-QR'))

def test_nw_seq3_seq4():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    a = NeedlemanWunsch(sub_matrix_file="substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    assert ( a.align(seq3, seq4) == (17.0, 'MAVHQLIRRP', 'M---QLIRHP'))

test_nw_alignment()