
class NeedlemanWunsch:
    def __init__(self, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
        """
        Initializes the scoring system for the alignment.
        """
        pass

    def read_fasta(self, file_path):
        """
        Reads a FASTA file and extracts exactly two sequences.
        """
        pass

    def initialize_matrices(self, seq1, seq2):
        """
        Creates the scoring matrix (filled with zeros initially) 
        and applies the accumulating gap penalties to the first row and column.
        """
        pass

    def fill_matrix(self, seq1, seq2, score_matrix):
        """
        Iterates through the scoring matrix, calculating the max score 
        for each cell based on the diagonal (match/mismatch), left (gap), 
        and up (gap) values.
        
        score_matrix from the "initialize_matrices()" function
        """
        pass

    def traceback(self, seq1, seq2, score_matrix):
        """
        Starts at the bottom-right corner of the filled matrix and works 
        backwards to the top-left, building the aligned sequences 
        and inserting gaps ('-') where necessary.
        
        score_matrix from the "fill_matrix()" function
        """
        pass

    def align(self, fasta_file_path):
        """
        The main runner method.
        1. Calls read_fasta()
        2. Calls initialize_matrices()
        3. Calls fill_matrix()
        4. Calls traceback()
        5. Returns the final aligned sequences
        """
        pass