#!/usr/bin/env python3
import sys

class SequenceAligner:
    def __init__(self, match_score=1, mismatch_score=-1, gap_penalty=-1):
        self.match = match_score
        self.mismatch = mismatch_score
        self.gap = gap_penalty
        self.headers = []
        self.sequences = []
        self.substitution_matrix = None

    def load(self, filename):
        try:
            with open(filename, "r") as infile:
                line = infile.readline()
                sequence = ""

                if line == "":
                    raise ValueError(f"The file {filename} is empty.")
                    
                while line != "":
                    if line.startswith(">"):
                        self.headers.append(line.strip())
                        if sequence != "":
                            self.sequences.append(sequence)
                            sequence = ""
                    else:
                        sequence += line.strip().replace(" ", "")
                    line = infile.readline()
                
                if sequence != "":
                    self.sequences.append(sequence)
                
        except FileNotFoundError:
            print(f"Error: The file {filename} could not be found. Please check the path")
        except PermissionError:
            print(f"Error: You do not have permission to read the file {filename}.")
        except IsADirectoryError:
            print(f"Error: The path {filename} points to a directory, not a file")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def load_substitution_matrix(self, filepath):
        self.substitution_matrix = {}
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                
            amino_acids = []
            for line in lines:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if not amino_acids:
                    amino_acids = parts
                    for aa in amino_acids:
                        self.substitution_matrix[aa] = {}
                else:
                    current_aa = parts[0]
                    scores = [int(x) for x in parts[1:]]
                    for i, score in enumerate(scores):
                        col_aa = amino_acids[i]
                        self.substitution_matrix[current_aa][col_aa] = score
            print(f"Successfully loaded substitution matrix from {filepath}\n")
        except Exception as e:
            print(f"Error loading matrix: {e}\n")
            self.substitution_matrix = None

    def identify_sequence_type(self, sequence):
        seq = sequence.upper()
        unique_chars = set(seq)
        dna_chars = set("ACGTN")
        rna_chars = set("ACGUN")
        
        if unique_chars.issubset(dna_chars):
            return "DNA"
        elif unique_chars.issubset(rna_chars):
            return "RNA"
        else:
            protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
            if unique_chars.issubset(protein_chars):
                return "Protein"
            else:
                return "Junk/Invalid"

    def _calculate_score(self, char1, char2):
        if self.substitution_matrix:
            try:
                return self.substitution_matrix[char1.upper()][char2.upper()]
            except KeyError:
                return self.mismatch
                
        if char1 == char2:
            return self.match
        return self.mismatch

    def needleman_wunsch(self, seq1, seq2):
        n, m = len(seq1), len(seq2)
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
        
        for i in range(1, n + 1):
            score_matrix[i][0] = i * self.gap
        for j in range(1, m + 1):
            score_matrix[0][j] = j * self.gap
            
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score_diag = score_matrix[i - 1][j - 1] + self._calculate_score(seq1[i - 1], seq2[j - 1])
                score_up = score_matrix[i - 1][j] + self.gap
                score_left = score_matrix[i][j - 1] + self.gap
                score_matrix[i][j] = max(score_diag, score_up, score_left)
                
        align1, align2 = "", ""
        i, j = n, m
        
        while i > 0 and j > 0:
            current_score = score_matrix[i][j]
            score_diff = self._calculate_score(seq1[i - 1], seq2[j - 1])
                
            if current_score == score_matrix[i - 1][j - 1] + score_diff:
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                i -= 1
                j -= 1
            elif current_score == score_matrix[i - 1][j] + self.gap:
                align1 += seq1[i - 1]
                align2 += "-"
                i -= 1
            else:
                align1 += "-"
                align2 += seq2[j - 1]
                j -= 1
                
        while i > 0:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        while j > 0:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
            
        return align1[::-1], align2[::-1], score_matrix[n][m]

    def smith_waterman(self, seq1, seq2):
        n, m = len(seq1), len(seq2)
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
        
        max_score = 0
        max_pos = (0, 0)
            
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score_diag = score_matrix[i - 1][j - 1] + self._calculate_score(seq1[i - 1], seq2[j - 1])
                score_up = score_matrix[i - 1][j] + self.gap
                score_left = score_matrix[i][j - 1] + self.gap
                
                current_score = max(0, score_diag, score_up, score_left)
                score_matrix[i][j] = current_score
                
                if current_score > max_score:
                    max_score = current_score
                    max_pos = (i, j)
                
        align1, align2 = "", ""
        i, j = max_pos
        
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            current_score = score_matrix[i][j]
            score_diff = self._calculate_score(seq1[i - 1], seq2[j - 1])
                
            if current_score == score_matrix[i - 1][j - 1] + score_diff:
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                i -= 1
                j -= 1
            elif current_score == score_matrix[i - 1][j] + self.gap:
                align1 += seq1[i - 1]
                align2 += "-"
                i -= 1
            else:
                align1 += "-"
                align2 += seq2[j - 1]
                j -= 1
            
        return align1[::-1], align2[::-1], max_score

    def align(self, algorithm="needleman_wunsch"):
        if len(self.sequences) < 2:
            print("Error: The loaded data does not contain at least two sequences.")
            return None

        seq1, seq2 = self.sequences[0], self.sequences[1]
        header1, header2 = self.headers[0], self.headers[1]

        type1 = self.identify_sequence_type(seq1)
        type2 = self.identify_sequence_type(seq2)

        print(f"Sequence 1 ({type1}): {header1}")
        print(f"Sequence 2 ({type2}): {header2}")
        print("-" * 30)

        if type1 == "Junk/Invalid" or type2 == "Junk/Invalid":
            print("Error: Cannot align. One or both sequences contain invalid characters.")
            return None
            
        if type1 != type2:
            print(f"Warning: Aligning {type1} with {type2} is unusual, but proceeding...\n")
        else:
            print(f"Proceeding to align two {type1} sequences...\n")

        if algorithm == "needleman_wunsch":
            print("Running Needleman-Wunsch (Global Alignment)...")
            aligned_1, aligned_2, score = self.needleman_wunsch(seq1, seq2)
        elif algorithm == "smith_waterman":
            print("Running Smith-Waterman (Local Alignment)...")
            aligned_1, aligned_2, score = self.smith_waterman(seq1, seq2)
        else:
            print(f"Error: Unknown algorithm '{algorithm}'")
            return None

        self.print_alignment(seq1, seq2, aligned_1, aligned_2, score)
        return aligned_1, aligned_2, score

    def print_alignment(self, seq1, seq2, align1, align2, score):
        print(f"Alignment Score: {score}")
        print("-" * 30)
        print(align1)
        match_line = "".join("|" if a1 == a2 and a1 != "-" else " " for a1, a2 in zip(align1, align2))
        print(match_line)
        print(align2)
        print("-" * 30)


# --- Command Line Interface ---
if __name__ == "__main__":
    # Check if the user provided enough arguments
    if len(sys.argv) < 2:
        print("Usage: ./align.py <filename> [algorithm]")
        print("Algorithms: needleman_wunsch (default), smith_waterman")
        sys.exit(1)

    # First argument is the FASTA file
    fasta_file = sys.argv[1]
    
    # Second argument is the algorithm (defaulting to Needleman-Wunsch)
    algorithm_choice = "needleman_wunsch"
    if len(sys.argv) >= 3:
        # Convert to lowercase and replace hyphens with underscores to match method names
        algorithm_choice = sys.argv[2].lower().replace("-", "_")

    # Validate algorithm choice
    if algorithm_choice not in ["needleman_wunsch", "smith_waterman"]:
        print(f"Error: Unknown algorithm '{algorithm_choice}'.")
        print("Please choose either 'needleman_wunsch' or 'smith_waterman'.")
        sys.exit(1)

    # Initialize the aligner
    aligner = SequenceAligner(match_score=2, mismatch_score=-1, gap_penalty=-2)

    # Optional: If you want to automatically load BLOSUM62 for proteins, 
    # you can uncomment the line below and ensure the file is in the same directory.
    # aligner.load_substitution_matrix("blosum62.txt")

    # Load the file and execute the alignment
    aligner.load(fasta_file)
    aligner.align(algorithm=algorithm_choice)