#!/usr/bin/env python3
import sys
import argparse

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
            print(f"Error: The file '{filename}' could not be found.")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred while reading the file: {e}")
            sys.exit(1)

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
            print(f"[Info] Successfully loaded substitution matrix from '{filepath}'")
        except FileNotFoundError:
            print(f"Error: Substitution matrix file '{filepath}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading matrix: {e}")
            sys.exit(1)

    def identify_sequence_type(self, sequence):
        if not sequence:
            return "Empty"
            
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

    def align(self, algorithm="nw"):
        if len(self.sequences) < 2:
            print("Error: The loaded data must contain at least two sequences.")
            sys.exit(1)

        seq1, seq2 = self.sequences[0], self.sequences[1]
        header1, header2 = self.headers[0], self.headers[1]

        # 1. Length/Memory safety check (Warns if matrix cells > 16 million)
        if len(seq1) * len(seq2) > 16000000:
            print(f"[Warning] Sequences are very long ({len(seq1)} and {len(seq2)} chars).")
            print("This may consume a large amount of memory and take a long time to compute.\n")

        # 2. Sequence Type Validation
        type1 = self.identify_sequence_type(seq1)
        type2 = self.identify_sequence_type(seq2)

        print(f"Seq 1 ({type1}): {header1[:50]}... ({len(seq1)} bp/aa)")
        print(f"Seq 2 ({type2}): {header2[:50]}... ({len(seq2)} bp/aa)")
        print("-" * 50)

        if type1 in ["Junk/Invalid", "Empty"] or type2 in ["Junk/Invalid", "Empty"]:
            print("Error: Cannot align. One or both sequences are empty or contain invalid characters.")
            sys.exit(1)
            
        if type1 != type2:
            print(f"[Warning] Aligning {type1} with {type2} is unusual, but proceeding...\n")
        else:
            print(f"Proceeding to align two {type1} sequences...\n")

        # 3. Execution
        if algorithm == "nw":
            print("Running Needleman-Wunsch (Global Alignment)...")
            aligned_1, aligned_2, score = self.needleman_wunsch(seq1, seq2)
        elif algorithm == "sw":
            print("Running Smith-Waterman (Local Alignment)...")
            aligned_1, aligned_2, score = self.smith_waterman(seq1, seq2)

        self.print_alignment(seq1, seq2, aligned_1, aligned_2, score)
        return aligned_1, aligned_2, score

    def print_alignment(self, seq1, seq2, align1, align2, score, line_length=60):
        print(f"\nFinal Alignment Score: {score}")
        print("=" * 60)
        
        # Build the match line
        match_line = "".join("|" if a1 == a2 and a1 != "-" else " " for a1, a2 in zip(align1, align2))
        
        # Print in chunks for readability
        for i in range(0, len(align1), line_length):
            chunk1 = align1[i:i+line_length]
            chunk_match = match_line[i:i+line_length]
            chunk2 = align2[i:i+line_length]
            
            print(f"S1: {chunk1}")
            print(f"    {chunk_match}")
            print(f"S2: {chunk2}\n")
        print("=" * 60)


# --- Professional Command Line Interface ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pairwise Sequence Aligner (Needleman-Wunsch & Smith-Waterman)")
    
    # Required Arguments
    parser.add_argument("fasta", help="Path to the input FASTA file containing at least 2 sequences.")
    
    # Optional Algorithm Choice
    parser.add_argument("-a", "--algorithm", choices=['nw', 'sw'], default='nw', 
                        help="Algorithm to use: 'nw' (Global) or 'sw' (Local). Default is 'nw'.")
    
    # Optional Substitution Matrix
    parser.add_argument("-x", "--matrix", help="Path to a substitution matrix file (e.g., blosum62.txt) for proteins.")
    
    # Optional Custom Scoring Parameters
    parser.add_argument("--match", type=int, default=2, help="Score for a match (default: 2)")
    parser.add_argument("--mismatch", type=int, default=-1, help="Penalty for a mismatch (default: -1)")
    parser.add_argument("--gap", type=int, default=-2, help="Penalty for a gap (default: -2)")

    args = parser.parse_args()

    # Initialize the aligner with custom or default arguments
    aligner = SequenceAligner(match_score=args.match, mismatch_score=args.mismatch, gap_penalty=args.gap)

    # Load matrix if provided
    if args.matrix:
        aligner.load_substitution_matrix(args.matrix)

    # Load the FASTA file and execute
    aligner.load(args.fasta)
    aligner.align(algorithm=args.algorithm)