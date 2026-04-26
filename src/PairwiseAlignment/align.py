#!/usr/bin/env python3
import sys
import argparse
from typing import Tuple, List, Optional, Dict

#Set a safety limit to prevent the matrix from eating all the RAM
MAX_MATRIX_SIZE = 16_000_000

#Standard BLOSUM62 Substitution Matrix
DEFAULT_BLOSUM62 = """\
#  Matrix made by matblas from blosum62.iij
#  * is a dummy
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

class SequenceAligner:
    def __init__(self, match_score: int = 1, mismatch_score: int = -1, gap_penalty: int = -1) -> None:
        self.match: int = match_score
        self.mismatch: int = mismatch_score
        self.gap: int = gap_penalty
        self.headers: List[str] = []
        self.sequences: List[str] = []
        self.substitution_matrix: Optional[Dict[str, Dict[str, int]]] = None

    def load(self, filename: str) -> None:
        """Reads a FASTA file, extracts headers and sequences, and validates formatting."""
        #Reset empty lists to prevent adding to old data if load() is called multiple times
        self.headers = []
        self.sequences = []
        
        try:
            with open(filename, "r") as infile:
                sequence_parts: List[str] = []
                for line in infile:
                    line = line.strip()
                    if not line:
                        continue 
                    
                    # If we hit a new header, save previous sequence and start new
                    if line.startswith(">"):
                        self.headers.append(line)
                        if sequence_parts:
                            self.sequences.append("".join(sequence_parts))
                            sequence_parts = [] 
                    else:
                        if not self.headers:
                            raise ValueError("Sequence data found before the first FASTA header.")
                        
                        #Strip out any spaces in the sequence
                        sequence_parts.append(line.replace(" ", ""))
                
                #Grab the last sequence in the file
                if sequence_parts:
                    self.sequences.append("".join(sequence_parts))
            
            #We validate some stuff post loading
            if not self.headers:
                raise ValueError(f"No FASTA headers were found in {filename}.")
                
            if len(self.headers) != len(self.sequences):
                raise ValueError("Mismatch between the number of FASTA headers and sequences.")
                
            if len(self.sequences) != 2:
                raise ValueError(f"Input FASTA file must contain exactly two sequences. Found {len(self.sequences)}.")
                
        except FileNotFoundError:
            raise FileNotFoundError(f"The file {filename} could not be found.")
        except Exception as e:
            if isinstance(e, ValueError):
                raise
            raise RuntimeError(f"An unexpected error occurred while reading the file: {e}")

    def _parse_matrix_lines(self, lines: List[str]) -> None:
        """Helper to parse matrix format from either a file or a hardcoded string."""
        self.substitution_matrix = {}
        amino_acids: List[str] = []
        for line in lines:
            #Skip any empty or comment lines
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            
            if not amino_acids:
                amino_acids = parts
                #For every amino acid, create a dict in the substitution matrix
                for aa in amino_acids:
                    self.substitution_matrix[aa] = {}
            else:
                #For every amino acid grab the scores in its row
                current_aa = parts[0]
                scores = [int(x) for x in parts[1:]]
                for i, score in enumerate(scores):
                    #Add that to the substitution matrix
                    col_aa = amino_acids[i]
                    self.substitution_matrix[current_aa][col_aa] = score

    def load_substitution_matrix(self, filepath: str) -> None:
        """Loads a scoring matrix from a user-provided file."""
        try:
            with open(filepath, 'r') as f:
                self._parse_matrix_lines(f.readlines())
            print(f"[Info] Successfully loaded custom substitution matrix from {filepath}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Substitution matrix file {filepath} not found.")
        except Exception as e:
            raise RuntimeError(f"Error loading matrix: {e}")

    def load_default_blosum62(self) -> None:
        """Loads the hardcoded BLOSUM62 matrix."""
        self._parse_matrix_lines(DEFAULT_BLOSUM62.split("\n"))
        print("[Info] Applying default BLOSUM62 substitution matrix.")

    def identify_sequence_type(self, sequence: str) -> str:
        """Simple check to guess if we are looking at DNA, RNA, or Protein."""
        if not sequence:
            return "Empty"
            
        seq = sequence.upper()
        unique_chars = set(seq)
        dna_chars = set("ACGTN")
        rna_chars = set("ACGUN")
        
        #Check if the unique characters of the sequence are DNA or RNA
        if unique_chars.issubset(dna_chars):
            return "DNA"
        elif unique_chars.issubset(rna_chars):
            return "RNA"
        else:
            #Or it's a protein or junk
            #Specifically for proteins, there are extra, commonly found characters
            protein_chars = set("ACDEFGHIKLMNPQRSTVWYBZX*")
            if unique_chars.issubset(protein_chars):
                return "Protein"
            else:
                return "Junk/Invalid"

    def _calculate_score(self, char1: str, char2: str) -> int:
        """Helper function to get the score between two characters."""
        #If a substitution matrix is found, it compares them based on that
        if self.substitution_matrix:
            c1, c2 = char1.upper(), char2.upper()
            try:
                return self.substitution_matrix[c1][c2]
            except KeyError:
                raise ValueError(f"Character '{c1}' or '{c2}' not found in substitution matrix.")
        
        #If not, then they are scored as DNA/RNA
        if char1 == char2:
            return self.match
        return self.mismatch

    def needleman_wunsch(self, seq1: str, seq2: str) -> Tuple[str, str, int]:
        """Global alignment using Needleman-Wunsch algorithm."""
        n, m = len(seq1), len(seq2)
        
        #Initialising a scoring matrix
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
        
        #Fill the first row and column with gap penalties, as in aligning a sequence with an empty sequence
        #This would produce a very bad alignment, full of gap penalties
        for i in range(1, n + 1):
            score_matrix[i][0] = i * self.gap
        for j in range(1, m + 1):
            score_matrix[0][j] = j * self.gap
        
        #Going through every pair on the grid, calculating the best score
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score_diag = score_matrix[i - 1][j - 1] + self._calculate_score(seq1[i - 1], seq2[j - 1])
                score_up = score_matrix[i - 1][j] + self.gap
                score_left = score_matrix[i][j - 1] + self.gap
                
                score_matrix[i][j] = max(score_diag, score_up, score_left)
                
        align1, align2 = "", ""
        i, j = n, m
        
        #Going backwards now, until we reach 0,0
        while i > 0 and j > 0:
            current_score = score_matrix[i][j]
            score_diff = self._calculate_score(seq1[i - 1], seq2[j - 1])
            
            #It checks the diagonal, top and left cells, going backwards until 0,0
            #Getting the highest score possible for every case  
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

        #At the end of one of the sequences, if there are still nucleotides in the other sequence
        #Gaps are added 
        while i > 0:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
        while j > 0:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
        
        #Reverse the sequences, as they were written backwards
        return align1[::-1], align2[::-1], score_matrix[n][m]

    def smith_waterman(self, seq1: str, seq2: str) -> Tuple[str, str, int]:
        """Local alignment using Smith-Waterman algorithm."""
        n, m = len(seq1), len(seq2)
        #Initialising a scoring matrix, with zeroes, instead of gap penalties
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
        
        max_score = 0
        max_pos = (0, 0)
        
        #Going over every cell in the matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                #It finds the best score for all three directions
                score_diag = score_matrix[i - 1][j - 1] + self._calculate_score(seq1[i - 1], seq2[j - 1])
                score_up = score_matrix[i - 1][j] + self.gap
                score_left = score_matrix[i][j - 1] + self.gap
                
                #That score can, at worst, be a 0. No negatives
                current_score = max(0, score_diag, score_up, score_left)
                score_matrix[i][j] = current_score
                
                if current_score > max_score:
                    max_score = current_score
                    max_pos = (i, j)
                
        align1, align2 = "", ""
        i, j = max_pos
        
        #Same procedure as in Needleman-Wunsch, going backwards
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

        #Returning the sequences reversed 
        return align1[::-1], align2[::-1], max_score

    def align(self, algorithm: str = "nw") -> Tuple[str, str, int]:
        """Main execution function to route to the correct algorithm."""
        seq1, seq2 = self.sequences[0], self.sequences[1]
        header1, header2 = self.headers[0], self.headers[1]

        # Check matrix size to avoid memory crashes
        if len(seq1) * len(seq2) > MAX_MATRIX_SIZE:
            print(f"[Warning] Sequences are very long ({len(seq1)} and {len(seq2)} chars).")
            print("This may consume a large amount of memory and take a long time to compute.\n")
        
        #Identifying if the sequences are DNA/RNA/Protein
        type1 = self.identify_sequence_type(seq1)
        type2 = self.identify_sequence_type(seq2)

        print(f"Seq 1 ({type1}): {header1[:50]}... ({len(seq1)} bp/aa)")
        print(f"Seq 2 ({type2}): {header2[:50]}... ({len(seq2)} bp/aa)")
        print("-" * 50)

        if type1 in ["Junk/Invalid", "Empty"] or type2 in ["Junk/Invalid", "Empty"]:
            raise ValueError("Cannot align. One or both sequences are empty or contain invalid characters.")
            
        if type1 != type2:
            print(f"[Warning] Aligning {type1} with {type2} is unusual, but proceeding...\n")
        else:
            print(f"Proceeding to align two {type1} sequences...")

        #If the sequences are protein and there is no provided substitution matrix, loading the default one
        if type1 == "Protein" and type2 == "Protein" and not self.substitution_matrix:
            self.load_default_blosum62()
        print()

        if algorithm == "nw":
            print("Running Needleman-Wunsch (Global Alignment)...")
            aligned_1, aligned_2, score = self.needleman_wunsch(seq1, seq2)
        elif algorithm == "sw":
            print("Running Smith-Waterman (Local Alignment)...")
            aligned_1, aligned_2, score = self.smith_waterman(seq1, seq2)
        else:
            raise ValueError(f"Unknown algorithm selected: {algorithm}")

        self.print_alignment(aligned_1, aligned_2, score)
        return aligned_1, aligned_2, score

    def print_alignment(self, align1: str, align2: str, score: int, line_length: int = 60) -> None:
        """Formats and prints the final alignment in readable chunks."""
        print(f"\nFinal Alignment Score: {score}")
        print("=" * 60)
        
        #Printing a connection if they are the same, else nothing between the sequences
        match_line = "".join("|" if a1 == a2 and a1 != "-" else " " for a1, a2 in zip(align1, align2))
        
        for i in range(0, len(align1), line_length):
            chunk1 = align1[i:i+line_length]
            chunk_match = match_line[i:i+line_length]
            chunk2 = align2[i:i+line_length]
            
            print(f"S1: {chunk1}")
            print(f"    {chunk_match}")
            print(f"S2: {chunk2}\n")
        print("=" * 60)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pairwise Sequence Aligner (Needleman-Wunsch & Smith-Waterman)")
    
    #Required arguments
    parser.add_argument("fasta", help="Path to the input FASTA file containing exactly 2 sequences.")
    parser.add_argument("algorithm", nargs="?", choices=["nw", "sw"], default="nw", 
                        help="Alignment algorithm: nw (global) or sw (local). Default is nw.")
    
    #Optional flags
    parser.add_argument("-x", "--matrix", help="Path to a substitution matrix file (e.g., blosum62.txt) for proteins.")
    parser.add_argument("--match", type=int, default=2, help="Score for a match (default: 2)")
    parser.add_argument("--mismatch", type=int, default=-1, help="Penalty for a mismatch (default: -1)")
    parser.add_argument("--gap", type=int, default=-2, help="Penalty for a gap (default: -2)")

    args = parser.parse_args()
    aligner = SequenceAligner(match_score=args.match, mismatch_score=args.mismatch, gap_penalty=args.gap)

    try:
        if args.matrix:
            aligner.load_substitution_matrix(args.matrix)

        aligner.load(args.fasta)
        aligner.align(algorithm=args.algorithm)
        
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nAlignment cancelled by user.", file=sys.stderr)
        sys.exit(1)