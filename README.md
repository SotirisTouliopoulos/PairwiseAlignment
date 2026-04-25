
https://teaching.healthtech.dtu.dk/22118/index.php/Pairwise_alignment

Description
Aligning sequences is of great importance in bioinformatics. Many discoveries are based on finding sequences that align to each other. 
Evolution theory and phylogeny are based on sequence alignments. This project is about implementing a well-known algorithm for aligning two sequences, 
i.e. finding where they match in an optimal fashion.

You must choose to implement either:

Smith-Waterman alignment where the goal is to find the best local alignment of the two sequences given as input, i.e. the optimal alignment that covers most/best of both sequences.
Needleman-Wunsch alignment where the goal is to find the best global alignment of the two sequences given as input, i.e. the optimal alignment that covers all of at least one sequence.
Or both if you are cool :-)

Input and output
The input is just a fasta file with two sequences, that should be aligned.
The output should be the the best alignment with clear notation where it is in both sequence inputs.
Note: Pairwise alignment works for both DNA and protein sequences.

Examples of program execution:

```
align.py <fastafile>
align.py fastafile.fsa
```

https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
https://books.google.dk/books?id=bvY21DGa1OwC&pg=PA87&lpg=PA87&dq=smith-waterman&source=web&ots=gJTIZXMkqv&sig=AO0TtuhrNFaH0ZuKIg2TUXHmqww&hl=en#v=onepage&q=smith-waterman&f=false
https://en.wikipedia.org/wiki/BLOSUM
https://en.wikipedia.org/wiki/Point_accepted_mutation



# Installation

### Install package dependencies
```
pip install pytest
```

### Install package
```
pip install -e .
```

# Uninstall package
```
pip uninstall PairwiseAlignment
```
