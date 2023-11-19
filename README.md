# GenomeAssembly_DNA_CodeBreaker
## Introduction

This Python script, named DNA_CODE_Breaker.py, is designed for DNA sequence assembly. The assembly process can be performed with or without considering mismatches in the reads. The script utilizes a De Bruijn graph approach for exact reads assembly,
and a trie-based method for reads containing mismatches.

## Usage

The script is designed to be run from the command line with the following syntax:

```bash
python DNA_CODE_Breaker.py reads_file reference_file Mismatch_Flag



**Usage**
The script is designed to be run from the command line with the following syntax:

reads_file: Path to the file containing DNA reads.
reference_file: Path to the file containing the reference DNA sequence.
Mismatch_Flag: Specify whether to consider mismatches in the reads (Y for yes, N for no).

## De Bruijin use
Exact Reads Assembly (Mismatch_Flag = N)
If Mismatch_Flag is set to N, the script performs assembly of exact reads. The De Bruijn graph approach is used to construct the graph, identify source and sink nodes, find cycles, and assemble the genome sequence.

## Trie Use
Reads Assembly with Mismatches (Mismatch_Flag = Y)
If Mismatch_Flag is set to Y, the script considers mismatches in the reads. The script constructs a trie using the reads file and performs trie matching against the reference genome. The resulting assembly includes positions for each matched pattern.

## Output
The assembled genome sequence is written to the an output file _"assembledreads.txt"_. For the case of mismatches, positions of the matched patterns are indicated in the output.

Examples
Assemble exact reads:

bash
Copy code
python DNA_CODE.py exact_reads.txt reference_sequence.txt N assembled_genome.txt
Assemble reads with mismatches:

bash
Copy code
python DNA_CODE.py reads_with_mismatches.txt reference_sequence.txt Y assembled_genome.txt
Note
The Mismatch_Flag parameter only accepts Y (yes) or N (no). Any other input will result in an error message.
