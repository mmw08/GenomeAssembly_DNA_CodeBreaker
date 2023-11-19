# GenomeAssembly_DNA_CodeBreaker
## Introduction

This Python script, named DNA_CODE_Breaker.py, is designed for DNA sequence assembly. The assembly process can be performed with or without considering mismatches in the reads. The script utilizes a De Bruijn graph approach for exact reads assembly,
and a trie-based method for reads containing mismatches.

## Usage

The script is designed to be run from the command line with the following syntax:

```bash
python DNA_Code_Breaker.py reads_file reference_file Mismatch_Flag
```


**Usage**
The script is designed to be run from the command line with the following syntax:

reads_file: Path to the file containing DNA reads.<br>
reference_file: Path to the file containing the reference DNA sequence.<br>
Mismatch_Flag: Specify whether to consider mismatches in the reads (Y for yes, N for no).

## De Bruijin use
Exact Reads Assembly (Mismatch_Flag = N).<br>
If Mismatch_Flag is set to N, the script performs assembly of exact reads. The De Bruijn graph approach is used to construct the graph, identify source and sink nodes, find cycles, and assemble the genome sequence.<br>
**Example to run file:**
```bash
python DNA_Code_Breaker.py reads.txt reference.txt N
```

## Trie Use
Reads Assembly with Mismatches (Mismatch_Flag = Y).<br>
If Mismatch_Flag is set to Y, the script considers mismatches in the reads. The script constructs a trie using the reads file and performs trie matching against the reference genome. The resulting assembly includes positions for each matched pattern.<br>
**Example to run file:**
```bash
python DNA_Code_Breaker.py reads.txt reference.txt Y
```

## Output
The assembled genome sequence is written to the an output file **_assembledreads.txt_**. For the case of mismatches, positions of the matched patterns are indicated in the output.<br>
The Mismatch_Flag parameter only accepts Y (yes) or N (no). Any other input will result in an error message.
