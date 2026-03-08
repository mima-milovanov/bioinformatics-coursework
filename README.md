# DNA Sequence Analyzer

This project was developed for the final exam of the [Python for Genomic Data Science](https://www.coursera.org/specializations/genomic-data-science) course, part of the Genomic Data Science Specialization on Coursera.

## Overview

The tool processes multi-FASTA files and performs bioinformatic analysis to answer specific questions about DNA sequences. It covers four areas: sequence parsing and record counting, sequence length analysis, Open Reading Frame (ORF) identification, and repeated substring detection.

## Features

**Record counting** identifies the total number of sequences in a FASTA file by detecting headers starting with the `>` symbol.

**Sequence length analysis** computes the length of every sequence and reports the longest and shortest along with their identifiers.

**ORF identification** searches for Open Reading Frames starting with `ATG` and ending with `TAA`, `TAG`, or `TGA`. It supports all three forward reading frames and can find the longest ORF across the entire file or within a specific sequence. Starting positions are reported using 1-based indexing. A sequence is only considered a valid ORF if it contains both a start codon and an in-frame stop codon.

**Repeat detection** finds all substrings of length `n` that occur more than once, including overlapping repeats, and reports the most frequent one along with its occurrence count.

## Requirements

- Python 3.x
- No external libraries required

## Installation

```bash
git clone https://github.com/mima-milovanov/dna-sequence-analyzer.git
cd dna-sequence-analyzer
```

## Usage

```bash
python fasta_sequence_analysis.py <fasta_file> [frame] [n] [target_id]
```

| Argument | Description | Default |
|---|---|---|
| `fasta_file` | Path to input FASTA file | required |
| `frame` | Forward reading frame (1, 2, or 3) | 1 |
| `n` | Length of repeats to search for | 6 |
| `target_id` | Sequence identifier to analyze individually | None |

## Examples

```bash
# Basic usage
python fasta_sequence_analysis.py dna.fasta

# Find ORFs in reading frame 2, repeats of length 12
python fasta_sequence_analysis.py dna.fasta 2 12

# Find longest ORF for a specific sequence
python fasta_sequence_analysis.py dna.fasta 1 6 "gi|142022655|gb|EQ086233.1|16"
```

## Output Example

```
Number of records: 18

Sequence lengths:
gi|142022655|gb|EQ086233.1|16 4804
...
Longest sequence length:  4894
Longest sequence IDs:  ['gi|142022655|gb|EQ086233.1|255']
Shortest sequence length: 115
Shortest sequence IDs:  ['gi|142022655|gb|EQ086233.1|346']

Longest ORF in file (frame 1):
Sequence ID:  gi|142022655|gb|EQ086233.1|45
Start position:  385
Length:  2394

Most frequent repeat of length 12:
Repeat: CATTCGCCATTC
Frequency: 10
```
