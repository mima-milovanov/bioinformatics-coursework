"""
DNA FASTA Analysis Tool
Answers questions about sequences in a multi-FASTA file:
1. Number of records
2. Sequence lengths (longest/shortest)
3. Open Reading Frames (ORFs) in a given reading frame
4. Repeated substrings of length n

Usage: python fasta_sequence_analysis.py <fasta_file> [frame] [n] [target_id]
Example: python fasta_sequence_analysis.py dna2.fasta 3 6 "gi|142022655|gb|EQ086233.1|16"
"""

import sys
from collections import Counter


def parse_fasta(filename):
    """ Parse FASTA file"""
    sequences = []
    current_id = None
    seq_parts = []

    with open(filename) as f:
        for line in f:
            line = line.strip()  # uklanja whitespace i razmake
            if not line:  # ako je prazna linija preskoci tu iteraciju
                continue
            if line.startswith(">"):  # provera da li je to header, ako jeste onda je to nova seq
                if current_id:  # ako smo vec imali prethodnu seq, moramo da je sacuvamo
                    sequences.append((current_id, "".join(seq_parts)))
                current_id = line[1:].split()[0]  # od hedera >seq1 description here skida > pa deli string po razmacima ['seq1', 'description', 'here'] i uzima index 0
                seq_parts = []  # resetovanje liste seq, pravimo novu listu za sledecu seq
            else:
                seq_parts.append(line.upper())  # upper() osigurava da su svi nukleotidi: ATGC a ne atgc

        if current_id:
            sequences.append((current_id, "".join(seq_parts)))  # da se sacuva poslednja seq

    return sequences


def count_records(sequences):
    """ Counts FASTA records"""
    return len(sequences)


def sequence_lengths(sequences):
    """"Computes seq lengths (longest/shortest/their IDs)"""
    lengths = {seq_id: len(seq) for seq_id, seq in sequences}  # dictionary comprehension

    max_len = max(lengths.values())
    min_len = min(lengths.values())

    longest_seq = [seq_id for seq_id, length in lengths.items() if length == max_len]
    shortest_seq = [seq_id for seq_id, length in lengths.items() if length == min_len]

    return lengths, max_len, min_len, longest_seq, shortest_seq


START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs(sequence, frame):
    """Returns list of open reading frame - ORFs:
    (start_position, length, sequence)"""
    orfs = []
    i = frame - 1  # Pošto Python koristi indeksiranje od 0, a biološki okviri se broje od 1

    while i + 3 <= len(sequence):
        codon = sequence[i:i + 3]
        if codon == START_CODON:
            j = i + 3  # krece NAKON start kodona
            while j + 3 <= len(sequence):
                stop = sequence[j:j + 3]
                if stop in STOP_CODONS:
                    orf_seq = sequence[i:j + 3]
                    orf_len = j + 3 - i  # include stop codon
                    orfs.append((i + 1, orf_len, orf_seq))  # 1-based start
                    break
                j += 3
        i += 3

    return orfs


def longest_orf_in_file(sequences, frame):
    """Finds the longest open reading frame (ORF) across multiple sequences"""
    longest_len = 0
    longest_info = None

    for seq_id, seq in sequences:
        orfs = find_orfs(seq, frame)

        for start, length, orf_seq in orfs:
            if length > longest_len:
                longest_len = length
                longest_info = (seq_id, start, length, orf_seq)

    return longest_info


def longest_orf_in_sequence(sequence, frame):
    """Returns tuple of the longest ORF in one sequence"""
    orfs = find_orfs(sequence, frame)
    if not orfs:
        return None
    return max(orfs, key=lambda x: x[1])  # vrati najveci element iz liste orfs po duzini (x[1])


def find_repeats(sequences, n):
    """Counts occurrences of all substrings of length n across all sequences"""
    counts = Counter()  # Counter je specijalni recnik koji automatski broji pojavljivanja

    for _, seq in sequences:  # raspakuje tuple (seq_id, seq)
        for i in range(len(seq) - n + 1):  # prolazi kroz sve moguce pozicije
            repeat = seq[i:i + n]
            counts[repeat] += 1

    return counts


def most_frequent_repeat(sequences, n):
    """Finds the most frequently occurring substring of length n across all sequences"""
    counts = find_repeats(sequences, n)
    repeats_only = {k: v for k, v in counts.items() if v > 1}  # samo oni koji se pojavljuju vise od jednom
    if not repeats_only:
        return None
    repeat = max(repeats_only, key=repeats_only.get)
    return repeat, repeats_only[repeat]


def main():
    if len(sys.argv) < 2:
        print("Usage: python fasta_sequence_analysis.py <fasta_file> [frame] [n] [target_id]")
        return

    fasta_file = sys.argv[1]
    frame = int(sys.argv[2]) if len(sys.argv) > 2 else 1  # default frame = 1
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 6      # default n = 6
    target_id = sys.argv[4] if len(sys.argv) > 4 else None

    sequences = parse_fasta(fasta_file)

    # Q1 - Count FASTA records
    print("Number of records:", count_records(sequences))

    # Q2 - Sequence lengths
    lengths, max_len, min_len, longest_seq, shortest_seq = sequence_lengths(sequences)
    print("\nSequence lengths:")
    for seq_id, length in lengths.items():
        print(seq_id, length)

    print("\nLongest sequence length: ", max_len)
    print("Longest sequence IDs: ", longest_seq)
    print("\nShortest sequence length:", min_len)
    print("Shortest sequence IDs: ", shortest_seq)

    # Q3 - Longest ORF in file
    result = longest_orf_in_file(sequences, frame)
    if result:
        seq_id, start, length, seq = result
        print(f"\nLongest ORF in file (frame {frame}):")
        print("Sequence ID: ", seq_id)
        print("Start position: ", start)
        print("Length: ", length)

    # Q3 - Longest ORF for a specific sequence ID
    if target_id:
        for seq_id, seq in sequences:
            if seq_id == target_id:
                result = longest_orf_in_sequence(seq, frame)
                if result:
                    start, length, orf_seq = result
                    print(f"\nLongest ORF for {target_id} (frame {frame}):")
                    print("Start position:", start)
                    print("Length:", length)
                else:
                    print(f"No ORFs found for {target_id} in frame {frame}")
                break
        else:
            print(f"ID '{target_id}' not found in file")

    # Q4 - Most frequent repeat of length n
    result = most_frequent_repeat(sequences, n)
    if result:
        repeat, freq = result
        print(f"\nMost frequent repeat of length {n}:")
        print("Repeat:", repeat)
        print("Frequency:", freq)
    else:
        print(f"\nNo repeats of length {n} found.")

    # Koliko razlicitih sekvenci se pojavljuje Max puta
    counts = find_repeats(sequences, n)
    repeats_only = {k: v for k, v in counts.items() if v > 1}
    max_freq = max(repeats_only.values())
    seqs_with_max = [k for k, v in repeats_only.items() if v == max_freq]
    print(f"Number of sequences occurring max ({max_freq}) times: {len(seqs_with_max)}")


if __name__ == "__main__":
    main()