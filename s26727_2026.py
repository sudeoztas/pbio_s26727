# Album: s26727
# Date: 2026
# Description: Random DNA sequence generator in FASTA format with extended bioinformatics features.

import random
import csv
import matplotlib.pyplot as plt

IUPAC_MAP = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}

STANDARD_NUCLEOTIDES = ['A', 'C', 'G', 'T']


def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length using only standard nucleotides."""
    return ''.join(random.choices(STANDARD_NUCLEOTIDES, k=length))


def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary of sequence statistics.
    Keys: 'A', 'C', 'G', 'T' (float, %), 'GC' (float, %), 'gc_ratio_A' (float, %).
    Only counts standard nucleotides, ignores embedded name characters and IUPAC ambiguous chars.
    """
    seq = ''.join(c for c in sequence if c.upper() in 'ACGT')
    n = len(seq)
    stats = {}
    for nuc in STANDARD_NUCLEOTIDES:
        stats[nuc] = (seq.upper().count(nuc) / n) * 100
    stats['GC'] = stats['G'] + stats['C']
    stats['gc_ratio_A'] = stats['GC']
    return stats


def insert_name(sequence: str, name: str) -> str:
    """Inserts a name at a random position in the sequence. Name written in lowercase letters."""
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


def format_fasta_record(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a single FASTA record as a string without EOF marker.
    Header starts with '>', ID and description separated by a space.
    Sequence is broken into lines of exactly line_width characters.
    """
    header = f">{seq_id}"
    if description:
        header += f" {description}"
    lines = [header]
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])
    return '\n'.join(lines) + '\n'


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a complete single-record FASTA file content as a string.
    Appends the required # EOF_1 marker at the very end.
    """
    return format_fasta_record(seq_id, description, sequence, line_width) + "# EOF_1\n"


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Gets an integer from the user in a specified range.
    In case of an error, repeats the question without raising an exception.
    """
    while True:
        val = input(prompt)
        try:
            n = int(val)
            if min_val <= n <= max_val:
                return n
            raise ValueError
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


def validate_id(prompt: str) -> str:
    """Gets a sequence ID from the user.
    Validates that ID is non-empty and contains no whitespace characters.
    """
    while True:
        seq_id = input(prompt)
        if seq_id and ' ' not in seq_id and '\t' not in seq_id:
            return seq_id
        print("Error: ID cannot be empty or contain whitespace.")