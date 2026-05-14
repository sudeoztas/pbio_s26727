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


def generate_sequence_weighted(length: int, weights: dict) -> str:
    """Returns a random DNA sequence using custom nucleotide probabilities.
    weights is a dict with keys A, C, G, T and float percentage values.
    """
    w = [weights[n] for n in STANDARD_NUCLEOTIDES]
    return ''.join(random.choices(STANDARD_NUCLEOTIDES, weights=w, k=length))


def generate_sequence_iupac(length: int, ambiguous_prob: float) -> str:
    """Returns a random DNA sequence that may contain IUPAC ambiguous characters.
    Each position has ambiguous_prob probability of being an ambiguous character.
    Reference: https://www.bioinformatics.org/sms/iupac.html
    """
    ambiguous_chars = list(IUPAC_MAP.keys())
    result = []
    for _ in range(length):
        if random.random() < ambiguous_prob:
            result.append(random.choice(ambiguous_chars))
        else:
            result.append(random.choice(STANDARD_NUCLEOTIDES))
    return ''.join(result)


def calculate_stats(sequence: str) -> dict:
    seq = ''.join(c for c in sequence if c.upper() in 'ACGT')
    n = len(seq)
    stats = {}
    for nuc in STANDARD_NUCLEOTIDES:
        stats[nuc] = (seq.upper().count(nuc) / n) * 100
    stats['GC'] = stats['G'] + stats['C']
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


def get_weighted_distribution() -> dict:
    """Prompts the user to enter percentage for each nucleotide.
    Validates that all four values sum to exactly 100.
    Returns a dict with keys A, C, G, T.
    """
    while True:
        print("Enter nucleotide percentages (must sum to 100):")
        try:
            a = float(input("  A (%): "))
            c = float(input("  C (%): "))
            g = float(input("  G (%): "))
            t = float(input("  T (%): "))
            if abs(a + c + g + t - 100) < 0.01:
                return {'A': a, 'C': c, 'G': g, 'T': t}
            print("Error: percentages must sum to 100.")
        except ValueError:
            print("Error: enter numeric values.")


def find_motif(sequence: str, motif: str) -> list:
    """Searches for all occurrences of a motif in the sequence.
    Returns a list of 1-based positions (biological convention).
    """
    positions = []
    seq = sequence.upper()
    motif = motif.upper()
    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)
        start = pos + 1
    return positions


def complement(sequence: str) -> str:
    """Returns the complementary DNA strand (5' to 3' direction preserved).
    Also handles IUPAC ambiguous characters by complementing each symbol.
    """
    iupac_complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    }
    return ''.join(iupac_complement.get(c, c) for c in sequence)


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complementary DNA strand."""
    return complement(sequence)[::-1]


def transcribe(sequence: str) -> str:
    """Returns the mRNA sequence by replacing T with U (in silico transcription)."""
    return sequence.upper().replace('T', 'U')


def translate(sequence: str) -> str:
    """Translates a DNA sequence into an amino acid sequence using the standard codon table.
    Stops translation at the first stop codon (*).
    Reference: https://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    seq = sequence.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = codon_table.get(codon, '?')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def sliding_window_gc(sequence: str, window: int, step: int) -> list:
    """Calculates GC content in a sliding window across the sequence.
    Returns a list of dicts with keys: start_position, gc_content.
    """
    results = []
    seq = sequence.upper()
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i + window]
        gc = (w.count('G') + w.count('C')) / len(w) * 100
        results.append({'start_position': i + 1, 'gc_content': round(gc, 2)})
    return results


def save_csv(data: list, filename: str):
    """Saves sliding window GC data to a CSV file with headers start_position and gc_content."""
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['start_position', 'gc_content'])
        writer.writeheader()
        writer.writerows(data)


def plot_gc_content(data: list, filename: str):
    """Generates a line chart of GC content along the sequence from sliding window data.
    Saves the chart as a PNG file using matplotlib.
    """
    positions = [d['start_position'] for d in data]
    gc_values = [d['gc_content'] for d in data]
    plt.figure(figsize=(12, 5))
    plt.plot(positions, gc_values, color='steelblue', linewidth=1)
    plt.title('GC Content Along the Sequence (Sliding Window)')
    plt.xlabel('Position (nt)')
    plt.ylabel('GC Content (%)')
    plt.ylim(0, 100)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def find_orfs(sequence: str, min_length: int = 100) -> list:
    """Finds all open reading frames (ORFs) in all 3 reading frames.
    An ORF starts at ATG and ends at the nearest stop codon (TAA/TAG/TGA).
    Only ORFs with length >= min_length nucleotides are returned.
    Returns a list of dicts with keys: start, end, length, frame (all 1-based positions).
    """
    orfs = []
    seq = sequence.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i + 3]
            if codon == 'ATG':
                found_stop = False
                for j in range(i + 3, len(seq) - 2, 3):
                    stop = seq[j:j + 3]
                    if stop in stop_codons:
                        length = j + 3 - i
                        if length >= min_length:
                            orfs.append({
                                'start': i + 1,
                                'end': j + 3,
                                'length': length,
                                'frame': frame + 1
                            })
                        i = j + 3
                        found_stop = True
                        break
                if not found_stop:
                    i += 3
            else:
                i += 3
    return orfs


def batch_mode():
    count = validate_positive_int("How many sequences to generate? ", 1, 1000)
    length = validate_positive_int("Enter sequence length for each sequence: ")
    description = input("Enter a description (optional, same for all): ")
    name = input("Enter your name: ")
    filename = "multi_sequences.fasta"
    with open(filename, 'w') as f:
        for i in range(1, count + 1):
            seq_id = f"Seq_{i:03d}"
            seq = generate_sequence(length)
            seq_with_name = insert_name(seq, name)
            f.write(format_fasta_record(seq_id, description, seq_with_name))
    print(f"Saved {count} sequences to {filename}")


def validate_fasta_file(filepath: str):
    """Loads a user-supplied FASTA file and validates its format.
    Checks header presence, allowed characters, and reports all errors found.
    """
    allowed = set('ACGTUacgtu' + ''.join(IUPAC_MAP.keys()) + ''.join(k.lower() for k in IUPAC_MAP.keys()))
    errors = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: file '{filepath}' not found.")
        return

    if not lines or not lines[0].startswith('>'):
        errors.append("Line 1: Missing header line starting with '>'.")

    in_sequence = False
    line_widths = []
    for i, line in enumerate(lines, 1):
        line = line.rstrip('\n')
        if line.startswith('>'):
            in_sequence = True
            line_widths = []
            continue
        if line.startswith('#') or line == '':
            continue
        if in_sequence:
            invalid = [c for c in line if c not in allowed]
            if invalid:
                errors.append(f"Line {i}: Invalid characters found: {set(invalid)}")
            line_widths.append(len(line))

    if errors:
        print("Validation errors found:")
        for e in errors:
            print(f"  - {e}")
    else:
        print("FASTA file is valid. No errors found.")


def main():
    """Main function. Handles user interaction, calls all feature functions, and writes output.
    Supports single sequence mode and batch mode with all 11 extended bioinformatics features.
    """
    print("=== DNA FASTA Generator ===\n")
    print("1) Single sequence")
    print("2) Batch mode")
    mode = input("Choose mode (1/2): ").strip()

    if mode == '2':
        batch_mode()
        return

    print("\nSequence generation method:")
    print("  1) Standard random (equal distribution)")
    print("  2) Custom nucleotide distribution")
    print("  3) IUPAC mode (with ambiguous characters)")
    gen_mode = input("Choose (1/2/3): ").strip()

    weights = None
    iupac_prob = 0.0

    if gen_mode == '2':
        weights = get_weighted_distribution()
    elif gen_mode == '3':
        while True:
            try:
                iupac_prob = float(input("Probability of ambiguous character per position (0.0 - 1.0): "))
                if 0.0 <= iupac_prob <= 1.0:
                    break
                print("Error: value must be between 0.0 and 1.0.")
            except ValueError:
                print("Error: enter a numeric value.")

    length = validate_positive_int("Enter sequence length: ")
    seq_id = validate_id("Enter sequence ID: ")
    description = input("Enter a description of the sequence: ")
    name = input("Enter your name: ")

    if gen_mode == '2' and weights:
        sequence = generate_sequence_weighted(length, weights)
    elif gen_mode == '3':
        sequence = generate_sequence_iupac(length, iupac_prob)
    else:
        sequence = generate_sequence(length)

    stats = calculate_stats(sequence)

    do_motif = input("Search for a motif? (y/n): ").strip().lower() == 'y'
    if do_motif:
        motif = input("Enter motif: ").strip().upper()
        positions = find_motif(sequence, motif)
        if positions:
            print(f"Motif '{motif}' found at positions: {positions}")
        else:
            print(f"Motif '{motif}' not found.")

    do_complement = input("Generate complement & reverse complement? (y/n): ").strip().lower() == 'y'
    do_transcribe = input("Generate mRNA transcription? (y/n): ").strip().lower() == 'y'
    do_translate = input("Translate to protein? (y/n): ").strip().lower() == 'y'
    do_orfs = input("Find ORFs? (y/n): ").strip().lower() == 'y'
    do_sliding = input("Sliding window GC analysis + chart? (y/n): ").strip().lower() == 'y'
    do_validate = input("Validate an existing FASTA file? (y/n): ").strip().lower() == 'y'

    seq_with_name = insert_name(sequence, name)

    records = [format_fasta_record(seq_id, description, seq_with_name)]

    if do_complement:
        comp_seq = complement(sequence)
        rev_comp_seq = reverse_complement(sequence)
        records.append(format_fasta_record(seq_id + "_comp", "complementary strand", comp_seq))
        records.append(format_fasta_record(seq_id + "_revcomp", "reverse complementary strand", rev_comp_seq))

    if do_transcribe:
        mrna = transcribe(sequence)
        records.append(format_fasta_record(seq_id + "_mRNA", "mRNA transcription (T->U)", mrna))

    if do_translate:
        protein = translate(sequence)
        records.append(format_fasta_record(seq_id + "_protein", "translated protein sequence", protein))

    filename = f"{seq_id}.fasta"
    with open(filename, 'w') as f:
        for record in records:
            f.write(record)

    print(f"\nSequence saved to file: {filename}")

    print(f"\nSequence statistics (n={length}):")
    for nuc in STANDARD_NUCLEOTIDES:
        print(f"  {nuc}: {stats[nuc]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")

    if do_orfs:
        min_orf = validate_positive_int("Minimum ORF length (nt): ", 1, 100_000)
        orfs = find_orfs(sequence, min_orf)
        if orfs:
            print(f"\nFound {len(orfs)} ORF(s):")
            for orf in orfs:
                print(f"  Frame {orf['frame']}: start={orf['start']}, end={orf['end']}, length={orf['length']} nt")
        else:
            print("No ORFs found.")

    if do_sliding:
        window = validate_positive_int("Window size (nt): ", 1, length)
        step = validate_positive_int("Step size (nt): ", 1, length)
        gc_data = sliding_window_gc(sequence, window, step)
        csv_file = f"{seq_id}_gc_sliding.csv"
        save_csv(gc_data, csv_file)
        print(f"Sliding window GC data saved to {csv_file}")
        chart_file = f"{seq_id}_gc_chart.png"
        plot_gc_content(gc_data, chart_file)
        print(f"GC content chart saved to {chart_file}")

    if do_validate:
        fasta_path = input("Enter path to FASTA file to validate: ").strip()
        validate_fasta_file(fasta_path)


if __name__ == "__main__":
    main()
