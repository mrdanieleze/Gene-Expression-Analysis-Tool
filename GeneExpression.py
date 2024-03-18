from Bio.Seq import Seq

def protein_to_rna(protein_sequence, codon_table):
    rna_sequence = ''
    for amino_acid in protein_sequence:
        if amino_acid not in codon_table:
            raise ValueError(f"Invalid amino acid '{amino_acid}' in the protein sequence.")
        possible_codons = codon_table[amino_acid]
        # Choose the codon based on its frequency
        chosen_codon = max(possible_codons, key=lambda x: x[1])[0].replace('T', 'U')
        rna_sequence += chosen_codon
    return rna_sequence

codon_table = {
    'A': [('GCU', 19.7), ('GCC', 15.0), ('GCA', 15.2), ('GCG', 11.9)],
    'C': [('UGU', 16.8), ('UGC', 14.6)],
    'D': [('GAU', 21.9), ('GAC', 24.4)],
    'E': [('GAA', 33.2), ('GAG', 18.4)],
    'F': [('UUU', 30.5), ('UUC', 18.2)],
    'G': [('GGU', 21.3), ('GGC', 33.4), ('GGA', 9.2), ('GGG', 8.6)],
    'H': [('CAU', 15.8), ('CAC', 13.1)],
    'I': [('AUU', 30.5), ('AUC', 18.2), ('AUA', 3.7)],
    'K': [('AAA', 33.2), ('AAG', 12.1)],
    'L': [('UUA', 15.2), ('UUG', 11.9), ('CUU', 11.9), ('CUC', 10.5), ('CUA', 5.3), ('CUG', 46.9)],
    'M': [('AUG', 24.8)],
    'N': [('AAU', 21.9), ('AAC', 24.4)],
    'P': [('CCU', 8.4), ('CCC', 6.4), ('CCA', 6.6), ('CCG', 26.7)],
    'Q': [('CAA', 12.1), ('CAG', 27.7)],
    'R': [('CGU', 21.1), ('CGC', 26.0), ('CGA', 4.3), ('CGG', 4.1), ('AGA', 1.4), ('AGG', 1.6)],
    'S': [('UCU', 5.7), ('UCC', 5.5), ('UCA', 7.8), ('UCG', 8.0), ('AGU', 7.2), ('AGC', 16.6)],
    'T': [('ACU', 8.0), ('ACC', 22.8), ('ACA', 6.4), ('ACG', 11.5)],
    'V': [('GUU', 16.8), ('GUC', 11.7), ('GUA', 11.5), ('GUG', 26.4)],
    'W': [('UGG', 10.7)],
    'Y': [('UAU', 16.8), ('UAC', 14.6)],
    '*': [('UAA', 1.8), ('UAG', 0.0), ('UGA', 1.0)]  # Stop codons
}

def GeneExpression(sequence):
    # Check if the sequence contains both T and U
    if 'T' in sequence and 'U' in sequence:
        print("Warning: Sequence contains both T and U. Sequence is invalid.")
        return
    
    # Check if the sequence contains only A, G, C, and U
    if all(base in 'AGCU' for base in sequence):
        # Transcribe RNA to DNA
        dna_sequence = Seq(sequence).back_transcribe()
        # Translate RNA to protein
        protein_sequence = Seq(sequence).translate()
        return dna_sequence, sequence, protein_sequence
    
    # Check if the sequence contains only A, T, G, and C
    if all(base in 'ATGC' for base in sequence):
        # Transcribe DNA to RNA
        rna_sequence = Seq(sequence).transcribe()
        # Translate RNA to protein
        protein_sequence = Seq(rna_sequence).translate()
        return sequence, rna_sequence, protein_sequence
    
    # Check if the sequence contains other letters besides A, T, G, and C
    if any(base not in 'ATGC' for base in sequence):
        # Back-translate protein sequence to RNA
        rna_sequence = protein_to_rna(sequence, codon_table)
        # Back-transcribe RNA to DNA
        dna_sequence = Seq(rna_sequence).back_transcribe()
        return dna_sequence, rna_sequence, sequence

# Example usage:
sequence = "IHEARTGIRL"  
result = GeneExpression(sequence)

if result:
    print("DNA sequence:")
    print(result[0])
    print("RNA sequence:")
    print(result[1])
    print("Protein sequence:")
    print(result[2])
