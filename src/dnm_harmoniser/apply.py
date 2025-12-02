import dnm_harmoniser

def get_homopolymer_length(ref: str, alt: str) -> int:
    """
    Calculate the length of the longest homopolymer (repeated nucleotide) in REF or ALT.

    For deletions: checks the deleted sequence (REF - ALT)
    For insertions: checks the inserted sequence (ALT - REF)
    For SNVs: returns 1 (single nucleotide)

    Parameters
    ----------
    ref : str
        Reference allele sequence
    alt : str
        Alternate allele sequence

    Returns
    -------
    int
        Length of longest homopolymer run
    """
    if not isinstance(ref, str) or not isinstance(alt, str):
        return 0

    ref = ref.upper()
    alt = alt.upper()

    # For SNVs, homopolymer length is 1
    if len(ref) == 1 and len(alt) == 1:
        return 1

    # For indels, check the longer sequence (contains the repetitive region)
    seq_to_check = ref if len(ref) > len(alt) else alt

    if len(seq_to_check) == 0:
        return 0

    # Find longest run of repeated nucleotides
    max_length = 1
    current_length = 1
    current_base = seq_to_check[0]

    for base in seq_to_check[1:]:
        if base == current_base:
            current_length += 1
            max_length = max(max_length, current_length)
        else:
            current_base = base
            current_length = 1

    return max_length

data = dnm_harmoniser.VariantDataset.from_tsv("../all_trios_N1316_denovoCNN0.6.csv")
reference = dnm_harmoniser.VariantDataset.from_tsv("../decode_parages.txt")

filtered_variants = data.variants.copy()
filtered_ref = reference.variants.copy()

filtered_variants['length'] = filtered_variants['ALT'].str.len() - filtered_variants['REF'].str.len()
filtered_ref['length'] = filtered_ref['Alt'].str.len() - filtered_ref['Ref'].str.len()

filtered_variants = filtered_variants[filtered_variants['length'].abs() < 20]
filtered_ref = filtered_ref[filtered_ref['length'].abs() < 20]

filtered_variants['homopolymer_length'] = filtered_variants.apply(lambda row: get_homopolymer_length(row['REF'], row['ALT']),axis=1)
filtered_ref['homopolymer_length'] = filtered_ref.apply(lambda row: get_homopolymer_length(row['Ref'], row['Alt']),axis=1)

filtered_variants = filtered_variants[filtered_variants['homopolymer_length'] < 8]
filtered_ref = filtered_ref[filtered_ref['homopolymer_length'] < 8]

filtered_variants.to_csv("input_filtered.tsv", index=False, sep="\t")
filtered_ref.to_csv("reference_filtered.tsv", index=False, sep="\t")