'''
 # @ Author: Charles Brottier and Jacky Ame
 # @ Creation: 2021-11-21 
 # @ Last update: 2022-01-06
 # @ Python Version : 3
'''

from timer import Timer
from Burrows_wheeler_minimal import *
from tools_karkkainen_sanders import simple_kark_sort
import argparse

###################################

## Argparse & help message setup ##

###################################


def help_msg():
    """
    Customized help message in case of user failing to input a valid command.

    Returns:
        [help]: [description]
    """
    return '''
        To use this script, please respect the following format:

        SequencesToSPSS.py -i reads.fasta -g genome.fasta [-t] x [-k] y
        
        [-i (required): Input the fasta file containing your reads.]
        [-g (required): Input the fasta file containing your reference genome.]
        [-t (optional [2]): Input x = the solidity threshold you want to use.]
        [-k (optional [31]): Input y = the k-mer size you want to use.]
        
    '''


parser = argparse.ArgumentParser(
    add_help=False,
    usage=help_msg(),
    description="Compressed representation and indexing\
        of a set of k-mers extracted from sequencing data.")
parser.add_argument("-i",
                    "--reads",
                    metavar='',
                    help="fasta file containing a set of reads",
                    dest='reads',
                    type=str,
                    required=True)
parser.add_argument("-g",
                    "--genome",
                    metavar='',
                    help="fasta file containing a reference genome",
                    dest='genome',
                    type=str,
                    required=True)
parser.add_argument("-t",
                    "--solidity_threshold",
                    metavar='',
                    help="solidity threshold (kmers occurring less than t\
                        times are not extracted) (default=2)",
                    dest='solid_t',
                    type=int,
                    required=False,
                    default=2)
parser.add_argument("-k",
                    "--kmer_size",
                    metavar='',
                    help="kmer size (default=31)",
                    dest='k_size',
                    type=int,
                    required=False,
                    default=31)

args = parser.parse_args()

###################################

######## Utility Functions ########

###################################


def validate_dna_seq(seq: str) -> bool:
    """
    Checks DNA sequences for incorrect characters.
    Anything else than A,T,G and C is considered incorrect.

    Args:
        seq (str): a DNA read

    Raises:
        Exception: sequence is empty
        TypeError: sequence is not a string
        Exception: incorrect character detected

    Returns:
        bool: sequence validity
    """
    if seq == "":
        raise Exception("Sequence is empty.")
    if not isinstance(seq, str):
        raise TypeError("Sequence is not readable (must be a string).")

    for nuc in seq:  # check for unwanted character
        if nuc.upper() not in 'ATGC':
            raise Exception("Incorrect character detected.")
    return True


def rev_comp(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        seq (str): a DNA sequence

    Returns:
        str: the reverse complement of the DNA sequence
    """
    rev = str.maketrans("ACGT", "TGCA")
    return seq.translate(rev)[::-1]


def canonical(kmer: str) -> str:
    """
    Returns the canonical form of a k-mer.

    Args:
        kmer (str): k-mer

    Returns:
        str: canonical k-mer
    """
    return min(kmer, rev_comp(kmer))


def kmer_generation(seq: str, kmer_freq_dict: dict, k_size: int) -> dict:
    """
    Fills a dictionary with the canonical k-mers present in a DNA sequence,
    associated with their frequencies.

    Args:
        seq (str): source DNA sequence
        kmer_freq_dict (dict): target dictionary
        k_size (int): size of the k-mers

    Returns:
        kmer_freq_dict (dict): canonical k-mers and their frequencies
    """
    for i in range(0, len(seq) - k_size + 1):
        # selecting the canonical k-mers
        kmer = canonical(seq[i:i + k_size])

        # while counting them, adding them to the dictionary
        if kmer in kmer_freq_dict:
            kmer_freq_dict[kmer] += 1
        else:
            kmer_freq_dict[kmer] = 1

    return kmer_freq_dict


# Loading Functions


def fasta_to_kmer(fasta_file: str, k_size: int, solid_t: int) -> set:
    """
    Reads the DNA sequences in a fasta file and returns a set of all its
    canonical kmers (of size k), as well as their frequencies, provided
    they pass the solidity threshold (t).

    Args:
        fasta_file (str): path to a fasta file
        k_size (int): size of the kmers
        solid_t (int): solidity threshold

    Raises:
        Exception: when encountering an invalid DNA sequence

    Returns:
        solid_kmers (set): set of all the solid canonical kmers
    """

    kmer_freq_dict = {}
    with open(fasta_file, "r") as file:
        for line in file:
            if not line.startswith(">"):
                line = line.strip().upper()
                # each DNA sequence must be valid
                if validate_dna_seq(line):
                    # filling the dictionary with k-mers and their frequencies
                    kmer_freq_dict = (kmer_generation(line, kmer_freq_dict,
                                                      k_size))
                else:
                    raise Exception(
                        f"A sequence of {fasta_file} is not valid.")

    # selecting solid k-mers
    solid_kmers = set()
    for kmer, freq in list(kmer_freq_dict.items()):
        if freq >= solid_t:
            solid_kmers.add(kmer)

    return solid_kmers


def fasta_simple_extract(fasta_file: str) -> str:
    """
    Extracts the DNA reads from a fasta file and returns a single sequence.

    Args:
        fasta_file (str): path to a fasta file

    Raises:
        Exception: when encountering an invalid DNA sequence.

    Returns:
        str: concatenation of all DNA reads
    """
    seq = ""
    with open(fasta_file, "r") as file:
        for line in file:
            if not line.startswith(">"):
                line = line.strip().upper()
                if validate_dna_seq(line):
                    seq += line
                else:
                    raise Exception(
                        f"A sequence in {fasta_file} is not valid.")

    return seq


###################################

##### SPSS Building Functions #####

###################################


def kmer_forward_extension(extension: str, kmer_set: set, k_size: int):
    """
    Returns the forward extension (to the right) of a given k-mer.

    Args:
        start (str): starting k-mer to extend
        kmer_set (set): set of all canonical k-mers
        k_size (int): size of the k-mers

    Returns:
        extension (str): forward extension of the k-mer
        kmer_set (set): set of all canonical k-mers
    """
    can_extend = True

    # loop set up to stop when no k-mer allowing extension remains
    while can_extend:
        can_extend = False
        # setting the suffix for extension
        suffix = extension[-(k_size - 1):]
        for letter in 'ATGC':
            next_k = canonical(suffix + letter)

            # extending if possible and looping over
            if next_k in kmer_set:
                extension = extension + letter
                can_extend = True
                # removing the used kmer from the list
                kmer_set.remove(next_k)
                break

    return extension, kmer_set


def unitig_generation(start: str, kmer_set: set, k_size: int):
    """
    Returns the full extension to the right and to the left of a given k-mer,
    which is a maximal unitig, as well as the unused k-mers.

    Args:
        start (str): k-mer to extend
        kmer_set (set): set of all canonical k-mers
        k_size (int): size of the k-mers

    Returns:
        extension (str): full extension of the k-mer
        kmer_set (set): remaining canonical k-mers
    """
    # forward extension
    extension, kmer_set = kmer_forward_extension(start, kmer_set, k_size)
    # reverse extension
    extension = rev_comp(extension)
    extension, kmer_set = kmer_forward_extension(extension, kmer_set, k_size)

    return extension, kmer_set


def all_unitigs(kmer_set: set, k_size: int) -> set:
    """
    Greedily generates all maximal unitigs from a set of k-mers.

    Args:
        kmer_set (set): canonical k-mer set
        k_size (int): size of the k-mers

    Returns:
        unitigs (set): set of all maximal unitigs
    """
    unitigs = set()
    # generating unitigs as lons as there are k-mers remaining in the set
    while kmer_set != set():
        start = kmer_set.pop()
        extension, kmer_set = unitig_generation(start, kmer_set, k_size)
        unitigs.add(extension)

    return unitigs


def unitig_compression(unitigs: set, ref_dict: dict, anchor_len: int) -> set:
    """
    Compresses unitigs based on a suffix/prefix match.
    This means unitig X and Y will be merged into Z if they share an "anchor"
    such as: X = ATGCCC[00000] ; Y = [00000]ATACG ; Z= ATGCCC[00000]ATACG ;
    the anchor being [00000], of length k = 5, in this case.

    Args:
        unitigs (set): set of unitigs
        ref_dict (dict): dictionary to keep track of used/unused unitigs
        anchor_len (int): length of the anchor allowing merging.

    Returns:
        compressed (set): all merged unitigs with the specified anchor length
    """

    compressed = set()
    for u_1 in unitigs:
        for u_2 in unitigs:
            # u_1 and u_2 must be different and unused
            if u_1 != u_2 and ref_dict[u_1] == 0 and ref_dict[u_2] == 0:
                # merging in case of suffix/prefix match
                if u_2[-anchor_len:] == u_1[:anchor_len]:
                    u_c = u_2 + u_1[anchor_len:]  # u_1 and u_2 merged
                    compressed.add(u_c)
                    ref_dict[u_1] += 1  # u_1 and u_2 marked as used
                    ref_dict[u_2] += 1
    return compressed


def spss_generation(unitigs: set, k_size: int) -> str:
    """
    Maximally compresses all unitigs in a step-by-step approach that consists in
    merging all unitigs based on a decreasing anchor length at each step.
    The compressed unitigs are then concatenated to form the SPSS.
    The maximal anchor length has empirically proven to be|k-mer| - 2.

    Args:
        unitigs (set): set of unitigs
        k_size (int): size of the k-mers

    Returns:
        (str): maximally compressed unitigs merged into the SPSS
    """
    compressed = set()
    anchor_len = k_size - 2  # maximal sanchor length is |k-mer| - 2

    while anchor_len > 0:
        ref = dict()  # dictionary used to flag merged unitigs
        for seq in unitigs:
            ref[seq] = 0  # value of 0 means unused (not yet merged)

        # compressing
        compression = unitig_compression(unitigs, ref, anchor_len)
        for seq in compression:
            compressed.add(seq)

        # dealing with unused unitigs
        for key, val in ref.items():
            # unused unitigs must be kept
            if val == 0:
                compression.add(key)
            # previously unused unitigs that become used must be removed
            if val == 1 and key in compressed:
                compressed.remove(key)

        unitigs = compression
        anchor_len -= 1  # decreasing anchor length for next iteration

    return ''.join(compressed)


def spss_max_compression(spss: str) -> str:
    """
    Maximally compresses the SPSS by replacing all series of identical
    characters with their number of repetitions (e.g. 'AAAA' becomes 'A4').

    Args:
        spss (str): SPSS to compress

    Returns:
        str: maximally compressed string
    """
    idx = 0
    min_spss = ""

    # screening the sequence until the before last caracter
    while idx != len(spss):
        count = 1
        # detecting repeats
        while (idx < len(spss) - 1) and (spss[idx] == spss[idx + 1]):
            count = count + 1
            idx = idx + 1
        if count == 1:
            min_spss += str(spss[idx])
        # replacing character repeats by "character + its occurrences count"
        else:
            min_spss += str(spss[idx]) + str(count)
        idx = idx + 1

    return min_spss


"""
NOT IMPLEMENTED. Raw attempt to reverse the compression, but takes too long.
Works fine on small sequences (AAAAATTTTTGGGGGCCCCC), but decompressed sequence
length doesn't match the original one with an actual SPSS. Out of time to debug.

def spss_decompression(spss: str) -> str:
    '''
    Decompression of the maximally compressed SPSS.
    '''
    idx = 0
    restored = ""

    while idx < len(spss):
        for idx, c in enumerate(spss):
            if c.isnumeric():
                for _ in range(int(c)):
                    restored += spss[idx - 1]
            idx += idx

    return restored
"""

###################################

#### SPSS: Querying & Indexing ####

###################################


def naive_fn_count(genome_kmers: set, spss: str) -> int:
    """
    Returns the number of false negatives (all k-mers in the genome 
    but not in SPSS sequence).

    Args:
        genome_kmers (set): k-mers of the genome
        spss_kmers (set): k-mers of the SPSS sequence

    Returns:
        fn_count (int): number of false negatives
    """
    fn_count = 0
    for kmer in genome_kmers:
        if kmer not in spss and rev_comp(kmer) not in spss:
            fn_count += 1

    return fn_count


def naive_fp_count(spss_kmers: set, genome: str) -> int:
    """
    Counts the number of false positives (k-mers in the SPSS but not in the 
    genome).

    Args:
        genome_kmers (set): k-mers in the genome
        spss_kmers (set): k-mers in the SPSS

    Returns:
        fp_count (int): number of false positives
    """
    fp_count = 0
    for kmer in spss_kmers:
        if kmer not in genome and rev_comp(kmer) not in genome:
            fp_count += 1

    return fp_count


def bwt_indexation(string: str) -> Tuple[str, dict, list]:
    """
    Applies the Burrows-Wheeler Transform on the passed string and returns the
    BWT, the dictionary of characters occurrences, and the characters ranks.

    Args:
        string (str): [description]

    Returns:
        [str, dict, list]: BWT ; occurrences and ranks of the characters (resp.)
    """
    string += '$'

    suffix_ar = simple_kark_sort(string)
    bwt = get_bwt(string, suffix_ar)
    # assert len(bwt) == len(string)
    n_dict = get_n(bwt)
    ranks = get_ranks(bwt)

    return bwt, n_dict, ranks


def indexed_fn_count(spss: str, genome_kmers: set) -> int:
    """
    Counts the number of false negatives (k-mers in the genome but not in the 
    SPSS) using the BWT.

    Args:
        spss (str): the SPSS
        genome_kmers (set): k-mers in the genome

    Returns:
        fn_count (int): number of false negatives
    """
    # BWT on the SPSS
    bwt, n_dict, ranks = bwt_indexation(spss)

    # Counting
    fn_count = 0
    for kmer in genome_kmers:
        present = contains(kmer, bwt, n_dict, ranks)[0]
        if not present:
            present = contains(rev_comp(kmer), bwt, n_dict, ranks)[0]
            if not present:
                fn_count += 1

    return fn_count


def indexed_fp_count(genome: str, spss_kmers: set) -> int:
    """
    Counts the number of false positives (all k-mers in the SPSS but not in 
    the genome) using the BWT.

    Args:
        genome (str): complete reference genome + '$'
        spss_kmers (set): All k-mers in the SPSS

    Returns:
        fp_count (int): number of false positives
    """

    # BWT on the genome
    bwt, n_dict, ranks = bwt_indexation(genome)

    # Counting
    fp_count = 0
    for kmer in spss_kmers:
        present = contains(kmer, bwt, n_dict, ranks)[0]
        if not present:
            present = contains(rev_comp(kmer), bwt, n_dict, ranks)[0]
            if not present:
                fp_count += 1
    return fp_count


# ############################################################################ #
#                                     MAIN                                     #
# ############################################################################ #


def main(reads, genome, k_size, solid_t):
    """
    Sine qua non conditions 
    """
    if k_size < 1:
        raise ValueError("K-mer size must be positive.")
    if solid_t < 1:
        raise ValueError("Minimal solidity threshold is 1.")
    # ------------------------------------------------------------------------ #
    """
    Building the k-mers set
    """
    with Timer() as total_time:
        spss_kmers = fasta_to_kmer(reads, k_size, solid_t)
    total_time.print("OUT TIME_SELECTING_KMERS={}")
    print(f"Nb of SOLID k-mers={len(spss_kmers)}")
    # ------------------------------------------------------------------------ #
    """
    Building the SPSS
    """
    with Timer() as total_time:
        unitigs = all_unitigs(spss_kmers, k_size)
        spss = spss_generation(unitigs, k_size)
    total_time.print("OUT TIME_SPSS_CONSTRUCTION={}")
    print(f"OUT |SPSS(K)|={len(spss)}")
    print(f"OUT #SPSS(K)={len(unitigs)}")
    print(f"Compression={round(len(spss)/len(''.join(unitigs))*100,2)}%")
    # ------------------------------------------------------------------------ #
    """
    Checking for false negatives and false positives
    """
    ref_genome = fasta_simple_extract(genome)
    genome_kmers = {}
    genome_kmers = kmer_generation(ref_genome, genome_kmers, k_size)
    # ------------------------------------------------------------------------ #
    """
    False negatives count
    """

    # Naive

    with Timer() as total_time:
        fn_count_naive = naive_fn_count(genome_kmers, spss)
    total_time.print("OUT TIME_FN_WITHOUT_INDEX={}")

    # Indexed

    with Timer() as total_time:
        fn_count_idx = indexed_fn_count(spss, genome_kmers)
    total_time.print("OUT TIME_FN_WITH_INDEX={}")
    assert fn_count_naive == fn_count_idx
    print(f"#FN={fn_count_idx}")
    # ------------------------------------------------------------------------ #
    """
    False positives count
    """
    # Naive

    spss_kmers = fasta_to_kmer('fasta_files/ecoli_sample_reads.fasta', k_size,
                               solid_t)
    with Timer() as total_time:
        fp_count_naive = naive_fp_count(spss_kmers, ref_genome)  # 55 sec
    total_time.print("OUT TIME_FP_WITHOUT_INDEX={}")

    # Indexed

    with Timer() as total_time:
        fp_count_idx = indexed_fp_count(ref_genome, spss_kmers)  # 6 sec
    total_time.print("OUT TIME_FP_WITH_INDEX={}")
    assert fp_count_naive == fp_count_idx
    print(f"#FP={fp_count_idx}")
    # ------------------------------------------------------------------------ #
    """
    Maximum SPSS Compression
    """
    with Timer() as total_time:
        min_spss = spss_max_compression(spss)
    total_time.print("OUT TIME_SPSS_MAX_COMPRESSION={}")
    print(f"OUT |min SPSS(K)|={len(min_spss)}")
    print(
        f"Max compression={round(len(min_spss)/len(''.join(unitigs))*100,2)}%")


###################################

#########  Main execution  ########

###################################

if __name__ == "__main__":
    with Timer() as total_time:
        main(args.reads, args.genome, args.k_size, args.solid_t)
    total_time.print("Total run time = {} seconds")
