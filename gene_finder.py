# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement(9)
    'Invalid input'

    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'Invalid input'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complement = [get_complement(a) for a in dna] #makes a list of the complements
    reverse_complement = ''.join(complement[::-1]) #makes the list a string and reverses it
    return reverse_complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAAGTAGG")
    'ATGAAG'
    >>> rest_of_ORF("ATGAAGAAG")
    'ATGAAGAAG'
    >>> rest_of_ORF("ATGAGAATAGG")
    'ATGAGAATAGG'
    >>> rest_of_ORF("ATGCATGAATGTAGATAGATGTAATGCCC")
    'ATGCATGAATGTAGA'
    """
    stop_codons = ['TAG','TAA','TGA']
    triples = []
    array = []
    for i in range (0,len(dna),3):
        triples.append(dna[i:i+3])
    found_string = False
    for i in stop_codons:
        if i in triples:
            index = triples.index(i)
            array.append(index)
            found_string = True #if there is a stop codon sets variable to True

    if found_string:
        return dna[:min(array)*3] #returns the dna until the stop codon if a stop codon was found
    else:
        return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCCGTAAATGCGAA")
    ['ATGCATGAATGTAGA', 'ATGTGCCCG', 'ATGCGAA']
    """
    result = []
    while dna: #empty string of dna is False
        if dna[:3] == 'ATG': #checks if the first 3 makes a start codon
            frame = rest_of_ORF(dna)
            result.append(frame)
            dna = dna[len(frame):]
        else:
            dna = dna[3:]
    return result


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frame_1 = find_all_ORFs_oneframe(dna)
    dna = dna[1:]
    frame_2 = find_all_ORFs_oneframe(dna)
    dna = dna[1:]
    frame_3 = find_all_ORFs_oneframe(dna)
    result = frame_1+frame_2+frame_3
    return result

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_dna = get_reverse_complement(dna)
    frame_a = find_all_ORFs(dna)
    frame_b = find_all_ORFs(reverse_dna)
    result = frame_a + frame_b
    return result

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    frames = find_all_ORFs_both_strands(dna)
    longest = max(frames or ['']) #if there are no frames, gives a value of ''
    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    ORFs = []
    for trial in range(num_trials):
        shuffled_dna = shuffle_string(dna) #creates a variable for the shuffled dna
        longest_frame = len(longest_ORF(shuffled_dna)) #finds length of longest ORF of the shuffled dna
        ORFs.append(longest_frame)
    return max(ORFs) 

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    triples = []
    amino_list = []
    for i in range (0,len(dna),3):
        triples.append(dna[i:i+3]) #occupies the list called triples with dna split up in 3s
    for triplet in triples:
        if len(triplet) == 3:
            amino_acid = aa_table[triplet]
            amino_list.append(amino_acid)
    return ''.join(amino_list)

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    final_frame = []
    threshold_frames = []
    threshold = longest_ORF_noncoding(dna, 1500) #sets the threshold for length of dna
    all_frames = find_all_ORFs_both_strands(dna)
    for frame in all_frames:
        if len(frame) > threshold:
            threshold_frames.append(frame)#appends frames longer than the threshold to the list threshold_frames
    for i in range(0,len(threshold_frames)):
        final_frame.append(coding_strand_to_AA(threshold_frames[i]))#changes frames to amino acids and appends to final_frame
    return final_frame

if __name__ == "__main__":
    import doctest
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    #print dna
    print gene_finder(dna)
    #doctest.run_docstring_examples(longest_ORF,globals(), verbose=True)
