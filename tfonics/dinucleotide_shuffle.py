"""Based on altschulEriksonDinuclShuffle.py. Originally by P. Clote, Oct 2003.

http://clavius.bc.edu/clotelab/RNAdinucleotideShuffle/ShuffleCodeParts/altschulEriksonDinuclShuffle.txt

This is an implementation of the Altschul and Erickson dinucleotide shuffle, as described in
https://doi.org/10.1093/oxfordjournals.molbev.a040370

Reworked and refactored to conform with pep8, as well as some performance tweaks
"""

import random
from collections import defaultdict

# This can be reworked with a sequence object, which would be handy
# when generating multiple strings from the same sequence, and would clean up some
# passing around of parameters. We may have to deep copy the dinucleotide_sequence
# parameter before shuffling, depending on whether we care about preserving the
# original sequence.


def get_dinucleotide_sequence(section):
    """Count the number of occurences of single and pairs of nucleotides in
    a section of DNA.

    Input parameters:
    section - string of nucleoutides

    Returns:
    dinucleotide_sequence - what nucleotides Y followed immediately after nucleotide X?
                            i.e dinucleotide_sequence[X] = [Y1, Y2, Y3] if the pairs
                            X Y1, X Y2 and X Y3 occured in section
    """
    # Construct a dictionary of lists to represent the vertices in the graph
    dinucleotide_sequence = defaultdict(list)

    # Compute count and lists
    for i in range(len(section)-1):
        start = section[i]
        end = section[i+1]
        dinucleotide_sequence[start].append(end)

    return dinucleotide_sequence


def connected_to_last(edge_list, nucleotide_list, last_character):
    """Make sure that an entry in the edge list is connected to the last character
    of the section.

    Inputs:
    edge_list - list of edges chosen
    nucleotide_list - list of nucleotides that occur in the section
    last_character - last character in the section

    Returns:
    boolean: whether or not the edge list forms a fully connected graph, and that
    graph is connected to the last character.
    """
    # Dictionary of whether a nucleotide is connected to the last character
    connected = {}
    for nucelotide in nucleotide_list:
        connected[nucelotide] = False

    # The last character is trivially connected to itself
    connected[last_character] = True

    # Back-propagate from the last nucleotide to all others
    # The longest possible path in the graph is n-1
    for _ in range(len(nucleotide_list)-1):
        for (start, end) in edge_list:
            if connected[end]:
                connected[start] = True

    return all(connected.values())


def pick_edges(section, nucleotide_list, dinucleotide_sequence, rng):
    """ Randomly pick a set of edges from the graph of dinucleotide connections,
    ensuring that the edges are all connected directly or indirectly to the last
    character in the section.

    Parameters:
    section - string of nucleotides to be shuffled
    nucleotide_list - list of unique nucleotides appearing in section
    dinucleotide_sequence - dictionary of lists of next-neigbhour nucleotide pairs

    Returns:
    edge_list - list of chosen graph edges
    """
    last_character = section[-1]

    while True:
        edge_list = []
        for start in nucleotide_list:
            if start != last_character:
                end = rng.choice(dinucleotide_sequence[start])
                edge_list.append((start, end))

        if connected_to_last(edge_list, nucleotide_list, last_character):
            break
    return edge_list


def get_nucleotide_list(section):
    """Find all unique characters in section, and make sure it
    only contains valid nucleotides.

    Parameters:
    section - string of nucleotides

    Returns:
    list of nucleotides found in section

    Throws:
    AssertionError if string contains letters other than A,C,G and T
    """
    nucleotide_set = set(section)
    invalid_nucleotides = nucleotide_set - {"A", "C", "G", "T"}
    assert not invalid_nucleotides, (
        "Input string contained non-nucleotide letters: {}"
        .format(invalid_nucleotides)
    )

    # The order of the nucleotide set is non-deterministic, so sort it
    # to ensure reproducibility
    return sorted(list(nucleotide_set))


def dinucleotide_shuffle(section, rng=None):
    """ Construct a new random eulerian path throught he nucleotide string by
    - Counting all the edges in the nucleotide connection graph
    - Randomly picking an edge emerging from each vertex in the graph except for
      one corresponding to the nucleotide in the section.
    - Ensure that the choices above are all connected directly or indirectly
      to the final character in the section, otherwise choose again
    - Remove these edges from the graph
    - Shuffle the remaining edges
    - Reinsert the removed edges
    - Traverse the graph to obtain a new, shuffled string

    Inputs:
    - section - string of nucleotides to be shuffled
    - rng - optional: instance of random.Random to be used for shuffling
    Outputs:
    shuffled string of nucleotides
    """

    if not rng:
        rng = random.Random()

    # Convert to uppercase before continuing
    section = section.upper()
    nucleotide_list = get_nucleotide_list(section)

    # dinucleotide sequence: { "A": ["A", "C", "A", "T",...], "T": [...]}
    dinucleotide_sequence = get_dinucleotide_sequence(section)
    edge_list = pick_edges(section, nucleotide_list, dinucleotide_sequence, rng)

    # move the edges in edge_list to the end of the vertex list, shuffle all other edges.
    for (start, end) in edge_list:
        dinucleotide_sequence[start].remove(end)
    for nucleotide in nucleotide_list:
        rng.shuffle(dinucleotide_sequence[nucleotide])
    for (start, end) in edge_list:
        dinucleotide_sequence[start].append(end)

    # construct the eulerian path
    path = section[0]
    previous_character = section[0]

    for _ in range(len(section)-1):
        current_character = dinucleotide_sequence[previous_character][0]
        path += current_character
        del dinucleotide_sequence[previous_character][0]
        previous_character = current_character

    return path
