import re
from typing import List, Tuple

from graphtagger.seqOps import revComp


def count_continuous_runs(dna_string: str) -> list:
    """
    Count the length of all continuous runs of a DNA base in a string.

    Parameters:
    - dna_string (str): The input DNA string.

    Returns:
    - list: A list of tuples where each tuple contains a character and the number
            of times it occurred consecutively in the input DNA string.
    """

    # Check if the input string is empty
    if not dna_string:
        return []

    # Initialize variables
    current_base = dna_string[0]
    current_count = 1
    result = []

    # Iterate through the DNA string starting from the second character
    for base in dna_string[1:]:
        if base == current_base:
            # If the current base is the same as the previous one, increment the count
            current_count += 1
        else:
            # If the current base is different, append the tuple and reset the count
            result.append((current_base, current_count))
            current_base = base
            current_count = 1

    # Append the last tuple after the loop
    result.append((current_base, current_count))

    return result


def construct_regex_pattern(motif_tuples: List[Tuple[str, int]]) -> str:
    """
    Construct a regex pattern to match sequences with runs of characters
    that differ by plus or minus one compared to the original input sequence.

    Parameters:
    - motif_tuples (List[Tuple[str, int]]): List of tuples where each tuple contains
                                            a character and the number of times it
                                            occurred consecutively.

    Returns:
    - str: The constructed regex pattern as a raw string.
    """

    pattern_parts = []

    for char, count in motif_tuples:
        if count == 1:
            pattern_parts.append(
                re.escape(char)
            )  # If count is 1, just escape the character
        else:
            # If count is greater than 1, add a range allowing for plus or minus one
            pattern_parts.append(rf"{re.escape(char)}{{{count-1},{count+1}}}")

    return rf"{''.join(pattern_parts)}"


def find_repeats_of_motif(text: str, motif: str) -> List[Tuple[int, int]]:
    """
    Find repeats of a motif in a given text.

    Args:
        text (str): The input text to search for motifs.
        motif (str): The motif to search for in the text.

    Returns:
        List[Tuple[int, int]]: A list of tuples containing start and end coordinates of motif matches.
    """
    pattern = re.compile(f"({motif})+")
    matches = pattern.finditer(text)

    coordinates = []
    for match in matches:
        start = match.start()
        end = match.end()
        coordinates.append((start, end))

    return coordinates


def get_flexi_motifs(motif: str) -> Tuple[str, str]:
    # Count consecutive runs of characters in a string
    motif_fwd_runs = count_continuous_runs(motif)
    # Repeat fro reverse motif
    motif_rev_runs = count_continuous_runs(revComp(motif))

    # Return the regex pattern for the fwd and rev motifs that allows +/- 1 bp variation
    # in any nucleotide run >= 2 bp
    return (
        construct_regex_pattern(motif_fwd_runs),
        construct_regex_pattern(motif_rev_runs),
    )
