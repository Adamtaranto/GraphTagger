from graphtagger.logs import init_logging

from typing import List, Tuple
import argparse
import logging
import os.path
import re
import sys

# GAF format: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
# GFA Path format: https://gfa-spec.github.io/GFA-spec/GFA1.html#required-fields-5
# Bandage paths: https://github.com/rrwick/Bandage/wiki/Graph-paths

# TODO: add default outfile naming
# TODO: Support input GFA file + append new paths to file end

def count_symbols(input_string: str) -> int:
    """
    Count the number of instances of ">" or "<" in a string.

    Args:
        input_string (str): The input string.

    Returns:
        int: The number of instances of ">" or "<" in the string.
    """
    count = 0
    
    count = input_string.count('>') + input_string.count('<')
    
    return count


def format_path_string(input_string: str) -> str:
    """
    Format a string containing numbers preceded by ">" or "<" into a comma-delimited string.

    Args:
        input_string (str): The input string containing numbers and their preceding characters.

    Returns:
        str: The formatted comma-delimited string.
    """
    # Use regular expression to find pairs of numbers and their preceding characters
    matches = re.findall(r'([<>])([^<>]*)', input_string)

    # Convert the matches to tuples with ">" replaced by "+" and "<" replaced by "-"
    result_tuples: List[Tuple[str, int]] = [
        (char.replace(">", "+").replace("<", "-"), name) for char, name in matches
    ]

    # Convert the list of tuples to a string
    result_string = ", ".join([f"{num}{char}" for char, num in result_tuples])

    return result_string


def process_gaf_file(input_file: str, output_file: str) -> None:
    """
    Process a GAF file, perform checks and transformations, and write results to an output file.

    Args:
        input_file (str): The path to the input GAF file.
        output_file (str): The path to the output file.
    """
    # Initialize counters
    path_counter = 0
    line_counter = 0
    
    # Open input and output files
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            # Skip lines starting with "#"
            if line.startswith("#"):
                continue
            # Count alignment line
            line_counter += 1
            
            # Split line into columns
            columns = line.strip().split("\t")

            # Check if line has at least 12 columns
            if len(columns) >= 12:
                # Extract string in $6
                path_string = columns[5]
                
                # Check if the string has at least 2 numbers
                if count_symbols(path_string) >= 2:
                    
                    # Format the string
                    formatted_string = format_path_string(path_string)
                    
                    # Extract the read name from the query name column
                    # Split on whitespace in case there is nanopore metadata in the name
                    read_name = columns[0].strip().split()[0]
                    
                    # Increment path counter
                    path_counter += 1
                    
                    # Write new tab-delimited line to output file
                    outfile.write(
                        f"P\tPath_{path_counter}:{read_name}\t{formatted_string}\t*\n"
                    )
            else:
                logging.warning(f'Skipping malformed line:\n{line}')
    # Summary            
    logging.info(f'Read {line_counter} GAF alignments.')
    logging.info(f'Extracted {path_counter} paths with > 1 segment.')
                        


def getArgs():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Convert GAF to GFA paths.", prog="gaf2path"
    )
    # parser.add_argument(
    #    "-i",
    #    "--input_gfa",
    #    required=True,
    #    help="Path to the input GFA file (can be gzipped).",
    # )
    parser.add_argument(
        "-g",
        "--input_gaf",
        required=True,
        help="Path to the input GAF file containing path information.",
    )
    parser.add_argument(
        "-o",
        "--output_gfa",
        default="output.paths.gfa",
        help="Path to the output GFA file. If not provided, "
        "the output file will have the same basename as the input with the 'paths.gfa' extension.",
    )
    # Parse command line arguments
    return parser.parse_args()


def main():
    # Set up logging
    init_logging()

    # Parse command line arguments
    args = getArgs()

    # Check if GAF exists
    if not os.path.exists(args.input_gaf):
        logging.error(f"Input file '{args.input_gaf}' does not exist.")
        sys.exit(1)
        
    # Make default outfile name
    
    # Generate GFA path lines from GAF file
    process_gaf_file(
        args.input_gaf,
        args.output_gfa,
    )


if __name__ == "__main__":
    main()

# """
# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start (0-based; closed)
# 4	int	Query end (0-based; open)
# 5	char	Strand relative to the path: + or -
# 6	string	Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
# 7	int	Path length
# 8	int	Start position on the path (0-based)
# 9	int	End position on the path (0-based)
# 10	int	Number of residue matches
# 11	int	Alignment block length
# 12	int	Mapping quality (0-255; 255 for missing)
# """
# 
# """
# Column	Field	Type	Regexp	Description
# 1	RecordType	Character	P	Record type
# 2	PathName	String	[!-)+-<>-~][!-~]*	Path name
# 3	SegmentNames	String	[!-)+-<>-~][!-~]*	A comma-separated list of segment names and orientations
# 4	Overlaps	String	\*\|([0-9]+[MIDNSHPX=])+	Optional comma-separated list of CIGAR strings
# """

# P  path_name   segment_path_list   overlaps
# P  path1   12+,33-,4+,66+  *
