# Authors: Roland Faure + Adam Taranto

"""
A small python script to add depth tags (DP) to an assembly.
Input: Takes an assembly (fasta or gfa, may be gzipped) and sequencing reads (fastq/fasta). 
Maps reads using minimap2, and calculates approximate mean depth per sequence.
Outputs an assembly in the same format as input with depth information in DP tags.
"""

from graphtagger.logs import init_logging
from graphtagger.utils import (
    are_tools_available,
    is_valid_fasta_file,
    is_valid_gfa_file,
)

import argparse
import gzip
import logging
import os
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add depth tags (DP) to an assembly", prog="map2tag"
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Input assembly (fasta or gfa). May be gzipped.",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--reads",
        help="Reads used to generate the assembly (fasta or fastq, can be gzipped)",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix for output file. Output will be written in same format as input.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use. Default: [1]",
        default=1,
    )
    parser.add_argument(
        "-x",
        "--preset",
        choices=["map-pb", "map-ont", "map-hifi", "sr"],
        help="Minimap2 preset to use. Default: [map-ont]",
        default="map-ont",
    )
    parser.add_argument(
        "--minimap2",
        help="Custom path to minimap2 executable [minimap2]",
        default="minimap2",
    )
    args = parser.parse_args()
    return args


def validate_input(args, temp_dir):
    # Check if minimap is available
    if args.minimap2 == "minimap2":
        # Check if Minimap2 is available on the path
        are_tools_available(["minimap2"])
    else:
        # Check custom minimap2 location if set.
        check_minimap2_availability(args.minimap2)

    # Check the input reads
    if not os.path.isfile(args.reads):
        logging.error(f"Input reads not found: {args.reads}")
        sys.exit(1)

    # Kill if asm file does not exist
    if not os.path.isfile(args.input):
        logging.error(f"Input assembly not found {args.input}")
        sys.exit(1)

    # Determine if the input file is gzipped
    is_gzipped = args.input.endswith(".gz")

    # If the assembly is in gfa format, convert it to fasta
    if is_valid_gfa_file(args.input, silent=True):
        is_gfa = True
        fasta_path = mktempfasta(args.input, temp_dir, is_gzipped)
        # log is gfa
    elif is_valid_fasta_file(args.input):
        # check if vailid fasta
        fasta_path = args.input
        # log is fasta
        is_gfa = False
    else:
        logging.error("Input is not valid GFA or FASTA.")
        sys.exit(1)

    # If output prefix is not provided, derive the output filename based on input file
    if args.prefix is None:
        if is_gzipped:
            # Strip .gz
            file_root, base_ext = os.path.splitext(args.input)
            # Strip extension
            outfile = os.path.splitext(file_root)[0]
        else:
            # Strip extension (i.e. .fa or .gfa)
            outfile = os.path.splitext(args.input)[0]
    else:
        # Use prefix if provided
        outfile = args.prefix

    # Add suffix and same extension as input
    if is_gfa:
        outfile = outfile + ".DP_tags.gfa"
    else:
        outfile = outfile + ".DP_tags.fa"

    return fasta_path, is_gfa, is_gzipped, outfile


def check_minimap2_availability(executable_path: str) -> bool:
    """
    Check if the minimap2 executable is available.

    Args:
        executable_path (str): Path to the minimap2 executable.

    Returns:
        bool: True if minimap2 is available, False otherwise.
    """
    # Check if the file exists
    if not os.path.isfile(executable_path):
        logging.error(f"File does not exist at locataion: {executable_path}")
        return False

    # Check if the command runs successfully
    try:
        subprocess.run(
            [executable_path, "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        logging.info(f"Custom minimap2 location seems good: {executable_path}")
        return True
    except subprocess.CalledProcessError:
        logging.error(
            f"Custom minimap2 location fails to run: {executable_path} --version"
        )
        return False


def get_depth(
    fasta: str,
    mmPath: str,
    reads: str,
    threads: int = 1,
    preset: str = "map-ont",
    temp_dir: str = None,
):
    if not temp_dir:
        temp_dir = os.getcwd()
    # Init empty dict for contig depths
    contig_depths = {}
    # Set temp PAF location
    paffile = temp_dir + "/reads2fasta.tmp.paf"
    # Prepare mm2 cmd
    command = [
        mmPath,
        "-t",
        str(threads),
        "--secondary=no",
        "-x",
        preset,
        fasta,
        reads,
    ]
    # Run the command
    try:
        logging.info(f"Call: {' '.join(command)}")        
        process = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True
            )

        # Print stderr
        logging.info("=== Minimap2 STDERR ===")
        logging.info(process.stderr)
        
        # Write the result to the output file
        with open(paffile, "w") as output:
            output.write(process.stdout)

        # Check outfile in correct location
        if os.path.isfile(paffile):
            logging.info(
                f"minimap2 completed successfully. Output written to: {paffile}"
            )
        else:
            logging.error(f"Failed to create alignment file: {paffile}")

    # Catch mm2 errors
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running minimap2: {e}")
        logging.debug("=== STDERR ===")
        logging.debug(e.stderr)
        sys.exit(1)

    # Now we can parse the PAF file and determine the coverage of each contig
    logging.info("Parsing minimap2 output.")
    with open(paffile, "r") as f:
        for line in f:
            ls = line.split("\t")
            if len(ls) < 11:
                continue
            contig = ls[5]
            length_of_contig = int(ls[6])
            start = int(ls[7])
            end = int(ls[8])
            if contig not in contig_depths:
                contig_depths[contig] = 0.0
            contig_depths[contig] += float(end - start) / float(length_of_contig)

    return contig_depths


def mktempfasta(input_gfa: str, temp_dir: str = None, is_gzipped: bool = False):
    if not temp_dir:
        temp_dir = os.getcwd()
    # Location to write temp fasta
    out_fasta = f"{temp_dir}/temp.fasta"
    logging.info(f"Converting gfa to fasta file: {out_fasta}")

    # Open temp fasta location for writing
    with open(out_fasta, "w") as output_file:
        # Open the input file accordingly (regular or gzipped)
        logging.info(f"Reading seq records from: {input_gfa}")
        with gzip.open(input_gfa, "rt") if is_gzipped else open(
            input_gfa, "r"
        ) as input_file:
            # Read gfa lines
            for line in input_file:
                if line.startswith("S"):
                    line = line.split("\t")
                    output_file.write(">" + line[1] + "\n")
                    output_file.write(line[2] + "\n")
    # Return path to temp fasta
    return out_fasta


def write_tags_to_file(infile, outfile, is_gfa, is_gzipped, depthDict):
    # Note: If no reads mapped, contig may be missing from depthDict.
    # Update tags in a GFA
    if is_gfa:
        with open(outfile, "w") as output_file:
            logging.info(f"Writing updated tags to gfa: {outfile}")
            # Open the input file accordingly (regular or gzipped)
            logging.info(f"Reading seq records from: {infile}")
            with gzip.open(infile, "rt") if is_gzipped else open(
                infile, "r"
            ) as input_file:
                for line in input_file:
                    if line.startswith("S"):
                        fields = line.strip("\n").split("\t")
                        contig_name = fields[1]
                        # This contig exists in depthDict
                        if contig_name in depthDict:
                            # Check for "DP" tags in optional columns
                            dp_tags = [
                                tag for tag in fields[3:] if tag.startswith("DP")
                            ]
                            # If existing DP tag found
                            if dp_tags:
                                # If "DP" tags exist, update the first one found
                                fields[
                                    fields.index(dp_tags[0])
                                ] = f"DP:f:{depthDict[contig_name]:.2f}"
                            else:
                                # If no "DP" tags exist, add a new one at the end
                                fields.append(f"DP:f:{depthDict[contig_name]:.2f}")
                        else:  # If contig not in dict, set DP to 0.0
                            # Check for "DP" tags in optional columns
                            dp_tags = [
                                tag for tag in fields[3:] if tag.startswith("DP")
                            ]
                            # If existing DP tag found
                            if dp_tags:
                                # If "DP" tags exist, update the first one found
                                fields[fields.index(dp_tags[0])] = "DP:f:0.0"
                            else:
                                # If no "DP" tags exist, add a new one at the end
                                fields.append("DP:f:0.0")
                        # Write the line to output
                        output_file.write("\t".join(fields) + "\n")
                    else:
                        # Write non-segment line unchanged
                        output_file.write(line)

    else:  # Write tags to fasta
        logging.info(f"Writing updated tags to fasta: {outfile}")
        with open(outfile, "w") as output_file:
            # Open the input file accordingly (regular or gzipped)
            logging.info(f"Reading seq records from: {infile}")
            with gzip.open(infile, "rt") if is_gzipped else open(
                infile, "r"
            ) as input_file:
                # Read lines
                for line in input_file:
                    if line.startswith(">"):
                        contig_name = line[1:].rstrip("\n").split(" ")[0]
                        if contig_name in depthDict:
                            depth = depthDict[contig_name]
                        else:
                            logging.warning(f"No data found for sequence: {contig_name}")
                            depth = 0.0
                        output_file.write(
                            line.rstrip("\n") + " DP:f:" + f"{depth:.2f}" + "\n"
                        )
                    else:
                        output_file.write(line)

def main():
    # Set up logging
    init_logging()
    # Load args
    args = parse_args()

    # Open tempdir
    # Note: use delete=False for testing.
    with tempfile.TemporaryDirectory(prefix="temp_", dir=os.getcwd(), delete=True) as temp_dir:
        logging.info(f"Open temp dir: {temp_dir}")
        # Check input files exist, correct format, and mm is accessible
        fasta_path, is_gfa, is_gzipped, outfile = validate_input(args, temp_dir)
        # Map reads to fasta and get depth
        depthDict = get_depth(
            fasta_path, args.minimap2, args.reads, args.threads, args.preset, temp_dir
        )
        # Write tags to output
        write_tags_to_file(args.input, outfile, is_gfa, is_gzipped, depthDict)
        # log complete
        logging.info("Finished!")


if __name__ == "__main__":
    main()
