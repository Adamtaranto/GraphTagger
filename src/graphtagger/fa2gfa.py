from typing import Optional
import argparse
import gzip
import hashlib
import logging
import os

from Bio import SeqIO

from graphtagger.logs import init_logging
from graphtagger.utils import is_valid_fasta_file


# gfa2fa
# Support headers longer than 80 chars, no tags
# awk '/^S/{print ">"$2; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' test.gfa > out.fa

# Allow for long seq names > 80 chars + include tags in header.
# awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' test.gfa > out.fa


def getArgs():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Convert FASTA file to GFA format.", prog="fa2gfa"
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        required=True,
        help="Path to the input FASTA file (can be gzipped).",
    )
    parser.add_argument(
        "-o",
        "--output_gfa",
        default=None,
        help="Path to the output GFA file. If not provided, "
        "the output file will have the same basename as the input with the '.gfa' extension.",
    )
    parser.add_argument(
        "-s",
        "--calc_hash",
        default=False,
        action="store_true",
        help="If set, calculate new SH tags from sha256 hash of sequence.",
    )

    # Parse command line arguments
    return parser.parse_args()


def convert_fasta_to_gfa(
    input_fasta: str, output_gfa: Optional[str] = None, calcHash: bool = False
) -> None:
    """
    Convert FASTA file to GFA format.

    Args:
        input_fasta (str): Path to the input FASTA file (can be gzipped).
        output_gfa (str, optional): Path to the output GFA file. If not provided,
            the output file will have the same basename as the input with the ".gfa" extension.
    """

    # Validate the input file
    if not is_valid_fasta_file(input_fasta):
        return

    # Determine if the input file is gzipped
    is_gzipped = input_fasta.endswith(".gz")

    # If output_gfa is not provided, derive the output filename based on input
    if output_gfa is None:
        if is_gzipped:
            file_root, base_ext = os.path.splitext(input_fasta)
            output_gfa = os.path.splitext(file_root)[0] + ".gfa"
        else:
            output_gfa = os.path.splitext(input_fasta)[0] + ".gfa"

    # Open the output GFA file for writing
    logging.info(f"Writing gfa to file: {output_gfa}")
    with open(output_gfa, "w") as output_file:
        # Open the input file accordingly (regular or gzipped)
        logging.info(f"Reading seq records from: {input_fasta}")
        with gzip.open(input_fasta, "rt") if is_gzipped else open(
            input_fasta, "r"
        ) as input_file:
            seq_count = 0
            # Read records from the FASTA file one at a time
            for seq_record in SeqIO.parse(input_file, "fasta"):
                # Add new SH tag to segment_tags dict if calcHash is True
                if calcHash:
                    # Calculate sha256 hash of the sequence
                    checksum = hashlib.sha256(str(seq_record.seq).encode()).hexdigest()
                    # Write the information from the current seqRecord to the output file in gfa segment line format
                    output_file.write(
                        f"S\t{seq_record.id}\t{str(seq_record.seq)}\tLN:i:{len(seq_record)}\tSH:H:{checksum}\n"
                    )
                else:
                    # Write the information from the current seqRecord to the output file in gfa segment line format
                    output_file.write(
                        f"S\t{seq_record.id}\t{str(seq_record.seq)}\tLN:i:{len(seq_record)}\n"
                    )
                seq_count += 1

    logging.info(f"Converted {seq_count} records to gfa.")


def main():
    # Set up logging
    init_logging()

    # Parse command line arguments
    args = getArgs()

    # Convert FASTA to GFA
    convert_fasta_to_gfa(args.input_fasta, args.output_gfa, args.calc_hash)


if __name__ == "__main__":
    main()
