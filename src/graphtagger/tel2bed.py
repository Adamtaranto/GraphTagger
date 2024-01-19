from graphtagger.logs import CustomFormatter
from graphtagger.utils import is_valid_fasta_file
from graphtagger.motifs import get_flexi_motifs, find_repeats_of_motif
from graphtagger.seqOps import revComp

from Bio import SeqIO

from typing import List, Tuple, Optional
import argparse
import gzip
import logging
import os.path
import sys

# TASK: Choose either "-" strand or rev comp "name" for output bed.
# TASH: Support GFA as input format.

def process_fasta(
    input_file: str, output_file: str, motif: str, minrep: Optional[int] = 1
) -> None:
    """
    Process a FASTA file, search for motif repeats, and write results to a BED file.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output BED file.
        motif (str): The motif to search for in each sequence.
    """

    # Validate the input file
    if not is_valid_fasta_file(input_file):
        return

    # Determine if the input file is gzipped
    is_gzipped = input_file.endswith(".gz")

    # If output_bed is not provided, derive the output filename based on input
    if output_file is None:
        if is_gzipped:
            file_root, base_ext = os.path.splitext(input_file)
            output_file = os.path.splitext(file_root)[0] + ".bed"
        else:
            output_file = os.path.splitext(input_file)[0] + ".bed"

    # Min length of sequential pattern matchs to report
    minreplen = len(motif) * minrep

    # Set default score value
    score = "."

    # Generate flexible fwd/rev regex patterns from base motif
    motifs = get_flexi_motifs(motif)
    revMotif = revComp(motif)

    # Open output bed file
    logging.info(f"Writing gfa to file: {output_file}")
    with open(output_file, "w") as output_bed:
        output_bed.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")

        logging.info(f"Reading seq records from: {input_file}")
        with gzip.open(input_file, "rt") if is_gzipped else open(
            input_file, "r"
        ) as input_handle:
            seq_count = 0
            pos_seq_count = 0
            total_hits = 0
            for record in SeqIO.parse(input_handle, "fasta"):
                logging.info(f"Searching sequence: {record.id}")

                # Search for forward orientation of flexi motif
                hits = find_repeats_of_motif(str(record.seq), motifs[0])
                strand = "+"
                fwd_count = 0
                for start, end in hits:
                    if int(end) - int(start) > minreplen:
                        score = int(end) - int(start)
                        output_bed.write(
                            f"{record.id}\t{start}\t{end}\t{motif}\t{score}\t{strand}\n"
                        )
                        fwd_count += 1
                if fwd_count > 0:
                    logging.info(f"Fwd motif runs found: {fwd_count}")

                # Search for reverse orientation of flexi motif
                hits = find_repeats_of_motif(str(record.seq), motifs[1])
                strand = "-"
                rev_count = 0
                for start, end in hits:
                    if int(end) - int(start) > minreplen:
                        score = int(end) - int(start)
                        output_bed.write(
                            f"{record.id}\t{start}\t{end}\t{revMotif}\t{score}\t{strand}\n"
                        )
                        rev_count += 1
                if rev_count > 0:
                    logging.info(f"Rev motif runs found: {rev_count}")

                # Increment count of motif positive sequences if any hits found
                if fwd_count or rev_count:
                    pos_seq_count += 1
                    total_hits += fwd_count + rev_count

                # Increment seq counter
                seq_count += 1

    logging.info(
        f"Screened {seq_count} seq records for motifs: Fwd={motifs[0]} and Rev={motifs[1]}."
    )
    logging.info(f"Found {total_hits} in {pos_seq_count} sequences.")


def getArgs():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Quick annotation of telomeric repeat runs in fasta file.", prog="tel2bed"
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        required=True,
        help="Path to the input FASTA file (can be gzipped).",
    )
    parser.add_argument(
        "-o",
        "--output_bed",
        default=None,
        help="Path to the output BED file. If not provided, "
        "the output file will have the same basename as the input with the '.bed' extension.",
    )
    parser.add_argument(
        "-m",
        "--motif",
        required=True,
        help="Telomeric motif to annotate in genome. By default tel2bed will search in fwd and rev orientations and will allow +/- 1bp flexibility in the pattern at any polynucleotide run.",
    )
    parser.add_argument(
        "-r",
        "--min_repeats",
        default=3,
        type=int,
        help="Minimum number of sequential pattern matches required for a hit to be reported. Default: 3",
    )

    # Parse command line arguments
    return parser.parse_args()


def main():
    # Set up logging
    fmt = "%(asctime)s | %(levelname)8s | %(module)s:%(lineno)s:%(funcName)20s() | %(message)s"
    handler_sh = logging.StreamHandler(sys.stdout)
    handler_sh.setFormatter(CustomFormatter(fmt))
    logging.basicConfig(format=fmt, level=logging.INFO, handlers=[handler_sh])

    # Parse command line arguments
    args = getArgs()

    # Convert FASTA to GFA
    process_fasta(
        args.input_fasta, args.output_bed, args.motif, minrep=args.min_repeats
    )


if __name__ == "__main__":
    main()
