from graphtagger.logs import CustomFormatter
from graphtagger.utils import is_valid_gfa_file

from collections import defaultdict
from typing import Dict, Tuple
import argparse
import csv
import gzip
import hashlib
import logging
import os.path
import sys

# TASK: Report total number of segments with updated tags.
# TASK: Support tab-delimited input.

def load_tags_from_csv(csv_file: str) -> defaultdict:
    """
    Load tag information from a CSV file into a nested dictionary.

    Parameters:
    - csv_file (str): Path to the CSV file. Columns: Sequence name, Tag, Type, Value.

    Returns:
    - defaultdict: A nested dictionary where keys are sequence names, and values are dictionaries
            with keys as tag names and values as tuples (Type, Value).
    """

    if not os.path.exists(csv_file):
        logging.error(f"Input file '{csv_file}' does not exist.")
        exit(1)

    # Use defaultdict to simplify dictionary construction
    tag_dict = defaultdict(dict)
    tag_counts = defaultdict(int)

    logging.info(f"Loading new tag data from: {csv_file}")
    with open(csv_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for row in reader:
            # Skip lines starting with "#" and check for missing values
            if not row or row[0].startswith("#"):
                continue
            elif len(row) != 4 or "" in row:
                # Log an error and skip the record
                logging.warning(f"Skipping record due to missing values - {row}")
                continue

            # Load line data
            sequence_name = row[0]
            tag = row[1]
            tag_type = row[2]
            value = row[3]

            # Construct the tuple (Type, Value)
            tag_info = (tag_type, value)

            # Check if the tag already exists for the sequence_name
            if tag in tag_dict[sequence_name]:
                # Log a warning and do not overwrite the existing tag_info
                logging.warning(
                    f"Tag '{tag}' already loaded for sequence '{sequence_name}'. Skipping csv duplicate."
                )
            else:
                # Add the tag info (TYPE,VALUE) to the internal TAG dictionary for SEQUENCE
                tag_dict[sequence_name][tag] = tag_info
                # Track instances of TAG:TYPE combinations
                tag_counts[f"{tag}:{tag_type}"] += 1

    # Log statistics
    logging.info(f"Number of segments with tag data loaded from csv: {len(tag_dict)}")

    logging.info("Unique combinations of TAG:TYPE and their occurrences in csv:")
    for tag_type_combination, count in tag_counts.items():
        logging.info(f"\"{tag_type_combination}\" x {count}")

    return tag_dict


def format_tags(sub_dict: Dict[str, Tuple[str, str]]) -> str:
    """
    Format the sub_dict into a string of tab-delimited tags.

    Parameters:
    - sub_dict (Dict[str, Tuple[str, str]]): The sub-dictionary for one sequence_name.

    Returns:
    - str: The formatted tags string.
    """
    formatted_tags = []
    for tag_name, (tag_type, value) in sub_dict.items():
        formatted_tag = f"{tag_name}:{tag_type}:{value}"
        formatted_tags.append(formatted_tag)
    return "\t".join(formatted_tags)


def update_tags(
    tag_dict1: Dict[str, Dict[str, Tuple[str, str]]],
    tag_dict2: Dict[str, Dict[str, Tuple[str, str]]],
    preserve: bool = True,
) -> Dict[str, Dict[str, Tuple[str, str]]]:
    """
    Update tag data for each sequence in the second tag_dict using data from the first tag_dict.

    Parameters:
    - tag_dict1 (Dict[str, Dict[str, Tuple[str, str]]]): The first tag dictionary.
    - tag_dict2 (Dict[str, Dict[str, Tuple[str, str]]]): The second tag dictionary.
    - preserve (bool): If True, preserve existing values in tag_dict2; if False, overwrite existing values.

    Returns:
    - Dict[str, Dict[str, Tuple[str, str]]]: The updated tag_dict2.
    """

    for sequence_name in tag_dict2:
        # Check if the sequence exists in the first dictionary
        if sequence_name in tag_dict1:
            tags1 = tag_dict1[sequence_name]
            tags2 = tag_dict2[sequence_name]

            for tag_name, (type1, value1) in tags1.items():
                # Check if the tag exists in the second dictionary
                if tag_name in tags2:
                    if not preserve:
                        # Overwrite the value in the second dictionary
                        tags2[tag_name] = (type1, value1)
                        logging.info(
                            f"Overwriting value for tag '{tag_name}' in sequence '{sequence_name}'"
                        )
                    else:
                        logging.info(
                            f"Preserve=True: Retain existing tag '{tag_name}' in sequence '{sequence_name}'"
                        )
                else:
                    # Add the new tag to the second dictionary
                    tags2[tag_name] = (type1, value1)
                    logging.info(
                        f"Adding new tag '{tag_name}' to sequence '{sequence_name}'"
                    )

    return tag_dict2


def update_gfa_tags(
    input_file: str,
    output_file: str,
    tag_dict: Dict[str, Dict[str, Tuple[str, str]]],
    preserve: bool = True,
    calcLen: bool = False,
    calcHash: bool = False,
):
    """
    Read a GFA file, update the Segment lines with information from tag_dict, and write to an output file.

    Parameters:
    - input_file (str): Path to the input GFA file.
    - output_file (str): Path to the output file.
    - tag_dict (Dict[str, Dict[str, Tuple[str, str]]]): The tag dictionary for updating Segment lines.
    - preserve (bool): If True, preserve existing values in tag_dict2; if False, overwrite existing values.
    - calcLen (bool): If True, add a new tag to the segment_tags dict with tag_name = "LN",
                     tag_type = "i", and value = len(dna_sequence).
    """

    # Determine if the input file is gzipped
    is_gzipped = input_file.endswith(".gz")

    # If output_gfa is not provided, derive the output filename based on input
    if output_file is None:
        if is_gzipped:
            file_root, base_ext = os.path.splitext(input_file)
            output_file = os.path.splitext(file_root)[0] + ".tagged.gfa"
        else:
            output_file = os.path.splitext(input_file)[0] + ".tagged.gfa"

    # Open the output GFA file for writing
    logging.info(f"Writing updated gfa to file: {output_file}")
    with open(output_file, "w") as outfile:
        # Open the input file accordingly (regular or gzipped)
        logging.info(f"Reading seq records from: {input_file}")
        with gzip.open(input_file, "rt") if is_gzipped else open(
            input_file, "r"
        ) as infile:
            # Process each line from the input gfa
            for line in infile:
                line = line.strip().split("\t")

                # Check if the line is a Segment line (starts with 'S')
                if line[0] == "S":
                    sequence_name = line[1]
                    dna_sequence = line[2]

                    # Initialize segment_tags as a defaultdict of dicts
                    segment_tags = defaultdict(dict)

                    # Load tags from Column 4 onwards into the segment_tags dict
                    for tag_info in line[3:]:
                        # Warn if whitespace found in tag
                        if " " in tag_info:
                            logging.warning(
                                f"Possible malformed tag, contains whitespace, tags bust be tab-separated: \"{tag_info}\""
                            )
                        # Split tag on ":"
                        tag_parts = tag_info.split(":")
                        if len(tag_parts) == 3:
                            tag_name, tag_type, value = tag_parts
                            # Raise warning if duplicate tag name exists in input gfa
                            if tag_name in segment_tags[sequence_name]:
                                logging.warning(
                                    f"Pre-existing duplicate tag {tag_name} on segment line {sequence_name}"
                                )
                            # Add new tag to segment dict, or overwrite duplicate
                            segment_tags[sequence_name][tag_name] = (tag_type, value)
                        else:
                            # Log a warning and skip the tag_info item
                            logging.warning(f"Skipping malformed tag_info: {tag_info}")

                    # Add a new tag to the segment_tags dict if calcLen is True
                    if calcLen:
                        segment_tags[sequence_name]["LN"] = (
                            "i",
                            str(len(dna_sequence)),
                        )
                        
                    # Add new SH tag to segment_tags dict if calcHash is True
                    if calcHash:
                        checksum = hashlib.sha256(str(dna_sequence).encode()).hexdigest()
                        segment_tags[sequence_name]["SH"] = (
                            "H",
                            str(checksum),
                        )

                    # If no tags have been found for this segment, init an empty tag_dict
                    # This gives update_tags() somewhere to copy the csv_tags into
                    if sequence_name not in segment_tags:
                        segment_tags[sequence_name] = {}

                    # Update segment_tags using the tag_dict
                    updated_tags = update_tags(
                        tag_dict, segment_tags, preserve=preserve
                    )

                    # Format the updated tags into a string
                    formatted_tags = format_tags(updated_tags[sequence_name])

                    # Append the formatted tags to the end of the line
                    updated_line = f"{line[0]}\t{sequence_name}\t{dna_sequence}\t{formatted_tags}\n"
                    outfile.write(updated_line)
                else:
                    # If not a Segment line, write the line to the output file unchanged
                    outfile.write("\t".join(line) + "\n")

    logging.info("Finished updating tags.")


def getArgs():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Add tags to GFA from csv file.", prog="csv2tag"
    )
    parser.add_argument(
        "-i",
        "--input_gfa",
        required=True,
        help="Path to the input GFA file (can be gzipped).",
    )
    parser.add_argument(
        "-c",
        "--input_csv",
        required=True,
        help="Path to the input CSV file containing tag information. Format = [NAME,TAG,TYPE,VALUE].",
    )
    parser.add_argument(
        "-o",
        "--output_gfa",
        default=None,
        help="Path to the output GFA file. If not provided, "
        "the output file will have the same basename as the input with the 'tagged.gfa' extension.",
    )
    parser.add_argument(
        "-p",
        "--preserve_tags",
        default=False,
        action="store_true",
        help="If set, preserve pre-existing tags.",
    )
    parser.add_argument(
        "-l",
        "--calc_len",
        default=False,
        action="store_true",
        help="If set, calculate new LN tags from length of sequence.",
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


def main():
    # Set up logging
    fmt = "%(asctime)s | %(levelname)8s | %(module)s:%(lineno)s:%(funcName)20s() | %(message)s"
    handler_sh = logging.StreamHandler(sys.stdout)
    handler_sh.setFormatter(CustomFormatter(fmt))
    logging.basicConfig(format=fmt, level=logging.INFO, handlers=[handler_sh])

    # Parse command line arguments
    args = getArgs()

    # Validate the input file
    if not is_valid_gfa_file(args.input_gfa):
        return

    # Load CSV tags
    csv_tags = load_tags_from_csv(args.input_csv)

    # Update gfa tags with tags from csv file
    update_gfa_tags(
        args.input_gfa,
        args.output_gfa,
        csv_tags,
        preserve=args.preserve_tags,
        calcLen=args.calc_len,
        calcHash=args.calc_hash,
    )


if __name__ == "__main__":
    main()
