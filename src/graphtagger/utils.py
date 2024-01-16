import logging
import os.path

def is_valid_fasta_file(input_fasta: str) -> bool:
    """
    Check if the input file is a valid FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file.

    Returns:
        bool: True if the file is a valid FASTA file, False otherwise.
    """
    valid_extensions = [".fa", ".fasta", ".fna"]

    if not os.path.exists(input_fasta):
        logging.error(f"Input file '{input_fasta}' does not exist.")
        return False

    # Extract the file extension, considering .gz if present
    file_root, base_ext = os.path.splitext(input_fasta)
    _, upstream_ext = os.path.splitext(file_root)

    # Check if the file is gzipped
    is_gzipped = base_ext.lower() == ".gz"

    if is_gzipped and upstream_ext.lower() not in valid_extensions:
        logging.error(
            f"Invalid file extension for '{input_fasta}'. Supported extensions are {', '.join(valid_extensions)}. These may be followed by .gz"
        )
        return False

    # Check if the file extension is valid
    elif not is_gzipped and base_ext.lower() not in valid_extensions:
        logging.error(
            f"Invalid file extension for '{input_fasta}'. Supported extensions are {', '.join(valid_extensions)}."
        )
        return False

    return True


def is_valid_gfa_file(input_gfa: str) -> bool:
    """
    Check if the input file is a valid GFA file.

    Args:
        input_gfa (str): Path to the input GFA file.

    Returns:
        bool: True if the file is a valid GFA file, False otherwise.
    """
    valid_extensions = [".gfa"]

    if not os.path.exists(input_gfa):
        logging.error(f"Input file '{input_gfa}' does not exist.")
        return False

    # Extract the file extension, considering .gz if present
    file_root, base_ext = os.path.splitext(input_gfa)
    _, upstream_ext = os.path.splitext(file_root)

    # Check if the file is gzipped
    is_gzipped = base_ext.lower() == ".gz"

    if is_gzipped and upstream_ext.lower() not in valid_extensions:
        logging.error(
            f"Invalid file extension for '{input_gfa}'. Supported extensions are {', '.join(valid_extensions)}. These may be followed by .gz"
        )
        return False

    # Check if the file extension is valid
    elif not is_gzipped and base_ext.lower() not in valid_extensions:
        logging.error(
            f"Invalid file extension for '{input_gfa}'. Supported extensions are {', '.join(valid_extensions)}."
        )
        return False

    return True