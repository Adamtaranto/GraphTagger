import hashlib
import re
import logging

# Based on functions from BiopythonBio/SeqIO/GfaIO.py
# By @michaelfm1211


def _check_tag_values(name, seq, tags):
    """Check a segment line's tags for inconsistencies (PRIVATE)."""
    for tag in tags:
        if tag[:2] == "LN":
            if int(tag[5:]) != len(seq):
                logging.warning(
                    f"Segment line '{name}' has incorrect length. Expected {len(seq)} but got {tag[5:]}."
                )
        elif tag[:2] == "SH":
            # SHA256 checksum
            checksum = hashlib.sha256(str(seq).encode()).hexdigest()
            if checksum.upper() != tag[5:]:
                logging.warning(
                    f"Segment line '{name}' has incorrect checksum. Expected {checksum} but got {tag[5:]}.",
                )


def _tags_to_annotations(tags):
    """Build an annotations dictionary from a list of tags (PRIVATE)."""
    annotations = {}
    for tag in tags:
        parts = tag.split(":")
        if len(parts) < 3:
            logging.error(f"Segment line has invalid tag: {tag}.")
        # Tag name must only be two alphanumeric characters
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9]", parts[0]) is None:
            logging.warning(
                f"Tag has invalid name: {parts[0]}. Are tags tab delimited?",
            )
        # Reassemble downstream splits in the value in case it contained ":" characters.
        parts[2] = ":".join(parts[2:])  # tag value may contain : characters

        # Add annotation to dict
        annotations[parts[0]] = (parts[1], parts[2])

        # Check type of the tag and raise warning on a mismatch. These RegExs
        # are part of the 1.0 standard.
        if parts[1] not in "AifZJHB":
            logging.warning(f"Tag has invalid type: {parts[1]}")
        elif parts[1] == "A" and re.fullmatch(r"[!-~]", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected printable character, got {parts[2]}."
            )
        elif parts[1] == "i" and re.fullmatch(r"[-+]?[0-9]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected signed integer, got {parts[2]}."
            )
        elif (
            parts[1] == "f"
            and re.fullmatch(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", parts[2])
            is None
        ):
            logging.warning(f"Tag has incorrect type. Expected float, got {parts[2]}.")
        elif parts[1] == "Z" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected printable string, got {parts[2]}."
            )
        elif parts[1] == "J" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected JSON excluding new-line and tab characters, got {parts[2]}."
            )
        elif parts[1] == "H" and re.fullmatch(r"[0-9A-Fa-f]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected byte array in hex format, got {parts[2]}."
            )
        elif (
            parts[1] == "B"
            and re.fullmatch(
                r"[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+", parts[2]
            )
            is None
        ):
            logging.warning(
                f"Tag has incorrect type. Expected array of integers or floats, got {parts[2]}.",
            )
    return annotations


def _validate_tags(tags):
    """Build an annotations dictionary from a list of tags (PRIVATE)."""
    # Init status checker
    tag_is_vaild = True

    for tag in tags:
        parts = tag.split(":")
        if len(parts) < 3:
            logging.error(f"Segment line has invalid tag: {tag}.")
            tag_is_vaild = False
        # Tag name must only be two alphanumeric characters
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9]", parts[0]) is None:
            logging.warning(
                f"Tag has invalid name: {parts[0]}. Are tags tab delimited?",
            )
            tag_is_vaild = False
        # Reassemble downstream splits in the value in case it contained ":" characters.
        parts[2] = ":".join(parts[2:])  # tag value may contain : characters

        # Log error if ":" in value of non-JSON tag
        if ":" in parts[2] and parts[1] != "J":
            logging.error(
                f"Non-JSON-type tag contains character ':' in value field. This is probably an error!: {tag}"
            )
            tag_is_vaild = False
        # Check type of the tag and raise warning on a mismatch. These RegExs
        # are part of the 1.0 standard.
        if parts[1] not in "AifZJHB":
            logging.warning(f"Tag has invalid type: {parts[1]}")
            tag_is_vaild = False
        elif parts[1] == "A" and re.fullmatch(r"[!-~]", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected printable character, got {parts[2]}."
            )
            tag_is_vaild = False
        elif parts[1] == "i" and re.fullmatch(r"[-+]?[0-9]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected signed integer, got {parts[2]}."
            )
            tag_is_vaild = False
        elif (
            parts[1] == "f"
            and re.fullmatch(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", parts[2])
            is None
        ):
            logging.warning(f"Tag has incorrect type. Expected float, got {parts[2]}.")
            tag_is_vaild = False
        elif parts[1] == "Z" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected printable string, got {parts[2]}."
            )
            tag_is_vaild = False
        elif parts[1] == "J" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected JSON excluding new-line and tab characters, got {parts[2]}."
            )
            tag_is_vaild = False
        elif parts[1] == "H" and re.fullmatch(r"[0-9A-Fa-f]+", parts[2]) is None:
            logging.warning(
                f"Tag has incorrect type. Expected byte array in hex format, got {parts[2]}."
            )
            tag_is_vaild = False
        elif (
            parts[1] == "B"
            and re.fullmatch(
                r"[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+", parts[2]
            )
            is None
        ):
            logging.warning(
                f"Tag has incorrect type. Expected array of integers or floats, got {parts[2]}.",
            )
            tag_is_vaild = False

    return tag_is_vaild
