import argparse
import logging

# Read coverage values from mosdepth summary txt file
# `{prefix}.mosdepth.summary.txt`
# chrom length bases mean min max
# chr13 115169878 8508921 0.07 0 8833
# chr13_region 13390 10009503 747.54 0 8833
# chr17 81195210 5713210 0.07 0 8689

# Extract col four header and check if text (print if median or mean)
# Extract col four value per segment name
 
# Update gfa.

def getArgs():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Add tags to GFA from mosDepth summary file.", prog="mos2tag"
    )
    parser.add_argument(
        "-i",
        "--input_gfa",
        required=True,
        help="Path to the input GFA file (can be gzipped).",
    )
    parser.add_argument(
        "-m",
        "--input_mos",
        required=True,
        help="Path to the input mosDepth summary file containing tag information.",
    )
    parser.add_argument(
        "-o",
        "--output_gfa",
        default=None,
        help="Path to the output GFA file. If not provided, "
        "the output file will have the same basename as the input with the 'tagged.gfa' extension.",
    )
    return parser.parse_args()

def main():
    # Set up logging
    logging.basicConfig(
        level=0,
        format="%(asctime)s:%(levelname)s:%(module)s:%(message)s",
    )

    # Parse command line arguments
    args = getArgs()

    print("This tool is under development.")
    # Validate the input file
    

    # Load mosDepth summary data
    

    # Update gfa tags with tags mosDepth summary file


if __name__ == "__main__":
    main()