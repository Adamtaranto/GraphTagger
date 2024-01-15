# author: Roland Faure
# small python script to take an assembly (fasta or gfa) and output an assembly in the same format with depth information

import sys
import argparse
import gzip
import re
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Add depth information to an assembly")
    parser.add_argument(
        "-i", "--input", help="Input assembly (fasta or gfa)", required=True
    )
    parser.add_argument(
        "-r",
        "--reads",
        help="Reads used to generate the assembly (fasta or fastq, can be gzipped)",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="Output assembly (same format as input)", required=True
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use [1]",
        required=False,
        default=1,
    )
    parser.add_argument(
        "-x",
        "--preset",
        help="Minimap2 preset to use [map-ont]",
        required=False,
        default="map-ont",
    )
    parser.add_argument(
        "--minimap2",
        help="path to minimap2 executable [minimap2]",
        required=False,
        default="minimap2",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    minimap = args.minimap2
    if minimap == "":
        minimap = "minimap2"
    assembly = args.input
    reads = args.reads
    output = args.output
    is_gfa = False

    # check the minimap2 dependency
    command = args.minimap2 + " --version"
    if os.system(command) != 0:
        print(
            "ERROR: minimap2 not found with path: ",
            args.minimap2,
            "\nPlease specify the path to minimap2 with the --minimap2 option",
        )
        sys.exit(1)

    # check the input assembly
    if not os.path.isfile(assembly):
        print("ERROR: input assembly not found")
        sys.exit(1)

    # check the input reads
    if not os.path.isfile(reads):
        print("ERROR: input reads not found")
        sys.exit(1)

    # if the assembly is in gfa format, convert it to fasta
    if assembly.endswith(".gfa"):
        is_gfa = True
        fasta = assembly + ".fasta"
        print("Converting gfa to fasta file ", fasta)
        f = open(assembly, "r")
        g = open(fasta, "w")
        for line in f:
            if line.startswith("S"):
                line = line.split("\t")
                g.write(">" + line[1] + "\n")
                g.write(line[2] + "\n")
        f.close()
        assembly = fasta

    # now we can run minimap2
    print("Running minimap2")
    paffile = output + ".tmp.paf"
    command = (
        minimap
        + " -t "
        + str(args.threads)
        + " --secondary=no -x "
        + args.preset
        + " "
        + assembly
        + " "
        + reads
        + " > "
        + paffile
    )
    if os.system(command) != 0:
        print("ERROR: minimap2 failed while running: ", command)
        sys.exit(1)

    # now we can parse the paf file and determine the coverage of each contig
    print("Parsing minimap2 output")
    contig2depth = {}
    f = open(paffile, "r")
    for line in f:
        ls = line.split("\t")
        if len(ls) < 11:
            continue
        contig = ls[5]
        length_of_contig = int(ls[6])
        start = int(ls[7])
        end = int(ls[8])
        if contig not in contig2depth:
            contig2depth[contig] = 0
        contig2depth[contig] += float(end - start) / float(length_of_contig)

    # now we can write the output
    print("Writing output")
    if is_gfa:
        f = open(args.input, "r")
        g = open(output, "w")
        for line in f:
            if line.startswith("S"):
                depth = 0
                contig = line.split("\t")[1]
                if contig in contig2depth:
                    depth = contig2depth[contig]
                g.write(line.rstrip("\n").rstrip("\t") + "\tDP:f:" + str(depth) + "\n")
            else:
                g.write(line)
        f.close()
        g.close()
    else:
        f = open(assembly, "r")
        g = open(output, "w")
        for line in f:
            if line.startswith(">"):
                contig = line[1:].rstrip("\n").split(" ")[0]
                depth = 0
                if contig in contig2depth:
                    depth = contig2depth[contig]
                g.write(line.rstrip("\n") + " DP:f:" + str(depth) + "\n")
            else:
                g.write(line)
        f.close()
        g.close()

    # remove the temporary paf file
    os.remove(paffile)

    print("Done")


if __name__ == "__main__":
    main()
