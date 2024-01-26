# GraphTagger
Collection of python tools to edit tags in assembly graph (.gfa) files. 

Sequences in GFA files are called Segments and may have meta-data stored in tags with the format `TAG:TYPE:VALUE`. See the [GFA-Spec](https://gfa-spec.github.io/GFA-spec/).


For example, mean read depth of 100 stored as a float would be: `DP:f:100`. You can use `csv2tag` to add or update arbritrary tags.

**New features**

Got an tag related task in need of automation? [Suggest a feature](https://github.com/Adamtaranto/GraphTagger/issues).

## Installation

GraphTagger requires Python >= v3.8

Install directly from this git repository.

```bash
pip install git+https://github.com/Adamtaranto/GraphTagger.git
```

Or clone and install locally.

```bash
git clone https://github.com/Adamtaranto/GraphTagger.git && cd GraphTagger
pip install -e .
```

## Tools 

### Tools for updating segment tags

**csv2tag**: Adds or updates gfa segment tags from a CSV. 
CSV lines have the format `Name`,`Tag`,`Type`,`Value`.

**map2tag**: Uses Minimap2 to map reads to the assembly and approximates the coverage of each contig. Updates the `DP` depth tag.

### Working with Paths
[DEV] **gaf2path**: Generate GFA format Path lines from GraphAligner GAF output.

### Feature annotation tools

**tel2bed**: Quick annotation of telomeric repeat runs in FASTA file, outputs BED format. This is useful for visualising potential chromosome ends in [Bandage-NG](https://github.com/asl/BandageNG). 

### Converting between GFA and FASTA

#### Converting from FASTA to GFA

**fa2gfa**: Convert FASTA to GFA format. Produces Segment records for each FASTA record. Adds `LN` tags.

#### Converting from GFA to FASTA

To convert an assembly graph into a FASTA of contigs you can use one of these `awk` one-liners.

If you want to preserve the GFA segment tags (i.e. segment length `LN`, or depth `DP`) in the FASTA description, option #2 will append them to the FASTA headers.

```bash
# Convert GFA to FASTA
# Supports sequence names longer than the 80 char sequence wrap length.
awk '/^S/{print ">"$2; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' in.gfa > out.fa

# Include GFA tags in FASTA header line.
awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' in.gfa > out.fa
```

## Usage

### csv2tag

Add tags to GFA from csv file.

INPUT_GFA: Path to the input GFA file (can be gzipped).   
INPUT_CSV: must have format = [NAME,TAG,TYPE,VALUE]

Options:

- --preserve_tags:   
If set, preserve pre-existing tags from gfa file.

- --calc_len:        
If set, calculate new LN tags from length of sequence.

```bash
csv2tag -i input.gfa -c new_tags.csv -o output.gfa
```

### map2tag

Approximate coverage from mapped long reads.   

```bash
map2tag -i input.gfa -r nanopore_reads.fq.gz -p outdir/tagged_graph -t 4 -x map-ont
# Output: outdir/tagged_graph.gfa
```
### gaf2path

Convert graph-alignments (GAF) to GFA path-lines.

INPUT_GAF: GAF file containing path information. Produced by [GraphAligner](https://github.com/maickrau/GraphAligner).   

OUTPUT_GFA: GFA formatted path-lines. Can be manually appended to the end of the original GFA graph used to produce the GAF file.

```bash
# Align long-reads to GFA 
GraphAligner --threads 8 --multimap-score-fraction 1 -x vg -f reads.fq -g input_graph.gfa -a longreads_aligned_on_gfa.gaf

# Convert GAF to GFA path-lines
gaf2path -g longreads_aligned_on_gfa.gaf -o pathlines.gfa

# Append pathlines to end of original graph
cat pathlines.gfa >> input_graph.gfa
```

### fa2gfa

Convert a FASTA file to GFA format.

INPUT_FASTA: Path to the input FASTA file (can be gzipped).

OUTPUT_GFA: Path to the output GFA file. If not provided, the output will be same as input but with the '.gfa' extension.

```bash
fa2gfa -i assembly.fa.gz -o output assembly.gfa
```

### tel2bed

Quick annotation of telomeric repeat runs in fasta file.

INPUT_FASTA: Path to the input FASTA file (can be gzipped).

Options:

- -o,--output_bed:  
Path to the output BED file. If not provided, the output file will have the same basename as the input with the '.bed' extension.
- -m,--motif:  
Telomeric motif to annotate in genome. By default tel2bed will search in fwd and rev orientations and will allow +/- 1bp flexibility in the pattern at any polynucleotide run.

- -r,--min_repeats:   
Minimum number of sequential pattern matches required for a hit to be reported. Default: 3

```bash
tel2bed -i contigs.fa -m TTAGGG -r 3
```

## Example: Add depth tags from MosDepth.

Convert [MosDepth](https://github.com/brentp/mosdepth) summary.txt file into CSV and use csv2tag to update `DP` tags in a corresponding `GFA`.


```bash
MINIASM_FASTA="path/to/miniasm_assembly.fa"
TEMP=/tmp

# Index asm
bwa-mem2 index -p minasm_idx $MINIASM_FASTA

# Map illumina data to asm
# Use smart pairing option -p to pass paired and unpaired reads
(seqtk mergepe paired_R?.fastq.gz; zcat unpaired.fastq.gz) | bwa-mem2 mem -p -t 8 -o aligned_reads.sam minasm_idx -

# Sort and write to bam
samtools view -u aligned_reads.sam | samtools sort -@ 8 -T $TEMP -o aligned_reads.sorted.bam 
# Delete sam
rm aligned_reads.sam 

# Index the sorted bam
samtools index -@ 8 aligned_reads.sorted.bam

# Calc median depth
mkdir mosdepth
mosdepth -t 8 --use-median --no-per-base -F 1796 --by 500 mosdepth/miniasm-asm aligned_reads.sorted.bam 
```

Inspect the mosdepth summary file.

```bash
cat mosdepth/miniasm-asm.mosdepth.summary.txt 
```

output:
```
chrom	length	bases	mean	min	max
utg000001l	1258785	128682030	102.23	0	422
utg000001l_region	1258785	128682030	102.23	0	422
utg000002l	1176747	122097652	103.76	0	314
utg000002l_region	1176747	122097652	103.76	0	314
utg000003l	872965	88294465	101.14	0	434
utg000003l_region	872965	88294465	101.14	0	434
utg000004c	81755	115821276	1416.69	0	3204
utg000004c_region	81755	115821276	1416.69	0	3204
...
total	34741114	3752939904	108.03	0	15782
total_region	34741114	3752939904	108.03	0	15782
```

Convert the `*.mosdepth.summary.txt` file to CSV
```bash
# Exclude lines that begin with "chrom   length" or "total"
# Exclude lines that contain "region"
# Write NAME,TAG,TYPE,VALUE in comma delimited format
awk -F'\t' 'BEGIN {OFS=","} !/^chrom\tlength|^total|region/ {print $1, "DP", "f", $4}' mosdepth/miniasm-asm.mosdepth.summary.txt > new_tags.csv
```

Update tags with csv2tag

```bash
csv2tag -i miniasm_assembly.gfa -c new_tags.csv -o miniasm_assembly_with_tags.gfa
```
