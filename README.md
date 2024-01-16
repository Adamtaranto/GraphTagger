# GraphTagger
Collection of small python tools to edit segment tags in assembly graph (.gfa) files. 

Sequences in GFA files are called Segments and may have meta-data stored in tags with the format `TAG:TYPE:VALUE`. For example, mean read depth of 100 stored as a float would be: `DP:f:100`.

See the [GFA-Spec](https://gfa-spec.github.io/GFA-spec/). 


## Tools 

**map2tag**: Uses Minimap2 to map reads to the assembly and approximates the coverage of each contig. Updates the `DP` depth tag.

**csv2tag**: Adds or updates gfa segment tags from a CSV. 
CSV lines have the format `Segment`,`Tag`,`Type`,`Value`.

[DEV] **mos2tag**: Reads depth info from a [MosDepth](https://github.com/brentp/mosdepth) summary.txt file and updates `DP` tags in a corresponding `GFA`.

**fa2gfa**: Convert FASTA to GFA format. Produces Segment records for each FASTA record. Adds `LN` tags.

**tel2bed**: Quick annotation of telomeric repeat runs in BED format. This is useful for visualising potential chromosome ends in [Bandage-NG](https://github.com/asl/BandageNG)

### Converting from GFA to FASTA

**Support long sequence names**
```bash
awk '/^S/{print ">"$2; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' in.gfa > out.fa
```

**Include GFA tags in FASTA header line**

If you also want to preserve the GFA segment tags (i.e. segment length `LN`, or depth `DP`) in the FASTA description you can append them to the headers:

```bash
awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' in.gfa > out.fa
```

## Usage

### map2tag

Approximate coverage from mapped long reads

```bash
map2tag -i input.gfa -r nanopore_reads.fq.gz -o output.gfa -t 4 -x map-ont
```

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