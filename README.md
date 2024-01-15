# GraphTagger
Collection of small python tools to edit segment tags in assembly graph (.gfa) files. 

Sequences in GFA files are called Segments and may have meta-data stored in tags with the format `TAG:TYPE:VALUE`. For example, mean read depth of 100 stored as a float would be: `DP:f:100`.

See the [GFA-Spec](https://gfa-spec.github.io/GFA-spec/). 


## Tools 

**map2gfa**: Uses Minimap2 to map reads to the assembly and approximates the coverage of each contig. Updates the `DP` depth tag.

[DEV] **csv2gfa**: Adds or updates gfa segment tags from a CSV. 
Lines have the format `Segment`,`Tag`,`Type`,`Value`.

[DEV] **mos2gfa**: Reads depth info from a [MosDepth](https://github.com/brentp/mosdepth) summary.txt file and updates `DP` tags in a corresponding `GFA`.

**fa2gfa**: Convert FASTA to GFA format. Produces Segment records for each FASTA record. Adds `LN` tags.

**tel2bed**: Quick annotation of telomeric repeat runs. Input: FASTA file, Telomeric motif (i.e. TTAGGG), min repeat count; Output: BED file.

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

```
usage: map2gfa [-h] -i INPUT -r READS -o OUTPUT [-t THREADS]
                           [-x PRESET] [--minimap2 MINIMAP2]

Add depth information to an assembly

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input assembly (fasta or gfa)
  -r READS, --reads READS
                        Reads used to generate the assembly (fasta or fastq,
                        can be gzipped)
  -o OUTPUT, --output OUTPUT
                        Output assembly (same format as input)
  -t THREADS, --threads THREADS
                        Number of threads to use [1]
  -x PRESET, --preset PRESET
                        Minimap2 preset to use [map-ont]
  --minimap2 MINIMAP2   path to minimap2 executable [minimap2]

```
