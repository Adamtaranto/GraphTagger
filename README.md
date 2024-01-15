# GraphTagger
Collection of small python tools to edit segment tags in assembly graph (.gfa) files. 

Sequences in GFA files are called Segments and may have meta-data stored in tags with the format `TAG:TYPE:VALUE`. For example, mean read depth of 100 stored as a float would be: `DP:f:100`.

See the [GFA-Spec](https://gfa-spec.github.io/GFA-spec/). 


## Tools 

**map2gfa**: Uses Minimap2 to map reads to the assembly and approximates the coverage of each contig. Updates the `DP` depth tag.

**csv2gfa**: Adds or updates gfa segment tags from a CSV. 
Lines have the format `Segment`,`Tag`,`Type`,`Value`.

**mos2gfa**: Reads depth info from a [MosDepth](https://github.com/brentp/mosdepth) summary.txt file and updates `DP` tags in a corresponding `GFA`.

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
