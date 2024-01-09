# compute_coverage
Small python script to compute coverage of contigs from reads. compute_coverage uses minimap2 to map the reads to the assembly and deduces the coverage of each contig.

## Usage

```
usage: compute_coverage.py [-h] -i INPUT -r READS -o OUTPUT [-t THREADS]
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
