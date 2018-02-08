## Scripts for: *Long read sequencing reveals Poxvirus evolution through rapid homogenization of gene arrays*

Python code used to characterize copy number and single nucleotide variation in an evolving population of vaccinia virus.

This code is archived on Zenodo. [![DOI](https://zenodo.org/badge/116744275.svg)](https://zenodo.org/badge/latestdoi/116744275)


### Dependencies

`python 2.7`

`pysam (v0.11.2.2)`

`numpy (v1.12.0)`

`pyfaidx (v0.4.9.2)`

`seaborn (v0.7.1)`

### Installation

```
git clone https://github.com/tomsasani/vacv-ont-manuscript.git
cd vacv-ont-manuscript
python setup.py install
```

### Usage

All scripts are run as subcommands from the main `vacv-plot` module. To view all available subcommands, simply run:

```
$ python -m vacv-plot

vacv-plot version: 0.0.1

commands:
summary       : text summary of ONT data
condensed     : linked boxes representing allele frequencies within copies
hist          : stacked bars representing copy number distributions
sb            : stacked bars representing frequencies of allele combinations
```

### Reproducing figures

To reproduce figures in the manuscript that were generated from ONT sequencing data, first download relevant ONT FASTQ files from the SRA (accession forthcoming).

Then, align FASTQ using BWA-MEM. For example:

```
bwa mem -x ont2d -t 8 $REF $FASTQ | sambamba view -S -f bam /dev/stdin | sambamba sort -o $OUT.BAM /dev/stdin
```

More details regarding each subcommand are given below, but commands for reproducing figures are as follows.

**Figure 2C**: `python -m vacv-plot hist --ref $REF --bams $P5.bam $P6.bam $P7.bam $P10.bam $P15.bam $P20.bam`

**Figure 3**: `python -m vacv-plot sb --ref $REF --bams $P5.bam $P6.bam $P7.bam $P10.bam $P15.bam $P20.bam`

**Figure 4**: `python -m vacv-plot condensed --ref $REF --bam $P5.bam` (repeat for P10, P15, and P20 BAM files)

**Figure 5**: `python -m vacv-plot sb --ref $REF --bams $P15.MOI.1.bam $P15.MOI.0.1.bam $P15.MOI.0.01.bam $P15.MOI.0.001.bam`

**Figure 6A**: `python -m vacv-plot hist --ref $REF --bams $P10.bam $P15.bam $P15.BHK.bam`

**Figure 6C**: `python -m vacv-plot sb --ref $REF --bams $P10.bam $P15.bam $P15.BHK.bam`

#### `summary`

```
python -m vacv-plot summary --ref reference> \
                            --bam sorted BAM file
```

This will print a summary of K3L CNV and K3L<sup>H47R</sup> allele frequencies in a given population. Allele frequencies are calculated from BAM pileups.

For example, output might look like this for alignments from P5...
```
locus	variant_site	ref	alt	af	coverage (# reads)	coverage (depth)	cn_distrib
K3L   30490   T    C    0.1366    1190  2663    1:787 | 2:206 | 3:122 | 4:45 | 5:17 | 6:8 | 7:3 | 8:1 | 9:1 | 10:0 | 11:0 | 12:0 | 13:0 | 14:0 | 15:0 
```

#### `hist`

Produces a stacked bar histogram of copy number within one or more populations of vaccinia. As seen in Figure 2C and 6A.

```
python -m vacv-plot hist --ref      reference \
                         --bams     sorted BAM file(s) \
                         -cn1       lower bound on copy number to plot (default = 1) \
                         -cn2       upper bound on copy number to plot (default = 5) \
                         -png       output plot as PNG (by default, output EPS) \
                         -o         name of output file (by default, "hist")
```

#### `condensed`

Produces "condensed" plots of K3L<sup>His47Arg</sup> allele frequency in each K3L copy within a population of vaccinia. As seen in Figure 3.

```
python -m vacv-plot condensed --ref      reference \
                              --bams     sorted BAM file(s) \
                              -rand      plot a uniform distribution of H47R alleles in the population \
                              -cn        plot genomes with up to this copy number (default = 5) \
                              -png       output plot as PNG (by default, output EPS) \
                              -o         name of output file (by default, "hist")
```

#### `sb`

Produces a stacked bar representation of homogenous and "mixed" WT/K3L<sup>His47Arg</sup> genomes within one or more populations of vaccinia. As seen in Figure 3, 5, and 6C.

```
python -m vacv-plot sb --ref      reference \
                       --bams     sorted BAM file(s) \
                       -rand      plot a uniform distribution of H47R alleles in the population \
                       -cn        plot genomes with up to this copy number (default = 5) \
                       -png       output plot as PNG (by default, output EPS) \
                       -o         name of output file (by default, "hist")
```
