## Scripts for: Long read sequencing reveals poxvirus evolution through rapid homogenization of gene arrays

Python code used to characterize copy number and single nucleotide variation in an evolving population of vaccinia virus.

This code is archived on Zenodo. [![DOI](https://zenodo.org/badge/116744275.svg)](https://zenodo.org/badge/latestdoi/116744275)

### Dependencies

`python 2.7`

`pysam (v0.12.0.1)`

`numpy (v1.12.0)`

`pyfaidx (v0.4.9.2)`

`seaborn (v0.7.1)`

### Installation

```
git clone https://github.com/tomsasani/vacv-ont-manuscript.git
cd vacv-ont-manuscript
```

### Usage

To get information about the arguments passed into a particular script, simply run:

```
$ python $script.py -h
```

Scripts produce various plots and summary statistics, as described below:

```
scripts:
summary.py       : text summary of ONT data
condensed.py     : linked boxes representing allele frequencies within copies
hist.py          : stacked bars representing copy number distributions
sb.py            : stacked bars representing frequencies of allele combinations
```

### Reproducing figures

To reproduce figures in the manuscript that were generated from ONT sequencing data, first download relevant ONT FASTQ files from the SRA (accession SRP128569) or from [Zenodo](DOI:10.5281/zenodo.1169394).

Then, align FASTQ using BWA-MEM. For example:

```
bwa mem -x ont2d -t $THREADS $REF $FASTQ | sambamba view -S -f bam /dev/stdin | sambamba sort -o $OUT.BAM /dev/stdin
```

Alternatively, download BAM files (aligned using the above command) from [Zenodo](DOI:10.5281/zenodo.1169394).

More details regarding each subcommand are given below, but commands for reproducing figures are as follows.

**Figure 2C**: `python hist.py --ref $REF --bams p5.bam p6.bam p7.bam p10.bam p15.bam p20.bam`

**Figure 3**: `python sb.py --ref $REF --bams p5.bam p6.bam p7.bam p10.bam p15.bam p20.bam`

**Figure 3-figure supplement 2**: `python simulate_sequencing_errors.py --ref $REF --bam p15.r7.bam --chemistry r7` (repeat for p15.bam and p15.MOI.0.1.bam, which were sequenced with R9 and R9.4 chemistries, respectively)

**Figure 3-figure supplement 3**: `python array_combinations.py --ref $REF --bams p10.bam p15.bam p20.bam -cn 3` 

**Figure 4**: `python condensed.py --ref $REF --bam p5.bam` (repeat for p10, p15, and p20 BAM files)

**Figure 6**: `python sb.py --ref $REF --bams p15.MOI.1.bam p15.MOI.0.1.bam p15.MOI.0.01.bam p15.MOI.0.001.bam`

**Figure 6A**: `python hist.py --ref $REF --bams p10.bam p15.bam p15.bhk.bam`

**Figure 6C**: `python sb.py --ref $REF --bams p10.bam p15.bam p15.bhk.bam`

#### `summary.py`

```
python summary.py --ref reference> \
                  --bam sorted and indexed BAM file
```

This will print a summary of K3L CNV and K3L<sup>H47R</sup> allele frequencies in a given population. Allele frequencies are calculated from BAM pileups.

For example, output might look like this for alignments from P5...
```
locus	variant_site	ref	alt	af	coverage (# reads)	coverage (depth)	cn_distrib
K3L   30490   T    C    0.1366    1190  2663    1:787 | 2:206 | 3:122 | 4:45 | 5:17 | 6:8 | 7:3 | 8:1 | 9:1 | 10:0 | 11:0 | 12:0 | 13:0 | 14:0 | 15:0 
```

#### `hist.py`

Produces a stacked bar histogram of copy number within one or more populations of vaccinia. As seen in Figure 2C and 6A.

```
python hist.py --ref      reference \
               --bams     sorted and indexed BAM file(s) \
               -cn1       lower bound on copy number to plot (default = 1) \
               -cn2       upper bound on copy number to plot (default = 5) \
               -png       output plot as PNG (by default, output EPS) \
               -o         name of output file (by default, "hist")
```

#### `condensed.py`

Produces "condensed" plots of K3L<sup>His47Arg</sup> allele frequency in each K3L copy within a population of vaccinia. As seen in Figure 3.

```
python condensed.py --ref      reference \
                    --bams     sorted and indexed BAM file(s) \
                    -rand      plot a random distribution of H47R alleles in the population \
                    -cn        plot genomes with up to this copy number (default = 5) \
                    -png       output plot as PNG (by default, output EPS) \
                    -o         name of output file (by default, "condensed")
```

#### `sb.py`

Produces a stacked bar representation of homogenous and "mixed" WT/K3L<sup>His47Arg</sup> arrays within one or more populations of vaccinia. As seen in Figure 3, 5, and 6C.

```
python sb.py --ref      reference \
             --bams     sorted and indexed BAM file(s) \
             -rand      plot a uniform distribution of H47R alleles in the population \
             -cn        plot genomes with up to this copy number (default = 5) \
             -png       output plot as PNG (by default, output EPS) \
             -o         name of output file (by default, "stackedbar")
```

#### `array_combinations.py`

Produces a plot that illustrates the abundance of K3L arrays with every possible combination of K3L<sup>WT</sup>/K3L<sup>His47Arg</sup> copies (for a given copy number). As seen in Figure 3-figure supplement 3.

```
python array_combinations.py --ref      reference \
                             --bams     sorted and indexed BAM file(s) \
                             -cn        plot genomes with this copy number (default = 3) \
                             -png       output plot as PNG (by default, output EPS) \
                             -o         name of output file (by default, "array-combinations")
```

#### `simulate_sequencing_errors.py`

Iterates over a population of K3L arrays and converts mixed arrays to homogenous arrays. Then, introduces "sequencing errors" by switching K3L<sup>WT</sup> and K3L<sup>His47Arg</sup> alleles.

```
python simulate_sequencing_errors.py --ref        reference \
                                     --bam        sorted and indexed BAM file(s) \
                                     --chemistry  the flowcell chemistry used to generate reads in BAM [r7, r9, r9.4] \
                                     -png         output plot as PNG (by default, output EPS) \
                                     -o           name of output file (by default, "sim-errors")
```
