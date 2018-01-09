import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from extract_genomes import extract_genomes, get_refseq
from collections import defaultdict
from operator import add
import seaborn as sns

def run(args):
    # Define the color gradient and figure formatting.
    cs = dict(zip(range(args.cn2 - args.cn1 + 2), 
            sns.color_palette('Reds', args.cn2 - args.cn1 + 2)))
    sns.set_style('white')
    ind = np.arange(len(args.bams))
    width = 0.9

    refseq = get_refseq(args.ref) 
    ax = plt.gca()

    sample_names = []
    all_cn_proportions = defaultdict(list)
    for bam in args.bams:
        sample_names.append(bam)
        # Get a list of vaccinia genomes.
        genomes = extract_genomes(bam, refseq, af=0, allele_filter='soft').genomes
        for cn in range(args.cn1, args.cn2 + 1):
            # Count the number of genomes of the specified copy number.
            cn_total = sum([1 if len(genome) == cn else 0 for genome in genomes])
            # Calculate the proportion of genomes of that copy number
            # in the population overall.
            cn_proportion = cn_total / float(len(genomes))
            # Add proportions to a list keyed on the BAM file being analyzed.
            all_cn_proportions[bam].append(cn_proportion)
    # Since we're plotting stacked bars, keep a running tally of the
    # "bottom" value we're adding the next bar onto.
    running_total = [0. for i in range(len(args.bams))]
    # Now, instead of looping over samples, loop over each copy number,
    # and get the proportions of genomes in each sample that match the
    # copy number of interest.
    for cn in range(args.cn1, args.cn2 + 1):
        idx = cn - args.cn1
        cn_props = [all_cn_proportions[bam][idx] for bam in all_cn_proportions]
        label = "CN {}".format(cn)
        color = cs[cn - args.cn1]
        plt.bar(ind, cn_props, width, bottom=running_total, 
                label=label, color=color, lw=0.25, edgecolor='w')
        # Keep track of the cumulative proportions being added to
        # the stacked bar plots.
        running_total = map(add, running_total, cn_props)
    # If the proportion of genomes in the specified range of copy numbers
    # doesn't add to 1, fill in the rest of the stacked bar plot.
    if any([x < 1.0 for x in running_total]):
        remaining_proportions = [1.0 - x for x in running_total]
        label = "CN > {}".format(args.cn2 + 1)
        color = cs[args.cn2 + 1 - args.cn1]
        plt.bar(ind, remaining_proportions, width, bottom=running_total, 
                label=label, color=color, lw=0.25, edgecolor='w')
    plt.ylabel("Proportion of reads")
    plt.xlabel('Passage')
    plt.legend()
    ax.set_xticks(ind)
    sns.despine()
    if args.png:
        plt.savefig(args.o + '.png', format='png', bbox_inches='tight')
    else:
        plt.savefig(args.o + '.eps', format='eps', bbox_inches='tight')

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bams", required=True, nargs='+', help='Path to sorted BAM file(s).')
    p.add_argument("--ref", required=True, help='Path to FASTA reference genome.')
    p.add_argument("-o", help='Name of output plot.', default='hist')
    p.add_argument("-png", help='Output as PNG.', action='store_true')
    p.add_argument("-cn1", type=int, help='Lower bound on copy number to plot', default=1)
    p.add_argument("-cn2", type=int, help='Upper bound on copy number to plot', default=5)
    run(p.parse_args(argv))
if __name__ == "__main__":
    main(sys.argv[1:])
