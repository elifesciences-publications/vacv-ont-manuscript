import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from extract_genomes import extract_genomes, get_refseq, create_random_genomes
from collections import defaultdict
import argparse
import seaborn as sns

def run(args):
    # Figure formatting.
    sns.set_style('white')
    n_samps = len(args.bams)
    x_len, y_len = 5 * n_samps, 2 * n_samps
    f, axarr = plt.subplots(nrows=1, ncols=n_samps, figsize=(x_len, y_len), sharey=True)
    # Define the color schema.
    cs = dict(zip(range(5), sns.color_palette('Blues', 5)))
    bar_width = 1

    allele_combos = defaultdict()
    refseq = get_refseq(args.ref)
    for bam in args.bams:
        sample_name = bam
        # Extract the full set of VACV genomes and allele combinations 
        # in the population.
        genomes = extract_genomes(bam, refseq, allele_filter='hard').genomes
        vals, af = get_vals(genomes, args.cn)
        # Generate random population of genomes if specified.
        if args.rand:
            rand_genomes = create_random_genomes(genomes, max_cn=args.cn, af=af)
            vals, af = get_vals(rand_genomes, args.cn)
        allele_combos[sample_name] = vals

    # Generating a standalone figure legend.
    wt_patch = mpatches.Patch(color=cs[0], label='homogenous WT')
    mix_patch = mpatches.Patch(color=cs[2], label='mixed alleles')
    mut_patch = mpatches.Patch(color=cs[4], label='homogenous H47R')
    arr = 0
    # Iterate over every passage's set of VACV genomes.
    for sample in allele_combos:
        values = allele_combos[sample]
        cn_range = range(1, len(values) + 1)
        pos = 1
        sep = 0
        for cn in cn_range:
            cn_vals = [x[1] for x in values if x[0] == cn][0]
            axarr[arr].bar(pos + sep,
                    cn_vals[0],
                    width=bar_width,
                    color = cs[0]
                    )
            axarr[arr].bar(pos + sep, 
                    cn_vals[1],
                    bottom=cn_vals[0],
                    width=bar_width,
                    color = cs[2], 
                    )
            axarr[arr].bar(pos + sep,
                    cn_vals[2],
                    bottom=cn_vals[1] + cn_vals[0],
                    width=bar_width,
                    color = cs[4]
                    )
            sep += 0.25
            pos += 0.75
        axarr[arr].set_title(sample.upper())
        axarr[arr].set_xlim([0.5, int(args.cn) + 0.5])
        arr += 1
    plt.legend(handles=[wt_patch, mix_patch, mut_patch], loc='upper right')
    sns.despine(right=True, top=True)
    axarr[0].set_ylabel('Proportion of Reads')
    axarr[n_samps / 2].set_xlabel('Copy Number')
    if args.png:
        plt.savefig(args.o + '.png', format='png', bbox_inches='tight')
    else:
        plt.savefig(args.o + '.eps', format='eps', bbox_inches='tight')

def count(genomes, cn):
    """
    Used in conjuction with `get_vals`. 
    Counts the proportions of all_wt, all_mut,
    or mixed genomes in a given set of `data`,
    up to a copy number of `cn`.
    """
    all_wt, all_mut, mixed = 0, 0, 0
    for genome in genomes:
        if len(genome) != cn: continue
        # homogenous WT
        if all([i == 0 for i in genome]): all_wt += 1 
        # homogenous H47R
        elif all([i == 1 for i in genome]): all_mut += 1
        # "mixed" -- requires genome to be longer than 1 copy
        elif len(set(genome)) > 1 and all([i in (0, 1) for i in genome]):
            mixed += 1
    return all_wt, mixed, all_mut 

def get_vals(genomes, max_cn):
    """
    Method to get the frequency of each genome with a given
    CN in the population. Also, get the frequency of each combination
    of alleles for each CN.
    """
    combo_lists = []
    for cn in range(1, max_cn + 1):
        # Get individual CN's frequency in genome.
        cn_freq = sum([1 if len(x) == cn else 0 for x in genomes]) / float(len(genomes))
        # Get frequencies of allele combinations.
        wt, mix, mut = 0, 0, 0
        cn_genomes = [x for x in genomes if len(x) == cn]
        # Proportions of each genome type (out of genomes of that size).
        combo_list = count(cn_genomes, cn)
        # Adjust these proportions by the total proportion of genomes in the population
        # of the given copy number.
        combo_list = [x / float(len(genomes)) for x in combo_list]
        combo_lists.append((cn, combo_list))
    true_af = float(sum([sum(g) for g in genomes])) / sum([len(g) for g in genomes])
    return combo_lists, true_af

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bams", required=True, nargs='*', help='Path to sorted BAM file.')
    p.add_argument("--ref", required=True, help='Path to FASTA reference genome.')
    p.add_argument("-cn", type=int, help='Plot genomes with up to this many copies of the specified region. (default = 5)', default=5)
    p.add_argument("-o", help='Name of output plot.', default='stackedbar')
    p.add_argument("-rand", help="Plot expected random distribution of alleles", action='store_true')
    p.add_argument("-png", help='Output as PNG.', action='store_true')
    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
