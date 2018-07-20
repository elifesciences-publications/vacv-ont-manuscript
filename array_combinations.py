import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import pandas as pd
import scipy.stats as ss
import numpy as np
import random
from collections import defaultdict, OrderedDict
from extract_arrays import extract_arrays, get_refseq, create_random_arrays
import argparse
from itertools import product

def get_counts(arrays, cart_prods, cn=2):
    """
    function to count every possible combination
    of WT and His47Arg copies in a population of K3L arrays
    """
    all_combos = OrderedDict()
    for combo in cart_prods:
        count = sum([1 if combo == ''.join([str(x) for x in a]) else 0 for a in arrays])
        all_combos[combo] = count
    return all_combos

def run(args):
    if len(args.bams) == 1:
        print >> sys.stderr, "ERROR: please specify more than 1 BAM"
        sys.exit()
    # setting up the figure object
    sns.set_style('ticks')
    f, axarr = plt.subplots(len(args.bams), 1, figsize=(8,10))
    pal = sns.color_palette('Blues', len(args.bams))

    refseq = get_refseq(args.ref)        
    for i,bam in enumerate(args.bams):
       
        name = bam.split('/')[-1].split('.')[0].upper()

        arrays = extract_arrays(bam, refseq, copy_filter='hard').arrays
        af = extract_arrays(bam, refseq, copy_filter='hard').af
        # limit analysis to arrays of the specified copy number
        filtered_arrays = [tuple(a) for a in arrays if len(a) == args.cn]
        # determine all possible combinations of alleles in a array of `cn` copy
        # number using cartesian products
        cart_prods = [''.join([str(x) for x in g]) for g in product([0,1], repeat=args.cn)]
        cart_prods = sorted(cart_prods, key=lambda x: x.count('1'))
        # count every instance of these combinations in the sequenced data
        all_combos = get_counts(filtered_arrays, cart_prods, cn=args.cn)
        # count every array of the specified copy number
        total_cn_arrays = float(sum(all_combos.values()))
        # count the total number of arrays with mixed alleles
        total_mixed = float(sum(all_combos.values()[1:-1]))
        frac_mixed = total_mixed / total_cn_arrays

        x = range(len(cart_prods))
        # get the fraction of each allele combination in the sequence data
        y = [_ / total_cn_arrays for _ in all_combos.values()]

        axarr[i].plot(x, y, color=pal[i], marker='o', label='observed')
        axarr[i].text((len(x) - 1) / 4., (np.max(y) - np.min(y)) / 2. + 0.05, 
                'Mixed array fraction {}'.format(round(sum(y[1:-1]), 2)))
        axarr[i].axvline(1, color='r', ls=':')
        axarr[i].axvline(len(x) - 2, color='r', ls=':')

        axarr[i].tick_params(axis='y', labelsize=12.5, color='k')
        axarr[i].tick_params(axis='x', labelsize=12.5, color='k')
        axarr[i].spines['left'].set_color('k')
        axarr[i].spines['bottom'].set_color('k')


        axarr[i].legend()
        # figure/axis formatting
        if i == len(args.bams) - 1:
            axarr[i].set_xticks(x)
            axarr[i].set_xticklabels(['-'.join(list(c)) for c in cart_prods])
            axarr[i].set_xlabel('Allele combination (0 = $K3L^{WT}$, 1 = $K3L^{His47Arg}$)')
        else:
            axarr[i].get_xaxis().set_visible(False)
        axarr[i].set_title(name)
        for tick in axarr[i].get_xticklabels():
            tick.set_rotation(45)
        if i == 1:
            axarr[i].set_ylabel("Proportion of arrays")
        sns.despine(ax=axarr[i], trim=True)
    if args.png: 
        plt.savefig(args.o + '.png', format='png', bbox_inches='tight')
    else:
        plt.savefig(args.o + '.eps', format='eps', bbox_inches='tight')

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bams", required = True, help='Path to sorted BAM files.', nargs='*')
    p.add_argument("--ref", required = True, help='Path to FASTA reference genome.')
    p.add_argument("-cn", type=int, default=3, help='Plot arrays with this many copies of K3L. (default = 3)')
    p.add_argument("-o", help='Name of output plot.', default='array-combinations')
    p.add_argument('-png', help='Output as png.', action='store_true')
    run(p.parse_args(argv))

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
