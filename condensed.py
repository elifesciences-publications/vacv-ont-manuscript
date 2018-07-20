import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle
from extract_arrays import extract_arrays, get_refseq, create_random_arrays

def plot_rects(ax, x, i, h47r_p, wt_p, shift):
    """
    Plot linked rectangles (i.e., copies).
    """
    ax.add_patch(
        Rectangle(
        (x + shift, i),
            h47r_p,
            0.5,
        facecolor='darkblue',
        edgecolor='k',
        lw=1))

    ax.add_patch(
        Rectangle(
        (x + h47r_p + shift, i),
            wt_p,
            0.5,
        facecolor='w',
        edgecolor='k',
        lw=1))

def run(args):
   
    # Figure formatting.
    sns.set_style("white")
    f, ax = plt.subplots()
    sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
    ax.set_xlim([1, (args.cn * 1.5)])
    ax.set_ylim([0.5, (args.cn * 1.2)])
    ax.set(yticklabels=[], xticklabels=[])

    refseq = get_refseq(args.ref)        
    data = extract_arrays(args.bam, refseq, copy_filter='hard')
    arrays = data.arrays
    af = float(sum([sum(a) for a in arrays])) / sum([len(a) for a in arrays])
    # Create a population with randomly distributed H47R alleles if specified.
    if args.rand:
        arrays = create_random_arrays(arrays, max_cn=args.cn, af=af)

    for cn in range(1, args.cn + 1):
        matrix = [array for array in arrays if len(array) == cn]
        if args.mixed_only:
            if cn == 1: continue
            matrix = [array for array in matrix if len(set(array)) == 2]
        matrix = sorted(matrix)
        m = np.array(matrix)
        # Total number of arrays with the specified copy number.
        total_arrays = float(m.shape[0])
        # Among all arrays of the given copy number, how many H47R alleles 
        # are observed in the first copy, second copy, etc.? 
        h47r_copies = np.sum(m, axis=0)
        # Convert these numbers to proportions at each copy.
        h47r_proportions = [x / total_arrays for x in h47r_copies]
        copy_count = 1
        for h47r_p in h47r_proportions:
            shift = copy_count * 0.2
            wt_p = 1 - h47r_p
            plot_rects(ax, copy_count, cn, h47r_p, wt_p, shift)
            # Overly specific parameters to plot lines in between rectangles 
            # (i.e., copies) and summary text for each copy number.
            if copy_count < cn:
                plt.plot([copy_count + h47r_p + wt_p + shift,
                          copy_count + h47r_p + wt_p + shift + 0.2],
                          [cn + 0.25, cn + 0.25], '-', c='k', lw=0.75)
            if copy_count == cn:
                plt.text(copy_count + h47r_p + wt_p + shift + 0.2,
                         cn + 0.15, str(int(total_arrays)), style='italic')
            copy_count += 1

    if args.png: 
        plt.savefig(args.o + '.png', dpi=300)
    else:
        plt.savefig(args.o + '.eps')

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required = True, help='Path to sorted BAM file.')
    p.add_argument("--ref", required = True, help='Path to FASTA reference genome.')
    p.add_argument("-mixed_only", help="Only plot mixed arrays", action="store_true")
    p.add_argument("-rand", help='Plot random distribution of alleles', action='store_true')
    p.add_argument("-cn", type=int, help='Plot K3L arrays with up to this many copies of the specified region. (default = 5)', default=5)
    p.add_argument("-png", help='Output as PNG.', action='store_true')
    p.add_argument("-o", help='Name of output plot.', default='condensed')
    run(p.parse_args(argv))

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
