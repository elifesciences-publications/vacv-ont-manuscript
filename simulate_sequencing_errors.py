import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import numpy as np
import seaborn as sns
from extract_arrays import extract_arrays, get_refseq

def run(args):
    sns.set_style('ticks')

    # read in data
    refseq = get_refseq(args.ref)
    arrays = extract_arrays(args.bam, refseq, copy_filter='hard').arrays

    switch_dict = {1:0, 0:1}
    # dictionary of "switch probabilities" (i.e., error rates)
    # for T-C or C-T changes at K3L-His47Arg
    switch_probs = { 'r7':  {1: 0.015, 0: 0.023},
                     'r9':  {1: 0.024, 0: 0.023},
                     'r94': {1: 0.026, 0: 0.005} }
    # use the switch probabilities for the particular sequencing
    # chemistry used to make the input BAM
    switch_probs = switch_probs[args.chemistry]
    all_homogeneous = []
    true_mixed, true_homogeneous = 0, 0
    # loop over every array in the population, 
    # and convert mixed arrays into homogenous arrays
    for array in arrays:
        # single-copy arrays can't be "mixed"
        if len(array) < 2: continue
        # if the array is mixed, convert it to a homogenous array
        if len(set(array)) == 2:
            true_mixed += 1
            # if the array contains equal numbers of K3L-His47Arg
            # and K3L-WT copies, randomly decide whether it gets converted
            # into a homogenous WT or His47Arg array
            if array.count(1) - array.count(0) == 0:
                switch = random.random() < 0.5
                if switch:
                    all_homogeneous.append([0 for a in range(len(array))])
                else:
                    all_homogeneous.append([1 for a in range(len(array))])
            # else, assign the array to the homogeneity of its
            # predominant allele
            elif array.count(1) - array.count(0) > 0:
                all_homogeneous.append([1 for a in range(len(array))])
            else:
                all_homogeneous.append([0 for a in range(len(array))])
        # if the array is already homogeneous, keep it as such
        else:
            true_homogeneous += 1
            all_homogeneous.append(array)
    simulation_results = []
    sims = 1000
    for sim in range(sims):
        # loop over every one of the converted homogenous arrays,
        # and introduce "switches" (i.e., T-C or C-T) changes to 
        # mimic sequencing errors at the K3L-His47Arg site
        new_arrays = []
        for array in all_homogeneous:
            new_array = []
            for copy in array:
                if random.random() < switch_probs[copy]:
                    new_array.append(switch_dict[copy])
                else: new_array.append(copy)
            new_arrays.append(new_array)
        # record the numbers of mixed and homogenous arrays in the
        # randomly "switched" population of arrays
        sim_mixed, sim_homogeneous = 0, 0
        for array in new_arrays:
            if len(set(array)) == 2:
                sim_mixed += 1
            else: sim_homogeneous += 1
        # record the fraction of mixed arrays in the simulated population
        simulation_results.append(sim_mixed / float(sim_mixed + sim_homogeneous))
    
    sns.kdeplot(np.array(simulation_results), label="Expected fraction of mixed arrays")
    plt.axvline(true_mixed / float(true_mixed + true_homogeneous), color='r', 
                    label="Observed fraction of mixed arrays", ls='--')
    plt.ylabel("Kernel density")
    plt.xlabel("Fraction of mixed arrays in population")
    plt.legend()
    sns.despine(trim=True)
    if args.png:
        plt.savefig(args.o + '_' + args.chemistry + '.png', dpi=200, bbox_inches='tight')
    else:
        plt.savefig(args.o + '_' + args.chemistry + '.eps', bbox_inches='tight')
            
def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True, help='Path to sorted BAM file.')
    p.add_argument("--ref", required=True, help='Path to FASTA reference genome.')
    p.add_argument("--chemistry", required=True, help='Flowcell chemistry used to generate reads in BAM. One of "r94", "r9", or "r7."')
    p.add_argument("-o", help='Name of output plot.', default='sim-errors')
    p.add_argument("-png", help='Output as PNG.', action='store_true')
    run(p.parse_args(argv))

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
