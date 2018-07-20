from extract_arrays import extract_arrays, get_refseq

def run(args):
    bam, ref = args.bam, args.ref
    header = ["locus", "variant_site", "ref", "alt", "af", "coverage (# reads)", 
                "coverage (depth)", "cn_distrib"]
    print ('\t'.join(header))

    refseq = get_refseq(ref)        
    arrays, af, ref, alt, depth = extract_arrays(bam, refseq, copy_filter='soft')
    copy_num = []
    for n in range(args.cn):
        copy = "%s:%s" % (str(n + 1), str([len(a) for a in arrays].count(n + 1)))
        copy_num.append(copy)
    info = ['K3L', '30490', ref, alt, str(af), str(len(arrays)), str(depth)]
    print ('\t'.join(info)), '\t', (' | '.join(copy_num))

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required = True, help='Path to sorted BAM file.')
    p.add_argument("--ref", required = True, help='Path to FASTA reference genome.')
    p.add_argument('-cn', type=int, default=15, help='Print summary for arrays up to this copy number (default = 15)')
    run(p.parse_args(argv))

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
