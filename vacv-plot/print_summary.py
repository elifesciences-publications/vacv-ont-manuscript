from extract_genomes import extract_genomes, get_refseq

def run(args):
    bam, ref = args.bam, args.ref
    header = ["locus", "variant_site", "ref", "alt", "af", "coverage (# reads)", "coverage (depth)", "cn_distrib"]
    print ('\t'.join(header))

    refseq = get_refseq(ref)        
    genomes, af, ref, alt, depth = extract_genomes(bam, refseq, allele_filter='soft')
    copy_num = []
    for n in range(1,16):
        copy = "%s:%s" % (str(n), str([len(x) for x in genomes].count(n)))
        copy_num.append(copy)
    info = ['K3L', '30490', ref, alt, str(af), str(len(genomes)), str(depth)]
    print ('\t'.join(info)), '\t', (' | '.join(copy_num))

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required = True, help='Path to sorted BAM file.')
    p.add_argument("--ref", required = True, help='Path to FASTA reference genome.')
    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
