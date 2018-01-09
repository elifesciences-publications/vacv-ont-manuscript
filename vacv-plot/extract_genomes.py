from collections import defaultdict, namedtuple
from pyfaidx import Fasta
import pysam
import random

def get_refseq(ref):
    refseq = None
    fa = Fasta(ref)
    for genome_id in fa.keys():
        genome = genome_id
        refseq = str(fa[genome_id])
    return refseq 

def site_is_polymorphic(refseq, pileupcolumn, af):
    """
    Quick and dirty method to determine if a site 
    within the user-inputted region is 'polymorphic' in the BAM.
    """
    ref_base = refseq[pileupcolumn.pos]
    alt_bases = []
    # Start counting the number of 'ref' pileup bases.
    # Any 'alt' bases will be added to the list above.
    read_ref_bases = 0
    for pileupread in pileupcolumn.pileups:
        if pileupread.alignment.mapping_quality < 20 or pileupread.query_position is None:
            continue
        # Store gaps with respect to the reference.
        if pileupread.is_del or pileupread.is_refskip:
            alt_bases.append('-')
            continue
        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
        if read_base == ref_base:
            read_ref_bases += 1
        elif read_base != ref_base:
            alt_bases.append(read_base)
    bases = ['A', 'T', 'C', 'G']
    # Count the frequency of each alternate base in the pileup.
    alt_count = [alt_bases.count(nuc) for nuc in bases]
    depth = pileupcolumn.n
    # Get the alternate base with the highest frequency in the pileup.
    max_nuc = bases[alt_count.index(max(alt_count))]
    # Check if the pileup site has an 'alt' allele frequency greater than
    # the user inputted limit, and is thus "polymorphic."
    if max(alt_count) / float(sum(alt_count) + read_ref_bases) > af:
        af = max(alt_count) / float(sum(alt_count) + read_ref_bases)
        return True, ref_base, max_nuc, af, depth
    else:
        return False, ref_base, max_nuc, None, depth

def get_clip(ct):
    """
    Ientify the number of bases that are soft or
    hard-clipped prior to a supplementary alignment of
    a query sequence (i.e., ONT read) using the pysam `cigar_tuple`.
    We use the number of bases clipped prior to an alignment as a proxy
    for that alignment's relative order in the original query, and therefore,
    its relative position in the viral genome.
    """
    MATCH = 0
    SOFT_CLIP = 4
    HARD_CLIP = 5

    clip = '-'
    for entry in ct:
        # Make sure we're looking at the first CIGAR entry in the CIGAR string.
        if ct.index(entry) == 0:
            if entry[0] in (SOFT_CLIP, HARD_CLIP): clip = entry[1]
            elif entry[0] == MATCH: clip = 0
        else: continue
    return clip

def create_coordinate_dict(bam, chromosome='gi|335317|gb|M35027.1|VACCG', 
                                left_coord=30363, right_coord=30629):
    """
    Catalog the start and end coordinates of each alignment
    within a query sequence (i.e., sequencing read). We use this
    to store the absolute start and end of each read with respect
    to the reference, for later filtering.
    """
    coords = defaultdict(list)
    samfile = pysam.Samfile(bam, "rb")
    for pileupread in samfile.fetch(chromosome, left_coord, right_coord):
        ref_end, ref_start = pileupread.reference_end, pileupread.reference_start
        q_name = str(pileupread.query_name)
        # Add the start and end coordinates of the alignment to the `coords` dictionary.
        coords[q_name].append(ref_start)
        coords[q_name].append(ref_end)
    return coords

def create_random_genomes(genomes, max_cn=5, af=0.5):
    """
    Method to generate a population of uniformly mutated
    viruses, such that the copy number distribution of the
    population is the same as the input population, and the
    allele frequency of the mutated population is the same
    as the input.
    """
    cn_freqs = []
    max_len = max([len(x) for x in genomes])
    for cn in range(1, max_len):
        # Get individual CN's frequency in genome
        cn_freq = sum([1 if len(x) == cn else 0 for x in genomes]) / float(len(genomes))
        cn_freqs.append(cn_freq)
    # Initialize population size, using observed frequencies of each
    # copy number in the sequence data.
    init = 500000
    pop_sizes = [int(round(init * x)) for x in cn_freqs]
    population = []
    # Initialize each genome in the population as a list 
    # of wild-type (0) genomes.
    for cn in range(1, max_len):
        lst = [[0] * cn for _ in range(pop_sizes[cn - 1])]
        population.extend(lst)
    # Set a limit on the total number of H47R alleles that could
    # be observed in the population, given our sequenced data.
    total_copies = sum([x * (pop_sizes.index(x) + 1) for x in pop_sizes])
    limit = int(round(total_copies * af))
    x = 0
    # Select random genomes and random sites within those genomes,
    # and mutate those sites to H47R, until we reach the observed
    # allele frequency in that population.
    while x < limit:
        i = random.randint(0, len(population) - 1)
        ii = random.randint(0, len(population[i]) - 1)
        probability = random.random()
        if probability < af and population[i][ii] == 0:
            population[i][ii] = 1
            x += 1
    return population 

def extract_genomes(bam, refseq, chromosome='gi|335317|gb|M35027.1|VACCG', 
                            left_coord=30300, right_coord=30800, gene_strand='-',
                            af=0.05, pos=30489, slop=150, allele_filter=None):
    """
    Identify the presence or absence of a variant site in every alignment 
    within all query sequences.
    """
    allele_dict = defaultdict(list)
    strand_dict = defaultdict(list)
    sa_tags = defaultdict(int)
    # Create a dictionary of the starts and ends of each alignment from each query sequence.
    coordinate_dict = create_coordinate_dict(bam, chromosome, left_coord, right_coord)
    samfile = pysam.Samfile(bam, "rb")
    # Iterate through the pileup columns within all reads that map to the locus.
    for pileupcolumn in samfile.pileup(chromosome, left_coord, right_coord):
        # Rather than examining every pileup column within the reads, only look at
        # columns within the portion of the reads that map to the locus.
        if pileupcolumn.pos != pos: continue
        if pileupcolumn.pos < left_coord: continue
        elif pileupcolumn.pos > right_coord: break
        # Get the ref/alt bases, as well as hacky estimates of allele frequency
        # and depth at the site.
        polymorphic, ref_base, alt_base, af, depth = site_is_polymorphic(refseq, pileupcolumn, af)
        # Iterate through each read that maps to the pileup column.
        for pileupread in pileupcolumn.pileups:
            mapq = pileupread.alignment.mapping_quality
            low_qual = False
            if mapq < 20: low_qual = True
            q_name = pileupread.alignment.query_name
            # Extract the number of supplementary alignments associated with
            # the query sequence, which we'll use as a proxy for copy number.
            try: sa_cn = len(pileupread.alignment.get_tag("SA").split(";"))
            except KeyError: sa_cn = 1
            sa_tags[q_name] = sa_cn
            del_skip, short = False, False
            # Get the soft/hard clip value from the 'left' end of the alignment.
            clip_value = get_clip(pileupread.alignment.cigartuples)
            # Catalog alignments that have deletions or gaps (with respect to the reference)
            # at the pileup site. We may want to filter them out in downstream analyses.
            if pileupread.is_del or pileupread.is_refskip:
                del_skip = True
            ref_end = pileupread.alignment.reference_end
            ref_start = pileupread.alignment.reference_start
            # Check that the individual alignment maps fully to at least the ends of the locus.
            # If the alignment doesn't map fully, add a placeholder entry to the `allele_dict`. 
            # As before, we need to catalog its existence (if we want to filter later).
            if ref_start > left_coord or ref_end < right_coord: 
                short = True
            q_pos = pileupread.query_position
            q_seq = pileupread.alignment.query_sequence
            if q_pos is not None:
                read_base = q_seq[q_pos]
            else: read_base = '-'
            # Get the strandedness of the alignment (and thus, the query).
            if pileupread.alignment.is_reverse: strand = '-'
            else: strand = "+"
            strand_dict[q_name].append(strand)
            # Double check that the alignment's parent query sequence is in our
            # dictionary of every query's coordinates. We'll be checking that the
            # terminal ends of the query map uniquely outside of the region later.
            # Allele calling.
            allele = '--'
            if read_base == ref_base: allele = 0
            elif read_base == alt_base: allele = 1
            # If the alignment has some "imperfection" (is missing a nucleotide at the mutated site,
            # doesn't align fully to the locus, low mapping quality, etc.), store it, 
            # but store its "allele" as a value  we can filter on later.
            if (del_skip and short) or low_qual: allele = '---'
            elif del_skip: allele = '--'
            elif short: allele = '-'
            # Check that the full query maps outside of the putative duplication.
            # If not, skip it. 
            if q_name not in coordinate_dict: continue
            if min(coordinate_dict[q_name]) > left_coord - slop or max(coordinate_dict[q_name]) < right_coord + slop:
                continue
            allele_dict[q_name].append((allele, clip_value, ref_base, alt_base))
    # Now, we go back through the dictionary of alleles, and sort each genome 5' to 3'.
    # Additionally, we remove any genomes with "imperfect" K3L alignments.
    # Going back through the dictionary is time-consuming, but in some cases
    # we want to catalog every allele in the genome of interest before we can apply filters.
    sorted_genomes = []
    for genome in allele_dict:
        # Sort each query's set of alignments by each alignment's clipping value.
        # A larger clipping value indicates that the alignment is "later" in a query.
        # Since K3L is on the negative strand, we want to reverse sort by clipping value.
        all_clips = [x[1] for x in allele_dict[genome]]
        contains_bad_clip = any([x == '-' for x in all_clips])
        if contains_bad_clip: continue
        ordered_genome = sorted(allele_dict[genome], key=lambda tup: tup[1], reverse=True)
        ordered_alleles = [x[0] for x in ordered_genome]
        # Check that we haven't catalogued more copies than there are SA.
        if len(ordered_alleles) > sa_tags[genome]: continue
        # Setting "hard" filter to be any K3L copy with an indel at position 30490,
        # a copy that doesn't align across the whole locus, or both.
        if allele_filter == 'hard':
            if any([i in ['-', '--', '---'] for i in ordered_alleles]): continue
        # The "soft" filter, used only to calculate gross copy number in the population,
        # includes alignments with indels at position 30490, since these alignments are 
        # otherwise indicative of a full K3L copy.
        elif allele_filter == 'soft' and any([i in ['-', '---'] for i in ordered_alleles]):
            continue
        sorted_genomes.append(ordered_alleles)

    results = namedtuple('Results', ['genomes', 'af', 'ref', 'alt', 'depth'])
    return results(sorted_genomes, af, ref_base, alt_base, depth)
