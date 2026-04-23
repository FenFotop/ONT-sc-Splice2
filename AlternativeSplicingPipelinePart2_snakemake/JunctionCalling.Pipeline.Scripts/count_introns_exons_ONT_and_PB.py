#Python script to grab all junctions from a bam file of aligned reads

# usage: 
# python count_introns_exons_ONT.py input.bam output_prefix

import pysam
import sys
import collections
from tqdm import tqdm

def robust_get_tag(read, tag="ts"):
    try: 
        return(read.get_tag(tag))
    except KeyError: 
        return(".")

BAM_CREF_SKIP = 3
match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
# Parse library_type argument and set tags accordingly
if len(sys.argv) > 3:
    library_type = sys.argv[3]
else:
    library_type = ""

if library_type == "ONT":
    cell_barcode_tag = "BC"
    umi_tag = "U8"
else:
    cell_barcode_tag = "CB"
    umi_tag = "XM"
library_end = sys.argv[4] if len(sys.argv) > 4 else "3prime"

def find_introns_single_cell(read_iterator, strandtag = "GS", cell_barcode_tag = cell_barcode_tag, umi_tag = umi_tag, total_reads=None):
    intron_counts = {}
    exon_counts = {}
    cell_barcodes = set()
    umis = collections.Counter()
        # Wrap iterator with tqdm for progress bar
    if total_reads is not None:
        iterator = tqdm(read_iterator, total=total_reads, desc="Processing reads")
    else:
        iterator = tqdm(read_iterator, desc="Processing reads")
    for r in iterator: # for each read
        base_position = r.pos # start position of read (probably not v trustworthy)
        chrom_strand = (r.reference_name, "-" if r.is_reverse else "+")
        cell_barcode = robust_get_tag(r, cell_barcode_tag)
        umi = robust_get_tag(r, umi_tag)
        umis[umi] += 1 # keep track of reads per umi
        if cell_barcode == ".": continue # unassigned read
        cell_barcodes.add(cell_barcode) # useful for outputting data later
        first = True
        exon_start = base_position
        for op, nt in r.cigartuples: # for each part of the CIGAR string
            if op in match_or_deletion: # exonic (or could be intron retention)
                base_position += nt
            elif op == BAM_CREF_SKIP:
                if not first: 
                    exon = (exon_start, base_position)
                    if not chrom_strand in exon_counts: 
                        exon_counts[chrom_strand] = {}
                    if not exon in exon_counts[chrom_strand]: 
                        exon_counts[chrom_strand][exon] = collections.Counter()
                    exon_counts[chrom_strand][exon][cell_barcode] += 1
                if first: # don't count first exon (don't trust start)
                    first = False
                junc_start = base_position
                base_position += nt
                if not chrom_strand in intron_counts: 
                    intron_counts[chrom_strand] = {}
                junc = (junc_start, base_position)
                if not junc in intron_counts[chrom_strand]: 
                    intron_counts[chrom_strand][junc] = collections.Counter()
                intron_counts[chrom_strand][junc][cell_barcode] += 1
                exon_start = base_position
    return(intron_counts,exon_counts,cell_barcodes,umis)
##MDS samples:
input_file_path = sys.argv[1]    
samfile = pysam.AlignmentFile(input_file_path, "rb")
# Get total number of reads for progress bar if possible
try:
    total_reads = samfile.mapped + samfile.unmapped
except AttributeError:
    total_reads = None
intron_counts,exon_counts,cell_barcodes,umis = find_introns_single_cell(samfile.fetch(), strandtag = "GS", total_reads=total_reads)
for chrom_strand,r in intron_counts.items(): 
    print(chrom_strand,len(r))

import gzip
import numpy as np
import os
from Bio.Seq import Seq

if library_type == "PB" and library_end == "3prime":
    cb_map = {str(Seq(cb).reverse_complement()): cb for cb in cell_barcodes}
else:
    cb_map = {cb: cb for cb in cell_barcodes}

cell_barcodes_rc = list(cb_map.keys())


def write_output(res_sc, output_file_path):
    with gzip.open(output_file_path, "wt", encoding='utf-8') as f:
        f.write("chrom\tstart\tend\tstrand\t" + "\t".join(cell_barcodes_rc) + "\n")
        for chrom_strand, r in res_sc.items():
            for junc, counts in r.items():
                row = "%s\t%i\t%i\t%s\t%s\n" % (
                    chrom_strand[0], junc[0], junc[1], chrom_strand[1],
                    "\t".join(str(counts[cb_map[cb_rc]]) for cb_rc in cell_barcodes_rc)
                )
                f.write(row)

write_output(exon_counts, sys.argv[2] + "_exons.tsv.gz")
write_output(intron_counts, sys.argv[2] + "_introns.tsv.gz")

umi_counts = list(umis.values())
print("Average value for the number of reads per umi:")
print(np.mean(umi_counts)) # knowles: 1.21 , Paulina: 1.611843008688349 					## Average value for the number of reads per umi's 
print("Median value for the number of reads per umi:")
print(np.median(umi_counts)) # knowles: 1.0, Paulina : 1.0   							## Median value for the number of reads per umi
print("Percent of umi's that have 2 supporting reads:")
print(np.sum(np.array(umi_counts)==2) / len(umi_counts)) # knowles: 0.139, Paulina: 0.16633508489002077  	## Percent of umi's that have 2 supporting reads
