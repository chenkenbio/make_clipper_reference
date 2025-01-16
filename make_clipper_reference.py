#!/usr/bin/env python3
"""
Author: Ken Chen (https://github.com/chenkenbio)
Date: 2025-01-15
"""

import argparse
import os
from tqdm import tqdm
import sys
import math
import random

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--gtf", required=True, help="GTF file")
    p.add_argument("-s", "--species", required=True, help="Species name")
    p.add_argument("-o", "--outdir", default=".", help="Output directory")
    # p.add_argument('--seed', type=int, default=2020)
    return p

def random_string(n: int=10) -> str:
    return "".join(random.choices("abcdefghijklmnopqrstuvwxyz", k=n))

def main():
    args = get_args().parse_args()

    TMPDIR = os.path.expanduser(os.environ.get("TMPDIR", "~/tmp"))

    tmp_genepred = os.path.join(TMPDIR, random_string(10) + ".genePred")

    # gtf2genepred
    cmd = f"gtfToGenePred -allErrors -genePredExt {args.gtf} {tmp_genepred}"
    if not os.path.exists(tmp_genepred):
        os.system(cmd)
    genes = dict()
    with open(tmp_genepred) as f:
        for line in tqdm(f):
            line = line.strip().split("\t")
            # genePredExt format
            chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_count, exon_starts, exon_ends = line[1:10]
            tx_start, tx_end = int(tx_start), int(tx_end)
            gene_id = line[11]
            exon_starts = [int(x) for x in exon_starts.rstrip(',').split(",") if x]
            exon_ends = [int(x) for x in exon_ends.rstrip(',').split(",") if x]
            # exons = list(zip(exon_starts.split(","), exon_ends.split(",")))
            exons = list(zip(exon_starts, exon_ends))
            if gene_id not in genes:
                genes[gene_id] = [chrom, tx_start, tx_end, strand, exons]
            else:
                assert genes[gene_id][0] == chrom
                assert genes[gene_id][3] == strand
                genes[gene_id][1] = min(genes[gene_id][1], tx_start)
                genes[gene_id][2] = max(genes[gene_id][2], tx_end)
                assert genes[gene_id][1] < genes[gene_id][2], f"{genes[gene_id][1]} < {genes[gene_id][2]} for {gene_id}"
                genes[gene_id][4].extend(exons)
    structure_fn = os.path.join(args.outdir, f"{args.species}.AS.STRUCTURE.COMPILED.gff")
    exon_fn = os.path.join(args.outdir, f"{args.species}_exons.bed")
    gene_fn = os.path.join(args.outdir, f"{args.species}_genes.bed")

    def _merge_intervals(intervals, to_length=False):
        intervals = sorted(intervals, key=lambda x: x[0])
        merged = []
        for i in intervals:
            if len(merged) == 0 or merged[-1][1] < i[0]:
                merged.append([i[0], i[1]])
            else:
                merged[-1][1] = max(merged[-1][1], i[1])
        if to_length:
            return sum([int(x[1]) - int(x[0]) for x in merged])
        return merged


    with open(structure_fn, "w") as structure_fp, open(exon_fn, "w") as exon_fp, open(gene_fn, "w") as gene_fp:
        for gene_id, (chrom, tx_start, tx_end, strand, exons) in genes.items():
            exons = sorted(exons, key=lambda x: (int(x[0]), int(x[1])))
            gene_fp.write(f"{chrom}\t{tx_start}\t{tx_end}\t{gene_id}\t0\t{strand}\n")
            premrna_length = int(tx_end) - int(tx_start)
            mrna_length = _merge_intervals(exons, to_length=True)
            # chr1    AS_STRUCTURE    gene    157784  157887  .       -       .       gene_id=ENSG00000222623.1;mrna_length=104;premrna_length=104
            structure_fp.write(f"{chrom}\tAS_STRUCTURE\tgene\t{tx_start + 1}\t{tx_end}\t.\t{strand}\t.\tgene_id={gene_id};mrna_length={mrna_length};premrna_length={premrna_length}\n")
            for i, (exon_start, exon_end) in enumerate(exons):
                exon_fp.write(f"{chrom}\t{exon_start}\t{exon_end}\t{gene_id}.{i}\t0\t{strand}\n")

    os.remove(tmp_genepred)


if __name__ == "__main__":
    args = get_args().parse_args()
    main()
