#!/usr/bin/env python3
"""Resolve missing (./.) CNV genotypes to hom-ref (0/0) where per-sample mosdepth-derived
depth ratios support it. See modules/local/depth_genotype_refine/main.nf for the CLI contract.

Ratio = (mean depth in the candidate CNV region) / (that sample's genome-wide mean depth), from
MOSDEPTH_RATIO. A ratio close to 1.0 means "no depth deviation here", i.e. this sample doesn't
carry the CNV -- but only missing genotypes are ever touched; existing calls are never overwritten
and no new CNVs are called from depth alone.
"""
import argparse
import sys

import pysam


def load_ratio_bed(path):
    """Parse a "<sample>.ratio.bed" file into {(chrom, start, end): ratio}."""
    ratios = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom, start, end, ratio = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if ratio == "NA":
                continue
            try:
                ratios[(chrom, start, end)] = float(ratio)
            except ValueError:
                continue
    return ratios


def sample_name_from_path(path):
    # Matches MOSDEPTH_RATIO's output naming: "<sample>.ratio.bed"
    basename = path.rsplit("/", 1)[-1]
    suffix = ".ratio.bed"
    if not basename.endswith(suffix):
        raise ValueError(f"Expected a '*.ratio.bed' file, got: {path}")
    return basename[: -len(suffix)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-in", required=True)
    parser.add_argument("--vcf-out", required=True)
    parser.add_argument("--low", type=float, required=True, help="Lower bound of the hom-ref ratio band")
    parser.add_argument("--high", type=float, required=True, help="Upper bound of the hom-ref ratio band")
    parser.add_argument("ratio_beds", nargs="+", help="One or more <sample>.ratio.bed files")
    args = parser.parse_args()

    ratios_by_sample = {sample_name_from_path(p): load_ratio_bed(p) for p in args.ratio_beds}

    vcf_in = pysam.VariantFile(args.vcf_in)
    vcf_out = pysam.VariantFile(args.vcf_out, "wz", header=vcf_in.header)

    refined = 0
    for rec in vcf_in:
        # pysam treats INFO/END as a reserved attribute, exposed via rec.stop rather than
        # rec.info -- it round-trips to the output untouched as long as the record itself isn't
        # rebuilt from scratch, which we don't do here (only sample GT values are modified).
        key = (rec.chrom, rec.pos - 1, rec.stop)

        for sample_name, sample_data in rec.samples.items():
            gt = sample_data.get("GT")
            if gt is None or any(allele is not None for allele in gt):
                continue  # only ever touch genotypes that are fully missing (./.)

            ratio = ratios_by_sample.get(sample_name, {}).get(key)
            if ratio is not None and args.low <= ratio <= args.high:
                sample_data["GT"] = (0, 0)
                sample_data.phased = False
                refined += 1

        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    print(f"Refined {refined} missing genotype(s) to hom-ref.", file=sys.stderr)


if __name__ == "__main__":
    main()
