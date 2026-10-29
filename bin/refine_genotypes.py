#!/usr/bin/env python3
"""Resolve missing (./.) CNV genotypes to hom-ref where per-sample mosdepth-derived depth ratios
support it. See modules/local/depth_genotype_refine/main.nf for the CLI contract.

Ratio = (mean depth in the candidate CNV region) / (that sample's genome-wide mean depth), from
MOSDEPTH_RATIO. A ratio of (expected CN / 2) means "no depth deviation here", i.e. this sample
doesn't carry the CNV -- but only missing genotypes are ever touched; existing calls are never
overwritten and no new CNVs are called from depth alone.

A resolved genotype also gets FORMAT/CN set to the expected CN, since slivar reads CN, not GT:
a 0/0 with CN "." is still an unknown sample to it. The expected CN is the one slivar's rules use
(get_expected_cn() in assets/cnv-slivar-functions.js): 2 on autosomes and PAR; on non-PAR chrX 1
for a male and 2 for a female; on non-PAR chrY 1 for a male. chrM, records straddling a PAR
boundary, unknown sex on non-PAR chrX/chrY (no --ped, or sex 0) and a female on non-PAR chrY have
no expected CN and are left missing. The ratio is divided by (expected CN / 2) before the band
check, so a normal male's chrX (ratio ~0.5) falls in the band, and the GT written has that many
alleles (0 for a male's non-PAR chrX/chrY, 0/0 elsewhere).

Records shorter than --min-size (END - POS) are left missing: depth over a few kb is too noisy to
tell a normal sample from a one-copy loss the caller missed. On the full-genome NA12878 trio, 16
of the 28 `denovo` calls refinement produced were on parents' losses that DRAGEN's joint caller
found, read as normal by depth; 13 of those 16 under 5 kb (see docs/journal.md).
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


# GRCh38 pseudoautosomal regions, 1-based inclusive, as CNV_PAR_X / CNV_PAR_Y in
# assets/cnv-slivar-functions.js (keep in sync)
PAR = {
    "X": [(10001, 2781479), (155701383, 156030895)],
    "Y": [(10001, 2781479), (56887903, 57217415)],
}


def load_ped_sexes(path):
    """Return {sample: "male" | "female"} from a PED's sex column (1/2; anything else is unknown)."""
    sexes = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 5 or line.startswith("#"):
                continue
            sex = {"1": "male", "2": "female"}.get(fields[4])
            if sex:
                sexes[fields[1]] = sex
    return sexes


def expected_cn(chrom, start1, end1, sex):
    """slivar's get_expected_cn(): the sample's normal CN over [start1, end1], or None."""
    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    if chrom in ("M", "MT"):
        return None
    if chrom not in ("X", "Y"):
        return 2
    pars = PAR[chrom]
    if any(start1 >= lo and end1 <= hi for lo, hi in pars):
        return 2
    if any(start1 <= hi and end1 >= lo for lo, hi in pars):
        return None  # straddles a PAR boundary: no single normal CN
    if sex == "male":
        return 1
    if sex == "female" and chrom == "X":
        return 2
    return None


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
    parser.add_argument("--min-size", type=int, default=0,
                        help="Leave records shorter than this (END - POS, bp) missing")
    parser.add_argument("--ped", help="PED whose sex column gives each VCF sample's sex (for chrX/chrY)")
    parser.add_argument("ratio_beds", nargs="+", help="One or more <sample>.ratio.bed files")
    args = parser.parse_args()

    ratios_by_sample = {sample_name_from_path(p): load_ratio_bed(p) for p in args.ratio_beds}
    sexes = load_ped_sexes(args.ped) if args.ped else {}

    vcf_in = pysam.VariantFile(args.vcf_in)
    vcf_out = pysam.VariantFile(args.vcf_out, "wz", header=vcf_in.header)
    has_cn = "CN" in vcf_in.header.formats

    refined = 0
    for rec in vcf_in:
        # pysam treats INFO/END as a reserved attribute, exposed via rec.stop rather than
        # rec.info -- it round-trips to the output untouched as long as the record itself isn't
        # rebuilt from scratch, which we don't do here (only sample GT/CN values are modified).
        key = (rec.chrom, rec.pos - 1, rec.stop)
        # start1 as slivar's variant.start + 1, end1 as its variant.stop
        start1, end1 = rec.pos, rec.stop
        if rec.stop - rec.pos < args.min_size:
            vcf_out.write(rec)
            continue

        for sample_name, sample_data in rec.samples.items():
            gt = sample_data.get("GT")
            if gt is None or any(allele is not None for allele in gt):
                continue  # only ever touch genotypes that are fully missing (./.)

            normal_cn = expected_cn(rec.chrom, start1, end1, sexes.get(sample_name))
            ratio = ratios_by_sample.get(sample_name, {}).get(key)
            if normal_cn is None or ratio is None:
                continue
            if args.low <= ratio / (normal_cn / 2) <= args.high:
                sample_data["GT"] = (0,) * normal_cn
                sample_data.phased = False
                if has_cn:
                    sample_data["CN"] = normal_cn
                refined += 1

        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    print(f"Refined {refined} missing genotype(s) to hom-ref, with their expected CN.", file=sys.stderr)


if __name__ == "__main__":
    main()
