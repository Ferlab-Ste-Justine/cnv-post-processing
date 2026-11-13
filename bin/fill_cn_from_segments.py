#!/usr/bin/env python3
"""Fill the missing (./.) genotypes of a family CNV VCF with each sample's own copy number, read from
its original per-sample DRAGEN VCF. See modules/local/fill_cn_from_segments/main.nf for the CLI
contract.

DRAGEN's per-sample CNV VCF covers the genome: REF blocks (ALT ".", CN 2, or 1 on a male's non-PAR
chrX/chrY) between the sample's CNV calls. NORMALIZE_CNV drops the REF blocks, so after the family
merge and collapse a sample with no call in a record is ./. with no FORMAT/CN, which slivar reads as
unknown. Here, for each such sample, every segment of its original VCF that overlaps the record (REF
blocks and its own calls, whatever their FILTER) is counted, by CN value, in bases of the record it
covers. When one CN covers at least --min-cover of the record, that CN is written:
  - the sample's expected CN (normal): GT becomes hom-ref, 0/0, or 0 on a male's non-PAR chrX/chrY
  - any other CN: FORMAT/CN only, GT stays ./. -- usually the sample's own call, which the collapse
    left in another record. slivar reads CN only; GT stays the sample's genotype for this record's
    ALT, which no call of its own supports, so FILTER_CNV_CARRIERS and Exomiser don't take it for a
    carrier.
The expected CN is slivar's get_expected_cn(), with the sex from the family's PED (expected_cn() in
refine_genotypes.py). Without one (chrM, a record straddling a PAR boundary, a non-PAR chrX/chrY
record and no sex), the CN is still written, and GT stays ./.

Existing genotypes are never changed, and no CNV is called here. Records shorter than --min-size
(END - POS) are left missing: under 5 kb the per-sample caller misses many losses the joint caller
finds, and a parent read as normal there turns the child's inherited CNV into `denovo`
(docs/journal.md).
"""
import argparse
import collections
import os
import sys

import pysam

# expected_cn() and load_ped_sexes() are depth refinement's, so the two fills and slivar's rules
# agree on the expected CN.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from refine_genotypes import expected_cn, load_ped_sexes  # noqa: E402


def load_segments(path):
    """Return (sample, {chrom: [(pos, end, cn)]}) from a single-sample VCF: every record with a CN,
    REF blocks and calls alike, as 1-based [POS, END]."""
    with pysam.VariantFile(path) as vcf:
        samples = list(vcf.header.samples)
        if len(samples) != 1:
            sys.exit(f"fill_cn_from_segments.py: {path} has {len(samples)} samples, expected the one of a per-sample VCF")
        if "CN" not in vcf.header.formats:
            sys.exit(f"fill_cn_from_segments.py: {path} declares no FORMAT/CN")
        segments = collections.defaultdict(set)  # a set: a split multiallelic site repeats its interval
        for rec in vcf:
            cn = rec.samples[0].get("CN")
            if cn is not None:
                segments[rec.chrom].add((rec.pos, rec.stop, cn))
    return samples[0], {chrom: sorted(segs) for chrom, segs in segments.items()}


def cn_cover(segments, chrom, pos, end):
    """{CN: bases of [pos, end] it covers} over the sample's segments."""
    cover = collections.Counter()
    for seg_pos, seg_end, cn in segments.get(chrom, ()):
        if seg_pos > end:
            break
        overlap = min(seg_end, end) - max(seg_pos, pos) + 1
        if overlap > 0:
            cover[cn] += overlap
    return cover


def covering_cn(segments, chrom, pos, end, min_cover):
    """The CN covering at least min_cover of [pos, end], or None (none does, or two tie)."""
    ranked = cn_cover(segments, chrom, pos, end).most_common(2)
    if not ranked or ranked[0][1] / (end - pos + 1) < min_cover:
        return None
    if len(ranked) == 2 and ranked[1][1] == ranked[0][1]:
        return None
    return ranked[0][0]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--vcf-in", required=True)
    parser.add_argument("--vcf-out", required=True)
    parser.add_argument("--min-size", type=int, default=0,
                        help="Leave records shorter than this (END - POS, bp) missing")
    parser.add_argument("--min-cover", type=float, required=True,
                        help="Share of a record one CN must cover in the sample's segments")
    parser.add_argument("--ped", help="PED whose sex column gives each VCF sample's sex (for chrX/chrY)")
    parser.add_argument("segment_vcfs", nargs="+", help="Each family member's original per-sample VCF")
    args = parser.parse_args()

    segments_by_sample = {}
    for path in args.segment_vcfs:
        sample, segments = load_segments(path)
        if sample in segments_by_sample:
            sys.exit(f"fill_cn_from_segments.py: two per-sample VCFs for sample {sample}")
        segments_by_sample[sample] = segments
    sexes = load_ped_sexes(args.ped) if args.ped else {}

    vcf_in = pysam.VariantFile(args.vcf_in)
    # A sample's segments are found by its VCF name, in the per-sample VCF's own header.
    no_segments = [name for name in vcf_in.header.samples if name not in segments_by_sample]
    if no_segments:
        sys.exit(
            f"fill_cn_from_segments.py: no per-sample VCF for sample(s) {', '.join(no_segments)} of {args.vcf_in} "
            f"(got: {', '.join(sorted(segments_by_sample))})"
        )
    if "CN" not in vcf_in.header.formats:
        sys.exit(f"fill_cn_from_segments.py: {args.vcf_in} declares no FORMAT/CN")
    vcf_out = pysam.VariantFile(args.vcf_out, "wz", header=vcf_in.header)

    filled = collections.Counter()
    for rec in vcf_in:
        if rec.stop - rec.pos < args.min_size:
            vcf_out.write(rec)
            continue
        for sample_name, sample_data in rec.samples.items():
            gt = sample_data.get("GT")
            if gt is not None and any(allele is not None for allele in gt):
                continue  # only ever fill genotypes that are fully missing (./.)

            cn = covering_cn(segments_by_sample[sample_name], rec.chrom, rec.pos, rec.stop, args.min_cover)
            if cn is None:
                continue
            sample_data["CN"] = cn
            normal_cn = expected_cn(rec.chrom, rec.pos, rec.stop, sexes.get(sample_name))
            if normal_cn is not None and cn == normal_cn:
                sample_data["GT"] = (0,) * normal_cn
                sample_data.phased = False
                filled["normal"] += 1
            else:
                filled["other"] += 1

        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    print(
        f"Filled {filled['normal'] + filled['other']} missing genotype(s) from the samples' own segments: "
        f"{filled['normal']} hom-ref with the expected CN, {filled['other']} with another CN (GT left missing).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
