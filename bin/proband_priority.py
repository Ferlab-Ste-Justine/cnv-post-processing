#!/usr/bin/env python3
"""Make truvari collapse keep the proband's CNV as each cluster's representative.

truvari collapse has no option to prefer a sample; its `--keep maxqual` keeps the highest-QUAL
record of each cluster, and falls back to `--keep first` (leftmost POS) on ties. So:

  mark     Save each record's QUAL in INFO/ORIG_QUAL, then set QUAL to 1 on records where a
           proband carries the CNV (non-reference GT) and 0 on all others. Run truvari collapse
           with --keep maxqual on the result: the proband's call wins every cluster it's in, and
           clusters without one keep the leftmost call, exactly as with the default --keep first.
  restore  Put the original QUAL back from INFO/ORIG_QUAL and drop that field, in place, on
           truvari's collapsed and removed outputs.

See modules/local/truvari_collapse/main.nf for the CLI contract.
"""
import argparse
import os
import sys

import pysam

ORIG_QUAL = "ORIG_QUAL"


def carries(rec, sample):
    gt = rec.samples[sample].get("GT")
    return gt is not None and any(allele not in (None, 0) for allele in gt)


def mark(vcf_in_path, vcf_out_path, probands):
    vcf_in = pysam.VariantFile(vcf_in_path)
    vcf_samples = list(vcf_in.header.samples)
    present = [p for p in probands if p in vcf_samples]
    missing = [p for p in probands if p not in vcf_samples]
    if missing:
        print(f"Proband(s) not in the VCF, ignored: {', '.join(missing)}", file=sys.stderr)

    vcf_in.header.info.add(ORIG_QUAL, 1, "Float", "QUAL before proband prioritization (temporary)")
    vcf_out = pysam.VariantFile(vcf_out_path, "wz", header=vcf_in.header)

    marked = 0
    for rec in vcf_in:
        if rec.qual is not None:
            rec.info[ORIG_QUAL] = rec.qual
        if any(carries(rec, p) for p in present):
            rec.qual = 1
            marked += 1
        else:
            rec.qual = 0
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    pysam.tabix_index(vcf_out_path, preset="vcf", force=True)
    print(f"Marked {marked} proband record(s) for {', '.join(present) or 'no proband'}.", file=sys.stderr)


def restore(vcf_path):
    vcf_in = pysam.VariantFile(vcf_path)
    header = vcf_in.header.copy()
    if ORIG_QUAL in header.info:
        header.info.remove_header(ORIG_QUAL)
    tmp_path = f"{vcf_path}.restore.tmp.vcf.gz"
    vcf_out = pysam.VariantFile(tmp_path, "wz", header=header)

    for rec in vcf_in:
        rec.qual = rec.info.get(ORIG_QUAL)
        if ORIG_QUAL in rec.info:
            del rec.info[ORIG_QUAL]
        rec.translate(header)
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    os.replace(tmp_path, vcf_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_mark = sub.add_parser("mark", help="Rewrite QUAL so the proband's records rank first")
    p_mark.add_argument("--vcf-in", required=True)
    p_mark.add_argument("--vcf-out", required=True, help="bgzipped and tabix-indexed")
    p_mark.add_argument("--proband", action="append", required=True, help="VCF sample name; repeatable")

    p_restore = sub.add_parser("restore", help="Restore QUAL from INFO/ORIG_QUAL, in place")
    p_restore.add_argument("vcfs", nargs="+")

    args = parser.parse_args()
    if args.command == "mark":
        mark(args.vcf_in, args.vcf_out, args.proband)
    else:
        for vcf in args.vcfs:
            restore(vcf)


if __name__ == "__main__":
    main()
