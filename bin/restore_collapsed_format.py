#!/usr/bin/env python3
"""Give back their FORMAT to the samples truvari's cohort collapse left empty. See
modules/local/restore_collapsed_format/main.nf for the CLI contract.

When truvari collapse folds calls into the one it keeps, it copies an absorbed call's FORMAT into
each sample that's empty in the kept record, but only when that call's GT reads as set to
truvari: `truvari.get_gt()` gives GT.NON for DRAGEN's `./1` (a gain of unknown zygosity, CN 3+)
and GT.UNK for a haploid `1` (a male's chrX), the same as `./.`. Those absorbed calls are dropped
with their FORMAT (CN, OID...), and the sample then reads as having no call there.

truvari writes every absorbed call to its removed VCF, with INFO/MatchId set to the kept
record's INFO/CollapseId. For each kept record, this copies into every sample with no GT allele
set the FORMAT of that sample's absorbed call, the one with the highest TruScore when there are
several (truvari's own order: most similar first), and adds them to INFO/NumConsolidated.
A sample is restored from an absorbed call where that sample has at least one GT allele set.
"""
import argparse
import sys
from collections import defaultdict

import pysam


def has_call(sample):
    return any(allele is not None for allele in (sample.get("GT") or ()))


def first_value(value):
    return value[0] if isinstance(value, tuple) else value


def absorbed_by_kept_record(removed_path):
    """Return {kept CollapseId: [(TruScore, file order, {sample: (FORMAT dict, phased)})]}."""
    absorbed = defaultdict(list)
    with pysam.VariantFile(removed_path) as removed:
        for order, rec in enumerate(removed):
            match_id = first_value(rec.info.get("MatchId"))
            if match_id is None:
                continue
            calls = {
                name: ({key: sample[key] for key in rec.format}, sample.phased)
                for name, sample in rec.samples.items()
                if has_call(sample)
            }
            if calls:
                score = rec.info.get("TruScore")
                absorbed[match_id].append((score if score is not None else -1, order, calls))
    for candidates in absorbed.values():
        candidates.sort(key=lambda c: (-c[0], c[1]))
    return absorbed


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--collapsed-vcf", required=True, help="truvari collapse -o")
    parser.add_argument("--removed-vcf", required=True, help="truvari collapse -c")
    parser.add_argument("--vcf-out", required=True)
    args = parser.parse_args()

    absorbed = absorbed_by_kept_record(args.removed_vcf)
    restored = defaultdict(int)
    with pysam.VariantFile(args.collapsed_vcf) as collapsed, \
            pysam.VariantFile(args.vcf_out, "wz", header=collapsed.header) as out:
        if "NumConsolidated" not in collapsed.header.info:
            sys.exit(f"{args.collapsed_vcf} has no INFO/NumConsolidated: is it truvari collapse output?")
        for rec in collapsed:
            candidates = absorbed.get(first_value(rec.info.get("CollapseId")), [])
            n_restored = 0
            for name, sample in rec.samples.items():
                if has_call(sample):
                    continue
                source = next((calls[name] for _, _, calls in candidates if name in calls), None)
                if source is None:
                    continue
                fmt, phased = source
                for key, value in fmt.items():
                    sample[key] = value
                sample.phased = phased
                restored[name] += 1
                n_restored += 1
            if n_restored:
                rec.info["NumConsolidated"] = (rec.info.get("NumConsolidated") or 0) + n_restored
            out.write(rec)

    total = sum(restored.values())
    detail = ", ".join(f"{name} {n}" for name, n in sorted(restored.items()))
    print(f"Restored {total} absorbed call(s) truvari left out" + (f": {detail}." if detail else "."))


if __name__ == "__main__":
    main()
