#!/usr/bin/env python3
"""Give back their FORMAT to the samples truvari's cohort collapse left empty, and a record of
their own to the calls it folded into a record already holding another call of their sample. See
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

A record has one column per sample, so when truvari folds two calls of one sample into the same
kept record, that column holds one of them and the other is in no record. With the cohort's
thresholds, that takes one call starting where the kept call does, at 70% of its length, and the
other starting where it ends (docs/journal.md). Each absorbed call of a sample (an ALT allele in
its GT) whose FORMAT/OID isn't that sample's in the kept record, once restored, is written after
that record as a record of its own, the way truvari writes a call it collapses with nothing: its
own site and INFO (minus truvari's match keys: TruScore, MatchId...), the FORMAT of each sample
it was lost for, the other samples missing (a call truvari copied into the kept record stays
only there).
"""
import argparse
import sys
from collections import defaultdict

import pysam


def has_call(sample):
    return any(allele is not None for allele in (sample.get("GT") or ()))


def is_alt_call(fmt):
    return 1 in (fmt.get("GT") or ())


def first_value(value):
    return value[0] if isinstance(value, tuple) else value


def absorbed_by_kept_record(removed_path):
    """Return {kept CollapseId: [(TruScore, file order, record, {sample: (FORMAT dict, phased)})]}."""
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
                absorbed[match_id].append((score if score is not None else -1, order, rec, calls))
    for candidates in absorbed.values():
        candidates.sort(key=lambda c: (-c[0], c[1]))
    return absorbed


def clear_sample(header, sample):
    """Set each FORMAT value of the sample missing, in its shape (GT keeps its ploidy). A string
    is written ".": pysam writes an unset one empty."""
    for key, value in list(sample.items()):
        if isinstance(value, tuple):
            sample[key] = (None,) * len(value)
        elif header.formats[key].type == "String":
            sample[key] = "."
        else:
            sample[key] = None


def own_record(header, absorbed, lost):
    """The absorbed record as truvari writes a call it collapses with nothing, in the output
    header: its own site and INFO, with only the samples in lost keeping their FORMAT."""
    new = absorbed.copy()
    # truvari's match keys are declared in the removed VCF's header only
    for key in list(new.info):
        if key not in header.info:
            del new.info[key]
    new.translate(header)
    for name, sample in new.samples.items():
        if name not in lost and has_call(sample):
            clear_sample(header, sample)
    return new


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--collapsed-vcf", required=True, help="truvari collapse -o")
    parser.add_argument("--removed-vcf", required=True, help="truvari collapse -c")
    parser.add_argument("--vcf-out", required=True)
    args = parser.parse_args()

    absorbed = absorbed_by_kept_record(args.removed_vcf)
    restored = defaultdict(int)
    own = []
    with pysam.VariantFile(args.collapsed_vcf) as collapsed, \
            pysam.VariantFile(args.vcf_out, "wz", header=collapsed.header) as out:
        if "NumConsolidated" not in collapsed.header.info:
            sys.exit(f"{args.collapsed_vcf} has no INFO/NumConsolidated: is it truvari collapse output?")
        if "OID" not in collapsed.header.formats:
            sys.exit(f"{args.collapsed_vcf} has no FORMAT/OID (set by NORMALIZE_CNV's ext.original_id): "
                     "can't tell which call a sample's column holds")
        for rec in collapsed:
            candidates = absorbed.get(first_value(rec.info.get("CollapseId")), [])
            n_restored = 0
            for name, sample in rec.samples.items():
                if has_call(sample):
                    continue
                source = next((calls[name] for _, _, _, calls in candidates if name in calls), None)
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

            for _, _, absorbed_rec, calls in sorted(candidates, key=lambda c: c[1]):
                lost = sorted(
                    name for name, (fmt, _) in calls.items()
                    if is_alt_call(fmt) and fmt.get("OID") != rec.samples[name].get("OID")
                )
                if lost:
                    out.write(own_record(out.header, absorbed_rec, lost))
                    own.extend(f"{name} {absorbed_rec.id}" for name in lost)

    total = sum(restored.values())
    detail = ", ".join(f"{name} {n}" for name, n in sorted(restored.items()))
    print(f"Restored {total} absorbed call(s) truvari left out" + (f": {detail}." if detail else "."))
    print(f"Gave {len(own)} absorbed call(s) a record of their own, the kept record holding another call of "
          "their sample" + (f": {', '.join(own)}." if own else "."))


if __name__ == "__main__":
    main()
