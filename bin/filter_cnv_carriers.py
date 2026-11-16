#!/usr/bin/env python3
"""Set each CNV record's FILTER from its carriers' FORMAT/FT and, with --pass-only, drop the
records that don't pass. See modules/local/filter_cnv_carriers/main.nf for the CLI contract.

A carrier is a sample whose GT has an ALT allele (0/1, 1/1, ./1, haploid 1). FORMAT/FT holds each
sample's own filters: DRAGEN's, for the sample's segment, in a jointly-called VCF; its call's
FILTER on the per-sample route, copied by bin/add_sample_filter.awk. A record's FILTER becomes:

- PASS when at least one carrier's FT is PASS;
- otherwise the filters its carriers' FT name, e.g. cnvLength;cnvQual;
- CarrierFT when no carrier's FT names one (no carrier, or FT unset).

The FILTER this replaces doesn't say whether the CNV passed. A DRAGEN joint VCF's is PASS as soon
as any sample passes, even one without the CNV, whose FT then describes its normal segment, and
SampleFT otherwise. A proband VCF's (JOIN_PROBAND_CNV) is the proband's own call's, though a
parent's call in the same family record may pass. Exomiser (analysisMode PASS_ONLY) reads FILTER
only. Records are written in the input's order.
"""
import argparse
import sys

import pysam

FALLBACK = "CarrierFT"
FALLBACK_DESCRIPTION = "No carrier of the CNV has FORMAT/FT PASS, and none names a filter in it"
OTHER_DESCRIPTION = "Named in the FORMAT/FT of a carrier of the CNV"


def ft_codes(value):
    """FORMAT/FT -> its filter names, [] when unset."""
    if value is None:
        return []
    if isinstance(value, tuple):
        value = ";".join(v for v in value if v)
    return [code for code in value.split(";") if code and code != "."]


def is_carrier(sample):
    return any(allele is not None and allele > 0 for allele in (sample.get("GT") or ()))


def carrier_filter(rec):
    """The record's new FILTER, as a list of filter names."""
    codes = []
    for sample in rec.samples.values():
        if not is_carrier(sample):
            continue
        ft = ft_codes(sample.get("FT"))
        if ft == ["PASS"]:
            return ["PASS"]
        codes.extend(code for code in ft if code not in codes and code != "PASS")
    return codes or [FALLBACK]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf-in", required=True)
    parser.add_argument("--vcf-out", required=True, help="Output path, ending in .vcf.gz (a .tbi is written next to it)")
    parser.add_argument("--pass-only", action="store_true", help="Drop the records whose new FILTER isn't PASS")
    args = parser.parse_args()

    if not args.vcf_out.endswith(".vcf.gz"):
        sys.exit("filter_cnv_carriers.py: --vcf-out must end in .vcf.gz")

    # First pass: the filter names to declare.
    with pysam.VariantFile(args.vcf_in) as vcf:
        if "FT" not in vcf.header.formats:
            sys.exit(
                f"filter_cnv_carriers.py: {args.vcf_in} has no FORMAT/FT (DRAGEN's in a joint VCF, set by "
                "NORMALIZE_CNV's ext.sample_filter on the per-sample route): can't tell whether a carrier's call passed"
            )
        used = set()
        for rec in vcf:
            used.update(carrier_filter(rec))

    total = n_pass = kept = 0
    plain_out = args.vcf_out[: -len(".gz")]
    with pysam.VariantFile(args.vcf_in) as vcf:
        header = vcf.header.copy()
        for code in sorted(used - set(header.filters)):
            header.filters.add(code, None, None, FALLBACK_DESCRIPTION if code == FALLBACK else OTHER_DESCRIPTION)
        with pysam.VariantFile(plain_out, "w", header=header) as out:
            for rec in vcf:
                total += 1
                codes = carrier_filter(rec)
                n_pass += codes == ["PASS"]
                if args.pass_only and codes != ["PASS"]:
                    continue
                rec.translate(out.header)
                rec.filter.clear()
                for code in codes:
                    rec.filter.add(code)
                out.write(rec)
                kept += 1
    pysam.tabix_index(plain_out, preset="vcf", force=True)

    print(
        f"FILTER set from the carriers' FORMAT/FT: {n_pass} of {total} record(s) PASS. "
        + (f"Kept those {kept} (--pass-only)." if args.pass_only else "Kept all."),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
