#!/usr/bin/env python3
"""Join the family's FORMAT onto a proband's own CNV calls, through FORMAT/OID. See
modules/local/join_proband_cnv/main.nf for the CLI contract.

The proband VCF holds the proband's own calls: normalized (SVTYPE DEL/DUP), collapsed per sample
and annotated by VEP, so its CSQ describes each call's own interval. The family VCF (merged,
collapsed across the family, depth-refined) has one record per cluster of calls: truvari keeps
one call per cluster, often a parent's, but each sample's FORMAT/OID still names its own call
(bin/add_original_id.awk). For every proband call, this writes one record built from:

- the proband's call: CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO (with CSQ);
- the family record whose proband OID names that call: FORMAT and every sample's column;
- INFO/COHORT_ID, COHORT_POS and COHORT_END: that family record (COHORT_ID percent-encoded like
  OID, since ":" isn't allowed in an INFO value either), and COHORT_NumCollapsed /
  COHORT_NumConsolidated when it has them (the per-sample collapse writes keys of the same name).

The rest of the family record's INFO describes the call truvari kept, not the proband's, and is
left out. slivar runs on the output.

A family record is found by (OID, ALT): a split multiallelic site repeats its ID on each allele.
A proband call with no family record, or in two, is an error: the per-sample segments don't
overlap, so each call should be in exactly one family record. Family records naming a proband
call that isn't in the proband VCF only get a warning: VEP drops records on contigs its cache
doesn't have, and a split multiallelic family site repeats the proband's OID on its other allele.
Records are written in the proband VCF's order.
"""
import argparse
import re
import sys

import pysam

COHORT_INFO_HEADERS = [
    '##INFO=<ID=COHORT_ID,Number=1,Type=String,Description="ID of the family record this proband call was joined to (\':\' percent-encoded as %3A)">',
    '##INFO=<ID=COHORT_POS,Number=1,Type=Integer,Description="POS of the family record this proband call was joined to">',
    '##INFO=<ID=COHORT_END,Number=1,Type=Integer,Description="END of the family record this proband call was joined to">',
]

# Family INFO keys copied as COHORT_<key>, with their header line, when the family VCF has them.
COHORT_INFO_COPIED = {
    "NumCollapsed": '##INFO=<ID=COHORT_NumCollapsed,Number=1,Type=Integer,Description="Number of calls truvari collapsed into the family record this proband call was joined to">',
    "NumConsolidated": '##INFO=<ID=COHORT_NumConsolidated,Number=1,Type=Integer,Description="Number of samples truvari consolidated into the family record this proband call was joined to">',
}


def decode_oid(oid):
    """Undo add_original_id.awk's percent-encoding (":" -> %3A, "%" -> %25), in one pass."""
    return re.sub(r"%(25|3A)", lambda m: "%" if m.group(1) == "25" else ":", oid)


def encode_id(record_id):
    """add_original_id.awk's encoding, for an ID written as an INFO value (COHORT_ID)."""
    return record_id.replace("%", "%25").replace(":", "%3A")


def open_text(path):
    return pysam.BGZFile(path, "r") if path.endswith(".gz") else open(path, "rb")


def read_vcf(path):
    """Return (## header lines, #CHROM line, records as column lists)."""
    meta, chrom_line, records = [], None, []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.decode().rstrip("\n")
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                chrom_line = line
            elif line:
                records.append(line.split("\t"))
    if chrom_line is None:
        sys.exit(f"join_proband_cnv.py: no #CHROM line in {path}")
    return meta, chrom_line, records


def parse_info(info):
    """INFO column -> ordered {key: value}, value None for a flag."""
    if info == ".":
        return {}
    entries = {}
    for item in info.split(";"):
        key, sep, value = item.partition("=")
        entries[key] = value if sep else None
    return entries


def format_info(entries):
    if not entries:
        return "."
    return ";".join(key if value is None else f"{key}={value}" for key, value in entries.items())


def record_end(cols, info):
    if "END" in info:
        return info["END"]
    return str(int(cols[1]) + len(cols[3]) - 1)


def sample_oid(cols, col):
    """The sample's FORMAT/OID in this record, "." when unset. A sample's trailing FORMAT values
    can be left out, so a short column means "." too."""
    keys = cols[8].split(":")
    if "OID" not in keys or col >= len(cols):
        return "."
    values = cols[col].split(":")
    index = keys.index("OID")
    return values[index] if index < len(values) else "."


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--family-vcf", required=True, help="The family VCF, with FORMAT/OID (depth-refined)")
    parser.add_argument("--proband-vcf", required=True, help="The proband's own VEP-annotated per-sample VCF")
    parser.add_argument("--proband", required=True, help="The proband's sample name in both VCFs")
    parser.add_argument("--vcf-out", required=True, help="Output path, ending in .vcf.gz (a .tbi is written next to it)")
    args = parser.parse_args()

    if not args.vcf_out.endswith(".vcf.gz"):
        sys.exit("join_proband_cnv.py: --vcf-out must end in .vcf.gz")

    family_meta, family_chrom_line, family_records = read_vcf(args.family_vcf)
    proband_meta, proband_chrom_line, proband_records = read_vcf(args.proband_vcf)

    samples = family_chrom_line.split("\t")[9:]
    if args.proband not in samples:
        sys.exit(f"join_proband_cnv.py: proband '{args.proband}' isn't a sample of {args.family_vcf} ({', '.join(samples)})")
    proband_col = 9 + samples.index(args.proband)
    if not any(line.startswith("##FORMAT=<ID=OID,") for line in family_meta):
        sys.exit(f"join_proband_cnv.py: {args.family_vcf} has no FORMAT/OID (set by NORMALIZE_CNV's ext.original_id)")
    proband_samples = proband_chrom_line.split("\t")[9:]
    if proband_samples != [args.proband]:
        sys.exit(
            f"join_proband_cnv.py: {args.proband_vcf} should hold {args.proband}'s calls only, but its samples are "
            f"{', '.join(proband_samples) or 'none'}. Is it that sample's own VCF?"
        )

    # The family records by the proband call their OID names, with its ALT.
    family_by_call, no_call = {}, 0
    for cols in family_records:
        oid = sample_oid(cols, proband_col)
        if oid == ".":
            no_call += 1
            continue
        key = (decode_oid(oid), cols[4])
        if key in family_by_call:
            first = family_by_call[key]
            sys.exit(
                f"join_proband_cnv.py: {args.proband}'s call {key[0]} ({key[1]}) is in two records of {args.family_vcf}, "
                f"at {first[0]}:{first[1]} and {cols[0]}:{cols[1]}. truvari's cohort collapse should keep each call in one record."
            )
        family_by_call[key] = cols

    joined, seen = [], set()
    for cols in proband_records:
        key = (cols[2], cols[4])
        if key in seen:
            sys.exit(f"join_proband_cnv.py: call {key[0]} ({key[1]}) is on more than one record of {args.proband_vcf}")
        seen.add(key)
        family = family_by_call.get(key)
        if family is None:
            sys.exit(
                f"join_proband_cnv.py: {args.proband}'s call {key[0]} ({key[1]}) at {cols[0]}:{cols[1]} is in no record of "
                f"{args.family_vcf}: no family record has it as {args.proband}'s FORMAT/OID. Either {args.proband_vcf} isn't "
                f"{args.proband}'s VCF from this run, or a step after the per-sample collapse lost the call (truvari's "
                "cohort collapse, RESTORE_COLLAPSED_FORMAT)."
            )

        family_info = parse_info(family[7])
        info = parse_info(cols[7])
        info["COHORT_ID"] = encode_id(family[2])
        info["COHORT_POS"] = family[1]
        info["COHORT_END"] = record_end(family, family_info)
        for name in COHORT_INFO_COPIED:
            if name in family_info:
                info[f"COHORT_{name}"] = family_info[name]
        joined.append(cols[:7] + [format_info(info)] + family[8:])

    family_info_ids = {m.group(1) for m in (re.match(r"^##INFO=<ID=([^,>]+)", line) for line in family_meta) if m}
    out_meta = (
        [line for line in proband_meta if not line.startswith("##FORMAT=")]
        + [line for line in family_meta if line.startswith("##FORMAT=")]
        + COHORT_INFO_HEADERS
        + [header for name, header in COHORT_INFO_COPIED.items() if name in family_info_ids]
    )

    plain_out = args.vcf_out[: -len(".gz")]
    with open(plain_out, "w") as out:
        out.write("\n".join(out_meta + [family_chrom_line]) + "\n")
        for cols in joined:
            out.write("\t".join(cols) + "\n")
    pysam.tabix_index(plain_out, preset="vcf", force=True)

    unmatched = len(family_by_call) - len(joined)
    print(
        f"Joined {len(joined)} {args.proband} call(s) to the family's FORMAT; "
        f"{no_call} family record(s) without a {args.proband} call left out.",
        file=sys.stderr,
    )
    if unmatched:
        print(
            f"WARNING: {unmatched} family record(s) name a {args.proband} call (FORMAT/OID and ALT) that isn't in "
            f"{args.proband_vcf}, and are left out. VEP drops records on contigs its cache doesn't have, and a split "
            f"multiallelic family site repeats {args.proband}'s OID on its other allele.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
