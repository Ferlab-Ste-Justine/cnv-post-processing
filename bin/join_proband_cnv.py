#!/usr/bin/env python3
"""Join a proband's own CNV calls to the family's slivar VCF, through FORMAT/OID. See
modules/local/join_proband_cnv/main.nf for the CLI contract.

truvari's cohort collapse keeps one call per cluster, often a parent's, so a slivar record's ID,
POS, END, QUAL, FILTER and caller INFO describe that call, not the proband's. Each sample's
FORMAT/OID still names its own call (bin/add_original_id.awk), so for every slivar record where
the proband has a call, this writes one record built from:

- the proband's own call, from its caller VCF: CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO;
- the slivar record: the INFO keys the proband's call doesn't have (slivar's tags, VEP's CSQ,
  truvari's NumCollapsed...), and FORMAT with every family sample;
- INFO/COHORT_ID, COHORT_POS and COHORT_END: the slivar record it comes from (COHORT_ID
  percent-encoded like OID, since ":" isn't allowed in an INFO value either).

Records where the proband has no call (OID ".") are left out.
"""
import argparse
import re
import sys

import pysam

COHORT_INFO_HEADERS = [
    '##INFO=<ID=COHORT_ID,Number=1,Type=String,Description="ID of the family slivar record this proband call was joined from (\':\' percent-encoded as %3A)">',
    '##INFO=<ID=COHORT_POS,Number=1,Type=Integer,Description="POS of the family slivar record this proband call was joined from">',
    '##INFO=<ID=COHORT_END,Number=1,Type=Integer,Description="END of the family slivar record this proband call was joined from">',
]

# Header lines taken from the proband's caller VCF when the slivar VCF doesn't define them: its
# INFO, FILTER and ALT values are written as is. FORMAT comes from slivar, so its lines aren't.
MERGED_HEADER_TYPES = ("INFO", "FILTER", "ALT", "contig")

HEADER_ID = re.compile(r"^##(\w+)=<ID=([^,>]+)")


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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--slivar-vcf", required=True, help="The family's slivar VCF")
    parser.add_argument("--proband-vcf", required=True, help="The proband's own caller VCF (the samplesheet's vcf)")
    parser.add_argument("--proband", required=True, help="The proband's sample name in both VCFs")
    parser.add_argument("--vcf-out", required=True, help="Output path, ending in .vcf.gz (a .tbi is written next to it)")
    args = parser.parse_args()

    if not args.vcf_out.endswith(".vcf.gz"):
        sys.exit("join_proband_cnv.py: --vcf-out must end in .vcf.gz")

    slivar_meta, slivar_chrom_line, slivar_records = read_vcf(args.slivar_vcf)
    proband_meta, _, proband_records = read_vcf(args.proband_vcf)

    samples = slivar_chrom_line.split("\t")[9:]
    if args.proband not in samples:
        sys.exit(f"join_proband_cnv.py: proband '{args.proband}' isn't a sample of {args.slivar_vcf} ({', '.join(samples)})")
    proband_col = 9 + samples.index(args.proband)
    if not any(line.startswith("##FORMAT=<ID=OID,") for line in slivar_meta):
        sys.exit(f"join_proband_cnv.py: {args.slivar_vcf} has no FORMAT/OID (set by NORMALIZE_CNV's ext.original_id)")

    # The proband's calls by ID. ALT="." records are the caller's ref blocks, which NORMALIZE_CNV
    # drops, so no OID can point to them.
    calls = {}
    for cols in proband_records:
        if cols[4] == ".":
            continue
        if cols[2] in calls:
            sys.exit(f"join_proband_cnv.py: ID '{cols[2]}' is on more than one record of {args.proband_vcf}")
        calls[cols[2]] = cols

    defined = {m.groups() for m in map(HEADER_ID.match, slivar_meta) if m}
    extra_meta = []
    for line in proband_meta:
        m = HEADER_ID.match(line)
        if m and m.group(1) in MERGED_HEADER_TYPES and m.groups() not in defined:
            extra_meta.append(line)
            defined.add(m.groups())
    out_meta = slivar_meta + extra_meta + COHORT_INFO_HEADERS
    contig_order = {m.group(2): i for i, m in enumerate(filter(None, map(HEADER_ID.match, out_meta))) if m.group(1) == "contig"}

    joined, no_call = [], 0
    for cols in slivar_records:
        keys = cols[8].split(":")
        values = cols[proband_col].split(":")
        oid_index = keys.index("OID") if "OID" in keys else None
        # A sample's trailing FORMAT values can be left out, so a short column means OID ".".
        oid = values[oid_index] if oid_index is not None and oid_index < len(values) else "."
        if oid == ".":
            no_call += 1
            continue
        original_id = decode_oid(oid)
        call = calls.get(original_id)
        if call is None:
            sys.exit(
                f"join_proband_cnv.py: {args.proband}'s OID '{original_id}' at {cols[0]}:{cols[1]} isn't in {args.proband_vcf}. "
                "Is it that sample's own caller VCF?"
            )

        slivar_info = parse_info(cols[7])
        info = parse_info(call[7])
        for key, value in slivar_info.items():
            info.setdefault(key, value)
        info["COHORT_ID"] = encode_id(cols[2])
        info["COHORT_POS"] = cols[1]
        info["COHORT_END"] = record_end(cols, slivar_info)

        # The call's own ALT, except on a multiallelic call that NORMALIZE_CNV split: there the
        # slivar record holds the one allele it's about.
        alt = cols[4] if "," in call[4] else call[4]
        joined.append(call[:4] + [alt, call[5], call[6], format_info(info)] + cols[8:])

    # The call's own POS can come before an earlier slivar record's.
    joined.sort(key=lambda c: (contig_order.get(c[0], len(contig_order)), c[0], int(c[1])))

    plain_out = args.vcf_out[: -len(".gz")]
    with open(plain_out, "w") as out:
        out.write("\n".join(out_meta + [slivar_chrom_line]) + "\n")
        for cols in joined:
            out.write("\t".join(cols) + "\n")
    pysam.tabix_index(plain_out, preset="vcf", force=True)

    print(
        f"Joined {len(joined)} slivar record(s) to {args.proband}'s own calls; "
        f"{no_call} record(s) without a {args.proband} call left out.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
