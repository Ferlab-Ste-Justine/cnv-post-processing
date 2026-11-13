#!/usr/bin/awk -f
# Copy each record's FILTER into the sample's FORMAT/FT, for a single-sample VCF read from stdin
# (see modules/local/normalize_cnv/main.nf). A jointly-called DRAGEN VCF gives each sample its own
# filters in FORMAT/FT; this gives a per-sample VCF's calls the same field, which bcftools merge
# and truvari collapse carry into the family VCF like OID, so FILTER_CNV_CARRIERS reads one field
# on both routes. FILTER's ";"-separated filter names are valid in FT as they are.
BEGIN { FS = OFS = "\t" }

/^##FORMAT=<ID=FT,/ {
    print "add_sample_filter.awk: the VCF already has a FORMAT/FT" > "/dev/stderr"
    exit 1
}

/^##/ { print; next }

/^#CHROM/ {
    if (NF != 10) {
        print "add_sample_filter.awk: expected exactly one sample column, got " (NF - 9) > "/dev/stderr"
        exit 1
    }
    print "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter: the call's FILTER, 'PASS' when all filters passed\">"
    print
    next
}

{
    $9 = $9 ":FT"
    $10 = $10 ":" $7
    print
}
