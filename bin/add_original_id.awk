#!/usr/bin/awk -f
# Copy each record's ID into the sample's FORMAT/OID, for a single-sample VCF read from stdin
# (see modules/local/normalize_cnv/main.nf). truvari collapse keeps one ID per cluster, but it
# copies FORMAT fields from every call it folds in, so OID still names each sample's own call in
# the cohort VCF. ":" separates FORMAT values, and DRAGEN IDs contain it
# (DRAGEN:LOSS:chr22:17287644-17294633), so OID percent-encodes it as %3A, as VCF 4.3 requires
# ("%" itself becomes %25 first).
BEGIN { FS = OFS = "\t" }

/^##/ { print; next }

/^#CHROM/ {
    if (NF != 10) {
        print "add_original_id.awk: expected exactly one sample column, got " (NF - 9) > "/dev/stderr"
        exit 1
    }
    print "##FORMAT=<ID=OID,Number=1,Type=String,Description=\"ID of this sample's own caller record, kept through the truvari collapses (':' percent-encoded as %3A)\">"
    print
    next
}

{
    oid = $3
    gsub(/%/, "%25", oid)
    gsub(/:/, "%3A", oid)
    $9 = $9 ":OID"
    $10 = $10 ":" oid
    print
}
