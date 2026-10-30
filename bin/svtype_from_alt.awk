#!/usr/bin/awk -f
# Set INFO/SVTYPE from a symbolic <DEL> or <DUP> ALT, for a VCF read from stdin whose multiallelic
# sites are already split (see modules/local/normalize_cnv/main.nf). DRAGEN writes SVTYPE=CNV on
# losses and gains alike, and truvari reads its variant type from SVTYPE before the ALT: CNV
# isn't one of its types, so every call was "UNK" and a loss could collapse with a gain, leaving
# the absorbed sample's GT pointing at the wrong ALT. Other ALTs and SVTYPEs are left as they are.
BEGIN { FS = OFS = "\t" }

/^#/ { print; next }

{
    if ($5 == "<DEL>") type = "DEL"
    else if ($5 == "<DUP>") type = "DUP"
    else { print; next }
    if ($8 ~ /(^|;)SVTYPE=[^;]*/) sub(/(^|;)SVTYPE=[^;]*/, (substr($8, 1, 7) == "SVTYPE=" ? "" : ";") "SVTYPE=" type, $8)
    else if ($8 == ".") $8 = "SVTYPE=" type
    else $8 = $8 ";SVTYPE=" type
    print
}
