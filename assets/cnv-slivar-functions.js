// cnv-slivar-functions.js  (ES5 for slivar's duktape)
//
// Table-driven CNV mode-of-inheritance and parent-of-origin rules for SLIVAR_EXPR.
// The inheritance logic lives in one table you can read top-to-bottom:
//
//   TRIO_TABLE    one affected carrier child + its two parents  -> a "child call"
//
// Everything below the table is plumbing that turns copy numbers into the
// tokens the table uses. To change a rule, edit a table row, then re-run
// cnv-slivar-moi-harness.js and look at what changed in the TSV.
//
// The call describes one child and its parents only. What the rest of the
// family looks like (is the transmitting parent affected, do sibs carry the
// CNV) is left to whoever reads the tags, with the PED in hand.
//
// Call shape — like the SNV pipeline, mode of inheritance is a --family-expr per tag:
// has_moi_call(fam, variant, "<tag>") is true when at least one affected child of the
// family gets that call from cnv_moi(), so a family with several affected children can
// carry several MoI tags. slivar values a family tag with the family's affected members
// (or the family ID when nobody is affected), not with the child the call is about: the
// per-sample po_* tags (--sample-expr, valued with the children they resolve) say which.
//
// slivar expr --vcf in.vcf.gz --ped fam.ped --js cnv-slivar-functions.js \
//   --family-expr 'denovo:has_moi_call(fam, variant, "denovo")' \
//   --family-expr 'denovo_suspect:has_moi_call(fam, variant, "denovo_suspect")' \
//   --family-expr 'inherited_maternal:has_moi_call(fam, variant, "inherited_maternal")' \
//   --family-expr 'inherited_paternal:has_moi_call(fam, variant, "inherited_paternal")' \
//   --family-expr 'recessive:has_moi_call(fam, variant, "recessive")' \
//   --family-expr 'recessive_lowconf:has_moi_call(fam, variant, "recessive_lowconf")' \
//   --family-expr 'possible_denovo:has_moi_call(fam, variant, "possible_denovo")' \
//   --family-expr 'possible_inherited_maternal:has_moi_call(fam, variant, "possible_inherited_maternal")' \
//   --family-expr 'possible_inherited_paternal:has_moi_call(fam, variant, "possible_inherited_paternal")' \
//   --family-expr 'inheritance_unresolved:has_moi_call(fam, variant, "inheritance_unresolved")' \
//   --family-expr 'mendelian_inconsistent:has_moi_call(fam, variant, "mendelian_inconsistent")' \
//   --family-expr 'unknown_cn:has_moi_call(fam, variant, "unknown_cn")' \
//   --family-expr 'x_linked_denovo:has_moi_call(fam, variant, "x_linked_denovo")' \
//   --family-expr 'x_linked_recessive:has_moi_call(fam, variant, "x_linked_recessive")' \
//   --sample-expr 'po_mother:cn_po_origin(sample, variant) == "maternal"' \
//   --sample-expr 'po_father:cn_po_origin(sample, variant) == "paternal"' \
//   --sample-expr 'po_both:cn_po_origin(sample, variant) == "both"' \
//   --sample-expr 'po_ambiguous:cn_po_origin(sample, variant) == "ambiguous"' \
//   --sample-expr 'po_denovo:cn_po_origin(sample, variant) == "denovo"' \
//   --sample-expr 'po_possible_denovo:cn_po_origin(sample, variant) == "possible_denovo"' \
//   --sample-expr 'po_possible_mother:cn_po_origin(sample, variant) == "possible_maternal"' \
//   --sample-expr 'po_possible_father:cn_po_origin(sample, variant) == "possible_paternal"'
//
// Vocabulary
//   dev      = CN - expected CN. Expected CN is 2 on autosomes and PAR, 1 for a
//              male on non-PAR chrX / chrY. So dev -1 = one copy lost, +1 = one
//              copy gained, -2 = both copies lost, +2 = two copies gained.
//   carrier  = a sample with dev != 0.
//   carrier child = an AFFECTED carrier with at least one parent in the family:
//              the only samples cnv_moi() gives a call.
//   copy     = what a parent passes to the child: one of its two haplotypes (CN 3
//              is [2,1]: it passes 2 or 1), or a father's single X / Y.
//   hit      = a passed copy that changed the child's way: a deleted copy for a
//              loss, an extra copy for a gain (copy_effect()).
//   passes nothing = a parent who gives this child no copy of this chromosome:
//              a father's X to a son, any mother's Y. Such a parent counts as
//              normal, and their CN isn't needed (does_parent_pass_nothing()).
//   missing parent = not in the PED, or in the VCF without a usable CN: it could
//              have passed 0, 1 or 2 copies, like any parent with CN 0 to 3.
//   mendelian_inconsistent = no way for the parents to transmit the child's CN
//              (see is_mendelian_possible() below) -> sample swap, miscall, missed call.
//   inheritance_unresolved = transmission is possible, but who passed what can't be
//              told (e.g. kid 3, mom 3, dad 3; or a parent with CN 4) -> needs a human.
//
// ---------------------------------------------------------------------------
// Map: what happens when slivar calls has_moi_call(fam, variant, tag)
// Once per family and MoI tag (--family-expr). It asks cnv_moi() for each member's
// call and is true if any of them is <tag>; cnv_moi() itself works on one child.
// [name] = shared building block, see the box at the bottom.
// (cn_po_origin has its own map below: slivar calls it once per sample via --sample-expr.)
// ---------------------------------------------------------------------------
//
// has_moi_call(fam, variant, tag)
//   for each s in fam: cnv_moi(s, variant) == tag ? ─────────────► true
//   none ─────────────────────────────────────────────────────────► false
//
// cnv_moi(s, variant)
// │
// ├─① s.affected && ([has_parent(s,"mom")] || [has_parent(s,"dad")]) ?
// │        no (unaffected, or a founder) ─────────────────────────► null
// │
// ├─② is_modellable(s, variant)
// │     [get_expected_cn(s)] < 0, or < 0 for one of get_needed_parents(s) ?
// │        yes ────────────────────────────────────────────────────► null
// │     get_needed_parents(s): the parents in the family, minus one who
// │     passes nothing here (does_parent_pass_nothing)
// │
// ├─③ is_cn_missing(s, variant)
// │       ![has_cn(s)]                                          ─► "unknown_cn"
// │       [get_cn_dev(s)] == 0 ?          not a carrier: parents not needed
// │       ![has_cn(p)] for every one of get_needed_parents(s)   ─► "unknown_cn"
// │       (only one of them without CN: that parent counts as missing, like in a duo)
// │
// ├─④ [get_cn_dev(s)] == 0 ?  not a carrier ─────────────────────► null
// │
// └─⑤ trio_call(s, variant), then get_x_linked_call(): on non-PAR chrX, the X_LINKED_CALLS name
//       │
//       │  kd     = [get_cn_dev(kid)]                                 // How far is the kid from normal. -1 = 1 copy lost, 1 = 1 copy gained.
//       │  kidTok = populate_kid_tokens(kd)                 → ["any","±1",...] // returns array of strings that is then used
//       │  outcome = [get_transmission_outcome(kid)]       → "new" | "maternal" | "both" | "duo_maybe_maternal" | ...
//       │  row    = get_trio_row(kid): the first TRIO_TABLE row matching (kidTok, outcome),
//       │           the same row cn_po_origin() reads
//       │
//       │  the row's call ≠ "no_pattern" ─────────────────────► that call
//       │  the row's call = "no_pattern" (outcome "undecided") ─► is_mendelian_possible(kid)
//       │                             yes → "inheritance_unresolved" (these parents could have produced
//       │                                                             this child's CN, but who passed what
//       │                                                             can't be told)
//       │                             no  → "mendelian_inconsistent" (they couldn't have)
//       │
//       └ is_mendelian_possible(kid)
//            ├ [has_cn(kid)]   no → true (can't judge)
//            ├ fromMom = transmissions(kid,"mom")
//            ├ fromDad = transmissions(kid,"dad")
//            │    ├ [get_parent_transmission(kid, which)].haps when it can tell
//            │    ├ diploid CN >= 4 (or MENDEL_ALLOW_RARE_SPLITS) → range(0, p.CN)
//            │    └ parent missing or no [has_cn] → range(0, kid.CN)
//            └ some a ∈ fromMom, b ∈ fromDad with a + b = kid.CN ?
//
// ---------------------------------------------------------------------------
// Map: what happens when slivar calls cn_po_origin(sample, variant)
// Once per sample (--sample-expr): a transmission label for every child, SNV-style.
// Each po_* tag tests one result, and slivar writes the IDs of the samples it's true
// for: po_mother "maternal", po_father "paternal", po_both "both", po_ambiguous
// "ambiguous", po_denovo "denovo", po_possible_denovo "possible_denovo",
// po_possible_mother "possible_maternal", po_possible_father "possible_paternal" (the
// SNV pipeline's names). It reads the same TRIO_TABLE row as cnv_moi(), so the MoI
// call and the parent of origin always agree.
// ---------------------------------------------------------------------------
//
// cn_po_origin(s, variant)
// │
// ├─① exp = [get_expected_cn(s)]
// │     exp < 0  (chrM, PAR-straddling record, unknown sex on non-PAR X/Y,
// │              female on non-PAR Y) ────────────────────────────────► null
// │     ![has_cn(s)] ─────────────────────────────────────────────────► null
// │     s.CN == exp  (not a carrier) ─────────────────────────────────► null
// │
// └─② get_trio_row(s)[3]: the parent-of-origin column of the child's TRIO_TABLE row,
//     found from get_transmission_outcome(s, variant):
//       │
//       ├ m = get_parent_transmission(s, "mom")    what mom can pass, and a normal copy's value
//       ├ d = get_parent_transmission(s, "dad")    same for dad
//       │     ├ PAR-straddling record                      → null (not modelled). When a record crosses a PAR boundary, return null.
//       │     ├ non-PAR X: dad → son (or unknown sex)      → {haps: [0], normal: 0}  passes nothing
//       │     ├ non-PAR Y: dad → daughter, or mom → anyone → {haps: [0], normal: 0}  passes nothing
//       │     │    (these two hold even if that parent is absent / has no CN:
//       │     │     open question, ** see docs/journal.md) **
//       │     ├ parent missing, or no [has_cn]             → null
//       │     ├ non-PAR X/Y: the father's single copy      → {haps: [dad.CN], normal: 1}
//       │     │    (his X to a daughter, his Y to a son)
//       │     └ diploid: CN_HAPLOTYPES[CN], CN 0-3         → {haps: [a, b], normal: 1}
//       │                CN >= 4 (several splits fit)      → null
//       │
//       ├ both parents known (m and d): every way, mom copy + dad copy == s.CN (get_way_kinds),
//       │   named by which copies are hits (copy_effect)
//       │     no way, both can only pass normal copies ──────────► "new"
//       │     the ways disagree ─────────────────────────────────► "ambiguous"
//       │     only mom's copy a hit / only dad's / both ─────────► "maternal" / "paternal" / "both"
//       ├ one known, the other missing (get_duo_outcome): the missing parent passes 0, 1 or 2
//       │     the known parent only passes normal copies, or only
//       │     its normal copy fits ──────────────────────────────► "duo_new"
//       │     its hit fits, the missing parent's copy normal ────► "duo_maternal" / "duo_paternal"
//       │     its hit and the missing parent's hit ──────────────► "duo_maternal_two" / "duo_paternal_two"
//       │     its hit, or its normal copy and the missing hit ───► "duo_maybe_maternal" / "duo_maybe_paternal"
//       │     (SNV / Radiant: a het child with a het mother and an unknown father is
//       │      POSSIBLE_MOTHER; a hom-alt child with a het mother is MOTHER)
//       └ otherwise (CN >= 4, no ordinary way, copies changed in opposite directions) ► "undecided"
//
//   Examples (autosome unless noted; mom / dad copies in brackets; outcome → parent of origin):
//     kid 1, mom 1 [1,0], dad 2 [1,1]       only 0 + 1            maternal  → "maternal"   (harness T4)
//     kid 1, mom 1 [1,0], dad 1 [1,0]       0 + 1 or 1 + 0        ambiguous → "ambiguous"  (T10)
//     kid 0, mom 1 [1,0], dad 1 [1,0]       only 0 + 0            both      → "both"       (T6)
//     kid 0, mom 0 [0,0], dad 1 [1,0]       only 0 + 0            both      → "both"       (T27)
//     kid 1, mom 2 [1,1], dad 2 [1,1]       no way, both normal   new       → "denovo"     (T1)
//     kid 0, mom 2 [1,1], dad 1 [1,0]       no way, dad carrier   undecided → null         (T19)
//     duo: kid 1, mom 2 [1,1], no father    her normal copy       duo_new   → "possible_denovo"   (D3)
//     duo: kid 1, mom 1 [1,0], no father    her 0 + 1, or 1 + 0   duo_maybe_maternal → "possible_maternal" (D1)
//     duo: kid 0, mom 1 [1,0], no father    her 0 + his 0         duo_maternal_two   → "maternal" (D2)
//     duo: kid 1, mom 0 [0,0], no father    her 0 + his 1         duo_maternal       → "maternal" (D5)
//     non-PAR X: son 0, mom 1 [1,0]         dad passes nothing    maternal  → "maternal"   (X2)
//     non-PAR X: daughter 3, mom 2 [1,1], dad 2 [2]   1 + 2      paternal  → "paternal"   (X8)
//
//
// ┌─ shared building blocks ───────────────────────────────────────────────────────────────────┐
// │                                                                                            │
// │  get_cn_dev(s, variant)                            "how far from normal?"  → number | null │
// │    ├ get_expected_cn(s, variant)                   2 / 1 / -1 (can't model)                │
// │    │    ├ fix_cnv_chrom_name(variant)              "chr7" → "7"                            │
// │    │    ├ is_cnv_in_par(chrom, start, end)         uses CNV_PAR_X / CNV_PAR_Y              │
// │    │    ├ is_cnv_straddles_par(chrom, start, end)  crosses a PAR boundary → -1             │
// │    │    └ s.sex                                    on non-PAR X / Y                        │
// │    └ has_cn(s)                                     integer CN ≥ 0 ?                        │
// │                                                                                            │
// │  has_parent(s, "mom"|"dad")                        "mom" in s  (reading s.mom throws       │
// │                                                     in slivar when it isn't set)           │
// │                                                                                            │
// │  get_transmission_outcome(kid, variant)            "who passed what?" → TRIO_TABLE outcome │
// │    ├ get_parent_transmission(kid, "mom"|"dad")     copies it can pass:                     │
// │    │    diploid → CN_HAPLOTYPES[CN] (CN >= 4 → can't tell);                                │
// │    │    non-PAR X/Y father → [his CN] to a daughter (X) / son (Y),                         │
// │    │    [0] to the other sex; mother on Y → [0]                                            │
// │    ├ missing parent (get_duo_outcome) → 0, 1 or 2                                          │
// │    └ get_way_kinds(): every pair adding up to kid.CN, by its hits (copy_effect)            │
// │  get_trio_row(kid, variant): TRIO_TABLE row → the MoI call (trio_call) and                 │
// │    the parent of origin (cn_po_origin), so the two always agree.                           │
// └────────────────────────────────────────────────────────────────────────────────────────────┘

/**
 * A sample object, as slivar hands it to the JavaScript.
 *
 * Reading an attribute slivar didn't set THROWS ("unknown attribute:mom"), it
 * does not return undefined like a plain JS object: slivar then logs a
 * javascript error and writes the record with no tag. So never read s.mom /
 * s.dad without has_parent() first.
 *
 * @typedef {Object} Sample
 * @property {string}   id        sample ID from the PED
 * @property {boolean}  affected  phenotype from the PED
 * @property {string}   sex       "male" | "female" | "unknown"
 * @property {number}   CN        FORMAT/CN; -2147483648 when missing (".": htslib's int32 missing value)
 * @property {Sample}   [mom]     only set when the mother is in the PED and the VCF
 * @property {Sample}   [dad]     only set when the father is in the PED and the VCF
 * @property {Sample[]} kids      this sample's children in the family (may be empty)
 */

/**
 * A variant (VCF record) object, as slivar hands it to the JavaScript.
 * Same rule as samples: a missing INFO field throws, so the CNV end is read
 * from variant.stop (always set; slivar takes it from INFO/END when the
 * record has one), never from INFO.END directly.
 * @typedef {Object} Variant
 * @property {string} CHROM  "chr7", "7", "chrX", ...
 * @property {number} start  0-based start (POS - 1)
 * @property {number} stop   1-based end (INFO/END when present)
 */

// ===========================================================================
// TABLE 1 — TRIO_TABLE: one carrier child, and what its parents passed
// ===========================================================================
// One row gives both answers for a child: the MoI call (cnv_moi(), affected children
// only) and the parent of origin (cn_po_origin(), every child), so the two can't
// disagree. The row is the first one matching the child's dev and its transmission
// outcome (get_transmission_outcome(): every way one copy from each parent can add
// up to the child's CN).
// kid column (the child's dev):
//   "±1"    dev is -1 or +1            (one copy lost or gained)
//   "±2+"   |dev| >= 2                 (CN 0 or CN >= 4 on an autosome)
//   "even"  dev is -2, +2, -4, ...     (could be one copy from each parent)
//   "any"
// outcome column. A "hit" is a passed copy that changed the child's way: a deleted
// copy for a loss, an extra copy for a gain.
//   Both parents known (a parent who passes nothing counts as known and normal):
//     "new"          no way to make the child's CN, and both parents can only pass
//                    normal copies
//     "maternal"     every way: only the mother's copy is a hit
//     "paternal"     every way: only the father's copy is a hit
//     "both"         every way: both copies are hits
//     "ambiguous"    the ways disagree (e.g. kid 1, mom 1, dad 1: 0 + 1 or 1 + 0)
//   One parent known, the other missing (not in the PED, or no CN): the missing
//   parent could have passed 0, 1 or 2 copies, like any parent with CN 0 to 3.
//     "duo_new"            the known parent only passed a normal copy: the CNV is
//                          new, or from the missing parent
//     "duo_maternal"       the mother surely passed it; the father's copy was normal
//     "duo_maternal_two"   two hits: the mother's, and the father's
//     "duo_maybe_maternal" the mother's hit, or her normal copy and the father's hit
//     (and the same three for the father)
//   "undecided"    a parent with CN >= 4 (several splits), no way with ordinary
//                  copies, or copies changed in opposite directions: the catch-all
//                  row, where is_mendelian_possible() decides
// parent-of-origin column: what cn_po_origin() returns (null: no po_* tag).
// First matching row wins.
// prettier-ignore
var TRIO_TABLE = [
  // kid     outcome                 child call                     parent of origin      harness rows
  ["±1",     "new",                  "denovo",                      "denovo"],            // T1 T2 X3 X6 X7 X12 Y3 Y4, Q4 Q5
  ["±2+",    "new",                  "denovo_suspect",              "denovo"],            // T3 T26 (P5: two hits at one locus is implausible)
  ["any",    "maternal",             "inherited_maternal",          "maternal"],          // T4 T7 T18 T23 T24 X2 X9 X13, Q9 kid
  ["any",    "paternal",             "inherited_paternal",          "paternal"],          // T5 X8 Y2 Y5, Q3 Q10, Q9 sib
  ["even",   "both",                 "recessive",                   "both"],              // T6 T9 T27 X10, Q1 Q2 (T9 = CN4 from two CN3, see P7)
  ["any",    "both",                 "inheritance_unresolved",      "both"],
  ["any",    "ambiguous",            "inheritance_unresolved",      "ambiguous"],         // T8 T10
  ["any",    "duo_new",              "possible_denovo",             "possible_denovo"],   // D3 D7 T11
  ["any",    "duo_maternal",         "inherited_maternal",          "maternal"],          // D5
  ["any",    "duo_paternal",         "inherited_paternal",          "paternal"],          // X11
  ["even",   "duo_maternal_two",     "recessive_lowconf",           "maternal"],          // D2 D6
  ["even",   "duo_paternal_two",     "recessive_lowconf",           "paternal"],
  ["any",    "duo_maybe_maternal",   "possible_inherited_maternal", "possible_maternal"], // D1
  ["any",    "duo_maybe_paternal",   "possible_inherited_paternal", "possible_paternal"], // D4
  ["any",    "any",                  "no_pattern",                  null]                 // -> inheritance_unresolved (T16 T17) or mendelian_inconsistent (T14 T19 T25)
];

// On non-PAR chrX, a de novo call is x_linked_denovo, for a son or a daughter (the CN and
// the sex say the rest), and a daughter with both X copies hit is x_linked_recessive; other
// calls, PAR and the other chromosomes keep the names above (get_x_linked_call). The SNV
// pipeline's other X-linked names (Radiant's transmission modes) read dominant or recessive
// off the zygosity and who is affected: they don't fit a son's gain (CN 2: his normal X plus
// an extra copy, X12 X13), and the calls don't look at who is affected. So a child with one
// inherited X copy keeps inherited_*, a son (X2) like a daughter (X9). denovo_suspect and
// the duo calls keep their names too. See docs/journal.md.
// prettier-ignore
var X_LINKED_CALLS = {
  male:   { denovo: "x_linked_denovo" },                                 // X3 X12: his only X changed
  female: { denovo: "x_linked_denovo",                                   // X7
            recessive: "x_linked_recessive" }                            // X10: both of her X copies hit
};

// ===========================================================================
// The classifier — the order of the checks is the order you read here
// ===========================================================================
/**
 * Does at least one child of this family get this MoI call on this record?
 * Entry point for every MoI --family-expr (see the call shape at the top of
 * the file): slivar calls it once per family, record and tag. It only asks
 * cnv_moi() about each member, so a family with several affected children can
 * get several tags (Q9: inherited_maternal for one child, inherited_paternal
 * for the other), and slivar values each tag with the family's affected
 * members.
 * @param {Sample[]} fam      all samples of one family
 * @param {Variant}  variant  the CNV record
 * @param {string}   call     the tag being tested, e.g. "denovo"
 * @returns {boolean}
 * @example
 * // T4: affected kid CN1, unaffected mom CN1, dad CN2
 * has_moi_call(fam, variant, "inherited_maternal")  // → true  (the kid's call)
 * has_moi_call(fam, variant, "denovo")             // → false
 * @example
 * // Q9: mom 1, dad 3; affected kid CN1 (from mom), affected sib CN3 (from dad)
 * has_moi_call(fam, variant, "inherited_maternal")  // → true  (the kid)
 * has_moi_call(fam, variant, "inherited_paternal")  // → true  (the sib)
 * // T13: nobody affected
 * has_moi_call(fam, variant, "inherited_maternal")  // → false (no call for unaffected children)
 */
function has_moi_call(fam, variant, call) {
  for (var i = 0; i < fam.length; i++) if (cnv_moi(fam[i], variant) === call) return true;
  return false;
}

/**
 * Mode-of-inheritance call for one child on one CNV record: the classifier
 * behind every MoI tag. slivar doesn't call it directly; has_moi_call() asks
 * it about each family member (see the maps at the top of the file). A child
 * gets at most one call per record.
 *
 * Only affected children get a call, and the call looks at the child and its
 * parents only: whether the transmitting parent is affected, or what the sibs
 * carry, is left to whoever reads the tags.
 *
 * Checks, in order:
 *    1) an affected sample with at least one parent in the family;
 *    2) can the record be modelled for this child (is_modellable);
 *    3) is a CN we need missing (is_cn_missing) -> "unknown_cn";
 *    4) does the child carry the CNV (dev != 0);
 *    5) its TRIO_TABLE call (trio_call), with its X-linked name on non-PAR chrX
 *       (get_x_linked_call).
 *
 * for each CNV record:
 *    for each family:
 *       for each of the 14 MoI --family-expr:  ← slivar
 *          has_moi_call(fam, variant, tag)
 *             for each member s:
 *                cnv_moi(s, variant)           ← steps ①–⑤
 *
 * @param {Sample}  s        any sample; only affected children get a call
 * @param {Variant} variant  the CNV record
 * @returns {?string} one call, e.g. "denovo", "inherited_maternal", "unknown_cn",
 *   "inheritance_unresolved"; null when the sample gets no tag (unaffected,
 *   founder, not a carrier, unmodellable record)
 * @example
 * // affected kid CN1, unaffected mom CN1, dad CN2, chr7  (harness T4; see assets/moi_scenarios.md)
 * cnv_moi(kid, variant)  // → "inherited_maternal"
 * cnv_moi(mom, variant)  // → null (a founder)
 * @example
 * // Each example is one child on one record; CNs of affected kid / mom / dad, chr7 unless noted.
 * // kid 1, mom 2, dad 2                         both parents normal       (T1)
 * cnv_moi(kid, variant)  // → "denovo"
 * // kid 1, mom 2, AFFECTED dad 1                the call ignores who is affected  (T5)
 * cnv_moi(kid, variant)  // → "inherited_paternal"
 * // kid 0, mom 1, dad 1                         one deleted copy from each parent  (T6)
 * cnv_moi(kid, variant)  // → "recessive"
 * // kid 3, mom 3, dad 3                         inherited, but from which parent?  (T8)
 * cnv_moi(kid, variant)  // → "inheritance_unresolved"
 * // kid 0, mom 0, dad 2                         dad can't pass a deleted copy  (T14)
 * cnv_moi(kid, variant)  // → "mendelian_inconsistent"
 * // kid 1, mom 2, dad CN missing                read as a duo without a father  (T11)
 * cnv_moi(kid, variant)  // → "possible_denovo"
 * // duo: kid 1, mom 2, no father in the PED     (D3)
 * cnv_moi(kid, variant)  // → "possible_denovo"
 * // non-PAR chrX: son 0, carrier mom 1, dad 1   his only X is hers  (X2)
 * cnv_moi(son, variant)  // → "inherited_maternal"
 * // non-PAR chrX: daughter 1, mom 2, dad 1      de novo in a daughter  (X7)
 * cnv_moi(daughter, variant)  // → "x_linked_denovo"
 * // non-PAR chrY: son 2, mom (no Y CN), dad 2   his Y is his father's  (Y2)
 * cnv_moi(son, variant)  // → "inherited_paternal"
 * // singleton kid 1 (S1), or any child on chrM (M1)
 * cnv_moi(kid, variant)  // → null
 */
function cnv_moi(s, variant) {
  // for the T4 example, s is the kid:
  // {
  //   id:       "kid",
  //   affected: true,
  //   sex:      "female",
  //   CN:       1,
  //   kids:     [],          // she has no children of her own
  //   mom:      { id: "mom", affected: false, sex: "female", CN: 1, kids: [ <the kid> ] },
  //   dad:      { id: "dad", affected: false, sex: "male",   CN: 2, kids: [ <the kid> ] }
  // }
  //
  // variant:
  // {
  //   CHROM: "chr7",
  //   start: 1000000,                          // 0-based, i.e. POS − 1
  //   stop:  1200000,                          // 1-based END
  //   INFO:  { END: 1200000, SVTYPE: "DEL" }
  // }
  if (!s.affected || !(has_parent(s, "mom") || has_parent(s, "dad"))) return null; // affected children only
  if (!is_modellable(s, variant)) return null; // chrM, PAR-straddling record, unknown sex on X/Y, female on Y (P1)
  if (is_cn_missing(s, variant)) return "unknown_cn"; // a CN we need is missing (P15)
  if (get_cn_dev(s, variant) === 0) return null; // not a carrier
  return get_x_linked_call(s, variant, trio_call(s, variant));
}

/**
 * TRIO_TABLE call for one affected carrier child (step ⑤ of the map): the
 * child call (3rd column) of its row (get_trio_row), the row cn_po_origin()
 * reads too. The catch-all "no_pattern" row is resolved here with
 * is_mendelian_possible(). Exported for the harness / debugging.
 *
 * The call describes this child and its parents only (who is affected among
 * them is not looked at). A parent who passes nothing here counts as normal
 * (does_parent_pass_nothing).
 *
 * @param {Sample}  kid      affected carrier with at least one parent in the family
 * @param {Variant} variant
 * @returns {string} the child call, e.g. "denovo", "inherited_maternal",
 *   "recessive", "inheritance_unresolved", "mendelian_inconsistent"
 * @example
 * // kid CN1, mom CN1, dad CN2, chr7  (harness T4; called from cnv_moi(kid, variant))
 * //
 * // kid:
 * // {
 * //   id:       "kid",
 * //   affected: true,
 * //   sex:      "female",
 * //   CN:       1,
 * //   kids:     [],          // she has no children of her own
 * //   mom:      { id: "mom", affected: false, sex: "female", CN: 1, kids: [ <the kid> ] },
 * //   dad:      { id: "dad", affected: false, sex: "male",   CN: 2, kids: [ <the kid> ] }
 * // }
 * //
 * // variant:
 * // {
 * //   CHROM: "chr7",
 * //   start: 1000000,                          // 0-based, i.e. POS − 1
 * //   stop:  1200000,                          // 1-based END
 * //   INFO:  { END: 1200000, SVTYPE: "DEL" }
 * // }
 * //
 * // kid tokens ["any","±1"]; mom [1,0] + dad [1,1]: only 0 + 1, mom's copy is the hit: outcome "maternal"
 * trio_call(kid, variant)  // → "inherited_maternal"
 * @example
 * // kid / mom / dad CNs, chr7 unless noted; the outcome in brackets
 * // kid 1, mom 2, dad 2             [new]                                       (T1)
 * trio_call(kid, variant)  // → "denovo"
 * // kid 0, mom 2, dad 2             [new], two copies changed                   (T3)
 * trio_call(kid, variant)  // → "denovo_suspect"
 * // kid 0, mom 1, dad 1             [both]: only 0 + 0                          (T6)
 * trio_call(kid, variant)  // → "recessive"
 * // kid 0, mom 0, dad 1             [both]: only 0 + 0 (her [0,0] has no other copy)  (T27)
 * trio_call(kid, variant)  // → "recessive"
 * // kid 1, mom 1, dad 3             [maternal]: dad's CN3 passes 1 or 2, only 0 + 1 fits  (Q9)
 * trio_call(kid, variant)  // → "inherited_maternal"
 * // duo: kid 1, mom 1, no father    [duo_maybe_maternal]: her 0 + his 1, or her 1 + his 0  (D1)
 * trio_call(kid, variant)  // → "possible_inherited_maternal"
 * // duo: kid 1, mom 0, no father    [duo_maternal]: her 0 + his 1             (D5)
 * trio_call(kid, variant)  // → "inherited_maternal"
 * // kid 3, mom 3, dad 3             [ambiguous]: 2 + 1 or 1 + 2                 (T8)
 * trio_call(kid, variant)  // → "inheritance_unresolved"
 * // kid 0, mom 0, dad 2             [undecided]: no way, not transmittable      (T14)
 * trio_call(kid, variant)  // → "mendelian_inconsistent"
 * // non-PAR chrX: son 0, mom 1, dad 1   [maternal]: dad passes no X to a son  (X2)
 * trio_call(kid, variant)  // → "inherited_maternal"
 * // non-PAR chrY: son 2, mom (no Y CN), dad 1   [new]: no Y from a mother  (Y3)
 * trio_call(kid, variant)  // → "denovo"
 */
function trio_call(kid, variant) {
  var row = get_trio_row(kid, variant); // e.g. T4: ["any", "maternal", "inherited_maternal", "maternal"]
  // Return the value of row[2], i.e. the 3rd value of the lookup table (the child call).
  //      Except if row[2] is "no_pattern" (the catch-all row: the outcome is "undecided"). If that is the
  //      case, we'll go with the is_mendelian_possible() which assess if:
  //          "mendelian_inconsistent" (none of our patterns fit this child. Could these parents have produced
  //                                    this child's CN at all -> No)
  //          "inheritance_unresolved" (none of our patterns fit this child. Could these parents have produced
  //                                    this child's CN at all -> Yes)
  if (row[2] !== "no_pattern") return row[2];
  return is_mendelian_possible(kid, variant) ? "inheritance_unresolved" : "mendelian_inconsistent";
}

/**
 * The TRIO_TABLE row of one carrier child: the first row whose kid column the
 * child's dev satisfies and whose outcome column is the child's transmission
 * outcome (get_transmission_outcome). trio_call() reads its child call (3rd
 * column) and cn_po_origin() its parent of origin (4th column), so the two
 * always agree. Exported for the harness / debugging.
 * @param {Sample}  kid      a carrier child
 * @param {Variant} variant
 * @returns {Array} a TRIO_TABLE row; the catch-all row when nothing more specific fits
 * @example
 * // kid CN1, mom CN1, dad CN2 (T4): kid tokens ["any", "±1"], outcome "maternal"
 * get_trio_row(kid, variant)  // → ["any", "maternal", "inherited_maternal", "maternal"]
 * // kid 3, mom 4, dad 3 (T16): outcome "undecided" (mom's CN4 can't be split)
 * get_trio_row(kid, variant)  // → ["any", "any", "no_pattern", null]
 */
function get_trio_row(kid, variant) {
  var kidTok = populate_kid_tokens(get_cn_dev(kid, variant)); // kid gets an array of strings. Ex could be ["any", "±2+", "even"] or ["any", "±1"] . Always contains "any" at first pos. "any" can be seen as a placeholder that triggers when we found no matches in the TRIO_TABLE.
  var outcome = get_transmission_outcome(kid, variant); // one single string value. ex: "new", "maternal", "both", "duo_maybe_maternal" or "undecided"
  for (var i = 0; i < TRIO_TABLE.length; i++) {
    var row = TRIO_TABLE[i];
    // excerpt from the TRIO table.
    // kid      outcome      child call             parent of origin   harness rows
    // ["±1",   "new",       "denovo",              "denovo"],         // T1 T2 X3 X6 X7 X12 Y3 Y4, Q4 Q5
    // ["±2+",  "new",       "denovo_suspect",      "denovo"],         // T3 T26 (P5: two hits at one locus is implausible)
    // ["any",  "maternal",  "inherited_maternal",  "maternal"],       // T4 T7 T18 T23 T24 X2 X9 X13, Q9 kid
    if (
      // ex: is "±1" found in ["any", "±1"]
      kidTok.indexOf(row[0]) >= 0 && // is the row's kid value in the kidTok's list? again kid's list could be  ["any", "±1"]. yes ->
      (row[1] === "any" || row[1] === outcome) // row accepts any outcome, or exactly this one? If both yes, this is the child's row.
    )
      return row;
  }
}

/**
 * The name a child's call takes on non-PAR chrX (X_LINKED_CALLS), by the
 * child's sex; everywhere else, and for calls the table doesn't list, the call
 * itself. Only reached for modellable records, so the child's sex is known.
 * @param {Sample}  kid
 * @param {Variant} variant
 * @param {string}  call     trio_call()'s result
 * @returns {string}
 * @example
 * get_x_linked_call(<son>, <non-PAR chrX>, "denovo")                 // → "x_linked_denovo"
 * get_x_linked_call(<daughter>, <non-PAR chrX>, "denovo")            // → "x_linked_denovo"
 * get_x_linked_call(<daughter>, <non-PAR chrX>, "recessive")         // → "x_linked_recessive"
 * get_x_linked_call(<son>, <non-PAR chrX>, "inherited_maternal")      // → "inherited_maternal" (unchanged)
 * get_x_linked_call(<son>, <chrX PAR1>, "denovo")                    // → "denovo" (PAR: autosomal names)
 * get_x_linked_call(<kid>, <chr7>, "recessive")                       // → "recessive"
 */
function get_x_linked_call(kid, variant, call) {
  if (fix_cnv_chrom_name(variant) !== "X" || is_cnv_in_par("X", variant.start + 1, variant.stop)) return call;
  var names = X_LINKED_CALLS[kid.sex];
  return names && names[call] ? names[call] : call;
}

// ---------------------------------------------------------------------------
// Tokens for the table
// ---------------------------------------------------------------------------
/**
 * Every TRIO_TABLE kid-column value that a child's dev satisfies. A row's kid
 * column matches when its value is in this list, so "any" always matches.
 * "±1" and "±2+" are exclusive; "even" comes on top of "±2+" for an even dev.
 * @param {number} d  the child's dev (CN - expected CN), never 0 here
 * @returns {string[]} always starts with "any"; then "±1", "±2+", "even" as they apply
 * @example
 * populate_kid_tokens(-1)  // → ["any", "±1"]
 * populate_kid_tokens(-2)  // → ["any", "±2+", "even"]
 * @example
 * populate_kid_tokens(+1)  // → ["any", "±1"]           (CN 3 on an autosome)
 * populate_kid_tokens(+2)  // → ["any", "±2+", "even"]  (CN 4: could be one extra copy from each parent)
 * populate_kid_tokens(+3)  // → ["any", "±2+"]          (CN 5: odd, so not "even")
 */
function populate_kid_tokens(d) {
  // every kid-column value this dev satisfies
  var t = ["any"];
  if (Math.abs(d) === 1) t.push("±1");
  if (Math.abs(d) >= 2) t.push("±2+");
  if (d !== 0 && d % 2 === 0) t.push("even");
  return t;
}
// ---------------------------------------------------------------------------
// Child helpers
// ---------------------------------------------------------------------------
/**
 * Is a CN we need for this child missing? The child's own CN always; its
 * parents' only when the child carries the CNV, and only the parents whose
 * copies matter here (get_needed_parents): a missing CN on a parent who passes
 * nothing, or on anyone else in the family, doesn't block the call (P15).
 * A needed parent without CN only blocks it when no needed parent has one:
 * otherwise that parent counts as missing (get_transmission_outcome), as for a
 * duo, the way the SNV pipeline and Radiant treat a parent's no-call.
 * @param {Sample}  kid      an affected child, with at least one parent in the family
 * @param {Variant} variant
 * @returns {boolean} true if the child has no usable CN, or is a carrier none of
 *   whose needed parents has a usable CN
 * @example
 * // kid CN1, mom CN2, dad CN missing (T11): dad read as absent
 * is_cn_missing(kid, variant)  // → false
 * // kid CN1, mom CN missing, dad CN missing
 * is_cn_missing(kid, variant)  // → true
 * @example
 * // kid CN missing, parents CN2 (T20)
 * is_cn_missing(kid, variant)  // → true
 * // kid CN2 (not a carrier), mom CN2, dad CN missing
 * is_cn_missing(kid, variant)  // → false (a non-carrier's parents aren't needed)
 * // non-PAR chrY: son CN2, mother (no CN on Y), dad CN2 (Y2)
 * is_cn_missing(son, variant)  // → false (a mother passes no Y: her CN isn't needed)
 */
function is_cn_missing(kid, variant) {
  if (!has_cn(kid)) return true;
  if (get_cn_dev(kid, variant) === 0) return false;
  var parents = get_needed_parents(kid, variant);
  for (var i = 0; i < parents.length; i++) if (has_cn(parents[i])) return false;
  return parents.length > 0;
}

/**
 * Can we say what the normal CN is on this record, for this child and for every
 * parent whose copies matter (get_needed_parents)? If not, the child gets no tag.
 * Not modellable: chrM; a chrX/chrY record straddling a PAR boundary; non-PAR
 * chrX/chrY with a child (or a needed parent) of unknown sex; a daughter on
 * non-PAR chrY.
 * @param {Sample}  kid      an affected child, with at least one parent in the family
 * @param {Variant} variant
 * @returns {boolean} false if get_expected_cn() is -1 for the child or a needed parent
 * @example
 * is_modellable(<kid of a trio>, <chr7 record>)                        // → true
 * is_modellable(<kid of a trio>, <chrM record>)                        // → false
 * is_modellable(<daughter, father's sex unknown>, <non-PAR chrX>)      // → false
 * is_modellable(<son, father's sex unknown>, <non-PAR chrX>)           // → true (he passes no X to a son)
 * is_modellable(<daughter, father's sex unknown>, <chrX PAR1>)         // → true (2 copies whatever the sex)
 * is_modellable(<daughter>, <non-PAR chrY>)                            // → false (no Y)
 * is_modellable(<son>, <non-PAR chrY>)                                 // → true (the mother passes no Y)
 * is_modellable(<kid of a trio>, <chrX record crossing the PAR1 end>)  // → false
 */
function is_modellable(kid, variant) {
  if (get_expected_cn(kid, variant) < 0) return false;
  var parents = get_needed_parents(kid, variant);
  for (var i = 0; i < parents.length; i++) if (get_expected_cn(parents[i], variant) < 0) return false;
  return true;
}

/**
 * The parents whose copies matter for this child on this record: those in the
 * family, minus a parent who passes nothing here (does_parent_pass_nothing).
 * @param {Sample}  kid
 * @param {Variant} variant
 * @returns {Sample[]} zero, one or two of kid.mom / kid.dad
 * @example
 * get_needed_parents(<kid of a trio>, <chr7>)          // → [mom, dad]
 * get_needed_parents(<duo kid, mother only>, <chr7>)   // → [mom]
 * get_needed_parents(<son>, <non-PAR chrX>)            // → [mom] (his X is hers)
 * get_needed_parents(<daughter>, <non-PAR chrX>)       // → [mom, dad]
 * get_needed_parents(<son>, <non-PAR chrY>)            // → [dad] (his Y is his father's)
 */
function get_needed_parents(kid, variant) {
  var out = [];
  if (has_parent(kid, "mom") && !does_parent_pass_nothing(kid, "mom", variant)) out.push(kid.mom);
  if (has_parent(kid, "dad") && !does_parent_pass_nothing(kid, "dad", variant)) out.push(kid.dad);
  return out;
}

/**
 * Does this parent give this child no copy at all of this chromosome? A father
 * passes no X to a son and a mother passes no Y, so on non-PAR chrX / chrY that
 * parent counts as passing nothing, whatever their own CN and even when they
 * aren't in the PED (open question, see docs/journal.md). Read straight from
 * get_parent_transmission(), which returns {haps: [0], normal: 0} for exactly
 * these parents, so the MoI calls and the parent of origin can't disagree.
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {Variant} variant
 * @returns {boolean}
 * @example
 * does_parent_pass_nothing(<son>, "dad", <non-PAR chrX>)       // → true
 * does_parent_pass_nothing(<daughter>, "dad", <non-PAR chrX>)  // → false (her father's X is hers)
 * does_parent_pass_nothing(<son>, "mom", <non-PAR chrY>)       // → true
 * does_parent_pass_nothing(<son>, "dad", <chrX PAR1>)          // → false (PAR: both parents pass a copy)
 * does_parent_pass_nothing(<kid>, "mom", <chr7>)               // → false
 */
function does_parent_pass_nothing(kid, which, variant) {
  var t = get_parent_transmission(kid, which, variant);
  return t !== null && t.normal === 0;
}

/**
 * The only safe way to ask whether a parent is there: reading s.mom / s.dad on a
 * sample without one throws in slivar (see the Sample typedef), so test with
 * `in` instead of reading the attribute.
 * @param {Sample} s
 * @param {string} which  "mom" | "dad"
 * @returns {boolean} true if that parent is in the family (s.mom / s.dad is set)
 * @example
 * has_parent(<duo kid, mother only>, "mom")  // → true
 * has_parent(<duo kid, mother only>, "dad")  // → false (and kid.dad would throw in slivar)
 * has_parent(<the mother, a founder>, "mom") // → false
 */
function has_parent(s, which) {
  return which in s;
}

// ===========================================================================
// Genome model: expected copy number (unchanged from the proposed file, P1)
// GRCh38 PAR, 1-based inclusive. GRCh37: PAR1 60001-2699520, PAR2 154931044-155260560.
// The pipeline refuses to run slivar unless vep_genome is GRCh38 (validateInputSamplesheet in
// subworkflows/local/utils_nfcore_cnv_post_processing_pipeline/main.nf).
// ===========================================================================
var CNV_PAR_X = [
  [10001, 2781479],
  [155701383, 156030895],
];
var CNV_PAR_Y = [
  [10001, 2781479],
  [56887903, 57217415],
];

/**
 * The record's chromosome without a "chr" prefix, so the rules work with both
 * naming styles.
 * @param {Variant} variant
 * @returns {string}
 * @example
 * fix_cnv_chrom_name(<CHROM "chr7">)  // → "7"
 * fix_cnv_chrom_name(<CHROM "chrX">)  // → "X"
 * fix_cnv_chrom_name(<CHROM "X">)     // → "X"
 * fix_cnv_chrom_name(<CHROM "chrM">)  // → "M"
 */
function fix_cnv_chrom_name(variant) {
  var c = variant.CHROM;
  return c.indexOf("chr") === 0 ? c.substring(3) : c;
}

/**
 * Does the CNV lie entirely inside a pseudoautosomal region (GRCh38, CNV_PAR_X /
 * CNV_PAR_Y)? There everyone has two copies, so the autosomal rules apply.
 * Always false off chrX / chrY.
 * @param {string} chrom   from fix_cnv_chrom_name()
 * @param {number} start1  1-based start
 * @param {number} end1    1-based end
 * @returns {boolean} true if the CNV lies entirely inside a pseudoautosomal region
 * @example
 * is_cnv_in_par("X", 1000000, 1010000)      // → true  (inside PAR1, 10,001-2,781,479)
 * is_cnv_in_par("Y", 56900000, 56950000)    // → true  (inside chrY PAR2)
 * is_cnv_in_par("X", 2780000, 2790000)      // → false (crosses the PAR1 end: see is_cnv_straddles_par)
 * is_cnv_in_par("X", 50000001, 50200000)    // → false (non-PAR)
 * is_cnv_in_par("7", 1000001, 1200000)      // → false (not chrX / chrY)
 */
function is_cnv_in_par(chrom, start1, end1) {
  // the CNV must lie entirely inside a PAR
  var pars = chrom === "X" ? CNV_PAR_X : chrom === "Y" ? CNV_PAR_Y : [];
  for (var i = 0; i < pars.length; i++) {
    if (start1 >= pars[i][0] && end1 <= pars[i][1]) return true;
  }
  return false;
}

/**
 * Does the CNV cross a PAR boundary (partly inside a PAR, partly outside)? Such a
 * record has no single normal CN for a male (2 in the PAR part, 1 outside), nor a
 * single copy a father passes on, so it isn't modelled at all.
 * @param {string} chrom   from fix_cnv_chrom_name()
 * @param {number} start1  1-based start
 * @param {number} end1    1-based end
 * @returns {boolean}
 * @example
 * is_cnv_straddles_par("X", 2780000, 2790000)      // → true  (crosses the PAR1 end, 2,781,479)
 * is_cnv_straddles_par("X", 155700000, 155710000)  // → true  (crosses the PAR2 start, 155,701,383)
 * is_cnv_straddles_par("X", 1000000, 1010000)      // → false (entirely inside PAR1)
 * is_cnv_straddles_par("X", 50000001, 50200000)    // → false (entirely outside)
 */
function is_cnv_straddles_par(chrom, start1, end1) {
  var pars = chrom === "X" ? CNV_PAR_X : chrom === "Y" ? CNV_PAR_Y : [];
  for (var i = 0; i < pars.length; i++) {
    var overlaps = start1 <= pars[i][1] && end1 >= pars[i][0];
    var inside = start1 >= pars[i][0] && end1 <= pars[i][1];
    if (overlaps && !inside) return true;
  }
  return false;
}

/**
 * Normal (reference) copy number for this sample on this record (P1). Every
 * rule works on the deviation from it (get_cn_dev), never on a hard-coded 2.
 * Whether a record is PAR is decided on its whole start..stop interval.
 * @param {Sample}  s
 * @param {Variant} variant
 * @returns {number} 2 on autosomes and PAR; 1 for a male on non-PAR X or Y;
 *   2 for a female on non-PAR X; -1 when it cannot be modelled (chrM, a
 *   record straddling a PAR boundary, unknown sex on non-PAR X/Y, female on
 *   non-PAR Y)
 * @example
 * get_expected_cn(boy, <chrX non-PAR record>)  // → 1
 * get_expected_cn(boy, <chr7 record>)          // → 2
 * @example
 * get_expected_cn(girl, <chrX non-PAR>)        // → 2
 * get_expected_cn(boy,  <chrX PAR1>)           // → 2  (his X and his Y both carry PAR1)
 * get_expected_cn(boy,  <chrY non-PAR>)        // → 1
 * get_expected_cn(girl, <chrY non-PAR>)        // → -1 (no Y)
 * get_expected_cn(<sex unknown>, <chrX non-PAR>)  // → -1
 * get_expected_cn(<sex unknown>, <chrX PAR1>)     // → 2  (sex doesn't matter in a PAR)
 * get_expected_cn(anyone, <chrM>)              // → -1
 * get_expected_cn(anyone, <chrX record crossing the PAR1 end>)  // → -1
 */
function get_expected_cn(s, variant) {
  var chrom = fix_cnv_chrom_name(variant);
  if (chrom === "M" || chrom === "MT") return -1;
  if (chrom !== "X" && chrom !== "Y") return 2;
  var start1 = variant.start + 1;
  var end1 = variant.stop; // not INFO.END: it throws when missing (Variant typedef)
  if (is_cnv_in_par(chrom, start1, end1)) return 2;
  if (is_cnv_straddles_par(chrom, start1, end1)) return -1; // no single normal CN: not modelled, for anyone
  if (s.sex !== "male" && s.sex !== "female") return -1;
  if (chrom === "X") return s.sex === "male" ? 1 : 2;
  return s.sex === "male" ? 1 : -1;
}

/**
 * Does this sample have a usable CN? It must be a non-negative integer; slivar's
 * missing sentinel (negative) and fractional CN (mosaic callers) count as no CN (P12).
 * @param {Sample} s
 * @returns {boolean}
 * @example
 * has_cn(<CN 2>)     // → true
 * has_cn(<CN 0>)     // → true  (a homozygous deletion is a real value)
 * has_cn(<CN -2147483648>)  // → false (what slivar passes when FORMAT/CN is "."). htslib's 32-bit "missing" marker, which slivar passes to the script when FORMAT/CN is .
 * has_cn(<CN 1.5>)   // → false (fractional, e.g. a mosaic call)
 */
function has_cn(s) {
  return typeof s.CN === "number" && s.CN >= 0 && s.CN === Math.round(s.CN);
}

/**
 * How far is this sample from its normal CN? The value every rule works with.
 * @param {Sample}  s
 * @param {Variant} variant
 * @returns {?number} CN - expected CN (-1 one copy lost, +1 one gained, ...),
 *   or null when there is no usable CN or the sample cannot be modelled
 * @example
 * get_cn_dev(<CN 1 on chr7>, v)             // → -1
 * get_cn_dev(<boy, CN 1 on non-PAR X>, v)   // → 0
 * @example
 * get_cn_dev(<CN 3 on chr7>, v)             // → +1
 * get_cn_dev(<boy, CN 0 on non-PAR X>, v)   // → -1 (his only X is deleted)
 * get_cn_dev(<girl, CN 1 on non-PAR X>, v)  // → -1
 * get_cn_dev(<CN missing on chr7>, v)       // → null
 * get_cn_dev(<CN 0 on chrM>, v)             // → null (not modelled)
 */
function get_cn_dev(s, variant) {
  var exp = get_expected_cn(s, variant);
  if (exp < 0 || !has_cn(s)) return null;
  return s.CN - exp;
}

// ===========================================================================
// Transmission: what each parent can pass on, and which parent the child's CNV
// came from. One engine, and one TRIO_TABLE row, for both the MoI call (trio_call)
// and the per-sample po_* tags (cn_po_origin), so the two can never disagree.
// ===========================================================================
// A diploid parent passes one of its two haplotypes. Its CN is split into the
// haplotype copy numbers it implies; CN >= 4 has several splits ([2,2], [3,1],
// [4,0]), so it can't be decided. On non-PAR chrX / chrY a father has a single
// copy, which he passes whole (X to a daughter, Y to a son), and some parents
// pass nothing at all (a father's X to a son, any mother's Y).
var CN_HAPLOTYPES = { 0: [0, 0], 1: [1, 0], 2: [1, 1], 3: [2, 1] };

/**
 * What this parent can pass to this child on this record. One half of the
 * transmission engine (with get_transmission_outcome), behind both the MoI
 * call and the per-sample po_* tags.
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {Variant} variant
 * @returns {?{haps: number[], normal: number}} the copy numbers the parent can
 *   actually pass, and the value a normal transmitted copy has (1, or 0 when that parent
 *   passes nothing here); null when it can't be known (parent absent or without
 *   CN, or a diploid CN >= 4)
 * @example
 * get_parent_transmission(<daughter>, "dad", <non-PAR chrX, dad CN 2>)   // → {haps: [2], normal: 1}
 *                                                                    // Dad has CN2 on a non-PAR X chromo.
 *                                                                    // He can therefore pass 2 copies, but
 *                                                                    // the normal is 1 copy (because chrx
 *                                                                    // non-par chr).
 * get_parent_transmission(<son>, "dad", <non-PAR chrX>)                  // → {haps: [0], normal: 0}
 * get_parent_transmission(<kid>, "mom", <chr7, mom CN 3>)                // → {haps: [2, 1], normal: 1}
 * @example
 * get_parent_transmission(<kid>, "mom", <chr7, mom CN 0>)        // → {haps: [0, 0], normal: 1} (always a deleted copy)
 * get_parent_transmission(<kid>, "mom", <chr7, mom CN 4>)        // → null ([2,2], [3,1] or [4,0]?)
 * get_parent_transmission(<duo kid>, "dad", <chr7, no father>)   // → null
 * get_parent_transmission(<son>, "dad", <non-PAR chrY, dad CN 2>)   // → {haps: [2], normal: 1}
 * get_parent_transmission(<son>, "mom", <non-PAR chrY, no mom CN>)  // → {haps: [0], normal: 0} (no Y from a mother)
 * get_parent_transmission(<son>, "dad", <chrX record crossing the PAR1 end>)  // → null (not modelled)
 */
function get_parent_transmission(kid, which, variant) {
  var chrom = fix_cnv_chrom_name(variant);
  var sexChrom = chrom === "X" || chrom === "Y";
  if (sexChrom && is_cnv_straddles_par(chrom, variant.start + 1, variant.stop)) return null; // not modelled (get_expected_cn)
  var hemi = sexChrom && !is_cnv_in_par(chrom, variant.start + 1, variant.stop);
  if (hemi) {
    // i.e. one gene on one chromo
    // Parents who pass nothing here, whatever their own CN (or whether they're in the PED).
    // Open question for the team: "any value" instead when the parent is absent / has no CN?
    // See docs/journal.md ("A parent who passes nothing on chrX/chrY").
    if (which === "dad" && (chrom === "X" ? kid.sex !== "female" : kid.sex !== "male")) return { haps: [0], normal: 0 };
    if (which === "mom" && chrom === "Y") return { haps: [0], normal: 0 };
  }
  if (!has_parent(kid, which)) return null;
  var p = kid[which];
  if (!has_cn(p)) return null;
  if (hemi && which === "dad") return { haps: [p.CN], normal: 1 }; // his single X (to a daughter) or Y (to a son)
  var h = CN_HAPLOTYPES[p.CN];
  return h ? { haps: h, normal: 1 } : null;
}

/**
 * How one passed copy compares with that parent's normal copy, seen from the
 * child: "normal", a "hit" (it changed the child's way: a deleted copy for a
 * loss, an extra copy for a gain), or "opposite" (it changed the other way).
 * @param {number} copy    the copy number passed
 * @param {number} normal  a normal copy's value for this parent (1, or 0 when it passes nothing)
 * @param {number} kd      the child's dev
 * @returns {string} "normal" | "hit" | "opposite"
 * @example
 * copy_effect(0, 1, -1)  // → "hit"       (a deleted copy, for a child who lost one)
 * copy_effect(2, 1, -1)  // → "opposite"  (an extra copy, for a child who lost one)
 * copy_effect(0, 0, -1)  // → "normal"    (a father's "nothing" to a son on non-PAR chrX)
 */
function copy_effect(copy, normal, kd) {
  if (copy === normal) return "normal";
  return (copy - normal) * kd > 0 ? "hit" : "opposite";
}

/**
 * Can this parent only pass normal copies here? A CN 2 parent ([1,1]), or one
 * who passes nothing (a father's X to a son, any mother's Y).
 * @param {{haps: number[], normal: number}} t  get_parent_transmission()'s result
 * @returns {boolean}
 * @example
 * passes_only_normal({haps: [1, 1], normal: 1})  // → true
 * passes_only_normal({haps: [1, 0], normal: 1})  // → false
 */
function passes_only_normal(t) {
  for (var i = 0; i < t.haps.length; i++) if (t.haps[i] !== t.normal) return false;
  return true;
}

/**
 * Every way one copy from the first parent plus one from the second adds up to
 * the child's CN, sorted into kinds by which copies are hits (copy_effect):
 * "first" (only the first parent's), "second", "both", or "opposite" (a copy
 * changed the other way). Two normal copies can't make a carrier.
 * @param {Sample}   kid
 * @param {number}   kd           the child's dev
 * @param {number[]} firstCopies  what the first parent can pass
 * @param {number}   firstNormal  its normal copy's value
 * @param {number[]} secondCopies what the second parent can pass
 * @param {number}   secondNormal its normal copy's value
 * @returns {Object} the kinds found, e.g. {first: true, second: true}; {} when no way fits
 * @example
 * // kid CN1: mom [1,0], dad [1,1] → only 0 + 1
 * get_way_kinds(kid, -1, [1, 0], 1, [1, 1], 1)  // → {first: true}
 * // kid CN1: mom [1,0], dad [1,0] → 0 + 1 or 1 + 0
 * get_way_kinds(kid, -1, [1, 0], 1, [1, 0], 1)  // → {first: true, second: true}
 */
function get_way_kinds(kid, kd, firstCopies, firstNormal, secondCopies, secondNormal) {
  var kinds = {};
  for (var i = 0; i < firstCopies.length; i++) {
    for (var j = 0; j < secondCopies.length; j++) {
      if (firstCopies[i] + secondCopies[j] !== kid.CN) continue; // this pair doesn't make the child's CN
      var a = copy_effect(firstCopies[i], firstNormal, kd),
        b = copy_effect(secondCopies[j], secondNormal, kd);
      if (a === "opposite" || b === "opposite") kinds.opposite = true;
      else if (a === "hit" && b === "hit") kinds.both = true;
      else if (a === "hit") kinds.first = true;
      else if (b === "hit") kinds.second = true;
    }
  }
  return kinds;
}

/**
 * The transmission outcome of one carrier child: TRIO_TABLE's outcome column
 * (see the legend there). From what each parent can pass
 * (get_parent_transmission), every way one copy from each adds up to the
 * child's CN (get_way_kinds). With one parent missing (not in the PED, or no
 * CN), get_duo_outcome() lets that parent pass any ordinary copy.
 * @param {Sample}  kid      a carrier child
 * @param {Variant} variant
 * @returns {string} "new" | "maternal" | "paternal" | "both" | "ambiguous" |
 *   "duo_new" | "duo_maternal" | "duo_paternal" | "duo_maternal_two" |
 *   "duo_paternal_two" | "duo_maybe_maternal" | "duo_maybe_paternal" | "undecided"
 * @example
 * // autosome unless noted; mom / dad copies in brackets
 * get_transmission_outcome(<kid 1, mom 1 [1,0], dad 2 [1,1]>)   // → "maternal"  (only 0 + 1; T4)
 * get_transmission_outcome(<kid 1, mom 1 [1,0], dad 1 [1,0]>)   // → "ambiguous" (0 + 1 or 1 + 0; T10)
 * get_transmission_outcome(<kid 0, mom 1 [1,0], dad 1 [1,0]>)   // → "both"      (only 0 + 0; T6)
 * get_transmission_outcome(<kid 1, mom 2 [1,1], dad 2 [1,1]>)   // → "new"       (no way, both normal; T1)
 * get_transmission_outcome(<kid 0, mom 2 [1,1], dad 1 [1,0]>)   // → "undecided" (no way, dad carries a CNV; T19)
 * get_transmission_outcome(<kid 3, mom 4, dad 1>)               // → "undecided" (mom's CN4 can't be split; T17)
 * get_transmission_outcome(<son 0, mom 1 [1,0]>, <non-PAR chrX>) // → "maternal" (dad passes no X to a son; X2)
 * get_transmission_outcome(<kid 1, mom 1 [1,0], no father>)     // → "duo_maybe_maternal" (D1)
 */
function get_transmission_outcome(kid, variant) {
  var kd = get_cn_dev(kid, variant);
  if (!kd) return "undecided"; // no CN, not modellable, or not a carrier: the callers check that first
  // A founder has no parent in the family: nothing to tell (without this, a male founder on
  // non-PAR chrX would see his "father" as known, since a father passes no X to a son).
  if (!has_parent(kid, "mom") && !has_parent(kid, "dad")) return "undecided";
  var m = get_parent_transmission(kid, "mom", variant); // m object is : (m.haps, m.normal)
  var d = get_parent_transmission(kid, "dad", variant); // d object is same as above
  if (m && d) {
    var k = get_way_kinds(kid, kd, m.haps, m.normal, d.haps, d.normal);
    var n = Object.keys(k).length;
    // if n == 0, means that we did not find any possible ways of transmission: new if neither parent
    // could have passed a hit, otherwise the Mendelian check decides.
    if (n === 0) return passes_only_normal(m) && passes_only_normal(d) ? "new" : "undecided";
    if (n > 1) return "ambiguous"; // if n > 1, the ways of transmission disagree.
    if (k.first) return "maternal";
    if (k.second) return "paternal";
    if (k.both) return "both";
    return "undecided"; // copies changed in opposite directions
  }
  if (m && is_parent_missing(kid, "dad", variant)) return get_duo_outcome(kid, kd, m, "maternal");
  if (d && is_parent_missing(kid, "mom", variant)) return get_duo_outcome(kid, kd, d, "paternal");
  return "undecided"; // a parent with CN >= 4, or no parent known
}

/**
 * The outcome for a child with one parent known and the other missing (not in
 * the PED, or no CN). The missing parent could have passed an ordinary copy:
 * 0, 1 or 2, like any parent with CN 0 to 3 (anything up to the child's CN with
 * MENDEL_ALLOW_RARE_SPLITS). This follows the SNV pipeline and Radiant (a het
 * child with a het mother and an unknown father is POSSIBLE_MOTHER, a hom-alt
 * one MOTHER), with a CNV carrier in place of a het.
 * @param {Sample} kid
 * @param {number} kd     the child's dev
 * @param {{haps: number[], normal: number}} known  what the known parent can pass
 * @param {string} side   "maternal" (the mother is known) | "paternal"
 * @returns {string} "duo_new" | "duo_<side>" | "duo_<side>_two" | "duo_maybe_<side>" | "undecided"
 * @example
 * // autosome, the mother known ([copies]), no father
 * get_duo_outcome(<kid 1>, -1, <mom 2 [1,1]>, "maternal")  // → "duo_new"            (D3)
 * get_duo_outcome(<kid 1>, -1, <mom 3 [2,1]>, "maternal")  // → "duo_new"            (only her 1 fits; D7)
 * get_duo_outcome(<kid 1>, -1, <mom 0 [0,0]>, "maternal")  // → "duo_maternal"       (her 0 + his 1; D5)
 * get_duo_outcome(<kid 0>, -2, <mom 1 [1,0]>, "maternal")  // → "duo_maternal_two"   (her 0 + his 0; D2)
 * get_duo_outcome(<kid 1>, -1, <mom 1 [1,0]>, "maternal")  // → "duo_maybe_maternal" (her 0 + his 1, or her 1 + his 0; D1)
 * get_duo_outcome(<kid 3>, +1, <mom 0 [0,0]>, "maternal")  // → "undecided"          (his copy would need 3)
 */
function get_duo_outcome(kid, kd, known, side) {
  if (passes_only_normal(known)) return "duo_new"; // its normal copy: the CNV is new, or from the missing parent
  var missing = MENDEL_ALLOW_RARE_SPLITS ? range(0, kid.CN) : [0, 1, 2];
  var k = get_way_kinds(kid, kd, known.haps, known.normal, missing, 1);
  if (k.opposite) return "undecided";
  if (k.first && k.second && !k.both) return "duo_maybe_" + side; // its hit, or its normal copy and the missing parent's hit
  if (Object.keys(k).length !== 1) return "undecided"; // no ordinary way: the Mendelian check decides
  if (k.first) return "duo_" + side; // its hit, and the missing parent's normal copy
  if (k.both) return "duo_" + side + "_two"; // its hit, and the missing parent's hit
  return "duo_new"; // only its normal copy fits: new, or from the missing parent
}

/**
 * Is this parent missing for this child: not in the PED, or in it without a
 * usable CN? A parent who passes nothing here isn't missing (it counts as
 * known and normal), and neither is one with CN 4 or more (known, but its
 * copies can't be told).
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {Variant} variant
 * @returns {boolean}
 */
function is_parent_missing(kid, which, variant) {
  var st = get_parent_status(kid, which, variant);
  return st === "absent" || st === "unknown";
}

/**
 * Parent of origin of one child's CNV: a transmission label for every child,
 * SNV-style. Entry point for the --sample-expr po_* lines. slivar calls it once
 * per sample and record (founders included: they have no parents, so they get
 * null), and values each po_* tag with the IDs of the samples whose result
 * matches. Affected status plays no part.
 *
 * It is the parent-of-origin column (4th) of the child's TRIO_TABLE row
 * (get_trio_row), the row whose child call cnv_moi() gives, so the MoI call and
 * the parent of origin always agree. "maternal" / "paternal" / "both" /
 * "ambiguous" come from which copies are hits; "denovo" needs both parents
 * known and normal, "possible_denovo" one of them missing; "possible_maternal"
 * / "possible_paternal" (the SNV pipeline's po_possible_mother /
 * po_possible_father) are for a duo whose known parent may have passed it. A
 * parent carrying a CNV that still can't explain the child's (a Mendelian
 * error) gets null: unlike the SNV pipeline's po_denovo, that isn't called de
 * novo.
 *
 * @param {Sample}  s        the child
 * @param {Variant} variant
 * @returns {?string} "maternal" | "paternal" | "both" | "ambiguous" | "denovo"
 *   | "possible_denovo" | "possible_maternal" | "possible_paternal", or null when
 *   it does not apply (not a carrier, no CN, unmodellable record, founder) or
 *   can't be told
 * @example
 * // kid CN1, mom CN1, dad CN2  → only mom 0 + dad 1 fits
 * cn_po_origin(kid, variant)  // → "maternal"
 * // kid CN1, mom CN1, dad CN1  → mom 0 + dad 1, or mom 1 + dad 0
 * cn_po_origin(kid, variant)  // → "ambiguous"
 * // non-PAR chrX: son CN0, mom CN1  → his only X is her deleted copy
 * cn_po_origin(son, variant)  // → "maternal"
 * @example
 * // chr7: kid CN0, mom CN1, dad CN1                  (T6)
 * cn_po_origin(kid, variant)  // → "both"
 * // chr7: UNAFFECTED kid CN1, mom CN1, dad CN2       (T13: affected status doesn't matter)
 * cn_po_origin(kid, variant)  // → "maternal"
 * // chr7: kid CN1, parents CN2                       (T1: no parent's copy is abnormal)
 * cn_po_origin(kid, variant)  // → "denovo"
 * // chr7 duo: kid CN1, mom CN2, no father            (D3: new, or from the father)
 * cn_po_origin(kid, variant)  // → "possible_denovo"
 * // chr7: kid CN0, mom CN2, dad CN1                  (T19: dad carries a CNV, yet can't explain it)
 * cn_po_origin(kid, variant)  // → null
 * // chr7 duo: kid CN1, mom CN1, no father            (D1: mom or the father?)
 * cn_po_origin(kid, variant)  // → "possible_maternal"
 * // chr7 duo: kid CN0, mom CN1, no father            (D2: one deleted copy is surely hers)
 * cn_po_origin(kid, variant)  // → "maternal"
 * // non-PAR chrX: daughter CN3, mom CN2, dad CN2     (X8: her father's duplicated X)
 * cn_po_origin(kid, variant)  // → "paternal"
 * // non-PAR chrY: son CN2, dad CN2                   (Y2: the mother isn't needed)
 * cn_po_origin(son, variant)  // → "paternal"
 * // non-PAR chrY: son CN2, dad CN1                   (Y3: the mother passes no Y)
 * cn_po_origin(son, variant)  // → "denovo"
 * // non-PAR chrX: son CN1 (normal for him), or a record crossing the PAR1 end  (X1, X4); a founder
 * cn_po_origin(s, variant)    // → null
 */
function cn_po_origin(s, variant) {
  var exp = get_expected_cn(s, variant);
  if (exp < 0 || !has_cn(s) || s.CN === exp) return null;
  return get_trio_row(s, variant)[3];
}

/**
 * What we know about one parent of this child on this record.
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {Variant} variant
 * @returns {string} "nothing" (passes no copy here: does_parent_pass_nothing),
 *   "absent" (not in the PED), "unknown" (no usable CN, or not modellable),
 *   "normal" (dev 0) or "carrier" (dev ≠ 0)
 * @example
 * get_parent_status(<son>, "dad", <non-PAR chrX>)         // → "nothing"
 * get_parent_status(<duo kid, mother only>, "dad", <chr7>) // → "absent"
 * get_parent_status(<kid>, "mom", <chr7, mom CN missing>)  // → "unknown"
 * get_parent_status(<kid>, "mom", <chr7, mom CN 2>)        // → "normal"
 * get_parent_status(<kid>, "dad", <chr7, dad CN 1>)        // → "carrier"
 */
function get_parent_status(kid, which, variant) {
  if (does_parent_pass_nothing(kid, which, variant)) return "nothing";
  if (!has_parent(kid, which)) return "absent";
  var d = get_cn_dev(kid[which], variant);
  if (d === null) return "unknown";
  return d === 0 ? "normal" : "carrier";
}

// ===========================================================================
// Is the child's CN transmittable at all?  (mendelian_inconsistent vs unresolved)
// ===========================================================================
// Same parent rules as get_parent_transmission(), but permissive where it can't
// decide: CN >= 4 allows any split ([2,2], [3,1], [4,0] -> passes 0..4), and an
// absent parent or one without CN could have passed anything.
// Set to true to also allow the rare hidden splits (CN2 = [2,0], CN3 = [3,0], ...).
// Then a diploid parent can pass anything from 0 to its CN, and only "child has
// more copies than both parents together" is still impossible — no deletion ever is.
var MENDEL_ALLOW_RARE_SPLITS = false;
/**
 * Could this child's CN have been transmitted by its parents?
 * The real question is "can we prove it impossible?": with no CN, or a
 * missing parent, nothing can be proven, so the answer is true.
 * trio_call() asks it when only the catch-all "no_pattern" row matches.
 * @param {Sample}  kid      the child (a founder would always come out true)
 * @param {Variant} variant
 * @returns {boolean} true unless no copy number mom can pass plus one dad can
 *   pass adds up to kid.CN (see transmissions)
 * @example
 * // kid CN0, mom CN0 ([0,0]), dad CN2 ([1,1]): 0 + 1 ≠ 0  (T14)
 * is_mendelian_possible(kid, variant)  // → false
 * @example
 * // chr7: kid CN1, mom CN1, dad CN1: 0 + 1 = 1           (T10)
 * is_mendelian_possible(kid, variant)  // → true
 * // chr7: kid CN0, mom CN2 ([1,1]), dad CN1: mom has no deleted copy to pass  (T19)
 * is_mendelian_possible(kid, variant)  // → false
 * // chr7: kid CN3, mom CN4 (any split: 0..4), dad CN3: e.g. 2 + 1  (T16)
 * is_mendelian_possible(kid, variant)  // → true
 * // chr7: kid CN1, parents CN2: a de novo isn't transmittable  (T1; why it only runs for the catch-all row)
 * is_mendelian_possible(kid, variant)  // → false
 * // a founder, or a kid with no CN
 * is_mendelian_possible(s, variant)    // → true (nothing to prove)
 * // chr7 duo: kid CN0, mom CN2, no father (he could have passed anything, but she can't pass 0)
 * is_mendelian_possible(kid, variant)  // → false
 * // non-PAR chrX: son CN2, mom CN0, no father (a father passes no X to a son anyway)
 * is_mendelian_possible(son, variant)  // → false (see the open question in docs/journal.md)
 */
function is_mendelian_possible(kid, variant) {
  if (!has_cn(kid)) return true; // With no CN there's nothing to test, so it can't prove anything is wrong, and it answers "possible" (true).
  var fromMom = transmissions(kid, "mom", variant);
  var fromDad = transmissions(kid, "dad", variant);
  for (var i = 0; i < fromMom.length; i++)
    for (var j = 0; j < fromDad.length; j++) if (fromMom[i] + fromDad[j] === kid.CN) return true;
  return false;
}

/**
 * Copy numbers this parent can pass to this child on this record, for
 * is_mendelian_possible(): get_parent_transmission(), with permissive fallbacks.
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {Variant} variant
 * @returns {number[]}
 *   parent passes nothing here       → [0] (a father's X to a son, any mother's Y)
 *   non-PAR X/Y, the father's copy   → [dad.CN] (X to a daughter, Y to a son)
 *   CN 0-3                           → the two halves in CN_HAPLOTYPES (CN2 → [1,1])
 *   CN >= 4, or MENDEL_ALLOW_RARE_SPLITS → 0..parent CN
 *   absent parent or no CN           → 0..kid.CN (could have passed anything)
 * @example
 * transmissions(kid, "dad", <chr7, dad CN 2>)               // → [1, 1]
 * transmissions(kid, "mom", <chr7, mom CN 4>)               // → [0, 1, 2, 3, 4]
 * transmissions(<kid CN 3>, "dad", <chr7, no father>)       // → [0, 1, 2, 3]
 * transmissions(<son>, "dad", <non-PAR chrX, no father>)    // → [0]
 * transmissions(<daughter>, "dad", <non-PAR chrX, dad CN 2>)  // → [2]
 */
function transmissions(kid, which, variant) {
  var t = get_parent_transmission(kid, which, variant);
  var p = has_parent(kid, which) ? kid[which] : null; // kid[which] alone throws for a founder in slivar
  if (t === null) return p !== null && has_cn(p) ? range(0, p.CN) : range(0, kid.CN);
  if (MENDEL_ALLOW_RARE_SPLITS && t.haps.length === 2) return range(0, p.CN);
  return t.haps;
}

/**
 * The integers from a to b, both included (empty when b < a).
 * @param {number} a
 * @param {number} b
 * @returns {number[]}
 * @example
 * range(0, 3)  // → [0, 1, 2, 3]
 * range(2, 2)  // → [2]
 * range(3, 1)  // → []
 */
function range(a, b) {
  var out = [];
  for (var i = a; i <= b; i++) out.push(i);
  return out;
}

if (typeof module !== "undefined")
  module.exports = {
    TRIO_TABLE: TRIO_TABLE,
    X_LINKED_CALLS: X_LINKED_CALLS,
    cnv_moi: cnv_moi,
    has_moi_call: has_moi_call,
    trio_call: trio_call,
    cn_po_origin: cn_po_origin,
    get_trio_row: get_trio_row,
    get_transmission_outcome: get_transmission_outcome,
    get_expected_cn: get_expected_cn,
    get_cn_dev: get_cn_dev,
    has_cn: has_cn,
  };
