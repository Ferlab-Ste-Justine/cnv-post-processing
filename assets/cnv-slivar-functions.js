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
//   --family-expr 'de_novo:has_moi_call(fam, variant, "de_novo")' \
//   --family-expr 'de_novo_suspect:has_moi_call(fam, variant, "de_novo_suspect")' \
//   --family-expr 'inherited_maternal:has_moi_call(fam, variant, "inherited_maternal")' \
//   --family-expr 'inherited_paternal:has_moi_call(fam, variant, "inherited_paternal")' \
//   --family-expr 'recessive:has_moi_call(fam, variant, "recessive")' \
//   --family-expr 'recessive_lowconf:has_moi_call(fam, variant, "recessive_lowconf")' \
//   --family-expr 'possible_de_novo:has_moi_call(fam, variant, "possible_de_novo")' \
//   --family-expr 'possible_inherited_maternal:has_moi_call(fam, variant, "possible_inherited_maternal")' \
//   --family-expr 'possible_inherited_paternal:has_moi_call(fam, variant, "possible_inherited_paternal")' \
//   --family-expr 'inheritance_unresolved:has_moi_call(fam, variant, "inheritance_unresolved")' \
//   --family-expr 'mendelian_inconsistent:has_moi_call(fam, variant, "mendelian_inconsistent")' \
//   --family-expr 'unknown_cn:has_moi_call(fam, variant, "unknown_cn")' \
//   --family-expr 'x_linked_recessive_de_novo:has_moi_call(fam, variant, "x_linked_recessive_de_novo")' \
//   --family-expr 'x_linked_dominant_de_novo:has_moi_call(fam, variant, "x_linked_dominant_de_novo")' \
//   --family-expr 'x_linked_recessive:has_moi_call(fam, variant, "x_linked_recessive")' \
//   --sample-expr 'po_mother:cn_po_origin(sample, variant) == "maternal"' \
//   --sample-expr 'po_father:cn_po_origin(sample, variant) == "paternal"' \
//   --sample-expr 'po_both:cn_po_origin(sample, variant) == "both"' \
//   --sample-expr 'po_ambiguous:cn_po_origin(sample, variant) == "ambiguous"' \
//   --sample-expr 'po_denovo:cn_po_origin(sample, variant) == "denovo"' \
//   --sample-expr 'po_possible_denovo:cn_po_origin(sample, variant) == "possible_denovo"'
//
// Vocabulary
//   dev      = CN - expected CN. Expected CN is 2 on autosomes and PAR, 1 for a
//              male on non-PAR chrX / chrY. So dev -1 = one copy lost, +1 = one
//              copy gained, -2 = both copies lost, +2 = two copies gained.
//   carrier  = a sample with dev != 0.
//   carrier child = an AFFECTED carrier with at least one parent in the family:
//              the only samples cnv_moi() gives a call.
//   passes nothing = a parent who gives this child no copy of this chromosome:
//              a father's X to a son, any mother's Y. Such a parent is read as
//              "normal", and their CN isn't needed (does_parent_pass_nothing()).
//   mendelian_inconsistent = no way for the parents to transmit the child's CN
//              (see is_mendelian_possible() below) -> sample swap, miscall, missed call.
//   inheritance_unresolved = transmission is possible, but no table row describes it
//              (e.g. kid 3, mom 4, dad 3) -> needs a human.
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
// │       (only one of them without CN: that parent counts as "absent", like a duo)
// │
// ├─④ [get_cn_dev(s)] == 0 ?  not a carrier ─────────────────────► null
// │
// └─⑤ trio_call(s, variant), then get_x_linked_call(): on non-PAR chrX, the X_LINKED_CALLS name
//       │
//       │  kd     = [get_cn_dev(kid)]                                 // How far is the kid from normal. -1 = 1 copy lost, 1 = 1 copy gained.
//       │  kidTok = populate_kid_tokens(kd)                 → ["any","±1",...] // returns array of strings that is then used
//       │  momTok = does_parent_pass_nothing(kid,"mom") ? "normal"
//       │             : populate_parent_token(kid,"mom",kd)
//       │                 ├ [has_parent] no, or no [has_cn] → "absent"
//       │                 └ [get_cn_dev(mom)] → "normal" | "same" | "half" | "other"
//       │  dadTok = the same for "dad"
//       │  origin = [get_transmission_origin(kid)]          → "maternal" | "paternal" | "both" | "ambiguous" | null
//       │
//       │  scan TRIO_TABLE, first row matching (kidTok, momTok, dadTok, origin)
//       │    call ≠ "no_pattern" ─────────────────────► that call
//       │    call = "no_pattern" ─► is_mendelian_possible(kid)
//       │                             yes → "inheritance_unresolved" (no pattern fits, but these parents
//       │                                                             could have produced this child's CN)
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
// "ambiguous", po_denovo "denovo", po_possible_denovo "possible_denovo". Independent
// of cnv_moi(), but the same [get_transmission_origin] fills TRIO_TABLE's origin
// column, so the two always agree.
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
// ├─② get_transmission_origin(s, variant)
// │     │
// │     ├ m = get_parent_transmission(s, "mom")    what mom can pass, and a normal copy's value
// │     ├ d = get_parent_transmission(s, "dad")    same for dad
// │     │     ├ PAR-straddling record                      → null (not modelled). When a record crosses a PAR boundary, return null.
// │     │     ├ non-PAR X: dad → son (or unknown sex)      → {haps: [0], normal: 0}  passes nothing
// │     │     ├ non-PAR Y: dad → daughter, or mom → anyone → {haps: [0], normal: 0}  passes nothing
// │     │     │    (these two hold even if that parent is absent / has no CN:
// │     │     │     open question, ** see docs/journal.md) **
// │     │     ├ parent missing, or no [has_cn]             → null
// │     │     ├ non-PAR X/Y: the father's single copy      → {haps: [dad.CN], normal: 1}
// │     │     │    (his X to a daughter, his Y to a son)
// │     │     └ diploid: CN_HAPLOTYPES[CN], CN 0-3         → {haps: [a, b], normal: 1}
// │     │                CN >= 4 (several splits fit)      → null
// │     │     m or d null ────────────────────────────────────────────► ③
// │     │
// │     ├ ways = the distinct (mom copy, dad copy) pairs with mom copy + dad copy == s.CN
// │     │     0 ways  (de novo, or a Mendelian error) ────────────────► ③
// │     │     2+ ways ────────────────────────────────────────────────► "ambiguous"
// │     │
// │     └ 1 way: which transmitted copy differs from that parent's normal one?
// │           mom's only ─────────────────────────────────────────────► "maternal"
// │           dad's only ─────────────────────────────────────────────► "paternal"
// │           both ───────────────────────────────────────────────────► "both"
// │     (any origin found is the answer)
// │
// └─③ get_denovo_origin(s, variant)     only when ② found nothing
//       mom, dad = get_parent_status(s, "mom" | "dad")
//         "nothing" (passes no copy here) | "absent" | "unknown" (no CN) | "normal" | "carrier"
//       both known and normal ("normal" or "nothing") ─────────────► "denovo"
//       one known and normal, the other "absent" or "unknown" ──────► "possible_denovo"
//       otherwise (a carrier parent: Mendelian error; or too unknown) ► null
//
//   Examples (autosome unless noted; mom / dad copies in brackets):
//     kid 1, mom 1 [1,0], dad 2 [1,1]       only 0 + 1          → "maternal"   (harness T4)
//     kid 1, mom 1 [1,0], dad 1 [1,0]       0 + 1 or 1 + 0      → "ambiguous"  (T10)
//     kid 0, mom 1 [1,0], dad 1 [1,0]       only 0 + 0          → "both"       (T6)
//     kid 1, mom 2 [1,1], dad 2 [1,1]       no way, both normal → "denovo"     (T1)
//     duo: kid 1, mom 2 [1,1], no father    no way, dad unknown → "possible_denovo"  (D3)
//     kid 0, mom 2 [1,1], dad 1 [1,0]       no way, dad carrier → null         (T19)
//     non-PAR X: son 0, mom 1 [1,0]         dad passes nothing  → "maternal"   (X2)
//     non-PAR X: daughter 3, mom 2 [1,1], dad 2 [2]   1 + 2     → "paternal"   (X8)
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
// │  get_transmission_origin(kid, variant)             "which parent's copy is abnormal?"      │
// │    ├ get_parent_transmission(kid, "mom"|"dad")     copies it can pass:                     │
// │    │    diploid → CN_HAPLOTYPES[CN] (CN >= 4 → can't tell);                                │
// │    │    non-PAR X/Y father → [his CN] to a daughter (X) / son (Y),                         │
// │    │    [0] to the other sex; mother on Y → [0]                                            │
// │    └ the distinct (mom copy, dad copy) pairs that add up to kid.CN:                        │
// │         0 → null, 2+ → "ambiguous", 1 → maternal / paternal / both                         │
// │  Also behind cn_po_origin() (--sample-expr), so MoI and PO agree.                          │
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
// TABLE 1 — TRIO_TABLE: one affected carrier child and its parents
// ===========================================================================
// kid column (the child's dev):
//   "±1"    dev is -1 or +1            (one copy lost or gained)
//   "±2+"   |dev| >= 2                 (CN 0 or CN >= 4 on an autosome)
//   "even"  dev is -2, +2, -4, ...     (could be one copy from each parent)
//   "any"
// mom / dad columns (the parent's dev, read RELATIVE to the kid):
//   "normal"  parent dev = 0
//   "same"    parent dev = kid dev            (parent has the same CNV)
//   "half"    parent dev = kid dev / 2        (parent is a carrier of one copy)
//   "absent"  parent not in the PED / VCF, or without a usable CN (as SNV/Radiant
//             treat a parent's no-call); is_cn_missing() gives unknown_cn instead
//             when no needed parent has a CN
//   "any"
// A parent who passes nothing here (a father's X to a son, any mother's Y) is read as
// "normal", even when absent from the PED or without a CN.
// origin column (get_transmission_origin(): which parent's transmitted copy is abnormal,
// from every way the parents' copies can add up to the child's CN; the same
// computation cn_po_origin() reports per sample, so the two always agree):
//   "maternal" / "paternal"   exactly one way fits, and only that parent's copy is abnormal
//   "any"                     no requirement (de novo and duo rows rely on the columns before)
// First matching row wins.
// prettier-ignore
var TRIO_TABLE = [
  // kid      mom        dad        origin       child call                      harness rows
  ["±1",     "normal",  "normal",  "any",       "de_novo"],                     // T1 T2 X3
  ["±2+",    "normal",  "normal",  "any",       "de_novo_suspect"],             // T3 T26 (P5: two hits at one locus is implausible)
  ["any",    "any",     "any",     "maternal",  "inherited_maternal"],          // T4 T7 T23 T24 X2, Q9 kid
  ["any",    "any",     "any",     "paternal",  "inherited_paternal"],          // T5, Q9 sib
  ["even",   "half",    "half",    "any",       "recessive"],                   // T6 T9 (T9 = CN4 from two CN3, see P7)
  ["any",    "normal",  "absent",  "any",       "possible_de_novo"],            // D3
  ["any",    "absent",  "normal",  "any",       "possible_de_novo"],
  ["±1",     "same",    "absent",  "any",       "possible_inherited_maternal"], // D1
  ["±1",     "absent",  "same",    "any",       "possible_inherited_paternal"],
  ["even",   "half",    "absent",  "any",       "recessive_lowconf"],           // D2
  ["even",   "absent",  "half",    "any",       "recessive_lowconf"],
  ["any",    "any",     "any",     "any",       "no_pattern"]                   // -> inheritance_unresolved (T8 T10) or mendelian_inconsistent (T14 T25)
];

// On non-PAR chrX, a child's call takes the SNV pipeline's X-linked name (Radiant's
// transmission modes), by the child's sex; other calls, PAR and the other chromosomes
// keep the names above (get_x_linked_call). A daughter with one inherited copy keeps
// inherited_*: dominant or carrier would depend on who is affected, which the calls
// don't look at. de_novo_suspect and the duo calls keep their names too.
// prettier-ignore
var X_LINKED_CALLS = {
  male:   { de_novo: "x_linked_recessive_de_novo",                       // X3: his only X changed
            inherited_maternal: "x_linked_recessive" },                  // X2: his only X is his mother's
  female: { de_novo: "x_linked_dominant_de_novo",                        // X7
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
 * @param {string}   call     the tag being tested, e.g. "de_novo"
 * @returns {boolean}
 * @example
 * // T4: affected kid CN1, unaffected mom CN1, dad CN2
 * has_moi_call(fam, variant, "inherited_maternal")  // → true  (the kid's call)
 * has_moi_call(fam, variant, "de_novo")             // → false
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
 *       for each of the 15 MoI --family-expr:  ← slivar
 *          has_moi_call(fam, variant, tag)
 *             for each member s:
 *                cnv_moi(s, variant)           ← steps ①–⑤
 *
 * @param {Sample}  s        any sample; only affected children get a call
 * @param {Variant} variant  the CNV record
 * @returns {?string} one call, e.g. "de_novo", "inherited_maternal", "unknown_cn",
 *   "inheritance_unresolved"; null when the sample gets no tag (unaffected,
 *   founder, not a carrier, unmodellable record)
 * @example
 * // affected kid CN1, unaffected mom CN1, dad CN2, chr7  (harness T4; see assets/moi_scenarios.md)
 * cnv_moi(kid, variant)  // → "inherited_maternal"
 * cnv_moi(mom, variant)  // → null (a founder)
 * @example
 * // Each example is one child on one record; CNs of affected kid / mom / dad, chr7 unless noted.
 * // kid 1, mom 2, dad 2                         both parents normal       (T1)
 * cnv_moi(kid, variant)  // → "de_novo"
 * // kid 1, mom 2, AFFECTED dad 1                the call ignores who is affected  (T5)
 * cnv_moi(kid, variant)  // → "inherited_paternal"
 * // kid 0, mom 1, dad 1                         one deleted copy from each parent  (T6)
 * cnv_moi(kid, variant)  // → "recessive"
 * // kid 3, mom 3, dad 3                         inherited, but from which parent?  (T8)
 * cnv_moi(kid, variant)  // → "inheritance_unresolved"
 * // kid 0, mom 0, dad 2                         dad can't pass a deleted copy  (T14)
 * cnv_moi(kid, variant)  // → "mendelian_inconsistent"
 * // kid 1, mom 2, dad CN missing                read as a duo without a father  (T11)
 * cnv_moi(kid, variant)  // → "possible_de_novo"
 * // duo: kid 1, mom 2, no father in the PED     (D3)
 * cnv_moi(kid, variant)  // → "possible_de_novo"
 * // non-PAR chrX: son 0, carrier mom 1, dad 1   his only X is hers  (X2)
 * cnv_moi(son, variant)  // → "x_linked_recessive"
 * // non-PAR chrX: daughter 1, mom 2, dad 1      de novo in a daughter  (X7)
 * cnv_moi(daughter, variant)  // → "x_linked_dominant_de_novo"
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
 * TRIO_TABLE call for one affected carrier child (step ⑤ of the map).
 *
 * Turns the child's and each parent's CN into tokens (populate_kid_tokens,
 * populate_parent_token) and works out the origin (get_transmission_origin), then
 * returns the call (5th column) of the first matching TRIO_TABLE row. The
 * catch-all "no_pattern" row is resolved here with is_mendelian_possible().
 * Exported for the harness / debugging.
 *
 * The call describes this child and its parents only (who is affected among
 * them is not looked at). A parent who passes nothing here is read as "normal"
 * (does_parent_pass_nothing).
 *
 * @param {Sample}  kid      affected carrier with at least one parent in the family
 * @param {Variant} variant
 * @returns {string} the child call, e.g. "de_novo", "inherited_maternal",
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
 * // tokens: kid ["any","±1"], mom "same", dad "normal"; origin "maternal"
 * trio_call(kid, variant)  // → "inherited_maternal"
 * @example
 * // kid / mom / dad CNs, chr7 unless noted; the TRIO_TABLE row that matches in brackets
 * // kid 1, mom 2, dad 2             [±1, normal, normal]                 (T1)
 * trio_call(kid, variant)  // → "de_novo"
 * // kid 0, mom 2, dad 2             [±2+, normal, normal]                (T3)
 * trio_call(kid, variant)  // → "de_novo_suspect"
 * // kid 0, mom 1, dad 1             [even, half, half]                   (T6)
 * trio_call(kid, variant)  // → "recessive"
 * // kid 1, mom 1, dad 3             origin "maternal": dad's CN3 passes 1 or 2, only 1 fits  (Q9)
 * trio_call(kid, variant)  // → "inherited_maternal"
 * // duo: kid 1, mom 1, no father    [±1, same, absent]                   (D1)
 * trio_call(kid, variant)  // → "possible_inherited_maternal"
 * // kid 3, mom 3, dad 3             no_pattern, but transmittable        (T8)
 * trio_call(kid, variant)  // → "inheritance_unresolved"
 * // kid 0, mom 0, dad 2             no_pattern, not transmittable        (T14)
 * trio_call(kid, variant)  // → "mendelian_inconsistent"
 * // non-PAR chrX: son 0, mom 1, dad 1   dad read as "normal" (he passes no X to a son)  (X2)
 * trio_call(kid, variant)  // → "inherited_maternal"
 * // non-PAR chrY: son 2, mom (no Y CN), dad 1   mom read as "normal" (no Y from a mother)  (Y3)
 * trio_call(kid, variant)  // → "de_novo"
 */
function trio_call(kid, variant) {
  var kd = get_cn_dev(kid, variant);
  var kidTok = populate_kid_tokens(kd); // kid gets an array of strings. Ex could be ["any", "±2+", "even"] or ["any", "±1"] . Always contains "any" at first pos. "any" can be seen as a placeholder that triggers when we found no matches in the TRIO_TABLE.
  // each parent (momTok and dadTok gets one single string value. ex: "normal", "same", "half", "absent" or "other"
  var momTok = does_parent_pass_nothing(kid, "mom", variant)
    ? "normal"
    : populate_parent_token(kid, "mom", kd, variant);
  var dadTok = does_parent_pass_nothing(kid, "dad", variant)
    ? "normal"
    : populate_parent_token(kid, "dad", kd, variant);
  var origin = get_transmission_origin(kid, variant); // "maternal" | "paternal" | "both" | "ambiguous" | null -- same answer as cn_po_origin()
  for (var i = 0; i < TRIO_TABLE.length; i++) {
    var row = TRIO_TABLE[i];
    // excerpt from the TRIO table.
    // kid       mom        dad        origin       child call                      harness rows
    // ["±1",    "normal",  "normal",  "any",       "de_novo"],                     // T1 T2 X3
    // ["±2+",   "normal",  "normal",  "any",       "de_novo_suspect"],             // T3   (P5: two hits at one locus is implausible)
    //["any",    "any",     "any",     "maternal",  "inherited_maternal"],          // T4 T7 X2, Q9 kid
    if (
      // ex: is "+/-1 found in ['any', '+1']"
      kidTok.indexOf(row[0]) >= 0 && // is the row's kid value in the kidTok's list? again kid's list could be  ["any", "±1"]. yes ->
      (row[1] === "any" || row[1] === momTok) && // row accepts any mom, or exactly this mom? yes ->
      (row[2] === "any" || row[2] === dadTok) && // row accepts any dad, or exactly this dad? yes ->
      (row[3] === "any" || row[3] === origin)
    ) {
      // row accepts any origin, or exactly this one? If all 4 yes, return the value of row[4] i.e. the 5th value of the lookup table.
      //      Except if row[4] of the lookup table is "no_pattern". If that is the case, we'll go with the is_mendelian_possible() which assess if:
      //          "mendelian_inconsistent" (none of our patterns fit this child. Could these parents have produced
      //                                    this child's CN at all -> No)
      //          "inheritance_unresolved" (none of our patterns fit this child. Could these parents have produced
      //                                    this child's CN at all -> Yes)
      if (row[4] !== "no_pattern") return row[4]; // So basically returns the value of the 5th column of TRIO_TABLE's current row unless that value is "no_pattern"
      return is_mendelian_possible(kid, variant) ? "inheritance_unresolved" : "mendelian_inconsistent";
    }
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
 * get_x_linked_call(<son>, <non-PAR chrX>, "de_novo")                 // → "x_linked_recessive_de_novo"
 * get_x_linked_call(<daughter>, <non-PAR chrX>, "de_novo")            // → "x_linked_dominant_de_novo"
 * get_x_linked_call(<son>, <non-PAR chrX>, "inherited_maternal")      // → "x_linked_recessive"
 * get_x_linked_call(<daughter>, <non-PAR chrX>, "inherited_maternal") // → "inherited_maternal" (unchanged)
 * get_x_linked_call(<son>, <chrX PAR1>, "de_novo")                    // → "de_novo" (PAR: autosomal names)
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
/**
 * How one parent compares with the child: the TRIO_TABLE mom / dad column value.
 * Checked in this order: absent (not in the family, or no usable CN), normal,
 * same, half, other. "same" means the parent carries the child's CNV as is (a
 * one-copy carrier passing it on); "half" means the parent carries half the
 * child's deviation (a recessive carrier); anything else is "other".
 * @param {Sample}  kid
 * @param {string}  which    "mom" | "dad"
 * @param {number}  kd       the child's dev
 * @param {Variant} variant
 * @returns {string} "absent" | "normal" (dev 0) | "same" (dev = kd)
 *   | "half" (dev × 2 = kd) | "other"
 * @example
 * // kid CN0 (kd = -2), mom CN1 (dev -1)
 * populate_parent_token(kid, "mom", -2, variant)  // → "half"
 * @example
 * // kid CN1 (kd = -1), chr7:
 * populate_parent_token(kid, "mom", -1, variant)  // mom CN2 (dev 0)  → "normal"
 * populate_parent_token(kid, "mom", -1, variant)  // mom CN1 (dev -1) → "same"
 * populate_parent_token(kid, "mom", -1, variant)  // mom CN0 (dev -2) → "other"
 * populate_parent_token(kid, "dad", -1, variant)  // duo, no father   → "absent"
 * populate_parent_token(kid, "dad", -1, variant)  // dad CN missing   → "absent"  (T11)
 * // non-PAR chrX: daughter CN3 (kd = +1), dad CN2 (dev +1 against his expected 1)
 * populate_parent_token(kid, "dad", 1, variant)   // → "same"
 */
function populate_parent_token(kid, which, kd, variant) {
  if (!has_parent(kid, which) || !has_cn(kid[which])) return "absent"; // no CN: read as a duo (is_cn_missing)
  var pd = get_cn_dev(kid[which], variant);
  if (pd === 0) return "normal";
  if (pd === kd) return "same";
  if (pd * 2 === kd) return "half";
  return "other";
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
 * otherwise that parent counts as "absent" (populate_parent_token), as for a
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
// came from. One engine for both the TRIO_TABLE origin column (trio_call) and
// the per-sample po_* tags (cn_po_origin), so the two can never disagree.
// ===========================================================================
// A diploid parent passes one of its two haplotypes. Its CN is split into the
// haplotype copy numbers it implies; CN >= 4 has several splits ([2,2], [3,1],
// [4,0]), so it can't be decided. On non-PAR chrX / chrY a father has a single
// copy, which he passes whole (X to a daughter, Y to a son), and some parents
// pass nothing at all (a father's X to a son, any mother's Y).
var CN_HAPLOTYPES = { 0: [0, 0], 1: [1, 0], 2: [1, 1], 3: [2, 1] };

/**
 * What this parent can pass to this child on this record. One half of the
 * transmission engine (with get_transmission_origin), behind both TRIO_TABLE's
 * origin column and the per-sample po_* tags.
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
 * Which parent's transmitted copy is abnormal, from every distinct way the
 * parents' copies (get_parent_transmission) can add up to the child's CN:
 *   0 ways -> null       (not explained by transmission: de novo, or a Mendelian error)
 *   2+     -> "ambiguous"
 *   1      -> "maternal" / "paternal" / "both", by which transmitted copy differs
 *             from that parent's normal one
 * @param {Sample}  kid
 * @param {Variant} variant
 * @returns {?string} "maternal" | "paternal" | "both" | "ambiguous" | null
 * @example
 * // chr7: kid CN1, mom CN1 [1,0], dad CN2 [1,1]  → only mom 0 + dad 1 fits
 * get_transmission_origin(kid, variant)  // → "maternal"
 * // non-PAR chrX: daughter CN3, dad CN2 [2], mom CN2 [1,1]  → only mom 1 + dad 2 fits
 * get_transmission_origin(kid, variant)  // → "paternal"
 * @example
 * // chr7: kid CN1, mom CN1 [1,0], dad CN1 [1,0]  → 0 + 1 or 1 + 0  (T10)
 * get_transmission_origin(kid, variant)  // → "ambiguous"
 * // chr7: kid CN0, mom CN1 [1,0], dad CN1 [1,0]  → only 0 + 0  (T6)
 * get_transmission_origin(kid, variant)  // → "both"
 * // chr7: kid CN1, mom CN1 [1,0], dad CN3 [2,1]  → only 0 + 1 (dad's CN3 is a separate dup)  (Q9)
 * get_transmission_origin(kid, variant)  // → "maternal"
 * // chr7: kid CN1, parents CN2 [1,1]             → no way (de novo)  (T1)
 * get_transmission_origin(kid, variant)  // → null
 * // chr7: kid CN3, mom CN4, dad CN1              → mom's CN4 can't be split  (T17)
 * get_transmission_origin(kid, variant)  // → null
 * // non-PAR chrX: son CN0, mom CN1 [1,0]          → his only X is her deleted copy  (X2)
 * get_transmission_origin(kid, variant)  // → "maternal"
 */
function get_transmission_origin(kid, variant) {
  if (!has_cn(kid)) return null;
  var m = get_parent_transmission(kid, "mom", variant); // m object is : (m.haps, m.normal)
  var d = get_parent_transmission(kid, "dad", variant); // d object is same as above
  if (!m || !d) return null;
  var ways = {},
    n = 0,
    fromMom,
    fromDad;
  // ex: kid.CN = 1, m.haps = [1, 0] and d.haps = [1, 0]
  for (var i = 0; i < m.haps.length; i++) {
    for (var j = 0; j < d.haps.length; j++) {
      if (m.haps[i] + d.haps[j] !== kid.CN) continue; // stays in the loop, but goes to the next j immediately
      var key = m.haps[i] + "," + d.haps[j]; // distinct by value: [1,1] gives the same way twice
      if (!ways[key]) {
        ways[key] = true;
        n++;
        fromMom = m.haps[i];
        fromDad = d.haps[j];
      }
    }
  }
  if (n === 0) return null; // if n == 0, means that we did not find any possible ways of transmission.
  if (n > 1) return "ambiguous"; // if n > 1, means that there is more than one ways of transmission.
  var momAbnormal = fromMom !== m.normal,
    dadAbnormal = fromDad !== d.normal;
  if (momAbnormal && dadAbnormal) return "both";
  if (momAbnormal) return "maternal";
  if (dadAbnormal) return "paternal";
  return null; // both copies normal: the child isn't a carrier
}

/**
 * Parent of origin of one child's CNV: a transmission label for every child,
 * SNV-style. Entry point for the --sample-expr po_* lines. slivar calls it once
 * per sample and record (founders included: they have no parents, so they get
 * null), and values each po_* tag with the IDs of the samples whose result
 * matches. Affected status plays no part.
 *
 * First, which parent's transmitted copy is abnormal (get_transmission_origin,
 * the same computation as TRIO_TABLE's origin column): "maternal", "paternal",
 * "both" or "ambiguous". When no combination of the parents' copies explains
 * the child's CN, get_denovo_origin() decides between "denovo" (both parents
 * known and normal) and "possible_denovo" (one parent normal, the other
 * unknown). A parent carrying a CNV that still can't explain the child's
 * (a Mendelian error) gets null: unlike the SNV pipeline's po_denovo, that
 * isn't called de novo.
 *
 * @param {Sample}  s        the child
 * @param {Variant} variant
 * @returns {?string} "maternal" | "paternal" | "both" | "ambiguous" | "denovo"
 *   | "possible_denovo", or null when it does not apply (not a carrier, no CN,
 *   unmodellable record, founder) or can't be told
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
 * cn_po_origin(kid, variant)  // → null
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
  var origin = get_transmission_origin(s, variant);
  if (origin !== null) return origin;
  return get_denovo_origin(s, variant);
}

/**
 * For a carrier child whose CN no combination of the parents' copies explains:
 * is the CNV new? Both parents known and normal → "denovo"; one known and
 * normal, the other unknown (absent or without CN) → "possible_denovo". A parent
 * who passes nothing here (a father's X to a son, any mother's Y) counts as
 * known and normal, as in TRIO_TABLE, so this always agrees with the de_novo /
 * possible_de_novo calls. A child with no parent at all in the family gets null.
 * @param {Sample}  kid      a carrier child (cn_po_origin checked it)
 * @param {Variant} variant
 * @returns {?string} "denovo" | "possible_denovo" | null
 * @example
 * get_denovo_origin(<kid CN1, mom CN2, dad CN2>, <chr7>)          // → "denovo"
 * get_denovo_origin(<son CN0, mom CN2, no father>, <non-PAR chrX>) // → "denovo" (a father passes no X)
 * get_denovo_origin(<kid CN1, mom CN2, no father>, <chr7>)        // → "possible_denovo"
 * get_denovo_origin(<kid CN1, mom CN2, dad CN missing>, <chr7>)   // → "possible_denovo"
 * get_denovo_origin(<kid CN0, mom CN2, dad CN1>, <chr7>)          // → null (dad carries a CNV)
 * get_denovo_origin(<son CN2, no mother, dad CN1>, <non-PAR chrX>) // → "possible_denovo" (his X is his mother's)
 * get_denovo_origin(<singleton son CN0>, <non-PAR chrX>)           // → null (no parent in the family)
 */
function get_denovo_origin(kid, variant) {
  if (!has_parent(kid, "mom") && !has_parent(kid, "dad")) return null;
  var mom = get_parent_status(kid, "mom", variant),
    dad = get_parent_status(kid, "dad", variant);
  var known = function (st) {
    return st === "normal" || st === "nothing";
  };
  var unknown = function (st) {
    return st === "absent" || st === "unknown";
  };
  if (known(mom) && known(dad)) return "denovo";
  if ((known(mom) && unknown(dad)) || (known(dad) && unknown(mom))) return "possible_denovo";
  return null;
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
    get_transmission_origin: get_transmission_origin,
    get_expected_cn: get_expected_cn,
    get_cn_dev: get_cn_dev,
    has_cn: has_cn,
  };
