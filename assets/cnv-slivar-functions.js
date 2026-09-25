// CNV mode-of-inheritance classification for `slivar expr --family-expr`.
//
// Mirrors the architecture of the SNV sibling pipeline's assets/slivar-functions.js as closely as
// CN-based genotypes allow: per-sample segregation predicates (segregating_*(s), each judging one
// family member's own CN against affected-status and slivar's own s.mom/s.dad/s.kids, which it
// populates from the PED) combined via fam.every(...), plus family-level structural gates
// (has_aff_parent, no_parents_in_fam) that decide whether the family's shape provides enough
// evidence for a given MoI to be assessed at all. This generalizes to duos/singletons/sibships the
// same way the SNV model does, not just trios.
//
// Differs from the SNV model in two ways specific to CN data (no direct analogs there):
//   - No hq1()-style AD/GQ/AB quality gate -- CN itself is already the output of depth-based
//     calling + refinement, so has_cn(s) (a resolved value) is the only per-sample quality check.
//   - Recessive needs its own CN-specific definition -- see cn_segregating_recessive below.
//
// Every tag function's INFO writeback happens through slivar's own tag:expr mechanism (each
// --family-expr becomes its own INFO/<tag>=<sample-list> field), confirmed to be the only
// annotation mechanism that persists here (see git history / PR discussion for the earlier
// --trio, single-string-label design this replaced, which relied on a manual
// variant.INFO[...] = ... write that does not actually persist).
//
// Missing FORMAT/CN reads back as htslib's raw BCF int32-missing sentinel in slivar's JS binding,
// not JS undefined/null -- confirmed directly. `< 0` catches that without hardcoding the
// sentinel's exact value; real copy numbers are never negative.
//
// Called like this:
// slivar expr \
//    --vcf variants.trio_NA12878_91_92_joint.cnv.vep.vcf.gz \
//    --ped trio_NA12878_91_92.ped \
//    --js cnv-slivar-functions.js \
//    --family-expr 'de_novo_candidate:moi_cnv_denovo(fam)' \
//    --family-expr 'dominant_inherited:moi_cnv_dominant(fam)' \
//    --family-expr 'recessive_candidate:moi_cnv_recessive(fam)' \
//    --family-expr 'candidate:moi_cnv_candidate(fam)' \
//    --family-expr 'ambiguous:moi_cnv_ambiguous(fam)' \
//    --family-expr 'unknown_cn:moi_cnv_unknown_cn(fam)' \
//    --family-expr 'po_mother:po_mother(fam)' \
//    --family-expr 'po_father:po_father(fam)' \
//    --family-expr 'po_ambiguous:po_ambiguous(fam)' \
//    -o trio_NA12878_91_92_joint.cnv.slivar.vcf.gz

function has_cn(sample) {
  return sample.CN >= 0;
}

function cn_deviation(sample) {
  return sample.CN - 2;
}

// ---------------------------------------------------------------------------
// Family-level structural gates (unchanged in spirit from the SNV model -- these only look at
// affected-status and parent/child relationships, nothing SNV-specific)
// ---------------------------------------------------------------------------

// True iff any family member is BOTH affected AND has children in the family (i.e. an affected
// parent we observe transmitting to the next generation). Gates `dominant` so it doesn't also
// fire on de novo configurations (kid affected+carrier, parents unaffected+non-carrier).
function has_aff_parent(fam) {
  for (var i = 0; i < fam.length; i++) {
    var s = fam[i];
    if (s.affected && s.kids && s.kids.length > 0) return true;
  }
  return false;
}

// True iff no family member has children in the family (i.e. no parents are in the PED) --
// covers both literal solos and sibship-only families. Both provide zero segregation evidence
// for recessive (nothing to confirm carrier parents against).
function no_parents_in_fam(fam) {
  for (var i = 0; i < fam.length; i++) {
    if (fam[i].kids && fam[i].kids.length > 0) return false;
  }
  return true;
}

// True iff every affected family member has both parents present in the family. Recessive can
// still fire without this (see moi_cnv_recessive), but its low-confidence moi_cnv_candidate flag
// distinguishes a fully-confirmed biallelic origin from a plausible-but-unconfirmed one.
function all_affected_have_both_parents(fam) {
  for (var i = 0; i < fam.length; i++) {
    var s = fam[i];
    if (s.affected && (!("mom" in s) || !("dad" in s))) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// Per-sample segregation predicates (called via fam.every(...))
// ---------------------------------------------------------------------------

// --- de novo ---
// Unaffected => confirmed non-carrier (CN==2).
// Affected   => a real carrier (CN!=2), AND both parents present in the family -- de novo isn't
//               provable otherwise (mirrors the SNV model's own self-gate to trios).
function cn_segregating_denovo(s) {
  if (!has_cn(s)) return false;
  if (!s.affected) return s.CN == 2;
  if (s.CN == 2) return false;
  return "mom" in s && "dad" in s;
}

// --- dominant ---
// dominant needs to see the trait passed down from an affected parent (vertical transmission, only one bad copy needed)
// Unaffected => confirmed non-carrier (CN==2).
// Affected   => a real carrier (CN!=2).
// Note: same per-sample shape as de novo -- what distinguishes the two is the family-level gate
// (moi_cnv_dominant requires has_aff_parent(fam), i.e. an affected parent transmitting the
// variant across a generation; a plain trio with only the kid affected never satisfies that gate,
// so it falls to de novo instead, exactly like the SNV model).
function cn_segregating_dominant(s) {
  if (!has_cn(s)) return false;
  if (s.affected) return s.CN != 2;
  return s.CN == 2;
}

// --- recessive ---
// recessive needs to see two unaffected parents each contributing half a "dose" that adds up to a full hit in the affected child (horizontal convergence, two bad copies needed)
// Affected => a real carrier (CN!=2) -- a "double dose" event (e.g. CN=0 full deletion, or a
//             symmetric CN=4 for a duplication scenario).
// Unaffected parent => carries exactly HALF of every affected child's deviation, same direction
//             (e.g. CN=1 when the affected child is CN=0) -- the CN analog of a heterozygous
//             carrier. Checked via s.kids, since (unlike SNV's absolute het/hom_ref states) what
//             counts as "carrier" here depends on the specific affected child's own CN.
// Unaffected non-parent (sibling) => conservatively required to be a confirmed non-carrier
//             (CN==2). The SNV model permits a sibling to be either a carrier or non-carrier
//             without checking a specific magnitude, because SNV het/hom_ref are absolute states;
//             a CNV sibling's "carrier" CN can only be judged against a specific affected child's
//             deviation the way a parent's can, which needs its own design if ever generalized
//             beyond trios/duos with a sibling -- deliberately conservative for this first pass.
function cn_segregating_recessive(s) {
  if (!has_cn(s)) return false;
  if (s.affected) return cn_deviation(s) != 0;
  if (s.kids && s.kids.length > 0) {
    for (var i = 0; i < s.kids.length; i++) {
      var kid = s.kids[i];
      if (!kid.affected) continue;
      if (!has_cn(kid)) return false;
      var kid_dev = cn_deviation(kid);
      if (kid_dev == 0) return false;
      if (cn_deviation(s) * 2 != kid_dev) return false;
    }
    return true;
  }
  return s.CN == 2;
}

// ---------------------------------------------------------------------------
// Mode-of-inheritance tag functions (one per --family-expr tag in conf/modules.config)
// ---------------------------------------------------------------------------

function moi_cnv_denovo(fam) {
  return fam.every(cn_segregating_denovo);
}

function moi_cnv_dominant(fam) {
  return has_aff_parent(fam) && fam.every(cn_segregating_dominant);
}

function moi_cnv_recessive(fam) {
  return !no_parents_in_fam(fam) && fam.every(cn_segregating_recessive);
}

// Low-confidence companion to moi_cnv_recessive: fires alongside it when family structure isn't
// complete enough to fully confirm biallelic origin (not every affected member has both parents
// present -- e.g. a duo with only one parent in the PED).
function moi_cnv_candidate(fam) {
  return moi_cnv_recessive(fam) && !all_affected_have_both_parents(fam);
}

// True when there's a real, known-CN variant in an affected family member, but the pattern
// doesn't match de novo / dominant / recessive -- e.g. one parent's CN matches neither 2 nor the
// expected carrier deviation.
function moi_cnv_ambiguous(fam) {
  if (moi_cnv_unknown_cn(fam)) return false;
  for (var i = 0; i < fam.length; i++) {
    var s = fam[i];
    if (s.affected && has_cn(s) && s.CN != 2) {
      if (moi_cnv_denovo(fam) || moi_cnv_dominant(fam) || moi_cnv_recessive(fam)) return false;
      return true;
    }
  }
  return false;
}

// True when there isn't enough information to classify: an affected family member's own CN is
// missing, or an affected member carries a real variant but a family member needed to resolve its
// origin (parent, for de novo/dominant/recessive) has missing CN.
function moi_cnv_unknown_cn(fam) {
  var any_affected_carrier = false;
  for (var i = 0; i < fam.length; i++) {
    var s = fam[i];
    if (!s.affected) continue;
    if (!has_cn(s)) return true;
    if (s.CN != 2) any_affected_carrier = true;
  }
  if (!any_affected_carrier) return false;
  for (var i = 0; i < fam.length; i++) {
    if (!has_cn(fam[i])) return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
// Confident parent-of-origin (autosomal CNVs only)
// ---------------------------------------------------------------------------
//
// For a kid with an abnormal CN (CN!=2), determines -- ONLY when the arithmetic is unambiguous --
// whether the variant haplotype can be confidently attributed to the maternal haplotype, the
// paternal haplotype, or both (biparental). This is deliberately narrow: most parent-of-origin
// inference from total CN alone is NOT resolvable (each observed CN is a SUM of two haplotype
// copy-numbers we never see directly), so the default is "not confident", and only the specific
// arithmetic pattern below earns a tag.
//
// Autosomal only -- out of scope: X/Y special-casing, de novo parent-of-origin (which chromosome
// a NEW mutation landed on), true UPD detection. moi_cnv_upd-style Mendelian-violation checks
// would be a separate function; this one assumes ordinary biparental transmission and only asks
// *which* parent's haplotype was transmitted.

// Decomposes a parent's own observed total CN into the haplotype pair it unambiguously implies,
// or null if it doesn't decompose unambiguously. CN=4 or higher is self-ambiguous at the
// single-person level (e.g. 4 could be [3,1] triplication+normal, or [2,2] double-duplication --
// indistinguishable from a single depth measurement), so only CN 0-3 are covered.
function cn_haplotypes(cn) {
  if (cn == 0) return [0, 0];
  if (cn == 1) return [1, 0];
  if (cn == 2) return [1, 1];
  if (cn == 3) return [2, 1];
  return null;
}

// Given the kid's own CN and both parents' haplotype pairs (see cn_haplotypes above), finds every
// mom/dad haplotype combination (kid gets exactly one haplotype from each parent) that reproduces
// the kid's observed CN, and resolves confidence from how many combinations match:
//   0 matches -> null       (inconsistent with Mendelian transmission -- de novo/bad-data
//                            territory already covered by moi_cnv_denovo/moi_cnv_ambiguous, not
//                            this function's job to flag)
//   1 match   -> the single combination pins down exactly which haplotype came from which parent:
//                "maternal" (only the transmitted maternal haplotype is abnormal), "paternal"
//                (only the transmitted paternal haplotype is abnormal), or "both" (both
//                transmitted haplotypes are abnormal -- a real, correct outcome, not an error)
//   2+ matches -> "ambiguous" (e.g. both parents are [2,1] simple-duplication carriers and the
//                kid is also CN=3 -- either parent could have transmitted the duplicated
//                haplotype)
//
// Combos are deduplicated by VALUE, not array position, before counting: a homozygous parent's
// own haplotype pair has two equal entries (CN=2 -> [1,1], CN=0 -> [0,0]), which would otherwise
// make the same biological explanation appear twice in the cartesian product below and falsely
// push a genuinely single-match case (e.g. mom CN=2 + dad CN=3 -> kid CN=3, the ordinary
// "one normal parent, one carrier parent" case) into the "ambiguous" branch.
function cn_po_resolve(kid_cn, mom_hap, dad_hap) {
  var combos = [
    [mom_hap[0], dad_hap[0]],
    [mom_hap[0], dad_hap[1]],
    [mom_hap[1], dad_hap[0]],
    [mom_hap[1], dad_hap[1]],
  ];
  var matches = [];
  var seen = {};
  for (var i = 0; i < combos.length; i++) {
    if (combos[i][0] + combos[i][1] != kid_cn) continue;
    var key = combos[i][0] + "," + combos[i][1];
    if (seen[key]) continue;
    seen[key] = true;
    matches.push(combos[i]);
  }
  if (matches.length == 0) return null;
  if (matches.length > 1) return "ambiguous";
  var transmitted_maternal = matches[0][0];
  var transmitted_paternal = matches[0][1];
  if (transmitted_maternal != 1 && transmitted_paternal == 1) return "maternal";
  if (transmitted_maternal == 1 && transmitted_paternal != 1) return "paternal";
  if (transmitted_maternal != 1 && transmitted_paternal != 1) return "both";
  // Both transmitted haplotypes normal -- no variant to attribute. Not normally reached since
  // the per-sample predicates below only call this on kids with CN!=2 to begin with.
  return null;
}

// Resolves a single sample's own parent-of-origin outcome, or null when this predicate doesn't
// apply to it at all (not a carrier kid, or missing/non-trio/uncallable parents) -- shared by the
// per-sample predicates and the family-level existential gates below so the same "does this even
// apply" logic isn't repeated four times.
function cn_po_origin(s) {
  if (!has_cn(s) || s.CN == 2) return null;
  if (!("mom" in s) || !("dad" in s)) return null;
  if (!has_cn(s.mom) || !has_cn(s.dad)) return null;
  var mom_hap = cn_haplotypes(s.mom.CN);
  var dad_hap = cn_haplotypes(s.dad.CN);
  if (!mom_hap || !dad_hap) return null;
  return cn_po_resolve(s.CN, mom_hap, dad_hap);
}

// Per-sample parent-of-origin predicates (called via fam.every(...), same composition as the MoI
// predicates above). Each only judges whether ITS OWN sample is consistent with the named origin
// -- a sample this predicate doesn't apply to (cn_po_origin(s) == null) vacuously passes true so
// it never blocks fam.every(...) for other, qualifying family members (e.g. siblings, or the
// parents themselves), matching how has_cn()/("mom" in s) gate other per-sample checks in this
// file. Because vacuous passes make fam.every(...) alone satisfiable by a family with no
// resolvable member at all, po_mother/po_father/po_ambiguous below additionally require
// fam.some(...) member to have actually resolved to that origin -- mirroring how moi_cnv_dominant
// and moi_cnv_recessive pair fam.every(...) with a family-level structural gate
// (has_aff_parent/no_parents_in_fam) rather than relying on fam.every(...) alone.
function cn_segregating_po_maternal(s) {
  var origin = cn_po_origin(s);
  return origin == null || origin == "maternal" || origin == "both";
}

function cn_segregating_po_paternal(s) {
  var origin = cn_po_origin(s);
  return origin == null || origin == "paternal" || origin == "both";
}

// Distinct from "not evaluated" (cn_segregating_po_maternal/paternal both vacuously true when
// this predicate's gates aren't met): this one is true only when the family WAS evaluable but the
// arithmetic genuinely couldn't resolve a single origin (2+ matching combinations), so a reviewer
// can tell "we looked, but couldn't resolve it" apart from "this case doesn't apply here at all".
function cn_segregating_po_ambiguous(s) {
  var origin = cn_po_origin(s);
  return origin == null || origin == "ambiguous";
}

function po_mother(fam) {
  return (
    fam.every(cn_segregating_po_maternal) &&
    fam.some(function (s) {
      var origin = cn_po_origin(s);
      return origin == "maternal" || origin == "both";
    })
  );
}

function po_father(fam) {
  return (
    fam.every(cn_segregating_po_paternal) &&
    fam.some(function (s) {
      var origin = cn_po_origin(s);
      return origin == "paternal" || origin == "both";
    })
  );
}

function po_ambiguous(fam) {
  return (
    fam.every(cn_segregating_po_ambiguous) &&
    fam.some(function (s) {
      return cn_po_origin(s) == "ambiguous";
    })
  );
}
