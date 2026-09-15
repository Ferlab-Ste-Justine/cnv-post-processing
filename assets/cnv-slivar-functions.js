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

function has_cn(sample) {
    return sample.CN >= 0
}

function cn_deviation(sample) {
    return sample.CN - 2
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
        var s = fam[i]
        if (s.affected && s.kids && s.kids.length > 0) return true
    }
    return false
}

// True iff no family member has children in the family (i.e. no parents are in the PED) --
// covers both literal solos and sibship-only families. Both provide zero segregation evidence
// for recessive (nothing to confirm carrier parents against).
function no_parents_in_fam(fam) {
    for (var i = 0; i < fam.length; i++) {
        if (fam[i].kids && fam[i].kids.length > 0) return false
    }
    return true
}

// True iff every affected family member has both parents present in the family. Recessive can
// still fire without this (see moi_cnv_recessive), but its low-confidence moi_cnv_candidate flag
// distinguishes a fully-confirmed biallelic origin from a plausible-but-unconfirmed one.
function all_affected_have_both_parents(fam) {
    for (var i = 0; i < fam.length; i++) {
        var s = fam[i]
        if (s.affected && (!("mom" in s) || !("dad" in s))) return false
    }
    return true
}

// ---------------------------------------------------------------------------
// Per-sample segregation predicates (called via fam.every(...))
// ---------------------------------------------------------------------------

// --- de novo ---
// Unaffected => confirmed non-carrier (CN==2).
// Affected   => a real carrier (CN!=2), AND both parents present in the family -- de novo isn't
//               provable otherwise (mirrors the SNV model's own self-gate to trios).
function cn_segregating_denovo(s) {
    if (!has_cn(s)) return false
    if (!s.affected) return s.CN == 2
    if (s.CN == 2) return false
    return ("mom" in s) && ("dad" in s)
}

// --- dominant ---
// Unaffected => confirmed non-carrier (CN==2).
// Affected   => a real carrier (CN!=2).
// Note: same per-sample shape as de novo -- what distinguishes the two is the family-level gate
// (moi_cnv_dominant requires has_aff_parent(fam), i.e. an affected parent transmitting the
// variant across a generation; a plain trio with only the kid affected never satisfies that gate,
// so it falls to de novo instead, exactly like the SNV model).
function cn_segregating_dominant(s) {
    if (!has_cn(s)) return false
    if (s.affected) return s.CN != 2
    return s.CN == 2
}

// --- recessive ---
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
    if (!has_cn(s)) return false
    if (s.affected) return cn_deviation(s) != 0
    if (s.kids && s.kids.length > 0) {
        for (var i = 0; i < s.kids.length; i++) {
            var kid = s.kids[i]
            if (!kid.affected) continue
            if (!has_cn(kid)) return false
            var kid_dev = cn_deviation(kid)
            if (kid_dev == 0) return false
            if (cn_deviation(s) * 2 != kid_dev) return false
        }
        return true
    }
    return s.CN == 2
}

// ---------------------------------------------------------------------------
// Mode-of-inheritance tag functions (one per --family-expr tag in conf/modules.config)
// ---------------------------------------------------------------------------

function moi_cnv_denovo(fam) { return fam.every(cn_segregating_denovo) }

function moi_cnv_dominant(fam) { return has_aff_parent(fam) && fam.every(cn_segregating_dominant) }

function moi_cnv_recessive(fam) { return !no_parents_in_fam(fam) && fam.every(cn_segregating_recessive) }

// Low-confidence companion to moi_cnv_recessive: fires alongside it when family structure isn't
// complete enough to fully confirm biallelic origin (not every affected member has both parents
// present -- e.g. a duo with only one parent in the PED).
function moi_cnv_candidate(fam) {
    return moi_cnv_recessive(fam) && !all_affected_have_both_parents(fam)
}

// True when there's a real, known-CN variant in an affected family member, but the pattern
// doesn't match de novo / dominant / recessive -- e.g. one parent's CN matches neither 2 nor the
// expected carrier deviation.
function moi_cnv_ambiguous(fam) {
    if (moi_cnv_unknown_cn(fam)) return false
    for (var i = 0; i < fam.length; i++) {
        var s = fam[i]
        if (s.affected && has_cn(s) && s.CN != 2) {
            if (moi_cnv_denovo(fam) || moi_cnv_dominant(fam) || moi_cnv_recessive(fam)) return false
            return true
        }
    }
    return false
}

// True when there isn't enough information to classify: an affected family member's own CN is
// missing, or an affected member carries a real variant but a family member needed to resolve its
// origin (parent, for de novo/dominant/recessive) has missing CN.
function moi_cnv_unknown_cn(fam) {
    var any_affected_carrier = false
    for (var i = 0; i < fam.length; i++) {
        var s = fam[i]
        if (!s.affected) continue
        if (!has_cn(s)) return true
        if (s.CN != 2) any_affected_carrier = true
    }
    if (!any_affected_carrier) return false
    for (var i = 0; i < fam.length; i++) {
        if (!has_cn(fam[i])) return true
    }
    return false
}
