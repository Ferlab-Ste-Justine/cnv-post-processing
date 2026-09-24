#!/usr/bin/env node
//
// Unit tests for assets/cnv-slivar-functions.js -- the slivar `--js`/`--family-expr` classification
// functions run by the SLIVAR_EXPR module. This exercises the classification LOGIC directly, as
// plain sample objects (the same shape slivar's own JS engine hands each function via `fam`), so it
// runs instantly with plain `node` and needs neither Docker nor a real slivar binary.
//
// This deliberately does NOT replace an actual nf-test run of the SLIVAR_EXPR module (see
// modules/local/slivar_expr/tests/main.nf.test) -- that's the only thing that verifies slivar's own
// VCF/PED parsing, --family-expr wiring, and JS-engine (goja) compatibility end to end. This file's
// job is narrower and cheaper: catch a logic regression in the classification functions themselves
// before it ever reaches a container. It exists because exactly this kind of bug slipped through
// once already -- the first cut of the parent-of-origin arithmetic (cn_po_resolve) double-counted a
// homozygous parent's own symmetric haplotype pair and silently misclassified the single most common
// real case (one carrier parent + one reference parent) as "ambiguous" instead of resolving it. That
// exact scenario is test 7/8 below.
//
// Run with: node assets/tests/cnv_slivar_functions.test.js
// Wired into: scripts/run-test-suite.sh (runs before the Docker-dependent nf-test suite).

// No 'use strict' here, deliberately: the functions under test are loaded via a top-level direct
// eval() below, and strict-mode eval() runs in its own scope, so its function declarations
// wouldn't leak into this file's scope the way non-strict "direct eval" does.

const fs = require("fs");
const path = require("path");

const SRC_PATH = path.join(__dirname, "..", "cnv-slivar-functions.js");
const src = fs.readFileSync(SRC_PATH, "utf8");
// eslint-disable-next-line no-eval -- loading the real functions under test, not arbitrary input
eval(src);

// ---------------------------------------------------------------------------
// Toy family builders
// ---------------------------------------------------------------------------

function sample(cn, affected) {
  return { CN: cn, affected: !!affected };
}

// Links a trio the way slivar's own PED-derived `fam` objects are linked: each parent's `.kids`
// lists their children, each kid's `.mom`/`.dad` point back at the parent objects.
function trio(mom, dad, kid) {
  mom.kids = [kid];
  dad.kids = [kid];
  kid.mom = mom;
  kid.dad = dad;
  return [mom, dad, kid];
}

// A duo (one parent + kid) -- e.g. the other parent genuinely isn't in the pedigree at all, not
// just missing a CN call.
function duo(mom, kid) {
  mom.kids = [kid];
  kid.mom = mom;
  return [mom, kid];
}

const ALL_TAGS = [
  "de_novo_candidate",
  "dominant_inherited",
  "recessive_candidate",
  "candidate",
  "ambiguous",
  "unknown_cn",
  "parent_of_origin_maternal",
  "parent_of_origin_paternal",
  "parent_of_origin_ambiguous",
];

function evalTags(fam) {
  return {
    de_novo_candidate: moi_cnv_denovo(fam),
    dominant_inherited: moi_cnv_dominant(fam),
    recessive_candidate: moi_cnv_recessive(fam),
    candidate: moi_cnv_candidate(fam),
    ambiguous: moi_cnv_ambiguous(fam),
    unknown_cn: moi_cnv_unknown_cn(fam),
    parent_of_origin_maternal: moi_cnv_po_maternal(fam),
    parent_of_origin_paternal: moi_cnv_po_paternal(fam),
    parent_of_origin_ambiguous: moi_cnv_po_ambiguous(fam),
  };
}

// ---------------------------------------------------------------------------
// Toy scenarios -- one real-world-shaped family per tag (or tag combination), each independently
// hand-traced against the real function bodies and cross-checked here. `expectTrue` lists every tag
// that must be true; every other tag in ALL_TAGS is asserted false, so a scenario also documents
// (and guards) which OTHER tags a real reviewer would expect to stay quiet.
// ---------------------------------------------------------------------------

const SCENARIOS = [
  {
    name: "1. De novo: unaffected CN=2 parents, affected CN=3 kid",
    fam: trio(sample(2, false), sample(2, false), sample(3, true)),
    expectTrue: ["de_novo_candidate"],
  },
  {
    name: "2. Dominant, inherited from mom: affected CN=3 mom, unaffected CN=2 dad, affected CN=3 kid",
    fam: trio(sample(3, true), sample(2, false), sample(3, true)),
    // A dominant-inherited CNV is *also* a resolvable parent-of-origin case -- these two tags
    // describe the same real event from two different angles, not a conflict.
    expectTrue: ["dominant_inherited", "parent_of_origin_maternal"],
  },
  {
    name: "3. Recessive (complete trio): unaffected CN=1 carrier parents, affected CN=0 kid",
    fam: trio(sample(1, false), sample(1, false), sample(0, true)),
    // Both parents independently transmitted their deleted haplotype -- a real, correct
    // "both" parent-of-origin outcome alongside recessive segregation.
    expectTrue: ["recessive_candidate", "parent_of_origin_maternal", "parent_of_origin_paternal"],
  },
  {
    name: "4. Candidate (recessive pattern, incomplete trio): unaffected CN=1 carrier mom, affected CN=0 kid, no dad in the pedigree at all",
    fam: duo(sample(1, false), sample(0, true)),
    // moi_cnv_candidate is defined as moi_cnv_recessive(fam) && !all_affected_have_both_parents(fam)
    // -- it can only ever fire alongside recessive_candidate, never instead of it.
    expectTrue: ["recessive_candidate", "candidate"],
  },
  {
    name: "5. Ambiguous: affected CN=0 kid, but neither parent's CN matches a de novo, dominant, or recessive pattern",
    fam: trio(sample(2, false), sample(1, false), sample(0, true)),
    expectTrue: ["ambiguous"],
  },
  {
    name: "6. Unknown CN: affected kid has no callable CN at all",
    fam: trio(sample(2, false), sample(2, false), sample(-1, true)),
    expectTrue: ["unknown_cn"],
  },
  {
    name: "7. Parent-of-origin, maternal only (all unaffected -- PO is independent of affected status): CN=1 carrier mom, CN=2 dad, CN=1 kid",
    fam: trio(sample(1, false), sample(2, false), sample(1, false)),
    expectTrue: ["parent_of_origin_maternal"],
  },
  {
    name: "8. Parent-of-origin, paternal only: CN=2 mom, CN=1 carrier dad, CN=1 kid",
    fam: trio(sample(2, false), sample(1, false), sample(1, false)),
    expectTrue: ["parent_of_origin_paternal"],
  },
  {
    name: '9. Parent-of-origin, genuinely ambiguous: both parents CN=3 carriers, kid also CN=3 (the tied "middle" case -- either parent could have transmitted)',
    fam: trio(sample(3, false), sample(3, false), sample(3, false)),
    expectTrue: ["parent_of_origin_ambiguous"],
  },
];

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

let failures = 0;

for (const scenario of SCENARIOS) {
  const actual = evalTags(scenario.fam);
  const wantTrue = new Set(scenario.expectTrue);
  const mismatches = [];
  for (const tag of ALL_TAGS) {
    const want = wantTrue.has(tag);
    if (actual[tag] !== want) {
      mismatches.push(`${tag}: expected ${want}, got ${actual[tag]}`);
    }
  }
  if (mismatches.length === 0) {
    console.log(`PASS  ${scenario.name}`);
  } else {
    failures++;
    console.log(`FAIL  ${scenario.name}`);
    for (const m of mismatches) console.log(`        ${m}`);
  }
}

console.log("");
if (failures === 0) {
  console.log(`All ${SCENARIOS.length} scenarios passed.`);
  process.exit(0);
} else {
  console.log(`${failures} of ${SCENARIOS.length} scenario(s) failed.`);
  process.exit(1);
}
