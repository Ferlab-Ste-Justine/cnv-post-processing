// node assets/cnv-slivar-moi-harness.js            run all checks
// node assets/cnv-slivar-moi-harness.js --accept   ... and save the results as the new reference
// (from any directory: every file it reads or writes lives next to this script)
//
// Test bench for cnv-slivar-functions.js (the SLIVAR_EXPR rules) only.
// Builds made-up families shaped like slivar's sample objects and writes:
//   moi_scenarios.tsv       one row per named scenario: the family (region, each member's CN),
//                           MoI call per child, parent of origin
//   moi_scenarios.md        the same table as Markdown
//   trio_lookup_table.tsv   TRIO_TABLE expanded to concrete autosomal CNs (kid, mom, dad -> call)
//
// Reference ("golden file"): moi_scenarios.expected.tsv
//   Every run compares moi_scenarios.tsv with it and lists the scenarios whose
//   results changed. Each change is a behaviour change you should be able to
//   explain. When you agree with the new results, re-run with --accept.
//
// Random families: 200000 random pedigrees are run to check that nothing
// crashes, that cnv_moi() only ever returns a known call (or null), that every
// inherited_* and de novo call agrees with the child's parent of origin, and that has_moi_call()
// (the --family-expr entry point) tags the family with exactly its children's calls.
//
// Wiring: conf/modules.config must have one slivar expression per call and parent of origin.
//
// Exit code 1 if anything changed versus the reference or a check failed.
var T = require("./cnv-slivar-functions.js");
var fs = require("fs");
var path = require("path");
function here(f) {
  return path.join(__dirname, f);
} // all files live next to this script, whatever the cwd
var MISSING = -2147483648; // what slivar passes for FORMAT/CN "." (htslib's bcf_int32_missing, INT32_MIN)
var EXPECTED = "moi_scenarios.expected.tsv";
var ACCEPT = process.argv.indexOf("--accept") >= 0;
var failed = false;

// ---------------------------------------------------------------------------
// Fake slivar objects
// ---------------------------------------------------------------------------
// slivar's samples and variants THROW when a rule reads an attribute slivar didn't
// set ("unknown attribute:mom", "unknown attribute:END"), where a plain JS object
// would return undefined; slivar then writes the record with no tag. strict()
// makes the fakes behave the same way, so such a read fails here, not only in slivar.
function strict(obj, what) {
  return new Proxy(obj, {
    get: function (target, prop) {
      if (typeof prop === "string" && !(prop in target))
        throw new Error("unknown attribute:" + prop + " (" + what + ")");
      return target[prop];
    },
  });
}
function mk(id, cn, affected, sex) {
  return strict(
    { id: id, CN: cn === null ? MISSING : cn, affected: !!affected, sex: sex || "female", kids: [] },
    "sample " + id,
  );
}
function link(kid, mom, dad) {
  if (mom) {
    kid.mom = mom;
    mom.kids.push(kid);
  }
  if (dad) {
    kid.dad = dad;
    dad.kids.push(kid);
  }
}
// spec: { kid: [CN, affected], kidsex, mom: [CN, affected] | null, dad: [CN, affected] | null,
//         sibs: [[CN, affected, sex], ...] }     CN null = missing
function family(spec) {
  var fam = [];
  var mom = spec.mom ? mk("mom", spec.mom[0], spec.mom[1]) : null;
  var dad = spec.dad ? mk("dad", spec.dad[0], spec.dad[1], "male") : null;
  var kid = mk("kid", spec.kid[0], spec.kid[1], spec.kidsex);
  link(kid, mom, dad);
  fam.push(kid);
  (spec.sibs || []).forEach(function (sb, i) {
    var s = mk("sib" + (i + 1), sb[0], sb[1], sb[2]);
    link(s, mom, dad);
    fam.push(s);
  });
  if (mom) fam.push(mom);
  if (dad) fam.push(dad);
  return fam;
}
// start is 0-based, stop the 1-based END, as slivar sets them.
function variant(chrom, start, end) {
  return strict(
    { CHROM: chrom, POS: start + 1, start: start, stop: end, INFO: strict({ END: end, SVTYPE: "DEL" }, "INFO") },
    "variant",
  );
}
var AUTO = variant("chr7", 1000000, 1200000);
var XNP = variant("chrX", 50000000, 50200000);
var XSTRADDLE = variant("chrX", 2779999, 2790000); // crosses the PAR1 end (2,781,479)
var XPAR = variant("chrX", 999999, 1010000); // entirely inside PAR1
var YNP = variant("chrY", 10660601, 10693753);
var MT = variant("chrM", 99, 1000);
// A record with no INFO/END (e.g. the slivar module-test fixtures): slivar sets stop
// to POS, and reading INFO.END throws.
var AUTO_NO_END = strict(
  { CHROM: "chr7", POS: 1000001, start: 1000000, stop: 1000001, INFO: strict({ SVTYPE: "DEL" }, "INFO") },
  "variant",
);

// Each child's MoI call (cnv_moi), as id=call. slivar tags the family with each of these calls
// (has_moi_call, --family-expr), valued with the family's affected members.
function moi(fam, v) {
  return (
    fam
      .map(function (s) {
        var call = T.cnv_moi(s, v);
        return call === null ? null : s.id + "=" + call;
      })
      .filter(function (x) {
        return x !== null;
      })
      .join(" ") || "-"
  );
}
function po(fam, v) {
  return (
    fam
      .filter(function (s) {
        return "mom" in s || "dad" in s;
      })
      .map(function (s) {
        return s.id + "=" + T.cn_po_origin(s, v);
      })
      .join(" ") || "-"
  );
}

// ---------------------------------------------------------------------------
// Named scenarios — add new ones at the end of their group
// Prefixes: T trio, X chrX, Y chrY, M chrM, D duo, S singleton, Q two children
// ---------------------------------------------------------------------------
var scenarios = [
  ["T1  de novo het deletion (trio)", AUTO, { kid: [1, 1], mom: [2, 0], dad: [2, 0] }],
  ["T2  de novo duplication (trio)", AUTO, { kid: [3, 1], mom: [2, 0], dad: [2, 0] }],
  ["T3  de novo CN=0 (implausible)", AUTO, { kid: [0, 1], mom: [2, 0], dad: [2, 0] }],
  ["T4  het del inherited from UNAFFECTED mom", AUTO, { kid: [1, 1], mom: [1, 0], dad: [2, 0] }],
  ["T5  het del inherited from AFFECTED dad", AUTO, { kid: [1, 1], mom: [2, 0], dad: [1, 1] }],
  ["T6  biallelic deletion (CN=0), carrier parents", AUTO, { kid: [0, 1], mom: [1, 0], dad: [1, 0] }],
  ["T7  dup inherited from unaffected mom", AUTO, { kid: [3, 1], mom: [3, 0], dad: [2, 0] }],
  ["T8  dup, both parents CN=3", AUTO, { kid: [3, 1], mom: [3, 0], dad: [3, 0] }],
  ["T9  CN=4 kid, both parents CN=3", AUTO, { kid: [4, 1], mom: [3, 0], dad: [3, 0] }],
  ["T10 het del kid, both parents CN=1", AUTO, { kid: [1, 1], mom: [1, 0], dad: [1, 0] }],
  ["T11 dad CN missing", AUTO, { kid: [1, 1], mom: [2, 0], dad: [null, 0] }],
  ["T12 kid CN=2, mom CN=1", AUTO, { kid: [2, 1], mom: [1, 0], dad: [2, 0] }],
  ["T13 unaffected kid CN=1", AUTO, { kid: [1, 0], mom: [1, 0], dad: [2, 0] }],
  ["T14 kid CN=0, mom CN=0, dad CN=2 (impossible)", AUTO, { kid: [0, 1], mom: [0, 0], dad: [2, 0] }],
  ["T15 AFFECTED dad CN=1, affected kid CN=2 (not segregating)", AUTO, { kid: [2, 1], mom: [2, 0], dad: [1, 1] }],
  ["T16 kid CN=3, mom CN=4, dad CN=3 (possible, origin unknown)", AUTO, { kid: [3, 1], mom: [4, 0], dad: [3, 0] }],
  ["T17 kid CN=3, mom CN=4, dad CN=1 (possible: maternal, or both)", AUTO, { kid: [3, 1], mom: [4, 0], dad: [1, 0] }],
  // T18-T22: the scenarios of the retired assets/tests/cnv_slivar_functions.test.js not covered above
  ["T18 dup inherited from AFFECTED mom", AUTO, { kid: [3, 1], mom: [3, 1], dad: [2, 0] }],
  ["T19 kid CN=0, mom CN=2, dad CN=1 (impossible)", AUTO, { kid: [0, 1], mom: [2, 0], dad: [1, 0] }],
  ["T20 affected kid CN missing", AUTO, { kid: [null, 1], mom: [2, 0], dad: [2, 0] }],
  ["T21 unaffected kid CN=1, dad CN=1 (paternal origin)", AUTO, { kid: [1, 0], mom: [2, 0], dad: [1, 0] }],
  ["T22 all unaffected, all CN=3 (origin ambiguous)", AUTO, { kid: [3, 0], mom: [3, 0], dad: [3, 0] }],
  // T23-T26: a reviewer's cases. The origin comes from every way the parents' copies can add up to
  // the child's CN, not from comparing their deviations, so a parent with a different CNV is fine.
  [
    "T23 kid CN=1, mom CN=0, dad CN=2 (mom can only pass a deleted copy)",
    AUTO,
    { kid: [1, 1], mom: [0, 0], dad: [2, 0] },
  ],
  ["T24 kid CN=3, mom CN=3, dad CN=1 (dad passed his normal copy)", AUTO, { kid: [3, 1], mom: [3, 0], dad: [1, 0] }],
  // At most 2 + 1: only a hidden [3,0] split of mom's CN 3 could give 4 (MENDEL_ALLOW_RARE_SPLITS).
  ["T25 kid CN=4, mom CN=3, dad CN=2 (impossible)", AUTO, { kid: [4, 1], mom: [3, 0], dad: [2, 0] }],
  ["T26 de novo CN=4 (implausible, like T3)", AUTO, { kid: [4, 1], mom: [2, 0], dad: [2, 0] }],
  ["X1  chrX MALE kid CN=1, dad CN=1 (normal)", XNP, { kid: [1, 1], kidsex: "male", mom: [2, 0], dad: [1, 0] }],
  [
    "X2  chrX male kid CN=0, carrier mom CN=1, dad CN=1",
    XNP,
    { kid: [0, 1], kidsex: "male", mom: [1, 0], dad: [1, 0] },
  ],
  [
    "X3  chrX male kid CN=0, mom CN=2, dad CN=1 (de novo)",
    XNP,
    { kid: [0, 1], kidsex: "male", mom: [2, 0], dad: [1, 0] },
  ],
  // A record crossing a PAR boundary has no single normal CN (a man has 2 copies in the PAR
  // part, 1 outside it), so it isn't modelled: no tag, no parent of origin, for anyone.
  [
    "X4  chrX record straddling the PAR1 boundary, son CN=3, parents CN=2 (not modelled)",
    XSTRADDLE,
    { kid: [3, 1], kidsex: "male", mom: [2, 0], dad: [2, 0] },
  ],
  [
    "X5  chrX record straddling the PAR1 boundary, daughter CN=1, mom CN=1, dad CN=2 (not modelled)",
    XSTRADDLE,
    { kid: [1, 1], kidsex: "female", mom: [1, 0], dad: [2, 0] },
  ],
  [
    "X6  chrX PAR1, son CN=3, parents CN=2 (PAR: autosomal rules)",
    XPAR,
    { kid: [3, 1], kidsex: "male", mom: [2, 0], dad: [2, 0] },
  ],
  // Daughters on non-PAR chrX: her father's single X always goes to her.
  [
    "X7  chrX daughter CN=1, mom CN=2, dad CN=1 (de novo)",
    XNP,
    { kid: [1, 1], kidsex: "female", mom: [2, 0], dad: [1, 0] },
  ],
  [
    "X8  chrX daughter CN=3, mom CN=2, dad CN=2 (father's duplicated X)",
    XNP,
    { kid: [3, 1], kidsex: "female", mom: [2, 0], dad: [2, 0] },
  ],
  [
    "X9  chrX daughter CN=1, mom CN=1, dad CN=1 (mother's deletion)",
    XNP,
    { kid: [1, 1], kidsex: "female", mom: [1, 0], dad: [1, 0] },
  ],
  [
    "X10 chrX daughter CN=0, carrier mom CN=1, dad CN=0 (both X copies hit)",
    XNP,
    { kid: [0, 1], kidsex: "female", mom: [1, 0], dad: [0, 0] },
  ],
  // chrY: a daughter has no Y, so she isn't modelled; a son's Y is his father's, and his
  // mother (no CN there) passes nothing.
  [
    "Y1  chrY daughter, no CN on Y (not modelled)",
    YNP,
    { kid: [null, 1], kidsex: "female", mom: [null, 0], dad: [2, 0] },
  ],
  [
    "Y2  chrY son CN=2, dad CN=2 (father's duplicated Y)",
    YNP,
    { kid: [2, 1], kidsex: "male", mom: [null, 0], dad: [2, 0] },
  ],
  [
    "Y3  chrY son CN=2, dad CN=1 (de novo duplication)",
    YNP,
    { kid: [2, 1], kidsex: "male", mom: [null, 0], dad: [1, 0] },
  ],
  ["Y4  chrY son CN=0, dad CN=1 (de novo deletion)", YNP, { kid: [0, 1], kidsex: "male", mom: [null, 0], dad: [1, 0] }],
  [
    "Y5  chrY son CN=0, dad CN=0 (father's deleted Y)",
    YNP,
    { kid: [0, 1], kidsex: "male", mom: [null, 0], dad: [0, 0] },
  ],
  ["M1  chrM son CN=0, parents CN=2 (not modelled)", MT, { kid: [0, 1], kidsex: "male", mom: [2, 0], dad: [2, 0] }],
  ["D1  duo kid CN=1, mom CN=1", AUTO, { kid: [1, 1], mom: [1, 0], dad: null }],
  ["D2  duo kid CN=0, mom CN=1", AUTO, { kid: [0, 1], mom: [1, 0], dad: null }],
  ["D3  duo kid CN=1, mom CN=2", AUTO, { kid: [1, 1], mom: [2, 0], dad: null }],
  ["D4  duo kid CN=1, AFFECTED dad CN=1", AUTO, { kid: [1, 1], mom: null, dad: [1, 1] }],
  ["S1  singleton kid CN=1", AUTO, { kid: [1, 1], mom: null, dad: null }],
  [
    "Q1  quad: recessive + UNAFFECTED CARRIER sib CN=1",
    AUTO,
    { kid: [0, 1], mom: [1, 0], dad: [1, 0], sibs: [[1, 0]] },
  ],
  ["Q2  quad: recessive + unaffected sib CN=2", AUTO, { kid: [0, 1], mom: [1, 0], dad: [1, 0], sibs: [[2, 0]] }],
  [
    "Q3  quad: kid CN=1 from AFFECTED dad CN=1, unaffected sib CN=1",
    AUTO,
    { kid: [1, 1], mom: [2, 0], dad: [1, 1], sibs: [[1, 0]] },
  ],
  ["Q4  quad: de novo, unaffected sib CN=2", AUTO, { kid: [1, 1], mom: [2, 0], dad: [2, 0], sibs: [[2, 0]] }],
  ["Q5  quad: two affected sibs CN=1, parents CN=2", AUTO, { kid: [1, 1], mom: [2, 0], dad: [2, 0], sibs: [[1, 1]] }],
  [
    "Q8  quad: AFFECTED kid CN=2, unaffected sib CN=1 from mom",
    AUTO,
    { kid: [2, 1], mom: [1, 0], dad: [2, 0], sibs: [[1, 0]] },
  ],
  ["Q9  quad: mom 1 dad 3; kid 1 (mat), sib 3 (pat)", AUTO, { kid: [1, 1], mom: [1, 0], dad: [3, 0], sibs: [[3, 1]] }],
  [
    "Q10 quad: two affected sibs CN=1, both from AFFECTED dad CN=1",
    AUTO,
    { kid: [1, 1], mom: [2, 0], dad: [1, 1], sibs: [[1, 1]] },
  ],
];

// ---------------------------------------------------------------------------
// 1. Scenario table (TSV + Markdown)
// ---------------------------------------------------------------------------
// The family columns are built from each scenario's own spec, so they always show what was tested.
var REGIONS = [
  [AUTO, "chr7"],
  [XNP, "chrX non-PAR"],
  [XSTRADDLE, "chrX across the PAR1 end"],
  [XPAR, "chrX PAR1"],
  [YNP, "chrY non-PAR"],
  [MT, "chrM"],
  [AUTO_NO_END, "chr7, no INFO/END"],
];
function region(v) {
  for (var i = 0; i < REGIONS.length; i++) if (REGIONS[i][0] === v) return REGIONS[i][1];
  return v.CHROM;
}
// One member's cell: CN ("." = missing), " aff" when affected, "son " / "daughter " when the
// spec sets the sex; "absent" = not in the PED. E.g. "son 0 aff", "1", ".", "absent".
function member(m, sex) {
  if (!m) return "absent";
  var who = sex === "male" ? "son " : sex === "female" ? "daughter " : "";
  return who + (m[0] === null ? "." : m[0]) + (m[1] ? " aff" : "");
}
var header = [
  "scenario",
  "region",
  "kid",
  "mom",
  "dad",
  "sib",
  "MoI call per child (cnv_moi)",
  "parent of origin (sample-expr)",
];
var MOI = header.indexOf("MoI call per child (cnv_moi)");
var rows = scenarios.map(function (sc) {
  var spec = sc[2],
    fam = family(spec),
    v = sc[1];
  return [
    sc[0],
    region(v),
    member(spec.kid, spec.kidsex),
    member(spec.mom),
    member(spec.dad),
    (spec.sibs || [])
      .map(function (sb) {
        return member(sb, sb[2]);
      })
      .join(", "),
    moi(fam, v),
    po(fam, v),
  ];
});
function to_tsv(rs) {
  return (
    [header]
      .concat(rs)
      .map(function (r) {
        return r.join("\t");
      })
      .join("\n") + "\n"
  );
}
fs.writeFileSync(here("moi_scenarios.tsv"), to_tsv(rows));

function md_cell(c) {
  return String(c)
    .replace(/\|/g, "\\|")
    .replace(/\s{2,}/g, " ");
}
var md = [
  "# CNV MoI scenarios",
  "",
  "Generated by `node cnv-slivar-moi-harness.js` — do not edit by hand.",
  "",
  "Scenario prefixes: **T** trio, **X** chrX, **Y** chrY, **M** chrM, **D** duo, **S** singleton, **Q** two children. `-` = no tag.",
  "",
  "Family columns: each member's CN (`.` = missing), `aff` = affected, `absent` = not in the PED; `son` / `daughter` where the scenario sets the sex. " +
    "Normal CN: 2 on autosomes and in the PARs; on non-PAR chrX 1 for a son and 2 for a daughter; on non-PAR chrY 1 for a son (a daughter has no Y).",
  "",
  "| " + header.join(" | ") + " |",
  "|" +
    header
      .map(function () {
        return "---";
      })
      .join("|") +
    "|",
];
rows.forEach(function (r) {
  var c = r.map(md_cell);
  c[0] = "**" + c[0].replace(/^(\S+)\s+/, "$1** ");
  c[MOI] = c[MOI].replace(/=(\S+)/g, "=`$1`");
  md.push("| " + c.join(" | ") + " |");
});
fs.writeFileSync(here("moi_scenarios.md"), md.join("\n") + "\n");

// ---------------------------------------------------------------------------
// 2. Compare with the reference (moi_scenarios.expected.tsv)
// ---------------------------------------------------------------------------
function by_id(tsv) {
  // { "T4": row, ... } keyed on the scenario ID
  var out = {};
  tsv
    .trim()
    .split("\n")
    .slice(1)
    .forEach(function (line) {
      var r = line.split("\t");
      out[r[0].split(/\s+/)[0]] = r;
    });
  return out;
}
console.log("Compared with " + EXPECTED + ":");
if (!fs.existsSync(here(EXPECTED))) {
  console.log("  no reference yet — run with --accept to create it from these results");
} else {
  var exp = by_id(fs.readFileSync(here(EXPECTED), "utf8")),
    now = by_id(to_tsv(rows));
  var changes = [];
  Object.keys(now).forEach(function (id) {
    if (!exp[id]) {
      changes.push("  + " + id + " (new scenario): " + now[id][MOI]);
      return;
    }
    for (var c = 1; c < header.length; c++) {
      if (exp[id][c] !== now[id][c])
        changes.push("  ~ " + id + " " + header[c] + ":  " + exp[id][c] + "  ->  " + now[id][c]);
    }
  });
  Object.keys(exp).forEach(function (id) {
    if (!now[id]) changes.push("  - " + id + " (scenario removed)");
  });
  if (changes.length === 0) console.log("  no changes");
  else {
    console.log(changes.join("\n"));
    if (!ACCEPT) failed = true;
  }
}
if (ACCEPT) {
  fs.writeFileSync(here(EXPECTED), to_tsv(rows));
  console.log("  --accept: saved these results as the new " + EXPECTED);
}

// ---------------------------------------------------------------------------
// 3. TRIO_TABLE expanded to concrete autosomal CNs
// ---------------------------------------------------------------------------
var CNS = [0, 1, 3, 4],
  PCNS = [0, 1, 2, 3, 4, "absent"]; // kid CN 2 is not a carrier
var look = [["kid_CN", "mom_CN", "dad_CN", "child_call", "parent_of_origin"]];
CNS.forEach(function (k) {
  PCNS.forEach(function (m) {
    PCNS.forEach(function (d) {
      var fam = family({ kid: [k, 1], mom: m === "absent" ? null : [m, 0], dad: d === "absent" ? null : [d, 0] });
      look.push([k, m, d, T.trio_call(fam[0], AUTO), String(T.cn_po_origin(fam[0], AUTO))]);
    });
  });
});
fs.writeFileSync(
  here("trio_lookup_table.tsv"),
  look
    .map(function (r) {
      return r.join("\t");
    })
    .join("\n") + "\n",
);

// ---------------------------------------------------------------------------
// 4. Random families: no crash, only known tags, MoI and parent of origin agree
// ---------------------------------------------------------------------------
var KNOWN = { unknown_cn: true, inheritance_unresolved: true, mendelian_inconsistent: true };
T.TRIO_TABLE.forEach(function (row) {
  if (row[4] !== "no_pattern") KNOWN[row[4]] = true;
});
["male", "female"].forEach(function (sex) {
  Object.keys(T.X_LINKED_CALLS[sex]).forEach(function (call) {
    KNOWN[T.X_LINKED_CALLS[sex][call]] = true;
  });
});
var KNOWN_PO = {
  null: true,
  maternal: true,
  paternal: true,
  both: true,
  ambiguous: true,
  denovo: true,
  possible_denovo: true,
  possible_maternal: true,
  possible_paternal: true,
};

// Fixed seed: same families every run. mulberry32, in exact 32-bit integer maths:
// a plain (seed * 1103515245 + 12345) LCG overflows JS's 2^53 exact range, so its
// low bits -- the ones seed % n reads -- were always 0 (rnd(2) never returned 1).
var seed = 42;
function rnd(n) {
  seed = (seed + 0x6d2b79f5) | 0;
  var t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) % n;
}
function pick(a) {
  return a[rnd(a.length)];
}
var CNPOOL = [0, 1, 1, 2, 2, 2, 2, 3, 3, 4, null];
var N = 200000,
  counts = {},
  po_counts = {},
  problems = {};
function problem(msg, spec, v) {
  if (!problems[msg]) problems[msg] = v.CHROM + " " + JSON.stringify(spec);
}

for (var it = 0; it < N; it++) {
  var v = pick([AUTO, AUTO, XNP, XNP, YNP, XSTRADDLE, AUTO_NO_END]);
  var spec = {
    kid: [pick(CNPOOL), 1],
    kidsex: pick(["male", "female"]),
    mom: rnd(5) ? [pick(CNPOOL), rnd(4) === 0 ? 1 : 0] : null,
    dad: rnd(5) ? [pick(CNPOOL), rnd(4) === 0 ? 1 : 0] : null,
    sibs: rnd(2) ? [[pick(CNPOOL), rnd(2), pick(["male", "female"])]] : [],
  };
  var fam = family(spec);
  try {
    fam.forEach(function (s) {
      var po = T.cn_po_origin(s, v);
      if (!KNOWN_PO[po]) problem("unknown parent of origin: " + po, spec, v);
      var call = T.cnv_moi(s, v);
      if (!("mom" in s || "dad" in s)) {
        if (call !== null) problem("a founder got a call: " + call, spec, v);
        return;
      }
      if (call !== null && !KNOWN[call]) problem("unknown call: " + call, spec, v);
      counts[call || "-"] = (counts[call || "-"] || 0) + 1;
      var d = T.get_cn_dev(s, v);
      if (!(s.affected && d !== null && d !== 0)) return;
      if (typeof call !== "string") {
        // every affected carrier child must get a call ...
        problem("an affected carrier child got no call", spec, v);
        return;
      }
      if (call === "inherited_maternal" && po !== "maternal")
        problem("inherited_maternal but parent of origin " + po, spec, v);
      if (call === "inherited_paternal" && po !== "paternal")
        problem("inherited_paternal but parent of origin " + po, spec, v);
      var denovo_call = /denovo/.test(call) && !/possible/.test(call); // denovo, denovo_suspect, x_linked_*_denovo
      if (call === "x_linked_recessive" && po !== (s.sex === "male" ? "maternal" : "both"))
        problem("x_linked_recessive (" + s.sex + ") but parent of origin " + po, spec, v);
      if (denovo_call !== (po === "denovo")) problem(call + " but parent of origin " + po, spec, v);
      if (call === "possible_denovo" && po !== "possible_denovo")
        problem("possible_denovo but parent of origin " + po, spec, v);
      // a duo call names the known parent; its parent of origin must name the same one
      var side = !("mom" in s) || !T.has_cn(s.mom) ? "paternal" : "maternal";
      if (/^possible_inherited_|^recessive_lowconf$/.test(call) && po !== side && po !== "possible_" + side)
        problem(call + " but parent of origin " + po, spec, v);
      po_counts[v.CHROM + " " + po] = (po_counts[v.CHROM + " " + po] || 0) + 1; // ... and agree with its parent of origin
    });
    // the --family-expr entry point tags the family with exactly its children's calls
    var calls = fam.map(function (s) {
      return T.cnv_moi(s, v);
    });
    Object.keys(KNOWN).forEach(function (c) {
      if (T.has_moi_call(fam, v, c) !== calls.indexOf(c) >= 0)
        problem("has_moi_call(" + c + ") disagrees with the children's calls", spec, v);
    });
  } catch (e) {
    problem("crash: " + e.message, spec, v);
    continue;
  }
}

console.log("\nRandom families (" + N + "):");
var probs = Object.keys(problems);
if (probs.length === 0)
  console.log(
    "  no crashes, no unknown calls, every affected carrier child got a call that agrees with its parent of origin, has_moi_call() matches the calls",
  );
else {
  failed = true;
  probs.forEach(function (p) {
    console.log("  PROBLEM " + p + "\n    e.g. " + problems[p]);
  });
}
console.log('  MoI calls seen, per child ("-" = no call):');
Object.keys(counts)
  .sort(function (a, b) {
    return counts[b] - counts[a];
  })
  .forEach(function (t) {
    console.log("    " + String(counts[t]).padStart(7) + "  " + t);
  });
console.log("  carrier children's parent of origin, by chromosome:");
Object.keys(po_counts)
  .sort()
  .forEach(function (t) {
    console.log("    " + String(po_counts[t]).padStart(7) + "  " + t);
  });

// ---------------------------------------------------------------------------
// 5. The slivar wiring: conf/modules.config has an expression for every call
// ---------------------------------------------------------------------------
// A call cnv_moi() can return but SLIVAR_EXPR never asks for is silently never written.
var config = fs.readFileSync(here("../conf/modules.config"), "utf8");
var wired = { moi: {}, po: {} };
config.replace(/--family-expr '(\w+):has_moi_call\(fam, variant, \\"(\w+)\\"\)'/g, function (m, tag, call) {
  if (tag !== call) problem_wiring("tag " + tag + " tests the call " + call);
  wired.moi[call] = true;
});
config.replace(/--sample-expr '(\w+):cn_po_origin\(sample, variant\) == \\"(\w+)\\"'/g, function (m, tag, po) {
  wired.po[po] = true;
});
var wiring = [];
function problem_wiring(msg) {
  wiring.push(msg);
}
Object.keys(KNOWN).forEach(function (call) {
  if (!wired.moi[call]) problem_wiring("no --family-expr for the call " + call);
});
Object.keys(wired.moi).forEach(function (call) {
  if (!KNOWN[call]) problem_wiring("--family-expr for an unknown call " + call);
});
Object.keys(KNOWN_PO).forEach(function (po) {
  if (po !== "null" && !wired.po[po]) problem_wiring("no --sample-expr for the parent of origin " + po);
});
Object.keys(wired.po).forEach(function (po) {
  if (!KNOWN_PO[po]) problem_wiring("--sample-expr for an unknown parent of origin " + po);
});
console.log("\nconf/modules.config (SLIVAR_EXPR):");
if (wiring.length === 0)
  console.log(
    "  " +
      Object.keys(wired.moi).length +
      " MoI and " +
      Object.keys(wired.po).length +
      " parent-of-origin expressions, one per value",
  );
else {
  failed = true;
  wiring.forEach(function (p) {
    console.log("  PROBLEM " + p);
  });
}

console.log("\nwrote moi_scenarios.tsv, moi_scenarios.md, trio_lookup_table.tsv");
if (failed) {
  console.log("\nFAILED — see above");
  process.exitCode = 1;
}
