// node assets/cnv-slivar-moi-harness.js            run all checks
// node assets/cnv-slivar-moi-harness.js --accept   ... and save the results as the new reference
// (from any directory: every file it reads or writes lives next to this script)
//
// Test bench for cnv-slivar-functions.js (the SLIVAR_EXPR rules) only.
// Builds made-up families shaped like slivar's sample objects and writes:
//   moi_scenarios.tsv       one row per named scenario: tag, child calls, parent of origin
//   moi_scenarios.md        the same table as Markdown
//   trio_lookup_table.tsv   TRIO_TABLE expanded to concrete autosomal CNs (kid, mom, dad -> call)
//
// Reference ("golden file"): moi_scenarios.expected.tsv
//   Every run compares moi_scenarios.tsv with it and lists the scenarios whose
//   results changed. Each change is a behaviour change you should be able to
//   explain. When you agree with the new results, re-run with --accept.
//
// Random families: 200000 random pedigrees are run to check that nothing
// crashes and that cnv_moi() only ever returns a known tag (or null).
//
// Exit code 1 if anything changed versus the reference or a check failed.
var T = require("./cnv-slivar-functions.js");
var fs = require("fs");
var path = require("path");
function here(f) {
  return path.join(__dirname, f);
} // all files live next to this script, whatever the cwd
var MISSING = -127; // slivar's sentinel for a missing FORMAT int
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
// A record with no INFO/END (e.g. the slivar module-test fixtures): slivar sets stop
// to POS, and reading INFO.END throws.
var AUTO_NO_END = strict(
  { CHROM: "chr7", POS: 1000001, start: 1000000, stop: 1000001, INFO: strict({ SVTYPE: "DEL" }, "INFO") },
  "variant",
);

function child_calls(fam, v) {
  return (
    fam
      .filter(function (s) {
        var d = T.cn_dev(s, v);
        return s.affected && ("mom" in s || "dad" in s) && d !== null && d !== 0;
      })
      .map(function (s) {
        return s.id + "=" + T.trio_call(s, v);
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
// Prefixes: T trio, X chrX, D duo, S singleton, Q two children
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
  ["T17 kid CN=3, mom CN=4, dad CN=1 (possible, maternal)", AUTO, { kid: [3, 1], mom: [4, 0], dad: [1, 0] }],
  // T18-T22: the scenarios of the retired assets/tests/cnv_slivar_functions.test.js not covered above
  ["T18 dup inherited from AFFECTED mom", AUTO, { kid: [3, 1], mom: [3, 1], dad: [2, 0] }],
  ["T19 kid CN=0, mom CN=2, dad CN=1 (impossible)", AUTO, { kid: [0, 1], mom: [2, 0], dad: [1, 0] }],
  ["T20 affected kid CN missing", AUTO, { kid: [null, 1], mom: [2, 0], dad: [2, 0] }],
  ["T21 unaffected kid CN=1, dad CN=1 (paternal origin)", AUTO, { kid: [1, 0], mom: [2, 0], dad: [1, 0] }],
  ["T22 all unaffected, all CN=3 (origin ambiguous)", AUTO, { kid: [3, 0], mom: [3, 0], dad: [3, 0] }],
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
  ["Q3  quad: dominant + unaffected sib CN=1", AUTO, { kid: [1, 1], mom: [2, 0], dad: [1, 1], sibs: [[1, 0]] }],
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
var header = ["scenario", "tag", "child calls (TRIO_TABLE)", "parent of origin (sample-expr)"];
var rows = scenarios.map(function (sc) {
  var fam = family(sc[2]),
    v = sc[1];
  return [sc[0], T.cnv_moi(fam, v) || "-", child_calls(fam, v), po(fam, v)];
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
  "Scenario prefixes: **T** trio, **X** chrX, **D** duo, **S** singleton, **Q** two children. `-` = no tag.",
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
  c[1] = "`" + c[1] + "`";
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
      changes.push("  + " + id + " (new scenario): " + now[id][1]);
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
var KNOWN = { "-": true, unknown_cn: true, inheritance_unresolved: true, mendelian_inconsistent: true };
T.FAMILY_TABLE.forEach(function (row) {
  KNOWN[row[0]] = true;
});
var KNOWN_PO = { null: true, maternal: true, paternal: true, both: true, ambiguous: true };
var YNP = variant("chrY", 10660601, 10693753);

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
  var v = pick([AUTO, AUTO, XNP, XNP, YNP, AUTO_NO_END]);
  var spec = {
    kid: [pick(CNPOOL), 1],
    kidsex: pick(["male", "female"]),
    mom: rnd(5) ? [pick(CNPOOL), rnd(4) === 0 ? 1 : 0] : null,
    dad: rnd(5) ? [pick(CNPOOL), rnd(4) === 0 ? 1 : 0] : null,
    sibs: rnd(2) ? [[pick(CNPOOL), rnd(2), pick(["male", "female"])]] : [],
  };
  var fam = family(spec),
    tag;
  try {
    tag = T.cnv_moi(fam, v) || "-";
    fam.forEach(function (s) {
      var po = T.cn_po_origin(s, v);
      if (!KNOWN_PO[po]) problem("unknown parent of origin: " + po, spec, v);
      var d = T.cn_dev(s, v);
      if (!(s.affected && ("mom" in s || "dad" in s) && d !== null && d !== 0)) return;
      var call = T.trio_call(s, v); // every carrier child must get a real call ...
      if (typeof call !== "string") {
        problem("trio_call() returned no call", spec, v);
        return;
      }
      if (call === "inherited_maternal" && po !== "maternal")
        problem("inherited_maternal but parent of origin " + po, spec, v);
      if (call === "inherited_paternal" && po !== "paternal")
        problem("inherited_paternal but parent of origin " + po, spec, v);
      po_counts[v.CHROM + " " + po] = (po_counts[v.CHROM + " " + po] || 0) + 1; // ... and agree with its parent of origin
    });
  } catch (e) {
    problem("crash: " + e.message, spec, v);
    continue;
  }
  if (!KNOWN[tag]) problem("unknown tag: " + tag, spec, v);
  counts[tag] = (counts[tag] || 0) + 1;
}

console.log("\nRandom families (" + N + "):");
var probs = Object.keys(problems);
if (probs.length === 0)
  console.log("  no crashes, no unknown tags, every carrier child got a call that agrees with its parent of origin");
else {
  failed = true;
  probs.forEach(function (p) {
    console.log("  PROBLEM " + p + "\n    e.g. " + problems[p]);
  });
}
console.log("  tags seen:");
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

console.log("\nwrote moi_scenarios.tsv, moi_scenarios.md, trio_lookup_table.tsv");
if (failed) {
  console.log("\nFAILED — see above");
  process.exitCode = 1;
}
