# Journal

Running notes on decisions, open questions and follow-ups for this pipeline.

## Decisions

### MoI tags are per child, with no family-level tag

_Decided 2026-10-09, BIOINFO-214, with the team._

**What:** `cnv_moi()` in `assets/cnv-slivar-functions.js` returns the `TRIO_TABLE` call of each affected carrier child, from the child and its parents only. Each MoI tag is a `--family-expr`, as in the SNV pipeline: `has_moi_call()` is true when at least one child of the family gets that call, and slivar values the tag with the family's affected members (2026-10-09; the same day's first version used one `--sample-expr` per tag, valued with the children themselves, and gave the same tags on the real trio). The parent-of-origin tags stay `--sample-expr`, so they name the child. The family-level step that combined the children's calls into one family tag (`FAMILY_TABLE`) is gone.

**Why:** simpler rules, which a geneticist can review from the trio lookup table alone. Interpreting the family (penetrance, segregation, who is affected) is left to whoever reads the tags, who has the PED and the phenotypes anyway.

**What it drops:**

- `dominant`: a child who inherits the CNV from an affected parent is now `inherited_maternal` / `inherited_paternal` like any other; the parent's status is in the PED.
- `parental_mosaicism_candidate`: two affected sibs with the same de novo CNV now each get `de_novo`.
- The segregation check across the family: an unaffected carrier sib, or an affected sib without the CNV, no longer changes a child's call, and an affected parent carrying a CNV their affected child doesn't carry gets no tag (harness T15). The open question of whether a non-penetrant carrier sib should block `dominant` goes away with it.
- The family-wide Mendelian check: only each carrier child's own CN is checked against their parents.

**What it adds:** families with several affected children can get a different call per child (harness Q9: the family carries both tags, and the `po_*` tags say which child is which), and a son's chrY gets a call: a mother passes no Y, so her missing CN there no longer blocks him (Y2 to Y5). On the real NA12878 trio, all 1,259 records get the same tags as before.

### Parent of origin is a transmission label for every child

_Decided 2026-10-10, BIOINFO-214._

**What:** besides `po_mother`, `po_father`, `po_both` and `po_ambiguous`, `cn_po_origin()` now returns `denovo` when no parent's copy explains the child's CN and both parents are known and normal (a parent who passes nothing on chrX/chrY counts as normal), and `possible_denovo` when one parent is normal and the other unknown (not in the PED, or no CN). So, as in the SNV pipeline, every child carrying a CNV gets a per-child transmission label when one can be given, affected or not.

**Why:** since the mode-of-inheritance tags became `--family-expr`, they list the family's affected members; the `po_*` tags are the only per-child record. Without `po_denovo`, a de novo CNV said nothing about which child it was in, and nothing at all for an unaffected child.

**Difference from SNV:** the SNV pipeline's `po_denovo` also covers any allele neither parent explains. Here a parent carrying a CNV that still can't explain the child's (harness T14, T19) gets no `po_*` tag: that's `mendelian_inconsistent`, not de novo.

**Consequence:** de novo is reported twice for an affected child (`de_novo` on the family, `po_denovo` on the child). The harness checks that the two always agree (`de_novo` / `de_novo_suspect` ⇔ `po_denovo`, `possible_de_novo` ⇒ `po_possible_denovo`). Not added: SNV's `po_possible_mother` / `po_possible_father` (a duo whose present parent carries the CNV, harness D1) and `po_unknown`.

### X-linked names on non-PAR chrX, as in the SNV pipeline

_Decided 2026-10-10, BIOINFO-214._

**What:** on non-PAR chrX, the clear-cut calls take the SNV pipeline's X-linked names (Radiant's transmission modes), by the child's sex, in place of the generic ones: a son's `de_novo` is `x_linked_recessive_de_novo`, a daughter's `x_linked_dominant_de_novo`; a son's `inherited_maternal` (his only X is his mother's) and a daughter's `recessive` (both of her X copies hit) are `x_linked_recessive`. PAR records keep the autosomal names, and chrY is unchanged (SNV has no Y-linked tags).

**Not renamed:** a daughter's single inherited copy keeps `inherited_maternal` / `inherited_paternal`: in SNV, `x_linked_dominant` needs the transmitting parent to be affected, and the CNV calls don't look at affected status, so dominant or carrier can't be told from the CNs. `de_novo_suspect` and the duo calls (`possible_*`, `recessive_lowconf`) keep their names too. The rest of the tags still use their own names, not Radiant's: a full rename to Radiant values would be its own change, and several CNV calls have no Radiant equivalent.

**Where:** `X_LINKED_CALLS`, next to `TRIO_TABLE` in `assets/cnv-slivar-functions.js`. The harness now also checks that `conf/modules.config` has a slivar expression for every call and parent of origin: the first version of this change forgot the three new lines, and the module test caught it.

### "denovo", not "de_novo", in every tag

_Decided 2026-10-20, BIOINFO-214, with the team._

**What:** the de novo mode-of-inheritance tags are spelled like the parent-of-origin ones (`po_denovo`, `po_possible_denovo`): `de_novo` → `denovo`, `de_novo_suspect` → `denovo_suspect`, `possible_de_novo` → `possible_denovo`, `x_linked_recessive_de_novo` → `x_linked_recessive_denovo`, `x_linked_dominant_de_novo` → `x_linked_dominant_denovo`. The decisions above keep their names as they were when written.

**Why:** one spelling across tags. The SNV pipeline is to be renamed the same way, separately.

**Consequence:** the two X-linked names are no longer Radiant's `transmission_mode` values verbatim (Radiant spells them `x_linked_recessive_de_novo` / `x_linked_dominant_de_novo`), so reading them into Radiant needs a mapping, or Radiant renames too.

## Open questions

### A parent who passes nothing on chrX/chrY: fixed at zero, or "any value" when absent?

_Noted 2026-10-05, BIOINFO-214. To discuss with the team._

**Context:** in `assets/cnv-slivar-functions.js`, `get_parent_transmission()` says what a parent can pass to a child on a given record. Three things build on it: `get_transmission_origin()`, behind the `TRIO_TABLE` origin column and the parent-of-origin tags (`po_mother`, `po_father`, ...), `transmissions()`, behind the `mendelian_inconsistent` check, and `does_parent_pass_nothing()`, behind how `TRIO_TABLE` reads that parent and whether their CN is needed. On non-PAR chrX/chrY some parents pass nothing at all, biologically: a father passes no X to a son, and a mother passes no Y.

**Current behaviour (since the transmission-engine change):** such a parent counts as passing nothing (`[0]`) even when they are absent from the PED or have no CN. Before, an absent or CN-less parent always counted as "could have passed any value" (`0..kid CN`), whoever they were.

**What it changes:**

- Mendelian check: a son on non-PAR X whose father isn't in the PED is now judged on his mother's copies alone. Example: son CN 2, mother CN 0, no father → was `inheritance_unresolved`, now `mendelian_inconsistent`. Since MoI tags became per child (2026-10-09), the same goes for a father in the PED without a CN: he used to give `unknown_cn`, and his CN is no longer needed.
- MoI calls: `TRIO_TABLE` reads such a parent as "normal", and their CN isn't needed. So a mother–son duo on chrX, or a father–son duo on chrY, can get `de_novo` rather than `possible_de_novo`, and a son on chrY gets a call although DRAGEN gives his mother no CN there.
- Parent of origin: a son's X origin resolves from his mother alone (a mother–son duo gets `po_mother`), and a son's Y origin from his father alone.

**If we switch to "any value" instead:**

- The example above goes back to `inheritance_unresolved`.
- A mother–son duo on chrX no longer gets a parent of origin.
- Parent of origin on chrY never resolves for anyone, and a son's chrY MoI call would need a rule for his mother: DRAGEN reports no CN for females on chrY, so she would always count as unknown. If we want "any value", it should probably apply to an absent parent only, not to a parent with no CN on a chromosome they don't carry.

**Where to change it:** the first `if (hemi)` block in `get_parent_transmission()`; `transmissions()`, `get_transmission_origin()` and `does_parent_pass_nothing()` follow from it.

## TODO

### Have a geneticist review the slivar MoI scenarios

_Noted 2026-10-06, BIOINFO-214._

**Why:** `assets/moi_scenarios.expected.tsv`, the golden file of the slivar rules harness, records what the code does, not what's right. It was produced with `node assets/cnv-slivar-moi-harness.js --accept`, and the harness and `ci-slivar-harness.yml` only check that the code keeps matching it. It only becomes a spec once a geneticist has gone through it and agreed. That review is part of the development work too.

**What to hand over:**

- `assets/moi_scenarios.md`: the 51 named scenarios (trios, duos, a singleton, families with two children, chrX/chrY/chrM), each with its family (region, each member's CN and affected status), the MoI call per child and the parent of origin.
- `assets/trio_lookup_table.tsv`: every autosomal kid/mom/dad CN combination with its child call and parent of origin, all 144 rows. For example, its first row says kid CN 0, mom CN 0, dad CN 0 → `inheritance_unresolved`, origin `both`. A geneticist can check calls like that at a glance. These calls are the MoI tags the pipeline writes for an affected child (see Decisions: there is no family-level tag).
- The open questions above.

**After the review:**

- A result they disagree with is either a rule change in `assets/cnv-slivar-functions.js` or a new open question.
- A case they bring up that the table doesn't cover becomes a new named scenario in the harness, which also makes it a test.

### Replace `commit_lint.yml` with the PR title check

_Noted 2026-09-28, BIOINFO-214._

**Current state:** both checks run on this repo.

- `.github/workflows/ci-pr-title-lint.yml` (job `PR Title Lint Check`) checks that the PR title follows `<type>: <TICKET-123> <description>`. It's copied unchanged from Post-processing-Pipeline, the fleet reference, where it replaces commit linting. The reasoning is in that repo's `docs/journal.md` ("Check the PR title, not every commit message"). In short, with squash merging the PR title becomes the commit message on `main`, so it's the one message worth checking.
- `.github/workflows/commit_lint.yml` (job `Commit Lint Check`) still runs `Ferlab-Ste-Justine/action-commit-lint` on every push. It's left in place on purpose, pending the checks below.

**Why it's pending:**

- `main`'s branch protection requires the `Commit Lint Check` status (checked through GitHub's public API on 2026-09-28; enforced for non-admins). Removing `commit_lint.yml` before that rule changes would block every PR, waiting for a check that never reports.
- Every PR merged so far (#1 to #12, up to 2025-08-01) was merged with a merge commit, not squashed. `action-commit-lint` relies on that: it checks every commit back to the last `Merge pull request #` commit. The repo's current merge settings (allowed methods, default squash commit message) aren't visible without admin rights, so whether squash merging is already enabled hasn't been confirmed.

**To do, with a repository admin:**

1. Confirm which merge methods are enabled, and switch to squash merging as the lab guideline says.
2. Set Settings → General → Pull Requests → "Default commit message" for squash merging to "Pull request title" (or "Pull request title and commit details"), so the checked title is what lands on `main`.
3. In the `main` branch protection rule, replace the required check `Commit Lint Check` with `PR Title Lint Check`.
4. Then delete `.github/workflows/commit_lint.yml`, and update `CLAUDE.md` and this entry.
