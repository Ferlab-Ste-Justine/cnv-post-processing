# Journal

Running notes on decisions, open questions and follow-ups for this pipeline.

## Open questions

### A parent who passes nothing on chrX/chrY: fixed at zero, or "any value" when absent?

_Noted 2026-10-05, BIOINFO-214. To discuss with the team._

**Context:** in `assets/cnv-slivar-functions.js`, `get_parent_transmission()` says what a parent can pass to a child on a given record. Two things build on it: `get_transmission_origin()`, behind the `TRIO_TABLE` origin column and the parent-of-origin tags (`po_mother`, `po_father`, ...), and `transmissions()`, behind the `mendelian_inconsistent` check. On non-PAR chrX/chrY some parents pass nothing at all, biologically: a father passes no X to a son, and a mother passes no Y.

**Current behaviour (since the transmission-engine change):** such a parent counts as passing nothing (`[0]`) even when they are absent from the PED or have no CN. Before, an absent or CN-less parent always counted as "could have passed any value" (`0..kid CN`), whoever they were.

**What it changes:**

- Mendelian check: a son on non-PAR X whose father isn't in the PED is now judged on his mother's copies alone. Example: son CN 2, mother CN 0, no father → was `inheritance_unresolved`, now `mendelian_inconsistent`. (A father in the PED but without a CN hits `unknown_cn` first, so that case doesn't change.)
- Parent of origin: a son's X origin resolves from his mother alone (a mother–son duo gets `po_mother`), and a son's Y origin from his father alone.

**If we switch to "any value" instead:**

- The example above goes back to `inheritance_unresolved`.
- A mother–son duo on chrX no longer gets a parent of origin.
- Parent of origin on chrY never resolves for anyone: DRAGEN reports no CN for females on chrY, so the mother would always count as unknown. If we want "any value", it should probably apply to an absent parent only, not to a parent with no CN on a chromosome they don't carry.

**Where to change it:** the first `if (hemi)` block in `get_parent_transmission()`; `transmissions()` and `get_transmission_origin()` follow from it.

## TODO

### Have a geneticist review the slivar MoI scenarios

_Noted 2026-10-06, BIOINFO-214._

**Why:** `assets/moi_scenarios.expected.tsv`, the golden file of the slivar rules harness, records what the code does, not what's right. It was produced with `node assets/cnv-slivar-moi-harness.js --accept`, and the harness and `ci-slivar-harness.yml` only check that the code keeps matching it. It only becomes a spec once a geneticist has gone through it and agreed. That review is part of the development work too.

**What to hand over:**

- `assets/moi_scenarios.md`: the 47 named scenarios (trios, duos, a singleton, families with two children, chrX/chrY/chrM), each with its family tag, child calls and parent of origin.
- `assets/trio_lookup_table.tsv`: every autosomal kid/mom/dad CN combination with its child call and parent of origin, all 144 rows. For example, its first row says kid CN 0, mom CN 0, dad CN 0 → `inheritance_unresolved`, origin `both`. A geneticist can check calls like that at a glance. These are child calls (`TRIO_TABLE`), before `FAMILY_TABLE` uses who is affected to turn them into the family tag.
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
