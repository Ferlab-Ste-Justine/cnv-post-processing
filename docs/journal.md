# Journal

Running notes on decisions, open questions and follow-ups for this pipeline.

## TODO

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
