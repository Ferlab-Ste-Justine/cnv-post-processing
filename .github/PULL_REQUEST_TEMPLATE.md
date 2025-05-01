<!--
# Ferlab-Ste-Justine/cnv-post-processing pull request

Many thanks for contributing to Ferlab-Ste-Justine/cnv-post-processing!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the main branch.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/tree/main/.github/CONTRIBUTING.md)
-->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
- [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/tree/main/.github/CONTRIBUTING.md))
- [ ] Make sure your code lints (`nf-core pipelines lint --release`).
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
- [ ] Check for unexpected warnings in debug mode (`nextflow run . -profile debug,test,docker --outdir <OUTDIR>`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] Reference Data Documentation in `docs/reference_data.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
