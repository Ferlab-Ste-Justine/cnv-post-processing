# `Ferlab-Ste-Justine/cnv-post-processing`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving Ferlab-Ste-Justine/cnv-post-processing.

We manage the required tasks for `Ferlab-Ste-Justine/cnv-post-processing` using JIRA. Please ensure that your changes are linked to a JIRA ticket for proper context and tracking.

Every branch or commit must include the corresponding JIRA ticket number to ensure traceability and alignment with Ferlab development's best practices ([see Ferlab Developer Handbook](https://www.notion.so/ferlab/Developer-Handbook-ca9d689d8aca4412a78eafa2dfa0f8a8))

> [!NOTE]
> If you need help using or modifying Ferlab-Ste-Justine/cnv-post-processing then the best place to ask is on the #bioinfo Slack channel.

## Contribution workflow

If you'd like to write some code for Ferlab-Ste-Justine/cnv-post-processing, the standard workflow is as follows:

1. **Create a feature branch**

- Follow the branch naming conventions outlined in the [Ferlab Developer Handbook](https://www.notion.so/ferlab/Developer-Handbook-ca9d689d8aca4412a78eafa2dfa0f8a8).

2. **Make necessary changes / additions**

- Implement your changes or additions on the feature branch.
- Use the `nf-core pipelines schema build` command to add any new parameters to the pipeline JSON schema.
- Update the documentation to reflect your changes, including any new parameters, processes, outputs or reference data.
- Run the `nf-core pipelines lint` command to ensure your changes do not introduce any new warnings or failures.
- Test your changes thoroughly. Execute the pipeline locally using the docker executor to ensure it completes successfully. If applicable, incorporate `nf-test` unit tests to validate your changes.

3. **Submit a pull request**

- Submit a pull request against the `main` branch.
- Ensure your commit messages follow the conventions outlined in the [Ferlab Developer Handbook](https://www.notion.so/ferlab/Developer-Handbook-ca9d689d8aca4412a78eafa2dfa0f8a8).
- If your branch contains many commits, we strongly recommend squashing them into a single commit before submitting the pull request.
- Include a detailed description in your pull request. This should summarize the changes made, reference the JIRA ticket, and describe the manual tests you performed to validate your changes. Providing this context will help reviewers understand and evaluate your work more efficiently.

4. **Update the changelog**

- Unless your changes are purely refactoring with no functional impact, update the `CHANGELOG.md` file to document your modifications.
- Clearly describe the change, addition, or fix you made, ensuring it is concise and easy to understand.
- Include the pull request link next to the description to provide additional context and traceability.
- Follow the formatting conventions already used in the `CHANGELOG` to maintain consistency.

5. **Code review and approval**

- At least one approval from a codeowner is required to merge the pull request.
- Address any feedback provided during the review process.

Once your pull request is approved and all checks pass, you will be authorized to merge your feature branch into the `main` branch.

If you're not used to this kind of workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

You have the option to test your changes locally by running the pipeline. For receiving warnings about process selectors and other `debug` information, it is recommended to use the debug profile. Ex:

```bash
nf-test test --profile debug,test,docker --verbose
nextflow run main.nf --profile debug,test,docker --verbose
```

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

Here are the tests that are run:

### Commit linter check

Ensure that the commit messages follow the ferlab naming convention as outlined in the [Developer Handbook](https://www.notion.so/ferlab/Developer-Handbook-ca9d689d8aca4412a78eafa2dfa0f8a8).

### Lint tests

Ensure that the nf-core lint command succeeds and report no failed test.

### Pipeline test

Verify that the pipeline can be executed successfully on a minimalist test dataset. This ensures that the core functionality of the pipeline remains intact.

## Pipeline contribution conventions

To make the `Ferlab-Ste-Justine/cnv-post-processing` code and processing logic more understandable for new contributors and to ensure quality, we semi-standardize the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test command in `.github/workflow/ci.yml`.
9. If applicable, add a description of the output files to `docs/output.md`.
10. Update the pipeline description in the `README.md`. You may need to update images as well.

### Default values

Parameters should be initialized / defined with default values within the `params` scope in `nextflow.config`.

Once there, use `nf-core pipelines schema build` to add to `nextflow_schema.json`.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. A nf-core standard set of labels that should be followed where possible can be seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config), which has the default process as a single core-process, and then different levels of multi-core configurations for increasingly large memory requirements defined with standardized labels.

The process resources can be passed on to the tool dynamically within the process with the `${task.cpus}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Nextflow version bumping

You may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

You may need to manually update the nextflow versions specified in the github workflow file `ci.yml`.

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).
