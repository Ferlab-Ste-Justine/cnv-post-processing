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

_Decided 2026-10-10, BIOINFO-214. Partly replaced on 2026-10-20: see "One X-linked de novo name, and a son's inherited X keeps inherited_maternal" below._

**What:** on non-PAR chrX, the clear-cut calls take the SNV pipeline's X-linked names (Radiant's transmission modes), by the child's sex, in place of the generic ones: a son's `de_novo` is `x_linked_recessive_de_novo`, a daughter's `x_linked_dominant_de_novo`; a son's `inherited_maternal` (his only X is his mother's) and a daughter's `recessive` (both of her X copies hit) are `x_linked_recessive`. PAR records keep the autosomal names, and chrY is unchanged (SNV has no Y-linked tags).

**Not renamed:** a daughter's single inherited copy keeps `inherited_maternal` / `inherited_paternal`: in SNV, `x_linked_dominant` needs the transmitting parent to be affected, and the CNV calls don't look at affected status, so dominant or carrier can't be told from the CNs. `de_novo_suspect` and the duo calls (`possible_*`, `recessive_lowconf`) keep their names too. The rest of the tags still use their own names, not Radiant's: a full rename to Radiant values would be its own change, and several CNV calls have no Radiant equivalent.

**Where:** `X_LINKED_CALLS`, next to `TRIO_TABLE` in `assets/cnv-slivar-functions.js`. The harness now also checks that `conf/modules.config` has a slivar expression for every call and parent of origin: the first version of this change forgot the three new lines, and the module test caught it.

### "denovo", not "de_novo", in every tag

_Decided 2026-10-20, BIOINFO-214, with the team._

**What:** the de novo mode-of-inheritance tags are spelled like the parent-of-origin ones (`po_denovo`, `po_possible_denovo`): `de_novo` → `denovo`, `de_novo_suspect` → `denovo_suspect`, `possible_de_novo` → `possible_denovo`, `x_linked_recessive_de_novo` → `x_linked_recessive_denovo`, `x_linked_dominant_de_novo` → `x_linked_dominant_denovo`. The decisions above keep their names as they were when written.

**Why:** one spelling across tags. The SNV pipeline is to be renamed the same way, separately.

**Consequence:** the two X-linked names are no longer Radiant's `transmission_mode` values verbatim (Radiant spells them `x_linked_recessive_de_novo` / `x_linked_dominant_de_novo`), so reading them into Radiant needs a mapping, or Radiant renames too.

### One computation for the MoI call and the parent of origin

_Decided 2026-10-20, BIOINFO-214, after a colleague's review._

**What:** `TRIO_TABLE` is keyed on the child's deviation and a transmission outcome (`get_transmission_outcome()`): every way one copy from each parent adds up to the child's CN, each passed copy named normal, a hit or opposite. The same row gives the MoI call and the parent of origin. The parent labels (`same`, `half`, `other`), which compared totals, are gone.

**Why:** the labels and the haplotype enumeration could disagree: a parent homozygous for a deletion (CN 0) is neither `same` nor `half`, so the call fell to `inheritance_unresolved` while the parent of origin was known. On the real NA12878 trio that was 139 of 2,532 records, all with a CN 0 parent; they are now `recessive` / `x_linked_recessive`.

**Choices made:**

- A missing parent (not in the PED, or no CN) passes 0, 1 or 2 copies, like any parent with CN 0 to 3. A way that needs more only counts for the Mendelian check (`inheritance_unresolved`), unless `MENDEL_ALLOW_RARE_SPLITS` is on.
- A copy changed in the opposite direction (a deletion passed to a child with a gain) explains nothing: such a child is `undecided`.
- A parent with CN 4 or more still can't be split: its child is `undecided`, and `is_mendelian_possible()` decides between `inheritance_unresolved` and `mendelian_inconsistent` (T16, T17).

### One X-linked de novo name, and a son's inherited X keeps inherited_maternal

_Decided 2026-10-20, BIOINFO-214, after a colleague's review._

**What:** on non-PAR chrX, a de novo CNV is `x_linked_denovo` for a son or a daughter, in place of `x_linked_recessive_denovo` (son) and `x_linked_dominant_denovo` (daughter). A son's CNV inherited from his mother keeps `inherited_maternal` instead of becoming `x_linked_recessive`, so `x_linked_recessive` is only a daughter with both X copies hit. 14 mode-of-inheritance tags instead of 15.

**Why:** the SNV pipeline's names (Radiant's `transmission_mode`) read dominant or recessive off the zygosity and who is affected. Radiant's lookup gives a hemizygous son (1/1) with ref parents `x_linked_recessive_de_novo`, a het daughter `x_linked_dominant_de_novo`, and a son `x_linked_recessive` only when his het mother is unaffected. For CNVs that breaks twice:

- A son's gain (CN 2 on non-PAR X: his normal X plus an extra copy) matches no SNV genotype, and "recessive" says nothing true about it (harness X12, de novo; X13, from his mother). This was the colleague's point.
- The calls don't look at who is affected, so a son who inherits a deletion from an affected mother got `x_linked_recessive` where the SNV side wouldn't. That is why a daughter's single inherited copy already kept `inherited_*` (2026-10-10).

The CN and the sex already say what the names inferred.

**Consequence:** the X-linked tags now leave the SNV names. Since the `denovo` rename, Radiant needs a mapping to read them anyway. On the real NA12878 trio, the daughter's 5 de novo records on non-PAR chrX change from `x_linked_dominant_denovo` to `x_linked_denovo`; nothing else changes.

### Each sample keeps its own call ID in FORMAT/OID

_Decided 2026-10-28, BIOINFO-214, after a colleague's suggestion._

**What:** `NORMALIZE_CNV` copies each per-sample record's `ID` into that sample's `FORMAT/OID` (`bin/add_original_id.awk`, set by `ext.original_id`). The collapses still keep one `ID` per cluster, the leftmost call's, but `OID` names each sample's own call in every family-level VCF. `:` is percent-encoded as `%3A`, since it separates FORMAT values and DRAGEN IDs contain it. Joint-called VCFs don't get `OID`: they skip the collapses, so their `ID` is already the caller's.

**Why:** after `TRUVARI_COLLAPSE_COHORT`, a record's `ID` is one sample's call, often a parent's, and nothing in the record said which call each other sample contributed. The link existed only in truvari's removed VCF (`MatchId` against the kept record's `CollapseId`), which nothing downstream reads. truvari keeps one `ID`, but it copies all FORMAT fields of each call it folds in (into samples that are still ref or missing), so a per-sample `OID` travels with the genotype. (Except for calls whose GT truvari reads as unset, `./1` and haploid `1`: `RESTORE_COLLAPSED_FORMAT` puts those back, see below.) On the real NA12878 trio the published slivar VCF is unchanged apart from `OID`; at chr22:17.29M, for example, the record's `ID` is NA12891's call while NA12878's and NA12892's `OID` both name `17288730-17294633`.

**Rejected:** making the proband's call each cluster's representative (QUAL rewrite + `--keep maxqual`, briefly on this branch). It changed which coordinates a record shows but still lost the other samples' IDs; `OID` keeps all of them without changing the collapse.

**Not covered:** a sample with two calls in one cluster (e.g. a fragmented CNV) keeps only the first one's `OID`; the second survives only in the removed VCF. The collapse membership itself is unchanged: a call is still compared with the cluster's representative only, and DRAGEN's `SVTYPE=CNV` still lets a DEL and a DUP at the same locus collapse together.

### Each proband gets its own calls back, joined to the slivar VCF through OID

_Decided 2026-10-28, BIOINFO-214._

**What:** `JOIN_PROBAND_CNV` (`bin/join_proband_cnv.py`) runs after slivar, once per affected child in the PED (`pedAffectedChildren`). For each slivar record where the proband's `OID` isn't `.`, it writes the proband's own call from its input VCF (ID, coordinates, QUAL, FILTER, INFO), plus the INFO keys only the slivar record has and the family's FORMAT. `COHORT_ID`/`COHORT_POS`/`COHORT_END` point back to the slivar record. Published as `slivar/<familyId>.<proband>.cnv.proband.vcf.gz`.

**Why:** the slivar record shows the call truvari kept, often a parent's, so the proband's own coordinates and caller metrics were hidden. Downstream (Radiant) wanted to join the proband's input VCF to the slivar VCF on `ID` = decoded `OID`; this does it once, in the pipeline. The join key is the proband's own column, since DRAGEN gives identical calls in different samples the same ID (18 of the trio's chr22 IDs are in more than one sample).

**Rejected:** rewriting the slivar record's ID and coordinates to the proband's call. Each member of a cluster was matched against the kept call only, so a parent's call can fail to match the proband's call it would then be shown with; and VEP, depth refinement and slivar would have to rerun on the new interval. Also rejected: matching on start-end, since the record's `POS`/`END` aren't the proband's, and DRAGEN's `POS` is one base before the start in its ID.

**Accepted:** `CSQ` and the slivar tags describe the truvari cluster, computed on the interval of the call kept to represent it, which may be a parent's. That's intended: the proband's call is in the cluster because it matched that call, and the slivar tags come from each sample's copy number, not from the interval. The one gap is `CSQ`: with `--pctovl 0.3`, the proband's interval can extend well beyond the kept call, and a gene overlapped only by that extra part isn't listed. Radiant has to be told.

### Calls truvari drops from the cohort collapse are copied back from its removed VCF

_Decided 2026-10-29, BIOINFO-214, from the first full-genome NA12878 trio run._

**What:** `RESTORE_COLLAPSED_FORMAT` (`bin/restore_collapsed_format.py`) runs between `TRUVARI_COLLAPSE_COHORT` and its sort. For each kept record, every sample with no GT allele set gets the FORMAT of its call among the absorbed ones (truvari's removed VCF, `MatchId` = the kept record's `CollapseId`), the highest `TruScore` first, and `NumConsolidated` counts it.

**Why:** truvari's `collapse_into_entry` only copies FORMAT from an absorbed call whose GT `truvari.get_gt()` reads as set. DRAGEN writes `./1` for a gain of unknown zygosity (CN 3+), which it reads as `GT.NON`, and `1` for a haploid call, which it reads as `GT.UNK`. Such a call absorbed into another sample's was dropped with its CN and OID: the sample looked like it had no call, slivar tagged it `unknown_cn`, and the proband VCF left the call out. On the full-genome NA12878 trio: 36 calls (NA12878 9, NA12891 15 with 2 haploid chrX, NA12892 12), exactly the 726 - 717 gap between NA12878's input calls and its proband VCF; the fix restores all 36, record count unchanged.

**Rejected:** rewriting GT before the collapse (`./1` -> `0/1`): it changes the GT Exomiser reads, and there's no honest diploid rewrite of a haploid `1`. Swapping GT out for truvari and back: fragile.

### Depth refinement writes the expected CN, sex-aware on chrX/chrY

_Decided 2026-10-29, BIOINFO-214, from the first full-genome NA12878 trio run._

**What:** when `bin/refine_genotypes.py` resolves a missing genotype to hom-ref, it also sets `FORMAT/CN` to the sample's expected CN, computed like `get_expected_cn()` in `assets/cnv-slivar-functions.js` (same GRCh38 PARs; sex from the family's PED, now an input of `DEPTH_GENOTYPE_REFINE`). The depth ratio is divided by (expected CN / 2) before the 0.8-1.2 band check, and the GT has as many alleles as the expected CN. No expected CN (chrM, PAR-straddling, unknown sex on non-PAR X/Y, a female's Y): left missing.

**Why:** slivar reads CN only. A `0/0:.` sample was still unknown to it, so on the full-genome trio 656 of 1212 records were `unknown_cn`, 207 `possible_inherited_*`, and `denovo` (both parents known and normal) never fired, against 0 `unknown_cn` and 43 `denovo` on the same trio's DRAGEN joint VCF, which gives every sample a CN. The band was also wrong for a male's non-PAR chrX: his ratio there is half a female's (median 0.34 for NA12891 against 0.55 and 0.59), so only 5% of his chrX regions fell in the band.

**Rejected:** forcing every `.` CN to 2. A `./.` that depth didn't confirm as normal may be a CNV the caller missed in a parent; writing 2 there would turn the child's inherited CNV into a `denovo`. Also rejected: having slivar read `0/0` with no CN as normal, which would make the rules depend on GT.

**Not covered:** whether a child with no call of its own should get any MoI tag at all (it gets `unknown_cn` when its CN is missing; see the open question "A child with no call of its own: `unknown_cn`, or no MoI tag?"); the 0.8-1.2 band itself is still an unvalidated default.

### Depth refinement only on records of 5 kb or more

_Decided 2026-10-29, BIOINFO-214, from the first full-genome NA12878 trio run._

**What:** `DEPTH_GENOTYPE_REFINE` leaves records shorter than 5 kb (END - POS) missing (`ext.min_size`, default 5000; `--min-size` in `bin/refine_genotypes.py`). Samples with no call of their own there stay `./.` with no CN, so slivar gives the child `possible_*` or `unknown_cn` instead of `denovo`.

**Why:** once refinement wrote CN, 28 records got `denovo`, both parents read as CN 2 by depth (ratios 0.81-1.13). Against the trio's DRAGEN joint VCF, only 10 were `denovo` there too; 8 were inherited (a parent CN 1), 8 had all three samples CN 1, 2 had no match. 21 of the 28 were under 5 kb, and 13 of the 16 that disagree with the joint VCF. Two measurements on the same run say small records are where both the caller and depth fail:

- Depth ratio by the sample's own call (autosomes): a real one-copy loss reads 0.54-0.59 at the median, but 5-7% of them fall in the 0.8-1.2 band. Samples with no call read a median 0.76 under 5 kb against 1.04 at 5-50 kb, so many small "no call" regions are probably losses the per-sample caller missed.
- truvari bench of each sample's calls in the final VCF against the joint VCF (`--pctseq 0 --pctovl 0.5 --pctsize 0.5 --refdist 2000`): precision 0.97, recall 0.68-0.75 over all sizes; 84% of the missed joint calls are under 5 kb, 85% of them losses. At 5 kb and more: precision 0.89-0.90, recall 0.87, F1 0.88. The final VCF scores within 4 calls of the raw per-sample DRAGEN VCFs, so this is the per-sample caller against the joint caller, not something the pipeline adds or loses.

**Rejected:** tightening the band (e.g. 0.9-1.1): most of the misread parents were at 0.85-0.98, and it would leave many truly normal parents missing. Marking depth-inferred CN (e.g. `FORMAT/CNSRC`) and downgrading `denovo` built on it to `denovo_suspect`: still possible on top of this, but it changes the slivar rules and needs a geneticist's review.

**Note:** DRAGEN's per-sample caller sets `cnvLength` on every call under 10 kb, while its joint caller lets small calls PASS, so a PASS-only filter downstream (Radiant, Exomiser) never sees our sub-10 kb calls. At 10 kb and more, PASS only, our calls score precision 0.92-0.95, recall 0.76-0.79 against the joint VCF.

### Family-mode Exomiser runs on the proband VCF

_Decided 2026-10-30, BIOINFO-214._

**What:** `EXOMISER_WORKFLOW` now runs after `JOIN_PROBAND_CNV`, on the proband VCF of the child the `familyPheno` phenopacket is about (`proband.subject.id`, read by `phenopacketProbandId()`). Its reports are `<familyId>.<proband>.exomiser.*`. A family with no such VCF (`DRAGEN_JOINT`, no `familyPed`, or a phenopacket proband who isn't a sequenced affected PED child) runs on the family VCF as before, `<familyId>.exomiser.*`, and `--exomiser_start_from_vep` only applies to those families.

**Why:** on the family VCF, Exomiser scored the truvari-kept call, often a parent's, so a gene overlapped only by the proband's own interval could be missed: the same gap as `CSQ` in the previous decision. The proband VCF has the proband's own calls and intervals, and every family member's FORMAT, so the pedigree analysis is unchanged. On the full-genome NA12878 trio, Exomiser read all 726 of NA12878's calls (1237 variant-gene pairs: a CNV counts once per gene it overlaps) with the trio pedigree. Its reports were empty, as on the family VCF: with the test phenotype data, the test analysis file's `priorityScoreFilter` keeps 2 of 19,171 genes. Rerun by hand without that filter: 261 pairs from the 151 `PASS` calls passed the failed-variant filter (`PASS_ONLY`; the other 575 calls are `cnvLength`), 18 passed the inheritance filter, and the 5 CNVs in Exomiser's VCF carry the proband's own call IDs.

**Rejected:** one Exomiser run per affected child: the family phenopacket describes one proband's phenotypes (the Phenopacket schema's `Family.proband` isn't repeated), and another child's calls would be scored against them. The other affected children still count through the pedigree, in Exomiser's segregation filters.

**Not done (yet):** a report for each affected child, from the `pheno` column of its own row (ignored on family rows today): Exomiser on that child's proband VCF with its own phenopacket, the family PED for segregation (`--ped` with `--sample`, to be checked). Worth it for affected sibs with different phenotypes, much less for identical twins. Decided 2026-10-30 to keep one run per family for now.

**Accepted:** this differs from the SNV pipeline, which runs Exomiser on the family VCF, in parallel with slivar. Exomiser's scores stay in its own reports and aren't written into the proband VCF; doing that would be a further step joining on the proband's call IDs.

### A parent who passes nothing on chrX/chrY passes nothing even when absent or without CN

_Noted 2026-10-05 as an open question; decided 2026-10-30, BIOINFO-214, from review v2 (R2, 2026-10-23)._

**What:** in `get_parent_transmission()` (`assets/cnv-slivar-functions.js`), on non-PAR chrX a father passes nothing to a son and, on non-PAR chrY, a mother passes nothing to anyone, whether or not that parent is in the PED or has a CN. `does_parent_pass_nothing()`, `get_transmission_outcome()` and `is_mendelian_possible()` follow from it: such a parent reads as known and normal, and their CN is never needed. Harness scenarios X14-X16 lock it in: a mother-son duo on chrX is judged on the mother alone (son CN 0, mother CN 1, father CN missing: `inherited_maternal` / `po_mother`; son CN 0, mother CN 2, no father: `x_linked_denovo` / `po_denovo`; son CN 2, mother CN 0, no father: `mendelian_inconsistent`).

**Why:** it's biology, not missing data. A son's X can only be his mother's and his Y his father's, so a father's X and a mother's Y never reach him, whatever the PED says. A mother-son duo is a complete trio for non-PAR chrX.

**Rejected:** treating such a parent as "could have passed anything" when absent or without CN. A mother-son duo on chrX would get `possible_*` calls and no parent of origin, as if the X could have come from the father, and since DRAGEN gives mothers no CN on chrY, every son's chrY call would become unresolved.

**Caveats:** this, like the whole expected-CN model, trusts the PED's sex column: a wrong sex makes every non-PAR chrX/chrY call of that family wrong (see the TODO "Check the PED sex against the data"). It also assumes 46,XX or 46,XY children; XXY, XYY or X0 break the expected-CN model everywhere, not just this rule.

### SVTYPE is set from the ALT, so truvari keeps losses and gains apart

_Decided 2026-10-30, BIOINFO-214._

**What:** `NORMALIZE_CNV` (both paths) sets `INFO/SVTYPE` to `DEL` or `DUP` from the split `<DEL>` / `<DUP>` ALT (`bin/svtype_from_alt.awk`), replacing DRAGEN's `SVTYPE=CNV`. Other ALTs are left as they are.

**Why:** truvari 5.4.0 takes a record's type from `INFO/SVTYPE` before the ALT (`VariantRecord.var_type()`), and `CNV` isn't one of its types, so every DRAGEN call was `UNK` and `same_type()` always passed: in effect `--typeignore`. DRAGEN's `GAIN` / `LOSS` are only in the ID, which truvari never reads. On the full-genome trio, 4 of the 168 calls the cohort collapse absorbed had joined a call of the other type; the absorbed sample's FORMAT was kept, so its GT then pointed at the wrong ALT (chr14:73530575: the father's CN 3 gain as `./1` on NA12878's `<DEL>`). slivar was right (it reads CN), but anything reading GT (Exomiser's segregation, Radiant) saw him as a deletion carrier, and VEP described the record as the loss only. Tested both ways in `modules/local/truvari_collapse/tests` (toy_loss_gain).

**Rejected:** putting `CNV` back after the collapse, for VEP and Exomiser: they gave the same output on the trio's joint VCF with `DEL` / `DUP` (same VEP consequences, Exomiser read every record), and `DEL` / `DUP` is the usual SVTYPE for CNVs.

**Accepted:** a sample whose call of the other type now has its own record is `./.` on the first one, and depth refinement doesn't fill it in (a gain's ratio is outside the hom-ref band). At chr14:73530575 NA12878 goes from `inherited_maternal` to `possible_inherited_maternal`; it's the only tag change on the trio. The proband VCF keeps the proband's own INFO from its input VCF, so its `SVTYPE` is still `CNV`.

## Open questions

### A child with no call of its own: `unknown_cn`, or no MoI tag?

_Noted 2026-10-29, BIOINFO-214, from the full-genome NA12878 trio run. To decide with the geneticists and Radiant._

**Context:** `cnv_moi()` checks for a missing CN before it checks whether the child carries anything:

```js
if (is_cn_missing(s, variant)) return "unknown_cn"; // a CN we need is missing (P15)
if (get_cn_dev(s, variant) === 0) return null; // not a carrier
```

So a child with no CN gets `unknown_cn`, whether or not it has a call. A family record exists because some sample has a CNV there; when that's only a parent, the child has no call of its own (`./.`) and gets a CN only if depth refinement confirms it's normal, which it never does under 5 kb (see "Depth refinement only on records of 5 kb or more").

**On the full-genome trio** (after the 5 kb floor), 538 of 1212 records are `unknown_cn`, all for NA12878:

| NA12878 at the record             | Parents' CN                  | Records                                  |
| --------------------------------- | ---------------------------- | ---------------------------------------- |
| no call of its own (`./.`, no CN) | at least one parent has a CN | 383 (321 under 5 kb, 63 at 5 kb or more) |
| its own call                      | both parents missing         | 155 (141 under 5 kb, 14 at 5 kb or more) |

The 383 are mostly a parent's CNV the child didn't inherit. `unknown_cn` is true of them (the child's CN is unknown) but it's a third of the family VCF's records, tagged "unknown" where the proband has nothing. The proband VCF (`JOIN_PROBAND_CNV`) only keeps records where the proband has a call, so there it's already the 155 meaningful ones only.

**Options:**

1. Keep as is: honest, and the proband VCF already hides the 383.
2. No MoI tag for a child with no CN (`return null` instead of `unknown_cn` when the child's own CN is missing); `unknown_cn` then only means "the child carries a CNV but the parents' CN we need is missing" (the 155). A child carrying a CNV the caller missed loses nothing, since there's no call to classify. Changes the meaning of a tag: one line in `cnv_moi()`, scenario P15 in `assets/cnv-slivar-moi-harness.js`, and the golden file re-accepted.
3. A separate tag (e.g. `child_no_call`) for "a parent's CNV, the child's status unknown"; the parents' own genotypes already say as much.

Leaning towards 2, pending that discussion.

## TODO

### Check the PED sex against the data

_Noted 2026-10-30, BIOINFO-214, from review v2 (R2)._

The slivar rules and depth refinement take each sample's sex from the PED, and nothing checks it. A wrong sex silently makes every non-PAR chrX/chrY call of the family wrong. Compare the PED sex with the data before slivar (DRAGEN's own sex/ploidy estimate, or chrX depth or heterozygosity), and stop on a mismatch, like `validatePedAgainstVcf`.

### Have a geneticist review the slivar MoI scenarios

_Noted 2026-10-06, BIOINFO-214._

**Why:** `assets/moi_scenarios.expected.tsv`, the golden file of the slivar rules harness, records what the code does, not what's right. It was produced with `node assets/cnv-slivar-moi-harness.js --accept`, and the harness and `ci-slivar-harness.yml` only check that the code keeps matching it. It only becomes a spec once a geneticist has gone through it and agreed. That review is part of the development work too.

**What to hand over:**

- `assets/moi_scenarios.md`: the 62 named scenarios (trios, duos, a singleton, families with two children, chrX/chrY/chrM), each with its family (region, each member's CN and affected status), the MoI call per child and the parent of origin.
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
