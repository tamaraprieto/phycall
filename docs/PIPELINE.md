# Pipeline — step-by-step run order

Every `*.PTA.sh` script is a SLURM job that takes a **config file** as its first
argument (some take an extra argument such as a chromosome, window size or depth):

```bash
sbatch [--array=1-N] src/<Script>.PTA.sh /path/to/Config.<PATIENT>.txt [extra-arg]
```

The config file is read by `src/ReadConfig.sh`, which exports the variables the scripts
rely on (`ORIDIR`, `WORKDIR`, `RESDIR`, `REF`, `IDLIST`, `SAMPLELIST`,
`PLATFORM`, `LIBRARY`, `SCRIPTDIR`, ...). Run `src/ReadConfig.sh --help` for the full
list of config keys.

Each script sources it with a bare `source ReadConfig.sh "$1"`, so `src/` must be on
your `$PATH` before launching jobs:

```bash
export PATH="$PWD/src:$PATH"
```

(Alternatively, edit the scripts to source it by absolute path,
`source "$SCRIPTDIR/ReadConfig.sh" "$1"`.)

---

## A. Read processing

| # | Step | Script | Example |
|---|------|--------|---------|
| 1 | Adapter trimming (Cutadapt) | `src/CutAdapt.PTA.sh` | `sbatch --array=1-356%120 src/CutAdapt.PTA.sh Config.4295.txt` |
| 2 | Map reads (BWA-MEM) | `src/BWA.PTA.sh` | `sbatch --array=1-356%50 src/BWA.PTA.sh Config.4295.txt` |
| 3 | Sort (Picard SortSam) | `src/SortSam.PTA.sh` | `sbatch --array=1-356%115 src/SortSam.PTA.sh Config.4295.txt` |
| 4 | Merge read groups + mark duplicates | `src/MergeWithMarkDuplicates.PTA.sh` | `sbatch --array=1-89 src/MergeWithMarkDuplicates.PTA.sh Config.4295.txt` |
| 5 | Cross-sample contamination | `src/CalculateContaminationSingleCell.PTA.sh` | `sbatch --array=1-89 src/CalculateContaminationSingleCell.PTA.sh Config.4295.txt` |
| 6 | Base quality score recalibration (BQSR I + II) | `src/BaseRecalibratorI.PTA.sh`, `src/BaseRecalibratorII.PTA.sh` | `sbatch --array=1-89%45 src/BaseRecalibratorI.PTA.sh Config.4295.txt` |

## B. Germline variant calling

| # | Step | Script |
|---|------|--------|
| 7 | HaplotypeCaller (per-sample GVCF) | `src/HaplotypeCaller.PTA.sh` |
| 8 | Joint genotyping (per interval — faster) | `src/GenotypeGVCFs.intervals.PTA.sh` (slow alt: `src/GenotypeGVCFs.PTA.sh`) |
| 9 | Concatenate interval VCFs | `src/ConcatenateVCFs.PTA.sh` |
| 10 | VQSR filtering | `src/VQSR.PTA.sh` |
| 11 | Somatic recall at invariable sites | `src/RecallSomatic.PTA.sh` |
| 12 | Allelic counts (recall / het sites) | `src/CollectAllelicCounts.PTA.sh` |

## C. Allelic imbalance / ADO

| # | Step | Script |
|---|------|--------|
| 13 | Allelic counts at het sites (chr21) | `src/CollectAllelicCounts.PTA.sh` |
| 14 | ADO vs depth — PTA H1 cell, in silico downsample | `src/DowsampleH1a.PTA.sh`, `src/DowsampleH1b.PTA.sh` |
| 15 | Amplification error at invariable sites | `src/ReadErrorsAtInvariableSites.sh` |

## D. CNV / LOH detection (allelic-imbalance method)

| # | Step | Script |
|---|------|--------|
| 16 | Phase germline het SNPs (Shapeit v4) | `src/ShapeitAI.PTA.sh` → `src/MergeShapeit.PTA.sh` |
| 17 | Barcode BAMs for Chisel | `src/PrepareChisel.PTA.sh` |
| 18 | Chisel A/B allele counts (5 kb / 50 kb) | `src/Chisel.PTA.sh` |
| 19 | Build per-window BAF count tables | `src/BAFprep.PTA.sh` |
| 20 | BAF profiles + haplotype-switch minimization | `src/BAF.PTA.sh` → `src/BAF.LOH.reduceswitches.py` |
| 21 | Classify LOH regions | `src/LOH.PTA.sh` |
| 22 | Mirrored-BAF & AUC validation | `src/BAF_AUC.py` |

## E. Somatic variant calling

| # | Step | Script |
|---|------|--------|
| 23 | Mutect2 multisample, matched normal (per interval) | `src/Mutect2MultiSample.intervals.PTA.sh` |
| 24 | FilterMutectCalls | `src/Mutect2MultiSample.Filter.intervals.PTA.sh` |
| 25 | Concatenate somatic VCFs | `src/ConcatenateVCFs.mutect.PTA.sh` |
| 26 | Intersect somatic ∩ germline → phylogeny SNV set | `src/IntersectSomaticGermline.PTA.sh` |
| 27 | Annotate somatic variants (ANNOVAR) | `src/AnnotateSomaticVariants.PTA.sh` |



## H. Phylogeny & mutation mapping

| # | Step | Script |
|---|------|--------|
| 28 | CellPhy ML tree (GT10+FO+E, 100 bootstraps) | `src/Cellphy.PTA.sh` |
| 29 | Map mutations to tree (CellPhy/treemut) | `src/Cellphy.mutationmapping.PTA.sh` |
| 30 | Rescue filtered-out Mutect2 variants by tree congruence (in vitro only) | `src/PreMutationMapping.PTA.sh` |
| 31 | Mutation mapping (treemut) | `src/MutationMapping.PTA.sh` → `src/R/MutationMappingCellPhyTreeMut.R`, `src/R/Treemut.R` |
| 32 | Mutation mapping, second pass | `src/Mutationmapping2.PTA.sh`; chrX: `src/Mutationmapping2.chrX.PTA.sh` |
| 33 | Per-bootstrap mapping + ultrametric trees (path analysis) | `src/MutationMappingForPathBootstrap.PTA.sh`, `src/ForPathCreateUltrametric.PTA.sh` |
| 34 | Map exome / CNV / rearrangements onto tree | `src/Mutationmapping.ExomeCNV.PTA.sh` → `src/R/MutationMappingLastStep.R` |
| 35 | Write VCF of mapped sites + branch labels | `src/GetVcfAfterMutationMapping.PTA.sh` |

## I. Signatures

| # | Step | Script |
|---|------|--------|
| 36 | Trinucleotide context per sample | `src/GetTriNucContextFromVCFsPerSample.PTA.sh` |
| 37 | Trinucleotide context per branch (internal/external) | `src/GetTriNucContextFromVCFsPerSamplePerBranch.PTA.sh` |

## J. Benchmarking & accuracy

| # | Step | Script |
|---|------|--------|
| 41 | Ground-truth SNV set (TP/FP) | `src/CreateGroundTruthSetOfSNVs.PTA.sh` |
| 42 | Specificity — invariable / non-variable sites | `src/InvariableSites.PTA.a.sh`, `src/MeasureAccuracy.PTA.sh` |
| 43 | Specificity — Invitro | `src/InvariableSites.PTA.invitro.sh`, `src/MeasureAccuracy.PTA.invitro.sh` |

## F. Exome analysis

| # | Step | Script |
|---|------|--------|
| 44 | Exome Mutect2 (optimized) | `src/Mutect2MultiSample.exome.optimized.PTA.sh`, `…Filter…`, `src/ConcatenateVCFs.mutect.exome.optimized.PTA.sh` |
| 45 | Annotate exome calls | `src/AnnotateSomaticVariants.exome.PTA.sh` 
| 46 | Extra exome filtering | `src/ExtraFilterExome.PTA.sh` → `src/R/ExtraFilterExome.R` |
| 47 | VAFs at exome calls (WGS & WES) | `src/CollectAllelicCountsExome.PTA.sh`, `…WGSinExome…`, `…WESinExome…`, `…HetSNPsinExome…` |
| 48 | Intersect WGS ∩ WES | `src/IntersectWGSWES.PTA.sh` |
| — | Stand-alone exome.pta workflow | `src/exome.pta/` (`CalculateContaminationSingleCell.sh`, `InsertSizeExome.sh`, `mgatk.sh`, `Mutect2MultiSample.sh`, `Mutect2FilterMultiSample.sh`, `Cellphy.sh`) |


## K. Other analyses

| Step | Script |
|------|--------|
| Immunoglobulin (MiXCR) | `src/MiXCR.PTA.sh`, `src/MiXCR.exome.PTA.sh` |
| Copy-number (Ginkgo) | `src/PrepareGinkgo.PTA.sh`, `src/Ginkgo.PTA.sh` 

## L. Illumina precision benchmarking

Produces `variants_table_cells.rds` — the precision input for Figure 4D / Figure S6C.
All scripts live in `src/illumina-precision/`.


| # | Step | Script |
|---|------|--------|
| L1 | Expand somatic BEDs + intersect het-SNPs (external cells A1–H2) | `src/illumina-precision/GetCalls.sh` |
| L2 | Expand somatic BEDs + intersect het-SNPs (internal branches) | `src/illumina-precision/internal/GetCalls.sh` |
| L3 | Count phased reads + classify TP/FP (external cells — array job) | `src/illumina-precision/CountandClassifyIlluminaReads.sh` |
| L4 | Count phased reads + classify TP/FP (internal branches — array job) | `src/illumina-precision/internal/CountandClassifyReads.sh` |
| L5 | Concatenate all classified files → `ClassifiedVariantsTable.out` | `src/illumina-precision/GetClassifiedVariantsTable.sh` |
| L6 | Join classifications + coverage → `variants_table_cells.rds` | `src/illumina-precision/BuildVariantsTableCells.R` |
