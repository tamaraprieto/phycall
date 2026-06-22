# Phycall

Somatic evolution of B-cell acute lymphoblastic leukaemia (B-ALL) and a DLD-1 in vitro clonal-evolution model, reconstructed from primary template amplification (PTA) single-cell whole-genome sequencing.
This repository holds the analysis pipeline and the code that reproduces manuscript Figures 4, 5, S5 and S6.

<img src="pipeline_graphic_abstract.png" alt="Description" width="900">


---
 
## Start here

| If you want to ... | Read |
|:----------------|:-----|
| Run the pipeline | [`docs/PIPELINE.md`](docs/PIPELINE.md) |
| Reproduce a figure | [`docs/FIGURES.md`](docs/FIGURES.md) |

---


## Analyses included

1. **Variant Calling and Analysis** — read processing (trim, map, dedup, BQSR),
   germline calling (HaplotypeCaller → GenotypeGVCFs → VQSR), somatic calling
   (Mutect2 multisample, matched normal), intersection of the two, and ANNOVAR
   annotation. → `docs/PIPELINE.md` A, B, E.

2. **Phylogenetic Reconstruction** — maximum-likelihood trees with CellPhy
   (GT10+FO+E, 100 bootstraps), mutation mapping with `treemut`, time calibration with
   `rtreefit`. → `docs/PIPELINE.md` H.

3. **Signature analysis** — 96-channel trinucleotide context, COSMIC signature fitting
   with deconstructSigs. → `docs/PIPELINE.md` I.


4. **Copy Number Variation Analysis** — Shapeit phasing → Chisel allele counts → BAF →
   LOH classification with haplotype-switch minimization. → `docs/PIPELINE.md` D.

5. **In Vitro Evolution Validation** — the DLD-1 clonal-evolution experiment
   (`Invitro` sample set).
   → `docs/PIPELINE.md`.

---

## Author

Tamara Prieto, Landau Lab, NYGC
