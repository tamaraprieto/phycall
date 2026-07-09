# Reproducing the manuscript figures

Figures 4, 5, S5, S6 (and optional Figure CNVs). 
---

## Panel → notebook/chunk mapping

| Figure | Panel | Notebook | Chunk / script |
|--------|-------|----------|----------------|
| 4 | A | `Figures.BALL.PTA.Rmd` | `# Fig. 4A` — `ESTIMATED TREE` |
| 4 | B | `Figures.BALL.PTA.Rmd` | `# Fig. 4B` — `heatmap similarity` |
| 4 | C | `Figures.BALL.PTA.Rmd` | `# Fig. 4C` — `Correlation branch lengths and sampling times` |
| 4 | D | `Figures.BALL.PTA.Rmd` | `# Fig. 4D` — `Precision/Recall invitro` |
| 4 | E | `Figures.BALL.PTA.Rmd` | `# Fig. 4E` — `Signature invitro tree` |
| 4 | F | `Figures.BALL.PTA.Rmd` | `# Fig. 4F` — `Deletion at chr13:41,875,000-42,400,000` ¹ |
| 4 | G | `src/BAF_AUC.py` | — (output to `src/Invitro_hg38/`) |
| 5 | A | illustrator | not code |
| 5 | B | `Figures.BALL.PTA.Rmd` | `number of mutations per branch` |
| 5 | C | `Figures.BALL.PTA.Rmd` | mutational-pattern chunk |
| 5 | D | `Figures.BALL.PTA.Rmd` | signatures chunk |
| 5 | E | `Figures.BALL.PTA.Rmd` | `PCA` chunk |
| 5 | F | `Figures.BALL.PTA.Rmd` | `patient phylogeny 4295` |
| 5 | G | `Figures.BALL.PTA.Rmd` | `patient phylogenies small` |
| 5 | H | external repo — https://github.com/jzinno/scABC | JAK2-mutant clone fitness analysis (patient A) ² |

| 5 | I | `src/PlotTrees.R` | `RunPhylodyn()` / `BNPR()` section ³ |
| 5 | J | `analysis/LOH.analysis.Rmd` | `complexheatmap` chunk |
| 5 | K | `Figures.BALL.PTA.Rmd` | `PhylogeneticSignalPatient` |
| 5 | L | `Figures.BALL.PTA.Rmd` | `445.zeb.cd10` |
| 5 | M | `Figures.BALL.PTA.Rmd` | `4084.kras.cd10` |
| S5 | a | `Figures.BALL.PTA.Rmd` | `# Figure S5` — `mad pta-mda` |
| S5 | b | `Figures.BALL.PTA.Rmd` | `depth-breadth relationship` |
| S5 | c–d | `Figures.BALL.PTA.Rmd` | `allelic imbalance` |
| S5 | e | `Figures.BALL.PTA.Rmd` | `ado` |
| S5 | f | `Figures.BALL.PTA.Rmd` | `recall by coverage` |
| S5 | g | `Figures.BALL.PTA.Rmd` | erroneous-reads chunk |
| S6 | A | illustrator | not code |
| S6 | B | `Figures.BALL.PTA.Rmd` | `# Fig. S6B` — `patient phylogeny` chunks |
| S6 | C | `Figures.BALL.PTA.Rmd` | `# Fig S6C` — `F1.invitro.comparison` |
| S6 | D | `Figures.BALL.PTA.Rmd` | `# Fig S6D` — `ComplexHeatmap.Exome` |
| S6 | E–F | `Figures.BALL.PTA.Rmd` | `# Fig 5E, Fig S6E-F` — `PCA` chunk |
| S6 | G | `Figures.BALL.PTA.Rmd` | effective-population-size chunk |
| CNVs | — | `analysis/LOH.analysis.Rmd` | `Plot states next to the phylogeny`, `Filter consecutive and by band` |

---

## Notes

**¹ Fig 4F — mirrored BAF input.** Run `src/generate_chr13_mirrored_baf.py` on the cluster first to produce `region_1_chr13_5k_imbalance.csv` from `dict_chr13.pkl` in `RESULTS/BAF/BAF_5000/`. The CSV covers chr13:40,000,000–43,000,000 in 5 kb windows; the notebook plots the sub-region 41,400,000–42,900,000.

**² Fig 5H — clonal fitness.** `analysis/ABC` (JAK2-mutant clone fitness analysis, patient A)

**³ Fig 5I — phylodynamics.** `src/PlotTrees.R` reads pre-computed ultrametric RDS files from `src/ForPathCreateUltrametric.PTA.sh`. Cluster paths:
- `timetrees/TreeMut4295_Time_Tree_Null.RDS`
- `Downloads/TreeMut417_Time_Tree_Null.RDS`
- `Downloads/TreeMut445_Time_Tree_Null.RDS`
