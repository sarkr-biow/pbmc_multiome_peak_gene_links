# Do enhancers regulate their nearest gene?

Testing how often the "nearest gene" rule agrees with a correlation-based
peak-to-gene assignment, using paired RNA and chromatin accessibility measured
in the same single cells.

## Question

When we find an accessible regulatory region in the genome, we usually don't know
which gene it controls. The standard shortcut is to assume it regulates the
closest gene. This analysis asks how often that shortcut gives the same answer as
a method that uses the data itself: correlating peak accessibility against gene
expression across thousands of cells.

## Data

Public 10x Genomics PBMC 10k Multiome dataset (paired scRNA-seq and scATAC-seq
from the same nuclei), aligned to hg38.

`https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/`

## Workflow

1. **Load** RNA and ATAC into a single Seurat object sharing one set of cell
   barcodes.
2. **Quality control** on TSS enrichment, nucleosome signal, and per-cell counts
   in both assays.
3. **Process each modality separately.** SCTransform + PCA for RNA; TF-IDF + SVD
   (LSI) for ATAC. LSI rather than PCA for ATAC because the data is sparse and
   close to binary, which breaks PCA's assumptions.
4. **Joint embedding with WNN,** which learns a per-cell weighting between the two
   modalities rather than blending them at a fixed ratio. LSI component 1 is
   excluded because it tracks sequencing depth rather than biology.
5. **Cell type annotation** by marker expression, read off a per-cluster dot plot
   (`figures/03a_marker_dotplot_clusters.png`).
6. **Peak-to-gene linkage** with `LinkPeaks()`, correlating peak accessibility
   against gene expression within a 500 kb window, correcting for GC content,
   peak width, and overall accessibility.
7. **Nearest-gene comparison.** For each linked peak, find the closest
   protein-coding TSS and ask whether it matches the gene the correlation chose.

## Results

11,909 cells were loaded and 11,070 passed quality filtering. 22 clusters were
annotated into 9 cell types: CD14 Mono, CD4 T, CD8 T, CD16 Mono, NK, B, DC, Other T, Unknown.

`LinkPeaks` was run on 3,000 variable genes and returned **7,304 significant
peak-gene links** across 2,481 testable genes.

### Nearest gene vs. linked gene

| Nearest-gene search space | Agreement |
|---|---|
| All 19,919 protein-coding genes | **40%** (2,923 / 7,304) |
| Only the 2,481 genes `LinkPeaks` could test | **63%** (4,571 / 7,304) |

**Distal peaks only (>10 kb from the linked TSS):** 6,042 links, 32% agreement
against all protein-coding genes and 59% against testable genes only.

### Distance distribution

`figures/04_link_distances.png` shows a sharp enrichment of links within roughly
50 kb of the TSS, decaying out to about 150 kb, then a flat tail extending all
the way to the 500 kb search boundary.

### Worked example

`figures/06_MS4A1_links.png` shows the MS4A gene cluster on chromosome 11. Peaks
linked to MS4A1 span roughly 200 kb and include regions near MS4A6E and MS4A7,
reaching across several intervening MS4A family genes. Accessibility across the
region is highest in B cells, matching where MS4A1 is expressed.

## Interpretation

Under the broadest comparison, the nearest protein-coding gene agreed with the
correlation-based link 40% of the time. **That number overstates the case, and
the reason matters.**

**The two searches were not symmetric.** `LinkPeaks` ran on the variable gene set
only, so it could never select a nearest gene outside that set - but the
nearest-gene search covered all 19,919 protein-coding genes. Restricting both
sides to the same gene universe raises agreement from 40% to 63%. Roughly half
the original disagreement was an artifact of that mismatch rather than biology.

**The remaining 37% is the more credible finding:** cases where the correlation
selected a gene other than the nearest one it could have selected.

**The long-distance tail is probably contaminated.** Real enhancer-promoter
interactions become less likely with distance, so a flat distribution out to the
edge of the search window is more consistent with false positives than with
genuine long-range regulation. Highly expressed housekeeping genes are a likely
source: RPL22, a ribosomal protein gene, was linked to peaks over 240 kb away,
which is more plausibly explained by co-variation with overall cell state than by
a regulatory relationship.

**Correlation is not contact.** Two features can co-vary without ever touching -
both may be responding to the same upstream signal. This analysis cannot
distinguish those cases, and that is the central limitation of the approach.

## Limitations

- Only variable genes were submitted to `LinkPeaks`. Both nearest-gene
  comparisons are reported because of this asymmetry.
- Only protein-coding genes were used as nearest-gene candidates. Including
  non-coding genes would change the baseline.
- Links are candidates, not validated regulatory relationships. No physical
  contact data was used to confirm any of them.
- The flat long-distance component of the distance distribution implies an
  unquantified false positive rate.
- Cell types are assigned from marker expression rather than reference mapping,
  so they are candidate identities. Clusters that could not be confidently
  assigned are labelled "Unknown". Platelet markers (PPBP, PF4) were not
  detected, so platelets are not identified here.
- One dataset, one tissue. Nothing here has been shown to generalize.

## Why the question matters

Most disease-associated genetic variants fall outside genes, in regulatory
regions, so interpreting them depends on knowing which gene each region controls.
Linear proximity is a weak proxy for that, and correlation across cells is a
better but still indirect one - it measures co-variation, not physical
interaction. Methods that measure chromatin contact directly, such as HiCAR
(Wei et al., *Molecular Cell*, 2022) and scHiCAR (Wei et al., *Nature
Biotechnology*, 2026), which profiles RNA, accessibility, and chromosome
conformation in the same cell, exist because both proxies fall short. This
analysis is an attempt to see where that shortfall shows up in real data.

## Repository contents

- `multiome_analysis.R` - the full analysis
- `peak_gene_links.csv` - every link with its nearest gene, nearest testable
  gene, and distance to the linked TSS

**Figures**

| File | Shows |
|---|---|
| `01_quality_control.png` | QC metrics: per-cell counts, TSS enrichment, nucleosome signal |
| `02_wnn_umap.png` | Joint RNA + ATAC UMAP, annotated |
| `03_marker_genes.png` | Marker expression on the joint UMAP |
| `03a_marker_dotplot_clusters.png` | Marker expression per cluster - the basis for annotation |
| `03b_marker_dotplot_labelled.png` | The same markers after labelling, as a check |
| `04_link_distances.png` | Distance from each peak to its linked gene TSS |
| `05_nearest_gene_accuracy.png` | Nearest-gene agreement under both search spaces |
| `06_MS4A1_links.png` | Accessibility, expression, and links across the MS4A cluster |

Both dot plots are included deliberately: `03a` shows the per-cluster evidence
the annotation was based on, and `03b` confirms the labels match the markers
after assignment.

## Reproducing

Requires R with Signac, Seurat, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38,
and ggplot2. The script downloads its own input data (~10 GB) and creates the
`data/` and `figures/` directories. `LinkPeaks` is the slow step and takes
1-3 hours. Roughly 16 GB of RAM is recommended.

Run under R 4.6.0 with Seurat 5.5.1 and Signac 1.17.1.

## Credit

The processing workflow - QC, LSI, WNN, and `LinkPeaks` - follows the Signac
joint RNA and ATAC vignette from the Stuart Lab. The nearest-gene comparison,
the testable-gene control, the distance analysis, and the interpretation are my
own additions.

Tools: Signac, Seurat, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38, ggplot2.
