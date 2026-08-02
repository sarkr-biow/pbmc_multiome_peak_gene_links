library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

dir.create("~/multiome-links", showWarnings = FALSE)
setwd("~/multiome-links")
dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

base_url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/"

files_to_get <- c(
  "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
)

options(timeout = 7200)

for (f in files_to_get) {
  if (!file.exists(file.path("data", f))) {
    download.file(paste0(base_url, f), file.path("data", f), mode = "wb")
  }
}

counts <- Read10X_h5("data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
frag_path <- "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag_path,
  annotation = annotation
)

n_cells_before <- ncol(pbmc)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

qc_plot <- VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

ggsave("figures/01_quality_control.png", qc_plot, width = 16, height = 5, dpi = 200)

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

DefaultAssay(pbmc) <- "RNA"

pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)


pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

pbmc <- FindClusters(
  pbmc,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)


saveRDS(pbmc, "data/pbmc_processed.rds")

pbmc <- readRDS("data/pbmc_processed.rds")

genes_to_test <- VariableFeatures(pbmc, assay = "SCT")
cat("Genes submitted to LinkPeaks:", length(genes_to_test), "\n")

DefaultAssay(pbmc) <- "ATAC"

pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = genes_to_test
)

saveRDS(pbmc, "data/pbmc_with_links.rds")

all_genes <- genes(EnsDb.Hsapiens.v86)

all_genes <- all_genes[all_genes$gene_biotype == "protein_coding"]

gene_table <- data.frame(
  gene_name = all_genes$gene_name,
  chr       = paste0("chr", as.character(seqnames(all_genes))),
  start     = start(all_genes),
  end       = end(all_genes),
  strand    = as.character(strand(all_genes)),
  stringsAsFactors = FALSE
)

head(gene_table)
nrow(gene_table)

gene_table$tss <- 0

for (i in 1:nrow(gene_table)) {
  
  if (gene_table$strand[i] == "+") {
    gene_table$tss[i] <- gene_table$start[i]
  } else {
    gene_table$tss[i] <- gene_table$end[i]
  }
  
}

normal_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
gene_table <- gene_table[gene_table$chr %in% normal_chromosomes, ]

link_table <- as.data.frame(Links(pbmc))


link_table$nearest_gene <- NA
link_table$distance_to_linked_gene <- NA

for (i in 1:nrow(link_table)) {
  peak_name  <- link_table$peak[i]
  pieces     <- strsplit(peak_name, "-")[[1]]
  peak_chr   <- pieces[1]
  peak_start <- as.numeric(pieces[2])
  peak_end   <- as.numeric(pieces[3])
  peak_middle <- (peak_start + peak_end) / 2
  
  
  genes_here <- gene_table[gene_table$chr == peak_chr, ]
  
  if (nrow(genes_here) == 0) {
    next
  }
  
  smallest_distance <- abs(genes_here$tss[1] - peak_middle)
  closest_gene_name <- genes_here$gene_name[1]
  
  for (j in 1:nrow(genes_here)) {
    
    this_distance <- abs(genes_here$tss[j] - peak_middle)
    
    if (this_distance < smallest_distance) {
      smallest_distance <- this_distance
      closest_gene_name <- genes_here$gene_name[j]
    }
    
  }
  
  link_table$nearest_gene[i] <- closest_gene_name
  
  
  linked_gene_name <- link_table$gene[i]
  
  linked_row <- genes_here[genes_here$gene_name == linked_gene_name, ]
  
  if (nrow(linked_row) > 0) {
    link_table$distance_to_linked_gene[i] <- abs(linked_row$tss[1] - peak_middle)
  }
  
}

link_table$nearest_was_right <- link_table$nearest_gene == link_table$gene
mean(link_table$nearest_was_right, na.rm = TRUE)


dist_plot <- ggplot(link_table, aes(x = distance_to_linked_gene / 1000)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  labs(
    x = "Distance from peak to linked gene TSS (kb)",
    y = "Number of links",
    title = "Most linked peaks are not sitting on the gene"
  ) +
  theme_minimal()

ggsave("figures/04_link_distances.png", dist_plot, width = 7, height = 5, dpi = 200)



wrong_ones <- link_table[
  !is.na(link_table$nearest_was_right) &
    link_table$nearest_was_right == FALSE &
    link_table$distance_to_linked_gene > 50000,
]


Idents(pbmc) <- pbmc$seurat_clusters

DefaultAssay(pbmc) <- "SCT"

marker_list <- c(
  "CD3D", "CD3E",           # T cells
  "IL7R", "CCR7",           # CD4 T
  "CD8A", "CD8B",           # CD8 T
  "NKG7", "GNLY", "KLRD1",  # NK / cytotoxic
  "MS4A1", "CD79A",         # B
  "LYZ", "CD14",            # CD14 monocytes
  "FCGR3A", "MS4A7",        # CD16 monocytes
  "FCER1A", "CST3"          # dendritic cells
)

dot <- DotPlot(pbmc, features = marker_list) + RotatedAxis()
ggsave("figures/03a_marker_dotplot_clusters.png", dot, width = 13, height = 8, dpi = 200)

testable_genes <- intersect(genes_to_test, gene_table$gene_name)
length(testable_genes)

gene_table_testable <- gene_table[gene_table$gene_name %in% testable_genes, ]

link_table$nearest_testable_gene <- NA

for (i in 1:nrow(link_table)) {
  
  pieces      <- strsplit(link_table$peak[i], "-")[[1]]
  peak_chr    <- pieces[1]
  peak_middle <- (as.numeric(pieces[2]) + as.numeric(pieces[3])) / 2
  
  genes_here <- gene_table_testable[gene_table_testable$chr == peak_chr, ]
  if (nrow(genes_here) == 0) next
  
  smallest_distance <- abs(genes_here$tss[1] - peak_middle)
  closest_gene_name <- genes_here$gene_name[1]
  
  for (j in 1:nrow(genes_here)) {
    this_distance <- abs(genes_here$tss[j] - peak_middle)
    if (this_distance < smallest_distance) {
      smallest_distance <- this_distance
      closest_gene_name <- genes_here$gene_name[j]
    }
  }
  
  link_table$nearest_testable_gene[i] <- closest_gene_name
}

link_table$nearest_testable_was_right <-
link_table$nearest_testable_gene == link_table$gene
  
cat("Nearest gene correct (all protein-coding):",
      mean(link_table$nearest_was_right, na.rm = TRUE), "\n")
cat("Nearest gene correct (testable genes only):",
      mean(link_table$nearest_testable_was_right, na.rm = TRUE), "\n")
  
distal_links <- link_table[
    !is.na(link_table$distance_to_linked_gene) &
      link_table$distance_to_linked_gene > 10000,
    ]
  
cat("Distal links (>10 kb):", nrow(distal_links), "\n")
cat("Distal, all protein-coding:",
    mean(distal_links$nearest_was_right, na.rm = TRUE), "\n")
cat("Distal, testable genes only:",
    mean(distal_links$nearest_testable_was_right, na.rm = TRUE), "\n")

number_right <- sum(link_table$nearest_was_right, na.rm = TRUE)
number_wrong <- sum(!link_table$nearest_was_right, na.rm = TRUE)

summary_table <- data.frame(
  comparison = rep(c("All protein-coding genes", "Testable genes only"), each = 2),
  answer     = rep(c("Nearest gene was right", "Nearest gene was wrong"), 2),
  count      = c(
    sum(link_table$nearest_was_right,           na.rm = TRUE),
    sum(!link_table$nearest_was_right,          na.rm = TRUE),
    sum(link_table$nearest_testable_was_right,  na.rm = TRUE),
    sum(!link_table$nearest_testable_was_right, na.rm = TRUE)
  )
)

bar_plot <- ggplot(summary_table, aes(x = answer, y = count, fill = answer)) +
  geom_col() +
  facet_wrap(~ comparison) +
  labs(x = NULL, y = "Number of peak-gene links") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/05_nearest_gene_accuracy.png", bar_plot, width = 9, height = 5, dpi = 200)

write.csv(link_table, "peak_gene_links.csv", row.names = FALSE)

Idents(pbmc) <- pbmc$seurat_clusters

pbmc <- RenameIdents(pbmc,
                     "0"  = "CD14 Mono",
                     "1"  = "CD4 T",
                     "2"  = "CD8 T",
                     "3"  = "CD4 T",
                     "4"  = "CD4 T",
                     "5"  = "CD16 Mono",
                     "6"  = "NK",
                     "7"  = "CD14 Mono",
                     "8"  = "CD8 T",
                     "9"  = "B",
                     "10" = "B",
                     "11" = "CD8 T",
                     "12" = "CD8 T",
                     "13" = "CD4 T",
                     "14" = "CD4 T",
                     "15" = "DC",
                     "16" = "Other T",
                     "17" = "Unknown",
                     "18" = "B",
                     "19" = "Unknown",
                     "20" = "Unknown",
                     "21" = "Unknown"
)

pbmc$cell_type <- Idents(pbmc)
table(pbmc$cell_type)

DefaultAssay(pbmc) <- "SCT"

umap_labelled <- DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave("figures/02_wnn_umap.png", umap_labelled, width = 7, height = 6, dpi = 200)

marker_plot <- FeaturePlot(
  pbmc,
  features = c("MS4A1", "CD3D", "CD8A", "LYZ", "NKG7", "FCGR3A"),
  reduction = "wnn.umap",
  ncol = 3
)

ggsave("figures/03_marker_genes.png", marker_plot, width = 12, height = 8, dpi = 200)

dot_check <- DotPlot(pbmc, features = c("CD3D", "CD8A", "NKG7", "MS4A1", "LYZ", "FCGR3A", "FCER1A")) +
  RotatedAxis()

ggsave("figures/03b_marker_dotplot_labelled.png", dot_check, width = 9, height = 6, dpi = 200)

print(link_table[link_table$gene == "MS4A1",
                 c("peak", "gene", "score", "distance_to_linked_gene")])

DefaultAssay(pbmc) <- "ATAC"

cov_plot <- CoveragePlot(
  object            = pbmc,
  region            = "MS4A1",
  features          = "MS4A1",
  expression.assay  = "SCT",
  idents            = c("B", "CD4 T", "CD8 T", "NK", "CD14 Mono", "CD16 Mono", "DC"),
  extend.upstream   = 150000,
  extend.downstream = 150000
)

ggsave("figures/06_MS4A1_links.png", cov_plot, width = 10, height = 7, dpi = 200)

print(head(wrong_ones[, c("peak", "gene", "nearest_gene", "distance_to_linked_gene")]))

saveRDS(pbmc, "data/pbmc_final.rds")
write.csv(link_table, "peak_gene_links.csv", row.names = FALSE)

sessionInfo()
