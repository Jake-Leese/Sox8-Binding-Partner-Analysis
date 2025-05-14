# 250425 script to perform differential expression analyses between HH6 cell types
# written using alexthiery-schelper-archr_dev_macs2-schelper-0.3.5.img container

.libPaths("/R/libs/AT_ArcR_macs2")
setwd('/data/Sox8_binding_partner_analysis/scRNAseq_objects')

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load and check the Seurat object
HHall_ectoderm <- readRDS("seurat_label_transfer.RDS")

# Subset HH6 cells only
HH6_ectoderm <- HHall_ectoderm[,HHall_ectoderm$stage %in% "HH6"]
unique(HH6_ectoderm$scHelper_cell_type) # Checking cell type labels

# Add another cel meta.data column to label cells as either neural or placode
HH6_ectoderm$ectoderm_type <- NA

# Assign "neural" based on scHelper_cell_type
HH6_ectoderm$ectoderm_type[HH6_ectoderm$scHelper_cell_type %in% c("aNP", "iNP", "eCN")] <- "neural"

# Assign "placode" based on scHelper_cell_type
HH6_ectoderm$ectoderm_type[HH6_ectoderm$scHelper_cell_type %in% c("PPR", "pPPR", "aPPR")] <- "placode"

# Perform a differential expression aalysis between neural and placodal cell types at HH6
HH6_markers_placode_vs_neural <- FindMarkers(
  HH6_ectoderm,
  group.by = "ectoderm_type",
  ident.1 = "placode",
  ident.2 = "neural",
  assay = "integrated",
  min.pct = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)

# add gene names as a column
HH6_markers_placode_vs_neural$gene <- rownames(HH6_markers_placode_vs_neural)

# Add significance column for colouring points that are beyond set logFC and pval thresholds
HH6_markers_placode_vs_neural <- HH6_markers_placode_vs_neural %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up in placode",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Up in neural",
      TRUE ~ "Not significant"
    )
  )

# Idenitify top DEGs for PPR and Neural to label
PPR_top_pval_genes <- HH6_markers_placode_vs_neural %>% filter(avg_log2FC > 0.25) %>% arrange(p_val_adj) %>% slice_head(n = 10)
PPR_top_fc_genes <- HH6_markers_placode_vs_neural %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10)
neural_top_pval_genes <- HH6_markers_placode_vs_neural %>% filter(avg_log2FC < -0.25) %>% arrange(p_val_adj) %>% slice_head(n = 10)
neural_top_fc_genes <- HH6_markers_placode_vs_neural %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% slice_head(n = 10)

# Combine and remove duplicates
top_PPR_genes <- bind_rows(PPR_top_pval_genes, PPR_top_fc_genes) %>% distinct(gene) %>% pull(gene)
top_neural_genes <- bind_rows(neural_top_pval_genes, neural_top_fc_genes) %>% distinct(gene) %>% pull(gene)

# Add label column and give top PPR and neural genes TRUE values for labelling in ggplot
HH6_markers_placode_vs_neural$label <- HH6_markers_placode_vs_neural$gene %in% c(top_PPR_genes, top_neural_genes, "BMI1", "EYA2", "GATA3", "JARID2", "TFAP2C") # Labelling differentially expressed chromatin remodelling factors also

# Plot volcano plotwith top placodal and neural genes labelled
ggplot(HH6_markers_placode_vs_neural, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = subset(HH6_markers_placode_vs_neural, label == TRUE),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,  # Allow as many labels as needed
    force = 1,           # Push labels away from each other
    segment.color = "black",  # Color of leader lines
    segment.size = 0.3        # Thickness of leader lines
  ) +
  scale_color_manual(values = c("Up in neural" = "steelblue", "Up in placode" = "firebrick"),
                     name = NULL) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: HH6 Placode vs Neural",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(text = element_text(size = 14))

# Save gene lists for all PPR and neural DEGs
PPR_DEGs <- HH6_markers_placode_vs_neural %>% filter(significance == "Up in placode")
Neural_DEGs <- HH6_markers_placode_vs_neural %>% filter(significance == "Up in neural")

write.csv(PPR_DEGs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH6_PPR_DEGs.csv")
write.csv(Neural_DEGs, "/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH6_neural_DEGs.csv")

writeLines(PPR_DEGs$gene , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH6_PPR_DEGs.txt")
writeLines(Neural_DEGs$gene , "/data/Sox8_binding_partner_analysis/scRNAseq_objects/DE_analyses/HH6_neural_DEGs.txt")


