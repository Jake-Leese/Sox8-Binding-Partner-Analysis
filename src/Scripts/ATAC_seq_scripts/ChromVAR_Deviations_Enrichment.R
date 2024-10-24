ArchRProject <- addMotifAnnotations(ss8_Sox8_ArchRProj, name = "Motif", motifPWMs = motifList, cutOff = 1e-05, force = T)

ArchRProject <- addBgdPeaks(ArchRProject)

ArchRProject <- addDeviationsMatrix(
  ArchRProject,
  peakAnnotation = "Motif",
  force = TRUE)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = ss8_Sox8_ArchRProj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp_NC <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 2,
    nudge_x = 2,
    color = "black") + 
  theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  theme(legend.position = "none")

ggUp

motifsdown <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = ss8_Sox8_ArchRProj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsdown)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp_NC <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 2,
    nudge_x = 2,
    color = "black") + 
  theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  theme(legend.position = "none")

ggUp_NC

