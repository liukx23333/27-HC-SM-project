###############################################################################
# SGMS1 KO IEL single-cell RNA-seq analysis
#这段只用了wt做不同subset的基因分析
###############################################################################

set.seed(20260729)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(svglite)
  library(ggsignif)
  library(rstatix)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(biomaRt)
  library(clusterProfiler)
  library(enrichplot)
})

# ---- Project paths -----------------------------------------------------------
project_dir <- path.expand(
  "~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/sgms1 ko iel"
)
figure_dir <- file.path(project_dir, "results", "figures")
table_dir  <- file.path(project_dir, "results", "tables")
object_dir <- file.path(project_dir, "results", "objects")
log_dir    <- file.path(project_dir, "results", "logs")

invisible(lapply(
  c(project_dir, figure_dir, table_dir, object_dir, log_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))


# ---- Reusable publication settings ------------------------------------------
sample_colors <- c(iwt = "#4F86C6", iko = "#D64065")

celltype_colors <- c(
  # CD8aa+ TCRab/gd：从原橙色 #F8845C 扩展
  "Memory-like CD8aa+ TCRab/gd T cells" =
    "#91D0B3",
  
  "Effector-like CD8aa+ TCRab/gd T cells" =
    "#3E9B72",
  
  # CD4：从原绿色 #5EB887 扩展
  "Memory-like CD4 T cells" =
    "#F8845C",
  
  "Effector-like CD4 T cells" =
    "#F6AA8D",
  "CD8ab+ T cells" =
    "#D64065",
  # Naive：保留原黄色
  "Naive T Cells" =
    "#F5D983",
  
  # Cycling：保留原红色
  "Cycling T Cells" =
    "#4F86C6",
  
  # Treg：在原红、蓝配色基础上增加紫色
  "Treg cells" =
    "#8E6BBE"
  
)

theme_publication <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}
theme_set(theme_publication())

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}


theme_publication <- function(
    base_size = 9,
    base_family = "sans"
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}

theme_set(theme_publication())

save_plot <- function(
    plot,
    filename,
    width,
    height,
    units = "in",
    dpi = 600,
    bg = "white"
) {
  stem <- file.path(
    figure_dir,
    safe_name(filename)
  )
  
  # PDF矢量图
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = units,
    bg = bg,
    limitsize = FALSE
  )
  
  # SVG矢量图
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height,
    bg = bg
  )
  print(plot)
  grDevices::dev.off()
  
  # PNG位图
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  
  invisible(plot)
}

save_table <- function(x, filename, row.names = FALSE) {
  write.csv(
    x, file.path(table_dir, filename),
    row.names = row.names, fileEncoding = "UTF-8"
  )
}

# ---- Read, merge, QC ---------------------------------------------------------
iwt_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/sgms1 ko iel/IWT/filtered_feature_bc_matrix")
iko_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/sgms1 ko iel/IKO/filtered_feature_bc_matrix")

iwt <- CreateSeuratObject(iwt_counts, project = "iwt")
iko <- CreateSeuratObject(iko_counts, project = "iko")

iel <- merge(iwt,iko,
             add.cell.ids = c("iwt","iko"),
             project = "total") 

iel$sample <- factor(iel$orig.ident, levels = c("iwt", "iko"))

# Mouse mitochondrial genes usually begin "mt-"; ignore.case covers "Mt-".
iel[["percent.mt"]] <- PercentageFeatureSet(
  iel, pattern = "^mt-", assay = "RNA"
)

p_qc_before <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_before)
save_plot(p_qc_before, "01_QC_before_filtering", 10.5, 3.5)

iel <- subset(
  iel,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA < 20000 &
    percent.mt < 5
)

p_qc_after <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_after)
save_plot(p_qc_after, "02_QC_after_filtering", 10.5, 3.5)

qc_summary <- iel[[]] |>
  group_by(sample) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
save_table(qc_summary, "QC_summary_after_filtering.csv")

# ---- Normalization, Harmony, clustering -------------------------------------
iel <- NormalizeData(iel, verbose = FALSE)
iel <- FindVariableFeatures(iel, nfeatures = 2000, verbose = FALSE)
iel <- ScaleData(iel, features = VariableFeatures(iel), verbose = FALSE)
iel <- RunPCA(iel, features = VariableFeatures(iel), verbose = FALSE)

p_elbow <- ElbowPlot(iel, ndims = 50) + theme_publication()
save_plot(p_elbow, "03_PCA_elbow", 5.5, 4)

# Harmony is run once. Downstream neighbors and UMAP explicitly use harmony.
iel <- RunHarmony(
  iel, group.by.vars = "sample", reduction = "pca",
  reduction.save = "harmony", assay.use = "RNA", verbose = FALSE
)

dims_main <- 1:20
iel <- FindNeighbors(iel, reduction = "harmony", dims = dims_main, verbose = FALSE)
iel <- FindClusters(iel, resolution = 0.3, random.seed = 20260729, verbose = FALSE)
iel <- RunUMAP(
  iel, reduction = "harmony", dims = dims_main,
  seed.use = 20260729, verbose = FALSE
)

p_clusters <- DimPlot(
  iel, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_clusters)
save_plot(p_clusters, "04_UMAP_clusters", 6, 5)

p_samples <- DimPlot(
  iel, reduction = "umap", group.by = "sample",
  cols = sample_colors, raster = FALSE
) + theme_publication()
print(p_samples)
save_plot(p_samples, "05_UMAP_samples", 6, 5)

p_split <- DimPlot(
  iel, reduction = "umap", group.by = "seurat_clusters",
  split.by = "sample", label = TRUE, repel = TRUE, raster = FALSE
) + theme_publication()
print(p_split)
save_plot(p_split, "06_UMAP_clusters_split_by_sample", 10, 5)
#############################################################
# 统计每个样本、每个 cluster 的细胞数
cluster_counts <- as.data.frame(
  table(
    sample = iel$sample,
    cluster = iel$seurat_clusters
  )
)

colnames(cluster_counts) <- c(
  "sample",
  "cluster",
  "cell_count"
)

# 计算每个样本内部各 cluster 的细胞比例
cluster_ratio <- cluster_counts |>
  group_by(sample) |>
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells
  ) |>
  ungroup()

# 保存完整的细胞数和比例
save_table(
  cluster_ratio,
  "cluster_counts_and_proportions.csv"
)

# 分别提取 IWT 和 IKO
ratio_iwt <- cluster_ratio |>
  filter(sample == "iwt") |>
  select(
    cluster,
    iwt_count = cell_count,
    iwt_proportion = proportion
  )

ratio_iko <- cluster_ratio |>
  filter(sample == "iko") |>
  select(
    cluster,
    iko_count = cell_count,
    iko_proportion = proportion
  )

# 合并 IWT 和 IKO 的结果
cluster_change <- full_join(
  ratio_iwt,
  ratio_iko,
  by = "cluster"
)

# 如果某个 cluster 只存在于一个样本中，将缺失值设为0
cluster_change[is.na(cluster_change)] <- 0

# 计算比例差异
cluster_change <- cluster_change |>
  mutate(
    # 正数表示KO中比例升高，负数表示KO中比例降低
    proportion_change = iko_proportion - iwt_proportion,
    
    # 转换成百分点，便于阅读
    percentage_point_change = proportion_change * 100,
    
    # 绝对变化幅度，用于寻找变化最大的cluster
    absolute_change = abs(proportion_change),
    
    # 比例变化倍数；加入极小值避免除以0
    log2_proportion_fold_change = log2(
      (iko_proportion + 1e-6) /
        (iwt_proportion + 1e-6)
    ),
    
    direction = case_when(
      proportion_change > 0 ~ "Higher in IKO",
      proportion_change < 0 ~ "Lower in IKO",
      TRUE ~ "No change"
    )
  ) |>
  arrange(desc(absolute_change)) |>
  mutate(
    change_rank = row_number()
  )

# 保存变化结果
save_table(
  cluster_change,
  "cluster_proportion_changes_IKO_vs_IWT.csv"
)

# 在控制台显示变化最大的cluster
print(
  cluster_change |>
    select(
      change_rank,
      cluster,
      iwt_count,
      iko_count,
      iwt_proportion,
      iko_proportion,
      percentage_point_change,
      log2_proportion_fold_change,
      direction
    )
)

# Cluster比例对比图
p_cluster_ratio <- ggplot(
  cluster_ratio,
  aes(
    x = cluster,
    y = proportion,
    fill = sample
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72
  ) +
  scale_fill_manual(values = sample_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Cluster proportions in IWT and IKO",
    x = "Cluster",
    y = "Cell proportion",
    fill = "Sample"
  ) +
  theme_publication()

print(p_cluster_ratio)

save_plot(
  p_cluster_ratio,
  "06A_cluster_proportions_IWT_vs_IKO",
  width = 8,
  height = 5
)


# 按变化幅度排序
cluster_change$cluster_ordered <- factor(
  cluster_change$cluster,
  levels = rev(cluster_change$cluster)
)

p_cluster_change <- ggplot(
  cluster_change,
  aes(
    x = cluster_ordered,
    y = percentage_point_change,
    fill = direction
  )
) +
  geom_col(width = 0.72) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Higher in IKO" = "#D64065",
      "Lower in IKO" = "#4F86C6",
      "No change" = "#D2DAE2"
    )
  ) +
  labs(
    title = "Cluster proportion changes: IKO vs IWT",
    x = "Cluster",
    y = "Change in IKO relative to IWT\n(percentage points)",
    fill = NULL
  ) +
  theme_publication()

print(p_cluster_change)

save_plot(
  p_cluster_change,
  "06B_cluster_proportion_changes_IKO_vs_IWT",
  width = 7,
  height = 5.5
)


# =============================================================================
# 计算每个 cluster 的特征基因
# =============================================================================

# 建立单独对象，避免改变主对象iel
iel_markers <- iel

DefaultAssay(iel_markers) <- "RNA"

# Seurat v5中，差异分析前合并RNA assay的不同layers
if (inherits(iel_markers[["RNA"]], "Assay5")) {
  iel_markers <- JoinLayers(
    iel_markers,
    assay = "RNA"
  )
}

# 将细胞身份设置为cluster
Idents(iel_markers) <- "seurat_clusters"

# 计算每个cluster相对于其他全部细胞的上调特征基因
cluster_markers <- FindAllMarkers(
  object = iel_markers,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# 判断Seurat版本使用的fold-change列名
fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cluster_markers)
)[1]

if (is.na(fc_column)) {
  stop("FindAllMarkers结果中没有找到fold-change列。")
}

# 按cluster、显著性和fold change排序
cluster_markers <- cluster_markers |>
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]])
  )

# 保存所有cluster的全部marker
save_table(
  cluster_markers,
  "all_cluster_marker_genes.csv"
)


top15_cluster_markers <- cluster_markers |>
  filter(
    p_val_adj < 0.05
  ) |>
  group_by(cluster) |>
  slice_max(
    order_by = .data[[fc_column]],
    n = 15,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(
    cluster,
    desc(.data[[fc_column]])
  )

save_table(
  top15_cluster_markers,
  "top15_marker_genes_per_cluster.csv"
)

#==============================================================

marker_panel <- c(
  "Trac", "Tcrg-V7", "Tcrg-C1", "Trdc","Cd8a","Cd4","Cd8b1","Il17a","Foxp3","Tcf7","Ms4a4b","Il7r","Xcl1","Gm8369","Cd244a","Klf3","Ccl4","Cd200r1","Themis",
  "Cd6","Gzmk","Ccl5","Gzmb","Gzma","Ifng","Mki67","Ccnb2","Stmn1","Foxp3","Cd40lg","Maf","Ctla4","Icos",
  "Tyrobp","Fcer1g","Klrd1","Cd200r2","Frmd5","Cd160","Kit","Id3",
  "Emid1","Pclaf","S1pr1","Lef1","Cd19","Klrb1c","Ncr1","Csf1r","Adgre1","Mpeg1","C1qa","Lyz2"
)
marker_panel <- intersect(marker_panel, rownames(iel))
p_markers <- DotPlot(iel, features = marker_panel) +
  coord_flip() + RotatedAxis() + theme_publication()
print(p_markers)
save_plot(p_markers, "07_cluster_marker_dotplot", 10, 8)

#260729 19:00rds
saveRDS(iel,file = "260807 sgms1 ko iel dimplot.rds")
# ---- Major cell-type annotation ---------------------------------------------
# IMPORTANT: verify these cluster-to-cell-type mappings from marker plots before
# biological interpretation. Named assignment avoids fragile numeric column use.
major_map <- c(
  "0"  = "Memory-like CD8aa+ TCRab/gd T cells",
  "1"  = "CD8ab+ T cells",
  "2"  = "Effector-like CD8aa+ TCRab/gd T cells",
  "3"  = "Effector-like CD4 T cells",
  "4"  = "Effector-like CD8aa+ TCRab/gd T cells",
  "5"  = "Naive T Cells",
  "6"  = "Memory-like CD4 T cells",
  "7"  = "Effector-like CD8aa+ TCRab/gd T cells",
  "8"  = "Cycling T Cells",
  "9"  = "Treg cells",
  "10" = "CD8ab+ T cells",
  "11" = "Monocyte",
  "12" = "Naive T Cells"
)
iel$celltype <- unname(major_map[as.character(iel$seurat_clusters)])


saveRDS(iel,file = "260803 sgms1 ko iel dimplot.rds")
iel_orig<-iel
iel<-iel_orig
iel <- subset(
  iel,
  subset = celltype != "Monocyte"
)

Idents(iel) <- "celltype"
celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd T cells",
  "Effector-like CD8aa+ TCRab/gd T cells",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T Cells",
  "Treg cells",
  "Cycling T Cells"
)

# 如果已经删除Monocyte，这里不再加入Monocyte
iel$celltype <- factor(
  as.character(iel$celltype),
  levels = celltype_order
)

Idents(iel) <- "celltype"

p_celltypes <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",cols = celltype_colors,
  repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_celltypes)
save_plot(p_celltypes, "08_UMAP_major_celltypes", 7, 5.5)

p_celltypes_split <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",
  split.by = "sample", cols = celltype_colors, repel = TRUE, raster = FALSE
) + theme_publication()
save_plot(p_celltypes_split, "09_UMAP_celltypes_split_by_sample", 11, 5.5)

cell_counts <- as.data.frame(table(
  celltype = iel$celltype, sample = iel$sample
)) |>
  filter(Freq > 0) |>
  group_by(sample) |>
  mutate(ratio = Freq / sum(Freq)) |>
  ungroup()
print(cell_counts)
save_table(cell_counts, "major_celltype_counts_and_ratios.csv")

p_ratio <- ggplot(cell_counts, aes(sample, ratio, fill = celltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell proportion", fill = "Cell type") +
  theme_publication()
print(p_ratio)
save_plot(p_ratio, "10_major_celltype_proportions", 6.5, 5)



# ---- Sgms1 and sphingomyelin-synthesis analysis in WT -----------------------
objWT<-subset(x = iel, subset = orig.ident == "iwt")
objWT1<-subset(x = objWT, idents=c("Memory-like CD8aa+ TCRab/gd T cells",
                                   "Effector-like CD8aa+ TCRab/gd T cells",
                                   "CD8ab+ T cells",
                                   "Memory-like CD4 T cells",
                                   "Effector-like CD4 T cells"))



expr_sgms1 <- FetchData(
  objWT1,
  vars = c(
    "Sgms1",
    "celltype",
    "sample"
  ),
  layer = "data"
)



expr_sgms1$celltype <- factor(
  as.character(expr_sgms1$celltype),
  levels = celltype_order
)


expr_sgms1_top25 <- expr_sgms1 |>
  group_by(celltype) |>
  arrange(
    desc(Sgms1),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()




target_group <- "Effector-like CD8aa+ TCRab/gd T cells"

comparisons <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd T cells"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)

p_sgms1 <- ggplot(
  expr_sgms1_top25,
  aes(
    x = celltype,
    y = Sgms1,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Sgms1",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none"
  )

print(p_sgms1)

save_plot(
  p_sgms1,
  "Sgms1_top25_Effector_CD8aa_comparison",
  width = 10,
  height = 7
)

expr_osbp <- FetchData(
  objWT1,
  vars = c(
    "Osbp",
    "celltype",
    "sample"
  ),
  layer = "data"
)

expr_osbp$celltype <- factor(
  as.character(expr_osbp$celltype),
  levels = celltype_order
)

expr_osbp_top25 <- expr_osbp |>
  group_by(celltype) |>
  arrange(
    desc(Osbp),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()

target_group <- "Effector-like CD8aa+ TCRab/gd T cells"

comparisons_osbp <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd T cells"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)

p_osbp <- ggplot(
  expr_osbp_top25,
  aes(
    x = celltype,
    y = Osbp,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = comparisons_osbp,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Osbp",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none"
  )

print(p_osbp)

save_plot(
  p_osbp,
  "Osbp_top25_Effector_CD8aa_comparison",
  width = 10,
  height = 7
)


expr_insig1 <- FetchData(
  objWT1,
  vars = c(
    "Insig1",
    "celltype",
    "sample"
  ),
  layer = "data"
)

expr_insig1$celltype <- factor(
  as.character(expr_insig1$celltype),
  levels = celltype_order
)

expr_insig1_top25 <- expr_insig1 |>
  group_by(celltype) |>
  arrange(
    desc(Insig1),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()

target_group <- "Effector-like CD8aa+ TCRab/gd T cells"

comparisons_insig1 <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd T cells"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)

p_insig1 <- ggplot(
  expr_insig1_top25,
  aes(
    x = celltype,
    y = Insig1,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = comparisons_insig1,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Insig1",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none"
  )

print(p_insig1)

save_plot(
  p_insig1,
  "Insig1_top25_Effector_CD8aa_comparison",
  width = 10,
  height = 7
)

expr_npc1 <- FetchData(
  objWT1,
  vars = c(
    "Npc1",
    "celltype",
    "sample"
  ),
  layer = "data"
)

expr_npc1$celltype <- factor(
  as.character(expr_npc1$celltype),
  levels = celltype_order
)

expr_npc1_top25 <- expr_npc1 |>
  group_by(celltype) |>
  arrange(
    desc(Npc1),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()

target_group <- "Effector-like CD8aa+ TCRab/gd T cells"

comparisons_npc1 <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd T cells"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)

p_npc1 <- ggplot(
  expr_npc1_top25,
  aes(
    x = celltype,
    y = Npc1,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = comparisons_npc1,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Npc1",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none"
  )

print(p_npc1)

save_plot(
  p_npc1,
  "Npc1_top25_Effector_CD8aa_comparison",
  width = 10,
  height = 7
)


sm_genes <- c(
  "Sgms1", "Sgms2", "Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2", "Degs1", "Kdsr"
)

objWT1<- AddModuleScore(
  objWT1, features = list(sm_genes), name = "SM_synthesis",
  seed = 20260729
)
sm_table <- objWT1[[]] |>
  transmute(celltype = celltype, SM_synthesis_score = SM_synthesis1)
save_table(sm_table, "WT_SM_synthesis_module_scores.csv")

p_sm_score <- ggplot(
  sm_table, aes(celltype, SM_synthesis_score, fill = celltype)
) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.65) +
  stat_compare_means(
    comparisons = comparisons, method = "wilcox.test", size = 3
  ) +
  scale_fill_manual(
    values = celltype_colors) +
  labs(title = "SM synthesis module score", x = NULL, y = "Module score") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
print(p_sm_score)
save_plot(p_sm_score, "16_WT_SM_synthesis_module_score", 7, 5.5)



chol_acquisition_genes <- c(
  "Srebf2",
  "Acat2",
  "Hmgcs1",
  "Hmgcr",
  "Sqle",
  "Mbtps1",
  "Mbtps2",
  "Lss",
  "Fdps",
  "Fdft1",
  "Mvd",
  "Mvk",
  "Pmvk",
  "Ggps1",
  "Cyp51",
  "Tm7sf2",
  "Msmo1",
  "Nsdhl",
  "Hsd17b7",
  "Ebp",
  "Sc5d",
  "Dhcr7",
  "Dhcr24",
  "Ldlr",
  "Npc1",
  "Npc2"
)


objWT1 <- AddModuleScore(
  objWT1,
  features = list(chol_acquisition_genes),
  name = "Chol_acquisition",
  seed = 20260729
)

# 提取评分及相关细胞信息
chol_table <- objWT1[[]] |>
  transmute(
    celltype = factor(
      as.character(celltype),
      levels = celltype_order
    ),
    sample = sample,
    Chol_acquisition_score = Chol_acquisition1
  )

# 保存每个细胞的模块评分
save_table(
  chol_table,
  "WT_cholesterol_acquisition_module_scores.csv"
)

# 设置目标细胞群
target_group <- "Effector-like CD8aa+ TCRab/gd T cells"


celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd T cells",
  "Effector-like CD8aa+ TCRab/gd T cells",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T Cells",
  "Treg cells",
  "Cycling T Cells"
)

objWT1$celltype <- factor(
  as.character(objWT1$celltype),
  levels = celltype_order
)
# 设置统计学比较
chol_comparisons <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd T cells"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)

# 绘制模块评分
p_chol_score <- ggplot(
  chol_table,
  aes(
    x = celltype,
    y = Chol_acquisition_score,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.65
  ) +
  stat_compare_means(
    comparisons = chol_comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 3,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Cholesterol acquisition module score",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_chol_score)

save_plot(
  p_chol_score,
  "WT_cholesterol_acquisition_module_score",
  width = 7,
  height = 5.5
)



################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd T cells",
  "Effector-like CD8aa+ TCRab/gd T cells",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T Cells",
  "Treg cells",
  "Cycling T Cells"
)

marker_genes <- c(
  "Trac",
  "Tcrg-C1",
  "Cd8a",
  "Cd4",
  "Cd8b1",
  "Tcf7",
  "Gzmb",
  "S1pr1",
  "Mki67",
  "Foxp3"
)

p_feature_umap <- FeaturePlot(
  object = iel,
  features = marker_genes,
  reduction = "umap",
  
  cols = c(
    "#D3D3D3",
    "#004C99"
  ),
  
  # 强化阳性细胞颜色
  min.cutoff = 0,
  max.cutoff = "q50",
  
  order = TRUE,
  pt.size = 0.6,
  
  ncol = 4,
  raster = TRUE,
  raster.dpi = c(300, 300),
  combine = TRUE
)
p_feature_umap <- p_feature_umap &
  theme_publication(
    base_size = 9,
    base_family = "sans"
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(
      face = "italic",
      size = 10,
      hjust = 0.5,
      colour = "black"
    ),
    
    legend.position = "none"
  )
p_feature_umap

save_plot(
  plot = p_feature_umap,
  filename = "UMAP_T_cell_marker_FeaturePlot",
  width = 10,
  height = 7.5,
  dpi = 600,
  bg = "white"
)
gene_colors <- c(
  "Trac"    = "#4C78A8",  # 蓝色
  "Tcrg-C1" = "#72A98F",  # 灰绿色
  "Cd8a"    = "#D98C95",  # 粉色
  "Cd4"     = "#E9C95B",  # 黄色
  "Cd8b1"   = "#C9365B",  # 玫红
  "Tcf7"    = "#F27658",  # 珊瑚橙
  "Gzmb"    = "#70B7C4",  # 青蓝
  "S1pr1"   = "#B6A6A0",  # 暖灰
  "Mki67"   = "#E99662",  # 浅橙
  "Il17a"   = "#8E6BBE",  # 紫色
  "Foxp3"   = "#58A65C"   # 绿色
)

stacked_violin_plot <- function(
    object,
    features,
    celltype_col = "celltype",
    celltype_order,
    assay = "RNA"
) {
  DefaultAssay(object) <- assay
  
  # 检查基因是否存在
  missing_features <- setdiff(features, rownames(object))
  
  if (length(missing_features) > 0) {
    stop(
      "以下基因未在当前 assay 中找到：",
      paste(missing_features, collapse = ", ")
    )
  }
  
  plot_data <- FetchData(
    object = object,
    vars = c(celltype_col, features)
  )
  
  plot_data$celltype <- factor(
    plot_data[[celltype_col]],
    levels = celltype_order
  )
  
  plot_data <- plot_data |>
    filter(!is.na(celltype)) |>
    select(celltype, all_of(features)) |>
    pivot_longer(
      cols = all_of(features),
      names_to = "feature",
      values_to = "expression"
    )
  
  # 控制基因从上至下的顺序
  plot_data$feature <- factor(
    plot_data$feature,
    levels = features
  )
  
  ggplot(
    plot_data,
    aes(
      x = celltype,
      y = expression,
      fill = feature,
      group = interaction(feature, celltype)
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "#D8D8D8",
      linewidth = 0.4
    ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      width = 0.86,
      adjust = 1,
      colour = "#555555",
      linewidth = 0.3
    ) +
    facet_grid(
      rows = vars(feature),
      scales = "free_y",
      switch = "y"
    ) +
    scale_x_discrete(
      limits = celltype_order,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = gene_colors,
      drop = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # 细胞群对应的纵向参考线
      panel.grid.major.x = element_line(
        colour = "#E5E5E5",
        linewidth = 0.4
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      
      axis.text.x = element_text(
        colour = "#333333",
        angle = 48,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
      axis.text.y = element_blank(),
      
      # 在左侧显示基因名称
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        colour = "#333333",
        angle = 0,
        hjust = 1,
        face = "italic",
        size = 9
      ),
      
      panel.spacing.y = grid::unit(0, "pt"),
      legend.position = "none",
      
      plot.margin = margin(
        t = 15,
        r = 15,
        b = 15,
        l = 15
      )
    )
}

p_stacked_violin <- stacked_violin_plot(
  object = iel,
  features = marker_genes,
  celltype_col = "celltype",
  celltype_order = celltype_order,
  assay = "RNA"
)
p_stacked_violin 
save_plot(
  plot = p_stacked_violin,
  filename = "stacked_violin_T_cell_markers",
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)
###############################################################################
# SGMS1 KO IEL single-cell RNA-seq analysis
#这部分是wt vs ko
###############################################################################


library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(data.table)


iwt<-Read10X(data.dir ="/single cell/Result-X101SC24094483-Z01-J001-B1-1.20241024/Result-X101SC24094483-Z01-J001-B1-1/2.Summary/IWT/filtered_feature_bc_matrix")
iko<-Read10X(data.dir ="/Result-X101SC24094483-Z01-J001-B1-1.20241024/Result-X101SC24094483-Z01-J001-B1-1/2.Summary/IKO/filtered_feature_bc_matrix")
siwt <- CreateSeuratObject(counts = iwt, project = "iwt")
siko <- CreateSeuratObject(counts = iko, project = "iko")
objBRC <- merge(siwt,siko,
                add.cell.ids = c("iwt","iko"),
                project = "total") 

objBRC[["percent.MT"]] <- PercentageFeatureSet(objBRC, pattern = "^Mt")
objBRC_bp <- objBRC
objBRC2=objBRC_bp
objBRC2 <- subset(objBRC2,subset = nFeature_RNA > 200 & nFeature_RNA <6000 & nCount_RNA<20000 & percent.MT <5)
VlnPlot(objBRC2, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
objBRC2 <- NormalizeData(objBRC2) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
objBRC2 <- RunPCA(objBRC2, features = VariableFeatures(object = objBRC2))


objBRC2=IntegrateLayers(
  object = objBRC2, method = HarmonyIntegration, 
  orig.reduction ="pca", new.reduction = "harmony",
  verbose = TRUE
)

objBRC2 <- RunPCA(objBRC2, features = VariableFeatures(object = objBRC2))
ElbowPlot(objBRC2, ndims = ncol(Embeddings(objBRC2, "pca")))#13

DimHeatmap(objBRC, dims = c(1:30), fast = T)

objBRC2<- FindNeighbors(objBRC2, dims = 1:20)

objBRC2<- FindClusters(objBRC2, resolution = 0.3)
objBRC2 <- RunUMAP(objBRC2, dim = 1:20)
DimPlot(objBRC, reduction = 'umap',label=T)
DimPlot(objBRC, reduction = 'umap', group.by ="orig.ident")
DimPlot(objFIN1, split.by = 'orig.ident')

Idents(objBRC)="seurat_clusters"

DotPlot(objBRC2,features = c("Trac","Tcrg-V7","Tcrg-C1","Trbv16","Cd8a","Cd8b1","Cd4","Foxp3","Il17a","Gzmb","Tcf7","Itgae","Tyrobp","Ccr7","S1pr1","Lef1","Rora","Cd5","Cd7","Gzma","Stmn1","Hist1h2ae","Glp1r","Cd19"))+coord_flip()+RotatedAxis()

meta <- objBRC@meta.data
objSUB2 = subset(x = objBRC2, idents = 3)
ElbowPlot(objSUB2, ndims = ncol(Embeddings(objSUB2, "pca")))#20
objSUB2<- FindNeighbors(objSUB2, dims = 1:20)

objSUB2<- FindClusters(objSUB2, resolution = 1)
objSUB2 <- RunUMAP(objSUB2, dim = 1:20)
DotPlot(objSUB2,features = c("Trac","Tcrg-C2","Cd8a","Cd8b1","Cd4","Foxp3","Il17a","Gzmb","Tcf7","Itgae","Tyrobp","Ccr7","S1pr1","Lef1","Rora","Cd5","Cd7","Gzma"))+coord_flip()+RotatedAxis()
DimPlot(objSUB2, split.by = 'orig.ident', label = T)
VlnPlot(objSUB2, features = c("Trac","Tcrg-V1","Cd8a","Cd8b1","Cd4","Tcf7","Gzmb"),pt.size = 0, ncol = 2)

metasub <- objSUB2@meta.data
metasub$celltype <- 'null'
metasub[which(metasub$seurat_clusters %in% c(2)),9] <- 'CD8aa+ CD4+ TCRab cells'
metasub[which(metasub$seurat_clusters %in% c(5)),9] <- 'Th17 cells'
metasub[which(metasub$seurat_clusters %in% c(0,1,3,4,6,7)),9] <- 'CD4 T cells'

objSUB2@meta.data <- metasub
Idents(objSUB2) <- 'celltype'
DimPlot(objSUB2, reduction = 'umap', label = T)
metamain <- objBRC2@meta.data
metamain$celltype <- 'null'
metamain[which(metamain$seurat_clusters %in% c(1)),7] <- 'CD8+ T cells'
metamain[which(metamain$seurat_clusters %in% c(0)),7] <- 'Memory CD8aa+ TCRab/gd T cells'
metamain[which(metamain$seurat_clusters %in% c(2,4,7)),7] <- 'Effector CD8aa+ TCRab/gd T cells'
metamain[which(metamain$seurat_clusters %in% c(5)),7] <- 'Naive T cells'
metamain[which(metamain$seurat_clusters %in% c(9)),7] <- 'Treg cells'
metamain[which(metamain$seurat_clusters %in% c(10)),7] <- 'ISG+ CD8 T Cells'
metamain[which(metamain$seurat_clusters %in% c(8)),7] <- 'Cycling T Cells'
metamain[which(metamain$seurat_clusters %in% c(6)),7] <- 'Memory CD4 T Cells'
objBRC2@meta.data <- metamain
Idents(objBRC2) <- 'celltype'
DimPlot(objBRC2, reduction = 'umap', label = T)

Idents(objSUB2)="Celltype_sub"
Idents(objBRC2)="Celltype_major"

objBRC2@meta.data$cellid <- rownames(objBRC2@meta.data)
objSUB2@meta.data$cellid <- rownames(objSUB2@meta.data)

metadata<- left_join(x=objBRC2@meta.data,y= objSUB2@meta.data, by= "cellid")
metadata$merge_celltype <- ifelse(is.na(metadata$celltype.y), metadata$celltype.x,metadata$celltype.y)
objBRC2@meta.data$merge_celltype <- metadata$merge_celltype
Idents(objBRC2) <- "merge_celltype"
Idents(objBRC2) <- "seurat_clusters"

objFIN = subset(x = objBRC2, idents = c(0,1,2,3,4,5,6,7,8,9,10))
Idents(objFIN) <- "merge_celltype"
DimPlot(objFIN, split.by = 'orig.ident')
objBRC<-objFIN 

meta <- objBRC@meta.data
meta$newcelltype <- 'null'

meta[which(meta$merge_celltype %in% c('Memory CD8aa+ TCRab/gd T cells','Effector CD8aa+ TCRab/gd T cells')),10] <- 'CD8aa+ TCRab/gd T cells'
meta[which(meta$merge_celltype %in% c('Memory CD4 T Cells','CD4 T cells','CD8aa+ CD4+ TCRab cells','Th17 cells','Treg cells')),10] <- 'CD4 T cells'
meta[which(meta$merge_celltype %in% c('Naive T cells')),10] <- 'Naive T cells'
meta[which(meta$merge_celltype %in% c("ISG+ CD8 T Cells",'CD8+ T cells')),10] <- 'CD8 T cells'
meta[which(meta$merge_celltype %in% c('Cycling T Cells')),10] <- 'Cycling T cells'

library(tidydr)

objBRC@meta.data <- meta
Idents(objBRC) <- "merge_celltype"

DimPlot(objBRC,cols = paried,split.by = "orig.ident",pt.size = 0.2)+theme(panel.grid=element_blank())
DimPlot(objBRC)

table(Idents(objBRC), objBRC$orig.ident)


objSUBN = subset(x = objBRC, idents = c("CD8aa+ TCRab/gd T cells"))

ElbowPlot(objSUBN, ndims = ncol(Embeddings(objSUBN, "pca")))
objSUBN<- FindNeighbors(objSUBN, dims = 1:13)

objSUBN<- FindClusters(objSUBN, resolution = 1)
objSUBN <- RunUMAP(objSUBN, dim = 1:13)

DimPlot(objSUBN, split.by = 'orig.ident')
VlnPlot(objSUBN,features = c("Trac","Tcrg-C2","Cd8a","Ccl4","Foxp3","Il17a","Gzmb","Tcf7","Prdm1","Lag3"))+coord_flip()+RotatedAxis()

meta1 <- objSUBN@meta.data
meta1$celltype <- 'null'

meta1[which(meta1$seurat_clusters %in% c(1,2,7,8,6)),7] <- 'Gzmlow CD8aa+ TCRab/gd T cells'
meta1[which(meta1$seurat_clusters %in% c(0,3,5,4,9)),7] <- 'Gzmb+ CD8aa+ TCRab/gd T cells'
objSUBN@meta.data <- meta1
Idents(objSUBN) <- "celltype"
DimPlot(objSUBN, reduction = 'umap')
FeaturePlot(objSUBN,feature="Gzmb")


DimPlot(objSUBN,cols = paried,pt.size = 0.1)+theme(panel.grid=element_blank())
FeaturePlot(objSUBN,features = c("Gzmb"),pt.size = 0.5,cols = c("lightgrey","#4f86c6"))
FeaturePlot(objSUBN,features = c("Gzma"),pt.size = 0.5,cols = c("lightgrey","#4f86c6"))


table(Idents(objSUBN), objSUBN$orig.ident)

Cellratio <- prop.table(table(Idents(objSUBN), objSUBN$orig.ident), margin = 2)

Cellratio
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.grid =element_blank(),panel.border = element_blank(),axis.line = element_blank(),legend.key = element_blank(),legend.background = element_blank())+
  scale_fill_manual(values=paried)


objSUBN$cluster_group<-paste(objSUBN$celltype,objSUBN$orig.ident)
Idents(objSUBN)<-"cluster_group"
DimPlot(objSUBN)
ENIEL <- FindMarkers(objSUBN, ident.1 = "Gzmb+ CD8aa+ TCRab/gd T cells iko", ident.2 = "Gzmb+ CD8aa+ TCRab/gd T cells iwt", verbose = FALSE)
data<-ENIEL
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0.58 #设置log2FC的阈值
data$Sig = ifelse(data$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
data = data.frame(data)
write.csv(data, file = "eniel2.csv", row.names = FALSE)
table(data$Sig) #查看数据统计情况


data$genename<-row.names(data)
data$label<-NA
data$label[which(data$genename=="Gzma")]<-"Gzma"
data$label[which(data$genename=="Ccl4")]<-"Ccl4"
data$label[which(data$genename=="Tnfrsf9")]<-"Tnfrsf9"
data$label[which(data$genename=="Itgad")]<-"Itgad"
data$label[which(data$genename=="Lag3")]<-"Lag3"
data$label[which(data$genename=="Cdk8")]<-"Cdk8"
data$label[which(data$genename=="Il34")]<-"Il34"
data$label[which(data$genename=="Ptpn5")]<-"Ptpn5"

data$label[which(data$genename=="Lef1")]<-"Lef1"
data$label[which(data$genename=="Id3")]<-"Id3"

data$label[which(data$genename=="Jun")]<-"Jun"
data$label[which(data$genename=="Nfkbia")]<-"Nfkbia"
data$label[which(data$genename=="Xcl1")]<-"Xcl1"
data$label[which(data$genename=="Bcl3")]<-"Bcl3"
data$label[which(data$genename=="Cd200r2")]<-"Cd200r2"
data$label[which(data$genename=="Klra1")]<-"Klra1"
data$label[which(data$genename=="Klra3")]<-"Klra3"
data$label[which(data$genename=="Klra7")]<-"Klra7"
data$label[which(data$genename=="Havcr2")]<-"Havcr2"
data$label[which(data$genename=="Snx9")]<-"Snx9"

library(ggplot2)
library(ggrepel)


###绘图——基础火山图###
p1 <- ggplot(data, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=1, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#4f86c6", "#d2dae2","#d64065")) + xlim(c(-1.5, 1.5))  + #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.7) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.7) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  theme_bw()+ylim(0,65)+theme(
    panel.border = element_blank(),  panel.background = element_blank(),panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

library(ggrepel)
p1
###添加基因名标记###
p1+geom_label_repel(data = data, aes(label = label),
                    size = 4,fontface="italic",                      # 设置标签
                    show.legend = FALSE,                # 不显示图例
                    max.overlaps = Inf, min.segment.length=0,point.padding=0)               # 设置标签重叠的最大次数


###################################################################
Idents(objBRC) <- "newcelltype"
objSUB48 = subset(x = objBRC, idents = c("CD4 T cells",'CD8 T cells'))
ElbowPlot(objSUB48, ndims = ncol(Embeddings(objSUB4, "pca")))
objSUB48<- FindNeighbors(objSUB4, dims = 1:18)

objSUB48<- FindClusters(objSUB4, resolution = 1)
objSUB48 <- RunUMAP(objSUB4, dim = 1:18)
DotPlot(objSUB48,features = c("Trac","Trbc1","Cd4","Cd8b1","Cd8a","Il17a","Foxp3","Lag3",
                              "Gzmb","Gzma","Il2rb","Ifng","Klrk1","Hcst",
                              "Tcf7","Ccr7","Lef1","S1pr1","Il7r"))+coord_flip()+RotatedAxis()

FeaturePlot(objSUB48,features = c("Gzmb","Gzma","Ccl4","Ifng","Tnf","Tbx21","Eomes","Tcf7","Il7r","Ccr7","Sell"))

meta2 <- objSUB48@meta.data
meta2$celltype <- 'null'
Idents(objSUB48) <- "seurat_clusters"

meta2[which(meta2$seurat_clusters %in% c(2)),7] <- 'Tcf1+ CD8 T cells'
meta2[which(meta2$seurat_clusters %in% c(0,3,6,9,12)),7] <- 'Gzmb+ CD8 T cells'
meta2[which(meta2$seurat_clusters %in% c(4,5,11,13,10)),7] <- 'Tcf1+ CD4 T cells'
meta2[which(meta2$seurat_clusters %in% c(1)),7] <- 'Gzmb+ CD4 T cells'
meta2[which(meta2$seurat_clusters %in% c(7)),7] <- 'Th17 cells'
meta2[which(meta2$seurat_clusters %in% c(8)),7] <- 'Treg cells'
objSUB48@meta.data <- meta2
Idents(objSUB48) <- "celltype"
DimPlot(objSUB48, reduction = 'umap')

paried=c("#f8845c","#5eb887","#f5d983","#d88c8f","#4f86c6","#d64065")
paried=c("#FADAB7","#d88c8f")

DimPlot(objSUB48,cols = paried,pt.size = 0.1)+theme(panel.grid=element_blank())
FeaturePlot(objSUB48,features = c("Gzmb"),cols = c("lightgrey","#4f86c6"))
FeaturePlot(objSUB48,features = c("Cd8b1"),cols = c("lightgrey","#4f86c6"))
FeaturePlot(objSUB48,features = c("Tcf7"),cols = c("lightgrey","#4f86c6"))
FeaturePlot(objSUB48,features = c("Gzma"),cols = c("lightgrey","#4f86c6"))

FeaturePlot(objSUB48,features = c("Il17a"),cols = c("lightgrey","#4f86c6"))
FeaturePlot(objSUB48,features = c("Foxp3"),cols = c("lightgrey","#4f86c6"))
table(Idents(objSUB48), objSUB48$orig.ident)

Cellratio <- prop.table(table(Idents(objSUB48), objSUB48$orig.ident), margin = 2)

Cellratio
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.grid =element_blank(),panel.border = element_blank(),axis.line = element_blank(),legend.key = element_blank(),legend.background = element_blank())


# 计算每个样本中各细胞类型的比例
ratio <- objSUBN@meta.data %>%
  # 按样本(orig.ident)和细胞类型(celltype_new_3)分组
  group_by(orig.ident, celltype) %>%
  # 计算每组细胞数量
  summarise(n = n()) %>%
  # 计算相对频率(当前组细胞数/同一样本总细胞数)
  mutate(relative_freq = n / sum(n))

ratio$celltype <- factor(
  ratio$celltype,
  levels = c(
    'Gzmb+ CD8 T cells','Tcf1+ CD8 T cells','Gzmb+ CD4 T cells', 'Tcf1+ CD4 T cells', 'Th17 cells',
    'Treg cells'
  )
)




# 设置样本的因子水平(控制绘图顺序)
ratio$orig.ident <- factor(ratio$orig.ident, levels = c('iwt', 'iko'))

# 创建组合图(堆叠柱状图 + 冲积图)
ggplot(ratio, 
       aes(x = orig.ident,          # X轴：样本分组
           y = relative_freq,        # Y轴：相对频率
           fill = celltype,    # 填充颜色：细胞类型
           label = celltype,   # 标签文本：细胞类型
           stratum = celltype, # 冲积图分层：细胞类型
           alluvium = celltype # 冲积图连接：细胞类型
       )) +
  # 绘制堆叠柱状图
  geom_col(width = 0.7) +
  # 添加流动曲线(展示比例变化)
  
  # 添加冲积图层(连接不同样本间的相同细胞类型)
  scale_fill_manual(values = c( "#d88c8f","#5eb887","#f5d983",  "#4f86c6" ,"#f8845c", "#d64065"))+
  # 设置坐标轴和图例标签
  labs(x = 'Sample', y = 'Ratio', fill = 'Celltype') +
  # 使用经典主题(无背景网格)
  theme_classic() +
  # 自定义主题细节
  theme(
    panel.grid =element_blank(),panel.border = element_blank(),axis.line = element_blank(),legend.key = element_blank(),legend.background = element_blank())+ 
  geom_flow(
    curve_type = 'sine',  # 使用正弦曲线
    alpha = 0.001,          # 设置透明度
    knot.pos = 0.5        # 节点位置(0.5=中间)
  ) +
  # 添加冲积图层(连接不同样本间的相同细胞类型)
  geom_alluvium(width = 0.6) +theme_classic()



obj$cluster_group<-paste(objSUBN$celltype,objSUBN$orig.ident)
Idents(objBRC)<-"cluster_group"
DimPlot(objBRC)
CD4 <- FindMarkers(objBRC, ident.1 = "CD4 T cells iko", ident.2 = "CD4 T cells iwt", verbose = FALSE)
CD8 <- FindMarkers(objBRC, ident.1 = "CD8 T cells iko", ident.2 = "CD8 T cells iwt", verbose = FALSE)
data<-CD8
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0.58 #设置log2FC的阈值
data$Sig = ifelse(data$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
data = data.frame(data)
data$genename<-row.names(data)
write.csv(data, file = "CD8.csv", row.names = FALSE)
table(data$Sig) #查看数据统计情况


data$genename<-row.names(data)
data$label<-NA
data$label[which(data$genename=="Gzma")]<-"Gzma"
data$label[which(data$genename=="Ccl4")]<-"Ccl4"
data$label[which(data$genename=="Cxcr6")]<-"Cxcr6"
data$label[which(data$genename=="Itgae")]<-"Itgae"
data$label[which(data$genename=="Cd38")]<-"Cd38"
data$label[which(data$genename=="Ccr5")]<-"Ccr5"
data$label[which(data$genename=="Cd244a")]<-"Cd244a"

library(ggplot2)
library(ggrepel)


###绘图——基础火山图###
p1 <- ggplot(data, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=1, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#4f86c6", "#d2dae2","#d64065")) + xlim(c(-1.5, 1.5))  + #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.7) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.7) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  theme_bw()+ylim(0,65)+theme(
    panel.border = element_blank(),  panel.background = element_blank(),panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
p1+geom_label_repel(data = data, aes(label = label),
                    size = 4,fontface="italic",                      # 设置标签
                    show.legend = FALSE,                # 不显示图例
                    max.overlaps = Inf, min.segment.length=0,point.padding=0)               # 设置标签重叠的最大次数
FeaturePlot(objSUB48,"Slamf7")


data<-CD8
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0.58 #设置log2FC的阈值
data$Sig = ifelse(data$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                  ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
data = data.frame(data)
data$genename<-row.names(data)
write.csv(data, file = "CD8.csv", row.names = FALSE)
table(data$Sig) #查看数据统计情况


data$genename<-row.names(data)
data$label<-NA
data$label[which(data$genename=="Lat2")]<-"Lat2"
data$label[which(data$genename=="Ccl4")]<-"Ccl4"
data$label[which(data$genename=="Lag3")]<-"Lag3"
data$label[which(data$genename=="Lgals3")]<-"Lgals3"
data$label[which(data$genename=="Ly6a")]<-"Ly6a"
data$label[which(data$genename=="Cd200r2")]<-"Cd200r2"
data$label[which(data$genename=="Tnfrsf9")]<-"Tnfrsf9"
data$label[which(data$genename=="Ifit1")]<-"Ifit1"
data$label[which(data$genename=="Ifit3")]<-"Ifit3"
data$label[which(data$genename=="Ifi213")]<-"Ifi213"

library(ggplot2)
library(ggrepel)


###绘图——基础火山图###
p1 <- ggplot(data, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=1, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#4f86c6", "#d2dae2","#d64065")) + xlim(c(-1.5, 1.5))  + #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.7) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.7) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  theme_bw()+ylim(0,65)+theme(
    panel.border = element_blank(),  panel.background = element_blank(),panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
p1+geom_label_repel(data = data, aes(label = label),
                    size = 4,fontface="italic",                      # 设置标签
                    show.legend = FALSE,                # 不显示图例
                    max.overlaps = Inf, min.segment.length=0,point.padding=0,force = 60)               # 设置标签重叠的最大次数


#############################################################################GO was done online by DAVID
library(org.Mm.eg.db)
library(biomaRt)

# 连接小鼠基因组数据库
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 基因名 → Ensembl ID
symbols <- CD8$genename
result <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                filters = "external_gene_name", 
                values = symbols, 
                mart = mart)

colnames(result) <- c("genename", "ensembl_id")
CD8<- merge(CD8, result, by = "genename", all.x = TRUE)

nielid$Sig = ifelse(nielid$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                      abs(nielid$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                    ifelse(nielid$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')

write.csv(CD8, file = "CD8id.csv", row.names = FALSE)
goniel<-read.delim("niel go.txt")


colnames(result) <- c("symbol", "ensembl_id")
nielid <- merge(data, result, by = "symbol", all.x = TRUE)

nielid$Sig = ifelse(nielid$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                      abs(nielid$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                    ifelse(nielid$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')

write.csv(nielid, file = "nielid.csv", row.names = FALSE)
data<-CD8
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0.58 #设置log2FC的阈值
data$Sig = ifelse(data$p_val_adj < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    abs(data$avg_log2FC) >= cut_off_log2FC)  #abs绝对值
goCD4<-read.delim("CD4 go.txt")
goCD8<-read.delim("CD8 go.txt")

goniel<- goniel %>% mutate(segment = "Natural")  # 添加来源标识

goCD4<- goCD4 %>% mutate(segment = "CD4") 
goCD8<- goCD8 %>% mutate(segment = "CD8")



go.N<-goniel[which(goniel$PValue<0.05),]
go.4<-goCD4[which(goCD4$PValue<0.05),]
go.8<-goCD8[which(goCD8$PValue<0.05),]



top5go_N <- head(go.N[order(go.N$PValue), ], 6)
top5go_4 <- head(go.4[order(go.4$PValue),], 6)
top5go_8 <- head(go.8[order(go.8$PValue),], 6)

top5go_N1 <- head(top5go_N[order(top5go_N$PValue,decreasing = T), ], 6)
top5go_41 <- head(top5go_4[order(top5go_4$PValue,decreasing = T), ], 6)
top5go_81 <- head(top5go_8[order(top5go_8$PValue,decreasing = T), ], 6)
go_up <- bind_rows(top10go_D, top10go_J, top10go_I, top10go_C) 


plotN<- ggplot(top5go_N1,aes(x=factor(segment, level=c("Natural")),y =factor(Term,levels=top5go_N1$Term) ))+
  geom_point(aes(size=Count,color=-log10(PValue)))+
  scale_color_gradient(low = "#6baedb",high = "#d64065")+theme_bw() +theme(axis.text.y=element_text(size=10, color = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size=10))+scale_size_continuous(range=c(3,9))
plotN

plot4<- ggplot(top5go_41,aes(x=factor(segment, level=c("CD4")),y =factor(Term,levels=top5go_41$Term)))+
  geom_point(aes(size=Count,color=-log10(PValue)))+
  scale_color_gradient(low = "#6baedb",high = "#d64065")+theme_bw() +theme(axis.text.y=element_text(size=10, color = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size=10))+scale_size_continuous(range=c(3,9))
plot4

plot8<- ggplot(top5go_81,aes(x=factor(segment, level=c("CD8")),y =factor(Term,levels=top5go_81$Term)))+
  geom_point(aes(size=Count,color=-log10(PValue)))+
  scale_color_gradient(low = "#6baedb",high = "#d64065")+theme_bw() +theme(axis.text.y=element_text(size=10, color = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size=10))+scale_size_continuous(range=c(3,9))
plot8
###############################################################################
iel vs lpl single cell
###############################################################################

set.seed(20260805)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(svglite)
  library(rstatix)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(ggsignif)
})

# ---- Project paths -----------------------------------------------------------
project_dir <- path.expand(
  "~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/iel lpl"
)
figure_dir <- file.path(project_dir, "results", "figures")
table_dir  <- file.path(project_dir, "results", "tables")
object_dir <- file.path(project_dir, "results", "objects")
log_dir    <- file.path(project_dir, "results", "logs")

invisible(lapply(
  c(project_dir, figure_dir, table_dir, object_dir, log_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))


# ---- Reusable publication settings ------------------------------------------
sample_colors <- c(iel = "#4F86C6", lpl = "#D64065")




celltype_colors <- c(
  # Unconventional T：绿色系
  "Memory-like unconventional T cells" =
    "#91D0B3",
  
  "Effector-like unconventional T cells" =
    "#3E9B72",
  
  # gdT17：黄绿色
  "gdT17 cells" =
    "#A8C95F",
  
  # NKT：蓝绿色
  "NKT cells" =
    "#35A6A6",
  
  # Naive/CD4 T：浅橙—橙红色系
  "Naive T cells" =
    "#F5D983",
  
  "CD4 T cells" =
    "#F6AA8D",
  
  "Th17 cells" =
    "#F07850",
  
  # Treg：紫色
  "Treg cells" =
    "#8E6BBE",
  
  "CD8ab+ T cells" =
    "#D64065",
  
  # Cycling：蓝色
  "Cycling T cells" =
    "#4F86C6"
)

theme_publication <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}
theme_set(theme_publication())

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}


theme_publication <- function(
    base_size = 9,
    base_family = "sans"
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}

theme_set(theme_publication())

save_plot <- function(
    plot,
    filename,
    width,
    height,
    units = "in",
    dpi = 600,
    bg = "white"
) {
  stem <- file.path(
    figure_dir,
    safe_name(filename)
  )
  
  # PDF矢量图
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = units,
    bg = bg,
    limitsize = FALSE
  )
  
  # SVG矢量图
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height,
    bg = bg
  )
  print(plot)
  grDevices::dev.off()
  
  # PNG位图
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  
  invisible(plot)
}

save_table <- function(x, filename, row.names = FALSE) {
  write.csv(
    x, file.path(table_dir, filename),
    row.names = row.names, fileEncoding = "UTF-8"
  )
}

# ---- Read, merge, QC ---------------------------------------------------------
iel_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/iel lpl/IIEL_wt/outs")
lpl_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/iel lpl/ILPL_wt/outs")

iel <- CreateSeuratObject(iel_counts, project = "iel")
lpl <- CreateSeuratObject(lpl_counts, project = "lpl")

iel <- merge(iel,lpl,
             add.cell.ids = c("iel","lpl"),
             project = "total") 

iel$sample <- factor(iel$orig.ident, levels = c("iel", "lpl"))

# Mouse mitochondrial genes usually begin "mt-"; ignore.case covers "Mt-".
iel[["percent.mt"]] <- PercentageFeatureSet(
  iel, pattern = "^mt-", assay = "RNA"
)

p_qc_before <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_before)
save_plot(p_qc_before, "01_QC_before_filtering", 10.5, 3.5)

iel <- subset(
  iel,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    nCount_RNA < 12000 &
    percent.mt < 5
)


saveRDS(iel, file.path(object_dir, "iellpl.rds"), compress = FALSE)

p_qc_after <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_after)
save_plot(p_qc_after, "02_QC_after_filtering", 10.5, 3.5)

qc_summary <- iel[[]] |>
  group_by(sample) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
save_table(qc_summary, "QC_summary_after_filtering.csv")

# ---- Normalization, Harmony, clustering -------------------------------------
iel <- NormalizeData(iel, verbose = FALSE)
iel <- FindVariableFeatures(iel, nfeatures = 2000, verbose = FALSE)
iel <- ScaleData(iel, features = VariableFeatures(iel), verbose = FALSE)
iel <- RunPCA(iel, features = VariableFeatures(iel), verbose = FALSE)

p_elbow <- ElbowPlot(iel, ndims = 50) + theme_publication()
print(p_elbow)
save_plot(p_elbow, "03_PCA_elbow", 5.5, 4)


dims_main <- 1:30
iel <- FindNeighbors(iel, dims = dims_main, verbose = FALSE)
iel <- FindClusters(iel, resolution = 1, random.seed = 20260805, verbose = FALSE)
iel <- RunUMAP(
  iel, dims = dims_main,
  seed.use = 20260805, verbose = FALSE
)



iel <- RunHarmony(
  iel, group.by.vars = "sample", reduction = "pca",
  reduction.save = "harmony", assay.use = "RNA", verbose = FALSE
)

iel <- FindNeighbors(iel, reduction = "harmony", dims = dims_main, verbose = FALSE)
iel <- FindClusters(iel, resolution = 1, random.seed = 20260805, verbose = FALSE)
iel <- RunUMAP(
  ielh, reduction = "harmony", dims = dims_main,
  seed.use = 20260805, verbose = FALSE
)



p_clusters <- DimPlot(
  iel, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_clusters)
save_plot(p_clusters, "04_UMAP_clusters", 6, 5)

p_samples <- DimPlot(
  iel, reduction = "umap", group.by = "sample",
  cols = sample_colors, raster = FALSE
) + theme_publication()
print(p_samples)
save_plot(p_samples, "05_UMAP_samples", 6, 5)

p_split <- DimPlot(
  iel, reduction = "umap", group.by = "seurat_clusters",
  split.by = "sample", label = TRUE, repel = TRUE, raster = FALSE
) + theme_publication()
print(p_split)
save_plot(p_split, "06_UMAP_clusters_split_by_sample", 10, 5)
#############################################################

saveRDS(iel, file.path(object_dir, "iel lpl_unannotated.rds"), compress = FALSE)
# =============================================================================
# 计算每个 cluster 的特征基因
# =============================================================================

# 建立单独对象，避免改变主对象iel
iel_markers <- iel

DefaultAssay(iel_markers) <- "RNA"

# Seurat v5中，差异分析前合并RNA assay的不同layers
if (inherits(iel_markers[["RNA"]], "Assay5")) {
  iel_markers <- JoinLayers(
    iel_markers,
    assay = "RNA"
  )
}

# 将细胞身份设置为cluster
Idents(iel_markers) <- "seurat_clusters"

# 计算每个cluster相对于其他全部细胞的上调特征基因
cluster_markers <- FindAllMarkers(
  object = iel_markers,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# 判断Seurat版本使用的fold-change列名
fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cluster_markers)
)[1]

if (is.na(fc_column)) {
  stop("FindAllMarkers结果中没有找到fold-change列。")
}

# 按cluster、显著性和fold change排序
cluster_markers <- cluster_markers |>
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]])
  )

# 保存所有cluster的全部marker
save_table(
  cluster_markers,
  "all_cluster_marker_genes.csv"
)


top15_cluster_markers <- cluster_markers |>
  filter(
    p_val_adj < 0.05
  ) |>
  group_by(cluster) |>
  slice_max(
    order_by = .data[[fc_column]],
    n = 15,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(
    cluster,
    desc(.data[[fc_column]])
  )

save_table(
  top15_cluster_markers,
  "top15_marker_genes_per_cluster.csv"
)

#==============================================================

marker_panel <- c(
  "Trac", "Trdc","Cd3d", "Trgc1", "Trgc2","Cd8a", "Cd4", "Cd8b1", "Klrb1b",
  "Nkg7", "Klrk1", "Klrb1c", "Ncr1", "Zbtb16", "Klrd1", "Il17a", "Il1r1",
  "Ccr6", "Il23r", "Foxp3", "Il10", "Sell", "Ccr7", "Cd44", "Mki67", "Stmn1",
  "Top2a", "Tbx21", "Ifng", "Gata3", "Il4", "Il5", "Il13", "Tcf7", "Ms4a4b",
  "Il7r", "Xcl1", "Cd244a", "Klf3", "Ccl4", "Cd200r1", "Themis",
  "Cd6", "Gzmk", "Ccl5", "Gzmb", "Gzma", "Cd40lg", "Maf", "Ctla4", "Icos",
  "Tyrobp", "Fcer1g", "Cd200r2", "Cd160", "Kit", "Id3",
  "Pclaf", "S1pr1", "Lef1", "Cd19", "C1qa",
  "Lyz2"
)

marker_panel <- intersect(marker_panel, rownames(iel))
p_markers <- DotPlot(iel, features = marker_panel) +
  coord_flip() + RotatedAxis() + theme_publication()
print(p_markers)
save_plot(p_markers, "07_cluster_marker_dotplot", 10, 8)

#260807 rds
saveRDS(iel,file = "260807 iel lpldimplot.rds")

#0805 end
# ---- Major cell-type annotation ---------------------------------------------
# IMPORTANT: verify these cluster-to-cell-type mappings from marker plots before
# biological interpretation. Named assignment avoids fragile numeric column use.
major_map <- c(
  "0"  = "Effector-like unconventional T cells",
  "1"  = "CD4 T cells",
  "2"  = "Effector-like unconventional T cells",
  "3"  = "Effector-like unconventional T cells",
  "4"  = "CD8ab+ T cells",
  "5"  = "Memory-like unconventional T cells",
  "6"  = "CD8ab+ T cells",
  "7"  = "Memory-like unconventional T cells",
  "8"  = "Naive T cells",
  "9"  = "Effector-like unconventional T cells",
  "10" = "Treg cells",
  "11" = "Th17 cells",
  "12" = "Cycling T cells",
  "13" = "CD8ab+ T cells",
  "14" = "NKT cells",
  "15" = "gdT17 cells",
  "16" = "CD8ab+ T cells",
  "17" = "NKT cells",
  "18" = "MHC-II+ myeloid/APC",
  "19" = "Effector-like unconventional T cells",
  "20" = "Plasma cells"
)
iel$celltype <- unname(major_map[as.character(iel$seurat_clusters)])

meta<-iel@meta.data
saveRDS(iel,file = "260810 sgms1 iel lpl dimplot.rds")
iel<-iel_b
iel <- subset(
  iel,
  subset = celltype != "MHC-II+ myeloid/APC"
)

iel <- subset(
  iel,
  subset = celltype != "Plasma cells"
)

Idents(iel) <- "celltype"
celltype_order <- c(
  "Memory-like unconventional T cells",
  "Effector-like unconventional T cells",
  "gdT17 cells",
  "NKT cells",
  "Naive T cells",
  "CD4 T cells",
  "Th17 cells",
  "Treg cells",
  "CD8ab+ T cells",
  "Cycling T cells"
)
iel$celltype <- factor(
  as.character(iel$celltype),
  levels = celltype_order
)


p_celltypes <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",cols = celltype_colors,
  repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_celltypes)
save_plot(p_celltypes, "08_UMAP_major_celltypes", 7, 5.5)

p_celltypes_split <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",
  split.by = "sample", cols = celltype_colors, repel = TRUE, raster = FALSE
) + theme_publication()
save_plot(p_celltypes_split, "09_UMAP_celltypes_split_by_sample", 11, 5.5)
print(p_celltypes_split)

saveRDS(iel,file = "260926 sgms1 iel lpl dimplot.rds")




obj_compare <- iel

DefaultAssay(obj_compare) <- "RNA"

if (inherits(obj_compare[["RNA"]], "Assay5")) {
  obj_compare <- JoinLayers(
    obj_compare,
    assay = "RNA"
  )
}

obj_compare$source <- factor(
  tolower(as.character(obj_compare$sample)),
  levels = c("iel", "lpl")
)

# 检查两组细胞数
print(table(obj_compare$source, useNA = "always"))

expr_sgms1 <- FetchData(
  obj_compare,
  vars = c(
    "Sgms1",
    "source"
  ),
  layer = "data",
  clean = FALSE
)

expr_sgms1 <- expr_sgms1 |>
  mutate(
    source = factor(
      as.character(source),
      levels = c("iel", "lpl")
    )
  ) |>
  filter(
    !is.na(source),
    !is.na(Sgms1)
  )

# 保留原有算法：
# 在整个IEL pool和整个LPL pool内分别选取Sgms1最高的25%细胞
expr_sgms1_top25 <- expr_sgms1 |>
  group_by(source) |>
  arrange(
    desc(Sgms1),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()


sgms1_statistics <- expr_sgms1_top25 |>
  wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    ),
    y.position = 3.35
  )

p_sgms1_top25 <- ggplot(
  expr_sgms1_top25,
  aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.3
  ) +
  stat_pvalue_manual(
    sgms1_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = c(
      "iel" = "#4F86C6",
      "lpl" = "#D64065"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "iel" = "IEL",
      "lpl" = "LPL"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 3.5, by = 0.5),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    ylim = c(1, 3.5)
  ) +
  labs(
    title = "Top 25% Sgms1-expressing cells",
    subtitle = "Top 25% selected separately from pooled IEL and pooled LPL",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

print(p_sgms1_top25)

source_colors <- c(
  "iel" = "#4F86C6",
  "lpl" = "#D64065"
)



save_plot(
  p_sgms1_top25,
  "SGMS1_pooled_IEL_vs_LPL_top25",
  width = 5,
  height = 5
)


sm_genes <- c(
  "Sgms1", "Sgms2",
  "Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2", "Degs1", "Kdsr"
)

set.seed(20260729)
obj_compare<-iel
obj_compare <- AddModuleScore(
  obj_compare,
  features = list(sm_genes),
  assay = "RNA",
  name = "SM_synthesis",
  seed = 20260729
)

sm_table <- obj_compare[[]] |>
  transmute(
    source = factor(
      tolower(as.character(sample)),
      levels = c("iel", "lpl")
    ),
    SM_synthesis_score = SM_synthesis1
  ) |>
  filter(
    !is.na(source),
    !is.na(SM_synthesis_score)
  )

sm_summary <- sm_table |>
  group_by(source) |>
  summarise(
    n_cells = n(),
    mean_score = mean(SM_synthesis_score),
    median_score = median(SM_synthesis_score),
    sd_score = sd(SM_synthesis_score),
    .groups = "drop"
  )

print(sm_summary)

# 不再 group_by(celltype)，只进行一次 IEL vs LPL 比较
sm_statistics <- sm_table |>
  wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  add_xy_position(
    x = "source"
  ) |>
  mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

print(sm_statistics)
source_colors <- c(
  "iel" = "#4F86C6",
  "lpl" = "#D64065"
)

p_sm_score_pooled <- ggplot(
  sm_table,
  aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  stat_pvalue_manual(
    sm_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = source_colors
  ) +
  scale_x_discrete(
    labels = c(
      "iel" = "IEL",
      "lpl" = "LPL"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "SM synthesis module score",
    subtitle = "Pooled IEL versus pooled LPL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

print(p_sm_score_pooled)

save_plot(
  p_sm_score_pooled,
  "SM_synthesis_pooled_IEL_vs_LPL",
  width = 5,
  height = 5
)

save_table(
  sm_summary,
  "SM_synthesis_pooled_IEL_vs_LPL_summary.csv"
)



saveRDS(iel, file.path(object_dir, "iel_lpl_annotated_20260810.rds"), compress = FALSE)

chol_acquisition_genes <- c(
  "Srebf2",
  "Acat2",
  "Hmgcs1",
  "Hmgcr",
  "Sqle",
  "Mbtps1",
  "Mbtps2",
  "Lss",
  "Fdps",
  "Fdft1",
  "Mvd",
  "Mvk",
  "Pmvk",
  "Ggps1",
  "Cyp51",
  "Tm7sf2",
  "Msmo1",
  "Nsdhl",
  "Hsd17b7",
  "Ebp",
  "Sc5d",
  "Dhcr7",
  "Dhcr24",
  "Ldlr",
  "Npc1",
  "Npc2"
)



obj_compare <- AddModuleScore(
  obj_compare,
  features = list(chol_acquisition_genes),
  assay = "RNA",
  name = "Chol_acquisition",
  seed = 20260729
)

chol_table <- obj_compare[[]] |>
  transmute(
    source = factor(
      tolower(as.character(sample)),
      levels = c("iel", "lpl")
    ),
    Chol_acquisition_score = Chol_acquisition1
  ) |>
  filter(
    !is.na(source),
    !is.na(Chol_acquisition_score)
  )

chol_summary <- chol_table |>
  group_by(source) |>
  summarise(
    n_cells = n(),
    mean_score = mean(Chol_acquisition_score),
    median_score = median(Chol_acquisition_score),
    sd_score = sd(Chol_acquisition_score),
    .groups = "drop"
  )

print(chol_summary)

# 不按 celltype 分组，只进行一次 pooled IEL vs pooled LPL 比较
chol_statistics <- chol_table |>
  wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  add_xy_position(
    x = "source"
  ) |>
  mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

print(chol_statistics)

source_colors <- c(
  "iel" = "#4F86C6",
  "lpl" = "#D64065"
)

p_chol_score_pooled <- ggplot(
  chol_table,
  aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  stat_pvalue_manual(
    chol_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = source_colors
  ) +
  scale_x_discrete(
    labels = c(
      "iel" = "IEL",
      "lpl" = "LPL"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Pooled IEL versus pooled LPL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

print(p_chol_score_pooled)

save_plot(
  p_chol_score_pooled,
  "Chol_acquisition_pooled_IEL_vs_LPL",
  width = 5,
  height = 5
)

save_table(
  chol_summary,
  "Chol_acquisition_pooled_IEL_vs_LPL_summary.csv"
)
# =============================================================================
# Simple explicit version: four cell types x three measurements
#
# Run after the updated cell-type annotation in the main script.
# IMPORTANT: remove or comment out the old line: iel <- iel_b
#
# There are no analysis loops and no general analysis/plotting function below.
# Each cell type has its own data, statistics, plot, and output object.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(svglite)
})

# ---- Simple export helper ----------------------------------------------------
# This is the only helper function. It only saves a finished plot.
save_plot_simple <- function(
    plot,
    filename,
    width = 5,
    height = 5
) {
  stem <- file.path(
    figure_dir,
    filename
  )
  
  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()
  
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height
  )
  print(plot)
  grDevices::dev.off()
  
  ggplot2::ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

# ---- Check the current annotation before analysis ----------------------------
print(
  table(
    iel$celltype,
    useNA = "ifany"
  )
)

required_celltypes <- c(
  "Memory-like unconventional T cells",
  "Effector-like unconventional T cells",
  "CD4 T cells",
  "CD8ab+ T cells"
)

missing_celltypes <- setdiff(
  required_celltypes,
  unique(as.character(iel$celltype))
)

if (length(missing_celltypes) > 0) {
  stop(
    "Missing cell types: ",
    paste(missing_celltypes, collapse = ", "),
    ". Check whether iel <- iel_b overwrote the updated annotation."
  )
}

# ---- Prepare object ----------------------------------------------------------
obj_compare <- iel
Seurat::DefaultAssay(obj_compare) <- "RNA"

if (inherits(obj_compare[["RNA"]], "Assay5")) {
  obj_compare <- Seurat::JoinLayers(
    obj_compare,
    assay = "RNA"
  )
}

obj_compare$source <- factor(
  tolower(as.character(obj_compare$sample)),
  levels = c("iel", "lpl")
)

print(
  table(
    obj_compare$celltype,
    obj_compare$source,
    useNA = "ifany"
  )
)

# ---- Colours ----------------------------------------------------------------
niel_memory_colors <- c(
  iel = "#91D0B3",
  lpl = "#C5E8D7"
)

niel_effector_colors <- c(
  iel = "#3E9B72",
  lpl = "#8CCDB1"
)

cd4_colors <- c(
  iel = "#F6AA8D",
  lpl = "#FAD0C0"
)

cd8_colors <- c(
  iel = "#D64065",
  lpl = "#EC9AAF"
)


# =============================================================================
# PART 1: Sgms1 top 25%
# =============================================================================

expr_sgms1 <- Seurat::FetchData(
  obj_compare,
  vars = c(
    "Sgms1",
    "source",
    "celltype"
  ),
  layer = "data",
  clean = FALSE
) |>
  dplyr::filter(
    !is.na(source),
    !is.na(celltype),
    !is.na(Sgms1)
  ) |>
  dplyr::mutate(
    source = factor(
      as.character(source),
      levels = c("iel", "lpl")
    )
  )

# ---- Memory-like nIEL: Sgms1 top 25% ---------------------------------------------
sgms1_niel_memory_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(celltype) == "Memory-like unconventional T cells"
  )

print(table(sgms1_niel_memory_data$source, useNA = "ifany"))

sgms1_niel_memory_top25 <- sgms1_niel_memory_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.25)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_niel_memory_statistics <- sgms1_niel_memory_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_niel_memory <- ggplot2::ggplot(
  sgms1_niel_memory_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 25% Sgms1-expressing cells",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_niel_memory)

save_plot_simple(
  p_sgms1_niel_memory,
  "SGMS1_top25_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: Sgms1 top 25% ---------------------------------------------
sgms1_niel_effector_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(celltype) == "Effector-like unconventional T cells"
  )

print(table(sgms1_niel_effector_data$source, useNA = "ifany"))

sgms1_niel_effector_top25 <- sgms1_niel_effector_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.25)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_niel_effector_statistics <- sgms1_niel_effector_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_niel_effector <- ggplot2::ggplot(
  sgms1_niel_effector_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 25% Sgms1-expressing cells",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_niel_effector)

save_plot_simple(
  p_sgms1_niel_effector,
  "SGMS1_top25_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD4 T cells: Sgms1 top 25% ---------------------------------------------
sgms1_cd4_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(celltype) == "CD4 T cells"
  )

print(table(sgms1_cd4_data$source, useNA = "ifany"))

sgms1_cd4_top25 <- sgms1_cd4_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.25)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_cd4_statistics <- sgms1_cd4_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_cd4 <- ggplot2::ggplot(
  sgms1_cd4_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_cd4_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 25% Sgms1-expressing cells",
    subtitle = "CD4 T cells",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_cd4)

save_plot_simple(
  p_sgms1_cd4,
  "SGMS1_top25_cd4_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD8ab+ T cells: Sgms1 top 25% ---------------------------------------------
sgms1_cd8_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(celltype) == "CD8ab+ T cells"
  )

print(table(sgms1_cd8_data$source, useNA = "ifany"))

sgms1_cd8_top25 <- sgms1_cd8_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.25)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_cd8_statistics <- sgms1_cd8_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_cd8 <- ggplot2::ggplot(
  sgms1_cd8_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 25% Sgms1-expressing cells",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_cd8)

save_plot_simple(
  p_sgms1_cd8,
  "SGMS1_top25_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)

# =============================================================================
# PART 2: SM synthesis module score
# =============================================================================

sm_genes <- c(
  "Sgms1", "Sgms2",
  "Cers1", "Cers2", "Cers3",
  "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2",
  "Degs1", "Kdsr"
)

set.seed(20260729)

obj_compare <- Seurat::AddModuleScore(
  obj_compare,
  features = list(sm_genes),
  assay = "RNA",
  name = "SM_synthesis",
  seed = 20260729
)

sm_all <- obj_compare[[]] |>
  dplyr::transmute(
    source = factor(
      tolower(as.character(source)),
      levels = c("iel", "lpl")
    ),
    celltype = as.character(celltype),
    SM_synthesis_score = SM_synthesis1
  ) |>
  dplyr::filter(
    !is.na(source),
    !is.na(celltype),
    !is.na(SM_synthesis_score)
  )

# ---- Memory-like nIEL: SM synthesis module score -------------------------------------------------
sm_niel_memory_data <- sm_all |>
  dplyr::filter(
    as.character(celltype) == "Memory-like unconventional T cells"
  )

print(table(sm_niel_memory_data$source, useNA = "ifany"))

sm_niel_memory_statistics <- sm_niel_memory_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_niel_memory <- ggplot2::ggplot(
  sm_niel_memory_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_niel_memory)

save_plot_simple(
  p_sm_niel_memory,
  "SM_synthesis_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: SM synthesis module score -------------------------------------------------
sm_niel_effector_data <- sm_all |>
  dplyr::filter(
    as.character(celltype) == "Effector-like unconventional T cells"
  )

print(table(sm_niel_effector_data$source, useNA = "ifany"))

sm_niel_effector_statistics <- sm_niel_effector_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_niel_effector <- ggplot2::ggplot(
  sm_niel_effector_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_niel_effector)

save_plot_simple(
  p_sm_niel_effector,
  "SM_synthesis_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD4 T cells: SM synthesis module score -------------------------------------------------
sm_cd4_data <- sm_all |>
  dplyr::filter(
    as.character(celltype) == "CD4 T cells"
  )

print(table(sm_cd4_data$source, useNA = "ifany"))

sm_cd4_statistics <- sm_cd4_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_cd4 <- ggplot2::ggplot(
  sm_cd4_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_cd4_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_cd4)

save_plot_simple(
  p_sm_cd4,
  "SM_synthesis_cd4_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD8ab+ T cells: SM synthesis module score -------------------------------------------------
sm_cd8_data <- sm_all |>
  dplyr::filter(
    as.character(celltype) == "CD8ab+ T cells"
  )

print(table(sm_cd8_data$source, useNA = "ifany"))

sm_cd8_statistics <- sm_cd8_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_cd8 <- ggplot2::ggplot(
  sm_cd8_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_cd8)

save_plot_simple(
  p_sm_cd8,
  "SM_synthesis_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)

# =============================================================================
# PART 3: Cholesterol acquisition module score
# =============================================================================

chol_acquisition_genes <- c(
  "Srebf2", "Acat2", "Hmgcs1", "Hmgcr", "Sqle",
  "Mbtps1", "Mbtps2", "Lss", "Fdps", "Fdft1",
  "Mvd", "Mvk", "Pmvk", "Ggps1", "Cyp51",
  "Tm7sf2", "Msmo1", "Nsdhl", "Hsd17b7", "Ebp",
  "Sc5d", "Dhcr7", "Dhcr24", "Ldlr", "Npc1", "Npc2"
)

obj_compare <- Seurat::AddModuleScore(
  obj_compare,
  features = list(chol_acquisition_genes),
  assay = "RNA",
  name = "Chol_acquisition",
  seed = 20260729
)

chol_all <- obj_compare[[]] |>
  dplyr::transmute(
    source = factor(
      tolower(as.character(source)),
      levels = c("iel", "lpl")
    ),
    celltype = as.character(celltype),
    Chol_acquisition_score = Chol_acquisition1
  ) |>
  dplyr::filter(
    !is.na(source),
    !is.na(celltype),
    !is.na(Chol_acquisition_score)
  )

# ---- Memory-like nIEL: Cholesterol acquisition module score -------------------------------------------------
chol_niel_memory_data <- chol_all |>
  dplyr::filter(
    as.character(celltype) == "Memory-like unconventional T cells"
  )

print(table(chol_niel_memory_data$source, useNA = "ifany"))

chol_niel_memory_statistics <- chol_niel_memory_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_niel_memory <- ggplot2::ggplot(
  chol_niel_memory_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_niel_memory)

save_plot_simple(
  p_chol_niel_memory,
  "Chol_acquisition_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: Cholesterol acquisition module score -------------------------------------------------
chol_niel_effector_data <- chol_all |>
  dplyr::filter(
    as.character(celltype) == "Effector-like unconventional T cells"
  )

print(table(chol_niel_effector_data$source, useNA = "ifany"))

chol_niel_effector_statistics <- chol_niel_effector_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_niel_effector <- ggplot2::ggplot(
  chol_niel_effector_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_niel_effector)

save_plot_simple(
  p_chol_niel_effector,
  "Chol_acquisition_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD4 T cells: Cholesterol acquisition module score -------------------------------------------------
chol_cd4_data <- chol_all |>
  dplyr::filter(
    as.character(celltype) == "CD4 T cells"
  )

print(table(chol_cd4_data$source, useNA = "ifany"))

chol_cd4_statistics <- chol_cd4_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_cd4 <- ggplot2::ggplot(
  chol_cd4_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_cd4_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_cd4)

save_plot_simple(
  p_chol_cd4,
  "Chol_acquisition_cd4_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- CD8ab+ T cells: Cholesterol acquisition module score -------------------------------------------------
chol_cd8_data <- chol_all |>
  dplyr::filter(
    as.character(celltype) == "CD8ab+ T cells"
  )

print(table(chol_cd8_data$source, useNA = "ifany"))

chol_cd8_statistics <- chol_cd8_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_cd8 <- ggplot2::ggplot(
  chol_cd8_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_cd8)

save_plot_simple(
  p_chol_cd8,
  "Chol_acquisition_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)


###############################################################################
#jem paper iel vs lpl
###############################################################################


set.seed(20260805)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(svglite)
  library(rstatix)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(ggsignif)
})

# ---- Project paths -----------------------------------------------------------
project_dir <- path.expand(
  "~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/jem single cell"
)
figure_dir <- file.path(project_dir, "results", "figures")
table_dir  <- file.path(project_dir, "results", "tables")
object_dir <- file.path(project_dir, "results", "objects")
log_dir    <- file.path(project_dir, "results", "logs")

invisible(lapply(
  c(project_dir, figure_dir, table_dir, object_dir, log_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))


# ---- Reusable publication settings ------------------------------------------
sample_colors <- c(iel = "#4F86C6", lpl = "#D64065")

celltype_colors <- c(
  # CD8aa+ nIEL：绿色系
  "Memory-like CD8aa+ TCRab/gd nIEL" =
    "#91D0B3",
  
  "Effector-like CD8aa+ TCRab/gd nIEL" =
    "#3E9B72",
  
  "CD8ab+ T cells" =
    "#D64065",
  
  # Naive/CD4/Th17：拉开明度和色相
  "Naive T cells" =
    "#F2D56B",  # 浅黄色
  
  "Memory-like CD4 T cells" =
    "#F3A45B",  # 浅橙色
  
  "Effector-like CD4 T cells" =
    "#E66B2E",  # 深橙色
  
  "Th17 cells" =
    "#A63D2D",  # 砖红色
  
  # Treg：紫色
  "Treg cells" =
    "#8E6BBE",
  
  # Cycling：蓝色
  "Cycling T cells" =
    "#4F86C6"
)



theme_publication <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}
theme_set(theme_publication())

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}


theme_publication <- function(
    base_size = 9,
    base_family = "sans"
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}

theme_set(theme_publication())

save_plot <- function(
    plot,
    filename,
    width,
    height,
    units = "in",
    dpi = 600,
    bg = "white"
) {
  stem <- file.path(
    figure_dir,
    safe_name(filename)
  )
  
  # PDF矢量图
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = units,
    bg = bg,
    limitsize = FALSE
  )
  
  # SVG矢量图
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height,
    bg = bg
  )
  print(plot)
  grDevices::dev.off()
  
  # PNG位图
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  
  invisible(plot)
}

save_table <- function(x, filename, row.names = FALSE) {
  write.csv(
    x, file.path(table_dir, filename),
    row.names = row.names, fileEncoding = "UTF-8"
  )
}

# ---- Read, merge, QC ---------------------------------------------------------
iel_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/jem single cell/GSE221006_4samples")

iel <- CreateSeuratObject(iel_counts, project = "iel")



# Mouse mitochondrial genes usually begin "mt-"; ignore.case covers "Mt-".
iel[["percent.mt"]] <- PercentageFeatureSet(
  iel, pattern = "^mt-", assay = "RNA"
)

p_qc_before <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_before)
save_plot(p_qc_before, "01_QC_before_filtering", 10.5, 3.5)

iel <- subset(
  iel,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    nCount_RNA < 20000 &
    percent.mt < 10
)


p_qc_after <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_after)
save_plot(p_qc_after, "02_QC_after_filtering", 10.5, 3.5)

qc_summary <- iel[[]] |>
  group_by(sample) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
save_table(qc_summary, "QC_summary_after_filtering.csv")

# ---- Normalization, Harmony, clustering -------------------------------------
iel <- NormalizeData(iel, verbose = FALSE)
iel <- FindVariableFeatures(iel, nfeatures = 2000, verbose = FALSE)
iel <- ScaleData(iel, features = VariableFeatures(iel), verbose = FALSE)
iel <- RunPCA(iel, features = VariableFeatures(iel), verbose = FALSE)

p_elbow <- ElbowPlot(iel, ndims = 50) + theme_publication()
print(p_elbow)
save_plot(p_elbow, "03_PCA_elbow", 5.5, 4)


dims_main <- 1:30
iel <- FindNeighbors(iel, dims = dims_main, verbose = FALSE)
iel <- FindClusters(iel, resolution = 0.8, random.seed = 20260805, verbose = FALSE)
iel <- RunUMAP(
  iel, dims = dims_main,
  seed.use = 20260805, verbose = FALSE
)


p_clusters <- DimPlot(
  iel, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_clusters)
save_plot(p_clusters, "04_UMAP_clusters", 6, 5)


#############################################################

saveRDS(iel, file.path(object_dir, "iel _jem260916.rds"), compress = FALSE)

marker_panel <- c(
  "Trac", "Trdc","Cd8a","Cd4","Cd8b1","Il17a","Foxp3","Tcf7","Ms4a4b","Il7r","Xcl1","Gm8369","Cd244a","Klf3","Ccl4","Cd200r1","Themis",
  "Cd6","Gzmk","Ccl5","Gzmb","Gzma","Ifng","Mki67","Ccnb2","Stmn1","Foxp3","Cd40lg","Maf","Ctla4","Icos",
  "Tyrobp","Fcer1g","Klrd1","Cd200r2","Frmd5","Cd160","Kit","Id3",
  "Emid1","Pclaf","S1pr1","Lef1","Cd19","Klrb1c","Ncr1","Csf1r","Adgre1","Mpeg1","C1qa","Lyz2"
)
marker_panel <- intersect(marker_panel, rownames(iel))
p_markers <- DotPlot(iel, features = marker_panel) +
  coord_flip() + RotatedAxis() + theme_publication()
print(p_markers)
save_plot(p_markers, "07_cluster_marker_dotplot", 10, 8)


major_map <- c(
  "0"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "1"  = "NT",
  "2"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "3"  = "Memory-like CD8aa+ TCRab/gd nIEL",
  "4"  = "NT",
  "5"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "6"  = "Memory-like CD8aa+ TCRab/gd nIEL",
  "7"  = "Effector-like CD4 T cells",
  "8"  = "CD8ab+ T cells",
  "9"  = "Memory-like CD4 T cells",
  "10" = "CD8ab+ T cells",
  "11" = "NT",
  "12" = "NT",
  "13" = "Naive T cells",
  "14" = "NT",
  "15" = "Effector-like CD8aa+ TCRab/gd nIEL",
  "16" = "NT",
  "17" = "NT",
  "18" = "NT",
  "19" = "NT",
  "20" = "Treg cells",
  "21" = "Th17 cells",
  "22" = "NT",
  "23" = "NT",
  "24" = "Effector-like CD8aa+ TCRab/gd nIEL",
  "25" = "NT",
  "26" = "Cycling T cells",
  "27" = "NT",
  "28" = "NT",
  "29" = "NT"
)

iel$celltype <- unname(major_map[as.character(iel$seurat_clusters)])


saveRDS(iel,file = "260818 jem iel dimplot.rds")

sample_map <- c(
  "1" = "Chow_IntraEpi_1",
  "2" = "HFHS_IntraEpi_1",
  "3" = "Chow_LaminaPropria_1",
  "4" = "HFHS_LaminaPropria_1"
)

barcode_suffix <- sub(
  pattern = ".*-",
  replacement = "",
  x = colnames(iel)
)

iel$sample <- unname(sample_map[barcode_suffix])

# 检查是否存在无法识别的 barcode
if (anyNA(iel$sample)) {
  stop("发现无法匹配到样本的 barcode，请检查 barcode 后缀。")
}

# 保存饮食和组织区室信息，方便后续分析
iel$diet <- ifelse(
  grepl("^Chow", iel$sample),
  "Chow",
  "HFHS"
)

iel$compartment <- ifelse(
  grepl("IntraEpi", iel$sample),
  "IntraEpi",
  "LaminaPropria"
)

# 检查各样本细胞数
print(table(iel$sample))

# ------------------------------------------------------------
# 2. 提取 Chow_IntraEpi_1
# ------------------------------------------------------------

chow_iel_cells <- colnames(iel)[
  iel$sample == "Chow_IntraEpi_1"
]

Chow_IntraEpi_1 <- subset(
  x = iel,
  cells = chow_iel_cells
)

# ------------------------------------------------------------
# 3. 保留并恢复原来的细胞类型身份
# ------------------------------------------------------------

Idents(Chow_IntraEpi_1) <- "celltype"

# 保留原来的细胞类型顺序
celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T cells",
  "Th17 cells",
  "Treg cells",
  "Cycling T cells"
)

existing_celltypes <- intersect(
  celltype_order,
  unique(as.character(Chow_IntraEpi_1$celltype))
)

Chow_IntraEpi_1$celltype <- factor(
  as.character(Chow_IntraEpi_1$celltype),
  levels = existing_celltypes
)

Idents(Chow_IntraEpi_1) <- "celltype"

# ------------------------------------------------------------
# 4. 核对提取结果
# ------------------------------------------------------------



# 使用原来的 UMAP 坐标展示，不重新聚类或重新注释
p_chow_iel <- DimPlot(
  Chow_IntraEpi_1,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  NoLegend() +
  theme_publication()

print(p_chow_iel)

save_plot(
  p_chow_iel,
  "Chow_IntraEpi_1_UMAP_original_annotation",
  width = 7,
  height = 5.5
)


iel <- subset(
  iel,
  subset = celltype != "NT"
)

Idents(iel) <- "celltype"
celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T cells",
  "Th17 cells",
  "Treg cells",
  "Cycling T cells"
)

# 如果已经删除Monocyte，这里不再加入Monocyte
iel$celltype <- factor(
  as.character(iel$celltype),
  levels = celltype_order
)

Idents(iel) <- "celltype"

p_celltypes <- DimPlot(Chow_IntraEpi_1, reduction = "umap", group.by = "celltype",cols = celltype_colors,
                       repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()

print(p_celltypes)
save_plot(p_celltypes, "08_UMAP_major_celltypes", 7, 5.5)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

selected_celltypes <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells"
)

# subset(..., idents = ...)依赖当前Idents
Idents(p_chow_iel) <- "celltype"

objWT1 <- subset(
  x = Chow_IntraEpi_1,
  idents = selected_celltypes
)

expr_sgms1 <- FetchData(
  objWT1,
  vars = c(
    "Sgms1",
    "celltype"
  ),
  layer = "data"
)

plot_celltype_order <- celltype_order[
  celltype_order %in% selected_celltypes
]

expr_sgms1$celltype <- factor(
  as.character(expr_sgms1$celltype),
  levels = plot_celltype_order
)


expr_sgms1_top25 <- expr_sgms1 |>
  filter(!is.na(celltype)) |>
  group_by(celltype) |>
  arrange(
    desc(Sgms1),
    .by_group = TRUE
  ) |>
  mutate(
    expression_rank = row_number(),
    total_cells_in_group = n(),
    top25_cell_number = ceiling(
      total_cells_in_group * 0.25
    )
  ) |>
  filter(
    expression_rank <= top25_cell_number
  ) |>
  ungroup()

target_group <-
  "Effector-like CD8aa+ TCRab/gd nIEL"

comparisons <- list(
  c(
    target_group,
    "Memory-like CD8aa+ TCRab/gd nIEL"
  ),
  c(
    target_group,
    "CD8ab+ T cells"
  ),
  c(
    target_group,
    "Memory-like CD4 T cells"
  ),
  c(
    target_group,
    "Effector-like CD4 T cells"
  )
)


plot_colors <- celltype_colors[plot_celltype_order]

p_sgms1 <- ggplot(
  expr_sgms1_top25,
  aes(
    x = celltype,
    y = Sgms1,
    fill = celltype
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4,
    step.increase = 0.08,
    hide.ns = FALSE
  ) +
  scale_fill_manual(
    values = plot_colors,
    breaks = plot_celltype_order
  ) +
  scale_x_discrete(
    limits = plot_celltype_order
  ) +
  labs(
    title = "Sgms1",
    x = NULL,
    y = "Normalized expression"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none"
  )

print(p_sgms1)

save_plot(
  p_sgms1,
  "Sgms1_top35_Effector_CD8aa_comparison",
  width = 10,
  height = 7
)



sm_genes <- c(
  "Sgms1", "Sgms2", "Cers1", "Cers2", "Cers3", "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2", "Degs1", "Kdsr"
)

objWT1<- AddModuleScore(
  objWT1, features = list(sm_genes), name = "SM_synthesis",
  seed = 20260729
)
sm_table <- objWT1[[]] |>
  transmute(celltype = celltype, SM_synthesis_score = SM_synthesis1)
save_table(sm_table, "WT_SM_synthesis_module_scores.csv")

p_sm_score <- ggplot(
  sm_table, aes(celltype, SM_synthesis_score, fill = celltype)
) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.65) +
  stat_compare_means(
    comparisons = comparisons, method = "wilcox.test", size = 3
  ) +
  scale_fill_manual(
    values = celltype_colors) +
  labs(title = "SM synthesis module score", x = NULL, y = "Module score") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
print(p_sm_score)
save_plot(p_sm_score, "16_WT_SM_synthesis_module_score", 7, 5.5)


# =============================================================================
# JEM normal Chow IEL versus normal Chow LPL: five groups x three measurements
#
# Run after cell-type annotation in the JEM script.
# Barcode suffix 1 = Chow IEL; suffix 3 = Chow LPL; suffix 2/4 are excluded.
#
# There are no analysis loops and no general analysis/plotting function below.
# Each cell type has its own data, statistics, plot, and output object.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(svglite)
})

# ---- Simple export helper ----------------------------------------------------
# This is the only helper function. It only saves a finished plot.
save_plot_simple <- function(
    plot,
    filename,
    width = 5,
    height = 5
) {
  stem <- file.path(
    figure_dir,
    filename
  )
  
  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()
  
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height
  )
  print(plot)
  grDevices::dev.off()
  
  ggplot2::ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

# =============================================================================
# PART 0: identify normal IEL and LPL from barcode suffix
# =============================================================================

obj_compare <- iel
Seurat::DefaultAssay(obj_compare) <- "RNA"

if (inherits(obj_compare[["RNA"]], "Assay5")) {
  obj_compare <- Seurat::JoinLayers(
    obj_compare,
    assay = "RNA"
  )
}

obj_compare$barcode <- colnames(obj_compare)

obj_compare$barcode_suffix <- sub(
  pattern = ".*-",
  replacement = "",
  x = obj_compare$barcode
)

print(
  table(
    obj_compare$barcode_suffix,
    useNA = "ifany"
  )
)

# Only the two normal Chow samples receive a source label.
# -1 = normal IEL; -3 = normal LPL.
# HFHS suffixes -2 and -4 stay NA and are excluded.
obj_compare$source <- NA_character_

obj_compare$source[
  obj_compare$barcode_suffix == "1"
] <- "iel"

obj_compare$source[
  obj_compare$barcode_suffix == "3"
] <- "lpl"

print(
  table(
    obj_compare$source,
    useNA = "ifany"
  )
)

normal_cells <- colnames(obj_compare)[
  !is.na(obj_compare$source)
]

obj_compare <- subset(
  obj_compare,
  cells = normal_cells
)

obj_compare$source <- factor(
  obj_compare$source,
  levels = c("iel", "lpl")
)

# Keep every requested group exactly as originally annotated.
# Memory-like and effector-like CD4 cells remain separate.
obj_compare$analysis_group <- as.character(
  obj_compare$celltype
)

required_groups <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "CD8ab+ T cells"
)

missing_groups <- setdiff(
  required_groups,
  unique(obj_compare$analysis_group)
)

if (length(missing_groups) > 0) {
  stop(
    "Missing analysis groups: ",
    paste(missing_groups, collapse = ", ")
  )
}

print(
  table(
    obj_compare$analysis_group,
    obj_compare$source,
    useNA = "ifany"
  )
)

# ---- Colours ----------------------------------------------------------------
niel_memory_colors <- c(
  iel = "#91D0B3",
  lpl = "#C5E8D7"
)

niel_effector_colors <- c(
  iel = "#3E9B72",
  lpl = "#8CCDB1"
)

cd4_memory_colors <- c(
  iel = "#F3A45B",
  lpl = "#F8CFA8"
)

cd4_effector_colors <- c(
  iel = "#D96B3B",
  lpl = "#F2AD8C"
)

cd8_colors <- c(
  iel = "#D64065",
  lpl = "#EC9AAF"
)


# =============================================================================
# PART 1: Sgms1 top 25%
# =============================================================================

expr_sgms1 <- Seurat::FetchData(
  obj_compare,
  vars = c(
    "Sgms1",
    "source",
    "analysis_group"
  ),
  layer = "data",
  clean = FALSE
) |>
  dplyr::filter(
    !is.na(source),
    !is.na(analysis_group),
    !is.na(Sgms1)
  ) |>
  dplyr::mutate(
    source = factor(
      as.character(source),
      levels = c("iel", "lpl")
    )
  )

# ---- Memory-like nIEL: Sgms1 top 25% ---------------------------------------------
sgms1_niel_memory_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD8aa+ TCRab/gd nIEL"
  )

print(table(sgms1_niel_memory_data$source, useNA = "ifany"))

sgms1_niel_memory_top25 <- sgms1_niel_memory_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.35)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_niel_memory_statistics <- sgms1_niel_memory_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_niel_memory <- ggplot2::ggplot(
  sgms1_niel_memory_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 35% Sgms1-expressing cells",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_niel_memory)

save_plot_simple(
  p_sgms1_niel_memory,
  "JEM_Chow_SGMS1_top25_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: Sgms1 top 25% ---------------------------------------------
sgms1_niel_effector_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD8aa+ TCRab/gd nIEL"
  )

print(table(sgms1_niel_effector_data$source, useNA = "ifany"))

sgms1_niel_effector_top25 <- sgms1_niel_effector_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.35)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_niel_effector_statistics <- sgms1_niel_effector_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_niel_effector <- ggplot2::ggplot(
  sgms1_niel_effector_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 35% Sgms1-expressing cells",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_niel_effector)

save_plot_simple(
  p_sgms1_niel_effector,
  "JEM_Chow_SGMS1_top35_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Memory-like CD4 T cells: Sgms1 top 25% ---------------------------------
sgms1_cd4_memory_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD4 T cells"
  )

print(table(sgms1_cd4_memory_data$source, useNA = "ifany"))

sgms1_cd4_memory_top25 <- sgms1_cd4_memory_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.35)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_cd4_memory_statistics <- sgms1_cd4_memory_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_cd4_memory <- ggplot2::ggplot(
  sgms1_cd4_memory_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_cd4_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 35% Sgms1-expressing cells",
    subtitle = "Memory-like CD4 T cells",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_cd4_memory)

save_plot_simple(
  p_sgms1_cd4_memory,
  "JEM_Chow_SGMS1_top25_cd4_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- Effector-like CD4 T cells: Sgms1 top 25% -------------------------------
sgms1_cd4_effector_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD4 T cells"
  )

print(table(sgms1_cd4_effector_data$source, useNA = "ifany"))

sgms1_cd4_effector_top25 <- sgms1_cd4_effector_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.35)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_cd4_effector_statistics <- sgms1_cd4_effector_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_cd4_effector <- ggplot2::ggplot(
  sgms1_cd4_effector_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_cd4_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 35% Sgms1-expressing cells",
    subtitle = "Effector-like CD4 T cells",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_cd4_effector)

save_plot_simple(
  p_sgms1_cd4_effector,
  "JEM_Chow_SGMS1_top25_cd4_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- CD8ab+ T cells: Sgms1 top 25% ---------------------------------------------
sgms1_cd8_data <- expr_sgms1 |>
  dplyr::filter(
    as.character(analysis_group) == "CD8ab+ T cells"
  )

print(table(sgms1_cd8_data$source, useNA = "ifany"))

sgms1_cd8_top25 <- sgms1_cd8_data |>
  dplyr::group_by(source) |>
  dplyr::arrange(
    dplyr::desc(Sgms1),
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    expression_rank = dplyr::row_number(),
    total_cells_in_group = dplyr::n(),
    top25_cell_number = ceiling(total_cells_in_group * 0.35)
  ) |>
  dplyr::filter(
    expression_rank <= top25_cell_number
  ) |>
  dplyr::ungroup()

sgms1_cd8_statistics <- sgms1_cd8_top25 |>
  rstatix::wilcox_test(
    Sgms1 ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sgms1_cd8 <- ggplot2::ggplot(
  sgms1_cd8_top25,
  ggplot2::aes(
    x = source,
    y = Sgms1,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sgms1_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Top 35% Sgms1-expressing cells",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Normalized Sgms1 expression"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sgms1_cd8)

save_plot_simple(
  p_sgms1_cd8,
  "JEM_Chow_SGMS1_top25_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)

# =============================================================================
# PART 2: SM synthesis module score
# =============================================================================

sm_genes <- c(
  "Sgms1", "Sgms2",
  "Cers1", "Cers2", "Cers3",
  "Cers4", "Cers5", "Cers6",
  "Sptlc1", "Sptlc2",
  "Degs1", "Kdsr"
)

set.seed(20260729)

obj_compare <- Seurat::AddModuleScore(
  obj_compare,
  features = list(sm_genes),
  assay = "RNA",
  name = "SM_synthesis",
  seed = 20260729
)

sm_all <- obj_compare[[]] |>
  dplyr::transmute(
    source = factor(
      tolower(as.character(source)),
      levels = c("iel", "lpl")
    ),
    analysis_group = as.character(analysis_group),
    SM_synthesis_score = SM_synthesis1
  ) |>
  dplyr::filter(
    !is.na(source),
    !is.na(analysis_group),
    !is.na(SM_synthesis_score)
  )

# ---- Memory-like nIEL: SM synthesis module score -------------------------------------------------
sm_niel_memory_data <- sm_all |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD8aa+ TCRab/gd nIEL"
  )

print(table(sm_niel_memory_data$source, useNA = "ifany"))

sm_niel_memory_statistics <- sm_niel_memory_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_niel_memory <- ggplot2::ggplot(
  sm_niel_memory_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_niel_memory)

save_plot_simple(
  p_sm_niel_memory,
  "JEM_Chow_SM_synthesis_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: SM synthesis module score -------------------------------------------------
sm_niel_effector_data <- sm_all |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD8aa+ TCRab/gd nIEL"
  )

print(table(sm_niel_effector_data$source, useNA = "ifany"))

sm_niel_effector_statistics <- sm_niel_effector_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_niel_effector <- ggplot2::ggplot(
  sm_niel_effector_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_niel_effector)

save_plot_simple(
  p_sm_niel_effector,
  "JEM_Chow_SM_synthesis_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Memory-like CD4 T cells: SM synthesis module score ---------------------
sm_cd4_memory_data <- sm_all |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD4 T cells"
  )

print(table(sm_cd4_memory_data$source, useNA = "ifany"))

sm_cd4_memory_statistics <- sm_cd4_memory_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_cd4_memory <- ggplot2::ggplot(
  sm_cd4_memory_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_cd4_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Memory-like CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_cd4_memory)

save_plot_simple(
  p_sm_cd4_memory,
  "JEM_Chow_SM_synthesis_cd4_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- Effector-like CD4 T cells: SM synthesis module score -------------------
sm_cd4_effector_data <- sm_all |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD4 T cells"
  )

print(table(sm_cd4_effector_data$source, useNA = "ifany"))

sm_cd4_effector_statistics <- sm_cd4_effector_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_cd4_effector <- ggplot2::ggplot(
  sm_cd4_effector_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_cd4_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "Effector-like CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_cd4_effector)

save_plot_simple(
  p_sm_cd4_effector,
  "JEM_Chow_SM_synthesis_cd4_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- CD8ab+ T cells: SM synthesis module score -------------------------------------------------
sm_cd8_data <- sm_all |>
  dplyr::filter(
    as.character(analysis_group) == "CD8ab+ T cells"
  )

print(table(sm_cd8_data$source, useNA = "ifany"))

sm_cd8_statistics <- sm_cd8_data |>
  rstatix::wilcox_test(
    SM_synthesis_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_sm_cd8 <- ggplot2::ggplot(
  sm_cd8_data,
  ggplot2::aes(
    x = source,
    y = SM_synthesis_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    sm_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "SM synthesis module score",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_sm_cd8)

save_plot_simple(
  p_sm_cd8,
  "JEM_Chow_SM_synthesis_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)

# =============================================================================
# PART 3: Cholesterol acquisition module score
# =============================================================================

chol_acquisition_genes <- c(
  "Srebf2", "Acat2", "Hmgcs1", "Hmgcr", "Sqle",
  "Mbtps1", "Mbtps2", "Lss", "Fdps", "Fdft1",
  "Mvd", "Mvk", "Pmvk", "Ggps1", "Cyp51",
  "Tm7sf2", "Msmo1", "Nsdhl", "Hsd17b7", "Ebp",
  "Sc5d", "Dhcr7", "Dhcr24", "Ldlr", "Npc1", "Npc2"
)

obj_compare <- Seurat::AddModuleScore(
  obj_compare,
  features = list(chol_acquisition_genes),
  assay = "RNA",
  name = "Chol_acquisition",
  seed = 20260729
)

chol_all <- obj_compare[[]] |>
  dplyr::transmute(
    source = factor(
      tolower(as.character(source)),
      levels = c("iel", "lpl")
    ),
    analysis_group = as.character(analysis_group),
    Chol_acquisition_score = Chol_acquisition1
  ) |>
  dplyr::filter(
    !is.na(source),
    !is.na(analysis_group),
    !is.na(Chol_acquisition_score)
  )

# ---- Memory-like nIEL: Cholesterol acquisition module score -------------------------------------------------
chol_niel_memory_data <- chol_all |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD8aa+ TCRab/gd nIEL"
  )

print(table(chol_niel_memory_data$source, useNA = "ifany"))

chol_niel_memory_statistics <- chol_niel_memory_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_niel_memory <- ggplot2::ggplot(
  chol_niel_memory_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_niel_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Memory-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_niel_memory)

save_plot_simple(
  p_chol_niel_memory,
  "JEM_Chow_Chol_acquisition_niel_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Effector-like nIEL: Cholesterol acquisition module score -------------------------------------------------
chol_niel_effector_data <- chol_all |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD8aa+ TCRab/gd nIEL"
  )

print(table(chol_niel_effector_data$source, useNA = "ifany"))

chol_niel_effector_statistics <- chol_niel_effector_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_niel_effector <- ggplot2::ggplot(
  chol_niel_effector_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_niel_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = niel_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Effector-like nIEL",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_niel_effector)

save_plot_simple(
  p_chol_niel_effector,
  "JEM_Chow_Chol_acquisition_niel_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)


# ---- Memory-like CD4 T cells: Cholesterol acquisition module score ----------
chol_cd4_memory_data <- chol_all |>
  dplyr::filter(
    as.character(analysis_group) == "Memory-like CD4 T cells"
  )

print(table(chol_cd4_memory_data$source, useNA = "ifany"))

chol_cd4_memory_statistics <- chol_cd4_memory_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_cd4_memory <- ggplot2::ggplot(
  chol_cd4_memory_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_cd4_memory_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_memory_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Memory-like CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_cd4_memory)

save_plot_simple(
  p_chol_cd4_memory,
  "JEM_Chow_Chol_acquisition_cd4_memory_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- Effector-like CD4 T cells: Cholesterol acquisition module score --------
chol_cd4_effector_data <- chol_all |>
  dplyr::filter(
    as.character(analysis_group) == "Effector-like CD4 T cells"
  )

print(table(chol_cd4_effector_data$source, useNA = "ifany"))

chol_cd4_effector_statistics <- chol_cd4_effector_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_cd4_effector <- ggplot2::ggplot(
  chol_cd4_effector_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_cd4_effector_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd4_effector_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "Effector-like CD4 T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_cd4_effector)

save_plot_simple(
  p_chol_cd4_effector,
  "JEM_Chow_Chol_acquisition_cd4_effector_IEL_vs_LPL",
  width = 5,
  height = 5
)




# ---- CD8ab+ T cells: Cholesterol acquisition module score -------------------------------------------------
chol_cd8_data <- chol_all |>
  dplyr::filter(
    as.character(analysis_group) == "CD8ab+ T cells"
  )

print(table(chol_cd8_data$source, useNA = "ifany"))

chol_cd8_statistics <- chol_cd8_data |>
  rstatix::wilcox_test(
    Chol_acquisition_score ~ source,
    exact = FALSE
  ) |>
  rstatix::add_xy_position(
    x = "source"
  ) |>
  dplyr::mutate(
    p_label = paste0(
      "p = ",
      format.pval(
        p,
        digits = 3,
        eps = 1e-300
      )
    )
  )

p_chol_cd8 <- ggplot2::ggplot(
  chol_cd8_data,
  ggplot2::aes(
    x = source,
    y = Chol_acquisition_score,
    fill = source
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.75,
    linewidth = 0.35
  ) +
  ggpubr::stat_pvalue_manual(
    chol_cd8_statistics,
    label = "p_label",
    hide.ns = FALSE,
    tip.length = 0.01
  ) +
  ggplot2::scale_fill_manual(
    values = cd8_colors
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      iel = "IEL",
      lpl = "LPL"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  ggplot2::labs(
    title = "Cholesterol acquisition module score",
    subtitle = "CD8ab+ T cells",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold")
  )

print(p_chol_cd8)

save_plot_simple(
  p_chol_cd8,
  "JEM_Chow_Chol_acquisition_cd8_IEL_vs_LPL",
  width = 5,
  height = 5
)


#---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "Memory-like CD4 T cells",
  "Effector-like CD4 T cells",
  "Naive T cells",
  "Th17 cells",
  "Treg cells",
  "Cycling T cells"
)

marker_genes <- c(
  "Trac",
  "Tcrg-C1",
  "Cd8a",
  "Cd4",
  "Cd8b1",
  "Tcf7",
  "Gzmb",
  "S1pr1",
  "Mki67",
  "Foxp3",
  "Il17a"
)

p_feature_umap <- FeaturePlot(
  object = iel,
  features = marker_genes,
  reduction = "umap",
  
  cols = c(
    "#D3D3D3",
    "#004C99"
  ),
  
  # 强化阳性细胞颜色
  min.cutoff = 0,
  max.cutoff = "q50",
  
  order = TRUE,
  pt.size = 0.6,
  
  ncol = 4,
  raster = TRUE,
  raster.dpi = c(300, 300),
  combine = TRUE
)
p_feature_umap <- p_feature_umap &
  theme_publication(
    base_size = 9,
    base_family = "sans"
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(
      face = "italic",
      size = 10,
      hjust = 0.5,
      colour = "black"
    ),
    
    legend.position = "none"
  )
p_feature_umap

save_plot(
  plot = p_feature_umap,
  filename = "UMAP_T_cell_marker_FeaturePlot",
  width = 10,
  height = 7.5,
  dpi = 600,
  bg = "white"
)


gene_colors <- c(
  "Trac"    = "#4C78A8",  # 蓝色
  "Tcrg-C1" = "#72A98F",  # 灰绿色
  "Cd8a"    = "#D98C95",  # 粉色
  "Cd4"     = "#E9C95B",  # 黄色
  "Cd8b1"   = "#C9365B",  # 玫红
  "Tcf7"    = "#F27658",  # 珊瑚橙
  "Gzmb"    = "#70B7C4",  # 青蓝
  "S1pr1"   = "#B6A6A0",  # 暖灰
  "Mki67"   = "#E99662",  # 浅橙
  "Il17a"   = "#8E6BBE",  # 紫色
  "Foxp3"   = "#58A65C"   # 绿色
)

stacked_violin_plot <- function(
    object,
    features,
    celltype_col = "celltype",
    celltype_order,
    assay = "RNA"
) {
  DefaultAssay(object) <- assay
  
  # 检查基因是否存在
  missing_features <- setdiff(features, rownames(object))
  
  if (length(missing_features) > 0) {
    stop(
      "以下基因未在当前 assay 中找到：",
      paste(missing_features, collapse = ", ")
    )
  }
  
  plot_data <- FetchData(
    object = object,
    vars = c(celltype_col, features)
  )
  
  plot_data$celltype <- factor(
    plot_data[[celltype_col]],
    levels = celltype_order
  )
  
  plot_data <- plot_data |>
    filter(!is.na(celltype)) |>
    select(celltype, all_of(features)) |>
    pivot_longer(
      cols = all_of(features),
      names_to = "feature",
      values_to = "expression"
    )
  
  # 控制基因从上至下的顺序
  plot_data$feature <- factor(
    plot_data$feature,
    levels = features
  )
  
  ggplot(
    plot_data,
    aes(
      x = celltype,
      y = expression,
      fill = feature,
      group = interaction(feature, celltype)
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "#D8D8D8",
      linewidth = 0.4
    ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      width = 0.86,
      adjust = 1,
      colour = "#555555",
      linewidth = 0.3
    ) +
    facet_grid(
      rows = vars(feature),
      scales = "free_y",
      switch = "y"
    ) +
    scale_x_discrete(
      limits = celltype_order,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = gene_colors,
      drop = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # 细胞群对应的纵向参考线
      panel.grid.major.x = element_line(
        colour = "#E5E5E5",
        linewidth = 0.4
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      
      axis.text.x = element_text(
        colour = "#333333",
        angle = 48,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
      axis.text.y = element_blank(),
      
      # 在左侧显示基因名称
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        colour = "#333333",
        angle = 0,
        hjust = 1,
        face = "italic",
        size = 9
      ),
      
      panel.spacing.y = grid::unit(0, "pt"),
      legend.position = "none",
      
      plot.margin = margin(
        t = 15,
        r = 15,
        b = 15,
        l = 15
      )
    )
}

p_stacked_violin <- stacked_violin_plot(
  object = iel,
  features = marker_genes,
  celltype_col = "celltype",
  celltype_order = celltype_order,
  assay = "RNA"
)
p_stacked_violin 
save_plot(
  plot = p_stacked_violin,
  filename = "stacked_violin_T_cell_markers",
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)

###############################################################################
# cyp27a1 KO IEL single-cell RNA-seq analysis

###############################################################################

set.seed(20260729)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(svglite)
  library(ggsignif)
  library(rstatix)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(biomaRt)
  library(clusterProfiler)
  library(enrichplot)
})

# ---- Project paths -----------------------------------------------------------
project_dir <- path.expand(
  "~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/cyp27a1 ko iel ileum only"
)
figure_dir <- file.path(project_dir, "results", "figures")
table_dir  <- file.path(project_dir, "results", "tables")
object_dir <- file.path(project_dir, "results", "objects")
log_dir    <- file.path(project_dir, "results", "logs")

invisible(lapply(
  c(project_dir, figure_dir, table_dir, object_dir, log_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))


# ---- Reusable publication settings ------------------------------------------
sample_colors <- c(iwt = "#4F86C6", iko = "#D64065")


celltype_colors <- c(
  # CD8aa+ TCRab/gd nIEL：绿色系
  "Memory-like CD8aa+ TCRab/gd nIEL" =
    "#91D0B3",
  
  "Effector-like CD8aa+ TCRab/gd nIEL" =
    "#3E9B72",
  
  "CD8ab+ T cells" =
    "#D64065",
  
  # CD4 谱系：橙色系
  "CD4 T cells" =
    "#F8845C",
  
  "Th17 cells" =
    "#E5A34A",
  
  # Naive：柔和黄色
  "Naive T cells" =
    "#F5D983",
  
  # Treg：紫色
  "Treg cells" =
    "#8E6BBE",
  
  # Cycling：蓝色
  "Cycling T cells" =
    "#4F86C6"
)


theme_publication <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}
theme_set(theme_publication())

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}


theme_publication <- function(
    base_size = 9,
    base_family = "sans"
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}

theme_set(theme_publication())

save_plot <- function(
    plot,
    filename,
    width,
    height,
    units = "in",
    dpi = 600,
    bg = "white"
) {
  stem <- file.path(
    figure_dir,
    safe_name(filename)
  )
  
  # PDF矢量图
  ggsave(
    filename = paste0(stem, ".pdf"),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = units,
    bg = bg,
    limitsize = FALSE
  )
  
  # SVG矢量图
  svglite::svglite(
    file = paste0(stem, ".svg"),
    width = width,
    height = height,
    bg = bg
  )
  print(plot)
  grDevices::dev.off()
  
  # PNG位图
  ggsave(
    filename = paste0(stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  
  invisible(plot)
}

save_table <- function(x, filename, row.names = FALSE) {
  write.csv(
    x, file.path(table_dir, filename),
    row.names = row.names, fileEncoding = "UTF-8"
  )
}

# ---- Read, merge, QC ---------------------------------------------------------
iwt_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/cyp27a1 ko iel ileum only/IIEL_wt/outs")
iko_counts<-Read10X(data.dir ="~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/cyp27a1 ko iel ileum only/IIEL_ko/outs")

iwt <- CreateSeuratObject(iwt_counts, project = "iwt")
iko <- CreateSeuratObject(iko_counts, project = "iko")

iel <- merge(iwt,iko,
             add.cell.ids = c("iwt","iko"),
             project = "total") 

iel$sample <- factor(iel$orig.ident, levels = c("iwt", "iko"))

# Mouse mitochondrial genes usually begin "mt-"; ignore.case covers "Mt-".
iel[["percent.mt"]] <- PercentageFeatureSet(
  iel, pattern = "^mt-", assay = "RNA"
)

p_qc_before <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_before)
save_plot(p_qc_before, "01_QC_before_filtering", 10.5, 3.5)

iel <- subset(
  iel,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    nCount_RNA < 10000 &
    percent.mt < 5
)

p_qc_after <- VlnPlot(
  iel,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample", ncol = 3, pt.size = 0
) & scale_fill_manual(values = sample_colors) & theme_publication()
print(p_qc_after)
save_plot(p_qc_after, "02_QC_after_filtering", 10.5, 3.5)

qc_summary <- iel[[]] |>
  group_by(sample) |>
  summarise(
    cells = n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )
save_table(qc_summary, "QC_summary_after_filtering.csv")

# ---- Normalization, Harmony, clustering -------------------------------------
iel <- NormalizeData(iel, verbose = FALSE)
iel <- FindVariableFeatures(iel, nfeatures = 2000, verbose = FALSE)
iel <- ScaleData(iel, features = VariableFeatures(iel), verbose = FALSE)
iel <- RunPCA(iel, features = VariableFeatures(iel), verbose = FALSE)

p_elbow <- ElbowPlot(iel, ndims = 50) + theme_publication()
print(p_elbow)
save_plot(p_elbow, "03_PCA_elbow", 5.5, 4)

# Harmony is run once. Downstream neighbors and UMAP explicitly use harmony.
iel <- RunHarmony(
  iel, group.by.vars = "sample", reduction = "pca",
  reduction.save = "harmony", assay.use = "RNA", verbose = FALSE
)

dims_main <- 1:30
iel <- FindNeighbors(iel, reduction = "harmony", dims = dims_main, verbose = FALSE)
iel <- FindClusters(iel, resolution = 0.5, random.seed = 20260729, verbose = FALSE)
iel <- RunUMAP(
  iel, reduction = "harmony", dims = dims_main,
  seed.use = 20260729, verbose = FALSE
)

p_clusters <- DimPlot(
  iel, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_clusters)
save_plot(p_clusters, "04_UMAP_clusters", 6, 5)

p_samples <- DimPlot(
  iel, reduction = "umap", group.by = "sample",
  cols = sample_colors, raster = FALSE
) + theme_publication()
print(p_samples)
save_plot(p_samples, "05_UMAP_samples", 6, 5)

p_split <- DimPlot(
  iel, reduction = "umap", group.by = "seurat_clusters",
  split.by = "sample", label = TRUE, repel = TRUE, raster = FALSE
) + theme_publication()
print(p_split)
save_plot(p_split, "06_UMAP_clusters_split_by_sample", 10, 5)
#############################################################


# =============================================================================
# 计算每个 cluster 的特征基因
# =============================================================================

# 建立单独对象，避免改变主对象iel
iel_markers <- iel

DefaultAssay(iel_markers) <- "RNA"

# Seurat v5中，差异分析前合并RNA assay的不同layers
if (inherits(iel_markers[["RNA"]], "Assay5")) {
  iel_markers <- JoinLayers(
    iel_markers,
    assay = "RNA"
  )
}

# 将细胞身份设置为cluster
Idents(iel_markers) <- "seurat_clusters"

# 计算每个cluster相对于其他全部细胞的上调特征基因
cluster_markers <- FindAllMarkers(
  object = iel_markers,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# 判断Seurat版本使用的fold-change列名
fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cluster_markers)
)[1]

if (is.na(fc_column)) {
  stop("FindAllMarkers结果中没有找到fold-change列。")
}

# 按cluster、显著性和fold change排序
cluster_markers <- cluster_markers |>
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]])
  )

# 保存所有cluster的全部marker
save_table(
  cluster_markers,
  "all_cluster_marker_genes.csv"
)


top15_cluster_markers <- cluster_markers |>
  filter(
    p_val_adj < 0.05
  ) |>
  group_by(cluster) |>
  slice_max(
    order_by = .data[[fc_column]],
    n = 15,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(
    cluster,
    desc(.data[[fc_column]])
  )

save_table(
  top15_cluster_markers,
  "top15_marker_genes_per_cluster.csv"
)

#==============================================================

marker_panel <- c(
  "Trac", "Trdc","Cd8a","Cd4","Cd8b1","Il17a","Foxp3","Tcf7","Ms4a4b","Il7r","Xcl1","Gm8369","Cd244a","Klf3","Ccl4","Cd200r1","Themis",
  "Cd6","Gzmk","Ccl5","Gzmb","Gzma","Ifng","Mki67","Ccnb2","Stmn1","Foxp3","Cd40lg","Maf","Ctla4","Icos",
  "Tyrobp","Fcer1g","Klrd1","Cd200r2","Frmd5","Cd160","Kit","Id3",
  "Emid1","Pclaf","S1pr1","Lef1","Cd19","Klrb1c","Ncr1","Csf1r","Adgre1","Mpeg1","C1qa","Lyz2"
)
marker_panel <- intersect(marker_panel, rownames(iel))
p_markers <- DotPlot(iel, features = marker_panel) +
  coord_flip() + RotatedAxis() + theme_publication()
print(p_markers)
save_plot(p_markers, "07_cluster_marker_dotplot", 10, 8)

#260729 19:00rds
saveRDS(iel,file = "260810 cyp27a1 ko iel.rds")
# ---- Major cell-type annotation ---------------------------------------------
# IMPORTANT: verify these cluster-to-cell-type mappings from marker plots before
# biological interpretation. Named assignment avoids fragile numeric column use.

major_map <- c(
  "0"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "1"  = "Memory-like CD8aa+ TCRab/gd nIEL",
  "2"  = "CD8ab+ T cells",
  "3"  = "CD8ab+ T cells",
  "4"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "5"  = "CD4 T cells",
  "6"  = "Effector-like CD8aa+ TCRab/gd nIEL",
  "7"  = "Memory-like CD8aa+ TCRab/gd nIEL",
  "8"  = "Cycling T cells",
  "9"  = "Naive T cells",
  "10" = "Th17 cells",
  "11" = "Treg cells",
  "12" = "Monocyte-like / T-monocyte doublets"
)
iel$celltype <- unname(major_map[as.character(iel$seurat_clusters)])


saveRDS(iel,file = "260803 sgms1 ko iel dimplot.rds")
iel_orig<-iel
iel<-iel_orig
iel <- subset(
  iel,
  subset = celltype != "Monocyte-like / T-monocyte doublets"
)

Idents(iel) <- "celltype"
celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "CD4 T cells",
  "Naive T cells",
  "Th17 cells",
  "Treg cells",
  "Cycling T cells"
)

# 如果已经删除Monocyte，这里不再加入Monocyte
iel$celltype <- factor(
  as.character(iel$celltype),
  levels = celltype_order
)

Idents(iel) <- "celltype"

p_celltypes <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",cols = celltype_colors,
  repel = TRUE, raster = FALSE
) + NoLegend() + theme_publication()
print(p_celltypes)
save_plot(p_celltypes, "08_UMAP_major_celltypes", 7, 5.5)

p_celltypes_split <- DimPlot(
  iel, reduction = "umap", group.by = "celltype",
  split.by = "sample", cols = celltype_colors, repel = TRUE, raster = FALSE
) + theme_publication()
save_plot(p_celltypes_split, "09_UMAP_celltypes_split_by_sample", 11, 5.5)

cell_counts <- as.data.frame(table(
  celltype = iel$celltype, sample = iel$sample
)) |>
  filter(Freq > 0) |>
  group_by(sample) |>
  mutate(ratio = Freq / sum(Freq)) |>
  ungroup()
print(cell_counts)
save_table(cell_counts, "major_celltype_counts_and_ratios.csv")

p_ratio <- ggplot(cell_counts, aes(sample, ratio, fill = celltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell proportion", fill = "Cell type") +
  theme_publication()
print(p_ratio)
save_plot(p_ratio, "10_major_celltype_proportions", 6.5, 5)


saveRDS(iel,file = "260917 cyp27a1 iel ko dimplot.rds")




HCdegs <- read.csv("~/Life in Lab/大吉大利/27 & sm & intestine/single cell Summary/cyp27a1 ko iel ileum only/YCS_27HCvsCtrl_20230214.csv")


HCDown <- HCdegs[which(HCdegs$log2FoldChange <= -3 & HCdegs$padj <= 0.05),]
HCDown <- HCDown$hgnc_symbol
HCDown <- HCDown[which(!duplicated(HCDown))]


HCDown_mm <- homologene(HCDown, inTax = 9606, outTax = 10090)

iel <- AddModuleScore(iel,features = list(HCDown_mm$"10090"),
                      name = "HCDown",
                      seed = 20260729
)

hcdown_table <- iel[[]] |>
  transmute(
    celltype = factor(
      as.character(celltype),
      levels = celltype_order
    ),
    source = factor(
      tolower(as.character(sample)),
      levels = c("iwt", "iko")
    ),
    HCDown_score = HCDown1
  ) |>
  filter(
    !is.na(celltype),
    !is.na(source)
  )

hcdown_statistics <- hcdown_table |>
  group_by(celltype) |>
  wilcox_test(
    HCDown_score ~ source,
    exact = FALSE
  ) |>
  adjust_pvalue(
    method = "BH"
  ) |>
  add_significance(
    p.col = "p.adj",
    output.col = "p.adj.signif"
  ) |>
  add_xy_position(
    x = "source",
    scales = "free_y"
  )

source_colors <- c(
  "iwt" = "#D55E00",
  "iko" = "#0072B2"
)

p_hcdown_score <- ggplot(
  hcdown_table,
  aes(
    x = source,
    y = HCDown_score,
    fill = source
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.65
  ) +
  facet_wrap(
    ~ celltype,
    ncol = 4,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = source_colors
  ) +
  stat_pvalue_manual(
    hcdown_statistics,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    title = "HCDown module score",
    subtitle = "wt versus ko",
    x = NULL,
    y = "Module score"
  ) +
  theme_publication() +
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

print(p_hcdown_score)

save_table(
  hcdown_table,
  "HCDown_wt_vs_ko_scores.csv"
)

save_table(
  hcdown_statistics,
  "HCDown_wt_vs_ko_statistics.csv"
)

save_plot(
  p_hcdown_score,
  "HCDown_wt_vs_ko",
  width = 10,
  height = 8
)






#----------------------------deg and go-----------------------
iel$cluster_group <- paste(
  iel$celltype,
  iel$sample
)

Idents(iel) <- "cluster_group"

DimPlot(
  iel,
  reduction = "umap",
  group.by = "cluster_group"
)


# DEG筛选阈值
cut_off_FDR <- 0.1
cut_off_log2FC <- 0.25
iel<- JoinLayers(
  iel,
  assay = "RNA"
)

# =============================================================================
# Memory-like CD8aa+ TCRab/gd nIEL
# =============================================================================

Memory_CD8aa <- FindMarkers(
  iel,
  ident.1 = "Memory-like CD8aa+ TCRab/gd nIEL iko",
  ident.2 = "Memory-like CD8aa+ TCRab/gd nIEL iwt",
  verbose = FALSE
)

Memory_CD8aa$Sig <- ifelse(
  Memory_CD8aa$p_val_adj < cut_off_FDR &
    abs(Memory_CD8aa$avg_log2FC) >= cut_off_log2FC,
  ifelse(
    Memory_CD8aa$avg_log2FC > cut_off_log2FC,
    "Up",
    "Down"
  ),
  "no"
)

Memory_CD8aa$genename <- rownames(Memory_CD8aa)

write.csv(
  Memory_CD8aa,
  file = file.path(
    table_dir,
    "Memory_CD8aa_KO_vs_WT_DEG.csv"
  ),
  row.names = FALSE
)

table(Memory_CD8aa$Sig)


# =============================================================================
# Effector-like CD8aa+ TCRab/gd nIEL
# =============================================================================

Effector_CD8aa <- FindMarkers(
  iel,
  ident.1 = "Effector-like CD8aa+ TCRab/gd nIEL iko",
  ident.2 = "Effector-like CD8aa+ TCRab/gd nIEL iwt",
  verbose = FALSE
)

Effector_CD8aa$Sig <- ifelse(
  Effector_CD8aa$p_val_adj < cut_off_FDR &
    abs(Effector_CD8aa$avg_log2FC) >= cut_off_log2FC,
  ifelse(
    Effector_CD8aa$avg_log2FC > cut_off_log2FC,
    "Up",
    "Down"
  ),
  "no"
)

Effector_CD8aa$genename <- rownames(Effector_CD8aa)

write.csv(
  Effector_CD8aa,
  file = file.path(
    table_dir,
    "Effector_CD8aa_KO_vs_WT_DEG.csv"
  ),
  row.names = FALSE
)

table(Effector_CD8aa$Sig)


# =============================================================================
# Memory-like CD8ab+ T cells
# =============================================================================

CD8ab <- FindMarkers(
  iel,
  ident.1 = "CD8ab+ T cells iko",
  ident.2 = "CD8ab+ T cells iwt",
  verbose = FALSE
)

CD8ab$Sig <- ifelse(
  CD8ab$p_val_adj < cut_off_FDR &
    abs(CD8ab$avg_log2FC) >= cut_off_log2FC,
  ifelse(
    CD8ab$avg_log2FC > cut_off_log2FC,
    "Up",
    "Down"
  ),
  "no"
)

CD8ab$genename <- rownames(CD8ab)

write.csv(
  CD8ab,
  file = file.path(
    table_dir,
    "CD8ab_KO_vs_WT_DEG.csv"
  ),
  row.names = FALSE
)

table(CD8ab$Sig)




# =============================================================================
# CD4 T cells
# =============================================================================

CD4 <- FindMarkers(
  iel,
  ident.1 = "CD4 T cells iko",
  ident.2 = "CD4 T cells iwt",
  verbose = FALSE
)

CD4$Sig <- ifelse(
  CD4$p_val_adj < cut_off_FDR &
    abs(CD4$avg_log2FC) >= cut_off_log2FC,
  ifelse(
    CD4$avg_log2FC > cut_off_log2FC,
    "Up",
    "Down"
  ),
  "no"
)

CD4$genename <- rownames(CD4)

write.csv(
  CD4,
  file = file.path(
    table_dir,
    "CD4_KO_vs_WT_DEG.csv"
  ),
  row.names = FALSE
)

table(CD4$Sig)

save.image(file="20260917 deg.RData")

deg_list <- list(
  Memory_CD8aa   = Memory_CD8aa,
  Effector_CD8aa = Effector_CD8aa,
  CD8ab   = CD8ab,
  CD4= CD4
)


run_go <- function(
    deg_data,
    celltype_name,
    direction
) {
  
  go_genes <- unique(
    deg_data$genename[
      deg_data$Sig == direction
    ]
  )
  
  # 将进入DEG分析的基因作为背景
  background_genes <- unique(
    deg_data$genename
  )
  
  go_result <- enrichGO(
    gene = go_genes,
    universe = background_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  # 去除语义高度相似的重复GO条目
  go_result_simple <- simplify(
    go_result,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  go_table <- as.data.frame(go_result_simple) |>
    mutate(
      celltype = celltype_name,
      direction = direction
    )
  
  write.csv(
    go_table,
    file = file.path(
      table_dir,
      paste0(
        celltype_name,
        "_",
        direction,
        "_GO_BP_simplified.csv"
      )
    ),
    row.names = FALSE
  )
  
  go_table
}
# KO上调基因
Memory_CD8aa_GO_Up <- run_go(
  Memory_CD8aa, "Memory_CD8aa", "Up"
)

Effector_CD8aa_GO_Up <- run_go(
  Effector_CD8aa, "Effector_CD8aa", "Up"
)

Memory_CD8ab_GO_Up <- run_go(
  Memory_CD8ab, "Memory_CD8ab", "Up"
)

Effector_CD8ab_GO_Up <- run_go(
  Effector_CD8ab, "Effector_CD8ab", "Up"
)

CD4_GO_Up <- run_go(
  CD4, "CD4", "Up"
)


# KO下调基因
Memory_CD8aa_GO_Down <- run_go(
  Memory_CD8aa, "Memory_CD8aa", "Down"
)

Effector_CD8aa_GO_Down <- run_go(
  Effector_CD8aa, "Effector_CD8aa", "Down"
)

Memory_CD8ab_GO_Down <- run_go(
  Memory_CD8ab, "Memory_CD8ab", "Down"
)

Effector_CD8ab_GO_Down <- run_go(
  Effector_CD8ab, "Effector_CD8ab", "Down"
)

CD4_GO_Down <- run_go(
  CD4, "CD4", "Down"
)



go_up <- bind_rows(
  Memory_CD8aa_GO_Up,
  Effector_CD8aa_GO_Up,
  Memory_CD8ab_GO_Up,
  Effector_CD8ab_GO_Up,
  CD4_GO_Up
) |>
  filter(p.adjust < 0.05) |>
  group_by(celltype) |>
  slice_min(
    order_by = p.adjust,
    n = 15,
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    celltype = factor(
      celltype,
      levels = c(
        "Memory_CD8aa",
        "Effector_CD8aa",
        "Memory_CD8ab",
        "Effector_CD8ab",
        "CD4"
      )
    ),
    log10_FDR = -log10(p.adjust)
  )

write.csv(
  go_up,
  file.path(
    table_dir,
    "Major_celltypes_Up_top10_overall_GO_BP.csv"
  ),
  row.names = FALSE
)
order_go_by_similarity <- function(go_data) {
  
  # 每个GO条目保留唯一的ID和名称
  go_terms <- go_data |>
    distinct(ID, Description)
  
  # 只有一个GO条目时无需聚类
  if (nrow(go_terms) <= 1) {
    go_data$Description <- factor(
      go_data$Description,
      levels = unique(go_data$Description)
    )
    
    return(go_data)
  }
  
  # 计算GO条目之间的语义相似度
  go_similarity <- mgoSim(
    go_terms$ID,
    go_terms$ID,
    semData = mouse_go_data,
    measure = "Wang",
    combine = NULL
  )
  
  go_similarity[is.na(go_similarity)] <- 0
  diag(go_similarity) <- 1
  
  # 相似度转换成距离并进行层次聚类
  go_cluster <- hclust(
    as.dist(1 - go_similarity),
    method = "average"
  )
  
  ordered_ids <- go_terms$ID[
    go_cluster$order
  ]
  
  ordered_descriptions <- go_terms$Description[
    match(
      ordered_ids,
      go_terms$ID
    )
  ]
  
  # rev()使聚类顺序从图的上方向下显示
  go_data |>
    mutate(
      Description = factor(
        Description,
        levels = rev(unique(ordered_descriptions))
      )
    )
}

library(GOSemSim)

mouse_go_data <- godata(
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  computeIC = TRUE
)

go_up <- order_go_by_similarity(go_up)

p_go_up <- ggplot(
  go_up,
  aes(
    x = celltype,
    y = Description
  )
) +
  geom_point(
    aes(
      size = Count,
      color = log10_FDR
    )
  ) +
  scale_color_gradient(
    low = "#6BAEDB",
    high = "#D64065",
    name = expression(-log[10]("FDR"))
  ) +
  scale_size_continuous(
    range = c(3, 9),
    name = "Gene count"
  ) +
  scale_x_discrete(
    labels = c(
      "Memory_CD8aa" = "Memory CD8aa",
      "Effector_CD8aa" = "Effector CD8aa",
      "Memory_CD8ab" = "Memory CD8ab",
      "Effector_CD8ab" = "Effector CD8ab",
      "CD4" = "CD4 T cells"
    )
  ) +
  labs(
    title = "Overall GO enrichment of genes upregulated in KO",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    )
  )

print(p_go_up)

save_plot(
  p_go_up,
  "Major_celltypes_Up_top10_overall_GO_BP",
  width = 8,
  height = 10
)


go_down_all <- bind_rows(
  Memory_CD8aa_GO_Down,
  Effector_CD8aa_GO_Down,
  Memory_CD8ab_GO_Down,
  Effector_CD8ab_GO_Down,
  CD4_GO_Down
)


# ---- 每群选择FDR最小的前10个GO条目 -----------------------------------------

go_down <- go_down_all |>
  filter(
    p.adjust < 0.05
  ) |>
  group_by(celltype) |>
  slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    celltype = factor(
      celltype,
      levels = c(
        "Memory_CD8aa",
        "Effector_CD8aa",
        "Memory_CD8ab",
        "Effector_CD8ab",
        "CD4"
      )
    ),
    log10_FDR = -log10(p.adjust)
  )


# ---- 根据GO语义相似性排列条目 -----------------------------------------------

go_down <- order_go_by_similarity(
  go_down
)


# ---- 保存用于绘图的GO结果 ----------------------------------------------------

write.csv(
  go_down,
  file = file.path(
    table_dir,
    "Major_celltypes_Down_top10_semantic_order_GO_BP.csv"
  ),
  row.names = FALSE
)


# ---- 绘制GO气泡图 ------------------------------------------------------------

p_go_down <- ggplot(
  go_down,
  aes(
    x = celltype,
    y = Description
  )
) +
  geom_point(
    aes(
      size = Count,
      color = log10_FDR
    )
  ) +
  scale_color_gradient(
    low = "#6BAEDB",
    high = "#D64065",
    name = expression(-log[10]("FDR"))
  ) +
  scale_size_continuous(
    range = c(3, 9),
    name = "Gene count"
  ) +
  scale_x_discrete(
    labels = c(
      "Memory_CD8aa" = "Memory CD8aa",
      "Effector_CD8aa" = "Effector CD8aa",
      "Memory_CD8ab" = "Memory CD8ab",
      "Effector_CD8ab" = "Effector CD8ab",
      "CD4" = "CD4 T cells"
    )
  ) +
  labs(
    title = "Overall GO enrichment of genes downregulated in KO",
    subtitle = "GO terms ordered by semantic similarity",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )

print(p_go_down)


# ---- 保存图片 ---------------------------------------------------------------

save_plot(
  p_go_down,
  "Major_celltypes_Down_top10_semantic_order_GO_BP",
  width = 12,
  height = 10
)

#-----------------------------vocanol---------------------------
cut_off_FDR <- 0.05
cut_off_log2FC <- 0.25

volcano_colors <- c(
  "Up" = "#D64065",
  "Down" = "#3977AF",
  "no" = "grey80"
)


Memory_CD8aa$genelist <- NA

# 早期活化
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Jun")] <- "Jun"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Junb")] <- "Junb"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Jund")] <- "Jund"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Fos")] <- "Fos"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Fosb")] <- "Fosb"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Cd69")] <- "Cd69"

# T-cell信号和免疫活化
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Cd6")] <- "Cd6"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Lat")] <- "Lat"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Zap70")] <- "Zap70"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Ctla4")] <- "Ctla4"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Dtx1")] <- "Dtx1"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Nfkbid")] <- "Nfkbid"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Thy1")] <- "Thy1"

# 组蛋白/染色质及细胞结构
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "H4c18")] <- "H4c18"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "H4c2")] <- "H4c2"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "H2aj")] <- "H2aj"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Tuba1a")] <- "Tuba1a"

# 代表性protein-folding基因
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Hspa1a")] <- "Hspa1a"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Hsph1")] <- "Hsph1"

# 下调：代谢、转运及蛋白合成
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Slc16a11")] <- "Slc16a11"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Slc14a1")] <- "Slc14a1"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Slc7a1")] <- "Slc7a1"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Enox2")] <- "Enox2"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Impact")] <- "Impact"
Memory_CD8aa$genelist[which(Memory_CD8aa$genename == "Rpl41")] <- "Rpl41"

Effector_CD8aa$genelist <- NA

# 早期活化
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Jun")] <- "Jun"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Junb")] <- "Junb"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Jund")] <- "Jund"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Fos")] <- "Fos"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Fosb")] <- "Fosb"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Cd69")] <- "Cd69"

# T-cell活化、细胞因子和迁移
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Il2")] <- "Il2"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Il7r")] <- "Il7r"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Il27ra")] <- "Il27ra"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Ccr5")] <- "Ccr5"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Cd6")] <- "Cd6"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Lat")] <- "Lat"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Ctla4")] <- "Ctla4"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Tnfaip3")] <- "Tnfaip3"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Nfkbid")] <- "Nfkbid"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Lgals9")] <- "Lgals9"

# 组蛋白/染色质
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "H4c18")] <- "H4c18"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "H2aj")] <- "H2aj"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Tuba1a")] <- "Tuba1a"

# 代表性protein-folding基因
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Hspa1a")] <- "Hspa1a"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Hsph1")] <- "Hsph1"

# 下调：代谢及跨膜转运
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Slc16a11")] <- "Slc16a11"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Slc14a1")] <- "Slc14a1"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Slc7a1")] <- "Slc7a1"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Slc22a15")] <- "Slc22a15"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Slc3a2")] <- "Slc3a2"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Acaca")] <- "Acaca"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Cyp4v3")] <- "Cyp4v3"
Effector_CD8aa$genelist[which(Effector_CD8aa$genename == "Enox2")] <- "Enox2"


Memory_CD8ab$genelist <- NA

# 早期活化
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Jun")] <- "Jun"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Junb")] <- "Junb"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Jund")] <- "Jund"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Fos")] <- "Fos"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Fosb")] <- "Fosb"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Nr4a2")] <- "Nr4a2"

# TCR信号和免疫活化
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Zap70")] <- "Zap70"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Lat")] <- "Lat"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Ptpn6")] <- "Ptpn6"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Rhoh")] <- "Rhoh"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Cd6")] <- "Cd6"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Thy1")] <- "Thy1"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Tnfaip3")] <- "Tnfaip3"

# 细胞增殖/细胞结构
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Tuba1a")] <- "Tuba1a"

# 代表性protein-folding基因
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Hspa1a")] <- "Hspa1a"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Hsph1")] <- "Hsph1"

# 下调：氨基酸转运、代谢和蛋白合成
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Slc38a2")] <- "Slc38a2"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Slc9a9")] <- "Slc9a9"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Impact")] <- "Impact"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Cpq")] <- "Cpq"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Rps21")] <- "Rps21"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Rplp2")] <- "Rplp2"
Memory_CD8ab$genelist[which(Memory_CD8ab$genename == "Rpl37a")] <- "Rpl37a"


Effector_CD8ab$genelist <- NA

# 早期活化
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Jun")] <- "Jun"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Junb")] <- "Junb"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Jund")] <- "Jund"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Fos")] <- "Fos"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Fosb")] <- "Fosb"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Cd69")] <- "Cd69"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Atf4")] <- "Atf4"

# T-cell相关
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Cd3d")] <- "Cd3d"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Cd8a")] <- "Cd8a"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Cd8b1")] <- "Cd8b1"

# 组蛋白
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "H2aj")] <- "H2aj"

# 代表性protein-folding基因
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Hspa1a")] <- "Hspa1a"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Hspa1b")] <- "Hspa1b"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Hsph1")] <- "Hsph1"

# 下调：代谢、转运及蛋白合成
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Acoxl")] <- "Acoxl"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Slc3a2")] <- "Slc3a2"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Rpl41")] <- "Rpl41"
Effector_CD8ab$genelist[which(Effector_CD8ab$genename == "Rps2")] <- "Rps2"

CD4$genelist <- NA

# 早期活化
CD4$genelist[which(CD4$genename == "Jun")] <- "Jun"
CD4$genelist[which(CD4$genename == "Junb")] <- "Junb"
CD4$genelist[which(CD4$genename == "Jund")] <- "Jund"
CD4$genelist[which(CD4$genename == "Fos")] <- "Fos"
CD4$genelist[which(CD4$genename == "Fosb")] <- "Fosb"

# T-cell activation和免疫调控
CD4$genelist[which(CD4$genename == "Cd3d")] <- "Cd3d"
CD4$genelist[which(CD4$genename == "Cd6")] <- "Cd6"
CD4$genelist[which(CD4$genename == "Ctla4")] <- "Ctla4"
CD4$genelist[which(CD4$genename == "Thy1")] <- "Thy1"
CD4$genelist[which(CD4$genename == "Tnfaip3")] <- "Tnfaip3"
CD4$genelist[which(CD4$genename == "Ccl5")] <- "Ccl5"

# 代表性protein-folding基因
CD4$genelist[which(CD4$genename == "Hspa1a")] <- "Hspa1a"
CD4$genelist[which(CD4$genename == "Hsph1")] <- "Hsph1"

# 下调：氨基酸转运、翻译和代谢调节
CD4$genelist[which(CD4$genename == "Slc38a1")] <- "Slc38a1"
CD4$genelist[which(CD4$genename == "Slc38a2")] <- "Slc38a2"
CD4$genelist[which(CD4$genename == "Impact")] <- "Impact"
CD4$genelist[which(CD4$genename == "Rpl39l")] <- "Rpl39l"
CD4$genelist[which(CD4$genename == "Rps2")] <- "Rps2"
CD4$genelist[which(CD4$genename == "Rpl41")] <- "Rpl41"







#----------------------------------------------------------
plot_manual_volcano <- function(
    sig,
    plot_title,
    file_name
) {
  
  sig <- sig |>
    mutate(
      Sig = factor(
        Sig,
        levels = c("no", "Down", "Up")
      ),
      minus_log10_FDR = -log10(
        pmax(
          p_val_adj,
          .Machine$double.xmin
        )
      )
    )
  
  # 上调和下调标签分开
  label_up <- sig |>
    filter(
      !is.na(genelist),
      avg_log2FC > 0
    )
  
  label_down <- sig |>
    filter(
      !is.na(genelist),
      avg_log2FC < 0
    )
  
  # 计算x轴范围，为标签预留空间
  x_limit <- max(
    abs(sig$avg_log2FC),
    na.rm = TRUE
  )
  
  p <- ggplot(
    sig,
    aes(
      x = avg_log2FC,
      y = minus_log10_FDR
    )
  ) +
    
    # 所有基因点；不再给标注点增加黑色描边
    geom_point(
      aes(color = Sig),
      size = 1.5,
      alpha = 0.72
    ) +
    
    # fold change阈值线
    geom_vline(
      xintercept = c(
        -cut_off_log2FC,
        cut_off_log2FC
      ),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey55"
    ) +
    
    # FDR阈值线
    geom_hline(
      yintercept = -log10(cut_off_FDR),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey55"
    ) +
    
    # -------------------------------------------------------------------------
  # 上调基因标签
  # 白色背景能够防止连接线穿过文字
  # -------------------------------------------------------------------------
  geom_label_repel(
    data = label_up,
    aes(label = genelist),
    direction = "y",
    nudge_x = 0.35,
    hjust = 0,
    size = 3.1,
    color = "#9E1B32",
    
    # 标签背景：白色、无边框
    fill = scales::alpha(
      "white",
      0.90
    ),
    label.size = 0,
    label.padding = unit(
      0.10,
      "lines"
    ),
    
    # 标签排斥设置
    box.padding = 0.45,
    point.padding = 0.25,
    force = 2,
    force_pull = 0.4,
    max.overlaps = Inf,
    max.time = 5,
    
    # 连接线：更细、更浅
    min.segment.length = 0.05,
    segment.color = scales::alpha(
      "#D64065",
      0.45
    ),
    segment.linewidth = 0.25,
    segment.curvature = 0.08,
    segment.ncp = 3,
    segment.angle = 20,
    
    seed = 20260729
  ) +
    
    # -------------------------------------------------------------------------
  # 下调基因标签
  # -------------------------------------------------------------------------
  geom_label_repel(
    data = label_down,
    aes(label = genelist),
    direction = "y",
    nudge_x = -0.35,
    hjust = 1,
    size = 3.1,
    color = "#174A73",
    
    # 标签背景：白色、无边框
    fill = scales::alpha(
      "white",
      0.90
    ),
    label.size = 0,
    label.padding = unit(
      0.10,
      "lines"
    ),
    
    # 标签排斥设置
    box.padding = 0.45,
    point.padding = 0.25,
    force = 2,
    force_pull = 0.4,
    max.overlaps = Inf,
    max.time = 5,
    
    # 连接线：更细、更浅
    min.segment.length = 0.05,
    segment.color = scales::alpha(
      "#3977AF",
      0.45
    ),
    segment.linewidth = 0.25,
    segment.curvature = -0.08,
    segment.ncp = 3,
    segment.angle = 20,
    
    seed = 20260729
  ) +
    
    scale_color_manual(
      values = c(
        "Up" = "#D64065",
        "Down" = "#3977AF",
        "no" = "grey82"
      ),
      breaks = c(
        "Up",
        "Down",
        "no"
      ),
      labels = c(
        "Up in KO",
        "Down in KO",
        "Not significant"
      ),
      drop = FALSE
    ) +
    
    scale_x_continuous(
      limits = c(
        -x_limit * 1.28,
        x_limit * 1.28
      ),
      expand = expansion(
        mult = c(0.02, 0.02)
      )
    ) +
    
    labs(
      title = plot_title,
      subtitle = "KO versus WT",
      x = expression(log[2]("fold change (KO/WT)")),
      y = expression(-log[10]("adjusted P value")),
      color = NULL
    ) +
    
    theme_publication() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.position = "top",
      plot.margin = margin(
        t = 8,
        r = 20,
        b = 8,
        l = 20
      )
    )
  
  print(p)
  
  save_plot(
    p,
    file_name,
    width = 8,
    height = 7
  )
  
  p
}
p_volcano_Memory_CD8aa <- plot_manual_volcano(
  Memory_CD8aa,
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Volcano_Memory_CD8aa_KO_vs_WT"
)

p_volcano_Effector_CD8aa <- plot_manual_volcano(
  Effector_CD8aa,
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "Volcano_Effector_CD8aa_KO_vs_WT"
)

p_volcano_Memory_CD8ab <- plot_manual_volcano(
  Memory_CD8ab,
  "Memory-like CD8ab+ T cells",
  "Volcano_Memory_CD8ab_KO_vs_WT"
)

p_volcano_Effector_CD8ab <- plot_manual_volcano(
  Effector_CD8ab,
  "Effector-like CD8ab+ T cells",
  "Volcano_Effector_CD8ab_KO_vs_WT"
)

p_volcano_CD4 <- plot_manual_volcano(
  CD4,
  "CD4 T cells",
  "Volcano_CD4_KO_vs_WT"
)



#------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

celltype_order <- c(
  "Memory-like CD8aa+ TCRab/gd nIEL",
  "Effector-like CD8aa+ TCRab/gd nIEL",
  "CD8ab+ T cells",
  "CD4 T cells",
  "Naive T cells",
  "Th17 cells",
  "Treg cells",
  "Cycling T cells"
)

marker_genes <- c(
  "Trac",
  "Trdc",
  "Cd8a",
  "Cd4",
  "Cd8b1",
  "Tcf7",
  "Gzmb",
  "S1pr1",
  "Mki67",
  "Foxp3",
  "Il17a"
)

p_feature_umap <- FeaturePlot(
  object = iel,
  features = marker_genes,
  reduction = "umap",
  
  cols = c(
    "#D3D3D3",
    "#004C99"
  ),
  
  # 强化阳性细胞颜色
  min.cutoff = 0,
  max.cutoff = "q50",
  
  order = TRUE,
  pt.size = 0.6,
  
  ncol = 4,
  raster = TRUE,
  raster.dpi = c(300, 300),
  combine = TRUE
)
p_feature_umap <- p_feature_umap &
  theme_publication(
    base_size = 9,
    base_family = "sans"
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(
      face = "italic",
      size = 10,
      hjust = 0.5,
      colour = "black"
    ),
    
    legend.position = "none"
  )
p_feature_umap

save_plot(
  plot = p_feature_umap,
  filename = "UMAP_T_cell_marker_FeaturePlot",
  width = 10,
  height = 7.5,
  dpi = 600,
  bg = "white"
)


gene_colors <- c(
  "Trac"    = "#4C78A8",  # 蓝色
  "Tcrg-C1" = "#72A98F",  # 灰绿色
  "Cd8a"    = "#D98C95",  # 粉色
  "Cd4"     = "#E9C95B",  # 黄色
  "Cd8b1"   = "#C9365B",  # 玫红
  "Tcf7"    = "#F27658",  # 珊瑚橙
  "Gzmb"    = "#70B7C4",  # 青蓝
  "S1pr1"   = "#B6A6A0",  # 暖灰
  "Mki67"   = "#E99662",  # 浅橙
  "Il17a"   = "#8E6BBE",  # 紫色
  "Foxp3"   = "#58A65C"   # 绿色
)

stacked_violin_plot <- function(
    object,
    features,
    celltype_col = "celltype",
    celltype_order,
    assay = "RNA"
) {
  DefaultAssay(object) <- assay
  
  # 检查基因是否存在
  missing_features <- setdiff(features, rownames(object))
  
  if (length(missing_features) > 0) {
    stop(
      "以下基因未在当前 assay 中找到：",
      paste(missing_features, collapse = ", ")
    )
  }
  
  plot_data <- FetchData(
    object = object,
    vars = c(celltype_col, features)
  )
  
  plot_data$celltype <- factor(
    plot_data[[celltype_col]],
    levels = celltype_order
  )
  
  plot_data <- plot_data |>
    filter(!is.na(celltype)) |>
    select(celltype, all_of(features)) |>
    pivot_longer(
      cols = all_of(features),
      names_to = "feature",
      values_to = "expression"
    )
  
  # 控制基因从上至下的顺序
  plot_data$feature <- factor(
    plot_data$feature,
    levels = features
  )
  
  ggplot(
    plot_data,
    aes(
      x = celltype,
      y = expression,
      fill = feature,
      group = interaction(feature, celltype)
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "#D8D8D8",
      linewidth = 0.4
    ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      width = 0.86,
      adjust = 1,
      colour = "#555555",
      linewidth = 0.3
    ) +
    facet_grid(
      rows = vars(feature),
      scales = "free_y",
      switch = "y"
    ) +
    scale_x_discrete(
      limits = celltype_order,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = gene_colors,
      drop = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # 细胞群对应的纵向参考线
      panel.grid.major.x = element_line(
        colour = "#E5E5E5",
        linewidth = 0.4
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      
      axis.text.x = element_text(
        colour = "#333333",
        angle = 48,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
      axis.text.y = element_blank(),
      
      # 在左侧显示基因名称
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        colour = "#333333",
        angle = 0,
        hjust = 1,
        face = "italic",
        size = 9
      ),
      
      panel.spacing.y = grid::unit(0, "pt"),
      legend.position = "none",
      
      plot.margin = margin(
        t = 15,
        r = 15,
        b = 15,
        l = 15
      )
    )
}

p_stacked_violin <- stacked_violin_plot(
  object = iel,
  features = marker_genes,
  celltype_col = "celltype",
  celltype_order = celltype_order,
  assay = "RNA"
)
p_stacked_violin 
save_plot(
  plot = p_stacked_violin,
  filename = "stacked_violin_T_cell_markers",
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)

