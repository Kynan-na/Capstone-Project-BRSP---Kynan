library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)
library(hgu133a.db)

library(clusterProfiler)
library(org.Hs.eg.db) # Database anotasi untuk manusia (Homo sapiens)
library(enrichplot)
library(ggplot2)

# Unduh dataset microarray dari GEO
gset <- getGEO("GSE122063", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Ekstraksi Data Ekspresi
ex <- exprs(gset)

# Cek dan aplikasikan log2 transformasi jika diperlukan
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) 
}
exprs(gset) <- ex

# Ekstraksi dan pengelompokan metadata (phenotype)
pData_df <- pData(gset)
groups <- rep("Other", nrow(pData_df))

# Beri label klasifikasi berdasarkan kondisi pasien di metadata
groups[grepl("Alzheimer", pData_df$characteristics_ch1, ignore.case = TRUE)] <- "Alzheimer"
groups[grepl("Control", pData_df$characteristics_ch1, ignore.case = TRUE)] <- "Control"

# Filter data: hanya gunakan sampel Alzheimer dan Control
valid_samples <- groups %in% c("Alzheimer", "Control")
gset_filtered <- gset[, valid_samples]
groups_filtered <- groups[valid_samples]

# Buat factor dan design matrix
f <- factor(groups_filtered, levels = c("Control", "Alzheimer"))
design <- model.matrix(~ 0 + f)
colnames(design) <- levels(f)

# Buat matriks kontras untuk perbandingan Alzheimer vs Control
contrast.matrix <- makeContrasts(Alzheimer_vs_Control = Alzheimer - Control, levels = design)

# Jalankan model statistik (limma pipeline)
fit <- lmFit(gset_filtered, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ==========================================
# Visualisasi Hasil Analisis
# ==========================================

# 1. Volcano Plot
# Ambil seluruh populasi gen untuk plotting global
top_genes <- topTable(fit2, adjust="fdr", sort.by="P", number=Inf)
ex_filtered <- exprs(gset_filtered)

# Klasifikasi gen (Up-regulated, Down-regulated, Not Significant)
top_genes$Expression <- "NO"
top_genes$Expression[top_genes$logFC > 1 & top_genes$adj.P.Val < 0.01] <- "Up-regulated"
top_genes$Expression[top_genes$logFC < -1 & top_genes$adj.P.Val < 0.01] <- "Down-regulated"

volcano_plot <- ggplot(top_genes, aes(x = logFC, y = -log10(P.Value), color = Expression)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs Alzheimer vs Control", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

print(volcano_plot)

# 2. Heatmap
# Ambil ID probe dari 50 gen paling signifikan
top_genes_50 <- topTable(fit2, adjust="fdr", sort.by="P", number=50)
top_gene_ids <- rownames(top_genes_50)

# Potong matriks ekspresi khusus untuk 50 gen teratas
heatmap_data <- ex_filtered[top_gene_ids, ]

# Setup anotasi kolom heatmap
annotation_col <- data.frame(
  Diagnosis = factor(groups_filtered, levels = c("Control", "Alzheimer"))
)
rownames(annotation_col) <- colnames(heatmap_data)
ann_colors <- list(
  Diagnosis = c(Control = "navy", Alzheimer = "firebrick")
)

pheatmap(heatmap_data, 
         scale = "row",                   # Normalisasi nilai Z-score per baris
         annotation_col = annotation_col, # Tambahkan pita grup di atas heatmap
         annotation_colors = ann_colors,  # Terapkan custom warna
         show_rownames = FALSE,           # Sembunyikan ID gen
         show_colnames = FALSE,           # Sembunyikan nama sampel
         cluster_cols = TRUE,             # Clustering kolom (sampel)
         cluster_rows = TRUE,             # Clustering baris (gen)
         main = "Heatmap of Top 50 DEGs (Alzheimer vs Control)")


# ==========================================
# Anotasi Probe ke Gene Symbol
# ==========================================

# Gabungkan data fitur anotasi langsung ke objek fit2
fit2$genes <- fData(gset_filtered)

# Ambil hasil lengkap dengan kolom anotasi bawaan
all_degs <- topTable(fit2, adjust="fdr", sort.by="P", number=Inf)

# Standarisasi penamaan kolom
colnames(all_degs)[colnames(all_degs) == "GENE_SYMBOL"] <- "Symbol"
colnames(all_degs)[colnames(all_degs) == "LOCUSLINK_ID"] <- "Entrez_ID"

# Pembersihan Data Anotasi
# a. Hapus baris tanpa Gene Symbol
annotated_degs <- all_degs[all_degs$Symbol != "" & !is.na(all_degs$Symbol), ]

# b. Hapus duplikasi Symbol, simpan yang memiliki p-value terkecil
annotated_degs <- annotated_degs[order(annotated_degs$adj.P.Val), ]
final_degs <- annotated_degs[!duplicated(annotated_degs$Symbol), ]


# ==========================================
# Functional Enrichment Analysis (GO & KEGG)
# ==========================================

# Ekstrak Entrez ID untuk gen target (signifikan) dan universe (background)
sig_genes <- subset(final_degs, adj.P.Val < 0.05 & abs(logFC) > 1)
sig_genes_clean <- sig_genes[!is.na(sig_genes$Entrez_ID), ]
sig_entrez <- as.character(sig_genes_clean$Entrez_ID[!is.na(sig_genes$Entrez_ID)])
universe_entrez <- as.character(final_degs$Entrez_ID[!is.na(final_degs$Entrez_ID)])

# 1. Eksekusi Analisis Gene Ontology (BP)
go_results <- enrichGO(gene          = sig_entrez,
                       universe      = universe_entrez, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",            
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

dotplot_go <- dotplot(go_results, showCategory = 15, title = "Top 15 Enriched GO Terms (Biological Process)")
print(dotplot_go)

# 2. Eksekusi Analisis KEGG Pathway
kegg_results <- enrichKEGG(gene         = sig_entrez,
                           organism     = 'hsa', 
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01)

# Konversi Entrez ID ke Symbol untuk label visualisasi
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

barplot_kegg <- barplot(kegg_results, showCategory = 8, title = "Top 8 Enriched KEGG Pathways")
print(barplot_kegg)

# 3. Visualisasi Network (Cnetplot)
# Siapkan vektor fold_changes berformat numerik murni
fold_changes <- as.numeric(sig_genes$logFC) 
names(fold_changes) <- sig_genes$Symbol
fold_changes <- fold_changes[!is.na(fold_changes)]

# Render Cnetplot
go_cnetplot <- cnetplot(go_results, 
                        foldChange = fold_changes, 
                        showCategory = 3)

print(go_cnetplot)