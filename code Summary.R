####################################################GF IEC analysis

rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(tidyverse)
library(stringr)
library(DESeq2)
library(sva)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(edgeR)
library(msigdbr)
library(enrichplot)
library(fgsea)
library(stringr)
library(org.Mm.eg.db)
library(pheatmap)
library(ggvenn)
library(RDAVIDWebService)



data <- read.csv("D spf-gf.csv", check.names = FALSE)
names(data)[1] <- "ENG"
colSums(data[,c(2:25)])

Mm <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl",
                       host = "https://www.ensembl.org")
geneinfo <- getBM(attributes = c("ensembl_gene_id","transcript_length","mgi_symbol"), mart = Mm)

data <- as.data.frame(data)
rownames(data) <- data$ENG
data$symbol <- geneinfo[match(data$ENG, geneinfo$ensembl_gene_id),3]
data <- data[which(rowSums(data[,c(2:25)]) > 0),]

###########

condition <- c(rep(c("SPF-C"), each = 3),rep(c("SPF-D"), each = 3),rep(c("SPF-J"), each = 3),rep(c("SPF-I"), each = 3),rep(c("GF-C"), each = 3),rep(c("GF-D"), each = 3),rep(c("GF-J"), each = 3),rep(c("GF-I"), each = 3))
colData <- data.frame(row.names = colnames(data[,c(2:25)]),condition)
summary(data[,2:25])
# 检查counts是否为整数（DESeq2要求输入为整数）

any(apply(data[,2:25], 2, function(x) !all(x == floor(x))))

integer_counts <- round(data[,2:25])

dds <- DESeqDataSetFromMatrix(countData =integer_counts, colData = colData, design = ~condition)
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = c("condition"), returnData = F)

dds <- DESeq(dds)

############
data_TPM <- data
data_TPM$transcript_length <- ''

test <- geneinfo
test <- test[!duplicated(geneinfo$ensembl_gene_id),]

for (i in 1:78317) {
  test[i,2] <- max(c(geneinfo[which(test[i,1] == geneinfo$ensembl_gene_id),'transcript_length']))
}

data_TPM$transcript_length <- test[match(data_TPM$ENG, test$ensembl_gene_id),2]

data_TPM <- na.omit(data_TPM)

#### TPM ####

tl <- as.numeric(data_TPM$transcript_length)
countdata <- data_TPM[,c(2:25)]

for (i in 1:24){
  n <- countdata[,i]/tl
  m <- sum(n)
  t <- n/m*1000000
  countdata[,i] <- t
}
data_TPM[,2:25] <- countdata
colSums(data_TPM[,c(2:25)])




########


########

DEG_up <- function(x){x$ENG <- rownames(x)
y <- x[which(x$padj <= 0.05 & x$log2FoldChange >= 1),]
z <- merge(y,data[,c(1,26)],by.x = "ENG",by.y = "ENG")
z <- z[order(z$log2FoldChange, decreasing = T),]
return(z)
}

DEG_down <- function(x){
  x$ENG <- rownames(x)
  y <- x[which(x$padj <= 0.05 & x$log2FoldChange <=  -1),]
  z <- merge(y,data[,c(1,26)],by.x = "ENG",by.y = "ENG")
  z <- z[order(z$log2FoldChange, decreasing = F),]
  return(z)
}


########### 
Duo <-  data.frame(results(dds, contrast = c("condition","GF-D","SPF-D")))
up_D <- DEG_up(x)
down_D <- DEG_down(x)

Jej <-  data.frame(results(dds, contrast = c("condition","GF-J","SPF-J")))
up_J <- DEG_up(x)
down_J <- DEG_down(x)

Ile <-  data.frame(results(dds, contrast = c("condition","GF-I","SPF-I")))
up_I <- DEG_up(x)
down_I <- DEG_down(x)

Colon <-  data.frame(results(dds, contrast = c("condition","GF-C","SPF-C")))
up_C <- DEG_up(x)
down_C <- DEG_down(x)
############################################GO analysis was done online by DAVID
goDup<-read.delim("david go d up.txt")
goJup<-read.delim("david go j up.txt")
goIup<-read.delim("david go i up.txt")
goCup<-read.delim("david go c up.txt")


goDup<- goDup %>% mutate(segment = "Duodeunum")  # 添加来源标识

goJup<- goJup %>% mutate(segment = "Jejunum") 
goIup<- goIup %>% mutate(segment = "Ileum")
goCup<- goCup %>% mutate(segment = "Colon")


go.D<-goDup[which(goDup$FDR<0.05),]
go.J<-goJup[which(goJup$FDR<0.05),]
go.I<-goIup[which(goIup$FDR<0.05),]
go.C<-goCup[which(goCup$FDR<0.05),]



top10go_D <- head(go.D[order(go.D$FDR), ], 3)
top10go_J <- head(go.J[order(go.J$FDR), ], 3)
top10go_I <- head(go.I[order(go.I$FDR), ], 3)
top10go_C <- head(go.C[order(go.C$FDR), ], 3)


go_up <- bind_rows(top10go_D, top10go_J, top10go_I, top10go_C) 


plot<- ggplot(go_up,aes(x=factor(segment, level=c('Duodeunum','Jejunum','Ileum','Colon')),Term))+
  geom_point(aes(size=Count,color=-log10(FDR)))+
  scale_color_gradient(low = "#6baedb",high = "#d64065")+theme_bw() +theme(axis.text.y=element_text(size=10, color = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size=10))+scale_size_continuous(range=c(3,9))
plot
#######################################################################################
Duo<- Duo %>% mutate(segment = "Duodeunum")  # 添加来源标识
Jej<- Jej %>% mutate(segment = "Jejunum") 
Ile<- Ile %>% mutate(segment = "Ileum")
Colon<- Colon %>% mutate(segment = "Colon")
dt <- bind_rows(Duo, Jej, Ile, Colon) 




Mm <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl",
                       host = "https://www.ensembl.org")
geneinfo <- getBM(attributes = c("ensembl_gene_id","transcript_length","mgi_symbol"), mart = Mm)


dt$symbol <- geneinfo[match(data$ENG, geneinfo$ensembl_gene_id),3]



up_D<- up_D %>% mutate(segment = "Duodeunum")  # 添加来源标识
up_J<- up_J %>% mutate(segment = "Jejunum") 
up_I<- up_I %>% mutate(segment = "Ileum")
up_C<- up_C %>% mutate(segment = "Colon")

down_D<- down_D %>% mutate(segment = "Duodeunum")  # 添加来源标识
down_J<- down_J %>% mutate(segment = "Jejunum") 
down_I<- down_I %>% mutate(segment = "Ileum")
down_C<- down_C %>% mutate(segment = "Colon")


up_D<- up_D %>% mutate(direction = "up")  # 添加来源标识
up_J<- up_J %>% mutate(direction = "up") 
up_I<- up_I %>% mutate(direction = "up")
up_C<- up_C %>% mutate(direction = "up")

down_D<- down_D %>% mutate(direction = "down")  # 添加来源标识
down_J<- down_J %>% mutate(direction = "down") 
down_I<- down_I %>% mutate(direction = "down")
down_C<- down_C %>% mutate(direction = "down")

sig<- rbind(up_D, up_J,up_I,up_C,down_D,down_J,down_I,down_C)


sig$segment <- factor(sig$segment, levels = unique(sig$segment))

sig$direction <- factor(sig$direction, levels = unique(sig$direction))



sig$genelist<-NA
sig$genelist[which(sig$symbol=="Dbp")]<-"Dbp"
sig$genelist[which(sig$symbol=="Tef")]<-"Tef"
sig$genelist[which(sig$symbol=="Per1")]<-"Per1"
sig$genelist[which(sig$symbol=="Per3")]<-"Per3"
sig$genelist[which(sig$symbol=="Per2")]<-"Per2"
sig$genelist[which(sig$symbol=="Ciart")]<-"Ciart"
sig$genelist[which(sig$symbol=="Nr1d1")]<-"Nr1d1"
sig$genelist[which(sig$symbol=="Nr1d2")]<-"Nr1d2"
sig$genelist[which(sig$symbol=="Rorc")]<-"Rorc"
sig$genelist[which(sig$symbol=="Cry2")]<-"Cry2"
sig$genelist[which(sig$symbol=="Hlf")]<-"Hlf"
sig$genelist[which(sig$symbol=="Cipc")]<-"Cipc"

sig$genelist[which(sig$symbol=="Muc2")]<-"Muc2"
sig$genelist[which(sig$symbol=="Muc4")]<-"Muc4"
sig$genelist[which(sig$symbol=="Muc6")]<-"Muc6"
sig$genelist[which(sig$symbol=="Muc13")]<-"Muc13"
sig$genelist[which(sig$symbol=="Muc20")]<-"Muc20"
sig$genelist[which(sig$symbol=="Reg3a")]<-"Reg3a"
sig$genelist[which(sig$symbol=="Reg3b")]<-"Reg3b"
sig$genelist[which(sig$symbol=="Reg3g")]<-"Reg3g"
sig$genelist[which(sig$symbol=="Lgr5")]<-"Lgr5"
sig$genelist[which(sig$symbol=="Ang4")]<-"Ang4"

#创建自定义色板：
mycol<- c("#7FBBC9", "#9A8CC3", "#96BC95", "#cda5a1", "#d64065","#517FFF")

p1
p<- ggplot() +
  geom_point(data = dt,
             aes(x = log2FoldChange, y = -log10(padj)),colour = 'grey',size=1) +
  coord_flip()

p1<- p + facet_grid(. ~ segment) #一行多列;如果需要一列多行则var ~ .

p2<- p1 + 
  geom_point(data = sig,
             aes(x = log2FoldChange, y = -log10(padj),
                 color= segment))
p3



#参数美化：
p3<- p2 +
  geom_vline(xintercept = c(-1, 1), 
             linewidth= 0.5, color = "grey50", lty = 'dashed')+ #添加阈值线
  scale_color_manual(values = mycol)+ #更改配色
  theme_bw()+
  theme(legend.position = 'none', #去掉图例
        panel.grid = element_blank(), #去掉背景网格
        axis.text = element_text(size = 10), #坐标轴标签大小
        axis.text.x = element_text(angle = 45, vjust = 0.8), #x轴标签旋转
        strip.text.x = element_text(size = 10, face = 'bold') #加粗分面标题
  )

p3 +geom_text_repel(data= sig,
                    aes(log2FoldChange, y = -log10(padj), 
                        label= genelist,color = direction),
                    fontface= 'italic',segment.size=0.4,
                    size= 7,nudge_y=120,direction="y",point.padding=0,box.padding=2,force = 5,segment.linetype = 3,max.overlaps=Inf,seed = 222,segment.alpha=0.5)+xlim(-10,10)+ylim(0,150)

################################################################################bulk seq analysis for 27-HC treated Jurkat cells
library(Seurat)
library(scRepertoire)
library(tidyverse)
library(patchwork)
library(harmony)
library(immunarch)
library(data.table)
library(patchwork)
library(anndata)
library(SeuratDisk)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(ComplexHeatmap)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(stringr)
library(reshape2)
library(data.table)
library(paletteer)
library(ggcorrplot)
library(RColorBrewer)
library(scRepertoire)
library(future)
library(glmGamPoi)
library(presto)
library(msigdbr)
library(DESeq2)
library(KEGGREST)
library(patchwork)

filename <- list.files(path = "BULK_LKX_250417/Data/")
filename <- str_extract(string = filename,pattern = ".*(?=Reads)")
filepath <- list.files(path = "BULK_LKX_250417/Data/",full.names = T)
for (i in 1:length(filepath)) {
  if (i == 1) {
    temp <- fread(filepath[i],verbose = T,skip = 4,header = F)
    temp <- temp[,1:2]
    counts <- temp
  }else{
    temp <- fread(filepath[i],verbose = T,skip = 4,header = F)
    counts <- bind_cols(counts,temp[,2])
  }
}
colnames(counts) <- c("Gene",filename)
counts <- counts[which(!duplicated(counts$Gene)),]
counts <- as.data.frame(counts)
rownames(counts) <- counts$Gene
counts <- counts[,-1]

load("GeneinfoENSG.RData")

coldata <- as.data.frame(colnames(counts))
coldata$Group <- rep(c("27HC","Ctrl"),each = 3)
colnames(coldata)[1] <- "names"
rownames(coldata) <- coldata$names

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Group)
dds <- DESeq(dds)
Res <- results(dds,tidy = T,contrast = c("Group","27HC","Ctrl"))
Resdata <- Res[which(Res$log2FoldChange != "NA"),]
Resdata <- Resdata[order(Resdata$log2FoldChange,decreasing = T),]
gseadata <- Resdata$log2FoldChange
names(gseadata) <- Resdata$row

#SGMSgenes <- c("SGMS1","SGMS2","CERS1","CERS2","CERS3","CERS4","CERS5","CERS6","SPTLC1","SPTLC2","DEGS1","KDSR")
#SGMSgenes <- read.gmt("GMTs/CholAcquisition_YCS.gmt")[,2]
SGMSgenes <- read.gmt("GMTs/KEGG_SPHINGOLIPID_METABOLISM.v2025.1.Hs.gmt")[,2]
SGMSgenes <- geneinfo$ensembl_gene_id[which(geneinfo$hgnc_symbol %in% SGMSgenes)]
SGMSgenes <- SGMSgenes[which(!duplicated(SGMSgenes))]
gmtgene <- as.data.frame(matrix(nrow = length(SGMSgenes),ncol = 2))
colnames(gmtgene) <- c("term","gene")
gmtgene$term <- "Sphingomylin Synthesis Associate Genes"
gmtgene$gene <- SGMSgenes

gseares <- GSEA(gseadata,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
P1 <- gseaplot2(gseares,c(1:length(gseares@result$ID)),
                pvalue_table = T,color = "black",base_size = 13,
                subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2))

SGMSgenes <- read.gmt("GMTs/CholAcquisition_YCS.gmt")[,2]
SGMSgenes <- geneinfo$ensembl_gene_id[which(geneinfo$hgnc_symbol %in% SGMSgenes)]
SGMSgenes <- SGMSgenes[which(!duplicated(SGMSgenes))]
gmtgene <- as.data.frame(matrix(nrow = length(SGMSgenes),ncol = 2))
colnames(gmtgene) <- c("term","gene")
gmtgene$term <- "Cholesterol Acquisition Genes"
gmtgene$gene <- SGMSgenes

gseares <- GSEA(gseadata,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
P2 <- gseaplot2(gseares,c(1:length(gseares@result$ID)),
                pvalue_table = T,color = "black",base_size = 13,
                subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2))

pdf("LKX_Bulk_GSEA_20250829.pdf",width = 7,height = 5)
P1
P2
dev.off()
############################################################single cell analysis
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
  theme(panel.grid =element_blank(),panel.border = element_blank(),axis.line = element_blank(),legend.key = element_blank(),legend.background = element_blank())+
  scale_fill_manual(values=paried)


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


################################sgms1 level analysis

library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(ggpubr)

objWT<-subset(x = objFIN, subset = orig.ident == "iwt")
objWT3<-subset(x = objWT, idents=c('CD8aa+ TCRab/gd T cells','CD8 T cells',"CD4 T cells"))

expr_sgms1 <- FetchData(objWT3, vars = c("Sgms1", "newcelltype")) # 使用标准化后的表达数据

comparisons <- list(
  c("CD8aa+ TCRab/gd T cells", "CD8 T cells"), 
  c("CD8aa+ TCRab/gd T cells", "CD4 T cells"),c("CD8 T cells","CD4 T cells"))

celllist <- c("CD8aa+ TCRab/gd T cells", "CD4 T cells","CD8 T cells")
for (i in 1:length(celllist)) {
  expr_sgms1_niel <- expr_sgms1[which(expr_sgms1$newcelltype == celllist[i]),]
  expr_sgms1_niel <- expr_sgms1_niel[order(expr_sgms1_niel$Sgms1,decreasing = T),]
  expr_sgms1_niel <- expr_sgms1_niel[1:round(nrow(expr_sgms1_niel)/4),]
  assign(paste("FQ",celllist[i],sep = "_"),expr_sgms1_niel)
}
new_exp <- rbind(`FQ_CD4 T cells`,`FQ_CD8 T cells`,`FQ_CD8aa+ TCRab/gd T cells`)

length(expr_sgms1_niel[which(expr_sgms1_niel$Sgms1 == 0),])
paried=c("#4f86c6","#5eb887","#d64065","#f8845c")
ggplot(new_exp,aes(`newcelltype`,`Sgms1`)) +
  geom_violin(aes(fill = `newcelltype`)) +
  scale_fill_manual(values=paried)+
  ylab("Sgms1") + xlab("cell type") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif",  size = 5  ) +
  geom_boxplot(width = 0.1,position = position_dodge(0.9))+theme_bw() +
  theme(
    plot.title = element_text(size = 20,hjust = 0.5),  # 标题字号
    axis.text.x = element_text(            # x轴文本设置
      hjust = 0.5, vjust = 0.5,            # 居中显示
      colour = "black", size = 10,         # 黑色/14号字
      face = 'bold'                        # 加粗
    ), 
    axis.text.y = element_text(size = 12), # y轴文本字号
    axis.title.y = element_text(           # y轴标题
      size = 14, face = 'bold'             # 14号/加粗
    ), 
    panel.border = element_blank(),        # 移除面板边框
    axis.line = element_line(              # 坐标轴线设置
      colour = "black", size = 1           # 黑色/粗细为1
    ),
    legend.text = element_text(            # 图例文本
      colour = "black", size = 12          # 黑色/14号字
    ),
    legend.title = element_text(           # 图例标题
      colour = "black", size = 12,         # 黑色/14号字
      face = 'bold'                        # 加粗
    ),
    panel.grid.major = element_blank(),   # 移除主网格线
    panel.grid.minor = element_blank()     # 移除次网格线
    ,legend.position = "none") +
  # 设置坐标轴标签：y轴标签文字，x轴留空
  ylab("Expression Level") + xlab("") +
  # 限制y轴范围(0-0.5)
  ylim(1, 3) +ggtitle("Sgms1")
# 添加统计检验：组间比较


SMsyn<-c("Sgms1","Sgms2","Cers1","Cers2","Cers3","Cers4","Cers5","Cers6","Sptlc1","Sptlc2","Degs1","Kdsr")
DotPlot(objWT3,features = SMtotal)+RotatedAxis()

aqesmsyn <- AddModuleScore(objWT3, features =list(SMsyn),name = "SMsynthesis")
VlnPlot(aqesmsyn,feature="SMsynthesis1",pt.size=0,cols = paried)+ylim(0,3)+geom_boxplot()
paried=c("#4f86c6","#5eb887","#d64065","#f8845c")


meta_data <- aqesmsyn@meta.data
new_table <- data.frame(newcelltype = meta_data$newcelltype, SMsynthesis1 = meta_data$SMsynthesis1)

ggplot(new_table ,aes(newcelltype,SMsynthesis1)) +
  geom_violin(aes(fill = `newcelltype`)) +
  scale_fill_manual(values=paried)+
  ylab("SMsynthesis") + xlab("cell type") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif",  size = 5  ) +
  geom_boxplot(width = 0.1,position = position_dodge(0.9))+theme_bw() +
  theme(
    plot.title = element_text(size = 20,hjust = 0.5),  # 标题字号
    axis.text.x = element_text(            # x轴文本设置
      hjust = 0.5, vjust = 0.5,            # 居中显示
      colour = "black", size = 10,         # 黑色/14号字
      face = 'bold'                        # 加粗
    ), 
    axis.text.y = element_text(size = 12), # y轴文本字号
    axis.title.y = element_text(           # y轴标题
      size = 14, face = 'bold'             # 14号/加粗
    ), 
    panel.border = element_blank(),        # 移除面板边框
    axis.line = element_line(              # 坐标轴线设置
      colour = "black", size = 1           # 黑色/粗细为1
    ),
    legend.text = element_text(            # 图例文本
      colour = "black", size = 12          # 黑色/14号字
    ),
    legend.title = element_text(           # 图例标题
      colour = "black", size = 12,         # 黑色/14号字
      face = 'bold'                        # 加粗
    ),
    panel.grid.major = element_blank(),   # 移除主网格线
    panel.grid.minor = element_blank()     # 移除次网格线
    ,legend.position = "none") +
  # 设置坐标轴标签：y轴标签文字，x轴留空
  ylab("Expression Level") + xlab("") +ggtitle("Syn")
# 限制y轴范围(0-0.5)