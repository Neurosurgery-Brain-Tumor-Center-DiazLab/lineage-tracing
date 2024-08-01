library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(dplyr)
library(plyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(escape)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dittoSeq)
library(sleepwalk)
library(DESeq2)
library(hexbin)
library(harmony)
library(BPCells)
library(DoubletFinder)
library(copykat)
library(infercnv)
library(scCustomize)
library(CellChat)
library(SeuratWrappers)
library(scATOMIC)
library(DiagrammeRsvg)
library(rsvg)
library(MAST)
library(openxlsx)
library(HGNChelper)
library(SCEVAN)
library(Rclusterpp)
library(cowplot)
library(ggbeeswarm)
library(Nebulosa)
library(ggrepel)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Azimuth)#infercnvとどこかで競合
library(future.apply)
library(scGSVA)
library(SCPA)
library(msigdbr)
library(stats)
library(multcomp)
library(broom)
library(ggpubr)
library(rstatix)
library(ggsignif)
library(ggtree) #LinRace入れた際に入っているのではないかと。
library(treeio) #LinRace入れた際に入っているのではないかと。
library(ggtreeExtra)
library(packcircles)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(ReactomePA)
library(biomaRt)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(monocle3)#もしかしたらharmonyと競合
library(CytoTRACE2)
library(slingshot)

library(TraceQC)

# #TraceQC, CRISPR barcodeの解析に必要
# #install.packages("knitr")#knitrが壊れてたら再インストールの必要あり
# library(devtools)
# devtools::install_github("LiuzLab/TraceQC")

# #ggtree知らない間に入っていた。
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")

# #ggtreeExtra
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtreeExtra")


# #packcircles, circle packing plot作りに必要
# install.packages("packcircles")

# # org.Mm.eg.db
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Mm.eg.db")

# # Orthology.eg.db、Human to Mouseできるっぽい。
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Orthology.eg.db")

# # Monocle3、非線形な発展も捉えられるらしい、計算重い、公式ホームページBiocManagerで色々入れろとあるが既に入っている模様
# library(devtools)
# devtools::install_github('cole-trapnell-lab/monocle3')

# # CytoTRACE2
# library(devtools)
# devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")

# # slingshot、連続性がある場合強いらしい、計算軽い
# # if (!require("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # 
# # BiocManager::install("slingshot")
#こっちの方が公式
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("kstreet13/slingshot")




setwd("/Volumes/TRANSCEND/DAISYdata1")
cellranger.dir <- "/Volumes/TRANSCEND/DAISYdata1"

sample_name <- "DAISY0" #少ないので一つ一つ確認して設定。
sample_name <- "DAISY13"

#### 各SeuratObject > DoubletFinderまで#.rdsの方が.Rdata（object名など環境も保存するらしい）よりファイルサイズが小さいので変更 ####
#setwdを全部変えてください。

load.file <- function(sample_name, treatment) { #sample_name, tumorvshealty, treatment, ponsvsthalamus
  cellranger.dir <- "/Volumes/TRANSCEND/DAISYdata1"
  
  setwd_path <- file.path(cellranger.dir,sample_name)
  setwd(setwd_path)
  
  file_path1 <- file.path(cellranger.dir,sample_name,"filtered_feature_bc_matrix")
  #pathはfunction中で指定、./DAISY0/filterd_feature_bc_matrixを読み込む
  
  rm(SeuratObject)
  
  CountData <- Read10X(data.dir = file_path1)
  #SeuratObject <- CreateSeuratObject(counts = CountData, project = "Nazarians", min.cells = 3, min.features = 200)
  #min.cells = 0にしないとintegrateする際にずれるらしい。
  SeuratObject <- CreateSeuratObject(counts = CountData, project = sample_name, min.cells = 0, min.features = 200)
  
  SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^mt-") #人だと"^MT-"
  SeuratObject[["percent.ribo"]] <- PercentageFeatureSet(SeuratObject, pattern = "^Rp[sl]") #人だと"^RP[SL]"
  
  plot1 <- VlnPlot(SeuratObject, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), cols = "lightblue", ncol = 4)
  plot2 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "black")
  plot3 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black")
  plot4 <- plot2+plot3
  ggsave(paste0(sample_name,".PreVlnplot.png"), plot= plot1, units="cm", width=15, height=10, dpi=200)
  ggsave(paste0(sample_name,".PreScatter.png"), plot= plot4, units="cm", width=25, height=10, dpi=200)
  
  #SeuratObject <- subset(SeuratObject, subset = nCount_RNA < 150000 & 200 < nFeature_RNA & nFeature_RNA < 20000 & percent.mt < 40)#Nazarian snRNA
  SeuratObject <- subset(SeuratObject, subset = nCount_RNA < 20000 & 200 < nFeature_RNA & nFeature_RNA < 5000 & percent.mt < 10)
  
  set.seed(605)
  
  #metadataをfunction(1,2,3,4)で付与
  SeuratObject@meta.data$SampleID <- sample_name
  #SeuratObject@meta.data$TumorvsHealthy <- tumorvshealty
  SeuratObject@meta.data$Treatment <- treatment
  #SeuratObject@meta.data$PonsvsThalamus <- ponsvsthalamus
  
  #正規化、Scaling、PCA
  SeuratObject <- NormalizeData(SeuratObject)
  SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)
  SeuratObject <- ScaleData(SeuratObject)
  SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(object = SeuratObject))
  print(SeuratObject[["pca"]], dims = 1:30, nfeatures = 10)
  
  # #次元決定 PC40くらいまで攻めた方が細分化
  # SeuratObject <- JackStraw(SeuratObject, num.replicate = 200)
  # SeuratObject <- ScoreJackStraw(SeuratObject, dims = 1:20)
  # JackStrawPlot(SeuratObject, dims = 1:20)
  plot5 <- ElbowPlot(SeuratObject, ndims = 50)
  ggsave(paste0(sample_name,".Elbow.png"), plot= plot5, units="cm", width=10, height=10, dpi=200, bg = "white")
  
  # UMAPで次元削減, PC40で以下計算、resolutionの値上げるとCluster数増える
  SeuratObject <- FindNeighbors(SeuratObject, dims = 1:20)
  SeuratObject <- FindClusters(SeuratObject, resolution = 1)
  # SeuratObject@meta.data$seurat_clustersがmeta.dataに追加
  SeuratObject <- RunUMAP(SeuratObject, dims = 1:20)
  plot6 <- DimPlot(SeuratObject, reduction = "umap") + ggtitle(sample_name) #+ plot_annotation(title = sample_name)
  ggsave(paste0(sample_name,".umap.png"), plot= plot6, units="cm", width=15, height=15, dpi=200)
  
  table(SeuratObject@meta.data$seurat_clusters)
  
  #save (SeuratObject, file = paste0(sample_name,".Rdata"))
  #return (SeuratObject) 
  
  # # DoubletFinder (pK Identificaton, ground-truth : CITE-seqとかで明確にdoubletの場合のみ)
  # sweep.res.list_obj.rejoined <- paramSweep(obj.rejoined, PCs = 1:20, sct = FALSE)
  # gt.calls <- obj.rejoined@meta.data[rownames(sweep.res.list_obj.rejoined[[1]]), "GT"]
  # ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
  # sweep.stats_obj.rejoined <- summarizeSweep(sweep.res.list_obj.rejoined, GT = TRUE, GT.calls = gt.calls)
  # bcmvn_obj.rejoined <- find.pK(sweep.stats_obj.rejoined)
  
  # DoubletFinder (pK Identificaton, no ground-truth)
  sweep.res.list_SeuratObject <- paramSweep(SeuratObject, PCs = 1:40, sct = FALSE)
  sweep.stats_SeuratObject <- summarizeSweep(sweep.res.list_SeuratObject, GT = FALSE)
  bcmvn_SeuratObject <- find.pK(sweep.stats_SeuratObject)
  
  # pKmax valueを選択
  ggplot(bcmvn_SeuratObject, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()
  pKmax <- bcmvn_SeuratObject %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
  pKmax <- as.numeric(as.character(pKmax[[1]]))
  
  plot7 <- ggplot(bcmvn_SeuratObject, aes(pK, BCmetric, group=1)) + geom_point() + geom_line() + theme_classic()
  ggsave(paste0(sample_name,".DoubletFinder.png"), plot= plot7, units="cm", width=20, height=5, dpi=200)
  
  #SeuratObject@meta.data #でClusteringResultsの名前選択
  annotations <- SeuratObject@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- SeuratObject@meta.data$例えばClusteringResults
  nExp_poi <- round(0.075*nrow(SeuratObject@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  SeuratObject <- doubletFinder(SeuratObject, PCs = 1:40, pN = 0.25, pK = pKmax, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  #再利用の場合は以下だが再利用とは？
  #SeuratObject <- doubletFinder(SeuratObject, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  
  meta.data <- colnames(SeuratObject@meta.data) #metadataのカラム名全部取り出して
  DF.classification <- meta.data[grep("DF.classification", colnames(SeuratObject@meta.data))] #DF.classification含まれるのを抽出して、DF.classificationとして格納
  #"DF.classification"になったDF.classification_0.25_pKmax_nExpを新しいメタデータ名DoubletFinderで追加
  SeuratObject@meta.data$DoubletFinder <- SeuratObject@meta.data$"DF.classification"
  
  #いちいち0.25_0.09_nExpのvalue入れなくてもwhichでワイルドカード的な扱い
  #SeuratObject.rmdoublet <- subset(SeuratObject, cells=rownames(SeuratObject@meta.data)[which(SeuratObject@meta.data$DF.classification == "Singlet")])
  
  
  #Cell Cycle Tirosh et al, 2015からSeuratで標準装備,#The following features are not present in the object出るけど無視
  #SeuratObject <- CellCycleScoring(SeuratObject, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  SeuratObject <- CellCycleScoring(SeuratObject, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)#, set.ident = TRUE)
  
  setwd("/Volumes/TRANSCEND/DAISYdata1")
  saveRDS(SeuratObject, file = paste0(sample_name,".rds"))
  
  file_path2 <- file.path(cellranger.dir,paste0(sample_name,".rds"))
  #SeuratObject <- readRDS(file_path2)
  return (SeuratObject) 
}

DAISY0 <- load.file("DAISY0","non")
DAISY13 <- load.file("DAISY13","Dox")

DAISY13$Phase

View(SeuratObject)
SeuratObject@assays$RNA

DimPlot(SeuratObject, reduction = "umap", group.by = "DoubletFinder")


#Cell Cycle Tirosh et al, 2015からSeuratで標準装備,#The following features are not present in the object出るけど無視
SeuratObject <- CellCycleScoring(SeuratObject, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
DimPlot(SeuratObject, reduction = "umap")#, group.by = "DoubletFinder")
FeaturePlot(SeuratObject, feature=c("Ick","Cdk20","Cxcl10","Cd274","Sphk1"), ncol=6)








#### reload.Rdata関数としてDoubletFinder後のrdsをsample_nameがSeuratObject名になったファイルを作る ####
RereadRDS <- function(sample_name){
  rm(SeuratObject)
  cellranger.dir <- "/Volumes/TRANSCEND/DAISYdata1"
  file_path2 <- file.path(cellranger.dir,paste0(sample_name,".rds"))
  SeuratObject <- readRDS(file_path2)
  assign(sample_name, SeuratObject, envir = .GlobalEnv) #assignを使わないと.GlobalEnv（どこからでもアクセス可能)に名前反映されない。
  return(sample_name)
}

setwd("/Volumes/TRANSCEND/DAISYdata1")
cellranger.dir <- "/Volumes/TRANSCEND/DAISYdata1"

sample_name_list <- c("DAISY0","DAISY13")
lapply(sample_name_list, RereadRDS)
colnames(DAISY0)

sample_name_list <- c("DAISY0_BC","DAISY13_BC")
lapply(sample_name_list, RereadRDS)
RereadRDS("DAISY13_BC")








####SeuratObjectの10xbarcodeの取得####
# Seuratオブジェクトからセル名を取得
tenxbarcode_with1 <- colnames(DAISY0)
# "-1" を削除
tenxbarcode <- sub("-1$", "", tenxbarcode_with1)
# テキストファイルに保存
write.table(tenxbarcode, file="DAISY0.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#writeLines(tenxbarcode, con="DAISY0_unix.txt", sep="\n") # Unix形式の改行コードで書き込む、seqkit動かない問題には意味なし


# Seuratオブジェクトからセル名を取得,# "-1" を削除
tenxbarcode <- colnames(DAISY13) %>% sub("-1$", "", .) #. は直前のobjectを意味する
# テキストファイルに保存
write.table(tenxbarcode, file="DAISY13.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#writeLines(tenxbarcode, con="DAISY13_unix.txt", sep="\n") # Unix形式の改行コードで書き込む、seqkit動かない問題には意味なし

tenxbarcode <- readLines("DAISY0.txt")
head(tenxbarcode) # ファイルの最初の数行を確認する

####ここからlinuxとpythonを使ってsequence解析を行います。####



####histgramでクローン数表示とDAISYの変異を長さが短くなることで表示####

barcode0 <- read.delim("/Users/masahirookada/max0.txt", header=TRUE, sep="\t")
barcode13 <- read.delim("/Users/masahirookada/max13.txt", header=TRUE, sep="\t")
clones0 <- table(barcode0$static)
clones13 <- table(barcode13$static)
mean(as.numeric(clones0))
mean(as.numeric(clones13))
df0 <- data.frame(clone=names(clones0), count=as.numeric(clones0), dataset="barcode0")
df13 <- data.frame(clone=names(clones13), count=as.numeric(clones13), dataset="barcode13")
df <- rbind(df0, df13)
#View(df)
ggplot(df0, aes(x=count)) +
  geom_histogram(binwidth=5, fill="blue", color="black", alpha=0.7) +
  labs(title="Clonal size", x="Counts", y="Frequency") +
  theme_minimal()

ggplot(df0, aes(x=count)) +
  geom_histogram(aes(y=..density..), binwidth=1, fill="blue", color="black", alpha=0.7) +
  geom_density(color="red", size=0.3) +
  labs(title="Clonal size", x="Counts", y="Density") +
  theme_minimal()

ggplot(df, aes(x=count, fill=dataset, color=dataset)) +
  geom_histogram(aes(y=..density..), binwidth=1, position="identity", alpha=0.7) +
  #geom_density(size=1, alpha=0.7) +
  scale_color_manual(values=c("gray10","orangered1"))+
  scale_fill_manual(values=c("blue2","lightpink1"))+
  labs(title="Clonal size", x="Counts", y="Density") +
  coord_cartesian(xlim = c(1, 3300))+ #, ylim = c(0, 0.02))+ ##xlim = c(1, 3300)
  #scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

mutate0 <- nchar(barcode0$mutate)
mutate13 <- nchar(barcode13$mutate)
hist(mutate0)
ggplot(data = data.frame(length = mutate13), aes(x = length)) +
  geom_histogram(aes(y=..density..), binwidth = 5, fill = "brown", color = "black") + #"lightslateblue" "brown"
  labs(title = "DAISY mutation", x = "DAISY length [bp]", y = "Frequency")+
  scale_x_continuous(expand = c(0,0),breaks=seq(0,150,10),limits=c(0,150))+
  theme_classic()




####ggtreeを使ってtree図の可視化####
# .newick形式(3-2njmodified.py)または.treefile形式(iqtree2)のファイルを読み込む。（nhx形式ならファイル注釈、nodeごとにannotationもついてくる）
tree <- read.tree("/Users/masahirookada/daisyprocessdata/test33.newick")
tree <- read.tree("/Users/masahirookada/daisyprocessdata/BCoutputATGC13/BC33.nexus.treefile")

tree <- read.tree("/Users/masahirookada/daisyprocessdata/test11.newick")
tree <- read.tree("/Users/masahirookada/daisyprocessdata/BCoutputATGC13/BC33.nexus.treefile")

tree <- read.tree("/Users/masahirookada/test_trees.nwk")
tree <- read.tree("/Users/masahirookada/testnj_trees.nwk")
tree <- read.tree("/Users/masahirookada/testnj1_trees.nwk")
tree <- read.tree("/Users/masahirookada/testnj2_trees.nwk")

# ツリーの可視化
ggtree(tree, size=0.5, color="purple",layout="roundrect") + geom_tiplab()#+
  # geom_strip('t1', 't3', barsize=2, color='red', label="associated taxa", offset.text=.1)+
  # geom_cladelab(node=3, label="another clade")+
  # geom_hilight(node=1, fill="steelblue", alpha=.6) +
  # geom_hilight(node=5, fill="darkgreen", alpha=.6) +
  # geom_hilight(node=3, fill='darkgreen', type="rect")
ggtree(tree, size=0.5, color="purple",branch.length="none",layout="roundrect") + geom_tiplab()
ggtree(tree, size=0.5, color="purple", ladderize=FALSE) + geom_tiplab()

ggtree(tree, layout="equal_angle", branch.length = 'none')+ geom_tiplab()

ggtree(tree, size=0.1, layout = "circular")
ggtree(tree, size=0.5, color="purple", layout = "circular", branch.length='none') + geom_tiplab() +
  geom_rootpoint(color="black", size=1) #+
  #geom_nodepoint(color="grey", size=0.01, alpha=0.8) +
  #geom_tippoint(color="red",size=0.01, alpha=0.5)#+ geom_tiplab()

## module2で色付ける
metadata <- data.frame(tenxbarcode = sub("-1$", "", rownames(DAISY13_BC@meta.data)),
                       module1 = DAISY13_BC$module1,module2 = DAISY13_BC$module2)
leafdata <- data.frame(tenxbarcode = tree$tip.label, stringsAsFactors = FALSE)
leafdata <- merge(leafdata, metadata, by.x = "tenxbarcode", by.y = "tenxbarcode", all.x = TRUE)
leafdata <- leafdata %>% mutate(color = case_when(module2 == "MES" ~ "firebrick1",
                                                  module2 == "AC" ~ "aquamarine1",
                                                  module2 == "OPC" ~ "lightskyblue1",
                                                  module2 == "NPC" ~ "khaki1",
                                                  TRUE ~ "grey90"))
# ggtreeプロットにノード情報を追加して色分け
ggtree(tree, size=0.5, color="purple",layout="circular", branch.length = 'none') %<+% leafdata +
  geom_point(aes(color = color), size = 4) +  # ノードの色分け
  geom_tiplab(aes(label = label), size = 1, align = TRUE) +  # ノードラベルの表示
  scale_color_identity()

## CytoTRACE2_Scoreのvalueで色付ける
metadata <- data.frame(tenxbarcode = sub("-1$", "", rownames(DAISY13_BC@meta.data)),
                       module1 = DAISY13_BC$module1,module2 = DAISY13_BC$module2,
                       CytoTRACE2_Score = DAISY13_BC$CytoTRACE2_Score, Phase = DAISY13_BC$Phase)
leafdata <- data.frame(tenxbarcode = tree$tip.label, stringsAsFactors = FALSE)
leafdata <- merge(leafdata, metadata, by.x = "tenxbarcode", by.y = "tenxbarcode", all.x = TRUE)

## CytoTRACE2_Scoreの値に基づいて色を設定
pal <- viridis(n = 256, option = "B") #"A": magma,"B": inferno,"C": plasma,"D": viridis,"E": cividis
leafdata <- leafdata %>%
  mutate(color = scales::col_numeric(palette = pal, domain = range(CytoTRACE2_Score, na.rm = TRUE))(CytoTRACE2_Score))

## Phaseに基づいて色を設定
leafdata <- leafdata %>% mutate(color = case_when(Phase == "G1" ~ "grey20",
                                                  Phase == "S" ~ "orangered1",
                                                  Phase == "G2M" ~ "grey90"))

# ggtreeプロットにノード情報を追加して色分け
ggtree(tree, size=0.5, color="purple",layout="roundrect", branch.length = 'none') %<+% leafdata +
  geom_point(aes(color = color), size = 4) +  # ノードの色分け
  geom_tiplab(aes(label = label), size = 2, align = TRUE) +  # ノードラベルの表示
  scale_color_identity()


## CytoTRACE2_Scoreのvalueで色付ける、heatmapの表示場所が変
library(ggtreeExtra)
library(ggnewscale)

metadata <- data.frame(tenxbarcode = sub("-1$", "", rownames(DAISY13_BC@meta.data)),
                       module1 = DAISY13_BC$module1,module2 = DAISY13_BC$module2,
                       CytoTRACE2_Score = DAISY13_BC$CytoTRACE2_Score, Phase = DAISY13_BC$Phase)
pal <- viridis(n = 256, option = "B") #"A": magma,"B": inferno,"C": plasma,"D": viridis,"E": cividis
metadata <- metadata %>%
  mutate(CytoTRACE2_Score_color = scales::col_numeric(palette = pal, domain = range(CytoTRACE2_Score, na.rm = TRUE))(CytoTRACE2_Score)) %>%
  mutate(Phase_color = case_when(Phase == "G1" ~ "grey20",
                                 Phase == "S" ~ "orangered1",
                                 Phase == "G2M" ~ "grey90")) %>%
  mutate(module2_color = case_when(module2 == "MES" ~ "firebrick1",
                                   module2 == "AC" ~ "aquamarine1",
                                   module2 == "OPC" ~ "lightskyblue1",
                                   module2 == "NPC" ~ "khaki1",
                                   TRUE ~ "grey90"))
leafdata <- data.frame(tenxbarcode = tree$tip.label, stringsAsFactors = FALSE)
leafdata <- merge(leafdata, metadata, by.x = "tenxbarcode", by.y = "tenxbarcode", all.x = TRUE)

# ggtree(tree, size=0.5, color="purple",layout="roundrect", branch.length = 'none') %<+% leafdata +
#   geom_tiplab(aes(label = label), size = 2, align = TRUE) +# ノードラベルの表示
#   geom_rect(aes(xmin = 1, xmax = 2, ymin = y - 0.5, ymax = y + 0.5, fill = Phase_color), color = NA) +
#   geom_rect(aes(xmin = 2, xmax = 3, ymin = y - 0.5, ymax = y + 0.5, fill = module2_color), color = NA) +
#   geom_rect(aes(xmin = 3, xmax = 4, ymin = y - 0.5, ymax = y + 0.5, fill = CytoTRACE2_Score_color), color = NA) +
#   scale_fill_identity() +
#   theme(legend.position = "right") +
#   labs(fill = "Legend") +
#   guides(fill = guide_colorbar(barwidth = 10, barheight = 1)) +
#   scale_x_continuous(breaks = 1:3, labels = c("Phase", "module2", "CytoTRACE2_Score"))

# Prepare the data for gheatmap
heatmap_data <- leafdata %>%
  dplyr::select(tenxbarcode, CytoTRACE2_Score_color, Phase_color, module2_color) %>%
  pivot_longer(cols = c(CytoTRACE2_Score_color, Phase_color, module2_color), names_to = "variable", values_to = "color") %>%
  pivot_wider(names_from = variable, values_from = color)  %>% column_to_rownames('tenxbarcode')
# # Convert color columns to factors for the legend
# heatmap_data$Phase_color <- factor(heatmap_data$Phase_color, levels = c("grey20", "orangered1", "grey90"), labels = c("G1", "S", "G2M"))
# heatmap_data$module2_color <- factor(heatmap_data$module2_color, levels = c("firebrick1", "aquamarine1", "lightskyblue1", "khaki1"), labels = c("MES", "AC", "OPC", "NPC"))
# Basic tree plot
p <- ggtree(tree, size=0.5, color="purple", layout="roundrect", branch.length = 'none') %<+% leafdata +
  geom_tiplab(aes(label = label), size = 2, align = TRUE)
# Adding heatmap, offset treeとheatmapの間、colorは枠線、
p <- gheatmap(p, heatmap_data, offset = 1.25, width = 0.1, colnames = FALSE, colnames_angle = 90, font.size = 5, color = "black")+
  scale_fill_identity() 
p

# Legend作るためだけにComplexHeatmapを使う
library(ComplexHeatmap)
library(circlize)
#Legend for CytoTRACE2_Score
# labels <- format(seq(min(metadata$CytoTRACE2_Score), max(metadata$CytoTRACE2_Score), length.out = 4), digits = 1, nsmall = 1)#digits有効数字、nsmall小数点の桁
# # Round the labels to 1 decimal place、round 0で整数、round 1で小数点１桁
# rounded_labels <- round(seq(min(metadata$CytoTRACE2_Score), max(metadata$CytoTRACE2_Score), length.out = 4), 1)
# # Format the labels to ensure they have 1 decimal place
# labels <- format(rounded_labels, nsmall = 1)
#まとめると
labels <-format(round(seq(min(metadata$CytoTRACE2_Score), max(metadata$CytoTRACE2_Score), length.out = 4), 1), nsmall = 1)
legend_CytoTRACE2_Score <- Legend(at = seq(min(metadata$CytoTRACE2_Score), max(metadata$CytoTRACE2_Score), length.out = 4),
                                  labels = labels,
                                  labels_gp = gpar(fontsize = 10),
                                  legend_gp = gpar(fill = pal),
                                  col_fun = scales::col_numeric(palette = pal, domain = range(metadata$CytoTRACE2_Score, na.rm = TRUE)),
                                  title = "CytoTRACE2 Score")
#draw(legend_CytoTRACE2_Score)
#Legend for Phase
legend_Phase <- Legend(at = c("G1", "S", "G2M"),
                       legend_gp = gpar(fill = c("grey20", "orangered1", "grey90")),
                       title = "Phase")
#draw(legend_Phase)
#Legend for Module2
legend_Module2 <- Legend(at = c("MES", "AC", "OPC", "NPC"),
                         legend_gp = gpar(fill = c("firebrick1", "aquamarine1", "lightskyblue1", "khaki1")),
                         title = "Module2")
#draw(legend_Module2)
#Combine legends using packLegend、いちいちクリアするのが手間なら
combined_legend <- packLegend(legend_CytoTRACE2_Score, legend_Phase, legend_Module2, direction = "vertical")#"horizontal"
draw(combined_legend)



# #分岐点がなぜか色ついてしまう
# ggtree(tree, size=0.5, color="purple",layout="circular", branch.length = 'none') %<+% leafdata +
#   geom_point(aes(color = color), size = 4) +  # ノードの色分け
#   geom_tiplab(aes(label = label), size = 1, align = TRUE)+
#   scale_color_manual(values = c("firebrick1" = "firebrick1",
#                               "aquamarine1" = "aquamarine1",
#                               "lightskyblue1" = "lightskyblue1",
#                               "khaki1" = "khaki1",
#                               "grey90" = "grey90"),
#                    name = "Module",
#                    labels = c("MES", "AC", "OPC", "NPC", "Other"),
#                    guide = guide_legend(title.position = "top", ncol = 1,
#                                         override.aes = list(size = 4))) +
#   theme(legend.position = "right")  # 右側に凡例を表示





#TraceQCだが使えるか怪しい
input_file ="/Users/masahirookada/mutate13cut.txt"
ref_file="/Users/masahirookada/Desktop/DAISYdata1/FLDAISYtraceqc.txt"
output_file = "./aligned_reads.txt"

ref <- parse_ref_file(ref_file)
plot_construct(ref,chr_per_row=120,chr_size=3)

# #pythonの方がglobal_alignmentから文字状態行列を作るの簡単。
# sequence_alignment_for_10x(input_file = input_file, ref_file = ref_file, output_file = output_file,
#                            match = 2, mismatch= -2, gapopen = -6, gapextension = -0.1, ncores = 10)
# sequence_permutation(ref_file,out="./seq_permutate.txt")
# 
# input_file <- system.file("extdata/test_data/raw_sequence","carlin_example.bam",package="TraceQC")
# View(input_file)


# # fastaファイルを読み込む
# fasta <- readDNAStringSet("/Users/masahirookada/Desktop/DAISYdata1/240528_LH00132_0314_B22MWC5LT3/Common_SB28_DAISY_DAY0_sub1_A1_S120_L007_merged_TAGTTACGCCAAGCTTGAATTC.fasta")
# 
# 
# # # fastaファイルから配列を読み込まないとエラー
# # ref <- readDNAStringSet("/Users/masahirookada/Desktop/DAISYdata1/FLDAISY.txt")
# # seq <- readDNAStringSet("/Users/masahirookada/mutate0.txt")
# 
# # 普通のtxtファイルとして読み込み
# ref <- readLines("/Users/masahirookada/Desktop/DAISYdata1/FLDAISYtraceqc.txt")[1] # 1行目の配列
# #ref <- DNAString(ref)#色がオシャレになるだけかと
# 
# seq <- readLines("/Users/masahirookada/mutate13cut.txt") # 各行に1配列
# #seq <- lapply(seq, DNAString)#色がオシャレになるだけかと
# 
# # Needleman-Wunsch (グローバルアラインメント、全長が一致しなくても強制的に）
# global_alignments <- lapply(seq, function(x) {
#   pairwiseAlignment(ref, x, type = "global")
# })
# 
# # # Smith-Waterman (ローカルアラインメント、最も似ているところを抽出)
# # local_alignments <- lapply(seq, function(x) {
# #   pairwiseAlignment(ref, x, type = "local")
# # })
# 
# # アラインメント結果のスコアを表示
# global_scores <- sapply(global_alignments, score)
# # local_scores <- sapply(local_alignments, score)
# 
# # スコアをデータフレームにまとめる
# alignment_scores <- data.frame(Global = global_scores, Local = local_scores)
# 
# # スコアをファイルに書き出す
# write.csv(alignment_scores, "alignment_scores.csv", row.names = FALSE)
# 
# 
# library(reticulate)
# use_python("/Users/masahirookada/miniforge3/envs/r-reticulate/bin/python", required=T)
# repl_python()






#### Barcodeのメタデータを合体させる####
BC <- read.delim("/Users/masahirookada/daisyprocessdata/BCoutputATGC13/BCall.txt", header=TRUE, sep=",") #sep=","で読み込む
# 必要なカラムを抽出
#BC <- BC %>% select(tenx, BCnumber)#何かと競合しているっぽい
BC <- BC %>% dplyr::select(tenx, BCnumber)
# tenxカラムに-1を追加
BC$tenx <- paste0(BC$tenx, "-1")
rownames(BC) <- BC$tenx
BC$tenx <- NULL 
#AddMetaDataでAAACCAAAGCCACATA-1ごとにBCnumber追加できる。
DAISY13 <- AddMetaData(object = DAISY13, metadata = BC)

saveRDS(DAISY13, file = paste0("DAISY13_BC",".rds"))

# # BCnumberの文字列から"BC"を削除して数値に変換し、数値の順番通りにファクターを作成
# BCorder <- unique(as.numeric(gsub("BC", "", DAISY13$BCnumber, fixed = TRUE)))
# BCorder <- sort(BCorder)
# # BCnumberを数値の順番通りにファクターに変換
# DAISY13$BCnumber <- factor(DAISY13$BCnumber, levels = paste0("BC", as.character(BCorder)))


BC <- read.delim("/Users/masahirookada/daisyprocessdata/BCoutputATGC0/BCall.txt", header=TRUE, sep=",") #sep=","で読み込む
# 必要なカラムを抽出
#BC <- BC %>% select(tenx, BCnumber)#何かと競合しているっぽい
BC <- BC %>% dplyr::select(tenx, BCnumber)
# tenxカラムに-1を追加
BC$tenx <- paste0(BC$tenx, "-1")
rownames(BC) <- BC$tenx
BC$tenx <- NULL 
#AddMetaDataでAAACCAAAGCCACATA-1ごとにBCnumber追加できる。
DAISY0 <- AddMetaData(object = DAISY0, metadata = BC)

# # BCnumberの文字列から"BC"を削除して数値に変換し、数値の順番通りにファクターを作成
# BCorder <- unique(as.numeric(gsub("BC", "", DAISY13$BCnumber, fixed = TRUE)))
# BCorder <- sort(BCorder)
# # BCnumberを数値の順番通りにファクターに変換
# DAISY13$BCnumber <- factor(DAISY13$BCnumber, levels = paste0("BC", as.character(BCorder)))

saveRDS(DAISY0, file = paste0("DAISY0_BC",".rds"))













####Circle packing plot####
# クラスター情報の抽出
cluster_info <- as.data.frame(table(DAISY13$BCnumber))
colnames(cluster_info) <- c("cluster", "size")
# サークルのパッキング計算
packing <- circleProgressiveLayout(cluster_info$size, sizetype = 'area') #"radius" or "area"
# データにパッキング結果をマージ
cluster_info <- cbind(cluster_info, packing)
# サークルの頂点の数を変更
dat.gg <- circleLayoutVertices(packing, npoints = 12)

# ラベルのサイズを調整するための関数を定義
# # サークルのサイズが大きいものから小さいものへと変更する
# adjust_label_size <- function(size, max_size, min_size) {
#   scaled_size <- scales::rescale(size, to = c(min_size, max_size))
#   return(scaled_size)
# }
# # プロットの作成
# ggplot() +
#   geom_polygon(data = dat.gg, aes(x, y, group = id, fill = as.factor(id)), color = "black", size=1, alpha = 0.8) +
#   #geom_text(data = cluster_info, aes(x, y, label = cluster), size = 3, color="navyblue") +
#   geom_text(data = cluster_info, aes(x, y, label = cluster, size = size), color="navyblue",
#             size = adjust_label_size(cluster_info$size, max_size = 15, min_size = 0.5)) + # サイズを調整
#   scale_fill_manual(values = rainbow(length(unique(cluster_info$cluster)))) + # 凡例のタイトルを変更
#   theme_void() +
#   theme(legend.position = "none") +
#   #theme(legend.position = "right", text = element_text(size = 12, family = "Arial")) +
#   coord_equal()

ggplot() +
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill = as.factor(id)), color = "black", linewidth=0.8, alpha = 0.8) +
  #scale_fill_distiller(palette = "RdPu", direction = 1 ) + #連続値じゃないと無理っぽい。
  scale_fill_manual(values = rainbow(length(unique(cluster_info$cluster))))+
  geom_text(data = cluster_info, aes(x, y, label = cluster, size = size), color="navyblue")+
  scale_size_continuous(range = c(0.5,15)) + #BC0は(0.5,1)でBC13は(0.5,15)
  theme_void() +
  theme(legend.position = "none") +
  #theme(legend.position = "right", text = element_text(size = 12, family = "Arial")) +
  coord_equal()



####Seuratの解析####
DAISY13 <- SetIdent(DAISY13, value = "BCnumber")

table(DAISY13$BCnumber)
#BCnumberメタデータがNAの場合、'nonBC'に置き換える
DAISY13@meta.data$BCnumber[is.na(DAISY13@meta.data$BCnumber)] <- "nonBC"
View(DAISY13@meta.data)

DimPlot_scCustom(DAISY13, reduction = "umap", group.by = c("seurat_clusters","Phase","DoubletFinder"))#+NoLegend()
DimPlot_scCustom(DAISY13, reduction = "umap", group.by = c("BCnumber"))+NoLegend()+ggtitle("Barcode clusters")
Cluster_Highlight_Plot(DAISY13, reduction = "umap", cluster_name = c("BC1","BC2","BC3","BC4","BC5","BC6","BC7","BC8"), pt.size = 0.1,
                       highlight_color = inferno(10),
                       background_color = "gainsboro")

DF <- data.frame(Phase = DAISY13@meta.data$Phase, BCnumber = DAISY13@meta.data$BCnumber, seurat_clusters = DAISY13@meta.data$seurat_clusters)

# BCnumberごとにPhaseの割合を計算
phase_percent <- DF %>%
  group_by(BCnumber, Phase) %>%
  tally() %>%  # tally()を使用して行数を数える
  mutate(percent = n / sum(n) * 100)
# BCnumberごとに全部の和が100%になるか確認
phase_percent <- phase_percent %>%
  group_by(BCnumber) %>%
  mutate(total_percent = sum(percent))
# グラフ作成
BCarbitrary <- paste0("BC", 1:50)
phase_percent <- phase_percent %>% filter(BCnumber %in% BCarbitrary) %>% mutate(BCnumber = factor(BCnumber, levels = BCarbitrary))
ggplot(phase_percent, aes(x = BCnumber, y = percent, fill = Phase)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~BCnumber, scales = "free") +
  labs(title = "Cell Cycle per Barcode Clone", x = "BCnumber", y = "Percentage") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # x軸のラベルを回転させて見やすくする

# BCnumberごとにseurat_clustersの割合を計算
seurat_clusters_percent <- DF %>%
  group_by(BCnumber, seurat_clusters) %>%
  tally() %>%  # tally()を使用して行数を数える
  mutate(percent = n / sum(n) * 100)
# BCnumberごとに全部の和が100%になるか確認
seurat_clusters_percent <- seurat_clusters_percent %>%
  group_by(BCnumber) %>%
  mutate(total_percent = sum(percent))
# グラフ作成
BCarbitrary <- paste0("BC", 1:50)
seurat_clusters_percent <- seurat_clusters_percent %>% filter(BCnumber %in% BCarbitrary) %>% mutate(BCnumber = factor(BCnumber, levels = BCarbitrary))
ggplot(seurat_clusters_percent, aes(x = BCnumber, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~BCnumber, scales = "free") +
  labs(title = "seurat_clusters per Barcode clone", x = "BCnumber", y = "Percentage") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # x軸のラベルを回転させて見やすくする


####Pseudobulk analysis####
pseudo <- AggregateExpression(DAISY13, assays = "RNA", return.seurat = T, group.by = c("BCnumber"))
head(Cells(pseudo))
tail(Cells(pseudo))

Idents(pseudo) <- "BCnumber"

DEG_BC <- FindMarkers(object = pseudo, ident.1 = c(paste0("BC", 1:10)), 
                      ident.2 = c(paste0("BC", 11:446)), test.use = "DESeq2")
head(DEG_BC, n = 100)

# DEG_BC <- FindMarkers(object = DAISY13, ident.1 = c(paste0("BC", 1:5)), 
#                       ident.2 = c(paste0("BC", 100:446)), test.use = "MAST")

# 簡単なvolcanoプロットを作成
ggplot(DEG_BC, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.1) +
  theme_classic() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(data = specific_pathway1, shape = 21, cex = 5, fill = "orangered2", color = "black", stroke = 0.3) +
  # #特定の値
  # geom_text_repel(data = subset(female_male_corr, group == "male"), aes(label = pathway), 
  #                 box.padding = unit(0.1, "lines"),
  #                 point.padding = unit(0.3, "lines"), color= "navy",parse = TRUE, size = 0.7)+
  # #特定の値
  # geom_text_repel(data = subset(female_male_corr, female_corr > -0.3 & male_corr < -0.5), aes(label = pathway), 
  #                 box.padding = unit(0.1, "lines"),
  #                 point.padding = unit(0.3, "lines"), color= "navy",parse = TRUE, size = 0.7)+
  #名前で指定
  geom_text_repel(data = subset(female_male_corr, grepl("hedge|cell_cycle", ignore.case = T, pathway)), aes(label = pathway), 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"), color="grey20",parse = TRUE, size = 2)

# 簡単なGO解析
significant <-DEG_BC[DEG_BC$p_val_adj<0.001,]
up <- significant[significant$avg_log2FC>0.3,]
down <- significant[significant$avg_log2FC<-0.5,]

# is.DEG <- as.logical(up$FDR<0.001)
# DEG.names <- rownames(up)[is.DEG]
DEG.names <- rownames(up)
View(DEG.names)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library(ReactomePA)
#keytypes(org.Mm.eg.db)
DEG.names.entrez <- bitr(DEG.names, fromType="SYMBOL", toType="ENTREZID",OrgDb=org.Mm.eg.db)
View(DEG.names.entrez)
ego <- enrichGO(gene          = DEG.names.entrez$ENTREZID,
                keyType       = 'ENTREZID',
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #MF: Molecular Function CC: Cellular Component BP: Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego))
dotplot(ego, showCategory=15)#900*400

gene_IDs <- DEG.names.entrez$ENTREZID
#View(up)
#logFC <- up$logFC
#names(logFC) <- DEG.names.entrez$ENTREZID
eReactome <- enrichPathway(gene=gene_IDs,organism="mouse",pvalueCutoff=0.05, readable=T)
dotplot(eReactome, showCategory=10)
emapplot(eReactome)
cnetplot(eReactome, categorySize="pvalue",foldChange=logFC,showCategory = 8)




####MetaProgram解析、AddModuleScoreよりもEscapeの方が正しい気がする####
#An integrative model of cellular states, plasticity and genetics for glioblastoma Cell(2019)のTableS2からコピペ

#### Excelからコピペして、("apple","banana","cat")作る方法 ####
#上から順にMES2,MES1,AC,OPC,NPC1,NPC2
data <- "HILPDA	ADM	DDIT3	NDRG1	HERPUD1	DNAJB9	TRIB3	ENO2	AKAP12	SQSTM1	MT1X	ATF3	NAMPT	NRN1	SLC2A1	BNIP3	LGALS3	INSIG2	IGFBP3	PPP1R15A	VIM	PLOD2	GBE1	SLC2A3	FTL	WARS	ERO1L	XPOT	HSPA5	GDF15	ANXA2	EPAS1	LDHA	P4HA1	SERTAD1	PFKP	PGK1	EGLN3	SLC6A6	CA9	BNIP3L	RPL21	TRAM1	UFM1	ASNS	GOLT1B	ANGPTL4	SLC39A14	CDKN1A	HSPA9"
data <- "CHI3L1	ANXA2	ANXA1	CD44	VIM	MT2A	C1S	NAMPT	EFEMP1	C1R	SOD2	IFITM3	TIMP1	SPP1	A2M	S100A11	MT1X	S100A10	FN1	LGALS1	S100A16	CLIC1	MGST1	RCAN1	TAGLN2	NPC2	SERPING1	C8orf4	EMP1	APOE	CTSB	C3	LGALS3	MT1E	EMP3	SERPINA3	ACTN1	PRDX6	IGFBP7	SERPINE1	PLP2	MGP	CLIC4	GFPT2	GSN	NNMT	TUBA1C	GJA1	TNFRSF1A	WWTR1"
data <- "CST3	S100B	SLC1A3	HEPN1	HOPX	MT3	SPARCL1	MLC1	GFAP	FABP7	BCAN	PON2	METTL7B	SPARC	GATM	RAMP1	PMP2	AQP4	DBI	EDNRB	PTPRZ1	CLU	PMP22	ATP1A2	S100A16	HEY1	PCDHGC3	TTYH1	NDRG2	PRCP	ATP1B2	AGT	PLTP	GPM6B	F3	RAB31	PPAP2B	ANXA5	TSPAN7"
data <- "BCAN	PLP1	GPR17	FIBIN	LHFPL3	OLIG1	PSAT1	SCRG1	OMG	APOD	SIRT2	TNR	THY1	PHYHIPL	SOX2-OT	NKAIN4	LPPR1	PTPRZ1	VCAN	DBI	PMP2	CNP	TNS3	LIMA1	CA10	PCDHGC3	CNTN1	SCD5	P2RX7	CADM2	TTYH1	FGF12	TMEM206	NEU4	FXYD6	RNF13	RTKN	GPM6B	LMF1	ALCAM	PGRMC1	HRASLS	BCAS1	RAB31	PLLP	FABP5	NLGN3	SERINC5	EPB41L2	GPR37L1"
data <- "DLL3	DLL1	SOX4	TUBB3	HES6	TAGLN3	NEU4	MARCKSL1	CD24	STMN1	TCF12	BEX1	OLIG1	MAP2	FXYD6	PTPRS	MLLT11	NPPA	BCAN	MEST	ASCL1	BTG2	DCX	NXPH1	HN1	PFN2	SCG3	MYT1	CHD7	GPR56	TUBA1A	PCBP4	ETV1	SHD	TNR	AMOTL2	DBN1	HIP1	ABAT	ELAVL4	LMF1	GRIK2	SERINC5	TSPAN13	ELMO1	GLCCI1	SEZ6L	LRRN1	SEZ6	SOX11"
data <- "STMN2	CD24	RND3	HMP19	TUBB3	MIAT	DCX	NSG1	ELAVL4	MLLT11	DLX6-AS1	SOX11	NREP	FNBP1L	TAGLN3	STMN4	DLX5	SOX4	MAP1B	RBFOX2	IGFBPL1	STMN1	HN1	TMEM161B-AS1	DPYSL3	SEPT3	PKIA	ATP1B1	DYNC1I1	CD200	SNAP25	PAK3	NDRG4	KIF5A	UCHL1	ENO2	KIF5C	DDAH2	TUBB2A	LBH	LOC150568	TCF4	GNG3	NFIB	DPYSL5	CRABP1	DBN1	NFIX	CEP170	BLCAP"

# ペーストしたデータを行ごとに分割し、各行をリストに格納する
lines <- strsplit(data, "\n")[[1]]
# 各行をタブで分割してリストに格納する
words_list <- lapply(lines, function(line) unlist(strsplit(line, "\t")))
# 各単語を引用符で囲み、カンマで区切る
quoted_lines <- lapply(words_list, function(words) paste0('("', paste(words, collapse = '","'), '")'))
# 各行を改行で結合する
quoted_data <- paste(quoted_lines, collapse = "\n")
# 結果を表示する
cat(quoted_data)

####GeneListの作成####
feature_list_human <- list(
  MES2 = c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT","NRN1","SLC2A1","BNIP3","LGALS3","INSIG2","IGFBP3","PPP1R15A","VIM","PLOD2","GBE1","SLC2A3","FTL","WARS","ERO1L","XPOT","HSPA5","GDF15","ANXA2","EPAS1","LDHA","P4HA1","SERTAD1","PFKP","PGK1","EGLN3","SLC6A6","CA9","BNIP3L","RPL21","TRAM1","UFM1","ASNS","GOLT1B","ANGPTL4","SLC39A14","CDKN1A","HSPA9"),
  MES1 = c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2","IFITM3","TIMP1","SPP1","A2M","S100A11","MT1X","S100A10","FN1","LGALS1","S100A16","CLIC1","MGST1","RCAN1","TAGLN2","NPC2","SERPING1","C8orf4","EMP1","APOE","CTSB","C3","LGALS3","MT1E","EMP3","SERPINA3","ACTN1","PRDX6","IGFBP7","SERPINE1","PLP2","MGP","CLIC4","GFPT2","GSN","NNMT","TUBA1C","GJA1","TNFRSF1A","WWTR1"),
  AC = c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1","GFAP","FABP7","BCAN","PON2","METTL7B","SPARC","GATM","RAMP1","PMP2","AQP4","DBI","EDNRB","PTPRZ1","CLU","PMP22","ATP1A2","S100A16","HEY1","PCDHGC3","TTYH1","NDRG2","PRCP","ATP1B2","AGT","PLTP","GPM6B","F3","RAB31","PPAP2B","ANXA5","TSPAN7"),
  OPC = c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4","LPPR1","PTPRZ1","VCAN","DBI","PMP2","CNP","TNS3","LIMA1","CA10","PCDHGC3","CNTN1","SCD5","P2RX7","CADM2","TTYH1","FGF12","TMEM206","NEU4","FXYD6","RNF13","RTKN","GPM6B","LMF1","ALCAM","PGRMC1","HRASLS","BCAS1","RAB31","PLLP","FABP5","NLGN3","SERINC5","EPB41L2","GPR37L1"),
  NPC1 = c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST","ASCL1","BTG2","DCX","NXPH1","HN1","PFN2","SCG3","MYT1","CHD7","GPR56","TUBA1A","PCBP4","ETV1","SHD","TNR","AMOTL2","DBN1","HIP1","ABAT","ELAVL4","LMF1","GRIK2","SERINC5","TSPAN13","ELMO1","GLCCI1","SEZ6L","LRRN1","SEZ6","SOX11"),
  NPC2 = c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4","DLX5","SOX4","MAP1B","RBFOX2","IGFBPL1","STMN1","HN1","TMEM161B-AS1","DPYSL3","SEPT3","PKIA","ATP1B1","DYNC1I1","CD200","SNAP25","PAK3","NDRG4","KIF5A","UCHL1","ENO2","KIF5C","DDAH2","TUBB2A","LBH","LOC150568","TCF4","GNG3","NFIB","DPYSL5","CRABP1","DBN1","NFIX","CEP170","BLCAP")
)

# 最初だけ大文字、残り小文字になるらしいがこれはfake
# feature_list_mouse <- lapply(feature_list_fake_mouse, function(x) {
#   str_to_title(tolower(x))
# })

# # biomaRtを使ってヒトの遺伝子リストをマウスの遺伝子リストに変換する関数、hostが正しく設定されてないらしく動かない。
# convert_human_to_mouse <- function(human_genes) {
#   # Ensemblデータベースに接続
#   ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#   mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#   # ヒトの遺伝子名から対応するマウスの遺伝子名を取得
#   result <- getLDS(attributes = c("hgnc_symbol"), 
#                    filters = "hgnc_symbol", 
#                    values = human_genes, 
#                    mart = ensembl, 
#                    attributesL = c("mgi_symbol"), 
#                    martL = mouse, 
#                    uniqueRows = TRUE)
#   # 結果をリストに変換
#   mouse_genes <- setNames(result$MGI.symbol, result$HGNC.symbol)
#   # 入力リストの順序に従って出力リストを作成し、対応するマウスの遺伝子がない場合はNAを設定
#   output <- sapply(human_genes, function(x) ifelse(x %in% names(mouse_genes), mouse_genes[[x]], NA))
#   return(output)
# }
# 
# MES2 <- c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT","NRN1","SLC2A1","BNIP3","LGALS3","INSIG2","IGFBP3","PPP1R15A","VIM","PLOD2","GBE1","SLC2A3","FTL","WARS","ERO1L","XPOT","HSPA5","GDF15","ANXA2","EPAS1","LDHA","P4HA1","SERTAD1","PFKP","PGK1","EGLN3","SLC6A6","CA9","BNIP3L","RPL21","TRAM1","UFM1","ASNS","GOLT1B","ANGPTL4","SLC39A14","CDKN1A","HSPA9")
# mMES2 <- convert_human_to_mouse(MES2)

human2mouse <- function(genes){
  ensemblegenes <- mapIds(org.Hs.eg.db, genes, "ENTREZID","SYMBOL") #Human SYMBOL to ENTREZID
  mapped <- select(Orthology.eg.db, keys=ensemblegenes, columns="Mus.musculus",keytype="Homo.sapiens") #Human ENTREZID to Mouse ENTREZID
  mapped$MouseSymbol <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID") # Mouse ENTREZID to Mouse SYMBOL
  return(unname(unlist(mapped$MouseSymbol))) #勝手にNA消えてるっぽい。SYMBOLだけ抽出するにはunname, unlistが大事
}

MES2 <- c("HILPDA","ADM","DDIT3","NDRG1","HERPUD1","DNAJB9","TRIB3","ENO2","AKAP12","SQSTM1","MT1X","ATF3","NAMPT","NRN1","SLC2A1","BNIP3","LGALS3","INSIG2","IGFBP3","PPP1R15A","VIM","PLOD2","GBE1","SLC2A3","FTL","WARS","ERO1L","XPOT","HSPA5","GDF15","ANXA2","EPAS1","LDHA","P4HA1","SERTAD1","PFKP","PGK1","EGLN3","SLC6A6","CA9","BNIP3L","RPL21","TRAM1","UFM1","ASNS","GOLT1B","ANGPTL4","SLC39A14","CDKN1A","HSPA9")
MES1 <- c("CHI3L1","ANXA2","ANXA1","CD44","VIM","MT2A","C1S","NAMPT","EFEMP1","C1R","SOD2","IFITM3","TIMP1","SPP1","A2M","S100A11","MT1X","S100A10","FN1","LGALS1","S100A16","CLIC1","MGST1","RCAN1","TAGLN2","NPC2","SERPING1","C8orf4","EMP1","APOE","CTSB","C3","LGALS3","MT1E","EMP3","SERPINA3","ACTN1","PRDX6","IGFBP7","SERPINE1","PLP2","MGP","CLIC4","GFPT2","GSN","NNMT","TUBA1C","GJA1","TNFRSF1A","WWTR1")
AC <- c("CST3","S100B","SLC1A3","HEPN1","HOPX","MT3","SPARCL1","MLC1","GFAP","FABP7","BCAN","PON2","METTL7B","SPARC","GATM","RAMP1","PMP2","AQP4","DBI","EDNRB","PTPRZ1","CLU","PMP22","ATP1A2","S100A16","HEY1","PCDHGC3","TTYH1","NDRG2","PRCP","ATP1B2","AGT","PLTP","GPM6B","F3","RAB31","PPAP2B","ANXA5","TSPAN7")
OPC <- c("BCAN","PLP1","GPR17","FIBIN","LHFPL3","OLIG1","PSAT1","SCRG1","OMG","APOD","SIRT2","TNR","THY1","PHYHIPL","SOX2-OT","NKAIN4","LPPR1","PTPRZ1","VCAN","DBI","PMP2","CNP","TNS3","LIMA1","CA10","PCDHGC3","CNTN1","SCD5","P2RX7","CADM2","TTYH1","FGF12","TMEM206","NEU4","FXYD6","RNF13","RTKN","GPM6B","LMF1","ALCAM","PGRMC1","HRASLS","BCAS1","RAB31","PLLP","FABP5","NLGN3","SERINC5","EPB41L2","GPR37L1")
NPC1 <- c("DLL3","DLL1","SOX4","TUBB3","HES6","TAGLN3","NEU4","MARCKSL1","CD24","STMN1","TCF12","BEX1","OLIG1","MAP2","FXYD6","PTPRS","MLLT11","NPPA","BCAN","MEST","ASCL1","BTG2","DCX","NXPH1","HN1","PFN2","SCG3","MYT1","CHD7","GPR56","TUBA1A","PCBP4","ETV1","SHD","TNR","AMOTL2","DBN1","HIP1","ABAT","ELAVL4","LMF1","GRIK2","SERINC5","TSPAN13","ELMO1","GLCCI1","SEZ6L","LRRN1","SEZ6","SOX11")
NPC2 <- c("STMN2","CD24","RND3","HMP19","TUBB3","MIAT","DCX","NSG1","ELAVL4","MLLT11","DLX6-AS1","SOX11","NREP","FNBP1L","TAGLN3","STMN4","DLX5","SOX4","MAP1B","RBFOX2","IGFBPL1","STMN1","HN1","TMEM161B-AS1","DPYSL3","SEPT3","PKIA","ATP1B1","DYNC1I1","CD200","SNAP25","PAK3","NDRG4","KIF5A","UCHL1","ENO2","KIF5C","DDAH2","TUBB2A","LBH","LOC150568","TCF4","GNG3","NFIB","DPYSL5","CRABP1","DBN1","NFIX","CEP170","BLCAP")
mMES2 <- human2mouse(MES2)
mMES1 <- human2mouse(MES1)
mAC <- human2mouse(AC)
mOPC <- human2mouse(OPC)
mNPC1 <- human2mouse(NPC1)
mNPC2 <- human2mouse(NPC2)

feature_list_mouse <- list(MES2=mMES2, MES1=mMES1, AC=mAC, OPC=mOPC, NPC1=mNPC1, NPC2=mNPC2)

#AddModuleScore速いけど案外適当らしい
DAISY13 <- AddModuleScore(DAISY13, features=feature_list_mouse, name = c("MES2_","MES1_","AC_","OPC_","NPC1_","NPC2_"), assays = "RNA", ctrl = 10, seed = 613)
DAISY13$MES2_1
DAISY13$MES1_2
DAISY13$AC_3
DAISY13$OPC_4
DAISY13$NPC1_5
DAISY13$NPC2_6

#Escape遅いけど、ちゃんと計算しているらしい
feature_list_mouse <- list(MES2=mMES2, MES1=mMES1, AC=mAC, OPC=mOPC, NPC1=mNPC1, NPC2=mNPC2)
GS.hallmark <- getGeneSets(library = "H")

#escape score計算,group = 細胞数, cores = 並列処理, ssGSEA.norm = 正規化で0-1、デフォルトだとssGSEAでenrichment analysis。
escape.DAISY13 <- enrichIt(DAISY13@assays$RNA$counts, #v5から、指定先変更
                           gene.sets = feature_list_mouse, cores = 8, min.size = 3, ssGSEA.norm = TRUE)
#escape.DAISY13re <- escape.matrix(DAISY13@assays$RNA$counts, #v5から、指定先変更
#                                  gene.sets = feature_list_mouse, method = "ssGSEA", cores = 8, min.size = 3, ssGSEA.norm = TRUE)

# 各行で最大値を持つカラム名を取得する
max_col <- apply(escape.DAISY13, 1, function(row) {
  colnames(escape.DAISY13)[which.max(row)]
})
# 新しいカラムを追加する
escape.DAISY13$module1 <- max_col
# 新しいカラムを追加し、条件に基づいてカラム名を変更する
escape.DAISY13$module2 <- apply(escape.DAISY13, 1, function(row) {
  max_col <- colnames(escape.DAISY13)[which.max(row)]
  if (grepl("MES", max_col)) {
    return("MES")
  } else if (grepl("NPC", max_col)) {
    return("NPC")
  } else {
    return(max_col)
  }
})

DAISY13 <- AddMetaData(DAISY13, metadata = escape.DAISY13)
DAISY13@meta.data

saveRDS(DAISY13, file = paste0("DAISY13_BC",".rds"))

module1_order = c("AC","MES1","MES2","OPC","NPC1","NPC2")
DAISY13$module1 <- factor(x = DAISY13$module1, levels = module1_order)#factor指定すると並び替えられる。

module2_order = c("AC","MES","OPC","NPC")
DAISY13$module2 <- factor(x = DAISY13$module2, levels = module2_order)#factor指定すると並び替えられる。


DF <- data.frame(module1 = DAISY13@meta.data$module1, BCnumber = DAISY13@meta.data$BCnumber, seurat_clusters = DAISY13@meta.data$seurat_clusters)

# BCnumberごとにmoduleの割合を計算
module1_percent <- DF %>%
  group_by(BCnumber, module1) %>%
  tally() %>%  # tally()を使用して行数を数える
  mutate(percent = n / sum(n) * 100)
# BCnumberごとに全部の和が100%になるか確認
module1_percent <- module1_percent %>%
  group_by(BCnumber) %>%
  mutate(total_percent = sum(percent))
# グラフ作成
BCarbitrary <- paste0("BC", 1:50)
module1_percent <- module1_percent %>% filter(BCnumber %in% BCarbitrary) %>% mutate(BCnumber = factor(BCnumber, levels = BCarbitrary))
ggplot(module1_percent, aes(x = BCnumber, y = percent, fill = module1)) +
  geom_bar(stat = "identity") +
  #facet_wrap(~BCnumber, scales = "free") +
  labs(title = "Module-like per Barcode Clone", x = "BCnumber", y = "Percentage") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # x軸のラベルを回転させて見やすくする







#### Monocle3 ####
set.seed(605)
DAISY13_BC <- SetIdent(DAISY13_BC, value = "BCnumber")
assign("data", DAISY13_BC, envir = .GlobalEnv) 

## Monocle3での解析
#Converting seuratobject to celldataset object for Monocle3
cds <- as.cell_data_set(data)
#Get feature/gene metadata
fData(cds)$gene_short_name <- rownames(fData(cds))

cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)#でPCAの確認

cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds)#UMAP

plot_cells(cds)
plot_cells(cds, color_cells_by="Phase") #"module2"
plot_cells(cds, color_cells_by="BCnumber")
plot_cells(cds, genes=c("Cd274", "Clpp"))

cds <- cluster_cells(cds)
#cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") #"cluster"

cds <- learn_graph(cds)
plot_cells(cds,color_cells_by = "cluster",
           label_groups_by_cluster=TRUE, group_label_size=3,
           label_leaves=FALSE, label_branch_points=FALSE)

cds <- order_cells(cds)#manualで選ぶ

# # a helper function to identify the root principal points:
# get_earliest_principal_node <- function(cds, time_bin="130-170"){
#   cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
#   
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }
# cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=TRUE,
           label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)


## Seuratから情報横領する場合
#Converting seuratobject to celldataset object for Monocle3
cds <- as.cell_data_set(data)
#Get feature/gene metadata
fData(cds)$gene_short_name <- rownames(fData(cds))
# #Get cell metadata
# head(colData(cds))
# #Get counts
# head(counts(cds))

#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
#Assign cluster information
list.cluster <- data@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
#Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE, group_label_size = 5)# + theme(legend.position = "right")

#cds <- learn_graph(cds, use_partition = TRUE) #use_partion=FALSEで全細胞を解析、defaultはTRUEつまり何らかのcluster
cds <- learn_graph(cds, use_partition = FALSE) #use_partion=FALSEで全細胞を解析、defaultはTRUEつまり何らかのcluster
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
           label_branch_points = TRUE, label_roots = TRUE, label_leaves = FALSE, group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 1]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()


BCarbitrary <- paste0("BC", 1:50)
data.pseudo <- data.pseudo %>% filter(BCnumber %in% BCarbitrary) %>% mutate(BCnumber = factor(BCnumber, levels = BCarbitrary))
ggplot(data.pseudo, aes(monocle3_pseudotime, BCnumber, fill = BCnumber)) + geom_boxplot() +NoLegend()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(BCnumber, monocle3_pseudotime), fill = BCnumber)) + geom_boxplot()

#Find genes that change as a function of pseudotime
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()


#### CytoTRACE2 ####
set.seed(727)
DAISY13_BC <- SetIdent(DAISY13_BC, value = "BCnumber")
assign("data", DAISY13_BC, envir = .GlobalEnv) 

#expression_data <- data@assays$RNA$counts #Raw count数でも大丈夫
expression_data <- data[["RNA"]]$counts #Raw count数でも大丈夫
#expression_data <- data[["RNA"]]$data #Normalized count数でも大丈夫

cytotrace2_result <- cytotrace2(expression_data) #デフォルトはマウス、species = "human"

DAISY13_BC <- AddMetaData(DAISY13_BC, metadata = cytotrace2_result)
DAISY13_BC@meta.data

saveRDS(DAISY13_BC, file = paste0("DAISY13_BC",".rds"))

DimPlot_scCustom(DAISY13_BC, reduction = "umap", group.by = c("CytoTRACE2_Potency"))#+NoLegend()
#FeaturePlot(DAISY13_BC, feature=c("CytoTRACE2_Score","CytoTRACE2_Relative"), ncol=2)
pal <- viridis(n = 256, option = "B") #"A": magma,"B": inferno,"C": plasma,"D": viridis,"E": cividis
FeaturePlot_scCustom(DAISY13_BC, feature=c("CytoTRACE2_Score","CytoTRACE2_Relative"), colors_use = pal)


DAISY13_BC@meta.data$BCnumber[is.na(DAISY13_BC@meta.data$BCnumber)] <- "nonBC"
table(DAISY13_BC$BCnumber)
DAISY13_BC_33 <- subset(DAISY13_BC, subset = BCnumber %in% c("BC33")) #%in%左側のオブジェクトが右側のオブジェクトの中に存在するかどうかを確認する

plot <- FeaturePlot_scCustom(DAISY13_BC_33, feature=c("CytoTRACE2_Score","CytoTRACE2_Relative"), colors_use = pal, pt.size=3, label = TRUE)
plot <- FeaturePlot_scCustom(DAISY13_BC_33, feature=c("CytoTRACE2_Score"), colors_use = pal, pt.size=3, label = FALSE)
LabelPoints(plot = plot, points = TopCells(object = DAISY13_BC_33[["umap"]], ncells=100), repel = TRUE, size=2) #10xbarcodeを書き込む

#こんな10xbarcodeのラベルの仕方もある
# Extract coordinates for barcodes
umap_coords <- Embeddings(DAISY13_BC_33, "umap")
barcode_labels <- rownames(umap_coords)
plot + geom_text(data = umap_coords, aes(x = umap_1, y = umap_2, label = barcode_labels), size = 2, vjust = 1, hjust = 1)








#### slingshot ####



