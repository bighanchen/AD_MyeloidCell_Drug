# author: Xiaowei Chen (chenxiaowei@ibp.ac.cn)
# import packages
library(Seurat)
library(VISION)
library(data.table)
library(DoubletFinder)
library(tidyverse)
library(sva)

#GSE143758
mtx <- fread("Public_scRNAseq/GSE143758/GSE143758_Admouse_Hippocampus_7m_AllNuclei_UMIcounts.txt", data.table = FALSE)
rownames(mtx) <- mtx[,1]
mtx <- mtx[, -1]

#DoubletFinder
hippo_pro <- CreateSeuratObject(mtx, project = "SeuratProject", assay = "RNA",
  min.cells = 10, min.features = 400, names.field = 1,
  names.delim = "+", meta.data = NULL)
 
hippo_pro <- SCTransform(hippo_pro)
hippo_pro <- RunPCA(hippo_pro, features = VariableFeatures(object = hippo_pro))
pc.num=1:15
hippo_pro <- RunUMAP(hippo_pro, dims=pc.num)
hippo_pro <- FindNeighbors(hippo_pro, dims = pc.num) %>% FindClusters(resolution = 1)

sweep.res.list <- paramSweep_v3(hippo_pro, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate = 0.012
homotypic.prop <- modelHomotypic(hippo_pro$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(hippo_pro))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

hippo_pro_dF <- doubletFinder_v3(hippo_pro, PCs = pc.num, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
mtx_dF <- mtx[,rownames(hippo_pro_dF@meta.data[hippo_pro_dF@meta.data$DF.classifications_0.25_0.21_619 =='Singlet',])]

#CCA
hippo_dF <- CreateSeuratObject(mtx_dF, project = "SeuratProject", assay = "RNA",
  min.cells = 10, min.features = 400, names.field = 1,
  names.delim = "+", meta.data = NULL)
 
batch <- c()
status <- c()
i=1
for (item in colnames(mtx_dF)){    
    sp1 <- strsplit(item,'_')
    bat <- strsplit(sp1[[1]][1],'\\.')[[1]][2]
    if (bat == 'Batch1'){
        batch[i]=1
    } else if(bat=='Batch3'){
        batch[i]=3
    }
    
    sp2 <- strsplit(sp1[[1]][2],'-')
    if (sp2[[1]][2] == 'WT' | sp2[[1]][1] == 'Wt'){
        status[i] <- 'WT'
    }else if(sp2[[1]][2] == 'AD'| sp2[[1]][1] == 'Untreated'){
        status[i] <- 'AD'
    }
    i=i+1    
}

hippo_dF@meta.data[,'batch'] <- batch
hippo_dF@meta.data[,'status'] <- status
hippo_dF.list <- SplitObject(hippo_dF, split.by = "batch")

# normalize and identify variable features for each dataset independently
hippo_dF.list <- lapply(X = hippo_dF.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features_integ <- SelectIntegrationFeatures(object.list = hippo_dF.list)

anchors <- FindIntegrationAnchors(object.list = hippo_dF.list, anchor.features = features_integ)

# this command creates an 'integrated' data assay
hippo_dF_integ <- IntegrateData(anchorset = anchors)

DefaultAssay(hippo_dF_integ) <- "integrated"

hippo_dF_integ <- ScaleData(hippo_dF_integ, verbose = FALSE)
hippo_dF_integ <- RunPCA(hippo_dF_integ, npcs = 30, verbose = FALSE)
hippo_dF_integ <- RunUMAP(hippo_dF_integ, reduction = "pca", dims = 1:30)
hippo_dF_integ <- FindNeighbors(hippo_dF_integ, reduction = "pca", dims = 1:30)
hippo_dF_integ <- FindClusters(hippo_dF_integ, resolution = 0.5)


new.cluster.ids <- c(
                          'ExN.DG',#0 Stxbp6
                          'Astrocytes',#1 Slc1a3
                          'ExN.CA1.1',#2
                          'ExN.CA3.1',#3
                          'ExN.Sub',#4 
                          'ExN.CA1.2',#5
                          'GABAergic.1',#6
                          'ExN.CA1.3',#7
                          'Ex.Neuron',#8 
                          'GABAergic.2',#9
                          'Oligodendrocyte precursor cells',#10 Pdgfra
                          'Endothelial',#11 Cldn5
                          'Ex.Neuron',#12 
                          'Oligodendrocytes',#13 Plp1
                          'GABAergic.3',#14
                          'Ependymal',#15 Ccdc153
                          'ExN.CA3.2',#16 
                          'Microglia',#17 Hmha1
                          'Microglia',#18 Hmha1
                          'ExN.CA3.3',#19
                          'Pericytes',#20 Pdgfrb, Vtn
                          'ExN.IEGs'  #21
)
names(new.cluster.ids) <- levels(hippo_dF_integ)
hippo_dF_integ <- RenameIdents(hippo_dF_integ, new.cluster.ids)
DimPlot(hippo_dF_integ, reduction = "umap", label = T)

hippo_dF_integ@active.ident <- factor(hippo_dF_integ@active.ident, 
                            levels=c('Astrocytes', 
                                     'Ependymal',
                                     'Endothelial', 
                                     'Oligodendrocytes', 
                                     'GABAergic.1', 
                                     'GABAergic.2', 
                                     'GABAergic.3', 
                                     'Oligodendrocyte precursor cells', 
                                     'Microglia', 
                                     'ExN.DG',
                                     'Pericytes',
                                     'Ex.Neuron',
                                     'ExN.CA3.1',
                                     'ExN.CA3.2',
                                     'ExN.CA3.3',
                                     'ExN.Sub',
                                     'ExN.CA1.1',
                                     'ExN.CA1.2',
                                     'ExN.CA1.3',
                                     "ExN.IEGs"))
features <- c('Homer1','Ndst4','Pex5l','Col6a1','Tshz2','Mgat4c', 'Slc6a13','Stxbp6','Hmha1','Pdgfra','Gad2','Plp1','Cldn5','Rarres2','Gfap','Slc1a3')
DefaultAssay(hippo_dF_integ) <- "RNA"
DotPlot(hippo_dF_integ, features = features) + RotatedAxis()

#Signature Scores
load('Signature_repository/all_deg_sig_list.bin')
DefaultAssay(hippo_dF_integ) <- 'RNA'
hippo_dF_integ_SigScore_all_deg <- Vision(hippo_dF_integ, signatures = all_deg_sig_list, assay='RNA')
hippo_dF_integ_SigScore_all_deg <- analyze(hippo_dF_integ_SigScore_all_deg)

for (signature_set in colnames(hippo_dF_integ_SigScore_all_deg@SigScores)) {
  
  sigScores <- hippo_dF_integ_SigScore_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  hippo_dF_integ@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(hippo_dF_integ@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

batch <- c()
status <- c()
i=1
for (item in rownames(hippo_dF_integ@meta.data)){    
    sp1 <- strsplit(item,'_')
    bat <- strsplit(sp1[[1]][1],'\\.')[[1]][2]
    batch[i] <- bat
    sp2 <- strsplit(sp1[[1]][2],'-')
    if (sp2[[1]][2] == 'WT' | sp2[[1]][1] == 'Wt'){
        status[i] <- 'WT'
    }else if(sp2[[1]][2] == 'AD'| sp2[[1]][1] == 'Untreated'){
        status[i] <- 'AD'
    }
    i=i+1    
}

hippo_dF_integ@meta.data[,'status'] <-status
hippo_dF_integ$CellType <- Idents(hippo_dF_integ)
all_deg_celltype_vln <- hippo_dF_integ@meta.data

all_deg_celltype_vln_df <- data.frame(score=hippo_dF_integ@meta.data$all_deg_sig,
                                      celltype=factor(hippo_dF_integ@meta.data$CellType,levels = c('Microglia','Endothelial','Astrocytes','Pericytes','Ependymal','Oligodendrocytes','Oligodendrocyte precursor cells','ExN.CA3.2','ExN.DG','GABAergic.2','GABAergic.3','Ex.Neuron','ExN.Sub','ExN.CA1.1','ExN.CA3.3','GABAergic.1','ExN.CA1.2','ExN.CA3.1','ExN.CA1.3','ExN.IEGs')),
                                      status=factor(hippo_dF_integ@meta.data$status,levels = c("WT","AD")))

options(repr.plot.width=28, repr.plot.height=7)
ggplot(na.omit(all_deg_celltype_vln_df), aes(x=celltype, y=score , fill = status)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1))

load('Signature_repository/region_common_all_deg_sig_list.bin')
DefaultAssay(hippo_dF_integ) <- 'RNA'
hippo_dF_integ_SigScore_region_all_deg <- Vision(hippo_dF_integ, signatures = region_common_all_deg_sig_list, assay='RNA')
hippo_dF_integ_SigScore_region_all_deg <- analyze(hippo_dF_integ_SigScore_region_all_deg)
for (signature_set in colnames(hippo_dF_integ_SigScore_region_all_deg@SigScores)) {
  
  sigScores <- hippo_dF_integ_SigScore_region_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  hippo_dF_integ@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(hippo_dF_integ@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

region_all_deg_vln_df <- data.frame(score=hippo_dF_integ@meta.data$region_common_all_deg_sig,
                                      celltype=factor(hippo_dF_integ@meta.data$CellType,levels = c('Microglia','Endothelial','Astrocytes','Oligodendrocytes','Ependymal','Pericytes','Oligodendrocyte precursor cells','ExN.DG','ExN.CA1.1','GABAergic.3','Ex.Neuron','GABAergic.2','ExN.CA1.2','ExN.Sub','ExN.CA3.1','GABAergic.1','ExN.CA1.3','ExN.CA3.2','ExN.CA3.3','ExN.IEGs')),
                                      status=factor(hippo_dF_integ@meta.data$status,levels = c("WT","AD")))
options(repr.plot.width=28, repr.plot.height=7)
ggplot(na.omit(region_all_deg_vln_df), aes(x=celltype, y=score , fill = status)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1))






