# author: Xiaowei Chen (chenxiaowei@ibp.ac.cn)
# import packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(VISION)
library(plyr)
library(patchwork)
library(topGO)
library(clusterProfiler)
library(qs)
library(DropletUtils)
library(SeuratDisk)

#Prepare the data
region_anno <- read.csv('STAligner/region.csv',header = FALSE)
colnames(region_anno) <- c('subregion','region')

st_ad_list <- c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')

spot_obj_list <- list()

for (samp in st_ad_list) {
    
    #load data
    st <- Load10X_Spatial(data.dir = samp,filename = 'filtered_feature_bc_matrix.h5')
    # normalization
    DefaultAssay(st) <- "Spatial"
    st <- SCTransform(st, assay = "Spatial", verbose = TRUE)
    DefaultAssay(st) <- 'SCT'
    st <- RunPCA(object = st, verbose = T)
    st <- RunUMAP(st, dims = 1:30, verbose = T)
        
    #read spot coordinates from space ranger output
    tissue_positions <- read.csv(paste0(samp,'/spatial/tissue_positions_list.csv'),header=FALSE)
    colnames(tissue_positions)<-c('barcode','covered','row','column','x','y')
    tissue_positions <- unite(tissue_positions, 'row_col', row, column, sep = '_', remove = FALSE)
    
    #read Abeta intensity measurement
    lower_abeta_intensity <- read.csv(paste0(samp,'/',samp,'_lower_abeta_intensity.csv'),header =TRUE)
    #read spot coordinates from Imaris 
    lower_position <- read.csv(paste0(samp,'/',samp,'_lower_position.csv'),header =TRUE)
    #get Abeta intensity measurement for each spot
    lower_abeta_intensity_position <- full_join(lower_position,lower_abeta_intensity,by = 'ID')
    lower_abeta_intensity_position <- unite(lower_abeta_intensity_position, 'row_col', row, col, sep = '_', remove = FALSE)
    lower_abeta_intensity_position <- left_join(lower_abeta_intensity_position,tissue_positions,by = 'row_col')
    lower_abeta <- lower_abeta_intensity_position[lower_abeta_intensity_position$covered == 1,c('barcode','Mean')]
    colnames(lower_abeta) <- c('barcode','mean_lower')
    
    #read Abeta intensity measurement
    upper_abeta_intensity <- read.csv(paste0(samp,'/',samp,'_upper_abeta_intensity.csv'),header =TRUE)
    #read spot coordinates from Imaris
    upper_position <- read.csv(paste0(samp,'/',samp,'_upper_position.csv'),header =TRUE)
    #get Abeta intensity measurement for each spot
    upper_abeta_intensity_position <- full_join(upper_position,upper_abeta_intensity,by = 'ID')
    upper_abeta_intensity_position <- unite(upper_abeta_intensity_position, 'row_col', row, col, sep = '_', remove = FALSE)
    upper_abeta_intensity_position <- left_join(upper_abeta_intensity_position,tissue_positions,by = 'row_col')
    upper_abeta <- upper_abeta_intensity_position[upper_abeta_intensity_position$covered == 1,c('barcode','Mean')]
    colnames(upper_abeta) <- c('barcode','mean_upper')
    
    mean_abeta <- inner_join(lower_abeta,upper_abeta,by = 'barcode')
    mean_abeta <- mutate(mean_abeta, average = (mean_lower + mean_upper) / 2)
    
    #get intensity mean
    mean_meta <- data.frame(barcode = mean_abeta[,c('barcode')],mean = mean_abeta[,c('average')])
    # remove rows containing na, nan, inf, -inf
    mean_meta <- mean_meta[!is.na(mean_meta$mean)&!is.nan(mean_meta$mean)&!is.infinite(mean_meta$mean),]
    
    region_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(region_meta) <- c("barcode","subregion")
    region_meta$region <- mapvalues(x = region_meta$subregion, from = region_anno$subregion, to = region_anno$region)
    region_meta <- region_meta[region_meta$subregion != 'VL',]
    
    spot_meta <- inner_join(mean_meta,region_meta,by = 'barcode')
    
    rownames(spot_meta)<-spot_meta$barcode
    spot_meta <- spot_meta[,c('subregion','region','mean'),drop = FALSE]
    spot_meta$subregion <- factor(spot_meta$subregion, levels = c("L1","L2/3","L4","L5","L6","PIR","Fxs+Int","opt","CC","CA-1","CA-2","DG","VL","CTXsp-1","CTXsp-2","STRD-1","STRD-2","SI","sAMY","RT","TH-1","TH-2","TH-3","DORsm","HY"))
    spot_meta$region <- factor(spot_meta$region, levels = c('cerebral_cortex','fiber_tracts','hippocampal_region','ventricular_systems','cerebral_nuclei','thalamus','hypothalamus'))
   
    # add annotation to metadata
    st <- AddMetaData(
        object = st,
        metadata = spot_meta
    )
    
    spot_obj_list[[samp]] <- st
} 

#Signature scores of genes affected by both abeta and aging: Global

ST_3M_AD_1_2 <- spot_obj_list[['ST_3M_AD_1_2']]
ST_6M_AD_2_1 <- spot_obj_list[['ST_6M_AD_2_1']]
ST_15M_AD_2_1 <- spot_obj_list[['ST_15M_AD_2_1']]

load('Signature_repository/all_deg_sig_list.bin')
DefaultAssay(ST_3M_AD_1_2) <- 'SCT'
ST_3M_AD_1_2_SigScore_all_deg <- Vision(ST_3M_AD_1_2, signatures = all_deg_sig_list, assay='SCT')
ST_3M_AD_1_2_SigScore_all_deg <- analyze(ST_3M_AD_1_2_SigScore_all_deg)
for (signature_set in colnames(ST_3M_AD_1_2_SigScore_all_deg@SigScores)) {
  
  sigScores <- ST_3M_AD_1_2_SigScore_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_3M_AD_1_2@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_3M_AD_1_2@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_6M_AD_2_1) <- 'SCT'
ST_6M_AD_2_1_SigScore_all_deg <- Vision(ST_6M_AD_2_1, signatures = all_deg_sig_list, assay='SCT')
ST_6M_AD_2_1_SigScore_all_deg <- analyze(ST_6M_AD_2_1_SigScore_all_deg)
for (signature_set in colnames(ST_6M_AD_2_1_SigScore_all_deg@SigScores)) {
  
  sigScores <- ST_6M_AD_2_1_SigScore_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_6M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_6M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_15M_AD_2_1) <- 'SCT'
ST_15M_AD_2_1_SigScore_all_deg <- Vision(ST_15M_AD_2_1, signatures = all_deg_sig_list, assay='SCT')
ST_15M_AD_2_1_SigScore_all_deg <- analyze(ST_15M_AD_2_1_SigScore_all_deg)
for (signature_set in colnames(ST_15M_AD_2_1_SigScore_all_deg@SigScores)) {
  
  sigScores <- ST_15M_AD_2_1_SigScore_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_15M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_15M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=21, repr.plot.height=7)

min_aggr_all_deg <- min(min(ST_3M_AD_1_2$all_deg_sig,na.rm = TRUE),min(ST_6M_AD_2_1$all_deg_sig,na.rm = TRUE),min(ST_15M_AD_2_1$all_deg_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(ST_3M_AD_1_2$all_deg_sig,na.rm = TRUE),max(ST_6M_AD_2_1$all_deg_sig,na.rm = TRUE),max(ST_15M_AD_2_1$all_deg_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(ST_3M_AD_1_2, features = 'all_deg_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_6M_AD_2_1, features = 'all_deg_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_15M_AD_2_1, features = 'all_deg_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )

all_deg_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$all_deg_sig,ST_6M_AD_2_1@meta.data$all_deg_sig,ST_15M_AD_2_1@meta.data$all_deg_sig),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))
							 
							 options(repr.plot.width=4, repr.plot.height=7)
ggplot(all_deg_vln_df, aes(x=section, y=score)) + 
geom_violin(trim=FALSE)+
stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.05 )

all_deg_vln_df_3Mvs6M <- t.test(ST_3M_AD_1_2@meta.data$all_deg_sig,ST_6M_AD_2_1@meta.data$all_deg_sig)
all_deg_vln_df_6Mvs15M <- t.test(ST_6M_AD_2_1@meta.data$all_deg_sig,ST_15M_AD_2_1@meta.data$all_deg_sig)

all_deg_region_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$all_deg_sig,ST_6M_AD_2_1@meta.data$all_deg_sig,ST_15M_AD_2_1@meta.data$all_deg_sig),
                                    region=factor(c(ST_3M_AD_1_2@meta.data$region,ST_6M_AD_2_1@meta.data$region,ST_15M_AD_2_1@meta.data$region),levels = c('fiber_tracts','hippocampal_region','thalamus','cerebral_nuclei','cerebral_cortex','ventricular_systems','hypothalamus')),
                                    subregion=factor(c(ST_3M_AD_1_2@meta.data$subregion,ST_6M_AD_2_1@meta.data$subregion,ST_15M_AD_2_1@meta.data$subregion),levels = c("opt","CC","Fxs+Int","CA-2","SI","DG","TH-2","CTXsp-1","DORsm","L1","STRD-1","L6","L5","RT","CTXsp-2","TH-1","PIR","TH-3","sAMY","HY","CA-1","STRD-2","L2/3","L4","VL")),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))
							 

options(repr.plot.width=28, repr.plot.height=14)
ggplot(na.omit(all_deg_region_vln_df), aes(x=region, y=score , fill = section)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1)),aes(group=1))

options(repr.plot.width=28, repr.plot.height=5)
ggplot(na.omit(all_deg_region_vln_df), aes(x=subregion, y=score, fill = section)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) 



#GO enrichment analysis
all_deg_global <- read.csv('Signature_repository/all_deg_global.csv',row.names = 1, header = TRUE)
all_deg_global_up <- all_deg_global[all_deg_global$logFC_3m>0,]
all_deg_global_down <- all_deg_global[all_deg_global$logFC_3m<0,]
all_deg_global_up_gene_id=bitr(all_deg_global_up$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 
all_deg_global_down_gene_id=bitr(all_deg_global_down$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

all_deg_global_up_GO_BP <- enrichGO(gene          = all_deg_global_up_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)

all_deg_global_up_GO_MF <- enrichGO(gene          = all_deg_global_up_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)

all_deg_global_up_GO_CC <- enrichGO(gene          = all_deg_global_up_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)

all_deg_global_up_GO_BP_simplify <- simplify(all_deg_global_up_GO_BP,cutoff = 0.5)
write.csv(all_deg_global_up_GO_BP_simplify@result,'all_deg_global_up_GO_BP_simplify.csv',quote=FALSE)

all_deg_global_up_GO_BP_simplify_plot <- all_deg_global_up_GO_BP_simplify
all_deg_global_up_GO_BP_simplify_plot@result <- all_deg_global_up_GO_BP_simplify_plot@result[c(1,2,3,4,5,13,23,34),]
options(repr.plot.width=7, repr.plot.height=7)

barplot(all_deg_global_up_GO_BP_simplify_plot, title="all_deg_global_up_GO_BP")+
theme(panel.background = element_blank(),panel.grid = element_blank() )

all_deg_global_up_GO_MF_simplify <- simplify(all_deg_global_up_GO_MF,cutoff = 0.5)
write.csv(all_deg_global_up_GO_MF_simplify@result,'all_deg_global_up_GO_MF_simplify.csv',quote=FALSE)

all_deg_global_up_GO_CC_simplify <- simplify(all_deg_global_up_GO_CC,cutoff = 0.5)
write.csv(all_deg_global_up_GO_CC_simplify@result,'all_deg_global_up_GO_CC_simplify.csv',quote=FALSE)

all_deg_global_down_GO_BP <- enrichGO(gene          = all_deg_global_down_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)

all_deg_global_down_GO_MF <- enrichGO(gene          = all_deg_global_down_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)

all_deg_global_down_GO_CC <- enrichGO(gene          = all_deg_global_down_gene_id$ENTREZID,
                #universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.01,
                #qvalueCutoff  = 0.05,
        readable      = TRUE)
all_deg_global_down_GO_BP_simplify <- simplify(all_deg_global_down_GO_BP,cutoff = 0.5)
write.csv(all_deg_global_down_GO_BP_simplify@result,'all_deg_global_down_GO_BP_simplify.csv',quote=FALSE)


all_deg_global_down_GO_BP_simplify_plot <- all_deg_global_down_GO_BP_simplify
all_deg_global_down_GO_BP_simplify_plot@result <- all_deg_global_down_GO_BP_simplify_plot@result[c(1,2,3,4,5,6,21,22),]
options(repr.plot.width=7, repr.plot.height=7)

barplot(all_deg_global_down_GO_BP_simplify_plot, title="all_deg_global_down_GO_BP")+
theme(panel.background = element_blank(),panel.grid = element_blank() )

#Signature scores of genes affected by both abeta and aging: Regional
load('Signature_repository/region_common_all_deg_sig_list.bin')
DefaultAssay(ST_3M_AD_1_2) <- 'SCT'
ST_3M_AD_1_2_SigScore_region_all_deg <- Vision(ST_3M_AD_1_2, signatures = region_common_all_deg_sig_list, assay='SCT')
ST_3M_AD_1_2_SigScore_region_all_deg <- analyze(ST_3M_AD_1_2_SigScore_region_all_deg)
for (signature_set in colnames(ST_3M_AD_1_2_SigScore_region_all_deg@SigScores)) {
  
  sigScores <- ST_3M_AD_1_2_SigScore_region_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_3M_AD_1_2@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_3M_AD_1_2@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_6M_AD_2_1) <- 'SCT'
ST_6M_AD_2_1_SigScore_region_all_deg <- Vision(ST_6M_AD_2_1, signatures = region_common_all_deg_sig_list, assay='SCT')
ST_6M_AD_2_1_SigScore_region_all_deg <- analyze(ST_6M_AD_2_1_SigScore_region_all_deg)
for (signature_set in colnames(ST_6M_AD_2_1_SigScore_region_all_deg@SigScores)) {
  
  sigScores <- ST_6M_AD_2_1_SigScore_region_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_6M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_6M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_15M_AD_2_1) <- 'SCT'
ST_15M_AD_2_1_SigScore_region_all_deg <- Vision(ST_15M_AD_2_1, signatures = region_common_all_deg_sig_list, assay='SCT')
ST_15M_AD_2_1_SigScore_region_all_deg <- analyze(ST_15M_AD_2_1_SigScore_region_all_deg)
for (signature_set in colnames(ST_15M_AD_2_1_SigScore_region_all_deg@SigScores)) {
  
  sigScores <- ST_15M_AD_2_1_SigScore_region_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_15M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_15M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=21, repr.plot.height=7)

min_aggr_all_deg <- min(min(ST_3M_AD_1_2$region_common_all_deg_sig,na.rm = TRUE),min(ST_6M_AD_2_1$region_common_all_deg_sig,na.rm = TRUE),min(ST_15M_AD_2_1$region_common_all_deg_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(ST_3M_AD_1_2$region_common_all_deg_sig,na.rm = TRUE),max(ST_6M_AD_2_1$region_common_all_deg_sig,na.rm = TRUE),max(ST_15M_AD_2_1$region_common_all_deg_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(ST_3M_AD_1_2, features = 'region_common_all_deg_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_6M_AD_2_1, features = 'region_common_all_deg_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_15M_AD_2_1, features = 'region_common_all_deg_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )
region_all_deg_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$region_common_all_deg_sig,ST_6M_AD_2_1@meta.data$region_common_all_deg_sig,ST_15M_AD_2_1@meta.data$region_common_all_deg_sig),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))
							 
options(repr.plot.width=4, repr.plot.height=7)
ggplot(region_all_deg_vln_df, aes(x=section, y=score)) + 
geom_violin(trim=FALSE)+
stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.05 )

region_all_deg_region_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$region_common_all_deg_sig,ST_6M_AD_2_1@meta.data$region_common_all_deg_sig,ST_15M_AD_2_1@meta.data$region_common_all_deg_sig),
                                    region=factor(c(ST_3M_AD_1_2@meta.data$region,ST_6M_AD_2_1@meta.data$region,ST_15M_AD_2_1@meta.data$region),levels = c('fiber_tracts','thalamus','hippocampal_region','cerebral_nuclei','cerebral_cortex','hypothalamus')),
                                    subregion=factor(c(ST_3M_AD_1_2@meta.data$subregion,ST_6M_AD_2_1@meta.data$subregion,ST_15M_AD_2_1@meta.data$subregion),levels = c("opt","CC","Fxs+Int","DORsm","CA-2","TH-2","SI","DG","CTXsp-1","RT","L6","L5","CTXsp-2","STRD-1","L1","TH-3","TH-1","PIR","CA-1","STRD-2","HY","sAMY","L2/3","L4")),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))


options(repr.plot.width=28, repr.plot.height=14)
ggplot(na.omit(region_all_deg_region_vln_df), aes(x=region, y=score , fill = section)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1))

options(repr.plot.width=28, repr.plot.height=5)
ggplot(na.omit(region_all_deg_region_vln_df), aes(x=subregion, y=score, fill = section)) + 
geom_boxplot(position=position_dodge(1),outlier.shape = NA) 




#Signature scores of CAS

CAS <- read.csv("Signature_repository/CAS.list",header=FALSE)
colnames(CAS) <- c("gene","sign")

CAS_sig_list <- c()
signature_vector <- c(CAS$sign)
names(signature_vector) <- c(CAS$gene)
signature_vector <- signature_vector[!(duplicated(names(signature_vector)) | duplicated(names(signature_vector),fromLast=TRUE))]
CAS_sig_list <- c(CAS_sig_list,createGeneSignature(name = "CAS_sig", sigData = signature_vector))
save(CAS_sig_list, file='Signature_repository/CAS_sig_list.bin')							 

DefaultAssay(ST_3M_AD_1_2) <- 'SCT'
ST_3M_AD_1_2_SigScore_cas <- Vision(ST_3M_AD_1_2, signatures = CAS_sig_list, assay='SCT')
ST_3M_AD_1_2_SigScore_cas <- analyze(ST_3M_AD_1_2_SigScore_cas)
for (signature_set in colnames(ST_3M_AD_1_2_SigScore_cas@SigScores)) {
  
  sigScores <- ST_3M_AD_1_2_SigScore_cas@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_3M_AD_1_2@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_3M_AD_1_2@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_6M_AD_2_1) <- 'SCT'
ST_6M_AD_2_1_SigScore_cas <- Vision(ST_6M_AD_2_1, signatures = CAS_sig_list, assay='SCT')
ST_6M_AD_2_1_SigScore_cas <- analyze(ST_6M_AD_2_1_SigScore_cas)
for (signature_set in colnames(ST_6M_AD_2_1_SigScore_cas@SigScores)) {
  
  sigScores <- ST_6M_AD_2_1_SigScore_cas@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_6M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_6M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

DefaultAssay(ST_15M_AD_2_1) <- 'SCT'
ST_15M_AD_2_1_SigScore_cas <- Vision(ST_15M_AD_2_1, signatures = CAS_sig_list, assay='SCT')
ST_15M_AD_2_1_SigScore_cas <- analyze(ST_15M_AD_2_1_SigScore_cas)
for (signature_set in colnames(ST_15M_AD_2_1_SigScore_cas@SigScores)) {
  
  sigScores <- ST_15M_AD_2_1_SigScore_cas@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  ST_15M_AD_2_1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(ST_15M_AD_2_1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=21, repr.plot.height=7)

min_aggr_all_deg <- min(min(ST_3M_AD_1_2$CAS_sig,na.rm = TRUE),min(ST_6M_AD_2_1$CAS_sig,na.rm = TRUE),min(ST_15M_AD_2_1$CAS_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(ST_3M_AD_1_2$CAS_sig,na.rm = TRUE),max(ST_6M_AD_2_1$CAS_sig,na.rm = TRUE),max(ST_15M_AD_2_1$CAS_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(ST_3M_AD_1_2, features = 'CAS_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_6M_AD_2_1, features = 'CAS_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg)),
     SpatialFeaturePlot(ST_15M_AD_2_1, features = 'CAS_sig')+
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )
cas_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$CAS_sig,ST_6M_AD_2_1@meta.data$CAS_sig,ST_15M_AD_2_1@meta.data$CAS_sig),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))
							 
options(repr.plot.width=4, repr.plot.height=7)
ggplot(cas_vln_df, aes(x=section, y=score)) + 
geom_violin(trim=FALSE)+
stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.05 )

cas_region_vln_df <- data.frame(score=c(ST_3M_AD_1_2@meta.data$CAS_sig,ST_6M_AD_2_1@meta.data$CAS_sig,ST_15M_AD_2_1@meta.data$CAS_sig),
                                    region=factor(c(ST_3M_AD_1_2@meta.data$region,ST_6M_AD_2_1@meta.data$region,ST_15M_AD_2_1@meta.data$region),levels = c('fiber_tracts','thalamus','hippocampal_region','cerebral_nuclei','cerebral_cortex','hypothalamus')),
                                    subregion=factor(c(ST_3M_AD_1_2@meta.data$subregion,ST_6M_AD_2_1@meta.data$subregion,ST_15M_AD_2_1@meta.data$subregion),levels = c("opt","CC","Fxs+Int","DORsm","CA-2","TH-2","SI","DG","CTXsp-1","RT","L6","L5","STRD-1","CTXsp-2","TH-1","L1","TH-3","PIR","CA-1","STRD-2","HY","sAMY","L2/3","L4")),
                             section=factor(c(rep('ST_3M_AD_1_2',dim(ST_3M_AD_1_2@meta.data)[1]),rep('ST_6M_AD_2_1',dim(ST_6M_AD_2_1@meta.data)[1]),rep('ST_15M_AD_2_1',dim(ST_15M_AD_2_1@meta.data)[1])),levels = c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')))


options(repr.plot.width=28, repr.plot.height=14)
ggplot(na.omit(cas_region_vln_df), aes(x=region, y=score , fill = section)) + 
geom_boxplot(position=position_dodge(1), outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1))

options(repr.plot.width=28, repr.plot.height=5)
ggplot(na.omit(cas_region_vln_df), aes(x=subregion, y=score, fill = section)) + 
geom_boxplot(position=position_dodge(1), outlier.shape = NA) #+
#geom_jitter() #+
#stat_summary(fun=median, geom="smooth",aes(group=1))

#Public ST APP/PS1-Signature scores

#Liu S, Li H, Shen Y, Zhu W et al. Moxibustion improves hypothalamus Aqp4 polarization in APP/PS1 mice: Evidence from spatial transcriptomics. Front Aging Neurosci 2023;15:1069155. PMID: 36819717

APP_PS1 <- Load10X_Spatial(data.dir = 'Public_ST/APP_PS1',filename = 'filtered_feature_bc_matrix.h5')

DefaultAssay(APP_PS1) <- "Spatial"
APP_PS1 <- SCTransform(APP_PS1, assay = "Spatial", verbose = TRUE)
DefaultAssay(APP_PS1) <- 'SCT'
APP_PS1 <- RunPCA(object = APP_PS1, verbose = T)
APP_PS1 <- RunUMAP(APP_PS1, dims = 1:30, verbose = T)

DefaultAssay(APP_PS1) <- 'SCT'
APP_PS1_SigScore_all_deg <- Vision(APP_PS1, signatures = all_deg_sig_list, assay='SCT')
APP_PS1_SigScore_all_deg <- analyze(APP_PS1_SigScore_all_deg)
for (signature_set in colnames(APP_PS1_SigScore_all_deg@SigScores)) {
  
  sigScores <- APP_PS1_SigScore_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  APP_PS1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(APP_PS1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=7, repr.plot.height=7)

min_aggr_all_deg <- min(min(APP_PS1$all_deg_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(APP_PS1$all_deg_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(APP_PS1, features = 'all_deg_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )

DefaultAssay(APP_PS1) <- 'SCT'
APP_PS1_SigScore_region_all_deg <- Vision(APP_PS1, signatures = region_common_all_deg_sig_list, assay='SCT')
APP_PS1_SigScore_region_all_deg <- analyze(APP_PS1_SigScore_region_all_deg)
for (signature_set in colnames(APP_PS1_SigScore_region_all_deg@SigScores)) {
  
  sigScores <- APP_PS1_SigScore_region_all_deg@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  APP_PS1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(APP_PS1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=7, repr.plot.height=7)

min_aggr_all_deg <- min(min(APP_PS1$region_common_all_deg_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(APP_PS1$region_common_all_deg_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(APP_PS1, features = 'region_common_all_deg_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )					 

DefaultAssay(APP_PS1) <- 'SCT'
APP_PS1_SigScore_cas <- Vision(APP_PS1, signatures = CAS_sig_list, assay='SCT')
APP_PS1_SigScore_cas <- analyze(APP_PS1_SigScore_cas)
for (signature_set in colnames(APP_PS1_SigScore_cas@SigScores)) {
  
  sigScores <- APP_PS1_SigScore_cas@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  APP_PS1@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(APP_PS1@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

options(repr.plot.width=7, repr.plot.height=7)

min_aggr_all_deg <- min(min(APP_PS1$CAS_sig,na.rm = TRUE))

max_aggr_all_deg <- max(max(APP_PS1$CAS_sig,na.rm = TRUE))

wrap_plots(
     SpatialFeaturePlot(APP_PS1, features = 'CAS_sig')+ 
     scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(min_aggr_all_deg,max_aggr_all_deg))
     )





