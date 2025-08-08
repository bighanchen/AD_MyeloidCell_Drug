# author: Xiaowei Chen (chenxiaowei@ibp.ac.cn)
# import packages
library(Seurat)
library(tidyr)
library(dplyr)
library(plyr)
library(edgeR)
library(VISION)
library(ggplot2)
library(VISION)

#Find differentially expressed genes between AD and WT brains
region_anno <- read.csv('STAligner/region.csv',header = FALSE)

colnames(region_anno) <- c('subregion','region')

st_3m_list <- c('ST_3M_AD_1_1','ST_3M_AD_1_2','ST_3M_WT_1_1','ST_3M_WT_1_2')
st_3m_expr_list <- list()
st_3m_meta_list <- list()

for (samp in st_3m_list) {
    st_3m_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(st_3m_meta) <- c("barcode","subregion")    
    st_3m_meta$region <- mapvalues(x = st_3m_meta$subregion, from = region_anno$subregion, to = region_anno$region)
    st_3m_meta <- st_3m_meta[st_3m_meta$subregion != 'VL',]
        
    st <- Load10X_Spatial(data.dir = samp,filename = "filtered_feature_bc_matrix.h5")
    expr <-as.data.frame(st@assays$Spatial@counts[,st_3m_meta$barcode])
    colnames(expr) <- paste(colnames(expr),strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep="_")
    expr$gene <- rownames(expr)
    st_3m_expr_list[[samp]] <- expr
    
    st_3m_meta$barcode <- paste(st_3m_meta$barcode,strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    samp_group <- strsplit(samp,'_')[[1]][3]
    st_3m_meta$group <- rep(samp_group,dim(st_3m_meta)[1])
    st_3m_meta_list[[samp]] <- st_3m_meta
}

st_3m_expr_joined <- join_all(st_3m_expr_list,by = "gene",type = "inner")

rownames(st_3m_expr_joined) <- st_3m_expr_joined$gene

st_3m_expr_joined <- st_3m_expr_joined[ , -which(colnames(st_3m_expr_joined) %in% c("gene"))]

st_3m_meta_joined <- rbind(st_3m_meta_list[['ST_3M_AD_1_1']],
                               st_3m_meta_list[['ST_3M_AD_1_2']],
                               st_3m_meta_list[['ST_3M_WT_1_1']],
                               st_3m_meta_list[['ST_3M_WT_1_2']])

deg_3m_list <- list()
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    group <- st_3m_meta_joined[st_3m_meta_joined$region == region, 'group']
    spot_DGEList <- DGEList(counts=st_3m_expr_joined[,st_3m_meta_joined[st_3m_meta_joined$region == region, 'barcode']],group = group)
    keep <- rowSums(cpm(spot_DGEList) > 0 ) > 0
    spot_DGEList <- spot_DGEList[keep, , keep.lib.sizes=FALSE]
    spot_DGEList <- calcNormFactors(spot_DGEList)
    
    design <- model.matrix(~0+group, data=spot_DGEList$samples)    
    disp <- estimateDisp(spot_DGEList, design, robust = TRUE)
    fit <- glmQLFit(disp, design, robust = TRUE) 
    
    qlf_ADvsWT <- glmQLFTest(fit, contrast=c(1,-1))
    toptags <- topTags(qlf_ADvsWT, n=Inf, adjust.method = "BH", p.value = 0.05)
    toptags_df <- toptags$table
    toptags_df$gene <- rownames(toptags_df)
    deg_3m_list[[region]] <- toptags_df
}

st_6m_list <- c('ST_6M_AD_2_1','ST_6M_AD_2_2','ST_6M_WT_2_1','ST_6M_WT_2_2')
st_6m_expr_list <- list()
st_6m_meta_list <- list()

for (samp in st_6m_list) {
    st_6m_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(st_6m_meta) <- c("barcode","subregion")    
    st_6m_meta$region <- mapvalues(x = st_6m_meta$subregion, from = region_anno$subregion, to = region_anno$region)
    st_6m_meta <- st_6m_meta[st_6m_meta$subregion != 'VL',]
        
    st <- Load10X_Spatial(data.dir = samp,filename = "filtered_feature_bc_matrix.h5")
    expr <-as.data.frame(st@assays$Spatial@counts[,st_6m_meta$barcode])
    colnames(expr) <- paste(colnames(expr),strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep="_")
    expr$gene <- rownames(expr)
    st_6m_expr_list[[samp]] <- expr
    
    st_6m_meta$barcode <- paste(st_6m_meta$barcode,strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    samp_group <- strsplit(samp,'_')[[1]][3]
    st_6m_meta$group <- rep(samp_group,dim(st_6m_meta)[1])
    st_6m_meta_list[[samp]] <- st_6m_meta
}

st_6m_expr_joined <- join_all(st_6m_expr_list,by = "gene",type = "inner")

rownames(st_6m_expr_joined) <- st_6m_expr_joined$gene

st_6m_expr_joined <- st_6m_expr_joined[ , -which(colnames(st_6m_expr_joined) %in% c("gene"))]

st_6m_meta_joined <- rbind(st_6m_meta_list[['ST_6M_AD_2_1']],
                               st_6m_meta_list[['ST_6M_AD_2_2']],
                               st_6m_meta_list[['ST_6M_WT_2_1']],
                               st_6m_meta_list[['ST_6M_WT_2_2']])

deg_6m_list <- list()
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    group <- st_6m_meta_joined[st_6m_meta_joined$region == region, 'group']
    spot_DGEList <- DGEList(counts=st_6m_expr_joined[,st_6m_meta_joined[st_6m_meta_joined$region == region, 'barcode']],group = group)
    keep <- rowSums(cpm(spot_DGEList) > 0 ) > 0
    spot_DGEList <- spot_DGEList[keep, , keep.lib.sizes=FALSE]
    spot_DGEList <- calcNormFactors(spot_DGEList)
    
    design <- model.matrix(~0+group, data=spot_DGEList$samples)    
    disp <- estimateDisp(spot_DGEList, design, robust = TRUE)
    fit <- glmQLFit(disp, design, robust = TRUE) 
    
    qlf_ADvsWT <- glmQLFTest(fit, contrast=c(1,-1))
    toptags <- topTags(qlf_ADvsWT, n=Inf, adjust.method = "BH", p.value = 0.05)
    toptags_df <- toptags$table
    toptags_df$gene <- rownames(toptags_df)
    deg_6m_list[[region]] <- toptags_df
}


st_15m_list <- c('ST_15M_AD_2_1','ST_15M_AD_2_2','ST_15M_WT_2_1','ST_15M_WT_2_2')
st_15m_expr_list <- list()
st_15m_meta_list <- list()

for (samp in st_15m_list) {
    st_15m_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(st_15m_meta) <- c("barcode","subregion")    
    st_15m_meta$region <- mapvalues(x = st_15m_meta$subregion, from = region_anno$subregion, to = region_anno$region)
    st_15m_meta <- st_15m_meta[st_15m_meta$subregion != 'VL',]
    
    st <- Load10X_Spatial(data.dir = samp,filename = "filtered_feature_bc_matrix.h5")
    expr <-as.data.frame(st@assays$Spatial@counts[,st_15m_meta$barcode])
    colnames(expr) <- paste(colnames(expr),strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep="_")
    expr$gene <- rownames(expr)
    st_15m_expr_list[[samp]] <- expr
    
    st_15m_meta$barcode <- paste(st_15m_meta$barcode,strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    samp_group <- strsplit(samp,'_')[[1]][3]
    st_15m_meta$group <- rep(samp_group,dim(st_15m_meta)[1])
    st_15m_meta_list[[samp]] <- st_15m_meta
}

st_15m_expr_joined <- join_all(st_15m_expr_list,by = "gene",type = "inner")

rownames(st_15m_expr_joined) <- st_15m_expr_joined$gene

st_15m_expr_joined <- st_15m_expr_joined[ , -which(colnames(st_15m_expr_joined) %in% c("gene"))]

st_15m_meta_joined <- rbind(st_15m_meta_list[['ST_15M_AD_2_1']],
                               st_15m_meta_list[['ST_15M_AD_2_2']],
                               st_15m_meta_list[['ST_15M_WT_2_1']],
                               st_15m_meta_list[['ST_15M_WT_2_2']])
deg_15m_list <- list()
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    group <- st_15m_meta_joined[st_15m_meta_joined$region == region, 'group']
    spot_DGEList <- DGEList(counts=st_15m_expr_joined[,st_15m_meta_joined[st_15m_meta_joined$region == region, 'barcode']],group = group)
    keep <- rowSums(cpm(spot_DGEList) > 0 ) > 0
    spot_DGEList <- spot_DGEList[keep, , keep.lib.sizes=FALSE]
    spot_DGEList <- calcNormFactors(spot_DGEList)
    
    design <- model.matrix(~0+group, data=spot_DGEList$samples)    
    disp <- estimateDisp(spot_DGEList, design, robust = TRUE)
    fit <- glmQLFit(disp, design, robust = TRUE) 
    
    qlf_ADvsWT <- glmQLFTest(fit, contrast=c(1,-1))
    toptags <- topTags(qlf_ADvsWT, n=Inf, adjust.method = "BH", p.value = 0.05)
    toptags_df <- toptags$table
    toptags_df$gene <- rownames(toptags_df)
    deg_15m_list[[region]] <- toptags_df
}

save(deg_3m_list, file='ST_DysregulatedGenes_Regional/deg_3m_list.bin')
save(deg_6m_list, file='ST_DysregulatedGenes_Regional/deg_6m_list.bin')
save(deg_15m_list, file='ST_DysregulatedGenes_Regional/deg_15m_list.bin')

deg_inner_list <- list()
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    deg_3m <- deg_3m_list[[region]]
    colnames(deg_3m) <- c(paste0(colnames(deg_3m)[1:5],'_3m'),'gene')
    deg_6m <- deg_6m_list[[region]]
    colnames(deg_6m) <- c(paste0(colnames(deg_6m)[1:5],'_6m'),'gene')
    deg_15m <- deg_15m_list[[region]]
    colnames(deg_15m) <- c(paste0(colnames(deg_15m)[1:5],'_15m'),'gene')
    
    deg_inner <- inner_join(deg_3m,deg_6m,by = 'gene') %>% inner_join(deg_15m,by='gene')
    deg_inner <- deg_inner[(deg_inner$logFC_3m >0 & deg_inner$logFC_6m >0 & deg_inner$logFC_15m >0) | (deg_inner$logFC_3m <0 & deg_inner$logFC_6m <0 & deg_inner$logFC_15m <0),]
    deg_inner_list[[region]] <- deg_inner
}

#Identify dysregulated genes in each region of 5XFAD mouse model brain
#Prepare data
st_ad_list <- c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')

region_anno <- read.csv('STAligner/region.csv',header = FALSE)

colnames(region_anno) <- c('subregion','region')

spot_meta_list <- list()

spot_expr_list <- list()

for (samp in st_ad_list) {
    
    #i <- which(st_ad_list==samp)
    
    #load data
    st <- Load10X_Spatial(data.dir = samp,filename = 'filtered_feature_bc_matrix.h5')
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
    
    #get natural logarithms log(intensity mean)
    mean_abeta_log <- data.frame(barcode = mean_abeta[,c('barcode')],logMean = log(mean_abeta[,c('average')]))
    # remove rows containing na, nan, inf, -inf
    mean_abeta_log <- mean_abeta_log[!is.na(mean_abeta_log$logMean)&!is.nan(mean_abeta_log$logMean)&!is.infinite(mean_abeta_log$logMean),]
    
    STAligner <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(STAligner) <- c("barcode","subregion")
    
    STAligner$region <- mapvalues(x = STAligner$subregion, from = region_anno$subregion, to = region_anno$region)
    STAligner <- STAligner[STAligner$subregion != 'VL',]
    
    spot_meta <- inner_join(mean_abeta_log,STAligner,by = 'barcode')
    
    #spot_meta <- na.omit(spot_meta)
    
    # get subset of sct data matrix
    spot_expr_mtx <- st@assays$Spatial@counts[,spot_meta$barcode]
    
    spot_expr <- as.data.frame(spot_expr_mtx)
    
    colnames(spot_expr) <- paste(colnames(spot_expr), strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][5],sep = '_')
    
    spot_expr$gene <- rownames(spot_expr)
    
    spot_expr_list[[samp]] <- spot_expr
    
    spot_meta$barcode <- paste(spot_meta$barcode,strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][5],sep = '_')
    
    samp_age <- as.integer(strsplit(strsplit(samp,'_')[[1]][2],'M')[[1]][1])
    
    spot_meta$age <- rep(samp_age,dim(spot_meta)[1])
    
    spot_meta_list[[samp]] <- spot_meta
} 

spot_expr_joined <- join_all(spot_expr_list,by = "gene",type = "inner")
rownames(spot_expr_joined) <- spot_expr_joined$gene
spot_expr_joined <- spot_expr_joined[,-which(names(spot_expr_joined) == "gene")]

spot_meta_joined <- rbind(spot_meta_list[['ST_3M_AD_1_2']],spot_meta_list[['ST_6M_AD_2_1']],spot_meta_list[['ST_15M_AD_2_1']])

#Full model containing abeta and age 
toptags_abeta_list <- list()

toptags_age_list <- list()

toptags_abeta_age_list <- list()

toptags_all_list <- list()

for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    
    spot_abeta <- spot_meta_joined[spot_meta_joined$region == region, 'logMean']
    
    spot_age   <- spot_meta_joined[spot_meta_joined$region == region, 'age']
    
    spot_DGEList <- DGEList(counts=spot_expr_joined[,spot_meta_joined[spot_meta_joined$region == region, 'barcode']])
    
    keep <- rowSums(cpm(spot_DGEList) > 0 ) > 0
    spot_DGEList <- spot_DGEList[keep, , keep.lib.sizes=FALSE]
    spot_DGEList <- calcNormFactors(spot_DGEList)
    
    design <- model.matrix(~spot_abeta+spot_age+spot_abeta:spot_age)
    
    disp <- estimateDisp(spot_DGEList, design, robust = TRUE)
    
    fit <- glmQLFit(disp, design, robust = TRUE)
    
    qlf_abeta <- glmQLFTest(fit,coef = 2)    
    toptags_abeta <- topTags(qlf_abeta, n=Inf, adjust.method = "BH", p.value = 0.05)    
    toptags_abeta_list[[region]] <- toptags_abeta$table
    
    qlf_age <- glmQLFTest(fit,coef = 3)    
    toptags_age <- topTags(qlf_age, n=Inf, adjust.method = "BH", p.value = 0.05)    
    toptags_age_list[[region]] <- toptags_age$table
    
    qlf_abeta_age <- glmQLFTest(fit,coef = 4)    
    toptags_abeta_age <- topTags(qlf_abeta_age, n=Inf, adjust.method = "BH", p.value = 0.05)    
    toptags_abeta_age_list[[region]] <- toptags_abeta_age$table
    
    qlf_all <- glmQLFTest(fit,coef = 2:4)    
    toptags_all <- topTags(qlf_all, n=Inf, adjust.method = "BH", p.value = 0.05) 
    toptags_all_list[[region]] <- toptags_all$table    
}

save(toptags_abeta_list, file='ST_DysregulatedGenes_Regional/toptags_abeta_list.bin')
save(toptags_age_list, file='ST_DysregulatedGenes_Regional/toptags_age_list.bin')
save(toptags_abeta_age_list, file='ST_DysregulatedGenes_Regional/toptags_abeta_age_list.bin')
save(toptags_all_list, file='ST_DysregulatedGenes_Regional/toptags_all_list.bin')

toptags_all_deg_inner_list <- list()
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    toptags_all <- toptags_all_list[[region]]
    toptags_all$gene <- rownames(toptags_all)
    toptags_all_deg_inner <- inner_join(toptags_all,deg_inner_list[[region]],by = 'gene')
    toptags_all_deg_inner_list[[region]] <- toptags_all_deg_inner
}

region_hist <- c()
up_hist <-c()
down_hist <- c()
deg_hist <- c()
i=1
for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    region_hist[i] <- region
    up_hist[i] <- sum(toptags_all_deg_inner_list[[region]]$logFC_3m>0, na.rm = TRUE)
    down_hist[i] <- sum(toptags_all_deg_inner_list[[region]]$logFC_3m<0, na.rm = TRUE)
    deg_hist[i] <- up_hist[i]+ down_hist[i]
    i=i+1
}
region_all_deg_hist <- data.frame(region = region_hist, Up_regulated = up_hist, Down_regulated = down_hist, Total = deg_hist)
region_all_deg_hist <- region_all_deg_hist[order(region_all_deg_hist$Total,decreasing = TRUE),]
region_all_deg_hist$region <- factor(region_all_deg_hist$region,levels = region_all_deg_hist$region)
region_all_deg_hist <- region_all_deg_hist %>% pivot_longer(c('Up_regulated','Down_regulated'),names_to = 'group', values_to = 'count')
region_all_deg_hist$group <- factor(region_all_deg_hist$group,levels = c('Up_regulated','Down_regulated'))

up_regulated_list <-list()
down_regulated_list <- list()

for (region in c('hippocampal_region','fiber_tracts','cerebral_nuclei','thalamus','hypothalamus','cerebral_cortex')){
    up_gene <- toptags_all_deg_inner_list[[region]][toptags_all_deg_inner_list[[region]]$logFC_3m>0,'gene']
    up_region <- rep(1,sum(toptags_all_deg_inner_list[[region]]$logFC_3m>0))
    up_regulated <- data.frame(gene = up_gene, reg = up_region)
    colnames(up_regulated) <- c('gene',region)
    down_gene <- toptags_all_deg_inner_list[[region]][toptags_all_deg_inner_list[[region]]$logFC_3m<0,'gene']
    down_region <- rep(1,sum(toptags_all_deg_inner_list[[region]]$logFC_3m<0))
    down_regulated <- data.frame(gene = down_gene, reg = down_region)
    colnames(down_regulated) <- c('gene',region)
    
    up_regulated_list[[region]] <- up_regulated
    down_regulated_list[[region]] <- down_regulated

}

up_regulated_full_join <- join_all(up_regulated_list,by = 'gene', type = 'full')
down_regulated_full_join <- join_all(down_regulated_list,by = 'gene', type = 'full')
down_regulated_full_join_common <- down_regulated_full_join[rowSums(down_regulated_full_join[2:7],na.rm = TRUE)>3,]
up_regulated_full_join_common <- up_regulated_full_join[rowSums(up_regulated_full_join[2:7],na.rm = TRUE)>3,]

region_common_all_deg_sig_vec <- c(rep(1,dim(up_regulated_full_join_common)[1]),rep(-1,dim(down_regulated_full_join_common)[1]))
names(region_common_all_deg_sig_vec) <- c(up_regulated_full_join_common$gene,down_regulated_full_join_common$gene)
region_common_all_deg_sig_vec <- region_common_all_deg_sig_vec[!(duplicated(names(region_common_all_deg_sig_vec)) | duplicated(names(region_common_all_deg_sig_vec),fromLast=TRUE))]
region_common_all_deg_sig_list <- c()
region_common_all_deg_sig <- createGeneSignature(name = "region_common_all_deg_sig", sigData = region_common_all_deg_sig_vec)
region_common_all_deg_sig_list <- c(region_common_all_deg_sig_list,region_common_all_deg_sig)
save(region_common_all_deg_sig_list, file='Signature_repository/region_common_all_deg_sig_list.bin')

dysregulated_full_join_common <- rbind(up_regulated_full_join_common,down_regulated_full_join_common)
write.csv(dysregulated_full_join_common,'Signature_repository/dysregulated_gene_common_regions.csv',quote = FALSE)

