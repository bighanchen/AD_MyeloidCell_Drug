# import packages
library(Seurat)
library(tidyr)
library(dplyr)
library(plyr)
library(edgeR)
library(VISION)
library(stringr)
library(ggplot2

#Find the differentially expressed genes between AD and WT brains
st_list <- c('ST_3M_AD_1_1','ST_3M_AD_1_2','ST_6M_AD_2_1','ST_6M_AD_2_2','ST_15M_AD_2_1','ST_15M_AD_2_2',
             'ST_3M_WT_1_1','ST_3M_WT_1_2','ST_6M_WT_2_1','ST_6M_WT_2_2','ST_15M_WT_2_1','ST_15M_WT_2_2')

st_expr_list <- list()
             
for (samp in st_list){
    
    st_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(st_meta) <- c("barcode","subregion")
    st_meta <- st_meta[st_meta$subregion != 'VL',]
    
    
    st <- Load10X_Spatial(data.dir = samp,filename = 'filtered_feature_bc_matrix.h5')
    expr <-as.data.frame(st@assays$Spatial@counts[,st_meta$barcode])
         
    colnames(expr) <- paste(colnames(expr),strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep="_")
    expr$gene <- rownames(expr)
    st_expr_list[[samp]] <- expr
}

ad_wt_expr <- join_all(st_expr_list,by = "gene",type = "inner")

rownames(ad_wt_expr) <- ad_wt_expr$gene

ad_wt_expr <- ad_wt_expr[ , -which(colnames(ad_wt_expr) %in% c("gene"))]

status_deg <- c(rep("AD",dim(st_expr_list[['ST_3M_AD_1_1']])[2]+dim(st_expr_list[['ST_3M_AD_1_2']])[2]+dim(st_expr_list[['ST_6M_AD_2_1']])[2]+dim(st_expr_list[['ST_6M_AD_2_2']])[2]+dim(st_expr_list[['ST_15M_AD_2_1']])[2]+dim(st_expr_list[['ST_15M_AD_2_2']])[2]-6),
                rep("WT",dim(st_expr_list[['ST_3M_WT_1_1']])[2]+dim(st_expr_list[['ST_3M_WT_1_2']])[2]+dim(st_expr_list[['ST_6M_WT_2_1']])[2]+dim(st_expr_list[['ST_6M_WT_2_2']])[2]+dim(st_expr_list[['ST_15M_WT_2_1']])[2]+dim(st_expr_list[['ST_15M_WT_2_2']])[2]-6)
                )
status_deg <-factor(status_deg,levels=c("WT","AD"))
age_deg <- c(rep("3m",dim(st_expr_list[['ST_3M_AD_1_1']])[2]+dim(st_expr_list[['ST_3M_AD_1_2']])[2]-2),
             rep("6m",dim(st_expr_list[['ST_6M_AD_2_1']])[2]+dim(st_expr_list[['ST_6M_AD_2_2']])[2]-2),
             rep("15m",dim(st_expr_list[['ST_15M_AD_2_1']])[2]+dim(st_expr_list[['ST_15M_AD_2_2']])[2]-2),
             rep("3m",dim(st_expr_list[['ST_3M_WT_1_1']])[2]+dim(st_expr_list[['ST_3M_WT_1_2']])[2]-2),
             rep("6m",dim(st_expr_list[['ST_6M_WT_2_1']])[2]+dim(st_expr_list[['ST_6M_WT_2_2']])[2]-2),
             rep("15m",dim(st_expr_list[['ST_15M_WT_2_1']])[2]+dim(st_expr_list[['ST_15M_WT_2_2']])[2]-2)
             )
age_deg <- factor(age_deg,levels=c("3m","6m","15m"))

group_deg=factor(paste(status_deg,age_deg,sep="."))


ad_wt_expr_DGEList <- DGEList(counts=as.matrix(ad_wt_expr),group = group_deg)

keep <- rowSums(cpm(ad_wt_expr_DGEList) > 0 ) > 0
ad_wt_expr_DGEList <- ad_wt_expr_DGEList[keep, , keep.lib.sizes=FALSE]
ad_wt_expr_DGEList <- calcNormFactors(ad_wt_expr_DGEList)

design_deg <- model.matrix(~0+group, data=ad_wt_expr_DGEList$samples)


disp_deg <- estimateDisp(ad_wt_expr_DGEList, design_deg, robust = TRUE)
fit_deg <- glmQLFit(disp_deg, design_deg, robust = TRUE) 
qlf_AD15mVSWT15m <- glmQLFTest(fit_deg, contrast=c(1,0,0,-1,0,0))
qlf_AD6mVSWT6m <- glmQLFTest(fit_deg, contrast=c(0,0,1,0,0,-1))
qlf_AD3mVSWT3m <- glmQLFTest(fit_deg, contrast=c(0,1,0,0,-1,0))


toptags_deg_3m <- topTags(qlf_AD3mVSWT3m, n=Inf, adjust.method = "BH", p.value = 0.05)
toptags_deg_3m_df <- toptags_deg_3m$table
toptags_deg_3m_df$gene <- rownames(toptags_deg_3m_df)
colnames(toptags_deg_3m_df) <- c("logFC_3m","logCPM_3m","F_3m","PValue_3m","FDR_3m","gene")

toptags_deg_6m <- topTags(qlf_AD6mVSWT6m, n=Inf, adjust.method = "BH", p.value = 0.05)
toptags_deg_6m_df <- toptags_deg_6m$table
toptags_deg_6m_df$gene <- rownames(toptags_deg_6m_df)
colnames(toptags_deg_6m_df) <- c("logFC_6m","logCPM_6m","F_6m","PValue_6m","FDR_6m","gene")

toptags_deg_15m <- topTags(qlf_AD15mVSWT15m, n=Inf, adjust.method = "BH", p.value = 0.05)
toptags_deg_15m_df <- toptags_deg_15m$table
toptags_deg_15m_df$gene <- rownames(toptags_deg_15m_df)
colnames(toptags_deg_15m_df) <- c("logFC_15m","logCPM_15m","F_15m","PValue_15m","FDR_15m","gene")

toptags_deg_inner <- join_all(list(toptags_deg_3m_df,toptags_deg_6m_df,toptags_deg_15m_df),by = "gene",type = "inner")
toptags_deg_inner_cons <- toptags_deg_inner[(toptags_deg_inner$logFC_3m >0 & toptags_deg_inner$logFC_6m >0 & toptags_deg_inner$logFC_15m >0) | (toptags_deg_inner$logFC_3m <0 & toptags_deg_inner$logFC_6m <0 & toptags_deg_inner$logFC_15m <0),]

#Identify dysregulated genes in brain of 5XFAD mouse model

#Prepare data
st_ad_list <- c('ST_3M_AD_1_2','ST_6M_AD_2_1','ST_15M_AD_2_1')

spot_meta_list <- list()

spot_expr_list <- list()

for (samp in st_ad_list) {
    
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
    mean_meta_log <- data.frame(barcode = mean_abeta[,c('barcode')],logMean = log(mean_abeta[,c('average')]))
    # remove rows containing na, nan, inf, -inf
    mean_meta_log <- mean_meta_log[!is.na(mean_meta_log$logMean)&!is.nan(mean_meta_log$logMean)&!is.infinite(mean_meta_log$logMean),]
    
    region_meta <- read.csv(paste0('STAligner/',samp,'.csv'),header=TRUE)
    colnames(region_meta) <- c("barcode","subregion")
    region_meta <- region_meta[region_meta$subregion != 'VL',]
    
    spot_meta <- inner_join(mean_meta_log,region_meta,by = 'barcode')
    
    # get subset of sct data matrix
    spot_expr_mtx <- st@assays$Spatial@counts[,spot_meta$barcode]
    
    spot_expr <- as.data.frame(spot_expr_mtx)
    
    colnames(spot_expr) <- paste(colnames(spot_expr), strsplit(samp,'_')[[1]][2],sep = '_')
    
    spot_expr$gene <- rownames(spot_expr)
    
    spot_expr_list[[samp]] <- spot_expr
    
    spot_meta$barcode <- paste(spot_meta$barcode,strsplit(samp,'_')[[1]][2],sep = '_')
    
    samp_age <- as.integer(strsplit(strsplit(samp,'_')[[1]][2],'M')[[1]][1])
    
    spot_meta$age <- rep(samp_age,dim(spot_meta)[1])
    
    spot_meta_list[[samp]] <- spot_meta
} 

spot_expr_joined <- join_all(spot_expr_list,by = "gene",type = "inner")
rownames(spot_expr_joined) <- spot_expr_joined$gene
spot_expr_joined <- spot_expr_joined[,-which(names(spot_expr_joined) == "gene")]

spot_meta_joined <- rbind(spot_meta_list[['ST_3M_AD_1_2']],
                          spot_meta_list[['ST_6M_AD_2_1']],                               
                          spot_meta_list[['ST_15M_AD_2_1']])
						  
#Full model containing abeta and age
abeta <- spot_meta_joined$logMean
age <- c(rep(3,dim(spot_meta_list[['ST_3M_AD_1_2']])[1]),rep(6,dim(spot_meta_list[['ST_6M_AD_2_1']])[1]),rep(15,dim(spot_meta_list[['ST_15M_AD_2_1']])[1]))

# find plaque associated genes
spot_DGEList <- DGEList(counts=spot_expr_joined)

#keep <- filterByExpr(subspot_DGEList,design_pag)
keep <- rowSums(cpm(spot_DGEList) > 0 ) > 0
spot_DGEList <- spot_DGEList[keep, , keep.lib.sizes=FALSE]
spot_DGEList <- calcNormFactors(spot_DGEList)

design <- model.matrix(~abeta+age+abeta:age)
disp <- estimateDisp(spot_DGEList, design, robust = TRUE)
fit <- glmQLFit(disp, design, robust = TRUE)

qlf_abeta <- glmQLFTest(fit,coef = 2)
toptags_abeta <- topTags(qlf_abeta, n=Inf, adjust.method = "BH", p.value = 0.05)
qlf_age <- glmQLFTest(fit,coef = 3)
toptags_age <- topTags(qlf_age, n=Inf, adjust.method = "BH", p.value = 0.05)
qlf_abeta_age <- glmQLFTest(fit,coef = 4)
toptags_abeta_age <- topTags(qlf_abeta_age, n=Inf, adjust.method = "BH", p.value = 0.05)
qlf_all <- glmQLFTest(fit,coef = 2:4)
toptags_all <- topTags(qlf_all, n=Inf, adjust.method = "BH", p.value = 0.05)

save(toptags_abeta, file='ST_DysregulatedGenes_Global/toptags_abeta.bin')
save(toptags_age, file='ST_DysregulatedGenes_Global/toptags_age.bin')
save(toptags_abeta_age, file='ST_DysregulatedGenes_Global/toptags_abeta_age.bin')
save(toptags_all, file='ST_DysregulatedGenes_Global/toptags_all.bin')

#Genes affected by both abeta and aging
toptags_all_df <- toptags_all$table
toptags_all_df$gene <- rownames(toptags_all_df)
toptags_all_deg_cons <- inner_join(toptags_all_df,toptags_deg_inner_cons,by = "gene")

toptags_all_deg_cons_1e5 <- toptags_all_deg_cons[toptags_all_deg_cons$FDR<0.00001 & toptags_all_deg_cons$FDR_3m<0.00001 & toptags_all_deg_cons$FDR_6m<0.00001 & toptags_all_deg_cons$FDR_15m<0.00001,]
all_deg_sig_vec <- sign(toptags_all_deg_cons_1e5$logFC_3m)
names(all_deg_sig_vec) <- toptags_all_deg_cons_1e5$gene
all_deg_sig_vec <- all_deg_sig_vec[!(duplicated(names(all_deg_sig_vec)) | duplicated(names(all_deg_sig_vec),fromLast=TRUE))]
all_deg_sig_list <- c()
all_deg_sig <- createGeneSignature(name = "all_deg_sig", sigData = all_deg_sig_vec)
all_deg_sig_list <- c(all_deg_sig_list,all_deg_sig)
save(all_deg_sig_list, file='Signature_repository/all_deg_sig_list.bin')
write.csv(toptags_all_deg_cons[toptags_all_deg_cons$FDR<0.00001,],'Signature_repository/all_global.csv',quote = FALSE)
write.csv(toptags_all_deg_cons[toptags_all_deg_cons$FDR<0.00001 & toptags_all_deg_cons$FDR_3m<0.00001 & toptags_all_deg_cons$FDR_6m<0.00001 & toptags_all_deg_cons$FDR_15m<0.00001,],'Signature_repository/all_deg_global.csv',quote = FALSE)

#Genes affected by abeta, no matter other factors
toptags_abeta_df <- toptags_abeta$table
toptags_abeta_df$gene <- rownames(toptags_abeta_df)
toptags_abeta_deg_cons <- inner_join(toptags_abeta_df,toptags_deg_inner_cons,by = "gene")

toptags_abeta_deg_cons_1e5 <- toptags_abeta_deg_cons[toptags_abeta_deg_cons$FDR<0.00001 & toptags_abeta_deg_cons$FDR_3m<0.00001 & toptags_abeta_deg_cons$FDR_6m<0.00001 & toptags_abeta_deg_cons$FDR_15m<0.00001,]
abeta_deg_sig_vec <- sign(toptags_abeta_deg_cons_1e5$logFC_3m)
names(abeta_deg_sig_vec) <- toptags_abeta_deg_cons_1e5$gene
abeta_deg_sig_vec <- abeta_deg_sig_vec[!(duplicated(names(abeta_deg_sig_vec)) | duplicated(names(abeta_deg_sig_vec),fromLast=TRUE))]
abeta_deg_sig_list <- c()
abeta_deg_sig <- createGeneSignature(name = "abeta_deg_sig", sigData = abeta_deg_sig_vec)
abeta_deg_sig_list <- c(abeta_deg_sig_list,abeta_deg_sig)
save(abeta_deg_sig_list, file='Signature_repository/abeta_deg_sig_list.bin')

write.csv(toptags_abeta_deg_cons_1e5,'ST_DysregulatedGenes_Global/toptags_abeta_deg_cons_1e5.csv',quote=FALSE)

#Genes affected by age, no matter other factors
toptags_age_df <- toptags_age$table
toptags_age_df$gene <- rownames(toptags_age_df)
toptags_age_deg_cons <- inner_join(toptags_age_df,toptags_deg_inner_cons,by = "gene")

toptags_age_deg_cons_1e5 <- toptags_age_deg_cons[toptags_age_deg_cons$FDR<0.00001 & toptags_age_deg_cons$FDR_3m<0.00001 & toptags_age_deg_cons$FDR_6m<0.00001 & toptags_age_deg_cons$FDR_15m<0.00001,]
age_deg_sig_vec <- sign(toptags_age_deg_cons_1e5$logFC_3m)
names(age_deg_sig_vec) <- toptags_age_deg_cons_1e5$gene
age_deg_sig_vec <- age_deg_sig_vec[!(duplicated(names(age_deg_sig_vec)) | duplicated(names(age_deg_sig_vec),fromLast=TRUE))]
age_deg_sig_list <- c()
age_deg_sig <- createGeneSignature(name = "age_deg_sig", sigData = age_deg_sig_vec)
age_deg_sig_list <- c(age_deg_sig_list,age_deg_sig)
save(age_deg_sig_list, file='Signature_repository/age_deg_sig_list.bin')
write.csv(toptags_age_deg_cons_1e5,'ST_DysregulatedGenes_Global/toptags_age_deg_cons_1e5.csv',quote=FALSE)

#Genes affected by abeta:age interaction
toptags_abeta_age_df <- toptags_abeta_age$table
toptags_abeta_age_df$gene <- rownames(toptags_abeta_age_df)
toptags_abeta_age_deg_cons <- inner_join(toptags_abeta_age_df,toptags_deg_inner_cons,by = "gene")

toptags_abeta_age_deg_cons_1e5 <- toptags_abeta_age_deg_cons[toptags_abeta_age_deg_cons$FDR<0.00001 & toptags_abeta_age_deg_cons$FDR_3m<0.00001 & toptags_abeta_age_deg_cons$FDR_6m<0.00001 & toptags_abeta_age_deg_cons$FDR_15m<0.00001,]
abeta_age_deg_sig_vec <- sign(toptags_abeta_age_deg_cons_1e5$logFC_3m)
names(abeta_age_deg_sig_vec) <- toptags_abeta_age_deg_cons_1e5$gene
abeta_age_deg_sig_vec <- abeta_age_deg_sig_vec[!(duplicated(names(abeta_age_deg_sig_vec)) | duplicated(names(abeta_age_deg_sig_vec),fromLast=TRUE))]
abeta_age_deg_sig_list <- c()
abeta_age_deg_sig <- createGeneSignature(name = "abeta_age_deg_sig", sigData = abeta_age_deg_sig_vec)
abeta_age_deg_sig_list <- c(abeta_age_deg_sig_list,abeta_age_deg_sig)
save(abeta_age_deg_sig_list, file='Signature_repository/abeta_age_deg_sig_list.bin')
write.csv(toptags_abeta_age_deg_cons_1e5,'ST_DysregulatedGenes_Global/toptags_abeta_age_deg_cons_1e5.csv',quote=FALSE)

#Genes affected by abeta only
toptags_abeta_only <- anti_join(toptags_abeta_df, toptags_age_df, by = "gene") %>%
                      anti_join(toptags_abeta_age_df, by = "gene")
toptags_abeta_only_deg_cons <- inner_join(toptags_abeta_only,toptags_deg_inner_cons,by = "gene")
toptags_abeta_only_deg_cons_1e5 <- toptags_abeta_only_deg_cons[toptags_abeta_only_deg_cons$FDR<0.00001 & toptags_abeta_only_deg_cons$FDR_3m<0.00001 & toptags_abeta_only_deg_cons$FDR_6m<0.00001 & toptags_abeta_only_deg_cons$FDR_15m<0.00001,]
abeta_only_deg_sig_vec <- sign(toptags_abeta_only_deg_cons_1e5$logFC_3m)
names(abeta_only_deg_sig_vec) <- toptags_abeta_only_deg_cons_1e5$gene
abeta_only_deg_sig_vec <- abeta_only_deg_sig_vec[!(duplicated(names(abeta_only_deg_sig_vec)) | duplicated(names(abeta_only_deg_sig_vec),fromLast=TRUE))]
abeta_only_deg_sig_list <- c()
abeta_only_deg_sig <- createGeneSignature(name = "abeta_only_deg_sig", sigData = abeta_only_deg_sig_vec)
abeta_only_deg_sig_list <- c(abeta_only_deg_sig_list,abeta_only_deg_sig)
save(abeta_only_deg_sig_list, file='Signature_repository/abeta_only_deg_sig_list.bin')

#Genes affected by age only
toptags_age_only <- anti_join(toptags_age_df, toptags_abeta_df, by = "gene") %>%
             anti_join(toptags_abeta_age_df, by = "gene")
toptags_age_only_deg_cons <- inner_join(toptags_age_only,toptags_deg_inner_cons,by = "gene")

toptags_age_only_deg_cons_1e5 <- toptags_age_only_deg_cons[toptags_age_only_deg_cons$FDR<0.00001 & toptags_age_only_deg_cons$FDR_3m<0.00001 & toptags_age_only_deg_cons$FDR_6m<0.00001 & toptags_age_only_deg_cons$FDR_15m<0.00001,]
age_only_deg_sig_vec <- sign(toptags_age_only_deg_cons_1e5$logFC_3m)
names(age_only_deg_sig_vec) <- toptags_age_only_deg_cons_1e5$gene
age_only_deg_sig_vec <- age_only_deg_sig_vec[!(duplicated(names(age_only_deg_sig_vec)) | duplicated(names(age_only_deg_sig_vec),fromLast=TRUE))]
age_only_deg_sig_list <- c()
age_only_deg_sig <- createGeneSignature(name = "age_only_deg_sig", sigData = age_only_deg_sig_vec)
age_only_deg_sig_list <- c(age_only_deg_sig_list,age_only_deg_sig)
save(age_only_deg_sig_list, file='Signature_repository/age_only_deg_sig_list.bin')



