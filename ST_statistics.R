# author: Xiaowei Chen (chenxiaowei@ibp.ac.cn)
#import packages
library(Seurat)
library(tidyr)
library(dplyr)
library(plyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Prepare data
st_ad_list <- c('ST_3M_AD_1_1','ST_3M_AD_1_2','ST_6M_AD_2_1','ST_6M_AD_2_2','ST_15M_AD_2_1','ST_15M_AD_2_2')

abeta_meta_list <- list()

for (samp in st_ad_list) {
        
    #load data
    #st <- Load10X_Spatial(data.dir = samp,filename = 'filtered_feature_bc_matrix.h5')
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
    lower_abeta <- lower_abeta_intensity_position[lower_abeta_intensity_position$covered == 1,c('barcode','Mean','StdDev')]
    colnames(lower_abeta) <- c('barcode','mean_lower','stddev_lower')
    
    #read Abeta intensity measurement
    upper_abeta_intensity <- read.csv(paste0(samp,'/',samp,'_upper_abeta_intensity.csv'),header =TRUE)
    
    #read spot coordinates from Imaris
    upper_position <- read.csv(paste0(samp,'/',samp,'_upper_position.csv'),header =TRUE)
    #get Abeta intensity measurement for each spot
    upper_abeta_intensity_position <- full_join(upper_position,upper_abeta_intensity,by = 'ID')
    upper_abeta_intensity_position <- unite(upper_abeta_intensity_position, 'row_col', row, col, sep = '_', remove = FALSE)
    upper_abeta_intensity_position <- left_join(upper_abeta_intensity_position,tissue_positions,by = 'row_col')
    upper_abeta <- upper_abeta_intensity_position[upper_abeta_intensity_position$covered == 1,c('barcode','Mean','StdDev')]
    colnames(upper_abeta) <- c('barcode','mean_upper','stddev_upper')
    
    abeta_meta <- inner_join(lower_abeta,upper_abeta,by = 'barcode')
    abeta_meta <- na.omit(abeta_meta)
    abeta_meta$barcode <- paste(abeta_meta$barcode,strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    abeta_meta$group <- paste(strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    samp_age <- as.integer(strsplit(strsplit(samp,'_')[[1]][2],'M')[[1]][1])
    
    abeta_meta$age <- rep(samp_age,dim(abeta_meta)[1])
    abeta_meta_list[[samp]] <- abeta_meta
} 

st_list <- c('ST_3M_AD_1_1','ST_3M_AD_1_2','ST_6M_AD_2_1','ST_6M_AD_2_2','ST_15M_AD_2_1','ST_15M_AD_2_2',
             'ST_3M_WT_1_1','ST_3M_WT_1_2','ST_6M_WT_2_1','ST_6M_WT_2_2','ST_15M_WT_2_1','ST_15M_WT_2_2')

spot_meta_list <- list()

for (samp in st_list) {
        
    #load data
    st <- Load10X_Spatial(data.dir = samp,filename = 'filtered_feature_bc_matrix.h5')
         
    # get spot meta.data
    spot_meta <- st@meta.data
    spot_meta$barcode <- rownames(spot_meta)
    spot_meta$orig.ident <- NULL
    
    spot_meta$barcode <- paste(spot_meta$barcode,strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    
    spot_meta$status <- strsplit(samp,'_')[[1]][3]
    
    spot_meta$group <- paste(strsplit(samp,'_')[[1]][2],strsplit(samp,'_')[[1]][3],strsplit(samp,'_')[[1]][5],sep = '_')
    
    samp_age <- as.integer(strsplit(strsplit(samp,'_')[[1]][2],'M')[[1]][1])
    
    spot_meta$age <- rep(samp_age,dim(spot_meta)[1])
    
    spot_meta_list[[samp]] <- spot_meta
} 

abeta_meta_joined <- rbind(abeta_meta_list[['ST_3M_AD_1_1']],
                           abeta_meta_list[['ST_3M_AD_1_2']],
                           abeta_meta_list[['ST_6M_AD_2_1']],
                           abeta_meta_list[['ST_6M_AD_2_2']],
                           abeta_meta_list[['ST_15M_AD_2_1']],
                           abeta_meta_list[['ST_15M_AD_2_2']])

spot_meta_joined <- rbind(spot_meta_list[['ST_3M_AD_1_1']],
                          spot_meta_list[['ST_3M_AD_1_2']],
                          spot_meta_list[['ST_6M_AD_2_1']],
                          spot_meta_list[['ST_6M_AD_2_2']],
                          spot_meta_list[['ST_15M_AD_2_1']],
                          spot_meta_list[['ST_15M_AD_2_2']],
                          spot_meta_list[['ST_3M_WT_1_1']],
                          spot_meta_list[['ST_3M_WT_1_2']],
                          spot_meta_list[['ST_6M_WT_2_1']],
                          spot_meta_list[['ST_6M_WT_2_2']],
                          spot_meta_list[['ST_15M_WT_2_1']],
                          spot_meta_list[['ST_15M_WT_2_2']])
						  
# Spot statistics
spot_meta_joined$group <- factor(spot_meta_joined$group,levels = c('3M_AD_1','3M_AD_2','3M_WT_1','3M_WT_2','6M_AD_1','6M_AD_2','6M_WT_1','6M_WT_2','15M_AD_1','15M_AD_2','15M_WT_1','15M_WT_2'))						  
spot_count <- table(spot_meta_joined$group)
spot_count_df <- as.data.frame(spot_count)
colnames(spot_count_df) <- c('samp','spot_count')
spot_count_df$age <-factor(c('3m','3m','3m','3m','6m','6m','6m','6m','15m','15m','15m','15m'),levels = c('3m','6m','15m'))
spot_count_df$status <-factor(c('AD_1','AD_2','WT_1','WT_2','AD_1','AD_2','WT_1','WT_2','AD_1','AD_2','WT_1','WT_2'),levels = c('AD_1','AD_2','WT_1','WT_2'))
ggplot(spot_count_df,mapping = aes(x = age,y = spot_count,fill= status))+
geom_bar(stat = 'identity', position ='dodge')

ncount_status_p <- ggplot(spot_meta_joined, aes(x= status, y=nCount_Spatial))+
  geom_violin(aes(fill = factor(status)))
ncount_age_p <- ggplot(spot_meta_joined, aes(x= factor(age), y=nCount_Spatial))+
  geom_violin(aes(fill = factor(age)))
nfeature_status_p <- ggplot(spot_meta_joined, aes(x= status, y=nFeature_Spatial))+
  geom_violin(aes(fill = factor(status)))
nfeature_age_p <- ggplot(spot_meta_joined, aes(x= factor(age), y=nFeature_Spatial))+
  geom_violin(aes(fill = factor(age)))  
options(repr.plot.width=9, repr.plot.height=6)
(ncount_status_p | ncount_age_p) / (nfeature_status_p | nfeature_age_p) 

# Abeta statistics
abeta_count <- table(abeta_meta_joined$group)
abeta_count_df <- as.data.frame(abeta_count)
colnames(abeta_count_df) <- c('samp','abeta_count')
abeta_spot_count <- left_join(abeta_count_df,spot_count_df,by = 'samp')

abeta_spot_count_long <- abeta_spot_count %>% pivot_longer(2:3,names_to = "type",values_to = "count")
ggplot(abeta_spot_count_long, aes(fill= type,x= status,y=count)) +  
geom_bar(position = "dodge", stat = "identity") + 
facet_grid(.~age)

#Abeta distribution and correlation
abeta_meta_mean <- abeta_meta_joined[,c('barcode','mean_lower','mean_upper','group')] %>% pivot_longer(c('mean_lower','mean_upper'),names_to = 'mean', values_to = 'count')
abeta_meta_mean$group <- factor(abeta_meta_mean$group,levels = c('3M_AD_1','3M_AD_2','6M_AD_1','6M_AD_2','15M_AD_1','15M_AD_2'))
options(repr.plot.width=18, repr.plot.height=3)
ggplot(abeta_meta_mean, aes(x= factor(mean), y=log(count)))+
geom_violin(aes(fill = factor(mean)))+
facet_grid(.~group)

abeta_meta_stddev <- abeta_meta_joined[,c('barcode','stddev_lower','stddev_upper','group')] %>% pivot_longer(c('stddev_lower','stddev_upper'),names_to = 'stddev', values_to = 'count')
abeta_meta_stddev$group <- factor(abeta_meta_stddev$group,levels = c('3M_AD_1','3M_AD_2','6M_AD_1','6M_AD_2','15M_AD_1','15M_AD_2'))
options(repr.plot.width=18, repr.plot.height=3)
ggplot(abeta_meta_stddev, aes(x= factor(stddev), y=log(count)))+
geom_violin(aes(fill = factor(stddev)))+
facet_grid(.~group)

samp_list <- c('3M_AD_1','3M_AD_2','6M_AD_1','6M_AD_2','15M_AD_1','15M_AD_2')
mean_corre_plot_list <- list()
stddev_corre_plot_list <- list()
for (samp in samp_list) {
    pm <- ggplot(data = abeta_meta_joined[abeta_meta_joined$group==samp,], mapping = aes(x = mean_lower, y = mean_upper))+
    geom_point()+geom_smooth(method="lm")+
    scale_x_log10(expand=c(0,0),limits = c(min(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_upper > 0,'mean_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_lower > 0,'mean_lower']), max(abeta_meta_joined[abeta_meta_joined$group==samp  & abeta_meta_joined$mean_upper > 0,'mean_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_lower > 0,'mean_lower']))) + 
    scale_y_log10(expand=c(0,0),limits = c(min(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_upper > 0,'mean_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_lower > 0,'mean_lower']), max(abeta_meta_joined[abeta_meta_joined$group==samp  & abeta_meta_joined$mean_upper > 0,'mean_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$mean_lower > 0,'mean_lower']))) + 
    #scale_x_log10(expand=c(0,0),limits = c(0.001,175)) + 
    #scale_y_log10(expand=c(0,0),limits = c(0.001,175)) +

    coord_fixed(ratio=1)
    mean_corre_plot_list[[samp]] <- pm
    
    ps <- ggplot(data = abeta_meta_joined[abeta_meta_joined$group==samp,], mapping = aes(x = stddev_lower, y = stddev_upper))+
    geom_point()+geom_smooth(method="lm")+
    scale_x_log10(expand=c(0,0),limits = c(min(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_upper > 0,'stddev_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_lower > 0,'stddev_lower']), max(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_upper > 0,'stddev_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_lower > 0,'stddev_lower']))) + 
    scale_y_log10(expand=c(0,0),limits = c(min(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_upper > 0,'stddev_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_lower > 0,'stddev_lower']), max(abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_upper > 0,'stddev_upper'],abeta_meta_joined[abeta_meta_joined$group==samp & abeta_meta_joined$stddev_lower > 0,'stddev_lower']))) + 
    #scale_x_log10(expand=c(0,0),limits = c(0.001,175)) + 
    #scale_y_log10(expand=c(0,0),limits = c(0.001,175)) +
    coord_fixed(ratio=1)
    stddev_corre_plot_list[[samp]] <- ps
}
options(repr.plot.width=18, repr.plot.height=3)
mean_corre_plot_list[[1]]|mean_corre_plot_list[[2]]|mean_corre_plot_list[[3]]|mean_corre_plot_list[[4]]|mean_corre_plot_list[[5]]|mean_corre_plot_list[[6]]
options(repr.plot.width=18, repr.plot.height=3)
stddev_corre_plot_list[[1]]|stddev_corre_plot_list[[2]]|stddev_corre_plot_list[[3]]|stddev_corre_plot_list[[4]]|stddev_corre_plot_list[[5]]|stddev_corre_plot_list[[6]]

#Abeta feature plots
ST_3M_AD_1_2 <- Load10X_Spatial(data.dir = "./ST_3M_AD_1_2/",filename = "filtered_feature_bc_matrix.h5")
ST_3M_AD_1_2_metadata <- data.frame(barcode=rownames(ST_3M_AD_1_2@meta.data))
ST_3M_AD_1_2_abeta <- abeta_meta_list[['ST_3M_AD_1_2']]
for (i in 1:dim(ST_3M_AD_1_2_abeta)[1]){
    ST_3M_AD_1_2_abeta[i,1] <- strsplit(ST_3M_AD_1_2_abeta[i,1],'_')[[1]][1]
}
ST_3M_AD_1_2_metadata_abeta <- left_join(ST_3M_AD_1_2_metadata,ST_3M_AD_1_2_abeta,by = "barcode")
rownames(ST_3M_AD_1_2_metadata_abeta)<-ST_3M_AD_1_2_metadata_abeta$barcode
ST_3M_AD_1_2_metadata_abeta <- ST_3M_AD_1_2_metadata_abeta[,-1]

ST_3M_AD_1_2 <- AddMetaData(
  object = ST_3M_AD_1_2,
  metadata = ST_3M_AD_1_2_metadata_abeta[,c('mean_lower','stddev_lower','mean_upper','stddev_upper')]
)

max_limit <- max(na.omit(ST_3M_AD_1_2@meta.data$mean_lower),na.omit(ST_3M_AD_1_2@meta.data$mean_upper))
ST_3M_AD_1_2_mean_lower_p <- SpatialFeaturePlot(ST_3M_AD_1_2, features = 'mean_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_3M_AD_1_2_mean_upper_p <- SpatialFeaturePlot(ST_3M_AD_1_2, features = 'mean_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)
ST_3M_AD_1_2_mean_lower_p | ST_3M_AD_1_2_mean_upper_p

max_limit <- max(na.omit(ST_3M_AD_1_2@meta.data$stddev_lower),na.omit(ST_3M_AD_1_2@meta.data$stddev_upper))
ST_3M_AD_1_2_stddev_lower_p <- SpatialFeaturePlot(ST_3M_AD_1_2, features = 'stddev_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_3M_AD_1_2_stddev_upper_p <- SpatialFeaturePlot(ST_3M_AD_1_2, features = 'stddev_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)
ST_3M_AD_1_2_stddev_lower_p | ST_3M_AD_1_2_stddev_upper_p

# 6 months
ST_6M_AD_2_1 <- Load10X_Spatial(data.dir = "./ST_6M_AD_2_1/",filename = "filtered_feature_bc_matrix.h5")
ST_6M_AD_2_1_metadata <- data.frame(barcode=rownames(ST_6M_AD_2_1@meta.data))
ST_6M_AD_2_1_abeta <- abeta_meta_list[['ST_6M_AD_2_1']]
for (i in 1:dim(ST_6M_AD_2_1_abeta)[1]){
    ST_6M_AD_2_1_abeta[i,1] <- strsplit(ST_6M_AD_2_1_abeta[i,1],'_')[[1]][1]
}
ST_6M_AD_2_1_metadata_abeta <- left_join(ST_6M_AD_2_1_metadata,ST_6M_AD_2_1_abeta,by = "barcode")
rownames(ST_6M_AD_2_1_metadata_abeta)<-ST_6M_AD_2_1_metadata_abeta$barcode
ST_6M_AD_2_1_metadata_abeta <- ST_6M_AD_2_1_metadata_abeta[,-1]

ST_6M_AD_2_1 <- AddMetaData(
  object = ST_6M_AD_2_1,
  metadata = ST_6M_AD_2_1_metadata_abeta[,c('mean_lower','stddev_lower','mean_upper','stddev_upper')]
)

max_limit <- max(na.omit(ST_6M_AD_2_1@meta.data$mean_lower),na.omit(ST_6M_AD_2_1@meta.data$mean_upper))
ST_6M_AD_2_1_mean_lower_p <- SpatialFeaturePlot(ST_6M_AD_2_1, features = 'mean_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_6M_AD_2_1_mean_upper_p <- SpatialFeaturePlot(ST_6M_AD_2_1, features = 'mean_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)
ST_6M_AD_2_1_mean_lower_p | ST_6M_AD_2_1_mean_upper_p

max_limit <- max(na.omit(ST_6M_AD_2_1@meta.data$stddev_lower),na.omit(ST_6M_AD_2_1@meta.data$stddev_upper))
ST_6M_AD_2_1_stddev_lower_p <- SpatialFeaturePlot(ST_6M_AD_2_1, features = 'stddev_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_6M_AD_2_1_stddev_upper_p <- SpatialFeaturePlot(ST_6M_AD_2_1, features = 'stddev_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)
ST_6M_AD_2_1_stddev_lower_p | ST_6M_AD_2_1_stddev_upper_p

# 15 months
ST_15M_AD_2_1 <- Load10X_Spatial(data.dir = "./ST_15M_AD_2_1/",filename = "filtered_feature_bc_matrix.h5")
ST_15M_AD_2_1_metadata <- data.frame(barcode=rownames(ST_15M_AD_2_1@meta.data))
ST_15M_AD_2_1_abeta <- abeta_meta_list[['ST_15M_AD_2_1']]
for (i in 1:dim(ST_15M_AD_2_1_abeta)[1]){
    ST_15M_AD_2_1_abeta[i,1] <- strsplit(ST_15M_AD_2_1_abeta[i,1],'_')[[1]][1]
}
ST_15M_AD_2_1_metadata_abeta <- left_join(ST_15M_AD_2_1_metadata,ST_15M_AD_2_1_abeta,by = "barcode")
rownames(ST_15M_AD_2_1_metadata_abeta)<-ST_15M_AD_2_1_metadata_abeta$barcode
ST_15M_AD_2_1_metadata_abeta <- ST_15M_AD_2_1_metadata_abeta[,-1]

ST_15M_AD_2_1 <- AddMetaData(
  object = ST_15M_AD_2_1,
  metadata = ST_15M_AD_2_1_metadata_abeta[,c('mean_lower','stddev_lower','mean_upper','stddev_upper')]
)

max_limit <- max(na.omit(ST_15M_AD_2_1@meta.data$mean_lower),na.omit(ST_15M_AD_2_1@meta.data$mean_upper))
ST_15M_AD_2_1_mean_lower_p <- SpatialFeaturePlot(ST_15M_AD_2_1, features = 'mean_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_15M_AD_2_1_mean_upper_p <- SpatialFeaturePlot(ST_15M_AD_2_1, features = 'mean_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)
ST_15M_AD_2_1_mean_lower_p | ST_15M_AD_2_1_mean_upper_p

max_limit <- max(na.omit(ST_15M_AD_2_1@meta.data$stddev_lower),na.omit(ST_15M_AD_2_1@meta.data$stddev_upper))
ST_15M_AD_2_1_stddev_lower_p <- SpatialFeaturePlot(ST_15M_AD_2_1, features = 'stddev_lower') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))
ST_15M_AD_2_1_stddev_upper_p <- SpatialFeaturePlot(ST_15M_AD_2_1, features = 'stddev_upper') +  
scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100),limits = c(0,max_limit))

options(repr.plot.width=14, repr.plot.height=7)

ST_15M_AD_2_1_stddev_lower_p | ST_15M_AD_2_1_stddev_upper_p
