pkg <- c('ppcor','tibble','randomForest','circlize','ComplexHeatmap','pROC','tidyr','ggforce','Boruta','nlme',
         'gridExtra','aod',
         'dplyr','ggplot2','phyloseq','RColorBrewer','ggpubr','reshape','scales','made4','metagMisc','readr')
sapply(pkg, require, character = TRUE)
setwd('/Users/m216453/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Mayo_project/Grover/FinalSubmit/')

## Loading data
load('data/Data.wk.RData')
source('code/Stats.R')
source('code/DailyCode.R')

data.obj$meta.dat$PA_c <- gsub('HV','Healthy volunteers',data.obj$meta.dat$PA_c)
data.obj$meta.dat$PA_c <- gsub('High','High PA PI-IBS',data.obj$meta.dat$PA_c)
data.obj$meta.dat$PA_c <- gsub('Low','Low PA PI-IBS',data.obj$meta.dat$PA_c)
human_meta <- data.obj$meta.dat

#===========Alpha Diversity of Human samples ================
otu <- data.obj$abund.list[['Species']] 
rownames(otu) <- gsub('.*s__','',rownames(otu))
otutable <- otu_table(otu, taxa_are_rows = TRUE) 
taxa <- data.frame(taxa = (data.obj$abund.list[['Species']]%>% rownames(.))) %>% 
  separate(taxa, into <- c('Phylum','Species'),sep <- ';') %>% 
  mutate(taxon = paste0('taxon',1:nrow(.)), name = gsub('s__','', Species)) %>% 
  column_to_rownames('name')
phy <- merge_phyloseq(otu_table(otu, taxa_are_rows = TRUE),tax_table(as.matrix(taxa)), sample_data(data.obj$meta.dat))

Observed <- estimate_richness(phyloseq_standardize_otu_abundance(phy, method = 'pa'), measures=c("Observed")) %>% rownames_to_column('SampleID')
alpha.obj <- generate_alpha_diversity(data.obj,rarefy=FALSE, depth=NULL, iter.no=5, measures=c('Shannon', 'InvSimpson'), seed=123)
alpha <- alpha.obj %>% rownames_to_column('SampleID') %>% inner_join(Observed) %>%
  inner_join(human_meta %>% rownames_to_column('SampleID') %>% dplyr::select(c('SampleID','PA_c','BMI_c'))) %>% 
  dplyr::filter(PA_c != 'PInoIBS')  %>% droplevels() %>% column_to_rownames('SampleID')  
alpha$PA_c <- as.factor(alpha$PA_c)

## linear regression on Alpha diversity
lm.obj <- lm(alpha$Shannon ~ alpha$BMI_c + alpha$PA_c)
obs0 <- prmatrix(summary(lm.obj)$coefficients)
lm.obj <- lm(alpha1$InvSimpson ~ alpha1$BMI_c + alpha1$PA_c)
obs1 <- prmatrix(summary(lm.obj)$coefficients)
lm.obj <- lm(alpha1$Observed ~ alpha1$BMI_c + alpha1$PA_c)
obs2 <- prmatrix(summary(lm.obj)$coefficients)

# Violin plot of alpha diversity
ggviolin(alpha, x = 'PA_c', y = 'Shannon', fill = 'PA_c',bxp.errorbar = TRUE,outlier.shape = NA,
         palette = c(brewer.pal(8,'Dark2')[5],'#f79646',brewer.pal(8,'Paired')[2]),
         add = c("none"), add.params = list(fill = "white")) + 
  geom_boxplot(data = alpha, aes(x = PA_c,y = Shannon), width =0.2, color = 'black', outlier.shape = NA) +
  geom_point(data = alpha, aes(x = PA_c,y = Shannon,fill = PA_c), 
             position=position_jitterdodge(jitter.width =0.2), size =1, color = 'black') +
  theme_bw() +
  theme(axis.text = element_text(color="black", size = 18),
        axis.title = element_text(color="black", size = 18),
        axis.text.x = element_text(color="black", size = 18, angle = 20, hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(2.5,4.3)) + 
  labs(x='',y='Shannon') +
  guides(fill=FALSE) +
  geom_bracket(xmin = c("Healthy volunteers", "High PA PI-IBS"), xmax = c("High PA PI-IBS", "Low PA PI-IBS"),
               y.position = c(4.1, 4.2),size = 1,label.size = 7,
               label = c(paste0('*'), paste0('**')))  
ggsave(path = 'result/','Fig1c_Shannon_jitter.pdf', width = 5, height =5, dpi = 100)

ggviolin(alpha, x = 'PA_c', y = 'InvSimpson', fill = 'PA_c',bxp.errorbar = TRUE,
         palette = c(brewer.pal(8,'Dark2')[5],'#f79646',brewer.pal(8,'Paired')[2]),
         add = "none", add.params = list(fill = "white"))+ 
  geom_boxplot(data = alpha, aes(x = PA_c,y = InvSimpson), width =0.2, color = 'black', outlier.shape = NA) +
  geom_point(data = alpha, aes(x = PA_c,y = InvSimpson,fill = PA_c), 
             position=position_jitterdodge(jitter.width =0.2), size =1, color = 'black') +
  theme_bw()+
  theme(axis.text = element_text(color="black", size = 18),
        axis.title = element_text(color="black", size = 18),
        axis.text.x = element_text(color="black", size = 18, angle = 20, hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())+
  labs(x='',y='InvSimpson') +
  guides(fill=FALSE) +
  scale_y_continuous(limits = c(5,30)) + 
  geom_bracket(xmin = 2, xmax = 3,size = 1,label.size = 7,
               y.position = c(29), label = paste0('*'))
ggsave(path = 'result/','Fig1c_InvSimpson_jitter.pdf', width = 5, height =5, dpi = 100)

ggviolin(alpha, x = 'PA_c', y = 'Observed', fill = 'PA_c',bxp.errorbar = TRUE,
         palette = c(brewer.pal(8,'Dark2')[5],'#f79646',brewer.pal(8,'Paired')[2]),
         add = "boxplot", add.params = list(fill = "white"))+ 
  geom_boxplot(data = alpha1, aes(x = PA_c,y = Observed), width =0.2, color = 'black', outlier.shape = NA) +
  geom_point(data = alpha1, aes(x = PA_c,y = Observed,fill = PA_c), 
             position=position_jitterdodge(jitter.width =0.2), size =1, color = 'black') +
  
  theme_bw()+
  theme(axis.text = element_text(color="black", size = 18),
        axis.title = element_text(color="black", size = 18),
        axis.text.x = element_text(color="black", size = 18, angle = 20, hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())+
  labs(x='',y='Observed') +
  scale_y_continuous(limits = c(0,92)) + 
  guides(fill=FALSE) +
  geom_bracket(xmin = c("Healthy volunteers", "High PA PI-IBS"), xmax = c("High PA PI-IBS", "Low PA PI-IBS"),
               y.position = c(86, 88),size = 1,label.size =7,
               label = c('**', '***'))
ggsave(path = 'result/','Fig1c_Observed_jitter.pdf', width = 5, height =5, dpi = 100)


#===========PCoA of human samples ===============
dist.temp <- dist.obj$BC_gen
obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
PCs <- cbind(obj$points[, 1], obj$points[, 2]) %>% as.data.frame()
PCs <- merge(PCs,human_meta, by = 0)%>% dplyr::rename(SampleID = Row.names)
head(PCs)
xlab <- paste0('PC1,', pve[1], '%')
ylab <- paste0('PC2,', pve[2], '%')
col <- c('Healthy volunteers' <- brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'<-'#f79646','Low PA PI-IBS'<-brewer.pal(8,'Paired')[2])
PCs <- PCs %>% filter(SampleID != 'S50') # S50 is the outlier, for visualization purpose, we exclude it

ggplot(PCs) +
  geom_point(aes(x = V1, y = V2,fill = PA_c, label = SampleID), size = 4, shape = 21) +
  geom_mark_ellipse(aes(x = V1, y = V2, color= PA_c),expand = unit(0, "mm"))+
  theme_bw() +
  scale_x_continuous(limits = c(-0.6,0.6))+
  scale_y_continuous(limits = c(-0.4,0.4))+
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  labs(fill = '') +
  xlab(xlab) + ylab(ylab) + 
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.title =  element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())+
  guides(color = FALSE)
ggsave(path = 'result/','Fig1b.pdf', width = 8, height =5, dpi = 100)

#===========Taxonomy plot of human samples =========
level = 'Phylum'
for (level in c('Phylum','Class','Family')){
  data.level <- data.obj$abund.list[[level]] 
  rownames(data.level) <- gsub('.*c__|.*p__|.*f__','', rownames(data.level))
  data.level <- apply(data.level, 2, function(x) x/sum(x))
  data.level <-  as.data.frame(t(data.level)) %>% rownames_to_column('SampleID')
  meta.dat <- human_meta %>% rownames_to_column('SampleID')%>% dplyr::filter(PA_c %in% c("High PA PI-IBS","Healthy volunteers")) %>% dplyr::select(c('SampleID','PA_c'))
  data.level <- inner_join(data.level, meta.dat) %>% column_to_rownames('SampleID') %>% group_by(PA_c) %>% 
    dplyr::summarise_each(funs(mean)) %>% column_to_rownames('PA_c') %>% t() %>% as.data.frame() %>% rownames_to_column(level)
  data.level1 <- data.level %>% dplyr::filter(`Healthy volunteers` > 0.01| `High PA PI-IBS` > 0.01) 
  data.level2 <- data.level1 %>% add_row(tibble_row(`Healthy volunteers` = 1-sum(data.level1$`Healthy volunteers`), 
                                                    `High PA PI-IBS` = 1-sum(data.level1$`High PA PI-IBS`)))
  data.level2[,level][nrow(data.level2)] <- 'Others'
  data.level <- data.level2 %>% melt()
  sort <- data.level %>% filter(data.level$variable == 'Healthy volunteers')
  sort <- sort[order(sort$value),][[level]]
  data.level[[level]] <- factor(data.level[[level]], levels = sort)

  p1 <- ggplot(data.level, aes(x= variable, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black", size = 0.2) +
    scale_fill_brewer(palette = 'Paired') +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() + 
    theme(text = element_text(size = 20, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 20, angle = 30, vjust = 0.9, hjust = 0.8),
          axis.text.y = element_text(color="black", size = 20),
          legend.text=element_text(size=20, face = 'italic'),
          legend.title=element_text(size=20),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'))
  legend <- cowplot::get_legend(p1)
  ggarrange(legend)
  ggsave(filename = paste0('result/','Fig1d_',level,'legend.pdf'),width = 3.5, height = 4.5, dpi = 500)
  graphics.off()

  p2 <- ggplot(class, aes(x= variable, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black", size = 0.2) +
    scale_fill_brewer(palette = 'Paired') +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() + 
    # scale_x_discrete(expand = c(0.5, 1)) +
    theme(text = element_text(size = 20, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 20, angle = 30, vjust = 0.9, hjust = 0.8),
          axis.text.y = element_text(color="black", size = 20),
          legend.text=element_text(size=20, face = 'italic'),
          legend.title=element_text(size=20),
          strip.text.x = element_text(angle=0, size = 20),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'),
          legend.position = 'none')
  ggsave(paste0('result/','Fig1d_',level,'main.pdf'),width = 3.5, height = 4.5, dpi = 500)
}


#===========Heatmap of human samples ========================
#  Fig1e #q < 0.2 on heatmap, correlation r = PA(all_samples), correlation q 
## Spearman correlation(select taxa based on the result of spearman q < 0.2, all these species will be plotted on the figure)
spr <- read.csv('data/Fig1f.csv') %>%  dplyr::rename(Species= X) %>% 
  dplyr::select(c('Species','Qvalue','SpearmanCorr')) %>% dplyr::filter(Qvalue <= 0.2) %>% 
  mutate(Species=gsub('.*;s__','',Species),grp = ifelse(Qvalue <= 0.1,'#a31415', 'black')) %>% 
  na.omit()

## Differential abundance result
fig1e <- read.csv('data/Fig1e.csv') %>% dplyr::rename(Species= X) %>% dplyr::select(c("Species","Qvalue","coef_")) %>% 
  mutate(Species = gsub('.*;s__','',Species)) %>% dplyr::filter(Species %in% c(spr$Species));dim(fig1e)
diff <- c(setdiff(spr$Species,fig1e$Species), setdiff(fig1e$Species,spr$Species) )
add <- c()
for(i in 1:length(diff)){
  add <- rbind(add, c(diff[i], 1, NA))
}
add <- as.data.frame(add)
colnames(add) <- colnames(fig1e)
fig1e <- rbind(fig1e, add) %>% column_to_rownames('Species')
fig1e <- fig1e[spr$Species,];dim(fig1e)

submet <- human_meta %>% rownames_to_column('SampleID')%>% filter(PA_c %in% c('High PA PI-IBS','Healthy volunteers')) %>% dplyr::select(c('SampleID','PA_c'))
Species <- data.obj$abund.list$Species 
rownames(Species) <- gsub('.*;s__','',rownames(Species))
Species <- apply(Species, 2, function(x) x/sum(x)) %>% as.data.frame() %>% 
  rownames_to_column('Species') %>% filter(Species %in% spr$Species) %>% 
  column_to_rownames('Species') %>% t() %>% as.data.frame()%>% rownames_to_column('SampleID') %>% 
  inner_join(submet) 
Species <- Species[order(Species$PA_c),];dim(Species);Species[1:5,1:5]
type1 <- as.vector(Species$PA_c)

mat <- Species %>% dplyr::select(-PA_c) %>% as.data.frame();mat[1:5,1:5]
rownames(mat) <- NULL
mat <- as.data.frame(mat) %>% column_to_rownames('SampleID')
mat <- mat[,spr$Species] %>% t() %>% sqrt()

ha = HeatmapAnnotation(HA = type1, annotation_height = 1,border = TRUE, 
                       col = list(HA = c('Healthy volunteers' = brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'='#f79646')),
                       # annotation_name_side = 'right',
                       simple_anno_size = unit(0.3, "cm"), annotation_label ='',
                       show_legend= FALSE)
pvalue_col_fun = colorRamp2(c(-0.6, 0, 0.6), c("#377EB8", "white", "#a31415"))#colorRamp2(c(-0.8, 0, 0.8), c(brewer.pal(11,'RdBu')[c(11,6,1)]))
heatmap_legend_param = list(title = 'Expression Level', 
                            at = seq(0, 1, 0.1), 
                            direction = 'horizontal', 
                            labels = c(0, '', '', '', '',0.5,'','','','', 1), 
                            color_bar = 'continuous', title_gp = gpar(fontsize = 15), 
                            legend_gp = gpar(fontsize = 20, fontfamily="Arial"), legend_height = unit(2, 'in'), 
                            legend_width = unit(2, 'in'), title_position="topcenter")

list <- Heatmap(spr$SpearmanCorr, name = "correlation", 
                col = pvalue_col_fun, show_row_dend = FALSE,border = TRUE,
                heatmap_legend_param = list(direction = "horizontal")) +
  Heatmap(mat, name = "sqrt(proportion)", heatmap_legend_param = list(direction = "horizontal"),
          col = colorRamp2(c(0, 0.2, max(mat)), brewer.pal(9, 'Purples')[c(1,5,9)]),
          top_annotation = ha, column_gap = unit(1, "mm"), 
          border = TRUE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(fig1e$Qvalue[i] <= 0.1) {
              grid.circle(x = x, y = y, r = unit(1/20,'cm'), gp = gpar(fill = 'black', col = 'black'))
            }
          },
          column_split = factor(c(rep("Healthy volunteers",21), rep("High PA PI-IBS", 12))),
          rect_gp = gpar(col= "white"),
          column_order= c(1:33), show_column_names = FALSE,
          show_row_dend = FALSE,
          row_names_gp = gpar(col = spr$grp, fontface = 'italic'),
          row_names_max_width = max_text_width(rownames(spr), gp = gpar(fontsize = 18))
  ) 

lgd_list = list(
  Legend(labels = c(''), title = "q < 0.1", type = "points", pch = 20, 
         legend_gp = gpar(col = 'black'))
)

pdf('result/Fig1f.pdf', width = 8, height = 8)
draw(list, merge_legend = TRUE, padding = unit(c(1, 1, 1,8), "cm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     annotation_legend_list = lgd_list, 
     use_raster = TRUE)
dev.off()




fig1e <- read.csv('data/Fig1e.csv') %>% filter(Qvalue <= 0.1) %>% dplyr::rename(Species= X) 
submet <- human_meta %>% rownames_to_column('SampleID')%>% filter(PA_c %in% c('High PA PI-IBS','Healthy volunteers')) %>% dplyr::select(c('SampleID','PA_c'))
Species <- data.obj$abund.list$Species 
Species.prop <- t(t(Species) /colSums(Species))
sub <- Species.prop[fig1e$Species,] %>% t() %>% merge(submet %>% column_to_rownames('SampleID'),by =0 ) %>% column_to_rownames('Row.names')
subm <- melt(sub)
subm$variable <- gsub('.*;s__','',subm$variable)
col <- c('Healthy volunteers' = brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'='#f79646')
ggplot(subm, aes(x = reorder(variable, -value), y = value,fill = PA_c)) + 
  stat_boxplot(geom ='errorbar', lwd = 1) +
  geom_boxplot(outlier.shape = NA, lwd = 1) + 
  scale_color_manual(values = col) +
  theme_bw() +
  scale_fill_manual(values = c('white','white')) + #
  geom_point(position=position_jitterdodge(jitter.width =0.1), size =1,aes(color = PA_c)) + 
  scale_y_continuous(trans = sqrt_trans(),
                     # limits = c(0,0.05),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2)))+
  labs(x = '', y = 'Proportion', fill = '')+
  theme(axis.text.x = element_text(color="black", size = 20, angle = 90,vjust = 0.35, hjust = 1, face = 'italic'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        text = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank())

Species <- apply(data.obj$abund.list$Species, 2, function(x) x/sum(x)) %>% as.data.frame() %>% 
  rownames_to_column('Species') %>% dplyr::filter(Species %in% fig1e$Species) %>% 
  column_to_rownames('Species') %>% as.matrix() %>% t() %>% as.data.frame()%>% rownames_to_column('SampleID') %>%
  inner_join(submet) %>% column_to_rownames('SampleID') %>% 
  melt() %>% mutate(variable=gsub('.*\\;','',variable), value = sqrt(value))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                        varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}
Species.m <- data_summary(Species, varname="value", groupnames=c("variable", "PA_c"))
Species.m$variable <- gsub('s__','',Species.m$variable)
Species.m$PA_c <- factor(Species.m$PA_c, levels = c('Healthy volunteers', 'High PA PI-IBS'))
col <- c('Healthy volunteers' = brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'='#f79646')

ggplot(Species.m, aes(x = reorder(variable, -value), y = value,  fill = PA_c)) +
  geom_bar(stat="identity", color="black", lwd =1,position=position_dodge()) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.5,lwd =1,position=position_dodge(.9)) +
  geom_point(data = subm, position=position_jitterdodge(jitter.width =0.2), size =0.8, color = 'black') + 
  scale_fill_manual(values = col) +
  theme_bw() +
  labs(x = '', y = 'Proportion', fill = '')+
  theme(axis.text.x = element_text(color="black", size = 20, angle = 90,vjust = 0.35, hjust = 1, face = 'italic'),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        text = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank()) + 
  guides(fill = guide_legend(show = FALSE))
ggsave(path = 'result/',paste0('Fig1g_barplot',level, '.pdf'), width = 10, height = 8, dpi = 300)



fig1e <- read.csv('data/Fig1f.csv') %>% dplyr::rename(Species= X) #%>% filter(Qvalue.perm < 0.2)
submet <- human_meta %>% rownames_to_column('SampleID')%>% filter(PA_c %in% c('High PA PI-IBS','Healthy volunteers')) %>% dplyr::select(c('SampleID','PA_c','PA'))
Species <- data.obj$abund.list$Species 
Species <- apply(Species, 2, function(x) x/sum(x)) %>% as.data.frame() %>% rownames_to_column('Species') %>% filter(Species %in% fig1e$Species)
Species$Species <- gsub('.*s__','',Species$Species)
Species <- Species[,c('Species',submet$SampleID)] %>% 
  column_to_rownames('Species') %>% as.matrix() %>% t() %>% 
  as.data.frame() %>% rownames_to_column('SampleID') %>% inner_join(submet) %>% 
  column_to_rownames('SampleID') %>% mutate(PA =- log(PA)) %>% droplevels()
Full1 <- c("Alistipes_putredinis",'Ruminococcus_bromii','Subdoligranulum_unclassified')
fig1e$Species <- gsub('.*s__','',fig1e$Species)
Full0 <- fig1e %>% filter(Species %in% Full1);head(Full0)
Full0$Species <- gsub('_unclassified','',Full0$Species)
scatter <- Species %>% dplyr::select(c('PA','PA_c',Full1)) %>% melt(id <-c('PA','PA_c'))
colnames(scatter)[2] <- 'PA Status'
scatter$variable <- gsub('_unclassified','',scatter$variable)
Full1 <- gsub('_unclassified','',Full1)

col <- c('Healthy volunteers'<-brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'<-'#f79646','Low PA PI-IBS'<-brewer.pal(8,'Paired')[2])
variables <- c()
for(i in 1:length(Full1)){
  q = format(round(Full0$Qvalue[i],3), scientific = F);r = format(round(Full0$SpearmanCorr[i],2), scientific = F)
  variables = c(variables,paste0('r=',r,',q=',q))
}

variables_name = c()
for(i in 1:length(Full1)){
  variables_name = c(variables_name,Full1[i])
}

ann_text <- data.frame(PA = c(8,8,8),value = c(0.15,0.12,0.15),lab = variables,
                       variable = variables_name) # BH 
graphics.off()
ggplot(scatter, aes(x = PA,y = value)) +
  facet_wrap(~variable, scales = 'free',nrow= 3, strip.position = "left") +
  geom_point(aes(color =`PA Status`,fill =`PA Status`), cex =5, shape = 17) +
  geom_smooth(method=lm, color = 'grey50')+
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = NULL, x = 'log(PA)', color = '') +
  guides(fill=guide_legend("")) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(text = element_text(size = 24, color = 'black'),
        axis.text = element_text(size = 24, color = 'black'),
        axis.title = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        strip.text = element_text(size = 24, color = 'black', face = 'italic'),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        strip.placement = "outside") +
  geom_text(
    data    = ann_text, size = 9,
    mapping = aes(x = -Inf, y = -Inf, label = lab),
    hjust   = -0.2,
    vjust   = -10
  )

ggsave(path = 'result/', 'Fig1f.pdf', width = 9, height =12, dpi = 100)


#===========Random forest prediction ===================
GPr  <- transform_sample_counts(phy, function(x) x / sum(x) )
GPfr <- filter_taxa(GPr, function(x) sum(x > 0.002) > (0.1*length(x)), TRUE) # same filter as random forest
otu <- t(otu_table(GPfr)@.Data) #extract the otu table from phyloseq object
tax <- as.data.frame(tax_table(GPfr)@.Data)#extract the taxonomy information
tax$label = gsub('s__','',tax$Species)

# ---- Boruta feature selection
b = as.data.frame(otu) %>% rownames_to_column('SampleID') %>% inner_join(human_meta %>% rownames_to_column('SampleID')%>%dplyr::select(c('SampleID',"PA_c2"))) %>% column_to_rownames('SampleID') %>% droplevels()
set.seed(222222)
boruta.wk6_train <- Boruta(PA_c2 ~., data = b, doTrace = 2)
boruta.wk6 <- TentativeRoughFix(boruta.wk6_train)
print(boruta.wk6)
par(mar = c(15, 4, 2, 2))

# traditional plot
bt <- getSelectedAttributes(boruta.wk6, withTentative = F)
plot(boruta.wk6, xlab = "", xaxt = "n", ylab='Importance z-score')
lz <-lapply(1:ncol(boruta.wk6$ImpHistory),function(i)
  boruta.wk6$ImpHistory[is.finite(boruta.wk6$ImpHistory[,i]),i])
name <- colnames(boruta.wk6$ImpHistory)
name <- gsub('`','',name)
names(lz) <- name
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.wk6$ImpHistory), cex.axis = 1)

# customized plot
dd <- as.data.frame(boruta.wk6$ImpHistory)
dd1 <- dd %>% dplyr::select(-c('shadowMean', 'shadowMin','shadowMax'))
mns <-apply(dd, 2, function(x) median(x))
apply(dd, 2, function(x) median(x))
dd <- dd[,order(mns)]
or <- as.vector(colnames(dd))
df <- as.data.frame(boruta.wk6$ImpHistory) %>% melt()
df$variable <- factor(df$variable,levels = names(Labels),ordered = TRUE)
select <- getSelectedAttributes(boruta.wk6, withTentative = F)
df$col <- 'red'
shadow <- unique(df$variable)[grep('^shadow', unique(df$variable))]
df[(df$variable %in% select)& !(df$variable %in% shadow),'col'] ='green'
df[(df$variable %in% shadow),'col'] = 'blue'
cols = c('blue'='grey30','red'='grey80','green'=brewer.pal(9,'Set1')[1])
ggplot(df, aes(x = variable, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(aes(fill = col),outlier.size = 0.8, outlier.colour = 'grey60') +theme_bw() +
  scale_fill_manual(values = cols)+
  labs(y='Importance z-score', x = '')+
  theme(text = element_text(size = 20, color = "black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length=unit(0.1, "cm"),
        axis.text.x = element_text(color="black", size = 18, angle = 90, hjust = 1, vjust = 0.25, face = 'italic'),
        axis.text.y = element_text(color="black", size = 20),
        strip.text.x = element_text(angle=0, size = 20),
        legend.position = "none")
ggsave(path = 'result/', 'Fig1h.pdf', width =17, height =10, dpi = 100)



# random forest prediction
group = 'PA_c2'
Sp = data.obj$abund.list[[level]]
Sp = apply(Sp, 2, function(x) x/sum(x))
Sp = as.data.frame(Sp)
Sp$taxon = paste0(level,seq(1:nrow(Sp)))
Sp_name = Sp %>% dplyr::select(taxon) %>% rownames_to_column(level)
Sp_name$Species = gsub('.*s__','',Sp_name$Species)
bt1 = Sp_name[Sp_name$Species %in% bt,]$taxon
rownames(Sp) = Sp$taxon 
Sp = Sp %>% dplyr::select(-taxon)
human_meta = data.obj$meta.dat;dim(human_meta)

#--- use	PIIBS vs HV
PIIBS.meta = human_meta %>% rownames_to_column('SampleID') %>% 
  dplyr::select(c('SampleID',group,'Phenotype')) %>% filter(Phenotype != 'PInoIBS') %>% dplyr::select(-Phenotype) ;dim(PIIBS.meta)
PIIBS.meta$PA_c2
PIIBS.Sp = t(Sp[,colnames(Sp) %in% (PIIBS.meta$SampleID)]) %>% as.data.frame();dim(PIIBS.Sp)
idx <- apply(PIIBS.Sp, 2,function(x){sum(x > 0.002) > (nrow(PIIBS.Sp)*0.1)})
PIIBS.Sp = as.data.frame(PIIBS.Sp[,idx]) %>% rownames_to_column('SampleID');dim(PIIBS.Sp)
rowSums(PIIBS.Sp[,-1]);dim(PIIBS.Sp)

RF_Sp = PIIBS.meta %>% inner_join(PIIBS.Sp) #%>% dplyr::select(c('SampleID','PA_c2',bt1)) %>% droplevels()


group_names = RF_Sp$SampleID
prob.list <- accuracy.list <- list()
set.seed(12111)
for (i in 1:500){
  g = sample(group_names, 1, replace = T)
  if (g=='S24') {cat(g,'\n')}
  train_data <- RF_Sp[RF_Sp['SampleID']!=g, ][,-1]
  test_data <- RF_Sp[RF_Sp['SampleID']==g, ][,-1]
  formula = as.formula(paste0(group ,"~ ."))
  data.rf <- randomForest(formula, data=train_data, ntree=501, importance=TRUE, proximities=TRUE)
  data.pred <- predict(data.rf, test_data, type = 'prob')
  accuracy.list[[paste(g)]] <- cbind(accuracy.list[[paste(g)]], data.rf$importance[,3])
  prob.list[[paste(g)]] <- cbind(prob.list[[paste(g)]], data.pred[,1])
}


probs <- c()
accuracys <- c()
labels <- c()
for (g in group_names){
  # mean probability
  accuracy <- rowMeans(accuracy.list[[paste(g)]])
  prob <- rowMeans(prob.list[[paste(g)]])
  cat(length(accuracy))
  cat(length(prob))
  cat(' ')
  probs <- c(probs, prob)
  accuracys <- cbind(accuracys, accuracy)
  # label
  label <- as.vector(RF_Sp[RF_Sp['SampleID']==g, ][, group])
  labels <- c(labels, label)
}


labels <- as.numeric(as.factor(labels))
labels.list <- probs.list <- list()
labels.list[[level]] <- labels
probs.list[[level]] <- probs
roc = pROC::roc(labels,probs)
df <- cbind(roc$sensitivities*100,roc$specificities*100, sens.ci[,c(1,3)] *100) %>% as.data.frame()
colnames(df) <- c('sensitivities','specificities','se.low','se.high')
size = 20
ggplot(df, aes(x = specificities,y = sensitivities))+
  geom_ribbon(aes(ymin = se.low, ymax = se.high), fill = "grey") +
  geom_path(aes(y=se.low), color = 'grey50') +
  geom_path(aes(y=se.high), color = 'grey50') +
  geom_path()+
  theme_bw()+
  xlim(c(100,0))+
  geom_abline(slope=1,intercept=100,color='grey')+
  xlab('Specificity (%)') +
  ylab('Sensitivity (%)')+
  annotate(geom='text', x=45, y=5, label=paste0("AUC:", round(ci(roc)[2]*100,1),'% (',round(ci(roc)[1]*100,1),'%-',round(ci(roc)[3]*100,1),'%)'), 
           size =size/3)+
  theme(axis.text = element_text(color="black", size = size),
        axis.title = element_text(color="black", size = size),
        legend.position = 'none',      
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color="black", hjust = 0.5),
        title =element_text(size=size)) 
ggsave(path = 'result/', 'Fig1i.pdf', width =5, height =5, dpi = 100)

##==========Supplemental Figure 1- Major taxa differences observed between low and high PA PI-IBS patients and taxa correlations with stool PA (Human Data from Fig 1)=======
# 1.Higher level taxonomy for high and low PA PI-IBS patients: Phylum, class, family (boxplots) low PA (n=14) 
for (level in c('Phylum','Class','Family')){
  class = data.obj$abund.list[[level]] 
  rownames(class) = gsub('.*c__|.*p__|.*f__','', rownames(class))
  class = apply(class, 2, function(x) x/sum(x))
  class =  class %>% t()%>%as.data.frame() %>% rownames_to_column('SampleID')
  meta.dat = human_meta %>% rownames_to_column('SampleID')%>% 
    dplyr::filter(PA_c %in% c('High PA PI-IBS','Low PA PI-IBS')) %>% dplyr::select(c('SampleID','PA_c'))
  class = inner_join(class, meta.dat) %>% column_to_rownames('SampleID') %>% group_by(PA_c) %>% 
    dplyr::summarise_each(funs(mean)) %>% column_to_rownames('PA_c') %>% t() %>% as.data.frame() %>% rownames_to_column(level)
  class1 = class %>% dplyr::filter(`Low PA PI-IBS` > 0.01| `High PA PI-IBS` > 0.01) 
  class2 = class1 %>% add_row(tibble_row(`Low PA PI-IBS` = 1-sum(class1$`Low PA PI-IBS`), `High PA PI-IBS` = 1-sum(class1$`High PA PI-IBS`)))
  class2[,level][nrow(class2)] = 'Others'
  class = class2 %>% melt()
  sort <- (class[class$variable == 'Low PA PI-IBS',] %>% arrange(desc(value)))[,level]
  class[[level]] = factor(class[[level]], levels = sort)
  cols <- c(RColorBrewer::brewer.pal(12, 'Paired'),RColorBrewer::brewer.pal(8, 'Set2'))
  p1 = ggplot(class, aes(x= variable, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black", size = 0.2) +
    scale_fill_manual(values = cols) +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() + 
    theme(text = element_text(size = 20, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 20, angle = 30, vjust = 0.9, hjust = 0.8),
          axis.text.y = element_text(color="black", size = 20),
          legend.text=element_text(size=20, face = 'italic'),
          legend.title=element_text(size=20),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'))
  legend <- cowplot::get_legend(p1)
  ggarrange(legend)
  ggsave(filename = paste0('result/','FigS1_',level,'legend.pdf'),width = 3.5, height = 4.5, dpi = 500)
  graphics.off()
  p2 = ggplot(class, aes(x= variable, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black", size = 0.2) +
    scale_fill_manual(values = cols) +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() + 
    # scale_x_discrete(expand = c(0.5, 1)) +
    theme(text = element_text(size = 20, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 20, angle = 30, vjust = 0.9, hjust = 0.8),
          axis.text.y = element_text(color="black", size = 20),
          legend.text=element_text(size=20, face = 'italic'),
          legend.title=element_text(size=20),
          strip.text.x = element_text(angle=0, size = 20),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'),
          legend.position = 'none')
  ggsave(paste0('result/','FigS1_',level,'main.pdf'),width = 3.5, height = 4.5, dpi = 500)
}


##==========Supplemental Figure 2-Correlation scatterplots between PA and taxa abundance (Human Data from Fig 1)============
# 1.Bacterial taxa identified as differentially abundant (q<0.1) are plotted against the logPA of patient stool. 
# Correlation coefficients and q-values come from comparisons within the larger cohort to assess the relationship 
# between PA and taxa abundance (n=21 HV, 12 high PA PI-IBS, 14 low PA PI-IBS, colored as in paper)  (JC/LY). 
# (Similar plots as to what was included in the main figure of the manuscript)
# r = spearman correlation; q = spearman q value
fig1e = read.csv('data/Fig1f.csv') %>% dplyr::rename(Species= X) %>% dplyr::filter(Qvalue <= 0.1)
submet = human_meta %>% rownames_to_column('SampleID')%>% filter(PA_c %in% c('High PA PI-IBS','Healthy volunteers')) %>% dplyr::select(c('SampleID','PA_c','PA'))
Species = data.obj$abund.list$Species 
Species = apply(Species, 2, function(x) x/sum(x)) %>% as.data.frame() %>% rownames_to_column('Species') %>% filter(Species %in% fig1e$Species)
Species$Species = gsub('.*s__','',Species$Species)
Species = Species[,c('Species',submet$SampleID)] %>% 
  column_to_rownames('Species') %>% as.matrix() %>% t() %>% 
  as.data.frame() %>% rownames_to_column('SampleID') %>% inner_join(submet) %>% 
  column_to_rownames('SampleID') %>% mutate(PA = log(PA)) %>% droplevels()
# combine part 
fig1e$Species = gsub('.*s__','',fig1e$Species)
Full1 = fig1e$Species[!fig1e$Species %in% c("Alistipes_putredinis",'Ruminococcus_bromii','Subdoligranulum_unclassified')]
Full0 = fig1e %>% dplyr::filter(Species %in% Full1);head(Full0)
scatter = Species %>% dplyr::select(c('PA','PA_c',Full1)) %>% melt(id =c('PA','PA_c'))
colnames(scatter)[2] = 'PA Status'
scatter$variable = gsub('_unclassified','',scatter$variable)
Full1 = gsub('_unclassified','',Full1)
col = c('Healthy volunteers'=brewer.pal(8,'Dark2')[5], 'High PA PI-IBS'='#f79646','Low PA PI-IBS'=brewer.pal(8,'Paired')[2])
variables = variables_name = c()
for(i in 1:length(Full1)){
  q = format(round(Full0$Qvalue[i],3), scientific = F);r = format(round(Full0$SpearmanCorr[i],2), scientific = F)
  variables = c(variables,paste0('r=',r,',q=',q))
  variables_name = c(variables_name,Full1[i])
}

ann_text <- data.frame(PA = c(rep(8,11)),
                       value = c(0.001,0.15,0.025,0.005,0.12,0.01, 0.08,0.002, 0.01, 0.1, 0.15),
                       lab = variables,variable = variables_name) 
ggplot(scatter, aes(x = PA,y = value)) +
  facet_wrap(~variable, scales = 'free',nrow= 4, strip.position = "left") +
  geom_point(aes(color =`PA Status`,fill =`PA Status`), cex =5, shape = 17) +
  geom_smooth(method=lm, color = 'grey50')+
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = NULL, x = 'log(PA)', color = '') +
  guides(fill=guide_legend("")) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(text = element_text(size = 24, color = 'black'),
        axis.text = element_text(size = 24, color = 'black'),
        axis.title = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.position = 'bottom',
        strip.text = element_text(size = 24, color = 'black', face = 'italic'),
        strip.background = element_blank(),
        strip.placement = "outside") +
  geom_text(
    data    = ann_text, size = 9,
    mapping = aes(x = -Inf, y = -Inf, label = lab),
    hjust   = -0.2,
    vjust   = -10)
ggsave(path = 'result/', 'FigS2.pdf', width =20, height =20, dpi = 300)
















##===== load data ======
tmp <- load('data/primary.RData')
meta <- as.data.frame(as.matrix(sample_data(primary)))
primary.t <- primary %>% transform_sample_counts(function(x) x / sum(x))
primary.6wks <- primary %>% subset_samples(obj1 %in% '6 wks')
meta.6wks <- as.data.frame(as.matrix(sample_data(primary.6wks)))

meta$PAStatus_ <- gsub("Healthy\\ volunteers\\ \\(Low\\ PA\\)","Healthy volunteers",meta$PAStatus_)
meta$PAStatus_ <- gsub("PI\\-IBS\\ \\(High\\ PA\\)",'High PA PI-IBS',meta$PAStatus_)
meta$PAStatus_ <- gsub('PI\\-IBS\\ \\(Low\\ PA\\)','Low PA PI-IBS',meta$PAStatus_)

meta.6wks$PAStatus_ <- gsub("Healthy\\ volunteers\\ \\(Low\\ PA\\)","Healthy volunteers",meta.6wks$PAStatus_)
meta.6wks$PAStatus_ <- gsub("PI\\-IBS\\ \\(High\\ PA\\)",'High PA PI-IBS',meta.6wks$PAStatus_)
meta.6wks$PAStatus_ <- gsub('PI\\-IBS\\ \\(Low\\ PA\\)','Low PA PI-IBS',meta.6wks$PAStatus_)


##======= PCoA of fecal and mice samples ============ 
sample_data(primary)$ecllipse <- as.character(meta$PatientID)
sample_data(primary)[sample_data(primary)$obj1 == 'Fecal Slurry', 'ecllipse'] = ''
sample_data(primary) <- sample_data(as.data.frame(as.matrix(primary@sam_data)) %>% mutate(patientid = PatientID))# %>% unite('label',PAStatus_:PatientID))
col <- c(brewer.pal(12, 'Set3'),brewer.pal(12, 'Paired'))
pslog <- transform_sample_counts(primary, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
eng <- as.data.frame(out.pcoa.log$vectors)[,c(1:2)] %>% rownames_to_column('SampleID') %>% 
  inner_join(as.data.frame(as.matrix(sample_data(primary))))
pc1 <- round((out.pcoa.log$values$Eigenvalues[1])/sum(out.pcoa.log$values$Eigenvalues),2)
pc2 <- round((out.pcoa.log$values$Eigenvalues[2])/sum(out.pcoa.log$values$Eigenvalues),2)

p1 <- ggplot(eng) +
  # add text and ellipse for mice part
  ggforce::geom_mark_ellipse(data = subset(eng, ecllipse != ""),
                             aes(x = Axis.1, y = Axis.2, fill=PAStatus_,
                                 group= ecllipse, label = PatientID)) +
  theme_bw() +
  geom_text_repel(aes(x = Axis.1, y = Axis.2, label = patientid), vjust =0.4,hjust =0.5, size = 4, fontface = "bold") + 
  geom_point(aes(x = Axis.1, y = Axis.2,color = factor(PAStatus_), shape = factor(obj1), fill = factor(PAStatus_)), size =3)+
  scale_color_manual(values = c('black','black', 'black')) +
  scale_shape_manual(values = c(24, 21))+
  scale_fill_manual(values = c(brewer.pal(8,'Dark2')[5], '#f79646',brewer.pal(8,'Paired')[2]))+
  labs(fill = 'Mice', shape = 'Donor') +
  xlab(paste0('PC1, ',(paste0(pc1*100,'%')))) + ylab(paste0('PC2, ',paste0(pc2*100, '%'))) + 
  theme(text = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        axis.text = element_text(size = 24, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.position = 'none')

data = eng %>% filter(obj1 =='Fecal Slurry')
l1 = ggplot(data, aes(x=Axis.1, y = Axis.2))+ 
  geom_point(aes(fill=factor(PAStatus_)), colour="black",pch=21, size=4) +
  theme_bw() +
  scale_fill_manual(values = c(brewer.pal(8,'Dark2')[5], '#f79646',brewer.pal(8,'Paired')[2])) +
  theme(legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(fill = 'Donors')

l1= get_legend(l1)
grid.draw(l1)

data = eng %>% filter(obj1 !='Fecal Slurry')
l2 = ggplot(data, aes(x=Axis.1, y = Axis.2))+ 
  geom_point(aes(fill=PAStatus_), pch=24, size=4) +
  theme_bw() +
  scale_fill_manual(values = c(brewer.pal(8,'Dark2')[5], '#f79646',brewer.pal(8,'Paired')[2])) +
  theme(legend.text=element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=24))+
  labs(fill = 'Mice')
l2= get_legend(l2)
grid.draw(l2)

ll = ggarrange(NULL,l1, l2,NULL, nrow = 4, ncol = 1)
pp =  ggarrange(p1,ll,  widths = c(2, 0.7),nrow = 1, ncol = 2)

ggsave(path = 'result/','Fig4a.pdf', width = 15, height = 10, dpi = 300)



##======= subject level distance permutation: fecal slurry VS  6 wks ============
dist.phy <- primary %>% phyloseq_standardize_otu_abundance(method = "pa")
dist.mat <- as.matrix(as.dist(phyloseq::distance(dist.phy, method = 'bray')))
grp.name <- 'obj1'
subject <- 'PatientID'
df <- as.data.frame(as.matrix(dist.phy@sam_data))
identical(df$SampleID, rownames(dist.mat))
df$obj1 <- as.factor(df$obj1)
df$PatientID <- as.factor(df$PatientID)
grp <- df[, grp.name]
IDs <- df[, subject]
sbj.levels <- levels(df$PatientID)
sbj.nlevels <- nlevels(df$PatientID)
grp.levels <- levels(df$obj1)
grp.nlevels <- nlevels(df$obj1)
dist.mat <- as.dist(phyloseq::distance(dist.phy, method = 'bray'))
dist.mat <- as.matrix(dist.mat)
df <- as.data.frame(as.matrix(dist.phy@sam_data))
# patientID and unique patientID
sbj.levels <- unique(df$PatientID)
patientID <- as.character(df$PatientID)
# calculate the orginial dist12, dist13
dist12.list <- dist13.list <- NULL
for (j in 1 : length(sbj.levels)) {
  fecal <- which(df$obj1 == grp.levels[1] & df$PatientID == sbj.levels[j])
  mice <- which(df$obj1 == grp.levels[2] & df$PatientID == sbj.levels[j]) 
  mice_not <- which(df$obj1 == grp.levels[2] & df$PatientID != sbj.levels[j]) 
  
  ind1 = fecal
  ind2 = mice
  ind3 = mice_not
  
  dist12 <- dist.mat[ind1, ind2]
  dist13 <- dist.mat[ind1, ind3]
  
  dist12.list <- c(dist12.list, mean(dist12))
  dist13.list <- c(dist13.list, mean(dist13))
}

stat.obj <- mean(dist12.list) - mean(dist13.list)

# permutation 
stat.perm <- NULL # the mean distance different for each iteration
for (i in 1:999){# Shuffle by patientID, namely, sbj.level for fecal 
  patientID.p <- patientID
  patientID.fecal <- patientID[df$obj1 == grp.levels[2]] 
  temp <- factor(patientID.fecal)
  levels(temp) <- sample(levels(temp)) ### shuffle the levels of patientID associated with fecal
  patientID.fecal.p <- as.character(temp)
  patientID.p[df$obj1 == grp.levels[2]]  <- patientID.fecal.p # change the patientID with shuffle patientID
  
  # Shuffle by patientID, namely, sbj.level for mice
  patientID.mice <- patientID[df$obj1 == grp.levels[1]] #grp.level[1] = 6wks
  temp <- factor(patientID.mice)
  levels(temp) <- sample(levels(temp)) ### shuffle the levels of patientID associated with mice
  patientID.mice.p <- as.character(temp)
  patientID.p[df$obj1 == grp.levels[1]] <- patientID.mice.p # change the patientID with shuffle patientID
  
  ## calculate the dist12, dist13
  dist12.list <- dist13.list <- NULL
  for (j in 1 : length(sbj.levels)) {
    fecal <- which(df$obj1 == grp.levels[2] & patientID.p == sbj.levels[j])
    mice <- which(df$obj1 == grp.levels[1] & patientID.p == sbj.levels[j]) 
    mice_not <- which(df$obj1 == grp.levels[1] & patientID.p != sbj.levels[j]) 
    
    ind1 = fecal
    ind2 = mice
    ind3 = mice_not
    
    dist12 <- dist.mat[ind1, ind2]
    dist13 <- dist.mat[ind1, ind3]
    
    dist12.list <- c(dist12.list, mean(dist12))
    dist13.list <- c(dist13.list, mean(dist13))
  }
  stat.perm <- c(stat.perm, mean(dist12.list) - mean(dist13.list))
}
# calculate pvalue
pv <- NULL
if(stat.obj > 0){
  pv <- mean(c(stat.perm >= stat.obj, TRUE))
}else{
  pv <- mean(c(stat.perm <= stat.obj, TRUE))
}

pv  

# generate plot
dist12.list <- dist13.list <- NULL
for (j in 1 : length(sbj.levels)) {
  fecal <- which(df$obj1 == grp.levels[1] & df$PatientID == sbj.levels[j])
  mice <- which(df$obj1 == grp.levels[2] & df$PatientID == sbj.levels[j]) 
  mice_not <- which(df$obj1 == grp.levels[2] & df$PatientID != sbj.levels[j]) 
  
  ind1 = fecal
  ind2 = mice
  ind3 = mice_not
  
  dist12 <- dist.mat[ind1, ind2]
  dist13 <- dist.mat[ind1, ind3]
  
  dist12.list[[j]] <- dist12
  dist13.list[[j]] <- dist13
}

dist12.df <- as.data.frame(unlist(dist12.list)) %>% mutate(group = 'within donor')
dist13.df <- as.data.frame(unlist(dist13.list)) %>% mutate(group = 'between donor')
colnames(dist12.df) <- colnames(dist13.df) <- c('dist','group')
dist.df <- rbind(dist12.df,dist13.df)

ggplot(dist.df, aes(x= group,  y = dist)) +
  stat_boxplot(geom='errorbar', linetype=1, width=0.15)+
  geom_boxplot(width=0.4,outlier.size = 0.5,lwd=0.7) +
  theme_bw()+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(color="black", size = 20),
        axis.text.y = element_text(color="black", size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        axis.title = element_text(color="black", size = 20))+
  labs(x='',y='Bray-Curtis Distance')+
  guides(fill=FALSE) + 
  annotate("text", x = 2.3, y = 0.2,label = paste0('p = ',pv), size = 7)
ggsave(path = 'result/','Fig4b.pdf', width = 6, height = 5, dpi = 300)





##======= Alpha diversity ==========
alpha_div <- estimate_richness(primary.6wks, measures = "Observed") %>% rownames_to_column('SampleID')
df <- as.data.frame(as.matrix(sample_data(primary.6wks))) %>% dplyr::select(c("SampleID","PatientID","PAStatus_","obj1")) %>% inner_join(alpha_div)
df$PAStatus_ <- as.factor(df$PAStatus_)
df$obj1 <- as.factor(df$obj1)
df <- within(df, PAStatus_ <- relevel(PAStatus_, ref = 2))
lm.obj <- lme(Observed ~ PAStatus_, random = ~ 1 | PatientID, df)
summary(lm.obj)
rownames(as.data.frame(coef(summary(lm.obj))))[3]
low_healthy <- round(as.data.frame(coef(summary(lm.obj)))[3,5],3)
p1 <- round(summary(lm.obj)$tTable['PAStatus_Healthy volunteers','p-value'],3)
p2 <- round(summary(lm.obj)$tTable['PAStatus_Low PA PI-IBS','p-value'],3)
pp1 <- ifelse(p1 <= 0.001,'***',ifelse(p1<= 0.01,'**','*'))
pp2 <- ifelse(p2 <= 0.001,'***',ifelse(p2<= 0.01,'**','*'))
df$PAStatus_ <- factor(df$PAStatus_, levels=c("Healthy volunteers",'High PA PI-IBS', 'Low PA PI-IBS'))
ggviolin(df, x = 'PAStatus_', y = 'Observed', fill = 'PAStatus_',bxp.errorbar = TRUE,
         palette = c(brewer.pal(8,'Dark2')[5], '#f79646',brewer.pal(8,'Paired')[2]),
         add = "boxplot", add.params = list(fill = "white"))+ 
  geom_point(data = df, aes(x = PAStatus_,y = Observed,fill = PAStatus_), position=position_jitterdodge(jitter.width =0.3), size =0.9, color = 'black') +
  theme_bw()+
  theme(axis.text = element_text(color="black", size = 24),
        axis.text.x = element_text(color="black", size = 24,angle = 20, hjust = 0.5, vjust = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        axis.title = element_text(color="black", size = 24))+
  labs(x='',y='Observed')+
  scale_y_continuous(limits = c(80, 235)) +
  guides(fill=FALSE) + # healthy to high = 0.001; low to high = 0.006; low to healthy = 0.382
  geom_bracket(xmin = c("Healthy volunteers", "High PA PI-IBS"), xmax = c("High PA PI-IBS", "Low PA PI-IBS"),
               y.position = c(225, 220), size = 1,label.size = 10,bracket.shorten = 0,
               label = c(paste0(pp1),
                         paste0(pp2)))
ggsave(path = 'result/','Fig4c_jitter.pdf', width = 7, height = 6, dpi = 300)



##======= PCoA of mice samples ==========
pslog <- transform_sample_counts(primary.6wks, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
eng <- as.data.frame(out.pcoa.log$vectors)[,c(1:2)] %>% rownames_to_column('SampleID') %>% inner_join(meta.6wks %>% dplyr::select(c('SampleID','PAStatus_')))
pc1 <- round((out.pcoa.log$values$Eigenvalues[1])/sum(out.pcoa.log$values$Eigenvalues),3)
pc2 <- round((out.pcoa.log$values$Eigenvalues[2])/sum(out.pcoa.log$values$Eigenvalues),3)

col.fill = c("High PA PI-IBS" = '#f79646',"Low PA PI-IBS" = brewer.pal(8,'Paired')[2],'Healthy volunteers' = brewer.pal(8,'Dark2')[5])
ggplot(eng) +
  geom_point(aes(x = Axis.1, y = Axis.2,fill = PAStatus_), size = 4, shape = 21)+
  geom_mark_ellipse(aes(x = Axis.1, y = Axis.2, color= PAStatus_))+
  theme_bw() +
  scale_fill_manual(values = col.fill)+
  scale_color_manual(values = col.fill)+
  scale_x_continuous(limits = c(min(eng$Axis.1) * 1.2,max(eng$Axis.1) * 1.2))+
  scale_y_continuous(limits = c(min(eng$Axis.2) * 1.2,max(eng$Axis.2) * 1.2))+
  labs(fill = '') +
  xlab(paste0('PC1, ',(paste0(pc1*100,'%')))) + ylab(paste0('PC2, ',pc2*100,'%')) + 
  theme(axis.text = element_text(color="black", size = 24),
        axis.text.x = element_text(color="black", size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        axis.title =  element_text(size = 24, color = "black"),
        legend.text = element_text(color="black", size = 24))+
  guides(color = FALSE)
ggsave(path = 'result/','Fig4d.pdf', width = 10, height = 6, dpi = 300)


##======= taxonomy barplot(phylum, class, family) ===========
primary.6wks.t <- transform_sample_counts(primary.6wks, function(x) x / sum(x))
meta.6wks <- as.data.frame(as.matrix(sample_data(primary.6wks.t)))
otu_df <- tibble::rownames_to_column(as.data.frame(otu_table(primary.6wks.t)@.Data),'OTU')
tax_df <- tibble::rownames_to_column(as.data.frame(tax_table(primary.6wks.t)@.Data),'OTU')
merge <- merge(otu_df,tax_df,'OTU')
df <- merge[,c('OTU',"Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain",colnames(primary.6wks@otu_table@.Data))]
levels <- rank_names(primary.6wks.t)

means <- list()
getwd()
color1 <- c(brewer.pal(12,'Paired'), brewer.pal(12,'Set3'),brewer.pal(8,'Accent'),brewer.pal(8,'Dark2'))
level <- levels[5]
for (level in levels[c(2:3,5)]) {
  data <- merge %>% 
    dplyr::select(-c(levels[(levels) != level],'OTU')) %>% 
    group_by(!!as.name(level)) %>% 
    dplyr::summarise_each(funs(sum)) %>%
    gather(SeqsID, value, -!!as.name(level)) %>% spread(!!as.name(level),value) %>%
    inner_join((meta.6wks %>% dplyr::select(PAStatus_, SeqsID)), by='SeqsID') %>%
    column_to_rownames('SeqsID')

  data <- aggregate(. ~ PAStatus_, data, function(x) mean(x[!is.na(x)])) %>%
    unite('id',PAStatus_, sep = '.') %>% 
    gather(!!as.name(level), value, -id) %>% spread(id, value)

  data1  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) > 0, ]
  data2  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) <= 0, ]
  
  data2 <- data2 %>% dplyr::select(-c(level))
  data2 <- as.data.frame(t(colSums(data2)))
  rownames(data2) <- 'Others'
  data2 <- data2 %>% rownames_to_column(level)
  bardata <- rbind(data1, data2) %>% melt(.) %>% dplyr::rename(PAStatus_ = variable) %>% droplevels() %>% na.omit()
  bardata$PAStatus_ <- factor(bardata$PAStatus_, levels = unique(bardata$PAStatus_))

  sort <- bardata %>% filter(PAStatus_ == bardata$PAStatus_[2])
  sort <- sort[order(sort$value),][[level]]
  bardata[[level]] = factor(bardata[[level]], levels = sort)
  p1 = ggplot(bardata, aes(x= PAStatus_, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black") +
    scale_fill_brewer(palette = 'Paired') +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'))
  legend <- cowplot::get_legend(p1)

  grid = grid.arrange(legend)
  ggsave(plot = grid,filename = paste0('result/','Fig4e_',level,'legend.pdf'),width = 3.5, height = 4.5, dpi = 500)
  graphics.off()
  
  p1 = ggplot(bardata, aes(x= PAStatus_, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.4, colour="black", size = 0.25) +
    scale_fill_brewer(palette = 'Paired') +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'),
          legend.position = 'none')
  ggsave(path = 'result/',paste0('Fig4e',level, '.pdf'), width = 3.5, height = 4.5, dpi = 300)
  
  
}


##======= KEGG heatmap ======= 
load('data/healthy_high_KEGG.RData')
pdf(file = 'result/FigS6d.pdf', width = 8, height = 8)
heatplot(high_healthy,method = "ward.D2",dend = "row",
         margins=c(5,21),keysize=0.5,key.par = list(cex=0.5),
         cexRow=0.8, cexCol=1)
dev.off()

##======= Differential abundance taxa =======
library(readr)
level2 <- read.csv("data/PAStatus_PI-IBS (Low PA)PI-IBS (High PA)AsBaseStrain_Fig3e_0.002.csv")
level1 <- read.csv("data/PAStatus_Healthy volunteers (Low PA)PI-IBS (High PA)AsBaseStrain_Fig3e_0.002.csv")
max(level2[level2$adjust.p<= 0.1,'adjust.p'])
# -------base =  PI-IBS(High PA), VS  Healthy volunters
max <- max(level1$fc.level1)
min <- min(level1$fc.level1)
head(level1)
ggplot(level1, aes(y=log10P, x=fc.level1,color = color, alpha = color)) +
  geom_point(size =3,pch = 17) +
  scale_color_manual(values=c("#999999", "#C32148"))+
  scale_alpha_manual(values = c(0.5, 1))+
  scale_x_continuous(limits = c(-max(abs(min), max),max(abs(min), max)))+
  geom_hline(yintercept = -log10(0.1), linetype="dashed",size=0.5) +
  geom_text(aes(label=sig.label),hjust=0.45, vjust=-1,size = 3, color = 'grey10', fontface = 3) +
  theme_bw() + 
  theme(axis.title = element_text(color="black", size = 20),
        axis.text = element_text(color="black", size = 20),
        legend.position = "none",
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10))+
  labs(y=expression(-log[10]*P),x='Effect size')+
  guides(fill=FALSE)
ggsave(path = 'result/', 'Fig4f.pdf', width = 8, height = 8, dpi = 300)
# -------base =  PI-IBS(High PA), VS  Low
head(level2)
max <- max(level2$fc.level2)
min <- min(level2$fc.level2)
ggplot(level2, aes(y=log10P, x=fc.level2,color = color, alpha = color)) +
  geom_point(size = 5,pch = 17) +
  scale_color_manual(values=c("#999999", "#C32148"))+
  scale_alpha_manual(values = c(0.5, 1))+
  scale_x_continuous(limits = c(-max(abs(min), max),max(abs(min), max)))+
  scale_y_continuous(limits = c(0, 15))+
  geom_hline(yintercept = -log10(0.1), linetype="dashed",size=0.5)+
  geom_text(aes(label=sig.label),hjust=0.45, vjust=-1,size = 3, color = 'grey10', fontface = 3) +
  theme_bw()+ 
  theme(axis.text = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10))+
  labs(y=expression(-log[10]*P),x='Effect size')+
  guides(fill=FALSE)
ggsave(path = 'result/', 'Fig4g.pdf', width = 8, height = 8, dpi = 300)

##======= Differential abundance KEGG ========
load('data/join_kegg.modules.RData')
kegg <- apply(join_kegg.modules,2,function(x){x/sum(x)})
# filter
kegg <- kegg[rowSums(kegg!=0) > ncol(kegg) *0.1, , drop=FALSE]	
kegg <- kegg[rowMaxs(kegg) > 0.002, , drop=FALSE]
# asin sqrt transform
join_kegg.modules.r <- asin(sqrt(kegg)) %>%
  as.data.frame(.) %>%
  rownames_to_column('kegg.name') %>%
  mutate(keggid = paste0('kegg',1:nrow(.)))
join_kegg.modules.r0 <- join_kegg.modules.r %>% 
  dplyr::select(-kegg.name) %>% 
  column_to_rownames('keggid') %>% as.matrix()
kegg_id <- join_kegg.modules.r %>% dplyr::select(kegg.name,keggid)

kegg_join <- as.data.frame(t(join_kegg.modules.r0)) %>% rownames_to_column('SampleID')%>% 
  inner_join(meta.6wks %>% dplyr::select(SampleID,PAStatus_, PatientID)) %>% 
  mutate(keggid = paste0('kegg',1:nrow(.))) %>% 
  dplyr::select(-keggid)

kegg_join$PAStatus_ <- as.factor(kegg_join$PAStatus_)
kegg_join$PAStatus_ <- relevel(kegg_join$PAStatus_, ref = 'High PA PI-IBS')
kegg_join$PatientID <- as.factor(kegg_join$PatientID)
data <- kegg_join %>% column_to_rownames('SampleID')
head(data)
pvs.level1 <-  pvs.level2 <-  fcs.level1 <- fcs.level2 <- c()
for (i in 1:(ncol(data)-2)) {
  tryCatch({
    m1.nb <- NULL
    m1.nb <- lme(as.formula(paste0(colnames(data)[i], ' ~ PAStatus_')), random = (~ 1 | PatientID), data)
    pv.level1 <- wald.test(b = fixed.effects(m1.nb), Sigma = vcov(m1.nb), Terms = 2)$result$chi2['P']
    names(pv.level1) <- colnames(data)[i]
    pvs.level1 <- c(pvs.level1, pv.level1)
    
    pv.level2 <- wald.test(b = fixed.effects(m1.nb), Sigma = vcov(m1.nb), Terms = 3)$result$chi2['P']
    names(pv.level2) <- colnames(data)[i]
    pvs.level2 <- c(pvs.level2, pv.level2)
    
    name.obj <- names(coef(m1.nb))[2:3]
    coef.nb <- fixed.effects(m1.nb)
    fc.level1 <- coef.nb['PAStatus_Healthy volunteers']
    fc.level2 <- coef.nb['PAStatus_Low PA PI-IBS']
    names(fc.level1) <- colnames(data)[i]
    names(fc.level2) <- colnames(data)[i]
    
    fcs.level1 <- c(fcs.level1,fc.level1)
    fcs.level2 <- c(fcs.level2,fc.level2)
    
  }, error=function(e){cat(paste0(colnames(data)[i]," ERROR :"),conditionMessage(e), "\n")})
}


pvs.level1.df <- 
  data.frame((pvs.level1)) %>% 
  rownames_to_column() %>% 
  dplyr::select(c(1:2)) %>%
  rename_at(colnames(.), list( ~ c('taxon','p-value.level1'))) %>% 
  mutate(taxon = gsub('.taxon.*','',taxon)) %>% 
  mutate(adjust.p = p.adjust(`p-value.level1`, method = 'fdr')) %>%
  dplyr::select(-`p-value.level1`)


pvs.level2.df <- 
  data.frame((pvs.level2)) %>% 
  rownames_to_column() %>% 
  dplyr::select(c(1:2)) %>%
  rename_at(colnames(.), list( ~ c('taxon','p-value.level2'))) %>% 
  mutate(taxon = gsub('.taxon.*','',taxon)) %>% 
  mutate(adjust.p = p.adjust(`p-value.level2`, method = 'fdr')) %>% 
  dplyr::select(-`p-value.level2`)

fcs.level1.df <- 
  data.frame((fcs.level1)) %>% 
  rownames_to_column() %>% 
  dplyr::select(c(1:2)) %>%
  rename_at(colnames(.), list( ~ c('taxon','fc.level1'))) %>% 
  mutate(taxon = gsub('.PAStatus.*','',taxon))

fcs.level2.df <- 
  data.frame((fcs.level2)) %>% 
  rownames_to_column() %>% 
  dplyr::select(c(1:2)) %>%
  rename_at(colnames(.), list( ~ c('taxon','fc.level2'))) %>% 
  mutate(taxon = gsub('.PAStatus.*','',taxon))

level1 <- inner_join(pvs.level1.df, fcs.level1.df) %>% 
  mutate(color = ifelse(adjust.p < 0.1, "sig", "nonsig")) %>% 
  mutate(log10P = -log10(adjust.p)) %>%
  mutate(keggid = taxon) %>%
  inner_join(kegg_id %>% dplyr::select(c(kegg.name,keggid))) 
level1 <- level1 %>% mutate(name = ifelse(adjust.p < 0.1, kegg.name, ""))
level2 <- inner_join(pvs.level2.df, fcs.level2.df) %>% 
  mutate(color = ifelse(adjust.p < 0.1, "sig", "nonsig")) %>% 
  mutate(log10P = -log10(adjust.p)) %>%
  mutate(keggid = taxon) %>%
  inner_join(kegg_id %>% dplyr::select(c(kegg.name,keggid))) 
level2 <- level2 %>% mutate(name = ifelse(adjust.p < 0.1, kegg.name, ""))

base <- levels(data$PAStatus_)[!levels(data$PAStatus_) %in% gsub('PAStatus_','',name.obj)]

name_level1 <- c(name.obj[1], base)
name_level1 <- gsub('PAStatus_','',name_level1)
name_level2 <- c(name.obj[2], base)
name_level2 <- gsub('PAStatus_','',name_level2)

max <- max(level1$fc.level1)
min <- min(level1$fc.level1)
ggplot(level1, aes(y=log10P, x=fc.level1,color = color, alpha = color)) +
  geom_point(size = 5,pch = 17) +
  scale_color_manual(values=c("#999999", "#C32148"))+
  scale_alpha_manual(values = c(0.5, 1))+
  scale_x_continuous(limits = c(-max, max))+
  # scale_y_continuous(limits = c(0, 15))+
  geom_hline(yintercept = -log10(0.1), linetype="dashed",size=0.5)+
  geom_text(aes(label=name),hjust=0.45, vjust=-1,size = 3, color = 'grey10', fontface = 3) +
  theme_bw()+ 
  theme(text = element_text(size = 20),
        axis.text = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        plot.title = element_text(size = 10))+
  labs(y=expression(-log[10]*P),x='Effect size')+
  guides(fill=FALSE) 
ggsave(path = 'result/', 'Fig4i.pdf', width = 8, height = 8, dpi = 300)



##======= random forest prediction =========
levels <- rank_names(primary.6wks)
level <- levels[8]
wk6 <- primary.6wks %>% 
  subset_samples(PAStatus_ %in% c("High PA PI-IBS","Low PA PI-IBS")) %>% 
  tax_glom(level) 
wk6_ <- wk6 %>% transform_sample_counts(function(x) x / sum(x))
wk6_ <- filter_taxa(wk6_, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
wk6_ <- filter_taxa(wk6_, function(x) max(x) > 0.002, TRUE)
load('data/otu.name.RData')
df <- as.data.frame(wk6_@otu_table) %>% rownames_to_column('taxon') %>% 
  inner_join(otu.name %>% rownames_to_column('taxon') %>% dplyr::select(c('taxon'))) %>%
  gather(SampleID, value, -taxon) %>% 
  spread(taxon,value) 
df <- df %>% inner_join(as.data.frame(as.matrix(wk6_@sam_data)) %>% dplyr::select(c('SampleID','PAStatus_','PatientID')))%>% column_to_rownames('SampleID')
df$PAStatus_ <- as.factor(df$PAStatus_)
group_names <- as.vector(unique(df$PatientID))

prob.list <- accuracy.list <- NULL
if (level == 'Strain'){
  seed = 1121
}else{
  seed = 1122
}
set.seed(seed)# strain seed = 1121
for (i in 1:500){
  g = sample(group_names, 1, replace = TRUE)
  train_data <- df[df['PatientID']!=g, ][,-(ncol(df))]
  test_data <- df[df['PatientID']==g, ][,-(ncol(df))]
  data.rf <- randomForest(PAStatus_~ ., data=train_data, ntree=501, importance=TRUE, proximities=TRUE)
  data.pred <- predict(data.rf, test_data[,1:ncol(test_data)-1], type = 'prob')
  accuracy.list[[paste(g)]] <- cbind(accuracy.list[[paste(g)]], data.rf$importance[,3])
  prob.list[[paste(g)]] <- cbind(prob.list[[paste(g)]], data.pred[,1])
}
probs <- c()
accuracys <- c()
labels <- c()
for (g in group_names){
  # mean probability
  accuracy <- rowMeans(accuracy.list[[paste(g)]])
  prob <- rowMeans(prob.list[[paste(g)]])
  cat(length(accuracy))
  cat(length(prob))
  cat(' ')
  probs <- c(probs, prob)
  accuracys <- cbind(accuracys, accuracy)
  
  # label
  label <- as.vector(df[df['PatientID']==g, ][, 'PAStatus_'])
  labels <- c(labels, label)
}
Importance <- as.data.frame(as.matrix(rowMeans(accuracys))) %>% rownames_to_column('taxon') %>%
  inner_join(otu.name %>% rownames_to_column('taxon') %>% dplyr::select(c('taxon',level)))
Importance1 <- Importance %>% filter(V1 >= quantile(Importance$V1,0.50))
Importance1 <- Importance1[order(-Importance1$V1),]
Importance1$cutoff <- c(rep('red',10),rep('gray',nrow(Importance1)-10))
Importance1[Importance1$cutoff == 'gray','cutoff1'] ='#00BFC4'
Importance1[Importance1$cutoff == 'red','cutoff1'] ='#F8766D'
values = c('#F8766D','#00BFC4') 
Importance1 = Importance1[order(-Importance1$V1),] %>% head(10)

ggplot(Importance1, aes(x = reorder(get(level), V1), y = V1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = values) +
  ylab("Mean Decrease in Accuracy (Importance)")+
  xlab('') +
  theme_bw() +
  coord_flip() + 
  theme(axis.text.x = element_text(color="black", size = 20, angle =90, vjust = 0.5, hjust =1),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') 
ggsave(paste0('result/','Fig4h.pdf'),width = 12, height = 6, dpi = 100)



## ROC curve
labels <- as.numeric(as.factor(labels))
labels.list <- probs.list <- list()
labels.list[[level]] <- labels
probs.list[[level]] <- probs
roc <- pROC::roc(labels,probs)
sens.ci <- ci.se(roc,specificities = roc$specificities * ifelse(roc$percent,100, 1))
df <- cbind(roc$sensitivities*100,roc$specificities*100, sens.ci[,c(1,3)] *100) %>% as.data.frame()
colnames(df) <- c('sensitivities','specificities','se.low','se.high')
size = 20

ggplot(df, aes(x = specificities,y = sensitivities))+
  geom_ribbon(aes(ymin = se.low, ymax = se.high), fill = "grey") +
  geom_path(aes(y=se.low), color = 'grey50') +
  geom_path(aes(y=se.high), color = 'grey50') +
  geom_path()+
  theme_bw()+
  xlim(c(100,0))+
  geom_abline(slope=1,intercept=100,color='grey')+
  xlab('Specificity (%)') +
  ylab('Sensitivity (%)')+
  annotate(geom='text', x=45, y=5, label=paste0("AUC:", round(ci(roc)[2]*100,1),'% (',round(ci(roc)[1]*100,1),'%-',round(ci(roc)[3]*100,1),'%)'), 
           size =size/3)+
  theme(axis.text = element_text(color="black", size = size),
        axis.title = element_text(color="black", size = size),
        legend.position = 'none',        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        plot.title = element_text(color="black", hjust = 0.5),
        title =element_text(size=size)) 
ggsave(path = 'result/', 'SupplementaryS6C.pdf', width =5, height =5, dpi = 300)



##======= Donor taxonomy plot ===========
## Humanized mice represent donor human stools with taxa differences observed between low and high PA patients and predictive of PA status (Data used for Fig. 3 Mouse and Human slurry comparisons)
# 1. Comparisons between high, low and healthy slurries and the microbiomes of the recipient mice (phyla, class, family) Bar-charts (JC/LY)
otu_df <- tibble::rownames_to_column(as.data.frame(otu_table(primary.t)@.Data),'OTU')
tax_df <- tibble::rownames_to_column(as.data.frame(tax_table(primary.t)@.Data),'OTU')
merge <- merge(otu_df,tax_df,'OTU')
levels <- rank_names(primary.t)
means <- list()
getwd()
color1 <- c(brewer.pal(12,'Paired'), brewer.pal(12,'Set3'),brewer.pal(8,'Accent'),brewer.pal(8,'Dark2'))
level <- levels[5]
for (level in levels[c(2:3,5)]) {
  data <- merge %>% 
    dplyr::select(-c(levels[(levels) != level],'OTU')) %>% 
    group_by(!!as.name(level)) %>% 
    dplyr::summarise_each(funs(sum)) %>%
    gather(SeqsID, value, -!!as.name(level)) %>% spread(!!as.name(level),value) %>%
    inner_join((meta %>% dplyr::select(PAStatus_, SeqsID)), by='SeqsID') %>%
    column_to_rownames('SeqsID')
  data <- aggregate(. ~ PAStatus_, data, function(x) mean(x[!is.na(x)])) %>%
    unite('id',PAStatus_, sep = '.') %>%
    gather(!!as.name(level), value, -id) %>% spread(id, value)

  data1  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) > 0, ]
  data2  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) <= 0, ]
  
  data2 <- data2 %>% dplyr::select(-c(level))
  data2 <- as.data.frame(t(colSums(data2)))
  rownames(data2) <- 'Others'
  data2 <- data2 %>% rownames_to_column(level)
  bardata <- rbind(data1, data2) %>% melt(.) %>% dplyr::rename(PAStatus_ = variable) %>% droplevels() %>% na.omit()
  bardata$PAStatus_ <- factor(bardata$PAStatus_, levels = unique(bardata$PAStatus_))

  sort <- bardata %>% filter(PAStatus_ == bardata$PAStatus_[2])
  sort <- sort[order(sort$value),][[level]]
  bardata[[level]] = factor(bardata[[level]], levels = sort)
  cols <- c(RColorBrewer::brewer.pal(12, 'Paired'),RColorBrewer::brewer.pal(8, 'Set2'))
  
  p1 = ggplot(bardata, aes(x= PAStatus_, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black") +
    scale_fill_manual(values = cols) +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'))
  legend <- cowplot::get_legend(p1)
  # pdf(paste0('result/','Fig3d_',level,'legend.pdf'),width = 3.5, height = 4.5)
  
  grid = grid.arrange(legend)
  ggsave(plot = grid,filename = paste0('result/','FigS5_',level,'legend.pdf'),width = 3.5, height = 5.5, dpi = 500)
  graphics.off()
  
  p1 = ggplot(bardata, aes(x= PAStatus_, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.4, colour="black", size = 0.25) +
    scale_fill_manual(values = cols) +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'),
          legend.position = 'none')
  ggsave(paste0('result/','FigS5_',level,'.pdf'),width = 6, height = 5, dpi = 500)
  
}
















#============ fmt alpha diveristy ========
load('data/FMT.RData')
fmt_control <- FMT %>% subset_samples(FMT1 %in% c('control mice','post FMT humanized mice'))
x <- estimate_richness(fmt_control, measures = c("Observed"))
df <- meta_FMT %>% dplyr::select(c("SampleID",'FMT1')) %>% filter(FMT1 %in% c('control mice','post FMT humanized mice')) ;table(df$FMT1)
x$SampleID <- rownames(x)
df_ <- inner_join(df, x) 
df_$FMT1 <- as.factor(df_$FMT1)
df_ <- within(df_, grp <- relevel(FMT1, ref = 1))
lm.obj <- lm(Observed ~ grp, df_)
fmt1 <- as.data.frame(coef(summary(lm.obj)))

df_$grp <- factor(df_$grp, levels=c('post FMT humanized mice','control mice'))
df_$grp <- gsub('post FMT humanized mice','Post-FMT',df_$grp)
df_$grp <- gsub('control mice','Post-control gavage',df_$grp)

ggviolin(df_, x = 'grp', y = 'Observed', fill = 'grp',bxp.errorbar = TRUE, remove = T,
         palette = c(`Post-control gavage`=brewer.pal(8,'Paired')[2], `Post-FMT`='#f79646'),
         add = "boxplot", add.params = list(fill = "white"))+ 
  geom_point(data = df_, aes(x = grp,y = Observed,fill = grp), position=position_jitterdodge(jitter.width =0.2), size =2, color = 'black') +
  theme_bw()+
  labs(x='',y='Observed') +
  scale_y_continuous(limits = c(min(df_$Observed) * 0.88,max(df_$Observed) * 1.11)) +
  guides(fill=FALSE) + 
  theme(axis.text.x = element_text(color="black", size = 24),
        axis.text.y = element_text(color="black", size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        axis.title = element_text(color="black", size = 24)) +
  geom_bracket(xmin = 1, xmax = 2,size = 1,label.size = 9,
               #tip.length = 0,
               y.position = c(max(df_$Observed) * 1.1),
               label = '*')
ggsave(path = 'result/','Fig5c_jitter.pdf', width = 7, height = 6.5, dpi = 300)


##======== PCoA of FMT ==============
pslog <- transform_sample_counts(FMT, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
eng <- as.data.frame(out.pcoa.log$vectors)[,c(1:2)] %>% rownames_to_column('SampleID') %>% inner_join(meta_FMT %>% dplyr::select(c('SampleID','FMT1')))
head(eng)
pc1 <- round((out.pcoa.log$values$Eigenvalues[1])/sum(out.pcoa.log$values$Eigenvalues),2)
pc2 <- round((out.pcoa.log$values$Eigenvalues[2])/sum(out.pcoa.log$values$Eigenvalues),2)
table(eng$FMT1)
eng$FMT1 <- sub('recipient humanized mice','Baseline (High PA)',eng$FMT1)
eng$FMT1 <-sub('post FMT humanized mice','Post-FMT',eng$FMT1)
eng$FMT1 <-sub('control mice','Post-control gavage',eng$FMT1)
eng$FMT1 <-sub('donor human','Healthy volunteers',eng$FMT1)
eng$FMT1 <-sub('recipient human','High PA PI-IBS',eng$FMT1)

eng$shape = 'Mice'
eng[eng$FMT1 %in% c('Healthy volunteers','High PA PI-IBS'),'shape'] = 'Human'
eng[eng$FMT1 %in% c('Post-control gavage'),'shape'] = 'control'
table(eng$shape);

eng$color = 'Mice'
eng[eng$FMT1 == 'Post-FMT','color'] = '#f79646'
eng[eng$FMT1 == 'Baseline (High PA)','color'] = brewer.pal(8,'Paired')[2]
eng[eng$FMT1 == 'Post-control gavage','color'] = brewer.pal(8,'Paired')[2]
eng[eng$FMT1 == 'Healthy volunteers','color'] = brewer.pal(8,'Dark2')[5]
eng[eng$FMT1 == 'High PA PI-IBS','color'] = brewer.pal(8,'Paired')[2]

p1 = ggplot(eng) +
  geom_point(aes(x = Axis.1, y = Axis.2,color = color, shape = shape, fill = color), size = 4)+
  theme_bw() +
  scale_color_manual(values = c(`#f79646`="#f79646",`#1F78B4`="#1F78B4",`#66A61E`="#66A61E")) +
  scale_shape_manual(values = c(Human=19,Mice=24, control = 2)) +
  scale_fill_manual(values = c(`#f79646`="#f79646",`#1F78B4`="#1F78B4",`#66A61E`="#66A61E")) +
  labs(fill = 'Mice', shape = 'Donor') +
  xlab(paste0('PC1, ',(paste0(pc1*100,'%')))) + 
  ylab(paste0('PC2, ',paste0(pc2*100, '%'))) + 
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 20, color = "black"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none')
l1 = ggplot(eng %>% filter(shape =='Human'))+ 
  geom_point(aes(x = Axis.1, y = Axis.2,shape = shape, fill=FMT1, color = FMT1), size=4, shape = 19) +
  theme_bw() +
  scale_color_manual(values = c(`High PA PI-IBS` = '#1F78B4', `Healthy volunteers` = '#66A61E'))+
  scale_fill_manual(values = c(`High PA PI-IBS` = '#1F78B4', `Healthy volunteers` = '#66A61E'))+
  # scale_shape_manual(values = c(19)) +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  labs(color = 'Human') +
  guides(shape = 'none',fill = 'none')
l1= get_legend(l1)
grid.draw(l1,recording=F)

col.fill = c("Post-FMT" = '#f79646',"Baseline (High PA)" = brewer.pal(8,'Paired')[2],'Post-control gavage' = 'white')
col.shape = c("Post-FMT" = '#f79646',"Baseline (High PA)" =brewer.pal(8,'Paired')[2],'Post-control gavage' = brewer.pal(8,'Paired')[2])

l2 = ggplot(eng %>% filter(shape !='Human') %>% dplyr::rename(Mice = FMT1))+ 
  geom_point(aes(x = Axis.1, y = Axis.2, shape = Mice, fill=Mice, color = Mice), size=4, shape = 24) +
  theme_bw()+
  scale_color_manual(values = col.shape)+
  scale_fill_manual(values = col.fill)+
  # scale_shape_manual(values = c(24,24,24)) +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))

l2= get_legend(l2)
grid.draw(l2,recording=F)

ll = ggarrange(NULL,l1, l2,NULL, nrow = 4, ncol = 1)
pp =  ggarrange(p1,ll,  widths = c(2, 1),nrow = 1, ncol = 2)
ggsave(path = 'result/','Fig5d.pdf', width =10, height = 6, dpi = 300)



##======== Differential abundance ========
taxon_level <- 'Species'
load('data/FMT_diff.RData')
pdf(paste0('result/Fig5e_Species.pdf'), width = 8, height = 8)
heatplot(data,scale="none",method = "ward.D2",zlim = c(-3,3),
         margins=c(10,15),keysize=0.6,key.par = list(cex=0.5),
         labRow=as.expression(lapply(rownames(data), function(a) bquote(italic(.(a))))),
         cexRow=1, cexCol=0.5)
dev.off()

# boxplot
box_bar = data %>% rownames_to_column(taxon_level) %>% melt() %>% group_by(!!as.name(taxon_level), variable) %>% as.data.frame()
box_bar$variable = gsub('.S.*','',box_bar$variable)
head(box_bar);table(box_bar$variable)
box_bar$variable = gsub('\\..*','',box_bar$variable)
idx <- grep('Post-control gavage', box_bar$variable)
sort <- box_bar[idx,]
sort <- aggregate(value ~ Species, data = sort, function(x) median(x))
sort <- as.data.frame(sort[order(sort$value, decreasing = T),])[,1]
box_bar[,taxon_level] = factor(box_bar$Species, levels = unique(sort))
head(box_bar)
ggplot(box_bar, aes(x = !!as.name(taxon_level), y = value, fill = variable)) +
  geom_boxplot(width=0.8,outlier.size = 0.4,lwd=0.4) +
  scale_y_continuous(trans = sqrt_trans(),
                     # limits = c(0,0.05),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) +
  theme_bw() +
  # scale_fill_brewer(palette = 'Dark2') +
  scale_fill_manual(values = c("Post-FMT" = '#f79646','Post-control gavage' = brewer.pal(8,'Paired')[2]))+
  theme(axis.text.x = element_text(color="black", size = 18, angle = 90, hjust = 1, vjust = 0.25, face = 'italic'),
        axis.text.y = element_text(color="black", size = 18),
        axis.title = element_text(color="black", size = 18),
        legend.text=element_text(size=18),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=18))+
  labs(x='',y='Proportion',fill = '')
ggsave(path = 'result/','Fig5g.pdf', width =15, height = 8, dpi = 300)

## ====== KEGG volcano plot ======
level1 <- read.csv("data/Fig4f.csv")
level1$name[level1$adjust.p<= 0.1] <- level1$kegg.name[level1$adjust.p<= 0.1]
max <- max(level1$fc.level1)
min <- min(level1$fc.level1)
ggplot(level1, aes(y=log10P, x=fc.level1,color = color, alpha = color)) +
  geom_point(size = 5,pch = 17) +
  scale_color_manual(values=c("#999999", "#C32148"))+
  scale_fill_manual(values=c("#999999", "#C32148"))+
  scale_alpha_manual(values = c(0.5, 1))+
  scale_x_continuous(limits = c(min, max))+
  geom_hline(yintercept = -log10(0.1), linetype="dashed",size=0.5)+
  geom_text(aes(label=name),hjust=0.45, vjust=1,size =3, color = 'grey10') +
  theme_bw()+ 
  theme(axis.text = element_text(color="black", size = 20),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        axis.title = element_text(color="black", size = 20),legend.position = "none")+
  labs(y=expression(-log[10]*P),x='Effect size')+
  guides(fill=FALSE)
ggsave('result/Fig5f.pdf', width = 8, height = 7, dpi = 100)



##========= barplot of taxonomy ==========
t_rary <- transform_sample_counts(fmt_control, function(x) x/sum(x))
otu <- as.data.frame(otu_table(t_rary)@.Data)
phy.0 <- as.data.frame(tax_table(t_rary)@.Data)
otu_df <- tibble::rownames_to_column(otu,'OTU')
tax_df <- tibble::rownames_to_column(phy.0,'OTU')
merge <- merge(otu_df,tax_df,'OTU')
levels <- rank_names(t_rary)
meta.dat <- meta_FMT %>% filter(FMT1 %in% c('control mice','post FMT humanized mice'))
color1 <- c(brewer.pal(12,'Paired'), brewer.pal(12,'Set3'),brewer.pal(8,'Accent'),brewer.pal(8,'Dark2'))
level <- levels[2]
for (level in levels[c(2:3,5)]) {
  data <- merge %>% 
    dplyr::select(-c(levels[(levels) != level],'OTU')) %>% 
    group_by(!!as.name(level)) %>% 
    dplyr::summarise_each(funs(sum)) %>%
    gather(SeqsID, value, -!!as.name(level)) %>% spread(!!as.name(level),value) %>%
    inner_join((meta.dat %>% dplyr::select(FMT1, SeqsID)), by='SeqsID') %>%
    column_to_rownames('SeqsID') %>% na.omit()
  data <- aggregate(. ~ FMT1, data, function(x) mean(x[!is.na(x)]))
  data <- data %>%
    unite('id',FMT1, sep = '.') %>% 
    gather(!!as.name(level), value, -id) %>% spread(id, value)
  data1  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) > 0, ]
  data2  <- data[rowSums((data %>% dplyr::select(-c(level))) >= 0.01) <= 0, ]
  
  data2 <- data2 %>% dplyr::select(-c(level))
  data2 <- as.data.frame(t(colSums(data2)))
  rownames(data2) <- 'Others'
  data2 <- data2 %>% rownames_to_column(level)
  bardata <- rbind(data1, data2) %>% melt(.) %>% dplyr::rename(FMT1 = variable) %>% droplevels() %>% na.omit()
  bardata$FMT1 <- factor(bardata$FMT1, levels = unique(bardata$FMT1))
  
  sort <- bardata %>% filter(FMT1 == bardata$FMT1[2])
  sort <- sort[order(sort$value),][[level]]
  bardata[[level]] = factor(bardata[[level]], levels = sort)
  bardata$FMT1 <- gsub('post FMT humanized mice','Post-FMT',bardata$FMT1)
  bardata$FMT1 <- gsub('control mice','Post-control gavage',bardata$FMT1)
  
  p1 = ggplot(bardata, aes(x= FMT1, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.5, colour="black") +
    scale_fill_brewer(palette = 'Paired') +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'))
  legend <- cowplot::get_legend(p1)
  grid = grid.arrange(legend)
  ggsave(plot = grid,filename = paste0('result/','FigS7_',level,'legend.pdf'),width = 3.5, height = 4.5, dpi = 500)
  graphics.off()
  p1 = ggplot(bardata, aes(x= FMT1, y=value, fill = get(level))) +
    geom_bar(stat="identity", width=0.4, colour="black", size = 0.25) +
    # scale_fill_brewer(palette = 'Paired') +
    scale_fill_manual(values = c(brewer.pal(11, 'Paired'),brewer.pal(8, 'Set1'))) +
    labs(x = '', y = "Proportion", fill = level)+
    theme_classic() +
    theme(text = element_text(size = 16, color = "black"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(0.1, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(color="black", size = 16, angle = 30, vjust = 0.7, hjust = 0.6),
          axis.text.y = element_text(color="black", size = 16),
          legend.text=element_text(size=16, face = 'italic'),
          legend.title=element_text(size=16),
          strip.text.x = element_text(angle=0, size = 16),
          panel.spacing = unit(0.02, "cm"),
          strip.background = element_rect(size = 0.5, color = 'black',fill = 'white'),
          legend.position = 'none')
  ggsave(paste0('result/','FigS7_',level,'.pdf'),width = 6, height = 5, dpi = 500)
  
}


load('data/FMT_t.RData')
## 2. Correlation scatterplots between PA and taxa abundance post FMT for taxa that were differentially abundant between FMT and control treated mice after the second FMT (JC/LY)
Fig4e <- read_delim("data/Fig4e.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>% dplyr::filter(adjust.p <= 0.1) %>% dplyr::select(c('taxon','Species')) %>% column_to_rownames('taxon')
sp <- as.data.frame(as.matrix(FMT_t@otu_table))[rownames(Fig4e),]
tb <- merge(sp,Fig4e, by = 0) %>% column_to_rownames('Species') %>% dplyr::select(-Row.names) %>% t() %>% as.data.frame() %>% rownames_to_column('SampleID')
met <- meta_FMT %>% filter(FMT1 =='post FMT humanized mice') %>% dplyr::mutate(PA = log(as.numeric(as.character(PA)))) %>% dplyr::select(c('SampleID','FMT1','PA'))
tb1 <- as.data.frame(tb) %>% inner_join(as.data.frame(met)) %>% dplyr::select(-FMT1) %>% column_to_rownames('SampleID')
names <- colnames(tb1)[colnames(tb1) != 'PA']

p.vec <- r.vec <- NULL
for(name in names){
  df = tb1[,c(name, 'PA'),drop =F]
  res <- cor.test(df[,name], df$PA, method = 'spearman')
  p <- res$p.value
  names(p) <- name
  r <- res$estimate
  names(r) <- name
  p.vec <- c(p.vec, p)
  r.vec <- c(r.vec, r)
}
q.vec <- p.adjust(p.vec, method = 'fdr')

spearman.PostFMT.logPA <- as.data.frame(cbind(p.vec, q.vec, r.vec))
colnames(spearman.PostFMT.logPA) = c('Pvalue','Qvalue','R')

## extract significant differential taxa between post-FMT and post control gavage
Fig4e.sig <- read_delim("data/Fig4e.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>% filter(adjust.p <= 0.1) %>%
  dplyr::select(c('taxon','Species')) %>% column_to_rownames('taxon')
spearman <- spearman.PostFMT.logPA %>% dplyr::filter(Qvalue <= 0.1)
tb <- as.data.frame(as.matrix(FMT_t@otu_table)) %>% 
  merge(as.data.frame(as.matrix(FMT_t@tax_table))[,'Species',drop =F],by = 0) %>% 
  dplyr::select(-Row.names) %>% column_to_rownames('Species') %>% 
  t() %>% as.data.frame() 
met <- meta_FMT %>% filter(FMT1 =='post FMT humanized mice') %>% dplyr::mutate(PA = log(as.numeric(as.character(PA)))) %>% dplyr::select(c('SampleID','FMT1','PA'))
tb1 <- as.data.frame(tb%>% rownames_to_column('SampleID')) %>% inner_join(as.data.frame(met)) %>% dplyr::select(-FMT1) %>% column_to_rownames('SampleID')
tb2 <- tb1[,c(rownames(spearman),'PA')]
names <- colnames(tb2)[colnames(tb2) !='PA']

variables = c()
for(i in 1:nrow(spearman)){
  q = format(round(spearman$Qvalue[i],3), scientific = F);r = format(round(spearman$R[i],2), scientific = F)
  variables = c(variables,paste0('r=',r,',q=',q))
}

ann_text <- data.frame(PA = c(rep(8,length(names))),value = c(rep(0.01, length(names))),lab = variables,
                       variable = rownames(spearman))
names <- colnames(tb2)[colnames(tb2) != 'PA']
for(name in names){
  scatter = tb2[,c(name,'PA'),drop=F] %>% rownames_to_column('SampleID')
  cat(name,'\n')
  p = ggplot(scatter, aes(x = PA,y = get(name))) +
    geom_point(cex =5, shape = 17, color = '#f79646') +
    geom_smooth(method=lm, color = 'grey50') +
    labs(y = name, x = 'log(PA)', color = '') +
    guides(fill=guide_legend("")) +
    theme_classic() +
    theme(text = element_text(size = 24, color = 'black'),
          axis.text = element_text(size = 24, color = 'black'),
          axis.title = element_text(size = 24, color = 'black'),
          axis.title.y = element_text(size = 24, color = 'black'),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24),
          legend.position = 'bottom',
          strip.text = element_text(size = 24, color = 'black', face = 'italic'),
          strip.background = element_blank(),
          strip.placement = "outside") +
    geom_text(
      data  = ann_text[ann_text$variable ==name,,drop=F], size = 9,
      mapping = aes(x = -Inf, y = -Inf, label = lab),
      hjust   = -0.2,
      vjust   = -10
    )
  ggsave(path = 'result/', paste0('FigS8_',name,'.pdf'), width =6, height =6, dpi = 100)
}






