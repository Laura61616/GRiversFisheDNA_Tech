# READNE: summary dataset, 20.11.2024
# Created by Yan Zhang (Nanjing University, China)
# ALL RIGHTS RESERVED
wd <- '../data/'
setwd(wd)

fun_path <- '../code/00_functions.R' ## CHANGE ME
source(fun_path)
useCors <- c('#2EC4B6', '#0A2845', '#E4A9CE', '#452E68','#4F5DCE', '#F4D35E')
##Please change the work space in the 00_functions.R
#############################################################
####local species list verus primer DB taxonomy
#############################################################
info <- read.csv('00PrimerList.csv', row.names = 1)
head(info)

data_info <- read.csv('00DataList.csv')
data_info <- unique(data_info[,c(1,4)])
rownames(data_info) <- data_info$DataID

all_record <- readRDS('Dataset_Fish_Recored_C.RDS')

IDs = rownames(info)

result <- NULL
for (ID in IDs) {
  primer = info[ID,1]
  DataID = info[ID,2]

  global_taxonomy <- read.table(paste0('db_taxonomy/', primer,'_taxonomy.txt'), sep = '\t',header = F)
  global_fish_taxonomy <- global_taxonomy[grep('Actinopteri', global_taxonomy$V2),]

  global_fish_spes <- gsub('.*;', '',global_fish_taxonomy$V2)
  global_fish_gens <- gsub(' .*', '',global_fish_spes)
  
  ###species list: Basin
  spelist_basin <- unlist(all_record[[DataID]]$B)
  spelist_basin_gen <- gsub(' .*','',spelist_basin)
  ###species list: Country
  spelist_country <- unlist(all_record[[DataID]]$C)
  spelist_country_gen <- gsub(' .*','',spelist_country)
  ###species list: Ecoregion
  spelist_ecoregion <- unlist(all_record[[DataID]]$E)
  spelist_ecoregion_gen <- gsub(' .*','',spelist_ecoregion)
  ###summary
  add <- data.frame(
    RawSeqID = ID,
    DataID = DataID,
    Continent = data_info[DataID,2],
    Primer = primer,
    Barcode = info[ID,7],
    Level = rep(c('Basin','Country','Ecoregion'), each = 2),
    Resolution = rep(c('Species', 'Genus'),3),
    Total_count = rep(c(length(spelist_basin),length(spelist_country),length(spelist_ecoregion)),each=2),
    Identified_count = c(sum(spelist_basin%in%global_fish_spes),
                         sum(spelist_basin_gen%in%global_fish_gens),
                         sum(spelist_country%in%global_fish_spes),
                         sum(spelist_country_gen%in%global_fish_gens),
                         sum(spelist_ecoregion%in%global_fish_spes),
                         sum(spelist_ecoregion_gen%in%global_fish_gens)))
  add$Prop = add$Identified_count/add$Total_count
  result <- rbind.data.frame(result, add)
}
write.table(result,'Res_FishRecord_primer_db.csv', row.names = F, sep = ',', quote = F)
##plotting
result <- read.csv('Res_FishRecord_primer_db.csv')
head(result)
p1 <- ggplot(result, aes(Primer, Prop, fill = Primer, color = Primer))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Database Coverage of local species')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1

p2 <- ggplot(result, aes(Barcode, Prop, fill = Barcode, color = Barcode))+
  facet_grid(Level~Resolution)+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_boxplot(fill = 'Transparent')+
  ylab('Database Coverage of local species')+
  geom_jitter()+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2

dt <- unique(result[,c(3,6:10)])
p3 <- ggplot(dt, aes(Continent, Prop, fill = Continent, color = Continent))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Database Coverage of local species')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3

pdf('../image/2_DatabaseCoverage_primer.pdf', width = 10, height = 6)
p1
dev.off()
pdf('../image/2_DatabaseCoverage_barcode.pdf', width = 6, height = 8)
p2
dev.off()
pdf('../image/2_DatabaseCoverage_continent.pdf', width = 8, height = 8)
p3
dev.off()
##significance
#---normality and homogeneity of variances
require(car)
shapiro_result <- result %>%
  group_by(Level, Resolution) %>%
  summarise(shapiro_test = list(shapiro.test(Prop)), .groups = "drop") %>%
  mutate(p_value = sapply(shapiro_test, function(x) x$p.value))

# Levene's Test between continets
levene_result_1 <- result %>%
  group_by(Level, Resolution) %>%
  do(leveneTest(Prop ~ Continent, data = .)[1,])
# Levene's Test between barcode
levene_result_2 <- result %>%
  group_by(Level, Resolution) %>%
  do(leveneTest(Prop ~ Barcode, data = .)[1,])
res <- data.frame(shapiro_result[,c(1,2,4)], levene_result_1[,c(3:5)], levene_result_2[,c(3:5)])
colnames(res) <- c('Level','Barcode','shapiro_p','levene_df_conti','levene_F_conti','levene_p_conti','levene_df_barcode','levene_F_barcode','levene_p_barcode')
write.table(res, 'Res_BarcodeCoverage_normality_test.csv', row.names = F, sep = ',', quote = F)
#---overall differences
sig1 <- rbind.data.frame(result %>%
                           group_by(Level, Resolution)%>%
                           do(tidy(kruskal.test(Prop ~ Barcode, data = .))) %>%
                           slice(1)%>%
                           mutate(group = 'Barcode'),
                         result %>%
                           group_by(Level, Resolution)%>%
                           do(tidy(kruskal.test(Prop ~ Continent, data = .))) %>%
                           slice(1)%>%
                           mutate(group = 'Continent'))
#---pairwise differences
sig2 <- rbind.data.frame(result %>%
                           mutate(Barcode = as.factor(Barcode)) %>% 
                           group_by(Level, Resolution) %>%
                           do({
                             dunnTest(Prop ~ Barcode, data = ., method = "bonferroni")$res
                           }) %>%
                           bind_rows() %>%
                           mutate(group = 'Barcode'),
                         result %>%
                           mutate(Continent = as.factor(Continent)) %>% 
                           group_by(Level, Resolution) %>%
                           do({
                             dunnTest(Prop ~ Continent, data = ., method = "bonferroni")$res
                           }) %>%
                           bind_rows() %>%
                           mutate(group = 'Continent'))

write.table(sig1, 'Res_BarcodeCoverage_sig_kruskall.csv', row.names = F, sep = ',', quote = F)
write.table(sig2, 'Res_BarcodeCoverage_sig_pairwise_dunn.csv', row.names = F, sep = ',', quote = F)

##add significant result to the plot
library(multcompView)
sig_res <- result %>%
  mutate(Barcode = as.factor(Barcode)) %>% 
  group_by(Level, Resolution) %>%
  do({
    dunn_result <- dunnTest(Prop ~ Continent, data = ., method = "bonferroni")$res
    p <- dunn_result$P.unadj
    names(p) <- gsub(' ','',dunn_result$Comparison)
    ls <- multcompLetters2(Prop ~ Continent, p, data = ., reversed = T, compare="<", threshold = 0.05)
    data.frame(Continent = names(ls$Letters),Letters = ls$Letters)
  }) %>%
  bind_rows()

Prop_mean <- result %>%
  mutate(Barcode = as.factor(Barcode)) %>% 
  group_by(Level, Resolution, Continent) %>%
  summarise(mean = mean(Prop))
  
dt <- unique(result[,c(3,6:10)])
dt2 <- dt %>%
  left_join(sig_res, by = c("Level", "Resolution", "Continent"))
p3_2 <- ggplot(dt2, aes(Continent, Prop, fill = Continent, color = Continent)) +
  facet_grid(Level ~ Resolution) +
  geom_boxplot(fill = 'Transparent') +
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2) +
  geom_jitter() +
  ylim(0,1.1)+
  stat_summary(aes(label = Letters), fun.y = "max", geom = "text", 
               position = position_dodge(width = 0.75), size = 4, vjust = -0.5, color = 'black') +
  ylab('Database Coverage of Local Species') +
  scale_color_manual(values = useCors) +
  scale_fill_manual(values = useCors) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3_2
pdf('../image/2_DatabaseCoverage_continent2.pdf', width = 8, height = 8)
p3_2
dev.off()
#############################################################
####global assignment verus local species list
#############################################################
info <- read.csv('00PrimerList.csv', row.names = 1)
head(info)

data_info <- read.csv('00DataList.csv')
data_info <- unique(data_info[,c(1,4)])
rownames(data_info) <- data_info$DataID

spe_lists <- readRDS('Dataset_Fish_Recored_C.RDS')
names(spe_lists)

ASV_tables <- readRDS('2_Format_SeperateDataset_ASVs_Stan.RDS')
names(ASV_tables)

data_ids <- names(ASV_tables)

res <- NULL
for (i in 1:length(data_ids)) {
  seqID <- gsub('.*_','',data_ids[i])
  dataID <- gsub('_.*','',data_ids[i])
  
  cat('---------------------\n')
  cat(paste0('----processing ', seqID, ': ', dataID,'...\n'))
  cat('---------------------\n')
  # get fish ASV_ids
  taxa_tb <- as.data.frame(ASV_tables[[data_ids[[i]]]]$Taxa)

  fish_taxa_tb <- taxa_tb[taxa_tb$class=='Actinopteri',]
  fish_taxa_gens <- sapply(strsplit(fish_taxa_tb$species,' '), '[', 1)
  fish_taxa_spes <- paste(fish_taxa_gens, sapply(strsplit(fish_taxa_tb$species,' '), '[', 2))
  fish_taxa_gens[fish_taxa_tb$genus=='Unclassified'] <-
    fish_taxa_spes[fish_taxa_tb$genus=='Unclassified'] <- 'Unclassified'
  remove(taxa_tb)
  # get asv table
  fish_asv_tb <- ASV_tables[[data_ids[[i]]]]$ASV[rownames(fish_taxa_tb),]
  # get record
  b_spes <- spe_lists[[dataID]]$B
  c_spes <- spe_lists[[dataID]]$C
  e_spes <- spe_lists[[dataID]]$E
  
  b_gens <- gsub(' .*','',b_spes)
  c_gens <- gsub(' .*','',c_spes)
  e_gens <- gsub(' .*','',e_spes)
  
  ## 2.1 check how namy original assignment have record
  add <- data.frame(
    RawSeqID = seqID,
    DataID = dataID,
    Continent = data_info[dataID,2],
    Primer = info[seqID,1],
    Barcode = info[seqID,7],
    Level = rep(c('Basin','Country','Ecoregion','New','Poor'), each = 2),
    Resolution = rep(c('Species','Genus'),5),
    Total_ASV = rep(length(fish_taxa_spes),10),
    Identified_ASV = c(sum(fish_taxa_spes %in% b_spes),
                       sum(fish_taxa_gens %in% b_gens),
                       sum(fish_taxa_spes %in% c_spes),
                       sum(fish_taxa_gens %in% c_gens),
                       sum(fish_taxa_spes %in% e_spes),
                       sum(fish_taxa_gens %in% e_gens),
                       sum((fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_spes %in% e_spes))),
                       sum((fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_gens %in% e_gens))),
                       sum(fish_taxa_tb$genus=='Unclassified'),
                       sum(fish_taxa_tb$genus=='Unclassified')
                       ),
    Total_Taxa = rep(c(length(unique(fish_taxa_spes)),length(unique(fish_taxa_gens))),5),
    Identified_Taxa = c(sum(unique(fish_taxa_spes) %in% b_spes),
                       sum(unique(fish_taxa_gens) %in% b_gens),
                       sum(unique(fish_taxa_spes) %in% c_spes),
                       sum(unique(fish_taxa_gens) %in% c_gens),
                       sum(unique(fish_taxa_spes) %in% e_spes),
                       sum(unique(fish_taxa_gens) %in% e_gens),
                       sum((unique(fish_taxa_spes)!='Unclassified')&(!(unique(fish_taxa_spes) %in% e_spes))),
                       sum((unique(fish_taxa_spes)!='Unclassified')&(!(unique(fish_taxa_gens) %in% e_gens))),
                       sum(unique(fish_taxa_spes)=='Unclassified'),
                       sum(unique(fish_taxa_gens)=='Unclassified')
                       ),
    Total_Seq = sum(fish_asv_tb),
    Median_Seq = median(colSums(fish_asv_tb)),
    Mean_Seq = mean(colSums(fish_asv_tb)),
    
    Identified_reads = c(sum(fish_asv_tb[fish_taxa_spes %in% b_spes,]),
                         sum(fish_asv_tb[fish_taxa_gens %in% b_gens,]),
                         sum(fish_asv_tb[fish_taxa_spes %in% c_spes,]),
                         sum(fish_asv_tb[fish_taxa_gens %in% c_gens,]),
                         sum(fish_asv_tb[fish_taxa_spes %in% e_spes,]),
                         sum(fish_asv_tb[fish_taxa_gens %in% e_gens,]),
                         sum(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_spes %in% e_spes)),]),
                         sum(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_gens %in% e_gens)),]),
                         sum(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',]),
                         sum(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',])
    ),
    Median_reads_Prop = c(median(colSums(fish_asv_tb[fish_taxa_spes %in% b_spes,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_gens %in% b_gens,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_spes %in% c_spes,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_gens %in% c_gens,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_spes %in% e_spes,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_gens %in% e_gens,])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_spes %in% e_spes)),])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_gens %in% e_gens)),])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',])/colSums(fish_asv_tb)),
                          median(colSums(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',])/colSums(fish_asv_tb))
    ),
    Mean_reads_Prop = c(mean(colSums(fish_asv_tb[fish_taxa_spes %in% b_spes,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_gens %in% b_gens,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_spes %in% c_spes,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_gens %in% c_gens,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_spes %in% e_spes,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_gens %in% e_gens,])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_spes %in% e_spes)),])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[(fish_taxa_tb$genus!='Unclassified')&(!(fish_taxa_gens %in% e_gens)),])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',])/colSums(fish_asv_tb)),
                        mean(colSums(fish_asv_tb[fish_taxa_tb$genus=='Unclassified',])/colSums(fish_asv_tb))
    )
      
    )
  res <- rbind.data.frame(res, add)
  
}
res$Identified_ASV_prop <- res$Identified_ASV/res$Total_ASV
res$Identified_Taxa_prop <- res$Identified_Taxa/res$Total_Taxa
res$Identified_Reads_prop <- res$Identified_reads/res$Total_Seq

write.table(res,'Res_FishRecord_global_assignment.csv', row.names = F, sep = ',', quote = F)
##plotting ASVs
res <- read.csv('Res_FishRecord_global_assignment.csv')
head(res)
p1 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Primer, Identified_ASV_prop, fill = Primer, color = Primer))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of ASVs assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1

p2 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Barcode, Identified_ASV_prop, fill = Barcode, color = Barcode))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of ASVs assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2

p3 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Continent, Identified_ASV_prop, fill = Continent, color = Continent))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of ASVs assigned to local taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3

##plotting Taxa
head(res)
p4 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Primer, Identified_Taxa_prop, fill = Primer, color = Primer))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Species assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p4

p5 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Barcode, Identified_Taxa_prop, fill = Barcode, color = Barcode))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Species assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p5

p6 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Continent, Identified_Taxa_prop, fill = Continent, color = Continent))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Species assigned to local taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p6
##plotting read
head(res)
p7 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Primer, Identified_Reads_prop, fill = Primer, color = Primer))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Reads assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p7

p8 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Barcode, Identified_Reads_prop, fill = Barcode, color = Barcode))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Reads assigned to local taxa')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p8

p9 <- ggplot(subset(res, Level %in% c('Basin','Country','Ecoregion')), aes(Continent, Identified_Reads_prop, fill = Continent, color = Continent))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Reads assigned to local taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p9
##---
pdf('../image/2_GlobalAssign_local_primer.pdf', width = 10, height = 6)
p1;p4;p7
dev.off()
pdf('../image/2_GlobalAssign_local_barcode.pdf', width = 6, height = 8)
p2;p5;p8
dev.off()
pdf('../image/2_GlobalAssign_local_continent.pdf', width = 8, height = 8)
p3;p6;p9
dev.off()
##plotting exotic
library(dplyr)
dt <- res %>%
  filter((Level == 'Ecoregion')&(Resolution == 'Genus'))%>%
  mutate(exotic_ASV = 1-Identified_ASV_prop,
         exotic_Taxa = 1-Identified_Taxa_prop,
         exotic_Reads = 1-Identified_Reads_prop)

dt <- res %>%
  filter((Level == 'New')&(Resolution == 'Genus'))%>%
  mutate(exotic_ASV = Identified_ASV_prop,
         exotic_Taxa = Identified_Taxa_prop,
         exotic_Reads = Identified_Reads_prop)
p7 <- ggplot(dt, aes(Continent, exotic_ASV, fill = Continent, color = Continent))+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.2, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% ASVs assigned to exotic taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p7
p8 <- ggplot(dt, aes(Continent, exotic_Taxa, fill = Continent, color = Continent))+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.2, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% Species assigned to exotic taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p8

p9 <- ggplot(dt, aes(Continent, exotic_Reads, fill = Continent, color = Continent))+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.2, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% Reads assigned to exotic taxa')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p9
pdf('../image/2_GlobalAssign_exotic_continent.pdf', width = 12, height = 4)
ggarrange(p7,p8,p9,ncol = 3, align = 'hv')
dev.off()

##significance
res <- subset(res, Level %in% c("Basin","Country","Ecoregion"))
#---normality and homogeneity of variances
require(car)
shapiro_result <- res %>%
  group_by(Level, Resolution) %>%
  summarise(shapiro_test = list(shapiro.test(Identified_Taxa_prop)), .groups = "drop") %>%
  mutate(p_value = sapply(shapiro_test, function(x) x$p.value))

# Levene's Test between continets
levene_result_1 <- res %>%
  group_by(Level, Resolution) %>%
  do(leveneTest(Identified_Taxa_prop ~ Continent, data = .)[1,])
# Levene's Test between barcode
levene_result_2 <- res %>%
  group_by(Level, Resolution) %>%
  do(leveneTest(Identified_Taxa_prop ~ Barcode, data = .)[1,])
normal_res <- data.frame(shapiro_result[,c(1,2,4)], levene_result_1[,c(3:5)], levene_result_2[,c(3:5)])
colnames(normal_res) <- c('Level','Barcode','shapiro_p','levene_df_conti','levene_F_conti','levene_p_conti','levene_df_barcode','levene_F_barcode','levene_p_barcode')
write.table(normal_res, 'Res_LocalSpeciesCoverage_GlobalAssign_normality_test.csv', row.names = F, sep = ',', quote = F)
#---overall differences
sig1 <- rbind.data.frame(res %>%
                           group_by(Level, Resolution)%>%
                           do(tidy(kruskal.test(Identified_Taxa_prop ~ Barcode, data = .))) %>%
                           slice(1)%>%
                           mutate(group = 'Barcode'),
                         res %>%
                           group_by(Level, Resolution)%>%
                           do(tidy(kruskal.test(Identified_Taxa_prop ~ Continent, data = .))) %>%
                           slice(1)%>%
                           mutate(group = 'Continent'))
#---pairwise differences
sig2 <- rbind.data.frame(res %>%
                           mutate(Barcode = as.factor(Barcode)) %>% 
                           group_by(Level, Resolution) %>%
                           do({
                             dunnTest(Identified_Taxa_prop ~ Barcode, data = ., method = "bonferroni")$res
                           }) %>%
                           bind_rows() %>%
                           mutate(group = 'Barcode'),
                         res %>%
                           mutate(Continent = as.factor(Continent)) %>% 
                           group_by(Level, Resolution) %>%
                           do({
                             dunnTest(Identified_Taxa_prop ~ Continent, data = ., method = "bonferroni")$res
                           }) %>%
                           bind_rows() %>%
                           mutate(group = 'Continent'))

write.table(sig1, 'Res_LocalSpeciesCoverage_GlobalAssign_sig_kruskall.csv', row.names = F, sep = ',', quote = F)
write.table(sig2, 'Res_LocalSpeciesCoverage_GlobalAssign_sig_pairwise_dunn.csv', row.names = F, sep = ',', quote = F)

##add significant result to the plot
library(multcompView)
sig_res <- res %>%
  mutate(Barcode = as.factor(Barcode)) %>% 
  group_by(Level, Resolution) %>%
  do({
    dunn_result <- dunnTest(Identified_Taxa_prop ~ Continent, data = ., method = "bonferroni")$res
    p <- dunn_result$P.unadj
    names(p) <- gsub(' ','',dunn_result$Comparison)
    ls <- multcompLetters2(Identified_Taxa_prop ~ Continent, p, data = ., reversed = T, compare="<", threshold = 0.05)
    data.frame(Continent = names(ls$Letters),Letters = ls$Letters)
  }) %>%
  bind_rows()

dt2 <- res %>%
  left_join(sig_res, by = c("Level", "Resolution", "Continent"))
p6_2 <-  ggplot(dt2, aes(Continent, Identified_Taxa_prop, fill = Continent, color = Continent))+
  facet_grid(Level~Resolution)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of Species assigned to local taxa')+
  ylim(0,1.1)+
  stat_summary(aes(label = Letters), fun.y = "max", geom = "text", 
               position = position_dodge(width = 0.75), size = 4, vjust = -0.5, color = 'black') +
  scale_color_manual(values = useCors) +
  scale_fill_manual(values = useCors) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p6_2

pdf('../image/2_GlobalAssign_local_continent2.pdf', width = 8, height = 8)
p6_2
dev.off()
