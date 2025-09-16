# READNE: summary dataset, 20.11.2024
# Created by Yan Zhang (Nanjing University, China)
# ALL RIGHTS RESERVED
wd <- '../data/'
setwd(wd)

fun_path <- '../code/00_functions.R' ## CHANGE ME
source(fun_path)
useCors <- c('#2EC4B6', '#0A2845', '#E4A9CE', '#452E68','#4F5DCE', '#F4D35E')
#############################################################
####compare corrected Taxa with original Taxa
#############################################################
info <- read.csv('00PrimerList.csv')
tmp <- read.csv('00DataList.csv')
tmp <- unique(tmp[,c(1,4)])
info <- merge.data.frame(info, tmp, by = 'DataID')
rownames(info) <- info$RawSeqID

old_assign <- readRDS('1_Original_SpeTableStan.RDS')

####basin
new_assign <- readRDS('2_Format_SeperateDataset_Spetb_LocalAssign_Basin.RDS')
records <- readRDS('Dataset_Fish_Recored_C.RDS')

names(new_assign)
names(old_assign)
names(records)

use_data_ids <- names(new_assign)
use_data_ids <- use_data_ids[gsub('_.*','',use_data_ids)%in%names(old_assign)]
res <- NULL
for (i in 1:length(use_data_ids)) {
  seqID <- gsub('.*_','',use_data_ids[i])
  dataID <- gsub('_.*','',use_data_ids[i])
  
  cat('---------------------\n')
  cat(paste0('----processing ', dataID, ': ', seqID,'...\n'))
  cat('---------------------\n')
  ## new_spe_tb
  new_spe_tb <- new_assign[[use_data_ids[i]]]
  new_spes <- rownames(new_spe_tb)
  ## old_spe_tb
  old_spe_tb <- old_assign[[dataID]]
  ids <- intersect(colnames(new_spe_tb), colnames(old_spe_tb))
  if (length(ids) < 2) {
    print('check!')
    break
  }
  old_spes <- rownames(old_spe_tb)[rowSums(old_spe_tb[,ids]) > 0]
  ## species record
  species_list <- records[[dataID]]$B
  
  add <- data.frame(
    ID = use_data_ids[i],
    DataID = info[seqID,'DataID'],
    primer = info[seqID,'N_PN'],
    barcode = info[seqID,'Target.region'],
    continent = info[seqID,'Continent'],
    type = 'basin',
    n.spes.common = length(intersect(new_spes, old_spes)),
    n.spes.new = length(setdiff(new_spes, old_spes)),
    n.spes.old = length(setdiff(old_spes, new_spes)),
    n.spes.old.local = length(intersect(setdiff(old_spes, new_spes), species_list)),
    n.spes.old.notlocal = length(setdiff(setdiff(old_spes, new_spes), species_list)),
    
    n.seqs = sum(new_spe_tb),
    mean.seqs = mean(colSums(new_spe_tb), na.rm = T),
    median.seqs = median(colSums(new_spe_tb), na.rm = T),
    prop.reads.common = sum(new_spe_tb[intersect(new_spes, old_spes),])/sum(new_spe_tb)
  )
  add_spes <- intersect(setdiff(old_spes, new_spes), species_list)
  global_taxonomy <- read.table(paste0('db_taxonomy/', add$primer,'_taxonomy.txt'), sep = '\t',header = F)
  global_fish_taxonomy <- global_taxonomy[grep('Actinopteri', global_taxonomy$V2),]
  global_fish_spes <- unique(gsub('.*;', '',global_fish_taxonomy$V2))
  add$n.spes.old.local.in.db <- sum(add_spes%in%global_fish_spes)
  res <- rbind.data.frame(res, add)
}
write.table(res,'2_Format_SeperateDataset_Spetb_compare.csv',
            row.names = F, sep = ',', quote = F)
#############################################################
####barplot-2
#############################################################
res <- read.csv('2_Format_SeperateDataset_Spetb_compare.csv')
head(res)
info <- read.csv('0_SummaryDataset_C.csv')
info_cleaned <- info %>%
  group_by(DataID) %>%
  summarise(
    ReferDB = first(unique(ReferDB)),
    Year = max(Year),
    .groups = "drop" )
data <- merge.data.frame(res, info_cleaned, by = 'DataID')
head(data)

dt <- melt(data[,c(1:8,10,11,17,18)], id.vars = c('DataID','ID','primer','barcode','continent','type','ReferDB','Year'))
head(dt)

dt$ID <- factor(dt$ID,levels = data[order(rowSums(data[,7:8]), decreasing = T),'ID'])
dt <- dt[order(dt$ID),]
dt$mark <- dt$continent; 
dt[dt$variable=='n.spes.old.notlocal','mark'] <- 'aa'; 
dt[dt$variable=='n.spes.old.local','mark'] <- 'ab'; 
dt[dt$variable=='n.spes.new','mark'] <- 'ac'

p1 <- ggplot(subset(dt, ReferDB=='Global_db'), aes(x = ID, y = value, fill = mark)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('grey','black','grey40',useCors),
                    labels = c(
                      "aa" = "Unique other spe in orginal",
                      "ab" = "Unique local spe in orginal",
                      "ac" = 'Unique in re-classified',
                      'Africa' = 'Africa',
                      'Asia' = 'Asia',
                      'Europe' = 'Europe',
                      'North_America' = 'North_America',
                      'Oceania' = 'Oceania',
                      'South_America' = 'South_America'),
                    breaks = c('aa','ab','ac','Africa','Asia','Europe','North_America','Oceania','South_America')) + 
  # scale_y_continuous(limits = c(0,30),breaks = c(0,5,10,15,20,25), labels = c(0,25,100,225,400,625))+
  ylim(0,250)+
  ylab("Number of species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p1

p2 <- ggplot(subset(dt, ReferDB=='Custom_db'), aes(x = ID, y = value, fill = mark)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('grey','black','grey40',useCors),
                    labels = c(
                      "aa" = "Unique other spe in orginal",
                      "ab" = "Unique local spe in orginal",
                      "ac" = 'Unique in re-classified',
                      'Africa' = 'Africa',
                      'Asia' = 'Asia',
                      'Europe' = 'Europe',
                      'North_America' = 'North_America',
                      'Oceania' = 'Oceania',
                      'South_America' = 'South_America'),
                    breaks = c('aa','ab','ac','Africa','Asia','Europe','North_America','Oceania','South_America')) + 
  ylab("Number of species") + 
  ylim(0,250)+
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p2

g1 <- ggarrange(p1,p2,ncol = 2, align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
g1

#---
dt5 <- melt(data[,c(1:6,10,16,17,18)], id.vars = c('DataID','ID','primer','barcode','continent','type','ReferDB','Year'))
head(dt5)
dt5$ID <- factor(dt5$ID,levels = data[order(rowSums(data[,7:8]), decreasing = T),'ID'])
p3 <- ggplot(subset(dt5, ReferDB == 'Global_db'), aes(x = ID, y = sqrt(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = sqrt(10), linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = c('black','red'))+
  scale_y_continuous(limits = c(0,11),breaks = c(0,3,6,9), labels = c(0,9,36,81))+
  ylab("Number of species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p3
p4 <- ggplot(subset(dt5, ReferDB == 'Custom_db'), aes(x = ID, y = sqrt(value), fill = variable)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = sqrt(10), linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = c('black','red'))+
  scale_y_continuous(limits = c(0,11),breaks = c(0,3,6,9), labels = c(0,9,36,81))+
  ylab("Number of species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p4
g2 <- ggarrange(p3,p4,ncol = 2, align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
g2
#---
dt6 <- data
dt6$ID <- factor(dt6$ID,levels = data[order(rowSums(data[,7:8]), decreasing = T),'ID'])
p5 <- ggplot(subset(dt6, ReferDB == 'Global_db'), aes(x = ID, y = prop.reads.common, fill = continent)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = useCors)+
  ylim(0,1)+
  ylab("% Reads of overlapped species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p5
p6 <- ggplot(subset(dt6, ReferDB == 'Custom_db'), aes(x = ID, y = prop.reads.common, fill = continent)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = useCors)+
  ylim(0,1)+
  ylab("% Reads of overlapped species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p6
g3 <- ggarrange(p5,p6,ncol = 2, align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
g3
#---
data$prop.spe.new <- data$n.spes.common/rowSums(data[,7:8])
dt6 <- data
dt6$ID <- factor(dt6$ID,levels = data[order(rowSums(data[,7:8]), decreasing = T),'ID'])
p7 <- ggplot(subset(dt6, ReferDB == 'Global_db'), aes(x = ID, y = prop.spe.new, fill = continent)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = useCors)+
  ylim(0,1)+
  ylab("% Number of overlapped species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p7
p8 <- ggplot(subset(dt6, ReferDB == 'Custom_db'), aes(x = ID, y = prop.spe.new, fill = continent)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey90')+
  scale_fill_manual(values = useCors)+
  ylim(0,1)+
  ylab("% Number of overlapped species") + 
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8))
p8

#---
pdf('../image/3_Compare_Species_Identity_barplot.pdf', width = 12, height = 12)
ggarrange(p1,p2,p5,p6,p3,p4,ncol = 2,nrow = 3, align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
ggarrange(p1,p2,p7,p8,p3,p4,ncol = 2,nrow = 3, align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
dev.off()
#############################################################
####boxplot
#############################################################
res <- read.csv('2_Format_SeperateDataset_Spetb_compare.csv')
head(res)
info <- read.csv('0_SummaryDataset_C.csv')
info_cleaned <- info %>%
  group_by(DataID) %>%
  summarise(
    ReferDB = first(unique(ReferDB)),
    Year = max(Year),
    .groups = "drop" )
data <- merge.data.frame(res, info_cleaned, by = 'DataID')
head(data)

data$prop_overlap_re <- data$n.spes.common/rowSums(data[,c(7,8)])
data$prop_overlap_ori <- data$n.spes.common/rowSums(data[,c(7,9)])
data$ReferDB <- factor(data$ReferDB, levels = c('Global_db','Custom_db'))

p1 <- ggplot(data, aes(x = continent, y =prop_overlap_re , fill = continent, color = continent)) +
  facet_wrap(~ReferDB)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% overlapped species\n(Reanalyzed assignments)')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1
p2 <- ggplot(data, aes(x = continent, y = prop_overlap_ori , fill = continent, color = continent)) +
  facet_wrap(~ReferDB)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% overlapped species\n(Original assignments)')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2

p3 <- ggplot(data, aes(x = continent, y = sqrt(n.spes.old.local) , fill = continent, color = continent)) +
  facet_wrap(~ReferDB)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = sqrt(10), linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of extra local species\nin Original assignments')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  scale_y_continuous(breaks = c(0,3,6,9), labels = c(0,9,36,81))+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3

p4 <- ggplot(subset(data, DataID != 'D013'), aes(x = continent, y =n.spes.old.local , fill = continent, color = continent)) +
  facet_wrap(~ReferDB)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('Number of extra local species\nin Original assignments')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p4

p5 <- ggplot(data, aes(x = continent, y = prop.reads.common , fill = continent, color = continent)) +
  facet_wrap(~ReferDB)+
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% Reads of overlapped species\nin Reanalyzed assignments')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p5
pdf('../image/3_Compare_Species_Identity_boxplot_1.pdf', width = 8, height = 4)
p1;p2;p3;p4;p5
dev.off()

p <- ggarrange(p2,p1,p3, nrow = 3, align = 'hv')
pdf('../image/3_Compare_Species_Identity_boxplot_2.pdf', width = 7, height = 12)
p
dev.off()
##
check <- data[(data$n.spes.old.local > 10)|(data$n.spes.old.local>data$n.spes.common),]
head(check)
coverage <- read.csv('Res_FishRecord_primer_db.csv')
head(coverage)
coverage <- subset(coverage, (Level == 'Basin')&(Resolution == 'Species'))
rownames(coverage) <- coverage$RawSeqID
check$db_coverage <- coverage[gsub('.*_','',check$ID),'Prop'] 
head(check)

write.table(check, 'Res_Compare_SpeciesIdentity_extraLocal_10.csv',
            row.names = F, sep = ',', quote = F)
#### not refer_db
p1 <- ggplot(data, aes(x = continent, y =prop_overlap_re , fill = continent, color = continent)) +
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% overlapped species (Re-classified)')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1
p2 <- ggplot(data, aes(x = continent, y =prop_overlap_ori , fill = continent, color = continent)) +
  geom_boxplot(fill = 'Transparent')+
  geom_hline(yintercept = 0.5, linetype = 'dashed', lwd = 0.2)+
  geom_jitter()+
  ylab('% overlapped species (Original)')+xlab('')+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2
#############################################################
####regression verus year, sequencing depth, referDB
#############################################################
res <- read.csv('2_Format_SeperateDataset_Spetb_compare.csv')
head(res)
info <- read.csv('0_SummaryDataset_C.csv')
info_cleaned <- info %>%
  group_by(DataID) %>%
  summarise(
    ReferDB = first(unique(ReferDB)),
    Year = max(Year),
    .groups = "drop" )
data <- merge.data.frame(res, info_cleaned, by = 'DataID')
data$prop.spes.new <- data$n.spes.common/rowSums(data[,7:8])
data$SR_Re <- rowSums(data[,7:8])
head(data)

summary(gam(prop.reads.common ~ log10(median.seqs+1):ReferDB, data = data))
summary(gam(prop.reads.common ~ Year:ReferDB, data = data))
summary(gam(prop.spes.new ~ log10(median.seqs+1):ReferDB, data = data))
summary(gam(prop.spes.new ~ Year:ReferDB, data = data))

mod1 <- gam(prop.reads.common ~ log10(median.seqs+1) + Year + ReferDB + SR_Re, data = data)
mod2 <- gam(prop.spes.new ~ log10(median.seqs+1) + Year + ReferDB + SR_Re, data = data)
sink('2_Format_SeperateDataset_Spetb_compare_gamModel.log')
summary(mod1)
summary(mod2)
sink()

p1 <- ggplot(data, aes(x = log10(median.seqs+1), y = prop.reads.common, size = Year)) +
  geom_point(aes(color = continent))+
  scale_color_manual(values = useCors)+
  scale_size_continuous(breaks = 2013:2023)+
  scale_x_continuous(breaks = 0:5, labels = c(1,10,100,1000,10000,100000))+
  geom_smooth(mapping = aes(x = log10(median.seqs+1), y = prop.reads.common),
              color = 'grey60',fill = 'grey60',method = 'glm',formula = y ~ x,level = 0.95,lwd = 1,alpha = 0.1)+
  xlab('Median of sequencing depth') + ylab('% Reads of overlapped species')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1

p2 <- ggplot(data, aes(x = log10(median.seqs+1), y = prop.spes.new, size = Year)) +
  geom_point(aes(color = continent))+
  scale_color_manual(values = useCors)+
  scale_size_continuous(breaks = 2013:2023)+
  scale_x_continuous(breaks = 0:5, labels = c(1,10,100,1000,10000,100000))+
  geom_smooth(mapping = aes(x = log10(median.seqs+1), y = prop.spes.new),
              color = 'grey60',fill = 'grey60',method = 'glm',formula = y ~ x,level = 0.95,lwd = 1,alpha = 0.1)+
  xlab('Median of sequencing depth') + ylab('% Number of overlapped species')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2

pdf('../image/3_Compare_Species_Identity_proportion_curve.pdf', width = 14, height = 5)
ggarrange(p1,p2,ncol = 2, align = 'hv', common.legend = T, legend = 'right',labels = c('a','b'))
dev.off()
#############################################################
####significance
#############################################################
res <- read.csv('2_Format_SeperateDataset_Spetb_compare.csv')
head(res)
info <- read.csv('0_SummaryDataset_C.csv')
data <- merge.data.frame(res, unique(info[,c('DataID','ReferDB')]), by = 'DataID')
head(data)

data <- subset(data, type=='basin')
data$prop_overlap_re <- data$n.spes.common/rowSums(data[,c(7,8)])
data$prop_overlap_ori <- data$n.spes.common/rowSums(data[,c(7,9)])
##check normality and homogeneity of variance
y_vars <- c("prop_overlap_re", "prop_overlap_ori", "n.spes.old.local")

res <- NULL
for (y in y_vars) {
  require(car)
  shapiro_result <- shapiro.test(data[[y]])$p.value
  
  # Levene's Test between continets
  levene_result_1 <- data %>%
    do(leveneTest(as.formula(paste(y, "~ continent")), data = .)[1,])
  # Levene's Test between global_db & custom_db
  levene_result_2 <- data %>%
    do(leveneTest(as.formula(paste(y, "~ ReferDB")), data = data)[1,])
  
  # Append results to the results data frame
  add <- data.frame(y,shapiro_result, 
             levene_result_1, 
             levene_result_2)
  colnames(add) <- c('y','shapiro_p','levene_df_conti','levene_F_conti','levene_p_conti','levene_df_ReferDB','levene_F_barcode','levene_p_ReferDB')
  res <- rbind(res, add)
}
rownames(res) <- NULL
print(res)
write.table(res, 'Res_Compare_SpeciesIdentity_normality_check.csv', row.names = F, sep = ',',quote = F)
##compare between continents
sig_res_overall <- rbind.data.frame(data %>%
                                      kruskal.test(prop_overlap_re ~ continent, data = .) %>%
                                      tidy()%>%
                                      mutate(y = "prop_overlap_re")%>%
                                      slice(1),
                                    data %>%
                                      kruskal.test(prop_overlap_ori ~ continent, data = .) %>%
                                      tidy()%>%
                                      mutate(y = "prop_overlap_ori")%>%
                                      slice(1),
                                    data %>%
                                      kruskal.test(n.spes.old.local ~ continent, data = .) %>%
                                      tidy()%>%
                                      mutate(y = "n.spes.old.local")%>%
                                      slice(1))

sig_posthoc <- rbind(dunnTest(prop_overlap_re ~ continent, data = data, method = "bonferroni")$res %>%
                       as.data.frame() %>%
                       mutate(y = "prop_overlap_re"),
                     dunnTest(prop_overlap_ori ~ continent, data = data, method = "bonferroni")$res %>%
                       as.data.frame() %>%
                       mutate(y = "prop_overlap_ori"),
                     dunnTest(n.spes.old.local ~ continent, data = data, method = "bonferroni")$res %>%
                       as.data.frame() %>%
                       mutate(y = "n.spes.old.local"))
write.table(sig_res_overall, 'Res_Compare_SpeciesIdentity_sig_conti.csv', row.names = F, sep = ',',quote = F)
write.table(sig_posthoc, 'Res_Compare_SpeciesIdentity_sig_conti_posthoc.csv', row.names = F, sep = ',',quote = F)
##compare between ReferDB
sig_res_overall <- rbind(data %>%
                 pairwise_wilcox_test(prop_overlap_re ~ ReferDB, p.adjust.method = "bonferroni"),
               data %>%
                 pairwise_wilcox_test(prop_overlap_ori ~ ReferDB, p.adjust.method = "bonferroni"),
               data %>%
                 pairwise_wilcox_test(n.spes.old.local ~ ReferDB, p.adjust.method = "bonferroni"))
write.table(sig_res_overall, 'Res_Compare_SpeciesIdentity_sig_ReferDB.csv', row.names = F, sep = ',',quote = F)
##comparison between two groups
sig_3 <- rbind(data %>%
                 group_by(ReferDB) %>%
                 pairwise_wilcox_test(prop_overlap_re ~ continent, p.adjust.method = "bonferroni"),
               data %>%
                 group_by(ReferDB) %>%
                 pairwise_wilcox_test(prop_overlap_ori ~ continent, p.adjust.method = "bonferroni"),
               data %>%
                 group_by(ReferDB) %>%
                 pairwise_wilcox_test(n.spes.old.local ~ continent, p.adjust.method = "bonferroni"))
write.table(sig_3, 'Res_Compare_SpeciesIdentity_sig_Conti_by_ReferDB.csv', row.names = F, sep = ',',quote = F)
sig_3 <- rbind(
  data %>%
    group_by(ReferDB) %>%
    do(tidy(aov(prop_overlap_re ~ continent, data = .))) %>%    
    mutate(y = "prop_overlap_re")%>%
    slice(1,3),
  data %>%
    group_by(ReferDB) %>%
    do(tidy(aov(prop_overlap_ori ~ continent, data = .))) %>%    
    mutate(y = "prop_overlap_ori")%>%
    slice(1,3),
  data %>%
    group_by(ReferDB) %>%
    do(tidy(aov(n.spes.old.local ~ continent, data = .))) %>%    
    mutate(y = "n.spes.old.local")%>%
    slice(1,3))
write.table(sig_3, 'Res_Compare_SpeciesIdentity_sig_Conti_by_ReferDB_anova.csv', row.names = F, sep = ',',quote = F)
#############---------************-----------################
#############################################################
####NMDS---prepare data
#############################################################
old_assign <- readRDS('1_Original_SpeTableStan.RDS')
global_assign <- readRDS('2_Format_SeperateDataset_Spetb_LocalAssign_Global.RDS')
basin_assign <- readRDS('2_Format_SeperateDataset_Spetb_LocalAssign_Basin.RDS')

names(global_assign)
names(old_assign)

use_data_ids <- names(global_assign)
use_data_ids <- use_data_ids[gsub('_.*','',use_data_ids)%in%names(old_assign)]
all_spes_new <- Reduce(union, lapply(global_assign[use_data_ids], rownames))
all_spes_old <- Reduce(union, lapply(old_assign[gsub('_.*','',use_data_ids)], rownames))
all_spes <- union(all_spes_new, all_spes_old)
## 2061 species: 1413 new and 1084 old
res <- as.data.frame(matrix(0,nrow = length(all_spes), ncol = length(use_data_ids)*3,
                            dimnames = list(all_spes, 
                                              c(paste(use_data_ids, 'original', sep = '-'),
                                                paste(use_data_ids, 'basin', sep = '-'),
                                                paste(use_data_ids, 'global', sep = '-')))))
for (i in 1:length(use_data_ids)) {
  seqID <- gsub('.*_','',use_data_ids[i])
  dataID <- gsub('_.*','',use_data_ids[i])
  
  cat('---------------------\n')
  cat(paste0('----processing ', dataID, ': ', seqID,'...\n'))
  cat('---------------------\n')
 
  ## old_spe_tb
  old_spe_tb <- old_assign[[dataID]]
  old_spes <- rownames(old_spe_tb)[rowSums(old_spe_tb) > 0]
  res[old_spes,paste(use_data_ids[i],'original',sep = '-')] <- 1
  ## global_spe_tb
  global_spe_tb <- global_assign[[use_data_ids[i]]]
  global_spes <- rownames(global_spe_tb)
  res[global_spes,paste(use_data_ids[i],'global',sep = '-')] <- 1
  ## basin_spe_tb
  if (use_data_ids[i]%in%names(basin_assign)) {
    basin_spe_tb <- basin_assign[[use_data_ids[i]]]
    basin_spes <- rownames(basin_spe_tb)
    res[basin_spes,paste(use_data_ids[i],'basin',sep = '-')] <- 1
  }
  
}
summary(rowSums(res)>0)
summary(colSums(res)>0)

res <- res[rowSums(res)>0,colSums(res)>0]
summary(rowSums(res))
summary(colSums(res))

write.table(res, 'Res_Format_SeperateDataset_Spetb_toNMDS.csv',
            sep = ',', quote = F)
#############################################################
####NMDS---running
#############################################################
info <- read.csv('00PrimerList.csv')
tmp <- read.csv('00DataList.csv')
tmp <- unique(tmp[,c(1,4)])
info <- merge.data.frame(info, tmp, by = 'DataID')

res <- read.csv('Res_Format_SeperateDataset_Spetb_toNMDS.csv')
colnames(res)
rownames(info)

data_info <- data.frame(
  ID = colnames(res),
  DataID = gsub('_.*','', colnames(res)),
  RawSeqID = gsub('\\..*','',gsub('.*_','', colnames(res))),
  type = gsub('.*\\.','', colnames(res))
)
data_info <- merge.data.frame(data_info, info[,c(2,3,8,14)], by = 'RawSeqID')
rownames(data_info) <- data_info$ID

## PCoA plotting
res.d <- metaMDS(t(res), distance = 'jaccard',trymax = 500, autotransform = FALSE, noshare = TRUE)
print(res.d)

stress <- res.d$stress
tp <- as.data.frame(res.d$points); 
colnames(tp) <- c('x','y')

tp <- data.frame(tp,
                 data_info[rownames(tp),])
p1 <- ggplot(tp, aes(x = x, y = y, shape =  type))+
  geom_point(mapping = aes(color = Continent),
             size = 2.5, 
             show.legend = T, 
             alpha = 0.8)+
  # stat_ellipse(mapping = aes(fill = Continent),
  #              type="norm",
  #              geom="polygon",
  #              alpha=0.2,color=NA, level = 0.95)+
  scale_shape_manual(values = c(15:17))+
  scale_color_manual(values = useCors)+
  xlab('NMDS1')+ylab('NMDS2')+
  annotate('text', x = -3.5, y = 2.5, label = paste0('Stress = ', round(stress,2)), size = 4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 legend.position = 'right')
p1

pdf('../image/3_Compare_Community_nmds.pdf', width = 6, height = 4)
p1
dev.off()

#############################################################
####PERMANOVA
#############################################################
info <- read.csv('00PrimerList.csv')
tmp <- read.csv('00DataList.csv')
tmp <- unique(tmp[,c(1,4)])
info <- merge.data.frame(info, tmp, by = 'DataID')

res <- read.csv('2_UseData/01_TechPaper_V01/Res_Format_SeperateDataset_Spetb_toNMDS.csv')
colnames(res)
rownames(info)

data_info <- data.frame(
  ID = colnames(res),
  DataID = gsub('_.*','', colnames(res)),
  RawSeqID = gsub('\\..*','',gsub('.*_','', colnames(res))),
  type = gsub('.*\\.','', colnames(res))
)
data_info <- merge.data.frame(data_info, info[,c(2,3,8,14)], by = 'RawSeqID')
rownames(data_info) <- data_info$ID

summary(colnames(res)%in%rownames(data_info))
data_info <- data_info[colnames(res),]

dist_matrix <- vegdist(t(res), method = "jaccard")
adonis_res_1 <- adonis2(dist_matrix ~ Continent, data = data_info, permutations = 999)
adonis_res_2 <- adonis2(dist_matrix ~ type, data = data_info, permutations = 999)
anosim_res_1 <- anosim(dist_matrix, grouping = data_info$Continent, permutations = 999)
anosim_res_2 <- anosim(dist_matrix, grouping = data_info$type, permutations = 999)

sink('2_UseData/01_TechPaper_V01/Res_Compare_SpeciesIdentity_nmds_sig.txt')
print('PERMANOVA for continents...')
print(adonis_res_1)
print('-------------------------------')
print('PERMANOVA for Type...')
print(adonis_res_2)
print('-------------------------------')
print('adonis2 for continents...')
print(anosim_res_1)
print('-------------------------------')
print('adonis2 for Type...')
print(anosim_res_2)
print('-------------------------------')
sink()
