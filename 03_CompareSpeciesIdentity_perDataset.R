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
  new_spe_tb_stan <- sweep(new_spe_tb, 2, colSums(new_spe_tb), '/')
  ## old_spe_tb
  old_spe_tb <- old_assign[[dataID]]
  ids <- intersect(colnames(new_spe_tb), colnames(old_spe_tb))
  if (length(ids) < 2) {
    print('check!')
    break
  }
  old_spes <- rownames(old_spe_tb)[rowSums(old_spe_tb[,ids]) > 0]
  ## species record
  add <- data.frame(
    ID = use_data_ids[i],
    DataID = info[seqID,'DataID'],
    primer = info[seqID,'N_PN'],
    barcode = info[seqID,'Target.region'],
    continent = info[seqID,'Continent'],
    type = 'basin',
    tot.n.seqs = sum(new_spe_tb),
    tot.mean.seqs = mean(colSums(new_spe_tb), na.rm = T),
    tot.median.seqs = median(colSums(new_spe_tb), na.rm = T),
    
    Species = new_spes,
    overlapped = new_spes%in%old_spes,
    
    prop.n.seqs = rowSums(new_spe_tb[new_spes,])/sum(new_spe_tb),
    prop.mean.seqs = rowMeans(new_spe_tb_stan[new_spes,], na.rm = T),
    prop.median.seqs = apply(new_spe_tb_stan[new_spes,], 1, function(x) median(x, na.rm = T))
  )
  rownames(add) <- NULL
  res <- rbind.data.frame(res, add)
}

res$type2 <- case_when(res$overlapped ~ 'Overlapped',
                       TRUE ~ 'Unique to reanalyzed')
write.table(res,'2_Format_SeperateDataset_Spetb_compare_perDataset.csv', row.names = F, sep = ',', quote = F)

#############################################################
####boxplot
#############################################################
tmp <- read.csv('2_Format_SeperateDataset_Spetb_compare.csv')

res <- read.csv('2_Format_SeperateDataset_Spetb_compare_perDataset.csv')
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

# data$ReferDB <- factor(data$ReferDB, levels = c('Global_db','Custom_db'))
data$type2 <- factor(data$type2, levels = c('Overlapped','Unique to reanalyzed'))
data$ID <- factor(data$ID, levels = tmp[order(rowSums(tmp[,7:8]), decreasing = T),'ID'])

p1 <- ggplot(subset(data, ReferDB == 'Global_db'), aes(x = ID, y = sqrt(prop.n.seqs), fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = type2), method = "wilcox.test",label = "p.signif",
                     label.y = 1.1, bracket.size = 0.5, tip.length = 0.01,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns") ) ) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9), labels = c(0,0.09,0.36,0.81))+
  ylab("Proportion of Reads") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p1

p2 <- ggplot(subset(data, ReferDB == 'Custom_db'), aes(x = ID, y = sqrt(prop.n.seqs), fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = type2), method = "wilcox.test",label = "p.signif",
                     label.y = 1.1, bracket.size = 0.5, tip.length = 0.01,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns") ) ) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9), labels = c(0,0.09,0.36,0.81))+
  ylab("Proportion of Reads") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p2

p3 <- ggplot(subset(data, ReferDB == 'Global_db'), aes(x = ID, y = sqrt(prop.mean.seqs), fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = type2), method = "wilcox.test",label = "p.signif",
                     label.y = 1.1, bracket.size = 0.5, tip.length = 0.01,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns") ) ) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9), labels = c(0,0.09,0.36,0.81))+
  ylab("Mean relative abundance") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p3

p4 <- ggplot(subset(data, ReferDB == 'Custom_db'), aes(x = ID, y = sqrt(prop.mean.seqs), fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  stat_compare_means(aes(group = type2), method = "wilcox.test",label = "p.signif",
                     label.y = 1.1, bracket.size = 0.5, tip.length = 0.01,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns") ) ) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9), labels = c(0,0.09,0.36,0.81))+
  ylab("Mean relative abundance") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p4

p5 <- ggplot(subset(data, ReferDB == 'Global_db'), aes(x = ID, y = prop.n.seqs*tot.n.seqs, fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  ylab("Number of sequences") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p5

p6 <- ggplot(subset(data, ReferDB == 'Custom_db'), aes(x = ID, y = sqrt(prop.n.seqs*tot.n.seqs), fill = type2)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("Overlapped" = "#1f77b4", "Unique to reanalyzed" = "#ff7f0e")) +
  # scale_y_continuous(breaks = c(0,0.3,0.6,0.9), labels = c(0,0.09,0.36,0.81))+
  ylab("Number of sequences") +
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
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )
p6

pdf('../image/3_Compare_Species_Identity_boxplot_perDataset.pdf', width = 15, height = 12)
ggarrange(p1,p2,p3,p4,ncol = 2,nrow=2,align = 'hv', widths = c(2/3,1/3), common.legend = T, legend = 'right')
dev.off()
##
library(glmmTMB)
sink('3_compare_Species_Identity_boxplot_perDataset_glmm.log')
mod1 <- glmmTMB(sqrt(prop.n.seqs) ~ type2 + Year + ReferDB + (1|ID), data)
mod2 <- glmmTMB(sqrt(prop.mean.seqs) ~ type2 + Year + ReferDB + (1|ID), data)
summary(mod1)
summary(mod2)

MuMIn::r.squaredGLMM(mod1)
MuMIn::r.squaredGLMM(mod2)

sink()
