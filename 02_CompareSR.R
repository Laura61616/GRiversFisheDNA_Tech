# READNE: summary dataset, 20.11.2024
# Created by Yan Zhang (Nanjing University, China)
# ALL RIGHTS RESERVED
wd <- '../data/'
setwd(wd)

fun_path <- '../code/00_functions.R' ## CHANGE ME
source(fun_path)
useCors <- c('#2EC4B6', '#0A2845', '#E4A9CE', '#452E68','#4F5DCE', '#F4D35E')
#############################################################
####compare corrected SR with original ones--all datasets
#############################################################
data <- read.csv('2_Format_SeperateDataset_SR.csv')
head(data)
data <- data[(!is.na(data$SR))&(!is.na(data$nASVs_G)),]
data <- as.data.frame(data %>%
                        group_by(Continent) %>%
                        mutate(mark = as.numeric(factor(DataID))) %>%
                        ungroup())
data$SR <- sqrt(data$SR);
data$nSpes_B <- sqrt(data$nSpes_B); data$nSpes_C <- sqrt(data$nSpes_C);
data$nSpes_E <- sqrt(data$nSpes_E); data$nSpes_G <- sqrt(data$nSpes_G)

# basin lev
b <- lmer(nSpes_B ~ SR + (1|DataID) + (1|Primer), data = data)
cond_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "conditional", nboot = 10,
                  max_level = 1)
marg_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "marginal", nboot = 10,
                  max_level = 1)
fit <- summary(b)

lab <- paste0('Fixed+Random = '
            , round(unlist(cond_r2$R2[2,2])*100,2)
            , '%\n'
            , lme_fit_lab(fit, 'SR', unlist(marg_r2$R2[2,2])))
x = min(data$SR, na.rm = T) +
  (max(data$SR, na.rm = T)-min(data$SR, na.rm = T))/4
y = max(data$nSpes_B, na.rm = T) -
  (max(data$nSpes_B, na.rm = T)-min(data$nSpes_B,na.rm = T))/10
p1 <- ggplot(data, aes(x = SR, y = nSpes_B)) +
  geom_point(aes(color = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15), labels = c(0,9,36,81,144,225))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_B),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  annotate('text',
           x = x, y = y,label = lab,
           size = 3) +
  ylab('Reclassified SR (Basin)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1
# country lev
b <- lmer(nSpes_C ~ SR + (1|DataID) + (1|Primer), data = data)
cond_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "conditional", nboot = 10,
                  max_level = 1)
marg_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "marginal", nboot = 10,
                  max_level = 1)
fit <- summary(b)

lab <- paste0('Fixed+Random = '
              , round(unlist(cond_r2$R2[2,2])*100,2)
              , '%\n'
              , lme_fit_lab(fit, 'SR', unlist(marg_r2$R2[2,2])))
x = min(data$SR, na.rm = T) +
  (max(data$SR, na.rm = T)-min(data$SR, na.rm = T))/4
y = max(data$nSpes_C, na.rm = T) -
  (max(data$nSpes_C, na.rm = T)-min(data$nSpes_C,na.rm = T))/10
p2 <- ggplot(data, aes(x = SR, y = nSpes_C)) +
  geom_point(aes(color = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15), labels = c(0,9,36,81,144,225))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_C),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  annotate('text',
           x = x, y = y,label = lab,
           size = 3) +
  ylab('Reclassified SR (Country)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2
# ecoregion lev
b <- lmer(nSpes_E ~ SR + (1|DataID) + (1|Primer), data = data)
cond_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "conditional", nboot = 10,
                  max_level = 1)
marg_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "marginal", nboot = 10,
                  max_level = 1)
fit <- summary(b)

lab <- paste0('Fixed+Random = '
              , round(unlist(cond_r2$R2[2,2])*100,2)
              , '%\n'
              , lme_fit_lab(fit, 'SR', unlist(marg_r2$R2[2,2])))
x = min(data$SR, na.rm = T) +
  (max(data$SR, na.rm = T)-min(data$SR, na.rm = T))/4
y = max(data$nSpes_E, na.rm = T) -
  (max(data$nSpes_E, na.rm = T)-min(data$nSpes_E,na.rm = T))/10
p3 <- ggplot(data, aes(x = SR, y = nSpes_E)) +
  geom_point(aes(color = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15), labels = c(0,9,36,81,144,225))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_E),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  annotate('text',
           x = x, y = y,label = lab,
           size = 3) +
  ylab('Reclassified SR (Ecoregion)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3
# global lev
b <- lmer(nSpes_G ~ SR + (1|DataID) + (1|Primer), data = data)
cond_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "conditional", nboot = 10,
                  max_level = 1)
marg_r2 <- partR2(b, partvars = 'SR', data = data,
                  R2_type = "marginal", nboot = 10,
                  max_level = 1)
fit <- summary(b)

lab <- paste0('Fixed+Random = '
              , round(unlist(cond_r2$R2[2,2])*100,2)
              , '%\n'
              , lme_fit_lab(fit, 'SR', unlist(marg_r2$R2[2,2])))
x = min(data$SR, na.rm = T) +
  (max(data$SR, na.rm = T)-min(data$SR, na.rm = T))/4
y = max(data$nSpes_G, na.rm = T) -
  (max(data$nSpes_G, na.rm = T)-min(data$nSpes_G,na.rm = T))/10
p4 <- ggplot(data, aes(x = SR, y = nSpes_G)) +
  geom_point(aes(color = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15), labels = c(0,9,36,81,144,225))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_G),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  annotate('text',
           x = x, y = y,label = lab,
           size = 3) +
  ylab('Reclassified SR (Global)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p4
pdf('../image/2_ReClassify_CorrectedTaxa_AllSR_compare_sqrt.pdf', width = 11,height = 8)
ggarrange(p4,p3,p2,p1,ncol = 2, nrow = 2, align = 'hv', common.legend = T, legend = 'right')
dev.off()
#############################################################
####compare corrected SR with original ones--all datasets
#############################################################
info <- read.csv('0_SummaryDataset_ItemID.csv')
data <- read.csv('2_Format_SeperateDataset_SR.csv')
data <- merge.data.frame(data, info[,c('ItemID','ReferDB')], by = 'ItemID')
head(data)
data <- data[(!is.na(data$SR))&(!is.na(data$nASVs_G)),]
data <- as.data.frame(data %>%
                        group_by(Continent) %>%
                        mutate(mark = as.numeric(factor(DataID))) %>%
                        ungroup())
data$SR <- sqrt(data$SR);
data$nSpes_B <- sqrt(data$nSpes_B); 
data$nSpes_C <- sqrt(data$nSpes_C); 
data$nSpes_E <- sqrt(data$nSpes_E); 
data$nSpes_G <- sqrt(data$nSpes_G)
# basin lev
p1 <- ggplot(data, aes(x = SR, y = nSpes_B)) +
  geom_point(aes(color = Continent, fill = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_B, linetype = ReferDB),
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1, inherit.aes = F)+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_B),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  ylab('Reclassified SR (Basin)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p1
# country lev
p2 <- ggplot(data, aes(x = SR, y = nSpes_C)) +
  geom_point(aes(color = Continent, fill = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_C, linetype = ReferDB),
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1, inherit.aes = F)+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_C),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  ylab('Reclassified SR (Country)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p2
# ecoregion lev
p3 <- ggplot(data, aes(x = SR, y = nSpes_E)) +
  geom_point(aes(color = Continent, fill = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_E, linetype = ReferDB),
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1, inherit.aes = F)+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_E),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  ylab('Reclassified SR (Ecoregion)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p3
# global lev
p4 <- ggplot(data, aes(x = SR, y = nSpes_G)) +
  geom_point(aes(color = Continent, fill = Continent, shape = as.character(mark)))+
  scale_color_manual(values = useCors)+
  scale_fill_manual(values = useCors)+
  scale_shape_manual(values = c(18,9,17:10))+
  scale_x_continuous(breaks = c(0,3,6,9,12), labels = c(0,9,36,81,144))+
  scale_y_continuous(breaks = c(0,3,6,9,12,15), labels = c(0,9,36,81,144,225))+
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey')+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_G, linetype = ReferDB),
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1, inherit.aes = F)+
  geom_smooth(mapping = aes(x = SR,
                            y = nSpes_G),
              color = 'grey60',
              fill = 'grey60',
              method = 'glm',
              formula = y ~ x,
              level = 0.95,
              lwd = 1,
              alpha = 0.1)+
  ylab('Reclassified SR (Global)') + xlab('Original SR')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))
p4

pdf('../image/2_ReClassify_CorrectedTaxa_AllSR_compare_sqrt.pdf', width = 11,height = 8)
ggarrange(p4,p3,p2,p1,ncol = 2, nrow = 2, align = 'hv', common.legend = T, legend = 'right')
dev.off()
#############################################################
####compare corrected SR with original ones--all datasets--LME
#############################################################
info <- read.csv('0_SummaryDataset_ItemID.csv')
data <- read.csv('2_Format_SeperateDataset_SR.csv')
data <- merge.data.frame(data, info[,c('ItemID','ReferDB')], by = 'ItemID')
head(data)
data <- data[(!is.na(data$SR))&(!is.na(data$nASVs_G)),]
data$SR <- sqrt(data$SR);
data$nSpes_B <- sqrt(data$nSpes_B); data$nSpes_C <- sqrt(data$nSpes_C);
data$nSpes_E <- sqrt(data$nSpes_E); data$nSpes_G <- sqrt(data$nSpes_G)
# 
dt <- melt(data[,c(12,6,11,22,18:21)], id.vars = c('ID','Continent','ReferDB','SR'))

res <- NULL
for (i in unique(dt$variable)) {
  dt_use <- subset(dt, variable == i)
  b <- lmer(value ~ SR + (1|ID), data = dt_use)
  cond_r2 <- partR2(b, partvars = 'SR', data = dt_use,
                    R2_type = "conditional", nboot = 10,
                    max_level = 1)
  marg_r2 <- partR2(b, partvars = 'SR', data = dt_use,
                    R2_type = "marginal", nboot = 10,
                    max_level = 1)
  fit <- summary(b)
  
  add1 <- data.frame(
    scale = i,
    type = 'All Datasets',
    t(fit$coefficients[2,]),
    cond_r2 = cond_r2$R2[2,2],
    marg_r2 = marg_r2$R2[2,2])
  #global_db
  b <- lmer(value ~ SR + (1|ID), data = subset(dt_use, ReferDB == 'Global_db'))
  cond_r2 <- partR2(b, partvars = 'SR', data = subset(dt_use, ReferDB == 'Global_db'),
                    R2_type = "conditional", nboot = 10,
                    max_level = 1)
  marg_r2 <- partR2(b, partvars = 'SR', data = subset(dt_use, ReferDB == 'Global_db'),
                    R2_type = "marginal", nboot = 10,
                    max_level = 1)
  fit <- summary(b)
  
  add2 <- data.frame(
    scale = i,
    type = 'Global_db',
    t(fit$coefficients[2,]),
    cond_r2 = cond_r2$R2[2,2],
    marg_r2 = marg_r2$R2[2,2])
  #custom_db
  b <- lmer(value ~ SR + (1|ID), data = subset(dt_use, ReferDB == 'Custom_db'))
  cond_r2 <- partR2(b, partvars = 'SR', data = subset(dt_use, ReferDB == 'Custom_db'),
                    R2_type = "conditional", nboot = 10,
                    max_level = 1)
  marg_r2 <- partR2(b, partvars = 'SR', data = subset(dt_use, ReferDB == 'Custom_db'),
                    R2_type = "marginal", nboot = 10,
                    max_level = 1)
  fit <- summary(b)
  
  add3 <- data.frame(
    scale = i,
    type = 'Custom_db',
    t(fit$coefficients[2,]),
    cond_r2 = cond_r2$R2[2,2],
    marg_r2 = marg_r2$R2[2,2])
  ##combine
  res <- rbind.data.frame(res,add1, add2, add3)
}
colnames(res) <- c('Scale','Type','Estimate', 'STD', 'df','t','p','cond_R2','marg_R2')
write.table(res, 'Res_Compare_SpeciesRichness_LME.csv', row.names = F, sep = ',', quote = F)
