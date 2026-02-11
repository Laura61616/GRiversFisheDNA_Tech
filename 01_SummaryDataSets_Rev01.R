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
####A01: distribution of sites
#############################################################
##calculating the centroids of datasets
Allsites <- read.csv('1_AllSites_Coordinants.csv')
sf_data <- Allsites %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

centroids <- sf_data %>%
  group_split(DataID) %>% 
  lapply(function(group) { 
    st_centroid(st_combine(group$geometry))
  }) %>%
  do.call(rbind, .)

centroids <- data.frame(
  DataID = unique(Allsites$DataID),
  geometry = centroids) %>%
  st_as_sf(crs = 4326)
##number of primer used
info <- read.csv('0_SummaryDataset_C.csv')
tmp <- info %>%
  group_by(DataID,Continent)%>%
  summarise(nPrimer = length(unique(Primer)))
tmp <- as.data.frame(tmp)
rownames(tmp) <- tmp$DataID
centroids$Continent <- tmp$Continent
centroids$nPrimer <- tmp$nPrimer
centroids <- centroids[,c(1,3,4,2)]
st_write(centroids, 'GIS//0_AllDatasets_Coordinants.shp')

sf_data$nPrimer <- tmp[sf_data$DataID,2]
sf_data <- sf_data[,c(1:9,11,10)]
st_write(sf_data, 'GIS/0_AllSites_Coordinants.shp')

##plotting dataset + sites map
useCors_dark <- c('#1C6B8D', '#061A2D', '#7A4B72', '#2F1A3A', '#3B3F91', '#D69C1D')

AllBasins <- st_read('GIS/Basin042017_3119.shp')
UseBasins <- st_read('GIS/00_S04_UseBasins_V03.shp')
pos <- st_read('GIS/0_AllSites_Coordinants.shp')
tmp=unique(as.data.frame(pos)[,c(1,5)])
tmp=tmp[order(tmp$Contnnt),]
tmp = as.data.frame(tmp%>%group_by(Contnnt)%>%mutate(mark = 1:n()))
rownames(tmp) = tmp$DataID
pos$mark <- tmp[pos$DataID,3]
pos2 <- st_read('GIS/0_AllDatasets_Coordinants.shp')
{
  mapworld <- borders("world",
                      fill = 'grey90',
                      colour = 'grey90')
  mp <- ggplot() + 
    mapworld +
    labs(x ='',y="")+
    theme_bw()+
    theme(legend.position = 'none',
          legend.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank())
}
p <- mp+
  geom_sf(aes(geometry = `geometry`),
          data = AllBasins,
          color = 'white',
          fill = 'grey75') +
  geom_sf(aes(geometry = `geometry`, fill = Continent ),
          data = pos2, alpha = 0.8, color = 'grey20', shape = 21, size = 3)+
  geom_sf(aes(geometry = `geometry`, color = Contnnt ),
          data = pos, alpha = 0.8, shape = 3, size = 1, stroke = 0.2, show.legend = F)+
  scale_fill_manual(values = useCors)+
  scale_color_manual(values = useCors_dark)
p
pdf('../image/1_AllDatasets_sites_map.pdf', width = 15, height = 5)
print(p)
dev.off()
#############################################################
####A02: Summary datasets: Data Types
#############################################################
data_list <- read.csv('0_SummaryDataset_C.csv')
data_list <- unique(data_list[,c(1,4,6:8)])
data_list <- data_list%>%
  mutate(type = case_when(RawSeq & SpeMatrix ~ 'RawSeq & SpeTable',
                   RawSeq & Ori_SR & !SpeMatrix ~ 'RawSeq & Richness',
                   RawSeq & !(Ori_SR | SpeMatrix) ~ 'RawSeq only',
                   SpeMatrix & !RawSeq ~ 'SpeTable only',
                   Ori_SR & !(RawSeq | SpeMatrix) ~ 'Richness only'
                   ))
sum_dt_type <- data_list%>%
  group_by(type, Continent)%>%
  summarise(n = n())
sum_dt_type <- rbind.data.frame(sum_dt_type,
                                data.frame('type' = c(rep('RawSeq & Richness',5),rep('RawSeq only',3),rep('Richness only',5)),
                                           'Continent' = c('Africa','Asia','North_America','South_America','Oceania',
                                                           'Africa','Oceania','Europe',
                                                           'Africa','Asia','North_America','Oceania','Oceania'),
                                           'n' = 0))
sum_dt_type$type <- factor(sum_dt_type$type
                           , levels = c('RawSeq & SpeTable', 'RawSeq & Richness',
                                        'RawSeq only', 'SpeTable only', 'Richness only'))
p1 <- ggplot(sum_dt_type, 
            aes(x = type, y = n, fill = Continent))+
  # geom_col(width = 0.7) +
  geom_col(width = 0.7, position = 'dodge2') +
  geom_text(aes(label = n), 
            position = position_dodge2(width = 0.7), 
            vjust = 0.1,hjust = 0.5,size=3) +
  scale_fill_manual(values = useCors) +
  scale_color_manual(values = useCors) +
  labs(y ='Number of datasets', x = '')+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p1
pdf('../image/1_DataTypes_summary.pdf',
    width = 6, height = 5)
p1
dev.off()
#############################################################
####A02: summary barcode and refer_db used
#############################################################
data <- read.csv('0_SummaryDataset_C.csv')
data <- unique(data[,c(1,4,11,12,18)])
data
##summary barcode
dt1 <- data %>%
  group_by(BarcodeRegion, Continent) %>%
  summarise(n = n())
dt1
dt1 <- rbind.data.frame(dt1,
                        data.frame('BarcodeRegion' = c('16S','16S','16S','COI','COI','COI','CytB','CytB','CytB','CytB'),
                                   'Continent' = c('Africa','Asia','South_America',
                                                   'Africa','Oceania','South_America',
                                                   'Africa','Asia','South_America','Oceania'),
                                   n = rep(0,10)))
p2 <- ggplot(dt1, 
             aes(x = BarcodeRegion, y = n, fill = Continent))+
  geom_col(width = 0.7, position = 'dodge2') +
  geom_text(aes(label = n), 
            position = position_dodge2(width = 0.7), 
            vjust = 0.1,hjust = 0.5,size=3) +
  scale_fill_manual(values = useCors) +
  scale_color_manual(values = useCors) +
  labs(y ='Number of datasets', x = '')+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p2
## summary ReferDB
dt2 <- data %>%
  group_by(Continent, ReferDB) %>%
  summarise(n = length(unique(DataID)))%>%
  group_by(Continent)%>%
  mutate(prop = n/sum(n))
dt2 <- rbind.data.frame(dt2,dt2[1,])
dt2[nrow(dt2),'ReferDB'] <- 'Custom_db';dt2[nrow(dt2),3:4] <- 0
dt2$mark <- paste(dt2$Continent, dt2$ReferDB)
dt2[dt2$ReferDB=='Global_db','mark'] <- 'Global_db'
dt2 <- dt2[c(1,12,2:11),]

dt2$ymax = dt2$prop
dt2[dt2$ReferDB=='Global_db','ymax'] <- 1
dt2$ymin = 0
dt2[dt2$ReferDB=='Global_db','ymin'] <- dt2[dt2$ReferDB=='Custom_db','ymax']

useCor2 <- c(useCors, 'grey')
marks <- c("Africa Custom_db", "Asia Custom_db","Europe Custom_db", "North_America Custom_db", "Oceania Custom_db","South_America Custom_db"
           ,"Global_db")

p3 <- ggplot(dt2, 
             aes(x = ReferDB, y = n, fill = Continent))+
  geom_col(width = 0.7, position = 'dodge2') +
  geom_text(aes(label = n), 
            position = position_dodge2(width = 0.7), 
            vjust = 0.5, hjust = 0.5, size = 3) +
  scale_fill_manual(values = useCors) +
  scale_color_manual(values = useCors) +
  labs(y ='Number of datasets', x = '')+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p3

pdf('../image/1_Barcode_and_db_summary.pdf', width = 6, height = 5)
ggarrange(p2,p3,ncol = 2, align = 'hv',widths = c(2/3,1/3))
dev.off()

pdf('../image/1_Type_and_Barcode_and_db_summary.pdf', width = 12, height = 5)
ggarrange(p1,p2,p3,ncol = 3, align = 'hv',widths = c(4/7, 3/7, 2/7))
dev.off()
#############################################################
####A03: summary number of sites
#############################################################
data <- read.csv('1_AllSites_Coordinants.csv')
head(data)
dt <- data %>%
  group_by(DataID, Continent) %>%
  summarise(Tsite = length(N_SiteID))
dt <- dt[order(dt$Tsite),]
dt$ID <- as.factor(1:nrow(dt))
p1 <- ggplot(dt, 
             aes(x = ID, y = Tsite, fill = Continent))+
  geom_col(width = 0.5) +
  scale_fill_manual(values = useCors) +
  labs(x ='Unique datasets', y = "Number of sites")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p1

dt2 <- dt
dt2 <- dt2 %>%
  mutate('site_cate' = 
           case_when(Tsite > 20 ~ "Exceed 20",
                     (Tsite <= 20)&(Tsite > 10) ~ "11 to 20",
                     (Tsite <= 10)&(Tsite > 3) ~ "4 to 10",
                     Tsite == 3 ~ "3",
                     Tsite == 2 ~ "2",
                     Tsite == 1 ~ "1" ))
dt2$site_cate <- factor(
  dt2$site_cate,
  levels = unique(dt2$site_cate, fromLast = T))
dt3 <- dt2 %>%
  group_by(Continent, site_cate) %>%
  summarise(Tstudy = n())

p2 <- ggplot(dt3, aes(x = site_cate, y = Tstudy, fill = Continent))+
  coord_flip()+
  geom_col(position = 'stack') +
  scale_fill_manual(values = useCors) +
  labs(x ='', y = "Number of datasets")+
  theme_bw()+
  theme(legend.position = 'right',
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p2
## save
pdf('../image/1_SiteNumber_data_1.pdf', width = 6, height = 3)
p1
dev.off()

pdf('../image/1_SiteNumber_data_2.pdf', width = 6, height = 3)
p2
dev.off()
#############################################################
####A03: summary number of species
#############################################################
data <- read.csv('0_SummaryDataset_ItemID.csv')
head(data)
dt <- unique(data[,c('DataID','Continent','Tspe')])
dt <- dt[order(dt$Tspe),]
dt <- na.omit(dt)
dt$ID <- as.factor(1:nrow(dt))
p1 <- ggplot(dt, 
             aes(x = ID, y = Tspe, fill = Continent))+
  geom_col(width = 0.5) +
  scale_fill_manual(values = useCors) +
  labs(x ='Unique datasets', y = "Number of species")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.75),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p1

dt2 <- dt
dt2 <- dt2 %>%
  mutate('site_cate' = 
           case_when(Tspe > 20 ~ "Exceed 20",
                     (Tspe <= 20)&(Tspe > 10) ~ "11 to 20",
                     (Tspe <= 10)&(Tspe > 5) ~ "6 to 10",
                     (Tspe <= 5)&(Tspe > 1) ~ "2 to 5"))
dt2$site_cate <- factor(
  dt2$site_cate,
  levels = unique(dt2$site_cate, fromLast = T))
dt3 <- dt2 %>%
  group_by(Continent, site_cate) %>%
  summarise(Tstudy = n())

p2 <- ggplot(dt3, aes(x = site_cate, y = Tstudy, fill = Continent))+
  coord_flip()+
  geom_col(position = 'stack') +
  scale_fill_manual(values = useCors) +
  labs(x ='', y = "Number of datasets")+
  theme_bw()+
  theme(legend.position = 'right',
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
p2
## save
pdf('../image/1_SpeciesNumber_data_1.pdf', width = 6, height = 3)
p1
dev.off()

pdf('../image/1_SpeciesNumber_data_2.pdf', width = 6, height = 3)
p2
dev.off()

#############################################################
####Rev01: representatives of the reanalyzed subset
#############################################################
data_list <- read.csv('0_SummaryDataset_ItemID.csv')
head(data_list)
info <- unique(data_list[,c("DataID","Continent","Year","Primer","BarcodeRegion","BarcodeLength","SamVolume","SamReps","TotVolume","ReferDB","RawSeqID","ReAssigned")])
head(info)

info$type <- case_when(info$ReAssigned %in% c("", 'no') ~ 'unreanalyzable', TRUE ~ 'reanalyzable')##93 items
summary(as.factor(info$type))#63 reanalyzable and 30 un-reanalyzable

##year of sampling
data1 <- info
data1$Year <- as.numeric(data1$Year)
data1 <- na.omit(data1)

wilcox_test_1 <- wilcox.test(Year ~ type, data = data1)
p1 <- ggplot(data1, aes(type, Year, fill = type, color = type))+
  geom_boxplot(alpha = 0.5)+
  geom_jitter(size = 2)+
  ylab('Year of sampling')+xlab("")+
  scale_color_aaas()+
  scale_fill_aaas()+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))+ 
  geom_signif(comparisons = list(c("reanalyzable", "unreanalyzable")), 
              map_signif_level = TRUE, 
              annotations = paste("Wilcoxon p = ", sprintf('%.3f', wilcox_test_1$p.value), sep = ''))
p1

##sampling volume
wilcox_test_2 <- wilcox.test(log10(TotVolume) ~ type, data = info)

p2 <- ggplot(info, aes(type, log10(TotVolume), fill = type, color = type))+
  geom_boxplot(alpha = 0.5)+
  geom_jitter(size = 2)+
  scale_y_continuous(breaks = c(-1,0,1,2), labels = c(0.1,1,10,100))+
  ylab('Sampling volume per site (L)')+xlab("")+
  scale_color_aaas()+
  scale_fill_aaas()+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))+ 
  geom_signif(comparisons = list(c("reanalyzable", "unreanalyzable")), 
              map_signif_level = TRUE, 
              annotations = paste("Wilcoxon p = ", sprintf('%.3f', wilcox_test_2$p.value)))
p2

##continent
info2 <- unique(info[c('DataID', 'type', 'Continent')])
data2 <- info2 %>%
  group_by(type, Continent)%>%
  summarise(n = n())%>%
  ungroup()

tab <- table(info$type,info$Continent)
fisher_1 <- fisher.test(tab)

p3 <- ggplot(as.data.frame(data2), aes(x = type, y = n))+
  geom_col(aes(fill = Continent, color = Continent), position = "stack") +
  scale_fill_manual(values = useCors) +
  scale_color_manual(values = useCors) +
  ylab('Number of datasets')+xlab("")+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))+
  geom_signif(comparisons = list(c("reanalyzable", "unreanalyzable")), 
              map_signif_level = TRUE, 
              annotations = paste("Fisher's exact test p = ", sprintf('%.3f', fisher_1$p.value), sep = ''))
p3

##BarcodeRegion
info3 <- unique(info[c('DataID', 'type', 'BarcodeRegion')])
data3 <- info3 %>%
  group_by(type, BarcodeRegion)%>%
  summarise(n = n())%>%
  ungroup()

tab <- table(info$type, info$BarcodeRegion)
fisher_2 <- fisher.test(tab)

p4 <- ggplot(as.data.frame(data3), aes(x = type, y = n))+
  geom_col(aes(fill = BarcodeRegion, color = BarcodeRegion), position = "stack") +
  ylab('Number of datasets')+xlab("")+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45,hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'grey'))+
  geom_signif(comparisons = list(c("reanalyzable", "unreanalyzable")), 
              map_signif_level = TRUE, 
              annotations = paste("Fisher's exact test p = ", sprintf('%.3f', fisher_2$p.value), sep = ''))
p4
pdf('../image/1_AllDatasets_representative_new.pdf', width = 8, height = 13)
ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2, align = 'hv', labels = letters[1:4])
dev.off()
