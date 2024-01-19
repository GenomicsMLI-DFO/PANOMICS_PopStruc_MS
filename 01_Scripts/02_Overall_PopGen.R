# Info --------------------------------------------------------------------

# Population genetics analysis of Pandalus borealis
#  at global spatial scale
#  1513 individuals
#  14331 snps with LD r2 < 0.5 over all individuals
#
# Audrey Bourret
# 2022-04-05
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(ggrepel)
library(parallel)


library("ggVennDiagram")
library(eulerr)
library(RColorBrewer)

library(multcompView)

library(adegenet)
library(hierfstat)
library(dartR)


library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")
library(robust)

`%nin%` = Negate(`%in%`)

library(QuickPop)


library(raster)
library(corrplot)
library(vegan)

#Libraries that we will need
library(codep)
library(adespatial)
library(adegraphics)
library(ape)
library(car)


library(raster)
library(rgdal)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)

library(terra)
library(tidyterra)

library(rEEMSplots)


#detach("package:ggVennDiagram", character.only = T)
library(sp)
library(rworldmap)
library(rworldxtra)


# Data --------------------------------------------------------------------

pop.data <- read_csv("00_Data/00_FIleInfos/Projet_Infos.csv") 

region.df <- tibble(Gen_ZONE = pop.data$Gen_ZONE %>% unique()) %>% 
  mutate(RegionAssesment = ifelse(Gen_ZONE %in% c("SFA-13", "SFA-14", "SFA-15", "SFA-16"), "Scotian",
                           ifelse(Gen_ZONE %in% c("SFA-8", "SFA-9", "SFA-10", "SFA-12"), "GSL",
                           ifelse(Gen_ZONE %in% c("SFA-4","SFA-5", "SFA-6", "SFA-7"), "NFL",
                           ifelse(Gen_ZONE %in% c("SFA-0", "WAZ", "EAZ"), "Northern",
                           ifelse(Gen_ZONE %in% "NAFO-3M", "NAFO-3M",NA))))))

region.df

pop.data <- pop.data %>% left_join(region.df) %>% 
              mutate(Gen_ZONE = ifelse(Gen_ZONE == "SFA-0", "EAZ", Gen_ZONE),
                                Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                                       "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                                       "SFA-9", "SFA-10", "SFA-12",
                                                                       "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                                       "NAFO-3M")),
                                Gen_ZONE_FG = factor(Gen_ZONE_FG, levels = c("SFA-0-1",
                                                                             paste(rep(levels(Gen_ZONE), each = 10), 1:10 , sep = "-"))),
                                Gen_ZONE_FG = factor(Gen_ZONE_FG),
                                Transcripto = ifelse(Gen_ZONE_FG %in% c("SFA-6-3", "SFA-8-4", "SFA-12-2", "SFA-12-3", "SFA-15-1"), "Yes", "No")
                                )


pop.too.small <-c("SFA-0-1", "WAZ-3", "SFA-7-7")

# Genetic Data ------------------------------------------------------------

load("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.adegenet.Rdata")

gl.final
gi.final

# Removing NAFO
gl.sfa <- gl.final[pop(gl.final) %nin%  c("NAFO-3M-1", "NAFO-3M-2", "NAFO-3M-3", "NAFO-3M")]
gi.sfa <- gi.final[pop(gi.final) %nin% c("NAFO-3M-1", "NAFO-3M-2", "NAFO-3M-3", "NAFO-3M")]
pop(gl.sfa) %>% table()

# Create subset with pop predefined 
gl.final.Gen_ZONE_FG <- gl.final
pop(gl.final.Gen_ZONE_FG) <- data.frame(ID_GQ = indNames(gl.final)) %>% left_join(pop.data) %>% pull(Gen_ZONE_FG)

gl.final.Gen_ZONE <- gl.final
pop(gl.final.Gen_ZONE) <- data.frame(ID_GQ = indNames(gl.final)) %>% left_join(pop.data) %>% pull(Gen_ZONE)

gl.final.Region <- gl.final
pop(gl.final.Region) <- data.frame(ID_GQ = indNames(gl.final)) %>% left_join(pop.data) %>% pull(RegionAssesment)

gl.sfa.Gen_ZONE_FG <- gl.sfa
pop(gl.sfa.Gen_ZONE_FG) <- data.frame(ID_GQ = indNames(gl.sfa)) %>% left_join(pop.data) %>% pull(Gen_ZONE_FG)

gl.sfa.Gen_ZONE <- gl.sfa
pop(gl.sfa.Gen_ZONE) <- data.frame(ID_GQ = indNames(gl.sfa)) %>% left_join(pop.data) %>% pull(Gen_ZONE)

gl.sfa.Region <- gl.sfa
pop(gl.sfa.Region) <- data.frame(ID_GQ = indNames(gl.sfa)) %>% left_join(pop.data) %>% pull(RegionAssesment)

# GI version

gi.final.Gen_ZONE_FG <- gi.final
pop(gi.final.Gen_ZONE_FG) <- data.frame(ID_GQ = indNames(gi.final)) %>% left_join(pop.data) %>% pull(Gen_ZONE_FG)

gi.final.Gen_ZONE <- gi.final
pop(gi.final.Gen_ZONE) <- data.frame(ID_GQ = indNames(gi.final)) %>% left_join(pop.data) %>% pull(Gen_ZONE)

gi.final.Region <- gi.final
pop(gi.final.Region) <- data.frame(ID_GQ = indNames(gi.final)) %>% left_join(pop.data) %>% pull(RegionAssesment)

# Env data from raster directly -------------------------------------------

## Loading the climatic rasters
ras_all <- stack(list.files("00_Data/99_SIG/", pattern = ".tif", full.names = T))
names(ras_all) <- names(ras_all) %>% str_remove("mon_")
names(ras_all) 

ras_current <- ras_all[[str_subset(names(ras_all), "2075", negate = T)]]
names(ras_current)

ras_future <- ras_all[[str_subset(names(ras_all), "2075", negate = F)]]
names(ras_future) <- names(ras_future) %>% str_remove("_2075")

remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}

ras_current <- remove.NAs.stack(ras_current)
ras_future  <- remove.NAs.stack(ras_future)

Env.data <- data.frame(Gen_ZONE_FG = unique(pop(gi.final.Gen_ZONE_FG))) %>% 
            left_join(pop.data %>% dplyr::select(Gen_ZONE_FG, Lat, Long) %>% distinct(.keep_all = T))

final.env.names <- data.frame(ID = c("Sbtm.ann","SSS.ann","SST.Larval","Tbtm.Summer","Tbtm.Winter"),
                              new.ID = c("BS-a", "SSS-a", "SST-l","BT-s", "BT-w"))

Env.ras <- data.frame(raster::extract(ras_current, Env.data[, c("Long", "Lat")]))
Env.r  <- cor(Env.ras, method = "pearson")

dimnames(Env.r)[[1]] <- data.frame(ID = dimnames(Env.r)[[1]]) %>% 
  left_join(final.env.names) %>% pull(new.ID)

dimnames(Env.r)[[2]] <- data.frame(ID = dimnames(Env.r)[[2]]) %>% 
  left_join(final.env.names) %>% pull(new.ID)

gg.r.env <- ggcorrplot::ggcorrplot(Env.r,   hc.order = TRUE, type = "full", lab = TRUE) +
            theme_bw(base_size = 8) + theme(axis.title = element_blank(),
                                            strip.background = element_rect(fill="white"),
                                             legend.position = "bottom") +
  facet_wrap(~"Spearman's correlations")

gg.r.env 


row.names(Env.ras) <- Env.data$Gen_ZONE_FG

Env.data <- bind_cols(Env.data, Env.ras)

# Do a PCA here
# PCA 
pca.env <- vegan::rda(Env.ras, scale = T)

res.env <- as.data.frame(scores(pca.env , display="sites", scaling=1)) %>% 
  mutate(Gen_ZONE_FG = dimnames(scores(pca.env, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct) %>% 
  mutate(Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                "SFA-9", "SFA-10", "SFA-12",
                                                "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                "NAFO-3M")),
         RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M")))

arrow.env.df <- data.frame(name = scores(pca.env, display="species", choices=1, scaling=1) %>% row.names() ,
                           x0 = 0,
                           y0 = 0,
                           xmax = scores(pca.env, display="species", choices=1, scaling=1) %>% as.vector(),
                           ymax = scores(pca.env, display="species", choices=2, scaling=1) %>% as.vector())  
arrow.factor = 1/4

perc <- round(100*(summary(pca.env)$cont$importance[2, 1:2]), 1)

gg.pca.env <-  res.env %>% 
  ggplot(aes(x =PC1, y = PC2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(fill = RegionAssesment), pch = 21, cex = 2, alpha = 0.5)+
  scale_fill_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta")
  )+ 
  
  geom_segment(data = arrow.env.df, aes(x = x0, y = y0, xend = xmax * arrow.factor, yend= ymax * arrow.factor),
               arrow = arrow(), size = 0.5)+
  geom_label_repel(data = arrow.env.df %>% left_join(final.env.names, by = c("name" = "ID")) , aes(label = new.ID, x = xmax *arrow.factor, y = ymax * arrow.factor), size = 2, max.overlaps = 20
  ) +
  #  geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-|NAFO-")), size = 3, max.overlaps = 20
  #  ) +
  labs(x = paste0("PC1 (", perc[1], "%)"), 
       y = paste0("PC2 (", perc[2], "%)")) + 
  #annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS.northern)$adj.r.squared %>% as.numeric() %>%  round(digits = 3)),
  #          x = -0.2, y = 0.25, vjust = "inward", hjust = "inward") +
  facet_wrap(~"Principale component analysis") +
  theme_bw(base_size = 8) + theme(legend.position = "none",
                                  strip.background = element_rect(fill="white"))

gg.pca.env

gg.env <- ggpubr::ggarrange(gg.r.env,
                            gg.pca.env,
                            labels = LETTERS,
                            widths = c(2,2))
gg.env

# Then check future variations
Env.ras.future <- data.frame(raster::extract(ras_future, Env.data[, c("Long", "Lat")]))
Env.ras.future$Gen_ZONE_FG <- Env.data$Gen_ZONE_FG
Env.ras$Gen_ZONE_FG <- Env.data$Gen_ZONE_FG


# Create a df.anomaly

env.change.df <- Env.ras %>% pivot_longer(-Gen_ZONE_FG) %>% 
  left_join(Env.ras.future %>% pivot_longer(-Gen_ZONE_FG), 
            by = c("Gen_ZONE_FG", "name"), 
            suffix = c(".current", ".future"))  %>% 
  left_join(final.env.names, by = c("name" = "ID")) %>% 
  mutate(pred.change = value.future - value.current) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE_FG, Gen_ZONE, RegionAssesment) %>% distinct(.keep_all = T)) %>% 
  mutate(Region = ifelse(Gen_ZONE %in% c("SFA-4", "SFA-5", "SFA-5", "SFA-7"), "NFL", 
                         ifelse(Gen_ZONE == "NAFO-3M", "3M",  RegionAssesment)),
         Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
        new.ID = factor(new.ID, levels = c("BT-w", "BT-s", "SST-l", "BS-a", "SSS-a"))) 

env.change.gg2 <- env.change.df %>% ggplot(aes(x = Region, y = value.current)) +
  geom_boxplot(col = "darkgray") +
  geom_jitter(aes(col = Region), alpha = 0.5, height = 0)+
  facet_wrap(~new.ID, scale = "free", nrow = 1) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))  + 
  scale_x_discrete(labels = c("ARC", "NFL", "GSL", "SS", "FLC")) +
  labs(x= "", y = "Current value (°C or ppm)") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", strip.background = element_rect(fill="white"))
env.change.gg2 

 env.change.gg.big2 <- ggarrange(gg.env,
                                env.change.gg2,
                                nrow = 2, ncol = 1, align = "hv", labels = c("", LETTERS[3:10]),
                                heights = c(2,1)
 )
 
 env.change.gg.big2
 
ggsave(plot = env.change.gg.big2,
       filename = "02_Results/01_Overall_PopGen/Raw_EnvData_20231017.png",
      height = 6, width = 8, units = "in", bg = "white")


# Identify outliers #1 - Bayescan ----------------------------------------------------------------

# Outlier can be loaded without reanalyses

?gl2bayescan

# Data conversion
gl.final.Gen_ZONE_FG@other$loc.metrics.flags$monomorphs <- FALSE
gl.final.Gen_ZONE_FG@ploidy <- rep(as.integer(2), nInd(gl.final.Gen_ZONE_FG))

gl2bayescan(gl.final.Gen_ZONE_FG, outfile = "00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.final.Gen_ZONE_FG.txt", outpath = ".")


gl.final.Gen_ZONE@other$loc.metrics.flags$monomorphs <- FALSE
gl.final.Gen_ZONE@ploidy <- rep(as.integer(2), nInd(gl.final.Gen_ZONE))

gl2bayescan(gl.final.Gen_ZONE, outfile = "00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.final.Gen_ZONE.txt", outpath = ".")

gl.sfa.Gen_ZONE_FG@other$loc.metrics.flags$monomorphs <- FALSE
gl.sfa.Gen_ZONE_FG@ploidy <- rep(as.integer(2), nInd(gl.sfa.Gen_ZONE_FG))

gl2bayescan(gl.sfa.Gen_ZONE_FG, outfile = "00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.sfa.Gen_ZONE_FG.txt", outpath = ".")


gl.sfa.Gen_ZONE@other$loc.metrics.flags$monomorphs <- FALSE
gl.sfa.Gen_ZONE@ploidy <- rep(as.integer(2), nInd(gl.sfa.Gen_ZONE))

gl2bayescan(gl.sfa.Gen_ZONE, outfile = "00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.sfa.Gen_ZONE.txt", outpath = ".")

# Run bayescan

cmd <- paste("00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.final.Gen_ZONE_FG.txt",
             "-snp",
             "-od", "./02_Results/01_Overall_PopGen/01_Bayescan/", #Output directory file
             "-threads", 20
)

cmd

#A <- system2("bayescan", cmd, stdout = T, stderr = T)
A

cmd <- paste("00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.final.Gen_ZONE.txt",
             "-snp",
             "-od", "./02_Results/01_Overall_PopGen/01_Bayescan/", #Output directory file
             "-threads", 20
)

cmd

#A <- system2("bayescan", cmd, stdout = T, stderr = T)
A

cmd <- paste("00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.sfa.Gen_ZONE_FG.txt",
             "-snp",
             "-od", "./02_Results/01_Overall_PopGen/01_Bayescan/", #Output directory file
             "-threads", 20
)

cmd

#A <- system2("bayescan", cmd, stdout = T, stderr = T)
A


cmd <- paste("00_Data/06_Filtering/06g_UniqueFinal/bayescan_gl.sfa.Gen_ZONE.txt",
             "-snp",
             "-od", "./02_Results/01_Overall_PopGen/01_Bayescan/", #Output directory file
             "-threads", 20
)

cmd

#A <- system2("bayescan", cmd, stdout = T, stderr = T)
A



source("~/Documents/Programs/BayeScan/R_functions/plot_R.r")

plot_bayescan("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.final.Gen_ZONE_FG_fst.txt", FDR = 0.05/15000, add_text= F)
plot_bayescan("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.final.Gen_ZONE_fst.txt", FDR = 0.05/15000, add_text = F)
plot_bayescan("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.sfa.Gen_ZONE_FG_fst.txt", FDR = 0.05)
plot_bayescan("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.sfa.Gen_ZONE_fst.txt", FDR = 0.05)

bayes.final.Gen_ZONE_FG <- read.table("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.final.Gen_ZONE_FG_fst.txt")
bayes.final.Gen_ZONE <- read.table("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.final.Gen_ZONE_fst.txt")
bayes.sfa.Gen_ZONE_FG <- read.table("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.sfa.Gen_ZONE_FG_fst.txt")
bayes.sfa.Gen_ZONE <- read.table("02_Results/01_Overall_PopGen/01_Bayescan/bayescan_gl.sfa.Gen_ZONE_fst.txt")


add.bayes.names <- function(fst, gl){
  loc <- data.frame(LOC = adegenet::locNames(gl)) %>% 
    mutate(LOC.num = row.names(.))
  
  fst %>% mutate(LOC.num = row.names(.)) %>% left_join(loc)
  
}


hist(bayes.final.Gen_ZONE$fst)

# The estimated alpha coefficient indicating the strength and direction of selection. 
# A positive value of alpha suggests diversifying selection, whereas negative values 
# suggest balancing or purifying selection.
#

bayes.final.Gen_ZONE_FG <- add.bayes.names(bayes.final.Gen_ZONE_FG, gl.final.Gen_ZONE_FG)
bayes.final.Gen_ZONE    <- add.bayes.names(bayes.final.Gen_ZONE, gl.final.Gen_ZONE)
bayes.sfa.Gen_ZONE_FG   <- add.bayes.names(bayes.sfa.Gen_ZONE_FG, gl.sfa.Gen_ZONE_FG)
bayes.sfa.Gen_ZONE      <- add.bayes.names(bayes.sfa.Gen_ZONE, gl.sfa.Gen_ZONE)

bayes.final.Gen_ZONE_FG %>% filter(LOC == "1464007:59:+")# %>% pull(LOC)

bayes.final.Gen_ZONE %>% filter(qval < 0.05) %>%   pull(alpha) %>% hist() 


bayes.final.Gen_ZONE %>% mutate(Cat = ifelse(qval <= 3.5e-06, "Top outliers",
                                      ifelse(qval < 0.05, "All outliers", "Neutral"))) %>% 
  
  ggplot() +
  geom_point(aes(x=LOC.num, y=-log10(qval), col = Cat), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
 # geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 1) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


test <- read_table("00_Data/06_Filtering/06a_r75_MAF05/populations.hapstats.tsv")
head(test)

bayes.outliers <- list(
  bayes.final.site = bayes.final.Gen_ZONE_FG  %>% filter(qval < 0.05) %>% pull(LOC) , 
  bayes.final.sfa = bayes.final.Gen_ZONE  %>% filter(qval < 0.05) %>% pull(LOC) , 
  bayes.sfa.site = bayes.sfa.Gen_ZONE_FG  %>% filter(qval < 0.05) %>% pull(LOC) , 
  bayes.sfa.sfa = bayes.sfa.Gen_ZONE  %>% filter(qval < 0.05) %>% pull(LOC) ,
  bayes.final.site.high = bayes.final.Gen_ZONE_FG  %>% filter(qval <= 3.5e-06) %>% pull(LOC) , 
  bayes.final.sfa.high = bayes.final.Gen_ZONE  %>% filter(qval <= 3.5e-06) %>% pull(LOC) , 
  bayes.sfa.site.high = bayes.sfa.Gen_ZONE_FG  %>% filter(qval <= 3.5e-06) %>% pull(LOC) , 
  bayes.sfa.sfa.high = bayes.sfa.Gen_ZONE  %>% filter(qval <= 3.5e-06) %>% pull(LOC) 
)


# Default plot
venn.bayes <- ggVennDiagram::ggVennDiagram(bayes.outliers[c(1,2,5,6)], 
                                           label_alpha = 0#,
                                            #category.names = c("Overall\n(site)", "Overall\n(area)", "Without 3M\n(site)", "Without 3M\n(area)")
                                           ) +
                  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,400)) +
                  ggplot2::scale_color_manual(values = rep("darkgray", 4)) +
                  scale_x_continuous(limits = c(-0.01, 1.01)) +
                  scale_y_continuous(limits = c(.1, .9))# +
                  #labs(title = "BayesScan") #+ theme_bw()
venn.bayes

plot(venn(bayes.outliers[c(1,2,5,6)]),
     quantities = T,
     fill = brewer.pal(n = 4, name = "Set1"),
     alpha = 0.7, edges = F,
     adjust_labels = F,
     main = "")

bayes.outliers.cat <- list(
  "BayesScan - all" = unique( bayes.final.Gen_ZONE_FG  %>% filter(qval < 0.05) %>% pull(LOC) , 
                        bayes.final.Gen_ZONE  %>% filter(qval < 0.05) %>% pull(LOC)) ,
  "BayesScan - high" = unique(bayes.final.Gen_ZONE_FG  %>% filter(qval <= 3.5e-06) %>% pull(LOC) , 
                              bayes.final.Gen_ZONE  %>% filter(qval <= 3.5e-06) %>% pull(LOC) )
)


plot(euler(bayes.outliers.cat, shape = "circle"),
     quantities = T,
     fill = brewer.pal(n = 4, name = "Dark2"),
     alpha = 0.7, edges = F,
     adjust_labels = T)


ggsave(plot = venn.bayes,
       filename = "02_Results/01_Overall_PopGen/01_Bayescan/Venn.allComparions.bayescan.png",
       height = 4.5, width = 6, units = "in", bg = "white")

venn.bayes <- ggVennDiagram::ggVennDiagram(bayes.outliers.high, 
                                           label_alpha = 0,
                                           category.names = c("Overall\n(site)", "Overall\n(area)", "Without 3M\n(site)", "Without 3M\n(area)")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,400)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 4)) +
  scale_x_continuous(limits = c(-0.01, 1.01)) +
  scale_y_continuous(limits = c(.1, .9))# +
#labs(title = "BayesScan") #+ theme_bw()


venn.bayes


# Identify outliers #2 - PCAdapt -----------------------------------------------------------------

# Convertion to plink .bed format

# Read .bed in PCAadapt
pcadapt.final.genotype  <- read.pcadapt("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.bed",
                                      type = "bed")


pcadapt.final.snp <- read.delim("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.bim",
                                header = F) %>% pull(V2)


pcadapt.sfa.genotype  <- read.pcadapt("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1399ind.n13HW.DP.r5.single.sfa.recode.bed",
                                        type = "bed")


pcadapt.sfa.snp <- read.delim("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1399ind.n13HW.DP.r5.single.sfa.recode.bim",
                                header = F) %>% pull(V2)

table(pcadapt.final.snp  == pcadapt.sfa.snp)

# Run pcadapt

K.init <- 10

pcadapt.final <- pcadapt(pcadapt.final.genotype, K =K.init, min.maf = 0.01)
pcadapt.sfa   <- pcadapt(pcadapt.sfa.genotype, K = K.init, min.maf = 0.01)

pcadapt.final$maf %>% min()
pcadapt.final$maf %>% min()

# Check screeplot

plot(pcadapt.final, option = "screeplot") 
plot(pcadapt.sfa, option = "screeplot")

# Check structure

plot(pcadapt.final, option = "scores") 
plot(pcadapt.sfa, option = "scores") 


# K = 2 pour sfa et K = 3 pour final

pcadapt.sfa.k2   <- pcadapt(pcadapt.sfa.genotype , K = 2,  min.maf = 0.01)
pcadapt.final.k3 <- pcadapt(pcadapt.final.genotype , K = 3,  min.maf = 0.01)

plot(pcadapt.sfa.k2, option = "manhattan")
plot(pcadapt.final.k3, option = "manhattan")

hist(pcadapt.sfa.k2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
hist(pcadapt.final.k3$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(pcadapt.sfa.k2, option = "qqplot")
plot(pcadapt.final.k3, option = "qqplot")

# Statistics
#x$pvalues 
alpha <- 0.05

qval.sfa.k2 <- qvalue::qvalue(pcadapt.sfa.k2$pvalues)$qvalues
outliers.sfa.k2 <- which(qval.sfa.k2 < alpha)
outliers.high.sfa.k2 <- which(qval.sfa.k2 < 3.5e-6)


qval.final.k3 <- qvalue::qvalue(pcadapt.final.k3$pvalues)$qvalues
outliers.final.k3 <-  which(qval.final.k3 < alpha)#, which(is.na(qval.final.k3))) %>% unique()
outliers.high.final.k3 <- which(qval.final.k3 < 3.5e-6) 

pcadapt.outliers <- list(

  pcadapt.final =pcadapt.final.snp[outliers.final.k3], 

  pcadapt.sfa =pcadapt.sfa.snp[outliers.sfa.k2],

  pcadapt.final.high =pcadapt.final.snp[outliers.high.final.k3], 
  
  pcadapt.sfa.high =pcadapt.sfa.snp[outliers.high.sfa.k2]
)



pcadapt.outliers.cat <- list(
  
  "PCAdapt - all" =pcadapt.final.snp[outliers.final.k3], 
  "PCAdapt - high" = pcadapt.final.snp[outliers.high.final.k3]

  )


length(outliers.final.k3)

plot(euler(pcadapt.outliers.cat, shape = "circle"),
     quantities = T,
     fill = brewer.pal(n = 4, name = "Dark2"),
     alpha = 0.7, edges = F,
     adjust_labels = T)


plot(euler(c(bayes.outliers.cat[1], pcadapt.outliers.cat[1]), shape = "circle"),
     quantities = T,
     fill = brewer.pal(n = 4, name = "Dark2"),
     alpha = 0.7, edges = F,
     adjust_labels = T)


plot(euler(c(bayes.outliers.cat[2], pcadapt.outliers.cat[2]), shape = "circle"),
     quantities = T,
     fill = brewer.pal(n = 4, name = "Dark2"),
     alpha = 0.7, edges = F,
     adjust_labels = T)


venn.pcadapt <- ggVennDiagram::ggVennDiagram(pcadapt.outliers, 
                                           label_alpha = 0#,
                                           #category.names = c("Overall", "Without 3M")
                                           ) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 4)) #+
  #scale_x_continuous(limits = c(-0.01, 1.01)) +
  #scale_y_continuous(limits = c(.1, .9))# +
#labs(title = "BayesScan") #+ theme_bw()


venn.pcadapt

ggsave(plot = venn.pcadapt,
       filename = "02_Results/01_Overall_PopGen/02_PCAdapt/Venn.allComparions.pcadapt.png",
       height = 4.5, width = 6, units = "in", bg = "white")


venn.pcadapt <- ggVennDiagram::ggVennDiagram(pcadapt.outliers.high, 
                                             label_alpha = 0,
                                             category.names = c("Overall", "Without 3M")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 4)) #+
#scale_x_continuous(limits = c(-0.01, 1.01)) +
#scale_y_continuous(limits = c(.1, .9))# +
#labs(title = "BayesScan") #+ theme_bw()


venn.pcadapt


# Identify outliers - JOIN ------------------------------------------------

# Default plot
venn.overall <- ggVennDiagram(c(pcadapt.outliers[c(1)], bayes.outliers[c(1,2)]), 
                              label_alpha = 0,
                              category.names = c("PCAdapt", "BayeScan\n(area)", "BayeScan\n(site)")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-5, 9)) +
  scale_y_continuous(limits = c(-5, 9))

#theme_bw()
venn.overall

venn.overall.high <- ggVennDiagram(c(pcadapt.outliers.high[c(1)], bayes.outliers.high[c(1,2)]), 
                              label_alpha = 0,
                              category.names = c("PCAdapt", "BayeScan\n(area)", "BayeScan\n(site)")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-5, 9)) +
  scale_y_continuous(limits = c(-5, 9))

#theme_bw()
venn.overall.high


venn.sfa <- ggVennDiagram(c(pcadapt.outliers[c(2)], bayes.outliers[c(3,4)]), 
              label_alpha = 0,
              category.names = c("PCAdapt", "BayeScan\n(area)", "BayeScan\n(site)")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-5, 9)) + 
  scale_y_continuous(limits = c(-5, 9))

venn.sfa.high <- ggVennDiagram(c(pcadapt.outliers.high[c(2)], bayes.outliers.high[c(3,4)]), 
                          label_alpha = 0,
                          category.names = c("PCAdapt", "BayeScan\n(area)", "BayeScan\n(site)")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-5, 9)) + 
  scale_y_continuous(limits = c(-5, 9))

venn.sfa.high 


putative.outliers.final <- c(pcadapt.outliers$pcadapt.final, 
                             bayes.outliers$bayes.final.site, 
                             bayes.outliers$bayes.final.sfa
                             ) %>% unique()

putative.outliers.final %>% length()

putative.outliers.final.high <- c(pcadapt.outliers$pcadapt.final.high, 
                             bayes.outliers$bayes.final.site.high, 
                             bayes.outliers$bayes.final.sfa.high
) %>% unique()

putative.outliers.final.high %>% length()

pcadapt.outliers$pcadapt.final.high %>% unique() %>% length()
c(bayes.outliers$bayes.final.site.high, bayes.outliers$bayes.final.sfa.high) %>% unique() %>% length()

putative.outliers.sfa <- c(pcadapt.outliers$pcadapt.sfa, 
                             bayes.outliers$bayes.sfa.site, 
                             bayes.outliers$bayes.sfa.sfa
) %>% unique()

putative.outliers.sfa.high <- c(pcadapt.outliers$pcadapt.sfa.high, 
                           bayes.outliers$bayes.sfa.site.high, 
                           bayes.outliers$bayes.sfa.sfa.high
) %>% unique()


venn.compa <- ggVennDiagram(list(putative.outliers.final , putative.outliers.sfa), 
                          label_alpha = 0,
                          category.names = c("Overall", "Without 3M")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-4.1, 8.1)) 

venn.compa

venn.compa.high <- ggVennDiagram(list(putative.outliers.final.high , putative.outliers.sfa.high), 
                            label_alpha = 0,
                            category.names = c("Overall", "Without 3M")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "red", limits = c(0,800)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10))+
  scale_x_continuous(limits = c(-4.1, 8.1)) 

venn.compa.high

venn.panel <- ggpubr::ggarrange(venn.compa,
                  ggpubr::ggarrange(venn.overall + theme(legend.position = "none"), 
                                    venn.sfa + theme(legend.position = "none"), 
                                    labels = c("B", "C"),
                                    nrow = 1, ncol = 2),
                  labels = c("A", ""),
                  common.legend = F,
                  nrow = 2, ncol = 1)
venn.panel

ggsave(plot = venn.panel,
       filename = "02_Results/01_Overall_PopGen/Venn.OutlierSNPs.png",
       height = 6, width = 6, units = "in", bg = "white")


length(putative.outliers.final)
length(putative.outliers.sfa)

#save(list = c("putative.outliers.final", "putative.outliers.sfa",
#                 "pcadapt.outliers", "bayes.outliers",
#              "putative.outliers.final.high", "putative.outliers.sfa.high",
#              "pcadapt.outliers.high", "bayes.outliers.high"),
#     file = "02_Results/01_Overall_PopGen/OutlierList.Rdata")


# Load outlier loci -------------------------------------------------------

load("02_Results/01_Overall_PopGen/OutlierList.Rdata")

# Proxy of neutral variation : 
putative.neutral.final <- locNames(gl.final) %>% str_subset(paste(putative.outliers.final, collapse = "|"), negate = T)
putative.neutral.final  %>% length()

putative.PcAdaptOutliers.final <- c(pcadapt.outliers$pcadapt.final) %>% unique()


putative.BayesOutliers.final <- c(bayes.outliers$bayes.final.site, 
                                  bayes.outliers$bayes.final.sfa) %>% unique()

# Allele frequency --------------------------------------------------------

# Allele frequency

# Perform a RDA
gp.final.Gen_ZONE_FG <- genind2genpop(gi.final.Gen_ZONE_FG[pop(gi.final.Gen_ZONE_FG) %nin% pop.too.small,])

# freq for all alleles
freq.allele.final.Gen_ZONE_FG <- adegenet::makefreq(gp.final.Gen_ZONE_FG )
freq.allele.final.Gen_ZONE_FG[1:5, 1:5]

# MAF
evens <- function(x) subset(x, x %% 2 == 0)

freq.MAF.final.Gen_ZONE_FG <- freq.allele.final.Gen_ZONE_FG[, evens(1:dim(freq.allele.final.Gen_ZONE_FG)[2]) ]

dim(freq.MAF.final.Gen_ZONE_FG )# %>% str()
freq.MAF.final.Gen_ZONE_FG[1:5, 1:5]

# Replace NA with most common genotype
#tab.imp.final <- apply(tab.final, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

#sum(is.na(tab.final))
sum(is.na(freq.MAF.final.Gen_ZONE_FG))


pfm.all.pca <- rda(freq.MAF.final.Gen_ZONE_FG, scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pfm.all.pca, type = "barplot", npcs=10, main="PCA Eigenvalues")
plot(pfm.all.pca)

# Performed a RDA using all SNPs

# Proxy of neutral variation : 
putative.neutral.final <- locNames(gl.final) %>% str_subset(paste(putative.outliers.final, collapse = "|"), negate = T)
putative.neutral.final  %>% length()

gp.neutral.snp <- genind2genpop(gi.final.Gen_ZONE_FG[pop(gi.final.Gen_ZONE_FG) %nin% pop.too.small,  loc = putative.neutral.final ])
gp.neutral.snp

pf.neutral.snp<- makefreq(gp.neutral.snp)
pfm.neutral.snp <- pf.neutral.snp[, evens(1:dim(pf.neutral.snp)[2]) ]

pf.neutral.snp[1:5, 1:5]
pfm.neutral.snp[1:5, 1:5]

pfm.neutral.pca <- rda(pfm.neutral.snp, scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pfm.neutral.pca, type = "barplot", npcs=10, main="PCA Eigenvalues")


# Neutral population structure table
PCs <- scores(pfm.neutral.pca, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = row.names(pfm.neutral.snp ), PCs)
colnames(PopStruct) <- c("Gen_ZONE_FG", "PC1", "PC2", "PC3")


# PCA functions -----------------------------------------------------------

plot_pca_eig <- list(geom_bar(stat = "identity"), 
                     geom_line(),
                     scale_fill_viridis_c(),
                     labs(y = "% variance", title = NULL, x = "PC axis"),
                     theme_bw(),
                     theme(axis.text.x = element_blank(), 
                           panel.grid = element_blank(),
                           axis.ticks.length = unit(0, "in"),
                           legend.position = "none"))

plot_pca <-  list(geom_hline(yintercept = 0),
                  geom_vline(xintercept = 0),
                  geom_point(alpha = 0.5, size = 2, pch = 21),
                  theme_bw())


# PCA - ind level --------------------------------------------------------------------

pca.final.14331snps  <- glPca(gl.final , center = TRUE, scale = FALSE,  
                              parallel = TRUE, n.core = 20, nf = 1000)

pca.final.Neutralsnps   <- glPca(gl.final[, locNames(gl.final) %nin% putative.outliers.final] , center = TRUE, scale = FALSE,  
                              parallel = TRUE, n.core = 20, nf = 1000)

pca.final.OutlierHighsnps  <- glPca(gl.final[, locNames(gl.final) %in% putative.outliers.final.high] , center = TRUE, scale = FALSE,  
                              parallel = TRUE, n.core = 20, nf = 1000)

pca.final.Outliersnps  <- glPca(gl.final[, locNames(gl.final) %in% putative.outliers.final] , center = TRUE, scale = FALSE,  
                                parallel = TRUE, n.core = 20, nf = 1000)

pca.final.PcAdaptOutliersnps  <- glPca(gl.final[, locNames(gl.final) %in% putative.PcAdaptOutliers.final] , center = TRUE, scale = FALSE,  
                                parallel = TRUE, n.core = 20, nf = 1000)

pca.final.BayesOutliersnps  <- glPca(gl.final[, locNames(gl.final) %in%  putative.BayesOutliers.final] , center = TRUE, scale = FALSE,  
                                parallel = TRUE, n.core = 20, nf = 1000)


#save(list = c("pca.final.14331snps", "pca.final.Neutralsnps", "pca.final.Outliersnps", 
#                "pca.final.BayesOutliersnps",  "pca.final.PcAdaptOutliersnps"   
#            ),
#     file = here("02_Results/01_Overall_PopGen/03_PCA/PCA.Rdata"))

load(here("02_Results/01_Overall_PopGen/03_PCA/PCA.Rdata"))

## Functions

## Eig var

final.14331snps.var <- pca_var(pca.final.14331snps, nInd(gl.final)-1) %>% 
  ggplot(aes(x =axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig
final.14331snps.var 

final.Neutralsnps.var <- pca_var(pca.final.Neutralsnps, nInd(gl.final)-1) %>% 
  ggplot(aes(x =axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig
final.Neutralsnps.var 

final.Outliersnps.var <- pca_var(pca.final.Outliersnps, length(putative.outliers.final)-1) %>% 
  ggplot(aes(x =axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig
final.Outliersnps.var 

# PCA
gPCA.final.14331snps <- pca.final.14331snps %>% QuickPop::pca_scoretable(naxe = 3) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M"))) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = RegionAssesment)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.5, size = 0.5) +  
  stat_ellipse(aes(col = RegionAssesment), cex = 1.1)+
  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"), 
                      labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+    
  #annotate("text",  x=Inf, y = Inf, label = paste("All snps:",  nLoc(gl.final)), vjust=1, hjust=1) +
  facet_wrap(~"Complete SNP panel") +  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
       x = paste0("PC1 (", QuickPop::pca_var(pca.final.14331snps)$p.eig[1] %>% round(3) *100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(pca.final.14331snps)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill="white"))
gPCA.final.14331snps

gPCA.final.Outliersnps <- pca.final.Outliersnps %>% QuickPop::pca_scoretable(naxe = 3) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M"))) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = RegionAssesment)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +

  geom_point(alpha = 0.5, size = 0.5) +  
  stat_ellipse(aes(col = RegionAssesment), cex = 1.1)+
  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+ 
  #annotate("text",  x=-Inf, y = Inf, label = paste("Outlier snps:",  length(putative.outliers.final)), vjust=1, hjust=0) +
  facet_wrap(~"Outlier SNP panel") + 
  labs(#title = paste("Outlier snps:",  length(putative.outliers.final)),
       x = paste0("PC1 (", QuickPop::pca_var(pca.final.Outliersnps)$p.eig[1] %>% round(3) *100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(pca.final.Outliersnps)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill="white"))
gPCA.final.Outliersnps

gPCA.final.Neutralsnps <- pca.final.Neutralsnps %>% QuickPop::pca_scoretable(naxe = 3) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  mutate(RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M"))) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = RegionAssesment)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +

  geom_point(alpha = 0.5, size = 0.5) +  
  stat_ellipse(aes(col = RegionAssesment), cex = 1.1)+
  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+ 
  facet_wrap(~"Neutral SNP panel") +  
  #annotate("text",  x=Inf, y = Inf, label =paste("Neutral snps:",  nLoc(gl.final) - length(putative.outliers.final)), vjust=1, hjust=1) +
  
  labs(#title = paste("Neutral snps:",  nLoc(gl.final) - length(putative.outliers.final)),
       x = paste0("PC1 (", QuickPop::pca_var(pca.final.Neutralsnps)$p.eig[1] %>% round(3) *100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(pca.final.Neutralsnps)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill="white"))
gPCA.final.Neutralsnps

# PCA - sampling location -------------------------------------------------

# Allele frequency

# Perform a RDA

# MAF
evens <- function(x) subset(x, x %% 2 == 0)

# ALL snps

gp.all.snp <- genind2genpop(gi.final.Gen_ZONE_FG[pop(gi.final.Gen_ZONE_FG) %nin% pop.too.small ])
gp.all.snp

pf.all.snp<- makefreq(gp.all.snp)
pfm.all.snp <- pf.all.snp[, evens(1:dim(pf.all.snp)[2]) ]

pf.all.snp[1:5, 1:5]
pfm.all.snp[1:5, 1:5]

pfm.all.pca <- rda(pfm.all.snp, scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pfm.all.pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

all.var <- data.frame(axis = 1:length(pfm.all.pca$CA$eig),
           eig = pfm.all.pca$CA$eig) %>% mutate(p.eig = eig / sum(eig))%>% 
  ggplot(aes(x =axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig

res.all <- as.data.frame(scores(pfm.all.pca, display="sites", scaling=1, choices = c(1,2,3,4))) %>% 
  mutate(Gen_ZONE_FG = dimnames(scores(pfm.all.pca, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct) %>% 
  mutate(Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                "SFA-9", "SFA-10", "SFA-12",
                                                "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                "NAFO-3M")),
         RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M")))
#res2.scotian <- dplyr::select(res.scotian, -c(RegionAssesment, Gen_ZONE, Gen_ZONE_FG))  

perc.all <- round(100*(summary(pfm.all.pca)$cont$importance[2, 1:2]), 1)


gg.pca.all <-  res.all %>% 
  ggplot(aes(x =PC1, y = PC2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(fill = RegionAssesment), pch = 21, cex = 2, alpha = 0.5)+
  scale_fill_manual(name = "", values = c("black","blue", "darkorange","red", "magenta"),
                    labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  
  # geom_segment(data = arrow.env.df, aes(x = x0, y = y0, xend = xmax * arrow.factor, yend= ymax * arrow.factor),
  #              arrow = arrow(), size = 0.5)+
  # geom_label_repel(data = arrow.env.df , aes(label = name, x = xmax *arrow.factor, y = ymax * arrow.factor), size = 3, max.overlaps = 20
  #  ) +
  geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-")), size = 3, max.overlaps = 1
  ) +
#  annotate("text",  x=Inf, y = -Inf, label = paste("All snps:",  nLoc(gl.final) ), vjust=0, hjust=1) +
  facet_wrap(~"Complete SNP panel") +  
  labs(x = paste0("PC1 (", perc.all[1], "%)"), 
       y = paste0("PC2 (", perc.all[2], "%)")) + 
  #annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS.northern)$adj.r.squared %>% as.numeric() %>%  round(digits = 3)),
  #          x = -0.2, y = 0.25, vjust = "inward", hjust = "inward") +
  
  theme_bw(base_size = 11)  +
  theme(strip.background = element_rect(fill="white"))

gg.pca.all

# Proxy of neutral variation : 
putative.neutral.final <- locNames(gl.final) %>% str_subset(paste(putative.outliers.final, collapse = "|"), negate = T)
putative.neutral.final  %>% length()


gp.neutral.snp <- genind2genpop(gi.final.Gen_ZONE_FG[pop(gi.final.Gen_ZONE_FG) %nin% pop.too.small,  loc = putative.neutral.final ])
gp.neutral.snp

pf.neutral.snp<- makefreq(gp.neutral.snp)
pfm.neutral.snp <- pf.neutral.snp[, evens(1:dim(pf.neutral.snp)[2]) ]

pf.neutral.snp[1:5, 1:5]
pfm.neutral.snp[1:5, 1:5]

pfm.neutral.pca <- rda(pfm.neutral.snp, scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pfm.neutral.pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

res.neutral <- as.data.frame(scores(pfm.neutral.pca, display="sites", scaling=1)) %>% 
  mutate(Gen_ZONE_FG = dimnames(scores(pfm.neutral.pca, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct) %>% 
  mutate(Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                "SFA-9", "SFA-10", "SFA-12",
                                                "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                "NAFO-3M")),
         RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M")))
#res2.scotian <- dplyr::select(res.scotian, -c(RegionAssesment, Gen_ZONE, Gen_ZONE_FG))  

perc.neutral <- round(100*(summary(pfm.neutral.pca)$cont$importance[2, 1:2]), 1)

gg.pca.neutral <-  res.neutral %>% 
  ggplot(aes(x =PC1, y = PC2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(fill = RegionAssesment), pch = 21, cex = 2, alpha = 0.5)+
  scale_fill_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+ 
  
  geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-")), size = 3, max.overlaps = 1
  ) +
  #annotate("text",  x=Inf, y = -Inf, label = paste("Neutral snps:",  nLoc(gl.final) -  length(putative.outliers.final)), vjust=0, hjust=1) +
  facet_wrap(~"Neutral SNP panel") +  
  labs(x = paste0("PC1 (", perc.neutral[1], "%)"), 
       y = paste0("PC2 (", perc.neutral[2], "%)")) + 
  #annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS.northern)$adj.r.squared %>% as.numeric() %>%  round(digits = 3)),
  #          x = -0.2, y = 0.25, vjust = "inward", hjust = "inward") +
  
  theme_bw(base_size = 11)  +
  theme(strip.background = element_rect(fill="white"))

gg.pca.neutral

# Outliers

gp.outlier.snp <- genind2genpop(gi.final.Gen_ZONE_FG[pop(gi.final.Gen_ZONE_FG) %nin% pop.too.small,  loc = putative.outliers.final ])
gp.outlier.snp

pf.outlier.snp<- makefreq(gp.outlier.snp)
pfm.outlier.snp <- pf.outlier.snp[, evens(1:dim(pf.outlier.snp)[2]) ]

pf.outlier.snp[1:5, 1:5]
pfm.outlier.snp[1:5, 1:5]

pfm.outlier.pca <- rda(pfm.outlier.snp, scale=T) # PCA in vegan uses the rda() call without any predictors
screeplot(pfm.outlier.pca, type = "barplot", npcs=10, main="PCA Eigenvalues")


res.outlier <- as.data.frame(scores(pfm.outlier.pca, display="sites", scaling=1)) %>% 
  mutate(Gen_ZONE_FG = dimnames(scores(pfm.outlier.pca, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct) %>% 
  mutate(Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                "SFA-9", "SFA-10", "SFA-12",
                                                "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                "NAFO-3M")),
         RegionAssesment = factor(RegionAssesment, levels = c("Northern", "NFL", "GSL", "Scotian", "NAFO-3M")))
#res2.scotian <- dplyr::select(res.scotian, -c(RegionAssesment, Gen_ZONE, Gen_ZONE_FG))  

perc.outlier <- round(100*(summary(pfm.outlier.pca)$cont$importance[2, 1:2]), 1)

gg.pca.outlier <-  res.outlier %>% 
  ggplot(aes(x =PC1, y = PC2)) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes(fill = RegionAssesment), pch = 21, cex = 2, alpha = 0.5)+
  scale_fill_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+ 
  
  # geom_segment(data = arrow.env.df, aes(x = x0, y = y0, xend = xmax * arrow.factor, yend= ymax * arrow.factor),
  #              arrow = arrow(), size = 0.5)+
  # geom_label_repel(data = arrow.env.df , aes(label = name, x = xmax *arrow.factor, y = ymax * arrow.factor), size = 3, max.overlaps = 20
  #  ) +
  geom_text_repel(aes(label = Gen_ZONE_FG %>% str_remove("SFA-")), size = 3, max.overlaps = 1
  ) +
  #annotate("text",  x=Inf, y = -Inf, label = paste("Outlier snps:",    length(putative.outliers.final)), vjust=0, hjust=1) +
  facet_wrap(~"Outlier SNP panel") + 
  labs(x = paste0("PC1 (", perc.outlier[1], "%)"), 
       y = paste0("PC2 (", perc.outlier[2], "%)")) + 
  #annotate("text", label = paste("R²-adjusted =", RsquareAdj(rdaS.northern)$adj.r.squared %>% as.numeric() %>%  round(digits = 3)),
  #          x = -0.2, y = 0.25, vjust = "inward", hjust = "inward") +
  
  theme_bw(base_size = 11)  +
  theme(strip.background = element_rect(fill="white"))

gg.pca.outlier

# FST ---------------------------------------------------------------------

final.14331snps_region.FST <- dartR::gl.fst.pop(gl.final.Region, nboots = 999, percent = 95, nclusters = 20)
final.14331snps_area.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE, nboots = 999, percent = 95, nclusters = 20)
final.14331snps_site.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE_FG, nboots = 999, percent = 95, nclusters = 20)

final.Outliersnps_region.FST <- dartR::gl.fst.pop(gl.final.Region[, locNames(gl.final.Gen_ZONE) %in% putative.outliers.final] , nboots = 999, percent = 95, nclusters = 20)
final.Outliersnps_area.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE[, locNames(gl.final.Gen_ZONE) %in% putative.outliers.final] , nboots = 999, percent = 95, nclusters = 20)
final.Outliersnps_site.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE_FG[, locNames(gl.final.Gen_ZONE_FG) %in% putative.outliers.final], nboots = 999, percent = 95, nclusters = 20)

final.Neutralsnps_region.FST <- dartR::gl.fst.pop(gl.final.Region[, locNames(gl.final.Gen_ZONE) %nin% putative.outliers.final] , nboots = 999, percent = 95, nclusters = 20)
final.Neutralsnps_area.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE[, locNames(gl.final.Gen_ZONE) %nin% putative.outliers.final] , nboots = 999, percent = 95, nclusters = 20)
final.Neutralsnps_site.FST <- dartR::gl.fst.pop(gl.final.Gen_ZONE_FG[, locNames(gl.final.Gen_ZONE_FG) %nin% putative.outliers.final], nboots = 999, percent = 95, nclusters = 20)

#save(list = c("final.14331snps_area.FST", "final.14331snps_site.FST", "final.14331snps_region.FST",
#              "final.Outliersnps_area.FST", "final.Outliersnps_site.FST", "final.Outliersnps_region.FST",
#              "final.Neutralsnps_area.FST", "final.Neutralsnps_site.FST", "final.Neutralsnps_region.FST"#,
#              ),
#     file = file.path(here::here(), "02_Results/01_Overall_PopGen/04_FST/Fst.All.Rdata"))

load(file.path(here::here(), "02_Results/01_Overall_PopGen/04_FST/Fst.All.Rdata"))


# Compared

# Functions

table.fst <- function(fst){
  res <-  fst$Bootstraps %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")
  
  return(res)
  
}

heat.fst <- function(fst){
  res <- bind_rows(table.fst(fst),
                   table.fst(fst) %>% mutate(Pop1.int = Population2,
                                             Pop2.int = Population1,
                                             Population1 = Pop1.int,
                                             Population2 = Pop2.int) %>% 
                     dplyr::select(-c(Pop1.int, Pop2.int))
  )
  
  return(res)         
  
}



heat.graph.panel <- function(heat){
  graph <-  heat %>%  
    ggplot(aes(x=Population1, y=Population2, fill=Fst, shapp = Sign)) + 
    geom_tile(colour=  "gray") +
    geom_point(aes(x = Population1, y=Population2, pch = Sign), size = 3)+
    #scale_y_discrete(limits=rev) +
    #scale_x_discrete(limits=rev) +
    scale_fill_gradient(low = "ivory1", high = "red3")+
    # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
    scale_shape_manual(values = c("","*"), guide = "none") +
    #scale_fill_distiller(palette = "Spectral") +
    labs(x = NULL, y = NULL) +
    #theme_bw() +
    facet_grid(SFA2~SFA1, space = "free", scale = "free") +
    theme_minimal() + 
    theme(#axis.text.x = element_blank(),
      strip.text = element_text(angle = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "cm"),
      panel.border = element_rect(fill = NA, colour = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
      legend.position = "right",
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))
  print(graph)
  return(graph)
}

# Results - Region

final.14331snps_region.FST %>% table.fst() %>% 
  summarise(Mean.NAFO = mean(Fst),
            sd.NAFO = sd(Fst),
            Min.NFAO = min(Fst),
            Max.NAFO = max(Fst),
            Max.Pvalue = round(max(`p-value`),))


g0.fst <- final.14331snps_region.FST %>% heat.fst() %>%  #%>% left_join(pop.data %>% dplyr::select(Population1 = RegionAssesment, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  #left_join(pop.data %>% dplyr::select(Population2 = RegionAssesment, SFA2 = RegionAssesment)  %>% distinct()) %>% #pull(SFA1) %>% unique()
  mutate(Population1 = ifelse(Population1 == "NAFO-3M", "FLC",
                      ifelse(Population1== "Northern", "ARC",
                      ifelse(Population1 == "Scotian", "SS",       
                             Population1))),
         Population2 = ifelse(Population2 == "NAFO-3M", "FLC",
                      ifelse(Population2== "Northern", "ARC",
                             ifelse(Population2 == "Scotian", "SS",  Population2))),
         Population1 = factor(Population1, levels =  c("ARC", "NFL", "GSL", "SS", "FLC")),
         Population2 = factor(Population2, levels = c("ARC", "NFL", "GSL", "SS", "FLC")),
        #SFA1 = factor(SFA1, levels = c("FC", "Scotian", "GSL", "NFL", "Artic")),
         #SFA2= factor(SFA2, levels = c("Artic", "NFL", "GSL", "Scotian", "FC")),
         Sign = ifelse(`p-value` <= 0.05, "*", "")) %>% 
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) + 
  geom_tile(colour=  "gray") +
 # geom_point(aes(x = Population1, y=Population2), size = 3)+
  geom_text(aes(label = paste0(format(round(Fst, dig = 4 ), scientific=F), Sign))) +
  #scale_y_discrete(limits=rev) +
  #scale_x_discrete(limits=rev) +
  scale_fill_gradient(low = "ivory1", high = "red3")+
  # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
  scale_shape_manual(values = c("","*"), guide = "none") +
  #scale_fill_distiller(palette = "Spectral") +
  labs(x = NULL, y = NULL) +
  #theme_bw() +
  #facet_grid(SFA2~SFA1, space = "free", scale = "free") +
  theme_minimal() + 
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(angle = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))

#g0.fst <-  g1.fst + scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev)  +    scale_fill_viridis_c(option = "C")
g0.fst


ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "04_FST", "Fst.byRegion.png"), 
       plot = g0.fst,
       height = 4, width = 4, units = "in")

# Results - Overall

final.14331snps_area.FST %>% table.fst() %>% mutate(Comp = paste(Population1, Population1)) %>% filter(str_detect(Comp, "NAFO-3M")) %>% 
  summarise(Mean.NAFO = mean(Fst),
            sd.NAFO = sd(Fst),
            Min.NFAO = min(Fst),
            Max.NAFO = max(Fst))

final.14331snps_area.FST %>% table.fst() %>% 
  filter(`p-value`< 0.05) %>% arrange(Fst) #%>% nrow()

g1.fst <- final.14331snps_area.FST %>% heat.fst() %>% left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE, SFA2 = RegionAssesment)  %>% distinct()) %>% #pull(SFA1) %>% unique()
  mutate(SFA1 = ifelse(SFA1 == "NAFO-3M", "FLC",
                ifelse(SFA1== "Northern", "ARC", 
                       ifelse(SFA1 == "Scotian", "SS",  
                       SFA1))),
         
         SFA2 = ifelse(SFA2 == "NAFO-3M", "FLC",
                ifelse(SFA2== "Northern", "ARC",ifelse(SFA2 == "Scotian", "SS",  SFA2))),
         Population1 = factor(Population1, levels = pop.data %>% pull(Gen_ZONE) %>% levels()),
         Population2 = factor(Population2, levels =pop.data %>% pull(Gen_ZONE) %>% levels()),
         SFA2 = factor(SFA2, levels = c("FLC", "SS", "GSL", "NFL", "ARC")),
         SFA1= factor(SFA1, levels = c("ARC", "NFL", "GSL", "SS", "FLC")),
         Sign = ifelse(`p-value` <= 0.05, "*", "")) %>% 
 # heat.graph.panel()# + theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
                    #         strip.text.y = element_text(angle = 0, vjust= 0.5, hjust = 0))
#g1.fst <-  g1.fst + scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev)  +    scale_fill_viridis_c(option = "C")
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) + 
  geom_tile(colour=  "gray") +
  # geom_point(aes(x = Population1, y=Population2), size = 3)+
  geom_text(aes(label = paste0(format(round(Fst, dig = 4 ), scientific=F), Sign)), cex = 3) +
  #scale_y_discrete(limits=rev) +
  #scale_x_discrete(limits=rev) +
  scale_fill_gradient(low = "ivory1", high = "red3")+
  # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
  scale_shape_manual(values = c("","*"), guide = "none") +
  #scale_fill_distiller(palette = "Spectral") +
  labs(x = NULL, y = NULL) +
  #theme_bw() +
  facet_grid(SFA2~SFA1, space = "free", scale = "free") +
  theme_minimal() + 
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(angle = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))
  
  g1.fst

  
  ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "04_FST", "Fst.overall_numb.png"), 
         plot = g1.fst,
         height =9, width = 9, units = "in")  
  
# By region

final.14331snps_site.FST %>% table.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = RegionAssesment)  %>% distinct()) %>% 
  dplyr::filter(SFA1 == SFA2) %>% 
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst))

# SFA-16-1

fst.df.16 <- final.14331snps_site.FST %>% table.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = RegionAssesment)  %>% distinct()) %>% 
  dplyr::filter(SFA1 == SFA2,
                SFA1 == "Scotian") %>% 
  mutate(Group = ifelse(Population1 == "SFA-16-1" | Population2 == "SFA-16-1", "Focal", "Other"))
 
fst.df.16 %>% group_by(Group) %>%  summarise(mean = mean(Fst),
                                             sd = sd(Fst)) 
 
nrow(fst.df.16)

t.test(fst.df.16 %>% dplyr::filter(Group == "Focal") %>% pull(Fst),
       fst.df.16 %>% dplyr::filter(Group == "Other") %>% pull(Fst)#,
       #alternative = "greater"
       )

# SFA-5-3

fst.df.53 <- final.14331snps_site.FST %>% table.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = RegionAssesment)  %>% distinct()) %>% 
  dplyr::filter(SFA1 == SFA2,
                SFA1 == "NFL") %>% 
  filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  
  mutate(Group = ifelse(Population1 == "SFA-5-3" | Population2 == "SFA-5-3", "Focal", "Other"))

fst.df.53 %>% group_by(Group) %>%  summarise(mean = mean(Fst),
                                             sd = sd(Fst)) 


t.test(fst.df.53 %>% dplyr::filter(Group == "Focal") %>% pull(Fst),
       fst.df.53 %>% dplyr::filter(Group == "Other") %>% pull(Fst)#,
       #alternative = "greater"
)


# SFA-E2

fst.df.E2 <- final.14331snps_site.FST %>% table.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = RegionAssesment)  %>% distinct()) %>% 
  dplyr::filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  
  dplyr::filter(SFA1 == SFA2,
                SFA1 == "Northern") %>% 
  filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  
  mutate(Group = ifelse(Population1 == "EAZ-2" | Population2 == "EAZ-2", "Focal", "Other"))

fst.df.E2 %>% group_by(Group) %>%  summarise(mean = mean(Fst),
                                             sd = sd(Fst)) 

t.test(fst.df.E2 %>% dplyr::filter(Group == "Focal") %>% pull(Fst),
       fst.df.E2 %>% dplyr::filter(Group == "Other") %>% pull(Fst)#,
       #alternative = "greater"
)

# Laurentian chenal


fst.df.chenal <- final.14331snps_site.FST %>% table.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = RegionAssesment) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = RegionAssesment)  %>% distinct()) %>% 
  dplyr::filter(#SFA1 == SFA2,
                SFA1 %in% c("Scotian", "GSL"),
                SFA2 %in% c("Scotian", "GSL"),
                Population1 %in% c("SFA-13-3", "SFA-13-5") | Population2 %in% c("SFA-13-3", "SFA-13-5") ,
                Population1 != "SFA-16-1", Population2 != "SFA-16-1"
                )%>% #View() 
  mutate(Group = ifelse(Population1 %in% c("SFA-13-3", "SFA-13-5") & Population2 %in% c("SFA-13-3", "SFA-13-5"), "Intra",
                        ifelse(SFA1 == "GSL" |  SFA2 == "GSL", "Focal", "Other")))

fst.df.chenal %>% group_by(Group) %>%  summarise(mean = mean(Fst),
                                             sd = sd(Fst)) 

nrow(fst.df.16)

t.test(fst.df.chenal %>% dplyr::filter(Group == "Focal") %>% pull(Fst),
       fst.df.chenal %>% dplyr::filter(Group == "Other") %>% pull(Fst)#,
       #alternative = "greater"
)



g3 <- final.14331snps_site.FST %>% heat.fst() %>%  mutate(Population1 = factor(Population1),
                                             Population2 = factor(Population2),
                                             SFA1 = sapply(str_split(Population1, "-"), `[`, 2), 
                                             SFA2 = sapply(str_split(Population2, "-"), `[`, 2),
                                             #SFA1 = ifelse(SFA1 %in% c("14", "15"), "14-15", SFA1),
                                             SFA1 = factor(SFA1, levels = c("16", "15", "14", "13")),
                                             #SFA2 = ifelse(SFA2 %in% c("14", "15"), "14-15", SFA2),
                                             SFA2= factor(SFA2, levels = c("13", "14", "15", "16")),
                                             Sign = ifelse(`p-value` <= 0.05, "Sign.", "Non-Sign."),
                                             Fst = ifelse(Fst < 0,0, Fst))  %>% 
  filter(!is.na(SFA1), !is.na(SFA2)) %>% 
  
  heat.graph.panel() + labs(title = "SS")
g3 <- g3 +  scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev) 
g3

g4 <- final.14331snps_site.FST  %>% heat.fst() %>%  mutate(Population1 = factor(Population1),#, levels = pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
                                         Population2 = factor(Population2),#, levels =pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
                                         SFA1 = sapply(str_split(Population1, "-"), `[`, 2), 
                                         SFA2 = sapply(str_split(Population2, "-"), `[`, 2),
                                         #SFA1 = ifelse(SFA1 %in% c("14", "15"), "14-15", SFA1),
                                         SFA1 = factor(SFA1, levels = c("12", "10", "9", "8")),
                                         #SFA2 = ifelse(SFA2 %in% c("14", "15"), "14-15", SFA2),
                                         SFA2= factor(SFA2, levels = c("8", "9", "10", "12")),
                                         Sign = ifelse(`p-value` <= 0.05, "Sign.", "Non-Sign."),
                                         Fst = ifelse(Fst < 0,0, Fst))  %>% 
  filter(!is.na(SFA1), !is.na(SFA2)) %>% 
  heat.graph.panel() + labs(title = "GSL")
g4 

g5.a <- final.14331snps_site.FST %>% heat.fst() %>%  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = Gen_ZONE) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = Gen_ZONE)  %>% distinct()) %>% 
  
  
  mutate(Population1 = factor(Population1),#, levels = pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         Population2 = factor(Population2),#, levels =pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         SFA1 =  SFA1 %>% str_remove("SFA-"), 
         SFA2 =  SFA2 %>% str_remove("SFA-"), 
         #SFA2 = sapply(str_split(Population2, "-"), `[`, 2),
         #SFA1 = ifelse(SFA1 %in% c("14", "15"), "14-15", SFA1),
         SFA1 = factor(SFA1, levels = c("EAZ", "WAZ", "4")),
         #SFA2 = ifelse(SFA2 %in% c("14", "15"), "14-15", SFA2),
         SFA2= factor(SFA2, levels = c("4", "WAZ", "EAZ")),
         Sign = ifelse(`p-value` <= 0.05, "Sign.", "Non-Sign."),
         Fst = ifelse(Fst < 0,0, Fst))  %>% 
  filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  filter(!is.na(SFA1), !is.na(SFA2)) %>% 
  heat.graph.panel() + labs(title = "ARC")
g5.a 

g5.b <- final.14331snps_site.FST %>% heat.fst() %>%  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = Gen_ZONE) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = Gen_ZONE)  %>% distinct()) %>% 
  
  
  mutate(Population1 = factor(Population1),#, levels = pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         Population2 = factor(Population2),#, levels =pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         SFA1 =  SFA1 %>% str_remove("SFA-"), 
         SFA2 =  SFA2 %>% str_remove("SFA-"), 
         #SFA2 = sapply(str_split(Population2, "-"), `[`, 2),
         #SFA1 = ifelse(SFA1 %in% c("14", "15"), "14-15", SFA1),
         SFA1 = factor(SFA1, levels = c("5", "6", "7")),
         #SFA2 = ifelse(SFA2 %in% c("14", "15"), "14-15", SFA2),
         SFA2= factor(SFA2, levels = c("7", "6", "5")),
         Sign = ifelse(`p-value` <= 0.05, "Sign.", "Non-Sign."),
         Fst = ifelse(Fst < 0,0, Fst))  %>% 
  filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  filter(!is.na(SFA1), !is.na(SFA2)) %>% 
  heat.graph.panel() + labs(title = "NFL")
g5.b  

g6 <- final.14331snps_site.FST %>% heat.fst() %>%  left_join(pop.data %>% dplyr::select(Population1 = Gen_ZONE_FG, SFA1 = Gen_ZONE) %>% distinct()  ) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Gen_ZONE_FG, SFA2 = Gen_ZONE)  %>% distinct()) %>% 
  
  
  mutate(Population1 = factor(Population1),#, levels = pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         Population2 = factor(Population2),#, levels =pop.data %>% pull(Gen_ZONE_FG) %>% levels()),
         SFA1 =  SFA1 %>% str_remove("NAFO-"), 
         SFA2 =  SFA2 %>% str_remove("NAFO-"), 
         #SFA2 = sapply(str_split(Population2, "-"), `[`, 2),
         #SFA1 = ifelse(SFA1 %in% c("14", "15"), "14-15", SFA1),
         SFA1 = factor(SFA1, levels = c("3M")),
         #SFA2 = ifelse(SFA2 %in% c("14", "15"), "14-15", SFA2),
         SFA2= factor(SFA2, levels = c("3M")),
         Sign = ifelse(`p-value` <= 0.05, "Sign.", "Non-Sign."),
         Fst = ifelse(Fst < 0,0, Fst))  %>% 
  filter(Population1  %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3"), Population2 %nin% c("SFA-0-1", "SFA-7-7", "WAZ-3")) %>% 
  filter(!is.na(SFA1), !is.na(SFA2)) %>% 
  heat.graph.panel() + labs(title = "FLC")
g6  


graph.area.Fst <- ggpubr::ggarrange(#g1 + scale_fill_gradient(low = "ivory1", high = "red3", limits = c(-0.001, 0.025)), 
  ggpubr::ggarrange(
                                    g3  + scale_fill_gradient(low = "white", high = "blue", limits = c(0,0.0038)) + theme(legend.position = "none"),
                                    g5.b + scale_fill_gradient(low = "white", high = "blue", limits = c(0,0.0038))+ theme(legend.position = "none"),
                                    nrow = 1, ncol = 2,  labels = LETTERS[1:2], align = "hv"),
  ggpubr::ggarrange(
                                    
                                    
                                    g5.a  + scale_fill_gradient(low = "white", high = "blue", limits = c(0,0.0038)) + theme(legend.position = "none"),
                                    g4 + scale_fill_gradient(low = "white", high = "blue", limits = c(0,0.0038)) + theme(legend.position = "none") ,
                                    g6  + scale_fill_gradient(low = "white", high = "blue", limits = c(0,0.0038)) ,
                                    nrow = 1, ncol = 3,  labels = LETTERS[3:5], align = "hv"),
                               nrow = 2, ncol = 1, common.legend = T, heights = c(2,1))


graph.area.Fst 

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "04_FST",  "Fst.byArea_20230627.png"), 
       plot = graph.area.Fst,
       height = 8, width = 8, units = "in")   


# He vs Ho ----------------------------------------------------------------

basic_stat_Region <- hierfstat::basic.stats(gi.final.Region, diploid = TRUE)
basic_stat_Gen_ZONE <- hierfstat::basic.stats(gi.final.Gen_ZONE, diploid = TRUE)
basic_stat_Gen_ZONE_FG <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG, diploid = TRUE)

#save(list = c("basic_stat_Region", "basic_stat_Gen_ZONE", "basic_stat_Gen_ZONE_FG"), 
#     file = file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/HeHo.RData"))

load(file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/HeHo.RData"))

# Homemade function to extract important infos
extract.basic.stats <- function(basic_stat_ALL){
  data.frame(ID = dimnames(basic_stat_ALL$Fis)[[2]],
           Fis = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE),
           Fis.IC95.low = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
           Fis.IC95.high = apply(basic_stat_ALL$Fis, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2]),
           
           He = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE),
           He.IC95.low = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
           He.IC95.high = apply(basic_stat_ALL$Hs, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2]),
           
           Ho = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE),
           Ho.IC95.low = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[1]),
           Ho.IC95.high = apply(basic_stat_ALL$Ho, MARGIN = 2, FUN = function(x) confint(lm(x ~ 1), level=0.95)[2])
           
           )
}

Hstat_Region      <- extract.basic.stats(basic_stat_Region)
Hstat_Gen_ZONE    <- extract.basic.stats(basic_stat_Gen_ZONE)
Hstat_Gen_ZONE_FG <- extract.basic.stats(basic_stat_Gen_ZONE_FG)

write_csv(Hstat_Region , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Region.csv"))
write_csv(Hstat_Gen_ZONE , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Gen_ZONE.csv"))
write_csv(Hstat_Gen_ZONE_FG , file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/Hstat_Gen_ZONE_FG.csv"))

# Stats

# 16-1 

df.he.16 <- rbind(Hstat_Gen_ZONE_FG %>% dplyr::select(ID, value = Ho, IC95.high = Ho.IC95.high, IC95.low =  Ho.IC95.low ) %>%
        mutate(name = "Ho"),
      Hstat_Gen_ZONE_FG %>% dplyr::select(ID, value = He, IC95.high = He.IC95.high, IC95.low =  He.IC95.low ) %>%
        mutate(name = "He") ) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct(),
            by = c("ID" = "Gen_ZONE_FG")) %>% 
  dplyr::filter(RegionAssesment == "Scotian",
                name == "He") %>% 
  mutate(Group = ifelse(ID == "SFA-16-1", "Focal", "Other"))

df.he.16 %>% group_by(Group) %>% summarise(mean = mean(value),
                                           sd= sd(value),
                                           maxCI = max(IC95.high),
                                           minCI = min (IC95.low ))


# 5-3 

df.he.53 <- rbind(Hstat_Gen_ZONE_FG %>% dplyr::select(ID, value = Ho, IC95.high = Ho.IC95.high, IC95.low =  Ho.IC95.low ) %>%
                    mutate(name = "Ho"),
                  Hstat_Gen_ZONE_FG %>% dplyr::select(ID, value = He, IC95.high = He.IC95.high, IC95.low =  He.IC95.low ) %>%
                    mutate(name = "He") ) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct(),
            by = c("ID" = "Gen_ZONE_FG")) %>% 
  dplyr::filter(RegionAssesment == "NFL",
                name == "He",
                ID != "SFA-7-7") %>% 
  mutate(Group = ifelse(ID == "SFA-5-3", "Focal", "Other"))


df.he.53 %>% group_by(Group) %>% summarise(mean = mean(value),
                                           sd= sd(value),
                                           maxCI = max(IC95.high),
                                           minCI = min (IC95.low ))


# Give more weigth

#basic_stat_Region <- hierfstat::basic.stats(gi.final.Region, diploid = TRUE)
basic_stat_Gen_ZONE_FG.neutral <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = putative.neutral.final], diploid = TRUE)
basic_stat_Gen_ZONE_FG.outlier <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = putative.outliers.final], diploid = TRUE)
basic_stat_Gen_ZONE_FG.adaptative <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = outliers_rdadapt_env.cor], diploid = TRUE)

#basic_stat_Gen_ZONE_FG <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG, diploid = TRUE)

#save(list = c("basic_stat_Gen_ZONE_FG.neutral", "basic_stat_Gen_ZONE_FG.outlier", "basic_stat_Gen_ZONE_FG.adaptative" ), 
#     file = file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/HeHo_neutral_outlier.RData"))

load(file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/HeHo_neutral_outlier.RData"))


Hstat_Gen_ZONE.comp <- bind_rows(extract.basic.stats(basic_stat_Gen_ZONE_FG.neutral)  %>% mutate(variable = "Neutral"),
                                 extract.basic.stats(basic_stat_Gen_ZONE_FG.outlier)  %>% mutate(variable = "Outlier"),
                                 extract.basic.stats(basic_stat_Gen_ZONE_FG.adaptative)  %>% mutate(variable = "Adaptative"),
                                 extract.basic.stats(basic_stat_Gen_ZONE_FG)  %>% mutate(variable = "Full")
                                 
                                 
) %>% left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct(.keep_all = T), by = c("ID" = "Gen_ZONE_FG")) %>% 
  mutate(Region = ifelse(Gen_ZONE %in% c("SFA-4", "SFA-5", "SFA-5", "SFA-7"), "NFL", 
                         ifelse(Gen_ZONE == "NAFO-3M", "3M",  RegionAssesment)),
         Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
         
         variable = factor(variable, levels = c("Full", "Neutral", "Outlier", "Adaptative"))) 

Hstat_Gen_ZONE.comp 


# One-way anova

# Compute the analysis of variance
res.H.aov <- lm(Ho ~ Region + variable + Region:variable, data = Hstat_Gen_ZONE.comp)
# Summary of the analysis
summary(res.H.aov)
anova(res.H.aov)

# Tukey tests
compute.H.stat <- function(df, snps){
  
  res.aov <- aov(He ~ Region, data = df %>% filter(variable == snps))
  tukey <- TukeyHSD(res.aov)
  
  # compact letter display
  cld <- multcompLetters4(res.aov, tukey)
  cld <- as.data.frame.list(cld$Region) 
  
  res <- cld %>% dplyr::select(Letters) %>% 
    mutate(Region = row.names(cld),
           variable = snps)
  return(res)
}

compute.H.stat(Hstat_Gen_ZONE.comp, "Full")


stat.H.letters <- bind_rows(compute.H.stat( Hstat_Gen_ZONE.comp, "Full"),
                            compute.H.stat( Hstat_Gen_ZONE.comp, "Neutral"),
                            compute.H.stat(Hstat_Gen_ZONE.comp, "Outlier")) 

stat.H.letters <- stat.H.letters %>% left_join(Hstat_Gen_ZONE.comp %>% group_by(Region, variable) %>% 
                                                 mutate(maxH = max(He) ) %>% 
                                                 dplyr::select(maxH) %>% 
                                                 distinct()) %>% 
  mutate (Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
          
          variable = factor(variable, levels = c("Neutral", "Outlier", "Full"))) 


stat.H.letters

gg.Hstat.full <- Hstat_Gen_ZONE.comp  %>%
 filter(variable == "Full") %>% 
  ggplot(aes(x = Region, y = He, col = Region)) +
  geom_jitter(height = 0, alpha = 0.5) +
  geom_boxplot(col = "gray10", fill = NA) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  scale_x_discrete(labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+
  labs(y = "He", x = "") +
  # Add space above my axis
  geom_blank(data = stat.H.letters  %>%  filter(variable == "Full"), aes(y = maxH + (maxH *0.05)))  +
  #  Add stats
  geom_text(data = stat.H.letters %>%  filter(variable == "Full") , aes(label = Letters, x = Region, y = maxH), col = "darkgray",  vjust = -1.5) +
  facet_wrap(~"Complete SNP panel") + 
  labs(y= "Expected heterozygosity") +
  #  facet_grid(variable ~ ., scale = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        #axis.ticks = element_blank(),
        #axis.line = element_line(colour = "grey50"),
        #panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        strip.background = element_rect(fill = "white")
        #plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )
gg.Hstat.full


gg.Hstat.outlier <- Hstat_Gen_ZONE.comp  %>%
  filter(variable == "Outlier") %>% 
  ggplot(aes(x = Region, y = He, col = Region)) +
  geom_jitter(height = 0, alpha = 0.5) +
  geom_boxplot(col = "gray10", fill = NA) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  scale_x_discrete(labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+
  labs(y = "He", x = "") +
  # Add space above my axis
  geom_blank(data = stat.H.letters  %>%  filter(variable == "Outlier"), aes(y = maxH + (maxH *0.05)))  +
  #  Add stats
  geom_text(data = stat.H.letters %>%  filter(variable == "Outlier") , aes(label = Letters, x = Region, y = maxH), col = "darkgray",  vjust = -1.5) +
  facet_wrap(~"Outlier SNP panel") + 
  labs(y= "Expected heterozygosity") +
  #  facet_grid(variable ~ ., scale = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        #axis.ticks = element_blank(),
        #axis.line = element_line(colour = "grey50"),
        #panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        strip.background = element_rect(fill = "white")
        #plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

gg.Hstat.outlier


gg.Hstat.neutral <- Hstat_Gen_ZONE.comp  %>%
  filter(variable == "Neutral") %>% 
  ggplot(aes(x = Region, y = He, col = Region)) +

  geom_jitter(height = 0, alpha = 0.5) +
  geom_boxplot(col = "gray10", fill = NA) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  scale_x_discrete(labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+
  labs(y = "He", x = "") +
  # Add space above my axis
  geom_blank(data = stat.H.letters  %>%  filter(variable == "Neutral"), aes(y = maxH + (maxH *0.05)))  +
  #  Add stats
  geom_text(data = stat.H.letters %>%  filter(variable == "Neutral") , aes(label = Letters, x = Region, y = maxH), col = "darkgray",  vjust = -1.5) +
  facet_wrap(~"Neutral SNP panel") + 
  labs(y= "Expected heterozygosity") +
  #  facet_grid(variable ~ ., scale = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        #axis.ticks = element_blank(),
        #axis.line = element_line(colour = "grey50"),
        #panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        strip.background = element_rect(fill = "white")
        #plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )
#, 
#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gg.Hstat.neutral


# Combined Pop struct figures ---------------------------------------------

gg.popstruc.1 <- ggpubr::ggarrange(#gPCA.final.14331snps + theme(legend.position = "none"),
                              gPCA.final.Neutralsnps+ theme(legend.position = "none") +
                                scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                                                   labels = c("ARC", "NFL", "GSL", "SS", "FLC")),
                              gPCA.final.Outliersnps+ theme(legend.position = "none"),
                              #gg.pca.all + theme(legend.position = "none"),
                              gg.pca.neutral + theme(legend.position = "none"), 
                              gg.pca.outlier + theme(legend.position = "none"),
                              #gg.Hstat.full,
                              gg.Hstat.neutral,
                              gg.Hstat.outlier +   labs(y= ""),
                              nrow = 3, ncol = 2, common.legend = T, align = "hv" ,
                              labels = LETTERS
                              
)
gg.popstruc.1

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "PopStruct_NeutralvsOutliers_20231006.png"), 
       plot = gg.popstruc.1,
       height = 8, width = 6, units = "in", bg = "white")   


gg.popstruc.2 <- ggpubr::ggarrange(gPCA.final.14331snps + theme(legend.position = "none") +
                                     scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                                                        labels = c("ARC", "NFL", "GSL", "SS", "FLC")),
  #gPCA.final.Neutralsnps+ theme(legend.position = "none"),
  #gPCA.final.Outliersnps+ theme(legend.position = "none"),
  gg.pca.all + theme(legend.position = "none"),
  #gg.pca.neutral + theme(legend.position = "none"), 
  #gg.pca.outlier + theme(legend.position = "none"),
  gg.Hstat.full,
  #gg.Hstat.neutral,
  #gg.Hstat.outlier,
  nrow = 1, ncol = 3, common.legend = T, align = "hv" ,
  labels = LETTERS
  
)
gg.popstruc.2

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "PopStruct_All_20231006.png"), 
       plot = gg.popstruc.2,
       height = 3.5, width = 8, units = "in", bg = "white")   

plot(table.fst(final.Neutralsnps_area.FST)$Fst, table.fst(final.Outliersnps_area.FST)$Fst)

table.fst(final.Neutralsnps_area.FST)$Fst %>% hist()

cor.test(table.fst(final.Neutralsnps_area.FST)$Fst, table.fst(final.Outliersnps_area.FST)$Fst,
         method = c("spearman"))

mean(table.fst(final.Outliersnps_area.FST)$Fst / table.fst(final.Neutralsnps_area.FST)$Fst)
mean(table.fst(final.Outliersnps_area.FST)$Fst)
mean(table.fst(final.Neutralsnps_area.FST)$Fst)


# RDA ---------------------------------------------------------------------

# Compile variables

env.summary.scale <- bind_cols(Env.data[,1], data.frame(scale(Env.data[,-1], center=TRUE, scale=TRUE)))
names(env.summary.scale)[1] <- "Gen_ZONE_FG"
env.summary.scale

## Standardization of the variables
## Recovering scaling coefficients
scale_env <- attr(scale(Env.data[,-1], center=TRUE, scale=TRUE), 'scaled:scale')
center_env <- attr(scale(Env.data[,-1], center=TRUE, scale=TRUE), 'scaled:center')

sapply(Env.data[,-1], min) 
sapply(Env.data[,-1], max) 

Variables <- data.frame(Gen_ZONE_FG = rownames(freq.MAF.final.Gen_ZONE_FG)) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE_FG, Lat, Long) %>% distinct()) %>% 
  left_join(env.summary.scale %>% dplyr::select(-c(Lat, Long))) %>% 
  left_join(PopStruct)

Variables %>% head()


## Environmental model only to check if my choice are significant

## Null model
RDA0 <- rda(freq.MAF.final.Gen_ZONE_FG ~ 1,  Variables) 

## Full model
#RDAfull <- rda(freq.MAF.final.Gen_ZONE_FG ~  T_0m_larv + S_0m_yearly + T_Bott_Winter, Variables)
RDAfull <- rda(freq.MAF.final.Gen_ZONE_FG ~  Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval, Variables)


#Explained variance by the global db-RDA model. Adjusted R2 accounts for the number of variables 
RsquareAdj(RDAfull)

#db-RDA Global model probability
anova(RDAfull, perm=999)
#p=0.001 

RDAfull$CCA$eig

screeplot(RDAfull )
plot(RDAfull , scaling = 3)


vegan::vif.cca(RDAfull)

mod <-  ordiR2step(RDA0, Pin = 0.01, scope = RDAfull)                

#Pin = 0.01, R2permutations = 1000, R2scope = T
#summary table with selected variables               
mod$anova

RsquareAdj(mod)

# Check MEM instead of lat / long


#Visualising the links longer than the treashold
#s.label(Coor.sites, nb = attr(dbmem.sites, "listw"))

#Visualising the 16 mem. Can create more spatial structure type when MEM are selected together in the same analysis
#The 1rst MEMs are large spatial scales, the last MEMs are small spatial scales
#s.value(Coor.sites, dbmem.sites[,1:4])
##MEM can be used in stats like any other explanatory variables

RDA0 <- rda(freq.MAF.final.Gen_ZONE_FG ~ 1,  MEM.df) 

## Full model
RDAfull <- rda(freq.MAF.final.Gen_ZONE_FG ~ . , MEM.df )

mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
mod

# MEM1 + MEM3

## Full model for variance partitioning
pRDAfull <- rda(freq.MAF.final.Gen_ZONE_FG ~ PC1 + PC2 + PC3 + Long + Lat + Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval,  Variables)
RsquareAdj(pRDAfull)
#anova(pRDAfull)

## Pure climate model
pRDAclim <- rda(freq.MAF.final.Gen_ZONE_FG ~  Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval + Condition(Long + Lat + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAclim)
#anova(pRDAclim)
pRDAclim$anova

## Pure neutral population structure model  
pRDAstruct <- rda(freq.MAF.final.Gen_ZONE_FG ~ PC1 + PC2 + PC3 + Condition(Long + Lat +   Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval),  Variables)
RsquareAdj(pRDAstruct)
#anova(pRDAstruct)

##Pure geography model 
pRDAgeog <- rda(freq.MAF.final.Gen_ZONE_FG ~ Long + Lat + Condition( Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval + PC1 + PC2 + PC3 ),  Variables)
RsquareAdj(pRDAgeog)

anova(pRDAgeog)


rda.summary <- data.frame(models = c("full", "env", "structure", "geography", "confounded", "unexplained", "total"))

rda.summary$inertia <- c(pRDAfull$CCA$tot.chi, pRDAclim$CCA$tot.chi, pRDAstruct$CCA$tot.chi, pRDAgeog$CCA$tot.chi,
                         pRDAfull$CCA$tot.chi -  pRDAclim$CCA$tot.chi - pRDAstruct$CCA$tot.chi - pRDAgeog$CCA$tot.chi,
                         pRDAclim$tot.chi - pRDAfull$CCA$tot.chi, pRDAclim$tot.chi
                         )

rda.summary$r2 <- c(RsquareAdj(pRDAfull)$adj.r.squared, RsquareAdj(pRDAclim)$adj.r.squared, RsquareAdj(pRDAstruct)$adj.r.squared, RsquareAdj(pRDAgeog)$adj.r.squared, NA, NA, NA)

rda.summary$p <- c(anova(pRDAfull)$`Pr(>F)`[1], anova(pRDAclim)$`Pr(>F)`[1], anova(pRDAstruct)$`Pr(>F)`[1], anova(pRDAgeog)$`Pr(>F)`[1], NA, NA, NA)

rda.summary$PropVariance <- c(1, pRDAclim$CCA$tot.chi/pRDAfull$CCA$tot.chi, pRDAstruct$CCA$tot.chi/pRDAfull$CCA$tot.chi, pRDAgeog$CCA$tot.chi/pRDAfull$CCA$tot.chi, 
                              (pRDAfull$CCA$tot.chi -  pRDAclim$CCA$tot.chi - pRDAstruct$CCA$tot.chi - pRDAgeog$CCA$tot.chi) / pRDAfull$CCA$tot.chi, NA, NA
                              )

rda.summary$PropVarianceTotal <- c(pRDAfull$CCA$tot.chi / pRDAclim$tot.chi , pRDAclim$CCA$tot.chi/pRDAfull$tot.chi, pRDAstruct$CCA$tot.chi/pRDAfull$tot.chi, pRDAgeog$CCA$tot.chi/pRDAfull$tot.chi, 
                                   (pRDAfull$CCA$tot.chi -  pRDAclim$CCA$tot.chi - pRDAstruct$CCA$tot.chi - pRDAgeog$CCA$tot.chi) / pRDAfull$tot.chi,
                                   (pRDAclim$tot.chi - pRDAfull$CCA$tot.chi) / pRDAclim$tot.ch, 1
)

rda.summary
write_csv(rda.summary, file = file.path(here::here(), "02_Results/01_Overall_PopGen/", "08_RDA", "PartialRDAsummary.csv"))

# RDA - ENV one by one ----------------------------------------------------

## Full model for variance partitioning
pRDAfull <- rda(freq.MAF.final.Gen_ZONE_FG ~ PC1 + PC2 + PC3 + Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval,  Variables)
RsquareAdj(pRDAfull)

pRDAclim <- rda(freq.MAF.final.Gen_ZONE_FG ~  Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval + Condition(PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAclim)
pRDAclim.anova <- anova(pRDAclim)

## Pure climate model
pRDAclim.Tbtm.Winter <- rda(freq.MAF.final.Gen_ZONE_FG ~  Tbtm.Winter  + Condition(PC1 + PC2 + PC3+ Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval),  Variables)
RsquareAdj(pRDAclim.Tbtm.Winter)
pRDAclim.Tbtm.Winter.anova <- anova(pRDAclim.Tbtm.Winter)

pRDAclim.Tbtm.Summer <- rda(freq.MAF.final.Gen_ZONE_FG ~  Tbtm.Summer  + Condition(PC1 + PC2 + PC3+ Tbtm.Winter + Sbtm.ann +  SSS.ann + SST.Larval),  Variables)
RsquareAdj(pRDAclim.Tbtm.Summer)
pRDAclim.Tbtm.Summer.anova <- anova(pRDAclim.Tbtm.Summer)

pRDAclim.Sbtm.ann <- rda(freq.MAF.final.Gen_ZONE_FG ~  Sbtm.ann  + Condition(PC1 + PC2 + PC3+ Tbtm.Summer + Tbtm.Summer  +  SSS.ann + SST.Larval),  Variables)
RsquareAdj(pRDAclim.Sbtm.ann)
pRDAclim.Sbtm.ann.anova <- anova(pRDAclim.Sbtm.ann)

pRDAclim.SSS.ann <- rda(freq.MAF.final.Gen_ZONE_FG ~ SSS.ann  + Condition(PC1 + PC2 + PC3+ Tbtm.Summer + Tbtm.Summer  +   Sbtm.ann +  SST.Larval),  Variables)
RsquareAdj(pRDAclim.SSS.ann)
pRDAclim.SSS.ann.anova <- anova(pRDAclim.SSS.ann)

pRDAclim.SST.Larval <- rda(freq.MAF.final.Gen_ZONE_FG ~ SST.Larval  + Condition(PC1 + PC2 + PC3+ Tbtm.Summer + Tbtm.Summer  +   Sbtm.ann + SSS.ann),  Variables)
RsquareAdj(pRDAclim.SST.Larval)
pRDAclim.SST.Larval.anova <- anova(pRDAclim.SST.Larval)



rda.env.summary <- data.frame(models = c("env", "Tbtm.Winter", "Tbtm.Summer", "Sbtm.ann", "SSS.ann", "SST.Larval",  "confounded",  "total"))

rda.env.summary$inertia <- c(pRDAclim$CCA$tot.chi, 
                             pRDAclim.Tbtm.Winter$CCA$tot.chi, pRDAclim.Tbtm.Summer$CCA$tot.chi, pRDAclim.Sbtm.ann$CCA$tot.chi, pRDAclim.SSS.ann$CCA$tot.chi, pRDAclim.SST.Larval$CCA$tot.chi,
                        pRDAclim$CCA$tot.chi -  pRDAclim.Tbtm.Winter$CCA$tot.chi - pRDAclim.Tbtm.Summer$CCA$tot.chi - pRDAclim.Sbtm.ann$CCA$tot.chi - pRDAclim.SSS.ann$CCA$tot.chi - pRDAclim.SST.Larval$CCA$tot.chi,
                       pRDAclim$tot.chi
                        
)

rda.env.summary$r2 <- c( RsquareAdj(pRDAclim)$adj.r.squared, RsquareAdj(pRDAclim.Tbtm.Winter)$adj.r.squared, RsquareAdj(pRDAclim.Tbtm.Summer)$adj.r.squared, RsquareAdj(pRDAclim.Sbtm.ann)$adj.r.squared, RsquareAdj(pRDAclim.SSS.ann)$adj.r.squared,RsquareAdj(pRDAclim.SST.Larval)$adj.r.squared, NA, NA)

rda.env.summary$p <- c(pRDAclim.anova$`Pr(>F)`[1], pRDAclim.Tbtm.Winter.anova$`Pr(>F)`[1], pRDAclim.Tbtm.Summer.anova$`Pr(>F)`[1], pRDAclim.Sbtm.ann.anova$`Pr(>F)`[1],  pRDAclim.SSS.ann.anova$`Pr(>F)`[1],  pRDAclim.SST.Larval.anova$`Pr(>F)`[1], NA, NA)

rda.env.summary$PropVariance <- c(1 , pRDAclim.Tbtm.Winter$CCA$tot.chi/pRDAclim$CCA$tot.chi,  pRDAclim.Tbtm.Summer$CCA$tot.chi/pRDAclim$CCA$tot.chi,  pRDAclim.Sbtm.ann$CCA$tot.chi/pRDAclim$CCA$tot.chi, pRDAclim.SSS.ann$CCA$tot.chi/pRDAclim$CCA$tot.chi, pRDAclim.SST.Larval$CCA$tot.chi/pRDAclim$CCA$tot.chi, 
                              (pRDAclim$CCA$tot.chi -  pRDAclim.Tbtm.Winter$CCA$tot.chi - pRDAclim.Tbtm.Summer$CCA$tot.chi - pRDAclim.Sbtm.ann$CCA$tot.chi- pRDAclim.SSS.ann$CCA$tot.chi - pRDAclim.SST.Larval$CCA$tot.chi  ) / pRDAclim$CCA$tot.chi, NA
)

rda.env.summary$PropVarianceTotal <- c(pRDAclim$CCA$tot.chi / pRDAclim$tot.chi , pRDAclim.Tbtm.Winter$CCA$tot.chi/pRDAclim$tot.chi,  pRDAclim.Tbtm.Summer$CCA$tot.chi/pRDAclim$tot.chi,  pRDAclim.Sbtm.ann$CCA$tot.chi/pRDAclim$tot.chi,pRDAclim.SSS.ann$CCA$tot.chi/pRDAclim$tot.chi, pRDAclim.SST.Larval$CCA$tot.chi/pRDAclim$tot.chi,     
                                       (pRDAclim$CCA$tot.chi -  pRDAclim.Tbtm.Winter$CCA$tot.chi - pRDAclim.Tbtm.Summer$CCA$tot.chi - pRDAclim.Sbtm.ann$CCA$tot.chi- pRDAclim.SSS.ann$CCA$tot.chi - pRDAclim.SST.Larval$CCA$tot.chi  ) / pRDAclim$tot.chi, 1
)

write_csv(rda.env.summary, file = file.path(here::here(), "02_Results/01_Overall_PopGen/", "08_RDA", "PartialRDAsummary_ENVonly.csv"))


# Genotype-environement association : identifying loci under selec --------

RDA_env.constrained   <- rda(freq.MAF.final.Gen_ZONE_FG ~   Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval + Condition(PC1 + PC2 + PC3),  Variables)


plot(RDA_env.constrained, scaling = 2)

screeplot(RDA_env.constrained, main="Eigenvalues of constrained axes")



source("01_Scripts/Functions/radapt.R")

rdadapt_env.constrained   <-rdadapt(RDA_env.constrained, 2)

## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env.constrained$p.values)

## Identifying the loci that are below the p-value threshold
outliers.constrained <- data.frame(Loci = colnames(freq.MAF.final.Gen_ZONE_FG)[which(rdadapt_env.constrained$p.values<thres_env)], 
                          p.value = rdadapt_env.constrained$p.values[which(rdadapt_env.constrained$p.values<thres_env)])
outliers.constrained    
outliers.constrained  %>% nrow()

## List of outlier names
outliers_rdadapt_env <- as.character(outliers.constrained$Loci)

RDA_env <- RDA_env.constrained

screeplot(RDA_env, main="Eigenvalues of constrained axes")

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores) %>% 
              mutate(type = ifelse(names %in% outliers_rdadapt_env, "Outlier", "Neutral"))

TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) %>% 
  mutate(ID = row.names(scores(RDA_env, choices=c(1,2), display="bp")) ) %>% left_join(final.env.names)
TAB_var  # pull the biplot scores

## Biplot of RDA loci and variables scores
gbiplot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
 # geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5) +
 ggrepel::geom_label_repel(data = TAB_var, aes(label = new.ID, x=1.1*RDA1, y=1.1*RDA2), size = 2, max.overlaps = 20
  ) +
  
  labs(x= paste0("RDA 1 (", round(eigenvals(RDA_env)[1]/sum(eigenvals(RDA_env)[1:5]),3) * 100 , "%)"), 
       y= paste0("RDA 2 (", round(eigenvals(RDA_env)[2]/sum(eigenvals(RDA_env)[1:5]),3) * 100 , "%)")) +
  
         facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.8)),
        strip.background = element_rect(fill = "white"))
gbiplot


# Enriched RDA ------------------------------------------------------------

RDA_outliers <- rda(freq.MAF.final.Gen_ZONE_FG[, outliers_rdadapt_env] ~ Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval ,  Variables)
screeplot(RDA_outliers, main="Eigenvalues of constrained axes")

round(eigenvals(RDA_outliers)[1]/sum(eigenvals(RDA_outliers)[1:5]),2)
eigenvals(RDA_outliers)[2]/sum(eigenvals(RDA_outliers))

## RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp")) %>% 
  mutate(ID = row.names(scores(RDA_outliers, choices=c(1,2), display="bp")) ) %>% left_join(final.env.names)
TAB_var  # pull the biplot scores


RsquareAdj(RDA_outliers)$adj.r.squared
anova.cca(RDA_outliers, step=9999, by="axis")

gbiplot.enr <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*4, y=RDA2*4), colour = "#F9A242FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
#  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5) +
  ggrepel::geom_label_repel(data = TAB_var, aes(label = new.ID, x=1.1*RDA1, y=1.1*RDA2), size = 2, force_pull = 100, max.overlaps = 20
  ) +
  # geom_path(data = data.frame(x1= seq(-0.9,0.9,0.15), y2 = rep(-0.9,13)), 
  #                           aes(x = x1, y=y2, colour = x1),
  #                            size = 1, arrow=arrow(length = unit(0.02, "npc")))+
  # geom_path(data = data.frame(x1= seq(-0.9,0.9,0.2), y2 = rep(-0.9,10)), 
  #           aes(y = x1, x=y2, colour = x1),
  #           size = 1, arrow=arrow(length = unit(0.02, "npc")))+
  # #geom_segment(aes(xend=0.9, yend=-1, x=-0.8, y=-0.9), colour=after_stat(value), size=2, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  # scale_color_viridis_c(option = "A", direction = -1) +
  #   
 labs(x= paste0("RDA 1 (", round(eigenvals(RDA_outliers)[1]/sum(eigenvals(RDA_outliers)[1:5]),3) * 100 , "%)"), 
      y= paste0("RDA 2 (", round(eigenvals(RDA_outliers)[2]/sum(eigenvals(RDA_outliers)[1:5]),3) * 100 , "%)")

 ) +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        strip.background = element_rect(fill = "white"))
gbiplot.enr


# Adaptative index - Constrained ------------------------------------------


## Loading the climatic rasters

newproj  = "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
  
## Species range shapefile (download information at beginning of tutorial)
range <- readOGR("./00_Data/99_SIG/range2.shp") 
crs(range) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras_current,  range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj_future <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras_future,  range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

admin <- ne_countries(scale = "medium", returnclass = "sf")

## Vectorization of the climatic rasters for ggplot
RDA_proj_current <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj_current <- lapply(RDA_proj_current, function(x) projectRaster(from = x, crs =  newproj))
RDA_proj_current <- lapply(RDA_proj_current, function(x) rasterToPoints(x))

RDA_proj_future <- list(res_RDA_proj_future$RDA1, res_RDA_proj_future$RDA2)
RDA_proj_future <- lapply(RDA_proj_future, function(x) projectRaster(from = x, crs =  newproj))
RDA_proj_future <- lapply(RDA_proj_future, function(x) rasterToPoints(x))


## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA_current <- as.data.frame(do.call(rbind, RDA_proj_current[1:2]))
colnames(TAB_RDA_current)[3] <- "value"
TAB_RDA_current$variable <- factor(c(rep("RDA1", nrow(RDA_proj_current[[1]])), rep("RDA2", nrow(RDA_proj_current[[2]]))), levels = c("RDA1","RDA2"))

TAB_RDA_current %>% head()# <- TAB_RDA %>% mutate(value.cor = )

pop.rda <- sf::st_as_sf(Variables, coords = c("Long", "Lat"), crs = 4326)

TAB_RDA_future <- as.data.frame(do.call(rbind, RDA_proj_future[1:2]))
colnames(TAB_RDA_future)[3] <- "value"
TAB_RDA_future$variable <- factor(c(rep("RDA1", nrow(RDA_proj_future[[1]])), rep("RDA2", nrow(RDA_proj_future[[2]]))), levels = c("RDA1","RDA2"))
TAB_RDA_future$model <- "Future"
TAB_RDA_current$model <- "Present"


# SFA

SFA.shp <- terra::vect("00_Data/99_SIG/SFA.shp")
Region.shp <-  terra::vect("00_Data/99_SIG/Regions_MS_corrected.shp")


#SFA.shp <- rgdal::readOGR("00_Data/01_SIG/SFA.shp")

# reproject data
SFA2.shp <- spTransform(SFA.shp,
                        "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

plot(SFA.shp,
     main = "Shapefile imported into R - SFA",
     axes = TRUE,
     border = "blue")

SFA.shp

SFA_df <- fortify(SFA2.shp)

bind_rows(TAB_RDA_current, TAB_RDA_future) %>%  ggplot(aes(x = value)) +
  geom_density(aes(col = model))+
  facet_grid(~variable)


future.RDA1 <- bind_rows(TAB_RDA_current, TAB_RDA_future) %>% filter(value < 8,
                                   variable == "RDA1") %>%   
  mutate(model = factor(model, levels = c("Present", "Future"))) %>% 
  ggplot() + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  #geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(-3, 4, length.out=10), include.lowest = T))) + 
  scale_fill_binned(n.breaks=12, type ="viridis", alpha = 0.8, direction = -1, option = "A") +
  #scale_fill_viridis_b(alpha = 0.8, direction = -1, option = "A", nbreaks = 10) +
  #scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  geom_sf(data = admin, fill=NA, size=0.1) +  
  geom_sf(data = pop.rda, cex = 0.5) +

#  # SFA
#  geom_path(data = SFA_df, 
#            aes(x = long, y = lat, group = group),
#            color = "gray25", size = .2) +
  
    
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +
  
  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "black",
          lwd = 0.5, fill = NA) +  
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  
  geom_sf_label(data = Region.shp, aes(label = layer),
               # col = rep(c("black","magenta", "darkorange","blue","red"),2),
                fill = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-300000),
                nudge_y = c(-900000,650000,100000,800000,220000)) +
  
  #coord_sf(xlim = c(-70, -43), ylim = c(43, 68), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA1")) +
  facet_wrap(~model) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        strip.background = element_rect(fill = "white"))
future.RDA1

splitFacet <- function(x){
  facet_vars <- names(x$facet$params$facets)         # 1
  x$facet    <- ggplot2::ggplot()$facet              # 2
  datasets   <- split(x$data, x$data[facet_vars])    # 3
  new_plots  <- lapply(datasets,function(new_data) { # 4
    x$data <- new_data
    x})
}    

new_future.RDA1 <- splitFacet(future.RDA1)
length(new_future.RDA1) # [1] 12
new_future.RDA1[[2]] 


future.RDA1.density <-  bind_rows(TAB_RDA_current, TAB_RDA_future) %>% filter(value < 8,
                                                                             variable == "RDA1") %>% 
  mutate(model = factor(model, levels = c("Future", "Present"))) %>% 
  ggplot(aes(x = value, y = model  ))+
 # ggridges::geom_density_ridges(scale = 100) +
  ggridges::geom_density_ridges_gradient(scale = 10000, aes(fill =stat(x), lty = model) )+# geom_density() +
   scale_fill_viridis_b(name = "RDA1", alpha = 0.5, direction = -1, option = "A", n.breaks=12, guide = "none") +
  scale_linetype_manual(values = c("solid", "dotted"), limits = c("Present", "Future"), name = "")+
  #geom_vline(xintercept = quantile(TAB_RDA_future %>% dplyr::filter( variable == "RDA1") %>% pull(value), c(0.025, 0.975) ),
  #           col = "gray10", lty = 5) +
  facet_wrap(~"Distribution of RDA1 scores") +
  labs(x=  "RDA1", y = "Density") +
   theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(5.5, 25, 5.5, 28), "points"),
        strip.background = element_rect(fill = "white")
        #legend.position = "none"
        )

future.RDA1.density

future.RDA2 <- bind_rows(TAB_RDA_current, TAB_RDA_future) %>% filter(value < 5,
                                                             variable == "RDA2") %>%   
  mutate(model = factor(model, levels = c("Present", "Future"))) %>% 
  ggplot() + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  #geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(-3, 4, length.out=10), include.lowest = T))) + 
  scale_fill_binned(n.breaks=9, type ="viridis", alpha = 0.8, direction = -1, option = "A") +
  #scale_fill_viridis_b(alpha = 0.8, direction = -1, option = "A", nbreaks = 10) +
  #scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  geom_sf(data = admin, fill=NA, size=0.1) +  
  geom_sf(data = pop.rda, cex = 0.5) +
  
  #  # SFA
  #  geom_path(data = SFA_df, 
  #            aes(x = long, y = lat, group = group),
  #            color = "gray25", size = .2) +
  
  
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +
  
  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "black",
          lwd = 0.5, fill = NA) +  

    geom_sf_label(data = Region.shp, aes(label = layer),
               # col = rep(c("black","magenta", "darkorange","blue","red"),2),
                fill = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-300000),
                nudge_y = c(-900000,650000,100000,800000,220000)) +
  
  
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA2")) +
  facet_wrap(~model) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        strip.background = element_rect(fill = "white"))
future.RDA2

new_future.RDA2 <- splitFacet(future.RDA2)
length(new_future.RDA2) # [1] 12
new_future.RDA2[[1]] 

future.RDA2.density <-  bind_rows(TAB_RDA_current, TAB_RDA_future) %>% filter(value < 5,
                                                                              variable == "RDA2") %>% 
  mutate(model = factor(model, levels = c("Future", "Present"))) %>% 
  ggplot(aes(x = value, y = model  ))+
  # ggridges::geom_density_ridges(scale = 100) +
  ggridges::geom_density_ridges_gradient(scale = 10000, aes(fill =stat(x), lty = model) )+# geom_density() +
  scale_fill_viridis_b(name = "RDA1", alpha = 0.5, direction = -1, option = "A", n.breaks=9, guide = "none") +
  scale_linetype_manual(values = c("solid", "dotted"), limits = c("Present", "Future"), name = "")+
  
  #geom_vline(xintercept = quantile(TAB_RDA_future %>% dplyr::filter( variable == "RDA2") %>% pull(value), c(0.025, 0.975) ),
  #           col = "gray10", lty = 5) +
  facet_wrap(~"Distribution of RDA2 scores") +
  labs(x=  "RDA2", y = "Density") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(5.5, 25, 5.5, 28), "points"),
        strip.background = element_rect(fill = "white")
        #legend.position = "none"
  )

future.RDA2.density

venn.out <-  ggVennDiagram(list( putative.BayesOutliers.final, as.character(outliers_rdadapt_env.cor),  putative.PcAdaptOutliers.final), 
                           label_alpha = 0,
                           set_size = 3,
                           label_size = 3,
                           label = "count",
                           category.names = c("BayeScan","pRDA",  "pcadapt")) +
  ggplot2::scale_fill_gradient(name = "N snps", low="white",high = "pink", limits = c(0,650)) +
  ggplot2::scale_color_manual(values = rep("darkgray", 10)) +
  scale_x_continuous(limits = c(-100, 1100)) +
  facet_wrap(~"Outliers vs adaptive SNPs") + theme_bw(base_size = 11) +
  theme(axis.ticks = element_blank(),
        #axis.line = element_blanline(colour = "grey50"),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(),
        axis.text =  element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")#,
#panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
#plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
venn.out

gg.BIG.adapt2 <- ggarrange(ggRDA <- ggarrange( gbiplot + theme(legend.position = "none"),  
                                               venn.out, 
                                               gbiplot.enr ,
                                              nrow = 3, ncol = 1, 
                                              align = "hv", labels = LETTERS),
                          
                                    ggarrange(new_future.RDA1[[1]] + facet_wrap(~"RDA1"), 
                                              new_future.RDA2[[1]] + facet_wrap(~"RDA2"), 
                                              nrow = 2, ncol = 1, heights = c(1,1), 
                                              align = "hv", labels = LETTERS[4]),
                                    

                          nrow = 1, ncol = 2, widths = c(2,3)
)

gg.BIG.adapt2

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "Adaptative.landscape_BIG3.png"), 
       plot =gg.BIG.adapt2, bg = "white",
       height = 8, width = 6.5, units = "in") 


# SVG ---------------------------------------------------------------------

FreqData <- freq.MAF.final.Gen_ZONE_FG
dimnames(FreqData)[[2]] <- str_remove(dimnames(FreqData)[[2]], "[.][:digit:]")
  
outliers_rdadapt_env.cor <- str_remove(outliers_rdadapt_env, "[.][:digit:]")

length( putative.outliers.final.high %>% unique())

adapt.stat <- function(FreqData, snps){

  mean_q <- apply(1-FreqData[,snps], 2, mean)

# Population adaptive index (Bonin et al. 2007) 
PAI <- lapply(1:nrow(FreqData), function(x) sum(abs((1-FreqData[x, snps])-mean_q), na.rm=T))
names(PAI) <- row.names(FreqData)

SVG <- lapply(lapply(1:nrow(FreqData), function(x) FreqData[x, snps]*(1-FreqData[x, snps])), function(y) mean(unlist(y)))
names(SVG) <- row.names(FreqData)


TAB <- data.frame(ID = names(SVG), PAI = unlist(PAI), SVG = unlist(SVG)) 

return(TAB)
  
}

table(putative.neutral.final %in% dimnames(FreqData)[[2]])

snps<- putative.neutral.final

SVG.res <- bind_rows(
adapt.stat(FreqData, outliers_rdadapt_env.cor) %>% mutate(variable = "Adaptative"),
adapt.stat(FreqData, putative.neutral.final) %>% mutate(variable = "Neutral"),
adapt.stat(FreqData, putative.outliers.final) %>% mutate(variable = "Outlier"),
adapt.stat(FreqData,locNames(gl.final)) %>% mutate(variable = "Full")
) %>% left_join(pop.data %>% dplyr::select(ID = Gen_ZONE_FG, Gen_ZONE, Region = RegionAssesment) %>% distinct(.keep_all = T))  %>% 
  mutate(Region = ifelse(Gen_ZONE %in% c("SFA-4", "SFA-5", "SFA-5", "SFA-7"), "NFL", 
                         ifelse(Gen_ZONE == "NAFO-3M", "3M",  Region)),
         Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
         
         variable = factor(variable, levels = c("Neutral", "Outlier", "Adaptative", "Full"))) 


# One-way anova

# Compute the analysis of variance
res.aov <- lm(SVG ~ Region + variable + Region:variable, data = SVG.res)
# Summary of the analysis
summary(res.aov)
anova(res.aov)

compute.SVG.stat <- function(df, snps){
  
  res.aov <- aov(SVG ~ Region, data = df %>% filter(variable == snps))
  tukey <- TukeyHSD(res.aov)
  
  # compact letter display
  cld <- multcompLetters4(res.aov, tukey)
  cld <- as.data.frame.list(cld$Region) 
  
  res <- cld %>% dplyr::select(Letters) %>% 
         mutate(Region = row.names(cld),
                variable = snps)
  return(res)
}

stat.letters <- bind_rows(compute.SVG.stat( SVG.res, "Adaptative"),
                          compute.SVG.stat( SVG.res, "Neutral"),
                          compute.SVG.stat( SVG.res, "Outlier")) 

stat.letters <- stat.letters %>% left_join(SVG.res %>% group_by(Region, variable) %>% 
                                           mutate(maxSVG = max(SVG) ) %>% 
                                          dplyr::select(maxSVG) %>% 
                                            distinct()) %>% 
  mutate (Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
         
         variable = factor(variable, levels = c("Neutral", "Outlier", "Adaptative"))) 


stat.letters


gg.SVG.adapt <- SVG.res %>%
  dplyr::filter(variable =="Adaptative") %>% # left_join(stat.letters) %>% 
  #  mutate(Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
  
  #        variable = factor(variable, levels = c("Neutral", "Outlier", "Adaptative"))) %>% 
  ggplot(aes(x = Region, y = SVG, col = Region)) +
  geom_boxplot(col = "darkgray") +
  geom_jitter(height = 0, alpha = 0.5) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  scale_x_discrete(labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+
  labs(y = "aSGV", x = "") +
  # Add space above my axis
  geom_blank(data = stat.letters %>%
               dplyr::filter(variable =="Adaptative"), aes(y = maxSVG + (maxSVG *0.05)))  +
  #  Add stats
  geom_text(data = stat.letters %>%
              dplyr::filter(variable =="Adaptative"), aes(label = Letters, x = Region, y = maxSVG), col = "darkgray",  vjust = -1.5) +
  facet_wrap(~"aSGV") +
  #facet_grid(variable ~ ., scale = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        #axis.ticks = element_blank(),
        #axis.line = element_line(colour = "grey50"),
        #panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_line(linetype = "dashed"),
        strip.background = element_rect(fill = "white")
  )
#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gg.SVG.adapt


##

#basic_stat_Region <- hierfstat::basic.stats(gi.final.Region, diploid = TRUE)
basic_stat_Gen_ZONE_FG.neutral <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = putative.neutral.final], diploid = TRUE)
basic_stat_Gen_ZONE_FG.outlier <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = putative.outliers.final], diploid = TRUE)
basic_stat_Gen_ZONE_FG.adaptative <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG[loc = outliers_rdadapt_env.cor], diploid = TRUE)

#basic_stat_Gen_ZONE_FG <- hierfstat::basic.stats(gi.final.Gen_ZONE_FG, diploid = TRUE)

#save(list = c("basic_stat_Gen_ZONE_FG.neutral", "basic_stat_Gen_ZONE_FG.outlier", "basic_stat_Gen_ZONE_FG.adaptative"), 
#     file = file.path(here::here(), "02_Results/01_Overall_PopGen/00_BasicStats/HeHo_panel.RData"))

load("02_Results/01_Overall_PopGen/00_BasicStats/HeHo_panel.RData")

Hstat_Gen_ZONE.comp <- bind_rows(extract.basic.stats(basic_stat_Gen_ZONE_FG.neutral)  %>% mutate(variable = "Neutral"),
                                 extract.basic.stats(basic_stat_Gen_ZONE_FG.outlier)  %>% mutate(variable = "Outlier"),
                                 extract.basic.stats(basic_stat_Gen_ZONE_FG.adaptative)  %>% mutate(variable = "Adaptative")
                                 
                                 
) %>% left_join(pop.data %>% dplyr::select(Gen_ZONE, Gen_ZONE_FG, RegionAssesment) %>% distinct(.keep_all = T), by = c("ID" = "Gen_ZONE_FG")) %>% 
  mutate(Region = ifelse(Gen_ZONE %in% c("SFA-4", "SFA-5", "SFA-5", "SFA-7"), "NFL", 
                  ifelse(Gen_ZONE == "NAFO-3M", "3M",  RegionAssesment)),
         Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
         
         variable = factor(variable, levels = c("Neutral", "Outlier", "Adaptative"))) 

Hstat_Gen_ZONE.comp 


# One-way anova

# Compute the analysis of variance
res.H.aov <- lm(Ho ~ Region + variable + Region:variable, data = Hstat_Gen_ZONE.comp)
# Summary of the analysis
summary(res.H.aov)
anova(res.H.aov)

# Tukey tests
compute.H.stat <- function(df, snps){
  
  res.aov <- aov(He ~ Region, data = df %>% filter(variable == snps))
  tukey <- TukeyHSD(res.aov)
  
  # compact letter display
  cld <- multcompLetters4(res.aov, tukey)
  cld <- as.data.frame.list(cld$Region) 
  
  res <- cld %>% dplyr::select(Letters) %>% 
    mutate(Region = row.names(cld),
           variable = snps)
  return(res)
}

compute.H.stat(Hstat_Gen_ZONE.comp, "Adaptative")


stat.H.letters <- bind_rows(compute.H.stat( Hstat_Gen_ZONE.comp, "Adaptative"),
                          compute.H.stat( Hstat_Gen_ZONE.comp, "Neutral"),
                          compute.H.stat(Hstat_Gen_ZONE.comp, "Outlier")) 

stat.H.letters <- stat.H.letters %>% left_join(Hstat_Gen_ZONE.comp %>% group_by(Region, variable) %>% 
                                             mutate(maxH = max(He) ) %>% 
                                             dplyr::select(maxH) %>% 
                                             distinct()) %>% 
  mutate (Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M")),
          
          variable = factor(variable, levels = c("Neutral", "Outlier", "Adaptative"))) 


stat.H.letters


Hstat_Gen_ZONE.comp %>% left_join(SVG.res %>% dplyr::select(ID, PAI, SVG, variable)) %>% 
  dplyr::select(Ho, SVG)  %>% 
  dplyr::filter(!is.na(Ho), !is.na(SVG))  %>% cor()

cor.test(Hstat_Gen_ZONE.comp %>% left_join(SVG.res %>% dplyr::select(ID, PAI, SVG, variable)) %>% 
           dplyr::select(Ho, SVG)  %>% 
           dplyr::filter(!is.na(Ho), !is.na(SVG)) %>% pull(Ho), Hstat_Gen_ZONE.comp %>% left_join(SVG.res %>% dplyr::select(ID, PAI, SVG, variable)) %>% 
           dplyr::select(Ho, SVG)  %>% 
           dplyr::filter(!is.na(Ho), !is.na(SVG)) %>% pull(SVG))

cor.test(Hstat_Gen_ZONE.comp %>% left_join(SVG.res %>% dplyr::select(ID, PAI, SVG, variable)) %>% 
           dplyr::select(He, SVG)  %>% 
           dplyr::filter(!is.na(He), !is.na(SVG)) %>% pull(He), Hstat_Gen_ZONE.comp %>% left_join(SVG.res %>% dplyr::select(ID, PAI, SVG, variable)) %>% 
           dplyr::select(He, SVG)  %>% 
           dplyr::filter(!is.na(He), !is.na(SVG)) %>% pull(SVG))


# Genetic offset - Constrained ----------------------------------------------------------
#library(RColorBrewer)

## Function to predict genomic offset from a RDA model
source("./01_Scripts/Functions/genomic_offset.R")
## Running the function for 2050 and 2080

RDA_outliers <- rda(freq.MAF.final.Gen_ZONE_FG[, outliers_rdadapt_env] ~ SST.Larval + Tbtm.Summer + Tbtm.Winter + SSS.ann,  Variables)


res_RDA_proj75 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_current, env_fut = ras_future, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
#res_RDA_proj2050 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_6190, env_fut = ras_2050, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

## Table global genetic offset predicted for 2050 and 2080
RDA_proj_offset <- data.frame(rbind(rasterToPoints( projectRaster(from = res_RDA_proj75$Proj_offset_global, crs =  newproj))), 
                              Date = c(rep("2075", nrow(rasterToPoints(projectRaster(from = res_RDA_proj75$Proj_offset_global, crs =  newproj))))))

hist(RDA_proj_offset$Global_offset)

## Projecting genomic offset on a map
#colors <- c(colorRampPalette(brewer.pal(11, "Spectral")[6:5])(2), colorRampPalette(brewer.pal(11, "Spectral")[4:3])(2), colorRampPalette(brewer.pal(11, "Spectral")[2:1])(3))
gg.offset.map <- ggplot(data = RDA_proj_offset) + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
 # geom_raster(aes(x = x, y = y, fill = cut(Global_offset, breaks=7, include.lowest = T)), alpha = 1) + 
#  scale_fill_manual(values = colors, labels = c("1-2","2-3","3-4","4-5","5-6","6-7","7-8"), guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_raster(aes(x = x, y = y, fill = Global_offset)) +
  #geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(-3, 4, length.out=10), include.lowest = T))) + 
  scale_fill_binned(n.breaks=10, type ="viridis", alpha = 0.7, option = "H",
                    name = "Global offset") +
  
#  scale_fill_distiller(palette = "Spectral", limits = c(0,1.5)) +
  
  geom_sf(data = admin, fill=NA, size=0.1) +  
  geom_sf(data = pop.rda, cex = 0.5) +  


  
    #  # SFA
  #  geom_path(data = SFA_df, 
  #            aes(x = long, y = lat, group = group),
  #            color = "gray25", size = .2) +
  
  
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +

  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "black",
          lwd = 0.5, fill = NA) +  
 # geom_sf(data = pop.rda %>% dplyr::filter(Gen_ZONE_FG %in% c("SFA-13-2", "SFA-15-2", "SFA-16-1")),
#          pch = 24, fill = viridis::viridis(n = 10)[10], cex = 10 ) +
  
  
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  
  geom_sf_label(data = Region.shp, aes(label = layer),
                col = c("black","magenta", "darkorange","blue","red"),
                fill = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-300000),
                nudge_y = c(-900000,650000,100000,800000,220000)) +
  
  
  xlab("Longitude") + ylab("Latitude") +
  #facet_grid(~ Date) +
  facet_wrap(~"Global genomic offset") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(),
        strip.background = element_rect(fill = "white"))
  
gg.offset.map



#ggResponse %>%   aplot::insert_right(gg.SVG.adapt +  theme_minimal(), width =.3) 
 
# Extract the offset by point

Offset.point <- data.frame(raster::extract(res_RDA_proj75$Proj_offset_global, Env.data[,c("Long", "Lat")]))
Offset.point$Gen_ZONE_FG <- Env.data$Gen_ZONE_FG

Offset.point <- Offset.point %>% left_join(pop.data %>% dplyr::select(Gen_ZONE_FG, Gen_ZONE, Region = RegionAssesment) %>% distinct(.keep_all = T))  %>% 
  mutate(Region = ifelse(Gen_ZONE %in% c("SFA-4", "SFA-5", "SFA-5", "SFA-7"), "NFL", 
                         ifelse(Gen_ZONE == "NAFO-3M", "3M",  Region)),
         Region = factor(Region, levels = c("Northern", "NFL", "GSL", "Scotian", "3M"))) 


names(Offset.point)[1] <- "Offset"

Offset.point %>% group_by(Region) %>% summarise(Mean = mean(Offset), sd = sd(Offset)) %>% arrange(Mean)


cor.test(Offset.point %>% left_join(adapt.stat(FreqData, outliers_rdadapt_env.cor),
                           by = c("Gen_ZONE_FG" = "ID")) %>% dplyr::select(Offset, SVG)  %>% 
              dplyr::filter(!is.na(Offset), !is.na(SVG)) %>% pull(SVG),
         Offset.point %>% left_join(adapt.stat(FreqData, outliers_rdadapt_env.cor),
                                   by = c("Gen_ZONE_FG" = "ID")) %>% dplyr::select(Offset, SVG)  %>% 
           dplyr::filter(!is.na(Offset), !is.na(SVG)) %>% pull(Offset)  
         
)


ggResponse <- Offset.point %>% left_join(adapt.stat(FreqData, outliers_rdadapt_env.cor),
                           by = c("Gen_ZONE_FG" = "ID")) %>% 
                 left_join( env.change.df %>% dplyr::filter(name == "SST.Larval")) %>% 
  ggplot(aes(x = SVG, y = Offset )) + 
 geom_smooth(method = "lm", col = "gray10", se = F) +
  geom_point(aes(col = Region, fill = Region), alpha = 0.5, size = 2) +
  scale_color_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"),
                     labels = c("ARC", "NFL", "GSL", "SS", "FLC"))+ 
  
  labs(x = "aSGV", y = "Genomic offset") +
  facet_wrap(~"Correlation") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))
ggResponse


gg.evol <- ggpubr::ggarrange(gg.offset.map + theme(legend.position = "bottom",
                                                   legend.title = element_blank(),
                                                   legend.text = element_text(angle = 90, hjust = 1)), 
                             ggarrange(
                               gg.SVG.adapt,
                               ggResponse,
                               nrow = 2, ncol = 1,
                               align = "hv", labels = c("", "C")),
                             labels = LETTERS,
                             widths = c(3,2) )

gg.evol

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "Fig6_evolution_20231017.png"), 
       plot = gg.evol, bg = "white",
       height = 6, width = 7, units = "in") 


# EEMS --------------------------------------------------------------------
source("./01_Scripts/Functions/EEMS.R")

# Create the different datasets

# Genetics

all.diff <- bed2diffs_v2(tab(gl.final, NA.method = c("asis")))
sort(round(eigen(all.diff)$values, digits = 2)) %>% tail()

neutral.diff <-  bed2diffs_v2(tab(gl.final[, putative.neutral.final], NA.method = c("asis")))
sort(round(eigen(neutral.diff)$values, digits = 2)) %>% tail()

outlier.diff <-  bed2diffs_v2(tab(gl.final[, putative.outliers.final], NA.method = c("asis")))
sort(round(eigen(outlier.diff)$values, digits = 2)) %>% tail()

#write.table(all.diff, file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "all.diffs"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(neutral.diff, file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "neutral.diffs"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(outlier.diff, file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "outlier.diffs"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Coordinates

eems.data <- data.frame(ID_GQ = indNames(gl.final)) %>% 
  left_join(pop.data %>% dplyr::select(ID_GQ, Long, Lat))

#write.table(eems.data %>% select(Long, Lat),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "all.coord"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(eems.data %>% select(Long, Lat),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "neutral.coord"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(eems.data %>% select(Long, Lat),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "outlier.coord"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Polygon habitat

poly.df <- read.csv(file = file.path("00_Data/99_SIG/EEMS_poly_points.csv"))
poly.df <- poly.df %>% arrange(order)

#write.table(poly.df %>% select(X, Y),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "all.outer"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(poly.df %>% select(X, Y),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "neutral.outer"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

#write.table(poly.df %>% select(X, Y),  file.path(here::here(), "02_Results/01_Overall_PopGen/06_EEMS/", "outlier.outer"), 
#            col.names = FALSE, row.names = FALSE, quote = FALSE)

length(putative.outliers.final)
length(putative.neutral.final)


# Create all the parameter files

for(i in 1:3){
  
  for(j in c(400, 600, 800, 1000))
    
    # Write the parameters
    cat(paste("datapath =", file.path(here::here(),"./02_Results/01_Overall_PopGen/06_EEMS/all")),
        paste("mcmcpath =", file.path(here::here(),"./02_Results/01_Overall_PopGen/06_EEMS", paste0("all_nDemes", j, "-chain",i))),
        paste("nIndiv =", 1513),   
        paste("nSites =", 14331), # all loci
       # paste("nSites =", 1552),   # all outliers
       # paste("nSites =", 12779), # neutrals
        paste("nDemes =",j),
        "diploid = true",
       # "numMCMCIter = 2000000",
       # "numBurnIter = 1000000",
       # "numThinIter = 9999",
       # NEW PARAMETERS
       "numMCMCIter = 3000000",
       "numBurnIter = 1000000",
       "numThinIter = 9999",
       # Variance for migration parameters
       paste("mEffctProposalS2 =", 0.3), # default = 0.1
       paste("mSeedsProposalS2 =", 0.01), # default = 0.01       
       paste("mrateMuProposalS2 =", 0.005), # default = 0.01           
       # Variance for diversity parameters
       paste( "qEffctProposalS2 =", 0.002),  # default = 0.001 
       paste("qSeedsProposalS2 =", 0.1), 
       # Variance df
       # paste( "dfProposalS2 =", 0.1),
       
        sep = "\n",
        file = file.path(here::here(),"./02_Results/01_Overall_PopGen/06_EEMS", paste0("all_params-nDemes", j, "-chain",i,".ini"))
    )
  
}


param.ini <- list.files(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/"),
                        pattern = ".ini", full.names = T) %>% 
                        str_subset("all_params|neutral_params|outlier_params") %>% 
                        str_subset("nDemes600|nDemes1000")

param.ini

eems.path <- file.path("/home/genobiwan/Documents/Programs/eems/runeems_snps/src/runeems_snps")

mclapply(param.ini , # %>% str_subset("outlier") %>% str_subset("1200"),
         FUN = function(x){
           
           seed <- str_sub(x, nchar(x)-4, nchar(x)-4) 
           
           cmd <- paste("--params", x,
                        "--seed", paste0(seed,23))
           
           A <- system2(eems.path, cmd, stdout = T, stderr = T)
           A
           
           cat(file = x %>% str_replace(".ini", ".log"),
               cmd, "\n\n",
               A, # what to put in my file
               append= T, sep = "\n")
         },
         mc.cores = 18
)          

# EEMS plot ---------------------------------------------------------------

## Part 2: Generate graphics
projection_none <- "+proj=longlat +datum=WGS84"
projection_out <- "+proj=lcc +datum=WGS84 +lon_0=-90 +lat_1=33 +lat_2=45"


for(p in c("all", "neutral", "outlier")){

#for(i in c(1)){
 
i <- 1:3 # mix the results from the 3 chains  
   
  for(j in c(800)){
    
    mcmcpath = file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS", paste0(p,"_nDemes", j, "-chain",i))
    plotpath = file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Plot", paste0(p,"_nDemes", j))
    
#    plotpath = file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Plot", paste0(p,"_nDemes", j, "-chain",i))
    
    print(plotpath)
    suppressMessages( # because take too much times
      eems.plots(mcmcpath , #= "/home/ian/Software/eems/runeems_snps/data/Crab_Neutral_NoPLBa/",
                 plotpath, # = paste0("/home/ian/Software/eems/runeems_snps/data/Crab_Neutral_NoPLBa/", "1000-demes_withMap"),
                 longlat = TRUE,
                 add.grid = T,
                 lwd.map = 1,
                 add.title = T,
                 #add.outline = T,
                 projection.in = projection_none,
                 add.colbar = T,
                 projection.out = projection_none,
                 #projection.out = projection_out,
                 add.map = T,
                 col.map = "black",
                 add.demes = T)
    )
  }
#}

}

# Extract the raster only

# the polygon

crop_extent <- rgdal::readOGR("00_Data/99_SIG/range2.shp")

plot(crop_extent,
     main = "Shapefile imported into R - crop extent",
     axes = TRUE,
     border = "blue")

for(p in c("all", "neutral", "outlier")){
  
  #for(i in c(1)){
   
  i <- 1:3
   
    for(j in c(800)){
      
      mcmcpath.temp <- file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS", paste0(p,"_nDemes", j, "-chain",i))
       
      print(mcmcpath.temp)
            
      # Values
      ras.m <- extract.EEMS.raster( mcmcpath.temp, is.mrates = TRUE)
      ras.q <- extract.EEMS.raster( mcmcpath.temp, is.mrates = FALSE)

      raster::crs(ras.m) <- "+proj=longlat +datum=WGS84 +no_defs"
      raster::crs(ras.q) <- "+proj=longlat +datum=WGS84 +no_defs"
            
      ras.mask.m <- mask(ras.m, crop_extent)
      ras.mask.q <- mask(ras.q, crop_extent)
      
      writeRaster(ras.mask.m, file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0(p,"_nDemes", j, "_m.tif")),options=c('TFW=YES'),overwrite=TRUE)
      writeRaster(ras.mask.q, file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0(p,"_nDemes", j, "_q.tif")),options=c('TFW=YES'), overwrite=TRUE)
    
      # Probabilities
      ras.prob.m <- extract.EEMS.prob.raster( mcmcpath.temp, is.mrates = TRUE)
      ras.prob.q <- extract.EEMS.prob.raster( mcmcpath.temp, is.mrates = FALSE)
      
      raster::crs(ras.prob.m) <- "+proj=longlat +datum=WGS84 +no_defs"
      raster::crs(ras.prob.q) <- "+proj=longlat +datum=WGS84 +no_defs"
      
      ras.prob.mask.m <- mask(ras.prob.m, crop_extent)
      ras.prob.mask.q <- mask(ras.prob.q, crop_extent)
      
      writeRaster(ras.prob.mask.m, file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0(p,"_nDemes", j, "_m_prob.tif")),options=c('TFW=YES'),overwrite=TRUE)
      writeRaster(ras.prob.mask.q, file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0(p,"_nDemes", j, "_q_prob.tif")),options=c('TFW=YES'), overwrite=TRUE)
      
      
            
     }
#  }
  
}


# The polygon

convert.ras <-function(file, newproj  = "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "){
  
  ras <- terra::rast(file)
  terra::crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs"
  ras.proj <- terra::project(ras, "epsg:4269")
  return(ras)
}


all.ras.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("all","_nDemes", 800, "_m.tif")))
all.ras.prob.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("all","_nDemes", 800, "_m_prob.tif")))

neutral.ras.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("neutral","_nDemes", 800, "_m.tif")))
neutral.ras.prob.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("neutral","_nDemes", 800, "_m_prob.tif")))

outlier.ras.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("outlier","_nDemes", 800, "_m.tif")))
outlier.ras.prob.m_df <- convert.ras(file.path(here::here(),"02_Results/01_Overall_PopGen/06_EEMS/Raster", paste0("outlier","_nDemes", 800, "_m_prob.tif")))

plot(all.ras.prob.m_df )

all.ras.prob.contour <- terra::as.contour(all.ras.prob.m_df, nlevels=20)
summary(all.ras.prob.contour)
table(all.ras.prob.contour$level)
all.ras.prob.contour$level
all.ras.prob.contour <- all.ras.prob.contour[c(1,19),]
plot(all.ras.prob.contour, "level")

all.ras.prob.contour$level
all.ras.prob.contour$PP <- c( "P{log(m) < 0} = 0.95", "P{log(m) > 0} = 0.95")


eems.data.demes <-  bind_cols(eems.data, read_table("02_Results/01_Overall_PopGen/06_EEMS/all_nDemes800-chain1/ipmap.txt", col_names = F))
eems.data.demes


demes.data <- read.table(file = file.path("02_Results/01_Overall_PopGen/06_EEMS/all_nDemes800-chain1/demes.txt"))
names(demes.data) <-  c("Long", "Lat")
demes.data$X1 <- 1:nrow(demes.data)
demes<- sf::st_as_sf(demes.data, coords = c("Long", "Lat"), crs ="+proj=longlat")

eems.data.demes <- eems.data.demes %>% left_join(demes.data, by = "X1") #%>% head()


pop.eems <- sf::st_as_sf(eems.data.demes %>% group_by(Long.y, Lat.y) %>% summarise(N = n()), coords = c("Long.y", "Lat.y"), crs = "+proj=longlat")

admin <- ne_countries(scale = "medium", returnclass = "sf", continent = c("north america"))
admin <- terra::vect(rnaturalearth::ne_countries(scale = "medium", continent = c("north america"), returnclass	  = "sf"))

# SFA

SFA.shp <- terra::vect("00_Data/99_SIG/SFA.shp")
Region.shp <-  terra::vect("00_Data/99_SIG/Regions_MS_corrected.shp")


# Try a graph

library(tidyterra)

devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)

gg.eems.m.all <- ggplot() +
  geom_spatraster(data = all.ras.m_df) +

  # SFA
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +
  
  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
         color = "black",
          lwd = 0.5, fill = NA) +  

  geom_spatraster_contour(data = all.ras.prob.m_df, 
                          aes(color = factor(after_stat(level))), 
                          breaks = c(0.05, 0.95), cex = 1, alpha = 1 , na.rm = T
                          ) +
#  geom_spatraster_contour(data = all.ras.prob.m_df, 
#                          colour = c("dodgerblue"), 
#                          breaks = c(0.95), cex = 1.5, alpha =1 , na.rm = T) +
  # EEMS raster
  #geom_sf(data = all.ras.m_df , aes(x = x, y = y, fill = value)) + 
  
  scale_fill_gradient2(name = "log(m)", low = "red", mid = "lightyellow", 
                       high = "dodgerblue", limits = c(-3,3),  na.value=NA) +
    

  scale_color_manual(name = "Posterior Probabilities", values = c("red", "dodgerblue"),
  label = c( "P{log(m)<0}=0.95", "P{log(m)>0}=0.95" )) +
  
  # Country  
  geom_sf(data = admin, fill="gray55", colour = NA, size=0.5, alpha = 1) +  

    # Demes
  geom_sf(data = demes, alpha = 1/5, pch = 4, size = 1) +  

  # Pop
  geom_sf(data = pop.eems) +
  
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "gray25",
          lwd = 0.5, fill = NA) +  
  geom_sf_label(data = Region.shp, aes(label = layer),
                fill = c("black","magenta", "darkorange","blue","red"),
                col = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-150000),
                nudge_y = c(-900000,650000,100000,800000,-50000)) +
  # Map limits
  scale_size_continuous( range = c(0.5, 3)) +
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  #coord_sf( crs = sf::st_crs("EPSG:6622"))+
  # Others
   xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~"Complete SNP panel") +
   theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
        )+
  guides(
    colour = guide_legend(order = 2),
    fill = guide_colorbar(order = 1),
    size = guide_legend(title = "Sampling site (N)", order = 3))

gg.eems.m.all

gg.eems.m.neutral <-ggplot() +
  geom_spatraster(data = neutral.ras.m_df) +
  
  # SFA
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +
  
  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "black",
          lwd = 0.5, fill = NA) +  
  
  geom_spatraster_contour(data = neutral.ras.prob.m_df, 
                          aes(color = factor(after_stat(level))), 
                          breaks = c(0.05, 0.95), cex = 1, alpha = 1 , na.rm = T
  ) +
  #  geom_spatraster_contour(data = all.ras.prob.m_df, 
  #                          colour = c("dodgerblue"), 
  #                          breaks = c(0.95), cex = 1.5, alpha =1 , na.rm = T) +
  # EEMS raster
  #geom_sf(data = all.ras.m_df , aes(x = x, y = y, fill = value)) + 
  
  scale_fill_gradient2(name = "log(m)", low = "red", mid = "lightyellow", 
                       high = "dodgerblue", limits = c(-3,3),  na.value=NA) +
  
  
  scale_color_manual(name = "Posterior Probabilities", values = c("red", "dodgerblue"),
                     label = c( "P{log(m)<0}=0.95", "P{log(m)>0}=0.95" )) +
  
  # Country  
  geom_sf(data = admin, fill="gray55", colour = NA, size=0.5, alpha = 1) +  
  
  # Demes
  geom_sf(data = demes, alpha = 1/5, pch = 4, size = 1) +  
  
  # Pop
  geom_sf(data = pop.eems) +
  
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "gray25",
          lwd = 0.5, fill = NA) +  
  geom_sf_label(data = Region.shp, aes(label = layer),
                fill = c("black","magenta", "darkorange","blue","red"),
                col = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-150000),
                nudge_y = c(-900000,650000,100000,800000,-50000)) +
  scale_size_continuous( range = c(0.5, 3)) +
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  #coord_sf( crs = sf::st_crs("EPSG:6622"))+
  # Others
  xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~"Neutral SNP panel") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )+
  guides(
    colour = guide_legend(order = 2),
    fill = guide_colorbar(order = 1),
    size = guide_legend(title = "Sampling site (N)", order = 3))
gg.eems.m.neutral 

gg.eems.m.outlier <-ggplot() +
  geom_spatraster(data = outlier.ras.m_df) +
  
  # SFA
  geom_sf(data = SFA.shp,  color = "gray25", lwd = 0.1, fill = NA) +
  
  # REGION
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "black",
          lwd = 0.5, fill = NA) +  
  
  geom_spatraster_contour(data = outlier.ras.prob.m_df, 
                          aes(color = factor(after_stat(level))), 
                          breaks = c(0.05, 0.95), cex = 1, alpha = 1 , na.rm = T
  ) +
  #  geom_spatraster_contour(data = all.ras.prob.m_df, 
  #                          colour = c("dodgerblue"), 
  #                          breaks = c(0.95), cex = 1.5, alpha =1 , na.rm = T) +
  # EEMS raster
  #geom_sf(data = all.ras.m_df , aes(x = x, y = y, fill = value)) + 
  
  scale_fill_gradient2(name = "log(m)", low = "red", mid = "lightyellow", 
                       high = "dodgerblue", limits = c(-3,3),  na.value=NA) +
  
  
  scale_color_manual(name = "Posterior Probabilities", values = c("red", "dodgerblue"),
                     label = c( "P{log(m)<0}=0.95", "P{log(m)>0}=0.95" )) +
  
  # Country  
  geom_sf(data = admin, fill="gray55", colour = NA, size=0.5, alpha = 1) +  
  
  # Demes
  geom_sf(data = demes, alpha = 1/5, pch = 4, size = 1) +  
  
  # Pop
  geom_sf(data = pop.eems) +
  
  geom_sf(data = Region.shp, # color = c("black","magenta", "darkorange","blue","red"),
          color = "gray25",
          lwd = 0.5, fill = NA) +  
  geom_sf_label(data = Region.shp, aes(label = layer),
                fill = c("black","magenta", "darkorange","blue","red"),
                col = "white",
                alpha = 0.75,
                nudge_x = c(400000,-150000,-300000,-200000,-150000),
                nudge_y = c(-900000,650000,100000,800000,-50000)) +
  # Map limits
  scale_size_continuous( range = c(0.5, 3)) +
  coord_sf(xlim = c(-20000, 1850000), ylim = c(-50000, 2650000), crs = sf::st_crs("EPSG:6622")) +
  #coord_sf( crs = sf::st_crs("EPSG:6622"))+
  # Others
  xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~"Outlier SNP panel") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )+
  guides(
    colour = guide_legend(order = 2),
    fill = guide_colorbar(order = 1),
    size = guide_legend(title = "Sampling site (N)", order = 3))

gg.eems.m.outlier 

gg.eems.m1 <- ggarrange(gg.eems.m.neutral, gg.eems.m.outlier ,
                         nrow = 1, ncol = 2, labels = LETTERS, common.legend = T,
                       legend = "right")
gg.eems.m1

ggsave(plot = gg.eems.m1 + theme(plot.margin = margin(0,5,0,0, "mm")),
       filename = "02_Results/01_Overall_PopGen/EEMS_migration_20231006.png",
       height = 6, width = 10, units = "in", bg = "white")


ggsave(plot = gg.eems.m.all,
       filename = "02_Results/01_Overall_PopGen/EEMS_migration_ALL_20231006.png",
       height = 6, width = 6, units = "in", bg = "white")

# Admixture - ALL ---------------------------------------------------------------

current.wd <- here::here()

# Original SNP

all.bed.file <- file.path(here::here(), "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.bed" )
all.fam.file <- all.bed.file %>% str_replace(".bed", ".fam")
all.fam <- read.table(all.fam.file)

for(k in 1:20){

print(k)  
  
setwd("./02_Results/01_Overall_PopGen/10_Admixture/")  

log.path <- paste0("populations.14331snps.1513ind.k",k, ".log")

if(file.exists(log.path)){
cat("\n This run was already done!\n")
  
}
# Run only if it doesn't exist yet
if(!file.exists(log.path)){

  cmd <- paste("--cv", # to perform cross-validation in the log file 
               all.bed.file,
               k, # the number of K
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("populations.14331snps.1513ind.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")

}
  
setwd(current.wd)
    
}

# Cross-validation results:

CV.all.res <- data.frame(k = 1:20,
                     CV = NA,
                    stringsAsFactors = F)


for(i in 1:nrow(CV.all.res)){
  # Which k
  k <- CV.all.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_Overall_PopGen/10_Admixture/", paste0("populations.14331snps.1513ind.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.all.res[i, "CV"] <- CV
  
}

CV.all.res$CV <- as.numeric(as.character(CV.all.res$CV))


gg.CV.all <- CV.all.res %>% CV.plot + facet_wrap(~"Complete SNP panel")
gg.CV.all


gg.CV.all <- CV.all.res %>% mutate(CV = as.numeric(as.character(CV)),
              color = ifelse(k == 4, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  facet_wrap(~"Complete SNPs panel") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) 

#ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen" ,"10_Admixture", "Admixture.global.CV.png"), 
#       plot = gg.CV,
#       height = 3.5, width = 4, units = "in")   



Q.res <-  read.table(file.path(here::here(), "02_Results/01_Overall_PopGen/10_Admixture/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.4.Q"))

head(Q.res)

Q.res <- cbind(all.fam$V1, Q.res)

names(Q.res) <- c("ID_GQ", "Q1", "Q2", "Q3", "Q4")

gg.admix.k4 <- Q.res %>% pivot_longer(cols = c("Q1", "Q2", "Q3", "Q4"), names_to = "Group", values_to = "Q") %>% 
 mutate(Group = factor(Group, levels = c("Q1","Q3", "Q2", "Q4"))) %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_bar(stat="identity") +
  facet_grid(.~ Gen_ZONE, space = "free", scale = "free", switch = "x") +
 scale_fill_manual(values = c( "magenta", "orange",  "deepskyblue", "blue")) +
  labs(y="Membership probability") +
  theme_minimal( base_size = 11) + 
   theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        #panel.spacing = unit(0, "cm"),
       # panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        
        panel.spacing.x = unit(0.2, "lines"),
        
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )

gg.admix.k4


ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen" ,"10_Admixture", "Admixture.global.k4.byArea.png"), 
       plot = gg.admix.k4,
       height = 3.5, width = 12, units = "in")   




Q.pie <- Q.res %>% pivot_longer(cols = c("Q1", "Q2", "Q3", "Q4"), names_to = "Group", values_to = "Q") %>% 
  mutate(Group = factor(Group, levels = c("Q1","Q3", "Q2", "Q4"))) %>% 
  left_join(pop.data) %>% 
  group_by(Gen_ZONE_FG, Group) %>% 
  summarise(value = mean(Q),
            N = n(),
            Lat = mean(Lat),
            Long = mean(Long))

library(scatterpie)
gg.map <- ggplot() + 
   geom_sf(data = admin, fill=NA, size=0.1) +
geom_scatterpie(aes(x=Long, y=Lat, group = Gen_ZONE_FG, r = N/60), 
                data = Q.pie, cols = "Group",  long_format = T) +
  coord_sf(xlim = c(-70, -43), ylim = c(43, 68), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

  

gg.map + geom_scatterpie(aes(x=Long, y=Lat, group = Gen_ZONE_FG, r = N/60), 
                data = Q.pie, cols = "Group",  long_format = T) 



ggpubr::ggarrange(  ggarrange(gPCA.final.14331snps + theme(legend.position = "none"),
                              gg.pca.all + theme(legend.position = "none"), 
                              het.gg.area + theme(legend.position = "bottom"),
                              g1.fst,  
                             
                              nrow =1, ncol = 4, labels = LETTERS,
                              widths = c(1,1,1,2)), 
                    gg.admix.k4, nrow=2, ncol = 1, labels = c("", "E"))


gg.pop.struct <- ggpubr::ggarrange(  ggarrange(gPCA.final.14331snps + theme(legend.position = "none"
                                                                           ) ,
                              gPCA.final.Neutralsnps+ theme(legend.position = "none"), 
                              gPCA.final.OutlierHighsnps + theme(legend.position = "none"),
                              het.gg.area + theme(legend.position = "right"),
                              #g1.fst,  
                              
                              nrow =1, ncol = 4, labels = LETTERS,
                              widths = c(1,1,1,2)), 
                    gg.admix.k4, nrow=2, ncol = 1, labels = c("", "E"),
                    heights = c(2,3))

gg.pop.struct

ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen" , "Fig2.popstruct.png"), 
       plot = gg.pop.struct,
       height = 5, width = 12, units = "in")   

gg.pop.struct 


# Admixture - Outlier ------ ----------------------------------------------
# Outliers SNPs

write.csv(putative.outliers.final, file.path( "02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.csv"), 
          row.names = F, quote = F)

#write.csv(SCAFFOLD.info.final, file.path(filter.ref.path,"06g_UniqueFinal", "Scaffold.info.n13.DP.r5.unique.final.csv"), 
#          row.names = F, quote = F)

# CREATE VCF WITH UNIQUE

cmd <- paste("--vcf", file.path(here::here(), "00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.vcf"), 
             "--recode",
             "--snps", file.path( "02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.csv"),
             "--out", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final")
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)


# Save as plink tped too for pcadapt

cmd1 <- paste("--vcf", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode.vcf"), 
              #"--recode",
              "--plink-tped",
              "--out", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode.tfam"), 
               "--tped", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode")
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a


bed.file <- file.path(here::here(), "./02_Results/01_Overall_PopGen/10_Admixture/", "putative.outlier.final.recode.bed" )
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)

for(k in 17:20){
  
  print(k)  
  
  setwd("./02_Results/01_Overall_PopGen/10_Admixture/Outliers")  
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("putative.outlier.final.recode.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(current.wd)
  
}


# Cross-validation results:

CV.outliers.res <- data.frame(k = 1:20,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.outliers.res)){
  # Which k
  k <- CV.outliers.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_Overall_PopGen/10_Admixture/Outliers", paste0("putative.outlier.final.recode.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.outliers.res[i, "CV"] <- CV
  
}

CV.outliers.res$CV <- as.numeric(as.character(CV.outliers.res$CV))

plot(CV.outliers.res$CV)

CV.outliers.res %>% arrange(CV)

CV.outliers.res %>% ggplot(aes(x = factor(k), y = CV)) + 
  geom_point() +
  theme_bw()

CV.outliers.res %>% arrange(CV) %>% head()

gg.CV.outlier <- CV.outliers.res %>% CV.plot + facet_wrap(~"Outlier SNP panel")
gg.CV.outlier


k <- 4

Q.outliers.res <-  read.table(file.path(here::here(), "02_Results/01_Overall_PopGen/10_Admixture/Outliers", paste0("putative.outlier.final.recode.",k,".Q")))

head(Q.outliers.res)

Q.outliers.res <- cbind(fam$V1, Q.outliers.res)

names(Q.outliers.res) <- c("ID_GQ", paste0("Q", 1:k))

gg.admix.k4.outlier <- Q.outliers.res %>% pivot_longer(cols =  paste0("Q", 1:k), names_to = "Group", values_to = "Q") %>% 
  mutate(Group = factor(Group, levels = paste0("Q", 1:k))) %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_bar(stat="identity") +
  facet_grid(.~ Gen_ZONE, space = "free", scale = "free", switch = "x") +
  scale_fill_manual(values = c( "magenta", "orange",  "deepskyblue", "blue")) +
  labs(y="Membership probability") +
  theme_minimal( base_size = 11) + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        #panel.spacing = unit(0, "cm"),
        # panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        
        panel.spacing.x = unit(0.2, "lines"),
        
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )


gg.admix.k4.outlier

# Admixture - Neutral -----------------------------------------------------


#write.csv(putative.neutral.final, file.path( "02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.csv"), 
#          row.names = F, quote = F)

#write.csv(SCAFFOLD.info.final, file.path(filter.ref.path,"06g_UniqueFinal", "Scaffold.info.n13.DP.r5.unique.final.csv"), 
#          row.names = F, quote = F)

# CREATE VCF WITH UNIQUE


cmd <- paste("--vcf", file.path(here::here(), "00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.vcf"), 
             "--recode",
             "--snps", file.path( "02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.csv"),
             
             "--out", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final")
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)


# Save as plink tped

cmd1 <- paste("--vcf", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode.vcf"), 
              #"--recode",
              "--plink-tped",
              "--out", file.path("02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode.tfam"), 
               "--tped", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode")
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a


neutral.bed.file <- file.path(here::here(), "./02_Results/01_Overall_PopGen/10_Admixture/", "putative.neutral.final.recode.bed" )
neutral.fam.file <- neutral.bed.file %>% str_replace(".bed", ".fam")
neutral.fam <- read.table(neutral.fam.file)

for(k in 11:20){
  
  print(k)  
  
  setwd("./02_Results/01_Overall_PopGen/10_Admixture/Neutrals")  
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               neutral.bed.file,
               k, # the number of K
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("putative.Neutral.final.recode.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(current.wd)
  
}


# Cross-validation results:

CV.neutral.res <- data.frame(k = 1:20,
                              CV = NA,
                              stringsAsFactors = F)


for(i in 1:nrow(CV.neutral.res)){
  # Which k
  k <- CV.neutral.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_Overall_PopGen/10_Admixture/Neutrals", paste0("putative.Neutral.final.recode.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.neutral.res[i, "CV"] <- CV
  
}

CV.neutral.res$CV <- as.numeric(as.character(CV.neutral.res$CV))

gg.CV.neutral <- CV.neutral.res %>% CV.plot + facet_wrap(~"Neutral SNP panel")
gg.CV.neutral


k <- 3

Q.neutral.res <-  read.table(file.path(here::here(), "02_Results/01_Overall_PopGen/10_Admixture/Neutrals", paste0("putative.neutral.final.recode.",k,".Q")))

head(Q.neutral.res)

Q.neutral.res <- cbind(neutral.fam$V1, Q.neutral.res)

names(Q.neutral.res) <- c("ID_GQ", paste0("Q", 1:k))

gg.admix.k3.neutral <- Q.neutral.res %>% pivot_longer(cols =  paste0("Q", 1:k), names_to = "Group", values_to = "Q") %>% 
  mutate(Group = factor(Group, levels = paste0("Q", 1:k))) %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_bar(stat="identity") +
  facet_grid(.~ Gen_ZONE, space = "free", scale = "free", switch = "x") +
  scale_fill_manual(values = c( "magenta", "orange",  "deepskyblue")) +
  labs(y="Membership probability") +
  theme_minimal( base_size = 11) + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        #panel.spacing = unit(0, "cm"),
        # panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        
        panel.spacing.x = unit(0.2, "lines"),
        
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )



# Admixture - Figures -----------------------------------------------------


gg.admix.3 <-  ggarrange(ggarrange(gg.CV.all, gg.CV.neutral, gg.CV.outlier, 
                                   nrow = 1, ncol = 3),
                         gg.admix.k4 + ggtitle("Complete SNPs panel (K = 4)"), 
                         gg.admix.k3.neutral + ggtitle("Neutral SNP panel (K = 3)"),
                         gg.admix.k4.outlier + ggtitle("Outlier SNP panel (K = 4)"),
                         nrow = 4, ncol = 1, labels = LETTERS, heights = c(3,4,4,4)
)

gg.admix.3



ggsave(filename = file.path(here::here(), "02_Results", "01_Overall_PopGen", "PopStruct_Admixture_Complete_Neutral_outlier_20231017.png"), 
       plot = gg.admix.3,
       height = 10, width = 12, units = "in", bg = "white")   

