# Info --------------------------------------------------------------------

# Filtering pipeline (including STACKS version >2 population)
# NS data - ONLY PB
#
# Audrey Bourret
# 2021-06-30

# Library -----------------------------------------------------------------
library(parallel) # to detect the number of cores
library(tidyverse)

library(vcfR)
library(adegenet)
library(hierfstat)

library(here)

`%nin%` = Negate(`%in%`)

# Data --------------------------------------------------------------------

stacks.ref.path <- "00_Data/05b_Stacks.ref/"
filter.ref.path <- "00_Data/06_Filtering/"

pop.data <- read_csv("00_Data/00_FIleInfos/Projet_Infos.csv")
                                
pop.data %>% head()

# What is the sample size
pop.data %>% pull(Cat_Sample) %>% unique()
pop.data %>% filter(Cat_Sample == "Sample") %>% 
             group_by(Gen_ZONE_FG) %>% 
             summarise(N = n()) #%>% write_csv("clipboard")
pop.data %>% filter(Cat_Sample == "Duplicate") %>% nrow() 

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

# Coverage ----------------------------------------------------------------
cov.data <- read_csv(file.path("02_Results/00_SNP_panel/05b_Stacks.ref/","AllIndividuals_RefGen_NreadsNloci.csv"))

cov.data <- cov.data %>% left_join(pop.data %>% select(sample = ID_GQ, Cat_Sample, Plaque_ID))

graph1.0 <- cov.data %>% mutate(RawRatio = n_used_fw_reads/(Raw/2),
                    LowRawRatio = ifelse(RawRatio < .1, "Low", "Normal" )) %>% 
             ggplot(aes(x = Raw/2, y = n_used_fw_reads, col = Espece, shape = LowRawRatio)) +
             geom_point() +
             theme_bw()
graph1.0


graph1.1 <- cov.data %>%  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = Gen_ZONE)) +
  geom_histogram() +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage") + 
  theme_bw() +
  theme(legend.position = "none")

graph1.1

graph1.2 <- cov.data %>%   ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Gen_ZONE))) +
  geom_point() +
  scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 80000, lty = "dashed", col = "darkgray") +
  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")

graph1.2

graph1.2b <- cov.data %>%   ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Plaque_ID))) +
  geom_point() +
  scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 80000, lty = "dashed", col = "darkgray") +
  facet_wrap(~Plaque_ID, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")

graph1.2b

graph1.3 <- cov.data %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = Gen_ZONE, col = as.factor(Gen_ZONE))) +
  geom_boxplot(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")

graph1.3

# List of species to use for next steps 
# I will use ONLY sample (no duplicate), and also I will perform 
# some stats by Gen_ZONE
ID.coverage <- cov.data %>% filter(mean_cov_ns >= 5,
                                   Cat_Sample == "Sample") %>% 
                            select(sample, Gen_ZONE) %>% 
                            mutate(Gen_ZONE = ifelse(Gen_ZONE == "SFA-0", "EAZ", Gen_ZONE %>% str_remove("-")))    

ID.coverage %>% group_by(Gen_ZONE) %>% summarise(N = n())


write.table(ID.coverage, 
            file = file.path(filter.ref.path, "popmap.coverage_5x.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)



# Filtering --------------------------------------------

# Parameters
r.value       <- 0.75 # Minimum within pop
R.value       <- 0.75 # Minimum overall
maf.value     <- 0.05 # Overall MAF
n.pop         <- ID.coverage %>% pull(Gen_ZONE) %>% unique() %>% length()

#maf.pop.value <- 0.05

# Filtering #1 : -r and -MAF ---------------------------------------------

cmd <- paste("-P", stacks.ref.path, 
             "-M", file.path(filter.ref.path, "popmap.coverage_5x.txt"),
             "--out-path",  file.path(filter.ref.path, "06a_r75_MAF05"),
             "-t", 8,
             "-r", r.value, #              
             "-R", R.value, # 
             "--min-maf", maf.value,
             "--min-populations", n.pop,
             #"--smooth",
             #"--write-single-snp",
             "--vcf",
             "--plink"
)

A <- system2("populations", cmd, stdout=T, stderr=T)
A

# General check : Ho and Fis, He vs Ho outliers ---------------------------------

# Load the post-filtration statistic table

filter.stat <- read.delim(file.path(filter.ref.path, "06a_r75_MAF05", "populations.sumstats.tsv"), 
                          skip=n.pop, sep = "\t", header = T )
names(filter.stat)[1] <- "Locus.ID"

nrow(filter.stat) / 15

summary(filter.stat)
head(filter.stat)

# N locus
filter.stat$Locus.ID %>% unique() %>% length()

# For fun, the distribution of the number of SNPs, by scaffold and locus

filter.stat %>% group_by(Locus.ID) %>% summarise(Nsnps = n() / 15) %>% 
                group_by(Nsnps) %>% summarise(Nloc = n()) %>% View()

filter.stat %>% group_by(Chr) %>% summarise(Nloc = length(Locus.ID %>% unique())  ) %>% 
                group_by(Nloc) %>% summarise(Nscaffold = n()) %>% View()

# Check the distribution of Fis and Ho

filter.stat %>% filter(Locus.ID == "814647")
#scaffold84112,15249,f87087Z15249

filter.stat %>% ggplot(aes(x = Obs.Het, fill=Pop.ID)) +
  geom_histogram() +
  theme_bw()

filter.stat %>% ggplot(aes(x = Fis)) +
  geom_histogram() +
  theme_bw()


filter.stat %>% ggplot(aes(x = Obs.Het, y = Fis, col = Pop.ID)) +
  #geom_point(alpha = 1/100) +
  geom_point()+
  geom_vline(xintercept = 0.5, col = "red", lty = "dashed")+
  geom_hline(yintercept = c(-0.3), col = "red", lty = "dashed") +
  theme_bw()

filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het, col = Fis)) +
  geom_point() +
  geom_vline(xintercept = 0.5, col = "red", lty = "dashed")+
  scale_colour_gradientn(colours=rainbow(4)) +
  geom_abline(slope = 1 ) +
  #geom_hline(yintercept = c(-0.3,0.3), col = "red", lty = "dashed") +
  theme_bw()


 # Filtration #2 : Missing data -----------------------------------------
 
 # # Then more filtering with VCF tools
 
 # Verif : individual and SNPs with > 30% missing data
 
 list.files(file.path(filter.ref.path, "06a_r75_MAF05"))
 
 vcf.path <- file.path(filter.ref.path, "06a_r75_MAF05", "populations.snps.vcf")
 
 cmd1 <- paste("--vcf", file.path(current.wd, vcf.path), 
               #"--out",  "MAX2.NArm",
               #" --max-alleles", 2,
               #"--max-missing", 0.80,
               #"--missing-site",
               "--missing-indv"
               #"--kept-sites"
               #"--recode"
 )
 
 # cmd2 <- paste("--vcf", file.path(current.wd, vcf.path), 
 #               #"--out",  "MAX2.NArm",
 #               #"--max-alleles", 2,
 #               #"--max-missing", 0.70,
 #               "--missing-site"
 #               #"--missing-indv",
 #               #"--kept-sites"
 #               #"--recode"
 # )
 
 
 setwd(file.path(filter.ref.path, "06b_MissingData")) 
 
 A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)
# A2 <- system2("vcftools", cmd2, stdout=T, stderr=T)
 
 A1
# A2
 
 
 cat(file = "populations.filt_GLOBAL_5x_b_Individuals_wMissing_A.log",
     "\n", cmd1, "\n",
     A1, # what to put in my file
     append= F, sep = "\n")
 

 setwd(current.wd)
 
 # Many ind with more 30% missing value
 
 imiss <- read.delim(file.path(filter.ref.path, "06b_MissingData", "out.imiss"), skip=0, sep = "\t", header = T )
 imiss %>% head()
 imiss %>% filter(F_MISS >.3)
 
 imiss %>% mutate(SP = str_sub(INDV,1,2)) %>%  ggplot(aes(x=F_MISS, fill = SP)) + geom_histogram()
 
 
graph2.0 <- imiss %>% left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
   mutate(Gen_ZONE = ifelse(Gen_ZONE == "SFA-0", "EAZ", Gen_ZONE)) %>%     
   ggplot(aes(x = Gen_ZONE, y = F_MISS)) + 
    geom_boxplot(col = "blue") +
    geom_jitter(height = 0, alpha = 0.2, col = "black") +
    geom_hline(yintercept = 0.25, lty = "dashed", col = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

good.ID <- imiss %>% filter(F_MISS <= .25) %>% pull(INDV) %>% as.character()
length(good.ID)
 
pop.data %>% filter(ID_GQ %in% good.ID) %>% #select(ID_GQ, Gen_ZONE) %>% 
   group_by(Gen_ZONE) %>% summarise(N = n())
 

 cmd2 <- paste("--vcf", vcf.path, 
               "--recode",
               paste("--indv", good.ID, collapse = " "),
               "--out",  vcf.path %>% str_replace("06a_r75_MAF05", "06b_MissingData") %>% 
                                      str_replace("populations.snps.vcf", paste0("populations.snps.", length(good.ID), "ind.vcf"))
 )
 
 cmd2
 
 A2 <- system2("vcftools", cmd2, stdout=T, stderr=T)
 A2 %>% tail()
 
 cat(file = file.path(filter.ref.path, "06b_MissingData","VCFtools_RemoveIndMissing.log"),
     "\n", cmd2, "\n",
     A2, # what to put in my file
     append= F, sep = "\n")
 
 list.files(file.path(filter.ref.path, "06b_MissingData"))
 
 vcf.path <- file.path(filter.ref.path, "06b_MissingData", "populations.snps.1658ind.vcf.recode.vcf")

 file.exists(vcf.path)
  
 cmd3 <- paste("--vcf", file.path(current.wd,  vcf.path ), 
               #"--out",  "MAX2.NArm",
               #"--max-alleles", 2,
               #"--max-missing", 0.70,
               "--missing-site"
               #"--missing-indv",
               #"--kept-sites"
               #"--recode"
 )
 
 cmd3
 
 setwd(file.path(filter.ref.path, "06b_MissingData")) 
 
 A3 <- system2("vcftools", cmd3, stdout=T, stderr=T)
 A3
 
  cat(file = "ppopulations.filt_GLOBAL_5x_b_Loci_wMissing_B.log",
      "\n", cmd3, "\n",
     A3, # what to put in my file
    append= F, sep = "\n")
  
 
 setwd(current.wd)
 
 # Change the initial 25% missing to 10% missing
 
 lmiss <- read.delim(file.path(filter.ref.path, "06b_MissingData", "out.lmiss"), header = T, sep = "\t" )
 head(lmiss %>% arrange(desc(F_MISS)))
 
graph2.1 <- lmiss %>% ggplot(aes(x=F_MISS)) + 
                      geom_histogram() + 
                      geom_vline(xintercept = 0.10, lty = "dashed", col = "red") +
                      theme_bw()

 
graph2.1  

good.LOC <- lmiss %>% dplyr::filter(F_MISS < .10) %>% dplyr::select(CHR,POS)
 
# N snps
nrow(good.LOC)
 
# N loci


 cmd4 <- paste("--vcf", file.path(filter.ref.path, "06b_MissingData", "populations.snps.1658ind.vcf.recode.vcf"), 
               "--recode",
               "--max-missing", "0.9",
               "--out", file.path(current.wd, filter.ref.path, "06b_MissingData", paste0("populations.", nrow(good.LOC),"snps", ".1658ind"))
 )
 
 cmd4
 
 
 A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
 A4
 
 cat(file = file.path(filter.ref.path, "06b_MissingData","populations.filt_GLOBAL_5x_b_Remove_Loci.log"),
     "\n", cmd4, "\n",
     A4, # what to put in my file
     append= F, sep = "\n")
 

# Filtration #3 : HW outliers ---------------------------------------

# This part can be long ... (data convertion to hierfstat)

list.files(file.path(filter.ref.path, "06b_MissingData"))
 
# Load the VCF file

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path, "06b_MissingData", "populations.36207snps.1658ind.recode.vcf"))

head(vcf.data)

vcf_field_names(vcf.data , tag = "FORMAT")

# Extract raw info from vcf file
gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))
head(gt.tidy)

# Proablement plus tard, sur filtration sur les individus

gt.meta <- gt.tidy %>% group_by(Indiv) %>%  summarise(Nsnps = length(gt_GT[!is.na(gt_GT)]),
                                                      DP = mean(gt_DP, na.rm=T)) %>%  
  left_join(pop.data, by = c("Indiv" = "ID_GQ"))
  
head(gt.meta)  

#gt.meta %>% filter(str_detect(Indiv, "EtOH")) %>% pull(Indiv) %>% View()
#gt.meta %>% filter(str_detect(Indiv, "1644|1645")) %>% View()

#gt.meta %>% filter(Loci < 3000) %>% View()

gt.meta %>%   ggplot(aes(x = Gen_ZONE, y = Nsnps)) +
  #geom_violin(fill="#C0C0C0", adjust = 1, scale = "count", trim = T) +
 # geom_jitter(height = 0, alpha = 1/5) +
  geom_jitter(height = 0, alpha = 1/5) +
  geom_boxplot(alpha = 0) +
  #  stat_summary(fun.data=mean_sdl, geom = "pointrange", color = "black")+
  #scale_y_continuous(trans = scales::log2_trans(), breaks = c(1,10,100,1000, 10000))+
  labs(y = "N snps", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

gt.meta %>%   ggplot(aes(x = DP, y = Nsnps, col = Plaque_ID)) +
geom_point()+
  #  stat_summary(fun.data=mean_sdl, geom = "pointrange", color = "black")+
  #scale_y_continuous(trans = scales::log2_trans(), breaks = c(1,10,100,1000, 10000))+
  labs(y = "N snps", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Conversion to various R formats
#gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

new.names <- data.frame(ID_GQ = indNames(gi.data))
new.names <- new.names %>% left_join(pop.data) %>% 
                           mutate(Gen_ZONE = ifelse(Gen_ZONE == "SFA-0", "EAZ", Gen_ZONE))
new.names %>% head()



# add names and info on library
#pop(gl.data) <- new.names$Gen_ZONE
pop(gi.data) <- new.names$Gen_ZONE

new.names$Gen_ZONE %>% table()
hf.data <- genind2hierfstat(gi.data) 

hw.res <- pegas::hw.test(gi.data)
hw.res %>% str()

hw.res %>% summary()
hw.res %>% summary()
hw.res %>% as_tibble() %>% filter(Pr.exact <= 0.05) %>% nrow()

14878/36207

# Computing diversity (Ho, He)

div <- summary(gi.data)

str(div)


# Try to do a purr and mapped
start_time <- Sys.time()

HW.pop <- tibble(Pop = pop(gi.data) %>% levels())

HW.pop <- HW.pop %>% mutate(# Return the number of individuals
                            nInd = map_int(Pop, function(d){
                                        nInd(gi.data[pop = d])
                                        }),
                            nLoc = map_int(Pop, function(d){
                              nLoc(gi.data[pop = d])
                            }),
                            # Perform HW on a subset
                             HW = map(Pop, function(d){
                                  gi <- gi.data[pop = d]
                                  #gi <- gi.data[pop = d, loc = 1:10]
                                  res <- pegas::hw.test(gi)
                                  res <- as.data.frame(res)
                                  return(res)
                             })
                          ) 

end_time <- Sys.time()
end_time - start_time

HW.pop$HW


HW.pop <- HW.pop %>% mutate(#How many loci <= 0.5
                            NP05 = map_int(HW, function(d){
                            d %>%  filter(Pr.exact <= 0.05) %>% nrow()
                            }),
                            PP05 = NP05/nLoc,
                            LociP05 = map(HW, function(d){
                            d %>% filter(Pr.exact <= 0.05) %>% row.names()  
                            })
                          )



HW.pop %>% pull(LociP05) %>% unlist() %>% table() %>% table() 
HW.pop %>% pull(LociP05) %>% unlist() %>% table() %>% length()
HW.pop %>% pull(LociP05) %>% unlist() %>% table() %>% hist(main = "N pop (max 15) in which a SNPs is not HW")

# New table with join info

HW.pop %>% pull(LociP05) %>% unlist() %>% table() %>% dimnames()

HW.loc <- tibble(ID = dimnames(HW.pop %>% pull(LociP05) %>% unlist() %>% table())[[1]],
                 Npop = HW.pop %>% pull(LociP05) %>% unlist() %>% table())

HW.loc


div.graph <- data.frame(ID = names(div$Hobs),
                          Hobs = div$Hobs,
                        Hexp = div$Hexp)


div.graph <- div.graph %>% left_join(HW.loc) %>% 
                           mutate(Npop = ifelse(is.na(Npop), 0, Npop))

div.graph %>% head()

#save(list = c("hw.res", "HW.pop", "div.graph"),
#     file = file.path(filter.ref.path, "06c_HW","SNPs_GLOBAL_PB_5x_HW.RES.data"))

load(file.path(filter.ref.path, "06c_HW","SNPs_GLOBAL_PB_5x_HW.RES.data"))

graph3.0 <-  div.graph %>% ggplot(aes(x = Hobs, y = Hexp, col = Npop)) +
                  geom_point() +
                  scale_colour_distiller(palette = "Spectral") +
                  geom_vline(xintercept = 0.5) +
                  labs(title = "Ho vs He considering HW disequillibirum in N pop")+ 
                  theme_bw()
graph3.0

graph3.1 <- div.graph %>% filter(Npop <=2) %>% 
  ggplot(aes(x = Hobs, y = Hexp, col = Npop)) +
  geom_point() +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.5) +
  theme_bw()
graph3.1


# Add info om loci

div.graph <- div.graph %>% mutate(Loc = sapply(str_split(ID, ":"), `[`, 1)) 

Nspns.tab <- div.graph %>% group_by(Loc) %>% 
                           summarise(N = n(),
                           meanNpop = mean(Npop))                    
Nspns.tab %>% ggplot(aes(x = N, y = meanNpop)) + geom_point(alpha = 1/10)
Nspns.tab %>% ggplot(aes(x = N, fill = factor(round(meanNpop,0)))) + geom_histogram()

# Keep 1 snp / rad loci
list.n2.unique <- div.graph %>% filter(Npop <=2)  %>% distinct(Loc, .keep_all = T) 
list.n13.unique <- div.graph %>% filter(Npop <=13)  %>% distinct(Loc, .keep_all = T) 

# Keep ALL
list.n2.all  <- div.graph %>% filter(Npop <=2) 
list.n13.all <- div.graph %>% filter(Npop <=13) 

nrow(list.n2.unique)
nrow(list.n13.unique)

nrow(list.n2.all)
nrow(list.n13.all)


write.csv(list.n2.unique %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.n2.unique.csv"), 
          row.names = F, quote = F)

write.csv(list.n13.unique %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.n13.unique.csv"), 
          row.names = F, quote = F)

write.csv(list.n2.all %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.n2.all.csv"), 
          row.names = F, quote = F)

write.csv(list.n13.all %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.n13.all.csv"), 
          row.names = F, quote = F)


#setwd(current.wd)

# CREATE VCF WITH UNIQUE

cmd1 <- paste("--vcf", file.path(filter.ref.path, "06b_MissingData", "populations.36207snps.1658ind.recode.vcf"), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06c_HW", "Loc.n2.unique.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06c_HW", paste0("populations.",list.n2.unique %>% nrow(),"snps.1658ind.n2HW.single"))
)

cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)
A1

cat(file = file.path(current.wd,filter.ref.path,"06c_HW","VCFtools_SnpWhiteList.n2HW.log"),
    "\n", cmd1, "\n",
    A1, # what to put in my file
    append= F, sep = "\n")


cmd2 <- paste("--vcf", file.path(filter.ref.path, "06b_MissingData", "populations.36207snps.1658ind.recode.vcf"), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06c_HW", "Loc.n13.unique.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06c_HW", paste0("populations.",list.n13.unique %>% nrow(),"snps.1658ind.n13HW.single"))
)

cmd2

A2 <- system2("vcftools", cmd2, stdout=T, stderr=T)
A2

cat(file = file.path(current.wd,filter.ref.path,"06c_HW","VCFtools_SnpWhiteList.n13HW.log"),
    "\n", cmd2, "\n",
    A2, # what to put in my file
    append= F, sep = "\n")

## NOW BY KEEPING ALL

cmd3 <- paste("--vcf", file.path(filter.ref.path, "06b_MissingData", "populations.36207snps.1658ind.recode.vcf"), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06c_HW", "Loc.n2.all.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06c_HW", paste0("populations.",list.n2.all %>% nrow(),"snps.1658ind.n2HW.all"))
)

cmd3

A3 <- system2("vcftools", cmd3, stdout=T, stderr=T)
A3

cat(file = file.path(current.wd,filter.ref.path,"06c_HW","VCFtools_SnpWhiteList.n2HW.log"),
    "\n", cmd3, "\n",
    A3, # what to put in my file
    append= T, sep = "\n")


cmd4 <- paste("--vcf", file.path(filter.ref.path, "06b_MissingData", "populations.36207snps.1658ind.recode.vcf"), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06c_HW", "Loc.n13.all.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06c_HW", paste0("populations.",list.n13.all %>% nrow(),"snps.1658ind.n13HW.all"))
)

cmd4

A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
A4

cat(file = file.path(current.wd,filter.ref.path,"06c_HW","VCFtools_SnpWhiteList.n13HW.log"),
    "\n", cmd4, "\n",
    A4, # what to put in my file
    append= T, sep = "\n")


# General Check #4 : Relatedness ----------------------------------------------

file.path(current.wd,filter.ref.path,"06c_HW") %>% list.files()

# In VCF tool, check for duplicated ind.

cmd <- paste("--vcf", file.path(current.wd, filter.ref.path,"06c_HW", "populations.17093snps.1658ind.n13HW.single.recode.vcf"), 
             #sub_indv,
             # --depth,
             "--relatedness2"
             
)


setwd(file.path(filter.ref.path,"06d_Relatedness")) 

A <- system2("vcftools", cmd, stdout=T, stderr=T)
A

cat(file = file.path("VCFtools.Relatedness.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

related2.all <- read.delim(file.path(filter.ref.path,"06d_Relatedness","out.relatedness2"), header = T, sep = "\t" )

related2.all %>% head()

length(related2.all$INDV1 %>% unique)

# Remove t=identical individuals
related2.ok <- related2.all %>% mutate(iden = ifelse(INDV1 == INDV2, 1, 0)) %>% # Remove identique
  filter(iden == 0) %>% 
  left_join(pop.data, by = c("INDV1" = "ID_GQ"))

related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% arrange(desc(RELATEDNESS_PHI)) %>% View()

related2.ok %>% arrange(desc(RELATEDNESS_PHI)) %>% head()

c(related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% pull(INDV1),
  related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% pull(INDV2)
) %>% table()

graph1 <- related2.ok %>% 
  #filter(RELATEDNESS_PHI > -20) %>% 
  ggplot(aes(x =  RELATEDNESS_PHI)) +
  labs(x = "Relatedness coefficient", y = "N observations", title = "Relatedness coefficient distribution")+
  geom_histogram() +
  geom_vline(xintercept = c(1/2, 1/4, 1/8, 1/16), lty = "dashed")+
  annotate("text", 
           x = c(1/2-0.005, 1/4-0.005, 1/8-0.005, 1/16-0.005), y = 1000, 
           label = c("Individual-self", "Siblings / Parent-offspring", "Half-siblings / Grandparent-grandchild", "First cousins"), 
           angle = 90, hjust = 0, vjust = 0) +
  theme_bw() #+
#  theme(axis.text.x = element_text(angle = 60, hjust = 1))

graph1

ggsave(filename = file.path(filter.ref.path,"06d_Relatedness", "Relatedness.png"),
       plot = graph1,
       width = 7, height = 5, units = "in")


# Filtration #5 : Too much depth ----------------------------------------------------------

vcf.path <- ("00_Data/06_Filtering/06c_HW/populations.34983snps.1658ind.n13HW.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

# Remove loci with too much depth


gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))
head(gt.tidy)

vcf.fix <- as.data.frame(vcf.data@fix) %>% mutate(Key = 1:nrow(.))

# group by snps
gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>% 
                                        left_join(vcf.fix %>% select(Key, ID))


graphDP1 <- gt.key %>% 
  ggplot(aes(x = medianDP, y = meanDP, col = maxDP)) +
  geom_jitter(alpha = 0.5)+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_color_viridis_c() +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graphDP1

graphDP2 <- gt.key %>% 
  ggplot(aes(x = medianDP)) +
 geom_histogram()+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graphDP2

LOC.DP <- gt.key %>% filter(meanDP > 25) %>% pull(ID)

LOC.N13 <- read.csv(file.path(filter.ref.path,"06c_HW", "Loc.n13.all.csv"))

LOC.N13.DP <- LOC.N13 %>% filter(ID %nin% LOC.DP)

nrow(LOC.N13.DP) + length(LOC.DP) == nrow(LOC.N13)

write.csv(LOC.N13.DP, file.path(filter.ref.path,"06e_DP", "Loc.n13.all.DP.csv"), 
          row.names = F, quote = F)


# Do the filtration

cmd5 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06e_DP", "Loc.n13.all.DP.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06e_DP", paste0("populations.",LOC.N13.DP %>% nrow(),"snps.1658ind.n13HW.all"))
)

cmd5

A5 <- system2("vcftools", cmd5, stdout=T, stderr=T)
A5

cat(file = file.path(current.wd,filter.ref.path,"06e_DP","VCFtools_SnpWhiteList.n13HW.DP.log"),
    "\n", cmd5, "\n",
    A5, # what to put in my file
    append= F, sep = "\n")

# #5 - Compute LD ---------------------------------------------------------

vcf.path <- ("00_Data/06_Filtering/06e_DP/populations.34485snps.1658ind.n13HW.DP.all.recode.vcf")
file.exists(vcf.path)

vcf.data <- vcfR::read.vcfR(vcf.path)
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# Plink all LD r 

# Transform.vcf in tped
cmd1a <- paste("--vcf", file.path(vcf.path), 
               #"--recode",
               "--plink-tped",
               "--out",  file.path(vcf.path) %>% str_remove(".recode.vcf") )


cmd1a

A1 <- system2("vcftools", cmd1a, stdout=T, stderr=T)


# Compute LD

cmd2a <- paste(#"--file", "./test.plink.bed",
  "--tfam", "00_Data/06_Filtering/06e_DP/populations.34485snps.1658ind.n13HW.DP.all.tfam", 
  "--tped", "00_Data/06_Filtering/06e_DP/populations.34485snps.1658ind.n13HW.DP.all.tped", 
  # "--allow-extra-chr",
  #"--make-bed",
  "--r2 inter-chr",
  "--ld-window-r2", 0.2,
  # "--indep-pairwise 50000kb 10 0.8", 
  "--out", "00_Data/06_Filtering/06f_LD/LD.34485snps.1658ind.n13HW.DP.all.r2"
)

#A2 <- system2("plink", cmd2a, stdout=T, stderr=T)

# Load results

LD.all <- read.table("00_Data/06_Filtering/06f_LD/LD.34485snps.1658ind.n13HW.DP.all.r2.ld", header = T)
nrow(LD.all)
head(LD.all)
hist(LD.all$R2)

# Transformer le tout en beau jeu de données :
SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  dplyr::select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1),
         scaffold.length = sapply(str_split(CHROM, ","), `[`,2))

SCAFFOLD.info

LD.all.info <- LD.all %>% select(SNP_A, SNP_B, R2) %>% 
  left_join(SCAFFOLD.info %>% select(SNP_A = ID, SCAFFOLD_A = scaffold, RADloc_A = RADloc)) %>% 
  left_join(SCAFFOLD.info %>% select(SNP_B = ID, SCAFFOLD_B = scaffold, RADloc_B = RADloc)) %>% 
  mutate(Cat = ifelse(SCAFFOLD_A == SCAFFOLD_B, ifelse(RADloc_A == RADloc_B, "Intra-RADloc", "Intra-Scaffold"), "Inter-Scaffold"))


LD.all.info %>% ggplot(aes(x = R2, fill = Cat)) +
  geom_density(alpha = 0.5) +
  theme_bw()

LD.all.info %>%  ggplot(aes(x = R2, fill = Cat)) +
  geom_histogram() +
  facet_grid(Cat ~ ., scale = "free_y") +
  theme_bw()

LD.all.info %>% filter(Cat == "Intra-Scaffold") %>% group_by(SCAFFOLD_A) %>% 
  summarise(N = n()) %>% arrange(desc(N))

LD.all.info %>% filter(Cat == "Intra-RADloc") %>% group_by(SCAFFOLD_A) %>% 
  summarise(N = n()) %>% arrange(desc(N))


## RADIATOR

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

library(magrittr)

radiator.LD.modif <- function(plink.LD, gi, r2){
  
  # Create metadata
  na.int <- na.gi.count(gi)
  names(na.int)
  
  markers.missing <- tibble::tibble(
    MARKERS = names(na.int),
    MISSING_PROP = na.int
  ) 
  
  ld.tibble <- plink.LD %>% filter(R2 > r2)
  
  markers.ld.list <- tibble::tibble(
    MARKERS = c(ld.tibble$SNP_A, ld.tibble$SNP_B)) %>%
    dplyr::count(MARKERS) %>%
    # dplyr::group_by(MARKERS) %>%
    # dplyr::tally(.) %>%
    # dplyr::ungroup(.) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::distinct(MARKERS) %>%
    purrr::flatten_chr(.)
  
  # we want to have them ordered from highest to lowest hence the approach above...
  geno.stats <- dplyr::filter(markers.missing, MARKERS %in% markers.ld.list)
  whitelist.markers <- blacklist.markers <- tibble::tibble(MARKERS = character(0))
  
  for (i in markers.ld.list) {
    # i <- markers.ld.list[1]
    dups <- dplyr::filter(ld.tibble, SNP_A %in% i | SNP_B %in% i)
    dups <- sort(unique(c(dups$SNP_A, dups$SNP_B)))
    
    # find all duplicates associated with the network
    new.dups <- 0L
    while(length(new.dups) > 0) {
      new.dups <- dplyr::filter(ld.tibble, SNP_A %in% dups | SNP_B %in% dups)
      new.dups <- sort(unique(c(new.dups$SNP_A, new.dups$SNP_A)))
      new.dups <- purrr::keep(.x = new.dups, .p = !new.dups %in% dups)
      if (length(new.dups) > 0) {
        dups <- c(dups, new.dups)
      }
    }
    dups <- tibble::tibble(MARKERS = dups)
    
    if (nrow(blacklist.markers) > 0) {
      dups <- dplyr::filter(dups, !MARKERS %in% blacklist.markers$MARKERS)
    }
    
    if (nrow(dups) > 0) {
      wm <- dups %>%
        dplyr::left_join(geno.stats, by = "MARKERS") %>%
        dplyr::mutate(MISSING_PROP = ifelse(is.na(MISSING_PROP), 1, MISSING_PROP)) %>% 
        dplyr::filter(MISSING_PROP == min(MISSING_PROP, na.rm = T)) %>%
        dplyr::sample_n(tbl = ., size = 1) %>% # make sure only 1 is selected
        dplyr::select(MARKERS)
      
      if (nrow(wm) > 0) whitelist.markers %<>% dplyr::bind_rows(wm)
      
      bm <- dplyr::filter(dups, !MARKERS %in% wm$MARKERS) %>%
        dplyr::select(MARKERS)
      
      if (nrow(bm) > 0) blacklist.markers %<>% dplyr::bind_rows(bm)
    }
  }
  
  dups <- bm <- wm <- i <- new.dups <- NULL
  
  blacklist.markers <- dplyr::distinct(blacklist.markers, MARKERS)
  
  return(blacklist.markers %>% pull(MARKERS))
  
}

#LD.r.8 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.8)
#LD.r.7 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.7)
#LD.r.6 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.6)
LD.r.5 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.5)
#LD.r.4 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.4)
#LD.r.3 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.3)
#LD.r.2 <- radiator.LD.modif(plink.LD = LD.all, gi = gi.data, r2 = 0.2)



LOC.N13.DP<- read.csv(file.path(filter.ref.path,"06e_DP", "Loc.n13.all.DP.csv"))

LOC.N13.DP.r5 <- LOC.N13.DP %>% filter(ID %nin% LD.r.5)

nrow(LOC.N13.DP.r5) + length(LD.r.5) == nrow(LOC.N13.DP)

write.csv(LOC.N13.DP.r5, file.path(filter.ref.path,"06f_LD", "Loc.n13.all.DP.r5.csv"), 
          row.names = F, quote = F)


# Do the filtration

cmd6 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"06f_LD", "Loc.n13.all.DP.r5.csv"),
              "--out", file.path(current.wd,filter.ref.path,"06f_LD", paste0("populations.",LOC.N13.DP.r5 %>% nrow(),"snps.1658ind.n13HW.DP.r5.all"))
)

cmd6

A6 <- system2("vcftools", cmd6, stdout=T, stderr=T)
A6

cat(file = file.path(current.wd,filter.ref.path,"06f_LD","VCFtools_SnpWhiteList.n13HW.DP.r5.log"),
    "\n", cmd6, "\n",
    A6, # what to put in my file
    append= F, sep = "\n")



# Final dataset -----------------------------------------------------------

vcf.path <- ("00_Data/06_Filtering/06f_LD/populations.20136snps.1658ind.n13HW.DP.r5.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

#gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

SCAFFOLD.info

data.frame(ID_GQ = indNames(gi.data)) %>% left_join(pop.data) %>% 
  group_by(Transcripto) %>% summarise(N = n())

ID.final <-  data.frame(ID_GQ = indNames(gi.data)) %>% 
  left_join(pop.data) %>% 
  filter(Gen_ZONE_FG %nin% c("SFA-6-3", "SFA-8-4", "SFA-12-2", "SFA-12-3", "SFA-15-1")) %>% 
  pull(ID_GQ)

length(ID.final)

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

# Function to create a list of loci, from a genind object

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
  # Create vectors for each loci
  MAF.res <- adegenet::minorAllele(gi)
  NA.res  <- na.gi.count(gi)
  
  # Filter by threshold
  MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
  cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
  
  NA.loc <- names(NA.res[NA.res <= NA.trs])
  cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
  
  # LOCI with both conditions
  LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
  LOCI.res %>% length()
  
  cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
  
  return(LOCI.res)
}


LOC.MAF05.NA10 <- filter.MAF.NA(gi.data[indNames(gi.data) %in% ID.final,] , MAF.trs = 0.05, NA.trs = 0.10)

# Then
SCAFFOLD.info.final <-  SCAFFOLD.info %>%  filter(ID %in% LOC.MAF05.NA10) %>% 
  arrange(scaffold, POS) %>% 
  distinct(RADloc, .keep_all = T)
  


write.csv(SCAFFOLD.info.final%>% select(ID), file.path(filter.ref.path,"06g_UniqueFinal", "Loc.n13.DP.r5.unique.final.csv"), 
          row.names = F, quote = F)

write.csv(SCAFFOLD.info.final, file.path(filter.ref.path,"06g_UniqueFinal", "Scaffold.info.n13.DP.r5.unique.final.csv"), 
          row.names = F, quote = F)


# CREATE VCF WITH UNIQUE

vcf.path

cmd <- paste("--vcf", file.path(filter.ref.path, "06f_LD", "populations.20136snps.1658ind.n13HW.DP.r5.all.recode.vcf"), 
              "--recode",

              paste("--indv",ID.final, collapse = " "),
              "--snps", file.path(current.wd,filter.ref.path, "06g_UniqueFinal", "Loc.n13.DP.r5.unique.final.csv"),
              
              "--out", file.path(current.wd,filter.ref.path, "06g_UniqueFinal", paste0("populations.", SCAFFOLD.info.final %>% nrow(),"snps.",length(ID.final),"ind.n13HW.DP.r5.single.final"))
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)

# A2


cat(file =  file.path(filter.ref.path, "06g_UniqueFinal", "VCF.filter.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Reload, and save as Rdata

vcf.path <- ("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

gl.final  <- vcfR::vcfR2genlight(vcf.data) 
gi.final  <- vcfR::vcfR2genind(vcf.data) 


pop(gl.final) <- data.frame(ID_GQ = indNames(gl.final)) %>% 
                   left_join(pop.data) %>% 
                   pull(Gen_ZONE_FG)

pop(gi.final) <- data.frame(ID_GQ = indNames(gi.final)) %>% 
  left_join(pop.data) %>% 
  pull(Gen_ZONE_FG)

save(list = c("gl.final", "gi.final"),
     file = "00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.adegenet.Rdata")

# Save as plink tped too for pcadapt

cmd1 <- paste("--vcf", vcf.path, 
               #"--recode",
               "--plink-tped",
               "--out",  vcf.path %>% str_remove(".vcf"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.tfam", 
               "--tped", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.tped", 
               "--make-bed", 
               "--out", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode"
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a

### idem with SFA

gl.sfa <- gl.final[pop(gl.final) %nin%  c("NAFO-3M-1", "NAFO-3M-2", "NAFO-3M-3", "NAFO-3M")]
ID.sfa <- indNames(gl.sfa)

cmd <- paste("--vcf", vcf.path, 
             "--recode",
             
             paste("--indv",ID.sfa, collapse = " "),
             "--out", vcf.path %>% str_remove(".recode.vcf") %>% 
                                   str_replace("1513", length(ID.sfa) %>% as.character() ) %>% 
                                   str_replace(".final", ".sfa"))

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)


cat(file =  file.path(filter.ref.path, "06g_UniqueFinal", "VCF.filter.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= T, sep = "\n")


"00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.recode.vcf"


cmd1 <- paste("--vcf", vcf.sfa.path, 
              #"--recode",
              "--plink-tped",
              "--out",  vcf.sfa.path %>% str_remove(".vcf"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1399ind.n13HW.DP.r5.single.sfa.recode.tfam", 
               "--tped", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1399ind.n13HW.DP.r5.single.sfa.recode.tped", 
               "--make-bed", 
               "--out", "./00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1399ind.n13HW.DP.r5.single.sfa.recode"
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a
