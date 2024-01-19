# Info --------------------------------------------------------------------

# RAD-seq pipeline with STACKS version >2, up-to Gstacks
# NS data
# Trimmomatics + no rad check on one side (cut 5 pb) 
# With newly created ref genome
#
# Audrey Bourret
# 2021-01-12
#


# Library -----------------------------------------------------------------

library(parallel)
library(tidyverse)

`%nin%` = Negate(`%in%`)

# Add python environment (containing multiqc) path available in R
Sys.setenv(PATH = paste(c("/home/genobiwan/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("multiqc", "--help")

# Data --------------------------------------------------------------------

pop.info <- read.table(file.path("00_Data/00_FIleInfos/", "popmap.txt"))
names(pop.info) <- c("Sample", "POP")
pop.info


pop.data <- read_csv("00_Data/00_FIleInfos/Projet_Infos.csv") 

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

# FastQC ------------------------------------------------------------------

# Took around 30 min by plate (NS and HI)
# Now with an overwrite function

fastqc <- function(folder.in, folder.out, overwrite = F, nthread) {
#@folder.in : where to find the data
#@folder.out : where to put the data
#@overwrite : if T, will remove everything within the folder  
    
  if(get_os() %in% c("os","linux")){ # to run only on os and not windows
    
    cat("Performing a FastQC analysis\n")
    
    # Remove old files or not (addition June 2020)
    if(isTRUE(overwrite)){
      file.remove(list.files(folder.out, full.name = T, pattern ="fastqc"))
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% 
        str_subset("R1.fastq|R2.fastq") %>% unique()
    }
    
    if(isTRUE(overwrite == F)){
      previous.analysis <- list.files(folder.out, pattern = ".html") %>% str_replace("_fastqc.html", ".fastq.gz")
      
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% unique() %>% 
        str_subset(paste(previous.analysis, collapse = "|"), negate = T) %>% 
        str_subset("R1.fastq|R2.fastq")
      
    }
    
    cat("Results could be find here:", folder.out ,"\n")
    
    mclapply(files.to.use,
             FUN = function(x){
               
               cmd <- paste("--outdir", folder.out, x, 
                            "-t", 1)
               system2("fastqc", cmd) 
              
             } ,
             mc.cores = nthread
    )
    

  } else {cat("Cannot perform FastQC on windows yet -- sorry!!")}
} # End of my function

# Test a multiqc function

multiqc <- function(folder.in){
  # Multi QC aggregation - run pretty fast ...
  for(s in c("R1", "R2")){
    print(s)  
    cmd <- paste(paste(list.files(folder.in, full.names = T) %>% 
                         #str_subset(l) %>% 
                         str_subset(paste0("_",s)) %>% 
                         str_subset(".zip"), collapse = " "),
                 "--outdir", file.path(folder.in, "MultiQC_report"),
                 "--filename", paste0("multiqc_report_", s, ".html"),
                 "-f" # to rewrite on previous data
    )
    
    system2("multiqc", cmd)
    
  } 
  
}


# Run FastQC - change for overwrite T for the first time to or it will not work
fastqc(folder.in = "./00_Data/01_Raw", folder.out = "./02_Results/01_FastQC/01_Raw", overwrite = T, nthread = 20)


# Multi QC aggregation - run pretty fast ...
multiqc(folder.in = get.value("./02_Results/00_SNP_panel/01_FastQC/01_Raw"))

# CHECK on the R2 side that the adapter "CGG" is present. 
# IF not, the PART2 of trimmomatic allows to remove 5 pb on this side

# Trimmomatics ------------------------------------------------------------

# In paired-end
# To remove the Illumina adapter and to cut 3pb in R2

# Check that you can reach the program
trimmomatic.path <- "~/Documents/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
system2("java", paste("-jar", trimmomatic.path, "PE", "-version"), stdout=T, stderr=T) 

# Check which files you want to use

files.to.use <- list.files( "./00_Data/01_Raw", full.names = T) %>% 
                           str_subset("R1") %>% 
                           str_subset(".md5", negate = T) 
files.to.use

#start_time3 <- Sys.time()

### BENCHMARK ###
# On 4 files of 1000 reads
# 80 cores in Trimmomatics: 1.249837 secs
# 8 cores in Trimmomatics: 1.179123 secs 
# 4 cores in Trimmomatics / parallele in R 4 cores: 0.3734035 secs
#  1 core in Trimmomatics / parallele in R 4 cores: 0.3527803 secs

# With 60 cores, memory is full, so 20 seems a better compromise
# checked in Terminal with "free -m" or "top"

# STEP 1 - REMOVE 3 PB on R2 (must do it first because of misalignment)

  mclapply(files.to.use,
         FUN = function(x){
           
           # Names of the important files
           file1 <- x
           file2 <- file1 %>% str_replace("_R1", "_R2")
           
           file2.out <- file2 %>% str_replace ("./00_Data/01_Raw","./00_Data/02_Trimmomatic") %>% 
                                  str_replace("_R2.fastq.gz", "_R2_HC3.fastq.gz")
           
           logout <- file1 %>% str_replace("./00_Data/01_Raw","./02_Results/00_SNP_panel/02_Trimmomatic") %>% 
                               str_replace("_R1.fastq.gz", ".log")
           
           # The command
           cmd1 <- paste("-jar",
                         trimmomatic.path, "SE",
                         "-threads", 1,
                         #"-trimlog", logout %>% str_replace(".log", "_HC.log"),
                         file2, file2.out,
                         #"ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
                         "HEADCROP:3",
                         sep = " ")
           
           A1 <- system2("java", cmd1, stdout=T, stderr=T) # to R console
           A1
           
           # save a file log
           cat(file = logout,
               #"STEP1 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
               "STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
               append= F, sep = "\n\n")
           
         } ,
         mc.cores = 20
)


  
  # STEP 2 - REMOVE Illumina adaptors
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Names of the important files
             file1 <- x
             file2 <- file1 %>% str_replace("./00_Data/01_Raw","./00_Data/02_Trimmomatic") %>% 
                                str_replace("_R1.fastq.gz", "_R2_HC3.fastq.gz")
             
             fileout <- file1 %>%  str_replace("./00_Data/01_Raw","./00_Data/02_Trimmomatic") %>% 
                                   str_replace("_R1.fastq.gz", ".fastq.gz")   
             
             logout <- file1 %>% str_replace("./00_Data/01_Raw", "./02_Results/00_SNP_panel/02_Trimmomatic") %>% 
                                 str_replace("_R1.fastq.gz", ".log")
             
             # The command

             cmd1 <- paste("-jar",
                           trimmomatic.path, "PE",
                           "-threads", 1,
                           #"-trimlog", logout,
                           file1, file2,
                           "-baseout", fileout,
                           "ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
                           #"HEADCROP:1-99",
                           sep = " ")
             
             A1 <- system2("java", cmd1, stdout=T, stderr=T) # to R console
             A1
             
             # save a file log
             cat(file = logout,
                 "STEP2 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
                 #"STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
                 append= T, sep = "\n\n")
             
           } ,
           mc.cores = 20
  )
  

# SOME POST ANALYSIS FILE MANIPULATION

# Rename the files

old.name <-list.files("./00_Data/02_Trimmomatic", full.names = T)
old.name

new.name <- old.name %>% 
  str_replace("_1P.fastq", "_R1.fastq") %>% 
  str_replace("_2P.fastq", "_R2.fastq") 
  
new.name

file.rename(from = old.name,
            to = new.name)

# Then remove ALL unecessary files

file.to.remove <- list.files("./00_Data/02_Trimmomatic", full.names = T, pattern = "1U.fastq.gz|2U.fastq.gz|HC3.fastq.gz")
sum(file.size(file.to.remove))/  1024 / 1024 / 1024

file.remove(file.to.remove)

# Run FastQC - change for overwrite T first the first time to or it will not work
# THIS PART SHOULD BE BENCHMARKED TOO

fastqc(folder.in = "./00_Data/02_Trimmomatic", folder.out = "./02_Results/00_SNP_panel/01_FastQC/02_Trimmomatic", 
       overwrite = T, nthread = 20)

multiqc("./02_Results/00_SNP_panel/01_FastQC/02_Trimmomatic")


# Demultiplex -------------------------------------------------------------

# TAKE CARE, all individuals with similar names will be collapse
length(pop.info$Sample)
length(pop.info$Sample %>% unique())

#Maybe keep plate number in their name to avoid ecraser a file  

# Compter 5 heures par plaque

sub.dir <- c("NS.1413.004", "NS.1416.001", "NS.1416.002")

# Before running, check that the right barcode file is found

# If you specify the same output directory for two different
# process_radtags runs, the second run will overwrite identical filenames
# from the first run. Output the data into separate directories, and then
# concatenate the shared samples together after the runs complete

#files.to.use <-  list.files("./00_Data/02_Trimmomatic", full.names = T, pattern = "_R1") %>% 
#  str_subset("NS.1416.001")
#


for(i in sub.dir){
  
  files.to.use <-  list.files("./00_Data/02_Trimmomatic", full.names = T, pattern = "_R1") %>% 
    str_subset(i) #%>% str_subset("P02|P03")
  
  # files.to.use
  
  cat(paste("\nWorking with the run:", i),
      paste("There are" , length(files.to.use) , "plates within this run"),
      sep= "\n")
  
  # Créer un sous-dossier s'il n'existe pas
  if(file.exists(file.path("./00_Data/03a_Demultiplex", i)) == F) {
    
    cat(paste("\nCreating a new directory:", file.path("./00_Data/03a_Demultiplex", i)),
        sep= "\n")
    
    dir.create(file.path("./00_Data/03a_Demultiplex", i), recursive = T)}
  
  
# Parallel version of process_radtag 

  
mclapply(files.to.use,
         FUN = function(x){
           
           # Files
           file1 <- x
           file2 <- file1 %>% str_replace("_R1", "_R2")
           
           barcode <- list.files("./00_Data/00_FileInfos", full.names = T, pattern = "barcodes.txt") %>% 
             str_subset(x %>% str_remove("./00_Data/02_Trimmomatic") %>% 
                          str_remove("_R1.fastq.gz") %>% 
                          str_remove(i) %>%
                          str_remove("[.][A-Z][:digit:][:digit:][:digit:][.]Panomics-") %>% 
                          str_remove("/")) 
           
           plate <- barcode %>% str_remove("./00_Data/00_FileInfos") %>% 
                                str_remove("_barcodes.txt") %>% 
                                str_remove("/")
           
           # Créer un sous-dossier s'il n'existe pas
           if(file.exists(file.path("./00_Data/03a_Demultiplex", i, plate)) == F) {
             
             #cat(paste("\nCreating a new directory:", file.path("./00_Data/03a_Demultiplex", i)),
             #     sep= "\n")
             
             dir.create(file.path("./00_Data/03a_Demultiplex", i, plate), recursive = T)}
     
           # Command
           cmd <- paste("--paired",
                        "-1", file1,
                        "-2", file2,
                        "-o", file.path("./00_Data/03a_Demultiplex", i, plate),
                        "--inline_null",   
                        "-b", barcode,
                        "--renz_1", "pstI", # Check RAD site on R1 
                        #"--renz_2", "mspI", # CGG on R2 
                        
                        # Some test when adapter are not working well on one site
                        #"--adapter-1","TGCAG",
                        #"--adapter-2","GGG",
                        #"--adapter_mm",3,
                        #"--disable-rad-check",
                        
                        "-E", "phred33",
                        "--filter-illumina",
                        "-c", # clean
                        "-r", #rescue
                        "-q", #check quality
                        "-t", 135,  #truncate at 150 - 6 (restriction site) - 8 (max barcode) = 136 (all reads within sample must be the same length)
                        "-i", "gzfastq"
           )
           
           A <- system2("process_radtags", cmd, stdout=T, stderr=T)
           A
           # save a log file 
           log.file <- file1 %>% str_replace("./00_Data/02_Trimmomatic", "./02_Results/00_SNP_panel/03_Demultiplex") %>% 
             str_replace(".fastq.gz","_summary.log") %>% 
             str_remove("_R1")
           
           cat(file = log.file,
               cmd, "\n\n",
               A, # what to put in my file
               append= F, sep = "\n")
           
           # Detailed log file
           file.rename(from = list.files(file.path("./00_Data/03a_Demultiplex", i, plate), full.names = T, pattern = ".log"),
                       to = log.file %>% str_replace("summary", "detailed")
                       )
           
         } ,
         mc.cores = 20
)

gc()

}


# Extract data for each summary_detailed
# You must change the str_subset code that is the index for each project

for(x in list.files( "./02_Results/00_SNP_panel/03_Demultiplex", pattern = "detailed", full.names = T)){
  
  data <- readLines(x) %>% 
    str_subset("Pb_|PB_|Dp_|Pm_") %>% 
    str_split(pattern = "\t")
  
  cat("\n",length(data), "samples retrieved in", x %>% str_remove("./03_Results/03_Demultiplex") %>% str_remove("/"))
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  # If there is a RUN column, it's because some individuals come from more than one sample
  names(data) <- c("Barcode", "Filename", "Total", "NoRadTag", "LowQuality", "Retained")
  
  data <- data %>% mutate(Total = as.numeric(as.character(Total)),
                          NoRadTag = as.numeric(as.character(NoRadTag)),
                          LowQuality = as.numeric(as.character(LowQuality)),
                          Retained = as.numeric(as.character(Retained)),
                          Run = x %>% str_remove("./03_Results/03_Demultiplex") %>% 
                            str_remove("_log_detailed.txt")  %>% 
                            str_remove("/")
  )
  write_csv(data, x %>% str_replace("detailed.log", "Nreads.csv"))
  
}

# Create one big log file 

Nreads.data <- data.frame()

for(x in list.files( "./03_Results/03_Demultiplex", pattern = "Nreads", full.names = T) %>% str_subset(".csv") %>% str_subset("NS")){
  data.int <- read_csv(x)
  
  Nreads.data <- bind_rows(Nreads.data, data.int)
  
}

Nreads.data$Removed <- Nreads.data$Total - Nreads.data$Retained  
Nreads.data$RunSeq <- paste(sapply(str_split(Nreads.data$Run, "[.]"),`[`,2),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,3),
                            sep = ".")

Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ"))

head(Nreads.data)

write_csv(Nreads.data, file.path("./02_Results/00_SNP_panel/03_Demultiplex","AllIndividuals_Nreads.csv"))

Nreads.data <- read_csv(file.path("./02_Results/00_SNP_panel/03_Demultiplex","AllIndividuals_Nreads.csv"))

# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)

Nreads.data$Filename %>% unique() %>% length()


# Stats
Nreads.data %>% group_by(Filename) %>% summarise(Retained = sum(Retained)) %>%
                group_by() %>% 
                summarise(mean = mean(Retained)/2,
                          sd = sd(Retained)/2,
                          min = min(Retained)/2,
                          max = max(Retained)/2)
  
sum(Nreads.data$Retained)/2
sum(Nreads.data$Total)

Nreads.data %>% 
  mutate(Espece = ifelse(Espece %in% c("Pb", "PB"), "Pb", Espece )) %>% 
  group_by(Filename, Espece, Plaque, Gen_ZONE) %>% summarise(Retained = sum(Retained)) %>%
  summary()


#colMeans("Nreads.data", col = 6)

graph0.0 <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(RunSeq, Plaque_ID, Cat) %>% summarise(N = sum(N)) %>% 
  ggplot(aes(x = factor(RunSeq), y = N, fill= Cat)) + 
  geom_bar(stat= "identity") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples") +
  facet_wrap(~ Plaque_ID, scale = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
graph0.0

graph0.1 <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, Plaque_ID, Cat) %>% summarise(N = sum(N)) %>% 
  ggplot(aes(x = Filename, y = N, fill= Cat)) + 
  geom_bar(stat= "identity") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples") +
  facet_wrap(~ Plaque_ID, scale = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")
graph0.1

graph0.2 <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, Espece, Plaque_ID, Cat) %>% summarise(N = sum(N)) %>% 
  filter(Cat == "Retained") %>% 
  ggplot(aes(x = N, fill= Espece)) + 
  geom_histogram() +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(x = "N reads")+
  #facet_wrap(~ Plaque_ID, scale = "free_x") +
  theme_bw() + 
  theme(#axis.text.x = element_blank(),
        legend.position = "bottom")
graph0.2

graph0.3 <- Nreads.data %>% 
  mutate(Espece = ifelse(Espece %in% c("Pb", "PB"), "Pb", Espece )) %>% 
  group_by(Filename, Espece, Plaque, Gen_ZONE) %>% summarise(Retained = sum(Retained)) %>% 
  ggplot(aes(x = Gen_ZONE, y = Retained )) + 
  geom_boxplot() +
  geom_jitter(aes(col= factor(Plaque))) +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000),
  #                   labels = c(1:10))+
  labs(y = "N reads", x = "Samples") +
  facet_grid(. ~ Espece, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")
graph0.3

graph0.4 <- Nreads.data %>% 
  mutate(Espece = ifelse(Espece %in% c("Pb", "PB"), "Pb", Espece )) %>% 
  group_by(Espece, Plaque, Gen_ZONE, Barcode.y, Filename) %>% summarise(Retained = sum(Retained)) %>% 
  ggplot(aes(x = Barcode.y, y = Retained )) + 
  geom_violin() +
 # geom_jitter() +
 # facet_wrap(~Barcode.y)
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000),
  #                   labels = c(1:10))+
  labs(y = "N reads ", x = "Barcodes") +
 # facet_grid(. ~ Espece, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "none")
graph0.4


graph0.5 <- Nreads.data %>% 
  mutate(Espece = ifelse(Espece %in% c("Pb", "PB"), "Pb", Espece )) %>% 
  group_by(Filename, Espece, Plaque, Gen_ZONE) %>% summarise(Retained = sum(Retained)) %>%
#  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_blank())
graph0.5

pdf(file.path("./03_Results/03_Demultiplex","BasicGraph_PostDemulti.pdf"))
  print(graph0.0)
  print(graph0.1)
  print(graph0.2)
  print(graph0.3)
  print(graph0.4)
  print(graph0.5)
dev.off()


## Create one file with the 3 files by individual
sub.dir <- c("NS.1413.004", "NS.1416.001", "NS.1416.002")

plates <- file.path("./00_Data/03a_Demultiplex", sub.dir[1]) %>% list.files()

for(p in plates){

files.to.use <- list.files(file.path("./00_Data/03a_Demultiplex", sub.dir[1], p), full.names = T)

cat(paste("\nWorking with the plate:", p),
    paste("There are" , length(files.to.use) , "files to process"),
    sep= "\n")

mclapply(files.to.use,
         FUN = function(x){

           # Files
           file1 <-  x
           file2 <- file1 %>% str_replace(sub.dir[1], sub.dir[2])
           file3 <- file1 %>% str_replace(sub.dir[1], sub.dir[3])
           
           file.join <- file1 %>% str_remove(file.path(sub.dir[1], p))# %>% 
                                  #str_remove("[:digit:][:digit:][:digit:][:digit:][.][:digit:][:digit:][:digit:][.][A-Z][:digit:][:digit:][:digit:][.]Panomics-" )

           # Cat command - this is so simple !
           cmd <- paste(file1, file2, file3, ">", file.join)

           system2("cat", cmd, stdout=T, stderr=T)
           
         } ,
         mc.cores = 20
)

}


list.files(file.path("./00_Data/03a_Demultiplex"), pattern = ".fq.gz") %>% length()

# Create the list of individuals for testing parameters :

hist(Nreads.data$Retained)
#summary(Nreads.data$Retained.by.sample)
summary(Nreads.data)
Nreads.data %>% pull(Cat_Sample) %>% unique()

test.ID <- Nreads.data %>% filter(Cat_Sample == "Sample",
                                  Espece %in% "Pb") %>%  
           group_by(Filename, POP) %>% 
           summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>% 
  group_by(POP) %>% 
  sample_n(1) %>% pull(Filename)  

test.ID  # Avec 8 pop, garder 2 ind / pop

 write.table(pop.info %>% select(Sample, POP) %>% 
                          filter(Sample %in% test.ID) %>% 
                          mutate(POP = "noPOP"), 
             file = file.path("./00_Data/00_FileInfos", "popmap.test_samples.txt"),
             quote = FALSE, sep = "\t",
             row.names = F, col.names = F)

 
 # BWA - index the reference genome ----------------------------------------
 
 # Done with the newly created ref genome
 
 # Check that 
 A <- system2("bwa", "", stdout=T, stderr=T)
 
 list.files("./00_Data/07_REF_Genome")
 
 cmd <- paste("index",
              "-p",  file.path("./00_Data/07_REF_Genome","Draft2.Genome"), 
              "-a", "bwtsw",
              file.path("./00_Data/07_REF_Genome","pilon_assembly_flye_16Sequel2HiFi.renamed.fasta.scaff_s90_c5_l0_d0_e30000_r0.05.scaffolds.fasta")
 )
 
 
 A <- system2("bwa", cmd, stdout=T, stderr=T)
 
 # save a log file 
 
 cat(file = file.path("./00_Data/07_REF_Genome", "Draft2.Genome.log" ),
     cmd, "\n\n",
     A, # what to put in my file
     append= F, sep = "\n")
 
 
 # BWA - Align reads to the reference genome -------------------------------
 
 # Can start with this
  # Peut valoir la peine d'updater le code pour utiliser la fonction parallele
 # Peut-être faire rouler sur un plus petit subset au début, question de s'assurer que ça fonctionne
 # 2-3 min par séquence
 
# ATTENTION - ICI POUR PAIRED 

 demulti.files <- list.files("./00_Data/03a_Demultiplex", pattern = ".1.fq.gz", full.names = T) %>% str_subset(".rem.", negate = T) #%>%  str_subset(TEST.ID)
 
 
  mclapply(demulti.files,
          FUN = function(x){
          # How I will rename all this : 
          file.R1 <- x
          file.R2 <- x %>% str_replace(".1.fq.gz", ".2.fq.gz")
          file.bam <- x %>% str_replace("./00_Data/03a_Demultiplex", "./00_Data/03b_Align") %>% 
                            str_replace(".1.fq.gz", ".bam")
          file.sort.bam <- file.bam %>% str_replace(".bam", ".sorted.bam")
          stat.tsv <- file.bam %>% str_replace(".bam", ".stat.tsv")
          
          # DO THE ALIGMENT  
          cmd1 <- paste("mem",
                        "-t", 32,
                        "-M",

                        file.path("./00_Data/07_REF_Genome","Draft2.Genome"), # the index ref genome
                        file.R1,
                        file.R2,
                        "2> /dev/null",
                        "| samtools", "view", "-Sb", 
                        "--threads", 15, # Number of additional threads
                        "-o", file.bam#,
          )# the file
          
          A <- system2("bwa", cmd1, stdout=T, stderr=T)
          
          # Sort the file 
          cmd2 <- paste("sort",
                        "--threads", 15,
                        #"-n",
                        "-O", "BAM",
                        # the new file (sorted)
                        "-o", file.sort.bam,
                        # the bam file to sort
                        file.bam
                        
          )
          
          A <- system2("samtools", cmd2, stdout=T, stderr=T)
          
          # Compute stats
          
          cmd3 <- paste("flagstat",
                        "--threads", 15,
                        "-O", "tsv",
                        file.sort.bam,
                        ">",
                        stat.tsv
          )
          
          A <- system2("samtools", cmd3, stdout=T, stderr=T)
          
          },
          mc.cores = 2
 ) 
 

  # Remove unsorted files that are unecessary 
  
  files.to.remove <- list.files("./00_Data/03b_Align", pattern = ".bam", full.names = T)  %>% 
    str_subset(pattern = "sorted", negate = T)
  
  length(files.to.remove)
  files.to.remove[1:20]
  
  for(x in files.to.remove){
    
    file.remove(x)
  }  
   


 
 # Compute the aligned reads -----------------------------------------------
 
 map.res <- data.frame(ID = character(),
                       total = numeric(),
                       secondary = numeric(),
                       supplementary = numeric(),
                       duplicates = numeric(),
                       mapped = numeric(),
                       mapped_perc = numeric(),
                       paired = numeric(),
                       read1 = numeric(),
                       read2 = numeric(),
                       properly_paired = numeric(),
                       properly_paired_perc = numeric(),
                       twith_itself_mate_mapped = numeric(),
                       singletons = numeric(),
                       singletons_perc = numeric(),
                       twith_mate_mapped_diff_chr = numeric(),
                       twith_mate_mapped_diff_chr_HMQ = numeric(),
                       #Nmappedprim = numeric(),
                        stringsAsFactors = F)
 
 
  files.to.use <- list.files("./00_Data/03b_Align", pattern = ".stat.tsv", full.names = T)  
  
for(x in seq_along(files.to.use)){
  
   ID <-   files.to.use[x] %>% str_remove("./00_Data/03b_Align") %>% 
                               str_remove(".stat.tsv") %>% 
                               str_remove("/")
  
   temp <-  sapply(str_split(readLines(files.to.use[x]), "\t"), `[`,1)  
   map.res[x,] <- c(ID, temp %>% str_remove("%"))
   
   
}  

for(x in 2:ncol(map.res)){
  map.res[,x] <- as.numeric(as.character(map.res[,x]))
  
}  
    
nrow(map.res)
head(map.res)  
summary(map.res)  
  
 map.res <- map.res %>%  left_join(pop.data, by = c("ID" = "ID_GQ"))

 
write_csv(map.res, file = "./02_Results/00_SNP_panel/03b_Align/Mapping_res.csv")
  
map.res <- read_csv("./02_Results/00_SNP_panel/03b_Align/Mapping_res.csv")
 
graph1.0 <- map.res %>% #filter(ID != "Dp_2066") %>% 
             ggplot(aes(x = mapped_perc, fill =Espece_revision)) + 
             geom_density(alpha = 0.5) +
             #facet_grid(Espece ~ . , scales = "free")
             labs(x = "Percentage of reads mapped", y = "Density") +
             theme_bw()  
 
graph1.0

graph1.1 <-  map.res %>%  ggplot(aes(y = mapped_perc, x = Espece_revision, group = Espece_revision)) + 
   geom_boxplot(col = "black", fill = "white") +
   #geom_boxplot(alpha = 0.5) +
   geom_jitter(aes(col = Espece_revision),height = 0, width = 0.5, alpha = 0.5) +
   geom_hline(yintercept = 95, lty = "dotted") +
   #facet_grid(Espece ~ . , scales = "free")
   labs(y = "Percentage of reads mapped", x = "Species") +
   theme_bw() + 
   theme(legend.position = "none") 

graph1.1 


 
graph1.2 <- map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
   ggplot(aes(x=total, y=mapped_perc, col=Espece_revision)) +
   geom_point(alpha = 0.5)+
   labs(y = "Percentage of reads mapped", x = "N read total") +
   theme_bw() 
 
graph1.2
 
graph1.3 <- map.res %>% filter(Espece_revision %in% c("Pb"))  %>% 
   ggplot(aes(y = mapped_perc, x = Gen_ZONE)) + 
   geom_boxplot(alpha = 0.5) +
   geom_hline(yintercept = 95, lty = "dashed") +
   geom_jitter(aes(col = total),height = 0, alpha = 0.75) +
   scale_color_distiller(palette = "Spectral") +
   #facet_wrap(~ Gen_ZONE)
   labs(y = "Percentage of reads mapped", x = "SFA") +
   theme_bw()  +
   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

graph1.3 
 

graph1.4 <- map.res %>% filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped, col = Plaque_ID)) +
  geom_density() +
  geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
             lty = "dashed" )+
  #facet_wrap(~ Gen_ZONE)+
  labs(x = "N reads mapped") +
  theme_bw()
graph1.4

graph1.5 <- map.res %>% filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped)) +
  geom_density() +
  geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
             lty = "dashed" )+
  labs(x = "N reads mapped") +
  facet_wrap(~ Gen_ZONE)+
  theme_bw()
graph1.5
  
graph1.6 <- map.res %>% filter(Espece_revision %in% c("Pb"))  %>% 
            mutate(Cat_cov = ifelse(mapped_perc >= 95, "OK > 95%", "Low < 95%")) %>%    
  ggplot(aes(x = mapped, col = Cat_cov)) +
  geom_density() +
  geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
             lty = "dashed" )+
  labs(x = "N reads mapped") +
  #facet_wrap(~ Gen_ZONE)+
  theme_bw()
graph1.6

map.res %>% filter(Espece_revision %in% c("Pb"),
                   mapped_perc >= 95,
                   Cat_Sample == "Sample"
                   #mapped > 1000000
                   )  %>%
  group_by(Cat_Sample, Gen_ZONE) %>% summarise(N = n())



pdf(file.path("./02_Results/00_SNP_panel/03b_Align/","BasicGraph_PostAlign_2021-06-22.pdf"))
  print(graph1.0)
  print(graph1.1)
  print(graph1.2)
  print(graph1.3)
  print(graph1.4)
  print(graph1.5)
  print(graph1.6)
dev.off()

# Compute stat

map.res %>% filter(Espece_revision %in% c("Pb"), mapped_perc >= 80) %>% summarise(mean = mean(mapped_perc),
                                                                                         sd = sd(mapped_perc))


## Overall popmap
write.table(map.res %>% filter(Espece_revision %in% c("Pb"),
                               mapped_perc >= 95) %>% select(ID, Espece_revision), 
            file = file.path("./00_Data/00_FileInfos", "popmap_good_align.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


 map.res %>% group_by(Espece_revision) %>% 
             summarise(MedianNread = median(total),
                       MedianPerc = median(mapped_perc),
                       MinPerc = min(mapped_perc),
                       MaxPerc = max(mapped_perc))
 
 map.res %>% filter(Cat_Sample == "Sample",
                    mapped_perc < 80) %>% 
             group_by(Espece_revision, Gen_ZONE_FG) %>% 
            summarise(N = n())
 
 map.res %>% filter(Cat_Sample == "Sample",
                    between(mapped_perc, 8, 95)) %>% 
   group_by(Espece, Gen_ZONE_FG) %>% 
   summarise(N = n())
 
map.res %>% filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=singletons, col=Espece_revision)) +
  geom_point()
 

plot(map.res[,c("total", "secondary", "mapped", "paired", "properly_paired", "singletons", "twith_itself_mate_mapped", "twith_mate_mapped_diff_chr")])

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=mapped, col=Espece_revision)) +
  geom_point()
 
map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=secondary, col=Plaque_ID)) +
  geom_smooth()+
  #geom_point() +
  facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=secondary, y=singletons, col=Plaque_ID)) +
  #geom_smooth()+
  geom_point() +
  facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=singletons, col=Espece_revision)) +
  geom_point()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=twith_itself_mate_mapped, y=singletons, col=Espece_revision)) +
  geom_point()


 map.res %>% filter(Cat_Sample == "Sample",
                    between(Permapped, 0.95, 0.97)) %>% 
   group_by(Espece, Gen_ZONE) %>% 
   summarise(N = n()) %>% arrange(desc(N)) #%>% pull(N) %>% sum()
 
hist(map.res$total) 
 
quantile(map.res$mapped, c(0.01,0.05,0.1))# %>% summary()
  
# map.res %>%  %>% summarise(Mean = mean(Permapped))
 

# Stacks - Testing parameters - Ref-genome --------------------------------


cmd <- paste("--samples", "./00_Data/03b_Align",
             "--popmap", file.path("./00_Data/00_FileInfos", "popmap.test_samples.txt"),
             "-o", "./00_Data/04b_Test.ref",

             "-T", 24,
             "-X", "\"populations:-r 0.80 --genepop\"",
             "-X", "\"gstacks:-S .sorted.bam\"" # en espérant que ça fonctionne ...
)

A <- system2("ref_map.pl", cmd, stdout=T, stderr=T)
A

# Stats


# Extract statistics

  cmd1 <- paste(file.path("./00_Data/04b_Test.ref", "populations.log.distribs"),
               "samples_per_loc_postfilters")
  
  res1 <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data.int <- res1[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
  names(data.int) <- c("n_samples", "n_loci")
  
  n_loci <- as.numeric(as.character(data.int$n_loci)) %>%  sum() 
  
  # N SNP / locus
  
  cmd2 <- paste(file.path("./00_Data/04b_Test.ref", "populations.log.distribs"),
               "snps_per_loc_postfilters")
  
  res2 <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data.int <- res2[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
  names(data.int) <- c("n_snps", "n_loci")
  
  data.int$n_snps <- as.numeric(as.character(data.int$n_snps))
  data.int$n_loci <- as.numeric(as.character(data.int$n_loci))  
 
  no_div <- data.int %>% filter(n_snps == 0) %>% pull(n_loci)  %>% as.numeric()
  n_loci_poly <- n_loci - no_div

 cat("\nThere is", n_loci_poly, "polymorphic loci (r80) out of", n_loci, "loci")


# Gstakcs - Ref - discovered SNPs -----------------------------------------

cmd <- paste("-I", "./00_Data/03b_Align",
             "-M", file.path("./00_Data/00_FileInfos", "popmap_good_align.txt"),
             "-O", "./00_Data/05b_Stacks.ref",
             "-S", ".sorted.bam",
             "-t", 8)

A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path("./02_Results/00_SNP_panel/05b_Stacks.ref", "gstacks.ref.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")



# Check the distribution of ...


# N SNP / locus

cmd <- paste(file.path("./00_Data/05b_Stacks.ref", "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
names(data) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")

graph2.1 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = POP)) +
  geom_histogram() +
  geom_vline(xintercept = 10, lty = "dashed", col = "darkgray") +
  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage") + 
  theme_bw() +
  theme(legend.position = "none")

graph2.1

ggsave(filename = file.path("./02_Results/00_SNP_panel/05b_Stacks.ref", "CoverageBySFA_June2021.png"),
       plot = graph2.1,
       width = 7.5, height = 4, units = "in")

graph2.2 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Gen_ZONE))) +
  geom_point() +
  scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 80000, lty = "dashed", col = "darkgray") +
  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")

graph2.2

ggsave(filename = file.path("./02_Results/00_SNP_panel/05b_Stacks.ref", "CoverageVSloci_June2021.png"),
       plot = graph2.2,
       width = 7, height = 4, units = "in")

graph2.3 <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
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

graph2.3

ggsave(filename = file.path("./02_Results/00_SNP_panel/05b_Stacks.ref", "Coverag_June2021.png"),
       plot = graph2.3,
       width = 7, height = 4, units = "in")

write_csv(data, file.path(get.value("stacks.ref.log"),"AllIndividuals_RefGen_NreadsNloci.csv"))
