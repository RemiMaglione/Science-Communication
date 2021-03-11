# load library
library(dada2)
library(ggplot2)

##### RUN1 #####
path <- "./run1" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
head(fnFs)
length(fnFs)

fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
head(fnRs)
length(fnRs)

#Extract sample names
sample.names<-as.vector(sapply(basename(fnFs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-1)],collapse="_"))}))

sample.namesR<-as.vector(sapply(basename(fnRs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-1)],collapse="_"))}))

which(!sample.names%in%sample.namesR)


## Plot quality profile of forward reads
plot12sample <- sample(1:length(fnFs), 12) 
fQual.gg<-plotQualityProfile(fnFs[plot12sample])

fQual.gg<-fQual.gg+
 geom_vline(xintercept=200)+
  geom_vline(xintercept=210)+# best
  geom_vline(xintercept=220)+
  geom_vline(xintercept=294)+
  geom_hline(yintercept = 30)

# Plot quality profile of reverse reads
rQual.gg<-plotQualityProfile(fnRs[plot12sample])
rQual.gg<-rQual.gg+
  geom_vline(xintercept=190)+
  geom_vline(xintercept=200)+# best
  geom_vline(xintercept=210)+
  geom_vline(xintercept=220)+
  geom_vline(xintercept=230)+
  geom_vline(xintercept=294)+
  geom_hline(yintercept = 30)
multiplot(fQual.gg, rQual.gg)

###qualMaxEE
#fastq.r1 <- read.csv(file="./fastqc/16S-RM2016-Mock-ML150-Q30_R1.fastq_109_fastqc/fastq.qual.csv", sep = "\t", header = TRUE)
#fastq.r2 <- read.csv(file="./fastqc/16S-RM2016-Mock-ML150-Q30_R2.fastq_109_fastqc/fastq.qual.csv", sep = "\t", header = TRUE)

#qualMaxEEplot(fastq.r1 = fastq.r1, fastq.r2 = fastq.r2)

# set filtered file folder path
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(16,19), truncLen = c(200, 200),
                     maxEE=c(2,3), 
                     maxN=0, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=6, 
                     verbose=TRUE)

filtFs <- list.files(filt_path, pattern="_F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filt_path, pattern="_R_filt.fastq.gz", full.names = TRUE)

filtFs<- sort(filtFs, method="radix")
filtRs<- sort(filtRs, method="radix")

sample.names<-as.vector(sapply(basename(filtFs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-2)],collapse="_"))}))

sample.namesR<-as.vector(sapply(basename(filtRs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-2)],collapse="_"))}))

#View(data.frame(F=sample.names, R=sample.namesR))

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

#Only for debugging purpose
#for (i in 1:nrow(as.data.frame(sample.names))) {
#  if (!(sample.names[i] %in% sample.namesR[i])) {
#    print(paste0("i=",i))
#    print(paste0("sample.names=",sample.names[i]))
#    print(paste0("sample.namesR=",sample.namesR[i]))
#  }
#}
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### BIG DATA workflow
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=5e8, multithread=3, randomize = TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=5e8, multithread=3, randomize = TRUE)
# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=3, pool = "pseudo", verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=3, pool = "pseudo", verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE, minOverlap = 30)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
saveRDS(mergers, "./mergers220_200_run1_pPool.rds")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "./seqtab220_200_run1_pPool.rds")
seqtab <- readRDS("./seqtab220_200_run1_pPool.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtabCollapse <- collapseNoMismatch(seqtab, orderBy = "abundance", minOverlap = 250, vec = TRUE, verbose = FALSE)
saveRDS(seqtabCollapse, "./seqtabCollapse220_200_run1_pPool_minOver_294.rds")
seqtabCollapse <- readRDS("./seqtabCollapse220_200_run1_pPool_minOver_294.rds")
dim(seqtabCollapse)


#### RUN2 ####
path <- "./run2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
head(fnFs)
length(fnFs)

fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
head(fnRs)
length(fnRs)

#Extract sample names
sample.names<-as.vector(sapply(basename(fnFs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-1)],collapse="_"))}))

sample.namesR<-as.vector(sapply(basename(fnRs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-1)],collapse="_"))}))

which(!sample.names%in%sample.namesR)


## Plot quality profile of forward reads
plot12sample <- sample(1:length(fnFs), 12) 
fQual.gg<-plotQualityProfile(fnFs[plot12sample])

fQual.gg<-fQual.gg+
  geom_vline(xintercept=200)+
  geom_vline(xintercept=220)+# best
  geom_vline(xintercept=220)+
  geom_vline(xintercept=294)+
  geom_hline(yintercept = 30)

# Plot quality profile of reverse reads
rQual.gg<-plotQualityProfile(fnRs[plot12sample])
rQual.gg<-rQual.gg+
  geom_vline(xintercept=190)+
  geom_vline(xintercept=200)+# best
  geom_vline(xintercept=200)+
  geom_vline(xintercept=220)+
  geom_vline(xintercept=220)+
  geom_vline(xintercept=294)+
  geom_hline(yintercept = 30)
multiplot(fQual.gg, rQual.gg)


###qualMaxEE
#fastq.r1 <- read.csv(file="./fastqc/16S-RM2016-Mock-ML150-Q30_R1.fastq_109_fastqc/fastq.qual.csv", sep = "\t", header = TRUE)
#fastq.r2 <- read.csv(file="./fastqc/16S-RM2016-Mock-ML150-Q30_R2.fastq_109_fastqc/fastq.qual.csv", sep = "\t", header = TRUE)

#qualMaxEEplot(fastq.r1 = fastq.r1, fastq.r2 = fastq.r2)

# set filtered file folder path
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(16,19), truncLen = c(220, 200),
                     maxEE=c(2,3), 
                     maxN=0, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=6, 
                     verbose=TRUE)

filtFs <- list.files(filt_path, pattern="_F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filt_path, pattern="_R_filt.fastq.gz", full.names = TRUE)

filtFs<- sort(filtFs, method="radix")
filtRs<- sort(filtRs, method="radix")

sample.names<-as.vector(sapply(basename(filtFs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-2)],collapse="_"))}))

sample.namesR<-as.vector(sapply(basename(filtRs),function(x){
  nx<-unlist(strsplit(x,split="_"))
  return(paste(nx[1:(length(nx)-2)],collapse="_"))}))

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

#Only for debugging purpose
#for (i in 1:nrow(as.data.frame(sample.names))) {
#  if (!(sample.names[i] %in% sample.namesR[i])) {
#    print(paste0("i=",i))
#    print(paste0("sample.names=",sample.names[i]))
#    print(paste0("sample.namesR=",sample.namesR[i]))
#  }
#}
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### BIG DATA workflow
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=5e8, multithread=3, randomize = TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=5e8, multithread=3, randomize = TRUE)
# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=3, pool = "pseudo", verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=3, pool = "pseudo", verbose = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE, minOverlap = 30)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
saveRDS(mergers, "./mergers220_200_run2_pPool.rds")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "./seqtab220_200_run2_pPool.rds")
seqtab <- readRDS("./seqtab220_200_run2_pPool.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtabCollapse <- collapseNoMismatch(seqtab, orderBy = "abundance", minOverlap = 250, vec = TRUE, verbose = FALSE)
saveRDS(seqtabCollapse, "./seqtabCollapse220_200_run2_pPool_minOver_294.rds")
seqtabCollapse <- readRDS("./seqtabCollapse220_200_run2_pPool_minOver_294.rds")
dim(seqtabCollapse)

### All RUN ###
#Input seqtab (from collapseNoMismatch) of each run
st1 <- readRDS("./seqtabCollapse220_200_run1_pPool_minOver_294.rds")
st2 <- readRDS("./seqtabCollapse220_200_run2_pPool_minOver_294.rds")
st.all <- mergeSequenceTables(st1, st2)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="pooled", multithread=6, verbose=TRUE)
saveRDS(seqtab.nochim, "./seqtabCollapse.220_200_2016_allRun_pPool_minOver_294_nochim.rds")
seqtab.nochim <- readRDS("./seqtabCollapse.220_200_2016_allRun_pPool_minOver_294_nochim.rds")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(st.all)

# identify taxonomy
taxa <- assignTaxonomy(seqtab.nochim,"./silva128/silva_nr_v128_train_set.fa.gz", minBoot = 50, tryRC = TRUE,
                       outputBootstraps = FALSE, 
                       taxLevels = c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species"), 
                       multithread = 6,
                       verbose = TRUE)

saveRDS(taxa, "./taxa.2016.220_200_allRun_minOver_294.rds")
taxa <- readRDS("./taxa.2016.220_200_allRun_minOver_294.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# exact species matching
taxa.sp <- addSpecies(taxa, "./silva128/silva_species_assignment_v128.fa.gz", tryRC = TRUE, allowMultiple=TRUE, verbose=FALSE)
saveRDS(taxa.sp, "./taxa.2016.sp_220_200_allRun.rds")
taxa.sp <- readRDS("./taxa.2016.sp_220_200_allRun.rds")

# inspect the taxonomic assignments
taxa.print <- taxa.sp # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
remove(taxa.print)
