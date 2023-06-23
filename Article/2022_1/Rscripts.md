---
title: "2022_1"
output: html_notebook
---

# load library
```{r lib}
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(picante)
library(gplots)
library(RColorBrewer)
library(lme4)
library(emmeans)
library(grid)
```

# loading Function
```{r makeTransparent}
makeTransparent <- function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}
```


### https://github.com/RemiMaglione/r-scripts/blob/master/scripts4Dada2/mostAbundASV.R
```{r mostAbundASV}
mostAbundASV <- function(otu_comm, topASV = 1, ASVseqColName = "MostAbundASVseq") { 
  mostAbundASV.sample <- data.frame("MostAbundASVseq"=character()) #create am empty data frame
  for (i in 1:nrow(otu_comm)) { #iterate on community_taxo matrix rows (=sample)
    tmp.row <- rownames(otu_comm[i,, drop=FALSE]) # temporary store the row names (=sample name)
    tmp.mostAbundASV <- rownames(data.frame(otu_comm[i, ])[order(-data.frame(otu_comm[i, ])),,drop=FALSE])[topASV]
    # Take only the 1rst ([1] at the end), most abundant (order(-...)) ASV (rowname) on the selected sample (i)
    # From the row of the comm taxo matrix (data.frame(otus_comm[i,])[..,,drop=FALSE"important to keep the row name during the process"]). 
    tmp.mostASVAbund <- (data.frame(otu_comm[i, ])[order(-data.frame(otu_comm[i, ])),,drop=FALSE])[topASV,]
    mostAbundASV.tmp <- data.frame("MostAbundASV"=tmp.mostAbundASV, "Abund"=tmp.mostASVAbund)
    rownames(mostAbundASV.tmp) <- tmp.row
    mostAbundASV.sample <- rbind(mostAbundASV.sample, mostAbundASV.tmp) #use rbind to add the last row to the final data frame
  }
  colnames(mostAbundASV.sample) <- c(ASVseqColName, "Abund")
  return(mostAbundASV.sample) 
}
```


```{r taxocomm}
##### taxocomm ####
#' Create taxonomy-aggregated community data.frame
#'
#' @param comm community format data.frame (samples in rows, taxa in columns, rows and columns have informative names)
#' @param taxo a data frame with rows in same order as columns of \code{comm} and taxonomic information for each row provided in columns
#' @param rank Taxonomic rank to aggregate at (a column name in \code{taxo})
#' @return A community data.frame where columns are taxonomic groups
#' @export
taxocomm <- function(comm, taxo, rank) {
  comm.taxo <- aggregate(t(comm), by=list(taxonrank=taxo[,rank]), sum)
  rownames(comm.taxo) <- comm.taxo[,1]
  return(t(comm.taxo[,-1]))
}
```


```{r mostAbundTaxa}
mostAbundTaxa <- function(otu_comm, topTaxa = 1, taxaLevel = "MostAbundTaxa") { 
  mostAbundTaxa.sample <- data.frame("MostAbundTaxa"=character()) #create am empty data frame
  for (i in 1:nrow(otu_comm)) { #iterate on community_taxo matrix rows (=sample)
    tmp.row <- rownames(otu_comm[i,, drop=FALSE]) # temporary store the row names (=sample name)
    tmp.mostAbundtaxa <- rownames(data.frame(otu_comm[i, ])[order(-data.frame(otu_comm[i, ])),,drop=FALSE])[topTaxa]
    # Take only the 1rst ([1] at the end), most abundant (order(-...)) taxa (rowname) on the selected sample (i)
    # From the row of the comm taxo matrix (data.frame(otus_comm[i,])[..,,drop=FALSE"important to keep the row name during the process"]). 
    tmp.mostTaxaAbund <- (data.frame(otu_comm[i, ])[order(-data.frame(otu_comm[i, ])),,drop=FALSE])[topTaxa,]
    mostAbundTaxa.tmp <- data.frame("MostAbundTaxa"=tmp.mostAbundtaxa, "Abund"=tmp.mostTaxaAbund)
    rownames(mostAbundTaxa.tmp) <- tmp.row
    mostAbundTaxa.sample <- rbind(mostAbundTaxa.sample, mostAbundTaxa.tmp) #use rbind to add the last row to the final data frame
  }
  colnames(mostAbundTaxa.sample) <- c(taxaLevel, "Abund")
  return(mostAbundTaxa.sample) 
}
```


```{r metaTree2df}
metaTree2df <- function(metaTree) {
  metaTree.gp <- metaTree %>% group_by(Class, date) %>% summarise(test = n())
  metaTree.gp.df <- as.data.frame(tidyr::pivot_wider(data = metaTree.gp, names_from = Class, values_from = test))
  rownames(metaTree.gp.df) <- metaTree.gp.df$date
  metaTree.gp.df <- metaTree.gp.df[,-1]
  metaTree.gp.df[is.na(metaTree.gp.df)] <- 0
  return(metaTree.gp.df)
}
```


# 2016
## input data
```{r}
mock.2016.tax <- readRDS(file="./taxa.2016.sp_220_200_mock.rds")
mock.2016.com <- readRDS(file="./seqtabCollapse.220_200_2016_mock_pPool_minOver_250_nochim.rds")
mock2016.meta <- read.csv(file = "./mock2016_metadata.csv", header = T)
```

## data processing
```{r}
rownames(mock2016.meta) <- mock2016.meta$sample
both1 <- intersect(rownames(mock2016.meta), rownames(mock.2016.com))
mock.2016.com <- mock.2016.com[both1,]
mock2016.meta <- mock2016.meta[both1,]

mock.2016.com<-mock.2016.com[which(apply(mock.2016.com, 1, sum)>0),]

abund <- apply(mock.2016.com,2,sum)
mock.2016.tax <- cbind(mock.2016.tax, abund)
taxonames <- rownames(mock.2016.tax)

# replace NAs with 'unidentified'
Mock_2016.taxo <-mock.2016.tax
Mock_2016.taxo[is.na(mock.2016.tax)] <- "unidentified"
Mock_2016.taxo <- as.data.frame(Mock_2016.taxo)
rownames(Mock_2016.taxo) <- taxonames
```


## Build phylogenetic tree
```{r}
mock.2016.com.mA_ASV <-mostAbundASV(mock.2016.com, 1) # extract most abundant ASV
mock.2016.com.mA_ASV.seq <- mock.2016.com.mA_ASV$MostAbundASVseq
names(mock.2016.com.mA_ASV.seq) <- rownames(mock.2016.com.mA_ASV)
seqList2fasta(mock.2016.com.mA_ASV.seq, file="mock.2016.topASV.fasta", write2file=TRUE, repName=FALSE, seqName = TRUE) #write seq asv as fasta

#################pynast#################
#align_seqs.py -i mock.2016.topASV.fasta
###############fastTree 2###############
#FastTree -nt mock.2016.topASV_aligned.fasta > mock.2016.topASV_aligned.tree
########################################
mock.2016.topASV.tree <- read.tree(file = "./pynast_aligned/mock.2016.topASV_aligned.tree") #import tree
```

## Tree metadata
```{r}
mock.2016.com.cla <- taxocomm(mock.2016.com, Mock_2016.taxo, "Class") #Class
mock.2016.com.cla.mostAbund <- mostAbundTaxa(mock.2016.com.cla, 1) #Class

mock2016.metaTree <- cbind(data.frame(id=rownames(mock2016.meta)), mock2016.meta)
both <- intersect(rownames(mock2016.metaTree), rownames(mock.2016.com.cla.mostAbund))
mock2016.metaTree <- mock2016.metaTree[both,]
mock2016.metaTree$Class <- mock.2016.com.cla.mostAbund[both,"MostAbundTaxa"]
```

## Plot tree
### Figure 2
```{r}
colDate <- c("Early"='black', "Mid"='seagreen3', "Late"='skyblue2', "Pre"='red',
             "Rye"='#D55E00', "Squash"='#999999')

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

mock.2016.topASV.ggtree <- ggtree(mock.2016.topASV.tree, layout='fan',branch.length="none", open.angle=90) %<+% mock2016.metaTree + geom_tippoint(aes(color=Class)) +
 scale_color_manual(values=cbbPalette) + ggtitle("2016") +
  theme(plot.title = element_text(hjust = 0.5))

mock2016.metaTree$date <- as.character(mock2016.metaTree$date)
mock2016.metaTree[which(mock2016.metaTree[,"date"]=="Date 0"),"date"] <- c("Pre")
mock2016.metaTree[which(mock2016.metaTree[,"date"]=="Date 1"),"date"] <- c("Early")
mock2016.metaTree[which(mock2016.metaTree[,"date"]=="Date 2"),"date"] <- c("Mid")
mock2016.metaTree[which(mock2016.metaTree[,"date"]=="Date 3"),"date"] <- c("Late")

mock2016.metaTree.date <- data.frame(date=mock2016.metaTree$date, row.names = rownames(mock2016.metaTree))

mock.2016.topASV.ggtree2 <- gheatmap(mock.2016.topASV.ggtree, mock2016.metaTree.date, offset = 1, color=NULL, 
         colnames_position="top", 
         colnames_angle=0, colnames_offset_y = 1, 
         hjust=0, font.size=4, width = 0.05) +
  scale_fill_manual(values=colDate)

mock2016.metaTree.material <- data.frame(material=mock2016.metaTree$material, row.names = rownames(mock2016.metaTree))

mock.2016.topASV.ggtree3 <- gheatmap(mock.2016.topASV.ggtree2, mock2016.metaTree.material, offset = 3, color=NULL, 
         colnames_position="top", 
         colnames_angle=0, colnames_offset_y = 1, 
         hjust=0, font.size=4, width = 0.05) +
  scale_fill_manual(values=colDate)

mock.2016.topASV.ggtree3
```


## phylo analysis
```{r}
mock2016.metaTree.gp <- mock2016.metaTree %>% group_by(Class, date) %>% summarise(test = n())
mock2016.metaTree.gp.df <- as.data.frame(tidyr::pivot_wider(data = mock2016.metaTree.gp, names_from = Class, values_from = test))
rownames(mock2016.metaTree.gp.df) <- mock2016.metaTree.gp.df$date
mock2016.metaTree.gp.df <- mock2016.metaTree.gp.df[,-1]

mock2016.metaTree.gp.df[is.na(mock2016.metaTree.gp.df)] <- 0
mock2016.metaTree.gp.df
```

## Standardized effect size of mean pairwise or mean nearest taxon distances in communities
```{r}
mock2016.metaTree.idda <- data.frame(id=mock2016.metaTree$id, date=mock2016.metaTree$date)
mock2016.metaTree.idda$test <- c(1)
mock2016.metaTree.idda.df <- as.data.frame(tidyr::pivot_wider(data = mock2016.metaTree.idda, names_from = id, values_from = test))
rownames(mock2016.metaTree.idda.df) <- mock2016.metaTree.idda.df$date
mock2016.metaTree.idda.df <- mock2016.metaTree.idda.df[,-1]
mock2016.metaTree.idda.df[is.na(mock2016.metaTree.idda.df)] <- 0

mock.2016.topASV.tree.dist <- cophenetic(mock.2016.topASV.tree)

mock.2016.topASV.sesmpd <- ses.mpd(mock2016.metaTree.idda.df, mock.2016.topASV.tree.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
mock.2016.topASV.sesmntd <- ses.mntd(mock2016.metaTree.idda.df, mock.2016.topASV.tree.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
```

## output results of mean pairwise distances
```{r}
mock.2016.topASV.sesmpd
```
##  output results of mean nearest taxon distances
```{r}
mock.2016.topASV.sesmntd
```

# 2017
## Input data
```{r}
mock.2017.tax <- read.table(file = "./mock2017.taxo", header = T)
mock2017.meta <- read.csv(file = "./mock2017.meta.csv", header = T)
rownames(mock2017.meta)<-mock2017.meta$id
rownames(mock.2017.tax) <-mock.2017.tax$sample

mock2017.meta <- mock2017.meta[which(!mock2017.meta$material%in%c("Rainwater", "Air")),]

both <- intersect(rownames(mock2017.meta), rownames(mock.2017.tax))

mock.2017.tax <- mock.2017.tax[both,]

mock.2017.tax.seq <- mock.2017.tax$seq

names(mock.2017.tax.seq) <- mock.2017.tax$sample
seqList2fasta(mock.2017.tax.seq, file="./2017/mock.2017.fasta", write2file=TRUE, repName=FALSE, seqName = TRUE) #write seq asv as fasta
#################pynast#################
#align_seqs.py -i mock.2017.fasta
###############fastTree 2###############
#cd pynast_aligned/
#FastTree -nt mock.2017_aligned.fasta > mock.2017_aligned.tree
########################################
mock.2017.tree <- read.tree(file = "./2017/pynast_aligned/mock.2017_aligned.tree") #import tree
```

```{r}
both2 <- intersect(rownames(mock2017.meta), mock.2017.tree$tip.label)
mock2017.metaTree <- mock2017.meta[both2,]

mock2017.merge <- merge(mock2017.metaTree, mock.2017.tax, by = 0)
mock2017.metaTree <- data.frame(mock2017.merge[,-1], row.names = mock2017.merge$id)
```


## Standardized effect size of mean pairwise or mean nearest taxon distances in communities
```{r}
mock2017.metaTree.idda <- data.frame(id=mock2017.metaTree$id, date=mock2017.metaTree$date)
mock2017.metaTree.idda$test <- c(1)
mock2017.metaTree.idda.df <- as.data.frame(tidyr::pivot_wider(data = mock2017.metaTree.idda, names_from = id, values_from = test))
rownames(mock2017.metaTree.idda.df) <- mock2017.metaTree.idda.df$date
mock2017.metaTree.idda.df <- mock2017.metaTree.idda.df[,-1]
mock2017.metaTree.idda.df[is.na(mock2017.metaTree.idda.df)] <- 0

mock2017.metaTree.dist <- cophenetic(mock.2017.tree)

mock.2017.sesmpd <- ses.mpd(mock2017.metaTree.idda.df, mock2017.metaTree.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
mock.2017.sesmntd <- ses.mntd(mock2017.metaTree.idda.df, mock2017.metaTree.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
```

## output results of mean pairwise distances
```{r}
mock.2017.sesmpd
```


##  output results of mean nearest taxon distances
```{r}
mock.2017.sesmntd
```

## Plot tree
### Figure 2
```{r}
colDate <- c("Early"='black', "Mid"='seagreen3', "Late"='skyblue2', "Pre"='red',
             "Rye"='#D55E00', "Squash"='#999999')
             #"Rainwater"='darkorange', "Air"='slateblue1')

mock2017.metaTree$date <- as.character(mock2017.metaTree$date)
mock2017.metaTree[which(mock2017.metaTree[,"date"]=="Date 0"),"date"] <- c("Pre")
mock2017.metaTree[which(mock2017.metaTree[,"date"]=="Date 1"),"date"] <- c("Early")
mock2017.metaTree[which(mock2017.metaTree[,"date"]=="Date 2"),"date"] <- c("Mid")
mock2017.metaTree[which(mock2017.metaTree[,"date"]=="Date 3"),"date"] <- c("Late")

mock.2017.ggtree <- ggtree(mock.2017.tree,  layout='fan', branch.length="none", open.angle=90) %<+% mock2017.metaTree + geom_tippoint(aes(color=Class)) +
 scale_color_manual(values=cbbPalette) + ggtitle("2017") +
  theme(plot.title = element_text(hjust = 0.5))

mock2017.metaTree.date <- data.frame(date=mock2017.metaTree$date, row.names = rownames(mock2017.metaTree))

mock.2017.ggtree2 <- gheatmap(mock.2017.ggtree, mock2017.metaTree.date, offset = 1, color=NULL, 
         colnames_position="top", 
         colnames_angle=0, colnames_offset_y = 1, 
         hjust=0, font.size=4, width = 0.05) +
  scale_fill_manual(values=colDate)

mock2017.metaTree.material <- data.frame(material=mock2017.metaTree$material, row.names = rownames(mock2017.metaTree))

mock.2017.ggtree3 <- gheatmap(mock.2017.ggtree2, mock2017.metaTree.material, offset = 4, color=NULL, 
         colnames_position="top", 
         colnames_angle=0, colnames_offset_y = 1, 
         hjust=0, font.size=4, width = 0.05) +
  scale_fill_manual(values=colDate)

mock.2017.ggtree3
```

## Figure 2 : final plot
```{r}
multiplot(mock.2016.topASV.ggtree3, mock.2017.ggtree3, cols = 2)
```

## Chi-Square test
### on class count for all sampling date
#### 2016
```{r}
mock2016.metaTree
mock2016.metaTree.df <-  metaTree2df(mock2016.metaTree)
chisq.test(apply(mock2016.metaTree.df, 2, sum))
```

#### 2017
```{r}
mock2017.metaTree.df <-  metaTree2df(mock2017.metaTree)
chisq.test(apply(mock2017.metaTree.df, 2, sum))
```

## Per sample type (squash leaves or rye material)
### 2016
```{r}
mock2016.metaTree.rye <- mock2016.metaTree[which(mock2016.metaTree$material=="Rye"),]
mock2016.metaTree.squ <- mock2016.metaTree[which(mock2016.metaTree$material=="Squash"),]

mock2016.metaTree.rye.df <-metaTree2df(mock2016.metaTree.rye)
mock2016.metaTree.squ.df <-metaTree2df(mock2016.metaTree.squ)
mock2016.metaTree.rye.df <- mock2016.metaTree.rye.df[c("Pre", "Early", "Mid", "Late"),]
mock2016.metaTree.squ.df <- mock2016.metaTree.squ.df[c("Early", "Mid", "Late"),]
```

### 2017
```{r}
mock2017.metaTree.rye <- mock2017.metaTree[which(mock2017.metaTree$material=="Rye"),]
mock2017.metaTree.squ <- mock2017.metaTree[which(mock2017.metaTree$material=="Squash"),]

mock2017.metaTree.rye.df <-metaTree2df(mock2017.metaTree.rye)
mock2017.metaTree.squ.df <-metaTree2df(mock2017.metaTree.squ)
mock2017.metaTree.rye.df <- mock2017.metaTree.rye.df[c("Pre", "Early", "Mid", "Late"),]
mock2017.metaTree.squ.df <- mock2017.metaTree.squ.df[c("Early", "Mid", "Late"),]
```

## Gammaproteobacteria count between sampling date
### Rye:2016&2017
```{r}
chisq.test(mock2017.metaTree.rye.df[,"Gammaproteobacteria"]+mock2016.metaTree.rye.df[,"Gammaproteobacteria"])
#chisq.test(mockAllYear.metaTree.rye.df.rel.final[,"Gammaproteobacteria"])
```

### Squash:2016&2017
```{r}
chisq.test(mock2017.metaTree.squ.df[,"Gammaproteobacteria"]+mock2016.metaTree.squ.df[,"Gammaproteobacteria"])
#chisq.test(mockAllYear.metaTree.squ.df.rel.final[,"Gammaproteobacteria"])
```

## Actinobacteria count between sampling date
### Rye:2016&2017
```{r}
chisq.test(mock2016.metaTree.rye.df[,"Actinobacteria"]+mock2017.metaTree.rye.df[,"Actinobacteria"])
```

### Squash: 2016&2017
```{r}
chisq.test(mock2016.metaTree.squ.df[,"Actinobacteria"]+mock2017.metaTree.squ.df[,"Actinobacteria"])
```

## Alphaproteobacteria count between sampling date
### Rye:2016&2017
```{r}
chisq.test(mock2016.metaTree.rye.df[,"Alphaproteobacteria"]+mock2017.metaTree.rye.df[,"Alphaproteobacteria"])
```

### Squash: 2016&2017
```{r}
chisq.test(mock2016.metaTree.squ.df[,"Alphaproteobacteria"]+mock2017.metaTree.squ.df[,"Alphaproteobacteria"])
```

# Figure 3
```{r}
MyPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#000000")
MyPalette2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#C0C0C0", "#0072B2")

mockAllYear.metaTree.squ.df <- bind_rows(mock2016.metaTree.squ.df %>% add_rownames("Season_Sampling"), 
          mock2017.metaTree.squ.df %>% add_rownames("Season_Sampling")) %>% 
    group_by(Season_Sampling)  %>% 
    replace(is.na(.), 0) %>%
    summarise_all(sum)

mockAllYear.metaTree.squ.df <- as.data.frame(mockAllYear.metaTree.squ.df)
rownames(mockAllYear.metaTree.squ.df) <- mockAllYear.metaTree.squ.df$Season_Sampling
mockAllYear.metaTree.squ.df <- mockAllYear.metaTree.squ.df[c("Early", "Mid", "Late"),]

mockAllYear.metaTree.rye.df <- bind_rows(mock2016.metaTree.rye.df %>% add_rownames("Season_Sampling"), 
          mock2017.metaTree.rye.df %>% add_rownames("Season_Sampling")) %>% 
    group_by(Season_Sampling)  %>% 
    replace(is.na(.), 0) %>% 
    summarise_all(sum)

mockAllYear.metaTree.rye.df <- as.data.frame(mockAllYear.metaTree.rye.df)
rownames(mockAllYear.metaTree.rye.df) <- mockAllYear.metaTree.rye.df$Season_Sampling
mockAllYear.metaTree.rye.df <- mockAllYear.metaTree.rye.df[c("Pre", "Early", "Mid", "Late"),]


mockAllYear.metaTree.squ.ggdf <- mockAllYear.metaTree.squ.df %>% tidyr::pivot_longer(!Season_Sampling, names_to = "class", values_to = "count") %>% mutate(material="Squash") %>%  group_by(Season_Sampling) %>% 
  mutate(count_rel = count / sum(count) * 100)

mockAllYear.metaTree.squ.df.rel <- mockAllYear.metaTree.squ.ggdf[,c("Season_Sampling","class", "count_rel")]
mockAllYear.metaTree.squ.df.rel.final <- mockAllYear.metaTree.squ.df.rel %>% tidyr::pivot_wider(names_from = "class", values_from = "count_rel") %>% summarise_all(sum)

mockAllYear.metaTree.rye.ggdf <- mockAllYear.metaTree.rye.df %>% tidyr::pivot_longer(!Season_Sampling, names_to = "class", values_to = "count") %>% mutate(material="Rye") %>%  group_by(Season_Sampling) %>% 
  mutate(count_rel = count / sum(count) * 100)

mockAllYear.metaTree.rye.df.rel <- mockAllYear.metaTree.rye.ggdf[,c("Season_Sampling","class", "count_rel")]
mockAllYear.metaTree.rye.df.rel.final <- mockAllYear.metaTree.rye.df.rel %>% tidyr::pivot_wider(names_from = "class", values_from = "count_rel") %>% summarise_all(sum)

mockAllYear.metaTree.squ.ggdf.gg<- ggplot(mockAllYear.metaTree.squ.ggdf, aes(x=Season_Sampling, y=count_rel, fill=class)) + geom_bar(stat = "identity") + ggtitle("Squash") +
 scale_fill_manual(values=MyPalette2) + ylab("Relative class abundance (%)") + xlab("Sampling Season") + theme_light() + scale_x_discrete(limits=c( "Early", "Mid", "Late"))  + labs(fill='Taxonomic class')

mockAllYear.metaTree.rye.ggdf.gg <- ggplot(mockAllYear.metaTree.rye.ggdf, aes(x=Season_Sampling, y=count_rel, fill=class)) + geom_bar(stat = "identity") + ggtitle("Rye") +
 scale_fill_manual(values=MyPalette) + ylab("Relative class abundance (%)") + xlab("Sampling Season") + theme_light() + scale_x_discrete(limits=c("Pre", "Early", "Mid", "Late")) + labs(fill='Taxonomic class')

multiplot(mockAllYear.metaTree.rye.ggdf.gg, mockAllYear.metaTree.squ.ggdf.gg, cols = 1) 
```

# Figure 4
## Data input
```{r}
comp19.dish <- read.csv(file = "./Inhibition de croissance des bactéries phytopathogènes Champion_2016_2017 (last).csv", sep = ",", header = T, na.strings = "NA")
```

## Data processing
```{r}
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
antago.list <- c("MC16-285","MP17-005", "MP17-115", "MP17-308")
antago.df <- data.frame(id=antago.list, col=cbbPalette[1:4])

for (i in 1:nrow(comp19.dish)){
  ifelse(comp19.dish[i, "id"]%in%antago.df$id,
         yes= comp19.dish[i, "col"] <- as.character(antago.df[which(antago.df$id==as.character(comp19.dish[i, "id"])), "col"]),
         no= comp19.dish[i, "col"] <- "grey90")
}
df3 <- merge(comp19.dish, antago.df, by="id")
```

## ANOVA
```{r anova}
MC16_285.aov <- anova(lm(P.s.courge.Per~xdate, data = df3[which(df3$id=="MC16-285"),]))
MP17_005.aov <- anova(lm(P.s.courge.Per~xdate, data = df3[which(df3$id=="MP17-005"),]))
MP17_115.aov <- anova(lm(P.s.courge.Per~date, data = df3[which(df3$id=="MP17-115"),]))
MP17_308.aov <- anova(lm(P.s.courge.Per~xdate, data = df3[which(df3$id=="MP17-308"),]))

#Manually export ANOVA results
df3$aovPval <- NA
df3[which(df3$id=="MC16-285"),"aovPval"] <- "p<0.0001"
df3[which(df3$id=="MP17-005"),"aovPval"] <- "p=0.0136"
df3[which(df3$id=="MP17-115"),"aovPval"] <- "p=0.9294"
df3[which(df3$id=="MP17-308"),"aovPval"] <- "p<0.0001"
```

## plot data and linear regression
```{r}
figure4.gg<-ggplot(df3, aes(x=xdate, y=P.s.courge.Per))+ 
  #geom_line(data=df2, colour="grey90") +
  geom_point(data=df3) +   geom_smooth(method="lm", se=T) +
  #geom_point(data=df3)
  #scale_color_manual(values = cbbPalette)  + 
  facet_grid(id~.) + theme_bw() + ylim(0,0.75) +
  #scale_x_discrete(labels=c("1d", "2d", "3d", "4d", "7d")) +
  xlab(expression(paste("Day count after ",italic("P. syringae"), " inoculation"))) + ylab("Pathogen inhibition (%)") + scale_x_continuous(labels = as.character(df3$xdate), breaks = df3$xdate) + geom_text(mapping =aes(x=6.5, y=0.65, label=aovPval))

svg(filename = "figure_4.svg")
figure4.gg
dev.off()

```


## Figure 5
```{r}
y_title <- expression(paste(italic("P. syringae"), " symptoms count"))

ggplot(data = test1.pseudo, aes(x=antagonist, y=pseudo_spot, fill=antagonist)) + geom_jitter(aes(color=antagonist), ) + geom_violin(alpha=0.3) +  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="red", linetype = 1, show.legend = FALSE) + theme_bw() + ylab(y_title) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_fill_discrete(name="Antagonist", breaks=c("EAU", "MC16-285", "MP17-005", "MP17-115", "MP17-308", "MIX"), labels=c("Water", "MC16-285", "MP17-005", "MP17-115", "MP17-308", "Mix")) + scale_color_discrete(name="Antagonist", breaks=c("EAU", "MC16-285", "MP17-005", "MP17-115", "MP17-308", "MIX"), labels=c("Water", "MC16-285", "MP17-005", "MP17-115", "MP17-308", "Mix")) 
```


# Figure 6
## Data input and processing
```{r}
test1 <- read.table(file = "./test1.csv", sep = ",", header = T)
test1.pseudo <- test1[which(test1[,"treatment"]=="B17-149"),]
test1.pseudo$rep <- factor(test1.pseudo$rep)
test1.pseudo$antagonist <- factor(test1.pseudo$antagonist)

test1.pseudo$antagonist <- factor(test1.pseudo$antagonist, levels = c("EAU", "MC16-285", "MP17-005", "MP17-115", "MP17-308", "MIX"))

test1.pseudo[test1.pseudo$antagonist=="EAU", "test"] <- c("neg")
test1.pseudo[test1.pseudo$antagonist=="MIX", "test"] <- c("mix")
test1.pseudo[!test1.pseudo$antagonist%in%c("EAU", "MIX"), "test"] <- c("antagonist")

test1.pseudo.wd <- test1.pseudo %>% pivot_longer(cols = c("Feuille.A","Feuille.B","Feuille.C","Feuille.D"), names_to = "leaf", values_to = "pseudo_count" )

test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.A", "Leaf"] <- "A"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.B", "Leaf"] <- "B"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.C", "Leaf"] <- "C"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.D", "Leaf"] <- "D"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.A", "Leafx"] <- "4"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.B", "Leafx"] <- "3"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.C", "Leafx"] <- "2"
test1.pseudo.wd[test1.pseudo.wd$leaf=="Feuille.D", "Leafx"] <- "1"

test1.pseudo.wd[test1.pseudo.wd$antagonist=="EAU"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "ab"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="EAU"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="EAU"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "bc"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="EAU"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "C"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MIX"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MIX"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "ab"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MIX"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "b"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MIX"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "c"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MC16-285"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MC16-285"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MC16-285"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MC16-285"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "b"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-005"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-005"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "ab"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-005"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "abc"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-005"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "c"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-115"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-115"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-115"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "b"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-115"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "b"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-308"&test1.pseudo.wd$leaf=="Feuille.A","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-308"&test1.pseudo.wd$leaf=="Feuille.B","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-308"&test1.pseudo.wd$leaf=="Feuille.C","tukeyHSD_glmer"] <- "a"
test1.pseudo.wd[test1.pseudo.wd$antagonist=="MP17-308"&test1.pseudo.wd$leaf=="Feuille.D","tukeyHSD_glmer"] <- "a"

test1.pseudo.wd$antagonist2 <- test1.pseudo.wd$antagonist
test1.pseudo.wd$antagonist2 <- as.character(test1.pseudo.wd$antagonist)
test1.pseudo.wd[test1.pseudo.wd$antagonist2=="EAU","antagonist2"] <- "Water"
```

## Normalization for figure 6.B
```{r}
test1.pseudo.wd.mean <-test1.pseudo.wd %>% group_by(antagonist, Leaf) %>% mutate(pseudo_count_mean = mean(pseudo_count))
leafA.mean <- test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$antagonist2=="Water"&test1.pseudo.wd.mean$Leaf=="A"),"pseudo_count_mean"][[1]][1]
leafB.mean <- test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$antagonist2=="Water"&test1.pseudo.wd.mean$Leaf=="B"),"pseudo_count_mean"][[1]][1]
leafC.mean <- test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$antagonist2=="Water"&test1.pseudo.wd.mean$Leaf=="C"),"pseudo_count_mean"][[1]][1]
leafD.mean <- test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$antagonist2=="Water"&test1.pseudo.wd.mean$Leaf=="D"),"pseudo_count_mean"][[1]][1]
test1.pseudo.wd.mean$pseudo_count_normPlus <- NA
test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="A"),"pseudo_count_normPlus"] <- (test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="A"),"pseudo_count"]+1)/(leafA.mean+1)
test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="B"),"pseudo_count_normPlus"] <- (test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="B"),"pseudo_count"]+1)/(leafB.mean+1)
test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="C"),"pseudo_count_normPlus"] <- (test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="C"),"pseudo_count"]+1)/(leafC.mean+1)
test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="D"),"pseudo_count_normPlus"] <- (test1.pseudo.wd.mean[which(test1.pseudo.wd.mean$Leaf=="D"),"pseudo_count"]+1)/(leafD.mean+1)
```

## Linear mixed modeling
```{r}
pseudo_count_normPlus.lmer2 <- lmerTest::lmer(log2(pseudo_count_normPlus) ~ 0 + antagonist * leaf + (1|rep_row) + (1|rep_col) , data = test1.pseudo.wd.mean)
pseudo_count_normPlus.lmer.emm <- emmeans(pseudo_count_normPlus.lmer2, list(pairwise ~ antagonist * leaf), adjust = "tukey")

myList <- c("EAU Feuille.A - (MC16-285 Feuille.A)",
"EAU Feuille.A - (MP17-005 Feuille.A)" ,
"EAU Feuille.A - (MP17-115 Feuille.A)" ,
"EAU Feuille.A - (MP17-308 Feuille.A)",
"EAU Feuille.A - MIX Feuille.A",
"EAU Feuille.B - (MC16-285 Feuille.B)",
"EAU Feuille.B - (MP17-005 Feuille.B)" ,
"EAU Feuille.B - (MP17-115 Feuille.B)" ,
"EAU Feuille.B - (MP17-308 Feuille.B)",
"EAU Feuille.B - MIX Feuille.B",
"EAU Feuille.C - (MC16-285 Feuille.C)",
"EAU Feuille.C - (MP17-005 Feuille.C)" ,
"EAU Feuille.C - (MP17-115 Feuille.C)" ,
"EAU Feuille.C - (MP17-308 Feuille.C)",
"EAU Feuille.C - MIX Feuille.C",
"EAU Feuille.D - (MC16-285 Feuille.D)",
"EAU Feuille.D - (MP17-005 Feuille.D)" ,
"EAU Feuille.D - (MP17-115 Feuille.D)" ,
"EAU Feuille.D - (MP17-308 Feuille.D)",
"EAU Feuille.D - MIX Feuille.D")

pseudo_count_normPlus.lmer.emm$`pairwise differences of antagonist, leaf`[which(pseudo_count_normPlus.lmer.emm$`pairwise differences of antagonist, leaf`@levels$`1`%in%myList)]

test1.pseudo.wd.mean$normPval <- NA
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285"&test1.pseudo.wd.mean$leaf=="Feuille.A","normPval"] <- "**"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285"&test1.pseudo.wd.mean$leaf=="Feuille.B","normPval"] <- "**"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285"&test1.pseudo.wd.mean$leaf=="Feuille.C","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285"&test1.pseudo.wd.mean$leaf=="Feuille.D","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005"&test1.pseudo.wd.mean$leaf=="Feuille.A","normPval"] <- "***"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005"&test1.pseudo.wd.mean$leaf=="Feuille.B","normPval"] <- "***"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005"&test1.pseudo.wd.mean$leaf=="Feuille.C","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005"&test1.pseudo.wd.mean$leaf=="Feuille.D","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115"&test1.pseudo.wd.mean$leaf=="Feuille.A","normPval"] <- "*"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115"&test1.pseudo.wd.mean$leaf=="Feuille.B","normPval"] <- "***"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115"&test1.pseudo.wd.mean$leaf=="Feuille.C","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115"&test1.pseudo.wd.mean$leaf=="Feuille.D","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308"&test1.pseudo.wd.mean$leaf=="Feuille.A","normPval"] <- "#"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308"&test1.pseudo.wd.mean$leaf=="Feuille.B","normPval"] <- "***"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308"&test1.pseudo.wd.mean$leaf=="Feuille.C","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308"&test1.pseudo.wd.mean$leaf=="Feuille.D","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX"&test1.pseudo.wd.mean$leaf=="Feuille.A","normPval"] <- "#"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX"&test1.pseudo.wd.mean$leaf=="Feuille.B","normPval"] <- "**"
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX"&test1.pseudo.wd.mean$leaf=="Feuille.C","normPval"] <- ""
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX"&test1.pseudo.wd.mean$leaf=="Feuille.D","normPval"] <- ""

test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285","normPvalPos"] <- max(log2(test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MC16-285","pseudo_count_normPlus"]))
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005","normPvalPos"] <- max(log2(test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-005","pseudo_count_normPlus"]))
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115","normPvalPos"] <- max(log2(test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-115","pseudo_count_normPlus"]))
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308","normPvalPos"] <- max(log2(test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MP17-308","pseudo_count_normPlus"]))
test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX","normPvalPos"] <- max(log2(test1.pseudo.wd.mean[test1.pseudo.wd.mean$antagonist=="MIX","pseudo_count_normPlus"]))
```

### Figure 6.A
```{r}
leaf.gg4 <- test1.pseudo.wd %>%
ggplot(aes(x=Leafx, y=log10(pseudo_count))) + geom_violin() + geom_jitter(alpha=0.25) +
facet_grid(antagonist2~.) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, col="grey30", linetype = 1) +
theme_bw() +    ylim(0,max(log10(test1.pseudo.wd$pseudo_count))+0.1) + geom_text(mapping =aes(y=max(log10(test1.pseudo.wd$pseudo_count))+0.05, label=tukeyHSD_glmer)) + ylab(y_title3)  + xlab( "Leaf from the top at inoculation (young to old)") + ggtitle("A.")
```

### Figure 6.B
```{r}
y_title6 <- expression(paste("Normalized ",italic("P. syringae"), " symptoms count log10(relative to mean(Water))"))
leaf.gg6c <- test1.pseudo.wd.mean[which(!test1.pseudo.wd.mean$antagonist2=="Water"),] %>%
ggplot(aes(x=Leafx, y=log2(pseudo_count_normPlus))) + geom_violin() + geom_jitter(alpha=0.25) +
facet_grid(antagonist2~., scales ="free") +
stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, col="grey30", linetype = 1) +
theme_bw()+ ylab(y_title6) + geom_hline(yintercept = 0, color="black") + xlab( "Leaf from the top at inoculation (young to old)") + ggtitle("B.") +
geom_text(mapping =aes(y=normPvalPos+0.1, label=normPval))
```

## Plotting final figure 6
```{r}
multiplot(leaf.gg4, leaf.gg6c, cols = 2)
```


# Supplemental Figure 5
```{r}
y_titleS5A <- expression(paste("Day count after pathogen inoculation"))

gg5SA <- ggplot(df3, aes(x=xdate, y=X.h.laitue.Per))+ 
  #geom_line(data=df2, colour="grey90") +
  geom_point(data=df3) +   geom_smooth(method="lm", se=T) +
  #geom_point(data=df3)
  #scale_color_manual(values = cbbPalette)  + 
  facet_grid(id~.) + theme_bw() + ylim(0,0.75) + ggtitle("A") + 
  #scale_x_discrete(labels=c("1d", "2d", "3d", "4d", "7d")) +
  xlab(y_titleS5A) + ylab("Pathogen inhibition (%)") + scale_x_continuous(labels = as.character(df3$xdate), breaks = df3$xdate)

y_titleS5B <- expression(paste("Day count after pathogen inoculation"))

gg5SB <- ggplot(df3, aes(x=xdate, y=X.c.chou.Per))+ 
  #geom_line(data=df2, colour="grey90") +
  geom_point(data=df3) +   geom_smooth(method="lm", se=T) +
  #geom_point(data=df3)
  #scale_color_manual(values = cbbPalette)  + 
  facet_grid(id~.) + theme_bw() + ylim(0,0.75) + ggtitle("B") + 
  #scale_x_discrete(labels=c("1d", "2d", "3d", "4d", "7d")) +
  xlab(y_titleS5B) + ylab("Pathogen inhibition (%)") + scale_x_continuous(labels = as.character(df3$xdate), breaks = df3$xdate)


y_titleS5C <- expression(paste("Day count after pathogen inoculation"))

gg5SC <- ggplot(df3, aes(x=xdate, y=P.s.haricot.Per))+ 
  #geom_line(data=df2, colour="grey90") +
  geom_point(data=df3) +   geom_smooth(method="lm", se=T) +
  #geom_point(data=df3)
  #scale_color_manual(values = cbbPalette)  + 
  facet_grid(id~.) + theme_bw() + ylim(0,0.75) + ggtitle("C") + 
  #scale_x_discrete(labels=c("1d", "2d", "3d", "4d", "7d")) +
  xlab(y_titleS5C) + ylab("Pathogen inhibition (%)") + scale_x_continuous(labels = as.character(df3$xdate), breaks = df3$xdate)
gg5SC

svg(filename = "supplemental_figure_5.svg", height = 7.5)
multiplot(gg5SA, gg5SC, gg5SB, cols = 2)
dev.off()
```


# Supplemental Figure 6
```{r}
test1.dry <- test1.pseudo[,c("antagonist", "plant_dry_weight")]
test1.dry$antagonist2 <- as.character(test1.dry$antagonist)
test1.dry[which(test1.dry$antagonist2=="EAU"), "antagonist2"] <- "Water"
test1.dry$antagonist2 <- as.factor(test1.dry$antagonist2)

ggplot(data = test1.dry, aes(x=antagonist2, y=plant_dry_weight)) + geom_jitter() + geom_violin(alpha=0.3) +  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="black", linetype = 1) + xlab("Antagonist") + ylab("Plant dry weight (g)") + ylim (0, max(test1$plant_dry_weight)) 
```


# Supplemental Figure 7
```{r}
chamTrt <- read.table(file="champion_treatment.tsv.txt", sep = "\t", header=T)
chamTrt.df <- chamTrt %>% group_by(year) %>%count(treatment)

chamTrt.ddf <- chamTrt.df %>% group_by(year) %>% mutate(count_rel = n / sum(n) * 100)

chamTrt.ddf$year <- as.factor(chamTrt.ddf$year)
colvec <- c("#FFA500", "#800000", "#0000FF", "#006400")
  
chamTrt.df.gg <- ggplot(chamTrt.ddf, aes(x=year, y=count_rel, fill=treatment, label=count_rel)) + geom_bar(stat = "identity") + ggtitle("A")+
 scale_fill_manual(values=colvec) + ylab("Relative strain isolation frequency") + xlab("Sampling year") + theme_light() + labs(fill="Treatment") + geom_text(size = 4, position = position_stack(vjust = 0.5), color="#FFFFFFFF", fontface ="bold")
chamTrt.df.gg


chamTrt[chamTrt$treatment=="RCC", "trt.mod"] <- "Cover Crop"
chamTrt[chamTrt$treatment=="CT-RCC", "trt.mod"] <- "Cover Crop"
chamTrt[chamTrt$treatment=="PC", "trt.mod"] <- "Non Cover Crop"
chamTrt[chamTrt$treatment=="BS", "trt.mod"] <- "Non Cover Crop"
chamTrt$year <- as.factor(chamTrt$year)

chamTrt.df2 <- chamTrt %>% group_by(year) %>%count(trt.mod)
chamTrt.ddf2 <- chamTrt.df2 %>% group_by(year) %>% mutate(count_rel = n / sum(n) * 100)

chamTrt.ddf2.gg <- ggplot(chamTrt.ddf2, aes(x=year, y=count_rel, fill=trt.mod, label=count_rel)) + geom_bar(stat = "identity") + ggtitle("B")+
 scale_fill_manual(values=MyPalette[2:1]) + ylab("Relative strain isolation frequency") + xlab("Sampling year") + theme_light() + labs(fill="Treatment category") + geom_text(size = 4, position = position_stack(vjust = 0.5), color="#FFFFFFFF", fontface ="bold")
chamTrt.ddf2.gg

svg(filename = "supplmental_figure_7.svg", height = 5, width = 7)
multiplot(chamTrt.df.gg, chamTrt.ddf2.gg, cols = 2)
dev.off()
```

```{r}
chisq.test(chamTrt.ddf2[which(chamTrt.ddf2$year==2016),"n"])
```

```{r}
chisq.test(chamTrt.ddf[which(as.data.frame(chamTrt.ddf[which(chamTrt.ddf$year==2016),"treatment"])$treatment%in%c("CT-RCC", "RCC")),"n"])
```

```{r}
chisq.test(chamTrt.ddf[which(chamTrt.ddf$year==2017),"n"])
```
