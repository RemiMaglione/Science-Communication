---
title: "2022_1"
output: html_notebook
---

###load library
```{r}
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(dplyr)
library(picante)
library(gplots)
library(RColorBrewer)
```

###loading Function
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


###import function
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


#Main
## input
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


##Build phylogenetic tree
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

##Tree metadata
```{r}
mock.2016.com.cla <- taxocomm(mock.2016.com, Mock_2016.taxo, "Class") #Class
mock.2016.com.cla.mostAbund <- mostAbundTaxa(mock.2016.com.cla, 1) #Class

mock2016.metaTree <- cbind(data.frame(id=rownames(mock2016.meta)), mock2016.meta)
both <- intersect(rownames(mock2016.metaTree), rownames(mock.2016.com.cla.mostAbund))
mock2016.metaTree <- mock2016.metaTree[both,]
mock2016.metaTree$Class <- mock.2016.com.cla.mostAbund[both,"MostAbundTaxa"]
```

## Plot tree
```{r}
colDate <- c("Early"='black', "Mid"='seagreen3', "Late"='skyblue2', "Pre"='red',
             "Rye"='#D55E00', "Squash"='#999999')

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

mock.2016.topASV.ggtree <- ggtree(mock.2016.topASV.tree, layout='fan',branch.length="none", open.angle=90) %<+% mock2016.metaTree + geom_tippoint(aes(color=Class)) +
 scale_color_manual(values=cbbPalette) 

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


### phylo analysis
```{r}
mock2016.metaTree.gp <- mock2016.metaTree %>% group_by(Class, date) %>% summarise(test = n())
mock2016.metaTree.gp.df <- as.data.frame(tidyr::pivot_wider(data = mock2016.metaTree.gp, names_from = Class, values_from = test))
rownames(mock2016.metaTree.gp.df) <- mock2016.metaTree.gp.df$date
mock2016.metaTree.gp.df <- mock2016.metaTree.gp.df[,-1]

mock2016.metaTree.gp.df[is.na(mock2016.metaTree.gp.df)] <- 0
mock2016.metaTree.gp.df
```

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

```{r}
mock.2016.topASV.sesmpd
```

```{r}
mock.2016.topASV.sesmntd
```
