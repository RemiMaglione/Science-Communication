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

## 2017
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

### Plot tree
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
 scale_color_manual(values=cbbPalette) 

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

```{r}
metaTree2tdyLonger <- function(metaTree) {
metaTree.df.tdy <- metaTree %>% group_by(Class, date) %>% summarise(test = n())
metaTree.df.tdy.wd <- tidyr::pivot_wider(data = metaTree.df.tdy, names_from = Class, values_from = test)
metaTree.df.tdy.lg <- metaTree.df.tdy.wd %>% tidyr::pivot_longer(!date, names_to = "class", values_to = "count")
metaTree.df.tdy.lg[is.na(metaTree.df.tdy.lg)] <- 0
return(metaTree.df.tdy.lg)}

mock2017.metaTree.rye.tlg <- metaTree2tdyLonger(mock2017.metaTree.rye)
mock2017.metaTree.squ.tlg <- metaTree2tdyLonger(mock2017.metaTree.squ)
```

#Chi-Square test
##on class count for all sampling date
###2016
```{r}
mock2016.metaTree
mock2016.metaTree.df <-  metaTree2df(mock2016.metaTree)
chisq.test(apply(mock2016.metaTree.df, 2, sum))
```

###2017
```{r}
mock2017.metaTree.df <-  metaTree2df(mock2017.metaTree)
chisq.test(apply(mock2017.metaTree.df, 2, sum))
```

```{r}
mock2016.metaTree.rye <- mock2016.metaTree[which(mock2016.metaTree$material=="Rye"),]
mock2016.metaTree.squ <- mock2016.metaTree[which(mock2016.metaTree$material=="Squash"),]

mock2016.metaTree.rye.df <-metaTree2df(mock2016.metaTree.rye)
mock2016.metaTree.squ.df <-metaTree2df(mock2016.metaTree.squ)
mock2016.metaTree.rye.df <- mock2016.metaTree.rye.df[c("Pre", "Early", "Mid", "Late"),]
mock2016.metaTree.squ.df <- mock2016.metaTree.squ.df[c("Early", "Mid", "Late"),]
```


```{r}
mock2017.metaTree.rye <- mock2017.metaTree[which(mock2017.metaTree$material=="Rye"),]
mock2017.metaTree.squ <- mock2017.metaTree[which(mock2017.metaTree$material=="Squash"),]

mock2017.metaTree.rye.df <-metaTree2df(mock2017.metaTree.rye)
mock2017.metaTree.squ.df <-metaTree2df(mock2017.metaTree.squ)
mock2017.metaTree.rye.df <- mock2017.metaTree.rye.df[c("Pre", "Early", "Mid", "Late"),]
mock2017.metaTree.squ.df <- mock2017.metaTree.squ.df[c("Early", "Mid", "Late"),]
```

##Gammaproteobacteria count between sampling date
###Rye:2016&2017
```{r}
chisq.test(mock2017.metaTree.rye.df[,"Gammaproteobacteria"]+mock2016.metaTree.rye.df[,"Gammaproteobacteria"])
#chisq.test(mockAllYear.metaTree.rye.df.rel.final[,"Gammaproteobacteria"])
```

###Squash:2016&2017
```{r}
chisq.test(mock2017.metaTree.squ.df[,"Gammaproteobacteria"]+mock2016.metaTree.squ.df[,"Gammaproteobacteria"])
#chisq.test(mockAllYear.metaTree.squ.df.rel.final[,"Gammaproteobacteria"])
```

##Actinobacteria count between sampling date
###Rye:2016&2017
```{r}
chisq.test(mock2016.metaTree.rye.df[,"Actinobacteria"]+mock2017.metaTree.rye.df[,"Actinobacteria"])
```
###Squash: 2016&2017
```{r}
chisq.test(mock2016.metaTree.squ.df[,"Actinobacteria"]+mock2017.metaTree.squ.df[,"Actinobacteria"])
```

##Alphaproteobacteria count between sampling date
###Rye:2016&2017
```{r}
chisq.test(mock2016.metaTree.rye.df[,"Alphaproteobacteria"]+mock2017.metaTree.rye.df[,"Alphaproteobacteria"])
```

###Squash: 2016&2017
```{r}
chisq.test(mock2016.metaTree.squ.df[,"Alphaproteobacteria"]+mock2017.metaTree.squ.df[,"Alphaproteobacteria"])
```

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

#Supplemental Figure 5
```{r}
g0 <- ggplot(df1, aes(y=P.s.courge.Per, x=xdate))+  geom_boxplot(aes(fill=xdate, colour=col)) + scale_fill_identity() +
  #scale_fill_manual(values = cbbPalette,name="Day after pathogen inoculation", labels=c("1d", "2d", "3d", "4d", "7d")) + 
  theme_bw() + 
  #scale_color_manual(values="grey70") + 
  facet_grid(.~id) + ylim(c(0,1)) + ylab("% pathogen inhibition") +  theme(legend.position="bottom") 

g.leg <-g_legend(g0)
grid.draw(g.leg)

g1<-ggplot(df3, aes(y=P.s.courge.Per*100, x=xdate))+  geom_boxplot(aes(fill=date)) + 
  scale_fill_manual(values = cbbPalette,name="Day after\nPathogen\ninoculation", labels=c("1d", "2d", "3d", "4d", "7d")) + 
  theme_bw() + scale_color_manual(values="grey70") + facet_grid(.~id) + ylab("Pathogen inhibition (%)") + xlab("Day after\nPathogen\ninoculation") +#theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  guides(fill=FALSE) + ggtitle(expression(paste("A.\nPathogen = ",italic("P. syringae"), " ; Host = ", italic("C. pepo")))) + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

g2<-ggplot(df3, aes(y=X.h.laitue.Per*100, x=xdate))+  geom_boxplot(aes(fill=date)) + 
  scale_fill_manual(values = cbbPalette,name="Day after\nPathogen\ninoculation", labels=c("1d", "2d", "3d", "4d", "7d")) +
  theme_bw() + scale_color_manual(values="grey70") + facet_grid(.~id) + ylab("Pathogen inhibition (%)") +   xlab("Day after\nPathogen\ninoculation") +#theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  guides(fill=FALSE)  + ggtitle(expression(paste("C.\nPathogen = ",italic("X. hortorum"), " ; Host = ", italic("L. sativa")))) + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

g3<-ggplot(df3, aes(y=X.c.chou.Per*100, x=xdate))+  geom_boxplot(aes(fill=date)) + 
  scale_fill_manual(values = cbbPalette,name="Day after\nPathogen\ninoculation", labels=c("1d", "2d", "3d", "4d", "7d")) +
  theme_bw() + scale_color_manual(values="grey70") + facet_grid(.~id)  + ylab("Pathogen inhibition (%)") + xlab("Day after\nPathogen\ninoculation") + #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill=FALSE) + ggtitle(expression(paste("B.\nPathogen = ",italic("X. campestris"), " ; Host = ", italic("B. oleracea")))) + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

g4<-ggplot(df3, aes(y=P.s.haricot.Per*100, x=xdate))+  geom_boxplot(aes(fill=date)) + 
  scale_fill_manual(values = cbbPalette,name="Day after\nPathogen\ninoculation", labels=c("1d", "2d", "3d", "4d", "7d")) + 
  theme_bw() + scale_color_manual(values="grey70") + facet_grid(.~id) + ylab("Pathogen inhibition (%)") +  xlab("Day after\nPathogen\ninoculation") + #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  guides(fill=FALSE)  +  ggtitle(expression(paste("D.\nPathogen = ",italic("X. syringae"), " ; Host = ", italic("P. vulgaris")))) + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))  +  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

multiplot(g1, g2, g3, g4, cols = 2)
```
