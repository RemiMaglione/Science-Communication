---
title: "Script final"
output: html_notebook 
---

```{r, echo=F, eval=F}
library(ggplot2)
library(ggsignif)
library(lme4)
library(nlme)
library(emmeans)
library(lmerTest)
library(pscl)
library(MASS)
library(picante)
library(vegan)
library(gridExtra)
library(CoDaSeq)
library(grid)
library(phyloseq)
library(zCompositions)
library(DESeq2)
library(dplyr)
library(ggtree)
library(gridBase)
library(foreach)
library(doParallel)
library(tidyr)
library(multcomp)
```


############################################################################################################
#Cover cropping reduced P. syringae abundance on squash leaves and improved fruit health and marketability.#
############################################################################################################
## 2016
## INPUT DATA
```{r input 2016 Fig1, echo=F, eval=F}
metadata2016<- read.csv(file = "./metadata_2016.csv")

metadata.courge2016 <- metadata2016[which(metadata2016$type=="Courge"),]
metadata.courge2016$trt <- factor(metadata.courge2016$trt)
metadata.courge2016 <- within(metadata.courge2016, trt <- relevel(trt, ref = "Sol_nu"))
metadata.courge2016$date <- factor(metadata.courge2016$date)
metadata.courge2016$block <- factor(metadata.courge2016$block)
metadata.courge2016$pseudo_rawcount <- as.numeric(metadata.courge2016$pseudo_rawcount)
metadata.courge2016$vol_saline <- with(metadata.courge2016,as.numeric(as.character(vol_saline)))

metadata.courge2016$count_rif_zero <- ifelse(metadata.courge2016$count_rif==0.1, 
                                             yes=metadata.courge2016$count_rif_zero<-c(0),
                                             no=metadata.courge2016$count_rif_zero<- metadata.courge2016$count_rif)

metadata.courge2016$pseudo_rawcount_zero <- with(metadata.courge2016, round((count_rif_zero*as.numeric(vol_saline))/(poids*0.01*dil_rif)))

metadata.courge2016$log_pseudo_zero <- log10(metadata.courge2016$pseudo_rawcount_zero)

metadata.courge2016$log_pseudo_zero[metadata.courge2016$log_pseudo_zero==-Inf] <- 0

metadata.courge2016$log_pseudo_zero_int <- as.integer(round(metadata.courge2016$log_pseudo_zero))

metadata.courge2016.d1 <- metadata.courge2016[metadata.courge2016$date=="Date 1",]
metadata.courge2016.d2 <- metadata.courge2016[metadata.courge2016$date=="Date 2",]
metadata.courge2016.d3 <- metadata.courge2016[metadata.courge2016$date=="Date 3",]

date_name <- c(`Date 1` = "Early season",
               `Date 2` = "Mid season",
               `Date 3` = "Late season")
```

##Generalized linear Mixed model by date
```{r model lmer 2016}
metadata.courge2016.d1.lmer1 <- lmerTest::lmer(log_pseudo_zero_int ~ trt + (1|block), data = metadata.courge2016.d1)
metadata.courge2016.d2.lmer1 <- lmerTest::lmer(log_pseudo_zero_int ~ trt + (1|block), data = metadata.courge2016.d2)
metadata.courge2016.d3.lmer1 <- lmerTest::lmer(log_pseudo_zero_int ~ trt + (1|block), data = metadata.courge2016.d3)

summary(metadata.courge2016.d1.lmer1)
summary(metadata.courge2016.d2.lmer1)
summary(metadata.courge2016.d3.lmer1)
```
##### No effect of blocking variable -> use linear model instead
```{r model fig1 2016}
metadata.courge2016.d1.lm <- lm(log_pseudo_zero_int ~ trt, data = metadata.courge2016.d1)
metadata.courge2016.d2.lm <- lm(log_pseudo_zero_int ~ trt, data = metadata.courge2016.d2)
metadata.courge2016.d3.lm <- lm(log_pseudo_zero_int ~ trt, data = metadata.courge2016.d3)

metadata.courge2016.d1.lm.aov <- aov(metadata.courge2016.d1.lm)
metadata.courge2016.d2.lm.aov <- aov(metadata.courge2016.d2.lm)
metadata.courge2016.d3.lm.aov <- aov(metadata.courge2016.d3.lm)

metadata.courge2016.d1.lm.tk <- TukeyHSD(metadata.courge2016.d1.lm.aov)
metadata.courge2016.d2.lm.tk <- TukeyHSD(metadata.courge2016.d2.lm.aov)
metadata.courge2016.d3.lm.tk <- TukeyHSD(metadata.courge2016.d3.lm.aov)

metadata.courge2016.d1.lm.tk
metadata.courge2016.d2.lm.tk
metadata.courge2016.d3.lm.tk
```

|                    trt                    |signif group (date1)|signif group (date2)|signif group (date3)|
| ----------------------------------------- |:------------------:|:------------------:|:------------------:|
| `r unique(metadata.courge2016.d1$trt)[1]` |            c       |          a         |          a         |
| `r unique(metadata.courge2016.d1$trt)[2]` |            c       |          a         |          a         |
| `r unique(metadata.courge2016.d1$trt)[3]` |          ab        |          a         |          a         |
| `r unique(metadata.courge2016.d1$trt)[4]` |           bc       |          a         |          a         |


##Figure 1 (top panel : 2016)
```{r violin fig1 2016}
label2016.pos <- max(metadata.courge2016$log_pseudo_zero) + 5*max(metadata.courge2016$log_pseudo_zero)/100

trt2016 <- c("Sol_nu","Plastique", "Seigle_Bru", "Seigle")
pos2016 <- c(label2016.pos, label2016.pos, label2016.pos, label2016.pos)
date1_2016 <- c("Date 1", "Date 1", "Date 1", "Date 1")
date2_2016 <- c("Date 2", "Date 2", "Date 2", "Date 2")
date3_2016 <- c("Date 3", "Date 3", "Date 3", "Date 3")
signif.2016.d1 <- c("c", "c", "bc", "ab")
signif.2016.d2 <- c("a", "a", "a", "a")
signif.2016.d3 <- c("a", "a", "a", "a")

label2016.df <- data.frame(trt=c(trt2016, trt2016, trt2016), 
                       y=c(pos2016, pos2016, pos2016), 
                       label=c(signif.2016.d1, signif.2016.d2, signif.2016.d3),
                       date= c(date1_2016, date2_2016, date3_2016))

y_title <- expression(paste(italic("P.syringae"), " count (log10)"))

pseudo2016.ggv <- ggplot(data=metadata.courge2016, aes(x=trt, y=log_pseudo_zero, fill=trt)) +
  geom_violin() + geom_jitter(aes(alpha=I(0.3)))+ facet_grid(~date, labeller = as_labeller(date_name)) + 
  ylab(y_title)+ 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_fill_manual(values = c(colvec[4], colvec[1], colvec[2], colvec[3]), name=c("Legend"),labels=c("BS","PC","RCC","CT-RCC")) +
  geom_text(data=label2016.df, aes(x=trt, y=y, label=label))
  
pseudo2016.ggv

ggsave(filename = "pseudo2016_final.svg", path = "._article1/", plot = pseudo2016.ggv, width = 8.5, height = 3, units = "in")
```

## 2017
## INUPUT DATA
```{r input fig1 2017, echo=F, eval=F}
metadata2017 <- read.csv(file = "./metadata.courge2017", sep = "\t", header = T)

metadata2017 <- metadata2017[which(metadata2017$type=="courge"),]
metadata2017$trt <- factor(metadata2017$trt)
metadata2017 <- within(metadata2017, trt <- relevel(trt, ref = "Sol_nu"))
metadata2017$date <- factor(metadata2017$date)
metadata2017$block <- factor(metadata2017$block)
metadata2017$log_cfu_rif_zero <- as.numeric(metadata2017$log_cfu_rif_zero)

metadata2017.d1 <- metadata2017[metadata2017$date=="Date 1",]
metadata2017.d2 <- metadata2017[metadata2017$date=="Date 2",]
metadata2017.d3 <- metadata2017[metadata2017$date=="Date 3",]

date_name <- c(`Date 1` = "Early season",
               `Date 2` = "Mid season",
               `Date 3` = "Late season")
```

## Generalized linear Mixed model by date
```{r lmer model 2017}
metadata2017.d1.lmer1 <- lmerTest::lmer(log_cfu_rif_zero ~ trt + (1|block), data = metadata2017.d1)
metadata2017.d2.lmer1 <- lmerTest::lmer(log_cfu_rif_zero ~ trt + (1|block), data = metadata2017.d2)
metadata2017.d3.lmer1 <- lmerTest::lmer(log_cfu_rif_zero ~ trt + (1|block), data = metadata2017.d3)

summary(metadata2017.d1.lmer1)
summary(metadata2017.d2.lmer1)
summary(metadata2017.d3.lmer1)
```
```{r}
emmeans(metadata2017.d1.lmer1, list(pairwise ~ trt), adjust = "tukey")
emmeans(metadata2017.d2.lmer1, list(pairwise ~ trt), adjust = "tukey")
emmeans(metadata2017.d3.lmer1, list(pairwise ~ trt), adjust = "tukey")
```

|              trt                   |signif group (date1)|signif group (date2)|signif group (date3)|
| ---------------------------------- |:------------------:|:------------------:|:------------------:|
| `r unique(metadata2017.d1$trt)[2]` |            c       |             d      |           b        |
| `r unique(metadata2017.d1$trt)[4]` |           b        |            c       |          ab        |
| `r unique(metadata2017.d1$trt)[1]` |          a         |           bc       |          ab        |
| `r unique(metadata2017.d1$trt)[3]` |          a         |          ab        |          a         |

##Figure 1 (bottom panel : 2017)
```{r fig1 2017}
trt2017 <- c("Sol_nu","Plastique", "Seigle_bru", "Seigle")
pos2017 <- c(label.pos, label.pos, label.pos, label.pos)
date1_2017 <- c("Date 1", "Date 1", "Date 1", "Date 1")
date2_2017 <- c("Date 2", "Date 2", "Date 2", "Date 2")
date3_2017 <- c("Date 3", "Date 3", "Date 3", "Date 3")
signif.2017.d1 <- c("c", "b", "a", "a")
signif.2017.d2 <- c("d", "c", "bc", "ab")
signif.2017.d3 <- c("b", "ab", "ab", "a")

label2017.df <- data.frame(trt=c(trt2017, trt2017, trt2017), 
                       y=c(pos2017, pos2017, pos2017), 
                       label=c(signif.2017.d1, signif.2017.d2, signif.2017.d3),
                       date= c(date1_2017, date2_2017, date3_2017))

y_title <- expression(paste(italic("P.syringae"), " count (log10)"))

pseudo2017.ggv <- ggplot(data=metadata2017, aes(x=trt, y=log_cfu_rif_zero, fill=trt)) +
  geom_violin() + geom_jitter(aes(alpha=I(0.3)))+ facet_grid(~date, labeller = as_labeller(date_name)) + 
  ylab(y_title)+ 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_fill_manual(values = c(colvec[4], colvec[1], colvec[2], colvec[3]), name=c("Legend"),labels=c("BS","PC","RCC","CT-RCC")) +
  geom_text(data=label2017.df, aes(x=trt, y=y, label=label))

pseudo2017.ggv
  
ggsave(filename = "pseudo2017_final.svg", path = "./", plot = pseudo2017.ggv, width = 8.5, height = 3, units = "in")
```


##Table 1: Proportion of squash fruit (mean +- standard deviation) with no P. syringae symptoms and marketable fruits with no damage for the two-growing seasons 2016 and 2017. 
### 2016 section
#### INPUT DATA
```{r}
courgeAgro2016.noNA <- read.table(file = './courgeAgro2016.noNA.csv', sep=",", header = T)

courgeAgro2016.noNA <- courgeAgro2016.noNA %>% 
  mutate(HV.Pseudo.cat=replace(HV.Pseudo, HV.Pseudo<10, 0)) %>%
  mutate(HV.Pseudo.cat=replace(HV.Pseudo.cat, HV.Pseudo.cat>=10, 1))

courgeAgro2016.noNA <- courgeAgro2016.noNA %>% 
  mutate(HV.Pénétrant.cat=replace(HV.Pénétrant, HV.Pénétrant<5, 0)) %>%
  mutate(HV.Pénétrant.cat=replace(HV.Pénétrant.cat, HV.Pénétrant.cat>=5, 1))

courgeAgro2016.noNA <- courgeAgro2016.noNA %>% 
  mutate(Pseudo.Cicatrice.cat=replace(Pseudo.Cicatrice, Pseudo.Cicatrice<10, 0)) %>%
  mutate(Pseudo.Cicatrice.cat=replace(Pseudo.Cicatrice.cat, Pseudo.Cicatrice.cat>=10, 1))

courgeAgro2016.noNA <- courgeAgro2016.noNA %>% 
  mutate(Black.Rot.cat=replace(Black.Rot, Black.Rot<5, 0)) %>%
  mutate(Black.Rot.cat=replace(Black.Rot.cat, Black.Rot.cat>=5, 1))

courgeAgro2016.noNA$marketable <- 0
courgeAgro2016.noNA <- courgeAgro2016.noNA %>%
  mutate(marketable=replace(marketable, rowSums(.[12:15])>0, 1))

courgeAgro2016.noNA$pseudo <- 0
courgeAgro2016.noNA <- courgeAgro2016.noNA %>%
  mutate(pseudo=replace(pseudo, rowSums(.[c(7:11)])>4, 1))

courgeSaine2016.avg <- courgeAgro2016.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(HV.Pseudo.cat.avg = mean(HV.Pseudo.cat))

courgeSaine2016.avg <- courgeAgro2016.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(HV.Pénétrant.cat.avg = mean(HV.Pénétrant.cat))

courgeSaine2016.avg <- courgeAgro2016.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(Pseudo.Cicatrice.cat.avg = mean(Pseudo.Cicatrice.cat))

courgeSaine2016.avg <- courgeAgro2016.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(Black.Rot.cat.avg = mean(Black.Rot.cat))

courgeSaine2016.avg <- courgeAgro2016.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(marketable.avg = mean(marketable))

courgeAgro2016.noNA.trt <- courgeAgro2016.noNA %>%
  group_by(Rep1, Traitement) %>%
  count(Traitement)
```

### Proportion of squash fruit without P. syringae symptoms (%) : mean and sd
```{r}
courgeAgro2016.noNA.trt <- courgeAgro2016.noNA %>%
  group_by(Rep1, Traitement) %>%
  count(Traitement)


courgeAgro2016.noNA.pseudo <- courgeAgro2016.noNA %>%
  group_by(Rep1, Traitement) %>%
  filter(pseudo==0) %>%
  count(Count=pseudo)

courgeAgro2016.noNA.pseudo$prop <- courgeAgro2016.noNA.pseudo[,"n"]/courgeAgro2016.noNA.trt[,"n"]
courgeAgro2016.noNA.pseudo.sub <- data.frame(Rep1=courgeAgro2016.noNA.pseudo$Rep1,
                                              Traitement=courgeAgro2016.noNA.pseudo$Traitement,
                                              prop=as.numeric(unlist(round(courgeAgro2016.noNA.pseudo[,"n"]/courgeAgro2016.noNA.trt[,"n"],digits=2))))

courgeAgro2016.noNA.pseudo.sub %>%
  group_by(Traitement) %>%
  summarise(avg=round(mean(prop)*100,digits=3), sd=round(sd(prop)*100,digits=3))
```

### squash fruit without P. syringae symptoms : generalized mixed model
```{r}
pseudo2016.glmer <- glmer(pseudo ~ Traitement + (1|Rep1), data = courgeAgro2016.noNA, family = "binomial")
summary(glht(pseudo2016.glmer, mcp(Traitement="Tukey")))
```
|   trt   |signif group(pseudo)|
| ------- |:------------------:|
| Plastic |         ab         | 
| Rye     |         a          |
| Rye+Gly |         ab         |
| Soil    |          b         |

### Proportion of marketable squash fruit with no damages (%) : mean and sd
```{r}
courgeAgro2016.noNA.marketable <- courgeAgro2016.noNA %>%
  group_by(Rep1, Traitement) %>%
  filter(marketable==0) %>%
  count(Count=marketable)

courgeAgro2016.noNA.marketable$prop <- courgeAgro2016.noNA.marketable[,"n"]/courgeAgro2016.noNA.trt[,"n"]
courgeAgro2016.noNA.marketable.sub <- data.frame(Rep1=courgeAgro2016.noNA.marketable$Rep1,
                                              Traitement=courgeAgro2016.noNA.marketable$Traitement,
                                              prop=as.numeric(unlist(round(courgeAgro2016.noNA.marketable[,"n"]/courgeAgro2016.noNA.trt[,"n"],digits=2))))

courgeAgro2016.noNA.marketable.sub %>%
  group_by(Traitement) %>%
  summarise(avg=round(mean(prop)*100,digits=3), sd=round(sd(prop)*100,digits=3))
```

### marketable squash fruit : generalized linear model 
```{r}
marketable.glm <- glm(marketable ~ relevel(Traitement, ref = "Soil"), data = courgeAgro2016.noNA, family = "binomial")
TukeyHSD(aov(marketable.glm))
```

|   trt   |signif group(market)|
| ------- |:------------------:|
| Plastic |         a          | 
| Rye     |          b         |
| Rye+Gly |         ab         |
| Soil    |         a          |

### 2017 section
#### INPUT DATA
```{r}
courgeAgro2017.noNA <- read.table(file = './courgeAgro2017.noNA.csv', sep=",", header = T)

courgeAgro2017.noNA <- courgeAgro2017.noNA %>% 
  mutate(HV.Pseudo.cat=replace(HV.Pseudo, HV.Pseudo<10, 0)) %>%
  mutate(HV.Pseudo.cat=replace(HV.Pseudo.cat, HV.Pseudo.cat>=10, 1))

courgeAgro2017.noNA <- courgeAgro2017.noNA %>% 
  mutate(HV.Pénétrant.cat=replace(HV.Pénétrant, HV.Pénétrant<5, 0)) %>%
  mutate(HV.Pénétrant.cat=replace(HV.Pénétrant.cat, HV.Pénétrant.cat>=5, 1))

courgeAgro2017.noNA <- courgeAgro2017.noNA %>% 
  mutate(Pseudo.Cicatrice.cat=replace(Pseudo.Cicatrice, Pseudo.Cicatrice<10, 0)) %>%
  mutate(Pseudo.Cicatrice.cat=replace(Pseudo.Cicatrice.cat, Pseudo.Cicatrice.cat>=10, 1))

courgeAgro2017.noNA <- courgeAgro2017.noNA %>% 
  mutate(Black.Rot.cat=replace(Black.Rot, Black.Rot<5, 0)) %>%
  mutate(Black.Rot.cat=replace(Black.Rot.cat, Black.Rot.cat>=5, 1))

courgeAgro2017.noNA$marketable <- 0
courgeAgro2017.noNA <- courgeAgro2017.noNA %>%
  mutate(marketable=replace(marketable, rowSums(.[12:15])>0, 1))

courgeAgro2017.noNA$pseudo <- 0
courgeAgro2017.noNA <- courgeAgro2017.noNA %>%
  mutate(pseudo=replace(pseudo, rowSums(.[c(7:11)])>4, 1))

courgeSaine2017.avg <- courgeAgro2017.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(HV.Pseudo.cat.avg = mean(HV.Pseudo.cat))

courgeSaine2017.avg <- courgeAgro2017.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(HV.Pénétrant.cat.avg = mean(HV.Pénétrant.cat))

courgeSaine2017.avg <- courgeAgro2017.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(Pseudo.Cicatrice.cat.avg = mean(Pseudo.Cicatrice.cat))

courgeSaine2017.avg <- courgeAgro2017.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(Black.Rot.cat.avg = mean(Black.Rot.cat))

courgeSaine2017.avg <- courgeAgro2017.noNA %>%
  group_by(X.Par, Traitement) %>%
  summarise(marketable.avg = mean(marketable))
```

### Proportion of squash fruit without P. syringae symptoms (%) : meand and sd
```{r 2017 symptoms}
courgeAgro2017.noNA.trt <- courgeAgro2017.noNA %>%
  group_by(Rep1, Traitement) %>%
  count(Traitement)


courgeAgro2017.noNA.pseudo <- courgeAgro2017.noNA %>%
  group_by(Rep1, Traitement) %>%
  filter(pseudo==0) %>%
  count(Count=pseudo)

courgeAgro2017.noNA.pseudo$prop <- courgeAgro2017.noNA.pseudo[,"n"]/courgeAgro2017.noNA.trt[,"n"]
courgeAgro2017.noNA.pseudo.sub <- data.frame(Rep1=courgeAgro2017.noNA.pseudo$Rep1,
                                              Traitement=courgeAgro2017.noNA.pseudo$Traitement,
                                              prop=as.numeric(unlist(round(courgeAgro2017.noNA.pseudo[,"n"]/courgeAgro2017.noNA.trt[,"n"],digits=2))))

courgeAgro2017.noNA.pseudo.sub %>%
  group_by(Traitement) %>%
  summarise(avg=round(mean(prop)*100,digits=3), sd=round(sd(prop)*100,digits=3))
```

### squash fruit without P. syringae symptoms : generalized mixed model 
```{r}
pseudo2017.glmer <- glmer(pseudo ~ Traitement + (1|Rep1), data = courgeAgro2017.noNA, family = "binomial")
summary(glht(pseudo2017.glmer, mcp(Traitement="Tukey")))
```
|   trt   |signif group (date1)|
| ------- |:------------------:|
| Plastic |         ab         | 
| Rye     |         ab         |
| Rye+Gly |         a          |
| Soil    |          b         |


### Proportion of marketable squash fruit with no damages (%)  : meand and sd
```{r 2017 market}
courgeAgro2017.noNA.trt <- courgeAgro2017.noNA %>%
  group_by(Rep1, Traitement) %>%
  count(Traitement)


courgeAgro2017.noNA.marketable <- courgeAgro2017.noNA %>%
  group_by(Rep1, Traitement) %>%
  filter(marketable==0) %>%
  count(Count=marketable)

courgeAgro2017.noNA.marketable$prop <- courgeAgro2017.noNA.marketable[,"n"]/courgeAgro2017.noNA.trt[,"n"]
courgeAgro2017.noNA.marketable.sub <- data.frame(Rep1=courgeAgro2017.noNA.marketable$Rep1,
                                              Traitement=courgeAgro2017.noNA.marketable$Traitement,
                                              prop=as.numeric(unlist(round(courgeAgro2017.noNA.marketable[,"n"]/courgeAgro2017.noNA.trt[,"n"],digits=2))))

courgeAgro2017.noNA.marketable.sub %>%
  group_by(Traitement) %>%
  summarise(avg=round(mean(prop)*100,digits=3), sd=round(sd(prop)*100,digits=3))
```

### Proportion of marketable squash fruit with no damages (%) : generalized mixed model 
```{r}
marketable2017.glmer <- glmer(marketable ~ Traitement + (1|Rep1), data = courgeAgro2017.noNA, family = "binomial")
summary(glht(marketable2017.glmer, mcp(Traitement="Tukey")))
```
|   trt   |signif group (date1)|
| ------- |:------------------:|
| Plastic |         b          | 
| Rye     |        a           |
| Rye+Gly |        ab          |
| Soil    |         b          |


##Supplemental Table 1. Number of P. syringae colony forming units (CFUs) recovered from squash leaves grown in different cover cropping practices#
### 2016 section
```{r}
#date 1
###RCC
d1.RCC <- round(10^with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)
d1.RCC.sd <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle",], sd(log_pseudo_zero)),digits = 2)
d1.RCC.mean <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)

###CT-RCC
d1.CT_RCC <- round(10^with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
d1.CT_RCC.sd <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle_Bru",], sd(log_pseudo_zero)),digits = 2)
d1.CT_RCC.mean <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
###PC
d1.PC <- round(10^with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
d1.PC.sd <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Plastique",], sd(log_pseudo_zero)),digits = 2)
d1.PC.mean <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
###NS
d1.NS <- round(10^with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)
d1.NS.sd <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Sol_nu",], sd(log_pseudo_zero)),digits = 2)
d1.NS.mean <- round(with(metadata.courge2016.d1[metadata.courge2016.d1$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)

#date 2
###RCC
d2.RCC <- round(10^with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)
d2.RCC.sd <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle",], sd(log_pseudo_zero)),digits = 2)
d2.RCC.mean <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)
###CT-RCC
d2.CT_RCC <- round(10^with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
d2.CT_RCC.sd <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle_Bru",], sd(log_pseudo_zero)),digits = 2)
d2.CT_RCC.mean <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
###PC
d2.PC <- round(10^with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
d2.PC.sd <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Plastique",], sd(log_pseudo_zero)),digits = 2)
d2.PC.mean <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
###NS
d2.NS <- round(10^with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)
d2.NS.sd <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Sol_nu",], sd(log_pseudo_zero)),digits = 2)
d2.NS.mean <- round(with(metadata.courge2016.d2[metadata.courge2016.d2$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)

#date 3 
###RCC
d3.RCC <- round(10^with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)
d3.RCC.sd <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle",], sd(log_pseudo_zero)),digits = 2)
d3.RCC.mean <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle",], mean(log_pseudo_zero)),digits = 2)
###CT-RCC
d3.CT_RCC <- round(10^with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
d3.CT_RCC.sd <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle_Bru",], sd(log_pseudo_zero)),digits = 2)
d3.CT_RCC.mean <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Seigle_Bru",], mean(log_pseudo_zero)),digits = 2)
###PC
d3.PC <- round(10^with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
d3.PC.sd <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Plastique",], sd(log_pseudo_zero)),digits = 2)
d3.PC.mean <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Plastique",], mean(log_pseudo_zero)),digits = 2)
###NS
d3.NS <- round(10^with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)
d3.NS.sd <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Sol_nu",], sd(log_pseudo_zero)),digits = 2)
d3.NS.mean <- round(with(metadata.courge2016.d3[metadata.courge2016.d3$trt=="Sol_nu",], mean(log_pseudo_zero)),digits = 2)
```

| date  |      mean RCC (log)           |         mean CT-RCC (log)           |         mean PC (log)       |        mean NS (log)        |
| ----- |:-----------------------------:|:-----------------------------------:|:---------------------------:|:---------------------------:|
| Date 1|`r d1.RCC.mean` ± `r d1.RCC.sd`|`r d1.CT_RCC.mean` ± `r d1.CT_RCC.sd`|`r d1.PC.mean` ± `r d1.PC.sd`|`r d1.NS.mean` ± `r d1.NS.sd`|
| Date 2|`r d2.RCC.mean` ± `r d2.RCC.sd`|`r d2.CT_RCC.mean` ± `r d2.CT_RCC.sd`|`r d2.PC.mean` ± `r d2.PC.sd`|`r d2.NS.mean` ± `r d2.NS.sd`|
| Date 3|`r d3.RCC.mean` ± `r d3.RCC.sd`|`r d3.CT_RCC.mean` ± `r d3.CT_RCC.sd`|`r d3.PC.mean` ± `r d3.PC.sd`|`r d3.NS.mean` ± `r d3.NS.sd`|

###2017 Section
```{r}
#date 1
###RCC
d1.RCC2017 <- round(10^with(metadata2017.d1[metadata2017.d1$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)
d1.RCC2017.sd <- round(with(metadata2017.d1[metadata2017.d1$trt=="Seigle",], sd(log_cfu_rif_zero)),digits = 2)
d1.RCC2017.mean <- round(with(metadata2017.d1[metadata2017.d1$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)

###CT-RCC
d1.CT_RCC2017 <- round(10^with(metadata2017.d1[metadata2017.d1$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
d1.CT_RCC2017.sd <- round(with(metadata2017.d1[metadata2017.d1$trt=="Seigle_bru",], sd(log_cfu_rif_zero)),digits = 2)
d1.CT_RCC2017.mean <- round(with(metadata2017.d1[metadata2017.d1$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
###PC
d1.PC2017 <- round(10^with(metadata2017.d1[metadata2017.d1$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
d1.PC2017.sd <- round(with(metadata2017.d1[metadata2017.d1$trt=="Plastique",], sd(log_cfu_rif_zero)),digits = 2)
d1.PC2017.mean <- round(with(metadata2017.d1[metadata2017.d1$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
###NS
d1.NS2017 <- round(10^with(metadata2017.d1[metadata2017.d1$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)
d1.NS2017.sd <- round(with(metadata2017.d1[metadata2017.d1$trt=="Sol_nu",], sd(log_cfu_rif_zero)),digits = 2)
d1.NS2017.mean <- round(with(metadata2017.d1[metadata2017.d1$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)

#date 2
###RCC
d2.RCC2017 <- round(10^with(metadata2017.d2[metadata2017.d2$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)
d2.RCC2017.sd <- round(with(metadata2017.d2[metadata2017.d2$trt=="Seigle",], sd(log_cfu_rif_zero)),digits = 2)
d2.RCC2017.mean <- round(with(metadata2017.d2[metadata2017.d2$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)
###CT-RCC
d2.CT_RCC2017 <- round(10^with(metadata2017.d2[metadata2017.d2$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
d2.CT_RCC2017.sd <- round(with(metadata2017.d2[metadata2017.d2$trt=="Seigle_bru",], sd(log_cfu_rif_zero)),digits = 2)
d2.CT_RCC2017.mean <- round(with(metadata2017.d2[metadata2017.d2$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
###PC
d2.PC2017 <- round(10^with(metadata2017.d2[metadata2017.d2$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
d2.PC2017.sd <- round(with(metadata2017.d2[metadata2017.d2$trt=="Plastique",], sd(log_cfu_rif_zero)),digits = 2)
d2.PC2017.mean <- round(with(metadata2017.d2[metadata2017.d2$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
###NS
d2.NS2017 <- round(10^with(metadata2017.d2[metadata2017.d2$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)
d2.NS2017.sd <- round(with(metadata2017.d2[metadata2017.d2$trt=="Sol_nu",], sd(log_cfu_rif_zero)),digits = 2)
d2.NS2017.mean <- round(with(metadata2017.d2[metadata2017.d2$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)

#date 3 
###RCC
d3.RCC2017 <- round(10^with(metadata2017.d3[metadata2017.d3$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)
d3.RCC2017.sd <- round(with(metadata2017.d3[metadata2017.d3$trt=="Seigle",], sd(log_cfu_rif_zero)),digits = 2)
d3.RCC2017.mean <- round(with(metadata2017.d3[metadata2017.d3$trt=="Seigle",], mean(log_cfu_rif_zero)),digits = 2)
###CT-RCC
d3.CT_RCC2017 <- round(10^with(metadata2017.d3[metadata2017.d3$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
d3.CT_RCC2017.sd <- round(with(metadata2017.d3[metadata2017.d3$trt=="Seigle_bru",], sd(log_cfu_rif_zero)),digits = 2)
d3.CT_RCC2017.mean <- round(with(metadata2017.d3[metadata2017.d3$trt=="Seigle_bru",], mean(log_cfu_rif_zero)),digits = 2)
###PC
d3.PC2017 <- round(10^with(metadata2017.d3[metadata2017.d3$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
d3.PC2017.sd <- round(with(metadata2017.d3[metadata2017.d3$trt=="Plastique",], sd(log_cfu_rif_zero)),digits = 2)
d3.PC2017.mean <- round(with(metadata2017.d3[metadata2017.d3$trt=="Plastique",], mean(log_cfu_rif_zero)),digits = 2)
###NS
d3.NS2017 <- round(10^with(metadata2017.d3[metadata2017.d3$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)
d3.NS2017.sd <- round(with(metadata2017.d3[metadata2017.d3$trt=="Sol_nu",], sd(log_cfu_rif_zero)),digits = 2)
d3.NS2017.mean <- round(with(metadata2017.d3[metadata2017.d3$trt=="Sol_nu",], mean(log_cfu_rif_zero)),digits = 2)
```

| date  |      mean RCC (log)           |         mean CT-RCC (log)           |         mean PC (log)       |        mean NS (log)        |
| ----- |:-----------------------------:|:-----------------------------------:|:---------------------------:|:---------------------------:|
| Date 1|`r d1.RCC2017.mean` ± `r d1.RCC2017.sd`|`r d1.CT_RCC2017.mean` ± `r d1.CT_RCC2017.sd`|`r d1.PC2017.mean` ± `r d1.PC2017.sd`|`r d1.NS2017.mean` ± `r d1.NS2017.sd`|
| Date 2|`r d2.RCC2017.mean` ± `r d2.RCC2017.sd`|`r d2.CT_RCC2017.mean` ± `r d2.CT_RCC2017.sd`|`r d2.PC2017.mean` ± `r d2.PC2017.sd`|`r d2.NS2017.mean` ± `r d2.NS2017.sd`|
| Date 3|`r d3.RCC2017.mean` ± `r d3.RCC2017.sd`|`r d3.CT_RCC2017.mean` ± `r d3.CT_RCC2017.sd`|`r d3.PC2017.mean` ± `r d3.PC2017.sd`|`r d3.NS2017.mean` ± `r d3.NS2017.sd`|


##Supplemental Table 2
### 2016 section
| date          |   CT-RCC/RCC       |     PC/RCC     |   NS/RCC       |    PC/CT-RCC      |     NS/CT-RCC     |   NS/PC       |
| ------------- |:------------------:|:--------------:|:--------------:|:-----------------:|:-----------------:|:-------------:|
| Date 1        |`r d1.CT_RCC/d1.RCC`|`r d1.PC/d1.RCC`|`r d1.NS/d1.RCC`|`r d1.PC/d1.CT_RCC`|`r d1.NS/d1.CT_RCC`|`r d1.NS/d1.PC`|
| Date 2        |`r d2.CT_RCC/d2.RCC`|`r d2.PC/d2.RCC`|`r d2.NS/d2.RCC`|`r d2.PC/d2.CT_RCC`|`r d2.NS/d2.CT_RCC`|`r d2.NS/d2.PC`|
| Date 3        |`r d3.CT_RCC/d3.RCC`|`r d3.PC/d3.RCC`|`r d3.NS/d3.RCC`|`r d3.PC/d3.CT_RCC`|`r d3.NS/d3.CT_RCC`|`r d3.NS/d3.PC`|

### 2017 section
| date          |   CT-RCC/RCC       |     PC/RCC     |   NS/RCC       |    PC/CT-RCC      |     NS/CT-RCC     |   NS/PC       |
| ------------- |:------------------:|:--------------:|:--------------:|:-----------------:|:-----------------:|:-------------:|
| Date 1        |`r d1.CT_RCC2017/d1.RCC2017`|`r d1.PC2017/d1.RCC2017`|`r d1.NS2017/d1.RCC2017`|`r d1.PC2017/d1.CT_RCC2017`|`r d1.NS2017/d1.CT_RCC2017`|`r d1.NS2017/d1.PC2017`|
| Date 2        |`r d2.CT_RCC2017/d2.RCC2017`|`r d2.PC2017/d2.RCC2017`|`r d2.NS2017/d2.RCC2017`|`r d2.PC2017/d2.CT_RCC2017`|`r d2.NS2017/d2.CT_RCC2017`|`r d2.NS2017/d2.PC2017`|
| Date 3        |`r d3.CT_RCC2017/d3.RCC2017`|`r d3.PC2017/d3.RCC2017`|`r d3.NS2017/d3.RCC2017`|`r d3.PC2017/d3.CT_RCC2017`|`r d3.NS2017/d3.CT_RCC2017`|`r d3.NS2017/d3.PC2017`|


####################################################################################
# Phyllosphere microbial communities differed between sampling dates and treatments#
####################################################################################

##INPUT & PROCESSING COMMUNITY DATA
###2016
```{r comm taxo metadata 2016}
asvs2016.comm <- readRDS(file="./m1_cut/seqtabCollapse.220_200_2016_allRun_pPool_minOver_294_nochim.rds")
asvs2016.taxo <- readRDS(file="./taxa.2016.sp_220_200_allRun.rds")
metadata2016 <- read.table(file = "./metadata_2016.csv", sep = ",", header = T)
```

```{r matching metadata to comm and taxo matrix 2016, message = FALSE}
asvs2016.comm <- asvs2016.comm[apply(asvs2016.comm,1,sum)>0,]
both <- intersect(rownames(metadata2016), rownames(asvs2016.comm))
asvs.comm.courge2016 <- asvs2016.comm[both,]
asvs.comm.courge2016 <- asvs.comm.courge2016[,apply(asvs.comm.courge2016,2,sum)>1] # Only keep column with more than 1 seq
metadata.courge2016 <- metadata2016[both,]
asvs.taxo.courge2016 <- asvs2016.taxo[colnames(asvs.comm.courge2016),]
```

###2017
```{r comm taxo metadata 2017}
asvs.comm <- readRDS(file="./seqtabCollapse230_210_courge2017_allRun_pPool_minOver_250_nochim.rds")
asvs.taxo <- readRDS(file="./taxa.2017.sp230_210_courge2017_allRun.rds")
metadata <- read.csv(file = "./metadata.courge2017", comment.char = "", sep = "\t", check.names = FALSE, quote="\"",
                     na.strings=c("NA","NaN", " "))
rownames(metadata) <- metadata[,1]
```

```{r matching metadata to comm and taxo matrix 2017, message = FALSE}
both <- intersect(rownames(metadata), rownames(asvs.comm))
asvs.comm.2017 <- asvs.comm[both,]
asvs.comm.2017 <- asvs.comm.2017[,apply(asvs.comm.2017,2,sum)>1] # Only keep column with more than 1 seq
metadata.2017 <- metadata[both,]
asvs.taxo.2017 <- asvs.taxo[colnames(asvs.comm.2017),]
```

##(line 175) Supplemental figure 3
###2016
```{r clr and PCA 2016 for control}
asvs.comm.courge2016.f<-t(codaSeq.filter(asvs.comm.courge2016, min.reads=min(apply(asvs.comm.courge2016, 1, sum))-1, min.prop=0.00001, min.occurrence=0.005, samples.by.row=T))
asvs.comm.courge2016.f.n0 <- as.matrix(cmultRepl(asvs.comm.courge2016.f, method="CZM", label=0))

asvs.comm.courge2016.f.n0.clr <- codaSeq.clr(asvs.comm.courge2016.f.n0)

myGreyCol <- c("grey30", "grey50", "grey70", "grey10")
myGreyColAlpha <- makeTransparent(myGreyCol, alpha = 0.3)

par(mfrow=c(1,2))
#plot1
asvs.comm.courge2016.f.n0.clr.meta <- metadata.courge2016[rownames(asvs.comm.courge2016.f.n0.clr),]
pcx <- prcomp(asvs.comm.courge2016.f.n0.clr)
# plot a PCA
# calculate percent variance explained for the axis labels
pc1 <- round(pcx$sdev[1]^2/sum(pcx$sdev^2),2)
pc2 <- round(pcx$sdev[2]^2/sum(pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")

ordiplot(pcx, display="sites", type = 'points', main = "2016", xlab=xlab, ylab=ylab)

with(asvs.comm.courge2016.f.n0.clr.meta, ordiellipse(pcx, 
                                                   group = date, kind = "sd", 
                                                   label = TRUE,
                                                   col = myGreyCol))
with(asvs.comm.courge2016.f.n0.clr.meta, ordispider(pcx, 
                                                   group = date,
                                                   label = TRUE,
                                                   col = myGreyColAlpha))
```

###2017
```{r clr and PCA 2017 for control}
asvs.comm.2017.f<-t(codaSeq.filter(asvs.comm.2017, min.reads=min(apply(asvs.comm.2017, 1, sum))-1, min.prop=0.00001, min.occurrence=0.005, samples.by.row=T))
asvs.comm.2017.f.n0 <- as.matrix(cmultRepl(asvs.comm.2017.f, method="CZM", label=0))

asvs.comm.2017.f.n0.clr <- codaSeq.clr(asvs.comm.2017.f.n0)

myGreyCol <- c("grey30", "grey50", "grey70", "grey10")
myGreyColAlpha <- makeTransparent(myGreyCol, alpha = 0.3)

par(mfrow=c(1,2))
#plot1
asvs.comm.2017.f.n0.clr.meta <- metadata.2017[rownames(asvs.comm.2017.f.n0.clr),]
pcx <- prcomp(asvs.comm.2017.f.n0.clr)
# plot a PCA
# calculate percent variance explained for the axis labels
pc1 <- round(pcx$sdev[1]^2/sum(pcx$sdev^2),2)
pc2 <- round(pcx$sdev[2]^2/sum(pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")

ordiplot(pcx, display="sites", type = 'points', xlab=xlab, ylab=ylab)

with(asvs.comm.2017.f.n0.clr.meta, ordiellipse(pcx, 
                                                   group = date, kind = "sd", 
                                                   label = FALSE,
                                                   col = myGreyCol))
with(asvs.comm.2017.f.n0.clr.meta, ordispider(pcx, 
                                                   group = date,
                                                   label = TRUE,
                                                   col = myGreyColAlpha))
```

## Data processing Suite
###2016
#### Metada Processing
```{r metadata processing 2016}
metadata.courge2016 <- metadata.courge2016[metadata.courge2016$type=="Courge",]
metadata.courge2016 <- metadata.courge2016[which((!metadata.courge2016$rep_seq=="2")|is.na(!metadata.courge2016$rep_seq=="2")),]
metadata.courge2016$type <- factor(metadata.courge2016$type)
rownames(metadata.courge2016) <- metadata.courge2016[,1]
both <- intersect(rownames(metadata.courge2016), rownames(asvs.comm.courge2016))
asvs.comm.courge2016 <- asvs.comm.courge2016[both,]
asvs.comm.courge2016 <- asvs.comm.courge2016[,apply(asvs.comm.courge2016,2,sum)>1] # Only keep column with more than 1 seq
metadata.courge2016 <- metadata.courge2016[both,]
asvs.taxo.courge2016 <- asvs.taxo.courge2016[colnames(asvs.comm.courge2016),]
```

###2017
#### Metada Processing
```{r matching metadata to comm and taxo matrix, message = FALSE}
metadata.courge2017 <- metadata.2017[which(metadata.2017$type=="courge"),]
metadata.courge2017$type <- factor(metadata.courge2017$type)
both <- intersect(rownames(metadata.courge2017), rownames(asvs.comm.2017))
asvs.comm.courge2017 <- asvs.comm.2017[both,]
asvs.comm.courge2017 <- asvs.comm.courge2017[,apply(asvs.comm.courge2017,2,sum)>1] # Only keep column with more than 1 seq
metadata.courge2017 <- metadata.2017[both,]
asvs.taxo.courge2017 <- asvs.taxo[colnames(asvs.comm.courge2017),]
```

##(line 204) Supplemental Figure 4. Rarefaction curves (ASVs versus number of sequences/sample) for each squash sample#
###2016 : top panel
```{r rarecurve, message = FALSE, eval=FALSE}
rarecurve(asvs.comm.courge2016, step = 100)
```

###2017 : bottom panel
```{r}
rarecurve(asvs.comm.courge2017, step = 100)
```

##Rarefy data - keep only samples with 5k reads
####2016
```{r rarefy 2016, message = FALSE}
asvs.comm.courge2016.5k <- asvs.comm.courge2016[which(apply(asvs.comm.courge2016,1,sum)>=5000),]
asvs.comm.courge2016.5k <- asvs.comm.courge2016.5k[,apply(asvs.comm.courge2016.5k,2,sum)>1]
asvs.taxo.courge2016.5k <- asvs.taxo.courge2016[colnames(asvs.comm.courge2016.5k),]
metadata.courge2016.5k <- metadata.courge2016[rownames(asvs.comm.courge2016.5k),]
dim(asvs.comm.courge2016.5k) #matrix dimention before rarefaction
asvs.comm.courge2016.rare5k <- rrarefy(asvs.comm.courge2016.5k,sample=5000) 
asvs.comm.courge2016.rare5k <- asvs.comm.courge2016.rare5k[,apply(asvs.comm.courge2016.rare5k,2,sum)>1]
metadata.courge2016.rare5k <- metadata.courge2016.5k[rownames(asvs.comm.courge2016.rare5k),]
asvs.taxo.courge2016.rare5k <- asvs.taxo.courge2016.5k[colnames(asvs.comm.courge2016.rare5k),]
dim(asvs.comm.courge2016.rare5k) #matrix dimention after rarefaction
remove(asvs.comm.courge2016.5k)
remove(asvs.taxo.courge2016.5k)
remove(metadata.courge2016.5k)

# Fix factor class for metadata after rarefaction
metadata.courge2016.rare5k$date <- factor(metadata.courge2016.rare5k$date)
metadata.courge2016.rare5k$trt <- factor(metadata.courge2016.rare5k$trt)
metadata.courge2016.rare5k$type <- factor(metadata.courge2016.rare5k$type)
metadata.courge2016.rare5k$block <- factor(metadata.courge2016.rare5k$block)
```

###2017
```{r rarefy 2017, message = FALSE}
asvs.comm.courge2017.5k <- asvs.comm.courge2017[which(apply(asvs.comm.courge2017,1,sum)>=5000),]
asvs.comm.courge2017.5k <- asvs.comm.courge2017.5k[,apply(asvs.comm.courge2017.5k,2,sum)>1]
asvs.taxo.courge2017.5k <- asvs.taxo.courge2017[colnames(asvs.comm.courge2017.5k),]
metadata.courge2017.5k <- metadata.courge2017[rownames(asvs.comm.courge2017.5k),]
dim(asvs.comm.courge2017.5k) #matrix dimention before rarefaction

asvs.comm.courge2017.rare5k <- rrarefy(asvs.comm.courge2017.5k,sample=5000)
asvs.comm.courge2017.rare5k <- asvs.comm.courge2017.rare5k[,apply(asvs.comm.courge2017.rare5k,2,sum)>1]
metadata.courge2017.rare5k <- metadata.courge2017.5k[rownames(asvs.comm.courge2017.rare5k),]
asvs.taxo.courge2017.rare5k <- asvs.taxo.courge2017.5k[colnames(asvs.comm.courge2017.rare5k),]

dim(asvs.comm.courge2017.rare5k) #matrix dimention after rarefaction
remove(asvs.comm.courge2017.5k)
remove(asvs.taxo.courge2017.5k)
remove(metadata.courge2017.5k)

# Fix factor class for metadata after rarefaction
metadata.courge2017.rare5k$date <- factor(metadata.courge2017.rare5k$date)
metadata.courge2017.rare5k$trt <- factor(metadata.courge2017.rare5k$trt)
metadata.courge2017.rare5k$type <- factor(metadata.courge2017.rare5k$type)
metadata.courge2017.rare5k$block <- factor(metadata.courge2017.rare5k$block)
```

##NMDS for all date and treatments
####2016
```{r nmds 2016, message = FALSE}
asvs.comm.courge2016.rare5k.mds <- metaMDS(decostand(asvs.comm.courge2016.rare5k,method="hellinger"))
```

###2017
```{r nmds 2017, message = FALSE}
asvs.comm.courge2017.rare5k.mds <- metaMDS(decostand(asvs.comm.courge2017.rare5k,method="hellinger"))
```

##(line 270) Supplemental Figure 8: ordination of squash phyllosphere community during the growing season of 2016 and 2017 at the 3 sampling dates#
###2016
```{r ordiplot comm vs date 2016}
par(mfrow=c(1,2))

myGreyCol <- c("grey30", "grey50", "grey70")
myGreyColAlpha <- makeTransparent(myGreyCol, alpha = 0.3)
                    

ordiplot(asvs.comm.courge2016.rare5k.mds, display="sites", type = 'n', main = "2016")

with(metadata.courge2016.rare5k, ordiellipse(asvs.comm.courge2016.rare5k.mds, 
                                                   group = date, kind = "sd", 
                                                   label = TRUE,
                                                   col = myGreyCol))

with(metadata.courge2016.rare5k, points(asvs.comm.courge2016.rare5k.mds,
                                                display = "sites", col = myGreyCol[as.numeric(metadata.courge2016.rare5k$date)],
                                                pch = 21, bg = myGreyCol[as.numeric(metadata.courge2016.rare5k$date)]))

with(metadata.courge2016.rare5k, ordispider(asvs.comm.courge2016.rare5k.mds, 
                                                   group = date,
                                                   label = TRUE,
                                                   col = myGreyColAlpha))

plot.new()
with(metadata.courge2016.rare5k, legend("topleft", legend = c("Date 1","Date 2", "Date 3"), bty = "n",
                                    pt.cex=2, cex=1.5, col = myGreyCol,
                                    pch = 21, pt.bg = myGreyCol))
```

###2017
```{r ordiplot comm vs date 2017}
par(mfrow=c(1,2))

myGreyCol <- c("grey30", "grey50", "grey70")
myGreyColAlpha <- makeTransparent(myGreyCol, alpha = 0.3)
                    

ordiplot(asvs.comm.courge2017.rare5k.mds, display="sites", type = 'n', main = "2017")

with(metadata.courge2017.rare5k, ordiellipse(asvs.comm.courge2017.rare5k.mds, 
                                                   group = date, kind = "sd", 
                                                   label = TRUE,
                                                   col = myGreyCol))

with(metadata.courge2017.rare5k, points(asvs.comm.courge2017.rare5k.mds,
                                                display = "sites", col = myGreyCol[as.numeric(metadata.courge2017.rare5k$date)],
                                                pch = 21, bg = myGreyCol[as.numeric(metadata.courge2017.rare5k$date)]))

with(metadata.courge2017.rare5k, ordispider(asvs.comm.courge2017.rare5k.mds, 
                                                   group = date,
                                                   label = TRUE,
                                                   col = myGreyColAlpha))
```


##PERMANOVA on rarefied data : trt vs date
###(line 270) 2016
```{r permanova ~ block/trt * date 2016}
adonis(decostand(asvs.comm.courge2016.rare5k, method="hellinger") ~ trt * date, strata=metadata.courge2016.rare5k$block, data=metadata.courge2016.rare5k)
```

###(line 271) 2017
```{r permanova ~ block/trt * date 2017}
adonis(decostand(asvs.comm.courge2017.rare5k, method="hellinger") ~ trt *date, strata = metadata.courge2017.rare5k$block, data=metadata.courge2017.rare5k)
```

## Modeling community diversity as function of treatment
###(Supp line 48) 2016
```{r linear mixed model of 2016 community diversity }
metadata.courge2016.rare5k.lmer <- lmerTest::lmer(alphaDiv ~ trt + (1|block), data = metadata.courge2016.rare5k)
summary(metadata.courge2016.rare5k.lmer)
```
##### no effect of blocking variable in 2016, perform linear model instead
```{r}
trt.allDate.div.aov <- aov(vegan::diversity(asvs.comm.courge2016.rare5k) ~ trt, data=metadata.courge2016.rare5k)
summary(trt.allDate.div.aov)
TukeyHSD(trt.allDate.div.aov, ordered = TRUE)
```

### (Supp line 49) 2017
```{r linear mixed model of 2017 community diversity}
metadata.courge2017.rare5k.lmer <- lmerTest::lmer(alphaDiv ~ trt + (1|block), data = metadata.courge2017.rare5k)
emmeans(metadata.courge2017.rare5k.lmer, list(pairwise ~ trt), adjust = "tukey")
```

##(line 268) supplemental Figure 6 
###2016 top left figure
```{r diversity plot by treatment during 2016 growing season}
metadata.courge2016.rare5k$alphaDiv <- c(vegan::diversity(asvs.comm.courge2016.rare5k))
ggplot(metadata.courge2016.rare5k, aes(x=trt, y=alphaDiv))+ geom_jitter(aes(alpha=0.5))+geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2016") + theme_bw() + xlab(NULL) + ylab("Alpha Diversity (Shannon)") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="red", linetype = 1) +  scale_x_discrete(labels = c('PC','RCC', 'CT-RCC', 'BS')) 
```

###2017 top right figure
```{r diversity plot by treatment during 2017 growing season}
metadata.courge2017.rare5k$alphaDiv <- c(vegan::diversity(asvs.comm.courge2017.rare5k))
ggplot(metadata.courge2017.rare5k, aes(x=trt, y=alphaDiv))+ geom_jitter(alpha=0.5) +geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2017") + theme_bw() + xlab(NULL) + ylab("Alpha Diversity (Shannon)") + 
  scale_x_discrete(labels = c('PC','RCC', 'CT-RCC', 'BS')) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="red", linetype = 1)
```

### 2016 bottom left figure
```{r diversity plot by date during 2016 growing season}
ggplot(metadata.courge2016.rare5k, aes(x=date, y=alphaDiv))+ geom_jitter(aes(alpha=0.5)) +geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2016") + theme_bw() + xlab("Season sampling") +  ylab("Alpha Diversity (Shannon)")+ stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.4, col="red", linetype = 1) + scale_x_discrete(labels= c("Early", "Mid", "Late"))
```

### 2017 bottom right figure
```{r diversity plot by date during 2017 growing season}
ggplot(metadata.courge2017.rare5k, aes(x=date, y=alphaDiv))+ geom_jitter(aes(alpha=0.5)) +geom_violin(alpha=0.7) + guides(alpha=FALSE) +
  ggtitle("2017") + theme_bw() + xlab(NULL) + ylab("Alpha Diversity (Shannon)") + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="red", linetype = 1) + scale_x_discrete(labels= c("Early", "Mid", "Late"))
```

## Modeling community diversity as function of sampling date
### 2016(Supp line 51)
```{r diversity courge2016 * treatment ANOVA }
date.allTrt.div.lm <- lm(vegan::diversity(asvs.comm.courge2016.rare5k) ~  date , data=metadata.courge2016.rare5k)
anova(date.allTrt.div.lm)
```
### 2017(Supp line 51)
```{r diversity courge2017 * treatment ANOVA }
date.allTrt.div.lm <- lm(vegan::diversity(asvs.comm.courge2017.rare5k) ~ date , data=metadata.courge2017.rare5k)
anova(date.allTrt.div.lm)
```

## Tukey post-hoc test on community diversity as function of sampling date model
### 2016 (Supp line 53)
```{r div date aov & TukeyHSD 2016}
date.Trt.div.aov <- aov(vegan::diversity(asvs.comm.courge2016.rare5k) ~ date, data=metadata.courge2016.rare5k)
TukeyHSD(date.Trt.div.aov, ordered = TRUE)
```
### 2017 (Supp line 53)
```{r div date aov & TukeyHSD 2017}
date.Trt.div.aov <- aov(vegan::diversity(asvs.comm.courge2017.rare5k) ~ date, data=metadata.courge2017.rare5k)
TukeyHSD(date.Trt.div.aov, ordered = TRUE)
```

##(line 268) supplemental Figure 7
### 2016 top left figure
```{r specnumber plot by date with all trt 2016}
metadata.courge2016.rare5k$specnumber <- specnumber(asvs.comm.courge2016.rare5k)
ggplot(metadata.courge2016.rare5k, aes(x=trt, y=specnumber))+ geom_jitter(aes(alpha=0.5)) +geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2016") + theme_bw() + xlab(NULL) + ylab("ASVs number")+ stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.45, col="red", linetype = 1) +  scale_x_discrete(labels = c('PC','RCC', 'CT-RCC', 'BS'))
```

### 2017 top right figure
```{r specnumber plot by date with all trt 2017}
metadata.courge2017.rare5k$specnumber <- specnumber(asvs.comm.courge2017.rare5k)
ggplot(metadata.courge2017.rare5k, aes(x=trt, y=specnumber))+ geom_jitter(alpha=0.5) +geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2017") + theme_bw() + xlab(NULL) + ylab("ASVs number") + 
  scale_x_discrete(labels = c('PC','RCC', 'CT-RCC', 'BS'))  + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5, col="red", linetype = 1)
```

##(Supp line 59)
### 2016
```{r spcnumber 2016 trt aov & TukeyHSD}
spcnmb.aov <- aov(specnumber(asvs.comm.courge2016.rare5k) ~ trt, data=metadata.courge2016.rare5k)
TukeyHSD(spcnmb.aov, ordered = TRUE)
```

### 2017
```{r spcnumber 2017 trt aov & TukeyHSD}
metadata.courge2017.rare5k.lmer <- lmerTest::lmer(specnumber ~ trt + (1|block), data = metadata.courge2017.rare5k)
emmeans(metadata.courge2017.rare5k.lmer, list(pairwise ~ trt), adjust = "tukey")
```

### 2016 bottom left figure
```{r}
ggplot(metadata.courge2016.rare5k, aes(x=date, y=specnumber))+ geom_jitter(aes(alpha=0.5)) +geom_violin(alpha=0.7) + guides(alpha=FALSE) + ggtitle("2016") + theme_bw() + xlab("Season sampling") + ylab("ASVs number")+ stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.45, col="red", linetype = 1) 
```

### 2017 bottom right figure
```{r}
ggplot(metadata.courge2017.rare5k, aes(x=date, y=specnumber))+ geom_jitter(aes(alpha=0.5)) +geom_violin(alpha=0.7) + guides(alpha=FALSE) +
  ggtitle("2017") + theme_bw() + xlab(NULL) + ylab("ASVs number") + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.3, col="red", linetype = 1) + scale_x_discrete(labels= c("Early", "Mid", "Late"))
```

##(Supp line 59)
### 2016
```{r spc 2016 date aov & TukeyHSD}
spcnmb.date.aov <- aov(specnumber(asvs.comm.courge2016.rare5k) ~ date, data=metadata.courge2016.rare5k)
anova(spcnmb.date.aov)
TukeyHSD(spcnmb.date.aov, ordered = TRUE)
```

### 2017
```{r spc 2017 date lmer & TukeyHSD}
metadata.courge2017.rare5k.lmer <- lmerTest::lmer(specnumber ~ date + (1|block), data = metadata.courge2017.rare5k)
emmeans(metadata.courge2017.rare5k.lmer, list(pairwise ~ date), adjust = "tukey")
```


#############################################
#Effect of treatments on community diversity#
#############################################

##Split matrix per season sampling
###2016
```{r split courge2016 based on Date 1, message = FALSE}
metadata.courge2016.rare5k.D1 <- metadata.courge2016.rare5k[metadata.courge2016.rare5k$date == 'Date 1',]
metadata.courge2016.rare5k.D1$date <- factor(metadata.courge2016.rare5k.D1$date)
metadata.courge2016.rare5k.D1$trt <- factor(metadata.courge2016.rare5k.D1$trt)
metadata.courge2016.rare5k.D1$type <- factor(metadata.courge2016.rare5k.D1$type)
metadata.courge2016.rare5k.D1$parcelle <- factor(metadata.courge2016.rare5k.D1$parcelle)
metadata.courge2016.rare5k.D1$block <- factor(metadata.courge2016.rare5k.D1$block)
metadata.courge2016.rare5k.D1$sub_ech <- factor(metadata.courge2016.rare5k.D1$sub_ech)
metadata.courge2016.rare5k.D1$run <- factor(metadata.courge2016.rare5k.D1$run)
metadata.courge2016.rare5k.D1$poids <- as.numeric(metadata.courge2016.rare5k.D1$poids)
metadata.courge2016.rare5k.D1$bact_count <- as.numeric(metadata.courge2016.rare5k.D1$bact_count)
metadata.courge2016.rare5k.D1$pseudo_count <- as.numeric(metadata.courge2016.rare5k.D1$pseudo_count)
asvs.comm.courge2016.rare5k.D1 <- asvs.comm.courge2016.rare5k[rownames(metadata.courge2016.rare5k.D1),]
asvs.taxo.courge2016.rare5k.D1 <- asvs.taxo.courge2016.rare5k[colnames(asvs.comm.courge2016.rare5k.D1),]
```
```{r split courge2016 based on Date 2, message = FALSE}
metadata.courge2016.rare5k.D2 <- metadata.courge2016.rare5k[metadata.courge2016.rare5k$date == 'Date 2',]
metadata.courge2016.rare5k.D2$date <- factor(metadata.courge2016.rare5k.D2$date)
metadata.courge2016.rare5k.D2$trt <- factor(metadata.courge2016.rare5k.D2$trt)
metadata.courge2016.rare5k.D2$type <- factor(metadata.courge2016.rare5k.D2$type)
metadata.courge2016.rare5k.D2$parcelle <- factor(metadata.courge2016.rare5k.D2$parcelle)
metadata.courge2016.rare5k.D2$block <- factor(metadata.courge2016.rare5k.D2$block)
metadata.courge2016.rare5k.D2$sub_ech <- factor(metadata.courge2016.rare5k.D2$sub_ech)
metadata.courge2016.rare5k.D2$run <- factor(metadata.courge2016.rare5k.D2$run)
metadata.courge2016.rare5k.D2$poids <- as.numeric(metadata.courge2016.rare5k.D2$poids)
metadata.courge2016.rare5k.D2$bact_count <- as.numeric(metadata.courge2016.rare5k.D2$bact_count)
metadata.courge2016.rare5k.D2$pseudo_count <- as.numeric(metadata.courge2016.rare5k.D2$pseudo_count)
asvs.comm.courge2016.rare5k.D2 <- asvs.comm.courge2016.rare5k[rownames(metadata.courge2016.rare5k.D2),]
asvs.taxo.courge2016.rare5k.D2 <- asvs.taxo.courge2016.rare5k[colnames(asvs.comm.courge2016.rare5k.D2),]
```
```{r split courge2016 based on Date 3, message = FALSE}
metadata.courge2016.rare5k.D3 <- metadata.courge2016.rare5k[metadata.courge2016.rare5k$date == 'Date 3',]
metadata.courge2016.rare5k.D3$date <- factor(metadata.courge2016.rare5k.D3$date)
metadata.courge2016.rare5k.D3$trt <- factor(metadata.courge2016.rare5k.D3$trt)
metadata.courge2016.rare5k.D3$type <- factor(metadata.courge2016.rare5k.D3$type)
metadata.courge2016.rare5k.D3$parcelle <- factor(metadata.courge2016.rare5k.D3$parcelle)
metadata.courge2016.rare5k.D3$block <- factor(metadata.courge2016.rare5k.D3$block)
metadata.courge2016.rare5k.D3$sub_ech <- factor(metadata.courge2016.rare5k.D3$sub_ech)
metadata.courge2016.rare5k.D3$run <- factor(metadata.courge2016.rare5k.D3$run)
metadata.courge2016.rare5k.D3$poids <- as.numeric(metadata.courge2016.rare5k.D3$poids)
metadata.courge2016.rare5k.D3$bact_count <- as.numeric(metadata.courge2016.rare5k.D3$bact_count)
metadata.courge2016.rare5k.D3$pseudo_count <- as.numeric(metadata.courge2016.rare5k.D3$pseudo_count)
asvs.comm.courge2016.rare5k.D3 <- asvs.comm.courge2016.rare5k[rownames(metadata.courge2016.rare5k.D3),]
asvs.taxo.courge2016.rare5k.D3 <- asvs.taxo.courge2016.rare5k[colnames(asvs.comm.courge2016.rare5k.D3),]
```
###2017
```{r split courge2017 based on Date 1, message = FALSE}
metadata.courge2017.rare5k.D1 <- metadata.courge2017.rare5k[metadata.courge2017.rare5k$date == 'Date 1',]
metadata.courge2017.rare5k.D1$date <- factor(metadata.courge2017.rare5k.D1$date)
metadata.courge2017.rare5k.D1$trt <- factor(metadata.courge2017.rare5k.D1$trt)
metadata.courge2017.rare5k.D1$type <- factor(metadata.courge2017.rare5k.D1$type)
metadata.courge2017.rare5k.D1$parcelle <- factor(metadata.courge2017.rare5k.D1$parcelle)
metadata.courge2017.rare5k.D1$block <- factor(metadata.courge2017.rare5k.D1$block)
metadata.courge2017.rare5k.D1$sub_ech <- factor(metadata.courge2017.rare5k.D1$sub_ech)
metadata.courge2017.rare5k.D1$run <- factor(metadata.courge2017.rare5k.D1$run)
metadata.courge2017.rare5k.D1$poids <- as.numeric(metadata.courge2017.rare5k.D1$poids)
metadata.courge2017.rare5k.D1$log_cfu_tsa <- as.numeric(metadata.courge2017.rare5k.D1$log_cfu_tsa)
metadata.courge2017.rare5k.D1$log_cfu_rif <- as.numeric(metadata.courge2017.rare5k.D1$log_cfu_rif)
asvs.comm.courge2017.rare5k.D1 <- asvs.comm.courge2017.rare5k[rownames(metadata.courge2017.rare5k.D1),]
asvs.comm.courge2017.rare5k.D1 <- asvs.comm.courge2017.rare5k.D1[,apply(asvs.comm.courge2017.rare5k.D1,2,sum)>0]
asvs.taxo.courge2017.rare5k.D1 <- asvs.taxo.courge2017.rare5k[colnames(asvs.comm.courge2017.rare5k.D1),]
```
```{r split courge2017 based on Date 2, message = FALSE}
metadata.courge2017.rare5k.D2 <- metadata.courge2017.rare5k[metadata.courge2017.rare5k$date == 'Date 2',]
metadata.courge2017.rare5k.D2$date <- factor(metadata.courge2017.rare5k.D2$date)
metadata.courge2017.rare5k.D2$trt <- factor(metadata.courge2017.rare5k.D2$trt)
metadata.courge2017.rare5k.D2$type <- factor(metadata.courge2017.rare5k.D2$type)
metadata.courge2017.rare5k.D2$parcelle <- factor(metadata.courge2017.rare5k.D2$parcelle)
metadata.courge2017.rare5k.D2$block <- factor(metadata.courge2017.rare5k.D2$block)
metadata.courge2017.rare5k.D2$sub_ech <- factor(metadata.courge2017.rare5k.D2$sub_ech)
metadata.courge2017.rare5k.D2$run <- factor(metadata.courge2017.rare5k.D2$run)
metadata.courge2017.rare5k.D2$poids <- as.numeric(metadata.courge2017.rare5k.D2$poids)
metadata.courge2017.rare5k.D2$log_cfu_tsa <- as.numeric(metadata.courge2017.rare5k.D2$log_cfu_tsa)
metadata.courge2017.rare5k.D2$log_cfu_rif <- as.numeric(metadata.courge2017.rare5k.D2$log_cfu_rif)
asvs.comm.courge2017.rare5k.D2 <- asvs.comm.courge2017.rare5k[rownames(metadata.courge2017.rare5k.D2),]
asvs.comm.courge2017.rare5k.D2 <- asvs.comm.courge2017.rare5k.D2[,apply(asvs.comm.courge2017.rare5k.D2,2,sum)>0]
asvs.taxo.courge2017.rare5k.D2 <- asvs.taxo.courge2017.rare5k[colnames(asvs.comm.courge2017.rare5k.D2),]
```
```{r split courge2017 based on Date 3, message = FALSE}
metadata.courge2017.rare5k.D3 <- metadata.courge2017.rare5k[metadata.courge2017.rare5k$date == 'Date 3',]
metadata.courge2017.rare5k.D3$date <- factor(metadata.courge2017.rare5k.D3$date)
metadata.courge2017.rare5k.D3$trt <- factor(metadata.courge2017.rare5k.D3$trt)
metadata.courge2017.rare5k.D3$type <- factor(metadata.courge2017.rare5k.D3$type)
metadata.courge2017.rare5k.D3$parcelle <- factor(metadata.courge2017.rare5k.D3$parcelle)
metadata.courge2017.rare5k.D3$block <- factor(metadata.courge2017.rare5k.D3$block)
metadata.courge2017.rare5k.D3$sub_ech <- factor(metadata.courge2017.rare5k.D3$sub_ech)
metadata.courge2017.rare5k.D3$run <- factor(metadata.courge2017.rare5k.D3$run)
metadata.courge2017.rare5k.D3$poids <- as.numeric(metadata.courge2017.rare5k.D3$poids)
metadata.courge2017.rare5k.D3$log_cfu_tsa <- as.numeric(metadata.courge2017.rare5k.D3$log_cfu_tsa)
metadata.courge2017.rare5k.D3$log_cfu_rif <- as.numeric(metadata.courge2017.rare5k.D3$log_cfu_rif)
asvs.comm.courge2017.rare5k.D3 <- asvs.comm.courge2017.rare5k[rownames(metadata.courge2017.rare5k.D3),]
asvs.comm.courge2017.rare5k.D3 <- asvs.comm.courge2017.rare5k.D3[,apply(asvs.comm.courge2017.rare5k.D3,2,sum)>0]
asvs.taxo.courge2017.rare5k.D3 <- asvs.taxo.courge2017.rare5k[colnames(asvs.comm.courge2017.rare5k.D3),]
```

##Perform NMDS
###2016
```{r nmds courge2016, message = FALSE}
asvs.comm.courge2016.rare5k.D1.mds <- metaMDS(decostand(asvs.comm.courge2016.rare5k.D1,method="hellinger"))
asvs.comm.courge2016.rare5k.D2.mds <- metaMDS(decostand(asvs.comm.courge2016.rare5k.D2,method="hellinger"))
asvs.comm.courge2016.rare5k.D3.mds <- metaMDS(decostand(asvs.comm.courge2016.rare5k.D3,method="hellinger"))
```
###2017
```{r nmds courge2017, message = FALSE}
asvs.comm.courge2017.rare5k.D1.mds <- metaMDS(decostand(asvs.comm.courge2017.rare5k.D1,method="hellinger"))
asvs.comm.courge2017.rare5k.D2.mds <- metaMDS(decostand(asvs.comm.courge2017.rare5k.D2,method="hellinger"))
asvs.comm.courge2017.rare5k.D3.mds <- metaMDS(decostand(asvs.comm.courge2017.rare5k.D3,method="hellinger"))
```

## (line 279) Linear model : alpha diversity ~ trt
###2016
```{r linear model alphaDiv per date, echo=FALSE}
multiplot(ggplot(metadata.courge2016.rare5k.re.D1, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2016.rare5k.re.D1$alphaDiv))), 
          ggplot(metadata.courge2016.rare5k.re.D2, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2016.rare5k.re.D2$alphaDiv))),
          ggplot(metadata.courge2016.rare5k.re.D3, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2016.rare5k.re.D3$alphaDiv))))


metadata.courge2016.rare5k.D1.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2016.rare5k.re.D1)
metadata.courge2016.rare5k.D2.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2016.rare5k.re.D2)
metadata.courge2016.rare5k.D3.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2016.rare5k.re.D3)

par(mfrow=c(1,3))
qqnorm(resid(metadata.courge2016.rare5k.D1.ad.lm), main = "Date 1")
qqline(resid(metadata.courge2016.rare5k.D1.ad.lm))

qqnorm(resid(metadata.courge2016.rare5k.D2.ad.lm), main = "Date 2")
qqline(resid(metadata.courge2016.rare5k.D2.ad.lm))

qqnorm(resid(metadata.courge2016.rare5k.D3.ad.lm), main = "Date 3")
qqline(resid(metadata.courge2016.rare5k.D3.ad.lm))

summary(metadata.courge2016.rare5k.D1.ad.lm)
summary(metadata.courge2016.rare5k.D2.ad.lm)
summary(metadata.courge2016.rare5k.D3.ad.lm)
```

#### (line 281) TukeyHSD on Linear model 
```{r}
metadata.courge2016.rare5k.D1.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2016.rare5k.D1.ad.lm))
metadata.courge2016.rare5k.D2.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2016.rare5k.D2.ad.lm))
metadata.courge2016.rare5k.D3.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2016.rare5k.D3.ad.lm))

metadata.courge2016.rare5k.D1.ad.lm.aov.tk
metadata.courge2016.rare5k.D2.ad.lm.aov.tk
metadata.courge2016.rare5k.D3.ad.lm.aov.tk
```

|                    trt                    |signif group (date1)|signif group (date2)|signif group (date3)|
| ----------------------------------------- |:------------------:|:------------------:|:------------------:|
| `r unique(metadata.courge2016.d1$trt)[1]` |            b       |          a         |          a         |#Sol_nu
| `r unique(metadata.courge2016.d1$trt)[2]` |            b       |          a         |          a         |#Plastique
| `r unique(metadata.courge2016.d1$trt)[3]` |           a        |          a         |          a         |#Seigle
| `r unique(metadata.courge2016.d1$trt)[4]` |           a        |          a         |          a         |#Seigle_Bru

###2017
```{r linear alphaDiv per date}
multiplot(ggplot(metadata.courge2017.rare5k.re.D1, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2017.rare5k.re.D1$alphaDiv))), 
          ggplot(metadata.courge2017.rare5k.re.D2, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2017.rare5k.re.D2$alphaDiv))),
          ggplot(metadata.courge2017.rare5k.re.D3, aes(alphaDiv))+ geom_density()+xlim(c(0, max(metadata.courge2017.rare5k.re.D3$alphaDiv))))


metadata.courge2017.rare5k.D1.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2017.rare5k.re.D1)
metadata.courge2017.rare5k.D2.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2017.rare5k.re.D2)
metadata.courge2017.rare5k.D3.ad.lm <- lm(alphaDiv ~ trt, data = metadata.courge2017.rare5k.re.D3)

par(mfrow=c(1,3))
qqnorm(resid(metadata.courge2017.rare5k.D1.ad.lm), main = "Date 1")
qqline(resid(metadata.courge2017.rare5k.D1.ad.lm))

qqnorm(resid(metadata.courge2017.rare5k.D2.ad.lm), main = "Date 2")
qqline(resid(metadata.courge2017.rare5k.D2.ad.lm))

qqnorm(resid(metadata.courge2017.rare5k.D3.ad.lm), main = "Date 3")
qqline(resid(metadata.courge2017.rare5k.D3.ad.lm))

summary(metadata.courge2017.rare5k.D1.ad.lm)
summary(metadata.courge2017.rare5k.D2.ad.lm)
summary(metadata.courge2017.rare5k.D3.ad.lm)
```


```{r}
metadata.courge2017.rare5k.D1.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2017.rare5k.D1.ad.lm))
metadata.courge2017.rare5k.D2.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2017.rare5k.D2.ad.lm))
metadata.courge2017.rare5k.D3.ad.lm.aov.tk <- TukeyHSD(aov(metadata.courge2017.rare5k.D3.ad.lm))

metadata.courge2017.rare5k.D1.ad.lm.aov.tk
metadata.courge2017.rare5k.D2.ad.lm.aov.tk
metadata.courge2017.rare5k.D3.ad.lm.aov.tk
```

|                       trt                        |signif group (date1)|signif group (date2) |signif group (date3)|
| ------------------------------------------------ |:------------------:|:-------------------:|:------------------:|
| `r unique(metadata.courge2017.rare5k.D1$trt)[1]` |           a        |          ab         |          a         |#Seigle_bru
| `r unique(metadata.courge2017.rare5k.D1$trt)[2]` |           a        |            c        |          a         |#Sol_nu
| `r unique(metadata.courge2017.rare5k.D1$trt)[3]` |           a        |          ab         |          a         |#Seigle
| `r unique(metadata.courge2017.rare5k.D1$trt)[4]` |            b       |           bc        |          a         |#Plastique

## (line 279)Figure 2
###2016 top panel
```{r alphaDiv per date 2016}
date_name <- c(`Date 1` = "Early season",
               `Date 2` = "Mid season",
               `Date 3` = "Late season")


metadata.courge2016.rare5k$alphaDiv <- vegan::diversity(asvs.comm.courge2016.rare5k)

metadata.courge2016.rare5k.D1$alphaDiv <- vegan::diversity(asvs.comm.courge2016.rare5k.D1)
metadata.courge2016.rare5k.D2$alphaDiv <- vegan::diversity(asvs.comm.courge2016.rare5k.D2)
metadata.courge2016.rare5k.D3$alphaDiv <- vegan::diversity(asvs.comm.courge2016.rare5k.D3)

metadata.courge2016.rare5k.re.D1 <- within(metadata.courge2016.rare5k.D1, trt <- relevel(trt, ref = "Sol_nu"))
metadata.courge2016.rare5k.re.D2 <- within(metadata.courge2016.rare5k.D2, trt <- relevel(trt, ref = "Sol_nu"))
metadata.courge2016.rare5k.re.D3 <- within(metadata.courge2016.rare5k.D3, trt <- relevel(trt, ref = "Sol_nu"))

violin_legend <- ggplot(metadata.courge2016.rare5k.re.D1, aes(x=trt, y=alphaDiv, fill=trt)) + geom_violin(alpha=0.7)+ scale_fill_manual(values = c("#0000FF", "#006400", "#800000", "#FFA500"), name="Treatment", labels=c("PC","RCC", "CT-RCC", "NS")) +  guides(alpha=FALSE)


g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- g_legend(violin_legend) 

grid.newpage()
grid.draw(legend) 


labelDiv.pos <- max(metadata.courge2016.rare5k$alphaDiv) + 5*(max(metadata.courge2016.rare5k$alphaDiv))/100

trt2016 <- c("Sol_nu","Plastique", "Seigle_Bru", "Seigle")
posDiv2016 <- c(labelDiv.pos, labelDiv.pos, labelDiv.pos, labelDiv.pos)
date1_2016 <- c("Date 1", "Date 1", "Date 1", "Date 1")
date2_2016 <- c("Date 2", "Date 2", "Date 2", "Date 2")
date3_2016 <- c("Date 3", "Date 3", "Date 3", "Date 3")
signifDiv.2016.d1 <- c("b", "b", "a", "a")
signifDiv.2016.d2 <- c("a", "a", "a", "a")
signifDiv.2016.d3 <- c("a", "a", "a", "a")

labelDiv2016.df <- data.frame(trt=c(trt2016, trt2016, trt2016), 
                       y=c(posDiv2016, posDiv2016, posDiv2016), 
                       label=c(signifDiv.2016.d1, signifDiv.2016.d2, signifDiv.2016.d3),
                       date= c(date1_2016, date2_2016, date3_2016))

alphadiv2016_final.gg<-ggplot(metadata.courge2016.rare5k, aes(x=trt, y=alphaDiv, fill=trt)) + geom_jitter(aes(alpha=0.5)) + facet_grid(.~date,  labeller = as_labeller(date_name)) + geom_violin(alpha=0.7) +
  ggtitle("2016") + theme_bw() + xlab(NULL) + ylab("Alpha Diversity (shannon Index)") + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.3, col="red", linetype = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("#0000FF", "#006400", "#800000", "#FFA500"), labels=c("PC","RCC", "CT-RCC", "BS"))+ 
  guides(alpha=F, fill=F) + geom_text(data=labelDiv2016.df, aes(x=trt, y=y, label=label)) 
alphadiv2016_final.gg
```

###2017 bottom panel
```{r alphaDiv per date 2017}
metadata.courge2017.rare5k.re.D1 <- within(metadata.courge2017.rare5k.D1, trt <- relevel(trt, ref = "Sol_nu"))
metadata.courge2017.rare5k.re.D2 <- within(metadata.courge2017.rare5k.D2, trt <- relevel(trt, ref = "Sol_nu"))
metadata.courge2017.rare5k.re.D3 <- within(metadata.courge2017.rare5k.D3, trt <- relevel(trt, ref = "Sol_nu"))

date_name <- c(`Date 1` = "Early season",
               `Date 2` = "Mid season",
               `Date 3` = "Late season")

violin_legend <- ggplot(metadata.courge2017.rare5k, aes(x=trt, y=alphaDiv, fill=trt)) + geom_violin(alpha=0.7)+ scale_fill_manual(values = c("#0000FF", "#006400", "#800000", "#FFA500"), name="Treatment", labels=c("PC","RCC", "CT-RCC", "BS")) +  guides(alpha=FALSE)


g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- g_legend(violin_legend) 

grid.newpage()
grid.draw(legend) 

labelDiv2017.pos <- max(metadata.courge2017.rare5k$alphaDiv) + 5*(max(metadata.courge2017.rare5k$alphaDiv))/100

trt2017 <- c("Seigle_bru","Sol_nu", "Seigle", "Plastique")
posDiv2017 <- c(labelDiv2017.pos, labelDiv2017.pos, labelDiv2017.pos, labelDiv2017.pos)
date1_2017 <- c("Date 1", "Date 1", "Date 1", "Date 1")
date2_2017 <- c("Date 2", "Date 2", "Date 2", "Date 2")
date3_2017 <- c("Date 3", "Date 3", "Date 3", "Date 3")
signifDiv.2017.d1 <- c("a", "a", "a", "b")
signifDiv.2017.d2 <- c("ab", "c", "ab", "bc")
signifDiv.2017.d3 <- c("a", "a", "a", "a")

labelDiv2017.df <- data.frame(trt=c(trt2017, trt2017, trt2017), 
                       y=c(posDiv2017, posDiv2017, posDiv2017), 
                       label=c(signifDiv.2017.d1, signifDiv.2017.d2, signifDiv.2017.d3),
                       date= c(date1_2017, date2_2017, date3_2017))


alphadiv2017_final.gg <- ggplot(metadata.courge2017.rare5k, aes(x=trt, y=alphaDiv, fill=trt))+ geom_jitter(aes(alpha=0.5))+ facet_grid(.~date, labeller = as_labeller(date_name))  +geom_violin(alpha=0.7)  +
  ggtitle("2017") + theme_bw() + xlab(NULL) + ylab("Alpha Diversity (Shannon Index)") + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.3, col="red", linetype = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("#0000FF", "#006400", "#800000", "#FFA500"), labels=c("PL","RCC", "CT-RCC", "BS"))+ 
  guides(alpha=F, fill=F) + geom_text(data=labelDiv2017.df, aes(x=trt, y=y, label=label)) 

alphadiv2017_final.gg
```

## (line 291) PERMANOVA
###2016
```{r permanova 2016 per date}
asvs.comm.courge2016.rare5k.D1.perm <-adonis(decostand(asvs.comm.courge2016.rare5k.D1, method="hellinger") ~ trt, strata=metadata.courge2016.rare5k.D1$block, data=metadata.courge2016.rare5k.D1, permutations = 9999)

paste("Date 1","trt:R2 =",sprintf("%.3f",asvs.comm.courge2016.rare5k.D1.perm$aov.tab$R2[1]), 
      "P =", asvs.comm.courge2016.rare5k.D1.perm$aov.tab$`Pr(>F)`[1], collapse = " ")

asvs.comm.courge2016.rare5k.D2.perm<-adonis(decostand(asvs.comm.courge2016.rare5k.D2, method="hellinger") ~ trt, strata=metadata.courge2016.rare5k.D2$block, data=metadata.courge2016.rare5k.D2, permutations = 9999)

paste("Date 2","trt:R2 =",sprintf("%.3f",asvs.comm.courge2016.rare5k.D2.perm$aov.tab$R2[1]), 
      "P =", asvs.comm.courge2016.rare5k.D2.perm$aov.tab$`Pr(>F)`[1], collapse = " ")

asvs.comm.courge2016.rare5k.D3.perm <- adonis(decostand(asvs.comm.courge2016.rare5k.D3, method="hellinger") ~ trt, strata=metadata.courge2016.rare5k.D3$block, data=metadata.courge2016.rare5k.D3, permutations = 9999)

paste("Date 3","trt:R2 =",sprintf("%.3f",asvs.comm.courge2016.rare5k.D3.perm$aov.tab$R2[1]), 
      "P =", asvs.comm.courge2016.rare5k.D3.perm$aov.tab$`Pr(>F)`[1], collapse = " ")

```

###2017
```{r permanova 2017 per date, message = FALSE}
asvs.comm.courge2017.rare5k.D1.perm <-adonis(decostand(asvs.comm.courge2017.rare5k.D1, method="hellinger") ~ trt, strata=metadata.courge2017.rare5k.D1$block, data=metadata.courge2017.rare5k.D1)

paste("Date 1","trt:R2 =",sprintf("%.3f",asvs.comm.courge2017.rare5k.D1.perm$aov.tab$R2[1]), 
       "P =", asvs.comm.courge2017.rare5k.D1.perm$aov.tab$`Pr(>F)`[1], collapse = " ")


asvs.comm.courge2017.rare5k.D2.perm<-adonis(decostand(asvs.comm.courge2017.rare5k.D2, method="hellinger") ~ trt, strata=metadata.courge2017.rare5k.D2$block, data=metadata.courge2017.rare5k.D2)

paste("Date 2","trt:R2 =",sprintf("%.3f",asvs.comm.courge2017.rare5k.D2.perm$aov.tab$R2[1]), 
       "P =", asvs.comm.courge2017.rare5k.D2.perm$aov.tab$`Pr(>F)`[1], collapse = " ")


asvs.comm.courge2017.rare5k.D3.perm <- adonis(decostand(asvs.comm.courge2017.rare5k.D3, method="hellinger") ~ trt, strata=metadata.courge2017.rare5k.D3$block, data=metadata.courge2017.rare5k.D3)

paste("Date 3","trt:R2 =",sprintf("%.3f",asvs.comm.courge2017.rare5k.D3.perm$aov.tab$R2[1]), 
       "P =", asvs.comm.courge2017.rare5k.D3.perm$aov.tab$`Pr(>F)`[1], collapse = " ")
```

## (line 292) Figure 3
###2016 top panel
```{r Figure 3 2016}
legend.lvl <- data.frame(freq=metadata.courge2016.rare5k[order(metadata.courge2016.rare5k$pseudo_count), "pseudo_count"])
legend.lvl.smry<-summary(legend.lvl)
legend.lvl <- data.frame(freq=c(as.numeric(strsplit(legend.lvl.smry, split = ":")[[1]][2]),
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[2]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[3]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[4]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[5]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[6]][2])
                                )
                         )
par(mfrow=c(1,4)) 
colvec <- c("#0000FF", "#006400", "#800000", "#FFA500")
#date1
ordiplot(asvs.comm.courge2016.rare5k.D1.mds, display="sites", type = 'n', main = "Early season 2016")

with(metadata.courge2016.rare5k.D1, ordispider(asvs.comm.courge2016.rare5k.D1.mds, 
                                                   group = trt, label = FALSE,
                                                   col = c("#e1e1ff", "#bfffbf", "#ffd0d0",  "#ffecc9")
                                               )
     )
with(metadata.courge2016.rare5k.D1, ordiellipse(asvs.comm.courge2016.rare5k.D1.mds, 
                                                   group = trt, kind = "sd",
                                                   col = c("#7a7aff", "#92ff92", "#ff8b8b",  "#ffd992")
                                                )
     )
with(metadata.courge2016.rare5k.D1, points(asvs.comm.courge2016.rare5k.D1.mds,
                                       display = "sites",
                                       pch = 21, 
                                       col = alpha(colvec[trt], 0.5),
                                       bg = alpha(colvec[trt], 0.5),
                                       cex=pseudo_count/2))
#date2
ordiplot(asvs.comm.courge2016.rare5k.D2.mds, display="sites", type = 'n', main = "Mid season 2016")

with(metadata.courge2016.rare5k.D2, ordispider(asvs.comm.courge2016.rare5k.D2.mds, 
                                                   group = trt, label = FALSE,
                                                   col = c("#e1e1ff", "#bfffbf", "#ffd0d0",  "#ffecc9")
                                               )
     )
with(metadata.courge2016.rare5k.D2, ordiellipse(asvs.comm.courge2016.rare5k.D2.mds, 
                                                   group = trt, kind = "sd",
                                                   col = c("#7a7aff","#92ff92", "#ff8b8b",  "#ffd992")
                                                )
     )
with(metadata.courge2016.rare5k.D2, points(asvs.comm.courge2016.rare5k.D2.mds,
                                       display = "sites",
                                       pch = 21, 
                                       col = alpha(colvec[trt], 0.5),
                                       bg = alpha(colvec[trt], 0.5),
                                       cex=pseudo_count/2))

#date3
ordiplot(asvs.comm.courge2016.rare5k.D3.mds, display="sites", type = 'n', main = "Late season 2016")

with(metadata.courge2016.rare5k.D3, ordispider(asvs.comm.courge2016.rare5k.D3.mds, 
                                                   group = trt, label = FALSE,
                                                   col = c("#e1e1ff", "#bfffbf", "#ffd0d0",  "#ffecc9")
                                               )
     )
with(metadata.courge2016.rare5k.D3, ordiellipse(asvs.comm.courge2016.rare5k.D3.mds, 
                                                   group = trt, kind = "sd",
                                                   col = c("#7a7aff", "#92ff92", "#ff8b8b",  "#ffd992")
                                                )
     )
with(metadata.courge2016.rare5k.D3, points(asvs.comm.courge2016.rare5k.D3.mds,
                                       display = "sites",
                                       pch = 21, 
                                       col = alpha(colvec[trt], 0.5),
                                       bg = alpha(colvec[trt], 0.5),
                                       cex=pseudo_count/2)
     )
plot.new()
with(legend.lvl, legend("topleft", legend = factor(freq), 
                        bty = "n", pch = 21,  
                        pt.cex=freq/2,
                        cex=2,
                        pt.bg = "black")
     )
with(metadata.courge2016.rare5k.D3, legend("bottomleft", legend = c("PC", "RCC", "CT-RCC", "BS"), 
                        bty = "n", pch = 21, 
                        pt.cex=1.5,
                        cex=1.2,
                        pt.bg = alpha(colvec, 0.5)
                        )
     )
```

###2017 bottom panel
```{r}
legend.lvl <- data.frame(freq=metadata.courge2017.rare5k[order(metadata.courge2017.rare5k$log_cfu_rif), "log_cfu_rif"])
legend.lvl.smry<-summary(legend.lvl)
legend.lvl <- data.frame(freq=c(as.numeric(strsplit(legend.lvl.smry, split = ":")[[1]][2]),
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[2]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[3]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[4]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[5]][2]), 
                                as.numeric(strsplit(legend.lvl.smry, split = ":")[[6]][2])
                                )
                         )

par(mfrow=c(1,4)) 
colvec <- c("#0000FF", #blue
            "#006400", #green
            "#800000", #red
            "#FFA500") #yellow
#date1
ordiplot(asvs.comm.courge2017.rare5k.D1.mds, display="sites", type = 'n', main = "Early season 2017")

with(metadata.courge2017.rare5k.D1, ordispider(asvs.comm.courge2017.rare5k.D1.mds, 
                                               group = trt, label = FALSE,
                                               col = c("#e1e1ff", "#bfffbf", "#ffd0d0", "#ffecc9")
                                               )
     )
with(metadata.courge2017.rare5k.D1, ordiellipse(asvs.comm.courge2017.rare5k.D1.mds, 
                                                group = trt, kind = "sd",
                                                col = c("#7a7aff", "#92ff92", "#ff8b8b", "#ffd992")
                                                )
     )
with(metadata.courge2017.rare5k.D1, points(asvs.comm.courge2017.rare5k.D1.mds,
                                           display = "sites",
                                           pch = 21, 
                                           col = alpha(colvec[trt], 0.5),
                                           bg = alpha(colvec[trt], 0.5),
                                           cex=log_cfu_rif/1.5))
#date2
ordiplot(asvs.comm.courge2017.rare5k.D2.mds, display="sites", type = 'n', main = "Mid season 2017")

with(metadata.courge2017.rare5k.D2, ordispider(asvs.comm.courge2017.rare5k.D2.mds, 
                                               group = trt, label = FALSE,
                                               col = c("#e1e1ff", "#bfffbf","#ffd0d0",  "#ffecc9")
                                               )
     )
with(metadata.courge2017.rare5k.D2, ordiellipse(asvs.comm.courge2017.rare5k.D2.mds, 
                                                group = trt, kind = "sd",
                                                col = c("#7a7aff","#92ff92", "#ff8b8b",  "#ffd992")
                                                )
     )
with(metadata.courge2017.rare5k.D2, points(asvs.comm.courge2017.rare5k.D2.mds,
                                           display = "sites",
                                           pch = 21, 
                                           col = alpha(colvec[trt], 0.5),
                                           bg = alpha(colvec[trt], 0.5),
                                           cex=log_cfu_rif/1.5))
#date3
ordiplot(asvs.comm.courge2017.rare5k.D3.mds, display="sites", type = 'n', main = "Late season 2017")

with(metadata.courge2017.rare5k.D3, ordispider(asvs.comm.courge2017.rare5k.D3.mds, 
                                               group = trt, label = FALSE,
                                               col = c("#e1e1ff", "#bfffbf","#ffd0d0",  "#ffecc9")
                                               )
     )

with(metadata.courge2017.rare5k.D3, ordiellipse(asvs.comm.courge2017.rare5k.D3.mds, 
                                                group = trt, kind = "sd",
                                                col = c("#7a7aff", "#92ff92","#ff8b8b",  "#ffd992")
                                                )
     )

with(metadata.courge2017.rare5k.D3, points(asvs.comm.courge2017.rare5k.D3.mds,
                                           display = "sites",
                                           pch = 21, 
                                           col = alpha(colvec[trt], 0.5),
                                           bg = alpha(colvec[trt], 0.5),
                                           cex=log_cfu_rif/1.5))
plot.new()
with(legend.lvl, legend("topleft", legend = factor(freq), 
                        bty = "n", pch = 21,  
                        pt.cex=freq/1.5,
                        cex=2,
                        pt.bg = "black")
     )
with(metadata.courge2017.rare5k.D3, legend("bottomleft", legend = c("PC", "RCC", "CT-RCC", "BS"), 
                                           bty = "n", pch = 21, 
                                           pt.cex=1.5,
                                           cex=1.2,
                                           pt.bg = alpha(colvec, 0.5)
                                           )
     )
```

## (line 302) MMDS scores against treatment model and Supplemental Figure 9 
### 2016
```{r extract axis scores 2016}
metadata.courge2016.rare5k.D1$NMDS1 <- vegan::scores(asvs.comm.courge2016.rare5k.D1.mds)[,1]
metadata.courge2016.rare5k.D1$NMDS2 <- vegan::scores(asvs.comm.courge2016.rare5k.D1.mds)[,2]

metadata.courge2016.rare5k.D2$NMDS1 <- vegan::scores(asvs.comm.courge2016.rare5k.D2.mds)[,1]
metadata.courge2016.rare5k.D2$NMDS2 <- vegan::scores(asvs.comm.courge2016.rare5k.D2.mds)[,2]

metadata.courge2016.rare5k.D3$NMDS1 <- vegan::scores(asvs.comm.courge2016.rare5k.D3.mds)[,1]
metadata.courge2016.rare5k.D3$NMDS2 <- vegan::scores(asvs.comm.courge2016.rare5k.D3.mds)[,2]
```
```{r modelling axis score as a function of treatment 2016}
courge2016.rare5k.D1.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2016.rare5k.D1)
courge2016.rare5k.D1.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2016.rare5k.D1)

courge2016.rare5k.D2.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2016.rare5k.D2)
courge2016.rare5k.D2.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2016.rare5k.D2)

courge2016.rare5k.D3.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2016.rare5k.D3)
courge2016.rare5k.D3.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2016.rare5k.D3)

summary(courge2016.rare5k.D1.NMDS1.lmer)
summary(courge2016.rare5k.D1.NMDS2.lmer)

summary(courge2016.rare5k.D2.NMDS1.lmer)
summary(courge2016.rare5k.D2.NMDS2.lmer)

summary(courge2016.rare5k.D3.NMDS1.lmer)
summary(courge2016.rare5k.D3.NMDS2.lmer)
```
```{r emmeans plot of MMDS scores against treatment model results 2016}
g1 <- plot(emmeans(courge2016.rare5k.D1.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 1 2016")
g2 <- plot(emmeans(courge2016.rare5k.D1.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 1 2016")

g3 <- plot(emmeans(courge2016.rare5k.D2.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 2 2016")
g4 <- plot(emmeans(courge2016.rare5k.D2.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 2 2016")

g5 <- plot(emmeans(courge2016.rare5k.D3.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 3 2016")
g6 <- plot(emmeans(courge2016.rare5k.D3.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 3 2016")

multiplot(g1, g3, g5, g2, g4, g6, cols = 3)
```

###2017
```{r extract axis scores 2017}
metadata.courge2017.rare5k.D1$NMDS1 <- vegan::scores(asvs.comm.courge2017.rare5k.D1.mds)[,1]
metadata.courge2017.rare5k.D1$NMDS2 <- vegan::scores(asvs.comm.courge2017.rare5k.D1.mds)[,2]

metadata.courge2017.rare5k.D2$NMDS1 <- vegan::scores(asvs.comm.courge2017.rare5k.D2.mds)[,1]
metadata.courge2017.rare5k.D2$NMDS2 <- vegan::scores(asvs.comm.courge2017.rare5k.D2.mds)[,2]

metadata.courge2017.rare5k.D3$NMDS1 <- vegan::scores(asvs.comm.courge2017.rare5k.D3.mds)[,1]
metadata.courge2017.rare5k.D3$NMDS2 <- vegan::scores(asvs.comm.courge2017.rare5k.D3.mds)[,2]
```
```{r modelling axis score as a function of treatment 2017}
courge2017.rare5k.D1.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2017.rare5k.D1)
courge2017.rare5k.D1.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2017.rare5k.D1)

courge2017.rare5k.D2.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2017.rare5k.D2)
courge2017.rare5k.D2.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2017.rare5k.D2)

courge2017.rare5k.D3.NMDS1.lmer <- lmerTest::lmer(NMDS1 ~ trt + (1|block), data = metadata.courge2017.rare5k.D3)
courge2017.rare5k.D3.NMDS2.lmer <- lmerTest::lmer(NMDS2 ~ trt + (1|block), data = metadata.courge2017.rare5k.D3)

summary(courge2017.rare5k.D1.NMDS1.lmer)
summary(courge2017.rare5k.D1.NMDS2.lmer)

summary(courge2017.rare5k.D2.NMDS1.lmer)
summary(courge2017.rare5k.D2.NMDS2.lmer)

summary(courge2017.rare5k.D3.NMDS1.lmer)
summary(courge2017.rare5k.D3.NMDS2.lmer)
```
```{r emmeans plot of MMDS scores against treatment model results 2017}
g71 <- plot(emmeans(courge2017.rare5k.D1.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 1 2017")
g72 <- plot(emmeans(courge2017.rare5k.D1.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 1 2017")

g73 <- plot(emmeans(courge2017.rare5k.D2.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 2 2017")
g74 <- plot(emmeans(courge2017.rare5k.D2.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 2 2017")

g75 <- plot(emmeans(courge2017.rare5k.D3.NMDS1.lmer, "trt"))+ggtitle("NMDS1 trt effect Date 3 2017")
g76 <- plot(emmeans(courge2017.rare5k.D3.NMDS2.lmer, "trt"))+ggtitle("NMDS2 trt effect Date 3 2017")

multiplot(g71, g72, g73, g74, g75, g76, cols = 3)
```


#######################################################################################################################################################
#Differentially abundant taxa among treatments at each sampling date: Sphingomonas and Methylobacterium were more abundant with cover crop treatments.#
#######################################################################################################################################################

##NA2taxa function
```{r NA2taxa}
NA2taxa<-function(dat.taxo, oneNa = F){sapply(1:length(dat.taxo$Genus), function(x){
  ifelse(test = isTRUE(oneNa),
         yes = ifelse(test = is.na(dat.taxo[x, "Genus"]), 
                       yes = ifelse(test = is.na(dat.taxo[x, "Family"]),
                                    yes=ifelse(test = is.na(dat.taxo[x, "Order"]),
                                               yes = ifelse(test = is.na(dat.taxo[x, "Class"]),
                                                            yes= paste(dat.taxo[x, "Phylum"], "NA", sep = "."),
                                                            no = paste(dat.taxo[x, "Class"],"NA", sep = ".")),
                                               no = paste(dat.taxo[x, "Order"],"NA", sep = ".")), 
                                    no=paste(dat.taxo[x, "Family"],"NA", sep = ".")),
                       no=paste(dat.taxo[x, "Genus"])),
         
         no= ifelse(test = is.na(dat.taxo[x, "Genus"]), 
                   yes = ifelse(test = is.na(dat.taxo[x, "Family"]),
                                yes=ifelse(test = is.na(dat.taxo[x, "Order"]),
                                           yes = ifelse(test = is.na(dat.taxo[x, "Class"]),
                                                        yes= paste(dat.taxo[x, "Phylum"], 
                                                                   dat.taxo[x, "Class"], 
                                                                   dat.taxo[x, "Order"], 
                                                                   dat.taxo[x, "Family"], 
                                                                   dat.taxo[x, "Genus"], sep = "."),
                                                        no = paste(dat.taxo[x, "Class"], 
                                                                   dat.taxo[x, "Order"], 
                                                                   dat.taxo[x, "Family"], 
                                                                   dat.taxo[x, "Genus"], sep = ".")),
                                           no = paste(dat.taxo[x, "Order"], 
                                                      dat.taxo[x, "Family"], 
                                                      dat.taxo[x, "Genus"], sep = ".")), 
                                no=paste(dat.taxo[x, "Family"], 
                                         dat.taxo[x, "Genus"], sep = ".")),
                   no=paste(dat.taxo[x, "Genus"])))
  })
}
```

##Tree2Heatmap function
```{r Tree2Heatmap}
tree2heatmap <-function(dataF, tree) {
  d = fortify(tree)
  d = subset(d, isTip)
  dataF <- dataF[with(d, label[order(y, decreasing = FALSE)]),]
  dataF$tip.label <- NA2taxa(dat.taxo = dataF, oneNa = T)
  t.ggtree <- ggtree(tree, branch.length = "none")
  
  dataF.mx <- as.matrix(dataF[,c(1:6)])
  
  nA<-c(1:length(dataF.mx[,1])) # Define a vector with row count
  
  #Min heat color
  col2rgb("#FFFFBF")
  rmin <- 1
  gmin <- 1
  bmin <- 191/255
  
  #Max Heat color for Postivie value
  rPmax <- 213/255
  gPmax <- 94/255
  bPmax <- 0
  
  #Max Heat color for Negative value
  col2rgb("#0072B2")
  rNmax <- 0
  gNmax <- 144/255
  bNmax <- 178/255  
  
  #viewport solution at : https://stackoverflow.com/questions/14124373/combine-base-and-ggplot-graphics-in-r-figure-window#14125565
  t.ggtree
  vp.right <- viewport(just =c("right","bottom"), height = unit(1, "npc"), width = unit(0.3, "npc"), x=0.3, y=0.015)
  #plot.new()
  #print(t.ggtree, vp=vp.right)
  
  par(bty="n",mar=c(0,0.2,0,0))
  
  plot(10,10, #Plot axis
       xlim=c(0, length(dataF.mx[1,])+0.5), #x axis as the number of column (=contrast)
       ylim=c(-3, length(dataF.mx[,1])), #y axis as the number of row (=DAA)
       xaxt="n", 
       yaxt="n", 
       xlab=NA, ylab=NA) # remove the actual output of axis
  
  par(bty="n", mar=c(0,0,0,0), 
      new=TRUE)
  print(t.ggtree, vp=vp.right)
  
  
  void<-sapply(nA,function(x){ #iterate of the DAA (=row)
    vx<-dataF.mx[x,] #Define vx as a vector of all column for the current line
    vNAx<-ifelse(is.na(vx)==T,0,vx) #Create a Vector with NA value as equal to zero
    
    #Create all polygon coordinate for the entire row
    x1 <- c(1:length(vx))-0.5 #define 1rst ppoint of the polygon
    x2 <- c(1:length(vx))+0.5 #define second point of the polygon
    y1 <- seq(x,x,length.out=length(vx))-0.5 #define 3rd point of the polygon
    y2 <- seq(x,x,length.out=length(vx))+0.5 #define 4th point of the polygon
    
    #Plot text for the current row
    par(fig=c(0.6,1,0,1), bty="n",mar=c(0,0,0,0), new=T)
    text(x=0, y=x, paste(dataF[x, "tip.label"]), cex = 0.5, pos = 4, family="Palatino",font = 1)
    
    
    par(fig=c(0.21,0.66,0,1),
        bty="n",
        #new=T,
        mar=c(0,4,0,0.2))
    
    void<-sapply(c(1:length(vNAx)), function(i){ #iterate over the contrast
      ifelse(vNAx[i]==0,
             yes = col <- "grey",
             no = ifelse(vNAx[i]>0,
                         yes = col <- rgb(red = rmin-(rmin-rPmax)*((vNAx[i]-1)/(max(dataF.mx, na.rm = T)-1)),
                                          green = gmin-(gmin-gPmax)*((vNAx[i]-1)/(max(dataF.mx, na.rm = T)-1)),
                                          blue = bmin-(bmin-bPmax)*((vNAx[i]-1)/(max(dataF.mx, na.rm = T)-1))),
                         
                         no =  col <- rgb(red = rmin-(rmin-rNmax)*((vNAx[i]+1)/(min(dataF.mx, na.rm = T)+1)),
                                          green = gmin-(gmin-gNmax)*((vNAx[i]+1)/(min(dataF.mx, na.rm = T)+1)),
                                          blue = bmin-(bmin-bNmax)*((vNAx[i]+1)/(min(dataF.mx, na.rm = T)+1)))
             )
      )
      
      polygon(x=c(x1[i], x2[i], x2[i], x1[i]), 
              y=c(y1[i], y1[i], y2[i], y2[i]),
              col= col,
              border = "grey", lwd=1)
      
    })
    
    
  })
  
  ### Legend
  cS <- c(-length(dataF.mx[1,]):-1,0:length(dataF.mx[1,]))
  
  void<-sapply(1:length(cS),function(x){
    ifelse(cS[x]>0, 
           yes = col <- rgb(red = rmin-(rmin-rPmax)*(cS[x]/max(cS)),
                            green = gmin-(gmin-gPmax)*(cS[x]/max(cS)),
                            blue = bmin-(bmin-bPmax)*(cS[x]/max(cS))),
           no =  col <- rgb(red = rmin-(rmin-rNmax)*(cS[x]/min(cS)),
                            green = gmin-(gmin-gNmax)*(cS[x]/min(cS)),
                            blue = bmin-(bmin-bNmax)*(cS[x]/min(cS)))
    )
    aDj=0.25
    Mn <- 1/2+aDj #Min scale
    Mx <- length(dataF.mx[1,])+aDj #Max Scale
    dMnMx <- Mx-Mn #distance Min to Max Scale
    
    vMn <- min(cS)  #Min value
    vMx <- max(cS)  #Max Value
    dV <- dist(c(vMn,vMx)) #Distance max: distance min to max value
    
    dMnX <- dist(c(vMn, cS[x])) #Current Distance: distance minValue to current value
    dR <- dMnX/dV #Distance ratio between current distance and max distance
    
    xCf <- Mn+dMnMx*dR #Final coordonate : applying distance-value-ratio to distance scale
    
    polygon(x=c(xCf-aDj, xCf+aDj, xCf+aDj, xCf-aDj), 
            y=c(-3,-3,-2,-2),
            col=col,
            border = "grey", lwd=1)
    })
  }
```

## 2016 DeSeq2
###Define contrast
```{r define contrast 2016}
contrast_test <- c("Seigle", "Seigle", "Seigle", "Seigle_Bru", "Seigle_Bru", "Plastique")
contrast_ref <- c("Sol_nu", "Seigle_Bru", "Plastique", "Plastique",  "Sol_nu", "Sol_nu")
contrast2016 <- data.frame(contrast_test=contrast_test, contrast_ref=contrast_ref)
```

###Deseq function
```{r Deseq function 2016}
deseqMe <- function(deseqDat, ps, test="LRT", alpha = 0.01, fold=1, contrast_test, contrast_ref){
comm.dds <- DESeq(deseqDat, test=test, reduced=~block,
                  sfType="poscounts", minmu=1e-6, minRep=Inf, parallel = T)
comm.dds.res <- results(comm.dds, contrast = c("trt", #group
                                               as.character(contrast_test), #tested
                                               as.character(contrast_ref))) #reference
comm.dds.padj <- comm.dds.res[which(comm.dds.res$padj < alpha), ]
comm.dds.padj.fold <- comm.dds.padj[which( comm.dds.padj$log2FoldChange > fold | comm.dds.padj$log2FoldChange < (-fold)),]
comm.dds.padj.fold.final <- cbind(as(comm.dds.padj.fold, "data.frame"),
                                  as(tax_table(ps)[rownames(comm.dds.padj.fold), ], "matrix"))
# Phylum order
x = tapply(comm.dds.padj.fold.final$log2FoldChange, comm.dds.padj.fold.final$Phylum, function(x) max(x))
x = sort(x, TRUE)
comm.dds.padj.fold.final$Phylum = factor(as.character(comm.dds.padj.fold.final$Phylum), levels=names(x))
# Genus order
x = tapply(comm.dds.padj.fold.final$log2FoldChange, comm.dds.padj.fold.final$Genus, function(x) max(x))
x = sort(x, TRUE)
comm.dds.padj.fold.final$Genus = factor(as.character(comm.dds.padj.fold.final$Genus), levels=names(x))
return(comm.dds.padj.fold.final)
}
```

### CodaSeq Filters
```{r filters 2016}
dim(asvs.comm.courge2016)
#filter
asvs.comm.courge2016.f<-t(codaSeq.filter(asvs.comm.courge2016, min.reads=1000, min.prop=0.00001, min.occurrence=0.005, samples.by.row=T))
dim(asvs.comm.courge2016.f)
both <- intersect(rownames(asvs.comm.courge2016.f), rownames(metadata.courge2016))
asvs.comm.courge2016.f <- asvs.comm.courge2016.f[both,]
metadata.courge2016.f <- metadata.courge2016[both,]
metadata.courge2016.f$trt <- factor(metadata.courge2016.f$trt)
metadata.courge2016.f$block <- factor(metadata.courge2016.f$block)
asvs.taxo.courge2016.f <- asvs.taxo.courge2016[colnames(asvs.comm.courge2016.f),]
```

###Zero imputation with log2(n_czm + 1-min(n_czm)) 
```{r Zero imputation 2016}
asvs.comm.courge2016.f.n0 <- as.matrix(cmultRepl(asvs.comm.courge2016.f, method="CZM", label=0, output = "p-counts"))
dim(asvs.comm.courge2016.f.n0)
which(asvs.comm.courge2016.f.n0<0)

both <- intersect(rownames(asvs.comm.courge2016.f.n0), rownames(metadata.courge2016))
asvs.comm.courge2016.f.n0 <- asvs.comm.courge2016.f.n0[both,]
metadata.courge2016.f.n0 <- metadata.courge2016[both,]
asvs.taxo.courge2016.f.n0 <- asvs.taxo.courge2016[colnames(asvs.comm.courge2016.f.n0),]
```

###adjust value
```{r adjust value 2016}
asvs.comm.courge2016.f.n0.adjusted <- asvs.comm.courge2016.f.n0+(1-min(asvs.comm.courge2016.f.n0))
```

###split by date
```{r split by date 2016}
metadata.courge2016.f.D1 <- metadata.courge2016.f[which(metadata.courge2016.f$date=="Date 1"),]
metadata.courge2016.f.D1$date <- factor(metadata.courge2016.f.D1$date)
metadata.courge2016.f.D1$block <- factor(metadata.courge2016.f.D1$block)
metadata.courge2016.f.D1$trt <- factor(metadata.courge2016.f.D1$trt)
asvs.comm.courge2016.f.n0.D1 <- asvs.comm.courge2016.f.n0.adjusted[rownames(metadata.courge2016.f.D1),]
asvs.comm.courge2016.f.n0.D1 <- asvs.comm.courge2016.f.n0.D1[,which(apply(asvs.comm.courge2016.f.n0.D1, 2, sum)>0)]
asvs.taxo.courge2016.f.n0.D1 <- asvs.taxo.courge2016.f.n0[colnames(asvs.comm.courge2016.f.n0.D1),]

metadata.courge2016.f.D2 <- metadata.courge2016.f[which(metadata.courge2016.f$date=="Date 2"),]
metadata.courge2016.f.D2$date <- factor(metadata.courge2016.f.D2$date)
metadata.courge2016.f.D2$block <- factor(metadata.courge2016.f.D2$block)
metadata.courge2016.f.D2$trt <- factor(metadata.courge2016.f.D2$trt)
asvs.comm.courge2016.f.n0.D2 <- asvs.comm.courge2016.f.n0.adjusted[rownames(metadata.courge2016.f.D2),]
asvs.comm.courge2016.f.n0.D2 <- asvs.comm.courge2016.f.n0.D2[,which(apply(asvs.comm.courge2016.f.n0.D2, 2, sum)>0)]
asvs.taxo.courge2016.f.n0.D2 <- asvs.taxo.courge2016.f.n0[colnames(asvs.comm.courge2016.f.n0.D2),]

metadata.courge2016.f.D3 <- metadata.courge2016.f[which(metadata.courge2016.f$date=="Date 3"),]
metadata.courge2016.f.D3$date <- factor(metadata.courge2016.f.D3$date)
metadata.courge2016.f.D3$block <- factor(metadata.courge2016.f.D3$block)
metadata.courge2016.f.D3$trt <- factor(metadata.courge2016.f.D3$trt)
asvs.comm.courge2016.f.n0.D3 <- asvs.comm.courge2016.f.n0.adjusted[rownames(metadata.courge2016.f.D3),]
asvs.comm.courge2016.f.n0.D3 <- asvs.comm.courge2016.f.n0.D3[,which(apply(asvs.comm.courge2016.f.n0.D3, 2, sum)>0)]
asvs.taxo.courge2016.f.n0.D3 <- asvs.taxo.courge2016.f.n0[colnames(asvs.comm.courge2016.f.n0.D3),]
```

###Phyloseq object
```{r Phyloseq object 2016}
asvs.taxo.courge2016.f.n0.D1.ps <- phyloseq(otu_table(asvs.comm.courge2016.f.n0.D1, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2016.f.n0.D1),
                                    sample_data(metadata.courge2016.f.D1))

asvs.taxo.courge2016.f.n0.D2.ps <- phyloseq(otu_table(asvs.comm.courge2016.f.n0.D2, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2016.f.n0.D2),
                                    sample_data(metadata.courge2016.f.D2))

asvs.taxo.courge2016.f.n0.D3.ps <- phyloseq(otu_table(asvs.comm.courge2016.f.n0.D3, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2016.f.n0.D3),
                                    sample_data(metadata.courge2016.f.D3))
```

###Physlodeq 2 deseq2 data
```{r Physlodeq 2 deseq2 data 2016}
asvs.taxo.courge2016.f.n0.D1.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2016.f.n0.D1.ps, design = ~ block + trt)

asvs.taxo.courge2016.f.n0.D2.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2016.f.n0.D2.ps, design = ~ block + trt)

asvs.taxo.courge2016.f.n0.D3.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2016.f.n0.D3.ps, design = ~ block + trt)
```

###Deseq2
```{r Deseq2 2016, message=FALSE}
asvs.taxo.courge2016.f.n0.D1.dds.final <- apply(contrast2016, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2016.f.n0.D1.ps.dat,
                 ps = asvs.taxo.courge2016.f.n0.D1.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})

asvs.taxo.courge2016.f.n0.D2.dds.final <- apply(contrast2016, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2016.f.n0.D2.ps.dat,
                 ps = asvs.taxo.courge2016.f.n0.D2.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})

asvs.taxo.courge2016.f.n0.D3.dds.final <- apply(contrast2016, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2016.f.n0.D3.ps.dat,
                 ps = asvs.taxo.courge2016.f.n0.D3.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})
```

## DATE 1 2016
### Extract Deseq2 results 
```{r extract Deseq2 results date 1 2016, message=FALSE}
#Date1
asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016 <- list(data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2016.f.n0.D1.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2016.f.n0.D1.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2016.f.n0.D1.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2016.f.n0.D1.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2016.f.n0.D1.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D1.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2016.f.n0.D1.dds.final[[6]][,"log2FoldChange"]))

D1.l2F.2016 <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[1]], 
                                   asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[2]], by="asv"),
                  asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[3]], by = "asv"),
              asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[4]], by = "asv"),
            asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[5]], by = "asv"),
          asvs.taxo.courge2016.f.n0.D1.dds.final.l2F.2016[[6]], by = "asv")

D1.l2F.2016$mean <- apply(asvs.taxo.courge2016.f.n0.D1.ps@otu_table@.Data[,D1.l2F.2016$asv],2, mean) # Extract mean from original comm matrix
D1.l2F.2016$abund <- apply(asvs.taxo.courge2016.f.n0.D1.ps@otu_table@.Data[,D1.l2F.2016$asv],2, sum) # Extract mean from original comm matrix
D1.l2F.2016 <- D1.l2F.2016[order(D1.l2F.2016$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D1.l2F.2016.final <- cbind(as(D1.l2F.2016, "data.frame"), as(tax_table(asvs.taxo.courge2016.f.n0.D1.ps)[D1.l2F.2016$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D1.l2F.2016.final <- D1.l2F.2016.final[which(!D1.l2F.2016.final$Class=="Chloroplast"),] # Remove chloroplast entry

#D1.l2F.2016.final <- D1.l2F.2016.final[which(!D1.l2F.2016.final$asv=="CGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAGCCATGCAGCACCTGTGTCCACGTCCCGAAGGAAGGAAACCGTCTCCGGTAACCGTCGTGGCATGTCAAGAGCTGGTAAGGTTCTGCGCGTT#GCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTTAATCTTGCGACCGTACTCCCCAGGCGGGAAGCTTAATGCGTTGGCTGCGCCACCGAGTGATAAATCACCCGACGGCTAGCTTCCATCGTTTACGGCGTGGACTAC"),] #remove weird #sequence

D1.l2F.2016.final.T100 <- D1.l2F.2016.final[1:100,] #Extract Top 100 DAA

D1.l2F.2016.final.T100 <- D1.l2F.2016.final.T100[order(D1.l2F.2016.final.T100$Genus),]#Order on Genus Name
D1.l2F.2016.final.T100 <- D1.l2F.2016.final.T100[order(D1.l2F.2016.final.T100$Family),]#Order on Genus Name
D1.l2F.2016.final.T100 <- D1.l2F.2016.final.T100[order(D1.l2F.2016.final.T100$Class),]#Order on Genus Name
D1.l2F.2016.final.T100 <- D1.l2F.2016.final.T100[order(D1.l2F.2016.final.T100$Order),]#Order on Genus Name
D1.l2F.2016.final.T100 <- D1.l2F.2016.final.T100[order(D1.l2F.2016.final.T100$Phylum),]#Order on Genus Name

D1.l2F.2016.final.T100 <- data.frame(D1.l2F.2016.final.T100[,2:length(names(D1.l2F.2016.final.T100))], row.names = D1.l2F.2016.final.T100$asv) #put asv as rownames

```

### Phylogenetic tree Date 1
```{r pynast to fastTree Date 1 2016}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D1.l2F.2016.final.T100), file = "./D1.l2F.2016.final.T100.fna", repName = FALSE)

#Import pynast + fastTree results
D1.l2F.2016.final.T100.tree <- read_tree("./D1.l2F.2016.final.T100_aligned_pfiltered.tre")
```

## DATE 2 2016 
### Extract Deseq2 results 
```{r dds results to log2foldchange date 2 2016}
#Date2
asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016 <- list(data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2016.f.n0.D2.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2016.f.n0.D2.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2016.f.n0.D2.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2016.f.n0.D2.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2016.f.n0.D2.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D2.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2016.f.n0.D2.dds.final[[6]][,"log2FoldChange"]))
D2.l2F.2016 <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[1]], 
                                   asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[2]], by="asv"),
                  asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[3]], by = "asv"),
              asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[4]], by = "asv"),
            asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[5]], by = "asv"),
          asvs.taxo.courge2016.f.n0.D2.dds.final.l2F.2016[[6]], by = "asv")


D2.l2F.2016$mean <- apply(asvs.taxo.courge2016.f.n0.D2.ps@otu_table@.Data[,D2.l2F.2016$asv],2, mean) # Extract mean from original comm matrix
D2.l2F.2016$abund <- apply(asvs.taxo.courge2016.f.n0.D2.ps@otu_table@.Data[,D2.l2F.2016$asv],2, sum) # Extract mean from original comm matrix
D2.l2F.2016 <- D2.l2F.2016[order(D2.l2F.2016$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D2.l2F.2016.final <- cbind(as(D2.l2F.2016, "data.frame"), as(tax_table(asvs.taxo.courge2016.f.n0.D2.ps)[D2.l2F.2016$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D2.l2F.2016.final <- D2.l2F.2016.final[which(!D2.l2F.2016.final$Class=="Chloroplast"),] # Remove chloroplast entry

D2.l2F.2016.final.T100 <- D2.l2F.2016.final[1:100,] #Extract Top 100 DAA

D2.l2F.2016.final.T100 <- D2.l2F.2016.final.T100[order(D2.l2F.2016.final.T100$Genus),]#Order on Genus Name
D2.l2F.2016.final.T100 <- D2.l2F.2016.final.T100[order(D2.l2F.2016.final.T100$Family),]#Order on Genus Name
D2.l2F.2016.final.T100 <- D2.l2F.2016.final.T100[order(D2.l2F.2016.final.T100$Class),]#Order on Genus Name
D2.l2F.2016.final.T100 <- D2.l2F.2016.final.T100[order(D2.l2F.2016.final.T100$Order),]#Order on Genus Name
D2.l2F.2016.final.T100 <- D2.l2F.2016.final.T100[order(D2.l2F.2016.final.T100$Phylum),]#Order on Genus Name

D2.l2F.2016.final.T100 <- data.frame(D2.l2F.2016.final.T100[,2:length(names(D2.l2F.2016.final.T100))], row.names = D2.l2F.2016.final.T100$asv) #put asv as rownames
```

### Phylogenetic tree Date 2
```{r pynast to fastTree Date 2 2016}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D2.l2F.2016.final.T100), file = "._article1/D2.l2F.2016.final.T100.fna", repName = FALSE)

#Import pynast + fastTree results
D2.l2F.2016.final.T100.tree <- read_tree("./D2.l2F.2016.final_aligned_pfiltered_mega.tre")
```

##Date 3
### Extract Deseq2 results
```{r extract Deseq2 results date 3 2016, message=FALSE}
#Date3
asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016 <- list(data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2016.f.n0.D3.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2016.f.n0.D3.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2016.f.n0.D3.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2016.f.n0.D3.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2016.f.n0.D3.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2016.f.n0.D3.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2016.f.n0.D3.dds.final[[6]][,"log2FoldChange"]))
D3.l2F.2016 <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[1]], 
                                   asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[2]], by="asv"),
                  asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[3]], by = "asv"),
              asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[4]], by = "asv"),
            asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[5]], by = "asv"),
          asvs.taxo.courge2016.f.n0.D3.dds.final.l2F.2016[[6]], by = "asv")

D3.l2F.2016$mean <- apply(asvs.taxo.courge2016.f.n0.D3.ps@otu_table@.Data[,D3.l2F.2016$asv],2, mean) # Extract mean from original comm matrix
D3.l2F.2016$abund <- apply(asvs.taxo.courge2016.f.n0.D3.ps@otu_table@.Data[,D3.l2F.2016$asv],2, sum) # Extract mean from original comm matrix
D3.l2F.2016 <- D3.l2F.2016[order(D3.l2F.2016$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D3.l2F.2016.final <- cbind(as(D3.l2F.2016, "data.frame"), as(tax_table(asvs.taxo.courge2016.f.n0.D3.ps)[D3.l2F.2016$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D3.l2F.2016.final <- D3.l2F.2016.final[which(!D3.l2F.2016.final$Class=="Chloroplast"),] # Remove chloroplast entry

D3.l2F.2016.final.T100 <- D3.l2F.2016.final[1:38,] #Extract Top 100 DAA

D3.l2F.2016.final.T100 <- D3.l2F.2016.final.T100[order(D3.l2F.2016.final.T100$Genus),]#Order on Genus Name
D3.l2F.2016.final.T100 <- D3.l2F.2016.final.T100[order(D3.l2F.2016.final.T100$Family),]#Order on Genus Name
D3.l2F.2016.final.T100 <- D3.l2F.2016.final.T100[order(D3.l2F.2016.final.T100$Class),]#Order on Genus Name
D3.l2F.2016.final.T100 <- D3.l2F.2016.final.T100[order(D3.l2F.2016.final.T100$Order),]#Order on Genus Name
D3.l2F.2016.final.T100 <- D3.l2F.2016.final.T100[order(D3.l2F.2016.final.T100$Phylum),]#Order on Genus Name

D3.l2F.2016.final.T100 <- data.frame(D3.l2F.2016.final.T100[,2:length(names(D3.l2F.2016.final.T100))], row.names = D3.l2F.2016.final.T100$asv) #put asv as rownames
```

### Phylogenetic tree Date 3
```{r pynast to fastTree Date 3 2016}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D3.l2F.2016.final.T100), file = "._article1/D3.l2F.2016.final.T100.fna", repName = FALSE)

#Import pynast + fastTree results
D3.l2F.2016.final.T100.tree <- read_tree("./D3.l2F.2016.final.T100_aligned.tre")
```

## (line 323) Supplemental Figure 10
### left panel
```{r tree2heatmap date 1 2016}
tree2heatmap(dataF = D1.l2F.2016.final.T100, tree = D1.l2F.2016.final.T100.tree)
```

### middle panel
```{r tree2heatmap date 2 2016}
tree2heatmap(dataF = D2.l2F.2016.final.T100, tree = D2.l2F.2016.final.T100.tree)
```

### right panel
```{r tree2heatmap date 3 2016}
tree2heatmap(dataF = D3.l2F.2016.final.T100, tree = D3.l2F.2016.final.T100.tree)
```

##2017
###Define contrast
```{r define contrast 2017}
contrast_test <- c("Seigle", "Seigle", "Seigle", "Seigle_bru", "Seigle_bru", "Plastique")
contrast_ref <- c("Sol_nu", "Seigle_bru", "Plastique", "Plastique",  "Sol_nu", "Sol_nu")
contrast2017 <- data.frame(contrast_test=contrast_test, contrast_ref=contrast_ref)
```

###Deseq function
```{r Deseq function 2017}
deseqMe <- function(deseqDat, ps, test="LRT", alpha = 0.01, fold=1, contrast_test, contrast_ref){
comm.dds <- DESeq(deseqDat, test=test, reduced=~block,
                  sfType="poscounts", minmu=1e-6, minRep=Inf, parallel = T)
comm.dds.res <- results(comm.dds, contrast = c("trt", #group
                                               contrast_test, #tested
                                               contrast_ref)) #reference
comm.dds.padj <- comm.dds.res[which(comm.dds.res$padj < alpha), ]
comm.dds.padj.fold <- comm.dds.padj[which( comm.dds.padj$log2FoldChange > fold | comm.dds.padj$log2FoldChange < (-fold)),]
comm.dds.padj.fold.final <- cbind(as(comm.dds.padj.fold, "data.frame"),
                                  as(tax_table(ps)[rownames(comm.dds.padj.fold), ], "matrix"))
# Phylum order
x = tapply(comm.dds.padj.fold.final$log2FoldChange, comm.dds.padj.fold.final$Phylum, function(x) max(x))
x = sort(x, TRUE)
comm.dds.padj.fold.final$Phylum = factor(as.character(comm.dds.padj.fold.final$Phylum), levels=names(x))
# Genus order
x = tapply(comm.dds.padj.fold.final$log2FoldChange, comm.dds.padj.fold.final$Genus, function(x) max(x))
x = sort(x, TRUE)
comm.dds.padj.fold.final$Genus = factor(as.character(comm.dds.padj.fold.final$Genus), levels=names(x))
return(comm.dds.padj.fold.final)
}
```

### CodaSeq Filters
```{r filters 2017}
dim(asvs.comm.courge2017)
#filter
asvs.comm.courge2017.f<-t(codaSeq.filter(asvs.comm.courge2017, min.reads=1000, min.prop=0.00001, min.occurrence=0.005, samples.by.row=T))
dim(asvs.comm.courge2017.f)
both <- intersect(rownames(asvs.comm.courge2017.f), rownames(metadata.courge2017))
asvs.comm.courge2017.f <- asvs.comm.courge2017.f[both,]
metadata.courge2017.f <- metadata.courge2017[both,]
metadata.courge2017.f$trt <- factor(metadata.courge2017.f$trt)
metadata.courge2017.f$block <- factor(metadata.courge2017.f$block)
asvs.taxo.courge2017.f <- asvs.taxo.courge2017[colnames(asvs.comm.courge2017.f),]
```

###Zero imputation with log2(n_czm + 1-min(n_czm)) 
```{r Zero imputation 2017}
asvs.comm.courge2017.f.n0 <- as.matrix(cmultRepl(asvs.comm.courge2017.f, method="CZM", label=0, output = "p-counts"))
dim(asvs.comm.courge2017.f.n0)
which(asvs.comm.courge2017.f.n0<0)

both <- intersect(rownames(asvs.comm.courge2017.f.n0), rownames(metadata.courge2017))
asvs.comm.courge2017.f.n0 <- asvs.comm.courge2017.f.n0[both,]
metadata.courge2017.f.n0 <- metadata.courge2017[both,]
asvs.taxo.courge2017.f.n0 <- asvs.taxo.courge2017[colnames(asvs.comm.courge2017.f.n0),]
```

###adjust value
```{r adjust value 2017}
asvs.comm.courge2017.f.n0.adjusted <- asvs.comm.courge2017.f.n0+(1-min(asvs.comm.courge2017.f.n0))
```

###split by date
```{r split by date 2017}
metadata.courge2017.f.D1 <- metadata.courge2017.f[which(metadata.courge2017.f$date=="Date 1"),]
metadata.courge2017.f.D1$date <- factor(metadata.courge2017.f.D1$date)
metadata.courge2017.f.D1$block <- factor(metadata.courge2017.f.D1$block)
metadata.courge2017.f.D1$trt <- factor(metadata.courge2017.f.D1$trt)
asvs.comm.courge2017.f.n0.D1 <- asvs.comm.courge2017.f.n0.adjusted[rownames(metadata.courge2017.f.D1),]
asvs.comm.courge2017.f.n0.D1 <- asvs.comm.courge2017.f.n0.D1[,which(apply(asvs.comm.courge2017.f.n0.D1, 2, sum)>0)]
asvs.taxo.courge2017.f.n0.D1 <- asvs.taxo.courge2017.f.n0[colnames(asvs.comm.courge2017.f.n0.D1),]

metadata.courge2017.f.D2 <- metadata.courge2017.f[which(metadata.courge2017.f$date=="Date 2"),]
metadata.courge2017.f.D2$date <- factor(metadata.courge2017.f.D2$date)
metadata.courge2017.f.D2$block <- factor(metadata.courge2017.f.D2$block)
metadata.courge2017.f.D2$trt <- factor(metadata.courge2017.f.D2$trt)
asvs.comm.courge2017.f.n0.D2 <- asvs.comm.courge2017.f.n0.adjusted[rownames(metadata.courge2017.f.D2),]
asvs.comm.courge2017.f.n0.D2 <- asvs.comm.courge2017.f.n0.D2[,which(apply(asvs.comm.courge2017.f.n0.D2, 2, sum)>0)]
asvs.taxo.courge2017.f.n0.D2 <- asvs.taxo.courge2017.f.n0[colnames(asvs.comm.courge2017.f.n0.D2),]

metadata.courge2017.f.D3 <- metadata.courge2017.f[which(metadata.courge2017.f$date=="Date 3"),]
metadata.courge2017.f.D3$date <- factor(metadata.courge2017.f.D3$date)
metadata.courge2017.f.D3$block <- factor(metadata.courge2017.f.D3$block)
metadata.courge2017.f.D3$trt <- factor(metadata.courge2017.f.D3$trt)
asvs.comm.courge2017.f.n0.D3 <- asvs.comm.courge2017.f.n0.adjusted[rownames(metadata.courge2017.f.D3),]
asvs.comm.courge2017.f.n0.D3 <- asvs.comm.courge2017.f.n0.D3[,which(apply(asvs.comm.courge2017.f.n0.D3, 2, sum)>0)]
asvs.taxo.courge2017.f.n0.D3 <- asvs.taxo.courge2017.f.n0[colnames(asvs.comm.courge2017.f.n0.D3),]
```

###Phyloseq object
```{r Phyloseq object 2017}
asvs.taxo.courge2017.f.n0.D1.ps <- phyloseq(otu_table(asvs.comm.courge2017.f.n0.D1, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2017.f.n0.D1),
                                    sample_data(metadata.courge2017.f.D1))

asvs.taxo.courge2017.f.n0.D2.ps <- phyloseq(otu_table(asvs.comm.courge2017.f.n0.D2, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2017.f.n0.D2),
                                    sample_data(metadata.courge2017.f.D2))

asvs.taxo.courge2017.f.n0.D3.ps <- phyloseq(otu_table(asvs.comm.courge2017.f.n0.D3, taxa_are_rows = F), 
                                    tax_table(asvs.taxo.courge2017.f.n0.D3),
                                    sample_data(metadata.courge2017.f.D3))
```

###Physlodeq 2 deseq2 data
```{r Physlodeq 2 deseq2 2017}
asvs.taxo.courge2017.f.n0.D1.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2017.f.n0.D1.ps, design = ~ block + trt)

asvs.taxo.courge2017.f.n0.D2.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2017.f.n0.D2.ps, design = ~ block + trt)

asvs.taxo.courge2017.f.n0.D3.ps.dat <- phyloseq_to_deseq2(asvs.taxo.courge2017.f.n0.D3.ps, design = ~ block + trt)
```

###Deseq2
```{r Deseq2 2017, message=FALSE}
asvs.taxo.courge2017.f.n0.D1.dds.final <- apply(contrast2017, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2017.f.n0.D1.ps.dat,
                 ps = asvs.taxo.courge2017.f.n0.D1.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})

asvs.taxo.courge2017.f.n0.D2.dds.final <- apply(contrast2017, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2017.f.n0.D2.ps.dat,
                 ps = asvs.taxo.courge2017.f.n0.D2.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})

asvs.taxo.courge2017.f.n0.D3.dds.final <- apply(contrast2017, 1, function(x){ 
  return(deseqMe(deseqDat = asvs.taxo.courge2017.f.n0.D3.ps.dat,
                 ps = asvs.taxo.courge2017.f.n0.D3.ps,
                 contrast_test = x[1],
                 contrast_ref = x[2]))
})
```

## Date 1 2017
###Extract Deseq2 results
```{r extract Deseq2 results date 1 2017, message=FALSE}
#Date1
asvs.taxo.courge2017.f.n0.D1.dds.final.l2F <- list(data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2017.f.n0.D1.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2017.f.n0.D1.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2017.f.n0.D1.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2017.f.n0.D1.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2017.f.n0.D1.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D1.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2017.f.n0.D1.dds.final[[6]][,"log2FoldChange"]))

D1.l2F <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[1]], 
                                   asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[2]], by="asv"),
                  asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[3]], by = "asv"),
              asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[4]], by = "asv"),
            asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[5]], by = "asv"),
          asvs.taxo.courge2017.f.n0.D1.dds.final.l2F[[6]], by = "asv")

D1.l2F$mean <- apply(asvs.taxo.courge2017.f.n0.D1.ps@otu_table@.Data[,D1.l2F$asv],2, mean) # Extract mean from original comm matrix
D1.l2F$abund <- apply(asvs.taxo.courge2017.f.n0.D1.ps@otu_table@.Data[,D1.l2F$asv],2, sum) # Extract mean from original comm matrix
D1.l2F <- D1.l2F[order(D1.l2F$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D1.l2F.final <- cbind(as(D1.l2F, "data.frame"), as(tax_table(asvs.taxo.courge2017.f.n0.D1.ps)[D1.l2F$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D1.l2F.final <- D1.l2F.final[which(!D1.l2F.final$Class=="Chloroplast"),] # Remove chloroplast entry

D1.l2F.final <- D1.l2F.final[which(!D1.l2F.final$asv=="CGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAGCCATGCAGCACCTGTGTCCACGTCCCGAAGGAAGGAAACCGTCTCCGGTAACCGTCGTGGCATGTCAAGAGCTGGTAAGGTTCTGCGCGTTGCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTTAATCTTGCGACCGTACTCCCCAGGCGGGAAGCTTAATGCGTTGGCTGCGCCACCGAGTGATAAATCACCCGACGGCTAGCTTCCATCGTTTACGGCGTGGACTAC"),] #remove weird sequence

D1.l2F.final.T100 <- D1.l2F.final[1:100,] #Extract Top 100 DAA

D1.l2F.final.T100 <- D1.l2F.final.T100[order(D1.l2F.final.T100$Genus),]#Order on Genus Name
D1.l2F.final.T100 <- D1.l2F.final.T100[order(D1.l2F.final.T100$Family),]#Order on Genus Name
D1.l2F.final.T100 <- D1.l2F.final.T100[order(D1.l2F.final.T100$Class),]#Order on Genus Name
D1.l2F.final.T100 <- D1.l2F.final.T100[order(D1.l2F.final.T100$Order),]#Order on Genus Name
D1.l2F.final.T100 <- D1.l2F.final.T100[order(D1.l2F.final.T100$Phylum),]#Order on Genus Name

D1.l2F.final.T100 <- data.frame(D1.l2F.final.T100[,2:length(names(D1.l2F.final.T100))], row.names = D1.l2F.final.T100$asv) #put asv as rownames

```

### Phylogenetic tree Date 1
```{r pynast to fastTree Date 1}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D1.l2F.final.T100), file = "D1.l2F.final.T100.fna", repName = FALSE)

#import pynast + FastTree result
D1.l2F.final.T100.tree <- read_tree("./D1.l2F.final.T100_aligned_pfiltered.tre")
```

## Date 2 2017
###Extract Deseq2 results
```{r extract Deseq2 results date 2 2017}
#Date2
asvs.taxo.courge2017.f.n0.D2.dds.final.l2F <- list(data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2017.f.n0.D2.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2017.f.n0.D2.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2017.f.n0.D2.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2017.f.n0.D2.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2017.f.n0.D2.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D2.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2017.f.n0.D2.dds.final[[6]][,"log2FoldChange"]))
D2.l2F <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[1]], 
                                   asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[2]], by="asv"),
                  asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[3]], by = "asv"),
              asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[4]], by = "asv"),
            asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[5]], by = "asv"),
          asvs.taxo.courge2017.f.n0.D2.dds.final.l2F[[6]], by = "asv")


D2.l2F$mean <- apply(asvs.taxo.courge2017.f.n0.D2.ps@otu_table@.Data[,D2.l2F$asv],2, mean) # Extract mean from original comm matrix
D2.l2F$abund <- apply(asvs.taxo.courge2017.f.n0.D2.ps@otu_table@.Data[,D2.l2F$asv],2, sum) # Extract mean from original comm matrix
D2.l2F <- D2.l2F[order(D2.l2F$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D2.l2F.final <- cbind(as(D2.l2F, "data.frame"), as(tax_table(asvs.taxo.courge2017.f.n0.D2.ps)[D2.l2F$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D2.l2F.final <- D2.l2F.final[which(!D2.l2F.final$Class=="Chloroplast"),] # Remove chloroplast entry

D2.l2F.final.T100 <- D2.l2F.final[1:100,] #Extract Top 100 DAA

D2.l2F.final.T100 <- D2.l2F.final.T100[order(D2.l2F.final.T100$Genus),]#Order on Genus Name
D2.l2F.final.T100 <- D2.l2F.final.T100[order(D2.l2F.final.T100$Family),]#Order on Genus Name
D2.l2F.final.T100 <- D2.l2F.final.T100[order(D2.l2F.final.T100$Class),]#Order on Genus Name
D2.l2F.final.T100 <- D2.l2F.final.T100[order(D2.l2F.final.T100$Order),]#Order on Genus Name
D2.l2F.final.T100 <- D2.l2F.final.T100[order(D2.l2F.final.T100$Phylum),]#Order on Genus Name

D2.l2F.final.T100 <- data.frame(D2.l2F.final.T100[,2:length(names(D2.l2F.final.T100))], row.names = D2.l2F.final.T100$asv) #put asv as rownames

```

### Phylogenetic tree Date 2
```{r pynast to fastTree Date 2}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D2.l2F.final.T100), file = "D2.l2F.final.T100.fna", repName = FALSE)

#import pynast + FastTree result
D2.l2F.final.T100.tree <- read_tree("./D2.l2F.final_aligned_pfiltered_mega.tre")
```

## Date 3 2017
###Extract Deseq2 results
```{r extract Deseq2 results date 3 2017}
#Date3
asvs.taxo.courge2017.f.n0.D3.dds.final.l2F <- list(data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[1]]),
                                                              RyeVsSoil=asvs.taxo.courge2017.f.n0.D3.dds.final[[1]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[2]]),
                                                              RyeVsCTR=asvs.taxo.courge2017.f.n0.D3.dds.final[[2]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[3]]),
                                                              RyeVsPlas=asvs.taxo.courge2017.f.n0.D3.dds.final[[3]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[4]]),
                                                              CTRVsPlas=asvs.taxo.courge2017.f.n0.D3.dds.final[[4]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[5]]),
                                                              CTRVsSoil=asvs.taxo.courge2017.f.n0.D3.dds.final[[5]][,"log2FoldChange"]),
                                                   data.frame(asv=rownames(asvs.taxo.courge2017.f.n0.D3.dds.final[[6]]),
                                                              PlasVsSoil=asvs.taxo.courge2017.f.n0.D3.dds.final[[6]][,"log2FoldChange"]))
D3.l2F <-  dplyr::full_join(
            dplyr::full_join(
              dplyr::full_join(
                dplyr::full_join(
                  dplyr::full_join(asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[1]], 
                                   asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[2]], by="asv"),
                  asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[3]], by = "asv"),
              asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[4]], by = "asv"),
            asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[5]], by = "asv"),
          asvs.taxo.courge2017.f.n0.D3.dds.final.l2F[[6]], by = "asv")

D3.l2F$mean <- apply(asvs.taxo.courge2017.f.n0.D3.ps@otu_table@.Data[,D3.l2F$asv],2, mean) # Extract mean from original comm matrix
D3.l2F$abund <- apply(asvs.taxo.courge2017.f.n0.D3.ps@otu_table@.Data[,D3.l2F$asv],2, sum) # Extract mean from original comm matrix
D3.l2F <- D3.l2F[order(D3.l2F$mean, decreasing = T, na.last = T),] #order DAA matrix based on original mean

D3.l2F.final <- cbind(as(D3.l2F, "data.frame"), as(tax_table(asvs.taxo.courge2017.f.n0.D3.ps)[D3.l2F$asv, ], "matrix")) # Bind taxonomic inofrmation to the DAA matrix

D3.l2F.final <- D3.l2F.final[which(!D3.l2F.final$Class=="Chloroplast"),] # Remove chloroplast entry

D3.l2F.final.T100 <- D3.l2F.final[1:100,] #Extract Top 100 DAA

D3.l2F.final.T100 <- D3.l2F.final.T100[order(D3.l2F.final.T100$Genus),]#Order on Genus Name
D3.l2F.final.T100 <- D3.l2F.final.T100[order(D3.l2F.final.T100$Family),]#Order on Genus Name
D3.l2F.final.T100 <- D3.l2F.final.T100[order(D3.l2F.final.T100$Class),]#Order on Genus Name
D3.l2F.final.T100 <- D3.l2F.final.T100[order(D3.l2F.final.T100$Order),]#Order on Genus Name
D3.l2F.final.T100 <- D3.l2F.final.T100[order(D3.l2F.final.T100$Phylum),]#Order on Genus Name

D3.l2F.final.T100 <- data.frame(D3.l2F.final.T100[,2:length(names(D3.l2F.final.T100))], row.names = D3.l2F.final.T100$asv) #put asv as rownames
```

### Phylogenetic tree Date 3
```{r pynast to fastTree Date 3}
#For pynast + FastTree
void <- ASV2fasta(asv_table = t(D3.l2F.final.T100), file = "D3.l2F.final.T100.fna", repName = FALSE)

#import pynast + FastTree result
D3.l2F.final.T100.tree <- read_tree("./D3.l2F.final_aligned_pfiltered_mega.tre")
```

## (line 323) Figure 4
### left panel
```{r tree2heatmap date 1 2017}
tree2heatmap(dataF = D1.l2F.final.T100, tree = D1.l2F.final.T100.tree)
```

### middle panel
```{r tree2heatmap date 2 2017}
tree2heatmap(dataF = D2.l2F.final.T100, tree = D2.l2F.final.T100.tree)
```

### right panel
```{r tree2heatmap date 3 2017}
tree2heatmap(dataF = D3.l2F.final.T100, tree = D3.l2F.final.T100.tree)
```

#########
#METHODS#
#########
## Effect of rarefaction over analysis and results
###2016
```{r}
# Defining Dataframe
rare.df2016.final <- NULL
#dimVSdate.df2016.final <- NULL
#asvVSdate.df2016.final <- NULL


###Initiate rarefy loop
for (i in seq(1000,10000,1000)) {
  #Parallelize
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #to not overload your computer
  registerDoParallel(cl)
  
  ###################################
  # Defining DESIRED Dataframe
  rare.df2016 <- NULL   #Permanova 
  #dimVSdate.df2016 <- NULL  #Sample Stat
  #asvVSdate.df2016 <- NULL #ASV number
  ###################################
  
  ###################################
  ###Keep only samples with the current loop minimum of sequences (i)
  asvs.comm.courge2016.temp <- asvs.comm.courge2016[which(apply(asvs.comm.courge2016,1,sum)>=i),]
  asvs.comm.courge2016.temp <- asvs.comm.courge2016.temp[,apply(asvs.comm.courge2016.temp,2,sum)>1]
  #asvs.taxo.courge2016.temp <- asvs.taxo.courge2016[colnames(asvs.comm.courge2016.temp),]
  metadata.courge2016.temp <- metadata.courge2016[rownames(asvs.comm.courge2016.temp),]
  ###################################
  
  ###################################
  #### CHANGE THE LOOP OUTPUT FOR THE DESIRED Dataframe (see above)
  ###Rarefy with 100 iteration
  rare.df2016 <- foreach(j=1:100, .combine=rbind, .packages = "vegan") %dopar% { 
    asvs.comm.courge2016.rare.temp <- rrarefy(asvs.comm.courge2016.temp, sample=i)
    asvs.comm.courge2016.rare.temp <- asvs.comm.courge2016.rare.temp[,apply(asvs.comm.courge2016.rare.temp,2,sum)>1]
    metadata.courge2016.rare.temp <- metadata.courge2016.temp[rownames(asvs.comm.courge2016.rare.temp),]
    #asvs.taxo.courge2016.rare.temp <- asvs.taxo.courge2016.temp[colnames(asvs.comm.courge2016.rare.temp),]
    ###Adjust metadata
    metadata.courge2016.rare.temp$date <- factor(metadata.courge2016.rare.temp$date)
    metadata.courge2016.rare.temp$trt <- factor(metadata.courge2016.rare.temp$trt)
    metadata.courge2016.rare.temp$type <- factor(metadata.courge2016.rare.temp$type)
    metadata.courge2016.rare.temp$block <- factor(metadata.courge2016.rare.temp$block)
    metadata.courge2016.rare.temp$sub_ech <- factor(metadata.courge2016.rare.temp$sub_ech)
    #split by date
    ###Date1
    metadata.courge2016.rare.temp.D1 <- metadata.courge2016.rare.temp[metadata.courge2016.rare.temp$date == 'Date 1',]
    metadata.courge2016.rare.temp.D1$date <- factor(metadata.courge2016.rare.temp.D1$date)
    metadata.courge2016.rare.temp.D1$trt <- factor(metadata.courge2016.rare.temp.D1$trt)
    metadata.courge2016.rare.temp.D1$type <- factor(metadata.courge2016.rare.temp.D1$type)
    metadata.courge2016.rare.temp.D1$block <- factor(metadata.courge2016.rare.temp.D1$block)
    metadata.courge2016.rare.temp.D1$sub_ech <- factor(metadata.courge2016.rare.temp.D1$sub_ech)
    metadata.courge2016.rare.temp.D1$run <- factor(metadata.courge2016.rare.temp.D1$run)
    asvs.comm.courge2016.rare.temp.D1 <- asvs.comm.courge2016.rare.temp[rownames(metadata.courge2016.rare.temp.D1),]
    asvs.comm.courge2016.rare.temp.D1 <- asvs.comm.courge2016.rare.temp.D1[,apply(asvs.comm.courge2016.rare.temp.D1,2,sum)>0]
    ###Date2
    metadata.courge2016.rare.temp.D2 <- metadata.courge2016.rare.temp[metadata.courge2016.rare.temp$date == 'Date 2',]
    metadata.courge2016.rare.temp.D2$date <- factor(metadata.courge2016.rare.temp.D2$date)
    metadata.courge2016.rare.temp.D2$trt <- factor(metadata.courge2016.rare.temp.D2$trt)
    metadata.courge2016.rare.temp.D2$type <- factor(metadata.courge2016.rare.temp.D2$type)
    metadata.courge2016.rare.temp.D2$block <- factor(metadata.courge2016.rare.temp.D2$block)
    metadata.courge2016.rare.temp.D2$sub_ech <- factor(metadata.courge2016.rare.temp.D2$sub_ech)
    metadata.courge2016.rare.temp.D2$run <- factor(metadata.courge2016.rare.temp.D2$run)
    asvs.comm.courge2016.rare.temp.D2 <- asvs.comm.courge2016.rare.temp[rownames(metadata.courge2016.rare.temp.D2),]
    asvs.comm.courge2016.rare.temp.D2 <- asvs.comm.courge2016.rare.temp.D2[,apply(asvs.comm.courge2016.rare.temp.D2,2,sum)>0]
    ###Date3
    metadata.courge2016.rare.temp.D3 <- metadata.courge2016.rare.temp[metadata.courge2016.rare.temp$date == 'Date 3',]
    metadata.courge2016.rare.temp.D3$date <- factor(metadata.courge2016.rare.temp.D3$date)
    metadata.courge2016.rare.temp.D3$trt <- factor(metadata.courge2016.rare.temp.D3$trt)
    metadata.courge2016.rare.temp.D3$type <- factor(metadata.courge2016.rare.temp.D3$type)
    metadata.courge2016.rare.temp.D3$block <- factor(metadata.courge2016.rare.temp.D3$block)
    metadata.courge2016.rare.temp.D3$sub_ech <- factor(metadata.courge2016.rare.temp.D3$sub_ech)
    metadata.courge2016.rare.temp.D3$run <- factor(metadata.courge2016.rare.temp.D3$run)
    asvs.comm.courge2016.rare.temp.D3 <- asvs.comm.courge2016.rare.temp[rownames(metadata.courge2016.rare.temp.D3),]
    asvs.comm.courge2016.rare.temp.D3 <- asvs.comm.courge2016.rare.temp.D3[,apply(asvs.comm.courge2016.rare.temp.D3,2,sum)>0]
    ####################################################
    
    #####################################################
    #Dim : Sample Stat
    #data.frame("Sample.nb.D1"= dim(asvs.comm.courge2016.rare.temp.D1)[1], 
    #           "ASVs.nb.D1"= dim(asvs.comm.courge2016.rare.temp.D1)[2],
    #           "Sample.nb.D2"= dim(asvs.comm.courge2016.rare.temp.D2)[1], 
    #           "ASVs.nb.D2"= dim(asvs.comm.courge2016.rare.temp.D2)[2],
    #           "Sample.nb.D3"= dim(asvs.comm.courge2016.rare.temp.D3)[1], 
    #           "ASVs.nb.D3"= dim(asvs.comm.courge2016.rare.temp.D3)[2])
    #####################################################
    
    #####################################################
    #Permanova
    permRare.tempD1.obj <-adonis(decostand(asvs.comm.courge2016.rare.temp.D1, method="hellinger") ~ trt, 
                                 strata=metadata.courge2016.rare.temp.D1$block, data=metadata.courge2016.rare.temp.D1, permutations = 9999)
    permRare.tempD2.obj <-adonis(decostand(asvs.comm.courge2016.rare.temp.D2, method="hellinger") ~ trt, 
                                 strata=metadata.courge2016.rare.temp.D2$block, data=metadata.courge2016.rare.temp.D2, permutations = 9999)
    permRare.tempD3.obj <-adonis(decostand(asvs.comm.courge2016.rare.temp.D3, method="hellinger") ~ trt, 
                                 strata=metadata.courge2016.rare.temp.D3$block, data=metadata.courge2016.rare.temp.D3, permutations = 9999)
    
    data.frame("R2.D1"=permRare.tempD1.obj$aov.tab$R2[1],
               "P.D1" = permRare.tempD1.obj$aov.tab$"Pr(>F)"[1],
               "R2.D2"= permRare.tempD2.obj$aov.tab$R2[1],
               "P.D2" = permRare.tempD2.obj$aov.tab$"Pr(>F)"[1],
               "R2.D3"= permRare.tempD3.obj$aov.tab$R2[1],
               "P.D3" = permRare.tempD3.obj$aov.tab$"Pr(>F)"[1], 
               "iteration" = paste("iteration_", j))
    #####################################################
  }
  stopCluster(cl)
  
  ####################################################
  ### permanova
  rare.df2016$rare <- factor(i)
  rare.df2016.final <- rbind(rare.df2016.final, rare.df2016)
  #
  ####################################################
  ### Dim
  #dimVSdate.df2016$rare <- factor(i)
  #dimVSdate.df2016.final <- rbind(dimVSdate.df2016.final, data.frame("Sample.nb.D1"= dimVSdate.df2016$Sample.nb.D1, "ASVs.nb.D1"= dimVSdate.df2016$ASVs.nb.D1,
  #                                                           "Sample.nb.D2"= dimVSdate.df2016$Sample.nb.D2, "ASVs.nb.D2"= dimVSdate.df2016$ASVs.nb.D2,
  #                                                           "Sample.nb.D3"= dimVSdate.df2016$Sample.nb.D3, "ASVs.nb.D3"= dimVSdate.df2016$ASVs.nb.D3,
  #                                                           "rare" = dimVSdate.df2016$rare))      
  ####################################################
  
  
  #### LOOP SENSOR ####
  print(paste("rare ", i))
  #####################
}
stopCluster(cl)
##################### END OF LOOP ###############################

dimVSdate.df2016.gg <- rbind(data.frame(sample=dimVSdate.df2016.final$Sample.nb.D1, ASV=dimVSdate.df2016.final$ASVs.nb.D1, group="Early", rare=dimVSdate.df2016.final$rare),
                         data.frame(sample=dimVSdate.df2016.final$Sample.nb.D2, ASV=dimVSdate.df2016.final$ASVs.nb.D2, group="Mid", rare=dimVSdate.df2016.final$rare),
                         data.frame(sample=dimVSdate.df2016.final$Sample.nb.D3, ASV=dimVSdate.df2016.final$ASVs.nb.D3, group="Late", rare=dimVSdate.df2016.final$rare))

rare.df2016.final.gg <-rbind(data.frame(R=rare.df2016.final$R2.D1, P=rare.df2016.final$P.D1, rare=rare.df2016.final$rare, group="Early"), 
                         data.frame(R=rare.df2016.final$R2.D2, P=rare.df2016.final$P.D2, rare=rare.df2016.final$rare, group="Mid"),
                         data.frame(R=rare.df2016.final$R2.D3, P=rare.df2016.final$P.D3, rare=rare.df2016.final$rare, group="Late"))

#1000 to 10 000 seq rarefaction
rare10k.df2016.final.gg <- rare.df2016.final.gg[which(as.numeric(as.character(rare.df2016.final.gg$rare))<10001),]

p8d2016 <-ggplot(rare10k.df2016.final.gg, aes(x=factor(rare), y=R, col=group, fill=group))+ geom_violin() + ylab("R2") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_vline(xintercept = 5, linetype="dotted") +
  ggtitle("A") + guides(col=F, fill=FALSE) + ylim(c(0, 0.254))

p9d2016 <-ggplot(rare10k.df2016.final.gg, aes(x=factor(rare), y=P, col=group, fill=group))+ geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_vline(xintercept = 5, linetype="dotted")+
  ggtitle("B") + guides(col=F, fill=FALSE) + ylim(c(0, 0.01))

p10d2016 <- ggplot(dimVSdate.df2016.gg, aes(x=factor(rare), y=sample, col=group, fill=group))+ geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_vline(xintercept = 5, linetype="dotted")+ ylim(c(50, 72)) + ggtitle("D") + guides(col=F, fill=F) 
```

###2017
```{r}
# Defining Dataframe
rare.df2017.final <- NULL
#dimVSdate.df2017.final <- NULL
#asvVSdate.df2017.final <- NULL


###Initiate rarefy loop
for (i in seq(1000,10000,1000)) {
  #Parallelize
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #to not overload your computer
  registerDoParallel(cl)
  
  ###################################
  # Defining DESIRED Dataframe
  rare.df2017 <- NULL   #Permanova 
  #dimVSdate.df2017 <- NULL  #Sample Stat
  #asvVSdate.df2017 <- NULL #ASV number
  ###################################
  
  ###################################
  ###Keep only samples with the current loop minimum of sequences (i)
  asvs.comm.courge2017.temp <- asvs.comm.courge2017[which(apply(asvs.comm.courge2017,1,sum)>=i),]
  asvs.comm.courge2017.temp <- asvs.comm.courge2017.temp[,apply(asvs.comm.courge2017.temp,2,sum)>1]
  #asvs.taxo.courge2017.temp <- asvs.taxo.courge2017[colnames(asvs.comm.courge2017.temp),]
  metadata.courge2017.temp <- metadata.courge2017[rownames(asvs.comm.courge2017.temp),]
  ###################################
  
  ###################################
  #### CHANGE THE LOOP OUTPUT FOR THE DESIRED Dataframe (see above)
  ###Rarefy with 100 iteration
  rare.df2017 <- foreach(j=1:100, .combine=rbind, .packages = "vegan") %dopar% { 
    asvs.comm.courge2017.rare.temp <- rrarefy(asvs.comm.courge2017.temp, sample=i)
    asvs.comm.courge2017.rare.temp <- asvs.comm.courge2017.rare.temp[,apply(asvs.comm.courge2017.rare.temp,2,sum)>1]
    metadata.courge2017.rare.temp <- metadata.courge2017.temp[rownames(asvs.comm.courge2017.rare.temp),]
    #asvs.taxo.courge2017.rare.temp <- asvs.taxo.courge2017.temp[colnames(asvs.comm.courge2017.rare.temp),]
    ###Adjust metadata
    metadata.courge2017.rare.temp$date <- factor(metadata.courge2017.rare.temp$date)
    metadata.courge2017.rare.temp$trt <- factor(metadata.courge2017.rare.temp$trt)
    metadata.courge2017.rare.temp$type <- factor(metadata.courge2017.rare.temp$type)
    metadata.courge2017.rare.temp$block <- factor(metadata.courge2017.rare.temp$block)
    metadata.courge2017.rare.temp$sub_ech <- factor(metadata.courge2017.rare.temp$sub_ech)
    #split by date
    ###Date1
    metadata.courge2017.rare.temp.D1 <- metadata.courge2017.rare.temp[metadata.courge2017.rare.temp$date == 'Date 1',]
    metadata.courge2017.rare.temp.D1$date <- factor(metadata.courge2017.rare.temp.D1$date)
    metadata.courge2017.rare.temp.D1$trt <- factor(metadata.courge2017.rare.temp.D1$trt)
    metadata.courge2017.rare.temp.D1$type <- factor(metadata.courge2017.rare.temp.D1$type)
    metadata.courge2017.rare.temp.D1$block <- factor(metadata.courge2017.rare.temp.D1$block)
    metadata.courge2017.rare.temp.D1$sub_ech <- factor(metadata.courge2017.rare.temp.D1$sub_ech)
    metadata.courge2017.rare.temp.D1$run <- factor(metadata.courge2017.rare.temp.D1$run)
    asvs.comm.courge2017.rare.temp.D1 <- asvs.comm.courge2017.rare.temp[rownames(metadata.courge2017.rare.temp.D1),]
    asvs.comm.courge2017.rare.temp.D1 <- asvs.comm.courge2017.rare.temp.D1[,apply(asvs.comm.courge2017.rare.temp.D1,2,sum)>0]
    ###Date2
    metadata.courge2017.rare.temp.D2 <- metadata.courge2017.rare.temp[metadata.courge2017.rare.temp$date == 'Date 2',]
    metadata.courge2017.rare.temp.D2$date <- factor(metadata.courge2017.rare.temp.D2$date)
    metadata.courge2017.rare.temp.D2$trt <- factor(metadata.courge2017.rare.temp.D2$trt)
    metadata.courge2017.rare.temp.D2$type <- factor(metadata.courge2017.rare.temp.D2$type)
    metadata.courge2017.rare.temp.D2$block <- factor(metadata.courge2017.rare.temp.D2$block)
    metadata.courge2017.rare.temp.D2$sub_ech <- factor(metadata.courge2017.rare.temp.D2$sub_ech)
    metadata.courge2017.rare.temp.D2$run <- factor(metadata.courge2017.rare.temp.D2$run)
    asvs.comm.courge2017.rare.temp.D2 <- asvs.comm.courge2017.rare.temp[rownames(metadata.courge2017.rare.temp.D2),]
    asvs.comm.courge2017.rare.temp.D2 <- asvs.comm.courge2017.rare.temp.D2[,apply(asvs.comm.courge2017.rare.temp.D2,2,sum)>0]
    ###Date3
    metadata.courge2017.rare.temp.D3 <- metadata.courge2017.rare.temp[metadata.courge2017.rare.temp$date == 'Date 3',]
    metadata.courge2017.rare.temp.D3$date <- factor(metadata.courge2017.rare.temp.D3$date)
    metadata.courge2017.rare.temp.D3$trt <- factor(metadata.courge2017.rare.temp.D3$trt)
    metadata.courge2017.rare.temp.D3$type <- factor(metadata.courge2017.rare.temp.D3$type)
    metadata.courge2017.rare.temp.D3$block <- factor(metadata.courge2017.rare.temp.D3$block)
    metadata.courge2017.rare.temp.D3$sub_ech <- factor(metadata.courge2017.rare.temp.D3$sub_ech)
    metadata.courge2017.rare.temp.D3$run <- factor(metadata.courge2017.rare.temp.D3$run)
    asvs.comm.courge2017.rare.temp.D3 <- asvs.comm.courge2017.rare.temp[rownames(metadata.courge2017.rare.temp.D3),]
    asvs.comm.courge2017.rare.temp.D3 <- asvs.comm.courge2017.rare.temp.D3[,apply(asvs.comm.courge2017.rare.temp.D3,2,sum)>0]
    ####################################################
    
    #####################################################
    #Dim : Sample Stat
    #data.frame("Sample.nb.D1"= dim(asvs.comm.courge2017.rare.temp.D1)[1], 
    #           "ASVs.nb.D1"= dim(asvs.comm.courge2017.rare.temp.D1)[2],
    #           "Sample.nb.D2"= dim(asvs.comm.courge2017.rare.temp.D2)[1], 
    #           "ASVs.nb.D2"= dim(asvs.comm.courge2017.rare.temp.D2)[2],
    #           "Sample.nb.D3"= dim(asvs.comm.courge2017.rare.temp.D3)[1], 
    #           "ASVs.nb.D3"= dim(asvs.comm.courge2017.rare.temp.D3)[2])
    #####################################################
    
    #####################################################
    #Permanova
    permRare.tempD1.obj <-adonis(decostand(asvs.comm.courge2017.rare.temp.D1, method="hellinger") ~ trt, 
                                 strata=metadata.courge2017.rare.temp.D1$block, data=metadata.courge2017.rare.temp.D1, permutations = 9999)
    permRare.tempD2.obj <-adonis(decostand(asvs.comm.courge2017.rare.temp.D2, method="hellinger") ~ trt, 
                                 strata=metadata.courge2017.rare.temp.D2$block, data=metadata.courge2017.rare.temp.D2, permutations = 9999)
    permRare.tempD3.obj <-adonis(decostand(asvs.comm.courge2017.rare.temp.D3, method="hellinger") ~ trt, 
                                 strata=metadata.courge2017.rare.temp.D3$block, data=metadata.courge2017.rare.temp.D3, permutations = 9999)
    
    data.frame("R2.D1"=permRare.tempD1.obj$aov.tab$R2[1],
               "P.D1" = permRare.tempD1.obj$aov.tab$"Pr(>F)"[1],
               "R2.D2"= permRare.tempD2.obj$aov.tab$R2[1],
               "P.D2" = permRare.tempD2.obj$aov.tab$"Pr(>F)"[1],
               "R2.D3"= permRare.tempD3.obj$aov.tab$R2[1],
               "P.D3" = permRare.tempD3.obj$aov.tab$"Pr(>F)"[1], 
               "iteration" = paste("iteration_", j))
    #####################################################
  }
  stopCluster(cl)
  
  ####################################################
  ### permanova
  rare.df2017$rare <- factor(i)
  rare.df2017.final <- rbind(rare.df2017.final, rare.df2017)
  #
  ####################################################
  ### Dim
  #dimVSdate.df2017$rare <- factor(i)
  #dimVSdate.df2017.final <- rbind(dimVSdate.df2017.final, data.frame("Sample.nb.D1"= dimVSdate.df2017$Sample.nb.D1, "ASVs.nb.D1"= dimVSdate.df2017$ASVs.nb.D1,
  #                                                           "Sample.nb.D2"= dimVSdate.df2017$Sample.nb.D2, "ASVs.nb.D2"= dimVSdate.df2017$ASVs.nb.D2,
  #                                                           "Sample.nb.D3"= dimVSdate.df2017$Sample.nb.D3, "ASVs.nb.D3"= dimVSdate.df2017$ASVs.nb.D3,
  #                                                           "rare" = dimVSdate.df2017$rare))      
  ####################################################
  
  
  #### LOOP SENSOR ####
  print(paste("rare ", i))
  #####################
}
stopCluster(cl)
##################### END OF LOOP ###############################

dimVSdate.df2017.gg <- rbind(data.frame(sample=dimVSdate.df2017.final$Sample.nb.D1, ASV=dimVSdate.df2017.final$ASVs.nb.D1, group="Early", rare=dimVSdate.df2017.final$rare),
                             data.frame(sample=dimVSdate.df2017.final$Sample.nb.D2, ASV=dimVSdate.df2017.final$ASVs.nb.D2, group="Mid", rare=dimVSdate.df2017.final$rare),
                             data.frame(sample=dimVSdate.df2017.final$Sample.nb.D3, ASV=dimVSdate.df2017.final$ASVs.nb.D3, group="Late", rare=dimVSdate.df2017.final$rare))

rare.df2017.final.gg <-rbind(data.frame(R=rare.df2017.final$R2.D1, P=rare.df2017.final$P.D1, rare=rare.df2017.final$rare, group="Early"), 
                             data.frame(R=rare.df2017.final$R2.D2, P=rare.df2017.final$P.D2, rare=rare.df2017.final$rare, group="Mid"),
                             data.frame(R=rare.df2017.final$R2.D3, P=rare.df2017.final$P.D3, rare=rare.df2017.final$rare, group="Late"))

#1000 to 10 000 seq rarefaction
rare10k.df2017.final.gg <- rare.df2017.final.gg[which(as.numeric(as.character(rare.df2017.final.gg$rare))<10001),]

p8d2017 <-ggplot(rare10k.df2017.final.gg, aes(x=factor(rare), y=R, col=group, fill=group))+ geom_violin() + ylab("R2") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_vline(xintercept = 5, linetype="dotted") +
  ggtitle("A") + guides(col=F, fill=FALSE) + ylim(c(0, 0.254))

p9d2017 <-ggplot(rare10k.df2017.final.gg, aes(x=factor(rare), y=P, col=group, fill=group))+ geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_vline(xintercept = 5, linetype="dotted")+
  ggtitle("B") + guides(col=F, fill=FALSE) + ylim(c(0, 0.01))

p10d2017 <- ggplot(dimVSdate.df2017.gg, aes(x=factor(rare), y=sample, col=group, fill=group))+ geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_vline(xintercept = 5, linetype="dotted")+ ylim(c(50, 72)) + ggtitle("D") + guides(col=F, fill=F) 
```

## (line 204) Supplemental Figure 5
```{r}
multiplot(p8d2016, p8d2017, p9d2016, p9d2017, p10d2016, p10d2017, cols = 3)
```
