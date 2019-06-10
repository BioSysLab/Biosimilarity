library(readxl)
library(tidyverse)
library(GeneExpressionSignature)

### read input data

instances_2 <- read_excel("cmap_instances_02.xls")
instances_2 <- instances_2 %>%
  filter(!is.na(cmap_name))
perturbations_2<- read.delim("rankMatrix.txt", header=TRUE, sep="\t")

### fix dimensions and names

pert_new <- perturbations_2[,-1]
pert_new <- pert_new[,-6101]
pert_new <- as.matrix(pert_new)
colnames(pert_new) <- gsub(pattern = "X",replacement = "",x = colnames(pert_new))
rownames(pert_new) <- perturbations_2[,1]

### find indices of cell lines

MCF7_idx <-  which(grepl(pattern = "MCF7",x = instances_2$cell2) == TRUE)
HL60_idx <-  which(grepl(pattern = "HL60",x = instances_2$cell2) == TRUE)
PC3_idx <-  which(grepl(pattern = "PC3",x = instances_2$cell2) == TRUE)
SKMEL5_idx <-  which(grepl(pattern = "SKMEL5",x = instances_2$cell2) == TRUE) ### SKMEL will be left out, only 18 instances

### split into 3 matrices, one for each cell

MCF7_pert <- pert_new[,MCF7_idx]
HL60_pert <- pert_new[,HL60_idx]
PC3_pert <- pert_new[,PC3_idx]

### filter MCF7
### keep only 1 array
MCF7_f1 <- instances_2[MCF7_idx,] %>% filter(array3 == "HT_HG-U133A")

### lost due to f1

lost_MCF7_f1 <- anti_join(instances_2[MCF7_idx,],MCF7_f1)

### remove 12 hours, all were removed from previous filter

MCF7_f2 <- MCF7_f1 %>% filter(`duration (h)` != "12")

### create bins for concentrations

MCF7_f2$bins <- 0

MCF7_f2 <- MCF7_f2 %>% 
  mutate(bins = if_else(`concentration (M)`>=10^(-8) & `concentration (M)` < 5*10^(-7),true = 10,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-7) & `concentration (M)` < 5*10^(-6),true = 9,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-6) & `concentration (M)` < 5*10^(-5),true = 8,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-5) & `concentration (M)` < 5*10^(-4),true = 7,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-4) & `concentration (M)` < 5*10^(-3),true = 6,false = bins))

MCF7_test <- MCF7_f2 %>%
  group_by(cmap_name) %>% nest(bins)

counts <- c(0,0,0,0,0)
names(counts) <- c("10","9","8","7","6")

for (i in 1:nrow(MCF7_test)) {
  
  if (10 %in% unlist(MCF7_test$data[[i]])) {
    counts[1] <- counts[1] + 1
  }
  if (9 %in% unlist(MCF7_test$data[[i]])) {
    counts[2] <- counts[2] + 1
  }
  if (8 %in% unlist(MCF7_test$data[[i]])) {
    counts[3] <- counts[3] + 1
  }
  if (7 %in% unlist(MCF7_test$data[[i]])) {
    counts[4] <- counts[4] + 1
  }
  if (6 %in% unlist(MCF7_test$data[[i]])) {
    counts[5] <- counts[5] + 1
  }
}

counts <- as.data.frame(counts) %>%
  rownames_to_column("bins")
counts$bins <- as.numeric(counts$bins)
MCF7_f2 <- left_join(MCF7_f2,counts)

MCF7_f2 <- MCF7_f2 %>% group_by(cmap_name) %>%
  filter(counts == max(counts)) %>% ungroup()


#### filter the rank matrix for mcf7 to have the intances that we kept

mcf7_idx <- which(colnames(MCF7_pert) %in% MCF7_f2$instance_id)
MCF7_pert <- MCF7_pert[,mcf7_idx]
########################## MCF7 DISTANCES #################################### 

pheno2 <- as.matrix(as.factor(MCF7_f2$cmap_name)[MCF7_f2$instance_id %in% colnames(MCF7_pert)])
rownames(pheno2) <- colnames(MCF7_pert)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno2))
exprset <- new("ExpressionSet",exprs=MCF7_pert,phenoData=pheno_new)
a <- RankMerging(exprset, MergingDistance = "Spearman")
ds_MCF7 <- ScoreGSEA(a,250,"avg")
saveRDS(ds_MCF7,"MCF7_GSEA_dist.RDS")
saveRDS(MCF7_f2,"MCF7_2616_instances_for_1211_unique_compounds.rds")


### filter HL60
### keep both arrays here, because HG-U133A 344 instances, HT_HG-U133A 885 instances
HL60_f1 <- instances_2[HL60_idx,]

### nothing is lost due to f1

### all have 6h duration

### create bins for concentrations

HL60_f1$bins <- 0

HL60_f2 <- HL60_f1 %>% 
  mutate(bins = if_else(`concentration (M)`>=10^(-8) & `concentration (M)` < 5*10^(-7),true = 10,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-7) & `concentration (M)` < 5*10^(-6),true = 9,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-6) & `concentration (M)` < 5*10^(-5),true = 8,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-5) & `concentration (M)` < 5*10^(-4),true = 7,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-4) & `concentration (M)` < 5*10^(-3),true = 6,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-3) & `concentration (M)` < 5*10^(-2),true = 5,false = bins))

HL60_test <- HL60_f2 %>%
  group_by(cmap_name) %>% nest(bins)

counts <- c(0,0,0,0,0,0)
names(counts) <- c("10","9","8","7","6","5")

for (i in 1:nrow(HL60_test)) {
  
  if (10 %in% unlist(HL60_test$data[[i]])) {
    counts[1] <- counts[1] + 1
  }
  if (9 %in% unlist(HL60_test$data[[i]])) {
    counts[2] <- counts[2] + 1
  }
  if (8 %in% unlist(HL60_test$data[[i]])) {
    counts[3] <- counts[3] + 1
  }
  if (7 %in% unlist(HL60_test$data[[i]])) {
    counts[4] <- counts[4] + 1
  }
  if (6 %in% unlist(HL60_test$data[[i]])) {
    counts[5] <- counts[5] + 1
  }
  if (5 %in% unlist(HL60_test$data[[i]])) {
    counts[6] <- counts[6] + 1
  }
}

counts <- as.data.frame(counts) %>%
  rownames_to_column("bins")
counts$bins <- as.numeric(counts$bins)
HL60_f2 <- left_join(HL60_f2,counts)

HL60_f2 <- HL60_f2 %>% group_by(cmap_name) %>%
  filter(counts == max(counts)) %>% ungroup()

### extra filter, to check drugs that are in both arrays. We will keep the ones that are in the "HT_HG-U133A"
HL_60_HT_HG_U133A <- HL60_f2[HL60_f2$array3=="HT_HG-U133A",]
HL_60_HG_U133A <- HL60_f2[HL60_f2$array3=="HG-U133A",]
to_remove <- which(HL_60_HG_U133A$cmap_name%in%HL_60_HT_HG_U133A$cmap_name)
HL_60_HG_U133A <- HL_60_HG_U133A[-to_remove,]

HL60_f2 <-rbind(HL_60_HT_HG_U133A,HL_60_HG_U133A)
#### filter the rank matrix for HL60 to have the intances that we kept

hl60_idx <- which(colnames(HL60_pert) %in% HL60_f2$instance_id)
HL60_pert <- HL60_pert[,hl60_idx]
########################## HL60 DISTANCES #################################### 

pheno2 <- as.matrix(as.factor(HL60_f2$cmap_name)[HL60_f2$instance_id %in% colnames(HL60_pert)])
rownames(pheno2) <- colnames(HL60_pert)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno2))
exprset <- new("ExpressionSet",exprs=HL60_pert,phenoData=pheno_new)
a <- RankMerging(exprset, MergingDistance = "Spearman")
ds_HL60 <- ScoreGSEA(a,250,"avg")

saveRDS(ds_HL60,"HL60_GSEA_dist.RDS")
saveRDS(HL60_f2,"HL60_1159_instances_for_1078_unique_compounds.rds")


### filter PC3
### keep only 1 array, HT_HG-U133A has 1617 instances, while HG-U133A has 124 instances

PC3_f1 <- instances_2[PC3_idx,] %>% filter(array3 == "HT_HG-U133A")

### lost due to f1

lost_PC3_f1 <- anti_join(instances_2[PC3_idx,],PC3_f1)

### all are 6h duration

### create bins for concentrations

PC3_f1$bins <- 0

PC3_f2 <- PC3_f1 %>% 
  mutate(bins = if_else(`concentration (M)`>=10^(-8) & `concentration (M)` < 5*10^(-7),true = 10,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-7) & `concentration (M)` < 5*10^(-6),true = 9,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-6) & `concentration (M)` < 5*10^(-5),true = 8,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-5) & `concentration (M)` < 5*10^(-4),true = 7,false = bins),
         bins = if_else(`concentration (M)`>=5*10^(-4) & `concentration (M)` < 5*10^(-3),true = 6,false = bins))

PC3_test <- PC3_f2 %>%
  group_by(cmap_name) %>% nest(bins)

counts <- c(0,0,0,0,0)
names(counts) <- c("10","9","8","7","6")

for (i in 1:nrow(PC3_test)) {
  
  if (10 %in% unlist(PC3_test$data[[i]])) {
    counts[1] <- counts[1] + 1
  }
  if (9 %in% unlist(PC3_test$data[[i]])) {
    counts[2] <- counts[2] + 1
  }
  if (8 %in% unlist(PC3_test$data[[i]])) {
    counts[3] <- counts[3] + 1
  }
  if (7 %in% unlist(PC3_test$data[[i]])) {
    counts[4] <- counts[4] + 1
  }
  if (6 %in% unlist(PC3_test$data[[i]])) {
    counts[5] <- counts[5] + 1
  }
}

counts <- as.data.frame(counts) %>%
  rownames_to_column("bins")
counts$bins <- as.numeric(counts$bins)
PC3_f2 <- left_join(PC3_f2,counts)

PC3_f2 <- PC3_f2 %>% group_by(cmap_name) %>%
  filter(counts == max(counts)) %>% ungroup()


#### filter the rank matrix for pc3 to have the intances that we kept

pc3_idx <- which(colnames(PC3_pert) %in% PC3_f2$instance_id)
PC3_pert <- PC3_pert[,pc3_idx]
########################## PC3 DISTANCES #################################### 

pheno2 <- as.matrix(as.factor(PC3_f2$cmap_name)[PC3_f2$instance_id %in% colnames(PC3_pert)])
rownames(pheno2) <- colnames(PC3_pert)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno2))
exprset <- new("ExpressionSet",exprs=PC3_pert,phenoData=pheno_new)
a <- RankMerging(exprset, MergingDistance = "Spearman")
ds_PC3 <- ScoreGSEA(a,250,"avg")
saveRDS(ds_PC3,"PC3_GSEA_dist.RDS")
saveRDS(PC3_f2,"PC3_1543_instances_for_1161_unique_compounds.rds")
