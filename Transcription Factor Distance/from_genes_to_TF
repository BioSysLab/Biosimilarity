library(tidyverse)
library(viper)
library(hgu133a.db)
library(hthgu133a.db)
library(CARNIVAL)

## weighting function, is used to get an aggregated perturbation column for the drugs that are replicated

weighted_mat <- function(m, clipLowWt=TRUE, lowThreshWT=0.1, clipLowCC=TRUE, lowThreshCC=0, metric="avg") {
  cc = cor(m, method="spearman") # pair-wise correlations
  if(clipLowCC) { # Trim low values
    cc[cc<lowThreshCC] = lowThreshCC
  }
  
  wt = 0.5 * colSums(cc) # Per-sample values
  if(clipLowWt) { # Trim low values
    wt[wt<lowThreshWT] = lowThreshWT
  }
  
  # Normalize the weights
  if(metric=="avg") {
    sumWT = sum(abs(wt))
  } else {
    sumWT = sqrt(sum(wt^2))
  }
  normWT = wt / sumWT  # Normalized weights
  
  # Return the scaled input values
  return(m %*% normWT)
}


## instance matrices of every cell line

MCF7_genes <- readRDS("MCF7_2616_instances_for_1211_unique_compounds.rds")
HL60_genes <- readRDS("HL60_1159_instances_for_1078_unique_compounds.rds")
PC3_genes <- readRDS("PC3_1543_instances_for_1161_unique_compounds.rds")

## total perturbations for arrays HG-U133A and HT_HG-U133A

## array 1 is the HG-U133A

array_1_tot_pert <- readRDS("array_1_tot_pert.RDS")

## array 2 is the HT_HG-U133A

array_2_tot_pert <- readRDS("array_2_tot_pert.RDS")

ar_1 <- as.data.frame(array_1_tot_pert)

ar_2<- as.data.frame(array_2_tot_pert)

## MCF7 perturbations, according to array 2, since this is the only array at MCF7

MCF7_genes$perturbation_scan_id <- paste0(MCF7_genes$perturbation_scan_id,".CEL","")
MCF7_genes$perturbation_scan_id <- str_remove(MCF7_genes$perturbation_scan_id,"'")
MCF7_pert_idx <- which(colnames(ar_2)%in%MCF7_genes$perturbation_scan_id)
MCF7_pert <- ar_2[,MCF7_pert_idx]

## build a perturbation matrix with unique perturbations using the replicates' correlations to get a 
## final perturbation matrix

MCF7_final_pert <- matrix(,ncol = length(unique(MCF7_genes$cmap_name)), nrow = nrow(ar_2))
colnames(MCF7_final_pert) <- unique(MCF7_genes$cmap_name)
for (i in 1:ncol(MCF7_final_pert)){
  drug <- unique(MCF7_genes$cmap_name)[i]
  replicates <- which(MCF7_genes$cmap_name==as.character(drug))
  files <- MCF7_genes$perturbation_scan_id[replicates]
  pert_idx <- which(colnames(MCF7_pert)%in%files)
  if (length(replicates)>1){
    merged_pert <- weighted_mat(m = as.matrix(MCF7_pert[,pert_idx]))
    MCF7_final_pert[,i] <- merged_pert
  } else{
    MCF7_final_pert[,i] <- MCF7_pert[,pert_idx]
  }
}

saveRDS(MCF7_final_pert,"MCF7_final_pert.RDS")

## HL60 perturbations. HL60 cell line contains both array 1 and array 2

HL60_genes$perturbation_scan_id <- paste0(HL60_genes$perturbation_scan_id,".CEL","")
HL60_genes$perturbation_scan_id <- str_remove(HL60_genes$perturbation_scan_id,"'")

## perturbations regarding array 1

HL60_pert_ar1_idx <- which(colnames(ar_1)%in%HL60_genes$perturbation_scan_id)
HL60_pert_ar1 <- ar_1[,HL60_pert_ar1_idx]
HL60_inst_ar1 <- HL60_genes[HL60_genes$array3=="HG-U133A",]

## build a perturbation matrix with unique perturbations using the replicates' correlations to get a 
## final perturbation matrix

HL60_final_pert_ar1 <- matrix(,ncol = length(unique(HL60_inst_ar1$cmap_name)), nrow = nrow(ar_1))
colnames(HL60_final_pert_ar1) <- unique(HL60_inst_ar1$cmap_name)
for (i in 1:ncol(HL60_final_pert_ar1)){
  drug <- unique(HL60_inst_ar1$cmap_name)[i]
  replicates <- which(HL60_inst_ar1$cmap_name==as.character(drug))
  files <- HL60_inst_ar1$perturbation_scan_id[replicates]
  pert_idx <- which(colnames(HL60_pert_ar1)%in%files)
  if (length(replicates)>1){
    merged_pert <- weighted_mat(m = as.matrix(HL60_pert_ar1[,pert_idx]))
    HL60_final_pert_ar1[,i] <- merged_pert
  } else{
    HL60_final_pert_ar1[,i] <- HL60_pert_ar1[,pert_idx]
  }
}

saveRDS(HL60_final_pert_ar1,"HL60_final_pert_ar1.RDS")

## perturbations regarding array 2

HL60_pert_ar2_idx <- which(colnames(ar_2)%in%HL60_genes$perturbation_scan_id)
HL60_pert_ar2 <- ar_2[,HL60_pert_ar2_idx]
HL60_inst_ar2 <- HL60_genes[HL60_genes$array3=="HT_HG-U133A",]

## build a perturbation matrix with unique perturbations using the replicates' correlations to get a 
## final perturbation matrix

HL60_final_pert_ar2 <- matrix(,ncol = length(unique(HL60_inst_ar2$cmap_name)), nrow = nrow(ar_2))
colnames(HL60_final_pert_ar2) <- unique(HL60_inst_ar2$cmap_name)
for (i in 1:ncol(HL60_final_pert_ar2)){
  drug <- unique(HL60_inst_ar2$cmap_name)[i]
  replicates <- which(HL60_inst_ar2$cmap_name==as.character(drug))
  files <- HL60_inst_ar2$perturbation_scan_id[replicates]
  pert_idx <- which(colnames(HL60_pert_ar2)%in%files)
  if (length(replicates)>1){
    merged_pert <- weighted_mat(m = as.matrix(HL60_pert_ar2[,pert_idx]))
    HL60_final_pert_ar2[,i] <- merged_pert
  } else{
    HL60_final_pert_ar2[,i] <- HL60_pert_ar2[,pert_idx]
  }
}

saveRDS(HL60_final_pert_ar2,"HL60_final_pert_ar2.RDS")


## PC3 perturbations, according to array 2, since this is the only array at PC3

PC3_genes$perturbation_scan_id <- paste0(PC3_genes$perturbation_scan_id,".CEL","")
PC3_genes$perturbation_scan_id <- str_remove(PC3_genes$perturbation_scan_id,"'")
PC3_pert_idx <- which(colnames(ar_2)%in%PC3_genes$perturbation_scan_id)
PC3_pert <- ar_2[,PC3_pert_idx]

## build a perturbation matrix with unique perturbations using the replicates' correlations to get a 
## final perturbation matrix

PC3_final_pert <- matrix(,ncol = length(unique(PC3_genes$cmap_name)), nrow = nrow(ar_2))
colnames(PC3_final_pert) <- unique(PC3_genes$cmap_name)
for (i in 1:ncol(PC3_final_pert)){
  drug <- unique(PC3_genes$cmap_name)[i]
  replicates <- which(PC3_genes$cmap_name==as.character(drug))
  files <- PC3_genes$perturbation_scan_id[replicates]
  pert_idx <- which(colnames(PC3_pert)%in%files)
  if (length(replicates)>1){
    merged_pert <- weighted_mat(m = as.matrix(PC3_pert[,pert_idx]))
    PC3_final_pert[,i] <- merged_pert
  } else{
    PC3_final_pert[,i] <- PC3_pert[,pert_idx]
  }
}

saveRDS(PC3_final_pert,"PC3_final_pert.RDS")


## Apply Viper to MCF7 to get the Transcription Factors

MCF7_pert_v <- MCF7_final_pert
rownames(MCF7_pert_v) <- rownames(MCF7_pert)
probes <- rownames(MCF7_pert_v)
anno <- AnnotationDbi::select(hthgu133a.db,
                              keys = probes,
                              columns = c("SYMBOL", "GENENAME","PROBEID"),
                              keytype = "PROBEID")
anno <- anno[,c("SYMBOL","PROBEID")]
anno <- anno%>%
  group_by(SYMBOL)%>%
  mutate(times1=n_distinct(PROBEID))%>%
  ungroup()%>%
  group_by(PROBEID)%>%
  mutate(times2=n_distinct(SYMBOL))%>%
  ungroup()

anno <- anno %>%
  filter(!is.na(SYMBOL))%>%
  filter(times1==1)%>%
  filter(times2==1)

symbol_sel <- which(rownames(MCF7_pert_v)%in%anno$PROBEID)

MCF7_pert_v <- MCF7_pert_v[symbol_sel,]

rownames(MCF7_pert_v) <- anno$SYMBOL

saveRDS(MCF7_pert_v,"MCF7_pert_for_viper.RDS")

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons

viper_MCF7 <- runDoRothEA(MCF7_pert_v, regulon=viper_regulon, confidence_level=c('A','B','C'))
saveRDS(viper_MCF7,"viper_MCF7.RDS")


## Apply Viper to HL60 to get the Transcription Factors 
## (We have two arrays here, HT_HG-U133A & HG-U133A)

## v1 for HG-U133A

HL60_genes_v1 <- HL60_inst_ar1
HL60_pert_v1 <- HL60_final_pert_ar1
rownames(HL60_pert_v1) <- rownames(HL60_pert_ar1)
probes <- rownames(HL60_pert_v1)
anno <- AnnotationDbi::select(hgu133a.db,
                              keys = probes,
                              columns = c("SYMBOL", "GENENAME","PROBEID"),
                              keytype = "PROBEID")
anno <- anno[,c("SYMBOL","PROBEID")]
anno <- anno%>%
  group_by(SYMBOL)%>%
  mutate(times1=n_distinct(PROBEID))%>%
  ungroup()%>%
  group_by(PROBEID)%>%
  mutate(times2=n_distinct(SYMBOL))%>%
  ungroup()

anno <- anno %>%
  filter(!is.na(SYMBOL))%>%
  filter(times1==1)%>%
  filter(times2==1)

symbol_sel <- which(rownames(HL60_pert_v1)%in%anno$PROBEID)

HL60_pert_v1 <- HL60_pert_v1[symbol_sel,]

rownames(HL60_pert_v1) <- anno$SYMBOL

saveRDS(HL60_pert_v1,"HL60_pert_v1_for_viper.RDS")

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons

viper_HL60_v1 <- runDoRothEA(HL60_pert_v1, regulon=viper_regulon, confidence_level=c('A','B','C'))
saveRDS(viper_HL60_v1,"viper_HL60_v1.RDS")

## v2 for HT_HG-U133A

HL60_genes_v2 <- HL60_inst_ar2
HL60_pert_v2 <- HL60_final_pert_ar2
rownames(HL60_pert_v2) <- rownames(HL60_pert_ar2)
probes <- rownames(HL60_pert_v2)
anno <- AnnotationDbi::select(hthgu133a.db,
                              keys = probes,
                              columns = c("SYMBOL", "GENENAME","PROBEID"),
                              keytype = "PROBEID")
anno <- anno[,c("SYMBOL","PROBEID")]
anno <- anno%>%
  group_by(SYMBOL)%>%
  mutate(times1=n_distinct(PROBEID))%>%
  ungroup()%>%
  group_by(PROBEID)%>%
  mutate(times2=n_distinct(SYMBOL))%>%
  ungroup()

anno <- anno %>%
  filter(!is.na(SYMBOL))%>%
  filter(times1==1)%>%
  filter(times2==1)

symbol_sel <- which(rownames(HL60_pert_v2)%in%anno$PROBEID)

HL60_pert_v2 <- HL60_pert_v2[symbol_sel,]

rownames(HL60_pert_v2) <- anno$SYMBOL

saveRDS(HL60_pert_v2,"HL60_pert_v2_for_viper.RDS")

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons

viper_HL60_v2 <- runDoRothEA(HL60_pert_v2, regulon=viper_regulon, confidence_level=c('A','B','C'))
saveRDS(viper_HL60_v2,"viper_HL60_v2.RDS")

## complete HL60 Viper, bind both of the calculated TFs matrices (v1 and v2)

viper_HL60 <- cbind(viper_HL60_v1,viper_HL60_v2)
saveRDS(viper_HL60,"viper_HL60_full.RDS")

## Apply Viper to PC3 to get the Transcription Factors

PC3_pert_v <- PC3_final_pert
rownames(PC3_pert_v) <- rownames(PC3_pert)
probes <- rownames(PC3_pert_v)
anno <- AnnotationDbi::select(hthgu133a.db,
                              keys = probes,
                              columns = c("SYMBOL", "GENENAME","PROBEID"),
                              keytype = "PROBEID")
anno <- anno[,c("SYMBOL","PROBEID")]
anno <- anno%>%
  group_by(SYMBOL)%>%
  mutate(times1=n_distinct(PROBEID))%>%
  ungroup()%>%
  group_by(PROBEID)%>%
  mutate(times2=n_distinct(SYMBOL))%>%
  ungroup()

anno <- anno %>%
  filter(!is.na(SYMBOL))%>%
  filter(times1==1)%>%
  filter(times2==1)

symbol_sel <- which(rownames(PC3_pert_v)%in%anno$PROBEID)

PC3_pert_v <- PC3_pert_v[symbol_sel,]

rownames(PC3_pert_v) <- anno$SYMBOL

saveRDS(PC3_pert_v,"PC3_pert_for_viper.RDS")

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons

viper_PC3 <- runDoRothEA(PC3_pert_v, regulon=viper_regulon, confidence_level=c('A','B','C'))
saveRDS(viper_PC3,"viper_PC3.RDS")
