library(rdist)
library(tidyverse)
library(dplyr)
library(GeneExpressionSignature)

## Compute Distance with GSEA score for MCF7

viper_MCF7 <- readRDS("viper_MCF7.RDS")

## make a Rank Matrix for the MCF7 Viper results, rank(-viper_MCF7), because we need descending, while rank()
## provides ascending order

viper_MCF7_r <- apply(X = -viper_MCF7, FUN = rank, MARGIN = 2)

pheno <- as.matrix(colnames(viper_MCF7_r))
rownames(pheno) <- colnames(viper_MCF7_r)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno))
exprset <- new("ExpressionSet",exprs=as.matrix(viper_MCF7_r),phenoData=pheno_new)

viper_MCF7_ds <- ScoreGSEA(exprset,20,"avg")
saveRDS(viper_MCF7_ds,"viper_MCF7_GSEA_distance.RDS")


## Compute Distance with GSEA score for HL60

viper_HL60 <- readRDS("viper_HL60_full.RDS")

## because the viper dataset came from a matrix binding, according to the arrays, the compound names and columns
## positions have to changed accordingly to the ones in the gene GSEA, to be aligned in the same positions  

HL60_gene_GSEA <- readRDS("HL60_GSEA_dist.RDS")
comp_name_order <- colnames(HL60_gene_GSEA)

viper_HL60_2 <- viper_HL60[,comp_name_order]

## make a Rank Matrix for the HL60 Viper results, rank(-viper_HL60), because we need descending, while rank()
## provides ascending order

viper_HL60_r <- apply(X = -viper_HL60_2, FUN = rank, MARGIN = 2)

pheno <- as.matrix(colnames(viper_HL60_r))
rownames(pheno) <- colnames(viper_HL60_r)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno))
exprset <- new("ExpressionSet",exprs=as.matrix(viper_HL60_r),phenoData=pheno_new)

viper_HL60_ds <- ScoreGSEA(exprset,20,"avg")
saveRDS(viper_HL60_ds,"viper_HL60_GSEA_distance.RDS")

## Compute Distance with GSEA score for PC3

viper_PC3 <- readRDS("viper_PC3.RDS")

## make a Rank Matrix for the PC3 Viper results, rank(-viper_PC3), because we need descending, while rank()
## provides ascending order

viper_PC3_r <- apply(X = -viper_PC3, FUN = rank, MARGIN = 2)

pheno <- as.matrix(colnames(viper_PC3_r))
rownames(pheno) <- colnames(viper_PC3_r)
pheno_new <- new("AnnotatedDataFrame",data=as.data.frame(pheno))
exprset <- new("ExpressionSet",exprs=as.matrix(viper_PC3_r),phenoData=pheno_new)

viper_PC3_ds <- ScoreGSEA(exprset,20,"avg")
saveRDS(viper_PC3_ds,"viper_PC3_GSEA_distance.RDS")
