## Calculate the Euclidean distances for every Viper cell line result

## Euclidean Distance for MCF7

viper_MCF7 <- readRDS("viper_MCF7.RDS")
viper_MCF7_eucl_distance <- as.matrix(dist(x=t(viper_MCF7), method = "euclidean"))
saveRDS(viper_MCF7_eucl_distance,"viper_MCF7_eucl_dist.rds")

## Euclidean Distance for HL60

viper_HL60 <- readRDS("viper_HL60_full.RDS")
viper_HL60_eucl_distance <- as.matrix(dist(x=t(viper_HL60), method = "euclidean"))
saveRDS(viper_HL60_eucl_distance,"viper_HL60_eucl_dist.rds")

## Euclidean Distance for PC3

viper_PC3 <- readRDS("viper_PC3.RDS")
viper_PC3_eucl_distance <- as.matrix(dist(x=t(viper_PC3), method = "euclidean"))
saveRDS(viper_PC3_eucl_distance,"viper_PC3_eucl_dist.rds")
