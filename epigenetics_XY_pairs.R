#!/usr/bin/env Rscript

library('glmm')
library('effects')
library(data.table)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(lmerTest)
library(lme4)
library('MuMIn')
library(grafify)
library(emmeans)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(lmtest)

#------------------------------
# annotation file
#------------------------------
gff3 <- fread("genes.gff3", h=T)
gff3


#---------------------------------
# XY pairs sRNA mapping by strata
#---------------------------------
sRNA_data_raw_Y_males <- fread("y_v4_raw_gene_mapping_inds.tsv")
sRNA_data_tpm_Y_males <- fread("y_v4_tpm_gene_mapping_inds.tsv")
sRNA_data_Y_males <- merge(sRNA_data_raw_Y_males, sRNA_data_tpm_Y_males)
sRNA_data_Y_males[, chr := rep('Y', nrow(sRNA_data_Y_males))]
setcolorder(sRNA_data_Y_males, c("gene", "individual", "counts", "TPM", "tissue", "sex", "path", "chr"))
names(sRNA_data_Y_males) <- c("geneY", "individual", "sRNA_count", "TPM", "tissue", "sex", "path", "chr")
sRNA_data_Y_males[, geneY := strsplit(geneY, '.', fixed=T)[[1]][1], by=1:nrow(sRNA_data_Y_males)]
sRNA_data_Y_males

# whole genome small RNA mapping (this file does not contain Y sRNA mapping)
sRNA_data <- fread("v4_regression_input_sRNA_gene_mapping.tsv")
names(sRNA_data) <- c("GeneName", "individual", "sRNA_count", "TPM", "tissue", "sex", "path")
sRNA_data[, chr := strsplit(GeneName, '_', fixed=T)[[1]][2], by=1:nrow(sRNA_data)]
table(sRNA_data$chr)
sRNA_data[is.na(chr)]
sRNA_data[, GeneName := strsplit(GeneName, '.', fixed=T)[[1]][1], by=1:nrow(sRNA_data)]
sRNA_data
table(sRNA_data$chr)

# XY pairs
# 951 XY pair information from less strict definition (no need to have 1-1 ortholog with outgroups)
XYpairs <- fread('XY_pairs.txt')
# consider genes starting above 330Mb as PAR --> 910 gene pairs remaining our of 951
XYpairs <- XYpairs[end < 330000000]
table(XYpairs$strata)
XYpairs[, GeneName := strsplit(geneX, '.', fixed=T)[[1]][1], by=1:nrow(XYpairs)]
XYpairs[, geneY := strsplit(geneY, '.', fixed=T)[[1]][1], by=1:nrow(XYpairs)]
table(XYpairs$strata)

# extract X genes data
X_sRNA_data <- merge(sRNA_data, XYpairs, by=c("GeneName"))
X_sRNA_data <- unique(X_sRNA_data)
table(table(X_sRNA_data$GeneName))
table(X_sRNA_data$chr.x)
X_sRNA_data[, c("chr.x", "chr.y") := NULL]
X_sRNA_data[, chr := rep("X", nrow(X_sRNA_data))]
names(X_sRNA_data)
X_sRNA_data

# extract Y genes data
Y_sRNA_data <- merge(sRNA_data_Y_males, XYpairs, by=c("geneY"))
Y_sRNA_data <- unique(Y_sRNA_data)
table(table(Y_sRNA_data$GeneName))
Y_sRNA_data[, c("chr.x", "chr.y") := NULL]
Y_sRNA_data[, chr := rep("Y", nrow(Y_sRNA_data))]
setcolorder(Y_sRNA_data, names(X_sRNA_data))
Y_sRNA_data

# merge X and Y data
XY_sRNA_data <- rbindlist(list(X_sRNA_data, Y_sRNA_data))
XY_sRNA_data[, Chr_sex := paste(chr, ' in ', sex, 's', sep=''), by=1:nrow(XY_sRNA_data)]
XY_sRNA_data[, strata2 := gsub('s', 'Stratum ', gsub('S', 'Stratum ', gsub('a', '', gsub('b', '', gsub('c', '', strata))))), by=1:nrow(XY_sRNA_data)]
XY_sRNA_data
table(XY_sRNA_data$strata)
table(XY_sRNA_data$strata2)

# divide sRNA by 2 for X in females
XY_sRNA_data[, corrected_TPM := if (sex == "female") {TPM/2} else {TPM}, by=1:nrow(XY_sRNA_data)]
XY_sRNA_data[, corrected_count := if (sex == "female") {round(sRNA_count/2)} else {sRNA_count}, by=1:nrow(XY_sRNA_data)]

fwrite(XY_sRNA_data, "XY_sRNA_data.txt", sep="\t")
XY_sRNA_data <- fread("XY_sRNA_data.txt", h=T)

# PAR genes start > 330000000 pb
gff3 <- fread("genes.gff3", h=T)
gff3 <- gff3[(chr == "chrX")&(min > 330000000)]
gff3
sRNA_data <- fread("v4_regression_input_sRNA_gene_mapping.tsv")
names(sRNA_data) <- c("GeneName", "individual", "sRNA_count", "TPM", "tissue", "sex", "path")
sRNA_data[, GeneName := strsplit(GeneName, '.', fixed=T)[[1]][1], by=1:nrow(sRNA_data)]
sRNA_data
sRNA_data_PAR <- merge(sRNA_data, gff3, by = "GeneName")
length(table(sRNA_data_PAR$GeneName))
sRNA_data_PAR[, strata := rep("PAR", nrow(sRNA_data_PAR))]
sRNA_data_PAR[, geneX := GeneName]
sRNA_data_PAR[, geneY := GeneName]
names(sRNA_data_PAR) <- c("GeneName", "individual", "sRNA_count", "TPM", "tissue", "sex", "path", "chr", "start", "end", "strand", "strata", "geneX", "geneY")
sRNA_data_PAR[, chr := gsub('chr', '', chr), by=1:nrow(sRNA_data_PAR)]
sRNA_data_PAR[, Chr_sex := paste(chr, ' in ', sex, 's', sep=''), by=1:nrow(sRNA_data_PAR)]
names(sRNA_data_PAR)
sRNA_data_PAR[, strand := NULL]
sRNA_data_PAR[, strata2 :=  strata]
sRNA_data_PAR[, corrected_TPM := TPM/2]
sRNA_data_PAR[, corrected_count := sRNA_count/2]
setcolorder(sRNA_data_PAR, names(XY_sRNA_data))
XY_sRNA_data <- rbindlist(list(sRNA_data_PAR, XY_sRNA_data))

# graph of sRNA mapping by strata
# ptgs
XY_sRNA_data$strata2 <- factor(XY_sRNA_data$strata2, levels = c("Stratum 1", "Stratum 2", "Stratum 3", "PAR"))  
jpeg(f="Figures/sRNA_XY_pairs_ptgs.jpeg", width=8, height=8, units = "in", res=500, pointsize = 18, quality=100)
ggplot(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), data=XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='ptgs')]) + 
  theme_bw(base_size=20) +
  #  geom_jitter(aes(x=Chr_sex, y=corrected_TPM), size = 0.05, color = "grey") +
  geom_violin(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), show.legend = T, linewidth = 0, alpha=0.7) +
  geom_boxplot(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), position = position_dodge(width=0.9), outlier.shape = 1, outlier.size = 0.1, outlier.fill = "black", alpha=0.2, show.legend = F, notch=TRUE) +
  stat_summary(aes(x=Chr_sex, y=corrected_TPM, group=Chr_sex), position = position_dodge(width=0.9), geom = "point", fun.y = "mean", col = "red", size = 2, shape = 4, show.legend = F) +
  coord_cartesian(ylim =  c(0, 120)) +
  scale_fill_manual(values=c("brown2", "#006633", "#333399")) +
  labs(x="Chromosome and sex", y="sRNA mapping (TPM)", fill="Chr_sex") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  guides(fill="none") +
  facet_grid(rows = vars(tissue), cols = vars(strata2), scales = "free")
dev.off()
# rddm
jpeg(f="Figures/sRNA_XY_pairs_rddm.jpeg", width=8, height=8, units = "in", res=500, pointsize = 18, quality=100)
ggplot(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), data=XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='rddm')]) + 
  theme_bw(base_size=20) +
  #  geom_jitter(aes(x=Chr_sex, y=corrected_TPM), size = 0.05, color = "grey") +
  geom_violin(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), show.legend = T, linewidth = 0, alpha=0.7) +
  geom_boxplot(aes(x=Chr_sex, y=corrected_TPM, fill=Chr_sex), position = position_dodge(width=0.9), outlier.shape = 1, outlier.size = 0.1, outlier.fill = "black", alpha=0.2, show.legend = F, notch=TRUE) +
  stat_summary(aes(x=Chr_sex, y=corrected_TPM, group=Chr_sex), position = position_dodge(width=0.9), geom = "point", fun.y = "mean", col = "red", size = 2, shape = 4, show.legend = F) +
  coord_cartesian(ylim =  c(0, 250)) +
  scale_fill_manual(values=c("brown2", "#006633", "#333399")) +
  labs(x="Chromosome and sex", y="sRNA mapping (TPM)", fill="Chr_sex") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  guides(fill="none") +
  facet_grid(rows = vars(tissue), cols = vars(strata2), scales = "free")
dev.off()


# glms rddm flowers
# with zero inflation and log transformation log1p(x)=log(1+|x|)
model_sRNA_XYpairs = glmmTMB(log1p(corrected_TPM) ~ strata2 * Chr_sex + (1|GeneName), dispformula =~ 1, ziformula =~ 1, family=gaussian(), data = XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='rddm')&(tissue == "flowers")])
res = simulateResiduals(model_sRNA_XYpairs)
jpeg(f="glmer_output/residuals_sRNA_XY_genes_rddm_flowers.jpeg", width=13, height=8, units = "in", res=500, pointsize = 18, quality=100)
plot(res) 
dev.off()
summary(model_sRNA_XYpairs)
glmmTMB:::Anova.glmmTMB(model_sRNA_XYpairs)
contrasts <- emmeans(model_sRNA_XYpairs, list(pairwise ~ strata2 * Chr_sex))
write.csv(contrasts$`pairwise differences of strata2, Chr_sex`, 'glmer_output/contrasts_sRNA_XY_genes_rddm_flowers.txt')


# glms rddm leaves
# with zero inflation and log transformation log1p(x)=log(1+|x|)
model_sRNA_XYpairs = glmmTMB(log1p(corrected_TPM) ~ strata2 * Chr_sex + (1|GeneName), dispformula =~ 1, ziformula =~ 1, family=gaussian(), data = XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='rddm')&(tissue == "leaves")])
res = simulateResiduals(model_sRNA_XYpairs)
jpeg(f="glmer_output/residuals_sRNA_XY_genes_rddm_leaves.jpeg", width=13, height=8, units = "in", res=500, pointsize = 18, quality=100)
plot(res) 
dev.off()
summary(model_sRNA_XYpairs)
glmmTMB:::Anova.glmmTMB(model_sRNA_XYpairs)
contrasts <- emmeans(model_sRNA_XYpairs, list(pairwise ~ strata2 * Chr_sex))
write.csv(contrasts$`pairwise differences of strata2, Chr_sex`, 'glmer_output/contrasts_sRNA_XY_genes_rddm_leaves.txt')


# glms ptgs flowers
# with zero inflation and log transformation log1p(x)=log(1+|x|)
model_sRNA_XYpairs = glmmTMB(log1p(corrected_TPM) ~ strata2 * Chr_sex + (1|GeneName), dispformula =~ 1, ziformula =~ 1, family=gaussian(), data = XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='ptgs')&(tissue == "flowers")])
res = simulateResiduals(model_sRNA_XYpairs)
jpeg(f="glmer_output/residuals_sRNA_XY_genes_ptgs_flowers.jpeg", width=13, height=8, units = "in", res=500, pointsize = 18, quality=100)
plot(res) 
dev.off()
summary(model_sRNA_XYpairs)
glmmTMB:::Anova.glmmTMB(model_sRNA_XYpairs)
contrasts <- emmeans(model_sRNA_XYpairs, list(pairwise ~ strata2 * Chr_sex))
write.csv(contrasts$`pairwise differences of strata2, Chr_sex`, 'glmer_output/contrasts_sRNA_XY_genes_ptgs_flowers.txt')


# glms ptgs leaves
# with zero inflation and log transformation log1p(x)=log(1+|x|)
model_sRNA_XYpairs = glmmTMB(log1p(corrected_TPM) ~ strata2 * Chr_sex + (1|GeneName), dispformula =~ 1, ziformula =~ 1, family=gaussian(), data = XY_sRNA_data[(!is.na(strata2))&(strata2!='')&(path =='ptgs')&(tissue == "leaves")])
res = simulateResiduals(model_sRNA_XYpairs)
jpeg(f="glmer_output/residuals_sRNA_XY_genes_ptgs_leaves.jpeg", width=13, height=8, units = "in", res=500, pointsize = 18, quality=100)
plot(res) 
dev.off()
summary(model_sRNA_XYpairs)
glmmTMB:::Anova.glmmTMB(model_sRNA_XYpairs)
contrasts <- emmeans(model_sRNA_XYpairs, list(pairwise ~ strata2 * Chr_sex))
write.csv(contrasts$`pairwise differences of strata2, Chr_sex`, 'glmer_output/contrasts_sRNA_XY_genes_ptgs_leaves.txt')




