# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
if (!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")
BiocManager::install("GEOquery")
if (!requireNamespace("limma", quietly = TRUE))
  install.packages("limma")
BiocManager::install("limma")
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
BiocManager::install("pheatmap")

library(BiocManager)
library(forcats)
library(ggrepel)
library(stringr)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(survminer)
library(GEOquery)
library(BiocGenerics)
library(limma)
library(umap)
library(parallel)
library(stats)
library(base)
library(org.Hs.eg.db)
rm(list = ls())
## change my_id to be the dataset that you want.
setwd("C:/Users/ESPIRAL/Desktop/FCFM-2021/verano2021/Nueva carpeta (4)")
my_id <- "GSE33126"
gse <- getGEO(my_id)
length(gse)
gse <- gse[[1]]
gse
## if more than one dataset is present, you can analyse the other dataset by changing the number inside the [[...]]
## e.g. gse <- gse[[2]]
pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data
pData(gse)$data_processing[1]
## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

library(dplyr)

sampleInfo <- pData(gse)
sampleInfo

## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns

sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)
sampleInfo

library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 

## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)

## If not, force the rownames to match the columns
#rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
		annotation_col=sampleInfo)

library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient)))+geom_point()+geom_text_repel()

		### CODE ONLY FOR DEMONSTRATION ONLY

		### lets' say are outliers are samples 1,2 and 3
		## replace 1,2,3 with the outliers in your dataset
		outlier_samples <- c(1,2,3)

		gse <- gse[,-outlier_samples]

		#exportando los datos
		library(readr)
		full_output <- cbind(fData(gse),exprs(gse))
		write_csv(full_output, path="gse_full_output.csv")
		features <- fData(gse)
		View(features)
		### Look at the features data frame and decide the names of the columns you want to keep
		features <- select(features,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
		full_output <- cbind(features,exprs(gse))
		write_csv(full_output, path="gse_full_output.csv")

#expresión diferencial
library(limma)
sampleInfo$group
dise<- model.matrix(~0 + sampleInfo$group)
dise
## the column names are a bit ugly, so we will rename
colnames(dise) <- c("Normal","Tumour")
summary(exprs(gse))

gse

## calculate median expression level
cutoff <- median(exprs(gse))
cutoff
## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff
dim(is_expressed)
## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

keep
## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]
dim(gse)
dim(dise)
	#pesos
	#aw<-arrayWeights(exprs(gse), design)
	#aw
	#fit <- lmFit(exprs(gse), dise, weigths=aw)
fit <- lmFit(exprs(gse), dise)
head(fit$coefficients)

contrasts <- makeContrasts(Tumour - Normal, levels=dise)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2)

topTable(fit2)
topTable(fit2, coef=1)
### to see the results of the second contrast (if it exists)
## topTable(fit2, coef=2)
decideTests(fit2)
table(decideTests(fit2))
## calculate relative array weights
aw <- arrayWeights(exprs(gse),dise)
aw

fit <- lmFit(exprs(gse), dise,weights = aw)
contrasts <- makeContrasts(Tumour - Normal, levels=dise)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

anno <- fData(gse)
anno

anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")

## Get the results for particular gene of interest
filter(full_results, Symbol == "SMOX")

## Get results for genes with TP53 in the name
filter(full_results, grepl("TP53", Symbol))

## Get results for one chromosome
filter(full_results, Chromosome==20)
p_cutoff <- 0.05
fc_cutoff <- 1

filter(full_results, adj.P.Val < 0.05, abs(logFC) > 1)

library(readr)
filter(full_results, adj.P.Val < 0.05, abs(logFC) > 1) %>%
  write_csv(path="filtered_de_results.csv")


### example code
library(pheatmap)
pheatmap(gene_matrix)

## Use to top 20 genes for illustration

topN <- 20
##
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)
gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Symbol) 
## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix, labels_row = gene_names)

pheatmap(gene_matrix, labels_row = gene_names, scale="row")

my_genes <- c("HIG2", "CA1","ETV4","FOXA1")
ids_of_interest <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(ID)

gene_names <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(Symbol)
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix,labels_row = gene_names,scale="row")

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

my_genes <- c("HIG2", "CA1","ETV4","FOXA1")

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              columns=c("ENSEMBL","GO"),
                              keys=my_genes,
                              keytype = "SYMBOL")

anno

anno <- AnnotationDbi::select(org.Hs.eg.db,
                              columns="SYMBOL",
                              keys="GO:0006338",
                              keytype="GO")

my_genes <- pull(anno, SYMBOL)
ids_of_interest <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(ID)

gene_names <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(Symbol)

gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row", cex=0.8)


#survival analysis
rm(list = ls())
library(GEOquery)
gse <- getGEO('GSE7390')[[1]]
View(pData(gse))
library(dplyr)
s_data <- pData(gse) %>% 
  dplyr::select(geo_accession, contains("characteristics"), -characteristics_ch1.2, -characteristics_ch1, -characteristics_ch1.1)

s_data <- s_data %>% 
  dplyr::rename(hospital = characteristics_ch1.3,
         age = characteristics_ch1.4,
         size = characteristics_ch1.5,
         surgery_type = characteristics_ch1.6,
         histtype = characteristics_ch1.7,
         angioinv = characteristics_ch1.8,
         lymp_infil = characteristics_ch1.9,
         node = characteristics_ch1.10,
         grade = characteristics_ch1.11,
         er = characteristics_ch1.12,
         t.rfs = characteristics_ch1.13,
         e.rfs = characteristics_ch1.14,
         t.os = characteristics_ch1.15,
         e.os = characteristics_ch1.16,
         t.dmfs = characteristics_ch1.17,
         e.dmfs =characteristics_ch1.18,
         t.tdm = characteristics_ch1.19,
         e.tdm = characteristics_ch1.20,
         risksg = characteristics_ch1.21,
         npi = characteristics_ch1.22,
         risknpi = characteristics_ch1.23,
         aol_os_10yr = characteristics_ch1.24,
         risk_aol = characteristics_ch1.25,
         veridex_risk = characteristics_ch1.26)

s_data <- s_data %>% 
  mutate(hospital = gsub("hospital: ", "", hospital, fixed=TRUE),
         age = as.numeric(gsub("age: ","", age, fixed=TRUE)),
         size = as.numeric(gsub("size: ","", size, fixed=TRUE)),
         surgery_type = gsub("Surgery_type: ","", surgery_type, fixed=TRUE),
         histtype = gsub("Histtype: ","", histtype, fixed=TRUE),
         angioinv = gsub("Angioinv: ","", angioinv, fixed=TRUE),
         lymp_infil = gsub("Lymp_infil: ","", lymp_infil, fixed=TRUE),
         node = gsub("node: ","", node, fixed=TRUE),
         grade = gsub("grade: ","", grade, fixed=TRUE),
         er = gsub("er: ","", er, fixed=TRUE),
         t.rfs = as.numeric(gsub("t.rfs: ","", t.rfs, fixed=TRUE)),
         e.rfs = as.numeric(gsub("e.rfs: ","", e.rfs, fixed=TRUE)),
         t.os = as.numeric(gsub("t.os: ","", t.os, fixed=TRUE)),
         e.os = as.numeric(gsub("e.os: ","", e.os, fixed=TRUE)),
         t.dmfs = as.numeric(gsub("t.dmfs: ","", t.dmfs, fixed=TRUE)),
         e.dmfs = as.numeric(gsub("e.dmfs: ","", e.dmfs, fixed=TRUE)),
         t.tdm = as.numeric(gsub("t.tdm: ","", t.tdm, fixed=TRUE)),
         e.tdm = as.numeric(gsub("e.tdm: ","", e.tdm, fixed=TRUE)),
         risksg = gsub("risksg: ","", risksg, fixed=TRUE),
         npi = as.numeric(gsub("NPI: ","", npi, fixed=TRUE)),
         risknpi = gsub("risknpi: ","", risknpi, fixed=TRUE),
         aol_os_10yr = as.numeric(gsub("AOL_os_10y: ","", aol_os_10yr, fixed=TRUE)),
         risk_aol = gsub("risk_AOL: ","", risk_aol, fixed=TRUE),
         veridex_risk = gsub("veridex_risk: ","", veridex_risk, fixed=TRUE)
         )
library(survival)
library(survminer)

fit <- survfit(Surv(t.os, e.os) ~ er, data = s_data)
ggsurvplot(fit, data = s_data,pval = TRUE)

exprs(gse)[1:5,1:5]
features <- fData(gse)
View(features)
features <- dplyr::select(features, ID, `Gene Symbol`)
features <- tidyr::separate(features,`Gene Symbol`, into = c("Symbol_1","Symbol_2"), sep=" /// ")

e_data <- exprs(gse) %>% 
  data.frame %>% 
  tibble::rownames_to_column("ID") %>% 
  left_join(features)

e_data

e_data <- filter(e_data,Symbol_1 %in% "ESR1") %>% 
  dplyr::select(-Symbol_1, -Symbol_2)
e_data

dim(s_data)
dim(e_data)

e_data <- e_data %>% tidyr::gather(geo_accession, Expression,-ID) %>% 
  tidyr::spread(ID, Expression)
e_data
## TO-DO: replace with pivot_wider and pivot_longer at some point

s_data <- left_join(s_data, e_data)
ggplot(s_data, aes(x = er, y =`205225_at`)) + geom_boxplot()
s_data <- mutate(s_data, Group = ifelse(`205225_at` > quantile(`205225_at`,0.25), "High","Low"))
fit <- survfit(Surv(t.os, e.os) ~ Group, data = s_data)
ggsurvplot(fit, data = s_data,pval = TRUE)

cutoff <- 0.5
s_data <- mutate(s_data, Group = ifelse(`205225_at` > quantile(`205225_at`,cutoff), "High","Low"))

fit <- survfit(Surv(t.os, e.os) ~ Group, data = s_data)
ggsurvplot(fit, data = s_data,pval = TRUE)













# load series and platform data from GEO
gset <- getGEO("GSE163530", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "undefined"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
LogC <- (qx[5] > 100) ||(qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
	ex[which(ex <= 0)] <- NaN
  	exprs(gset) <- log2(ex)
	}

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("undefined"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE163530", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE163530", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE163530")