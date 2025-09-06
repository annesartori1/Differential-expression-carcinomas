#Anne Sartori
#BS 831
#Final Project

setwd("C:\\Users\\annes\\Dropbox (Personal)\\Personal\\Biostats\\BS831 Genomics\\Final Project")
getwd()

#rm(list = ls()) #Clear Global Environment

library(dplyr)
library(ggplot2)
library(Biobase)
library(DESeq2)

# load raw data and explore it. exprs data 34422 genes by 120 samples
#pData 120 samples by 25 variables
HNC <- readRDS("HNSC_htseq_raw_counts_AEvsG1vsG3.rds")
dim(HNC)
head(exprs(HNC))
str(exprs(HNC))
head(pData(HNC))
dim(pData(HNC)) #has data on smoking and HPV status and alcohol consn 120x 25
head(fData(HNC)) #gene symbols are here in hgnc_symbol
names(pData(HNC))
summary(pData(HNC))
table(pData(HNC)$patient.number_pack_years_smoked, useNA="ifany")
table(pData(HNC)$patient.tumor_tissue_site) #all head and neck
table(pData(HNC)$patient.frequency_of_alcohol_consumption, useNA="ifany") #about 47 nonmissing
table(pData(HNC)$grade, useNA="ifany")#There are no missing, but I want to remove AE
table(pData(HNC)$grade, pData(HNC)$patient.frequency_of_alcohol_consumption, useNA=
        "ifany")
table(pData(HNC)$patient.hpv_status_by_ish_testing, pData(HNC)$patient.hpv_status_by_p16_testing,useNA="ifany")
#We don't seem to really have the HPV dtaa. There 7 samples with data for either HPV variable, and only 1 positive
#106/120 that have missing on both HPV variables--so they are basing the claims on 14 samples or have more data than these
table(pData(HNC)$grade, useNA="ifany")

###Preprocessing####

#Look at the expression data
hist(exprs(HNC))
exprs(HNC)[1:5,1:5]

# Subset HNC to exclude samples where grade == "AE"
HNC_filtered <- HNC[, pData(HNC)$grade != "AE"]
pData(HNC_filtered)$grade <- droplevels(pData(HNC_filtered)$grade)
dim(exprs(HNC_filtered)) #34422 by 80: this worked.
table(HNC_filtered$grade)
hist(exprs(HNC)) #unsurprisingly, visually the same

# Remove genes with zero counts
genes_zeros <- rowSums(exprs(HNC_filtered)) == 0
table(genes_zeros)

HNC_filtered <- HNC_filtered[!genes_zeros, ]
dim(HNC_filtered) #32208 x 80
hist(exprs(HNC_filtered))
exprs(HNC_filtered)[1:5,50:80]#There are a couple in the 4th row in later cols that are not zero, so the zeros in the first 5 cols don't indicate a problem with filtering

#Calculate MADs
mads <- apply(exprs(HNC_filtered), 1, mad)
length(mads)

# order MAD values
mads_sort <- sort(mads, decreasing = TRUE)
head(mads_sort)

# subset top 5000
mads_sub <- mads_sort[1:5000]
head(mads_sub)
length(mads_sub)

# subset data object
# using 'names()' to select names of genes
HNC_filtered <- HNC_filtered[names(mads_sub), ]
dim(HNC_filtered)

#SAVE a version of HNC_filtered for clustering and classification before filtering more for differential analysis
cdata <- HNC_filtered 

#### Differential expression analyses using DESeq2 ###


# Format data for DESeq2
#Remove data with missing on covariates
# Logical vector: TRUE for samples where both smoking and alcohol are not NA
keep2 <- !is.na(pData(HNC_filtered)$patient.frequency_of_alcohol_consumption) & !is.na(pData(HNC_filtered)$patient.number_pack_years_smoked)


# Subset the ExpressionSet to keep only those samples
HNC_filtered <- HNC_filtered[, keep2]
dim(HNC_filtered) #Now I have only 18 samples 
cols_keep <- c(
  "patient.frequency_of_alcohol_consumption", 
  "patient.number_pack_years_smoked", 
   "grade"
)

col_data <- as.data.frame(pData(HNC_filtered)[, cols_keep])
col_data$grade <- factor(col_data$grade, levels = c("g1", "g3")) #g1 is "control," or reference level

# define formula
# where phenotype variable of interest is last, so DESeq2::results() extracts it by default
fm <- "~ patient.frequency_of_alcohol_consumption + patient.number_pack_years_smoked  + grade"

# remove any rows with NAs from col_data (required by DESeq2)
nas_remove <- apply(col_data, 1, anyNA)
table(nas_remove)

exprs_deseq <- exprs(HNC_filtered)[, !nas_remove]
col_data <- col_data[!nas_remove, ]
dim(exprs_deseq) #This did not remove any more
dim(col_data)


dds <- DESeqDataSetFromMatrix(
  countData = exprs_deseq, 
  colData = col_data, 
  design = formula(fm)
)

# run DESeq2
dds <- DESeq(dds)

# extract results
res_deseq <- results(dds)
head(res_deseq)

# number of significant genes
#(I think there was only 1 significant gene with smoking*alcohol interaction--come back to this)
table(res_deseq$padj < 0.05)
#top significant genes
res_deseq_sort <- res_deseq[order(res_deseq$padj), ]
head(res_deseq_sort)
#gene symbols
fData(HNC_filtered)[rownames(res_deseq_sort)[1:3], ]
# MA plot
DESeq2::plotMA(res_deseq, alpha = 0.05)

#ggplot  MA plot of the results. Significance
#value matches the threshold above: adjusted p-value is <0.05.

library(ggplot2)

# Convert results to data frame
res_df <- as.data.frame(res_deseq)
res_df$gene <- rownames(res_deseq)
str(res_df)
dim(res_df)

# Define significance

res_df$significant <- ifelse(
  is.na(res_df$padj), 
  "NA", 
  ifelse(res_df$padj < 0.05, "Significant", "Not significant")
)

# Plot
ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_x_continuous(trans = "log10") + 
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean Expression (baseMean, log scale)",
    y = "Log2 Fold Change"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")

### Gene set enrichment analysis ####
remotes::install_github("lmweber/hypeR")
BiocManager::install("msigdbr")
install.packages("msigdbr")
library(hypeR)
library(msigdbr)
?msigdbr::msigdbr


#We connect to the database on the internet and download the database we want
hallmark <- msigdb_gsets(species = "Homo sapiens", db_species = "HS", collection = "H")

# check object
hallmark
str(hallmark)
class(hallmark)
?gsets

# gene sets
str(hallmark$genesets)
length(hallmark$genesets)  # contains 50 gene sets





#########Alt version below
HALLMARK <- msigdb::msigdb_gsets(
  species = "Homo sapiens",
  db_species = "HS",
  collection = "H"
)
library(msigdbr)

# Get Hallmark gene sets for Homo sapiens
HALLMARK <- msigdbr(species = "Homo sapiens", category = "H")
HALLMARK
str(HALLMARK)
class(HALLMARK)

#Get descriptions of the gene sets
descriptions <- unique(HALLMARK[, c("gs_name", "gs_description")])
print(descriptions)

#Try this first with Epithelial
# extract gene set of interest
gs_emt <- subset(HALLMARK, gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

# check output
str(gs_emt)
head(gs_emt)
length(gs_emt)

dim(HNC_filtered)
# extract gene symbol list from expression set
head(fData(HNC_filtered))
genes_eset <- fData(HNC_filtered)$hgnc_symbol
head(genes_eset)

# convert to upper case to avoid mismatches due to upper / lower case
gs_upper <- toupper(gs_emt$gene_symbol)
genes_eset_upper <- toupper(genes_eset)

head(gs_upper)
head(genes_eset_upper)

# calculate overlap
# 141
sum(gs_upper %in% genes_eset_upper)
sum(genes_eset_upper %in% gs_upper)

# calculate indices of matching genes in the two gene lists
# after converting gene symbols to upper case to avoid mismatches due to upper / lower case
gene_idx <- match(gs_upper, genes_eset_upper)

gene_idx
length(gene_idx)
head(genes_eset_upper[gene_idx]) 
# remove any missing / NAs
gene_idx <- gene_idx[!is.na(gene_idx)]

gene_idx
length(gene_idx) 

# Remove NAs and rank by stat, which are Wald statistics, follow a standard Normal under the null
res <- res_deseq[!is.na(res_deseq$stat), ]
Wald_stats <- res$stat
names(Wald_stats) <- rownames(res)

# Sort in decreasing order (most upregulated first)
gene_stats <- sort(gene_stats, decreasing = TRUE)
head(gene_stats)
length(gene_stats)
# Use in ksGenescore
# look at the ranks of the Wald statistics

# Rank the genes (higher stat = higher rank)
gene_ranks <- rank(Wald_stats)
# show first few ranks
head(sort(gene_ranks))
length(gene_ranks)         # Should be the same as gene_idx
summary(gene_ranks)        # Should have integers between 1 and 5000
any(is.na(gene_ranks))     # Should be FALSE
gene_ranks_vec <- as.numeric(gene_ranks)
matched_ranks <- gene_ranks_vec[gene_idx]
# load ksGenescore() function: Calculates the K-S gene scores
source("ksGenescore.R")
# calculate K-S gene scores and generate plot using ksGenescore()
ksGenescore(n.x = 5000, y = matched_ranks, do.plot = TRUE)


#Now prepare the data to look at all the genesets.
# Assuming res.df is the DESeq2 result dataframe
# Assuming fData(HNC_filtered)$hgnc_symbol contains gene symbols

#Convert rownames to column in res.df and get gene symbols from fData
# res_df has a column named 'gene' with gene IDs
# and 'fData(HNC_filtered)$hgnc_symbol' contains gene symbols

# Perform the inner join by the 'gene' column
res_with_symbols <- dplyr::inner_join(
  fData(HNC_filtered) |> 
    tibble::rownames_to_column("gene") |>  # Convert rownames to "gene" column (assuming rownames are gene IDs)
    dplyr::select(gene, hgnc_symbol),  # Keep only gene_id and hgnc_symbol
  
  res_df,  # DESeq2 results (which already has 'gene' column)
  by = "gene"  # Join by the 'gene' column
)
str(res_with_symbols)
# Now 'res_with_symbols' contains both DESeq2 results and gene symbols
#Professor Weber has code to extract ranks by symbols. I did that above for the top few.
str(HALLMARK)
##Enrichment for whole Hallmark set of sets
ranked_signature <- res_with_symbols |> 
  dplyr::arrange(stat) |> 
  dplyr::select(hgnc_symbol, stat) |> 
  tibble::deframe()

#These are Wald stats
head(ranked_signature)  # most negative: up-regulated in G3 (second category)
tail(ranked_signature)  # most positive: up-regulated in G1 (first category)


ks_obj <- hypeR::hypeR(
  signature = ranked_signature, 
  genesets = hallmark, 
  test = "kstest",
  background = nrow(HNC) #how many tests are we starting with: #rows before filtering
)

# look at results
str(ks_obj)
class(ks_obj) #info stored in ks_obj
?hyp
str(ks_obj$data)
ks_obj$data[1:10, ]
ks_results <- ks_obj$data
significant_sets <- ks_results[ks_results$fdr < 0.05, ]

# View the top significant sets
head(significant_sets)
# Sort by FDR
significant_sets <- significant_sets[order(significant_sets$fdr), ]

# Save to CSV
write.csv(significant_sets, "significant_enriched_sets.csv", row.names = FALSE)


# dot plots (enrichment dot plot)
hypeR::hyp_dots(ks_obj, top = 20) #top 20 of the hallmark lists ranked in significance by FDR
#the one we were interested in at the start is significant, but not the most significant

##### CLUSTERING (and after that, Classification) #####
library(dplyr)
library(ggplot2)
library(matrixStats)
library(caret)
library(cba)
library(naivebayes)
library(mclust)
library(glmnet)
# from Bioconductor
library(Biobase)
library(ComplexHeatmap)

# load helper functions
source("hcopt.R")
source("misc.R")

#Note that I use the *full group of g1 and g3 samples* for the clustering,
#not just the 18 samples that were left after I removed missing on my covariates
#for DESeq2 earlier

#This filtered full group is in cdata.

# set up heatmap annotation
# color formatted as a list
table(cdata$grade)

dds_full <- DESeqDataSetFromMatrix(
  countData = exprs(cdata), 
  colData = pData(cdata),
  design = ~ grade
)


#Apply variance-stabilizing transformation to output of DESeq2 analysis
vst_data <- vst(dds_full, blind = TRUE)
#Extract transformed expression matrix
vst_mat <- assay(vst_data)

all(colnames(vst_mat) == rownames(pData(cdata)))

hm_ann2 <- HeatmapAnnotation(
  grade = pData(cdata)$grade,
  col = list(grade = c("g1" = "pink", "g3" = "lightgray"))
)

hc_col <- hcopt(dist(t(exprs(cdata))), method = "ward.D")
hc_row <- hcopt(as.dist(1 - cor(t(exprs(cdata)))), method = "ward.D")

Heatmap(
  vst_mat, 
  name = "expression", 
  # col = BS831::colGradient(c("blue", "white", "darkred"), 9),
  top_annotation = hm_ann2, 
  cluster_rows = hc_row,  # can also set to TRUE for default clustering
  cluster_columns = hc_col,  # can also set to TRUE for default clustering
  column_split = 2, 
  row_split = 3, 
  show_parent_dend_line = TRUE, 
  row_title = "", 
  show_column_names = FALSE, 
  show_row_names = FALSE
)

#Classification

#Fit an ElasticNet model using glmnet
#install.packages("glmnet")
library(caret)
library(caret)
library(glmnet)
library(Biobase)

# Prepare expression data
exprs_df <- as.data.frame(t(exprs(cdata)))
gene_class <- factor(pData(cdata)$grade, levels = c("g3", "g1"))  # g3 as positive class
exprs_df$Class <- gene_class

# Split data into training and test sets
set.seed(42)
train_index <- createDataPartition(exprs_df$Class, p = 0.75, list = FALSE)
train_data <- exprs_df[train_index, ]
test_data  <- exprs_df[-train_index, ]

# Train control
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Grid of alpha (mix of L1/L2) and lambda (penalty strength)
tune_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.25),      # Ridge to Lasso
  lambda = 10^seq(-4, 1, length = 20)
)

# Train with glmnet
enet_model <- train(
  Class ~ .,
  data = train_data,
  method = "glmnet",
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = tune_grid,
  preProcess = c("center", "scale")
)

# View results
enet_model$bestTune
plot(enet_model)

## Extract final glmnet model
final_model <- enet_model$finalModel

# Coefficients at best alpha/lambda
coefs <- coef(final_model, s = enet_model$bestTune$lambda)

# Convert to a readable data frame
selected_genes <- as.matrix(coefs)
# Drop the intercept row
selected_genes <- selected_genes[rownames(selected_genes) != "(Intercept)", , drop = FALSE]
# Filter out any zeroes
selected_genes <- selected_genes[selected_genes[,1] != 0, , drop = FALSE]

# View selected gene names and coefficients
print(selected_genes)


## Predictions on training data
train_predictions <- predict(enet_model, newdata = train_data)

# Confusion matrix for training data
train_confusion <- confusionMatrix(train_predictions, train_data$Class)

# Print training accuracy
train_accuracy <- train_confusion$overall['Accuracy']
print(train_accuracy)

# Evaluate on held-out test set
pred_probs <- predict(enet_model, newdata = test_data, type = "prob")
pred_class <- predict(enet_model, newdata = test_data)

confusionMatrix(pred_class, test_data$Class)
roc_curve <- pROC::roc(test_data$Class, pred_probs[, "g3"])
plot(roc_curve)


# Extract gene indices corresponding to non-zero coefficients
selected_gene_indices <- rownames(selected_genes)
head(selected_gene_indices)
length(selected_gene_indices)#10 because I took out intercept
head(fData(cdata)$hgnc_symbol)

# Get gene symbols 
str(fData(cdata))

# Loop through selected gene IDs and print matching HGNC symbols
for (gene_id in selected_gene_ids) {
  match_index <- which(fData(cdata)$ensembl_gene_id == gene_id)
  if (length(match_index) > 0) {
    symbol <- fData(cdata)$hgnc_symbol[match_index]
    cat(gene_id, ":", symbol, "\n")
  } else {
    cat(gene_id, ": not found in fData(cdata)$ensembl_gene_id\n")
  }
}


#fData(HNC)[fData(HNC)$ == "ENSG00000166819", "hgnc_symbol"]
#fData(HNC)[fData(HNC)$ensembl_gene_id == "ENSG00000167588", "hgnc_symbol"]
#fData(HNC)[fData(HNC)$ensembl_gene_id == "ENSG00000184530", "hgnc_symbol"]


