######################## DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS ########################

################## NCBI DATASET (GSE224356) ##############################

########### Using raw counts from featurecounts for analysis ######################

############### Set up your workin directory ################################

#-------------------------------------------------------------------------------

#################### STEP-1: LOAD LIBRARIES & PACKAGES #########################

#-------------------------------------------------------------------------------


library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(org.Hs.eg.db)


#-------------------------------------------------------------------------------

######################### STEP-2: LOAD COUNTS & SAMPLE DATA ####################

#-------------------------------------------------------------------------------


############ Load counts data ##########


counts <- data.matrix(read.csv("counts/counts.csv", 
                     header = TRUE, 
                     sep = ",",
                     row.names = 1))#### DESEQ2 expects geneIDs as rownames


class(counts) ##### counts should be a matrix 

head(counts) ##### check the first few enteries of the matrix


######## Load the sample information table ###########

coldata <- read.csv("Sample_Condition.csv", row.names = 1)

head(coldata)

######## Check if counts & coldata are consistent in terms of sample order

###### The sample order is not same
###### Rearrange to make sample order consistent

all(rownames(coldata) %in% colnames(counts)) ### Check if all samples are present in counts & coldata

all(rownames(coldata) == colnames(counts)) ##### Check the sample order

counts <- counts[, rownames(coldata)] #### Rearrange the rows of counts as per coldata

all(rownames(coldata) == colnames(counts)) ### Check the sample order again

#### The rows of coldata matches the columns of counts

#-------------------------------------------------------------------------------

####################### STEP-3: CREATE DESeqDataSetFromMatrix #################

#-------------------------------------------------------------------------------

######## Define the factors from coldata

coldata$Condition <- factor(coldata$Condition,
                            levels = c("Tumor", "Normal"))

levels(coldata$Condition) #### Verify the factor levels 

######## The design is based on Condition (Tumor or Normal)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Condition)
dds

#-------------------------------------------------------------------------------

############################# STEP-4: PRE-FILTERING ############################

#-------------------------------------------------------------------------------

##### This step is optional

###### Saves computational resources by removing low counts genes

###### Keep the rows that have counts of at least 10 

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds

#-------------------------------------------------------------------------------

########################### STEP-5: RUN DESEQ2 FOR DE ANALYSIS #################

#-------------------------------------------------------------------------------

dds <- DESeq(dds)

res <- results(dds)

res

#-------------------------------------------------------------------------------

############################ STEP-6: COMPARE EXPRESSION ########################

#-------------------------------------------------------------------------------

########### Extract DE Results based on Condition (Tumor vs Normal)

res <- results(dds, contrast=c("Condition","Tumor","Normal")) 


######## Order the results table by adj_p_value

resOrdered <- res[order(res$pvalue),]

summary(res)

sum(res$padj < 0.1, na.rm=TRUE) ### how many adj_p_vales were less than 0.1


res05 <- results(dds, alpha=0.05) #### Select DE genes with a cutoff 0.05

summary(res05)

res_sig <- res05[which(res05$padj < 0.05), ] ##### setting adj_p_value < 0.05

summary(res_sig)

#-------------------------------------------------------------------------------

############################ STEP-7: PLOTTING THE RESULTS ######################

#-------------------------------------------------------------------------------


######## MA Plot

plotMA(res)


############# Volcano Plot

EnhancedVolcano(res_sig, lab=rownames(res_sig), x='log2FoldChange', y='pvalue')

######### Convert the geneIDs to gene_names for better understanding


gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(res_sig),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")


res_sig$gene_name <- gene_symbols

view(res_sig)

############ Volcano Plot 

EnhancedVolcano(res_sig, lab=res_sig$gene_name, x='log2FoldChange', y='pvalue')

############ Upregulated & Downregulated genes 

up <- res_sig[res_sig$log2FoldChange > 1, ]

nrow(up)

down <- res_sig[res_sig$log2FoldChange < -1, ]

nrow(down)

######## Save the results into a CSV

write.csv(as.data.frame(res_sig), file = "Results/DE_genes.csv")
write.csv(as.data.frame(up), file = "Results/Upregulated_genes.csv")
write.csv(as.data.frame(down), file = "Results/Downregulated_genes.csv")


