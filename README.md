# Differential Expression Analysis Using DESEQ2

The goal of this project is to identify differentially expressed genes in papillary thyroid carcinoma (PTC). For this purpose RNA-Seq data (raw FASTQ reads) was downloaded from NCBI GEO with accession [GSE224356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224356). The dataset consisted of 6 samples; 3 carcinoma tissue and 3 paired normal tissue.

## Analysis Workflow

1.  **Data Import**

-   Raw FASTQ reads were obtained for all samples

2.  **Quality Check**

-   Tool: fastqc
-   Quality of raw reads was checked.

3.  **Quality Control**

-   Tool: fastp
-   All-in-one pre-processing was performed, including adapter trimming and removal of low quality bases.

4.  **Alignment**

-   Tool: HISAT2
-   The pre-processed reads were aligned with reference genome (hg38) using default parameters and BAM files were obtained.

5.  **Post-alignment Processing**

-   Tool: Markduplicates & Rmdup
-   The PCR duplicates were first marked and removed from BAM files.

6.  **Transcript Quantification**

-   Tool: featurecounts
-   The aligned reads were assembled into transcripts and quantified. The final counts matrix was used to perform DE analysis.

7.  **DE Analysis**

-   Tool: DESEQ2
-   Counts matrix was imported into RStudio, DESeqDataSet object was created, low count genes were removed and final DE analysis was done. The DE genes (padj \< 0.05 and \|log2FC\| ≥ 1) lists were obtained and upregulated and downregulated genes were identified.

8.  **Visualization**

-   Tools: MA Plot, Volcano Plot.
-   MA and enhanced volcano plots were used to view DE genes.

## Results

A total of 2752 genes were known to be DE. Among these, 1466 genes were upregulated and 903 genes were downregulated. Many of these DE genes are known to be involved in pathogenesis of PTC and act as biomarkers of the disease. Some of these known PTC-associated genes include S100A6, COL1A1, DHRS3, COL3A1, ZAP70, TIMP1, and SERPINA1. These genes contribute towards PTC and are known to be diagnostic and therapeutic biomarkers of the disease.

## References

-   [Identification of Differentially Expressed Genes in Papillary Thyroid Cancers](https://pmc.ncbi.nlm.nih.gov/articles/PMC2649849/)
-   [Identification of Biomarkers Based on Differentially Expressed Genes in Papillary Thyroid Carcinoma](https://www.nature.com/articles/s41598-018-28299-9)
-   [Network Analyses of Integrated Differentially Expressed Genes in Papillary Thyroid Carcinoma to Identify Characteristic Genes](https://www.mdpi.com/2073-4425/10/1/45)
