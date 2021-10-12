# README

This script consists of 3 functions: **Normalization**, **Biomart query creation** and **Gene enrichment analysis**

## Normalization
Data manipulation function requires a raw counts data with 2 groups to work with, these groups can be different conditions with different number of samples. The default normalization method is 'TMM', however it can be changed as a parameter as well to TMMwsp, RLE or upperquartile as a variable. The normalization and data manipulation function gives a quasilike log lineer model for statistical result table as a result to detect differences between conditions.

*Before starting the second function the ensembl_gene_id of the result table must be created as a data frame since it will be required in the gene enrichment analysis.*

## Biomart Query
Biomart function searches and lists the gene IDs of significantly meaningful regulated genes from the statistical results.It draws the complete information of the desired organism dataset from Ensembl, like *Homo sapiens* or *Mus musculus* and then aligns this information with the normalized counts data from the first function.

## Gene Enrichment Analysis
This function takes the normalized counts data and biomart results and divides the differently regulated genes between groups, and performs gene enrichment of these genes. The integer number is infiniteby default but it can be restricted to any desired number as a variable, and the p value can be used as default which is 0.05 but if test requires results with higher statistical significance it can be changes as a variable.

*After completing the gene enrichment analysis the results can be written as a table, but if desired a dotplot can be drawn to visualize the results for the GO terms. Also the result enrichment file can be used in different analysis.*
