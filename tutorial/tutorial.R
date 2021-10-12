source("Documents/TestRepo/miniproject.R")


#raw counts data for data normalization and statistical calculations, it should include ensembl gene_id as the first column and count values of the samples in the other columns 
data <- read.csv("Documents/Alotaibi-RNAseq-Histone-Variant-09.09.2021/Rawdata/H2afv.counts.csv",
                 header = T,
                 row.names = 1)


#groups can consist sample replicas or different datasets with the same condition 
group1 <- c("H2afv1","H2afv2")
group2 <- c("Scr1","Scr2")


#data frames required for the functions first data frame comes from the normalization function and the second one create gene_id list from the first data frame
df <- DGE.analysis(data, group1 , group2, method)
df2 <- rownames_to_column(df$table, "ensembl_gene_id_version")

#dataset variable is needed for mart function and it should be given according to the organism of interest for the analysis
dataset <- "mmusculus_gene_ensembl"

#data frame for the mart function 
df3 <- mart.func(df2, dataset)

# number variable is an integer and if there is a desired limit for the most significantly meaningful differently regulated genes, 
# the number variable can be changed according to the requirement of the analysis. 
number <- Inf

# p value gives the statistical significance and if desired it can be lowered than the default 0.05,
p <- 0.05


# lastly organism of interest should be given for gene enrichment analysis
organism <- "org.Mm.eg.db"


# data frame of the final data with the information of gene enrichment
# dotplot visualization of the desired gene enrichment sets for 25 categories, dotplot can be altered according to gene enrichment results

df4 <- go.func(df,df2,df3,number,p,organism)

dotplot(df4,
        showCategory = 25,
        font.size = 10,
        title = "Dotplot of Enriched Pathways"
)
