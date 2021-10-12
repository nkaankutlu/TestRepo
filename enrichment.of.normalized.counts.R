library("edgeR")
library("biomaRt")
library("clusterProfiler")
library("org.Mm.eg.db")
library("tibble")
library("dplyr")

DGE.analysis <- function(rawdata, group1, group2){
  
  df <- data
  
  ####creating a DGElist for manipulation and normalization of the count values data####

  group <- factor(c(group1, group2))
  y <- DGEList(as.matrix(data), group = group)
  keep <- rowSums(cpm(y) > 1) >= 2
  dgeList <- y[keep, , keep.lib.size=FALSE]
  dgeList <- calcNormFactors(dgeList, method="TMM")
  dgeList$samples
  levels(dgeList$samples$group)
  
  ####performing GLM calculations####
  design <- model.matrix(~group, data = dgeList$samples)
  colnames(design) <- levels(dgeList$samples$group)
  dgeList <- estimateDisp.DGEList(dgeList, design)
  fit <- glmQLFit(dgeList, design)
  qlf <- glmQLFTest(fit, coef = 2)
  df.qlf <-rownames_to_column(qlf$table, "ensembl_gene_id_version")
  return(qlf)
}

mart.func <- function(df2, dataset){
  ####finding the gene IDs of significantly meaningful regulated genes####
  listMarts()
  ensembl <- useMart("ensembl", dataset= dataset, host = "asia.ensembl.org")
  datasets <- listDatasets(mart= ensembl)
  
  ensembl = useDataset(dataset, mart= ensembl)
  
  
  mart <- useEnsembl(biomart = "ensembl",
                     dataset= dataset,
                     mirror = "asia")
  
  BM <- getBM(attributes= c('ensembl_gene_id_version', 'external_gene_name'), 
              filters = "ensembl_gene_id_version", 
              values = df2[,c(1)], 
              mart = mart)
   return(BM)
}

go.func <- function(df, df2, df3, number, p, organism){
  
  merge <- merge(df2, df3, by="ensembl_gene_id_version", all=TRUE)
  
  #### gene set enrichment analysis of differently regulated genes####
  
  topTags <- topTags(df, n=number, adjust.method="fdr", p.value = p, sort.by = "logFC")
  df.topTags <- as.data.frame(topTags)
  df.topTags <-rownames_to_column(df.topTags, "ensembl_gene_id_version")
  mergetoptags <- merge(df.topTags, df3, by= "ensembl_gene_id_version", all.x = TRUE)
  genenames <- mergetoptags$external_gene_name
  
  
  GOenrich <- enrichGO(gene = genenames,
                       OrgDb = organism,
                       ont = "BP",
                       keyType = "SYMBOL"
  )
  
 
  return(GOenrich)
}
