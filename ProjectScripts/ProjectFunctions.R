packageF("biomaRt")
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"), mart = ensembl)


if(!"ermineR" %in% rownames(installed.packages())){
  install_github("PavlidisLab/ermineR", force = T)
}
library(ermineR)


GenericHumanAnno <-  GetAnnoFiles("Generic_human")


DESeq2RUN <- function(data, Meta, model){
  DESeqDS <- DESeqDataSetFromMatrix(countData = data, colData = Meta, design = model)
  DESeqOut <- DESeq(DESeqDS)
  return(DESeqOut)
}

GetDESeq2Results <- function(DESeqOut, coef, alpha = 0.05, indepFilter = TRUE){
  DEresults <- results(DESeqOut, name = coef, alpha = alpha, format = "DataFrame", independentFiltering = indepFilter)
  DEresults$GeneSymbol <- geneNames$hgnc_symbol[match(rownames(DEresults), geneNames$ensembl_gene_id)]
  DEresults$EnsemblID <- rownames(DEresults)
  DEresults %<>% data.frame %>% filter(GeneSymbol != "")
  return(DEresults)
}


GetOneSidedPval <- function(ResultsObj, adjust = "BH", logFCcol = "log2FoldChange", GeneCol = "GeneSymbol", pvalCol = "pvalue"){
  DESeqResultsDF <- data.frame(ResultsObj)
  #Just for now - remove the duplicated genes (5 at this point) 
  DESeqResultsDF <- DESeqResultsDF[!duplicated(DESeqResultsDF$GeneSymbol),]
  DESeqResultsDF$DownPval <- apply(DESeqResultsDF %>% select_(.dots = c(logFCcol, pvalCol)), 1, function(x){
    if(x[1] < 0){
      x[2]/2
    } else {
      1-x[2]/2
    }
  })
  DESeqResultsDF$DownPvalAdj <- p.adjust(DESeqResultsDF$DownPval, "BH")
  DESeqResultsDF$UpPval <- apply(DESeqResultsDF %>% select_(.dots = c(logFCcol, pvalCol)), 1, function(x){
    if(x[1] > 0){
      x[2]/2
    } else {
      1-x[2]/2
    }
  })
  DESeqResultsDF$UpPvalAdj <- p.adjust(DESeqResultsDF$UpPval, "BH")
  rownames(DESeqResultsDF) <- as.character(DESeqResultsDF[[GeneCol]])
  return(DESeqResultsDF)
}

CompareResults <- function(data1, data2, name1 = "NBB", name2 = "Norway", colorCol = "Cohort", Title){
  Signif_1 <- data1 %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
  Signif_2 <- data2 %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
  
  UnionSignif <- merge(data1 %>% data.frame %>% filter(GeneSymbol %in% c(Signif_1, Signif_2)),
                       data2 %>% data.frame %>% filter(GeneSymbol %in% c(Signif_1, Signif_2)), by = "GeneSymbol",
                       all.x=T, all.y = T, suffixes = c(paste0("_", name1),
                                                        paste0("_", name2)))
  
  UnionSignif$Cohort <- sapply(UnionSignif$GeneSymbol, function(x){
    x = as.character(x)
    if(x %in% Signif_1 & x %in% Signif_2){
      "Both"
    } else if(x %in% Signif_1){
      name1
    } else if(x %in% Signif_2){
      name2
    }
  }) %>% factor
  
  CompareCol <- grep("log2FoldChange", names(UnionSignif), value = T)
  ggplot(UnionSignif, aes_string(CompareCol[1], CompareCol[2], color = colorCol)) +
    theme_classic() +
    labs(title = Title) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
}

CompareResultsAll <- function(data1, data2, name1 = "NBB", name2 = "Norway", colorCol = "Cohort", Title){
  Signif_1 <- data1 %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
  Signif_2 <- data2 %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
  
  UnionSignif <- merge(data1 %>% data.frame,
                       data2 %>% data.frame, by = "GeneSymbol",
                       all.x=T, all.y = T, suffixes = c(paste0("_", name1),
                                                        paste0("_", name2)))
  
  UnionSignif$Cohort <- sapply(UnionSignif$GeneSymbol, function(x){
    x = as.character(x)
    if(x %in% Signif_1 & x %in% Signif_2){
      "Both"
    } else if(x %in% Signif_1){
      name1
    } else if(x %in% Signif_2){
      name2
    } else {
      "NS"
    }
  }) %>% factor
  
  CompareCol <- grep("log2FoldChange", names(UnionSignif), value = T)
  ggplot(UnionSignif, aes_string(CompareCol[1], CompareCol[2])) +
    theme_classic() +
    labs(title = Title) +
    geom_point(color = "grey", alpha = 0.5) +
    geom_point(data = UnionSignif %>% filter(Cohort != "NS") %>% droplevels(),
               aes_string(CompareCol[1], CompareCol[2],  color = "Cohort"),
               inherit.aes = FALSE) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_density2d()
}

GetAdjCountDESeq <- function(dds, Gene,  adjCov){
  GeneRow <- which(rownames(dds) == Gene)
  Mu <- log2(t(t(assays(dds)$mu[GeneRow,])/sizeFactors(dds)))
  Counts <-  log2(t(t(assays(dds)$counts[GeneRow,])/sizeFactors(dds)))
  Resid <- Counts - Mu
  corrFactor <- sizeFactors(dds)
  Coef <- coef(dds)[GeneRow,]
  mod <- attr(dds, "modelMatrix")
  modAdj <-mod
  for(cov in as.character(adjCov$Cov)){
    adjType = adjCov %>% filter(Cov == cov) %>% .$adjType %>% as.character()
    if(adjType == "mean"){
      modAdj[,cov] <- mean(modAdj[,cov], na.rm=T)
      
    } else if (adjType == "base"){
      modAdj[,cov] <- 0
    }
  }
  AdjValue <- (modAdj %*% Coef) + Resid
  return(AdjValue)
}
