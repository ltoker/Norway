packageF("DESeq2")
if(!"ermineR" %in% rownames(installed.packages())){
  install_github("oganm/ermineR", force = T)
}
library(ermineR)

if(length(list.files(pattern = "Generic_human_noParents.an.txt")) == 0){
  download.file("http://chibi.ubc.ca/microannots/Generic_human_noParents.an.txt.gz", destfile="Generic_human_noParents.an.txt")
} else {
  warning("using existing GPL file")
}

GenericHumanAnno <-  read.table("Generic_human_noParents.an.txt", header = TRUE, sep = "\t", quote = "")
studyFinal$Cortex$Metadata$age_years <- as.numeric(studyFinal$Cortex$Metadata$age_years)

MetaCovar = "sex + age_years + pm_time_min + rin + Batch"
Model = as.formula(paste0("~Profile + cohort +", MetaCovar))

DESeq2RUN <- function(data, Meta, model){
  DESeqDS <- DESeqDataSetFromMatrix(countData = data, colData = Meta, design = model)
  DESeqDS = estimateSizeFactors(DESeqDS)
  DESeqOut <- DESeq(DESeqDS)
  return(DESeqOut)
}

GetDESeq2Results <- function(DESeqOut, coef, alpha = 0.05, indepFilter = TRUE){
  DEresults <- results(DESeqOut, name = coef, alpha = alpha, format = "DataFrame", independentFiltering = indepFilter)
  DEresults$GeneSymbol <- geneNames$hgnc_symbol[match(rownames(DEresults), geneNames$ensembl_gene_id)]
  DEresults$EnsemblID <- rownames(DEresults)
  return(DEresults)
}

GetOneSidedPval <- function(ResultsObj, adjust = "BH"){
  DESeqResultsDF <- data.frame(ResultsObj)
  #Just for now - remove the duplicated genes (5 at this point) 
  DESeqResultsDF <- DESeqResultsDF[!duplicated(DESeqResultsDF$GeneSymbol),]
  DESeqResultsDF$DownPval <- apply(DESeqResultsDF %>% select(log2FoldChange, pvalue), 1, function(x){
    if(x[1] < 0){
      x[2]/2
    } else {
      1-x[2]
    }
  })
  DESeqResultsDF$DownPvalAdj <- p.adjust(DESeqResultsDF$DownPval, "BH")
  DESeqResultsDF$UpPval <- apply(DESeqResultsDF %>% select(log2FoldChange, pvalue), 1, function(x){
    if(x[1] > 0){
      x[2]/2
    } else {
      1-x[2]
    }
  })
  DESeqResultsDF$UpPvalAdj <- p.adjust(DESeqResultsDF$UpPval, "BH")
  rownames(DESeqResultsDF) <- as.character(DESeqResultsDF$GeneSymbol)
  return(DESeqResultsDF)
}

#### Create DESeq compatible object
CountsDESeq <- countMatrixFiltered
CountsDESeq[,-1] <- apply(CountsDESeq[-1], c(1, 2), function(x) as.integer(round(x, digits = 0)))

#Remove the samples identified as outliers
CountsDESeq %<>% select_(.dots = c("genes", as.character(studyFinal$Cortex$Metadata$RNA2)))

#Filter low expressed genes based on sex genes and excluding genes with low variance 
CountsDESeq <- CountsDESeq[CountsDESeq$genes %in% studyFinal$Cortex$aned_high$ensemblID,]

rownames(CountsDESeq) <- as.character(CountsDESeq$genes)

#Run DESeq2 on both cohorts, withou MGP adjustment
DESeqOutBoth <- DESeq2RUN(data =  CountsDESeq[-1], Meta = studyFinal$Cortex$Metadata, model = Model)
DESeqResultsBoth <- GetDESeq2Results(DESeqOutBoth, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Both <- GetOneSidedPval(ResultsObj = DESeqResultsBoth)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_Both$GeneSymbol) %>% droplevels
EnrichListPDdown_Both <- gsr(scores = DESeqResultsDF_Both, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_Both <- gsr(scores = DESeqResultsDF_Both, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")


##### Repeat for each cohort separately
## NBB
Model2 = as.formula(paste0("~Profile + ", MetaCovar))
MetaNBB <- studyFinal$Cortex$Metadata %>% filter(cohort == "NBB") %>% droplevels()
CountsNBB = CountsDESeq %>% select_(.dots = as.character(MetaNBB$RNA2))

DESeqOut_NBB <- DESeq2RUN(data =  CountsNBB, Meta = MetaNBB, model = Model2)
DESeqResults_NBB <- GetDESeq2Results(DESeqOut_NBB, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_NBB <- GetOneSidedPval(ResultsObj = DESeqResults_NBB)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_NBB$GeneSymbol) %>% droplevels
EnrichListPDdown_NBB <- gsr(scores = DESeqResultsDF_NBB, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_NBB <- gsr(scores = DESeqResultsDF_NBB, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")

##Norway
MetaNorway <- studyFinal$Cortex$Metadata %>% filter(cohort == "Norway") %>% droplevels()
CountsNorway = CountsDESeq %>% select_(.dots = as.character(MetaNorway$RNA2))

DESeqOut_Norway <- DESeq2RUN(data =  CountsNorway, Meta = MetaNorway, model = Model2)
DESeqResults_Norway <- GetDESeq2Results(DESeqOut_Norway, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Norway <- GetOneSidedPval(ResultsObj = DESeqResults_Norway)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_Norway$GeneSymbol) %>% droplevels
EnrichListPDdown_Norway <- gsr(scores = DESeqResultsDF_Norway, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_Norway <- gsr(scores = DESeqResultsDF_Norway, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")


#Repeat but correcting just for Oligo MGP
OligModel <- "Oligo_Genes"

Model2Olig <- as.formula(paste0("~Profile + ", MetaCovar, " + ", OligModel))

#Norway
MetaNorwayOlig <- MetaTemp %>% filter(cohort == "Norway") %>% droplevels()

DESeqOut_Norway_Olig <- DESeq2RUN(data =  CountsNorway, Meta = MetaNorwayOlig, model = ~Profile + sex + age_years + pm_time_min + rin + Batch + Oligo_Genes)
DESeqResults_Norway_Olig <- GetDESeq2Results(DESeqOut_Norway_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Norway_Olig <- GetOneSidedPval(ResultsObj = DESeqResults_Norway_Olig)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_Norway_Olig$GeneSymbol) %>% droplevels
EnrichListPDdown_Norway_Olig <- gsr(scores = DESeqResultsDF_Norway_Olig, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_Norway_Olig <- gsr(scores = DESeqResultsDF_Norway_Olig, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")

#NBB
MetaNBBOlig <- MetaTemp %>% filter(cohort == "NBB") %>% droplevels()

DESeqOut_NBB_Olig <- DESeq2RUN(data =  CountsNBB, Meta = MetaNBBOlig, model = ~Profile + sex + age_years + pm_time_min + rin + Batch + Oligo_Genes)
DESeqResults_NBB_Olig <- GetDESeq2Results(DESeqOut_NBB_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_NBB_Olig <- GetOneSidedPval(ResultsObj = DESeqResults_NBB_Olig)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_NBB_Olig$GeneSymbol) %>% droplevels
EnrichListPDdown_NBB_Olig <- gsr(scores = DESeqResultsDF_NBB_Olig, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_NBB_Olig <- gsr(scores = DESeqResultsDF_NBB_Olig, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")

#Compare results from both cohorts before Olig adjustment
SignifNBB <- DESeqResults_NBB %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
SignifNorway <- DESeqResults_Norway %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()

UnionSignif <- merge(DESeqResults_NBB %>% data.frame %>% filter(GeneSymbol %in% c(SignifNBB, SignifNorway)),
                     DESeqResults_Norway %>% data.frame %>% filter(GeneSymbol %in% c(SignifNBB, SignifNorway)), by = "GeneSymbol",
                     all.x=T, all.y = T, suffixes = c("_NBB", "_Norway"))

UnionSignif$Cohort <- sapply(UnionSignif$GeneSymbol, function(x){
  x = as.character(x)
  if(x %in% SignifNBB){
    "NBB"
  } else if(x %in% SignifNorway){
    "Norway"
  }
}) %>% factor

ggplot(UnionSignif, aes(log2FoldChange_Norway, log2FoldChange_NBB, color = Cohort)) +
  theme_classic() +
  labs(title = "Before Oligo correction") +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

#Compare results from both cohorts after Olig adjustment
SignifNBB_adjOlig <- DESeqResults_NBB_Olig %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()
SignifNorway_adjOlig <- DESeqResults_Norway_Olig %>% data.frame() %>% filter(padj < 0.05) %>% .$GeneSymbol %>% as.character()

UnionSignif_adjOlig <- merge(DESeqResults_NBB_Olig %>% data.frame %>% filter(GeneSymbol %in% c(SignifNBB_adj, SignifNorway_adj)),
                         DESeqResults_Norway_Olig %>% data.frame %>% filter(GeneSymbol %in% c(SignifNBB_adj, SignifNorway_adj)), by = "GeneSymbol",
                         all.x=T, all.y = T, suffixes = c("_NBB", "_Norway"))

UnionSignif_adj$Cohort <- sapply(UnionSignif_adj$GeneSymbol, function(x){
  x = as.character(x)
  if(x %in% SignifNBB_adj){
    "NBB"
  } else if(x %in% SignifNorway_adj){
    "Norway"
  }
}) %>% factor

ggplot(UnionSignif_adj, aes(log2FoldChange_Norway, log2FoldChange_NBB, color = Cohort)) +
  theme_classic() +
  labs(title = "After Oligo correction") +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

#Rerun both cohorts, adjusting form Olig
ModelOlig <- as.formula(paste0("~Profile + cohort + ", MetaCovar, " + ", OligModel))

DESeqOutBoth_Olig <- DESeq2RUN(data =  CountsDESeq[-1], Meta = MetaTemp, model = ModelOlig)
DESeqResultsBoth_Olig <- GetDESeq2Results(DESeqOutBoth_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)
DESeqResultsBoth_Olig_LFC <- lfcShrink(DESeqOutBoth_Olig, coef="Profile_PD_vs_Cont", alpha = 0.05, type = "normal")
DESeqResultsBoth_LFC <- lfcShrink(DESeqOutBoth, coef="Profile_PD_vs_Cont", alpha = 0.05, type = "normal")

DESeqResultsDF_Both_Olig <- GetOneSidedPval(ResultsObj = DESeqResultsBoth_Olig)
ExpressedGenes <- GenericHumanAnno %>% filter(ProbeName %in% DESeqResultsDF_Both_Olig$GeneSymbol) %>% droplevels
EnrichListPDdown_Both_Olig <- gsr(scores = DESeqResultsDF_Both_Olig, scoreColumn = "DownPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
EnrichListPDup_Both_Olig <- gsr(scores = DESeqResultsDF_Both_Olig, scoreColumn = "UpPvalAdj", bigIsBetter = F, logTrans = T, annotation = ExpressedGenes, aspects = "B")
