source("SetUp.R")
packageF("biomaRt")
packageF("tidyr")
packageF("parallel")
name = "Parkome"
resultsPath = "GeneralResults"

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}


ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "uniprot_gn_symbol", "ensembl_gene_id", "gene_biotype"), mart = ensembl)
geneNames$hgnc_symbol <- apply(geneNames, 1, function(gene){
  if(gene[1] == ""){
    gene[2]
  } else {
    gene[1]
  }
})
MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]

Meta <- read.table(paste0(name, "/meta/mapping_ids.csv"), header = T, sep = "\t")
Meta$Series_sample_id <- Meta$RNA2
Meta2 <- read.table(paste0(name, "/meta/clinical_data.csv"), header = T, sep = "\t")

Metadata <- merge(Meta, Meta2 %>% select(-condition), by.x = "RNA2", by.y = "sample_id", all.x = F, all.y = F)
Metadata <- Metadata[!grepl("child", Metadata$condition, ignore.case = T),]
write.table(Metadata, paste0(name, "/meta/Metadata.tsv"), sep = "\t", row.names = F)

Metadata$Profile <- sapply(Metadata$condition, function(x){
  x <- as.character(x)
  if(x == "Control"){
    "Cont"
  } else {
    "PD"
  }
})

Metadata %<>% arrange(Profile)
Metadata$Profile <- factor(Metadata$Profile, levels = c("Cont", "PD"))

Metadata$CommonName <- sapply(levels(Metadata$Profile), function(x){
  paste0(x, "_", 1:nrow(Metadata %>% filter(Profile == x)))
}, simplify = FALSE) %>% unlist

Metadata$OrgRegion = factor("Cortex")
Metadata %<>% mutate(NeuExpRegion = OrgRegion,
                     Filename = RNA2,
                     Series_sample_id = RNA2,
                     Study = "Parkome")


#Get the count matrix and filter mitochondrial genes
countMatrix <- read.csv(paste0(name, "/data/countMatrix.genes"), header=TRUE, sep = "\t")

CohortData <- sapply(levels(Metadata$cohort), function(Cohort){
  subMeta = Metadata %>% filter(cohort == Cohort) %>% droplevels()
  subExp = countMatrix %>% select(c("genes", as.character(subMeta$RNA2)))
  names(subExp)[-1] <- subMeta$CommonName[match(names(subExp)[-1], subMeta$RNA2)] %>% as.character()
  SbjCor <- cor(subExp %>% select(matches("_")))
 
  #Remove genes no variance
  VarGenes <- apply(subExp %>% select(matches("_")), 1, sd)
  subExp <- subExp[VarGenes > 0.01,]
  list(Metadata = subMeta, countMatrix = subExp)
}, simplify = FALSE)


CohortData <- lapply(CohortData, function(Cohort){
  countMatrix = Cohort$countMatrix
  Metadata = Cohort$Metadata
  CountSum <- apply(countMatrix %>% select(matches("_")), 2, sum)
  
  MitoCountSum <- apply(countMatrix %>% filter(genes %in% MitoGenes$ensembl_gene_id) %>% select(matches("_")), 2, sum)
  
  MitoCountFiltered <- countMatrix %>% filter(!genes %in% MitoGenes$ensembl_gene_id)
  
  MitoFiltCountSum = apply(MitoCountFiltered[-1], 2, sum)
  
  TopFiveProportion <- sapply(names(CountSum), function(sbj){
    SubMatrix = countMatrix %>% select_(.dots = c("genes", sbj))
    names(SubMatrix)[2] <- "Counts"
    TopFive = SubMatrix %>% arrange(desc(Counts)) %>% head(5)
    TopFive %<>%  mutate(Proportion = Counts/CountSum[sbj])
    Genes <- geneNames[match(TopFive$genes, geneNames$ensembl_gene_id),]  %>% select(-ensembl_gene_id)
    Genes$Filename = sbj
    temp <- cbind(Genes, TopFive)
    names(temp)[names(temp) == "genes"] <- "ensemblID"
    temp
  }, simplify = FALSE) %>% rbindlist()
  
  TopFiveProportionNoMT <- sapply(names(CountSum), function(sbj){
    SubMatrix = MitoCountFiltered %>% select_(.dots = c("genes", sbj))
    names(SubMatrix)[2] <- "Counts"
    TopFive = SubMatrix %>% arrange(desc(Counts)) %>% head(5)
    TopFive %<>%  mutate(Proportion = Counts/MitoFiltCountSum[sbj])
    Genes <- geneNames[match(TopFive$genes, geneNames$ensembl_gene_id),] %>% select(-ensembl_gene_id)
    Genes$Filename = sbj
    temp <- cbind(Genes, TopFive)
    names(temp)[names(temp) == "genes"] <- "ensemblID"
    temp
  }, simplify = FALSE) %>% rbindlist() %>% data.frame()
  
  TopFiveProportionNoMT <- merge(TopFiveProportionNoMT, geneNames, by.x = "ensemblID", by.y = "ensembl_gene_id", all.x = T, all.y = F)
  
  TopFiveSum <- TopFiveProportionNoMT %>% group_by(Filename) %>%
    summarise(TotProp = sum(Proportion)) %>%
    data.frame %>% arrange(TotProp)
  
  
  TopFiveGeneFreq <- TopFiveProportionNoMT %>% group_by(ensemblID) %>%
    summarise(n = n()) %>%
    data.frame
  
  TopFiveGeneFreq %<>% mutate(ensemblID2 = paste0(ensemblID, " (", n, ")"))
  TopFiveGeneFreq <- merge(TopFiveGeneFreq, geneNames[!duplicated(geneNames$ensembl_gene_id),], by.x = "ensemblID", by.y = "ensembl_gene_id", all.x = T, all.y = F)
  TopFiveGeneFreq %<>% arrange(n)
  
  TopFiveProportionNoMT$ensemblID2 <- TopFiveGeneFreq$ensemblID2[match(TopFiveProportionNoMT$ensemblID, TopFiveGeneFreq$ensemblID)]
  
  TopFiveProportionNoMT$Filename <- factor(TopFiveProportionNoMT$Filename, levels = TopFiveSum$Filename)
  TopFiveProportionNoMT$ensemblID2 <- factor(TopFiveProportionNoMT$ensemblID2, levels = rev(as.character(TopFiveGeneFreq$ensemblID2)))
  
  Plot <- ggplot(TopFiveProportionNoMT, aes(Filename, Proportion, fill = ensemblID2)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_blank()) +
    labs(y = "Proportion of reads", x = "Sample", title = "Top five genes with the highest read count") + 
    scale_fill_manual(values = c("brown1", "goldenrod3", "aquamarine4",
                                 "darkred", "cornflowerblue", "darksalmon", gray.colors(6)),
                      name = "GeneID (n)") +
    geom_bar(stat = "identity")
  print(Plot)  
  
  #Get the top expressed genes in each sample and the sum of their counts
  TopGenes <- sapply(grep("_", names(MitoCountFiltered), value = T), function(smpl){
    TopFive <- MitoCountFiltered %>% arrange_(.dots = smpl) %>% tail(5)
    TopFive %<>% select_(.dots = c("genes", smpl))
    names(TopFive)[2] <- "Counts"
    TopFive %>% arrange(desc(Counts))
  }, simplify = FALSE)
  
  TopFiveCount <- lapply(TopGenes, function(smpl){
    sum(smpl$Counts)
  }) %>% unlist
  
  TopFiveNames <- lapply(TopGenes, function(smpl){
    smpl$genes %>% as.character
  }) %>% unlist %>% data.frame()
  
  names(TopFiveNames) <- "ensemblID"
  TopFiveNames %<>% mutate(GeneSymbol = geneNames$hgnc_symbol[match(TopFiveNames$ensemblID, geneNames$ensembl_gene_id)],
                           GeneType = geneNames$gene_biotype[match(TopFiveNames$ensemblID, geneNames$ensembl_gene_id)])
  
  TopFiveTable <- TopFiveNames %>% group_by(ensemblID) %>% summarise(n = n()) %>% data.frame() %>% arrange(desc(n))
  TopFiveTable %<>% mutate(GeneSymbol = geneNames$hgnc_symbol[match(TopFiveTable$ensemblID, geneNames$ensembl_gene_id)],
                           GeneType = geneNames$gene_biotype[match(TopFiveTable$ensemblID, geneNames$ensembl_gene_id)])
  
  #Get the common genes with the highest count in all the samples (and the ribosomal gene..)
  CommonTopGenes <- TopFiveTable[TopFiveTable$n > 0.5*length(grep("_", names(countMatrix))) | TopFiveTable$ensemblID == "ENSG00000226958",]
  CommonTopGenesSum <- apply(MitoCountFiltered %>% filter(genes %in% CommonTopGenes$ensemblID) %>% select(matches("_")), 2, sum)
  TopFiveSum <- TopFiveProportionNoMT %>% group_by(Filename) %>%
    summarise(TotProp = sum(Proportion)) %>%
    data.frame %>% arrange(TotProp)
  
  countMatrixFiltered <- MitoCountFiltered %>% filter(!genes %in% as.character(CommonTopGenes$ensemblID)) %>% droplevels()
  ZeroCount <- apply(countMatrixFiltered[-1], 1, function(x){
    sum(as.numeric(x)==0)
  })
  
  countMatrixFiltered <- countMatrixFiltered[ZeroCount < 0.8*ncol(countMatrixFiltered),]
  rownames(countMatrixFiltered) <- countMatrixFiltered$genes
  #Create log2 CPM matrix after removal of mitochondria-encoded genes
  cpmMatrixFiltered <- Count2CPM(countMatrixFiltered[,-1]) %>% data.frame()
  cpmMatrixFiltered <- apply(cpmMatrixFiltered, c(1,2), function(x) log2(x+1)) %>% data.frame()
  cpmMatrixFiltered <- cbind(as.character(countMatrixFiltered$genes), cpmMatrixFiltered)
  colnames(cpmMatrixFiltered)[1] <- "genes"
  
  cpmCountFiltered <- apply(cpmMatrixFiltered[,-1], 2, sum)
  
  #Add gene symbols
  GeneSymbolAll <- data.frame(GeneSymbol = geneNames$hgnc_symbol[match(rownames(cpmMatrixFiltered), geneNames$ensembl_gene_id)],
                              Probe = rownames(cpmMatrixFiltered),
                              ensemblID = rownames(cpmMatrixFiltered))
  
  
  ExpDataCPM <- cbind(GeneSymbolAll, cpmMatrixFiltered[-1])
  
  ExpDataCPM <- ExpDataCPM[!is.na(ExpDataCPM$GeneSymbol),]
  ExpDataCPM <- ExpDataCPM[!ExpDataCPM$GeneSymbol == "",]
  
  #Remove duplicated gene symbols
  ExpDataCPM %<>% arrange(desc(.[,4]))
  
  ExpDataCPM %<>% filter(!duplicated(GeneSymbol))
  
  list(Metadata = Metadata, aned = ExpDataCPM, SmplCor = Cohort$SbjCor, countMatrix = countMatrixFiltered)
})


studyFinal <- lapply(CohortData, function(Cohort) {
  Metadata = Cohort$Metadata
  Metadata$age_years <- as.numeric(Metadata$age_years)
  Metadata$PMI_hours <- as.numeric(Metadata$pm_time_min)/60
  Metadata$Batch <- factor(Metadata$Batch)
  aned = Cohort$aned
  RunComBat = FALSE
  source(paste0(GenScriptPath, "pre-proccessRNAseq.R"), local=T)
  output <- list(aned_high, aned_good, aned_low, MaxNoise,
                 exclude_samples_low, exclude_samples_high, exclude_samples, Metadata_org, Metadata)
  names(output) <- c("aned_high", "aned_good", "aned_low", "NoiseThreshold",
                     "exclude_samples_low", "exclude_samples_high", "exclude_samples", "Metadata_org", "Metadata")
  output
})


studyFinal <- lapply(studyFinal, function(cohData) {
  cohData$Metadata <- GeneSex(cohData$aned_good, Metadata=cohData$Metadata)
  return(cohData)
})

missmatched <- lapply(studyFinal, function(x){
  meta <- x$Metadata
  names(meta) <- tolower(names(meta))
  meta$sex <- sapply(meta$sex, function(sex){
    if(grepl("^male|^man|^m$", tolower(sex))){
      "M"
    } else if(grepl("female|^wom|w|^f$", tolower(sex))){
      "F"
    }
  })
  meta$commonname[meta$sex != meta$biogender]
})

#Create gender HeatMaps
datas <- datasGenerate(c("XIST", "KDM5D", "RPS4Y1"))

#Estimating cell type proportions
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

#PCA analysis without the missmatched samples
PCA_results <- as.list(names(studyFinal))
names(PCA_results) <- names(studyFinal)

PCA_results <- mclapply(PCA_results, function(x){
  region = studyFinal[[x]]$Metadata$NeuExpRegion %>% unique
  CellType_genes <- GetMarkers(region)
  
  #Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
  if(region == "Cortex"){
    CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]
    temp <- read.table("HumanNeuronVSglia.tsv", sep = "\t", header = T)
    CellType_genes$NeuronAll_Genes <- temp$HugoName %>% as.character()
  }
  
  aned_high <- studyFinal[[x]]$aned_high
  aned_high <- aned_high[,!names(aned_high) %in% missmatched[[x]]]
  
  #bootstrap with replacement the samples in each group to ensure equal number of samples/group (90% of the samples in the smaller group)
  groups <- sapply(grep("_", names(aned_high), value=T),
                   function(x) gsub("_.*", "", x)) %>% table
  MinGrp <- round(0.9*min(groups))
  AllSamples <- sapply(names(groups), function(grp){
    grep(grp, names(aned_high), value = T)
  }, simplify = FALSE)
  results <- list()
  for(i in c(1:100)){
    BootSamples <- lapply(AllSamples, function(grp){
      grp[sample(1:length(grp), MinGrp,replace = FALSE)]
    }) %>% unlist %>% as.character
    
    aned_highSub <- aned_high %>% select_(.dots = c("GeneSymbol", BootSamples))
    results[[i]] <- PCA_genes_All_based(dataset_id=x,
                                        dataset=aned_highSub,
                                        CellType_genes=CellType_genes,
                                        NoiseThershold = studyFinal[[x]]$NoiseThreshold)
  }
  
  return(results)
}, mc.cores = length(PCA_results))

PCA_resultsMean <- lapply(PCA_results, function(region){
  MeanPCA <- sapply(names(region[[1]]$modified), function(celltype){
    temp <- data.frame(CommonName = names(region[[1]]$modified[[celltype]]$x[,1]),
                       Rot = region[[1]]$modified[[celltype]]$x[,1])
    for(i in 2:length(region)){
      temp <- merge(temp, data.frame(CommonName = names(region[[i]]$modified[[celltype]]$x[,1]),
                                     Rot = region[[i]]$modified[[celltype]]$x[,1]), by = "CommonName", all = TRUE)
    }
    names(temp)[2:ncol(temp)] <- paste0("Rot", c(1:c(ncol(temp)-1)))
    temp$MeanRot <- rowMeans(temp[-1], na.rm = T)
    temp
  }, simplify=FALSE)
})

#Add estimation to Metadata 
for(study in names(studyFinal)){
  studyFinal[[study]]$Metadata %<>% select(-matches("_Genes"))
  estimates <- lapply(PCA_resultsMean[[study]], function(cells){
    temp <- cells %>% select(MeanRot)
    rownames(temp) <- cells$CommonName
    temp
  }) %>% do.call(cbind, .)
  names(estimates) <- names(PCA_resultsMean[[study]])
  studyFinal[[study]]$Metadata <- merge(studyFinal[[study]]$Metadata,
                                        estimates,
                                        by.x="CommonName",
                                        by.y="row.names",
                                        all.x=TRUE,
                                        sort = F)
  studyFinal[[study]]$Metadata$Profile <- as.factor(studyFinal[[study]]$Metadata$Profile)
  studyFinal[[study]]$Metadata$Profile <- relevel(studyFinal[[study]]$Metadata$Profile, ref="Cont")
}

#Add the count matrix, remove outlier samples and low expressed genes based on previous steps
for(Cohort in names(studyFinal)){
  temp <-CohortData[[Cohort]]$countMatrix
  temp <- temp[temp$genes %in% studyFinal[[Cohort]]$aned_high$ensemblID,]
  temp[,-1] <- apply(temp[-1], c(1, 2), function(x) as.integer(round(x, digits = 0)))
  temp %<>% select("genes", grep("_", names(studyFinal[[Cohort]]$aned_high), value = T))
  studyFinal[[Cohort]]$countMatrix <- temp 
}

#Print MGP plots
# sapply(names(studyFinal), function(stdName){
#   if(!stdName %in% list.dirs(name, full.names = FALSE)){
#     dir.create(paste0(name,"/", stdName))
#   }
#   meta <- studyFinal[[stdName]]$Metadata
#   meta <- meta[apply(meta, 2, function(x) sum(!is.na(x))) > 0]
#   names(meta) <- sapply(names(meta), function(x) gsub(" ", "_", x))
#   sapply(grep("_Genes", names(meta), value = TRUE), function(mgp){
#     temp <- PlotPCggplot(data=meta, CellVar = mgp,
#                          name = paste(name, stdName, sep = "-"), txtSize = 16, pValSize = )
#     ggsave(paste0(name,"/", stdName, "/", gsub("Genes", "MGP", mgp), ".pdf"),
#            plot = temp, width = 12, height = 8, units = "in",
#            dpi=300)
#   })
# })

rm(AllsampleData, AllNames, CountData, countMatrix, cpmMatrix, cpmMatrixFiltered, datas, ensembl,
   estimates, ExpDataCPM, geneNames, GeneSymbolAll, Meta, Meta2, missmatched,
   MitoCountFiltered, RegionData, TopFiveCount, TopFiveProportion, TopFiveProportionNoMT, TopGenes,
   CommonTopGenes, CommonTopGenesSum, cpmCount, cpmCountAsIs, MaxSignal, MitoCountSum,
   MitoCountFiltered)

save.image(paste0(GeneralResultsPath, name, ".RData"))
save(studyFinal, file = paste0(GeneralResultsPath, "studyFinal", name, ".rda"))
save(PCA_results, file = paste0(GeneralResultsPath, "PCAresults", name, ".Rda"))

ggplot(studyFinal$Norway$Metadata, aes(Profile,Oligo_Genes)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(width = 0.2)

#Run DE analysis
source(paste0(ProjScriptPath, "DESeqAnalysis.R"))
