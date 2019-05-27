source("SetUp.R")
packageF("biomaRt")
packageF("parallel")
name = "Parkome"
resultsPath = "GeneralResults"

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

GeneTPM <- function(transData, mitoGenes = MitoGenes, rmGenes){
  transData %<>% filter(!gene_id %in% c(mitoGenes$ensembl_gene_id, rmGenes))
  transData %<>% mutate(RPK = 1000*expected_count/effective_length)
  transData$RPK[is.infinite(transData$RPK)|is.na(transData$RPK)] <- 0
  TotRPK = sum(transData$RPK)
  transData %<>% mutate(TPM2 = (10^6)*RPK/TotRPK)
  GeneData <- transData %>% group_by(gene_id) %>% summarise(geneTPM = sum(TPM2)) %>% data.frame()
  return(GeneData)
}

ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"), mart = ensembl)

MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]

Metadata <- read.table("mapping_ids.csv", header = T, sep = "\t")

Meta <- read.table("mapping_ids.csv", header = T, sep = "\t")
Meta$Series_sample_id <- Meta$RNA2
Meta2 <- read.table("clinical_data.csv", header = T, sep = "\t")

Metadata <- merge(Meta, Meta2 %>% select(-condition), by.x = "RNA2", by.y = "sample_id", all.x = F, all.y = F)
Metadata <- Metadata[!grepl("child", Metadata$condition, ignore.case = T),]
write.table(Metadata," Metadata.tsv", sep = "\t", row.names = F)

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
                     #Batch = lane,
                     Study = "Parkome")

#Get the count matrix and filter mitochondrial genes
countMatrix <- read.csv(paste0(name, "/data/countMatrix.genes"), header=TRUE, sep = "\t")

MaxSignal <- apply(countMatrix, 1, max)
#countMatrix <- countMatrix[MaxSignal > 5,]


CountSum <- apply(countMatrix %>% select(matches("SL")), 2, sum)

MitoCountSum <- apply(countMatrix %>% filter(genes %in% MitoGenes$ensembl_gene_id) %>% select(matches("SL")), 2, sum)

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
TopFiveProportionNoMT$ensemblID2 <- factor(TopFiveProportionNoMT$ensemblID2, levels = TopFiveGeneFreq$ensemblID2)

ggplot(TopFiveProportionNoMT, aes(Filename, Proportion, fill = ensemblID2)) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_blank()) +
  labs(y = "Proportion of reads", x = "Sample", title = "Top five genes with the highest read count") + 
  scale_fill_manual(values = c(gray.colors(9), "brown1", "goldenrod3", "aquamarine4",
                               "darkred", "cornflowerblue", "darksalmon" ),
                    name = "GeneID (n)") +
  geom_bar(stat = "identity")

#Get the top expressed genes in each sample and the sum of their counts
TopGenes <- sapply(grep("^SL", names(MitoCountFiltered), value = T), function(smpl){
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
CommonTopGenes <- TopFiveTable[TopFiveTable$n > 0.5*length(grep("SL", names(countMatrix))) | TopFiveTable$ensemblID == "ENSG00000226958",]
CommonTopGenesSum <- apply(MitoCountFiltered %>% filter(genes %in% CommonTopGenes$ensemblID) %>% select(matches("SL")), 2, sum)
TopFiveSum <- TopFiveProportionNoMT %>% group_by(Filename) %>%
  summarise(TotProp = sum(Proportion)) %>%
  data.frame %>% arrange(TotProp)

if(length(list.files(path = paste0(name, "/data"), pattern = "isoforms")) > 0){ # When the data was proccessed from BAM/FASTQ files 
  #Recreate CPM matrix after exclusion of mitochondrial genes and the top expressed genes
  AllSampleData <- sapply(list.files(path = paste0(name, "/data"), pattern = "isoforms"), function(x) strsplit(x,"\\.")[[1]][1]) 
  names(AllSampleData) <- AllSampleData
  AllSampleData %<>% as.list
  
  AllSampleData <- mclapply(AllSampleData, function(smpl){
    read.csv(paste0(name, "/data/", smpl, ".isoforms.results") , header = T, sep = "\t", quote = "")
  }, mc.cores = 30)
  
  
  #Get the effective lengths of the genes
  GeneLengthList <- lapply(AllSampleData, function(smpl){
    smpl %<>% mutate(WeighLength = effective_length*IsoPct/100)
    smpl %>% group_by(gene_id) %>% summarise(geneLength = sum(WeighLength)) %>% data.frame
  })
  
  GeneLengthDF <- sapply(names(GeneLengthList), function(x){
    GeneLengthList[[x]][2]
  }) %>% do.call(cbind, .) %>% data.frame
  
  
  colnames(GeneLengthDF) <- sapply(colnames(GeneLengthDF), function(x) strsplit(x, "\\.gene")[[1]][1]) %>%
    gsub("PEC_ASD_UCLA_.*_mRNA_HiSeq2.00_", "", .)
  rownames(GeneLengthDF) <- levels(AllSampleData[[1]]$gene_id)
  
  
  GeneSymbolAll <- data.frame(GeneSymbol = geneNames$hgnc_symbol[match(rownames(GeneLengthDF), geneNames$ensembl_gene_id)],
                              Probe = geneNames$hgnc_symbol[match(rownames(GeneLengthDF), geneNames$ensembl_gene_id)],
                              GeneType = geneNames$gene_biotype[match(rownames(GeneLengthDF), geneNames$ensembl_gene_id)],
                              ensemblID = rownames(GeneLengthDF))
  GeneLengthDF <- cbind(GeneSymbolAll, GeneLengthDF)
  GeneLengthDF$SD <- apply(GeneLengthDF %>% select(matches("GSM|_")), 1, function(gene){
    gene <- gene[gene > 0]
    round(sd(gene), digits = 2)
  })
  GeneLengthDF$Max <-  apply(GeneLengthDF %>% select(matches("GSM|_")), 1, function(x){
    max(x) 
  })
  
  GeneLengthDF %<>% filter(!ensemblID %in% MitoGenes$ensembl_gene_id) %>% droplevels()
  GeneLengthDF$GeneType <- as.character(GeneLengthDF$GeneType)
  GeneLengthDF$GeneType[is.na(GeneLengthDF$GeneType)] <- "NotAnnotated"
  GeneLengthDF$GeneType <- factor(GeneLengthDF$GeneType)
  
  #Identify short genes which are likely to correspond to non-coding genes. This is mostly for ribosomal depletion protocols
  GenesRM <- GeneLengthDF[GeneLengthDF$Max < 150 & GeneLengthDF$Max > 0,] 
  GeneRMType <- GenesRM %>% group_by(GeneType) %>% summarise(n = n()) %>% data.frame()
  
  GeneRMcount <- apply(MitoCountFiltered %>% filter(genes %in% as.character(GenesRM$ensemblID)) %>% select(matches("GSM|_")), 2, sum)
  
  GeneRMTypeCount <- sapply(GeneRMType$GeneType, function(type){
    subGenesRM <- GenesRM %>% filter(GeneType == type) %>% droplevels() 
    apply(MitoCountFiltered %>% filter(genes %in% as.character(subGenesRM$ensemblID)) %>% select(matches("GSM|_")), 2, sum)
  }, simplify = FALSE)
  names(GeneRMTypeCount) <- GeneRMType$GeneType
  
  GeneRMTypeCount %<>% do.call(rbind, .) %>% data.frame()
  
  GeneRMTypePercent <- sapply(1:length(CountSum), function(i){
    Percent = 100*GeneRMTypeCount[,i]/(CountSum[i] - MitoCountSum[i])
    sapply(Percent, function(x) round(x, digits = 2))
  }, simplify = T) %>% data.frame
  
  names(GeneRMTypePercent) <- names(GeneRMTypeCount)
  GeneRMTypePercent <- cbind(rownames(GeneRMTypeCount), GeneRMTypePercent)
  names(GeneRMTypePercent)[1] <- "GeneType"
  
  #Repeat for all genes
  GeneTypeTable <- GeneLengthDF %>% filter(Max > 0) %>% group_by(GeneType) %>% summarise(n = n()) %>% data.frame()
  
  GeneTypeCount <- sapply(GeneTypeTable$GeneType, function(type){
    subGenes <- GeneLengthDF %>% filter(Max > 0, GeneType == type) 
    apply(countMatrix %>% filter(genes %in% subGenes$ensemblID) %>% select(matches("GSM|_")), 2, sum)
  }, simplify = FALSE)
  
  names(GeneTypeCount) <- GeneTypeTable$GeneType
  
  GeneTypeCount %<>% do.call(rbind, .) %>% data.frame()
  
  GeneTypeRatio <- sapply(1:length(CountSum), function(i){
    Percent = GeneTypeCount[,i]/(CountSum[i] - MitoCountSum[i])
    sapply(Percent, function(x) round(x, digits = 2))
  }, simplify = T) %>% data.frame
  
  names(GeneTypeRatio) <- names(GeneTypeCount)
  GeneTypeRatio <- cbind(rownames(GeneTypeCount), GeneTypeRatio)
  names(GeneTypeRatio)[1] <- "GeneType"
  
  #Create filtered count matrix, after removal of the high impact genes and short reads (the effective library)
  countMatrixFiltered <- MitoCountFiltered %>% filter(!genes %in% c(as.character(GenesRM$ensemblID), as.character(CommonTopGenes$ensemblID))) %>% droplevels()
  EffectiveLibCount <- apply(countMatrixFiltered[,-1], 2, sum)
  
  #Get Proportion of differnt gene types
  SummaryTable <- data.frame(Sample = names(CountSum),
                             TotalCount = CountSum,
                             MTgeneCount = MitoCountSum,
                             TopGenesCount = CommonTopGenesSum,
                             ShortLengthCount = GeneRMcount,
                             EffectiveCount = EffectiveLibCount,
                             ProtCoding = unlist(GeneTypeCount[rownames(GeneTypeCount) == "protein_coding",]))
  
  SummaryTable %<>% mutate(MT_Ratio = round(MTgeneCount/TotalCount, digits = 2),
                           TopGenesRatioAll = round(TopGenesCount/TotalCount, digits = 2),
                           TopGenesRatio = round(TopGenesCount/(TotalCount-MTgeneCount), digits = 2),
                           ShortRatioAll = round(ShortLengthCount/TotalCount, digits = 2),
                           ShortRatio = round(ShortLengthCount/(TotalCount-MTgeneCount), digits = 2),
                           EffectiveRatioAll = round(EffectiveCount/TotalCount, digits = 2),
                           EffectiveRatio = round(EffectiveCount/(TotalCount-MTgeneCount), digits = 2),
                           ProtCodingRatioAll = round(ProtCoding/TotalCount, digits = 2),
                           ProtCodingRatio = round(ProtCoding/(TotalCount - MTgeneCount), digits = 2))
  
  
  temp <- melt(SummaryTable %>% select(matches("Sampl|MT_R|Ratio")), id.vars = "Sample", variable.name = "Type", value.name = "Proportion")
  temp$Type <- factor(temp$Type, levels =  unique(c(grep("MT_Ratio|All", unique(temp$Type), value = T),
                                                    grep("Ratio$", unique(temp$Type), value = T))))
  temp$Group <- Metadata$Profile[match(as.character(temp$Sample), Metadata$Series_sample_id)]
  temp$ReadsType <- sapply(as.character(temp$Type), function(x){
    if(grepl("MT|All", x)){
      "All Reads"
    } else {
      "MT reads excluded"
    }
  })
  
  levels(temp$Type) <- sapply(levels(temp$Type), function(x) gsub("Ratio|All", "", x))
  levels(temp$Type) <- sapply(levels(temp$Type), function(x) gsub("Effective", "EffectiveLibrary", x))
  levels(temp$Type) <- sapply(levels(temp$Type), function(x) gsub("MT_", "Mitochondrial", x))
  
  
  ggplot(temp[!grepl("Eff|Prot", temp$Type),], aes(Type, Proportion)) +
    theme_classic() +
    labs(x = "", y = "Proportion of counts") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(outlier.shape =  NA) +
    geom_jitter(width = 0.2, size = 1, aes(color = Group)) +
    facet_wrap(~ReadsType, scales = "free")
  
  ggplot(temp[grepl("Eff|Prot", temp$Type),], aes(Type, Proportion)) +
    theme_classic() +
    labs(x = "", y = "Proportion of counts") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_boxplot(outlier.shape =  NA) +
    geom_jitter(width = 0.2, size = 1, aes(color = Group)) +
    facet_wrap(~ReadsType, scales = "free")
} else {
  countMatrixFiltered <- MitoCountFiltered %>% filter(!genes %in% as.character(CommonTopGenes$ensemblID)) %>% droplevels()
}

#Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrix <- Count2CPM(MitoCountFiltered[,-1]) %>% data.frame()
cpmCountAsIs <- apply(cpmMatrix, 2, sum)
cpmMatrix <- apply(cpmMatrix, c(1,2), function(x) log2(x+1)) %>% data.frame()
cpmMatrix <- cbind(as.character(MitoCountFiltered$genes), cpmMatrix)
colnames(cpmMatrix)[1] <- "genes"

cpmCount <- apply(cpmMatrix[,-1], 2, sum)

cpmMatrixFiltered <- cpmMatrix
countMatrixFiltered <- MitoCountFiltered

ZeroCount <- apply(cpmMatrixFiltered[-1], 1, function(x){
  sum(as.numeric(x)==0)
})

cpmMatrixFiltered <- cpmMatrixFiltered[ZeroCount < 0.8*ncol(cpmMatrixFiltered),]
rownames(cpmMatrixFiltered) <- cpmMatrixFiltered$genes

countMatrixFiltered <- countMatrixFiltered[ZeroCount < 0.8*ncol(countMatrixFiltered),]
rownames(countMatrixFiltered) <- countMatrixFiltered$genes

#Add gene symbols
GeneSymbolAll <- data.frame(GeneSymbol = geneNames$hgnc_symbol[match(rownames(cpmMatrixFiltered), geneNames$ensembl_gene_id)],
                            Probe = rownames(cpmMatrixFiltered),
                            ensemblID = rownames(cpmMatrixFiltered))


ExpDataCPM <- cbind(GeneSymbolAll, cpmMatrixFiltered[-1])

ExpDataCPM <- ExpDataCPM[!is.na(ExpDataCPM$GeneSymbol),]
ExpDataCPM <- ExpDataCPM[!ExpDataCPM$GeneSymbol == "",]

RegionData <- sapply(levels(Metadata$OrgRegion), function(region){
  subMeta = Metadata %>% filter(OrgRegion == region)
  Reps <- table(as.character(subMeta$CommonName))
  subExp = ExpDataCPM %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID", as.character(subMeta$RNA2)))
  #Remove replicates with the lower correlation to other samples
  SbjCor <- cor(subExp %>% select(matches("SL")))
  diag(SbjCor) <- NA
  MedianCor <- apply(SbjCor, 2, function(x) median(x, na.rm = T)) %>% sort(decreasing = T)
  SampleNames <- data.frame(FileName = names(MedianCor),
                            CommonName = subMeta$CommonName[match(names(MedianCor), subMeta$RNA2)])
  SampleRM <- SampleNames[duplicated(SampleNames$CommonName),]
  
  subMeta <- subMeta %>%  filter(!RNA2 %in% as.character(SampleRM$FileName))
  
  subExp <- subExp[,!names(subExp) %in% as.character(SampleRM$FileName)]
  names(subExp)[-c(1:3)] <- subMeta$CommonName[match(names(subExp)[-c(1:3)], subMeta$RNA2)] %>% as.character()
  
  #Remove genes no variance
  VarGenes <- apply(subExp %>% select(matches("_")), 1, sd)
  subExp <- subExp[VarGenes > 0.01,]
  
  list(Metadata = subMeta, aned = subExp, SmplCor = SbjCor)
}, simplify = FALSE)


studyFinal <- lapply(RegionData, function(region) {
  Metadata = region$Metadata
  aned = region$aned
  source(paste0(GenScriptPath, "pre-proccessRNAseq.R"), local=T)
  output <- list(aned_high, aned_good, aned_low, MaxNoise,
                 exclude_samples_low, exclude_samples_high, exclude_samples, Metadata_org, Metadata)
  names(output) <- c("aned_high", "aned_good", "aned_low", "NoiseThreshold",
                     "exclude_samples_low", "exclude_samples_high", "exclude_samples", "Metadata_org", "Metadata")
  output
})


studyFinal <- lapply(studyFinal, function(regData) {
  regData$Metadata <- GeneSex(regData$aned_good, Metadata=regData$Metadata)
  return(regData)
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

HeatMapGen2(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath, missmatched=missmatched, save=1)


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
    temp <- cells$MeanRot
    names(temp) <- cells$CommonName
    temp
  }) %>% List2df
  studyFinal[[study]]$Metadata <- merge(studyFinal[[study]]$Metadata,
                                        estimates,
                                        by.x="CommonName",
                                        by.y="row.names",
                                        all.x=TRUE)
  studyFinal[[study]]$Metadata$Profile <- as.factor(studyFinal[[study]]$Metadata$Profile)
  studyFinal[[study]]$Metadata$Profile <- relevel(studyFinal[[study]]$Metadata$Profile, ref="Cont")
}

#Print MGP plots
sapply(names(studyFinal), function(stdName){
  if(!stdName %in% list.dirs(name, full.names = FALSE)){
    dir.create(paste0(name,"/", stdName))
  }
  meta <- studyFinal[[stdName]]$Metadata
  meta <- meta[apply(meta, 2, function(x) sum(!is.na(x))) > 0]
  names(meta) <- sapply(names(meta), function(x) gsub(" ", "_", x))
  sapply(grep("_Genes", names(meta), value = TRUE), function(mgp){
    temp <- PlotPCggplot(data=meta, CellVar = mgp,
                         name = paste(name, stdName, sep = "-"), txtSize = 16, pValSize = )
    ggsave(paste0(name,"/", stdName, "/", gsub("Genes", "MGP", mgp), ".pdf"),
           plot = temp, width = 12, height = 8, units = "in",
           dpi=300)
  })
})

rm(AllsampleData, AllNames, CountData, countMatrix, cpmMatrix, cpmMatrixFiltered, datas, ensembl,
   estimates, ExpDataCPM, geneNames, GeneSymbolAll, Meta, Meta2, missmatched,
   MitoCountFiltered, RegionData, TopFiveCount, TopFiveProportion, TopFiveProportionNoMT, TopGenes,
   CommonTopGenes, CommonTopGenesSum, cpmCount, cpmCountAsIs, MaxSignal, MitoCountSum,
   MitoCountFiltered)

save.image(paste0(GeneralResultsPath, name, ".RData"))
save(studyFinal, file = paste0(GeneralResultsPath, "studyFinal", name, ".rda"))
save(PCA_results, file = paste0(GeneralResultsPath, "PCAresults", name, ".Rda"))

#Run DE analysis
#source(paste0(ProjScriptPath, "DESeqAnalysis.R"))