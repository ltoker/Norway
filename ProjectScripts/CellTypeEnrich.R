TasicDataOligiDowb <- read.table("~/CellularProportionsPsychiatry/GeneralResults/TasicExpCellsScaled2.txt",
                                 sep = "\t", header = T, quote = '"', comment = "!")
rownames(TasicDataOligiDowb) <- TasicDataOligiDowb$Probe

TasicDataOligiUp <- read.table("~/CellularProportionsPsychiatry/GeneralResults/TasicExpCellsScaled.txt",
                                 sep = "\t", header = T, quote = '"')
rownames(TasicDataOligiUp) <- TasicDataOligiUp$Probe


NeuroExpressoDataOligiDowb <- read.table("~/CellularProportionsPsychiatry/GeneralResults/NeuroExpCellsScaled2.txt",
                                 sep = "\t", header = T, quote = '"', comment = "!")
rownames(NeuroExpressoDataOligiDowb) <- NeuroExpressoDataOligiDowb$Probe

NeuroExpressoDataOligiUp <- read.table("~/CellularProportionsPsychiatry/GeneralResults/NeuroExpCellsScaled.txt",
                               sep = "\t", header = T, quote = '"')
rownames(NeuroExpressoDataOligiUp) <- NeuroExpressoDataOligiUp$Probe

#Get the annotation platforms from Gemma
GPL339 = GetAnnoFiles("GPL339")
MouseGeneralGenes <- GetAnnoFiles("Generic_mouse")
HumanGeneralGenes <- GetAnnoFiles("Generic_human")

aspects = c("C", "B", "M")
TasicOligoEnrichDown <- gsr(scores = TasicDataOligiDowb, scoreColumn = "Oligodendrocyte", bigIsBetter = T, logTrans = F, annotation = MouseGeneralGenes, aspects = aspects)
TasicOligoEnrichUp <- gsr(scores = TasicDataOligiUp, scoreColumn = "Oligodendrocyte", bigIsBetter = T, logTrans = F, annotation = MouseGeneralGenes, aspects = aspects)

NeuExprOligoEnrichDown <- gsr(scores = NeuroExpressoDataOligiDowb, scoreColumn = "Oligodendrocyte", bigIsBetter = T, logTrans = F, annotation = GPL339, aspects = aspects)
NeuExprOligoEnrichUp <- gsr(scores = NeuroExpressoDataOligiUp, scoreColumn = "Oligodendrocyte", bigIsBetter = T, logTrans = F, annotation = GPL339, aspects = aspects)

source("ProjectScripts/HumanSingleCell.R")
#Repeat for Darmais data
DarmanisMetaFiltered <- DarmanisMeta %>% filter(!cellType %in% c("hybrid", "fetal_quiescent", "fetal_replicating")) %>% droplevels()
DarmanisExpFiltered <- DarmanisExp[colnames(DarmanisExp) %in% DarmanisMetaFiltered$GSM]

DarmanisScaled <- scale(DarmanisExpFiltered %>% t) %>% t %>% data.frame()
rownames(DarmanisScaled) <- rownames(DarmanisExpFiltered)

DarmanisExpCellsScaled <- sapply(unique(DarmanisMetaFiltered$cellType), function(celltype){
  samples <- DarmanisMetaFiltered %>% filter(cellType == celltype) %>% .$GSM %>% make.names()
  expr <- DarmanisScaled %>%  select_(.dots = samples)
  rowMeans(expr)
}) %>% data.frame

rownames(DarmanisExpCellsScaled) <-  rownames(DarmanisScaled)
DarmanisExpCellsScaled2 <- data.frame(apply(DarmanisExpCellsScaled, 2, function(x) {-1*x}))

DarmanisScaledDF <- DarmanisExpCellsScaled %>% mutate(GeneSymbol = rownames(.))

DarmanisOligoEnrichDown <- gsr(scores = DarmanisExpCellsScaled2, scoreColumn = "oligodendrocytes", bigIsBetter = T, logTrans = F, annotation = HumanGeneralGenes, aspects = aspects)
DarmanisOligoEnrichUp <- gsr(scores = DarmanisExpCellsScaled, scoreColumn = "oligodendrocytes", bigIsBetter = T, logTrans = F, annotation = HumanGeneralGenes, aspects = aspects)

DarmanisNeuronEnrichDown <- gsr(scores = DarmanisExpCellsScaled2, scoreColumn = "neurons", bigIsBetter = T, logTrans = F, annotation = HumanGeneralGenes, aspects = aspects)
DarmanisNeuronEnrichUp <- gsr(scores = DarmanisExpCellsScaled, scoreColumn = "neurons", bigIsBetter = T, logTrans = F, annotation = HumanGeneralGenes, aspects = aspects)


NEneuronVSglia <- read.table("NEnuronsVSglia.tsv", header = T, sep = "\t")
NEneuronVSglia %<>% filter(logFC > 3.2)
NEneuronVSglia$HugoName <- sapply(as.character(NEneuronVSglia$Symbol), function(x){
  temp <- mouse2human(x)$humanGene
  if(length(temp) == 0){
    NA
  } else if(length(temp) > 1){
    paste0(temp, collapse = "|")
  } else {
    temp
  }
}) %>% unlist

NEneuronVSglia$HumanNeurEnich <- sapply(NEneuronVSglia$HugoName, function(x){
  if(is.na(x)){
    temp = NA
  } else {
    temp <- DarmanisScaledDF %>% filter(GeneSymbol == x) %>% .$neurons
    if(length(temp) == 0 ){
      temp = NA
    }
  }
  temp
}) %>% unlist

NEneuronVSglia %<>% filter(HumanNeurEnich > 0.2)
write.table(NEneuronVSglia, "HumanNeuronVSglia", sep = "\t", row.names = FALSE, col.names = TRUE)



