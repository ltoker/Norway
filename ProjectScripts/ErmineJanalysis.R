#Run ErmineJ

#Getting the NeuroExpresso cortical cell expression profiles
if(length(list.files(pattern = "n_expressoMeta.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/n_expressoSamples.rda?raw=true", destfile = "n_expressoMeta.rda")
} else {
  warning("using existing file")
}

load("n_expressoMeta.rda")

#Removes samples not included in NeuroExpresso
n_expressoSamples %<>% filter(!is.na(.$PyramidalDeep))

if(length(list.files(pattern = "n_expressoExpr.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/n_expressoExpr.rda?raw=true", destfile = "n_expressoExpr.rda")
} else {
  warning("using existing file")
}

load("n_expressoExpr.rda")

if(length(list.files(pattern = "regionHierarchy.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/regionHierarchy.rda?raw=true", destfile = "regionHierarchy.rda")
} else {
  warning("using existing file")
}
load("regionHierarchy.rda")

#Scale the expression for ErmineJ analysis
n_expressoExpr<- n_expressoExpr[!grepl("\\|", n_expressoExpr$Gene.Symbol),]
NeuroExpScaled <- scale(n_expressoExpr[names(n_expressoExpr) %in% n_expressoSamples$sampleName] %>% as.matrix %>% t) %>% t %>% data.frame()
rownames(NeuroExpScaled) <- as.character(n_expressoExpr$Probe)

CortexCells <- n_expressoSamples[grepl("Cortex", n_expressoSamples$Region),] %>% .$ShinyNames %>% unique

NeuroExpCellsScaled <- sapply(CortexCells, function(celltype){
  samples <- n_expressoSamples %>% filter(ShinyNames == celltype) %>% .$sampleName
  expr <- NeuroExpScaled %>%  select_(.dots = samples)
  rowMeans(expr)
}) %>% data.frame() %>% mutate(Probe = rownames(.))

names(NeuroExpCellsScaled) <- sapply(names(NeuroExpCellsScaled), function(x) gsub("\\.", "", x))
NeuroExpCellsScaled %<>% select(Probe, Astrocyte,  Microglia, Oligodendrocyte,
                                FSBasketG42, MartinottiGIN,  VIPRelnG30, PyramidalGlt25d2, PyramidalS100a10, PyramidalCrtThalamic )

rownames(NeuroExpCellsScaled) <- as.character(NeuroExpCellsScaled$Probe)

# Download annotation file
if(length(list.files(pattern = "GPL339.gz")) == 0){
  download.file(paste0("http://chibi.ubc.ca/microannots/", "GPL339",  "_noParents.an.txt.gz"), destfile="GPL339.gz")
} else {
  warning("using existing GPL file")
}

GPL339 <- read.table("GPL339.gz", header = TRUE, sep = "\t", quote = "")

write.table(GPL339 %>% filter(ProbeName %in% n_expressoExpr$Probe), "NeuExprProbes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

if(!"ermineR" %in% rownames(installed.packages())){
  install_github("oganm/ermineR", force = T)
}
library(ermineR)

CommonGeneAnnoFile <- GPL339 %>% filter(ProbeName %in% n_expressoExpr$Probe) %>% droplevels
EnrichList <- sapply(names(NeuroExpCellsScaled)[-1], function(celltype){
  gsr(scores = NeuroExpCellsScaled, scoreColumn = celltype, bigIsBetter = T, logTrans = F, annotation = CommonGeneAnnoFile, aspects = "B")
}, simplify = FALSE)

#Repeat for Tasic data
if(length(list.files(pattern = "TasicMouseExp.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/TasicMouseExp.rda?raw=true", destfile = "TasicMouseExp.rda")
}  else {
  warning("using existing file")
}
load("TasicMouseExp.rda")

if(length(list.files(pattern = "TasicMouseMeta.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/TasicMouseMeta.rda?raw=true", destfile = "TasicMouseMeta.rda")
}  else {
  warning("using existing file")
}
load("TasicMouseMeta.rda")

if(length(list.files(pattern = "TasicNEmatch.rda")) == 0){
  download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/meltedSingleCells.rda?raw=true", destfile = "TasicNEmatch.rda")
}  else {
  warning("using existing file")
}
load("TasicNEmatch.rda")

TasicMouseMeta$ShinyNames <- meltedSingleCells$ShinyNames[match(TasicMouseMeta$primary_type, meltedSingleCells$sampleName)]

TasicExpScaled <- scale(TasicMouseExp %>% as.matrix %>% t) %>% t %>% data.frame()

TasicExpCellsScaled <- sapply(unique(TasicMouseMeta$ShinyNames), function(celltype){
  samples <- TasicMouseMeta %>% filter(ShinyNames == celltype) %>% .$sample_title %>% make.names()
  expr <- TasicExpScaled %>%  select_(.dots = samples)
  rowMeans(expr)
}) %>% data.frame() %>% mutate(Probe = rownames(.))

names(TasicExpCellsScaled) <- sapply(names(TasicExpCellsScaled), function(x) gsub("\\.| ", "", x))
TasicExpCellsScaled <- TasicExpCellsScaled[-2]
TasicExpCellsScaled %<>% select(ncol(TasicExpCellsScaled), 1:c(ncol(TasicExpCellsScaled)-1))
rownames(TasicExpCellsScaled) <- as.character(TasicExpCellsScaled$Probe)

write.table(TasicExpCellsScaled, "GeneralResults/TasicExpCellsScaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

if(length(list.files(pattern = "Generic_mouse_noParents.an.txt")) == 0){
  download.file("http://chibi.ubc.ca/microannots/Generic_mouse_noParents.an.txt.gz", destfile="Generic_mouse_noParents.an.txt")
} else {
  warning("using existing GPL file")
}
GenericMouseAnno <-  read.table("Generic_mouse_noParents.an.txt", header = TRUE, sep = "\t", quote = "")

CommonGeneAnnoFileTasic <- GenericMouseAnno %>% filter(ProbeName %in% TasicExpCellsScaled$Probe) %>% droplevels

EnrichListTasic <- sapply(names(TasicExpCellsScaled)[-1], function(celltype){
  gsr(scores = TasicExpCellsScaled, scoreColumn = celltype, bigIsBetter = T, logTrans = F, annotation = CommonGeneAnnoFileTasic, aspects = "B")
}, simplify = FALSE)