source("SetUp.R")
packageF("DESeq2")

name = "Parkome"
load(paste0(GeneralResultsPath, "Parkome.RData"))

studyFinal$Cortex$Metadata$age_years <- as.numeric(studyFinal$Cortex$Metadata$age_years)
studyFinal$Cortex$Metadata$PMI_hours <- as.numeric(studyFinal$Cortex$Metadata$pm_time_min)/60
studyFinal$Cortex$Metadata$Batch <- factor(studyFinal$Cortex$Metadata$Batch)

MetaCovar = "sex + age_years + rin + PMI_hours + Batch"
Model = as.formula(paste0("~Profile + cohort +", MetaCovar))

source("ProjectScripts/ProjectFunctions.R")

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
EnrichListPDdown_Both <- gsr(scores = DESeqResultsDF_Both, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_Both <- gsr(scores = DESeqResultsDF_Both, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))


##### Repeat for each cohort separately
## NBB
Model2 = as.formula(paste0("~Profile + ", MetaCovar))
MetaNBB <- studyFinal$Cortex$Metadata %>% filter(cohort == "NBB") %>% droplevels()
CountsNBB = CountsDESeq %>% select_(.dots = as.character(MetaNBB$RNA2))

DESeqOut_NBB <- DESeq2RUN(data =  CountsNBB, Meta = MetaNBB, model = Model2)
DESeqResults_NBB <- GetDESeq2Results(DESeqOut_NBB, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_NBB <- GetOneSidedPval(ResultsObj = DESeqResults_NBB)
EnrichListPDdown_NBB <- gsr(scores = DESeqResultsDF_NBB, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_NBB <- gsr(scores = DESeqResultsDF_NBB, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))


##Norway
MetaNorway <- studyFinal$Cortex$Metadata %>% filter(cohort == "Norway") %>% droplevels()
CountsNorway = CountsDESeq %>% select_(.dots = as.character(MetaNorway$RNA2))

DESeqOut_Norway <- DESeq2RUN(data =  CountsNorway, Meta = MetaNorway, model = Model2)
DESeqResults_Norway <- GetDESeq2Results(DESeqOut_Norway, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Norway <- GetOneSidedPval(ResultsObj = DESeqResults_Norway)
EnrichListPDdown_Norway <- gsr(scores = DESeqResultsDF_Norway, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_Norway <- gsr(scores = DESeqResultsDF_Norway, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))


MetaTemp = studyFinal$Cortex$Metadata
#Repeat but correcting just for Oligo MGP
OligModel <- "Oligo_Genes"

Model2Olig <- as.formula(paste0("~Profile + ", MetaCovar, " + ", OligModel))

#Norway
MetaNorwayOlig <- MetaTemp %>% filter(cohort == "Norway") %>% droplevels()

DESeqOut_Norway_Olig <- DESeq2RUN(data =  CountsNorway, Meta = MetaNorwayOlig, model = Model2Olig)
DESeqResults_Norway_Olig <- GetDESeq2Results(DESeqOut_Norway_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Norway_Olig <- GetOneSidedPval(ResultsObj = DESeqResults_Norway_Olig)
EnrichListPDdown_Norway_Olig <- gsr(scores = DESeqResultsDF_Norway_Olig, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_Norway_Olig <- gsr(scores = DESeqResultsDF_Norway_Olig, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))

#NBB
MetaNBBOlig <- MetaTemp %>% filter(cohort == "NBB") %>% droplevels()

DESeqOut_NBB_Olig <- DESeq2RUN(data =  CountsNBB, Meta = MetaNBBOlig, model = Model2Olig)
DESeqResults_NBB_Olig <- GetDESeq2Results(DESeqOut_NBB_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_NBB_Olig <- GetOneSidedPval(ResultsObj = DESeqResults_NBB_Olig)
EnrichListPDdown_NBB_Olig <- gsr(scores = DESeqResultsDF_NBB_Olig, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_NBB_Olig <- gsr(scores = DESeqResultsDF_NBB_Olig, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))


#Repeat but correcting just for General neuron MGP
NeuronModel <- "NeuronAll_Genes"

Model2Neuron <- as.formula(paste0("~Profile + ", MetaCovar, " + ", NeuronModel))

#Norway
MetaNorwayNeuron <- MetaTemp %>% filter(cohort == "Norway") %>% droplevels()

DESeqOut_Norway_Neuron <- DESeq2RUN(data =  CountsNorway, Meta = MetaNorwayNeuron, model = Model2Neuron)
DESeqResults_Norway_Neuron <- GetDESeq2Results(DESeqOut_Norway_Neuron, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_Norway_Neuron <- GetOneSidedPval(ResultsObj = DESeqResults_Norway_Neuron)
EnrichListPDdown_Norway_Neuron <- gsr(scores = DESeqResultsDF_Norway_Neuron, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_Norway_Neuron <- gsr(scores = DESeqResultsDF_Norway_Neuron, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))

#NBB
MetaNBBNeuron <- MetaTemp %>% filter(cohort == "NBB") %>% droplevels()

DESeqOut_NBB_Neuron <- DESeq2RUN(data =  CountsNBB, Meta = MetaNBBNeuron, model = Model2Neuron)
DESeqResults_NBB_Neuron <- GetDESeq2Results(DESeqOut_NBB_Neuron, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)

DESeqResultsDF_NBB_Neuron <- GetOneSidedPval(ResultsObj = DESeqResults_NBB_Neuron)
EnrichListPDdown_NBB_Neuron <- gsr(scores = DESeqResultsDF_NBB_Neuron, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_NBB_Neuron <- gsr(scores = DESeqResultsDF_NBB_Neuron, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))

#Compare results from both cohorts before Olig adjustment
CompareResultsAll(data1 = DESeqResults_NBB, data2 = DESeqResults_Norway, name1 = "NBB",
               name2 = "Norway", colorCol = "Cohort",
               Title = "Before Oligo correction, no pH correction")


#Compare results from both cohorts after Olig adjustment
CompareResultsAll(data1 = DESeqResults_NBB_Olig, data2 = DESeqResults_Norway_Olig, name1 = "NBB",
               name2 = "Norway", colorCol = "Cohort",
               Title = "After Oligo correction")



#Rerun both cohorts, adjusting for Olig
ModelOlig <- as.formula(paste0("~Profile + cohort + ", MetaCovar, " + ", OligModel))

DESeqOutBoth_Olig <- DESeq2RUN(data =  CountsDESeq[-1], Meta = MetaTemp, model = ModelOlig)
DESeqResultsBoth_Olig <- GetDESeq2Results(DESeqOutBoth_Olig, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)


DESeqResultsDF_Both_Olig <- GetOneSidedPval(ResultsObj = DESeqResultsBoth_Olig)
EnrichListPDdown_Both_Olig <- gsr(scores = DESeqResultsDF_Both_Olig, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))
EnrichListPDup_Both_Olig <- gsr(scores = DESeqResultsDF_Both_Olig, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = c("B", "M", "C"))

rm(datas, cpmCountFiltered, cpmCountFiltered, ensembl, estimates, ExpDataAll, ExpDataCPM)
save.image(file = paste0(GeneralResultsPath, "DESeqAnalysisParkome.Rdata"))


AdjCovar <- data.frame(Cov = c("sex_M_vs_F", 
                               "Batch_2_vs_1", "Batch_3_vs_1", "Batch_4_vs_1",
                               "age_years","PMI_hours", "rin",
                               "Oligo_Genes"),
                       adjType = c(rep("base", 4), rep("mean", 4)))
PTPRH_PV <- plotCounts(DESeqOut_Norway_Olig, gene = "ENSG00000080031", intgroup = "condition", returnData = T) %>% data.frame
PTPRH_PV %<>% mutate(AdjPTPRH_Olig = GetAdjCountDESeq(dds = DESeqOut_Norway_Olig, Gene = "ENSG00000080031", adjCov = AdjCovar),
                     AdjPTPRH = GetAdjCountDESeq(dds = DESeqOut_Norway, Gene = "ENSG00000080031", adjCov = AdjCovar %>% filter(Cov != "Oligo_Genes")),
                     count = log2(count),
                     RNA2 = rownames(attr(DESeqOut_Norway_Olig, "modelMatrix")),
                     Cohort = "PV")

AdjCovar <- data.frame(Cov = c("sex_M_vs_F", "Batch_1_vs_0",
                               "Batch_2_vs_0", "Batch_3_vs_0", "Batch_4_vs_0",
                               "age_years","PMI_hours", "rin",
                               "Oligo_Genes"),
                       adjType = c(rep("base", 5), rep("mean", 4)))
PTPRH_NBB <- plotCounts(DESeqOut_NBB_Olig, gene = "ENSG00000080031", intgroup = "condition", returnData = T) %>% data.frame
PTPRH_NBB %<>% mutate(AdjPTPRH_Olig = GetAdjCountDESeq(dds = DESeqOut_NBB_Olig, Gene = "ENSG00000080031", adjCov = AdjCovar),
                     AdjPTPRH = GetAdjCountDESeq(dds = DESeqOut_NBB, Gene = "ENSG00000080031", adjCov = AdjCovar %>% filter(Cov != "Oligo_Genes")),
                     count = log2(count),
                     RNA2 = rownames(attr(DESeqOut_NBB_Olig, "modelMatrix")),
                     Cohort = "NBB")
PTPRH_both <- rbind(PTPRH_PV, PTPRH_NBB) %>% data.frame() %>% gather(key = "CountType", value = "PTPRHcount", AdjPTPRH_Olig, AdjPTPRH)
PTPRH_both$condition <- relevel(PTPRH_both$condition,ref = "Control")
levels(PTPRH_both$condition) <- c("Cont", "PD")
PTPRH_both$CountType <- factor(PTPRH_both$CountType, levels = c("AdjPTPRH", "AdjPTPRH_Olig"))
levels(PTPRH_both$CountType) <- c("NotOligoAdjusted", "OligoAdjusted")

Plot <- ggplot(PTPRH_both, aes(condition, PTPRHcount, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "log2(PTPRH)") +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.5) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_grid(CountType~Cohort, scales = "free_y")
ggsave("PTPRHgeneLevels.pdf", plot = Plot, device = "pdf", width =  7, height =5, dpi = 300, useDingbats=F)
write.table(PTPRH_both, file = "PTPRH_bothCohorts.tsv", row.names = F, col.names = T, sep = "\t")
