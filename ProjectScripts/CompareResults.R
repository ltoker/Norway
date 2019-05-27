packageF("pROC")

CompareResults <- function(data1, data2, name1 = "NBB", name2 = "Norway", colorCol = "Cohort", Title, FcCol1 = "log2FoldChange", FcCol2 = "log2FoldChange"){
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

CompareResultsAll <- function(data1, data2, name1 = "NBB", name2 = "Norway", colorCol = "Cohort",
                              pThresh = 0.05,Title, FcCol1 = "log2FoldChange", FcCol2 = "log2FoldChange"){
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
  }) %>% factor(levels = c("NS", name1, name2, "Both"))
  CheckSig <- function(x, pThresh){
    temp <- if(is.na(x)){
      NA
    } else {
      if(x < pThresh){
        "Yes"
      } else {
        "No"
      }
    }
    return(temp)
  }

  UnionSignif[[paste0(name1, "Up")]] <- sapply(UnionSignif[[paste0("UpPvalAdj_", name1)]], function(x) {CheckSig(x, pThresh)})
  UnionSignif[[paste0(name1, "Down")]] <- sapply(UnionSignif[[paste0("DownPvalAdj_", name1)]], function(x) {CheckSig(x, pThresh)})
  UnionSignif[[paste0(name2, "Up")]] <- sapply(UnionSignif[[paste0("UpPvalAdj_", name2)]], function(x) {CheckSig(x, pThresh)})
  UnionSignif[[paste0(name2, "Down")]] <- sapply(UnionSignif[[paste0("DownPvalAdj_", name2)]], function(x) {CheckSig(x, pThresh)})
  
  CompareCol <- grep(paste0(FcCol1, "|", FcCol2), names(UnionSignif), value = T)
  signifData <- UnionSignif %>% filter(Cohort != "NS")
  signifData$Cohort <- factor(signifData$Cohort, levels = c(name1, name2, "Both"))
  
  p <- ggplot(UnionSignif, aes_string(CompareCol[1], CompareCol[2])) +
    theme_classic() +
    labs(title = Title) +
    ylim(-2,2) +
    xlim(-2,2) +
    geom_point(color = "grey", alpha = 0.5) +
    geom_point(data = signifData,
               aes_string(CompareCol[1], CompareCol[2],
                          color = "Cohort"),inherit.aes = FALSE) +
    scale_color_discrete(drop=FALSE) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
  return(list(data = UnionSignif,
              GenePlot = p))
}

PlotROC <- function(data, response, predictor, title = NULL){
  par(mfrow = c(1,2))
  data %<>% select_(.dots = c(response, predictor))
  data$Predictor <- data[[predictor]]
  if(is.null(title)){
    title = paste0(predictor, " vs ", response)
  }
  formString <- paste0(response, "~ Predictor")
  ROC <- roc(formula(formString), data, levels = c("No", "Yes"), na.rm = T, 
             direction = "<", plot = TRUE, print.auc = T,
             ylim = c(0,1),
             print.auc.x = 0.4, main = title)
  plot(precision ~ recall, t(coords(ROC, "all", ret = c("recall", "precision"))), type="l")
  par(mfrow = c(1,1))
}

GetROC <- function(data, response, predictor, title = NULL){
  formString <- paste0(response, "~" , predictor)
  ROC <- roc(formula(formString), data, levels = c("Yes", "No"), na.rm = T)
             
  return(ROC)
}

PlotCompar <- function(data, title){
  mfrow(c(1,2))
  plot(data$GenePlot)
  
}
Parkome <-new.env()
GSE20295 <- new.env()

load("GeneralResults/DESeqAnalysisParkome.Rdata", envir = Parkome)
load("GeneralResults/GSE20295.RData", envir = GSE20295)



DESeqResultsDF_Parkome_Both_Olig <- get("DESeqResultsDF_Both_Olig", envir = Parkome)
DESeqResultsDF_Parkome_Both <- get("DESeqResultsDF_Both", envir = Parkome)
DESeqResultsDF_Parkome_Norway <- get("DESeqResultsDF_Norway", envir = Parkome)
DESeqResultsDF_Parkome_Norway_Olig <- get("DESeqResultsDF_Norway_Olig", envir = Parkome)
DESeqResultsDF_Parkome_NBB <- get("DESeqResultsDF_NBB", envir = Parkome)
DESeqResultsDF_Parkome_NBB_Olig <- get("DESeqResultsDF_NBB_Olig", envir = Parkome)



GSE68719_NBB <- CompareResultsAll(data1 = DESeqResultsDF_Both,
                                  data2 = DESeqResultsDF_Parkome_NBB,
                                  name1 = "GSE68719",
                                  name2 = "NBB", colorCol = "Cohort",
                                  Title = "GSE68719 vs NBB, before MGP correction")


GSE68719_ParkVest <- CompareResultsAll(data1 = DESeqResultsDF_Both,
                                       data2 = DESeqResultsDF_Parkome_Norway,
                                       name1 = "GSE68719",
                                       name2 = "ParkVest", colorCol = "Cohort",
                                       Title = "GSE68719 vs ParkVest, before MGP correction")
                                                  

GSE68719_Parkome <- CompareResultsAll(data1 = DESeqResultsDF_Both,
                                      data2 = DESeqResultsDF_Parkome_Both,
                                      name1 = "GSE68719",
                                      name2 = "Parkome", colorCol = "Cohort",
                                      Title = "GSE68719 vs Parkome, before oligo MGP correction")

#After correction
GSE68719_NBB_after <- CompareResultsAll(data1 = DESeqResultsDF_Both_Oligo,
                                        data2 = DESeqResultsDF_Parkome_NBB_Olig,
                                        name1 = "GSE68719",
                                        name2 = "NBB", colorCol = "Cohort",
                                        Title = "GSE68719 vs NBB, after oligo MGP correction")

GSE68719_NBB_afterMicroglia <- CompareResultsAll(data1 = DESeqResultsDF_Both_Microglia,
                                        data2 = DESeqResultsDF_Parkome_NBB_Olig,
                                        name1 = "GSE68719",
                                        name2 = "NBB", colorCol = "Cohort",
                                        Title = "GSE68719 vs NBB, after microglia MGP correction")
                        

GSE68719_ParkVest_afterMicroglia <- CompareResultsAll(data1 = DESeqResultsDF_Both_Microglia,
                                             data2 = DESeqResultsDF_Parkome_Norway_Olig,
                                             name1 = "GSE68719",
                                             name2 = "ParkVest", colorCol = "Cohort",
                                             Title = "GSE68719 vs ParkVest, after microglia MGP correction")

GSE68719_ParkVest_after <- CompareResultsAll(data1 = DESeqResultsDF_Both_Oligo,
                                             data2 = DESeqResultsDF_Parkome_Norway_Olig,
                                             name1 = "GSE68719",
                                             name2 = "ParkVest", colorCol = "Cohort",
                                             Title = "GSE68719 vs ParkVest, after oligo MGP correction")

GSE68719_Parkome_after <- CompareResultsAll(data1 = DESeqResultsDF_Both_Oligo,
                                            data2 = DESeqResultsDF_Parkome_Both_Olig,
                                            name1 = "GSE68719",
                                            name2 = "Parkome", colorCol = "Cohort",
                                            Title = "GSE68719 vs Parkome, after oligo MGP correction")

GSE68719_NBB$ROCdown <- GetROC(GSE68719_NBB$data,
                               response = "GSE68719Down",
                               predictor = "DownPvalAdj_NBB",
                               title = "GSE68719Down, NBB")
GSE68719_NBB$ROCup <- GetROC(GSE68719_NBB$data,
                             response = "GSE68719Up",
                             predictor = "UpPvalAdj_NBB",
                             title = "GSE68719Up, NBB")
GSE68719_ParkVest$ROCdown <- GetROC(GSE68719_ParkVest$data,
                                    response = "GSE68719Down",
                                    predictor = "DownPvalAdj_ParkVest",
                                    title ="GSE68719Down,ParkVest")
GSE68719_ParkVest$ROCup <- GetROC(GSE68719_ParkVest$data,
                                  response = "GSE68719Up",
                                  predictor = "UpPvalAdj_ParkVest", title = "GSE68719Up,ParkVest")
GSE68719_Parkome$ROCdown <- GetROC(GSE68719_Parkome$data,
                                   response = "GSE68719Down",
                                   predictor = "DownPvalAdj_Parkome", title = "GSE68719Down, Parkome")
GSE68719_Parkome$ROCup <- GetROC(GSE68719_Parkome$data,
                                 response = "GSE68719Up",
                                 predictor = "UpPvalAdj_Parkome",
                                 title = "GSE68719Up, Parkome")

GSE68719_NBB_after$ROCdown <- GetROC(GSE68719_NBB_after$data,
                                     response = "GSE68719Down",
                                     predictor = "DownPvalAdj_NBB",
                                     title = "State: GSE68719Down, NBB (After)")
GSE68719_NBB_after$ROCup <- GetROC(GSE68719_NBB_after$data,
                                   response = "GSE68719Up",
                                   predictor = "UpPvalAdj_NBB",
                                   title = "State: GSE68719Up, NBB (After)")
GSE68719_ParkVest_after$ROCdown <- GetROC(GSE68719_ParkVest_after$data,
                                          response = "GSE68719Down",
                                          predictor = "DownPvalAdj_ParkVest",
                                          title = "GSE68719Down, ParkVest (After)")
GSE68719_ParkVest_after$ROCup <- GetROC(GSE68719_ParkVest_after$data,
                                        response = "GSE68719Up",
                                        predictor = "UpPvalAdj_ParkVest",
                                        title = "GSE68719Up, ParkVest (After)")
GSE68719_Parkome_after$ROCdown <- GetROC(GSE68719_Parkome_after$data,
                                         response = "GSE68719Down",
                                         predictor = "DownPvalAdj_Parkome",
                                         title = "GSE68719Down, Parkome (After)")
GSE68719_Parkome_after$ROCup <- GetROC(GSE68719_Parkome_after$data,
                                       response = "GSE68719Up",
                                       predictor = "UpPvalAdj_Parkome",
                                       title = "GSE68719Up, Parkome (After)")


rm(list = ls(pat = "^Top|Matrix|^Mito|Data|ensembl|^Count|countM|Common|Expressed|^PCA|Model|Zero|Covar|cpmC|Parkome|Generic|missmatched|GeneSymbolAll|estimates|geneNames|mouseMark|study|subCol"))
#save.image(file = paste0(GeneralResultsPath, "CompareResultsParkome_GSE68719.Rdata"))
