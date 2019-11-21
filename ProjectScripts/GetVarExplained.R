lapply(PCA_results$GSE68719[[1]]$All, function(cell){
  temp <- cell %>% summary
  if(length(temp) > 3){
    temp$importance %>% .[2,1:3] 
  }
  }) %>% do.call(rbind,.)



lapply(studyFinal, function(Region){
  LogicData <- apply(Region$Metadata, 2, function(x){
    if(sum(is.na(x)) == length(x)){
      FALSE
      } else {
        TRUE
        }
    })
  data <- Region$Metadata[,LogicData]
  data$Profile <- relevel(data$Profile, ref = "PD")
  Cells <- grep("_Genes", names(data), value = TRUE)
  sapply(Cells, function(Cell){
    temp <- wilcox.test(as.formula(paste0(Cell, "~Profile")), data = data, conf.int = TRUE)
    data.frame(CellType = Cell,
               MinEst = signif(temp$conf.int[1], digits = 2),
               MaxEst = signif(temp$conf.int[2], digits = 2),
               Estimate = signif(temp$estimate, digits = 2),
               pVAl = signif(temp$`p.value`, digits = 2))
  }, simplify = F) %>% rbindlist()
})
