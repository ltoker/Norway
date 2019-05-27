source(paste0(GenScriptPath,"graphic_functions.R")) #networks and heatmaps
source(paste0(GenScriptPath,"correlation_functions.R")) #corrlations and distributions functions

OutSamples <- function(data){
  Q1 = data %>% summary %>% .["1st Qu."]
  Q3 = data %>% summary %>% .["3rd Qu."]
  IQR = Q3-Q1
  Min = Q1 - IQR
  Max = Q3 + IQR
  return(names(data)[data < Min])
}

Metadata %<>% droplevels()
if(length(ls(pat="resultsPath", envir = .GlobalEnv)) == 0){
  print("No results path specified, using working directory")
  resultsPath = paste0(getwd(), "/")
}
if(length(ls(pat="study")) == 0){
  study = ""
}
if(length(ls(pat="region")) == 0){
  region = ""
}

aned_good <- aned[complete.cases(aned),]

#Remove batch effect
source(paste0(GenScriptPath,"ComBat_script.R"), local=T)

#Get max signal for noise data based on expression of known non-expressed genes
Fgene <- grep("XIST", aned_good$GeneSymbol, value=TRUE, ignore.case = T) %>% unique
Mgene <- grep("KDM5D|RPS4Y1", aned_good$GeneSymbol, value=TRUE, ignore.case = T)

if(length(c(Fgene, Mgene)) == 0){
  if(nrow(aned_good) > 20000){
    print("no gender genes detected, setting noise threshold to median")
    MaxNoise = aned_good[-c(1:3)] %>% as.matrix() %>% median(na.rm =T)
    } else {
      print("no gender genes detected, assuming filtered data,setting noise threshold to Zero")
      MaxNoise = 0
      }
  } else {
    if(length(unique(Metadata$Sex)) == 2 | length(unique(Metadata$Sex)) == 0){
      GeneGender <- GeneSex(aned=aned_good,
                            Metadata=data.frame(CommonName = names(aned_good[sapply(aned_good,
                                                                                    function(x) is.numeric(x))])))
      if("Gender genes disagree, cannot decide about biological gender" %in% names(warnings())){
        print("Can't use Sex genes to determine noise threshold, setting noise threshold to median")
        MaxNoise = aned_good[-c(1:3)] %>% median(is.na =T)
      } else {
        Noise <- sapply(c(Fgene, Mgene), function(gene){
          if(gene %in% Fgene){
            gender = "M"
          } else if (gene %in% Mgene) {
            gender = "F"
          }
          #Detect samples with potential quality issues and remove them from noise calculation
          if(length(Fgene) > 0 & length(Mgene) > 0){
            FemaleExp = aned_good %>% filter(GeneSymbol %in% Fgene) %>% select(-c(1:3)) %>% apply(2, mean)
            MaleExp = aned_good %>% filter(GeneSymbol %in% Mgene) %>% select(-c(1:3)) %>% apply(2, mean)
            SexDiff  = MaleExp - FemaleExp
            RmSample = names(SexDiff)[SexDiff > -3 & SexDiff < 2]
            if(length(RmSample) > 0){
              print(paste("Suspected data quality in:", paste(RmSample, collapse = ", ")))
            }
          } else{
            RmSample = "none"
          }
          temp <- aned_good %>% filter(GeneSymbol == gene) %>%
            select_(.dots=GeneGender %>% filter(BioGender == gender, !CommonName %in% RmSample) %>% .$CommonName %>% as.character)  %>% unlist
          #quantile(temp, 0.9)
        })  %>% unlist
        
        MaxNoise <- quantile(Noise, 0.95)
        if(is.na(MaxNoise)){
          print("Can't use sex genes, setting noise threshold to median")
          MaxNoise = aned_good[-c(1:3)] %>% as.matrix() %>% median(is.na =T)
        }
      }
      
    } else if(length(unique(Metadata$Sex)) == 1){
      print("Only one Meta-sex, setting noise threshold to median")
      MaxNoise = aned_good[-c(1:3)] %>% as.matrix() %>% median(is.na =T)
    }
    
  }
    


print(paste("Noise threshold:", MaxNoise))

#Define genes above and below noise threshold
ProbeSum <- apply(aned_good %>% select(matches("GSM|_")), 1, function(x) quantile(x, 0.95) > MaxNoise)
aned_high <- aned_good[ProbeSum,]
aned_low <- aned_good[!ProbeSum,]

#print("Getting the heatmaps")
#Create Z-score heatmaps for low and high signals
samples <- 4:ncol(aned_good)
CorSamplesGood <- cor(aned_good %>% select(matches("GSM|_")), method = "spearman")
diag(CorSamplesGood) <- NA
MedianCorGood <- apply(CorSamplesGood, 2, function(x) median(x, na.rm = T)) %>% sort()

png(paste0(name, "/", name, "SampleCorAllExp.png"))
boxplot(data.frame(CorSamplesGood) %>% select_(.dots = names(MedianCorGood)),outline = F, las = 2)
dev.off()

CorSamplesHigh <- cor(aned_good %>% select(matches("GSM|_")), method = "spearman")
diag(CorSamplesHigh) <- NA
MedianCorHigh <- apply(CorSamplesHigh, 2, function(x) median(x, na.rm = T)) %>% sort()

png(paste0(name, "/", name, "SampleCorHighExp.png"))
boxplot(data.frame(CorSamplesHigh) %>% select_(.dots = names(MedianCorHigh)), outline = F, las = 2)
dev.off()
exclude_samples_high <- OutSamples(MedianCorHigh)

if(nrow(aned_low) > 3) {
  CorSamplesLow <- cor(aned_good %>% select(matches("GSM|_")), method = "spearman")
  diag(CorSamplesLow) <- NA
  MedianCorLow <- apply(CorSamplesLow, 2, function(x) median(x, na.rm = T)) %>% sort()
  png(paste0(name, "/", name, "SampleCorLowExp.png"))
  boxplot(data.frame(CorSamplesLow) %>% select_(.dots = names(MedianCorLow)), outline = F, las = 2)
  dev.off()
  exclude_samples_low <- OutSamples(MedianCorLow)
} else {
  exclude_samples_low <- NULL
}

exclude_samples <- union(exclude_samples_high, exclude_samples_low)
Metadata_org <- Metadata

if(ncol(aned_good) > 20){
  aned_good_org <- aned_good
  aned_good <- aned_good[,!names(aned_good) %in% exclude_samples]
  
  Metadata <- Metadata[!Metadata$CommonName %in% exclude_samples,]
} else {
  if(length(exclude_samples) > 0){
    warning(paste0(paste0(exclude_samples, collapse = ", "), "are suspected outliers. Less than 20 samples, no outliers removed")
)  }
}

#Define genes above and below noise threshold
ProbeSum <- apply(aned_good[,-c(1:3)], 1, function(x) quantile(x, 0.95) > MaxNoise)
aned_high <- aned_good[ProbeSum,]
aned_low <- aned_good[!ProbeSum,]

#Multiple-probset treatment
aned_high <- aned_high[names(rev(sort(apply(aned_high[,-c(1:3)], 1, sd)))),]
aned_high <- aned_high[!duplicated(aned_high[,2]),]

#print("Working on the last histograms")
#screen(3) ; hist(cor(t(aned_low[,-c(1:3)])), main="Low signals end", freq=F)
#screen(4) ; hist(cor(t(aned_high[,-c(1:3)])), main="High signals end", freq=F)

#close.screen(all.screens=T)
closeDev()
