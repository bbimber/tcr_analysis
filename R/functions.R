library(lattice)
library(gplots)
library(ggplot2)
library(reshape)
library(dplyr)
library(Rlabkey)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(data.table)
#library(plyr)
library(dtplyr)

pullTCRMetaFromLabKey <- function(requireGeneFile = FALSE){
  df <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Internal/Bimber/145", 
    schemaName="lists", 
    queryName="TCR_Datasets", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c('ReadsetId','ReadsetId/application','ReadsetId/status','ReadsetId/workbook','ReadsetId/totalForwardReads','StimId','StimId/AnimalId','StimId/Date','StimId/Peptide','StimId/Treatment','Population','Replicate','Cells','StimId/ActivatedFreq','CellClass','StimId/Background','ReadsetId/numCDR3s','ReadsetId/distinctLoci','ReadsetId/numTcrRuns','SequenceComments','GroupId','GeneTable','Activated','Metadata/geneCountFiles','Metadata/estimatedLibrarySize', 'Key'), 
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=makeFilter(c('ReadsetId', 'NON_BLANK', ''))
  )
  
  names(df)[names(df)=="stimid_animalid"] <- "AnimalId"
  df$AnimalId <- as.factor(df$AnimalId)
  
  names(df)[names(df)=="readsetid"] <- "ReadsetId"
  df$ReadsetId <- as.factor(df$ReadsetId)
  
  names(df)[names(df)=="metadata_estimatedlibrarysize"] <- "EstimatedLibrarySize"
  df$EstimatedLibrarySize <- as.numeric(df$EstimatedLibrarySize)
  
  names(df)[names(df)=="readsetid_totalforwardreads"] <- "TotalReads"
  df$TotalReads <- as.integer(df$TotalReads)
  
  names(df)[names(df)=="stimid_peptide"] <- "Peptide"
  df$Peptide <- as.factor(df$Peptide)
  
  names(df)[names(df)=="cellclass"] <- "Cellclass"
  df$Cellclass <- as.factor(df$Cellclass)
  
  names(df)[names(df)=="readsetid_numcdr3s"] <- "NumCDR3s"
  df$NumCDR3s <- as.integer(df$NumCDR3s)
  
  names(df)[names(df)=="population"] <- "Population"
  df$Population <- as.factor(df$Population)
  
  names(df)[names(df)=="stimid"] <- "StimId"
  df$StimId <- as.factor(df$StimId)
  
  names(df)[names(df)=="stimid_treatment"] <- "Treatment"
  df$Treatment <- as.factor(df$Treatment)
  
  names(df)[names(df)=="readsetid_distinctloci"] <- "DistinctLoci"
  df$DistinctLoci <- as.factor(df$DistinctLoci)
  
  names(df)[names(df)=="readsetid_numtcrruns"] <- "RunTCRRuns"
  df$RunTCRRuns <- as.integer(df$RunTCRRuns)

    names(df)[names(df)=="metadata_genecountfiles"] <- "OutputFileId"
  df$OutputFileId <- as.factor(df$OutputFileId)
  
  names(df)[names(df)=="readsetid_application"] <- "Application"
  df$Application <- as.factor(df$Application)
  
  names(df)[names(df)=="stimid_date"] <- "SampleDate"
  names(df)[names(df)=="key"] <- "DatasetId"
  
  names(df)[names(df)=="readsetid_workbook"] <- "WorkbookId"
  df$WorkbookId <- as.factor(df$WorkbookId)
  
  if (requireGeneFile){
    df <- df[!is.na(df$OutputFileId),]  
  }
  
  df$Status <- df$readsetid_status
  df <- df[is.na(df$Status),]
  
  # df$ResponseType <- c('Other')
  # df$ResponseType[df$Peptide %in% c('Gag120', 'Gag69')] <- c('Supertope')
  # df$ResponseType[df$Peptide %in% c('Gag16', 'Gag18', 'Gag23', 'Gag30', 'Gag33', 'Gag50', 'Gag65', 'Gag97', 'Gag109', 'Gag119')] <- c('Other E-Restricted')
  # df$ResponseType[df$Peptide %in% c('Gag27', 'Gag37', 'QI9', 'LF8')] <- c('Other Class II-Restricted')
  # df$ResponseType[df$Peptide %in% c('RL8', 'RL9', 'RL10')] <- c('Conventional')
  # df$ResponseType <- as.factor(df$ResponseType)
  
  df <- df %>% 
    group_by(StimId, Population, replicate, Application) %>% 
    mutate(totalCellsForGroup = sum(cells))
  
  return(df)
}

pullTCRResultsFromLabKey <- function(colFilter = NULL){
  df <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Internal/Bimber/", 
    schemaName="assay.TCRdb.TCRdb", 
    queryName="Data", 
    viewName="", 
    colSelect=c('sampleName', 'subjectId','date','libraryId','libraryId/label','analysisId', 'analysisId/readset', 'locus','vHit','dHit','jHit','cHit','vGene', 'jGene','CDR3','count','fraction', 'disabled'), 
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=colFilter
  )
  
  df <- df[!is.na(df$cdr3),]
  df <- df[is.na(df$disabled),]
  df$LocusCDR3 <- paste0(df$locus, ":", df$cdr3)
  df$LocusCDR3 <- as.factor(df$LocusCDR3)
  
  df$LocusCDR3WithUsage <- paste0(df$locus, ":", df$cdr3,":",df$vgene,"-",df$jgene)
  df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  
  dateFormat <- date_format(format ="%Y-%m-%d")
  df$dateFormatted <- dateFormat(df$date)
  
  df$locus <- as.factor(df$locus)
  
  return(df)
}

pullResultsFromLabKey <- function(resultFilter = NULL){
  results <- pullTCRResultsFromLabKey(colFilter = resultFilter)  
  meta <- pullTCRMetaFromLabKey(requireGeneFile = FALSE)
  
  results <- merge(results, meta, by.x=c('analysisid_readset'), by.y=c('ReadsetId'))
  results <- qualityFilterResults(results)
  results <- groupSingleCellData(results)
  results <- filterBasedOnCellCount(results)
  
  return(results)
}

qualityFilterResults <- function(df = NULL, minLibrarySize=2000000, minReadCount=500000){
  print(paste0('performing quality filter'))
  
  r <- nrow(df)
  df <- df[df$EstimatedLibrarySize > minLibrarySize,]
  print(paste0('rows dropped for library size: ', (r - nrow(df))))
  
  r <- nrow(df)
  df <- df[df$TotalReads > minReadCount,]
  print(paste0('rows dropped for read count: ', (r - nrow(df))))
  
  return(df)
}

groupSingleCellData <- function(df = NULL, lowFreqThreshold = 0.05, minTcrReads = 50){
  #print(nrow(df))
  df <- df %>% 
    group_by(StimId, Population, replicate, Application) %>% 
    mutate(totalTcrReadsForGroup = sum(count))
  
  df <- df %>% 
    group_by(StimId, Population, replicate, Application, locus) %>% 
    mutate(totalTcrReadsForGroupAndLocus = sum(count))

  #bulk
  dfb <- df[df$Application != 'RNA-seq, Single Cell',]
  dfb$percentage = dfb$count / dfb$totalTcrReadsForGroup
  dfb$percentageForLocus = dfb$count / dfb$totalTcrReadsForGroupAndLocus
  dfb <- dfb[dfb$totalTcrReadsForGroup >= minTcrReads,]
  dfb$totalCellsForCDR3 <- c(0)

  columns <- c('StimId', 'Population', 'Application', 'AnimalId', 'SampleDate', 'Peptide', 'Treatment', 'locus', 'vhit', 'jhit', 'cdr3', 'count', 'date', 'dateFormatted', 'LocusCDR3', 'LocusCDR3WithUsage', 'Cellclass', 'replicate', 'totalCellsForGroup', 'cells', 'percentage', 'percentageForLocus', 'totalCellsForCDR3')
  dfb <- dfb[,columns]
  print(paste0('total bulk rows: ', nrow(dfb)))
  
  #single cell
  dfs <- df[df$Application == 'RNA-seq, Single Cell',]
  print(paste0('starting single cell rows: ', nrow(dfs)))
  dfs <- dfs %>%
    group_by(StimId, Population, Application, AnimalId, SampleDate, Peptide, Treatment, locus, vhit, jhit, cdr3, date, dateFormatted, LocusCDR3, LocusCDR3WithUsage, Cellclass, replicate, totalCellsForGroup) %>%
    summarize(totalCellsForCDR3 = sum(cells))
  dfs$cells <- c(1)
  dfs$count <- dfs$totalCellsForCDR3
  
  dfs$percentage <- dfs$totalCellsForCDR3 / dfs$totalCellsForGroup
  dfs$percentageForLocus <- dfs$totalCellsForCDR3 / dfs$totalCellsForGroup
  
  print(paste0('total single cell rows: ', nrow(dfs)))
  
  df <- rbind(dfb[columns], dfs[columns])
  print(paste0('combined: ', nrow(df)))
  
  df <- df %>% group_by(LocusCDR3, AnimalId) %>% mutate(maxPercentageInAnimal = max(percentage)) 
  df <- df %>% group_by(LocusCDR3, Peptide) %>% mutate(maxPercentageInPeptide = max(percentage)) 
  
  if (lowFreqThreshold){
    df$LocusCDR3 <- as.character(df$LocusCDR3)
    df$LocusCDR3[df$maxPercentageInAnimal < lowFreqThreshold] <- c('Low Freq.')
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    
    df$LocusCDR3WithUsage <- as.character(df$LocusCDR3WithUsage)
    df$LocusCDR3WithUsage[df$maxPercentageInAnimal < lowFreqThreshold] <- c('Low Freq.')
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  }
  
  df <- df %>% group_by(StimId, Population, replicate, Application) %>% mutate(totalTcrReadsForGroupAfterLowFreq = sum(count)) 
  
  #TODO
  df$percentageAfterLowFreq <- df$count / df$totalTcrReadsForGroupAfterLowFreq
  
  df <- arrange(df, LocusCDR3)
  
  return(df)
}

filterBasedOnCellCount <- function(df = NULL){
  startRows <- nrow(df)
  df$minPctToInclude <- c(0)
  
  if (length(df$Application == 'RNA-seq, Single Cell') > 0){
    df$minPctToInclude[df$Application == 'RNA-seq, Single Cell'] <- c(0.1)
    df$minPctToInclude[df$Application == 'RNA-seq, Single Cell'] <- (1 / df$cells[df$Application == 'RNA-seq, Single Cell']) 
  }

    #df <- df[df$percentageForLocus >= df$minPctToInclude | df$percentage >= (df$minPctToInclude / 4),]
  print(paste0('rows filtered: ', (startRows - nrow(df))))
  
  return(df)
}

makePlot <- function(locus, df, pal = 'Reds', bulkOnly = FALSE, tnfPosOnly = FALSE, showLegend = FALSE, yMax=1.05, yField='percentage', xField='Label', fillField='LocusCDR3'){
  locusData <- df[df$locus == locus,]
  
  if (bulkOnly){
    locusData <- locusData[locusData$Application == 'RNA-seq',]
  }
  
  if (tnfPosOnly){
    locusData <- locusData[locusData$Population == 'TNF-Pos',]
  }
  
  locusData$LocusCDR3 <- as.factor(as.character(locusData$LocusCDR3))
  totalCDR3 = length(unique(locusData$LocusCDR3))
  print(paste(locus, ', total CDR3', totalCDR3))
  
  if (nrow(locusData) == 0){
    return(ggplot(data.frame()) +geom_blank())
  }
  
  getPalette = colorRampPalette(brewer.pal(9, pal))
  
  colorsDT <-  data.table(LocusCDR3=levels(locusData$LocusCDR3), Color=getPalette(totalCDR3))
  locusData <- data.table(locusData)
  locusData <- merge(locusData, colorsDT, by = c("LocusCDR3"), all = TRUE)
  locusData$Color[locusData$LocusCDR3 == 'Low Freq.'] <- c('#D3D3D3')
  locusData$Color <- as.factor(locusData$Color)
  
  colorValues <- getPalette(totalCDR3)
  #colorValues <- rev(colorValues)
  colorValues[[1]] <- '#D3D3D3'  #grey
  
  #w <- as.numeric(difftime(max(df$SampleDate), min(df$SampleDate), units = c('days')))
  #w <- w / length(unique(df$SampleDate))
  #w <- w / 2
  
  if (showLegend){
    lp = 'right'
  }
  else {
    lp = 'none'
  }
  
  P1<-ggplot(locusData) +
    aes_string(x = xField, y = yField, order = 'percentage', width = 0.5) +
    geom_bar(aes_string(fill = fillField), stat = 'identity', colour="black") +
    scale_fill_manual(values = colorValues) +
    theme_bw() +
    theme(
      legend.position=lp,
      legend.text=element_text(size=12,color="black"),
      legend.title = element_blank(),
      plot.margin = unit(c(0.3,0.3,0.3,0.3),"in"),
      plot.title = element_text(size=24,vjust=2, hjust = 0.5),
      axis.title.x = element_text(size=20,vjust=-1),
      axis.title.y = element_text(size=20,vjust=0),
      plot.background = element_rect(color = 'white'),
      axis.text = element_text(size=16,color="black", vjust = 0),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.text.x = element_text(size = 14,face="bold"),
      strip.text.y = element_text(size = 14,face="bold")
    ) +
    ylim(c(0,yMax)) +
    guides(fill=guide_legend(ncol=1)) +
    xlab('') +
    labs(title = paste0(locus)) +
    #scales = "free_y"
    facet_grid(Population ~ dateFormatted, scales = 'fixed') +
    ylab('')
  
  return(P1)
}

makePlotGroup <- function(df, subDir, basename, tnfPosOnly=FALSE, includeTRB=FALSE, includeTRG=FALSE, pal='Reds', showLegend=FALSE, bulkOnly=FALSE, yMax=NA, yField='percentage', xField='Label', fillField='LocusCDR3'){
  plots = list()
  i = 1
  
  for (locus in sort(unique(df$locus))){
    P1 <- makePlot(locus, df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly, yField=yField, xField=xField, fillField=fillField)
    plots[[i]] <- P1
    i=i+1
  }
  
  suffix = ''
  if (bulkOnly){
    suffix = paste0(suffix, '.tnf')
  }
  
  if (tnfPosOnly){
    suffix = paste0(suffix, '.bulk')
  }
  
  w <- 1400
  #if (showLegend){
  #  w <- 1800
  #}
  
  if (!file.exists(subDir)){
    dir.create(file.path('./', subDir))
  }
  
  fn <- paste0('./',subDir,'/', basename, suffix)
  
  write.table(df, file=paste0(fn, '.txt'), sep='\t', row.names = FALSE)
  
  png(paste0(fn,'.png'), height = 1000, width = w)
  do.call(grid.arrange, plots, c(ncol=2))
  dev.off()
  
  if (includeTRB){
    png(paste0(fn,'.TRB.png'), height = 500, width = 700)
    P_B <- makePlot('TRB', df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly,yField='percentageForLocus')
    P_B = P_B + labs(title = paste0(basename, " TRB"))
    print(P_B)
    dev.off()
  }
  
  if (includeTRG){
    png(paste0(fn,'.TRG.png'), height = 500, width = 700)
    P_G <- makePlot('TRG', df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly,yField='percentageForLocus')
    P_G = P_G + labs(title = paste0(basename, " TRG"))
    print(P_G)
    dev.off() 
  }
}