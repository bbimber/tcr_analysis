library(lattice)
library(gplots)
library(ggplot2)
library(reshape)
#Note: if loaded, this must load prior to dplyr.
#library(plyr)
library(dplyr)
library(Rlabkey)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(data.table)

library(dtplyr)

pullTCRMetaFromLabKey <- function(requireGeneFile = FALSE, filterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE){
  df <- rbind(pullFullTranscriptMetaFromLabKey(filterClauses = filterClauses, replicateAsSuffix = replicateAsSuffix), pullEnrichedMetaFromLabKey(filterClauses = filterClauses, separateEnriched = separateEnriched, replicateAsSuffix = replicateAsSuffix, skipEnriched = skipEnriched))
  
  if (requireGeneFile){
    df <- df[!is.na(df$OutputFileId),]  
  }
  
  df <- df[is.na(df$Status),]
  
  df$ResponseType <- c('Other')
  df$ResponseType[df$Peptide %in% c('Gag120', 'Gag69')] <- c('MHC-E Supertope')
  df$ResponseType[df$Peptide %in% c('Gag16', 'Gag18', 'Gag23', 'Gag30', 'Gag33', 'Gag50', 'Gag65', 'Gag97', 'Gag109', 'Gag119')] <- c('Other E-Restricted')
  df$ResponseType[df$Peptide %in% c('Gag27/28', 'Gag37')] <- c('Other Class II-Restricted')
  df$ResponseType[df$Peptide %in% c('RL8', 'RL9', 'RL10', 'QI9', 'LF8')] <- c('Conventional')
  df$ResponseType <- as.factor(df$ResponseType)
  
  df <- df %>% 
    group_by(SeqDataType, StimId, Population, Replicate, IsSingleCell) %>% 
    mutate(totalCellsForGroup = sum(Cells))
  
  return(df)
}

reLevelFactor <- function(f, ordered){
  for (l in rev(ordered)){
    if (l %in% levels(f)){
      f <- relevel(f, l)
    } else {
      print(paste0('unknown level: ', l))
    }
  }
  
  return(f)
}

pullEnrichedMetaFromLabKey <- function(filterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE){
  clauses <- list(
    c('enrichedReadsetId', 'NON_BLANK', ''), 
    c('enrichedReadsetId/totalFiles', 'GT', 0)
  )
  
  if (!is.null(filterClauses)){
    clauses = append(clauses, filterClauses)
  }
  
  filter <- do.call(makeFilter, clauses)
  if (length(clauses) > 0){
    print('TCR enriched meta filter')
    print(filter)
  }

  df <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Internal/Bimber/", 
    schemaName="tcrdb", 
    queryName="cdnas", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c(
      'enrichedReadsetId','enrichedReadsetId/status','enrichedReadsetId/workbook','enrichedReadsetId/totalForwardReads','enrichedReadsetId/numCDR3s','enrichedReadsetId/distinctLoci','enrichedReadsetId/numTcrRuns',
      'sortId', 'sortId/stimId','sortId/stimId/animalId','sortId/stimId/date','sortId/stimId/stim','sortId/stimId/treatment','sortId/stimId/activated','sortId/stimId/background',
      'sortId/population','sortId/replicate','sortId/cells', "rowid"
    ),
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=filter
  )
  
  df <- doSharedColumnRename(df)
  
  names(df)[names(df)=="enrichedreadsetid"] <- "ReadsetId"
  df$ReadsetId <- as.integer(df$ReadsetId)
  
  names(df)[names(df)=="enrichedreadsetid_totalforwardreads"] <- "TotalReads"
  df$TotalReads <- as.integer(df$TotalReads)
  
  names(df)[names(df)=="enrichedreadsetid_numcdr3s"] <- "NumCDR3s"
  df$NumCDR3s <- as.integer(df$NumCDR3s)
  
  names(df)[names(df)=="enrichedreadsetid_distinctloci"] <- "DistinctLoci"
  df$DistinctLoci <- as.factor(df$DistinctLoci)
  
  names(df)[names(df)=="enrichedreadsetid_numtcrruns"] <- "RunTCRRuns"
  df$RunTCRRuns <- as.integer(df$RunTCRRuns)
  
  names(df)[names(df)=="enrichedreadsetid_workbook"] <- "WorkbookId"
  df$WorkbookId <- as.factor(df$WorkbookId)
  
  names(df)[names(df)=="enrichedreadsetid_status"] <- "Status"
  df$Status <- as.factor(df$Status)
  
  df$SeqDataType <- 'TCR_Enriched'
  df$SeqDataType <- as.factor(df$SeqDataType)
  
  df$Label <- as.character(df$Peptide)
  
  if (replicateAsSuffix){
    df$Label[!is.na(df$Replicate)] <- paste0(df$Label[!is.na(df$Replicate)], '_', df$Replicate[!is.na(df$Replicate)])  
  }
  
  df$Label[df$IsSingleCell] <- paste0(df$Label[df$IsSingleCell], '**')
  if (separateEnriched){
    df$Label <- paste0(df$Label, "+")
    
  }
  df$Label <- as.factor(df$Label)
  
  print(paste0('total TCR enriched rows: ' , nrow(df)))
  if (skipEnriched){
    print('dropping all TCR enriched rows')
    df <- df[c(FALSE),]  
  }

  
  return(df)
}

pullFullTranscriptMetaFromLabKey <- function(filterClauses = NULL, replicateAsSuffix = TRUE){
  clauses = list(
    c('readsetId', 'NON_BLANK', ''), 
    c('readsetId/totalFiles', 'GT', 0)
  )
  
  if (!is.null(filterClauses)){
    clauses = append(clauses, filterClauses)
  }
  
  filter <- do.call(makeFilter, clauses)
  if (length(clauses) > 0){
    print('Full transcript meta filter')
    print(filter)
  }
  
  df <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Internal/Bimber/", 
    schemaName="tcrdb", 
    queryName="cdnas", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c(
      'readsetId','readsetId/status','readsetId/workbook','readsetId/totalForwardReads','readsetId/numCDR3s','readsetId/distinctLoci','readsetId/numTcrRuns',
      'sortId', 'sortId/stimId','sortId/stimId/animalId','sortId/stimId/date','sortId/stimId/stim','sortId/stimId/treatment','sortId/stimId/activated','sortId/stimId/background',
      'sortId/population','sortId/replicate','sortId/cells', "rowid"
    ),
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=filter
  )
  
  df <- doSharedColumnRename(df)
  
  names(df)[names(df)=="readsetid"] <- "ReadsetId"
  df$ReadsetId <- as.integer(df$ReadsetId)
  
  names(df)[names(df)=="readsetid_totalforwardreads"] <- "TotalReads"
  df$TotalReads <- as.integer(df$TotalReads)
  
  names(df)[names(df)=="readsetid_numcdr3s"] <- "NumCDR3s"
  df$NumCDR3s <- as.integer(df$NumCDR3s)
  
  names(df)[names(df)=="readsetid_distinctloci"] <- "DistinctLoci"
  df$DistinctLoci <- as.factor(df$DistinctLoci)
  
  names(df)[names(df)=="readsetid_numtcrruns"] <- "RunTCRRuns"
  df$RunTCRRuns <- as.integer(df$RunTCRRuns)
  
  names(df)[names(df)=="readsetid_workbook"] <- "WorkbookId"
  df$WorkbookId <- as.factor(df$WorkbookId)
  
  names(df)[names(df)=="readsetid_status"] <- "Status"
  df$Status <- as.factor(df$Status)
  
  df$SeqDataType <- 'Full_Transcriptome'
  df$SeqDataType <- as.factor(df$SeqDataType)
  
  df$Label <- as.character(df$Peptide)
  
  if (replicateAsSuffix){
    df$Label[!is.na(df$Replicate)] <- paste0(df$Label[!is.na(df$Replicate)], '_', df$Replicate[!is.na(df$Replicate)])  
  }
  
  df$Label[df$IsSingleCell] <- paste0(df$Label[df$IsSingleCell], '**')
  
  df$Label <- as.factor(df$Label)
  
  print(paste0('total full transcriptome rows: ' , nrow(df)))
  
  return(df)
}

doSharedColumnRename <- function(df){
  names(df)[names(df)=="sortid_stimid_animalid"] <- "AnimalId"
  df$AnimalId <- as.factor(df$AnimalId)
  
  names(df)[names(df)=="sortid_stimid_date"] <- "SampleDate"
  names(df)[names(df)=="rowid"] <- "DatasetId"
  
  names(df)[names(df)=="sortid_cells"] <- "Cells"
  df$Cells <- as.integer(df$Cells)
  
  names(df)[names(df)=="sortid_replicate"] <- "Replicate"
  df$Replicate <- as.factor(df$Replicate)
  
  names(df)[names(df)=="sortid"] <- "SortId"
  df$SortId <- as.factor(df$SortId)
  
  df$IsSingleCell <- df['Cells'] == 1
  
  #TODO: consider if this is the best approach
  df$Replicate[df$IsSingleCell] <- c(NA)
  
  names(df)[names(df)=="sortid_stimid_treatment"] <- "Treatment"
  df$Treatment <- as.factor(df$Treatment)
  
  names(df)[names(df)=="sortid_stimid_stim"] <- "Peptide"
  df$Peptide <- as.factor(df$Peptide)
  
  #names(df)[names(df)=="cellclass"] <- "Cellclass"
  #df$Cellclass <- as.factor(df$Cellclass)
  
  names(df)[names(df)=="sortid_population"] <- "Population"
  df$Population <- as.factor(df$Population)
  
  names(df)[names(df)=="sortid_stimid"] <- "StimId"
  df$StimId <- as.factor(df$StimId)
  
  return(df)
}

pullTCRResultsFromLabKey <- function(colFilterClauses = NULL){
  colFilter = NULL
  if (!is.null(colFilterClauses)){
    colFilter = do.call(makeFilter, colFilterClauses)
  }
  
  if (!is.null(colFilter)){
    print('assay result filter')
    print(colFilter)

  }
  
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
  
  print(paste0('total assay rows (unfiltered): ', nrow(df)))
  
  df <- df[!is.na(df$cdr3),]
  df <- df[is.na(df$disabled),]
  df$LocusCDR3 <- paste0(df$locus, ":", df$cdr3)
  df$LocusCDR3 <- as.factor(df$LocusCDR3)
  
  df$LocusCDR3WithUsage <- paste0(df$locus, ":", df$cdr3,":",df$vgene,"-",df$jgene)
  df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  
  dateFormat <- date_format(format ="%Y-%m-%d")
  df$dateFormatted <- dateFormat(df$date)
  
  df$locus <- as.factor(df$locus)
  
  print(paste0('total assay rows: ', nrow(df)))
  
  return(df)
}

pullResultsFromLabKey <- function(resultFilterClauses = NULL, metaFilterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE){
  results <- pullTCRResultsFromLabKey(colFilterClauses = resultFilterClauses)
  meta <- pullTCRMetaFromLabKey(requireGeneFile = FALSE, filterClauses = metaFilterClauses, separateEnriched = separateEnriched, replicateAsSuffix = replicateAsSuffix, skipEnriched = skipEnriched)
  meta$Label <- reorder(meta$Label)
  
  results <- merge(results, meta, by.x=c('analysisid_readset'), by.y=c('ReadsetId')) #, all.x=FALSE, all.y=FALSE
  print(paste0('total merged rows (prefilter): ', nrow(results)))

  print(paste0('total merged rows: ', nrow(results)))
  
  results <- results %>% 
    group_by(SeqDataType, StimId, Population, Replicate, IsSingleCell, locus) %>% 
    mutate(totalCellsForGroupAndLocus = sum(Cells))
  
  print(paste0('total merged rows after group: ', nrow(results)))
  
  results <- qualityFilterResults(results)
  results <- groupSingleCellData(results)
  results <- filterBasedOnCellCount(results)
  
  return(list(results=results, meta = meta))
}

qualityFilterResults <- function(df = NULL, minLibrarySize=2000000, minReadCount=500000, minReadCountForEnriched=5000){
  print(paste0('performing quality filter, initial rows: ', nrow(df)))
  
  r <- nrow(df)
  #TODO: restore this
  #df <- df[df$EstimatedLibrarySize > minLibrarySize,]
  #print(paste0('rows dropped for library size: ', (r - nrow(df))))

  r <- nrow(df)
  #toDrop <- unique(df[!(df$SeqDataType == 'Full_Transcriptome' & df$TotalReads >= minReadCount),c('AnimalId', 'SampleDate', 'Peptide', 'Population', 'TotalReads')])
  #write.table(toDrop, file = 'toDrop.txt', sep = '\t', row.names = FALSE, quote = FALSE)
  
  df <- df[!(df$SeqDataType == 'Full_Transcriptome' & df$TotalReads < minReadCount),]
  print(paste0('total RNA rows dropped for read count: ', (r - nrow(df)), ', remaining: ', nrow(df)))

  r <- nrow(df)
  df <- df[!(df$SeqDataType == 'TCR_Enriched' & df$TotalReads < minReadCountForEnriched),]
  print(paste0('total TCR enriched rows dropped for read count: ', (r - nrow(df)), ', remaining: ', nrow(df)))
  
  return(df)
}

groupSingleCellData <- function(df = NULL, lowFreqThreshold = 0.05, minTcrReads = 50){
  #print(nrow(df))
  
  df <- df %>% 
    group_by(SeqDataType, StimId, Population, Replicate, IsSingleCell) %>% 
    mutate(totalTcrReadsForGroup = sum(count))
  
  df <- df %>% 
    group_by(SeqDataType, StimId, Population, Replicate, IsSingleCell, locus) %>% 
    mutate(totalTcrReadsForGroupAndLocus = sum(count))

  #bulk
  dfb <- df[!df$IsSingleCell,]
  dfb$percentage = dfb$count / dfb$totalTcrReadsForGroup
  dfb$percentageForLocus = dfb$count / dfb$totalTcrReadsForGroupAndLocus
  dfb <- dfb[dfb$totalTcrReadsForGroup >= minTcrReads,]
  dfb$totalCellsForCDR3 <- c(0)
  columns <- c('SeqDataType', 'StimId', 'Population', 'IsSingleCell', 'Label', 'AnimalId', 'SampleDate', 'Peptide', 'Treatment', 'locus', 'vhit', 'jhit', 'cdr3', 'count', 'date', 'dateFormatted', 'LocusCDR3', 'LocusCDR3WithUsage', 'Replicate', 'totalCellsForGroup', 'totalCellsForGroupAndLocus', 'Cells', 'percentage', 'percentageForLocus', 'totalCellsForCDR3')
  dfb <- dfb[,columns]
  print(paste0('total bulk rows: ', nrow(dfb)))
  
  #single cell
  dfs <- df[df$IsSingleCell,]
  print(paste0('starting single cell rows: ', nrow(dfs)))
  if (nrow(dfs) > 0){
    dfs <- dfs %>%
      group_by(SeqDataType, StimId, Population, IsSingleCell, Label, AnimalId, SampleDate, Peptide, Treatment, locus, vhit, jhit, cdr3, date, dateFormatted, LocusCDR3, LocusCDR3WithUsage, Replicate, totalCellsForGroup, totalCellsForGroupAndLocus) %>%
      summarize(totalCellsForCDR3 = sum(Cells))
    dfs$Cells <- c(1)
    dfs$count <- dfs$totalCellsForCDR3
    
    dfs$percentage <- dfs$totalCellsForCDR3 / dfs$totalCellsForGroup
    dfs$percentageForLocus <- dfs$totalCellsForCDR3 / dfs$totalCellsForGroupAndLocus
  } else {
    dfs$totalCellsForCDR3 <- numeric()
    dfs$Cells <- integer()
    dfs$count <- dfs$totalCellsForCDR3
    dfs$percentage <- numeric()
    dfs$percentageForLocus <- numeric()
  }
  
  print(paste0('total single cell rows: ', nrow(dfs)))
  
  df <- rbind(dfb[columns], dfs[columns])
  print(paste0('combined: ', nrow(df)))
  
  df <- df %>% group_by(LocusCDR3, AnimalId) %>% mutate(maxPercentageInAnimal = max(percentage)) 
  df <- df %>% filter(Peptide != 'SEB') %>% group_by(LocusCDR3, AnimalId, Population) %>% mutate(maxPercentageInAnimalPopulation = max(percentage)) 
  df <- df %>% group_by(LocusCDR3, Peptide) %>% mutate(maxPercentageInPeptide = max(percentage)) 
  
  if (lowFreqThreshold){
    df$LocusCDR3 <- as.character(df$LocusCDR3)
    df$LocusCDR3[df$maxPercentageInAnimal < lowFreqThreshold] <- c('Low Freq.')
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    
    df$LocusCDR3WithUsage <- as.character(df$LocusCDR3WithUsage)
    df$LocusCDR3WithUsage[df$maxPercentageInAnimal < lowFreqThreshold] <- c('Low Freq.')
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  }
  
  df <- df %>% group_by(StimId, Population, Replicate, IsSingleCell) %>% mutate(totalTcrReadsForGroupAfterLowFreq = sum(count)) 
  
  #TODO
  df$percentageAfterLowFreq <- df$count / df$totalTcrReadsForGroupAfterLowFreq
  
  df <- arrange(df, LocusCDR3)
  
  return(df)
}

filterBasedOnCellCount <- function(df = NULL){
  startRows <- nrow(df)
  df$minPctToInclude <- c(0)
  
  if (length(df$IsSingleCell == TRUE) > 0){
    df$minPctToInclude[df$IsSingleCell] <- c(0.1)
    df$minPctToInclude[df$IsSingleCell] <- (1 / df$Cells[df$IsSingleCell]) 
  }

  print(paste0('rows filtered by cell count: ', (startRows - nrow(df))))
  
  return(df)
}

makePlot <- function(locus, df, pal = 'Reds', bulkOnly = FALSE, tnfPosOnly = FALSE, showLegend = FALSE, yMax=1.05, yField='percentage', xField='Label', fillField='LocusCDR3', facet1 = 'Population', facet2 = 'dateFormatted', doFacet=TRUE){
  locusData <- df[df$locus == locus,]
  
  if (bulkOnly){
    locusData <- locusData[!locusData$IsSingleCell,]
  }
  
  if (tnfPosOnly){
    locusData <- locusData[locusData$Population == 'TNF-Pos',]
  }
  
  locusData[fillField] <- as.factor(as.character(locusData[[fillField]]))
  totalCDR3 = length(unique(locusData[[fillField]]))
  print(paste(locus, ', total CDR3', totalCDR3))
  
  if (nrow(locusData) == 0){
    return(ggplot(data.frame()) +geom_blank())
  }
  
  getPalette = colorRampPalette(brewer.pal(9, pal))
  
  colorsDT <-  data.table(fillField=levels(locusData[[fillField]]), Color=getPalette(totalCDR3))
  locusData <- data.table(locusData)
  locusData <- merge(locusData, colorsDT, by.x = c(fillField), by.y=c('fillField'), all = TRUE)
  locusData$Color[locusData[[fillField]] == 'Low Freq.'] <- c('#D3D3D3')
  locusData$Color <- as.factor(locusData$Color)
  
  colorValues <- getPalette(totalCDR3)
  #colorValues <- rev(colorValues)
  colorValues[[1]] <- '#D3D3D3'  #grey
  
  facetFormula <- as.formula(paste0(facet1, ' ~ ', facet2))
  
  N <- locusData[, c(facet2, xField), with = FALSE]
  names(N) <- c('GroupField', 'CountField')

  N <- N %>% 
    group_by(GroupField) %>% 
    summarise(V1 = n_distinct(CountField))
  
  names(N) <- c(facet2, 'V1')
  
  N$scaleFactor = N$V1 / max(N$V1)
  
  locusData = merge(locusData, N, by.x = c(facet2), by.y = c(facet2), all.x = TRUE)
  if (showLegend){
    lp = 'right'
  }
  else {
    lp = 'none'
  }
  
  P1<-ggplot(locusData) +
    aes_string(x = xField, y = yField, order = 'percentage', width = '0.5*scaleFactor') +
    geom_bar(aes_string(fill = fillField), stat = 'identity', colour="black") +
    scale_fill_manual(values = colorValues) +
    theme_bw() +   #base_family = "Arial-Black"
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
    ylab('')
  
  if (doFacet){
    P1 <- P1 + facet_grid(facetFormula, scales = 'free_x')
  }
  
  return(P1)
}

makePlotGroup <- function(df, subDir, basename, tnfPosOnly=FALSE, includeTRB=FALSE, includeTRG=FALSE, pal='Reds', showLegend=FALSE, bulkOnly=FALSE, yMax=NA, yField='percentage', xField='Label', fillField='LocusCDR3', replicateAsSuffix = TRUE, allLocusWidth=1800, singleLocusWidth=1000, asEPS = FALSE, doFacet=TRUE){
  plots = list()
  i = 1
  
  facet1 <- 'Population'
  facet2 <- 'dateFormatted'

  for (locus in sort(unique(df$locus))){
    P1 <- makePlot(locus, df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly, yField=yField, xField=xField, fillField=fillField, facet1=facet1, facet2=facet2, doFacet=doFacet)
    plots[[i]] <- P1
    i=i+1
  }
  
  suffix = ''
  if (bulkOnly){
    suffix = paste0(suffix, '.bulk')
  }
  
  if (tnfPosOnly){
    suffix = paste0(suffix, '.tnf')
  }
  
  if (!file.exists(subDir)){
    dir.create(file.path('./', subDir))
  }
  
  fn <- paste0('./',subDir,'/', basename, suffix)
  
  write.table(df, file=paste0(fn, '.txt'), sep='\t', row.names = FALSE)
  
  if (asEPS) {
    postscript(paste0(fn,'.eps'), height = 1000, width = allLocusWidth)  
  } else {
    png(paste0(fn,'.png'), height = 1000, width = allLocusWidth)  
  }
  
  do.call(grid.arrange, plots, c(ncol=2))
  dev.off()
  
  if (includeTRB){
    if (asEPS) {
      postscript(paste0(fn,'.TRB.eps'), height = 500, width = singleLocusWidth)
    } else {
      png(paste0(fn,'.TRB.png'), height = 500, width = singleLocusWidth)
    }
    P_B <- makePlot('TRB', df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly,yField='percentageForLocus', xField=xField, fillField=fillField, facet1=facet1, facet2=facet2, doFacet=doFacet)
    P_B = P_B + labs(title = paste0(basename, " TRB"))
    print(P_B)
    dev.off()
  }
  
  if (includeTRG){
    if (asEPS) {
      postscript(paste0(fn,'.TRG.eps'), height = 500, width = singleLocusWidth)
    } else {
      png(paste0(fn,'.TRG.png'), height = 500, width = singleLocusWidth)
    }
    P_G <- makePlot('TRG', df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly,yField='percentageForLocus', xField=xField, fillField=fillField, facet1=facet1, facet2=facet2, doFacet=doFacet)
    P_G = P_G + labs(title = paste0(basename, " TRG"))
    print(P_G)
    dev.off() 
  }
  
  if (includeTRB){
    return(P_B)
  }
}