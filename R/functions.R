library(lattice)
library(gplots)
library(ggplot2)
library(data.table)
#Note: if loaded, this must load prior to dplyr.
#library(plyr)
library(dplyr)
library(Rlabkey)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(naturalsort)
library(stringi)
library(plotly)

library(dtplyr)

pullTCRMetaFromLabKey <- function(requireGeneFile = FALSE, filterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE, collapseReplicates = FALSE, savePrefix, overwrite = T, folderPath = "/Labs/Bimber/"){
  saveFile <- paste0(savePrefix, '.metadata.txt')
  if (!overwrite && file.exists(saveFile)) {
    print('reading metadata from file')
    df <- read.table(saveFile, sep = '\t', header = T)
  } else {
    dfT <- pullFullTranscriptomeMetaFromLabKey(filterClauses = filterClauses, replicateAsSuffix = replicateAsSuffix, collapseReplicates = collapseReplicates, folderPath = folderPath)
    dfE <- pullEnrichedMetaFromLabKey(filterClauses = filterClauses, separateEnriched = separateEnriched, replicateAsSuffix = replicateAsSuffix, collapseReplicates = collapseReplicates, skipEnriched = skipEnriched, folderPath = folderPath)
    
    if (nrow(dfT) > 0 && nrow(dfE) > 0) {
      df <- rbind(dfT, dfE)
    } else if (nrow(dfT) > 0) {
      df <- dfT  
    } else if (nrow(dfE) > 0) {
      df <- dfE
    }
    write.table(df, file = saveFile, row.names = F, quote = F, sep = '\t')
  }
  
  if (requireGeneFile){
    df <- df[!is.na(df$OutputFileId),]  
  }
  
  toDrop <- sum(!is.na(df$Status))
  if (toDrop > 0) {
    print(paste0('Total rows dropped due to readset status: ', toDrop))
    df <- df[is.na(df$Status),]
  }
  
  toDrop <- sum(!is.na(df$cDNA_Status))
  if (toDrop > 0) {
    print(paste0('Total rows dropped due to cDNA status: ', toDrop))
    print(paste0('status: ', paste0(unique(df$cDNA_Status), collapse = ';')))
    df <- df[is.na(df$cDNA_Status),]
  }
  
  if (nrow(df) > 0) {
    df$ResponseType <- c('Other')
    df$ResponseType[df$Peptide %in% c('Gag120', 'Gag69')] <- c('MHC-E Supertope')
    df$ResponseType[df$Peptide %in% c('Gag53', 'Gag73')] <- c('MHC-II Supertope')
    df$ResponseType[df$Peptide %in% c('Gag16', 'Gag18', 'Gag23', 'Gag30', 'Gag33', 'Gag50', 'Gag65', 'Gag97', 'Gag109', 'Gag119')] <- c('Other E-Restricted')
    df$ResponseType[df$Peptide %in% c('Gag27/28', 'Gag37')] <- c('Other Class II-Restricted')
    df$ResponseType[df$Peptide %in% c('RL8', 'RL9', 'RL10', 'QI9', 'LF8')] <- c('Conventional')
    df$ResponseType <- as.factor(df$ResponseType)
  } else {
    df$ResponseType <- character()
  }
  
  df <- df %>% 
    group_by(SeqDataType, SampleId, Population, Replicate, IsSingleCell) %>% 
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

pullEnrichedMetaFromLabKey <- function(filterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE, collapseReplicates = FALSE, folderPath = "/Labs/Bimber/"){
  print('downloading metadata for TCR enriched datasets')
  clauses <- list(
    c('tcrreadsetid', 'NON_BLANK', ''), 
    c('tcrreadsetid/totalFiles', 'GT', 0)
  )
  
  if (!is.null(filterClauses) && !is.na(filterClauses)){
    clauses = append(clauses, filterClauses)
  }
  
  filter <- do.call(makeFilter, clauses)
  if (length(clauses) > 0){
    print('TCR enriched meta filter')
    print(filter)
  }

  df <- suppressWarnings(labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath=folderPath, 
    schemaName="singlecell", 
    queryName="cdna_libraries", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c(
      'tcrreadsetid','tcrreadsetid/status','tcrreadsetid/workbook','tcrreadsetid/totalForwardReads','tcrreadsetid/numCDR3s','tcrreadsetid/distinctLoci','tcrreadsetid/numTcrRuns',
      'sortId', 'sortId/sampleId','sortId/sampleId/subjectId','sortId/sampleId/sampledate','sortId/sampleId/stim','sortId/sampleId/assayType',
      'sortId/population','sortId/replicate','sortId/cells', "rowid", "status", "tcrreadsetid/librarytype"
    ),
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=filter
  ))
  
  df <- doSharedColumnRename(df)
  
  names(df)[names(df)=="tcrreadsetid"] <- "ReadsetId"
  df$ReadsetId <- as.integer(df$ReadsetId)
  
  names(df)[names(df)=="tcrreadsetid_totalforwardreads"] <- "TotalReads"
  df$TotalReads <- as.integer(df$TotalReads)
  
  names(df)[names(df)=="tcrreadsetid_numcdr3s"] <- "NumCDR3s"
  df$NumCDR3s <- as.integer(df$NumCDR3s)
  
  names(df)[names(df)=="tcrreadsetid_distinctloci"] <- "DistinctLoci"
  df$DistinctLoci <- as.factor(df$DistinctLoci)
  
  names(df)[names(df)=="tcrreadsetid_numtcrruns"] <- "RunTCRRuns"
  df$RunTCRRuns <- as.integer(df$RunTCRRuns)
  
  names(df)[names(df)=="tcrreadsetid_workbook"] <- "WorkbookId"
  df$WorkbookId <- as.factor(df$WorkbookId)
  
  names(df)[names(df)=="status"] <- "cDNA_Status"
  df$cDNA_Status <- as.factor(df$cDNA_Status)
  
  names(df)[names(df)=="tcrreadsetid_status"] <- "Status"
  df$Status <- as.factor(df$Status)

  if (nrow(df) == 0){
    df$SeqDataType <- character()
  } else {
    df$SeqDataType <- c('TCR_Enriched')
    df$SeqDataType[grep(x = df$tcrreadsetid_librarytype, pattern = 'VDJ')] <- c('10x_VDJ')
    
    df$IsSingleCell[df$SeqDataType == '10x_VDJ'] <- F
  }
  
  df$Label <- as.character(df$Peptide)
  
  df <- df[!(names(df) %in% c('tcrreadsetid_librarytype'))]
  
  if (nrow(df)> 0 && collapseReplicates){
    print('collapsing replicates')
    replicateAsSuffix <- FALSE
    df$Replicate <- c(NA)
  } else if (replicateAsSuffix){
    print('adding replicate to label')
    df$Label[!is.na(df$Replicate)] <- paste0(df$Label[!is.na(df$Replicate)], '_', df$Replicate[!is.na(df$Replicate)])  
  }
  
  if (nrow(df) > 0) {
    df$Label[df$IsSingleCell] <- paste0(df$Label[df$IsSingleCell], '**')
    if (separateEnriched){
      df$Label <- paste0(df$Label, "+")
      
    }
  }
  
  df$Label <- naturalfactor(df$Label)
  
  print(paste0('total TCR enriched rows: ' , nrow(df)))
  print(paste0('total TCR enriched rows (non-10x): ' , sum(df$SeqDataType == 'TCR_Enriched')))
  print(paste0('total TCR enriched rows (10x): ' , sum(df$SeqDataType == '10x_VDJ')))
  if (skipEnriched){
    print('dropping all TCR enriched rows')
    df <- df[df$SeqDataType == 'TCR_Enriched',]  
  }

  return(df)
}

pullFullTranscriptomeMetaFromLabKey <- function(filterClauses = NULL, replicateAsSuffix = TRUE, collapseReplicates = FALSE, folderPath = "/Labs/Bimber/"){
  print('downloading metadata for whole transcriptome datasets')
  clauses = list(
    c('readsetId', 'NON_BLANK', ''), 
    c('readsetId/totalFiles', 'GT', 0),
    c('readsetId/librarytype', 'NOT_EQUAL_OR_NULL', '10x 5\' GEX')
  )
  
  if (!is.null(filterClauses) && !is.na(filterClauses)){
    clauses = append(clauses, filterClauses)
  }
  
  filter <- do.call(makeFilter, clauses)
  if (length(clauses) > 0){
    print('Full transcript meta filter')
    print(filter)
  }
  
  df <- suppressWarnings(labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath=folderPath, 
    schemaName="singlecell", 
    queryName="cdna_libraries", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c(
      'readsetId','readsetId/status','readsetId/workbook','readsetId/totalForwardReads','readsetId/numCDR3s','readsetId/distinctLoci','readsetId/numTcrRuns',
      'sortId', 'sortId/sampleId','sortId/sampleId/subjectId','sortId/sampleId/sampledate','sortId/sampleId/stim','sortId/sampleId/assayType',
      'sortId/population','sortId/replicate','sortId/cells', "rowid", "status"
    ),
    containerFilter=NULL,
    colNameOpt='rname',
    colFilter=filter
  ))
  
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
  
  names(df)[names(df)=="status"] <- "cDNA_Status"
  df$cDNA_Status <- as.factor(df$cDNA_Status)
  
  names(df)[names(df)=="readsetid_status"] <- "Status"
  df$Status <- as.factor(df$Status)
  
  if (nrow(df) > 0){
    df$SeqDataType <- 'Full_Transcriptome'
    df$SeqDataType <- as.factor(df$SeqDataType)
  }
  
  df$Peptide <- naturalfactor(df$Peptide)
  
  df$Label <- as.character(df$Peptide)
  
  if (collapseReplicates){
    replicateAsSuffix <- FALSE
    df$Replicate <- c(NA)
  } else if (replicateAsSuffix){
    df$Label[!is.na(df$Replicate)] <- paste0(df$Label[!is.na(df$Replicate)], '_', df$Replicate[!is.na(df$Replicate)])  
  }
  
  df$Label[df$IsSingleCell] <- paste0(df$Label[df$IsSingleCell], '**')
  
  df$Label <- naturalfactor(df$Label)
  
  print(paste0('total full transcriptome rows: ' , nrow(df)))
  
  return(df)
}

doSharedColumnRename <- function(df){
  names(df)[names(df)=="sortid_sampleid_subjectid"] <- "SubjectId"
  df$SubjectId <- as.factor(df$SubjectId)
  
  names(df)[names(df)=="sortid_sampleid_sampledate"] <- "SampleDate"
  names(df)[names(df)=="rowid"] <- "DatasetId"
  
  names(df)[names(df)=="sortid_cells"] <- "Cells"
  df$Cells <- as.integer(df$Cells)
  
  names(df)[names(df)=="sortid_replicate"] <- "Replicate"
  df$Replicate <- as.character(df$Replicate)
  
  names(df)[names(df)=="sortid"] <- "SortId"
  df$SortId <- as.integer(df$SortId)
  
  if (nrow(df) == 0) {
    df$IsSingleCell <- logical()
  } else {
    df$IsSingleCell <- df['Cells'] == 1
  }
  
  if (sum(is.na(df$IsSingleCell)) > 0 ){
    stop('Some rows are NA for IsSingleCell!')
  }
  
  #TODO: consider if this is the best approach
  df$Replicate[df$IsSingleCell] <- c(NA)
  
  names(df)[names(df)=="sortid_sampleid_assaytype"] <- "AssayType"
  df$AssayType <- as.factor(df$AssayType)
  
  names(df)[names(df)=="sortid_sampleid_stim"] <- "Peptide"
  df$Peptide <- naturalfactor(df$Peptide)
  
  #names(df)[names(df)=="cellclass"] <- "Cellclass"
  #df$Cellclass <- as.factor(df$Cellclass)
  
  names(df)[names(df)=="sortid_population"] <- "Population"
  df$Population <- naturalfactor(df$Population)
  
  names(df)[names(df)=="sortid_sampleid"] <- "SampleId"
  df$SampleId <- as.integer(df$SampleId)
  
  return(df)
}

pullTCRResultsFromLabKey <- function(colFilterClauses = NULL, savePrefix, overwrite = T, folderPath = "/Labs/Bimber/"){
  colFilter = NULL
  if (!is.null(colFilterClauses) && !is.na(colFilterClauses)){
    colFilter = do.call(makeFilter, colFilterClauses)
  }
  
  saveFile <- paste0(savePrefix, '.raw.results.txt')
  if (!overwrite && file.exists(saveFile)) {
    print('reading assay data from file')
    df <- read.table(saveFile, sep = '\t', header = T)
  } else {
    df <- suppressWarnings(labkey.selectRows(
      baseUrl="https://prime-seq.ohsu.edu", 
      folderPath=folderPath, 
      schemaName="assay.TCRdb.TCRdb", 
      queryName="Data", 
      viewName="", 
      colSelect=c('sampleName', 'cdna', 'calculatedPopulation', 'libraryId','libraryId/label','analysisId', 'analysisId/readset', 'locus','vHit','dHit','jHit','cHit','vGene', 'jGene','CDR3','count','fraction', 'disabled'), 
      containerFilter=NULL,
      colNameOpt='rname',
      colFilter=colFilter
    ))
    
    write.table(df, file = saveFile, row.names = F, quote = F, sep = '\t')
  }
  
  print(paste0('total assay rows (unfiltered): ', nrow(df)))

    #Drop rows without results, which can be stored in the DB to track failed samples 
  df <- df[!is.na(df$cdr3),]
  df <- df[is.na(df$disabled),]
  
  if (sum(is.na(df$locus)) > 0) {
    stop('there are rows without locus')
  }
  
  if (nrow(df) > 0) {
    df$LocusCDR3 <- paste0(df$locus, ":", df$cdr3)
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    
    df$LocusCDR3WithUsage <- paste0(df$locus, ":", df$cdr3,":",df$vgene,"-",df$jgene)
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)

    df$locus <- as.factor(df$locus)
  } else {
    df$LocusCDR3 <- character()
    df$LocusCDR3WithUsage <- character()
  }
  
  print(paste0('total assay rows: ', nrow(df)))
  
  return(df)
}

pullResultsFromLabKey <- function(resultFilterClauses = NULL, metaFilterClauses = NULL, separateEnriched = FALSE, replicateAsSuffix = TRUE, skipEnriched = FALSE, collapseReplicates = FALSE, lowFreqThreshold = 0.05, backgroundThreshold = 0.5, savePrefix = NULL, overwrite = T, folderPath = "/Labs/Bimber/"){
  if (is.null(savePrefix)) {
    stop('Must supply a save prefix')
  }
  
  results <- pullTCRResultsFromLabKey(colFilterClauses = resultFilterClauses, savePrefix = savePrefix, overwrite = overwrite, folderPath = folderPath)
  origResults <- nrow(results)

  meta <- pullTCRMetaFromLabKey(requireGeneFile = FALSE, filterClauses = metaFilterClauses, separateEnriched = separateEnriched, replicateAsSuffix = replicateAsSuffix, skipEnriched = skipEnriched, collapseReplicates = collapseReplicates, savePrefix = savePrefix, overwrite = overwrite, folderPath = folderPath)
  
  if (nrow(results[is.na(results$cdna),]) > 0) {
    stop('Exptected all result rows to have cDNA!')
  }
  
  #Note: because cell hashing requires a one-to-many readset to cDNA relationship and b/c legacy data did not store the cDNA ID in the assay, split this:
  results <- merge(results[!is.na(results$cdna),], meta, by.x=c('cdna'), by.y=c('DatasetId')) 
  results <- results[ , -which(names(results) %in% c("ReadsetId"))]
  names(results)[names(results) == 'cdna'] <- 'DatasetId'
  print(paste0('Rows with metadata identified by cDNA: ', nrow(results)))

  if (nrow(results) == 0) {
    warning('No results found')
    return(NULL)
  }
  
  print(paste0('total merged rows (prefilter): ', nrow(results)))

  # This is to allow TCR results to be imported for subsets, such as based on transcription:
  if (sum(!is.na(results$calculatedpopulation)) > 0) {
    results$Population <- as.character(results$Population)
    results$Population[!is.na(results$calculatedpopulation)] <- as.character(results$calculatedpopulation[!is.na(results$calculatedpopulation)])
    results$Population <- as.factor(results$Population)
  }
  
  # This will occur if the filters supplied for  results and metadata are not the same
  if (origResults != nrow(results)){
    print(paste0('Total results are equal before/after merge: ', origResults, '/', nrow(results)))
  }
  
  print(paste0('total merged rows: ', nrow(results)))
  
  results <- qualityFilterResults(results, savePrefix = savePrefix)
  print(paste0('total rows after quality filter: ', nrow(results)))
  
  results <- results %>% 
    group_by(SeqDataType, SampleId, Population, Replicate, IsSingleCell, locus) %>% 
    mutate(totalCellsForGroupAndLocusBulk = sum(Cells), totalCellsForGroupAndLocusSS = n_distinct(analysisid))
  
  results <- results %>% 
    group_by(SeqDataType, SampleId, Population, Replicate, IsSingleCell, locus, analysisid) %>% 
    mutate(totalCDR3ForCellAndLocus = n())
  
  print(paste0('total merged rows after group: ', nrow(results)))

  dateFormat <- date_format(format ="%Y-%m-%d")
  results$dateFormatted <- dateFormat(as.Date(results$SampleDate))
  
  results <- groupSingleCellData(results, lowFreqThreshold = lowFreqThreshold, backgroundThreshold = backgroundThreshold)
  if (!all(is.na(results)) && nrow(results) > 0) {
    results <- filterBasedOnCellCount(results)
  } else {
    print('No rows in results')
    return(NA)
  }
  
  # Add clone names:
  labelDf <- suppressWarnings(labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath=folderPath, 
    schemaName="tcrdb", 
    queryName="clones", 
    viewName="", 
    showHidden=TRUE,
    colSelect=c('cloneName','chain','cdr3','displayname'),
    colNameOpt='rname'
  ))
  colnames(labelDf)[colnames(labelDf)=="chain"] <- "locus"
  
  labelDf$LabelCol <- coalesce(labelDf$displayname, labelDf$clonename)
  
  labelsGroup <- labelDf %>% 
    group_by(locus, cdr3) %>% 
    summarize(CloneName = toString(sort(unique(LabelCol))))
  
  results <- merge(results, labelsGroup, by.x = c('locus', 'cdr3'), by.y = c('locus', 'cdr3'), all.x = TRUE, all.y = FALSE)
  
  results$LocusCDR3 <- as.character(results$LocusCDR3)
  results$LocusCDR3[!is.na(results$CloneName) & !results$IsFiltered] <- results$CloneName[!is.na(results$CloneName) & !results$IsFiltered]
  results$LocusCDR3 <- as.factor(results$LocusCDR3)
  
  return(list(results=results, meta = meta))
}

qualityFilterResults <- function(df = NULL, minLibrarySize=2000000, minReadCount=150000, minReadCountForSingle=150000, minReadCountForEnriched=5000, savePrefix = NULL){
  print(paste0('performing quality filter, initial rows: ', nrow(df)))
  
  dropCols <- c('SeqDataType', 'SubjectId', 'SampleDate', 'Peptide', 'Population', 'TotalReads', 'IsSingleCell', 'count')
  toDrop <- unique(df[(df$SeqDataType %in% c('Full_Transcriptome', '10x_VDJ') & !df$IsSingleCell & df$TotalReads < minReadCount),dropCols])
  
  r <- nrow(df)
  df <- df[!(df$SeqDataType %in% c('Full_Transcriptome', '10x_VDJ') & !df$IsSingleCell & df$TotalReads < minReadCount),]
  print(paste0('total bulk RNA rows dropped for read count: ', (r - nrow(df)), ', remaining: ', nrow(df)))

  r <- nrow(df)
  toDrop <- rbind(toDrop, unique(df[(df$SeqDataType %in% c('Full_Transcriptome', '10x_VDJ') & df$IsSingleCell & df$TotalReads < minReadCountForSingle), dropCols]))
  df <- df[!(df$SeqDataType %in% c('Full_Transcriptome', '10x_VDJ') & df$IsSingleCell & df$TotalReads < minReadCountForSingle),]
  print(paste0('total single cell RNA rows dropped for read count: ', (r - nrow(df)), ', remaining: ', nrow(df)))

  r <- nrow(df)
  toDrop <- rbind(toDrop, unique(df[(df$SeqDataType == 'TCR_Enriched' & df$TotalReads < minReadCountForEnriched), dropCols]))
  df <- df[!(df$SeqDataType == 'TCR_Enriched' & df$TotalReads < minReadCountForEnriched),]
  print(paste0('total TCR enriched rows dropped for read count: ', (r - nrow(df)), ', remaining: ', nrow(df)))
  
  write.table(toDrop, file = paste0(savePrefix, '.DropForReadCount.txt'), sep = '\t', row.names = FALSE, quote = FALSE)

  return(df)
}

groupSingleCellData <- function(df = NULL, lowFreqThreshold = 0.05, minTcrReads = 10, backgroundThreshold = 0.5, expandThreshold = 12){
  df <- df %>% 
    group_by(SeqDataType, SampleId, Population, Replicate, IsSingleCell) %>% 
    mutate(totalTcrReadsForGroup = sum(count))
  
  df <- df %>% 
    group_by(SeqDataType, SampleId, Population, Replicate, IsSingleCell, locus) %>% 
    mutate(totalTcrReadsForGroupAndLocus = sum(count))

  #bulk
  dfb <- df[!df$IsSingleCell,]
  columns <- c('DatasetId', 'SeqDataType', 'SampleId', 'Population', 'IsSingleCell', 'Label', 'SubjectId', 'SampleDate', 'Peptide', 'AssayType', 'locus', 'vhit', 'jhit', 'cdr3', 'count', 'dateFormatted', 'LocusCDR3', 'LocusCDR3WithUsage', 'Replicate', 'totalCellsForGroup', 'totalCellsForGroupAndLocus', 'Cells', 'percentage', 'percentageForLocus', 'totalCellsForCDR3')
  if (nrow(dfb) > 1) {
    dfb$totalCellsForGroupAndLocus <- dfb$totalCellsForGroupAndLocusBulk
  
    #this allows replicates to collapse:
    dfb <- dfb %>%
      group_by(DatasetId, SeqDataType, SampleId, Population, IsSingleCell, Label, SubjectId, SampleDate, Peptide, AssayType, locus, vhit, jhit, cdr3, dateFormatted, LocusCDR3, LocusCDR3WithUsage, Replicate, totalCellsForGroup, totalCellsForGroupAndLocus, totalTcrReadsForGroup, totalTcrReadsForGroupAndLocus, Cells) %>%
      summarize(count = sum(count))
    dfb$percentage = dfb$count / dfb$totalTcrReadsForGroup
    dfb$percentageForLocus = dfb$count / dfb$totalTcrReadsForGroupAndLocus
    r <- nrow(dfb)
    dfb <- dfb[dfb$totalTcrReadsForGroup >= minTcrReads,]
    print(paste0('total bulk rows filtered due to low TCR reads: ', (r - nrow(dfb))))
    
    dfb$totalCellsForCDR3 <- c(0)
    dfb <- as.data.frame(dfb[,columns])
    print(paste0('total bulk rows: ', nrow(dfb)))
  } else {
    print('No bulk rows')  
    dfb$totalCellsForCDR3 <- numeric()
    dfb$Cells <- integer()
    dfb$count <- integer()
    dfb$percentage <- numeric()
    dfb$percentageForLocus <- numeric()
    dfb$totalCellsForGroupAndLocus <- integer()
    
  }
  
  #single cell
  dfs <- df[df$IsSingleCell,]
  print(paste0('starting single cell rows: ', nrow(dfs)))
  if (nrow(dfs) > 0){
    dfs$totalCellsForGroupAndLocus <- dfs$totalCellsForGroupAndLocusSS
    dfs$ScaledCells <- dfs$Cells / dfs$totalCDR3ForCellAndLocus
    dfs$DatasetId <- paste0(dfs$SeqDataType, '-', dfs$SampleId)
    dfs <- dfs %>%
      group_by(DatasetId, SeqDataType, SampleId, Population, IsSingleCell, Label, SubjectId, SampleDate, Peptide, AssayType, locus, vhit, jhit, cdr3, dateFormatted, LocusCDR3, LocusCDR3WithUsage, Replicate, totalCellsForGroup, totalCellsForGroupAndLocus) %>%
      summarize(totalCellsForCDR3 = sum(ScaledCells))
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
    dfs$totalCellsForGroupAndLocus <- integer()
  }
  
  print(paste0('total single cell rows after grouping: ', nrow(dfs)))

  if (nrow(dfb) > 0 && nrow(dfs) > 0) {
    dfs <- as.data.frame(dfs[,columns])
    df <- rbind(dfb[columns], dfs[columns])  
  } else if (nrow(dfb) > 0) {
    df <- dfb[columns]
  } else {
    dfs <- as.data.frame(dfs[,columns])
  }
  
  print(paste0('combined: ', nrow(df)))
  df <- df %>% group_by(LocusCDR3, SubjectId) %>% mutate(maxPercentageInAnimal = max(percentage)) 
  df <- df %>% group_by(LocusCDR3, Peptide) %>% mutate(maxPercentageInPeptide = max(percentage)) 

  #track the highest frequency in a matched negative sample:
  negativePopulations <- list(
    c('TNF-Neg', 'TNF-Pos'), 
    c('Tetramer-Neg', 'Tetramer-Pos'),
    c('Tet-', 'Tet+')
  )
  
  negativePopulationNames <- character()
  for (pops in negativePopulations) {
    negativePopulationNames <- append(negativePopulationNames, pops[[1]])
  }
  negativePopulationNames <- negativePopulationNames[negativePopulationNames %in% as.character(df$Population)]
  
  df$IsFiltered <- c(FALSE)
  df$ratioToBackground <- c(0)
  if (length(negativePopulationNames) > 0) {
    negdf <- df %>% filter(Population %in% negativePopulationNames) %>% group_by(LocusCDR3, SampleId, Population) %>% summarize(maxPercentageInMatchedNeg = max(percentage)) 
    for (pops in negativePopulations) {
      p1 <- pops[[1]]
      p2 <- pops[[2]]
      negdf$Population <- as.character(negdf$Population)
      negdf$Population[negdf$Population == p1] <- c(p2)
      negdf$Population <- as.factor(negdf$Population)
    }
  
    df <- merge(df, negdf[c('LocusCDR3', 'SampleId', 'Population', 'maxPercentageInMatchedNeg')], by = c('LocusCDR3', 'SampleId', 'Population'), all.x = TRUE, suffix = c("", ""))
    rm(negdf)
    
    if (nrow(df) == 0) {
      return(df)
    }
  } else {
    print('no matching negatives found, skipping')
    df$maxPercentageInMatchedNeg <- 0
  }
  
  if (length(df$ratioToBackground[!is.na(df$maxPercentageInMatchedNeg)]) > 0){
    df$ratioToBackground[!is.na(df$maxPercentageInMatchedNeg)] <- (df$maxPercentageInMatchedNeg[!is.na(df$maxPercentageInMatchedNeg)] / df$percentage[!is.na(df$maxPercentageInMatchedNeg)])  
  } else {
    print(paste0('no rows have a matched negative sample'))
  }

  df$LocusCDR3Orig <- df$LocusCDR3
  if (lowFreqThreshold){
    print(paste0('Low freq threshold: ', lowFreqThreshold))
    
    df$LocusCDR3 <- as.character(df$LocusCDR3)
    df$LocusCDR3[df$maxPercentageInAnimal < lowFreqThreshold & is.na(df$CloneName)] <- c('Low Freq.')
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    df$IsFiltered[df$maxPercentageInAnimal < lowFreqThreshold] <- c(TRUE)
    
    df$LocusCDR3WithUsage <- as.character(df$LocusCDR3WithUsage)
    df$LocusCDR3WithUsage[df$maxPercentageInAnimal < lowFreqThreshold & is.na(df$CloneName)] <- c('Low Freq.')
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  }
  
  if (!is.na(backgroundThreshold)){
    total <- length(df$LocusCDR3[df$ratioToBackground > backgroundThreshold])
    if (total > 0){
      print(paste0('total filtered due to high background: ', total))
    }

    df$LocusCDR3 <- as.character(df$LocusCDR3)
    df$LocusCDR3[df$ratioToBackground > backgroundThreshold] <- c('Background')
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    df$IsFiltered[df$ratioToBackground > backgroundThreshold] <- c(TRUE)

    df$LocusCDR3WithUsage <- as.character(df$LocusCDR3WithUsage)
    df$LocusCDR3WithUsage[df$ratioToBackground > backgroundThreshold] <- c('Background')
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  }
  
  if (expandThreshold) {
    df$ExpandThreshold <- df$percentageForLocus / df$percentage
    
    total <- length(df$LocusCDR3[df$ExpandThreshold > expandThreshold])
    if (total > 0){
      print(paste0('total filtered due to low reads in locus: ', total))
      print(paste0('max expand: ', max(df$ExpandThreshold)))
    }
    
    df$LocusCDR3 <- as.character(df$LocusCDR3)
    df$LocusCDR3[df$ExpandThreshold > expandThreshold] <- c('Low Reads for Locus')
    df$LocusCDR3 <- as.factor(df$LocusCDR3)
    df$IsFiltered[df$ExpandThreshold > expandThreshold] <- c(TRUE)
    
    df$LocusCDR3WithUsage <- as.character(df$LocusCDR3WithUsage)
    df$LocusCDR3WithUsage[df$ExpandThreshold > expandThreshold] <- c('Low Reads for Locus')
    df$LocusCDR3WithUsage <- as.factor(df$LocusCDR3WithUsage)
  }
  
  df <- df %>% group_by(SampleId, Population, Replicate, IsSingleCell) %>% mutate(totalTcrReadsForGroupAfterLowFreq = sum(count)) 
  
  #TODO
  df$percentageAfterLowFreq <- df$count / df$totalTcrReadsForGroupAfterLowFreq
  
  df <- arrange(df, LocusCDR3)
  
  return(df)
}

filterBasedOnCellCount <- function(df = NULL){
  startRows <- nrow(df)
  df$minPctToInclude <- c(0)
  if (sum(df$IsSingleCell) > 0) {
    df$minPctToInclude[df$IsSingleCell] <- c(0.1)
    df$minPctToInclude[df$IsSingleCell] <- (1 / df$Cells[df$IsSingleCell]) 
  }

  print(paste0('rows filtered by cell count: ', (startRows - nrow(df))))
  
  return(df)
}

makePlot <- function(locus, df, pal = 'Reds', bulkOnly = FALSE, tnfPosOnly = FALSE, showLegend = FALSE, yMax=1.05, yField='percentage', xField='Label', fillField='LocusCDR3', facet1 = 'Population', facet2 = 'dateFormatted', doFacet=TRUE, lowFreqThreshold = 0.05, reverseColors = F, colorSteps = NA, locusDataFile = NA, repeatGrouping = F){
  locusData <- df[df$locus == locus,]
  if (lowFreqThreshold){
    print(paste0('Plot low freq threshold: ', lowFreqThreshold))
    locusData[[fillField]] <- as.character(locusData[[fillField]])
    print(paste0('dropped: ', length(locusData[[fillField]][locusData[[yField]] < lowFreqThreshold])))
    locusData[[fillField]][locusData[[yField]] < lowFreqThreshold] <- c('Low Freq.')
    locusData[[fillField]] <- naturalfactor(locusData[[fillField]])
  }
  
  if (bulkOnly){
    locusData <- locusData[!locusData$IsSingleCell,]
  }
  
  if (tnfPosOnly){
    locusData <- locusData[locusData$Population == 'TNF-Pos',]
  }
  
  locusData[fillField] <- naturalfactor(as.character(locusData[[fillField]]))
  totalCDR3 = length(unique(locusData[[fillField]]))
  #print(paste(locus, ', total CDR3:', totalCDR3))
  
  if (nrow(locusData) == 0){
    return(ggplot(data.frame()) +geom_blank())
  }
  
  if (is.na(colorSteps)) {
    colorSteps <- min(9, totalCDR3)
    colorSteps <- max(3, colorSteps)
  }
  
  #print(paste0('total color steps: ', colorSteps))
  getPalette = colorRampPalette(brewer.pal(colorSteps, pal))
  #colorsDT <-  data.table(fillField=levels(locusData[[fillField]]), Color=getPalette(totalCDR3))
  #locusData <- merge(locusData, colorsDT, by.x = c(fillField), by.y=c('fillField'), all = TRUE)
  #locusData$Color[locusData[[fillField]] == 'Low Freq.'] <- c('#D3D3D3')
  #locusData$Color[locusData[[fillField]] == 'Background'] <- c('#000000')
  #locusData$Color <- as.factor(locusData$Color)
  
  extraColors <- character()
  if ('Background' %in% locusData[[fillField]]) {
    extraColors <- c(extraColors, '#000000')  #black
  }
  

  if ('Low Freq.' %in% locusData[[fillField]]) {
    extraColors <- c(extraColors, '#D3D3D3')  #grey  
  }
  
  if ('Low Reads for Locus' %in% locusData[[fillField]]) {
    extraColors <- c(extraColors, '#808080')  #grey  
  }
  
  colorValues <- getPalette(totalCDR3 - length(extraColors))
  if (reverseColors) {
    colorValues <- rev(colorValues)
  }
  colorValues <- c(extraColors, colorValues)
  #write.table(file = 'colorValues.csv', colorValues)

  if (doFacet){
    N <- data.table(locusData)
    N <- N[, c(facet2, xField), with = FALSE]
    names(N) <- c('GroupField', 'CountField')
    
    N <- as.data.frame(N) %>% 
      group_by(GroupField) %>% 
      summarise(V1 = n_distinct(CountField))
    
    names(N) <- c(facet2, 'V1')
    
    N$scaleFactor = N$V1 / max(N$V1)
    locusData = merge(locusData, N, by.x = c(facet2), by.y = c(facet2), all.x = TRUE)
    widthStr <- '0.5' #*scaleFactor'
    
    #locusData[[facet2]] <- naturalfactor(locusData[[facet2]], ordered = T)
    #locusData[[xField]] <- naturalfactor(locusData[[xField]], ordered = T)
  } else {
    locusData$scaleFactor <- c(1)
    widthStr <- '0.5'
  }
  
  if (showLegend){
    lp = 'right'
  }
  else {
    lp = 'none'
  }
  
  if (!is.na(locusDataFile)) {
    write.table(locusData, file = locusDataFile, sep= '\t', quote = F, row.names = F)
  }
  
  if (repeatGrouping) {
    print('regrouping data')
    groupFields <- c(fillField, xField, 'LocusCDR3Orig')
    if (doFacet && !is.null(facet1) && !is.na(facet1) && facet1 != '.') {
      groupFields <- c(groupFields, facet1)
    }
    
    if (doFacet && !is.null(facet2) && !is.na(facet2) && facet2 != '.') {
      groupFields <- c(groupFields, facet2)
    }
    
    locusData <- locusData %>% group_by_at(groupFields) %>% summarise_at(c(yField, 'count', 'percentage'), sum)
    print(str(locusData))
    locusData$hoverText <- paste(
      fillField, ': ', locusData[[fillField]],
      '<br>', xField, ': ', locusData[[xField]],
      '<br>', yField, ': ', locusData[[yField]],
      '<br>Cell/Read #: ', locusData$count
    )
    
  } else {
    locusData$hoverText <- paste(
      fillField, ': ', locusData[[fillField]],
      '<br>', xField, ': ', locusData[[xField]],
      '<br>', yField, ': ', locusData[[yField]],
      '<br>Cell/Read #: ', locusData$count,
      '<br>Total Cell/Read #: ', locusData$totalTcrReadsForGroupAfterLowFreq,
      '<br>Cells For Group: ', locusData$totalCellsForGroup, 
      '<br>Cells For Group And Locus: ', locusData$totalCellsForGroupAndLocus, 
      '<br>Cells For CDR3: ', locusData$totalCellsForCDR3,
      '<br>cDNA: ', locusData$DatasetId
      )
  }
  
  locusData[[fillField]] <- forcats::fct_relevel(locusData[[fillField]], 'Low Freq.')
  #locusData[[xField]] <- as.character(locusData[[xField]])
  P1<-ggplot(locusData) +
    aes_string(x = xField, y = yField, order = 'percentage', width = widthStr, text = 'hoverText') +
    geom_bar(aes_string(fill = fillField), stat = 'identity', colour="black", position = 'fill') +
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
    facetFormula <- as.formula(paste0(facet1, ' ~ ', facet2))
    P1 <- P1 + facet_grid(facetFormula, scales = 'free_x', space = 'free_x')
  }
  
  return(P1)
}

makePlotGroup <- function(df, subDir, basename, tnfPosOnly=FALSE, locusPlots = c('TRB'), pal='Reds', showLegend=FALSE, bulkOnly=FALSE, yMax=NA, yField='percentage', xField='Label', fillField='LocusCDR3', replicateAsSuffix = TRUE, allLocusWidth=1800, singleLocusWidth=1000, singleLocusHeight=500, asEPS = FALSE, doFacet=TRUE, facet1 = 'Population', facet2 = 'dateFormatted', lowFreqThreshold = 0.05, plotTitle = NULL, reverseColors = F, doPrint = T, repeatGrouping = F){
  plots = list()
  i = 1

  for (locus in sort(unique(df$locus))){
    P1 <- makePlot(locus, df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly, yField=yField, xField=xField, fillField=fillField, facet1=facet1, facet2=facet2, doFacet=doFacet, lowFreqThreshold = lowFreqThreshold, reverseColors = reverseColors, repeatGrouping = repeatGrouping)
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

  returnPlots <- list()
  for (locus in locusPlots){
    print(paste0('Making plot for locus: ', locus))
    if (asEPS) {
      postscript(paste0(fn,'.',locus,'.eps'), height = singleLocusHeight, width = singleLocusWidth)
    } else {
      png(paste0(fn,'.', locus, '.png'), height = singleLocusHeight, width = singleLocusWidth)
    }

    locusDataFile <- paste0('./',subDir,'/', basename, '.', locus, '.txt')
    print(locusDataFile)
    P_B <- makePlot(locus, df, tnfPosOnly=tnfPosOnly, pal=pal, showLegend=showLegend, yMax=yMax, bulkOnly=bulkOnly,yField='percentageForLocus', xField=xField, fillField=fillField, facet1=facet1, facet2=facet2, doFacet=doFacet, lowFreqThreshold=lowFreqThreshold, reverseColors = reverseColors, locusDataFile = locusDataFile)
    if (is.null(plotTitle)) {
      P_B = P_B + labs(title = paste0(basename, " ", locus))
    } else if (is.na(plotTitle)) {
      P_B = P_B + labs(title = NULL)
    } else {
      P_B = P_B + labs(title = plotTitle)
    }
    
    print(P_B)
    dev.off()
    
    if (doPrint) {
      print(ggplotly(P_B, tooltip = "hoverText"))
    }
    
    returnPlots[locus] <- list(P_B)
  }

  return(returnPlots)
}

getIDsForGroup <- function(projectName) {
  groupMembers <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Labs/Bimber", 
    schemaName="laboratory", 
    queryName="project_usage", 
    colSelect = "subjectid",
    colFilter=makeFilter(c("project", "EQUALS", projectName)),
    containerFilter=NULL, 
    colNameOpt="rname"
  )
  
  groupIDs <- sort(unique(groupMembers$subjectid))
  print(groupIDs)
  
  return(groupIDs)
}

getWorkbooksWithData <- function(groupIDs) {
  workbooks <- labkey.selectRows(
    baseUrl="https://prime-seq.ohsu.edu", 
    folderPath="/Labs/Bimber", 
    schemaName="singlecell", 
    queryName="cdna_libraries", 
    colSelect = "workbook/workbookid",
    colFilter=makeFilter(c("sortId/sampleId/subjectId", "CONTAINS_ONE_OF", paste0(groupIDs, collapse = ';'))),
    containerFilter=NULL, 
    colNameOpt="rname"
  )

  workbookIds <- sort(unique(workbooks$workbook_workbookid))
  print(workbookIds)
  
  return(workbookIds)
}