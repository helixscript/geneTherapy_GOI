library(dplyr)
library(readr)
library(GenomicRanges)
library(parallel)
library(data.table)
library(ggplot2)
library(ggrepel)
library(grDevices)
library(RColorBrewer)
library(kableExtra)
library(geneRxCluster)
library(UpSetR)

set.seed(1)
source('lib.R')

# Read in the configuration file.
config <- yaml::read_yaml('config.yml')

if(! dir.exists(config$outputDirPath)) dir.create(config$outputDirPath)
if(! dir.exists(config$outputDirPath)) stop('Error - can not create output directory.')

# Read intSite data and limit to sites near genes.
d <- read_tsv(config$inputDataPath, progress = FALSE, col_types = cols())

# Convert patient report format to expected format.
if('chromosome' %in% names(d))  d<- rename(d, 'chromosome' = 'seqnames')
if('internalSampleID' %in% names(d)) d <- rename(d, 'internalSampleID' = 'GTSP')
if('subject' %in% names(d)) d <- rename(d, 'subject' = 'patient')
if('position' %in% names(d)){
  d$start <- d$position
  d$end <- d$position
  d$position <- NULL
}
if(! 'posid' %in% names(d)) d$posid <- paste0(d$seqnames, d$strand, d$start)

if(! 'timePointDays' %in% names(d)){
  tps <- expandTimePoints(d$timePoint)
  d <- bind_cols(d, tps)
  if(any(is.na(d$timePointDays))){
    message('Warning - some time points could not be convereted to days.')
    d <- d[which(is.na(d$timePointDays)) * -1,]
  }
}



# Limit sites to study those within config$maxDistNearestGene nt of a gene boundary
d <- subset(d, abs(nearestFeatureDist) <= config$maxDistNearestGene)

# Since we are studying genes, expand instances where a site intersects more than one gene.
#    chr7+543431  ABC1,XYZ1  ->
#    chr7+543431  ABC1
#    chr7+543431  XYZ1

d$nearestFeature <- toupper(gsub('\\s', '', d$nearestFeature))
d <- mutate(d, nearestFeature = strsplit(nearestFeature, ",")) %>% tidyr::unnest(nearestFeature)


# Limit sites to those from samples with at least config$minSampleAbund estimated cells.
sampleTotalAbundance <- group_by(d, GTSP) %>%
                        summarise(totalAbund = sum(estAbund)) %>%
                        ungroup() %>%
                        filter(totalAbund >= config$minSampleAbund)

d <- subset(d, GTSP %in% sampleTotalAbundance$GTSP)


# Limit sites to those with nearest genes seen in at least config$minGeneSubjects patients.
genePatients <- group_by(d, nearestFeature) %>%
                summarise(nPatients = n_distinct(patient)) %>%
                ungroup() %>%
                filter(nPatients >= config$minGeneSubjects)

d <- subset(d, nearestFeature %in% genePatients$nearestFeature)


# Create a gene-centric table of maximal abundance and relative abundance values.
# This table will be joined to the main data table (d).

geneAbundances <- group_by(subset(d, d$timePointDays > config$earlyVsLateCutoffDays), nearestFeature) %>%
                  summarise(maxAbund = max(estAbund),
                            maxRelAbund = max(relAbund)) %>%
                  ungroup()

# Load external data files.
HGNC     <- loadHGNCfull()
HGNCprev <- loadHGNCprev()

# Use the retired HGNC gene ids to update gene symbols to current symbols.
d$nearestFeature <- updateGeneSymbols(d$nearestFeature)


# Separate the data into early and late data sets.
a <- d[d$timePointDays <= config$earlyVsLateCutoffDays,]
b <- d[d$timePointDays >  config$earlyVsLateCutoffDays,]


# ScanStats
#-------------------------------------------------------------------------------

# Create GRanges for early and late sites.
ag <- makeGRangesFromDataFrame(a, keep.extra.columns = TRUE)
bg <- makeGRangesFromDataFrame(b, keep.extra.columns = TRUE)

# Remove duplicate sites from each group across all subjects.
ag <- ag[! duplicated(a$posid)]
bg <- bg[! duplicated(b$posid)]

# Interpret value1 and value2 from the clusterSource flag.
scanStatsDF <- data.frame(scanStats(ag, bg))
scanStatsDF <- bind_rows(lapply(split(scanStatsDF, 1:nrow(scanStatsDF)), function(x){
                 if(x$clusterSource == 'A'){
                   x$earlySites = max(c(x$value1, x$value2))
                   x$lateSites = min(c(x$value1, x$value2))
                 } else {
                   x$earlySites = min(c(x$value1, x$value2))
                   x$lateSites = max(c(x$value1, x$value2))
                 }
  
                 x
                })) %>% select(seqnames, start, end, width, earlySites, lateSites, target.min)

# Clean up scanStats result data.frame.
names(scanStatsDF) <- c('chromosome', 'start', 'end', 'width', 'earlySites', 'lateSites', 'target.min')
scanStatsDF$target.min <- round(scanStatsDF$target.min, 3)

# Create a data.frame corresponding to the data used for the scanStats analysis to draw nearest genes from.
g <- data.frame(c(ag, bg))

# Look within each cluster and tally the number of nearest genes for final report.
scanStatsDF <- bind_rows(lapply(split(scanStatsDF, 1:nrow(scanStatsDF)), function(x){
  tab <- sort(table(subset(g, seqnames == x$chromosome & start >= x$start & start <= x$end)$nearestFeature), decreasing = TRUE)
  x$nearestGenes <- paste0(paste0(names(tab), ' (', tab, ')'), collapse = ', ')
  x <- mutate(x, cluster = paste0(x$chromosome, ':', x$start, '-', x$end), .before = 'chromosome')
  x$url <- paste0(config$scanStatsBEDfileURL, "&position=", x$chromosome, "%3A", x$start-10, "-", x$end+10)
  select(x, -chromosome, start, end)
}))

# Create UCSF BED files for early and late sites.
ag <- group_by(a, patient, posid) %>% dplyr::slice_max(estAbund, with_ties = FALSE) %>% ungroup()
ag$siteLabel <- gsub('\\s', '_', paste0(ag$patient, ':', ag$start))
createIntUCSCTrack(ag, title = paste0('sites <= ', config$earlyVsLateCutoffDays, ' days'), outputFile = file.path(config$outputDirPath, 'early.ucscTrack'), siteLabel = 'siteLabel')

bg <- group_by(b, patient, posid) %>% dplyr::slice_max(estAbund, with_ties = FALSE) %>% ungroup()
bg$siteLabel <- gsub('\\s', '_', paste0(bg$patient, ':', bg$start))
createIntUCSCTrack(bg, title = paste0('sites > ', config$earlyVsLateCutoffDays, ' days'), outputFile = file.path(config$outputDirPath, 'late.ucscTrack'), siteLabel = 'siteLabel')

# Merge tracks together into a single track file.
# (!) This file needs to be moved to a web server and accessible via the URL
#     defined in config$scanStatsBEDfileURL.

system(paste0('cat ', file.path(config$outputDirPath, 'early.ucscTrack'), ' ', file.path(config$outputDirPath, 'late.ucscTrack'), ' > ', file.path(config$outputDirPath, 'tracks.ucsc')))


# End ScanStats
#-------------------------------------------------------------------------------


# Count unique sites in each time point grouping.
a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)

genesToSubjects <- group_by(d, nearestFeature) %>% 
                   summarise(subjects = n_distinct(patient),
                             totalSites = n_distinct(posid)) %>% 
                   ungroup()


# Build a table describing genes in later timepoints.
lateGeneData <- subset(b, timePointDays >= config$longitudinal_minTimeDays) %>%
                group_by(nearestFeature) %>%
                  summarise(latestTimePointDays = sort(unique(timePointDays), decreasing = TRUE)[1],
                  longitudinalSubjects = n_distinct(patient),
                  longitudinalSites = n_distinct(posid),
                  longitudinalTimePoints = n_distinct(timePointDays)) %>%
                ungroup()

# For early and late data groups (a and b), create tables of site counts.
ag <- group_by(a, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()
bg <- group_by(b, nearestFeature) %>% summarise(nSites = n_distinct(posid)) %>% ungroup()            

              
# Update table names to reflect their source since they will be jointed together
# otherwise they identical column names would be appended with .x and .y suffixes. 
names(ag) <- paste0(names(ag), '_a')
names(bg) <- paste0(names(bg), '_b')              


# Count TOTAL sites in early and late groupings since they will be constants used
# in the following operation.
a_totalSites <- n_distinct(a$posid)
b_totalSites <- n_distinct(b$posid)


# Join the early and late tables together by gene symbols.
# Then, for each gene, determine change in integration frequency and the significance
# of that change using Fisher's Exact tests. After each gene test is completed,
# join pre-built tables to the result.

k <- left_join(ag, bg, by = c('nearestFeature_a' = 'nearestFeature_b')) %>% 
     tidyr::drop_na() %>%
     group_by(nearestFeature_a) %>%
       mutate(gene = nearestFeature_a,
              pVal = fisher.test(matrix(c(a_totalSites - nSites_a, 
                                          b_totalSites - nSites_b, 
                                          nSites_a , 
                                          nSites_b), 
                                        ncol = 2, byrow = FALSE))$p.value,
              earlyCount = nSites_a,
              lateCount  = nSites_b,
              earlyFreq  = nSites_a / a_totalSites,
              lateFreq   = nSites_b / b_totalSites,
              percentChange = ((lateFreq - earlyFreq) / earlyFreq) * 100) %>%
     ungroup() %>% 
     select(-nearestFeature_a) %>%
     left_join(genesToSubjects, by = c('gene' = 'nearestFeature')) %>%
     left_join(geneAbundances, by = c('gene' = 'nearestFeature')) %>%
     left_join(lateGeneData, by = c('gene' = 'nearestFeature')) %>% 
     select(-nSites_a, -nSites_b) %>%
     as.data.table()
 

# We will study the result both with and without correction for multiple comparisons. 
# Here we save the uncorrected pValues (raw), correct for multiple comparisons (adj)
# Then set pVal to NA. pVal will be toggled latter between raw and adj values.
k$pVal.raw <- k$pVal
k$pVal.adj <- p.adjust(k$pVal, method = 'BH')
k$pVal     <- NA


# Read in oncogene lists.
oncoGeneLists <- list('cosmic'     = unique(updateGeneSymbols(unique(readLines(config$COSMIC_oncogene_table)))),
                      'cosmic_tsg' = unique(updateGeneSymbols(unique(readLines(config$COSMIC_tsg_table)))),
                      'allOnco'    = unique(updateGeneSymbols(unique(readLines(config$allOnco_oncogene_table)))))

# Test if gene symbols are in the included oncogene lists.
k$oncoGene <- ifelse(k$gene %in% unique(oncoGeneLists[['cosmic']]), 
                     ifelse(k$gene %in% oncoGeneLists$cosmic_tsg, 'TSG', 'Yes'), '')


# First, set the working pValues to the BH corrected pValues.
k$pVal <- k$pVal.adj      # Set working pValue
k <- setCategories(k)     # Set categories based on working pValue
pValAdjVolcanoPlot <- volcanoPlot(k)      # Build volvano plot.
ggsave(file.path(config$outputDirPath, 'pValAdjVolcanotPlot.pdf'), pValAdjVolcanoPlot, units = 'in', width = 10, height = 15)
readr::write_tsv(k, file.path(config$outputDirPath, 'geneData_corrected_pVals.tsv'))


# Run a series of enrichment tests on 3rd party gene lists.
correctedGeneListTests <- geneListTests() 
correctedGeneListTests$sigMarks <- ifelse(correctedGeneListTests$pval <= 0.001, '***',
                                            ifelse(correctedGeneListTests$pval <= 0.01, '**',
                                                   ifelse(correctedGeneListTests$pval <= 0.05, '*', 'n.s.')))
correctedGeneListTests$Enriched <-  paste (ifelse(correctedGeneListTests$Enriched, 'Enriched', 'Depleted'), correctedGeneListTests$sigMarks)
correctedGeneListTests <- select(arrange(correctedGeneListTests, OncoGeneList), -pval, -sigMarks) 


# Second, set the working pValues to the uncorrected pValues.
k$pVal <- k$pVal.raw
k <- setCategories(k)
pValVolcanoPlot <- volcanoPlot(k)
ggsave(file.path(config$outputDirPath, 'pValRawVolcanotPlot.pdf'), pValVolcanoPlot, units = 'in', width = 10, height = 12)
readr::write_tsv(k, file.path(config$outputDirPath, 'geneData_uncorrected_pVals.tsv'))

# Run a series of enrichment tests on 3rd party gene lists for uncorrected values.
uncorrectedGeneListTests <- geneListTests()
uncorrectedGeneListTests$sigMarks <- ifelse(uncorrectedGeneListTests$pval <= 0.001, '***',
                                              ifelse(uncorrectedGeneListTests$pval <= 0.01, '**',
                                                      ifelse(uncorrectedGeneListTests$pval <= 0.05, '*', 'n.s.')))
uncorrectedGeneListTests$Enriched <-  paste (ifelse(uncorrectedGeneListTests$Enriched, 'Enriched', 'Depleted'), uncorrectedGeneListTests$sigMarks)
uncorrectedGeneListTests <- select(arrange(uncorrectedGeneListTests, OncoGeneList), -pval, -sigMarks) 

upSetPlot <- buildUpSetPlot()

rmarkdown::render('report.Rmd',
                  output_file = file.path(config$outputDirPath, 'report.pdf'),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = config$reportTitle,
                                'author' = config$reportAuthor))

invisible(unlink('report_cache', recursive = TRUE))
invisible(unlink('report_files', recursive = TRUE))
invisible(unlink(file.path(config$outputDirPath, 'report_files'), recursive = TRUE))