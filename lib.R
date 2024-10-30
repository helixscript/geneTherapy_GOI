
ppNum <- function (n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

loadHGNCprev <- function(){
  HGNC <- readr::read_tsv(config$HGNC_symbols, col_types = cols())
  
  HGNC <- select(HGNC, `Approved symbol`, `Previous symbols`)
  names(HGNC) <- c('symbol', 'prev')
  
  HGNC$symbol <- toupper(HGNC$symbol)
  
  HGNC$prev <- toupper(gsub('\\s', '', HGNC$prev))
  
  HGNC <- distinct(HGNC)
  
  # Create a table that only has alternative symbols.
  mutate(HGNC, prev = strsplit(prev, ",")) %>%
  tidyr::unnest(prev) %>%
  tidyr::drop_na() %>%
  distinct() %>%
  as.data.table()
}


loadHGNCfull <- function(){
  HGNC <- readr::read_tsv(config$HGNC_symbols, col_types = cols())
  
  HGNC <- select(HGNC, `Approved symbol`, `Previous symbols`, Status, `Alias symbols`)
  names(HGNC) <- c('symbol', 'prevSymbols', 'status', 'aliasSymbols')
  
  HGNC$symbol <- toupper(HGNC$symbol)
  HGNC$prevSymbols <- toupper(gsub('\\s', '', HGNC$prevSymbols))
  HGNC$aliasSymbols <- toupper(gsub('\\s', '', HGNC$aliasSymbols))

  as.data.table(distinct(HGNC))
}


buildUpSetPlot <- function(){
  z <- bind_rows(lapply(split(k, 1:nrow(k)), function(x){
         if(x$flag == 'none') return(tibble())
         o <- tibble('gene' = x$gene, 'Depleated' = 0, 'Enriched' = 0, 'Abundant' = 0, 'Longintudinal' = 0)
         if(grepl('D', x$flag)) o$Depleated <- 1
         if(grepl('E', x$flag)) o$Enriched <- 1
         if(grepl('A', x$flag)) o$Abundant <- 1
         if(grepl('L', x$flag)) o$Longintudinal <- 1
         o
  })) %>% as.data.frame()
  rownames(z) <- z$gene
  z$gene <- NULL
        
  upset(z, order.by="freq", text.scale = 1.2, point.size = 3)
}


geneListTests <- function(){
  bind_rows(lapply(c('D', 'E', 'A', 'L'), function(x){
    
    # Subset current gene list.
    a <- k[grepl(x, k$flag),]
    b <- k[! grepl(x, k$flag),]
    
    bind_rows(mapply(function(genes, label){
      tibble('Class' = x, 
                  'Genes in Class' = n_distinct(a$gene),
                  OncoGeneList = label,
                  'Percent genes in list' = sprintf("%.1f%%", (sum(a$gene %in% genes) / n_distinct(a$gene))*100),
                  Enriched = (sum(a$gene %in% genes) / n_distinct(a$gene)) > (sum(b$gene %in% genes) / n_distinct(b$gene)),
                  pval = fisher.test(matrix(c(sum(a$gene %in% genes),
                                              sum(! a$gene %in% genes),
                                              sum(b$gene %in% genes),
                                              sum(! b$gene %in% genes)), ncol = 2, byrow = FALSE))$p.val)
    }, oncoGeneLists, names(oncoGeneLists), SIMPLIFY = FALSE))
  }))
}

logValue <- function(x) log2(abs(x)) * ifelse(x < 0, -1, 1)


updateGeneSymbols <- function(g){
  g <- toupper(gsub('\\s', '', g))
  g2 <- g
  for(gene in unique(g)){
    if(gene %in% HGNCprev$prev & ! gene %in% HGNCprev$symbol){
      o <- HGNCprev[prev == gene]
      if(nrow(o) == 1){
        message(gene, '-> ', o$symbol)
        g2[g2 == gene] <- o$symbol
      }
    }
  }
  
  if('geneNameReassignment_table' %in% names(config)){
    if(file.exists(config$geneNameReassignment_table)){
      o <- readr::read_tsv(config$geneNameReassignment_table, col_names = TRUE)
      if(any(o$from %in% g2)){
        o <- o[o$from %in% g2,]
        
        invisible(lapply(split(o, 1:nrow(o)), function(oo){
          g2[which(g2 == oo$from)] <<- oo$to
        }))
        
      }
    }
  }
  
  message(sprintf("%.2f%%", (sum(g != g2) / length(g))*100), ' gene symbols updated.')
  g2
}

setCategories <- function(k){
  depletedGenes <- subset(k, percentChange <= 0 & pVal <= 0.05)
  enrichedGenes <- subset(k, percentChange > 0 & pVal <= 0.05)
  
  abundantGenes <- subset(k, maxAbund >= config$minAbundCatThreshold)

  k$flag <- ''
  k[k$percentChange <= 0 & 
    k$gene %in% depletedGenes$gene,]$flag <- paste0(k[k$percentChange <= 0 & 
                                                        k$gene %in% depletedGenes$gene,]$flag, 'D')
  
  k[k$percentChange > 0 & 
    k$gene %in% enrichedGenes$gene,]$flag <- paste0(k[k$percentChange > 0 & 
                                                      k$gene %in% enrichedGenes$gene,]$flag, 'E')
  
  k[k$gene %in% abundantGenes$gene,]$flag <- paste0(k[k$gene %in% abundantGenes$gene,]$flag, 'A')
  
  # Latter timepoints already baked into longitudinalSubjects and longitudinalSites columns.
  
  k[k$longitudinalSubjects >= config$longitudinal_minNumSubjects & 
    k$longitudinalSites >= config$longitudinal_minNumSites,]$flag <- paste0(k[k$longitudinalSubjects >= config$longitudinal_minNumSubjects & 
                                                                            k$longitudinalSites >= config$longitudinal_minNumSites,]$flag, 'L')
  k[k$flag == '',]$flag <- 'none'
  k
}

volcanoPlot <- function(k){
  k$volcanoPlot_x <- logValue(k$percentChange)
  k$volcanoPlot_y <- log(1 / k$pVal)
  
  k$geneLabel <- k$gene
  o <- split(k, k$percentChange >= 0)
  
  o[[1]] <- arrange(o[[1]], desc(volcanoPlot_y))
  o[[1]][(config$volcanoPlot_numTopGeneLabels+1):nrow(o[[1]]),]$geneLabel <- 'none'
  
  o[[2]] <- arrange(o[[2]], desc(volcanoPlot_y))
  o[[2]][(config$volcanoPlot_numTopGeneLabels+1):nrow(o[[2]]),]$geneLabel <- 'none'
  
  k <- bind_rows(o)
  
  flagLevels <- unique(k$flag)
  
  flagLevels <- flagLevels[flagLevels != 'none']
  flagLevels <- c('none', flagLevels[order(nchar(flagLevels))])
  
  k$flag <- factor(k$flag, levels = flagLevels)
  
  colors <- c('gray70', brewer.pal(12, 'Paired')) # grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(k$flag) - 1))
  
  jitter_pos <- position_jitter(width = 0.35, height = 0.35, seed = 1)
  
  k <- bind_rows(list(subset(k, flag == 'none'), subset(k, flag != 'none')))
  
  k$geneLabel <- ifelse(k$geneLabel == 'none', '', k$geneLabel)
  
  k$oncoGene <- ifelse(k$oncoGene == '', 'No', k$oncoGene)
  
  # Reorder rows so that interesting points will be displayed on the top of the data point stack.
  nk <- nrow(k)
  k <- bind_rows(subset(k, oncoGene == 'No' & oncoGene != 'TSG' & flag != 'EAL'), 
                 subset(k, oncoGene == 'Yes' & oncoGene != 'TSG' & flag != 'EAL'), 
                 subset(k, oncoGene != 'TSG' & flag == 'EAL'), 
                 subset(k, oncoGene == 'TSG'))
  
  if(nrow(k) != nk) stop('Error - reordering of volcano plot data points failed.')
  
  ggplot(subset(k, abs(percentChange) >= 2), aes(volcanoPlot_x, volcanoPlot_y, label = geneLabel, fill = flag, shape = oncoGene)) + 
       scale_fill_manual(name = 'Class', values = colors, labels = flagLevels, drop = FALSE) +
       scale_x_continuous(limits = c(min(k$volcanoPlot_x - 3), max(k$volcanoPlot_x))) +
       scale_y_continuous(limits = c(-1, max(k$volcanoPlot_y+3))) +
       scale_shape_manual(values = c(21, 24, 22)) +
       geom_hline(yintercept = log(1/0.05), color = 'black', linetype = 'dashed') +
       geom_jitter(position = jitter_pos, size = 2, alpha = 0.9) +  
       geom_text_repel(position = jitter_pos, size = 3.25, point.size = 2, direction = 'both', seed = 1,
                       max.overlaps = Inf,  max.time = 10, max.iter = 50000) +
       labs(x = 'log2(percent gene integration frequency change)', y = 'log(1/p-value)') +
       theme(legend.position = "bottom",
             legend.key=element_blank(),
             text = element_text(size = 12),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")) +
       guides(fill = guide_legend(override.aes = list(size=4, shape=21), title.position = "top"),
              shape = guide_legend(override.aes = list(size=4), title.position = "top"))
}



#' Run geneRxCluster scan statistics.
#' Bioinformatics. 2014 Jun 1;30(11):1493-500. 
#' 
#' @param  kvals: The interpretation of kvals = 2L:5L is that genomic windows will be drawn for every consecutive groups of integration sites with group sizes of 2, 3, 4, and 5 integration sites. 
#' @return Data frame of identified clusters with counts and statistics where target.min is the smallest nominal False Discoveries Expected for each cluster.  
#'
#' @export
scanStats <- function(gr1, gr2, gr1.label = 'A', gr2.label = 'B', kvals = '15L:30L', nperm = 100,
                      cutpt.tail.expr = 'critVal.target(k, n, target = 5, posdiff = x)',
                      cutpt.filter.expr = 'as.double(apply(x, 2, median, na.rm = TRUE))'
                      #cutpt.filter.expr = paste0('apply(x, 2, quantile, probs = ', p, ', na.rm = TRUE)')
){
  library(GenomicRanges)
  library(geneRxCluster)
  
  gr1$clusterSource <- TRUE
  gr2$clusterSource <- FALSE
  
  gr3 <- c(gr1, gr2)
  
  df <- GenomicRanges::as.data.frame(gr3)
  df <- df[order(df$seqnames, df$start),]
  df$seqnames <- droplevels(df$seqnames)
  row.names(df) <- NULL
  df <- df[, c('seqnames', 'start', 'clusterSource')]
  
  comm <- paste0('gRxCluster(df$seqnames, df$start, df$clusterSource, ',  kvals, ', nperm=', nperm, ', cutpt.tail.expr=', cutpt.tail.expr, ', cutpt.filter.expr=', cutpt.filter.expr, ')')
  scan <- eval(parse(text = comm))
  
  if(length(scan)==0) return(GRanges())
  
  gr1.overlap <- GenomicRanges::findOverlaps(gr1, scan, ignore.strand=TRUE)
  gr2.overlap <- GenomicRanges::findOverlaps(gr2, scan, ignore.strand=TRUE)
  
  if(length(gr1.overlap) > 0 & length(gr2.overlap) == 0){
    scan$clusterSource <- gr1.label
    return(scan)
  }
  
  if(length(gr1.overlap) == 0 & length(gr2.overlap) > 0){
    scan$clusterSource <- gr2.label
    return(scan)
  }
  
  gr1.overlap <- stats::aggregate(queryHits ~ subjectHits, data=as.data.frame(gr1.overlap), FUN=length)
  gr2.overlap <- stats::aggregate(queryHits ~ subjectHits, data=as.data.frame(gr2.overlap), FUN=length)
  
  gr1.list <- list()
  gr1.list[gr1.overlap[,1]] <- gr1.overlap[,2]
  
  gr2.list <- list()
  gr2.list[gr2.overlap[,1]] <- gr2.overlap[,2]
  
  # index error protection
  gr1.list[seq(length(gr1.list)+1, (length(gr1.list) + length(scan)-length(gr1.list) + 1))] <- 0
  gr2.list[seq(length(gr2.list)+1, (length(gr2.list) + length(scan)-length(gr2.list) + 1))] <- 0
  
  scan$clusterSource <- '?'
  
  for(i in 1:length(scan))
  {
    if(! is.null(gr1.list[[i]]) & is.null(gr2.list[[i]])){
      scan[i]$clusterSource <- gr1.label
    } else if (is.null(gr1.list[[i]]) & ! is.null(gr2.list[[i]])){
      scan[i]$clusterSource <- gr2.label
    } else if (gr1.list[[i]] > gr2.list[[i]]){
      scan[i]$clusterSource <- gr1.label
    } else {
      scan[i]$clusterSource <- gr2.label
    }
  }
  
  scan
}


#' Run a number of combination of settings for geneRxCluster scan statistics.
#' Bioinformatics. 2014 Jun 1;30(11):1493-500. 
#' 
#' @return Data frame of clusters statistics from each combination of paramteres tested. 
#'
#' @export
tuneScanStatistics <- function(A, B, CPUs=25, 
                               minWindow=5, maxWindow=100, windowStep=5, 
                               minProb=0.75, maxProb=0.95, probStep=0.025, 
                               minTarget=5, maxTarget=100, targetStep=5){
  
  t <- do.call(rbind, lapply(minWindow:(maxWindow - minWindow), function(i){
    r <- data.frame()
    
    for(i2 in (minWindow+windowStep):maxWindow){
      if(i2-i>=windowStep){
        for(p in seq(minProb, maxProb, by=probStep)){
          for(t in seq(minTarget, maxTarget, by=targetStep)){
            kvals <- paste0(i, 'L:', i2, 'L')
            cutpt.filter.expr <- paste0('apply(x, 2, quantile, probs = ', p, ', na.rm = TRUE)')
            cutpt.tail.expr   <- paste0('critVal.target(k, n, target = ', t, ', posdiff = x)')
            
            r <- rbind(r, data.frame(kvals=kvals,
                                     cutpt.filter.expr=cutpt.filter.expr,
                                     cutpt.tail.expr=cutpt.tail.expr))
          }
        }
      }
    }
    
    r
  }))
  
  t$s <- ceiling(seq_along(1:nrow(t))/(nrow(t)/CPUs))
  save(list = ls(all.names = TRUE), file='tuneData.RData', envir = environment())
  
  cluster <- parallel::makeCluster(CPUs)
  t2 <- do.call(rbind, parallel::parLapply(cluster, split(t, t$s), function(x){
    #t2 <- do.call(rbind, lapply(split(t, t$s), function(x){  
    load('tuneData.RData')
    
    do.call(rbind, lapply(1:nrow(x), function(i){
      message('chunk ', x$s[1], ' parameter row ', i)
      
      tryCatch({
        scan <- gt23::scanStats(A, B,
                                kvals=x[i,]$kvals, 
                                nperm='1000',
                                cutpt.filter.expr=x[i,]$cutpt.filter.expr,
                                cutpt.tail.expr=x[i,]$cutpt.tail.expr)
        
        data.frame(kvals=x[i,]$kvals,
                   cutpt.filter.expr=x[i,]$cutpt.filter.expr,
                   cutpt.tail.expr=x[i,]$cutpt.tail.expr,
                   FDR=gRxSummary(scan)$FDR,
                   totalClusters=gRxSummary(scan)$Clusters_Discovered,
                   Aclusters=length(subset(scan, clusterSource=='A')),
                   Bclusters=length(subset(scan, clusterSource=='B')),
                   avgClusterWidth=(sum(width(scan)) / length(scan)))
      }, 
      error = function(cond) {
        return(data.frame())
      })
    }))
  }))
  stopCluster(cluster)
  t2
}


createIntUCSCTrack <- function(d, abundCuts = c(5,10,50), 
                               posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                               negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                               title = 'intSites', outputFile = 'track.ucsc', visibility = 1,
                               padSite = 0, siteLabel = NA){
  
  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')
  
  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')
  
  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not defined or can not be found in your data.')
  
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  
  # Convert Hex color codes to RGB color codes. 
  # col2rgb() returns a matrix, here we collapse the columns into comma delimited strings.
  #   grDevices::col2rgb(posColors)
  #         [,1] [,2] [,3] [,4]
  #   red    140  103   66   29
  #   green  157  104   52    0
  #   blue   255  227  199  171
  
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s",
                       title, title, visibility)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}


expandTimePoints <- function(tps){
  d <- tibble::tibble(tp = sub('_', '.', tps))
  d$n <- 1:nrow(d)
  
  d$timePointType <- stringr::str_match(base::toupper(d$tp), '[DMY]')
  d$timePointType[which(is.na(d$timePointType))] <- 'X'
  
  d <- dplyr::bind_rows(lapply(split(d, d$timePointType), function(x){
    n <- as.numeric(stringr::str_match(x$tp, '[\\d\\.]+')) * ifelse(grepl('\\-', x$tp), -1, 1)
    
    if(x$timePointType[1] == 'D'){
      x$timePointMonths <- base::round(n / 30.4167, digits = 0)
      x$timePointDays   <- base::round(n, digits = 0)
    } else if(x$timePointType[1] == 'M'){
      x$timePointMonths <- base::round(n, digits = 0)
      x$timePointDays   <- base::round(n * 30.4167, digits = 0)
    } else if(x$timePointType[1] == 'Y'){
      x$timePointMonths <- base::round(n * 12, digits = 0)
      x$timePointDays   <- base::round(n * 365, digits = 0)
    } else {
      message('Warning - could not determine date unit for: ', paste0(unique(x$timePoint), collapse = ', '))
      x$timePointMonths <- n
      x$timePointDays   <- n 
    }
    x
  }))
  
  data.frame(dplyr::arrange(d, n) %>% dplyr::select(timePointMonths, timePointDays))
}