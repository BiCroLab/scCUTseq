#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2021-05-12
## Script to perform cnv calling from a list of binned read counts

## Load/install packages
packages = c("data.table", "argparser", "DNAcopy", "ParDNAcopy", "pbapply", "gtools", "randomForest", "copynumber", "ASCAT.sc", "aCGH", "matrixStats")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Filter bincounts with the use of a blacklist and correct counts for GC content")
parser = add_argument(parser, "--counts", help = "Path to counts", nargs = 1)
parser = add_argument(parser, "--bins", help = "Path to bin file", nargs = 1)
parser = add_argument(parser, "--blacklist", help = "Path to blacklist file", nargs = 1)
parser = add_argument(parser, "--gc", help = "Path to GC content file", nargs = 1)
parser = add_argument(parser, "--binsize", help = "Binsize to run with", nargs = 1, type = "numeric")
parser = add_argument(parser, "--normseg", help = "Normalize additional segments", nargs = 1, type = "character")
parser = add_argument(parser, "--removethreshold", help = "Remove samples with fewer reads than this (0 to not remove any samples)", nargs = 1, type = "numeric", default = 3e5)
parser = add_argument(parser, "--segmentation", help = "Type of segmentation (single or joint)", nargs = 1, type = "character")
parser = add_argument(parser, "--alpha", help = "alpha parameter for segmantation with DNACopy", nargs = 1, type = "numeric")
parser = add_argument(parser, "--prune", help = "undo.prune parameter for segmentation with DNAcopy", nargs = 1, type = "numeric")
parser = add_argument(parser, "--penalty", help = "Penalty value (gamma) for joint segmentation with mpcf from copynumber", nargs = 1, type = "numeric", default = 10)
parser = add_argument(parser, "--type", help = "Type of sequencing 'single' or 'bulk'", nargs = 1, type = "character")
parser = add_argument(parser, "--randomforest", help = "Path to randomforest model for single cell classification", 
                      default = NULL, nargs = 1, type = "character")
parser = add_argument(parser, "--rfthreshold", help = "threshold of randomforest model for single cell classification", 
                      nargs = 1, type = "character")
parser = add_argument(parser, "--minploidy", help = "Min ploidy of sample (single-cell only)", 
                      nargs = 1, type = "integer", default = 1.5)
parser = add_argument(parser, "--maxploidy", help = "Max ploidy of sample (single-cell only)", 
                      nargs = 1, type = "integer", default = 6)
parser = add_argument(parser, "--minpurity", help = "Min purity of sample (1 for single cell / pure cell lines)", 
                      nargs = 1, type = "integer", default = 1.5)
parser = add_argument(parser, "--maxpurity", help = "Max purity of sample", 
                      nargs = 1, type = "integer", default = 6)
parser = add_argument(parser, "--sex", help = "Sex of sample (either male or female)", nargs = 1, type = "character")
parser = add_argument(parser, "--threads", help = "Number of threads to use", nargs = 1, type = "integer")
parser = add_argument(parser, "--output", help = "Path to output file", nargs = 1)
argv = parse_args(parser)

# argv = list()
# argv$counts = "/mnt/AchTeraD/data/scCUTseq/prostate/P3_subs/cnv/500000/all-counts.tsv.gz"
# argv$gc = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copynumber-pipeline/cnv-calling/files/hg19/GC_variable_500000_150_bwa"
# argv$blacklist = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copynumber-pipeline/cnv-calling/files/hg19/hg19-blacklist.v2_adjusted.bed"
# argv$bins = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copynumber-pipeline/cnv-calling/files/hg19/variable_500000_150_bwa.bed"
# argv$norm = "/mnt/AchTeraD/Documents/Projects/scCUTseq/copynumber-pipeline/cnv-calling/files/hg19/normalize.tsv"
# argv$binsize = 500000
# argv$alpha = 0.0001
# argv$undo.prune = 0.05
# argv$penalty = 6
# argv$type = "single"
# argv$segmentation = "joint"
# argv$minploidy = 1.7
# argv$maxploidy = 3.5
# argv$minpurity = 1
# argv$maxpurity = 1
# argv$sex = "male"
# argv$threads = 32
# argv$output = ""
# argv$randomforest = "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/randomforest.rds"
# argv$rfthreshold = 0.3
# argv$removethreshold = 300000

# Check input parameters
if(!file.exists(argv$counts)) {
  stop("Error: --type argument has to be either 'single' or 'bulk'")
}
if(!file.exists(argv$gc)) {
  stop("Error: --gc file does not exist")
}
if(!file.exists(argv$blacklist)) {
  stop("Error: --blacklist file does not exist")
}
if(!file.exists(argv$bins)) {
  stop("Error: --bins file does not exist")
}
if(!argv$type %in% c("single", "bulk")) {
  stop("Error: --type argument has to be either 'single' or 'bulk'")
}

# Make list of all data that needs to be included in final output
out = list()
out$counts = fread(argv$counts)
out$removethreshold = argv$removethreshold
# Remove low read samples
out$counts = out$counts[, colnames(out$counts)[colSums(out$counts, na.rm = T) >= out$removethreshold], with = F]
out$gc = fread(argv$gc)
out$blacklist = fread(argv$blacklist)
out$bins = fread(argv$bins)
out$penalty = argv$penalty
out$alpha = argv$alpha
out$undo.prune = argv$prune
out$segmentation_type = argv$segmentation

# Set parameters
type = argv$type
binsize = argv$binsize
samples = colnames(out$counts)
minploidy = argv$minploidy
maxploidy = argv$maxploidy
minpurity = argv$minpurity
maxpurity = argv$maxpurity
threads = argv$threads

if(argv$sex == "female") {
  ismale = FALSE
} else if(argv$sex == "male") {
  ismale = TRUE
} else {
  stop("Error: --sex must be either 'male' or 'female'")
}


output = argv$output
if(!is.na(argv$randomforest)) rf = readRDS(argv$randomforest)
rfthreshold = argv$rfthreshold



## Define functions
# LOWESS regression for GC normalization
lowess.gc = function(x, y) {
  low = lowess(x, log(y), f=0.05);
  z = approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}

# Set names
setnames(out$bins, c("chr", "start", "end"))
setnames(out$blacklist, c("chr", "start", "end", "reason"))

# Only keep large blacklisted regions
out$blacklist = out$blacklist[(end - start) > 0.2 * binsize,]

# Get overlaps with blacklist
bins = copy(out$bins) # Make copy, do not assign by reference
blacklist = copy(out$blacklist) # Make copy, do not assign by reference
setkey(bins, chr, start, end)
setkey(blacklist, chr, start, end)

# Get the overlaps and retain bins that do not overlap
overlaps = foverlaps(bins, blacklist)
# Reorder
overlaps = overlaps[mixedorder(overlaps$chr)]
overlaps = unique(overlaps, by = c("chr", "i.start", "i.end"))
indices = which(is.na(overlaps$start))

# Filter reads by indices
out$counts = out$counts[indices]
out$gc = out$gc[indices]
out$bins = out$bins[indices]

# Add pseudocount to readcounts
out$counts = out$counts + 1

# Normalize by mean and GC correction
message("Running GC normalization...")
out$counts_gc = out$counts[, pblapply(.SD, function(sample) {
  lowess.gc(out$gc$V1, (sample + 1 / mean(sample + 1)))
}, cl = threads)]

# Add 0.001 to 0 ratios and make log ratios
# out$counts_gc[out$counts_gc == 0] = 1e-3
out$counts_lrr = log2(out$counts_gc)

# Normalize additional segments if requested
if (file_test("-f", argv$norm)) {
  message("Running Segment normalization...")
  normal = fread(argv$norm)
  setnames(normal, c("chr", "start", "end", "adjust"))
  normal[, chr := as.character(chr)]
  setkey(normal, chr, start, end)
  setkey(out$bins, chr, start, end)
  
  # Get indices to normalize
  toCorrect = foverlaps(normal, out$bins)
  toCorrect = merge(out$bins, toCorrect[, .(chr, start, end, adjust)], by = c("chr", "start", "end"), all.x = T)
  toCorrect[is.na(adjust), adjust := 0]
  
  # Get unique
  toCorrect = unique(toCorrect, by = c("chr", "start", "end"))
  toCorrect = toCorrect[mixedorder(toCorrect$chr)]
  
  # Reorder bins to mixedorder
  out$bins = out$bins[mixedorder(out$bins$chr)]
  
  # Normalize the segments
  out$counts_lrr = out$counts_lrr - toCorrect$adjust
  out$counts_gc = 2^out$counts_lrr
} else {
  message("Warning, normalization file does not exist, skipping...")
}


if (out$segmentation_type == "single") {
  message("Segmenting profiles using CBS...")
  # Make CNA object, do smoothing and segmentation
  cna = CNA(out$counts_lrr, out$bins$chr, out$bins$start, data.type = "logratio", sampleid = colnames(out$counts_lrr)) 
  cna_smooth = smooth.CNA(cna)
  
  
  cna_segment = parSegment(cna_smooth, alpha = out$alpha, min.width = 5, undo.splits = "prune", undo.prune = out$undo.prune,
                           njobs = threads, distrib = "Rparallel")
  
  
  # Restructure
  out$segments_long = data.table(cna_segment$output)
  out$segments_long = out$segments_long[mixedorder(out$segments_long$chrom)]
  setorder(out$segments_long, ID)
  out$segments_long[, ID := gsub("X", "", ID)]
  
  # Add segments to output list
  segments_rep = out$segments_long[, .(chr = rep(chrom, num.mark), start = rep(loc.start, num.mark), 
                                       seg.mean = rep(seg.mean, num.mark)), by = .(ID)]
  #segments_rep = segments_rep[mixedorder(segments_rep$chr)]
  #setorder(segments_rep, ID)
  segments_rep[, bin_num := rep(seq_len(nrow(out$counts_lrr)), ncol(out$counts_lrr))]
  out$segments = dcast(segments_rep, bin_num ~ ID, value.var = "seg.mean")[, -1]
  
  # Reorder out$segments to match sample orders of previous analysis
  setcolorder(out$segments, colnames(out$counts_lrr))
}
if (out$segmentation_type == "joint") {
  message("Segmenting profiles using multipcf...")
  # Multipcf
  mpcf_input = as.data.frame(cbind(out$bins[, 1:2], out$counts_lrr))
  mpcf = multipcf(mpcf_input, gamma = out$penalty, normalize = FALSE, verbose = FALSE)
  
  # Restructure
  out$segments_long = melt(data.table(mpcf), 
                           id.vars = c("chrom", "start.pos", "end.pos", "arm", "n.probes"), variable.name = "ID")
  setnames(out$segments_long, c("chrom", "start.pos", "end.pos", "arm", "num.mark", "ID", "seg.mean"))
  mpcf_rep = out$segments_long[, .(chr = rep(chrom, num.mark), 
                                   start = rep(start.pos, num.mark), 
                                   end = rep(end.pos, num.mark), 
                                   seg.mean = rep(seg.mean, num.mark)), by = ID]
  
  # Set order
  #setorder(mpcf_rep, ID, chr, start)
  mpcf_rep[, bin_num := seq_len(.N), by = ID]
  out$segments = dcast(mpcf_rep, bin_num ~ ID, value.var = "seg.mean")[, -1]
  
  # Reorder out$segments to match sample orders of previous analysis
  setcolorder(out$segments, colnames(out$counts_lrr))
}


#### If including mergeLevels need to add aCGH package and need to adjust following segments to include merged segments instead of non-merged segments
# MergeLevels
message("Merging similar segments...")
segments_merged = pblapply(colnames(out$segments), function(cell) {
  # # Run mergeLevels and transform to data.table
  # merged = mergeLevels(out$counts_lrr[[cell]],
  #                      out$segments[[cell]],
  #                      pv.thres = 1e-10,
  #                      verbose = 0)$vecMerged |> data.table()
  
  # Run mergeLevels per chromosome and transform to data.table
  merged = lapply(unique(out$bins$chr), function(chr) {
    indices = which(out$bins$chr == chr)
    mergeLevels(out$counts_lrr[[cell]][indices],
                out$segments[[cell]][indices],
                pv.thres = 1e-5,
                thresMax = .3,
                verbose = 0)$vecMerged |> data.table()
  })
  merged = rbindlist(merged)
  
  
  setnames(merged, cell)
  return(merged)
  
}, cl = threads)
out$segments_merged = do.call(cbind, segments_merged)

# # Modify log ratio segments to contain median normalized read counts per bin
# median_segments = pblapply(colnames(out$counts_gc), function(cell) {
#   # Subset segments_long to only have correct sample
#   dt = data.table(num = out$segments_long[ID == cell]$num.mark)
#   # Calculate start and end bin for segments
#   dt[, end := cumsum(num)]
#   dt[, start := data.table::shift(end, fill = 0) + 1]
# 
#   medians = sapply(1:nrow(dt), function(x) {
#     return(median(out$counts_gc[[cell]][dt[x, start]:dt[x, end]]))
#   })
#   return(data.table(rep(medians, dt$num)))
# }, cl = threads)
# out$segments_read = do.call(cbind, median_segments)
# setnames(out$segments_read, colnames(out$counts_gc))

# Modify log ratio segments to contain median normalized read counts per bin 
#### IN CASE OF MERGELEVELS
median_segments = pblapply(colnames(out$counts_gc), function(cell) {
  # Get split points
  splits = c(1, which(out$segments_merged[[cell]] != data.table::shift(out$segments_merged[[cell]])),
             nrow(out$segments_merged) + 1)
  
  dt = data.table(num = diff(splits))
  dt[, end := cumsum(num)]
  dt[, start := data.table::shift(end, fill = 0) + 1]
  
  
  
  medians = sapply(seq_len(nrow(dt)), function(x) {
    return(median(out$counts_gc[[cell]][dt[x, start]:dt[x, end]]))
  })
  final = data.table(rep(medians, dt$num))
  setnames(final, cell)
  
  return(final)
}, cl = threads)

out$segments_read = do.call(cbind, median_segments)
setnames(out$segments_read, colnames(out$counts_gc))

# Renormalize
# out$segments_read = data.table(sweep(out$segments_read, 2, colMeans(out$segments_read), '/'))

# Perform integer copy number calling if it's single cell data
if (type == "single") {
  message("Calculating integer copy numbers...")
  # # Initialize CN inference variables for SoS method -- Code adapted from ginkgo pipeline
  # n_cells = ncol(out$segments_read)
  # ploidy = rbind(c(0,0), c(0,0))
  # CNgrid = seq(minploidy, maxploidy, by = 0.05)
  # n_ploidy = length(CNgrid)  # Number of ploidy tests during CN inference
  # CNmult = matrix(0, n_ploidy, n_cells)
  # CNerror = matrix(0, n_ploidy, n_cells)
  # outerColsums = matrix(0, n_ploidy, n_cells)
  # CN = numeric(n_cells)
  # 
  # # Loop through cells and get integer copy number
  # res = pbsapply(1:n_cells, function(cell) {
  #   outerRaw = as.matrix(out$segments_read[, cell, with = F]) %o% CNgrid
  #   outerRound = round(outerRaw)
  #   outerDiff = (outerRaw - outerRound) ^ 2
  #   outerColsums[, cell] <<- colSums(outerDiff, na.rm = FALSE, dims = 1) ## Assignment outside scope
  #   CNmult[, cell] <<- CNgrid[order(outerColsums[, cell])] ## Assignment outside scope
  #   CNerror[, cell] <<- round(sort(outerColsums[, cell]), digits=2) ## Assignment outside scope
  #   
  #   # Define CN multiplier and assign to CN vector
  #   CN[cell] <<- CNmult[1, cell] ## Assignment outside scope
  #   # return integer copy number profile
  #   return(out$segments_read[, cell, with = F] * CN[cell])
  # })
  
  ploidies = seq(minploidy, maxploidy, .01)
  purities = seq(minpurity, maxpurity, .01)
  
  # # Loop through segments and get lowest ploidy error
  # res_ploidies = pblapply(unique(out$segments_long$ID), function(cell) {
  #   dist_mat = buildDistanceMatrix(out$segments_long[ID == cell, seg.mean], out$segments_long[ID == cell, num.mark], purs = purities, ploidies = ploidies, maxTumourPhi = 8, ismale = ismale)
  #   mins = arrayInd(which.min(dist_mat), dim(dist_mat))
  #   ploidy = ploidies[mins[2]]
  #   purity = purities[mins[1]]
  # 
  #   return(data.table(sample = cell,
  #                     ploidy = ploidy,
  #                     purity = purity))
  # }, cl = threads)
  
  res_ploidies = pblapply(colnames(out$segments_merged), function(cell) {
    dist_mat = buildDistanceMatrix(out$segments_merged[[cell]], 
                                   rep(1, nrow(out$segments_merged)), 
                                   purs = purities, ploidies = ploidies, maxTumourPhi = 8, ismale = ismale)
    mins = arrayInd(which.min(dist_mat), dim(dist_mat))
    ploidy = ploidies[mins[2]]
    purity = purities[mins[1]]
    
    return(data.table(sample = cell,
                      ploidy = ploidy,
                      purity = purity))
  }, cl = threads)
  
  out$ploidies = rbindlist(res_ploidies)
  # 
  # # Setnames to CN
  # names(CN) = colnames(out$segments_read)
  # out$cn = CN
  
  # Scale counts_lrr
  res = pblapply(colnames(out$counts_lrr), function(cell) {
    out$counts_lrr[[cell]] * out$ploidies[sample == cell, ploidy]
    }, cl = threads)
  out$counts_lrr_scaled = data.table(do.call(cbind, res))
  setnames(out$counts_lrr_scaled, colnames(out$counts_lrr))
  
  # Save to out
  res = pblapply(colnames(out$segments_read), function(cell) {
    out$segments_read[[cell]] * out$ploidies[sample == cell, ploidy]
    }, cl = threads)
  out$segments_scaled = data.table(do.call(cbind, res))
  setnames(out$segments_scaled, colnames(out$segments_read))
  out$copynumber = round(out$segments_scaled)
  
  # Write stats
  message("Calculating statistics...")
  
  suppressWarnings({
    out$stats = data.table(sample = samples)
    out$stats[, total_reads := colSums(out$counts[, ..samples])]
    out$stats[, mean_reads := colMeans(out$counts[, ..samples])]
    out$stats[, median_reads := out$counts[, ..samples] |> as.matrix() |> colMedians()]
    out$stats[, spikiness := out$counts[, ..samples] |> as.matrix() |> diff() |> abs() |> colSums() / total_reads]
    out$stats[, non_integerness := (out$copynumber[, ..samples] - out$segments_scaled[, ..samples]) |> abs() |> as.matrix() |> colMedians()]
    out$stats[, bin_to_medians := (out$segments_scaled[, ..samples] - out$counts_lrr_scaled[, ..samples]) |> abs() |> as.matrix() |> colMedians()]
    out$stats[, bin_to_integer := (out$copynumber[, ..samples] - out$counts_lrr_scaled[, ..samples]) |> abs() |> as.matrix() |> colMedians()]
    out$stats[, coef_variation := (out$counts_lrr[, ..samples] |> as.matrix() |> colSds(na.rm = TRUE)) / colMeans(out$counts_lrr[, ..samples])]
    
    # Calculate autocorrelation
    out$stats[, autocorrelation := 
                pbsapply(sample, function(x) {
                  tail(acf(out$counts_lrr_scaled[[x]], 1, na.action = na.pass, type = "correlation", plot = FALSE)$acf, 1)
                }, cl = threads)]
  
    out$stats[, mean_absolute_deviation := out$counts_lrr[, ..samples] |> as.matrix() |> colMads(constant = 1)]
    out$stats[, mean_variance := out$copynumber[, ..samples] |> var() |> colMeans()]
              
    # Calculate halfiness
    halfiness = (-log2(abs(pmin(abs(out$segments_scaled[, ..samples] - out$counts_lrr_scaled[, ..samples]), 0.499) - 0.5))) - 1
    out$stats[, total_halfiness := colSums(halfiness, na.rm = TRUE)]
    out$stats[, scaled_halfiness := colSums(halfiness / (out$counts_lrr_scaled[, ..samples] + 1), na.rm = TRUE)]
    out$stats[, breakpoints := out$segments_long[, 1:3] |> unique() |> nrow()]
    out$stats[, mean_cn := colMeans(out$copynumber[, ..samples])]
    out$stats[, mode_cn := sapply(samples, function(x) names(sort(-table(out$copynumber[[x]]))[1]))]  
    
    # Run randomforest classification
    prediction = as.data.table(predict(rf, out$stats, type="prob"))
    out$stats[, classifier_prediction := ifelse(prediction$good >= rfthreshold, "good", "bad")]
  })
            
} else {
  # Write stats
  message("Calculating statistics...")
  
  # Write stats for bulk
  suppressWarnings({
    out$stats = data.table(sample = samples)
    out$stats[, total_reads := colSums(out$counts[, ..samples])]
    out$stats[, mean_reads := colMeans(out$counts[, ..samples])]
    out$stats[, median_reads := out$counts[, ..samples] |> as.matrix() |> colMedians()]
    out$stats[, spikiness := out$counts[, ..samples] |> as.matrix() |> diff() |> abs() |> colSums() / total_reads]
    out$stats[, coef_variation := (out$counts_lrr[, ..samples] |> as.matrix() |> colSds(na.rm = TRUE)) / colMeans(out$counts_lrr[, ..samples])]
    out$stats[, mean_absolute_deviation := out$counts_lrr[, ..samples] |> as.matrix() |> colMads(constant = 1)]
    out$stats[, breakpoints := nrow(unique(out$segments_long[, 1:3]))]
  })
}

# Write output
message("Writing output...")
saveRDS(out, output)
message("Finished.")
