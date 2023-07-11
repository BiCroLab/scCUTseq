## Author: Luuk Harbers
## Date: 2021-05-16
## Script for plotting bulk or single-cell copynumber profiles

## Load/install packages
packages = c("data.table", "argparser", "RColorBrewer", "scales", "ggdendro", "cowplot", "ggplot2", "tidyverse", "pbapply")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Load in cnv results and plot copy number profiles for bulk or single cell CUTseq")
parser = add_argument(parser, "--rds", help = "Path to rds file of cnv results", nargs = 1)
parser = add_argument(parser, "--runtype", help = "Runtype (either 'single' or 'bulk')", nargs = 1)
parser = add_argument(parser, "--threads", help = "Number of threads to use", nargs = 1, type = "integer")
parser = add_argument(parser, "--outdir", help = "Path to output directory", nargs = 1)
argv = parse_args(parser)

# argv = list()
# argv$rds = "/mnt/AchTeraD/data/scCUTseq/prostate/P3/cnv/500000/cnv.rds"
# argv$runtype = "single"
# argv$outdir = "/mnt/AchTeraD/data/scCUTseq/prostate/P3/cnv/500000/plots/"
# argv$threads = 20

# Set parameters
threads = argv$threads

# Set cowplot theme as default theme
theme_set(theme_cowplot())

# # 
# # Make outdir if not exists
# if(!dir.exists(argv$outdir)) {
#   dir.create(argv$outdir)
#   dir.create(paste0(argv$outdir, "/profiles"))
#   dir.create(paste0(argv$outdir, "/genomewide"))
# }

# Read in rds
cnv = readRDS(argv$rds)

# Prepare bins for plotting (same for single cell and bulk)
bins = cnv$bins
bins[, bin := seq_along(chr)]
bins[, end_cum := cumsum((end - start) + 1)]
bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

# Colors for single-cell sequencing
if (argv$runtype == 'single') {
  colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
             "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
  names(colors) = c(as.character(0:10), "10+")
} else {
  colors = c(brewer.pal(3, "Set1")[1:2], "#c1c1c1")
  names(colors) = c("Amplification", "Deletion", "Neutral")
}

# Get sample names
samples = colnames(cnv$counts_gc)

# Plot individual profiles
if (argv$runtype == 'single'){
  cat("Plotting individual profiles...\n")
  pblapply(samples, function(id) {
    dt = cbind(bins, cnv$copynumber[[id]], cnv$counts_gc[[id]] * cnv$ploidies[sample == id, ploidy])
    setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "cn", "raw"))
    
    dt[, cn := ifelse(cn < 11, cn, 11)]
    dt[, col := ifelse(cn < 11, as.character(cn), "10+")]
    dt[, col := factor(col, levels = c(as.character(0:10), "10+"))]
    
    # Sample stats
    stat_string = paste0(id,
                         " | reads: ", cnv$stats[sample == id, total_reads],
                         " | avg reads/bin: ", round(cnv$stats[sample == id, mean_reads]),
                         " | spikiness: ", round(cnv$stats[sample == id, spikiness], 3),
                         " | RF classifier: ", cnv$stats[sample == id, classifier_prediction],
                         " | segmentation type: ", cnv$segmentation_type)
    
    # save plot
    plot = ggplot(dt, aes(x = bin)) +
      geom_point(aes(y = raw, color = col), size = 0.7) +
      geom_point(aes(y = cn), size = 1) +
      scale_color_manual(values = colors, drop = F) +
      scale_y_continuous(labels=comma_format(accuracy = 1), breaks = pretty_breaks(6), limits = c(0, 12)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(y = "Copy Number", x = "", subtitle = stat_string) +
      geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
      geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", inherit.aes = F) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # Save plot
    ggsave(filename = paste0(argv$outdir, "/profiles/", id, ".png"), plot = plot,
           width = 14, height = 7, units = "in", dpi = 300)
  }, cl = threads)
  
  # Plot heatmap
  dt = cbind(bins, cnv$copynumber)
  
  # Make dendrogram
  cat("Calculate distances and clustering samples...\n")
  hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
  dhc = as.dendrogram(hc)

  # Rectangular lines
  ddata = dendro_data(dhc, type = "rectangle")

  # Plot Dendrogram
  dendro = ggplot(ggdendro::segment(ddata)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0.004, 0.004)) +
    theme_dendro()
  
  # Prepare for heatmap
  dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
  dt_melt[, value := factor(value)]
  dt_melt[as.numeric(value) > 10, value := "10+"]
  dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
  
  # Set sample order
  dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
  
  # Plot heatmap 20K+
  heatmap = ggplot(dt_melt) +
    geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 2) +
    coord_flip() +
    scale_color_manual(values = colors, drop = F) +
    labs(color = "Copy Number") + 
    scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
    geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank())
  
  # Combine plots and save
  #combined = plot_grid(dendro, heatmap,  align = "h", rel_widths = c(0.2, 2), ncol = 2)
  cat("Plotting Genomewide heatmap...\n")
  ggsave(filename = paste0(argv$outdir, "/genomewide/genomewideheatmap.png"), plot = heatmap,
         width = 20, height = 28, dpi = 900)
} else {
  pblapply(samples, function(id) {
    cat("Plotting individual profiles...\n")
    dt = cbind(bins, cnv$segments[[id]], cnv$counts_lrr[[id]])
    setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "segment", "raw"))
    
    # Assign colors
    dt[segment > log2(2.5/2), col := "Amplification"]
    dt[segment < log2(1.5/2), col := "Deletion"]
    dt[is.na(col), col := "Neutral"]
    

    # Sample stats
    stat_string = paste0(id,
                         " | reads: ", cnv$stats[sample == id, total_reads],
                         " | avg reads/bin: ", round(cnv$stats[sample == id, mean_reads]),
                         " | spikiness: ", round(cnv$stats[sample == id, spikiness], 3),
                         " | segmentation type: ", cnv$segmentation_type)
    
    # save plot
    plot = ggplot(dt, aes(x = bin)) +
      geom_point(aes(y = raw, color = col), size = 0.7) +
      geom_point(aes(y = segment), size = 1) +
      scale_color_manual(values = colors, drop = F) +
      scale_y_continuous(limits = c(-2, 4), labels=comma_format(accuracy = 1), breaks = pretty_breaks(6)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(y = "Log2 ratio", x = "", subtitle = stat_string) +
      geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
      geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", inherit.aes = F) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # Save plot
    ggsave(filename = paste0(argv$outdir, "/profiles/", id, ".png"), plot = plot,
           width = 14, height = 7, units = "in", dpi = 300)
  }, cl = threads)
  # Plot heatmap
  dt = cbind(bins, cnv$segments)
  
  # Make dendrogram
  cat("Calculate distances and cluster all samples...\n")
  hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
  dhc = as.dendrogram(hc)
  
  # Rectangular lines
  ddata = dendro_data(dhc, type = "rectangle")
  
  # Plot Dendrogram
  dendro = ggplot(ggdendro::segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() + 
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    theme_dendro()
  
  # Prepare for heatmap
  dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
  
  # Set sample order
  dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
  
  # Plot heatmap 20K+
  heatmap = ggplot(dt_melt) +
    geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 2) +
    coord_flip() +
    scale_color_gradient2(midpoint = 0, low = "blue", high = "red", limits = c(-2, 4), oob = scales::squish) +
    labs(color = "Log2 ratio") + 
    scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
    geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank())
  
  # Combine plots and save
  combined = plot_grid(dendro, heatmap,  align = "h", rel_widths = c(0.2, 2), ncol = 2)
  cat("Plotting Genomewide heatmap...\n")
  ggsave(filename = paste0(argv$outdir, "/genomewide/genomewideheatmap.png"), plot = combined,
         width = 16, height = 4, dpi = 900)
}
