#!/usr/bin/Rscript

library(ggplot2)


# Rscript test_rscript.R --input ../methylation_record_gene_flank_bin20_binfl100.txt 
#   --label gene --nbin 20 --output CG_meth_gene60_flank60 --flank --nbin_flank 100

# Rscript test_rscript.R --input ../methylation_record_promoter2k_bin20.txt,../methylation_record_exon_bin20.txt,../methylation_record_intron_bin20.txt 
#   --label promoter,exon,intron --nbin 20 --output CG_meth_byFeature_nbin20


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  flanking <- FALSE
  fig <- "meth_output"
  fig_height = 10 
  fig_width = 15
  lim <- 0
  #grouping <- "sample"
  
  for(i in 1:length(args)){
    if(args[i] == "--help"){
      cat("USAGE: Rscript MethByBin_Plot.R --input <file1.tab,file2.tab,...> --nbin <INT> --label <lab1,lab2,...>", 
"[--output <outFileName> --flank --nbin_flank <INT> --width --height] \n",
"Required arguments: \n --input \t list of input files separated by ',' and no space\n --nbin \t number of bins to be considered for feature \n --label \t list of labels of thr input files separated by ',' and no space \n",
"Optional arguments: \n --output \t default 'meth_output' \n --flank \t if flanking regions are to be considered \n --nbin_flank \t [required if --flank] number of bins of flanking region \n", 
"--width \t figure width in cm, default 15 \n --height \t figure height in cm, default 10 \n --ylim \t upper limit of y axis. Default ajusted according to max values \n",
"--group \t select group do plot. options: 'sample'[default] or 'coord'\n")
    }
    if(args[i] == "--input"){
      meth <- args[i+1]
      if(grepl(",", meth)) {
        meth <- unlist(strsplit(meth, ","))
        }
    }
    
    if(args[i] == "--flank") {
      flanking  <- TRUE
    }
    if(args[i] == "--nbin_flank") {
      nbin_flank  <- as.numeric(args[i+1])
      adjustXaxis <- 1000/nbin_flank
    }
    if(args[i] == "--nbin"){
      nbin  <- as.numeric(args[i+1])
    }
    if(args[i] == "--output"){
      fig = args[i+1]
    }
    if(args[i] == "--height"){
      fig_height <- as.numeric(args[i+1])
    }
    if(args[i] == "--width"){
      fig_width  <- as.numeric(args[i+1])
    }
    if(args[i] == "--label"){
      lab  <- args[i+1]
      if(grepl(",", lab)) {
        lab <- unlist(strsplit(lab, ","))
        }
    }
    if(args[i] == "--ylim"){
      lim  <- c(0, as.numeric(args[i+1]))
    }
    if(args[i] == "--group"){
      grouping  <- args[i+1]
    }
    
  }
  #print(lim)
  if (!"--input" %in% args) {stop("no input file")}
  if (!"--group" %in% args) {stop("select group to plot means")}
  meth_table <- data.frame()
  for (file in meth) {
    temp_tab <- read.table(file, head = T) 
    meth_table <- rbind(meth_table, temp_tab)
  }
  meth_table <- meth_table[order(meth_table[,"bin"]),]
  min <- min(meth_table$bin)  ## by default the lowest value is -19
  max <- abs(max(meth_table$bin))
  if (grouping == "coord") {
    meth_table$interact <- interaction(meth_table$file, meth_table$coordID)
    p <- ggplot(meth_table, aes(x=bin, y=methylation_level, group=coordID, fill = sample)) +
      geom_line(aes(group=interact),colour="grey75") +
      stat_summary(aes(colour=coordID),fun.y = "mean", geom = "line") +
      expand_limits(y = lim) +
      theme(legend.title=element_blank()) +
      scale_colour_manual(values = c("#7FCDBB","#1D91C0", "#225EA8", "#081D58")) +
      ylab("Methylation level")
    if (flanking) {
      flank1 <- paste( -(min -1)/adjustXaxis, "kb", sep = " ")
      p <- p + scale_x_continuous(breaks=c(min, 1, nbin, max), 
                                  labels=c(paste("-",flank1), "TSS", "TTS", flank1)) +
        xlab(lab) +
        geom_vline(xintercept = c(1, max + min - 1), linetype = "dashed")
    } else {
      p = p + facet_wrap(~feature, scales = "fixed") +
        scale_x_continuous(breaks=c(min, max), 
                           labels=c("5'", "3'")) + xlab("")
    }
    ggsave(paste(fig, ".png"), p, height = fig_height, width = fig_width, unit="cm", dpi=300)
    #png(paste(fig, ".png"), height = fig_height, width = fig_width, units = "cm", res=300)
    #p
    #dev.off()
  } else {
    p <- ggplot(meth_table, aes(x=bin, y=methylation_level, group=sample, fill = sample)) +
      geom_line(aes(group=file),colour="grey75") +
      stat_summary(aes(colour=sample),fun.y = "mean", geom = "line") +
      expand_limits(y = lim) +
      theme(legend.title=element_blank()) +
      ylab("Methylation level")
    if (flanking) {
      flank1 <- paste( -(min -1)/adjustXaxis, "kb", sep = " ")
      p <- p + scale_x_continuous(breaks=c(min, 1, nbin, max), 
                                  labels=c(paste("-",flank1), "TSS", "TTS", flank1)) +
        xlab(lab) +
        geom_vline(xintercept = c(1, max + min - 1), linetype = "dashed")
    } else {
      p = p + facet_wrap(~feature, scales = "fixed") +
        scale_x_continuous(breaks=c(min, max), 
                           labels=c("5'", "3'")) + xlab("")
    }
    ggsave(paste(fig, ".png"), p, height = fig_height, width = fig_width, unit="cm", dpi=300)
    #png(paste(fig, ".png"), height = fig_height, width = fig_width, units = "cm", res=300)
    #p
    #dev.off()
    
  }

  
}

main()



