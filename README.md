# parsing-BSseqData

This module includes python script and an R script that result from an adaptation of MethOverRegion module from ViewBS (Huang X et al. 2018. Bioinformatics), to generate metaplots for methylation frequency over given genomic features, normalized for a given bin length.

This adaptation was created to generate similar metaplots, for feature and flanking regions, but also includes the option of ploting several features at in the same metaplot (e.g. promoter, exon, intron and intergenic regions).


MethByBin.py (python2.7) takes as input: 
- the genome feature coordinates (use gff2tab.pl https://github.com/gpertea/gscripts/blob/master/gff2tab.pl ) 
- a .txt file containing the path to each Genome-wide cytosine methylation report (from Bismark) and the corresponding sample name, separated by tab and one entry per line: e.g. \<path/to/file.bgz> \<sampleName>
- REQUIRES tabix and argparse

```
usage: MethByBin.py [-h] --featCoord str [-listCoord] --feature str --sample
                    str [--bin INT] [--out STR] [--min_depth INT]
                    [--max_depth INT] [-flanks] [--fl_length INT]
                    [--bin_fl INT]


Get methylation frequencies per bin for a given genomic feature (and flanking
regions)

optional arguments:
  -h, --help       show this help message and exit
  --featCoord str  [REQUIRED] file name, or list of paths to files (one per
                   line),containing the coordinates to the feature to be
                   considered. - File containing a list of paths (use argument
                   -listCoord) should have the format: <path> <sampleID>-
                   File(s) with coordinates should have the following format
                   <chr> <feature> <stt> <end> <strand> <id> <product> (e.g.
                   output of gff2tab.pl)
  -listCoord       Use this argument if --featCoord contains a list of paths
                   for coordinates. May be usefull when different groups of
                   the same features
  --feature str    name of the feature being considered. e.g. transcript,
                   exon, intron, ...
  --sample str     [REQUIRED] list of files and corresponding sample names as:
                   <path/to/file.bgz> <sampleName>
  --bin INT        the feature will be divided in this many bins. If flanks =
                   False, features with size < bin number will not be
                   considered. Default 60.
  --out STR        Output file name
  --min_depth INT  minimum mC base coverage to be considered. Default: 5
  --max_depth INT  maximum mC base coverage to be considered. Default: 1000000
  -flanks          if selected, the analysis will include the regions flanking
                   the feature (flanks = True)
  --fl_length INT  if flanks, size in bp of the flanking region to be
                   considered. Default: 2000.
  --bin_fl INT     if flanks, the flanking regions will be divided in this
                   many bins. Default 100.


```

The output of this script will be a methylation report by bin, for all samples, in the corresponding format:

\<file>	\<sample>	\<feature> \<coordID>	\<bin>	\<C_number>	\<T_number>	\<methylation_level>
  
This output is then used as input in MethByBin_Plot.R to plot the data (requires ggplot2 package):
```
Usage: Rscript MethByBin_Plot.R --input <file1.tab,file2.tab,...> --nbin <INT> --label <lab1,lab2,...> 
               [--output <outFileName>] [--flank --nbin_flank <INT>] [--width] [--height] 

Required arguments: 
  --input       list of input files separated by ',' and no space
  --nbin        number of bins to be considered for feature 
  --label       list of labels of thr input files separated by ',' and no space
  --group       select group to plot means. Options: 'sample' or 'coord'

Optional arguments:
  --output      default 'meth_output' 
  --flank       if flanking regions are to be considered 
  --nbin_flank  [required if --flank] number of bins of flanking region 
  --width       figure width in cm, default 15 
  --height      figure height in cm, default 10 
```

If methylation report includes flanking regions, option *--flank* should be used together with *--nbin_flank <INT>*.
  Bin number in *--nbin* (and *nbin_flank*) should match the number of bins considered for MethByBin.py. 
