
# General workflow for whole genome methylation analysis:

<br/>

## Get the coordinates for the feature of interest, depending on the plots to generate:

[gff2tab](https://github.com/gpertea/gscripts/blob/master/gff2tab.pl) script can be used to get gene, exon and intron coordinates:


```{bash eval=FALSE, include=TRUE}

gff2tab.pl -T genomeSequence-annot.gff > genomeSequence-transcriptGFFCoord.tab
gff2tab.pl -E -U genomeSequence-annot.gff > genomeSequence-exonUniqueGFFCoord.tab
gff2tab.pl -I -U genomeSequence-annot.gff > genomeSequence-intronUniqueGFFCoord.tab

```


For promoter regions I did a costum sript to collect the coordinates, based on a given size: [gff2promoter](https://github.com/pedro-mb/python_parsers/blob/master/gff2promoter.py)

```{bash eval=FALSE, include=TRUE}

python gff2promoter.py --gff genomeSequence-annot.gff --size 2000 --out genomeSequence-promoter2kGeneGFFCoord.tab 

```


To get intergenic regions I used Bedtools. 'bedtools complement' requires a gff file only with gene rows and a tab separated file (*scaffoldSizes.tsv*) containing all scaffoldIDs and corresponding size in bp:

For example:


  scaffold1  &nbsp; 22220034 <br/>
  scaffold2  &nbsp; 30300303 <br/>
  ... <br/>
<br/>

```{bash eval=FALSE, include=TRUE}
awk '$3== "gene"' genomeSequence-annot.gff > genomeSequence-gene.gff

/data/bin/bedtools2-master/bin/bedtools complement -i genomeSequence-gene.gff -g GenericData/scaffoldSizes.tsv > genomeSequence-intergenicGFFcoord.tab


```

the file with intergenic coordinates need to be further edited to get the same format as the files for the other features. Back then I did this in R ... but I believe that are faster ways to get it:

```{r eval=FALSE, include=TRUE}
getwd()
coord <- read.table("../metilation/genomeSequence-intergenicGFFcoord.tab")
feat <- rep("intergenic", nrow(coord))
strand <- rep(".", nrow(coord))
coord2<- data.frame(coord$V1, feat, coord$V2, coord$V3, strand)

write.table(coord2, file = "genomeSequence-intergenicGFFcoord_edit.tab", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

```

<br/>

## Convert CX_report files and create an index file: 


```{bash eval=FALSE, include=TRUE}

# if report is in '.gz' format, it needs to be converted to '.bgz' before creating the index with tabix. Use the following command or skip to the next

gunzip -c sample1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz | bgzip  > sample1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz

tabix -p vcf sample1_bismark_bt2_pe.deduplicated.CpG_report.txt.bgz


# if you have multiple files 

for i in $(ls -1 */*.txt.gz); do prefix=${i%.gz}; echo $prefix; gunzip -c $i | bgzip  > $prefix.bgz; done

for i in $(ls -1 */*.txt.bgz); do echo $i; tabix -p vcf $i; done

# an index file with extension .tbi is generated 

```

<br/>

## Create a text file containing all BGZ files 

Add full PATH and the corresponding sample group, in the same line separated by tab: 

For example (listSampleBGZgiles.txt):

  path/to/BGZfiles/tissueA_rep1.bgz &nbsp; tissueA <br/>
  path/to/BGZfiles/tissueA_rep2.bgz &nbsp; tissueA <br/>
  path/to/BGZfiles/tissueB_rep1.bgz &nbsp; tissueB <br/>
  path/to/BGZfiles/tissueB_rep2.bgz &nbsp; tissueB <br/>
  ... <br/>

<br/>

## Example1: For a combined plot with promoter, exon and intron features and all tissues/samples: **

for the *MethByBin.py*, the pytabix module is required [link](https://pypi.org/project/pytabix/) . <br/>
for the *MethByBin_Plot.R*, ggplot2 is required.

```{bash eval=FALSE, include=TRUE}

# generate a report containg methylation frequency by genomic "bins" (a genomic window of 20 bp by default, but can be changed) for each .bgz file
# this report will be saved with the name specified in -out
# if necessary, the context can be specified (CG, CHG or CHH)


python MethByBin.py --featCoord genomeSequence-promoter2kGeneGFFCoord.tab --sample listSampleBGZgiles.txt --feature promoter --bin 20 --context CG --out methylation_record_promoter2k_bin20.txt

python MethByBin.py --featCoord genomeSequence-intronUniqueGFFCoord.tab --sample listSampleBGZfiles.txt --feature intron --out methylation_record_intron20.txt --bin 20 --context CG

python MethByBin.py --featCoord genomeSequence-exonUniqueGFFCoord.tab --sample listSampleBGZfiles.txt --feature exon --out methylation_record_exon20.txt --bin 20 --context CG

# use the multiple output files as input for the R script that will generate the plots

Rscript MethByBin_Plot.R --input methylation_record_promoter2k_bin20.txt,methylation_record_exon20.txt,methylation_record_intron20.txt --label promoter,exon,intron --nbin 20 --output Example1 --group sample
```
![Example 1](https://github.com/pedro-mb/parsing-BSseqData/blob/master/tutorial/Example1.png)
<br/>

## Example2: For gene body plot with flanking regions and all tissues/samples:

use the *transcriptGFFCoord.tab* file..The length of the bins can be changed, but they need to be the same in the python and R scripts

```{bash eval=FALSE, include=TRUE}

python MethByBin.py --featCoord genomeSequence-transcriptGFFCoord.tab --sample listSampleBGZfiles.txt --feature gene --out methylation_record_gene_flank_bin20_binfl100.txt --bin 20 --context CG -flanks --bin_fl 100

Rscript MethByBin_Plot.R --input methylation_record_gene_flank_bin20_binfl100.txt --label gene --nbin 20 --output Example2 --flank --nbin_flank 100 --group sample

```
![Example 2](https://github.com/pedro-mb/parsing-BSseqData/blob/master/tutorial/Example2.png)

<br/>

## Example3: Using a list of different files containing feature coordinates:

Generate a plot for the same tissue/sample (including replicates) and different groups of genes. For example, define groups of genes according to their expression level in that tissue (in 4 quartiles: Q4-Highly expressed, Q3-moderatly expressed, Q2-low expressed, Q1-very low expressed/inactive)

A) filter the *genomeSequence-transcriptGFFCoord.tab* to get coordinate files only for the genes of interest from each quartil; <br/>
B) create a text file (e.g. *quartil_transcript_coord.txt*) containing two columns (tab separared) with the path for these new files and the group ID (Q1, ...):
For example:

  path/to/coord/gene_group1.gff.tab &nbsp; Q1 <br/>
  path/to/coord/gene_group2.gff.tab &nbsp; Q2 <br/>
  path/to/coord/gene_group3.gff.tab &nbsp; Q3 <br/>
  path/to/coord/gene_group4.gff.tab &nbsp; Q4 <br/>
  ... <br/>

<br/> 

To run the *MethByBin.py*, use the previous file as input for argument --featCoord and add the option -listCord 

For the --sample argument use only a list of bgz files (biol. replicates) for one sample at a time, since the MethByBin_Plot.R is not prepared to deal with multiple factors (sample and gene groups)

To run the *MethByBin_Plot.R*, use --group coord



```{bash eval=FALSE, include=TRUE}
python MethByBin.py --featCoord quartil_transcript_coord.txt --sample listTissueABGZgiles.txt --feature transcript --context CG --out methylation_record_quartile_gene20_flank100.txt -flanks --bin 20 -listCoord


Rscript MethByBin_Plot.R --input methylation_record_quartile_gene20_flank100.txt --label gene --nbin 20 --output Example3 --flank --nbin_flank 100 --group coord
```

![Example 3](https://github.com/pedro-mb/parsing-BSseqData/blob/master/tutorial/Example3.png)
<br/>

MethByBin_Plot.R uses ggplot and the script can be edited to change general display (geom_line(), stat_summary()), ...) or to change line colours.
By default, plots show mean methylation frequency per bin (if different replicates are used for the same tissue/sample) in coloured lines, and methylation frequency per bin for each replicate in  grey lines (to check sample dispersion). To plot mean and confidence interval (in grey) add "-cl" option to the command


