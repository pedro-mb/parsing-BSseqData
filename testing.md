<br/>

General workflow for whole genome methylation analysis:
<br/>

**- First step, get the coordinates for the feature of interest, depending on the plots to generate: **

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

![a test](https://github.com/pedro-mb/parsing-BSseqData/blob/master/Example3.png)
To get intergenic regions I used Bedtools. 'bedtools complement' requires a gff file only with gene rows and a tab separated file (*scaffoldSizes.tsv*) containing all scaffoldIDs and corresponding size in bp:

For example:


  scaffold1  &nbsp; 22220034 <br/>
  scaffold2  &nbsp; 30300303 <br/>
  ... <br/>
<br/>
