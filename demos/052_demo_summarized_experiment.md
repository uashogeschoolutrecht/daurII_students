---
title: "daur2_lesson5_b"
author: "Marc A.T. Teunis, Alyanne de Haan"
date: '2020-05-19'
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---









## Introduction: SummarizedExperiment class, rows and columns

Based on the SummarizedExperiment class vignette.

### A little reminder of the ExpressionSet
In the previous lesson we encountered the Expressionset objects in GEO. Expression sets are a datacontainer for high throughput array data from expression experiments (in a matrix called _exprs_) and the related metadata. As you will have realised, you can access:

  + (always) expression data  with _exprs()_ ( one row per feature (gene or probe) on the chip and one column per sample. 
  
  + phenotypic data (information about the samples) with _pdata()_ (one row per sample, and one column per covariate)
  
  + and feature data (information about genes or probe sets) with _fdata()_ (as with the expression data, one row per feature, one column per variable)

Each data point in the exprs matrix reflects how that particular gene or feature was expressed in that particular sample. 

### But there is more
But the Expressionset is not the only possible datacontainer format in bioconductor.

Sometimes you want to share or analyse data from an experiment with multiple types of observations, such as DNA mutations and abundance of RNA and proteins. For sets of assays with the same information across all rows (e.g., genes or genomic ranges) the datacontainer is called `SummarizedExperiment`. SummarizedExperiment objects can therefore contain multiple matrices with assay data. 

But it can also be used to store data from 1 assay. SummarizedExperiment are based on Expressionsets, but can contain GenomicRanges. `SummarizedExperiment` is therefore often used to store data from sequencing-based experiments (with GenomicRanges in the rows) while `Expressionsets` are more often used for array-based experiments (with features in the rows). In both datacontainers, the columns represent samples.


### Structure 

The structure of an SummarizedExperiment object looks like this:

<div class="figure" style="text-align: center">
<img src="D:/r_projects/daurII_students/images/summarizedexperiment.png" alt="..." width="50%" />
<p class="caption">...</p>
</div>

`assay()`, `assays()`: A matrix-like or list of matrix-like objects of identical dimension

  + matrix-like: implements dim(), dimnames(), and 2-dimensional [, [<- methods.
  
  + rows: genes, genomic coordinates, etc.
  
  + columns: samples, cells, etc.
  
`colData()`: Annotations on each column, as a DataFrame.

  + E.g., description of each sample
  + So this is phenotype (sample) information, like pData for ExpressionSets. 
  
`rowData()` and / or `rowRanges(`): Annotations on each row.

  + E.g., rowRanges(): coordinates of gene / exons in transcripts / etc.
  
  + E.g., rowData():P-values and log-fold change of each gene after differential expresssion analysis.
  
`metadata()`: List of unstructured metadata describing the overall content of the object.

You can subset SummarizedExperiment objects like you would subset R: matrices. So suppose your SummarizedExperiment object is called _se_, you can make a subset containing only the first 10 rows with _se[1:10,]_.


As most tutorials, we will use the example RNA sequencing dataset `Airway` to illustrate the use of SummarizedExperiment (paper PMID: 24926665, GEO: GSE52778, remember you learned how to use this information with rentrez and geoquery in a previous lesson.)

__Example (1)__

```r
BiocManager::install("airway")
```

```r
library(airway)
data(airway)
```

This looks a lot like an Expressionset. You can see that there are 8 samples and 64102 features (genes). But what is this about? Let's check PubMed as we have done before, and use XML package (install.packages("XML") first!) to extract the abstract. Read at least the first 4 sentences of the abstract and make sure you understand what this experiment was investigating. 

__Example (2)__

```r
library(rentrez)
library(XML)

fetch.pubmed <- entrez_fetch(db = "pubmed", id = 24926665,
                      rettype = "xml", parsed = T)
abstract = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x) xmlValue(xmlChildren(x)$Abstract))
```

We can check the phenotypic data:
__Example (3)__

```r
colData(airway)
```

```
## DataFrame with 8 rows and 9 columns
##            SampleName     cell      dex    albut        Run avgLength
##              <factor> <factor> <factor> <factor>   <factor> <integer>
## SRR1039508 GSM1275862   N61311    untrt    untrt SRR1039508       126
## SRR1039509 GSM1275863   N61311      trt    untrt SRR1039509       126
## SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126
## SRR1039513 GSM1275867  N052611      trt    untrt SRR1039513        87
## SRR1039516 GSM1275870  N080611    untrt    untrt SRR1039516       120
## SRR1039517 GSM1275871  N080611      trt    untrt SRR1039517       126
## SRR1039520 GSM1275874  N061011    untrt    untrt SRR1039520       101
## SRR1039521 GSM1275875  N061011      trt    untrt SRR1039521        98
##            Experiment    Sample    BioSample
##              <factor>  <factor>     <factor>
## SRR1039508  SRX384345 SRS508568 SAMN02422669
## SRR1039509  SRX384346 SRS508567 SAMN02422675
## SRR1039512  SRX384349 SRS508571 SAMN02422678
## SRR1039513  SRX384350 SRS508572 SAMN02422670
## SRR1039516  SRX384353 SRS508575 SAMN02422682
## SRR1039517  SRX384354 SRS508576 SAMN02422673
## SRR1039520  SRX384357 SRS508579 SAMN02422683
## SRR1039521  SRX384358 SRS508580 SAMN02422677
```

And make subsets based on colData as well. For instance, suppose you only want to analyse data from the treated condition:
__Example (4)__

```r
airway[, airway$dex == "trt"]
```

```
## class: RangedSummarizedExperiment 
## dim: 64102 4 
## metadata(1): ''
## assays(1): counts
## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
## rowData names(0):
## colnames(4): SRR1039509 SRR1039513 SRR1039517 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

Colnames() is similar to sampleNames from ExpressionSet and rownames() are like featureNames.

__Example (4)__

```r
head(rownames(airway))
```

```
## [1] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
## [5] "ENSG00000000460" "ENSG00000000938"
```

And the actual measurements data are found with assay() and assays(). This particular dataset has only one assay, and it is called "counts" and it contains RNA sequencing count data:
__Example (5)__

```r
assayNames(airway)
```

```
## [1] "counts"
```


Let's print the first few rows:
__Example (6)__

```r
head(assays(airway)$counts,10)
```

```
##                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
## ENSG00000000003        679        448        873        408       1138
## ENSG00000000005          0          0          0          0          0
## ENSG00000000419        467        515        621        365        587
## ENSG00000000457        260        211        263        164        245
## ENSG00000000460         60         55         40         35         78
## ENSG00000000938          0          0          2          0          1
## ENSG00000000971       3251       3679       6177       4252       6721
## ENSG00000001036       1433       1062       1733        881       1424
## ENSG00000001084        519        380        595        493        820
## ENSG00000001167        394        236        464        175        658
##                 SRR1039517 SRR1039520 SRR1039521
## ENSG00000000003       1047        770        572
## ENSG00000000005          0          0          0
## ENSG00000000419        799        417        508
## ENSG00000000457        331        233        229
## ENSG00000000460         63         76         60
## ENSG00000000938          0          0          0
## ENSG00000000971      11027       5176       7995
## ENSG00000001036       1439       1359       1109
## ENSG00000001084        714        696        704
## ENSG00000001167        584        360        269
```

### Granges

So far so good. This could just have been an ExpressionSet. But SummarizedExperiments have the option to store Genomic Ranges lists for each feature (i.e. row). So in other words: in the assay matrix in the airway dataset, each row contains data from a gene. For each gene, we also get a list (GRanges list) containing the exons of the gene. You can access this with _rowRanges()_. 


You can see where the first gene in the assay matrix is located: on chromosome x. Also, it has  17 exons:

__Example (7)__

```r
rowRanges(airway)[1,]
```

```
## GRangesList object of length 1:
## $ENSG00000000003
## GRanges object with 17 ranges and 2 metadata columns:
##        seqnames            ranges strand |   exon_id       exon_name
##           <Rle>         <IRanges>  <Rle> | <integer>     <character>
##    [1]        X 99883667-99884983      - |    667145 ENSE00001459322
##    [2]        X 99885756-99885863      - |    667146 ENSE00000868868
##    [3]        X 99887482-99887565      - |    667147 ENSE00000401072
##    [4]        X 99887538-99887565      - |    667148 ENSE00001849132
##    [5]        X 99888402-99888536      - |    667149 ENSE00003554016
##    ...      ...               ...    ... .       ...             ...
##   [13]        X 99890555-99890743      - |    667156 ENSE00003512331
##   [14]        X 99891188-99891686      - |    667158 ENSE00001886883
##   [15]        X 99891605-99891803      - |    667159 ENSE00001855382
##   [16]        X 99891790-99892101      - |    667160 ENSE00001863395
##   [17]        X 99894942-99894988      - |    667161 ENSE00001828996
##   -------
##   seqinfo: 722 sequences (1 circular) from an unspecified genome
```


But you can also do it the other way around, and look for the data containing information about a certain region of interest. In this example the interval 100,000 to 110,000 of chromosome 1:

__Example (8)__

```r
# define region of interest:
roi <- GRanges(seqnames="1", ranges=100000:1100000)
# new SummarizedExperiment object containing 74 rows:
roi_se <- subsetByOverlaps(airway, roi) 
# print the first 6 rows:
head(assays(roi_se)$counts,6) 
```

```
##                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
## ENSG00000131591        111        114        137         64        151
## ENSG00000177757          2          5          8          1          7
## ENSG00000185097          0          0          0          0          0
## ENSG00000187583         19          9         13          6         17
## ENSG00000187608        319        234        214         90        199
## ENSG00000187634        136         68         98          6        101
##                 SRR1039517 SRR1039520 SRR1039521
## ENSG00000131591        199        130        122
## ENSG00000177757          6          3          4
## ENSG00000185097          1          1          0
## ENSG00000187583         19          9         14
## ENSG00000187608        228        143        174
## ENSG00000187634         44         54         18
```

The GRanges object can contain  information about genes on the standard chromosomes, as well as genes or haplotypes of genes on haplotype chromosomes and mitochondrial. So that is why the seqnames column of a GRanges object may contain something other than numbers between 1 and 22, x and y.

vragen ongeveer hier



