Capstone Exam - ABDS, Module DAUR2/BIOC
================
\<Your Name & student\_id\>
2020-05-20

# Directions and tips

**PLEASE PROVIDE THE FOLLOWING IN THE YAML HEADER ABOVE BEFORE
SUBMITTING THIS EXAM**

1.  The GEO accesion number and GSM sample numbers as a subtitle
2.  Your Full Name, as registered at the Hogeschool Utrecht in the
    author field
3.  Your HU student\_id in the author field

**Write an Rmd file, containing code chunks and answer all the questions
(or steps) required to complete the exam below. Also include a rendered
version of this Rmd in your submission. Submitting both files can be
doen before the deadline 21 June, 2019 in CANVAS under assignment
“BIOC\_EXAM”**

**The total number of points that can be gained in this exam is 24
(excluding the bonus question. The grade will be calculated from:**

\((\frac{your score + bonus.score}{26}) * 10\)

## TIPS

You probably will need at least to use the following packages and
functions:

``` 
    tidyverse
     readr
     forcats
     ggplot2
     dplyr
     tibble
    GEOquery
      getGEO()
      getGEOSuppFiles()
    biomaRt
      getBM()
      listAttributes()
    DESeq2
      <many functions>
    heatmap.2 / pheatmap
    here
      here()
    MutationalPatterns
      
```

*Mind that these are suggestions and that this list is not complete, you
will need to figure out which packages you will need additionally.*

You are allowed to Google around for code and reuse code available on
the internet. Please provide adequate reference to code if you take it
from a published source such as blogs, kernels, publications, vignettes
etc. Code form fora such as StackOverflow usually do not need reference.

# Packages

``` r
library(tidyverse)
## add additional packages here:--->
```

# Setting the project root

We need an object for the root of the project or use the `{here}`
package instead.

``` r
if (!require("rprojroot")) install.packages("rprojroot")
library("rprojroot")
root <- find_root_file(criterion = is_rstudio_project)
root

## or alternatively use the `{here}` package
```

# Introduction to the exam

In this Exam we request you to perform a full analyis of a given
Next-Generation Differential Gene Expression Experiment Dataset. The
dataset will be provided to you as a Gene Expression Omnibus (GEO)
accession number.

You are required to write an anlaysis R-Markdown (Rmd) script that
addresses all of the questions (steps) below. An additional bonus
question on the `{MutationalPatterns}` package is also included in the
exam.

The Exam Rmd script and it’s rendered html version (knit) need to be
uploaded before the deadline:

**22 june 2019, 12:00 a.m. (noon) hrs CEST**

Submissions after this deadline will be considered invalid and will be
graded as ‘NA’ in OSIRIS, so make sure you submit on time.

# Steps in the required analysis

To pass this Exam you need to complete all parts of the ANALYSIS STEPS
below:

**In short:** *You need to address the following in your analysis:*

  - Download the data (including raw counts) of the associated GEO
    accession id
  - Choose the appropriate design formula for your dataset: For this you
    need to analyze the structure of the dataset and the information
    about the performed experiment. You need to discover which
    experimental groups were used. Try answering the question: “What
    makes sense to compare in this experiment?” and choose your design
    (formula) accordingly
  - Generate all exploratory graphs and visualizations such as
    dispersion or MA plots
  - Generate a result table that shows the results of your analysis
  - Generate four annotated (gene\_names/symbols) heatmaps (two showing
    the top 100 significant (up and down in a seperate heatmap) and the
    top 20 significant (up and down in a seperate heatmap)
  - Find the homologous gene sequences for the top 5 genes for Orang
    Utan
  - Perform a multiple sequence alignment where you compare the sequence
    of the top 5 genes to the 5 homologous genes from the Orang Utan
  - Find the protein sequences for the homologous sequences detected

For a succesfull submission of this exam execute the following steps in
your Rmd script:

# *STEP 1: Download the data* (4 points)

Download the provided dataset (raw counts - matrix) and any additional
data you need to complete the steps below. You will also be required to
access other data sources for example via `{GEOquery}`, `{Organism.db}`
or `{Annotation.dbi}` packages. The data required has different types:
RNA sequencing raw counts, DNA sequences, Genomic information,
annotations etc.

# *STEP2: Generate a `SummarizedExperiment` object* (4 points)

The `SummarizedExperiment` class R object must at least include the
following items:

    - `AssayData`: A raw counts matrix, annotated with feature ids (rownames)
    - `PhenoData`/`ColData`: An annotatedDataFrame with column information, corresponding to the AssayData matrix  
    - `ExperimentData`/Metadata: An object (MIAME) with metainformation about the Experiment
    - `FeatureData`: A dataframe with linkage between the feature ids (usually ensembl_ids or entrez_ids) and additional gene annotions such as symbol (gene_name).

## Construct DESeqDataSet from counts and coldata

If you have the counts and coldata in two seperate files or objects you
can construct a DESeqDataset by using the function
`DESeq2::DESeqDataSetFromMatrix()`

You need to make sure that the colnames of the counts matrix are equal
to the rownames of the coldata

Example

    DESeq2::DESeqDataSetFromMatrix(countData = ...,
    colData = ..., design = ~ cell + dex)

# *STEP3: Differential Gene Expression Profile (DE)* (8 points)

Using the workflow explained during the course, run a complete
Differential Gene Expression analysis on the supplied dataset.

See: lesson 6 “rna\_seq workflow” theory, the demo and the bioconductor
workflow:
<https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html>

**In short: You need to address the following in your DE analysis:**

  - Choose the appropriate design formula for your dataset: For this you
    need to analyze the structure of the dataset and the information
    about the performed experiment. You need to discover which
    experimental groups were used. Try answering the question: “What
    makes sense to compare in this experiment?” and choose your design
    (formula) accordingly.
  - Generate all exploratory graphs and visualizations such as
    dispersion or MA/pca plots.
  - Generate a result table that shows the results of your analysis.

# *STEP 4: Heatmaps* (4 points)

  - Generate four annotated (gene\_names/symbols) heatmaps (two showing
    the top 100 significant (up and down in a seperate heatmap) and the
    top 20 significant (up and down in a seperate heatmap). Use
    e.g. `{biomaRt}` or `{biomartr}` or `{Annotation.dbi}` for
    annotating the result-table. Think about one-to-many or many-to-one
    mappings. Remove rows if neccessary.

**REMEMBER** The top genes can be found by ranking the genes on the
basis of the value for `p.adjusted` in the results rable from you
Differential Expression analysis.

# *STEP 5: Sequence alignment* (4 points)

  - Find the `coding sequence` for the top 5 genes, using the
    `{biomaRt}` or `{biomartr}` package.
  - Find `homologous` sequences for your gene for the `Orang utan` using
    `{biomaRt}`
  - Perform alignments for the top 5 genes (DNA sequences) and
    homologous genes to Orang utan. If no genes were found for Orang
    utan, choose a species of your own choice. (for example chinchilla,
    lama, horse, cow, sheep, chimpanzee)

# *STEP 6: Conclusions* (2 POINTS)

  - Write a short conclusive paragraph on your differential expression
    analysis. What steps do you think are needed next?

# *7; Bonus QUESTION* (4 points)

During the lesson on functions and iterations you learned how to write a
function and perform loop in R.

Write a function that takes a single GEO accession number and some GSM
sample ids as input and creates the full analysis report for that
function. Test your function with the GEO dataset you analyzed for the
exam. Think about any additional arguments to your function you might
need and implement them in your function. Include the function
definition and proof that it works in your exam submission.
