Basic GO Analysis with mogo
================
Dan MacLean
23 June, 2021

## mogo

`mogo` is a package for generating GO analysis using *M*.oryzae genes
only. You can provide the main function with the list of genes of
interest and it does (most of) the rest\!

## Installation

To install you’ll need the `devtools` R package. From the R console in
RStudio type

``` r
install.packages("devtools")
```

Once `devtools` is installed you can use that to install `mogo`

``` r
devtools::install_github("TeamMacLean/mogo")
```

## Preparation

The first step is to load our gene expression data file and get a simple
list of genes for GO analysis

### Data loading

First read in your gene expression file. You can do that with
`read_csv()`. Minimally it should contain the gene ID, the log fold
change and the \(p\)-value.

``` r
library(readr)
library(dplyr)
library(here)

gene_expression <- read_csv(here("inst", "extdata","sample_gene_expression.csv"))
gene_expression
#> # A tibble: 132 x 3
#>    gene_id   log2fc  p.adj
#>    <chr>      <dbl>  <dbl>
#>  1 MGG_16981  3.53  0.291 
#>  2 MGG_04557  1.51  0.217 
#>  3 MGG_03894 -7.03  0.255 
#>  4 MGG_09290 -7.19  0.377 
#>  5 MGG_01448  2.11  0.174 
#>  6 MGG_09589 -3.04  0.183 
#>  7 MGG_03476 -0.341 0.199 
#>  8 MGG_16442  3.60  0.0434
#>  9 MGG_01048 -3.36  0.309 
#> 10 MGG_14883 -4.93  0.157 
#> # … with 122 more rows
```

### Filter as required

We can now filter the genes to select only the ones with e.g
\(p <= 0.05\), using `filter()`

``` r
filtered_gene_expression <- filter(gene_expression, p.adj <= 0.05)
filtered_gene_expression
#> # A tibble: 16 x 3
#>    gene_id   log2fc    p.adj
#>    <chr>      <dbl>    <dbl>
#>  1 MGG_16442  3.60  0.0434  
#>  2 MGG_14022  3.51  0.00180 
#>  3 MGG_00427  0.462 0.0174  
#>  4 MGG_10367 -1.30  0.0215  
#>  5 MGG_00015  0.539 0.0398  
#>  6 MGG_01295 -4.53  0.0360  
#>  7 MGG_05721 -8.60  0.000341
#>  8 MGG_13781  7.58  0.0353  
#>  9 MGG_14886 -0.674 0.00469 
#> 10 MGG_14863 -1.63  0.0322  
#> 11 MGG_01492 -6.62  0.0437  
#> 12 MGG_04901 -2.27  0.0440  
#> 13 MGG_04562 -1.67  0.0439  
#> 14 MGG_00559  0.574 0.0324  
#> 15 MGG_14897  5.88  0.00736 
#> 16 MGG_03526 -0.932 0.0440
```

### Extract the gene id column

We can now extract the `gene_id` column using the `$` syntax

``` r
gene_ids <- filtered_gene_expression$gene_id
gene_ids
#>  [1] "MGG_16442" "MGG_14022" "MGG_00427" "MGG_10367" "MGG_00015" "MGG_01295"
#>  [7] "MGG_05721" "MGG_13781" "MGG_14886" "MGG_14863" "MGG_01492" "MGG_04901"
#> [13] "MGG_04562" "MGG_00559" "MGG_14897" "MGG_03526"
```

## GO enrichment

The GO enrichment is done in the `mogo` package. Load that and use the
`do_enrich()` function, passing it the vector of `gene_ids` to calculate
the enrichment.

``` r
library(mogo)
enrich <- do_enrich(gene_ids)
enrich
#> #
#> # over-representation test
#> #
#> #...@organism     UNKNOWN 
#> #...@ontology     UNKNOWN 
#> #...@gene     chr [1:16] "MGG_16442" "MGG_14022" "MGG_00427" "MGG_10367" "MGG_00015" ...
#> #...pvalues adjusted by 'BH' with cutoff <0.05 
#> #...5 enriched terms found
#> 'data.frame':    5 obs. of  9 variables:
#>  $ ID         : chr  "GO:0032259" "GO:0008168" "GO:0008171" "GO:0008033" ...
#>  $ Description: chr  "The process in which a methyl group is covalently attached to a molecule." "Catalysis of the transfer of a methyl group to an acceptor molecule." "Catalysis of the transfer of a methyl group to the oxygen atom of an acceptor molecule." "The process in which a pre-tRNA molecule is converted to a mature tRNA, ready for addition of an aminoacyl group." ...
#>  $ GeneRatio  : chr  "16/16" "16/16" "3/16" "2/16" ...
#>  $ BgRatio    : chr  "133/10093" "136/10093" "18/10093" "55/10093" ...
#>  $ pvalue     : num  3.27e-31 4.77e-31 2.63e-06 3.33e-03 6.27e-03
#>  $ p.adjust   : num  6.21e-30 6.21e-30 2.28e-05 2.17e-02 3.26e-02
#>  $ qvalue     : num  3.27e-30 3.27e-30 1.20e-05 1.14e-02 1.72e-02
#>  $ geneID     : chr  "MGG_16442/MGG_14022/MGG_00427/MGG_10367/MGG_00015/MGG_01295/MGG_05721/MGG_13781/MGG_14886/MGG_14863/MGG_01492/M"| __truncated__ "MGG_16442/MGG_14022/MGG_00427/MGG_10367/MGG_00015/MGG_01295/MGG_05721/MGG_13781/MGG_14886/MGG_14863/MGG_01492/M"| __truncated__ "MGG_14022/MGG_00427/MGG_00015" "MGG_04562/MGG_00559" ...
#>  $ Count      : int  16 16 3 2 2
#> #...Citation
#>   Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
#>   clusterProfiler: an R package for comparing biological themes among
#>   gene clusters. OMICS: A Journal of Integrative Biology
#>   2012, 16(5):284-287
```

As you can see the `enrich` object has a lot of information in it. The
result table can be extracted using `as.data.frame()` to convert the
`enrich` to a data.frame (note `glimpse` is a helpful function for
printing out big dataframes).

``` r
result_table <- as.data.frame(enrich)
glimpse(result_table)
#> Rows: 5
#> Columns: 9
#> $ ID          <chr> "GO:0032259", "GO:0008168", "GO:0008171", "GO:0008033", "G…
#> $ Description <chr> "The process in which a methyl group is covalently attache…
#> $ GeneRatio   <chr> "16/16", "16/16", "3/16", "2/16", "2/16"
#> $ BgRatio     <chr> "133/10093", "136/10093", "18/10093", "55/10093", "76/1009…
#> $ pvalue      <dbl> 3.269250e-31, 4.773264e-31, 2.629052e-06, 3.331689e-03, 6.…
#> $ p.adjust    <dbl> 6.205243e-30, 6.205243e-30, 2.278511e-05, 2.165598e-02, 3.…
#> $ qvalue      <dbl> 3.265917e-30, 3.265917e-30, 1.199217e-05, 1.139788e-02, 1.…
#> $ geneID      <chr> "MGG_16442/MGG_14022/MGG_00427/MGG_10367/MGG_00015/MGG_012…
#> $ Count       <int> 16, 16, 3, 2, 2
```

You can save the result to an excel-compatible csv file with,
`write_csv()`

``` r
write_csv(result_table, "my_GO_results.csv")
```

## Plotting

The `enrich` object is from the package `ClusterProfiler` and can be
used directly in most of it’s plot types. See them at the
`ClusterProfiler` page
<http://yulab-smu.top/clusterProfiler-book/chapter12.html#bar-plot>

### Barplot

Works ok, but the text can make it a bit unwieldy

``` r
library(clusterProfiler)
#> Warning: package 'clusterProfiler' was built under R version 4.0.3
#> clusterProfiler v3.18.1  For help: https://guangchuangyu.github.io/software/clusterProfiler
#> 
#> If you use clusterProfiler in published research, please cite:
#> Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
#> 
#> Attaching package: 'clusterProfiler'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(enrichplot)
#> Warning: package 'enrichplot' was built under R version 4.0.3
barplot(enrich, showCategory=5)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Dotplot

Similar

``` r
dotplot(enrich, showCategory=5)
#> wrong orderBy parameter; set to default `orderBy = "x"`
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### UpSetPlot

These are like a really sophisticated Venn/Euler diagram.
<https://jku-vds-lab.at/tools/upset/>

``` r
upsetplot(enrich)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### GoPlot

Plots from the package `GoPlot` such as bubble plots can be made, but
need you to convert the `enrich` object to a DAVID compatible object,
and make an accessory data.frame of the expression information. Use
`enricher_to_david()` for the first part, and the relevant columns of
the gene expression data you created at the beginning for the second
part.

``` r
david <- enricher_to_david(enrich)
#> 
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   GeneID = col_character(),
#>   `GO term definition` = col_character(),
#>   `GO term accession` = col_character(),
#>   `GO domain` = col_character()
#> )
expr_info <- data.frame(
  ID = filtered_gene_expression$gene_id,
  logFC = filtered_gene_expression$log2fc
)
```

Then you can convert that to the data format needed for `GoPlot`

``` r
library(GOplot)
#> Loading required package: ggplot2
#> Loading required package: ggdendro
#> Loading required package: gridExtra
#> 
#> Attaching package: 'gridExtra'
#> The following object is masked from 'package:dplyr':
#> 
#>     combine
#> Loading required package: RColorBrewer
circ <- circle_dat(david, expr_info)
```

### Go Bar

A different sort of GO barchart

``` r
GOBar(subset(circ, category == 'BP'))
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### Go Bubble

A plot with bubbles

``` r
GOBubble(circ, labels=3)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->