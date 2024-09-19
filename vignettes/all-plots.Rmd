---
title: "Structure of the all_plots list"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Structure of the all_plots list}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r opts, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

It's only necessary to understand `all_plots` if you intend to modify the automatically generated
plots. 

Before the report is rendered, the plot objects are stored in a list. This allows the user to adjust
any of the plots and re-run `arnaq.create.report()` or `export_all_svg()` to include the
updated version in the output. As `ggplot2` plots, they can be latered using any of the usual 
`ggplot2` grammar. Alternatively, or in combination with such adjustments, the functions that create
these plots can be used with different parameters. Plots can be turned by altering the
`arnaq.plot.flags` object and rerunning one of the above functions.

The `_ps` suffix is a 
convention to indicate that the list entry contains a list of unnamed subplots, usually for 
clarity where there are more samples than can comfortably fir on one plot.

The full structure of the object is:

```{text}
all_plots
├── $barplots_ps
│   ├── [[1]]
│   ├── [[...]]
│   └── [[i]]
├── $barplots2_ps
│   ├── [[1]]
│   ├── [[...]]
│   └── [[i]]
├── $violin1
├── $violin2
├── $normplot
├── $hcluster_ps
│   ├── [[1]]
│   ├── [[2]]
│   ├── [[3]]
│   └── [[4]]
├── $pca_ids_ps
│   ├── [[1]]
│   ├── [[...]]
│   └── [[i]]
├── $pca_group_pss
│   ├── [[1]]
│   │   ├── [[1]]
│   │   ├── [[...]]
│   │   └── [[i]]
│   ├── [[...]]
│   │   ├── [[1]]
│   │   ├── [[...]]
│   │   └── [[i]]
│   └── [[j]]
│       ├── [[1]]
│       ├── [[...]]
│       └── [[i]]
├── $complexity_plot
├── $complexity_plot_2
├── $detected_plot
├── $detected_plot_2
├── $ercc_mix1_plot
├── $ercc_mix2_plot
├── $ercc_fc_ps
│   ├── [[1]]
│   ├── [[...]]
│   └── [[i]]
├── $biotype_ps
│   ├── [[1]]
│   ├── [[...]]
│   └── [[i]]
├── $biotype_pie_plot
└── $read_scatterplots_ps
    ├── [[1]]
    ├── [[...]]
    └── [[i]]
```

# Plot descriptions

Not all plots will be present in every run. Inclusion is based on factors such as what data is 
available, how many samples are in the set, and the parameters supplied to `arnaq()`.

## `barplots_ps` and `barplots2_ps`

*Included when read summary data (as from featureCounts) is given in `resources.yml`.*

These are lists of bar plots showing the read assignment summaries. In `barplots_ps`, they are 
displayed as counts; in `barplots2_ps` they are displayed as proportions.

How the samples are split up across multiple plots is controlled by the `max.samples.per.page` 
parameter to `arnaq()`. Note that even when all the samples fit in a single plot, these entries
are still lists of length 1.

## `violin1` and `violin2`

*Included when additional sample metrics (as from Picard) are given in `resources.yml`.*

These show additional metrics in the form of violin plots, split by group. In this case, only one 
group is used, rather than every group as in other plots. The group used is the first named in the 
`treat.groups` parameter.

The metrics are split across two plots, for clarity, based on scale. Count-based metrics are in 
`violin1`. Those based on proportions are in `violin2`.

## `normplot`

*Included when variance-stabilising normalisation is performed.*

This is the standard `meanSdPlot` in the `DESeq2` package.

## `hcluster_ps`

*Always included.*

This is a list of four plots showing how the samples cluster. Two types of plot are show; a 
sample-to-sample similarity heatmap and a hiercarchical clustering dendrogram. Two similarity 
metrics are considered; Euclidean distance based on the normalised read counts and Spearman's rank 
correlation of the same.

## `pca_ids_ps`

*Always included.*

These show the PCA plot with samples labelled by their names. Each requested component from the 2nd
onwards is plotted against component 1. The number of components considered is controlled by 
`pca.depth` and the default is to consider every component that accounts for 10% of the total 
variability or more.

Even if there is only a single plot, this is still a list of length 1.

## `pca_group_pss`

*Always included.*

This is a list of lists of PCA plots.

Each top level list represents one of the groups specified in `pca.groups`. In the case where this 
is only a single group, the list will be of length 1.

The sublists have the same format as `pca_ids_ps`, above, with one entry per depth of PCA. The 
samples are coloured by group membership.

## `complexity_plot` and `detected_plot`

*Included when the number of samples is <= 40.*

The complexity plot shows what proportion of the total reads are assigned to the top N genes, as 
N varies. The detected genes plot show how many genes are detected at different thresholds. Each 
line is one sample.

In these versions, every sample is shown.

## `complexty_plot_2` and `detected_plot_2`

*Included when the number of samples is > 20.*

As above, but in this case only the worst 8 samples for each plot are shown in colour. The remaining
samples in the set are coloured grey, behind.

## `ercc_mix1_plot` and `ercc_mix2_plot`

*Included when `ERCC` is `TRUE` and the there is an `ERCC` column in the sample metadata.*

These two plots show the observed cpm for ERCC spike-in mixes against the expected.

## `ercc_fc_ps`

*Included when `ERCC` is `TRUE` and at least one pair of samples has been given in `ERCC.pairs`.*

A list of plots showing observed ERCC spike-in fold-changes between two samples with different 
mixes against the known fold-changes.

## `biotype_ps`

*Included always.*

A list of barplots showing the counts of biotypes for the observed reads, as given in the `.gtf` 
support file and mapped to a smaller nub mer of categories by the biotype conversion table. Where 
the number of samples exceeds `max.samples.per.page`, this is split across multiple plots for 
clarity, as with `barplots_ps`. Similarly, if all the samples fit on one plot, the list will be of
length 1.

## `biotype_pie`

*Included always.*

A pie chart showing the biotypes of observed genes, not with respect to read abundance.

## `read_scatterplots_ps`

*Included when the length of `scatter.pairs` is at least 1.*

A list of plots showing sample-to-sample comparisons, one point per gene.