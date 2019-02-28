# Sincell
### R/Bioconductor package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data

[Bioconductor version: Release (3.8)](http://bioconductor.org/packages/release/bioc/html/sincell.html)

Cell differentiation processes are achieved through a continuum of hierarchical intermediate cell-states that might be captured by single-cell RNA seq. Existing computational approaches for the assessment of cell-state hierarchies from single-cell data might be formalized under a general workflow composed of i) a metric to assess cell-to-cell similarities (combined or not with a dimensionality reduction step), and ii) a graph-building algorithm (optionally making use of a cells-clustering step). Sincell R package implements a methodological toolbox allowing flexible workflows under such framework. Furthermore, Sincell contributes new algorithms to provide cell-state hierarchies with statistical support while accounting for stochastic factors in single-cell RNA seq. Graphical representations and functional association tests are provided to interpret hierarchies.

Author: Miguel Julia <migueljuliamolina at gmail.com>, Amalio Telenti <atelenti at jcvi.org>, Antonio Rausell <antonio.rausell at isb-sib.ch>

Maintainer: Miguel Julia <migueljuliamolina at gmail.com>, Antonio Rausell <antonio.rausell at isb-sib.ch>

Citation (from within R, enter `citation("sincell")`):

```
Juliá M, Telenti A, Rausell A (2014). “Sincell: R package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data.” bioRxiv. http://dx.doi.org/.
```

## Installation
To install this package, start R (version "3.5") and enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sincell", version = "3.8")
```

For older versions of R, please refer to the appropriate [Bioconductor release](http://bioconductor.org/about/release-announcements/). 

## Documentation
To view documentation for the version of this package installed in your system, start R and enter:

```
browseVignettes("sincell")
```

| Document      | Script        | Title              |
| :-------------: |:-------------:| :------------------|
| [PDF](http://bioconductor.org/packages/release/bioc/vignettes/sincell/inst/doc/sincell-vignette.pdf) | [R Script](http://bioconductor.org/packages/release/bioc/vignettes/sincell/inst/doc/sincell-vignette.R) | Sincell: Analysis of cell state hierarchies from single-cell RNA-seq |
| [PDF](http://bioconductor.org/packages/release/bioc/manuals/sincell/man/sincell.pdf) | - | Reference Manual |
| [Text](http://bioconductor.org/packages/release/bioc/news/sincell/NEWS) | - | NEWS |

## [Download stats](http://bioconductor.org/packages/stats/bioc/sincell/)



