# Missing-Value-Imputation-for-Mixed-Data-via-Gaussian-Copula

This repository provides codes to reproduce the experiment results in paper: Yuxuan Zhao, and Madeleine Udell. ["Missing Value Imputation for Mixed Data via Gaussian Copula."](https://dl.acm.org/doi/abs/10.1145/3394486.3403106?casa_token=y5buF87Ip1kAAAAA:LNesorhIMdx6ZuXxXz8UdTuwPuJrG2q-CWphz1mxnI-s6KA6FrnCD6KghK9t_9UXkKE4-z_BzhsShA). The imputation algorithm introduced in the paper is available in both [R](https://github.com/udellgroup/mixedgcImp) and [Python](https://github.com/udellgroup/online_mixed_gc_imp).

## Software loading

In our experiments, we use our copula EM algorithm in the [R pacakge](https://github.com/udellgroup/mixedgcImp). For other algorithms used in our paper, we adopt the following implementation: [sbgcop](https://cran.r-project.org/web/packages/sbgcop/index.html), [missForest](https://cran.r-project.org/web/packages/missForest/index.html), imputeFAMD in R package [missMDA](https://cran.r-project.org/web/packages/missMDA/index.html), [xPCA](https://gitlab.com/xpca/xpcar), [softImpute](https://cran.r-project.org/web/packages/softImpute/index.html) and [GLRM](https://github.com/madeleineudell/LowRankModels.jl). All codes are public available. To replicate our experimental results, make sure all invovled algorithms are installed.

## Datasets

We include a copy for each used real world dataset. They are also all public available: [GSS](https://gssdataexplorer.norc.org/), [MovieLens 1M](https://grouplens.org/datasets/movielens/1m/), [CAL500exp](http://slam.iis.sinica.edu.tw/demo/CAL500exp/), [GBSG](https://cran.r-project.org/web/packages/mfp/), [TIPS](http://ggobi.org/book/), and ESL and LEV at https://waikato.github.io/weka-wiki/datasets/ (item named "A gzip'ed tar containing ordinal, real-world datasets donated by Professor Arie Ben David (datasets-arie_ben_david.tar.gz, 11,348 Bytes)").

## Organization

Our provided codes are divided by figures/tables. While replicating the exact results may be time demanding, users can use a smaller number of repetitions. See each file for details. Currently, codes for replicating Figure 4 and Table 4 are provided. 

More to be updated...
