# Incorporating Animal Movement Into Distance Sampling

# Author Contributions Checklist Form

## Data

### Abstract

Line transect survey of Pan-tropical Spotted Dolphins

### Availability

Data are provided by the Southwest Fisheries Science Center (SWFSC), National Oceanic and Atmospheric Administration (NOAA).

### Description

The data is provided in the supplementary materials.

## Code

### Abstract

The code provided reproduces the methods applied in the paper including scripts to apply the methods to simulated data and to reproduce the simulation study results.

### Description

The code is provided in two ZIP files. One ZIP file contains an R package (moveds) to be installed. The second ZIP file uses this R package to show how the methods are applied and to reproduce the simulation studies. MIT License applies to all code.

### Optional Information

Code should run on any desktop computer with modest resources. The Description file within the R package details what dependencies must be installed to run the code.

## Instructions for Use

### Reproducibility 

Please refer to the README provided in the supplementary material.

### Replication 

The package included contains vignettes that explain the implementation of the methods and show how to apply it to simulated data. Below is the version of R, version of R packages, and other system-specific information that was used when producing the results in the paper:

```
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.3 LTS

Matrix products: default
BLAS:   /usr/local/lib/R/lib/libRblas.so
LAPACK: /usr/local/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.2.1  goftest_1.1-1  Distance_0.9.8 mrds_2.2.0     moveds_0.1.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2          magrittr_1.5        splines_3.6.1      
 [4] tidyselect_0.2.5    munsell_0.5.0       colorspace_1.4-1   
 [7] lattice_0.20-38     R6_2.4.0            rlang_0.4.0        
[10] Rsolnp_1.16         optimx_2018-7.10    dplyr_0.8.3        
[13] tools_3.6.1         parallel_3.6.1      grid_3.6.1         
[16] nlme_3.1-140        gtable_0.3.0        mgcv_1.8-28        
[19] withr_2.1.2         assertthat_0.2.1    lazyeval_0.2.2     
[22] tibble_2.1.3        numDeriv_2016.8-1.1 crayon_1.3.4       
[25] Matrix_1.2-17       purrr_0.3.2         glue_1.3.1         
[28] compiler_3.6.1      pillar_1.4.2        scales_1.0.0       
[31] truncnorm_1.0-8     pkgconfig_2.0.3 
```
