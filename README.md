Source code for manuscript "A Utility Approach to Individualized Optimal Dose Selection Using Biomarkers" by Pin Li, Jeremy M.G. Taylor, Spring Kong, Shruti Jolly and Matthew J. Schipper. For questions, comments or remarks about the code please contact Pin Li (pinli@umich.edu).

The code has been written using R with the following version information. The packages used are tmvtnorm, stats, ggplot2, glmnet, BB, numDeriv.

To reproduce the Fig 1 and Fig 2, run the file "Fig12.R"

To reproduce the Fig 3 and Table 1, run the file "Simulation.R" with different scenarios by changing the coefficients as in Table 1.

To reproduce the Table A1, run the file "TwoTox.R" with different scenarios by changing the coefficients as in Table A1.


R version 3.5.2 (2018-12-20)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggplot2_3.0.0     numDeriv_2016.8-1 BB_2014.10-1      glmnet_2.0-16    
 [5] foreach_1.4.4     tmvtnorm_1.4-10   gmm_1.6-2         sandwich_2.4-0   
 [9] Matrix_1.2-15     mvtnorm_1.0-10   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1       pillar_1.3.0     compiler_3.5.2   plyr_1.8.4      
 [5] bindr_0.1.1      iterators_1.0.10 tibble_1.4.2     gtable_0.2.0    
 [9] lattice_0.20-38  pkgconfig_2.0.1  rlang_0.3.4      bindrcpp_0.2.2  
[13] withr_2.1.2      dplyr_0.7.6      grid_3.5.2       tidyselect_0.2.4
[17] glue_1.3.0       R6_2.2.2         purrr_0.2.5      magrittr_1.5    
[21] scales_0.5.0     codetools_0.2-15 assertthat_0.2.0 colorspace_1.3-2
[25] quadprog_1.5-5   lazyeval_0.2.1   munsell_0.5.0    crayon_1.3.4    
[29] zoo_1.8-3       



