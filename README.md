Partitionning the metacommunity *CV*<sup>2</sup>
================

This document illustrates how to use the R functions attached to
Segrestin & Lepš (2022) <i>Journal of Ecology</i>. We propose a framework aiming at
disentangling the relative effects of population stability and different
types of synchronies on metacommunity stability.

## Data preparation

The data set should be a `data.frame`. Column names must include
`Community` and `Species`. The annual biomass of individual populations must
be organized in columns (one column per year). No extra column is
allowed.

The following dummy data set shows the required format for a fictitious
metacommunity of 4 communities (C1 to C4) including the biomass
productivity of 19 populations of 9 species (A, B, C, …, I) for 7 years
(Y1 to Y7)

``` r
print(data, row.names = F, right = F)
```

    ##  Community Species Y1 Y2 Y3 Y4 Y5 Y6 Y7
    ##  C1        A       10  8  7 12 15 13 10
    ##  C1        B       20 22 19 16 12 17 19
    ##  C1        C        5  7  4  4  6  2  8
    ##  C1        D        2  0  0  0  1  0  0
    ##  C1        E        0  0  0  0  2  0  0
    ##  C2        A       12 10 11 15 15  8  9
    ##  C2        B       18 20 21 15 17 20 11
    ##  C2        D        0  1  0  0  0  0  0
    ##  C2        F        7  6  3  1  8  4  9
    ##  C2        G        8  7  7  8  7  7  8
    ##  C2        H        0  0  0  0  1  0  0
    ##  C3        B       13 18 20 15  7 21 15
    ##  C3        C        9 12  7  5  9  7  6
    ##  C3        H       20 18 20 22 21 17 20
    ##  C4        A        7  5  6 15 17 12  9
    ##  C4        B       17 22 20 22 10 25 13
    ##  C4        F        3  7  0  0  1  4  7
    ##  C4        H        0  2  0  1  0  0  0
    ##  C4        I       12  7  8  6 10  9  4

This dummy data set can be downloaded [here](https://raw.githubusercontent.com/jsegrestin/MetacomStab/master/data/dummy_data.csv)

## Run the function

You can source the `R` functions from GitHub using the following code

``` r
require(devtools)
url <- "https://raw.githubusercontent.com/jsegrestin/MetacomStab/master/R/functions.R"
devtools::source_url(url)
```

Once the functions and the data set are loaded, you can run `cv2_decomp()` using the
prepared data set

``` r
cv2_decomp(data)
```

    ## 
    ## Decomposition of the metacommunity squared coefficient of variation
    ## See Segrestin & Leps (2022)
    ## 
    ## CV2 = 0.0026, Pop.var = 0.0059, Pop.sync = -0.0033
    ## 
    ##  Segrestin & Leps (2022)      Hammond et al. (2020)
    ##  Pop.sync[direct] = -0.003    Delta = 0.0186       
    ##  Pop.sync[intra] = 0.0068     Beta[MP] = 0.0112    
    ##  pop.sync[indirect] = -0.0082 Beta[CCi] = 0.0485   
    ##  pop.sync[no] = 0.001         Beta[CCno] = 0.0055

It returns the metacommunity *CV*<sup>2</sup>, *Pop.var*,
*Pop.sync* and the four indices of *Pop.sync*
describing population synchrony at different spatial and organizational
scales. The corresponding indices of asynchrony described in [Hammond
*et al.* (2020) *Ecosphere*](https://www.doi.org/10.1002/ecs2.3078)
are also given.

If you want to include a randomization procedure to estimate the range of
*pop.sync* under the hypothesis of independent fluctuations between
populations, use the argument `nrand`

``` r
cv2_decomp(data, nrand = 10000)
```

    ## 
    ## Decomposition of the metacommunity squared coefficient of variation
    ## See Segrestin & Leps (2022)
    ## 
    ## CV2 = 0.0026, Pop.var = 0.0059, Pop.sync = -0.0033 [-0.0044; 0.0059]
    ## 
    ##  Segrestin & Leps (2022)      Rand (95% CI)     Hammond et al. (2020)
    ##  Pop.sync[direct] = -0.003    [-0.0025; 0.0029] Delta = 0.0186       
    ##  Pop.sync[intra] = 0.0068     [-0.0024; 0.0031] Beta[MP] = 0.0112    
    ##  pop.sync[indirect] = -0.0082 [-0.004; 0.0045]  Beta[CCi] = 0.0485   
    ##  pop.sync[no] = 0.001         [-0.0015; 0.0014] Beta[CCno] = 0.0055

It now includes the 95% confidence interval of *Pop.sync*
values corresponding to independent fluctuations between populations.
Observed values out of the range are significantly different from zero. 
Positive or negative values represent prevalent synchrony or anti-synchrony between populations, respectively.
