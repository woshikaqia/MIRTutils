MIRTutils package
================

## Overview

This package provides a set of functions to calculate some commonly used
values when doing IRT analysis, especially with new item types such as
item bundles (testlets).

IRT model currently supported are 1-3 PL, GPCM, Rasch Testlet Model, and
a mix of these models. Values that can be computed include the
probability of correct response, item information, item score residuals,
person scores, person fit index, etc. It also provides a function for
multidimensional IRT data simulation.

In the context of the current package, test items can be either
standalone items (**SA**) or cluster items (**Cluster**). Both SA and
cluster items consist of ***assertions***. An assertion is the smallest
item scoring unit in the test. You can think of assertion as the same as
a traditional item that is typically modeled by a unidimensional IRT
model. An SA item typically includes 1 to a handful of assertions,
whereas a cluster item typically includes 6 or more assertions. An SA
item is more like a small set of traditional items assuming local
dependence among assertions, whereas a cluster is essentially a testlet.

The IRT model used is a special model similar to a Rasch Testlet model
but allows for SA items to load only on the overall dimension. When
there are only SA items, the model reduces to a regular unidimensional
IRT model, and all the values computed in this package fits to the
traditional unidimensional IRT paradigm. When there are only cluster
items, the model reduces to a Rasch testlet model. Unlike the
unidimensional model where theta is usually estimated by the univariate
MLE, EAP or MAP, and the multidimensional model where theta is usually
estimated by the multivarite MLE or EAP, most of the values computed by
this package involving testlets is based on the marginal maximum
likelihood estimation (MMLE) of theta. MMLE is a hybrid of MLE and EAP,
where the nuisance dimension of a testlet (cluster) is first integrated
out, and the resulting marginal likelihood is then maximized to find the
theta estimate.

## Installation

`devtools::install_github("woshikaqia/MIRTutils")`

## Github Repo

<https://github.com/woshikaqia/MIRTutils/>

## Examples

### Simulate Data

Simulate 7 examniees’ item responses to a test consist of 20 SA items
(each with only 1 assertion; 17 of them are dichotomous and 3 are
polytomous) and 7 Cluster items (with a total of 81 assertions).

``` r
library(MIRTutils)
data(example_SA_parm)
data(example_Cluster_parm)
sigma <- diag(c(1, sqrt(unique(example_Cluster_parm$cluster_var))))
mu <- rep(0, nrow(sigma))
thetas <- MASS::mvrnorm(7,mu,sigma)
thetas[,1] <- seq(-3,3,1) #overall dimension theta values
itmDat <- sim_data(thetas = thetas, SA_parm = example_SA_parm, Cluster_parm = example_Cluster_parm)
```

``` r
library(dplyr)
as_tibble(example_SA_parm)
```

    ## # A tibble: 20 × 7
    ##        a      b1    b2    b3     g ItemID AssertionID
    ##    <dbl>   <dbl> <dbl> <dbl> <dbl> <chr>  <chr>      
    ##  1   1    0.588     NA  NA       0 SA_i1  SA_i1_AS1  
    ##  2   1    0.423     NA  NA       0 SA_i2  SA_i2_AS1  
    ##  3   1   -1.06      NA  NA       0 SA_i3  SA_i3_AS1  
    ##  4   1   -0.156     NA  NA       0 SA_i4  SA_i4_AS1  
    ##  5   1   -0.416     NA  NA       0 SA_i5  SA_i5_AS1  
    ##  6   1   -0.135     NA  NA       0 SA_i6  SA_i6_AS1  
    ##  7   1   -1.14      NA  NA       0 SA_i7  SA_i7_AS1  
    ##  8   1   -0.573     NA  NA       0 SA_i8  SA_i8_AS1  
    ##  9   1   -0.859     NA  NA       0 SA_i9  SA_i9_AS1  
    ## 10   1   -0.505     NA  NA       0 SA_i10 SA_i10_AS1 
    ## 11   1    0.297     NA  NA       0 SA_i11 SA_i11_AS1 
    ## 12   1    0.404     NA  NA       0 SA_i12 SA_i12_AS1 
    ## 13   1   -0.386     NA  NA       0 SA_i13 SA_i13_AS1 
    ## 14   1   -0.540     NA  NA       0 SA_i14 SA_i14_AS1 
    ## 15   1   -1.45      NA  NA       0 SA_i15 SA_i15_AS1 
    ## 16   1    0.170     NA  NA       0 SA_i16 SA_i16_AS1 
    ## 17   1    0.0839    NA  NA       0 SA_i17 SA_i17_AS1 
    ## 18   1   -1          0   0.5    NA CR_i1  CR_i1_AS1  
    ## 19   1.5  0          1   1.5    NA CR_i2  CR_i2_AS1  
    ## 20   2    1          2  NA      NA CR_i3  CR_i3_AS1

``` r
as_tibble(example_Cluster_parm)
```

    ## # A tibble: 81 × 6
    ##        a       b cluster_var position ItemID AssertionID
    ##    <dbl>   <dbl>       <dbl>    <dbl> <chr>  <chr>      
    ##  1     1 -1.21         0.371        1 CL_i1  CL_i1_AS1  
    ##  2     1 -0.0959       0.371        1 CL_i1  CL_i1_AS2  
    ##  3     1  0.0881       0.371        1 CL_i1  CL_i1_AS3  
    ##  4     1  0.0188       0.371        1 CL_i1  CL_i1_AS4  
    ##  5     1  0.248        0.371        1 CL_i1  CL_i1_AS5  
    ##  6     1 -0.0965       0.371        1 CL_i1  CL_i1_AS6  
    ##  7     1 -0.546        0.371        1 CL_i1  CL_i1_AS7  
    ##  8     1 -0.727        0.371        1 CL_i1  CL_i1_AS8  
    ##  9     1  0.787        0.371        1 CL_i1  CL_i1_AS9  
    ## 10     1 -0.210        0.371        1 CL_i1  CL_i1_AS10 
    ## # ℹ 71 more rows

``` r
unique(example_Cluster_parm$ItemID)
```

    ## [1] "CL_i1" "CL_i2" "CL_i3" "CL_i4" "CL_i5" "CL_i6" "CL_i7"

``` r
as_tibble(itmDat)
```

    ## # A tibble: 7 × 101
    ##   SA_i1_AS1 SA_i2_AS1 SA_i3_AS1 SA_i4_AS1 SA_i5_AS1 SA_i6_AS1 SA_i7_AS1
    ##       <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
    ## 1         0         0         0         0         0         0         0
    ## 2         0         0         1         0         1         0         0
    ## 3         0         0         0         1         0         0         1
    ## 4         1         1         1         1         0         1         1
    ## 5         1         0         1         1         1         1         1
    ## 6         1         1         1         1         1         1         1
    ## 7         1         1         1         1         1         1         1
    ## # ℹ 94 more variables: SA_i8_AS1 <dbl>, SA_i9_AS1 <dbl>, SA_i10_AS1 <dbl>,
    ## #   SA_i11_AS1 <dbl>, SA_i12_AS1 <dbl>, SA_i13_AS1 <dbl>, SA_i14_AS1 <dbl>,
    ## #   SA_i15_AS1 <dbl>, SA_i16_AS1 <dbl>, SA_i17_AS1 <dbl>, CR_i1_AS1 <dbl>,
    ## #   CR_i2_AS1 <dbl>, CR_i3_AS1 <dbl>, CL_i1_AS1 <dbl>, CL_i1_AS2 <dbl>,
    ## #   CL_i1_AS3 <dbl>, CL_i1_AS4 <dbl>, CL_i1_AS5 <dbl>, CL_i1_AS6 <dbl>,
    ## #   CL_i1_AS7 <dbl>, CL_i1_AS8 <dbl>, CL_i1_AS9 <dbl>, CL_i1_AS10 <dbl>,
    ## #   CL_i1_AS11 <dbl>, CL_i1_AS12 <dbl>, CL_i1_AS13 <dbl>, CL_i1_AS14 <dbl>, …

### Scoring

…continued from the data simulation, this part scores the first
examinee. Note that the `scoring` function should be run one person at a
time. In general you can construct a loop to score all of the examinees,
or do parallel processing if you are dealing with extremely a large
dataset

``` r
SA_dat <- itmDat[,1:20]
Cluster_dat <- itmDat[,-1:-20]
# Score the first examinee
out_scoring <- scoring(SA_dat[1,], Cluster_dat[1,], example_SA_parm, example_Cluster_parm, n.nodes = 11, SE=TRUE)
(est_theta <- out_scoring$par)  # estimated theta
```

    ## [1] -3.081139

``` r
as.vector(out_scoring$SE)  # estimated standard error
```

    ## [1] 0.4445252

### Person Fit

…continued from the data simulation and scoring, this part computes the
person fit statistics for the first examinee.

- For traditional test (i.e., standalone items only), this is lz\*
  described in Snijders (2001) or Sinharay (2016).
- For Cluster items, this is lzt\* described in Lin et.al (2024), which
  extends the lz\* to the Rasch testlet model when theta is estimated by
  MMLE.

You can specify one ore more alpha levels using the `Alpha` argument for
flagging.

``` r
# Compute person fit for the first examinee
rst <- person_fit(est_theta, SA_dat[1,], Cluster_dat[1,], example_SA_parm, example_Cluster_parm, Alpha=c(0.01,0.05))
rst
```

    ##       theta raw        ll    exp_ll var_ll_org correction   var_ll        lz
    ## 1 -3.081139  14 -28.48542 -31.13043   46.16721   36.35753 9.809683 0.8445016
    ##   test_info flag_0.01 flag_0.05
    ## 1  4.265769         0         0

### Other Useful Values

You can use the `utility()` function to compute some useful values if
needed. The example below computes the **expected item raw score** and
**item information**.

``` r
# Compute person fit for the first examinee
data(example_SA_parm)
data(example_Cluster_parm)
theta <- seq(-3,3,1)
output <- utility(theta, example_SA_parm, example_Cluster_parm, what = c("escore","iteminfo"))
round(output$escore,3)
```

    ##   theta SA_i1_AS1 SA_i2_AS1 SA_i3_AS1 SA_i4_AS1 SA_i5_AS1 SA_i6_AS1 SA_i7_AS1
    ## 1    -3     0.027     0.032     0.126     0.055     0.070     0.054     0.134
    ## 2    -2     0.070     0.081     0.282     0.137     0.170     0.134     0.296
    ## 3    -1     0.170     0.194     0.516     0.301     0.358     0.296     0.534
    ## 4     0     0.357     0.396     0.743     0.539     0.603     0.534     0.757
    ## 5     1     0.602     0.640     0.887     0.761     0.805     0.757     0.894
    ## 6     2     0.804     0.829     0.955     0.896     0.918     0.894     0.958
    ## 7     3     0.918     0.929     0.983     0.959     0.968     0.958     0.984
    ##   SA_i8_AS1 SA_i9_AS1 SA_i10_AS1 SA_i11_AS1 SA_i12_AS1 SA_i13_AS1 SA_i14_AS1
    ## 1     0.081     0.105      0.076      0.036      0.032      0.068      0.079
    ## 2     0.194     0.242      0.183      0.091      0.083      0.166      0.188
    ## 3     0.395     0.465      0.379      0.215      0.197      0.351      0.387
    ## 4     0.639     0.702      0.624      0.426      0.400      0.595      0.632
    ## 5     0.828     0.865      0.818      0.669      0.645      0.800      0.823
    ## 6     0.929     0.946      0.924      0.846      0.831      0.916      0.927
    ## 7     0.973     0.979      0.971      0.937      0.931      0.967      0.972
    ##   SA_i15_AS1 SA_i16_AS1 SA_i17_AS1 CR_i1_AS1 CR_i2_AS1 CR_i3_AS1  CL_i1 CL_i2
    ## 1      0.175      0.040      0.044     0.131     0.011     0.000  1.105 2.008
    ## 2      0.365      0.102      0.111     0.337     0.048     0.002  2.582 3.446
    ## 3      0.610      0.237      0.253     0.809     0.199     0.018  5.271 5.085
    ## 4      0.810      0.458      0.479     1.620     0.675     0.123  8.998 6.574
    ## 5      0.920      0.696      0.714     2.385     1.639     0.595 12.731 7.678
    ## 6      0.969      0.862      0.872     2.770     2.558     1.405 15.427 8.361
    ## 7      0.988      0.944      0.949     2.917     2.895     1.877 16.904 8.722
    ##    CL_i3 CL_i4  CL_i5 CL_i6  CL_i7
    ## 1  2.887 0.369  1.401 0.705  4.796
    ## 2  4.677 0.877  2.754 1.507  6.571
    ## 3  6.689 1.820  4.773 2.732  8.160
    ## 4  8.601 3.182  7.294 4.148  9.366
    ## 5 10.177 4.714  9.895 5.387 10.152
    ## 6 11.336 6.088 12.077 6.224 10.598
    ## 7 12.103 7.059 13.552 6.671 10.824

``` r
round(output$iteminfo,4)
```

    ##   theta SA_i1_AS1 SA_i2_AS1 SA_i3_AS1 SA_i4_AS1 SA_i5_AS1 SA_i6_AS1 SA_i7_AS1
    ## 1    -3    0.0262    0.0306    0.1102    0.0520    0.0653    0.0510    0.1162
    ## 2    -2    0.0650    0.0748    0.2023    0.1179    0.1413    0.1161    0.2086
    ## 3    -1    0.1409    0.1565    0.2497    0.2103    0.2299    0.2085    0.2488
    ## 4     0    0.2296    0.2392    0.1908    0.2485    0.2395    0.2489    0.1840
    ## 5     1    0.2397    0.2303    0.1000    0.1821    0.1571    0.1841    0.0945
    ## 6     2    0.1575    0.1419    0.0426    0.0930    0.0752    0.0946    0.0399
    ## 7     3    0.0755    0.0656    0.0166    0.0392    0.0308    0.0399    0.0155
    ##   SA_i8_AS1 SA_i9_AS1 SA_i10_AS1 SA_i11_AS1 SA_i12_AS1 SA_i13_AS1 SA_i14_AS1
    ## 1    0.0745    0.0941     0.0704     0.0344     0.0311     0.0636     0.0725
    ## 2    0.1561    0.1835     0.1496     0.0830     0.0760     0.1385     0.1529
    ## 3    0.2389    0.2488     0.2353     0.1686     0.1583     0.2279     0.2372
    ## 4    0.2306    0.2090     0.2347     0.2446     0.2401     0.2409     0.2326
    ## 5    0.1423    0.1166     0.1487     0.2215     0.2291     0.1600     0.1454
    ## 6    0.0659    0.0513     0.0698     0.1303     0.1401     0.0771     0.0678
    ## 7    0.0266    0.0202     0.0283     0.0588     0.0646     0.0317     0.0274
    ##   SA_i15_AS1 SA_i16_AS1 SA_i17_AS1 CR_i1_AS1 CR_i2_AS1 CR_i3_AS1  CL_i1  CL_i2
    ## 1     0.1442     0.0387     0.0419    0.1266    0.0247    0.0013 0.7032 0.5083
    ## 2     0.2319     0.0920     0.0984    0.3108    0.1062    0.0099 1.1321 0.5761
    ## 3     0.2379     0.1807     0.1889    0.6559    0.4025    0.0713 1.4725 0.5808
    ## 4     0.1541     0.2482     0.2496    0.8907    1.0817    0.4498 1.5988 0.5277
    ## 5     0.0732     0.2115     0.2041    0.5742    1.6803    1.4709 1.4745 0.4279
    ## 6     0.0299     0.1191     0.1118    0.2322    0.8992    1.4709 1.1328 0.3010
    ## 7     0.0114     0.0526     0.0487    0.0843    0.2325    0.4498 0.7012 0.1788
    ##    CL_i3  CL_i4  CL_i5  CL_i6  CL_i7
    ## 1 0.3958 0.3106 0.4971 0.4624 0.3013
    ## 2 0.4296 0.6040 0.6373 0.7147 0.2991
    ## 3 0.4355 0.9235 0.7278 0.8775 0.2787
    ## 4 0.4191 1.1144 0.7640 0.8799 0.2415
    ## 5 0.3841 1.1173 0.7434 0.7283 0.1913
    ## 6 0.3316 0.9363 0.6588 0.4871 0.1358
    ## 7 0.2631 0.6274 0.5125 0.2611 0.0848

## References
 - Lin, Z., Jiang, T., Rijmen, F., & Van Wamelen, P. (2024). Asymptotically Correct Person Fit z-Statistics For the Rasch Testlet Model. _Psychometrika_,  https://doi.org/10.1007/s11336-024-09997-y
 - Sinharay, S. (2016). Asymptotically correct standardization of person-fit statistics beyond dichotomous items. _Psychometrika_, 81(4), 992-1013.
 - Snijders, T. A. (2001). Asymptotic null distribution of person fit statistics with estimated person parameter. _Psychometrika_, 66, 331-342.
