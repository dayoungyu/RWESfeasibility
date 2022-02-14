
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RWESfeasibility

The goal of RWESfeasibility is to run simulations for finding CIs of
extrapolated feasibility counts with weighted proportions.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dayoungyu/RWESfeasibility")
```

## Example

``` r
library(RWESfeasibility)

n_cohorts = 3
n_den = c(160, 250, 100)
n_rev = c(50, 50, 45)
n_case = c(32, 2, 19)
N_sim = 1000

sim_results = ci_sim(n_cohorts, n_den, n_rev, n_case, N_sim)
```

##### Summary table output

``` r
sim_results$summary_table
#> # A tibble: 4 x 2
#>   Statistic          Frequency                
#>   <fct>              <chr>                    
#> 1 Mean               154 (30.3%)              
#> 2 95% CI             135 (26.4%) - 174 (34.1%)
#> 3 90% probability N≥ 141 (27.7%)              
#> 4 75% probability N≥ 148 (29.1%)
```
