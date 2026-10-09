# Advanced Topics

## Introduction

Most `EpiExposure` workflows can be completed using the default model
specification introduced in previous sections

However, advanced applications often require additional flexibility
regarding:

- model engines;
- uncertainty estimation;
- spline specification;
- identifiability diagnostics;
- computational performance.

This section discusses these advanced topics and provides
recommendations for robust DLNM analyses.

## Alternative model engines

`EpiExposure` separates exposure-lag specification from model
estimation.

Consequently, the same DLNM structure can be fitted using different
regression engines.

## Supported families by modeling engine

The modeling engines available in
[`fit_epidlnm()`](https://tomazrg.github.io/EpiExposure/reference/fit_epidlnm.md)
and
[`find_bestfit()`](https://tomazrg.github.io/EpiExposure/reference/find_bestfit.md)
differ in the response families exposed through the `EpiExposure`
interface. The table below summarizes the current `EpiExposure` support.
It describes the package interface rather than every family accepted
natively by the underlying modeling packages.

| Engine    | Beta | Binomial       | Poisson | Gamma | Gaussian | Negative binomial |
|:----------|:----:|:---------------|:-------:|:-----:|:--------:|:-----------------:|
| `glm`     |  No  | Yes            |   Yes   |  Yes  |   Yes    |        No         |
| `glmmTMB` | Yes  | Yes            |   Yes   |  Yes  |   Yes    |     Yes, NB2      |
| `gam`     | Yes  | Yes            |   Yes   |  Yes  |   Yes    |     Yes, NB2      |
| `gamm`    |  No  | Yes            |   Yes   |  Yes  |   Yes    |        No         |
| `gls`     |  No  | No             |   No    |  No   |   Yes    |        No         |
| `spaMM`   | Yes  | Yes            |   Yes   |  Yes  |   Yes    |     Yes, NB2      |
| `brms`    | Yes  | Yes, Bernoulli |   Yes   |  Yes  |   Yes    |     Yes, NB2      |
| `INLA`    | Yes  | Yes            |   Yes   |  Yes  |   Yes    |     Yes, NB2      |
| `bdlnm`   | Yes  | Yes            |   Yes   |  Yes  |   Yes    |     Yes, NB2      |

**Note:** Negative-binomial support refers to the quadratic-variance NB2
parameterization. The table describes the families currently exposed
through the `EpiExposure` interface, not every family that may be
available natively in the underlying modeling packages. For INLA-backed
engines, availability also depends on the likelihoods provided by the
installed INLA version.

Because family specification follows the same
[`fit_epidlnm()`](https://tomazrg.github.io/EpiExposure/reference/fit_epidlnm.md)
interface across engines, the sections below do not repeat every
supported engine-family combination. Instead, each example highlights a
capability that differs meaningfully across engines, such as
conventional random effects, spatial autocorrelation, or Bayesian
inference.

## Reference DLNM design

The examples below use the same cross-basis specification whenever
possible. This allows the modeling engines to be compared without
repeatedly rebuilding the exposure-lag design.

``` r

library(EpiExposure)
```

``` r

data("poisson_data")

head(poisson_data)
#>   epi_id block time    tmean  wetness      rain y
#> 1      1   B01    0 25.32685 14.55902 11.629805 5
#> 2      1   B01    1 25.72697 13.35332  0.000000 5
#> 3      1   B01    2 24.88694 13.00895  8.036259 5
#> 4      1   B01    3 26.04031 13.22694  0.000000 5
#> 5      1   B01    4 25.94690 15.30767  6.709351 5
#> 6      1   B01    5 24.85839 17.41708  8.826920 5
```

``` r

cb <- define_exposures(
  data = poisson_data,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)

dat <- build_design(
  data = poisson_data,
  cb_templates = cb,
  groups = c(
    "epi_id",
    "block"
  )
)
```

``` r

dat_poisson <- prepare_response(
  data = dat,
  response = "y",
  family = "poisson"
)
```

## Fixed-effect reference model

A generalized linear model provides the simplest fixed-effect
specification. It is useful as a reference when random or spatial
dependence is not required.

``` r

fit_glm <- fit_epidlnm(
  data = dat_poisson,
  model_engine = "glm",
  family = "poisson",
  epiexposure_spec = attr(
    cb,
    "spec"
  )
)
```

``` r

summary(fit_glm)
#> 
#> Call:
#> stats::glm(formula = y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + 
#>     cb_tmean_3 + cb_tmean_4 + cb_tmean_5 + cb_tmean_6 + cb_rain_1 + 
#>     cb_rain_2 + cb_rain_3 + cb_rain_4 + cb_rain_5 + cb_rain_6 + 
#>     cb_wetness_1 + cb_wetness_2 + cb_wetness_3 + cb_wetness_4 + 
#>     cb_wetness_5 + cb_wetness_6, family = "poisson", data = dat_poisson)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -3.303652   1.132719  -2.917 0.003539 ** 
#> cb_tmean_1    0.090028   0.013837   6.506 7.71e-11 ***
#> cb_tmean_2    0.012584   0.025493   0.494 0.621554    
#> cb_tmean_3    0.140732   0.043795   3.213 0.001312 ** 
#> cb_tmean_4    0.270869   0.097282   2.784 0.005363 ** 
#> cb_tmean_5   -0.027014   0.016290  -1.658 0.097252 .  
#> cb_tmean_6    0.162029   0.039149   4.139 3.49e-05 ***
#> cb_rain_1     0.091313   0.025531   3.577 0.000348 ***
#> cb_rain_2    -0.259121   0.033323  -7.776 7.48e-15 ***
#> cb_rain_3     0.078788   0.015872   4.964 6.91e-07 ***
#> cb_rain_4    -0.064189   0.022854  -2.809 0.004975 ** 
#> cb_rain_5     0.001674   0.036894   0.045 0.963806    
#> cb_rain_6     0.035621   0.050073   0.711 0.476846    
#> cb_wetness_1  0.051027   0.003629  14.062  < 2e-16 ***
#> cb_wetness_2 -0.086751   0.010892  -7.965 1.65e-15 ***
#> cb_wetness_3  0.049067   0.014053   3.492 0.000480 ***
#> cb_wetness_4 -0.098276   0.036715  -2.677 0.007434 ** 
#> cb_wetness_5  0.059851   0.009419   6.355 2.09e-10 ***
#> cb_wetness_6 -0.084802   0.018470  -4.591 4.40e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 2866.42  on 519  degrees of freedom
#> Residual deviance:  955.26  on 501  degrees of freedom
#> AIC: 3133.2
#> 
#> Number of Fisher Scoring iterations: 4
```

## Conventional random effects

When epidemics are grouped within experimental blocks, observations from
the same block may share unmeasured conditions. A random intercept can
account for this hierarchical dependence.

``` r

fit_glmm <- fit_epidlnm(
  data = dat_poisson,
  model_engine = "glmmTMB",
  family = "poisson",
  random_effect = "block",
  epiexposure_spec = attr(
    cb,
    "spec"
  )
)
```

``` r

summary(fit_glmm)
#>  Family: poisson  ( log )
#> Formula:          
#> y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 +  
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 +  
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 +  
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 +  
#>     (1 | block)
#> Data: dat_poisson
#> 
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>    2920.3    3005.4   -1440.2    2880.3       500 
#> 
#> Random effects:
#> 
#> Conditional model:
#>  Groups Name        Variance Std.Dev.
#>  block  (Intercept) 0.04283  0.207   
#> Number of obs: 520, groups:  block, 20
#> 
#> Conditional model:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -2.574612   1.158675  -2.222 0.026281 *  
#> cb_tmean_1    0.081602   0.014101   5.787 7.17e-09 ***
#> cb_tmean_2   -0.008901   0.025844  -0.344 0.730544    
#> cb_tmean_3    0.115294   0.044915   2.567 0.010260 *  
#> cb_tmean_4    0.170722   0.098491   1.733 0.083031 .  
#> cb_tmean_5   -0.023288   0.017047  -1.366 0.171890    
#> cb_tmean_6    0.119456   0.039905   2.994 0.002758 ** 
#> cb_rain_1     0.129022   0.026462   4.876 1.08e-06 ***
#> cb_rain_2    -0.239983   0.033977  -7.063 1.63e-12 ***
#> cb_rain_3     0.056071   0.016278   3.445 0.000572 ***
#> cb_rain_4    -0.057666   0.023722  -2.431 0.015063 *  
#> cb_rain_5    -0.038092   0.037558  -1.014 0.310473    
#> cb_rain_6     0.060427   0.052017   1.162 0.245368    
#> cb_wetness_1  0.052540   0.003729  14.091  < 2e-16 ***
#> cb_wetness_2 -0.093637   0.011114  -8.425  < 2e-16 ***
#> cb_wetness_3  0.053720   0.014580   3.685 0.000229 ***
#> cb_wetness_4 -0.093537   0.037422  -2.499 0.012437 *  
#> cb_wetness_5  0.055911   0.009722   5.751 8.87e-09 ***
#> cb_wetness_6 -0.090618   0.018671  -4.853 1.21e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The term generated internally is equivalent to:

`(1 | block)`

This effect represents discrete group-level heterogeneity. It does not
model spatial correlation as a function of distance.

### GAM estimation

The `gam` engine can fit the DLNM design with `mgcv`. Beta and
negative-binomial models use `mgcv` extended families and are fitted
through the supported REML route by default.

``` r

data("nb2_data")

head(nb2_data)
#>   epi_id block time    tmean  wetness      rain y
#> 1      1   B01    0 28.51652 3.304398  9.234377 7
#> 2      1   B01    1 25.78990 5.600390 19.225955 7
#> 3      1   B01    2 26.86089 7.060567  0.000000 7
#> 4      1   B01    3 25.19746 5.064939  8.700396 7
#> 5      1   B01    4 25.80954 5.634410  0.000000 7
#> 6      1   B01    5 25.28982 6.813331  5.096005 7
```

``` r

dat_nb2 <- build_design(
  data = nb2_data,
  cb_templates = cb,
  groups = c(
    "epi_id",
    "block"
  )
)

dat_nb2 <- prepare_response(
  data = dat_nb2,
  response = "y",
  family = "negative_binomial"
)
```

``` r

fit_gam <- fit_epidlnm(
  data = dat_nb2,
  model_engine = "gam",
  family = "negative_binomial",
  epiexposure_spec = attr(
    cb,
    "spec"
  )
)
```

``` r

summary(fit_gam)
#> 
#> Family: Negative Binomial(3.312) 
#> Link function: log 
#> 
#> Formula:
#> y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + 
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + 
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + 
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6
#> 
#> Parametric coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.184097   2.303002  -0.080 0.936287    
#> cb_tmean_1    0.050136   0.028851   1.738 0.082253 .  
#> cb_tmean_2   -0.129826   0.052604  -2.468 0.013588 *  
#> cb_tmean_3    0.009794   0.089261   0.110 0.912633    
#> cb_tmean_4   -0.166707   0.200448  -0.832 0.405595    
#> cb_tmean_5   -0.052436   0.036715  -1.428 0.153242    
#> cb_tmean_6    0.012849   0.082609   0.156 0.876399    
#> cb_rain_1     0.198301   0.057344   3.458 0.000544 ***
#> cb_rain_2    -0.328951   0.071447  -4.604 4.14e-06 ***
#> cb_rain_3     0.027017   0.036851   0.733 0.463466    
#> cb_rain_4    -0.094688   0.050348  -1.881 0.060015 .  
#> cb_rain_5    -0.085088   0.085783  -0.992 0.321251    
#> cb_rain_6    -0.085474   0.106752  -0.801 0.423318    
#> cb_wetness_1  0.066901   0.006780   9.868  < 2e-16 ***
#> cb_wetness_2 -0.082650   0.021324  -3.876 0.000106 ***
#> cb_wetness_3  0.094385   0.025050   3.768 0.000165 ***
#> cb_wetness_4 -0.012930   0.072911  -0.177 0.859238    
#> cb_wetness_5  0.058337   0.016182   3.605 0.000312 ***
#> cb_wetness_6 -0.036181   0.044278  -0.817 0.413850    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> R-sq.(adj) =  0.467   Deviance explained =   55%
#> -REML = 1848.3  Scale est. = 1         n = 520
```

### GAMM estimation

The `gamm` engine combines the fixed DLNM design with a conventional
random intercept. In the current `EpiExposure` interface,
`random_effect` is required.

``` r

fit_gamm <- fit_epidlnm(
  data = dat_poisson,
  model_engine = "gamm",
  family = "poisson",
  random_effect = "block",
  epiexposure_spec = attr(
    cb,
    "spec"
  )
)
#> 
#>  Maximum number of PQL iterations:  20
```

The fitted object contains two native components:

``` r

summary(fit_gamm)
#> 
#> Family: poisson 
#> Link function: log 
#> 
#> Formula:
#> y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + 
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + 
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + 
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6
#> 
#> Parametric coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  -2.573106   1.180352  -2.180 0.029726 *  
#> cb_tmean_1    0.081602   0.014365   5.681 2.28e-08 ***
#> cb_tmean_2   -0.008901   0.026325  -0.338 0.735418    
#> cb_tmean_3    0.115297   0.045756   2.520 0.012051 *  
#> cb_tmean_4    0.170724   0.100318   1.702 0.089407 .  
#> cb_tmean_5   -0.023289   0.017367  -1.341 0.180525    
#> cb_tmean_6    0.119458   0.040642   2.939 0.003442 ** 
#> cb_rain_1     0.129018   0.026943   4.789 2.21e-06 ***
#> cb_rain_2    -0.239982   0.034612  -6.934 1.27e-11 ***
#> cb_rain_3     0.056073   0.016573   3.383 0.000772 ***
#> cb_rain_4    -0.057667   0.024167  -2.386 0.017396 *  
#> cb_rain_5    -0.038088   0.038251  -0.996 0.319851    
#> cb_rain_6     0.060423   0.052990   1.140 0.254717    
#> cb_wetness_1  0.052540   0.003798  13.832  < 2e-16 ***
#> cb_wetness_2 -0.093637   0.011321  -8.271 1.20e-15 ***
#> cb_wetness_3  0.053718   0.014853   3.617 0.000329 ***
#> cb_wetness_4 -0.093536   0.038125  -2.453 0.014490 *  
#> cb_wetness_5  0.055910   0.009904   5.645 2.77e-08 ***
#> cb_wetness_6 -0.090616   0.019022  -4.764 2.49e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> R-sq.(adj) =  0.631   
#>   Scale est. = 1         n = 520
```

``` r

summary(
  fit_gamm,
  component = "lme"
)
#> Linear mixed-effects model fit by maximum likelihood
#>   Data: data 
#>        AIC     BIC    logLik
#>   524.4704 609.547 -242.2352
#> 
#> Random effects:
#>  Formula: ~1 | block
#>         (Intercept) Residual
#> StdDev:   0.2066749        1
#> 
#> Variance function:
#>  Structure: fixed weights
#>  Formula: ~invwt 
#> Fixed effects:  list(fixed) 
#>                    Value Std.Error  DF   t-value p-value
#> X(Intercept)  -2.5731056 1.1803515 482 -2.179949  0.0297
#> Xcb_tmean_1    0.0816015 0.0143649 482  5.680623  0.0000
#> Xcb_tmean_2   -0.0089008 0.0263249 482 -0.338115  0.7354
#> Xcb_tmean_3    0.1152970 0.0457561 482  2.519817  0.0121
#> Xcb_tmean_4    0.1707239 0.1003175 482  1.701835  0.0894
#> Xcb_tmean_5   -0.0232888 0.0173666 482 -1.341008  0.1805
#> Xcb_tmean_6    0.1194585 0.0406424 482  2.939258  0.0034
#> Xcb_rain_1     0.1290176 0.0269427 482  4.788593  0.0000
#> Xcb_rain_2    -0.2399821 0.0346118 482 -6.933535  0.0000
#> Xcb_rain_3     0.0560726 0.0165732 482  3.383323  0.0008
#> Xcb_rain_4    -0.0576668 0.0241673 482 -2.386149  0.0174
#> Xcb_rain_5    -0.0380883 0.0382508 482 -0.995753  0.3199
#> Xcb_rain_6     0.0604231 0.0529900 482  1.140273  0.2547
#> Xcb_wetness_1  0.0525396 0.0037983 482 13.832359  0.0000
#> Xcb_wetness_2 -0.0936366 0.0113209 482 -8.271142  0.0000
#> Xcb_wetness_3  0.0537184 0.0148528 482  3.616711  0.0003
#> Xcb_wetness_4 -0.0935359 0.0381246 482 -2.453426  0.0145
#> Xcb_wetness_5  0.0559104 0.0099040 482  5.645222  0.0000
#> Xcb_wetness_6 -0.0906163 0.0190219 482 -4.763798  0.0000
#>  Correlation: 
#>               X(Int) Xcb_t_1 Xcb_t_2 Xcb_t_3 Xcb_t_4 Xcb_t_5 Xcb_t_6 Xcb_r_1
#> Xcb_tmean_1   -0.879                                                        
#> Xcb_tmean_2   -0.730  0.756                                                 
#> Xcb_tmean_3   -0.945  0.821   0.673                                         
#> Xcb_tmean_4   -0.866  0.910   0.863   0.773                                 
#> Xcb_tmean_5   -0.381  0.210   0.420   0.559   0.147                         
#> Xcb_tmean_6   -0.664  0.655   0.596   0.513   0.845  -0.279                 
#> Xcb_rain_1     0.012 -0.028  -0.016  -0.002  -0.021   0.022  -0.013         
#> Xcb_rain_2     0.032 -0.074  -0.071  -0.021  -0.077   0.021  -0.056  -0.076 
#> Xcb_rain_3    -0.063  0.022  -0.051   0.052  -0.001   0.036  -0.027  -0.662 
#> Xcb_rain_4    -0.004  0.061   0.057   0.001   0.057  -0.015   0.032   0.033 
#> Xcb_rain_5    -0.034  0.001  -0.048   0.035  -0.017   0.040  -0.041  -0.602 
#> Xcb_rain_6     0.005  0.052   0.031  -0.015   0.040  -0.033   0.028   0.050 
#> Xcb_wetness_1  0.002 -0.062  -0.036  -0.048  -0.047  -0.018  -0.016   0.054 
#> Xcb_wetness_2 -0.023 -0.011   0.002  -0.038  -0.014  -0.035  -0.023   0.002 
#> Xcb_wetness_3 -0.198 -0.043  -0.035   0.007  -0.027   0.023  -0.003  -0.045 
#> Xcb_wetness_4 -0.048 -0.047   0.000  -0.046  -0.052   0.019  -0.059   0.015 
#> Xcb_wetness_5 -0.122  0.006  -0.013   0.009   0.010  -0.008   0.026  -0.097 
#> Xcb_wetness_6 -0.033 -0.002  -0.014  -0.002  -0.043   0.016  -0.060   0.048 
#>               Xcb_r_2 Xcb_r_3 Xcb_r_4 Xcb_r_5 Xcb_r_6 Xcb_w_1 Xcb_w_2 Xcb_w_3
#> Xcb_tmean_1                                                                  
#> Xcb_tmean_2                                                                  
#> Xcb_tmean_3                                                                  
#> Xcb_tmean_4                                                                  
#> Xcb_tmean_5                                                                  
#> Xcb_tmean_6                                                                  
#> Xcb_rain_1                                                                   
#> Xcb_rain_2                                                                   
#> Xcb_rain_3     0.040                                                         
#> Xcb_rain_4    -0.605  -0.011                                                 
#> Xcb_rain_5     0.025   0.938   0.020                                         
#> Xcb_rain_6    -0.669  -0.010   0.866  -0.004                                 
#> Xcb_wetness_1  0.046  -0.102   0.058  -0.108   0.027                         
#> Xcb_wetness_2 -0.045   0.053  -0.011   0.072  -0.012  -0.109                 
#> Xcb_wetness_3  0.029   0.113  -0.097   0.108  -0.055  -0.149   0.122         
#> Xcb_wetness_4  0.048  -0.035  -0.045  -0.016  -0.036   0.381   0.470   0.106 
#> Xcb_wetness_5  0.017   0.153  -0.122   0.148  -0.073  -0.579   0.203   0.785 
#> Xcb_wetness_6  0.051  -0.076  -0.008  -0.061   0.007   0.358   0.128  -0.118 
#>               Xcb_w_4 Xcb_w_5
#> Xcb_tmean_1                  
#> Xcb_tmean_2                  
#> Xcb_tmean_3                  
#> Xcb_tmean_4                  
#> Xcb_tmean_5                  
#> Xcb_tmean_6                  
#> Xcb_rain_1                   
#> Xcb_rain_2                   
#> Xcb_rain_3                   
#> Xcb_rain_4                   
#> Xcb_rain_5                   
#> Xcb_rain_6                   
#> Xcb_wetness_1                
#> Xcb_wetness_2                
#> Xcb_wetness_3                
#> Xcb_wetness_4                
#> Xcb_wetness_5 -0.123         
#> Xcb_wetness_6  0.732  -0.433 
#> 
#> Standardized Within-Group Residuals:
#>         Min          Q1         Med          Q3         Max 
#> -2.81800137 -0.78615982 -0.07432645  0.70191859  5.10243548 
#> 
#> Number of Observations: 520
#> Number of Groups: 20
```

The default summary describes the fixed/population GAM component. The
`component = "lme"` option displays the mixed-model component.

For non-Gaussian responses,
[`mgcv::gamm()`](https://rdrr.io/pkg/mgcv/man/gamm.html) uses a
PQL-based fitting route. Binary responses should therefore be treated
cautiously.

## Spatially autocorrelated effects with `spaMM`

Consider a regional study of a foliar rice disease conducted at
spatially distributed field locations. Daily mean temperature, leaf
wetness duration, and rainfall are recorded during the 86 days preceding
the final assessment. The outcome is the final number of lesions per
plant.

The cross-basis terms estimate non-linear and delayed environmental
effects. A Matérn random effect accounts for residual geographic
similarity among nearby fields after those measured environmental
effects have been considered.

### Spatially dependent fields

``` r

data("st_poisson")

head(st_poisson)
#>   epi_id block time    x_coord    y_coord    tmean  wetness     rain y
#> 1      1   B01    0 -0.4210796 -0.1649767 24.06189 10.81782 0.000000 4
#> 2      1   B01    1 -0.4210796 -0.1649767 24.82615 11.15336 0.000000 4
#> 3      1   B01    2 -0.4210796 -0.1649767 26.11794 12.20862 0.000000 4
#> 4      1   B01    3 -0.4210796 -0.1649767 26.17299 11.70761 6.440981 4
#> 5      1   B01    4 -0.4210796 -0.1649767 26.39594 13.08992 0.000000 4
#> 6      1   B01    5 -0.4210796 -0.1649767 26.23450 14.78350 0.000000 4
```

``` r


cb_spatial <- define_exposures(
  data = st_poisson,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)

dat_spatial <- build_design(
  data = st_poisson,
  cb_templates = cb_spatial,
  groups = c(
    "epi_id",
    "x_coord",
    "y_coord"
  )
)

dat_spatial <- prepare_response(
  data = dat_spatial,
  response = "y",
  family = "poisson"
)
```

``` r

fit_spatial <- fit_epidlnm(
  data = dat_spatial,
  model_engine = "spamm",
  family = "poisson",
  random_effect = NULL,
  spatial_effect = c(
    "x_coord",
    "y_coord"
  ),
  spatial_structure = "matern",
  spatial_group = NULL,
  epiexposure_spec = attr(
    cb_spatial,
    "spec"
  )
)
```

``` r

summary(fit_spatial)
```

This specification generates:

`Matern(1 | x_coord + y_coord)`

The Matérn term is distinct from a conventional random intercept. It
represents autocorrelation that changes continuously with the distance
between locations.

### Independent spatial fields by year

If the disease process is expected to have a different spatial
realization in each growing season, `year` can define independent Matérn
fields:

``` r

data("spatial_year")

head(spatial_year)
#>      epi_id location_id year block block_id time    x_coord  y_coord    tmean
#> 1 L001_2022        L001 2022   B01 2022_B01    0 -0.6698125 0.622697 24.38006
#> 2 L001_2022        L001 2022   B01 2022_B01    1 -0.6698125 0.622697 24.53832
#> 3 L001_2022        L001 2022   B01 2022_B01    2 -0.6698125 0.622697 24.90186
#> 4 L001_2022        L001 2022   B01 2022_B01    3 -0.6698125 0.622697 26.04159
#> 5 L001_2022        L001 2022   B01 2022_B01    4 -0.6698125 0.622697 23.60107
#> 6 L001_2022        L001 2022   B01 2022_B01    5 -0.6698125 0.622697 25.79864
#>     wetness      rain y
#> 1  4.903199 12.904436 3
#> 2  6.781023  7.651138 3
#> 3 11.195548  0.000000 3
#> 4 10.530333  0.000000 3
#> 5 10.322681  0.000000 3
#> 6 12.433350  0.000000 3
```

``` r

dat_spatial_year <- build_design(
  data = spatial_year,
  cb_templates = cb,
  groups = c(
    "epi_id",
    "x_coord",
    "y_coord",
    "year",
    "location_id",
    "block_id"
  )
)

dat_spatial_year <- prepare_response(
  data = dat_spatial_year,
  response = "y",
  family = "poisson"
)
```

``` r

fit_spatial_year <- fit_epidlnm(
  data = dat_spatial_year,
  model_engine = "spamm",
  family = "poisson",
  random_effect = "block_id",
  spatial_effect = c(
    "x_coord",
    "y_coord"
  ),
  spatial_structure = "matern",
  spatial_group = "year",
  epiexposure_spec = attr(
    cb_spatial,
    "spec"
  )
)
```

``` r

summary(fit_spatial_year)
```

The resulting spatial term is:

`Matern(1 | x_coord + y_coord %in% year)`

In this model, `(1 | block_id)` represents conventional experimental
heterogeneity, whereas the Matérn term represents spatially
autocorrelated residual variation within each year.

### Independent spatial fields across epidemiological replicates

The grouping variable in a Matérn term can identify independent
repetitions of an entire epidemiological experiment. Consider a field
study conducted at 130 fixed locations. The same host, pathogen,
experimental layout, and observation period are used in four independent
repetitions.

Each repetition contains an 86-day exposure profile followed by one
final disease assessment at each location. Temperature, rainfall, and
leaf wetness can differ among repetitions, but the temporal positions
remain comparable: the most recent observation corresponds to lag zero,
and the oldest observation corresponds to lag 85.

Each combination of location and repetition represents a separate
epidemic:

`L001_Replicate_1`

`L001_Replicate_2`

`L001_Replicate_3`

`L001_Replicate_4`

The cross-basis terms use the exposure profile from each epidemic to
estimate the shared non-linear and delayed environmental associations.
However, the residual spatial pattern of disease is not required to
remain in the same locations across repetitions.

Independent residual spatial fields can be specified as:

`Matern(1 | x_coord + y_coord %in% epidemic_replicate)`

The model therefore allows each repetition to produce a different
realized map of residual disease intensity. For example, areas with high
residual disease in the first repetition do not need to remain high in
the second repetition.

The Matérn parameters remain shared among repetitions. Consequently, all
repetitions contribute jointly to estimation of the general magnitude,
smoothness, and spatial scale of residual dependence, while retaining
separate realized spatial fields.

``` r

data("spatial_replicate")

head(spatial_replicate)
#>             epi_id location_id epidemic_replicate block        block_id time
#> 1 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    0
#> 2 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    1
#> 3 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    2
#> 4 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    3
#> 5 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    4
#> 6 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    5
#>     x_coord    y_coord    tmean   wetness rain y
#> 1 0.7586476 -0.4480224 23.89557  6.798144    0 5
#> 2 0.7586476 -0.4480224 24.45162  7.355172    0 5
#> 3 0.7586476 -0.4480224 24.49880  8.565241    0 5
#> 4 0.7586476 -0.4480224 24.68348 12.894026    0 5
#> 5 0.7586476 -0.4480224 25.65530 14.958210    0 5
#> 6 0.7586476 -0.4480224 23.98108 13.466749    0 5
```

``` r


cb_replicate <- define_exposures(
  data = spatial_replicate,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)

dat_spatial_replicate <- build_design(
  data = spatial_replicate,
  cb_templates = cb_replicate,
  groups = c(
    "epi_id",
    "x_coord",
    "y_coord",
    "epidemic_replicate",
    "location_id",
    "block_id"
    
  )
)

dat_spatial_replicate <- prepare_response(
  data = dat_spatial_replicate,
  response = "y",
  family = "poisson"
)
```

``` r

fit_spatial_replicate <- fit_epidlnm(
  data = dat_spatial_replicate,
  model_engine = "spamm",
  family = "poisson",
  random_effect = "block_id",
  spatial_effect = c(
    "x_coord",
    "y_coord"
  ),
  spatial_structure = "matern",
  spatial_group = "epidemic_replicate",
  epiexposure_spec = attr(
    cb_replicate,
    "spec"
  )
)
```

``` r

summary(fit_spatial_replicate)
```

The conventional random intercept and the spatial term represent
different sources of dependence. The `block_id` effect accounts for
discrete experimental heterogeneity, whereas the Matérn term accounts
for continuous residual similarity among nearby locations within each
epidemiological repetition.

An epidemiological replicate must represent an independent repetition of
the complete disease process. Repeated assessments made during the
progression of a single epidemic should not be treated as independent
replicates.

## Gaussian generalized least squares

The `gls` engine is restricted to Gaussian outcomes with the identity
link and does not support the `EpiExposure` random-intercept interface.

``` r

data("gaussian_data")

head(gaussian_data)
#>   epi_id block time    tmean  wetness     rain       y
#> 1      1   B01    0 29.16636 11.58879  0.00000 12.4501
#> 2      1   B01    1 27.30294 12.71218  0.00000 12.4501
#> 3      1   B01    2 26.86088 12.74422  0.00000 12.4501
#> 4      1   B01    3 26.29974 11.92657 10.12782 12.4501
#> 5      1   B01    4 26.67772 13.41525  0.00000 12.4501
#> 6      1   B01    5 26.63974 15.59525  0.00000 12.4501
```

``` r

cb_gaussian <- define_exposures(
  data = gaussian_data,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)

dat_gaussian <- build_design(
  data = gaussian_data,
  cb_templates = cb_gaussian,
  groups = c(
    "epi_id"))

dat_gaussian <- prepare_response(
  data = dat_gaussian,
  response = "y",
  family = "gaussian"
)
```

``` r

fit_gls <- fit_epidlnm(
  data = dat_gaussian,
  model_engine = "gls",
  family = "gaussian",
  epiexposure_spec = attr(
    cb_gaussian,
    "spec"
  ),
  method = "REML"
)
```

``` r

summary(fit_gls)
```

**Note:** The INLA and B-DLNM commands in the following sections are
shown but not evaluated during vignette construction because these
models may require substantial computational time. Users can execute the
code interactively after installing the optional packages `INLA`,
`bdlnm`, and `sn`.

## MCMC estimation with `brms`

The `brms` engine provides posterior sampling for the fixed DLNM
coefficients and any requested hierarchical components.

``` r

fit_brms <- fit_epidlnm(
  data = dat_poisson,
  model_engine = "brms",
  family = "poisson",
  random_effect = "block",
  epiexposure_spec = attr(
    cb,
    "spec"
  ),
  chains = 4,
  iter = 2000,
  seed = 123
)
```

``` r

summary(fit_brms, eval=FALSE)
```

## Bayesian DLNM with `bdlnm`

The `bdlnm` engine provides a Bayesian implementation specifically
designed for distributed lag linear and non-linear models. Internally,
the native `bdlnm` package fits the model using Integrated Nested
Laplace Approximation through
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) and then draws
samples from the approximate posterior distribution using
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
【1-fe3db1】

Unlike the direct `inla` engine, which receives the epidemic-level
`cb_*` columns as conventional fixed predictors, the `bdlnm` engine
receives the original named `crossbasis` objects in the model formula.
The fitted object stores the underlying INLA model, the basis objects,
posterior coefficient draws, and summaries of those draws. 【1-fe3db1】

In `EpiExposure`, this engine is integrated into the same prediction and
uncertainty contract used by the other engines. Posterior fixed-effect
draws are retained for draw-by-draw propagation through
[`summarise_effects()`](https://tomazrg.github.io/EpiExposure/reference/summarise_effects.md),
[`predict_outcomes()`](https://tomazrg.github.io/EpiExposure/reference/predict_outcomes.md),
and other downstream functions.

``` r

cb_bdlnm <- define_exposures(
  data = poisson_data,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)


dat_bdlnm <- build_design(
  data = poisson_data,
  cb_templates = cb_bdlnm,
  groups = c(
    "epi_id",
    "block"
  )
)

dat_bdlnm <- prepare_response(
  data = dat_bdlnm,
  response = "y",
  family = "poisson"
)

bdlnm_bases <- attr(
  dat_bdlnm,
  "epiexposure_bdlnm_basis_objects",
  exact = TRUE
)
```

The model can then be fitted through the harmonized `EpiExposure`
interface:

``` r

fit_bdlnm <- fit_epidlnm(
  data = dat_bdlnm,
  model_engine = "bdlnm",
  family = "poisson",
  random_effect = "block",
  random_effect_prior = list(
    prec = list(prior = "pc.prec",param = c(1,0.01))),
  epiexposure_spec = attr(dat_bdlnm, "epiexposure_spec"),
  sample.arg = list(n = 1000, seed = 123))
```

``` r

summary(fit_bdlnm)
```

``` r

fit_bdlnm$model
#fit_bdlnm$coefficients
fit_bdlnm$coefficients.summary


fit_bdlnm$model$summary.random
fit_bdlnm$model$summary.hyperpar

fit_bdlnm$model$summary.random$.epiexposure_re_index

head(fit_bdlnm$model$marginals.hyperpar)
```

The `random_effect = "block"` argument is translated internally to an
independent latent effect of the form

`f(block, model = "iid")`

for INLA-backed engines. When supplied, `random_effect_prior` is passed
directly to the `hyper` argument of the latent effect specification,
allowing custom hyperprior definitions for the random-effect precision.
For example,

`f(block, model = "iid", hyper = random_effect_prior)`

The `sample.arg` list controls posterior coefficient sampling performed
by the native `bdlnm` engine. A nonzero seed provides reproducible
posterior sampling. The native package uses 1,000 posterior samples by
default, although the number of samples can be changed according to the
required precision and available computational resources.

For deterministic `EpiExposure` predictions, posterior means of the
fixed/population coefficients define the central expected response. When
uncertainty is requested, the stored joint posterior coefficient draws
are propagated through the inverse link one draw at a time.

The resulting intervals describe uncertainty in the population-level
expected response. They do not include observation-level noise or
group-specific random effects, even when random effects are included in
the fitted model.

For both `model_engine = "inla"` and `model_engine = "bdlnm"`, random-
effect hyperparameters can be inspected through the underlying INLA fit
(e.g., `summary.hyperpar`), allowing users to evaluate posterior
distributions of random-effect precision parameters.

## Direct INLA versus B-DLNM

Both engines rely on INLA for Bayesian inference and support the same
likelihood families, random effects, and hyperprior specifications. The
main difference lies in how the DLNM basis is represented and what is
returned to the user.

- `model_engine = "inla"` fits the epidemic-level `cb_*` columns
  directly as fixed predictors through
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html). The
  resulting object contains the standard INLA outputs, including
  posterior summaries, marginal distributions, hyperparameter estimates,
  and random-effect summaries.

- `model_engine = "bdlnm"` passes named `crossbasis` objects to the
  native `bdlnm` interface. Internally, the model is fitted through
  INLA, but the original DLNM basis representation is preserved and
  posterior coefficient samples are generated automatically for
  simulation-based uncertainty propagation.

``` r

fit_inla <- fit_epidlnm(
  data = dat_poisson,
  model_engine = "inla",
  family = "poisson",
  random_effect = "block",
  random_effect_prior = list(
      prec = list(
      prior = "pc.prec",
      param = c(1,0.01))),
  epiexposure_spec = attr(cb, "spec"))
```

#### Outputs shared by INLA and B-DLNM

Fixed-effect posterior summaries:

``` r

fit_inla$summary.fixed

fit_bdlnm$model$summary.fixed
```

Random-effect posterior summaries:

``` r

fit_inla$summary.random

fit_bdlnm$model$summary.random
```

Posterior distributions of hyperparameters:

``` r

fit_inla$marginals.hyperpar

fit_bdlnm$model$marginals.hyperpar
```

Marginal likelihood:

``` r

fit_inla$mlik

fit_bdlnm$model$mlik
```

For both engines, random-effect hyperparameters can be inspected through
the INLA summaries. For example:

``` r

fit_inla$summary.hyperpar

fit_bdlnm$model$summary.hyperpar
```

### Outputs specific to B-DLNM

``` r

names(fit_bdlnm)
```

``` r

names(fit_inla)
```

Original DLNM basis objects:

``` r

fit_bdlnm$basis
```

Posterior coefficient draws:

``` r


dim(fit_bdlnm$coefficients)

head(fit_bdlnm$coefficients)
```

Posterior coefficient summaries:

``` r


fit_bdlnm$coefficients.summary
```

#### Model objects

``` r

class(fit_inla)

class(fit_bdlnm)

class(fit_bdlnm$model)
```

The direct INLA engine is useful when users want a conventional INLA
representation of the epidemic-level design matrix and direct access to
standard INLA outputs. The B-DLNM engine is useful when users want to
retain the native DLNM representation and work directly with posterior
coefficient samples, fitted basis objects, and Bayesian DLNM workflows.
