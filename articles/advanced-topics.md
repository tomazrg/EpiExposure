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

poisson_data = readxl::read_xlsx("data/poisson_data.xlsx")
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

data_nb2 = readxl::read_xlsx("data/nb2_data.xlsx")
```

``` r

dat_nb2 <- build_design(
  data = data_nb2,
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

The cross-basis terms estimate nonlinear and delayed weather effects. A
Matérn random effect accounts for residual geographic similarity among
nearby fields after those measured weather effects have been considered.

### Dependence spatial fields

``` r

spatial_data = readxl::read_xlsx("data/spamm_st_poisson.xlsx")
```

``` r


cb_spatial <- define_exposures(
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

dat_spatial <- build_design(
  data = spatial_data,
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
#> formula: y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + 
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + 
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + 
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 + 
#>     Matern(1 | x_coord + y_coord)
#> Estimation of corrPars and lambda by ML (p_v approximation of logL).
#> Estimation of fixed effects by ML (p_v approximation of logL).
#> Estimation of lambda by 'outer' ML, maximizing logL.
#> family: poisson( link = log ) 
#>  ------------ Fixed effects (beta) ------------
#>              Estimate Cond. SE t-value
#> (Intercept)  -3.16606 2.010242 -1.5750
#> cb_tmean_1    0.10014 0.023029  4.3485
#> cb_tmean_2   -0.05371 0.040190 -1.3364
#> cb_tmean_3    0.17495 0.081186  2.1549
#> cb_tmean_4   -0.03074 0.160071 -0.1921
#> cb_tmean_5    0.02089 0.035865  0.5825
#> cb_tmean_6    0.01820 0.074983  0.2427
#> cb_rain_1     0.19496 0.040978  4.7577
#> cb_rain_2    -0.23400 0.048868 -4.7884
#> cb_rain_3     0.08277 0.026439  3.1306
#> cb_rain_4    -0.15720 0.035905 -4.3784
#> cb_rain_5     0.01401 0.061227  0.2288
#> cb_rain_6    -0.06629 0.075269 -0.8808
#> cb_wetness_1  0.06302 0.005466 11.5307
#> cb_wetness_2 -0.11132 0.014985 -7.4290
#> cb_wetness_3  0.03777 0.020279  1.8624
#> cb_wetness_4 -0.14240 0.047853 -2.9757
#> cb_wetness_5  0.04277 0.013694  3.1233
#> cb_wetness_6 -0.09258 0.030615 -3.0241
#>  --------------- Random effects ---------------
#> Family: gaussian( link = identity ) 
#>                    --- Correlation parameters:
#>      1.nu     1.rho 
#> 0.3202374 0.1502092 
#>            --- Variance parameters ('lambda'):
#> lambda = var(u) for u ~ Gaussian; 
#>    x_coord +.  :  0.1207  
#> # of obs: 520; # of groups: x_coord +., 520 
#>  ------------- Likelihood values  -------------
#>                         logLik
#> logL       (p_v(h)): -1498.842
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

spatial_year = readxl::read_xlsx("data/sim_data_spatial_year.xlsx")
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
#> formula: y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + 
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + 
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + 
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 + 
#>     (1 | block_id) + Matern(1 | x_coord + y_coord %in% year)
#> Estimation of corrPars and lambda by ML (p_v approximation of logL).
#> Estimation of fixed effects by ML (p_v approximation of logL).
#> Estimation of lambda by 'outer' ML, maximizing logL.
#> family: poisson( link = log ) 
#>  ------------ Fixed effects (beta) ------------
#>               Estimate Cond. SE t-value
#> (Intercept)  -1.986851 1.543010 -1.2876
#> cb_tmean_1    0.078318 0.018085  4.3307
#> cb_tmean_2   -0.095705 0.030529 -3.1348
#> cb_tmean_3    0.084204 0.062337  1.3508
#> cb_tmean_4   -0.153358 0.124726 -1.2296
#> cb_tmean_5   -0.003494 0.032335 -0.1081
#> cb_tmean_6   -0.045205 0.059690 -0.7573
#> cb_rain_1     0.156362 0.036365  4.2998
#> cb_rain_2    -0.199387 0.048135 -4.1422
#> cb_rain_3     0.082163 0.025598  3.2098
#> cb_rain_4    -0.134819 0.031945 -4.2203
#> cb_rain_5    -0.024302 0.057439 -0.4231
#> cb_rain_6    -0.105148 0.064188 -1.6381
#> cb_wetness_1  0.064872 0.006663  9.7366
#> cb_wetness_2 -0.098978 0.017189 -5.7583
#> cb_wetness_3  0.092460 0.032176  2.8736
#> cb_wetness_4 -0.125711 0.057193 -2.1980
#> cb_wetness_5  0.076935 0.017248  4.4604
#> cb_wetness_6 -0.091845 0.030094 -3.0519
#>  --------------- Random effects ---------------
#> Family: gaussian( link = identity ) 
#>                    --- Correlation parameters:
#>       2.nu      2.rho 
#> 0.39937143 0.04894806 
#>            --- Variance parameters ('lambda'):
#> lambda = var(u) for u ~ Gaussian; 
#>    block_id  :  0.02595 
#>    x_coord +.  :  0.1532  
#> # of obs: 520; # of groups: block_id, 40; x_coord +., 520 
#>  ------------- Likelihood values  -------------
#>                         logLik
#> logL       (p_v(h)): -1516.401
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

Each repetition contains an 86-day exposure history followed by one
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

The cross-basis terms use the exposure history from each epidemic to
estimate the shared nonlinear and delayed environmental associations.
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

spatial_replicate = readxl::read_xlsx("data/sim_data_spatial_replicate.xlsx")
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
#> formula: y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + 
#>     cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + 
#>     cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + 
#>     cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 + 
#>     (1 | block_id) + Matern(1 | x_coord + y_coord %in% epidemic_replicate)
#> Estimation of corrPars and lambda by ML (p_v approximation of logL).
#> Estimation of fixed effects by ML (p_v approximation of logL).
#> Estimation of lambda by 'outer' ML, maximizing logL.
#> family: poisson( link = log ) 
#>  ------------ Fixed effects (beta) ------------
#>               Estimate Cond. SE  t-value
#> (Intercept)  -2.337437 2.179094 -1.07266
#> cb_tmean_1    0.071969 0.023917  3.00910
#> cb_tmean_2   -0.100318 0.041184 -2.43584
#> cb_tmean_3    0.117159 0.092526  1.26623
#> cb_tmean_4   -0.115241 0.158802 -0.72569
#> cb_tmean_5   -0.035818 0.046221 -0.77492
#> cb_tmean_6    0.050521 0.074130  0.68152
#> cb_rain_1     0.188871 0.041777  4.52093
#> cb_rain_2    -0.262496 0.052716 -4.97942
#> cb_rain_3     0.069414 0.026027  2.66702
#> cb_rain_4    -0.109453 0.035869 -3.05148
#> cb_rain_5    -0.022291 0.059042 -0.37754
#> cb_rain_6    -0.006898 0.075571 -0.09128
#> cb_wetness_1  0.070358 0.006307 11.15474
#> cb_wetness_2 -0.085368 0.016292 -5.24002
#> cb_wetness_3  0.060046 0.027462  2.18650
#> cb_wetness_4 -0.058927 0.054821 -1.07491
#> cb_wetness_5  0.063947 0.016237  3.93825
#> cb_wetness_6 -0.099298 0.032569 -3.04886
#>  --------------- Random effects ---------------
#> Family: gaussian( link = identity ) 
#>                    --- Correlation parameters:
#>      2.nu     2.rho 
#> 0.6645319 0.2009695 
#>            --- Variance parameters ('lambda'):
#> lambda = var(u) for u ~ Gaussian; 
#>    block_id  :  0.03897 
#>    x_coord +.  :  0.129  
#> # of obs: 520; # of groups: block_id, 40; x_coord +., 520 
#>  ------------- Likelihood values  -------------
#>                         logLik
#> logL       (p_v(h)): -1572.785
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

gaussian_data = readxl::read_xlsx("data/gaussian_data.xlsx")
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
#> Generalized least squares fit by REML
#>   Model: y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 +      cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 +      cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 +      cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 
#>   Data: dat_gaussian 
#>        AIC     BIC    logLik
#>   2204.988 2289.32 -1082.494
#> 
#> Coefficients:
#>                   Value Std.Error   t-value p-value
#> (Intercept)   2.7834500  5.207711  0.534486  0.5932
#> cb_tmean_1    0.1894659  0.072089  2.628230  0.0088
#> cb_tmean_2   -0.4654679  0.135706 -3.429981  0.0007
#> cb_tmean_3    0.0508730  0.205349  0.247739  0.8044
#> cb_tmean_4   -0.6100407  0.515812 -1.182681  0.2375
#> cb_tmean_5   -0.2300612  0.122648 -1.875786  0.0613
#> cb_tmean_6    0.0891854  0.256516  0.347680  0.7282
#> cb_rain_1     0.5863313  0.161256  3.636023  0.0003
#> cb_rain_2    -0.4958344  0.200520 -2.472745  0.0137
#> cb_rain_3     0.3609263  0.111125  3.247928  0.0012
#> cb_rain_4    -0.5409101  0.142304 -3.801081  0.0002
#> cb_rain_5     0.2930206  0.262242  1.117366  0.2644
#> cb_rain_6    -0.3416745  0.308615 -1.107122  0.2688
#> cb_wetness_1  0.1750801  0.020310  8.620565  0.0000
#> cb_wetness_2 -0.2704239  0.061895 -4.369079  0.0000
#> cb_wetness_3  0.2737794  0.075205  3.640428  0.0003
#> cb_wetness_4 -0.1655927  0.192532 -0.860078  0.3902
#> cb_wetness_5  0.2266449  0.048390  4.683718  0.0000
#> cb_wetness_6 -0.1054772  0.128835 -0.818702  0.4133
#> 
#>  Correlation: 
#>              (Intr) cb_t_1 cb_t_2 cb_t_3 cb_t_4 cb_t_5 cb_t_6 cb_r_1 cb_r_2
#> cb_tmean_1   -0.816                                                        
#> cb_tmean_2   -0.655  0.636                                                 
#> cb_tmean_3   -0.883  0.663  0.547                                          
#> cb_tmean_4   -0.842  0.867  0.782  0.621                                   
#> cb_tmean_5   -0.116 -0.072  0.273  0.461 -0.166                            
#> cb_tmean_6   -0.579  0.519  0.423  0.261  0.793 -0.613                     
#> cb_rain_1     0.039 -0.050 -0.041 -0.020 -0.042 -0.005 -0.009              
#> cb_rain_2     0.051 -0.032 -0.105  0.005 -0.079  0.026 -0.084 -0.048       
#> cb_rain_3    -0.102  0.054  0.058  0.077  0.056  0.048  0.026 -0.666 -0.011
#> cb_rain_4    -0.019  0.003  0.013 -0.017  0.026 -0.045  0.044  0.090 -0.541
#> cb_rain_5    -0.082  0.051  0.042  0.068  0.051  0.037  0.026 -0.605 -0.004
#> cb_rain_6    -0.021  0.005  0.033 -0.020  0.033 -0.035  0.048  0.100 -0.607
#> cb_wetness_1  0.003  0.028 -0.007 -0.042  0.019 -0.032  0.010 -0.036 -0.012
#> cb_wetness_2  0.030 -0.022 -0.075 -0.054 -0.018 -0.092  0.030 -0.047  0.016
#> cb_wetness_3 -0.147 -0.102 -0.032 -0.083 -0.054 -0.043  0.017 -0.024 -0.074
#> cb_wetness_4  0.047 -0.103 -0.071 -0.165 -0.068 -0.083 -0.001  0.004  0.014
#> cb_wetness_5 -0.060 -0.105 -0.051 -0.042 -0.061 -0.021  0.010  0.008 -0.011
#> cb_wetness_6  0.013 -0.054  0.004 -0.070 -0.033 -0.001 -0.022  0.014  0.020
#>              cb_r_3 cb_r_4 cb_r_5 cb_r_6 cb_w_1 cb_w_2 cb_w_3 cb_w_4 cb_w_5
#> cb_tmean_1                                                                 
#> cb_tmean_2                                                                 
#> cb_tmean_3                                                                 
#> cb_tmean_4                                                                 
#> cb_tmean_5                                                                 
#> cb_tmean_6                                                                 
#> cb_rain_1                                                                  
#> cb_rain_2                                                                  
#> cb_rain_3                                                                  
#> cb_rain_4    -0.113                                                        
#> cb_rain_5     0.944 -0.089                                                 
#> cb_rain_6    -0.124  0.862 -0.103                                          
#> cb_wetness_1  0.065 -0.038  0.060 -0.045                                   
#> cb_wetness_2  0.032  0.057  0.001  0.041 -0.208                            
#> cb_wetness_3  0.032  0.061  0.061  0.098 -0.242 -0.033                     
#> cb_wetness_4 -0.001 -0.011 -0.014 -0.004  0.300  0.187  0.045              
#> cb_wetness_5  0.013  0.033  0.046  0.071 -0.569  0.189  0.750 -0.178       
#> cb_wetness_6 -0.018 -0.034 -0.030 -0.052  0.120 -0.146 -0.044  0.722 -0.353
#> 
#> Standardized residuals:
#>          Min           Q1          Med           Q3          Max 
#> -2.492063118 -0.638964653 -0.009302727  0.664222592  3.699250013 
#> 
#> Residual standard error: 1.865893 
#> Degrees of freedom: 520 total; 501 residual
```

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

summary(fit_brms)
#>  Family: poisson 
#>   Links: mu = log 
#> Formula: y_model ~ 1 + cb_tmean_1 + cb_tmean_2 + cb_tmean_3 + cb_tmean_4 + cb_tmean_5 + cb_tmean_6 + cb_rain_1 + cb_rain_2 + cb_rain_3 + cb_rain_4 + cb_rain_5 + cb_rain_6 + cb_wetness_1 + cb_wetness_2 + cb_wetness_3 + cb_wetness_4 + cb_wetness_5 + cb_wetness_6 + (1 | block) 
#>    Data: structure(list(epi_id = c(1, 2, 3, 4, 5, 6, 7, 8,  (Number of observations: 520) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~block (Number of levels: 20) 
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(Intercept)     0.23      0.04     0.16     0.33 1.00      768     1705
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept       -2.61      1.17    -4.89    -0.38 1.00     1338     1906
#> cb_tmean_1       0.08      0.01     0.05     0.11 1.00     1308     1650
#> cb_tmean_2      -0.01      0.03    -0.06     0.04 1.00     1438     2059
#> cb_tmean_3       0.12      0.05     0.03     0.20 1.00     1441     2130
#> cb_tmean_4       0.17      0.10    -0.03     0.37 1.00     1238     1535
#> cb_tmean_5      -0.02      0.02    -0.06     0.01 1.00     2616     2758
#> cb_tmean_6       0.12      0.04     0.04     0.20 1.00     1696     2045
#> cb_rain_1        0.13      0.03     0.08     0.18 1.00     2616     2529
#> cb_rain_2       -0.24      0.03    -0.31    -0.18 1.00     2235     2867
#> cb_rain_3        0.06      0.02     0.02     0.09 1.00     2159     2658
#> cb_rain_4       -0.06      0.02    -0.10    -0.01 1.00     2043     2806
#> cb_rain_5       -0.04      0.04    -0.11     0.03 1.00     2152     2500
#> cb_rain_6        0.06      0.05    -0.04     0.16 1.00     1934     2283
#> cb_wetness_1     0.05      0.00     0.05     0.06 1.00     3395     2625
#> cb_wetness_2    -0.09      0.01    -0.12    -0.07 1.00     3159     3047
#> cb_wetness_3     0.05      0.01     0.03     0.08 1.00     2618     2816
#> cb_wetness_4    -0.09      0.04    -0.17    -0.02 1.00     2197     2263
#> cb_wetness_5     0.06      0.01     0.04     0.07 1.00     2180     2550
#> cb_wetness_6    -0.09      0.02    -0.13    -0.05 1.00     2345     2649
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

## Bayesian DLNM with `bdlnm`

The `bdlnm` engine provides a Bayesian implementation specifically
designed for distributed lag linear and nonlinear models. Internally,
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
#>                      Length Class  Mode   
#> model                   53  inla   list   
#> basis                    3  -none- list   
#> coefficients         19000  -none- numeric
#> coefficients.summary   114  -none- numeric
```

``` r

fit_bdlnm$model
#> Time used:
#>   Pre = 1.78, Running = 1.75, Post = 0.52, Total = 4.05
#fit_bdlnm$coefficients
fit_bdlnm$coefficients.summary
#>                       mean          sd  0.025quant     0.5quant   0.975quant
#> (Intercept)   -2.612308666 1.181011217 -4.86173293 -2.659073564 -0.271741415
#> tmean_v1.l1    0.081827689 0.014285124  0.05468580  0.081841641  0.108500115
#> tmean_v1.l2   -0.008881804 0.026473706 -0.05979620 -0.009364572  0.044742272
#> tmean_v2.l1    0.116775335 0.046154672  0.02792660  0.118181920  0.209262892
#> tmean_v2.l2    0.170686691 0.097874377 -0.01795049  0.173623565  0.358115565
#> tmean_v3.l1   -0.023144634 0.017116497 -0.05412196 -0.022838542  0.008978175
#> tmean_v3.l2    0.119642910 0.038671644  0.04575193  0.120948587  0.197361979
#> rain_v1.l1     0.128096901 0.026269755  0.07581237  0.129003281  0.178838059
#> rain_v1.l2    -0.239686348 0.033790575 -0.30519545 -0.239825436 -0.169897806
#> rain_v2.l1     0.056297062 0.016763660  0.02603588  0.055803307  0.090203720
#> rain_v2.l2    -0.058826311 0.024053725 -0.10389328 -0.059182823 -0.012161076
#> rain_v3.l1    -0.038708686 0.038858203 -0.11258207 -0.038844164  0.038305274
#> rain_v3.l2     0.058778348 0.053433229 -0.04509105  0.056894884  0.164445332
#> wetness_v1.l1  0.052609959 0.003723103  0.04538249  0.052604149  0.059713221
#> wetness_v1.l2 -0.093896682 0.011069018 -0.11651535 -0.093612535 -0.073576562
#> wetness_v2.l1  0.053719862 0.014373251  0.02434727  0.053718821  0.082410374
#> wetness_v2.l2 -0.093244212 0.038430918 -0.16970272 -0.093340193 -0.017618942
#> wetness_v3.l1  0.055673487 0.009673492  0.03697027  0.055688038  0.073729222
#> wetness_v3.l2 -0.090573738 0.018771943 -0.12753240 -0.090908305 -0.056623194
#>                      mode
#> (Intercept)   -2.76660074
#> tmean_v1.l1    0.08230571
#> tmean_v1.l2   -0.01506542
#> tmean_v2.l1    0.12123077
#> tmean_v2.l2    0.17477959
#> tmean_v3.l1   -0.02269479
#> tmean_v3.l2    0.12843803
#> rain_v1.l1     0.13628989
#> rain_v1.l2    -0.24148249
#> rain_v2.l1     0.05473600
#> rain_v2.l2    -0.05831061
#> rain_v3.l1    -0.03565000
#> rain_v3.l2     0.03976396
#> wetness_v1.l1  0.05214164
#> wetness_v1.l2 -0.09609658
#> wetness_v2.l1  0.05363214
#> wetness_v2.l2 -0.10424887
#> wetness_v3.l1  0.05508213
#> wetness_v3.l2 -0.09485818


fit_bdlnm$model$summary.random
#> $.epiexposure_re_index
#>    ID         mean         sd   0.025quant     0.5quant  0.975quant
#> 1   1 -0.073388600 0.07401964 -0.219230268 -0.073289776  0.07188893
#> 2   2  0.428869338 0.06779024  0.297485255  0.428177628  0.56421789
#> 3   3 -0.458404662 0.08169741 -0.621827940 -0.457389857 -0.30078336
#> 4   4  0.314021054 0.06715646  0.183353933  0.313502297  0.44765874
#> 5   5  0.121811376 0.07259087 -0.020045878  0.121497678  0.26546539
#> 6   6  0.007532425 0.07655867 -0.142698675  0.007433813  0.15832881
#> 7   7 -0.039792085 0.07263410 -0.182653273 -0.039789699  0.10305691
#> 8   8  0.159194341 0.07213827  0.018453096  0.158800956  0.30218902
#> 9   9 -0.187911275 0.07984201 -0.345930061 -0.187531395 -0.03206174
#> 10 10 -0.102041861 0.08013782 -0.260104513 -0.101850109  0.05492633
#> 11 11 -0.319360601 0.07959719 -0.477726471 -0.318687590 -0.16483988
#> 12 12  0.148699710 0.07219596  0.007838944  0.148308530  0.29180133
#> 13 13 -0.070269438 0.07550691 -0.218992968 -0.070178362  0.07793448
#> 14 14 -0.216943963 0.07863385 -0.372770103 -0.216504362 -0.06362880
#> 15 15  0.067999068 0.07306866 -0.075133978  0.067802882  0.21225578
#> 16 16  0.083026631 0.07243394 -0.058666284  0.082760793  0.22624234
#> 17 17  0.080406530 0.07133697 -0.059314168  0.080198346  0.22131952
#> 18 18  0.071200642 0.07351787 -0.072888772  0.071032322  0.21625386
#> 19 19 -0.093189357 0.07683668 -0.244633056 -0.093055729  0.05749175
#> 20 20  0.078548921 0.07661603 -0.071499291  0.078349588  0.22973814
#>            mode          kld
#> 1  -0.073291040 1.653910e-08
#> 2   0.428191219 3.488174e-08
#> 3  -0.457401128 4.078086e-08
#> 4   0.313512456 2.913033e-08
#> 5   0.121503415 1.968406e-08
#> 6   0.007435662 1.446849e-08
#> 7  -0.039789324 1.727782e-08
#> 8   0.158808209 2.165640e-08
#> 9  -0.187535834 1.667722e-08
#> 10 -0.101852209 1.337840e-08
#> 11 -0.318695621 2.573571e-08
#> 12  0.148315763 2.158594e-08
#> 13 -0.070179398 1.539789e-08
#> 14 -0.216509789 1.888167e-08
#> 15  0.067806478 1.772473e-08
#> 16  0.082765763 1.916018e-08
#> 17  0.080202269 1.928677e-08
#> 18  0.071035304 1.703447e-08
#> 19 -0.093057267 1.483206e-08
#> 20  0.078352811 1.500913e-08
fit_bdlnm$model$summary.hyperpar
#>                                         mean       sd 0.025quant 0.5quant
#> Precision for .epiexposure_re_index 22.38976 7.788687   10.19738  21.3789
#>                                     0.975quant     mode
#> Precision for .epiexposure_re_index   40.41147 19.41154

fit_bdlnm$model$summary.random$.epiexposure_re_index
#>    ID         mean         sd   0.025quant     0.5quant  0.975quant
#> 1   1 -0.073388600 0.07401964 -0.219230268 -0.073289776  0.07188893
#> 2   2  0.428869338 0.06779024  0.297485255  0.428177628  0.56421789
#> 3   3 -0.458404662 0.08169741 -0.621827940 -0.457389857 -0.30078336
#> 4   4  0.314021054 0.06715646  0.183353933  0.313502297  0.44765874
#> 5   5  0.121811376 0.07259087 -0.020045878  0.121497678  0.26546539
#> 6   6  0.007532425 0.07655867 -0.142698675  0.007433813  0.15832881
#> 7   7 -0.039792085 0.07263410 -0.182653273 -0.039789699  0.10305691
#> 8   8  0.159194341 0.07213827  0.018453096  0.158800956  0.30218902
#> 9   9 -0.187911275 0.07984201 -0.345930061 -0.187531395 -0.03206174
#> 10 10 -0.102041861 0.08013782 -0.260104513 -0.101850109  0.05492633
#> 11 11 -0.319360601 0.07959719 -0.477726471 -0.318687590 -0.16483988
#> 12 12  0.148699710 0.07219596  0.007838944  0.148308530  0.29180133
#> 13 13 -0.070269438 0.07550691 -0.218992968 -0.070178362  0.07793448
#> 14 14 -0.216943963 0.07863385 -0.372770103 -0.216504362 -0.06362880
#> 15 15  0.067999068 0.07306866 -0.075133978  0.067802882  0.21225578
#> 16 16  0.083026631 0.07243394 -0.058666284  0.082760793  0.22624234
#> 17 17  0.080406530 0.07133697 -0.059314168  0.080198346  0.22131952
#> 18 18  0.071200642 0.07351787 -0.072888772  0.071032322  0.21625386
#> 19 19 -0.093189357 0.07683668 -0.244633056 -0.093055729  0.05749175
#> 20 20  0.078548921 0.07661603 -0.071499291  0.078349588  0.22973814
#>            mode          kld
#> 1  -0.073291040 1.653910e-08
#> 2   0.428191219 3.488174e-08
#> 3  -0.457401128 4.078086e-08
#> 4   0.313512456 2.913033e-08
#> 5   0.121503415 1.968406e-08
#> 6   0.007435662 1.446849e-08
#> 7  -0.039789324 1.727782e-08
#> 8   0.158808209 2.165640e-08
#> 9  -0.187535834 1.667722e-08
#> 10 -0.101852209 1.337840e-08
#> 11 -0.318695621 2.573571e-08
#> 12  0.148315763 2.158594e-08
#> 13 -0.070179398 1.539789e-08
#> 14 -0.216509789 1.888167e-08
#> 15  0.067806478 1.772473e-08
#> 16  0.082765763 1.916018e-08
#> 17  0.080202269 1.928677e-08
#> 18  0.071035304 1.703447e-08
#> 19 -0.093057267 1.483206e-08
#> 20  0.078352811 1.500913e-08

head(fit_bdlnm$model$marginals.hyperpar)
#> $`Precision for .epiexposure_re_index`
#>               x            y
#>  [1,]  5.262877 3.677687e-04
#>  [2,]  5.563902 6.037148e-04
#>  [3,]  6.443650 1.508954e-03
#>  [4,]  8.756236 7.228568e-03
#>  [5,] 10.197377 1.400626e-02
#>  [6,] 11.580572 2.239259e-02
#>  [7,] 13.354980 3.402146e-02
#>  [8,] 14.669751 4.197071e-02
#>  [9,] 15.785451 4.761113e-02
#> [10,] 16.794004 5.152970e-02
#> [11,] 17.741349 5.405639e-02
#> [12,] 18.655367 5.540598e-02
#> [13,] 19.106387 5.568861e-02
#> [14,] 19.555811 5.573030e-02
#> [15,] 20.006315 5.554376e-02
#> [16,] 20.458859 5.514031e-02
#> [17,] 20.640904 5.492043e-02
#> [18,] 20.823872 5.466775e-02
#> [19,] 20.915721 5.452929e-02
#> [20,] 21.007821 5.438283e-02
#> [21,] 21.192811 5.406623e-02
#> [22,] 21.378899 5.371853e-02
#> [23,] 21.566148 5.334032e-02
#> [24,] 21.754633 5.293218e-02
#> [25,] 21.849458 5.271685e-02
#> [26,] 21.944717 5.249402e-02
#> [27,] 22.136583 5.202605e-02
#> [28,] 22.330315 5.152871e-02
#> [29,] 22.823384 5.015987e-02
#> [30,] 23.330497 4.861673e-02
#> [31,] 23.855472 4.689995e-02
#> [32,] 24.400909 4.501459e-02
#> [33,] 25.569635 4.074642e-02
#> [34,] 26.879116 3.582296e-02
#> [35,] 28.398234 3.023788e-02
#> [36,] 30.252846 2.396381e-02
#> [37,] 32.720973 1.694342e-02
#> [38,] 36.675252 9.074608e-03
#> [39,] 40.411465 4.754373e-03
#> [40,] 45.123917 1.992539e-03
#> [41,] 56.251531 2.131189e-04
#> [42,] 66.783545 2.297470e-05
#> [43,] 74.596908 4.545596e-06
#> attr(,"hyperid")
#> [1] "1001|.epiexposure_re_index"
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
#>                      mean          sd  0.025quant     0.5quant   0.975quant
#> (Intercept)  -2.591106208 1.158867213 -4.86350062 -2.591116679 -0.318652173
#> cb_tmean_1    0.081789429 0.014101386  0.05413782  0.081789494  0.109440665
#> cb_tmean_2   -0.008803491 0.025844186 -0.05948172 -0.008803368  0.041874047
#> cb_tmean_3    0.115609698 0.044917339  0.02753051  0.115609957  0.203687417
#> cb_tmean_4    0.170841534 0.098491876 -0.02229254  0.170841979  0.363973067
#> cb_tmean_5   -0.023434331 0.017048938 -0.05686528 -0.023434454  0.009997315
#> cb_tmean_6    0.119556167 0.039905670  0.04130404  0.119556540  0.197806167
#> cb_rain_1     0.129278390 0.026463744  0.07738752  0.129277719  0.181173080
#> cb_rain_2    -0.239973194 0.033977749 -0.30659970 -0.239973368 -0.173345697
#> cb_rain_3     0.055882237 0.016277830  0.02396230  0.055882493  0.087800719
#> cb_rain_4    -0.057647986 0.023724497 -0.10416900 -0.057648122 -0.011126200
#> cb_rain_5    -0.038490568 0.037558099 -0.11213916 -0.038490257  0.035156255
#> cb_rain_6     0.060488279 0.052021394 -0.04151905  0.060487752  0.162498607
#> cb_wetness_1  0.052652397 0.003728930  0.04534043  0.052652362  0.059964559
#> cb_wetness_2 -0.093616653 0.011114282 -0.11541091 -0.093616566 -0.071822893
#> cb_wetness_3  0.054014970 0.014580718  0.02542389  0.054014874  0.082606594
#> cb_wetness_4 -0.093549426 0.037424146 -0.16693418 -0.093549527 -0.020164091
#> cb_wetness_5  0.055838400 0.009722635  0.03677308  0.055838469  0.074903329
#> cb_wetness_6 -0.090635888 0.018671733 -0.12724937 -0.090635867 -0.054022531
#>                      mode          kld
#> (Intercept)  -2.591116730 5.510611e-11
#> cb_tmean_1    0.081789494 5.522464e-11
#> cb_tmean_2   -0.008803368 5.526328e-11
#> cb_tmean_3    0.115609958 5.520280e-11
#> cb_tmean_4    0.170841981 5.527258e-11
#> cb_tmean_5   -0.023434454 5.487164e-11
#> cb_tmean_6    0.119556542 5.526886e-11
#> cb_rain_1     0.129277715 5.536695e-11
#> cb_rain_2    -0.239973369 5.523186e-11
#> cb_rain_3     0.055882495 5.533254e-11
#> cb_rain_4    -0.057648122 5.501625e-11
#> cb_rain_5    -0.038490256 5.527607e-11
#> cb_rain_6     0.060487749 5.504803e-11
#> cb_wetness_1  0.052652362 5.518856e-11
#> cb_wetness_2 -0.093616565 5.523869e-11
#> cb_wetness_3  0.054014874 5.506714e-11
#> cb_wetness_4 -0.093549528 5.519009e-11
#> cb_wetness_5  0.055838470 5.514581e-11
#> cb_wetness_6 -0.090635867 5.524441e-11

fit_bdlnm$model$summary.fixed
#>                      mean          sd  0.025quant     0.5quant   0.975quant
#> (Intercept)   -2.59110621 1.158867216 -4.86350063 -2.591116683 -0.318652172
#> tmean_v1.l1    0.08178943 0.014101386  0.05413782  0.081789494  0.109440665
#> tmean_v1.l2   -0.00880349 0.025844186 -0.05948172 -0.008803368  0.041874047
#> tmean_v2.l1    0.11560970 0.044917339  0.02753051  0.115609957  0.203687417
#> tmean_v2.l2    0.17084153 0.098491876 -0.02229254  0.170841980  0.363973068
#> tmean_v3.l1   -0.02343433 0.017048938 -0.05686528 -0.023434454  0.009997315
#> tmean_v3.l2    0.11955617 0.039905670  0.04130404  0.119556541  0.197806167
#> rain_v1.l1     0.12927839 0.026463744  0.07738752  0.129277718  0.181173080
#> rain_v1.l2    -0.23997319 0.033977749 -0.30659970 -0.239973368 -0.173345697
#> rain_v2.l1     0.05588224 0.016277830  0.02396230  0.055882494  0.087800719
#> rain_v2.l2    -0.05764799 0.023724497 -0.10416900 -0.057648122 -0.011126200
#> rain_v3.l1    -0.03849057 0.037558099 -0.11213916 -0.038490257  0.035156256
#> rain_v3.l2     0.06048828 0.052021395 -0.04151905  0.060487752  0.162498607
#> wetness_v1.l1  0.05265240 0.003728930  0.04534043  0.052652362  0.059964559
#> wetness_v1.l2 -0.09361665 0.011114282 -0.11541091 -0.093616566 -0.071822893
#> wetness_v2.l1  0.05401497 0.014580718  0.02542389  0.054014874  0.082606594
#> wetness_v2.l2 -0.09354943 0.037424146 -0.16693418 -0.093549528 -0.020164092
#> wetness_v3.l1  0.05583840 0.009722635  0.03677308  0.055838469  0.074903329
#> wetness_v3.l2 -0.09063589 0.018671733 -0.12724937 -0.090635867 -0.054022531
#>                       mode          kld
#> (Intercept)   -2.591116734 5.510561e-11
#> tmean_v1.l1    0.081789494 5.522355e-11
#> tmean_v1.l2   -0.008803368 5.526328e-11
#> tmean_v2.l1    0.115609958 5.520280e-11
#> tmean_v2.l2    0.170841982 5.527276e-11
#> tmean_v3.l1   -0.023434454 5.487155e-11
#> tmean_v3.l2    0.119556542 5.526940e-11
#> rain_v1.l1     0.129277715 5.536571e-11
#> rain_v1.l2    -0.239973369 5.523487e-11
#> rain_v2.l1     0.055882495 5.533336e-11
#> rain_v2.l2    -0.057648122 5.501702e-11
#> rain_v3.l1    -0.038490256 5.527607e-11
#> rain_v3.l2     0.060487749 5.504811e-11
#> wetness_v1.l1  0.052652362 5.517297e-11
#> wetness_v1.l2 -0.093616565 5.523869e-11
#> wetness_v2.l1  0.054014874 5.506714e-11
#> wetness_v2.l2 -0.093549528 5.519040e-11
#> wetness_v3.l1  0.055838470 5.514811e-11
#> wetness_v3.l2 -0.090635867 5.524068e-11
```

Random-effect posterior summaries:

``` r

fit_inla$summary.random
#> $.epiexposure_re_index
#>    ID         mean         sd   0.025quant     0.5quant  0.975quant
#> 1   1 -0.073388600 0.07401964 -0.219230271 -0.073289776  0.07188894
#> 2   2  0.428869340 0.06779024  0.297485255  0.428177630  0.56421790
#> 3   3 -0.458404665 0.08169742 -0.621827945 -0.457389861 -0.30078336
#> 4   4  0.314021055 0.06715646  0.183353931  0.313502298  0.44765874
#> 5   5  0.121811377 0.07259087 -0.020045880  0.121497679  0.26546539
#> 6   6  0.007532425 0.07655867 -0.142698677  0.007433813  0.15832881
#> 7   7 -0.039792085 0.07263410 -0.182653275 -0.039789699  0.10305691
#> 8   8  0.159194342 0.07213828  0.018453095  0.158800957  0.30218903
#> 9   9 -0.187911276 0.07984201 -0.345930064 -0.187531397 -0.03206174
#> 10 10 -0.102041861 0.08013782 -0.260104516 -0.101850110  0.05492633
#> 11 11 -0.319360603 0.07959719 -0.477726475 -0.318687592 -0.16483988
#> 12 12  0.148699711 0.07219597  0.007838943  0.148308532  0.29180134
#> 13 13 -0.070269439 0.07550691 -0.218992970 -0.070178362  0.07793448
#> 14 14 -0.216943964 0.07863386 -0.372770107 -0.216504364 -0.06362880
#> 15 15  0.067999069 0.07306866 -0.075133980  0.067802883  0.21225578
#> 16 16  0.083026632 0.07243394 -0.058666286  0.082760794  0.22624234
#> 17 17  0.080406531 0.07133697 -0.059314170  0.080198347  0.22131952
#> 18 18  0.071200642 0.07351787 -0.072888775  0.071032322  0.21625386
#> 19 19 -0.093189357 0.07683668 -0.244633059 -0.093055730  0.05749175
#> 20 20  0.078548921 0.07661603 -0.071499293  0.078349588  0.22973814
#>            mode          kld
#> 1  -0.073291040 1.653910e-08
#> 2   0.428191221 3.488175e-08
#> 3  -0.457401131 4.078085e-08
#> 4   0.313512458 2.913034e-08
#> 5   0.121503416 1.968406e-08
#> 6   0.007435663 1.446849e-08
#> 7  -0.039789324 1.727783e-08
#> 8   0.158808210 2.165640e-08
#> 9  -0.187535835 1.667723e-08
#> 10 -0.101852209 1.337840e-08
#> 11 -0.318695623 2.573571e-08
#> 12  0.148315764 2.158594e-08
#> 13 -0.070179398 1.539789e-08
#> 14 -0.216509790 1.888167e-08
#> 15  0.067806479 1.772473e-08
#> 16  0.082765764 1.916019e-08
#> 17  0.080202270 1.928677e-08
#> 18  0.071035305 1.703448e-08
#> 19 -0.093057268 1.483207e-08
#> 20  0.078352812 1.500913e-08

fit_bdlnm$model$summary.random
#> $.epiexposure_re_index
#>    ID         mean         sd   0.025quant     0.5quant  0.975quant
#> 1   1 -0.073388600 0.07401964 -0.219230268 -0.073289776  0.07188893
#> 2   2  0.428869338 0.06779024  0.297485255  0.428177628  0.56421789
#> 3   3 -0.458404662 0.08169741 -0.621827940 -0.457389857 -0.30078336
#> 4   4  0.314021054 0.06715646  0.183353933  0.313502297  0.44765874
#> 5   5  0.121811376 0.07259087 -0.020045878  0.121497678  0.26546539
#> 6   6  0.007532425 0.07655867 -0.142698675  0.007433813  0.15832881
#> 7   7 -0.039792085 0.07263410 -0.182653273 -0.039789699  0.10305691
#> 8   8  0.159194341 0.07213827  0.018453096  0.158800956  0.30218902
#> 9   9 -0.187911275 0.07984201 -0.345930061 -0.187531395 -0.03206174
#> 10 10 -0.102041861 0.08013782 -0.260104513 -0.101850109  0.05492633
#> 11 11 -0.319360601 0.07959719 -0.477726471 -0.318687590 -0.16483988
#> 12 12  0.148699710 0.07219596  0.007838944  0.148308530  0.29180133
#> 13 13 -0.070269438 0.07550691 -0.218992968 -0.070178362  0.07793448
#> 14 14 -0.216943963 0.07863385 -0.372770103 -0.216504362 -0.06362880
#> 15 15  0.067999068 0.07306866 -0.075133978  0.067802882  0.21225578
#> 16 16  0.083026631 0.07243394 -0.058666284  0.082760793  0.22624234
#> 17 17  0.080406530 0.07133697 -0.059314168  0.080198346  0.22131952
#> 18 18  0.071200642 0.07351787 -0.072888772  0.071032322  0.21625386
#> 19 19 -0.093189357 0.07683668 -0.244633056 -0.093055729  0.05749175
#> 20 20  0.078548921 0.07661603 -0.071499291  0.078349588  0.22973814
#>            mode          kld
#> 1  -0.073291040 1.653910e-08
#> 2   0.428191219 3.488174e-08
#> 3  -0.457401128 4.078086e-08
#> 4   0.313512456 2.913033e-08
#> 5   0.121503415 1.968406e-08
#> 6   0.007435662 1.446849e-08
#> 7  -0.039789324 1.727782e-08
#> 8   0.158808209 2.165640e-08
#> 9  -0.187535834 1.667722e-08
#> 10 -0.101852209 1.337840e-08
#> 11 -0.318695621 2.573571e-08
#> 12  0.148315763 2.158594e-08
#> 13 -0.070179398 1.539789e-08
#> 14 -0.216509789 1.888167e-08
#> 15  0.067806478 1.772473e-08
#> 16  0.082765763 1.916018e-08
#> 17  0.080202269 1.928677e-08
#> 18  0.071035304 1.703447e-08
#> 19 -0.093057267 1.483206e-08
#> 20  0.078352811 1.500913e-08
```

Posterior distributions of hyperparameters:

``` r

fit_inla$marginals.hyperpar
#> $`Precision for .epiexposure_re_index`
#>               x            y
#>  [1,]  5.262894 3.677747e-04
#>  [2,]  5.563919 6.037264e-04
#>  [3,]  6.443657 1.508963e-03
#>  [4,]  8.756237 7.228572e-03
#>  [5,] 10.197377 1.400626e-02
#>  [6,] 11.580572 2.239259e-02
#>  [7,] 13.354980 3.402146e-02
#>  [8,] 14.669752 4.197071e-02
#>  [9,] 15.785451 4.761112e-02
#> [10,] 16.794004 5.152969e-02
#> [11,] 17.741349 5.405639e-02
#> [12,] 18.655368 5.540598e-02
#> [13,] 19.106388 5.568861e-02
#> [14,] 19.555812 5.573030e-02
#> [15,] 20.006315 5.554375e-02
#> [16,] 20.458860 5.514030e-02
#> [17,] 20.640905 5.492042e-02
#> [18,] 20.823873 5.466775e-02
#> [19,] 20.915721 5.452928e-02
#> [20,] 21.007822 5.438283e-02
#> [21,] 21.192811 5.406622e-02
#> [22,] 21.378900 5.371852e-02
#> [23,] 21.566148 5.334032e-02
#> [24,] 21.754634 5.293218e-02
#> [25,] 21.849458 5.271685e-02
#> [26,] 21.944718 5.249402e-02
#> [27,] 22.136584 5.202605e-02
#> [28,] 22.330316 5.152871e-02
#> [29,] 22.823385 5.015987e-02
#> [30,] 23.330498 4.861674e-02
#> [31,] 23.855473 4.689995e-02
#> [32,] 24.400910 4.501459e-02
#> [33,] 25.569635 4.074643e-02
#> [34,] 26.879117 3.582297e-02
#> [35,] 28.398234 3.023789e-02
#> [36,] 30.252846 2.396381e-02
#> [37,] 32.720973 1.694342e-02
#> [38,] 36.675252 9.074608e-03
#> [39,] 40.411465 4.754374e-03
#> [40,] 45.123917 1.992539e-03
#> [41,] 56.251528 2.131189e-04
#> [42,] 66.783527 2.297480e-05
#> [43,] 74.596781 4.545721e-06
#> attr(,"hyperid")
#> [1] "1001|.epiexposure_re_index"

fit_bdlnm$model$marginals.hyperpar
#> $`Precision for .epiexposure_re_index`
#>               x            y
#>  [1,]  5.262877 3.677687e-04
#>  [2,]  5.563902 6.037148e-04
#>  [3,]  6.443650 1.508954e-03
#>  [4,]  8.756236 7.228568e-03
#>  [5,] 10.197377 1.400626e-02
#>  [6,] 11.580572 2.239259e-02
#>  [7,] 13.354980 3.402146e-02
#>  [8,] 14.669751 4.197071e-02
#>  [9,] 15.785451 4.761113e-02
#> [10,] 16.794004 5.152970e-02
#> [11,] 17.741349 5.405639e-02
#> [12,] 18.655367 5.540598e-02
#> [13,] 19.106387 5.568861e-02
#> [14,] 19.555811 5.573030e-02
#> [15,] 20.006315 5.554376e-02
#> [16,] 20.458859 5.514031e-02
#> [17,] 20.640904 5.492043e-02
#> [18,] 20.823872 5.466775e-02
#> [19,] 20.915721 5.452929e-02
#> [20,] 21.007821 5.438283e-02
#> [21,] 21.192811 5.406623e-02
#> [22,] 21.378899 5.371853e-02
#> [23,] 21.566148 5.334032e-02
#> [24,] 21.754633 5.293218e-02
#> [25,] 21.849458 5.271685e-02
#> [26,] 21.944717 5.249402e-02
#> [27,] 22.136583 5.202605e-02
#> [28,] 22.330315 5.152871e-02
#> [29,] 22.823384 5.015987e-02
#> [30,] 23.330497 4.861673e-02
#> [31,] 23.855472 4.689995e-02
#> [32,] 24.400909 4.501459e-02
#> [33,] 25.569635 4.074642e-02
#> [34,] 26.879116 3.582296e-02
#> [35,] 28.398234 3.023788e-02
#> [36,] 30.252846 2.396381e-02
#> [37,] 32.720973 1.694342e-02
#> [38,] 36.675252 9.074608e-03
#> [39,] 40.411465 4.754373e-03
#> [40,] 45.123917 1.992539e-03
#> [41,] 56.251531 2.131189e-04
#> [42,] 66.783545 2.297470e-05
#> [43,] 74.596908 4.545596e-06
#> attr(,"hyperid")
#> [1] "1001|.epiexposure_re_index"
```

Marginal likelihood:

``` r

fit_inla$mlik
#>                                            [,1]
#> log marginal-likelihood (integration) -1583.137
#> log marginal-likelihood (Gaussian)    -1583.524

fit_bdlnm$model$mlik
#>                                            [,1]
#> log marginal-likelihood (integration) -1583.137
#> log marginal-likelihood (Gaussian)    -1583.524
```

For both engines, random-effect hyperparameters can be inspected through
the INLA summaries. For example:

``` r

fit_inla$summary.hyperpar
#>                                         mean       sd 0.025quant 0.5quant
#> Precision for .epiexposure_re_index 22.38976 7.788686   10.19738  21.3789
#>                                     0.975quant     mode
#> Precision for .epiexposure_re_index   40.41147 19.41154

fit_bdlnm$model$summary.hyperpar
#>                                         mean       sd 0.025quant 0.5quant
#> Precision for .epiexposure_re_index 22.38976 7.788687   10.19738  21.3789
#>                                     0.975quant     mode
#> Precision for .epiexposure_re_index   40.41147 19.41154
```

### Outputs specific to B-DLNM

``` r

names(fit_bdlnm)
#> [1] "model"                "basis"                "coefficients"        
#> [4] "coefficients.summary"
```

``` r

names(fit_inla)
#>  [1] "names.fixed"                 "summary.fixed"              
#>  [3] "marginals.fixed"             "summary.lincomb"            
#>  [5] "marginals.lincomb"           "size.lincomb"               
#>  [7] "summary.lincomb.derived"     "marginals.lincomb.derived"  
#>  [9] "size.lincomb.derived"        "mlik"                       
#> [11] "cpo"                         "gcpo"                       
#> [13] "po"                          "waic"                       
#> [15] "residuals"                   "model.random"               
#> [17] "summary.random"              "marginals.random"           
#> [19] "size.random"                 "summary.linear.predictor"   
#> [21] "marginals.linear.predictor"  "summary.fitted.values"      
#> [23] "marginals.fitted.values"     "size.linear.predictor"      
#> [25] "summary.hyperpar"            "marginals.hyperpar"         
#> [27] "internal.summary.hyperpar"   "internal.marginals.hyperpar"
#> [29] "offset.linear.predictor"     "model.spde2.blc"            
#> [31] "summary.spde2.blc"           "marginals.spde2.blc"        
#> [33] "size.spde2.blc"              "model.spde3.blc"            
#> [35] "summary.spde3.blc"           "marginals.spde3.blc"        
#> [37] "size.spde3.blc"              "logfile"                    
#> [39] "misc"                        "dic"                        
#> [41] "mode"                        "joint.hyper"                
#> [43] "nhyper"                      "version"                    
#> [45] "Q"                           "graph"                      
#> [47] "ok"                          "cpu.intern"                 
#> [49] "cpu.used"                    "all.hyper"                  
#> [51] ".args"                       "call"                       
#> [53] "model.matrix"
```

Original DLNM basis objects:

``` r

fit_bdlnm$basis
#> $tmean
#>        tmean_v1.l1 tmean_v1.l2 tmean_v2.l1  tmean_v2.l2  tmean_v3.l1
#>   [1,]   9.2720803   6.0999038    20.69975 -0.296999198  -8.55582100
#>   [2,]   7.3149810   3.8839791    22.12025  0.657264200 -11.49886510
#>   [3,]   5.3143838   4.9881901    22.01426  1.143664859 -12.13531743
#>   [4,]   5.5176082   4.3670748    22.34706  1.180775396 -11.80361340
#>   [5,]  16.6086278  -0.2180438    16.12271  0.490455892   4.27917021
#>   [6,]  13.8360139   4.3640033    19.19509 -0.506888326  -5.71928105
#>   [7,]  14.4294978   4.3417052    19.84170 -0.046449341  -8.12076616
#>   [8,]  -2.2600726   1.6432877    19.17995  4.872516607 -10.96980937
#>   [9,]  15.0654188   4.5265723    19.33312 -0.213168183  -7.15268938
#>  [10,]  18.3006232   1.7051605    17.39383  0.185257011  -2.83274104
#>  [11,]  12.2453946   5.0546480    20.19254 -0.454197104  -7.79487980
#>  [12,]  15.3830583   4.2361157    18.49756 -0.546306326  -3.87180673
#>  [13,]  11.7108131   4.6126807    20.26392  0.507963320  -9.22366906
#>  [14,]  14.3931286   4.6860440    19.98276 -0.155084529  -8.67736207
#>  [15,]   8.9151260   6.7641954    21.14676 -0.263104306  -9.96325001
#>  [16,]  10.6084838   7.5270445    20.19330 -0.209594869  -9.24997499
#>  [17,]  11.0351735   3.3972142    21.34780  0.767380708 -11.07656616
#>  [18,]   8.8361168   7.1281097    20.87420 -0.507707323  -8.87251772
#>  [19,]   2.7998105   4.0860062    22.48013  1.754510097 -12.65420015
#>  [20,]  16.8340811   3.3633177    18.34276 -0.075838451  -4.99206167
#>  [21,]  15.4120358   3.3205073    19.22443  0.225893308  -6.47861310
#>  [22,]  13.8357950   5.4635667    20.35957 -0.226837201  -9.61537242
#>  [23,]  -0.7877796   2.5983180    19.86823  3.712856269 -11.27555747
#>  [24,]   9.8873004   5.1312756    21.65734  0.229772478 -11.45524719
#>  [25,]  17.1563309   4.1574899    18.78468  0.113263102  -7.32552760
#>  [26,]  14.8930336   2.9836582    19.21570  0.331257679  -6.29028590
#>  [27,]   3.7255299   4.4241769    22.25154  1.590740146 -12.43626383
#>  [28,]   7.3458147   6.3015107    21.74190  0.502628386 -11.42792051
#>  [29,]   9.9179169   8.2706569    20.30230 -0.562576850  -8.66181904
#>  [30,]  17.0672181   2.0041954    17.67576 -0.144274100  -2.35771140
#>  [31,]   6.7440729   5.6800382    22.06046  0.497178335 -11.79717462
#>  [32,]   6.6616385   5.5506889    21.33138  0.477699245 -11.01087185
#>  [33,]   4.5916767   5.6381794    21.27418  1.549808589 -11.66857130
#>  [34,]   9.1308402   3.9134280    20.62331  0.572514878  -8.64361657
#>  [35,]  12.9245221   6.4095487    20.39941 -0.736568303  -8.66265638
#>  [36,]  14.1111691   4.1395461    20.25885  0.405408150  -9.84901546
#>  [37,]  14.0855949   4.7584687    19.58107 -0.020211598  -7.10362494
#>  [38,]  10.9006575   6.0211029    20.58343 -0.392761530  -8.69686485
#>  [39,]   2.3222556   5.1650656    22.30944  1.559855170 -12.37378185
#>  [40,]  14.2684745   6.6051624    19.46359 -0.487530693  -7.98898714
#>  [41,]  17.1745850  -0.1006841    16.20545  0.360501152   3.29484207
#>  [42,]   7.7160155   8.3510577    20.73816 -0.012053021 -10.02502851
#>  [43,]   6.8143028   6.9737238    20.15230  1.574463696 -10.12356971
#>  [44,]  16.1052837   2.6344313    18.10665 -0.259749291  -3.30498512
#>  [45,]   0.3773680   2.8417513    21.88730  3.057890880 -12.38702195
#>  [46,]  10.7100872   6.2656213    20.16981 -0.539580677  -8.16126716
#>  [47,]  17.0323899   2.7669844    18.28052  0.116535855  -4.84548777
#>  [48,]   7.6444634   5.8212632    21.14083  0.799967020 -10.93956318
#>  [49,]  14.6496376   4.9614738    18.95094 -0.615792847  -5.59617383
#>  [50,]   5.7850301   6.9217024    21.44369  0.669318593 -11.58464705
#>  [51,]  -0.8475368   2.6592485    20.25371  4.241794650 -11.48444044
#>  [52,]   3.2452038   6.4922485    21.52830  1.291441236 -11.99385302
#>  [53,]   4.2451980   5.9312840    21.71294  1.215338134 -11.66236100
#>  [54,]  13.3673572   4.8801562    19.26255 -0.139967199  -6.74684241
#>  [55,]   9.5942945   7.9701409    20.78158 -0.640588975  -9.27832082
#>  [56,]   3.6947468   6.1967159    21.07907  1.223942519 -11.50631345
#>  [57,]  11.3869360   6.1351313    20.14074  0.210525486  -9.63667506
#>  [58,]  15.1189693   4.8812951    19.28994 -0.238838118  -6.97939051
#>  [59,]   8.7927514   6.2946216    20.04045  0.091970491  -7.42933734
#>  [60,]  10.5757855   7.7931783    20.08031 -0.634655150  -7.34513704
#>  [61,]  10.3325752   6.1250471    20.59617  0.486186346 -10.32254039
#>  [62,]   9.9882627   4.5009904    21.49754  0.672627706 -11.27237660
#>  [63,]   7.6762885   6.5924713    20.65055  0.588925424 -10.40936389
#>  [64,]  15.5270247   5.8149306    19.10813 -0.778879740  -6.90274437
#>  [65,]   3.7263214   4.3956537    22.86958  0.936383631 -12.60452102
#>  [66,]   8.8207774   4.2616811    21.62602  0.820562738 -10.76001713
#>  [67,]  13.6015710   4.6112764    20.17738  0.350064908  -9.51009145
#>  [68,]  10.9694811   6.1049555    20.95608  0.044751821 -10.05391744
#>  [69,]  17.5245216   3.1592551    18.37645  0.129044062  -5.91104395
#>  [70,]  13.6156416   4.3745801    19.15357  0.146016938  -7.07385220
#>  [71,]  16.3201119   4.5511712    18.80229 -0.311220659  -6.87554302
#>  [72,]   7.3792690   6.9987350    21.39664  0.331485858 -11.12012958
#>  [73,]  12.6052504   5.3092003    19.41547 -0.215901455  -6.35819859
#>  [74,]  19.1714305   3.0963146    17.68196  0.294345533  -4.93637946
#>  [75,]  12.3879359   7.1151940    20.44765 -0.762521087  -9.07382570
#>  [76,]   9.1842988   6.5076611    20.90071  0.269877410 -10.70953456
#>  [77,]   5.4060235   5.5766460    21.20075  1.795471692 -11.55701859
#>  [78,]  17.3285110   3.5641489    18.82631  0.312337029  -7.22272062
#>  [79,]   7.8600337   5.7656149    21.75787  0.347736931 -11.54797291
#>  [80,]   1.2225477   4.3139697    21.41925  2.489347286 -12.09946064
#>  [81,]   3.1072786   5.5508212    21.42703  1.687189069 -11.49132264
#>  [82,]   1.0383572   3.9944495    21.67453  2.213117836 -12.22586994
#>  [83,]   5.2967949   5.4608154    21.89163  0.752747183 -11.77608635
#>  [84,]  19.2361557   2.9856378    17.33103  0.002656079  -3.63643656
#>  [85,]  18.9603024   0.2702727    16.32755  1.076457177   0.61038753
#>  [86,]   5.9973952   5.5706528    21.68784  0.506519990 -11.48334471
#>  [87,]  13.1939561   4.5430113    20.19322 -0.078752689  -8.81348546
#>  [88,]  13.6063200   5.2511363    19.05132 -0.872104468  -4.77944622
#>  [89,]   3.2219624   4.9241821    22.38650  1.847099103 -12.47887239
#>  [90,]  10.4444092   6.0670206    20.98997 -0.261376079 -10.38624333
#>  [91,]   1.2565735   2.4218718    21.77477  2.285219622 -12.34702464
#>  [92,]  13.9312107   3.3678122    18.55572 -0.759981469  -2.80709416
#>  [93,]  18.4049919   0.8891724    16.98698  0.816861557  -1.26815983
#>  [94,]  12.4868612   6.1038719    20.54078 -0.311756041  -9.66243918
#>  [95,]   7.5374783   6.6805024    20.47074  0.435983202 -10.05737271
#>  [96,]  10.4735950   5.8903613    20.62727  0.143194058  -9.52236405
#>  [97,]   3.3560731   3.7766934    21.92856  1.813604425 -12.00229890
#>  [98,]  12.9621349   7.1651365    20.33189 -0.885233504  -8.86172817
#>  [99,]  17.2106201   4.1306413    18.30677  0.066768749  -5.63126043
#> [100,]  14.3396426   5.8672238    19.17793 -0.514337287  -6.38589837
#> [101,]   3.5680897   3.8721577    21.87765  1.847459548 -11.95111650
#> [102,]  15.1938200   4.7879937    19.60587 -0.472890765  -7.85936381
#> [103,]  12.9096727   4.9177787    19.76257 -0.049254217  -7.60220132
#> [104,]  14.1664366   5.2826508    19.98838 -0.402535346  -8.81519120
#> [105,]  17.7543379   0.9404889    16.52373  0.318991574   1.34878311
#> [106,]   6.7323204   5.9049364    22.25935  0.507417126 -12.17772526
#> [107,]  11.9457880   5.5302588    20.77969 -0.166936353  -9.88103666
#> [108,]  15.0835299   3.8392274    19.13979 -0.320073612  -5.99401482
#> [109,]   7.7740192   7.0587810    21.74395 -0.044422746 -11.26775688
#> [110,]  14.2734091   5.3833037    19.80720 -0.655408277  -7.86785531
#> [111,]  16.3008476   3.3328098    18.78854  0.010266881  -6.06360542
#> [112,]   9.8159780   5.2177983    20.78652 -0.238865386  -9.25729912
#> [113,]  16.5182470   2.4805972    18.15156 -0.133868033  -3.63165992
#> [114,]  -1.1451716   2.1482520    20.81282  4.537680533 -11.89787356
#> [115,]  15.0873325   5.0471307    18.90426 -0.278291539  -5.75283310
#> [116,]   1.4653830   3.3598075    22.75351  2.265493620 -12.87869542
#> [117,]  11.4024441   6.1139299    20.86212 -0.482326226  -9.65523271
#> [118,]   0.5006046   2.7653090    22.04472  2.312140140 -12.45045560
#> [119,]  10.6654639   6.4647913    21.12966 -0.154719952 -10.91267300
#> [120,]   5.4965462   6.7747060    20.96616  0.917143111 -11.25215228
#> [121,]   7.4010387   6.7825246    22.02628  0.042308717 -11.61249440
#> [122,]  13.4521984   5.8953179    19.40831 -0.741474699  -6.31963552
#> [123,]  16.1185679   2.3095646    18.44979  0.210811483  -4.36749399
#> [124,]   6.8163569   6.1316707    21.12859  0.915018976 -11.09661864
#> [125,]   2.0945793   5.1458218    20.69815  2.133910105 -11.57554720
#> [126,]  13.1418490   4.7702265    19.31552 -0.259178495  -5.33513594
#> [127,]   8.5174844   8.8570310    20.57558 -0.313769293 -10.11594495
#> [128,]   4.4523754   6.0625767    21.98140  1.240780807 -12.13110415
#> [129,]  16.1422733   0.8479539    17.17346  0.461091378   0.75419375
#> [130,]   3.9026106   5.5632086    22.18464  1.214228199 -12.27188830
#> [131,]  -0.6766918   1.6880324    20.55074  3.043253552 -11.72019501
#> [132,]   8.2343664   2.7033989    22.12085  1.239742448 -11.74026944
#> [133,]   7.7972600   5.9747237    21.24792 -0.554547443  -9.42164533
#> [134,]   3.8611260   5.6505349    22.39721  1.028847232 -12.22300034
#> [135,]   5.5554630   5.1681253    22.05704  1.521024596 -12.07000649
#> [136,]   4.2844201   4.3822313    21.45357  1.653939241 -11.97337114
#> [137,]  12.5645741   5.8627716    20.41916 -0.341338140  -9.00792519
#> [138,]   2.8436957   4.7099278    22.69696  1.681586042 -12.82086641
#> [139,]  -2.1563522   2.1092542    19.49358  3.820460569 -11.12935397
#> [140,]  14.6509051   2.6455296    18.98162 -0.111035224  -5.42330107
#> [141,]   3.0796455   3.4325583    21.93722  1.605690183 -12.17264432
#> [142,]   0.3926410   2.9875236    20.86047  2.453637383 -11.76466363
#> [143,]  -1.0423726   1.8633639    20.38617  3.463404471 -11.61275488
#> [144,]   6.6040716   5.9530499    21.17102  0.639359543 -11.23349091
#> [145,]  13.6004018   5.8387479    20.40029 -0.404742913  -9.63955583
#> [146,]  11.1289643   7.2200821    20.21759 -0.424730354  -8.85137188
#> [147,]   1.7467875   4.2266122    21.27942  2.525021481 -11.88935869
#> [148,]   7.2529210   6.2462093    21.86439  0.632349077 -11.88093366
#> [149,]  15.2224454   4.5748797    18.86487 -0.332769857  -6.06167050
#> [150,]  -0.1789172   4.0161096    21.04784  2.758109829 -11.97885322
#> [151,]  16.1215578  -1.3214483    15.73554  0.882486505   6.45349283
#> [152,]  12.7945340   2.5963282    19.10994 -0.564067161  -3.54501605
#> [153,]   4.4086768   5.0513454    21.87720  1.177133271 -11.73966379
#> [154,]   7.1666276   4.8672288    21.79997  0.840591208 -11.69787006
#> [155,]   4.0833315   2.4167447    22.42984  1.748362857 -12.54684429
#> [156,]   1.8534443   4.9774601    21.98192  2.010844415 -12.32564985
#> [157,]  14.1044329   4.2220481    19.03633 -0.656969879  -4.81127359
#> [158,]   1.6574775   3.4862960    20.81773  2.581525301 -11.69722748
#> [159,]   8.6479572   6.5750652    21.07396  0.040071643 -10.77397809
#> [160,]   5.7138909   5.3625837    21.64644  1.255562986 -11.72330527
#> [161,]  17.1816456   3.5348498    18.34507 -0.121261038  -5.52021644
#> [162,]  11.6094942   4.5890648    20.80477  0.229152519  -9.60311787
#> [163,]  10.0078546   7.2554404    20.76129 -0.438175067  -9.96000185
#> [164,]  14.1564446   6.0095768    19.31892 -0.490895993  -6.90392946
#> [165,]  14.8510644   3.0219070    18.38892 -0.188697468  -2.79820791
#> [166,]   9.2085427   5.7061144    20.98015  0.387899617 -10.28120835
#> [167,]  13.5225022   6.5098322    20.10378 -0.781589196  -8.75404417
#> [168,]  15.0113265   4.0564119    19.55116  0.028027695  -7.69778964
#> [169,]  15.0095213   2.8781980    18.82188 -0.008319700  -4.48505572
#> [170,]  11.9094947   5.5610046    20.90914 -0.036734702 -10.08015624
#> [171,]  12.3380366   6.9726746    20.19958 -0.664236354  -8.47409199
#> [172,]  10.2875601   6.2819190    20.76317 -0.345022089  -8.95391784
#> [173,]   6.2042195   5.6942568    21.92942  0.496479737 -11.44963452
#> [174,]   3.6473906   5.1212548    20.50993  2.186336499 -11.30341446
#> [175,]   8.2405963   5.6355312    21.46756  0.658479662 -11.34125430
#> [176,]  17.2342452   1.7987154    17.54074 -0.155380496  -2.02252352
#> [177,]   0.8501792   4.7017660    20.64133  2.372300851 -11.67551601
#> [178,]   7.3416484   5.4047025    20.92212  1.638322384 -10.88807817
#> [179,]  12.7025408   6.5692761    19.63802 -0.465583064  -8.41540371
#> [180,]   9.2912503   4.4641728    21.36703  0.605666683 -10.61352893
#> [181,]   6.7935896   3.6958852    21.88733  0.923439090 -11.78677801
#> [182,]  15.1433515   3.9869682    18.39190 -0.878379898  -3.19634185
#> [183,]  12.5109927   5.7727315    19.55790 -0.100909057  -7.72519225
#> [184,]   8.3462840   6.2090954    21.30980  0.370419512 -11.09164213
#> [185,]   6.7249081   5.6078146    21.06810  1.399173285 -11.35422125
#> [186,]   2.5711466   4.5483163    21.61920  1.907308704 -11.95198520
#> [187,]  -1.6043477   1.8810382    21.94809  2.651149180 -12.51562341
#> [188,]  11.7480803   5.9186115    20.87541 -0.107679912  -9.97062741
#> [189,]  -0.1251348   2.9691801    20.64640  3.777264180 -11.76957606
#> [190,]   7.7865195   6.4807172    21.64115  0.050828664 -10.96815697
#> [191,]  11.1133419   7.1227977    20.81429 -0.478089035  -9.88295264
#> [192,]   9.5657197   6.5616726    20.88070 -0.015259322  -9.97175149
#> [193,]   2.9317969   4.7691644    20.78880  2.364040018 -11.66137540
#> [194,]  12.6678878   4.9298884    19.81565 -0.562649938  -7.01680796
#> [195,]  13.7683507   5.0511745    19.39472 -0.670541220  -6.55291478
#> [196,]   5.1237589   5.2133492    21.87869  1.218841091 -11.89195606
#> [197,]   5.0313774   6.7359877    21.22797  1.414965458 -11.58939768
#> [198,]   6.5763626   4.3480000    21.80415  1.004780630 -11.75585147
#> [199,]   9.9048303   6.5271346    21.43587 -0.398734024 -10.62027783
#> [200,]   8.5574529   6.8689470    21.52583  0.169713876 -11.05936720
#> [201,]   9.0614984   6.8900153    21.35630 -0.458827761  -9.81409153
#> [202,]   4.6851726   6.1672124    21.68196  1.105849414 -11.72742893
#> [203,]   8.4804430   5.5914196    20.69067  0.498529662 -10.21336404
#> [204,]   3.7641140   5.6190695    21.38546  1.335863246 -11.83280027
#> [205,]   1.9474171   3.6809205    22.39098  2.522665118 -12.71884778
#> [206,]   5.2512055   6.8005066    22.22688  0.272066643 -11.97081114
#> [207,]   6.5062962   6.0749384    21.28439  0.300196772 -11.19025929
#> [208,]   5.0578460   5.4068729    22.08367  1.072600252 -11.85108716
#> [209,]  13.5212296   4.8876170    19.38929 -0.344429947  -6.57065903
#> [210,]   8.3445643   7.5580633    21.27737 -0.025651740 -11.19372638
#> [211,]  15.3643547   5.2849223    19.56710 -0.413949912  -7.97133996
#> [212,]   7.5141910   5.5641919    21.06647  0.890894232 -11.28977083
#> [213,]  -1.8893847   0.8421338    16.18542  5.823295038  -9.25389901
#> [214,]   1.9755618   4.3158780    21.06782  2.528507307 -11.76429073
#> [215,]   8.3380582   6.6878867    21.30103 -0.207600461 -10.68190543
#> [216,]   8.2895017   7.6097714    21.25357 -0.468735873 -10.05000591
#> [217,]   3.6467864   4.9641832    21.61400  1.481796918 -11.99889256
#> [218,]  17.9873108  -1.0442544    15.67083  0.847275171   4.34327879
#> [219,]   8.9708950   6.2409012    21.16389  0.283976781 -10.57318510
#> [220,]  -0.6355265   3.4124643    20.76760  2.883915564 -11.75579429
#> [221,]   2.5468397   4.4177973    21.63387  1.548703814 -11.77527832
#> [222,]   7.3946110   6.3314235    21.53173  0.607056346 -11.06866701
#> [223,]  -0.8186487   2.8268874    21.30183  3.284315505 -12.09997960
#> [224,]   6.3131949   4.5630343    21.79188  1.292854051 -11.79905089
#> [225,]  11.7001904  -4.1261506    15.39609  1.271055650  13.31764048
#> [226,]   9.5480334   5.1693236    21.77920  0.356728011 -11.34827172
#> [227,]  11.6244636   6.8904994    19.99977 -0.170591287  -8.08434473
#> [228,]  10.8149961   6.5854296    20.93035 -0.316669472 -10.03319868
#> [229,]  14.5177324   4.1389138    19.59104 -0.026786735  -7.14290173
#> [230,]  14.4550807   4.2204916    19.50925 -0.344846752  -6.93171416
#> [231,]  14.9952283   3.9868606    18.77471 -0.386958158  -4.41147056
#> [232,]   9.9507708   5.4247027    21.16231  0.288545607 -10.41633848
#> [233,]  13.5188677   5.9987368    19.15921 -1.275525759  -5.24306137
#> [234,]   4.4610828   6.1072289    21.54517  1.145958246 -11.34917987
#> [235,]   3.2054982   3.3865880    22.86793  1.854613759 -12.91705664
#> [236,]   9.0934667   6.7509964    21.29668 -0.317556092  -9.96185025
#> [237,]   8.6983039   5.3473398    21.43281 -0.045939188 -10.08507819
#> [238,]  12.8338845   5.3736057    20.25071  0.295511674  -9.83693787
#> [239,]  10.8884985   7.9580178    20.31124 -0.601602982  -9.08637814
#> [240,]   0.4430152   3.5137360    21.60047  2.409412931 -12.18931457
#> [241,]  10.8059146   6.3180793    20.83929 -0.395856731  -9.59959132
#> [242,]   8.0741797   7.5367209    20.97589  0.110363909 -10.12613365
#> [243,]  17.5760404   2.5309912    17.42954 -0.334809967  -1.99690768
#> [244,]  15.6930538   2.7917354    18.73347 -0.260881404  -4.86464748
#> [245,]   8.6215070   5.7813202    21.75949  0.275070339 -10.98497716
#> [246,]  15.4211511   5.2822294    18.78716 -0.291462903  -6.07911343
#> [247,]  11.6467798   5.5832526    20.93632  0.148797666 -10.35824641
#> [248,]   9.7004350   6.2251892    21.22881 -0.264500366 -10.35276269
#> [249,]  13.7675388   5.2537889    19.89554 -0.611247163  -7.86156851
#> [250,]   6.3530167   6.4118472    21.67389  0.520646151 -11.65145206
#> [251,]  12.5448030   5.3150677    19.77963 -0.688832674  -6.65779784
#> [252,]   9.7536811   5.3411061    21.08209  0.297349909 -10.54409557
#> [253,]  19.6890083   2.9492711    17.17034  0.353515633  -3.66514108
#> [254,]  -0.2729401   1.5381791    20.45501  3.134181775 -11.66334151
#> [255,]   6.3268062   6.6437763    21.79166  0.184937106 -11.67343921
#> [256,]  17.3318498   3.4304158    18.01435 -0.328944682  -4.22926019
#> [257,]  16.8677830   4.3816796    18.73763 -0.059506074  -6.80325024
#> [258,]   9.8517677   6.7349927    20.44776 -0.063301796  -8.84770578
#> [259,]  10.8415691   7.2096061    20.98692 -0.411055085 -10.35691813
#> [260,]   6.5203498   5.7460078    22.12717  0.593407111 -11.72876465
#> [261,]  10.3316742   3.2218971    21.40260  0.952163341 -11.00380233
#> [262,]  14.1101024   5.3591485    19.73946 -0.411245936  -7.60968966
#> [263,]  10.7113440   7.1416921    19.92857 -0.771908135  -7.66828938
#> [264,]   8.0123863   4.6274413    21.90979  0.773030645 -11.61472607
#> [265,]  14.3565177   4.5045014    19.04930  0.074043898  -5.82585865
#> [266,]  14.3262218   6.1204362    19.61302 -0.669036184  -7.68591116
#> [267,]   4.5699392   4.8339044    22.12163  1.323264059 -11.88032939
#> [268,]   9.5537411   5.3634663    20.51324  1.075497503 -10.11912133
#> [269,]  11.9538661   5.6198060    20.28204 -0.475413256  -7.77254501
#> [270,]   5.7611808   5.6510078    22.00673  0.610226931 -11.07160506
#> [271,]  11.6166009   5.9096759    20.61777 -0.229378639 -10.01008448
#> [272,]   2.7020826   5.3229445    22.42988  1.678449179 -12.68907167
#> [273,]   8.8133403   4.5247297    22.16841  0.454970990 -11.66878678
#> [274,]  15.8901608   2.9560630    18.78360 -0.210287089  -5.56774300
#> [275,]  11.6296099   5.8346856    20.89623 -0.127337931 -10.10034124
#> [276,]  12.6132953   6.2283487    19.57273 -1.019776814  -6.73175491
#> [277,]  17.8982738   1.1964868    16.78032  0.191602850   0.10992629
#> [278,]   3.4780178   7.0142971    20.16467  1.131688774 -10.91752625
#> [279,]   1.1060852   3.7621698    20.97282  2.922946794 -11.89003891
#> [280,]  14.3456519   3.3771061    19.41167 -0.011007953  -6.08699600
#> [281,]  14.2708128   4.4286426    19.13484 -0.668008594  -5.36366774
#> [282,]   7.5476964   4.7953216    22.24607  0.620223175 -11.85530789
#> [283,]   7.9725019   6.5206726    21.17772  0.079309788 -10.59041589
#> [284,]   9.5819703   3.2670543    21.71457  0.908819062 -11.21989359
#> [285,]  12.3849812   4.2000702    19.32150 -0.364447363  -4.55828297
#> [286,]  -1.9575869   2.3109500    20.98552  3.471077130 -11.97570867
#> [287,]  16.3546509   1.1194218    17.32968  0.436275676  -0.28379919
#> [288,]  -2.8658436   1.4451616    19.35622  5.032329432 -11.07679878
#> [289,]  10.3119458   6.5887156    21.11391 -0.448028812  -9.57143237
#> [290,]   3.4005555   4.8256553    22.17407  1.746186241 -12.40949166
#> [291,]  16.4280564   3.2113300    18.56822 -0.121659454  -5.13604051
#> [292,]   4.5619541   6.1033414    21.80335  1.137517189 -11.95909866
#> [293,]  11.4053222   5.6992446    20.30731 -0.443422384  -8.30093921
#> [294,]   9.4299244   6.2008978    20.71281  0.108262175  -9.90670717
#> [295,]  13.9196293   5.9510779    19.79721 -0.856018340  -7.51273334
#> [296,]  13.5431014   5.5212491    19.87054 -0.596882068  -7.65266805
#> [297,]   9.1937419   7.0097044    20.20848  0.109608101  -9.31369033
#> [298,]   8.0160017   5.4961629    21.61087  0.365693296 -11.28128366
#> [299,]   6.8437882   4.6059828    22.27183  1.090256970 -12.07624140
#> [300,]  15.5087860   5.2824090    19.43461 -0.510285643  -7.79552621
#> [301,]   6.3641110   4.5476675    22.54285  0.703606662 -12.12596592
#> [302,]  11.7580414   6.9920033    20.35953 -0.825859550  -8.53162749
#> [303,]   2.4210085   5.2472945    22.85688  0.998096346 -12.70763664
#> [304,]  17.1230816   2.2920521    17.73304 -0.188466521  -2.80754814
#> [305,]  14.1423833   3.8115168    19.08999 -0.653959101  -4.71827102
#> [306,]   5.7833900   5.2220736    21.69699  0.825839959 -11.60019384
#> [307,]   9.5153611   7.2566237    21.04283 -0.101762105 -10.44601590
#> [308,]   2.1834680   5.5545273    21.53168  1.674560030 -11.96033912
#> [309,]   9.8023401   6.8719746    21.46899 -0.193126801 -11.10043192
#> [310,]   6.5885512   5.2753503    22.72778  0.557814656 -12.53835296
#> [311,]   3.3821689   5.1274762    21.82532  1.413335876 -12.13299797
#> [312,]  -1.2304062   2.9849728    20.99369  3.476136876 -11.97018690
#> [313,]  16.3120839   3.2122717    18.76615 -0.188787870  -5.76148144
#> [314,]   4.3934733   5.9269215    21.84625  0.836919761 -11.51170062
#> [315,]   3.4763725   6.0850429    21.16838  1.598605514 -11.25342108
#> [316,]   4.8583802   7.1882954    21.77903  0.858097553 -11.79194909
#> [317,]  16.2955485   4.3401017    18.57131 -0.394477405  -5.22821415
#> [318,]   9.2867636   7.0537591    21.36999 -0.312625360 -10.61377192
#> [319,]   0.3978228   2.8747734    22.20759  2.455337552 -12.63534639
#> [320,]  -0.7043640   3.0534958    20.27327  3.348265181 -11.53285995
#> [321,]   7.9620390   6.5752846    20.94485  0.290955007  -9.83980796
#> [322,]   5.4098287   6.8147071    21.38639  0.673017832 -11.49376241
#> [323,]  11.1241918   7.2524835    20.31064 -0.578298909  -8.79629511
#> [324,]   6.4340313   6.1647355    20.61267  1.900889495 -10.97949683
#> [325,]  15.1097292   5.0524504    19.41830 -0.043402021  -8.03175778
#> [326,]   7.4268105   7.0191105    21.91830 -0.216513815 -10.93993746
#> [327,]   7.9265701   4.8257538    22.01904  0.477270974 -11.54908538
#> [328,]  12.8220281   5.7675660    19.70559 -0.948521306  -6.63109481
#> [329,]  11.5987038   5.7288615    20.79760 -0.325949626  -9.22797562
#> [330,]  16.8778927   3.9298955    18.04867 -0.165946347  -4.50213687
#> [331,]  16.6066780  -0.8671473    16.04478  0.664335762   4.65438765
#> [332,]   5.9970167   7.7114265    21.56537  0.570687195 -11.58545220
#> [333,]  19.5275344   1.8250429    16.80056  0.693207590  -1.91445333
#> [334,]  15.6950169   4.2202823    19.33545 -0.115725379  -7.36638218
#> [335,]  10.5427533   6.7936302    20.51646 -0.012624308 -10.29109353
#> [336,]   6.7989649   6.2676245    21.01903  0.784199525 -11.00320507
#> [337,]  14.2614180   4.1799391    19.40200  0.002918540  -6.92677614
#> [338,]  11.3067258   5.5300008    20.92542 -0.047534686 -10.19390243
#> [339,]   5.3896485   4.0005504    22.63777  1.357865524 -12.71847646
#> [340,]  10.8193449   6.0565783    20.79763 -0.359915129  -8.99480018
#> [341,]   0.7013823   3.9920356    20.26255  3.113779559 -11.47135100
#> [342,]  11.1273542   5.0545409    20.59582  0.343156868  -9.96516041
#> [343,]   7.3760235   6.7859854    20.81324  0.284891458  -9.71078091
#> [344,]  10.6338885   5.8505095    20.75091 -0.418041630  -8.93274452
#> [345,]   5.9631209   6.3725906    21.29559  0.737269671 -11.39394871
#> [346,]  13.3408734   4.8183943    20.36688 -0.120942423  -9.04692996
#> [347,]  10.3725859   5.6284343    21.21747 -0.161404849 -10.36800140
#> [348,]  15.7257894  -0.7350107    17.10078  0.268219507   1.66417827
#> [349,]  10.3650005   7.0498551    20.18386 -0.785604910  -8.41581120
#> [350,]   7.3470374   5.0215925    21.79176  0.738642481 -11.42639467
#> [351,]  10.8516562   6.7251473    20.96108 -0.497449286  -9.41850279
#> [352,]   4.6889525   6.2380546    21.85468  0.815107216 -11.66013706
#> [353,]   3.8723132   5.9227528    21.74868  1.153025845 -11.98164989
#> [354,]  18.0257268   2.6379442    17.94711  0.027186749  -4.54209082
#> [355,]   4.9325727   5.8232334    22.65138  0.737398245 -12.60754690
#> [356,]   8.0726934   5.8045716    21.02283  0.207458275 -10.15021970
#> [357,]  11.6164438   5.1095604    20.12160  0.194927046  -8.45060006
#> [358,]  10.0922198   6.1535609    21.29617  0.214605606 -10.97533772
#> [359,]  15.3503694   4.5878376    18.73589 -0.889896466  -4.67557639
#> [360,]  12.5959024   5.6662102    20.24356 -0.476992596  -8.32654254
#> [361,]   8.3911358   4.4369414    21.19583  0.702687003 -10.43731802
#> [362,]  17.6482312   1.0895396    16.87353  0.083695120   0.13109953
#> [363,]   9.3192114   6.7700994    20.25928 -0.426090180  -8.40376431
#> [364,]  12.0655439   5.7507298    20.69718 -0.414276143  -9.38679054
#> [365,]   1.4290701   3.4586289    21.19897  3.010979350 -11.91225873
#> [366,]  10.2319118   5.4320986    21.07181  0.053815492  -9.77109109
#> [367,]  13.4552006   5.0016922    19.95238 -0.477666866  -7.74460509
#> [368,]  14.8520263   2.9043594    19.74211  0.801511214  -8.89794725
#> [369,]   8.1414116   5.8790869    21.34050  0.311394570 -10.70603661
#> [370,]  15.3093274   5.1562383    19.30381 -0.437374555  -7.42018314
#> [371,]   7.5242639   6.0398774    21.12295  0.334395752 -10.68705418
#> [372,]  16.6632374   2.9956956    18.18640  0.027469804  -3.90369017
#> [373,]   0.7879645   4.2070281    21.23259  2.659201529 -11.85585278
#> [374,]  14.8384219   4.4489251    19.87550  0.286420039  -8.83886108
#> [375,]   2.2861287   4.0980591    21.25101  2.745217522 -11.92829012
#> [376,]  14.3514833   5.7331873    19.17692 -0.532997592  -6.69073910
#> [377,]  12.5934919   5.4438603    20.18290  0.046341107  -8.70108628
#> [378,]  13.6491448   5.4889938    19.63595 -0.324141183  -8.16466329
#> [379,]  13.6681530   5.7743142    19.51770 -0.648543250  -6.94582342
#> [380,]   8.6666929   5.6853576    20.96582  0.820595353 -10.77123706
#> [381,]   8.2920347   7.2300749    21.41512 -0.484025048 -10.80286651
#> [382,]   2.5999710   3.7642071    20.77104  2.062500647 -11.11269878
#> [383,]   1.4667346   4.3809023    19.96803  3.314911067 -11.24071352
#> [384,]   8.9396876   4.6111019    20.65020  0.675847679 -10.23117619
#> [385,]   9.7689552   4.3411505    21.10656  0.986320616 -10.71977005
#> [386,]   2.5909848   4.7748754    22.27402  1.837242288 -12.48777584
#> [387,]  15.0560253   6.2774340    18.87119 -0.856518411  -5.57924928
#> [388,]  13.8957361   5.7338824    19.87692 -0.723708298  -8.09091459
#> [389,]  15.8178897   4.8799831    19.08186 -0.054036788  -7.32016548
#> [390,]   5.1531567   5.3261957    20.44703  2.028529993 -10.70372172
#> [391,]   9.0119391   7.8190056    20.85142 -0.237007078 -10.29950993
#> [392,]   3.0452590   4.4981479    21.82092  1.657556038 -12.17085107
#> [393,]   1.0078938   3.7106704    20.97908  2.320461842 -11.58831705
#> [394,]   4.0339725   3.4106597    22.20138  1.711559935 -12.45666255
#> [395,]  14.6728214   4.2266631    19.55280  0.023053333  -7.66688411
#> [396,]   7.6028405   4.4708163    22.23167  0.823105568 -11.98515585
#> [397,]   8.6698868   6.1216208    21.45280  0.057739494 -10.02307506
#> [398,]   6.2962616   7.6503806    21.29102  0.194048207 -10.68583500
#> [399,]   6.8220483   7.4825293    22.01331 -0.038093681 -11.66410142
#> [400,]   4.0713644   7.0099503    21.03666  1.201832081 -11.37042749
#> [401,]  12.3346208   7.4692227    20.15356 -1.026143042  -8.11986729
#> [402,]  -1.8077530   1.3551081    17.09205  4.879724970  -9.76530542
#> [403,]  15.6789129  -2.5322614    15.23718  0.979034806   8.91512421
#> [404,]  14.9220672   4.9702031    19.19232 -0.329808389  -6.30956790
#> [405,]  -1.0619768   2.8403408    20.19390  3.827080559 -11.50421328
#> [406,]  16.3691134  -0.1401553    16.98425  0.169550543   1.18594330
#> [407,]   5.6222064   6.1052864    22.10979  0.647943075 -11.91503121
#> [408,]  15.7459720   5.1310622    19.16198 -0.127473780  -7.65751752
#> [409,]  17.1203494   4.2651782    18.06568 -0.285030790  -4.81080368
#> [410,]  14.4971101   4.7713819    19.15254 -0.623800476  -5.38463610
#> [411,]   9.9883030   6.1341212    20.60956  0.176589507  -9.51333340
#> [412,]  11.3580064   7.0307284    20.48851 -0.685047000  -8.65291345
#> [413,]   4.4226108   5.1729928    21.84868  1.131843267 -11.98686883
#> [414,]   7.8627677   5.8912429    21.49221  0.325899126 -11.22032658
#> [415,]  12.0843236   5.7699470    20.21879 -0.296785907  -9.07584701
#> [416,]  14.7732589   5.4875392    19.29985 -0.357470532  -7.32948600
#> [417,]  10.5975351   7.0245102    20.97680  0.029951661 -10.76138648
#> [418,]  14.3372859  -1.0695625    16.81203  0.181876462   4.52400483
#> [419,]  12.0944913   5.8355606    20.61280 -0.073126084  -9.41108294
#> [420,]  10.2925203   6.6924900    20.65437 -0.306107347  -9.33830655
#> [421,]  17.2095072   1.2755673    16.98881  0.390549230   0.19445870
#> [422,]   2.3978979   4.3344818    21.47391  2.432603038 -12.06425717
#> [423,]  15.8431814   0.7337682    17.71704 -0.159335802  -0.96383406
#> [424,]  12.8946390   3.6080379    19.57859  0.234281954  -6.69715477
#> [425,]   8.2674335   6.8747915    20.74797  1.101171415 -10.82493433
#> [426,]   9.0574937   6.3521746    20.56511  0.530547191 -10.43953915
#> [427,]   5.0621979   6.1743366    21.54755  0.801481494 -11.56451722
#> [428,]   1.5618598   3.6179520    22.26983  2.137402470 -12.62280710
#> [429,]  15.3282549  -1.9268103    16.73637  0.696273141   3.59289694
#> [430,]  16.7236446   3.1815989    17.79097 -0.070208823  -2.59849513
#> [431,]   4.3211454   6.2988583    20.22372  1.595215995 -10.98052051
#> [432,]  15.7600376   5.2448509    18.71584 -0.691394206  -5.19422456
#> [433,]   9.4162472   6.1897406    21.02299 -0.778619181  -9.07500707
#> [434,]  11.9657417   5.3747911    20.51307  0.240387886  -9.44778893
#> [435,]   2.3611607   4.8041574    20.47556  2.592942116 -11.36188380
#> [436,]  10.9646965   6.4733315    20.31019 -0.323730672  -7.97965734
#> [437,]   1.0551580   3.9985347    21.66824  2.045084966 -11.99958603
#> [438,]  10.5470585   6.4633655    20.78776 -0.713789701  -8.92483258
#> [439,]   6.4216266   5.8131506    20.78170  1.551377149 -11.09178514
#> [440,]  10.6944118   5.7665055    21.01960  0.053915280 -10.47965197
#> [441,]  -0.6454678   2.7294430    21.00274  2.971097752 -11.91647000
#> [442,]   3.4640131   6.8096117    20.49563  2.110637066 -11.28092613
#> [443,]   9.8898126   7.0178490    20.70823  0.193901530 -10.28398110
#> [444,]   5.0481798   3.5081871    22.40638  1.282137562 -12.36388475
#> [445,]   5.8020075   7.0678192    21.91861  0.170288175 -11.65843592
#> [446,]  12.7554167   6.0802052    20.15974 -0.617358499  -8.45307932
#> [447,]   6.1332963   7.2689376    20.99975  0.388070222 -10.48214795
#> [448,]  11.7225429   4.8234619    20.18633  0.133078078  -8.28506831
#> [449,]   5.3344931   4.2296432    21.08527  1.485249051 -11.33821369
#> [450,]   2.0835592   2.5397438    22.04896  2.112085766 -12.39758815
#> [451,]  12.0521495   5.2163581    19.82150 -0.632770280  -5.97732069
#> [452,]   6.1878676   5.4743611    22.03676  1.241938000 -12.21547988
#> [453,]  -1.4864358   2.9805865    20.99212  3.572727504 -11.99555661
#> [454,]  16.6832860   4.1249918    18.60877 -0.307790843  -5.85336029
#> [455,]   9.0198428   5.4654575    21.54642  0.505551304 -11.27226894
#> [456,]   6.1370600   5.2412422    21.88981  0.648757623 -11.46413142
#> [457,]  14.0358508   5.6150492    20.15759 -0.110999335  -9.43179988
#> [458,]  15.8293800   2.4314217    17.45858 -0.066798930  -0.02729488
#> [459,]   3.7637940   3.9566652    22.45379  1.888558338 -12.60368535
#> [460,]  10.5922832   4.3566786    20.68540  0.377187269  -9.71236067
#> [461,]  15.8046084   4.1971369    19.14566  0.057678570  -6.93490891
#> [462,]  12.9571122   5.8696851    20.27737 -0.399010865  -8.99221664
#> [463,]  10.5307594   7.1175584    20.71824 -0.685674343  -8.66269022
#> [464,]  15.8362001   3.3195884    18.21574 -0.605497236  -3.00198904
#> [465,]  17.9692521   4.3541463    18.33883 -0.362900518  -6.06053861
#> [466,]  13.5282236   5.1382310    20.20209 -0.164628091  -8.67811119
#> [467,]  10.6936711   7.0205431    21.15291 -0.690483767 -10.12623617
#> [468,]   9.5656568   5.7606843    20.90417  0.267767386 -10.41688104
#> [469,]   8.6886401   5.0976435    21.30961  0.441590141  -9.92716755
#> [470,]  -1.8862152   0.8949068    17.79542  5.391983007 -10.16688665
#> [471,]  13.4114635   5.6933336    19.75855 -0.313364951  -8.29060371
#> [472,]  16.5492561   2.5506918    17.60607 -0.117634333  -1.49255501
#> [473,]  15.1511829   3.9529361    18.93890 -0.227224663  -5.94898430
#> [474,]  12.5208451   5.5678467    19.77019 -0.404488046  -7.35764662
#> [475,]   5.2708631   6.5945500    21.34405  0.555722700 -11.13973738
#> [476,]  13.5889873   4.6058377    20.42862 -0.134175613  -9.17777079
#> [477,]  12.1473917   6.3833223    19.66279 -0.554741070  -7.41694558
#> [478,]  13.5810752   4.6094161    19.08217 -0.623285714  -5.12242702
#> [479,]   5.0653011   6.4014555    20.73899  1.430427751 -11.20571078
#> [480,]  13.3922461   4.5915717    19.43375 -0.507959331  -6.03949446
#> [481,]   6.1384448   4.6576329    21.24737  1.133824573 -11.19533153
#> [482,]  10.1443583   4.4976779    21.11925  0.237534264  -9.62443365
#> [483,]  10.9607587   5.3730522    20.79395  0.315620352  -9.61011223
#> [484,]  10.2165665   6.0599463    21.22839 -0.517472844  -9.51500493
#> [485,]  -0.4892663   4.3847283    19.93835  3.622293047 -11.30513853
#> [486,]  12.8260467   2.8630024    20.67545  0.795513182  -9.46931481
#> [487,]   8.2889552   5.6792661    21.21300  0.447975786 -10.35703068
#> [488,]   7.2182383   6.9591627    21.45026  0.263830714 -11.42595980
#> [489,]  -1.3074667   0.8941275    19.00541  4.076764584 -10.83316474
#> [490,]   5.7623418   3.7749611    21.35869  2.125402672 -11.66904653
#> [491,]   7.1891650   6.0404103    21.85100  0.233424186 -11.32041312
#> [492,]  11.7097869   5.6690187    20.56938 -0.586321911  -8.52321732
#> [493,]  16.1138782   4.1938550    18.87280 -0.250041681  -6.21996806
#> [494,]  10.9047316   5.0534616    20.77284 -0.098709088  -8.94355885
#> [495,]  12.5183968   7.6213511    20.32202 -0.849134279  -9.23705490
#> [496,]   9.3985625   6.4869506    20.96010 -0.247905158  -9.92414726
#> [497,]   0.3921751   1.7405172    23.26763  2.080945191 -13.26463080
#> [498,]   1.7783407   3.3117971    21.86205  1.988238712 -12.27698621
#> [499,]   6.4926873   4.2084303    21.81486  0.783644037 -11.24945888
#> [500,]   8.7294253   5.7694763    21.45407  0.301666234 -10.45418917
#> [501,]   2.8714095   6.6968982    21.63423  0.929756753 -11.76547208
#> [502,]  14.5096908   4.8634148    19.25879 -0.857729516  -6.09226685
#> [503,]  14.3982309   2.4065033    18.02296 -0.272454514  -1.17218869
#> [504,]  10.8389013   6.2616044    21.00369 -0.166022973  -9.93779544
#> [505,]  -1.5555294   1.5568946    22.50372  2.484005877 -12.86204457
#> [506,]   9.4324207   7.4378682    20.96422 -0.836099890  -9.64209427
#> [507,]  15.9993426   4.0799824    19.00468  0.085702629  -6.75439029
#> [508,]   7.9509632   7.0549204    21.49905  0.035286060 -10.88259386
#> [509,]  14.7102365   3.7332764    19.13720 -0.398455039  -5.70133985
#> [510,]  12.9858891   3.7690782    20.19050  0.146887373  -8.08756797
#> [511,]   4.7068379   6.7106541    20.91349  0.847166659 -11.02511996
#> [512,]   8.5910651   4.8143750    21.84597  0.792504195 -11.78280215
#> [513,]   3.5137817   5.6371141    21.97113  1.114607223 -12.01745061
#> [514,]   7.6014850   6.8902065    21.34420  0.243677416 -10.84476465
#> [515,]   9.9717041   5.2820272    21.23517  0.429004837 -10.37722079
#> [516,]   9.7302570   7.1147539    21.16763 -0.145107749 -10.34968347
#> [517,]   8.7110950   7.9764031    21.32106 -0.780070659 -10.24385593
#> [518,]  14.0721540   5.7344365    19.45863 -0.510280147  -7.00100075
#> [519,]  14.9150431   4.1805928    19.51404  0.343495760  -8.36733284
#> [520,]   7.7654924   6.4948167    21.01858  0.581140362 -10.85053825
#>         tmean_v3.l2
#>   [1,]  2.663669257
#>   [2,]  0.412627877
#>   [3,] -0.274256660
#>   [4,] -0.116583536
#>   [5,]  5.380644408
#>   [6,]  4.141714105
#>   [7,]  2.151598241
#>   [8,] -2.783918840
#>   [9,]  2.735411937
#>  [10,]  4.136810907
#>  [11,]  3.193970514
#>  [12,]  4.192584145
#>  [13,]  0.842785997
#>  [14,]  2.000951116
#>  [15,]  2.128802879
#>  [16,]  1.609830566
#>  [17,]  0.176661018
#>  [18,]  2.684843474
#>  [19,] -0.887363656
#>  [20,]  3.310607759
#>  [21,]  2.125185473
#>  [22,]  1.568345954
#>  [23,] -2.057836181
#>  [24,]  0.759475765
#>  [25,]  1.746307861
#>  [26,]  2.239738445
#>  [27,] -0.741903956
#>  [28,]  0.466622682
#>  [29,]  2.558442412
#>  [30,]  5.069384092
#>  [31,]  0.324731730
#>  [32,]  0.729749591
#>  [33,] -0.556928336
#>  [34,]  1.506727330
#>  [35,]  2.615898716
#>  [36,]  0.728922011
#>  [37,]  1.983322106
#>  [38,]  2.373745059
#>  [39,] -0.553658799
#>  [40,]  2.329111176
#>  [41,]  5.706176651
#>  [42,]  1.678949041
#>  [43,]  0.018085204
#>  [44,]  4.977715420
#>  [45,] -1.667296571
#>  [46,]  3.091792162
#>  [47,]  3.201497541
#>  [48,]  0.343243327
#>  [49,]  4.119028967
#>  [50,]  0.194698611
#>  [51,] -2.369071790
#>  [52,] -0.387883503
#>  [53,] -0.147927707
#>  [54,]  2.743790500
#>  [55,]  2.348411412
#>  [56,] -0.111169549
#>  [57,]  1.068403980
#>  [58,]  2.423607632
#>  [59,]  2.247590279
#>  [60,]  2.797014680
#>  [61,]  0.647183674
#>  [62,]  0.026846856
#>  [63,]  0.600812363
#>  [64,]  3.387849610
#>  [65,] -0.129950689
#>  [66,]  0.540630121
#>  [67,]  0.673620602
#>  [68,]  1.160628727
#>  [69,]  2.701287385
#>  [70,]  2.599235382
#>  [71,]  2.890977940
#>  [72,]  0.782480547
#>  [73,]  2.990049045
#>  [74,]  1.996483295
#>  [75,]  2.206767697
#>  [76,]  0.677525327
#>  [77,] -0.704286090
#>  [78,]  1.408670633
#>  [79,]  0.562877337
#>  [80,] -1.306494061
#>  [81,] -0.312917678
#>  [82,] -1.113290037
#>  [83,]  0.319898443
#>  [84,]  3.231933879
#>  [85,]  2.488779420
#>  [86,]  0.503999473
#>  [87,]  2.307454045
#>  [88,]  4.718910364
#>  [89,] -0.871747705
#>  [90,]  1.754930100
#>  [91,] -1.249468698
#>  [92,]  6.468352730
#>  [93,]  2.686453498
#>  [94,]  1.870829626
#>  [95,]  1.055626485
#>  [96,]  1.274466347
#>  [97,] -0.735830206
#>  [98,]  2.591880238
#>  [99,]  1.989924709
#> [100,]  2.745120330
#> [101,] -0.704677721
#> [102,]  3.116731840
#> [103,]  2.198629487
#> [104,]  2.385233354
#> [105,]  4.554000207
#> [106,]  0.176134064
#> [107,]  1.585782131
#> [108,]  3.589042546
#> [109,]  1.070183064
#> [110,]  3.157428571
#> [111,]  2.975026712
#> [112,]  2.454552213
#> [113,]  4.471838944
#> [114,] -2.590947044
#> [115,]  2.515074571
#> [116,] -1.210267010
#> [117,]  2.147465312
#> [118,] -1.234173719
#> [119,]  1.079809789
#> [120,]  0.253761930
#> [121,]  0.888862203
#> [122,]  4.214991969
#> [123,]  3.335431538
#> [124,]  0.163902711
#> [125,] -1.010583649
#> [126,]  3.054412180
#> [127,]  1.690202732
#> [128,] -0.420494008
#> [129,]  4.196693810
#> [130,] -0.324094132
#> [131,] -1.740134449
#> [132,] -0.732282702
#> [133,]  2.657584872
#> [134,] -0.134021087
#> [135,] -0.556993785
#> [136,] -0.721174714
#> [137,]  1.920460318
#> [138,] -0.816073664
#> [139,] -2.160526996
#> [140,]  4.266503259
#> [141,] -0.673098083
#> [142,] -1.290468050
#> [143,] -1.964538111
#> [144,]  0.417910945
#> [145,]  1.867928770
#> [146,]  2.171383242
#> [147,] -1.265220132
#> [148,]  0.040029443
#> [149,]  3.267848528
#> [150,] -1.515415719
#> [151,]  5.229125701
#> [152,]  6.385241373
#> [153,] -0.232429963
#> [154,]  0.054051499
#> [155,] -0.859396861
#> [156,] -0.903007247
#> [157,]  4.774525911
#> [158,] -1.306365646
#> [159,]  1.133258574
#> [160,] -0.357382587
#> [161,]  3.220271684
#> [162,]  1.530236132
#> [163,]  2.031278928
#> [164,]  2.739015249
#> [165,]  4.318845266
#> [166,]  0.783119791
#> [167,]  2.694895334
#> [168,]  2.361251751
#> [169,]  3.515957630
#> [170,]  1.094243403
#> [171,]  2.399227812
#> [172,]  2.200298963
#> [173,]  0.727189873
#> [174,] -0.922728187
#> [175,]  0.282629594
#> [176,]  5.353429121
#> [177,] -1.212842404
#> [178,] -0.248914724
#> [179,]  2.396564614
#> [180,]  0.379706272
#> [181,] -0.051589329
#> [182,]  5.617441375
#> [183,]  2.106533362
#> [184,]  0.814089383
#> [185,] -0.229665007
#> [186,] -0.779471892
#> [187,] -1.492964055
#> [188,]  1.391155616
#> [189,] -2.136048662
#> [190,]  1.106241472
#> [191,]  1.740576046
#> [192,]  1.317042900
#> [193,] -1.174987059
#> [194,]  3.638147554
#> [195,]  3.966831399
#> [196,] -0.239145548
#> [197,] -0.391767754
#> [198,]  0.003761584
#> [199,]  1.688373096
#> [200,]  0.654961662
#> [201,]  2.134687860
#> [202,] -0.061847894
#> [203,]  0.821102645
#> [204,] -0.394687522
#> [205,] -1.393016308
#> [206,]  0.525079954
#> [207,]  0.672929697
#> [208,] -0.120004746
#> [209,]  3.189055452
#> [210,]  0.858289828
#> [211,]  2.326086301
#> [212,]  0.058500634
#> [213,] -3.325254948
#> [214,] -1.263139247
#> [215,]  1.485558333
#> [216,]  2.220344691
#> [217,] -0.536665795
#> [218,]  5.016241256
#> [219,]  0.732626750
#> [220,] -1.548411914
#> [221,] -0.393864082
#> [222,]  0.282792216
#> [223,] -1.801349929
#> [224,] -0.416809157
#> [225,]  7.275534150
#> [226,]  0.468231347
#> [227,]  1.536036527
#> [228,]  1.505009307
#> [229,]  2.303898345
#> [230,]  3.248901593
#> [231,]  3.796508454
#> [232,]  0.988907775
#> [233,]  5.425944412
#> [234,]  0.317800549
#> [235,] -0.977984228
#> [236,]  1.831610220
#> [237,]  1.813357724
#> [238,]  0.780077494
#> [239,]  2.367488615
#> [240,] -1.244315682
#> [241,]  1.852116431
#> [242,]  1.258081491
#> [243,]  5.128292893
#> [244,]  4.565078736
#> [245,]  0.587450528
#> [246,]  2.435790533
#> [247,]  0.885126531
#> [248,]  1.645100447
#> [249,]  3.387703835
#> [250,]  0.489124882
#> [251,]  4.061542457
#> [252,]  0.916600715
#> [253,]  1.960647088
#> [254,] -1.776890606
#> [255,]  0.752993459
#> [256,]  4.139204302
#> [257,]  2.155972134
#> [258,]  1.836884064
#> [259,]  1.408846235
#> [260,]  0.401895918
#> [261,]  0.166777879
#> [262,]  2.592986289
#> [263,]  3.315274298
#> [264,]  0.105937468
#> [265,]  2.022649029
#> [266,]  2.699986167
#> [267,] -0.275323381
#> [268,]  0.160946456
#> [269,]  2.886831653
#> [270,]  0.682711608
#> [271,]  1.661423277
#> [272,] -0.813222225
#> [273,]  0.481459866
#> [274,]  4.221806531
#> [275,]  1.472212977
#> [276,]  4.366266119
#> [277,]  4.751762976
#> [278,] -0.099076619
#> [279,] -1.605013209
#> [280,]  3.022042187
#> [281,]  4.452369288
#> [282,]  0.287644610
#> [283,]  1.324791967
#> [284,]  0.233199956
#> [285,]  3.907320785
#> [286,] -1.964669502
#> [287,]  4.075073558
#> [288,] -2.879330988
#> [289,]  2.165918644
#> [290,] -0.780541824
#> [291,]  3.565338373
#> [292,] -0.272511344
#> [293,]  3.151705687
#> [294,]  1.206119029
#> [295,]  3.492194580
#> [296,]  3.037922998
#> [297,]  1.697053355
#> [298,]  0.658753791
#> [299,] -0.194369548
#> [300,]  2.683146973
#> [301,]  0.285919417
#> [302,]  2.766428758
#> [303,] -0.186969394
#> [304,]  4.993162214
#> [305,]  5.043744432
#> [306,]  0.130217648
#> [307,]  1.136508732
#> [308,] -0.643211839
#> [309,]  1.039080529
#> [310,]  0.086034591
#> [311,] -0.485170093
#> [312,] -1.953638721
#> [313,]  3.755343111
#> [314,]  0.262461448
#> [315,] -0.268157170
#> [316,]  0.112092615
#> [317,]  3.323206777
#> [318,]  1.537453694
#> [319,] -1.348200863
#> [320,] -1.853133424
#> [321,]  1.375277187
#> [322,]  0.316573898
#> [323,]  2.343718291
#> [324,] -0.649543928
#> [325,]  1.904720394
#> [326,]  1.429580660
#> [327,]  0.568773975
#> [328,]  4.330679165
#> [329,]  2.114013280
#> [330,]  3.222877537
#> [331,]  5.493442278
#> [332,]  0.462841847
#> [333,]  1.997231403
#> [334,]  2.337016773
#> [335,]  1.114585001
#> [336,]  0.272763816
#> [337,]  2.701114319
#> [338,]  1.409787475
#> [339,] -0.596277739
#> [340,]  2.243096571
#> [341,] -1.709279607
#> [342,]  1.109503999
#> [343,]  1.397124177
#> [344,]  2.594406598
#> [345,]  0.167606909
#> [346,]  1.815337718
#> [347,]  1.497538098
#> [348,]  6.864010159
#> [349,]  3.033240874
#> [350,]  0.238353145
#> [351,]  1.897848568
#> [352,]  0.310037463
#> [353,] -0.310252682
#> [354,]  3.594286447
#> [355,] -0.109927324
#> [356,]  1.301990680
#> [357,]  1.260883409
#> [358,]  0.753256116
#> [359,]  4.965308340
#> [360,]  2.701663129
#> [361,]  0.659212271
#> [362,]  5.279031383
#> [363,]  2.828083940
#> [364,]  2.340530620
#> [365,] -1.566085569
#> [366,]  1.359236825
#> [367,]  3.151372446
#> [368,]  0.476389442
#> [369,]  1.005453227
#> [370,]  2.829643223
#> [371,]  0.790214601
#> [372,]  3.226018491
#> [373,] -1.277551949
#> [374,]  0.945701319
#> [375,] -1.458637049
#> [376,]  3.211832171
#> [377,]  1.520291119
#> [378,]  2.469372451
#> [379,]  3.314214213
#> [380,]  0.529717646
#> [381,]  1.777806574
#> [382,] -0.682694933
#> [383,] -1.776967879
#> [384,]  0.607481168
#> [385,]  0.224508226
#> [386,] -0.871289988
#> [387,]  3.495693737
#> [388,]  3.200392712
#> [389,]  1.852066934
#> [390,] -0.508182978
#> [391,]  1.275844309
#> [392,] -0.724010613
#> [393,] -1.047605905
#> [394,] -0.829513991
#> [395,]  2.090110986
#> [396,] -0.047100345
#> [397,]  1.481855208
#> [398,]  1.307687637
#> [399,]  0.838726781
#> [400,] -0.087167463
#> [401,]  3.255470211
#> [402,] -2.777933776
#> [403,]  6.383029563
#> [404,]  2.801047438
#> [405,] -2.134949411
#> [406,]  6.513604065
#> [407,]  0.304642892
#> [408,]  2.049552638
#> [409,]  3.385506941
#> [410,]  3.878329253
#> [411,]  1.293982922
#> [412,]  2.568939816
#> [413,] -0.190960042
#> [414,]  0.706341546
#> [415,]  2.268841939
#> [416,]  2.383370754
#> [417,]  0.951089135
#> [418,]  7.635126703
#> [419,]  1.190294203
#> [420,]  1.770766679
#> [421,]  3.904750422
#> [422,] -1.222944131
#> [423,]  6.700886221
#> [424,]  2.385962889
#> [425,] -0.075498263
#> [426,]  0.671031383
#> [427,]  0.219666822
#> [428,] -1.146338369
#> [429,]  6.713182474
#> [430,]  3.476662222
#> [431,] -0.383497497
#> [432,]  3.512410902
#> [433,]  3.559881350
#> [434,]  0.954714897
#> [435,] -1.237225441
#> [436,]  2.410858391
#> [437,] -0.881478164
#> [438,]  3.153495433
#> [439,] -0.432806499
#> [440,]  1.166995909
#> [441,] -1.613625738
#> [442,] -0.849129084
#> [443,]  0.873304312
#> [444,] -0.431444806
#> [445,]  0.734216468
#> [446,]  2.942998698
#> [447,]  1.038101238
#> [448,]  1.849303177
#> [449,] -0.317883027
#> [450,] -1.103281194
#> [451,]  3.710804294
#> [452,] -0.477951900
#> [453,] -2.023696664
#> [454,]  3.251384615
#> [455,]  0.489149113
#> [456,]  0.560477713
#> [457,]  1.190081747
#> [458,]  4.301674530
#> [459,] -0.923495093
#> [460,]  1.249836344
#> [461,]  1.748994814
#> [462,]  2.212405923
#> [463,]  2.603000063
#> [464,]  5.301273793
#> [465,]  3.020182971
#> [466,]  1.950882417
#> [467,]  2.265310822
#> [468,]  0.842267013
#> [469,]  1.135472475
#> [470,] -3.073729731
#> [471,]  2.604805432
#> [472,]  4.340024783
#> [473,]  3.331248118
#> [474,]  3.204966018
#> [475,]  0.651149382
#> [476,]  1.993755119
#> [477,]  3.151012263
#> [478,]  4.950835459
#> [479,] -0.246113209
#> [480,]  4.004805834
#> [481,] -0.053151739
#> [482,]  1.375704136
#> [483,]  0.916894658
#> [484,]  2.607574451
#> [485,] -1.970751630
#> [486,]  0.708545815
#> [487,]  0.719980787
#> [488,]  0.656174664
#> [489,] -2.311831694
#> [490,] -0.972368414
#> [491,]  0.816756192
#> [492,]  2.775805452
#> [493,]  2.974870060
#> [494,]  2.404970258
#> [495,]  2.394190074
#> [496,]  1.776786094
#> [497,] -1.153249230
#> [498,] -1.005477023
#> [499,]  0.317281710
#> [500,]  1.187621439
#> [501,]  0.093384436
#> [502,]  4.642327836
#> [503,]  5.376828137
#> [504,]  1.586796959
#> [505,] -1.415105029
#> [506,]  2.852787216
#> [507,]  1.934040071
#> [508,]  0.997093537
#> [509,]  4.214401909
#> [510,]  2.186835476
#> [511,]  0.409146387
#> [512,] -0.117689813
#> [513,] -0.136817153
#> [514,]  0.865646024
#> [515,]  0.772638365
#> [516,]  1.226632548
#> [517,]  2.204344467
#> [518,]  3.000429693
#> [519,]  1.120091555
#> [520,]  0.594708277
#> attr(,"df")
#> [1] 3 2
#> attr(,"range")
#> [1] 15 36
#> attr(,"lag")
#> [1]  0 85
#> attr(,"argvar")
#> attr(,"argvar")$fun
#> [1] "ns"
#> 
#> attr(,"argvar")$knots
#> [1] 23.80351 26.88520
#> 
#> attr(,"argvar")$intercept
#> [1] FALSE
#> 
#> attr(,"argvar")$Boundary.knots
#> [1] 15 36
#> 
#> attr(,"arglag")
#> attr(,"arglag")$fun
#> [1] "ns"
#> 
#> attr(,"arglag")$knots
#> numeric(0)
#> 
#> attr(,"arglag")$intercept
#> [1] TRUE
#> 
#> attr(,"arglag")$Boundary.knots
#> [1]  0 85
#> 
#> attr(,"class")
#> [1] "crossbasis" "matrix"    
#> 
#> $rain
#>         rain_v1.l1    rain_v1.l2 rain_v2.l1   rain_v2.l2 rain_v3.l1
#>   [1,]  1.61389594  2.631104e-01   8.535782  2.598576967 -3.0641914
#>   [2,]  0.66381000  2.244610e-01   9.317647  0.278047380 -2.3484429
#>   [3,]  1.34430978 -1.610633e-01   6.172373 -0.600013638 -1.5551083
#>   [4,]  1.07171014  9.813093e-02   8.041785  0.137612975 -3.0150779
#>   [5,]  1.14097739  7.234636e-01   6.858040  2.610365307 -2.5548606
#>   [6,]  0.52849201 -5.309371e-01  10.813620  1.716028168 -4.4820798
#>   [7,]  1.55457387  3.236603e-01  11.606426  1.459932350 -3.7490115
#>   [8,]  1.61217991 -3.106883e-01   9.168084  0.152286377 -2.6948152
#>   [9,]  1.99258955  6.198441e-02   9.066006  2.737267417 -2.5791367
#>  [10,]  0.62990786  7.009364e-01   7.636336  2.324968833 -2.7460341
#>  [11,]  1.00300585 -7.715530e-01   8.243855 -1.599695626 -3.0869961
#>  [12,]  0.76323243  7.293650e-01   7.706129  2.097159401 -2.9820286
#>  [13,]  0.16673923  3.809821e-01  10.066491 -0.695607928 -4.1379496
#>  [14,]  0.98750929  4.288169e-01   6.689666  1.308141251 -2.4177214
#>  [15,]  0.84822750  1.820475e-01   6.264129  0.471484705 -1.6417335
#>  [16,]  0.39433130  8.269901e-02   7.018000  1.446899235 -2.0009305
#>  [17,]  0.62349642  6.558574e-01   8.727608  1.455422207 -3.4417011
#>  [18,]  1.03632485  4.050714e-01   6.048743 -0.369525637 -1.8225831
#>  [19,]  0.64896966 -4.101067e-01   8.228547  0.009211242 -3.0871958
#>  [20,]  1.42129396  6.157997e-01   7.494445  2.507374222 -2.7576564
#>  [21,]  0.96930346  1.680071e-01   5.054699  2.059495280 -1.7517050
#>  [22,]  2.36412867  9.685225e-01   9.194740  1.164312581 -3.0771454
#>  [23,]  0.94921208  3.261234e-01   8.795670 -0.796930896 -2.4234610
#>  [24,]  1.01198565 -2.331677e-01   6.664225  1.314166177 -2.2860528
#>  [25,]  1.87495293 -2.310920e-01   6.268657  0.123652709 -1.6927591
#>  [26,]  0.82065740 -1.403732e-01   5.924495  0.999235180 -2.3139791
#>  [27,]  0.28064586  4.531402e-01   4.434266  0.679121520 -1.4337696
#>  [28,]  1.73583070  2.756799e-01  11.235679 -0.874085666 -3.1410319
#>  [29,]  1.13979210  1.362102e-01  10.382990 -0.735268770 -3.3184988
#>  [30,]  0.68706448 -8.031877e-03   8.426742  1.157587859 -3.1430531
#>  [31,]  0.97345839 -2.619748e-01   8.144696 -1.765323724 -2.8133465
#>  [32,]  0.22817292 -3.172326e-01  15.218608  2.542780104 -6.2576835
#>  [33,]  1.85843641 -2.700308e-01   8.158620 -1.168462508 -2.1339602
#>  [34,]  0.65730269  3.532336e-01   6.951718  1.063226810 -2.8550997
#>  [35,]  0.76879123  3.145490e-01   9.803270  0.618218317 -4.0131055
#>  [36,]  0.77092591  3.587540e-01   6.975915  2.152190165 -2.4148847
#>  [37,]  0.11484257 -4.280235e-02   7.739060  0.554791534 -3.2009233
#>  [38,]  1.09673109 -1.695713e-02   8.941663 -1.893835371 -2.7677763
#>  [39,]  1.05449163 -5.426839e-01   8.830709  3.148644992 -3.4456073
#>  [40,] -0.23508699  1.118329e-01   7.536622  0.184750801 -2.6528276
#>  [41,]  0.41841882 -4.066485e-01   6.801038 -0.223431036 -2.6590748
#>  [42,]  1.54207398  2.761138e-01  10.579521  1.114188272 -3.9055524
#>  [43,]  1.00005276  3.033173e-01   5.861766  0.093104274 -2.1588257
#>  [44,]  1.38118370  7.125840e-01   9.677534  0.795184072 -2.6122157
#>  [45,]  2.21314539  7.829270e-01   8.706249 -1.201904872 -1.9235036
#>  [46,] -0.23285140 -2.099412e-01   7.788568 -0.208465130 -2.5074914
#>  [47,] -0.11101350  8.106936e-02   9.637131 -1.607131618 -3.8251921
#>  [48,]  0.25542346  6.445687e-02   7.089200 -0.012092818 -1.9883575
#>  [49,]  1.03985662  8.936534e-01   8.884046 -0.208652451 -3.0332819
#>  [50,]  1.41696958  2.073683e-01   7.712809  1.351326476 -2.6011073
#>  [51,]  1.53293373  2.263615e-01  10.931574  0.907077074 -4.0766168
#>  [52,]  1.91312326  1.147739e+00   8.404044 -0.504006800 -2.3881874
#>  [53,]  2.17672304  4.920955e-01   8.407800 -0.553002165 -2.7944656
#>  [54,] -0.09938209 -4.060658e-01   7.933621  2.099568699 -2.8036261
#>  [55,]  1.95374718  3.120906e-01  11.451923  1.746342249 -3.7575487
#>  [56,]  0.73579798 -1.329973e-01   7.285780 -0.271097612 -2.7760403
#>  [57,]  1.15103427 -3.419503e-01   8.757348  0.871288408 -2.8934164
#>  [58,]  2.74613990  9.802095e-01  10.156525  2.469193670 -3.3051593
#>  [59,]  2.90812754 -2.052931e-01   9.676334 -0.572989197 -2.8896985
#>  [60,]  1.12562745  1.090522e+00   7.726039 -0.335781801 -2.1774115
#>  [61,]  1.20849414  3.578701e-01   9.982734  0.650674655 -2.5292279
#>  [62,]  1.62853799  3.755695e-01   7.919000  0.645685626 -2.4627176
#>  [63,]  1.38692319 -2.052908e-01  11.603538 -1.173451467 -3.8575076
#>  [64,]  1.86219397 -1.243918e+00   8.528913 -1.110197074 -2.2911452
#>  [65,]  0.68560414 -7.116060e-01   5.666426 -0.935499081 -1.7903593
#>  [66,]  0.01327362  9.193131e-01   8.322114 -0.077200180 -3.6963033
#>  [67,]  1.87596145  4.133933e-01   7.403957  1.613175679 -2.1951271
#>  [68,]  1.12359052  7.419587e-01   6.416401 -0.019818917 -2.3255793
#>  [69,]  1.58676780 -3.320972e-01   6.478228 -0.074602720 -1.9692098
#>  [70,]  1.53626217  1.725499e-01  12.833542  1.712628998 -4.1740904
#>  [71,]  1.57140349 -3.404690e-01   8.685400  1.636081167 -3.2138423
#>  [72,]  1.02072928 -1.050873e-01  10.180775  0.085085114 -3.5519067
#>  [73,]  0.53000168  7.747432e-01  11.243992  0.539151791 -4.0287258
#>  [74,]  0.83867217  2.238914e-01   5.959292  0.955147383 -2.4679916
#>  [75,]  1.00090554  6.891278e-01  14.858758  1.398197307 -5.7364203
#>  [76,]  1.55116161  8.921379e-01   8.726932  2.069789867 -2.7310047
#>  [77,]  1.80377260 -2.381268e-01  11.002655  0.302317465 -3.5036520
#>  [78,]  1.42057028 -3.885348e-01   8.825753 -0.283893752 -3.0884598
#>  [79,]  1.52368002  1.336029e+00  10.678668 -1.190578430 -3.8941446
#>  [80,]  0.60301812 -8.909974e-01   9.223577  1.151526321 -3.7324084
#>  [81,]  0.58377887 -1.240440e-01   7.937862  0.551174464 -3.1108975
#>  [82,]  1.00999949  3.735721e-01   7.756431  0.154606330 -2.9303658
#>  [83,]  1.77309622 -5.504337e-01   9.235588  2.362595237 -2.7129297
#>  [84,]  1.04736960  4.924176e-01   9.240703  0.002516658 -3.2988056
#>  [85,]  0.99626822 -4.783180e-01  10.168873  1.704238162 -3.5818338
#>  [86,]  2.64846600  1.912999e-01   8.839562  1.521324775 -2.2399197
#>  [87,]  1.31214176 -3.237450e-01   9.024238  0.878944717 -3.0039032
#>  [88,]  0.12242240  5.900050e-01   3.659170  1.099377360 -0.6724961
#>  [89,]  1.86213249  7.316544e-01  10.132264 -1.000801139 -3.2664065
#>  [90,]  1.64887684 -9.673150e-01  11.780189  2.684977073 -4.0669839
#>  [91,]  0.70033370  7.247003e-02   6.359833  1.966843182 -2.0807996
#>  [92,]  1.67970376 -6.254607e-01   7.397716  0.683030889 -1.6545490
#>  [93,]  1.01821464 -1.082894e+00   8.719060  2.417217149 -2.8633281
#>  [94,]  0.34603020  3.177105e-01   8.302586  0.535508482 -2.8819915
#>  [95,]  1.45261187 -7.362627e-03  14.046202  1.827933873 -3.5891837
#>  [96,]  0.89751181  1.900651e-02   7.418635  0.610728444 -2.7610503
#>  [97,]  0.87495016 -5.819174e-01   5.487962  0.095517748 -1.6307185
#>  [98,]  1.53409403 -6.893321e-02   6.710323  0.845267549 -2.0023743
#>  [99,] -0.34689915 -4.801922e-01   9.325574  2.563034629 -3.5884383
#> [100,]  1.62116368 -1.032242e+00   6.731417  0.062803142 -1.9029819
#> [101,]  0.74932026  6.079337e-01   7.822322  1.199021850 -3.1489631
#> [102,]  0.91004863  1.591242e-01   7.453312  2.252448559 -2.7081349
#> [103,]  1.81405018 -3.306643e-01   7.141610 -0.487604602 -1.6463324
#> [104,]  1.60342261 -3.500471e-01   9.507524  2.606415028 -3.4360096
#> [105,]  0.70356280  8.742565e-02   6.436422  0.487092326 -1.9856272
#> [106,]  1.25974479  3.752352e-01   5.727948  0.016447361 -1.2166372
#> [107,]  0.62623321  4.879701e-01   7.934547  0.206988790 -2.2565823
#> [108,]  0.39580684  6.310726e-02   6.078374  0.582011939 -2.3182977
#> [109,]  1.06973342 -1.407351e-01  10.872470  2.074910919 -3.6588770
#> [110,]  0.28715207 -4.965101e-02   8.089367  2.972415680 -3.3280760
#> [111,]  1.01793949 -7.019077e-01   7.689108  1.020735411 -2.8129019
#> [112,]  2.61483040  4.886516e-01   9.444219 -0.617145955 -2.8363696
#> [113,]  1.21588052 -1.699755e-01   4.852108  1.517217569 -1.2130189
#> [114,]  0.91662406 -6.232206e-01   6.955795 -0.019192610 -2.7951971
#> [115,]  0.75375358 -3.924438e-02   6.185438  2.250423585 -2.4865799
#> [116,]  0.67872839  8.768499e-01   6.926895 -0.647059168 -2.5578655
#> [117,]  0.36587992 -3.085692e-01   4.884134  0.739680834 -1.9785237
#> [118,]  0.49441153  8.846475e-01  10.328667 -0.109492431 -4.1823890
#> [119,]  0.30891711 -7.402094e-01   6.658075  2.867244005 -1.8441250
#> [120,]  0.84199105  3.633688e-01  11.035583 -0.844737737 -3.4968100
#> [121,]  0.86173548  2.523332e-01   3.575372  0.727306476 -1.0400567
#> [122,] -0.37997026  7.340399e-01   6.735477 -0.097391777 -2.6397227
#> [123,]  2.42032201  1.536740e-01  15.841480  2.216936689 -5.9390131
#> [124,]  0.76122759 -7.110155e-01  10.085963  0.132251984 -2.9787098
#> [125,]  0.31193790 -2.141689e-01   2.637892  0.213196526 -0.8757132
#> [126,]  0.83044508  2.385234e-01  10.307660  0.506613602 -3.9315185
#> [127,]  2.36279651 -7.698254e-01  14.117534  0.998999007 -4.3900272
#> [128,] -0.18660327  5.045643e-01   7.186186  0.897542492 -2.6977468
#> [129,]  0.45084952  1.828128e-01   3.234798  0.354725021 -1.1816613
#> [130,]  0.45294304 -1.862809e-01   6.720201  0.381880339 -2.5667035
#> [131,]  2.05791888  6.921392e-01  13.138760  1.439428567 -4.8079958
#> [132,]  1.41048918  6.272614e-02  10.426713 -0.109393813 -3.8898923
#> [133,]  0.55943489 -5.652095e-01   8.228497  0.853033214 -3.1466479
#> [134,]  1.18992806 -7.731652e-02   9.505226  0.143361929 -3.3862456
#> [135,]  1.14294447 -2.319454e-01   7.381366 -0.323458162 -2.5835916
#> [136,]  0.83912919 -1.549771e-01   6.577492 -0.796631104 -1.6143000
#> [137,]  1.23926349  1.801163e-01  11.052908 -0.429344367 -3.8903532
#> [138,]  0.72244241  9.785688e-01   5.106530 -1.015799723 -1.9445585
#> [139,]  0.02992166 -5.910753e-01   5.594714  0.836808345 -1.3220322
#> [140,] -0.42876729 -1.457832e-01   5.211108 -0.466464189 -2.3165259
#> [141,]  2.77591370  6.676812e-01  10.132794  0.260019581 -2.7669000
#> [142,]  0.67507429 -5.200159e-01  11.866010  1.191639071 -4.6100173
#> [143,]  0.33478669 -2.498145e-02   6.474298  0.779129507 -2.6438908
#> [144,] -0.29077737  2.415318e-01   9.223105  1.461795772 -3.9007137
#> [145,]  1.62100993 -7.752260e-01  14.200432  0.967462823 -4.8388483
#> [146,]  1.57494055  6.474924e-01   8.550999  0.955118145 -2.7766616
#> [147,]  2.18391444  4.101000e-01  10.607118 -0.352244127 -3.1853643
#> [148,]  1.07033339  7.510368e-02  10.228139  1.401241139 -3.7190132
#> [149,]  0.57588770  7.393344e-01  12.179880  0.699218737 -4.4501997
#> [150,]  0.24248438 -4.832327e-01   4.424378 -2.011636199 -1.5261634
#> [151,]  0.22414601  6.184691e-01   4.897449  1.177765037 -2.0675112
#> [152,]  0.40085325 -1.044368e-01  11.590494 -1.285608494 -3.6142708
#> [153,]  0.50628991 -6.752442e-01  12.043604  0.229971392 -4.5248039
#> [154,]  0.71241687  8.428693e-01   4.749215 -0.372666653 -1.4371731
#> [155,]  0.77495866  5.682710e-01  10.870130  1.344740555 -4.3134182
#> [156,]  0.39373026 -3.128835e-01   9.097423  1.438546459 -3.2732530
#> [157,]  2.15806274  1.267634e-01  15.081570  1.903627205 -4.1233492
#> [158,]  1.00186966 -5.552855e-01   8.395318 -1.380512298 -2.8002625
#> [159,]  0.94623871  7.525279e-02   5.532246  1.347887323 -2.1912159
#> [160,]  0.74789737  6.506898e-01   6.833790  0.447518050 -1.6228895
#> [161,]  2.14761419 -1.808886e-01   9.702603  0.085124551 -1.4750297
#> [162,]  0.05711187 -6.817896e-01   7.858727  0.935546742 -3.2732331
#> [163,]  1.05295515 -5.700463e-01   6.905186  1.383244713 -2.2570154
#> [164,]  1.42062520  1.604400e+00  10.878036  1.336081274 -4.0488885
#> [165,]  0.30775759 -5.095423e-01   4.962789  3.604957306 -1.9024717
#> [166,]  1.94107768  4.160641e-01  12.917864  2.096636433 -4.9384408
#> [167,]  1.07023068 -5.568019e-01   9.756235 -0.250038455 -3.0231671
#> [168,]  0.28704439 -2.695548e-01  10.496358  0.415409291 -3.8687146
#> [169,]  1.96470293 -2.019122e-01   9.073906  0.663118572 -2.9185113
#> [170,] -0.41379672  2.455457e-01   7.403318  1.329659757 -3.2656775
#> [171,]  1.22236198  2.933358e-01  10.651475  2.578349594 -3.2772900
#> [172,]  0.07365927  1.533844e-01   6.570093  2.251126922 -2.1851057
#> [173,]  0.48924521  1.204354e-01   8.826907  0.823886897 -3.3117605
#> [174,]  1.46127038 -1.259464e-01   8.058069 -1.636358376 -2.1818841
#> [175,]  0.76241947  1.308839e-01   8.353090 -1.184050624 -3.3740603
#> [176,]  1.29192911 -8.123312e-02   8.411763  1.550388765 -2.9630631
#> [177,]  0.46839467  1.326332e-01   6.209774 -1.471944827 -2.3582772
#> [178,]  3.43719577  9.111003e-01  14.039874  1.005814574 -3.7433835
#> [179,]  1.25594463  1.750923e-01  13.112961 -0.435360951 -5.1490917
#> [180,]  1.01365350 -1.040619e+00   7.581569 -0.311833181 -2.1343170
#> [181,]  1.38194696  4.187756e-01   6.978511  1.550085100 -2.0751128
#> [182,]  1.85620727  1.231882e+00  10.162262  1.774385168 -3.1904417
#> [183,]  0.43040612 -1.635746e-01   8.744092  0.737406899 -2.9965006
#> [184,]  1.09776221  8.536928e-01  10.635200  0.568208463 -3.5641817
#> [185,]  0.79944104 -2.055474e-01   7.221119 -1.192643831 -2.3911191
#> [186,]  1.56785285  8.486747e-01  11.295168  0.694543901 -4.4216724
#> [187,]  1.77717080  6.922966e-01  11.966965  2.467354081 -4.1388416
#> [188,]  0.38426617 -4.539695e-01   6.503759  0.443995979 -2.6824547
#> [189,]  0.68320746  9.061865e-01   7.626236  1.364847356 -2.5414582
#> [190,]  1.77547956  3.448004e-01   9.160575 -0.829340474 -2.9604424
#> [191,]  1.37677197  1.008870e+00  11.872801  1.730652090 -3.8880441
#> [192,]  1.03845891  1.108880e-01   6.064897  1.403589991 -2.2410239
#> [193,]  0.67138393  2.223570e-01   6.930809  0.238700291 -2.6138260
#> [194,]  0.72832061 -3.540097e-01   7.508715  1.806805661 -1.3971116
#> [195,]  2.04472116 -5.291750e-01  12.045093 -0.183789673 -4.2791118
#> [196,]  0.57405705  1.964588e-01  10.121307  1.024955105 -4.0177815
#> [197,]  1.85316692 -8.281554e-01  11.360339  1.469823916 -4.0103017
#> [198,]  2.05396978 -6.603863e-01  10.297948  2.245158986 -2.8278429
#> [199,]  1.13839536 -7.444149e-01   4.811078  0.258046888 -1.5324967
#> [200,]  2.37892247 -7.931115e-01  14.551051  1.456002623 -4.2116940
#> [201,]  0.61011834 -4.301255e-01   9.333181  0.051960869 -3.9042116
#> [202,]  0.92719999  1.504945e-01   7.605042 -0.957839846 -2.9918537
#> [203,]  2.04409758  7.006884e-01  12.659720 -0.997552604 -3.8352235
#> [204,]  0.61508187 -3.719685e-01  13.460308  1.153164830 -5.1824573
#> [205,] -0.40820886 -9.723241e-02   4.539834  1.070482136 -1.9630641
#> [206,]  1.56870436  7.281991e-01  11.213168  1.386928797 -3.5765115
#> [207,]  1.34596845  2.020603e-01   6.627725 -0.200276795 -2.2618826
#> [208,]  2.24596796 -1.833659e-01  11.390867  2.263276224 -3.4387923
#> [209,]  1.40060303 -3.679687e-01  11.140460  0.027946650 -3.4606634
#> [210,]  0.97938970  4.244857e-01   8.780449  0.408420407 -3.1228059
#> [211,]  1.34560621  6.624599e-01   6.468746  2.039349401 -2.2833041
#> [212,]  1.27258677  7.332283e-02   8.334366 -1.374982399 -2.0630928
#> [213,]  0.48002153 -3.296284e-01   5.143160 -0.315978092 -1.8245378
#> [214,]  1.03745528 -4.876497e-02   5.842632 -0.331763493 -2.0065308
#> [215,]  1.68478436  6.922751e-02   8.994459  1.847572410 -2.6830028
#> [216,]  0.60526107 -4.556485e-01   7.325831 -2.569041737 -2.1970433
#> [217,]  0.98492383  3.687182e-01   6.186049  0.350076204 -1.5328197
#> [218,]  0.96531202  1.471819e-01   8.162095  0.870040026 -2.5654896
#> [219,]  1.67832841  5.604567e-01   6.831701 -1.032407183 -2.2312306
#> [220,]  1.15537129  9.451675e-01   6.209800  1.742103508 -2.1274219
#> [221,]  1.87467156  4.005098e-01  11.566419  1.131503258 -3.9513453
#> [222,] -0.14348416 -1.692153e-02   3.612819  0.806214410 -1.0121716
#> [223,]  0.29393176  7.160355e-01   9.378026  1.334865661 -3.2531656
#> [224,]  1.21064694 -8.591564e-02   4.438098 -0.011973468 -1.3087174
#> [225,]  0.46812805 -7.608544e-01   6.799569 -0.017184330 -2.6913625
#> [226,] -0.65359166  3.668761e-02   9.548720 -0.567510679 -4.0874169
#> [227,]  1.01093422 -3.790122e-01  10.678945  0.339060290 -4.2828255
#> [228,]  1.05600078 -8.783265e-01   7.997296  0.119335038 -2.7476544
#> [229,]  0.62292887  4.377444e-01   5.007759  0.061538328 -1.6789438
#> [230,]  1.74276512  1.113689e+00   8.566327  0.603502332 -2.9463762
#> [231,]  0.57344131  2.805595e-01   6.689231  1.017920509 -2.2305298
#> [232,]  0.97142758 -4.769870e-01   8.667672  1.604698045 -3.3828826
#> [233,]  0.05556516 -7.017386e-01   5.659467  0.910284527 -2.3543566
#> [234,]  0.41375396 -1.142725e-01   8.323610  1.160718506 -3.4867289
#> [235,]  1.13321563  3.902313e-01   5.247184  0.676083746 -1.3513023
#> [236,]  0.78227190  4.839666e-01   8.515709 -1.019298386 -3.1362207
#> [237,] -0.17716446 -8.143400e-01   5.581438 -0.175224453 -2.2460814
#> [238,]  1.32929031  3.111593e-01  12.483440  0.954457523 -3.5451635
#> [239,]  0.57536003 -2.375402e-02   3.706848  1.739688062 -1.3542697
#> [240,]  0.42775157 -2.667009e-01   5.463262 -2.052185410 -1.5145050
#> [241,]  1.67713832  7.906781e-01  11.200439  2.484259276 -3.7200630
#> [242,] -0.05886693  1.047247e+00  10.846015 -0.761838187 -4.5677246
#> [243,]  0.79394310  6.326967e-02  10.929857  1.803871399 -4.3500399
#> [244,]  0.23125047  4.315033e-01   7.429581  0.455976136 -2.9148082
#> [245,]  0.54598101 -2.273451e-01   5.031964  0.815667263 -1.7878367
#> [246,]  0.87296331 -3.752401e-01   9.707247  0.057897640 -2.7769366
#> [247,]  0.58699075  6.219344e-01   7.727907  1.907144380 -2.8414069
#> [248,]  0.98048372  5.453175e-02  10.605986  1.696496006 -3.5372110
#> [249,]  0.99285173 -3.347515e-01   5.873830  1.016100995 -1.4199457
#> [250,] -0.19406386  6.974588e-01   7.157422  2.596259937 -3.0320171
#> [251,]  1.25469712  4.880697e-01   9.910442  0.985817363 -3.3774268
#> [252,]  2.00898349  2.668201e-01  12.314731  0.020191776 -4.4060794
#> [253,]  0.59768226  3.288627e-01   7.315325  2.122871021 -2.4540407
#> [254,]  0.30699888 -1.564161e-01   4.582720  1.762820371 -1.2661721
#> [255,]  0.23359335 -3.686488e-02   7.123536  0.101134557 -1.9790589
#> [256,]  1.92943542  2.430579e-01   7.430945  1.181267704 -2.1979245
#> [257,]  0.64661782  9.319584e-01   9.042913  0.412195132 -3.2262759
#> [258,]  1.56022781  2.831061e-01   6.133412 -0.083815671 -1.8124674
#> [259,]  1.96158441  1.228147e+00  10.669358  0.902283137 -3.8254891
#> [260,]  0.94682284  2.049542e-01   4.435166  0.779746165 -1.2434425
#> [261,]  0.53162297 -2.366974e-01   8.195832  0.912696173 -3.0579169
#> [262,]  0.48150748  4.793964e-01   7.163952  1.775411990 -2.6339500
#> [263,]  0.93118520  1.313594e-02   4.637890  1.633811508 -1.0985892
#> [264,] -0.04447794  9.607789e-01   8.824664  0.657875010 -3.7859489
#> [265,]  1.65986610 -2.290367e-01   8.213511  1.384078527 -2.6351502
#> [266,] -0.22691492 -2.020932e-01   7.715008  0.509470990 -3.3090103
#> [267,]  0.59061493 -4.568014e-01   4.042523  0.317178244 -1.1800545
#> [268,]  0.60520787 -3.856363e-01   7.114444  0.514644649 -2.4883024
#> [269,]  1.21156868  3.044878e-01  11.366361 -0.903250288 -3.8192592
#> [270,]  0.36187418  6.294421e-01   7.363100  0.698132657 -2.9195236
#> [271,]  1.07644373  2.846759e-01   7.752697  0.255349238 -3.0157483
#> [272,]  1.43370948  5.110590e-02  10.955729 -0.468604669 -3.7539248
#> [273,]  0.84805742  9.206490e-01   7.519994  0.245092259 -2.9713778
#> [274,]  2.41422809  4.507182e-01   7.956498  0.777143451 -1.7089980
#> [275,]  1.18462367 -9.262499e-02  10.483827 -0.179103396 -3.7623587
#> [276,]  1.62597977  1.359175e-01  11.773243  0.840428225 -4.2274840
#> [277,]  0.99403244 -6.418161e-02   9.231719  0.310472452 -3.1901628
#> [278,]  2.64858192 -3.163301e-01  12.470192 -0.092777940 -3.6690807
#> [279,]  0.23153136 -2.343747e-01   6.388769  2.331443538 -2.5761997
#> [280,]  0.49297404  4.627429e-01   7.571705 -0.559592834 -2.6025717
#> [281,]  1.94668268  3.031774e-01   7.388389  1.154742877 -2.0917810
#> [282,]  0.75970666  4.941492e-01   9.595023  2.674368767 -3.8071762
#> [283,]  0.47968123  2.127912e-02   7.115266 -0.383512708 -1.9199884
#> [284,]  1.40632756  1.187265e-02  10.175012  2.432642403 -3.6225606
#> [285,]  1.21060816 -4.262649e-01   5.067901 -1.547800875 -1.1730873
#> [286,]  1.22829005  1.811093e-01  11.941388  1.152514168 -4.4964979
#> [287,]  0.40950999 -2.963586e-01   5.320354  0.650339656 -1.6457470
#> [288,]  1.27516737  1.811612e-01  11.970338  1.133157588 -4.5241082
#> [289,]  1.51202639 -8.474553e-01  14.782904  0.005928183 -5.4008290
#> [290,]  1.57707131 -2.152382e-02  11.783729  1.422029652 -3.4009686
#> [291,]  0.87756854  6.820717e-01   6.312398  2.144128021 -2.2448453
#> [292,]  0.47882195  6.746915e-01   8.351265  1.967915203 -2.1635132
#> [293,]  0.33149858 -7.859155e-02   1.737903  0.199803670 -0.6369908
#> [294,]  1.61046872  1.238443e+00   8.618963  1.112827531 -3.1069579
#> [295,]  1.04317104 -2.712755e-01   6.080981  0.280819212 -2.0044976
#> [296,]  0.74526625 -4.280751e-01   7.616746 -0.713796084 -2.8276115
#> [297,]  2.07794761 -4.419822e-02   9.882549  1.052146873 -2.6364415
#> [298,]  0.89270210  4.763962e-01   6.135523  1.838867730 -2.4289140
#> [299,]  0.57344203 -3.119229e-01   7.895141  0.416480202 -2.8907841
#> [300,] -0.19500335  1.083910e+00   7.438773 -0.831684808 -3.2801135
#> [301,]  1.43359689  3.205327e-02   7.705790  1.681293709 -2.6474736
#> [302,]  1.43109421  3.850758e-01  13.244101 -0.986732392 -4.7774152
#> [303,]  0.66238334 -1.238213e-01   9.810036  2.578536566 -3.4639432
#> [304,]  1.20218217  4.655512e-01   4.982395  1.110461142 -0.6953820
#> [305,]  0.25318342 -9.250878e-01  10.319418  2.566707608 -3.6471973
#> [306,]  0.91438719  5.253523e-02   7.716017  0.299821938 -2.7785467
#> [307,]  0.83763585  3.996261e-01   7.451526  0.756737881 -2.8018687
#> [308,]  1.13961311 -9.700123e-01  11.636858  2.597432294 -3.5342745
#> [309,]  1.06873782  8.474282e-01   6.611498  2.906985482 -1.7445662
#> [310,]  1.28451095  4.161292e-01  13.346347  1.935156799 -5.1913661
#> [311,]  1.41723633 -3.441697e-01   7.373979 -0.338826520 -2.8457931
#> [312,]  0.10174886  2.908855e-01   6.835448  1.422049507 -2.8646996
#> [313,]  1.77257638  2.010786e-01   9.588275 -0.107113460 -3.3254120
#> [314,]  1.35188707  7.937523e-01   9.652571 -0.133300024 -3.1886043
#> [315,]  0.41747864  5.024073e-01   8.349916  2.208424935 -2.9389458
#> [316,]  2.36201001  1.038076e+00  11.209547  1.616564816 -3.4553789
#> [317,]  1.68154420  7.860201e-02  11.950484 -1.114781634 -4.3506613
#> [318,]  1.87086322  4.666186e-01  10.296195  1.099364540 -2.1193826
#> [319,]  1.16393050  1.459413e-01   7.538896  0.833982729 -2.6442955
#> [320,]  0.28624871 -5.499857e-01   7.946351  0.914371243 -2.8672287
#> [321,]  1.05990549  7.217754e-01   8.610590 -1.151517179 -2.8388311
#> [322,]  1.67025445  9.694436e-02   9.840816  1.590543048 -3.5111614
#> [323,]  1.96922841  2.546662e-01  10.814770  0.842524946 -3.2778304
#> [324,]  1.17364945  7.058845e-01  11.507099  1.470048422 -3.9776991
#> [325,]  2.81713651  9.500505e-02   8.971930 -0.675469873 -2.3669389
#> [326,]  0.34850268 -5.543866e-01  13.788451  0.999301750 -5.7807858
#> [327,]  0.30507672 -4.225715e-01   6.511637  1.308245967 -2.2083126
#> [328,]  0.69134632  5.925187e-02   6.842581  0.964426446 -2.4650038
#> [329,]  0.40469293 -1.107497e-02   9.705767  2.508931677 -3.7501568
#> [330,]  1.53436793 -6.799390e-02   9.339475  1.816313152 -2.7390260
#> [331,]  1.39567664  7.786354e-01  13.260030  1.546499556 -4.5402533
#> [332,]  1.44190565  1.031798e+00  10.332621  1.174635284 -3.2752269
#> [333,]  0.21546826 -4.522147e-01   8.131294  0.271442825 -3.2179176
#> [334,]  0.76595688  7.639316e-01   5.654922  1.147856752 -1.6358451
#> [335,]  1.27581023 -1.818718e-01  10.445429  2.776491354 -3.5056261
#> [336,]  0.68422239 -7.959311e-01   7.869637  0.810178864 -3.0574634
#> [337,]  0.72119835  1.910366e-02   6.031306  2.759221980 -2.0966942
#> [338,]  2.15783402  9.786693e-02  11.811307 -0.478685487 -3.6824556
#> [339,]  1.35819166  1.884972e-01   9.050565  0.592792007 -3.1256335
#> [340,]  0.87778350 -2.687657e-02   5.532610 -0.200672338 -1.7700278
#> [341,]  1.00532414 -2.871897e-02   9.801568  0.983570344 -3.6317625
#> [342,]  0.75877583 -2.971131e-01   4.348865 -0.992153556 -1.1845878
#> [343,]  1.75866052  1.934766e-01  10.416143  0.537981821 -3.0847002
#> [344,]  0.99995291 -2.046507e-01   8.787563 -0.809088686 -3.3170168
#> [345,]  1.44506880  4.632794e-01   8.285343  2.181623032 -3.0609213
#> [346,]  0.74312217  2.534767e-01   9.121914  0.742840855 -2.2284982
#> [347,]  0.49287271  3.583616e-04   7.432731 -0.019897202 -3.0579119
#> [348,]  1.41289352  2.299152e-01   6.946196  0.211957696 -2.5331178
#> [349,]  1.48618648  1.210987e-01   6.760299  0.874408914 -2.1670113
#> [350,]  0.69485102 -2.937339e-01  11.002034 -0.211292664 -4.0511014
#> [351,]  1.19750644  8.000828e-01  11.509173 -0.163467339 -3.6741655
#> [352,]  0.39612759 -1.445775e-01   7.090662  1.794751999 -2.2903858
#> [353,]  1.58176398  1.254186e-01  11.169673  1.587928962 -4.1536836
#> [354,]  1.38982974 -3.682295e-02   6.892075  1.643071061 -1.7676466
#> [355,]  1.57177186  7.192168e-01   6.872826 -0.955845412 -2.3159191
#> [356,]  0.85753543  7.325952e-03   4.358680 -0.211143689 -1.5094830
#> [357,]  0.25197387  1.826303e-01   9.627525  1.191268516 -3.5353436
#> [358,] -0.71881286  5.659530e-01   7.901381  0.736870468 -3.5100623
#> [359,]  0.14637054 -8.680576e-02  12.530101  2.476710556 -4.5442948
#> [360,]  1.39858075 -1.459438e-01   7.018404  1.117874928 -2.0615489
#> [361,]  0.75845306  4.259280e-01   4.604433 -0.669089698 -1.3711197
#> [362,]  0.80149030  6.802702e-01   6.500670  0.017433276 -2.2666697
#> [363,] -0.10593979 -1.118957e-01   7.501522  1.196878715 -2.9139080
#> [364,]  1.59702892 -1.853536e-01   9.175472  2.769692006 -3.3148520
#> [365,]  0.45159812  4.362832e-01   5.427344  2.827734876 -1.7351347
#> [366,]  1.24880845  3.200886e-01   7.883585  0.800497293 -1.9420739
#> [367,]  1.69538201  4.447592e-01  10.395640  2.573090589 -3.7998208
#> [368,]  1.21045451 -7.346079e-01   9.770412 -0.377975662 -3.1973764
#> [369,]  0.68080670 -4.411930e-01   8.960272  1.229751637 -3.0470464
#> [370,]  0.45053562  5.234585e-01   4.161340  0.732969796 -0.7410514
#> [371,]  2.09571715  4.087508e-01  10.883641  0.046342991 -3.4585605
#> [372,]  0.50294088  2.837938e-01   7.460598  4.124109356 -2.2909091
#> [373,]  0.90718909 -5.955900e-01  10.623973  0.176797121 -4.0732609
#> [374,]  0.63570703 -8.733054e-02  10.185701 -0.189453629 -3.9929841
#> [375,]  1.22132233  4.965863e-01  11.462978  1.892317050 -3.9680342
#> [376,]  0.58293247 -2.009650e-01   4.309350 -1.031814381 -0.8759877
#> [377,]  1.71025586  7.200502e-02   9.796041  2.441258156 -2.7590654
#> [378,]  1.70882261  7.940383e-01  13.607556  1.317136215 -4.1739335
#> [379,]  2.12726252 -6.880980e-01  11.289561 -0.312396549 -3.6922584
#> [380,]  1.60713051  8.940365e-01   9.699241  1.069063125 -2.9844972
#> [381,]  0.02296055 -1.663825e-02   9.075635  2.441606551 -3.1422480
#> [382,]  0.54944655  5.119333e-01   4.762282  0.106162648 -1.7318419
#> [383,]  0.71609581  4.726346e-01   6.174492  0.321055875 -2.0830774
#> [384,] -0.39864430 -2.483748e-01   2.891632 -0.391144867 -0.8153553
#> [385,]  0.91285764  1.225082e-01   5.366500  0.073874896 -1.3775575
#> [386,]  1.15382674 -3.658521e-01   9.908629 -0.736383202 -3.5396372
#> [387,]  1.39521667 -5.510531e-01  13.942440  0.971271925 -4.6317327
#> [388,]  0.72193712  6.177870e-01   6.060647 -0.010275498 -2.5107399
#> [389,]  3.20680974  2.232577e-01  10.592269  0.664511957 -2.9417512
#> [390,]  0.77129013  6.095408e-01  10.433095 -0.389573564 -3.2272631
#> [391,]  0.37093770 -1.096360e-01   8.382474  1.833315076 -2.8306769
#> [392,]  2.21185254  1.371765e-01  10.110774  1.060055823 -3.0509019
#> [393,]  1.99335536 -1.373826e-01   9.844948  0.925151928 -3.6956000
#> [394,]  0.33568107 -2.263940e-01   5.441522  0.215941071 -1.8650593
#> [395,] -0.67916221  3.867982e-01  10.923793 -0.057183556 -4.2012035
#> [396,]  2.09807985 -6.791459e-01   9.414709  3.082743155 -2.8767992
#> [397,]  1.10616905  5.436060e-01  12.273288  1.211069376 -4.5250291
#> [398,]  2.33418563 -2.262612e-01   9.375420  0.896909264 -2.5560984
#> [399,]  0.98156219 -1.254665e-01   8.980151  1.380717835 -1.9881825
#> [400,]  1.20066789 -4.998299e-01  10.451410  1.720820151 -3.3287485
#> [401,]  0.61064280  5.313276e-01   7.894011 -0.960771923 -2.7427832
#> [402,]  2.03318512 -5.335950e-01   9.169984  0.322706807 -2.6055557
#> [403,]  0.37785663  1.675804e-01   5.336303 -1.367558528 -2.1497676
#> [404,]  1.28931316 -4.471247e-01   8.053416  0.137186217 -1.4326881
#> [405,]  2.60845062  7.347885e-01   9.707472  2.624580812 -1.9652801
#> [406,]  0.73044973 -6.735287e-01  10.469410 -0.021950867 -4.1976229
#> [407,]  0.65513572 -9.689238e-01   9.268933  1.840174437 -3.1973643
#> [408,]  0.57157078 -1.488147e-01  11.727046  3.410832961 -4.3317227
#> [409,]  1.00649909  3.736259e-02   8.344349  0.371696426 -3.0659200
#> [410,]  0.42742240 -8.166046e-02   9.219964 -0.533867051 -3.3021853
#> [411,]  0.55487104  6.311269e-01   8.520908  1.350283491 -2.6187582
#> [412,]  0.65771663 -7.419739e-01   8.203736  1.088304066 -3.2202722
#> [413,]  0.73603326  1.547272e-01   6.970501  0.485535420 -2.0204517
#> [414,]  1.64487686 -1.105997e+00  10.517413  0.707181251 -3.6887420
#> [415,]  0.57629260  7.130674e-02  11.994811 -1.589297797 -4.1570604
#> [416,]  2.07590686 -7.191393e-01   8.286847 -0.214567434 -2.0314452
#> [417,]  1.15564976 -3.983431e-01   5.207820 -1.026195825 -1.6654800
#> [418,] -0.20725421  7.849754e-01   7.194948 -0.615862435 -2.8466678
#> [419,]  0.94639691  4.794755e-01   6.597617  1.574876647 -2.1716111
#> [420,]  1.56075217 -1.079509e-01  12.170736  2.177234015 -4.7421342
#> [421,]  0.93589439  5.231794e-01   7.802144 -0.086372492 -2.9675671
#> [422,]  1.77806584 -1.699430e-01   6.823188 -1.234852024 -1.3307455
#> [423,]  1.45863133 -8.745505e-01  11.711475  1.661621013 -3.6140603
#> [424,]  2.21447979  4.019009e-01   8.021623  0.754908232 -2.5473335
#> [425,]  1.88635700  4.350122e-02   9.190484  0.879759442 -3.2677638
#> [426,]  1.80932024 -3.985106e-01  10.317055  0.616029207 -3.1782912
#> [427,]  0.94155305  1.417316e-01   6.532707  1.239034058 -1.8303436
#> [428,]  0.88152182  2.354754e-01   7.826986  0.393020007 -3.0562326
#> [429,]  0.60858186 -2.237041e-01   6.346406  1.117078009 -2.4686316
#> [430,]  1.34944876 -8.148679e-01  11.945805  0.335529172 -4.2188019
#> [431,]  0.89399964  8.390119e-02   5.814954  0.775326123 -1.3863935
#> [432,]  1.35923681  5.754708e-01   9.852649  1.346203608 -2.5816562
#> [433,]  2.71699339  8.193906e-01  11.643601  1.620108240 -3.4576465
#> [434,]  1.93359733  5.967052e-01   9.201045  0.790485934 -2.5427114
#> [435,] -0.81587265  6.369568e-01   8.395886 -1.554368335 -3.5949789
#> [436,]  1.21688525  2.284390e-01  10.781103  0.641861920 -4.0920552
#> [437,] -0.21666953 -3.471311e-01  10.497094  3.039758669 -4.0972432
#> [438,]  1.06027250  3.710734e-02   7.288109  1.263223192 -2.6004526
#> [439,]  0.92024628  4.704327e-01   9.780747 -1.018578144 -3.7768652
#> [440,]  0.51136873 -8.956713e-02   3.815447 -0.802180328 -1.0294564
#> [441,]  1.58656190 -4.568527e-01   9.961643  0.354633363 -3.5003535
#> [442,]  1.40692527 -9.656842e-02   7.090291  0.972868481 -2.2545683
#> [443,]  1.34458782  5.929538e-01  11.734556  0.782635119 -4.2940980
#> [444,] -0.07232947  4.511600e-01   7.642309  0.576131080 -3.3092620
#> [445,]  1.38477225  1.314011e-02   6.615456  0.312065116 -2.3392504
#> [446,] -0.03273125 -4.876177e-01   5.635796  1.548436161 -2.4099197
#> [447,]  1.31029714  6.794452e-01   6.781268  1.874732652 -2.3159799
#> [448,]  1.67063247  1.566857e-01  12.466718 -0.184974157 -4.5686529
#> [449,]  1.85626645 -4.613177e-01  10.290364 -0.129697596 -2.9714776
#> [450,] -0.04784599  8.985820e-02  13.540122  2.803590607 -5.8370631
#> [451,]  1.31134723 -5.489108e-01  11.901209  2.061969144 -4.4924801
#> [452,]  0.94125008  3.582869e-01   6.150194  0.710186529 -2.2814494
#> [453,]  1.68627302  9.975770e-01   6.778458  3.422894208 -2.2320415
#> [454,]  1.21894629  5.274041e-01  10.476697  2.702200714 -3.6942569
#> [455,]  0.79625647  1.171686e-01   5.533446  0.917135072 -1.9516466
#> [456,]  1.57014070  3.436476e-01   7.650579 -0.023033231 -2.3268497
#> [457,]  0.29082529  9.463458e-02   6.365171  1.231609535 -2.3864878
#> [458,]  1.18175100  1.081261e-01   3.862783  1.411085670 -0.6717798
#> [459,]  1.32792408  9.418318e-02  10.328096  0.755912093 -3.6247156
#> [460,]  0.11568460  3.346669e-01  12.112420 -1.772880528 -4.8952385
#> [461,]  1.72252238  2.391364e-01   7.335728 -0.812652296 -1.2740551
#> [462,]  0.71381364  2.108312e-01   7.903119  2.208339904 -2.6841361
#> [463,]  0.17275546  6.079823e-01   5.274767  0.291302240 -2.1459323
#> [464,]  1.24761914 -2.994294e-01   8.867902  1.251468314 -3.3737535
#> [465,]  1.37312730  4.569573e-01  12.285867 -0.146240423 -4.5035018
#> [466,]  1.62959963 -7.451419e-01   6.565652 -0.737516332 -2.1660013
#> [467,]  0.67334400  6.836608e-01  12.699553  3.834743334 -5.1095026
#> [468,]  1.91687107 -6.119177e-01  10.694551  1.277470228 -3.7130062
#> [469,]  0.45313107 -1.586503e-01   6.938668  0.800079965 -2.5226629
#> [470,]  0.21772258  5.067253e-01   5.887641 -0.344606487 -2.2537754
#> [471,]  0.24229591  6.435065e-01   3.060611 -0.546317997 -1.1744788
#> [472,]  1.03909527 -2.128129e-01   5.339943  1.295739811 -1.6108318
#> [473,]  2.23192003  7.033952e-01  10.666261  1.782714799 -2.5542881
#> [474,]  1.95807220 -1.213504e-03  11.545162  0.824226036 -3.7751194
#> [475,]  0.92217121 -5.017716e-01   6.755582 -1.064427607 -2.0053871
#> [476,]  2.07944937  2.813764e-01   7.757836  2.137003051 -1.8814220
#> [477,]  1.39040974  1.013171e+00  11.178588  1.752410027 -4.1168870
#> [478,]  1.36116438 -1.843161e-01   9.913733 -0.385958920 -3.0125992
#> [479,]  1.07484669  5.125855e-01  10.762121  1.726736568 -4.1340141
#> [480,]  0.60914860  6.926645e-02   6.701696  1.165706707 -2.0234979
#> [481,]  2.01921641 -8.411014e-01   7.387036  0.588452372 -1.8973485
#> [482,]  0.97591444 -5.256645e-01   8.032173  0.070427452 -2.7987207
#> [483,]  0.78170271  1.252882e-02   7.933950  0.146742177 -3.2319445
#> [484,]  1.94687447 -7.076624e-01   9.389818  1.207093926 -3.0671255
#> [485,]  0.81752218  9.869771e-01  10.219698  1.705316516 -3.9664829
#> [486,]  1.85651006 -2.652654e-02   7.024862  0.871356545 -2.4855091
#> [487,]  0.43410588  7.040177e-01   6.395327 -1.114104334 -2.2205446
#> [488,]  1.11325372  3.371561e-01   9.083264  1.760399290 -3.2435338
#> [489,]  0.14327019 -2.612760e-01   4.975057  2.276354305 -1.5648624
#> [490,]  0.38793982 -5.344362e-01  10.003114 -0.432608216 -3.9453830
#> [491,]  0.15426592 -1.268304e+00   9.336141  0.025481652 -3.4510854
#> [492,]  0.43591661  5.174832e-01   5.391203  0.803419350 -2.0122028
#> [493,]  1.12485490 -1.325833e-01  13.375666  1.040356003 -4.6821431
#> [494,]  1.21691259 -4.157245e-01   8.142717  0.844160034 -1.9064774
#> [495,]  1.37692485  1.009355e+00  11.057601  0.513549415 -3.5093320
#> [496,]  1.75740967 -1.213261e+00  12.599097  1.391035105 -4.2758902
#> [497,]  0.96243204  4.045642e-01   6.008584  1.180535517 -2.0445352
#> [498,]  0.37044894 -9.690454e-01   6.123074 -0.689065248 -2.0310766
#> [499,]  1.32813549  5.826978e-01   5.493331  2.092809357 -1.8063241
#> [500,]  0.95240391  5.240030e-01   5.375696  1.311254926 -1.7500269
#> [501,]  1.09844564  7.027087e-01  12.164127  0.657162780 -4.1190493
#> [502,]  1.90740523  2.503732e-01   6.849124  2.899733985 -1.6333624
#> [503,]  1.14508817  3.457508e-01  10.216329  1.742342256 -3.3889118
#> [504,] -0.24516072  3.074796e-02   7.020373  0.144793033 -2.7902922
#> [505,]  0.56178669  8.445674e-01  10.970340  0.349045475 -3.8751587
#> [506,]  1.27903454  2.770675e-01   8.593485 -0.033016833 -3.1825334
#> [507,]  1.22113251  3.980791e-05  11.061909 -0.286465966 -3.8820271
#> [508,]  1.51856214  3.604312e-01   5.310311  1.111920094 -1.2096723
#> [509,] -0.55298066 -4.241228e-01   9.689785  1.290204479 -4.1270815
#> [510,]  0.32107518 -5.102280e-01   3.552676  1.859755200 -1.4397672
#> [511,]  0.61422095 -5.138276e-03   9.822202  0.719203012 -3.9103363
#> [512,]  1.17120903  2.414479e-01   7.850299  2.116049695 -2.6489659
#> [513,]  0.99872706 -1.944096e-02   4.343183  1.881524584 -1.2479248
#> [514,]  2.06376984  5.530839e-01   8.621882  1.106767539 -2.6369969
#> [515,]  1.87623150  7.043383e-01   7.768200  1.419427277 -2.4713343
#> [516,]  0.07328367 -2.757490e-01   4.116329  0.476725759 -1.8129108
#> [517,]  1.16219548  6.964812e-01   5.279605  0.070409458 -1.4170461
#> [518,]  0.73055588  7.464099e-01   7.020243  1.780516114 -2.5878330
#> [519,]  1.11979492  1.139145e+00  10.963575  1.270670594 -3.9178134
#> [520,]  0.28428552  7.518946e-01  14.265761  0.219001503 -5.1302780
#>           rain_v3.l2
#>   [1,] -1.0146053000
#>   [2,]  0.2767089470
#>   [3,]  0.0850995564
#>   [4,]  0.1120852161
#>   [5,] -0.9302672242
#>   [6,] -1.0365543702
#>   [7,] -0.5362411752
#>   [8,] -0.6165910018
#>   [9,] -0.5964281405
#>  [10,] -0.7495384165
#>  [11,]  0.2952896467
#>  [12,] -0.7892190367
#>  [13,]  0.3522360164
#>  [14,] -0.5336659374
#>  [15,]  0.1641858724
#>  [16,] -0.3996320205
#>  [17,] -0.4073175055
#>  [18,]  0.1386494869
#>  [19,]  0.0042274110
#>  [20,] -0.9365873107
#>  [21,] -0.6759766031
#>  [22,] -0.0867017383
#>  [23,] -0.2711313385
#>  [24,] -0.8496831282
#>  [25,] -0.1056925138
#>  [26,] -0.5579365524
#>  [27,] -0.0273732330
#>  [28,]  0.5626476770
#>  [29,]  0.3496240870
#>  [30,] -0.7475440830
#>  [31,]  1.0252459737
#>  [32,] -1.2150001573
#>  [33,]  0.4367474213
#>  [34,] -0.2330010891
#>  [35,] -0.2265476679
#>  [36,] -0.8026963215
#>  [37,] -0.1942544079
#>  [38,]  0.2368592594
#>  [39,] -1.5342209675
#>  [40,] -0.5598808086
#>  [41,]  0.1383521954
#>  [42,] -0.1958281846
#>  [43,]  0.0377824191
#>  [44,] -0.1567284980
#>  [45,]  0.2802016213
#>  [46,]  0.0474969145
#>  [47,]  0.7213211334
#>  [48,]  0.4865770790
#>  [49,]  0.3872001454
#>  [50,] -0.5212172036
#>  [51,] -0.2259143116
#>  [52,]  0.5924883593
#>  [53,]  0.5138656262
#>  [54,] -0.2544866029
#>  [55,] -0.3411148430
#>  [56,]  0.1367236283
#>  [57,]  0.0504359420
#>  [58,] -0.6869516463
#>  [59,] -0.0506906783
#>  [60,]  1.4256874828
#>  [61,]  0.5710579386
#>  [62,] -0.2156451372
#>  [63,]  0.5195770529
#>  [64,] -0.5309922160
#>  [65,]  0.4994415267
#>  [66,]  0.2450399544
#>  [67,] -0.4156254096
#>  [68,]  0.3229996639
#>  [69,] -0.4341384905
#>  [70,] -0.5231219384
#>  [71,] -0.8235439788
#>  [72,] -0.3973577314
#>  [73,]  0.2827974473
#>  [74,] -0.4398034608
#>  [75,]  0.1636623023
#>  [76,] -0.1138924933
#>  [77,] -0.5339962701
#>  [78,]  0.0635015327
#>  [79,]  1.0082309389
#>  [80,] -0.9629815661
#>  [81,] -0.2513025323
#>  [82,] -0.1119343773
#>  [83,] -0.9609247514
#>  [84,]  0.1405458559
#>  [85,] -0.8451516420
#>  [86,] -0.6793647328
#>  [87,] -0.3126418384
#>  [88,] -0.2165381507
#>  [89,]  0.5466703904
#>  [90,] -1.5765245253
#>  [91,] -0.6910501348
#>  [92,] -0.5296054981
#>  [93,] -1.2126822296
#>  [94,]  0.1264874554
#>  [95,] -0.5600534508
#>  [96,] -0.2768391794
#>  [97,] -0.2845016907
#>  [98,] -0.7923556355
#>  [99,] -1.6354068061
#> [100,] -0.5561048968
#> [101,] -0.3327434808
#> [102,] -0.5531127125
#> [103,]  0.2589276874
#> [104,] -1.3534319088
#> [105,] -0.3224188099
#> [106,]  0.9369350216
#> [107,] -0.2013963388
#> [108,] -0.0310659834
#> [109,] -0.6955865100
#> [110,] -1.2560514878
#> [111,] -0.6703981696
#> [112,]  0.7529307606
#> [113,] -0.6360454416
#> [114,] -0.0712863041
#> [115,] -0.9768407970
#> [116,]  0.4992425604
#> [117,] -0.4170408767
#> [118,]  0.3308477869
#> [119,] -1.5949537659
#> [120,]  0.2565585732
#> [121,] -0.1958168741
#> [122,]  0.3832932710
#> [123,] -0.9875047063
#> [124,] -0.1223505477
#> [125,]  0.0068491431
#> [126,]  0.1981335727
#> [127,] -1.4545709046
#> [128,] -0.3002885718
#> [129,] -0.1773515992
#> [130,] -0.3097327641
#> [131,]  0.0178582393
#> [132,]  0.5186994406
#> [133,] -0.5818077264
#> [134,]  0.0606226793
#> [135,] -0.1632448062
#> [136,]  0.5688568296
#> [137,]  0.3554161448
#> [138,]  0.8661838331
#> [139,]  0.0781434648
#> [140,]  0.0858382816
#> [141,]  0.2919463804
#> [142,] -0.5169562323
#> [143,] -0.3924477884
#> [144,] -0.3607259191
#> [145,] -0.3684953837
#> [146,] -0.0008561825
#> [147,]  0.3207925516
#> [148,] -0.4258088176
#> [149,]  0.2530274745
#> [150,]  0.6562402698
#> [151,] -0.3252963399
#> [152,]  0.4074267688
#> [153,] -0.4577673619
#> [154,]  0.8343175856
#> [155,] -0.1543093143
#> [156,] -0.6403809573
#> [157,] -0.5785010503
#> [158,]  0.0049056124
#> [159,] -0.4373276112
#> [160,]  0.0920879092
#> [161,] -0.5358039185
#> [162,] -0.6218148013
#> [163,] -0.2678267463
#> [164,]  0.0402687459
#> [165,] -1.6806106694
#> [166,] -0.7525466315
#> [167,] -0.3251794084
#> [168,] -0.2666635000
#> [169,] -0.3009949413
#> [170,] -0.4535092164
#> [171,] -0.6817792336
#> [172,] -0.6093874977
#> [173,] -0.0817518341
#> [174,]  0.5633118368
#> [175,]  0.4814931214
#> [176,] -0.4979755254
#> [177,]  0.7511163300
#> [178,] -0.4072916492
#> [179,]  0.3511021117
#> [180,] -0.4556622222
#> [181,] -0.3608687720
#> [182,]  0.0891176803
#> [183,] -0.1141585932
#> [184,]  0.0855129775
#> [185,] -0.0398192508
#> [186,]  0.1834279598
#> [187,] -0.5912969000
#> [188,] -0.2066219749
#> [189,] -0.0888385341
#> [190,]  0.7256313251
#> [191,] -0.8035350400
#> [192,] -0.4836962669
#> [193,]  0.0720480343
#> [194,] -0.1082587666
#> [195,] -0.4720423706
#> [196,] -0.5112296322
#> [197,] -1.0419494434
#> [198,] -1.0722192125
#> [199,] -0.2781924342
#> [200,] -0.7298848556
#> [201,] -0.2131483743
#> [202,]  0.5623254285
#> [203,]  0.9414450441
#> [204,] -0.6069610428
#> [205,] -0.4184906420
#> [206,] -0.5260371536
#> [207,]  0.2785360711
#> [208,] -1.1459601815
#> [209,] -0.3557825140
#> [210,]  0.3968760321
#> [211,] -0.6521965163
#> [212,]  0.4098542268
#> [213,]  0.1469716944
#> [214,]  0.0760444959
#> [215,] -0.5777144122
#> [216,]  0.4361240428
#> [217,]  0.3287092767
#> [218,] -0.7245330034
#> [219,]  0.7281972427
#> [220,] -0.2030511223
#> [221,] -0.2572487886
#> [222,] -0.0367713950
#> [223,] -0.5508172725
#> [224,] -0.0243384651
#> [225,] -0.2076859101
#> [226,]  0.4070828050
#> [227,] -0.2670716220
#> [228,] -0.2676964120
#> [229,] -0.1630404940
#> [230,]  0.0306441239
#> [231,] -0.2953236286
#> [232,] -0.9110260684
#> [233,] -0.6558410570
#> [234,] -0.6014122155
#> [235,] -0.5402654194
#> [236,]  0.7678161306
#> [237,] -0.3164934574
#> [238,] -0.6413990951
#> [239,] -0.6200251163
#> [240,]  0.1875005912
#> [241,] -1.0487163011
#> [242,]  0.8245525852
#> [243,] -0.7278083517
#> [244,] -0.2770352490
#> [245,] -0.4390723542
#> [246,] -0.1504646271
#> [247,] -0.4666944634
#> [248,] -0.3321764717
#> [249,] -1.1106327104
#> [250,] -0.8930300846
#> [251,] -0.4085696418
#> [252,]  0.3204002574
#> [253,] -0.4930624455
#> [254,] -0.2071144333
#> [255,]  0.4913781459
#> [256,] -0.5923887634
#> [257,] -0.3720274090
#> [258,]  0.0458658237
#> [259,]  0.4144482972
#> [260,] -0.1475580482
#> [261,] -0.5887903841
#> [262,] -0.6664111365
#> [263,] -1.3978900142
#> [264,] -0.0448620102
#> [265,] -0.8783962784
#> [266,] -0.0715598418
#> [267,] -0.4860027891
#> [268,] -0.3755217249
#> [269,] -0.0469819727
#> [270,] -0.1281685196
#> [271,]  0.0403596670
#> [272,]  0.4940918295
#> [273,]  0.1004728408
#> [274,] -0.1080040107
#> [275,] -0.0735316069
#> [276,]  0.1818042682
#> [277,] -0.5257711258
#> [278,] -0.6450204249
#> [279,] -1.2192756750
#> [280,]  0.9149875949
#> [281,] -0.4634201756
#> [282,] -1.1594803194
#> [283,] -0.5278889645
#> [284,] -0.8321031605
#> [285,]  0.3151377718
#> [286,] -0.2171719257
#> [287,] -0.4185623074
#> [288,] -0.2898486571
#> [289,] -0.3327101718
#> [290,] -0.9009374771
#> [291,] -0.3388145214
#> [292,]  0.0187280554
#> [293,] -0.0803742576
#> [294,] -0.0405584665
#> [295,] -0.2596904919
#> [296,]  0.1736535057
#> [297,] -0.1685337860
#> [298,] -0.7554978247
#> [299,] -0.6090169632
#> [300,]  0.5828235134
#> [301,] -0.7251970988
#> [302,]  0.5914926565
#> [303,] -0.9620588912
#> [304,] -0.6837451440
#> [305,] -1.7522071390
#> [306,] -0.0589709997
#> [307,] -0.2920687686
#> [308,] -1.3953142757
#> [309,] -1.6274345399
#> [310,] -0.6155316939
#> [311,]  0.1048655103
#> [312,] -0.5340161211
#> [313,] -0.0852866306
#> [314,]  0.4977310084
#> [315,] -0.6899107800
#> [316,] -0.1876487785
#> [317,]  0.2498355590
#> [318,] -0.0096422278
#> [319,] -0.1073883319
#> [320,] -0.5738523751
#> [321,]  0.7247516046
#> [322,] -0.7446658111
#> [323,] -0.3214074893
#> [324,] -0.1226393361
#> [325,]  0.5525206881
#> [326,] -0.4031692438
#> [327,] -0.4005963518
#> [328,] -0.7024027698
#> [329,] -0.6727870318
#> [330,] -0.4487776409
#> [331,] -0.1067751417
#> [332,]  0.3952392057
#> [333,] -0.5804170163
#> [334,]  0.1129850114
#> [335,] -1.3279212739
#> [336,] -0.6012439037
#> [337,] -1.1458119596
#> [338,] -0.0012072262
#> [339,] -0.0213213358
#> [340,]  0.0484477034
#> [341,] -0.0865161050
#> [342,]  0.1469372031
#> [343,] -0.0891093903
#> [344,]  0.4618617592
#> [345,] -0.7047893627
#> [346,] -0.4497626405
#> [347,]  0.1073283464
#> [348,] -0.1198912335
#> [349,] -0.3504789032
#> [350,]  0.1712833401
#> [351,]  0.6702069775
#> [352,] -0.4726748204
#> [353,] -0.4403241236
#> [354,] -0.7800096225
#> [355,]  0.7696533132
#> [356,] -0.1444085367
#> [357,] -0.2439147156
#> [358,] -0.1745224595
#> [359,] -0.9547773912
#> [360,] -0.3803543967
#> [361,]  0.7925376201
#> [362,]  0.3231902902
#> [363,] -0.4730283747
#> [364,] -1.1992397584
#> [365,] -0.7831723705
#> [366,] -0.9888201717
#> [367,] -1.0066185975
#> [368,] -0.2342077225
#> [369,] -0.9890180557
#> [370,] -0.5102124906
#> [371,] -0.0162884610
#> [372,] -2.4492852571
#> [373,] -0.5143092244
#> [374,] -0.0482651052
#> [375,] -0.2303028202
#> [376,] -0.1519235758
#> [377,] -0.8390402673
#> [378,] -0.0522494822
#> [379,] -0.3431881896
#> [380,] -0.2245782169
#> [381,] -0.1648470654
#> [382,]  0.0518432689
#> [383,]  0.0267792038
#> [384,]  0.0196863958
#> [385,]  0.1449207443
#> [386,]  0.4823379752
#> [387,] -0.9471782292
#> [388,]  0.2292955959
#> [389,] -0.2189994325
#> [390,]  0.1874208712
#> [391,] -0.9461459568
#> [392,] -0.8261495069
#> [393,] -0.3734653281
#> [394,] -0.2700175813
#> [395,] -0.3507965061
#> [396,] -1.3765645387
#> [397,]  0.1761583655
#> [398,] -0.5179337181
#> [399,] -0.5448112713
#> [400,] -1.2852919198
#> [401,]  0.5367975179
#> [402,] -0.4251427384
#> [403,]  0.7648416259
#> [404,] -0.9439808267
#> [405,] -1.2804902418
#> [406,] -0.2855544084
#> [407,] -1.4393002292
#> [408,] -1.7996031056
#> [409,] -0.1502754048
#> [410,]  0.4227751823
#> [411,]  0.1555528633
#> [412,] -0.4903789083
#> [413,] -0.3023399232
#> [414,] -0.5415447856
#> [415,]  0.8438426320
#> [416,] -0.3814149831
#> [417,]  0.0889398442
#> [418,]  0.4841212549
#> [419,] -0.3933696681
#> [420,] -1.2007015436
#> [421,]  0.1098637709
#> [422,]  0.5923493570
#> [423,] -1.0198043495
#> [424,] -0.2050593014
#> [425,] -0.3897260160
#> [426,] -0.2483499829
#> [427,] -0.6995334390
#> [428,] -0.2569381821
#> [429,] -0.5509036089
#> [430,] -0.6170393098
#> [431,] -0.0556232446
#> [432,] -0.0865751032
#> [433,]  0.1038814645
#> [434,] -0.0529145511
#> [435,]  0.7437587219
#> [436,] -0.0758669473
#> [437,] -1.3816449370
#> [438,] -0.4798458384
#> [439,]  0.7019740320
#> [440,] -0.1441781358
#> [441,] -0.2869018115
#> [442,] -0.4626520928
#> [443,] -0.2220824248
#> [444,] -0.1020591830
#> [445,] -0.1324211136
#> [446,] -0.7999913999
#> [447,] -0.7474174397
#> [448,]  0.2328805148
#> [449,] -0.0152224756
#> [450,] -1.1441947865
#> [451,] -1.3433341138
#> [452,] -0.1626975935
#> [453,] -0.8590753077
#> [454,] -1.0953062908
#> [455,] -0.2954360869
#> [456,]  0.6581346658
#> [457,] -0.4760083928
#> [458,] -0.5893953353
#> [459,] -0.6438825540
#> [460,]  0.9328126684
#> [461,]  0.3456163014
#> [462,] -0.6645764143
#> [463,]  0.1050438683
#> [464,] -0.7497901786
#> [465,]  0.3476046746
#> [466,]  0.5653861963
#> [467,] -1.1907482326
#> [468,] -0.7181124293
#> [469,] -0.6611656641
#> [470,]  0.5470056907
#> [471,]  0.3650651297
#> [472,] -0.9502582858
#> [473,]  0.0372055617
#> [474,] -0.6217030769
#> [475,]  0.3656589932
#> [476,] -0.7539907194
#> [477,] -0.2372431881
#> [478,]  0.4993480693
#> [479,] -0.6732882557
#> [480,] -0.3986210041
#> [481,] -0.8497880912
#> [482,] -0.1765288070
#> [483,]  0.0351993122
#> [484,] -0.9693702604
#> [485,] -0.4800853942
#> [486,] -0.4594851252
#> [487,]  1.2636557302
#> [488,] -0.4976283715
#> [489,] -0.4344557803
#> [490,]  0.0839517914
#> [491,] -0.6473230163
#> [492,] -0.2128173567
#> [493,] -0.9536845463
#> [494,] -0.9145626524
#> [495,]  0.8772639275
#> [496,] -1.2102771996
#> [497,] -0.4918702743
#> [498,] -0.0587696581
#> [499,] -0.6864822584
#> [500,] -0.2934595632
#> [501,]  0.0011372711
#> [502,] -0.8136387530
#> [503,] -0.2593691559
#> [504,]  0.1537571006
#> [505,]  0.3731537711
#> [506,] -0.0998489406
#> [507,]  0.1410799925
#> [508,] -0.1792960176
#> [509,] -0.5174754302
#> [510,] -0.9320177357
#> [511,] -0.3246611240
#> [512,] -0.4465846090
#> [513,] -0.3651306172
#> [514,] -0.3487912208
#> [515,] -0.1067710802
#> [516,] -0.2797867526
#> [517,]  0.3008880774
#> [518,] -0.3089512432
#> [519,] -0.4193687567
#> [520,] -0.1426524080
#> attr(,"df")
#> [1] 3 2
#> attr(,"range")
#> [1]  0 40
#> attr(,"lag")
#> [1]  0 85
#> attr(,"argvar")
#> attr(,"argvar")$fun
#> [1] "ns"
#> 
#> attr(,"argvar")$knots
#> [1] 0.06716823 0.53734584
#> 
#> attr(,"argvar")$intercept
#> [1] FALSE
#> 
#> attr(,"argvar")$Boundary.knots
#> [1]  0 40
#> 
#> attr(,"arglag")
#> attr(,"arglag")$fun
#> [1] "ns"
#> 
#> attr(,"arglag")$knots
#> numeric(0)
#> 
#> attr(,"arglag")$intercept
#> [1] TRUE
#> 
#> attr(,"arglag")$Boundary.knots
#> [1]  0 85
#> 
#> attr(,"class")
#> [1] "crossbasis" "matrix"    
#> 
#> $wetness
#>        wetness_v1.l1 wetness_v1.l2 wetness_v2.l1 wetness_v2.l2 wetness_v3.l1
#>   [1,]   11.12909284   3.981042098      21.06198    0.17721607   -9.28802903
#>   [2,]   10.63292490   3.564103406      20.82608    0.86063541   -9.02321825
#>   [3,]   -3.50670887   0.670818337      19.60646    2.75691687  -10.99731725
#>   [4,]   17.28727375   0.584494518      18.08473    1.37664591   -4.12116770
#>   [5,]    5.72161587  -0.630209992      22.06087    2.12289094  -11.75344038
#>   [6,]    1.99665425   2.927081165      21.22029    1.76173121  -11.31277318
#>   [7,]    9.88367946   3.830501574      21.62510    0.35317573  -10.25092443
#>   [8,]   13.80513678   2.253078228      19.85827    0.66921798   -7.91976968
#>   [9,]   -1.02031978   0.399499786      21.01392    3.56858456  -11.49145060
#>  [10,]    3.06365504  -0.363843814      23.55436    2.29839817  -13.04637246
#>  [11,]   17.23471394   1.958880251      18.20884    0.77453693   -4.14137857
#>  [12,]   -3.12248337  -0.567931162      16.55922    0.99171023   -9.30049675
#>  [13,]    6.35449823   2.372479772      22.98622    1.17349586  -12.35673835
#>  [14,]    9.35951542   0.888475437      21.31741    0.92459317  -10.18370223
#>  [15,]   14.46717363   2.156070653      18.35040    0.19842533   -2.29471859
#>  [16,]    8.32428200   3.814366115      20.61764    0.57060417   -9.92496527
#>  [17,]   -3.53601616   0.162410480      18.40730    4.31118140  -10.35016046
#>  [18,]   10.78393264   2.631727079      21.11585    1.11245134   -8.77774056
#>  [19,]    2.65870828   0.968252793      22.60555    1.79808084  -12.48158069
#>  [20,]   17.01546737   1.234105077      17.78046    1.52024616   -2.38864325
#>  [21,]   20.42548495   1.465727024      16.89038    1.53686268   -2.94486883
#>  [22,]   11.99030770   0.791356602      20.41589    1.63571450   -8.91100051
#>  [23,]    8.21975085   1.575233081      22.28280    1.27391529  -11.76193900
#>  [24,]   19.29313526   1.290078387      16.50765    0.98600716   -0.13020020
#>  [25,]   12.30750202   2.398898087      20.89369    1.06752521   -9.31634245
#>  [26,]   13.86170331   1.676867114      20.13214    1.67915810   -8.60507038
#>  [27,]   -4.71782129  -0.102106617      17.50556    3.61422585   -9.85326391
#>  [28,]   16.11975051   3.729677755      19.15959    0.64003925   -6.85985988
#>  [29,]    9.66738767   1.803312037      21.25722    0.82168252   -9.29196428
#>  [30,]   -1.01703043   2.045378624      22.24848    2.09344125  -12.43114201
#>  [31,]    1.31494607  -1.545792115      22.53046    2.55437242  -12.32016403
#>  [32,]    5.61114426   1.411283065      22.56266    1.46000230  -11.75045658
#>  [33,]   18.30818913   0.729955299      17.95235    1.35958812   -4.47401258
#>  [34,]    0.02045654   0.832052574      21.03268    2.43977452  -11.54390940
#>  [35,]    2.91937480   0.955894119      21.04400    2.35002879  -11.11203981
#>  [36,]    8.29989922   3.649148239      22.43334    1.17479692  -11.87579592
#>  [37,]   12.70937365   1.891531797      19.91429    2.01653977   -7.61226768
#>  [38,]   16.35450969   2.741570895      17.55621    0.66268469   -0.68231150
#>  [39,]   15.57011321  -0.734168149      16.12386    0.79135445    5.84640171
#>  [40,]    5.96550015   2.573669204      20.77457    1.19627833   -9.72971608
#>  [41,]   14.74534857   3.829500228      17.93152   -0.79722537   -0.63380869
#>  [42,]    4.99376200  -0.334194103      22.17324    1.07954600  -11.79046715
#>  [43,]   10.03977624   1.673114741      20.78921    1.56111689   -9.01731712
#>  [44,]   19.19042159   1.842671588      17.20457    0.37185071   -2.62545298
#>  [45,]    1.34433524   2.623073204      21.99378    1.86596178  -12.00579816
#>  [46,]    6.66409980   2.968199258      22.33789    0.45993999  -11.16273147
#>  [47,]   15.47480257   3.684434518      19.71645    0.38239767   -8.14833470
#>  [48,]   16.58812098   2.399457641      19.09491    0.62458986   -6.75286765
#>  [49,]   -2.38130325  -0.255425876      21.02606    2.28130808  -11.80249892
#>  [50,]    9.62868671   2.142564262      20.53843    0.99092625   -8.40885768
#>  [51,]   16.25921498   2.219487678      18.40171    0.73447535   -3.73271569
#>  [52,]   17.12273008   1.041777573      18.44285    1.81581082   -4.99155612
#>  [53,]    7.73203489   5.014454546      21.88036    0.05371877  -10.28267933
#>  [54,]   12.46954379   4.615229270      20.23717    0.74427505   -8.79680593
#>  [55,]    5.44265194   2.392258087      21.85462    2.44760540  -11.29189350
#>  [56,]   14.90061467   1.642231142      19.46926    1.09174440   -6.99801269
#>  [57,]   14.51689804   0.985361344      19.43211    1.73132776   -6.80021901
#>  [58,]    9.76058539   2.355116341      21.49549    1.79687714  -10.92011560
#>  [59,]   13.23980153   4.281145343      19.04958    0.47712235   -4.15656901
#>  [60,]    3.84378250   2.781969157      22.81661    1.59019551  -12.51823837
#>  [61,]   10.67077554   1.478852954      21.61074    1.60312086  -10.55011553
#>  [62,]   16.25577804   2.208591089      18.98083    0.76972423   -6.27798368
#>  [63,]   -3.37535531   0.046108868      19.57112    4.57782182  -11.00533825
#>  [64,]   13.80576706   1.330776537      20.54434    1.70058805   -9.19553050
#>  [65,]    7.86914029  -0.921833758      22.41490    2.14099853  -11.63549438
#>  [66,]    5.98504527   4.483237401      22.67505    0.92114706  -12.11814713
#>  [67,]    9.16661919  -0.320945223      21.45217    1.75521279  -10.46159914
#>  [68,]   11.25549664   3.432410307      19.71704    0.74033613   -4.88659266
#>  [69,]   10.58465839  -0.004959931      21.96253    1.88662913  -11.37239614
#>  [70,]    5.52757242   1.580224677      22.20651    1.14397945  -11.90498012
#>  [71,]    4.95088111   1.405255215      22.19782    1.06914955  -11.42736913
#>  [72,]   14.57196493   4.059424075      20.21307    0.32280499   -8.86379145
#>  [73,]    0.28446189   1.448955619      23.00504    1.54685064  -12.83429661
#>  [74,]    6.21220866   2.392351982      22.26468    1.64454690  -11.33594300
#>  [75,]   14.43281366   2.332342167      18.22583    0.42378156   -1.22253738
#>  [76,]   -0.23622893   1.060716184      20.82858    1.66421021  -11.50428435
#>  [77,]    5.44655393   0.924447973      21.71753    1.33209042  -11.64754229
#>  [78,]    4.79811680   3.803430767      22.91436    1.33238443  -12.38922184
#>  [79,]    9.43919369   1.451240773      21.70727    1.56114722   -9.53572685
#>  [80,]   15.01916482   2.620742902      19.04010    0.81792376   -5.83123737
#>  [81,]    8.61295011   2.152428551      21.53742    2.47441714  -11.32060646
#>  [82,]    8.20262505  -0.249018140      21.67921    2.64035947  -10.91920874
#>  [83,]   17.86800989  -0.595098156      15.98070    0.88386905    3.58122668
#>  [84,]    8.13677819   4.759754843      21.97097    0.40533165  -11.03024656
#>  [85,]    4.19838878   3.561704257      23.39049    1.32124165  -12.73979887
#>  [86,]    1.82976728   1.086800981      23.64921    2.11656913  -13.00785666
#>  [87,]   17.80441472   2.299193541      17.34009    1.15846890   -1.52716142
#>  [88,]   11.17374918   1.675526168      21.04785    1.39516231   -9.31102284
#>  [89,]    5.56818797   1.400093259      21.92026    1.60696610  -11.18079177
#>  [90,]   18.45979450   1.876107904      17.45082    1.57594709   -2.66243072
#>  [91,]   14.70424558   2.101796041      18.96854    1.06679048   -4.18731354
#>  [92,]    7.61529614   1.825744924      19.28717    2.98705061   -9.44187324
#>  [93,]   10.23153496   0.463503757      21.63553    1.75853068  -10.70934800
#>  [94,]   -2.86136419  -0.994397960      22.00067    2.07831228  -12.33813000
#>  [95,]   14.11030372   0.759045334      19.17838    0.66631305   -5.07620980
#>  [96,]   18.51751291   2.491512466      17.41529    0.94052982   -2.73379978
#>  [97,]   17.96560979   0.556054967      18.37685    2.06293188   -5.72606233
#>  [98,]    6.97890823   4.105022635      20.62613   -0.57660127   -8.85497996
#>  [99,]    4.51381161   2.251891889      22.53081    1.90050725  -11.94074247
#> [100,]   11.51940897   1.871487198      21.19121    1.23543346   -9.96783851
#> [101,]    6.27775724   3.925892757      22.44647    1.29368270  -11.67065785
#> [102,]   -1.24795730  -1.833160963      22.23173    2.90772970  -12.31045946
#> [103,]   -1.35371378   0.205463290      19.75639    1.28625815  -10.90054546
#> [104,]   13.77529372   1.119208861      20.76479    1.64519379   -9.91037393
#> [105,]   10.72767054   2.760718473      21.43775    0.87633348   -9.96885335
#> [106,]   16.02951989   4.830677845      18.11340    0.34924645   -3.65388993
#> [107,]   10.00595856   2.188197966      21.45262    1.19084475   -9.20610241
#> [108,]   -0.99976394  -0.731488658      22.17193    1.12607400  -12.34921762
#> [109,]   17.97033583   2.093491967      17.63232    0.62720696   -2.85705621
#> [110,]    1.99741859  -1.184126891      22.31284    2.66710432  -12.05191751
#> [111,]    0.86039644  -0.199959361      23.55866    1.57203483  -13.12998360
#> [112,]   -0.97141809   2.556199899      20.29694    2.51983273  -11.21965688
#> [113,]   10.91493845   2.161842635      20.68683    1.44177249   -7.60760449
#> [114,]    4.13576430   0.250209649      23.36313    2.09239424  -12.73929105
#> [115,]   11.35675725   3.941431478      20.75835    0.79262960   -9.06193131
#> [116,]   -1.26615753   3.312996619      21.53126    1.92621686  -11.73210009
#> [117,]    0.20951527  -1.948669803      23.56107    2.16919637  -13.04968741
#> [118,]   19.50033128   2.741793403      16.61537    0.91078495   -0.80613955
#> [119,]    8.37307769   2.555025757      21.48339    1.42258533  -10.44573328
#> [120,]   14.58457210   1.210832493      19.80100    1.29051666   -8.01991057
#> [121,]   15.79859096   0.761035657      18.74024    1.06185844   -4.41314590
#> [122,]   -1.72003888   1.073245848      20.24536    3.00903415  -11.34382801
#> [123,]   -1.68084194   2.889243763      20.30957    2.13725773  -11.34248955
#> [124,]   13.49994066   1.709104898      20.42958    1.52538754   -8.79806277
#> [125,]   11.33588652   3.141632334      20.42222    0.23859486   -8.83517158
#> [126,]    6.21560372   5.874987344      21.77166    0.70656182  -10.54302585
#> [127,]   15.53936037   2.075005940      18.48988    0.79062762   -4.59131471
#> [128,]   16.46636417   1.586233579      18.24464    1.72753903   -3.63386966
#> [129,]   12.54353076   2.927659560      21.11497    1.00262873  -10.22606530
#> [130,]    6.66285453   4.660482167      21.33880    0.68867726  -10.21951935
#> [131,]   10.13504905   2.531190164      21.04549    1.54829910   -9.87834913
#> [132,]   -3.76410164  -1.118833277      15.04491    4.37533822   -8.46505786
#> [133,]   -3.23944817  -0.018912952      23.47933    2.09167768  -13.20546136
#> [134,]    8.31752249   3.356215248      22.18291    0.63667635  -11.50525002
#> [135,]   -3.01731693  -1.133188921      18.21628    2.08451254  -10.20390188
#> [136,]    9.27629161  -0.016009433      22.20766    1.95954738  -11.45682973
#> [137,]    6.79676501   0.582745528      22.49519    1.84914739  -11.26576367
#> [138,]   11.09123369   1.297153872      21.39500    1.59140258  -10.18769092
#> [139,]   -1.70572702   1.591283194      22.39626    2.30771127  -12.51723900
#> [140,]   12.27183037   4.082186701      20.69150    0.32876947   -8.81317261
#> [141,]   -4.79580939  -0.277829632      15.66389    3.15613042   -8.81675267
#> [142,]   17.41130967   1.274645271      18.54104    0.85316489   -5.66392468
#> [143,]   -4.25320403  -0.177623197      18.84660    3.01024770  -10.60277330
#> [144,]    6.87536023   1.927261924      22.13403    0.53349692  -11.13390510
#> [145,]    7.32502057   1.236566000      22.94716    1.64119862  -12.20378855
#> [146,]   15.58930008   2.812371982      17.97472    0.19987756   -1.43771016
#> [147,]   15.05954255  -0.340229026      20.08864    2.15673554   -8.88244067
#> [148,]   11.75680810  -1.194395059      20.48017    2.49650276   -8.12088827
#> [149,]    2.04246906   0.800662299      19.23429    2.32809220  -10.52053009
#> [150,]   15.58596920   2.698980368      17.31309    0.46569138    1.24875120
#> [151,]    8.16351996   3.946598411      22.22892    1.17108701  -11.71150379
#> [152,]   -3.62051411   0.536001701      16.83044    1.69412538   -9.44874282
#> [153,]   17.32560838   2.229942542      17.98473    0.88972421   -3.42297038
#> [154,]   16.43002044   2.662935713      18.24044    0.43618978   -3.94183212
#> [155,]    0.94862472  -0.661493624      22.17879    2.92374277  -11.93881907
#> [156,]   17.08295288   1.926155608      18.41913    0.95841385   -4.92436436
#> [157,]   -0.42712297   3.230913158      21.93990    2.42200626  -12.15979895
#> [158,]    8.51366797   1.155796998      22.03970    1.95753119  -11.50626542
#> [159,]   16.90712691   2.604732267      18.01996    0.77070136   -3.13640208
#> [160,]   13.62271901  -0.373566217      20.39943    2.32094606   -8.75304048
#> [161,]    1.62089021  -1.399911557      22.71174    1.75082750  -12.57994207
#> [162,]   10.10422371   1.905081092      21.25791    0.83629977  -10.23589584
#> [163,]   15.30189814  -0.941007917      16.50286    1.06174701    4.72710164
#> [164,]    5.08169493   2.815571999      22.95969    1.53931422  -12.37985241
#> [165,]    3.78178039   1.261042277      22.72069    1.31232946  -12.40515943
#> [166,]   17.97861786   2.008622682      17.33055    1.06767019   -1.67359118
#> [167,]    6.27373004   3.308939932      21.21558    1.35114613  -10.93485487
#> [168,]   16.26076695  -0.530088483      19.05809    2.73279240   -6.41430970
#> [169,]    3.51338585   3.776926098      22.11406    1.11122831  -11.90036198
#> [170,]   10.36315614   2.937881748      20.74958    0.99945280   -7.69807664
#> [171,]   11.32875072  -0.261147697      21.31212    2.24400395  -10.00439481
#> [172,]   13.50399681   6.241955474      19.41751    0.41764066   -7.11621467
#> [173,]    4.84671216   2.115690148      23.65737    1.54912446  -13.01488299
#> [174,]    5.61735678   3.524900147      21.71050    1.62413923  -11.29886004
#> [175,]   11.63503335   1.697375098      20.87455    1.62114435   -9.54170606
#> [176,]    6.43664696   2.519053748      22.07493    1.16836882  -11.39482876
#> [177,]   -0.93005184   1.183298378      23.23624    2.43474491  -12.99487827
#> [178,]   -0.20086928   2.553248666      22.37709    1.69465605  -12.41153365
#> [179,]    2.16214692   0.757957222      23.23479    1.91989150  -12.81434080
#> [180,]    9.31146665   3.649276524      22.10905    0.99304980  -11.44839215
#> [181,]   -3.52086837  -0.579442379      19.91593    1.66305193  -11.20232274
#> [182,]    3.85692346   2.146392141      22.89718    1.13642033  -12.54622098
#> [183,]   17.22878812   2.352960576      16.69597    0.22110282    1.48073215
#> [184,]   12.31770420   2.128548561      20.45554    1.65786139   -8.91470885
#> [185,]   12.37827451   3.317418669      21.02998    0.94975632  -10.08494693
#> [186,]    8.53133056   3.520566210      21.35096    0.95499949   -9.63916288
#> [187,]    9.03387533   2.884404931      21.40281    0.85321796  -10.08056103
#> [188,]   14.96498767   4.420146789      18.70486   -0.42573586   -3.91803883
#> [189,]    3.68014939   1.452266997      22.77776    1.22588681  -12.21923046
#> [190,]    0.70207000  -0.209228343      22.59056    1.86850558  -12.46479618
#> [191,]    7.58242325   2.173994358      22.33603    1.61500286  -11.23586886
#> [192,]    3.46520719   4.135362606      23.51349    0.78261523  -12.64627844
#> [193,]    6.62120774   0.828447699      23.06629    1.76649009  -12.49255644
#> [194,]   -0.26205945   0.545418543      22.06966    1.85945076  -12.26685421
#> [195,]   17.47068847   1.440240016      18.45124    0.87968429   -6.05741857
#> [196,]   -4.89676638  -1.001965376      21.09957    1.68445956  -11.86518656
#> [197,]    8.36923640   3.341058581      22.04515    0.76615755  -11.05862606
#> [198,]    6.65566343   1.030516503      21.98666    1.70970573  -11.41528354
#> [199,]    5.26662520  -0.698232059      22.97820    2.33184683  -12.12818196
#> [200,]   19.24856281   2.517660211      17.26286    1.12703817   -2.93466530
#> [201,]    2.55414965   1.133305252      22.98850    2.51096806  -12.07198878
#> [202,]    8.04245516   3.623167219      21.61592    1.56203940  -11.11292047
#> [203,]   -3.43155157  -0.683886401      11.57318    4.03949892   -6.51126694
#> [204,]    7.43614966   4.511308382      22.66650    0.48886941  -11.85953781
#> [205,]    7.15467592   3.591574013      22.36591    1.30927650  -11.64225124
#> [206,]   17.71663933   2.964611591      17.98591    1.03570250   -3.83728791
#> [207,]   -2.52044401  -0.801190208      21.04493    2.81564018  -11.80804288
#> [208,]   13.21428406   4.508980591      20.40619   -0.02622914   -8.59443726
#> [209,]   -4.05712075  -1.140069091      14.46062    3.39066678   -8.13523536
#> [210,]   12.43965661   0.397100924      19.99175    1.70000972   -7.62020934
#> [211,]    9.31206716   4.925108611      21.14520    0.34571131   -9.52326110
#> [212,]    5.38613391  -0.642440630      22.94837    1.94186189  -12.30815006
#> [213,]    1.38492893  -1.500863931      21.63053    1.84397269  -11.84094104
#> [214,]   13.68168355   2.107622040      19.89610    1.13746898   -6.92177509
#> [215,]   11.78461070   0.524183676      20.32208    1.72555440   -8.54929492
#> [216,]    9.21559803   2.775277708      22.18191    1.11945595  -11.49050506
#> [217,]   -4.08468122   0.895933798      19.14692    2.66760701  -10.77502586
#> [218,]    9.19184840   1.437901250      21.37452    2.11128778   -9.53356084
#> [219,]   15.16917820   3.489646415      19.86558    0.69343477   -8.38883016
#> [220,]   17.97125035   2.022308740      17.07579    1.05208504   -0.66414730
#> [221,]    6.61482431   1.544968906      23.05733    1.54445578  -12.37372067
#> [222,]   17.93846359   1.924258545      17.52568    1.13509620   -2.32339176
#> [223,]    2.36988264   0.739966714      23.27170    1.57635589  -12.60657733
#> [224,]   10.92596973   1.981640431      21.07419    1.90572726   -9.56343153
#> [225,]   19.11284329   0.797759274      16.81506    1.37738767   -1.07569067
#> [226,]   -3.85392551  -0.392922257      15.03555    2.78227381   -8.46301368
#> [227,]    7.24813000   1.708587994      23.05717    1.38212535  -12.27648624
#> [228,]   -0.13669094   2.076475071      21.35439    2.22207575  -11.84725558
#> [229,]   11.97172041   0.708514408      20.91637    2.02421087   -9.19261373
#> [230,]   14.75897176   1.479848573      19.15218    0.14657298   -5.18263628
#> [231,]   16.80356700   1.981053898      17.64462    0.35693030   -1.42505260
#> [232,]    7.93318502   2.350176668      21.61668    1.85537032  -10.88736633
#> [233,]    4.95947166  -1.138011116      22.44047    2.16615553  -11.87843369
#> [234,]    5.19823341  -0.405447113      22.69987    1.98259951  -11.87714020
#> [235,]    9.26437508   1.480995550      20.71033    1.42996267   -7.79144136
#> [236,]   -1.31878346   3.173710190      21.46823    3.05704043  -11.92631437
#> [237,]   11.81242662   4.395115134      20.57411    1.08652294   -9.38207265
#> [238,]   -0.41504355  -0.203140312      22.49665    2.12580769  -12.52724716
#> [239,]   -0.52674189   0.327187269      19.25746    0.96308669  -10.73403064
#> [240,]   14.91747300   1.442750806      19.95951    1.12093642   -8.43754043
#> [241,]   -4.54443679  -0.455474367      15.54097    2.73522141   -8.74705967
#> [242,]    8.36590110   4.589674962      20.78259    0.27033168   -9.15533739
#> [243,]    6.88976138   0.606694434      22.06237    2.58897852  -10.80404713
#> [244,]    4.78901451   1.282784483      22.26180    1.09740814  -11.84136271
#> [245,]   -2.31174636   1.024259819      22.66732    2.37254649  -12.73564860
#> [246,]   -4.57836582   0.558555354      17.17431    1.88676896   -9.63862977
#> [247,]    5.27696572   1.474887740      23.18650    1.69556066  -12.64799770
#> [248,]   14.57663183   2.806633933      19.80721    0.83357603   -7.78896422
#> [249,]   19.04453365   1.165120365      17.65275    1.66888515   -4.17335248
#> [250,]    5.13619205   2.794673731      23.14886    1.49129129  -12.64935027
#> [251,]   12.30015155   1.992091554      20.88536    1.43759629   -9.44792910
#> [252,]   -2.95741524   1.772040134      21.28768    1.50166580  -11.94758763
#> [253,]    1.78206873   2.552926101      22.10199    2.46534758  -12.09199666
#> [254,]    0.49309934   0.897620796      22.16827    2.10066078  -12.27846064
#> [255,]   11.77645814  -0.139712071      21.22641    1.92720560  -10.32733176
#> [256,]    0.75822017   1.194828477      19.65111    1.91150676  -10.87053582
#> [257,]   -3.10430363   0.621850807      19.34229    4.95796475  -10.86049584
#> [258,]    1.57751944   0.716615516      22.21958    0.66930636  -12.07062550
#> [259,]    7.55165400   2.635207270      22.23967    1.17251697  -11.67317594
#> [260,]   17.11064052   2.837176699      18.19052    0.97215462   -3.93435201
#> [261,]    2.95201989   1.933849042      21.04057    2.58286898  -10.54375732
#> [262,]    1.22260635   0.701696482      22.90974    2.05690092  -12.65086000
#> [263,]    8.36180564   2.899252361      22.12590    0.74919680  -11.39595901
#> [264,]   -1.36398572   1.159077079      23.27889    2.23386831  -13.03697542
#> [265,]   16.24726890  -0.906549687      16.09687    1.34053846    5.12514297
#> [266,]    7.41651941   2.571112071      22.47097    1.12706447  -11.02502815
#> [267,]   -1.06713834   2.143295088      22.24236    2.29237086  -12.21614892
#> [268,]   11.30335015   1.050082456      21.64651    1.86435804  -10.63728797
#> [269,]    0.88628961  -0.780310369      23.27434    1.65786934  -12.84204156
#> [270,]   17.58561541   1.480952819      17.37655    1.22956465   -1.54285995
#> [271,]    1.89386989   0.899542473      23.28141    1.59906726  -12.85562464
#> [272,]   11.96598646   1.974447165      21.40487    1.24487319  -10.48398487
#> [273,]   10.57648289   2.812830933      21.56972    1.05928744  -10.56056889
#> [274,]   11.72069459   0.915913072      21.26528    1.76013633  -10.41429272
#> [275,]    2.21367983   1.451127345      22.84434    1.26450444  -12.63047958
#> [276,]   -2.46769640   0.923754295      21.02153    2.54870344  -11.81258949
#> [277,]   12.37105991   1.366031747      20.83747    1.60913488   -9.21057843
#> [278,]   -4.07278656  -0.462522205      14.39402    2.73023103   -8.09927854
#> [279,]   10.12988548   0.828193286      21.42418    2.15810744   -9.63865869
#> [280,]    8.48523062   2.471282722      22.25265    1.37592750  -11.42304327
#> [281,]   12.58748556   4.416657767      19.81495   -0.00926777   -6.52680200
#> [282,]   14.40629273   2.850159746      19.10044    0.75558657   -5.15573704
#> [283,]   16.20419400   0.901219071      19.24967    1.52957866   -7.21915642
#> [284,]   -2.61964061  -1.322886637      21.49578    1.14222829  -12.04825567
#> [285,]    8.52238196   3.048146658      20.70557    0.42488605   -9.51878029
#> [286,]   18.09386595   1.091768124      16.53923    1.27077593    1.08509263
#> [287,]   12.58386277   3.425561397      19.56479    0.40397765   -5.50679499
#> [288,]   16.12070006  -1.360394898      15.78293    0.57462485    6.41865691
#> [289,]   15.70407116  -0.122523357      18.22672    1.23083386   -3.12437274
#> [290,]   13.61464143   4.502616771      20.35283    0.45583897   -9.22101922
#> [291,]   -2.34972076   2.202359366      19.62073    3.05140667  -10.97945981
#> [292,]   -2.17605175   0.683551390      23.52329    2.17518882  -13.20018904
#> [293,]   15.69872479   1.093138881      18.90379    1.42833754   -5.57598136
#> [294,]   16.25220127   3.482844302      19.44434    0.33833217   -7.68029243
#> [295,]    4.52712796   2.969683477      19.65579    0.83360221   -9.97241837
#> [296,]   16.35250862   2.060830428      18.38988    1.36397589   -3.92159599
#> [297,]   19.20020697   2.164518629      16.42176    0.71098765    0.28219252
#> [298,]    9.24629471   2.414923599      21.14896    1.42773025  -10.27491065
#> [299,]   -4.39356631  -0.069001369      18.38413    0.62435576  -10.34733893
#> [300,]    9.60338133   4.344346717      21.47906    0.28112547  -10.32485838
#> [301,]   15.32696330   2.104246233      18.42931    0.36586630   -2.91852976
#> [302,]   -0.88715991   1.775323231      21.58344    1.65461528  -11.98254272
#> [303,]   14.70145021   1.593087329      19.19283    1.30114941   -5.73015162
#> [304,]   -1.03630527   0.308883592      22.28345    2.24028353  -12.33451225
#> [305,]    1.38276759   2.124125521      21.62833    0.60118915  -11.87548823
#> [306,]   -0.76516362   1.981405984      20.87669    0.52565565  -11.59032837
#> [307,]    4.41894411   3.999979302      22.67543    0.75993715  -11.93692950
#> [308,]    7.74316138   2.794348578      21.28131    0.72369529  -10.57946450
#> [309,]   14.34022414   2.286014875      17.98720    0.88929405   -0.10359924
#> [310,]    2.83558655   4.434833468      21.17950    1.91650885  -11.44785704
#> [311,]   -0.17828257   1.831894896      21.30004    1.89388017  -11.87987624
#> [312,]   14.69386436   2.714914977      20.00548    0.47471953   -7.90172834
#> [313,]    6.40668920   3.025912782      22.52399    1.14536659  -11.83431255
#> [314,]    1.42696311   2.540605574      18.95544    2.11090645  -10.35885344
#> [315,]   13.98428793   1.358101272      17.70796    0.91132402    1.69606554
#> [316,]    5.13379116   3.020232669      22.35006    1.30305583  -11.98913141
#> [317,]    0.39327924   0.656521623      20.90877    2.69038161  -11.48818878
#> [318,]   -2.94252078  -0.190425877      15.62112    5.24227158   -8.77030711
#> [319,]   12.16190308   4.171636794      20.94356    0.42616693   -9.60697831
#> [320,]   -3.46492473   0.404154884      21.22439    2.68124010  -11.93596759
#> [321,]   -0.21015026   0.819250891      19.33417    3.31016098  -10.73680703
#> [322,]   16.05569617   0.516341562      18.83662    1.39519278   -5.50809533
#> [323,]    9.85400682   4.361656174      20.09776   -0.24322979   -7.85251706
#> [324,]   -4.46872480  -0.561202660      17.35426    2.33148681   -9.76781230
#> [325,]   10.79643361  -0.067912694      20.39694    2.50336892   -5.95437550
#> [326,]    6.96691348   1.533648436      22.85480    1.49057085  -12.09674572
#> [327,]   10.09581916   2.002936256      21.54886    2.16397746  -10.59368374
#> [328,]    5.51466473   4.035480011      22.54877    1.22906038  -11.86815651
#> [329,]   12.10395266   1.843087696      20.86833    1.21271483   -9.39242465
#> [330,]   14.67761848   3.725540670      19.25627    0.34412738   -5.68243055
#> [331,]   15.39780466   3.271591804      19.80939    0.77824365   -8.57819133
#> [332,]   12.83372291   5.347200683      19.41600    1.18623727   -6.91461682
#> [333,]   12.12457154   3.755654957      20.19541    0.84895130   -7.28449045
#> [334,]    1.55485007   2.868263313      21.43941    1.48155122  -11.67366102
#> [335,]    2.39860572   2.502748363      22.35460    0.65002580  -11.95198901
#> [336,]    2.26373846   2.618457967      20.53363    2.04294428  -11.09350626
#> [337,]    8.36457225   2.659241777      22.25543    1.46102042  -11.20136195
#> [338,]    8.18314527   1.577319058      22.30588    1.49313076  -11.58632727
#> [339,]   11.67009920   3.027846995      20.65387    0.85858290   -8.74642342
#> [340,]    5.01635141   2.905114814      23.18724    1.21670352  -12.48519012
#> [341,]    0.20375885   2.299918921      22.57551    2.40023666  -12.53796759
#> [342,]   -1.01569281   0.401494965      18.49416    2.13515227  -10.36213064
#> [343,]   17.63274672   1.312909484      18.06884    0.87904101   -4.11262177
#> [344,]    6.83207847   1.431787011      22.64525    1.74367670  -11.89437107
#> [345,]   -3.10541002  -2.081396499      22.10178    0.36192152  -12.39993460
#> [346,]    1.64969400   2.280606225      22.52016    2.92875199  -12.42380922
#> [347,]    0.57955978   3.575624427      22.71711    2.31670880  -12.62033535
#> [348,]   10.12292739   3.600550505      21.47230    0.61182923  -10.04765559
#> [349,]    0.62002449   0.512748617      22.04859    1.20437313  -12.23427743
#> [350,]   16.84436100   0.301440312      18.80123    2.25425943   -5.86111795
#> [351,]    8.93922778   1.986222799      20.95349    1.96217398   -9.97068986
#> [352,]   16.09768901   0.345840341      18.11236    0.69092490   -2.55701027
#> [353,]   -3.97866059  -0.419332033      12.23597    1.25949273   -6.88689770
#> [354,]   11.69218645   1.807504767      20.72780    1.13231456   -9.24776565
#> [355,]   11.64478365   2.998924323      20.48047    1.16235343   -7.12541616
#> [356,]   15.75884681   3.095827895      18.27972    1.12764184   -3.55946399
#> [357,]    1.53567262   1.965921284      21.47134    2.43317928  -11.86967509
#> [358,]    3.81096142   2.343361851      23.12007    1.45514607  -12.53696098
#> [359,]    9.95368740  -1.950244153      21.30355    2.60170052  -10.34427277
#> [360,]    8.96160462   2.843181950      21.86185    1.10983168  -10.77939588
#> [361,]    9.53761657   0.647285181      21.74020    1.93344378  -10.09480937
#> [362,]   14.59372627   2.577829237      19.80721    0.85906448   -7.84716941
#> [363,]    6.39947651  -0.601033310      21.60523    1.95298113  -11.13277054
#> [364,]    4.05191814   4.267827531      21.90911    1.81304777  -11.78318795
#> [365,]   18.12092901   1.245828138      17.04739    1.51931463   -0.86428486
#> [366,]    2.74648087   2.294403973      22.11191    1.20160121  -11.68904389
#> [367,]   15.11138016   3.651541936      19.00113    0.75348416   -5.11801736
#> [368,]   12.01113797   1.390527587      20.78570    1.59117465   -9.02782292
#> [369,]    2.81163910   1.066635659      22.74659    2.05531920  -12.46256158
#> [370,]   14.88865891   1.861228345      19.62521    1.18893201   -7.30866665
#> [371,]   14.05342444   2.587615855      19.33211    0.53440974   -5.45987657
#> [372,]   -1.27165621  -1.897977474      18.29198    2.59199950  -10.18899167
#> [373,]    3.73859644   0.445853216      22.87785    2.07250237  -12.46061713
#> [374,]   -3.28397261   1.111292953      19.82488    2.66694460  -11.12851230
#> [375,]    7.75592681   2.786456871      22.53805    1.09079983  -11.60569231
#> [376,]    1.41743379   2.695787172      22.01289    1.34172265  -12.18496760
#> [377,]    4.48257670   1.848493885      22.24378    1.96045962  -11.85148786
#> [378,]    2.03510071  -2.047903339      21.97501    2.84695189  -11.14880149
#> [379,]    8.01587061   2.776634296      21.80092    1.18900324  -10.83966043
#> [380,]    8.19188329   2.433286603      21.72557    1.12815017  -10.69814539
#> [381,]   16.30459006   2.660734615      18.99006    0.95755938   -6.33402609
#> [382,]   12.74217109   4.417897355      19.94834    0.78168127   -6.79407455
#> [383,]   13.42026332   2.585096065      20.44234    0.81998989   -8.86529957
#> [384,]   12.21706812   0.376093329      20.83903    1.54163795   -9.51695216
#> [385,]    6.22041493   1.167075999      22.61130    1.42622755  -11.79390860
#> [386,]   13.48111457   3.496604034      20.59778    0.91259191   -9.44655509
#> [387,]    9.88978038   3.376199054      22.35535    0.89345339  -11.84300271
#> [388,]   10.21490002   2.956743199      21.23279    0.79412152   -8.43376499
#> [389,]    3.65945920   1.945301157      23.72345    1.46844498  -13.09720960
#> [390,]    5.26897526   1.712112341      22.75076    1.77484529  -12.13180889
#> [391,]    9.88067546   0.870032729      21.29299    1.48500095   -9.82033049
#> [392,]    4.54331890   3.777586038      22.64375    1.64068434  -12.08899228
#> [393,]    1.05074292  -0.329051201      21.72690    2.44393430  -11.87857036
#> [394,]    3.66760339   1.809841744      22.65346    1.94619241  -12.32698559
#> [395,]    5.14839514   0.975464541      22.97372    1.70150620  -12.27885486
#> [396,]   -2.36252692  -0.151978183      20.22080    3.59985001  -11.34922632
#> [397,]   16.53598465   2.054867582      15.83935    1.14576111    5.73696323
#> [398,]    5.60784412   0.480698072      22.44250    1.82945369  -11.09898657
#> [399,]   12.92229451   2.044725822      20.66061    1.27010470   -9.21181060
#> [400,]   11.64208413   4.153335629      20.29610    0.94298001   -9.03496712
#> [401,]   13.18977499   4.552522061      20.64584    0.14013738   -9.22686375
#> [402,]   18.15476436   1.140364859      15.84057    1.12960682    3.75630302
#> [403,]   -3.59165700  -0.492751210      10.79807    2.79351395   -6.07792225
#> [404,]   16.98421896   1.639916751      17.87073    1.07454722   -2.77962352
#> [405,]   11.40331192   3.740182146      20.51934   -0.01448401   -7.40019697
#> [406,]    3.45632322   6.272039830      21.39395    0.75619441  -10.47886913
#> [407,]   -4.08437533  -1.635365326      14.95875    4.90597471   -8.41862482
#> [408,]    4.88726469   5.031795940      22.05070    1.13947235  -11.59423949
#> [409,]   19.47087508   1.561921300      16.39142    0.68236595    0.08634134
#> [410,]   -2.58812884  -0.866960264      21.92219    1.78318443  -12.29783935
#> [411,]   10.58208147   3.779414702      21.78572    0.55754989  -10.92552243
#> [412,]    0.93995982   2.389994559      22.08341    2.13396342  -12.10290076
#> [413,]   12.53646465   0.841829390      20.88194    1.59323984   -9.83565366
#> [414,]    2.42697955   1.070654761      21.90323    1.40202926  -12.01193934
#> [415,]    6.29893640  -0.099913268      23.07307    2.14464612  -12.42118002
#> [416,]   13.79190769   2.240028771      20.24509    0.99222903   -8.43679466
#> [417,]   18.76925986   1.587256185      16.29521    1.03877011    1.28326873
#> [418,]    3.78425448   0.393486730      22.02076    2.70586790  -11.93203206
#> [419,]   11.75968652   2.987050501      20.59184    1.51821660   -9.55687014
#> [420,]   17.93722234   2.906748167      17.95823    0.73490000   -4.19187034
#> [421,]    9.50055668   2.929113002      21.13886    1.62394649   -9.62117589
#> [422,]    1.67120636  -0.256264005      21.89019    1.60025233  -11.30428461
#> [423,]    3.76708593   4.520596536      22.89376    1.81192359  -12.61724562
#> [424,]   -3.05825407  -0.169529334      23.02569    1.92169486  -12.95392967
#> [425,]    5.96821032   0.877484629      22.35066    2.38031301  -11.68147744
#> [426,]    7.24843032   2.973824582      21.81075    2.35769334  -11.36965366
#> [427,]   16.27783737   2.791089556      18.69733    0.54724023   -5.33834835
#> [428,]    3.04898401   3.192032395      22.11789    1.58909544  -12.02745951
#> [429,]   10.07675422   6.123684130      20.75821    0.54395951   -9.70962988
#> [430,]   10.65662692   2.439225283      21.17207    1.01880868   -9.60777036
#> [431,]   11.60836312   1.853062743      19.84361    1.00947230   -7.32315684
#> [432,]   -0.37605470   1.756606028      21.73136    0.88795532  -12.13195637
#> [433,]    5.22416583   0.600865504      21.28279    2.96599850  -10.37414121
#> [434,]   16.87592110   0.851588749      19.12883    1.63983513   -7.16021209
#> [435,]   -1.70154586  -0.277088307      23.55456    2.08529280  -13.20376075
#> [436,]   18.43166658   2.206021487      17.55826    1.04309921   -3.13990194
#> [437,]    0.96693455   2.238936123      22.50898    2.84317428  -12.51433351
#> [438,]   18.82942094   1.058882707      17.95095    1.68711563   -5.12527971
#> [439,]   12.18455739   2.659099878      20.39648    0.96763927   -7.60660144
#> [440,]   10.60905789   0.085172503      21.38702    1.96580417  -10.35135504
#> [441,]    1.88646988  -2.320536638      22.56300    1.17790651  -12.39125987
#> [442,]   17.07617817   1.231758466      18.45162    1.22181944   -4.94965973
#> [443,]   -2.44308931   0.900320677      22.28261    3.24571248  -12.48957220
#> [444,]   17.45353719   0.429995199      18.11634    1.30157044   -4.10425626
#> [445,]   10.92454437   2.025265794      20.60320    0.75370122   -9.82292987
#> [446,]    2.38657816   1.553854001      21.53530    2.50751538  -11.65319069
#> [447,]    1.59010838   0.580418512      23.05900    2.84072979  -12.76868035
#> [448,]   14.58260557   0.257330723      19.54320    2.15550583   -6.54014279
#> [449,]   -1.02913915   1.130888014      22.40728    3.17991761  -12.50037866
#> [450,]    8.75438753   2.512421004      21.19704    1.72515173   -9.61350674
#> [451,]    5.80258556   1.399548465      22.03553    1.25780086  -11.05328755
#> [452,]   12.49407595   0.453992462      21.15539    1.99241300  -10.27171421
#> [453,]    1.66846230   1.840979342      22.47856    2.58492364  -12.16967420
#> [454,]    4.00230434   3.171289822      22.82364    0.56964242  -11.97345703
#> [455,]    3.24186923   0.057547599      21.71884    2.60328562  -11.35635477
#> [456,]   10.80056988   0.119929744      20.91251    1.93283783   -9.81882570
#> [457,]    0.82193054   2.617444459      22.42461    1.73856696  -11.94088594
#> [458,]    4.15580531   2.421754629      22.98506    1.66316979  -12.33536198
#> [459,]   16.73616486  -1.122702195      16.15133    1.00109752    4.32366467
#> [460,]    5.19223346   4.367723031      21.30374    2.44691410  -11.34335256
#> [461,]    9.35281492   1.592310758      21.75681    1.54684109  -11.10977254
#> [462,]    8.77914929   2.272152325      22.39019    1.03235761  -11.80411876
#> [463,]   13.52233058   1.725269183      19.91136    1.42977151   -7.30644719
#> [464,]    0.52214011   0.532593252      23.59967    2.05119182  -13.14805178
#> [465,]    5.38101079   3.283875024      22.33531    2.07245427  -12.06772912
#> [466,]    5.56281314   2.381102212      22.49869    1.18010701  -12.03448964
#> [467,]    4.92268868   1.981446253      22.60119    1.75101596  -11.66968339
#> [468,]    9.18586136   3.543473408      21.78951    1.24056591  -10.75858319
#> [469,]   10.24680494   2.490835835      21.22808    1.11620724  -10.41040595
#> [470,]   15.83265508   2.285211667      19.52373    1.05335742   -7.52940217
#> [471,]    9.76046943   2.119798263      21.17310    1.51152932   -9.55283060
#> [472,]    3.00767154   2.986571616      21.82758    3.34069637  -11.93259634
#> [473,]    1.53476785   2.027860926      23.45563    1.97125528  -12.85349034
#> [474,]    3.88539259   4.429504018      20.71544    1.40766777  -10.59507632
#> [475,]    6.35825849   0.141541472      22.44740    1.96670514  -11.78338850
#> [476,]   15.82559676   1.165918953      18.56361    0.98397489   -4.23469951
#> [477,]   19.36978937   1.944275286      17.02645    1.02637428   -2.17513551
#> [478,]   15.06285309   3.077674490      19.10366    0.68540355   -5.59885949
#> [479,]    2.81606960   1.673651861      19.88235    1.30446200  -10.37860972
#> [480,]    6.18245191   3.773107014      22.13989    1.71091737  -11.78395834
#> [481,]   17.03941836   2.499611205      18.77674    0.63884920   -6.12450695
#> [482,]    6.69035338   1.552237278      22.29136    1.81584837  -11.81717299
#> [483,]    2.99764026   1.062381010      22.19779    3.35713966  -12.14721096
#> [484,]   -2.29639271  -0.318097602      22.45047    1.82479573  -12.57847413
#> [485,]   11.29431552   3.198734394      20.77966    0.89146798   -9.61987278
#> [486,]    3.45341053  -0.831707065      22.78343    2.55663580  -12.39712403
#> [487,]   11.53118788   1.656140936      20.16655    1.36406565   -7.19362110
#> [488,]   -1.20475932   3.619025546      20.18903    3.93507350  -11.21452278
#> [489,]    3.60332183   0.216159778      21.68271    1.33946275  -11.71075298
#> [490,]   17.23734168   2.556321856      19.08234    0.60852668   -7.31436595
#> [491,]   -5.04584259  -0.082337540      17.34437    1.06077999   -9.76103396
#> [492,]    0.01567990  -1.884164850      23.30620    1.94328316  -12.98960065
#> [493,]    1.52156871   3.705680163      22.36253    0.76828302  -11.62549551
#> [494,]    1.28196555   1.639068472      22.89178    2.04735685  -12.54282003
#> [495,]    9.12879234   0.675047936      22.13801    1.44294689  -11.06967159
#> [496,]    9.58498751   2.115009770      21.60194    0.95606100  -10.03607694
#> [497,]    1.00112703   0.857943933      20.37020    3.22132378  -11.23727231
#> [498,]    9.38761245   3.112198168      20.59349    0.73258719   -9.14964455
#> [499,]    4.18008816   4.185549888      22.64114    1.50383214  -12.16604230
#> [500,]   12.20536862   2.041870311      20.69844    1.22001082   -8.95423198
#> [501,]   11.33510460   3.983094508      21.03598    0.46930092   -9.66734962
#> [502,]   10.24113775   2.539368599      21.46532    1.70624384  -10.19103388
#> [503,]    6.82157708   2.783504687      20.41516    0.48979087  -10.37946653
#> [504,]   10.44386142   3.995338964      20.76593    0.61537229   -7.75927303
#> [505,]   12.13092568   6.218535471      20.51143    0.17277155   -9.16136429
#> [506,]    3.08558179   1.880081912      21.82397    3.06542564  -11.88537048
#> [507,]   10.55855751   3.079836837      21.37459    0.92077421  -10.11096104
#> [508,]    9.30122069   2.510249112      21.29428    1.15053942  -10.24897810
#> [509,]    0.07304677   2.885580079      22.02103    3.12680077  -12.26248924
#> [510,]   14.31515262   2.224741680      19.24019    1.29675037   -5.05977486
#> [511,]    5.69500111   4.746005863      21.12956    1.63313640  -10.87877088
#> [512,]    1.11240814   2.102161569      22.41514    2.43578053  -12.28270459
#> [513,]    2.27419863  -0.396389531      22.28657    1.81192706  -12.27153876
#> [514,]    9.13636708   1.129436940      21.51793    2.25389590   -9.57881479
#> [515,]   11.43386607  -0.464762279      21.00870    2.07796028   -9.75391364
#> [516,]    4.99754112   1.027013260      21.70480    2.22103534  -11.37634130
#> [517,]    6.64434010   5.633558205      21.58132    0.47386498  -10.70674508
#> [518,]    5.21493901   3.723680369      21.47315    2.08655171  -10.94381098
#> [519,]    1.30475187  -2.938260395      23.40691    1.85873264  -13.06402164
#> [520,]    7.60362303   1.455319829      21.95783    1.39217623  -10.66848625
#>        wetness_v3.l2
#>   [1,]   1.943324604
#>   [2,]   0.405054553
#>   [3,]  -1.511110258
#>   [4,]   0.693656246
#>   [5,]  -1.296217257
#>   [6,]  -0.529109818
#>   [7,]   0.983871843
#>   [8,]   1.046988348
#>   [9,]  -2.039827107
#>  [10,]  -1.384356480
#>  [11,]   1.493337472
#>  [12,]  -0.564946592
#>  [13,]  -0.379611917
#>  [14,]   0.402981443
#>  [15,]   3.930750599
#>  [16,]   0.459008205
#>  [17,]  -2.417739555
#>  [18,]  -0.255831452
#>  [19,]  -0.979833503
#>  [20,]  -0.418200774
#>  [21,]  -0.737246330
#>  [22,]  -0.950975331
#>  [23,]  -0.476594682
#>  [24,]   1.539224312
#>  [25,]  -0.093449928
#>  [26,]  -0.766537404
#>  [27,]  -2.034123108
#>  [28,]   0.441276047
#>  [29,]   1.578612842
#>  [30,]  -1.099732095
#>  [31,]  -1.519410526
#>  [32,]  -0.545286682
#>  [33,]   0.749418345
#>  [34,]  -1.244409434
#>  [35,]  -1.400042561
#>  [36,]  -0.289048686
#>  [37,]  -1.560867376
#>  [38,]   1.173229322
#>  [39,]   4.743686406
#>  [40,]   0.231071968
#>  [41,]   5.759114576
#>  [42,]  -0.313127738
#>  [43,]  -0.916569810
#>  [44,]   3.175021765
#>  [45,]  -0.927419631
#>  [46,]   1.023904274
#>  [47,]   1.173141204
#>  [48,]   1.553461710
#>  [49,]  -1.288762215
#>  [50,]   0.050879019
#>  [51,]   1.363560376
#>  [52,]  -1.367921977
#>  [53,]   1.717727668
#>  [54,]   0.995213372
#>  [55,]  -1.475815646
#>  [56,]   0.401195996
#>  [57,]  -0.931952908
#>  [58,]  -0.736989536
#>  [59,]   1.118420772
#>  [60,]  -0.656496761
#>  [61,]  -0.816903697
#>  [62,]   1.405903320
#>  [63,]  -2.573956083
#>  [64,]  -1.116538941
#>  [65,]  -1.145816300
#>  [66,]   0.033597085
#>  [67,]  -1.060718994
#>  [68,]   0.321893701
#>  [69,]  -1.020878726
#>  [70,]  -0.320542995
#>  [71,]   0.007165944
#>  [72,]   0.957284653
#>  [73,]  -0.811414063
#>  [74,]  -0.923216606
#>  [75,]   2.703273308
#>  [76,]  -0.823587230
#>  [77,]  -0.482776779
#>  [78,]  -0.435728481
#>  [79,]  -0.814845257
#>  [80,]   0.634578570
#>  [81,]  -1.227383595
#>  [82,]  -1.755050296
#>  [83,]   4.223127578
#>  [84,]   0.548793097
#>  [85,]  -0.446446913
#>  [86,]  -1.193998460
#>  [87,]  -0.335786091
#>  [88,]  -0.588456749
#>  [89,]  -0.762485347
#>  [90,]  -1.390650380
#>  [91,]   0.300738460
#>  [92,]  -1.637934358
#>  [93,]  -0.826676117
#>  [94,]  -1.190305376
#>  [95,]   3.148200182
#>  [96,]   0.410552212
#>  [97,]  -1.715711784
#>  [98,]   3.028078291
#>  [99,]  -1.023001664
#> [100,]  -0.460769445
#> [101,]   0.008945283
#> [102,]  -1.739508393
#> [103,]  -0.754276516
#> [104,]  -0.803674892
#> [105,]   0.487117665
#> [106,]   0.888636018
#> [107,]  -0.095285870
#> [108,]  -0.637210203
#> [109,]   1.971950685
#> [110,]  -1.658994911
#> [111,]  -0.867383575
#> [112,]  -1.202116093
#> [113,]  -1.095373045
#> [114,]  -1.158230583
#> [115,]   0.244825128
#> [116,]  -0.659763026
#> [117,]  -1.308917824
#> [118,]   0.067955608
#> [119,]  -0.481077473
#> [120,]   0.107477732
#> [121,]   1.847358868
#> [122,]  -1.648367754
#> [123,]  -1.103299205
#> [124,]  -0.593794874
#> [125,]   1.202825449
#> [126,]   0.985225054
#> [127,]   1.250183101
#> [128,]  -1.724618799
#> [129,]  -0.156039447
#> [130,]   0.784253063
#> [131,]  -0.474207066
#> [132,]  -2.464160068
#> [133,]  -1.175699430
#> [134,]   0.215841540
#> [135,]  -1.191073101
#> [136,]  -1.303251686
#> [137,]  -1.099751585
#> [138,]  -1.025138682
#> [139,]  -1.222470383
#> [140,]   1.694404254
#> [141,]  -1.776455125
#> [142,]   1.954887407
#> [143,]  -1.693836407
#> [144,]   1.072696853
#> [145,]  -0.662834942
#> [146,]   2.824830615
#> [147,]  -1.169603383
#> [148,]  -1.977250240
#> [149,]  -1.294644397
#> [150,]   1.866156983
#> [151,]  -0.131527629
#> [152,]  -0.924277730
#> [153,]   0.794645147
#> [154,]   2.127525342
#> [155,]  -1.807451539
#> [156,]   0.820763828
#> [157,]  -1.165279013
#> [158,]  -1.109321616
#> [159,]   0.934920262
#> [160,]  -1.886477604
#> [161,]  -1.019171358
#> [162,]   0.791325613
#> [163,]   3.967449782
#> [164,]  -0.684434341
#> [165,]  -0.580157472
#> [166,]   0.387204223
#> [167,]  -0.335018749
#> [168,]  -2.953745574
#> [169,]  -0.111209315
#> [170,]   0.228447564
#> [171,]  -1.728473109
#> [172,]   0.418756703
#> [173,]  -0.701136716
#> [174,]  -0.612587016
#> [175,]  -0.829631597
#> [176,]  -0.363355595
#> [177,]  -1.318288764
#> [178,]  -0.851240564
#> [179,]  -1.134153692
#> [180,]  -0.066539207
#> [181,]  -0.938834768
#> [182,]  -0.561636935
#> [183,]   3.310819389
#> [184,]  -0.865841148
#> [185,]  -0.003550336
#> [186,]   0.554407562
#> [187,]   0.103373229
#> [188,]   3.487142088
#> [189,]  -0.354539621
#> [190,]  -1.078949591
#> [191,]  -0.956749846
#> [192,]   0.141577288
#> [193,]  -0.993607769
#> [194,]  -1.045570054
#> [195,]   1.427913460
#> [196,]  -0.956121125
#> [197,]   0.559917315
#> [198,]  -0.962027574
#> [199,]  -1.557763652
#> [200,]  -0.463270725
#> [201,]  -1.636012592
#> [202,]  -0.265396798
#> [203,]  -2.270142144
#> [204,]   0.413967375
#> [205,]  -0.518440783
#> [206,]  -0.633676484
#> [207,]  -1.598376894
#> [208,]   2.519755306
#> [209,]  -1.910941545
#> [210,]  -1.135683463
#> [211,]   1.733319060
#> [212,]  -1.157523592
#> [213,]  -1.141675540
#> [214,]   0.075433596
#> [215,]  -1.243221536
#> [216,]  -0.083674969
#> [217,]  -1.498645755
#> [218,]  -1.545344648
#> [219,]   0.448858835
#> [220,]   0.403404043
#> [221,]  -0.680597774
#> [222,]   0.210927812
#> [223,]  -0.629440893
#> [224,]  -1.241926995
#> [225,]   0.669457943
#> [226,]  -1.565964987
#> [227,]  -0.432422411
#> [228,]  -1.099840204
#> [229,]  -1.601837181
#> [230,]   4.325413708
#> [231,]   3.099475534
#> [232,]  -0.883827195
#> [233,]  -1.485275003
#> [234,]  -1.138509755
#> [235,]  -0.743626692
#> [236,]  -1.555459591
#> [237,]   0.035104589
#> [238,]  -1.237104817
#> [239,]  -0.537952298
#> [240,]   0.607341884
#> [241,]  -1.539279111
#> [242,]   1.688752689
#> [243,]  -1.800320399
#> [244,]  -0.357651983
#> [245,]  -1.323367889
#> [246,]  -1.026594993
#> [247,]  -0.946362861
#> [248,]   0.574222537
#> [249,]  -0.889666181
#> [250,]  -0.691720652
#> [251,]  -0.618524151
#> [252,]  -0.799712032
#> [253,]  -1.241670163
#> [254,]  -1.096409384
#> [255,]  -1.049553808
#> [256,]  -1.058559211
#> [257,]  -2.763593717
#> [258,]  -0.008292590
#> [259,]  -0.434564126
#> [260,]  -0.187875764
#> [261,]  -1.042412973
#> [262,]  -1.177485992
#> [263,]   0.325055626
#> [264,]  -1.215933093
#> [265,]   2.893314327
#> [266,]   0.059827240
#> [267,]  -1.032364964
#> [268,]  -1.207198257
#> [269,]  -0.929784697
#> [270,]   0.357952450
#> [271,]  -0.835873341
#> [272,]  -0.211255372
#> [273,]  -0.351987769
#> [274,]  -0.890526103
#> [275,]  -0.565260041
#> [276,]  -1.416440281
#> [277,]  -0.871342104
#> [278,]  -1.533171420
#> [279,]  -1.673810355
#> [280,]  -0.307360811
#> [281,]   1.675661014
#> [282,]   1.005291804
#> [283,]  -0.184255344
#> [284,]  -0.661254937
#> [285,]   0.441868619
#> [286,]   0.684996345
#> [287,]   2.691307888
#> [288,]   6.349818127
#> [289,]   2.086391249
#> [290,]   0.622573194
#> [291,]  -1.644109532
#> [292,]  -1.206810642
#> [293,]  -0.109001449
#> [294,]   1.403326251
#> [295,]  -0.135109146
#> [296,]  -0.774521592
#> [297,]   1.513232096
#> [298,]  -0.680159300
#> [299,]  -0.351231768
#> [300,]   1.252780104
#> [301,]   3.019640079
#> [302,]  -0.803385267
#> [303,]  -0.095865950
#> [304,]  -1.263073292
#> [305,]  -0.040259763
#> [306,]  -0.135517904
#> [307,]   0.410358534
#> [308,]   0.387915735
#> [309,]   0.633309525
#> [310,]  -0.746464987
#> [311,]  -0.980062562
#> [312,]   1.772362587
#> [313,]  -0.102961877
#> [314,]  -0.929122913
#> [315,]   1.726707655
#> [316,]  -0.509729489
#> [317,]  -1.565133974
#> [318,]  -2.923960436
#> [319,]   1.034359829
#> [320,]  -1.501840013
#> [321,]  -1.812271664
#> [322,]   0.794609729
#> [323,]   2.389822062
#> [324,]  -1.312401410
#> [325,]  -2.941749241
#> [326,]  -0.744831868
#> [327,]  -1.626039457
#> [328,]  -0.207353832
#> [329,]   0.050816639
#> [330,]   1.304921628
#> [331,]   0.400402044
#> [332,]  -0.694482564
#> [333,]  -0.130554181
#> [334,]  -0.497881223
#> [335,]   0.139140812
#> [336,]  -0.835613530
#> [337,]  -0.831048359
#> [338,]  -0.753001368
#> [339,]   0.133025343
#> [340,]  -0.549704742
#> [341,]  -1.221475672
#> [342,]  -1.201707743
#> [343,]   1.896479736
#> [344,]  -0.863362145
#> [345,]  -0.225020899
#> [346,]  -1.487869563
#> [347,]  -1.097370974
#> [348,]   0.470004583
#> [349,]  -0.639206925
#> [350,]  -2.101043958
#> [351,]  -1.115062887
#> [352,]   3.748547137
#> [353,]  -0.709083264
#> [354,]   1.160696317
#> [355,]  -0.707107616
#> [356,]  -0.642103244
#> [357,]  -1.299735213
#> [358,]  -0.425851581
#> [359,]  -1.878129004
#> [360,]   0.338169694
#> [361,]  -1.491437143
#> [362,]   0.929963757
#> [363,]  -1.335707578
#> [364,]  -0.657822876
#> [365,]  -0.434036955
#> [366,]  -0.429442623
#> [367,]   0.091763782
#> [368,]  -0.517030967
#> [369,]  -1.091984004
#> [370,]  -0.020254327
#> [371,]   1.343667395
#> [372,]  -1.505599484
#> [373,]  -1.188881576
#> [374,]  -1.468210752
#> [375,]   0.036717626
#> [376,]  -0.556314171
#> [377,]  -0.761654820
#> [378,]  -2.214936715
#> [379,]   0.225502377
#> [380,]  -0.256021101
#> [381,]   0.358656185
#> [382,]  -0.073368481
#> [383,]   0.714848077
#> [384,]  -0.474319239
#> [385,]  -0.794949261
#> [386,]  -0.022658536
#> [387,]  -0.130224978
#> [388,]   0.368769262
#> [389,]  -0.688123443
#> [390,]  -1.061312961
#> [391,]  -0.488301443
#> [392,]  -0.666053627
#> [393,]  -1.437550013
#> [394,]  -0.997327513
#> [395,]  -0.858162794
#> [396,]  -2.009340531
#> [397,]   0.007641818
#> [398,]  -0.537231846
#> [399,]  -0.593316985
#> [400,]   0.248723607
#> [401,]   1.449432280
#> [402,]   1.184349614
#> [403,]  -1.572353506
#> [404,]   0.973442488
#> [405,]   3.614899111
#> [406,]   1.469535097
#> [407,]  -2.761947157
#> [408,]   0.149452277
#> [409,]   2.345341465
#> [410,]  -1.013979235
#> [411,]   0.419844381
#> [412,]  -0.969165215
#> [413,]  -0.891835423
#> [414,]  -0.780495387
#> [415,]  -1.227559391
#> [416,]   0.445129343
#> [417,]   0.996652367
#> [418,]  -1.615906134
#> [419,]  -0.335396300
#> [420,]   0.729380552
#> [421,]  -0.750079461
#> [422,]  -1.122161234
#> [423,]  -0.785518325
#> [424,]  -1.081285706
#> [425,]  -1.452969296
#> [426,]  -1.277643778
#> [427,]   1.945657628
#> [428,]  -0.663789508
#> [429,]   0.558718003
#> [430,]   0.733876204
#> [431,]   0.739793775
#> [432,]  -0.429041330
#> [433,]  -1.982900555
#> [434,]  -0.489879089
#> [435,]  -1.186566079
#> [436,]   0.219881168
#> [437,]  -1.588392119
#> [438,]  -0.859873407
#> [439,]   0.139394070
#> [440,]  -1.071870807
#> [441,]  -0.748598150
#> [442,]   0.644975360
#> [443,]  -1.781292504
#> [444,]   1.350648018
#> [445,]   0.018456977
#> [446,]  -1.465469832
#> [447,]  -1.583411460
#> [448,]  -1.836902033
#> [449,]  -1.796649038
#> [450,]  -1.023991340
#> [451,]  -0.064886049
#> [452,]  -1.304257443
#> [453,]  -1.280575373
#> [454,]   0.609616248
#> [455,]  -1.621271140
#> [456,]  -1.123568799
#> [457,]  -0.364252705
#> [458,]  -0.799837488
#> [459,]   4.427613660
#> [460,]  -1.082186765
#> [461,]  -0.718864409
#> [462,]  -0.297873612
#> [463,]  -0.282919792
#> [464,]  -1.154198086
#> [465,]  -0.842103484
#> [466,]  -0.394073223
#> [467,]  -1.065880015
#> [468,]  -0.319659605
#> [469,]  -0.192677159
#> [470,]   0.068806454
#> [471,]  -0.208234356
#> [472,]  -1.833460541
#> [473,]  -0.811098065
#> [474,]  -0.143371660
#> [475,]  -1.026073925
#> [476,]   1.551784057
#> [477,]   0.597854002
#> [478,]   0.961138550
#> [479,]  -0.844590294
#> [480,]  -0.436274138
#> [481,]   1.417353486
#> [482,]  -1.077132859
#> [483,]  -1.909390200
#> [484,]  -1.026871124
#> [485,]   0.443056925
#> [486,]  -1.581563903
#> [487,]  -0.717105122
#> [488,]  -2.059704107
#> [489,]  -0.749760489
#> [490,]   1.431825143
#> [491,]  -0.595337440
#> [492,]  -1.139433969
#> [493,]   0.405478105
#> [494,]  -1.076622351
#> [495,]  -0.241795979
#> [496,]   0.546608485
#> [497,]  -1.784862095
#> [498,]   0.701399374
#> [499,]  -0.427700617
#> [500,]  -0.087747502
#> [501,]   0.837501118
#> [502,]  -1.135529395
#> [503,]   0.056952809
#> [504,]   0.964890430
#> [505,]   1.037613496
#> [506,]  -1.716970583
#> [507,]   0.174490011
#> [508,]  -0.506015742
#> [509,]  -1.650969288
#> [510,]  -0.598439981
#> [511,]  -0.263539498
#> [512,]  -1.224913464
#> [513,]  -1.101757600
#> [514,]  -1.890463293
#> [515,]  -1.341966240
#> [516,]  -0.865847156
#> [517,]   1.079705753
#> [518,]  -1.093609394
#> [519,]  -1.089647625
#> [520,]  -0.624675630
#> attr(,"df")
#> [1] 3 2
#> attr(,"range")
#> [1]  0 24
#> attr(,"lag")
#> [1]  0 85
#> attr(,"argvar")
#> attr(,"argvar")$fun
#> [1] "ns"
#> 
#> attr(,"argvar")$knots
#> [1]  8.75715 12.16244
#> 
#> attr(,"argvar")$intercept
#> [1] FALSE
#> 
#> attr(,"argvar")$Boundary.knots
#> [1]  0 24
#> 
#> attr(,"arglag")
#> attr(,"arglag")$fun
#> [1] "ns"
#> 
#> attr(,"arglag")$knots
#> numeric(0)
#> 
#> attr(,"arglag")$intercept
#> [1] TRUE
#> 
#> attr(,"arglag")$Boundary.knots
#> [1]  0 85
#> 
#> attr(,"class")
#> [1] "crossbasis" "matrix"
```

Posterior coefficient draws:

``` r


dim(fit_bdlnm$coefficients)
#> [1]   19 1000

head(fit_bdlnm$coefficients)
#>                  sample1     sample2      sample3     sample4     sample5
#> (Intercept) -3.779096979 -2.08849885 -3.460286325 -1.63904393 -2.84145063
#> tmean_v1.l1  0.094815319  0.07586089  0.097962142  0.07905103  0.08405857
#> tmean_v1.l2 -0.008261747 -0.01318477  0.005568889 -0.02739723  0.01394428
#> tmean_v2.l1  0.137619978  0.08919768  0.143689440  0.07628214  0.12185522
#> tmean_v2.l2  0.261685086  0.15390735  0.303849446  0.11882935  0.21243926
#> tmean_v3.l1 -0.050784853 -0.04137848 -0.027766968 -0.03469468 -0.01633367
#>                  sample6     sample7     sample8     sample9      sample10
#> (Intercept) -3.797347857 -1.55588840 -2.93079969 -1.38250786 -1.9486317049
#> tmean_v1.l1  0.086459374  0.07302896  0.07594046  0.04982782  0.0860304029
#> tmean_v1.l2 -0.002355427 -0.01290876 -0.01835669 -0.04001902  0.0003176844
#> tmean_v2.l1  0.186405869  0.05704914  0.13747745  0.06483365  0.1031298995
#> tmean_v2.l2  0.177874248  0.17357299  0.15794983  0.03160396  0.1260306585
#> tmean_v3.l1  0.016206920 -0.03541640 -0.01716999 -0.03256766 -0.0128545351
#>                 sample11    sample12    sample13     sample14     sample15
#> (Intercept) -4.475315230 -4.80969093 -4.11260390 -2.486539018 -2.713984131
#> tmean_v1.l1  0.103744991  0.10304503  0.10354053  0.090854748  0.076708456
#> tmean_v1.l2  0.028266498  0.01074011  0.04037308  0.004182017  0.016636407
#> tmean_v2.l1  0.164275106  0.21904065  0.15341175  0.131774223  0.123438718
#> tmean_v2.l2  0.293908788  0.25571473  0.32748678  0.151555297  0.211323552
#> tmean_v3.l1 -0.009580847  0.01653628 -0.02125101  0.001549944 -0.007237924
#>                 sample16     sample17     sample18    sample19     sample20
#> (Intercept) -3.673304917 -2.758208915 -4.352563344 -3.27142814 -1.279481897
#> tmean_v1.l1  0.098010334  0.088334706  0.088775631  0.08499615  0.057216344
#> tmean_v1.l2  0.005196725 -0.008928737 -0.005310461 -0.01499579 -0.022712802
#> tmean_v2.l1  0.154056448  0.112412866  0.189827455  0.11936352  0.080872246
#> tmean_v2.l2  0.220093402  0.228973657  0.222412215  0.19859496  0.049322481
#> tmean_v3.l1 -0.017654457 -0.023504625 -0.020713024 -0.03780698 -0.008094865
#>                 sample21    sample22     sample23     sample24    sample25
#> (Intercept) -4.129204799 -2.08580783 -2.850003266 -3.771438228 -1.16638008
#> tmean_v1.l1  0.081052228  0.06946573  0.088399342  0.091050460  0.07190738
#> tmean_v1.l2  0.010990400 -0.02943281  0.004096972 -0.001981254 -0.01569092
#> tmean_v2.l1  0.178668408  0.10828385  0.157601551  0.177517012  0.06609011
#> tmean_v2.l2  0.191250737  0.12461100  0.223772930  0.171675854  0.12518726
#> tmean_v3.l1  0.006589239 -0.01911304 -0.008656836  0.015193596 -0.03182509
#>                 sample26    sample27    sample28    sample29    sample30
#> (Intercept) -2.657651503 -1.66796052 -2.81993271 -2.02152492 -0.80069553
#> tmean_v1.l1  0.080030450  0.06456126  0.09316526  0.07477004  0.05979767
#> tmean_v1.l2 -0.009550179 -0.02231559  0.01290745 -0.04180232 -0.05098432
#> tmean_v2.l1  0.109993611  0.09996979  0.09985550  0.07566076  0.04068960
#> tmean_v2.l2  0.206446923  0.08988498  0.26423822  0.10762365 -0.03402899
#> tmean_v3.l1 -0.036112738 -0.01152298 -0.05050839 -0.04840361 -0.03366453
#>                sample31    sample32     sample33    sample34    sample35
#> (Intercept) -3.22088354 -2.59332326 -3.727814822 -5.20498555 -4.16031016
#> tmean_v1.l1  0.10090004  0.08225127  0.081716426  0.11599873  0.09311281
#> tmean_v1.l2  0.01327027 -0.01991379  0.003357016  0.04158000  0.01640278
#> tmean_v2.l1  0.13539728  0.11795623  0.166907727  0.20935717  0.15440034
#> tmean_v2.l2  0.24885236  0.10642993  0.182982198  0.42675903  0.28534817
#> tmean_v3.l1 -0.02948983 -0.01698778  0.018429970 -0.02458307 -0.02711961
#>                sample36    sample37     sample38    sample39    sample40
#> (Intercept) -3.33321533 -3.18626623 -3.437078962 -1.93021369 -2.81606801
#> tmean_v1.l1  0.08990430  0.08378313  0.085467378  0.08298317  0.09006512
#> tmean_v1.l2 -0.01373356 -0.02217271 -0.005946555 -0.01523790 -0.02116935
#> tmean_v2.l1  0.13041667  0.12903210  0.126183914  0.07317161  0.13045047
#> tmean_v2.l2  0.24883093  0.14503895  0.169329615  0.16315101  0.15665439
#> tmean_v3.l1 -0.04139460 -0.02544903 -0.018604002 -0.03994608 -0.02478738
#>                 sample41      sample42    sample43     sample44      sample45
#> (Intercept) -3.611033003 -1.308005e+00 -3.02224408 -2.621309677 -4.5098632318
#> tmean_v1.l1  0.090271329  5.859741e-02  0.08869797  0.086271756  0.1060296944
#> tmean_v1.l2  0.024987237 -6.058383e-02  0.02499724 -0.001422973  0.0617011834
#> tmean_v2.l1  0.142646332  1.015597e-01  0.10808281  0.091250379  0.1651711247
#> tmean_v2.l2  0.224450039 -8.700444e-05  0.23846400  0.188538860  0.3846551589
#> tmean_v3.l1 -0.005132264 -1.008183e-02 -0.01242916 -0.023211181 -0.0008217666
#>                sample46    sample47    sample48     sample49    sample50
#> (Intercept) -0.24024339 -1.68722246 -2.34781975 -2.033424556 -4.60162543
#> tmean_v1.l1  0.05541917  0.06477729  0.07796775  0.067731018  0.10134737
#> tmean_v1.l2 -0.04807470 -0.01400435 -0.02247138  0.004314775  0.02417617
#> tmean_v2.l1  0.03902176  0.06766877  0.10324140  0.072022505  0.19366355
#> tmean_v2.l2 -0.06125891  0.09520606  0.14609310  0.137672802  0.33383458
#> tmean_v3.l1 -0.01041795 -0.04651488 -0.02897704 -0.028646285 -0.01920612
#>                sample51    sample52    sample53    sample54     sample55
#> (Intercept) -1.99065238 -1.98498070 -2.62086275 -4.31401519 -3.294605926
#> tmean_v1.l1  0.07205438  0.07729534  0.08757511  0.10126890  0.085226773
#> tmean_v1.l2 -0.01891057 -0.02577651  0.01782448  0.02303171  0.008941194
#> tmean_v2.l1  0.09677163  0.09903799  0.11633102  0.17375377  0.146935488
#> tmean_v2.l2  0.09642624  0.11762562  0.16445291  0.24474614  0.246848538
#> tmean_v3.l1 -0.01613753 -0.02721311 -0.01267040 -0.01213958 -0.015599815
#>                 sample56    sample57    sample58     sample59    sample60
#> (Intercept) -3.524612896 -3.85240704 -2.29368295 -2.160765810 -1.25334684
#> tmean_v1.l1  0.089946219  0.10241047  0.08401615  0.071116425  0.06586679
#> tmean_v1.l2 -0.009203407 -0.01671348  0.01038143 -0.005572887 -0.03252070
#> tmean_v2.l1  0.149460325  0.17562906  0.12793084  0.119439910  0.07592673
#> tmean_v2.l2  0.190566734  0.22822284  0.12759851  0.103035630  0.10835667
#> tmean_v3.l1 -0.013694268 -0.02184762  0.01577805 -0.012811178 -0.04614746
#>                sample61     sample62      sample63     sample64      sample65
#> (Intercept) -2.74565765 -3.284356151 -3.5678632332 -1.438819537 -3.8602445752
#> tmean_v1.l1  0.08024576  0.101580701  0.0957188648  0.075247979  0.0910730234
#> tmean_v1.l2 -0.01999759 -0.005169146  0.0150049320 -0.004890181  0.0006331908
#> tmean_v2.l1  0.13196572  0.121372425  0.1579305175  0.062716966  0.1524804114
#> tmean_v2.l2  0.16360593  0.260955249  0.2390463939  0.177294731  0.2388847325
#> tmean_v3.l1 -0.02475883 -0.035183379 -0.0001418616 -0.052614261 -0.0305327316
#>                sample66    sample67     sample68    sample69    sample70
#> (Intercept) -3.08480913 -2.65906885 -2.314788935 -2.58945176 -3.83499332
#> tmean_v1.l1  0.08559242  0.09257283  0.075337987  0.08145854  0.08163360
#> tmean_v1.l2 -0.01639263 -0.00149441 -0.003962374 -0.02692439  0.02471758
#> tmean_v2.l1  0.15081949  0.11059638  0.117021023  0.10778240  0.14931854
#> tmean_v2.l2  0.17562844  0.20099431  0.143767570  0.07600594  0.28206918
#> tmean_v3.l1 -0.01140864 -0.03375152 -0.006130039 -0.02240655 -0.02073504
#>                sample71    sample72     sample73    sample74     sample75
#> (Intercept) -1.15124735 -2.98853994 -2.132949535 -1.41735504 -0.989572347
#> tmean_v1.l1  0.06771340  0.07833748  0.075678639  0.06464121  0.061102163
#> tmean_v1.l2 -0.02134340  0.00444092  0.008436183 -0.02358660 -0.047719002
#> tmean_v2.l1  0.05799952  0.13486095  0.088597221  0.04606004  0.088799041
#> tmean_v2.l2  0.09889111  0.21568022  0.149518836  0.11843935 -0.003388174
#> tmean_v3.l1 -0.04382275 -0.01759275 -0.008292401 -0.05908038 -0.018286769
#>                sample76    sample77    sample78     sample79    sample80
#> (Intercept) -2.69100005 -2.53467342 -1.82942665 -2.318873091 -2.75039394
#> tmean_v1.l1  0.08761628  0.07588085  0.06948432  0.082047047  0.07610505
#> tmean_v1.l2  0.01209836 -0.01787112 -0.02039286 -0.008740031 -0.01466767
#> tmean_v2.l1  0.12397283  0.12094666  0.09317630  0.096720697  0.11651263
#> tmean_v2.l2  0.20743502  0.14451968  0.08134451  0.175665688  0.14696309
#> tmean_v3.l1 -0.01983145 -0.02333160 -0.02045959 -0.025252076 -0.02876636
#>                sample81     sample82    sample83     sample84    sample85
#> (Intercept) -3.52767683 -2.988814475 -2.07149415 -2.936799217 -2.74598118
#> tmean_v1.l1  0.08692746  0.082576273  0.07147816  0.072399513  0.08685359
#> tmean_v1.l2 -0.03386415 -0.008413759 -0.01572377 -0.031218472  0.01055456
#> tmean_v2.l1  0.14940287  0.140562061  0.08766007  0.156602571  0.10002665
#> tmean_v2.l2  0.11660179  0.212566246  0.14883036  0.109718131  0.23641554
#> tmean_v3.l1 -0.01994058 -0.025942113 -0.03632573 -0.006967335 -0.04744879
#>                 sample86     sample87    sample88    sample89     sample90
#> (Intercept) -2.370299895 -4.255515626 -3.71407637 -3.82745928 -3.197716647
#> tmean_v1.l1  0.081866943  0.104282368  0.08614971  0.09829317  0.097364949
#> tmean_v1.l2  0.009036393  0.014335867 -0.01412775  0.03569337  0.006214827
#> tmean_v2.l1  0.092051167  0.203009733  0.16228983  0.16966592  0.139730675
#> tmean_v2.l2  0.200493448  0.269115226  0.18608032  0.33346836  0.314116380
#> tmean_v3.l1 -0.022263254 -0.009027383 -0.01400429 -0.02215520 -0.031601053
#>                sample91    sample92    sample93     sample94    sample95
#> (Intercept) -2.46465010 -3.05118524 -5.32729035 -4.955503560 -4.50269183
#> tmean_v1.l1  0.09275380  0.07834050  0.10830453  0.102837217  0.09602142
#> tmean_v1.l2  0.02180129 -0.01801835  0.02534238  0.025392537  0.02244915
#> tmean_v2.l1  0.11634951  0.10965452  0.21372676  0.207715251  0.16950337
#> tmean_v2.l2  0.20683562  0.21045085  0.39080632  0.311021677  0.31160108
#> tmean_v3.l1 -0.02357994 -0.04915137 -0.02992451  0.002189059 -0.02605278
#>                 sample96    sample97    sample98     sample99   sample100
#> (Intercept) -2.274375544 -2.23819403 -1.94973076 -4.095334819 -1.37771621
#> tmean_v1.l1  0.071740379  0.06571244  0.08172779  0.096633393  0.07305242
#> tmean_v1.l2 -0.029397088 -0.03801979  0.01372066  0.005786678 -0.03930338
#> tmean_v2.l1  0.116817982  0.10022207  0.09619474  0.173247743  0.07452285
#> tmean_v2.l2  0.092321476  0.07864493  0.21156932  0.263758333  0.07421572
#> tmean_v3.l1 -0.005590808 -0.03064758 -0.02188884 -0.005295341 -0.02476475
#>               sample101    sample102   sample103    sample104   sample105
#> (Intercept) -0.80268772 -2.947752087 -1.27458556 -0.501390552 -2.04736436
#> tmean_v1.l1  0.06654215  0.079802104  0.06683148  0.055017222  0.07545730
#> tmean_v1.l2 -0.05654118  0.002586612 -0.04970806 -0.063090193 -0.01784039
#> tmean_v2.l1  0.05685734  0.127734858  0.06326825  0.027418017  0.07130767
#> tmean_v2.l2 -0.04099920  0.159683151  0.05392275 -0.008524464  0.19454387
#> tmean_v3.l1 -0.01863863 -0.014261718 -0.05184630 -0.045853362 -0.05338530
#>               sample106    sample107   sample108   sample109   sample110
#> (Intercept) -0.63177326 -1.941156867 -3.94878432 -5.05728471 -1.98093623
#> tmean_v1.l1  0.05680902  0.078178858  0.09011016  0.10132808  0.07143715
#> tmean_v1.l2 -0.03052090 -0.045283679  0.03319220  0.05460802 -0.04687311
#> tmean_v2.l1  0.03813446  0.123224194  0.16319340  0.17966906  0.07860109
#> tmean_v2.l2  0.04792184  0.073517043  0.27769363  0.39366801  0.10251460
#> tmean_v3.l1 -0.02980040 -0.007292874 -0.01067947 -0.01026005 -0.05259020
#>                sample111   sample112   sample113    sample114     sample115
#> (Intercept) -0.226805707 -3.51761648 -2.76979546 -2.836729879 -2.792706e+00
#> tmean_v1.l1  0.055342971  0.09297284  0.07878366  0.084732848  8.873727e-02
#> tmean_v1.l2 -0.057736516 -0.01046955  0.01578583 -0.008499294 -3.294198e-05
#> tmean_v2.l1  0.036598406  0.17745856  0.12748196  0.130042750  1.018350e-01
#> tmean_v2.l2 -0.004515145  0.21582942  0.21783695  0.214109344  2.044038e-01
#> tmean_v3.l1 -0.034166966 -0.00907173 -0.00280160 -0.022213075 -3.770416e-02
#>               sample116   sample117   sample118    sample119     sample120
#> (Intercept) -3.22028547 -2.86952746 -1.41205672 -1.164516057 -2.7273777397
#> tmean_v1.l1  0.07994551  0.08305957  0.07623541  0.062378798  0.0819393295
#> tmean_v1.l2 -0.03061102 -0.01676008 -0.03447778 -0.016374632  0.0003071041
#> tmean_v2.l1  0.14311468  0.14371155  0.05607026  0.087432718  0.1137261426
#> tmean_v2.l2  0.17644987  0.16456153  0.11870197  0.033665674  0.2006271648
#> tmean_v3.l1 -0.04649078 -0.03293830 -0.05578983 -0.001803153 -0.0249378897
#>               sample121    sample122   sample123    sample124   sample125
#> (Intercept) -4.08190759 -1.761772546 -3.81604114 -3.272359084 -3.01416537
#> tmean_v1.l1  0.09986363  0.079791707  0.10244907  0.093398441  0.09558498
#> tmean_v1.l2  0.03135951  0.013433601  0.01107021  0.031286460 -0.02989137
#> tmean_v2.l1  0.18462150  0.086136909  0.15916758  0.150786117  0.14914488
#> tmean_v2.l2  0.32757004  0.096328764  0.28119524  0.243815451  0.11027921
#> tmean_v3.l1 -0.00589173 -0.002165154 -0.01867516  0.008695073 -0.01214020
#>               sample126    sample127   sample128   sample129   sample130
#> (Intercept) -1.58302883 -3.798224360 -1.36078349 -1.08735958 -2.32636696
#> tmean_v1.l1  0.06901878  0.088450824  0.07583498  0.06916483  0.09128682
#> tmean_v1.l2 -0.05077525 -0.016512400 -0.03451069 -0.03684089 -0.02676429
#> tmean_v2.l1  0.06367033  0.172348402  0.06018975  0.05872450  0.10755362
#> tmean_v2.l2  0.06700607  0.178832042  0.09887979  0.09327509  0.14362496
#> tmean_v3.l1 -0.03567665  0.001034329 -0.04572207 -0.05132196 -0.03159552
#>               sample131    sample132   sample133   sample134   sample135
#> (Intercept) -0.95055009 -2.842065947 -0.07571771 -2.60811566 -3.55682372
#> tmean_v1.l1  0.06983730  0.080514040  0.05779559  0.08135872  0.09836128
#> tmean_v1.l2 -0.04353463 -0.007516182 -0.02032767 -0.01786233  0.01118496
#> tmean_v2.l1  0.06939697  0.141746780  0.02201095  0.13205015  0.15631938
#> tmean_v2.l2  0.06553224  0.198129493  0.02181917  0.13374006  0.31801341
#> tmean_v3.l1 -0.03387288 -0.017190404 -0.02095698 -0.01622357 -0.03383529
#>               sample136   sample137    sample138   sample139   sample140
#> (Intercept) -1.76326984 -3.40908712 -1.893327848 -3.47232213 -1.57813097
#> tmean_v1.l1  0.06689904  0.08367116  0.077946111  0.08040115  0.07216619
#> tmean_v1.l2 -0.02281506 -0.01872136 -0.032771811 -0.03266087 -0.04475815
#> tmean_v2.l1  0.08112087  0.12276875  0.113142889  0.14331540  0.06515410
#> tmean_v2.l2  0.11100242  0.22442303  0.084819616  0.13477938  0.07080539
#> tmean_v3.l1 -0.02739885 -0.03887508 -0.007364599 -0.02328721 -0.03365825
#>               sample141    sample142   sample143    sample144   sample145
#> (Intercept) -1.42868173 -4.592354425 -2.54069701 -3.349082963 -1.62965154
#> tmean_v1.l1  0.06447990  0.103209634  0.08794358  0.090853962  0.07740914
#> tmean_v1.l2 -0.01393286  0.032389070  0.03160293  0.030403280 -0.03341298
#> tmean_v2.l1  0.08527011  0.184223139  0.09612116  0.145030168  0.08641061
#> tmean_v2.l2  0.05836253  0.310166349  0.23645003  0.245021298  0.12187086
#> tmean_v3.l1 -0.00277967  0.001148334 -0.03113766 -0.004307507 -0.03486983
#>                sample146    sample147    sample148   sample149    sample150
#> (Intercept) -4.565146513 -3.006568524 -0.284950277 -3.86112419 -4.924864767
#> tmean_v1.l1  0.102843531  0.074744539  0.059012872  0.08513891  0.113602136
#> tmean_v1.l2  0.030049102 -0.016712544 -0.065603879  0.03878128  0.056801171
#> tmean_v2.l1  0.198206861  0.137083116  0.033905407  0.14382780  0.193432037
#> tmean_v2.l2  0.302718740  0.108119674  0.004343583  0.27655549  0.391532149
#> tmean_v3.l1 -0.003040919  0.007478121 -0.047974047 -0.01035856  0.007285137
#>                sample151    sample152    sample153    sample154    sample155
#> (Intercept) -2.194954044 -1.677297116 -3.638706634 -3.176193007 -3.181261551
#> tmean_v1.l1  0.076322635  0.072176016  0.084049221  0.095775664  0.080756726
#> tmean_v1.l2  0.009816218 -0.009238616  0.006988005  0.006468828 -0.025286379
#> tmean_v2.l1  0.091916340  0.089272532  0.145427737  0.120670685  0.154004392
#> tmean_v2.l2  0.160210575  0.152990191  0.203633979  0.282046472  0.119715657
#> tmean_v3.l1 -0.006587015 -0.034746719 -0.009834519 -0.049026286 -0.004710737
#>                sample156   sample157   sample158   sample159   sample160
#> (Intercept) -4.081507400 -1.77366634 -0.95960358 -2.87171515 -2.65844072
#> tmean_v1.l1  0.099878466  0.07136933  0.05699122  0.07885035  0.07295547
#> tmean_v1.l2  0.033834272 -0.03926620 -0.02763156 -0.00934365 -0.02144586
#> tmean_v2.l1  0.157924030  0.09052447  0.06725817  0.12682783  0.10378846
#> tmean_v2.l2  0.259779612  0.06479660  0.05121371  0.17829764  0.17574728
#> tmean_v3.l1 -0.006249957 -0.03280280 -0.02591815 -0.02035214 -0.03678486
#>               sample161   sample162    sample163    sample164   sample165
#> (Intercept) -1.84919904 -2.03492028 -3.712237219 -3.634666685 -2.50279302
#> tmean_v1.l1  0.08504777  0.07697193  0.094644437  0.088596484  0.07434121
#> tmean_v1.l2 -0.00408919 -0.04972421  0.014011364 -0.004487537 -0.04916801
#> tmean_v2.l1  0.07007811  0.09862877  0.169424721  0.153246635  0.11470572
#> tmean_v2.l2  0.19999547  0.11212154  0.229606080  0.237711101  0.06092762
#> tmean_v3.l1 -0.04101147 -0.05275599 -0.004758997 -0.019159005 -0.02581333
#>               sample166   sample167   sample168    sample169    sample170
#> (Intercept) -2.11008230 -1.66674380 -0.75337224 -0.190374595 -2.694681506
#> tmean_v1.l1  0.07461965  0.07479966  0.06112907  0.058363876  0.073079727
#> tmean_v1.l2 -0.03300214 -0.01427849 -0.03619301 -0.043332089 -0.047405759
#> tmean_v2.l1  0.08353168  0.07614899  0.04136681  0.007662706  0.158435197
#> tmean_v2.l2  0.14715819  0.14063531  0.05258484  0.017326404  0.111713909
#> tmean_v3.l1 -0.04512240 -0.04722554 -0.03844605 -0.047713421 -0.003762902
#>               sample171    sample172    sample173    sample174    sample175
#> (Intercept) -2.25631161 -2.882477390 -2.588368196 -2.840911807 -3.747354373
#> tmean_v1.l1  0.07209200  0.084655001  0.077075499  0.077409011  0.097455992
#> tmean_v1.l2 -0.05642695 -0.014155829  0.001076071 -0.006073988  0.025446239
#> tmean_v2.l1  0.11048040  0.126527567  0.131634294  0.111145059  0.168435932
#> tmean_v2.l2  0.05746268  0.125021601  0.200206199  0.196959691  0.234123757
#> tmean_v3.l1 -0.03332846 -0.007606943 -0.017725650 -0.026424692  0.003388883
#>               sample176    sample177   sample178   sample179   sample180
#> (Intercept) -2.86585826 -1.606819537 -1.47202091 -3.19100426 -1.94542628
#> tmean_v1.l1  0.08274040  0.074215096  0.08370774  0.10565699  0.06687109
#> tmean_v1.l2  0.01824347 -0.016327381 -0.01420120  0.02169300 -0.02633985
#> tmean_v2.l1  0.12397571  0.085344072  0.06699103  0.13691700  0.09822628
#> tmean_v2.l2  0.22352362  0.064953961  0.14633654  0.27087640  0.09550112
#> tmean_v3.l1 -0.01783821 -0.005520822 -0.04422771 -0.01922206 -0.01911774
#>               sample181    sample182   sample183    sample184    sample185
#> (Intercept) -2.01494484 -3.943965592 -3.06698692 -2.397485961 -2.860107288
#> tmean_v1.l1  0.06860288  0.087069754  0.08356354  0.081884045  0.092361047
#> tmean_v1.l2  0.02731468  0.010019023 -0.03330401 -0.002365899 -0.010247517
#> tmean_v2.l1  0.08546133  0.162938327  0.15792116  0.119305271  0.145036409
#> tmean_v2.l2  0.17828391  0.251726281  0.11536910  0.146582116  0.180975351
#> tmean_v3.l1 -0.01437721 -0.003233656 -0.01247651 -0.008821848 -0.002055731
#>               sample186   sample187   sample188   sample189   sample190
#> (Intercept) -4.07865799 -2.67014335 -1.68801088 -3.56075553 -1.73056413
#> tmean_v1.l1  0.10428176  0.07898203  0.07397676  0.10537875  0.07137175
#> tmean_v1.l2  0.02391261  0.01252470 -0.02096758  0.01859996 -0.03128123
#> tmean_v2.l1  0.17107425  0.09571771  0.08454302  0.12971978  0.07166032
#> tmean_v2.l2  0.29887276  0.18176156  0.11435071  0.32174546  0.06029766
#> tmean_v3.l1 -0.02595074 -0.01671757 -0.02746604 -0.04335093 -0.04054282
#>               sample191   sample192   sample193   sample194    sample195
#> (Intercept) -5.31521715 -1.56028885 -2.34936503 -3.27419826 -0.867037764
#> tmean_v1.l1  0.10813566  0.04533688  0.06551189  0.09637628  0.051309536
#> tmean_v1.l2 -0.01208122 -0.04103528 -0.02488992  0.03913427 -0.046700224
#> tmean_v2.l1  0.20786748  0.06614822  0.07452829  0.12568600  0.032499475
#> tmean_v2.l2  0.32037784  0.01193416  0.10610559  0.32047487 -0.004625752
#> tmean_v3.l1 -0.02684566 -0.01855041 -0.04683058 -0.02500606 -0.041198128
#>                sample196   sample197   sample198    sample199    sample200
#> (Intercept) -3.153238687 -2.30385720 -2.98616203 -3.318543364 -2.134663665
#> tmean_v1.l1  0.084833352  0.08016744  0.08032390  0.091370431  0.083223790
#> tmean_v1.l2  0.002761304 -0.05846974  0.01073375  0.002415158 -0.002755847
#> tmean_v2.l1  0.139321079  0.08672542  0.11917086  0.152151532  0.099656439
#> tmean_v2.l2  0.211474258  0.07520656  0.22933540  0.236362674  0.190902715
#> tmean_v3.l1 -0.018480620 -0.04984695 -0.02163703 -0.008652337 -0.023864960
#>               sample201    sample202    sample203   sample204   sample205
#> (Intercept) -3.71915077 -3.257550565 -4.288865357 -1.41747854 -3.64017888
#> tmean_v1.l1  0.09073077  0.084817472  0.094402153  0.06691471  0.09307455
#> tmean_v1.l2  0.01738866 -0.004084766  0.019429805 -0.01928848  0.00791251
#> tmean_v2.l1  0.14741925  0.121290762  0.195120741  0.08232474  0.10013250
#> tmean_v2.l2  0.28274548  0.212239882  0.253285822  0.06341939  0.27300689
#> tmean_v3.l1 -0.02585795 -0.028152470  0.008977732 -0.01174263 -0.04539233
#>               sample206    sample207   sample208   sample209    sample210
#> (Intercept) -1.06566208 -2.447070965 -1.81466917 -2.14275951 -2.343422316
#> tmean_v1.l1  0.06840736  0.089330036  0.07188559  0.07616915  0.087774239
#> tmean_v1.l2 -0.02665284  0.008263266 -0.01876821 -0.01685770  0.008459528
#> tmean_v2.l1  0.03285829  0.083410891  0.10528583  0.10458481  0.108136232
#> tmean_v2.l2  0.10117109  0.257292504  0.09294637  0.13765937  0.247017125
#> tmean_v3.l1 -0.05708467 -0.048882377 -0.01493845 -0.02301335 -0.027116113
#>               sample211   sample212   sample213    sample214   sample215
#> (Intercept) -4.71362567 -1.68040819 -3.97729927 -5.397584816 -2.07498240
#> tmean_v1.l1  0.09359313  0.07824508  0.08485536  0.105806828  0.08162437
#> tmean_v1.l2  0.02753654 -0.03078305  0.02323172  0.027835317  0.02735507
#> tmean_v2.l1  0.20209538  0.07071215  0.13774136  0.229173772  0.07565815
#> tmean_v2.l2  0.32173432  0.11941318  0.27288728  0.337950320  0.24344772
#> tmean_v3.l1 -0.01211697 -0.04727831 -0.02177397 -0.006836793 -0.03803404
#>                 sample216    sample217     sample218   sample219   sample220
#> (Intercept) -1.6984100930 -2.194822209 -3.5320976173 -2.04114652 -2.80489770
#> tmean_v1.l1  0.0734986292  0.086692113  0.0942554460  0.07915896  0.07959590
#> tmean_v1.l2 -0.0004918809  0.008742069  0.0006073774 -0.04063880 -0.01581401
#> tmean_v2.l1  0.1044052912  0.075748811  0.1343498135  0.11250604  0.11197696
#> tmean_v2.l2  0.1590593317  0.208382346  0.2435470228  0.07780140  0.17325564
#> tmean_v3.l1 -0.0248138452 -0.038700483 -0.0304204367 -0.01296712 -0.04224584
#>               sample221    sample222   sample223    sample224   sample225
#> (Intercept) -0.36518806 -1.135035851 -2.55908633 -3.502218069 -1.95548116
#> tmean_v1.l1  0.06190973  0.055249547  0.07850973  0.095878061  0.07423570
#> tmean_v1.l2 -0.03348651 -0.057088685 -0.00278381 -0.008545691 -0.03843149
#> tmean_v2.l1  0.02367690  0.059410358  0.10169996  0.147436567  0.08561301
#> tmean_v2.l2  0.03120317 -0.001870107  0.16429373  0.209860325  0.13127057
#> tmean_v3.l1 -0.04567207 -0.028600355 -0.02605267 -0.023981392 -0.04499736
#>                 sample226    sample227   sample228    sample229   sample230
#> (Intercept) -2.2912055920 -2.915120652 -2.50902368 -3.075973609 -3.08372517
#> tmean_v1.l1  0.0831463155  0.080072360  0.09551302  0.092529749  0.08891143
#> tmean_v1.l2  0.0119144527 -0.014330311  0.01932527 -0.009418508  0.01863593
#> tmean_v2.l1  0.1171352656  0.148612141  0.09813886  0.143202052  0.11265279
#> tmean_v2.l2  0.1920077782  0.123682335  0.23388348  0.209722748  0.25674794
#> tmean_v3.l1 -0.0008149802 -0.008865919 -0.02481832 -0.022425190 -0.03333979
#>                sample231    sample232   sample233   sample234   sample235
#> (Intercept) -2.852972691 -0.900649585 -1.35630381 -2.76919235 -3.61321755
#> tmean_v1.l1  0.084145970  0.066469087  0.07311717  0.09332670  0.09414138
#> tmean_v1.l2 -0.003101436  0.005845426 -0.03175081 -0.01387823 -0.00822173
#> tmean_v2.l1  0.134297779  0.033144878  0.07214466  0.13402155  0.14600186
#> tmean_v2.l2  0.246980903  0.063743678  0.08950376  0.18716860  0.27401476
#> tmean_v3.l1 -0.026628308 -0.028097762 -0.02899523 -0.02527241 -0.04491859
#>               sample236   sample237    sample238   sample239     sample240
#> (Intercept) -1.87686536 -2.77628260 -2.311354139 -1.59755186 -2.7478807431
#> tmean_v1.l1  0.07536307  0.07959953  0.069587612  0.05702612  0.0950641760
#> tmean_v1.l2  0.01869129 -0.02592026 -0.025772423 -0.06088487 -0.0004279487
#> tmean_v2.l1  0.07842609  0.12418539  0.118241808  0.04342732  0.0926487722
#> tmean_v2.l2  0.20726272  0.08116933  0.077402409  0.03268290  0.2513629460
#> tmean_v3.l1 -0.02426624 -0.02118117 -0.005893115 -0.05382117 -0.0532370931
#>               sample241    sample242   sample243   sample244   sample245
#> (Intercept) -3.16234522 -2.506237163 -2.85302865 -2.71762316 -1.60187631
#> tmean_v1.l1  0.08021464  0.079812462  0.09412934  0.08450771  0.07716354
#> tmean_v1.l2 -0.01272401 -0.003694373  0.03227498 -0.01332661 -0.01409603
#> tmean_v2.l1  0.11374441  0.116428296  0.10660735  0.15027927  0.09138384
#> tmean_v2.l2  0.17869700  0.188524589  0.28353284  0.14030590  0.15345010
#> tmean_v3.l1 -0.03979106 -0.013772385 -0.01844533  0.01112849 -0.02306387
#>                 sample246   sample247   sample248   sample249    sample250
#> (Intercept) -3.4166402324 -2.13665591 -3.29987060 -0.61144747 -3.166250899
#> tmean_v1.l1  0.0922982834  0.05916641  0.08862133  0.06019802  0.084799161
#> tmean_v1.l2 -0.0002962083 -0.02921323  0.05017570 -0.07256126  0.027728910
#> tmean_v2.l1  0.1216581343  0.09992365  0.14565120  0.03977142  0.134796938
#> tmean_v2.l2  0.2150971326  0.05012714  0.26313154 -0.04455051  0.255572095
#> tmean_v3.l1 -0.0285191639 -0.01061494 -0.00241291 -0.04237092 -0.005115722
#>               sample251   sample252   sample253    sample254   sample255
#> (Intercept) -1.30502052 -1.61347160 -1.11434147 -3.920928144 -5.13758089
#> tmean_v1.l1  0.05757972  0.07154131  0.07262623  0.084703564  0.11689244
#> tmean_v1.l2 -0.02842138 -0.01846810 -0.01981542  0.007993125  0.06103426
#> tmean_v2.l1  0.06922094  0.06877563  0.05695770  0.170191088  0.24841346
#> tmean_v2.l2  0.04312211  0.12956070  0.10816400  0.226387848  0.32379781
#> tmean_v3.l1 -0.02791238 -0.04407058 -0.02836962 -0.003993397  0.04476410
#>               sample256    sample257   sample258   sample259    sample260
#> (Intercept) -2.78510779 -1.624430009 -3.55767757 -3.30157787 -2.604251241
#> tmean_v1.l1  0.08531726  0.061792435  0.10211534  0.10639024  0.088961337
#> tmean_v1.l2 -0.03075471 -0.013212711  0.03498986  0.01265593  0.005650889
#> tmean_v2.l1  0.11250052  0.082406414  0.13525063  0.16279214  0.103714723
#> tmean_v2.l2  0.13044847  0.088146302  0.32545119  0.29808612  0.176082103
#> tmean_v3.l1 -0.03224886 -0.003632764 -0.03706221 -0.01891530 -0.029758815
#>               sample261   sample262    sample263   sample264   sample265
#> (Intercept) -2.81720228 -1.03630226 -1.344512510 -2.17040032 -2.75766325
#> tmean_v1.l1  0.08497739  0.06775702  0.072164675  0.08091494  0.08037797
#> tmean_v1.l2  0.00858170 -0.02583082 -0.006971272 -0.03395128 -0.01615718
#> tmean_v2.l1  0.11242406  0.03765747  0.084511308  0.08518277  0.13532600
#> tmean_v2.l2  0.22482946  0.13076250  0.107494732  0.09760972  0.16386370
#> tmean_v3.l1 -0.02198219 -0.05222414 -0.013778971 -0.03862242 -0.02318946
#>               sample266     sample267    sample268   sample269   sample270
#> (Intercept) -2.92396853 -1.8945560713 -4.821682787 -2.61916654 -3.65149039
#> tmean_v1.l1  0.08079424  0.0774197820  0.099367116  0.08117318  0.08225756
#> tmean_v1.l2 -0.03631512 -0.0003322761  0.034001567 -0.01861004 -0.02218090
#> tmean_v2.l1  0.12137544  0.0839461216  0.181765565  0.13953624  0.15612132
#> tmean_v2.l2  0.15279377  0.1906322898  0.299098718  0.14835824  0.20077324
#> tmean_v3.l1 -0.03213478 -0.0304660100 -0.005533117 -0.01228950 -0.01872251
#>               sample271    sample272    sample273     sample274   sample275
#> (Intercept) -2.29776647 -3.701844093 -4.860115881 -2.8942090352 -1.95796542
#> tmean_v1.l1  0.07620329  0.105373889  0.108405745  0.0814030770  0.07841281
#> tmean_v1.l2 -0.03349513  0.000786027  0.029282678  0.0135730865 -0.03836540
#> tmean_v2.l1  0.12327659  0.135983335  0.215157761  0.1395803326  0.07770008
#> tmean_v2.l2  0.11529552  0.286569048  0.314221548  0.2134264763  0.08508080
#> tmean_v3.l1 -0.02715306 -0.038549899  0.001050139  0.0005778607 -0.04057782
#>               sample276     sample277    sample278   sample279    sample280
#> (Intercept)  0.01589878 -3.5652735390 -3.421734424 -0.59158343 -2.655821492
#> tmean_v1.l1  0.05489848  0.0902121057  0.089265083  0.06070244  0.074690830
#> tmean_v1.l2 -0.02538682 -0.0003005896 -0.003779219 -0.04293099 -0.003216308
#> tmean_v2.l1  0.03719311  0.1624444396  0.143745847  0.02585591  0.107176146
#> tmean_v2.l2  0.01529408  0.2373304723  0.226702628 -0.01786236  0.195032582
#> tmean_v3.l1 -0.01503379 -0.0042754056 -0.020486290 -0.04017316 -0.019998221
#>                sample281     sample282    sample283   sample284    sample285
#> (Intercept) -3.617414484 -3.8594469525 -2.314232921 -1.73032668 -4.482325442
#> tmean_v1.l1  0.084661409  0.1042511396  0.094689978  0.06151283  0.105358479
#> tmean_v1.l2 -0.003837665  0.0008147676 -0.003047401 -0.06538850  0.000242814
#> tmean_v2.l1  0.157442551  0.1684797621  0.111940693  0.07256551  0.195965768
#> tmean_v2.l2  0.214564548  0.2366666776  0.218099068  0.05047650  0.254436125
#> tmean_v3.l1 -0.020050805  0.0015722113 -0.029647773 -0.06657029 -0.025204675
#>               sample286    sample287    sample288    sample289   sample290
#> (Intercept) -2.43580380 -3.492790944 -1.540573575 -1.832416287 -2.05852137
#> tmean_v1.l1  0.08805159  0.092779904  0.068368836  0.068566382  0.07920522
#> tmean_v1.l2 -0.02005754  0.002522615 -0.010631296  0.004102971 -0.01114228
#> tmean_v2.l1  0.10957508  0.112482273  0.102803996  0.101040608  0.09202176
#> tmean_v2.l2  0.19014323  0.259424236  0.076907029  0.146530051  0.16321961
#> tmean_v3.l1 -0.03493790 -0.053740985  0.008640956 -0.015480985 -0.02634194
#>               sample291   sample292   sample293   sample294   sample295
#> (Intercept) -1.38200870 -2.05917858 -2.87683622 -2.57312265 -0.98813189
#> tmean_v1.l1  0.07776487  0.06196454  0.08537081  0.07284649  0.06812941
#> tmean_v1.l2 -0.02139920 -0.03584251 -0.04916061 -0.01926486 -0.02147512
#> tmean_v2.l1  0.10079152  0.10183192  0.12405485  0.10523996  0.04159450
#> tmean_v2.l2  0.07741107  0.05679702  0.13278037  0.16447879  0.08775655
#> tmean_v3.l1 -0.01095230 -0.01103579 -0.03034082 -0.03763304 -0.03900209
#>                sample296   sample297   sample298   sample299    sample300
#> (Intercept) -3.892160279 -1.61118555 -3.32199924 -1.38283666 -3.518430342
#> tmean_v1.l1  0.092759703  0.07853593  0.08994475  0.07688402  0.083710475
#> tmean_v1.l2 -0.002329045 -0.02643088  0.03022749 -0.03801726  0.002196663
#> tmean_v2.l1  0.172290583  0.06992496  0.14635921  0.06923631  0.163245910
#> tmean_v2.l2  0.172518420  0.15988486  0.26149493  0.11883515  0.177387122
#> tmean_v3.l1  0.005122102 -0.04497772 -0.01406063 -0.05171132  0.005114716
#>                sample301   sample302    sample303   sample304     sample305
#> (Intercept) -3.003188120 -2.61887936 -3.265900609 -1.54518794 -3.3043401094
#> tmean_v1.l1  0.092496025  0.07719648  0.099395993  0.06939725  0.0891348885
#> tmean_v1.l2 -0.003350913 -0.01531905  0.008838444 -0.02682634  0.0005465786
#> tmean_v2.l1  0.110128704  0.11547459  0.134712195  0.09611185  0.1298114315
#> tmean_v2.l2  0.195405472  0.16536629  0.276796063  0.14001454  0.2141942124
#> tmean_v3.l1 -0.037992521 -0.03125677 -0.023940768 -0.03153908 -0.0213888194
#>                sample306   sample307   sample308   sample309   sample310
#> (Intercept) -4.735254964 -4.70420021  0.70947601 -1.54089161 -1.33344379
#> tmean_v1.l1  0.103979594  0.09928947  0.04709568  0.07202635  0.07579833
#> tmean_v1.l2  0.027103486  0.01449633 -0.04834511 -0.02034903 -0.01605196
#> tmean_v2.l1  0.206085687  0.19038106 -0.03214228  0.08721466  0.05992003
#> tmean_v2.l2  0.283886451  0.27140389 -0.02118821  0.12854611  0.17966228
#> tmean_v3.l1  0.003769796 -0.00941388 -0.05168909 -0.02911875 -0.04934977
#>               sample311   sample312   sample313    sample314   sample315
#> (Intercept) -3.28367039 -0.89959521 -0.85002251 -2.509560051 -2.64821375
#> tmean_v1.l1  0.09040105  0.05847305  0.05956900  0.078226532  0.08396001
#> tmean_v1.l2 -0.02899590 -0.03696850 -0.06363078 -0.006183555 -0.01870489
#> tmean_v2.l1  0.11376058  0.02258923  0.06868481  0.104387917  0.09437902
#> tmean_v2.l2  0.17668226  0.06458826 -0.03035144  0.150192363  0.17540855
#> tmean_v3.l1 -0.04300936 -0.05729002 -0.02576804 -0.023908425 -0.03852225
#>                sample316   sample317    sample318   sample319   sample320
#> (Intercept) -3.885985065 -4.16410050 -3.091037537 -1.62947654 -1.85480210
#> tmean_v1.l1  0.080450544  0.09515105  0.086890151  0.06771431  0.06422315
#> tmean_v1.l2 -0.011064625  0.02536842  0.005652432 -0.03482335 -0.03699289
#> tmean_v2.l1  0.175603277  0.16358969  0.125745346  0.08070410  0.09624664
#> tmean_v2.l2  0.172231702  0.32433930  0.186385880  0.05282334  0.08167202
#> tmean_v3.l1  0.007898981 -0.02982611 -0.019589819 -0.02671238 -0.02550781
#>               sample321    sample322   sample323   sample324    sample325
#> (Intercept) -1.86004396 -5.515364492 -0.84965054 -4.92479780 -2.766823160
#> tmean_v1.l1  0.07339052  0.103822069  0.06500448  0.11195994  0.077285218
#> tmean_v1.l2 -0.01624118  0.027938531 -0.02293164  0.06737964 -0.001447659
#> tmean_v2.l1  0.06852158  0.226559991  0.06369885  0.18121569  0.131153945
#> tmean_v2.l2  0.11560564  0.329723580  0.03915382  0.45571583  0.141704554
#> tmean_v3.l1 -0.03585852  0.007370055 -0.02301562 -0.02741888  0.006520745
#>               sample326    sample327    sample328    sample329   sample330
#> (Intercept) -1.26898401 -2.421686801 -3.079171339 -3.554835247 -0.74687531
#> tmean_v1.l1  0.07369727  0.081779619  0.099738335  0.087516759  0.05215005
#> tmean_v1.l2 -0.03338668 -0.004596269  0.008114039 -0.001562652 -0.03891091
#> tmean_v2.l1  0.07132697  0.124351719  0.123745029  0.160650590  0.03245709
#> tmean_v2.l2  0.07805539  0.133484686  0.193845728  0.151965438  0.06188357
#> tmean_v3.l1 -0.02221801 -0.008919447 -0.022710894  0.014844855 -0.04473858
#>                sample331    sample332    sample333   sample334    sample335
#> (Intercept) -2.113756943 -3.265416192 -4.643634133 -1.65358261 -4.575301116
#> tmean_v1.l1  0.088042517  0.087807122  0.087315548  0.06939319  0.097752142
#> tmean_v1.l2 -0.002789959 -0.009906265  0.008797590 -0.03342198  0.024574453
#> tmean_v2.l1  0.077364265  0.139303126  0.199150388  0.07550997  0.184965722
#> tmean_v2.l2  0.214779489  0.179734400  0.253307468  0.04451205  0.287258821
#> tmean_v3.l1 -0.044558491 -0.013693977 -0.000673334 -0.01844605 -0.002984981
#>                 sample336   sample337    sample338   sample339   sample340
#> (Intercept) -3.5203093166 -1.42965193 -3.279824018 -4.08699267 -2.47488575
#> tmean_v1.l1  0.0882958768  0.05558372  0.088304806  0.10071091  0.09049253
#> tmean_v1.l2 -0.0009778301 -0.05914446  0.009696734  0.04875447  0.02147652
#> tmean_v2.l1  0.1442659933  0.09510552  0.175141872  0.16577028  0.12068642
#> tmean_v2.l2  0.2159136353 -0.03330417  0.229014005  0.35038721  0.23885498
#> tmean_v3.l1 -0.0166978966 -0.01014100 -0.003176400 -0.01553526 -0.01257641
#>               sample341   sample342    sample343    sample344   sample345
#> (Intercept) -4.30149478 -3.21281406 -2.761751003 -2.223083021 -1.48667840
#> tmean_v1.l1  0.10203586  0.09133204  0.090276970  0.076305433  0.06211038
#> tmean_v1.l2 -0.03323590  0.02418702 -0.002601692 -0.023890240 -0.04087554
#> tmean_v2.l1  0.18891191  0.11852934  0.140829558  0.112276764  0.07852240
#> tmean_v2.l2  0.23386614  0.23331448  0.245461772  0.008916363  0.06044620
#> tmean_v3.l1 -0.02527177 -0.01892862 -0.017780393  0.006310629 -0.02574665
#>               sample346    sample347   sample348   sample349   sample350
#> (Intercept) -3.80261179 -2.873195130 -1.25215958 -2.72494837 -4.31131590
#> tmean_v1.l1  0.09896771  0.086357109  0.06850796  0.08062756  0.10863350
#> tmean_v1.l2 -0.01634431 -0.006573806 -0.02240855  0.01428622  0.01089741
#> tmean_v2.l1  0.16455692  0.140095026  0.07166442  0.07329937  0.18194068
#> tmean_v2.l2  0.22453241  0.153667476  0.10763477  0.25644969  0.29047930
#> tmean_v3.l1 -0.02128263 -0.005902295 -0.02561934 -0.04663334 -0.01792392
#>                 sample351   sample352    sample353    sample354   sample355
#> (Intercept) -1.5578592662 -3.72540524 -3.391617825 -2.427678185 -5.44093820
#> tmean_v1.l1  0.0623732306  0.09195611  0.091692745  0.071679888  0.10880928
#> tmean_v1.l2 -0.0397813400  0.02419451 -0.028160997 -0.008527424  0.05466528
#> tmean_v2.l1  0.0831864412  0.16066735  0.168973075  0.110821543  0.21685154
#> tmean_v2.l2  0.0005307081  0.30298386  0.140696070  0.153024361  0.40440821
#> tmean_v3.l1 -0.0211546168 -0.01562253 -0.008890721 -0.019822865 -0.01179513
#>               sample356   sample357    sample358   sample359   sample360
#> (Intercept) -3.18665542 -2.11847207 -2.569452120 -3.54515210 -3.70324454
#> tmean_v1.l1  0.08004907  0.08462793  0.083817510  0.09109186  0.10372484
#> tmean_v1.l2 -0.01267212 -0.01840108  0.008579787 -0.01768321  0.01544655
#> tmean_v2.l1  0.12682050  0.11350967  0.110841083  0.15240943  0.13936834
#> tmean_v2.l2  0.17779782  0.14907805  0.255357655  0.19248097  0.29856105
#> tmean_v3.l1 -0.02411850 -0.02598190 -0.036443502 -0.03291361 -0.03715504
#>                sample361   sample362   sample363    sample364   sample365
#> (Intercept) -3.580255370 -1.71056528 -3.64187770 -2.627317204 -3.61440962
#> tmean_v1.l1  0.088796143  0.06002953  0.08591759  0.075840716  0.10252433
#> tmean_v1.l2 -0.005887962 -0.02362074 -0.02263977  0.004876997  0.03184393
#> tmean_v2.l1  0.164483161  0.07588448  0.16911185  0.137235815  0.16042799
#> tmean_v2.l2  0.248432004  0.10270542  0.19230870  0.212407008  0.34062406
#> tmean_v3.l1 -0.022240496 -0.03071613 -0.01721652 -0.014905028 -0.01643257
#>               sample366   sample367    sample368   sample369   sample370
#> (Intercept) -0.69743429 -4.39787330 -1.936935801 -4.49614916 -3.61013075
#> tmean_v1.l1  0.07122421  0.10866069  0.074162336  0.10290612  0.09477536
#> tmean_v1.l2 -0.04208700  0.03308922 -0.024074852  0.02351845  0.03322792
#> tmean_v2.l1  0.04103128  0.18140852  0.113438658  0.15834799  0.17350056
#> tmean_v2.l2  0.08807001  0.30898314  0.084038449  0.36271936  0.31067109
#> tmean_v3.l1 -0.04473615  0.00155141 -0.005756648 -0.03836273 -0.01373654
#>                sample371    sample372   sample373   sample374   sample375
#> (Intercept) -2.934992497 -3.066120818 -3.16474506 -0.42000731 -1.70509295
#> tmean_v1.l1  0.080195674  0.077949490  0.09064908  0.05670066  0.06409306
#> tmean_v1.l2  0.001906367 -0.013588450 -0.01173355 -0.06654990 -0.06786909
#> tmean_v2.l1  0.127613005  0.157974738  0.13018561  0.05625597  0.08981011
#> tmean_v2.l2  0.156385416  0.163700557  0.23552077 -0.03900863  0.02185047
#> tmean_v3.l1 -0.018052778 -0.007677538 -0.02747080 -0.02811995 -0.03533957
#>               sample376    sample377    sample378     sample379    sample380
#> (Intercept) -1.86065213 -3.638167966 -3.140439242 -3.7574762127 -4.855193366
#> tmean_v1.l1  0.07426732  0.094678900  0.094747597  0.0967518403  0.094529796
#> tmean_v1.l2 -0.02014341  0.002607119 -0.004196839 -0.0008951848  0.021093435
#> tmean_v2.l1  0.10189883  0.153637365  0.138389704  0.1614559805  0.218364604
#> tmean_v2.l2  0.08978192  0.274767436  0.210878872  0.2258141219  0.311517141
#> tmean_v3.l1 -0.01865599 -0.030211479 -0.022514023 -0.0231416151  0.004974149
#>                sample381   sample382    sample383   sample384    sample385
#> (Intercept) -3.743417530 -3.45946424 -4.325285657 -2.23154559 -1.964278870
#> tmean_v1.l1  0.087699896  0.10451263  0.095774860  0.07745541  0.069918076
#> tmean_v1.l2  0.018231076  0.01244764 -0.005302957 -0.04727693 -0.004673942
#> tmean_v2.l1  0.170278488  0.17004452  0.164348989  0.09295572  0.105844928
#> tmean_v2.l2  0.245813014  0.28485896  0.307742114  0.11472265  0.101905693
#> tmean_v3.l1 -0.007795989 -0.01916265 -0.045955279 -0.03990653 -0.009541263
#>               sample386    sample387   sample388   sample389   sample390
#> (Intercept) -1.17986478 -3.178146805 -2.76926252 -1.71564831 -1.82483705
#> tmean_v1.l1  0.07836188  0.084366801  0.08545090  0.06553652  0.08809968
#> tmean_v1.l2 -0.06390830 -0.002523818 -0.01037945 -0.04216427 -0.02559208
#> tmean_v2.l1  0.04153609  0.123272872  0.11886045  0.08521363  0.04910094
#> tmean_v2.l2  0.06799910  0.253474561  0.19780202  0.03341279  0.21567304
#> tmean_v3.l1 -0.06847963 -0.037298518 -0.02878592 -0.01700775 -0.08392190
#>               sample391    sample392    sample393   sample394   sample395
#> (Intercept) -2.96553287 -1.950879162 -2.711736413 -1.22088353 -1.46914787
#> tmean_v1.l1  0.08969020  0.074831261  0.089698112  0.07069642  0.06541154
#> tmean_v1.l2 -0.01703977 -0.006907343  0.005925808 -0.02051504 -0.02209971
#> tmean_v2.l1  0.12855404  0.085289821  0.143318069  0.08025313  0.07082885
#> tmean_v2.l2  0.18914811  0.148212635  0.198703691  0.09541051  0.10108901
#> tmean_v3.l1 -0.03512428 -0.025159944 -0.008844282 -0.01849955 -0.02161155
#>               sample396   sample397     sample398   sample399    sample400
#> (Intercept) -2.59975504 -0.53345157 -2.3906214299 -2.98224246 -1.448389573
#> tmean_v1.l1  0.08371012  0.07233665  0.0713671220  0.09227936  0.056490182
#> tmean_v1.l2 -0.02641783 -0.05532744 -0.0001618407  0.01633223 -0.002242084
#> tmean_v2.l1  0.10710774  0.04890896  0.1012388897  0.10844914  0.084882717
#> tmean_v2.l2  0.15936335  0.00229288  0.1611600343  0.27350826  0.042021977
#> tmean_v3.l1 -0.03752881 -0.03367866 -0.0176227216 -0.04740752  0.014242234
#>                 sample401    sample402    sample403   sample404   sample405
#> (Intercept) -2.3928069371 -3.351244127 -2.950077515 -3.02885052 -2.11343048
#> tmean_v1.l1  0.0754294922  0.082662490  0.092867731  0.08634033  0.07321683
#> tmean_v1.l2  0.0241563283 -0.002899583 -0.009412832 -0.02077161 -0.01976623
#> tmean_v2.l1  0.1243364092  0.150635649  0.101323381  0.14896334  0.09684978
#> tmean_v2.l2  0.1757771569  0.168708302  0.252577073  0.17663916  0.17901389
#> tmean_v3.l1  0.0007591115 -0.013578372 -0.069053427 -0.02614054 -0.03085229
#>                sample406    sample407   sample408   sample409   sample410
#> (Intercept) -2.120096082 -2.548111398 -1.82285926 -0.98472129 -2.50583158
#> tmean_v1.l1  0.074259249  0.082083573  0.07967524  0.06112923  0.08311022
#> tmean_v1.l2 -0.002191734 -0.001252609 -0.02372228 -0.03297368 -0.01924369
#> tmean_v2.l1  0.101824411  0.129816325  0.08336581  0.05917850  0.10489421
#> tmean_v2.l2  0.151723186  0.145202361  0.17982284  0.04101847  0.11791167
#> tmean_v3.l1 -0.017702387 -0.013516469 -0.04966940 -0.03869813 -0.03476710
#>                sample411   sample412    sample413    sample414    sample415
#> (Intercept) -3.265066544 -2.67303988 -2.902387480 -1.844506370 -0.771668191
#> tmean_v1.l1  0.092296735  0.08026697  0.081106054  0.077926362  0.077172218
#> tmean_v1.l2 -0.009348882 -0.02035415 -0.023482189 -0.005883208 -0.009724613
#> tmean_v2.l1  0.130317162  0.11849458  0.116817469  0.095536852  0.029682384
#> tmean_v2.l2  0.261653912  0.18664287  0.121365059  0.137510982  0.111665641
#> tmean_v3.l1 -0.031029247 -0.04122832 -0.008630375 -0.017392806 -0.045608390
#>                sample416   sample417    sample418   sample419    sample420
#> (Intercept) -1.839901356 -2.72270261 -2.649807162 -2.79094152 -1.555914980
#> tmean_v1.l1  0.065299656  0.09014596  0.093903651  0.08924832  0.064069208
#> tmean_v1.l2 -0.005858366 -0.01663926  0.004968839  0.01216986 -0.030840367
#> tmean_v2.l1  0.101177428  0.13546955  0.133909980  0.12347549  0.066896460
#> tmean_v2.l2  0.124930523  0.14648792  0.239886020  0.24053078  0.012334691
#> tmean_v3.l1 -0.010420899 -0.01365152 -0.019825330 -0.02619519 -0.009868708
#>                sample421    sample422    sample423   sample424   sample425
#> (Intercept) -2.966490422 -2.597540668 -1.766624563 -0.29022433 -3.83044869
#> tmean_v1.l1  0.080322052  0.088403645  0.079711134  0.06333363  0.10454111
#> tmean_v1.l2 -0.016881637 -0.002244444 -0.005359393 -0.02217430  0.03086961
#> tmean_v2.l1  0.166717375  0.117207224  0.063479598  0.03928557  0.14801091
#> tmean_v2.l2  0.118181040  0.168519235  0.176935468  0.01403979  0.33293810
#> tmean_v3.l1  0.009141218 -0.018007315 -0.033867963 -0.02182546 -0.04020565
#>               sample426   sample427   sample428   sample429   sample430
#> (Intercept) -4.03087176 -2.05100724 -3.18074433 -2.77538342 -0.52724521
#> tmean_v1.l1  0.10136995  0.07197171  0.10098945  0.08176796  0.07050145
#> tmean_v1.l2  0.01911252 -0.02408735  0.04590246 -0.01150859 -0.01515117
#> tmean_v2.l1  0.14061178  0.08445563  0.12050668  0.11061477  0.02311596
#> tmean_v2.l2  0.32974702  0.09560528  0.37428160  0.19660101  0.12814294
#> tmean_v3.l1 -0.04845277 -0.03458566 -0.04791035 -0.04212328 -0.05807549
#>               sample431   sample432   sample433   sample434   sample435
#> (Intercept) -3.27960526 -1.19541768 -0.43542248 -2.43806694 -3.84624331
#> tmean_v1.l1  0.07942137  0.07104798  0.04827266  0.08181634  0.08974041
#> tmean_v1.l2 -0.03559310 -0.04459505 -0.04620699 -0.01324831  0.01245365
#> tmean_v2.l1  0.14010967  0.04690394  0.02963958  0.12707356  0.17569616
#> tmean_v2.l2  0.17518794  0.05157817 -0.02291605  0.15726837  0.23490954
#> tmean_v3.l1 -0.02972220 -0.05115781 -0.02282589 -0.02111021 -0.01261359
#>                 sample436   sample437   sample438    sample439    sample440
#> (Intercept) -2.4492064436 -3.93680433 -2.73846846 -1.586775953 -2.997384082
#> tmean_v1.l1  0.0841030139  0.09798259  0.09024016  0.062716279  0.095372782
#> tmean_v1.l2  0.0001409512  0.03681374  0.01943656  0.004382563  0.007720372
#> tmean_v2.l1  0.0934987169  0.15627826  0.11811367  0.065268865  0.148020734
#> tmean_v2.l2  0.1639281292  0.31247428  0.25486740  0.047414872  0.169868138
#> tmean_v3.l1 -0.0262823413 -0.01991615 -0.01739453 -0.004690239 -0.014584790
#>               sample441   sample442   sample443   sample444   sample445
#> (Intercept) -4.08288512 -2.77603140 -2.52552943 -1.42750171 -1.30106639
#> tmean_v1.l1  0.09816005  0.08714511  0.08263899  0.06430323  0.07723791
#> tmean_v1.l2  0.03193309  0.01604250 -0.04528737 -0.05862232 -0.03910266
#> tmean_v2.l1  0.17370424  0.12148312  0.13960151  0.08335446  0.06209080
#> tmean_v2.l2  0.31168883  0.22151590  0.11744821  0.04917040  0.08878778
#> tmean_v3.l1 -0.01408201 -0.02320406 -0.01982209 -0.04128799 -0.05111518
#>               sample446   sample447   sample448   sample449   sample450
#> (Intercept) -3.26923206 -3.60910298 -4.81511582 -0.77909052 -4.49505645
#> tmean_v1.l1  0.08858178  0.10166664  0.10899319  0.06950446  0.09768405
#> tmean_v1.l2 -0.02815501  0.02990822  0.05820704 -0.05577179  0.04633351
#> tmean_v2.l1  0.16569988  0.13939528  0.19285242  0.05574930  0.19180995
#> tmean_v2.l2  0.21870711  0.33228955  0.40695091  0.02636460  0.31776675
#> tmean_v3.l1 -0.02611520 -0.03257345 -0.01393683 -0.03487178  0.01311912
#>                sample451    sample452   sample453    sample454    sample455
#> (Intercept) -4.157951079 -2.472104801 -0.45504425 -5.696729807 -3.978188374
#> tmean_v1.l1  0.106288734  0.080386032  0.05496847  0.107411757  0.103343999
#> tmean_v1.l2  0.008878278 -0.015958430 -0.04088257  0.079537765 -0.013405971
#> tmean_v2.l1  0.187145198  0.139242036  0.05281097  0.211040898  0.188810942
#> tmean_v2.l2  0.286278574  0.154275587 -0.03242359  0.466099406  0.219772921
#> tmean_v3.l1 -0.011377327 -0.008883086 -0.01658192  0.005207856 -0.007802275
#>               sample456    sample457    sample458   sample459   sample460
#> (Intercept) -5.87446702 -1.991138456 -4.637222717 -1.68795890 -1.17520067
#> tmean_v1.l1  0.12238210  0.083794224  0.107210142  0.07029465  0.06569101
#> tmean_v1.l2  0.06084208 -0.007360501  0.007808954  0.01326971 -0.04394120
#> tmean_v2.l1  0.21016057  0.095676989  0.192015394  0.05193240  0.05486705
#> tmean_v2.l2  0.49097690  0.216703530  0.296121586  0.20025776  0.05412409
#> tmean_v3.l1 -0.03065011 -0.033282762 -0.018664823 -0.02424690 -0.03364992
#>               sample461    sample462   sample463    sample464    sample465
#> (Intercept) -1.51030598 -2.937623768 -2.63736036 -3.702847852 -2.264199272
#> tmean_v1.l1  0.07638668  0.080143135  0.08384103  0.108497007  0.068191671
#> tmean_v1.l2 -0.01701648 -0.009047151 -0.01468816  0.008361003 -0.042811596
#> tmean_v2.l1  0.07300862  0.136709886  0.11100506  0.149567471  0.126508689
#> tmean_v2.l2  0.09062984  0.185162658  0.17206990  0.254143486  0.083206008
#> tmean_v3.l1 -0.01877564 -0.022309430 -0.02455797 -0.026874321 -0.006560031
#>                sample466   sample467    sample468    sample469    sample470
#> (Intercept) -0.849616592 -0.84686707 -3.606154134 -2.074031428 -2.236049480
#> tmean_v1.l1  0.061923257  0.06449558  0.096583615  0.069937764  0.077233336
#> tmean_v1.l2 -0.051300467 -0.03429628  0.026934884 -0.015001696 -0.009380261
#> tmean_v2.l1  0.061746086  0.06670929  0.170919187  0.109802230  0.093606573
#> tmean_v2.l2  0.003190418  0.03136374  0.287159029  0.085520778  0.173674141
#> tmean_v3.l1 -0.032031958 -0.01724485 -0.003088975 -0.009926726 -0.035086357
#>               sample471   sample472    sample473   sample474    sample475
#> (Intercept) -2.46785546 -2.02948832 -2.862544315 -0.94206584 -4.200792631
#> tmean_v1.l1  0.07590828  0.08480108  0.082174307  0.04874442  0.095042968
#> tmean_v1.l2 -0.01914750  0.01056016 -0.001688972 -0.05708007  0.001242216
#> tmean_v2.l1  0.09571600  0.07789534  0.120211565  0.05409701  0.204385801
#> tmean_v2.l2  0.15227598  0.26046201  0.183966660 -0.01682540  0.240813838
#> tmean_v3.l1 -0.03202771 -0.04792945 -0.019447775 -0.02553897  0.006558340
#>               sample476   sample477   sample478   sample479    sample480
#> (Intercept) -3.25927717 -1.99878914 -2.63446386 -1.55272654 -2.897938891
#> tmean_v1.l1  0.08501825  0.07125470  0.08002624  0.07061554  0.100248132
#> tmean_v1.l2 -0.02502984 -0.03913465 -0.01427547 -0.03390024  0.007373342
#> tmean_v2.l1  0.13727197  0.10748254  0.11374009  0.05887931  0.118122031
#> tmean_v2.l2  0.14484240  0.06795740  0.15946602  0.07490680  0.237596820
#> tmean_v3.l1 -0.02182202 -0.02446822 -0.01323589 -0.04745119 -0.037807889
#>               sample481    sample482   sample483    sample484   sample485
#> (Intercept) -3.41669123 -2.216205061 -1.79592054 -2.368587273 -3.45582925
#> tmean_v1.l1  0.09375848  0.087989332  0.07177818  0.073974798  0.09462034
#> tmean_v1.l2  0.02359334 -0.007975325 -0.04069174 -0.021918745 -0.01959690
#> tmean_v2.l1  0.13798176  0.109504561  0.09492869  0.116682028  0.14680997
#> tmean_v2.l2  0.35012716  0.182184824  0.02502611  0.102050502  0.20216450
#> tmean_v3.l1 -0.02460968 -0.030496926 -0.00843769 -0.008934704 -0.03263343
#>               sample486   sample487   sample488    sample489   sample490
#> (Intercept) -1.47514184 -2.65154699 -1.79868628 -2.392500654 -3.00978517
#> tmean_v1.l1  0.06925655  0.07634880  0.07597778  0.081076030  0.07982467
#> tmean_v1.l2 -0.02161647 -0.02448157 -0.01156812  0.009709737  0.01023971
#> tmean_v2.l1  0.07422744  0.10604735  0.11249116  0.123040237  0.13471958
#> tmean_v2.l2  0.08743741  0.12026119  0.10324224  0.172608455  0.20249685
#> tmean_v3.l1 -0.02769620 -0.02666855 -0.00323974 -0.011733826 -0.01024647
#>               sample491   sample492    sample493   sample494    sample495
#> (Intercept) -1.44965705 -4.49968460 -3.478679810 -3.44971207 -3.417642968
#> tmean_v1.l1  0.06614333  0.10269263  0.091750892  0.09842852  0.093708264
#> tmean_v1.l2 -0.01374390  0.02495063  0.006437786 -0.02963142  0.003008611
#> tmean_v2.l1  0.07444119  0.19117766  0.151182035  0.12738199  0.147024467
#> tmean_v2.l2  0.11311915  0.23049793  0.237427060  0.23876337  0.255798229
#> tmean_v3.l1 -0.01739438  0.01142118 -0.011989349 -0.05304810 -0.028156231
#>               sample496   sample497     sample498     sample499    sample500
#> (Intercept) -3.39783956 -1.14841188 -3.1367361062 -3.4791775063 -2.521189558
#> tmean_v1.l1  0.09059235  0.07173004  0.0835392125  0.0895467486  0.080060468
#> tmean_v1.l2  0.05497012 -0.04913322 -0.0004994823 -0.0001500896 -0.002648006
#> tmean_v2.l1  0.12518914  0.05176733  0.1264754798  0.1453925404  0.121529397
#> tmean_v2.l2  0.30426392  0.01913608  0.2244511687  0.2134765853  0.167847288
#> tmean_v3.l1 -0.01305946 -0.03717515 -0.0250926686 -0.0250213598 -0.003830402
#>               sample501   sample502   sample503   sample504   sample505
#> (Intercept) -2.43893215 -3.63226727 -1.44782525 -1.12734734 -2.27655249
#> tmean_v1.l1  0.08830585  0.09293165  0.06474655  0.06981341  0.06951753
#> tmean_v1.l2 -0.02944542 -0.01108674 -0.02992732 -0.03165973 -0.04264417
#> tmean_v2.l1  0.10146484  0.14919007  0.06958652  0.07424273  0.09057184
#> tmean_v2.l2  0.17476367  0.25561204  0.08951417  0.05986685  0.10911915
#> tmean_v3.l1 -0.04540951 -0.03620201 -0.03247231 -0.02501672 -0.04581180
#>               sample506   sample507   sample508   sample509   sample510
#> (Intercept) -1.10502421 -2.94407429 -1.90741675 -5.29249150 -1.36288070
#> tmean_v1.l1  0.07369738  0.08178225  0.07754302  0.11452227  0.07453562
#> tmean_v1.l2 -0.02882826 -0.02485185 -0.03806399  0.03466558 -0.04668286
#> tmean_v2.l1  0.06962054  0.13353904  0.07896279  0.19987398  0.06406581
#> tmean_v2.l2  0.10545233  0.11834091  0.14886592  0.35809025  0.05206551
#> tmean_v3.l1 -0.03989399 -0.01321576 -0.06481208 -0.02508930 -0.04229067
#>               sample511   sample512   sample513    sample514     sample515
#> (Intercept) -4.29807801 -0.05426331 -1.86461348 -2.923964939  0.5164922287
#> tmean_v1.l1  0.11070999  0.05468702  0.08157278  0.081199775  0.0467645540
#> tmean_v1.l2  0.04970747 -0.04873620  0.01494411 -0.007885379 -0.0446834001
#> tmean_v2.l1  0.16323549  0.03596450  0.10213313  0.115179918  0.0157018083
#> tmean_v2.l2  0.38776663  0.01879014  0.18451819  0.159758192 -0.0005838052
#> tmean_v3.l1 -0.02455373 -0.03333284 -0.01461894 -0.029193105 -0.0449248425
#>                sample516    sample517    sample518   sample519   sample520
#> (Intercept) -3.846615409 -2.547824877 -2.678530821 -2.27100333 -1.41194116
#> tmean_v1.l1  0.105186756  0.080821702  0.078527145  0.07839416  0.06562279
#> tmean_v1.l2  0.005345388  0.028210053  0.004683499 -0.02399126 -0.02625008
#> tmean_v2.l1  0.176519186  0.110701748  0.115518551  0.09690739  0.07037180
#> tmean_v2.l2  0.238333015  0.218293640  0.205721553  0.17517930  0.06517231
#> tmean_v3.l1 -0.012488327 -0.003713945 -0.014695578 -0.05280838 -0.03274559
#>               sample521   sample522   sample523   sample524    sample525
#> (Intercept) -1.41325835 -0.99803718 -3.40977110 -1.95716097 -2.322629717
#> tmean_v1.l1  0.06285760  0.07011691  0.09186742  0.07565567  0.087101761
#> tmean_v1.l2 -0.05160464 -0.01987664 -0.03785402 -0.02361322 -0.006103913
#> tmean_v2.l1  0.07117420  0.07988184  0.14262497  0.09053885  0.116585993
#> tmean_v2.l2  0.02605555  0.06035277  0.15539800  0.08895475  0.185898981
#> tmean_v3.l1 -0.03201144 -0.01894685 -0.03635705 -0.02176780 -0.024526729
#>               sample526    sample527   sample528    sample529    sample530
#> (Intercept) -3.12163669 -2.489093302 -1.93749258 -4.476828042 -0.886318669
#> tmean_v1.l1  0.09467896  0.071714536  0.06936947  0.107040000  0.061715974
#> tmean_v1.l2 -0.01509343  0.002242196 -0.01689963  0.027705165 -0.040023569
#> tmean_v2.l1  0.13948150  0.117056753  0.07178748  0.193635985  0.058510009
#> tmean_v2.l2  0.19807440  0.213090551  0.13190177  0.331506926 -0.005524225
#> tmean_v3.l1 -0.01875391 -0.007961932 -0.03438182 -0.003964928 -0.018878830
#>               sample531   sample532   sample533   sample534   sample535
#> (Intercept) -3.10869378 -1.33380105  0.75823443 -0.40791145 -4.24800171
#> tmean_v1.l1  0.08400543  0.06105519  0.03420551  0.05313841  0.10080974
#> tmean_v1.l2  0.00759246 -0.02760656 -0.09604009 -0.04469869  0.04667806
#> tmean_v2.l1  0.13344860  0.06646384 -0.01440527  0.04099284  0.14484375
#> tmean_v2.l2  0.22192550  0.07860226 -0.12817467 -0.03983088  0.37014591
#> tmean_v3.l1 -0.01146740 -0.03864309 -0.05353086 -0.02076475 -0.01625871
#>               sample536   sample537   sample538   sample539    sample540
#> (Intercept) -0.16332741 -2.67127019 -3.91100905 -1.08937776 -3.287032915
#> tmean_v1.l1  0.05305581  0.08382481  0.09393934  0.06226615  0.086714290
#> tmean_v1.l2 -0.01514204 -0.02620475  0.01814675 -0.05390314 -0.015855984
#> tmean_v2.l1  0.01778685  0.13533032  0.15874447  0.07086272  0.154937798
#> tmean_v2.l2  0.06006936  0.15984035  0.28493322  0.00219004  0.150628039
#> tmean_v3.l1 -0.02952113 -0.02579301 -0.02309806 -0.02197556 -0.004239041
#>                sample541    sample542    sample543   sample544   sample545
#> (Intercept) -2.130488416 -2.190933991 -2.941971073 -2.02329226 -1.59249623
#> tmean_v1.l1  0.081319848  0.076717557  0.098782365  0.07013412  0.07800229
#> tmean_v1.l2 -0.003592607  0.010963986 -0.008562443 -0.02717174 -0.03933916
#> tmean_v2.l1  0.120812903  0.107250640  0.154954009  0.11490669  0.08421352
#> tmean_v2.l2  0.162588013  0.161437529  0.238139390  0.09328848  0.10665295
#> tmean_v3.l1 -0.005746460 -0.003787333 -0.014038331 -0.02215559 -0.03530373
#>                sample546    sample547   sample548    sample549    sample550
#> (Intercept) -3.773389178 -2.434183410 -4.78587844 -2.590161206 -3.215218339
#> tmean_v1.l1  0.090268331  0.083078493  0.10130007  0.083185874  0.086840393
#> tmean_v1.l2  0.023370405  0.005101769  0.02040015 -0.007354994  0.026399376
#> tmean_v2.l1  0.177678824  0.130403397  0.19744315  0.103580729  0.132016006
#> tmean_v2.l2  0.219536072  0.174090192  0.27349924  0.178260092  0.234907727
#> tmean_v3.l1  0.007239605 -0.006276962 -0.00832226 -0.038226196 -0.006459075
#>               sample551   sample552    sample553    sample554   sample555
#> (Intercept) -2.10758720 -2.28538673 -1.283937607 -4.158885727 -2.52027810
#> tmean_v1.l1  0.08672999  0.09074260  0.059803473  0.093854330  0.08811207
#> tmean_v1.l2 -0.01002653 -0.02506309 -0.009319843  0.035626601  0.02644199
#> tmean_v2.l1  0.06305091  0.08871790  0.079901600  0.182823172  0.11748169
#> tmean_v2.l2  0.23674419  0.17650379  0.089935936  0.274883951  0.26195073
#> tmean_v3.l1 -0.06625470 -0.04888194 -0.007415525  0.003811255 -0.01753686
#>               sample556   sample557   sample558    sample559   sample560
#> (Intercept) -3.59079825 -2.38235426 -2.12870223 -4.625015649 -0.19374347
#> tmean_v1.l1  0.08573204  0.08842541  0.08267692  0.111400570  0.06705600
#> tmean_v1.l2  0.01792761  0.00606245 -0.01009010  0.052233581 -0.03174651
#> tmean_v2.l1  0.15478384  0.10262718  0.12015351  0.199765255  0.02793964
#> tmean_v2.l2  0.24395045  0.20503505  0.13959707  0.366145580  0.08657991
#> tmean_v3.l1 -0.02287878 -0.02810918 -0.01484996 -0.001531185 -0.05354528
#>               sample561    sample562    sample563    sample564   sample565
#> (Intercept) -3.75064514 -3.895794333 -2.762054017 -3.597766757 -1.10477997
#> tmean_v1.l1  0.09429959  0.111841276  0.081433613  0.098154261  0.06787938
#> tmean_v1.l2  0.03209254  0.009926575 -0.009430781  0.008289672 -0.04365984
#> tmean_v2.l1  0.12596840  0.172278624  0.129711083  0.141125656  0.05727365
#> tmean_v2.l2  0.27232755  0.350292726  0.170381816  0.251318296  0.06381674
#> tmean_v3.l1 -0.03099240 -0.043521962 -0.009680649 -0.016631182 -0.04865864
#>                sample566   sample567   sample568   sample569   sample570
#> (Intercept) -2.483467394 -1.94510321 -3.75380733 -2.05140553 -3.43837691
#> tmean_v1.l1  0.073399063  0.07617257  0.10252339  0.07802930  0.08193672
#> tmean_v1.l2 -0.009792371 -0.01588055  0.03490279 -0.01605827 -0.03576851
#> tmean_v2.l1  0.124500619  0.09864962  0.14627226  0.08574434  0.13423702
#> tmean_v2.l2  0.096509356  0.16074834  0.33318708  0.18926281  0.20366980
#> tmean_v3.l1 -0.015052276 -0.02661750 -0.03951195 -0.05385855 -0.05004782
#>               sample571   sample572   sample573   sample574   sample575
#> (Intercept) -2.27138019 -1.87936516 -3.48046514 -1.58286826 -3.70829058
#> tmean_v1.l1  0.08114003  0.07296593  0.10196192  0.06159014  0.10209238
#> tmean_v1.l2 -0.03159230 -0.02435094  0.02451897 -0.03320451  0.02906159
#> tmean_v2.l1  0.12474607  0.10730019  0.13298540  0.06157854  0.15812152
#> tmean_v2.l2  0.10780721  0.13190170  0.33722298  0.07897416  0.32715066
#> tmean_v3.l1 -0.01902687 -0.02962480 -0.03739812 -0.03318417 -0.01916012
#>                sample576    sample577    sample578    sample579   sample580
#> (Intercept) -2.384707536 -4.020415963 -2.529532106 -3.576924575 -1.87403046
#> tmean_v1.l1  0.085323101  0.087799138  0.086734805  0.088111174  0.07406476
#> tmean_v1.l2 -0.003053555 -0.001692815  0.002340296  0.032926922 -0.02823313
#> tmean_v2.l1  0.093743913  0.186443565  0.108960864  0.169010757  0.08422127
#> tmean_v2.l2  0.167769265  0.219168751  0.148128060  0.255513505  0.12648434
#> tmean_v3.l1 -0.031455198  0.004909984 -0.006899390  0.008995454 -0.03213688
#>                sample581   sample582   sample583   sample584   sample585
#> (Intercept) -1.010898384 -3.52081659 -3.27003433 -0.47474651 -2.86566200
#> tmean_v1.l1  0.057545570  0.09933842  0.08853375  0.06400889  0.08514525
#> tmean_v1.l2 -0.038206984  0.01378703 -0.01159168 -0.03770233  0.01998819
#> tmean_v2.l1  0.073634835  0.14710992  0.12065525  0.04927742  0.11209061
#> tmean_v2.l2  0.006226431  0.26373056  0.22590284  0.01586855  0.21844044
#> tmean_v3.l1 -0.023429954 -0.02300986 -0.04883826 -0.01746455 -0.01595046
#>                sample586   sample587   sample588    sample589   sample590
#> (Intercept) -3.401103992 -1.96473675 -3.16076774 -2.888123657 -2.04303662
#> tmean_v1.l1  0.082996681  0.07138806  0.08153172  0.072818038  0.06720712
#> tmean_v1.l2 -0.009745941 -0.04358790 -0.02202736 -0.010864662 -0.01087491
#> tmean_v2.l1  0.148954256  0.11084469  0.12539425  0.167280910  0.09689192
#> tmean_v2.l2  0.196527076  0.06521221  0.11563492  0.134344384  0.12242371
#> tmean_v3.l1 -0.008335330 -0.01410555 -0.01476761  0.008851952 -0.02191709
#>                sample591   sample592    sample593    sample594   sample595
#> (Intercept) -2.003044354 -1.93127895 -3.614478563 -3.593481204 -3.31633607
#> tmean_v1.l1  0.081040602  0.07133204  0.087049225  0.083166098  0.09205078
#> tmean_v1.l2  0.004844867 -0.01435561  0.031876089 -0.006782969 -0.01145682
#> tmean_v2.l1  0.123916052  0.08191754  0.142914019  0.171353124  0.14352958
#> tmean_v2.l2  0.188365697  0.16799472  0.260595252  0.189797880  0.25379095
#> tmean_v3.l1 -0.011283275 -0.03725874 -0.006556828 -0.001257343 -0.03259921
#>                sample596   sample597   sample598   sample599    sample600
#> (Intercept) -2.631851912 -1.04851954 -3.34352984 -1.46035666 -2.503497557
#> tmean_v1.l1  0.083364186  0.04984615  0.08545646  0.06790834  0.083304594
#> tmean_v1.l2  0.009812454 -0.03818784 -0.03291050 -0.02338683 -0.009678781
#> tmean_v2.l1  0.107635641  0.06008103  0.14951544  0.09580250  0.120918802
#> tmean_v2.l2  0.182820165 -0.04238838  0.14755697  0.06421567  0.199837100
#> tmean_v3.l1 -0.024944521 -0.01186677 -0.02013845 -0.02274027 -0.024497541
#>                sample601   sample602   sample603    sample604     sample605
#> (Intercept) -3.166873325 -2.73032926 -3.12740576 -3.501140181 -2.7787243705
#> tmean_v1.l1  0.095198958  0.09233146  0.07873495  0.087942369  0.0815480898
#> tmean_v1.l2  0.023121586  0.01744382 -0.01044050  0.013827990 -0.0008766962
#> tmean_v2.l1  0.141162847  0.12529607  0.13912688  0.165879459  0.1220441431
#> tmean_v2.l2  0.247361413  0.21757337  0.17714154  0.218080088  0.1980131218
#> tmean_v3.l1 -0.007889898 -0.01307076 -0.02225568  0.003502735 -0.0156536667
#>               sample606   sample607    sample608   sample609   sample610
#> (Intercept) -2.93588405 -3.51607729 -2.990909695 -2.57549470 -1.53450294
#> tmean_v1.l1  0.10019218  0.09290791  0.082923096  0.08819631  0.05498399
#> tmean_v1.l2 -0.01591827  0.01530120 -0.007080901 -0.02246666 -0.04241147
#> tmean_v2.l1  0.14009707  0.12543270  0.135474583  0.12374799  0.08493833
#> tmean_v2.l2  0.24593935  0.27610220  0.163191494  0.14059801  0.01215654
#> tmean_v3.l1 -0.04591700 -0.03297702 -0.018061471 -0.02149413 -0.01322257
#>               sample611    sample612   sample613    sample614   sample615
#> (Intercept) -1.79256505 -4.617855450 -2.31223006 -2.100438883 -1.11051222
#> tmean_v1.l1  0.08955080  0.103618254  0.08899650  0.076431870  0.06946915
#> tmean_v1.l2 -0.05467085  0.036424128 -0.03384352  0.006805387 -0.02773127
#> tmean_v2.l1  0.10359925  0.208436979  0.13434854  0.098895362  0.05439538
#> tmean_v2.l2  0.11195994  0.293804539  0.11592910  0.168904112  0.03676440
#> tmean_v3.l1 -0.04348182 -0.001601736 -0.02169060 -0.015026509 -0.03257851
#>                sample616   sample617    sample618   sample619   sample620
#> (Intercept) -2.851046778 -3.10986877 -3.550956620 -2.51260695 -2.94886526
#> tmean_v1.l1  0.094575568  0.08444271  0.092934610  0.08723677  0.08774564
#> tmean_v1.l2  0.002691819  0.02280181  0.004904828 -0.01786028  0.01402539
#> tmean_v2.l1  0.119013006  0.11134173  0.167887230  0.11779235  0.11334640
#> tmean_v2.l2  0.251981960  0.24441541  0.286835647  0.16304118  0.26024296
#> tmean_v3.l1 -0.037880161 -0.01853610 -0.023381426 -0.03512044 -0.02421619
#>               sample621    sample622   sample623   sample624    sample625
#> (Intercept) -1.17891479 -3.161398506 -3.27773983 -3.72547947 -5.090945287
#> tmean_v1.l1  0.06657813  0.085602661  0.08362636  0.07745052  0.099044737
#> tmean_v1.l2 -0.08094599  0.004640483  0.01795582 -0.01818111  0.010885934
#> tmean_v2.l1  0.04649693  0.118332294  0.14001884  0.14697171  0.214960271
#> tmean_v2.l2  0.02120033  0.219516925  0.26202292  0.17582045  0.291483173
#> tmean_v3.l1 -0.06116326 -0.016767181 -0.01633427 -0.02977887  0.002290751
#>                sample626   sample627    sample628   sample629   sample630
#> (Intercept) -2.225773848 -2.41670792 -2.260043974 -2.55261645 -4.04088701
#> tmean_v1.l1  0.077660234  0.08042227  0.083959123  0.07768121  0.09192784
#> tmean_v1.l2  0.014640194  0.01499541 -0.023198589 -0.01935792  0.00639910
#> tmean_v2.l1  0.114806732  0.09188312  0.128431899  0.10401651  0.17696682
#> tmean_v2.l2  0.153347937  0.21274248  0.095622662  0.15670977  0.26170891
#> tmean_v3.l1  0.002355555 -0.02332004 -0.005643929 -0.03217643 -0.01175502
#>               sample631   sample632    sample633   sample634   sample635
#> (Intercept) -3.36301215 -2.33141378 -3.104559298 -1.53188839 -1.30549224
#> tmean_v1.l1  0.08822595  0.08561809  0.092042133  0.06297045  0.07622174
#> tmean_v1.l2  0.01302240 -0.02638911 -0.002136461 -0.01936336 -0.01485768
#> tmean_v2.l1  0.13688685  0.09772403  0.156519929  0.04905070  0.06983014
#> tmean_v2.l2  0.20188988  0.11642822  0.195631565  0.10875596  0.14468986
#> tmean_v3.l1 -0.01546407 -0.03166860 -0.012863472 -0.04519299 -0.06240012
#>               sample636   sample637    sample638   sample639   sample640
#> (Intercept) -3.17718348 -2.54864438 -2.743043429 -0.84370545 -1.48524884
#> tmean_v1.l1  0.08998120  0.07145385  0.088340530  0.05968945  0.07357803
#> tmean_v1.l2 -0.04321608 -0.01984168  0.008135092 -0.05386396 -0.01925639
#> tmean_v2.l1  0.13859091  0.12409408  0.141888038  0.04432773  0.08417561
#> tmean_v2.l2  0.17594112  0.13156911  0.242855909  0.03262799  0.13930371
#> tmean_v3.l1 -0.04429806 -0.00833352 -0.014270441 -0.05068324 -0.03289532
#>                sample641   sample642   sample643    sample644    sample645
#> (Intercept) -4.145821832 -2.80537521 -3.94382390 -3.800674003 -3.145491416
#> tmean_v1.l1  0.105338503  0.06927669  0.09984936  0.095789243  0.086173176
#> tmean_v1.l2  0.036860549 -0.05367231  0.03825984  0.001848811  0.025966591
#> tmean_v2.l1  0.173568228  0.11999828  0.15476249  0.149060926  0.126426628
#> tmean_v2.l2  0.318822391  0.07048834  0.34765643  0.272524473  0.203432372
#> tmean_v3.l1 -0.007495497 -0.01877624 -0.01776810 -0.027483249 -0.007403326
#>               sample646    sample647   sample648   sample649    sample650
#> (Intercept) -1.06308739 -3.459331910 -1.45621502 -2.61711827 -3.328107423
#> tmean_v1.l1  0.06397750  0.089417451  0.07145971  0.07086349  0.077779603
#> tmean_v1.l2 -0.05866419  0.007489478 -0.01444342 -0.04400636 -0.007334954
#> tmean_v2.l1  0.08967157  0.158919487  0.06820326  0.12079666  0.143186237
#> tmean_v2.l2 -0.01385355  0.255482243  0.09133069  0.09056011  0.147163714
#> tmean_v3.l1 -0.01774925 -0.029769798 -0.02051443 -0.01950027  0.005272944
#>               sample651   sample652    sample653     sample654    sample655
#> (Intercept) -2.49758016 -4.20207143 -2.768167925 -3.4838861140  0.084942953
#> tmean_v1.l1  0.08596264  0.10638995  0.085271463  0.0979814976  0.056910129
#> tmean_v1.l2  0.01682261  0.03436124  0.009912503  0.0193787570 -0.053654089
#> tmean_v2.l1  0.10787041  0.17652976  0.125844963  0.1521901338  0.016848518
#> tmean_v2.l2  0.21909508  0.35610650  0.276321936  0.2450814386  0.009441559
#> tmean_v3.l1 -0.01381673 -0.02133349 -0.022411699 -0.0001308189 -0.044534640
#>                sample656    sample657    sample658   sample659   sample660
#> (Intercept) -4.353673792 -4.030654506 -2.681576508 -1.57712161 -2.12668809
#> tmean_v1.l1  0.088255769  0.100067235  0.079056316  0.07171999  0.07085796
#> tmean_v1.l2 -0.001196185  0.006392858  0.026798689 -0.04664867 -0.02699024
#> tmean_v2.l1  0.180086588  0.187724290  0.129401982  0.06279270  0.10722422
#> tmean_v2.l2  0.219524936  0.252248781  0.196816751  0.08535877  0.05558082
#> tmean_v3.l1 -0.009164691 -0.003679254  0.007309216 -0.04354225 -0.01798158
#>               sample661    sample662   sample663    sample664    sample665
#> (Intercept) -3.21881657 -2.285314157 -3.39253571 -2.427381323 -4.003483867
#> tmean_v1.l1  0.08428616  0.078468498  0.09863450  0.074561210  0.081441348
#> tmean_v1.l2 -0.01841264 -0.015432067  0.03120214 -0.001713057 -0.015314155
#> tmean_v2.l1  0.14594428  0.121815040  0.13929837  0.100182626  0.182277843
#> tmean_v2.l2  0.15892310  0.093251978  0.28330948  0.196136254  0.212124246
#> tmean_v3.l1 -0.02905791  0.005016867 -0.02060248 -0.035996374  0.003280331
#>               sample666    sample667   sample668     sample669    sample670
#> (Intercept) -4.36192584 -2.011728493 -4.10874452 -4.3549937356 -4.172658775
#> tmean_v1.l1  0.11067196  0.064150305  0.08892035  0.0920748791  0.098525802
#> tmean_v1.l2  0.06124925 -0.008931930  0.01115223  0.0003765713 -0.001143387
#> tmean_v2.l1  0.15076668  0.096512135  0.18044454  0.1813438537  0.176132924
#> tmean_v2.l2  0.42636709  0.101906925  0.25058185  0.2359661982  0.285315350
#> tmean_v3.l1 -0.01321178 -0.003337033  0.00882608 -0.0049831898 -0.032272061
#>               sample671    sample672    sample673   sample674    sample675
#> (Intercept)  0.30598084 -4.264377951 -3.465656726 -2.41607168 -1.991482773
#> tmean_v1.l1  0.04894128  0.105792673  0.097282640  0.08065717  0.073271239
#> tmean_v1.l2 -0.08022932  0.037598073  0.008421587 -0.01337603 -0.005224331
#> tmean_v2.l1 -0.01145920  0.183470133  0.153705991  0.09211515  0.102326911
#> tmean_v2.l2 -0.05655525  0.304953505  0.280049630  0.19784605  0.100884621
#> tmean_v3.l1 -0.07003544  0.001814678 -0.033403257 -0.03910512 -0.014547642
#>               sample676    sample677   sample678   sample679   sample680
#> (Intercept) -3.05691517 -1.283839544 -2.04485502 -2.98212756 -2.00078885
#> tmean_v1.l1  0.07271274  0.066106486  0.07103620  0.08684149  0.08084400
#> tmean_v1.l2 -0.01785561 -0.006167534 -0.02803739 -0.00864355 -0.03450761
#> tmean_v2.l1  0.14053227  0.048306680  0.10361527  0.13646545  0.09663781
#> tmean_v2.l2  0.06215759  0.155445947  0.11519860  0.19449068  0.09539394
#> tmean_v3.l1  0.01791798 -0.047686458 -0.02635109 -0.02285039 -0.02830025
#>               sample681   sample682    sample683   sample684   sample685
#> (Intercept) -2.16446868 -0.77836925 -0.853258487 -3.85907913 -3.03981099
#> tmean_v1.l1  0.08319321  0.05524262  0.058324934  0.09392454  0.08154584
#> tmean_v1.l2 -0.03765995 -0.05304857 -0.039849103  0.03070053 -0.03120733
#> tmean_v2.l1  0.10451713  0.04288047  0.054303331  0.18987484  0.14460127
#> tmean_v2.l2  0.12063310 -0.01647490 -0.008165018  0.25692291  0.09497892
#> tmean_v3.l1 -0.04437152 -0.04543415 -0.024034726  0.02315815 -0.01031869
#>               sample686   sample687   sample688    sample689   sample690
#> (Intercept) -3.75867874 -3.26987505 -4.50157379 -3.166528214 -2.03837610
#> tmean_v1.l1  0.09570330  0.09161991  0.10758127  0.079274621  0.07453450
#> tmean_v1.l2  0.04203439 -0.02019662  0.04534376  0.009282668 -0.05113442
#> tmean_v2.l1  0.15900466  0.15326863  0.19676454  0.148835236  0.09935223
#> tmean_v2.l2  0.33754002  0.19035861  0.29195982  0.206446593  0.08825847
#> tmean_v3.l1 -0.01900403 -0.02178874  0.02060072 -0.013697559 -0.04023283
#>                sample691    sample692    sample693    sample694    sample695
#> (Intercept) -3.075577515 -2.988881291 -2.110132659 -3.476757520 -2.755066795
#> tmean_v1.l1  0.080462202  0.088366156  0.081394086  0.087591502  0.082252486
#> tmean_v1.l2 -0.004319920 -0.014575546 -0.005854796  0.005281485  0.025694939
#> tmean_v2.l1  0.140001281  0.141258911  0.108360453  0.154814759  0.126937092
#> tmean_v2.l2  0.158697968  0.142445355  0.178128742  0.180519728  0.229463773
#> tmean_v3.l1 -0.005512799 -0.003201552 -0.024558149 -0.013570088 -0.005576534
#>               sample696   sample697    sample698   sample699   sample700
#> (Intercept) -1.61071198 -0.84538645 -1.225006788 -3.50566083 -2.99116835
#> tmean_v1.l1  0.08038246  0.06451695  0.056089200  0.09285222  0.08302450
#> tmean_v1.l2 -0.03397353 -0.04774324 -0.059580898 -0.01485406 -0.03006942
#> tmean_v2.l1  0.06879184  0.05110146  0.089810774  0.12497912  0.14485782
#> tmean_v2.l2  0.13581980 -0.01164957 -0.037747406  0.18179732  0.13804240
#> tmean_v3.l1 -0.05929261 -0.02728920  0.008901788 -0.03536154 -0.01369640
#>               sample701   sample702   sample703    sample704   sample705
#> (Intercept) -1.40049081 -2.08910650 -3.33431846 -5.423440999 -1.30800946
#> tmean_v1.l1  0.07342374  0.08405186  0.08640925  0.108621326  0.06790648
#> tmean_v1.l2 -0.03953890 -0.02614816  0.01841403  0.042730081 -0.03074867
#> tmean_v2.l1  0.07805067  0.11798370  0.14621731  0.223672747  0.07905866
#> tmean_v2.l2  0.05071649  0.17338256  0.22458944  0.402832135  0.06059468
#> tmean_v3.l1 -0.03752810 -0.03230270 -0.02671435 -0.005574607 -0.03481893
#>               sample706    sample707   sample708    sample709   sample710
#> (Intercept) -2.37869890 -3.371708515 -1.55513245 -2.743215687 -1.04444256
#> tmean_v1.l1  0.07105137  0.088690603  0.07260617  0.088473374  0.07208029
#> tmean_v1.l2 -0.02766416  0.008922813 -0.01901987  0.003626158 -0.01854162
#> tmean_v2.l1  0.11901886  0.148815792  0.05982952  0.084552807  0.05405760
#> tmean_v2.l2  0.06847558  0.250410838  0.10880311  0.239299064  0.11230064
#> tmean_v3.l1  0.00345207 -0.017640884 -0.02543610 -0.042822936 -0.03892659
#>                sample711   sample712   sample713   sample714    sample715
#> (Intercept) -1.903182123 -5.29635638 -2.98632353 -5.14824902 -4.047883674
#> tmean_v1.l1  0.074208424  0.11390270  0.09091166  0.11072200  0.100625910
#> tmean_v1.l2  0.004073494  0.04765638  0.01506519  0.01743517  0.029354176
#> tmean_v2.l1  0.056105569  0.19723549  0.13546823  0.24840172  0.189043068
#> tmean_v2.l2  0.160844492  0.42187042  0.26614528  0.28388651  0.312286597
#> tmean_v3.l1 -0.041923819 -0.02889586 -0.02047951  0.01787677  0.006878819
#>               sample716   sample717   sample718    sample719     sample720
#> (Intercept) -3.75575178 -4.29264109 -2.39231618 -2.145068821 -2.7580865900
#> tmean_v1.l1  0.08569147  0.09468917  0.08506265  0.081360306  0.0817172261
#> tmean_v1.l2 -0.01595530  0.00524170  0.01524904 -0.007872861  0.0025334623
#> tmean_v2.l1  0.16657262  0.20926047  0.10429695  0.138948463  0.1530132923
#> tmean_v2.l2  0.20307915  0.18110337  0.21801967  0.110649847  0.1671294474
#> tmean_v3.l1 -0.01447907  0.01631713 -0.01799341 -0.009346181  0.0001484011
#>               sample721   sample722    sample723    sample724   sample725
#> (Intercept) -2.02585692 -3.78998437 -2.187992439 -3.790572824 -0.77695032
#> tmean_v1.l1  0.07448253  0.10401385  0.080659843  0.093427951  0.05429711
#> tmean_v1.l2 -0.01950030  0.01911011  0.002109383  0.001101831 -0.02712182
#> tmean_v2.l1  0.07151615  0.14315712  0.113005364  0.168280893  0.04795396
#> tmean_v2.l2  0.15059872  0.27862965  0.144847151  0.249370592  0.01987857
#> tmean_v3.l1 -0.04771661 -0.01493765 -0.014050288 -0.016121674 -0.01040495
#>               sample726    sample727    sample728   sample729    sample730
#> (Intercept) -4.45915019 -3.382552276 -2.448278443 -1.10547263 -3.688196614
#> tmean_v1.l1  0.10126737  0.088550767  0.067986560  0.05526845  0.095821272
#> tmean_v1.l2  0.02771397  0.024953785 -0.006388993 -0.03705750  0.016712547
#> tmean_v2.l1  0.16900824  0.155426870  0.111380493  0.05466132  0.166831250
#> tmean_v2.l2  0.31043063  0.217887352  0.120774925  0.06416566  0.272614351
#> tmean_v3.l1 -0.01679556  0.003031814 -0.010040834 -0.03612725 -0.002554989
#>               sample731    sample732    sample733   sample734   sample735
#> (Intercept) -1.25545786 -3.243494698 -3.663237364 -1.96290345 -1.69116082
#> tmean_v1.l1  0.05809485  0.094357831  0.095847821  0.07092908  0.06753063
#> tmean_v1.l2 -0.01255379 -0.006469271  0.007864019 -0.02275191  0.01001971
#> tmean_v2.l1  0.05096682  0.150207970  0.144530142  0.09101549  0.08549761
#> tmean_v2.l2  0.11122715  0.195255176  0.196124280  0.13331018  0.15064637
#> tmean_v3.l1 -0.03066587 -0.016219335  0.002677943 -0.02015480 -0.01450848
#>               sample736   sample737   sample738   sample739     sample740
#> (Intercept) -1.88653713 -1.87752412 -2.04653861 -0.82657210 -1.8292199665
#> tmean_v1.l1  0.06680131  0.07252039  0.07395100  0.07039719  0.0766372463
#> tmean_v1.l2 -0.03191182 -0.02833482 -0.02512142 -0.03780183  0.0008300401
#> tmean_v2.l1  0.08239426  0.07515159  0.08429847  0.02830710  0.0629542141
#> tmean_v2.l2  0.08592632  0.09696030  0.11841633  0.03692575  0.1847769898
#> tmean_v3.l1 -0.02949589 -0.04400877 -0.03637235 -0.04826146 -0.0488712406
#>               sample741   sample742   sample743    sample744    sample745
#> (Intercept)  0.01039347 -1.83195189 -3.20356022 -3.331543269 -3.575086682
#> tmean_v1.l1  0.05379150  0.06221688  0.09455280  0.080419478  0.095406949
#> tmean_v1.l2 -0.05462118 -0.05354938  0.02592862 -0.001779737  0.004981971
#> tmean_v2.l1  0.02880228  0.09264125  0.11125588  0.158654783  0.159847306
#> tmean_v2.l2 -0.02018028  0.05215882  0.28180261  0.210035238  0.230762159
#> tmean_v3.l1 -0.04409006 -0.02974270 -0.03976367 -0.009514096 -0.011986377
#>                sample746   sample747   sample748    sample749    sample750
#> (Intercept) -2.057057872 -2.66767624 -2.69258205 -2.701804745 -4.787259961
#> tmean_v1.l1  0.073206587  0.08288387  0.07531899  0.080492323  0.104546034
#> tmean_v1.l2 -0.009300549 -0.02707108 -0.01873367  0.002584534 -0.000516584
#> tmean_v2.l1  0.081203798  0.09911524  0.16466540  0.133988367  0.171571281
#> tmean_v2.l2  0.169421426  0.16302878  0.08078348  0.169232672  0.279397656
#> tmean_v3.l1 -0.039941984 -0.04790964  0.02266591 -0.010498463 -0.036033355
#>                sample751    sample752   sample753   sample754   sample755
#> (Intercept) -4.162776727 -0.630303658 -0.39937608 -2.04806267 -2.98312800
#> tmean_v1.l1  0.101822793  0.056424264  0.05974855  0.07224426  0.09672102
#> tmean_v1.l2  0.023529759 -0.055416887 -0.06289514 -0.02008978  0.02107097
#> tmean_v2.l1  0.184084345  0.047399407  0.02504782  0.11036998  0.12408775
#> tmean_v2.l2  0.286198121 -0.002666641 -0.01792280  0.11975714  0.25599689
#> tmean_v3.l1 -0.002625791 -0.031516424 -0.05638197 -0.02815253 -0.01810390
#>               sample756   sample757   sample758   sample759   sample760
#> (Intercept) -1.71931003 -3.33584521 -2.32955187 -3.32987618 -0.35054255
#> tmean_v1.l1  0.06799537  0.08776864  0.08433887  0.10220513  0.05441105
#> tmean_v1.l2 -0.03535576  0.01449593 -0.02073245 -0.01047975 -0.04934808
#> tmean_v2.l1  0.09072397  0.14486695  0.11890444  0.10865812  0.02088806
#> tmean_v2.l2  0.05651076  0.18078877  0.19140979  0.29146829 -0.02255652
#> tmean_v3.l1 -0.02832487 -0.01338788 -0.04453629 -0.06477159 -0.04195576
#>               sample761   sample762    sample763   sample764   sample765
#> (Intercept) -3.44664087 -1.48089275 -2.815141880 -4.12366734 -1.80888426
#> tmean_v1.l1  0.08973659  0.06812140  0.091336566  0.10630530  0.07856941
#> tmean_v1.l2  0.02148083 -0.02098982  0.006532624  0.01642418 -0.00507734
#> tmean_v2.l1  0.14258093  0.07422253  0.110731780  0.17710241  0.08375368
#> tmean_v2.l2  0.29692729  0.12257408  0.245328960  0.34375273  0.19183532
#> tmean_v3.l1 -0.01436751 -0.03601564 -0.034850179 -0.02632209 -0.03304015
#>               sample766   sample767   sample768   sample769   sample770
#> (Intercept) -2.65907828 -4.38188008 -1.83694303 -3.97084051 -0.29298284
#> tmean_v1.l1  0.07084567  0.10450107  0.07419918  0.10482788  0.05200240
#> tmean_v1.l2 -0.02040570  0.03463697 -0.04174769  0.01272138 -0.07901007
#> tmean_v2.l1  0.10988841  0.18313368  0.11490203  0.18361270  0.01441302
#> tmean_v2.l2  0.14206345  0.27375511  0.06883055  0.21777009 -0.09484729
#> tmean_v3.l1 -0.01973543  0.00207037 -0.01725944  0.01011835 -0.04848739
#>               sample771    sample772   sample773    sample774    sample775
#> (Intercept) -5.63222672 -3.410289298 -1.29304701 -2.579998377 -2.148854508
#> tmean_v1.l1  0.11546261  0.101610342  0.06943054  0.089509372  0.076815575
#> tmean_v1.l2  0.02327496  0.005953105 -0.01511015  0.006050778 -0.008810432
#> tmean_v2.l1  0.21618205  0.154245098  0.06944951  0.111869269  0.078386150
#> tmean_v2.l2  0.37781216  0.243802825  0.08520214  0.226334933  0.223305769
#> tmean_v3.l1 -0.02440433 -0.012721938 -0.04127983 -0.030861293 -0.046718411
#>                sample776   sample777   sample778   sample779   sample780
#> (Intercept) -3.592703651 -2.18237985 -0.34830039 -1.43091618 -1.85321194
#> tmean_v1.l1  0.090193894  0.07861734  0.07390058  0.07336512  0.06749886
#> tmean_v1.l2  0.008867864 -0.04034437 -0.03872870  0.01303568 -0.03871010
#> tmean_v2.l1  0.126017661  0.12616445  0.05476327  0.06863994  0.11476573
#> tmean_v2.l2  0.308142939  0.06617997  0.07710433  0.15199115  0.10326298
#> tmean_v3.l1 -0.048385880 -0.01056593 -0.03286589 -0.03060537 -0.02774333
#>                sample781    sample782   sample783   sample784    sample785
#> (Intercept) -2.413335922 -3.813150655 -2.11774993 -2.85276732 -3.214335039
#> tmean_v1.l1  0.070155482  0.094681838  0.08078631  0.08609874  0.085529388
#> tmean_v1.l2 -0.018649848  0.008228711 -0.01541301 -0.03993836  0.026114715
#> tmean_v2.l1  0.126414963  0.161274109  0.08793282  0.15888512  0.128050821
#> tmean_v2.l2  0.080282263  0.287501628  0.16362723  0.12750247  0.182736843
#> tmean_v3.l1 -0.002561343 -0.013902315 -0.03392027 -0.01766972  0.002627862
#>                sample786   sample787    sample788    sample789   sample790
#> (Intercept) -2.648499245 -1.09259366 -3.301262953 -2.883214408 -3.29143442
#> tmean_v1.l1  0.078856120  0.06226874  0.096137062  0.085699147  0.10226562
#> tmean_v1.l2 -0.009069593 -0.01613951 -0.007373738  0.007241624  0.03396458
#> tmean_v2.l1  0.165228244  0.07823668  0.179516276  0.151264725  0.14495437
#> tmean_v2.l2  0.138249840  0.06504600  0.213754683  0.202086170  0.28950734
#> tmean_v3.l1  0.003834364 -0.01880703 -0.002332273  0.005896027 -0.01209929
#>                sample791   sample792   sample793   sample794    sample795
#> (Intercept) -1.891998345 -2.03052444 -3.23159393 -2.35537930 -3.868084002
#> tmean_v1.l1  0.076663182  0.06878678  0.09447423  0.06094646  0.098231248
#> tmean_v1.l2 -0.002873052 -0.04207407  0.02254210 -0.02954966  0.038810374
#> tmean_v2.l1  0.096682699  0.08321086  0.11988863  0.08953413  0.164744889
#> tmean_v2.l2  0.152495337  0.09132536  0.25439458  0.12231508  0.312082660
#> tmean_v3.l1 -0.018918197 -0.04501551 -0.02904719 -0.03634293  0.006868744
#>               sample796    sample797    sample798     sample799   sample800
#> (Intercept) -2.92991692 -1.276234596 -2.722287165 -2.8079560927 -1.43016713
#> tmean_v1.l1  0.07537166  0.074266908  0.078583302  0.0969869747  0.05761518
#> tmean_v1.l2 -0.01597423 -0.009993188  0.004879882  0.0010263381 -0.02895880
#> tmean_v2.l1  0.16193453  0.076160816  0.094703546  0.1286636643  0.08414615
#> tmean_v2.l2  0.10823663  0.093941247  0.164086767  0.1737124263  0.06224252
#> tmean_v3.l1  0.00274459 -0.006743161 -0.018151905  0.0003705624 -0.02336492
#>                sample801    sample802   sample803    sample804    sample805
#> (Intercept) -3.701011045 -2.684256790 -2.44457820 -2.139734393 -2.773780493
#> tmean_v1.l1  0.095549786  0.097214136  0.08384265  0.077880688  0.083637756
#> tmean_v1.l2  0.025968817 -0.008077327 -0.01392264  0.009838002  0.025207423
#> tmean_v2.l1  0.165261765  0.104182948  0.10523940  0.091311175  0.134800864
#> tmean_v2.l2  0.261381209  0.188602188  0.15273902  0.235116834  0.232140568
#> tmean_v3.l1 -0.006069706 -0.030059311 -0.02282670 -0.024976868 -0.006184249
#>                sample806    sample807    sample808   sample809    sample810
#> (Intercept) -3.508155444 -2.846458492 -4.330528003 -3.11926599 -2.617253837
#> tmean_v1.l1  0.079823031  0.083821297  0.103256765  0.09292190  0.088633416
#> tmean_v1.l2  0.017047827 -0.007999355  0.016824687 -0.01606902 -0.001254002
#> tmean_v2.l1  0.172196805  0.132282451  0.188282008  0.12712958  0.155339933
#> tmean_v2.l2  0.202588912  0.216166532  0.241055230  0.23422108  0.187150871
#> tmean_v3.l1  0.008165795 -0.030283815 -0.001924623 -0.04201302  0.001891688
#>               sample811    sample812   sample813   sample814   sample815
#> (Intercept) -3.98755335 -3.148693756 -0.34086935 -1.31752086 -0.19832064
#> tmean_v1.l1  0.10179530  0.103039602  0.05412169  0.06616737  0.06235160
#> tmean_v1.l2  0.01734337  0.003661358 -0.03332192 -0.02133598 -0.05370416
#> tmean_v2.l1  0.15120667  0.160170814  0.03065650  0.04720707  0.02269554
#> tmean_v2.l2  0.34290719  0.274828043  0.03076094  0.08847729  0.03175712
#> tmean_v3.l1 -0.04857517 -0.022259599 -0.03821035 -0.04350942 -0.05849261
#>               sample816   sample817   sample818     sample819     sample820
#> (Intercept) -1.86304955 -3.36726434 -4.40338501  0.1015948368 -3.4104688916
#> tmean_v1.l1  0.08901092  0.09358454  0.09900157  0.0642236774  0.0931355156
#> tmean_v1.l2 -0.02216450  0.02047460  0.04102762 -0.0427034489  0.0002161071
#> tmean_v2.l1  0.07701173  0.15001928  0.16305651  0.0009528968  0.1460009347
#> tmean_v2.l2  0.16263429  0.27503230  0.35570140  0.0522723295  0.2425238593
#> tmean_v3.l1 -0.05876081 -0.01487987 -0.01424979 -0.0640760335 -0.0168317628
#>               sample821   sample822    sample823    sample824   sample825
#> (Intercept) -1.11031616 -2.09007325 -2.317394936 -2.858904699 -2.14727826
#> tmean_v1.l1  0.06331571  0.08271692  0.072347775  0.085919462  0.07153863
#> tmean_v1.l2 -0.05977601  0.02402381 -0.026670882 -0.003437751 -0.05949603
#> tmean_v2.l1  0.05558658  0.10207023  0.117688050  0.134228218  0.08144371
#> tmean_v2.l2  0.01441627  0.18062480  0.061926195  0.194198683  0.16065766
#> tmean_v3.l1 -0.04909231 -0.00114294 -0.009213194 -0.026055970 -0.06815673
#>               sample826   sample827    sample828    sample829   sample830
#> (Intercept) -2.70760010 -1.75911589 -2.133657066 -3.008300792 -2.26929571
#> tmean_v1.l1  0.08797984  0.06838772  0.079918506  0.086980004  0.08077331
#> tmean_v1.l2  0.00928045 -0.02660557 -0.001669047 -0.006143071 -0.04032671
#> tmean_v2.l1  0.12939611  0.08133660  0.090258526  0.127355724  0.11406980
#> tmean_v2.l2  0.17455627  0.07297403  0.167617762  0.189631687  0.15024550
#> tmean_v3.l1 -0.00361538 -0.02734835 -0.018555757 -0.021934267 -0.04723819
#>               sample831     sample832    sample833   sample834   sample835
#> (Intercept) -3.85388712 -0.8619717961 -2.648958208 -2.64740996 -1.53783701
#> tmean_v1.l1  0.09329027  0.0614580819  0.081239970  0.08339059  0.08113743
#> tmean_v1.l2 -0.02338388 -0.0678822758 -0.003307844 -0.01582332 -0.02735331
#> tmean_v2.l1  0.16692910  0.0504574118  0.115662332  0.10243886  0.07231809
#> tmean_v2.l2  0.19723400  0.0002205431  0.118370972  0.16614467  0.15997373
#> tmean_v3.l1 -0.02651345 -0.0483495540 -0.003663434 -0.02549406 -0.04858312
#>               sample836   sample837    sample838     sample839   sample840
#> (Intercept) -3.52622971 -0.26708234 -1.147190693 -1.5220264751 -4.45209691
#> tmean_v1.l1  0.09968092  0.05670428  0.073636497  0.0761478917  0.11654810
#> tmean_v1.l2  0.02864228 -0.01648460 -0.006525104 -0.0002689094  0.06241838
#> tmean_v2.l1  0.15080330  0.01781862  0.048137184  0.0853668530  0.16954583
#> tmean_v2.l2  0.32406419  0.03982477  0.126382286  0.1636561569  0.45022641
#> tmean_v3.l1 -0.02785254 -0.03577944 -0.039461908 -0.0237362426 -0.02156924
#>               sample841    sample842    sample843   sample844   sample845
#> (Intercept) -4.45772447 -3.623706938 -3.526314200 -1.95830834 -2.67668631
#> tmean_v1.l1  0.10070228  0.099180371  0.108813874  0.06826272  0.08164538
#> tmean_v1.l2  0.03599434  0.001971726  0.007637362 -0.05045776 -0.01128141
#> tmean_v2.l1  0.18415758  0.161837006  0.179619627  0.09673235  0.09872695
#> tmean_v2.l2  0.32321870  0.243594769  0.291503156  0.10176877  0.12048645
#> tmean_v3.l1 -0.01017842 -0.015443643 -0.019440166 -0.05403465 -0.02124251
#>                 sample846    sample847   sample848   sample849   sample850
#> (Intercept) -3.5229966813  0.041153813 -4.18738905  0.03819708 -4.76156247
#> tmean_v1.l1  0.0867220439  0.056997190  0.09555106  0.05837147  0.09498787
#> tmean_v1.l2 -0.0003409854 -0.063696163 -0.01482460 -0.05348622 -0.01516661
#> tmean_v2.l1  0.1559121420  0.034150066  0.15736509  0.03611706  0.21336533
#> tmean_v2.l2  0.2342786601 -0.006728809  0.27093222  0.03107502  0.20638614
#> tmean_v3.l1 -0.0034176704 -0.041615052 -0.04613277 -0.04370037 -0.00807923
#>               sample851   sample852    sample853    sample854   sample855
#> (Intercept) -2.78706060 -2.17299284 -2.839463590 -0.840263058 -2.25970524
#> tmean_v1.l1  0.09245219  0.07718859  0.087536584  0.066985424  0.07386416
#> tmean_v1.l2  0.01674220 -0.01713527  0.003806463 -0.046573746 -0.01126046
#> tmean_v2.l1  0.12256735  0.08732680  0.144212879  0.052078591  0.12530290
#> tmean_v2.l2  0.24076400  0.18832347  0.203131200  0.008998352  0.12239539
#> tmean_v3.l1 -0.02614747 -0.04237130 -0.018725525 -0.023904649 -0.01534855
#>                 sample856   sample857   sample858    sample859    sample860
#> (Intercept) -2.4296827679 -4.20758097 -2.30734493 -4.616948893 -1.518321194
#> tmean_v1.l1  0.0924195973  0.11393143  0.08398579  0.098949700  0.080596403
#> tmean_v1.l2  0.0002810488  0.04959737 -0.03019188  0.009357802 -0.003187192
#> tmean_v2.l1  0.0790039594  0.15786593  0.08899226  0.182322639  0.063073509
#> tmean_v2.l2  0.2573343215  0.33516081  0.11987288  0.307712728  0.178158187
#> tmean_v3.l1 -0.0688438557 -0.01856525 -0.02886357 -0.017580046 -0.053816520
#>               sample861   sample862   sample863   sample864   sample865
#> (Intercept)  0.81843048 -3.71659410 -1.36739516 -1.99647800 -4.24048293
#> tmean_v1.l1  0.03455907  0.09211049  0.07049221  0.08064827  0.10075578
#> tmean_v1.l2 -0.07009564 -0.01484225 -0.03497026 -0.02526625  0.02029400
#> tmean_v2.l1 -0.01608244  0.15878727  0.09941483  0.07536644  0.18791892
#> tmean_v2.l2 -0.11321990  0.19029251  0.06845576  0.16430903  0.33023812
#> tmean_v3.l1 -0.03847067 -0.02352060 -0.02014185 -0.05182819 -0.01436954
#>                sample866    sample867    sample868   sample869   sample870
#> (Intercept) -3.207051616 -3.040788835 -3.384206824 -1.64457107 -1.03002387
#> tmean_v1.l1  0.091234488  0.078538569  0.093914378  0.06674485  0.05886534
#> tmean_v1.l2  0.017303561  0.005908163  0.013664776 -0.02095998 -0.02019862
#> tmean_v2.l1  0.133056427  0.153246910  0.147512002  0.09405331  0.04907977
#> tmean_v2.l2  0.213068423  0.161589936  0.234303650  0.12013604  0.06441514
#> tmean_v3.l1 -0.008740589  0.008460698 -0.002931327 -0.02181957 -0.02391756
#>               sample871    sample872   sample873   sample874    sample875
#> (Intercept) -5.62011509 -1.591911279 -2.82824383 -3.60042957 -5.969433290
#> tmean_v1.l1  0.10928531  0.076236130  0.09126199  0.08999985  0.121182545
#> tmean_v1.l2  0.01045976 -0.001232745  0.01378248 -0.02118801  0.008952369
#> tmean_v2.l1  0.22513553  0.109285967  0.11714016  0.14011152  0.231210535
#> tmean_v2.l2  0.30999267  0.119456885  0.24108750  0.21845461  0.381396437
#> tmean_v3.l1 -0.01168073  0.001710053 -0.02716497 -0.04252997 -0.025593606
#>               sample876    sample877   sample878    sample879    sample880
#> (Intercept) -2.15194839 -3.569871460 -3.07432191 -5.432051390 -3.392201945
#> tmean_v1.l1  0.08420836  0.094215658  0.07709394  0.106315015  0.084399100
#> tmean_v1.l2 -0.01603636  0.005854559 -0.01090719  0.038171078  0.005087667
#> tmean_v2.l1  0.09800610  0.135010833  0.14280456  0.231930489  0.148044489
#> tmean_v2.l2  0.15501247  0.258916871  0.19610616  0.363507671  0.224823511
#> tmean_v3.l1 -0.05333821 -0.006574921 -0.01548303  0.009701459 -0.008338663
#>               sample881    sample882    sample883   sample884   sample885
#> (Intercept) -0.64708542 -2.923442545 -4.544564338 -4.06768437 -2.67976290
#> tmean_v1.l1  0.05473615  0.074311526  0.092876512  0.09526560  0.08402246
#> tmean_v1.l2 -0.06212055 -0.003791503 -0.018544983  0.01836863 -0.03443601
#> tmean_v2.l1  0.06767144  0.134583901  0.191885466  0.16411193  0.10391612
#> tmean_v2.l2 -0.02804104  0.188484162  0.211684138  0.27122129  0.14818796
#> tmean_v3.l1 -0.02695459 -0.014826104 -0.005161146 -0.01614720 -0.04693997
#>                sample886   sample887   sample888    sample889   sample890
#> (Intercept) -1.103760974 -1.69843819 -2.26223298 -0.842967493 -4.71743185
#> tmean_v1.l1  0.064812237  0.06945282  0.07562687  0.058639824  0.10364414
#> tmean_v1.l2 -0.056768200 -0.03223029 -0.03358750 -0.042775223 -0.00326363
#> tmean_v2.l1  0.066912416  0.06607221  0.12093786  0.041026908  0.17777756
#> tmean_v2.l2  0.009006114  0.10676476  0.10728011 -0.007215848  0.28973506
#> tmean_v3.l1 -0.027358641 -0.04607052 -0.02902886 -0.026594993 -0.04289937
#>               sample891   sample892   sample893   sample894     sample895
#> (Intercept) -0.51246165 -2.93162378 -4.10816266 -5.43233370 -3.9510437305
#> tmean_v1.l1  0.05485289  0.08091390  0.08856322  0.11409818  0.0938924241
#> tmean_v1.l2 -0.05139046 -0.02001154 -0.01121304  0.04474175  0.0002338085
#> tmean_v2.l1  0.04547956  0.13250834  0.16646277  0.21786225  0.1568291847
#> tmean_v2.l2 -0.01903042  0.14858987  0.24353650  0.38394108  0.2496176934
#> tmean_v3.l1 -0.02816316 -0.02145319 -0.02725186 -0.01420523 -0.0290553778
#>               sample896    sample897   sample898   sample899   sample900
#> (Intercept) -2.80507450 -1.589225794 -1.24598984 -1.43625622 -1.42293703
#> tmean_v1.l1  0.07268908  0.078830004  0.05975283  0.05749604  0.06491785
#> tmean_v1.l2 -0.01169329 -0.006151718 -0.03685288 -0.02578881 -0.03395315
#> tmean_v2.l1  0.13254708  0.093020879  0.07572295  0.07019071  0.06905798
#> tmean_v2.l2  0.17086227  0.163817322  0.00999694  0.10579561  0.09998181
#> tmean_v3.l1 -0.01372304 -0.039434083 -0.01822234 -0.03240384 -0.03443004
#>                sample901   sample902   sample903   sample904    sample905
#> (Intercept) -4.522246761 -1.99639094 -1.96917020 -0.25634643 -4.038473648
#> tmean_v1.l1  0.085884082  0.07090232  0.07142992  0.04461629  0.092599073
#> tmean_v1.l2 -0.002862393 -0.01645074 -0.03625357 -0.03014609  0.004221393
#> tmean_v2.l1  0.209750318  0.07270335  0.08677495  0.03532732  0.164351876
#> tmean_v2.l2  0.227667686  0.11249092  0.08070558  0.03123350  0.244690404
#> tmean_v3.l1 -0.004937294 -0.02562278 -0.03088672 -0.01370445 -0.012022681
#>                sample906    sample907    sample908   sample909   sample910
#> (Intercept) -3.208993973  0.848145950 -4.170652964 -1.23735077 -2.68876172
#> tmean_v1.l1  0.074172529  0.046428702  0.088276121  0.07415003  0.08862978
#> tmean_v1.l2 -0.006948267 -0.052719069  0.010571164 -0.03660890 -0.01981629
#> tmean_v2.l1  0.122025791 -0.008063137  0.202013493  0.06166058  0.12403345
#> tmean_v2.l2  0.145652723 -0.062502121  0.229329908  0.03318256  0.10297634
#> tmean_v3.l1 -0.031771871 -0.037349359  0.008912137 -0.02557072 -0.02166286
#>                sample911     sample912    sample913     sample914   sample915
#> (Intercept) -2.066440985 -2.8878339253 -2.665714569 -3.579215e+00 -1.85336077
#> tmean_v1.l1  0.079597953  0.0825398533  0.088277978  8.338605e-02  0.07744404
#> tmean_v1.l2  0.002817126  0.0003004605 -0.005039972  7.173735e-06 -0.01497673
#> tmean_v2.l1  0.120695111  0.1461318337  0.098437944  1.683993e-01  0.08073145
#> tmean_v2.l2  0.148880806  0.2223610620  0.189099139  1.703888e-01  0.14854912
#> tmean_v3.l1 -0.002471400 -0.0209949670 -0.030217240 -1.141782e-03 -0.03539714
#>               sample916   sample917   sample918   sample919    sample920
#> (Intercept) -2.02438803 -0.63644442 -4.35225458 -0.92617161 -4.122190269
#> tmean_v1.l1  0.07585788  0.05540923  0.10322644  0.05545311  0.090461453
#> tmean_v1.l2 -0.02022577 -0.04900028  0.04476277 -0.05908793 -0.006112375
#> tmean_v2.l1  0.09238595  0.02909943  0.17146901  0.04405005  0.184639205
#> tmean_v2.l2  0.19430601  0.03161258  0.35910291  0.05289080  0.193391473
#> tmean_v3.l1 -0.05372330 -0.05407919 -0.01758310 -0.05167603  0.003350028
#>               sample921   sample922   sample923   sample924   sample925
#> (Intercept) -1.26956805 -0.30301435 -2.01868450 -2.34600170 -2.76717963
#> tmean_v1.l1  0.06403960  0.05463835  0.07344352  0.08815119  0.07955555
#> tmean_v1.l2 -0.05619986 -0.06413038 -0.01127141 -0.01144526 -0.01710843
#> tmean_v2.l1  0.05497626  0.04067825  0.08746035  0.11710996  0.13786578
#> tmean_v2.l2  0.04906364 -0.05072079  0.18347611  0.16951810  0.17564551
#> tmean_v3.l1 -0.05709263 -0.01757357 -0.04126652 -0.02186534 -0.02201236
#>                sample926    sample927   sample928    sample929   sample930
#> (Intercept) -1.745119964 -2.022482952 -4.22482741 -2.816400840 -2.35199659
#> tmean_v1.l1  0.076334550  0.065024135  0.10178953  0.076794525  0.07629459
#> tmean_v1.l2 -0.004241409 -0.007987856 -0.01801637 -0.029581884  0.01054677
#> tmean_v2.l1  0.092932238  0.092897587  0.21238862  0.124016992  0.09758659
#> tmean_v2.l2  0.153339914  0.149714337  0.23095702  0.091206234  0.15916048
#> tmean_v3.l1 -0.012920233 -0.019026235 -0.01493416 -0.007654989 -0.01052044
#>               sample931   sample932    sample933   sample934   sample935
#> (Intercept) -4.10939143  0.75764038 -3.470341653 -2.05038094 -2.97947078
#> tmean_v1.l1  0.10058664  0.02973334  0.103304741  0.07891300  0.06936566
#> tmean_v1.l2  0.03349020 -0.06301673  0.008229092 -0.02387661 -0.01549463
#> tmean_v2.l1  0.14907175 -0.00875112  0.151828184  0.09175411  0.13876933
#> tmean_v2.l2  0.30132223 -0.11577547  0.274744211  0.16482865  0.09429521
#> tmean_v3.l1 -0.02952696 -0.02331125 -0.028476833 -0.04429825  0.01024263
#>               sample936    sample937   sample938   sample939   sample940
#> (Intercept) -2.29818508 -3.569679738 -2.66445870 -0.85612491 -2.26954132
#> tmean_v1.l1  0.08667734  0.092743927  0.08302217  0.06570689  0.08256223
#> tmean_v1.l2 -0.01111021 -0.005806569 -0.01214733 -0.06731918  0.00194664
#> tmean_v2.l1  0.11013298  0.153865473  0.11825088  0.03489119  0.07146217
#> tmean_v2.l2  0.14999322  0.244285139  0.18272847  0.03367522  0.14756613
#> tmean_v3.l1 -0.03980515 -0.024960458 -0.01767791 -0.06336391 -0.02531072
#>               sample941   sample942     sample943   sample944    sample945
#> (Intercept) -1.60809295 -3.33074312 -4.3885831864 -2.47832235 -0.143748130
#> tmean_v1.l1  0.07597274  0.08590094  0.0984576192  0.09223853  0.058689446
#> tmean_v1.l2 -0.02063263  0.02652179  0.0589885389 -0.01120055 -0.044606574
#> tmean_v2.l1  0.07737942  0.13391786  0.1897027728  0.09861551  0.009895964
#> tmean_v2.l2  0.15355387  0.25575245  0.3322090964  0.16615913 -0.009103748
#> tmean_v3.l1 -0.04215231 -0.01880563 -0.0003823895 -0.01470619 -0.043086991
#>               sample946    sample947  sample948   sample949   sample950
#> (Intercept) -3.51391786 -3.428575695 -3.7114678 -1.22712210 -2.80529126
#> tmean_v1.l1  0.08120370  0.087438869  0.1009601  0.06704018  0.07898778
#> tmean_v1.l2  0.01581507  0.047435911  0.0207053 -0.03414107 -0.02380161
#> tmean_v2.l1  0.15823429  0.153371832  0.1617619  0.07813540  0.10965661
#> tmean_v2.l2  0.11957565  0.277581778  0.2579130  0.04064615  0.15968709
#> tmean_v3.l1  0.01241471  0.009344604 -0.0179726 -0.02894725 -0.03023693
#>               sample951   sample952    sample953   sample954   sample955
#> (Intercept) -3.64179464 -4.54375673 -1.449753000 -2.43636629 -2.36714531
#> tmean_v1.l1  0.09586340  0.10570818  0.067568578  0.08364885  0.07981485
#> tmean_v1.l2  0.01331646  0.02761598 -0.037294767 -0.01916338 -0.01230941
#> tmean_v2.l1  0.14621116  0.20957454  0.100260078  0.12270246  0.10376201
#> tmean_v2.l2  0.27035657  0.33126691  0.049385619  0.12028895  0.15097755
#> tmean_v3.l1 -0.03025951 -0.01015259 -0.007092159 -0.01212972 -0.02023273
#>               sample956    sample957   sample958   sample959   sample960
#> (Intercept) -1.93507033 -3.711263777 -1.60416215 -3.08199373 -2.93955091
#> tmean_v1.l1  0.06725416  0.091096496  0.07966499  0.08618756  0.07984698
#> tmean_v1.l2 -0.03152334  0.016887545 -0.01573324 -0.03877811 -0.02037556
#> tmean_v2.l1  0.08561235  0.173497682  0.06377440  0.12194962  0.11015046
#> tmean_v2.l2  0.09333690  0.213769422  0.16904415  0.12271531  0.16860563
#> tmean_v3.l1 -0.03003503 -0.007790386 -0.05399711 -0.02712181 -0.04007825
#>               sample961   sample962   sample963   sample964    sample965
#> (Intercept) -1.08539021 -0.33271339 -2.52678054 -3.35245515 -2.069553695
#> tmean_v1.l1  0.07656260  0.06108666  0.07538388  0.09822344  0.076787797
#> tmean_v1.l2 -0.01825729 -0.04387114 -0.00172161  0.02307317  0.006145296
#> tmean_v2.l1  0.06809874  0.02190702  0.11312754  0.12980565  0.085105774
#> tmean_v2.l2  0.08685580 -0.01068416  0.17619562  0.29075504  0.187156827
#> tmean_v3.l1 -0.03876298 -0.02883705 -0.02964190 -0.03088410 -0.042097374
#>                sample966    sample967    sample968    sample969    sample970
#> (Intercept) -1.469670517  0.017152596 -0.774321677 -3.838548806 -3.334074139
#> tmean_v1.l1  0.055756020  0.068055313  0.065302979  0.091571509  0.081060523
#> tmean_v1.l2 -0.038168123 -0.064134553 -0.030361464 -0.006867957 -0.002506829
#> tmean_v2.l1  0.088786962  0.042243129  0.043531844  0.188699730  0.123609696
#> tmean_v2.l2 -0.007927579  0.008228037  0.008908907  0.198922820  0.207090069
#> tmean_v3.l1  0.007485815 -0.044132151 -0.029396610 -0.005552755 -0.033204296
#>                sample971   sample972   sample973   sample974   sample975
#> (Intercept) -2.772444439 -2.67144981 -1.49420652 -1.69864633 -3.86495851
#> tmean_v1.l1  0.079479459  0.09243760  0.07431226  0.06816487  0.10082279
#> tmean_v1.l2 -0.009169975 -0.02034834 -0.04004454 -0.01613373  0.02305540
#> tmean_v2.l1  0.125707240  0.12547911  0.07310322  0.05903795  0.15573412
#> tmean_v2.l2  0.139017540  0.17136597  0.13271336  0.14546415  0.33818913
#> tmean_v3.l1 -0.009779151 -0.02908781 -0.05191245 -0.04556908 -0.03672088
#>               sample976    sample977   sample978    sample979     sample980
#> (Intercept) -1.86828894 -0.531830403 -1.31764280 -3.020702170 -3.068429e+00
#> tmean_v1.l1  0.06704996  0.058457835  0.07076395  0.084245998  9.835593e-02
#> tmean_v1.l2 -0.05000182 -0.070160696 -0.05540107  0.004816869  1.340305e-05
#> tmean_v2.l1  0.09416952  0.052224482  0.06508186  0.133792844  1.733174e-01
#> tmean_v2.l2  0.05249417 -0.009651846  0.03042652  0.203743078  2.394502e-01
#> tmean_v3.l1 -0.03808720 -0.028587990 -0.04344607 -0.005039408  3.836748e-03
#>               sample981     sample982    sample983   sample984   sample985
#> (Intercept) -0.80546317 -3.5753266196 -2.004578654 -5.44900384 -3.43407293
#> tmean_v1.l1  0.08401210  0.0902106775  0.059904665  0.10558758  0.09130452
#> tmean_v1.l2 -0.02338995 -0.0002038959 -0.015521830  0.03415448  0.01951846
#> tmean_v2.l1  0.06418486  0.1625865107  0.095557904  0.22173060  0.10835804
#> tmean_v2.l2  0.10403164  0.2005712559  0.047572333  0.39085439  0.27112182
#> tmean_v3.l1 -0.03281003 -0.0059216723  0.005435755 -0.01068518 -0.03681456
#>                 sample986   sample987    sample988   sample989   sample990
#> (Intercept) -3.2939363124 -5.00150643 -2.824780987 -0.27186088 -2.03422160
#> tmean_v1.l1  0.0849532854  0.08280324  0.084815980  0.04289821  0.07395735
#> tmean_v1.l2  0.0007456022  0.01911763 -0.014533977 -0.06658685 -0.02982008
#> tmean_v2.l1  0.1316488758  0.21568600  0.157127841  0.01569090  0.08739473
#> tmean_v2.l2  0.2242009126  0.27046076  0.150225063 -0.05810534  0.11246188
#> tmean_v3.l1 -0.0392704447  0.01048075 -0.004148644 -0.03994517 -0.02625547
#>                sample991    sample992   sample993   sample994    sample995
#> (Intercept) -2.408556445 -2.605612426 -3.92041043 -2.89824446 -0.502014377
#> tmean_v1.l1  0.082840820  0.075178763  0.11400314  0.08557832  0.056296861
#> tmean_v1.l2 -0.007292709 -0.019831546  0.00811155 -0.01880723 -0.049531762
#> tmean_v2.l1  0.107761380  0.129429435  0.16515602  0.13781575  0.034867734
#> tmean_v2.l2  0.184716923  0.110656841  0.32334845  0.17076278 -0.007137611
#> tmean_v3.l1 -0.034554618 -0.009890857 -0.04386683 -0.02350569 -0.033762826
#>               sample996   sample997    sample998   sample999   sample1000
#> (Intercept) -0.34171588 -1.98405242 -3.050901946 -3.49539138 -3.663866058
#> tmean_v1.l1  0.05144153  0.07529620  0.086636776  0.08378304  0.087400330
#> tmean_v1.l2 -0.04065108 -0.01135425 -0.008792713 -0.03575424  0.006764306
#> tmean_v2.l1  0.03263626  0.05896222  0.135238920  0.14448274  0.157132817
#> tmean_v2.l2 -0.01511849  0.17407197  0.239778917  0.18511009  0.224957814
#> tmean_v3.l1 -0.02557461 -0.06804312 -0.028812107 -0.05122032 -0.007320438
```

Posterior coefficient summaries:

``` r


fit_bdlnm$coefficients.summary
#>                       mean          sd  0.025quant     0.5quant   0.975quant
#> (Intercept)   -2.612308666 1.181011217 -4.86173293 -2.659073564 -0.271741415
#> tmean_v1.l1    0.081827689 0.014285124  0.05468580  0.081841641  0.108500115
#> tmean_v1.l2   -0.008881804 0.026473706 -0.05979620 -0.009364572  0.044742272
#> tmean_v2.l1    0.116775335 0.046154672  0.02792660  0.118181920  0.209262892
#> tmean_v2.l2    0.170686691 0.097874377 -0.01795049  0.173623565  0.358115565
#> tmean_v3.l1   -0.023144634 0.017116497 -0.05412196 -0.022838542  0.008978175
#> tmean_v3.l2    0.119642910 0.038671644  0.04575193  0.120948587  0.197361979
#> rain_v1.l1     0.128096901 0.026269755  0.07581237  0.129003281  0.178838059
#> rain_v1.l2    -0.239686348 0.033790575 -0.30519545 -0.239825436 -0.169897806
#> rain_v2.l1     0.056297062 0.016763660  0.02603588  0.055803307  0.090203720
#> rain_v2.l2    -0.058826311 0.024053725 -0.10389328 -0.059182823 -0.012161076
#> rain_v3.l1    -0.038708686 0.038858203 -0.11258207 -0.038844164  0.038305274
#> rain_v3.l2     0.058778348 0.053433229 -0.04509105  0.056894884  0.164445332
#> wetness_v1.l1  0.052609959 0.003723103  0.04538249  0.052604149  0.059713221
#> wetness_v1.l2 -0.093896682 0.011069018 -0.11651535 -0.093612535 -0.073576562
#> wetness_v2.l1  0.053719862 0.014373251  0.02434727  0.053718821  0.082410374
#> wetness_v2.l2 -0.093244212 0.038430918 -0.16970272 -0.093340193 -0.017618942
#> wetness_v3.l1  0.055673487 0.009673492  0.03697027  0.055688038  0.073729222
#> wetness_v3.l2 -0.090573738 0.018771943 -0.12753240 -0.090908305 -0.056623194
#>                      mode
#> (Intercept)   -2.76660074
#> tmean_v1.l1    0.08230571
#> tmean_v1.l2   -0.01506542
#> tmean_v2.l1    0.12123077
#> tmean_v2.l2    0.17477959
#> tmean_v3.l1   -0.02269479
#> tmean_v3.l2    0.12843803
#> rain_v1.l1     0.13628989
#> rain_v1.l2    -0.24148249
#> rain_v2.l1     0.05473600
#> rain_v2.l2    -0.05831061
#> rain_v3.l1    -0.03565000
#> rain_v3.l2     0.03976396
#> wetness_v1.l1  0.05214164
#> wetness_v1.l2 -0.09609658
#> wetness_v2.l1  0.05363214
#> wetness_v2.l2 -0.10424887
#> wetness_v3.l1  0.05508213
#> wetness_v3.l2 -0.09485818
```

#### Model objects

``` r

class(fit_inla)
#> [1] "inla"

class(fit_bdlnm)
#> [1] "bdlnm"

class(fit_bdlnm$model)
#> [1] "inla"
```

The direct INLA engine is useful when users want a conventional INLA
representation of the epidemic-level design matrix and direct access to
standard INLA outputs. The B-DLNM engine is useful when users want to
retain the native DLNM representation and work directly with posterior
coefficient samples, fitted basis objects, and Bayesian DLNM workflows.
