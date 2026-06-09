# From statistical output to epidemiological interpretation

## Objective

Reduzir dimensionalidade (interpretabilidade)

``` r
red <- reduce_effects(fit)

plot(red)
```

👉 aqui você explica:

DLNM = alta dimensão → difícil interpretar redução permite visualizar:

Efeito acumulado ao longo do tempo

``` r
plot(red$type, red$effect,
     type="l", lwd=2,
     main="Reduced exposure-response relationship")
```
