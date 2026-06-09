# Understanding exposure-lag-response relationships

## Objective

`Visualizar a superfície completa`

``` r
surf <- predict_surface(fit)

plot(surf)
```

``` r
plot(surf, type="contour", main="Exposure–lag–response surface")
```

👉 Aqui você explica:

resposta depende de DOIS fatores:

intensidade da exposição momento (lag)

👉 conceito importante:

O efeito não é instantâneo — ele é distribuído ao longo do tempo.
