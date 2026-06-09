# Identifying critical exposure periods in epidemics

## Objective

Descobrir períodos críticos

``` r
crit <- identify_critical_lags(fit)
```

🧠 Discussão 👉 aqui está uma DAS PARTES MAIS IMPORTANTES Explique:

Nem todos os lags são igualmente importantes.

👉 insight epidemiológico:

janelas críticas períodos de susceptibilidade

``` r
barplot(crit$importance,
        names.arg = crit$lag,
        col = "tomato",
        main = "Critical lag importance")
```
