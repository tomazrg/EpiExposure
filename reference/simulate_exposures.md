# Simulate chronological exposure profiles

Generates one or more complete exposure profiles in chronological order,
from the earliest/oldest observation to the most recent observation.

## Usage

``` r
simulate_exposures(
  max_lag,
  n = 1,
  mode = c("profile", "pattern"),
  pattern = c("random_times", "alternating_times", "block_times"),
  fixed_value = NULL,
  n_times = NULL,
  time_range = NULL,
  alternating_values = NULL,
  block_times = NULL,
  background = list(dist = "normal", mean = 0, sd = 1),
  bounds = NULL,
  seed = NULL,
  cumulative = FALSE
)
```

## Arguments

- max_lag:

  Non-negative integer maximum retrospective lag. Every generated
  profile has exactly \`max_lag + 1\` chronological positions.

- n:

  Positive integer number of exposure profiles to simulate. Default is
  1.

- mode:

  Character. \`"profile"\` generates the background process only;
  \`"pattern"\` generates the background and then applies one explicit
  chronological pattern.

- pattern:

  Character pattern used when \`mode = "pattern"\`:

  \- \`"random_times"\`: assign \`fixed_value\` at \`n_times\` distinct
  chronological positions sampled without replacement from
  \`time_range\`; - \`"alternating_times"\`: recycle
  \`alternating_values\` through \`time_range\` in chronological
  order; - \`"block_times"\`: assign \`fixed_value\` to the complete
  inclusive interval defined by \`block_times\`.

- fixed_value:

  Finite numeric scalar used by \`"random_times"\` and
  \`"block_times"\`.

- n_times:

  Positive integer number of selected chronological positions for
  \`"random_times"\`.

- time_range:

  Optional unique integer vector of chronological times used by
  \`"random_times"\` or \`"alternating_times"\`. Values must be between
  0 and \`max_lag\`, inclusive. If \`NULL\`, every chronological
  position is eligible.

  \`time_range\` is not used by \`"block_times"\`; supplying it with
  that pattern is an error rather than being silently ignored.

- alternating_values:

  Numeric vector containing at least two finite values. Values are
  recycled across the selected \`time_range\` in chronological order
  when \`pattern = "alternating_times"\`.

- block_times:

  Unique integer vector of length two defining the first and final
  chronological positions of the inclusive block when \`pattern =
  "block_times"\`.

- background:

  Named list describing the background process. Supported distributions
  are:

  \`"normal"\`

  : \`list(dist = "normal", mean = ..., sd = ...)\`

  \`"empirical"\`

  : \`list(dist = "empirical", values = ...)\`

  \`"ar1"\`

  : \`list(dist = "ar1", mean = ..., sd = ..., phi = ...)\`

  \`"fixed"\`

  : \`list(dist = "fixed", value = ...)\`

  For the stationary Gaussian AR(1) generator, \`sd\` is the stationary
  marginal standard deviation and \`phi\` must lie strictly inside
  \`(-1, 1)\`.

  Unknown background fields are rejected explicitly to catch misspelled
  parameter names.

- bounds:

  Optional finite numeric vector \`c(lower, upper)\`. Background values
  are clipped to these bounds before the pattern and values are clipped
  again after the pattern.

  If \`cumulative = TRUE\`, \`bounds\` constrain the exposure
  \*\*increments\*\* before \`cumsum()\`; the final cumulative
  trajectory is not clipped again.

- seed:

  Optional finite integer seed controlling the complete set of \`n\`
  exposure simulations. The caller's random-number state is restored
  after the function returns.

- cumulative:

  Logical. If \`TRUE\`, applies \`cumsum()\` in chronological order
  after background generation, pattern assignment, and bounds.

## Value

An object of class \`"epiexposure_simulated_exposures"\` with:

- \`profiles\`:

  Canonical list containing exactly \`n\` simulated numeric
  chronological profiles, named \`"simulation_1"\`, ...,
  \`"simulation_n"\`.

- \`profile\`:

  Backward-compatible numeric alias when \`n = 1\`; otherwise \`NULL\`.

- \`simulation_data\`:

  Long data frame containing \`simulation\`, \`position\`, \`time\`,
  \`lag\`, and \`profile\`. This is a tabular view of the same simulated
  exposure profiles, not model-parameter samples.

- \`meta\`:

  Common simulation metadata, including \`max_lag\`, chronological
  order, \`n\`, background specification, and RNG contract.

- \`simulation_meta\`:

  List with simulation-specific pattern metadata, such as selected
  chronological positions and corresponding retrospective lags.

## Details

The first profile value corresponds internally to the maximum
retrospective lag and the final profile value corresponds internally to
lag 0:

“\` time_0 = earliest/oldest observation = lag max_lag time_max_lag =
most recent observation = lag 0 “\`

\`simulate_exposures()\` simulates exposure profiles only. It does not
draw model coefficients and therefore has no \`"summary"\`/\`"samples"\`
output switch. The number of simulated exposure profiles is controlled
exclusively by \`n\`.

\## One simulation contract

\`n\` is the only argument controlling how many exposure trajectories
are generated:

“\` x1 \<- simulate_exposures(..., n = 1) length(x1\$profiles) \# 1

x100 \<- simulate_exposures(..., n = 100) length(x100\$profiles) \# 100
“\`

There is no statistical summarization inside this function. If \`n =
100\`, all 100 simulated profiles are retained.

This is deliberately different from \`predict_outcomes()\` and
\`compare_predictions()\`, where \`"samples"\` refers to uncertainty
draws of fitted model parameters and \`"summary"\` summarizes those
parameter-draw predictions.

\## Direct use in \`compare_predictions()\`

The whole returned object can be used directly as one exposure
component:

“\` rain_A \<- simulate_exposures(..., n = 100)

scenario_A \<- list( tmean = 25, rain = rain_A, wetness = 10 ) “\`

\`compare_predictions()\` recognizes \`\$profiles\` as the canonical
collection. If other exposure variables in the same scenario also
contain 100 profiles, they are paired by simulation index:

“\` tmean profile 1 + rain profile 1 + wetness profile 1 tmean profile
2 + rain profile 2 + wetness profile 2 ... “\`

A variable with exactly one profile may be recycled across the common
scenario profile count. No Cartesian product is constructed.

\## Chronological versus retrospective indexing

All user-facing pattern arguments use chronological time:

“\` time = 0 -\> earliest/oldest observation time = max_lag -\> most
recent observation “\`

The corresponding retrospective lag is:

\$\$lag = max\\lag - time.\$\$

No profile reversal is performed after simulation.

\## Background generators

For \`background\$dist = "normal"\`, positions are independent Gaussian
draws.

For \`"empirical"\`, positions are sampled with replacement from the
supplied finite empirical values.

For \`"ar1"\`, the process is initialized from its stationary marginal
distribution and uses innovation SD

\$\$ sd\_{innovation} = sd\sqrt{1-\phi^2}. \$\$

Thus \`background\$sd\` is the stationary marginal SD rather than the
innovation SD.

For \`"fixed"\`, every background position equals \`background\$value\`.

\## Seeds across exposure variables

\`seed\` controls one call to \`simulate_exposures()\`. Reusing exactly
the same seed in separate stochastic calls for different exposure
variables can reuse the same underlying pseudo-random sequence and
unintentionally induce dependence among simulated exposures. Use
different seeds for independently simulated variables unless
synchronized random structure is intentional.

\## Pattern-specific arguments

Pattern arguments that do not belong to the selected mode/pattern are
rejected rather than silently ignored.
