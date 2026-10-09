# Compare chronological exposure profiles

Compares two or more exposure profiles, simulated scenarios, or
epidemic-specific exposure profiles while preserving the EpiExposure
chronological profile convention.

## Usage

``` r
compare_exposures(
  exposures = NULL,
  exposure1 = NULL,
  exposure2 = NULL,
  time = NULL,
  group = NULL,
  agg_fun = c("median", "mean"),
  mode = c("global", "timewise"),
  cor_method = c("pearson", "spearman", "kendall"),
  q = 0.8,
  eps = 1e-12
)
```

## Arguments

- exposures:

  Named list of numeric exposure profiles. Every element must represent
  the same exposure variable, in the same units, over the same
  chronological positions, and all profiles must have equal length.

  Profiles must be stored from the earliest/oldest observation to the
  most recent observation. A single numeric vector is still accepted for
  backward compatibility, but at least two profiles or aggregated groups
  are required for an actual comparison.

- exposure1:

  Legacy numeric exposure vector. It must be supplied together with
  \`exposure2\` when \`exposures = NULL\`.

- exposure2:

  Legacy numeric exposure vector. It must be supplied together with
  \`exposure1\` when \`exposures = NULL\`.

  If \`exposures\` is supplied, \`exposure1\` and \`exposure2\` must
  both be \`NULL\`; legacy arguments are not silently ignored.

- time:

  Optional numeric vector of chronological time positions. It must have
  the same length as every exposure profile, contain unique finite
  values, be strictly increasing, and be equally spaced.

  Equal spacing is required because EpiExposure profiles represent
  complete regular exposure profiles and several global metrics treat
  each temporal position equally.

  If \`NULL\`, \`0:(n - 1)\` is used.

- group:

  Optional group membership for the supplied profiles. When provided,
  profiles are first aggregated within each group at every time position
  using \`agg_fun\`, and pairwise comparisons are made between the
  resulting group-level profiles.

  An unnamed vector is matched by profile position. A named vector is
  safer: its names must match \`names(exposures)\` exactly and are used
  to align group membership to profiles before aggregation.

- agg_fun:

  Character. \`"median"\` or \`"mean"\`. Used only when \`group\` is
  supplied. The default is \`"median"\`.

- mode:

  Character. \`"global"\` or \`"timewise"\`.

  \`"global"\`

  : Returns one row per pairwise profile comparison with distance,
    magnitude, correlation, peak-timing, center-of-mass, and
    high-exposure overlap summaries.

  \`"timewise"\`

  : Returns one row per chronological position and pairwise comparison,
    retaining the complete exposure profiles.

- cor_method:

  Character correlation method used in \`mode = "global"\`. Available
  options are:

  \- \`"pearson"\` (default): linear association between the original
  exposure values; - \`"spearman"\`: monotonic association calculated
  from exposure ranks; - \`"kendall"\`: Kendall's rank correlation
  (tau), also measuring monotonic association and often useful for small
  samples or many non-linear monotonic relationships.

  The selected method changes only the \`corr\` summary. All distance,
  magnitude, timing, center-of-mass, and overlap metrics are unchanged.

  Correlation is returned as \`NA\` when either profile is constant
  because none of the supported correlation coefficients is defined in
  that case.

- q:

  Finite numeric scalar between 0 and 1. Each profile's own \`q\`
  quantile defines its high-exposure positions for the Jaccard-style
  \`overlap_above_q\` metric. The default is \`0.8\`.

- eps:

  Positive finite numeric tolerance. It is used only to decide when a
  time-wise ratio denominator or a non-negative profile's total
  center-of-mass weight is numerically indistinguishable from zero.

  \`eps\` is \*\*not added to exposure values\*\* and therefore never
  modifies a profile before comparison.

## Value

A data frame.

With \`mode = "timewise"\`, the output contains:

\- \`exposure1\`, \`exposure2\`: compared profile/group names; -
\`position\`: chronological position \`1, ..., n\`; - \`time\`:
chronological time supplied by the user; - \`lag\`: positional
retrospective lag, with lag 0 at the most recent observation; -
\`value1\`, \`value2\`; - \`diff = value1 - value2\`; - \`abs_diff =
abs(diff)\`; - \`ratio = value1 / value2\` when \`abs(value2) \> eps\`,
otherwise \`NA\`; - \`ratio_defined\`: whether that numerical ratio was
computed.

With \`mode = "global"\`, the output contains one row per pair and
includes:

\- discrete distances \`L1\`, \`L2\`, and \`Linf\`; - \`MAE\` and
\`RMSE\`, which are length-normalized profile discrepancies; -
individual and between-profile sums, means, and SDs; - the selected
Pearson, Spearman, or Kendall correlation when both profiles vary; -
peak values and peak timing; - explicit information about tied maxima; -
temporal centers of mass when the corresponding profile is non-negative
and has positive total mass; - profile-specific high-exposure thresholds
and Jaccard overlap.

Output attributes record chronological order, common time step,
comparison direction, grouping information, and the ratio/center-of-mass
contracts.

## Details

Each profile must be ordered from the \*\*oldest/earliest\*\* exposure
observation to the \*\*most recent\*\* observation. Consequently, for a
profile of length \`n\`, the final element corresponds to retrospective
lag 0 and the first element corresponds to positional lag \`n - 1\`.

\`compare_exposures()\` is descriptive. It compares exposure profiles
themselves and does not estimate DLNM effects, disease risks, or model
predictions.

\## Comparison direction

Pairwise differences are directional:

\$\$difference = exposure1 - exposure2.\$\$

Therefore positive \`diff\`, \`diff_sum\`, or \`diff_mean\` values
indicate larger exposure in the first named profile than in the second.

\`timing_shift_peak\` and \`center_of_mass_shift\` follow the same
direction:

\$\$time_1 - time_2.\$\$

Positive values therefore mean that the corresponding timing summary
occurs later in \`exposure1\`.

\## Discrete distance metrics

For paired profile values \\x\_{1t}\\ and \\x\_{2t}\\,

\$\$L1 = \sum_t \|x\_{1t} - x\_{2t}\|,\$\$

\$\$L2 = \sqrt{\sum_t (x\_{1t} - x\_{2t})^2},\$\$

and

\$\$Linf = \max_t \|x\_{1t} - x\_{2t}\|.\$\$

\`L1\` and \`L2\` increase with the number of represented time
positions. \`MAE\` and \`RMSE\` are therefore also returned when
comparisons across equal-unit profiles of different lengths are needed.

These distances retain the units/scaling of the exposure profile and
should only be compared across profiles of the same exposure variable
and units.

\## Time-wise ratios

Earlier EpiExposure code calculated

\$\$x_1 / (x_2 + eps),\$\$

which changed every denominator, including perfectly valid non-zero
values. The current implementation never alters observed exposure
values. It returns

\$\$x_1/x_2\$\$

only when \`abs(x2) \> eps\`; otherwise the ratio is \`NA\`.

A numerical ratio is not automatically scientifically meaningful. It
should only be interpreted as a relative exposure contrast when the
exposure is measured on a ratio scale with a meaningful zero. For
example, ratios of temperatures expressed in degrees Celsius generally
should not be interpreted as relative temperature effects.

\## Correlation method

\`cor_method\` controls only the global \`corr\` column:

\- Pearson measures linear co-variation of the original profile
values; - Spearman is Pearson correlation of the ranks and measures
monotonic association; - Kendall returns Kendall's tau, based on
concordant and discordant pairs.

Correlation describes similarity in temporal pattern, not similarity in
absolute exposure magnitude. For example, two profiles can have
correlation close to 1 while one is systematically much larger than the
other. Distance and magnitude metrics should therefore be interpreted
alongside \`corr\`.

No correlation is causal and no p-value is computed by this descriptive
function.

\## Peak timing and ties

\`peak_time1\` and \`peak_time2\` retain the first chronological
occurrence of the maximum for backward compatibility. The function
additionally returns first and last peak times, the midpoint of their
temporal range, the number of maximum positions, and indicators of tied
maxima. Thus a plateau or repeated equal maximum is never hidden by
\`which.max()\`.

\`timing_shift_peak\` uses the first peak times for backward
compatibility. \`timing_shift_peak_midpoint\` compares the midpoints of
the full peak-time ranges and is usually more informative when maxima
are tied.

\## Temporal center of mass

For a non-negative profile with positive total exposure, the temporal
center of mass is

\$\$COM = \frac{\sum_t time_t x_t}{\sum_t x_t}.\$\$

Earlier versions silently replaced negative exposure values by zero
before computing COM. That changes the supplied profile and is no longer
done.

If a profile contains any negative value, or its non-negative total mass
is not greater than \`eps\`, its COM is returned as \`NA\` and
\`com_defined\` is \`FALSE\`. Even for non-negative data, COM is
scientifically interpretable only when zero and exposure magnitude have
a meaningful weighting interpretation.

\## High-exposure overlap

Each profile obtains its own threshold:

\$\$Q_j(q) = quantile(x_j, q).\$\$

High-exposure positions satisfy \\x\_{jt} \>= Q_j(q)\\. The returned
\`overlap_above_q\` is the Jaccard index

\$\$ \frac{\|H_1 \cap H_2\|}{\|H_1 \cup H_2\|}. \$\$

This measures whether the two profiles experience their \*\*own
relatively high exposure\*\* at the same times. It is not an overlap
above one common absolute exposure threshold. Ties at the quantile
threshold can make the number of high positions larger than exactly
\`(1 - q) \* n\`; counts and thresholds are returned explicitly.

\## Group aggregation

With \`group\`, aggregation occurs independently at every chronological
position before any pairwise metric is calculated. The resulting
comparison therefore describes representative group-level profiles, not
a distribution of pairwise epidemic-level differences. No uncertainty
interval or inferential test is produced by this function.
