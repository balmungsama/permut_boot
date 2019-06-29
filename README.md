# permut&boot.R
Conduct permutation tests with bootstrapped confidence intervals

This package is designed to provide easy and user-friendly access to permuted differences tests and bootstrap-resampled confidence intervals.
the advantage of this package over others is the relative ease in which it can be implemented, with very few user-defined parameters necessary to conduct the statistical analysis. The other advantage is the inclusion of a one-sample permutation test function, which allows to an estimation of the p-value and reliability that the provided data is different from zero. This is accomplished by estimating a null binomial distrobition against which to compare the user-provided data.
