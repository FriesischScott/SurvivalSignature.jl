# SurvivalSignature.jl

![CI](https://github.com/friesischscott/SurvivalSignature.jl/actions/workflows/ci.yml/badge.svg) [![codecov](https://codecov.io/gh/FriesischScott/SurvivalSignature.jl/branch/master/graph/badge.svg?token=1967M26ATE)](https://codecov.io/gh/FriesischScott/SurvivalSignature.jl) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4306113.svg)](https://doi.org/10.5281/zenodo.4306113)


Julia package for the computation of survival signatures as introduced by [Coolen et al. (2013)](https://doi.org/10.1007/978-3-642-30662-4_8).

In addition to the regular analytical computations, this package contains a Monte Carlo simulation based algorithm to approximate the survival signature for systems where the computational demand for the standard approach is too high. If you use this package for a publication, please cite [Behrensdorf et al. (2021)](https://doi.org/10.1016/j.ress.2021.107935).

## References

Behrensdorf, J., Regenhardt, T.-E., Broggi, M., Beer, M. (2021) Numerically efficient computation of the survival signature for the reliability analysis of large networks, *Reliability Engineering & System Safety*, 107935, https://doi.org/10.1016/j.ress.2021.107935.

Coolen F.P.A., Coolen-Maturi T. (2013) Generalizing the Signature to Systems with Multiple Types of Components. *In: Zamojski W., Mazurkiewicz J., Sugier J., Walkowiak T., Kacprzyk J. (eds) Complex Systems and Dependability. Advances in Intelligent and Soft Computing*, 170. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-30662-4_8