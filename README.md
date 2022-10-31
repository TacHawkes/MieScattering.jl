# MieScattering.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/github/TacHawkes/MieScattering.jl?svg=true)](https://ci.appveyor.com/project/TacHawkes/MieScattering-jl)
[![Coverage](https://codecov.io/gh/TacHawkes/MieScattering.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/TacHawkes/MieScattering.jl)

MieScattering.jl is a pure-Julia package for calculation of Mie scattering. This package is basically a 1:1 port of the great [miepython](https://github.com/scottprahl/miepython) Python package written by Scott Prahl. I highly recommed to read his [documentation](https://miepython.readthedocs.io/en/latest/index.html) which also gives a lot of mathematical insights into the topic. This package basically follows the same syntax as `miepython` but with small tweaks to adapt to Julia specifics. This package gives comparable or even better (thanks to real multi-threading) performance than the JIT-version of `miepython` using `numba`.
The test suite has been adapted as well, so the whole package is numerically tested against Wiscombe's code and of course `miepython` itself.

## Documentation

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://tachawkes.github.io/MieScattering.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tachawkes.github.io/MieScattering.jl/stable)

## Installation

Install the package using the package manager:

```julia
] add https://github.com/TacHawkes/MieScattering.jl
```

## Alternatives

This is not the first Julia package for calculating Mie scattering. The whole purpose of this package is to have all features of `miepython` without having to call Python from Julia.

- [JLMie](https://github.com/Hinamoooon/jlmie)
- [MieScatter] (https://github.com/dronir/MieScatter.jl)