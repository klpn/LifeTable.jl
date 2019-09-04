# LifeTable

[![Build Status](https://travis-ci.org/klpn/LifeTable.jl.svg?branch=master)](https://travis-ci.org/klpn/LifeTable.jl)

The `LifeTable` module contains the function `periodlifetable`, which can be used to calculate life tables in accordance with the methods used by [Human Mortality Database](http://www.mortality.org/) (in particular, see the [Methods Protocol](http://www.mortality.org/Public/Docs/MethodsProtocol.pdf) 38--39). The function is called like `periodlifetable(inframe, sex, openend, intype)`, where `inframe` is a [DataFrame](https://github.com/JuliaStats/DataFrames.jl), `sex` may be given as `1` or `2` for males or females (the calculation of average numbers lived for those dying as infants differs between the sexes), and `intype` can be `"count"` (the default) or `"rate"`.

The first column in `inframe` is assumed to contain the start of the included age intervals. If `openend=true` (default), the last row is assumed to include data for an open interval at the end of life; if `openend=false` it is assumed to include a closed interval of the same size as the second last row. If `intype="count"`, the second column is assumed to be a vector of age-specific population at risk, and the third column is assumed to be a vector of age-specific death counts. If `intype="rate"`, the second column is assumed to contain age-specific death rates.

The function returns a new DataFrame with the life table.

| Column | Description
| ------ | -----------
| age | Start of age interval.
| m | Age-specific death rates.
| a | Expected numbers of years lived at a given age for someone dying at that age. 
| q | Age-specific death probabilities.
| p | Age-specific survival probabilities.
| l | Probability of survival to a given age.
| d | Distribution of deaths by age.
| ld | Number of years lived at a given age.
| t | Remaining number of years to live at a given age.
| e | Life expectancy at a given age.

The function `causelife(lifetable, causefreq)` takes a lifetable returned by `periodlifetable` and a vector of the same length as the number of rows in the life table, with the proportion of deaths at a given age caused by a specific disease (or another cause). It returns a DataFrame with the following columns:


| Column | Description
| ------ | -----------
| age | Start of age interval (same as the input life table).
|m | Age-specific death rates normalized to the frequency of the cause of death (the next column).
|f | The proportion of the population surviving to a given age that will die of the specific cause.

The DataFrame returned by `causelife` can be used as input to `periodlifetable` (with `intype="rate"`) to compute a new lifetable for the subpopulation dying of a specific cause.

The function `mortsurvfit(lifetable, numbdeaths, func, functype)` fits the `m` column (if `functype="rate"`) or the `l` column (if `functype="surv"`) against the `age` column in a DataFrame returned by `periodlifetable` with a Gompertz function (if `func="gompertz"`) or a two-parameter Weibull function (if `func="weibull"`). See e.g. [Juckett and Rosenberg (1993)](http://www.ncbi.nlm.nih.gov/pubmed/8377524) for a discussion of these different functions and their applications to human mortality and survival data. The columns are fitted using [LsqFit](https://github.com/JuliaOpt/LsqFit.jl), and the data points are weighted by the square root of the age-specific number of deaths (given in `numbdeaths`). All input columns are assumed to be DataArrays, and are converted to ordinary arrays in order to work with LsqFit. `mortsurvfit` returns a dictionary with the `func` and `functype` parameters and the `LsqFitResult` returned by the curve fitting.

`mortsurvparamsfit(msfits)` takes a list of dictionaries returned by `mortsurvfit` and performs linear regression (using [GLM](https://github.com/JuliaStats/GLM.jl)) on transformations of the coefficients in the list. This can be used to determine intersections between mortality or survival curves. A dictionary containing (1) a DataFrame with the original parameters and transformations and (2) the results from the regression is returned.
