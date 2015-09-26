# LifeTable
The `LifeTable` module contains the function `PeriodLifeTable`, which can be used to calculate life tables in accordance with the methods used by [Human Mortality Database](http://www.mortality.org/) (in particular, see the [Methods Protocol](http://www.mortality.org/Public/Docs/MethodsProtocol.pdf) 38--39). The function is called like `PeriodLifeTable(inframe, sex, intype)`, where `inframe` is a [DataFrame](https://github.com/JuliaStats/DataFrames.jl), `sex` may be given as `1` or `2` for males or females (the calculation of average numbers lived for those dying as infants differs between the sexes), and `intype` can be `"count"` (the default) or `"rate"`.

The first column in `inframe` is assumed to contain the start of the included age intervals (the last row thus includes data for an open interval at the end of life). If `intype="count"`, the second column is assumed to be a vector of age-specific population at risk, and the third column is assumed to be a vector of age-specific death counts. If `intype="rate"`, the second column is assumed to contain age-specific death rates.

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
