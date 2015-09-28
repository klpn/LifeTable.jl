using LifeTable
using Base.Test
using DataFrames

# Builds a lifetable for Swedish males 2014 (based on data from Statistics Sweden). 
mse14 = readtable("data/mse14.csv")
ltmse14 = PeriodLifeTable(mse14, 1)

@test_approx_eq_eps ltmse14[:e][1] 80.3453 1e-3

# Performs a regression on Weibull parameters for survival for Swedish 
# females dying of leukemia in age intervals from 30--34 to 80--84 years
# 1951--2012 (based on data from WHO Mortality Database).

leuk_msfs = Dict[]
fse5112tot3084 = readtable("data/fse5112tot3084.csv")
fse5112pop3084 = readtable("data/fse5112pop3084.csv")
fse5112leuk3084 = readtable("data/fse5112leuk3084.csv")

for i in 2:63
	dths = DataFrame(age = fse5112pop3084[1], p = fse5112pop3084[i], d = fse5112tot3084[i])
	lt = PeriodLifeTable(dths, 2, false)
	leuk_prop = fse5112leuk3084[i] ./ fse5112tot3084[i]
	leuk_life = CauseLife(lt, leuk_prop)
	leuk_lt = PeriodLifeTable(leuk_life, 2, false, "rate")
	leuk_msf = MortSurvFit(leuk_lt, fse5112tot3084[i], "weibull", "surv")
	push!(leuk_msfs, leuk_msf)
end

fit_leuk_msfs = MortSurvParamsFit(leuk_msfs)

@test_approx_eq_eps coef(fit_leuk_msfs["fit"])[2] -4.5 1e-1
