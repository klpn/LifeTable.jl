using LifeTable
using Test
using DataFrames
using CSV 
using GLM

# Builds a lifetable for Swedish males 2014 (based on data from Statistics Sweden). 
mse14 = CSV.read("data/mse14.csv")
ltmse14 = periodlifetable(mse14, 1)

@test isapprox(ltmse14[:e][1], 80.3453; atol=1e-3)

# Performs a regression on Weibull parameters for survival for Swedish 
# females dying of leukemia in age intervals from 30--34 to 80--84 years
# 1951--2012 (based on data from WHO Mortality Database).

leuk_msfs = Dict[]
fse5112tot3084 = CSV.read("data/fse5112tot3084.csv")
fse5112pop3084 = CSV.read("data/fse5112pop3084.csv")
fse5112leuk3084 = CSV.read("data/fse5112leuk3084.csv")

for i in 2:63
	dths = DataFrame(age = fse5112pop3084[1], p = fse5112pop3084[i], d = fse5112tot3084[i])
	lt = periodlifetable(dths, 2, false)
	leuk_prop = fse5112leuk3084[i] ./ fse5112tot3084[i]
	leuk_life = causelife(lt, leuk_prop)
	leuk_lt = periodlifetable(leuk_life, 2, false, "rate")
	leuk_msf = mortsurvfit(leuk_lt, fse5112tot3084[i], "weibull", "surv")
	push!(leuk_msfs, leuk_msf)
end

fit_leuk_msfs = mortsurvparamsfit(leuk_msfs)

@test isapprox(coef(fit_leuk_msfs["fit"])[2], -4.5; atol=1e-1)
