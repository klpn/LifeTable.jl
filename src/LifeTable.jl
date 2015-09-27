module LifeTable

using DataArrays, DataFrames, LsqFit, GLM 

export PeriodLifeTable, CauseLife, MortSurvFit, MortSurvParamsFit

function PeriodLifeTable(inframe, sex, intype="count")
	if sex == 1
		aint = 0.045
		acoef = 2.684
		aconst = 0.330
	elseif sex == 2
		aint = 0.053
		acoef = 2.8
		aconst = 0.350
	end

	nrows = size(inframe, 1)
	age = inframe[1]
	i = age[2:nrows] .- age[1:nrows-1]

	if intype == "count"
		m = inframe[3]./inframe[2]
	elseif intype == "rate"
		m = inframe[2]
	end
	
	if age[1] == 0
		if m[1] >= 0.107
			a0 = aconst
		else
			a0 = aint + acoef * m[1]
		end
	else
		a0 = 0.5
	end

	aend = 1 / m[nrows]
	aoth = fill(0.5, nrows-2)
	aclosed = [a0, aoth]
	a = [aclosed, aend]
	

	q(m, a, i) = (i*m) / (1+(1-a)*i*m)
	qclosed = map(q, m[1:nrows-1], aclosed, i)
	q = [qclosed, 1]
	p = 1 .- q

	l = cumprod([1, p[1:nrows-1]])
	dclosed = q[1:nrows-1] .* l[1:nrows-1]
	d = [dclosed, l[nrows]]

	ld(l, a, d, i) = (l - (1-a) * d) * i
	ldclosed = map(ld, l[1:nrows-1], aclosed, dclosed, i)
	ld = [ldclosed, l[nrows]*aend]

	revl = reverse(ld)
	tclosed = reverse(cumsum(revl))[1:nrows-1]
	t = [tclosed, ld[nrows]]
	e = t./l
	
	outframe = DataFrame(age = age, m = m, a = a, q = q, p = p, l = l, d = d,
	ld = ld, t = t, e = e)

	return outframe

end

function CauseLife(lifetable, causefreq)
	dca = causefreq .* lifetable[:d]
	l0ca = reverse(cumsum(reverse(dca)))
	fca = l0ca ./ lifetable[:l]
	mca = (causefreq .* lifetable[:m])
	mcacoh = mca ./ fca

	outframe = DataFrame(age = lifetable[:age], m = mcacoh, f = fca)
	return outframe
end

function MortSurvFit(lifetable, numbdeaths, func, functype)
	if func == "gompertz"
		ratemodel(x, p) = p[1] * exp(p[2].*x)
		survmodel(x, p) = exp(-(p[1]/p[2]) * (exp(p[2].*x)-1))
		p = [exp(-18), 0.1]
	elseif func == "weibull"
		ratemodel(x, p) = (p[1]/p[2]) .* (x./p[2]).^(p[1]-1)
		survmodel(x, p) = exp(-(x./p[2]).^p[1]) 
		p = [10.0, 80.0]
	end
	
	if functype == "rate"
		ycol = lifetable[:m]
		model = ratemodel
	elseif functype == "surv"
		ycol = lifetable[:l]
		model = survmodel
	end

	xarr = convert(Array{Int}, lifetable[:age])
	yarr = convert(Array{Float64}, ycol)
	warr = sqrt(convert(Array{Int}, numbdeaths))

	fit = curve_fit(model, xarr, yarr, warr, p)
	return ["fit" => fit, "func" => func, "functype" => functype]
end

function MortSurvParamsFit(msfits)
	p1arr = Float64[]
	p2arr = Float64[]

	for msfit in msfits
		push!(p1arr, msfit["fit"].param[1]) 
		push!(p2arr, msfit["fit"].param[2])
	end
	
	func = msfits[1]["func"]
	functype = msfits[1]["functype"]
	if func == "gompertz"
		ratetrans_x(p1, p2) = p2
		ratetrans_y(p1, p2) = log(p1)
		survtrans_x(p1, p2) = p2
		survtrans_y(p1, p2) = log(p1/p2)
	elseif func == "weibull"
		ratetrans_x(p1, p2) = p1 - 1
		ratetrans_y(p1, p2) = log(p1/p2) - (p1-1) * log(p2)
		survtrans_x(p1, p2) = p1
		survtrans_y(p1, p2) = -p1 * log(p2) 
	end
	
	if functype == "rate"
		trans_x = ratetrans_x
		trans_y = ratetrans_y
	elseif functype == "surv"
		trans_x = survtrans_x
		trans_y = survtrans_y
	end

	xarr = map(trans_x, p1arr, p2arr)
	yarr = map(trans_y, p1arr, p2arr)
	trans_params = DataFrame(p1 = p1arr, p2 = p2arr, X = xarr, Y = yarr)
	
	fit = glm(Y~X, trans_params, Normal(), IdentityLink())
	return ["trans_params" => trans_params, "fit" => fit]
end

end
