module LifeTable

using DataFrames, LsqFit, GLM 

export periodlifetable, causelife, mortsurvfit, predict, mortsurvparamsfit

function periodlifetable(inframe, sex, openend=true, intype="count")
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

	if openend
		i = age[2:nrows] .- age[1:nrows-1]
	else
		ifirst = age[2:nrows] .- age[1:nrows-1]
		iend = ifirst[nrows-1]
		i = [ifirst; iend]
	end

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
	
	q(m, a, i) = (i*m) / (1+(1-a)*i*m)
	ld(l, a, d, i) = (l - (1-a) * d) * i

	if openend
		aend = 1 / m[nrows]
		aoth = fill(0.5, nrows-2)
		aclosed = [a0; aoth]
		a = [aclosed; aend]
		qclosed = map(q, m[1:nrows-1], aclosed, i)
		q = [qclosed; 1]
	else
		aoth = fill(0.5, nrows-1)
		a = [a0; aoth]
		q = map(q, m, a, i)
	end

	p = 1 .- q
	l = cumprod([1; p[1:nrows-1]])
	d = q .* l

	if openend 
		ldclosed = map(ld, l[1:nrows-1], aclosed, d[1:nrows-1], i)
		ld = [ldclosed; l[nrows]*aend]
		revl = reverse(ld)
		tclosed = reverse(cumsum(revl))[1:nrows-1]
		t = [tclosed; ld[nrows]]
	else
		ld = map(ld, l, a, d, i)
		revl = reverse(ld)
		t = reverse(cumsum(revl))
	end

	e = t./l
	
	outframe = DataFrame(age = age, m = m, a = a, q = q, p = p, l = l, d = d,
	ld = ld, t = t, e = e)

	return outframe

end

function causelife(lifetable, causefreq)
	dca = causefreq .* lifetable[:d]
	l0ca = reverse(cumsum(reverse(dca)))
	fca = l0ca ./ lifetable[:l]
	mca = (causefreq .* lifetable[:m])
	mcacoh = mca ./ fca

	outframe = DataFrame(age = lifetable[:age], m = mcacoh, f = fca)
	return outframe
end

gomprate(x,p) = p[1] * exp.(p[2].*x)
gompsurv(x,p) = exp.(-(p[1]/p[2]) .* (exp.(p[2].*x)-1))
weibrate(x,p) = (p[1]/p[2]) .* (x./p[2]).^(p[1]-1)
weibsurv(x,p) = exp.(-(x./p[2]).^p[1])
# Functions for blocks with redundant components with exponentially distributed lifespans, 
# see Gavrilov and Gavrilova, Journal of Theoretical Biology (2001).
gavrblrate(x,p) = (p[1]*p[2]*exp.(-p[2].*x).*(1-exp.(-p[2].*x)).^(p[1]-1))./
	(1-(1-exp.(-p[2].*x)).^p[1])
gavrblsurv(x,p) = 1-(1-exp.(-p[2].*x)).^p[1]

models = Dict(
	"gompertz" => Dict("rate" => gomprate, "surv" => gompsurv, 
		"p" => [exp(-18); 0.1]),
	"weibull" => Dict("rate" => weibrate, "surv" => weibsurv,
		"p" => [10.0; 80.0]),
	"gavrblock" => Dict("rate" => gavrblrate, "surv" => gavrblsurv,
		"p" => [10.0; 0.03]),
	)

function mortsurvfit(lifetable, numbdeaths, func, functype)
	ycols = Dict("rate" => :m, "surv" => :l)

	xarr = convert(Array{Int}, lifetable[:age])
	yarr = convert(Array{Float64}, lifetable[ycols[functype]])
	warr = sqrt.(convert(Array{Int}, numbdeaths))
	model = models[func][functype]

	fit = curve_fit(model, xarr, yarr, warr, models[func]["p"])
	return Dict("fit" => fit, "func" => func, "functype" => functype, "model" => model)
end

predict(x, msfit) = msfit["model"](x, msfit["fit"].param)

gompratetr_x(p1, p2) = p2
gompratetr_y(p1, p2) = log(p1) 
gompsurvtr_x(p1, p2) = p2
gompsurvtr_y(p1, p2) = log(p1/p2) 
weibratetr_x(p1, p2) = p1 - 1
weibratetr_y(p1, p2) = log(p1/p2) - (p1-1) * log(p2)
weibsurvtr_x(p1, p2) = p1
weibsurvtr_y(p1, p2) = -p1 * log(p2) 

transmodels = Dict(
	"gompertz" =>  Dict("rate" => [gompratetr_x; gompratetr_y],
		"surv" => [gompsurvtr_x; gompsurvtr_y]),
	"weibull" =>  Dict("rate" => [weibratetr_x; weibratetr_y],
		"surv" => [weibsurvtr_x; weibsurvtr_y]),
	)
	
function mortsurvparamsfit(msfits)
	p1arr = Float64[]
	p2arr = Float64[]

	for msfit in msfits
		push!(p1arr, msfit["fit"].param[1]) 
		push!(p2arr, msfit["fit"].param[2])
	end
	
	func = msfits[1]["func"]
	functype = msfits[1]["functype"]

	trans_x = transmodels[func][functype][1]
	trans_y = transmodels[func][functype][2]

	xarr = map(trans_x, p1arr, p2arr)
	yarr = map(trans_y, p1arr, p2arr)
	trans_params = DataFrame(p1 = p1arr, p2 = p2arr, X = xarr, Y = yarr)
	
	fit = glm(@formula(Y~X), trans_params, Normal(), IdentityLink())
	return Dict("trans_params" => trans_params, "fit" => fit)
end

end
