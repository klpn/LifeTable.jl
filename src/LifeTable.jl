module LifeTable

using DataArrays, DataFrames

export PeriodLifeTable

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

end
