using DataArrays, DataFrames

function LifeTable(inframe, sex, intype="count")
	if sex == 1
		aint = 0.045
		acoef = 2.684
	elseif sex == 2
		aint = 0.053
		acoef = 2.8
	end

	nrows = size(inframe, 1)
	outframe = DataFrame()
	outframe[:age] = inframe[1]
	i = outframe[:age][2:nrows] .- outframe[:age][1:nrows-1]

	if intype == "count"
		outframe[:m] = inframe[3]./inframe[2]
	elseif intype == "rate"
		outframe[:m] = inframe[2]
	end
	
	a0 = aint + acoef * outframe[:m][1]
	aend = 1 / outframe[:m][nrows]
	aoth = fill(0.5, nrows-2)
	aclosed = [a0, aoth]
	outframe[:a] = [aclosed, aend]
	

	q(m, a, i) = (i*m) / (1+(1-a)*i*m)
	qclosed = map(q, outframe[:m][1:nrows-1], aclosed, i)
	outframe[:q] = [qclosed, 1]
	outframe[:p] = 1 .- outframe[:q]

	outframe[:l] = cumprod([1, outframe[:p][1:nrows-1]])
	dclosed = outframe[:q][1:nrows-1] .* outframe[:l][1:nrows-1]
	outframe[:d] = [dclosed, outframe[:l][nrows]]

	ld(l, a, d) = l - (1-a) * d
	ldclosed = map(ld, outframe[:l][1:nrows-1], aclosed, dclosed)
	outframe[:ld] = [ldclosed, outframe[:l][nrows]*aend]

	revl = reverse(outframe[:ld])
	tclosed = reverse(cumsum(revl))[1:nrows-1]
	outframe[:t] = [tclosed, outframe[:ld][nrows]]
	outframe[:e] = outframe[:t]./outframe[:l]

	return outframe

end
