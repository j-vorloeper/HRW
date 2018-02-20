# Author:
# Jürgen Vorloeper, Hochschule Ruhr West, 2018


import GLPK

"""
 simple interface to glpk
 solve f*x->max under constraints Ax<=b, x>=0,
 with A \in R^(n\times m), b\in R^n, f,lb,ub \in R^m

 #### Fields
 - 'f::Vector{Real}': xx
 - 'A::Matrix{Real}': xx
 - 'b::Vector{Real}': xx
"""
function glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
#function glpk(f::Array{Real,1}, A::Array{Real,2}, b::Array{Real,1})
#function glpk{T}(f::AbstractArray{T,1}, A::AbstractArray{T,2}, b::AbstractArray{T,1})
#function glpk(f::Array{Float64,1}, A::Array{Float64,2}, b::Array{Float64,1})

	eflag = 0

	n = size(A, 1)
	m = size(A, 2)


	lp = GLPK.Prob()
	GLPK.set_prob_name(lp, "sample")
	GLPK.set_obj_dir(lp, GLPK.MAX)


	GLPK.add_rows(lp, n)

	for i=1:n
		GLPK.set_row_bnds(lp, i, GLPK.UP, 0.0, b[i])
	end

	GLPK.add_cols(lp, m)

	for j=1:m
		GLPK.set_col_bnds(lp, j, GLPK.LO, 0.0, 0.0)
		GLPK.set_obj_coef(lp, j, f[j])
	end


	GLPK.load_matrix(lp, sparse(A))
	GLPK.simplex(lp)

	z = GLPK.get_obj_val(lp)

	x = Array{Float64}(m)

	for j=1:m
		x[j] = GLPK.get_col_prim(lp, j)
    end

	return x, z, eflag

end


# simple interface to glpk
# solve LP Ax (<=,==,>=) b, lb<= x <= ub, f*x->max
# with A \in R^(n\times m), b\in R^n, f,lb,ub \in R^m
# with A[i,:] <= b[i], if ctype[i] == "U"
# with A[i,:] == b[i], if ctype[i] == "S",
# with A[i,:] >= b[i], if ctype[i] == "L"
#
# with x[i] integer,    if vartype[i] == "I"
# with x[i] continuous, if vartype[i] == "C"
# with x[i] binary,     if vartype[i] == "B"
function glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64}, ctype::String, vartype::String)

	eflag = 0

	mip = false

	n = size(A, 1)
	m = size(A, 2)


	lp = GLPK.Prob()
	GLPK.set_prob_name(lp, "sample")
	GLPK.set_obj_dir(lp, GLPK.MAX)


	GLPK.add_rows(lp, n)

	for i = 1:n
		if ctype[i] == 'U'
			GLPK.set_row_bnds(lp, i, GLPK.UP, 0.0, b[i])
		end

		if ctype[i] == 'L'
			GLPK.set_row_bnds(lp, i, GLPK.LO, b[i], 0.0)
		end

		if ctype[i] == 'S'
			GLPK.set_row_bnds(lp, i, GLPK.FX, b[i], b[i])
		end
	end


	GLPK.add_cols(lp, m)

	for j = 1:m
		if abs(lb[j]-ub[j])<1e-12
			GLPK.set_col_bnds(lp, j, GLPK.FX, lb[j], ub[j])
		else
			GLPK.set_col_bnds(lp, j, GLPK.DB, lb[j], ub[j])
		end
		GLPK.set_obj_coef(lp, j, f[j])

		if vartype[j] == 'C'
			GLPK.set_col_kind(lp,j, GLPK.CV)
		end

		if vartype[j] == 'I'
			mip = true
			GLPK.set_col_kind(lp, j, GLPK.IV)
		end

		if vartype[j] == 'B'
			mip = true
			GLPK.set_col_kind(lp, j, GLPK.BV)
		end
	end

	GLPK.load_matrix(lp, sparse(A))
	v1 = GLPK.simplex(lp)
	if v1>0
		println("*** glpk: No solution found!")
		eflag = -1
	end

	if mip==true
		v2 = GLPK.intopt(lp)

		if v2>0
			println("*** glpk: No solution found!")
			eflag = -1
		end
	end

	if mip==true
		z = GLPK.mip_obj_val(lp)
	else
		z = GLPK.get_obj_val(lp)
	end


	x = Array{Float64}(m)

	if mip==true
		for j=1:m
			x[j] = GLPK.mip_col_val(lp, j)
		end
	else
		for j=1:m
			x[j] = GLPK.get_col_prim(lp, j)
		end
    end


	return x, z, eflag
end




# simple interface to glpk
# solve LP Ax <= b, lb<= x <= ub, f*x->max
# with A \in R^(n\times m), b\in R^n, f,lb,ub \in R^m
function glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})

	n = size(A, 1)
	m = size(A, 2)

	ctype = repeat("U", n)
	vartype = repeat("C",m)

	return glpk(f, A, b, lb, ub, ctype, vartype)
end
