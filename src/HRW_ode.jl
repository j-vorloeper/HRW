# Author:
# Jürgen Vorloeper, Hochschule Ruhr West, 2018


import NLsolve

function solve_dgl(methode::String, fh::Function, z0::Vector{Float64}, tspan::Vector{Float64}, h::Float64)

const n = round(Int64, (tspan[end]-tspan[1])/h)
zj = copy(z0)

t = zeros(n+1)

const N = length(z0)

z = Array{Array{Float64,1}}(n+1)

z[1] = z0

for k=2:length(z)
	z[k] = zeros(N)
end


if methode == "expl_Euler"
    for j=0:n-1
		tj = tspan[1] + j*h
        zj += h * fh(tj,zj)

        t[j+2] = tj+h
        z[j+2] = zj
    end
elseif methode == "verb_Euler"
    for j=0:n-1
        tj = tspan[1] + j * h
        k1 = zj + h/2*fh(tj,zj)
        zj += h * fh(tj+h/2, k1)

        t[j+2] = tj+h
        z[j+2] = zj
    end
elseif methode == "RK4"
    for j=0:n-1
        tj = tspan[1] + j * h
        k1 = fh(tj,zj)
        k2 = fh(tj+h/2,zj+h/2*k1)
        k3 = fh(tj+h/2,zj+h/2*k2)
        k4 = fh(tj+h,zj+h*k3)
        zj += h/6 * (k1+2*k2+2*k3+k4)

        t[j+2] = tj+h
        z[j+2] = zj
    end
elseif methode== "impl_Euler"
	for j=0:n-1
		tj = tspan[1] + j*h
		zj = z[j+1]

		function ff!(ffx::Vector{Float64},x::Vector{Float64})
			ffx[:] = (x-zj-h*fh(tj+h,x))[:]
		end

		SR = NLsolve.nlsolve(ff!, zj)
		zj = SR.zero

		t[j+2] = tj+h
		z[j+2] = zj
	end
end

return t, z

end
