Julia HRW module
=================




HRW.jl is written for teaching purposes at the Hochschule Ruhr West (University of Applied Sciences) at Mülheim an der Ruhr, Germany. It contains a wrapper for [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl) and simple ODE solvers.



## Installation

The package is unregistered. It can be installed with `Pkg.clone`.

```
julia> Pkg.clone("https://github.com/j-vorloeper/HRW.git")
```




## Documentation glpk-Interface

`glpk` solves (Mixed-Integer) Linear Programs.

```
x, z, eflag = glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64}, ctype::String, vartype::String)

# solves
# Ax <= b, f*x -> max
# lb <= x <= ub
#
# with A[i,:] <= b[i], if ctype[i] == "U"
# with A[i,:] == b[i], if ctype[i] == "S",
# with A[i,:] >= b[i], if ctype[i] == "L"
#
# with x[i] integer,    if vartype[i] == "I"
# with x[i] continuous, if vartype[i] == "C"
# with x[i] binary,     if vartype[i] == "B"
```

```
x, z, eflag = glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})

# solves
# Ax <= b, f*x -> max
# lb <= x <= ub
```

```
x, z, eflag = glpk(f::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})

# solves
# Ax <= b, f*x -> max
# 0 <= x
```

## Documentation ODE-Solvers

`solve_dgl` solves ordinary system of differential equations `y'(t) = f(t,y(t)), y(t_0) = y0`

```
solve_dgl(methode::String, fh::Function, z0::Vector{Float64}, tspan::Vector{Float64}, h::Float64)
```

The parameteras are as follows
- `methode`: if `expl_Euler` selects explicit Euler scheme, `verb_Euler` selects improved Euler scheme, `RK4` selects Runge-Kutta scheme, `impl_Euler` selects implicit Euler scheme,
- `fh`is function handle to right hand side
- `tspan` determines time range
- `h`is constant time stepping



## Project Status

The package is tested against Julia `0.5` on Windows.


## References

W. Dahmen and A. Reusken: "Numerik für Ingenieure und Naturwissenschaftler", Springer, 2008.
