


# Julia HRW module


HRW.jl is written by [Jürgen Vorloeper](<juergen.vorloeper@hs-ruhrwest.de>) for teaching purposes at the [Hochschule Ruhr West](https://www.hochschule-ruhr-west.de/). It contains a wrapper for [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl) and basic ODE solvers.



## Installation

The package is unregistered. It can be installed with `Pkg.clone`.

```
julia> Pkg.clone("https://github.com/j-vorloeper/HRW.git")
```



## Documentation glpk-Interface

`glpk` solves (Mixed-Integer) Linear Programs f&sdot;x&rarr;max under the constraints Ax (&le;,=,&ge;)b, lb&le;x&le;ub. Components of x are allowed to be real, integer or binary. It is a wrapper around [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl) and provides a very simple and convenient  interface to the [GNU Linear Programming Kit ](http://www.gnu.org/software/glpk) library.

The function
```
x, z, eflag = glpk(f, A, b, lb, ub, ctype, vartype)
```
solves f&sdot;x&rarr;max under the constraints Ax (&le;,=,&ge;)b, lb&le;x&le;ub, with
-  A<sub>i</sub>&sdot;x &le;b<sub>i</sub> if `ctype[i] == "U" `
-  A<sub>i</sub>&sdot;x =b<sub>i</sub> if `ctype[i] == "S" `
-  A<sub>i</sub>&sdot;x &ge;b<sub>i</sub> if `ctype[i] == "L" `

(where A<sub>i</sub> denotes the i-th row of A) and
- x<sub>i</sub>&isin; Z, if `vartype[i] == "I" `
- x<sub>i</sub>&isin; R, if `vartype[i] == "C" `
- x<sub>i</sub>&isin; {0,1}, if `vartype[i] == "B" `

Input arguments are of type `f::Vector{Float64}`, `A::Matrix{Float64}`, `b::Vector{Float64}`, `lb::Vector{Float64}`, `ub::Vector{Float64}`, `ctype::String`, `vartype::String`.

The function
```
x, z, eflag = glpk(f, A, b, lb, ub)
```
solves f&sdot;x&rarr;max under the constraints Ax &le;b, lb&le;x&le;ub with x&isin;R<sup>m</sup>.

The function
```
x, z, eflag = glpk(f, A, b)
```
solves f&sdot;x&rarr;max under the constraints Ax &le;b, 0&le;x&le; x&isin;R<sup>m</sup>

Return values are:
- `x::Vector{Float64}` optimal solution x<sup>&ast;</sup>&isin;R<sup>n</sup>,
- `z::Float64` is equal to f&sdot;x<sup>&ast;</sup>,
- `eflag::Int64` indicates an error, if its value is nonzero.

#### Example 1.1

Solve

4x<sub>1</sub> + 3x<sub>2</sub>&rarr;max

with constraits

<table border="0" rules="rows">
<tr>
<td>4x<sub>1</sub> + 2x<sub>2</sub></td> <td>&le;</td> <td>20</td> </tr>
<tr> <td>2x<sub>1</sub> + 4x<sub>2</sub></td> <td>&le;</td> <td>18</td> </tr>
<tr> <td>x<sub>1</sub> + x<sub>2</sub></td> <td>&le;</td> <td>6</td>
</tr>
</table>

and  0&le; x<sub>1</sub>, x<sub>2</sub> with the following code:

```
import HRW

A = [
    4.0 2.0;
    2.0 4.0;
    1.0 1.0
    ]
b = [20.0, 18.0, 6.0]
f = [4.0, 3.0]

x, z = HRW.glpk(f,A,b)

@printf("Solution: z = %4.2f  x1 = %4.2f  x2 = %4.2f\n", z, x[1], x[2])
```

#### Example 1.2

Solve

x<sub>1</sub> + x<sub>2</sub>  + x<sub>3</sub>   + x<sub>4</sub>  + x<sub>5</sub>   + x<sub>6</sub>   &rarr;min

with constraits

<table border="0" rules="rows">
<tr>
<td>x<sub>1</sub> + x<sub>6</sub></td> <td>&ge;</td> <td>3</td>
</tr>
<tr>
<td>x<sub>1</sub> + x<sub>2</sub></td> <td>&ge;</td> <td>8</td>
</tr>
<tr>
<td>x<sub>2</sub> + x<sub>3</sub></td> <td>&ge;</td> <td>19</td>
</tr>
<tr>
<td>x<sub>3</sub> + x<sub>4</sub></td> <td>&ge;</td> <td>8</td>
</tr>
<tr>
<td>x<sub>4</sub> + x<sub>5</sub></td> <td>&ge;</td> <td>14</td>
</tr>
<tr>
<td>x<sub>5</sub> + x<sub>6</sub></td> <td>&ge;</td> <td>5</td>
</tr>
</table>

where  x<sub>1</sub>,..., x<sub>6</sub> are required to be non-negative integers, with the following code:


```
import HRW

f = -1.0*ones(6)

A = [
    -1.0  0.0  0.0  0.0  0.0 -1.0;
    -1.0 -1.0  0.0  0.0  0.0  0.0;
     0.0 -1.0 -1.0  0.0  0.0  0.0;
     0.0  0.0 -1.0 -1.0  0.0  0.0;
     0.0  0.0  0.0 -1.0 -1.0  0.0;
     0.0  0.0  0.0  0.0 -1.0 -1.0
    ]

b = [-3.0, -8.0, -10.0, -8.0, -14.0, -5.0]

lb = zeros(6)
ub = Inf*ones(6)

ctype = repeat("U",6)
vartype = repeat("I",6)

x, z = HRW.glpk(f, A, b, lb, ub, ctype, vartype)

@printf("Solution: %4.2f\n", -z)
```

## Documentation ODE-Solvers

`solve_dgl` solves systems of explicit ordinary differential equations y'(t) = f(t,y(t)), y(t<sub>0</sub>) = y<sub>0</sub> within the interval [t<sub>0</sub>,T]. Here are y: R&rarr;R<sup>m</sup> and f: R&times;R<sup>m</sup>&rarr;R<sup>m</sup>.

```
t, y = solve_dgl(methode, fh, y0, tspan, h)
```

Input arguments are of type `methode::String`, `fh::Function`, `y0::Vector{Float64}`, `tspan::Vector{Float64}`, `h::Float64`.

The parameters are as follows:
<ul>
  <li> `methode` detemines different time stepping schemes:
    <ul>
      <li>`expl_Euler` first order explicit Euler scheme,</li>
      <li>`verb_Euler` second order (improved) Euler scheme,</li>
      <li>`RK4` classical Runge-Kutta scheme of order 4,</li>
      <li>`impl_Euler` first order implicit Euler scheme (used for stiff problems),</li>
    </ul>
  </li>
<li> `fh` function handle to right hand side `fval = f(t::Float64,y::Vector{Float64})::Vector{Float64}`,</li>
<li> `tspan` determines time range, `tspan[1]=t0`, `tspan[2]=T` final time,</li>
<li> `h` constant time stepping.</li>
</ul>


Return values are:
- `t::Vector{Float64}` vector of time steps,
- `z::Vector{Float64}` value of solution y at each time step.


#### Example 2.1

Let y be the solution of

y'(t) = y(t) / (t-0.5)<sup>0.5</sup>, y(1) = 1.

To compute an approximation of y(9), the  following code can be used:

```
import HRW

# right hand side of ODE
f(t::Float64,y::Vector{Float64}) = [y[1]/(sqrt(t-0.5))]

t, z = HRW.solve_dgl("RK4", f, [1.0], [1.0, 9.0], 0.01)

@printf("Solution: %4.2f\n", z[end])

```

#### Example 2.2

The following example solves the 1 dimensional heat equation

T<sub>t</sub>(x,t) - T<sub>xx</sub>(x,t) = 0, t&ge;0, x&isin; (0,1).

Initial conditions are T(x,0)=sin(&pi;x), boundary conditions are T(0,t)=T(1,t)=0.

The solution is computed using the Line-Method, that transforms the PDE using finite differences in space into a system of ODEs, see [DR] chapter 11.1.

```
function f(t::Float64, y::Vector{Float64})

const N = round(Int64,1.0/hx) - 1

if N==1
	return -1.0/hx^2*2*y[1]
end

if N==2
	return -1.0/hx^2*[2 -1; -1 2]*y
end

result = zeros(N)

result[1] = 2*y[1]-y[2]

for k=2:N-1
   result[k] = -y[k-1]+2*y[k]-y[k+1]
end

result[N] = -y[N-1]+2*y[N]

result = -result/hx^2

end



# main

const hx = 0.2
const ht = 0.001
const t0 = 0.0
const T  = 2.0

# initial condition
const y0 = sin(pi * collect(hx:hx:1-hx))

const tspan = [t0, T]

import HRW

method = "impl_Euler" # expl_Euler,...
t,z = HRW.solve_dgl(method, f, y0, tspan, ht)


# copy vectors to plot solution with correct
# boundary values in space
k = length(z)
l = length(z[1])

Z = zeros(k,l+2)

for i=1:k
   for j=1:l
      Z[i,j+1] = z[i][j]
   end
end

using PyPlot


X = linspace(0, 1, l+2)'
Y = t

surf = plot_surface(X,Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=false, cmap="coolwarm")
xlabel("x")
ylabel("t")
zlabel("T")
title("Solution of Heat Equation")


```


## Project Status

The package is tested against Julia `0.5` on Windows.


## References

[DR] W. Dahmen and A. Reusken: "Numerik für Ingenieure und Naturwissenschaftler", Springer, 2008.
