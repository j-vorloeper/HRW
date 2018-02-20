# Author:
# Jürgen Vorloeper, Hochschule Ruhr West, 2018

# code partially based on Julia package "Gaston"
# see: https://github.com/mbaz/Gaston.jl


# verify gnuplot is present on this system
if !success(`gnuplot --version`)
    error("*** Error: gnuplot is not available on this system.")
end



type GnuplotState

running::Bool  # true when gnuplot is already running
fid
figs::Array  # list of identifiers to all figures

function GnuplotState(running, fid, figs)

  new(running,fid,figs)

end

end



function GnuPlotOutput_Init()

global gnuplot_state

pin = Base.Pipe()
r = spawn(`gnuplot`, pin)

gnuplot_state.running = true
gnuplot_state.fid = (pin, r)

w = write(pin, string("set terminal wxt","\n"))
if !(w > 0)
  println("Something went wrong writing to gnuplot STDIN.")
  return
end


return (pin, r)

end



# close gnuplot pipe
function GnuPlotOutput_Close(x...)
    global gnuplot_state

    if gnuplot_state.running
        close(gnuplot_state.fid[1])
    end

    # reset gnuplot_state
    gnuplot_state.running = false
    gnuplot_state.figs = Any[]

    return 0
end



function GnuPlotOutput_Direct(x::Vector{Float64}, y::Vector{Float64})


global gnuplot_state

if gnuplot_state.running == false
  GnuPlotOutput_Init()
end

@assert(length(x) == length(y), "Vectors x and y must have the same number of elements")

s = string("set xlabel \"x\" ");
write(gnuplot_state.fid[1], string(s,"\n"))

s = string("set ylabel \"y\" ");
write(gnuplot_state.fid[1], string(s,"\n"))


s = string("plot [",x[1],":",x[end],"] '-' notitle with lines")
write(gnuplot_state.fid[1], string(s,"\n"))

for i=1:length(x)
  s = string(i, "  ", y[i])
  write(gnuplot_state.fid[1], string(s,"\n"))
end

write(gnuplot_state.fid[1],string("e","\n"))
flush(gnuplot_state.fid[1])


end



function GnuPlotOutput(x::Vector{Float64}, y::Vector{Float64}, fname::String, title::String)

if length(x)!=length(y)
			println("*** Error")
			return -1
end


fid = open(fname*".gnu","w")

if fid==0
  error("*** GnuPlotOutput: Unable to open file")
end


@printf(fid, "set title \"%s\" \n", title)
@printf(fid, "set xlabel \"x\"\n")
@printf(fid, "set ylabel \"y\"\n")
@printf(fid, "plot \"%s.dat\" notitle with lines\n", fname)
@printf(fid, "reset\n")
close(fid)


fid = fopen(fname*".dat", "w")
if fid==0
    error("*** GnuPlotOutput: Unable to open file")
end


for k=1:length(x)
			@printf(fid, "%10.8e %10.8e\n", x[k], y[k])
end

close(fid)

println("GnuPlotOutput: output in files '%s.dat + .gnu'\n", fname)
end
