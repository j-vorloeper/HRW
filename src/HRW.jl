# Autor:        Juergen Vorloeper
# Organization: Hochschule Ruhr West
# Date:         Summer 2017

module HRW


export glpk
export solve_dgl
export GnuPlotOutput_Direct, GnuPlotOutput, GnuPlotOutput_Close

# load files
include("HRW_glpk.jl")
include("HRW_ode.jl")
include("HRW_gnuplot.jl")


# set up global variables, that stores gnuplot's state
gnuplot_state = GnuplotState(false,0,Any[])


# Register function "GnuPlotOutput_Close" to be called,
# when gnuplot_state goes out of scope
finalizer(gnuplot_state,GnuPlotOutput_Close)

end
