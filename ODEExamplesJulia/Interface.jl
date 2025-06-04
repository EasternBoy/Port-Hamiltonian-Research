push!(LOAD_PATH, ".")
import Pkg
Pkg.activate(".")

using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, Plots

# abstract type myAbstract end

# include("ODEExamplesJulia/DoublePendulum.jl")
# include("ODEExamplesJulia/MassSpringDamper.jl")
# include("ODEExamplesJulia/MagBall.jl")
# include("ODEExamplesJulia/CartPole.jl")
include("RLCLadderNetworks.jl")
# include("ODEExamplesJulia/TodaLattice.jl")
# include("ODEExamplesJulia/ComposeSystem.jl")




n_cells = 2
R  = ones(n_cells)
G  = 10e-6*ones(n_cells)
L  = 2e-6*ones(n_cells)
C₀ = 1e-6*ones(n_cells)
V₀ = ones(n_cells)
initialState = [C₀; zeros(n_cells)]


RLClad = RLC_LadNet(n_cells, initialState, R, L, C₀, G, V₀)
sampling = 0.:0.2e-6:1e-5
saved_values = update!(RLClad, sampling)


ns = length(saved_values.t)
y = [saved_values.saveval[i][2][1] for i in 1:ns]
u = [saved_values.saveval[i][3][1] for i in 1:ns]


plot()
plot!(saved_values.t, y)
plot!(saved_values.t, u)