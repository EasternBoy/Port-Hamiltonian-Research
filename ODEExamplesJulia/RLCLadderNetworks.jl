# ***********************************************************************
# Chaturantabut et al. - 2016 - Structure-Preserving Model Reduction for Nonlinear Port-Hamiltonian 
# https://epubs.siam.org/doi/pdf/10.1137/15M1055085?casa_token=SucwcwUXkrkAAAAA:gv7hL3-epFoZWJ_eyJVFDd_DYQQzFm90HEp5Ux1FRySA41XtfuJT6D5DIkidr5kYzJJpW539z6of
# ***********************************************************************

# state s
# resistance R
# inductance L
# capacitance C
# conductance G
# initial voltage V₀
# input u
push!(LOAD_PATH, ".")
import Pkg
Pkg.activate(".")

using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, Plots, NPZ


using LinearAlgebra
using SparseArrays

mutable struct RLC_LadNet
    n_cells::Int64
    s::Vector{Float64}
    R::Vector{Float64}
    L::Vector{Float64}
    C₀::Vector{Float64}
    G::Vector{Float64}
    V₀::Vector{Float64}

    function RLC_LadNet(n_cells::Int64, s::Vector{Float64}, R::Vector{Float64}, L::Vector{Float64}, C₀::Vector{Float64}, G::Vector{Float64}, V₀::Vector{Float64})
        obj = new(n_cells, s, R, L, C₀, G, V₀)
        return obj
    end
end

function input(s::Vector{Float64}, t::Float64)
    σ = 0.5e-6 
    m = 3e-6
    #Gaussian Windows
    return [3*exp(-(t - m)^2/(2π*σ^2)), 0.]
end


function PHSystem(ds, s, config, t)
    n = config.n_cells
    u = input(s,t)
 
    C₀ = config.C₀
    L  = config.L
    G  = config.G
    R  = config.R
    V₀ = config.V₀

    S = spzeros(n,n)
    [S[i,i]   =  1. for i in 1:n]
    [S[i,i+1] = -1. for i in 1:n-1]
    Jₚ = [spzeros(n,n) S; -S' spzeros(n,n)]
    Rₚ = [spdiagm(G) spzeros(n,n); spzeros(n,n) spdiagm(R)]
    Gₚ = [spzeros(n)' 1 spzeros(n-1)'; spzeros(n-1)' 1 spzeros(n)']'

    ∇H = spzeros(2n)
    for i in 1:n
        ∇H[i]   = V₀[i]*(exp(s[i]/(C₀[i]*V₀[i]))-1)
        ∇H[i+n] = s[i+n]/L[i]
    end

    ds[1:2n] = Vector((Jₚ - Rₚ)*∇H + Gₚ*u)
end

function ODEcallback(s, t, integrator)
    config = integrator.p
    n = config.n_cells
    u = input(s,t)
 
    C₀ = config.C₀
    L  = config.L
    V₀ = config.V₀
    R  = config.R

    S  = spzeros(n,n)
    Jₚ = [spzeros(n,n) S; -S' spzeros(n,n)]
    Rₚ = [spdiagm(G) spzeros(n,n); spzeros(n,n) spdiagm(R)]
    Gₚ = [spzeros(n)' 1 spzeros(n-1)'; spzeros(n-1)' 1 spzeros(n)']'

    H  = 0
    ∇H = spzeros(2n)
    for i in 1:n
        H  = H + C₀[i]*V₀[i]^2*(exp(s[i]/(C₀[i]*V₀[i]))-1) - s[i]*V₀[i] + s[i+n]^2/(2*L[i])
        ∇H[i] = V₀[i]*(exp(s[i]/(C₀[i]*V₀[i]))-1)
        ∇H[i+n] = s[i+n]/L[i]
    end

    ds = Vector((Jₚ - Rₚ)*∇H + Gₚ*u)
    y  = Gₚ'*∇H

    return (s, ds, H, ∇H, y, u)
end


function update!(obj::RLC_LadNet, sampling)
    s₀           = obj.s
    prob         = ODEProblem(PHSystem, s₀, (sampling[1], sampling[end]), obj);
    saved_values = SavedValues(Float64, Tuple{Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}});
    cb           = SavingCallback(ODEcallback, saved_values, saveat = sampling);
    sol          = solve(prob, Tsit5(), callback = cb, abstol = 1e-12, reltol = 1e-12);
    obj.s        = sol.u[end];
    return saved_values
end


n_cells = 10
R  = ones(n_cells)
G  = 10e-6*ones(n_cells)
L  = 2e-6*ones(n_cells)
C₀ = 1e-6*ones(n_cells)
V₀ = ones(n_cells)
initialState = [C₀; zeros(n_cells)]


RLClad       = RLC_LadNet(n_cells, initialState, R, L, C₀, G, V₀)
sampling     = 0.:0.2e-6:1e-5
saved_values = update!(RLClad, sampling)


ns    = length(saved_values.t)

Xs    = reduce(hcat, [saved_values.saveval[i][1] for i in 1:ns])
Xdots = reduce(hcat, [saved_values.saveval[i][2] for i in 1:ns])
Hs    = [saved_values.saveval[i][3] for i in 1:ns]
∇Hs   = reduce(hcat, [saved_values.saveval[i][4] for i in 1:ns])
y     = reduce(hcat, [saved_values.saveval[i][5] for i in 1:ns])
u     = reduce(hcat, [saved_values.saveval[i][6] for i in 1:ns])


#Select basic vectors
#Reduce order r
r = 6
V, Λ, U = svd(Xs)
V = V[:,1:r]


# Establish training data with Gaussian windown input
npzwrite(os.path.join("..","data", "RLC_data_train.npz"), Dict(
    "Xs"     => V'*Xs,
    "Xdots"  => V'*Xdots,
    "Ys"     => y,
    "Hs"     => Hs,
    "gradHs" => ∇Hs,
    "Us"     => u,
    "Ts"     => saved_values.t,
    "Vtrans" => V,
))

plot()
plot!(sampling, y[1,:], label="y")
plot!(sampling, u[1,:], label="u")