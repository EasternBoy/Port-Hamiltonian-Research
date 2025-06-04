# Yuwei Geng, Data-Driven Reduced-Order Models for Port-Hamiltonian Systems with Operator Inference
# https://arxiv.org/pdf/2501.02183
# Toda lattice


push!(LOAD_PATH, ".")
import Pkg
Pkg.activate(".")

using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, Plots, NPZ


using LinearAlgebra
using SparseArrays




mutable struct TodaLattice
    n_cells::Int64
    s::Vector{Float64}
    γ::Vector{Float64}
    function TodaLattice(n_cells::Int64, s::Vector{Float64}, γ::Vector{Float64})
        return new(n_cells, s, γ)
    end
end


# mode    = 1 #Training
mode    = 2 #Testing
global pre_u = 0
function input(t::Float64; mode = mode)
    global pre_u
    if mode == 1
        if t%1 == 0
            pre_u = 0.5(rand() - rand())
        end
        return pre_u
    else
        res = 0
        if t%200 <= 5
            res = 0.1
        else
            res = 0
        end
        return res
    end
end

if mode == 1
    sampling     = 0.:0.1:100
else
    sampling     = 0.:0.1:100
end


function PHSystem(ds, s, config, t)
    N = config.n_cells
    u = input(t)
    γ = config.γ

 

    Jₚ = [spzeros(N,N) sparse(I, N, N); -sparse(I, N, N) spzeros(N,N)]
    Rₚ = [spzeros(N,N) spzeros(N,N); spzeros(N,N) spdiagm(γ)]
    Gₚ = [spzeros(N)' 1. spzeros(N-1)']'

    ∇H = spzeros(2N)
    for i in 1:N
        if i == 1
            ∇H[i]   = exp(s[i] - s[i+1]) - 1
        elseif i==N
            ∇H[i]   = exp(s[i]) - exp(s[i-1] - s[i])
        else
            ∇H[i]   = exp(s[i] - s[i+1]) - exp(s[i-1] - s[i])
        end
        ∇H[i+N] = s[i+N]
    end

    ds[1:2N] = Vector((Jₚ - Rₚ)*∇H + Gₚ*u)
end


function ODEcallback(s, t, integrator)
    config = integrator.p
    N = config.n_cells
    u = input(t)
    γ = config.γ

 
    Jₚ = [spzeros(N,N) sparse(I, N, N); -sparse(I, N, N) spzeros(N,N)]
    Rₚ = [spzeros(N,N) spzeros(N,N); spzeros(N,N) spdiagm(γ)]
    Gₚ = [spzeros(N)' 1. spzeros(N-1)']'


    H  = 0
    ∇H = spzeros(2N)
    for i in 1:N
        if i == 1
            ∇H[i]   = exp(s[i] - s[i+1]) - 1
            H  += (1/2)*s[i+N]^2 + exp(s[i] - s[i+1])
        elseif i == N
            ∇H[i]   = exp(s[i]) - exp(s[i-1] - s[i])
        else
            ∇H[i]   = exp(s[i] - s[i+1]) - exp(s[i-1] - s[i])
            H  += (1/2)*s[i+N]^2 + exp(s[i] - s[i+1])
        end
        ∇H[i+N] = s[i+N]
    end

    H += exp(s[N]) - s[1] - N
    ds = Vector((Jₚ - Rₚ)*∇H + Gₚ*u)
    y  = Gₚ'*∇H
    return s[1:end], ds, H, ∇H, y, u
end


function update!(obj::TodaLattice, sampling)
    s0           = obj.s
    prob         = ODEProblem(PHSystem, s0, (sampling[1], sampling[end]), obj);
    saved_values = SavedValues(Float64, Tuple{Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Float64, Float64});
    cb           = SavingCallback(ODEcallback, saved_values, saveat = sampling);
    sol          = solve(prob, Tsit5(), callback = cb, saveat = sampling);
    obj.s        = sol.u[end];
    return saved_values, sol
end


n_cells = 50
γ       = 0.2*ones(n_cells)
initial = zeros(2n_cells)

toda         = TodaLattice(n_cells, initial, γ)

saved_values, sol = update!(toda, sampling)



ns    = length(saved_values.t)

Xs    = reduce(hcat, [saved_values.saveval[i][1] for i in 1:ns])
Xdots = reduce(hcat, [saved_values.saveval[i][2] for i in 1:ns])
Hs    = [saved_values.saveval[i][3] for i in 1:ns]
∇Hs   = reduce(hcat, [saved_values.saveval[i][4] for i in 1:ns])
y     = [saved_values.saveval[i][5] for i in 1:ns]
u     = [saved_values.saveval[i][6] for i in 1:ns]


#Select basic vectors
#Reduce order r
r = 10
V, Λ, U = svd(Xs)
V = V[:,1:r]


if mode == 1
    npzwrite(os.path.join("..","data", "TodaLat_data_train.npz"), Dict(
        "Xs"     => Xs,
        "Xdots"  => Xdots,
        "Ys"     => y,
        "Hs"     => Hs,
        "gradHs" => ∇Hs,
        "Us"     => u,
        "Ts"     => saved_values.t,
        "Vtrans" => V))
else
    npzwrite(os.path.join("..","data", "TodaLat_data_test.npz"), Dict(
        "Xs"     => Xs,
        "Xdots"  => Xdots,
        "Ys"     => y,
        "Hs"     => Hs,
        "gradHs" => ∇Hs,
        "Us"     => u,
        "Ts"     => saved_values.t,
        "Vtrans" => V))
end


plot()
plot!(sampling, y, label="y")
plot!(sampling, u, label="u")
