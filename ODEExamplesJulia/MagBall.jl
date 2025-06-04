"""
# Thomas Beckers, Gaussian Process Port-Hamiltonian Systems: Bayesian Learning with Physics Prior CDC2022
# https://arxiv.org/pdf/2305.09017
# Arguments
- 'R': Resistant
- 'm': Mass of ball
- 'c': Drag parameter
- 'u': External input
"""

mutable struct MagBall
    s::Vector{Float64}
    R::Float64
    m::Float64
    c::Float64
    u::Float64

    function MagBall(s::Vector{Float64}, R::Float64, m::Float64, c::Float64, u::Float64)
        return new(s, R, m, c, u)
    end
end

function Hamiltonian!(ds, s, mb, t)
    R = mb.R
    m = mb.m
    c = mb.c
    u = mb.u

    J = [0  1  0;
        -1  0  0;
         0  0  0]
    R = [0 0           0;
         0 c*abs(s[2]) 0;
         0 0           1/R]
    G = [0; 0; 1]

    Lₛ = 1/(0.1 + s[1]^2)

    H   =  1/(2m)*s[2]^2 + 1/(2Lₛ)*s[3]^2
    ∇ₛH  = [s[1]*s[3]^2; 1/m*s[2]; 1/Lₛ*s[3]]
    ds[1:3]  = (J-R)*∇ₛH + G*u
    return s, ds, ∇ₛH, H, J, R, G
end

function PHSystem(obj::MagBall)
    return Hamiltonian!(zeros(4), obj.s, obj, 0.)
end

function update!(obj::MagBall, tspan::Tuple{Float64, Float64})
    s₀    = obj.s
    prob  = ODEProblem(Hamiltonian!, s₀, tspan, obj)
    sol   = solve(prob,Tsit5())
    obj.s = sol.u[end]
    return sol
end
