"""
# ***********************************************************************
# https://diamweb.ewi.tudelft.nl/~jeltsema/Book/PHSystemsTheory_printedbook.pdf
# https://underactuated.mit.edu/acrobot.html
# ***********************************************************************

# H = 1/2q̇ᵀM(q)q̇ + U(q)
# q̇ = ∂H/∂p (p = M(q)q̇)
# ṗ = -∂H/∂q + Bu = M(q)q̈ + Ṁ(q)q̇ = -C(q,q̇)q̇ + τᵧ(q) + Bu + Ṁ(q)q̇

# Arguments
- `s`: state variable of
- `m₁`: The input and output dimension of the system
- `m₂`: The amount of damping
- `ℓ`: Length of pendulum
- `f`: External force
"""

mutable struct CartPole
    s::Vector{Float64}
    m₁::Float64    
    m₂::Float64
    ℓ::Float64
    f::Float64

    function CartPole(init::Vector{Float64}, m₁::Float64, m₂::Float64, ℓ::Float64, f::Float64)
        return new(init, m₁, m₂, ℓ, f)
    end
end

function EulerLagrange!(dz, z, obj, t)
    g  = 9.8
    x  = z[1]
    θ  = z[2]
    fₓ = obj.f
    m₁ = obj.m₁
    m₂ = obj.m₂
    ℓ  = obj.ℓ

    dz[1] = z[3] 
    dz[2] = z[4]

    dz[3] = 1/(m₁+m₂*sin(θ)^2)*(fₓ + m₂*sin(θ)*(ℓ*dθ^2 + g*cos(θ)))
    dz[4] = 1/(ℓ*(m₁+m₂*sin(θ)^2))*(-fₓ*cos(θ) - m₂*ℓ*dθ^2*cos(θ)*sin(θ)-(m₁+m₂)*g*sin(θ))
end


function PHSystem(obj::CartPole)
    g  = 9.8
    x  = obj.s[1]
    θ  = obj.s[2]
    ẋ  = obj.s[3]
    dθ = obj.s[4]
    fₓ = obj.f

    m₁ = obj.m₁
    m₂ = obj.m₂
    ℓ  = obj.ℓ
    M  = [m₁+m₂       m₂*ℓ*cos(θ);
          m₂*ℓ*cos(θ) m₂*ℓ^2]
    
    q = [x;  θ]
    q̇ = [ẋ; dθ]
    p = M*q̇
    z = [q; p]

    B = [1;0]


    K = 1/2*dot(q̇,M,q̇)
    P = m₂*g*ℓ*cos(θ)
    H = K + P

    ∂M∂x = [0 0;
            0 0]
    ∂M∂θ = [0          -m₂*ℓ*sin(θ);
           -m₂*ℓ*sin(θ) 0]
    ∂H∂q = 1/2*[0; dot(q̇, ∂M∂θ, q̇)] + [0; -m₂*ℓ*sin(θ)]

    ṗ   = -∂H∂q + B*fₓ
    ż   = [q̇; ṗ]
    ∂H∂z = [inv(M)q̇ ;-∂H∂q]

    J = [zeros(2,2) I(2); 
        -I(2)       zeros(2,2)]
    R = zeros(4,4)
    G = [0; 0; 1; 0]

    return z, ż, ∂H∂z, H, J, R, G
end

function update!(obj::CartPole, tspan::Tuple{Float64, Float64})
    z0    = obj.s
    prob  = ODEProblem(EulerLagrange!, z0, tspan, obj)
    sol   = solve(prob,Tsit5())
    obj.s = sol.u[end]
    return sol
end
