"""
# Modified from https://github.com/Algopaul/PortHamiltonianBenchmarkSystems.jl.git
# https://www.sciencedirect.com/science/article/abs/pii/S0005109812002257
Structure-preserving tangential interpolation for model reduction of port-Hamiltonian systems
# Arguments
- `n_cells`: The number of masses. The system dimension is `2n_cells`
- `io_dim`: The input and output dimension of the system
- `c`: The amount of damping
- `m`: The weight of the masses
- `k`: The stiffness of the springs
# Outputs
- `config`: The configuration struct for the system. The system can subsequently be created with `construct_system(config)`
"""


struct MSDconfig
    n_cells::Int
    io_dim::Int
    c::Vector{Float64}
    m::Vector{Float64}
    k::Vector{Float64}
    function MSDconfig(n_cells::Int=10, io_dim::Int=2, c::Vector{Float64}=ones(10), m::Vector{Float64}=4ones(10), k::Vector{Float64}=4ones(10)) 
        @assert n_cells > 0        "number of cells must be positive"
        @assert io_dim  > 0        "number of inputs and outputs must be positive"
        @assert io_dim  <= n_cells "number of inputs and outputs must be less than or equal to the number of cells"
        return new(n_cells, io_dim, c, m, k)
    end
end


function construct_system(config::MSDconfig)
    n_cells = config.n_cells
    io_dim  = config.io_dim 
    c = config.c
    m = config.m
    k = config.k
    n = 2*n_cells
    # B is initialized as dense matrix. Since all results of transfer function
    # computations will lead to dense results.
    B = zeros(n, io_dim)
    [B[2 * i, i] = 1.0 for i = 1:io_dim]
    J = zeros(n, n)
    [J[i, i + 1] = 1.0 for i = 1:2:(n - 1)]
    J = J - J'
    # Set constants.
    R = msd_construct_R(n_cells, c)
    Q = msd_construct_Q(n_cells, k, m)
    return (J = J, R = R, Q = Q, B = B)
end

function msd_construct_R(n_cells, c)
    n = 2n_cells
    R = zeros(n, n)
    for (j, i) in enumerate(2:2:n)
        R[i, i] = c[j]
    end
    return R
end

function msd_construct_Q(n_cells, k, m)
    n = 2n_cells
    Q = zeros(n, n)
    n = size(Q, 1)
    for (j, i) in enumerate(1:2:(n - 3))
        msd_construct_Q_add_k_stencil(Q, k[j], i)
    end
    Q[end - 1, end - 1] += k[n_cells]
    for (j, i) in enumerate(2:2:n)
        Q[i, i] = 1 / m[j]
    end
    return Q
end

function msd_construct_Q_add_k_stencil(Q, k, i)
    @views Q[i:(i + 2), i:(i + 2)] .+= [k 0 -k; 0 0 0; -k 0 k]
end


# ẋ = (J-R)Qx + Gu
# y = GᵀQx
function PHSystem(dz, z, obj::MSDconfig, t)
    n = obj.n_cells

    J, R, Q, G = construct_system(obj)
    n, m = size(G)
    E = I(n)
    P = zeros(n, m)
    S = zeros(m, m)
    N = zeros(m, m)


    dz[1:n] = (J-R)*Q*z[1:n] + G*u
    return E, J, R, Q, G, P, S, N
end

"""
Setup in [Gugercin2012](https://doi.org/10.1016/j.automatica.2012.05.052)
# Outputs
"""
function gugercin_pH_msd_chain()
    n_cells = 50
    io_dim  = 2
    c = 1.0
    m = 4.0
    k = 4.0
    return MSDconfig(n_cells, io_dim, c*ones(n_cells), m*ones(n_cells), k*ones(n_cells))
end

"""
generate_MSD_plant(n_cells::Int)

Constructor for the configuration as control problem.
# Arguments
- `n_cells`: The number of masses. The system dimension is `2n_cells`.
# Outputs
- `A`: The system matrix of the plant.
- `B`: The input matrix of the plant, as a concatenation of `B=hcat(B1, B2)`, where `B1` and `B2` describe the disturbance input and the control input, respectively. `B2` has `nw` columns.
- `C`: The output matrix of the plant, as a concatenation of `C=vcat(C1, C2)`, where `C1` and `C2` describe the performance output and the measured output, respectively. `C2` has `nz` columns.
- `D`: The feedthrough matrix of the plant, as a concatenation of `D=vcat(hcat(D11, D12), hcat(D21, D22))`.
- dimension `nw`: dimension of the disturbance input.
- dimension `nz`: dimension of the performance output.
"""
function generate_MSD_plant(n_cells)
    config = MSDconfig(n_cells, 2, ones(n_cells), 4ones(n_cells), 4ones(n_cells))
    J, R, Q, B = construct_system(config)
    # Add disturbance inputs.
    Bw = zero(B)
    Bw[6, 1] = 1.0
    Bw[8, 2] = 1.0
    A = (J - R) * Q
    B = hcat(Bw, B)
    C = B' * Q
    D = zeros(size(B, 2), size(B, 2))
    nz, nw = 2, 2
    return (A, B, C, D, nz, nw)
end