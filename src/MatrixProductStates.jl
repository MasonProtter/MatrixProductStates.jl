# [[file:~/.julia/dev/MatrixProductStates/README.org::*Source%20Code][Source Code:1]]
module MatrixProductStates

using LinearAlgebra, TensorOperations, TensorCast, LowRankApprox

export *, /, ==, ≈, isequal, adjoint, getindex, randn
export MPS, MPO, Left, Right, compress

include("MPS.jl")
include("MPO.jl")
include("contraction.jl")
include("timeevolution.jl")
include("compression.jl")
# #---------------------------------------------------------------------
# # DMRG Algorithm

# function DMRG(ψ::MPS, H::MPO)
#     ψ = ψguess
#     converged = false
    
#     while !(converged)
#         right_sweep!(ψ)
#         left_sweep!(ψ)
#         converged = check_convergance(ψ)
#     end
# end

end
# Source Code:1 ends here
