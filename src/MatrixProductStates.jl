# Module Definition
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Module%20Definition][Module Definition:1]]
module MatrixProductStates

using LinearAlgebra, TensorOperations, TensorCast, LowRankApprox, Arpack, Strided

export *, /, ==, ≈, isequal, adjoint, getindex, randn
export MPS, MPO, left, right, compress, imag_time_evolution, rightcanonical, leftcanonical 
export ground_state

include("MPS.jl")
include("MPO.jl")
include("compression.jl")
include("contraction.jl")
include("timeevolution.jl")
include("groundstate.jl")
include("utils.jl")

end
# Module Definition:1 ends here
