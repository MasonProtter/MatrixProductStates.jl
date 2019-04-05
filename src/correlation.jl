# Correlation Functions
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Correlation%20Functions][Correlation Functions:1]]
"""
    two_point_correlator(i, j, L, op_i::Matrix, op_j::Matrix)

Create an MPO on `L` sites (with bond dimension 1) representing identity operators everywhere except
sites `i` and `j` where `op_i` and `op_j` are inserted instead. ie.

    𝟙 ⊗ 𝟙 ⊗ ... ⊗ op_i ⊗ 𝟙 ⊗ ... ⊗ op_j ⊗ 𝟙 ⊗ ... ⊗ 𝟙
"""
function two_point_correlator(i, j, L, op_i::Matrix, op_j::Matrix)
    d = size(op_i)[1]
    @assert (size(op_i) == (d, d)) && (size(op_j) == (d, d))
    @assert i in 1:L
    @assert j in 1:L

    op_i_tnsr = resphape(convert(Matrix{Complex{Float64}}, op_i), 1, 1, d, d)
    op_j_tnsr = resphape(convert(Matrix{Complex{Float64}}, op_j), 1, 1, d, d)

    id_tnsr   = reshape(typeof(op_i)(one(eltype(op_i_tnsr))I, d, d), 1, 1, d, d)

    tensors = map(1:L) do l
        O_tnsr = (l == i ? op_i_tnsr : 
                  l == j ? op_j_tnsr : 
                  id_tnsr)
    end 
    MPO{L,Complex{Float64}}(tensors)
end
# Correlation Functions:1 ends here
