# Matrix Product Operators
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Matrix%20Product%20Operators][Matrix Product Operators:1]]
"""
    MPO{L, T<:Number}

Matrix product operator on L sites. The `i`th tensor in the operator
has indices `[aⁱ⁻¹, aⁱ, σⁱ, σ′ⁱ]` where `(σⁱ, σ′ⁱ)` are the physical
indices and `(aⁱ⁻¹, aⁱ)` are bond indices.

A four site MPS would be diagrammatically represented

    σ¹          σ²          σ³          σ⁴
    |           |           |           | 
    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--•     
    |           |           |           | 
    σ′¹         σ′²         σ′³         σ′⁴


Note that `a⁰` and `aᴸ` must be of dimension 1.
"""
struct MPO{L, T<:Number}
    tensors::Vector{Array{T,4}}
end


"""
    MPO(W::Array{T,4}, L)
Create an `MPO` for `L` sites with all interior sites containing the tensor
`W`. The tensor is assumed to have the usual matrix-of-operators structure,
with the first two indices being the bond (matrix) dimension and the last two
indices being the physical (operator) dimension. The first and last sites only
use the last row and first column of `W`, respectively.

For example, the MPO form of the Hamiltonian for the TFIM is
constructed as with coupling `g` and length `L` is constructed as
follows:

    id = [1 0
          0 1]

    σᶻ = [1  0 
          0 -1]

    σˣ = [0 1
          1 0]

    σʸ = [0  -im
          im   0]

    W = zeros(3, 3, 2, 2)
    W[1, 1, :, :] = id
    W[2, 1, :, :] = σᶻ
    W[3, 1, :, :] = -g*σˣ
    W[3, 2, :, :] = -σᶻ
    W[3, 3, :, :] = id

returning 
 
    Ĥ::MPO = Ŵ¹ Ŵ² Ŵ³ ⋅⋅⋅ Ŵᴸ⁻¹ Wᴸ
"""
function MPO(W::Array{T,4}, L) where {T}
    L >= 2 || throw(DomainError(L, "At least 2 sites."))

    tensors = Vector{Array{T,4}}(undef, L)
    
    tensors[1] = W[end:end, :, :, :] # Row vector.
    for i in 2:(L-1)
        tensors[i] = W # Matrix
    end
    tensors[L] = W[:, 1:1, :, :] # Column vector.

    MPO{L,T}(tensors)
end

Base.:(==)(O::MPO, U::MPO) = O.tensors == U.tensors
Base.:(≈)(O::MPO, U::MPO)  = O.tensors ≈ U.tensors
Base.getindex(O::MPO, args...) = getindex(O.tensors, args...)

function Base.show(io::IO, ::MIME"text/plain", O::MPO{L, T}) where {L, T}
    d = length(O[2][1, 1, 1, :])
    bonddims = [size(O[i][:, :, 1, 1]) for i in 1:L]
    println(io, "Matrix product Operator on $L sites")
    _show_mpo_dims(io, L, d, bonddims)
end

function _show_mpo_dims(io::IO, L, d, bonddims)
    println(io, "  Physical dimension: $d")
    print(io, "  Bond dimensions:   ")
    if L > 8
        for i in 1:8
            print(io, bonddims[i], " × ")
        end
        print(io, " ... × ", bonddims[L])
    else
        for i in 1:(L-1)
            print(io, bonddims[i], " × ")
        end
        print(io, bonddims[L])
    end
end

function Base.show(io::IO, O::MPO{L, T}) where {L, T}
    print(io, "MPO on $L sites")
end
# Matrix Product Operators:1 ends here
