#---------------------------------------------------------------------
# MPS Compression

abstract type Direction end
struct Left  <: Direction end
struct Right <: Direction end

flip(::Left)  = Right()
flip(::Right) = Left()

function partial_svd(A::AbstractMatrix; rank::Integer=typemax(Int))
    U, S, V = svd(A)
    K = min(length(S), rank)
    @inbounds U[:, 1:K], S[1:K], V[:, 1:K]
end

function compress(ψ::MPS{L, T}, to_the::Left; Dcut=100) where {L, T}
    tensors = Array{T, 3}[]
    
    A = ψ[L]
    Dᴸ⁻¹, Dᴸ, d = size(A)
    @cast Am[aᴸ⁻¹, (σᴸ, aᴸ)] |= A[aᴸ⁻¹, aᴸ, σᴸ]
    
    U, S, V = partial_svd(Am, rank=Dcut)
    
    @cast B[aᴸ⁻¹, aᴸ, σᴸ] |= V'[aᴸ⁻¹, (σᴸ, aᴸ)] (σᴸ:d)
    push!(tensors, B)
    
    for i ∈ (L-1):-1:1
        A = ψ[i]
        Dⁱ⁻¹, Dⁱ, d = size(A)
        @tensor M[aⁱ⁻¹, aⁱ, σⁱ]    := A[aⁱ⁻¹, aⁱ′, σⁱ] * (U * Diagonal(S))[aⁱ′, aⁱ]
        @cast   Mm[aⁱ⁻¹, (σⁱ, aⁱ)] |= M[aⁱ⁻¹, aⁱ, σⁱ]
        
        U, S, V = partial_svd(Mm, rank=Dcut)
        
        @cast B[aⁱ⁻¹, aⁱ, σⁱ] |= V'[aⁱ⁻¹, (σⁱ, aⁱ)] (σⁱ:d)
        push!(tensors, B)
    end
    # @assert (length(U) == 1) && (length(S) == 1)
    # tensors[1] *= U[1]*S[1]
    MPS{L, T}(reverse(tensors)), Right()
end

# function compress(ψ::MPS{L, T}, Dcut=100, to_the::Right) where {L, T}
#     tensors = Array{T, 3}[]
    
#     B = ψ[L]
#     D⁰, D¹, d = size(B)
#     @cast Bm[(σ¹, a⁰), a¹] := B[a⁰, a¹, σ¹]
#     U, S, V = partial_svd(Bm, Dcut)
#     A = reshape(U, (1, length(S), d))
#     push!(tensors, A)
    
#     for i ∈ (L-1):-1:2
#         A = ψ[i]
#         Dⁱ⁻¹, Dⁱ, d = size(A)
#         @tensor M[aⁱ⁻¹, aⁱ, σⁱ]    := A[aⁱ⁻¹, aⁱ′, σⁱ] * (U * Diagonal(S))[aⁱ′, aⁱ]
#         @cast   Mm[aⁱ⁻¹, (aⁱ, σⁱ)] := M[aⁱ⁻¹, aⁱ, σⁱ]
#         U, S, V = partial_svd(Mm, rank=Dcut)
#         B = reshape(V', (length(S), Dⁱ, d))
#         push!(tensors, B)
#     end
#     MPS{L, T}(tensors), Left()
# end
