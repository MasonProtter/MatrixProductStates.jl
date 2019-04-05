# Compression
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Compression][Compression:1]]
function compress(ψ::MPS{L, T}, to_the::Right; Dcut::Int=typemax(Int)) where {L, T}
    tensors = Array{T, 3}[]
    
    B = ψ[1]
    d = length(B[1, 1, :])
    
    @cast Bm[(σ¹, a⁰), a¹] |= B[a⁰, a¹, σ¹]
    U, S, V = psvd(Bm, rank=Dcut)
    #S = S/√sum(S .^ 2)

    @cast A[a⁰, a¹, σ¹] |= U[(σ¹, a⁰), a¹] (σ¹:d)
    push!(tensors, A)
    
    for i ∈ 2:L
        B = ψ[i]
        d = length(B[1, 1, :])

        @tensor M[aⁱ⁻¹, aⁱ, σⁱ] := (Diagonal(S)*V')[aⁱ⁻¹, aⁱ⁻¹′] * B[aⁱ⁻¹′, aⁱ, σⁱ]
        @cast   Mm[(σⁱ, aⁱ⁻¹), aⁱ] |= M[aⁱ⁻¹, aⁱ, σⁱ]
        
        U, S, V = psvd(Mm, rank=Dcut)
        #S = S/√sum(S .^ 2)

        @cast A[aⁱ⁻¹, aⁱ, σⁱ] |= U[(σⁱ, aⁱ⁻¹), aⁱ] (σⁱ:d)
        push!(tensors, A)
    end
    MPS{L, T}(tensors), Left()
end

leftcanonical(ψ) = compress(ψ, right)[1]

function compress(ψ::MPS{L, T}, to_the::Left; Dcut::Int=typemax(Int)) where {L, T}
    tensors = Array{T, 3}[]
    
    A = ψ[L]
    d = length(A[1, 1, :])
    @cast Am[aᴸ⁻¹, (σᴸ, aᴸ)] |= A[aᴸ⁻¹, aᴸ, σᴸ]
    
    U, S, V = psvd(Am, rank=Dcut)
    #S = S/√sum(S .^ 2)    

    @cast B[aᴸ⁻¹, aᴸ, σᴸ] |= V'[aᴸ⁻¹, (σᴸ, aᴸ)] (σᴸ:d)
    push!(tensors, B)
    
    for i ∈ (L-1):-1:1
        A = ψ[i]
        d = length(A[1, 1, :])
        @tensor M[aⁱ⁻¹, aⁱ, σⁱ]    := A[aⁱ⁻¹, aⁱ′, σⁱ] * (U * Diagonal(S))[aⁱ′, aⁱ]
        @cast   Mm[aⁱ⁻¹, (σⁱ, aⁱ)] |= M[aⁱ⁻¹, aⁱ, σⁱ]
        
        U, S, V = psvd(Mm, rank=Dcut)
        #S = S/√sum(S .^ 2)

        @cast B[aⁱ⁻¹, aⁱ, σⁱ] |= V'[aⁱ⁻¹, (σⁱ, aⁱ)] (σⁱ:d)
        push!(tensors, B)
    end
    MPS{L, T}(reverse(tensors)), Right()
end

rightcanonical(ψ) = compress(ψ, left)[1]

compress(ψ; Dcut) = compress(ψ, left, Dcut=Dcut)[1]
# Compression:1 ends here
