

# #+HTML: </p>
# #+HTML: </details>

# #+HTML: <details><summary>MPS Contraction</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Matrix%20Product%20States][Matrix Product States:3]]
"""
    Base.:(*)(ψ′::Adjoint{T, MPS{L, T}}, ϕ::MPS{L, T}) where {L, T}
representing
    •--(b¹ b¹)--•--(b² b²)--•--(b³ b³)--•       
    |           |           |           | 
    σ′¹         σ′²         σ′³         σ′⁴
    σ′¹         σ′²         σ′³         σ′⁴
    |           |           |           | 
    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--•
"""
function Base.:(*)(ψ′::Adjoint{T, MPS{L, T}}, ϕ::MPS{L, T}) where {L, T}
    ψ = ψ′.parent

    M   = ϕ.tensors[1]
    M̃dg = dg(ψ.tensors[1])
    
    @tensor cont[b₁, a₁] := M̃dg[b₁, 1, σ₁] * M[1, a₁, σ₁]
    
    for i in 2:L-1
        M   = ϕ.tensors[i]
        M̃dg = dg(ψ.tensors[i])

        @tensor cont[bᵢ, aᵢ] := M̃dg[bᵢ, bᵢ₋₁, σᵢ] * cont[bᵢ₋₁, aᵢ₋₁] * M[aᵢ₋₁, aᵢ, σᵢ]
    end
    M   = ϕ.tensors[L]
    M̃dg = dg(ψ.tensors[L])
    
    @tensor M̃dg[1, bᴸ⁻¹, σᴸ] * cont[bᴸ⁻¹, aᴸ⁻¹] * M[aᴸ⁻¹, 1, σᴸ]
end
# Matrix Product States:3 ends here

 

# #+HTML: <details><summary>MPO Contraction</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Matrix%20Product%20Operators][Matrix Product Operators:2]]
"""
    Base.:(*)(O::MPO, ψ::MPS)
representing

    σ¹          σ²          σ³          σ⁴
    |           |           |           | 
    •--(b¹ b¹)--•--(b² b²)--•--(b³ b³)--•     
    |           |           |           | 
    σ′¹         σ′²         σ′³         σ′⁴
    σ′¹         σ′²         σ′³         σ′⁴
    |           |           |           | 
    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--•     
"""
function Base.:(*)(O::MPO{L, T}, ψ::MPS{L, T}) where {L, T}
    tensors = Array{T,3}[]
    for i in 1:L
        W = O.tensors[i]
        M = ψ.tensors[i]

        @reduce N[(bᵢ₋₁, aᵢ₋₁), (bᵢ, aᵢ), σᵢ] :=  sum(σ′ᵢ) W[bᵢ₋₁, bᵢ, σᵢ, σ′ᵢ] * M[aᵢ₋₁, aᵢ, σ′ᵢ]
        
        push!(tensors, N)
    end
    MPS{L, T}(tensors)
end


"""
    Base.:(*)(O1::MPO, O2::MPO)
representing

    σ¹          σ²          σ³          σ⁴
    |           |           |           | 
    •--(b¹ b¹)--•--(b² b²)--•--(b³ b³)--•     
    |           |           |           | 
    σ′′¹        σ′′²        σ′′³        σ′′⁴
    σ′′¹        σ′′²        σ′′³        σ′′⁴
    |           |           |           | 
    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--• 
    |           |           |           | 
    σ′¹         σ′²         σ′³         σ′⁴    
"""
function Base.:(*)(O1::MPO{L, T}, O2::MPO{L, T}) where {L, T}
    tensors = Array{T,4}[]
    for i in 1:L
        W1 = O1.tensors[i]
        W2 = O2.tensors[i]

        @reduce V[(bᵢ₋₁, aᵢ₋₁), (bᵢ, aᵢ), σᵢ, σ′ᵢ] :=  sum(σ′′ᵢ) W1[bᵢ₋₁, bᵢ, σᵢ, σ′′ᵢ] * W2[aᵢ₋₁, aᵢ, σ′′ᵢ, σ′ᵢ]
        
        push!(tensors, V)
    end
    MPO{L, T}(tensors)
end

"""
    Base.:(*)(ψ::Adjoint{T,MPS{L,T}}, O::MPO) where {L,T}
representing

    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--•       
    |           |           |           | 
    σ′¹         σ′²         σ′³         σ′⁴
    σ′¹         σ′²         σ′³         σ′⁴
    |           |           |           | 
    •--(b¹ b¹)--•--(b² b²)--•--(b³ b³)--•
    |           |           |           | 
    σ¹          σ²          σ³          σ⁴ 
"""
function Base.:(*)(ψ′::Adjoint{T,MPS{L,T}}, O::MPO{L, T}) where {L,T}
    ψ = ψ′.parent
    tensors = Array{T,3}[]
    Ws = dg.(reverse(O.tensors))
    for i in 1:L
        W = Ws[i]
        M = ψ.tensors[i]

        @reduce N[(bᵢ₋₁, aᵢ₋₁), (bᵢ, aᵢ), σᵢ] :=  sum(σ′ᵢ) W[bᵢ₋₁, bᵢ, σᵢ, σ′ᵢ] * M[aᵢ₋₁, aᵢ, σ′ᵢ]
        push!(tensors, N)
    end
    adjoint(MPS{L, T}(tensors))
end
# Matrix Product Operators:2 ends here
