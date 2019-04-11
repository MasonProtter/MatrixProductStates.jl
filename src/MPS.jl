# Matrix Product States
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Matrix%20Product%20States][Matrix Product States:1]]
"""
    MPS{L, T<:Number}

Matrix product state on L sites. 

The `i`th tensor in the state has indices `[aⁱ⁻¹, aⁱ, σⁱ]` where
`(aⁱ⁻¹, aⁱ)` are bond indices and `σⁱ` is the physical index.

A four site MPS would be diagrammatically represented

    σ¹          σ²          σ³          σ⁴
    |           |           |           | 
    •--(a¹ a¹)--•--(a² a²)--•--(a³ a³)--•     

Note that `a⁰` and `aᴸ` must be of dimension 1.
"""
struct MPS{L, T<:Number} 
    tensors::Vector{Array{T,3}}
end

Base.isequal(ψ::MPS, ϕ::MPS)     = (isequal(ψ.tensors, ϕ.tensors))
Base.isapprox(ψ::MPS, ϕ::MPS)   = isapprox(ψ.tensors, ϕ.tensors)

Base.eltype(::Type{MPS{L, T}}) where {L, T} = T

Base.length(::MPS{L, T}) where {L, T} = L

Base.size(::MPS{L, T}) where {L, T} = (L,)
Base.getindex(ψ::MPS, i::Int) = getindex(ψ.tensors, i)

Base.:(*)(ψ::MPS{L, T}, x::Number) where {L, T} = MPS{L,T}(ψ.tensors .* x)
Base.:(*)(x::Number, ψ::MPS) = ψ * x
Base.:(/)(ψ::MPS{L,T}, x::Number) where {L, T} = MPS{L,T}(ψ.tensors ./ x)
Base.copy(ψ::MPS{L, T}) where {L, T} = MPS{L,T}(copy.(ψ.tensors))

function Base.randn(::Type{MPS{L, T}}, D::Int, d::Int) where {L, T}
    tensors = [randn(1, D, d), [randn(D, D, d) for _ in 2:(L-1)]..., randn(D, 1, d)]
    MPS{L, T}(tensors) |> leftcanonical |> rightcanonical
end

"""
    MPS(vs::Vector{Vector})
Create an `MPS` representing a product state (all bonds have dimension 1),
where each site is described by the corresponding element of `vs`.
"""
function MPS(vs::Vector{Vector{T}}) where {T}
    L = length(vs)

    tensrs = Vector{Array{T,3}}(undef, L)
    for i in 1:L
        tensrs[i] = reshape(copy(vs[i]), 1, 1, :)
    end

    MPS{L,T}(tensrs)
end

"""
    MPS(v::Vector, L)
Create an `MPS` for `L` sites representing a uniform product state (all bonds
have dimension 1), where each site is described by `v`.
"""
MPS(v::Vector, L) = MPS([v for _ in 1:L])

function Base.show(io::IO, ::MIME"text/plain", ψ::MPS{L, T}) where {L, T}
    d = length(ψ.tensors[2][1, 1, :])
    bonddims = [size(ψ[i][:, :, 1]) for i in 1:L]
    println(io, "Matrix product state on $L sites")
    _show_mps_dims(io, L, d, bonddims)
end

function Base.show(ψ::MPS{L, T}) where {L, T}
    d = length(ψ.tensors[2][1, 1, :])
    bonddims = [size(ψ[i][:, :, 1]) for i in 1:L]
    println("Matrix product state on $L sites")
    _show_mps_dims(L, d, bonddims)
end

function _show_mps_dims(io::IO, L, d, bonddims)
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

function Base.show(io::IO, ψ::MPS{L, T}) where {L, T}
    print(io, "MPS on $L sites")
end
# Matrix Product States:1 ends here



# #+HTML: <details><summary>Adjoint MPS</summary>
# #+HTML: <p>


# [[file:~/.julia/dev/MatrixProductStates/README.org::*Matrix%20Product%20States][Matrix Product States:2]]
function Base.adjoint(ψ::MPS{L, T}) where {L,T}
    Adjoint{T, MPS{L, T}}(ψ)
end

function Base.show(io::IO, ::MIME"text/plain", ψ::Adjoint{T, MPS{L, T}}) where {L, T}
    d = length(ψ.parent[2][1, 1, :])
    bonddims = reverse([reverse(size(ψ.parent[i][:, :, 1])) for i in 1:L])
    println(io, "Adjoint matrix product state on $L sites")
    _show_mps_dims(io, L, d, bonddims)
end

function Base.show(io::IO, ψ::Adjoint{T, MPS{L, T}}) where {L, T}
    print(io, "Adjoint MPO on $L sites")t
end

Base.size(::Adjoint{T, MPS{L, T}}) where {L, T} = (1, L)

function Base.getindex(ψ::Adjoint{T, MPS{L, T}}, args...) where {L, T}
    out = getindex(reverse(ψ.parent.tensors), args...)
    permutedims(conj.(out), (2, 1, 3))
end

adjoint_tensors(ψ::MPS) = reverse(conj.(permutedims.(ψ.tensors, [(2, 1, 3)])))
# Matrix Product States:2 ends here
