# Iterative Ground State Search
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Iterative%20Ground%20State%20Search][Iterative Ground State Search:1]]
function R_exprs(ψ::MPS{L, T}, H::MPO{L, T}) where {L, T}
    Rs = Array{T, 3}[]
    B = ψ[L]
    W = H[L]
    @reduce R[bᴸ⁻¹, aᴸ⁻¹, aᴸ⁻¹′] := sum(σᴸ, σᴸ′, bᴸ, aᴸ,  aᴸ′) begin 
        (conj.(B))[aᴸ⁻¹, aᴸ, σᴸ] * W[bᴸ⁻¹, bᴸ, σᴸ, σᴸ′] * B[aᴸ⁻¹′, aᴸ′, σᴸ′]
    end
    push!(Rs, R)
    for i in (L-1):-1:2
        B = ψ[i]
        W = H[i]
        @reduce R[bⁱ⁻¹, aⁱ⁻¹, aⁱ⁻¹′] := sum(σⁱ, σⁱ′, bⁱ, aⁱ,  aⁱ′) begin 
            (conj.(B))[aⁱ⁻¹,aⁱ,σⁱ] * W[bⁱ⁻¹,bⁱ,σⁱ,σⁱ′] * B[aⁱ⁻¹′,aⁱ′,σⁱ′] * R[bⁱ,aⁱ,aⁱ′]
        end
        push!(Rs, R)
    end
    reverse(Rs)
end

function sweep!(::Right, ψ::MPS{L, T}, H::MPO{L, T}, R_exs, Dcut) where {L, T}
    L_exs = Array{T, 3}[]
    L_ex  = ones(T, 1, 1, 1)
    for l in 1:(L-1)
        M    = ψ[l]
        Dˡ⁻¹, Dˡ, d = size(M)
        W    = H[l]
        R_ex = R_exs[l]
        @cast  v[(σˡ, aˡ⁻¹, aˡ)] |= M[aˡ⁻¹, aˡ, σˡ]

        @reduce h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] |= sum(bˡ⁻¹, bˡ) begin
            L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R_ex[bˡ, aˡ, aˡ′]
        end

        λ, Φ = eigs(h, v0=v, nev=1, which=:SR)
        E = λ[1]
        v⁰ = Φ[:,1]

        @cast Mm[(σˡ, aˡ⁻¹), aˡ] |= v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
        U, S, V = psvd(Mm, rank=Dcut)
        @cast A[aˡ⁻¹, aˡ, σˡ] |= U[(σˡ, aˡ⁻¹), aˡ] (σˡ:d)

        @reduce L_ex[bˡ, aˡ, aˡ′] := sum(σˡ, σˡ′, bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′) begin
            L_ex[bˡ⁻¹,aˡ⁻¹,aˡ⁻¹′] * (conj.(A))[aˡ⁻¹,aˡ,σˡ] * W[bˡ⁻¹,bˡ,σˡ,σˡ′] * A[aˡ⁻¹′,aˡ′,σˡ′]
        end
        push!(L_exs, L_ex)

        SVp = Diagonal(S)*(V')
        Bp1 = ψ.tensors[l+1]
        @tensor Mp1[sⁱ⁻¹, aⁱ, σⁱ] := SVp[sⁱ⁻¹, aⁱ⁻¹] * Bp1[aⁱ⁻¹, aⁱ, σⁱ]
        ψ.tensors[l+1] = Mp1
    end
    return ψ, L_exs
end

# function sweep!(::Right, ψ::MPS{L, T}, H::MPO{L, T}, Rs, Dcut) where {L, T}
#     L_exs = Array{T, 3}[]
#     B = ψ[1]
#     D⁰, D¹, d = size(B)
#     W = H[1]
#     R = Rs[1]
#     @cast v[(σ¹, a⁰, a¹)] |= B[a⁰, a¹, σ¹]
#     @reduce h[(σ¹, a¹), (σ¹′, a¹′)] |= sum(b⁰, b¹) begin
#         W[b⁰, b¹, σ¹, σ¹′] * R[b¹, a¹, a¹′]
#     end
#     λ⁰, v⁰ = eigs(h, v0=v, nev=1, which=:SR)
#     @cast M[(σ¹, a⁰), a¹] |= v⁰[(σ¹, a⁰, a¹)] (a⁰:D⁰, a¹:D¹, σ¹:d)
#     U, S, V = psvd(M, rank=Dcut)
#     @cast   A[a⁰, a¹, σ¹]  |= U[(σ¹, a⁰), a¹] (σ¹:d)
#     @reduce L_ex[b¹, a¹, a¹′] := sum(σ¹, σ¹′, b⁰, a⁰, a⁰′) begin
#         (conj.(A))[a⁰, a¹, σ¹] * W[b⁰, b¹, σ¹, σ¹′] * A[a⁰′, a¹′, σ¹′]
#     end
#     push!(L_exs, L_ex)
#     ψ.tensors[1] = A # Mutate ψ

#     for l in 2:(L-1)
#         B = ψ[l]
#         W = H[l]
#         R = Rs[l]
#         SVp = Diagonal(S)*V'
#         @reduce M[aⁱ⁻¹, aⁱ, σⁱ]   := sum(aⁱ⁻¹′) SVp[aⁱ⁻¹, aⁱ⁻¹′] * B[aⁱ⁻¹′, aⁱ, σⁱ]
#         @cast   v[(σⁱ, aⁱ⁻¹, aⁱ)] |= M[aⁱ⁻¹, aⁱ, σⁱ]
#         @reduce h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] |= sum(bˡ⁻¹, bˡ) begin
#             L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R[bˡ, aˡ, aˡ′]
#         end
#         sizem1 = size(M)
#         #println((l=l, sizeB=size(B), sizeM=size(M), sizeR=size(R), sizeRp1=size(Rs[l-1])))
        
#         λ⁰, v⁰ = eigs(h, v0=v, nev=1, which=:SR)
#         Dˡ⁻¹, Dˡ, d = size(M)       
#         @cast M[aˡ⁻¹, aˡ, σˡ]    |= v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
#         @cast Mm[(σˡ, aˡ⁻¹), aˡ] |= M[aˡ⁻¹, aˡ, σˡ]

#         U, S, V = psvd(Mm, rank=Dcut)
#         #@cast A[aⁱ⁻¹, aⁱ, σⁱ] |= U[(σⁱ, aⁱ⁻¹), aⁱ] (σⁱ:d)
#         @cast A[aˡ⁻¹, aˡ, σˡ] |= U[(σˡ, aˡ⁻¹), aˡ] (σˡ:d)
#         println((sizeB =size(B), sizem1=sizem1, sizem2=size(M), sizeA=size(A), sizeU=size(U)))
#         @reduce L_ex[bˡ, aˡ, aˡ′] := sum(σˡ, σˡ′, bˡ⁻¹, aˡ⁻¹,  aˡ⁻¹′) begin
#             L_ex[bˡ⁻¹,aˡ⁻¹,aˡ⁻¹′] * (conj.(A))[aˡ⁻¹,aˡ,σˡ] * W[bˡ⁻¹,bˡ,σˡ,σˡ′] * A[aˡ⁻¹′,aˡ′,σˡ′]
#         end
#         push!(L_exs, L_ex)
#         ψ.tensors[l] = A
#     end
#     B = ψ.tensors[L]
#     @reduce M[aᴸ⁻¹, aᴸ, σᴸ] := sum(aᴸ⁻¹′) (Diagonal(S)*V')[aᴸ⁻¹, aᴸ⁻¹′] * B[aᴸ⁻¹′, aᴸ, σᴸ]
#     ψ.tensors[L] = M
#     ψ, L_exs, λ⁰[1]
# end

function sweep!(::Left, ψ::MPS{L, T}, H::MPO{L, T}, L_exs, Dcut) where {L, T}
    R_exs = Array{T, 3}[]
    A = ψ[L]
    Dᴸ⁻¹, Dᴸ, d = size(A)
    W = H[L]
    L_ex = L_exs[L-1]
    @cast v[(σᴸ, aᴸ⁻¹, aᴸ)] |= A[aᴸ⁻¹, aᴸ, σᴸ]
    @reduce h[(σᴸ, aᴸ⁻¹), (σᴸ′, aᴸ⁻¹′)] |= sum(bᴸ⁻¹, bᴸ) begin
        L_ex[bᴸ⁻¹, aᴸ⁻¹, aᴸ⁻¹′] * W[bᴸ⁻¹, bᴸ, σᴸ, σᴸ′]
    end
    λ⁰, v⁰ = eigs(h, v0=v, nev=1, which=:SR)
    @cast M[aᴸ⁻¹, (σᴸ, aᴸ)] |= v⁰[(σᴸ, aᴸ⁻¹, aᴸ)] (aᴸ⁻¹:Dᴸ⁻¹, aᴸ:Dᴸ, σᴸ:d)

    U, S, V = psvd(M, rank=Dcut)

    @cast B[aᴸ⁻¹, aᴸ, σᴸ] |= V'[aᴸ⁻¹, (σᴸ, aᴸ)] (σᴸ:d)

    @reduce R_ex[bᴸ⁻¹, aᴸ⁻¹, aᴸ⁻¹′] := sum(σᴸ, σᴸ′, bᴸ, aᴸ,  aᴸ′) begin 
        (conj.(B))[aᴸ⁻¹, aᴸ, σᴸ] * W[bᴸ⁻¹, bᴸ, σᴸ, σᴸ′] * B[aᴸ⁻¹′, aᴸ′, σᴸ′]
    end
    push!(R_exs, R_ex)
    for l in (L-1):-1:2
        A = ψ[l]
        Dᴸ⁻¹, Dᴸ, d = size(A)
        W = H[l]
        L_ex = L_exs[l-1]
        US = U * Diagonal(S)
        @reduce M[aˡ⁻¹, aˡ, σˡ]   := sum(aˡ′) A[aˡ⁻¹, aˡ′, σˡ] * US[aˡ′, aˡ]
        @cast   v[(σˡ, aˡ⁻¹, aˡ)] |= M[aˡ⁻¹, aˡ, σˡ]

        @reduce h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] |= sum(bˡ⁻¹, bˡ) begin
            L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R_ex[bˡ, aˡ, aˡ′]
        end
        
        λ⁰, v⁰ = eigs(h, v0=v, nev=1, which=:SR)
        Dˡ⁻¹, Dˡ, d = size(M)       
        @cast M[aˡ⁻¹, (σˡ, aˡ)] |= v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
        U, S, V = psvd(M, rank=Dcut)

        @cast B[aˡ⁻¹, aˡ, σˡ] |= V'[aˡ⁻¹, (σˡ, aˡ)] (σˡ:d)

        @reduce R_ex[bⁱ⁻¹, aⁱ⁻¹, aⁱ⁻¹′] := sum(σⁱ, σⁱ′, bⁱ, aⁱ,  aⁱ′) begin 
            (conj.(B))[aⁱ⁻¹,aⁱ,σⁱ] * W[bⁱ⁻¹,bⁱ,σⁱ,σⁱ′] * B[aⁱ⁻¹′,aⁱ′,σⁱ′] * R_ex[bⁱ,aⁱ,aⁱ′]
        end
        push!(R_exs, R_ex)
    end
    A = ψ.tensors[1]
    @reduce M[a⁰, a¹, σ¹] := sum(a¹′) A[a⁰, a¹′, σ¹] * (U * Diagonal(S))[a¹′, a¹]
    ψ.tensors[1] = M
    ψ, reverse(R_exs), λ⁰[1]
end

not(x) = ~x

function isconverged(ψ::MPS, H::MPO; ϵ=1e-2)
    ϕ = rightcanonical(ψ)
    realize(ϕ' * (H * H * ϕ) - (ϕ' * (H * ϕ))^2) < ϵ
end

function ground_state(ψ::MPS{L, T}, H::MPO{L, T}, Dcut) where {L, T}
    ψ = ψ |> copy
    R_exs = R_exprs(ψ, H)
    converged = false
    count     = 0
    while not(converged)
        #println(size.(ψ.tensors))
        ψ, L_exs = sweep!(right, ψ, H, R_exs, Dcut)
        #println(size.(ψ.tensors))
        ψ, R_exs, E₀′ = sweep!(left,  ψ, H, L_exs, Dcut)
        #println(size.(ψ.tensors))
        #println(size.(R_exs))

        count += 1
        @show count
        if isconverged(ψ, H)
            converged = true
        elseif count >= 100
            @warn "Did not converge in 200 iterations"
            break
        end
    end
    ψ, E₀′
end

realize(x::Number) = error("Unrecognized numerical type")
realize(x::Real) = x
function realize(x::Complex)
    abs(imag(x)) < 1e-13 || error("Non-zero imaginary component")
    real(x)
end
# Iterative Ground State Search:1 ends here
