# Iterative Ground State Search
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Iterative%20Ground%20State%20Search][Iterative Ground State Search:1]]
function R_exprs(ψ::MPS{L, T}, H::MPO{L, T}) where {L, T}
    R_exs = Array{T, 3}[]
    B = ψ[L]
    W = H[L]
    R_ex = ones(T, 1, 1, 1)
    for i in L:-1:2
        B = ψ[i]
        W = H[i]
        @tensor R_ex[bⁱ⁻¹, aⁱ⁻¹, aⁱ⁻¹′] := (conj.(B))[aⁱ⁻¹,aⁱ,σⁱ] * W[bⁱ⁻¹,bⁱ,σⁱ,σⁱ′] * B[aⁱ⁻¹′,aⁱ′,σⁱ′] * R_ex[bⁱ,aⁱ,aⁱ′]
        push!(R_exs, R_ex)
    end
    reverse(R_exs)
end

function sweep!(::Right, ψ::MPS{L, T}, H::MPO{L, T}, R_exs) where {L, T}
    L_exs = Array{T, 3}[]
    L_ex  = ones(T, 1, 1, 1)
    E = zero(T)
    for l in 1:(L-1)
        M    = ψ[l]
        Dˡ⁻¹, Dˡ, d = size(M)
        W    = H[l]
        R_ex = R_exs[l]

        @cast v[(σˡ, aˡ⁻¹, aˡ)] |= M[aˡ⁻¹, aˡ, σˡ]

        @reduce h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] |= sum(bˡ⁻¹, bˡ) begin
            L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R_ex[bˡ, aˡ, aˡ′]
        end strided
        h = collect(h)

        λ, Φ = eigs(h, v0=v, nev=1, which=:SR)
        E = λ[1]::T
        v⁰ = Φ[:,1]::Vector{T}

        @cast Mm[(σˡ, aˡ⁻¹), aˡ] := v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)

        U, S, V = svd(Mm)
        @cast A[aˡ⁻¹, aˡ, σˡ] |= U[(σˡ, aˡ⁻¹), aˡ] (σˡ:d, aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ)
        @cast A[aˡ⁻¹, aˡ, σˡ] |= U[(σˡ, aˡ⁻¹), aˡ] (σˡ:d, aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ)

        ψ.tensors[l] = A

        @tensor L_ex[bˡ, aˡ, aˡ′] := L_ex[bˡ⁻¹,aˡ⁻¹,aˡ⁻¹′] * (conj.(A))[aˡ⁻¹,aˡ,σˡ] * W[bˡ⁻¹,bˡ,σˡ,σˡ′] * A[aˡ⁻¹′,aˡ′,σˡ′]

        push!(L_exs, L_ex)

        SVp = Diagonal(S)*(V')
        Bp1 = ψ.tensors[l+1]
        @tensor Mp1[sⁱ⁻¹, aⁱ, σⁱ] := SVp[sⁱ⁻¹, aⁱ⁻¹] * Bp1[aⁱ⁻¹, aⁱ, σⁱ]
        ψ.tensors[l+1] = Mp1
    end
    ψ, L_exs, E
end

function sweep!(::Left, ψ::MPS{L, T}, H::MPO{L, T}, L_exs) where {L, T}
    R_exs = Array{T, 3}[]
    R_ex  = ones(T, 1, 1, 1)
    E = zero(T)
    for l in L:-1:2
        M = ψ[l]
        Dˡ⁻¹, Dˡ, d = size(M)
        W    = H[l]
        L_ex = L_exs[l-1]
        @cast v[(σˡ, aˡ⁻¹, aˡ)] |= M[aˡ⁻¹, aˡ, σˡ]

        @reduce h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] |= sum(bˡ⁻¹, bˡ) begin
            L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R_ex[bˡ, aˡ, aˡ′]
        end strided

        h = collect(h)

        λ, Φ = eigs(h, v0=v, nev=1, which=:SR)
        E = λ[1]::T
        v⁰ = Φ[:,1]::Vector{T}
        @cast Mm[aˡ⁻¹, (σˡ, aˡ)] |= v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
        U, S, V = svd(Mm)
        @cast B[aˡ⁻¹, aˡ, σˡ] |= V'[aˡ⁻¹, (σˡ, aˡ)] (σˡ:d)

        ψ.tensors[l] = B

        @tensor R_ex[bⁱ⁻¹, aⁱ⁻¹, aⁱ⁻¹′] := (conj.(B))[aⁱ⁻¹,aⁱ,σⁱ] * W[bⁱ⁻¹,bⁱ,σⁱ,σⁱ′] * B[aⁱ⁻¹′,aⁱ′,σⁱ′] * R_ex[bⁱ,aⁱ,aⁱ′]
        push!(R_exs, R_ex)

        US = U * Diagonal(S)
        Am1 = ψ.tensors[l-1]
        @tensor Mm1[aˡ⁻², sˡ⁻¹, σˡ⁻¹] :=  Am1[aˡ⁻², aˡ⁻¹′, σˡ⁻¹] * US[aˡ⁻¹′, sˡ⁻¹]
        ψ.tensors[l-1] = Mm1
    end
    ψ, R_exs, E
end

function ground_state(ψ::MPS{L, T}, H::MPO{L, T}; maxiter=40, quiet=false) where {L, T}
    ϕ = ψ |> copy
    
    quiet || println("Computing R expressions")

    R_exs = R_exprs(ψ, H)
    converged = false
    count     = 0
    E₀ = zero(T)
    while not(converged)
        quiet || println("Performing right sweep")
        ϕ, L_exs, _ = sweep!(right, ϕ, H, R_exs)

        quiet || println("Performing left sweep")
        ϕ, R_exs, E₀ = sweep!(left,  ϕ, H, L_exs)

        count += 1

        if iseigenstate(ϕ, H)
            quiet || println("Converged in $count iterations")
            converged = true
        elseif count >= maxiter
            @warn "Did not converge in $maxiter iterations"
            break
        end
    end
    ϕ, E₀
end

function iseigenstate(ψ::MPS, H::MPO; ϵ=1e-5)
    ϕ = rightcanonical(ψ)
    realize(ϕ' * (H * H * ϕ) - (ϕ' * (H * ϕ))^2) < ϵ
end
# Iterative Ground State Search:1 ends here
