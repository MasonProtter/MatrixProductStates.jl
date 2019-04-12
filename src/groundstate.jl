# Iterative Ground State Search
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Iterative%20Ground%20State%20Search][Iterative Ground State Search:1]]
function R_exprs(ψ::MPS{L, T}, H::MPO{L, T}) where {L, T}
    R_exs = Array{T, 3}[]
    R_ex = ones(T, 1, 1, 1)
    for l in L:-1:2
        R_ex = iterate_R_ex(ψ[l], H[l], R_ex) 
        push!(R_exs, R_ex)
    end
    reverse(R_exs)
end

# function preallocate_hs(ψ::MPS{L, T}) where {L, T}
#     h_tnsrs = map(ψ.tensors) do M
#         Dˡ⁻¹, Dˡ, d = size(M)
#         Array{T, 6}(undef, d, Dˡ⁻¹, Dˡ, d, Dˡ⁻¹, Dˡ)
#     end
# end


function sweep!(::Right, ψ::MPS{L, T}, H::MPO{L, T}, R_exs) where {L, T}
    L_exs = Array{T, 3}[]
    L_ex  = ones(T, 1, 1, 1)
    E = zero(T)
    for l in 1:(L-1)
        W = H[l]
        
        E, A, SVp = eigenproblem(right, ψ[l], L_ex, W, R_exs[l])
        ψ.tensors[l] = A

        L_ex = iterate_L_ex(A, W, L_ex)
        push!(L_exs, L_ex)

        Bp1 = ψ.tensors[l+1]
        @tensor Mp1[sⁱ⁻¹, aⁱ, σⁱ] := SVp[sⁱ⁻¹, aⁱ⁻¹] * Bp1[aⁱ⁻¹, aⁱ, σⁱ]
        ψ.tensors[l+1] = Mp1
    end
    L_exs, E
end

function sweep!(::Left, ψ::MPS{L, T}, H::MPO{L, T}, L_exs) where {L, T}
    R_exs = Array{T, 3}[]
    R_ex  = ones(T, 1, 1, 1)
    E = zero(T)
    for l in L:-1:2
        W = H[l]

        E, US, B = eigenproblem(left, ψ[l], L_exs[l-1], W, R_ex)
        ψ.tensors[l] = B

        R_ex = iterate_R_ex(B, W, R_ex) 
        push!(R_exs, R_ex)

        Am1 = ψ.tensors[l-1]
        @tensor Mm1[aˡ⁻², sˡ⁻¹, σˡ⁻¹] :=  Am1[aˡ⁻², aˡ⁻¹′, σˡ⁻¹] * US[aˡ⁻¹′, sˡ⁻¹]
        ψ.tensors[l-1] = Mm1
    end
    R_exs, E
end

function h_matrix(L_ex::Array{T,3}, W::Array{T,4}, R_ex::Array{T,3}) where {T}
    @tensor h[σˡ, aˡ⁻¹, aˡ, σˡ′, aˡ⁻¹′, aˡ′] := L_ex[bˡ⁻¹, aˡ⁻¹, aˡ⁻¹′] * W[bˡ⁻¹, bˡ, σˡ, σˡ′] * R_ex[bˡ, aˡ, aˡ′]
    @cast h[(σˡ, aˡ⁻¹, aˡ), (σˡ′, aˡ⁻¹′, aˡ′)] := h[σˡ, aˡ⁻¹, aˡ, σˡ′, aˡ⁻¹′, aˡ′]
end

function eigenproblem(dir::Direction, M::Array{T, 3}, L_ex::Array{T, 3}, W::Array{T, 4}, R_ex::Array{T, 3}) where {T}
    @cast v[(σˡ, aˡ⁻¹, aˡ)] |= M[aˡ⁻¹, aˡ, σˡ]
    
    h = h_matrix(L_ex, W, R_ex)
    λ, Φ = eigs(h, v0=v, nev=1, which=:SR)
    E  = λ[1]::T 
    v⁰ = (Φ[:,1])::Vector{T}

    (E, split_tensor(dir, v⁰, size(M))...)
end

function split_tensor(::Right, v⁰::Vector, (Dˡ⁻¹, Dˡ, d))
    @cast Mm[(σˡ, aˡ⁻¹), aˡ] := v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
    U, S, V = svd(Mm)
    @cast A[aˡ⁻¹, aˡ, σˡ] |= U[(σˡ, aˡ⁻¹), aˡ] (σˡ:d, aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ)
    A, Diagonal(S)*V'
end

function split_tensor(::Left, v⁰::Vector, (Dˡ⁻¹, Dˡ, d))
    @cast Mm[aˡ⁻¹, (σˡ, aˡ)] |= v⁰[(σˡ, aˡ⁻¹, aˡ)] (aˡ⁻¹:Dˡ⁻¹, aˡ:Dˡ, σˡ:d)
    U, S, V = svd(Mm)
    @cast B[aˡ⁻¹, aˡ, σˡ] |= V'[aˡ⁻¹, (σˡ, aˡ)] (σˡ:d)
    U*Diagonal(S), B
end

function iterate_R_ex(B, W, R_ex) where {T}
    @tensoropt R_ex′[bⁱ⁻¹, aⁱ⁻¹, aⁱ⁻¹′] := (conj.(B))[aⁱ⁻¹,aⁱ,σⁱ] * W[bⁱ⁻¹,bⁱ,σⁱ,σⁱ′] * B[aⁱ⁻¹′,aⁱ′,σⁱ′] * R_ex[bⁱ,aⁱ,aⁱ′]
end

function iterate_L_ex(A, W, L_ex) where {T}
    @tensoropt L_ex′[bˡ, aˡ, aˡ′] := L_ex[bˡ⁻¹,aˡ⁻¹,aˡ⁻¹′] * (conj.(A))[aˡ⁻¹,aˡ,σˡ] * W[bˡ⁻¹,bˡ,σˡ,σˡ′] * A[aˡ⁻¹′,aˡ′,σˡ′]
end


"""
    ground_state(ψ::MPS{L, T}, H::MPO{L, T}; maxiter=10, quiet=false, ϵ=1e-8) where {L, T}

Perform the finite system density matrix renormalization group
algorithm. First this will build up the R expressions, then do right
and left sweeps until either
 1) The state converges to an eigenstate `ϕ` such that
    ϕ' * H * H * ϕ ≈ (ϕ' * H * ϕ) 
to the requested tolerance `ϵ`
 2) The energy eigenvalue stops changing (possible signaling the algorithm is 
stuck in a local minimum)
 3) The number of full (right and left) sweeps exceeds `maxiter`. 

Setting `quiet=true` will suppress notifications about the algorithm's
progress but *not* warnings due to non-convergence.
"""
function ground_state(ψ::MPS{L, T}, H::MPO{L, T}; maxiter=10, quiet=false, ϵ=1e-8) where {L, T}
    ϕ = ψ |> copy

    quiet || println("Computing R expressions")
    R_exs = R_exprs(ψ, H)

    converged = false
    count     = 0
    E₀ = zero(T)
    enable_cache(maxsize=5*10^9)
    while not(converged)
        quiet || println("Performing right sweep")
        L_exs, E₀′ = sweep!(right, ϕ, H, R_exs)

        quiet || println("Performing left sweep")
        R_exs, E₀  = sweep!(left,  ϕ, H, L_exs)

        count += 1
        if iseigenstate(ϕ, H, ϵ=ϵ)
            quiet || println("Converged in $count iterations")
            converged = true
        elseif count > 1 && E₀ ≈ E₀′
                @warn """
Energy eigenvalue converged but state is not an eigenstate.
Consider either lowering your requested tolerance or 
implementing a warm-up algorithm to avoid local minima.
"""
            break
        elseif count >= maxiter
            @warn "Did not converge in $maxiter iterations"
            break
        end
    end
    clear_cache()
    ϕ, E₀
end


function iseigenstate(ψ::MPS, H::MPO; ϵ=1e-8)
    ϕ = rightcanonical(ψ)
    isapprox(ϕ' * (H * H * ϕ), (ϕ' * (H * ϕ))^2, rtol=ϵ)
end
# Iterative Ground State Search:1 ends here
