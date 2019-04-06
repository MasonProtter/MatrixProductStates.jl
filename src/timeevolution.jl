# Imaginary Time Evolution
# I don't think this works!
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Imaginary%20Time%20Evolution][Imaginary Time Evolution:1]]
# Fixme! this does not appear to find ground states!

function _MPO_handed_time_evolver(hs::Vector{Matrix{T}}, τ, L, d) where {T}
    tensors = Array{T, 4}[]
    for h in hs
        O = exp(-τ*h)
        @cast P[(σⁱ, σⁱ′), (σⁱ⁺¹, σⁱ⁺¹′)] |= O[(σⁱ, σⁱ⁺¹), (σⁱ′, σⁱ⁺¹′)] (σⁱ:d, σⁱ′:d)
        U, S, V = svd(P)

        @cast U[1, k, σⁱ, σⁱ′]     := U[(σⁱ, σⁱ′), k] * √(S[k])      (σⁱ:d)
        @cast Ū[k, 1, σⁱ⁺¹, σⁱ⁺¹′] := √(S[k]) * V'[k, (σⁱ⁺¹, σⁱ⁺¹′)] (σⁱ⁺¹:d)
        push!(tensors, U, Ū)
    end
    MPO{L, T}(tensors)
end

function MPO_time_evolvers(h1::Matrix, hi::Matrix, hL::Matrix, τ, L, d)
    if iseven(L)
        odd_hs  = [h1, [hi for _ in 3:2:(L-1)]...]
        even_hs = [[hi for i in 2:2:(L-1)]..., hL]
    else
        odd_hs  = [h1, [hi for _ in 3:2:(L-1)]..., hL]
        even_hs = [hi for i in 2:2:(L-1)]
    end
    
    Uodd  = _MPO_handed_time_evolver(odd_hs, τ, L, d)
    Ueven = _MPO_handed_time_evolver(even_hs, τ, L, d)
    Uodd, Ueven
end

function imag_time_evolution(ψ::MPS{L, T}, h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, 
                             β, N, Dcut) where {L, T}
    @warn "This probably still doesn't work!"
    τ = β/N
    d = length(ψ[1][1, 1, :])
    ϕ = ψ  # Ground state guess
    dir = left
    Uodd, Ueven = MPO_time_evolvers(h1, hi, hL, τ, L, d)
    for _ in 1:N
        ϕ1, dir = compress(Uodd  * ϕ,  dir, Dcut=Dcut)
        ϕ,  dir = compress(Ueven * ϕ1, dir, Dcut=Dcut)
        #ϕ,  dir = compress(Uodd  * ϕ2, dir, Dcut=Dcut)
    end
    ϕ
end
# Imaginary Time Evolution:1 ends here
