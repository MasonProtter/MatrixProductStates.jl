# Utils
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Utils][Utils:1]]
export ⊗, realize

A ⊗ B = kron(A, B)

realize(x::Number) = error("Unrecognized numerical type")
realize(x::Real) = x
function realize(x::Complex; ϵ=1e-13)
    abs(imag(x)) < ϵ || error("Non-zero imaginary component")
    real(x)
end

dg(M::Array{T, 4}) where {T} = permutedims(conj.(M), (2, 1, 3, 4))
dg(M::Array{T, 3}) where {T} = permutedims(conj.(M), (2, 1, 3))

not(x) = ~x

function iseigenstate(ψ::MPS, H::MPO; ϵ=1e-5)
    ϕ = rightcanonical(ψ)
    realize(ϕ' * (H * H * ϕ) - (ϕ' * (H * ϕ))^2) < ϵ
end
# Utils:1 ends here
