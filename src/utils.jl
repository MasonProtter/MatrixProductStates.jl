# Utils
# #+HTML: <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Utils][Utils:1]]
export ⊗, realize

abstract type Direction end

struct Left  <: Direction end # Often useful to dispatch on direction an algorithm is going
struct Right <: Direction end

const left  = Left()
const right = Right()

A ⊗ B = kron(A, B)

realize(x::Number) = error("Unrecognized numerical type")
realize(x::Real) = x
function realize(x::Complex; ϵ=1e-10)
    abs(imag(x)) < ϵ || error("Non-zero imaginary component, $(imag(x))")
    real(x)
end

dg(M::Array{T, 4}) where {T} = permutedims(conj.(M), (2, 1, 3, 4))
dg(M::Array{T, 3}) where {T} = permutedims(conj.(M), (2, 1, 3))

not(x) = ~x
# Utils:1 ends here
