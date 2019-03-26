#---------------------------------------------------------------------
# Imaginary Time Evolution Assuming only nearest neighbor interactions

"""
     MPO_odd_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
"""
function MPO_odd_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
    tensors = Array{T, 4}[]


    if iseven(L)
    end
end
    
"""
    MPO_even_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
"""
function MPO_even_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
    tensors = Array{T, 4}[]

    if iseven(L)
    end
end

"""
    MPO_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
"""
function MPO_time_evolver(h1::Matrix{T}, hi::Matrix{T}, hL::Matrix{T}, τ, L) where {T}
    Uoddhalf = MPO_odd_time_evolver(h1, hi, hL, τ/2, L)
    Ueven    = MPO_even_time_evolver(h1, hi, hL, τ, L)
    compress(Uoddhalf * Ueven * Uoddhalf)
end
