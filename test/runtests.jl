# Tests
# #+HTML:  <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Tests][Tests:1]]
using Test, MatrixProductStates

@testset "TFIM" begin
    g = 1.0; L = 7

    function H_TFIM(g, L)
        id = [1  0; 
              0  1]
        σˣ = [0  1; 
              1  0]
        σᶻ = [1  0; 
              0 -1]
        W_tnsr = zeros(Complex{Float64}, 3, 3, 2, 2)
        W_tnsr[1, 1, :, :] = id    
        W_tnsr[2, 1, :, :] = -σᶻ  
        W_tnsr[3, 1, :, :] = -g*σˣ
        W_tnsr[3, 2, :, :] = σᶻ   
        W_tnsr[3, 3, :, :] = id   

        return MPO(W_tnsr, L)
    end
    H = H_TFIM(g, L)
    ψ = randn(MPS{L, Float64}, 100, 2)

    @test (ψ' * H) * ψ ≈ ψ' * (H * ψ)
    
    ψ̃ = compress(ψ, left, Dcut=80)[1] # Note: no actual information is lost in this 
                                        # compression because of the small size of the chain

    @test              ψ̃'ψ̃ ≈ 1
    @test          ψ'ψ/ψ'ψ ≈ ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ̃))/ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ))/ψ̃'ψ

    ϕ, E₀ = ground_state(ψ, H, quiet=true)
    @test ϕ' * H * H * ϕ ≈ (ϕ'*H*ϕ)^2
end
# Tests:1 ends here
