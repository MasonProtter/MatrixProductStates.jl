# Tests
# #+HTML:  <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Tests][Tests:1]]
using Test, MatrixProductStates

@testset "TFIM" begin
    g = 1.0; L = 7

    id = [1.0 0.0; 0.0 1.0]
    σˣ = [0.0 1.0; 1.0 0.0]
    σᶻ = [1.0 0.0; 0.0 -1.0]

    H_tnsr = zeros(3, 3, 2, 2)
    H_tnsr[1, 1, :, :] = id
    H_tnsr[2, 1, :, :] = -σᶻ
    H_tnsr[3, 1, :, :] = -g*σˣ
    H_tnsr[3, 2, :, :] = σᶻ
    H_tnsr[3, 3, :, :] = id
    H = MPO(H_tnsr, L)

    ψ = randn(MPS{L, Float64}, 100, 2)

    @test (ψ' * H) * ψ ≈ ψ' * (H * ψ)
    
    ψ̃ = compress(ψ, left, Dcut=80)[1] # Note: no actual information is lost in this 
                                        # compression because of the small size of the chain

    @test              ψ̃'ψ̃ ≈ 1
    @test          ψ'ψ/ψ'ψ ≈ ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ̃))/ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ))/ψ̃'ψ

    

end
# Tests:1 ends here
