# Tests
# #+HTML:  <details><summary>Source</summary>
# #+HTML: <p>

# [[file:~/.julia/dev/MatrixProductStates/README.org::*Tests][Tests:1]]
using Test, MatrixProductStates, SparseArrays, Arpack

@testset "TFIM   " begin
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
    ψ = randn(MPS{L, Complex{Float64}}, 100, 2)
    
    ψ̃ = compress(ψ, left, Dcut=80)[1] # Note: no actual information is lost in this 
    # compression because of the small size of the chain

    @test              ψ̃'ψ̃ ≈ 1
    @test          ψ'ψ/ψ'ψ ≈ ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ̃))/ψ̃'ψ̃
    @test ((ψ'*(H*ψ))/ψ'ψ) ≈ (ψ̃' * (H * ψ))/ψ̃'ψ

    ϕ, E₀ = ground_state(ψ, H, quiet=true)
    @test ϕ' * H * H * ϕ ≈ (ϕ'*H*ϕ)^2
end

@testset "Hubbard" begin

    id = [1 0
          0 1]
    c  = [0 0
          1 0] #Anti commuting matrix
    c_up = c  ⊗ id
    c_dn = id ⊗ c
    id²  = id ⊗ id
    n_up = c_up' * c_up
    n_dn = c_dn' * c_dn

    P_up = (id² - 2c_up'*c_up) # Spin up parity operator
    P_dn = (id² - 2c_dn'*c_dn) # Spin down parity operator

    function H_hub(U, μ, L)
        W_tnsr = zeros(Complex{Float64}, 6, 6, 4, 4)
        W_tnsr[1, 1, :, :] = id²
        W_tnsr[2, 1, :, :] = c_up'
        W_tnsr[3, 1, :, :] = c_dn'
        W_tnsr[4, 1, :, :] = c_up
        W_tnsr[5, 1, :, :] = c_dn
        W_tnsr[6, 1, :, :] = U*(n_up * n_dn) - μ*(n_up + n_dn)
        W_tnsr[6, 2, :, :] =  c_up  * P_up  # Must multiply by the parity operator to get 
        W_tnsr[6, 3, :, :] =  c_dn  * P_dn  # correct off-site commutation relations!
        W_tnsr[6, 4, :, :] = -c_up' * P_up
        W_tnsr[6, 5, :, :] = -c_dn' * P_dn
        W_tnsr[6, 6, :, :] = id²
        MPO(W_tnsr, L)
    end

    function solve_hub(U, μ, L; retfull=true, quiet=true)
        H = H_hub(U, μ, L)
        ψ = randn(MPS{L, Complex{Float64}}, 100, 4)
        (ϕ, E₀), t, bytes = @timed ground_state(ψ, H, ϵ=1e-5, quiet=quiet)

        (ϕ=ϕ, E₀=E₀, H=H, t=t, Gbytes=bytes/1e9)
    end

    function Hub_ED(U, μ, L,)
        Û = U*(n_up * n_dn) - μ*(n_up + n_dn)
        c_dg_up(i) = foldl(⊗, sparse.([i==j ? c_up' : id² for j in 1:L]))
        cup(i)     = foldl(⊗, sparse.([i==j ? c_up  : id² for j in 1:L]))
        c_dg_dn(i) = foldl(⊗, sparse.([i==j ? c_dn' : id² for j in 1:L]))
        cdn(i)     = foldl(⊗, sparse.([i==j ? c_dn  : id² for j in 1:L]))
        Ûf(i)      = foldl(⊗, sparse.([i==j ? Û     : id² for j in 1:L]))
        function c_dg_c(i) 
            out = c_dg_up(i)*cup(i+1) + c_dg_dn(i)*cdn(i+1)
            out + out'
        end
        H = -sum(c_dg_c, 1:(L-1)) + sum(Ûf, 1:L)

        λ, ϕ = eigs(H, nev=1, which=:SR)
        (ϕ'H*ϕ)[]
    end

   
    U = 3.0; μ = -1.0; L = 4
    H = H_hub(U, μ, L)

    ϕ, E₀ = solve_hub(U, μ, L, retfull=true, quiet=true)
    @test ϕ' * H * H * ϕ ≈ (ϕ'*H*ϕ)^2  # Make sure energy is eigenvalue
    @test ϕ' * H * ϕ ≈ E₀              # make sure eigenvalue matches one produced by alogrithm
    @test ϕ' * H * ϕ ≈ Hub_ED(U, μ, L) # check against exact diagonalization
end
# Tests:1 ends here
