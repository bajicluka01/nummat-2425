using MKG
using Test
using SparseArrays

@testset "Ujemanje dimenzij" begin
    A = spzeros(3,3)
    b = [1, 1, 1]
    x = conj_grad_baseline(A, b)
    @test b == [1, 1, 1]
end

@testset "3x3 matrika" begin
    I = [1., 1, 2, 2, 2, 3, 3]
    J = [1., 2, 1, 2, 3, 2, 3]
    V = [2., -1, -1, 2, -1, -1, 2]
    A = sparse(I, J, V)
    b = [1, 1, 1]
    x, it = conj_grad_baseline(A, b)

    @test isapprox(x, [1.5, 2, 1.5], atol=10e-10)
    @test it < 1000
end

@testset "4x4 matrika" begin
    I = [1., 1, 2, 2, 2, 3, 3, 4]
    J = [1., 2, 1, 2, 3, 2, 3, 4]
    V = [2., -1, -1, 2, -1, -1, 2, 4]
    A = sparse(I, J, V)
    b = [3.5, -2, 1, 3]
    x, it = conj_grad_baseline(A, b)

    @test isapprox(x, [15/8, 1/4, 5/8, 3/4], atol=10e-10)
    @test it < 1000
end

# Primer iz https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
@testset "Nepopolni razcep Choleskega" begin
    I = [1., 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
    J = [1., 1, 2, 2, 3, 1, 3, 4, 1, 4, 5]
    V = [5., -2, 5, -2, 5, -2, -2, 5, -2, -2, 5]
    A = sparse(I, J, V)
    L = nep_chol(A)

    I_res = [1., 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5]
    J_res = [1., 2, 4, 5, 1, 2, 3, 4, 5, 2, 3, 4, 1, 2, 3, 4, 5, 1, 2, 4, 5]
    V_res = [5, -2, -2, -2, -2, 5, -2, 0.8, 0.8, -2, 5, -2, -2, 0.8, -2, 5, -2, -2, 0.8, -2, 5]
    res = sparse(I_res, J_res, V_res, 5, 5)

    @test isapprox(L*L', res, atol=10e-10)
end

@testset "MKG brez predpogojevanja" begin
    I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
    J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
    V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]

    A = sparse(I, J, V)
    b = [2, 3, -5, 1, 0.2]
    x, it = conj_grad_baseline(A, b)

    @test isapprox(x, [70/1279, 1880/1279, -6, 2, 13/3], atol=10e-10)
end

@testset "MKG s predpogojevanjem" begin
    I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
    J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
    V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]
    A = sparse(I, J, V)
    b = [2, 3, -5, 1, 0.2]
    L = nep_chol(A)
    x, it = conj_grad(A, b, L)

    @test isapprox(x, [70/1279, 1880/1279, -6, 2, 13/3], atol=10e-10)
end
