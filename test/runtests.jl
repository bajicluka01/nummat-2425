using Airy
using Test
using SpecialFunctions
using LinearAlgebra

# točne vrednosti za prvih 5 ničel, pridobljene z Wolfram Alpha
n1 = -2.3381074104597670384891972524467354406385401456723878524838544372136680027
n2 = -4.0879494441309706166369887014573910602247646991085297549841608760251219468
n3 = -5.5205598280955510591298555129312935737972142806175251048328875769574946413
n4 = -6.7867080900717589987802463844961769660538824773934936165235290935562365568
n5 = -7.9441335871208531231382805557982685321406743969722148086438542857164486402
n = [n1, n2, n3, n4, n5]

k=3000

@testset "Pravilno izračunanih prvih 5 ničel z bisekcijo" begin
    n_test = airy_k_nicel(5)
    for i=1:5
        @test isapprox(n[i], n_test[i], atol=10e-10)
    end
end

@testset "Pravilno izračunanih prvih 5 ničel s tangentno metodo" begin
    n_test = airy_k_nicel(5, metoda="tangentna")
    for i=1:5
        @test isapprox(n[i], n_test[i], atol=10e-10)
    end
end

@testset "Pravilno izračunanih prvih 5 ničel z metodo regula falsi" begin
    n_test = airy_k_nicel(5, metoda="regula")
    for i=1:5
        @test isapprox(n[i], n_test[i], atol=10e-10)
    end
end

@testset "Pravilno izračunane ničle na intervalu [-8, 0]" begin
    n_test = airy_nicle_na_intervalu(-8.0)
    @test length(n_test) == 5
    for i=1:5
        @test isapprox(n[i], n_test[i], atol=10e-10)
    end
end

@testset "Pravilno izračunanih prvih $k ničel z bisekcijo" begin
    n_test = airy_k_nicel(k)
    for i=1:k
        @test isapprox(airyai(n_test[i]), 0, atol=10e-10)
    end
end

@testset "Pravilno izračunanih prvih $k ničel s tangentno metodo" begin
    n_test = airy_k_nicel(k, metoda="tangentna")
    for i=1:k
        @test isapprox(airyai(n_test[i]), 0, atol=10e-10)
    end
end

@testset "Pravilno izračunanih prvih $k ničel z metodo regula falsi" begin
    n_test = airy_k_nicel(k, metoda="regula")
    for i=1:k
        @test isapprox(airyai(n_test[i]), 0, atol=10e-10)
    end
end
