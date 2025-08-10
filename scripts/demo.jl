#' # Metoda konjugiranih gradientov s predpogojevanjem
#' Avtor: Luka Bajić

#' ## Opis problema
#' Metoda konjugiranih gradientov je postopek za reševanje linearnega sistema enačb $Ax = b$,
#' ob predpostavki, da je matrika A pozitivno definitna.
#' 

#' ## Opis rešitve
#' ### Nepopolni razcep Choleskega

#' ### Metoda konjugiranih gradientov brez predpogojevanja

#' ### Metoda konjugiranih gradientov s predpogojevanjem


#' ## Primer uporabe


using MKG, Plots, LinearAlgebra, SparseArrays
tmp = [1 0 1; 0 2 0; 0 0 1]
Tv = sparse(tmp)
print(Tv)


I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]
A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
L = nep_chol(A)
x1, it1, res1 = conj_grad_baseline(A, b, vrniresid=true)
x2, it2, res2 = conj_grad(A, b, L, vrniresid=true)

plot(res1, label="brez predpogojevanja", title="Primerjava residualov")
plot!(res2, label="s predpogojevanjem", title="Primerjava residualov")

#' Kot vidimo

I = [1., 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
J = [1., 1, 2, 2, 3, 1, 3, 4, 1, 4, 5]
V = [5., -2, 5, -2, 5, -2, -2, 5, -2, -2, 5]
A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
L = nep_chol(A)
x, it = conj_grad(A, b, L, tol=0.1)

exp = [2/5, 19/25, -87/125, 51/625, 727/3125]
x




I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]

A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
x, it = conj_grad_baseline(A, b)
L = nep_chol(A)
x2, it = conj_grad(A, b, L)

isapprox(x, x2, atol=10e-10)
