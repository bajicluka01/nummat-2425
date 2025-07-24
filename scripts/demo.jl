#' # Metoda konjugiranih gradientov s predpogojevanjem

#' Metoda konjugiranih gradientov je postopek za reševanje linearnega sistema enačb $Ax = b$,
#' ob predpostavki, da je matrika A pozitivno definitna.

using MKG, Plots, LinearAlgebra, SparseArrays
tmp = [1 0 1; 0 2 0; 0 0 1]
Tv = sparse(tmp)
print(Tv)


I = [1, 1, 2, 2, 2, 3, 3]
J = [1, 2, 1, 2, 3, 2, 3]
V = [2, -1, -1, 2, -1, -1, 2]
A = sparse(I, J, V)
b = [1, 1, 1]
x, it = conj_grad(A, b)

x1, it1, res1 = conj_grad(A, b, false, true)
x2, it2, res2 = conj_grad(A, b, true, true)

plot(res1, label="brez predpogojevanja", title="Primerjava residualov")
plot!(res2, label="s predpogojevanjem", title="Primerjava residualov")


I = [1, 1, 2, 2, 2, 3, 3, 4]
J = [1, 2, 1, 2, 3, 2, 3, 4]
V = [2, -1, -1, 2, -1, -1, 2, 4]
A = sparse(I, J, V)
b = [3.5, -2, 1, 3]
x, it = conj_grad(A, b)


I = [1., 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
J = [1., 1, 2, 2, 3, 1, 3, 4, 1, 4, 5]
V = [5., -2, 5, -2, 5, -2, -2, 5, -2, -2, 5]
A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
x, it = conj_grad(A, b)
x, it = conj_grad(A, b)

x, it = conj_grad_basic(A, b)







I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]

A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
x, it = conj_grad_basic(A, b)

L = nep_chol(A)
x, it = conj_grad(A, b, L)