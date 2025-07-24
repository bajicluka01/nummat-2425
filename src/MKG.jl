module MKG
using LinearAlgebra
using SparseArrays
export conj_grad, conj_grad_basic, nep_chol

function nep_chol(A_::AbstractSparseMatrix)
    A = copy(A_)
    L = copy(A)

    for i=1:size(A,2)
        sum = 0
        for j=1:i-1
            sum = sum + L[i,j]^2
        end

        L[i,i] = sqrt(A[i,i]-sum)

        for j=i+1:size(A,2)
            # elemente nad diagonalo nastavimo na 0
            L[i,j] = 0

            for k=1:i-1
                A[j,i] = A[j,i] - L[j,k] * L[i,k]
            end

            if L[j,i] != 0
                L[j,i] = A[j,i] / L[i,i]
            end
        end
    end

    return dropzeros(L)
end

function obr_sub(U_::AbstractSparseMatrix, r_::AbstractVector)
    U = copy(U_)
    r = copy(r_)
    n = size(r,1)
    x = zeros(n)

    x[n] = r[n] / U[n, n]

    for i=n-1:-1:1
        s = r[i]
        for j=1+1:n
            s = s - U[i,j] * x[j]
        end
        x[i] = s / U[i,i]
    end

    return x
end


function conj_grad(A_::AbstractSparseMatrix, b_::AbstractVector, L_::AbstractSparseMatrix, vrniresid::Bool = false)
    A = copy(A_)
    b = copy(b_)
    L = copy(L_)

    print(obr_sub(L, b))

    # Začnemo s poljubno začetno vrednostjo x_0
    x_0 = zeros(size(b))

    # Nastavimo želeno toleranco
    eps = 1e-10

    # Preprečimo neskončno zanko, v primeru da metoda ne konvergira
    max_iter = 1000

    # Izračunamo začetni residual
    r_prev = b - A*x_0
    p = r_prev
    iter = 0
    x_out = x_0

    if vrniresid
        residuali = [norm(r_prev)]
    end

    while norm(r_prev)^2 > eps && iter < max_iter
        alpha_k = (r_prev' * r_prev) / (p' * A * p)
        x_out = x_out + alpha_k * p
        r_next = r_prev - alpha_k * A * p
        beta_k = (r_next' * r_next) / (r_prev' * r_prev)
        p = r_next + beta_k * p 
        r_prev = r_next
        iter = iter + 1

        if vrniresid
            push!(residuali, norm(r_prev))
        end
    end
        
    x_out = L' * x_out

    if vrniresid
        return x_out, iter, residuali
    end
    return x_out, iter
end

function conj_grad_basic(A_::AbstractSparseMatrix, b_::AbstractVector, vrniresid::Bool = false)
    A = copy(A_)
    b = copy(b_)

    # Začnemo s poljubno začetno vrednostjo x_0
    x_0 = zeros(size(b))

    # Nastavimo želeno toleranco
    eps = 1e-10

    # Preprečimo neskončno zanko, v primeru da metoda ne konvergira
    max_iter = 1000

    # Izračunamo začetni residual
    r_0 = b - A*x_0
    p = r_0
    iter = 0
    r_prev = r_0
    x_out = x_0

    if vrniresid
        residuali = [norm(r_prev)]
    end

    while norm(r_prev)^2 > eps && iter < max_iter
        alpha_k = (r_prev' * r_prev) / (p' * A * p)
        x_out = x_out + alpha_k * p
        r_next = r_prev - alpha_k * A * p
        beta_k = (r_next' * r_next) / (r_prev' * r_prev)
        p = r_next + beta_k * p 
        r_prev = r_next
        iter = iter + 1

        if vrniresid
            push!(residuali, norm(r_prev))
        end
    end

    if vrniresid
        return x_out, iter, residuali
    end
    return x_out, iter
end

end # module MKG
