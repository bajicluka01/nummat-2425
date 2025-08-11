module MKG
using LinearAlgebra
using SparseArrays
export conj_grad, conj_grad_baseline, nep_chol

"""
    L = nep_chol(A)

Izračunaj nepopolni razcep Choleskega matrike `A`.
"""
function nep_chol(A::AbstractSparseMatrix)
    I_, J_, V_ = findnz(A)
    # odsranimo elemente nad diagonalo
    for k=1:length(J_)
        if J_[k] > I_[k]
            V_[k] = 0
        end
    end
    L = dropzeros(sparse(I_, J_, V_))

    for i=1:size(A,2)
        L = dropzeros(L)
        I, J, V = findnz(L)

        sum = 0
        for k=1:length(J)
            if I[k] == i && J[k] < i
                sum = sum + V[k]*V[k]
            end
        end

        L[i,i] = sqrt(A[i,i]-sum)

        for j=i+1:size(A,2)
            sum2 = A[j,i]

            for k=1:i-1
                sum2 = sum2 - L[j,k]*L[i,k]
            end

            if L[j,i] != 0
                L[j,i] = sum2/L[i,i]
            end
        end
    end

    return dropzeros(L)
end

"""
    x_out, iter = conj_grad(A, b, L)

Izvedi metodo konjugiranih gradientov nad matriko `A` in vektorjem `b`, s predpogojevalcem `L`.
"""
function conj_grad(A::AbstractSparseMatrix, b::AbstractVector, L::AbstractSparseMatrix; tol::Float64=10e-10, vrniresid::Bool = false)
    function obr_sub(U, r::AbstractVector)
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

    function prema_sub(L::AbstractSparseMatrix, b::AbstractVector)
        n = size(b,1)
        y = zeros(n)

        y[1] = b[1] / L[1, 1]
        for i=2:n 
            s = b[i]
            for j=1:i-1
                s = s - L[i,j] * y[j]
            end
            y[i] = s / L[i,i]
        end
        return y
    end

    # vrne z, da velja Mz=r, brez da bi računali M (velja M=LL^T)
    function sistem(L::AbstractSparseMatrix, r::AbstractVector)
        # La = r 
        a = prema_sub(L, r)

        # L^Tz = a 
        z = obr_sub(L', a)
        return z
    end

    # Začnemo s poljubno začetno vrednostjo x_0
    x_out = zeros(size(b))

    # Preprečimo neskončno zanko, v primeru da metoda ne konvergira
    max_iter = 1000

    # Izračunamo začetni residual
    r_prev = b - A*x_out
    z_prev = sistem(L, r_prev)
    p = z_prev
    iter = 1

    if vrniresid
        residuali = [norm(r_prev)]
    end

    while iter < max_iter
        alpha_k = (r_prev' * z_prev) / (p' * A * p)
        x_out = x_out + alpha_k * p
        r_next = r_prev - alpha_k * A * p
        
        if vrniresid
            push!(residuali, norm(r_next))
        end

        if norm(r_next, Inf) < tol
            break
        end

        z_next = sistem(L, r_next)
        beta_k = (r_next' * z_next) / (r_prev' * z_prev)
        p = z_next + beta_k * p

        r_prev = r_next
        z_prev = z_next
        iter = iter + 1
    end

    if vrniresid
        return x_out, iter, residuali
    end
    return x_out, iter
end

"""
    x_out, iter = conj_grad(A, b)

Izvedi metodo konjugiranih gradientov nad matriko `A` in vektorjem `b`, brez predpogojevanja.
"""
function conj_grad_baseline(A::AbstractSparseMatrix, b::AbstractVector; tol::Float64=10e-10, vrniresid::Bool = false)
    # Začnemo s poljubno začetno vrednostjo x_0
    x_0 = zeros(size(b))

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

    while iter < max_iter
        alpha_k = (r_prev' * r_prev) / (p' * A * p)
        x_out = x_out + alpha_k * p
        r_next = r_prev - alpha_k * A * p
        beta_k = (r_next' * r_next) / (r_prev' * r_prev)
        p = r_next + beta_k * p 

        iter = iter + 1

        if vrniresid
            push!(residuali, norm(r_next))
        end

        if norm(r_next, Inf) < tol
            break
        end

        r_prev = r_next
    end

    if vrniresid
        return x_out, iter, residuali
    end
    return x_out, iter
end

end # module MKG
