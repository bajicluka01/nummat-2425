module Fresnel
using LinearAlgebra

export fresnel_cos

function fresnel_cos(x::Float64, n_vozlov::Integer=100)
    function aux(x::Float64, n::Integer)
        l, w = gauss_laguerre(n)

        y = pi*x*x/2

        sum_f = 0.0
        sum_g = 0.0
        for i=1:n
            # definiramo vmesne spremenljivke, da se izognemo dvakratnemu računanju
            c = l[i]/y
            c2 = c^2+1
            faktor_f = 1/(sqrt(c)*(c2))
            faktor_g = sqrt(c)/(c2)

            sum_f = sum_f + w[i] * faktor_f
            sum_g = sum_g + w[i] * faktor_g
        end

        faktor = 2/(pi*x*x*pi*sqrt(2))
        sum_f = sum_f*faktor
        sum_g = sum_g*faktor

        return sum_f, sum_g
    end

    function gauss_laguerre(n::Integer)
        # skonstruiramo tridiagonalno matriko
        d = zeros(n)
        p = zeros(n-1)
        for i=1:n-1
            p[i] = i 
            d[i] = 2.0 * i - 1
        end
        d[n] = 2.0 * n - 1
        M = Tridiagonal(p, d, p)

        # izračunamo vozle in uteži
        l, v = eigen(M)
        w = v[1,:].^2
        return l, w
    end

    mult = 1
    if x == 0.0
        return x 
    elseif x <= 0 # ker velja C(-x) = -C(x)
        mult = -1 
        x = abs(x)
    end

    # za dovolj majhne x-e uporabimo kar potenčno vrsto, ki konvergira že za
    # dovolj majhne n-je (zanemarljivo s stališča časovne kompleksnosti)
    if abs(x) <= 1.5
        res = 0.0
        for n=0:12
            st = (-1)^n * (0.5*pi)^(2*n) * x^(4n+1)
            im = factorial(big(2*n)) * (4*n+1)
            res = res + st/im
        end
    else # sicer uporabimo pomožne funkcije
        res_prev = 0
        #while true
        f, g = aux(x, n_vozlov)
        param = pi*x^2/2
        res = 0.5 + f * sin(param) - g * cos(param)

        #    if abs(res-res_prev)
        #end
    end

    return mult * res
end

end # module Fresnel
