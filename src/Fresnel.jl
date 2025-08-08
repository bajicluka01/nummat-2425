module Fresnel
using LinearAlgebra

export fresnel_cos

function fresnel_cos(x::Float64)
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
    else
        res = 0
    end

    return mult * res
end

end # module Fresnel
