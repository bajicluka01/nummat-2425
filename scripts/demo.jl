#' # Fresnelov integral

using Fresnel, Plots

x = 0.25
fresnel_cos(x)


xs = []
ys = []

for i=-1.5:0.1:1.5
    push!(xs, i)
    push!(ys, fresnel_cos(i))
end
plot(xs[:], ys[:], label="Fresnelov kosinus")




function fr(x)
    res = 0.0

    for n=0:21
        st = (-1)^n * (0.5*pi)^(2*n) * x^(4n+1)
        im = factorial(big(2*n)) * (4*n+1)
        tmp = st/im
        res = res + tmp
    end

    return res
end

fr(0.25)
