#' # Fresnelov integral

using Fresnel, Plots

x = 1.75
fresnel_cos(x)


xs = []
ys = []
for i=-10:0.05:10
    push!(xs, i)
    push!(ys, fresnel_cos(i))
end
plot(xs[:], ys[:], label="Fresnelov kosinus")

n = 1000
x = 2.0
c = fresnel_cos(x, n)
isapprox(c, 0.4728183813385479287835767887961982009647206121516617412533085818867, atol=5*10e-11)


c = fresnel_cos(x, 1000)
c = fresnel_cos(x, 2000)
c = fresnel_cos(x, 3000)
c = fresnel_cos(x, 4000)
c = fresnel_cos(x, 5000)
c = fresnel_cos(x, 6000)
c = fresnel_cos(x, 7000)
c = fresnel_cos(x, 8000)
c = fresnel_cos(x, 9000)
c = fresnel_cos(x, 10000)
