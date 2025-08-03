#' # Ničle Airyjeve funkcije
#' 
#' $\begin{bmatrix}A1 &A2 & A3\\A4 & A5 & A6 \end{bmatrix}$

using Airy, Plots
using SpecialFunctions

x = Float64(-3.5)
x2 = airyai(x)
x1, xs, ys, xs_, ys_ = airy_nicle_na_intervalu(x, 0.0005)
x1[1]
#start = 23381
#stop = 23384
start = 2338
stop = 2341
plot(xs_[start:stop], ys_[start:stop])
scatter!(xs, ys)

xs
x1


airy_k_nicel(10)



x = Float64(-8.5)
x2 = airyai(x)
x1 = airy_nicle_na_intervalu(x, 0.005)
x1[1]

airyai(x1[1])
abs(x1[1]-n1)

isapprox(x1[1], n1, atol=10e-4)
isapprox(x1[1], n1, atol=10e-6)
isapprox(x1[1], n1, atol=10e-10)

n1 = -2.3381074104597670384891972524467354406385401456723878524838544372136680027
     -2.3380074105267843

airyai(-6.7866080902074835)
airyai(-6.786608090112116)
airyai(-6.786608090016749)


n2 = -4.0879494441309706166369887014573910602247646991085297549841608760251219468
n3 = -5.5205598280955510591298555129312935737972142806175251048328875769574946413
n4 = -6.7867080900717589987802463844961769660538824773934936165235290935562365568
n5 = -7.9441335871208531231382805557982685321406743969722148086438542857164486402













n1 = -2.3381074104597670384891972524467354406385401456723878524838544372136680027
n2 = -4.0879494441309706166369887014573910602247646991085297549841608760251219468
n3 = -5.5205598280955510591298555129312935737972142806175251048328875769574946413
n4 = -6.7867080900717589987802463844961769660538824773934936165235290935562365568
n5 = -7.9441335871208531231382805557982685321406743969722148086438542857164486402
n = [n1, n2, n3, n4, n5]

n_test = airy_k_nicel(5, metoda="tangentna")
for i=1:5
    isapprox(n[i], n_test[i], atol=10e-10)
end


10e-10





#' test



a = 2.0
h = 0.001
i = 0.0

gamma23 = 1.3541179394264004169452880
gamma13 = 2.6789385347077476336556929

y_prej = [1/(3^(2/3) * gamma23), -1/(3^(1/3) * gamma13)]
x1 = airyai(i)
x2 = airy_premik(y_prej, i, h)[1]
isapprox(x1, x2)
while i <= a 
    x1 = airyai(i+h)
    x2 = airy_premik(y_prej, i, h)[1]
    #isapprox(x1, x2)
    println(abs(x1-x2))
    i = i + h
end 








a = -2.0
h = -0.001
i = 0.0

gamma23 = 1.3541179394264004169452880
gamma13 = 2.6789385347077476336556929

y_prej = [1/(3^(2/3) * gamma23), -1/(3^(1/3) * gamma13)]


x1 = airyai(i)
x2 = airy_premik(y_prej, i, h)[1]

while i >= a 
    x1 = airyai(i-h)
    x2 = airy_premik(y_prej, i, h)[1]
    #@test isapprox(x1, x2)
    #println(abs(x1-x2))
    if isapprox(x1, x2)
        println("OK")
    end
    i = i + h
end 













#backup
function airy_vrednost(x::Float64, h::Float64=0.01, vrnivmesne::Bool=false)
    vmesne = []
    gamma23 = 1.3541179394264004169452880
    gamma13 = 2.6789385347077476336556929

    y_prej = [1/(3^(2/3) * gamma23), -1/(3^(1/3) * gamma13)]

    if x < 0
        step = -h
    else
        step = h
    end

    i = 0.0
    while abs(i) <= abs(x)
        push!(vmesne, y_prej[1])
        y = exp(sigma(step, A1(i, step), A2(i, step))) * y_prej
        y_prej = y
        i = i + step
    end

    if vrnivmesne
        return y_prej[1], vmesne
    end

    return y_prej[1]
end
