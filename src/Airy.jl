module Airy
using LinearAlgebra
using Plots

export airy_k_nicel, airy_nicle_na_intervalu

function airy_premik(y_prej::Vector, x_prej::Float64, h::Float64=0.01)
    function A1(xk::Float64, h::Float64)
        # c = 1/2 - sqrt(3)/6
        c = 0.2113248654051871177454256097490212721761991
        tmp = xk + (c * h)
        return [0 1; tmp 0]
    end

    function A2(xk::Float64, h::Float64)
        # c = 1/2 + sqrt(3)/6
        c = 0.7886751345948128822545743902509787278238008
        tmp = xk + (c * h)
        return [0 1; tmp 0]
    end

    function sigma(h::Float64, A1::Matrix, A2::Matrix)
        # c = sqrt(3)/12
        c = 0.1443375672974064411272871951254893639119004
        return h/2 * (A1 + A2) - c * h*h * kom(A1, A2)
    end

    function kom(A::Matrix, B::Matrix)
        return A*B - B*A
    end

    y = exp(sigma(h, A1(x_prej, h), A2(x_prej, h))) * y_prej
    return y
end

"""
najde prvih k ničel
"""
function airy_k_nicel(k::Integer; h::Float64=0.01, metoda::String="bisekcija", max_korakov::Integer=10^12)
    nicle = []

    y_prej1 = 0.35502805388781723926006318600418317639797917419917724058332651030081004245
    y_prej2 = -0.25881940379280679840518356018920396347909113835493458221000181385610277267

    y_prej = [y_prej1, y_prej2]

    step = -h

    xs_total = [0.0]
    ys_total = [y_prej[1]]

    i = 0.0
    korak = 0
    while korak <= max_korakov
        y = airy_premik(y_prej, i, step)
        korak = korak + 1

        # najden interval na katerem se nahaja ničla
        if sign(y[1]) != sign(y_prej[1])
            # kličemo eno izmed metod, da dobimo ničlo z želeno natančnostjo
            if metoda == "tangentna"
                n, iter = tangentna(y_prej, i)
            elseif metoda == "regula"
                n, iter = regula_falsi(y, i+step, y_prej, i)
            else
                n, iter = bisekcija(y, i+step, y_prej, i)
            end

            push!(nicle, n)
        end

        push!(xs_total, i)
        push!(ys_total, y[1])

        y_prej = y
        i = i + step

        # če smo našli zahtevano število ničel, zaključimo
        if length(nicle) == k 
            return nicle
        end
    end
    # če smo presegli maksimalno število korakov, vseeno vrnemo vse ničle, ki smo jih našli
    # na tem mestu bi po potrebi lahko vrnili tudi napako
    return nicle
end

"""
najde vse ničle na intervalu [a, 0)
kjer je a negativno število
"""
function airy_nicle_na_intervalu(a::Float64; h::Float64=0.01, metoda::String="bisekcija")
    nicle = []
    xs = []
    ys = []

    # vnaprej točno izračunane vrednosti začetnih pogojev y(0) in y'(0)
    y_prej1 = 0.35502805388781723926006318600418317639797917419917724058332651030081004245
    y_prej2 = -0.25881940379280679840518356018920396347909113835493458221000181385610277267

    y_prej = [y_prej1, y_prej2]

    if a < 0
        step = -h
    else
        step = h
    end

    xs_total = [0.0]
    ys_total = [y_prej[1]]

    i = 0.0
    while abs(i) <= abs(a)
        y = airy_premik(y_prej, i, step)

        # najden interval na katerem se nahaja ničla
        if sign(y[1]) != sign(y_prej[1])
            # kličemo eno izmed metod, da dobimo ničlo z želeno natančnostjo
            if metoda == "tangentna"
                n, iter = tangentna(y_prej, i)
            elseif metoda == "regula"
                n, iter = regula_falsi(y, i+step, y_prej, i)
            else
                n, iter = bisekcija(y, i+step, y_prej, i)
            end

            push!(nicle, n)
        end

        push!(xs_total, i)
        push!(ys_total, y[1])

        y_prej = y
        i = i + step
    end

    return nicle#, xs, ys, xs_total, ys_total
end

function tangentna(y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=1000)
    i = 0
    x = 0
    for i=1:max_iter
        x = x_prej - y_prej[1]/y_prej[2]
        h = x-x_prej

        y = airy_premik(y_prej, x_prej, h)
        
        if abs(y[1]-y_prej[1]) < tol && abs(x-x_prej) < tol
            break
        end

        y_prej = y
        x_prej = x
    end

    return x, i
end

function regula_falsi(y::Vector, x::Float64, y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=1000)
    c = x
    i = 0
    for i=1:max_iter
        c = x_prej - y_prej[1] * (x-x_prej)/(y[1]-y_prej[1])

        y_c = airy_premik(y_prej, x_prej, c-x_prej)
        if y[1] * y_c[1] < 0
            y_prej = y_c
            x_prej = c
        else
            y = y_c
            x = c
        end

        # preverimo ali izračunana točka ustreza želeni toleranci
        if abs(y[1]-y_prej[1]) < tol && abs(x-x_prej) < tol
            break
        end

    end
    return c, i
end

function bisekcija(y::Vector, x::Float64, y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=1000)
    tmp = 1
    xs = [x, x_prej]
    ys = [y[1], y_prej[1]]
    x_m = 0
    i = 0
    for i=1:max_iter
        # izračunamo sredino intervala
        h = (x-x_prej)/2
        x_m = x_prej + h

        # izračunamo vrednost funkcije v srednji točki
        y_m = airy_premik(y_prej, x_prej, h)

        # preverimo na katerem od podintervalov se predznaka razlikujeta
        if sign(y[1]) != sign(y_m[1])
            y_prej = y_m
            x_prej = x_m
            push!(xs, x_m)
            push!(ys, y_m[1])
        else
            y = y_m
            x = x_m
            push!(xs, x_m)
            push!(ys, y_m[1])
        end

        # preverimo ali izračunana točka ustreza želeni toleranci
        if abs(y[1]-y_prej[1]) < tol && abs(x-x_prej) < tol
            break
        end

        tmp = tmp + 1
    end

    return x_m, i
end


end # module Airy
