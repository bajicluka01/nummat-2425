module Airy
using LinearAlgebra
using Plots

export airy_k_nicel, airy_nicle_na_intervalu

"""
    y = airy_premik(y_prej, x_prej, h)

Izračunaj vrednost Airyjeve funkcije za `x_prej`+`h` iz prejšnje vrednosti `y_prej`.
"""
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
    nicle = airy_k_nicel(k, h)

Izračunaj prvih `k` ničel Airyjeve funkcije, s premiki s korakom `h`.
"""
function airy_k_nicel(k::Integer; h::Float64=0.01, metoda::String="bisekcija", max_korakov::Integer=10^12, vrni_vmesne::Bool=false, vrni_povp_iter::Bool=false)
    nicle = []
    y_prej1 = 0.35502805388781723926006318600418317639797917419917724058332651030081004245
    y_prej2 = -0.25881940379280679840518356018920396347909113835493458221000181385610277267
    y_prej = [y_prej1, y_prej2]
    step = -h
    xs_total = [0.0]
    ys_total = [y_prej[1]]
    i = 0.0
    korak = 0
    total_iter = 0
    while korak <= max_korakov
        y = airy_premik(y_prej, i, step)
        korak = korak + 1

        # najden interval na katerem se nahaja ničla
        if sign(y[1]) != sign(y_prej[1])
            # kličemo eno izmed metod, da dobimo ničlo z želeno natančnostjo
            if metoda == "bisekcija"
                n, iter = bisekcija(y, i+step, y_prej, i)
            elseif metoda == "regula"
                n, iter = regula_falsi(y, i+step, y_prej, i)
            else # privzeta metoda je tangentna
                n, iter = tangentna(y_prej, i)
            end

            push!(nicle, n)
            total_iter = total_iter + iter
        end

        push!(xs_total, i)
        push!(ys_total, y[1])

        y_prej = y
        i = i + step

        # če smo našli zahtevano število ničel, zaključimo
        if length(nicle) == k 
            if vrni_vmesne
                if vrni_povp_iter
                    return nicle, xs_total, ys_total, total_iter/length(nicle)
                end
                return nicle, xs_total, ys_total
            else
                if vrni_povp_iter
                    return nicle, total_iter/length(nicle)
                end
                return nicle
            end
        end
    end
    # če smo presegli maksimalno število korakov, vseeno vrnemo vse ničle, ki smo jih našli
    # na tem mestu bi po potrebi lahko vrnili tudi napako
    if vrni_vmesne
        if vrni_povp_iter
            return nicle, xs_total, ys_total, total_iter/length(nicle)
        end
        return nicle, xs_total, ys_total
    else
        if vrni_povp_iter
            return nicle, total_iter/length(nicle)
        end
        return nicle
    end
end

"""
    nicle = airy_nicle_na_intervalu(a, h)

Izračunaj vse ničle Airyjeve funkcije na intervalu [`a`, 0], s premiki s korakom `h`.
"""
function airy_nicle_na_intervalu(a::Float64; h::Float64=0.01, metoda::String="bisekcija", vrni_vmesne::Bool=false)
    nicle = []

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
            if metoda == "bisekcija"
                n, iter = bisekcija(y, i+step, y_prej, i)
            elseif metoda == "regula"
                n, iter = regula_falsi(y, i+step, y_prej, i)
            else # privzeta metoda je tangentna
                n, iter = tangentna(y_prej, i)
            end

            push!(nicle, n)
        end

        push!(xs_total, i)
        push!(ys_total, y[1])

        y_prej = y
        i = i + step
    end

    if vrni_vmesne
        return nicle, xs_total, ys_total
    end
    return nicle
end

"""
    x, iter = tangentna(y_prej, x_prej)

Izračunaj ničlo s tangentno metodo iz točke `x_prej` in vrednosti `y_prej`.
"""
function tangentna(y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=100)
    iter = 0
    x = 0
    for _=1:max_iter
        x = x_prej - y_prej[1]/y_prej[2]
        h = x-x_prej

        y = airy_premik(y_prej, x_prej, h)
        
        iter = iter + 1

        if abs(y[1]-y_prej[1]) < tol && abs(x-x_prej) < tol
            break
        end

        y_prej = y
        x_prej = x
    end

    return x, iter
end

"""
    c, iter = regula_falsi(y, x, y_prej, x_prej)

Izračunaj ničlo z metodo regula falsi na intervalu [`x`, `x_prej`], kjer sta robni vrednosti [`y` in `y_prej`].
"""
function regula_falsi(y::Vector, x::Float64, y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=100)
    c = x
    iter = 0
    for _=1:max_iter
        c = x_prej - (y_prej[1]*(x-x_prej))/(y[1]-y_prej[1])
        y_c = airy_premik(y_prej, x_prej, c-x_prej)

        if y[1] * y_c[1] > 0
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

        iter = iter + 1
    end
    return c, iter
end

"""
    x_m, iter = bisekcija(y, x, y_prej, x_prej)

Izračunaj ničlo z bisekcijo na intervalu [`x`, `x_prej`], kjer sta robni vrednosti [`y` in `y_prej`].
"""
function bisekcija(y::Vector, x::Float64, y_prej::Vector, x_prej::Float64; tol::Float64=10e-11, max_iter::Integer=100)
    xs = [x, x_prej]
    ys = [y[1], y_prej[1]]
    x_m = 0
    iter = 0
    for _=1:max_iter
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

        iter = iter + 1

        # preverimo ali izračunana točka ustreza želeni toleranci
        if abs(y[1]-y_prej[1]) < tol && abs(x-x_prej) < tol
            break
        end
    end

    return x_m, iter
end

end # module Airy
