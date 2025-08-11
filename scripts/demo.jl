#' # Metoda konjugiranih gradientov s predpogojevanjem
#' Avtor: Luka Bajić

#' ## Opis problema
#' Metoda konjugiranih gradientov je postopek za reševanje linearnega sistema enačb $Ax = b$,
#' ob predpostavki, da je matrika A pozitivno definitna. V našem primeru se ukvarjamo z razpršenimi matrikami,
#' torej matrikami, ki imajo večino elementov ničelnih. Da izboljšamo hitrost konvergence metode konjugiranih
#' gradientov, uporabimo predpogojevanje z matriko $M=LL^T$, kjer spodnjetrikotno matriko $L$ dobimo z nepopolnim razcepom Choleskega. 

#' Ker imamo opravka z razpršenimi matrikami, je s stališča prostorske zahtevnosti smiselno, da jih hranimo
#' v posebni strukturi, ki hrani samo neničelne elemente in sestoji iz treh vektorjev: vektor 
#' vrednosti v matriki, vektor indeksov po vrsticah in vektor indeksov po stolpcih. Na primer za sledečo matriko:

#' $A=\begin{bmatrix}7 & 1.1 & 0 & 0 & 0 \\ 1.1 & 2 & 0 & 0 & 0 \\ 0 & 0 & 3 & 0 & 3 \\ 0 & 0 & 0 & 0.5 & 0 \\ 0 & 0 & 3 & 0 & 4.2 \end{bmatrix}$

#' hranimo vektorje

#' $I = \begin{bmatrix}1 & 1 & 2 & 2 & 3 & 3 & 4 & 5 & 5 \end{bmatrix}$
#' $J = \begin{bmatrix}1 & 2 & 1 & 2 & 3 & 5 & 4 & 3 & 5 \end{bmatrix}$
#' $V = \begin{bmatrix}7 & 1.1 & 1.1 & 2 & 3 & 3 & 0.5 & 3 & 4.2 \end{bmatrix}$

#' Uporabnost tega pristopa postane bolj očitna pri bistveno večjih matrikah in če je število ničelnih elementov dovolj veliko.

#' ## Opis rešitve
#' Implementacija ponuja tri metode: $\textit{nep\_chol}$, ki izračuna nepopolni razcep Choleskega
#' za podano matriko $A$, $\textit{conj\_grad\_baseline}$, ki izvede metodo konjugiranih gradientov
#' nad podano matriko $A$ in vektorjem desnih strani $b$ in vrne iskan $x$, in $\textit{conj\_grad}$, 
#' ki prav tako vrne $x$, le da metodo konjugiranih gradientov izvede s predpogojevanjem s podano 
#' spodnjetrikotno matriko $L$.

#' ### Nepopolni razcep Choleskega
#' 

#' ### Metoda konjugiranih gradientov brez predpogojevanja
#' Iščemo rešitev sistema $Ax=b$, kjer je matrika $A$ simetrična in pozitivno definitna. Za začetni približek $x_0$
#' vzamemo kar ničelni vektor in izračunamo $r_0 = b - Ax_0$. Nastavimo še $p_0=x_0$ in nato ponavljamo sledečo iteracijo:

#' $\alpha_k=\frac{r_k^Tr_k}{p_k^TAp_k}$

#' $x_{k+1}=x_k-\alpha_kp_k$

#' $r_{k+1}=r_k-\alpha_kAp_k$

#' $\beta_k=\frac{r_{k+1}^Tr_{k+1}}{r_k^Tr_k}$

#' $p_{k+1}=r_{k+1}+\beta_kp_k$

#' Postopek zaključimo kadar je norma vektorja $r_{k+1}$ manjša od želene tolerance.

#' ### Metoda konjugiranih gradientov s predpogojevanjem
#' Predpogojevanje izvajamo z matriko $M=LL^T$, kjer spodnjetrikotno matriko $L$ dobimo z nepopolnim razcepom Choleskega.
#' Rešiti moramo sistem $Mz=r$, vendar ker nočemo računati inverza $M^{-1}$ (iz vidika časovne zahtevnosti), lahko
#' najprej s premo substitucijo rešimo sistem $La=r$ in nato z obratno substitucijo rešimo sistem $L^Tz=a$. 
#' Tako prema kot obratna substitucija imata časovno zahtevnost $O(n^2)$, kar je bistveno hitrejše od računanja inverza.

#' Zgornji postopek moramo ponoviti v vsakem koraku iteracije, in sicer na sledeč način: najprej kot začetno vrednost $p_0$
#' nastavimo rešitev sistema $Mz_0=r_0$, korake pa izvajamo na sledeč način:

#' $\alpha_k=\frac{r_k^Tz_k}{p_k^TAp_k}$

#' $x_{k+1}=x_k-\alpha_kp_k$

#' $r_{k+1}=r_k-\alpha_kAp_k$

#' Rešimo sistem $Mz_{k+1}=r_{k+1}$.

#' $\beta_k=\frac{r_{k+1}^Tz_{k+1}}{r_k^Tz_k}$

#' $p_{k+1}=z_{k+1}+\beta_kp_k$

#' Zaustavitveni pogoj ostane enak kot pri osnovni metodi.

#' ## Primer uporabe
#' Spodnji izsek kode demonstrira razliko med metodo konjugiranih gradientov brez predpogojevanja in
#' s predpogojevanjem. 

using MKG, Plots, LinearAlgebra, SparseArrays

I = [1., 1, 2, 2, 3, 3, 4, 5, 5]
J = [1., 2, 1, 2, 3, 5, 4, 3, 5]
V = [7., 1.1, 1.1, 2, 3, 3, 0.5, 3, 4.2]
A = sparse(I, J, V)
b = [2, 3, -5, 1, 0.2]
L = nep_chol(A)
x1, it1, res1 = conj_grad_baseline(A, b, vrniresid=true)
x2, it2, res2 = conj_grad(A, b, L, vrniresid=true)

println("MKG brez predpogojevanja se zaustavi po $it1 korakih.")
println("MKG s predpogojevanjem se zaustavi po $it2 korakih.")

#' Kot vidimo že na relativno majhnem primeru, predpogojevanje bistveno izboljša konvergenco.

#' V spodnjem izseku si ogledamo še primer z bistveno večjo matriko, ki jo zgeneriramo psevdonaključno:
#' za vse diagonalne elemente zgeneriramo naključna števila iz nekega intervala, za preostale elemente pa 
#' z neko manjšo verjetnostjo (npr. 0.1) zgeneriramo manjša števila, da ohranimo dominantnost diagonale.
#' Vsako izmed teh števil zapišemo tako na indeks $(i,j)$ kot $(j,i)$, da ohranimo simetričnost. Ta postopek sicer ne
#' zagotavlja pozitivne definitnosti zgenerirane matrike, vendar se empirično izkaže kot dovolj dober za potrebe te demonstracije.

#' Z zastavico $\textit{vrniresid}$ nam metodi vrneta tabelo residualov (oziroma neskončnih norm vektorjev $r_{k+1}$)
#' na posameznem koraku iteracije. Residuale izrišemo s paketom Plots in tako dobimo vizualno primerjavo hitrosti konvergence
#' med metodama. 

n=700
I = [1.0]
J = [1.0]
V = [rand()*10.0]
for i=2:n
    push!(V, rand()*80.0)
    push!(I, i)
    push!(J, i)

    if rand() < 0.1
        r = rand()
        push!(V, r)
        push!(V, r)
        r1 = rand(1:n)
        r2 = rand(1:n)
        while r1==i || r2==i # poskrbimo, da ne prepišemo diagonale
            r1 = rand(1:n)
            r2 = rand(1:n)
        end
        push!(I, r1)
        push!(J, r2)
        push!(I, r2)
        push!(J, r1)
    end
end
A = sparse(I, J, V)
b = rand(n)
x1, it1, res1 = conj_grad_baseline(A, b, vrniresid=true, tol=10e-20)
L = nep_chol(A)
x2, it2, res2 = conj_grad(A, b, L, vrniresid=true, tol=10e-20)

plot(res1, label="brez predpogojevanja", title="Primerjava residualov")
plot!(res2, label="s predpogojevanjem", title="Primerjava residualov")





n=600
I = [1.0]
J = [1.0]
V = [rand()*10.0]
for i=2:n
    push!(V, rand()*80.0)
    push!(I, i)
    push!(J, i)

    if rand() < 0.1
        r = rand()
        push!(V, r)
        push!(V, r)
        r1 = rand(1:n)
        r2 = rand(1:n)
        while r1==i || r2==i # poskrbimo, da ne prepišemo diagonale
            r1 = rand(1:n)
            r2 = rand(1:n)
        end
        push!(I, r1)
        push!(J, r2)
        push!(I, r2)
        push!(J, r1)
    end
end
A = sparse(I, J, V)
b = rand(n)
x1, it1, res1 = conj_grad_baseline(A, b, vrniresid=true, tol=10e-20)
L2 = nep_chol(A)
L = test(A)
x2, it2, res2 = conj_grad(A, b, L, vrniresid=true, tol=10e-20)
isapprox(L, L2)

tmp = L2-L
for i=1:600
    for j=1:600
        if tmp[i,j] != 0
            println(tmp[i,j])
        end
    end
end

function test(A::AbstractSparseMatrix)
    L = copy(A)
    
    for i=1:size(A,2)
        L = dropzeros(L)
        I, J, V = findnz(L)
        sum = 0
        for k=1:length(J)
            if I[k] == i && J[k] < i #&& V[k] != 0
                println(i, " ", k, " " , sum, " " , J[k], " " , V[k])
                sum = sum + V[k]^2
            end
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

    for i = 1:size(A,2)
        for j = i+1:size(A,2)
            L[i,j] = 0
        end
    end

    return dropzeros(L)
end

A
I, J, V = findnz(A)
L2 = test(A)
L = nep_chol(A)
isapprox(L, L2, atol=10e-10)
l2 = L2*L2'
l = L*L'
isapprox(l2, l)



I = [1., 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
J = [1., 1, 2, 2, 3, 1, 3, 4, 1, 4, 5]
V = [5., -2, 5, -2, 5, -2, -2, 5, -2, -2, 5]
A = sparse(I, J, V)
L2 = test(A)
L = nep_chol(A)

isapprox(L2*L2', res, atol=10e-10)
isapprox(L*L', res, atol=10e-10)

I_res = [1., 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5]
J_res = [1., 2, 4, 5, 1, 2, 3, 4, 5, 2, 3, 4, 1, 2, 3, 4, 5, 1, 2, 4, 5]
V_res = [5, -2, -2, -2, -2, 5, -2, 0.8, 0.8, -2, 5, -2, -2, 0.8, -2, 5, -2, -2, 0.8, -2, 5]
res = sparse(I_res, J_res, V_res, 5, 5)






for j=i+1:size(A,2)
            # elemente nad diagonalo nastavimo na 0
            L[i,j] = 0

            sum = 0
            for k=1:length(J)
                if I[k] == i && J[k] < i 
                    if I[k] == j && J[k] < j
                    #if I[j] == j && J[j] < j
                        #sum = sum + V[k]*V[j]
                        A[j,i] = A[j,i] - V[k]*V[k]
                    end
                end
            end

            if L[j,i] != 0
                #L[j,i] = (A[j,i] - sum)/L[i,i]
                L[j,i] = A[j,i] / L[i,i]
            end
        end
