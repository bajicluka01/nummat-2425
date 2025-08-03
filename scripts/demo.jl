#' # Ničle Airyjeve funkcije
#' Avtor: Luka Bajić

#' ## Opis problema
#' Problem, ki ga rešujemo je sestavljen iz dveh delov: numeričnega reševanja diferencialne enačbe drugega reda
#' in iskanja ničel funkcije z metodami kot so bisekcija in regula falsi. Rezultati morajo biti natančni na deset decimalnih mest.

#' Želimo poiskati čimveč ničel Airyjeve funkcije, ki je dana z naslednjo diferencialno enačbo:

#' $Ai''(x)-xAi(x)=0$

#' pri začetnih pogojih

#' $Ai(0) = \frac{1}{3^{\frac{2}{3}}\Gamma(\frac{2}{3})}, Ai'(0) = -\frac{1}{3^{\frac{1}{3}}\Gamma(\frac{1}{3})}$

#' Vrednosti funkcije lahko računamo z uporabo Magnusove metode, pri kateri se z izbranim korakom $h$ premikamo od začetne vrednosti
#' v levo po $x$-osi (ker vemo, da ima funkcija $Ai$ vse ničle negativne). Premik izvajamo z naslednjo formulo:

#' $y_{k+1} = exp(\frac{h}{2}(A_1+A_2)-\frac{\sqrt{3}}{12}h^2[A_1,A_2])$

#' pri čemer je $A_{1,2}=A(x_k+(\frac{1}{2}\pm \frac{\sqrt{3}}{6})h)$ in $[A_1,A_2] = A_1A_2-A_2A_1$. 
#' Matriko $A$ pa dobimo tako, da zgornjo diferencialno enačbo drugega reda prevedemo na sistem diferencialnih enačb prvega reda,
#' in sicer: 

#' $Y'(x)=\begin{bmatrix}Ai'(x) \\ Ai''(x)\end{bmatrix} = \begin{bmatrix}Ai'(x) \\ xAi'(x)\end{bmatrix}$

#' kar lahko preoblikujemo v ustrezno obliko za Magnusovo metodo:

#' $Y'(x) = \begin{bmatrix}0 & 1 \\ x & 0\end{bmatrix}\begin{bmatrix}Ai(x) \\ Ai'(x)\end{bmatrix}$

#' kjer je iskana matrika torej $A=\begin{bmatrix}0 & 1 \\ x & 0\end{bmatrix}$

#' ## Opis rešitve
#' Za izračun vrednosti Airyjeve funkcije enostavno sledimo zgoraj navedenim formulam, pri čemer lahko določene vrednosti,
#' npr. vrednosti začetnih pogojev, vnaprej izračunamo z orodjem kot je Wolfram Alpha, ker je njihova vrednost neodvisna
#' od ostalih parametrov. Ostali izračuni se izvajajo v metodi $\text{airy\_premik}$, ki na podlagi prejšnjih vrednosti $y_k$, in $y_k'$
#' pri $x_k$, in parametra za korak $h$, ki ga lahko poljubno določimo, izračuna vrednosti funkcije pri $x_k+h$, 
#' torej $y_{k+1}$ in $y_{k+1}'$.

#' ### Iskanje ničel
#' Implementacija ponuja dve možnosti za iskanje ničel: uporabnik specificira interval $[a,0]$, na katerem 
#' metoda $\text{airy\_nicle\_na\_intervalu}$ najde vse ničle, ali pa specificira vrednost $k$, nakar 
#' metoda $\text{airy\_k\_nicel}$ vrne prvih $k$ ničel od koordinatnega izhodišča proti $-\infty$.

#' Razlika med metodama je samo v zaustavitvenem pogoju, postopek iskanja posamezne ničle ostaja enak, in sicer:
#' z zgoraj omenjenim postopkom se s korakom $h$ premikamo po funkciji in na vsakem koraku preverjamo ali se predznak vrednosti
#' funkcije razlikuje od predznaka vrednosti funkcije na prejšnjem koraku. Ker vemo, da je funkcija zvezna, razlika med prezdnakoma pomeni, 
#' da se na intervalu med $x_k$ in $x_{k+1}$ nahaja ničla. Ker je zahtevana natančnost na deset decimalk ($h$ pa je tipično bistveno večji 
#' in posledično ni dovolj natančen), se na tem mestu izvede ena izmed metod regula falsi, bisekcija ali tangentna metoda 
#' (izbiro podamo kot argument) za iskanje ničle.

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
