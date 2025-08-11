# 18.1.3 Metoda konjugiranih gradientov s predpogojevanjem

Avtor: Luka Bajić (<lb4129@student.uni-lj.si>)

## Opis problema

Repozitorij vsebuje paket za reševanje sistemov $Ax=b$, kjer je $A$ simetrična in pozitivno definitna, z uporabo metode konjugiranih gradientov. Izberemo lahko med metodo s predpogojevanjem ali brez, pri čemer je za predpogojevanje implementirana tudi pomožna funkcija, ki izračuna nepopolen razcep Choleskega.

## Navodila za uporabo kode

Repozitorij kloniramo z ukazom:

```
git clone https://github.com/bajicluka01/nummat-2425.git
```

in nato izberemo vejo dn1:

```
git checkout dn1
```

Glavni del kode se nahaja v datoteki `src/MKG.jl`. Gre za modul `MKG`, katerega lahko uvozimo in uporabljamo v ostalih skriptah. Primer uporabe je demonstriran v datoteki `scripts/demo.jl`, ki prav tako služi kreaciji poročila, ki se nahaja v datoteki `report/report.pdf`.

## Navodila za poganjanje testov 

Najprej v korenski mapi projekta odpremo Julia REPL in nato vstopimo v paketni način, kjer poženemo naslednji ukaz:

```
test MKG
```

## Navodila za kreiranje poročila

Najprej v korenski mapi projekta odpremo Julia REPL in nato poženemo naslednji ukaz:

```
include("scripts/makedocs.jl")
```
