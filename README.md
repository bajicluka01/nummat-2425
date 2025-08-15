# 18.2.2 Fresnelov integral

Avtor: Luka Bajić (<lb4129@student.uni-lj.si>)

## Opis problema

Repozitorij vsebuje paket za numerično računanje Fresnelovega kosinusa. Paket ponuja funkcijo `fresnel_cos`, pri kateri lahko izbiramo med dvema metodama: adaptivno Simpsonovo pravilo in (manj natančne) Gauss-Laguerrove kvadrature.

## Navodila za uporabo kode

Repozitorij kloniramo z ukazom:

```
git clone https://github.com/bajicluka01/nummat-2425.git
```

in nato izberemo vejo dn2:

```
git checkout dn2
```

Glavni del kode se nahaja v datoteki `src/Fresnel.jl`. Gre za modul `Fresnel`, katerega lahko uvozimo in uporabljamo v ostalih skriptah. Primer uporabe je demonstriran v datoteki `scripts/demo.jl`, ki prav tako služi kreaciji poročila, ki se nahaja v datoteki `report/report.pdf`.

## Navodila za poganjanje testov 

Najprej v korenski mapi projekta odpremo Julia REPL in nato vstopimo v paketni način, kjer poženemo naslednji ukaz:

```
test Fresnel
```

## Navodila za kreiranje poročila

Najprej v korenski mapi projekta odpremo Julia REPL in nato poženemo naslednji ukaz:

```
include("scripts/makedocs.jl")
```
