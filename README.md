# 18.3.1 Ničle Airyjeve funkcije

Avtor: Luka Bajić (<lb4129@student.uni-lj.si>)

## Opis problema

Repozitorij vsebuje paket za iskanje ničel Airyjeve funkcije. Paket ponuja dve funkciji: `airy_k_nicel`, ki najde prvih k ničel, začenši v koordinatnem izhodišču proti $-\infty$, in `airy_nicle_na_intervalu`, ki najde vse ničle na intervalu $[a,0]$. Za metodo za iskanje ničel lahko izbiramo med bisekcijo, regula falsi in tangentno metodo. 

## Navodila za uporabo kode

Repozitorij kloniramo z ukazom:

```
git clone https://github.com/bajicluka01/nummat-2425.git
```

in nato izberemo vejo dn3:

```
git checkout dn3
```

Glavni del kode se nahaja v datoteki `src/Airy.jl`. Gre za modul `Airy`, katerega lahko uvozimo in uporabljamo v ostalih skriptah. Primer uporabe je demonstriran v datoteki `scripts/demo.jl`, ki prav tako služi kreaciji poročila, ki se nahaja v datoteki `report/report.pdf`.

## Navodila za poganjanje testov 

Najprej v korenski mapi projekta odpremo Julia REPL in nato vstopimo v paketni način, kjer poženemo naslednji ukaz:

```
test Airy
```

## Navodila za kreiranje poročila

Najprej v korenski mapi projekta odpremo Julia REPL in nato poženemo naslednji ukaz:

```
include("scripts/makedocs.jl")
```
