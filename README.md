# Betting-Against-Beta
Kode og data til analyse og backtest i afhandlingen om Betting Against Beta

** Backtesting analyse kode **
-
I denne fil findes alt kode vedrørende de 49 industriporteføljer. I koden findes der 
- Beregning af BAB & SDBAB
- Beregning af P1-P10
- Performance variable for BAB, SDBAB, P1-P10
- Proposition 3
- Beregning af Transaktionsomkostninger
- Beregning af Bear/Bull marked
- Følsomhedsanalyse
- Diverse plots

Koden er sat op til at kunne køre ved først at uploade Market-filen og derefter FFM3-filen.
Koden er opsat til at køre efter en specifik konfiguration, følgende parametre skal ændres manuelt, hvis anden konfiguration ønskes testet

Log / ikke log
- Linje 55 & 58

Rullende afkast
- Linje 79

Shrinkage vægt
- Linje 140

Stokastisk dominans (Orden eller trim)
- Linje 218 & 363

Stokastisk dominans (dominansperiode)
- Linje 236 & 247

Stokastisk dominans (signifikansniveau)
- Linje 274 & 367

Transaktionsomkostninger (Låne- og transaktionsomkostinger)
- Linje 542 & 543

P1-P10 (Log / ikke log)
- Linje 1491 & 1492

Ved ændring i nogle af ovenstående parametre kan der forekomme ændringer i længderne af forskellige vektorere, som skal ændres manuelt.

** Empirisk analyse kode **
-
I denne fil findes alt kode vedrørende det empiriske data. I koden findes der 
- Beregning af BAB & SDBAB
- Beregning af P1-P10
- Performance variable for BAB, SDBAB, P1-P10
- Proposition 3
- Beregning af Transaktionsomkostninger
- Følsomhedsanalyse
- Diverse plots

Koden er sat op til at kunne køre ved først at uploade Market-filen og derefter FFM3-filen.
Koden er opsat til at køre efter en specifik konfiguration, følgende parametre skal ændres manuelt, hvis anden konfiguration ønskes testet

Log / ikke log
- Linje 46 & 49

Rullende afkast
- Linje 71

Shrinkage vægt
- Linje 131

Stokastisk dominans (Orden eller trim)
- Linje 212 & 334

Stokastisk dominans (dominansperiode)
- Linje 226 & 235

Stokastisk dominans (signifikansniveau)
- Linje 259 & 337

Transaktionsomkostninger (Låne- og transaktionsomkostinger)
- Linje 511 & 512

P1-P10 (Log / ikke log)
- Linje 953 & 954

Ved ændring i nogle af ovenstående parametre kan der forekomme ændringer i længderne af forskellige vektorere, som skal ændres manuelt.

** Kumulativ afkast - shrinkage faktor **
-
I denne fil kan man få et overblik over BAB-strategiens kumulative afkast med forskellige shrinkage faktorer 

