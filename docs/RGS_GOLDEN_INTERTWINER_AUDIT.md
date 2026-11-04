# RGS–PSU–Q6: audit dell'intertwiner aureo

## Oggetto

Si confrontano due insiemi finiti da 64 punti, senza identificare in anticipo
le etichette dei vertici:

1. il seme ricostruito X_Oh del modello RGS 64B;
2. la proiezione fisica P_parallel({0,1}^6 - 1/2) della costruzione aurea
   6D→3D.

La matrice di proiezione usa phi=(1+sqrt(5))/2 e
phi*=(1-sqrt(5))/2, con i due sottospazi fisico e interno ortonormali
e mutuamente ortogonali.

## Risultato numerico

L'audit Python src.rgs_psu_intertwiner.audit_rgs_golden_intertwiner calcola
covarianze, shell radiali e spettro completo delle 2016 distanze a coppie.
Per i dati canonici:

- covarianza Q6 proiettata: C_Q6 = 0.25 I_3;
- covarianza del seme RGS: C_RGS = 1.25 I_3;
- scala che sarebbe imposta da un'eventuale similitudine: sqrt(5);
- residuo RMS dello spettro delle distanze normalizzato: circa 0.069209259652.

Il confronto delle shell radiali dà molteplicità differenti:
RGS 8, 8, 24, 24; proiezione Q6 12, 20, 20, 12.

## Ostacolo affine: dimostrazione

Siano X,Y sottoinsiemi di R^3 i due cloud centrati, con la stessa
molteplicità uniforme di 64 punti. Se esistesse una corrispondenza affine
biunivoca y=A x+b, dopo il centramento si avrebbe

    C_Y = A C_X A^T.

Poiché C_X=c_X I e C_Y=c_Y I, segue

    A A^T = (c_Y/c_X) I.

Quindi A=s O, con O ortogonale: ogni equivalenza affine sarebbe
necessariamente una similitudine. Ma una similitudine conserva lo spettro
delle distanze a meno di un fattore comune; il residuo normalizzato non nullo
0.069209259652 lo esclude.

**Conclusione:** non esiste un'intertwiner affine biunivoco tra questi due
insiemi di punti, con la definizione canonica e i pesi uniformi adottati.
Non è un fallimento della decomposizione aurea 3+3: è un no-go per
l'identificazione diretta dei due cloud.

## Cosa resta aperto

- una mappa non affine motivata da una metrica/action RGS;
- un'intertwiner fibrico che includa le variabili PSU/Anti-Twin e non soltanto
  le coordinate fisiche 3D;
- una finestra di accettazione o un seed RGS alternativi, purché derivati
  da un principio variazionale e non scelti per adattare i dati;
- il confronto dei gruppi di simmetria come rappresentazioni, non solo come
  insiemi di punti;
- il limite continuo e la convergenza Mosco.

Questo audit non dimostra una teoria fisica unificata, né identifica la
proiezione aurea con il seme RGS. Fornisce un vincolo quantitativo che ogni
costruzione futura deve rispettare.

## Esecuzione

    python -m src.rgs_psu_intertwiner --output results/rgs_golden_intertwiner.json
    pytest -q tests/test_rgs_psu_intertwiner.py
