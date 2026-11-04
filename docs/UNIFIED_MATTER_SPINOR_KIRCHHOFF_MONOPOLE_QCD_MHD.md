# Materia, correnti spinoriali, monopoli, quark e MHD — raccordo unitario RGS–PSU

**Inserto della monografia canonica in 27 capitoli.** Non è un nuovo capitolo numerato. Collega i capitoli già esistenti su geometria RGS/TEGR, Q6–Anti-Twin, 2I, Clifford/Dirac, elettromagnetismo, QCD/quark, superfluidi e turbolenza. Le Fasi del repository sono milestone di verifica, non nuovi capitoli.

## 1. Azione unitaria come obiettivo

La teoria deve derivare geometria, campi di materia, correnti e fluidi da un'unica azione, mantenendo distinti i gruppi di gauge:
[
S_{m tot}=S_{m TEGR}[e,omega]+S_D[e,omega,A,W,G,psi]
+S_{m Maxwell}[A]+S_{m YM}[W,G]+S_{m fluid}
+S_{m int}.
]
Qui (e^a{}_mu) è il coframe, (omega^a{}_{bmu}) la connessione di spin, (A_mu) il campo U(1) elettromagnetico, (W_mu) e (G_mu) i campi non-abeliani nei settori in cui la rappresentazione li richiede, e (psi) i multipleti fermionici. I campi bosonici/condensati hanno termini cinetici e statistiche proprie: non sono identificati con i fermioni.

[
D_mupsi=left(partial_mu+rac14omega_{abmu}gamma^{ab}
-iqA_mu-ig_2W_mu^It^I-ig_3G_mu^AT^Aight)psi,
quad
S_D=int d^4x,eleft[rac i2arpsigamma^ae_a{}^mu
overleftrightarrow D_mupsi-marpsipsiight].
]
Il current elettromagnetico deriva dalla variazione rispetto ad (A_mu):
(j^mu_{m EM}=q e_a{}^muarpsigamma^apsi), sommato sulle specie con le cariche effettive.

## 2. Correnti spinoriali come correnti di Kirchhoff su Q6

La nuova implementazione del repository definisce sul Q6 canonico (64 vertici, 192 archi, 240 quadrati) uno spinore (psi_vinmathbb C^2) per vertice e link (U_{ij}in SU(2)) sugli archi orientati (i<j). Le trasformazioni locali sono
[
psi_imapsto G_ipsi_i,qquad U_{ij}mapsto G_jU_{ij}G_i^dagger.
]
Con
[
(D_Upsi)_{ij}=psi_j-U_{ij}psi_i,qquad H_U=D_U^dagger W D_U,
]
il flusso orientato è
[
J_{ij}=2w_{ij}operatorname{Im}(psi_i^dagger U_{ij}^daggerpsi_j).
]
Con la convenzione d'incidenza del repository (B_1[i,e]=-1), (B_1[j,e]=+1), l'equazione (idotpsi=H_Upsi) implica esattamente
[
oxed{dotho=B_1J,qquad ho_i=psi_i^daggerpsi_i.}
]
Questa è la legge di Kirchhoff discreta per la probabilità spinoriale: variazione della densità al nodo uguale al flusso netto entrante. Il current è invariante per trasformazioni locali; il trasporto su cammini aperti trasforma agli estremi e l'olonomia chiusa per coniugazione alla base.

**Distinzione obbligatoria:** il current implementato è un flusso di probabilità SU(2) su un grafo finito. Non è ancora un current elettrico in ampere. Per ottenere quello fisico occorre includere la carica U(1) corretta nel multipletto e derivare il current dalla stessa azione accoppiata. U(1) elettromagnetico, SU(2) di spin/gauge e SU(3) di colore non sono intercambiabili.

File aggiunti: `src/spinor_kirchhoff.py`, `tests/test_spinor_kirchhoff.py`, `examples/spinor_kirchhoff_audit.py`. Il JSON prodotto registra seed, tolleranze, SHA-256 degli input canonici, residui di covarianza dei cammini aperti, coniugazione delle olonomie chiuse, invarianza del current, hermiticità, intertwining del Laplaciano e continuità di Kirchhoff. Riproduzione:
```bash
PYTHONPATH=. python examples/spinor_kirchhoff_audit.py --seed 20261104 --gauge-seed 20261105
```

## 3. Monopolo di Dirac: U(1), fibrato e carica

Il monopolo si descrive con potenziali locali su due carte di (S^2), incollati da una funzione di transizione U(1). La quantizzazione è l'integralità del primo numero di Chern:
[
n=rac{q}{2pihbar}int_{S^2}Finmathbb Z.
]
In convenzioni gaussiane equivale a (qg=nhbar/2); in unità razionalizzate i fattori (4pi) cambiano. La convenzione va fissata prima di usare cariche numeriche.

Per il ponte PSU bisogna: (i) scegliere un 2-ciclo chiuso e orientato che circondi il difetto; (ii) definire fasi U(1) sui link e flussi sulle facce; (iii) costruire le due carte e le transizioni; (iv) calcolare il Chern number e verificarne l'indipendenza dalla posizione della stringa; (v) accoppiare (A_mu) al current carico.

Il complesso Möbius 12/18/6 è non orientabile e non è, da solo, la superficie (S^2) orientata richiesta dall'integrale standard di Chern. Serve una sfera discreta separata (per esempio un guscio icosaedrico) o un 2-ciclo chiuso orientato verificato nel complesso raffinato. Non si deduce carica magnetica dalla sola non-orientabilità. Il monopolo rimane un settore teorico condizionato, non una particella osservata.

## 4. Il ponte dei 120 gradi: 2I, A5 e tripletto SU(3)

Tenere separati tre fatti:
- il gruppo rotazionale icosaedrico è (A_5), ordine 60;
- il suo doppio rivestimento binario (2Isubset SU(2)) ha ordine 120;
- l'icosaedro ha rotazioni di 120° attorno agli assi ternari. L'ordine 120 del gruppo non è un angolo di 120°; il lift spinoriale usa il semiangolo.

I tre pesi della rappresentazione fondamentale di SU(3), disegnati nel piano di Cartan, formano un triangolo equilatero con separazioni angolari di 120°. Questo offre una motivazione geometrica per confrontare il settore ternario icosaedrico con un tripletto di colore, ma l'angolo non dimostra equivalenza di rappresentazioni.

Il test effettivo è costruire un embedding specifico di (A_5) e calcolare
[
mathrm{Hom}_{A_5}(V_{A_5},mathbb C^3)
={X:ho_3(g)X=Xho_{A_5}(g) orall gin A_5}.
]
Verificare dimensione, unicità, unitarietà, commutatori e compatibilità con i generatori di carica. Il repository ha già un report a 120 direzioni/Hart e moduli per (2I), rappresentazioni e intertwiners; il report delle 120 direzioni non va confuso con le rotazioni di 120°.

## 5. Cariche frazionarie: il settore elettrodebole

Con la convenzione (Q=T_3+Y), il doppietto sinistro di quark ha (Y=1/6):
[
Q(u_L)=+rac12+rac16=+rac23,qquad
Q(d_L)=-rac12+rac16=-rac13.
]
I singoletti destri hanno (T_3=0), (Y(u_R)=2/3), (Y(d_R)=-1/3). Con la convenzione (Q=T_3+Y/2), gli ipercarichi si raddoppiano. Il down quark ha carica (-1/3), non (+1/3); (1/3) può essere usato come modulo solo dichiarandolo.

Il materiale RGS archiviato contiene la matrice attivo/sterile
[
V_{m biverso}=phi^{-2}operatorname{diag}(2/3,2/3,1/3)
]
e il tensore simmetrico traceless del tripletto (T_{2u}),
[
H_{m aniso}=egin{pmatrix}
-0.0191&-0.0200&0.0164\
-0.0200&0.0139&0.0068\
0.0164&0.0068&0.0052
end{pmatrix},quad |H_{m aniso}|_Fsimeq0.0449.
]
Sono ingredienti candidati per l'intertwiner/charge-sector audit, ma (V_{m biverso}) non è la matrice delle cariche elettriche: agisce in uno spazio con definizione diversa. Per dimostrare che il modello ricava (2/3) e (-1/3), occorre costruire la mappa verso (SU(3)_c\times SU(2)_L\times U(1)_Y), recuperare il generatore (Q), derivare le rappresentazioni chirali e verificare cancellazione delle anomalie. Il tripletto (T_{2u}) non è automaticamente il tripletto fondamentale di colore.

## 6. MHD, superfluidi e settori bosonico/fermionico

Il campo bosonico condensato è descritto da un parametro d'ordine complesso e dai suoi modi collettivi; la materia fermionica da spinori e statistiche di Fermi. In (^{3})He l'accoppiamento attrattivo forma coppie di Cooper, che hanno comportamento collettivo bosonico. Questo collega i due settori senza identificarli.

Distinguere i sistemi:
- il grande vortice in (^{4})He è un superfluido neutro e genera una geometria acustica efficace per le eccitazioni; non è di per sé un esperimento di MHD elettromagnetica;
- per He II la descrizione intermedia appropriata è il modello a due fluidi/HVBK, con componente normale, superfluida, vortici quantizzati e mutual friction;
- MHD classica descrive fluidi conduttori/plasmi carichi e richiede Maxwell + legge di Ohm + equazione della quantità di moto;
- condensati carichi/superconduttori richiedono un ordine parametrico accoppiato al campo U(1).

Per un fluido conduttore, in MHD resistiva:
[
partial_tmathbf B=
abla	imes(mathbf v	imesmathbf B)+eta_m
abla^2mathbf B,
qquad 
ablacdotmathbf B=0.
]
Per He II si parte invece dalle equazioni a due fluidi/HVBK. Il ponte comune è la struttura variazionale e di conservazione: correnti di materia sorgono nei campi di gauge; stress di gauge e forze di Lorentz entrano nella quantità di moto del fluido; il coframe/connessione fissa i frame locali. Le leggi costitutive restano specifiche del sistema.

## 7. Sequenza e gate della teoria unificata

[
e	oomega	o T^a,R^{ab}	o D_{m Clifford}
	o j^mu_{m spinor}	o(A_mu,F_{mu
u},	ext{monopolo})
	o(	ext{Yang--Mills/QCD},	ext{MHD/HVBK})	o	ext{coarse graining}.
]
È una sequenza di dipendenze, non una lista di risultati già dimostrati. I gate sono: recuperare (Hinmathbb R^{192\times192}); chiarire (K_0,K_{tangent},K_{geometry,tangent}) e inverse/pseudoinverse; canonizzare incidenza Q6 e mappa PRISM; derivare connessione metricamente compatibile dal coframe; fissare Lorentz/Spin/Pin; costruire Dirac e i suoi currents; verificare Chern number U(1); calcolare gli intertwiners (2I/A_5)–SU(3) e le cariche; derivare HVBK o MHD secondo il sistema; infine dimostrare equidistribuzione, measured GH e Mosco (liminf + recovery) sulle raffinazioni. Solo dopo si può dichiarare un limite continuo TEGR+Dirac+Maxwell/Yang–Mills+fluidi.

## 8. Riferimenti esterni essenziali

- R. Heras, “Dirac quantisation condition: a comprehensive review”, *Contemporary Physics* 59 (2018), [doi:10.1080/00107514.2018.1527974](https://doi.org/10.1080/00107514.2018.1527974): patch, stringa di Dirac e quantizzazione.
- “Monopole charge quantization or why electromagnetism is a U(1)-gauge theory”, *Annals of Physics* 312 (2004), [doi:10.1016/j.aop.2004.02.001](https://doi.org/10.1016/j.aop.2004.02.001): fibrato U(1) e primo Chern.
- “Clifford Spinors and Root System Induction: H4 and the Grand Antiprism”, *Advances in Applied Clifford Algebras* (2021), [doi:10.1007/s00006-021-01139-2](https://doi.org/10.1007/s00006-021-01139-2): 2I, H4, Spin/Pin.
- Skrbek, Schmoranzer & Sreenivasan, “Phenomenology of transition to quantum turbulence in flows of superfluid helium”, *PNAS* 121 (2024), [doi:10.1073/pnas.2302256121](https://doi.org/10.1073/pnas.2302256121): quantized circulation, pairing fermionico e due fluidi.
- Dubrulle & Gibbon, *Phil. Trans. R. Soc. A* 380 (2022), [doi:10.1098/rsta.2021.0092](https://doi.org/10.1098/rsta.2021.0092), e Gibbon & Vincenzi, *J. Fluid Mech.* (2026): connessioni matematiche multifrattale–Navier–Stokes.
- D. Tong, *The Standard Model* lecture notes: generatore elettrodebole (Q=T_3+Y), con convenzione degli ipercarichi esplicita.

## 9. Registro epistemico

**Identità finite:** (B_1B_2=0); covarianza dei cammini e coniugazione delle olonomie; continuità spinoriale Q6 (dotho=B_1J); integrità del Chern number per un U(1)-bundle su una superficie chiusa orientata.

**Fisica consolidata:** cariche elettrodeboli dei quark, tripletto di colore SU(3), corrente spinoriale da azione gauge-invariante, quantizzazione di Dirac condizionata al monopolo, due fluidi per He II e MHD per fluidi conduttori.

**Ipotesi RGS da testare:** mappa T2u/2I verso tripletto SU(3), derivazione delle cariche dalla matrice attivo/sterile, connessione fra current Q6 e corrente elettrica, azione comune PSU–spinore–Maxwell–Yang–Mills–fluidi.

**Aperti:** monopoli fisici, mappa di rappresentazione e cancellazione anomalie, calibrazione SI dei currents, azione comune, limite Lorentziano/TEGR/Dirac e convergenza GH/Mosco.
