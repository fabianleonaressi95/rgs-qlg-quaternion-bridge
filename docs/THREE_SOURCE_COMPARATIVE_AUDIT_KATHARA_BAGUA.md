# Comparative source audit: Krystal Spiral, Fibonacci/Krystal claims, and Bagua trigram conventions

**Audit date:** 2026-11-09.  
**Method:** complete readable article text and embedded diagram descriptions were inspected for the three user-supplied pages. The Kathara 12-center/15-line transcription was cross-checked against the Dance for Life, Love, Freedom manual diagram (2002). This is a source/content audit, not endorsement of the sources' metaphysical claims.

## 1. Emerald Awakening — “Spirals of Creation”

URL: https://emerald24.org/spirals-of-creation/

### What the page actually specifies

- The Kathara diagram has 12 numbered centers and a central vertical 1–12 column; AzurA/CSA is placed at the midpoint of that column.
- The 1–12 column length and 6–7 width are stated to have ratio 2:1.
- The next grid is rotated around AzurA so that its center-9-to-AzurA ray overlays the previous grid's center-12-to-AzurA ray.
- At stage 2, the new center 6 overlays the stage-0 center 12.
- Each 45-degree stage scales lengths by sqrt(2); two stages rotate 90 degrees and scale by 2.
- The page distinguishes the Krystal/Kathara spiral from its Golden Mean and Fibonacci constructions. It gives the golden-mean rectangle construction and later discusses possible pivot placements of golden-mean/Fibonacci rectangles inside the Kathara grid.

### Exact mathematical content and comparison

The 45° clockwise similarity is
\[
T=\sqrt2R_{-\pi/4}=\begin{pmatrix}1&1\\-1&1\end{pmatrix},
\quad T^TT=2I,\quad T^2=2R_{-\pi/2}.
\]
This matches the prior square-similarity theorem T057. The overlay equations are
\[
Tp_9=p_{12},\qquad T^2p_6=p_{12},
\]
which imply \(p_9=Tp_6\). This is exactly the three-center relation already recorded in T057.

**Important correction to our previous interpretation:** the page says the ray \(9\to\mathrm{AzurA}\) overlays the prior ray \(12\to\mathrm{AzurA}\), not that a complete center coordinate is independently measured. The coordinate derivation is valid after choosing a Euclidean frame and using the stated 2:1 proportion and the stage similarity. It does not recover all 12 center coordinates or all line incidences by itself.

The page's claims about Veca-universe grids, energetic operation and cosmic structure are not mathematical consequences of the displayed similarity. The same page's golden-mean/Fibonacci constructions are distinct geometric constructions and should not be conflated with \(T\).

## 2. Les Ailes du Karma — “Spirale Fibonacci vs spirale Krystal”

URL: https://www.les-ailes-du-karma.com/blog/articles/spirale-fibonacci-vs-spirale-krystal

### Mathematical statements found

- It gives the Fibonacci recurrence \(F_{n+1}=F_n+F_{n-1}\), the sequence \(0,1,1,2,3,5,\ldots\), and the convergence of successive ratios to \(\varphi\).
- It describes a Kathara stage rule with eight 45° orientations per full turn and the 12-to-9 overlay.
- It gives the doubling sequence \(0,1,1,2,4,8,16,\ldots\) for the Krystal construction.
- It also says the next Krystal number is the sum of all previous values plus one.

### Internal consistency issue

The two numerical prescriptions are not identical as written. If \(K_0=0,K_1=1\) and
\[
K_{n+1}=\sum_{j=0}^{n}K_j+1,
\]
then the first terms are \(0,1,2,4,8,\ldots\), not \(0,1,1,2,4,\ldots\). The latter doubling sequence can be obtained with a different initialization/convention, but the page does not resolve the mismatch. We therefore retain only the independently specified geometric doubling \(2^m\) at 90° quadrant stages and \((\sqrt2)^n\) at 45° stages, not the ambiguous cumulative-sum recurrence.

The article's claims that Fibonacci/Golden Mean geometry is “anti-life”, vampiric, causes illness, or consumes biological energy are unsupported assertions, not mathematical statements or evidence about physiology. They are excluded from the RGS model.

## 3. FengshuiTell — Bagua / Eight Trigrams

URL: https://fengshuitell.com/es/blogs/guide/unveiling-the-mysteries-of-the-bagua-eight-trigrams-a-comprehensive-guide

### Useful content

The article correctly identifies eight basic trigram symbols and the combinatorial count \(8\times8=64\) ordered upper/lower combinations for hexagrams. It also describes yao as yin/yang lines and gives several traditional associations (Heaven, Earth, Fire, Thunder, Wind, Water, Mountain, Lake).

### Terminology and line-order corrections

The page repeatedly calls each three-line trigram a “hexagram”. This is terminologically incorrect: a trigram has three yao; a hexagram has six yao and can be formed by stacking an upper and a lower trigram.

Some line descriptions are incomplete or inconsistent: for example, the article describes Zhen, Kan, Gen and Dui using only upper/lower-line language without reliably specifying all three yao; some descriptions conflict with the actual three-line symbols. Do not use its prose as the source of the bit encoding.

For our explicit convention—bits ordered bottom-to-top, 0 = yin, 1 = yang—the eight exact triples are:
\[
\begin{array}{c|c}
000&\text{Kun }☷\\
001&\text{Gen }☶\\
010&\text{Kan }☵\\
011&\text{Xun }☴\\
100&\text{Zhen }☳\\
101&\text{Li }☲\\
110&\text{Dui }☱\\
111&\text{Qian }☰
\end{array}
\]
This is the convention already used in docs/TRIGRAM_Q3_WALSH_REFINEMENT_THEOREMS.md. The eight triples are \(V(Q_3)=\mathbb Z_2^3\); single-line flips give the 12 edges and the six coordinate-square faces. The 64 hexagrams are ordered pairs of trigrams, equivalently \(Q_3^{\rm lower}\times Q_3^{\rm upper}\cong Q_6\).

The Bagua article's cosmological, divinatory, health and Feng Shui applications are cultural/spiritual descriptions, not mathematical validation of Q3, Q6, Kathara, or RGS.

## 4. Kathara graph cross-check against the 2002 manual

The Emerald page's diagram gives the 12-center geometric arrangement and 2:1 proportions. The Dance for Life, Love, Freedom manual (2002), section “Structure of the Kathara Level-1 12-Tree Grid Core Scalar Template”, explicitly states 12 primary centers, 15 primary lines, three vertical line IDs 1–3, and central vertical line ID 1. The diagram-level endpoint transcription used by this repository is documented in docs/KATHARA_15_LINE_GRAPH_FACE_AUDIT.md.

The graph is a 12-vertex, 21-elementary-edge planar graph after expanding the three vertical polylines. Its chosen straight-line embedding has 10 bounded faces and one exterior face, satisfying V-E+F=2. These are faces of the transcribed planar graph. The source does not explicitly declare them as physical 2-cells, so a Kathara cellular complex still requires a separate modeling convention.

## 5. Comparison with the current RGS constructions

| Construction | Exact combinatorics | What is not identified |
|---|---|---|
| Kathara transcription | 12 vertices, 15 named primary lines, 21 elementary edges, 10 bounded planar faces | Not yet a chain complex with source-defined physical 2-cells |
| Trigram cube Q3 | 8 vertices, 12 edges, 6 square faces | No canonical map to 12 Kathara centers |
| Hexagram cube Q6 | 64 vertices, 192 edges, 240 square faces | No canonical map to the Kathara graph |
| Kathara stage similarity | T=sqrt(2)R(-45°), T²=2R(-90°) | Not a 3-channel transfer matrix |
| Tribonacci/supergolden transfer | Mpsi has polynomial lambda³-lambda²-1 | No incidence-derived map from Kathara graph to this matrix |

The graph and the trigram cube must remain separate objects until a map is explicitly defined and its incidence/orientation properties proved. In particular, 12 Kathara centers cannot be identified with 8 trigrams by relabeling alone.

## 6. Source references

- Emerald Awakening, “Spirals of Creation”: https://emerald24.org/spirals-of-creation/
- Anne Bocquet, “Spirale Fibonacci vs spirale Krystal”: https://www.les-ailes-du-karma.com/blog/articles/spirale-fibonacci-vs-spirale-krystal
- FengshuiTell, “Revelando los misterios del Bagua - Ocho trigramas”: https://fengshuitell.com/es/blogs/guide/unveiling-the-mysteries-of-the-bagua-eight-trigrams-a-comprehensive-guide
- Dance for Life, Love, Freedom (2002), section “Structure of the Kathara Level-1 12-Tree Grid Core Scalar Template”: https://fliphtml5.com/wcfyn/imae/2002.08_Dance_For_Life_Love_Freedom_Manual_original_25MB/84/
