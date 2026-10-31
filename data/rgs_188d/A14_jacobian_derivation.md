# A14 — Jacobiana RGS–PSU / QLG

## Definizione

Per i 64 siti $X_i=(x_i,y_i,z_i)$, $X\in\mathbb R^{64\times3}$, i vincoli sono

$$F(X)=\begin{pmatrix}\sum_i x_i\\\sum_i y_i\\\sum_i z_i\\\frac12(\sum_i\|X_i\|^2-R_0^2)\end{pmatrix}=0.$$

La Jacobiana rispetto a $X_{flat}\in\mathbb R^{192}$ è

$$J=\frac{\partial F}{\partial X_{flat}}=\begin{pmatrix}1&0&0&1&0&0&\cdots\\0&1&0&0&1&0&\cdots\\0&0&1&0&0&1&\cdots\\x_1&y_1&z_1&x_2&y_2&z_2&\cdots\end{pmatrix}\in\mathbb R^{4\times192}.$$

Il quarto vincolo usa $R_0^2=\sum_i\|X_i\|^2$ dell'embedding archiviato: è la varietà a scala fissata passante per $X$. La sua derivata è $X_{flat}^T$.

## Esito numerico

- rank(J)=4; dim ker(J)=188.
- Il test indipendente è $\|JQ\|_F$ e il confronto fra proiettori $QQ^T$ e $Q_{SVD}Q_{SVD}^T$.
- La base SVD delle normali è una scelta ortonormale; le etichette traslazione x/y/z e radiale appartengono alle righe di $J$, non a una scelta arbitraria di quattro vettori normali.

Il report JSON contiene i residui misurati.
