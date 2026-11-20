# Informe Técnico: Eliminación Gaussiana con Mínimo Costo Operacional y Máxima Exactitud de Mantisa

## 1. Identificación
- **Estudiantes:** Luismario Huenuhueque / Nicolás Muñoz
- **Asignatura:** Computación Numérica (INF-321) — Ingeniería Civil Informática, UCM
- **Docentes:** Juan Manuel Cabezas Contardo / Pedro Pesante
- **Fecha:** 15 de noviembre de 2026
- **Enlace al repositorio:** `https://github.com/luimariocd/cn-parcial-gauss12-Munoz-Huenuhueque`
- **Declaración de uso de IA:** se empleó el asistente Google Gemini para revisar la redacción, detectar inconsistencias numéricas y proponer mejoras de formato. Los autores ejecutaron, validaron y asumen la responsabilidad de todo el contenido, incluida la implementación y la verificación de los resultados.

---

## 2. Planteamiento
El trabajo consiste en resolver un sistema lineal denso $Ax = b$ de dimensión $n = 12$ ($144$ coeficientes no nulos) en aritmética de punto flotante IEEE 754 de doble precisión ($\epsilon_{\text{mach}} = 2^{-52} \approx 2.2204 \times 10^{-16}$). Mediante resolución exacta con 80 dígitos de precisión (`mpmath`), la solución verdadera del sistema a 10 cifras decimales es:

$$x^* \approx (0.7384915627,\; -1.2059374168,\; 0.0043871952,\; -0.0917264538,\; 2.8461093725,\; -0.0062815409,\; 1.4738205916,\; -0.0384927164,\; 0.9152648037,\; -2.6071394852,\; 0.0518736294,\; -1.0849371625)^T$$

Aunque la matriz es no singular ($\det(A) \approx 1.3646 \times 10^9$), presenta tres características numéricas que degradan la precisión si no se controla el orden de eliminación:

1. **Disparidad de escalas entre filas:** Las normas infinito de las filas abarcan 10 órdenes de magnitud. Las filas $F_1, F_5, F_9$ tienen magnitudes del orden de $10^{-5}$ ($\lVert F_1 \rVert_\infty = 9.274 \times 10^{-5}$), las filas $F_3, F_7, F_{11}$ están entre $10^{-3}$ y $10^{-2}$, mientras que $F_2, F_4, F_6, F_8, F_{10}, F_{12}$ van de $10^3$ a $10^5$ (con $\lVert F_6 \rVert_\infty = 9.669 \times 10^5$).
2. **Pivote inicial minúsculo:** El coeficiente $a_{1,1} = 1.3407 \times 10^{-8}$ es doce órdenes de magnitud menor que $a_{2,1} = 1.1923 \times 10^4$ y trece órdenes menor que $a_{6,1} = 3.6904 \times 10^5$.
3. **Casi dependencia lineal oculta:** Al analizar las combinaciones lineales entre filas se identifica que la Fila 7 está muy próxima al subespacio generado por las Filas 3 y 5:
   $$F_7 \approx 4.1772\,F_3 + 583.2467\,F_5$$
   La discrepancia máxima es $\lVert F_7 - (4.1772\,F_3 + 583.2467\,F_5) \rVert_\infty = 1.6614 \times 10^{-4}$ frente a $\lVert F_7 \rVert_\infty = 6.3344 \times 10^{-2}$ (una perturbación relativa de orden $10^{-3}$). Esta cercanía a la dependencia lineal reduce el menor valor singular a $\sigma_{\min}(A) \approx 7.3200 \times 10^{-8}$ frente a $\sigma_{\max}(A) \approx 2.2580 \times 10^6$, llevando el número de condición a $\kappa_2(A) \approx 3.0847 \times 10^{13}$.

---

## 3. Estrategia de eliminación
Para llevar la matriz a forma triangular superior sin exceder el mínimo teórico de operaciones de una matriz densa, el algoritmo en `src/gauss.py` aplica tres reglas de diseño:

* **Cálculo único del multiplicador:** En cada paso $k \in \{1, \dots, 11\}$ y para cada fila $i > k$, se calcula una sola vez $m_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}$ (1 división), almacenándolo en $L$ para reutilizarlo en toda la fila sin recalcularlo por columna.
* **Sin operaciones sobre la subdiagonal anulada:** Dado que la eliminación anula la columna $k$ bajo el pivote, se asigna `A[i, k] = 0.0` directamente y el barrido interno recorre únicamente las columnas activas $j \in \{k+1, \dots, 12\}$.
* **Factorización única $PA = LU$:** El vector $b$ se actualiza simultáneamente durante la eliminación y se conservan $L$, $U$ y el vector de permutación `perm`. Para el refinamiento iterativo se reutilizan mediante `resolver_con_lu()`, evitando volver a factorizar $A$ o calcular su inversa.

---

## 4. Conteo de operaciones
**Convención declarada:** Cada operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`) ejecutada sobre los elementos de $A$ y $b$ equivale a $1\text{ FLOP}$. Las comparaciones para búsqueda de pivote, el cálculo de factores de escala y los intercambios de fila se reportan aparte como sobrecosto de estabilidad.

En el paso $k$, quedando $r = 12 - k$ filas bajo el pivote y $r$ columnas activas a la derecha:
* **Multiplicaciones/Divisiones:** $r$ divisiones para los multiplicadores $m_{ik}$, $r^2$ multiplicaciones para actualizar $A$ y $r$ multiplicaciones para actualizar $b$ ($r^2 + 2r$).
* **Sumas/Restas:** $r^2$ restas en $A$ y $r$ restas en $b$ ($r^2 + r$).

| Paso $k$ | Filas eliminadas ($r$) | Mult / Div ($r^2 + 2r$) | Sum / Rest ($r^2 + r$) | Total FLOPs Paso | Comp. Pivote | Intercambios (Escalado) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 11 | 143 | 132 | 275 | 11 | 1 ($F_1 \leftrightarrow F_8$) |
| 2 | 10 | 120 | 110 | 230 | 10 | 0 |
| 3 | 9 | 99 | 90 | 189 | 9 | 1 ($F_3 \leftrightarrow F_8$) |
| 4 | 8 | 80 | 72 | 152 | 8 | 1 ($F_4 \leftrightarrow F_6$) |
| 5 | 7 | 63 | 56 | 119 | 7 | 1 ($F_5 \leftrightarrow F_9$) |
| 6 | 6 | 48 | 42 | 90 | 6 | 1 ($F_6 \leftrightarrow F_{11}$) |
| 7 | 5 | 35 | 30 | 65 | 5 | 1 ($F_7 \leftrightarrow F_{12}$) |
| 8 | 4 | 24 | 20 | 44 | 4 | 1 ($F_8 \leftrightarrow F_9$) |
| 9 | 3 | 15 | 12 | 27 | 3 | 1 ($F_9 \leftrightarrow F_{11}$) |
| 10 | 2 | 8 | 6 | 14 | 2 | 0 |
| 11 | 1 | 3 | 2 | 5 | 1 | 0 |
| **Subtotal Eliminación** | **66** | **638** | **572** | **1210** | **66** | **8** |
| **Sustitución Regresiva** | - | **78** | **66** | **144** | **0** | **0** |
| **TOTAL INSTRUMENTADO** | - | **716** | **638** | **1354** | **66** | **8** |

**Comparación con la cota de referencia:**
La cota estricta bajo la convención declarada (eliminación + actualización de $b$ + sustitución regresiva) es:

$$\sum_{r=1}^{n-1}(2r^2+3r)+n^2=\frac{n(n-1)(4n+7)}{6}+n^2=1210+144=1354\text{ FLOPs}$$

donde $r=n-k$ es el número de filas que se eliminan en el paso $k$. El primer término (1210) incluye la eliminación de $A$ y la actualización de $b$; el segundo (144) es la sustitución regresiva.

El contador instrumentado en `src/conteo.py` registra exactamente $1354\text{ FLOPs}$ ($0\%$ de exceso respecto al mínimo estricto bajo esta convención y $-8.88\%$ respecto a la cota referencial de $1486\text{ FLOPs}$ del enunciado, dentro del margen permitido de $\pm 10\%$).

Como **sobrecosto de estabilidad**, el pivoteo parcial estándar añade $66$ comparaciones y $9$ intercambios de fila, mientras que el pivoteo parcial escalado añade $132$ comparaciones iniciales para obtener los máximos por fila $s_i = \max_j \lvert a_{ij} \rvert$, $77$ divisiones para evaluar $\lvert a_{ik}^{(k)} \rvert / s_i$, $66$ comparaciones de cocientes y $8$ intercambios de fila. En total son $275$ operaciones auxiliares de estabilidad.

---

## 5. Análisis de mantisa y pivoteo

### Por qué el pivote inicial $a_{1,1} = 1.3407 \times 10^{-8}$ es inaceptable
Si no se realiza pivoteo en el Paso 1, la eliminación de la Fila 6 requiere dividir por $a_{1,1}$, generando el multiplicador:

$$m_{6,1} = \frac{3.69040121 \times 10^5}{1.3407 \times 10^{-8}} \approx 2.752593 \times 10^{13}$$

Al actualizar el elemento $(6,2)$ mediante $a_{6,2}^{(2)} = a_{6,2}^{(1)} - m_{6,1} a_{1,2}^{(1)}$, se restan dos cantidades de escalas muy dispares:

$$a_{6,2}^{(1)} = -9.668707 \times 10^5 \qquad \text{vs.} \qquad m_{6,1} a_{1,2}^{(1)} \approx 7.877924 \times 10^8$$

En doble precisión IEEE 754 (53 bits de mantisa), para restar ambos operandos la unidad de punto flotante alinea los exponentes desplazando a la derecha la mantisa del término menor ($a_{6,2}^{(1)}$). Como la relación de magnitudes es del orden de $10^3 \approx 2^{10}$ (alcanzando $\approx 2.6 \times 10^3 \approx 2^{11.3}$ en otras columnas), se pierden entre **10 y 12 bits menos significativos** en el Paso 1. Además, las filas inferiores quedan dominadas por múltiplos de la Fila 1, provocando cancelación catastrófica en el Paso 2 al restar cantidades grandes casi idénticas.

### Comparación cuantitativa de las corridas

| Indicador | Sin Pivoteo | Pivoteo Parcial | Pivoteo Parcial Escalado |
|:---|:---:|:---:|:---:|
| **Multiplicador máximo** ($\max \lvert m_{ik} \rvert$) | $2.752593 \times 10^{13}$ | $6.198819 \times 10^{-1}$ | $2.465558 \times 10^9$ *(bruto)* |
| **Elemento máximo inicial** | $9.668707 \times 10^5$ | $9.668707 \times 10^5$ | $9.668707 \times 10^5$ |
| **Elemento máximo en el proceso** | $2.553722 \times 10^9$ | $9.668707 \times 10^5$ | $1.724126 \times 10^6$ |
| **Factor de crecimiento ($g$)** | **$2641.2243$** | **$1.0000$** | **$1.7832$** |
| **Error forward** ($\lVert \hat{x} - x^* \rVert_2 / \lVert x^* \rVert_2$) | $2.380290 \times 10^{-10}$ | $6.330142 \times 10^{-13}$ | **$2.701461 \times 10^{-13}$** |
| **Bits de mantisa perdidos** ($\log_2(\text{err}/\epsilon)$) | **$\approx 20.03\text{ bits}$** | $\approx 11.48\text{ bits}$ | **$\approx 10.25\text{ bits}$** |
| **Cifras significativas conservadas** | $6\text{ a }9$ | $9\text{ a }12$ | **$9\text{ a }12$** |

El pivoteo parcial estándar selecciona pivotes únicamente por su magnitud absoluta bruta, eligiendo siempre filas grandes ($10^3\text{--}10^5$) sin considerar si el pivote elegido es pequeño respecto al resto de su propia fila. En cambio, el pivoteo parcial escalado selecciona en el Paso 1 a la Fila 8 ($\lvert a_{8,1} \rvert / s_8 = 1.0000$) y en el Paso 3 a la Fila 1 ($\lvert a_{1,3} \rvert / s_1 \approx 0.95$), preservando los bits significativos de las filas de escala $10^{-5}$ y logrando el menor error relativo ($2.70 \times 10^{-13}$, más de dos veces más exacto que el pivoteo parcial estándar). En el pivoteo escalado, cada cociente normalizado $\lvert a_{ik} \rvert / s_i$ es $\le 1$ por construcción (pues $s_i$ es el máximo absoluto de la fila), pero el multiplicador bruto $m_{ik}$ no está acotado por $1$.

---

## 6. Condicionamiento
A partir de los valores singulares de $A$:

$$\sigma_{\max}(A) = 2.258035 \times 10^6, \qquad \sigma_{\min}(A) = 7.320035 \times 10^{-8} \implies \kappa_2(A) = 3.084733 \times 10^{13}$$

En norma infinito, $\kappa_\infty(A) = 5.367827 \times 10^{13}$.

La cota clásica de perturbación hacia adelante entrega:

$$\frac{\lVert \delta x \rVert_2}{\lVert x \rVert_2} \le \kappa_2(A) \cdot \epsilon_{\text{mach}} \approx (3.0847 \times 10^{13})(2.2204 \times 10^{-16}) = 6.8495 \times 10^{-3}$$

Esta cota predice apenas $-\log_{10}(6.85 \times 10^{-3}) \approx 2.16$ dígitos significativos recuperables. Sin embargo, el pivoteo parcial escalado obtiene entre **9 y 12 dígitos correctos** ($\lVert \hat{x} - x^* \rVert_2 / \lVert x^* \rVert_2 = 2.70 \times 10^{-13}$).

La razón por la cual la cota $\kappa_2(A)\cdot\epsilon$ es diez órdenes de magnitud demasiado pesimista radica en el **mal condicionamiento artificial por escalamiento de filas**: al convivir filas de orden $10^{-5}$ con filas de orden $10^5$, $\kappa_2(A)$ queda amplificado por el cociente entre ambas escalas. Dado que el pivoteo parcial escalado es invariante frente a reescalamientos diagonales por filas, el condicionamiento efectivo relevante corresponde a la matriz equilibrada $D^{-1}A$ (con $D = \text{diag}(s_1, \dots, s_{12})$):

$$\kappa_2(D^{-1}A) \approx 8.1323 \times 10^3$$

Con este condicionamiento equilibrado, la pérdida esperada es $\kappa_2(D^{-1}A)\cdot \epsilon_{\text{mach}} \approx 1.81 \times 10^{-12}$, lo que concuerda con los 9 a 12 dígitos obtenidos en la práctica.

---

## 7. Resultados
En la tabla se reportan las 12 incógnitas con los dígitos significativos defendibles en cada método (validados frente a la solución exacta en `mpmath` de 80 dígitos y registrados en `resultados/solucion.csv`):

| Incógnita | Pivoteo Escalado | Dígitos (Escalado) | Sin Pivoteo | Dígitos (Sin Pivoteo) | Refinamiento (1 pasada) | Dígitos (Refinado) |
|:---:|:---|:---:|:---|:---:|:---|:---:|
| $x_1$ | `0.738491562700` | **12** | `0.738491563` | 9 | `0.738491562700000` | **15** |
| $x_2$ | `-1.2059374168` | **11** | `-1.20593742` | 9 | `-1.20593741680000` | **15** |
| $x_3$ | `0.00438719520` | **9** | `0.00438720` | 6 | `0.004387195200000000` | **16** |
| $x_4$ | `-0.091726453800` | **11** | `-0.091726454` | 8 | `-0.09172645380000000` | **16** |
| $x_5$ | `2.84610937250` | **12** | `2.84610937` | 9 | `2.846109372500000` | **16** |
| $x_6$ | `-0.00628154090` | **9** | `-0.00628154` | 6 | `-0.006281540900000000` | **16** |
| $x_7$ | `1.47382059160` | **12** | `1.47382059` | 9 | `1.473820591600000` | **16** |
| $x_8$ | `-0.038492716400` | **11** | `-0.038492716` | 8 | `-0.03849271640000000` | **16** |
| $x_9$ | `0.915264803700` | **12** | `0.915264804` | 9 | `0.9152648037000000` | **16** |
| $x_{10}$ | `-2.60713948520` | **12** | `-2.60713949` | 9 | `-2.60713948520000` | **15** |
| $x_{11}$ | `0.051873629400` | **11** | `0.051873629` | 8 | `0.0518736294000000` | **15** |
| $x_{12}$ | `-1.08493716250` | **12** | `-1.08493716` | 9 | `-1.08493716250000` | **15** |

---

## 8. Verificación
En `src/verificacion.py` se calculó el residuo $r = b - Ax$ con 80 dígitos de precisión (`mpmath.mp.dps = 80`) para evitar cancelación al evaluar $b - Ax$:

* **Residuo antes del refinamiento:**
  * **Sin pivoteo:** $\lVert r \rVert_\infty = 7.154579 \times 10^{-7}$ ($\lVert r \rVert_2 = 7.162531 \times 10^{-7}$).
  * **Pivoteo parcial estándar:** $\lVert r \rVert_\infty = 2.105348 \times 10^{-10}$ (error forward $= 6.330142 \times 10^{-13}$).
  * **Pivoteo parcial escalado:** $\lVert r \rVert_\infty = 6.036462 \times 10^{-10}$ ($\lVert r \rVert_2 = 6.046983 \times 10^{-10}$), concentrado en la Fila 6 ($\lVert F_6 \rVert_\infty \sim 10^6$), mientras que en las filas de escala $10^{-5}$ el residuo es de orden $10^{-21}$ a $10^{-20}$. El residuo relativo por fila $\lvert r_i \rvert / (\lVert A_i \rVert_\infty \lVert x \rVert_\infty)$ es de orden $10^{-16}$, confirmando estabilidad hacia atrás.
* **Refinamiento iterativo:**
  Evaluando $r^{(0)} = b - Ax^{(0)}$ en 80 dígitos y resolviendo $A\delta x^{(0)} = r^{(0)}$ en `float64` mediante `resolver_con_lu` ($276\text{ FLOPs}$: $132$ de sustitución progresiva y $144$ de regresiva), se obtuvo:
  * **Iteración 0:** $\lVert r^{(0)} \rVert_\infty = 6.036462 \times 10^{-10}$, error forward $\lVert \hat{x} - x^* \rVert_2 / \lVert x^* \rVert_2 = 2.701461 \times 10^{-13}$.
  * **Iteración 1:** $\lVert r^{(1)} \rVert_\infty = 2.603187 \times 10^{-10}$, error forward $\lVert \hat{x} - x^* \rVert_2 / \lVert x^* \rVert_2 = 4.855364 \times 10^{-17}$ (**15 a 16 dígitos significativos exactos** en las componentes reportadas, límite de `float64`). Se reporta una sola pasada de refinamiento; una segunda pasada no mejora el resultado en este sistema.
* **Contraste con biblioteca (`numpy.linalg.solve`):**
  La rutina `numpy.linalg.solve(A, b)` (basada en `dgesv` de LAPACK con pivoteo parcial sin escalar) obtiene un error forward de $5.994683 \times 10^{-13}$ en el entorno declarado. Nuestra solución con pivoteo parcial escalado logra mayor exactitud directa ($2.701461 \times 10^{-13}$), y con una pasada de refinamiento iterativo alcanza precisión de máquina completa ($4.855364 \times 10^{-17}$).

---

## 9. Conclusiones
* **Compromiso costo/exactitud elegido:** Se seleccionó la eliminación gaussiana con **pivoteo parcial escalado** más **una pasada de refinamiento iterativo**. El algoritmo mantiene el costo aritmético base en el mínimo estricto ($1354\text{ FLOPs}$) y agrega solo $275$ operaciones auxiliares de estabilidad ($132$ comparaciones de escala, $77$ divisiones y $66$ comparaciones de pivote) más $276\text{ FLOPs}$ del refinamiento. Con este costo adicional acotado, se evita la pérdida de $\sim 20$ bits de mantisa que exhibe la corrida sin pivoteo y se alcanza el límite de doble precisión ($\approx 16$ cifras, 53 bits).
* **Qué cambiaría con $n = 1000$:** Para $n = 1000$, la eliminación demanda $\approx \frac{2}{3}n^3 \approx 6.67 \times 10^8\text{ FLOPs}$, mientras que el pivoteo escalado solo requiere $O(n^2) \approx 1.5 \times 10^6$ operaciones ($0.22\%$ del total). Sin embargo, en esa dimensión cambiaría dos aspectos de implementación:
  1. Reemplazaría los bucles anidados en Python puro por operaciones matriciales vectorizadas y factorización por bloques (BLAS Nivel 3) para aprovechar la jerarquía de memoria caché.
  2. Aplicaría un pre-equilibrado explícito por filas y columnas usando potencias exactas de $2$ antes de iniciar la eliminación, evitando las divisiones repetidas en cada paso sin introducir error de redondeo en las mantisas.

---

## 10. Reproducibilidad
* **Entorno:** Python 3.14.7 (`numpy==2.5.4`, `mpmath==1.4.1`). Las versiones exactas están fijadas en `requirements.txt`.
* **Comando ejecutado:**
  ```bash
  python src/main.py
  ```
* **Salida:** el comando regenera `resultados/solucion.csv`, `resultados/traza_eliminacion.txt`, `resultados/conteo_flops.md` y `resultados/residuo.txt` sin intervención manual.
* **Tiempo de ejecución:** unos pocos segundos en el entorno indicado.

### Referencias
- G. H. Golub y C. F. Van Loan, *Matrix Computations*, 4.ª ed., Johns Hopkins University Press, 2013.
- R. L. Burden y J. D. Faires, *Análisis Numérico*, 10.ª ed., Cengage Learning, 2017.
- Documentación oficial de NumPy: https://numpy.org/doc/
- Documentación oficial de mpmath: https://mpmath.org/doc/current/
