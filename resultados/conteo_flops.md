# Conteo de Operaciones Elementales (FLOPs)

## Convención Declarada
- **1 FLOP** = 1 operación aritmética de punto flotante (`+`, `-`, `*`, `/`).
- En cada paso $k \in \{1, \dots, 11\}$ y para cada fila $i > k$:
  - Se calcula **una sola vez** el multiplicador $m_{ik} = a_{ik} / a_{kk}$ (1 división).
  - Se fija $a_{ik} = 0$ directamente sin operar la subdiagonal ya anulada.
  - Se actualizan las columnas $j \in \{k+1, \dots, n\}$ con $a_{ij} \leftarrow a_{ij} - m_{ik} a_{kj}$ (1 multiplicación y 1 resta por columna).
  - Se actualiza el vector derecho $b_i \leftarrow b_i - m_{ik} b_k$ (1 multiplicación y 1 resta).
- El costo de búsqueda de pivote, cálculo de factores de escala e intercambios se contabiliza aparte como **sobrecosto de estabilidad**.

## Desglose Paso a Paso de la Eliminación ($n = 12$)

| Paso $k$ | Filas eliminadas | Mult / Div | Sum / Rest | Total FLOPs Paso | Comp. Pivoteo (Parcial) | Intercambio (Escalado) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 11 | 143 | 132 | 275 | 11 | 1 |
| 2 | 10 | 120 | 110 | 230 | 10 | 0 |
| 3 | 9 | 99 | 90 | 189 | 9 | 1 |
| 4 | 8 | 80 | 72 | 152 | 8 | 1 |
| 5 | 7 | 63 | 56 | 119 | 7 | 1 |
| 6 | 6 | 48 | 42 | 90 | 6 | 1 |
| 7 | 5 | 35 | 30 | 65 | 5 | 1 |
| 8 | 4 | 24 | 20 | 44 | 4 | 1 |
| 9 | 3 | 15 | 12 | 27 | 3 | 1 |
| 10 | 2 | 8 | 6 | 14 | 2 | 0 |
| 11 | 1 | 3 | 2 | 5 | 1 | 0 |
| **Subtotal Eliminación** | **66** | **638** | **572** | **1210** | **66** | **8** |
| **Sustitución Regresiva** | - | **78** | **66** | **144** | 0 | 0 |
| **TOTAL GENERAL** | - | **716** | **638** | **1354** | **66** | **8** |

## Comparación de Estrategias y Cotas de Referencia

| Estrategia | Mult / Div | Sum / Rest | Total FLOPs Aritméticos | Cota de Referencia | Diferencia vs. 1486 | Sobrecosto de Estabilidad |
|:---|:---:|:---:|:---:|:---:|:---:|:---|
| Sin pivoteo | 716 | 638 | **1354** | 1486 | -132 (-8.88%) | 0 comparaciones, 0 intercambios |
| Pivoteo parcial | 716 | 638 | **1354** | 1486 | -132 (-8.88%) | 66 comparaciones, 9 intercambios |
| Pivoteo parcial escalado | 716 | 638 | **1354** | 1486 | -132 (-8.88%) | 132 comp. escala + 77 div. cocientes + 66 comp. pivote, 8 intercambios |

El total de 1354 FLOPs es el mínimo estricto bajo la convención declarada (eliminación + actualización de b + sustitución regresiva). La cota de 1486 es la referencia del enunciado; la diferencia se calcula respecto de ella.
