from dataclasses import dataclass
from typing import List, Optional
import numpy as np
from conteo import ContadorOperaciones


@dataclass
class ResultadoGauss:
    estrategia: str
    x: np.ndarray
    U: np.ndarray
    L: np.ndarray
    P: np.ndarray
    c_mod: np.ndarray
    contador: ContadorOperaciones
    max_multiplicador: float
    factor_crecimiento: float
    max_inicial: float
    max_proceso: float
    traza_pasos: List[str]
    perm: np.ndarray


def _formatear_matriz_aumentada(A: np.ndarray, b: np.ndarray, titulo: str) -> str:
    n = A.shape[0]
    lineas = [titulo]
    header = "Fila | " + " ".join(f"{'c'+str(j+1):>14s}" for j in range(n)) + f" | {'b':>15s}"
    lineas.append(header)
    lineas.append("-" * len(header))
    for i in range(n):
        fila_str = " ".join(f"{A[i, j]:14.6e}" for j in range(n))
        lineas.append(f"F{i+1:02d}  | {fila_str} | {b[i]:15.7e}")
    lineas.append("")
    return "\n".join(lineas)


def eliminar_y_resolver(
    A_in: np.ndarray,
    b_in: np.ndarray,
    estrategia: str = "escalado",
) -> ResultadoGauss:
    """
    Resuelve Ax = b por eliminación gaussiana con mínimo costo operacional.
    estrategia: 'sin_pivoteo', 'parcial', 'escalado'
    """
    A = A_in.astype(np.float64).copy()
    b = b_in.astype(np.float64).copy()
    n = A.shape[0]

    L = np.eye(n, dtype=np.float64)
    perm = np.arange(n)
    contador = ContadorOperaciones(estrategia=estrategia)

    max_inicial = float(np.max(np.abs(A)))
    max_proceso = max_inicial
    max_multiplicador = 0.0
    traza_pasos: List[str] = []

    traza_pasos.append(_formatear_matriz_aumentada(A, b, f"=== ESTADO INICIAL ({estrategia.upper()}) ==="))

    # Factores de escala por fila (solo una vez al inicio para pivoteo escalado)
    escalas: Optional[np.ndarray] = None
    if estrategia == "escalado":
        escalas = np.zeros(n, dtype=np.float64)
        for i in range(n):
            max_fila = abs(A[i, 0])
            for j in range(1, n):
                contador.escala_comparaciones += 1
                val = abs(A[i, j])
                if val > max_fila:
                    max_fila = val
            escalas[i] = max_fila

    # Eliminación hacia adelante (11 pasos)
    for k in range(n - 1):
        paso_cnt = contador.iniciar_paso(k + 1)

        # Selección de pivote según estrategia
        p = k
        if estrategia == "parcial":
            mejor_val = abs(A[k, k])
            for i in range(k + 1, n):
                paso_cnt.comparaciones += 1
                val = abs(A[i, k])
                if val > mejor_val:
                    mejor_val = val
                    p = i
        elif estrategia == "escalado":
            assert escalas is not None
            contador.escala_divisiones += 1
            mejor_cociente = abs(A[k, k]) / escalas[k]
            for i in range(k + 1, n):
                contador.escala_divisiones += 1
                paso_cnt.comparaciones += 1
                cociente = abs(A[i, k]) / escalas[i]
                if cociente > mejor_cociente:
                    mejor_cociente = cociente
                    p = i

        # Intercambio de filas si corresponde
        if p != k:
            A[[k, p], :] = A[[p, k], :]
            b[[k, p]] = b[[p, k]]
            perm[[k, p]] = perm[[p, k]]
            if k > 0:
                L[[k, p], :k] = L[[p, k], :k]
            if escalas is not None:
                escalas[[k, p]] = escalas[[p, k]]
            paso_cnt.intercambio = 1

        pivote = A[k, k]
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote nulo en el paso k={k+1}")

        # Barrido de filas inferiores
        for i in range(k + 1, n):
            # Un único multiplicador por fila (1 división)
            m_ik = A[i, k] / pivote
            paso_cnt.mult_div += 1

            L[i, k] = m_ik
            if abs(m_ik) > max_multiplicador:
                max_multiplicador = float(abs(m_ik))

            # La subdiagonal queda exactamente en cero sin gastar FLOPs
            A[i, k] = 0.0

            # Barrido de columnas desde k+1 en adelante
            for j in range(k + 1, n):
                A[i, j] = A[i, j] - m_ik * A[k, j]
                paso_cnt.mult_div += 1
                paso_cnt.sum_rest += 1
                if abs(A[i, j]) > max_proceso:
                    max_proceso = float(abs(A[i, j]))

            # Actualización del lado derecho b
            b[i] = b[i] - m_ik * b[k]
            paso_cnt.mult_div += 1
            paso_cnt.sum_rest += 1

        info_inter = f"Intercambio F{k+1} <-> F{p+1}" if p != k else "Sin intercambio"
        titulo_paso = (
            f"--- PASO {k+1:02d} ({estrategia}) | Fila original: F{perm[k]+1} | "
            f"Pivote = {pivote:.6e} | {info_inter} ---"
        )
        traza_pasos.append(_formatear_matriz_aumentada(A, b, titulo_paso))

    # Sustitución regresiva
    x = np.zeros(n, dtype=np.float64)
    x[n - 1] = b[n - 1] / A[n - 1, n - 1]
    contador.sust_mult_div += 1

    for i in range(n - 2, -1, -1):
        suma = 0.0
        for j in range(i + 1, n):
            prod = A[i, j] * x[j]
            contador.sust_mult_div += 1
            if j == i + 1:
                suma = prod
            else:
                suma = suma + prod
                contador.sust_sum_rest += 1
        x[i] = (b[i] - suma) / A[i, i]
        contador.sust_sum_rest += 1
        contador.sust_mult_div += 1

    P = np.eye(n, dtype=np.float64)[perm, :]
    factor_crecimiento = max_proceso / max_inicial

    return ResultadoGauss(
        estrategia=estrategia,
        x=x,
        U=A,
        L=L,
        P=P,
        c_mod=b,
        contador=contador,
        max_multiplicador=max_multiplicador,
        factor_crecimiento=factor_crecimiento,
        max_inicial=max_inicial,
        max_proceso=max_proceso,
        traza_pasos=traza_pasos,
        perm=perm,
    )


def resolver_con_lu(L: np.ndarray, U: np.ndarray, perm: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    """
    Resuelve Ax = b reutilizando los factores L, U y el vector de permutación perm
    sin volver a factorizar la matriz A.
    """
    n = L.shape[0]
    b_p = b_in[perm].astype(np.float64)

    # Sustitución progresiva: L y = P b (L tiene diagonal unitaria)
    y = np.zeros(n, dtype=np.float64)
    for i in range(n):
        suma = 0.0
        for j in range(i):
            suma += L[i, j] * y[j]
        y[i] = b_p[i] - suma

    # Sustitución regresiva: U x = y
    x = np.zeros(n, dtype=np.float64)
    for i in range(n - 1, -1, -1):
        suma = 0.0
        for j in range(i + 1, n):
            suma += U[i, j] * x[j]
        x[i] = (y[i] - suma) / U[i, i]

    return x