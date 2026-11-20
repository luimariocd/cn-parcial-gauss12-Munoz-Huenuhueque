from dataclasses import dataclass
from typing import List, Tuple
import numpy as np
import mpmath as mp
from gauss import resolver_con_lu

mp.mp.dps = 80


@dataclass
class AnalisisExactitud:
    x_exacta: List[mp.mpf]
    errores_relativos_comp: np.ndarray
    digitos_correctos_comp: List[int]
    error_relativo_norm2: float
    residuo_ext: List[mp.mpf]
    norma_inf_residuo: mp.mpf
    norma_2_residuo: mp.mpf


def leer_datos_texto(ruta_a: str, ruta_b: str) -> Tuple[List[List[str]], List[str]]:
    with open(ruta_a, "r", encoding="utf-8") as fa:
        filas_a = [line.strip().split(",") for line in fa if line.strip()]
    with open(ruta_b, "r", encoding="utf-8") as fb:
        filas_b = [line.strip() for line in fb if line.strip()]
    return filas_a, filas_b


def calcular_solucion_exacta_mp(filas_a_str: List[List[str]], filas_b_str: List[str]) -> List[mp.mpf]:
    n = len(filas_b_str)
    A_mp = mp.matrix([[mp.mpf(val) for val in fila] for fila in filas_a_str])
    b_mp = mp.matrix([mp.mpf(val) for val in filas_b_str])
    x_mp = mp.lu_solve(A_mp, b_mp)
    return [x_mp[i] for i in range(n)]


def calcular_residuo_extendido(
    filas_a_str: List[List[str]],
    filas_b_str: List[str],
    x_float: np.ndarray,
) -> Tuple[List[mp.mpf], mp.mpf, mp.mpf]:
    n = len(filas_b_str)
    A_mp = [[mp.mpf(val) for val in fila] for fila in filas_a_str]
    b_mp = [mp.mpf(val) for val in filas_b_str]
    x_mp = [mp.mpf(float(val)) for val in x_float]

    r_mp: List[mp.mpf] = []
    for i in range(n):
        ax_i = mp.fsum(A_mp[i][j] * x_mp[j] for j in range(n))
        r_mp.append(b_mp[i] - ax_i)

    norma_inf = max(abs(val) for val in r_mp)
    norma_2 = mp.sqrt(mp.fsum(val * val for val in r_mp))
    return r_mp, norma_inf, norma_2


def evaluar_exactitud(
    filas_a_str: List[List[str]],
    filas_b_str: List[str],
    x_calc: np.ndarray,
    x_exacta: List[mp.mpf],
) -> AnalisisExactitud:
    n = len(x_exacta)
    errores_comp = np.zeros(n, dtype=np.float64)
    digitos_comp: List[int] = []

    num_norm2 = mp.mpf("0")
    den_norm2 = mp.mpf("0")

    for i in range(n):
        xc = mp.mpf(float(x_calc[i]))
        xe = x_exacta[i]
        diff = abs(xc - xe)
        rel = diff / abs(xe)
        errores_comp[i] = float(rel)

        num_norm2 += diff * diff
        den_norm2 += xe * xe

        if rel == 0:
            digitos_comp.append(16)
        else:
            # Criterio conservador: se descuenta medio ulp decimal.
            # d cifras correctas <==> error relativo < 0.5 * 10^(-d).
            d = int(mp.floor(-mp.log10(2 * rel)))
            digitos_comp.append(max(0, min(16, d)))

    err_norm2 = float(mp.sqrt(num_norm2) / mp.sqrt(den_norm2))
    r_mp, n_inf, n_2 = calcular_residuo_extendido(filas_a_str, filas_b_str, x_calc)

    return AnalisisExactitud(
        x_exacta=x_exacta,
        errores_relativos_comp=errores_comp,
        digitos_correctos_comp=digitos_comp,
        error_relativo_norm2=err_norm2,
        residuo_ext=r_mp,
        norma_inf_residuo=n_inf,
        norma_2_residuo=n_2,
    )


def refinamiento_iterativo(
    filas_a_str: List[List[str]],
    filas_b_str: List[str],
    L: np.ndarray,
    U: np.ndarray,
    perm: np.ndarray,
    x_inicial: np.ndarray,
    pasadas: int = 2,
) -> List[Tuple[int, np.ndarray, mp.mpf]]:
    """
    Ejecuta refinamiento iterativo:
    1. Calcula r^(m) = b - A x^(m) en precisión extendida (mpmath 80 dígitos).
    2. Resuelve A d^(m) = r^(m) reutilizando L, U y perm en float64.
    3. Actualiza x^(m+1) = x^(m) + d^(m).
    """
    historial = []
    x_actual = x_inicial.copy()
    _, n_inf_0, _ = calcular_residuo_extendido(filas_a_str, filas_b_str, x_actual)
    historial.append((0, x_actual.copy(), n_inf_0))

    for p in range(1, pasadas + 1):
        r_mp, _, _ = calcular_residuo_extendido(filas_a_str, filas_b_str, x_actual)
        r_float = np.array([float(val) for val in r_mp], dtype=np.float64)
        delta_x = resolver_con_lu(L, U, perm, r_float)
        x_actual = x_actual + delta_x
        _, n_inf_p, _ = calcular_residuo_extendido(filas_a_str, filas_b_str, x_actual)
        historial.append((p, x_actual.copy(), n_inf_p))

    return historial