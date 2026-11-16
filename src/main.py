from pathlib import Path
import numpy as np
import mpmath as mp

from conteo import generar_markdown_conteo
from gauss import eliminar_y_resolver
from verificacion import (
    leer_datos_texto,
    calcular_solucion_exacta_mp,
    evaluar_exactitud,
    refinamiento_iterativo,
)

def formatear_digitos_defendibles(valor: float, digitos: int) -> str:
    if digitos <= 0:
        return "0"
    if valor == 0.0:
        return f"{0.0:.{digitos - 1}f}"

    # Todas las incógnitas tienen 1 dígito entero distinto de cero (orden 10^-1 a 10^0)
    abs_val = abs(valor)
    if abs_val >= 1.0:
        decimales = max(0, digitos - 1)
    else:
        # Ej: 0.8 tiene primer dígito significativo en la 1ra posición decimal
        orden = int(np.floor(np.log10(abs_val)))
        decimales = max(0, digitos - 1 - orden)
    return f"{valor:.{decimales}f}"


def main() -> None:
    raiz = Path(__file__).resolve().parent.parent
    ruta_a = raiz / "datos" / "matriz_a.csv"
    ruta_b = raiz / "datos" / "vector_b.csv"
    dir_res = raiz / "resultados"
    dir_res.mkdir(exist_ok=True)

    # 1. Lectura exacta de datos
    filas_a_str, filas_b_str = leer_datos_texto(str(ruta_a), str(ruta_b))
    A = np.array([[float(val) for val in fila] for fila in filas_a_str], dtype=np.float64)
    b = np.array([float(val) for val in filas_b_str], dtype=np.float64)

    # 2. Ejecución de las estrategias de eliminación
    res_sin = eliminar_y_resolver(A, b, estrategia="sin_pivoteo")
    res_parc = eliminar_y_resolver(A, b, estrategia="parcial")
    res_esc = eliminar_y_resolver(A, b, estrategia="escalado")

    # 3. Verificación en precisión extendida y control con biblioteca
    x_exacta_mp = calcular_solucion_exacta_mp(filas_a_str, filas_b_str)
    x_numpy = np.linalg.solve(A, b)

    eval_sin = evaluar_exactitud(filas_a_str, filas_b_str, res_sin.x, x_exacta_mp)
    eval_parc = evaluar_exactitud(filas_a_str, filas_b_str, res_parc.x, x_exacta_mp)
    eval_esc = evaluar_exactitud(filas_a_str, filas_b_str, res_esc.x, x_exacta_mp)
    eval_np = evaluar_exactitud(filas_a_str, filas_b_str, x_numpy, x_exacta_mp)

    # 4. Refinamiento iterativo sobre la solución con pivoteo escalado
    hist_ref = refinamiento_iterativo(
        filas_a_str, filas_b_str, res_esc.L, res_esc.U, res_esc.perm, res_esc.x, pasadas=2
    )
    x_refinado = hist_ref[-1][1]
    eval_ref = evaluar_exactitud(filas_a_str, filas_b_str, x_refinado, x_exacta_mp)

    # 5. Generar resultados/conteo_flops.md
    md_conteo = generar_markdown_conteo(res_sin.contador, res_parc.contador, res_esc.contador)
    (dir_res / "conteo_flops.md").write_text(md_conteo, encoding="utf-8")

    # 6. Generar resultados/traza_eliminacion.txt
    texto_traza = []
    texto_traza.append("================================================================================")
    texto_traza.append("TRAZA COMPLETA DE ELIMINACIÓN GAUSSIANA (11 PASOS)")
    texto_traza.append("================================================================================\n")
    texto_traza.append(">>> CORRIDA 1: PIVOTEO PARCIAL ESCALADO (ESTRATEGIA ÓPTIMA) <<<\n")
    texto_traza.extend(res_esc.traza_pasos)
    max_sub_esc = float(np.max(np.abs(np.tril(res_esc.U, -1))))
    texto_traza.append(f"Verificación triangular superior (máximo valor bajo la diagonal): {max_sub_esc:.1e}\n\n")

    texto_traza.append(">>> CORRIDA 2: SIN PIVOTEO (COMPARACIÓN DE PÉRDIDA DE MANTISA) <<<\n")
    texto_traza.extend(res_sin.traza_pasos)
    max_sub_sin = float(np.max(np.abs(np.tril(res_sin.U, -1))))
    texto_traza.append(f"Verificación triangular superior (máximo valor bajo la diagonal): {max_sub_sin:.1e}\n")

    (dir_res / "traza_eliminacion.txt").write_text("\n".join(texto_traza), encoding="utf-8")

    # 7. Generar resultados/solucion.csv
    lineas_sol = [
        "incognita,valor_pivoteo_escalado,digitos_correctos_escalado,valor_refinado,digitos_correctos_refinado,valor_sin_pivoteo,digitos_correctos_sin_pivoteo,error_relativo_escalado,error_relativo_sin_pivoteo"
    ]
    for i in range(12):
        d_esc = eval_esc.digitos_correctos_comp[i]
        d_ref = eval_ref.digitos_correctos_comp[i]
        d_sin = eval_sin.digitos_correctos_comp[i]

        val_esc_str = formatear_digitos_defendibles(float(res_esc.x[i]), d_esc)
        val_ref_str = formatear_digitos_defendibles(float(x_refinado[i]), d_ref)
        val_sin_str = formatear_digitos_defendibles(float(res_sin.x[i]), d_sin)

        lineas_sol.append(
            f"x{i+1},{val_esc_str},{d_esc},{val_ref_str},{d_ref},{val_sin_str},{d_sin},"
            f"{eval_esc.errores_relativos_comp[i]:.4e},{eval_sin.errores_relativos_comp[i]:.4e}"
        )
    (dir_res / "solucion.csv").write_text("\n".join(lineas_sol) + "\n", encoding="utf-8")

    # 8. Generar resultados/residuo.txt
    cond_2 = float(np.linalg.cond(A, 2))
    cond_inf = float(np.linalg.cond(A, np.inf))
    lineas_res = [
        "======================================================================",
        "VERIFICACIÓN EN PRECISIÓN EXTENDIDA (mpmath 80 dígitos) Y REFINAMIENTO",
        "======================================================================",
        f"Condicionamiento kappa_2(A)   : {cond_2:.6e}",
        f"Condicionamiento kappa_inf(A) : {cond_inf:.6e}",
        "",
        "--- 1. Residuo r = b - A*x (Pivoteo Parcial Escalado) ---",
    ]
    for i, r_val in enumerate(eval_esc.residuo_ext):
        lineas_res.append(f"r[{i+1:02d}] = {mp.nstr(r_val, 25)}")
    lineas_res.extend([
        f"Norma infinito ||r||_inf : {mp.nstr(eval_esc.norma_inf_residuo, 15)}",
        f"Norma 2        ||r||_2   : {mp.nstr(eval_esc.norma_2_residuo, 15)}",
        f"Error relativo ||dx||/||x||_2 : {eval_esc.error_relativo_norm2:.6e}",
        "",
        "--- 2. Residuo r = b - A*x (Sin Pivoteo) ---",
        f"Norma infinito ||r||_inf : {mp.nstr(eval_sin.norma_inf_residuo, 15)}",
        f"Norma 2        ||r||_2   : {mp.nstr(eval_sin.norma_2_residuo, 15)}",
        f"Error relativo ||dx||/||x||_2 : {eval_sin.error_relativo_norm2:.6e}",
        "",
        "--- 3. Residuo r = b - A*x (Pivoteo Parcial Estándar) ---",
        f"Norma infinito ||r||_inf : {mp.nstr(eval_parc.norma_inf_residuo, 15)}",
        f"Error relativo ||dx||/||x||_2 : {eval_parc.error_relativo_norm2:.6e}",
        "",
        "--- 4. Evolución del Refinamiento Iterativo (sobre Pivoteo Escalado) ---",
    ])
    for paso_ref, x_it, n_inf_it in hist_ref:
        ev_it = evaluar_exactitud(filas_a_str, filas_b_str, x_it, x_exacta_mp)
        lineas_res.append(
            f"Iteración {paso_ref}: ||r||_inf = {mp.nstr(n_inf_it, 12)} | "
            f"Error relativo ||dx||/||x||_2 = {ev_it.error_relativo_norm2:.6e}"
        )

    lineas_res.extend([
        "",
        "--- 5. Contraste contra rutina de biblioteca (numpy.linalg.solve) ---",
        f"Diferencia ||x_escalado - x_numpy||_inf : {np.max(np.abs(res_esc.x - x_numpy)):.6e}",
        f"Diferencia ||x_refinado - x_numpy||_inf : {np.max(np.abs(x_refinado - x_numpy)):.6e}",
        f"Error relativo numpy ||dx||/||x||_2     : {eval_np.error_relativo_norm2:.6e}",
    ])
    (dir_res / "residuo.txt").write_text("\n".join(lineas_res) + "\n", encoding="utf-8")

    print("Ejecución completada exitosamente. Archivos generados en 'resultados/'.")


if __name__ == "__main__":
    main()