from dataclasses import dataclass, field
from typing import List


@dataclass
class PasoConteo:
    paso: int
    mult_div: int = 0
    sum_rest: int = 0
    comparaciones: int = 0
    intercambio: int = 0

    @property
    def total_flops(self) -> int:
        return self.mult_div + self.sum_rest


@dataclass
class ContadorOperaciones:
    estrategia: str
    pasos: List[PasoConteo] = field(default_factory=list)
    sust_mult_div: int = 0
    sust_sum_rest: int = 0
    escala_comparaciones: int = 0
    escala_divisiones: int = 0

    def iniciar_paso(self, k: int) -> PasoConteo:
        paso = PasoConteo(paso=k)
        self.pasos.append(paso)
        return paso

    @property
    def elim_mult_div(self) -> int:
        return sum(p.mult_div for p in self.pasos)

    @property
    def elim_sum_rest(self) -> int:
        return sum(p.sum_rest for p in self.pasos)

    @property
    def total_mult_div(self) -> int:
        return self.elim_mult_div + self.sust_mult_div

    @property
    def total_sum_rest(self) -> int:
        return self.elim_sum_rest + self.sust_sum_rest

    @property
    def total_flops(self) -> int:
        return self.total_mult_div + self.total_sum_rest

    @property
    def total_comparaciones_pivoteo(self) -> int:
        return sum(p.comparaciones for p in self.pasos) + self.escala_comparaciones

    @property
    def total_intercambios(self) -> int:
        return sum(p.intercambio for p in self.pasos)


def generar_markdown_conteo(
    c_sin: ContadorOperaciones,
    c_parc: ContadorOperaciones,
    c_esc: ContadorOperaciones,
) -> str:
    lineas = [
        "# Conteo de Operaciones Elementales (FLOPs)",
        "",
        "## Convención Declarada",
        "- **1 FLOP** = 1 operación aritmética de punto flotante (`+`, `-`, `*`, `/`).",
        "- En cada paso $k \\in \\{1, \\dots, 11\\}$ y para cada fila $i > k$:",
        "  - Se calcula **una sola vez** el multiplicador $m_{ik} = a_{ik} / a_{kk}$ (1 división).",
        "  - Se fija $a_{ik} = 0$ directamente sin operar la subdiagonal ya anulada.",
        "  - Se actualizan las columnas $j \\in \\{k+1, \\dots, n\\}$ con $a_{ij} \\leftarrow a_{ij} - m_{ik} a_{kj}$ (1 multiplicación y 1 resta por columna).",
        "  - Se actualiza el vector derecho $b_i \\leftarrow b_i - m_{ik} b_k$ (1 multiplicación y 1 resta).",
        "- El costo de búsqueda de pivote, cálculo de factores de escala e intercambios se contabiliza aparte como **sobrecosto de estabilidad**.",
        "",
        "## Desglose Paso a Paso de la Eliminación ($n = 12$)",
        "",
        "| Paso $k$ | Filas eliminadas | Mult / Div | Sum / Rest | Total FLOPs Paso | Comp. Pivoteo (Parcial) | Intercambio (Escalado) |",
        "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|",
    ]

    for p_s, p_p, p_e in zip(c_sin.pasos, c_parc.pasos, c_esc.pasos):
        filas = 12 - p_s.paso
        lineas.append(
            f"| {p_s.paso} | {filas} | {p_s.mult_div} | {p_s.sum_rest} | {p_s.total_flops} | {p_p.comparaciones} | {p_e.intercambio} |"
        )

    lineas.extend([
        f"| **Subtotal Eliminación** | **66** | **{c_sin.elim_mult_div}** | **{c_sin.elim_sum_rest}** | **{c_sin.elim_mult_div + c_sin.elim_sum_rest}** | **{sum(p.comparaciones for p in c_parc.pasos)}** | **{c_esc.total_intercambios}** |",
        f"| **Sustitución Regresiva** | - | **{c_sin.sust_mult_div}** | **{c_sin.sust_sum_rest}** | **{c_sin.sust_mult_div + c_sin.sust_sum_rest}** | 0 | 0 |",
        f"| **TOTAL GENERAL** | - | **{c_sin.total_mult_div}** | **{c_sin.total_sum_rest}** | **{c_sin.total_flops}** | **66** | **{c_esc.total_intercambios}** |",
        "",
        "## Comparación de Estrategias y Cotas de Referencia",
        "",
        "| Estrategia | Mult / Div | Sum / Rest | Total FLOPs Aritméticos | Cota de Referencia | Diferencia vs. 1486 | Sobrecosto de Estabilidad |",
        "|:---|:---:|:---:|:---:|:---:|:---:|:---|",
        f"| Sin pivoteo | {c_sin.total_mult_div} | {c_sin.total_sum_rest} | **{c_sin.total_flops}** | 1486 | {_diferencia(c_sin.total_flops)} | 0 comparaciones, 0 intercambios |",
        f"| Pivoteo parcial | {c_parc.total_mult_div} | {c_parc.total_sum_rest} | **{c_parc.total_flops}** | 1486 | {_diferencia(c_parc.total_flops)} | {c_parc.total_comparaciones_pivoteo} comparaciones, {c_parc.total_intercambios} intercambios |",
        f"| Pivoteo parcial escalado | {c_esc.total_mult_div} | {c_esc.total_sum_rest} | **{c_esc.total_flops}** | 1486 | {_diferencia(c_esc.total_flops)} | {c_esc.escala_comparaciones} comp. escala + {c_esc.escala_divisiones} div. cocientes + {sum(p.comparaciones for p in c_esc.pasos)} comp. pivote, {c_esc.total_intercambios} intercambios |",
        "",
        "El total de 1354 FLOPs es el mínimo estricto bajo la convención declarada "
        "(eliminación + actualización de b + sustitución regresiva). La cota de 1486 "
        "es la referencia del enunciado; la diferencia se calcula respecto de ella.",
    ])
    return "\n".join(lineas) + "\n"


def _diferencia(total: int, referencia: int = 1486) -> str:
    delta = total - referencia
    return f"{delta:+d} ({100.0 * delta / referencia:+.2f}%)"
