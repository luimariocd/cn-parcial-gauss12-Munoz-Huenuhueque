# Eliminación Gaussiana 12×12: costo operacional mínimo y control de mantisa

Resuelve un sistema lineal denso de 12 incógnitas mediante eliminación gaussiana con pivoteo parcial escalado, conteo instrumentado de FLOPs y verificación del residuo en precisión extendida.

## Autoría y asignatura

- **Estudiantes:** Nicolás Muñoz, Luismario Huenuhueque
- **Asignatura:** Computación Numérica (INF-321) — Ingeniería Civil Informática, Universidad Católica del Maule (UCM)
- **Docentes:** Juan Manuel Cabezas Contardo / Pedro Pesante
- **Fecha:** Noviembre de 2026

## Requisitos

- **Python:** 3.14 (probado en Python 3.14.7).
- **Instalación de dependencias:**

  ```bash
  pip install -r requirements.txt
  ```

  Dependencias fijadas: `numpy==2.5.4` y `mpmath==1.4.1`.

## Ejecución

Desde la raíz del repositorio:

```bash
python src/main.py
```

El comando genera o actualiza automáticamente:

- `resultados/solucion.csv`
- `resultados/traza_eliminacion.txt`
- `resultados/conteo_flops.md`
- `resultados/residuo.txt`

## Estructura de carpetas

```text
cn-parcial-gauss12-Munoz-Huenuhueque/
├─ readme.md
├─ informe.md
├─ requirements.txt
├─ datos/
│  ├─ matriz_a.csv
│  └─ vector_b.csv
├─ src/
│  ├─ gauss.py
│  ├─ conteo.py
│  ├─ verificacion.py
│  └─ main.py
├─ resultados/
│  ├─ solucion.csv
│  ├─ traza_eliminacion.txt
│  ├─ conteo_flops.md
│  └─ residuo.txt
└─ .gitignore
```

## Resultado resumido

| Incógnita | Valor |
|:---:|---:|
| $x_1$ | 0.7384915627 |
| $x_2$ | -1.2059374168 |
| $x_3$ | 0.0043871952 |
| $x_4$ | -0.0917264538 |
| $x_5$ | 2.8461093725 |
| $x_6$ | -0.0062815409 |
| $x_7$ | 1.4738205916 |
| $x_8$ | -0.0384927164 |
| $x_9$ | 0.9152648037 |
| $x_{10}$ | -2.6071394852 |
| $x_{11}$ | 0.0518736294 |
| $x_{12}$ | -1.0849371625 |

- **Total de FLOPs:** 1354
- **Norma infinito del residuo:** $6.04 \times 10^{-10}$ (solución directa); $2.60 \times 10^{-10}$ en la solución refinada.
- **Error forward tras una pasada de refinamiento:** $4.86 \times 10^{-17}$.

## Informe

El análisis completo de costo, pivoteo, condicionamiento y verificación está en [informe.md](informe.md).

## Licencia y uso académico

Trabajo académico desarrollado para la asignatura Computación Numérica (INF-321) de la Universidad Católica del Maule. Su uso es educativo; se permite consultar y reutilizar el código citando la fuente.
