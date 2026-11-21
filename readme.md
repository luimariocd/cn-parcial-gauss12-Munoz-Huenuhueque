# Eliminación Gaussiana 12×12: costo operacional mínimo y control de mantisa

Este proyecto resuelve un sistema de ecuaciones lineales de 12 incógnitas utilizando el método de **eliminación gaussiana con pivoteo parcial escalado**.

El programa realiza los cálculos automáticamente, registra las operaciones efectuadas, contabiliza el costo computacional (FLOPs) y verifica la precisión de la solución obtenida mediante el cálculo del residuo en precisión extendida.

El objetivo es estudiar el comportamiento numérico del método y analizar su costo operacional.

## Autoría y asignatura

- **Estudiantes:** Nicolás Muñoz y Luismario Huenuhueque.
- **Asignatura:** Computación Numérica (INF-321).
- **Carrera:** Ingeniería Civil Informática.
- **Institución:** Universidad Católica del Maule (UCM).
- **Docentes:** Juan Manuel Cabezas Contardo / Pedro Pesante.
- **Fecha:** Noviembre de 2026.

---

## 1. Requisitos previos

Para utilizar este proyecto no es necesario tener conocimientos de programación, pero sí contar con los siguientes programas instalados:

- **Python 3.14:** es el lenguaje utilizado para ejecutar el programa. Se recomienda utilizar la versión 3.14.7, con la que se realizaron las pruebas.
- **pip:** herramienta que permite instalar las bibliotecas necesarias para ejecutar el proyecto. Generalmente viene incluida con Python.

### ¿Cómo instalar Python?

1. Ingresar al sitio oficial: https://www.python.org/downloads/
2. Descargar Python 3.14 para el sistema operativo correspondiente.
3. Ejecutar el instalador y completar la instalación.
4. En Windows, marcar la opción **Add Python to PATH** antes de continuar con la instalación, si aparece disponible.

Para comprobar que Python está instalado correctamente, abrir una terminal y ejecutar:

```bash
python --version
```

En Windows también se puede utilizar:

```bash
py --version
```

Si la instalación fue correcta, se mostrará la versión de Python instalada.

---

## 2. Descargar el proyecto

No es necesario instalar Git para obtener los archivos del proyecto.

1. Ingresar al repositorio del proyecto en GitHub.
2. Presionar el botón verde **Code**.
3. Seleccionar la opción **Download ZIP**.
4. Una vez descargado el archivo, descomprimirlo en una ubicación de preferencia.
5. Abrir la carpeta del proyecto.

Después de estos pasos, se tendrán disponibles todos los archivos necesarios para ejecutar el programa.

---

## 3. Instalar las dependencias

El proyecto utiliza dos bibliotecas externas de Python que permiten realizar los cálculos numéricos y las verificaciones correspondientes.

Las dependencias son:

| Biblioteca | Versión | Utilidad |
|---|---|---|
| NumPy | 2.5.4 | Operaciones y cálculos numéricos. |
| mpmath | 1.4.1 | Cálculos con precisión extendida. |

Para instalarlas, primero es necesario abrir una terminal dentro de la carpeta principal del proyecto.

**En Windows:**

1. Abrir la carpeta del proyecto.
2. Hacer clic en la barra de direcciones del Explorador de archivos.
3. Escribir `powershell` y presionar Enter.
4. Se abrirá una terminal ubicada directamente en la carpeta del proyecto.

Luego ejecutar:

```bash
py -3.14 -m pip install -r requirements.txt
```

**En Linux o macOS:**

Abrir una terminal, ingresar a la carpeta del proyecto y ejecutar:

```bash
python3.14 -m pip install -r requirements.txt
```

Este comando utiliza el archivo `requirements.txt` para instalar automáticamente las bibliotecas y versiones necesarias.

No es necesario instalar NumPy ni mpmath por separado.

---

## 4. Ejecutar el programa

Una vez instaladas las dependencias, el proyecto está listo para ejecutarse.

Desde la terminal abierta en la carpeta principal, ingresar el siguiente comando:

**Windows:**

```bash
py -3.14 src/main.py
```

**Linux o macOS:**

```bash
python3.14 src/main.py
```

Al ejecutar el comando, el programa realizará automáticamente el proceso de eliminación gaussiana, contabilizará las operaciones realizadas y verificará la solución obtenida.

No es necesario ingresar datos manualmente, ya que el programa utiliza las matrices y los vectores almacenados en la carpeta `datos`.

### Resultados generados

Al finalizar la ejecución, se generarán o actualizarán automáticamente los siguientes archivos dentro de la carpeta `resultados`:

| Archivo | Contenido |
|---|---|
| `solucion.csv` | Valores obtenidos para las 12 incógnitas del sistema. |
| `traza_eliminacion.txt` | Registro del procedimiento de eliminación gaussiana. |
| `conteo_flops.md` | Desglose del conteo de operaciones de punto flotante (FLOPs). |
| `residuo.txt` | Resultados de la verificación del residuo y la precisión de la solución. |

Estos archivos permiten revisar tanto la solución final como el procedimiento y las verificaciones realizadas.

---

## 5. Estructura del proyecto

El proyecto está organizado en diferentes carpetas para separar los datos de entrada, el código fuente y los resultados obtenidos.

```text
cn-parcial-gauss12-Munoz-Huenuhueque/
│
├── readme.md
├── informe.md
├── requirements.txt
├── .gitignore
│
├── datos/
│   ├── matriz_a.csv
│   └── vector_b.csv
│
├── src/
│   ├── gauss.py
│   ├── conteo.py
│   ├── verificacion.py
│   └── main.py
│
└── resultados/
    ├── solucion.csv
    ├── traza_eliminacion.txt
    ├── conteo_flops.md
    └── residuo.txt
```

### Descripción de las carpetas y archivos

- **`datos/`**: contiene los archivos CSV con la matriz de coeficientes y el vector de términos independientes que utiliza el programa.
- **`src/`**: contiene el código fuente del proyecto.
  - `main.py`: archivo principal desde el cual se ejecuta el programa.
  - `gauss.py`: implementa el método de eliminación gaussiana.
  - `conteo.py`: realiza el conteo de operaciones de punto flotante (FLOPs).
  - `verificacion.py`: realiza las verificaciones numéricas de la solución.
- **`resultados/`**: almacena los archivos generados durante la ejecución.
- **`requirements.txt`**: contiene las bibliotecas y versiones necesarias para ejecutar el proyecto.
- **`informe.md`**: contiene el desarrollo y análisis del trabajo.
- **`.gitignore`**: especifica los archivos y carpetas que no deben incluirse en el control de versiones.

---

## 6. Resultados obtenidos

A continuación, se presentan los valores obtenidos para las 12 incógnitas del sistema lineal:

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

### Indicadores numéricos

- **Total de FLOPs:** 1354.
- **Norma infinito del residuo (solución directa):** $6.04 \times 10^{-10}$.
- **Norma infinito del residuo (solución refinada):** $2.60 \times 10^{-10}$.
- **Error forward después de una pasada de refinamiento:** $4.86 \times 10^{-17}$.

---

## 7. Informe del proyecto

Para conocer en detalle el desarrollo del método, el análisis del costo operacional, el pivoteo, el condicionamiento y las verificaciones numéricas, se puede consultar el siguiente documento:

[**Ver informe completo del proyecto (informe.md)**](informe.md)

---

## 8. Licencia y uso académico

Trabajo académico desarrollado para la asignatura Computación Numérica (INF-321) de la carrera de Ingeniería Civil Informática de la Universidad Católica del Maule.

