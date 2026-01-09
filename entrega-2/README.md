# Entrega 2 — Guía 2 (Ejercicios 6 y 7)

Esta entrega corresponde a los **Ejercicios 6 y 7 de la Guía de Laboratorio 2** de la materia **Análisis Numérico**.

## Descripción del problema

Se estudia numéricamente la **ecuación del calor bidimensional**:

U_t = U_xx + U_yy,   (x,y) ∈ Ω = [0,2]²,   t ∈ (0,T_f)

con condiciones de borde de Dirichlet homogéneas y condición inicial gaussiana.

El problema se discretiza mediante **diferencias finitas** en una grilla espacial uniforme
y una grilla temporal, y se aproxima la solución usando distintos esquemas temporales.

Los objetivos principales son:
- Implementar y comparar distintos métodos numéricos.
- Analizar la **estabilidad** de los métodos.
- Comparar los **tiempos de ejecución** entre métodos implícitos.

## Métodos implementados

- Método explícito
- Método de Crank–Nicholson
- Método ADI (Alternating Direction Implicit)

## Contenido

- `DiferenciasFinitas.jl/`: librería propia con la implementación de los métodos.
- `scripts/`: scripts para ejecutar los experimentos pedidos.
- `figures/`: gráficos de la solución numérica y comparaciones.
- `Entrega2.pdf`: consigna original de la entrega.

## Lenguaje

- Julia
