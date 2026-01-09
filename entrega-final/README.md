# Entrega Final — FEM 2D

Esta entrega corresponde al **Trabajo Final**.

## Descripción del problema

Se estudia un **problema de reacción–difusión singularmente perturbado** en dos dimensiones,
definido en el dominio Ω = [0,1]²:

−ε² Δu + u = f,   (x,y) ∈ Ω  
u = 0,           (x,y) ∈ ∂Ω

Se elige el término fuente f de modo que la solución exacta u(x,y) sea conocida,
permitiendo analizar el comportamiento numérico para distintos valores del parámetro ε.

El problema se resuelve mediante el **método de elementos finitos** usando mallas triangulares
uniformes y funciones lineales a trozos.

Los objetivos principales son:
- Implementar la generación de mallas triangulares en 2D.
- Resolver el problema FEM para distintos valores de ε y tamaños de malla.
- Analizar el fenómeno de **capas límite**.
- Calcular errores en las normas L² y H¹ y estimar el **orden de convergencia**.

## Contenido

- `EntregaFinal_Comerso.jl`: implementación en Julia del método de elementos finitos.
- `entrega_final.pdf`: consigna original del trabajo final.

## Lenguaje

- Julia
