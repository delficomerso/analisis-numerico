### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ b3234e54-32f3-4727-b9ee-35ad7911c849
begin
	import Pkg
	#Pkg.activate("ruta/a/DiferenciasFinitas") <- acá poner el directorio a la carpeta DiferenciasFinitas
    Pkg.instantiate()
end

# ╔═╡ 6e339a23-bfa0-43b1-b8d8-2f52cc55b94e
using DiferenciasFinitas, Revise, Plots, LinearAlgebra

# ╔═╡ 22e36d32-3b47-42e9-98b9-dc9086edd5f4
begin
	h = 0.05        
	Tf = 0.1 #no muy grande para que no explote
	u0 = cond_inicial(h)    
end

# ╔═╡ bad1a2e1-8c69-4907-8f63-27d1ae7fa13c
md"""
#### Método Explícito

El método explícito para la ecuación del calor en $2$D es condicionalmente estable.  
La condición de estabilidad es

$r = \frac{\Delta t}{h^2} \le \frac{1}{4}$

Por ejemplo, si tomamos $h = 0.05$:  
- Con $\Delta t = 0.0001$ (cumpliendo la condición), la solución se difunde suavemente y se ve correcta.  
- Con $\Delta t = 0.01$ (violando la condición), la solución "explota".
"""

# ╔═╡ 08bd6b2d-a85e-4b84-ace4-a9ee08c30af3
begin
	dt_estable = 0.0001
	Nsteps_estable = Int(Tf / dt_estable)
	local u = copy(u0)
	for n in 1:Nsteps_estable
	    u = paso_explicito(u,h,dt_estable)
	end
	
	 heatmap(u, title="Explícito estable: Δt = $dt_estable, h = $h", color=:viridis)
end

# ╔═╡ b5ece3c0-c0f8-4bff-9f96-21b3c98cf889
begin
    dt_inestable = 0.01
    Nsteps_inestable = Int(Tf / dt_inestable)
    local u = copy(u0)
    for n in 1:Nsteps_inestable
        u = paso_explicito(u, h, dt_inestable)
    end
    heatmap(u, title="Explícito inestable: Δt = $dt_inestable, h = $h", color=:viridis)
end


# ╔═╡ f131f0ec-b6c6-4008-9b70-b41e367f46bc
md"""
#### Método Crank–Nicolson

El método de Crank–Nicolson es incondicionalmente estable, así que las soluciones no deberían "explotar" para ningún $\Delta t$.  
Sin embargo, la precisión y el comportamiento cualitativo de la solución sí dependen del tamaño del paso temporal.

Por ejemplo, tomando $h = 0.05$:

- Con $\Delta t = 0.0001$, la solución es precisa, suave y se parece a la reconstruida por el método explícito.  

- Con $\Delta t = 0.1$, se tiene $r = 40$, un valor muy grande. Aunque el método sigue siendo estable, aparecen oscilaciones no físicas y valores negativos debido a la falta de amortiguamiento numérico.  

Esto no contradice la estabilidad incondicional del método, significa que para pasos de tiempo grandes, deja de representar correctamente el comportamiento suavizante de la ecuación del calor.
"""


# ╔═╡ e4d98edb-882d-47dc-9ef2-d7c99067e026
#este tarda un poco en correr
begin
	dt_CN_small = 0.0001
	Nsteps_CN_small = Int(Tf / dt_CN_small)
    local u = copy(u0)
    for n in 1:Nsteps_CN_small
        u = paso_CN(u,h,dt_CN_small)
    end
    heatmap(u, title="Solución C-N, Δt = $dt_CN_small, h = $h", color=:viridis)
end

# ╔═╡ 04be86b6-f80c-4c1d-8bdd-ccd5de2ec2df
begin
	dt_CN_big = 0.1
	Nsteps_CN_big = Int(Tf / dt_CN_big)
    local u = copy(u0)
    for n in 1:Nsteps_CN_big
        u = paso_CN(u,h,dt_CN_big)
    end
    heatmap(u, title="Solución C-N, Δt = $dt_CN_big, h = $h", color=:viridis)
end

# ╔═╡ 0932d7cc-cbcc-4ad9-a565-43cf4458ee61
md"""
#### Método ADI (Alternating Direction Implicit)

El método ADI también es incondicionalmente estable.

Por ejemplo, tomando $h = 0.05$:

- Con $\Delta t = 0.0001$, la solución es suave y muy similar a la obtenida con el método explícito y Crank–Nicolson.

- Con $\Delta t = 0.1$, se tiene $r = 40$. El método sigue siendo estable, y aunque la solución comienza a mostrar ligeras deformaciones, no aparecen valores negativos ni oscilaciones fuertes como en Crank–Nicolson. Esto indica que el esquema ADI controla mejor los errores numéricos para pasos de tiempo grandes.

"""

# ╔═╡ b328d224-02fe-4f70-8e43-8baf3daab59b
begin
    dt_ADI_small = 0.0001
	Nsteps_ADI_small = Int(Tf / dt_ADI_small)
    local u = copy(u0)
    for n in 1:Nsteps_ADI_small
        u = paso_ADI(u,h,dt_ADI_small)
    end
    heatmap(u, title="Solución ADI, Δt = $dt_ADI_small, h = $h", color=:viridis)
end

# ╔═╡ ae12a6fb-3f74-47cf-8910-99b26758a9a4
begin
	dt_ADI_big = 0.1
	Nsteps_ADI_big = Int(Tf / dt_ADI_big)
    local u = copy(u0)
    for n in 1:Nsteps_ADI_big
        u = paso_ADI(u,h,dt_ADI_big)
    end
    heatmap(u, title="Solución ADI, Δt = $dt_ADI_big, h = $h", color=:viridis)
end

# ╔═╡ 1c423d5d-3215-4231-a474-826a0572c1de
md"""
#### Comparación del costo computacional

Al comparar el desempeño del método ADI con el de Crank–Nicolson, se observa una diferencia importante en el costo computacional.  
El método de Crank–Nicolson requiere resolver en cada iteración un sistema lineal bidimensional de gran tamaño, lo que implica un costo computacional alto, especialmente cuando la malla espacial es fina.  

En cambio, el método ADI divide cada paso temporal en dos subpasos unidimensionales: uno en la dirección $x$ y otro en la dirección $y$. En lugar de resolver un único sistema grande, el método resuelve dos sistemas tridiagonales, lo que reduce el tiempo de cómputo.  

En la práctica, ambos métodos poseen la misma estabilidad y orden de precisión, pero ADI resulta mucho más eficiente en términos de tiempo de ejecución.  
"""

# ╔═╡ 04c2665d-1967-4a1d-84ec-9f2c6cf9265e
begin
	dt = 0.001
	Nsteps = Int(Tf/dt)
    local u = copy(u0)
    println("Tiempo C-N:")
    @time for n in 1:Nsteps
        u = paso_CN(u,h,dt)
    end

    local u = copy(u0)
    println("Tiempo ADI:")
    @time for n in 1:Nsteps
        u = paso_ADI(u,h,dt)
    end
end

# ╔═╡ Cell order:
# ╠═b3234e54-32f3-4727-b9ee-35ad7911c849
# ╠═6e339a23-bfa0-43b1-b8d8-2f52cc55b94e
# ╠═22e36d32-3b47-42e9-98b9-dc9086edd5f4
# ╟─bad1a2e1-8c69-4907-8f63-27d1ae7fa13c
# ╠═08bd6b2d-a85e-4b84-ace4-a9ee08c30af3
# ╠═b5ece3c0-c0f8-4bff-9f96-21b3c98cf889
# ╟─f131f0ec-b6c6-4008-9b70-b41e367f46bc
# ╠═e4d98edb-882d-47dc-9ef2-d7c99067e026
# ╠═04be86b6-f80c-4c1d-8bdd-ccd5de2ec2df
# ╟─0932d7cc-cbcc-4ad9-a565-43cf4458ee61
# ╠═b328d224-02fe-4f70-8e43-8baf3daab59b
# ╠═ae12a6fb-3f74-47cf-8910-99b26758a9a4
# ╟─1c423d5d-3215-4231-a474-826a0572c1de
# ╠═04c2665d-1967-4a1d-84ec-9f2c6cf9265e
