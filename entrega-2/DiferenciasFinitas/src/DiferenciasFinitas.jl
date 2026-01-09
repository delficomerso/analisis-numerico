module DiferenciasFinitas

using LinearAlgebra

function meshgrid(x,y)
    X = [xi for yi in y,xi in x]
    Y = [yi for yi in y,xi in x]
    return X,Y
end

#condición inicial con condiciones de borde Dirichlet
function cond_inicial(h)
    xs = 0:h:2
    ys = 0:h:2
    X,Y = meshgrid(xs,ys)
    u0 = exp.(-10*((X.-0.5).^2 .+ (Y.-0.5).^2))
    u0[1,:] .= 0; u0[end,:] .= 0
    u0[:,1] .= 0; u0[:,end] .= 0
    return u0
end

#discretizaciones de 2da derivada (sin dividir por h^2 porque uso r=t/h^2 por afuera)
function δ2x(u,h,i,j)
    uxx = (u[i+1,j]-2*u[i,j]+u[i-1,j])
    return uxx
end

function δ2y(u,h,i,j)
    uyy = (u[i,j+1]-2*u[i,j]+u[i,j-1])
    return uyy
end

#método explícito (un paso)
function paso_explicito(u,h,dt)
    r = dt/h^2
    J = size(u,1)-1 #J+1 nodos totales
    u_new = copy(u)
    for i in 2:J #completo los nodos interiores
        for j in 2:J
            uxx = δ2x(u,h,i,j)
            uyy = δ2y(u,h,i,j)
            u_new[i,j] = u[i,j] + r*(uxx+uyy)
        end
    end
    u_new[1, :] .= 0; u_new[end, :] .= 0
    u_new[:, 1] .= 0; u_new[:, end] .= 0
    return u_new
end 

#Laplaciano 1D
function matriz_1D(N) #N es el número de nodos interiores
    return Tridiagonal(ones(N-1),-2*ones(N),ones(N-1))
end

#Laplaciano 2D
function matrices_2D(J,r)
    Dx = matriz_1D(J-1) #J-1 nodos interiores 
    Dy = matriz_1D(J-1)
    I1 = Matrix(I,J-1,J-1)
    L = kron(I1,Dx) + kron(Dy,I1)
    A = Matrix(I,(J-1)^2,(J-1)^2) - r/2*L
    B = Matrix(I,(J-1)^2,(J-1)^2) + r/2*L

    return A,B
end

#método Crank-Nicholson (un paso)
function paso_CN(u,h,dt)
    r = dt/h^2
    J = size(u,1)-1 #J+1 nodos totales
    A,B = matrices_2D(J,r)
    u_vec = vec(u[2:end-1,2:end-1])
    b = B*u_vec
    u_new_vec = A\b
    u_new = reshape(u_new_vec,J-1,J-1)
    u_full = zeros(J+1,J+1)
    u_full[2:end-1,2:end-1] .= u_new
    return u_full
end

function matrices_ADI(J,r)
    D = matriz_1D(J-1)
    Dx = I - r/2*D
    Dy = I - r/2*D
    return Dx,Dy
end

#método ADI (un paso)
function paso_ADI(u,h,dt)
    r = dt/h^2
    J = size(u,1)-1 #J+1 nodos totales
    u_half = copy(u)
    u_new = copy(u)
    Dx, Dy = matrices_ADI(J,r)
    for j in 2:J # nodos interiores
        b = [u[i,j] + (r/2)*δ2y(u,h,i,j) for i in 2:J]
        u_half[2:J,j] = Dx\b
    end
    for i in 2:J
        b = [u_half[i,j] + (r/2)*δ2x(u_half,h,i,j) for j in 2:J]
        u_new[i,2:J] = Dy\b
    end
    u_new[1,:] .= 0; u_new[end,:] .= 0
    u_new[:,1] .= 0; u_new[:,end] .= 0
    return u_new
end 

export meshgrid, cond_inicial, δ2x, δ2y, matriz_1D, matrices_2D, matrices_ADI, paso_explicito, paso_CN, paso_ADI

end # module DiferenciasFinitas
