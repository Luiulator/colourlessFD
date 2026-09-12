# src/diffusion.jl

"""
    lin_solve!(x::Matrix{T}, x0::Matrix{T}, a::T, c::T, b::Int, params::FluidParams{T}; iter::Int=20) where {T}

Resuelve el sistema lineal elíptico (I - a ∇²)x = x0 usando Gauss-Seidel.

En cada celda (i, j), la discretización por diferencias finitas es:
    c · x[i, j] = x0[i, j] + a · (x[i-1, j] + x[i+1, j] + x[i, j-1] + x[i, j+1])
"""
function lin_solve!(x::Matrix{T}, x0::Matrix{T}, a::T, c::T, b::Int, params::FluidParams{T}; iter::Int=20) where {T}
    nx, ny = params.nx, params.ny
    inv_c = one(T) / c

    for _ in 1:iter
        for j in 2:ny+1
            @inbounds for i in 2:nx+1
                neighbors_sum = x[i-1, j] + x[i+1, j] + x[i, j-1] + x[i, j+1]
                x[i, j] = (x0[i, j] + a * neighbors_sum) * inv_c
            end
        end
        # Aplicar condiciones de contorno tras cada barrido de Gauss-Seidel
        set_bnd!(x, b, params)
    end
end

"""
    diffuse!(x::Matrix{T}, x0::Matrix{T}, diff_rate::T, b::Int, params::FluidParams{T}; iter::Int=20) where {T}

Aplica la **Difusión Implícita** (ecuación del calor):
    ∂x/∂t = ν ∇²x  ==>  (I - ν Δt ∇²) x_(n+1) = x_n

### Parámetros:
- `x`: Matriz de destino donde se guarda el nuevo estado difundido.
- `x0`: Matriz de origen con los valores antes de la difusión.
- `diff_rate`: Coeficiente de difusión (viscosidad `visc` para velocidad, o `diff` para tinte/humo).
- `b`: Tipo de frontera (`0` escalar, `1` velocidad X, `2` velocidad Y).
- `iter`: Número de iteraciones de Gauss-Seidel (20 es el estándar en Stable Fluids).
"""
function diffuse!(x::Matrix{T}, x0::Matrix{T}, diff_rate::T, b::Int, params::FluidParams{T}; iter::Int=20) where {T}
    # Si no hay difusión, simplemente copiamos los datos
    if diff_rate <= zero(T)
        x .= x0
        return
    end

    # Constante a = diff_rate * dt / (dx * dx)
    a = params.dt * diff_rate / (params.dx * params.dx)
    c = one(T) + T(4) * a

    lin_solve!(x, x0, a, c, b, params; iter=iter)
end
