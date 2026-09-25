# src/advection.jl

"""
    advect!(d::Matrix{T}, d0::Matrix{T}, u::Matrix{T}, v::Matrix{T}, b::Int, params::FluidParams{T}) where {T}

Realiza la **Advección Semi-Lagrangiana** incondicionalmente estable (Método de Jos Stam).

### ¿Qué hace?
Transporta una cantidad `d0` (que puede ser la densidad o las propias componentes de velocidad `u`, `v`)
a través del campo de velocidades `(u, v)`, guardando el resultado en `d`.

### Algoritmo:
1. Para cada celda `(i, j)`, retrocede en el tiempo siguiendo la línea de corriente:
   `x_prev = i - dt * u[i, j] / dx`
   `y_prev = j - dt * v[i, j] / dx`
2. Clampea la posición para no salir del dominio de celdas interiores `[1.5, N + 0.5]`.
3. Interpola bilinealmente el valor de `d0` a partir de las 4 celdas vecinas más cercanas.
4. Aplica las condiciones de contorno con `set_bnd!(d, b, params)`.
"""
function advect!(d::Matrix{T}, d0::Matrix{T}, u::Matrix{T}, v::Matrix{T}, b::Int, params::FluidParams{T}) where {T}
    nx, ny = params.nx, params.ny
    dt0 = params.dt / params.dx  # Factor de escala temporal respecto al tamaño de celda

    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            # 1. Retroceder a lo largo de la velocidad para encontrar la posición origen
            x = T(i) - dt0 * u[i, j]
            y = T(j) - dt0 * v[i, j]

            # 2. Clampear para que no se salga de las celdas interiores
            x = clamp(x, T(1.5), T(nx) + T(0.5))
            y = clamp(y, T(1.5), T(ny) + T(0.5))

            # 3. Índices de las 4 celdas circundantes
            i0 = floor(Int, x)
            i1 = i0 + 1
            j0 = floor(Int, y)
            j1 = j0 + 1

            # Pesos de interpolación (fracciones de distancia)
            s1 = x - T(i0)
            s0 = one(T) - s1
            t1 = y - T(j0)
            t0 = one(T) - t1

            # 4. Interpolación bilineal
            d[i, j] = s0 * (t0 * d0[i0, j0] + t1 * d0[i0, j1]) +
                      s1 * (t0 * d0[i1, j0] + t1 * d0[i1, j1])
        end
    end

    # 5. Aplicar condiciones de contorno
    set_bnd!(d, b, params)
end
