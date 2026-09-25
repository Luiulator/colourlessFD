# src/boundary.jl

"""
    set_bnd!(x::Matrix{T}, b::Int, params::FluidParams{T}) where {T}

Aplica las condiciones de contorno de **pared sólida (Free-Slip / No-Penetration)**
en las celdas fantasma (*ghost cells*) que rodean la malla de simulación.

### Parámetro `b` (Tipo de campo):
- `b = 0`: **Escalar** (Densidad / Presión). Flujo nulo a través de la pared (Neumann homogéneo: ∂x/∂n = 0).
- `b = 1`: **Velocidad Horizontal `u`**. La componente normal se anula en las paredes verticales izquierda/derecha.
- `b = 2`: **Velocidad Vertical `v`**. La componente normal se anula en el suelo y en el techo.

### Geometría física de la pared:
1. **Paredes laterales (Izquierda `i=1` y Derecha `i=nx+2`)**:
   - Si `b == 1` (`u`), invertimos el signo (`-x`) para que la velocidad en la frontera física exacta (el punto medio entre celda 1 y 2) sea `0`.
   - Si `b == 0` o `b == 2`, copiamos el valor adyacente para permitir deslizamiento natural.

2. **Suelo y Techo (Fondo `j=1` y Techo `j=ny+2`)**:
   - Si `b == 2` (`v`), invertimos el signo (`-x`) para que el fluido no atraviese el suelo ni el techo.
   - Si `b == 0` o `b == 1`, copiamos el valor adyacente para que el fluido pueda fluir horizontalmente a lo largo del suelo.

3. **Esquinas**:
   - Se promedian con las dos celdas frontera adyacentes para evitar discontinuidades.
"""
function set_bnd!(x::Matrix{T}, b::Int, params::FluidParams{T}) where {T}
    nx, ny = params.nx, params.ny

    # 1. Paredes Izquierda (i = 1) y Derecha (i = nx + 2)
    @inbounds for j in 2:ny+1
        x[1, j]      = (b == 1) ? -x[2, j]      : x[2, j]   # Izquierda  
        x[nx+2, j]   = (b == 1) ? -x[nx+1, j]   : x[nx+1, j]# Derecha
    end

    # 2. Suelo (j = 1) y Techo (j = ny + 2)
    @inbounds for i in 2:nx+1
        x[i, 1]      = (b == 2) ? -x[i, 2]      : x[i, 2]
        x[i, ny+2]   = (b == 2) ? -x[i, ny+1]   : x[i, ny+1]
    end

    # 3. Las cuatro esquinas (promedio de las celdas adyacentes)
    @inbounds begin
        x[1, 1]          = T(0.5) * (x[2, 1]          + x[1, 2])
        x[1, ny+2]       = T(0.5) * (x[2, ny+2]       + x[1, ny+1])
        x[nx+2, 1]       = T(0.5) * (x[nx+1, 1]       + x[nx+2, 2])
        x[nx+2, ny+2]    = T(0.5) * (x[nx+1, ny+2]    + x[nx+2, ny+1])
    end
end
