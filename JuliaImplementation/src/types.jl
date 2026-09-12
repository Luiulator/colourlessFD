"""
Parámetros físicos y geométricos de la simulación del fluido.
Tipado con "T<:AbstractFloat".
"""
struct FluidParams{T<:AbstractFloat}
    nx::Int     # Número de celdas útiles en el eje horizontal (X)
    ny::Int     # Número de celdas útiles en el eje vertical (Y)
    dx::T       # Espaciado espacial entre celdas (Δx)
    dt::T       # Paso de tiempo entre iteraciones (Δt, ej. 0.016s ≈ 60 FPS)
    visc::T     # Viscosidad cinemática del fluido (ν). 0 = fluido no viscoso / ideal
    diff::T     # Coeficiente de difusión del tinte o humo
end

"""
Estructura principal que almacena todos los campos escalares y vectoriales de la malla 2D.

Usa un Double Buffering:
Cada propiedad tiene un array actual y uno temporal (`_prev`) para poder leer los valores
del instante anterior sin sobrescribirlos mientras se calculan los nuevos en el mismo paso.
"""
mutable struct FluidGrid2D{T<:AbstractFloat}
    # Campo de velocidad horizontal (eje X)
    u::Matrix{T}            # Velocidad X actual en cada celda
    u_prev::Matrix{T}       # Búfer temporal / velocidad X del paso anterior

    # Campo de velocidad vertical (eje Y)
    v::Matrix{T}            # Velocidad Y actual en cada celda
    v_prev::Matrix{T}       # Búfer temporal / velocidad Y del paso anterior

    # Campo escalar para visualización
    density::Matrix{T}      # Densidad actual de humo en cada celda
    density_prev::Matrix{T} # Búfer temporal de densidad

    # Búferes auxiliares para el cálculo de la presión
    p::Matrix{T}            # Campo de presión p(x, y) que compensa la divergencia
    div::Matrix{T}          # Divergencia del campo de velocidad intermedio (∇ · u*)

    # Parámetros físicos asociados a esta malla
    params::FluidParams{T}
end

"""
    FluidGrid2D([::Type{T}=Float32], nx::Int, ny::Int; dx=1.0f0, dt=0.016f0, visc=0.0f0, diff=0.0f0)
    FluidGrid2D{T}(nx::Int, ny::Int; dx=1.0, dt=0.016, visc=0.0, diff=0.0)

Constructor que inicializa la malla de fluidos preasignando toda la memoria necesaria.

¿Por qué `nx + 2` y `ny + 2`?
Se añaden 2 celdas extra por dimensión llamadas Celdas Fantasma:
- Las celdas interiores de simulación van de `2:nx+1` y `2:ny+1`.
- Los índices `1` y `nx+2` (en X), y `1` y `ny+2` (en Y) representan las fronteras o paredes
  del contenedor donde se aplican las condiciones de contorno.
"""
function FluidGrid2D{T}(nx::Int, ny::Int; dx=1.0, dt=0.016, visc=0.0, diff=0.0) where {T<:AbstractFloat}
    alloc() = zeros(T, nx + 2, ny + 2)
    params = FluidParams{T}(nx, ny, T(dx), T(dt), T(visc), T(diff))
    return FluidGrid2D{T}(
        alloc(), alloc(),  # u, u_prev
        alloc(), alloc(),  # v, v_prev
        alloc(), alloc(),  # density, density_prev
        alloc(), alloc(),  # p, div
        params
    )
end

function FluidGrid2D(::Type{T}, nx::Int, ny::Int; dx=1.0, dt=0.016, visc=0.0, diff=0.0) where {T<:AbstractFloat}
    return FluidGrid2D{T}(nx, ny; dx=dx, dt=dt, visc=visc, diff=diff)
end

function FluidGrid2D(nx::Int, ny::Int; dx=1.0f0, dt=0.016f0, visc=0.0f0, diff=0.0f0)
    return FluidGrid2D{Float32}(nx, ny; dx=dx, dt=dt, visc=visc, diff=diff)
end