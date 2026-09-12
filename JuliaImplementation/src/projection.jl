# src/projection.jl

"""
    project!(u::Matrix{T}, v::Matrix{T}, p::Matrix{T}, div::Matrix{T}, params::FluidParams{T}; iter::Int=20) where {T}

Algoritmo en 3 pasos (revisar Teorema de la Descomposición de Helmholtz-Hodge):

1. Cálculo de la Divergencia (∇ · w):
2. Resolución de la Ecuación de Poisson (∇²p = ∇ ⋅ w):
   Calcula la presión "p" necesaria en cada celda para equilibrar las compresiones/expansiones
   mediante el solver "lin_solve!" (Gauss-Seidel con a = 1, c = 4).
3. Sustracción del Gradiente de Presión:
   Corrige las velocidades restando la fuerza que empuja de zonas de alta presión a baja presión,
   garantizando que el campo final cumpla exactamente ∇ · u = 0.
"""
function project!(u::Matrix{T}, v::Matrix{T}, p::Matrix{T}, div::Matrix{T}, params::FluidParams{T}; iter::Int=20) where {T}
    nx, ny = params.nx, params.ny
    dx = params.dx

    # =========================================================================
    # Fase 1: Calcular la divergencia del campo de velocidad y poner p a cero
    # =========================================================================
    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            # Divergencia por diferencias finitas centrales
            div[i, j] = -T(0.5) * dx * (
                (u[i+1, j] - u[i-1, j]) +
                (v[i, j+1] - v[i, j-1])
            )
            # Inicializamos la presión estimada en cero
            p[i, j] = zero(T)
        end
    end

    # Aplicamos contornos a la divergencia y a la presión (b = 0 para escalares)
    set_bnd!(div, 0, params)
    set_bnd!(p, 0, params)

    # =========================================================================
    # Fase 2: Resolver la ecuación de Poisson de la presión: ∇²p = div
    # =========================================================================
    # Es decir: lin_solve! con a = 1, c = 4
    lin_solve!(p, div, one(T), T(4), 0, params; iter=iter)

    # =========================================================================
    # Fase 3: Restar el gradiente de presión a las velocidades (u = u* - ∇p)
    # =========================================================================
    inv_dx = T(0.5) / dx

    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            u[i, j] -= inv_dx * (p[i+1, j] - p[i-1, j])
            v[i, j] -= inv_dx * (p[i, j+1] - p[i, j-1])
        end
    end

    # Aplicamos condiciones de contorno finales al campo de velocidades resultante
    set_bnd!(u, 1, params)
    set_bnd!(v, 2, params)
end
