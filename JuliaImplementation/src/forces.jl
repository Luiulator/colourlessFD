"""
    add_buoyancy!(grid::FluidGrid2D{T}; gx=zero(T), gy=T(-9.81)) where {T}

Aplica el vector de aceleración/gravedad (gx, gy) de la IMU proporcional
a la concentración de densidad (humo/materia). Permite inclinar y agitar
el dispositivo en 2D haciendo que el gas reaccione a la orientación real.
"""
function add_buoyancy!(grid::FluidGrid2D{T}; gx=zero(T), gy=T(-9.81)) where {T}
    dt = grid.params.dt
    nx, ny = grid.params.nx, grid.params.ny
    gx_T = T(gx)
    gy_T = T(gy)

    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            d = grid.density[i, j]
            grid.u[i, j] += d * gx_T * dt
            grid.v[i, j] += d * gy_T * dt
        end
    end
end

"""
    dissipate_density!(grid::FluidGrid2D{T}; decay=T(0.995)) where {T}

Aplica una tasa de disipación / decaimiento a la densidad del humo en cada paso
para evitar que el lienzo se sature completamente con el tiempo.
"""
function dissipate_density!(grid::FluidGrid2D{T}; decay=T(0.995)) where {T}
    decay_T = T(decay)
    nx, ny = grid.params.nx, grid.params.ny

    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            grid.density[i, j] *= decay_T
        end
    end
end
