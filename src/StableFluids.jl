module StableFluids

export FluidParams, FluidGrid2D
export step!, add_buoyancy!, dissipate_density!, diffuse!, advect!, project!, lin_solve!, set_bnd!

include("types.jl")
include("forces.jl")
include("boundary.jl")
include("diffusion.jl")
include("advection.jl")
include("projection.jl")

"""
    step!(grid::FluidGrid2D)

Ejecuta un ciclo completo de Navier-Stokes incompresible (Stable Fluids):
1. Difusión de velocidad (si visc > 0)
2. Proyección de velocidad (∇ · u = 0)
3. Autoadvección de velocidad (u y v se transportan a sí mismas)
4. Proyección de velocidad final (∇ · u = 0)
5. Difusión de densidad (si diff > 0)
6. Advección de la densidad con el nuevo campo de velocidades

Aclaro un poco, en el paper original del 99 el orden era fuerzas -> advección -> difusión -> proyección.
En una evolución del método que publicó el mismo autor en 2003 Stam, Jos "Real Time Fluid Dynamics for Games"
establece la mejora que representa proyectar antes y después de la advección dado que para aplicar la advección
semi-lagrangiana, uno de los requisitos (que no tuvo en cuenta originalmente) es que el campo al que se la aplicas
sea solenoidal.

Como tanto la adición de fuerzas externas como la difusión introducen una divergencia al campo de velocidades, se
requiere ese paso de proyección previo a la advección semi-lagrangiana.
"""
function step!(grid::FluidGrid2D{T}; iter::Int=20) where {T}
    params = grid.params

    # 1. Difusión de la velocidad (viscosidad)
    if params.visc > zero(T)
        grid.u_prev .= grid.u
        grid.v_prev .= grid.v
        diffuse!(grid.u, grid.u_prev, params.visc, 1, params; iter=iter)
        diffuse!(grid.v, grid.v_prev, params.visc, 2, params; iter=iter)
    end

    # 2. Proyección antes de advectar (asegura campo solenoidal)
    project!(grid.u, grid.v, grid.p, grid.div, params; iter=iter)

    # 3. Autoadvección de la velocidad
    grid.u_prev .= grid.u
    grid.v_prev .= grid.v
    advect!(grid.u, grid.u_prev, grid.u_prev, grid.v_prev, 1, params)
    advect!(grid.v, grid.v_prev, grid.u_prev, grid.v_prev, 2, params)

    # 4. Proyección final (garantiza incompresibilidad ∇ · u = 0)
    project!(grid.u, grid.v, grid.p, grid.div, params; iter=iter)

    # 5. Difusión del escalar / densidad
    if params.diff > zero(T)
        grid.density_prev .= grid.density
        diffuse!(grid.density, grid.density_prev, params.diff, 0, params;  iter=iter)
    end

    # 6. Advección del escalar / densidad
    grid.density_prev .= grid.density
    advect!(grid.density, grid.density_prev, grid.u, grid.v, 0, params)
end

end # module StableFluids