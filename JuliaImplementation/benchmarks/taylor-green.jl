using Pkg
Pkg.activate(dirname(@__DIR__))

using StableFluids
using CairoMakie
using BenchmarkTools

function run_benchmark(nx, ny, iters)
    # Constantes para el constructor

    L = 1.0f0
    U_0 = 1.0f0
    k = pi / L

    nu = 0.01f0
    dx, dy = L / nx, L / ny

    # esto es un arreglillo para que el número de courant sea constante
    # C = 0.016f0 aunque vayamos escalando la malla. A ver si resuelve el problema del error

    dt = 0.1f0 * dx / U_0

    t_final = 1.0f0 / 30.0f0
    n_steps = round(Int, t_final / dt)

    error_L2 = zeros(n_steps)
    times = zeros(n_steps)
    print_freq = 20

    #en T-G no nos interesa la densidad, por eso diff=0
    grid = FluidGrid2D(nx, ny; dx=dx, dt=dt, visc=nu, diff=0.0f0)

    # Establecemos las condiciones iniciales del vórtice de T-G
    for j in 2:ny+1
        for i in 2:nx+1

            x = (i - 1.5f0) * dx
            y = (j - 1.5f0) * dy

            grid.u[i, j] = U_0 * sin(k * x) * cos(k * y)
            grid.v[i, j] = -U_0 * cos(k * x) * sin(k * y)

        end
    end

    set_bnd!(grid.u, 1, grid.params)
    set_bnd!(grid.v, 2, grid.params)

    # Funciones que guardan las soluciones analíticas del campo de velocidades
    u_exact(x, y, t) = U_0 * sin(k * x) * cos(k * y) * exp(-2 * nu * k^2 * t)
    v_exact(x, y, t) = -U_0 * cos(k * x) * sin(k * y) * exp(-2 * nu * k^2 * t)

    # Bucle principal:
    # Avanzamos pasos de la simulación y calculamos el error

    for step in 1:n_steps

        step!(grid; iter=iters)
        t = step * dt
        times[step] = t

        diffsq = 0.0f0

        for j in 2:ny+1
            @inbounds for i in 2:nx+1

                x = (i - 1.5f0) * dx
                y = (j - 1.5f0) * dy

                diffsq += (grid.u[i, j] - u_exact(x, y, t))^2 + (grid.v[i, j] - v_exact(x, y, t))^2

            end
        end

        diffsq = diffsq / (nx * ny)

        error_L2[step] = sqrt(diffsq)

        #if step % print_freq == 0
        #    println("Paso ", step, " | Error L2: ", error_L2[step])
        #end
    end

    return times, error_L2
end

mesh_refinement = 3
iter_refinement = 5

fig = Figure(size=(800, 500))

#= ax = Axis(fig[1, 1],
    title="Variación del Error con el número de iteraciones de Gauss-Siedel",
    xlabel="Tiempo simulado (s)",
    ylabel="Error L2 en velocidad",
    xgridvisible=true,
    ygridvisible=true
) =#

#= for mult in 1:mesh_refinement
    nx, ny = mult * 32, mult * 32
    #gs_iters = mult * 20   #Número de iteraciones incrementales

    time, error = run_benchmark(nx, ny, 50)

    lines!(ax, time, error, linewidth=2.5, label="$(nx)x$(ny) 50 iters GS")
end

for mult in 1:mesh_refinement
    nx, ny = mult * 32, mult * 32
    #gs_iters = mult * 20   #Número de iteraciones incrementales

    time, error = run_benchmark(nx, ny, 100)

    lines!(ax, time, error, linewidth=2.5, label="$(nx)x$(ny) 100 iters GS")
end
 =#

#=
for mult in 1:iter_refinement
    nx = 32
    ny = nx

    gs_iters = mult * 20

    time, error = run_benchmark(nx, ny, gs_iters)

    lines!(ax, time, error, linewidth=2.5, label="$(nx)x$(ny) $(gs_iters) iteraciones")

end
=#

# =============================================================================
# Cálculo de la Vorticidad Numérica (Rotacional 2D)
# =============================================================================

function compute_vorticity(grid::FluidGrid2D{T}) where {T}
    nx, ny = grid.params.nx, grid.params.ny
    dx = grid.params.dx
    inv_2dx = T(0.5) / dx

    omega = zeros(T, nx, ny)

    for j in 2:ny+1
        @inbounds for i in 2:nx+1
            # Diferencias centrales usando las celdas contiguas
            dv_dx = (grid.v[i+1, j] - grid.v[i-1, j]) * inv_2dx
            du_dy = (grid.u[i, j+1] - grid.u[i, j-1]) * inv_2dx

            # Guardamos en matriz indexada de 1 a nx, 1 a ny
            omega[i-1, j-1] = dv_dx - du_dy
        end
    end

    return omega
end

# =============================================================================
# Visualización y Comparación de Vorticidad (Numérica vs Analítica)
# =============================================================================

function compare_vorticity(nx=64, ny=64, iters=50; t_target=1.0f0 / 30.0f0)
    L = 1.0f0
    U_0 = 1.0f0
    k = Float32(pi) / L
    nu = 0.01f0
    dx, dy = L / nx, L / ny
    dt = 0.1f0 * dx / U_0

    n_steps = round(Int, t_target / dt)
    t = n_steps * dt

    grid = FluidGrid2D(nx, ny; dx=dx, dt=dt, visc=nu, diff=0.0f0)

    # Condiciones iniciales de Taylor-Green
    for j in 2:ny+1
        for i in 2:nx+1
            x = (i - 1.5f0) * dx
            y = (j - 1.5f0) * dy
            grid.u[i, j] = U_0 * sin(k * x) * cos(k * y)
            grid.v[i, j] = -U_0 * cos(k * x) * sin(k * y)
        end
    end
    set_bnd!(grid.u, 1, grid.params)
    set_bnd!(grid.v, 2, grid.params)

    # Avanzar simulación
    for step in 1:n_steps
        step!(grid; iter=iters)
    end

    # Coordenadas en los centros de celda
    xs = [(i - 0.5f0) * dx for i in 1:nx]
    ys = [(j - 0.5f0) * dy for j in 1:ny]

    # Vorticidad analítica
    ω_exact(x, y, t_val) = 2.0f0 * U_0 * k * sin(k * x) * sin(k * y) * exp(-2.0f0 * nu * k^2 * t_val)

    # Vorticidad numérica
    ω_num = compute_vorticity(grid)
    ω_ana = [ω_exact(x, y, t) for x in xs, y in ys]
    ω_diff = ω_num .- ω_ana

    # =========================================================================
    # Visualización con Makie (contourf)
    # =========================================================================
    fig_vort = Figure(size=(2*380, 2*380))
    cmap = :roma

    val_max = max(maximum(abs, ω_num), maximum(abs, ω_ana))
    vort_levels = range(-val_max, val_max, length=50)

    # 1. Vorticidad Numérica
    ax1 = Axis(fig_vort[1, 1],
        title="Vorticidad Numérica (t = $(round(t, digits=4))s)",
        xlabel="x", ylabel="y",
        aspect=DataAspect()
    )
    cf1 = contourf!(ax1, xs, ys, ω_num, colormap=cmap, levels=vort_levels)
    Colorbar(fig_vort[1, 2], cf1, label="ω_num")

    # 2. Vorticidad Analítica
    ax2 = Axis(fig_vort[1, 3],
        title="Vorticidad Analítica (t = $(round(t, digits=4))s)",
        xlabel="x", ylabel="y",
        aspect=DataAspect()
    )
    cf2 = contourf!(ax2, xs, ys, ω_ana, colormap=cmap, levels=vort_levels)
    Colorbar(fig_vort[1, 4], cf2, label="ω_exact")

    # 3. Diferencia (Error celda por celda)
    err_max = max(maximum(abs, ω_diff), 1.0f-6)
    err_levels = range(-err_max, err_max, length=25)
    ax3 = Axis(fig_vort[2, 1],
        title="Error (Numérica - Analítica)",
        xlabel="x", ylabel="y",
        aspect=DataAspect()
    )
    cf3 = contourf!(ax3, xs, ys, ω_diff, colormap=cmap, levels=err_levels)
    contour!(ax3, xs, ys, ω_diff, levels=30, linewidth=0.85 , color=:white)
    Colorbar(fig_vort[2, 2], cf3, label="Δω")

    # Guardar imagen y retornar figura
    save("taylor_green_vorticity.png", fig_vort)
    println("Visualización guardada en taylor_green_vorticity.png")

    return fig_vort
end

# Ejecutar la comparación
fig = compare_vorticity(64, 64, 50)
display(fig)