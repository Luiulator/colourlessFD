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

    t_final = 1.0f0/30.0f0
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

ax = Axis(fig[1, 1],
    title="Variación del Error con el número de iteraciones de Gauss-Siedel",
    xlabel="Tiempo simulado (s)",
    ylabel="Error L2 en velocidad",
    xgridvisible=true,
    ygridvisible=true
)

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

@benchmark run_benchmark(32, 32, 30) samples=500 seconds=180

#axislegend(ax)

#display(fig)