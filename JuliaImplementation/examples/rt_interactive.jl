using Pkg
Pkg.activate(dirname(@__DIR__))

using StableFluids
using GLMakie

function run_simulation()
    nx, ny = 300, 300
    grid = FluidGrid2D(nx, ny; dt=0.016f0, visc=1.0f0, diff=0.5f0)

    # Radios adaptativos proporcionales a la resolución de la malla
    min_dim = min(nx, ny)
    smoke_radius = max(1, round(Int, 0.03f0 * min_dim))   # ~3% del dominio (ej. 6 en 200x200, 2 en 50x50)
    brush_radius = max(1, round(Int, 0.025f0 * min_dim))  # ~2.5% del dominio (ej. 5 en 200x200, 1 en 50x50)

    # Inyección de humo inicial
    center_x, center_y = nx / 2, ny / 2
    for j in 2:ny+1
        for i in 2:nx+1
            dist2 = (i - 1 - center_x)^2 + (j - 1 - center_y)^2
            if dist2 <= smoke_radius^2
                grid.density[i, j] = 1.0f0
            end
        end
    end

    
    display_buffer = Observable(copy(grid.density[2:nx+1, 2:ny+1]))
    title_obs = Observable("FPS: --")
    
    fig = Figure(backgroundcolor=:black)
    ax = Axis(fig[1, 1], aspect=DataAspect(),
        title=title_obs,
        titlecolor=:white, titlesize=14, backgroundcolor=:black)
    hidedecorations!(ax)
    empty!(interactions(ax))

    heatmap!(ax, display_buffer, colormap=:grays, colorrange=(0.0f0, 1.0f0), interpolate=false)
    screen = display(fig)

    #=
    println("==========================================================")
    println("Simulación interactiva:")
    println(" • Click + Arrastre : Inyectar densidad e inducir velocidad/advección")
    println(" • A / D : Inclinar a la izquierda / derecha (IMU X)")
    println(" • W / S : Inclinar hacia arriba / abajo (IMU Y)")
    println(" • ESPACIO : Inyectar nueva nube de gas")
    println("==========================================================")
    =#

    prev_mouse = nothing
    last_time = time()
    frame_count = 0

    while isopen(screen)
        # 1. Lectura del teclado con WASD
        gx = 0.0f0
        gy = -15.9f0  # Gravedad natural hacia abajo en reposo



        if ispressed(fig, Mouse.left)
            mp = mouseposition(ax.scene)
            mx, my = Float32(mp[1]), Float32(mp[2])
            cx, cy = round(Int, mx), round(Int, my)

            # Calcular el desplazamiento del ratón para inducir advección
            if prev_mouse !== nothing
                dx_drag = mx - prev_mouse[1]
                dy_drag = my - prev_mouse[2]
            else
                dx_drag = 0.0f0
                dy_drag = 0.0f0
            end
            prev_mouse = (mx, my)

            force_scale = 15.0f0  # Escala de fuerza aplicada por el arrastre

            for j in max(2, cy - brush_radius + 1):min(ny + 1, cy + brush_radius + 1)
                for i in max(2, cx - brush_radius + 1):min(nx + 1, cx + brush_radius + 1)
                    if (i - 1 - cx)^2 + (j - 1 - cy)^2 <= brush_radius^2
                        grid.density[i, j] = 1.0f0
                        grid.u[i, j] += dx_drag * force_scale
                        grid.v[i, j] += dy_drag * force_scale
                    end
                end
            end
        else
            prev_mouse = nothing
        end

        if ispressed(fig, Makie.Keyboard.a)
            gx -= 35.0f0
        end
        if ispressed(fig, Makie.Keyboard.d)
            gx += 35.0f0
        end
        if ispressed(fig, Makie.Keyboard.w)
            gy = 25.0f0   # Poner boca arriba / acelerar hacia el techo
        end
        if ispressed(fig, Makie.Keyboard.s)
            gy = -35.0f0  # Fuerte aceleración hacia el suelo
        end

        # Inyectar con la barra espaciadora
        if ispressed(fig, Makie.Keyboard.space)
            for j in max(2, round(Int, center_y - smoke_radius + 1)):min(ny + 1, round(Int, center_y + smoke_radius + 1))
                for i in max(2, round(Int, center_x - smoke_radius + 1)):min(nx + 1, round(Int, center_x + smoke_radius + 1))
                    if ((i - 1 - center_x)^2 + (j - 1 - center_y)^2) <= smoke_radius^2
                        grid.density[i, j] = 4.0f0
                    end
                end
            end
        end

        # 2. Aplicar la aceleración vectorial de la IMU al gas
        add_buoyancy!(grid; gx=gx, gy=gy)

        # 3. Paso de simulación completo
        step!(grid)

        # 4. Actualizar visualización
        display_buffer[] = grid.density[2:nx+1, 2:ny+1]

        # 5. Contador de FPS
        frame_count += 1
        t_now = time()
        elapsed = t_now - last_time
        if elapsed >= 0.25
            fps = frame_count / elapsed
            title_obs[] = string("FPS: ", round(fps, digits=1))
            frame_count = 0
            last_time = t_now
        end

        sleep(0.001)
    end
end

# Ejecutar automáticamente al invocar el script
run_simulation()
