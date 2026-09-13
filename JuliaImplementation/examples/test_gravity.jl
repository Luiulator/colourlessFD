using Pkg
Pkg.activate(dirname(@__DIR__))

using StableFluids
using GLMakie

function run_simulation()
    # 1. Malla 32x128
    nx, ny = 32, 128
    grid = FluidGrid2D(nx, ny; dt=0.016f0, visc=1.0f0, diff=0.5f0)

    # Inyección de humo inicial
    center_x, center_y = 16, 64
    radius = 6
    for j in 2:ny+1
        for i in 2:nx+1
            dist2 = (i - 1 - center_x)^2 + (j - 1 - center_y)^2
            if dist2 <= radius^2
                grid.density[i, j] = 4.0f0
            end
        end
    end

    # 2. Buffer Observable para Makie (32x128)
    display_buffer = Observable(copy(grid.density[2:nx+1, 2:ny+1]))
    title_obs = Observable("IMU: W A S D | [Espacio]: Inyectar")

    # 3. Ventana estilo pantalla OLED
    fig = Figure(size=(280, 850), backgroundcolor=:black)
    ax = Axis(fig[1, 1], aspect=DataAspect(), 
              title=title_obs, 
              titlecolor=:white, titlesize=14, backgroundcolor=:black)
    hidedecorations!(ax)

    heatmap!(ax, display_buffer, colormap=:grays, colorrange=(0.0f0, 1.0f0), interpolate=false)
    screen = display(fig)

    println("==========================================================")
    println("Simulación interactiva WASD iniciada:")
    println(" • A / D : Inclinar a la izquierda / derecha (IMU X)")
    println(" • W / S : Inclinar hacia arriba / abajo (IMU Y)")
    println(" • ESPACIO : Inyectar nueva nube de gas")
    println("==========================================================")

    while isopen(screen)
        # 1. Lectura del teclado con WASD
        gx = 0.0f0
        gy = -15.0f0  # Gravedad natural hacia abajo en reposo

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
            for j in max(2, center_y - radius):min(ny + 1, center_y + radius)
                for i in max(2, center_x - radius):min(nx + 1, center_x + radius)
                    if ((i - 1 - center_x)^2 + (j - 1 - center_y)^2) <= radius^2
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

        sleep(0.016)
    end
end

# Ejecutar automáticamente al invocar el script
run_simulation()
