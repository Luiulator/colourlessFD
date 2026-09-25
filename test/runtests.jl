using Test
using StableFluids

@testset "StableFluids.jl Suite" begin

    @testset "Inicialización de Malla y Parámetros" begin
        # 1. Grid Float32 por defecto
        nx, ny = 32, 64
        grid = FluidGrid2D(nx, ny)
        @test grid isa FluidGrid2D{Float32}
        @test size(grid.u) == (nx + 2, ny + 2)
        @test size(grid.v) == (nx + 2, ny + 2)
        @test size(grid.density) == (nx + 2, ny + 2)
        @test size(grid.p) == (nx + 2, ny + 2)
        @test size(grid.div) == (nx + 2, ny + 2)
        @test grid.params.nx == nx
        @test grid.params.ny == ny
        @test grid.params.dx == 1.0f0
        @test grid.params.dt == 0.016f0
        @test all(grid.u .== 0.0f0)
        @test all(grid.density .== 0.0f0)

        # 2. Grid Float64 explícito
        grid64 = FluidGrid2D(Float64, 16, 16; dx=0.5, dt=0.01, visc=0.001, diff=0.002)
        @test grid64 isa FluidGrid2D{Float64}
        @test grid64.params.dx ≈ 0.5
        @test grid64.params.visc ≈ 0.001
    end

    @testset "Condiciones de Contorno (set_bnd!)" begin
        nx, ny = 10, 10
        grid = FluidGrid2D(nx, ny)
        params = grid.params

        # Matriz de prueba
        mat = zeros(Float32, nx + 2, ny + 2)
        mat[2:nx+1, 2:ny+1] .= 5.0f0

        # Escalar (b = 0, Neumann)
        set_bnd!(mat, 0, params)
        @test mat[1, 5] == mat[2, 5] == 5.0f0
        @test mat[nx+2, 5] == mat[nx+1, 5] == 5.0f0
        @test mat[5, 1] == mat[5, 2] == 5.0f0
        @test mat[5, ny+2] == mat[5, ny+1] == 5.0f0

        # Velocidad Horizontal u (b = 1, Dirichlet en X, Neumann en Y)
        set_bnd!(mat, 1, params)
        @test mat[1, 5] == -mat[2, 5]  # No penetración en pared vertical
        @test mat[nx+2, 5] == -mat[nx+1, 5]
        @test mat[5, 1] == mat[5, 2]   # Deslizamiento en suelo

        # Velocidad Vertical v (b = 2, Neumann en X, Dirichlet en Y)
        set_bnd!(mat, 2, params)
        @test mat[1, 5] == mat[2, 5]   # Deslizamiento en paredes verticales
        @test mat[5, 1] == -mat[5, 2]  # No penetración en suelo
        @test mat[5, ny+2] == -mat[5, ny+1] # No penetración en techo
    end

    @testset "Difusión Implícita (diffuse!)" begin
        nx, ny = 16, 16
        grid = FluidGrid2D(nx, ny; dt=0.1f0, diff=0.5f0)
        grid.density[8:9, 8:9] .= 10.0f0
        initial_mass = sum(grid.density[2:nx+1, 2:ny+1])

        grid.density_prev .= grid.density
        diffuse!(grid.density, grid.density_prev, grid.params.diff, 0, grid.params; iter=40)

        # La densidad debe haberse suavizado
        @test grid.density[8, 8] < 10.0f0
        @test grid.density[7, 8] > 0.0f0
        # Conservación aproximada de masa en celdas interiores con contorno Neumann
        diffused_mass = sum(grid.density[2:nx+1, 2:ny+1])
        @test isapprox(initial_mass, diffused_mass, rtol=0.05)
    end

    @testset "Advección Semi-Lagrangiana (advect!)" begin
        nx, ny = 20, 20
        grid = FluidGrid2D(nx, ny; dt=1.0f0, dx=1.0f0)
        
        # Flujo constante hacia la derecha u = 1, v = 0
        grid.u .= 1.0f0
        grid.v .= 0.0f0
        
        # Pulso inicial en x = 5
        grid.density[5, 10] = 5.0f0
        grid.density_prev .= grid.density

        advect!(grid.density, grid.density_prev, grid.u, grid.v, 0, grid.params)

        # Con u = 1 y dt = 1, el pulso debió moverse hacia la derecha (i = 6)
        @test grid.density[6, 10] > grid.density[5, 10]
    end

    @testset "Proyección de Helmholtz-Hodge (∇ · u = 0)" begin
        nx, ny = 16, 16
        grid = FluidGrid2D(nx, ny)
        
        # Introducir velocidades con divergencia conocida
        for j in 2:ny+1, i in 2:nx+1
            grid.u[i, j] = sin(Float32(i) / 2.0f0)
            grid.v[i, j] = cos(Float32(j) / 2.0f0)
        end
        set_bnd!(grid.u, 1, grid.params)
        set_bnd!(grid.v, 2, grid.params)

        calc_max_div(u, v) = maximum([
            abs(-0.5f0 * grid.params.dx * ((u[i+1, j] - u[i-1, j]) + (v[i, j+1] - v[i, j-1])))
            for i in 2:nx+1, j in 2:ny+1
        ])

        div_before = calc_max_div(grid.u, grid.v)

        # Proyectar
        project!(grid.u, grid.v, grid.p, grid.div, grid.params; iter=40)

        div_after = calc_max_div(grid.u, grid.v)

        # La divergencia máxima debe reducirse significativamente
        @test div_after < div_before
        @test all(isfinite.(grid.u))
        @test all(isfinite.(grid.v))
    end

    @testset "Fuerza de Gravedad / IMU (add_buoyancy!)" begin
        nx, ny = 16, 16
        grid = FluidGrid2D(nx, ny; dt=0.1f0)
        grid.density[8, 8] = 2.0f0
        
        add_buoyancy!(grid; gx=10.0f0, gy=-20.0f0)
        
        @test grid.u[8, 8] ≈ 2.0f0 * 10.0f0 * 0.1f0
        @test grid.v[8, 8] ≈ 2.0f0 * (-20.0f0) * 0.1f0
        @test grid.u[2, 2] == 0.0f0 # Sin densidad no hay aceleración
    end

    @testset "Disipación / Decaimiento de Humo (dissipate_density!)" begin
        nx, ny = 16, 16
        grid = FluidGrid2D(nx, ny)
        grid.density[8, 8] = 2.0f0

        dissipate_density!(grid; decay=0.9f0)
        @test grid.density[8, 8] ≈ 1.8f0

        # Celdas vacías siguen en cero
        @test grid.density[2, 2] == 0.0f0
    end

    @testset "Ciclo Completo y Estabilidad Numérica (step!)" begin
        nx, ny = 32, 32
        grid = FluidGrid2D(nx, ny; dt=0.016f0, visc=0.01f0, diff=0.01f0)
        
        # Inyección central
        grid.density[15:18, 15:18] .= 5.0f0
        grid.u[15:18, 15:18] .= 2.0f0
        grid.v[15:18, 15:18] .= 3.0f0

        # Simular 30 pasos de tiempo
        for step in 1:30
            add_buoyancy!(grid; gx=0.0f0, gy=-9.81f0)
            step!(grid)
        end

        # Verificar que no hay NaNs ni Infs (incondicionalmente estable)
        @test all(isfinite.(grid.u))
        @test all(isfinite.(grid.v))
        @test all(isfinite.(grid.density))
        @test all(isfinite.(grid.p))
        @test any(grid.density .> 0.0f0)
    end
end
