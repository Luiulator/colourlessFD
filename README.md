# StableFluids.jl

<div align="center">

[![CI](https://github.com/luinux/StableFluids.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/luinux/StableFluids.jl/actions/workflows/CI.yml)
[![Julia 1.9+](https://img.shields.io/badge/julia-v1.9%2B-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)

**Simulación incondicionalmente estable de fluidos incompresibles 2D en Julia pura.**  
*Basado en el algoritmo Stable Fluids de Jos Stam (SIGGRAPH 1999 & GDC 2003).*

</div>

---

## 🌟 Características Principales

- **Incondicionalmente Estable**: Emplea advección Semi-Lagrangiana, permitiendo pasos de tiempo ($\Delta t$) arbitrariamente grandes sin divergencia ni inestabilidades numéricas.
- **Campos Solenoidales Rigurosos ($\nabla \cdot \mathbf{u} = 0$)**: Implementa la mejora de Stam (2003) proyectando el campo de velocidad antes y después de la advección para preservar la incompresibilidad física.
- **Optimizado para Sistemas Embebidos y FPUs**: Tipado nativo en `Float32` por defecto (ideal para microcontroladores y hardware embebido como la FPU de ESP32) con soporte genérico para `Float64`.
- **Cero Asignaciones Dinámicas de Memoria**: La memoria se preasigna en la inicialización (`FluidGrid2D`), garantizando 0 *heap allocations* durante el bucle principal de simulación.
- **Condiciones de Contorno Robustas**: Celdas fantasma (*ghost cells*) con condiciones de *Free-Slip* / No-Penetración para paredes y esquinas.
- **Cero Dependencias Gráficas Obligatorias**: El núcleo del motor es puro Julia; las librerías gráficas (Makie) solo son requeridas por ejemplos específicos.

---

## 📐 Fundamentos Teóricos

El motor resuelve numéricamente las **Ecuaciones de Navier-Stokes para fluidos incompresibles**:

$$\frac{\partial \mathbf{u}}{\partial t} = -(\mathbf{u} \cdot \nabla)\mathbf{u} + \nu \nabla^2 \mathbf{u} - \frac{1}{\rho}\nabla p + \mathbf{f}$$

$$\nabla \cdot \mathbf{u} = 0$$

### El Ciclo de Simulación (`step!`)

```text
┌────────────────────────────────────────────────────────────────────────┐
│                        Ciclo de Navier-Stokes                          │
│                                                                        │
│   1. Difusión de velocidad (ν) ───► Solver Gauss-Seidel implícito      │
│   2. Proyección 1               ───► Descomposición de Helmholtz-Hodge │
│   3. Autoadvección              ───► Trazado Semi-Lagrangiano          │
│   4. Proyección 2 (Final)       ───► Garantiza ∇ · u = 0               │
│   5. Difusión de densidad       ───► Dispersión de humo/tinte          │
│   6. Advección de densidad      ───► Transporte escalar                │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 🚀 Instalación y Uso Rápido

### Instalación vía Julia Package Manager

```julia
using Pkg
Pkg.add(url="https://github.com/luinux/StableFluids.jl.git")
```

### Ejemplo Mínimo

```julia
using StableFluids

# 1. Crear una malla de 64x64 celdas
grid = FluidGrid2D(64, 64; dt=0.03f0, visc=0.01f0, diff=0.005f0)

# 2. Inyectar densidad (humo) y velocidad en el centro
grid.density[30:34, 30:34] .= 5.0f0
grid.v[30:34, 30:34] .= 10.0f0

# 3. Bucle de física
for frame in 1:200
    add_buoyancy!(grid; gx=0.0f0, gy=-9.81f0) # Aplicar gravedad / aceleración IMU
    step!(grid)                               # Resolver Navier-Stokes
    
    # grid.density, grid.u y grid.v están listos para renderizar o transmitir
end
```

---

## 🧪 Ejemplos (`examples/`)

### 🎮 Simulación Interactiva WASD + IMU (GLMakie)
Visualizador gráfico estilo display OLED con control de inclinación y gravedad IMU mediante las teclas `W`, `A`, `S`, `D` y `Espacio`:
```bash
julia --project=. examples/test_gravity.jl
```

---

## 🧩 Arquitectura de la Malla (`FluidGrid2D`)

Para evitar problemas de frontera sin coste computacional extra, el dominio físico de $N_x \times N_y$ se rodea de una capa de **celdas fantasma** (*ghost cells*):

```text
    j = ny+2  ┌─────────────────────────────────┐  (Pared Superior)
              │ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ │
              │ ░ ┌─────────────────────────┐ ░ │
              │ ░ │                         │ ░ │
              │ ░ │    Celdas Interiores    │ ░ │
              │ ░ │      (Simulación)       │ ░ │
              │ ░ │                         │ ░ │
              │ ░ └─────────────────────────┘ ░ │
              │ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ ░ │
    j = 1     └─────────────────────────────────┘  (Suelo)
              i = 1                          i = nx+2
```

- **Escalares (Densidad / Presión)**: Condición de Neumann ($\partial x / \partial n = 0$).
- **Velocidades ($u, v$)**: Componente normal invertida en paredes para forzar velocidad nula en la frontera física exacta (No-Penetración) y componente tangencial copiada (Free-Slip).

---

## 🔬 Ejecución de Pruebas Unitarias

La suite de pruebas automatizadas valida la inicialización, condiciones de contorno, difusión, advección, reducción de divergencia en proyección y estabilidad multi-paso:

```bash
julia --project=. test/runtests.jl
```

---

## 📚 Referencias

1. **Stam, Jos (1999)**. *"Stable Fluids"*. Proceedings of the 26th Annual Conference on Computer Graphics and Interactive Techniques (SIGGRAPH '99), pp. 121–128.
2. **Stam, Jos (2003)**. *"Real-Time Fluid Dynamics for Games"*. Proceedings of the Game Developer Conference (GDC '03).
3. **Chorin, A. J. (1968)**. *"Numerical solution of the Navier-Stokes equations"*. Mathematics of Computation, 22(104), pp. 745–762.

---

## 📄 Licencia

Este proyecto está distribuido bajo la licencia **MIT**. Consulta el archivo [LICENSE](LICENSE) para más detalles.
