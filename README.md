# [StableFluids.jl](https://github.com/luinux/StableFluids.jl)

A pure Julia implementation of Jos Stam's [Stable Fluids](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf) algorithm (1999/2003), widely used in computer graphics and games to simulate smoke, fire, and gaseous phenomena in real time. It is lightweight, unconditionally stable, and serves as the reference implementation and playground for embedded ports like [ESP32-stableFluids](https://github.com/Luiulator/ESP32-stableFluids).



## A few Highlights

This package solves the incompressible 2D Navier-Stokes equations for fluid flow:

**Continuity Equation (Incompressibility):**
$$ \nabla \cdot \mathbf{u} = 0 $$

**Momentum Equation:**
$$ \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = \nu \nabla^2 \mathbf{u} - \frac{1}{\rho} \nabla p + \mathbf{f} $$

**Dye / Smoke Transport:**
$$ \frac{\partial \rho}{\partial t} + (\mathbf{u} \cdot \nabla) \rho = \kappa \nabla^2 \rho + S $$

It operates via operator splitting, advancing each physical step independently within each time frame:
1. **External Forces:** Adding body forces (e.g. gravitational buoyancy, IMU tilt, or user mouse impulses) to the velocity field.
2. **Viscous Diffusion:** Implicit solver for fluid viscosity using Gauss-Seidel relaxation.
3. **Pressure Projection:** Solving the Poisson equation ($\nabla^2 p = \nabla \cdot \mathbf{u}^*$) via Helmholtz-Hodge decomposition to ensure a divergence-free velocity field.
4. **Semi-Lagrangian Advection:** Unconditionally stable back-tracing of streamlines with bilinear interpolation to transport velocities and scalar densities.
5. **Density Diffusion & Dissipation:** Transporting and gradually decaying smoke concentration to keep high contrast and prevent saturation.

Lastly, Stable Fluids is, as per its name, unconditionally stable. That means you won't get infinite velocities under any time step.

*Be advised!!!* As a tradeoff for unconditional stability, the semi-Lagrangian advection method introduces numerical diffusion (artificial smoothing). You should expect smooth, highly aesthetic fluid-like motion with swirling vortices, but do not use it for high-precision aerodynamic or DNS CFD simulations.

---

## Interactive Simulation

An interactive real-time simulation is available in [`JuliaImplementation/examples/rt_interactive.jl`](JuliaImplementation/examples/rt_interactive.jl):

- **Mouse Click & Drag:** Spawns smoke density and imparts directional velocity momentum along the drag vector.
- **W / A / S / D:** Controls 2D gravity / acceleration vector dynamically (simulating an IMU accelerometer tilt).
- **SPACE:** Injects a dense puff of smoke at the center.
- **Adaptive Brush:** Brush radius and spawn size scale dynamically to any grid resolution ($50 \times 50$, $200 \times 200$, $300 \times 300$, etc.).
- **Live FPS Counter:** Real-time performance monitoring displayed in the title bar.

To launch the interactive demo:
```bash
julia --project=JuliaImplementation JuliaImplementation/examples/rt_interactive.jl
```

---

## Running Tests

Run the test suite covering grid initialization, boundary conditions, implicit diffusion, advection, projection, and continuous dissipation:

```bash
julia --project=JuliaImplementation JuliaImplementation/test/runtests.jl
```

---

## Future Work

A list of optimizations and features I want to add to this:

- Multi-threaded **Red-Black Gauss-Seidel** solver with `@simd` vectorization for high-resolution real-time grids ($500 \times 500+$).
- GPU acceleration using `KernelAbstractions.jl` / `CUDA.jl`.
- **Vorticity Confinement** (Fedkiw et al.) to counteract numerical dissipation and preserve fine turbulent swirls.
- Extension to 3D grid volumes with volumetric rendering.