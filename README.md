# StableFluids.jl

A pure Julia implementation of Jos Stam's [Stable Fluids](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf) algorithm (1999/2003), widely used in computer graphics and games to simulate smoke, fire, and gaseous phenomena in real time. It is lightweight, unconditionally stable, and serves as the reference implementation and playground for embedded ports like [ESP32-stableFluids](https://github.com/Luiulator/ESP32-stableFluids).


## A few Highlights

This package solves the incompressible 2D Navier-Stokes equations for fluid flow:

$$ \nabla \cdot \mathbf{u} = 0 $$

$$ \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = \nu \nabla^2 \mathbf{u} - \frac{1}{\rho} \nabla p + \mathbf{f} $$

and also an equation for transporting the dye

$$ \frac{\partial \rho}{\partial t} + (\mathbf{u} \cdot \nabla) \rho = \kappa \nabla^2 \rho + S $$

It does so by breaking down each term of the sum and solving them separately, then summing up. Lastly, it applies a correction pressure so that we enforce that the fluid remains incompressible.

Lastly, Stable Fluids is, as per its name, uncondi1tionally stable. That means you won't get infinite velocities under any time step.

*Be advised!!!* As a tradeoff for unconditional stability, the semi-Lagrangian advection method introduces numerical diffusion (artificial smoothing). You should expect smooth, highly aesthetic fluid-like motion with swirling vortices, but do not use it for high-precision aerodynamic or CFD simulations.

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
- **Vorticity Confinement** to counteract numerical dissipation and preserve fine turbulent swirls.
- Extension to 3D grid volumes with volumetric rendering.