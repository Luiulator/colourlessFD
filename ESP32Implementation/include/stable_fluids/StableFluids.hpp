#pragma once

#include "advection.hpp"
#include "diffusion.hpp"
#include "projection.hpp"
#include "types.hpp"
#include <cstring>

namespace stable_fluids {
template <size_t NX, size_t NY>
void step(FluidGrid2D<NX, NY>& grid, int iter = 20) {
  FluidParams<NX, NY>& params = grid.params;

  if (params.visc > 0.0f) {
    std::memcpy(grid.u_prev, grid.u, sizeof(grid.u));
    std::memcpy(grid.v_prev, grid.v, sizeof(grid.v));
    diffuse(grid.u, grid.u_prev, params.visc, 1, params, iter);
    diffuse(grid.v, grid.v_prev, params.visc, 2, params, iter);
  }

  project(grid.u, grid.v, grid.p, grid.div, params, iter = iter);

  std::memcpy(grid.u_prev, grid.u, sizeof(grid.u));
  std::memcpy(grid.v_prev, grid.v, sizeof(grid.v));

  advect(grid.u, grid.u_prev, grid.u_prev, grid.v_prev, 1, params);
  advect(grid.v, grid.v_prev, grid.u_prev, grid.v_prev, 2, params);

  project(grid.u, grid.v, grid.p, grid.div, params, iter);

  if (params.diff > 0.0f) {
    std::memcpy(grid.density_prev, grid.density, sizeof(grid.density));
    diffuse(grid.density, grid.density_prev, params.diff, 0, params,
            iter = iter);
  }

  std::memcpy(grid.density_prev, grid.density, sizeof(grid.density));
  advect(grid.density, grid.density_prev, grid.u, grid.v, 0, params);
}

} // namespace stable_fluids