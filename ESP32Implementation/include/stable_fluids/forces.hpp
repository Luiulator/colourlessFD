#pragma once

#include "types.hpp"

namespace stable_fluids {
template <size_t NX, size_t NY>
void add_buoyancy(FluidGrid2D<NX, NY> &grid, float gx = 0.0f,
                  float gy = -9.81f) {
  float dt = grid.params.dt;

  for (size_t i = 1; i <= NX; i++) {
    for (size_t j = 1; j <= NY; j++) {

      float d = grid.density[i][j];
      grid.u[i][j] += d * gx * dt;
      grid.v[i][j] += d * gy * dt;
    }
  }
}

}