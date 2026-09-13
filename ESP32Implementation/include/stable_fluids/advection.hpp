#pragma once

#include "boundary.hpp"
#include "types.hpp"

// A partir de C++17 está la función std::clamp en la cabecera <algorithm>
#include <algorithm>

namespace stable_fluids {

template <size_t NX, size_t NY>
void advect(float (&d)[NX + 2][NY + 2], const float (&d0)[NX + 2][NY + 2],
            float (&u)[NX + 2][NY + 2], const float (&v)[NX + 2][NY + 2], int b,
            FluidParams<NX, NY> &params) {

  float dt0 = params.dt / params.dx;

  for (size_t i = 1; i <= NX; i++) {
    for (size_t j = 1; j <= NY; j++) {

      float x = i - dt0 * u[i][j];
      float y = j - dt0 * v[i][j];

      // std::clamp necesita que los 3 argumentos sean del mismo tipo, por lo
      // que hay que castear
      x = std::clamp(x, 0.5f, static_cast<float>(NX) + 0.5f);
      y = std::clamp(y, 0.5f, static_cast<float>(NY) + 0.5f);

      // De paso, hay que castear también los índices de las celdas
      size_t i0 = static_cast<size_t>(x);
      size_t i1 = i0 + 1;
      size_t j0 = static_cast<size_t>(y);
      size_t j1 = j0 + 1;

      // Y ahora devolverlos a float para hallar los pesos xd
      float s1 = x - static_cast<float>(i0);
      float s0 = 1.0f - s1;
      float t1 = y - static_cast<float>(j0);
      float t0 = 1.0f - t1;

      // Ahora sí, bilineal interpolation
      d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
                s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
    }
  }

  set_bnd<NX, NY>(b, d);
}
} // namespace stable_fluids