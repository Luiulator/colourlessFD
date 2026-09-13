#pragma once

#include "boundary.hpp"
#include "diffusion.hpp"
#include "types.hpp"

namespace stable_fluids {

template <size_t NX, size_t NY>
void project(float (&u)[NX + 2][NY + 2], float (&v)[NX + 2][NY + 2],
             float (&p)[NX + 2][NY + 2], float (&div)[NX + 2][NY + 2],
             FluidParams<NX, NY> &params, int iters = 20) {

  float dx = params.dx;

  for (size_t i = 1; i <= NX; i++) {
    for (size_t j = 1; j <= NY; j++) {

      div[i][j] = -0.5f * dx *
                  ((u[i + 1][j] - v[i - 1][j]) + (v[i][j + 1] - v[i][j - 1]));

      p[i][j] = 0.0f;
    }
  }

  set_bnd<NX, NY>(div, 0);
  set_bnd<NX, NY>(p, 0);

  /*
  void lin_solve(int b, float (&x)[NX+2][NY+2], const float (&x0)[NX+2][NY+2],
                   float a, float c, int iter = 20)

    b:= condicioón de contorno, en este caso b=0 es para escalares
    x:= matriz de salida que hay que resolver
    x0:= término independiente
    a:=peso de  los vecinos
    c:= número de vecinos (4 en 2D)
  */
  lin_solve<NX, NY>(0, p, div, 1.0f, 4.0f, 0, params, iters = iters);

  float inv_dx = 0.5f / dx;

  for (size_t i = 1; i <= NX; i++) {
    for (size_t j = 1; i <= NY; j++) {
      u[i][j] -= inv_dx * (p[i + 1][j] - p[i - 1][j]);
      v[i][j] -= inv_dx * (p[i][j + 1] - p[i][j - 1]);
    }
  }

  set_bnd<NX, NY>(u, 1);
  set_bnd<NX, NY>(v, 2);
};

} // namespace stable_fluids