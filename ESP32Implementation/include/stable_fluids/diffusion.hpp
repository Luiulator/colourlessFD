#pragma once

#include <cstring>
#include "types.hpp"
#include "boundary.hpp"

namespace stable_fluids {

    template <size_t NX, size_t NY>
    void lin_solve(int b, float (&x)[NX+2][NY+2], const float (&x0)[NX+2][NY+2],
                   float a, float c, int iter = 20) {

        float inv_c = 1.0f / c;

        for (int it = 0; it < iter; ++it) {
            for (size_t i = 1; i <= NX; ++i) {
                for (size_t j = 1; j <= NY; ++j) {
                    float neighbours_sum = x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1];
                    x[i][j] = (x0[i][j] + a * neighbours_sum) * inv_c;
                }
            }
            set_bnd<NX, NY>(b, x);
        }
    }

    template <size_t NX, size_t NY>
    void diffuse(float (&x)[NX+2][NY+2], const float (&x0)[NX+2][NY+2], float diff_rate,
                 int b, const FluidParams<NX, NY>& params, int iter = 20) {

        if (diff_rate <= 0.0f) {

            /* 
            C++ no me deja usar directamente x=x0. Para copiar todos los números es o hacer un bucle o usar std::memcpy
            No está usando memoria dinámica porque sizeof(x) calcula el tamaño en bytes que ocupa la matriz en tiempo de compilación
             */

            std::memcpy(x, x0, sizeof(x));
            return;
        }

        float a = params.dt * diff_rate / (params.dx * params.dx);
        float c = 1.0f + 4.0f * a;

        lin_solve<NX, NY>(b, x, x0, a, c, iter);
    }

} // namespace stable_fluids
