#pragma once

#include "types.hpp"

/* 
Por suerte para mí en C++ se pueden usar ternarios como en Julia
condicion ? valor_si_true : valor_si_false
*/

namespace stable_fluids{

    template <size_t NX, size_t NY>
    void set_bnd(int b, float (&x)[NX+2][NY+2]) {

        for (size_t j =1; j <= NY; j++) {
            x[0][j]     = (b==1)    ? -x[1][j]      : x[1][j];
            x[NX+1][j]  = (b==1)    ? -x[NX][j]     : x[NX][j];
        }

        for (size_t i=1; i <= NX; i++) {
            x[i][0]     = (b==2)    ? -x[i][1]      : x[i][1];
            x[i][NY+1]  = (b==2)    ? -x[i][NY]     : x[i][NY];
        }

        x[0][0]           = 0.5f * (x[1][0]          + x[0][1]);
        x[0][NY + 1]      = 0.5f * (x[1][NY + 1]     + x[0][NY]);
        x[NX + 1][0]      = 0.5f * (x[NX][0]         + x[NX + 1][1]);
        x[NX + 1][NY + 1] = 0.5f * (x[NX][NY + 1]    + x[NX + 1][NY]);
    }
}