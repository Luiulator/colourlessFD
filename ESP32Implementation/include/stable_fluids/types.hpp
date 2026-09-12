#pragma once

#include <cstddef>

namespace stable_fluids {

    /* 
    Aquí he hecho un template. Si este código lo quisiera correr en otro microcontrolador más potente o menos potente, tendría
    que ir cambiando el tamaño de la malla. Para aseguirar que solo tendré que hacer el cambio en un sitio, convierto el struct en un
    template. La idea es que FluidParams deja de ser un tipo de dato por sí solo, y ahora acepta NX y NY como parámetros para construir
    un tipo de dato a partir de ellos. P.Ej:

    FluidParams<32, 32> el compilador creará un tipo con nx=32, ny=32
    FluidParams<64, 64> hará lo mismo con nx=64, ny=64


    Esto es útil porque sabemos en tiempo de compilación qué espacio total ocupa la malla, es decir, el compilador podrá optimizar el código
    con loop unrolling.
     */


    /* 
    Otro recordatorio sobre el loop unrolling, porque tuve que refrescar. Cuando escribes un bucle for:

    for (int i = 0; i < 4; i++) {
        array[i] = 0.0f;
    }

    La CPU del ESP32 no solo ejecuta la asignación array[i] = 0.0f. En cada una de las 4 vueltas, la CPU está obligada a hacer "trabajo burocrático extra":

    - array[i] = 0.0f
    - i++
    - ¿i < 4?
    - Volver a saltar a la instrucción del principio del bucle

    El loop unrolling es básicamente que el compilador desenrolla el bucle eliminando todos los saltos y comprobaciones que puede.
    Lo que hace es traducir directamente ese código a:

    
    array[0] = 0.0f;
    array[1] = 0.0f;
    array[2] = 0.0f;
    array[3] = 0.0f;

    Te ahorra incrementos de contador, comparaciones, etc.
    */

    template <size_t NX, size_t NY>
    struct FluidParams {   
        static constexpr size_t nx = NX;
        static constexpr size_t ny = NY;
        float dx;
        float dt;
        float visc;
        float diff;

        FluidParams(float dx_ = 1.0f, float dt_ = 0.016f,
                float visc_ = 0.0f, float diff_ = 0.0f)
            :dx(dx_), dt(dt_), visc(visc_), diff(diff_) {}

    };

    template <size_t NX, size_t NY>
    struct FluidGrid2D {
        FluidParams<NX, NY> params;


        // La sentencia = {} le asigna directamente ceros en todas las posiciones a las matrices

        float u[NX+2][NY+2] = {};
        float u_prev[NX+2][NY+2] = {};
        float v[NX+2][NY+2] = {}; 
        float v_prev[NX+2][NY+2] = {};

        float density[NX+2][NY+2] = {};
        float density_prev[NX+2][NY+2] = {};

        float p[NX+2][NY+2] = {};
        float div[NX+2][NY+2] = {};

        FluidGrid2D() = default;


        /*
        & es para pasar por referencia, va escrito luego del tipo
        
        Recordatorio:
        
        FluidGrid2D(FluidParams p_params) haría una copia innecesaria en memoria
        FluidGrid2D(FluidParams& p_params) no haría la copia, pero el constructor podría modificar accidentalmente p_params
        
        Lo mejor es 
        
        FluidGrid2D(const FluidParams& p_params)
        */

        FluidGrid2D(const FluidParams<NX, NY>& p_params)
        : params(p_params) {}
    };
}