#include "../include/stable_fluids/StableFluids.hpp"
#include "../include/stable_fluids/forces.hpp"
#include <iostream>

int main() {
  // 1. Instanciamos una malla de 32x32
  stable_fluids::FluidGrid2D<32, 32> grid;

  // 2. Probamos las fuerzas y el paso de simulación
  stable_fluids::add_buoyancy(grid, 0.0f, -9.81f);
  stable_fluids::step(grid, 20);

  std::cout << "✅ ¡Todo compila y ejecuta correctamente!" << std::endl;
  return 0;
}
