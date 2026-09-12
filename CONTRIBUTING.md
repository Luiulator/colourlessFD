# Guía de Contribución a StableFluids.jl

¡Gracias por tu interés en contribuir a **StableFluids.jl**! Toda ayuda, desde corrección de errores y documentación hasta optimizaciones de rendimiento y nuevos ejemplos, es bienvenida.

---

## 🛠️ Entorno de Desarrollo Local

1. **Clonar el repositorio:**
   ```bash
   git clone https://github.com/tu-usuario/StableFluids.jl.git
   cd StableFluids.jl
   ```

2. **Abrir Julia en el entorno del paquete:**
   ```bash
   julia --project=.
   ```

3. **Ejecutar la suite de pruebas unitarias:**
   ```bash
   julia --project=. test/runtests.jl
   ```

---

## 📋 Directrices de Código

- **Estabilidad de tipos**: Asegúrate de que las funciones numéricas admitan el parámetro genérico `T<:AbstractFloat` y no generen conversiones implícitas no deseadas.
- **Cero asignaciones dinámicas**: El ciclo `step!(grid)` está diseñado para no realizar `malloc` ni crear matrices intermedias en caliente. Usa mutaciones con `@inbounds` y difusión in-place.
- **Documentación**: Mantén los docstrings en formato Markdown estándar de Julia explicando los fundamentos matemáticos y las condiciones de contorno.

---

## 🚀 Flujo de Pull Requests

1. Crea una rama descriptiva para tu cambio: `git checkout -b feature/mi-mejora` o `git checkout -b fix/mi-bug`.
2. Asegúrate de que todos los tests pasen (`julia --project=. test/runtests.jl`).
3. Añade nuevos tests si introduces nueva funcionalidad.
4. Abre un Pull Request detallando los cambios y las decisiones de diseño adoptadas.
