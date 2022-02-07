# lmfit2 - Levenberg-Marquardt minimization

A modification of the original *lmfit* library to simplify its use and allow some parallelization.

Just add `lmfit.c` and `lmfit.h` to your projects (don't forget the license). Or add this as a submodule or with `FetchContent` and link to `lmfit2`. Curve fitting is also included.

This version allows parallelization with OpenMP. If your `evaluate()` function is very expensive,
this may accelerate computation by calling it in parallel when computing the Jacobian. If the function is simple
there will be no gain, or even a performance penalty due to threading overhead.

To enable multithreading enable `LMFIT_OPENMP` in CMake configuration. Or if you don't use Cmake, enable OpenMP
in your compiler and define the macro `LMFIT_OPENMP` when compiling `lmfit.c`.

`main.cpp` is just a dumb usage example derived from the original project.
