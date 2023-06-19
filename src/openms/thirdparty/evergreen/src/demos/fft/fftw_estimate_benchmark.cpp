#include <iostream>
#include <cstring>
#include "fftw3.h"
#include "../../Utility/Clock.hpp"

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: fftw_benchmark <LOG_N>" << std::endl;
    return 1;
  }
  int log_n = atoi(argv[1]);
  int n = 1<<log_n;

  // To avoid allocating buffers, use existing memory (which in C++,
  // will usually be allocated with new); the alternative would be to
  // use fftw_malloc and then copy memory.
  fftw_complex*x = new fftw_complex[n];
  fftw_complex*y = new fftw_complex[n];
  // Initialize input data:
  for (int i=0; i<n; ++i) {
    x[i][0] = i;
    x[i][1] = i;
  }

  std::cout << n << " ";

  // Cold start with FFTW_ESTIMATE (good use-case for FFTs of unknown
  // size, no buffers needed):
  Clock c;
  const int shape[] = {n};
  auto plan = fftw_plan_dft(sizeof(shape)/sizeof(int), &shape[0], x, y, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  std::cout << c.tock() << std::endl;

  delete[] x;
  delete[] y;

  return 0;
}
