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

  // Actual inputs:
  fftw_complex*x = new fftw_complex[n];
  fftw_complex*y = new fftw_complex[n];
  // Initialize input data:
  for (int i=0; i<n; ++i) {
    x[i][0] = i;
    x[i][1] = i;
  }

  // Buffers (must be hard-coded in place to reuse plan)
  fftw_complex*in = (fftw_complex*)fftw_malloc(n*sizeof(fftw_complex));
  fftw_complex*out = (fftw_complex*)fftw_malloc(n*sizeof(fftw_complex));

  std::cout << n << " ";

  // Cold start:
  Clock c;
  const int shape[] = {n};
  memcpy(in, x, n*sizeof(fftw_complex));
  auto plan = fftw_plan_dft(sizeof(shape)/sizeof(int), &shape[0], in, out, FFTW_FORWARD, FFTW_MEASURE);
  fftw_execute(plan);
  memcpy(y, out, n*sizeof(fftw_complex));
  std::cout << c.tock() << " ";

  // Re-initialize the data:
  for (int i=0; i<n; ++i) {
    x[i][0] = i;
    x[i][1] = -i;
  }

  // Warm start:
  c.tick();
  memcpy(in, x, n*sizeof(fftw_complex));
  fftw_execute(plan);
  memcpy(y, out, n*sizeof(fftw_complex));
  std::cout << c.tock() << std::endl;
  fftw_destroy_plan(plan);

  fftw_free(in);
  fftw_free(out);

  delete[] x;
  delete[] y;

  return 0;
}
