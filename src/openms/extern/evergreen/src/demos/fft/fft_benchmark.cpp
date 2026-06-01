#include <iostream>
#include "../../FFT/FFT.hpp"
#include "../../Utility/Clock.hpp"

Clock c;

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: fft_benchmark <LOG_N>" << std::endl;
    return 1;
  }
  int log_n = atoi(argv[1]);
  unsigned long n = 1ul<<log_n;

  Tensor<cpx> x({n});
  for (unsigned long i=0; i<n; ++i)
    x[i] = cpx{double(i), double(i)};

  std::cout << n << " ";

  Clock c;
  // In-place FFT:
  // true, true arguments say to apply shuffling and to undo
  // transpositions. If the application was complex convolution, both
  // of these could be false to get the convolution result faster.
  apply_fft<DIF, true, true>(x);
  std::cout << c.tock() << std::endl;
}
