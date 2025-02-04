#include <iostream>
#include <cstring>
#include "../../Utility/Clock.hpp"
#include "../../Convolution/p_convolve.hpp"

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: conv_benchmark <LOG_N>" << std::endl;
    return 1;
  }
  int log_n = atoi(argv[1]);
  unsigned long n = 1ul<<log_n;

  // Actual inputs:
  Tensor<cpx> x({n});
  for (unsigned long i=0; i<n; ++i)
    x[i] = cpx{double(i),double(i)};

  Tensor<cpx> y({n});
  for (unsigned long i=0; i<n; ++i)
    y[i] = cpx{-double(i),-double(i)};

  std::cout << n << " ";

  Clock c;
  Tensor<cpx> z = naive_convolve(x, y);
  std::cout << c.tock() << std::endl;

  return 0;
}
