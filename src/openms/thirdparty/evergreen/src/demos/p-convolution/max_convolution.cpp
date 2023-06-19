#include <iostream>
#include <fstream>
#include "../../Convolution/p_convolve.hpp"
#include "../../Utility/Clock.hpp"

Clock c;

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "usage: max_conv <filename with n and x and y>" << std::endl;
    exit(1);
  }

  std::cout.precision(100);

  std::ifstream fin(argv[1]);

  unsigned long n;
  fin >> n;
  Tensor<double> x({n});
  for (unsigned long i=0; i<n; ++i)
    fin >> x[i];
  Tensor<double> y({n});
  for (unsigned long i=0; i<n; ++i)
    fin >> y[i];

  Clock c;
  auto z = naive_max_convolve(x,y);
  c.ptock();
  std::cout << z.flat() << std::endl;

  c.tick();
  auto z2 = numeric_p_convolve(x,y,std::numeric_limits<double>::infinity());
  c.ptock();
  std::cout << z2.flat() << std::endl;
}
