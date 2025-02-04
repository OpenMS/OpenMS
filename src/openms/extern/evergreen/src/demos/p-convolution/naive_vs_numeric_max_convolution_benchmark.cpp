#include <iostream>
#include "../../Convolution/p_convolve.hpp"
#include "../../Utility/Clock.hpp"

Clock c;

void init_data(Tensor<double> & x, Tensor<double> & y) {
  unsigned long k;
  for (k=0; k<x.flat_size(); ++k)
    x[k] = exp( - (k - 128.0)*(k - 128.0) / (100.0*100.0) );
  for (k=0; k<y.flat_size(); ++k)
    y[k] = x[k] + exp( - (k - 700.0)*(k - 700.0) / (10.0*10.0) );
}

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: convolution_benchmark <LOG_N>" << std::endl;
    return 1;
  }

  const unsigned int log_n = atoi(argv[1]);
  const unsigned long n = 1ul<<log_n;
    
  Tensor<double> x({n});
  Tensor<double> y({n});

  init_data(x,y);

  x.flat() /= sum( x.flat() );
  y.flat() /= sum( y.flat() );

  Clock c;

  c.tick();
  auto z = naive_max_convolve(x,y);
  std::cout << n << " " << c.tock() << " ";

  c.tick();
  auto z2 = numeric_p_convolve(x,y,std::numeric_limits<double>::infinity());
  std::cout << c.tock() << std::endl;
}
