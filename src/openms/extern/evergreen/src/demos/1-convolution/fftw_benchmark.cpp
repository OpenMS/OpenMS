#include <iostream>
#include <cstring>
#include "fftw3.h"
#include "../../Utility/Clock.hpp"

fftw_complex* convolve(fftw_complex*x, fftw_complex*y, int n) {
  // Buffers:
  fftw_complex *input = (fftw_complex*)fftw_malloc(2*n*sizeof(fftw_complex));
  fftw_complex *output = (fftw_complex*)fftw_malloc(2*n*sizeof(fftw_complex));
  fftw_complex *temp = (fftw_complex*)fftw_malloc(2*n*sizeof(fftw_complex));

  const int shape[] = {2*n};
  auto plan = fftw_plan_dft(sizeof(shape)/sizeof(int), &shape[0], input, output, FFTW_FORWARD, FFTW_ESTIMATE);
  
  // Zero pad x:
  for (int i=0; i<n; ++i) {
    input[i][0] = x[i][0];
    input[i][1] = x[i][1];
  }
  for (int i=n; i<2*n; ++i) {
    input[i][0] = 0;
    input[i][1] = 0;
  }
  
  // FFT zero padded x:
  fftw_execute(plan);

  // Copy FFT of zero padded x to temp:
  for (int i=0; i<2*n; ++i) {
    temp[i][0] = output[i][0];
    temp[i][1] = output[i][1];
  }
  
  // Zero pad y:
  for (int i=0; i<n; ++i) {
    input[i][0] = y[i][0];
    input[i][1] = y[i][1];
  }
  for (int i=n; i<2*n; ++i) {
    input[i][0] = 0;
    input[i][1] = 0;
  }
  // FFT zero padded y:
  fftw_execute(plan);

  // Multiply FFT results:
  for (int i=0; i<2*n; ++i) {
    const double r1 = output[i][0];
    const double i1 = output[i][1];

    const double r2 = temp[i][0];
    const double i2 = temp[i][1];

    input[i][0] = r1*r2 - i1*i2;

    // Conjugate inline:
    input[i][1] = -(i1*r2 + r1*i2);
  }

  // input contains conjugated FFT of result

  // Conjugate input and output to reuse plan (input conjugation is
  // already performed above):
  fftw_execute(plan);
  
  fftw_destroy_plan(plan);

  // Conjugate output and divide by 2*n:
  fftw_complex *z = new fftw_complex[2*n-1];
  double one_over_two_n = 1.0 / (2*n);
  for (int i=0; i<2*n-1; ++i) {
    // Multiplication is faster than division:
    z[i][0] = output[i][0] * one_over_two_n;
    // Conjugation inline:
    z[i][1] = output[i][1] * -one_over_two_n;
  }

  fftw_free(input);
  fftw_free(output);
  fftw_free(temp);
  
  return z;  
}

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: fftw_conv_benchmark <LOG_N>" << std::endl;
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

  // Initialize input data:
  for (int i=0; i<n; ++i) {
    y[i][0] = -i;
    y[i][1] = -i;
  }
  
  std::cout << n << " ";

  Clock c;
  fftw_complex*z = convolve(x, y, n);
  std::cout << c.tock() << std::endl;

  delete[] x;
  delete[] y;
  delete[] z;

  return 0;
}
