#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <cassert>
#include <iostream>

using namespace OpenMS;

void testSmoothData_basic()
{
  ModifiedSincSmoother smoother;

  std::vector<double> input{1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> output;

  smoother.smoothData(input, output, 3, 0.5);

  assert(output.size() == input.size());

  for (size_t i = 0; i < output.size(); ++i)
  {
    assert(std::abs(output[i] - 1.0) < 1e-10);
  }
}

void testSmooth_inplace()
{
  ModifiedSincSmoother smoother;
  std::vector<double> data{1.0, 2.0, 3.0, 4.0, 5.0};

  smoother.smooth(data, 2, 1); // degree = 2, half_kernel_size = 1

  for (size_t i = 0; i < data.size(); ++i)
  {
    std::cout << "smoothed[" << i << "] = " << data[i] << std::endl;
  }
}

void testInvalidWindowSize()
{
  ModifiedSincSmoother smoother;
  std::vector<double> input{1.0, 2.0, 3.0};
  std::vector<double> output;

  try
  {
    smoother.smoothData(input, output, 2, 0.5); // even window size
    assert(false); // should not reach here
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Caught expected exception: " << e.what() << std::endl;
  }
}

void testInvalidCutoff()
{
  ModifiedSincSmoother smoother;
  std::vector<double> input{1.0, 2.0, 3.0};
  std::vector<double> output;

  try
  {
    smoother.smoothData(input, output, 3, -1.0); // invalid cutoff
    assert(false);
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Caught expected exception: " << e.what() << std::endl;
  }
  try
  {
    smoother.smoothData(input, output, 3, 1.5); // invalid cutoff
    assert(false);
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Caught expected exception: " << e.what() << std::endl;
  }
  try
  {
    smoother.smoothData(input, output, 3, 0.0); // invalid cutoff
    assert(false);
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Caught expected exception: " << e.what() << std::endl;
  }
}

int main()
{
  std::cout << "Running ModifiedSincSmoother tests..." << std::endl;
  testSmoothData_basic();
  testSmooth_inplace();
  testInvalidWindowSize();
  testInvalidCutoff();
  std::cout << "All tests completed." << std::endl;
  return 0;
}
