#ifndef _CLOCK_HPP
#define _CLOCK_HPP

#include <iostream>
#include <iomanip>
#include <chrono>

class Clock {
protected:
  std::chrono::steady_clock::time_point startTime;
public:
  Clock() {
    tick();
  }
  void tick() {
    startTime = std::chrono::steady_clock::now();
  }
  float tock() {
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();;
    return float(std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()) / 1e6f;
  }
  // Print elapsed time with newline
  void ptock() {
    float elapsed = tock();
    std::cout << "Took " << elapsed << " seconds" << std::endl;
  }
};

#endif
