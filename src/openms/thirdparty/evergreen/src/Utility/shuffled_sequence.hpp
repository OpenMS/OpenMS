#ifndef _SHUFFLED_SEQUENCE_HPP
#define _SHUFFLED_SEQUENCE_HPP

inline std::vector<unsigned long> shuffled_sequence(const unsigned long N) {
  std::vector<unsigned long> result(N);
  for (unsigned long i=0; i<N; ++i)
    result[i] = i;

  // Swap with random index:
  for (unsigned long i=0; i<N; ++i)
    std::swap(result[i], result[rand() % N]);

  return result;
}

#endif
