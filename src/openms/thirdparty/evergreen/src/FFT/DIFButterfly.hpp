#ifndef _DIFBUTTERFLY_HPP
#define _DIFBUTTERFLY_HPP

#include "Twiddles.hpp"

template<unsigned long N>
class DIFButterfly {
public:
  // Can improve speed, but makes compilation very resource intensive
  // for larger N:

  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    // Butterfly, then multiply twiddles into second half of list:

    cpx twiddle = cpx{1.0, 0.0};
    for (unsigned long i=0; i<N/2; ++i) {
      cpx temp = data[(i + N/2)];
      data[(i + N/2)] = data[i] - temp;
      data[(i + N/2)] *= twiddle;

      data[i] += temp;

      // Compute next twiddle factor:
      Twiddles<N/2>::advance(twiddle);
    }
      
    DIFButterfly<N/2>::apply(data);
    DIFButterfly<N/2>::apply(data+N/2);
  }
};

template<>
class DIFButterfly<0ul> {
public:
  inline static void apply(cpx* __restrict  const) {
    // do nothing
  }
};

template<>
class DIFButterfly<1ul> {
public:
  inline static void apply(cpx* __restrict  const) {
    // do nothing
  }
};

// Same as DITButterfly<2>:
template<>
class DIFButterfly<2ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict  const data) {
    data[1] = data[0] - data[1];
    data[0] = data[0] + data[0] - data[1];
  }
};

// Same as DITButterfly<4>, but shuffled so [1] and [2] are swapped:
template<>
class DIFButterfly<4ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    // Note: index written even when index is zero for
    // clarity; however, all constant multiplications will be
    // computed at compile time:

    cpx t = data[2];
    data[2] = data[0] - t;
    data[0] += t;
    t = data[3];
    data[3] = cpx{data[1].i - t.i, t.r - data[1].r};
    data[1] += t;
    t = data[1];
    data[1] = data[0] - t;
    data[0] += t;
    t = data[3];
    data[3] = data[2] - t;
    data[2] += t;
  }
};

// No need to manually swap, because <4> is called manually, as it would be above:
template<>
class DIFButterfly<8ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {

    const double sqrt2Over2 = Twiddles<4>::sin();
    
    cpx temp = data[4];
    data[4] = data[0] - temp;
    data[0] += temp;
    
    cpx twiddle{sqrt2Over2, -sqrt2Over2};
    temp = data[5];
    data[5] = data[1] - temp;
    data[5] *= twiddle;
    data[1] += temp;
    
    twiddle = cpx{0.0, -1.0};
    temp = data[6];
    data[6] = data[2] - temp;
    data[6] *= twiddle;
    data[2] += temp;
    
    twiddle = cpx{-sqrt2Over2, -sqrt2Over2};
    temp = data[7];
    data[7] = data[3] - temp;
    data[7] *= twiddle;
    data[3] += temp;

    DIFButterfly<4ul>::apply(data);
    DIFButterfly<4ul>::apply(data+4);
  }
};

// No need to manually swap, because <8> is called manually, as it would be above:
template<>
class DIFButterfly<16ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {

    const double sqrt2Over2 = Twiddles<4>::sin();
    const double sinPiOver8 = Twiddles<8>::sin();
    const double cosPiOver8 = Twiddles<8>::cos();
    
    cpx temp = data[8];
    data[8] = data[0] - temp;
    data[0] += temp;
    
    cpx twiddle{cosPiOver8, -sinPiOver8};
    temp = data[9];
    data[9] = data[1] - temp;
    data[9] *= twiddle;
    data[1] += temp;
    
    twiddle = cpx{sqrt2Over2, -sqrt2Over2};
    temp = data[10];
    data[10] = data[2] - temp;
    data[10] *= twiddle;
    data[2] += temp;
    
    twiddle = cpx{sinPiOver8, -cosPiOver8};
    temp = data[11];
    data[11] = data[3] - temp;
    data[11] *= twiddle;
    data[3] += temp;

    twiddle = cpx{0.0, -1.0};
    temp = data[12];
    data[12] = data[4] - temp;
    data[12] *= twiddle;
    data[4] += temp;

    twiddle = cpx{-sinPiOver8, -cosPiOver8};
    temp = data[13];
    data[13] = data[5] - temp;
    data[13] *= twiddle;
    data[5] += temp;

    twiddle = cpx{-sqrt2Over2, -sqrt2Over2};
    temp = data[14];
    data[14] = data[6] - temp;
    data[14] *= twiddle;
    data[6] += temp;

    twiddle = cpx{-cosPiOver8, -sinPiOver8};
    temp = data[15];
    data[15] = data[7] - temp;
    data[15] *= twiddle;
    data[7] += temp;

    DIFButterfly<8ul>::apply(data);
    DIFButterfly<8ul>::apply(data+8);
  }
};

#endif
