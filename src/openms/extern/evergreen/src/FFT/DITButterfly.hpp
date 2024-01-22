#ifndef _DITBUTTERFLY_HPP
#define _DITBUTTERFLY_HPP

#include "Twiddles.hpp"

template<unsigned long N>
class DITButterfly {
public:
  // Can improve speed, but makes compilation very resource intensive
  // for larger N:

  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    DITButterfly<N/2>::apply(data);
    DITButterfly<N/2>::apply(data+N/2);

    cpx twiddle = cpx{1.0, 0.0};
    for (unsigned long i=0; i<N/2; ++i) {
      cpx temp = data[(i + N/2)]*twiddle;
      data[(i + N/2)] = data[i] - temp;
      data[i] += temp;

      // Compute next twiddle factor:
      Twiddles<N/2>::advance(twiddle);
    }
  }
};

template<>
class DITButterfly<0ul> {
public:
  inline static void apply(cpx* __restrict  const) {
    // do nothing
  }
};

template<>
class DITButterfly<1ul> {
public:
  inline static void apply(cpx* __restrict  const) {
    // do nothing
  }
};

template<>
class DITButterfly<2ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict  const data) {
    data[1] = data[0] - data[1];
    data[0] = data[0] + data[0] - data[1];
    // x = x + x - y 
    // equivalent to 
    // x <<= 1
    // x += y
    //    data[0] += data[0];
    //    data[0] -= data[1];
  }
};

template<>
class DITButterfly<4ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    // Note: index written even when index is zero for
    // clarity; however, all constant multiplications will be
    // computed at compile time:
    
    cpx t = data[1];
    data[1] = data[0] - t;
    data[0] += t;
    t = data[3];
    data[3] = cpx{data[2].i - t.i, t.r - data[2].r};
    data[2] += t;
    t = data[2];
    data[2] = data[0] - t;
    data[0] += t;
    t = data[3];
    data[3] = data[1] - t;
    data[1] += t;
  }
};

template<>
class DITButterfly<8ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    DITButterfly<4ul>::apply(data);
    DITButterfly<4ul>::apply(data+4);
    
    const double sqrt2Over2 = Twiddles<4>::sin();
    
    cpx temp = data[4];
    data[4] = data[0] - temp;
    data[0] += temp;
    
    cpx twiddle{sqrt2Over2, -sqrt2Over2};
    temp = data[(1 + 4)]*twiddle;
    data[(1 + 4)] = data[1] - temp;
    data[1] += temp;
    
    twiddle = cpx{0.0, -1.0};
    temp = data[(2 + 4)]*twiddle;
    data[(2 + 4)] = data[2] - temp;
    data[2] += temp;
    
    twiddle = cpx{-sqrt2Over2, -sqrt2Over2};
    temp = data[(3 + 4)]*twiddle;
    data[(3 + 4)] = data[3] - temp;
    data[3] += temp;
  }
};

template<>
class DITButterfly<16ul> {
public:
  //  __attribute__((always_inline))
  inline static void apply(cpx* __restrict const data) {
    DITButterfly<8ul>::apply(data);
    DITButterfly<8ul>::apply(data+8);

    const double sqrt2Over2 = Twiddles<4>::sin();
    const double sinPiOver8 = Twiddles<8>::sin();
    const double cosPiOver8 = Twiddles<8>::cos();

    cpx temp = data[8];
    data[8] = data[0] - temp;
    data[0] += temp;
    
    cpx twiddle{cosPiOver8, -sinPiOver8};
    temp = data[(1 + 8)]*twiddle;
    data[(1 + 8)] = data[1] - temp;
    data[1] += temp;
    
    twiddle = cpx{sqrt2Over2, -sqrt2Over2};
    temp = data[(2 + 8)]*twiddle;
    data[(2 + 8)] = data[2] - temp;
    data[2] += temp;
    
    twiddle = cpx{sinPiOver8, -cosPiOver8};
    temp = data[(3 + 8)]*twiddle;
    data[(3 + 8)] = data[3] - temp;
    data[3] += temp;

    twiddle = cpx{0, -1.0};
    temp = data[(4 + 8)]*twiddle;
    data[(4 + 8)] = data[4] - temp;
    data[4] += temp;

    twiddle = cpx{-sinPiOver8, -cosPiOver8};
    temp = data[(5 + 8)]*twiddle;
    data[(5 + 8)] = data[5] - temp;
    data[5] += temp;

    twiddle = cpx{-sqrt2Over2, -sqrt2Over2};
    temp = data[(6 + 8)]*twiddle;
    data[(6 + 8)] = data[6] - temp;
    data[6] += temp;

    twiddle = cpx{-cosPiOver8, -sinPiOver8};
    temp = data[(7 + 8)]*twiddle;
    data[(7 + 8)] = data[7] - temp;
    data[7] += temp;
  }
};

#endif
