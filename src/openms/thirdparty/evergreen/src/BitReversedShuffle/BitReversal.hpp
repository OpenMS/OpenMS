#ifndef _BITREVERSAL_HPP
#define _BITREVERSAL_HPP

template <unsigned char LOG_N>
class BitReversal {
protected:
  static const unsigned char reversed_byte_table[256];

  inline static int fast_log2(unsigned int i) {
    // Note: This is left in for reference, but it can be performed
    // more generally (i.e., LOG_N>=25) using the clz opcode and then
    // subtracting.
    static_assert(LOG_N < 25, "Fast logarithm by float casting only works for 24 bits or less; result for larger numbers of bits can be built from the total number of bits and either __builtin_clz or __builtin_clzl.");
    float f = i;
    unsigned int exponent = ((*(unsigned int *)&f >> 23) - 0x7f);
    return exponent;
  }

public:
  inline static unsigned int reverse_int_logical(unsigned int x) {
    // swap odd and even bits
    x = ((x >> 1) & 0x55555555) | ((x & 0x55555555) << 1);
    // swap consecutive pairs
    x = ((x >> 2) & 0x33333333) | ((x & 0x33333333) << 2);
    // swap nibbles ... 
    x = ((x >> 4) & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    // swap bytes
    x = ((x >> 8) & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    // swap 2-byte long pairs
    x = ( x >> 16             ) | ( x               << 16);
    return x;
  }
  
  inline static unsigned short reverse_short_byte_table(unsigned short x){
    unsigned char inByte0 = (x & 0xFF);
    unsigned char inByte1 = (x & 0xFF00) >> 8;
    return (reversed_byte_table[inByte0] << 8) | reversed_byte_table[inByte1];
  }

  inline static unsigned int reverse_int_byte_table(unsigned int x){
    unsigned char inByte0 = (x & 0xFF);
    unsigned char inByte1 = (x & 0xFF00) >> 8;
    unsigned char inByte2 = (x & 0xFF0000) >> 16;
    unsigned char inByte3 = (x & 0xFF000000) >> 24;
    return (reversed_byte_table[inByte0] << 24) | (reversed_byte_table[inByte1] << 16) | (reversed_byte_table[inByte2] << 8) | reversed_byte_table[inByte3];
  }

  inline static unsigned long reverse_long_byte_table(unsigned long x){
    unsigned char inByte0 = (x & 0xFF);
    unsigned char inByte1 = (x & 0xFF00) >> 8;
    unsigned char inByte2 = (x & 0xFF0000) >> 16;
    unsigned char inByte3 = (x & 0xFF000000) >> 24;
    unsigned char inByte4 = (x & 0xFF00000000ul) >> 32;
    unsigned char inByte5 = (x & 0xFF0000000000ul) >> 40;
    unsigned char inByte6 = (x & 0xFF000000000000ul) >> 48;
    unsigned char inByte7 = (x & 0xFF00000000000000ul) >> 56;
    return ((unsigned long)reversed_byte_table[inByte0] << 56) | ((unsigned long)reversed_byte_table[inByte1] << 48) | ((unsigned long)reversed_byte_table[inByte2] << 40) | ((unsigned long)reversed_byte_table[inByte3] << 32) | (reversed_byte_table[inByte4] << 24) | (reversed_byte_table[inByte5] << 16) | (reversed_byte_table[inByte6] << 8) | reversed_byte_table[inByte7];
  }

  inline static unsigned long reverse_bitwise(unsigned long x) {
    unsigned long maskFromLeft = 1<<LOG_N;
    unsigned long res = 0;
    unsigned int bitNum = LOG_N;
    while (maskFromLeft > 0) {
      unsigned char bit = (x & maskFromLeft) >> bitNum;
      res |= ( bit << (LOG_N-1-bitNum) );
      --bitNum;
      maskFromLeft >>= 1;
    }
    return res;
  }

  inline static unsigned long reverse_bytewise(unsigned long x) {
    // if (constexpr) statements should be eliminated by compiler to
    // choose correct case at compile time:

    if (LOG_N > sizeof(unsigned int)*8) {
      // Work with long reversal:
      
      // sizeof(unsigned long) * 8:
      //      const unsigned int bitsPerLong = sizeof(unsigned long)<<3;

      // Pure bit reversal of 1 will result in 1<<63; need to shift
      // right so that reversal of 1 yields LOG_N:
      return reverse_long_byte_table(x) >> (sizeof(unsigned long)*8 - LOG_N);
    }
    else if (LOG_N > sizeof(unsigned short)*8) {
      // Work with int reversal:
      
      // sizeof(unsigned int) * 8:
      //    const unsigned int bitsPerInt = sizeof(unsigned int)<<3;
      // Pure bit reversal of 1 will result in 1<<31; need to shift
      // right so that reversal of 1 yields LOG_N:
      return reverse_int_byte_table(x) >> (sizeof(unsigned int)*8 - LOG_N);
    }
    else if (LOG_N > sizeof(unsigned char)*8) {
      // Work with short int reversal:
      return reverse_short_byte_table(x) >> (sizeof(unsigned short)*8 - LOG_N);
    }
    // Work with char reversal:
    return reversed_byte_table[x] >> (sizeof(unsigned char)*8 - LOG_N);
  }

  inline static unsigned int reverse_bytewise(unsigned int x) {
    // To prevent unnecessary warnings about (desired) behavior when LOG_N == 0:
    if (LOG_N == 0)
      return x;
    return reverse_int_byte_table(x) >> (sizeof(unsigned int)*8 - LOG_N);
  }

  // Using XOR recurrence:
  inline static void advance_index_and_reversed(unsigned long & index, unsigned long & reversed) {
    unsigned long temp = index+1;
    unsigned long tail = ( index ^ temp );
    // tail is of the form 00...011...1

    index = temp;
    // create the reverse of tail, which is of form 11...100...0:
    auto shift = __builtin_clzl(tail);
    tail <<= shift;
    tail >>= ((sizeof(unsigned long)*8)-LOG_N);

    // xor reversed with reversed tail gives reversed of index+1:
    reversed ^= tail;
  }
};

template<unsigned char LOG_N>
const unsigned char BitReversal<LOG_N>::reversed_byte_table[256] = {0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2, 0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA, 0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE, 0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1, 0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5, 0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD, 0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB, 0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};

#endif
