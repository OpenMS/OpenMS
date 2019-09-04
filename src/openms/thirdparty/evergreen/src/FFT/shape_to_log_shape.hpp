#ifndef _SHAPE_TO_LOG_SHAPE_HPP
#define _SHAPE_TO_LOG_SHAPE_HPP

inline unsigned char integer_log2(unsigned long val) {
  // Note: that can be sped up using bit twiddling (cast to float,
  // etc.) or using the clz opcode. But it's still reasonably fast as
  // is:
  unsigned char res = round( log2(val) );
  #ifdef SHAPE_CHECK
  assert( (1ul<<res) == val);
  #endif
    
  return res;
}

inline Vector<unsigned char> shape_to_log_shape(const Vector<unsigned long> & shape) {
  Vector<unsigned char> log_shape(shape.size());
  for (unsigned char k=0; k<shape.size(); ++k)
    log_shape[k] = integer_log2(shape[k]);
  return log_shape;
}

inline unsigned long real_length_to_packed_length(unsigned long len) {
  if (len == 0)
    return 0;

  return len / 2 + 1;
}

inline unsigned long packed_length_to_real_length(unsigned long packed_len) {
  if (packed_len == 0)
    return 0;

  if (packed_len == 1)
    return 1;

  return (packed_len - 1)*2;
}

// Note: returns log_shape for equivalent complex FFT:
inline Vector<unsigned char> packed_shape_to_log_shape(const Vector<unsigned long> & packed_shape) {
  Vector<unsigned char> log_equiv_shape(packed_shape.size());
  unsigned char k;
  for (k=0; k<packed_shape.size()-1; ++k)
    log_equiv_shape[k] = integer_log2( packed_shape[k] );
  log_equiv_shape[k] = integer_log2( packed_length_to_real_length(packed_shape[k]) );

  return log_equiv_shape;
}

inline Vector<unsigned char> reversed_packed_shape_to_log_shape(const Vector<unsigned long> & packed_shape) {
  Vector<unsigned char> log_equiv_shape(packed_shape.size());
  unsigned char k=0;
  log_equiv_shape[k] = integer_log2( packed_length_to_real_length(packed_shape[k]) );
  for (k=1; k<packed_shape.size(); ++k)
    log_equiv_shape[k] = integer_log2( packed_shape[k] );

  return log_equiv_shape;
}

#endif
