#ifndef _VECTORTRIOT_HPP
#define _VECTORTRIOT_HPP

template <typename ...VECTORS>
void check_vector_pack_lengths(const VECTORS & ...args, unsigned long length) {
  #ifdef SHAPE_CHECK
  unsigned long sizes[] = { args.size()... };
  for (unsigned long s : sizes)
    assert(s >= length);
  #endif
}

// Note: Vectorizing functions also work on VectorView types; a
// common base class could be used, but that would require virtual
// functions, which would considerably slow the methods. So for now,
// this is performed with duck typing:

// Allows no modifications:
template <typename FUNCTION, typename ...VECTORS>
void for_each_vectors(FUNCTION function, unsigned long length, const VECTORS & ...args) {
  check_vector_pack_lengths<VECTORS...>(args..., length);
  for(unsigned long k=0; k<length; ++k) {
    function(args[k]...);
  }
}

// Allows modifications to all arguments:
template <typename FUNCTION, typename ...VECTORS>
void modify_vectors(FUNCTION function, unsigned long length, VECTORS & ...args) {
  check_vector_pack_lengths<VECTORS...>(args..., length);
  for(unsigned long k=0; k<length; ++k) {
    function(args[k]...);
  }
}

// Allows modifications only to dest:
template <typename FUNCTION, typename DEST_VECTOR, typename ...SOURCE_VECTORS>
void apply_vectors(FUNCTION function, unsigned long length, DEST_VECTOR & dest, const SOURCE_VECTORS & ...args) {
  check_vector_pack_lengths<DEST_VECTOR, SOURCE_VECTORS...>(dest, args..., length);
  for(unsigned long k=0; k<length; ++k) {
    function(dest[k], args[k]...);
  }
}

#endif
