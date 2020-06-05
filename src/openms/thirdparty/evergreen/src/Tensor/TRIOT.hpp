#ifndef _TRIOT_HPP
#define _TRIOT_HPP

#include "TemplateSearch.hpp"
#include "TensorUtils.hpp"

// TODO: Can a namespace casing like this improve compilation time
// by restricting lookup of template classes within this
// namespace?

// TODO: Would explicit template arguments enable faster compilation?
namespace TRIOT {
  //---------------------------------------------
  // For each:
  //---------------------------------------------

  template <unsigned char DIMENSION, unsigned char CURRENT>
  class ForEachFixedDimensionHelper {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      for (counter[CURRENT]=0; counter[CURRENT]<shape[CURRENT]; ++counter[CURRENT])
	TRIOT::ForEachFixedDimensionHelper<DIMENSION-1, CURRENT+1>::template apply<FUNCTION, TENSORS...>(counter, shape, function, args...);
    }
  };
  
  template <unsigned char CURRENT>
  class ForEachFixedDimensionHelper<1u, CURRENT> {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      for (counter[CURRENT]=0; counter[CURRENT]<shape[CURRENT]; ++counter[CURRENT])
	// For explicitly forcing the compiler to recognize that the
	// dimensionality is constant; this will be necessary on older
	// compilers:

	function(args[tuple_to_index_fixed_dimension<CURRENT+1>(counter, args.data_shape())]...);
	//  function(args[tuple_to_index(counter, args.data_shape(), CURRENT+1)]...);
    }
  };

    
  // for a tensor with dimension 0
  template <unsigned char CURRENT>
  class ForEachFixedDimensionHelper<0u, CURRENT> {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      // Do nothing
    }
  };

  template <unsigned char DIMENSION>
  class ForEachFixedDimension {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      unsigned long counter[DIMENSION];
      memset(counter, 0, DIMENSION*sizeof(unsigned long));
      TRIOT::ForEachFixedDimensionHelper<DIMENSION,0>::template apply<FUNCTION, TENSORS...>(counter, shape, function, args...);
    }
  };
  
  template<>
    class ForEachFixedDimension<0U>	{
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(const_tup_t shape, FUNCTION function, TENSORS & ...args) {
		// do nothing, so that memset is not called with size = 0 which is a GCC extension
	}
  };

  //---------------------------------------------
  // For each, with visible counter:
  //---------------------------------------------

  template <unsigned char DIMENSION, unsigned char CURRENT>
  class ForEachVisibleCounterFixedDimensionHelper {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      for (counter[CURRENT]=0; counter[CURRENT]<shape[CURRENT]; ++counter[CURRENT])
	ForEachVisibleCounterFixedDimensionHelper<DIMENSION-1, CURRENT+1>::template apply<FUNCTION, TENSORS...>(counter, shape, function, args...);
    }
  };
  
  template <unsigned char CURRENT>
  class ForEachVisibleCounterFixedDimensionHelper<1u, CURRENT> {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      for (counter[CURRENT]=0; counter[CURRENT]<shape[CURRENT]; ++counter[CURRENT])
	// Cast the counter to a const_tup_t pointer so that its
	// contents cannot be modified by function:

	// For explicitly forcing the compiler to recognize that the
	// dimensionality is constant; this will be necessary on older
	// compilers (which may not see this through constant
	// propagation):

	function(static_cast<const_tup_t>(counter), CURRENT+1, args[tuple_to_index_fixed_dimension<CURRENT+1>(counter, args.data_shape())]...);

	// function(static_cast<const_tup_t>(counter), CURRENT+1, args[tuple_to_index(counter, args.data_shape(), CURRENT+1)]...);
    }
  };

  template <unsigned char CURRENT>
  class ForEachVisibleCounterFixedDimensionHelper<0u, CURRENT> {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(tup_t counter, const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      // Do nothing
    }
  };

  template <unsigned char DIMENSION>
  class ForEachVisibleCounterFixedDimension {
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(const_tup_t shape, FUNCTION function, TENSORS & ...args) {
      unsigned long counter[DIMENSION];
      memset(counter, 0, DIMENSION*sizeof(unsigned long));
      ForEachVisibleCounterFixedDimensionHelper<DIMENSION,0>::template apply<FUNCTION, TENSORS...>(counter, shape, function, args...);
    }
  };
  
  template<>
    class ForEachVisibleCounterFixedDimension<0U>	{
  public:
    template <typename FUNCTION, typename ...TENSORS>
    inline static void apply(const_tup_t shape, FUNCTION function, TENSORS & ...args) {
		// do nothing, so that memset is not called with size = 0 which is a GCC extension
	}
  };
}

template <typename ...TENSORS>
void check_tensor_pack_bounds(const TENSORS & ...args, const Vector<unsigned long> & shape) {
  #ifdef SHAPE_CHECK
  // Verify same shapes:

  // TODO: this could be faster by using an array of references; C++
  // does not allow an array of references, but it would allow an
  // array of structs containing only the reference. 
  Vector<unsigned long> shapes[] = { args.view_shape()... };
  for (const Vector<unsigned long> & s : shapes) {
    // Check that all dimensions match:
    assert(s.size() == shape.size());

    // Check that iterating over shape is in bounds with respect to
    // the current view_shape:
    assert(s >= shape);
  }
  #endif
}

template <typename ...TENSORS>
void check_tensor_pack_bounds(const Vector<unsigned long> & shape) {
}

template <typename ...TENSORS>
Vector<unsigned long> bounding_shape(const TENSORS & ...args) {
  // Verify same shapes:
  Vector<unsigned long> shapes[] = { args.view_shape()... };
  Vector<unsigned long> result = shapes[0];
  for (const Vector<unsigned long> & s : shapes) {
    #ifdef SHAPE_CHECK
    // Check that all dimensions match:
    assert(s.size() == result.size());
    #endif

    for (unsigned int i=0; i<result.size(); ++i)
      result[i] = std::min(result[i], s[i]);
  }
  return result;
}

// Interface for external use (note: these functions also work with
// the tensor view types); rather than TENSOR types, they could be
// treated as TensorLike and WritableTensorLike, but for now duck
// typing is a bit simpler. This works because const & parameter
// typing allows rvalues and lvalues to be passed, and && typing would
// normally only allow rvalue references, but because of templating,
// the compiler can automatically choose the type as T& && --> T&.

// Allows no modifications:
template <typename FUNCTION, typename ...TENSORS>
void for_each_tensors(FUNCTION function, const Vector<unsigned long> & shape, const TENSORS & ...args) {
  check_tensor_pack_bounds<TENSORS...>(args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachFixedDimension>::apply(shape.size(), shape, function, args...);
}

template <typename FUNCTION, typename ...TENSORS>
void enumerate_for_each_tensors(FUNCTION function, const Vector<unsigned long> & shape, const TENSORS & ...args) {
  check_tensor_pack_bounds<TENSORS...>(args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachVisibleCounterFixedDimension>::apply(shape.size(), shape, function, args...);
}

// Allows modifications to all arguments:
template <typename FUNCTION, typename ...DEST_TENSORS>
void modify_tensors(FUNCTION function, const Vector<unsigned long> & shape, DEST_TENSORS && ...args) {
  check_tensor_pack_bounds<DEST_TENSORS...>(args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachFixedDimension>::apply(shape.size(), shape, function, args...);
}

template <typename FUNCTION, typename ...TENSORS>
void enumerate_modify_tensors(FUNCTION function, const Vector<unsigned long> & shape, TENSORS && ...args) {
  check_tensor_pack_bounds<TENSORS...>(args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachVisibleCounterFixedDimension>::apply(shape.size(), shape, function, args...);
}

// Allow modifications only to dest:
template <typename FUNCTION, typename DEST_TENSOR, typename ...SOURCE_TENSORS>
void apply_tensors(FUNCTION function, const Vector<unsigned long> & shape, DEST_TENSOR && dest, const SOURCE_TENSORS & ...source_args) {
  check_tensor_pack_bounds<DEST_TENSOR, SOURCE_TENSORS...>(dest, source_args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachFixedDimension>::apply(shape.size(), shape, function, dest, source_args...);
}

template <typename FUNCTION, typename DEST_TENSOR, typename ...SOURCE_TENSORS>
void enumerate_apply_tensors(FUNCTION function, const Vector<unsigned long> & shape, DEST_TENSOR && dest, const SOURCE_TENSORS & ...source_args) {
  check_tensor_pack_bounds<DEST_TENSOR, SOURCE_TENSORS...>(dest, source_args..., shape);
  LinearTemplateSearch<0u,MAX_TENSOR_DIMENSION,TRIOT::ForEachVisibleCounterFixedDimension>::apply(shape.size(), shape, function, dest, source_args...);
}

#endif
