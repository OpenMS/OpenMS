#ifndef _P_CONVOLVE_HPP
#define _P_CONVOLVE_HPP

// See Pfeuffer and Serang 2016 JMLR for details. Note that this
// implementation does not yet perform the linear regression
// correction from the reference above (which would further improve
// accuracy).

#include <array>
#include <set>
#include "naive_convolve.hpp"
#include "fft_convolve.hpp"

// Used in numeric p-convolution:
const double tau_denom = 1e-9;

// Empirically chosen runtime constant of FFT convolution:
const double FFT_CONV_RUNTIME_CONSTANT = 10.0;

// Speedup in runtime constant of naive when p=1 or p=inf relative to
// naive for arbitrary p (because they do not need to compute pow):

const double SPEEDUP_OF_NAIVE_P1_OR_PINF = 2.0;

// The reasoning behind possible choices of p^*_max is described in
// Pfeuffer and Serang 2016 JMLR (it is a heuristically reasonable
// compromise between accuracy and speed):
const double MAX_P_NUMERIC = log2(0.7) / log2(0.999) * 2;

inline Tensor<double> fft_p_convolve_to_p(const Tensor<double> & lhs, const Tensor<double> & rhs, double p_goal) {
  // Note: does not check if values are nonnegative (this is done one
  // time in numeric_p_convolve so as to not check repeatedly).
  
  Tensor<double> lhs_pow = lhs, rhs_pow = rhs;
  for (unsigned long k=0; k<lhs_pow.flat_size(); ++k)
    lhs_pow[k] = custom_pow(lhs_pow[k], p_goal);
  for (unsigned long k=0; k<rhs_pow.flat_size(); ++k)
    rhs_pow[k] = custom_pow(rhs_pow[k], p_goal);
  
  Tensor<double> res = fft_convolve(lhs_pow, rhs_pow);
  for (unsigned long k=0; k<rhs_pow.flat_size(); ++k)
    res[k] = fabs(res[k]);
  return res;
}

// Note: this iterative method for computing pow makes the algorithm
// in O(n log(n) |P|^2); however, this is still significantly faster
// than using pow.
inline double fast_pow_from_interleaved_p_index(double val, unsigned int p_ind) {
  for (unsigned int i=0; i<p_ind/2; ++i)
    val *= val;
  
  // Interleaved powers of 2: e.g., p=[1,1.5,2,3,4,6,8 ...], so
  // p_ind=5 corresponds to 6, which is not a power of 2, so perform
  // sqrt(val^3) for final step:

  // Faster version of if (p_ind % 2 == 1):
  if ( (p_ind & 1u) == 1)
    val = sqrt(val*val*val);
  return val;
}

inline Tensor<double> fft_p_convolve_to_p_from_p_index(const Tensor<double> & lhs, const Tensor<double> & rhs, unsigned int p_ind) {
  // Note: does not check if values are nonnegative (this is done one
  // time in numeric_p_convolve so as to not check repeatedly).

  // Note: could possibly pass these by reference to prevent multiple
  // allocations.
  Vector<unsigned long> conv_shape_doubles = padded_convolution_shape(lhs, rhs);
  Tensor<double> lhs_padded_doubles(conv_shape_doubles);
  Tensor<double> rhs_padded_doubles(conv_shape_doubles);

  // Embed in larger and take to p:
  apply_tensors([p_ind](double & lhs_padded, double lhs) {
      lhs_padded = fast_pow_from_interleaved_p_index(lhs, p_ind);
    },
    lhs.view_shape(),
    lhs_padded_doubles, lhs);

  apply_tensors([p_ind](double & rhs_padded, double rhs) {
      rhs_padded = fast_pow_from_interleaved_p_index(rhs, p_ind);
    },
    rhs.view_shape(),
    rhs_padded_doubles, rhs);

  Tensor<double> res = fft_convolve_already_padded_rvalue(std::move(lhs_padded_doubles), std::move(rhs_padded_doubles), lhs.data_shape() + rhs.data_shape() - 1ul);

  for (unsigned long k=0; k<res.flat_size(); ++k)
    res[k] = fabs(res[k]);
  return res;
}

inline Vector<double> interleaved_powers_of_2(unsigned int log_max_p) {
  Vector<double> result(2*log_max_p+1);

  // Length should be odd (final power of two will occur outside the loop):
  double val = 1.0;
  unsigned int k;
  for (k=0; k<(result.size()-1)/2; ++k) {
    result[2*k] = val;
    result[2*k+1] = val*1.5;
    val *= 2.0;
  }
  result[2*k] = val;
  return result;
}

inline double best_tau_fft_for_length(unsigned long flat_length) {
  // Note: Based on the empirical data (mentioned below), using a
  // statically defined tau of 1e-5 should be roughly large enough
  // for the largest supported FFT with length 2^32; however, it
  // will be unnecessarily strict
  //return 1e-5;
    
  // Based on empirical error data (using the maximum error from the
  // lowest and highest values in a list) on a large range of 1D
  // sizes. Note that there is likely a theoretical basis for
  // believing the error will grow logarithmically with the length, as
  // mentioned very briefly in Pfeuffer and Serang 2016.

  const double log_x1 = log(1e3);
  const double log_y1 = log(2e-12);
  const double min_tau = 1e-12;
  const double log_x2 = log(4e6);
  const double log_y2 = log(2e-8);
  const double slope = (log_y2-log_y1) / (log_x2-log_x1);
  const double bias = log_y1 - slope*log_x1;
  const double tau_val = exp(bias  + slope*log(double(flat_length)));
  // Using 15* expected optimal tau value to be conservative in case
  // the log vs. log linear curve is slightly jagged or is slightly
  // concave up (also the benchmark data used to fit the error used
  // complex FFT convolution, whereas this now uses real FFT for real
  // convolutions, which uses a recurrence that lowers the accuracy
  // slightly compared to using a cache of complex exponentials). If
  // necessary, more recent empirical data can be taken. Also, note
  // that this may be conservative since it uses the flat size for
  // multidimensional problems (whereas it is unknown but possible
  // that the numeric stability may be limited by the largest axis,
  // etc.):

  // Note: the FFT estimates of the norms should be monotonic whenever
  // the FFT is stable (i.e., the p-norm should be >= the 2p-norm). If
  // not, the p-norm estimate may work, but the quadratic and linear
  // estimates may be unstable.

  return std::max(15*tau_val, min_tau);
}

inline static double linear_projection(const std::array<double,2> & norms, double p1, double p2, double p_goal) {
  // assumes norms[0] is not too close to zero (single stable norm should never = 0)
  double delta = p2-p1;
  double root = norms[1] / norms[0];
  if ( fabs(root) < tau_denom )
    // when it's unstable to solve the linear (via division), return
    // the p-norm estimate using the highest stable p (note that the
    // sequence of p was chosen to not exceed p_goal):
    return custom_pow(norms[1], 1.0/p2);
  double alpha = custom_pow(root, 1.0/delta);
  double n = norms[0] / custom_pow(alpha, p1);
  return alpha * custom_pow(n, 1.0/p_goal);
}

inline double check_nan_call_linear_projection(double val, const std::array<double,4> & norms, double p1, double p2, double p_goal){
  if(std::isnan(val)){
    const std::array<double,2> lin_norms{{norms[2], norms[3]}};
    return linear_projection(lin_norms, p1, p2, p_goal);
  }
  return val;
}

inline double quadratic_projection(const std::array<double,4> & norms, double p1, double p2, double p_goal) {
  double delta = p2-p1;
  
  // Quadrataic coefficients from the null space of the matrix of norms:
  double c = norms[1]*norms[3] - norms[2]*norms[2];
  double b = norms[1]*norms[2] - norms[0]*norms[3];
  double a = norms[0]*norms[2] - norms[1]*norms[1];
  // Solve c + b*x + a*x^2 == 0.

  if ( fabs(a) > tau_denom ) {
    
    double disc = b*b - 4*a*c;

    if ( disc >= 0.0 ) {
      double root1 = (-b + sqrt(disc)) / (2*a);
      double root2 = (-b - sqrt(disc)) / (2*a);

      if ( root1 >= 0.0 && root2 >= 0.0 ) {
        // only perform when numerically stable:
        double alpha1 = custom_pow(root1, 1.0/delta);
        double alpha2 = custom_pow(root2, 1.0/delta);
        // Ensure alpha1 > alpha2:
        if ( alpha2 > alpha1 )
          std::swap(alpha1, alpha2);
	  
        double alpha1_up_p1 = custom_pow(alpha1,p1);
        double alpha1_up_p2 = custom_pow(alpha1,p2);
        
        double alpha2_up_p1 = custom_pow(alpha2,p1);
        double alpha2_up_p2 = custom_pow(alpha2,p2);
        
        double denom = alpha1_up_p2*alpha2_up_p1 - alpha1_up_p1*alpha2_up_p2;
        if ( fabs(denom) > tau_denom ) {
          double n1 = ( norms[1]*alpha2_up_p1 - norms[0]*alpha2_up_p2 ) / denom;
          double n2 = ( norms[0]*alpha1_up_p2 - norms[1]*alpha1_up_p1 ) / denom;

          if ( alpha1 > tau_denom )
            return check_nan_call_linear_projection(alpha1 * custom_pow( n1 + n2*custom_pow(alpha2/alpha1,p_goal), 1.0/p_goal ), norms, p1, p2, p_goal);
          return check_nan_call_linear_projection(custom_pow(n1 * custom_pow(alpha1,p_goal) + n2 * custom_pow(alpha2,p_goal), 1.0/p_goal), norms, p1, p2, p_goal);
        }
      }
    }
  }
  
  const std::array<double,2> lin_norms{{norms[2], norms[3]}};
  return linear_projection(lin_norms, p1, p2, p_goal);
}

inline void compute_quadratic_projections(const std::vector<Tensor<double> > & p_index_to_norms, const Vector<double> & all_p, double p_goal, Tensor<double> & result, const Tensor<bool> & solved, const Tensor<int> & highest_stable_p_index) {
  // Note: It may potentially be worth tranposing so that norms for a
  // given index are in cache order.

  // Fill in remaining entries in result:
  for (unsigned long i=0; i<result.flat_size(); ++i)
    if ( ! solved[i] ) {
      // Compute quadratic projection of p-norm:
      int highest_stable = highest_stable_p_index[i];
      double result_at_index;
      // Powers of 2 are even indices (highest in sequence of 4
      // evenly spaced points must be a power of 2):

      // Note: it may be more efficient to replace this if...else
      // ladder with a switch statement:

      // Use & 1 as a faster replacement for % 2:
      if (highest_stable >= 4 && (highest_stable & 1u) == 0) {
        // 5 points available {0,1,2,3,4,...} --> 4 evenly spaced points
        std::array<double, 4> norms{{p_index_to_norms[highest_stable-4][i],p_index_to_norms[highest_stable-2][i],p_index_to_norms[highest_stable-1][i],p_index_to_norms[highest_stable][i]}};
        result_at_index = quadratic_projection(norms, all_p[highest_stable-1], all_p[highest_stable], p_goal);
      }
      // Use & 1 as a faster replacement for % 2:
      else if (highest_stable >= 5 && (highest_stable & 1u) == 1) {
        // Use highest_stable-1, since you're dropping by 1 to reach the next power of 2:
        std::array<double, 4> norms{{p_index_to_norms[highest_stable-5][i],p_index_to_norms[highest_stable-3][i],p_index_to_norms[highest_stable-2][i],p_index_to_norms[highest_stable-1][i]}};

        result_at_index = quadratic_projection(norms, all_p[highest_stable-2], all_p[highest_stable-1], p_goal);
      }
      else if (highest_stable >= 1) {
        // 2 evenly points available
        std::array<double, 2> norms{{p_index_to_norms[highest_stable-1][i],p_index_to_norms[highest_stable][i]}};
        result_at_index = linear_projection(norms, all_p[highest_stable-1], all_p[highest_stable], p_goal);
      }
      else {
        // only one point stable
	// TODO: could divide by length of u vector to improve error:
	// (\| u \|_p^p / len(u))^(1/p)
        result_at_index = custom_pow(p_index_to_norms[highest_stable][i], 1.0/all_p[highest_stable]);
      }
      
      result[i] = result_at_index;
    }
}

// Note: This function is slow, but it should be called quite rarely. 
inline double naive_p_convolve_at_index(const Tensor<double> & lhs, const Tensor<double> & rhs, const Vector<unsigned long> & ind, double p_goal) {
  double max_val = 0.0;
  Vector<unsigned long> rhs_ind(ind.size());
  enumerate_for_each_tensors([&ind, &rhs_ind, &rhs, &max_val](const_tup_t lhs_tup, const unsigned char dim, double lhs_val){
      for (unsigned char i=0; i<dim; ++i)
	rhs_ind[i] = ind[i] - lhs_tup[i];

      // Note: Using TensorView once would be much faster:
      if (rhs_ind < rhs.data_shape())
	max_val = std::max(max_val, lhs_val*rhs[rhs_ind]);
    },
    lhs.data_shape(),
    lhs);

  if (max_val == 0.0)
    return max_val;

  // Divide by max_val before taking to power p_goal to preserve precision:
  double res = 0.0;
  enumerate_for_each_tensors([&ind, &rhs_ind, &rhs, max_val, &res, p_goal](const_tup_t lhs_tup, const unsigned char dim, double lhs_val){
      for (unsigned char i=0; i<dim; ++i)
	rhs_ind[i] = ind[i] - lhs_tup[i];

      // Note: Using TensorView once would be much faster:
      if (rhs_ind < rhs.data_shape())
	res += custom_pow(lhs_val*rhs[rhs_ind]/max_val, p_goal);
    },
    lhs.data_shape(),
    lhs);

  return max_val*custom_pow(res, 1.0/p_goal);
}

inline void perform_affine_correction(const Tensor<double> & lhs, const Tensor<double> & rhs, const double p_goal, const Tensor<int> & highest_stable_p_index, Tensor<double> & result) {
  // Note: for greater performance, this could be a bitset:
  std::set<int> used_p_indices;
  for (unsigned long i=0; i<result.flat_size(); ++i)
    used_p_indices.insert(highest_stable_p_index[i]);

  for (int p_ind : used_p_indices) {

    // Find min and max result values and their tuple indices:
    double min_res_in_contour = std::numeric_limits<double>::infinity();
    Vector<unsigned long> min_index(result.dimension());
    double max_res_in_contour = 0.0;
    Vector<unsigned long> max_index(result.dimension());
    enumerate_for_each_tensors([&min_res_in_contour, &min_index, &max_res_in_contour, &max_index, p_ind](const_tup_t tup, const unsigned char dim, double res_val, int res_p_ind) {

	if (res_p_ind == p_ind) {
	  if (res_val < min_res_in_contour) {
	    min_res_in_contour = res_val;
	    for (unsigned char i=0; i<dim; ++i)
	      min_index[i] = tup[i];
	  }
	  if (res_val > max_res_in_contour) {
	    max_res_in_contour = res_val;
	    for (unsigned char i=0; i<dim; ++i)
	      max_index[i] = tup[i];
	  }
	}
      },
      result.data_shape(),
      result, highest_stable_p_index);

    // Compute exact p-convolution at min_index and max_index:
    double exact_at_min_index = naive_p_convolve_at_index(lhs, rhs, min_index, p_goal);
    double exact_at_max_index = naive_p_convolve_at_index(lhs, rhs, max_index, p_goal);

    const double denom = max_res_in_contour - min_res_in_contour;
    if (denom > tau_denom) {
      const double contour_slope = (exact_at_max_index - exact_at_min_index) / denom;
      // new_est = contour_bias + contour_slope * old_est
      // new_est = exact_at_min_index + contour_slope * (old_est - min_res_in_contour)
      // contour_bias = exact_at_min_index - contour_slope*min_res_in_contour
      const double contour_bias = exact_at_min_index - contour_slope*min_res_in_contour;
      for (unsigned long i=0; i<result.flat_size(); ++i) {
	if (highest_stable_p_index[i] == p_ind) {
	  double & res_val = result[i];
	  res_val = contour_bias + res_val*contour_slope;
	}
      }
    }
  }
}


inline Tensor<double> numeric_p_convolve_helper(const Tensor<double> & lhs, const Tensor<double> & rhs, double max_p, double p_goal) {
  if (p_goal >= 1.0) {
    // Note: max_p should already be <= p_goal
    unsigned int log_max_p = log2(max_p);
    Vector<double> all_p = interleaved_powers_of_2( log_max_p );

    double tau_fft = best_tau_fft_for_length( std::max(lhs.flat_size(), rhs.flat_size()) );

    Vector<unsigned long> res_shape = lhs.data_shape() + rhs.data_shape() - 1ul;
    
    Tensor<int> highest_stable_p_index(res_shape);
    // Use std::vector as a container: Vector<T> is for numeric T types
    // only and does not call constructors, etc.:
    std::vector<Tensor<double> > p_index_to_norms(all_p.size());

    // Initialized to false, since default value is 0:
    Tensor<bool> solved(res_shape);

    Tensor<double> result;
  
    // Note: isinf may not work on every compiler when aggressive
    // optimizations are turned on:
    if ( ! std::isinf(p_goal) ) {
    
      /* p-norm is not included in powers of two sequence. therefore,
	 try p-norm directly here (you may get lucky and have all
	 results numerically stable with p, and so you only need a
	 single FFT) */

      // Note: can avoid ^p when p is a power of 2 here:
      result = fft_p_convolve_to_p(lhs, rhs, p_goal);
    
      for (unsigned long k=0; k<result.flat_size(); ++k)
	if ( result[k] > tau_fft ) {
	  result[k] = custom_pow(result[k], 1.0/p_goal);
	  solved[k] = true;
	}
    
      if (all(solved.flat()))
	return result;

      // If p_goal is in the power of two sequence, then add it to the
      // collection so that the convolution is not recomputed:
      // (Do not use custom_pow, because it may be less stable)
      if ( 1ul<<(unsigned int)log_max_p == p_goal )
	p_index_to_norms[p_index_to_norms.size()-1] = result;
    }
    else
      // result has not yet been allocated:
      result = Tensor<double>(res_shape);
  
    Tensor<bool> sufficient_for_projection = solved;

    // Note: this loop is a strong candidate for parallelization
    // (although it might require some massaging with OpenMP, since
    // this loop can return as soon as all indices have been
    // sufficiently computed)
    for (int p_ind=(int)all_p.size()-1; p_ind>=0; --p_ind) {

      // If p_goal is in all_p, it has already been added to p_index_to_norms:
      if ( all_p[p_ind] != p_goal ) {
	// The following replaces calling with pow:
	//      p_index_to_norms[p_ind] = fft_p_convolve_to_p(lhs, rhs, p_ind);      
	p_index_to_norms[p_ind] = fft_p_convolve_to_p_from_p_index(lhs, rhs, p_ind);
      }
    
      const Tensor<double> & row = p_index_to_norms[p_ind];
      for (unsigned long k=0; k<highest_stable_p_index.flat_size(); ++k)
	if ( ! solved[k] )
	  // Since it is going in descending order, stop assigning to
	  // highest_stable_p_index if it has already been assigned a
	  // nonzero value:
	  if ( highest_stable_p_index[k] == 0 && row[k] > tau_fft )
	    highest_stable_p_index[k] = p_ind;

      for (unsigned long k=0; k<highest_stable_p_index.flat_size(); ++k)
	// powers of 2 are the even indices (highest in sequence of 4
	// evenly spaced points must be a power of 2):
	// Use & 1 == 0 as a faster replacement for % 2 == 0:
	if ( ( (highest_stable_p_index[k] & 1u) == 0 && highest_stable_p_index[k] - p_ind >= 5 ) || highest_stable_p_index[k] - p_ind >= 6 ) {
	  sufficient_for_projection[k] = true;
	}

      // Check if enough contours of p have been solved at every index
      // so that the desired level of extrapolation; if so, break:

      // Note: A potentially large speedup may be possible here by
      // stopping once n - cardinality(sufficient_for_projection) is <
      // some constant; those remaining constant number of indices
      // could be solved each in O(n), potentially preventing a full
      // O(n log(n)) convolution from being computed for smaller
      // values of p.
      if (all(sufficient_for_projection.flat()))
	break;
    }
  
    compute_quadratic_projections(p_index_to_norms, all_p, p_goal, result, solved, highest_stable_p_index);

    perform_affine_correction(lhs, rhs, p_goal, highest_stable_p_index, result);

    return result;
  }
  else {
    #ifndef NUMERIC_CHECK
    // Note: Negative p_goal could be processed by taking 1.0/value
    // for each value and using positive p.
    assert(p_goal >= 0);
    #endif

    // p < 1

    // Assume that for p<1, raw p-convolution (without piecewise) is
    // numerically stable (it will decrease the dynamic range of the data).
    Tensor<double> result = fft_p_convolve_to_p(lhs, rhs, p_goal);
    for (unsigned long i=0; i<result.flat_size(); ++i)
      result[i] = custom_pow(result[i], 1.0/p_goal);
    return result;
  }
}

inline Tensor<double> numeric_p_convolve(const Tensor<double> & lhs, const Tensor<double> & rhs, double p_goal) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  #endif
  #ifdef NUMERIC_CHECK

  // Numeric p-convolution is only for nonnegative values; it may be
  // possible on all values by encoding the sign into complex numbers
  // and then performing a complex convolution instead of a real
  // convolution:
  assert( lhs.flat() >= 0.0 );
  assert( rhs.flat() >= 0.0 );

  // Note: When p_goal is negative, believe could use 1/lhs and 1/rhs with -p_goal:
  assert(p_goal > 0.0);
  #endif

  Vector<unsigned long> res_shape = lhs.data_shape() + rhs.data_shape() - 1ul;
  // Empirically, at a flat size of 64, the fft version becomes more
  // efficient than naive. But because of doubling all but the last
  // axis in real fft convolution (for zero padding), there is roughly
  // a 2^(dim-1) slowdown on fft version (using a linear runtime for
  // fft, which is realistic for small lengths where the log2 factor
  // should be negligable for an efficient implementation).
  unsigned long flat_size = flat_length(static_cast<unsigned long*>(res_shape), res_shape.size());

  double max_p = std::min(p_goal, MAX_P_NUMERIC);

  double approx_fft_conv_runtime = flat_size*log2(flat_size)*log2(max_p)*FFT_CONV_RUNTIME_CONSTANT;
  double approx_naive_conv_runtime = flat_size*flat_size;
  
  if (p_goal == 1.0) {
    if ( approx_fft_conv_runtime*SPEEDUP_OF_NAIVE_P1_OR_PINF > approx_naive_conv_runtime ){
      return naive_convolve(lhs, rhs);
    }
  }
  else if (std::isinf(p_goal)) {
    // When true max-convolution is wanted, then the cost of computing
    // naive version is slightly improved, because it does not need to
    // use pow function.
    if ( approx_fft_conv_runtime*SPEEDUP_OF_NAIVE_P1_OR_PINF > approx_naive_conv_runtime ){
      return naive_max_convolve(lhs, rhs);
    }
  }
  else {
    if ( flat_size*log2(flat_size)*log2(max_p)*FFT_CONV_RUNTIME_CONSTANT > flat_size*flat_size ){
      return naive_p_convolve(lhs, rhs, p_goal);
    }
  }
  
  double lhs_max = max(lhs.flat());
  double rhs_max = max(rhs.flat());

  // Don't divide by 0; if lhs_max or rhs_max == 0, return Tensor full of 0 values:
  if (lhs_max == 0.0 || rhs_max == 0.0)
    return Tensor<double>(lhs.data_shape() + rhs.data_shape() - 1ul);
  
  Tensor<double> lhs_prime = lhs;
  lhs_prime.flat() /= lhs_max;
  Tensor<double> rhs_prime = rhs;
  rhs_prime.flat() /= rhs_max;
  
  Tensor<double> res = numeric_p_convolve_helper(lhs_prime, rhs_prime, max_p, p_goal);
  
  res.flat() *= lhs_max*rhs_max;

  // Correct instability for values very close to zero: when the phase
  // can make the sign negative, force the value to be nonnegative:
  for (unsigned long k=0; k<res.flat_size(); ++k)
    res[k] = fabs(res[k]);

  return res;
}

#endif
