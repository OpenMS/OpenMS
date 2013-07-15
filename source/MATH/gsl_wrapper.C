/*
 * gsl_wrapper.c
 *
 *  Created on: Jul 10, 2013
 *      Author: Hans-Christian Ehrlich
 */

#include "OpenMS/MATH/gsl_wrapper.h"

#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

struct deprecated_gsl_bspline_deriv_workspace : public gsl_bspline_deriv_workspace {};
struct deprecated_gsl_bspline_workspace : public gsl_bspline_workspace {};
struct deprecated_gsl_fft_real_wavetable : public gsl_fft_real_wavetable {};
struct deprecated_gsl_fft_real_workspace : public gsl_fft_real_workspace {};
struct deprecated_gsl_interp : public gsl_interp {};
struct deprecated_gsl_interp_accel : public gsl_interp_accel {};
struct deprecated_gsl_interp_type : public gsl_interp_type {};
struct deprecated_gsl_matrix : public gsl_matrix {};
struct deprecated_gsl_multifit_fdfsolver : public gsl_multifit_fdfsolver {};
struct deprecated_gsl_multifit_fdfsolver_type : public gsl_multifit_fdfsolver_type {};
struct deprecated_gsl_multifit_function_fdf : public gsl_multifit_function_fdf {};
struct deprecated_gsl_multifit_linear_workspace : public gsl_multifit_linear_workspace {};
struct deprecated_gsl_permutation : public gsl_permutation {};
struct deprecated_gsl_ran_discrete_t : public gsl_ran_discrete_t {};
struct deprecated_gsl_rng : public gsl_rng {};
struct deprecated_gsl_rng_type : public gsl_rng_type {};
struct deprecated_gsl_spline : public gsl_spline {};
struct deprecated_gsl_vector : public gsl_vector {};
struct deprecated_gsl_vector_view : public gsl_vector_view {};

void
deprecated_wrapper_gsl_rng_default_seed_set( unsigned long int x )
{
	gsl_rng_default_seed = x;
}

const deprecated_gsl_rng_type *
deprecated_wrapper_gsl_rng_taus_get()
{
	return static_cast<const deprecated_gsl_rng_type*>(gsl_rng_taus);
}

const deprecated_gsl_rng_type * deprecated_wrapper_get_gsl_rng_mt19937()
{
	return static_cast<const deprecated_gsl_rng_type*>( gsl_rng_mt19937 );
}

const deprecated_gsl_multifit_fdfsolver_type * deprecated_wrapper_get_multifit_fdfsolver_lmsder()
{
	return static_cast<const deprecated_gsl_multifit_fdfsolver_type*>(gsl_multifit_fdfsolver_lmsder);
}

const deprecated_gsl_rng_type * deprecated_wrapper_get_gsl_rng_default()
{
	return static_cast<const deprecated_gsl_rng_type*>(gsl_rng_default);
}

const deprecated_gsl_vector*
deprecated_wrapper_gsl_vector_view_get_vector(deprecated_gsl_vector_view_ptr v)
{
	return static_cast<deprecated_gsl_vector*>( &v->vector );
}

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_x( deprecated_gsl_multifit_fdfsolver* s )
{
	return static_cast<deprecated_gsl_vector*>(s->x);
}

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_f( deprecated_gsl_multifit_fdfsolver* s )
{
	return static_cast<deprecated_gsl_vector*>(s->f);
}

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_dx( deprecated_gsl_multifit_fdfsolver* s )
{
	return static_cast<deprecated_gsl_vector*>(s->dx);
}

deprecated_gsl_matrix*
deprecated_wrapper_gsl_multifit_fdfsolver_get_J( deprecated_gsl_multifit_fdfsolver* s )
{
	return static_cast<deprecated_gsl_matrix*>(s->J);
}

size_t
deprecated_wrapper_gsl_matrix_get_size1( deprecated_gsl_matrix* M)
{
	return M->size1;
}

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_linear()
{
	return static_cast<const deprecated_gsl_interp_type*>( gsl_interp_linear );
}

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_polynomial()
{
	return static_cast<const deprecated_gsl_interp_type*>( gsl_interp_polynomial );
}

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_cspline()
{
	return static_cast<const deprecated_gsl_interp_type*>( gsl_interp_cspline );
}

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_akima()
{
	return static_cast<const deprecated_gsl_interp_type*>( gsl_interp_akima );
}

size_t
deprecated_wrapper_gsl_interp_type_get_min_size( const deprecated_gsl_interp_type * t )
{
	return t->min_size;
}

deprecated_gsl_vector *
deprecated_wrapper_gsl_bspline_workspace_get_knots( deprecated_gsl_bspline_workspace * w )
{
	return static_cast<deprecated_gsl_vector *>( w->knots );
}

size_t
deprecated_wrapper_gsl_vector_get_size( const deprecated_gsl_vector* v)
{
	return v->size;
}

double *
deprecated_wrapper_gsl_vector_get_data( const deprecated_gsl_vector* v)
{
	return v->data;
}

deprecated_gsl_multifit_function_fdf_ptr
deprecated_wrapper_gsl_multifit_fdfsolver_lmsder_new(
		int (* f) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f),
		int (* df) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_matrix * df),
		int (* fdf) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f, deprecated_gsl_matrix *df),
		size_t n,   /* number of functions */
		size_t p,   /* number of independent variables */
		void * params )
{
	deprecated_gsl_multifit_function_fdf* ptr
			= (deprecated_gsl_multifit_function_fdf*) malloc(sizeof(struct deprecated_gsl_multifit_function_fdf));
	deprecated_gsl_multifit_function_fdf_ptr a (
			ptr,
			std::ptr_fun( free ) );

	a->f = (int (*)(const gsl_vector* x, void* params, gsl_vector* f)) f;
	a->df = (int (*)(const gsl_vector* x, void* params, gsl_matrix* J)) df;
	a->fdf = (int (*)(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)) fdf;
	a->n = n;
	a->p = p;
	a->params = params;
	return a;
}

const char * deprecated_gsl_strerror (const int e)
{
	return gsl_strerror(e);
}


const char * deprecated_gsl_rng_name (const deprecated_gsl_rng * r)
{
	return gsl_rng_name (r);
}


const deprecated_gsl_rng_type * deprecated_gsl_rng_env_setup (void)
{
	return  static_cast<const deprecated_gsl_rng_type*>( gsl_rng_env_setup() );
}


double deprecated_gsl_blas_dnrm2 (const deprecated_gsl_vector * X)
{
	return gsl_blas_dnrm2 (X);
}


double deprecated_gsl_cdf_gaussian_P (const double x, const double sigma)
{
	return gsl_cdf_gaussian_P (x, sigma);
}


double deprecated_gsl_cdf_tdist_Pinv (const double P, const double nu)
{
	return gsl_cdf_tdist_Pinv (P, nu);
}


double deprecated_gsl_ran_cauchy (const deprecated_gsl_rng * r, const double a)
{
	return gsl_ran_cauchy (r, a);
}


double deprecated_gsl_ran_cauchy_pdf (const double x, const double a)
{
	return gsl_ran_cauchy_pdf (x, a);
}


double deprecated_gsl_ran_exponential (const deprecated_gsl_rng * r, const double mu)
{
	return gsl_ran_exponential (r, mu);
}


double deprecated_gsl_ran_exponential_pdf (const double x, const double mu)
{
	return gsl_ran_exponential_pdf (x, mu);
}


double deprecated_gsl_ran_flat (const deprecated_gsl_rng * r, const double a, const double b)
{
	return gsl_ran_flat (r, a, b);
}


double deprecated_gsl_ran_gaussian (const deprecated_gsl_rng * r, const double sigma)
{
	return gsl_ran_gaussian (r, sigma);
}


double deprecated_gsl_stats_absdev (const double data[], const size_t stride, const size_t n)
{
	return gsl_stats_absdev (data , stride, n);
}


double deprecated_gsl_stats_covariance (const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n)
{
	return gsl_stats_covariance (data1 , stride1, data2 , stride2, n);
}


double deprecated_gsl_stats_int_mean (const int data[], const size_t stride, const size_t n)
{
	return gsl_stats_int_mean (data , stride, n);
}


double deprecated_gsl_stats_int_sd (const int data[], const size_t stride, const size_t n)
{
	return gsl_stats_int_sd (data , stride, n);
}


double deprecated_gsl_stats_mean (const double data[], const size_t stride, const size_t n)
{
	return gsl_stats_mean (data , stride, n);
}


double deprecated_gsl_stats_median_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n)
{
	return gsl_stats_median_from_sorted_data (sorted_data , stride, n) ;
}

double deprecated_gsl_stats_quantile_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, const double f)
{
	return gsl_stats_quantile_from_sorted_data (sorted_data , stride, n, f) ;
}

double deprecated_gsl_stats_sd (const double data[], const size_t stride, const size_t n)
{
	return gsl_stats_sd (data , stride, n);
}


double deprecated_gsl_stats_sd_with_fixed_mean (const double data[], const size_t stride, const size_t n, const double mean)
{
	return gsl_stats_sd_with_fixed_mean (data , stride, n, mean);
}


double deprecated_gsl_stats_variance (const double data[], const size_t stride, const size_t n)
{
	return gsl_stats_variance (data , stride, n);
}


double deprecated_gsl_stats_variance_m (const double data[], const size_t stride, const size_t n, const double mean)
{
	return gsl_stats_variance_m (data , stride, n, mean);
}


deprecated_gsl_fft_real_wavetable * deprecated_gsl_fft_real_wavetable_alloc (size_t n)
{
	return static_cast<deprecated_gsl_fft_real_wavetable*>( gsl_fft_real_wavetable_alloc (n) );
}


deprecated_gsl_fft_real_workspace * deprecated_gsl_fft_real_workspace_alloc (size_t n)
{
	return static_cast<deprecated_gsl_fft_real_workspace*>( gsl_fft_real_workspace_alloc (n) );
}


deprecated_gsl_matrix *deprecated_gsl_matrix_alloc (const size_t n1, const size_t n2)
{
	return static_cast<deprecated_gsl_matrix*>( gsl_matrix_alloc(n1, n2) );
}


deprecated_gsl_matrix * deprecated_gsl_matrix_calloc (const size_t n1, const size_t n2)
{
	return static_cast<deprecated_gsl_matrix*>( gsl_matrix_calloc(n1, n2) );
}

double deprecated_gsl_matrix_get(const deprecated_gsl_matrix * m, const size_t i, const size_t j)
{
	return gsl_matrix_get(m, i, j);
}

void deprecated_gsl_multifit_fdfsolver_free (deprecated_gsl_multifit_fdfsolver * s)
{
	gsl_multifit_fdfsolver_free(s);
}


int deprecated_gsl_multifit_fdfsolver_iterate (deprecated_gsl_multifit_fdfsolver * s)
{
	return gsl_multifit_fdfsolver_iterate(s);
}


deprecated_gsl_multifit_linear_workspace * deprecated_gsl_multifit_linear_alloc (size_t n, size_t p)
{
	return static_cast<deprecated_gsl_multifit_linear_workspace*>( gsl_multifit_linear_alloc(n, p) );
}


void deprecated_gsl_multifit_linear_free (deprecated_gsl_multifit_linear_workspace * work)
{
	gsl_multifit_linear_free(work);
}


deprecated_gsl_permutation * deprecated_gsl_permutation_alloc (const size_t n)
{
	return static_cast<deprecated_gsl_permutation*>( gsl_permutation_alloc (n) );
}


deprecated_gsl_ran_discrete_t * deprecated_gsl_ran_discrete_preproc (size_t K, const double *P)
{
	return static_cast<deprecated_gsl_ran_discrete_t*>( gsl_ran_discrete_preproc (K, P) );
}


deprecated_gsl_rng * deprecated_gsl_rng_alloc (const deprecated_gsl_rng_type * T)
{
	return static_cast<deprecated_gsl_rng*>( gsl_rng_alloc (T) );
}


unsigned long int deprecated_gsl_rng_get (const deprecated_gsl_rng * r)
{
	return gsl_rng_get (r);
}


double deprecated_gsl_rng_uniform (const deprecated_gsl_rng * r)
{
	return gsl_rng_uniform (r);
}


deprecated_gsl_vector * deprecated_gsl_vector_alloc (const size_t n)
{
	return static_cast<deprecated_gsl_vector*>( gsl_vector_alloc (n) );
}


double deprecated_gsl_vector_get (const deprecated_gsl_vector * v, const size_t i)
{
	return gsl_vector_get (v, i);
}


void deprecated_gsl_vector_set (deprecated_gsl_vector * v, const size_t i, double x)
{
	gsl_vector_set (v, i, x);
}


deprecated_gsl_vector_view_ptr deprecated_gsl_vector_view_array (double *v, size_t n)
{
	deprecated_gsl_vector_view* ptr
				= (deprecated_gsl_vector_view*) malloc(sizeof(struct deprecated_gsl_vector_view));
	deprecated_gsl_vector_view_ptr b (
				ptr,
				std::ptr_fun( free ) );
	gsl_vector_view a = gsl_vector_view_array (v, n);
	b->vector = a.vector;

	return b;
}


int deprecated_gsl_linalg_LU_decomp (deprecated_gsl_matrix * A, deprecated_gsl_permutation * p, int *signum)
{
	return gsl_linalg_LU_decomp (A, p, signum);
}


int deprecated_gsl_matrix_add (deprecated_gsl_matrix * a, const deprecated_gsl_matrix * b)
{
	return gsl_matrix_add (a, b);
}


int deprecated_gsl_matrix_fprintf (FILE * stream, const deprecated_gsl_matrix * m, const char * format)
{
	return gsl_matrix_fprintf (stream, m, format);
}


int deprecated_gsl_matrix_scale (deprecated_gsl_matrix * a, const double x)
{
	return gsl_matrix_scale (a, x);
}


int deprecated_gsl_multifit_covar (const deprecated_gsl_matrix * J, double epsrel, deprecated_gsl_matrix * covar)
{
	return gsl_multifit_covar (J, epsrel, covar);
}


int deprecated_gsl_vector_memcpy (deprecated_gsl_vector * dest, const deprecated_gsl_vector * src)
{
	return gsl_vector_memcpy (dest, src);
}


size_t deprecated_gsl_ran_discrete (const deprecated_gsl_rng *r, const deprecated_gsl_ran_discrete_t *g)
{
	return gsl_ran_discrete (r, g);
}


unsigned int deprecated_gsl_ran_binomial (const deprecated_gsl_rng * r, double p, unsigned int n)
{
	return gsl_ran_binomial (r, p, n);
}


unsigned int deprecated_gsl_ran_poisson (const deprecated_gsl_rng * r, double mu)
{
	return gsl_ran_poisson (r, mu);
}


unsigned long int deprecated_gsl_rng_max (const deprecated_gsl_rng * r)
{
	return gsl_rng_max (r);
}


unsigned long int deprecated_gsl_rng_min (const deprecated_gsl_rng * r)
{
	return gsl_rng_min (r);
}


void deprecated_gsl_fft_real_wavetable_free (deprecated_gsl_fft_real_wavetable * wavetable)
{
	gsl_fft_real_wavetable_free (wavetable);
}


void deprecated_gsl_matrix_free (deprecated_gsl_matrix * m)
{
	gsl_matrix_free (m);
}


void deprecated_gsl_matrix_set_zero (deprecated_gsl_matrix * m)
{
	gsl_matrix_set_zero (m);
}

void deprecated_gsl_matrix_set(deprecated_gsl_matrix * m, const size_t i, const size_t j, const double x)
{
	gsl_matrix_set(m, i, j, x);
}

void deprecated_gsl_permutation_free (deprecated_gsl_permutation * p)
{
	gsl_permutation_free (p);
}

double deprecated_gsl_pow_int(double x, int n)
{
	return gsl_pow_int(x, n);
}


void deprecated_gsl_rng_free (deprecated_gsl_rng * r)
{
	gsl_rng_free (r);
}


void deprecated_gsl_rng_set (const deprecated_gsl_rng * r, unsigned long int seed)
{
	gsl_rng_set (r, seed);
}


void deprecated_gsl_sort (double * data, const size_t stride, const size_t n)
{
	gsl_sort (data, stride, n);
}


void deprecated_gsl_sort_index (size_t * p, const double * data, const size_t stride, const size_t n)
{
	gsl_sort_index (p, data, stride, n);
}


void deprecated_gsl_sort_vector (deprecated_gsl_vector * v)
{
	gsl_sort_vector (v);
}


void deprecated_gsl_vector_free (deprecated_gsl_vector * v)
{
	gsl_vector_free (v);
}


void deprecated_gsl_vector_minmax (const deprecated_gsl_vector * v, double * min_out, double * max_out)
{
	gsl_vector_minmax (v, min_out, max_out);
}


void deprecated_gsl_vector_set_zero (deprecated_gsl_vector * v)
{
	gsl_vector_set_zero (v);
}


int deprecated_gsl_linalg_SV_decomp (deprecated_gsl_matrix * A, deprecated_gsl_matrix * V, deprecated_gsl_vector * S, deprecated_gsl_vector * work)
{
	return gsl_linalg_SV_decomp(A, V, S, work);
}


deprecated_gsl_multifit_fdfsolver *
deprecated_gsl_multifit_fdfsolver_alloc (
		const deprecated_gsl_multifit_fdfsolver_type * T,
		size_t n,
		size_t p)
{
	return static_cast<deprecated_gsl_multifit_fdfsolver*>( gsl_multifit_fdfsolver_alloc(T, n, p) );
}


int
deprecated_gsl_multifit_fdfsolver_set (
		deprecated_gsl_multifit_fdfsolver * s,
		deprecated_gsl_multifit_function_fdf * fdf,
		const deprecated_gsl_vector * x)
{
	return gsl_multifit_fdfsolver_set(s, fdf, x);
}


int deprecated_gsl_multifit_linear (const deprecated_gsl_matrix * X, const deprecated_gsl_vector * y, deprecated_gsl_vector * c, deprecated_gsl_matrix * cov, double * chisq, deprecated_gsl_multifit_linear_workspace * work)
{
	return gsl_multifit_linear(X, y, c, cov, chisq, work);
}


int deprecated_gsl_multifit_linear_est (const deprecated_gsl_vector * x, const deprecated_gsl_vector * c, const deprecated_gsl_matrix * cov, double *y, double *y_err)
{
	return gsl_multifit_linear_est(x, c, cov, y, y_err);
}


int deprecated_gsl_multifit_wlinear (const deprecated_gsl_matrix * X, const deprecated_gsl_vector * w, const deprecated_gsl_vector * y, deprecated_gsl_vector * c, deprecated_gsl_matrix * cov, double * chisq, deprecated_gsl_multifit_linear_workspace * work)
{
	return gsl_multifit_wlinear(X, w, y, c, cov, chisq, work);
}


int deprecated_gsl_blas_dgemm (deprecated_gsl_CBLAS_TRANSPOSE_t TransA, deprecated_gsl_CBLAS_TRANSPOSE_t TransB, double alpha, const deprecated_gsl_matrix* A, const deprecated_gsl_matrix* B, double beta, deprecated_gsl_matrix* C)
{
	return gsl_blas_dgemm ((CBLAS_TRANSPOSE_t)TransA, (CBLAS_TRANSPOSE_t)TransB, alpha, A, B, beta, C);
}


int deprecated_gsl_fft_real_transform (double data[], const size_t stride, const size_t n, const deprecated_gsl_fft_real_wavetable * wavetable, deprecated_gsl_fft_real_workspace * work)
{
	return gsl_fft_real_transform (data , stride, n, wavetable, work);
}


int deprecated_gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, const size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq)
{
	return gsl_fit_linear (x, xstride, y, ystride, n, c0, c1, cov00, cov01, cov11, sumsq);
}


int deprecated_gsl_fit_mul (const double * x, const size_t xstride, const double * y, const size_t ystride, const size_t n, double * c1, double * cov11, double * sumsq)
{
	return gsl_fit_mul (x, xstride, y, ystride, n, c1, cov11, sumsq);
}


int deprecated_gsl_fit_wlinear (const double * x, const size_t xstride, const double * w, const size_t wstride, const double * y, const size_t ystride, const size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * chisq)
{
	return gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, n, c0, c1, cov00, cov01, cov11, chisq);
}


int deprecated_gsl_linalg_LU_solve (const deprecated_gsl_matrix * LU, const deprecated_gsl_permutation * p, const deprecated_gsl_vector * b, deprecated_gsl_vector * x)
{
	return gsl_linalg_LU_solve (LU, p, b, x);
}


int deprecated_gsl_multifit_test_delta (const deprecated_gsl_vector * dx, const deprecated_gsl_vector * x, double epsabs, double epsrel)
{
	return gsl_multifit_test_delta (dx, x, epsabs, epsrel);
}

double deprecated_gsl_spline_eval_deriv(const deprecated_gsl_spline * spline, double x, deprecated_gsl_interp_accel * a)
{
	return gsl_spline_eval_deriv(spline, x, a);
}

double deprecated_gsl_spline_eval(const deprecated_gsl_spline * spline, double x, deprecated_gsl_interp_accel * a)
{
	return gsl_spline_eval(spline, x, a);
}

void deprecated_gsl_spline_free(deprecated_gsl_spline * spline)
{
	gsl_spline_free(spline);
}
deprecated_gsl_interp_accel * deprecated_gsl_interp_accel_alloc(void)
{
	return static_cast<deprecated_gsl_interp_accel*>( gsl_interp_accel_alloc() );
}

void deprecated_gsl_interp_accel_free(deprecated_gsl_interp_accel * a)
{
	gsl_interp_accel_free(a);
}
deprecated_gsl_spline * deprecated_gsl_spline_alloc(const deprecated_gsl_interp_type * T, size_t size)
{
	return static_cast<deprecated_gsl_spline*>( gsl_spline_alloc(T, size) );
}

int deprecated_gsl_spline_init(deprecated_gsl_spline * spline, const double xa[], const double ya[], size_t size)
{
	return gsl_spline_init(spline, xa, ya, size);
}

deprecated_gsl_bspline_workspace* deprecated_gsl_bspline_alloc(const size_t k, const size_t nbreak)
{
  return static_cast<deprecated_gsl_bspline_workspace*>( gsl_bspline_alloc(k, nbreak) );
}

deprecated_gsl_bspline_deriv_workspace* deprecated_gsl_bspline_deriv_alloc(const size_t k)
{
  return  static_cast<deprecated_gsl_bspline_deriv_workspace*>( gsl_bspline_deriv_alloc(k) );
}

int deprecated_gsl_bspline_deriv_eval(const double x, const size_t nderiv, deprecated_gsl_matrix* dB, deprecated_gsl_bspline_workspace* w, deprecated_gsl_bspline_deriv_workspace * dw)
{
  return gsl_bspline_deriv_eval(x, nderiv, dB, w, dw);
}

void deprecated_gsl_bspline_deriv_free(deprecated_gsl_bspline_deriv_workspace* w)
{
  gsl_bspline_deriv_free(w);
}

int deprecated_gsl_bspline_eval(const double x, deprecated_gsl_vector* B, deprecated_gsl_bspline_workspace* w)
{
  return gsl_bspline_eval(x, B, w);
}

void deprecated_gsl_bspline_free(deprecated_gsl_bspline_workspace* w)
{
  gsl_bspline_free(w);
}

int deprecated_gsl_bspline_knots(const deprecated_gsl_vector* breakpts, deprecated_gsl_bspline_workspace * w)
{
  return gsl_bspline_knots(breakpts, w);
}

int deprecated_gsl_bspline_knots_uniform(const double a, const double b, deprecated_gsl_bspline_workspace* w)
{
  return gsl_bspline_knots_uniform(a, b, w);
}

size_t deprecated_gsl_bspline_ncoeffs(deprecated_gsl_bspline_workspace* w)
{
  return  gsl_bspline_ncoeffs(w);
}

deprecated_gsl_interp* deprecated_gsl_interp_alloc(const deprecated_gsl_interp_type* T, size_t n)
{
  return static_cast<deprecated_gsl_interp*>( gsl_interp_alloc(T, n) );
}

double deprecated_gsl_interp_eval(const deprecated_gsl_interp* obj, const double xa[], const double ya[], double x, deprecated_gsl_interp_accel* a)
{
  return gsl_interp_eval(obj, xa, ya, x, a);
}

void deprecated_gsl_interp_free(deprecated_gsl_interp* interp)
{
  gsl_interp_free(interp);
}

int deprecated_gsl_interp_init(deprecated_gsl_interp* obj, const double xa[], const double ya[], size_t size)
{
  return gsl_interp_init(obj, xa, ya, size);
}

void deprecated_gsl_ran_discrete_free(deprecated_gsl_ran_discrete_t *g)
{
	gsl_ran_discrete_free(g);
}

void  deprecated_gsl_fft_real_workspace_free (deprecated_gsl_fft_real_workspace * workspace)
{
	gsl_fft_real_workspace_free (workspace);
}

double  deprecated_gsl_sf_psi(const double x)
{
	return gsl_sf_psi(x);
}
