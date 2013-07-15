/*
* deprecated_gsl_wrapper.h
*
* Created on: Jul 9, 2013
* Author: Hans-Christian Ehrlich
*/

#ifndef GSL_WRAPPER_H_
#define GSL_WRAPPER_H_
#include <stdio.h>

#include <boost/shared_ptr.hpp>

struct deprecated_gsl_bspline_deriv_workspace;
struct deprecated_gsl_bspline_deriv_workspace;
struct deprecated_gsl_bspline_workspace;
struct deprecated_gsl_fft_real_wavetable;
struct deprecated_gsl_fft_real_workspace;
struct deprecated_gsl_interp;
struct deprecated_gsl_interp_accel;
struct deprecated_gsl_interp_type;
struct deprecated_gsl_matrix;
struct deprecated_gsl_multifit_function_fdf;
struct deprecated_gsl_multifit_fdfsolver;
struct deprecated_gsl_multifit_fdfsolver_type;
struct deprecated_gsl_multifit_linear_workspace;
struct deprecated_gsl_permutation;
struct deprecated_gsl_ran_discrete_t;
struct deprecated_gsl_rng;
struct deprecated_gsl_rng_type;
struct deprecated_gsl_spline;
struct deprecated_gsl_vector;
struct deprecated_gsl_vector_view;

typedef boost::shared_ptr<deprecated_gsl_multifit_function_fdf> deprecated_gsl_multifit_function_fdf_ptr;
typedef boost::shared_ptr<deprecated_gsl_vector_view> deprecated_gsl_vector_view_ptr;

#define deprecated_gsl_interp_linear gsl_interp_linear
#define deprecated_gsl_interp_polynomial gsl_interp_polynomial
#define deprecated_gsl_interp_cspline gsl_interp_cspline
#define deprecated_gsl_interp_akima gsl_interp_akima

/* copy and past from gsl (not needed parts removed)*/
enum {
  deprecated_gsl_SUCCESS  = 0,
  deprecated_gsl_FAILURE  = -1,
  deprecated_gsl_CONTINUE = -2,  /* iteration has not converged */
  deprecated_gsl_EDOM     = 1   /* input domain error, e.g sqrt(-1) */
} ;
enum deprecated_gsl_CBLAS_TRANSPOSE {
	deprecated_gsl_CblasNoTrans=111
};
typedef  enum deprecated_gsl_CBLAS_TRANSPOSE deprecated_gsl_CBLAS_TRANSPOSE_t;

void
deprecated_wrapper_gsl_rng_default_seed_set( unsigned long int x );

const deprecated_gsl_rng_type *
deprecated_wrapper_gsl_rng_taus_get();

const deprecated_gsl_rng_type *
deprecated_wrapper_get_gsl_rng_mt19937();

const deprecated_gsl_multifit_fdfsolver_type *
deprecated_wrapper_get_multifit_fdfsolver_lmsder();

const deprecated_gsl_rng_type *
deprecated_wrapper_get_gsl_rng_default();

const deprecated_gsl_vector*
deprecated_wrapper_gsl_vector_view_get_vector( deprecated_gsl_vector_view_ptr v );

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_x( deprecated_gsl_multifit_fdfsolver* s );

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_f( deprecated_gsl_multifit_fdfsolver* s );

deprecated_gsl_vector*
deprecated_wrapper_gsl_multifit_fdfsolver_get_dx( deprecated_gsl_multifit_fdfsolver* s );

deprecated_gsl_matrix*
deprecated_wrapper_gsl_multifit_fdfsolver_get_J( deprecated_gsl_multifit_fdfsolver* s );

size_t
deprecated_wrapper_gsl_matrix_get_size1( deprecated_gsl_matrix* M);

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_linear();

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_polynomial();

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_cspline();

const deprecated_gsl_interp_type *
deprecated_wrapper_get_gsl_interp_akima();

size_t
deprecated_wrapper_gsl_interp_type_get_min_size( const deprecated_gsl_interp_type * t );

deprecated_gsl_vector *
deprecated_wrapper_gsl_bspline_workspace_get_knots( deprecated_gsl_bspline_workspace * w );

size_t
deprecated_wrapper_gsl_vector_get_size( const deprecated_gsl_vector* v);

double *
deprecated_wrapper_gsl_vector_get_data( const deprecated_gsl_vector* v);


deprecated_gsl_multifit_function_fdf_ptr
deprecated_wrapper_gsl_multifit_fdfsolver_lmsder_new(
		int (* f) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f),
		int (* df) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_matrix * df),
		int (* fdf) (const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f, deprecated_gsl_matrix *df),
		size_t n,   /* number of functions */
		size_t p,   /* number of independent variables */
		void * params );

int deprecated_gsl_blas_dgemm (
		deprecated_gsl_CBLAS_TRANSPOSE_t TransA,
		deprecated_gsl_CBLAS_TRANSPOSE_t TransB,
		double alpha,
		const deprecated_gsl_matrix* A,
		const deprecated_gsl_matrix* B,
		double beta,
		deprecated_gsl_matrix* C);

double deprecated_gsl_blas_dnrm2 (const deprecated_gsl_vector* X);

double deprecated_gsl_cdf_gaussian_P (const double x, const double sigma);
double deprecated_gsl_cdf_tdist_Pinv (const double P, const double nu);

int deprecated_gsl_fft_real_transform (double data[], const size_t stride, const size_t n, const deprecated_gsl_fft_real_wavetable* wavetable, deprecated_gsl_fft_real_workspace* work);

deprecated_gsl_fft_real_wavetable* deprecated_gsl_fft_real_wavetable_alloc (size_t n);
void deprecated_gsl_fft_real_wavetable_free (deprecated_gsl_fft_real_wavetable* wavetable);
deprecated_gsl_fft_real_workspace* deprecated_gsl_fft_real_workspace_alloc (size_t n);

int deprecated_gsl_fit_linear (const double* x, const size_t xstride, const double* y, const size_t ystride, const size_t n, double* c0, double* c1, double* cov00, double* cov01, double* cov11, double* sumsq);
int deprecated_gsl_fit_mul (const double* x, const size_t xstride, const double* y, const size_t ystride, const size_t n, double* c1, double* cov11, double* sumsq);
int deprecated_gsl_fit_wlinear (const double* x, const size_t xstride, const double* w, const size_t wstride, const double* y, const size_t ystride, const size_t n, double* c0, double* c1, double* cov00, double* cov01, double* cov11, double* chisq);

void deprecated_gsl_interp_accel_free(deprecated_gsl_interp_accel* a);

int deprecated_gsl_linalg_LU_decomp (deprecated_gsl_matrix* A, deprecated_gsl_permutation* p, int* signum);
int deprecated_gsl_linalg_LU_solve (const deprecated_gsl_matrix* LU, const deprecated_gsl_permutation* p, const deprecated_gsl_vector* b, deprecated_gsl_vector* x);
int deprecated_gsl_linalg_SV_decomp (deprecated_gsl_matrix* A, deprecated_gsl_matrix* V, deprecated_gsl_vector* S, deprecated_gsl_vector* work);

int deprecated_gsl_matrix_add (deprecated_gsl_matrix* a, const deprecated_gsl_matrix* b);
deprecated_gsl_matrix* deprecated_gsl_matrix_alloc (const size_t n1, const size_t n2);
deprecated_gsl_matrix* deprecated_gsl_matrix_calloc (const size_t n1, const size_t n2);
int deprecated_gsl_matrix_fprintf (FILE* stream, const deprecated_gsl_matrix* m, const char* format);
void deprecated_gsl_matrix_free (deprecated_gsl_matrix* m);
double deprecated_gsl_matrix_get(const deprecated_gsl_matrix* m, const size_t i, const size_t j);
int deprecated_gsl_matrix_scale (deprecated_gsl_matrix* a, const double x);
void deprecated_gsl_matrix_set(deprecated_gsl_matrix* m, const size_t i, const size_t j, const double x);
void deprecated_gsl_matrix_set_zero (deprecated_gsl_matrix* m);

int
deprecated_gsl_multifit_covar (
		const deprecated_gsl_matrix* J,
		double epsrel,
		deprecated_gsl_matrix* covar);

deprecated_gsl_multifit_fdfsolver*
deprecated_gsl_multifit_fdfsolver_alloc (
		const deprecated_gsl_multifit_fdfsolver_type* T,
		size_t n,
		size_t p);

void
deprecated_gsl_multifit_fdfsolver_free (deprecated_gsl_multifit_fdfsolver* s);

int
deprecated_gsl_multifit_fdfsolver_iterate (deprecated_gsl_multifit_fdfsolver* s);

int
deprecated_gsl_multifit_fdfsolver_set (
		deprecated_gsl_multifit_fdfsolver* s,
		deprecated_gsl_multifit_function_fdf* fdf,
		const deprecated_gsl_vector* x);

int
deprecated_gsl_multifit_linear (
		const deprecated_gsl_matrix* X,
		const deprecated_gsl_vector* y,
		deprecated_gsl_vector* c,
		deprecated_gsl_matrix* cov,
		double* chisq,
		deprecated_gsl_multifit_linear_workspace* work);

deprecated_gsl_multifit_linear_workspace* deprecated_gsl_multifit_linear_alloc (size_t n, size_t p);
int deprecated_gsl_multifit_linear_est (const deprecated_gsl_vector* x, const deprecated_gsl_vector* c, const deprecated_gsl_matrix* cov, double* y, double* y_err);
void deprecated_gsl_multifit_linear_free (deprecated_gsl_multifit_linear_workspace* work);
int deprecated_gsl_multifit_test_delta (const deprecated_gsl_vector* dx, const deprecated_gsl_vector* x, double epsabs, double epsrel);
int deprecated_gsl_multifit_wlinear (const deprecated_gsl_matrix* X, const deprecated_gsl_vector* w, const deprecated_gsl_vector* y, deprecated_gsl_vector* c, deprecated_gsl_matrix* cov, double* chisq, deprecated_gsl_multifit_linear_workspace* work);

deprecated_gsl_permutation* deprecated_gsl_permutation_alloc (const size_t n);
void deprecated_gsl_permutation_free (deprecated_gsl_permutation* p);

double deprecated_gsl_pow_int(double x, int n);

unsigned int deprecated_gsl_ran_binomial (const deprecated_gsl_rng* r, double p, unsigned int n);
double deprecated_gsl_ran_cauchy (const deprecated_gsl_rng* r, const double a);
double deprecated_gsl_ran_cauchy_pdf (const double x, const double a);
size_t deprecated_gsl_ran_discrete (const deprecated_gsl_rng* r, const deprecated_gsl_ran_discrete_t* g);
deprecated_gsl_ran_discrete_t* deprecated_gsl_ran_discrete_preproc (size_t K, const double* P);
double deprecated_gsl_ran_exponential (const deprecated_gsl_rng* r, const double mu);
double deprecated_gsl_ran_exponential_pdf (const double x, const double mu);
double deprecated_gsl_ran_flat (const deprecated_gsl_rng* r, const double a, const double b);
double deprecated_gsl_ran_gaussian (const deprecated_gsl_rng* r, const double sigma);
unsigned int deprecated_gsl_ran_poisson (const deprecated_gsl_rng* r, double mu);

deprecated_gsl_rng* deprecated_gsl_rng_alloc (const deprecated_gsl_rng_type* T);
const deprecated_gsl_rng_type* deprecated_gsl_rng_env_setup (void);
void deprecated_gsl_rng_free (deprecated_gsl_rng* r);
unsigned long int deprecated_gsl_rng_get (const deprecated_gsl_rng* r);
unsigned long int deprecated_gsl_rng_max (const deprecated_gsl_rng* r);
unsigned long int deprecated_gsl_rng_min (const deprecated_gsl_rng* r);
const char* deprecated_gsl_rng_name (const deprecated_gsl_rng* r);
void deprecated_gsl_rng_set (const deprecated_gsl_rng* r, unsigned long int seed);
double deprecated_gsl_rng_uniform (const deprecated_gsl_rng* r);

void deprecated_gsl_sort (double* data, const size_t stride, const size_t n);
void deprecated_gsl_sort_index (size_t* p, const double* data, const size_t stride, const size_t n);
void deprecated_gsl_sort_vector (deprecated_gsl_vector* v);

double deprecated_gsl_spline_eval(const deprecated_gsl_spline* spline, double x, deprecated_gsl_interp_accel* a);
double deprecated_gsl_spline_eval_deriv(const deprecated_gsl_spline* spline, double x, deprecated_gsl_interp_accel* a);
void deprecated_gsl_spline_free(deprecated_gsl_spline* spline);
int deprecated_gsl_spline_init(deprecated_gsl_spline* spline, const double xa[], const double ya[], size_t size);

double deprecated_gsl_stats_absdev (const double data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_covariance (const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n);
double deprecated_gsl_stats_int_mean (const int data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_int_sd (const int data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_mean (const double data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_median_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_quantile_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, const double f);
double deprecated_gsl_stats_sd (const double data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_sd_with_fixed_mean (const double data[], const size_t stride, const size_t n, const double mean);
double deprecated_gsl_stats_variance (const double data[], const size_t stride, const size_t n);
double deprecated_gsl_stats_variance_m (const double data[], const size_t stride, const size_t n, const double mean);

const char* deprecated_gsl_strerror (const int e);

deprecated_gsl_vector* deprecated_gsl_vector_alloc (const size_t n);
void deprecated_gsl_vector_free (deprecated_gsl_vector* v);
double deprecated_gsl_vector_get (const deprecated_gsl_vector* v, const size_t i);
int deprecated_gsl_vector_memcpy (deprecated_gsl_vector* dest, const deprecated_gsl_vector* src);
void deprecated_gsl_vector_minmax (const deprecated_gsl_vector* v, double* min_out, double* max_out);
void deprecated_gsl_vector_set (deprecated_gsl_vector* v, const size_t i, double x);
void deprecated_gsl_vector_set_zero (deprecated_gsl_vector* v);
deprecated_gsl_vector_view_ptr deprecated_gsl_vector_view_array (double* v, size_t n);



deprecated_gsl_bspline_workspace* deprecated_gsl_bspline_alloc(const size_t k, const size_t nbreak);

deprecated_gsl_bspline_deriv_workspace* deprecated_gsl_bspline_deriv_alloc(const size_t k);

int deprecated_gsl_bspline_deriv_eval(const double x, const size_t nderiv, deprecated_gsl_matrix* dB, deprecated_gsl_bspline_workspace* w, deprecated_gsl_bspline_deriv_workspace * dw);

void deprecated_gsl_bspline_deriv_free(deprecated_gsl_bspline_deriv_workspace* w);

int deprecated_gsl_bspline_eval(const double x, deprecated_gsl_vector* B, deprecated_gsl_bspline_workspace* w);

void deprecated_gsl_bspline_free(deprecated_gsl_bspline_workspace* w);

int deprecated_gsl_bspline_knots(const deprecated_gsl_vector* breakpts, deprecated_gsl_bspline_workspace * w);

int deprecated_gsl_bspline_knots_uniform(const double a, const double b, deprecated_gsl_bspline_workspace* w);

size_t deprecated_gsl_bspline_ncoeffs(deprecated_gsl_bspline_workspace* w);

deprecated_gsl_interp_accel* deprecated_gsl_interp_accel_alloc(void);

deprecated_gsl_interp* deprecated_gsl_interp_alloc(const deprecated_gsl_interp_type* T, size_t n);

double deprecated_gsl_interp_eval(const deprecated_gsl_interp* obj, const double xa[], const double ya[], double x, deprecated_gsl_interp_accel* a);

void deprecated_gsl_interp_free(deprecated_gsl_interp* interp);

int deprecated_gsl_interp_init(deprecated_gsl_interp* obj, const double xa[], const double ya[], size_t size);

deprecated_gsl_spline* deprecated_gsl_spline_alloc(const deprecated_gsl_interp_type* T, size_t size);

void deprecated_gsl_ran_discrete_free(deprecated_gsl_ran_discrete_t *g);

void  deprecated_gsl_fft_real_workspace_free (deprecated_gsl_fft_real_workspace * workspace);

double  deprecated_gsl_sf_psi(const double x);

#endif /* deprecated_gsl_WRAPPER_H_*/
