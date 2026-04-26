/*
 * TRON: Trust Region Newton method for L2-loss / logistic regression.
 *
 * This file is taken (without source-level modifications) from
 * John Halloran's Percolator-TRON integration at
 * https://bitbucket.org/jthalloran/percolator_upgrade. Halloran's copy
 * is itself a near-verbatim copy of liblinear v2.11 (downloaded
 * April 2017 per Halloran & Rocke 2018) with the upstream BSD-3
 * file header stripped.
 *
 * Copyright (c) 2007-2015 The LIBLINEAR Project. All rights reserved.
 * Distributed under the BSD 3-Clause License — see ../license.txt
 * for the full text.
 *
 * Algorithmic references:
 *   Lin, Weng, Keerthi (2008), "Trust Region Newton Method for
 *     Large-Scale Logistic Regression", JMLR 9:627-650.
 *   Halloran & Rocke (2018), "A Matter of Time: Faster Percolator
 *     Analysis...", J. Proteome Res. 17(5):1978-1982.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <Eigen/Dense>
#include "tron.h"

#ifndef min
template <class T> static inline T min(T x,T y) { return (x<y)?x:y; }
#endif

#ifndef max
template <class T> static inline T max(T x,T y) { return (x>y)?x:y; }
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// Eigen-backed replacements for the four BLAS-1 calls TRON used to make
// against extern "C" {dnrm2_, ddot_, daxpy_, dscal_}. The signatures are
// kept simple (n + raw pointers) so the call sites below stay as close as
// possible to the upstream liblinear v2.11 source. Eigen::Map adds zero
// allocations / copies; under -O3 the body inlines to the same x86_64 SIMD
// the system BLAS would emit, but without the .so dependency.
namespace {

inline double tron_nrm2(int n, const double* x)
{
  return Eigen::Map<const Eigen::VectorXd>(x, n).norm();
}

inline double tron_dot(int n, const double* x, const double* y)
{
  return Eigen::Map<const Eigen::VectorXd>(x, n).dot(
         Eigen::Map<const Eigen::VectorXd>(y, n));
}

inline void tron_axpy(int n, double alpha, const double* x, double* y)
{
  Eigen::Map<Eigen::VectorXd>(y, n) +=
      alpha * Eigen::Map<const Eigen::VectorXd>(x, n);
}

inline void tron_scal(int n, double alpha, double* x)
{
  Eigen::Map<Eigen::VectorXd>(x, n) *= alpha;
}

}  // anonymous namespace

static void default_print(const char *buf)
{
	fputs(buf,stdout);
	fflush(stdout);
}

void TRON::info(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	(*tron_print_string)(buf);
}

TRON::TRON(const function *fun_obj, double eps, double eps_cg, int max_iter)
{
	this->fun_obj=const_cast<function *>(fun_obj);
	this->eps=eps;
	this->eps_cg=eps_cg;
	this->max_iter=max_iter;
	tron_print_string = default_print;
}

TRON::~TRON()
{
}

void TRON::tron(double *w)
{
	// Parameters for updating the iterates.
	double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;

	// Parameters for updating the trust region size delta.
	double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

	int n = fun_obj->get_nr_variable();
	int i, cg_iter;
	double delta, snorm, one=1.0;
	double alpha, f, fnew, prered, actred, gs;
	int search = 1, iter = 1, inc = 1;
	double *s = new double[n];
	double *r = new double[n];
	double *g = new double[n];

	// calculate gradient norm at w=0 for stopping condition.
	double *w0 = new double[n];
	for (i=0; i<n; i++)
		w0[i] = 0;
	fun_obj->fun(w0);
	fun_obj->grad(w0, g);
	double gnorm0 = tron_nrm2(n, g);
	delete [] w0;

	f = fun_obj->fun(w);
	fun_obj->grad(w, g);
	delta = tron_nrm2(n, g);
	double gnorm = delta;

	if (gnorm <= eps*gnorm0)
		search = 0;

	iter = 1;

	double *w_new = new double[n];
	bool reach_boundary;
	while (iter <= max_iter && search)
	{
		cg_iter = trcg(delta, g, s, r, &reach_boundary);

		memcpy(w_new, w, sizeof(double)*n);
		tron_axpy(n, one, s, w_new);

		gs = tron_dot(n, g, s);
		prered = -0.5*(gs-tron_dot(n, s, r));
		fnew = fun_obj->fun(w_new);

		// Compute the actual reduction.
		actred = f - fnew;

		// On the first iteration, adjust the initial step bound.
		snorm = tron_nrm2(n, s);
		if (iter == 1)
			delta = min(delta, snorm);

		// Compute prediction alpha*snorm of the step.
		if (fnew - f - gs <= 0)
			alpha = sigma3;
		else
			alpha = max(sigma1, -0.5*(gs/(fnew - f - gs)));

		// Update the trust region bound according to the ratio of actual to predicted reduction.
		if (actred < eta0*prered)
			delta = min(max(alpha, sigma1)*snorm, sigma2*delta);
		else if (actred < eta1*prered)
			delta = max(sigma1*delta, min(alpha*snorm, sigma2*delta));
		else if (actred < eta2*prered)
			delta = max(sigma1*delta, min(alpha*snorm, sigma3*delta));
		else
		{
			if (reach_boundary)
				delta = sigma3*delta;
			else
				delta = max(delta, min(alpha*snorm, sigma3*delta));
		}

		info("iter %2d act %5.3e pre %5.3e delta %5.3e f %5.3e |g| %5.3e CG %3d\n", iter, actred, prered, delta, f, gnorm, cg_iter);

		if (actred > eta0*prered)
		{
			iter++;
			memcpy(w, w_new, sizeof(double)*n);
			f = fnew;
			fun_obj->grad(w, g);

			gnorm = tron_nrm2(n, g);
			if (gnorm <= eps*gnorm0)
				break;
		}
		if (f < -1.0e+32)
		{
			info("WARNING: f < -1.0e+32\n");
			break;
		}
		if (prered <= 0)
		{
			info("WARNING: prered <= 0\n");
			break;
		}
		if (fabs(actred) <= 1.0e-12*fabs(f) &&
		    fabs(prered) <= 1.0e-12*fabs(f))
		{
			info("WARNING: actred and prered too small\n");
			break;
		}
	}

	delete[] g;
	delete[] r;
	delete[] w_new;
	delete[] s;
}

int TRON::trcg(double delta, double *g, double *s, double *r, bool *reach_boundary)
{
	int i, inc = 1;
	int n = fun_obj->get_nr_variable();
	double one = 1;
	double *d = new double[n];
	double *Hd = new double[n];
	double rTr, rnewTrnew, alpha, beta, cgtol;

	*reach_boundary = false;
	for (i=0; i<n; i++)
	{
		s[i] = 0;
		r[i] = -g[i];
		d[i] = r[i];
	}
	cgtol = eps_cg*tron_nrm2(n, g);

	int cg_iter = 0;
	rTr = tron_dot(n, r, r);
	while (1)
	{
		if (tron_nrm2(n, r) <= cgtol)
			break;
		cg_iter++;
		fun_obj->Hv(d, Hd);

		alpha = rTr/tron_dot(n, d, Hd);
		tron_axpy(n, alpha, d, s);
		if (tron_nrm2(n, s) > delta)
		{
			info("cg reaches trust region boundary\n");
			*reach_boundary = true;
			alpha = -alpha;
			tron_axpy(n, alpha, d, s);

			double std = tron_dot(n, s, d);
			double sts = tron_dot(n, s, s);
			double dtd = tron_dot(n, d, d);
			double dsq = delta*delta;
			double rad = sqrt(std*std + dtd*(dsq-sts));
			if (std >= 0)
				alpha = (dsq - sts)/(std + rad);
			else
				alpha = (rad - std)/dtd;
			tron_axpy(n, alpha, d, s);
			alpha = -alpha;
			tron_axpy(n, alpha, Hd, r);
			break;
		}
		alpha = -alpha;
		tron_axpy(n, alpha, Hd, r);
		rnewTrnew = tron_dot(n, r, r);
		beta = rnewTrnew/rTr;
		tron_scal(n, beta, d);
		tron_axpy(n, one, r, d);
		rTr = rnewTrnew;
	}

	delete[] d;
	delete[] Hd;

	return(cg_iter);
}

double TRON::norm_inf(int n, double *x)
{
	double dmax = fabs(x[0]);
	for (int i=1; i<n; i++)
		if (fabs(x[i]) >= dmax)
			dmax = fabs(x[i]);
	return(dmax);
}

void TRON::set_print_string(void (*print_string) (const char *buf))
{
	tron_print_string = print_string;
}

}}}  // namespace OpenMS::Internal::Percolator
