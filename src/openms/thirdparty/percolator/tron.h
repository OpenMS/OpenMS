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

#ifndef _TRON_H
#define _TRON_H

namespace OpenMS { namespace Internal { namespace Percolator {

class function
{
public:
	virtual double fun(double *w) = 0 ;
	virtual void grad(double *w, double *g) = 0 ;
	virtual void Hv(double *s, double *Hs) = 0 ;

	virtual int get_nr_variable(void) = 0 ;
	virtual ~function(void){}
};

class TRON
{
public:
	TRON(const function *fun_obj, double eps = 0.1, double eps_cg = 0.1, int max_iter = 1000);
	~TRON();

	void tron(double *w);
	void set_print_string(void (*i_print) (const char *buf));

private:
	int trcg(double delta, double *g, double *s, double *r, bool *reach_boundary);
	double norm_inf(int n, double *x);

	double eps;
	double eps_cg;
	int max_iter;
	function *fun_obj;
	void info(const char *fmt,...);
	void (*tron_print_string)(const char *buf);
};

}}}  // namespace OpenMS::Internal::Percolator

#endif
