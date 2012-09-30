// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENSWATH_ALGO__stats_COR__H__
#define OPENSWATH_ALGO__stats_COR__H__

/*! \file Corr.h cor {stats} R Documentation

\brief Correlation, Variance and Covariance (Matrices)


Description

var, cov and cor compute the variance of x and the covariance or correlation of x and y if these are vectors. If x and y are matrices then the covariances (or correlations) between the columns of x and the columns of y are computed.

cov2cor scales a covariance matrix into the corresponding correlation matrix efficiently.
Usage

- var(x, y = NULL, na.rm = FALSE, use)

- cov(x, y = NULL, use = "everything",
method = c("pearson", "kendall", "spearman"))

- cor(x, y = NULL, use = "everything",
method = c("pearson", "kendall", "spearman"))

cov2cor(V)

Arguments
x 	a numeric vector, matrix or data frame.
y 	NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).
na.rm 	logical. Should missing values be removed?
use 	an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
method 	a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
V 	symmetric numeric matrix, usually positive definite such as a covariance matrix.
Details

For cov and cor one must either give a matrix or data frame for x or give both x and y.

var is just another interface to cov, where na.rm is used to determine the default for use when that is unspecified. If na.rm is TRUE then the complete observations (rows) are used (use = "na.or.complete") to compute the variance. Otherwise, by default use = "everything".

If use is "everything", NAs will propagate conceptually, i.e., a resulting value will be NA whenever one of its contributing observations is NA.
If use is "all.obs", then the presence of missing observations will produce an error. If use is "complete.obs" then missing values are handled by casewise deletion (and if there are no complete cases, that gives an error).
"na.or.complete" is the same unless there are no complete cases, that gives NA. Finally, if use has the value "pairwise.complete.obs" then the correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables. This can result in covariance or correlation matrices which are not positive semi-definite, as well as NA entries if there are no complete pairs for that pair of variables. For cov and var, "pairwise.complete.obs" only works with the "pearson" method. Note that (the equivalent of) var(double(0), use=*) gives NA for use = "everything" and "na.or.complete", and gives an error in the other cases.

The denominator n - 1 is used which gives an unbiased estimator of the (co)variance for i.i.d. observations. These functions return NA when there is only one observation (whereas S-PLUS has been returning NaN), and fail if x has length zero.

For cor(), if method is "kendall" or "spearman", Kendall's tau or Spearman's rho statistic is used to estimate a rank-based measure of association. These are more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.
For cov(), a non-Pearson method is unusual but available for the sake of completeness. Note that "spearman" basically computes cor(R(x), R(y)) (or cov(.,.)) where R(u) := rank(u, na.last="keep"). In the case of missing values, the ranks are calculated depending on the value of use, either based on complete observations, or based on pairwise completeness with reranking for each pair.

Scaling a covariance matrix into a correlation one can be achieved in many ways, mathematically most appealing by multiplication with a diagonal matrix from left and right, or more efficiently by using sweep(.., FUN = "/") twice. The cov2cor function is even a bit more efficient, and provided mostly for didactical reasons.
Value

For r <- cor(*, use = "all.obs"), it is now guaranteed that all(r <= 1).
References

Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
See Also

cor.test for confidence intervals (and tests).

cov.wt for weighted covariance computation.

sd for standard deviation (vectors).

*/


//#include "StatsInclude/Base/Constants.h"
//#include "StatsInclude/Utils/Functors.h"
using namespace std;

#include <math.h>
#include <iterator>
/*! \file Stats/runmed.h

\brief Compute running medians of odd span. This is the �most robust� scatter plot smoothing possible.

*/

namespace OpenSwath
{

				template<typename TInputIterator, typename TInputIteratorY>
				typename std::iterator_traits<TInputIterator>::value_type cor_pearson(
					TInputIterator xBeg,
					TInputIterator xEnd,
					TInputIteratorY yBeg
					)
				{
					typedef typename std::iterator_traits<TInputIterator>::value_type value_type;
					value_type   m1, m2;
					value_type   s1, s2;
					value_type   corr;
					m1 = m2 = s1 = s2 = 0.0;
					corr = 0.0;
					ptrdiff_t n = std::distance(xBeg , xEnd );
					value_type nd= static_cast<value_type>(n);
					for(; xBeg != xEnd ;++xBeg,++yBeg)
					{
						corr += *xBeg * *yBeg;
						m1 += *xBeg;
						m2 += *yBeg;
						s1 += *xBeg* *xBeg;
						s2 += *yBeg* *yBeg;
					}
					m1 /= nd;
					m2 /= nd;
					s1 -= m1*m1*nd;
					s2 -= m2*m2*nd;

					if( s1 < 1.0e-12 || s2 < 1.0e-12 )
						return(0.0);
					else
					{
						corr -= m1*m2*(double)n;
						corr /= sqrt(s1*s2);
						return(corr);
					}
				}
}

#endif

////////11////////21////////31////////41////////51////////61////////71////////8
////////11////////21////////31////////41////////51////////61////////71////////8
