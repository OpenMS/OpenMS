// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_MATH_STATISTICS_GUMBELDISTRIBUTIONFITTER_H
#define OPENMS_MATH_STATISTICS_GUMBELDISTRIBUTIONFITTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <vector>

// gsl includes
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>


namespace OpenMS
{
	namespace Math
	{
	  /** 
	  	@brief Implements a fitter for the Gumbel distribution.
	
	    This class fits a Gumbel distribution to a number of data points.
	    The results as well as the initial guess are specified using the struct 
			GumbelDistributionFitResult.
	
	    The formula with the fitted parameters can be transformed into a
	    gnuplot formula using getGnuplotFormula() after fitting.
	
			The implementation is done using GSL fitting algorithms.
			
			@ingroup Math
		*/
		class OPENMS_DLLAPI GumbelDistributionFitter
		{
			public:
	
				/// struct to represent the parameters of a gumbel distribution
				struct GumbelDistributionFitResult
				{
					public:
						
						GumbelDistributionFitResult()
							: a(1.0),
								b(2.0)
						{
						}

						GumbelDistributionFitResult(const GumbelDistributionFitResult& rhs)
							: a(rhs.a),
								b(rhs.b)
						{
						}

						GumbelDistributionFitResult& operator = (const GumbelDistributionFitResult& rhs)
						{
							if (this != &rhs)
							{
								a = rhs.a;
								b = rhs.b;
							}
							return *this;
						}
									
						/// location parameter a
						double a;
	
						/// scale parameter b
						double b;
				};
		
				/// Default constructor
				GumbelDistributionFitter();
				/// Destructor
				virtual ~GumbelDistributionFitter();
				
				/// sets the gumbel distribution start parameters a and b for the fitting 
				void setInitialParameters(const GumbelDistributionFitResult& result);
	
				/** 
					@brief Fits a gumbel distribution to the given data points
	
					@param points Input parameter which represents the point used for the fitting
	
					@exception Exception::UnableToFit is thrown if fitting cannot be performed
				*/
				GumbelDistributionFitResult fit(std::vector<DPosition<2> >& points);
	
				/// returns the gnuplot formula of the fitted gumbel distribution
				const String& getGnuplotFormula() const;
				
			protected:
				
				static int gumbelDistributionFitterf_(const gsl_vector* x, void* params, gsl_vector* f);
	
				static int gumbelDistributionFitterdf_(const gsl_vector* x, void* params, gsl_matrix* J);
	
				static int gumbelDistributionFitterfdf_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);
	
				void printState_(size_t iter, gsl_multifit_fdfsolver* s);
				
				GumbelDistributionFitResult init_param_;
				
				String gnuplot_formula_;
	
			private:
				/// Copy constructor (not implemented)
				GumbelDistributionFitter(const GumbelDistributionFitter& rhs);
				/// assignment operator (not implemented)
				GumbelDistributionFitter& operator = (const GumbelDistributionFitter& rhs);
		};
	}
}

#endif

