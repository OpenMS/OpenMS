// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_MATH_STATISTICS_GAUSSFITTER_H
#define OPENMS_MATH_STATISTICS_GAUSSFITTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <vector>

// gsl includes
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

namespace OpenMS
{
	class GaussFitter
	{
		public:

			struct GaussFitResult
			{
				public:
					double A;
					double x0; 
					double sigma;
			};
		
			GaussFitter();

			GaussFitter(const GaussFitter&);

			virtual ~GaussFitter();

			GaussFitter& operator = (const GaussFitter&);

			GaussFitResult fit(std::vector<DPosition<2> >&);

			const GaussFitResult& getInitialParameters() const;

			void setInitialParameters(const GaussFitResult&);

			const String& getGnuplotFormula() const;
			
		protected:
			
			static int gauss_fitter_f_(const gsl_vector* x, void* params, gsl_vector* f);

			static int gauss_fitter_df_(const gsl_vector* x, void* params, gsl_matrix* J);

			static int gauss_fitter_fdf_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);
			
			void print_state_(size_t iter, gsl_multifit_fdfsolver * s);
			
			GaussFitResult init_param_;
			
			String gnuplot_formula_;
	};
}

#endif

