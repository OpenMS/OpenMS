// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch, Sven Nahnsen $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_IDDECOYPROBABILITY_H
#define OPENMS_ANALYSIS_ID_IDDECOYPROBABILITY_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/MATH/STATISTICS/GammaDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Write s.th. useful here

		here also, but with more detail
 		
  */
  class IDDecoyProbability : public DefaultParamHandler
  {
    public:

      /// Default constructor
      IDDecoyProbability();

			/// Copy constructor
			IDDecoyProbability(const IDDecoyProbability& rhs);

			/// Desctructor
			virtual ~IDDecoyProbability();

			/// assignment operator
			IDDecoyProbability& operator = (const IDDecoyProbability& rhs);

		  /// 
			void apply(	std::vector<PeptideIdentification>& prob_ids,
									const std::vector<PeptideIdentification>& fwd_ids, 
									const std::vector<PeptideIdentification>& rev_ids) throw (Exception::MissingInformation);


		protected:

			struct Transformation_
			{
			  double y_factor;
			  double x_factor;
			  double x_shift;
			  double x_max;
			  double y_max_bin;
			};

			void normalizeBins_(const std::vector<double>& scores, HashMap<double, double>& binned, Transformation_& trafo);

			double getProbability_(const GammaDistributionFitter::GammaDistributionFitResult& result_gamma,
														const Transformation_& gamma_trafo,
														const GaussFitter::GaussFitResult& result_gauss,
														const Transformation_& gauss_trafo,
														double score);
			
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDDECOYPROBABILITY_H

