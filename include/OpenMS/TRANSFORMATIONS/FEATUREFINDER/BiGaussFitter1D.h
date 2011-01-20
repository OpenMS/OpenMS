// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BIGAUSSFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BIGAUSSFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


namespace OpenMS
{
    /** 
      @brief Bigaussian distribution fitter (1-dim.) approximated using linear interpolation.
                 
      @htmlinclude OpenMS_BiGaussFitter1D.parameters
    */
    class OPENMS_DLLAPI BiGaussFitter1D
    : public MaxLikeliFitter1D
    {
        public:
	
            /// Default constructor
            BiGaussFitter1D();
        
            /// copy constructor
            BiGaussFitter1D(const BiGaussFitter1D& source);
        
            /// destructor
            virtual ~BiGaussFitter1D();
        
            /// assignment operator
            virtual BiGaussFitter1D& operator = (const BiGaussFitter1D& source);
    
            /// create new BiGaussModel object (function needed by Factory)
            static Fitter1D* create()
            {
              return new BiGaussFitter1D();
            }
            
            /// return interpolation model
            QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);
    
            /// name of the model (needed by Factory)
            static const String getProductName()
            {
              return "BiGaussFitter1D";
            }
		
        protected:
			
            /// statistics for first peak site
            Math::BasicStatistics<> statistics1_;
            /// statistics for second peak site
            Math::BasicStatistics<> statistics2_;
			
	  void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BIGAUSSFITTER1D_H
