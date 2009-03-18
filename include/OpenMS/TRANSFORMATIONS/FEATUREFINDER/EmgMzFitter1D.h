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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGMZFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGMZFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
    /** 
        @brief Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (GSL implementation) for parameter optimization.
                   
        @ref EmgFitter1D_Parameters are explained on a separate page.   
    
        @ingroup FeatureFinder
    */
    class EmgMzFitter1D : EmgFitter1D
    {
        public:
            
            /// Default constructor
            EmgMzFitter1D();
            
            /// copy constructor
            EmgMzFitter1D(const EmgMzFitter1D& source);
            
            /// destructor
            virtual ~EmgMzFitter1D();
            
            /// assignment operator
            virtual EmgMzFitter1D& operator = (const EmgMzFitter1D& source);
                          

						/// create new BiGaussModel object (function needed by Factory)
            static Fitter1D* create()
            {
                return new EmgMzFitter1D();
            }       
            
        protected:
          
          /// Compute start parameter
          virtual void setInitialParameters_(const RawDataArrayType& set);
    };

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGFITTER1D_H
