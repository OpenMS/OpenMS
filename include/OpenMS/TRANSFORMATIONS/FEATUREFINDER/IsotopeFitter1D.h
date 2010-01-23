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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
    /** 
      @brief Isotope distribution fitter (1-dim.) approximated using linear interpolation.
                 
      @htmlinclude OpenMS_IsotopeFitter1D.parameters                 
     */
    class OPENMS_DLLAPI IsotopeFitter1D
    : public MaxLikeliFitter1D
    {
        public:
	
            /// Default constructor
            IsotopeFitter1D();
        
            /// copy constructor
            IsotopeFitter1D(const IsotopeFitter1D& source);
        
            /// destructor
            virtual ~IsotopeFitter1D();
        
            /// assignment operator
            virtual IsotopeFitter1D& operator = (const IsotopeFitter1D& source);
      
            /// create new IsotopeFitter1D object (function needed by Factory)
            static Fitter1D* create()
            {
              return new IsotopeFitter1D();
            }
      
            /// name of the model (needed by Factory)
            static const String getProductName()
            {
              return "IsotopeFitter1D";
            }
            
            /// return interpolation model
            QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);
		
        protected:
	
            /// isotope charge
            CoordinateType charge_;
            /// standard derivation in isotope
            CoordinateType isotope_stdev_;
            /// maximum isotopic rank to be considered
          	Int max_isotope_;
      
            void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H
