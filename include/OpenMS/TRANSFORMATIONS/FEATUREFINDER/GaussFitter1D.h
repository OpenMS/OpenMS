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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
    /** 
      @brief Gaussian distribution fitter (1-dim.) approximated using linear interpolation.
                 
      @htmlinclude OpenMS_GaussFitter1D.parameters               
    */
    class OPENMS_DLLAPI GaussFitter1D
    : public MaxLikeliFitter1D
    {
        public:
	
          /// Default constructor
          GaussFitter1D();
      
          /// copy constructor
          GaussFitter1D(const GaussFitter1D& source);
      
          /// destructor
          virtual ~GaussFitter1D();
      
          /// assignment operator
          virtual GaussFitter1D& operator = (const GaussFitter1D& source);

          /// create new GaussFitter1D object (function needed by Factory)
          static Fitter1D* create()
          {
            return new GaussFitter1D();
          }
          
          /// return interpolation model
          QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);
  
          /// name of the model (needed by Factory)
          static const String getProductName()
          {
            return "GaussFitter1D";
          }
		
        protected:
            
          void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSFITTER1D_H
