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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LOGNORMALFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LOGNORMALFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
   
    /** 
        @brief LogNormal distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (GSL implementation) for parameter optimization.
        @ingroup FeatureFinder
    */
    class LogNormalFitter1D 
    : public LevMarqFitter1D
    {
        public:
            
            /// Default constructor
            LogNormalFitter1D();
            
            /// copy constructor
            LogNormalFitter1D(const LogNormalFitter1D& source);
            
            /// destructor
            virtual ~LogNormalFitter1D();
            
            /// assignment operator
            virtual LogNormalFitter1D& operator = (const LogNormalFitter1D& source);
            
            /// create new BiGaussModel object (function needed by Factory)
            static Fitter1D* create()
            {
                return new LogNormalFitter1D();
            }
            
            /// name of the model (needed by Factory)
            static const String getProductName()
            {
                return "LogNormalFitter1D";
            }
              
            /// return interpolation model
            QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);
            
	     protected:
          
          /// Helper struct (contains the size of an area and a raw data container)
          struct Data
          {
            typedef RawDataPoint1D RawDataPointType;
            typedef DPeakArray<RawDataPointType > RawDataArrayType;
                
            UInt n;
            RawDataArrayType set;
          };
         
          /// Compute start parameter
          void setInitialParameters_(const RawDataArrayType& set);

          /// Evaluation of the target function for nonlinear optimization
          static Int residual_(const gsl_vector* x, void* params, gsl_vector* f);
        
          /// Compute the Jacobian matrix, where each row of the matrix corresponds to a point in the data
          static Int jacobian_(const gsl_vector* x, void* params, gsl_matrix* J);
        
          /// Driver function for the evaluation of function and jacobian
          static Int evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);
        
          /** Diplay the intermediate state of the solution. The solver state contains 
              the vector s->x which is the current position, and the vector s->f with 
              corresponding function values */
          void printState_(Int iter, gsl_multifit_fdfsolver * s);  
            
          /// parameter of log normal function - ratio between h and the height at which w and s are computed
          CoordinateType r_;
          /// Parameter of emg - peak height
          CoordinateType height_;
          /// Parameter of emg - peak width
          CoordinateType width_;
          /// Parameter of emg - peak symmetry
          CoordinateType symmetry_;
          /// Parameter of emg - peak retention time
          CoordinateType retention_;
       
          void updateMembers_();
    };

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LOGNORMALFITTER1D_H
