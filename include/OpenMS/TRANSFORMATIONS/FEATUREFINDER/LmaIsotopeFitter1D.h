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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>



namespace OpenMS
{
    /**
      @brief Isotope distribution fitter (1-dim.) approximated using Levenberg-Marquardt algorithm (GSL implementation) for parameter optimization.

      @htmlinclude OpenMS_LmaIsotopeFitter1D.parameters
     */
    class OPENMS_DLLAPI LmaIsotopeFitter1D
    : public LevMarqFitter1D
    {
        public:

            enum Averagines{C=0,H,N,O,S,AVERAGINE_NUM};

            /// Default constructor
            LmaIsotopeFitter1D();

            /// copy constructor
            LmaIsotopeFitter1D(const LmaIsotopeFitter1D& source);

            /// destructor
            virtual ~LmaIsotopeFitter1D();

            /// assignment operator
            virtual LmaIsotopeFitter1D& operator = (const LmaIsotopeFitter1D& source);

            /// create new LmaIsotopeFitter1D object (function needed by Factory)
            static Fitter1D* create()
            {
              return new LmaIsotopeFitter1D();
            }

            /// name of the model (needed by Factory)
            static const String getProductName()
            {
              return "LmaIsotopeFitter1D";
            }

            /// return interpolation model
            QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);

	protected:

          /// Helper struct (contains the size of an area, a raw data container, the relative abundance of i-th isotopic peak and the distance between consecutive isotopic peaks)
          struct Data
          {
            typedef Peak1D PeakType;
            typedef std::vector<PeakType > RawDataArrayType;
            typedef std::vector < double > ContainerType;
            typedef Feature::CoordinateType CoordinateType;

            Size n;
            RawDataArrayType set;
            ContainerType isotopes_exact;
            CoordinateType isotope_distance;
            // bool mono_known;
            // CoordinateType monoisotopic_mz;
            CoordinateType isotopes_stdev;
            CoordinateType sigma;
          };

          /// Compute start parameter
          void setInitialParameters_();

          /// Evaluation of the target function for nonlinear optimization
          static Int residual_(const gsl_vector* x, void* params, gsl_vector* f);

          /// Compute the Jacobian matrix, where each row of the matrix corresponds to a point in the data
          static Int jacobian_(const gsl_vector* x, void* params, gsl_matrix* J);

          /// Driver function for the evaluation of function and jacobian
          static Int evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);

          /** Display the intermediate state of the solution. The solver state contains
              the vector s->x which is the current position, and the vector s->f with
              corresponding function values
          */
          void printState_(Int iter, gsl_multifit_fdfsolver * s);

          /// isotope charge
          UInt charge_;
          /// standard derivation in isotope
          CoordinateType isotope_stdev_;
          /// total intensity (area under curve)
          CoordinateType total_intensity_;
          /// monoisotopic mass
          CoordinateType monoisotopic_mz_;
          /// maximum isotopic rank to be considered
          Int max_isotope_;
          /// cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered
          DoubleReal trim_right_cutoff_;
          /// distance between consecutive isotopic peaks
          DoubleReal isotope_distance_;
          /// Centroid m/z (as opposed to monoisotopic m/z)
          CoordinateType mean_;
          /// number of an atom per Dalton of mass
          DoubleReal averagine_[AVERAGINE_NUM];
          /// relative abundance of i-th isotopic peak
          ContainerType isotopes_exact_;
          /// The position of the monoisotopic mass is known(=1) or unknown(=0).
          bool monoisotopic_mass_known_;

	  			void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H
