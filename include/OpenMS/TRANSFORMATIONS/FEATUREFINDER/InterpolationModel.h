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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

namespace OpenMS
{
  /** 
  	@brief Abstract class for 1D-models that are approximated using linear interpolation
		
		Model wrapping LinearInterpolation for speed-up in calculation of predicted intensities
		Derived classes have to implement setSamples()
		
    @htmlinclude OpenMS_InterpolationModel.parameters

		@ingroup FeatureFinder

	*/
    class OPENMS_DLLAPI InterpolationModel
    	: public BaseModel<1>
    {

      public:
			typedef DoubleReal IntensityType;
      typedef DPosition<1> PositionType;
			typedef DoubleReal CoordinateType;
			typedef Math::LinearInterpolation<DoubleReal> LinearInterpolation;

      /// Default constructor
      InterpolationModel()
				: BaseModel<1>(),
					interpolation_()
			{
				this->defaults_.setValue("interpolation_step",0.1,"Sampling rate for the interpolation of the model function ");
				this->defaults_.setValue("intensity_scaling",1.0,"Scaling factor used to adjust the model distribution to the intensities of the data");
			}

      /// copy constructor
      InterpolationModel(const InterpolationModel& source)
				: BaseModel<1>(source),
					interpolation_(source.interpolation_),
					interpolation_step_(source.interpolation_step_),
					scaling_(source.scaling_)
			{
			}

      /// destructor
      virtual ~InterpolationModel()
      {
      }

      /// assignment operator
      virtual InterpolationModel& operator = (const InterpolationModel& source)
			{
				if (&source ==this) return *this;

				BaseModel<1>::operator = (source);
				interpolation_step_ = source.interpolation_step_;
				interpolation_ = source.interpolation_;
				scaling_ = source.scaling_;

				return *this;
			}

      /// access model predicted intensity at position @p pos
      IntensityType getIntensity(const PositionType& pos) const
			{
				return interpolation_.value(pos[0]);
			}

			/// access model predicted intensity at position @p pos
      IntensityType getIntensity(CoordinateType coord) const
			{
				return interpolation_.value(coord);
			}

			/// Returns the interpolation class
			const LinearInterpolation& getInterpolation() const
			{
				return interpolation_;
			}

			/** @brief get the scaling for the model

				A scaling factor of @p scaling means that the area under the model equals
				@p scaling. Default is 1.
			*/
			CoordinateType getScalingFactor() const
			{
				return scaling_;
			}

			/** @brief set the offset of the model

				The whole model will be shifted to the new offset without being recomputed all over.
				Setting takes affect immediately.
			*/
			virtual void setOffset(CoordinateType offset)
			{
				interpolation_.setOffset(offset);
			}

			/// get reasonable set of samples from the model (i.e. for printing)
			void getSamples(SamplesType& cont) const
			{
				cont = SamplesType();
				BaseModel<1>::PeakType peak;
				for (Size i=0; i<interpolation_.getData().size(); ++i)
				{
					peak.setIntensity( interpolation_.getData()[i] );
					peak.getPosition()[0] = interpolation_.index2key(i);
					cont.push_back(peak);
				}
			}

			/// "center" of the model, particular definition (depends on the derived model)
			virtual CoordinateType getCenter() const
			{
				throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				return CoordinateType(); // we will never get here, but this avoids a warning
			}

			/// set sample/supporting points of interpolation wrt params.
			virtual void setSamples()
			{
				throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			};

			/**
				@brief Set the interpolation step for the linear interpolation of the model

				For setting to take affect, call setSamples().
			*/
			void setInterpolationStep(CoordinateType interpolation_step)
			{
				interpolation_step_ = interpolation_step;
				this->param_.setValue("interpolation_step",interpolation_step_);
			}

			void setScalingFactor(CoordinateType scaling)
			{
			  scaling_ = scaling;
				this->param_.setValue("intensity_scaling",scaling_);
			}

		protected:
			LinearInterpolation interpolation_;
			CoordinateType interpolation_step_;
			CoordinateType scaling_;

			void updateMembers_()
			{
				BaseModel<1>::updateMembers_();
				interpolation_step_ = this->param_.getValue("interpolation_step");
				scaling_ = this->param_.getValue("intensity_scaling");
			}
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
