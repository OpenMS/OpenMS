// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Early Pearly $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/KERNEL/KernelTraits.h>

namespace OpenMS
{
  /** @brief abstract class for 1D-models that are approximated using linear interpolation
	
			Model wrapping LinearInterpolation for speed-up in calculation of predicted intensities
			Derived classes have to implement setSamples()

			Parameters:
			<table>
			<tr><td>interpolation_step</td>
					<td>step size used to interpolate model</td></tr>
			<tr><td>intensity_scaling</td>
					<td>factor used to scale the calculated intensities</td></tr>
			<tr><td>cutoff</td>
					<td>peak with intensity below cutoff is not considered
							 to be part of the model</td></tr>
			</table>
			
			@ingroup FeatureFinder
			
	*/
  template < typename Traits = KernelTraits >
    class InterpolationModel
    : public BaseModel<1,Traits>
    {

      public:
			typedef typename DPeak<1,Traits>::IntensityType IntensityType;
      		typedef DPosition<1,Traits> PositionType;
			typedef typename PositionType::CoordinateType CoordinateType;
			typedef Math::LinearInterpolation<CoordinateType,IntensityType> LinearInterpolation;
			typedef typename LinearInterpolation::container_type ContainerType;
			typedef DPeakArray<1, DPeak<1,Traits> > SamplesType;

      /// standard constructor
      InterpolationModel()
				: BaseModel<1,Traits>(),
					interpolation_()
			{
				setScalingFactor(1.0);
				setInterpolationStep(0.1);

				this->defaults_.setValue("interpolation_step",0.1);
				this->defaults_.setValue("intensity_scaling",1.0);
			}

      /// copy constructor
      InterpolationModel(const InterpolationModel& source)
				: BaseModel<1,Traits>(source),
					interpolation_(source.interpolation_)
			{
				setScalingFactor(source.scaling_);
				setInterpolationStep(source.interpolation_step_);
			}

      /// destructor
      virtual ~InterpolationModel(){}

      /// assignment operator
      virtual InterpolationModel& operator = (const InterpolationModel& source)
			{
				BaseModel<1,Traits>::operator = (source);
				setInterpolationStep(source.interpolation_step_);
				interpolation_ = source.interpolation_;
				setScalingFactor(source.scaling_);
				return *this;
			}

      /// acess model predicted intensity at position @p pos
      IntensityType getIntensity(const PositionType& pos) const
			{
				return interpolation_.value(pos[0]);
			}
			
			/// acess model predicted intensity at position @p pos
      IntensityType getIntensity(const CoordinateType& coord) const
			{
				return interpolation_.value(coord);
			}

			/** @brief set the interpolation step for the linear interpolation of the model

				For setting to take affect, call setSamples() or setParm.
			*/
			void setInterpolationStep(CoordinateType interpolation_step)
			{
				interpolation_step_ = interpolation_step;
			}

			const LinearInterpolation& getInterpolation() const
			{
				return interpolation_;
			}

			LinearInterpolation& getInterpolation()
			{
				return interpolation_;
			}

			/** @brief set the scaling for the model

				A scaling factor of @p scaling means that the area under the model equals
				@p scaling. Default is 1.
				For setting to take affect, call setSamples() or setParm.
			*/
			void setScalingFactor(CoordinateType scaling)
			{
				scaling_ = scaling;
			}

			/** @brief get the scaling for the model

				A scaling factor of @p scaling means that the area under the model equals
				@p scaling. Default is 1.
			*/
			const CoordinateType& getScalingFactor() const
			{
				return scaling_;
			}

			virtual void setParam(const Param& param)
			{
				BaseModel<1,Traits>::setParam(param);
				interpolation_step_ = this->param_.getValue("interpolation_step");
				scaling_ = this->param_.getValue("intensity_scaling");
			}

    	virtual const Param& getParam() const
			{
				this->param_.setValue("interpolation_step",interpolation_step_);
				this->param_.setValue("intensity_scaling",scaling_);
				return BaseModel<1,Traits>::getParam();
			}

    	virtual Param& getParam()
			{
				this->param_.setValue("interpolation_step",interpolation_step_);
				this->param_.setValue("intensity_scaling",scaling_);
				return BaseModel<1,Traits>::getParam();
			}


			/** @brief set the offset of the model

				The whole model will be shifted to the new offset without being computing all over.
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
				DPeak<1,Traits> peak;
				for (Size i=0; i<interpolation_.getData().size(); ++i)
				{
					peak.setIntensity( interpolation_.getData()[i] );
					peak.getPosition()[0] = interpolation_.index2key(i);
					cont.push_back(peak);
				}
			}

			/// "center" of the model, particular definition (depends on the derived model)
			virtual const CoordinateType getCenter() const=0;

			/// set sample/supporting points of interpolation wrt params.
			virtual void setSamples() =0;

		protected:
			LinearInterpolation interpolation_;
			CoordinateType interpolation_step_;
			CoordinateType scaling_;
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_INTERPOLATIONMODEL_H
