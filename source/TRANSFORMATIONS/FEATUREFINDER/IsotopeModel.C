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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <numeric>

namespace OpenMS
{
    IsotopeModel::IsotopeModel()
		: InterpolationModel<>(),
			monoisotopic_mz_(0.0)
		{
			setName(getProductName());
			
			defaults_.setValue("averagines:C",0.0443);
			defaults_.setValue("averagines:H",0.0);
			defaults_.setValue("averagines:N",0.0037);
			defaults_.setValue("averagines:O",0.022);
			defaults_.setValue("averagines:S",0.0);
			defaults_.setValue("isotope:trim_right_cutoff",0.001);
			defaults_.setValue("isotope:maximum",1000000);
			defaults_.setValue("isotope:distance",1.000495);
			defaults_.setValue("isotope:stdev",0.1);
			defaults_.setValue("charge",1);
			defaults_.setValue("statistics:mean",0.0);
			
			defaultsToParam_();
		}

  	IsotopeModel::IsotopeModel(const IsotopeModel& source)
		: InterpolationModel<>(source)
		{
			setParam(source.mean_, source.charge_, source.isotope_stdev_);
			updateMembers_();
		}

    IsotopeModel::~IsotopeModel()
    {
    }

   	IsotopeModel& IsotopeModel::operator = (const IsotopeModel& source)
		{
			if (&source == this) return *this;
			
			InterpolationModel<>::operator = (source);
			setParam(source.mean_, source.charge_, source.isotope_stdev_);
			updateMembers_();
			
			return *this;
		}

		void IsotopeModel::setSamples()
		{
			/// MAGIC alert, num stdev for smooth table for normal distribution
			CoordinateType normal_widening_num_stdev = 4.;
			/// Actual width for values in the smooth table for normal distribution
			CoordinateType normal_widening_width = isotope_stdev_
																						* normal_widening_num_stdev;

			typedef std::vector < double > ContainerType;
			ContainerType isotopes_exact;
			CoordinateType mass = mean_ * charge_;

			int C_num = int( 0.5 + mass * averagine_[C]);
			int N_num = int( 0.5 + mass * averagine_[N]);
			int O_num = int( 0.5 + mass * averagine_[O]);
			int H_num = int( 0.5 + mass * averagine_[H]);
			int S_num = int( 0.5 + mass * averagine_[S]);

			String form("");
			if (C_num) form.append("C").append(String(C_num));
			if (H_num) form.append("H").append(String(H_num));
			if (N_num) form.append("N").append(String(N_num));
			if (O_num) form.append("O").append(String(O_num));
			if (S_num) form.append("S").append(String(S_num));

			EmpiricalFormula formula(form);
			typedef IsotopeDistribution::iterator IsoIter;
			IsotopeDistribution isotope_distribution = formula.getIsotopeDistribution(max_isotope_);
			isotope_distribution.setTrimRightCutoff(trim_right_cutoff_);
			isotope_distribution.trimRight();
			isotope_distribution.renormalize();

			// compute the average mass (-offset)
			CoordinateType isotopes_mean = 0;
			int i=0;
			for (	IsoIter iter = isotope_distribution.begin();
						iter != isotope_distribution.end(); ++iter,++i)
			{
				isotopes_exact.push_back(iter->second);
				isotopes_mean += iter->second*i;
			}
			isotopes_mean *= isotope_distance_ / charge_;
			// (Need not divide by sum of probabilities, which is 1.)

			// "stretch" the averagine isotope distribution
			size_t isotopes_exact_size = isotopes_exact.size();
			isotopes_exact.resize(size_t( (isotopes_exact_size-1)
							*isotope_distance_/interpolation_step_+1.6)); // round up a bit more

			for ( size_t i = isotopes_exact_size-1; i; --i )
			// we don't need to move the 0-th entry
			{
				isotopes_exact [size_t(CoordinateType(i)*
								isotope_distance_/interpolation_step_/charge_+0.5)]
				=	isotopes_exact [ i ];
				isotopes_exact [ i ] = 0;
			}

			// compute the normal distribution (to be added for widening the averagine isotope distribution)
			Math::BasicStatistics<> normal_widening_model;
			normal_widening_model.setSum  (1);
			normal_widening_model.setMean (0);
			normal_widening_model.setVariance(isotope_stdev_*isotope_stdev_);
			// fill a container with CoordinateType points
			ContainerType normal_widening_coordinate;
			for ( double coord = -normal_widening_width;
						coord <= normal_widening_width;
						coord += interpolation_step_
					)
			{
				normal_widening_coordinate.push_back(coord);
			}
			// compute normal approximation at these CoordinateType points
			ContainerType normal_widening;
			normal_widening_model.normalApproximation(normal_widening,normal_widening_coordinate );

			// fill linear interpolation
			const ContainerType& left = isotopes_exact;
			const ContainerType& right = normal_widening;
			ContainerType& result = interpolation_.getData();
	    result.clear ();

      int rMax = std::min ( int( left.size() + right.size() - 1 ),
														int(2*normal_widening_width/interpolation_step_*max_isotope_+1) );
      result.resize ( rMax, 0 );

      // we loop backwards because then the small products tend to come first
      // (for better numerics)
      for ( int i = left.size() - 1; i >= 0; --i )
			{
				for ( int j = std::min ( rMax - i, int ( right.size() ) ) - 1; j >= 0 ; --j )
				{
					result[i+j] += left[i] * right[j];
				}
      }

			interpolation_.setOffset(mean_-isotopes_mean-normal_widening_width);
			monoisotopic_mz_ = mean_-isotopes_mean;
			interpolation_.setScale ( interpolation_step_ );

			// scale data so that integral over distribution equals one
			// multiply sum by interpolation_step_ -> rectangular approximation of integral
			IntensityType factor = scaling_ / interpolation_step_ /
						std::accumulate ( result.begin(), result.end(), IntensityType(0) );

			for ( ContainerType::iterator iter = result.begin(); iter != result.end(); ++iter)
			{
				*iter *= factor;
			}
		}

		void IsotopeModel::setParam(CoordinateType mean, UnsignedInt charge, CoordinateType isotope_stdev)
		{
			charge_ = charge;
			isotope_stdev_ = isotope_stdev;
			mean_ = mean;

			param_.setValue("charge", static_cast<SignedInt>(charge_));
			param_.setValue("isotope:stdev",isotope_stdev_);
			param_.setValue("statistics:mean", mean_);
			
			setSamples();
		}

		void IsotopeModel::setOffset(double offset)
		{
			double diff = offset - getInterpolation().getOffset();
			mean_ += diff;
			monoisotopic_mz_ += diff;

			InterpolationModel<>::setOffset(offset);

			param_.setValue("statistics:mean", mean_);
		}

		const IsotopeModel::CoordinateType& IsotopeModel::getOffset()
		{
			return getInterpolation().getOffset();
		}

		UnsignedInt IsotopeModel::getCharge()
		{
			return charge_;
		}

		const IsotopeModel::CoordinateType IsotopeModel::getCenter() const
		{
			return monoisotopic_mz_;
		}

		void IsotopeModel::updateMembers_()
		{
			InterpolationModel<>::updateMembers_();

			charge_ = param_.getValue("charge");
			isotope_stdev_ = param_.getValue("isotope:stdev");
			mean_ = param_.getValue("statistics:mean");
			max_isotope_ = param_.getValue("isotope:maximum");
			trim_right_cutoff_ = param_.getValue("isotope:trim_right_cutoff");
			isotope_distance_ = param_.getValue("isotope:distance");

			averagine_[C] = param_.getValue("averagines:C");
			averagine_[H] = param_.getValue("averagines:H");
			averagine_[N] = param_.getValue("averagines:N");
			averagine_[O] = param_.getValue("averagines:O");
			averagine_[S] = param_.getValue("averagines:S");
			
			setSamples();
		}
}
