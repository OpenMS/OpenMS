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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{

  IsotopeModelGeneral::IsotopeModelGeneral()
  : IsotopeModel()
  {
    setName(getProductName());

    defaultsToParam_();
  }

  IsotopeModelGeneral::IsotopeModelGeneral(const IsotopeModelGeneral& source)
  : IsotopeModel(source)
  {
    setParameters( source.getParameters() );
    updateMembers_();
  }

  IsotopeModelGeneral::~IsotopeModelGeneral()
  {
  }

  IsotopeModelGeneral& IsotopeModelGeneral::operator = (const IsotopeModelGeneral& source)
  {
    if (&source == this) return *this;

    IsotopeModel::operator = (source);
    setParameters( source.getParameters() );
    updateMembers_();

    return *this;
  }

  void IsotopeModelGeneral::setSamples(EmpiricalFormula formula)
  {
    // MAGIC alert, num stdev for smooth table for normal distribution
    CoordinateType normal_widening_num_stdev = 4.;
    // Actual width for values in the smooth table for normal distribution
    CoordinateType normal_widening_width = isotope_stdev_ * normal_widening_num_stdev;

    //cout << "isotope_distance_ " << isotope_distance_ << endl;
    //isotope_distance_ = 0.96;
    //cout << "isotope_distance_ " << isotope_distance_ << endl;

    typedef std::vector < double > ContainerType;
    ContainerType isotopes_exact;

    typedef IsotopeDistribution::iterator IsoIter;
    IsotopeDistribution isotope_distribution = formula.getIsotopeDistribution(max_isotope_);
    //isotope_distribution.setTrimRightCutoff(trim_right_cutoff_);
    isotope_distribution.trimRight(trim_right_cutoff_);
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

    Int rMax = std::min( int( left.size() + right.size() - 1 ), int(2*normal_widening_width/interpolation_step_*max_isotope_+1) );
    result.resize(rMax,0);

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for ( SignedSize i = left.size() - 1; i >= 0; --i )
    {
      for ( SignedSize j = std::min<SignedSize>( rMax - i, right.size() ) - 1; j >= 0 ; --j )
      {
        result[i+j] += left[i] * right[j];
      }
    }

    interpolation_.setMapping(interpolation_step_, normal_widening_width / interpolation_step_, mean_ - isotopes_mean)	;
    monoisotopic_mz_ = mean_-isotopes_mean;

    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ / accumulate ( result.begin(), result.end(), IntensityType(0) );

    for ( ContainerType::iterator iter = result.begin(); iter != result.end(); ++iter)
    {
      *iter *= factor;
    }
  }

  void IsotopeModelGeneral::updateMembers_()
  {
    IsotopeModel::updateMembers_();
  }

}
