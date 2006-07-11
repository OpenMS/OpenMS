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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H
#define OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/Param.h>

namespace OpenMS
{
  /**
    @brief Linear Resampling of raw data.

	  This class can be used to generate uniform data from non-uniform raw data (e.g. ESI-TOF or MALDI-TOF experiments).
	  Therefore the intensity at every position x in the input raw data is spread to the two
	  adjacent resampling points.
	  This method preserves the area of the input signal and also the centroid position of a peak.
	  Therefore it is recommended for quantitation as well as for identification experiments.
	
	  @note Use this method only for high resoluted data (< 0.1 Th between two adjacent raw data points).
       		The resampling rate should be >= the accuracy.

    @todo write tests (Eva)

   */


  template <typename PeakType>
  class LinearResampler
  {

  public:

    ///\name Typedefs
    //@{
    ///
    typedef DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
    ///
    typedef typename MSExperiment< PeakType >::iterator SpectrumIterator;
    ///
    typedef typename MSExperiment< PeakType >::const_iterator SpectrumConstIterator;
    ///
    typedef MSSpectrum< PeakType > Spectrum;
    ///
    typedef typename Spectrum::const_iterator ConstPeakIterator;
    ///
    typedef typename Spectrum::iterator PeakIterator;

    //@}

    /**@brief
    */
    LinearResampler()
        : spacing_(0.05)
    {
      ms_exp_resampled = 0;
    }

    LinearResampler(const Param& parameters)
    {
      param_ = parameters;

      // if a parameter is missed in the param object the value should be substituted by a dv value
      DataValue dv = param_.getValue("ResamplingWidth");
      if (dv.isEmpty() || dv.toString() == "") spacing_ = 0.05;
      else spacing_ = (double)dv;

      ms_exp_resampled = 0;
    }

    /// Copy constructor.
    LinearResampler( LinearResampler const & lr )
        : param_(lr.param_),
        spacing_(lr.spacing_)
    {
      ms_exp_resampled = 0;
    }

    /// Destructor.
    ~LinearResampler()
    {}

    /// Assignment operator
    LinearResampler& operator= (const LinearResampler& source)
    {
      if (&source == this) return *this;

      // Note: you have to set the ms_exp_resampled by your own!
      ms_exp_resampled = 0;

      param_ = source.param_;
      spacing_ = source.spacing_;
      return *this;
    }


    /// Take the reference to the output data
    LinearResampler& operator()(MSExperiment< PeakType >& ms_exp)
    {
      ms_exp_resampled = &ms_exp;

      return *this;
    }


    // assumes that the resampled_first iterator points on a spectrum with proper size (containing peaks with zero intensity)
    void start(ConstPeakIterator first, ConstPeakIterator last, PeakIterator resampled_first)
    {
      double end_pos = (last-1)->getPos();
      double start_pos = first->getPos();
      int number_raw_points = distance(first,last);
      int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing_ + 1));

      // generate the resampled peaks at positions origin+i*spacing_
      for (int i=0; i < number_resampled_points; ++i)
      {
        (resampled_first+i)->getPos() = start_pos+ i*spacing_;
      }

      // spread the intensity h of the data point at position x to the left and right
      // adjacent resampled peaks
      double distance_left = 0.;
      double distance_right = 0.;
      int left_index = 0;
      int right_index = 0;

      for (int i=0; i < number_raw_points ; ++i)
      {
        left_index = (int)floor(((first+i)->getPos() - start_pos) / spacing_);
        right_index = left_index + 1;

        // compute the distance between x and the left adjacent resampled peak
        distance_left = fabs((first+i)->getPos() - (resampled_first + left_index)->getPos()) / spacing_;
        // compute the distance between x and the right adjacent resampled peak
        distance_right = fabs((first+i)->getPos() - (resampled_first + right_index)->getPos());

        // add the distance_right*h to the left resampled peak and distance_left*h to the right resampled peak
        (resampled_first + left_index)->getIntensity() += (first+i)->getIntensity()*distance_right / spacing_;
        (resampled_first + right_index)->getIntensity() += (first+i)->getIntensity()*distance_left;
      }
    }



    ///\name Accessors
    //@{
    /// Non-mutable access to the resampled data
    inline const MSExperiment<PeakType>* getResampledDataPointer() const { return ms_exp_resampled; }
    /// Mutable access to the resampled data
    inline MSExperiment<PeakType>* getResampledDataPointer() { return ms_exp_resampled; }
    /// Mutable access to the resampled data
    inline void setResampledDataPointer(MSExperiment<PeakType>& ms_exp) { ms_exp_resampled = &ms_exp; }

    /// Non-mutable access to spacing
    inline const double getSpacing() const { return spacing_; }
    /// Mutable access to the spacing
    inline double& getSpacing() { return spacing_; }
    /// Mutable access to the spacing
    inline void setSpacing(const double spacing) { spacing_ = spacing; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param)
    {
      param_ = param;

      // set the new values
      DataValue dv = param_.getValue("ResamplingWidth");
      if (!(dv.isEmpty() || dv.toString() == "")) spacing_ = (double)dv;
    }
    //@}


  protected:
    // ms_exp_resampled points to the resampled MSExperiment
    MSExperiment<PeakType>* ms_exp_resampled;

    // Parameter object
    Param param_;

    // spacing of the resampled data
    double spacing_;
  };

  template < typename PeakType >
  const MSExperiment< PeakType >&
  operator>>(const MSExperiment< PeakType >& ms_exp, LinearResampler< PeakType >& lr)
  {
    *static_cast<ExperimentalSettings*>(lr.getResampledDataPointer()) = ms_exp;

    typename LinearResampler< PeakType >::SpectrumConstIterator first_scan = ms_exp.begin();
    typename LinearResampler< PeakType >::SpectrumConstIterator last_scan = ms_exp.end();

    while (first_scan != last_scan)
    {
      typename MSSpectrum< PeakType >::const_iterator first_data_point = first_scan->begin();
      typename MSSpectrum< PeakType >::const_iterator last_data_point = first_scan->end();

      double end_pos = (last_data_point-1)->getPos();
      double start_pos = first_data_point->getPos();
      int number_resampled_points = (int)ceil((end_pos -start_pos) / lr.getSpacing() + 1);

      DPeakArrayNonPolymorphic<1,PeakType> data(number_resampled_points);

      MSSpectrum< PeakType > resampled_spec;
      resampled_spec.setContainer(data);

	    resampled_spec.setRetentionTime(first_scan->getRetentionTime(), first_scan->getRetentionTimeStart(), first_scan->getRetentionTimeStop());
	    resampled_spec.setMSLevel(first_scan->getMSLevel());
	    resampled_spec.setName(first_scan->getName());

      lr.start(first_data_point,last_data_point,resampled_spec.begin());

      lr.getResampledDataPointer()->std::vector< MSSpectrum< PeakType > >::push_back(resampled_spec);
      ++first_scan;
    }
    return *lr.getResampledDataPointer();
  }

} // namespace OpenMS

#endif // OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLER_H
