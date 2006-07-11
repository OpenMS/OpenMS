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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_SMOOTHING_DGAUSSFILTER_H
#define OPENMS_FILTERING_SMOOTHING_DGAUSSFILTER_H

#include <OpenMS/FILTERING/SMOOTHING/DSmoothFilter.h>

#include <OpenMS/FORMAT/Param.h>

#include <math.h>

//#define DEBUG_FILTERING
//#undef DEBUG_FILTERING

namespace OpenMS
{
  using namespace Math;

  /**
    @brief This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform raw data.

    Gaussian filters are important in many signal processing,
    image processing, and communication applications. These filters are characterized by narrow bandwidths,
    sharp cutoffs, and low passband ripple. A key feature of Gaussian filters is that the Fourier transform of a
    Gaussian is also a Gaussian, so the filter has the same response shape in both the time and frequency domains.
    The coefficients \f$ \emph{coeffs} \f$ of the Gaussian-window with length \f$ \emph{frameSize} \f$ are calculated
    from the gaussian distribution
    \f[ \emph{coeff}(x) = \frac{1}{\sigma \sqrt{2\pi}} e^{\frac{-x^2}{2\sigma^2}} \f]
    where \f$ x=[-\frac{frameSize}{2},...,\frac{frameSize}{2}] \f$ represents the window area and \f$ \sigma \f$
    is the standard derivation.

    @note The wider the kernel width the smoother the signal (the more detail information get lost!).
          Use a gaussian filter kernel which has approximately the same width as your mass peaks,
          whereas the gaussian peak width corresponds approximately to 8*sigma.
  */

  template <Size D, typename MapType = MSExperiment< DRawDataPoint<1> > >
  class DGaussFilter : public DSmoothFilter<D, MapType>
  {
  public:
    /** @name Type definitions
     */
    //@{
    ///
    typedef DSmoothFilter<D, MapType> BaseClass;
    typedef typename BaseClass::DimensionDescription DimensionDescription;
    typedef typename BaseClass::RawDataConstIterator RawDataConstIterator;
    typedef typename BaseClass::RawDataIterator RawDataIterator;
    typedef typename BaseClass::RawData RawData;
    ///
    using BaseClass::raw_filtered_;
    using BaseClass::mz_dim_;
    using BaseClass::rt_dim_;
    using BaseClass::coeffs_;
    ///
    //@}

    /** @name Constructors and Destructor
     */

    //@{
    inline DGaussFilter()
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }

      sigma_ = .1;
      spacing_ = 0.01;

      //compute the filter kernel coefficients
      init(sigma_,spacing_);
      raw_filtered_=0;
    }

    /// Class constructor setting the coefficients in the filter. Please note that the frame_size must be odd.
    inline DGaussFilter(const Param& parameters) throw (Exception::InvalidValue)
    {
      param_ = parameters;

      if (D == 1) // one dimensional picking
      {
        mz_dim_ = 0;
        rt_dim_ = -1;
      }
      else // due to our precondition, we know that D == 2
      {
        rt_dim_=0;
        mz_dim_=1;
      }

      // if a smoothing parameter is missed in the param object the value should be substituted by a dv value
      DataValue dv;
      float kernel_width = 0.;

      dv = param_.getValue("GaussianWidth");
      if (dv.isEmpty() || dv.toString() == "") kernel_width = .8;
      else kernel_width = (float)dv;

      std::cout << "KERNEL " << kernel_width << std::endl;

      if (kernel_width <= 0)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The kernel width should be greater than zero!",String(kernel_width));
      }

      // The standard deviation corresponds approximately kernel_width / 8
      sigma_ = kernel_width / 8.;

      //compute the filter kernel coefficients with at least 50 data points
			spacing_= 4*sigma_ / 50;
      init(sigma_,spacing_);

      raw_filtered_=0;
    }
    ///
    inline DGaussFilter(const DGaussFilter& g)
        : DSmoothFilter<D>(g),
        sigma_(g.sigma_),
        spacing_(g.spacing_),
        param_(g.param_)
    { }
    ///
    virtual ~DGaussFilter()
    {   }

    //@}

    /** @name Assignment
     */
    //@{
    ///
    inline DGaussFilter& operator=(const DGaussFilter& s)
    {
      param_ = s.param_;

      spacing_=s.spacing_;
      raw_filtered_=s.raw_filtered_;
      mz_dim_=s.mz_dim_;
      rt_dim_=s.rt_dim_;
      coeffs_=s.coeffs_;
      sigma_=s.sigma_;

      return *this;
    }
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to the sigma
    inline const double& getSigma() const { return sigma_; }
    /// Mutable access to the sigma
    inline void setSigma(const double sigma)
    {
      sigma_=sigma;
	  spacing_ = 4*sigma_ / 50;
      init(sigma_,spacing_);
    }

    /// Non-mutable access to the kernel width
    inline double getKernelWidth() const
    {
      return (sigma_ * 8.);
    }
    /// Mutable access to the kernel width
    inline void setKernelWidth(const double kernel_width)
    {
      sigma_ = kernel_width / 8.;
      init(sigma_,spacing_);
    }

    /// Non-mutable access to the spacing
    inline const double& getSpacing() const { return spacing_; }
    /// Mutable access to the spacing
    inline void setSpacing(const double spacing)
    {
      spacing_=spacing;
	  OPENMS_PRECONDITION((4*sigma_ > spacing), "You have to choose a smaller spacing for the kernel coefficients!" );
      init(sigma_,spacing_);
    }


    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) throw (Exception::InvalidValue)
    {
      param_ = param;
      DataValue dv = param_.getValue("GaussianWidth");

      if (!(dv.isEmpty() || dv.toString() == ""))
      {
        double kernel_width = 0.;
        kernel_width = (float)dv;

        if (kernel_width <= 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The kernel width should be greater than zero!",String(kernel_width));
        }

        sigma_ = kernel_width / 8.;
        init(sigma_,spacing_);
      }
    }
    //@}


    /** Build a gaussian distribution for the current spacing and standard deviation.
    We store the coefficiens of gaussian in the vector<double> coeffs_;

    We only need a finite amount of points since the gaussian distribution
    decays fast. We take 4*sigma (99.993666% of the area is within four standard deviations), since at that point the function
    has dropped to ~ -10^-4
    */
    void init(float sigma, float spacing)
    {
      sigma_= sigma;
      spacing_=spacing;

      int number_of_points_right = (int)(ceil(4*sigma_ / spacing_))+1;
      coeffs_.resize(number_of_points_right);
      coeffs_[0] = 1.0/(sigma_ * sqrt(2.0 * M_PI));

      for (int i=1; i < number_of_points_right; i++)
      {
        coeffs_[i] = gauss_(i*spacing_);
      }
#ifdef DEBUG_FILTERING
      std::cout << "Coeffs: " << std::endl;
      for (int i=0; i < number_of_points_right; i++)
      {
        std::cout << i*spacing_ << ' ' << coeffs_[i] << std::endl;
      }
#endif

    }

  protected:
    /// The standard derivation  \f$ \sigma \f$.
    double sigma_;
    /// The spacing of the pre-tabulated kernel coefficients
    double spacing_;
    /// Parameter object
    Param param_;
		

    // computes the value of the gaussian distribution (mean=0 and standard deviation=sigma) at position x
    double gauss_(double x)
    {
      return (1.0/(sigma_ * sqrt(2.0 * M_PI)) * exp(-(x*x) / (2 * sigma_ * sigma_)));
    }

    ///
    virtual void convolute_(RawDataConstIterator it_begin, RawDataConstIterator it_end, RawDataIterator new_raw_first)
    {
#ifdef DEBUG_FILTERING
      std::cout << "KernelWidth: " << 8*sigma_ << std::endl;
#endif
      RawDataConstIterator help = it_begin;
      while (help != it_end)
      {
        new_raw_first->setPosition(help->getPosition());
        new_raw_first->setIntensity(integrate_(help,it_begin,it_end));
        ++new_raw_first;
        ++help;
      }
    }

    double integrate_(RawDataConstIterator x, RawDataConstIterator first, RawDataConstIterator last)
    {
      double v = 0.;
      // norm the gaussian kernel area to one
      double norm = 0.;
      int middle = coeffs_.size();

      double start_pos = ((x->getPosition()[mz_dim_]-(middle*spacing_)) > first->getPosition()[mz_dim_]) ? (x->getPosition()[mz_dim_]-(middle*spacing_))
                         : first->getPosition()[mz_dim_];
      double end_pos = ((x->getPosition()[mz_dim_]+(middle*spacing_)) < (last-1)->getPosition()[mz_dim_]) ? (x->getPosition()[mz_dim_]+(middle*spacing_))
                       : (last-1)->getPosition()[mz_dim_];


      RawDataConstIterator help = x;
#ifdef DEBUG_FILTERING
      std::cout << "integrate from middle to start_pos "<< help->getPosition()[mz_dim_] << " until " << start_pos << std::endl;
#endif

      //integrate from middle to start_pos
      while ((help != first) && ((help-1)->getPosition()[mz_dim_] > start_pos))
      {
        // search for the corresponding datapoint of help in the gaussian (take the left most adjacent point)
        double distance_in_gaussian = fabs(x->getPosition()[mz_dim_] - help->getPosition()[mz_dim_]);
        unsigned int left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j=0; ((j<3) &&  (distance(first,help-j) >= 0)); ++j)
        {
          if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
          {
            left_position +=j;
            break;
          }
        }

        // interpolate between the left and right data points in the gaussian to get the true value at position distance_in_gaussian
        int right_position = left_position+1;
        double d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
        // check if the right data point in the gaussian exists
        double coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];
#ifdef DEBUG_FILTERING
        std::cout << "distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " right_position " << right_position << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at left_position "  <<  coeffs_[left_position] << std::endl;
        std::cout << "coeffs_ at right_position "  <<  coeffs_[right_position] << std::endl;
        std::cout << "interpolated value left " << coeffs_right << std::endl;
#endif


        // search for the corresponding datapoint for (help-1) in the gaussian (take the left most adjacent point)
        distance_in_gaussian = fabs(x->getPosition()[mz_dim_] - (help-1)->getPosition()[mz_dim_]);
        left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j=0; ((j<3) && (distance(first,help-j) >= 0)); ++j)
        {
          if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
          {
            left_position +=j;
            break;
          }
        }

        // start the interpolation for the true value in the gaussian
        right_position = left_position+1;
        d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
        double coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                            : coeffs_[left_position];
#ifdef DEBUG_FILTERING
        std::cout << " help-1 " << (help-1)->getPosition()[mz_dim_] << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " right_position " << right_position << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[left_position]<<std::endl;
        std::cout << "coeffs_ at right_position " <<   coeffs_[right_position]<<std::endl;
        std::cout << "interpolated value right " << coeffs_left << std::endl;

        std::cout << " intensity " << fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. << " * " << (help-1)->getIntensity() << " * " << coeffs_left <<" + " << (help)->getIntensity()<< "* " << coeffs_right
        << std::endl;
#endif


        norm += fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. * (coeffs_left + coeffs_right);

        v+= fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. * ((help-1)->getIntensity()*coeffs_left + help->getIntensity()*coeffs_right);
        --help;
      }


      //integrate from middle to end_pos
      help = x;
#ifdef DEBUG_FILTERING
      std::cout << "integrate from middle to endpos "<< (help)->getPosition()[mz_dim_] << " until " << end_pos << std::endl;
#endif
      while ((help != (last-1)) && ((help+1)->getPosition()[mz_dim_] < end_pos))
      {
        // search for the corresponding datapoint for help in the gaussian (take the left most adjacent point)
        double distance_in_gaussian = fabs(x->getPosition()[mz_dim_] - help->getPosition()[mz_dim_]);
        int left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j=0; ((j<3) && (distance(help+j,last-1) >= 0)); ++j)
        {
          if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
          {
            left_position +=j;
            break;
          }
        }
        // start the interpolation for the true value in the gaussian
        int right_position = left_position+1;
        double d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
        double coeffs_left= (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                            : coeffs_[left_position];

#ifdef DEBUG_FILTERING
        std::cout << " help " << (help)->getPosition()[mz_dim_] << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at right_position " <<  coeffs_[left_position]<<std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
        std::cout << "interpolated value left " << coeffs_left << std::endl;
#endif

        // search for the corresponding datapoint for (help+1) in the gaussian (take the left most adjacent point)
        distance_in_gaussian = fabs(x->getPosition()[mz_dim_] - (help+1)->getPosition()[mz_dim_]);
        left_position = (unsigned int)floor(distance_in_gaussian / spacing_);

        // search for the true left adjacent data point (because of rounding errors)
        for (int j=0; ((j<3) && (distance(help+j,last-1) >= 0)); ++j)
        {
          if (((left_position-j)*spacing_ <= distance_in_gaussian) && ((left_position-j+1)*spacing_ >= distance_in_gaussian))
          {
            left_position -= j;
            break;
          }

          if (((left_position+j)*spacing_ < distance_in_gaussian) && ((left_position+j+1)*spacing_ < distance_in_gaussian))
          {
            left_position +=j;
            break;
          }
        }

        // start the interpolation for the true value in the gaussian
        right_position = left_position+1;
        d = fabs((left_position*spacing_)-distance_in_gaussian) / spacing_;
        double coeffs_right = (right_position < middle) ? (1-d)*coeffs_[left_position]+d*coeffs_[right_position]
                              : coeffs_[left_position];
#ifdef DEBUG_FILTERING
        std::cout << " (help + 1) " << (help+1)->getPosition()[mz_dim_] << " distance_in_gaussian " << distance_in_gaussian << std::endl;
        std::cout << " left_position " << left_position << std::endl;
        std::cout << "coeffs_ at right_position " <<   coeffs_[left_position]<<std::endl;
        std::cout << "coeffs_ at left_position " <<  coeffs_[right_position]<<std::endl;
        std::cout << "interpolated value right " << coeffs_right << std::endl;

        std::cout << " intensity " <<  fabs(help->getPosition()[mz_dim_] - (help+1)->getPosition()[mz_dim_]) / 2.
        << " * " << help->getIntensity() << " * " << coeffs_left <<" + " << (help+1)->getIntensity()
        << "* " << coeffs_right
        << std::endl;
#endif

        norm += fabs((help-1)->getPosition()[mz_dim_]-help->getPosition()[mz_dim_]) / 2. * (coeffs_left + coeffs_right);

        v+= fabs(help->getPosition()[mz_dim_] - (help+1)->getPosition()[mz_dim_]) / 2. * (help->getIntensity()*coeffs_left + (help+1)->getIntensity()*coeffs_right);
        ++help;
      }

			if (v > 0)
				{
					return v / norm;
				}
			else
				{
					return 0;
				}
    }

  };

} // namespace OpenMS
#endif
