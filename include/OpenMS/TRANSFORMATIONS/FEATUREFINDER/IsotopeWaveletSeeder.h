// -*- C++: make; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff, Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/DATASTRUCTURES/HashMap.h>

#include <gsl/gsl_cdf.h>

#include <limits>

namespace OpenMS
{
	/**
		@brief Seeding module which uses a isotopic wavelet to find seeds.
		
		This seeder select interesting regions in the map using a wavelet funtion
		modelling the distribution of isotopic peak intensities.
		
		Parameters:
		
		<table>
		<tr><td></td><td></td><td>max_charge</td>
		<td>The mother wavelet is precomputed for different charge states. This
		is the maximum charge state considered.</td></tr>
		<tr><td></td><td></td><td>min_charge</td>
		<td>The mother wavelet is precomputed for different charge states. This
		is the minimum charge considered.</td></tr>
		<tr><td></td><td></td><td>intensity_factor</td>
		<td>Scores below the intensity of this point times this parameter are not
		considered for charge estimation.</td></tr>
		<tr><td></td><td></td><td>avg_intensity_factor</td>
		<td>influences the threshold for interesting points in the wavelet transform. </td></tr>
		<tr><td></td><td></td><td>min_samplingrate</td>
		<td>minimum sampling rate (e.g. step size for cwt), usally determined by the average m/z spacing</td></tr>
		<tr><td></td><td></td><td>mass_tolerance_right</td>
		<td>width of seed bounding box to the right. </td></tr>
		<tr><td></td><td></td><td>mass_tolerance_left</td>
		<td>width of seed bounding box to the right. </td></tr>
		<tr><td></td><td></td><td>scans_to_sumup</td>
		<td>number of scans used for alignment</td></tr>
		<tr><td></td><td></td><td>tolerance_scansum</td>
		<td>mass tolerance during point alignment</td></tr>
		</table>
		
		There are additional parameter : @see BaseSweepSeeder
		 
		@ref IsotopeWaveletSeeder_Parameters are explained on a separate page.

		@ingroup FeatureFinder
	*/
  class IsotopeWaveletSeeder
    : public BaseSweepSeeder
  {

  public:

		/// intensity of a peak
		typedef FeaFiTraits::IntensityType IntensityType;
		/// coordinate ( in rt or m/z )
		typedef FeaFiTraits::CoordinateType CoordinateType;
		/// score
		typedef DoubleReal ProbabilityType;

		/// a single MS spectrum
		typedef BaseSweepSeeder::SpectrumType SpectrumType;
		/// a peak
		typedef SpectrumType::PeakType PeakType;
		/// a container of peaks
		typedef SpectrumType::ContainerType ContainerType;
		/// a container of 2D peaks (only temporarily used)
		typedef MSSpectrum< Peak2D >::ContainerType TempContainerType;

		/// charge state estimate with associated score
		typedef BaseSweepSeeder::ScoredChargeType ScoredChargeType;
		/// m/z position in spectrum with charge estimate and score
		typedef BaseSweepSeeder::ScoredMZType ScoredMZType;
		/// container of scored m/z positions
		typedef BaseSweepSeeder::ScoredMZVector ScoredMZVector;

		/// The mother wavelets (one for each charge state)
		typedef std::vector<std::vector<double> > WaveletCollection;
		/// The Hash entry, stores pairs of scan number and list of charge scores
		typedef std::pair<std::list<UInt>, std::list<double> > DoubleList;
		/// Stores the charge states examined
		typedef std::list<UInt> ChargeVector;


    /// Default constructor
    IsotopeWaveletSeeder();

    /// destructor
    virtual ~IsotopeWaveletSeeder();

    /// Copy constructor
    IsotopeWaveletSeeder(const IsotopeWaveletSeeder& rhs);

    /// Assignment operator
    IsotopeWaveletSeeder& operator= (const IsotopeWaveletSeeder& rhs);

    static BaseSeeder* create()
    {
      return new IsotopeWaveletSeeder();
    }

    static const String getProductName()
    {
      return "IsotopeWaveletSeeder";
    }

  protected:

		/** @brief detects isotopic pattern

			@improvement Why is a pointer being used here?  It would be faster to reuse a static member.  Avoid reallocation, pwts can be cleared instead.  (Clemens aking Maintainer)
		*/
		ScoredMZVector detectIsotopicPattern_(SpectrumType& scan);

		/// keeps member and param entries in synchrony
  	virtual void updateMembers_();

		/// Computes m/z spacing of the LC-MS map
		void computeSpacings_();

		/// Precompute and store the gamma function (for the mother wavelet)
		void generateGammaValues_();

		/// Compute null variance
		void computeNullVariance_(const DPeakArray<PeakType >& cwt, const UInt charge );

		/// Compute local variance (in an interval) and test its significance
		ProbabilityType testLocalVariance_(const DPeakArray<PeakType >& cwt, const UInt& start, const UInt charge);

		/**
			@brief Computes the wavelet transform for several charges in nearly the same time.

			The working horse of the discrete-time continuous wavelet transform, but for several charges at the same time.
		 	Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 	this yields the same.
		*/
		void fastMultiCorrelate_(const SpectrumType& signal, std::vector<DPeakArray<PeakType > >* pwts);

		/**
				@brief Returns the lamba parameter of the mother wavelet

				The lambda parameter essentially influences the shape of the wavelet.
		 		Since isotope patterns depend on mass, the wavelet has to adapt its shape.
		 		For more insights look at the formula of the wavelet function.
		*/
		inline CoordinateType getLambda_(CoordinateType real_mass) const
		{
			return (0.035 + 0.000678*real_mass);
		}

		/// The wavelet (mother) function.
 		inline double phiRaw_(double t, double lambda, double a) throw ()
		{
			if (t>2*peak_cut_off_)	return(0);

			Int x0, x1; double f0, f1, fi, res=0;
			x0 = (Int) trunc ((t/a + 1)/min_spacing_);
			x1 = x0+1;
			if ((UInt) x1 < preComputedGamma_.size())
			{
				f0 = preComputedGamma_[x0];
				f1 = preComputedGamma_[x1];
				fi = (f0 + (f1-f0)/((x1-x0)*min_spacing_) * ((t/a+1)-x0*min_spacing_));
				res = (sin(2*M_PI*t/a) * exp(-lambda)) * ((pow(lambda,t/a)) / fi);
				return (res);
			}

			res = (sin(2*M_PI*t/a) * exp(-lambda) * (pow(lambda,t/a)) / tgamma ((t/a)+1));

			return (res);
		}

		UInt findNextMax_(const DPeakArray<PeakType >& cwt, const UInt index);

		/// Assigns scores to each charge state of a isotopic pattern
		ScoredMZVector identifyCharge_(std::vector<DPeakArray<PeakType > >& candidates, SpectrumType& scan);

		/// Interpolates between to data points
		inline double getInterpolatedValue_(double x0, double x, double x1, double f0, double f1) const
		{
			return (f0 + (f1-f0)/(x1-x0) * (x-x0));
		}

		/// Returns an index pair containing the mass/charge @p mz
		inline std::pair<Int, Int> getNearBys_(const SpectrumType& scan, double mz, UInt start=0)
		{
			for (UInt i=start; i<scan.getContainer().size(); ++i)
			{
				if (scan.getContainer()[i].getPos() > mz)
				{
					return std::make_pair(i-1, i);
				}
			}

			//not found
			return std::make_pair(-1,-1);
		}

		/// Returns the absolute mean of the intensities in @p signal
		double getAbsMean_(const DPeakArray<PeakType >& signal,
																		UInt startIndex,
																		UInt endIndex) const;

		/// Does not need much explanation.
		bool wavelet_initialized_;
		/// Number of isotopic peaks a wavelet should contain
		UInt peak_cut_off_;
		/// Length of the mother wavelet
		UInt waveletLength_;
		/// Average spacing in a MS scan
		CoordinateType avMZSpacing_;
		/// Minium spacing
		CoordinateType min_spacing_;
		/// Charge states being tested
		ChargeVector charges_;
		/// Stores the Gamme function
		HashMap<UInt, double> preComputedGamma_;
		/// Determines threshold for the minimum score of a peak
		IntensityType signal_avg_factor_;
		/// Determines threshold for cwt of a peak
		IntensityType cwt_avg_factor_;
		/// Tolerance for scan alignment
		CoordinateType tolerance_scansum_;
		/// variance on empty interval of cwt (used as null hypothesis)
		std::vector<IntensityType> null_var_;
		/// Number of samples for null hypothesis
		std::vector<UInt> n_null_;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H
