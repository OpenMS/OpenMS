// -*- C++: make; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff, Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

// GSL includes
#include <gsl/gsl_sf_gamma.h>
// #include <gsl/gsl_integration.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <hash_map.h>
#include <map.h>
#include <math.h>
#include <values.h>
#include <algorithm>
#include <vector>


namespace OpenMS
{

	/** @brief Seeding module for the peptide quantification algorithm in OpenMS.
	
			This seeder select interesting regions in the map using a wavelet funtion
			modelling the isotopic distribution.
				
			@ingroup FeatureFinder
		
	*/ 
  class IsotopeWaveletSeeder 
    : public BaseSeeder
  {

  public:
  
  	enum DimensionId
    {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    typedef FeaFiTraits::IntensityType IntensityType;
    typedef FeaFiTraits::CoordinateType CoordinateType;
    typedef KernelTraits::ProbabilityType ProbabilityType;
		typedef FeaFiTraits::MapType MapType;
		typedef MapType::PIterator PeakIterator;
		typedef MapType::PeakType PeakType;
		
		/// This comparator is used to compare
		typedef PeakType::NthPositionLess<0> MZless;
		
		/// The mother wavelets (one for each charge state)
		typedef std::vector<std::vector<double> > WaveletCollection; 
		///
		typedef std::pair<std::list<double>, std::list<double> > DoubleList;
		///
		typedef hash_multimap<unsigned int, DoubleList> SweepLineHash;
		///
		typedef std::multimap<unsigned int, DoubleList> SweepLineMap;	
		///
		typedef std::list<unsigned int> ChargeVector;

    /// standard constructor
    IsotopeWaveletSeeder();

    /// destructor 
    virtual ~IsotopeWaveletSeeder();

    /// return next seed 
    IndexSet nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new IsotopeWaveletSeeder();
    }

    static const String getName()
    {
      return "IsotopeWaveletSeeder";
    }
		
  protected:
  
		/// Computes m/z spacing of the LC-MS map
		void computeSpacings_();
		
		/// Precompute and store the gamma function
		void generateGammaValues_();
		
		/**
			@brief Computes the wavelet transform for several charges in nearly the same time.
			
			The working horse of the discrete-time continuous wavelet transform, but for several charges at the same time.
		 Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 this yields the same.		
		**/																			
		void fastMultiCorrelate(const DPeakArray<1, PeakType >& signal, 
	 																	std::vector<DPeakArray<1, PeakType > >* pwts, 
																		std::vector<double>* wt_thresholds);
		
		/** The lambda parameter essentially influences the shape of the wavelet.
		 		Since isotope patterns depend on mass, the wavelet has to adapt its shape. 
		 		For more insights look at the formula of the wavelet function. 
		**/
		inline CoordinateType getLambda (const CoordinateType realMass) const {	return (0.035 + 0.000678*realMass); }																	
		
		/// The wavelet (mother) function. 
 		inline double phiRaw (const double t, const double lambda, const double a) throw ()
		{	
			if (t>2*peak_cut_off_)	return(0);
			
			int x0, x1; double f0, f1, fi, res=0;
			x0 = (int) trunc ((t/a + 1)/min_spacing_);
			x1 = x0+1;
			if ((unsigned int) x1 < preComputedGamma_.size())
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
		
		void identifyCharge (const std::vector<DPeakArray<1, PeakType > >& candidates, 
																std::vector<double>* wt_thresholds, 
																const UnsignedInt scan, 
																const CoordinateType current_rt);
		
		/// Returns the interpolated value TODO: Check who is calling this function and where
		inline double getInterpolatedValue (const double x0, const double x, const double x1, 
			const double f0, const double f1) const throw ()
		{
			return (f0 + (f1-f0)/(x1-x0) * (x-x0));
		}

		/// Returns a bucket containing the mass/charge @p mz
		inline std::pair<int, int> getNearBys (const unsigned int scan, const double mz, const unsigned int start=0)
		{
				for (unsigned int i=start; i<traits_->getData()[scan].getContainer().size(); ++i)
				{
						if (traits_->getData()[scan].getContainer()[i].getPos() > mz)
							return (std::pair<int, int> (i-1, i));
				}

				//not found
				return (std::pair<int, int> (-1, -1));
		}	
		
		/// Returns the absolute mean of the intensities in this scan
		double getAbsMean (const DPeakArray<1, PeakType >& signal,
																	const unsigned int startIndex, 
																	const unsigned int endIndex) const;
		
		/// Removes patterns from hash occuring in less then rt_votes_cutoff_ scans
		void filterHashByRTVotes ();
										
		/// Does not need much explanation.
		bool is_initialized_;
		/// Number of isotopic peaks a wavelet should contain		
		UnsignedInt peak_cut_off_;
		/// Length of the mother wavelet
		UnsignedInt waveletLength_;
		/// Average spacing in a MS scan
		CoordinateType avMZSpacing_;
		/// Minium spacing
		CoordinateType min_spacing_;
		/// Minimum number of scans in which an isotopic pattern has to occur 
		UnsignedInt rt_votes_cutoff_;
		/// Charge states being tested
		ChargeVector charges_;
		/// Maximum charge state tested
		UnsignedInt max_charge_;
		//// Hash storing the detected regions
		SweepLineHash hash_;
		/// Stores the Gamme function
		hash_map<UnsignedInt, double> preComputedGamma_;
		/// Iterator pointing at the next region
		SweepLineHash::const_iterator hash_iter_;
		/// Determines threshold for the minimum score of a peak
		IntensityType intensity_factor_;
		/// Determines threshold for cwt of a peak
		IntensityType avg_intensity_factor_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H
