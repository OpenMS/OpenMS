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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <hash_map.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which uses a isotopic wavelet to find seeds.
	
		This seeder select interesting regions in the map using a wavelet funtion
		modelling the isotopic distribution.
		
		@todo Write test with more than one scan and test not only the intensity (Ole)
		@todo Document parameters (Ole)
		
		@ingroup FeatureFinder
	*/ 
  class IsotopeWaveletSeeder 
    : public BaseSeeder
  {

  public:
  
  	enum DimensionId
    {
      RT = DimensionDescription < LCMS_Tag >::RT,
      MZ = DimensionDescription < LCMS_Tag >::MZ
    };

    typedef FeaFiTraits::IntensityType IntensityType;
    typedef FeaFiTraits::CoordinateType CoordinateType;
    typedef KernelTraits::ProbabilityType ProbabilityType;
		typedef FeaFiTraits::MapType MapType;
		typedef MapType::PeakType PeakType;
		typedef MapType::SpectrumType SpectrumType;
		
		/// The mother wavelets (one for each charge state)
		typedef std::vector<std::vector<double> > WaveletCollection; 
		/// The Hash entry, stores pairs of scan number and list of charge scores
		typedef std::pair<std::list<UnsignedInt>, std::list<double> > DoubleList;
		/// The Hash. Maps mass bins to scans and charge scores
		typedef hash_multimap<UnsignedInt, DoubleList> SweepLineHash;
		/// Stores the charge states examined
		typedef std::list<UnsignedInt> ChargeVector;

    /// Default constructor
    IsotopeWaveletSeeder();

    /// destructor 
    virtual ~IsotopeWaveletSeeder();

    /// Copy constructor
    IsotopeWaveletSeeder(const IsotopeWaveletSeeder& rhs);
    
    /// Assignment operator
    IsotopeWaveletSeeder& operator= (const IsotopeWaveletSeeder& rhs);

    /// return next seed 
    IndexSet nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new IsotopeWaveletSeeder();
    }

    static const String getProductName()
    {
      return "IsotopeWaveletSeeder";
    }
		
  protected:
  	virtual void updateMembers_();
  	
		/// Computes m/z spacing of the LC-MS map
		void computeSpacings_();
		
		/// Precompute and store the gamma function (for the mother wavelet)
		void generateGammaValues_();
		
		/**
			@brief Computes the wavelet transform for several charges in nearly the same time.
			
			The working horse of the discrete-time continuous wavelet transform, but for several charges at the same time.
		 	Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 	this yields the same.		
		*/																			
		void fastMultiCorrelate(const SpectrumType& signal, std::vector<DPeakArray<1, PeakType > >* pwts, std::vector<double>* wt_thresholds);
		
		/** 
				@brief Returns the lamba parameter of the mother wavelet
				
				The lambda parameter essentially influences the shape of the wavelet.
		 		Since isotope patterns depend on mass, the wavelet has to adapt its shape. 
		 		For more insights look at the formula of the wavelet function. 
		*/
		inline CoordinateType getLambda (const CoordinateType realMass) const 
		{	
			return (0.035 + 0.000678*realMass); 
		}																	
		
		/// The wavelet (mother) function. 
 		inline double phiRaw (double t, double lambda, double a) throw ()
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
		
		/// Assigns scores to each charge state of a isotopic pattern
		void identifyCharge(const std::vector<DPeakArray<1, PeakType > >& candidates, std::vector<double>* wt_thresholds, UnsignedInt scan);
		
		/// Returns the interpolated value 
		inline double getInterpolatedValue (double x0, double x, double x1, double f0, double f1) const throw ()
		{
			return (f0 + (f1-f0)/(x1-x0) * (x-x0));
		}

		/// Returns a bucket containing the mass/charge @p mz
		inline std::pair<int, int> getNearBys(UnsignedInt scan, double mz, UnsignedInt start=0)
		{
			for (unsigned int i=start; i<traits_->getData()[scan].getContainer().size(); ++i)
			{
				if (traits_->getData()[scan].getContainer()[i].getPos() > mz)
				{
					return std::make_pair(i-1, i);
				}
			}

			//not found
			return std::make_pair(-1,-1);
		}	
		
		/// Returns the absolute mean of the intensities in @p signal
		double getAbsMean (const DPeakArray<1, PeakType >& signal,
																	UnsignedInt startIndex, 
																	UnsignedInt endIndex) const;
		
		/// Removes entries from hash occuring in less then rt_votes_cutoff_ scans
		void filterHashByRTVotes ();
		
	  /// Sums the intensities in adjacent scans
	  void sumUp_(SpectrumType& scan, UnsignedInt current_scan_index);
		
		///Aligns the two scans and increases intensities of peaks in @p scan if those peaks are present in @p neighbour
		void AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour);
										
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
		/// Determines distance of left box frame from monoisotopic bin
		CoordinateType mass_tolerance_right_;
		/// Determines distance of right box frame from monoisotopic bin
		CoordinateType mass_tolerance_left_;
		/// Tolerance for scan alignment
		CoordinateType tolerance_scansum_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELETSEEDER_H
