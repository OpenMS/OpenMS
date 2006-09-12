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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFINDERL_H

#include <iostream>
#include <fstream>
#include <strstream>
#include <list>
#include <hash_map.h>
#include <map.h>
#include <math.h>
#include <values.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
// #include <fftw3.h>
// #include <OpenMS/FORMAT/MzXMLFile.h>
// #include <OpenMS/FORMAT/DTA2DFile.h>
// #include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMedian.h>
// #include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>

#ifndef PRECISION_CUT_OFF
#define PRECISION_CUT_OFF 1e-10
#endif 

#ifndef LAPLACE_SMOOTH_EPSILON
#define LAPLACE_SMOOTH_EPSILON 1e-6
#endif

namespace OpenMS
{

/**
 * @brief A class to identify isotope Patterns in 2D mass spectra.
 * */
template <typename MapType>
class IsotopeFinder
{
	public:
				
		/** A simple matrix-like structure to "collect" several sampled wavelet functions in a single and simple 
	 	* data structure. */
		typedef std::vector<std::vector<double> > WaveletCollection; //Ugly, but simple and suited for this purposes 
		typedef std::vector<std::vector<double> > Matrix;

		typedef std::pair<std::list<double>, std::list<double> > DoubleList;
		typedef hash_multimap<unsigned int, DoubleList> SweepLineHash;
		typedef std::multimap<unsigned int, DoubleList> SweepLineMap;	

		/** The standard constructor. */
		IsotopeFinder () throw ();

		/** Extended constructor.
		 * @param experiment An experimental 2D mass spectrum
		 * @param integration_workspace A GSL (GNU scientific library) parameter for integration. There is no need to change
		 * this parameter.
		 * @param integration_epsilon A GSL (GNU scientific library) parameter for integration. There is no need to change
		 * this parameter.
		 * @param peak_cut_off The number of isotope peaks a wavelet should contain. Since the wavelet function has been
		 * tuned to resemble specific peak probabilies (dependent on the mass and a binomial distribution e.g.) you should
		 * usually not change this parameter. */
		IsotopeFinder (const MapType& experiment, const unsigned int integration_workspace=100, 
			double integration_epsilon=1e-6, unsigned int peak_cut_off=5, double score_cut_off=1e6, unsigned int 
			rt_votes_cut_off=1, double wtCutOff_=0) throw ();

		/* virtual*/ void initializeMe () throw ();

		/** The destructor */
		virtual ~IsotopeFinder () throw ();

		/** Reads in a DTA2D file. There is already a class named DTA2D, implemented exactly for the same purpose. 
		 * Unfortunately, reading-in a file by this class takes much longer than by this simple function here. */  		
		virtual void readTabFile (const string& filename) throw (Exception::FileNotFound);

		/** Computes a discrete-time continuous wavelet transform. 
		 * @param scanNumber The scan of your 2D map you are interested in (by index, not in RT).
		 * @param charge The charge of the pattern you are searching for. */
		virtual MSSpectrum<DRawDataPoint<2> > cwt (const unsigned int scanNumber, const unsigned int charge) throw ();

		/** Essentially the same function as cwt.
		 * Computes the wavelet transform for several charges in nearly the same time. */
		virtual std::vector<MSSpectrum<DRawDataPoint<2> > > cwtMulti (const unsigned int scanNumber, 
			const std::list<unsigned int>& charges) throw ();

		virtual void identifyCharge (std::vector<MSSpectrum<DRawDataPoint<2> > >& candidates) throw ();

		virtual SweepLineHash findFeatures (unsigned int startScan, unsigned int endScan, bool sweepLine=false) throw ();	

		
		/** The integration_workspace is a parameter of the GNU scientifc library (GSL). For more information see GSL's 
		 * documentation on integration. */
		inline unsigned int getIntegrationWorkSpace () const throw ()
			{ return (INTEGRATION_WORKSPACE); }
				
		/** The integration_workspace is a parameter of the GNU scientifc library (GSL). For more information see GSL's 
		 * documentation on integration. */
		inline void setIntegrationWorkSpace (const unsigned int integrationWorkspace) throw ()
			{ INTEGRATION_WORKSPACE = integrationWorkspace; }
		
		/** The integration_epsilon is a parameter of the GNU scientifc library (GSL). For more information see GSL's 
		 * documentation on integration. */
		inline double getIntegrationEpsilon () const throw ()
			{ return (INTEGRATION_EPSILON); }
		
		/** The integration_epsilon is a parameter of the GNU scientifc library (GSL). For more information see GSL's 
		 * documentation on integration. */
		inline void setIntegrationEpsilon (const double integrationEpsilon) throw ()
			{ INTEGRATION_EPSILON = integrationEpsilon; }

		/** The peak_cut parameter determines the number of isotope peaks a wavelet should contain. 
		 * Since the wavelet function has been tuned to resemble specific peak probabilies 
		 * (dependent on the mass and a binomial distribution e.g.) you should usually not change this parameter. */
		inline unsigned int getPeakCutOff () const throw ()
			{ return (PEAK_CUT_OFF); }
		
		/** The peak_cut parameter determines the number of isotope peaks a wavelet should contain. 
		 * Since the wavelet function has been tuned to resemble specific peak probabilies 
		 * (dependent on the mass and a binomial distribution e.g.) you should usually not change this parameter. */
		inline void setPeakCutOff (const unsigned int peakCutOff) throw ()
			{ PEAK_CUT_OFF = peakCutOff; }

		inline void setWtCutOff (const double wtCutOff) throw ()
			{ WT_CUT_OFF = wtCutOff; }		
		
		inline double getWtCutOff () throw ()
			{ return(WT_CUT_OFF); }

		inline double getScoreCutOff () const throw ()
			{ return (SCORE_CUT_OFF); }
		
		inline void setScoreCutOff (const double peakScoreOff) throw ()
			{ SCORE_CUT_OFF = peakScoreOff; }

		inline unsigned int getRTVotesCutOFF () const throw ()
			{ return (RT_VOTES_CUT_OFF); }
		
		inline void setRTVotesCutOff (const unsigned int RTVotesCutOff) throw ()
			{ RT_VOTES_CUT_OFF = RTVotesCutOff; }

		inline unsigned int getNumScans () const throw ()
			{ return (experiment_.size()); }	
			
		inline double getAvMZSpacing() const 
			{ return avMZSpacing_; }
		
		
		/** The common index operator. Returns the scan "index" (by index, not by RT).
		 * @todo Do we have to throw something like e.g. OutOfBoundsError here? */
		virtual inline MSSpectrum<DRawDataPoint<2> > operator[] (unsigned int index) const throw ()
			{ return (experiment_[index]); }


	protected:

		/** The working horse of the discrete-time continuous wavelet transform. 
		 * Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 * this yields the same. */
		virtual MSSpectrum<DRawDataPoint<2> > fastCorrelate (const MSSpectrum<DRawDataPoint<2> >& scan, 
			const unsigned int charge) throw ();
		
		/** The working horse of the discrete-time continuous wavelet transform, but for several charges at the same time.
		 * Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 * this yields the same. */
		virtual std::vector<MSSpectrum<DRawDataPoint<2> > > fastMultiCorrelate (const MSSpectrum<DRawDataPoint<2> >& scan, 						const std::list<unsigned int>& charges) throw ();

		/** The wavelet (mother) function. Exaclty the same as phiRaw, but needed for use with the GNU scientific library (GSL).
 		* The function has to be static, since we need a template-free pointer to this function.
 		* Note the different semantics of gamma between C++ and R. */
		static double gslPhiRaw (double x, void* params) throw ();

		virtual inline double phiRaw (const double t, const double lambda, const double a) throw ()
		{	
			int x0, x1; double f0, f1, fi;
			x0 = (int) trunc ((t/a + 1)/min_spacing_);
			x1 = x0+1;
			
			if ((unsigned int) x1 >= preComputedGamma_.size())
			{ return (0); }; 
			
			f0 = preComputedGamma_[x0];
			f1 = preComputedGamma_[x1];
			fi = (f0 + (f1-f0)/((x1-x0)*min_spacing_) * ((t/a+1)-x0*min_spacing_));
			double res = sin(2*M_PI*t/a) * exp(-lambda) * (pow(lambda,t)) / fi;
			//return (res);
			
			return (res<PRECISION_CUT_OFF ? 0:res); 
		}

		/** Estimates the average spacing for the m/z dimension. 
		 * Used internally to compute the circular convolution. */
		virtual double averageMZSpacing (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal) const throw ();

		/** The lambda parameter essentially influences the shape of the wavelet.
		 * Since isotope patterns depend on mass, the wavelet has the adapt its shape. 
		 * For more insights look at the formula of the wavelet function. */
		virtual inline double getLambda (const double realMass) const throw ()
		{	return (0.035 + 0.000678*realMass); }

		/** The wavelet (mother) function 
		 * Note the different semantics of gamma between C++ and R. */
		
		/** The wavelet function with dilation parameter.
		 * Note that this function can only be used if the m/z direction is equally spaced. Otherwise you have to sample
		 * directly from phiRaw (what is also done internally). 
		 * @param lambda See member function getLambda.
		 * @param charge The charge you are interested in (will be converted into the dilation parameter of the wavelet).
		 * @param resolution The *constant* m/z resolution.*/
		virtual std::vector<double> phiA (const double lambda, const double charge, const double resolution) throw ();

		/** ############################ We should throw something like OutOfBounds e.g. */
		double getMean (const MSSpectrum<DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		double getAbsMean (const MSSpectrum<DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		/** ############################ We should throw something like OutOfBounds e.g. */
		double getUpShiftedMoment (const MSSpectrum<DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int moment, const unsigned int endIndex) throw ();

		double getWaveletPotential (const MSSpectrum<DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int charge) const throw ();
		
		/** ############################ We should throw something like OutOfBounds e.g. 
		 *	############################ and an error if signal has length 1.
		 * */
		double getSd (const MSSpectrum<DRawDataPoint<2> >& signal, const double mean, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		double getAbsSd (const MSSpectrum<DRawDataPoint<2> >& signal, const double mean, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		void generateGammaValues () throw ();

		virtual void filterHashByRTVotes () throw ();
						
		/** The exerimental 2D mass spectrum. */
		MapType experiment_;
		
		/** Internal parameters. See their respective get and set functions for documentation. */
		unsigned int INTEGRATION_WORKSPACE;
		double INTEGRATION_EPSILON;  
		unsigned int PEAK_CUT_OFF;
		double SCORE_CUT_OFF;
		unsigned int RT_VOTES_CUT_OFF;		
		double WT_CUT_OFF;

		unsigned int waveletLength_; //Will be initialized with -1 and therefore it should be a signed integer	
		double avMZSpacing_;
		double min_spacing_;
		WaveletCollection* minSpacingPhis_;	
		/** The hash map for the sweep line algorithm */
		SweepLineHash hash_;
		hash_map<unsigned int, double> preComputedGamma_;
		

}; //class

			

template <typename MapType>
IsotopeFinder<MapType>::IsotopeFinder () throw () : INTEGRATION_WORKSPACE(100), INTEGRATION_EPSILON (1e-6), 
	PEAK_CUT_OFF (5), SCORE_CUT_OFF (1e6), RT_VOTES_CUT_OFF (1), WT_CUT_OFF (0), waveletLength_ (0), 
	avMZSpacing_ (0), minSpacingPhis_ (NULL) 
{
}


template <typename MapType>
IsotopeFinder<MapType>::~IsotopeFinder () throw ()
{
}


template <typename MapType>
IsotopeFinder<MapType>::IsotopeFinder (const MapType& experiment, const unsigned int integration_workspace, 
	double integration_epsilon, unsigned int peak_cut_off, double score_cut_off, unsigned int rt_votes_cut_off, 
	double wtCutOff) throw () 
	: INTEGRATION_WORKSPACE(integration_workspace), INTEGRATION_EPSILON (integration_epsilon), PEAK_CUT_OFF (peak_cut_off), 
	SCORE_CUT_OFF (score_cut_off), RT_VOTES_CUT_OFF (rt_votes_cut_off), WT_CUT_OFF (wtCutOff), waveletLength_(0), 
	avMZSpacing_(0), minSpacingPhis_ (NULL) 

{
	experiment_ = experiment;
	initializeMe();	
}


template <typename MapType>
void IsotopeFinder<MapType>::initializeMe () throw ()
{	
	//Since the signal might be unequally spaced, we have to sample the wavelet function for each translational step 
	//First, we have to estimate the average spacing in M/Z direction; otherwise we could sample too many or too few points
	//for the wavelet function.
	experiment_.updateRanges();
	//experiment_.sortSpectra(true);
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > signal; experiment_.get2DData(signal); 
	avMZSpacing_ = averageMZSpacing (signal);
	waveletLength_ = (int) (PEAK_CUT_OFF/avMZSpacing_);
 	std::vector<double> help (waveletLength_, 0);		
	delete (minSpacingPhis_);	
	std::list<unsigned int> charges;
	charges.push_back(1); charges.push_back(2); charges.push_back(3); charges.push_back(4); 
	minSpacingPhis_ = new WaveletCollection (charges.size(), help);			
	min_spacing_ = INT_MAX; double tmp=0;
	for (unsigned int i=0; i<signal.size()-1; ++i)
	{
		tmp = signal[i+1].getPosition().Y() - signal[i].getPosition().Y();
		if (fabs(tmp) < min_spacing_)
			min_spacing_ = tmp;
	}
	generateGammaValues();
}


template <typename MapType>
void IsotopeFinder<MapType>::readTabFile (const string& filename) throw (Exception::FileNotFound)
{
	typename MapType::SpectrumType spec, spec2;
  typename MapType::SpectrumType::PeakType p;

	TopHatFilter filter;
	
	std::ifstream ifile (filename.c_str());
 	if (!ifile)
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);				
	}	
	double c_x, old_x=-1, c_y, c_intensity; 
	while (ifile >> c_x)
	{
		ifile >> c_y;
		ifile >> c_intensity;

		//i.e. a new spectrum begins
		if (c_x != old_x && !spec.empty())
		{
			spec.setRetentionTime(old_x);
			filter.filter (spec.begin(), spec.end(), spec2);
			spec2.setRetentionTime(old_x);
			experiment_.push_back(spec2);
			spec.clear();	spec2.clear();	
		}
		
		p.setIntensity (c_intensity);
		p.getPosition()[0] = c_y;
		p.getPosition()[1] = c_x; old_x = c_x;
		
		spec.push_back(p);
	}
	
	ifile.close();		
	initializeMe();	
}


template <typename MapType>
MSSpectrum<DRawDataPoint<2> > IsotopeFinder<MapType>::cwt (const unsigned int scanNumber, const unsigned int charge) 
	throw ()
{
	return (fastCorrelate (experiment_[scanNumber], charge)); 
}


template <typename MapType>
std::vector<MSSpectrum<DRawDataPoint<2> > > IsotopeFinder<MapType>::cwtMulti (const unsigned int scanNumber, 
	const std::list<unsigned int>& charges) throw ()
{
	return (fastMultiCorrelate (experiment_[scanNumber], charges));
}




template <typename MapType>
typename IsotopeFinder<MapType>::SweepLineHash IsotopeFinder<MapType>::findFeatures (unsigned int startScan, unsigned int endScan, const bool sweepLine) throw ()
{
	hash_.clear();			
	experiment_.updateRanges(); //this needs identifyCharge
	std::list<unsigned int> charges;
	charges.push_back(1); charges.push_back(2); charges.push_back(3); charges.push_back(4); 
	std::vector<MSSpectrum<DRawDataPoint<2> > > pwts;
	std::vector<double> RTs (experiment_.size());
	std::cout << std::endl;
	for (unsigned int i=startScan; i<=endScan; ++i)
	{
		std::cout << "Spectrum " << i << " of " << endScan << std::endl; 
		pwts = cwtMulti (i, charges);
		identifyCharge (pwts);
		RTs[i] = experiment_[i].getRetentionTime();
	}	
	
	
	if (sweepLine)
		filterHashByRTVotes ();	
	
	//Do not move this further up!
	SweepLineHash::iterator iter_h;
	std::list<double>::iterator iter_s; double c_sum=0;
	for (iter_h=hash_.begin(); iter_h != hash_.end(); ++iter_h)
	{
		DoubleList& c_entry = iter_h->second;
		c_sum=0;
		for (iter_s=c_entry.second.begin(); iter_s != c_entry.second.end(); ++iter_s)
			c_sum += *iter_s;
		for (iter_s=c_entry.second.begin(); iter_s != c_entry.second.end(); ++iter_s)
			*iter_s /= c_sum;
	}	

	for (SweepLineHash::const_iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
	{
		if (iter->first != 0)
			std::cout << experiment_.getMin().Y()+(iter->first-1)*avMZSpacing_ << " : " 
			<< experiment_.getMin().Y()+iter->first*avMZSpacing_ << "\t [ ";
		else
			std::cout << experiment_.getMin().Y() << " : " << experiment_.getMin().Y()+iter->first*avMZSpacing_ << "\t [ ";
	 	for (std::list<double>::const_iterator iter_cl2=iter->second.first.begin(); 
			iter_cl2 != iter->second.first.end(); ++iter_cl2)
				std::cout	<< *iter_cl2 << " "; std::cout << "]  { ";	
		for (std::list<double>::const_iterator iter_l=iter->second.second.begin();
			iter_l != iter->second.second.end(); ++iter_l)
				std::cout << *iter_l << " ";
		std::cout << "}" << std::endl;
	}
// 	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > data;
// 	experiment_.get2DData(data);
// 	std::cout << std::endl << "Found " << hash_.size() << " potential isotope pattern(s) (outof " 
// 		<< data.size() << " possible points)." << std::endl;	

	return (hash_);
}


template <typename MapType>
void IsotopeFinder<MapType>::filterHashByRTVotes () throw ()
{
	std::cout << std::endl << "Hash size before filtering: " << hash_.size() << std::endl;	
	SweepLineHash::iterator iter_f, iter_b;
 	std::vector<SweepLineHash::iterator> toDelete;
	std::list<double>::iterator iter_l;

	for (SweepLineHash::iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
		if (iter->second.first.size() < RT_VOTES_CUT_OFF)
			toDelete.push_back(iter);

	for (unsigned int i=0; i<toDelete.size(); ++i)
		hash_.erase (toDelete[i]);	

	toDelete.clear();
	bool delete_it=true;
	for (SweepLineHash::iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
	{
		for (iter_l=iter->second.second.begin(); iter_l!=iter->second.second.end(); ++iter_l)		
			if (*iter_l >= SCORE_CUT_OFF)
			{
				delete_it = false;
				break;
			}
		if (delete_it)
			toDelete.push_back(iter);
		else
			delete_it = true;
	}

	for (unsigned int i=0; i<toDelete.size(); ++i)
		hash_.erase (toDelete[i]);
	
	std::cout << "Hash size after filtering: " << hash_.size() << std::endl << std::endl;	
	
}

				
				
template <typename MapType>
/*std::pair<typename IsotopeFinder<MapType>::WaveletCollection, typename IsotopeFinder<MapType>::WaveletCollection>*/ 
void IsotopeFinder<MapType>::identifyCharge (std::vector<MSSpectrum<DRawDataPoint<2> > >& candidates) throw ()
{	
	std::vector<double> int_mins (candidates[0].size(),INT_MIN), zeros (candidates[0].size(),0);
	WaveletCollection meansC (candidates.size(), int_mins);
	WaveletCollection sdsC (candidates.size(), zeros);
	WaveletCollection scoresC (candidates.size(), zeros);

	std::vector<unsigned int> start_indices, end_indices;
	//In order to determine the start and end indices, we first need to know the width of the region one should consider 
	//to estimate the mean and the sd of the pattern candidate. 
	//That region is defined by the position of the heighst amplitude +/- waveletLength_.
	
	typedef MSSpectrum<DRawDataPoint<2> >::ContainerType containerType;
	containerType::iterator iter; 
	unsigned int start_index, end_index, c_index, i_iter; //Helping variables
	double abs_mean_threshold, sd, seed_mz; 
	std::vector<bool> processed (candidates[0].size(), false);
	//double c_snr=0;
	DSignalToNoiseEstimatorMedian<1, containerType::iterator> sTn;
	for (unsigned int c=0; c<candidates.size(); ++c)		
	{
		processed = std::vector<bool> (candidates[0].size(), false); //Reset
		candidates[c].updateRanges();
		containerType c_candidate = candidates[c].getContainer();  
		//Ugly, but do not how to do this in a better (and easy) way
		for (unsigned int i=0; i<c_candidate.size(); ++i)
			c_candidate[i].setPosition(DPosition<2>(c_candidate[i].getPosition().X(), i));
			
		sort (c_candidate.begin(), 
					c_candidate.end(), 
					ReverseComparator< DRawDataPoint<2>::IntensityLess >::ReverseComparator() );
		
		//sort (c_candidate.begin(), c_candidate.end(), comparator); 
	
		i_iter=0;
		for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter, ++i_iter)
		{
			c_index = (int) (iter->getPosition().Y());			
			start_index = c_index-waveletLength_+1;
			end_index = c_index+waveletLength_-1;		

			//Catch impossible cases
			if (end_index >= candidates[c].size() || start_index > end_index)
				continue;			

			abs_mean_threshold = getAbsMean (candidates[c], start_index, end_index);
			sd = getAbsSd (candidates[c], abs_mean_threshold, start_index, end_index);

			//sTn.init (candidates[c].begin()+start_index, candidates[c].begin()+end_index);
			//c_snr = sTn.getSignalToNoise (candidates[c].begin()+c_index);
			//std::cout << c_snr << std::endl; 
			
			if (iter->getIntensity() <= abs_mean_threshold+2*sd || iter->getIntensity() < WT_CUT_OFF || processed[c_index])	
				continue;

			/*if (c_snr < 3)
			{
				std::cout << "SNR aborted! " << std::endl;
				continue;
			}*/

	
			//Mark as processed
			for (unsigned int z=start_index; z<=end_index; ++z)
				processed[z] = true;

			seed_mz=iter->getPosition().X();	
			
			meansC[c][c_index] = getUpShiftedMoment (candidates[c], 1, start_index, end_index)
 				/ (fabs(getMean (candidates[c], start_index, end_index))+LAPLACE_SMOOTH_EPSILON);
			sdsC[c][c_index] = getUpShiftedMoment (candidates[c], 2, start_index, end_index)
				/ (getUpShiftedMoment (candidates[c], 3, start_index, end_index)+LAPLACE_SMOOTH_EPSILON);
			scoresC[c][c_index] = meansC[c][c_index]*sdsC[c][c_index]*iter->getIntensity();
		}	
	}

	std::vector<int> positions (candidates.size(), 0); std::list<double> c_list;
	double c_mz, c_fill_mz=experiment_.getMin().Y()+avMZSpacing_; //The first right boundary of the hash cell
		
	unsigned int c_fill_index=0, count_finished_charges=0;
	DoubleList c_pair; std::list<double>::iterator iter_cl, iter_cl_hash; 
	double allZero; std::list<double> c_fill_list; std::pair<SweepLineHash::iterator, SweepLineHash::iterator> iter_hash;
	while (1) // ;-)
	{	
		//Termination criterion
		for (unsigned int c=0; c<candidates.size(); ++c)
		{
			if (positions[c] >= (int) candidates[c].size())
			{
				if (++count_finished_charges >= candidates.size())
					goto FINISHED_HASHING;
				positions[c] = -1;
			}
		}
						
		for (unsigned int c=0; c<candidates.size(); ++c)		
		{
			if (positions[c] < 0 || positions[c] >= (int) candidates[c].size())
				continue;
			c_mz = candidates[c].getContainer()[positions[c]].getPos();				
			while (c_mz <= c_fill_mz)
			{						
				c_list.push_back (scoresC[c][positions[c]]);
				
				if (++(positions[c]) >= (int) candidates[c].size())
					break;
				c_mz = candidates[c].getContainer()[positions[c]].getPos();				
			}
		}				

		allZero=true;	
		for (iter_cl=c_list.begin(); iter_cl!=c_list.end(); ++iter_cl)
			if (*iter_cl != 0)
				allZero=false;
		
		if (!c_list.empty() && !allZero)
		{
			iter_hash=hash_.equal_range(c_fill_index);
			while (iter_hash.first != iter_hash.second)		
			{
				//std::cout << "M/Z already in hash table ... " << std::endl;
				typename MapType::const_iterator beforeRT = experiment_.RTBegin(candidates[0].getRetentionTime());
				if (find(iter_hash.first->second.first.begin(), iter_hash.first->second.first.end(), 
					(--beforeRT)->getRetentionTime()) == iter_hash.first->second.first.end())
				{
					//i.e. there is no neighbouring entry before this retention time
					//i.e. we can treat this case as if no entry is present in the hash
					++iter_hash.first;
					continue;
				}
				c_fill_list = iter_hash.first->second.first; c_fill_list.push_back(candidates[0].getRetentionTime());	
				for (iter_cl = c_list.begin(), iter_cl_hash = iter_hash.first->second.second.begin(); iter_cl != c_list.end(); 
					++iter_cl, ++iter_cl_hash)
						*iter_cl += *iter_cl_hash;
				hash_.erase (iter_hash.first);
			  c_pair = DoubleList (c_fill_list, c_list);
				goto FINISH;
			}
			
			//std::cout << "New M/Z entry for hash table ..." << std::endl;
			//Store RT and the corresponding score
			c_fill_list.clear(); 
			c_fill_list.push_back(candidates[0].getRetentionTime());	
			c_pair = DoubleList (c_fill_list, c_list);

		FINISH:			
			hash_.insert (SweepLineHash::value_type(c_fill_index, c_pair));
		}

		c_list.clear();
	  ++c_fill_index; 
		c_fill_mz += avMZSpacing_;
	}

 FINISHED_HASHING:
	
	return;	
}	


template <typename MapType>
double  IsotopeFinder<MapType>::getWaveletPotential (const MSSpectrum<DRawDataPoint<2> >& signal, 
	const unsigned int startIndex, const unsigned int charge) const throw ()
{	
	double realMass = signal.getContainer()[startIndex].getPos() * charge;
	double lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet

	std::vector<double> phi (waveletLength_); double cumSpacing=0, cSpacing=0;
	for (unsigned int j=0; j<waveletLength_; ++j)
	{	
		phi[j] = phiRaw (cumSpacing, lambda, 1.0/(double)charge);
		cSpacing = signal.getContainer()[(startIndex+j+1)%signal.size()].getPos() 
			- signal.getContainer()[(startIndex+j)%signal.size()].getPos();
		//cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal.size()-1).
		//Since this case is only of theoretical interest (we do not expect any important signals at the very end of the 
		//spectrum), we simlply use the average spacing in that case.
		if (cSpacing < 0)
			cumSpacing += avMZSpacing_;
		else //the "normal" case
			cumSpacing += signal.getContainer()[(startIndex+j+1)%signal.size()].getPos() 
				- signal.getContainer()[(startIndex+j)%signal.size()].getPos() ; 
	}		

	double min=INT_MAX;
	for (unsigned int i=0; i<phi.size(); ++i)
		if (signal.getContainer()[i].getIntensity() < min)
				min = signal.getContainer()[i].getIntensity();
	min = fabs(min);
	
	double res=0, mean=0;
			
	for (unsigned int i=0; i<phi.size(); ++i)
		res += signal.getContainer()[i].getIntensity() + min;
	mean = (res/(double)(phi.size()));
	
	return (mean);
}


template <typename MapType>
double IsotopeFinder<MapType>::getMean (const MSSpectrum<DRawDataPoint<2> >& signal, 
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += signal.getContainer()[i].getIntensity();

	return (res/(double)(endIndex-startIndex+1));	
}

				
template <typename MapType>
double IsotopeFinder<MapType>::getAbsMean (const MSSpectrum<DRawDataPoint<2> >& signal, 
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
			res += fabs(signal.getContainer()[i].getIntensity());

	return (res/(double)(endIndex-startIndex+1));	
}


template <typename MapType>
double IsotopeFinder<MapType>::getUpShiftedMoment (const MSSpectrum<DRawDataPoint<2> >& signal, const unsigned int moment, 
	const unsigned int startIndex, const unsigned int endIndex) throw ()
{
	double tmp=0, min=INT_MAX; 
	for (unsigned int i=startIndex; i<endIndex; ++i)
		if (signal.getContainer()[i].getIntensity() < min)
			min = signal.getContainer()[i].getIntensity();
	
	if (min < 0)
		min *=-1;
	else min=0;
	
	double mean=0, res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += signal.getContainer()[i].getIntensity() + min;
	mean = (res/(double)(endIndex-startIndex+1));

	if (moment == 1)
		return (mean);

	if (moment == 2) //I.e. variance
	{
		res=0;
		for (unsigned int i=startIndex; i<endIndex; ++i)
		{
			tmp = signal.getContainer()[i].getIntensity() + min - mean;
			res += tmp*tmp;
		}
		return (sqrt(res/(double)(endIndex-startIndex)));
	}
	
	//I.e. skewness
	
	res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
	{
		tmp = signal.getContainer()[i].getIntensity() + min - mean;
		res += tmp*tmp*tmp;
	}

	return ((pow(fabs(res/(double)(endIndex-startIndex+1)), 1.0/3.0)));
}


template <typename MapType>
double IsotopeFinder<MapType>::getSd (const MSSpectrum<DRawDataPoint<2> >& signal, const double mean, 
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += (signal.getContainer()[i].getIntensity()-mean)*(signal.getContainer()[i].getIntensity()-mean);

	return (sqrt(res/(double)(endIndex-startIndex)));	
}


template <typename MapType>
double IsotopeFinder<MapType>::getAbsSd (const MSSpectrum<DRawDataPoint<2> >& signal, const double mean, 
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += (fabs(signal.getContainer()[i].getIntensity())-mean)*(fabs(signal.getContainer()[i].getIntensity())-mean);

	return (sqrt(res/(double)(endIndex-startIndex)));	
}



template <typename MapType>
MSSpectrum<DRawDataPoint<2> > IsotopeFinder<MapType>::fastCorrelate (const MSSpectrum<DRawDataPoint<2> >& scan, 
			const unsigned int charge) throw ()
{
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > res (scan.getContainer());
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > signal (scan.getContainer());
	double sum=0; unsigned int k=0;		
	
	std::vector<double> phi (waveletLength_); //the wavelet
	
	double cumSpacing=0, cSpacing=0, realMass=0, lambda=0; //helping variables
	for (unsigned int i=0; i<signal.size(); ++i)
	{		
		//Now, let's sample the wavelet
		cumSpacing=0;			
		realMass = signal[i].getPos() * charge;
		lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet

		for (unsigned int j=0; j<waveletLength_; ++j)
		{	
			phi[j] = phiRaw (cumSpacing, lambda, 1.0/(double)charge);
			cSpacing = signal[(i+j+1)%signal.size()].getPos() - signal[(i+j)%signal.size()].getPos();
			//cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal.size()-1).
			//Since this case is only of theoretical interest (we do not expect any important signals at the very end of the 
			//spectrum), we simlply use the average spacing in that case.
			if (cSpacing < 0)
				cumSpacing += avMZSpacing_;
			else //the "normal" case
				cumSpacing += signal[(i+j+1)%signal.size()].getPos() - signal[(i+j)%signal.size()].getPos() ; 
		}			
					
		sum=0;
		k=0;
		for (unsigned int j=i; j<signal.size() && k<phi.size(); ++j, ++k)
			sum += signal[j].getIntensity()*phi[k];

		for (unsigned int l=0; l<i && k<phi.size(); ++l, ++k)
			sum += signal[l].getIntensity()*phi[k];

		res[i].setIntensity(sum);
	}

	MSSpectrum<DRawDataPoint<2> > result (scan);
 	result.setContainer (res);
	
	return (result);
}


template <typename MapType>
std::vector<MSSpectrum<DRawDataPoint<2> > > IsotopeFinder<MapType>::fastMultiCorrelate 
	(const MSSpectrum<DRawDataPoint<2> >& scan, const std::list<unsigned int>& charges) throw ()
{					
	std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > > res (charges.size(), scan.getContainer());
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > signal (scan.getContainer());
	
	WaveletCollection phis (charges.size(), std::vector<double> (waveletLength_)); //all necessary wavelets (by rows)
	
	double cumSpacing=0, cSpacing=0, realMass=0, lambda=0; //helping variables
	std::list<unsigned int>::const_iterator charge_iter; unsigned int k=0; //helping variables

	for (unsigned int i=0; i<signal.size(); ++i)
	{		
		//Now, let's sample the wavelets
		for (charge_iter=charges.begin(), k=0; charge_iter!=charges.end(); ++charge_iter, ++k)
		{
			cumSpacing=0;				
			realMass = signal[i].getPos() * (*charge_iter);
			lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
			for (unsigned int j=0; j<waveletLength_; ++j)
			{	
				realMass = signal[(i+j+1)%signal.size()].getPos() * (*charge_iter);
				lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
				phis[k][j] = phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter));
				cSpacing = signal[(i+j+1)%signal.size()].getPos() - signal[(i+j)%signal.size()].getPos();
				//cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal.size()-1).
				//Since this case is only of theoretical interest (we do not expect any important signals at the very end of the 
				//spectrum), we simlply use the average spacing in that case.
				if (cSpacing < 0)
					cumSpacing += avMZSpacing_;
				else //The "normal" case
					cumSpacing += signal[(i+j+1)%signal.size()].getPos() - signal[(i+j)%signal.size()].getPos(); 
			}
		}			
					
		std::vector<double> sums (charges.size());
		k=0;
		for (unsigned int j=i; j<signal.size() && k<phis[0].size(); ++j, ++k)
		{	//Since all wavelet functions have the same length, we can simply use phis[0].size()
			
			for (unsigned int m=0; m<charges.size(); ++m)
				sums[m] += signal[j].getIntensity()*phis[m][k];
		}

		for (unsigned int l=0; l<i && k<phis[0].size(); ++l, ++k)
		{	//Since all wavelet functions have the same length, we can simply use phis[0].size()
	
			for (unsigned int m=0; m<charges.size(); ++m)
				sums[m] += signal[l].getIntensity()*phis[m][k];
		}

		//Store the current convolution result
		for (unsigned int m=0; m<charges.size(); ++m)
			res[m][i].setIntensity(sums[m]);
	}

	std::vector<MSSpectrum<DRawDataPoint<2> > > results (charges.size());
	MSSpectrum<DRawDataPoint<2> > result (scan);
	for (unsigned int i=0; i<charges.size(); ++i)
	{
 		result.setContainer (res[i]);
		result.setRetentionTime(scan.getRetentionTime());
		results[i] = result;
	}

	return (results);
}


template <typename MapType>
double IsotopeFinder<MapType>::averageMZSpacing (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal) 
	const throw ()
{
	double MZspacing=0, avMZspacing=0, leftOuts=0;
	for (unsigned int i=0; i<signal.size()-1; ++i)
	{			
		MZspacing = signal[i+1].getPosition().Y() - signal[i].getPosition().Y();
		if (signal[i+1].getPosition().X() != signal[i].getPosition().X()) //i.e. a new scan begins
		{
			++leftOuts;
			continue;			
		}
		avMZspacing += MZspacing;
	}
	avMZspacing /= (double)(signal.size()-leftOuts-1); //-1, since there are n data points and hence n-1 spacings
		
	return (avMZspacing);		
}


template <typename MapType>
double IsotopeFinder<MapType>::gslPhiRaw (double x, void* params) throw ()
{
	double* parameters = reinterpret_cast<double*> (params);
	double lambda = parameters[0];
	double a = parameters[1];
		
	return (sin(2*M_PI*x/a) * exp(-lambda) * (pow(lambda,x)) / tgamma((x/a)+1));
}


template <typename MapType>
void IsotopeFinder<MapType>::generateGammaValues () throw ()
{
	double query = 0; unsigned int counter=0;
	while (query <= 4*PEAK_CUT_OFF +1) //4 because of max_charge	
	{
		preComputedGamma_[counter] = tgamma (query);
		query += min_spacing_;			
		++counter;
	}
}
				

template <typename MapType>
std::vector<double> IsotopeFinder<MapType>::phiA (const double lambda, const double charge, const double resolution) 
	throw ()
{
	double a = 1.0/charge; //Convert the charge into the dilation coefficient of the wavelet
	std::vector<double > res (waveletLength_);	
	
	gsl_integration_workspace* GSL_WSP = gsl_integration_workspace_alloc (INTEGRATION_WORKSPACE);
	
	double iValue, absError;
	gsl_function F;
	F.function = &IsotopeFinder<MapType>::gslPhiRaw;
	double* parameters = new double [2];
	parameters[0] = lambda;
	parameters[1] = a;
	F.params = parameters;
	
	gsl_integration_qag (&F, 0, waveletLength_, 2*INTEGRATION_EPSILON, INTEGRATION_EPSILON, INTEGRATION_WORKSPACE-1, 
		GSL_INTEG_GAUSS15, GSL_WSP, &iValue, &absError);  
	gsl_integration_workspace_free (GSL_WSP);
	
	double mean = iValue / waveletLength_;	
	for (unsigned int i=0; i<waveletLength_; ++i)
		res[i] = phiRaw (i*resolution, lambda, a) - mean;

	delete[] (parameters);
	
	return (res);				
}

#endif 


} //namespace
