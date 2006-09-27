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
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFINDER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <hash_map.h>
#include <map.h>
#include <math.h>
#include <values.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMedian.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>

#ifndef LAPLACE_SMOOTH_EPSILON
#define LAPLACE_SMOOTH_EPSILON 1e-6
#endif

#ifndef DEBUG_ISOTOPE_FINDER
#define DEBUG_ISOTOPE_FINDER 0
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
			double integration_epsilon=1e-6, unsigned int peak_cut_off=5, double score_cut_off=0, 
			unsigned int rt_votes_cut_off=6	, double wtCutOff=0, double mz_interleave=3) throw ();

		virtual void initializeMe () throw ();

		/** The destructor */
		virtual ~IsotopeFinder () throw ();

		/** Reads in a DTA2D file. There is already a class named DTA2D, implemented exactly for the same purpose. 
		 * Unfortunately, reading-in a file by this class takes much longer than by this simple function here. */  		
		virtual void readTabFile (const std::string& filename) throw (Exception::FileNotFound);

		/** Essentially the same function as cwt.
		 * Computes the wavelet transform for several charges in nearly the same time. */
		virtual inline void cwtMulti (const unsigned int scanNumber, 
			const std::list<unsigned int>& charges, std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* pwts,
			std::vector<double>* wt_thresholds) throw ();

		virtual void identifyCharge (const std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >& candidates, 
			std::vector<double>* wt_thresholds, const unsigned int scan, const double RT) throw ();

		virtual SweepLineHash findFeatures (const unsigned int start_scan, const unsigned int end_scan, 
			const bool sweepLine=true) throw ();	

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
			{ PEAK_CUT_OFF = peakCutOff; initializeMe(); }

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

		/** The common index operator. Returns the scan "index" (by index, not by RT).
		 * @todo Do we have to throw something like e.g. OutOfBoundsError here? */
		virtual inline MSSpectrum<DRawDataPoint<2> > operator[] (unsigned int index) const throw ()
			{ return (experiment_[index]); }

		inline unsigned int getNumScans () const throw ()
    	{ return (experiment_.size()); }

		virtual void printMapEntry (SweepLineMap::iterator iter) throw ();

		virtual void printHashEntry (SweepLineHash::iterator iter) throw ();

		virtual void createGNUplot (unsigned int alignedTo, double mz, unsigned int charge,  
			const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >* signal, std::vector<double>* wavelet, 
			DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >* transform) throw ();
		
		inline unsigned int unique_merge (const std::list<double>& a, const std::list<double>& b,
			std::vector<double>::iterator& c_begin, std::vector<double>::iterator& c_end, std::vector<double>* res) throw ();
			
			double getAvMZSpacing() { return avMZSpacing_; }

	protected:

		/** The working horse of the discrete-time continuous wavelet transform, but for several charges at the same time.
		 * Note that you should compute a convolution instead of an correlation. Since we do not mirror the wavelet function
		 * this yields the same. */
		virtual void fastMultiCorrelate (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, 
			const std::list<unsigned int>& charges, std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* pwts,
			std::vector<double>* wt_thresholds) throw ();

		/** The wavelet (mother) function. Exaclty the same as phiRaw, but needed for use with the GNU scientific library (GSL).
 		* The function has to be static, since we need a template-free pointer to this function.
 		* Note the different semantics of gamma between C++ and R. */
		
		static double gslPhiRaw (double t, void* params) throw ();

		virtual inline double phiRaw (const double t, const double lambda, const double a) throw ()
		{	
			if (t>2*PEAK_CUT_OFF)	
				return(0);
			
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

		inline double getInterpolatedValue (const double x0, const double x, const double x1, 
			const double f0, const double f1) const throw ()
		{
			return (f0 + (f1-f0)/(x1-x0) * (x-x0));
		}

		inline std::pair<int, int> getNearBys (const unsigned int scan, const double mz, const unsigned int start=0)
		{
				for (unsigned int i=start; i<experiment_[scan].getContainer().size(); ++i)
				{
						if (experiment_[scan].getContainer()[i].getPos() > mz)
							return (std::pair<int, int> (i-1, i));
				}

				//not found
				return (std::pair<int, int> (-1, -1));
		}	

		virtual double phiRawInt (const double lambda, const double a) throw ();

		/** Estimates the average spacing for the m/z dimension. 
		 * Used internally to compute the circular convolution. */
		virtual double averageMZSpacing (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal) const throw ();
		virtual double averageRTSpacing (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal) const throw ();
		
		/** The lambda parameter essentially influences the shape of the wavelet.
		 * Since isotope patterns depend on mass, the wavelet has the adapt its shape. 
		 * For more insights look at the formula of the wavelet function. */
		virtual inline double getLambda (const double realMass) const throw ()
		{	return (0.035 + 0.000678*realMass); }

		virtual void filterHashByRTVotes () throw ();

		double getMean (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		double getAbsMean (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int endIndex) const throw ();

		double getUpShiftedMoment (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const unsigned int startIndex, 
			const unsigned int moment, const unsigned int endIndex) throw ();

		double getSd (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const double mean, 
			const unsigned int startIndex, const unsigned int endIndex) const throw ();

		double getAbsSd (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const double mean, 
			const unsigned int startIndex, const unsigned int endIndex) const throw ();

		void generateGammaValues () throw ();

		inline double getMZbyHashKey (const unsigned int key) const throw ()
		{ return (experiment_.getMin().Y()+key*avMZSpacing_); }

		virtual void prepareGNUplotFiles (const std::string& file="gnu.plot") throw ();

		/** The exerimental 2D mass spectrum. */
		MapType experiment_;
		
		/** Internal parameters. See their respective get and set functions for documentation. */
		unsigned int INTEGRATION_WORKSPACE;
		double INTEGRATION_EPSILON;  
		unsigned int PEAK_CUT_OFF;
		double SCORE_CUT_OFF;
		unsigned int RT_VOTES_CUT_OFF;
    double WT_CUT_OFF;
		unsigned int MZ_INTERLEAVE;
		
		unsigned int waveletLength_; //Will be initialized with -1 and therefore it should be a signed integer	
		double avMZSpacing_, avRTSpacing_, min_spacing_, max_spacing_, av_intens_, sd_intens_;	
		
		/** The hash map for the sweep line algorithm */
		SweepLineHash hash_;
		static hash_map<unsigned int, double> preComputedGamma_;

		std::list<double> mzsToGnuFiles_;
		unsigned int writtenGnuFiles_;
		
};  // end of class


template <typename MapType>
hash_map<unsigned int, double> IsotopeFinder<MapType>::preComputedGamma_;

inline bool comparator (const DRawDataPoint<2> a, const DRawDataPoint<2> b)
{
	return (a.getIntensity() > b.getIntensity());
}	


template <typename MapType>
IsotopeFinder<MapType>::IsotopeFinder () throw () : INTEGRATION_WORKSPACE(100), INTEGRATION_EPSILON (1e-6),
  PEAK_CUT_OFF (5), SCORE_CUT_OFF (0), RT_VOTES_CUT_OFF (6), WT_CUT_OFF (0), MZ_INTERLEAVE (2), waveletLength_ (0),
  avMZSpacing_ (0), avRTSpacing_ (0), min_spacing_ (0), max_spacing_ (1), av_intens_ (0), sd_intens_ (0), 
	writtenGnuFiles_ (0) 
{
}


template <typename MapType>
IsotopeFinder<MapType>::~IsotopeFinder () throw ()
{
}


template <typename MapType>
IsotopeFinder<MapType>::IsotopeFinder (const MapType& experiment, const unsigned int integration_workspace,
  double integration_epsilon, unsigned int peak_cut_off, double score_cut_off, unsigned int rt_votes_cut_off,
  double wtCutOff, double mz_interleave) throw () : INTEGRATION_WORKSPACE (integration_workspace), 
	INTEGRATION_EPSILON (integration_epsilon), PEAK_CUT_OFF (peak_cut_off), SCORE_CUT_OFF (score_cut_off), 
	RT_VOTES_CUT_OFF (rt_votes_cut_off), WT_CUT_OFF (wtCutOff), MZ_INTERLEAVE (mz_interleave),waveletLength_(0), 
	avMZSpacing_(0), min_spacing_ (0), max_spacing_ (1), av_intens_ (0), sd_intens_ (0), writtenGnuFiles_ (0)
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
	//experiment_.sortSpectra(true); Do we need this here?
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > signal; experiment_.get2DData(signal); 
	avMZSpacing_ = averageMZSpacing (signal);
	std::cout << "Average m/z spacing: " << avMZSpacing_ << std::endl;
	avRTSpacing_ = averageRTSpacing (signal);
	std::cout << "Average RT spacing: " << avRTSpacing_ << std::endl;
	waveletLength_ = (int) (PEAK_CUT_OFF/avMZSpacing_);	
	std::list<unsigned int> charges;
	charges.push_back(1); charges.push_back(2); 
	min_spacing_ = INT_MAX; max_spacing_= 0; av_intens_ =0; double tmp=0;
	for (unsigned int i=0; i<signal.size()-1; ++i)
	{
		tmp = signal[i+1].getPosition().Y() - signal[i].getPosition().Y();
		if (fabs(tmp) < min_spacing_)
			min_spacing_ = tmp;		
		if (fabs(tmp) > max_spacing_)
			max_spacing_ = tmp;
		av_intens_ += signal[i].getIntensity();
	}
	av_intens_ += signal[signal.size()-1].getIntensity();
	av_intens_ /= (double) signal.size();
	sd_intens_ = getAbsSd (signal, av_intens_, 0, signal.size()); 	
	std::cout << "Minimal m/z spacing: " << min_spacing_ << std::endl;
	std::cout << "Maximal m/z spacing: " << max_spacing_ << std::endl;
	std::cout << "Average intensity: " << av_intens_ << std::endl;
	std::cout << "Intensity sd: " << sd_intens_ << std::endl;
	generateGammaValues();
	prepareGNUplotFiles();
}


template <typename MapType>
void IsotopeFinder<MapType>::readTabFile (const std::string& filename) throw (Exception::FileNotFound)
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

	//Pushing the last scan
	spec.setRetentionTime(old_x);
	filter.filter (spec.begin(), spec.end(), spec2);
	spec2.setRetentionTime(old_x);
	experiment_.push_back(spec2);

	if (experiment_.empty())
	{
		std::cout << "Error: file is empty." << std::endl;
		exit(-1);
	}
	
	ifile.close();		
	initializeMe();	
}


template <typename MapType>
void IsotopeFinder<MapType>::cwtMulti (const unsigned int scanNumber, 
	const std::list<unsigned int>& charges, std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* pwts, 
	std::vector<double>* wt_thresholds) throw ()
{
	fastMultiCorrelate (experiment_[scanNumber].getContainer(), charges, pwts, wt_thresholds);
}


template <typename MapType>
typename IsotopeFinder<MapType>::SweepLineHash IsotopeFinder<MapType>::findFeatures 
	(const unsigned int start_scan, const unsigned int end_scan, const bool sweepLine) throw ()
{
	experiment_.updateRanges(); //this needs identifyCharge
	std::list<unsigned int> charges;
	charges.push_back(1); charges.push_back(2); 
	std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* pwts = NULL;
	std::vector<double>* wt_thresholds = NULL;
	for (unsigned int i=start_scan; i<=end_scan; ++i)
	{
		std::cout << "Spectrum " << i << " (" << experiment_[i].getRetentionTime() << ") of " << end_scan << std::endl; 
		pwts = new std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > > (charges.size(), experiment_[i].getContainer());
		wt_thresholds = new std::vector<double> (charges.size(), 0);
		cwtMulti (i, charges, pwts, wt_thresholds);
		identifyCharge (*pwts, wt_thresholds, i, experiment_[i].getRetentionTime());
		delete (pwts);
		delete (wt_thresholds);
	}

	if (sweepLine)
    filterHashByRTVotes ();

	return (hash_);
}
				
		
template <typename MapType>
unsigned int IsotopeFinder<MapType>::unique_merge (const std::list<double>& a, const std::list<double>& b, 
	std::vector<double>::iterator& c_begin, std::vector<double>::iterator& c_end, std::vector<double>* res) throw ()
{
	std::list<double>::const_iterator iter_a, iter_b;
  double old;
	iter_a=a.begin(); iter_b=b.begin();
	unsigned int occ=0;
	while (iter_a != a.end() && iter_b != b.end())
	{
			while (*iter_a < *iter_b)
			{
				(*res)[occ++] = *iter_a;
				if (++iter_a == a.end())
					break;
			}
			while (*iter_b < *iter_a)
			{
				(*res)[occ++] = *iter_b;
				if (++iter_b == b.end())
					break;
			}
			
			if (*iter_a == *iter_b)
				(*res)[occ++] = *iter_a;

			old = (*res)[occ-1];
			
			while (*iter_a == old)
				if (++iter_a == a.end())
					break;
		
			
			while (*iter_b == old)				
				if (++iter_b == b.end())
					break;

	}

	while (iter_a != a.end())
	{
		(*res)[occ++] = *iter_a;
		++iter_a;
	}
	while (iter_b != b.end())
	{
		(*res)[occ++] = *iter_b;
		++iter_b;
	}

	c_begin = res->begin();
	c_end = res->begin()+occ;
	return (occ);
}


template <typename MapType>
void IsotopeFinder<MapType>::printMapEntry (SweepLineMap::iterator iter) throw ()
{
	std::cout << getMZbyHashKey(iter->first) << "\t [ ";
	for (std::list<double>::const_iterator iter_cl2=iter->second.first.begin(); 
		iter_cl2 != iter->second.first.end(); ++iter_cl2)
	std::cout	<< *iter_cl2 << " "; std::cout << "]  { ";	
	for (std::list<double>::const_iterator iter_l=iter->second.second.begin();
		iter_l != iter->second.second.end(); ++iter_l)
	std::cout << *iter_l << " ";
	std::cout << "}" << std::endl;	
}
			

template <typename MapType>
void IsotopeFinder<MapType>::printHashEntry (SweepLineHash::iterator iter) throw ()
{
	std::cout << getMZbyHashKey(iter->first) << "\t [ ";
	for (std::list<double>::const_iterator iter_cl2=iter->second.first.begin(); 
		iter_cl2 != iter->second.first.end(); ++iter_cl2)
	std::cout	<< *iter_cl2 << " "; std::cout << "]  { ";	
	for (std::list<double>::const_iterator iter_l=iter->second.second.begin();
		iter_l != iter->second.second.end(); ++iter_l)
	std::cout << *iter_l << " ";
	std::cout << "}" << std::endl;	
}



template <typename MapType>
void IsotopeFinder<MapType>::filterHashByRTVotes () throw ()
{
	std::cout << "Hash size before filtering: " << hash_.size() << std::endl << std::endl;	
	
	SweepLineHash::iterator iter_f, iter_b;
 	std::vector<SweepLineHash::iterator> toDelete;
	std::list<double>::iterator iter_l;

	for (SweepLineHash::iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
	{
		if (iter->second.first.size() <= RT_VOTES_CUT_OFF)
			toDelete.push_back(iter);
	}

	for (unsigned int i=0; i<toDelete.size(); ++i)
		hash_.erase (toDelete[i]);	

	std::cout << "Hash size after filtering: " << hash_.size() << std::endl << std::endl;	
}


template <typename MapType>
/*std::pair<typename IsotopeFinder<MapType>::WaveletCollection, typename IsotopeFinder<MapType>::WaveletCollection>*/ 
void IsotopeFinder<MapType>::identifyCharge (const std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >& candidates, std::vector<double>* wt_thresholds, const unsigned int scan, const double RT) throw ()
{	
	std::vector<double> int_mins (candidates[0].size(),INT_MIN), zeros (candidates[0].size(),0);
	WaveletCollection scoresC (candidates.size(), zeros);

	std::vector<unsigned int> start_indices, end_indices;
	//In order to determine the start and end indices, we first need to know the width of the region one should consider 
	//to estimate the mean and the sd of the pattern candidate. 
	//That region is defined by the position of the heighst amplitude +/- waveletLength_.
	
	typedef MSSpectrum<DRawDataPoint<2> >::ContainerType containerType;
	containerType::iterator iter; 
	unsigned int start_index, end_index, c_index, i_iter; //Helping variables
	double seed_mz, c_check_point, c_val, c_av_intens;
	std::vector<bool> processed (candidates[0].size(), false);
	std::pair<int, int> c_between; 			
	int start, end, goto_left;

	for (unsigned int c=0; c<candidates.size(); ++c)		
	{
		processed = std::vector<bool> (candidates[0].size(), false); //Reset
		containerType c_candidate = candidates[c]; 
		
		//Ugly, but do not how to do this in a better (and easy) way
		for (unsigned int i=0; i<c_candidate.size(); ++i)
			c_candidate[i].setPosition(DPosition<2>(c_candidate[i].getPosition().X(), i));

		sort (c_candidate.begin(), c_candidate.end(), comparator); 
		c_av_intens = getAbsMean (candidates[c], 0, candidates[c].size());

		for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter)
		{
			if (iter->getIntensity() <= (*wt_thresholds)[c]*5*c_av_intens)
				break;
		}

		c_candidate.erase (iter, c_candidate.end());
		

		i_iter=0;
		for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter, ++i_iter)
		{
			c_index = (int) (iter->getPosition().Y());		
			
			if (processed[c_index])
        continue;		
			
			start_index = c_index-waveletLength_-1;
			end_index = c_index+waveletLength_+1;				
			seed_mz=iter->getPosition().X();

			//Catch impossible cases
			// Unsigned int < zero is always false !
      if (/*start_index < 0 ||*/ end_index >= candidates[c].size() || start_index > end_index)
        continue;			
	
			//Mark as processed
			for (unsigned int z=start_index; z<=end_index; ++z)
				processed[z] = true;	
			
			start=(-2*(PEAK_CUT_OFF-1))+1, end=(2*(PEAK_CUT_OFF-1))-1;
			goto_left = c_index - waveletLength_ - 1; 
			//std::cout << start << "\t" << end << std::endl;
			for (int v=start; v<=end; ++v)
			{
				c_check_point = seed_mz+v*0.5/((double)c+1);
				c_between = getNearBys (scan, c_check_point, goto_left);
				if (c_between.first < 0 || c_between.second < 0)
					break;
				c_val = getInterpolatedValue (candidates[c][c_between.first].getPos(), 
					c_check_point, 
					candidates[c][c_between.second].getPos(), 
					candidates[c][c_between.first].getIntensity(), 
					candidates[c][c_between.second].getIntensity());

 				if (fabs(c_val) < c_av_intens)
					continue;
				
				if (abs(v)%2 == 1) //i.e. whole
					scoresC[c][c_index] -= c_val;
				else //i.e. peak
					scoresC[c][c_index] += c_val; 
			}

			if (scoresC[c][c_index] <= 1.5*iter->getIntensity())
				scoresC[c][c_index] = 0;
		}	
	}

	//Now, since we computed all scores, we can hash all mz positions
 	
	unsigned int numOfCharges = candidates.size(), numOfMZPositions = candidates[0].size();
	//This is a vector telling us the next mz position in charge i we have to hash
	std::vector<unsigned int> positions (numOfCharges, 0); 
	
	unsigned int c_hash_key, count_finished_charges=0; double allZero; 
	DoubleList c_pair; 	std::list<double> c_list, c_fill_list; std::list<double>::iterator iter_cl, iter_cl_hash;
	std::pair<SweepLineHash::iterator, SweepLineHash::iterator> iter_hash;
	while (1) // ;-)
	{	
		//Termination criterion
		//Test for every charge ...
		for (unsigned int c=0; c<numOfCharges; ++c)
		{
			//... if we hashed already all possible mz coordinates
			if (positions[c] >= numOfMZPositions)
			{
				//if so ... goto FINISHED_HASHING
				if (++count_finished_charges >= numOfCharges)
					goto FINISHED_HASHING;
				positions[c] = -1;
			}
		}
		//End of Termination criterion
						
		
		//The hashing
		for (unsigned int c=0; c<numOfCharges; ++c)		
		{
			//if positions[c]<0 then all mz coordinates of this charge are hashed
			//this is also the case if the number of entries equals numOfMZPositions
			
			// unsigned int < 0 is always false !
			if (/*positions[c] < 0 ||*/ positions[c] >= numOfMZPositions)
				continue;
			
			//positions[c] also tells us the next candidate to hash
			c_list.push_back (scoresC[c][positions[c]++]);
		}				
	
		for (unsigned int c=0; c<numOfCharges-1; ++c)		
			if (positions[c+1] != positions[c])
			{
				std::cout << "Quadro Zack!" << std::endl;
				exit(-1);
			}


		//std::cout << experiment_[scan].getContainer()[positions[0]-1].getPos() << std::endl;	
		//std::cout << experiment_.getMin().Y() << std::endl;
		c_hash_key = (unsigned int) ((experiment_[scan].getContainer()[positions[0]-1].getPos() - experiment_.getMin().Y())
			/ avMZSpacing_);
		
		allZero=true;	
		for (iter_cl=c_list.begin(); iter_cl!=c_list.end(); ++iter_cl)
			if (*iter_cl != 0)
				allZero=false;
		
		if (!c_list.empty() && !allZero)
		{
			iter_hash=hash_.equal_range(c_hash_key);
			while (iter_hash.first != iter_hash.second)		
			{
				//std::cout << "M/Z already in hash table ... " << std::endl;

				if (scan != 0)	
					if (find(iter_hash.first->second.first.begin(), iter_hash.first->second.first.end(), 
						experiment_[scan-1].getRetentionTime()) == iter_hash.first->second.first.end())
					{
						//i.e. there is no neighbouring entry before this retention time
						//i.e. we can treat this case as if no entry is present in the hash
						++iter_hash.first; 
						continue;
					}
				c_fill_list = iter_hash.first->second.first; c_fill_list.push_back(RT);	
				c_fill_list.unique(); //It might be the case, that we have several votes for the same RT and MZ by different
				//charges => unique the list.
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
			c_fill_list.push_back(RT);	
			c_pair = DoubleList (c_fill_list, c_list);

		FINISH:			
			hash_.insert (SweepLineHash::value_type(c_hash_key, c_pair));
		}

		c_list.clear();
	}

 FINISHED_HASHING:
	
	return;	
}	


template <typename MapType>
double IsotopeFinder<MapType>::getMean (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal,
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += signal[i].getIntensity();

	return (res/(double)(endIndex-startIndex+1));	
}

				
template <typename MapType>
double IsotopeFinder<MapType>::getAbsMean (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal,
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
			res += fabs(signal[i].getIntensity());

	return (res/(double)(endIndex-startIndex+1));	
}


template <typename MapType>
double IsotopeFinder<MapType>::getUpShiftedMoment (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal,
	const unsigned int moment, const unsigned int startIndex, const unsigned int endIndex) throw ()
{
	double tmp=0, min=INT_MAX; 
	for (unsigned int i=startIndex; i<endIndex; ++i)
		if (signal[i].getIntensity() < min)
			min = signal[i].getIntensity();
	
	if (min < 0)
		min *=-1;
	else min=0;
	
	double mean=0, res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += signal[i].getIntensity() + min;
	mean = (res/(double)(endIndex-startIndex+1));

	
	if (moment == 1)
		return (mean);

	if (moment == 2) //I.e. variance
	{
		res=0;
		for (unsigned int i=startIndex; i<endIndex; ++i)
		{
			tmp = signal[i].getIntensity() + min - mean;
			res += tmp*tmp;
		}
		return (sqrt(res/(double)(endIndex-startIndex)));
	}
	
	//I.e. skewness
	
	res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
	{
		tmp = signal[i].getIntensity() + min - mean;
		res += tmp*tmp*tmp;
	}

	return ((pow(fabs(res/(double)(endIndex-startIndex)), 1.0/3.0)));
}


template <typename MapType>
double IsotopeFinder<MapType>::getSd (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const double mean, 
	const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
	double res=0;
	for (unsigned int i=startIndex; i<endIndex; ++i)
		res += (signal[i].getIntensity()-mean)*(signal[i].getIntensity()-mean);

	return (sqrt(res/(double)(endIndex-startIndex)));	
}


template <typename MapType>
void IsotopeFinder<MapType>::fastMultiCorrelate 
	(const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const std::list<unsigned int>& charges, 
	 std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* pwts, std::vector<double>* wt_thresholds) throw ()
{					
	std::vector<DPeakArrayNonPolymorphic<1, DRawDataPoint<2> > >* res = pwts;
	unsigned int signal_size = signal.size();
	
	WaveletCollection phis (charges.size(), std::vector<double> (waveletLength_)); //all necessary wavelets (by rows)
	
	double cumSpacing=0, cSpacing=0, realMass=0, lambda=0, w_sum=0, w_s_sum=0, max_w_monoi_intens=0.25, 
		align_offset, tmp_pos, tmp_pos1; //helping variables
	std::list<unsigned int>::const_iterator charge_iter; unsigned int k=0; //helping variables
	std::list<double>::iterator formzs;
	std::list<unsigned int> formzs_indices;
	std::list<WaveletCollection> back_phis;
	std::list<double> tmpMzsToGnuFiles = mzsToGnuFiles_;
					
	//double integration=0;
	double max=0;
	formzs=tmpMzsToGnuFiles.begin();
	for (unsigned int i=0; i<signal_size; ++i)
	{	
		//Now, let's sample the wavelets
		for (charge_iter=charges.begin(), k=0; charge_iter!=charges.end(); ++charge_iter, ++k)
		{
			cumSpacing=0; w_sum=0; w_s_sum=0;				
			realMass = signal[i].getPos() * (*charge_iter);
			lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
			//integration = phiRawInt (lambda, 1.0/(double)(*charge_iter));
			//std::cout << "Exact: " << integration << "\t" << lambda << "\t" << 1.0/(double)(*charge_iter) << std::endl;
			
			max_w_monoi_intens=0.25/(*charge_iter);

			//Align the maximum monoisotopic peak of the wavelet with some signal point
			unsigned int j=0; double last=0;
			while (cumSpacing < max_w_monoi_intens)
			{
				cSpacing = signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos();
			 	last=cumSpacing;	
				if (cSpacing < 0)
					cumSpacing += avMZSpacing_;
				else //The "normal" case
					cumSpacing += signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos(); 
				++j;
			}
			
			align_offset = max_w_monoi_intens-last;
			
			cumSpacing=align_offset;	
			for (unsigned int j=0; j<waveletLength_; ++j)
			{
				tmp_pos = signal[(i+j)%signal_size].getPos();	
				tmp_pos1 = signal[(i+j+1)%signal_size].getPos();
				
				realMass = tmp_pos1 * (*charge_iter);
				lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
				//std::cout << "phi before: " << phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter))
					//<< "\tphi after" << phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter))- integration << std::endl;
				//std::cout << "cumSpacing: " << cumSpacing << std::endl;
				phis[k][j] = phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter));
				w_sum += phis[k][j];
				w_s_sum += phis[k][j] * phis[k][j];
				//std::cout << "t: " << cumSpacing << "\t" << lambda << "\t" << 1.0/(double)(*charge_iter)
					//			<< "\t" << integration << "\t" << phis[k][j] << std::endl;
				cSpacing = tmp_pos1 - tmp_pos;
				//cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal_size-1).
				//Since this case is only of theoretical interest (we do not expect any important signals at the very end of the 
				//spectrum), we simlply use the average spacing in that case.
				if (cSpacing < 0)
					cumSpacing += avMZSpacing_;
				else //The "normal" case
					cumSpacing += tmp_pos1 - tmp_pos; 
			}
			max=-INT_MAX;
			for (unsigned int j=0; j<waveletLength_; ++j)
			{
				phis[k][j] -= (w_sum/(double)waveletLength_);
				if (phis[k][j] > max)
					max = phis[k][j];
			}
			for (unsigned int j=0; j<waveletLength_; ++j)
				phis[k][j] /= max;

			(*wt_thresholds)[k] = w_s_sum; 							
		}			
					
		std::vector<double> sums (charges.size());
		k=0;
		for (unsigned int j=i; j<signal_size && k<phis[0].size(); ++j, ++k)
		{	//Since all wavelet functions have the same length, we can simply use phis[0].size()
			
			for (unsigned int m=0; m<charges.size(); ++m)
			{
				sums[m] += signal[j].getIntensity()*phis[m][k];
			}
		}

		for (unsigned int l=0; l<i && k<phis[0].size(); ++l, ++k)
		{	//Since all wavelet functions have the same length, we can simply use phis[0].size()
	
			for (unsigned int m=0; m<charges.size(); ++m)
				sums[m] += signal[l].getIntensity()*phis[m][k];
		}

		//Store the current convolution result
		for (unsigned int m=0; m<charges.size(); ++m)
			(*res)[m][i].setIntensity(sums[m]);

		if (formzs == tmpMzsToGnuFiles.end())
			continue;
		
		//std::cout << "*formzs: " << *formzs <<	"\t signal[i].getPos() " << signal[i].getPos() << std::endl;
		if (*formzs == signal[i].getPos())
		{
			std::cout << "index+1: " << i+1 << " " << signal[i+1].getPos() << std::endl;
			formzs_indices.push_back (i);
			back_phis.push_back (phis);
		}
		while (*formzs <= signal[i].getPos())
		{
			tmpMzsToGnuFiles.erase (formzs);
			formzs = tmpMzsToGnuFiles.begin();
			if (formzs == tmpMzsToGnuFiles.end())
				break;
		}
		
	}

	std::list<unsigned int>::iterator iterf;	
	std::list<WaveletCollection>::iterator iterbp;

	for (unsigned int m=0; m<charges.size(); ++m)
	{
		for (formzs=mzsToGnuFiles_.begin(), iterf=formzs_indices.begin(), iterbp=back_phis.begin(); 
			formzs != mzsToGnuFiles_.end(); ++formzs, ++iterf, ++iterbp)		
			createGNUplot (*iterf, *formzs, m+1, &signal, &((*iterbp)[m]), &(*res)[m]);
	}
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
double IsotopeFinder<MapType>::averageRTSpacing (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal) 
	const throw ()
{
	double RTspacing=0, avRTspacing=0, counts=0;
	for (unsigned int i=0; i<signal.size()-1; ++i)
	{			
		if (signal[i+1].getPosition().X() != signal[i].getPosition().X()) //i.e. a new scan begins
		{
			RTspacing = signal[i+1].getPosition().X() - signal[i].getPosition().X();
			avRTspacing += RTspacing;
			++counts;
		}
	}

	if (counts==0)
		return (1); //should be neutral
	
	avRTspacing /= counts;		

	return (avRTspacing);		
}

				
template <typename MapType>
double IsotopeFinder<MapType>::getAbsSd (const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >& signal, const double mean,
  const unsigned int startIndex, const unsigned int endIndex) const throw ()
{
  double res=0, tmp;
  for (unsigned int i=startIndex; i<endIndex; ++i)
	{
		tmp = fabs(signal[i].getIntensity())-mean;
    res += tmp*tmp;
	}

  return (sqrt(res/(double)(endIndex-startIndex)));
}
				
				
template <typename MapType>
double IsotopeFinder<MapType>::gslPhiRaw (double t, void* params) throw ()
{
	double* parameters = reinterpret_cast<double*> (params);
	double lambda = parameters[0];
	double a = parameters[1];
	double min_spacing = parameters[2];
	
	int x0, x1; double f0, f1, fi;
	double res=0;
	x0 = (int) trunc ((t/a + 1)/min_spacing);
	x1 = x0+1;
	if ((unsigned int) x1 < preComputedGamma_.size())
	{
		f0 = preComputedGamma_[x0];
		f1 = preComputedGamma_[x1];
		fi = (f0 + (f1-f0)/((x1-x0)*min_spacing) * ((t/a+1)-x0*min_spacing));
		res = (sin(2*M_PI*t/a) * exp(-lambda)) * ((pow(lambda,t/a)) / fi);
		return (res);
	}
			
	res = (sin(2*M_PI*t/a) * exp(-lambda) * (pow(lambda,t/a)) / tgamma ((t/a)+1)); 
	return (res);
}


template <typename MapType>
void IsotopeFinder<MapType>::generateGammaValues () throw ()
{
	std::cout << "Precomputing the Gamma function ...";
	preComputedGamma_.clear();
	double query = 0; unsigned int counter=0;
	while (query <= 4*PEAK_CUT_OFF +1) //4 because of max_charge	
	{
		preComputedGamma_[counter] = tgamma (query);
		query += min_spacing_;	
		++counter;
	}
	std::cout << " done." << std::endl;
}



template <typename MapType>
double IsotopeFinder<MapType>::phiRawInt (const double lambda, const double a) throw ()
{
	gsl_function F;
	gsl_integration_workspace* GSL_WSP = gsl_integration_workspace_alloc (INTEGRATION_WORKSPACE);
	F.function = &IsotopeFinder<MapType>::gslPhiRaw;
	double* parameters = new double [3]; double iValue, absError;
	parameters[0] = lambda;
	parameters[1] = a;
	parameters[2] = min_spacing_;
	F.params = parameters;
	//unsigned int neval;
	
	std::cout << "Warning: you are using slow qag integration. " << std::endl;
	
	gsl_integration_qag (&F, 0, PEAK_CUT_OFF, INTEGRATION_EPSILON, 0.01, 
		INTEGRATION_WORKSPACE-1, GSL_INTEG_GAUSS61, GSL_WSP, &iValue, &absError);  
	gsl_integration_workspace_free (GSL_WSP);

	return (iValue);
}


template <typename MapType>
void IsotopeFinder<MapType>::prepareGNUplotFiles (const std::string& file) throw ()
{
	mzsToGnuFiles_.clear();			
	std::ifstream ifile (file.c_str());
	double cMZ;
	while (ifile >> cMZ)
		mzsToGnuFiles_.push_back(cMZ);	
	
	ifile.close();
}


template <typename MapType>
void IsotopeFinder<MapType>::createGNUplot (unsigned int alignedTo, double mz, unsigned int charge,  
	const DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >* signal, std::vector<double>* wavelet, 
	DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >* transform) throw ()
{
	std::ofstream ofile;
	
	if (signal != NULL)
	{
		std::stringstream str;
		str << "s" << mz << ".wt" << '\0';
		ofile.open (str.str().c_str());
		DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >::const_iterator iter;
		for (iter = signal->begin(); iter != signal->end(); ++iter)
				ofile << iter->getPos() << "\t" << iter->getIntensity() << std::endl;
		ofile.close();
	}
	
	if (wavelet != NULL)	
	{
		unsigned int i=0;
		std::stringstream str;

		str << "w" << mz << "_" << charge << ".wt" << '\0';
		ofile.open (str.str().c_str());
		std::vector<double>::iterator iter;
		std::cout << "alignedTo " << alignedTo << "\t" << "mz " << mz << std::endl;
		for (iter = wavelet->begin(), i=alignedTo; iter != wavelet->end(); ++iter, ++i)
			ofile << (*signal)[i].getPos() << "\t" << *iter << std::endl;
		ofile.close();
		str.clear();
	}	
	
	if (transform != NULL)				
	{	
		std::stringstream str;
		str << "t" << mz << "_" << charge << ".wt" << '\0';
		ofile.open (str.str().c_str());
		DPeakArrayNonPolymorphic<1, DRawDataPoint<2> >::iterator iter;
		for (iter = transform->begin(); iter != transform->end(); ++iter)
			ofile << iter->getPos() << "\t" << iter->getIntensity() << std::endl;
		ofile.close();
		str.clear();
	}
	
	++writtenGnuFiles_;
}



} //namespace

#endif
