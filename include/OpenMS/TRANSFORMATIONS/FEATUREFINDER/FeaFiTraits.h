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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H

#include <OpenMS/DATASTRUCTURES/ScanIndexMSExperiment.h>
#include <OpenMS/DATASTRUCTURES/IndexSet.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>

#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

namespace OpenMS
{
/**
	 @brief Traits class for the feature finding algorithm.
	 
	 This class is rather an "umbrella" for the different modules / steps of the algorithm
	 than a traits class in the traditional sense.
	
	 @ingroup FeatureFinder 	
**/
class FeaFiTraits
{

public:

    /// Defines the coordinates of peaks / features.
    enum DimensionId
    {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    /// Flag for each data point
    enum Flag { UNUSED, SEED, INSIDE_FEATURE };

    typedef MSExperimentExtern<DPeak<1> > MapType;

    typedef std::vector<Flag> FlagVector;

    typedef DRawDataPoint<2> PeakType;
    typedef DRawDataPoint<2>::IntensityType IntensityType;
    typedef DRawDataPoint<2>::CoordinateType CoordinateType;
    typedef DRawDataPoint<2>::PositionType PositionType;
    typedef DFeature<1>::ChargeType ChargeType;

    typedef PeakType::NthPositionLess< RT > RTless;
    typedef PeakType::NthPositionLess< MZ > MZless;

    typedef DFeatureMap<2> FeatureVector;
    typedef DFeature<2>::ConvexHullType ConvexHullType;
    typedef FeaFiModule::NoSuccessor NoSuccessor;

    /// standard constructor
    FeaFiTraits() {}

    /// destructor
    virtual ~FeaFiTraits() {}

    /// copy constructor
    FeaFiTraits(const FeaFiTraits& source)
            : map_(source.map_),
            flags_(source.flags_),
            scan_index_(source.scan_index_),
            features_(source.features_)
    {}

    /// assignment operator
    FeaFiTraits& operator = (const FeaFiTraits& source)
    {
        if (&source == this)
            return *this;

        map_             = source.map_;
        flags_             = source.flags_;
        scan_index_   = source.scan_index_;
        features_        = source.features_;

        return *this;
    }

    /// set internal data and update range information
    void setData(MapType& exp)
    {			
				if (exp.size() == 0)
				{
					std::cout << "No data provided. Aborting. " << std::endl;
					return;
				}
		
				std::cout << "Storing MSExperimentExtern " << std::endl;
				map_.setBufferSize( exp.getBufferSize() );
				map_.updateBuffer();
	
				// copy scanwise such that we can remove tandem spectra
				for (UnsignedInt i=0; i<exp.size(); ++i)
				{
					if (exp[i].getMSLevel() == 1) map_.push_back(exp[i]);
				}	
		 
				std::cout << "Updating range information. " << std::endl;
        // update range informations
        map_.updateRanges();				
				
				std::cout << "This map contains " << map_.size() << " scans  ";
				std::cout << "and " << map_.getSize() << " data points. " << std::endl;

				std::cout << "Setting flags. " << std::endl;
        // resize internal data structures
        flags_.reserve(map_.getSize());

				// set peak flags
        for (UnsignedInt i=0; i<map_.getSize(); ++i)
            flags_.push_back(FeaFiTraits::UNUSED);
				
				std::cout << "Initialising scan index DS. " << std::endl;					
				scan_index_.init ( map_.peakBegin(), map_.peakEnd() );
   }
		
		/// copy input data to external memory and update range information 
		/// NOTE: This is slow since all peaks are copied individually
    void setData(MSExperiment<DPeak<1> >& exp)
    {
				if (exp.size() == 0)
				{
					std::cout << "No data provided. Aborting. " << std::endl;
					return;
				}
		
				std::cout << "Storing MSExperiment " << std::endl;
			
				for (UnsignedInt i=0; i<exp.size(); ++i)
				{
					if (exp[i].getMSLevel() == 1) map_.push_back(exp[i]);
				}	
								
				std::cout << "Updating range information. " << std::endl;
        // update range informations
        map_.updateRanges();
				
				std::cout << "This map contains " << map_.size() << " scans ";
				std::cout << "and " << map_.getSize() << " data points. " << std::endl;

				std::cout << "Setting flags. " << std::endl;
        // resize internal data structures
        flags_.reserve(map_.getSize());

				// set peak flags
        for (UnsignedInt i=0; i<map_.getSize(); ++i)
            flags_.push_back(FeaFiTraits::UNUSED);
						
				if (map_.getSize() == 0)
				{
					std::cout << "No data provided. Aborting. " << std::endl;
					return;
				}
				std::cout << "Initialising scan index DS. " << std::endl;
        scan_index_.init ( map_.peakBegin(), map_.peakEnd() );
    }
			
		/// Mutable access to LC-MS map
		MapType& getData() { return map_; }
		/// Const access to LC-MS map
		const MapType& getData() const { return map_; }
		
    /// non-mutable access flag with index @p index .
    const Flag& getPeakFlag(const UnsignedInt index) const throw (Exception::IndexOverflow) {  return flags_.at(index); }
    /// mutable access flag with index @p index.
    Flag& getPeakFlag(const UnsignedInt index) throw (Exception::IndexOverflow) { return flags_.at(index); }
    /// access peak with index @p index.
   	PeakType getPeak(const UnsignedInt index) const throw (Exception::IndexOverflow) { return map_.getPeak(index);  }

    /// retrieve the number of peaks.
    const UnsignedInt getNumberOfPeaks()  { return map_.getSize(); }
		/// Retrieve index datastructure 
    const ScanIndexMSExperiment<MapType >& getScanIndex() { return scan_index_; }
	
    /// access intensity of peak with index @p index.
    const IntensityType& getPeakIntensity(const UnsignedInt index) const throw (Exception::IndexOverflow) { return map_.getPeak(index).getIntensity(); }
    /// access m/z of peak with index @p index .
    const CoordinateType& getPeakMz(const UnsignedInt index) const throw (Exception::IndexOverflow) { return map_.getPeak(index).getPosition()[MZ]; }
    /// access retention time of peak with index @p index.
    const CoordinateType& getPeakRt(const UnsignedInt index) const throw (Exception::IndexOverflow) { return map_.getPeak(index).getPosition()[RT]; }
		
    /// access scan number of peak with index @p index
    const UnsignedInt getPeakScanNr(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
        if (index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "SimpleFeaFiTraits::getScanNr()", index, map_.getSize());
        CoordinateType current_rt = getPeakRt(index);
        return scan_index_.getRank(current_rt);
    }

    /** @brief get index of next peak in m/z dimensio.

    \param index of the peak whose successor is requested
    \return index of the next peak 
    */
    UnsignedInt getNextMz(UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor)
    {
        if (index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index, map_.getSize());
        if (index >= (map_.getSize()-1) )
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);

        // check whether we walked out of the current scan i.e. the retention
        // time has changed
        if (getPeakRt(index) != getPeakRt(index+1))
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);

        // since we sorted by rt and then by m/z, the peak with the same rt but
        // the larger m/z is simply one step further in the peak vector
        return ++index;
    }

    /** @brief get index of previous peak in m/z dimension.

    \param index of the peak whose predecessor is requested
    \return index of the previous peak
    */
    UnsignedInt getPrevMz(UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor)
    {
        if (index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, map_.getSize());

        // if we are at the beginning of the peak vector, there will be no previous peak ;-)
        if (index == 0)
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);

        // check whether we walked out of the current scan i.e. the retention
        // time has changed (same problem as above in nextMz() )
        if (getPeakRt(index) != getPeakRt(index-1))
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);

        // same as above
        return --index;
    }

    /** @brief get index of next peak in retention time dimension.
     
    \param index of the peak whose successor is requested
    \return index of the next peak
    */
    UnsignedInt getNextRt(const UnsignedInt index) throw (Exception::IndexOverflow, NoSuccessor)
    {
        if (index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, map_.size());

        const PeakType p  = map_.getPeak(index);

        MapType::PIterator piter;
        try
        {
            piter = scan_index_.getNextRt(p.getPosition()[RT],p.getPosition()[MZ]);
        }
        catch (Exception::Base ex)
        {
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
        }
        
				UnsignedInt peak_index = piter.getPeakNumber();
				
        if (peak_index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, map_.size());

        return peak_index;
    }

    /** @brief get index of next peak in retiontion time dimension.

    \param index of the peak whose predecessor is requested
    \return index of the previous peak
    */
    UnsignedInt getPrevRt(const UnsignedInt index) throw (Exception::IndexOverflow, NoSuccessor)
    {
        if (index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevRt", index, map_.size());

        const PeakType p = getPeak(index);
        MapType::PIterator piter;
        try
        {
            piter = scan_index_.getPrevRt(p.getPosition()[RT],p.getPosition()[MZ]);
        }
        catch (Exception::Base ex)
        {
            throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevRt", index);
        }
  
        UnsignedInt peak_index = piter.getPeakNumber();
				
        if (peak_index>=map_.getSize())
            throw Exception::IndexOverflow(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index, map_.size());

        return peak_index;
    }

    /// run main loop
    const FeatureVector& run(const std::vector<BaseSeeder*>& seeders,
                             const std::vector<BaseExtender*>& extenders,
                             const std::vector<BaseModelFitter*>& fitters)
    {
        // Visualize seeds and extension in TOPPView:
        // get all Seeds and save corresponding peaks as "features"
        // get convex hull of the extension and use it for the "feature"

        // counts the number of features collected so far,
        // is needed for the gnuplot output.
#ifdef DEBUG_FEATUREFINDER
        int nr_feat = 0;
#endif

        if (map_.getSize() == 0)
        {
            std::cout << " No data provided! Aborting..." << std::endl;
            return features_;
        }

        if (seeders.size() == 0 ||
             extenders.size() == 0 ||
             fitters.size() == 0)
        {
            std::cout << " No modules set. Aborting..." << std::endl;
            return features_;
        }

        // gather information for fitting summary
        std::map<String,int> exception;									//count exceptions
        int no_exceptions = 0;
        std::map<String,int> mz_model;									//count used mz models
        std::map<float,int> mz_stdev;										//count used mz standard deviations
        std::vector<int> charge(10);											//count used charges
        double corr_mean=0.0, corr_max=0.0, corr_min=1.0; 	//boxplot for correlation

        StopWatch watch;
        unsigned int seed_count = 0;
        try
        {
            while (true)
            {
//                	UnsignedInt seed = seeders[0]->nextSeed();
								IndexSet seed_region = seeders[0]->nextSeed();

                watch.start();
                std::cout << "Extension ..." << std::endl;
                IndexSet peaks = extenders[0]->extend(seed_region);
                watch.stop();
                std::cout << "Time spent for extension: " << watch.getClockTime() << std::endl;
                watch.reset();
                ++seed_count;
								std::cout << "This is seed nr " << seed_count << std::endl;
                try
                {

                    watch.start();
                    features_.push_back(fitters[0]->fit(peaks));
                    watch.stop();
                    std::cout << "Time spent for fitting: " << watch.getClockTime() << std::endl;
                    watch.reset();

										#ifdef DEBUG_FEATUREFINDER
                    writeGnuPlotFile_(peaks,false,nr_feat++);
										#endif
                    // gather information for fitting summary
                    const DFeature<2>& f = features_[features_.size()-1];

                    float corr = f.getOverallQuality();
                    corr_mean += corr;
                    if (corr<corr_min)
                        corr_min = corr;
                    if (corr>corr_max)
                        corr_max = corr;

                    // count estimated charge states
                    unsigned int ch = f.getCharge();
                    if (ch>= charge.size())
                    {
                        charge.resize(ch);
                    }
                    charge.at(ch)++;

                    const Param& p = f.getModelDescription().getParam();
                    ++mz_model[ p.getValue("MZ") ];

                    DataValue dp = p.getValue("MZ:isotope:stdev");
                    if (dp != DataValue::EMPTY)
                    {
                        ++mz_stdev[p.getValue("MZ:isotope:stdev")];
                    }

                }
                catch( BaseModelFitter::UnableToFit ex)
                {
                    // set unused flag for all data points
                    for (IndexSet::ConstIterator it=peaks.begin(); it!=peaks.end(); ++it)
                    {
                        getPeakFlag(*it) = FeaFiTraits::UNUSED;
                    }
                    std::cout << " " << ex.what() << std::endl;
                    watch.stop();
                    std::cout << "Time spent for fitting: " << watch.getClockTime() << std::endl;
                    watch.reset();
                    ++no_exceptions;
                    ++exception[ex.getName()];
										
										#ifdef DEBUG_FEATUREFINDER
                    writeGnuPlotFile_(peaks,false,nr_feat++);
										#endif
                }

            } // end of while(true)
        }
        catch(NoSuccessor ex)
        { }

        // Print summary:
        Size size = features_.size();

        std::cout << size << " features were found. " << std::endl;

//         std::cout << "seed count " << seed_count << std::endl;

        std::cout << "FeatureFinder summary:\n"
        << "Correlation:\n\tminimum: " << corr_min << "\n\tmean: " << corr_mean/size
        << "\n\tmaximum: " << corr_max << std::endl;

        std::cout << "Exceptions:\n";
        for (std::map<String,int>::const_iterator it=exception.begin(); it!=exception.end(); ++it)
        {
            std::cout << "\t" << it->first << ": " << it->second*100/no_exceptions
            << "% (" << it->second << ")\n";
        }

        std::cout << "Chosen mz models:\n";
        for (std::map<String,int>::const_iterator it=mz_model.begin(); it!=mz_model.end(); ++it)
        {
            std::cout << "\t" << it->first << ": " << it->second*100/size
            << "% (" << it->second << ")\n";
        }

        std::cout << "Chosen mz stdevs:\n";
        for (std::map<float,int>::const_iterator it=mz_stdev.begin(); it!=mz_stdev.end(); ++it)
        {
            std::cout << "\t" << it->first << ": " << it->second*100/(size-charge[0])
            << "% (" << it->second << ")\n";
        }

        std::cout << "Charges:\n";
        for (unsigned int i=1; i<charge.size(); ++i)
            if (charge[i]!=0)
            {
                std::cout << "\t+" << i << ": " << charge[i]*100/(size-charge[0])
                << "% (" << charge[i] << ")\n";
            }

#ifdef DEBUG_FEATUREFINDER
        IndexSet empty;
        writeGnuPlotFile_(empty,true,nr_feat);
#endif

        return features_;

    } // end of run(seeders, extenders, fitters)

    ///@brief Calculate the convex hull of the peaks contained in @p set
    const ConvexHullType calculateConvexHull(const IndexSet& set )
    {
      ConvexHullType hull;
      
    	//convert to vector of DPosition
			std::vector< DPosition<2> > points;				
			for (IndexSet::const_iterator it = set.begin(); it!=set.end(); ++it)
      {
      	points.push_back(getPeak(*it).getPosition());
      }
      //calculate convex hull
      hull = points;
      return hull;
    }

protected:
   /// Writes gnuplot output (only for debugging purposes)
    void writeGnuPlotFile_(IndexSet peaks, bool last,int nr_feat)
    {
        // write feature + surrounding region to file
        if (!last)
        {
            String gp_fname("plot.gp");
            std::ofstream gpfile( gp_fname.c_str(), std::ios_base::app  );

            String file    = "region" + String(nr_feat);

						if (nr_feat == 0)
							gpfile << "splot \'" << file << "\' w i title \"\" " << std::endl;
						else	
              gpfile << "replot \'" << file << "\' w i title \"\" " << std::endl;

            std::ofstream myfile(file.c_str()); // data file
            IndexSet::const_iterator citer = peaks.begin();

            while (citer != peaks.end())
            {
                myfile << getPeakRt(*citer) << " " << getPeakMz(*citer) << " " << getPeakIntensity(*citer) << std::endl;
                citer++;
            }
            myfile.close();
            gpfile.close();

        }
        else
        {

            String gp_fname("plot.gp");

            std::ofstream gpfile( gp_fname.c_str() , std::ios_base::app );
            gpfile << "pause -1 \'Please hit enter to continue....\' " << std::endl;
            gpfile.close();
        }

    }

    /// Container for peak data
    MapType map_;

    /// Flags indicating whether a peak is unused, a seed or inside a feature region
    FlagVector flags_;

    /// stores reference to the scan numbers for each peak.
    ScanIndexMSExperiment<MapType, MapType::PIterator > scan_index_;

    /// The found features in the LC/MS map
    FeatureVector features_;
};

namespace Internal
{
/// Iterator adapter that makes operator*()
/// return intensity of the corresponding peak
struct IntensityIterator : IndexSet::const_iterator
{
    IntensityIterator ( IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
            : IndexSet::const_iterator(iter),
            traits_(traits)
    {}
    FeaFiTraits::IntensityType operator * () const throw()
    {
        return traits_->getPeakIntensity( IndexSet::const_iterator::operator *() );
    }
protected:
    FeaFiTraits const * traits_;
};

/// Iterator adapter that makes operator*()
/// return mz of the corresponding peak
struct MzIterator : IndexSet::const_iterator
{
    MzIterator ( IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
            : IndexSet::const_iterator(iter),
            traits_(traits)
    {}
    FeaFiTraits::CoordinateType operator * () const throw()
    {
        return traits_->getPeakMz( IndexSet::const_iterator::operator *() );
    }
protected:
    FeaFiTraits const * traits_;
};

/// Iterator adapter that makes operator*()
/// return retention time of the corresponding peak
struct RtIterator : IndexSet::const_iterator
{
    RtIterator ( IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
            : IndexSet::const_iterator(iter),
            traits_(traits)
    {}
    FeaFiTraits::CoordinateType operator * () const throw()
    {
        return traits_->getPeakRt( IndexSet::const_iterator::operator *() );
    }
protected:
    FeaFiTraits const * traits_;
};

} // namespace Internal

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
