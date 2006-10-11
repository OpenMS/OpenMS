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

#include <OpenMS/DATASTRUCTURES/ScanIndex.h>
#include <OpenMS/DATASTRUCTURES/IndexSet.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>

#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <OpenMS/KERNEL/MSExperiment.h>

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

    typedef std::vector<Flag> FlagVector;

    typedef DRawDataPoint<2> PeakType;
    typedef DRawDataPoint<2>::IntensityType IntensityType;
    typedef DRawDataPoint<2>::CoordinateType CoordinateType;
    typedef DRawDataPoint<2>::PositionType PositionType;
    typedef DFeature<2>::ChargeType ChargeType;

    typedef PeakType::NthPositionLess< RT > RTless;
    typedef PeakType::NthPositionLess< MZ > MZless;

    typedef DPeakArray<2, DRawDataPoint<2> > PeakVector;
    typedef DFeatureMap<2> FeatureVector;
    typedef DFeature<2>::ConvexHullType ConvexHullType;
    typedef FeaFiModule::NoSuccessor NoSuccessor;

    /// standard constructor
    FeaFiTraits();

    /// destructor
    virtual ~FeaFiTraits();

    /// copy constructor
    FeaFiTraits(const FeaFiTraits& source);

    /// assignment operator
    virtual FeaFiTraits& operator = (const FeaFiTraits& source);

    /// fill the internal data structure using an instance of MSExperiment
    void setData(MSExperiment<DPeak<1> >& exp);

    /// set iterator range as data for FeatureFinder
    template <class ConstPeakIterator>
    void setData(ConstPeakIterator begin, ConstPeakIterator end)
    {
			for (ConstPeakIterator it=begin; it!=end;++it)
			{
				addSinglePeak(*it);
			}
			// sorts the peak data
			sortData_();
    }

    void addSinglePeak(const DRawDataPoint<2>& peak)
    {
			peaks_.push_back(peak);
			flags_.push_back(UNUSED);
    }

    /// non-mutable acess flag with index @p index .
    const Flag& getPeakFlag(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
			return flags_.at(index);
    }
    /// mutable acess flag with index @p index.
    Flag& getPeakFlag(const UnsignedInt index) throw (Exception::IndexOverflow)
    {
			return flags_.at(index);
    }

    /// acess peak with index @p index.
    const PeakType& getPeak(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
			return peaks_.at(index);
    }
    /// retrieve the number of peaks.
    const UnsignedInt getNumberOfPeaks()
    {
			return peaks_.size();
    }
	
	PeakVector& getAllPeaks() 
	{
		return peaks_;
	}
	
	ScanIndex<PeakVector> getScanIndex() 
	{
		return scan_index_;
	}
	
	

    /// acess intensity of peak with index @p index.
    const IntensityType& getPeakIntensity(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
			return peaks_.at(index).getIntensity();
    }
    /// acess m/z of peak with index @p index .
    const CoordinateType& getPeakMz(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
			return peaks_.at(index).getPosition()[MZ];
    }
    /// acess retention time of peak with index @p index.
    const CoordinateType& getPeakRt(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
			return peaks_.at(index).getPosition()[RT];
    }
    /// returns signal/noise ration of peak with index @p index
//     const double& getPeakSN(const UnsignedInt index) const throw (Exception::IndexOverflow)
//     {
// 			return sn_ratios_.at(index);
//     }
    /// acess scan number of peak with index @p index
    const UnsignedInt getPeakScanNr(const UnsignedInt index) const throw (Exception::IndexOverflow);

    /** @brief get index of next peak in m/z dimensio.

		\param index of the peak whose successor is requested
		\return index of the next peak 
    */
    UnsignedInt getNextMz(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of previous peak in m/z dimension.

		\param index of the peak whose predecessor is requested
		\return index of the previous peak
    */
    UnsignedInt getPrevMz(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of next peak in retention time dimension.
     
		\param index of the peak whose successor is requested
		\return index of the next peak
    */
    UnsignedInt getNextRt(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /** @brief get index of next peak in retiontion time dimension.

		\param index of the peak whose predecessor is requested
		\return index of the previous peak
    */
    UnsignedInt getPrevRt(const UnsignedInt index) const throw (Exception::IndexOverflow, NoSuccessor);

    /// run main loop
    const FeatureVector& run(const std::vector<BaseSeeder*>& seeders,
                             const std::vector<BaseExtender*>& extenders,
                             const std::vector<BaseModelFitter*>& fitters);

    /** @brief Calculate the convex hull of the peaks contained in @p set

    Uses the gift wrap algorithm 
    */
    const ConvexHullType calculateConvexHull(const IndexSet& set);


	 protected:

    /** @brief We sort the peaks according to their position.

		In 1D m/z, in the 2D case m/z and rt. That is,
		the peaks are first sorted by their rt value
		and peaks with equal rt (i.e. scan index) are 
		then sorted by m/z. In addition,
		we initialise the vector of scan indizes
		in order to retrieve quickly the scan number of a peak.
    */
    void sortData_();

    /// @todo Remove. Only for debugging purposes
    void writeGnuPlotFile_(IndexSet peaks, bool last,int nr_feat);

    /// Calculate area of a triangle (needed for gift wrap algorithm)
    inline double triangleArea_(IndexSet::const_iterator it0, IndexSet::const_iterator it1, IndexSet::const_iterator it2)
    {
			// triangle area via determinant: x0*y1+x1*y2+x2*y0-x2*y1-x1*y0-x0*y2
			return getPeakMz(*it0)*getPeakRt(*it1) + getPeakMz(*it1)*getPeakRt(*it2) + getPeakMz(*it2)*getPeakRt(*it0)
				- getPeakMz(*it2)*getPeakRt(*it1) - getPeakMz(*it1)*getPeakRt(*it0) - getPeakMz(*it0)*getPeakRt(*it2);
    }

    /// vector of peaks
    PeakVector peaks_;

    /// Flags indicating whether a peak is unused, a seed or inside a feature region
    FlagVector flags_;

    /// stores reference to the scan numbers for each peak.
    ScanIndex<PeakVector> scan_index_;

    /// The (hopefully) found features in the LC/MS map
	FeatureVector features_;

    /// Stores a the signal / noise ratio for each peak
    std::vector<double> sn_ratios_;

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
		struct PeakIterator : IndexSet::const_iterator
		{
			PeakIterator ( IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
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
