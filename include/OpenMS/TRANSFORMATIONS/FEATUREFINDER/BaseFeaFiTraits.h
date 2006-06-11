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
// $Id: BaseFeaFiTraits.h,v 1.35 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEFEAFITRAITS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEFEAFITRAITS_H

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/DATASTRUCTURES/IndexSet.h>

#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <vector>
#include <ostream> 

namespace OpenMS
{

  /** @brief Abstract base class for FeatureFinder traits holding datastructures and main-loop 

   					Every derived class has to implement the static functions
      			"T* create()" and "const String getName()" (see FactoryProduct for details).
  */
  class BaseFeaFiTraits 
    : public FactoryProduct
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
    typedef std::vector<const Flag*> FlagRefVector;

    typedef DRawDataPoint<2> PeakType;
    typedef DRawDataPoint<2>::IntensityType IntensityType;
    typedef DRawDataPoint<2>::CoordinateType CoordinateType;
    typedef DRawDataPoint<2>::PositionType PositionType;
    typedef DFeature<2>::ChargeType ChargeType;

    typedef PeakType::NthPositionLess< RT > RTless;
    typedef PeakType::NthPositionLess< MZ > MZless;
    
    typedef DPeakArrayNonPolymorphic<2, DRawDataPoint<2> > PeakVector;
    typedef std::vector<const PeakType*> PeakRefVector;
    typedef DFeatureMap<2> FeatureVector;
		typedef DFeature<2>::ConvexHullType ConvexHullType;
    typedef FeaFiModule::NoSuccessor NoSuccessor;
    
    /// standard constructor 
    BaseFeaFiTraits();

    /// copy constructor 
    BaseFeaFiTraits(const BaseFeaFiTraits& source);

    /// destructor 
    virtual ~BaseFeaFiTraits();

    /// assignment operator 
    virtual BaseFeaFiTraits& operator = (const BaseFeaFiTraits& source);

    /// register all derived classes here
    static void registerChildren();

		/** @brief set extenders used in run

			Registers a vector of seeders in this instance
	 		of BaseFeaFiTraits. Moreover, each seeder receives
	 		a pointer to this traits class.
		*/
    void setSeeders(std::vector<BaseSeeder*> seeders);

		/** @brief set extenders used in run

			Registers a vector of extender classes in this instance
	   	of BaseFeaFiTraits. Moreover, each extender receives
	   	a pointer to this traits class.
		*/
    void setExtenders(std::vector<BaseExtender*> extenders);

    /** @brief set fitters used in run

			Registers a vector of extender classes in this instance
	   	of BaseFeaFiTraits. Moreover, each extender receives
	   	a pointer to this traits class.
		*/
    void setFitters(std::vector<BaseModelFitter*> fitters);

    /// set iterator range as data for FeatureFinder
    template <class ConstPeakIterator>
    void setData(ConstPeakIterator begin, ConstPeakIterator end) 
    {
    	  	
    	for (ConstPeakIterator it=begin; it!=end;it ++)
   		{
      	addSinglePeak(*it);
   		} 
   		
   		// sorts the peak data
   		sortData_();

    }
    
    /// sets the verbosity of the debug messages
    void setDebugLevel(UnsignedInt lvl); 
    
    /// get debug level
    const UnsignedInt& getDebugLevel() const; 
    
    /// get debug level
    UnsignedInt& getDebugLevel();
    
    /// sets instance number
    void setInstanceId(String instance); 
    
    /// returns the instance number (not mutable)
    const String& getInstanceId() const; 
    
    /// get instance number (mutable)
    String& getInstanceId();
            
    /// sets the stream to which the debug messages are written (default is std::cout)
    void setDebugStream(std::ostream* os);
    
    /// returns the debug stream
    std::ostream* getDebugStream();
           
    virtual void setData(const MSExperiment<DPeak<1> >& exp) = 0;
    
     /// add a single peak to internal datastructure  
    virtual void addSinglePeak(const DRawDataPoint<2>& peak) = 0;
    /// Non-mutable acess flag with Index @p index
    virtual const Flag& getPeakFlag(const UnsignedInt index) const throw (Exception::IndexOverflow) =0;
    /// Mutable acess flag with Index @p index 
    virtual Flag& getPeakFlag(const UnsignedInt index) throw (Exception::IndexOverflow) =0;

    /// acess range of flags through pointers
    virtual const FlagRefVector& getFlags(const IndexSet& index_set) throw (Exception::IndexOverflow) =0;

    /// Non-mutable acess to all flags 
    virtual const FlagVector& getAllFlags() const=0;
    /// Mutable acess to all flags 
    virtual FlagVector& getAllFlags()=0;

    /// acess peak with Index @p index 
    virtual const PeakType& getPeak(const UnsignedInt index) const throw (Exception::IndexOverflow) =0;
    /// acess range of peaks 
    virtual const PeakRefVector& getPeaks(const IndexSet& index_set) throw (Exception::IndexOverflow) =0;
    /// acess all peaks 
		virtual const PeakVector& getAllPeaks()=0;
    /// retrieve the number of peaks 
    virtual const UnsignedInt getNumberOfPeaks()=0;

    /// access intensity of peak with index @p index 
    virtual const IntensityType& getPeakIntensity(const UnsignedInt index) const throw(Exception::IndexOverflow) =0;
    /// access m/z of peak with index @p index
    virtual const CoordinateType& getPeakMz(const UnsignedInt index) const throw(Exception::IndexOverflow)=0;
    /// access retention time of peak with index @p index 
    virtual const CoordinateType& getPeakRt(const UnsignedInt index) const throw(Exception::IndexOverflow)=0;
		/// acess scan number of peak with index @p index 
    virtual const UnsignedInt getPeakScanNr(const UnsignedInt index) const throw (Exception::IndexOverflow)=0;

    /** @brief get index of next peak in m/z dimension
    
       	\param index of the peak whose successor is requested
      	\return index of the next peak 
    */
    virtual UnsignedInt getNextMz(UnsignedInt index) const throw(Exception::IndexOverflow, NoSuccessor)=0;

    /** @brief get index of previous peak in m/z dimension 
    
       \param index of the peak whose predecessor is requested
       \return index of the previous peak
    */
    virtual UnsignedInt getPrevMz(UnsignedInt index) const throw(Exception::IndexOverflow, NoSuccessor)=0;

    /** @brief get index of next peak in retention time dimension 
    
    		\param index of the peak whose successor is requested
       	\return index of the next peak
    */
    virtual UnsignedInt getNextRt(UnsignedInt index) const throw(Exception::IndexOverflow, NoSuccessor)=0;

    /** @brief get index of next peak in retiontion time dimension 
    
       \param index of the peak whose predecessor is requested
       \return index of the previous peak
    */
    virtual UnsignedInt getPrevRt(UnsignedInt index) const throw(Exception::IndexOverflow, NoSuccessor)=0;
    	
    /// run main loop using set seeders, extenders and fitters
    virtual const FeatureVector& run()=0;

		/** @brief Calculate the convex hull of the peaks contained in @p set
		  
				Uses the gift wrap algorithm 
		 
		 */
		const ConvexHullType& calculateConvexHull(const IndexSet& set);

  protected:
  	
  	/// sorts the peaks according to retention time and m/z 
  	virtual void sortData_() = 0;

		/// Calculate area of a triangle (needed for gift wrap algorithm)
		inline double triangleArea_(IndexSet::const_iterator it0, IndexSet::const_iterator it1, IndexSet::const_iterator it2)
		{
			// triangle area via determinant: x0*y1+x1*y2+x2*y0-x2*y1-x1*y0-x0*y2
			return getPeakMz(*it0)*getPeakRt(*it1) + getPeakMz(*it1)*getPeakRt(*it2) + getPeakMz(*it2)*getPeakRt(*it0)
				   - getPeakMz(*it2)*getPeakRt(*it1) - getPeakMz(*it1)*getPeakRt(*it0) - getPeakMz(*it0)*getPeakRt(*it2);
		}
  	
    std::vector<BaseSeeder*> seeders_;
    std::vector<BaseExtender*> extenders_;
    std::vector<BaseModelFitter*> fitters_;

    FeatureVector features_;
        
     /// debug level (0: very limited output, 1: verbose output)
    UnsignedInt debug_;
    
    /// stream used for debug messages
    std::ostream* debug_stream_;
    
    /// instance number (needed for debugging messages)
    String instance_;
    
   
	private:
		/// @brief tempory storage of the calculated convex hull only accessible via call to calculateConvexHull_(IndexSet)
		ConvexHullType convex_hull_;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEFEAFITRAITS_H
