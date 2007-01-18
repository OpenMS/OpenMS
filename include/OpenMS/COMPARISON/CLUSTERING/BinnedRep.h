// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch  $
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H
#define OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H

#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/BinnedSparseVector.h>


namespace OpenMS
{
  typedef double intensity;
  /**
	  @brief Binned Representation of a PeakSpectrum (stick spectrum)
	  
	  Bin Dimensions: <br>
	    BinSize = size of the bins <br>
	    BinSpread > 0 adds peaks to more bins <br>
	      expample: <br>
	      spread 0:  | =>   #      <br>
	      spread 1:  | =>  ###     <br>
	      spread 2:  | =? #####    <br>
	  since small binsizes tend to produce very sparse Vectors ( < 1%)
	  a sparse Vector is used
  
  */
  class BinnedRep
  {
  public:
    typedef BinnedSparseVector::const_iterator const_iterator;
    typedef BinnedSparseVector::iterator iterator;
		typedef BinnedSparseVector::const_iterator ConstIterator;
		typedef BinnedSparseVector::iterator Iterator;

		/** @name Constructors and destructors
		*/
   	//@{
		/// default constructor
    BinnedRep();

    /// detailed constructor with declaration of bin dimensions
    BinnedRep(const double, const uint = 0);

		/// detailed constructor with PeakSpectrum and bin settings
		BinnedRep(const PeakSpectrum& spec, double binsize = 1.0, uint binspread = 0);
		
    /// copy constructor
    BinnedRep(const BinnedRep& source);

    /// destructor
    virtual ~BinnedRep();
		//@}

		/** @name Accessors
		*/
		//@{
    /// assignment operator
    BinnedRep& operator = (const BinnedRep& source);

		/// returns the id of the spectrum
    uint id() const{ return id_;}

		/// returns the bin size of the spectrum
    double getBinSize() const{ return binsize_;}

		/// returns the maximum of the m/z values
    double max() const{return end_;}

		/// returns the minimum of the m/z values
    double min() const{return begin_;}

		/// returns the retention time
    double getRetention() const { return retention_;}

		/// returns the spreading of the bins
    unsigned int getBinSpread() const {return spread_;}

		/// returns the mass-to-charge ratio of the parent ion
    double getParentmz() const { return parent_m_z_;}

		/// returns the charge of the parent ion
    unsigned int getPrecursorPeakCharge() const {return precursorpeakcharge_;}

		/// returns the number of bins
    unsigned int size() const {return bins_.size();}

		/// converts it to a string
    String str() const;

		/// returns an iterator pointing at the first bin
    Iterator begin() { return bins_.begin();}

		/// returns an end iterator
    Iterator end() { return bins_.end();}

		/// returns an constant begin iterator
    ConstIterator begin() const{ return bins_.begin();}

		/// returns an constant end Iterator
    ConstIterator end() const{ return bins_.end();}

		/// access to a bin with bracket operator
    intensity operator[] (int) const;
    //@}
		//
    /// fill bins with stick spectrum
    friend void operator << (BinnedRep& bin_rep, const PeakSpectrum& spec);

    /// scale all peaks from 0 to 1
    void normalize();

  private:
    
    /// sparse vector containing the summed intensity <br>
    BinnedSparseVector bins_;

    /// size of the bins
    double binsize_;

		/// spreading of the bins
    unsigned int spread_;

		/// first m/z 
    double begin_;

		/// last m/z
    double end_;

    /// the spectrum id
    unsigned int id_;

		/// retention time of the spectrum
    double retention_;

		/// the m/z of the parent ion
    double parent_m_z_;

		/// charge of the parent ion
    unsigned int precursorpeakcharge_;

		/// clears all data
    void clear_();

  };
}

#endif //OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H
