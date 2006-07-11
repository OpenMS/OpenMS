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
// $Maintainer:  $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H
#define OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H

#include <vector>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/BinnedSparseVector.h>


namespace OpenMS
{
  typedef double intensity;
  /**
	  @brief Bin Representation of a 1D MSSpectrum (stick spectrum)
	  
	  Bin Dimensions: <br>
	    BinSize = size of the bins <br>
	    BinSpread > 0 adds peaks to more bins <br>
	      expample: <br>
	      spread 0:  | =>   # <br>
	      spread 1:  | =>  ### <br>
	      spread 2:  | =? ##### <br> <br>
	  since small binsizes tend to produce very sparse Vectors ( < 1%)<br>
	  a sparse Vector is used<br>
  
  */
  class BinnedRep
  {
  public:
    typedef BinnedSparseVector::const_iterator const_iterator;
    typedef BinnedSparseVector::iterator iterator;
    //typedef std::vector<double>::iterator iterator;
    //typedef std::vector<double>::const_iterator const_iterator;
    /** @brief standard constructor <br> */
    BinnedRep();

    /** @brief constructor with declaration of bin dimensions<br> */
    BinnedRep(const double, const uint = 0);

    /** @brief copy constructor <br> */
    BinnedRep(const BinnedRep& source);

    /** @brief destructor <br> */
    ~BinnedRep();

    /** @brief assignment operator <br> */
    BinnedRep& operator = (const BinnedRep& source);

    /** @name readonly accessors */
    //@[
    uint id() const{ return id_;}
    double getBinSize() const{ return binsize_;}
    double max() const{return end_;}
    double min() const{return begin_;}
    double getRetention() const { return retention_;}
    unsigned int getBinSpread() const {return spread_;}
    double getParentmz() const { return parent_m_z_;}
    unsigned int getPrecursorPeakCharge() const {return precursorpeakcharge_;}
    unsigned int size() const {return bins_.size();}
    String str() const ;
    //@}

    /** @name access to individual bins */
    //@{
    iterator begin() { return bins_.begin();}
    iterator end() { return bins_.end();}
    const_iterator begin() const{ return bins_.begin();}
    const_iterator end() const{ return bins_.end();}
    intensity operator[] (int) const ;
    //@}
    /** @brief fill bins with stick spectrum <br> */
    friend const MSSpectrum< DPeak<1> >& operator<<(BinnedRep&,const MSSpectrum< DPeak<1> >&);

    /** @brief scale all peaks from 0 to 1 */
    void normalize();

  private:
    //std::vector<intensity> bins_;
    /**
    sparse vector containing the summed intensity <br>
    */
    BinnedSparseVector bins_;

    /** @name bin dimensions */
    //@{
    double binsize_;
    unsigned int spread_;
    //@}

    double begin_;
    double end_;

    /** @name information about the source stick spectrum */
    //@{
    unsigned int id_;
    double retention_;
    double parent_m_z_;
    unsigned int precursorpeakcharge_;
    //@}

    void clear_();

  };
}

#endif //OPENMS_COMPARISON_CLUSTERING_BINNEDREP_H
