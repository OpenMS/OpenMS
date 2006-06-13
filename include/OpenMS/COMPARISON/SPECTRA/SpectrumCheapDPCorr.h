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
// $Id: SpectrumCheapDPCorr.h,v 1.5 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#ifndef OPENMS_DATASTRUCTURES_HASHMAP_H
	#include <OpenMS/DATASTRUCTURES/HashMap.h>
#endif

namespace OpenMS
{

  /**
  	@brief SpectrumCheapDPCorr calculates an optimal alignment on stick spectra
  
	  to save computing time only Peak Pairs that could get a score > 0 are considered
	  which Peak Pairs could get scores > 0 ? <br>
	  Peaks get a score depending on the difference in position and the heights of the peaks <br>
	  pairs with positions that differ more than some limit get score 0
	  
	  \param variation maximum difference in position (in percent of the current m/z)
	      note that big values of variation ( 1 being the maximum ) result in consideration
	      of all possible pairings which has a running time of O(n²)
    \param int_cnt: how the peak heights are used in the score<br>
      0 = product<br>
      1 = sqrt(product)<br>
      2 = sum <br>
      3 = agreeing intensity<br>
    \param keeppeaks
      keep peaks without alignment partner in the consensus spectrum<br>
  	
  	@todo correct int_cnt sum (is that really normalizing? it is not!) (Andreas)
  */
  class SpectrumCheapDPCorr : public CompareFunctor
  {
  public:
    /** @brief standard constructor <br>*/
    SpectrumCheapDPCorr();

    /** @brief copy constructor <br>*/
    SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source);

    /** @brief destructor <br>*/
    ~SpectrumCheapDPCorr();

    /** @brief assignment operator <br> */
    SpectrumCheapDPCorr& operator=(const SpectrumCheapDPCorr& source);

    static FactoryProduct* create() { return new SpectrumCheapDPCorr();}

		static const String getName()
		{
			return "SpectrumCheapDPCorr";
		}

    String info() const;

    double operator()(const ClusterSpectrum& csa ,const ClusterSpectrum& csb)const;

    /** @brief return consensus spectrum from last funtion call operator<br> */
    const MSSpectrum< DPeak<1> >& lastconsensus() const;

		HashMap<Size, Size> getPeakMap() const;

    /** @brief set weighting of <i>csb</i> for consensus from next function call operator<br> */
    void setFactor(double f);
  private:

    /** @brief O(n^2) dynamical programming <br> */
    double dynprog_(const MSSpectrum< DPeak<1> >& , const MSSpectrum< DPeak<1> >& , int, int, int, int ) const;

    /** @brief similarity of two peaks <br> */
    double comparepeaks_(double posa, double posb, double inta, double intb) const ;

    static const String info_;

    /** @brief consensus spectrum of the last comparison <br> */
    mutable MSSpectrum< DPeak<1> > lastconsensus_;

    /** @brief should peaks with no alignment partner be kept in the consensus? <br> */
    bool keeppeaks_;

    /** @brief weighting factor for the next consensus spectrum <br> */
    mutable double factor_;

		mutable HashMap<Size, Size> peak_map_;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
