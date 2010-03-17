// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow, Andreas Bertsch $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
#define OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>
#include <set>
#include <stack>

namespace OpenMS
{

	/**	
  	@brief SpectraMerger Bla
		
		@todo Add Logger compatibility (Andreas)

		@htmlinclude OpenMS_SpectraMerger.parameters

  */
  class OPENMS_DLLAPI SpectraMerger
    : public DefaultParamHandler
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectraMerger();

    /// copy constructor 
    SpectraMerger(const SpectraMerger& source);

    /// destructor
    virtual ~SpectraMerger();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectraMerger& operator=(const SpectraMerger& source);
		// @}

		///
		template <typename ExperimentType> void mergeSpectraBlockWise(ExperimentType& exp)
		{
			return;
		}

		/// merges spectra with similar precursors
		template <typename ExperimentType> void mergeSpectraPrecursors(ExperimentType& exp)
		{
			DoubleReal mz_tolerance(param_.getValue("precursor_method:mz_tolerance"));
			DoubleReal rt_tolerance(param_.getValue("precursor_method:rt_tolerance"));
			DoubleReal mz_binning_width(param_.getValue("mz_binning_width"));
			DoubleReal mz_binning_unit(param_.getValue("mz_binning_unit"));

			typedef typename ExperimentType::ConstIterator const_exp_iter;
			Map<Size, std::vector<Size> > spectra_by_idx;
			Size count1(0);
			for (const_exp_iter it1 = exp.begin(); it1 != exp.end(); ++it1, ++count1)
			{
				if (it1->getMSLevel() == 1)
				{
					continue;
				}
				DoubleReal rt1(it1->getRT());
				if (it1->getPrecursors().size() == 0)
				{
					std::cerr << "SpectrumMerger::mergeSpectraPrecursors(): no precursor defined at spectrum: RT=" << rt1 << ", skipping!" << std::endl;
					continue;
				}
				else if (it1->getPrecursors().size() > 1)
				{
					std::cerr << "SpectrumMerger::mergeSpectraPrecursors(): multiple precursors defined at spectrum RT=" << rt1 << ", using only first one!" << std::endl;
				}
				
				DoubleReal precursor_mz1(it1->getPrecursors().begin()->getMZ());

				Size count2(count1);
				for (const_exp_iter it2 = it1 + 1; it2 != exp.end(); ++it2, ++count2)
				{
					if (it2->getMSLevel() == 1)
					{
						continue;
					}
					DoubleReal rt2(it2->getRT());
					if (it1->getPrecursors().size() == 0)
        	{
						continue;
        	}
					DoubleReal precursor_mz2(it2->getPrecursors().begin()->getMZ());
					
					if (fabs(precursor_mz1 - precursor_mz2) < mz_tolerance && fabs(rt1 - rt2) < rt_tolerance)
					{
						spectra_by_idx[count1].push_back(count2);
					}
				}
			}
	
			// identify which spectra are merged 	
			std::set<Size> used_spectra;
			Map<Size, std::vector<Size> > spectra_to_merge;
			for (Map<Size, std::vector<Size> >::ConstIterator it = spectra_by_idx.begin(); it != spectra_by_idx.end(); ++it)
			{
				if (used_spectra.find(it->first) != used_spectra.end())
				{
					continue;
				}
				used_spectra.insert(it->first);
				std::stack<Size> steak;
				for (std::vector<Size>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
				{
					if (used_spectra.find(*sit) == used_spectra.end())
					{
						steak.push(*sit);
					}
				}
	
				spectra_to_merge[it->first] = std::vector<Size>();
				while (steak.size() != 0)
				{
					Size spec_idx = steak.top();
					spectra_to_merge[it->first].push_back(spec_idx);
					steak.pop();
					for (std::vector<Size>::const_iterator sit = spectra_by_idx[spec_idx].begin(); sit != spectra_by_idx[spec_idx].end(); ++sit)
					{
						if (used_spectra.find(*sit) == used_spectra.end())
						{
							steak.push(*sit);
						}
					}
				}
			}


			// merge spectra	
			for (Map<Size, std::vector<Size> >::ConstIterator it = spectra_to_merge.begin(); it != spectra_to_merge.end(); ++it)
			{
				PeakSpectrum all_peaks = exp[it->first];			
				for (std::vector<Size>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
				{
					for (PeakSpectrum::ConstIterator pit = exp[*sit].begin(); pit != exp[*sit].end(); ++pit)
					{
						all_peaks.push_back(*pit);
					}
				}
				all_peaks.sortByPosition();

  			PeakSpectrum consensus_spec;
		  	consensus_spec.setMSLevel(2);
		  	Peak1D old_peak = *all_peaks.begin();
		  	for (PeakSpectrum::ConstIterator it = (++consensus_spec.begin()); it != consensus_spec.end(); ++it)
		  	{
		    	if (fabs(old_peak.getMZ() - it->getMZ()) < mz_binning_width) // TODO use unit
		    	{
		      	old_peak.setIntensity(old_peak.getIntensity() + it->getIntensity());
		    	}
   		 		else
    			{
      			consensus_spec.push_back(old_peak);
     		 		old_peak = *it;
    			}
  			}
			}


			return;
		}


		void mergeSpectraBlockWisePeakMap(PeakMap& exp);

		/// merges spectra with similar precursors
		void mergeSpectraPrecursorsPeakMap(PeakMap& exp);
		// @}
	
  };
	
}
#endif //OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
