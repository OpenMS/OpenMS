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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
#define OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <vector>
#include <set>
#include <stack>

namespace OpenMS
{

	/**	
  	@brief Merges blocks of MS or MS2 spectra
		
    Parameter's are accessible via the DefaultParamHandler.

		@htmlinclude OpenMS_SpectraMerger.parameters

  */
  class OPENMS_DLLAPI SpectraMerger
    : public DefaultParamHandler
  {
  public:

    /// blocks of spectra (master-spectrum index to sacrifice-spectra(the ones being merged into the master-spectrum))
    typedef Map<Size, std::vector<Size> > MergeBlocks;

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

		// @name Merging functions
		// @{
		/// 
		template <typename MapType> void mergeSpectraBlockWise(MapType& exp)
		{
			IntList ms_levels = (IntList) (param_.getValue("block_method:ms_levels"));
      Int rt_block_size(param_.getValue("block_method:rt_block_size"));

      for (IntList::iterator it_mslevel = ms_levels.begin(); it_mslevel<ms_levels.end(); ++it_mslevel)
      {

        MergeBlocks spectra_to_merge;
        Size idx_block(0);
        SignedSize block_size_count(rt_block_size+1);
        Size idx_spectrum(0);
        for (typename MapType::const_iterator it1 = exp.begin(); it1 != exp.end(); ++it1)
			  {
				  if (Int(it1->getMSLevel()) == *it_mslevel)
				  {
            // block full
            if (++block_size_count >= rt_block_size)
            {
              block_size_count=0;
              idx_block = idx_spectrum;
            }
            else
            {
              spectra_to_merge[idx_block].push_back(idx_spectrum);
            }
          }

          ++idx_spectrum;
        }
        // check if last block had sacrifice spectra
        if (block_size_count==0)
        { //block just got initialized
          spectra_to_merge[idx_block] = std::vector<Size>();
        }

        // merge spectra, remove all old MS spectra and add new consensus spectra
        mergeSpectra_(exp, spectra_to_merge, *it_mslevel);
      }


      exp.sortSpectra();

			return;
		}

		/// merges spectra with similar precursors
		template <typename MapType> void mergeSpectraPrecursors(MapType& /*exp*/)
		{
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    /*
    idea: merge spectra with "similar" precursors to enhance S/N
    method: either leave #spectra identical and add neighbouring spectra if within deltas
            or do full clustering (single,complete linkage etc) and report only merged spectra
               with each spectrum being added to one cluster exlcusively

    untested and probably buggy code:


			DoubleReal mz_tolerance(param_.getValue("precursor_method:mz_tolerance"));
			DoubleReal rt_tolerance(param_.getValue("precursor_method:rt_tolerance"));


			typedef typename MapType::ConstIterator ConstExpIterator;
			typedef typename MapType::SpectrumType SpectrumType;
			Map<Size, std::vector<Size> > spectra_by_idx;
			Size count1(0);
      // iterate over spectra
			for (ConstExpIterator it1 = exp.begin(); it1 != exp.end(); ++it1, ++count1)
			{
        // only MS2 and above
				if (it1->getMSLevel() == 1)
				{
					continue;
				}
				DoubleReal rt1(it1->getRT());
				if (it1->getPrecursors().size() == 0)
				{
					LOG_DEBUG << "SpectrumMerger::mergeSpectraPrecursors(): no precursor defined at spectrum: RT=" << rt1 << ", skipping!" << std::endl;
					continue;
				}
				else if (it1->getPrecursors().size() > 1)
				{
					LOG_WARN << "SpectrumMerger::mergeSpectraPrecursors(): multiple precursors defined at spectrum RT=" << rt1 << ", using only first one!" << std::endl;
				}
				
				DoubleReal precursor_mz1(it1->getPrecursors().begin()->getMZ());

        // from current spectrum --> last spectrum
				Size count2(count1 + 1);
				for (ConstExpIterator it2 = it1 + 1; it2 != exp.end(); ++it2, ++count2)
				{
					if (it2->getMSLevel() == 1)	continue;
					if (it1->getPrecursors().size() == 0)	continue;
					
					DoubleReal rt2(it2->getRT());
					DoubleReal precursor_mz2(it2->getPrecursors().begin()->getMZ());
					
					if (fabs(precursor_mz1 - precursor_mz2) < mz_tolerance && fabs(rt1 - rt2) < rt_tolerance)
					{
						spectra_by_idx[count1].push_back(count2);
					}
				}
			}
	
			// identify which spectra are merged 	
			std::set<Size> used_spectra;
			MergeBlocks spectra_to_merge;
			for (MergeBlocks::ConstIterator it = spectra_by_idx.begin(); it != spectra_by_idx.end(); ++it)
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

      // merge spectra, remove all old MS2 spectra and add new consensus spectra
      mergeSpectra_(exp, spectra_to_merge, 2);

      exp.sortSpectra();
    
    */

			return;
		}


		// @}

    protected:
    
    /**
        @brief merges blocks of spectra of a certain level

        Merges spectra belonging to the same block, setting their MS level to @p ms_level.
        All old spectra of level @p ms_level are removed, and the new consensus spectra (one per block)
        are added.
        The resulting map is NOT sorted!

    */
    template <typename MapType>
    void mergeSpectra_(MapType& exp, const MergeBlocks& spectra_to_merge, const UInt ms_level)
    {
			DoubleReal mz_binning_width(param_.getValue("mz_binning_width"));
			String mz_binning_unit(param_.getValue("mz_binning_width_unit"));

      // merge spectra
			MapType merged_spectra;

			for (Map<Size, std::vector<Size> >::ConstIterator it = spectra_to_merge.begin(); it != spectra_to_merge.end(); ++it)
			{
        
        typename MapType::SpectrumType all_peaks = exp[it->first];			
        DoubleReal rt_average=all_peaks.getRT();

				for (std::vector<Size>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
				{
          rt_average+=exp[*sit].getRT();
					for (typename MapType::SpectrumType::ConstIterator pit = exp[*sit].begin(); pit != exp[*sit].end(); ++pit)
					{
						all_peaks.push_back(*pit);
					}
				}
				all_peaks.sortByPosition();
        rt_average/=it->second.size()+1;

  			typename MapType::SpectrumType consensus_spec;
        // todo: what about metainfo and precursor information?
		  	consensus_spec.setMSLevel(ms_level);
        consensus_spec.setRT(rt_average);

        if (all_peaks.size()==0) continue;
        else
        {
          typename MapType::PeakType old_peak = *all_peaks.begin();
          DoubleReal distance;
		  	  for (typename MapType::SpectrumType::ConstIterator it = (++all_peaks.begin()); it != all_peaks.end(); ++it)
		  	  {
            if (mz_binning_unit=="Da") distance=fabs(old_peak.getMZ() - it->getMZ());   //Da delta
            else distance= fabs(old_peak.getMZ() - it->getMZ())*1e6 / old_peak.getMZ(); //ppm delta

		    	  if (distance < mz_binning_width)
		    	  {
		      	  old_peak.setIntensity(old_peak.getIntensity() + it->getIntensity());
		    	  }
   		 		  else
    			  {
      			  consensus_spec.push_back(old_peak);
     		 		  old_peak = *it;
    			  }
  			  }
          consensus_spec.push_back(old_peak); // store last peak

				  merged_spectra.push_back(consensus_spec);
        }
			}

			// remove level "X" spectra and add consensus spectra
      exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<typename MapType::SpectrumType>(IntList::create(String(ms_level)), false)), exp.end());

			for (typename MapType::const_iterator it = merged_spectra.begin(); it != merged_spectra.end(); ++it)
			{
				exp.push_back(*it);
			}
    }
	
  };
	
}
#endif //OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
