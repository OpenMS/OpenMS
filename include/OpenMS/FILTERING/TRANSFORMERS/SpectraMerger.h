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
			Map<DoubleReal, std::vector<Size> > spectra_by_mz;
			Size count(0);
			for (const_exp_iter it = exp.begin(); it != exp.end(); ++it, ++count)
			{
				if (it->getMSLevel() == 1)
				{
					continue;
				}
				DoubleReal rt(it->getRT());
				if (it->getPrecursors().size() == 0)
				{
					std::cerr << "SpectrumMerger::mergeSpectraPrecursors(): no precursor defined at spectrum: RT=" << rt << ", skipping!" << std::endl;
				}
				else if (it->getPrecursors().size() > 1)
				{
					std::cerr << "SpectrumMerger::mergeSpectraPrecursors(): multiple precursors defined at spectrum RT=" << rt << ", using only first one!" << std::endl;
				}
				DoubleReal precursor_mz(it->getPrecursors().begin()->getMZ());
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
