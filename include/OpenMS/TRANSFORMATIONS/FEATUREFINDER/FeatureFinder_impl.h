// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_IMPL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_IMPL_H

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm_impl.h>

namespace OpenMS
{
	// This is documented in the declaration, see FeatureFinder.h
	template<class PeakType, class FeatureType>
	void FeatureFinder::run(const String& algorithm_name, MSExperiment<PeakType> const & input_map, FeatureMap<FeatureType> & features, const Param& param)
	{
		// Nothing to do if there is no data
		if (input_map.size()==0)
		{
			return;
		}
	
		// check input
		{	
			// We need updated ranges => check number of peaks
			if (algorithm_name != "mrm" && input_map.getSize()==0)
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__, "FeatureFinder needs updated ranges on input map. Aborting.");
			}

			// We need MS1 data only => check levels
			if (algorithm_name != "mrm" && (input_map.getMSLevels().size() != 1 || input_map.getMSLevels()[0] != 1 ))
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__, "FeatureFinder can only operate on MS level 1 data. Please do not use MS/MS data. Aborting.");
			}
			
			//Check if the peaks are sorted according to m/z
			for (Size s=0; s<input_map.size(); ++s)
			{
				if (input_map[s].size()==0) continue;
				DoubleReal last_pos = input_map[s][0].getMZ();
				for (Size p=1; p<input_map[s].size(); ++p)
				{
					if (input_map[s][p].getMZ()<last_pos)
					{
						throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__, "FeatureFinder can only operate spectra that contain the peaks ordered according to m/z. Aborting.");
					}
					last_pos = input_map[s][p].getMZ();
				}
			}
		}

		// initialize
		{
			// Resize peak flag vector
			flags_.resize(input_map.size());
			for (Size i=0; i<input_map.size(); ++i)
			{
				flags_[i].assign(input_map[i].size(), UNUSED);
			}
		}
		
		// do the work
		if (algorithm_name!="none")
		{
			FeatureFinderAlgorithm<PeakType, FeatureType>* algorithm = Factory<FeatureFinderAlgorithm<PeakType, FeatureType> >::create(algorithm_name);
			algorithm->setParameters(param);
			algorithm->setData(input_map,features,*this);
			algorithm->run();
			delete(algorithm);
		}
		
		//report RT apex spectrum index and native ID for each feature
		for (Size i=0; i<features.size(); ++i)
		{
			//index
			Size spectrum_index = input_map.RTBegin(features[i].getRT()) - input_map.begin();
			features[i].setMetaValue("spectrum_index",(UInt) spectrum_index);
			//native id
			String native_id = input_map[spectrum_index].getNativeID();
			features[i].setMetaValue("spectrum_native_id", native_id);
		}
	}

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_IMPL_H
