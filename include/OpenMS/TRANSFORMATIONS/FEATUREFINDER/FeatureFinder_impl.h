// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
	void FeatureFinder::run(const String& algorithm_name, MSExperiment<PeakType> & input_map, FeatureMap<FeatureType> & features, const Param& param, const FeatureMap<FeatureType>& seeds)
	{
		// Nothing to do if there is no data
		if ((algorithm_name != "mrm" && input_map.empty()) || (algorithm_name == "mrm" && input_map.getChromatograms().empty()))
		{
		  features.clear(true);
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
      if (!input_map.isSorted(true))
      {
        LOG_WARN << "Input map is not sorted by RT and m/z! This is done now, before applying the algorithm!" << std::endl;
        input_map.sortSpectra(true);
        input_map.sortChromatograms(true);
      }
			for (Size s=0; s<input_map.size(); ++s)
			{
				if (input_map[s].empty()) continue;
				if (input_map[s][0].getMZ()<0)
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__, "FeatureFinder can only operate on spectra that contain peaks with positive m/z values. Filter the data accordingly beforehand! Aborting.");
				}
			}
		}

		// initialize
		if (algorithm_name!="mrm" && algorithm_name!="centroided")
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
			algorithm->setSeeds(seeds);
			algorithm->run();
			delete(algorithm);
		}
		
		if (algorithm_name!="mrm") // mrm  works on chromatograms; the next section is only for conventional data
		{
			//report RT apex spectrum index and native ID for each feature
			for (Size i=0; i<features.size(); ++i)
			{
				//index
				Size spectrum_index = input_map.RTBegin(features[i].getRT()) - input_map.begin();
				features[i].setMetaValue("spectrum_index", spectrum_index);
				//native id
				if (spectrum_index < input_map.size())
				{
					String native_id = input_map[spectrum_index].getNativeID();
					features[i].setMetaValue("spectrum_native_id", native_id);
				}
				else
				{
					/// @todo that happens sometimes using IsotopeWaveletFeatureFinder (Rene, Marc, Andreas, Clemens)
					std::cerr << "FeatureFinderAlgorithm_impl, line=" << __LINE__ << "; FixMe this cannot be, but happens" << std::endl;
				}
			}
		}
	}

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_IMPL_H
