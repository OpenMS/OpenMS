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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_OPENSWATH_TEST_H
#define OPENMS_OPENSWATH_TEST_H

// #include "DataAccessImpl/OpenSWATH.h"
// #include "OpenSwath/DataAccess/DataStructures.h"
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

namespace OpenSWATH_Test
{

  using namespace OpenMS;

  // Above are the definitions using real chromatograms.
  // Below are the definitions using spectra objects to store the data.
  //
  // This is necessary as long as peak picking cannot be done on chromatograms
  // natively.  the classes needed at the moment are SavitzkyGolayFilter,
  // GaussFilter, PeakPickerHiRes -- SignalToNoiseEstimatorMedian seems to work
  // already.
#if 0
  typedef MSChromatogram<ChromatogramPeak> RichPeakChromatogram;
  typedef MRMTransitionGroup<MSChromatogram, ChromatogramPeak> MRMTransitionGroupType;
#else
  typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
  //typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, OpenSwath::LightTransition> MRMTransitionGroupType;
  typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, ReactionMonitoringTransition> MRMTransitionGroupType;
#endif

  void setup_MRMFeatureFinderScoring(MRMTransitionGroupType& transition_group, std::vector< RichPeakChromatogram >& picked_chroms)
  {

    // Load the chromatograms (mzML) and the meta-information (TraML)
    PeakMap exp;
    //OpenSwath::LightTargetedExperiment transitions;
    TargetedExperiment transitions;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), exp);

    {
      TargetedExperiment transition_exp_;
      TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
      //convertTargetedExp(transition_exp_, transitions);
      transitions = transition_exp_;
    }

    // add all the transitions to the peakgroup
    transition_group.setTransitionGroupID("mypeptide");
    transition_group.addTransition(transitions.getTransitions()[0], transitions.getTransitions()[2].getNativeID());
    transition_group.addTransition(transitions.getTransitions()[2], transitions.getTransitions()[0].getNativeID());
    transition_group.addTransition(transitions.getTransitions()[3], transitions.getTransitions()[4].getNativeID());

    // add all the chromatograms to the peakgroup
    {
      Chromatogram chromatogram_old = exp.getChromatograms()[1];
      RichPeakChromatogram chromatogram;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
      {   
        ChromatogramPeak peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      } 
      //cout << "Chromatogram 0 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
      chromatogram.setMetaValue("product_mz", 618.31);
      chromatogram.setNativeID(chromatogram_old.getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
    }
    {
      Chromatogram chromatogram_old = exp.getChromatograms()[0];
      RichPeakChromatogram chromatogram;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
      {   
        ChromatogramPeak peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      } 
      //cout << "Chromatogram 2 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
      chromatogram.setMetaValue("product_mz", 628.45);
      chromatogram.setNativeID(chromatogram_old.getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
    }
    {
      Chromatogram chromatogram_old = exp.getChromatograms()[4]; 
      RichPeakChromatogram chromatogram;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
      {   
        ChromatogramPeak peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      } 
      //cout << "Chromatogram 3 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
      chromatogram.setMetaValue("product_mz", 651.3);
      chromatogram.setNativeID(chromatogram_old.getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
    }

    // do peakpicking, e.g. find a peak at 3120 RT / 170 intensity in all the spectra
    for(Size k = 0; k < transition_group.getChromatograms().size(); k++)
    {
      RichPeakChromatogram picked_chrom;
      ChromatogramPeak peak;
      peak.setMZ(3120);
      peak.setIntensity(170);
      picked_chrom.push_back(peak);

      picked_chrom.getFloatDataArrays().clear();
      picked_chrom.getFloatDataArrays().resize(3);
      picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
      picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
      picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
      picked_chrom.getFloatDataArrays()[0].push_back(1000.0);
      picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
      picked_chrom.getFloatDataArrays()[2].push_back(3140.0);

      picked_chroms.push_back(picked_chrom);
    }

  }
}

#endif
