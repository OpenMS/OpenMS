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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <map>
#include <vector>

#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
#include <OpenMS/METADATA/ProteinHit.h>

using std::map;
using std::vector;

namespace OpenMS
{
  LabelFreeLabeler::LabelFreeLabeler()
    : BaseLabeler()
  {
  }

  LabelFreeLabeler::~LabelFreeLabeler()
  {
  }

  void LabelFreeLabeler::preCheck(Param & /* param */) const
  {
    // no specific requirements
  }

  // merge all channels into the first one
  // no further influence of the simulation process needed
  void LabelFreeLabeler::setUpHook(FeatureMapSimVector & features)
  {
    if(features.size() == 1) return;
    else
    {
      LOG_INFO << "Merging input FASTA files into one. Intensities will be summed up if duplicates occur.";
      FeatureMapSim final_map = mergeProteinIdentificationsMaps_(features);
      features.clear();
      features.push_back(final_map);
    }
  }

  /// Labeling between digestion and rt simulation
  void LabelFreeLabeler::postDigestHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between RT and Detectability
  void LabelFreeLabeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void LabelFreeLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void LabelFreeLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void LabelFreeLabeler::postRawMSHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  void LabelFreeLabeler::postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)
  {

  }
} // namespace OpenMS
