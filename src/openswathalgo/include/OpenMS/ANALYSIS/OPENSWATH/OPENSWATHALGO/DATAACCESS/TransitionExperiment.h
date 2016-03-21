// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_TRANSITIONEXPERIMENT_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_TRANSITIONEXPERIMENT_H

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  struct OPENSWATHALGO_DLLAPI LightTransition
  {
public:
    std::string transition_name;
    std::string peptide_ref;
    double library_intensity;
    double product_mz;
    double precursor_mz;
    int product_charge;
    bool decoy;
    bool detecting_transition;
    bool quantifying_transition;
    bool identifying_transition;

    int getProductChargeState() const
    {
      return product_charge;
    }

    std::string getNativeID() const
    {
      return transition_name;
    }

    std::string getPeptideRef() const
    {
      return peptide_ref;
    }

    double getLibraryIntensity() const
    {
      return library_intensity;
    }

    void setLibraryIntensity(double l)
    {
      library_intensity = l;
    }

    double getProductMZ() const
    {
      return product_mz;
    }

    double getPrecursorMZ() const
    {
      return precursor_mz;
    }

    void setDetectingTransition (bool d)
    {
      detecting_transition = d;
    }

    bool isDetectingTransition() const
    {
      return detecting_transition;
    }

    void setQuantifyingTransition (bool q)
    {
      quantifying_transition = q;
    }

    bool isQuantifyingTransition() const
    {
      return quantifying_transition;
    }

    void setIdentifyingTransition (bool i)
    {
      identifying_transition = i;
    }

    bool isIdentifyingTransition() const
    {
      return identifying_transition;
    }
  };

  struct OPENSWATHALGO_DLLAPI LightModification
  {
    int location;
    std::string unimod_id;
  };

  struct OPENSWATHALGO_DLLAPI LightPeptide
  {
    double rt;
    int charge;
    std::string sequence;
    std::vector<std::string> protein_refs;
    // Peptide group label (corresponds to MS:1000893, all peptides that are isotopic forms of the same peptide should be assigned the same peptide group label)
    std::string peptide_group_label;
    std::string id;

    int getChargeState() const
    {
      return charge;
    }

    std::vector<LightModification> modifications;
  };

  struct OPENSWATHALGO_DLLAPI LightProtein
  {
    std::string id;
    std::string sequence;
  };

  struct OPENSWATHALGO_DLLAPI LightTargetedExperiment
  {
    LightTargetedExperiment() : peptide_reference_map_dirty_(true) {}

    typedef LightTransition Transition;
    typedef LightPeptide Peptide;
    typedef LightProtein Protein;

    std::vector<LightTransition> transitions;
    std::vector<LightPeptide> peptides;
    std::vector<LightProtein> proteins;
    std::vector<LightTransition> & getTransitions()
    {
      return transitions;
    }

    std::vector<LightPeptide> & getPeptides()
    {
      return peptides;
    }

    std::vector<LightProtein> & getProteins()
    {
      return proteins;
    }

    const LightPeptide& getPeptideByRef(const std::string& ref)
    {
      if (peptide_reference_map_dirty_)
      {
        createPeptideReferenceMap_();
      }
      return *(peptide_reference_map_[ref]);
    }

  private:

    void createPeptideReferenceMap_()
    {
      for (size_t i = 0; i < getPeptides().size(); i++)
      {
        peptide_reference_map_[getPeptides()[i].id] = &getPeptides()[i];
      }
      peptide_reference_map_dirty_ = false;
    }

    bool peptide_reference_map_dirty_;
    std::map<std::string, LightPeptide*> peptide_reference_map_;

  };

} //end Namespace OpenSwath


#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_DATAACCESS_TRANSITIONEXPERIMENT_H
