// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  struct LightTransition
  {
public:

    LightTransition() :
      fragment_charge(0)
    {
    }

    std::string transition_name;
    std::string peptide_ref;
    double library_intensity;
    double product_mz;
    double precursor_mz;
    int fragment_charge;
    bool decoy;
    bool detecting_transition;
    bool quantifying_transition;
    bool identifying_transition;

    int getProductChargeState() const
    {
      return fragment_charge;
    }

    bool isProductChargeStateSet() const
    {
      return !(fragment_charge == 0);
    }

    std::string getNativeID() const
    {
      return transition_name;
    }

    std::string getPeptideRef() const
    {
      return peptide_ref;
    }

    std::string getCompoundRef() const
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

  struct LightModification
  {
    int location;
    int unimod_id;
  };

  // A compound is either a peptide or a metabolite
  struct LightCompound
  {

    LightCompound() :
      drift_time(-1),
      charge(0)
    {
    }

    double drift_time;
    double rt;
    int charge;
    std::string sequence;
    std::vector<std::string> protein_refs;
    // Peptide group label (corresponds to MS:1000893, all peptides that are isotopic forms of the same peptide should be assigned the same peptide group label)
    std::string peptide_group_label;
    std::string gene_name;
    std::string id;

    // for metabolites
    std::string sum_formula;
    std::string compound_name;

    // By convention, if there is no (metabolic) compound name, it is a peptide 
    bool isPeptide() const
    {
      return compound_name.empty();
    }

    void setChargeState(int ch)
    {
      charge = ch;
    }

    int getChargeState() const
    {
      return charge;
    }

    void setDriftTime(double d)
    {
      drift_time = d;
    }

    double getDriftTime() const
    {
      return drift_time;
    }

    std::vector<LightModification> modifications;
  };

  struct LightProtein
  {
    std::string id;
    std::string sequence;
  };

  struct LightTargetedExperiment
  {
    LightTargetedExperiment() : compound_reference_map_dirty_(true) {}

    typedef LightTransition Transition;
    typedef LightCompound Peptide;
    typedef LightCompound Compound;
    typedef LightProtein Protein;

    std::vector<LightTransition> transitions;
    std::vector<LightCompound> compounds;
    std::vector<LightProtein> proteins;
    std::vector<LightTransition> & getTransitions()
    {
      return transitions;
    }

    const std::vector<LightTransition> & getTransitions() const
    {
      return transitions;
    }

    std::vector<LightCompound> & getCompounds()
    {
      return compounds;
    }

    const std::vector<LightCompound> & getCompounds() const
    {
      return compounds;
    }

    std::vector<LightProtein> & getProteins()
    {
      return proteins;
    }

    const std::vector<LightProtein> & getProteins() const
    {
      return proteins;
    }

    // legacy
    const LightCompound& getPeptideByRef(const std::string& ref)
    {
      return getCompoundByRef(ref);
    }

    const LightCompound& getCompoundByRef(const std::string& ref)
    {
      if (compound_reference_map_dirty_)
      {
        createPeptideReferenceMap_();
      }
      return *(compound_reference_map_[ref]);
    }

  private:

    void createPeptideReferenceMap_()
    {
      for (size_t i = 0; i < getCompounds().size(); i++)
      {
        compound_reference_map_[getCompounds()[i].id] = &getCompounds()[i];
      }
      compound_reference_map_dirty_ = false;
    }

    // Map of compounds (peptides or metabolites)
    bool compound_reference_map_dirty_;
    std::map<std::string, LightCompound*> compound_reference_map_;

  };

} //end Namespace OpenSwath

