// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>
#include <map>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  enum class TransType
  {
    PEPTIDE, ///< Transition is a peptide transition
    NUCTIDE, ///< Transition is a nucleotide transition
    COMPOUND, ///< Transition is a metabolomics compound transition
    UNKNOWN ///< Transition type is unknown
  };

  struct LightTransition
  {
    std::string transition_name;
    std::string transition_ref;
    double library_intensity{};
    double product_mz{};
    double precursor_mz{};
    double precursor_im{-1};
    int fragment_charge{};
    bool decoy{};
    bool detecting_transition{};
    bool quantifying_transition{};
    bool identifying_transition{};
    TransType transition_type;

    int getProductChargeState() const
    {
      return fragment_charge;
    }

    bool isProductChargeStateSet() const
    {
      return !(fragment_charge == 0);
    }

    bool isPrecursorImSet() const
    {
      return !(precursor_im == -1);
    }

    std::string getNativeID() const
    {
      return transition_name;
    }

    std::string getTransRef() const
    {
      return transition_ref;
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

    double getPrecursorIM() const
    {
      return precursor_im;
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

    TransType getType() const
    {
      return transition_type;
    }
  };

  struct LightModification
  {
    int location;
    int unimod_id;
  };

  // A compound is either a peptide a nuctide or a metabolite
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

    // for oligos
    //TODO: Strictly speaking we should be able to use the same structure as for peptides, but we keep this separate for now
    std::vector<std::string> oligo_refs;
    std::string nuctide_group_label;
    

    // By convention, if there is no (metabolic) compound name, it is a peptide 
    bool isPeptide() const
    {
      return (compound_name.empty() && oligo_refs.empty());
    }

    bool isNuctide() const
    {
      return !oligo_refs.empty();
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

  struct LightOligo
  {
    std::string id;
    std::string sequence;
  };

  struct LightTargetedExperiment
  {
    LightTargetedExperiment() : compound_reference_map_dirty_(true) {}

    typedef LightTransition Transition;
    typedef LightCompound Peptide;
    typedef LightCompound Nuctide;
    typedef LightOligo Oligo;
    typedef LightCompound Compound;
    typedef LightProtein Protein;

    std::vector<LightTransition> transitions;
    std::vector<LightCompound> compounds;
    std::vector<LightProtein> proteins;
    std::vector<LightOligo> oligos;
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

    std::vector<LightOligo> & getOligos()
    {
      return oligos;
    }

    const std::vector<LightOligo> & getOligos() const
    {
      return oligos;
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
        createReferenceMap_();
      }
      return *(compound_reference_map_[ref]);
    }

  private:

    void createReferenceMap_()
    {
      for (size_t i = 0; i < getCompounds().size(); i++)
      {
        compound_reference_map_[getCompounds()[i].id] = &getCompounds()[i];
      }
      compound_reference_map_dirty_ = false;
    }

    // Map of compounds (peptides, nuctides, or metabolites)
    bool compound_reference_map_dirty_;
    std::map<std::string, LightCompound*> compound_reference_map_;

  };

} //end Namespace OpenSwath

