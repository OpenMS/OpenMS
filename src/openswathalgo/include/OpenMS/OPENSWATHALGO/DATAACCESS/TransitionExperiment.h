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
#include <cstdint>
#include <functional>
#include <limits>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{

  /// @brief Compact enum for fragment ion types (replaces string storage)
  /// Reduces memory from ~32 bytes (std::string) to 1 byte
  enum class FragmentIonType : uint8_t
  {
    Unknown = 0,
    AIon,      ///< a-ion
    BIon,      ///< b-ion
    CIon,      ///< c-ion
    XIon,      ///< x-ion
    YIon,      ///< y-ion
    ZIon,      ///< z-ion
    ZPrimeIon, ///< z'-ion (z prime)
    ZDotIon,   ///< z.-ion (z dot)
    Precursor, ///< Precursor ion
    BMinusH2O, ///< b-ion with water loss
    YMinusH2O, ///< y-ion with water loss
    BMinusNH3, ///< b-ion with ammonia loss
    YMinusNH3, ///< y-ion with ammonia loss
    Empty = 255 ///< No fragment type set
  };

  /// @brief Convert fragment ion type string to enum
  /// @param s Fragment type string (e.g. "y", "b", "prec")
  /// @return Corresponding FragmentIonType enum value
  inline FragmentIonType stringToFragmentIonType(const std::string& s)
  {
    if (s.empty()) return FragmentIonType::Empty;
    if (s == "a") return FragmentIonType::AIon;
    if (s == "b") return FragmentIonType::BIon;
    if (s == "c") return FragmentIonType::CIon;
    if (s == "x") return FragmentIonType::XIon;
    if (s == "y") return FragmentIonType::YIon;
    if (s == "z") return FragmentIonType::ZIon;
    if (s == "z'") return FragmentIonType::ZPrimeIon;
    if (s == "z.") return FragmentIonType::ZDotIon;
    if (s == "prec" || s == "precursor") return FragmentIonType::Precursor;
    if (s == "b-H2O" || s == "b-H20") return FragmentIonType::BMinusH2O;
    if (s == "y-H2O" || s == "y-H20") return FragmentIonType::YMinusH2O;
    if (s == "b-NH3") return FragmentIonType::BMinusNH3;
    if (s == "y-NH3") return FragmentIonType::YMinusNH3;
    return FragmentIonType::Unknown;
  }

  /// @brief Convert fragment ion type enum to string
  /// @param t FragmentIonType enum value
  /// @return Fragment type string (e.g. "y", "b", "prec")
  inline std::string fragmentIonTypeToString(FragmentIonType t)
  {
    switch (t)
    {
      case FragmentIonType::Empty: return "";
      case FragmentIonType::AIon: return "a";
      case FragmentIonType::BIon: return "b";
      case FragmentIonType::CIon: return "c";
      case FragmentIonType::XIon: return "x";
      case FragmentIonType::YIon: return "y";
      case FragmentIonType::ZIon: return "z";
      case FragmentIonType::ZPrimeIon: return "z'";
      case FragmentIonType::ZDotIon: return "z.";
      case FragmentIonType::Precursor: return "prec";
      case FragmentIonType::BMinusH2O: return "b-H2O";
      case FragmentIonType::YMinusH2O: return "y-H2O";
      case FragmentIonType::BMinusNH3: return "b-NH3";
      case FragmentIonType::YMinusNH3: return "y-NH3";
      case FragmentIonType::Unknown:
      default: return "unknown";
    }
  }

  /// @brief Packed boolean flags for transitions
  /// Reduces memory from 4 bytes (4 separate bools) to 1 byte
  struct TransitionFlags
  {
    uint8_t decoy : 1;
    uint8_t detecting : 1;
    uint8_t quantifying : 1;
    uint8_t identifying : 1;
    uint8_t reserved : 4;

    TransitionFlags() : decoy(0), detecting(1), quantifying(1), identifying(0), reserved(0) {}
  };

  struct LightTransition
  {
    std::string transition_name;
    std::string peptide_ref;
    double library_intensity{};
    double product_mz{};
    double precursor_mz{};
    double precursor_im{-1};
    int8_t fragment_charge{};                    ///< Fragment charge (compact: range typically 1-8)
    TransitionFlags flags{};                     ///< Packed boolean flags

    // Additional fields for roundtrip TSV/PQP I/O
    int16_t fragment_nr{-1};                     ///< Fragment ion ordinal (e.g. 7 for y7)
    FragmentIonType fragment_type{FragmentIonType::Empty}; ///< Fragment ion type enum
    std::vector<std::string> peptidoforms;       ///< Peptidoforms for IPF

    // Legacy bool accessors for API compatibility
    bool getDecoy() const { return flags.decoy; }
    void setDecoy(bool d) { flags.decoy = d; }

    // Fragment type string accessors for backward compatibility
    std::string getFragmentType() const { return fragmentIonTypeToString(fragment_type); }
    void setFragmentType(const std::string& s) { fragment_type = stringToFragmentIonType(s); }

    /// Reconstruct annotation from fragment_type, fragment_nr, and fragment_charge
    /// Returns string like "b8", "y6^2", "prec", or empty if no fragment info set
    std::string getAnnotation() const
    {
      if (fragment_type == FragmentIonType::Empty)
      {
        return "";
      }
      if (fragment_type == FragmentIonType::Precursor)
      {
        return "prec";
      }
      std::string result = fragmentIonTypeToString(fragment_type);
      if (fragment_nr >= 0)
      {
        result += std::to_string(fragment_nr);
      }
      if (fragment_charge > 0)
      {
        result += "^" + std::to_string(static_cast<int>(fragment_charge));
      }
      return result;
    }

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

    double getPrecursorIM() const
    {
      return precursor_im;
    }

    void setDetectingTransition (bool d)
    {
      flags.detecting = d;
    }

    bool isDetectingTransition() const
    {
      return flags.detecting;
    }

    void setQuantifyingTransition (bool q)
    {
      flags.quantifying = q;
    }

    bool isQuantifyingTransition() const
    {
      return flags.quantifying;
    }

    void setIdentifyingTransition (bool i)
    {
      flags.identifying = i;
    }

    bool isIdentifyingTransition() const
    {
      return flags.identifying;
    }

    /// Equality operator - compares transition_name (consistent with hash)
    bool operator==(const LightTransition& rhs) const
    {
      return transition_name == rhs.transition_name;
    }

    bool operator!=(const LightTransition& rhs) const
    {
      return !(*this == rhs);
    }
  };

  struct LightModification
  {
    int location;
    int unimod_id;

    /// Equality operator - compares location and unimod_id (consistent with hash)
    bool operator==(const LightModification& rhs) const
    {
      return location == rhs.location && unimod_id == rhs.unimod_id;
    }

    bool operator!=(const LightModification& rhs) const
    {
      return !(*this == rhs);
    }
  };

  // A compound is either a peptide or a metabolite
  struct LightCompound
  {

    LightCompound() :
      drift_time(-1),
      rt(std::numeric_limits<double>::quiet_NaN()),
      rt_start(std::numeric_limits<double>::quiet_NaN()),
      rt_end(std::numeric_limits<double>::quiet_NaN()),
      charge(0)
    {
    }

    double drift_time;
    /// Retention time (seconds). NaN if unset — `ChromatogramExtractor::prepare_coordinates`
    /// throws if a non-negative `rt_extraction_window` is requested but `rt` is NaN.
    double rt;
    /// Optional retention time range, used by ChromatogramExtractor::prepare_coordinates
    /// when called with rt_extraction_window=NaN. Both NaN means "no range encoded";
    /// otherwise rt_start..rt_end define the per-compound extraction window.
    /// Populated by OpenSwathDataAccessHelper::convertTargetedExp from a heavyweight
    /// TargetedExperiment::Peptide that has two retention-time entries.
    double rt_start;
    double rt_end;
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

    // Additional fields for roundtrip TSV/PQP I/O
    std::string label_type;                      ///< Label type (e.g. "heavy" or "light")
    std::string smiles;                          ///< SMILES representation (metabolomics)
    std::string adducts;                         ///< Adducts (metabolomics)

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

    /// Equality operator - compares id (consistent with hash)
    bool operator==(const LightCompound& rhs) const
    {
      return id == rhs.id;
    }

    bool operator!=(const LightCompound& rhs) const
    {
      return !(*this == rhs);
    }
  };

  struct LightProtein
  {
    std::string id;
    std::string sequence;

    // Additional fields for roundtrip TSV/PQP I/O
    std::string uniprot_id;                      ///< UniProt identifier

    /// Equality operator - compares id (consistent with hash)
    bool operator==(const LightProtein& rhs) const
    {
      return id == rhs.id;
    }

    bool operator!=(const LightProtein& rhs) const
    {
      return !(*this == rhs);
    }
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

// Hash function specializations for OpenSwath types
namespace std
{
  /**
   * @brief Hash function for LightTransition.
   *
   * Hashes based on the transition_name which serves as the unique identifier.
   * This enables use in std::unordered_map and std::unordered_set.
   */
  template<>
  struct hash<OpenSwath::LightTransition>
  {
    std::size_t operator()(const OpenSwath::LightTransition& t) const noexcept
    {
      return std::hash<std::string>{}(t.transition_name);
    }
  };

  /**
   * @brief Hash function for LightCompound.
   *
   * Hashes based on the id which serves as the unique identifier.
   * This enables use in std::unordered_map and std::unordered_set.
   */
  template<>
  struct hash<OpenSwath::LightCompound>
  {
    std::size_t operator()(const OpenSwath::LightCompound& c) const noexcept
    {
      return std::hash<std::string>{}(c.id);
    }
  };

  /**
   * @brief Hash function for LightProtein.
   *
   * Hashes based on the id which serves as the unique identifier.
   * This enables use in std::unordered_map and std::unordered_set.
   */
  template<>
  struct hash<OpenSwath::LightProtein>
  {
    std::size_t operator()(const OpenSwath::LightProtein& p) const noexcept
    {
      return std::hash<std::string>{}(p.id);
    }
  };

  /**
   * @brief Hash function for LightModification.
   *
   * Hashes based on location and unimod_id using a proper combining function
   * to reduce collision probability.
   * This enables use in std::unordered_map and std::unordered_set.
   */
  template<>
  struct hash<OpenSwath::LightModification>
  {
    std::size_t operator()(const OpenSwath::LightModification& m) const noexcept
    {
      std::size_t seed = std::hash<int>{}(m.location);
      // Use standard hash combine formula (from boost) for better distribution
      seed ^= std::hash<int>{}(m.unimod_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
} // namespace std

