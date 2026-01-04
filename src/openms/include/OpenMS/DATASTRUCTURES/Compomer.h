// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/OpenMSConfig.h>

#include <functional>
#include <map>
#include <vector>

namespace OpenMS
{

  class String;

  /**
    @brief Holds information on an edge connecting two features from a (putative) charge ladder
    
    A Compomer represents the chemical composition difference between two mass spectrometry features.
    It stores information about the adducts (ions, molecules, or atoms) that explain the mass and
    charge differences between these features.
    
    The Compomer has two sides:
    - LEFT side: adducts subtracted from the first feature
    - RIGHT side: adducts added to the first feature
    
    This model allows representing the relationship between features that correspond to the same
    analyte but with different adduct compositions or charge states.
    
    The Compomer maintains metadata such as:
    - Net charge (difference between right and left sides)
    - Total mass difference
    - Probability score of this adduct combination
    - Expected RT shift caused by the adducts
    
    This class is used extensively in the feature decharging and adduct annotation processes.
    
    @ingroup Datastructures
  */
  class OPENMS_DLLAPI Compomer
  {
public:
    /**
      @brief Enumeration for specifying which side of the compomer to operate on
      
      - LEFT: The left side (adducts subtracted from the first feature)
      - RIGHT: The right side (adducts added to the first feature)
      - BOTH: Both sides of the compomer
    */
    enum SIDE {LEFT, RIGHT, BOTH};

    /// Type definition for one side of a compomer (maps adduct labels to Adduct objects)
    typedef std::map<String, Adduct> CompomerSide;
    
    /**
      @brief Container for both sides of a compomer
      
      Vector with exactly two elements:
      - [0] = left side (adducts subtracted)
      - [1] = right side (adducts added)
    */
    typedef std::vector<CompomerSide> CompomerComponents;

    /**
      @brief Default Constructor
      
      Initializes an empty compomer with zero net charge, mass, and probability.
    */
    Compomer();

    /**
      @brief Constructor with net-charge, mass, and probability

      @param[in] net_charge Net charge of the compomer (right side - left side)
      @param[in] mass Mass difference represented by the compomer
      @param[in] log_p Log probability of this adduct combination
    */
    Compomer(Int net_charge, double mass, double log_p);

    /**
      @brief Copy constructor

      @param[in] p Source compomer to copy from
    */
    Compomer(const Compomer& p);

    /**
      @brief Assignment Operator

      @param[in] source Source compomer to assign from
      @return Reference to this object
    */
    Compomer& operator=(const Compomer& source);

    /**
      @brief Add an adduct to a specific side of the compomer

      Adds the specified amount of the adduct to the given side and
      updates the compomer's properties (net charge, mass, etc.).

      @param[in] a The adduct to add
      @param[in] side Which side to add the adduct to (0=LEFT, 1=RIGHT)
    */
    void add(const Adduct& a, UInt side);

    /**
      @brief Determines if two compomers conflict with each other

      Checks if these two compomers can coexist for one feature by examining
      if they have conflicting adduct compositions on the specified sides.

      @param[in] cmp The other Compomer to compare against
      @param[in] side_this Which side of this compomer to check (0=LEFT, 1=RIGHT)
      @param[in] side_other Which side of the other compomer to check (0=LEFT, 1=RIGHT)
      @return True if the compomers conflict (cannot coexist), false otherwise
     */
    bool isConflicting(const Compomer& cmp, UInt side_this, UInt side_other) const;

    /**
      @brief Set a unique identifier for this compomer

      @param[in] id The unique ID to assign
    */
    void setID(const Size& id);
    
    /**
      @brief Get the unique identifier of this compomer
      
      @return The unique ID of this compomer
    */
    const Size& getID() const;
    
    /**
      @brief Get both sides (left and right) of this compomer
      
      @return Reference to the compomer components (left and right sides)
    */
    const CompomerComponents& getComponent() const;

    /**
      @brief Get the net charge of this compomer
      
      The net charge is calculated as the difference between the right and left sides.
      
      @return Net charge value
    */
    const Int& getNetCharge() const;

    /**
      @brief Get the total mass difference represented by this compomer
      
      @return Mass difference in Da
    */
    const double& getMass() const;

    /**
      @brief Get the sum of positive charges in this compomer
      
      @return Total positive charges
    */
    const Int& getPositiveCharges() const;

    /**
      @brief Get the sum of negative charges in this compomer
      
      @return Total negative charges
    */
    const Int& getNegativeCharges() const;

    /**
      @brief Get the log probability of this adduct combination
      
      Higher values indicate more likely combinations.
      
      @return Log probability value
    */
    const double& getLogP() const;

    /**
      @brief Get the expected retention time shift caused by this compomer
      
      @return Expected RT shift value
    */
    const double& getRTShift() const;

    /**
      @brief Get a string representation of all adducts in this compomer
      
      @return String representation of adducts on both sides
    */
    String getAdductsAsString() const;

    /**
      @brief Get a string representation of adducts on a specific side

      @param[in] side Which side to get adducts for (LEFT, RIGHT, or BOTH)
      @return String representation of adducts on the specified side
    */
    String getAdductsAsString(UInt side) const;

    /**
      @brief Check if the compomer contains only a single adduct on the specified side

      @param[out] a Output parameter that will contain the adduct if found
      @param[in] side Which side to check (LEFT or RIGHT)
      @return True if only a single adduct is present on the specified side
    */
    bool isSingleAdduct(Adduct& a, const UInt side) const;

    /**
      @brief Remove all adducts of type @p a

      Remove ALL instances of the given adduct,
      BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
    **/
    Compomer removeAdduct(const Adduct& a) const;

    /**
      @brief Remove all adducts of type @p a from @p side (LEFT or RIGHT)

      remove ALL instances of the given adduct from the given side (LEFT or RIGHT),
      BUT use the given adducts parameters (charge, logp, mass etc) to update the compomers members
    */
    Compomer removeAdduct(const Adduct& a, const UInt side) const;

    /**
      @brief Returns the adduct labels from @p side (LEFT or RIGHT)

      Get a list of labels for the @p side
      (useful for assigning channels (light, heavy etc) to features).
    */
    StringList getLabels(const UInt side) const;


    /**
      @brief Add a complete set of adducts to a specific side of the compomer

      @param[in] add_side The set of adducts to add
      @param[in] side Which side to add the adducts to (LEFT or RIGHT)
    */
    void add(const CompomerSide& add_side, UInt side);

    /**
      @brief Comparison operator for sorting compomers

      Sorts compomers by (in order of importance):
      1. Net charge
      2. Mass
      3. Probability

      @param[in] c1 First compomer to compare
      @param[in] c2 Second compomer to compare
      @return True if c1 should be ordered before c2
    */
    friend OPENMS_DLLAPI bool operator<(const Compomer& c1, const Compomer& c2);

    /**
      @brief Output stream operator for printing compomer contents

      @param[in,out] os Output stream to write to
      @param[in] cmp Compomer to print
      @return Reference to the output stream
    */
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Compomer& cmp);

    /**
      @brief Equality comparison operator

      @param[in] a First compomer to compare
      @param[in] b Second compomer to compare
      @return True if the compomers are equal
    */
    friend OPENMS_DLLAPI bool operator==(const Compomer& a, const  Compomer& b);

private:

    CompomerComponents cmp_;   ///< Adducts of left and right side
    Int net_charge_;           ///< Net charge (right - left)
    double mass_;              ///< Net mass (right - left)
    Int pos_charges_;          ///< Sum of positive charges
    Int neg_charges_;          ///< Sum of negative charges
    double log_p_;             ///< Log probability of this adduct combination
    double rt_shift_;          ///< Expected net RT shift (-shift_leftside + shift_rightside)
    Size id_;                  ///< Unique identifier for this compomer

  }; // \Compomer

} // namespace OpenMS

// Hash function specialization for Compomer
// Note: Only hash fields used in operator== (cmp_, net_charge_, mass_, pos_charges_, neg_charges_, log_p_, id_)
// Do NOT hash rt_shift_ as it is not compared in operator==
namespace std
{
  template<>
  struct hash<OpenMS::Compomer>
  {
    std::size_t operator()(const OpenMS::Compomer& c) const noexcept
    {
      std::size_t seed = OpenMS::hash_int(c.getNetCharge());
      OpenMS::hash_combine(seed, OpenMS::hash_float(c.getMass()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(c.getPositiveCharges()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(c.getNegativeCharges()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(c.getLogP()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(c.getID()));

      // Hash the compomer components (vector<map<String, Adduct>>)
      const auto& components = c.getComponent();
      OpenMS::hash_combine(seed, OpenMS::hash_int(components.size()));
      for (const auto& side : components)
      {
        OpenMS::hash_combine(seed, OpenMS::hash_int(side.size()));
        for (const auto& [key, adduct] : side)
        {
          OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(key));
          OpenMS::hash_combine(seed, std::hash<OpenMS::Adduct>{}(adduct));
        }
      }
      return seed;
    }
  };
} // namespace std

