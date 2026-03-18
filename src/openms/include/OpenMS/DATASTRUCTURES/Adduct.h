// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <functional>

namespace OpenMS
{

  class OPENMS_DLLAPI Adduct
  {
public:

    typedef std::vector<Adduct> AdductsType;

    /// Default C'tor
    Adduct();

    /// C'tor with initial charge
    Adduct(Int charge);

    /// C'tor for all members
    Adduct(Int charge, Int amount, double singleMass, const String& formula, double log_prob, double rt_shift, const String& label = "");

    /// Increase amount of this adduct by factor @param[in] m
    Adduct operator*(const Int m) const;
    /// Add two adducts amount if they are equal (defined by equal formula)
    Adduct operator+(const Adduct& rhs);
    /// Add other adducts amount to *this (equal formula required!)
    void operator+=(const Adduct& rhs);


    /// Print the contents of an Adduct to a stream.
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Adduct& a);

    /// Comparator
    friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b);

    //@{ Accessors
    const Int& getCharge() const;

    void setCharge(const Int& charge);

    const Int& getAmount() const;
    void setAmount(const Int& amount);

    const double& getSingleMass() const;
    void setSingleMass(const double& singleMass);

    const double& getLogProb() const;
    void setLogProb(const double& log_prob);

    const String& getFormula() const;
    void setFormula(const String& formula);

    const double& getRTShift() const;
    const String& getLabel() const;

    // convert a ion string to adduct string with charge information (eg. ion_string = "Na1", charge = "1" --> "[M+Na]+")
    String toAdductString(const String& ion_string, const Int& charge);

    /// Convert to adduct string with explicit multiplier (e.g., mol_multiplier=2 -> "[2M+H]+")
    static String toAdductString(const String& ion_string, const Int& charge, Int mol_multiplier);
    //}

private:
    Int charge_; ///< usually +1
    Int amount_; ///< number of entities
    double singleMass_; ///< mass of a single entity
    double log_prob_; ///< log probability of observing a single entity of this adduct
    String formula_; ///< chemical formula (parsable by EmpiricalFormula)
    double rt_shift_; ///< RT shift induced by a single entity of this adduct (this is for adducts attached prior to ESI, e.g. labeling)
    String label_; ///< Label for this adduct (can be used to indicate heavy labels)

    String checkFormula_(const String& formula);

  };

} // namespace OpenMS

// Hash function specialization for Adduct
// Note: Only hash fields used in operator== (charge, amount, singleMass, log_prob, formula)
// Do NOT hash rt_shift_ and label_ as they are not compared in operator==
namespace std
{
  template<>
  struct hash<OpenMS::Adduct>
  {
    std::size_t operator()(const OpenMS::Adduct& a) const noexcept
    {
      std::size_t seed = OpenMS::hash_int(a.getCharge());
      OpenMS::hash_combine(seed, OpenMS::hash_int(a.getAmount()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(a.getSingleMass()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(a.getLogProb()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(a.getFormula()));
      return seed;
    }
  };
} // namespace std


