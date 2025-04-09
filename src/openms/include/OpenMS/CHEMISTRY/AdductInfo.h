// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS
{
  class OPENMS_DLLAPI AdductInfo
  {
  public:
    /**
      C'tor, to build a representation of an adduct.

      @param name Identifier as given in the Positive/Negative-Adducts file, e.g. 'M+2K-H;1+'
      @param adduct Formula of the adduct, e.g. '2K-H'
      @param charge The charge (must not be 0; can be negative), e.g. 1
      @param mol_multiplier Molecular multiplier, e.g. for charged dimers '2M+H;+1'

    **/
    AdductInfo(const String& name, const EmpiricalFormula& adduct, int charge, UInt mol_multiplier = 1);

    /// returns the neutral mass of the small molecule without adduct (creates monomer from nmer, decharges and removes the adduct (given m/z of [nM+Adduct]/|charge| returns mass of [M])
    double getNeutralMass(double observed_mz) const;

    /// returns the m/z of the small molecule with neutral mass @p neutral_mass if the adduct is added (given mass of [M] returns m/z of [nM+Adduct]/|charge|)
    double getMZ(double neutral_mass) const;

    /// returns the mass shift caused by this adduct if charges are compensated with protons
    double getMassShift(bool use_avg_mass = false) const;

    /// checks if an adduct (e.g.a 'M+2K-H;1+') is valid, i.e. if the losses (==negative amounts) can actually be lost by the compound given in @p db_entry.
    /// If the negative parts are present in @p db_entry, true is returned.
    bool isCompatible(const EmpiricalFormula& db_entry) const;

    /// get charge of adduct
    int getCharge() const;

    /// original string used for parsing
    const String& getName() const;

    /// sum formula of adduct itself. Useful for comparison with feature adduct annotation
    const EmpiricalFormula& getEmpiricalFormula() const;

    /// get molecular multiplier
    UInt getMolMultiplier() const;

    /// parse an adduct string containing a formula (must contain 'M') and charge, separated by ';'.
    /// e.g. M+H;1+
    /// 'M' can have multipliers, e.g. '2M + H;1+' (for a singly charged dimer)
    static AdductInfo parseAdductString(const String& adduct);

    /// equality operator
    bool operator==(const AdductInfo& other) const;

  private:
    /// members
    String name_; ///< arbitrary name, only used for error reporting
    EmpiricalFormula ef_; ///< Sum formula for the actual adduct e.g. 'H' in 2M+H;+1
    double mass_; ///< computed from ef_.getMonoWeight(), but stored explicitly for efficiency
    int charge_;  ///< negative or positive charge; must not be 0
    UInt mol_multiplier_; ///< Mol multiplier, e.g. 2 in 2M+H;+1
  };
}
