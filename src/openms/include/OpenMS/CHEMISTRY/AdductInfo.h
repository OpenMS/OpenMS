// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
    bool isCompatible(EmpiricalFormula db_entry) const;

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
