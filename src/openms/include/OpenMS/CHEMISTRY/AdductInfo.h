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
  /**
    @brief Mass-arithmetic helper for ionizing adducts (e.g. [M+H]+, [2M-H]-,
    [M+2K-H]+).

    An AdductInfo describes how a neutral compound @em M is observed as an ion
    in a mass spectrum: which atoms are added or lost, how many copies of @em M
    are clustered (n-mer multiplier), and what net charge state the resulting
    ion carries. The class provides:
      - parsing of the textual adduct notation used throughout OpenMS
        (see @ref AdductInfo_grammar below),
      - forward conversion @em M (neutral mass) → observed @em m/z
        (#getMZ), and the inverse @em m/z → @em M (#getNeutralMass),
      - the bare mass shift contributed by the adduct alone (#getMassShift,
        useful for charge/adduct grouping in feature deconvolution),
      - compatibility checks against a candidate database compound
        (#isCompatible).

    All conversions correctly account for electron mass: a positive charge
    means the ion is missing @c |charge| electrons relative to the neutral
    form, a negative charge means it carries @c |charge| extra electrons.
    The EmpiricalFormula stored internally is always the @em neutral atomic
    composition of the adduct part — the charge is tracked separately.

    @section AdductInfo_grammar Adduct string grammar

    Adduct strings are parsed by #parseAdductString and follow the form
    @verbatim
      [<n>]M[<sign><k><Formula>...]; <charge_magnitude><sign>
    @endverbatim
    where
      - @c [\<n\>]M is the molecular ion ("M" optionally prefixed by a small
        integer multiplier for n-mers, e.g. @c 2M for dimers),
      - each @c \<sign\>\<k\>\<Formula\> term adds (@c +) or subtracts (@c -) a
        stoichiometric amount of a neutral fragment, with optional integer
        coefficient @c \<k\> (default 1),
      - the trailing @c ;\<n\>\<sign\> gives the absolute charge value followed
        by a @c + or @c - sign for the polarity.

    Examples (all valid):
    @verbatim
      M+H;1+              ; [M+H]+
      M-H;1-              ; [M-H]-
      M+Na;1+             ; sodium adduct
      M+2K-H;1+           ; potassium adduct with proton loss
      2M+H;1+             ; protonated dimer
      2M+CH3CN+Na;1+      ; sodiated acetonitrile dimer adduct
      M;1+                ; bare positive molecular ion (rare)
    @endverbatim

    Whitespace is allowed and stripped during parsing. Mixed double signs
    such as @c ++ or @c +- inside the formula part are rejected. The
    semicolon between formula and charge is mandatory.

    @section AdductInfo_invariants Invariants

      - @c charge is non-zero (a charge of 0 throws on construction).
      - @c mol_multiplier is at least 1.
      - The stored EmpiricalFormula has @c getCharge() == 0; charge is held
        separately because EmpiricalFormula's own (de)protonation accounting
        is unsuitable for non-protic adducts like sodium or potassium.

    @section AdductInfo_example Example

    @code
      AdductInfo a = AdductInfo::parseAdductString("M+H;1+");
      double mz = a.getMZ(180.0634);   // glucose [M+H]+: ~181.0707
      double m  = a.getNeutralMass(mz); // round-trips to ~180.0634

      AdductInfo dimer = AdductInfo::parseAdductString("2M+Na;1+");
      double dimer_mz  = dimer.getMZ(180.0634); // sodiated glucose dimer
    @endcode

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI AdductInfo
  {
  public:
    /**
      @brief Construct an AdductInfo from already-parsed parts.

      Use #parseAdductString if you have the textual notation; this
      constructor is for callers that already have the components on hand
      (e.g. when building an adduct programmatically or restoring from
      a database).

      @param[in] name             Textual identifier (typically the original
                                  adduct string, used for diagnostics).
      @param[in] adduct           Neutral EmpiricalFormula of the adduct part
                                  alone (e.g. @c "H" for [M+H]+, or the empty
                                  formula for the bare [M]+). Must have
                                  @c getCharge() == 0; the charge of the ion
                                  is given separately via @p charge.
      @param[in] charge           Signed ion charge. Must be non-zero;
                                  positive for cations, negative for anions.
      @param[in] mol_multiplier   N-mer multiplier (1 = monomer, 2 = dimer,
                                  ...). Must be at least 1.
      @throws Exception::InvalidParameter if @p charge is 0, @p mol_multiplier
             is 0, or @p adduct already carries a non-zero charge.
    */
    AdductInfo(const String& name, const EmpiricalFormula& adduct, int charge, UInt mol_multiplier = 1);

    /**
      @brief Compute the neutral monomer mass M from an observed m/z.

      Inverse of #getMZ. Given an observed @em m/z assumed to come from
      this adduct, returns the neutral mass of one molecule @em M (i.e.
      with the n-mer collapsed and the adduct composition removed). The
      computation accounts for the electron mass that distinguishes a
      neutral atom from its ionized form.

      @param[in] observed_mz Observed m/z (must be > 0).
      @return Neutral monomer mass of @em M in Dalton.
    */
    double getNeutralMass(double observed_mz) const;

    /**
      @brief Compute the observed m/z of [nM+Adduct]^(charge) from a
      candidate neutral monomer mass.

      Forward direction; inverse of #getNeutralMass. Useful when scoring a
      database candidate against an observed peak.

      @param[in] neutral_mass Neutral mass of one monomer @em M in Dalton.
      @return Predicted m/z of the adduct ion.
    */
    double getMZ(double neutral_mass) const;

    /**
      @brief Net mass shift introduced by this adduct relative to the
      neutral monomer.

      Returns the mass difference between an observed ion @c [nM+Adduct]
      and @c n*M when expressed in @em proton-compensated form: a singly
      charged proton-adduct [M+H]+ yields a mass shift of zero (the
      "added" proton and the missing electron net out to one proton mass
      compensated by adding a hydrogen). Other adducts give a non-trivial
      shift; for example [M+Na]+ shifts by mass(Na) − mass(H).

      Primarily used by feature/charge deconvolution algorithms that
      group features by their mass differences.

      @param[in] use_avg_mass If true, use average atomic weights rather
                              than monoisotopic ones. Useful when comparing
                              against centroided low-resolution data.
      @return Mass shift in Dalton.
    */
    double getMassShift(bool use_avg_mass = false) const;

    /**
      @brief Check whether this adduct is physically compatible with a
      candidate compound.

      An adduct that removes atoms from the molecule (e.g. @c -H in
      [M-H]-) can only be observed if the candidate actually contains
      those atoms. This method returns true when every element subtracted
      by the adduct is present in sufficient quantity in @p db_entry.

      Plus-only adducts (no atoms subtracted) are trivially compatible
      with any compound and always return true.

      @param[in] db_entry EmpiricalFormula of the candidate compound to
                          test against.
      @return True iff @p db_entry could plausibly produce this adduct.
    */
    bool isCompatible(const EmpiricalFormula& db_entry) const;

    /**
      @brief Signed ion charge.
      @return Non-zero integer; positive for cations, negative for anions.
    */
    int getCharge() const;

    /**
      @brief Original textual identifier supplied to the constructor (or
      the parsed source string, when constructed via #parseAdductString).
      @return Const reference to the stored name; mainly for diagnostics.
    */
    const String& getName() const;

    /**
      @brief Neutral EmpiricalFormula of the adduct part alone
      (no @em M, no @em n-mer multiplier).

      Useful when comparing against per-feature adduct annotations.
      Charge is @em not encoded in the returned formula (it always has
      @c getCharge() == 0); use #getCharge for the ion's charge state.

      @return Const reference to the stored EmpiricalFormula.
    */
    const EmpiricalFormula& getEmpiricalFormula() const;

    /**
      @brief N-mer multiplier (1 for monomers, 2 for dimers, ...).
    */
    UInt getMolMultiplier() const;

    /**
      @brief Parse an adduct string into an AdductInfo.

      The accepted notation is described in @ref AdductInfo_grammar. The
      original string (after whitespace stripping) is stored as the
      AdductInfo's name and retrievable via #getName.

      @param[in] adduct Adduct string, e.g. @c "M+H;1+" or @c "2M+CH3CN+Na;1+".
      @return Parsed AdductInfo.
      @throws Exception::InvalidValue if the string is malformed (missing
             semicolon, missing charge sign, invalid operators, missing
             @em M, unexpected character).
      @throws Exception::ConversionError if a numeric component (charge,
              stoichiometric coefficient, n-mer multiplier) cannot be parsed
              as an integer.
      @throws Exception::ParseError if an adduct fragment formula cannot be
             parsed as an EmpiricalFormula.
      @throws Exception::InvalidParameter if the parsed adduct violates the
             constructor invariants (e.g. zero charge or a zero n-mer
             multiplier).
     */
    static AdductInfo parseAdductString(const String& adduct);

    /**
      @brief Equality on all four fields (name, formula, charge, n-mer
      multiplier). Two AdductInfos parsed from the same canonical string
      compare equal; two AdductInfos describing the same physical adduct
      but with different @em names do not.
    */
    bool operator==(const AdductInfo& other) const;

  private:
    String name_;             ///< original adduct string, kept for diagnostics
    EmpiricalFormula ef_;     ///< neutral formula of the adduct part (no M, no n-mer multiplier)
    double mass_;             ///< cached monoisotopic mass of @c ef_ (avoids recomputation)
    int charge_;              ///< signed ion charge; non-zero (positive = cation, negative = anion)
    UInt mol_multiplier_;     ///< n-mer multiplier (1 = monomer, 2 = dimer, ...)
  };
}
