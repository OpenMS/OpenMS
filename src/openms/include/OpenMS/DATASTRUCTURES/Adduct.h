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

  /**
    @brief One elementary adduct contribution to a feature's ionization state.

    An Adduct describes a single chemical "piece" that can be added to (or
    removed from) a neutral molecule during ionization — for example one
    sodium cation, one ammonia adduct, or one lost proton. It is the
    @em building-block used by feature/charge deconvolution algorithms
    (FeatureDeconvolution, MetaboliteFeatureDeconvolution, Compomer,
    QTClusterFinder); several Adducts combine to describe the complete
    ion observed for one feature.

    @section Adduct_fields Fields

      - @c formula        : chemical formula of one entity, e.g. @c "Na",
                            @c "H", @c "NH3", @c "C2H3N". Must parse as a
                            neutral EmpiricalFormula; explicit charge in the
                            formula triggers a warning at construction.
      - @c amount         : how many copies of @c formula contribute
                            (stoichiometric count; e.g. 2 in @c "M+2K-H").
                            Negative amounts are accepted but emit a warning.
      - @c singleMass     : monoisotopic mass of one entity (one copy of
                            @c formula). The full mass contribution of this
                            Adduct is @c amount @c * @c singleMass.
      - @c charge         : per-entity charge of the adduct (e.g. @c +1 for
                            Na+, @c -1 for a proton loss). The ion's total
                            charge is @c amount @c * @c charge.
      - @c log_prob       : log-prior on observing this adduct, used by
                            deconvolution algorithms when scoring competing
                            hypotheses.
      - @c rt_shift       : optional retention-time shift induced by one
                            entity (non-zero for adducts attached @em prior
                            to ESI, e.g. covalent labels; zero for
                            in-source adducts like @c Na or @c K that do
                            not affect chromatographic behavior).
      - @c label          : free-form label, typically used to mark heavy
                            isotopic labels.

    @section Adduct_vs_adductinfo Relationship to AdductInfo

    AdductInfo (in CHEMISTRY/) represents an entire @em adduct specification
    — the parsed form of a string such as @c "[M+2K-H]+" with the molecular
    entity @em M, an n-mer multiplier, and a net charge. Adduct represents
    only a single piece of that specification (e.g. just the @c "2K" part
    or the @c "-H" part). AdductInfo is used for matching observed @em m/z
    against database compounds; Adduct is used to build up @em hypotheses
    about which combinations of pieces might explain a feature.

    @section Adduct_arithmetic Arithmetic

    @c operator*(Int) scales the amount: @c a*3 returns a copy of @c a
    with @c amount tripled. @c operator+ and @c operator+= merge two
    Adducts that share the same @c formula, summing their amounts (used
    when accumulating identical pieces of an ion's adduct composition).
    Mixing adducts with different formulae throws.

    @section Adduct_equality Equality and hashing

    Equality compares @c charge, @c amount, @c singleMass, @c log_prob,
    and @c formula. @c rt_shift and @c label are deliberately excluded so
    that two Adducts describing the same physical entity but tagged with
    different labels still compare equal — this matches the std::hash
    specialization for Adduct (see bottom of this header).

    @ingroup Datastructures
  */
  class OPENMS_DLLAPI Adduct
  {
public:

    /// Convenience typedef for a list of Adducts (used by feature deconvolution).
    typedef std::vector<Adduct> AdductsType;

    /**
      @brief Default constructor. All numeric fields zero-initialised; formula
      and label empty. Use this only as a placeholder; populate fields via
      setters before use.
    */
    Adduct();

    /**
      @brief Construct with only the per-entity charge set.

      All other fields are zero/empty. @c amount, @c singleMass,
      @c log_prob, and @c formula can be populated via setters; @c rt_shift
      and @c label can only be set via the full constructor.

      @param[in] charge Per-entity charge (e.g. @c +1 for Na+).
    */
    Adduct(Int charge);

    /**
      @brief Construct with all fields populated.

      @param[in] charge       Per-entity charge (e.g. @c +1 for Na+).
      @param[in] amount       Stoichiometric count of this adduct in the ion
                              (e.g. @c 2 for @c "2K" in @c "M+2K-H"). A
                              negative value is accepted but emits a warning.
      @param[in] singleMass   Monoisotopic mass of one entity in Dalton.
      @param[in] formula      Chemical formula of one entity (e.g. @c "Na",
                              @c "NH3"). Must parse as a neutral
                              EmpiricalFormula; warnings are logged for
                              empty, charged, or ambiguous single-element
                              formulas (e.g. @c "C2" with no other element).
      @param[in] log_prob     Log-prior on observing this adduct (default 0).
      @param[in] rt_shift     RT shift contributed by one entity in seconds
                              (default 0; non-zero for pre-ESI labels).
      @param[in] label        Free-form label, typically heavy-isotope tag.
    */
    Adduct(Int charge, Int amount, double singleMass, const String& formula, double log_prob, double rt_shift, const String& label = "");

    /**
      @brief Multiply the stoichiometric @c amount by an integer factor.

      Returns a copy with @c amount scaled; other fields unchanged. Used
      when expanding @c m*Adduct hypotheses (e.g. promoting one @c Na to
      @c 3 @c Na).

      @param[in] m Integer multiplier.
      @return A new Adduct with the scaled amount.
    */
    Adduct operator*(const Int m) const;

    /**
      @brief Sum two Adducts that describe the same chemical entity.

      Both operands must share the same @c formula; only the @c amount
      fields are summed. Used when accumulating identical pieces of an
      ion's adduct composition.

      @param[in] rhs The other Adduct (must have identical @c formula).
      @return A new Adduct with the summed amount.
      @throws const char* C-string literal when @p rhs has a different
              @c formula than @c *this.
    */
    Adduct operator+(const Adduct& rhs);

    /**
      @brief In-place version of @c operator+. Same precondition
      (matching formulas).

      @param[in] rhs The other Adduct (must have identical @c formula).
      @throws const char* C-string literal when @p rhs has a different
              @c formula than @c *this.
    */
    void operator+=(const Adduct& rhs);


    /// Print all fields except @c rt_shift and @c label to @p os
    /// (for diagnostics; not a parseable serialization).
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Adduct& a);

    /**
      @brief Equality on @c charge, @c amount, @c singleMass, @c log_prob,
      and @c formula. Deliberately ignores @c rt_shift and @c label so
      that adducts differing only in labeling tag still compare equal.
    */
    friend OPENMS_DLLAPI bool operator==(const Adduct& a, const Adduct& b);

    //@{ Accessors

    /// @brief Per-entity charge; @c +1 for typical cation adducts (Na+, H+).
    const Int& getCharge() const;

    /**
      @brief Set the per-entity charge.
      @param[in] charge Per-entity charge.
    */
    void setCharge(const Int& charge);

    /// @brief Stoichiometric count of this adduct in the ion.
    const Int& getAmount() const;

    /**
      @brief Set the stoichiometric count. Logs a warning on negative
      values but accepts them.
      @param[in] amount Stoichiometric count (typically @c >= 0).
    */
    void setAmount(const Int& amount);

    /// @brief Monoisotopic mass of one entity (Dalton).
    const double& getSingleMass() const;

    /**
      @brief Set the monoisotopic mass of one entity.
      @param[in] singleMass Monoisotopic mass of one entity in Dalton.
    */
    void setSingleMass(const double& singleMass);

    /// @brief Log-prior on observing this adduct.
    const double& getLogProb() const;

    /**
      @brief Set the log-prior on observing this adduct.
      @param[in] log_prob Log-prior (negative for unlikely adducts; @c 0
                          for no prior preference).
    */
    void setLogProb(const double& log_prob);

    /// @brief Chemical formula of one entity (canonicalised by
    /// EmpiricalFormula, e.g. @c "H2O" — not @c "OH2").
    const String& getFormula() const;

    /**
      @brief Set the chemical formula of one entity. The input is parsed
      as an EmpiricalFormula and re-emitted in canonical order; warnings
      are logged for empty, explicitly charged, or ambiguous-single-
      element formulas.
      @param[in] formula EmpiricalFormula-parseable string (e.g. @c "Na",
                         @c "NH3", @c "C2H3N").
    */
    void setFormula(const String& formula);

    /// @brief RT shift induced by one entity (seconds). Non-zero only for
    /// adducts attached prior to ESI (e.g. covalent labels).
    const double& getRTShift() const;

    /// @brief Free-form label (typically a heavy-isotope tag).
    const String& getLabel() const;

    /**
      @brief Render @p ion_string + @p charge as a canonical adduct string,
      assuming a monomer (n=1).

      Convenience overload that calls the three-argument version with
      @c mol_multiplier = 1.

      @param[in] ion_string Element/formula portion of the adduct,
                            EmpiricalFormula-parseable (e.g. @c "Na1",
                            @c "NH3").
      @param[in] charge     Signed net charge of the resulting ion.
      @return Adduct string in the form @c "[M+Na]+", @c "[M-H]-", etc.
    */
    String toAdductString(const String& ion_string, const Int& charge);

    /**
      @brief Render @p ion_string + @p charge as a canonical adduct string,
      with an explicit n-mer multiplier.

      Element/coefficient terms are sorted alphabetically by element symbol
      so equivalent formulas always produce identical strings.

      Examples:
      @verbatim
        ("Na1", +1, 1) -> "[M+Na]+"
        ("Na1", +1, 2) -> "[2M+Na]+"
        ("H1",  -1, 1) -> "[M-H]-"
        ("NH3", +1, 1) -> "[M+3H+N]+"   (canonical sort of NH3)
      @endverbatim

      @param[in] ion_string     EmpiricalFormula-parseable element string.
      @param[in] charge         Signed net charge of the resulting ion.
      @param[in] mol_multiplier N-mer multiplier (@c 1 = monomer,
                                @c 2 = dimer, ...).
      @return Canonical adduct string.
    */
    static String toAdductString(const String& ion_string, const Int& charge, Int mol_multiplier);
    //}

private:
    Int charge_;        ///< per-entity charge (e.g. +1 for Na+)
    Int amount_;        ///< stoichiometric count of this entity
    double singleMass_; ///< monoisotopic mass of one entity (Da)
    double log_prob_;   ///< log-prior on observing this adduct
    String formula_;    ///< chemical formula of one entity (EmpiricalFormula-parseable, canonicalised)
    double rt_shift_;   ///< RT shift induced by one entity (s); non-zero only for pre-ESI labels
    String label_;      ///< free-form tag (e.g. heavy-isotope label name)

    /// Validate and canonicalise @p formula. Warns on empty, charged, or
    /// ambiguous single-element inputs. Returns the canonical form.
    String checkFormula_(const String& formula);

  };

} // namespace OpenMS

// Hash function specialization for Adduct.
//
// Note: only the fields compared by operator== (charge, amount, singleMass,
// log_prob, formula) feed into the hash. rt_shift and label are deliberately
// excluded so two Adducts that compare equal also hash equal.
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
