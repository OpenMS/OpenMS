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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <iostream>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Representation of a ribonucleotide (modified or unmodified)

      The available information is based on the Modomics database (http://modomics.genesilico.pl/modifications/).

      @see RibonucleotideDB
  */
  class OPENMS_DLLAPI Ribonucleotide
  {
    friend class RibonucleotideDB;
  public:

    enum TermSpecificityNuc
    {
      ANYWHERE,
      FIVE_PRIME,
      THREE_PRIME,
      NUMBER_OF_TERM_SPECIFICITY
    };

    /** @name Constructors
    */
    //@{
    /// Constructor
    Ribonucleotide(const String& name = "unknown ribonucleotide",
                   const String& code = ".",
                   const String& new_code = "",
                   const String& html_code = ".",
                   const EmpiricalFormula& formula = EmpiricalFormula(),
                   char origin = '.',
                   double mono_mass = 0.0,
                   double avg_mass = 0.0,
                   enum TermSpecificityNuc term_spec = ANYWHERE,
                   const EmpiricalFormula& baseloss_formula =
                   default_baseloss_);

    /// Copy constructor
    Ribonucleotide(const Ribonucleotide& ribo) = default;

    /// Destructor
    virtual ~Ribonucleotide();
    //@}

    /** @name Assignment
     */
    //@{
    /// assignment operator
    Ribonucleotide& operator=(const Ribonucleotide& ribo) = default;
    //@}


    /** @name Equality
     */
    //@{
    /// Equality operator
    bool operator==(const Ribonucleotide& ribonucleotide) const;
    //@}

    /** Accessors
     */
    //@{

    /// Return the short name
    const String getCode() const;

    /// Set the short name
    void setCode(const String& code);

    /// Get the name of the ribonucleotide
    const String getName() const;

    /// Set the name of the ribonucleotide
    void setName(const String& name);

    /// Get formula for the ribonucleotide
    const EmpiricalFormula getFormula() const;

    /// Set the empirical formula for the ribonucleotide
    void setFormula(const EmpiricalFormula& formula);

    /// Get the monoisotopic mass of the ribonucleotide
    double getMonoMass() const;

    /// Set the monoisotopic mass of the ribonucleotide
    void setMonoMass(double mono_mass);

    /// Set the average mass of the ribonucleotide
    double getAvgMass() const;

    /// Get the average mass of the ribonucleotide
    void setAvgMass(double avg_mass);

    /// Get the "new" (Modomics) code
    const String getNewCode() const;

    /// Set the "new" (Modomics) code
    void setNewCode(const String &new_code);

    /// ostream iterator to write the residue to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Ribonucleotide& ribo);

    /// Get the code of the unmodified base (e.g., "A", "C", ...)
    char getOrigin() const;

    /// Set the code of the unmodified base (e.g., "A", "C", ...)
    void setOrigin(char origin);

    /// Set the HTML (RNAMods) code
    String getHTMLCode() const;

    /// Get the HTML (RNAMods) code
    void setHTMLCode(const String& html_code);

    /// Get the terminal specificity
    enum TermSpecificityNuc getTermSpecificity() const;

    /// Set the terminal specificity
    void setTermSpecificity(enum TermSpecificityNuc term_spec);

    /// Get sum formula after loss of the nucleobase
    const EmpiricalFormula getBaselossFormula() const;

    /// Set the sum formula after loss of the nucleobase
    void setBaselossFormula(const EmpiricalFormula& formula);

    //@}

    /// Return true if this is a modified ribonucleotide and false otherwise
    bool isModified() const;

    /// Return whether this is an "ambiguous" modification (representing isobaric modifications on the base/ribose)
    bool isAmbiguous() const;

  protected:
    /// Default value for sum formula after nucleobase loss
    static const EmpiricalFormula default_baseloss_;

    String name_; ///< full name
    String code_; ///< short name
    String new_code_; ///< Modomics code
    String html_code_; ///< RNAMods code
    EmpiricalFormula formula_; ///< sum formula
    char origin_; ///< character of unmodified version of ribonucleotide
    double mono_mass_; ///< monoisotopic mass
    double avg_mass_; ///< average mass
    enum TermSpecificityNuc term_spec_; ///< terminal specificity
    EmpiricalFormula baseloss_formula_; ///< sum formula after loss of the nucleobase
  };

  /// Stream output operator
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Ribonucleotide& ribo);

  /// Dummy nucleotide used to represent 5' and 3' chain ends. Usually, just the phosphates.
  using RibonucleotideChainEnd = Ribonucleotide;

}
