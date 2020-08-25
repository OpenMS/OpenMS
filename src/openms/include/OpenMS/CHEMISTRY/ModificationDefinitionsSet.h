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
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CONCEPT/Types.h> // for "UInt"
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h> // for "StringList"
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <set>

namespace OpenMS
{
  /** @ingroup Chemistry

          @brief Representation of a set of modification definitions

          This class enhances the modification definitions as defined in the
          class ModificationDefinition into a set of definitions. This is also
          e.g. used as input parameters  in search engines.
  */
  class OPENMS_DLLAPI ModificationDefinitionsSet
  {
public:

    /** @name Constructor and Destructors
    */
    //@{
    /// default constructor
    ModificationDefinitionsSet();

    /// copy constructor
    ModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);

    /// detailed constructor with StringLists
    ModificationDefinitionsSet(const StringList& fixed_modifications, const StringList& variable_modifications);

    /// destructor
    virtual ~ModificationDefinitionsSet();
    //@}

    /** @name Accessors
    */
    //@{
    /// sets the maximal number of modifications allowed per peptide
    void setMaxModifications(Size max_mod);

    /// return the maximal number of modifications allowed per peptide
    Size getMaxModifications() const;

    /// returns the number of modifications stored in this set
    Size getNumberOfModifications() const;

    /// returns the number of fixed modifications stored in this set
    Size getNumberOfFixedModifications() const;

    /// returns the number of variable modifications stored in this set
    Size getNumberOfVariableModifications() const;

    /// adds a modification definition to the set
    void addModification(const ModificationDefinition& mod_def);

    /// sets the modification definitions
    void setModifications(const std::set<ModificationDefinition>& mod_defs);

    /** @brief set the modification definitions from a string

            The strings should contain a comma separated list of modifications. The names
            can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
            definitions and other specific definitions are given by the modifications themselves.
    */
    void setModifications(const String& fixed_modifications, const String& variable_modifications);

    /// same as above, but using StringList instead of comma separated strings
    void setModifications(const StringList& fixed_modifications, const StringList& variable_modifications);

    /// returns the stored modification definitions
    std::set<ModificationDefinition> getModifications() const;

    /// returns the stored fixed modification definitions
    const std::set<ModificationDefinition>& getFixedModifications() const;

    /// returns the stored variable modification definitions
    const std::set<ModificationDefinition>& getVariableModifications() const;

    /// returns only the names of the modifications stored in the set
    std::set<String> getModificationNames() const;

    /// populates the output lists with the modification names (use e.g. for ProteinIdentification::SearchParameters)
    void getModificationNames(StringList& fixed_modifications, StringList& variable_modifications) const;

    /// returns only the names of the fixed modifications
    std::set<String> getFixedModificationNames() const;

    /// returns only the names of the variable modifications
    std::set<String> getVariableModificationNames() const;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    ModificationDefinitionsSet& operator=(const ModificationDefinitionsSet& element);
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications
    bool isCompatible(const AASequence& peptide) const;

    /// equality operator
    bool operator==(const ModificationDefinitionsSet& rhs) const;

    /// inequality operator
    bool operator!=(const ModificationDefinitionsSet& rhs) const;
    //@}

    /**
       @brief Finds modifications in the set that match a given (delta) mass

       @param matches Matching modifications (output), sorted by mass error
       @param mass Query mass
       @param residue Query residue
       @param term_spec Query term specificity
       @param consider_variable Consider the variable modifications?
       @param consider_fixed Consider the fixed modifications?
       @param is_delta Is @p mass a delta mass (mass difference)?
       @param tolerance Numeric tolerance (in Da) for mass matching

       @throw Exception::IllegalArgument if both @p consider_variable and @p consider_fixed are false
    */
    void findMatches(std::multimap<double, ModificationDefinition>& matches, double mass, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY, bool consider_fixed = true, bool consider_variable = true, bool is_delta = true, double tolerance = 0.01) const;

    /// Infers the sets of defined modifications from the modifications present on peptide identifications
    void inferFromPeptides(const std::vector<PeptideIdentification>& peptides);

protected:

    std::set<ModificationDefinition> variable_mods_;

    std::set<ModificationDefinition> fixed_mods_;

    Size max_mods_per_peptide_;

    /// helper function for findMatches() - finds matching modifications in @p source and adds them to @p matches
    static void addMatches_(std::multimap<double, ModificationDefinition>& matches, double mass, const String& residue, ResidueModification::TermSpecificity term_spec, const std::set<ModificationDefinition>& source, bool is_delta, double tolerance);
 
  };

} // namespace OpenMS

