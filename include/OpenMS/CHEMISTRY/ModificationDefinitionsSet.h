// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MODIFICATIONDEFINITIONSSET_H
#define OPENMS_CHEMISTRY_MODIFICATIONDEFINITIONSSET_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

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
    ModificationDefinitionsSet(const ModificationDefinitionsSet & rhs);

    /// detailed constructor with comma separated list of modifications
    ModificationDefinitionsSet(const String & fixed_modifications, const String & variable_modifications = "");

    /// detailed constructor with StringLists
    ModificationDefinitionsSet(const StringList & fixed_modifications, const StringList & variable_modifications = StringList::create(""));

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
    void addModification(const ModificationDefinition & mod_def);

    /// sets the modification definitions
    void setModifications(const std::set<ModificationDefinition> & mod_defs);

    /** @brief set the modification definitions from a string

            The strings should contain a comma separated list of modifications. The names
            can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
            definitions and other specific definitions are given by the modifications themselves.
    */
    void setModifications(const String & fixed_modifications, const String & variable_modifications);

    /// same as above, but using StringList instead of comma separated strings
    void setModifications(const StringList & fixed_modifications, const StringList & variable_modifications);

    /// returns the stored modification definitions
    std::set<ModificationDefinition> getModifications() const;

    /// returns the stored fixed modification definitions
    const std::set<ModificationDefinition> & getFixedModifications() const;

    /// returns the stored variable modification definitions
    const std::set<ModificationDefinition> & getVariableModifications() const;

    /// return only the names of the modifications stored in the set
    std::set<String> getModificationNames() const;

    /// return only the names of the fixed modifications
    std::set<String> getFixedModificationNames() const;

    /// return only the names of the variable modifications
    std::set<String> getVariableModificationNames() const;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    ModificationDefinitionsSet & operator=(const ModificationDefinitionsSet & element);
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications
    bool isCompatible(const AASequence & peptide) const;

    /// equality operator
    bool operator==(const ModificationDefinitionsSet & rhs) const;

    /// inequality operator
    bool operator!=(const ModificationDefinitionsSet & rhs) const;
    //@}


protected:

    std::set<ModificationDefinition> variable_mods_;

    std::set<ModificationDefinition> fixed_mods_;

    Size max_mods_per_peptide_;
  };


} // namespace OpenMS

#endif
