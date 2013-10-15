// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_CHEMISTRY_MODIFICATIONDEFINITION_H
#define OPENMS_CHEMISTRY_MODIFICATIONDEFINITION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

namespace OpenMS
{
  /** @ingroup Chemistry

          @brief Representation of modification definition

          This class defines a modification type e.g. a input parameter of a search engine.
          The modification is defined using an unique name of the modification present
          in the modifications DB instance.
  */
  class OPENMS_DLLAPI ModificationDefinition
  {
public:

    /** @name Constructor and Destructors
    */
    //@{
    /// default constructor
    ModificationDefinition();

    /// copy constructor
    ModificationDefinition(const ModificationDefinition & rhs);

    /// detailed constructor specifying the modifications name
    explicit ModificationDefinition(const String & mod);

    /// destructor
    virtual ~ModificationDefinition();
    //@}

    /** @name Accessors
    */
    //@{
    /// sets the allowed position of the modification
    void setTermSpecificity(ResidueModification::Term_Specificity pos);

    /// returns the allowed position of the modification
    ResidueModification::Term_Specificity getTermSpecificity() const;

    /// sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
    void setFixedModification(bool fixed);

    /// returns if the modification if fixed true, else false
    bool isFixedModification() const;

    /// set the maximal number of occurences per peptide, unbound if 0
    void setMaxOccurences(UInt num);

    /// returns the maximal number of occurences per peptide
    UInt getMaxOccurences() const;

    /// returns the modification set
    String getModification() const;

    /// sets the modification, allowed are unique names provided by ModificationsDB
    void setModification(const String & modification);
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    ModificationDefinition & operator=(const ModificationDefinition & element);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const ModificationDefinition & rhs) const;

    /// inequality operator
    bool operator!=(const ModificationDefinition & rhs) const;

    /// less than operator for e.g. usage in maps; only mod FullIds are compared!
    bool operator<(const OpenMS::ModificationDefinition &) const;
    //@}

protected:

    /// allowed position
    ResidueModification::Term_Specificity term_spec_;

    /// the modification
    const ResidueModification * mod_;

    /// fixed (true) or variable (false)
    bool fixed_modification_;

    /// maximal number of occurences per peptide
    UInt max_occurences_;
  };

} // namespace OpenMS

#endif
