// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  class CVMappingTerm;
  
  /**
      @brief Representation of a CV Mapping rule used by CVMappings

      Representation of a controlled vocabulary mapping rule.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI CVMappingRule
  {
public:

    /// enum to specify the requirement level
    enum RequirementLevel
    {
      MUST = 0,
      SHOULD = 1,
      MAY = 2
    };

    /// enum to specify the combination operator
    enum CombinationsLogic
    {
      OR = 0,
      AND = 1,
      XOR = 2
    };

    /// Default constructor
    CVMappingRule();

    /// Copy constructor
    CVMappingRule(const CVMappingRule & rhs);

    /// Destructor
    virtual ~CVMappingRule();

    /// Assignment operator
    CVMappingRule & operator=(const CVMappingRule & rhs);

    /** @name Accessors
    */
    //@{
    /// sets the identifier of the rule
    void setIdentifier(const String & identifier);

    /// returns the identifier of the rule
    const String & getIdentifier() const;

    /// sets the path of the element, where this rule is allowed
    void setElementPath(const String & element_path);

    /// returns the path of the element, where this rule is allowed
    const String & getElementPath() const;

    /// sets the requirement level of this rule
    void setRequirementLevel(RequirementLevel level);

    /// returns the requirement level of this rule
    RequirementLevel getRequirementLevel() const;

    /// sets the combination operator of the rule
    void setCombinationsLogic(CombinationsLogic combinations_logic);

    /// returns the combinations operator of the rule
    CombinationsLogic getCombinationsLogic() const;

    /// sets the scope path of the rule
    void setScopePath(const String & path);

    /// returns the scope path of the rule
    const String & getScopePath() const;

    /// sets the terms which are allowed
    void setCVTerms(const std::vector<CVMappingTerm> & cv_terms);

    /// returns the allowed terms
    const std::vector<CVMappingTerm> & getCVTerms() const;

    /// adds a term to the allowed terms
    void addCVTerm(const CVMappingTerm & cv_terms);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVMappingRule & rhs) const;

    /// inequality operator
    bool operator!=(const CVMappingRule & rhs) const;
    //@}

protected:

    String identifier_;

    String element_path_;

    RequirementLevel requirement_level_;

    String scope_path_;

    CombinationsLogic combinations_logic_;

    std::vector<CVMappingTerm> cv_terms_;
  };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H
