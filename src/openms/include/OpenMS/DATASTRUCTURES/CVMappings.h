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

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGS_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGS_H

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <map>
#include <vector>

namespace OpenMS
{
  class CVMappingRule;
  class CVReference;
  
  /**
    @brief Representation of controlled vocabulary mapping rules (for PSI formats)

    This file serves as object for the controlled vocabulary term usage definitions
    used in CV-Mapping files. All the supported attributes supported in the
    mapping file are supported by this class.

    @ingroup Format
  */
  class OPENMS_DLLAPI CVMappings
  {
public:

    /// Default constructor
    CVMappings();

    /// Copy constructor
    CVMappings(const CVMappings& rhs);

    /// Destructor
    virtual ~CVMappings();

    /// Assignment operator
    CVMappings& operator=(const CVMappings& rhs);

    /** @name Accessors
    */
    //@{
    /// sets the mapping rules of the mapping file
    void setMappingRules(const std::vector<CVMappingRule>& cv_mapping_rules);

    /// returns the mapping rules
    const std::vector<CVMappingRule>& getMappingRules() const;

    /// adds a mapping rule
    void addMappingRule(const CVMappingRule& cv_mapping_rule);

    /// sets the CV references
    void setCVReferences(const std::vector<CVReference>& cv_references);

    /// returns the CV references
    const std::vector<CVReference>& getCVReferences() const;

    /// adds a CV reference
    void addCVReference(const CVReference& cv_reference);
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if a CV reference is given
    bool hasCVReference(const String& identifier);

    /// equality operator
    bool operator==(const CVMappings& rhs) const;

    /// inequality operator
    bool operator!=(const CVMappings& rhs) const;
    //@}

protected:

    std::vector<CVMappingRule> mapping_rules_;

    std::map<String, CVReference> cv_references_;

    std::vector<CVReference> cv_references_vector_;
  };
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGS_H
