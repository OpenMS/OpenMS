// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

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

