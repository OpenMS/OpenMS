// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Controlled Vocabulary Reference

      Reference to a controlled vocabulary, defined in the first section of a mapping file.

      @ingroup Datastructures
  */

  class OPENMS_DLLAPI CVReference
  {
public:

    /// Default constructor
    CVReference();

    /// Copy constructor
    CVReference(const CVReference& rhs);

    /// Destructor
    virtual ~CVReference();

    /// Assignment operator
    CVReference& operator=(const CVReference& rhs);

    /** @name Accessors
    */
    //@{
    /// sets the name of the CV reference
    void setName(const std::string& name);

    /// returns the name of the CV reference
    const std::string& getName() const;

    /// sets the CV identifier which is referenced
    void setIdentifier(const std::string& identifier);

    /// returns the CV identifier which is referenced
    const std::string& getIdentifier() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVReference& rhs) const;

    /// inequality operator
    bool operator!=(const CVReference& rhs) const;
    //@}


protected:

    std::string name_;

    std::string identifier_;
  };


} // namespace OpenMS

