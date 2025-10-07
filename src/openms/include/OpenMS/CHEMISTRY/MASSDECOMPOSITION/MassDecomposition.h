// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <map>

namespace OpenMS
{
  class String;

  /**
    @brief Class represents a decomposition of a mass into amino acids

    This class represents a mass decomposition into amino acids. A
    decomposition are amino acids given with frequencies which add
    up to a specific mass.
  */
  class OPENMS_DLLAPI MassDecomposition
  {
public:

    /**
      @name Constructors and destructors
    */
    //@{
    /// default constructor
    MassDecomposition();

    /// copy constructor
    MassDecomposition(const MassDecomposition& deco);

    /// constructor with String as parameter
    explicit MassDecomposition(const String& deco);
    //@}

    /**
      @name Operators and accessors
    */
    //@{
    /// assignment operator
    MassDecomposition& operator=(const MassDecomposition& rhs);

    /// adds the mass decomposition d to this object
    MassDecomposition& operator+=(const MassDecomposition& d);

    /// returns the decomposition as a string
    String toString() const;

    /// returns the decomposition as a string; instead of frequencies the amino acids are repeated
    String toExpandedString() const;

    /// adds this decomposition and the decomposition given and returns a new composition
    MassDecomposition operator+(const MassDecomposition& rhs) const;

    /// returns the max frequency of this composition
    Size getNumberOfMaxAA() const;
    //@}

    /**
      @name Predicates
    */
    //@{
    /// less than predicate
    bool operator<(const MassDecomposition& rhs) const;

    /// equality operator
    bool operator==(const String& deco) const;

    /// returns true if tag is contained in the mass decomposition
    bool containsTag(const String& tag) const;

    /// returns true if the mass decomposition if contained in this instance
    bool compatible(const MassDecomposition& deco) const;
    //@}

protected:
    std::map<char, Size> decomp_;
    Size number_of_max_aa_;
  };
} // namespace OpenMS

